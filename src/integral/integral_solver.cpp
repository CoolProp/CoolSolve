/**
 * @file integral_solver.cpp
 * @brief Time-marching orchestrator for the equation-based dynamic solver.
 *
 * Strategy: the integral equations `y = base + INTEGRAL(dydt, t, t0, tf)` are
 * *removed* from a working copy of the IR, producing an "algebraic subsystem"
 * in which the state variables `y` and the integration variable `t` have no
 * defining equation — they become free external values that the integrator
 * fixes each step. The algebraic `Solver` (reused unmodified) then resolves
 * the remaining equations (the derivative definitions `dydt = f(...)` and any
 * algebraic constraints) and the integrand values are read back as `f = dy/dt`.
 * See docs/integral_table.md §2 (Algorithm overview) and §3 (Architecture).
 */
#include "coolsolve/ees_string.h"
#include "coolsolve/integral/integral_solver.h"

#include "coolsolve/autodiff_node.h"
#include "coolsolve/parser.h"   // astToString
#include "coolsolve/variable_inference.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace coolsolve {

namespace {

std::string toLowerCopy(const std::string& s) { return eesFoldCase(s); }

/// The algebraic solve of an evaluation of the right-hand side failed (IntegralSolver::evaluateRHS).
struct AlgebraicSolveError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/// Does `expr` contain a top-level-or-nested INTEGRAL() call?
bool containsIntegralCall(const ExprPtr& expr) {
    if (!expr) return false;
    bool found = false;
    std::visit([&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, FunctionCall>) {
            std::string n = node.name;
            std::transform(n.begin(), n.end(), n.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (n == "integral") found = true;
        }
    }, expr->node);
    if (found) return true;
    std::visit([&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, UnaryOp>) {
            found = containsIntegralCall(node.operand);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            found = containsIntegralCall(node.left) || containsIntegralCall(node.right);
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            for (const auto& a : node.args) if ((found = containsIntegralCall(a))) return;
            for (const auto& na : node.namedArgs) if ((found = containsIntegralCall(na.second))) return;
        } else if constexpr (std::is_same_v<T, Variable>) {
            for (const auto& idx : node.indices) if ((found = containsIntegralCall(idx))) return;
        }
    }, expr->node);
    return found;
}

/// Variables of `ir` that (structurally) depend on one of `seeds`: the closure of the seeds
/// through the blocks of the analysis. A block depends on a seed when one of its inputs does
/// (an equation variable it does not determine); a loop (block of several variables) that
/// contains a dependent variable is dependent as a whole.
std::set<std::string, CaseInsensitiveLess> dependentVariables(
    const IR& ir, const StructuralAnalysisResult& analysis,
    const std::set<std::string, CaseInsensitiveLess>& seeds) {
    std::set<std::string, CaseInsensitiveLess> dep = seeds;
    const auto& eqs = ir.getEquations();
    bool changed = true;
    for (int pass = 0; changed && pass < 1000; ++pass) {
        changed = false;
        for (const auto& blk : analysis.blocks) {
            std::set<std::string, CaseInsensitiveLess> outputs(blk.variables.begin(), blk.variables.end());
            bool dynamic = false;
            if (blk.variables.size() > 1) {
                for (const auto& v : blk.variables)
                    if (dep.count(v)) { dynamic = true; break; }
            }
            for (int id : blk.equationIds) {
                if (dynamic) break;
                if (id < 0 || id >= static_cast<int>(eqs.size())) continue;
                for (const auto& v : eqs[id].variables)
                    if (!outputs.count(v) && dep.count(v)) { dynamic = true; break; }
            }
            if (!dynamic) continue;
            for (const auto& v : blk.variables)
                if (dep.insert(v).second) changed = true;
        }
    }
    return dep;
}

}  // namespace

// ============================================================================
// Construction — build the algebraic subsystem
// ============================================================================

IntegralSolver::IntegralSolver(const Program& program, const IR& ir,
                               const StructuralAnalysisResult& analysis,
                               const SolverOptions& options,
                               const LookupTableStore* lookupStore)
    : fullIr_(ir), solverOpts_(options), lookupStore_(lookupStore) {
    // The time-march loop reuses the algebraic `Solver` thousands of times
    // (several RHS evaluations per integration step). The per-block
    // `progressCallback` is designed for a single top-level solve (a handful
    // of block events); firing it inside every internal step floods the GUI's
    // SSE progress stream with tens of thousands of events and freezes the
    // browser. The parametric study follows the same convention — it does not
    // propagate `progressCallback` to its per-point solves (server.cpp ~1870).
    // The `cancelToken` is preserved so the Stop button stays responsive.
    solverOpts_.progressCallback = nullptr;

    problem_ = extractIntegralProblem(ir, analysis);

    // The $IntegralTable spec (if any), straight off the AST.
    IntegralTableSpec tableSpec;
    for (const auto& stmt : program.statements) {
        if (stmt->is<Directive>()) {
            const auto& d = stmt->as<Directive>();
            if (d.hasIntegralTableSpec) {
                tableSpec = d.integralTableSpec;
                break;
            }
        }
    }

    if (!problem_.hasIntegral) {
        problem_.valid = false;
        problem_.errorMessage = "no INTEGRAL() equations found";
        return;
    }
    if (!problem_.valid) return;   // solve() reports problem_.errorMessage

    // 1. Build the program of the algebraic subsystem: the equations declaring a user state
    //    (`y = base + INTEGRAL(...)`) are removed (the integrator owns `y`), the equations that
    //    hold INTEGRAL calls inside a larger expression are replaced by their copy with the
    //    auxiliary states, and the non-variable integrands get an equation of their own.
    //    The reduced system then pins the integrator-owned variables (states + integration
    //    variable) with "driver" equations (v = <literal>): they make it square (otherwise the
    //    free states/time fail the analyzer's squareness gate) and let the integrator set these
    //    values each step by mutating the shared NumberLiteral.
    //
    //    The plain additive form is only a state of its own when its base is constant along the
    //    trajectory (it is evaluated once, at t0). A base that depends on the integration variable
    //    or on a state through the algebraic equations is detected on the reduced system and the
    //    equation is then handled with an auxiliary state, like any other form.
    std::set<int> forceAuxiliary;
    for (int attempt = 0; attempt < 2; ++attempt) {
        if (attempt == 1) problem_ = extractIntegralProblem(ir, analysis, &forceAuxiliary);
        if (!problem_.valid) return;
        problem_.tableSpec = tableSpec;

        baseProgram_ = Program{};
        baseProgram_.sourceFilename = program.sourceFilename;
        std::map<std::pair<const Expression*, const Expression*>, const RewrittenEquation*> rewritten;
        for (const auto& rw : problem_.rewrittenEquations)
            rewritten[{rw.originalLhs.get(), rw.originalRhs.get()}] = &rw;
        size_t rewrittenFound = 0, declarationsDropped = 0;
        for (const auto& stmt : program.statements) {
            if (stmt->is<Equation>()) {
                const auto& eq = stmt->as<Equation>();
                auto it = rewritten.find({eq.lhs.get(), eq.rhs.get()});
                if (it != rewritten.end()) {
                    auto copy = std::make_shared<Statement>(*stmt);
                    copy->as<Equation>().lhs = it->second->lhs;
                    copy->as<Equation>().rhs = it->second->rhs;
                    baseProgram_.statements.push_back(copy);
                    ++rewrittenFound;
                    continue;
                }
                if (containsIntegralCall(eq.rhs) || containsIntegralCall(eq.lhs)) {
                    ++declarationsDropped;   // declaration of a user state
                    continue;
                }
            }
            baseProgram_.statements.push_back(stmt);
        }
        if (rewrittenFound != problem_.rewrittenEquations.size() ||
            declarationsDropped != problem_.integralEquationIds.size()) {
            problem_.valid = false;
            problem_.errorMessage =
                "INTEGRAL() calls inside a DUPLICATE loop are not supported: write the equations "
                "one by one (see docs/integral_table.md §7.1).";
            return;
        }
        for (const auto& [name, expr] : problem_.integrandEquations)
            baseProgram_.statements.push_back(makeEquation(makeVariable(name), expr));

        std::vector<std::string> driven;
        for (const auto& st : problem_.states) driven.push_back(st.name);
        driven.push_back(problem_.integrationVar);
        march_ = buildReduced({}, driven);

        // A user state whose base depends on t or on a state is not a plain initial value.
        forceAuxiliary.clear();
        if (attempt == 0 && march_->analysis->success) {
            std::set<std::string, CaseInsensitiveLess> seeds(driven.begin(), driven.end());
            const auto dynamic = dependentVariables(*march_->ir, *march_->analysis, seeds);
            for (const auto& st : problem_.states) {
                if (st.auxiliary) continue;
                for (const auto& v : st.baseVariables)
                    if (dynamic.count(v)) { forceAuxiliary.insert(st.equationId); break; }
            }
        }
        if (forceAuxiliary.empty()) break;
    }

    // The blocks the march solves at every evaluation of the right-hand side: those that depend on the
    // integration variable or on a state (the driver blocks included). The others are solved by the
    // complete solve at the initial point and keep their values.
    if (problem_.valid && march_ && march_->analysis->success) {
        std::set<std::string, CaseInsensitiveLess> seeds;
        for (const auto& st : problem_.states) seeds.insert(st.name);
        seeds.insert(problem_.integrationVar);
        const auto dynamic = dependentVariables(*march_->ir, *march_->analysis, seeds);
        const auto& blocks = march_->analysis->blocks;
        for (size_t b = 0; b < blocks.size(); ++b)
            if (std::any_of(blocks[b].variables.begin(), blocks[b].variables.end(),
                            [&](const std::string& v) { return dynamic.count(v) > 0; }))
                marchBlocks_.push_back(b);
    }
}

std::unique_ptr<IntegralSolver::ReducedSystem> IntegralSolver::buildReduced(
    const std::vector<StmtPtr>& extraStatements, const std::vector<std::string>& driven) const {
    auto sys = std::make_unique<ReducedSystem>();
    Program prog;
    prog.sourceFilename = baseProgram_.sourceFilename;
    prog.statements = baseProgram_.statements;
    for (const auto& st : extraStatements) prog.statements.push_back(st);
    for (const auto& name : driven) {
        ExprPtr lit = makeNumber(0.0);
        sys->drivers[name] = lit;
        prog.statements.push_back(makeEquation(makeVariable(name), lit));
    }

    // Build the reduced IR, run inference/initialisation on it, then copy the user's
    // guesses / loaded initials from the full IR.
    sys->ir = std::make_unique<IR>(IR::fromAST(prog));
    inferVariables(*sys->ir);
    initializeVariables(*sys->ir);
    for (const auto& [name, info] : sys->ir->getVariables()) {
        (void)info;
        if (const auto* full = fullIr_.getVariable(name)) {
            auto* mut = sys->ir->getVariableMutable(name);
            if (mut) {
                mut->guessValue = full->guessValue;
                mut->solutionValue = full->solutionValue;
                mut->lowerBound = full->lowerBound;
                mut->upperBound = full->upperBound;
            }
        }
    }

    // Structural analysis on the reduced (algebraic) system, and the algebraic solver
    // (reused at every time step).
    sys->analysis = std::make_unique<StructuralAnalysisResult>(StructuralAnalyzer::analyze(*sys->ir));
    sys->solver = std::make_unique<Solver>(*sys->ir, *sys->analysis, solverOpts_.coolpropConfig);
    if (lookupStore_) sys->solver->setLookupTableStore(lookupStore_);
    return sys;
}

// ============================================================================
// RHS evaluation and helpers
// ============================================================================

double IntegralSolver::evalBaseExpr(
    const ExprPtr& expr,
    const std::map<std::string, double, CaseInsensitiveLess>& vars) const {
    if (!expr) return 0.0;
    if (expr->is<NumberLiteral>()) return expr->as<NumberLiteral>().value;
    if (expr->is<Variable>()) {
        const auto& v = expr->as<Variable>();
        const std::string name = v.flattenedName();   // array elements are named `x[2]`
        auto it = vars.find(name);
        if (it != vars.end()) return it->second;
        // Fall back to the full IR's stored value (guesses / initials).
        const auto* info = fullIr_.getVariable(name);
        if (info) {
            if (info->solutionValue) return *info->solutionValue;
            if (info->guessValue) return *info->guessValue;
            return 1.0;
        }
        throw std::runtime_error("The initial value of a state refers to '" + name +
                                 "', which is not a variable of the model");
    }
    if (expr->is<UnaryOp>()) {
        const auto& u = expr->as<UnaryOp>();
        double v = evalBaseExpr(u.operand, vars);
        return (u.op == "-") ? -v : v;
    }
    if (expr->is<BinaryOp>()) {
        const auto& b = expr->as<BinaryOp>();
        double l = evalBaseExpr(b.left, vars), r = evalBaseExpr(b.right, vars);
        if (b.op == "+") return l + r;
        if (b.op == "-") return l - r;
        if (b.op == "*") return l * r;
        if (b.op == "/") return (r != 0.0) ? l / r : 0.0;
        if (b.op == "^") return std::pow(l, r);
    }
    if (expr->is<FunctionCall>()) {
        // Built-in math functions (abs, pi, sqrt, min, IF in both forms, ...) share the
        // evaluator of the algebraic solver; anything else throws "Unknown or unsupported
        // function" instead of silently evaluating to 0.
        std::vector<ADValue> args;
        for (const auto& a : expr->as<FunctionCall>().args)
            args.push_back(ADValue::constant(evalBaseExpr(a, vars), 0));
        return evaluateStandardFunction(expr->as<FunctionCall>().name, args).value;
    }
    return 0.0;
}

Eigen::VectorXd IntegralSolver::evaluateRHS(double t, const Eigen::VectorXd& y) {
    // Pin the integrator-owned variables: update the driver equations' literals
    // (so the explicit driver blocks resolve to these values) and set them as
    // algebraic guesses so any block reading them as an external sees the same.
    setDriverValue(*march_, problem_.integrationVar, t);
    march_->solver->setGuess(problem_.integrationVar, t);
    for (int i = 0; i < static_cast<int>(problem_.states.size()); ++i) {
        setDriverValue(*march_, problem_.states[i].name, y(i));
        march_->solver->setGuess(problem_.states[i].name, y(i));
    }

    std::string error;
    if (march_->solver->solveBlocks(marchBlocks_, solverOpts_, &error) != SolverStatus::Success) {
        throw AlgebraicSolveError("Algebraic solve failed at t=" + std::to_string(t) + ": " + error);
    }

    // f_i = derivative of state i, read from the algebraic solution.
    const auto& state = march_->solver->getEvaluator().getAllVariables();
    Eigen::VectorXd f(static_cast<int>(problem_.states.size()));
    for (int i = 0; i < f.size(); ++i) {
        const auto& sv = problem_.states[i];
        auto it = state.find(sv.integrandVar);
        if (it == state.end())
            throw std::runtime_error("The integrand '" + sv.integrandVar + "' of INTEGRAL() is not "
                                     "a variable of the model (nothing defines it)");
        f(i) = it->second;
    }
    return f;
}

double IntegralSolver::readDerivative(int /*stateIndex*/) const {
    // Reserved for future streaming use; evaluation happens in evaluateRHS.
    return 0.0;
}

void IntegralSolver::setDriverValue(ReducedSystem& sys, const std::string& name, double value) {
    auto it = sys.drivers.find(name);
    if (it != sys.drivers.end() && it->second && it->second->is<NumberLiteral>())
        it->second->as<NumberLiteral>().value = value;
}

// ============================================================================
// Symbolic limits: solve the model at the initial point and read them
// ============================================================================

bool IntegralSolver::resolveSymbolicInterval(double& t0, double& tf, double& step,
                                             std::string& error, DiagnosticCollector& diag) {
    // The initial point: t = lower limit, every state at its initial value. t is not driven by the
    // integrator here but defined by `t = <lower limit>`; each limit / step is a variable defined
    // by its expression, so the algebraic solver evaluates it like any other equation (property
    // calls, arrays, user functions...).
    std::vector<StmtPtr> extra;
    const auto& calls = problem_.calls;
    extra.push_back(makeEquation(makeVariable(problem_.integrationVar), makeVariable(calls[0].loVar)));
    for (const auto& c : calls) {
        extra.push_back(makeEquation(makeVariable(c.loVar), c.lowerExpr));
        extra.push_back(makeEquation(makeVariable(c.hiVar), c.upperExpr));
        if (c.stepExpr) extra.push_back(makeEquation(makeVariable(c.stepVar), c.stepExpr));
    }
    std::vector<std::string> driven;
    for (const auto& st : problem_.states) driven.push_back(st.name);
    auto sys = buildReduced(extra, driven);
    if (!sys->analysis->success) {
        error = "The algebraic equations of the model are not structurally solvable (" +
                sys->analysis->errorMessage + ")";
        return false;
    }

    // The limits and the step must not change along the trajectory: they cannot depend on the
    // integration variable or on a state (EES evaluates them once).
    std::set<std::string, CaseInsensitiveLess> seeds(driven.begin(), driven.end());
    seeds.insert(problem_.integrationVar);
    const auto dynamic = dependentVariables(*sys->ir, *sys->analysis, seeds);
    for (const auto& c : calls) {
        const struct { const std::string* var; const ExprPtr* expr; const char* what; } checks[] = {
            {&c.loVar, &c.lowerExpr, "lower limit"}, {&c.hiVar, &c.upperExpr, "upper limit"},
            {&c.stepVar, &c.stepExpr, "step"}};
        for (const auto& chk : checks) {
            if (chk.var->empty() || !dynamic.count(*chk.var)) continue;
            error = std::string("The ") + chk.what + " '" + astToString(*chk.expr) + "' of INTEGRAL() depends on the "
                    "integration variable '" + problem_.integrationVar + "' or on a state of the model: the "
                    "limits and the step must stay constant during the integration";
            return false;
        }
    }

    // Initial point: states at their guess (user states, as in the march) or 0 (auxiliary states).
    for (const auto& st : problem_.states) {
        double v = 0.0;
        if (!st.auxiliary) {
            const auto* info = fullIr_.getVariable(st.name);
            v = (info && info->guessValue) ? *info->guessValue : 1.0;
        }
        setDriverValue(*sys, st.name, v);
        sys->solver->setGuess(st.name, v);
    }
    SolveResult r = sys->solver->solve(solverOpts_, false);
    if (!r.success) {
        error = "Solving the model at the initial point to resolve the INTEGRAL limits failed: " +
                r.errorMessage;
        return false;
    }
    auto value = [&](const std::string& name, double& out) {
        auto it = r.variables.find(name);
        if (it == r.variables.end() || !std::isfinite(it->second)) return false;
        out = it->second;
        return true;
    };
    double lo0 = 0.0, hi0 = 0.0;
    if (!value(calls[0].loVar, lo0) || !value(calls[0].hiVar, hi0)) {
        error = "The limits of INTEGRAL() could not be evaluated (not finite)";
        return false;
    }
    t0 = lo0;
    tf = hi0;
    step = 0.0;
    bool stepTaken = false;
    for (const auto& c : calls) {
        double lo = 0.0, hi = 0.0;
        if (!value(c.loVar, lo) || !value(c.hiVar, hi)) {
            error = "The limits of INTEGRAL() could not be evaluated (not finite)";
            return false;
        }
        const double tolLo = 1e-9 * std::max(1.0, std::abs(lo0));
        const double tolHi = 1e-9 * std::max(1.0, std::abs(hi0));
        if (std::abs(lo - lo0) > tolLo || std::abs(hi - hi0) > tolHi) {
            error = "All INTEGRAL() calls must share the same [lo, hi] interval; found [" +
                    std::to_string(lo0) + ", " + std::to_string(hi0) + "] and [" +
                    std::to_string(lo) + ", " + std::to_string(hi) + "]";
            return false;
        }
        if (c.stepExpr) {
            double h = 0.0;
            if (!value(c.stepVar, h)) {
                error = "The step of INTEGRAL() could not be evaluated (not finite)";
                return false;
            }
            if (h < 0.0) {
                error = "INTEGRAL() step must be positive (0 = automatic step); got " + std::to_string(h);
                return false;
            }
            if (h > 0.0) {
                if (!stepTaken) { step = h; stepTaken = true; }
                else if (h != step)
                    diag.push(DiagnosticSeverity::Warning, "I002",
                              "INTEGRAL() calls use different steps (" + std::to_string(step) + " and " +
                              std::to_string(h) + "): the first explicit step is used for the whole integration",
                              "integral");
            }
        }
    }
    return true;
}

/// More rows than this over the integration interval is an output interval typed wrongly.
constexpr double kMaxTableRows = 1.0e6;

// ============================================================================
// Time-march loop
// ============================================================================

IntegralSolveResult IntegralSolver::solve(const IntegratorOptions& intOpt) {
    IntegralSolveResult result;
    result.problem = problem_;
    const bool fixedStepMethod =
        intOpt.method == IntegratorOptions::EulerExplicit ||
        intOpt.method == IntegratorOptions::EulerImplicit ||
        intOpt.method == IntegratorOptions::RK4;
    result.method = methodToString(intOpt.method) + (intOpt.richardson && fixedStepMethod ? " + Richardson" : "");
    if (!problem_.valid || problem_.states.empty()) {
        result.errorMessage = problem_.errorMessage.empty()
                                  ? std::string("invalid IntegralProblem")
                                  : problem_.errorMessage;
        return result;
    }
    if (!march_ || !march_->analysis || !march_->analysis->success) {
        result.errorMessage = "Algebraic subsystem analysis failed";
        if (march_ && march_->analysis && !march_->analysis->errorMessage.empty())
            result.errorMessage += ": " + march_->analysis->errorMessage;
        return result;
    }
    for (const auto& w : problem_.warnings)
        result.diagnostics.push(DiagnosticSeverity::Warning, "I002", w, "integral");

    // ---- Integration interval and step -------------------------------------
    // Constant limits and step (the common case) are folded at extraction; symbolic ones are
    // read from the model solved at the initial point.
    double t0 = problem_.lowerLimit;
    double tf = problem_.upperLimit;
    double modelStep = problem_.fixedStep;    // 5th INTEGRAL argument, 0 = automatic
    if (!problem_.limitsConstant || !problem_.stepConstant) {
        std::string error;
        if (!resolveSymbolicInterval(t0, tf, modelStep, error, result.diagnostics)) {
            result.errorMessage = error;
            return result;
        }
    }
    if (!(tf > t0)) {
        result.errorMessage = "Invalid integration interval [" + std::to_string(t0) + ", " +
                              std::to_string(tf) + "]: the upper limit must be greater than the lower limit";
        return result;
    }
    const double span = tf - t0;
    // Report the interval and step actually used (resolved when the limits are symbolic).
    result.problem.lowerLimit = t0;
    result.problem.upperLimit = tf;
    result.problem.fixedStep = modelStep;

    // Build the integrator (Richardson only for fixed-step methods).
    auto ig = wrapRichardson(createIntegrator(intOpt.method),
                             intOpt.richardson && fixedStepMethod);

    // Step size. An explicit step (5th INTEGRAL argument, else `integralFixedStep`) is honoured as
    // given: the number of steps follows from (hi - lo) / step and `integralMaxSteps` does not
    // apply. Without one, a fixed-step method takes exactly `integralMaxSteps` steps; RK45 adapts
    // and `integralMaxSteps` bounds the work (see below).
    const bool adaptive = (intOpt.method == IntegratorOptions::RK45);
    const int maxSteps = std::max(1, intOpt.maxSteps);
    double h = 0.0;
    long long plannedSteps = 0;   // fixed-step methods
    if (!adaptive) {
        const double explicitStep = modelStep > 0.0 ? modelStep : intOpt.fixedStep;
        if (explicitStep > 0.0) {
            const double n = span / explicitStep;
            if (!(n < 2.0e9)) {
                result.errorMessage = "The integration step " + std::to_string(explicitStep) +
                                      " would need more than 2e9 steps over [" + std::to_string(t0) +
                                      ", " + std::to_string(tf) + "]";
                return result;
            }
            plannedSteps = std::max<long long>(1, static_cast<long long>(std::ceil(n - 1e-9)));
            h = explicitStep;
        } else {
            plannedSteps = maxSteps;
            h = span / static_cast<double>(maxSteps);
        }
    } else {
        h = (modelStep > 0.0) ? modelStep : std::min(span / 10.0, span);
    }
    const double hMin = intOpt.minStep > 0.0 ? intOpt.minStep : span * 1e-9;
    const double hMax = intOpt.maxStep > 0.0 ? intOpt.maxStep : span;
    const double tolT = 1e-9 * span;   // what is left of the interval below this is not a step

    // ---- Set up the table --------------------------------------------------
    // Column names are matched to the model's variables without regard to case (EES names are
    // case-insensitive); the header uses the spelling of the model.
    auto canonical = [&](const std::string& name, bool& known) {
        const auto& vars = fullIr_.getVariables();
        auto it = vars.find(name);
        known = it != vars.end();
        return known ? it->first : name;
    };
    IntegralTableSpec spec = problem_.tableSpec.isPresent() ? problem_.tableSpec
                                                            : IntegralTableSpec{};
    std::vector<std::string> wanted;
    if (!spec.isPresent()) {
        spec.integrationVar = problem_.integrationVar;
        wanted.push_back(problem_.integrationVar);
        for (const auto& s : problem_.states)
            if (!s.tableName.empty()) wanted.push_back(s.tableName);
        for (const auto& d : problem_.derivativeNames)
            if (!problem_.isInternal(d)) wanted.push_back(d);
    } else if (spec.columns.empty()) {
        wanted.push_back(problem_.integrationVar);
    } else {
        wanted = spec.columns;
    }
    std::vector<std::string> columns;
    {
        std::set<std::string, CaseInsensitiveLess> seen;
        for (size_t i = 0; i < wanted.size(); ++i) {
            bool known = false;
            std::string name = canonical(wanted[i], known);
            if (i == 0 && !known) known = true;   // the integration variable column is always filled
            if (!seen.insert(name).second) {
                if (problem_.tableSpec.isPresent())
                    result.diagnostics.push(DiagnosticSeverity::Warning, "I001",
                                            "$IntegralTable: '" + wanted[i] + "' is listed twice; "
                                            "the second column is ignored", "integral");
                continue;
            }
            if (!known)
                result.diagnostics.push(DiagnosticSeverity::Warning, "I001",
                                        "$IntegralTable: '" + wanted[i] + "' is not a variable of the model; "
                                        "its column stays empty", "integral");
            columns.push_back(name);
        }
        if (problem_.tableSpec.isPresent() && !columns.empty() &&
            toLowerCopy(columns[0]) != toLowerCopy(problem_.integrationVar))
            result.diagnostics.push(DiagnosticSeverity::Warning, "I001",
                                    "$IntegralTable: the first name '" + columns[0] + "' is not the integration "
                                    "variable '" + problem_.integrationVar + "' of the INTEGRAL() calls; its column "
                                    "holds the integration variable", "integral");
    }
    result.table = IntegralTable(problem_.integrationVar);
    result.table.setColumns(columns);

    // Values of every tabulated column at (t, yCur): the algebraic subsystem is evaluated there so
    // states, derivatives and algebraic variables are consistent at this point.
    auto evalRow = [&](double t, const Eigen::VectorXd& yCur) {
        try {
            (void)evaluateRHS(t, yCur);
        } catch (...) {
            // Keep the algebraic state the failed re-solve reached.
        }
        // The integration variable and the states as the integrator holds them, the others from the state
        const auto& state = march_->solver->getEvaluator().getAllVariables();
        std::map<std::string, double> vals;
        for (const auto& col : result.table.columns()) {
            if (eesCompareFolded(col, problem_.integrationVar) == 0) {
                vals[col] = t;
                continue;
            }
            bool isState = false;
            for (int i = 0; i < static_cast<int>(problem_.states.size()) && !isState; ++i) {
                if (eesCompareFolded(col, problem_.states[i].name) == 0) {
                    vals[col] = yCur(i);
                    isState = true;
                }
            }
            if (isState) continue;
            auto it = state.find(col);
            if (it != state.end()) vals[col] = it->second;
        }
        return vals;
    };

    // ---- Initial state at t0 ----------------------------------------------
    Eigen::VectorXd y0(static_cast<int>(problem_.states.size()));
    // Provisional algebraic solve at t0: the user states keep their IR guesses, the auxiliary
    // states (integrals from t0 to t0) are 0.
    for (int i = 0; i < y0.size(); ++i) {
        const auto* info = fullIr_.getVariable(problem_.states[i].name);
        y0(i) = problem_.states[i].auxiliary ? 0.0
                                             : ((info && info->guessValue) ? *info->guessValue : 1.0);
    }
    SolveResult initRes;
    try {
        // Solve once at t0 to populate algebraic vars (e.g. y0 baseline).
        setDriverValue(*march_, problem_.integrationVar, t0);
        march_->solver->setGuess(problem_.integrationVar, t0);
        for (int i = 0; i < y0.size(); ++i) {
            setDriverValue(*march_, problem_.states[i].name, y0(i));
            march_->solver->setGuess(problem_.states[i].name, y0(i));
        }
        initRes = march_->solver->solve(solverOpts_, false);
    } catch (const std::exception& e) {
        result.errorMessage = std::string("Initial algebraic solve failed: ") + e.what();
        return result;
    }
    if (!initRes.success) {
        result.errorMessage = "Initial algebraic solve failed: " + initRes.errorMessage;
        return result;
    }
    // Output interval: `:n` of the directive (a number, or a variable holding a constant), else the
    // configuration key `integralOutputInterval`; 0 = a row at every step.
    double outInterval = solverOpts_.integralOutputInterval;
    std::string intervalSource = "integralOutputInterval";
    bool fromDirective = false;
    if (spec.hasOutputInterval) {
        outInterval = spec.outputInterval;
        intervalSource = "$IntegralTable: the output interval";
        fromDirective = true;
    } else if (!spec.outputIntervalVar.empty()) {
        auto it = initRes.variables.find(spec.outputIntervalVar);
        if (it != initRes.variables.end() && std::isfinite(it->second)) {
            outInterval = it->second;
            intervalSource = "$IntegralTable: the output interval '" + spec.outputIntervalVar + "'";
            fromDirective = true;
        } else {
            outInterval = 0.0;
            result.diagnostics.push(DiagnosticSeverity::Warning, "I001",
                                    "$IntegralTable: the output interval '" + spec.outputIntervalVar +
                                    "' is not a variable of the model with a finite value; "
                                    "a row is written at every step", "integral");
        }
    }
    if (!std::isfinite(outInterval) || outInterval < 0.0 || (fromDirective && outInterval == 0.0)) {
        result.errorMessage = intervalSource + " must be " + (fromDirective ? "positive" : "zero (every step) or positive") +
                              "; got " + std::to_string(outInterval);
        return result;
    }
    if (outInterval > 0.0 && span / outInterval > kMaxTableRows) {
        result.errorMessage = intervalSource + " " + std::to_string(outInterval) + " would write more than " +
                              std::to_string(static_cast<long long>(kMaxTableRows)) + " rows over [" +
                              std::to_string(t0) + ", " + std::to_string(tf) + "]";
        return result;
    }
    // Recover y(t0) from each state's base expression.
    try {
        for (int i = 0; i < y0.size(); ++i)
            y0(i) = evalBaseExpr(problem_.states[i].baseExpr, initRes.variables);
    } catch (const std::exception& e) {
        result.errorMessage = std::string("Initial value of a state: ") + e.what();
        return result;
    }
    std::optional<std::map<std::string, double>> rowAtT = evalRow(t0, y0);   // row at the start of the step, while known
    result.table.appendRow(t0, *rowAtT);

    // ---- March -------------------------------------------------------------
    Eigen::VectorXd y = y0;
    double t = t0;
    int taken = 0;
    int rejected = 0;
    // Rows at t0 + k*outInterval, independent of the steps (EES): a row between two steps is
    // interpolated linearly between the rows of the steps around it (`bracket`, two rows).
    long long nextRow = 1;
    IntegralTable bracket(problem_.integrationVar);
    bracket.setColumns(columns);
    RHSFunction rhs = [this](double tt, const Eigen::VectorXd& yy) {
        return evaluateRHS(tt, yy);
    };
    // Fixed-step methods take exactly `plannedSteps` steps. RK45 may reject steps: its budget is
    // 4 * integralMaxSteps attempted (accepted + rejected) steps.
    const long long maxAttempts = adaptive ? 4LL * maxSteps : plannedSteps;
    long long attempts = 0;

    try {
        while (t < tf - tolT && attempts < maxAttempts) {
            ++attempts;
            // Honour the Stop button promptly between steps (mirrors the
            // parametric study's per-point cancel check in server.cpp).
            if (solverOpts_.cancelToken &&
                solverOpts_.cancelToken->load(std::memory_order_relaxed)) {
                result.success = false;
                result.errorMessage = "Integration cancelled by user";
                result.totalSteps = taken;
                result.rejectedSteps = rejected;
                return result;
            }
            double hTry = h;
            if (t + hTry >= tf - tolT) hTry = tf - t;  // land exactly on tf
            const double tPrev = t;
            const Eigen::VectorXd yPrev = y;

            StepResult sr = ig->step(t, y, rhs, hTry, intOpt);

            if (!sr.accepted) {
                ++rejected;
                double hn = sr.nextStep > 0.0 ? sr.nextStep : hTry * 0.5;
                h = std::min(hMax, std::max(hMin, hn));
                if (h <= hMin) {
                    result.errorMessage = "Step size dropped below minimum (" +
                                          std::to_string(h) + ") at t=" + std::to_string(t);
                    result.totalSteps = taken;
                    result.rejectedSteps = rejected;
                    return result;
                }
                continue;
            }

            y = sr.yNew;
            t += sr.stepTaken;
            if (t >= tf - tolT) t = tf;   // what is left of the interval is rounding
            ++taken;
            result.acceptedStepSizes.push_back(sr.stepTaken);
            if (adaptive) {
                double hn = sr.nextStep > 0.0 ? sr.nextStep : hTry;
                h = std::min(hMax, std::max(hMin, hn));
            }

            // Rows of this step: at every step (no interval), else the output times it reached
            // (a time on the step's end is that row; the final point is always a row).
            const bool isLast = (t >= tf);
            bool onStep = isLast || outInterval <= 0.0;
            std::vector<double> between;   // output times strictly inside (tPrev, t)
            while (outInterval > 0.0) {
                const double tOut = t0 + static_cast<double>(nextRow) * outInterval;
                if (tOut >= tf - tolT || tOut > t + tolT) break;
                if (tOut >= t - tolT) onStep = true;
                else between.push_back(tOut);
                ++nextRow;
            }
            if (between.empty() && !onStep) {
                rowAtT.reset();   // no row at this step: the next step evaluates its own start
            } else if (between.empty()) {
                rowAtT = evalRow(t, y);
                result.table.appendRow(t, *rowAtT);
            } else {
                const auto rowPrev = rowAtT ? *rowAtT : evalRow(tPrev, yPrev);
                rowAtT = evalRow(t, y);
                bracket.clear();
                bracket.appendRow(tPrev, rowPrev);
                bracket.appendRow(t, *rowAtT);
                for (double tOut : between) {
                    std::map<std::string, double> vals;
                    for (const auto& col : columns) vals[col] = bracket.interpolate(col, tOut);
                    result.table.appendRow(tOut, vals);
                }
                if (onStep) result.table.appendRow(t, *rowAtT);
            }
        }
    } catch (const std::exception& e) {
        result.errorMessage = e.what();
        // An explicit method beyond its stability limit makes the states grow until an algebraic solve fails
        if (dynamic_cast<const AlgebraicSolveError*>(&e) && !adaptive && intOpt.method != IntegratorOptions::EulerImplicit) {
            std::ostringstream hint;
            hint << "\nThe step h = " << h << " may be too large for " << methodToString(intOpt.method)
                 << ": an explicit method with a fixed step is unstable above a limit set by the fastest dynamics of the "
                    "model, and the states then diverge. Reduce the step of INTEGRAL() (or integralFixedStep), or use "
                    "RK45 (adaptive step) or EulerImplicit";
            result.errorMessage += hint.str();
        }
        result.totalSteps = taken;
        result.rejectedSteps = rejected;
        // The table holds whatever was recorded before the failure.
        result.success = false;
        return result;
    }

    result.totalSteps = taken;
    result.rejectedSteps = rejected;
    if (t < tf) {
        // Never report a truncated trajectory as a success.
        result.success = false;
        result.errorMessage =
            "The integration stopped at " + problem_.integrationVar + " = " + std::to_string(t) +
            " before reaching the upper limit " + std::to_string(tf) + ": the budget of " +
            std::to_string(maxAttempts) + " attempted steps (4 x integralMaxSteps = 4 x " +
            std::to_string(maxSteps) + ") is exhausted. Increase integralMaxSteps, give a step in "
            "INTEGRAL(), or loosen integralRelTol / integralAbsTol";
        return result;
    }
    result.success = true;
    // Expose the final-step algebraic state (states + derivatives + algebraic
    // variables at t = tf) for downstream code (solution check, debug, .sol).
    // The internal variables of the INTEGRAL extraction (auxiliary states, integrand
    // definitions) are not part of the model: they are not reported.
    // (the last row evaluation solved the system at the upper limit: the state holds that point)
    result.algebraicResult.success = true;
    result.algebraicResult.status = SolverStatus::Success;
    result.algebraicResult.variables = march_->solver->getEvaluator().getAllVariables();
    result.algebraicResult.stringVariables = march_->solver->getEvaluator().getAllStringVariables();
    for (const auto& name : problem_.internalNames)
        result.algebraicResult.variables.erase(name);
    return result;
}

}  // namespace coolsolve
