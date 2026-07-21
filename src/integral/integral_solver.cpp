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
#include "coolsolve/integral/integral_solver.h"

#include "coolsolve/variable_inference.h"

#include <algorithm>
#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <stdexcept>
#include <string>

namespace coolsolve {

namespace {

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

}  // namespace

// ============================================================================
// Construction — build the algebraic subsystem
// ============================================================================

IntegralSolver::IntegralSolver(const Program& program, const IR& ir,
                               const StructuralAnalysisResult& analysis,
                               const SolverOptions& options)
    : fullIr_(ir), solverOpts_(options) {
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

    // Harvest the $IntegralTable spec (if any) straight off the AST.
    for (const auto& stmt : program.statements) {
        if (stmt->is<Directive>()) {
            const auto& d = stmt->as<Directive>();
            if (d.hasIntegralTableSpec) {
                problem_.tableSpec = d.integralTableSpec;
                break;
            }
        }
    }

    if (!problem_.hasIntegral) {
        problem_.valid = false;
        problem_.errorMessage = "no INTEGRAL() equations found";
        return;
    }

    // 1. Build a filtered Program that drops the integral-declaring equations,
    //    and add "driver" equations (v = <literal>) for each integrator-owned
    //    variable (states + integration var). The drivers make the reduced IR
    //    square (otherwise the free states/time fail the analyzer's squareness
    //    gate) and let the integrator pin these values each step by mutating
    //    the shared NumberLiteral.
    Program filtered;
    filtered.sourceFilename = program.sourceFilename;
    for (const auto& stmt : program.statements) {
        if (stmt->is<Equation>()) {
            const auto& eq = stmt->as<Equation>();
            if (containsIntegralCall(eq.rhs) || containsIntegralCall(eq.lhs)) continue;
        }
        filtered.statements.push_back(stmt);
    }
    auto addDriver = [&](const std::string& varName) {
        ExprPtr lit = makeNumber(0.0);
        driverLiterals_[varName] = lit;
        filtered.statements.push_back(makeEquation(makeVariable(varName), lit));
    };
    for (const auto& s : problem_.states) addDriver(s.name);
    addDriver(problem_.integrationVar);

    // 2. Build the reduced IR, run inference/initialisation on it, then copy
    //    the user's guesses / loaded initials from the full IR.
    reducedIr_ = std::make_unique<IR>(IR::fromAST(filtered));
    inferVariables(*reducedIr_);
    initializeVariables(*reducedIr_);
    for (const auto& [name, info] : reducedIr_->getVariables()) {
        (void)info;
        if (const auto* full = ir.getVariable(name)) {
            auto* mut = reducedIr_->getVariableMutable(name);
            if (mut) {
                mut->guessValue = full->guessValue;
                mut->solutionValue = full->solutionValue;
                mut->lowerBound = full->lowerBound;
                mut->upperBound = full->upperBound;
            }
        }
    }

    // 3. Structural analysis on the reduced (algebraic) system.
    reducedAnalysis_ = std::make_unique<StructuralAnalysisResult>(
        StructuralAnalyzer::analyze(*reducedIr_));

    // 4. The algebraic solver, reused at every time step.
    algebraicSolver_ = std::make_unique<Solver>(*reducedIr_, *reducedAnalysis_,
                                                options.coolpropConfig);
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
        auto it = vars.find(v.name);
        if (it != vars.end()) return it->second;
        // Fall back to the full IR's stored value (guesses / initials).
        const auto* info = fullIr_.getVariable(v.name);
        if (info) {
            if (info->solutionValue) return *info->solutionValue;
            if (info->guessValue) return *info->guessValue;
        }
        return 1.0;
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
        std::string n = expr->as<FunctionCall>().name;
        std::transform(n.begin(), n.end(), n.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (n == "pi") return M_PI;
        if (n == "abs") return std::abs(evalBaseExpr(expr->as<FunctionCall>().args.at(0), vars));
    }
    return 0.0;
}

Eigen::VectorXd IntegralSolver::evaluateRHS(double t, const Eigen::VectorXd& y) {
    // Pin the integrator-owned variables: update the driver equations' literals
    // (so the explicit driver blocks resolve to these values) and set them as
    // algebraic guesses so any block reading them as an external sees the same.
    setDriverValue(problem_.integrationVar, t);
    algebraicSolver_->setGuess(problem_.integrationVar, t);
    for (int i = 0; i < static_cast<int>(problem_.states.size()); ++i) {
        setDriverValue(problem_.states[i].name, y(i));
        algebraicSolver_->setGuess(problem_.states[i].name, y(i));
    }

    SolveResult r = algebraicSolver_->solve(solverOpts_, false);
    lastAlgebraic_ = r;  // refreshed for recordRow / readDerivative
    if (!r.success) {
        throw std::runtime_error("Algebraic solve failed at t=" + std::to_string(t) +
                                 ": " + r.errorMessage);
    }

    // f_i = derivative of state i, read from the algebraic solution.
    Eigen::VectorXd f(static_cast<int>(problem_.states.size()));
    for (int i = 0; i < f.size(); ++i) {
        const auto& sv = problem_.states[i];
        if (!sv.integrandVar.empty()) {
            auto it = r.variables.find(sv.integrandVar);
            f(i) = (it != r.variables.end()) ? it->second : 0.0;
        } else {
            f(i) = evalBaseExpr(sv.integrandExpr, r.variables);
        }
    }
    return f;
}

double IntegralSolver::readDerivative(int /*stateIndex*/) const {
    // Reserved for future streaming use; evaluation happens in evaluateRHS.
    return 0.0;
}

void IntegralSolver::setDriverValue(const std::string& name, double value) {
    auto it = driverLiterals_.find(name);
    if (it != driverLiterals_.end() && it->second && it->second->is<NumberLiteral>())
        it->second->as<NumberLiteral>().value = value;
}

// ============================================================================
// Time-march loop
// ============================================================================

IntegralSolveResult IntegralSolver::solve(const IntegratorOptions& intOpt) {
    IntegralSolveResult result;
    result.problem = problem_;
    if (!problem_.valid || problem_.states.empty()) {
        result.errorMessage = problem_.errorMessage.empty()
                                  ? std::string("invalid IntegralProblem")
                                  : problem_.errorMessage;
        return result;
    }
    if (!reducedAnalysis_ || !reducedAnalysis_->success) {
        result.errorMessage = "Algebraic subsystem analysis failed";
        return result;
    }

    // ---- Integration interval and initial state ----------------------------
    if (!problem_.limitsConstant) {
        result.errorMessage = "Non-constant integration limits are not yet supported "
                              "(resolve parameters before the INTEGRAL call).";
        return result;
    }
    const double t0 = problem_.lowerLimit;
    const double tf = problem_.upperLimit;
    if (tf <= t0) {
        result.errorMessage = "Invalid integration interval [" + std::to_string(t0) + ", " +
                              std::to_string(tf) + "]";
        return result;
    }

    // Build the integrator (Richardson only for fixed-step methods).
    const bool fixedStepMethod =
        intOpt.method == IntegratorOptions::EulerExplicit ||
        intOpt.method == IntegratorOptions::EulerImplicit ||
        intOpt.method == IntegratorOptions::RK4;
    auto ig = wrapRichardson(createIntegrator(intOpt.method),
                             intOpt.richardson && fixedStepMethod);

    // Determine step size.
    double h = 0.0;
    const bool adaptive = (intOpt.method == IntegratorOptions::RK45);
    if (!adaptive) {
        if (problem_.fixedStep > 0.0) h = problem_.fixedStep;
        else if (intOpt.fixedStep > 0.0) h = intOpt.fixedStep;
        else h = (tf - t0) / std::max(1, intOpt.maxSteps);
    } else {
        h = (problem_.fixedStep > 0.0) ? problem_.fixedStep
                                       : std::min((tf - t0) / 10.0, (tf - t0));
    }
    const double hMin =
        intOpt.minStep > 0.0 ? intOpt.minStep : (tf - t0) * 1e-9;
    const double hMax =
        intOpt.maxStep > 0.0 ? intOpt.maxStep : (tf - t0);

    // ---- Set up the table --------------------------------------------------
    IntegralTableSpec spec = problem_.tableSpec.isPresent() ? problem_.tableSpec
                                                            : IntegralTableSpec{};
    if (!spec.isPresent()) {
        spec.integrationVar = problem_.integrationVar;
        spec.columns.push_back(problem_.integrationVar);
        for (const auto& s : problem_.states) spec.columns.push_back(s.name);
        for (const auto& d : problem_.derivativeNames) spec.columns.push_back(d);
    } else if (spec.columns.empty()) {
        spec.columns.push_back(problem_.integrationVar);
    }
    result.table = IntegralTable(problem_.integrationVar);
    result.table.setColumns(spec.columns);
    const double outInterval = spec.outputInterval;

    auto recordRow = [&](double t) {
        // Evaluate the algebraic subsystem at (t, y) so every tabulated column
        // (states, derivatives, algebraic vars) is consistent at this point.
        const Eigen::VectorXd yCur = currentY_;  // set by the march loop below
        try {
            (void)evaluateRHS(t, yCur);
        } catch (...) {
            // Keep the previous algebraic state if the re-solve fails.
        }
        SolveResult last = lastAlgebraic_;  // refreshed by evaluateRHS
        std::map<std::string, double> vals;
        for (const auto& [n, v] : last.variables) vals[n] = v;
        for (int i = 0; i < static_cast<int>(problem_.states.size()); ++i)
            vals[problem_.states[i].name] = yCur(i);
        vals[problem_.integrationVar] = t;
        result.table.appendRow(t, vals);
    };

    // ---- Initial state at t0 ----------------------------------------------
    Eigen::VectorXd y0(static_cast<int>(problem_.states.size()));
    // Provisional algebraic solve at t0 (states keep their IR guesses).
    for (int i = 0; i < y0.size(); ++i) {
        const auto* info = fullIr_.getVariable(problem_.states[i].name);
        y0(i) = (info && info->guessValue) ? *info->guessValue : 1.0;
    }
    SolveResult initRes;
    try {
        // Solve once at t0 to populate algebraic vars (e.g. y0 baseline).
        setDriverValue(problem_.integrationVar, t0);
        algebraicSolver_->setGuess(problem_.integrationVar, t0);
        for (int i = 0; i < y0.size(); ++i) {
            setDriverValue(problem_.states[i].name, y0(i));
            algebraicSolver_->setGuess(problem_.states[i].name, y0(i));
        }
        initRes = algebraicSolver_->solve(solverOpts_, false);
    } catch (const std::exception& e) {
        result.errorMessage = std::string("Initial algebraic solve failed: ") + e.what();
        return result;
    }
    if (!initRes.success) {
        result.errorMessage = "Initial algebraic solve failed: " + initRes.errorMessage;
        return result;
    }
    // Recover y(t0) from each state's base expression.
    for (int i = 0; i < y0.size(); ++i)
        y0(i) = evalBaseExpr(problem_.states[i].baseExpr, initRes.variables);
    currentY_ = y0;
    lastAlgebraic_ = initRes;
    recordRow(t0);

    // ---- March -------------------------------------------------------------
    Eigen::VectorXd y = y0;
    double t = t0;
    int taken = 0;
    int rejected = 0;
    double nextRecordT = (outInterval > 0.0) ? t0 + outInterval : tf;
    const int maxSteps = std::max(1, intOpt.maxSteps);
    RHSFunction rhs = [this](double tt, const Eigen::VectorXd& yy) {
        return evaluateRHS(tt, yy);
    };

    try {
        for (int step = 0; step < maxSteps * 4 && t < tf - 1e-12; ++step) {
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
            if (t + hTry > tf) hTry = tf - t;  // land exactly on tf

            StepResult sr = ig->step(t, y, rhs, hTry, intOpt);
            lastAlgebraic_ = {};  // recordRow will refresh via evaluateRHS

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
            currentY_ = y;
            ++taken;
            result.acceptedStepSizes.push_back(sr.stepTaken);
            if (!adaptive) {
                h = hTry;  // fixed step
            } else {
                double hn = sr.nextStep > 0.0 ? sr.nextStep : hTry;
                h = std::min(hMax, std::max(hMin, hn));
            }

            // Record honouring the output interval (always record the final point).
            bool isLast = (t >= tf - 1e-9);
            if (outInterval <= 0.0 || isLast || t >= nextRecordT - 1e-9) {
                recordRow(t);
                if (!isLast && outInterval > 0.0)
                    nextRecordT = std::min(tf, t + outInterval);
            }
        }
    } catch (const std::exception& e) {
        result.errorMessage = e.what();
        result.totalSteps = taken;
        result.rejectedSteps = rejected;
        // The table holds whatever was recorded before the failure.
        result.success = false;
        return result;
    }

    result.totalSteps = taken;
    result.rejectedSteps = rejected;
    result.success = true;
    // Expose the final-step algebraic state (states + derivatives + algebraic
    // variables at t = tf) for downstream code (solution check, debug, .sol).
    result.algebraicResult = lastAlgebraic_;
    return result;
}

}  // namespace coolsolve
