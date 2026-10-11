/**
 * @file integral_extraction.cpp
 * @brief Extract an `IntegralProblem` from a parsed IR.
 *
 * Walks every equation looking for `INTEGRAL(integrand, var, lo, hi[, step])`
 * calls (anywhere in an equation: an INTEGRAL call is a *value*, the integral from
 * `lo` to the current value of `var`), classifies state vs algebraic variables,
 * validates structural preconditions, and runs a conservative high-index
 * rejection.  See `docs/integral_table.md` §1 (Mathematical model) and §3
 * (Architecture).
 *
 * Two kinds of state, one per call:
 *  - the plain additive form `y = base + INTEGRAL(...)` (base free of INTEGRAL and
 *    constant along the trajectory) makes `y` a state owned by the integrator;
 *  - any other form (`y = m*INTEGRAL(...)`, `-INTEGRAL(...)`, `INTEGRAL(...)/3.6E6`,
 *    two calls in one equation, an INTEGRAL on the left-hand side...) gets an
 *    *auxiliary* state per call (internal name, I(lo) = 0, dI/dt = integrand) and the
 *    equation, with the call replaced by that state, stays an ordinary algebraic
 *    equation.
 *
 * The initial state is *not* extracted here: at t = t0 the integral term has a
 * zero-width interval and evaluates to 0, so `y(t0)` falls out of the algebraic
 * solve at the first step (handled by `IntegralSolver`).
 */
#include "coolsolve/ees_string.h"
#include "coolsolve/integral/integral_problem.h"

#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <set>
#include <sstream>

namespace coolsolve {

namespace {

// ---- small helpers ---------------------------------------------------------

// Names are compared as EES does: without regard to case, accented letters included
std::string toLower(const std::string& s) { return eesFoldCase(s); }

/// Find the first `INTEGRAL(...)` FunctionCall in `expr` (depth-first).
/// Returns nullptr if none.
const FunctionCall* findIntegralCall(const ExprPtr& expr) {
    if (!expr) return nullptr;
    const FunctionCall* found = nullptr;
    std::visit([&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, FunctionCall>) {
            if (toLower(node.name) == "integral") found = &node;
        }
    }, expr->node);
    if (found) return found;
    std::visit([&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, UnaryOp>) {
            found = findIntegralCall(node.operand);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            found = findIntegralCall(node.left);
            if (!found) found = findIntegralCall(node.right);
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            for (const auto& a : node.args)
                if ((found = findIntegralCall(a))) return;
            for (const auto& na : node.namedArgs)
                if ((found = findIntegralCall(na.second))) return;
        } else if constexpr (std::is_same_v<T, Variable>) {
            for (const auto& idx : node.indices)
                if ((found = findIntegralCall(idx))) return;
        }
    }, expr->node);
    return found;
}

/// Collect the outermost INTEGRAL calls of `expr`, in depth-first (left to right) order.
/// A call nested in the arguments of another call is not collected (the extraction
/// rejects it separately).
void collectIntegralCalls(const ExprPtr& expr, std::vector<const FunctionCall*>& out) {
    if (!expr) return;
    if (expr->is<FunctionCall>() && toLower(expr->as<FunctionCall>().name) == "integral") {
        out.push_back(&expr->as<FunctionCall>());
        return;
    }
    std::visit([&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, UnaryOp>) {
            collectIntegralCalls(node.operand, out);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            collectIntegralCalls(node.left, out);
            collectIntegralCalls(node.right, out);
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            for (const auto& a : node.args) collectIntegralCalls(a, out);
            for (const auto& na : node.namedArgs) collectIntegralCalls(na.second, out);
        } else if constexpr (std::is_same_v<T, Variable>) {
            for (const auto& idx : node.indices) collectIntegralCalls(idx, out);
        }
    }, expr->node);
}

/// Names of the variables used by `expr` (array elements by their flattened name).
void collectVariableNames(const ExprPtr& expr, std::set<std::string, CaseInsensitiveLess>& out) {
    if (!expr) return;
    std::visit([&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Variable>) {
            out.insert(node.flattenedName());
            for (const auto& idx : node.indices) collectVariableNames(idx, out);
        } else if constexpr (std::is_same_v<T, UnaryOp>) {
            collectVariableNames(node.operand, out);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            collectVariableNames(node.left, out);
            collectVariableNames(node.right, out);
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            for (const auto& a : node.args) collectVariableNames(a, out);
            for (const auto& na : node.namedArgs) collectVariableNames(na.second, out);
        }
    }, expr->node);
}

/// Is `target` reached from `e` only through additions and subtractions, with a
/// net + sign (`a + INTEGRAL`, `INTEGRAL - b + c`, `a - (b - INTEGRAL)`)?  Anything else
/// (a factor, a division, `a - INTEGRAL`, a function argument...) is not the plain
/// additive state form.
bool reachedAdditively(const ExprPtr& e, const FunctionCall* target, int sign) {
    if (!e) return false;
    if (e->is<FunctionCall>()) return &e->as<FunctionCall>() == target && sign > 0;
    if (e->is<UnaryOp>()) {
        const auto& u = e->as<UnaryOp>();
        if (u.op == "-") return reachedAdditively(u.operand, target, -sign);
        if (u.op == "+") return reachedAdditively(u.operand, target, sign);
        return false;
    }
    if (e->is<BinaryOp>()) {
        const auto& b = e->as<BinaryOp>();
        if (b.op != "+" && b.op != "-") return false;
        std::vector<const FunctionCall*> inLeft;
        collectIntegralCalls(b.left, inLeft);
        if (!inLeft.empty()) return reachedAdditively(b.left, target, sign);
        return reachedAdditively(b.right, target, (b.op == "-") ? -sign : sign);
    }
    return false;
}

/// Clone an expression, applying `replace` to every node; if `replace` returns
/// non-null the cloned subtree is substituted (used to blank INTEGRAL → 0).
ExprPtr cloneExpr(const ExprPtr& e, const std::function<ExprPtr(const ExprPtr&)>& replace) {
    if (!e) return nullptr;
    if (auto r = replace(e)) return r;
    auto out = std::make_shared<Expression>(*e);
    std::visit([&](auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, UnaryOp>) {
            node.operand = cloneExpr(node.operand, replace);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            node.left = cloneExpr(node.left, replace);
            node.right = cloneExpr(node.right, replace);
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            for (auto& a : node.args) a = cloneExpr(a, replace);
            for (auto& na : node.namedArgs) na.second = cloneExpr(na.second, replace);
        } else if constexpr (std::is_same_v<T, Variable>) {
            for (auto& idx : node.indices) idx = cloneExpr(idx, replace);
        }
    }, out->node);
    return out;
}

/// Try to fold `expr` to a double. Recognises numbers, +/-/*/, ^, unary -,
/// and the `pi()` constant. Returns false on any Variable / StringLiteral.
bool evalConstantExpr(const ExprPtr& e, double& out) {
    if (!e) return false;
    if (e->is<NumberLiteral>()) { out = e->as<NumberLiteral>().value; return true; }
    if (e->is<UnaryOp>()) {
        const auto& u = e->as<UnaryOp>();
        double v;
        if (!evalConstantExpr(u.operand, v)) return false;
        out = (u.op == "-") ? -v : v;
        return true;
    }
    if (e->is<BinaryOp>()) {
        const auto& b = e->as<BinaryOp>();
        double l, r;
        if (!evalConstantExpr(b.left, l) || !evalConstantExpr(b.right, r)) return false;
        if (b.op == "+") { out = l + r; return true; }
        if (b.op == "-") { out = l - r; return true; }
        if (b.op == "*") { out = l * r; return true; }
        if (b.op == "/") { if (r == 0.0) return false; out = l / r; return true; }
        if (b.op == "^") { out = std::pow(l, r); return true; }
        return false;
    }
    if (e->is<FunctionCall>()) {
        if (toLower(e->as<FunctionCall>().name) == "pi") { out = M_PI; return true; }
        return false;
    }
    return false;  // Variable, StringLiteral
}

}  // namespace

// ============================================================================
// Public API
// ============================================================================

bool hasIntegral(const IR& ir) {
    for (const auto& eq : ir.getEquations())
        if (findIntegralCall(eq.rhs) || findIntegralCall(eq.lhs)) return true;
    return false;
}

IntegralProblem extractIntegralProblem(const IR& ir,
                                       const StructuralAnalysisResult& /*analysis*/,
                                       const std::set<int>* forceAuxiliary) {
    IntegralProblem prob;

    // Pull the `$IntegralTable` spec off the parsed AST, if present, by scanning
    // the IR equations' source is not enough — the spec lives on the Directive
    // statements in the Program.  The IR does not retain directives, so the
    // table spec is wired in by the runner (Phase 6) before solving.  Here we
    // only leave `tableSpec` default-constructed; extraction is structural.

    auto fail = [&](const std::string& msg) {
        prob.valid = false;
        prob.errorMessage = msg;
        return prob;
    };

    // -- 1. Walk equations, harvest every INTEGRAL call. --------------------
    struct EquationCalls {
        const EquationInfo* eq;
        std::vector<const FunctionCall*> calls;
    };
    std::vector<EquationCalls> work;
    std::set<std::string, CaseInsensitiveLess> integralLhsNames;  // LHS variables of those equations
    for (const auto& eq : ir.getEquations()) {
        EquationCalls ec{&eq, {}};
        collectIntegralCalls(eq.lhs, ec.calls);
        collectIntegralCalls(eq.rhs, ec.calls);
        if (ec.calls.empty()) continue;
        prob.hasIntegral = true;
        if (eq.lhs && eq.lhs->is<Variable>()) integralLhsNames.insert(eq.lhs->as<Variable>().flattenedName());
        work.push_back(std::move(ec));
    }

    // Internal variable names that cannot clash with the model's variables.
    std::set<std::string, CaseInsensitiveLess> usedInternal;
    auto internalName = [&](const std::string& stem, int k) {
        std::string name = stem + std::to_string(k);
        while (ir.getVariable(name) || usedInternal.count(name)) name += "_";
        usedInternal.insert(name);
        prob.internalNames.push_back(name);
        return name;
    };

    // -- 2. Validate each call and create its state. --------------------------
    int callIndex = 0;
    bool allLimitsConstant = true;
    bool allStepsConstant = true;
    bool firstExplicitStepTaken = false;
    double firstConstLo = 0.0, firstConstHi = 0.0;
    bool firstLimitsConstant = false;

    for (const auto& ec : work) {
        const EquationInfo& eq = *ec.eq;
        const std::string where = " (line " + std::to_string(eq.sourceLine) + ")";

        // Pass 1: validate the arguments of every call of this equation.
        for (const FunctionCall* ic : ec.calls) {
            const auto& args = ic->args;
            if (args.size() != 4 && args.size() != 5) {
                return fail("INTEGRAL() expects 4 or 5 arguments (integrand, var, lo, hi[, step]); "
                            "got " + std::to_string(args.size()) + where + ".");
            }
            if (!args[1] || !args[1]->is<Variable>() || !args[1]->as<Variable>().indices.empty()) {
                return fail("INTEGRAL() integration variable must be a scalar variable" + where + ".");
            }
            for (const auto& a : args) {
                if (findIntegralCall(a)) {
                    return fail("Nested INTEGRAL() calls are not supported in this version "
                                "(multi-variable integration). See docs/integral_table.md §7.1.");
                }
            }
        }

        // Is this the plain additive state form  y = base + INTEGRAL(...)  with a base that is
        // constant along the trajectory?  (one call, scalar variable on the left, the call reached
        // by additions only, base free of the integration variable and of any variable defined
        // through an INTEGRAL)
        bool userState = false;
        ExprPtr baseExpr;
        std::vector<std::string> baseVariableList;
        if (ec.calls.size() == 1 && eq.lhs && eq.lhs->is<Variable>() &&
            eq.lhs->as<Variable>().indices.empty() && eq.rhs &&
            reachedAdditively(eq.rhs, ec.calls[0], +1)) {
            baseExpr = cloneExpr(eq.rhs, [&](const ExprPtr& e) -> ExprPtr {
                if (e && e->is<FunctionCall>() && &e->as<FunctionCall>() == ec.calls[0])
                    return makeNumber(0.0, e->sourceLineNumber);
                return nullptr;
            });
            std::set<std::string, CaseInsensitiveLess> baseVars;
            collectVariableNames(baseExpr, baseVars);
            const std::string integVarHere = ec.calls[0]->args[1]->as<Variable>().name;
            userState = !baseVars.count(integVarHere);
            for (const auto& v : baseVars)
                if (integralLhsNames.count(v)) userState = false;   // base follows another state (or y itself)
            if (forceAuxiliary && forceAuxiliary->count(eq.id)) userState = false;
            if (userState) baseVariableList.assign(baseVars.begin(), baseVars.end());
        }

        // Pass 2: one state per call.
        std::vector<ExprPtr> replacements;    // auxiliary variable of each call, in order
        for (const FunctionCall* ic : ec.calls) {
            const auto& args = ic->args;
            const std::string integVar = args[1]->as<Variable>().name;

            StateVariable sv;
            sv.equationId = eq.id;
            sv.integrandExpr = args[0];
            if (args[0] && args[0]->is<Variable>()) {
                sv.integrandVar = args[0]->as<Variable>().flattenedName();
            } else {
                // Any other integrand expression is defined by an extra algebraic equation, so it is
                // evaluated by the algebraic solver (arrays, property calls, user functions...).
                sv.integrandVar = internalName("__integrand_", callIndex + 1);
                prob.integrandEquations.emplace_back(sv.integrandVar, args[0]);
            }
            if (userState) {
                sv.name = eq.lhs->as<Variable>().name;
                sv.baseExpr = baseExpr;
                sv.baseVariables = baseVariableList;
                sv.tableName = sv.name;
            } else {
                sv.name = internalName("__integral_", callIndex + 1);
                sv.baseExpr = makeNumber(0.0);
                sv.auxiliary = true;
                if (ec.calls.size() == 1 && eq.lhs && eq.lhs->is<Variable>())
                    sv.tableName = eq.lhs->as<Variable>().flattenedName();
                replacements.push_back(makeVariable(sv.name, eq.sourceLine));
            }

            // Limits / step: prefer constant folding, keep the expressions for the solver.
            double lo = 0.0, hi = 0.0, step = 0.0;
            const bool loOk = evalConstantExpr(args[2], lo);
            const bool hiOk = evalConstantExpr(args[3], hi);
            const bool hasStep = args.size() == 5;
            const bool stepOk = hasStep && evalConstantExpr(args[4], step);
            if (hasStep && stepOk && step < 0.0)
                return fail("INTEGRAL() step must be positive (0 = automatic step)" + where + ".");

            IntegralCall call;
            call.equationId = eq.id;
            call.stateIndex = static_cast<int>(prob.states.size());
            call.lowerExpr = args[2];
            call.upperExpr = args[3];
            if (hasStep) call.stepExpr = args[4];

            if (prob.calls.empty()) {
                prob.integrationVar = integVar;
                prob.lowerLimitExpr = args[2];
                prob.upperLimitExpr = args[3];
                prob.lowerLimit = lo;
                prob.upperLimit = hi;
                firstLimitsConstant = loOk && hiOk;
                firstConstLo = lo;
                firstConstHi = hi;
            } else {
                if (toLower(integVar) != toLower(prob.integrationVar)) {
                    return fail("All INTEGRAL() calls must share the same integration variable; "
                                "found '" + integVar + "' and '" + prob.integrationVar + "'.");
                }
                // Constant intervals are compared here; symbolic ones by the solver once the
                // limits are known (see IntegralSolver::solve).
                if (firstLimitsConstant && loOk && hiOk && (lo != firstConstLo || hi != firstConstHi)) {
                    return fail("All INTEGRAL() calls must share the same [lo, hi] interval.");
                }
            }
            if (!(loOk && hiOk)) allLimitsConstant = false;
            if (hasStep && !stepOk) allStepsConstant = false;
            // The first explicit step wins; a step of 0 means "automatic".
            if (hasStep && stepOk && step > 0.0) {
                if (!firstExplicitStepTaken) {
                    prob.fixedStep = step;
                    firstExplicitStepTaken = true;
                } else if (step != prob.fixedStep) {
                    prob.warnings.push_back(
                        "INTEGRAL() calls use different steps (" + std::to_string(prob.fixedStep) +
                        " and " + std::to_string(step) + "): the first explicit step is used for the whole "
                        "integration" + where + ".");
                }
            }

            prob.calls.push_back(call);
            prob.states.push_back(sv);
            ++callIndex;
        }

        if (userState) {
            prob.integralEquationIds.push_back(eq.id);
        } else {
            RewrittenEquation rw;
            rw.equationId = eq.id;
            rw.originalLhs = eq.lhs;
            rw.originalRhs = eq.rhs;
            size_t next = 0;
            auto replaceCall = [&](const ExprPtr& e) -> ExprPtr {
                if (e && e->is<FunctionCall>() && toLower(e->as<FunctionCall>().name) == "integral")
                    return replacements.at(next++);
                return nullptr;
            };
            // Calls are numbered in depth-first order, left side first, as collected above.
            rw.lhs = cloneExpr(eq.lhs, replaceCall);
            rw.rhs = cloneExpr(eq.rhs, replaceCall);
            prob.rewrittenEquations.push_back(std::move(rw));
            prob.auxEquationIds.push_back(eq.id);
        }
    }

    if (!prob.hasIntegral) {
        prob.valid = false;
        prob.errorMessage = "no INTEGRAL() equations found";
        return prob;
    }
    // The integration variable belongs to the integrator: no INTEGRAL equation may define it
    // (names are case-insensitive: `T = ... INTEGRAL(f, t, ...)` is the same variable).
    for (const auto& ec : work) {
        if (ec.eq->lhs && ec.eq->lhs->is<Variable>() &&
            toLower(ec.eq->lhs->as<Variable>().flattenedName()) == toLower(prob.integrationVar)) {
            return fail("An equation with an INTEGRAL() call defines '" + ec.eq->lhs->as<Variable>().name +
                        "', which is the integration variable '" + prob.integrationVar +
                        "' (variable names are not case-sensitive) (line " +
                        std::to_string(ec.eq->sourceLine) + ").");
        }
    }
    prob.limitsConstant = allLimitsConstant;
    prob.stepConstant = allStepsConstant;
    if (!allLimitsConstant || !allStepsConstant) {
        // The solver resolves symbolic limits/steps from the initial solve; it defines them as
        // internal variables `__int_lo_k`, `__int_hi_k`, `__int_step_k`.
        for (size_t k = 0; k < prob.calls.size(); ++k) {
            auto& c = prob.calls[k];
            c.loVar = internalName("__int_lo_", static_cast<int>(k) + 1);
            c.hiVar = internalName("__int_hi_", static_cast<int>(k) + 1);
            if (c.stepExpr) c.stepVar = internalName("__int_step_", static_cast<int>(k) + 1);
        }
    }

    // Convenience name lists (user states and the derivative variables).
    for (const auto& st : prob.states) {
        prob.stateNames.push_back(st.name);
        if (!st.integrandVar.empty()) prob.derivativeNames.push_back(st.integrandVar);
    }

    // -- 3. (nested integrals were rejected above) ---------------------------

    // -- 4. Classify variables: state / algebraic / integration var. ---------
    // (auxiliary states are internal variables, not IR variables: only the user states count here)
    std::set<std::string, CaseInsensitiveLess> stateSet;
    for (const auto& st : prob.states)
        if (!st.auxiliary) stateSet.insert(st.name);
    for (const auto& [name, info] : ir.getVariables()) {
        if (stateSet.count(name)) continue;          // state variable (integrator-owned)
        if (toLower(name) == toLower(prob.integrationVar)) continue;   // independent variable (integrator-owned)
        prob.algebraicVars.push_back(name);
    }

    // Algebraic equations = every equation that is not the declaration of a user state
    // (equations with auxiliary states are ordinary algebraic equations).
    std::set<int> integralSet(prob.integralEquationIds.begin(),
                              prob.integralEquationIds.end());
    for (const auto& eq : ir.getEquations())
        if (!integralSet.count(eq.id)) prob.algebraicEquationIds.push_back(eq.id);

    // -- 5. Conservative high-index detection (placeholder for Pantelides). --
    // A state variable may legitimately appear in its own integral equation and
    // in the equation that defines its derivative. If it appears in *another*
    // algebraic equation that does not reference its derivative, the model is
    // likely higher-than-index-1 and we cannot solve it without index reduction.
    std::set<std::string, CaseInsensitiveLess> derivSet(prob.derivativeNames.begin(),
                                                        prob.derivativeNames.end());
    for (const auto& eq : ir.getEquations()) {
        if (integralSet.count(eq.id)) continue;
        // Does this algebraic equation define a derivative? (its output or any
        // var is a derivative name) -> fine.
        bool definesDeriv = false;
        for (const auto& v : eq.variables)
            if (derivSet.count(v)) { definesDeriv = true; break; }
        if (definesDeriv) continue;
        // Does it constrain a state variable directly?
        for (const auto& v : eq.variables) {
            if (stateSet.count(v)) {
                prob.diagnostics.push_back(
                    "Potential high-index DAE: state variable '" + v +
                    "' is constrained by an algebraic equation (line " +
                    std::to_string(eq.sourceLine) +
                    ") that does not reference its derivative. "
                    "Index reduction is not yet supported; the solve may fail.");
                break;
            }
        }
    }

    // -- 6. Structural squareness note (algebraic subsystem). ---------------
    // After removing the K integral equations and K state variables, the
    // algebraic subsystem should be square. The integration variable (held by
    // the integrator) accounts for the +1 "equation" t = t_current.
    const long kStates = static_cast<long>(stateSet.size());
    const long algebraicUnknowns = static_cast<long>(prob.algebraicVars.size());
    const long algebraicEqs = static_cast<long>(prob.algebraicEquationIds.size());
    if (algebraicEqs != algebraicUnknowns) {
        std::ostringstream os;
        os << "Algebraic subsystem is not square: " << algebraicEqs
           << " equations vs " << algebraicUnknowns << " unknowns "
           << "(after removing " << kStates << " state variable(s)). "
           << "The model may be over- or under-constrained.";
        prob.diagnostics.push_back(os.str());
    }

    prob.valid = true;
    return prob;
}

const StateVariable* IntegralProblem::state(const std::string& name) const {
    for (const auto& s : states)
        if (s.name == name) return &s;
    return nullptr;
}

bool IntegralProblem::isInternal(const std::string& name) const {
    for (const auto& n : internalNames)
        if (toLower(n) == toLower(name)) return true;
    return false;
}

}  // namespace coolsolve
