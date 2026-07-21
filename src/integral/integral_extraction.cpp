/**
 * @file integral_extraction.cpp
 * @brief Extract an `IntegralProblem` from a parsed IR.
 *
 * Walks every equation looking for top-level `INTEGRAL(integrand, var, lo, hi[, step])`
 * calls, classifies state vs algebraic variables, validates structural
 * preconditions, and runs a conservative high-index rejection.  See
 * `docs/integral_table.md` §1 (Mathematical model) and §3 (Architecture).
 *
 * The initial state is *not* extracted here: at t = t0 the integral term has a
 * zero-width interval and evaluates to 0, so `y(t0)` falls out of the algebraic
 * solve at the first step (handled by `IntegralSolver`).
 */
#include "coolsolve/integral/integral_problem.h"

#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"

#include <algorithm>
#include <cctype>
#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <functional>
#include <set>
#include <sstream>

namespace coolsolve {

namespace {

// ---- small helpers ---------------------------------------------------------

std::string toLower(const std::string& s) {
    std::string o(s.size(), '\0');
    std::transform(s.begin(), s.end(), o.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return o;
}

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
        if (findIntegralCall(eq.rhs)) return true;
    return false;
}

IntegralProblem extractIntegralProblem(const IR& ir,
                                       const StructuralAnalysisResult& /*analysis*/) {
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

    // -- 1. Walk equations, harvest INTEGRAL declarations. -------------------
    for (const auto& eq : ir.getEquations()) {
        const FunctionCall* ic = findIntegralCall(eq.rhs);
        if (!ic) continue;
        prob.hasIntegral = true;

        // LHS must be a plain variable (the state).
        if (!eq.lhs || !eq.lhs->is<Variable>() || !eq.lhs->as<Variable>().indices.empty()) {
            return fail("INTEGRAL() must appear on the right-hand side of a scalar variable "
                        "equation (line " + std::to_string(eq.sourceLine) + ").");
        }
        std::string stateName = eq.lhs->as<Variable>().name;

        // Validate argument count: 4 (auto step) or 5 (fixed step).
        const auto& args = ic->args;
        if (args.size() != 4 && args.size() != 5) {
            return fail("INTEGRAL() expects 4 or 5 arguments (integrand, var, lo, hi[, step]); "
                        "got " + std::to_string(args.size()) + " (line " +
                        std::to_string(eq.sourceLine) + ").");
        }

        // 2nd argument must be a plain variable (the integration variable).
        if (!args[1] || !args[1]->is<Variable>() || !args[1]->as<Variable>().indices.empty()) {
            return fail("INTEGRAL() integration variable must be a scalar variable "
                        "(line " + std::to_string(eq.sourceLine) + ").");
        }
        std::string integVar = args[1]->as<Variable>().name;

        // Integrand: usually a variable (the derivative); record both ways.
        StateVariable sv;
        sv.name = stateName;
        sv.integrandExpr = args[0];
        if (args[0] && args[0]->is<Variable>() && args[0]->as<Variable>().indices.empty())
            sv.integrandVar = args[0]->as<Variable>().name;

        // Base expression = RHS with the INTEGRAL call blanked to 0; this is
        // the part that evaluates to y(t0) once the algebraic system is solved
        // at t = t0 (handled by IntegralSolver).
        sv.baseExpr = cloneExpr(eq.rhs, [](const ExprPtr& e) -> ExprPtr {
            if (e && e->is<FunctionCall>() && toLower(e->as<FunctionCall>().name) == "integral")
                return makeNumber(0.0, e->sourceLineNumber);
            return nullptr;
        });
        sv.equationId = eq.id;

        // Limits / step: prefer constant folding, keep expressions for Phase 5.
        double lo = 0.0, hi = 0.0, step = 0.0;
        bool loOk = evalConstantExpr(args[2], lo);
        bool hiOk = evalConstantExpr(args[3], hi);
        bool stepOk = (args.size() == 5) && evalConstantExpr(args[4], step);

        // -- 2. Validate consistency across all INTEGRAL equations. -----------
        if (prob.states.empty()) {
            prob.integrationVar = integVar;
            prob.lowerLimitExpr = args[2];
            prob.upperLimitExpr = args[3];
            prob.lowerLimit = lo;
            prob.upperLimit = hi;
            prob.limitsConstant = loOk && hiOk;
            if (args.size() == 5 && stepOk) prob.fixedStep = step;
        } else {
            if (integVar != prob.integrationVar) {
                return fail("All INTEGRAL() calls must share the same integration variable; "
                            "found '" + integVar + "' and '" + prob.integrationVar + "'.");
            }
            double eLo = 0.0, eHi = 0.0;
            bool eLoOk = evalConstantExpr(args[2], eLo);
            bool eHiOk = evalConstantExpr(args[3], eHi);
            if (prob.limitsConstant && eLoOk && eHiOk &&
                (eLo != prob.lowerLimit || eHi != prob.upperLimit)) {
                return fail("All INTEGRAL() calls must share the same [lo, hi] interval.");
            }
            if (args.size() == 5 && stepOk && prob.fixedStep == 0.0) {
                prob.fixedStep = step;  // adopt the first explicit step
            }
        }

        prob.states.push_back(sv);
        prob.integralEquationIds.push_back(eq.id);
    }

    if (!prob.hasIntegral) {
        prob.valid = false;
        prob.errorMessage = "no INTEGRAL() equations found";
        return prob;
    }

    // Convenience name lists.
    for (const auto& s : prob.states) {
        prob.stateNames.push_back(s.name);
        if (!s.integrandVar.empty()) prob.derivativeNames.push_back(s.integrandVar);
    }

    // -- 3. Reject nested integrals (multi-variable integration). ------------
    for (const auto& s : prob.states) {
        // The integrand expression must not itself contain an INTEGRAL call.
        if (findIntegralCall(s.integrandExpr)) {
            return fail("Nested INTEGRAL() calls are not supported in this version "
                        "(multi-variable integration). See docs/integral_table.md §7.1.");
        }
    }

    // -- 4. Classify variables: state / algebraic / integration var. ---------
    std::set<std::string, CaseInsensitiveLess> stateSet(prob.stateNames.begin(),
                                                        prob.stateNames.end());
    for (const auto& [name, info] : ir.getVariables()) {
        if (stateSet.count(name)) continue;          // state variable (integrator-owned)
        if (name == prob.integrationVar) continue;   // independent variable (integrator-owned)
        prob.algebraicVars.push_back(name);
    }

    // Algebraic equations = every equation that is not an integral declaration.
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
    const long kStates = static_cast<long>(prob.states.size());
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

}  // namespace coolsolve
