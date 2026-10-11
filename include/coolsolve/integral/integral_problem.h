#pragma once
/**
 * @file integral_problem.h
 * @brief `IntegralProblem` — the structural description of an equation-based
 *        dynamic model, extracted from the IR by `extractIntegralProblem()`.
 *
 * Classifies every equation/variable into the semi-explicit index-1 DAE form
 *     y' = f(t, y, z)    (state variables y, expressed via INTEGRAL)
 *     0  = g(t, y, z)    (algebraic variables z, solved each step)
 * and validates the structural preconditions.  See `docs/integral_table.md`
 * §1 (Mathematical model) and §3 (Architecture).
 */
#include "coolsolve/ast.h"
#include "coolsolve/integral/integral_table.h"
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace coolsolve {

class IR;                  // forward (full def in ir.h)
struct StructuralAnalysisResult;  // forward

/**
 * @brief One state variable of the dynamic system (one per INTEGRAL call).
 *
 * Two kinds:
 *  - a **user state**, declared by an equation of the plain additive form
 *        y = y0 + INTEGRAL(dydt, t, t0, tf)
 *    Here `name`="y", `baseExpr` is the non-integral part (y0) used to recover
 *    the initial state at t0, and the equation itself is removed from the
 *    algebraic subsystem (the integrator owns `y`);
 *  - an **auxiliary state** (`auxiliary == true`) for an INTEGRAL call that is
 *    part of a larger expression (`y = m*INTEGRAL(...)`, `-INTEGRAL(...)`, two calls
 *    in one equation...). `name` is an internal variable that replaces the call
 *    in the equation (`I(t0) = 0`, `dI/dt` = integrand); the equation stays in the
 *    algebraic subsystem as an ordinary equation.
 * `integrandVar` is always a variable name: the integrand when it is a variable,
 * otherwise an internal variable defined by an extra equation `integrandVar = expr`.
 */
struct StateVariable {
    std::string name;            ///< The state variable (LHS, or internal name when auxiliary).
    std::string integrandVar;    ///< Derivative variable name (1st INTEGRAL arg), e.g. "dydt".
    ExprPtr integrandExpr;       ///< The 1st INTEGRAL argument expression.
    ExprPtr baseExpr;            ///< Non-integral RHS part (evaluates to y(t0) at t=t0); 0 when auxiliary.
    int equationId = -1;         ///< Equation containing the INTEGRAL call.
    bool auxiliary = false;      ///< Internal state of an INTEGRAL embedded in a larger expression.
    std::vector<std::string> baseVariables;   ///< Variables of `baseExpr` (user states only).
    std::string tableName;       ///< Variable shown for this call in the default table (LHS), "" if none.
};

/**
 * @brief One INTEGRAL(integrand, var, lo, hi[, step]) call, with its limit expressions.
 */
struct IntegralCall {
    int equationId = -1;
    int stateIndex = -1;         ///< Index in `IntegralProblem::states`.
    ExprPtr lowerExpr, upperExpr;
    ExprPtr stepExpr;            ///< null when the call has 4 arguments.
    /// Internal variables holding the limits/step once resolved (only when some limit or step
    /// is not a constant): defined by the solver as `loVar = lowerExpr`, ...
    std::string loVar, hiVar, stepVar;
};

/**
 * @brief An equation containing INTEGRAL calls inside a larger expression: the
 *        original sides are kept for reference, `lhs`/`rhs` hold the copies in
 *        which every call is replaced by its auxiliary state variable.
 */
struct RewrittenEquation {
    int equationId = -1;
    ExprPtr originalLhs, originalRhs;   ///< The shared AST nodes (identify the statement).
    ExprPtr lhs, rhs;                   ///< Call-free copies.
};

/**
 * @brief Description of a dynamic model ready for `IntegralSolver`.
 *
 * The integration variable and the interval `[lo, hi]` are shared by all INTEGRAL calls of the
 * model. Limits and step may be constants (`limitsConstant`/`stepConstant`: folded here, the
 * historical zero-overhead path) or arbitrary expressions of the model's variables, resolved by
 * the solver at the initial point (they must not depend on the integration variable or on a state).
 */
struct IntegralProblem {
    bool hasIntegral = false;    ///< True if any equation contains an INTEGRAL() call.

    std::string integrationVar;       ///< Independent variable, e.g. "t".
    ExprPtr lowerLimitExpr;           ///< 3rd INTEGRAL argument.
    ExprPtr upperLimitExpr;           ///< 4th INTEGRAL argument.
    double lowerLimit = 0.0;          ///< Resolved constant limit (when `limitsConstant`).
    double upperLimit = 0.0;
    bool limitsConstant = false;      ///< True if both limits evaluated to constants.
    double fixedStep = 0.0;           ///< Optional fixed step (5th INTEGRAL argument); 0 = auto.
    /// True if every explicit step is a constant (or absent). When `!limitsConstant ||
    /// !stepConstant` the limits/step are resolved by the solver from the initial solve.
    bool stepConstant = true;

    std::vector<IntegralCall> calls;            ///< One entry per INTEGRAL call (same order as `states`).
    std::vector<RewrittenEquation> rewrittenEquations;   ///< Equations with auxiliary states.
    /// Extra equations `<internal variable> = <expression>` for the integrands that are not variables.
    std::vector<std::pair<std::string, ExprPtr>> integrandEquations;
    std::vector<std::string> warnings;          ///< User-visible warnings (reported by the solver).
    /// Internal variables created by the extraction (auxiliary states `__integral_<k>`, integrand
    /// definitions `__integrand_<k>`); they are never shown in user outputs. Names cannot clash
    /// with the model's variables (a numeric suffix is added when needed).
    std::vector<std::string> internalNames;
    bool isInternal(const std::string& name) const;

    std::vector<StateVariable> states;
    std::vector<std::string> stateNames;        ///< Convenience: names of `states`.
    std::vector<std::string> derivativeNames;   ///< The dydt variables (algebraic, computed/step).

    /// Variables solved algebraically each step (all non-state unknowns,
    /// including the derivative names).
    std::vector<std::string> algebraicVars;

    std::vector<int> integralEquationIds;   ///< The `y = ... integral(...)` equations of the user states.
    std::vector<int> auxEquationIds;        ///< Equations with INTEGRAL calls inside a larger expression.
    std::vector<int> algebraicEquationIds;  ///< Everything else (solved each step).

    IntegralTableSpec tableSpec;            ///< From `$IntegralTable` (may be absent).

    bool valid = false;                     ///< True when the model is structurally solvable.
    std::string errorMessage;               ///< First fatal problem (when `!valid`).
    std::vector<std::string> diagnostics;   ///< Non-fatal notes / high-index warnings.

    /// Convenience: state variable by name, or nullptr.
    const StateVariable* state(const std::string& name) const;
};

/**
 * @brief Cheap pre-check: does `ir` contain any INTEGRAL() call?
 *
 * Used by the runner (Phase 6) to dispatch to the integral path with zero
 * overhead on non-integral models.
 */
bool hasIntegral(const IR& ir);

/**
 * @brief Extract and validate the `IntegralProblem` from a parsed IR.
 *
 * Walks every equation, finds top-level INTEGRAL() calls, classifies state vs
 * algebraic variables, validates that all integrals share one integration
 * variable and interval, and runs a conservative high-index rejection
 * (placeholder for future Pantelides-style index reduction).
 *
 * `forceAuxiliary` lists equation ids that must be handled with an auxiliary state even
 * though they have the plain additive form (used by the solver when the base of such an
 * equation turns out to vary along the trajectory).
 */
IntegralProblem extractIntegralProblem(const IR& ir,
                                       const StructuralAnalysisResult& analysis,
                                       const std::set<int>* forceAuxiliary = nullptr);

}  // namespace coolsolve
