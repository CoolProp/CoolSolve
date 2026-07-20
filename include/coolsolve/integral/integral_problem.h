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
#include <string>
#include <vector>

namespace coolsolve {

struct IR;                 // forward (full def in ir.h)
struct StructuralAnalysisResult;  // forward

/**
 * @brief One state variable of the dynamic system.
 *
 * Declared by an equation of the form
 *     y = y0 + INTEGRAL(dydt, t, t0, tf)
 * Here `name`="y", `integrandVar`="dydt" (the derivative, an algebraic
 * variable computed each step), `baseExpr` is the non-integral part (y0)
 * used to recover the initial state at t0.
 */
struct StateVariable {
    std::string name;            ///< The state variable (LHS), e.g. "y".
    std::string integrandVar;    ///< Derivative variable name (1st INTEGRAL arg), e.g. "dydt".
    ExprPtr integrandExpr;       ///< The 1st INTEGRAL argument expression (usually a Variable).
    ExprPtr baseExpr;            ///< Non-integral RHS part (evaluates to y(t0) at t=t0).
    int equationId = -1;         ///< Equation declaring this state.
};

/**
 * @brief Description of a dynamic model ready for `IntegralSolver`.
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

    std::vector<StateVariable> states;
    std::vector<std::string> stateNames;        ///< Convenience: names of `states`.
    std::vector<std::string> derivativeNames;   ///< The dydt variables (algebraic, computed/step).

    /// Variables solved algebraically each step (all non-state unknowns,
    /// including the derivative names).
    std::vector<std::string> algebraicVars;

    std::vector<int> integralEquationIds;   ///< The `y = ... integral(...)` equations.
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
 */
IntegralProblem extractIntegralProblem(const IR& ir,
                                       const StructuralAnalysisResult& analysis);

}  // namespace coolsolve
