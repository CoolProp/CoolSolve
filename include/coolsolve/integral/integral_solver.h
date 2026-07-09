#pragma once
/**
 * @file integral_solver.h
 * @brief `IntegralSolver` — orchestrates the time-marching solve for an
 *        equation-based dynamic model.
 *
 * Reuses the algebraic `Solver` *unmodified* at every step: the integral
 * equations are stripped into a reduced "algebraic subsystem" IR, the state
 * variables `y` and the integration variable `t` are fixed as external values
 * each step, and the algebraic solve returns the derivatives `f = dy/dt`. An
 * `Integrator` (Phase 1) then advances the state. See
 * `docs/integral_table_plan.md` §Phase 5 and §3.4.
 */
#include "coolsolve/ast.h"
#include "coolsolve/integral/integral_problem.h"
#include "coolsolve/integral/integral_table.h"
#include "coolsolve/integral/integrator.h"
#include "coolsolve/ir.h"
#include "coolsolve/solver.h"
#include "coolsolve/structural_analysis.h"

#include <memory>
#include <string>
#include <vector>

namespace coolsolve {

/**
 * @brief Result of a dynamic solve, mirroring `SolveResult` for the algebraic
 *        solver plus the trajectory table.
 */
struct IntegralSolveResult {
    bool success = false;
    std::string errorMessage;

    IntegralProblem problem;                ///< Extracted dynamic structure.

    IntegralTable table;                    ///< Tabulated trajectory.
    int totalSteps = 0;                     ///< Accepted steps taken.
    int rejectedSteps = 0;                  ///< RK45 steps rejected (diagnostics).
    std::vector<double> acceptedStepSizes;  ///< Per-step `h` (diagnostics).

    /// Final-step algebraic result (variables, status, block results).
    SolveResult algebraicResult;
};

/**
 * @brief Time-marching solver for equation-based dynamic models.
 */
class IntegralSolver {
public:
    /**
     * @param program  Parsed AST (used to build the algebraic subsystem IR;
     *                  function/procedure definitions are preserved).
     * @param ir       Full IR (already inferred/initialised by the runner).
     * @param analysis Structural analysis of the full IR.
     * @param options  Algebraic solver options (reused per step).
     */
    IntegralSolver(const Program& program, const IR& ir,
                   const StructuralAnalysisResult& analysis,
                   const SolverOptions& options);

    /**
     * @brief March the system from `t0` to `tf` and tabulate the trajectory.
     *
     * The method, step size, and tolerances come from `intOpt`; the algebraic
     * solver options come from the constructor's `options`.
     */
    IntegralSolveResult solve(const IntegratorOptions& intOpt);

    /// Read-only access to the extracted problem (for diagnostics/debug output).
    const IntegralProblem& problem() const { return problem_; }

private:
    const IR& fullIr_;
    SolverOptions solverOpts_;              // copied (caller may pass a temporary)
    IntegralProblem problem_;

    std::unique_ptr<IR> reducedIr_;                          ///< algebraic subsystem
    std::unique_ptr<StructuralAnalysisResult> reducedAnalysis_;
    std::unique_ptr<Solver> algebraicSolver_;                ///< reused each step

    /// Evaluate the RHS f(t, y) by fixing (t, y) and solving the algebraic
    /// subsystem. Throws std::runtime_error on algebraic failure.
    Eigen::VectorXd evaluateRHS(double t, const Eigen::VectorXd& y);

    /// Read the i-th state's derivative from the last algebraic solve.
    double readDerivative(int stateIndex) const;

    /// Evaluate a state's base expression (initial value at t0) against the
    /// given variable map.
    double evalBaseExpr(const ExprPtr& expr,
                        const std::map<std::string, double, CaseInsensitiveLess>& vars) const;

    // March-loop scratch state (kept as members so recordRow can read them).
    Eigen::VectorXd currentY_;
    SolveResult lastAlgebraic_;

    /// Driver equations (`v = <literal>`) for integrator-owned variables
    /// (states + integration var). Each step the literal is updated to the
    /// integrator's current value so the algebraic solve sees them as fixed.
    std::map<std::string, ExprPtr> driverLiterals_;
    void setDriverValue(const std::string& name, double value);
};

}  // namespace coolsolve
