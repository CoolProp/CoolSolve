#pragma once
/**
 * @file integrator.h
 * @brief Abstract integrator interface and options for the equation-based
 *        dynamic solver (`INTEGRAL` / `$IntegralTable`).
 *
 * This module is intentionally decoupled from the EES parser/IR/evaluator so
 * that the numerical methods (Euler, RK4, Dormand-Prince RK45) can be unit
 * tested in isolation against analytical solutions (see
 * `tests/test_integrator_*.cpp`).  The orchestrator that wires these methods
 * to the algebraic `Solver` lives in `integral_solver.h`.
 *
 * References:
 *   - Hairer, Nørsett, Wanner, *Solving Ordinary Differential Equations I*
 *     (Nonstiff Problems), Springer.
 *   - Dormand & Prince (1980), "A family of embedded Runge-Kutta formulae",
 *     J. Comp. Appl. Math. 6(1):19-26.
 *   - Press et al., *Numerical Recipes* 3rd ed., Ch. 17.
 */

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <string>

namespace coolsolve {

/**
 * @brief Configuration for a single integrator instance.
 *
 * All fields are inert by default; a non-integral model never constructs an
 * integrator (zero overhead — see `docs/contributing.md §1.7`).
 */
struct IntegratorOptions {
    /** Available integration methods. */
    enum Method {
        EulerExplicit,  ///< Fixed-step explicit Euler (1st order).
        EulerImplicit,  ///< Fixed-step implicit Euler (1st order, A-stable).
        RK4,            ///< Fixed-step classic Runge-Kutta (4th order). Default.
        RK45            ///< Variable-step Dormand-Prince (embedded 4(5)). Adaptive.
    };

    Method method = RK4;          ///< Integration method (default: RK4).
    double fixedStep = 0.0;       ///< Requested step size. 0 => derive from `maxSteps`.
    int    maxSteps  = 1000;      ///< Upper bound on the number of steps.
    double relTol    = 1e-6;      ///< RK45 relative error control.
    double absTol    = 1e-9;      ///< RK45 absolute error floor.
    double minStep   = 0.0;       ///< Minimum step size for RK45 (0 = auto).
    double maxStep   = 0.0;       ///< Maximum step size for RK45 (0 = auto).
    bool   richardson = false;    ///< Richardson extrapolation (fixed-step methods only).
};

/**
 * @brief Right-hand side function `f(t, y)` returning `dy/dt`.
 *
 * A pure function with no IR/evaluator coupling, so the integrators can be
 * tested analytically.  In Phase 5 the `IntegralSolver` adapts the algebraic
 * solver into this signature.
 */
using RHSFunction = std::function<Eigen::VectorXd(double t, const Eigen::VectorXd& y)>;

/**
 * @brief Result of a single integration step.
 */
struct StepResult {
    Eigen::VectorXd yNew;          ///< State at `t + stepTaken`.
    double stepTaken     = 0.0;    ///< Step size actually used.
    double nextStep      = 0.0;    ///< Recommended next step (RK45 adaptive; = h for fixed).
    double errorEstimate = 0.0;    ///< Normalised local error (RK45); 0 for fixed-step.
    bool   accepted      = true;   ///< Whether the step satisfied the tolerance (RK45).
    int    rhsEvaluations = 0;     ///< Number of RHS evaluations consumed (diagnostics).
};

/**
 * @brief Abstract base class for one-step ODE integrators.
 *
 * Each concrete integrator advances the state from `(t, y)` by a proposed
 * step `h`.  For fixed-step methods `h` is used as-is; for the adaptive
 * `RK45` method `h` is the proposed step and `StepResult::stepTaken` / the
 * error estimate drive the outer step-size controller.
 */
class Integrator {
public:
    virtual ~Integrator() = default;

    /**
     * @brief Advance one step of size `h` from `(t, y)`.
     * @param t   Current independent variable.
     * @param y   Current state.
     * @param rhs Right-hand side `f(t, y)`.
     * @param h   Proposed step size.
     * @param opt Integrator options (tolerances etc.).
     */
    virtual StepResult step(double t, const Eigen::VectorXd& y,
                            const RHSFunction& rhs, double h,
                            const IntegratorOptions& opt) = 0;

    /// Human-readable method name (diagnostics / debug reports).
    virtual const char* name() const = 0;

    /**
     * @brief Global order of accuracy `p` (error ∝ h^p).
     *
     * Used by the Richardson extrapolation wrapper to pick the correct
     * combination factor 2^p.  Defaults to 1; concrete methods override.
     */
    virtual int order() const { return 1; }
};

/**
 * @brief Factory: construct the concrete integrator for `method`.
 *
 * Mirrors the `createSolver(SolverStrategy)` pattern used by the algebraic
 * solver pipeline (`src/solver.cpp`).
 */
std::unique_ptr<Integrator> createIntegrator(IntegratorOptions::Method method);

/**
 * @brief Wrap `base` so each step is refined by Richardson extrapolation.
 *
 * Applies `I ≈ (4·I(h/2) − I(h)) / 3` by taking one `h` step and two `h/2`
 * steps from the base integrator and combining the results.  Only valid for
 * fixed-step methods; reduces trunc error by 1–2 orders (EES docs,
 * `richardson_extrapolation.htm`).  Returns `base` unchanged if Richardson
 * is disabled.
 */
std::unique_ptr<Integrator> wrapRichardson(std::unique_ptr<Integrator> base,
                                           bool enable);

/// Convert a method enum to its canonical string name.
std::string methodToString(IntegratorOptions::Method m);

/**
 * @brief Parse an integration method name (case-insensitive).
 * @return true if parsing succeeded.
 */
bool parseIntegralMethod(const std::string& s, IntegratorOptions::Method& out);

}  // namespace coolsolve
