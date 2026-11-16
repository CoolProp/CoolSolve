#pragma once
/**
 * @file solver_common.h
 * @brief Shared utilities for all non-linear solvers.
 *
 * Contains helper functions used by Newton, TrustRegion, LM, and other
 * solvers to avoid code duplication across implementations.
 */

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string_view>

namespace coolsolve {

// ============================================================================
// Weighted residual (convergence tests)
// ============================================================================

/// Smallest relative weight of a residual: an equation of terms below it is met at |F| < tolerance * this
/// (1E-21 at the default tolerance): terms of 1E-23 are rounding noise of the iteration, not values to resolve.
constexpr double kMinResidualWeight = 1e-12;

/// Rounding level of an equation, in ulp of the scale S of its terms: the residual of a sum of tens of terms,
/// each from a CoolProp call or an AD product with a few ulp of error, cannot be resolved below it.
constexpr double kRoundingUlps = 64.0;

/// Relative weight w of an equation of finite scale S in a block of largest finite scale \p blockScale (`weightedResidual`).
inline double relativeResidualWeight(double scale, double blockScale) {
    return blockScale < 1.0 ? std::max(scale, kMinResidualWeight) : 1.0;
}

/// Rounding level 64 eps S of an equation of finite scale S (`kRoundingUlps`).
inline double residualRoundingLevel(double scale) {
    return kRoundingUlps * std::numeric_limits<double>::epsilon() * scale;
}

/**
 * @brief Weighted residual rho of an equation of residual \p f and scale S = \p scale (`EvaluationResult::scale`)
 *        in a block whose largest finite scale is \p blockScale (`blockResidualScale`; S itself for one equation).
 *
 * A convergence test `rho < thr` (thr = `tolerance`, `100 * tolerance`, `lsRelaxedTolerance`...) is exactly
 * |F| < max(thr * w, 64 eps S), with
 *   - w = max(S, 1E-12) in a block in which every equation has terms below 1 (blockScale < 1): the test is relative
 *     to the terms (an equation of terms 1E-20 is met at 1E-21, not at 1E-9, which any value below 1E-9 satisfies);
 *   - w = 1 otherwise, the absolute test for every equation of the block (an equation of tiny terms beside
 *     enthalpies of 1E6 is rounding noise of the block: examples `condenser_3zones`, library `CSL-0127`);
 *   - 64 eps S, the rounding level (`kRoundingUlps`), the same for every threshold: a residual within it gives
 *     rho = 0. It exceeds tolerance * w only for S > 7E4 (default tolerance).
 * A non-finite S (an infinite derivative, `sqrt` at 0) gives the absolute test rho = |F|; a NaN \p f gives NaN.
 */
inline double weightedResidual(double f, double scale, double blockScale) {
    const double r = std::abs(f);
    if (!std::isfinite(scale)) return r;
    if (r <= residualRoundingLevel(scale)) return 0.0;   // false for a NaN f
    return r / relativeResidualWeight(scale, blockScale);
}

/// The largest finite entry of \p scale (0 when there is none): the `blockScale` of `weightedResidual`.
template <class Scale>   // Eigen::VectorXd or std::vector<double>
double blockResidualScale(const Scale& scale) {
    double blockScale = 0.0;
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(scale.size()); ++i)
        if (std::isfinite(scale[i])) blockScale = std::max(blockScale, static_cast<double>(scale[i]));
    return blockScale;
}

/**
 * @brief ||rho||_inf with rho_i the weighted residual of equation i (`weightedResidual`); NaN if a residual is NaN.
 *
 * An empty or mismatching \p scale (a problem that gives none: no Jacobian has been evaluated yet, or a
 * hand-built problem) means the plain ||F||_inf.
 */
template <class Scale>   // Eigen::VectorXd or std::vector<double>
double weightedResidualNorm(const Eigen::VectorXd& F, const Scale& scale) {
    const bool weighted = static_cast<Eigen::Index>(scale.size()) == F.size();
    const double blockScale = weighted ? blockResidualScale(scale) : 0.0;
    double norm = 0.0;
    for (Eigen::Index i = 0; i < F.size(); ++i) {
        const double rho = weighted ? weightedResidual(F[i], scale[i], blockScale) : std::abs(F[i]);
        if (!(rho <= norm)) {
            if (std::isnan(rho)) return rho;
            norm = rho;
        }
    }
    return norm;
}

/**
 * @brief Distance to convergence ||F_i / a_i||_inf, with a_i = max(tolerance * w_i, 64 eps S_i) the residual at which
 *        equation i is met (`weightedResidual`; a_i = tolerance for a non-finite S): the `relativeTolerance` test
 *        compares it with its value at the start. A start of huge terms (a residual of 1E20 at the guess) then does
 *        not make a point that is still far from convergence look converged, as |F| / |F0| would.
 *        An empty or mismatching \p scale gives ||F||_inf / tolerance.
 */
template <class Scale>   // Eigen::VectorXd or std::vector<double>
double convergenceDistance(const Eigen::VectorXd& F, const Scale& scale, double tolerance) {
    const bool weighted = static_cast<Eigen::Index>(scale.size()) == F.size();
    const double blockScale = weighted ? blockResidualScale(scale) : 0.0;
    double norm = 0.0;
    for (Eigen::Index i = 0; i < F.size(); ++i) {
        double allowed = tolerance;
        if (weighted && std::isfinite(scale[i]))
            allowed = std::max(tolerance * relativeResidualWeight(scale[i], blockScale), residualRoundingLevel(scale[i]));
        const double d = std::abs(F[i]) / allowed;
        if (!(d <= norm)) {
            if (std::isnan(d)) return d;
            norm = d;
        }
    }
    return norm;
}

/**
 * @brief Compute automatic variable-scaling factors.
 *
 * Variables with widely different magnitudes (T~300, P~1e7, h~3e5) produce
 * ill-conditioned Jacobians.  Scaling each variable by its order of magnitude
 * improves convergence.
 *
 * @param x  Current iterate (unscaled)
 * @return   Vector of scaling factors (s_i ≈ 10^floor(log10|x_i|), clamped to [1e-6, 1e6]).
 */
inline Eigen::VectorXd computeScalingFactors(const Eigen::VectorXd& x) {
    const int n = x.size();
    Eigen::VectorXd scale(n);
    for (int i = 0; i < n; ++i) {
        double xi = std::abs(x(i));
        if (xi < 1.0) {
            scale(i) = 1.0;
        } else {
            double log10x = std::floor(std::log10(xi));
            scale(i) = std::pow(10.0, log10x);
            scale(i) = std::max(1e-6, std::min(scale(i), 1e6));
        }
    }
    return scale;
}

// ============================================================================
// Non-Monotone Merit History
// ============================================================================

/**
 * @brief Tracks recent merit function values for non-monotone line search.
 *
 * Standard (monotone) line search requires φ(x_{k+1}) < φ(x_k) at every step.
 * This can trap solvers in narrow curved valleys or near saddle points.
 *
 * Non-monotone line search (Grippo, Lampariello, Lucidi 1986) relaxes this
 * to: φ(x_{k+1}) < max(φ(x_{k-M+1}), ..., φ(x_k)) + sufficient decrease.
 * The solver can temporarily accept larger merit values to escape local traps,
 * while still guaranteeing global convergence.
 *
 * **Bounded reference** (`boundedRef`): the raw max can be arbitrarily larger
 * than the current merit when a single early bad iterate inflates the window.
 * This makes the Armijo condition trivially satisfied, allowing steps that
 * destroy recent progress.  `boundedRef(phi, R)` caps the reference at
 * `R × phi`, preventing the reference from exceeding `R` times the current
 * merit.  Default ratio R = 10.
 *
 * Usage in solver iteration loop:
 *   NonMonotoneHistory history(options.lsNonMonotoneMemory);
 *   // each iteration:
 *   double phi = 0.5 * F.squaredNorm();
 *   history.push(phi);
 *   double refPhi = history.boundedRef(phi);  // bounded non-monotone reference
 *
 * When memory == 1, this degenerates to monotone line search (maxValue = current phi).
 */
struct NonMonotoneHistory {
    explicit NonMonotoneHistory(int memory) : memory_(std::max(memory, 1)) {}

    /// Record a new merit function value.
    void push(double phi) {
        values_.push_back(phi);
        if (static_cast<int>(values_.size()) > memory_)
            values_.erase(values_.begin());
    }

    /// Maximum of the stored merit values (the non-monotone reference).
    double maxValue() const {
        if (values_.empty()) return 0.0;
        return *std::max_element(values_.begin(), values_.end());
    }

    /**
     * @brief Bounded non-monotone reference value.
     *
     * Returns min(max(history), currentPhi * maxRatio).  This prevents
     * the reference from being arbitrarily larger than the current merit,
     * which can cause the Armijo condition to admit steps that destroy
     * recent progress (observed when a single early bad iterate inflates
     * the history window far above the current merit).
     *
     * @param currentPhi  Current merit value φ(x_k)
     * @param maxRatio    Maximum allowed ratio refPhi / currentPhi (e.g. 10)
     * @return Bounded reference value for the non-monotone Armijo condition
     */
    double boundedRef(double currentPhi, double maxRatio = 10.0) const {
        double mv = maxValue();
        double cap = currentPhi * maxRatio;
        return std::min(mv, cap);
    }

    bool empty() const { return values_.empty(); }
    int size() const { return static_cast<int>(values_.size()); }

private:
    int memory_;
    std::vector<double> values_;
};

/**
 * @brief Acceptance test of a backtracking line-search step.
 *
 * Non-monotone Armijo condition (Grippo et al. 1986)
 *   φ(x + λd) ≤ refPhi + α λ ∇φ·d,  with ∇φ·d = -2 φ₀ for a Newton step,
 * plus two relaxations near a minimum, where φ only changes by rounding
 * noise and the strict test would reject every step.  Used by the Newton
 * line search and by the tearing outer iteration.
 *
 * @param phiNew  Merit at the trial point (non-finite → rejected)
 * @param phi0    Merit at the start of the step
 * @param refPhi  Non-monotone reference (`NonMonotoneHistory::boundedRef`)
 * @param lambda  Step length of the trial point
 * @param alpha   Armijo parameter (`SolverOptions::lsAlpha`)
 */
inline bool lineSearchAccepts(double phiNew, double phi0, double refPhi,
                              double lambda, double alpha) {
    if (!std::isfinite(phiNew)) return false;
    if (phiNew <= refPhi - alpha * lambda * 2.0 * phi0) return true;
    // Near-minimum relaxation
    if (phi0 < 0.05 && phiNew < phi0 * 2.0) return true;
    if (lambda < 0.01 && phi0 < 0.1 && phiNew <= phi0 * 1.5) return true;
    return false;
}

/**
 * @brief Record one iteration in a SolverTrace.
 *
 * Centralises the boilerplate that every solver repeats to push an
 * iteration record.
 */
struct SolverTrace;  // forward — full definition in solver.h

// ============================================================================
// Fatal error detection
// ============================================================================

/**
 * @brief Check if an exception message indicates a fatal evaluation error
 *        that no solver strategy can recover from.
 *
 * Fatal errors include unsupported functions or unknown fluids, where
 * iterating with different values won't help, and a `CALL ERROR` of the model
 * (EES halts the calculation as soon as any iterate reaches it; the message is
 * "CALL ERROR: <text>", see ExpressionEvaluator::evaluateProcedureCall).
 * A CoolProp failure is not fatal: it is recoverable (a penalty value).
 */
inline bool isFatalEvaluationError(std::string_view what) {
    return what.find("Unknown or unsupported function") != std::string_view::npos
        || what.find("Unknown fluid") != std::string_view::npos
        || what.find("CALL ERROR: ") != std::string_view::npos;
}

/**
 * @brief Rethrow the exception being handled if it is a fatal evaluation error.
 *
 * Call it first in every `catch` block that would otherwise treat an evaluation
 * failure as a recoverable point (backtrack, shrink the step, try another
 * start): a fatal error then leaves the solver as an exception, which
 * `Solver::runSolverStrategy` turns into `SolverStatus::EvaluationError` with
 * the message. Does nothing for any other exception.
 */
inline void rethrowIfFatalEvaluationError() {
    try {
        throw;
    } catch (const std::exception& e) {
        if (isFatalEvaluationError(e.what())) throw;
    } catch (...) {
    }
}

/**
 * @brief Give \p problem (a `NonLinearSolver::Problem`) its residual scale when it has none yet, by one
 *        evaluation with the Jacobian at \p x, for a solver that otherwise evaluates residuals only (KINSOL
 *        Picard and fixed point, BisectionND): its tests are then weighted as those of the other solvers.
 *        A failed evaluation leaves the scale empty (the absolute tests); a fatal evaluation error is rethrown.
 */
template <class Problem>
void ensureResidualScale(Problem& problem, const Eigen::VectorXd& x) {
    if (problem.residualScale.size() != 0) return;
    Eigen::VectorXd F(x.size());              // sized, as the solvers' own evaluations are
    Eigen::MatrixXd J(x.size(), x.size());
    try {
        problem.evaluate(x, F, J, true);
    } catch (...) {
        rethrowIfFatalEvaluationError();
    }
}

}  // namespace coolsolve
