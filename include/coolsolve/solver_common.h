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

namespace coolsolve {

// ============================================================================
// Weighted residual (convergence tests)
// ============================================================================

/// Smallest weight of a residual (`residualWeight`): an equation of terms below it is met at |F| < tolerance * this
/// (1E-21 at the default tolerance): terms of 1E-23 are rounding noise of the iteration, not values to resolve.
constexpr double kMinResidualWeight = 1e-12;

/**
 * @brief Weight w of the residual of an equation of scale S (`EvaluationResult::scale`) in a block whose
 *        largest scale is \p blockScale (the scale of the equation itself for a block of one equation).
 *
 * The convergence tests compare rho = |F| / w with `tolerance`, `100 * tolerance`, `lsRelaxedTolerance`...
 *   - a block in which every equation has terms below 1 (blockScale < 1): w = max(S, 1E-12), the test is
 *     relative to the terms. An equation of terms 1E-20 is met at a residual of 1E-21, not at 1E-9, which any
 *     value below 1E-9 satisfies;
 *   - a block with an equation of terms of order 1 or more: w = 1, today's absolute test for all its equations.
 *     The residual of an equation of tiny terms beside enthalpies of 1E6 is rounding noise of the block (a
 *     pressure drop of 1E-23 Pa, `MIN(a, b)` of two noisy 1E-8: examples `condenser_3zones`, library
 *     `CSL-0127`), and a relative test on it would never be met;
 *   - S > 64 eps / tolerance (7E4 at the default tolerance): w = 64 eps S / tolerance, so that the test is
 *     not below 64 ulp of the terms, a rounding level that the residual cannot reach.
 */
inline double residualWeight(double scale, double blockScale, double tolerance) {
    const double relativeWeight = blockScale < 1.0 ? std::min(1.0, std::max(scale, kMinResidualWeight)) : 1.0;
    const double roundingWeight = 64.0 * std::numeric_limits<double>::epsilon() * scale / tolerance;
    return std::max(relativeWeight, roundingWeight);
}

/**
 * @brief ||rho||_inf with rho_i = |F_i| / w_i (`residualWeight`); NaN if a residual is NaN.
 *
 * An empty or mismatching \p scale (a problem that gives none: no Jacobian has been evaluated yet, or a
 * hand-built problem) means w = 1: the plain ||F||_inf.
 */
template <class Scale>   // Eigen::VectorXd or std::vector<double>
double weightedResidualNorm(const Eigen::VectorXd& F, const Scale& scale, double tolerance) {
    const bool weighted = static_cast<Eigen::Index>(scale.size()) == F.size();
    double blockScale = 0.0;
    if (weighted)
        for (Eigen::Index i = 0; i < F.size(); ++i) blockScale = std::max(blockScale, static_cast<double>(scale[i]));
    double norm = 0.0;
    for (Eigen::Index i = 0; i < F.size(); ++i) {
        const double rho = weighted ? std::abs(F[i]) / residualWeight(scale[i], blockScale, tolerance) : std::abs(F[i]);
        if (!(rho <= norm)) {
            if (std::isnan(rho)) return rho;
            norm = rho;
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
// Equilibrated linear solve (Newton steps)
// ============================================================================

/// The power of two nearest below 1/m for a magnitude m > 0 (1 for 0 or a non-finite m): scaling by it is exact.
inline double inversePowerOfTwo(double m) {
    if (!(m > 0.0) || !std::isfinite(m)) return 1.0;
    int exponent = 0;
    std::frexp(m, &exponent);   // m = f 2^exponent, f in [0.5, 1)
    return std::ldexp(1.0, -exponent);
}

/**
 * @brief The solution of J dx = rhs by a column-pivoting QR factorisation of the equilibrated matrix R J C.
 *
 * R and C are diagonal: each row of J is scaled so that its largest entry is in [0.5, 1), then each column of the
 * result likewise, by powers of two (exact: the scaling adds no rounding; the row/column equilibration of LAPACK's
 * xGEEQU). A block whose equations have terms of very different sizes (an energy balance of 1E9 J/kmol beside mole
 * fractions of 1E-2, register CS-BUG-BLOCK-SCALE-SPREAD) has a Jacobian whose condition number is ~1E19 as it stands
 * and ~1E5 equilibrated: the QR, which treats the pivots below eps times its largest one as zero, then keeps every
 * direction of the small equations.
 * @param rank  when given, the rank the QR found for the equilibrated matrix
 */
inline Eigen::VectorXd solveEquilibrated(const Eigen::MatrixXd& J, const Eigen::VectorXd& rhs, Eigen::Index* rank = nullptr) {
    const Eigen::Index rows = J.rows(), cols = J.cols();
    Eigen::VectorXd r(rows), c(cols);
    for (Eigen::Index i = 0; i < rows; ++i) r(i) = inversePowerOfTwo(J.row(i).cwiseAbs().maxCoeff());
    Eigen::MatrixXd A = r.asDiagonal() * J;
    for (Eigen::Index j = 0; j < cols; ++j) c(j) = inversePowerOfTwo(A.col(j).cwiseAbs().maxCoeff());
    A = A * c.asDiagonal();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    if (rank) *rank = qr.rank();
    return c.cwiseProduct(qr.solve(r.cwiseProduct(rhs)));
}

/**
 * @brief The Newton step: the solution of J dx = rhs by the column-pivoting QR of J, as the solvers always took it;
 *        when that step leaves linear equations unmet that the convergence test counts, the solution of the
 *        equilibrated system if it meets them better.
 *
 * The plain step is kept when its linear residual J dx - rhs passes the convergence test of the block
 * (weightedResidualNorm with the residual scale and the tolerance of the solver, below \p tolerance): every block it
 * solves takes the same steps and iterations as before, a singular but consistent block (a family of solutions, of
 * which the equilibrated QR would pick another member) and the equations of rounding noise of a block (their weight,
 * residualWeight) included. Otherwise (CS-BUG-BLOCK-SCALE-SPREAD: the QR dropped a direction of the small equations, and
 * the step leaves their residuals untouched) solveEquilibrated is solved too, and the step whose linear residual is
 * smaller in that test is taken. The cost is one product J dx, and a second factorisation only in that case.
 * @param residualScale  the scale S of the residuals (`Problem::residualScale`; empty: the absolute test)
 * @param rank  when given, the rank the QR found for the matrix of the step taken (KINSOL's singularity test)
 */
template <class Scale>   // Eigen::VectorXd or std::vector<double>
Eigen::VectorXd solveNewtonStep(const Eigen::MatrixXd& J, const Eigen::VectorXd& rhs, const Scale& residualScale,
                                double tolerance, Eigen::Index* rank = nullptr) {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(J);
    Eigen::VectorXd dx = qr.solve(rhs);
    if (rank) *rank = qr.rank();
    auto unmet = [&](const Eigen::VectorXd& step) {   // the weighted linear residual (NaN or inf: not finite)
        return step.allFinite() ? weightedResidualNorm(Eigen::VectorXd(J * step - rhs), residualScale, tolerance)
                                : std::numeric_limits<double>::infinity();
    };
    const double plain = unmet(dx);
    if (plain < tolerance) return dx;
    Eigen::Index equilibratedRank = 0;
    Eigen::VectorXd equilibrated = solveEquilibrated(J, rhs, &equilibratedRank);
    if (!(unmet(equilibrated) < plain)) return dx;
    if (rank) *rank = equilibratedRank;
    return equilibrated;
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
inline bool isFatalEvaluationError(const std::string& what) {
    return what.find("Unknown or unsupported function") != std::string::npos
        || what.find("Unknown fluid") != std::string::npos
        || what.find("CALL ERROR: ") != std::string::npos;
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

}  // namespace coolsolve
