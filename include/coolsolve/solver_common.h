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

namespace coolsolve {

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
 * @brief Record one iteration in a SolverTrace.
 *
 * Centralises the boilerplate that every solver repeats to push an
 * iteration record.
 */
struct SolverTrace;  // forward — full definition in solver.h

}  // namespace coolsolve
