#pragma once
/**
 * @file solver_common.h
 * @brief Shared utilities for all non-linear solvers.
 *
 * Contains helper functions used by Newton, TrustRegion, LM, and other
 * solvers to avoid code duplication across implementations.
 */

#include <Eigen/Dense>

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

/**
 * @brief Record one iteration in a SolverTrace.
 *
 * Centralises the boilerplate that every solver repeats to push an
 * iteration record.
 */
struct SolverTrace;  // forward — full definition in solver.h

}  // namespace coolsolve
