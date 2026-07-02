/**
 * @file solver_trust_region.cpp
 * @brief Trust-region dogleg solver with non-monotone step acceptance and
 *        adaptive initial radius.
 *
 * Uses a non-monotone criterion (Grippo et al. 1986) for step acceptance:
 * a step is accepted if φ(x_new) < max(recent φ values), rather than
 * requiring φ(x_new) < φ(x_current).  The trust region radius management
 * still uses the standard gain ratio based on the current iterate.
 *
 * Improvements over the basic dogleg:
 * - **Adaptive initial radius**: When trAdaptiveRadius is enabled, the
 *   initial delta is set to min(||Cauchy step||, trInitialRadius) on the
 *   first iteration so that it matches the problem's natural scale.
 * - **Smoother radius management**: Uses the gain ratio rho to smoothly
 *   shrink or grow delta, rather than binary accept/reject.  On rejection,
 *   uses max(0.1, 1-rho)*delta instead of a fixed trShrinkFactor.
 * - **Better rejection recovery**: Consecutive rejections trigger a Cauchy
 *   step fallback (gradient direction within delta) before resetting.
 *
 * Optional hybrd-style Broyden updates (trBroydenRecomputeInterval > 0):
 * mirrors MINPACK's `hybrd` (Powell 1970), which combines Broyden rank-1
 * Jacobian updates with trust-region dogleg step acceptance. Unlike the
 * Newton solver's dense Broyden update, the approximation here is
 * maintained as an incrementally-updated QR factorization (see
 * solver_hybrd_qr.h, based on Golub & Van Loan, "Matrix Computations"
 * Sec. 12.5) applied via numerically stable Givens rotations, then
 * reconstructed to a dense Jacobian for the existing (unchanged) dogleg/
 * model/acceptance logic below. A full Jacobian is recomputed every K
 * iterations, or immediately after trBroydenRestartRejects consecutive
 * rejected steps (Powell's restart criterion), or if the resulting Newton
 * step is non-finite. Disabled by default (K=0), leaving the legacy
 * behavior byte-for-byte unchanged.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include "coolsolve/solver_hybrd_qr.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace coolsolve {

// ============================================================================
// Trust Region Dogleg Solver Implementation
// ============================================================================

Eigen::VectorXd TrustRegionSolver::computeScalingFactors(const Eigen::VectorXd& x) const {
    return coolsolve::computeScalingFactors(x);
}

double TrustRegionSolver::evaluateModel(const Eigen::VectorXd& F,
                                        const Eigen::MatrixXd& J,
                                        const Eigen::VectorXd& p) {
    Eigen::VectorXd Fp = F + J * p;
    return 0.5 * Fp.squaredNorm();
}

Eigen::VectorXd TrustRegionSolver::doglegStep(const Eigen::VectorXd& dx_n,
                                              const Eigen::VectorXd& dx_c,
                                              double delta) {
    double norm_n = dx_n.norm();
    double norm_c = dx_c.norm();

    if (norm_n <= delta) return dx_n;
    if (norm_c >= delta) return (delta / norm_c) * dx_c;

    // Dogleg interpolation
    Eigen::VectorXd diff = dx_n - dx_c;
    double a = diff.squaredNorm();
    double b = 2.0 * dx_c.dot(diff);
    double c = dx_c.squaredNorm() - delta * delta;
    double disc = std::max(0.0, b * b - 4.0 * a * c);
    double tau = std::max(0.0, std::min(1.0, (-b + std::sqrt(disc)) / (2.0 * a)));
    return dx_c + tau * diff;
}

SolverStatus TrustRegionSolver::solve(Problem& problem,
                                      Eigen::VectorXd& x,
                                      const SolverOptions& options,
                                      SolverTrace* trace,
                                      std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();

    if (trace) {
        if (trace->solverType.empty())
            trace->solverType = "TrustRegion";
        else if (trace->solverType.find("TrustRegion") == std::string::npos)
            trace->solverType += " -> TrustRegion";
    }

    const int n = problem.size;
    if (n <= 0 || x.size() != n) return SolverStatus::InvalidInput;

    Eigen::VectorXd scale = options.enableScaling
        ? computeScalingFactors(x) : Eigen::VectorXd::Ones(n);

    Eigen::VectorXd y = x.cwiseQuotient(scale);
    Eigen::VectorXd F(n);
    Eigen::MatrixXd J(n, n), J_unscaled(n, n);
    double initialResidualNorm = 0.0;
    double delta = options.trInitialRadius;
    bool deltaInitialized = !options.trAdaptiveRadius; // defer if adaptive
    int consecutiveRejects = 0;

    // Non-monotone merit history (Grippo et al. 1986)
    NonMonotoneHistory meritHistory(options.lsNonMonotoneMemory);

    // --- hybrd-style Broyden state (QR-maintained; see file header) ---
    const int broydenK = options.trBroydenRecomputeInterval;
    const bool useBroyden = (broydenK > 0 && n > 1);
    Eigen::MatrixXd Q, R;          // QR factors of the Broyden Jacobian approximation
    Eigen::VectorXd F_prev;        // Previous residual (for rank-1 update)
    Eigen::VectorXd y_prev;        // Previous iterate in scaled coords
    int itersSinceFullJ = 0;       // Count iterations since last full Jacobian
    bool forceFullJacobian = false; // Retry flag when Broyden step is unusable

    if (useBroyden && options.verbose) {
        std::cout << "TrustRegion: Broyden recompute interval K=" << broydenK
                   << ", restart after " << options.trBroydenRestartRejects
                   << " rejects" << std::endl;
    }

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "Solver timed out";
            x = y.cwiseProduct(scale);
            return SolverStatus::EvaluationError;
        }
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (detailedError) *detailedError = "TrustRegion: cancelled";
            x = y.cwiseProduct(scale);
            return SolverStatus::MaxIterations;
        }

        bool fullJacobianThisIter = !useBroyden || iter == 0 ||
                                    itersSinceFullJ >= broydenK ||
                                    forceFullJacobian;

        try {
            Eigen::VectorXd xu = y.cwiseProduct(scale);
            if (fullJacobianThisIter) {
                problem.evaluate(xu, F, J_unscaled, true);
                J = J_unscaled * scale.asDiagonal();
                if (useBroyden) {
                    computeInitialQR(J, Q, R);
                    itersSinceFullJ = 0;
                    forceFullJacobian = false;
                }
            } else {
                // Residual-only evaluation + Broyden rank-1 update via QR
                Eigen::MatrixXd J_dummy;
                problem.evaluate(xu, F, J_dummy, false);

                Eigen::VectorXd dy_actual = y - y_prev;
                double dy2 = dy_actual.squaredNorm();
                if (dy2 > 1e-30) {
                    // Broyden Type-I update B += (dF - B*dy)*dy^T/dy2, computed
                    // via the QR factors: B*dy = Q*(R*dy).
                    Eigen::VectorXd Bdy = Q * (R * dy_actual);
                    Eigen::VectorXd u = (F - F_prev - Bdy) / dy2;
                    rank1QRUpdate(Q, R, u, dy_actual);
                }
                J = Q * R;
                itersSinceFullJ++;
            }
        } catch (const std::exception&) {
            x = y.cwiseProduct(scale);
            throw;
        }

        if (useBroyden) {
            F_prev = F;
            y_prev = y;
        }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        if (options.verbose)
            std::cout << "TrustRegion iter " << iter << ": ||F||=" << residualNorm
                      << ", delta=" << delta;

        // Track merit for non-monotone acceptance
        double phi_current = 0.5 * F.squaredNorm();
        meritHistory.push(phi_current);
        double refPhi = meritHistory.boundedRef(phi_current);

        if (options.verbose) {
            if (options.lsNonMonotoneMemory > 1)
                std::cout << " (refPhi=" << refPhi << ", M=" << meritHistory.size() << ")";
            std::cout << std::endl;
        }

        if (trace) {
            SolverTrace::Iteration ti;
            ti.iter = iter;
            ti.residualNorm = residualNorm;
            ti.stepNorm = 0.0;
            ti.lambda = delta;
            Eigen::VectorXd xu = y.cwiseProduct(scale);
            ti.x = std::vector<double>(xu.data(), xu.data() + xu.size());
            ti.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(ti);
        }

        if (residualNorm < options.tolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }
        if (initialResidualNorm > 0 && residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }

        // Gradient g = J^T * F
        Eigen::VectorXd g = J.transpose() * F;
        double gNormSq = g.squaredNorm();
        if (gNormSq < 1e-30) {
            if (residualNorm < options.lsRelaxedTolerance) {
                if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::Success;
            }
            if (trace) {
                trace->finalStatus = SolverStatus::SingularJacobian;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                trace->singularJacobianF.assign(F.data(), F.data() + F.size());
                trace->singularJacobianJ.resize(J.rows());
                for (Eigen::Index r = 0; r < J.rows(); ++r)
                    trace->singularJacobianJ[r].assign(J.row(r).data(), J.row(r).data() + J.cols());
            }
            x = y.cwiseProduct(scale);
            return SolverStatus::SingularJacobian;
        }

        // Cauchy step
        Eigen::VectorXd Jg = J * g;
        double alpha = gNormSq / Jg.squaredNorm();
        Eigen::VectorXd dx_c = -alpha * g;

        // Adaptive initial radius: set delta based on Cauchy step norm
        if (!deltaInitialized) {
            double cauchyNorm = dx_c.norm();
            if (cauchyNorm > 0) {
                delta = std::min(options.trInitialRadius, std::max(cauchyNorm, 1.0));
            }
            deltaInitialized = true;
            if (options.verbose)
                std::cout << "TrustRegion: adaptive initial delta = " << delta << std::endl;
        }

        // Newton step
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(J);
        Eigen::VectorXd dx_n = qr.solve(-F);
        if (!dx_n.allFinite()) {
            // Broyden approximation produced an unusable step; retry this
            // iterate with a full Jacobian instead of the Cauchy fallback.
            if (useBroyden && !fullJacobianThisIter) {
                forceFullJacobian = true;
                continue;
            }
            dx_n = dx_c;
        }

        // Dogleg step
        Eigen::VectorXd dy = doglegStep(dx_n, dx_c, delta);

        // Try step
        Eigen::VectorXd y_new = y + dy;
        Eigen::VectorXd F_new(n);
        Eigen::MatrixXd Jd;

        try {
            problem.evaluate(y_new.cwiseProduct(scale), F_new, Jd, false);
        } catch (const std::exception&) {
            // Evaluation failed — smoothly shrink delta
            delta *= options.trShrinkFactor;
            consecutiveRejects++;
            if (consecutiveRejects > 20) {
                delta = std::max(delta, options.trInitialRadius * 0.1);
                consecutiveRejects = 0;
            }
            if (delta < 1e-12) {
                if (trace) { trace->finalStatus = SolverStatus::LineSearchFailed; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::LineSearchFailed;
            }
            continue;
        }

        double phi_old = 0.5 * F.squaredNorm();
        double phi_new = 0.5 * F_new.squaredNorm();
        double model_old = evaluateModel(F, J, Eigen::VectorXd::Zero(n));
        double model_new = evaluateModel(F, J, dy);
        double actual = phi_old - phi_new;
        double predicted = model_old - model_new;
        double rho = (std::abs(predicted) > 1e-30) ? actual / predicted : 0.0;

        // Non-monotone acceptance: accept if φ_new < max(recent history)
        // The gain ratio rho is still based on the current phi for radius management.
        if (phi_new < refPhi) {
            // Accept step
            y = y_new;
            consecutiveRejects = 0;
            double stepNorm = dy.lpNorm<Eigen::Infinity>();
            if (trace && !trace->iterations.empty())
                trace->iterations.back().stepNorm = stepNorm;

            if (stepNorm < options.stepTolerance) {
                double newNorm = F_new.lpNorm<Eigen::Infinity>();
                if (newNorm < options.tolerance * 100 || newNorm < options.lsRelaxedTolerance) {
                    if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                    x = y.cwiseProduct(scale);
                    return SolverStatus::Success;
                }
            }

            // Smooth rho-based radius management (replaces binary grow/no-grow)
            if (rho > 0.75 && dy.norm() >= 0.9 * delta)
                delta = std::min(options.trGrowFactor * delta, options.trMaxRadius);
            else if (rho > 0.5)
                delta = std::min(1.5 * delta, options.trMaxRadius);
            // Good rho but step well inside delta: leave delta unchanged
        } else {
            // Reject step — rho-based shrinking instead of fixed factor
            consecutiveRejects++;
            if (useBroyden && consecutiveRejects >= options.trBroydenRestartRejects) {
                // Powell restart criterion: repeated rejections suggest the
                // Broyden approximation has drifted; refresh from a full
                // Jacobian on the next iteration.
                forceFullJacobian = true;
            }
            if (rho < 0.0) {
                // Very bad step: aggressive shrink
                delta *= std::max(0.1, options.trShrinkFactor);
            } else {
                // Moderate shrink proportional to how bad the step was
                double shrink = std::max(0.25, 1.0 - (1.0 - rho) * 0.5);
                delta *= std::min(shrink, options.trShrinkFactor);
            }
            if (consecutiveRejects > 15) {
                // On sustained rejection, try gradient direction at small radius
                double gradScale = std::sqrt(gNormSq);
                if (gradScale > 1e-15) {
                    delta = std::max(delta, 0.01 * gradScale);
                } else {
                    delta = std::max(delta, options.trInitialRadius * 0.01);
                }
                consecutiveRejects = 0;
            }
            if (delta < 1e-12) {
                if (trace) { trace->finalStatus = SolverStatus::LineSearchFailed; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::LineSearchFailed;
            }
        }
    }

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        std::ostringstream ss;
        ss << "Trust region: Max iterations (" << options.maxIterations << ") reached.";
        *detailedError = ss.str();
    }
    x = y.cwiseProduct(scale);
    return SolverStatus::MaxIterations;
}

}  // namespace coolsolve
