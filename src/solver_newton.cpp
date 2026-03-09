/**
 * @file solver_newton.cpp
 * @brief Newton solver with non-monotone backtracking line search and optional
 *        Broyden quasi-Newton Jacobian reuse.
 *
 * Uses the Grippo-Lampariello-Lucidi (1986) non-monotone line search:
 * instead of requiring strict decrease at every step (standard Armijo),
 * the merit function is compared against the maximum over the last M
 * iterations.  This helps escape narrow curved valleys and saddle points
 * that trap monotone methods.  When M=1, behaviour is identical to the
 * classic monotone Armijo line search.
 *
 * When broydenRecomputeInterval > 0, the full Jacobian is computed only
 * every K iterations; intermediate iterations use a Broyden Type-I rank-1
 * update (Broyden 1965), saving O(n) CoolProp evaluations per iteration.
 * If a Broyden step fails the line search, the full Jacobian is
 * automatically recomputed and the step retried.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace coolsolve {

// ============================================================================
// Newton Solver Implementation
// ============================================================================

Eigen::VectorXd NewtonSolver::computeScalingFactors(const Eigen::VectorXd& x) const {
    return coolsolve::computeScalingFactors(x);
}

double NewtonSolver::lineSearch(Problem& problem,
                                const Eigen::VectorXd& x,
                                const Eigen::VectorXd& dx,
                                const Eigen::VectorXd& F,
                                const SolverOptions& options,
                                double refPhi) {
    double phi0 = 0.5 * F.squaredNorm();
    // Directional derivative: ∇φ·dx = -||F||² = -2φ₀ for the full Newton step
    double dphi0 = -2.0 * phi0;

    double lambda = 1.0;
    Eigen::VectorXd x_new(x.size());
    Eigen::VectorXd F_new(x.size());
    Eigen::MatrixXd J_dummy;

    for (int lsIter = 0; lsIter < options.lsMaxIterations; ++lsIter) {
        x_new = x + lambda * dx;
        try {
            problem.evaluate(x_new, F_new, J_dummy, false);
        } catch (const std::exception&) {
            lambda *= options.lsRho;
            continue;
        }

        double phi_new = 0.5 * F_new.squaredNorm();

        // Non-monotone Armijo condition (Grippo et al. 1986):
        //   φ(x + λd) ≤ refPhi + α λ ∇φ·d
        // where refPhi = max(φ_{k-M+1}, ..., φ_k) from the outer iteration.
        // When M=1, refPhi = φ₀ and this reduces to the standard monotone Armijo.
        if (phi_new <= refPhi + options.lsAlpha * lambda * dphi0) return lambda;
        // Near-minimum relaxation
        if (phi0 < 0.05 && phi_new < phi0 * 2.0) return lambda;
        if (lambda < 0.01 && phi0 < 0.1 && phi_new <= phi0 * 1.5) return lambda;

        lambda *= options.lsRho;
        if (lambda < options.lsMinStep) return 0.0;
    }
    return 0.0;
}

SolverStatus NewtonSolver::solve(Problem& problem,
                                 Eigen::VectorXd& x,
                                 const SolverOptions& options,
                                 SolverTrace* trace,
                                 std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();

    if (trace) {
        if (trace->solverType.empty())
            trace->solverType = "Newton";
        else if (trace->solverType.find("Newton") == std::string::npos)
            trace->solverType += " -> Newton";
    }

    const int n = problem.size;
    if (n <= 0 || x.size() != n) return SolverStatus::InvalidInput;

    // Scaling
    Eigen::VectorXd scale = options.enableScaling
        ? computeScalingFactors(x) : Eigen::VectorXd::Ones(n);
    if (options.verbose && options.enableScaling) {
        std::cout << "Newton: Scaling factors min=" << scale.minCoeff()
                  << ", max=" << scale.maxCoeff() << std::endl;
    }

    Eigen::VectorXd y = x.cwiseQuotient(scale);
    Eigen::VectorXd F(n), dy(n), x_unscaled(n);
    Eigen::MatrixXd J(n, n), J_unscaled(n, n);
    double initialResidualNorm = 0.0;

    // Non-monotone merit history (Grippo et al. 1986)
    NonMonotoneHistory meritHistory(options.lsNonMonotoneMemory);

    // --- Broyden quasi-Newton state ---
    const int broydenK = options.broydenRecomputeInterval;
    const bool useBroyden = (broydenK > 0 && n > 1);
    Eigen::MatrixXd B;            // Broyden Jacobian approximation (scaled coords)
    Eigen::VectorXd F_prev;       // Previous residual (for rank-1 update)
    Eigen::VectorXd y_prev;       // Previous iterate in scaled coords
    int itersSinceFullJ = 0;      // Count iterations since last full Jacobian
    bool forceFullJacobian = false; // Retry flag when Broyden step fails

    if (useBroyden && options.verbose) {
        std::cout << "Newton: Broyden recompute interval K=" << broydenK << std::endl;
    }

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "Solver timed out";
            return SolverStatus::EvaluationError;
        }
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (detailedError) *detailedError = "Newton: cancelled";
            return SolverStatus::MaxIterations;
        }

        // --- Evaluate F and J (full or Broyden update) ---
        bool fullJacobianThisIter = !useBroyden || iter == 0 ||
                                    itersSinceFullJ >= broydenK ||
                                    forceFullJacobian;

        try {
            x_unscaled = y.cwiseProduct(scale);
            if (fullJacobianThisIter) {
                // Full Jacobian evaluation
                problem.evaluate(x_unscaled, F, J_unscaled, true);
                J = J_unscaled * scale.asDiagonal();
                if (useBroyden) {
                    B = J;
                    itersSinceFullJ = 0;
                    forceFullJacobian = false;
                }
            } else {
                // Residual-only evaluation + Broyden rank-1 update
                Eigen::MatrixXd J_dummy;
                problem.evaluate(x_unscaled, F, J_dummy, false);

                // Broyden Type-I update: B += (dF - B*dy) * dy^T / (dy^T * dy)
                Eigen::VectorXd dF = F - F_prev;
                Eigen::VectorXd dy_actual = y - y_prev;
                double dy2 = dy_actual.squaredNorm();
                if (dy2 > 1e-30) {
                    B += ((dF - B * dy_actual) * dy_actual.transpose()) / dy2;
                }
                J = B;
                itersSinceFullJ++;
            }
        } catch (const std::exception&) {
            throw;
        }

        // Save state for Broyden update on next iteration
        if (useBroyden) {
            F_prev = F;
            y_prev = y;
        }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        // Track merit for non-monotone line search
        double phi = 0.5 * F.squaredNorm();
        meritHistory.push(phi);
        double refPhi = meritHistory.boundedRef(phi);

        if (options.verbose) {
            std::cout << "Newton iter " << iter << ": ||F||_inf = " << residualNorm;
            if (useBroyden)
                std::cout << (fullJacobianThisIter ? " [full J]" : " [Broyden]");
            if (options.lsNonMonotoneMemory > 1)
                std::cout << " (refPhi=" << refPhi << ", M=" << meritHistory.size() << ")";
            std::cout << std::endl;
        }

        // Record trace
        if (trace) {
            SolverTrace::Iteration ti;
            ti.iter = iter;
            ti.residualNorm = residualNorm;
            ti.stepNorm = 0.0;
            ti.lambda = 1.0;
            Eigen::VectorXd xu = y.cwiseProduct(scale);
            ti.x = std::vector<double>(xu.data(), xu.data() + xu.size());
            ti.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(ti);
        }

        // Convergence: absolute
        if (residualNorm < options.tolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }
        // Convergence: relative
        if (initialResidualNorm > 0 && residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }

        // Solve J*dy = -F
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(J);
        dy = qr.solve(-F);

        if (!dy.allFinite()) {
            // If using Broyden approximation, retry with full Jacobian
            if (useBroyden && !fullJacobianThisIter) {
                forceFullJacobian = true;
                continue;
            }
            if (trace) {
                trace->finalStatus = SolverStatus::SingularJacobian;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                trace->singularJacobianF.assign(F.data(), F.data() + F.size());
                trace->singularJacobianJ.resize(J.rows());
                for (Eigen::Index r = 0; r < J.rows(); ++r)
                    trace->singularJacobianJ[r].assign(J.row(r).data(), J.row(r).data() + J.cols());
            }
            return SolverStatus::SingularJacobian;
        }

        // Line search in scaled coordinates with non-monotone reference
        Problem lsProblem;
        lsProblem.size = n;
        lsProblem.evaluate = [&](const Eigen::VectorXd& yt, Eigen::VectorXd& Fo, Eigen::MatrixXd& Jdummy, bool) {
            problem.evaluate(yt.cwiseProduct(scale), Fo, Jdummy, false);
        };
        double lambda = lineSearch(lsProblem, y, dy, F, options, refPhi);

        if (lambda == 0.0) {
            // If Broyden step failed line search, retry with full Jacobian
            if (useBroyden && !fullJacobianThisIter) {
                forceFullJacobian = true;
                // Restore F_prev/y_prev so next Broyden update is correct
                if (!F_prev.size()) { F_prev = F; y_prev = y; }
                continue;
            }

            if (residualNorm < options.lsRelaxedTolerance) {
                if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::Success;
            }
            // Fallback: tiny step
            for (double fb : {0.1, 0.01, 0.001}) {
                Eigen::VectorXd yt = y + fb * dy;
                Eigen::VectorXd Ft(n);
                Eigen::MatrixXd Jd;
                try {
                    problem.evaluate(yt.cwiseProduct(scale), Ft, Jd, false);
                    if (0.5 * Ft.squaredNorm() < 0.5 * F.squaredNorm()) { lambda = fb; break; }
                } catch (...) { continue; }
            }
            if (lambda == 0.0) {
                if (trace) { trace->finalStatus = SolverStatus::LineSearchFailed; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                if (detailedError) {
                    std::ostringstream ss;
                    ss << "Line search failed at iteration " << iter << ".\n"
                       << "Residual norm ||F|| = " << residualNorm << "\n"
                       << "Newton step norm ||dy|| = " << dy.norm() << "\n";
                    *detailedError = ss.str();
                }
                x = y.cwiseProduct(scale);
                return SolverStatus::LineSearchFailed;
            }
        }

        double stepNorm = (lambda * dy).lpNorm<Eigen::Infinity>();
        y += lambda * dy;
        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = stepNorm;
            trace->iterations.back().lambda = lambda;
        }

        if (stepNorm < options.stepTolerance && residualNorm < options.tolerance * 100) {
            if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }
    }

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        std::ostringstream ss;
        ss << "Max iterations (" << options.maxIterations << ") reached. Last ||F|| = " << F.lpNorm<Eigen::Infinity>();
        *detailedError = ss.str();
    }
    x = y.cwiseProduct(scale);
    return SolverStatus::MaxIterations;
}

}  // namespace coolsolve
