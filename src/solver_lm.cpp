/**
 * @file solver_lm.cpp
 * @brief Levenberg-Marquardt solver with Nielsen's lambda adaptation,
 *        cumulative diagonal scaling, and optional geodesic acceleration.
 *
 * Uses a non-monotone criterion (Grippo et al. 1986) for step acceptance:
 * a step is accepted if φ(x_new) < max(recent φ values), rather than
 * requiring strict decrease.
 *
 * Improvements over the basic LM:
 * - **Nielsen's lambda adaptation**: Smoother lambda transitions using
 *   λ = λ * max(1/3, 1 - (2ρ-1)³) on acceptance, λ = λ*ν with ν doubling
 *   on consecutive rejections (Nielsen 1999, Madsen et al. 2004).
 * - **Cumulative diagonal scaling**: D_k = max(D_{k-1}, diag(J^T J)),
 *   preventing scale collapse when the Jacobian changes dramatically.
 * - **Geodesic acceleration** (Transtrum & Sethna 2012): Adds a second-order
 *   correction to the LM step by evaluating the directional second derivative
 *   of F along the velocity step.  Costs 1 extra residual evaluation per
 *   iteration but can halve the number of iterations on curved problems.
 *   Controlled by `lmGeodesicAcceleration` (default: true).
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace coolsolve {

// ============================================================================
// Levenberg-Marquardt Solver Implementation
// ============================================================================

Eigen::VectorXd LevenbergMarquardtSolver::computeScalingFactors(const Eigen::VectorXd& x) const {
    return coolsolve::computeScalingFactors(x);
}

SolverStatus LevenbergMarquardtSolver::solve(Problem& problem,
                                              Eigen::VectorXd& x,
                                              const SolverOptions& options,
                                              SolverTrace* trace,
                                              std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();

    if (trace) {
        if (trace->solverType.empty())
            trace->solverType = "LevenbergMarquardt";
        else if (trace->solverType.find("LevenbergMarquardt") == std::string::npos)
            trace->solverType += " -> LevenbergMarquardt";
    }

    const int n = problem.size;
    if (n <= 0 || x.size() != n) return SolverStatus::InvalidInput;

    Eigen::VectorXd scale = options.enableScaling
        ? computeScalingFactors(x) : Eigen::VectorXd::Ones(n);

    Eigen::VectorXd y = x.cwiseQuotient(scale);
    Eigen::VectorXd F(n), x_unscaled(n);
    Eigen::MatrixXd J(n, n), J_unscaled(n, n);
    double lambda = options.lmInitialLambda;
    double nu = 2.0;   // Nielsen's rejection multiplier (doubles on consecutive rejects)
    double initialResidualNorm = 0.0;

    // Cumulative diagonal scaling (Nielsen): D never decreases across iterations
    Eigen::VectorXd D_cumulative = Eigen::VectorXd::Constant(n, 1e-6);

    // Non-monotone merit history (Grippo et al. 1986)
    NonMonotoneHistory meritHistory(options.lsNonMonotoneMemory);

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "LM solver timed out";
            x = y.cwiseProduct(scale);
            return SolverStatus::EvaluationError;
        }
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (detailedError) *detailedError = "LevenbergMarquardt: cancelled";
            x = y.cwiseProduct(scale);
            return SolverStatus::MaxIterations;
        }

        try {
            x_unscaled = y.cwiseProduct(scale);
            problem.evaluate(x_unscaled, F, J_unscaled, true);
            J = J_unscaled * scale.asDiagonal();
        } catch (const std::exception&) { throw; }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        if (options.verbose)
            std::cout << "LM iter " << iter << ": ||F||=" << residualNorm
                      << ", lambda=" << lambda << ", nu=" << nu;

        // Track merit for non-monotone acceptance
        double phi_track = 0.5 * F.squaredNorm();
        meritHistory.push(phi_track);
        double refPhi = meritHistory.boundedRef(phi_track);

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
            ti.lambda = lambda;
            Eigen::VectorXd xu = y.cwiseProduct(scale);
            ti.x = std::vector<double>(xu.data(), xu.data() + xu.size());
            ti.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(ti);
        }

        // Convergence checks
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

        // Normal equations: (J^T J + lambda * D) dy = -J^T F
        Eigen::MatrixXd JtJ = J.transpose() * J;
        Eigen::VectorXd JtF = J.transpose() * F;

        // Cumulative Marquardt diagonal scaling (Nielsen): D_k = max(D_{k-1}, diag(J^T J))
        for (int i = 0; i < n; ++i)
            D_cumulative(i) = std::max(D_cumulative(i), JtJ(i, i));

        // Build damped system
        Eigen::MatrixXd A = JtJ;
        for (int i = 0; i < n; ++i)
            A(i, i) += lambda * std::max(D_cumulative(i), 1e-6);

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
        Eigen::VectorXd v = qr.solve(-JtF);  // "velocity" step

        if (!v.allFinite()) {
            lambda *= options.lmLambdaIncrease;
            nu = 2.0;
            if (lambda > options.lmMaxLambda) {
                if (trace) { trace->finalStatus = SolverStatus::SingularJacobian; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::SingularJacobian;
            }
            continue;
        }

        // --- Geodesic acceleration (Transtrum & Sethna 2012) ---
        // Compute second-order correction: a = -(J^TJ + λD)^{-1} J^T fvv
        // where fvv ≈ (F(x+hv) - F(x) - h*J*v) / (h²/2) is the directional
        // second derivative of F along v.
        Eigen::VectorXd dy = v;  // total step = v + 0.5*a
        if (options.lmGeodesicAcceleration && v.norm() > 1e-15) {
            double h = std::sqrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, y.norm()) /
                       std::max(1e-15, v.norm());
            Eigen::VectorXd y_pert = y + h * v;
            Eigen::VectorXd F_pert(n);
            Eigen::MatrixXd J_dummy;
            try {
                problem.evaluate(y_pert.cwiseProduct(scale), F_pert, J_dummy, false);
                // fvv = (F_pert - F - h*J*v) / (0.5*h*h)
                Eigen::VectorXd fvv = (F_pert - F - h * (J * v)) / (0.5 * h * h);
                Eigen::VectorXd a = qr.solve(-J.transpose() * fvv);
                if (a.allFinite()) {
                    // Only apply acceleration if it doesn't dominate the velocity step
                    // (avoids instability when the second-order term is unreliable)
                    double accelRatio = (0.5 * a).norm() / v.norm();
                    if (accelRatio < 1.0) {
                        dy = v + 0.5 * a;
                        if (options.verbose)
                            std::cout << "  LM geodesic: |a|/|v| = " << accelRatio << std::endl;
                    }
                }
            } catch (...) {
                // Perturbation failed — fall back to velocity-only step
            }
        }

        // Trial point
        Eigen::VectorXd y_new = y + dy;
        Eigen::VectorXd F_new(n);
        Eigen::MatrixXd Jd;
        double phi_old = 0.5 * F.squaredNorm();
        double phi_new;
        try {
            problem.evaluate(y_new.cwiseProduct(scale), F_new, Jd, false);
            phi_new = 0.5 * F_new.squaredNorm();
        } catch (...) {
            lambda *= options.lmLambdaIncrease;
            nu = 2.0;
            if (lambda > options.lmMaxLambda) {
                if (trace) { trace->finalStatus = SolverStatus::EvaluationError; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::EvaluationError;
            }
            continue;
        }

        // Gain ratio (based on current phi for lambda management)
        Eigen::VectorXd D = D_cumulative.cwiseMax(1e-6);
        double predicted = 0.5 * dy.dot(lambda * D.cwiseProduct(dy) - JtF);
        double actual = phi_old - phi_new;
        double rho = (std::abs(predicted) > 1e-30) ? actual / predicted : 0.0;

        // Non-monotone acceptance: accept if φ_new < max(recent history)
        if (phi_new < refPhi) {
            // Accept step
            y = y_new;
            double stepNorm = dy.lpNorm<Eigen::Infinity>();
            if (trace && !trace->iterations.empty())
                trace->iterations.back().stepNorm = stepNorm;

            // Nielsen's smooth lambda adaptation (Madsen et al. 2004):
            //   λ = λ * max(1/3, 1 - (2ρ - 1)³), ν = 2
            if (options.lmNielsenUpdate) {
                double tmp = 2.0 * rho - 1.0;
                double factor = std::max(1.0 / 3.0, 1.0 - tmp * tmp * tmp);
                lambda = std::max(lambda * factor, options.lmMinLambda);
                nu = 2.0;
            } else {
                // Legacy behavior
                if (rho > 0.75)
                    lambda = std::max(lambda * options.lmLambdaDecrease, options.lmMinLambda);
                else if (rho < 0.25)
                    lambda = std::min(lambda * options.lmLambdaIncrease, options.lmMaxLambda);
            }

            if (stepNorm < options.stepTolerance) {
                double nn = F_new.lpNorm<Eigen::Infinity>();
                if (nn < options.tolerance * 100 || nn < options.lsRelaxedTolerance) {
                    if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                    x = y.cwiseProduct(scale);
                    return SolverStatus::Success;
                }
            }
        } else {
            // Reject step — Nielsen's exponential increase: λ *= ν, ν *= 2
            if (options.lmNielsenUpdate) {
                lambda = std::min(lambda * nu, options.lmMaxLambda);
                nu *= 2.0;
            } else {
                lambda = std::min(lambda * options.lmLambdaIncrease, options.lmMaxLambda);
            }
        }
    }

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        std::ostringstream ss;
        ss << "Levenberg-Marquardt: Max iterations (" << options.maxIterations << ") reached.";
        *detailedError = ss.str();
    }
    x = y.cwiseProduct(scale);
    return SolverStatus::MaxIterations;
}

}  // namespace coolsolve
