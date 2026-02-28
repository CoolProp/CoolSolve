/**
 * @file solver_trust_region.cpp
 * @brief Trust-region dogleg solver.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
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
    int consecutiveRejects = 0;  // track consecutive rejections for adaptive reset

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "Solver timed out";
            x = y.cwiseProduct(scale);
            return SolverStatus::EvaluationError;
        }

        try {
            Eigen::VectorXd xu = y.cwiseProduct(scale);
            problem.evaluate(xu, F, J_unscaled, true);
            J = J_unscaled * scale.asDiagonal();
        } catch (const std::exception&) {
            x = y.cwiseProduct(scale);
            throw;
        }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        if (options.verbose)
            std::cout << "TrustRegion iter " << iter << ": ||F||=" << residualNorm
                      << ", delta=" << delta << std::endl;

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

        // Newton step
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(J);
        Eigen::VectorXd dx_n = qr.solve(-F);
        if (!dx_n.allFinite()) dx_n = dx_c;

        // Dogleg step
        Eigen::VectorXd dy = doglegStep(dx_n, dx_c, delta);

        // Try step
        Eigen::VectorXd y_new = y + dy;
        Eigen::VectorXd F_new(n);
        Eigen::MatrixXd Jd;

        try {
            problem.evaluate(y_new.cwiseProduct(scale), F_new, Jd, false);
        } catch (const std::exception&) {
            // Evaluation failed — smoothly shrink delta instead of hard reset
            delta *= options.trShrinkFactor;
            consecutiveRejects++;
            if (consecutiveRejects > 20) {
                // Gradual reset: try a moderate radius instead of full trInitialRadius
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

        if (actual > 0) {
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

            // Grow delta on good steps
            if (rho > 0.75 && dy.norm() >= 0.9 * delta)
                delta = std::min(options.trGrowFactor * delta, options.trMaxRadius);
        } else {
            // Reject step — smooth shrink
            consecutiveRejects++;
            delta *= options.trShrinkFactor;
            if (consecutiveRejects > 20) {
                delta = std::max(delta, options.trInitialRadius * 0.1);
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
