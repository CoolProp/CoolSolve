/**
 * @file solver_lm.cpp
 * @brief Levenberg-Marquardt solver (damped least-squares).
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
    double initialResidualNorm = 0.0;

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
                      << ", lambda=" << lambda << std::endl;

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

        // Solve (J^T J + lambda * D) dy = -J^T F  (Marquardt diagonal scaling)
        Eigen::MatrixXd JtJ = J.transpose() * J;
        Eigen::VectorXd JtF = J.transpose() * F;
        Eigen::VectorXd diag_JtJ = JtJ.diagonal();
        for (int i = 0; i < n; ++i)
            JtJ(i, i) += lambda * std::max(diag_JtJ(i), 1e-6);

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(JtJ);
        Eigen::VectorXd dy = qr.solve(-JtF);

        if (!dy.allFinite()) {
            lambda *= options.lmLambdaIncrease;
            if (lambda > options.lmMaxLambda) {
                if (trace) { trace->finalStatus = SolverStatus::SingularJacobian; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::SingularJacobian;
            }
            continue;
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
            if (lambda > options.lmMaxLambda) {
                if (trace) { trace->finalStatus = SolverStatus::EvaluationError; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                x = y.cwiseProduct(scale);
                return SolverStatus::EvaluationError;
            }
            continue;
        }

        // Gain ratio
        Eigen::VectorXd D = diag_JtJ.cwiseMax(1e-6);
        double predicted = 0.5 * dy.dot(lambda * D.cwiseProduct(dy) - JtF);
        double actual = phi_old - phi_new;
        double rho = (std::abs(predicted) > 1e-30) ? actual / predicted : 0.0;

        if (actual > 0) {
            // Accept
            y = y_new;
            double stepNorm = dy.lpNorm<Eigen::Infinity>();
            if (trace && !trace->iterations.empty())
                trace->iterations.back().stepNorm = stepNorm;

            if (rho > 0.75)
                lambda = std::max(lambda * options.lmLambdaDecrease, options.lmMinLambda);
            else if (rho < 0.25)
                lambda = std::min(lambda * options.lmLambdaIncrease, options.lmMaxLambda);

            if (stepNorm < options.stepTolerance) {
                double nn = F_new.lpNorm<Eigen::Infinity>();
                if (nn < options.tolerance * 100 || nn < options.lsRelaxedTolerance) {
                    if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                    x = y.cwiseProduct(scale);
                    return SolverStatus::Success;
                }
            }
        } else {
            lambda = std::min(lambda * options.lmLambdaIncrease, options.lmMaxLambda);
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
