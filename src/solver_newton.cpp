/**
 * @file solver_newton.cpp
 * @brief Newton solver with backtracking line search.
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
                                const SolverOptions& options) {
    double phi0 = 0.5 * F.squaredNorm();
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

        // Armijo condition
        if (phi_new <= phi0 + options.lsAlpha * lambda * dphi0) return lambda;
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

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "Solver timed out";
            return SolverStatus::EvaluationError;
        }

        // Evaluate F(x), J(x) in original coordinates, then scale Jacobian
        try {
            x_unscaled = y.cwiseProduct(scale);
            problem.evaluate(x_unscaled, F, J_unscaled, true);
            J = J_unscaled * scale.asDiagonal();
        } catch (const std::exception&) {
            throw;
        }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        if (options.verbose)
            std::cout << "Newton iter " << iter << ": ||F||_inf = " << residualNorm << std::endl;

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

        // Line search in scaled coordinates
        Problem lsProblem;
        lsProblem.size = n;
        lsProblem.evaluate = [&](const Eigen::VectorXd& yt, Eigen::VectorXd& Fo, Eigen::MatrixXd& Jdummy, bool) {
            problem.evaluate(yt.cwiseProduct(scale), Fo, Jdummy, false);
        };
        double lambda = lineSearch(lsProblem, y, dy, F, options);

        if (lambda == 0.0) {
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
