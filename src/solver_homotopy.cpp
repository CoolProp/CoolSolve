/**
 * @file solver_homotopy.cpp
 * @brief Homotopy continuation solver for global convergence.
 *
 * Constructs a homotopy H(x,t) = t·F(x) + (1-t)·(x - x0) and tracks
 * the solution curve from t=0 (trivial: x=x0) to t=1 (target: F(x)=0).
 *
 * Uses Newton as the corrector at each continuation step with adaptive
 * step-size control: dt increases when Newton converges quickly and
 * decreases when it fails or converges slowly.
 *
 * Theoretical guarantee (Chow-Mallet-Paret-Yorke, 1978): for almost all
 * starting points x0, the homotopy path is a smooth curve that reaches
 * a solution of F(x) = 0.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace coolsolve {

SolverStatus HomotopySolver::solve(Problem& problem,
                                    Eigen::VectorXd& x,
                                    const SolverOptions& options,
                                    SolverTrace* trace,
                                    std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();
    const int n = problem.size;

    if (trace) {
        if (trace->solverType.empty()) trace->solverType = "Homotopy";
        else trace->solverType += " -> Homotopy";
    }

    if (n <= 0 || x.size() != n) return SolverStatus::InvalidInput;

    // x0 is our starting point (the initial guess)
    Eigen::VectorXd x0 = x;

    // Adaptive continuation parameters
    double t = 0.0;
    double dt = 0.1;               // initial step size
    const double dtMin = 1e-6;     // minimum step size before giving up
    const double dtMax = 0.3;      // maximum step size
    const int maxSteps = 200;      // max continuation steps
    const int correctorMaxIter = 30; // Newton iterations per step

    // Corrector tolerances (tighter at higher t)
    auto correctorTol = [&](double t_val) {
        return options.tolerance * (1.0 + 100.0 * (1.0 - t_val));
    };

    for (int step = 0; step < maxSteps; ++step) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "Homotopy: timed out";
            break;
        }

        double t_new = std::min(t + dt, 1.0);

        // Build homotopy problem: H(x, t_new) = t_new * F(x) + (1 - t_new) * (x - x0) = 0
        NonLinearSolver::Problem hProblem;
        hProblem.size = n;
        hProblem.evaluate = [&](const Eigen::VectorXd& xv,
                                Eigen::VectorXd& H,
                                Eigen::MatrixXd& Jh,
                                bool computeJac) {
            Eigen::VectorXd F(n);
            Eigen::MatrixXd J(n, n);
            problem.evaluate(xv, F, J, computeJac);

            // H(x) = t_new * F(x) + (1 - t_new) * (x - x0)
            H = t_new * F + (1.0 - t_new) * (xv - x0);

            if (computeJac) {
                // dH/dx = t_new * J + (1 - t_new) * I
                Jh = t_new * J;
                for (int i = 0; i < n; ++i)
                    Jh(i, i) += (1.0 - t_new);
            }
        };

        // Solve H(x, t_new) = 0 using Newton with reduced iterations
        SolverOptions correctorOpts;
        correctorOpts.maxIterations = correctorMaxIter;
        correctorOpts.tolerance = correctorTol(t_new);
        correctorOpts.relativeTolerance = 1e-10;
        correctorOpts.stepTolerance = options.stepTolerance;
        correctorOpts.enableScaling = options.enableScaling;
        correctorOpts.lsAlpha = options.lsAlpha;
        correctorOpts.lsRho = options.lsRho;
        correctorOpts.lsMaxIterations = options.lsMaxIterations;
        correctorOpts.lsMinStep = options.lsMinStep;
        correctorOpts.lsRelaxedTolerance = correctorTol(t_new) * 10.0;
        correctorOpts.verbose = false;

        Eigen::VectorXd x_trial = x;  // warm-start from previous solution
        NewtonSolver corrector;
        SolverStatus corrStatus = corrector.solve(hProblem, x_trial, correctorOpts);

        if (trace) {
            SolverTrace::Iteration ti;
            ti.iter = step;
            // Evaluate actual F at current x for reporting
            Eigen::VectorXd Factual(n);
            Eigen::MatrixXd Jd;
            try { problem.evaluate(x_trial, Factual, Jd, false); }
            catch (...) { Factual.setConstant(1e30); }
            ti.residualNorm = Factual.lpNorm<Eigen::Infinity>();
            ti.stepNorm = dt;
            ti.lambda = t_new;
            ti.x = std::vector<double>(x_trial.data(), x_trial.data() + n);
            ti.residuals = std::vector<double>(Factual.data(), Factual.data() + n);
            ti.detail = "       [homotopy t=" + std::to_string(t_new) + ", dt=" + std::to_string(dt) + "]\n";
            trace->iterations.push_back(ti);
        }

        if (corrStatus == SolverStatus::Success) {
            // Corrector converged — accept step
            x = x_trial;
            t = t_new;

            if (t >= 1.0) {
                // Verify actual F(x) is small at t=1
                Eigen::VectorXd Fcheck(n);
                Eigen::MatrixXd Jd;
                try {
                    problem.evaluate(x, Fcheck, Jd, false);
                    if (Fcheck.lpNorm<Eigen::Infinity>() < options.tolerance) {
                        if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                        return SolverStatus::Success;
                    }
                    // Close but not converged — do a final Newton polish
                    Eigen::VectorXd x_save = x;
                    SolverOptions polishOpts = options;
                    polishOpts.maxIterations = 50;
                    NewtonSolver polish;
                    SolverStatus polishStatus = polish.solve(problem, x, polishOpts);
                    if (polishStatus == SolverStatus::Success) {
                        if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                        return SolverStatus::Success;
                    }
                    x = x_save;  // restore if polish diverged
                } catch (...) {}

                // Homotopy reached t=1 but final solution isn't precise enough
                // Return the best we have
                if (trace) { trace->finalStatus = SolverStatus::MaxIterations; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                return SolverStatus::MaxIterations;
            }

            // Increase step size on success (but cap at dtMax)
            dt = std::min(dt * 1.5, dtMax);
        } else {
            // Corrector failed — reduce step size and retry
            dt *= 0.5;
            if (dt < dtMin) {
                if (options.verbose)
                    std::cout << "Homotopy: step size too small (dt=" << dt
                              << ") at t=" << t << std::endl;
                // Try once more with LM as corrector before giving up
                Eigen::VectorXd x_lm = x;
                LevenbergMarquardtSolver lmCorrector;
                SolverOptions lmOpts = correctorOpts;
                lmOpts.maxIterations = 100;
                SolverStatus lmStatus = lmCorrector.solve(hProblem, x_lm, lmOpts);
                if (lmStatus == SolverStatus::Success) {
                    x = x_lm;
                    t = t_new;
                    dt = dtMin * 10;  // reset to a small but usable step
                    continue;
                }
                break;  // give up
            }
        }
    }

    // Reached here means homotopy didn't fully converge.
    // Try a final Newton polish on the actual problem from where we are.
    {
        Eigen::VectorXd x_save = x;
        NewtonSolver finalNewton;
        SolverOptions finalOpts = options;
        finalOpts.maxIterations = 50;
        SolverStatus finalStatus = finalNewton.solve(problem, x, finalOpts);
        if (finalStatus == SolverStatus::Success) {
            if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }
        x = x_save;  // restore if polish diverged
    }

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        std::ostringstream ss;
        ss << "Homotopy: did not converge. Last t=" << t << ", dt=" << dt;
        *detailedError = ss.str();
    }
    return SolverStatus::MaxIterations;
}

}  // namespace coolsolve
