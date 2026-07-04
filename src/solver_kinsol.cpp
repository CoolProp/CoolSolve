/**
 * @file solver_kinsol.cpp
 * @brief SUNDIALS-KINSOL-style nonlinear solver with three globalisation
 *        strategies: inexact Newton + Dennis-Schnabel line search (KIN_LINESEARCH),
 *        Picard / Richardson fixed-point iteration (KIN_PICARD), and
 *        Anderson-accelerated fixed-point iteration (KIN_FP).
 *
 * This is a faithful in-tree port of the SUNDIALS KINSOL globalisation modes.
 * No external library is linked: the algorithms are implemented directly on
 * top of Eigen, following the same conventions (variable scaling, traces,
 * non-monotone-free sufficient decrease) as the other CoolSolve solvers.
 *
 * References:
 *  - Hindmarsh, Brown, Lee, Serban, Shumaker, Woodward (2005): *SUNDIALS:
 *    Suite of Nonlinear and Differential/Algebraic Equation Solvers*,
 *    ACM Trans. Math. Softw. 31(3).  KINSOL is the nonlinear solver component.
 *  - Dennis, J. E., Schnabel, R. B. (1983/1996): *Numerical Methods for
 *    Unconstrained Optimization and Nonlinear Equations*, SIAM.  The line
 *    search is Algorithm A6.3.1mod (backtracking with quadratic then cubic
 *    interpolation and the Armijo sufficient-decrease test on φ = ½||F||²).
 *  - Anderson, D. G. (1965): *Iterative Procedures for Nonlinear Integral
 *    Equations*, J. ACM 12(4).  Origin of Anderson acceleration.
 *  - Walker, H. F., Ni, P. (2011): *Anderson Acceleration for Fixed-Point
 *    Iterations*, SIAM J. Numer. Anal. 49(4).  Reference for the
 *    least-squares depth-m Anderson implementation used in KIN_FP.
 *
 * The problem is posed as F(x) = 0.  Picard and FixedPoint modes solve the
 * equivalent fixed-point form G(x) = x − ωF(x) = 0 (ω = kinsolPicardOmega,
 * default 1).  Convergence is declared on ||F||_∞ < tolerance for all modes.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <Eigen/Dense>
#include <sstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

namespace coolsolve {

// ============================================================================
// KINSOL Solver Implementation
// ============================================================================

Eigen::VectorXd KINSOLSolver::computeScalingFactors(const Eigen::VectorXd& x) const {
    return coolsolve::computeScalingFactors(x);
}

// ----------------------------------------------------------------------------
// Dennis-Schnabel backtracking line search on φ(λ) = ½||F(y + λ·s)||²
// (Dennis & Schnabel 1983, Algorithm A6.3.1mod).
//
//   φ0    = φ(0)        = ½||F(y)||²
//   φ'0   = dφ/dλ|₀     = Fᵀ·(Js) = −||F||²  for the Newton step (Js = −F)
//
// Acceptance test (Armijo sufficient decrease):
//   φ(λ) ≤ φ0 + α·λ·φ'0,    α = kinsolLineSearchAlpha (default 1e-4).
// First backtrack uses quadratic interpolation of (0,φ0,φ'0) & (λ,φ(λ));
// subsequent backtracks use cubic interpolation through the two most recent
// points and the origin.  Each interpolated λ is safeguarded to lie in
// [0.1·λ_old, 0.5·λ_old] so progress is guaranteed.
// Returns the accepted λ (or 0 on failure) and writes F(y+λs) to F_out.
// ----------------------------------------------------------------------------

double KINSOLSolver::lineSearch(Problem& problem,
                                const Eigen::VectorXd& y,
                                const Eigen::VectorXd& s,
                                const Eigen::VectorXd& scale,
                                double phi0,
                                double dphi0,
                                const SolverOptions& options,
                                Eigen::VectorXd& F_out) {
    const int n = y.size();
    const double alpha = options.kinsolLineSearchAlpha;
    const int maxSteps = std::max(1, options.kinsolLineSearchMaxIters);
    const double minStep = std::max(options.lsMinStep, 1e-12);

    Eigen::VectorXd xTrial(n), FTrial(n);
    Eigen::MatrixXd Jdummy;

    double lambda = 1.0;                 // try the full Newton step first
    double lambdaPrev = lambda;          // previous (rejected) λ for cubic
    double phiPrev = 0.0;                // φ at λPrev
    bool havePrev = false;
    double bestLambda = 0.0, bestPhi = std::numeric_limits<double>::max();
    Eigen::VectorXd bestF = F_out;

    for (int step = 0; step < maxSteps; ++step) {
        if (lambda < minStep) break;
        xTrial = y + lambda * s;
        try {
            problem.evaluate(xTrial.cwiseProduct(scale), FTrial, Jdummy, false);
        } catch (const std::exception&) {
            // Evaluation failed (e.g. CoolProp penalty region): shrink and retry.
            lambdaPrev = lambda;
            havePrev = true;
            lambda *= 0.5;
            continue;
        }
        double phi = 0.5 * FTrial.squaredNorm();
        if (!std::isfinite(phi)) {
            lambdaPrev = lambda;
            havePrev = true;
            lambda *= 0.5;
            continue;
        }

        // Track best feasible point seen (for graceful failure).
        if (phi < bestPhi) { bestPhi = phi; bestLambda = lambda; bestF = FTrial; }

        // Armijo sufficient-decrease test.
        if (phi <= phi0 + alpha * lambda * dphi0) {
            F_out = FTrial;
            return lambda;
        }

        // Choose the next λ by interpolation.
        double lambdaNew;
        if (!havePrev) {
            // Quadratic interpolation through (0,φ0,φ'0) and (λ,φ):
            //   q(λ) = φ0 + φ'0·λ + c·λ²,   c = (φ − φ0 − φ'0·λ)/λ²
            //   minimiser: λ* = −φ'0·λ² / (2(φ − φ0 − φ'0·λ))
            double denom = 2.0 * (phi - phi0 - dphi0 * lambda);
            if (std::abs(denom) < 1e-300) {
                lambdaNew = 0.5 * lambda;
            } else {
                lambdaNew = -dphi0 * lambda * lambda / denom;
            }
        } else {
            // Cubic interpolation through (0,φ0,φ'0), (λPrev,φPrev), (λ,φ).
            // φ(λ) − φ0 − φ'0·λ = a·λ³ + b·λ²  ⇒  a·λ + b = u/λ²  with
            //   u = φ − φ0 − φ'0·λ.  Solve 2×2 for (a,b), then minimise.
            double u_prev = phiPrev - phi0 - dphi0 * lambdaPrev;
            double u_curr = phi - phi0 - dphi0 * lambda;
            double a_cubic, b_cubic;
            if (std::abs(lambdaPrev) < 1e-300 || std::abs(lambda) < 1e-300) {
                lambdaNew = 0.5 * lambda;
            } else {
                double r_prev = u_prev / (lambdaPrev * lambdaPrev);
                double r_curr = u_curr / (lambda * lambda);
                double ddiv = (lambdaPrev - lambda);
                if (std::abs(ddiv) < 1e-300) {
                    lambdaNew = 0.5 * lambda;
                } else {
                    a_cubic = (r_prev - r_curr) / ddiv;
                    b_cubic = r_prev - a_cubic * lambdaPrev;
                    // Minimise ½(3a·λ² + 2b·λ + φ'0) → λ = (−b + √(b²−3a·φ'0))/(3a)  if a≠0
                    double disc = b_cubic * b_cubic - 3.0 * a_cubic * dphi0;
                    if (std::abs(a_cubic) < 1e-300 || disc < 0.0) {
                        // Cubic unusable: fall back to quadratic using the origin
                        // and the current point.
                        double denom2 = 2.0 * u_curr;
                        lambdaNew = (std::abs(denom2) < 1e-300)
                                        ? 0.5 * lambda
                                        : -dphi0 * lambda * lambda / denom2;
                    } else {
                        lambdaNew = (-b_cubic + std::sqrt(disc)) / (3.0 * a_cubic);
                    }
                }
            }
        }

        // Safeguard: ensure λ strictly decreases and stays in [0.1·λ, 0.5·λ].
        if (!std::isfinite(lambdaNew) || lambdaNew > 0.5 * lambda || lambdaNew < 0.1 * lambda) {
            lambdaNew = 0.5 * lambda;
        }

        lambdaPrev = lambda;
        phiPrev = phi;
        havePrev = true;
        lambda = lambdaNew;
    }

    // No point satisfied the sufficient-decrease test.  Return the best point
    // if it at least improved on φ0 (helps the outer loop accept a usable step);
    // otherwise signal failure with λ = 0.
    if (bestLambda > 0.0 && bestPhi < phi0) {
        F_out = bestF;
        return bestLambda;
    }
    return 0.0;
}

// ----------------------------------------------------------------------------
// KIN_LINESEARCH: inexact Newton (exact direct linear solve) + line search.
// ----------------------------------------------------------------------------

SolverStatus KINSOLSolver::solveLineSearch(Problem& problem, Eigen::VectorXd& y,
                                           const Eigen::VectorXd& scale,
                                           const SolverOptions& options,
                                           SolverTrace* trace,
                                           std::string* detailedError) {
    const int n = y.size();
    Eigen::VectorXd F(n), dy(n);
    Eigen::MatrixXd J(n, n), J_unscaled(n, n);
    double initialResidualNorm = 0.0;
    auto startTime = std::chrono::high_resolution_clock::now();

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "KINSOL: timed out";
            return SolverStatus::EvaluationError;
        }
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (detailedError) *detailedError = "KINSOL: cancelled";
            return SolverStatus::MaxIterations;
        }

        Eigen::VectorXd xu = y.cwiseProduct(scale);
        try {
            problem.evaluate(xu, F, J_unscaled, true);
        } catch (const std::exception& e) {
            if (detailedError) *detailedError = std::string("KINSOL: ") + e.what();
            return SolverStatus::EvaluationError;
        }
        J = J_unscaled * scale.asDiagonal();

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        if (options.verbose) {
            std::cout << "KINSOL(LS) iter " << iter << ": ||F||_inf = " << residualNorm << std::endl;
        }
        if (trace) {
            SolverTrace::Iteration ti;
            ti.iter = iter;
            ti.residualNorm = residualNorm;
            ti.stepNorm = 0.0;
            ti.lambda = 1.0;
            Eigen::VectorXd xtr = y.cwiseProduct(scale);
            ti.x = std::vector<double>(xtr.data(), xtr.data() + xtr.size());
            ti.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(ti);
        }

        // Convergence: absolute.
        if (residualNorm < options.tolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }
        // Convergence: relative.
        if (initialResidualNorm > 0 &&
            residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }

        // Exact linear solve: J · dy = −F (inexact-Newton with η = 0, since
        // CoolSolve uses direct solvers rather than Krylov iterations).
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(J);
        dy = qr.solve(-F);
        if (!dy.allFinite() || qr.rank() < n) {
            if (trace) { trace->finalStatus = SolverStatus::SingularJacobian;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            if (detailedError) {
                std::ostringstream ss;
                ss << "KINSOL(LS): singular Jacobian at iteration " << iter
                   << " (rank " << qr.rank() << "/" << n << "), ||F|| = " << residualNorm;
                *detailedError = ss.str();
            }
            return SolverStatus::SingularJacobian;
        }

        // Line search.  φ'₀ = Fᵀ·(J·dy) = Fᵀ·(−F) = −||F||².
        // lineSearch works in scaled coordinates: it receives the real
        // `problem` (which takes unscaled x) and applies `scale` itself once
        // when evaluating each trial point.
        double phi0 = 0.5 * F.squaredNorm();
        double dphi0 = -F.squaredNorm();
        Eigen::VectorXd F_new = F;
        double lambda = lineSearch(problem, y, dy, scale, phi0, dphi0, options, F_new);

        if (lambda == 0.0) {
            // Line search could not find sufficient decrease.
            if (residualNorm < options.lsRelaxedTolerance) {
                if (trace) { trace->finalStatus = SolverStatus::Success;
                             trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                return SolverStatus::Success;
            }
            if (trace) { trace->finalStatus = SolverStatus::LineSearchFailed;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            if (detailedError) {
                std::ostringstream ss;
                ss << "KINSOL(LS): line search failed at iteration " << iter
                   << ", ||F|| = " << residualNorm;
                *detailedError = ss.str();
            }
            return SolverStatus::LineSearchFailed;
        }

        double stepNorm = (lambda * dy).lpNorm<Eigen::Infinity>();
        y += lambda * dy;
        F = F_new;
        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = stepNorm;
            trace->iterations.back().lambda = lambda;
        }

        if (stepNorm < options.stepTolerance && residualNorm < options.tolerance * 100) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }
    }

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations;
                 trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        std::ostringstream ss;
        ss << "KINSOL(LS): max iterations (" << options.maxIterations
           << ") reached, ||F|| = " << F.lpNorm<Eigen::Infinity>();
        *detailedError = ss.str();
    }
    return SolverStatus::MaxIterations;
}

// ----------------------------------------------------------------------------
// KIN_PICARD: fixed-point (Richardson) iteration  y ← y − ω·F(y).
// ----------------------------------------------------------------------------

SolverStatus KINSOLSolver::solvePicard(Problem& problem, Eigen::VectorXd& y,
                                       const Eigen::VectorXd& scale,
                                       const SolverOptions& options,
                                       SolverTrace* trace,
                                       std::string* detailedError) {
    const int n = y.size();
    const double omega = options.kinsolPicardOmega;
    Eigen::VectorXd F(n);
    double initialResidualNorm = 0.0;
    auto startTime = std::chrono::high_resolution_clock::now();

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "KINSOL(Picard): timed out";
            return SolverStatus::EvaluationError;
        }
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (detailedError) *detailedError = "KINSOL(Picard): cancelled";
            return SolverStatus::MaxIterations;
        }

        Eigen::VectorXd xu = y.cwiseProduct(scale);
        Eigen::MatrixXd Jdummy;
        try {
            problem.evaluate(xu, F, Jdummy, false);
        } catch (const std::exception& e) {
            if (detailedError) *detailedError = std::string("KINSOL(Picard): ") + e.what();
            return SolverStatus::EvaluationError;
        }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        if (options.verbose) {
            std::cout << "KINSOL(Picard) iter " << iter << ": ||F||_inf = " << residualNorm << std::endl;
        }
        if (trace) {
            SolverTrace::Iteration ti;
            ti.iter = iter;
            ti.residualNorm = residualNorm;
            ti.stepNorm = 0.0;
            ti.lambda = omega;
            Eigen::VectorXd xtr = y.cwiseProduct(scale);
            ti.x = std::vector<double>(xtr.data(), xtr.data() + xtr.size());
            ti.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(ti);
        }

        if (!std::isfinite(residualNorm)) {
            if (trace) { trace->finalStatus = SolverStatus::Diverged;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            if (detailedError) *detailedError = "KINSOL(Picard): residual became non-finite";
            return SolverStatus::Diverged;
        }
        if (residualNorm < options.tolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }
        if (initialResidualNorm > 0 &&
            residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }

        // Fixed-point (Richardson) update: y_{k+1} = y_k − ω·F(y_k).
        Eigen::VectorXd step = -omega * F;
        double stepNorm = step.lpNorm<Eigen::Infinity>();
        y += step;
        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = stepNorm;
        }

        if (stepNorm < options.stepTolerance && residualNorm < options.tolerance * 100) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }
    }

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations;
                 trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        std::ostringstream ss;
        ss << "KINSOL(Picard): max iterations (" << options.maxIterations
           << ") reached, ||F|| = " << F.lpNorm<Eigen::Infinity>();
        *detailedError = ss.str();
    }
    return SolverStatus::MaxIterations;
}

// ----------------------------------------------------------------------------
// KIN_FP: fixed-point iteration with Anderson acceleration (Anderson 1965,
// Walker & Ni 2011).  G(y) = y − ω·F(y), ω = kinsolPicardOmega.
//
// History of the last (m+1) iterates is kept.  At each step:
//   f_i = G(y_i) − y_i = −ω·F(y_i)             (fixed-point residual)
//   ΔF[:,j] = f_{j+1} − f_j,  j = 0..m_k−1
//   γ = argmin_γ || f_{m_k} − ΔF·γ ||₂          (least squares, col-pivoted QR)
//   β_0 = γ_0 ; β_j = γ_j − γ_{j−1} ; β_{m_k} = 1 − γ_{m_k−1}   (so Σβ = 1)
//   y_next = (1−θ)·y_{m_k} + θ·Σ_i β_i·G(y_i)
// ----------------------------------------------------------------------------

SolverStatus KINSOLSolver::solveFixedPoint(Problem& problem, Eigen::VectorXd& y,
                                           const Eigen::VectorXd& scale,
                                           const SolverOptions& options,
                                           SolverTrace* trace,
                                           std::string* detailedError) {
    const int n = y.size();
    const int depth = std::max(0, options.kinsolAndersonDepth);
    const double theta = options.kinsolAndersonRelaxation;
    const double omega = options.kinsolPicardOmega;

    auto startTime = std::chrono::high_resolution_clock::now();
    double initialResidualNorm = 0.0;

    // Evaluate G(y_0) and store y_0, G(y_0).
    auto evalG = [&](const Eigen::VectorXd& yv, Eigen::VectorXd& Gout, double& resNorm)
                     -> bool {
        Eigen::VectorXd xu = yv.cwiseProduct(scale);
        Eigen::VectorXd Fv(n);
        Eigen::MatrixXd Jdummy;
        try {
            problem.evaluate(xu, Fv, Jdummy, false);
        } catch (const std::exception&) {
            return false;
        }
        resNorm = Fv.lpNorm<Eigen::Infinity>();
        Gout = yv - omega * Fv;     // G(y) = y − ω·F(y)
        return true;
    };

    std::vector<Eigen::VectorXd> Ys;   // iterates y_i
    std::vector<Eigen::VectorXd> Gs;   // G(y_i)

    Eigen::VectorXd G0(n);
    double res0 = 0.0;
    if (!evalG(y, G0, res0)) {
        if (detailedError) *detailedError = "KINSOL(FP): initial evaluation failed";
        return SolverStatus::EvaluationError;
    }
    initialResidualNorm = res0;
    Ys.push_back(y);
    Gs.push_back(G0);

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "KINSOL(FP): timed out";
            return SolverStatus::EvaluationError;
        }
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (detailedError) *detailedError = "KINSOL(FP): cancelled";
            return SolverStatus::MaxIterations;
        }

        int k = static_cast<int>(Ys.size()) - 1;     // index of the current iterate
        double residualNorm = (Ys[k] - Gs[k]).lpNorm<Eigen::Infinity>() / std::max(omega, 1e-300);
        // (||F|| = ||y − G(y)||/ω ; for ω=1 this is just ||Ys−Gs||.)

        if (options.verbose) {
            std::cout << "KINSOL(FP) iter " << iter << ": ||F||_inf = " << residualNorm << std::endl;
        }
        if (trace) {
            SolverTrace::Iteration ti;
            ti.iter = iter;
            ti.residualNorm = residualNorm;
            ti.stepNorm = 0.0;
            ti.lambda = theta;
            Eigen::VectorXd xtr = Ys[k].cwiseProduct(scale);
            ti.x = std::vector<double>(xtr.data(), xtr.data() + xtr.size());
            // Residuals = ω·F = (y − G(y)).
            Eigen::VectorXd Fk = Ys[k] - Gs[k];
            ti.residuals = std::vector<double>(Fk.data(), Fk.data() + Fk.size());
            trace->iterations.push_back(ti);
        }

        if (!std::isfinite(residualNorm)) {
            if (trace) { trace->finalStatus = SolverStatus::Diverged;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            if (detailedError) *detailedError = "KINSOL(FP): residual became non-finite";
            return SolverStatus::Diverged;
        }
        if (residualNorm < options.tolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            y = Ys[k];
            return SolverStatus::Success;
        }
        if (initialResidualNorm > 0 &&
            residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) { trace->finalStatus = SolverStatus::Success;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            y = Ys[k];
            return SolverStatus::Success;
        }

        // Number of usable differences = min(current index k, depth).
        int mk = std::min(k, depth);
        Eigen::VectorXd yNext(n);
        bool usedAnderson = false;

        if (mk >= 1) {
            // Build ΔF (n × mk): column j = f_{j+1} − f_j = (G_{j+1}−G_j) − (Y_{j+1}−Y_j).
            Eigen::MatrixXd dF(n, mk);
            for (int j = 0; j < mk; ++j) {
                Eigen::VectorXd fNext = Gs[j + 1] - Ys[j + 1];   // f_{j+1} = −ω·F_{j+1}
                Eigen::VectorXd fCurr = Gs[j] - Ys[j];           // f_j
                dF.col(j) = fNext - fCurr;
            }
            Eigen::VectorXd fk = Gs[k] - Ys[k];                  // current residual

            // Least-squares: γ = argmin || fk − ΔF·γ ||.  Use col-pivoted QR
            // (rank-revealing) so rank-deficient ΔF degrades gracefully.
            Eigen::VectorXd gamma = dF.colPivHouseholderQr().solve(fk);
            if (gamma.allFinite()) {
                // β weights: β_0=γ_0, β_j=γ_j−γ_{j−1}, β_{mk}=1−γ_{mk−1}.
                Eigen::VectorXd yAnderson = Eigen::VectorXd::Zero(n);
                double beta0 = gamma(0);
                yAnderson += beta0 * Gs[0];
                for (int j = 1; j < mk; ++j) {
                    double bj = gamma(j) - gamma(j - 1);
                    yAnderson += bj * Gs[j];
                }
                double betaLast = 1.0 - gamma(mk - 1);
                yAnderson += betaLast * Gs[mk];     // Σ_i β_i·G(y_i), i = 0..mk
                // Damping: y_next = (1−θ)·y_current + θ·y_Anderson.
                yNext = (1.0 - theta) * Ys[k] + theta * yAnderson;
                usedAnderson = true;
            }
        }

        if (!usedAnderson) {
            // Plain fixed-point step (mk == 0 or LS failed): y_next = G(y_current).
            yNext = Gs[k];
        }

        // Evaluate G at the new iterate for the next iteration.
        Eigen::VectorXd Gnew(n);
        double resNew = 0.0;
        if (!evalG(yNext, Gnew, resNew)) {
            if (detailedError) *detailedError = "KINSOL(FP): evaluation failed during iteration";
            if (trace) { trace->finalStatus = SolverStatus::EvaluationError;
                         trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            y = Ys[k];
            return SolverStatus::EvaluationError;
        }

        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = (yNext - Ys[k]).lpNorm<Eigen::Infinity>();
        }

        // Append to history and trim to (depth + 1) most recent pairs.
        Ys.push_back(yNext);
        Gs.push_back(Gnew);
        while (static_cast<int>(Ys.size()) > depth + 1) {
            Ys.erase(Ys.begin());
            Gs.erase(Gs.begin());
        }
        y = yNext;
    }

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations;
                 trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        int k = static_cast<int>(Ys.size()) - 1;
        double rn = (Ys[k] - Gs[k]).lpNorm<Eigen::Infinity>() / std::max(omega, 1e-300);
        std::ostringstream ss;
        ss << "KINSOL(FP): max iterations (" << options.maxIterations << ") reached, ||F|| = " << rn;
        *detailedError = ss.str();
    }
    return SolverStatus::MaxIterations;
}

// ----------------------------------------------------------------------------
// Top-level dispatch.
// ----------------------------------------------------------------------------

SolverStatus KINSOLSolver::solve(Problem& problem,
                                 Eigen::VectorXd& x_guess,
                                 const SolverOptions& options,
                                 SolverTrace* trace,
                                 std::string* detailedError) {
    if (trace && trace->solverType.empty()) trace->solverType = "KINSOL";

    const int n = problem.size;
    if (n <= 0 || x_guess.size() != n) return SolverStatus::InvalidInput;

    Eigen::VectorXd scale = options.enableScaling
        ? computeScalingFactors(x_guess) : Eigen::VectorXd::Ones(n);
    if (options.verbose && options.enableScaling) {
        std::cout << "KINSOL: Scaling factors min=" << scale.minCoeff()
                  << ", max=" << scale.maxCoeff() << std::endl;
    }

    Eigen::VectorXd y = x_guess.cwiseQuotient(scale);   // work in scaled coordinates

    SolverStatus status;
    switch (options.kinsolGlobalStrategy) {
        case KinsolGlobalStrategy::LineSearch:
            status = solveLineSearch(problem, y, scale, options, trace, detailedError);
            break;
        case KinsolGlobalStrategy::Picard:
            status = solvePicard(problem, y, scale, options, trace, detailedError);
            break;
        case KinsolGlobalStrategy::FixedPoint:
            status = solveFixedPoint(problem, y, scale, options, trace, detailedError);
            break;
        default:
            if (detailedError) *detailedError = "KINSOL: unknown global strategy";
            return SolverStatus::InvalidInput;
    }

    // Write the final iterate back into the caller's vector (in unscaled
    // coordinates), so the pipeline can warm-start the next solver regardless
    // of whether this one converged.
    x_guess = y.cwiseProduct(scale);
    return status;
}

}  // namespace coolsolve
