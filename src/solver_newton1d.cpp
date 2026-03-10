/**
 * @file solver_newton1d.cpp
 * @brief Specialized 1D root-finding solver for size-1 implicit blocks.
 *
 * Standard Newton with line search can fail when the initial guess is very
 * far from the solution (e.g. enthalpy equations where default guess is 1
 * but solution is ~1e5).  Newton1D uses a multi-phase approach:
 *
 *   Phase 1 — Trust-region Newton with adaptive radius and sign-change
 *             (bracket) detection for bisection fallback.
 *   Phase 2 — Multi-probe exploration: if Phase 1 stalls, evaluates the
 *             residual at ~900 probe points (log-spaced ±1e8, plus values
 *             near every external variable × {0.5, 0.9, …, 2.0}).  All
 *             sign changes are scored by midpoint residual to prefer true
 *             roots over poles.
 *   Phase 3 — Bisection + Newton hybrid: once a bracket is found,
 *             alternates bisection and Newton steps within the bracket
 *             for fast, guaranteed convergence.
 *   Phase 4 — Final Newton polish from the best point found, with
 *             relaxed tolerance acceptance.
 *
 * Falls through to the standard pipeline if all phases fail.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <sstream>

namespace coolsolve {

// ============================================================================
// Newton1D Solver Implementation
// ============================================================================

SolverStatus Newton1DSolver::solve(
    Eval1D& eval,
    double& x,
    const std::map<std::string, double>& externalVars,
    const SolverOptions& options,
    SolverTrace* trace,
    std::string* outErrorMessage)
{
    double xCur = x;
    double radius = std::max(std::abs(xCur) * 2.0, 100.0);
    bool converged = false;

    // Track bracket for bisection fallback
    bool hasBracket = false;
    double xLo = 0, xHi = 0, fLo = 0, fHi = 0;

    // --------------------------------------------------------------------
    // Phase 1: Newton with trust-region limiting
    // --------------------------------------------------------------------
    // Uses fewer iterations since Phase 2 probing is more effective
    // for problems where the initial guess is far from the solution.
    int phase1MaxIter = std::min(options.maxIterations, 50);
    for (int iter = 0; iter < phase1MaxIter && !converged; ++iter) {
        double f, j;
        try {
            auto [fv, jv] = eval(xCur);
            f = fv; j = jv;
        } catch (const std::exception& e) {
            if (isFatalEvaluationError(e.what())) {
                if (outErrorMessage) *outErrorMessage = e.what();
                return SolverStatus::EvaluationError;
            }
            // Recoverable evaluation failure — try reducing x toward zero
            xCur *= 0.5;
            continue;
        } catch (...) {
            xCur *= 0.5;
            continue;
        }

        if (trace) {
            SolverTrace::Iteration traceIter;
            traceIter.iter = iter;
            traceIter.residualNorm = std::abs(f);
            traceIter.stepNorm = 0.0;
            traceIter.lambda = 1.0;
            traceIter.x = {xCur};
            traceIter.residuals = {f};
            trace->iterations.push_back(traceIter);
        }

        if (std::abs(f) < options.tolerance) {
            converged = true;
            break;
        }

        // Update bracket
        if (!hasBracket) {
            if (iter == 0) {
                xLo = xCur; fLo = f;
            } else {
                if (f * fLo < 0) {
                    hasBracket = true;
                    xHi = xCur; fHi = f;
                    if (xLo > xHi) { std::swap(xLo, xHi); std::swap(fLo, fHi); }
                } else {
                    xLo = xCur; fLo = f;
                }
            }
        } else {
            // Keep bracket tight
            if (f * fLo < 0) { xHi = xCur; fHi = f; }
            else { xLo = xCur; fLo = f; }
            if (xLo > xHi) { std::swap(xLo, xHi); std::swap(fLo, fHi); }
        }

        double dx;
        if (hasBracket) {
            // Bisection-Newton hybrid (Illinois/Dekker-style)
            double xNewton = (std::abs(j) > 1e-30) ? xCur - f / j : (xLo + xHi) / 2.0;
            if (xNewton > xLo && xNewton < xHi) {
                dx = xNewton - xCur;
            } else {
                dx = (xLo + xHi) / 2.0 - xCur;
            }
        } else if (std::abs(j) > 1e-30) {
            // Newton step with trust-region clamping
            dx = -f / j;
            if (std::abs(dx) > radius) {
                dx = (dx > 0 ? 1.0 : -1.0) * radius;
            }
            radius = std::max(radius, std::abs(dx) * 2.0);
        } else {
            // Zero derivative — explore by stepping
            dx = radius;
            radius *= 2.0;
        }

        xCur += dx;
        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = std::abs(dx);
        }
    }

    // --------------------------------------------------------------------
    // Phase 2: Multi-probe exploration to find a sign change
    // --------------------------------------------------------------------
    // If Newton didn't converge and no bracket found, scan both positive
    // and negative ranges on a log scale, plus intermediate values that
    // may cross sign near poles.
    if (!converged && !hasBracket) {
        std::vector<double> probes;

        // Dense scan around zero
        for (double v : {0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0}) {
            probes.push_back(v);
            probes.push_back(-v);
        }
        // Log-spaced scan
        for (double scale = 100; scale <= 1e8; scale *= 2.0) {
            probes.push_back(scale);
            probes.push_back(-scale);
        }
        // KEY STRATEGY: Probe near external variable values.
        // Solutions of equations are typically at scales related to their inputs.
        // Without this, narrow sign-change regions near poles are missed.
        for (const auto& [name, val] : externalVars) {
            if (!std::isfinite(val) || val == 0.0) continue;
            for (double frac : {0.5, 0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 1.5, 2.0}) {
                probes.push_back(val * frac);
                if (val > 0) probes.push_back(-val * frac);
            }
        }

        // Evaluate all probes and record results for bracket detection
        struct ProbeResult { double x; double f; bool valid; };
        std::vector<ProbeResult> results;
        double bestF = 1e30;
        double bestX = xCur;

        for (double probe : probes) {
            try {
                auto [f, j] = eval(probe);
                if (!std::isfinite(f)) continue;
                results.push_back({probe, f, true});
                if (std::abs(f) < bestF) {
                    bestF = std::abs(f);
                    bestX = probe;
                }
                if (std::abs(f) < options.tolerance) {
                    xCur = probe;
                    converged = true;
                    break;
                }
            } catch (...) {
                results.push_back({probe, 0.0, false});
            }
        }

        // Collect ALL sign changes, pick the best one (smallest midpoint |f|)
        // to avoid locking onto poles instead of roots.
        if (!converged) {
            std::sort(results.begin(), results.end(),
                      [](const ProbeResult& a, const ProbeResult& b) { return a.x < b.x; });

            struct Bracket { double lo, flo, hi, fhi, midF; };
            std::vector<Bracket> brackets;

            for (size_t i = 0; i + 1 < results.size(); ++i) {
                if (!results[i].valid || !results[i+1].valid) continue;
                if (results[i].f * results[i+1].f < 0) {
                    double lo = results[i].x, hi = results[i+1].x;
                    double flo_v = results[i].f, fhi_v = results[i+1].f;
                    double mid = (lo + hi) / 2.0;
                    double fmid = 1e30;
                    try {
                        auto [fm, jm] = eval(mid);
                        fmid = std::abs(fm);
                    } catch (...) {
                        fmid = 1e30; // Evaluation failed — likely a pole
                    }
                    brackets.push_back({lo, flo_v, hi, fhi_v, fmid});
                }
            }

            if (!brackets.empty()) {
                auto& best = *std::min_element(brackets.begin(), brackets.end(),
                    [](const Bracket& a, const Bracket& b) { return a.midF < b.midF; });
                hasBracket = true;
                xLo = best.lo; fLo = best.flo;
                xHi = best.hi; fHi = best.fhi;
                if (xLo > xHi) { std::swap(xLo, xHi); std::swap(fLo, fHi); }
            }
        }

        if (!converged && !hasBracket) {
            xCur = bestX;
        }
    }

    // --------------------------------------------------------------------
    // Phase 3: Bisection + Newton hybrid within bracket
    // --------------------------------------------------------------------
    if (!converged && hasBracket) {
        for (int iter = 0; iter < options.maxIterations && !converged; ++iter) {
            double xMid = (xLo + xHi) / 2.0;
            try {
                auto [f, j] = eval(xMid);

                if (trace) {
                    SolverTrace::Iteration traceIter;
                    traceIter.iter = 1000 + iter; // Distinguish bisection iterations
                    traceIter.residualNorm = std::abs(f);
                    traceIter.stepNorm = (xHi - xLo) / 2.0;
                    traceIter.lambda = 1.0;
                    traceIter.x = {xMid};
                    traceIter.residuals = {f};
                    traceIter.detail = "       [bisection] bracket: ["
                        + std::to_string(xLo) + ", " + std::to_string(xHi) + "]\n";
                    trace->iterations.push_back(traceIter);
                }

                if (std::abs(f) < options.tolerance) {
                    xCur = xMid;
                    converged = true;
                    break;
                }

                // Try Newton step within bracket
                if (std::abs(j) > 1e-30) {
                    double xNewton = xMid - f / j;
                    if (xNewton > xLo && xNewton < xHi) {
                        try {
                            auto [fn, jn] = eval(xNewton);
                            if (std::abs(fn) < std::abs(f)) {
                                if (fn * fLo < 0) { xHi = xNewton; fHi = fn; }
                                else { xLo = xNewton; fLo = fn; }
                                continue;
                            }
                        } catch (...) {}
                    }
                }

                // Fall back to bisection
                if (f * fLo < 0) { xHi = xMid; fHi = f; }
                else { xLo = xMid; fLo = f; }

                if (xHi - xLo < options.stepTolerance) {
                    xCur = (xLo + xHi) / 2.0;
                    converged = true;
                }
            } catch (...) {
                // Evaluation failed at midpoint — narrow bracket from the other side
                xHi = xMid;
            }
        }
    }

    // --------------------------------------------------------------------
    // Phase 4: Final Newton polish from best position
    // --------------------------------------------------------------------
    if (!converged) {
        for (int iter = 0; iter < 50 && !converged; ++iter) {
            try {
                auto [f, j] = eval(xCur);
                if (std::abs(f) < options.tolerance) { converged = true; break; }
                if (std::abs(f) < options.lsRelaxedTolerance) { converged = true; break; }
                if (std::abs(j) < 1e-30) break;
                double dx = -f / j;
                double maxStep = std::max(std::abs(xCur) * 2.0, 1e6);
                if (std::abs(dx) > maxStep) dx = (dx > 0 ? 1.0 : -1.0) * maxStep;
                xCur += dx;
            } catch (...) {
                break;
            }
        }
    }

    if (converged) {
        x = xCur;
        if (trace) {
            trace->finalStatus = SolverStatus::Success;
        }
        return SolverStatus::Success;
    }

    return SolverStatus::MaxIterations;
}

// ============================================================================
// Simplified Newton1D for use inside symbolic reduction
// ============================================================================

SolverStatus Newton1DSolver::solveSimple(
    Eval1D& eval,
    double& x,
    const SolverOptions& options)
{
    double xCur = x;
    for (int iter = 0; iter < options.maxIterations; ++iter) {
        try {
            auto [f, j] = eval(xCur);
            if (std::abs(f) < options.tolerance) {
                x = xCur;
                return SolverStatus::Success;
            }
            if (std::abs(j) < 1e-30) break;
            double dx = -f / j;
            double maxStep = std::max(std::abs(xCur) * 2.0, 1e6);
            if (std::abs(dx) > maxStep) dx = (dx > 0 ? 1.0 : -1.0) * maxStep;
            xCur += dx;
        } catch (...) {
            break;
        }
    }
    return SolverStatus::MaxIterations;
}

}  // namespace coolsolve
