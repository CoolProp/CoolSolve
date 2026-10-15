/**
 * @file solver_newton1d.cpp
 * @brief Specialized 1D root-finding solver for size-1 implicit blocks.
 *
 * Standard Newton with line search can fail when the initial guess is very
 * far from the solution (e.g. enthalpy equations where default guess is 1
 * but solution is ~1e5).  Newton1D uses a multi-phase approach:
 *
 *   Phase 1 — Trust-region Newton with adaptive radius and sign-change
 *             (bracket) detection for bisection fallback.  A Newton step
 *             that increases |F| is rejected and retried with half the step.
 *   Phase 2 — Multi-probe exploration: if Phase 1 stalls, evaluates the
 *             residual at ~900 probe points (log-spaced ±1e8, plus values
 *             near every external variable × {0.5, 0.9, …, 2.0}).  All
 *             sign changes are scored by midpoint residual to prefer true
 *             roots over poles.
 *   Phase 3 — Bisection + Newton hybrid: once a bracket is found,
 *             alternates bisection and Newton steps within the bracket
 *             for fast, guaranteed convergence.
 *   Phase 4 — Final Newton polish from the best point found (same step
 *             rejection as Phase 1), with relaxed tolerance acceptance.
 *
 * Falls through to the standard pipeline if all phases fail.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace coolsolve {

namespace {

/**
 * Step control of the Newton iterations of Phases 1 and 4: a Newton step that
 * increases |F| (or reaches a point where F is not finite) is rejected and
 * retried from the last accepted iterate with half the step.  Without it, the
 * step of an equation that is hyperbolic in its unknown (f = c - k/x from a
 * guess above the root) jumps over the pole at x = 0, |F| grows and the
 * iterate runs away.  A step that does not come from a Newton direction
 * (zero derivative, bisection) is never rejected.
 */
struct StepControl {
    /// Successive halvings of one step after which the Newton direction is given up.
    static constexpr int kMaxRejections = 30;

    double xAcc = 0.0;     ///< last accepted iterate
    double fAcc = std::numeric_limits<double>::infinity();  ///< F at the last accepted iterate
    bool newtonTrial = false;  ///< the current iterate is a Newton step taken from xAcc
    int rejections = 0;    ///< successive rejections of the current step

    /// If the iterate \p x with residual \p f, reached by a Newton step, increased |F|: replace
    /// \p x by the point half-way back to xAcc and return true.  A step that changes the sign of F
    /// is kept when \p keepSignChange (it brackets a root: Phase 1 bisects it).
    bool reject(double& x, double f, bool keepSignChange) {
        if (!newtonTrial) return false;
        if (std::isfinite(f) && (std::abs(f) <= std::abs(fAcc) || (keepSignChange && f * fAcc < 0))) return false;
        x = xAcc + 0.5 * (x - xAcc);
        ++rejections;
        return true;
    }

    /// Accept the iterate \p x (residual \p f); \p newtonStep tells whether the step taken from it is a Newton step.
    void accept(double x, double f, bool newtonStep) {
        xAcc = x;
        fAcc = f;
        newtonTrial = newtonStep;
        rejections = 0;
    }

    bool exhausted() const { return rejections > kMaxRejections; }
};

}  // namespace

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
    StepControl stepControl;

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
            traceIter.lambda = std::ldexp(1.0, -stepControl.rejections);  // share of the Newton step tried
            traceIter.x = {xCur};
            traceIter.residuals = {f};
            trace->iterations.push_back(traceIter);
        }

        if (std::abs(f) < options.tolerance) {
            converged = true;
            break;
        }

        // Step control: a Newton step that increased |F| is retried, shortened, from the
        // last accepted iterate (the trust radius shrinks with it).  A rejected point does
        // not enter the bracket: a retry across a pole would pair it with a root-free interval.
        if (!hasBracket && stepControl.reject(xCur, f, true)) {
            if (stepControl.exhausted()) {
                xCur = stepControl.xAcc;
                break;  // the Newton direction gives no descent: let Phase 2 explore
            }
            radius = std::abs(xCur - stepControl.xAcc);
            if (trace && !trace->iterations.empty()) {
                trace->iterations.back().stepNorm = radius;
            }
            continue;
        }

        // Update bracket (a pole where F is not finite is not a sign change)
        if (!std::isfinite(f)) {
            // nothing to record
        } else if (!hasBracket) {
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
        bool newtonStep = false;
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
            newtonStep = true;
        } else {
            // Zero derivative — explore by stepping
            dx = radius;
            radius *= 2.0;
        }

        stepControl.accept(xCur, f, newtonStep);
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
        StepControl polishControl;
        for (int iter = 0; iter < 50 && !converged; ++iter) {
            try {
                auto [f, j] = eval(xCur);
                if (std::abs(f) < options.tolerance) { converged = true; break; }
                if (polishControl.reject(xCur, f, false)) {
                    if (polishControl.exhausted()) break;
                    continue;
                }
                const bool newtonStep = std::abs(j) >= 1e-30;
                polishControl.accept(xCur, f, newtonStep);
                if (!newtonStep) break;
                double dx = -f / j;
                double maxStep = std::max(std::abs(xCur) * 2.0, 1e6);
                if (std::abs(dx) > maxStep) dx = (dx > 0 ? 1.0 : -1.0) * maxStep;
                xCur += dx;
            } catch (...) {
                break;
            }
        }
        // The polish cannot go further: accept the best iterate if its residual is small
        if (!converged && std::abs(polishControl.fAcc) < options.lsRelaxedTolerance) {
            xCur = polishControl.xAcc;
            converged = true;
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
