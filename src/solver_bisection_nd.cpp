/**
 * @file solver_bisection_nd.cpp
 * @brief Multi-dimensional bisection solver for small systems.
 *
 * Implements a sign-change–based bisection approach inspired by Vrahatis (1988).
 * Only applicable to systems with n ≤ MAX_BISECTION_DIM (typically 8-10).
 *
 * Strategy:
 *   1. Probe a structured grid of points around the initial guess to find
 *      a simplex whose vertices have distinct sign patterns in F(x).
 *   2. Repeatedly bisect the largest edge of the simplex, replacing the
 *      vertex whose sign pattern matches the midpoint's.
 *   3. Stop when the simplex diameter is below tolerance or F is small enough.
 *
 * This is a last-resort solver: expensive but guaranteed to converge if a
 * simplex with proper sign patterns is found.
 */
#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>

namespace coolsolve {

/// Default maximum dimension.  The effective limit is read from SolverOptions::bisectionNDMaxBlockSize.
static constexpr int DEFAULT_MAX_BISECTION_DIM = 8;

/**
 * @brief Compute the sign pattern vector of F(x): +1 or -1 per component.
 */
static std::vector<int> signPattern(const Eigen::VectorXd& F) {
    std::vector<int> s(F.size());
    for (int i = 0; i < F.size(); ++i)
        s[i] = (F(i) >= 0) ? 1 : -1;
    return s;
}

/**
 * @brief Try to solve F(x)=0 via multi-dimensional bisection.
 *
 * Works by maintaining n+1 vertices (a simplex in R^n) with distinct
 * sign patterns.  Each bisection step evaluates F at the midpoint of the
 * longest edge and replaces the vertex with matching sign pattern.
 */
SolverStatus BisectionNDSolver::solve(Problem& problem,
                                       Eigen::VectorXd& x,
                                       const SolverOptions& options,
                                       SolverTrace* trace,
                                       std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();
    const int n = problem.size;

    if (trace) {
        if (trace->solverType.empty()) trace->solverType = "BisectionND";
        else trace->solverType += " -> BisectionND";
    }

    // Guard: only for small systems.
    // The maximum block size is controlled by options.bisectionNDMaxBlockSize (default: 8).
    // Blocks exceeding this limit are skipped by the pipeline (InvalidInput is treated as
    // "not applicable" and the next solver in the pipeline is tried instead).
    const int maxDim = (options.bisectionNDMaxBlockSize > 0)
                       ? options.bisectionNDMaxBlockSize
                       : DEFAULT_MAX_BISECTION_DIM;
    if (n <= 0 || x.size() != n) {
        if (detailedError)
            *detailedError = "BisectionND: invalid problem size (n=" + std::to_string(n) + ")";
        return SolverStatus::InvalidInput;
    }
    if (n > maxDim) {
        if (detailedError)
            *detailedError = "BisectionND: block too large (n=" + std::to_string(n)
                           + " > bisectionNDMaxBlockSize=" + std::to_string(maxDim)
                           + "). Increase bisectionNDMaxBlockSize in coolsolve.conf to allow"
                             " larger blocks (warning: cost is exponential in n).";
        return SolverStatus::InvalidInput;
    }

    // Helper: safe evaluate returning success flag
    Eigen::VectorXd F(n);
    Eigen::MatrixXd Jdummy;
    auto safeEval = [&](const Eigen::VectorXd& xv, Eigen::VectorXd& Fout) -> bool {
        try {
            problem.evaluate(xv, Fout, Jdummy, false);
            return Fout.allFinite();
        } catch (...) { return false; }
    };

    // --- Phase 1: Build initial simplex with diverse sign patterns ---
    // We need n+1 vertices with distinct sign patterns.
    // Probe on a structured grid: x0 ± scale_i * e_i for each coordinate,
    // plus random perturbations if needed.

    struct Vertex { Eigen::VectorXd x; Eigen::VectorXd F; std::vector<int> signs; };
    std::vector<Vertex> simplex;
    // Keep a pool of all valid evaluated points for fallback
    std::vector<Vertex> pool;

    Eigen::VectorXd scale = coolsolve::computeScalingFactors(x);

    // Evaluate at the initial guess itself
    {
        Vertex v;
        v.x = x;
        v.F.resize(n);
        if (safeEval(v.x, v.F)) {
            v.signs = signPattern(v.F);
            if (v.F.lpNorm<Eigen::Infinity>() < options.tolerance) {
                // Already converged
                if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                return SolverStatus::Success;
            }
            pool.push_back(v);
        }
    }

    // Probe along each coordinate axis (wider range to cross sign boundaries)
    std::vector<double> probeFactors = {-5.0, -2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0, 5.0};
    for (int i = 0; i < n; ++i) {
        for (double f : probeFactors) {
            Vertex v;
            v.x = x;
            v.x(i) += f * std::max(scale(i), 1.0);
            v.F.resize(n);
            if (safeEval(v.x, v.F)) {
                v.signs = signPattern(v.F);
                if (v.F.lpNorm<Eigen::Infinity>() < options.tolerance) {
                    x = v.x;
                    if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                    return SolverStatus::Success;
                }
                pool.push_back(v);
            }
        }
    }

    // Probe at diagonal (multi-axis) combinations to find all 2^n sign patterns.
    // Axis-aligned probes can miss sign combinations that require simultaneous
    // offsets in multiple directions (e.g. [+,+] when both components must be >0
    // simultaneously). For n ≤ 4 we try corners at two radii {±2, ±5} × scale;
    // for n > 4 a representative subset is used to keep the cost bounded.
    {
        // Try corners at two radii so we can reach sign boundaries farther out
        for (double cornerRad : {2.0, 5.0}) {
            std::vector<std::vector<double>> corners;
            if (n <= 4) {
                // All 2^n corners of {-cornerRad, +cornerRad}^n
                int numCorners = 1 << n;  // 2^n
                for (int mask = 0; mask < numCorners; ++mask) {
                    std::vector<double> offsets(n);
                    for (int k = 0; k < n; ++k)
                        offsets[k] = ((mask >> k) & 1) ? cornerRad : -cornerRad;
                    corners.push_back(offsets);
                }
            } else {
                // For n>4, use a subset: all-positive, all-negative, and axis-flipped variants
                {
                    std::vector<double> o(n, cornerRad);  corners.push_back(o);
                    std::vector<double> m(n, -cornerRad); corners.push_back(m);
                }
                for (int k = 0; k < n; ++k) {
                    std::vector<double> o(n,  cornerRad); o[k] = -cornerRad; corners.push_back(o);
                    std::vector<double> m(n, -cornerRad); m[k] =  cornerRad; corners.push_back(m);
                }
            }
            for (const auto& offsets : corners) {
                Vertex v;
                v.x = x;
                for (int k = 0; k < n; ++k)
                    v.x(k) += offsets[k] * std::max(scale(k), 1.0);
                v.F.resize(n);
                if (safeEval(v.x, v.F)) {
                    v.signs = signPattern(v.F);
                    if (v.F.lpNorm<Eigen::Infinity>() < options.tolerance) {
                        x = v.x;
                        if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                        return SolverStatus::Success;
                    }
                    pool.push_back(v);
                }
            }
        }
    }

    // Also probe at random points in a ball around x
    std::mt19937 rng(42);
    std::normal_distribution<double> randn(0.0, 1.0);
    for (int probe = 0; probe < 20 * n; ++probe) {
        Vertex v;
        v.x = x;
        for (int i = 0; i < n; ++i)
            v.x(i) += randn(rng) * std::max(scale(i), 1.0);
        v.F.resize(n);
        if (safeEval(v.x, v.F)) {
            v.signs = signPattern(v.F);
            if (v.F.lpNorm<Eigen::Infinity>() < options.tolerance) {
                x = v.x;
                if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                return SolverStatus::Success;
            }
            pool.push_back(v);
        }
    }

    // Select n+1 vertices with the most distinct sign patterns
    // Greedy: pick first point, then iteratively pick the point with the
    // most different sign pattern from existing vertices.
    if (pool.empty()) {
        if (detailedError) *detailedError = "BisectionND: no evaluable point found";
        return SolverStatus::EvaluationError;
    }

    auto signDist = [](const std::vector<int>& a, const std::vector<int>& b) {
        int d = 0;
        for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) d++;
        return d;
    };

    // Select up to 2^n vertices with one representative per distinct sign pattern.
    // Having all 2^n sign-pattern combinations in the simplex (one vertex per pattern)
    // is the key requirement for the Vrahatis bisection to converge; with fewer patterns
    // the algorithm may cycle without making progress toward the root.
    // For each pattern we keep the vertex with the smallest ||F||.
    // Greedy: start from the vertex with smallest ||F||, then add one vertex per new
    // sign pattern (choosing the representative with the smallest ||F||).
    const int maxSimplexSize = 1 << n;   // 2^n

    // Pick best initial vertex (smallest |F|)
    simplex.push_back(*std::min_element(pool.begin(), pool.end(),
        [](const Vertex& a, const Vertex& b) {
            return a.F.lpNorm<Eigen::Infinity>() < b.F.lpNorm<Eigen::Infinity>();
        }));

    while (static_cast<int>(simplex.size()) < maxSimplexSize && !pool.empty()) {
        // For each pool point, compute min sign distance to existing simplex vertices.
        // We want to add one representative per ENTIRELY NEW sign pattern.
        int bestIdx = -1;
        int bestMinDist = -1;
        double bestNorm = 1e30;
        for (size_t pi = 0; pi < pool.size(); ++pi) {
            int minDist = n + 1;
            for (const auto& sv : simplex)
                minDist = std::min(minDist, signDist(pool[pi].signs, sv.signs));
            // Only add if this point introduces at least one new sign pattern (minDist>0)
            if (minDist > bestMinDist || (minDist == bestMinDist && pool[pi].F.lpNorm<Eigen::Infinity>() < bestNorm)) {
                bestMinDist = minDist;
                bestIdx = static_cast<int>(pi);
                bestNorm = pool[pi].F.lpNorm<Eigen::Infinity>();
            }
        }
        if (bestIdx < 0 || bestMinDist == 0) break;  // no more distinct sign patterns
        simplex.push_back(pool[bestIdx]);
        pool.erase(pool.begin() + bestIdx);
    }

    if (static_cast<int>(simplex.size()) < 2) {
        // Need at least 2 vertices with different sign patterns to bisect
        if (detailedError) *detailedError = "BisectionND: could not find diverse sign patterns";
        // Return best point found
        x = simplex[0].x;
        return SolverStatus::MaxIterations;
    }

    // --- Phase 2: Iterative bisection ---
    // Bisection converges at a linear rate (one bit of precision per step per dimension),
    // so it needs more iterations than quadratically converging Newton methods.
    // The effective iteration budget is maxIterations * bisectionNDIterFactor.
    const int bisMaxIter = static_cast<int>(
        std::max(1.0, options.maxIterations * std::max(1.0, options.bisectionNDIterFactor)));
    // Repeatedly: find the longest edge, bisect it, replace the vertex
    // whose sign pattern matches the midpoint.
    for (int iter = 0; iter < bisMaxIter; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "BisectionND: timed out";
            break;
        }
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (detailedError) *detailedError = "BisectionND: cancelled";
            break;
        }

        // Find vertex with smallest residual (candidate solution)
        int bestVtx = 0;
        double bestResidual = simplex[0].F.lpNorm<Eigen::Infinity>();
        for (size_t vi = 1; vi < simplex.size(); ++vi) {
            double r = simplex[vi].F.lpNorm<Eigen::Infinity>();
            if (r < bestResidual) { bestResidual = r; bestVtx = static_cast<int>(vi); }
        }

        if (trace) {
            SolverTrace::Iteration ti;
            ti.iter = iter;
            ti.residualNorm = bestResidual;
            ti.stepNorm = 0.0;
            ti.lambda = 1.0;
            ti.x = std::vector<double>(simplex[bestVtx].x.data(),
                                        simplex[bestVtx].x.data() + n);
            ti.residuals = std::vector<double>(simplex[bestVtx].F.data(),
                                               simplex[bestVtx].F.data() + n);
            trace->iterations.push_back(ti);
        }

        if (bestResidual < options.tolerance) {
            x = simplex[bestVtx].x;
            if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
            return SolverStatus::Success;
        }

        // Find the best edge to bisect: longest edge whose midpoint is not
        // already in the simplex (i.e., not too close to any existing vertex).
        // This avoids a degenerate fixed-point loop where the midpoint always
        // maps back to an existing vertex and the simplex never shrinks.
        struct EdgeCandidate { int i, j; double len; };
        std::vector<EdgeCandidate> edges;
        for (size_t i = 0; i < simplex.size(); ++i) {
            for (size_t j = i + 1; j < simplex.size(); ++j) {
                double len = (simplex[i].x - simplex[j].x).norm();
                edges.push_back({static_cast<int>(i), static_cast<int>(j), len});
            }
        }
        std::sort(edges.begin(), edges.end(), [](const EdgeCandidate& a, const EdgeCandidate& b) {
            return a.len > b.len;  // longest first
        });
        double maxLen = edges.empty() ? 0.0 : edges[0].len;

        // Check simplex diameter convergence
        if (maxLen < options.stepTolerance) {
            x = simplex[bestVtx].x;
            if (bestResidual < options.lsRelaxedTolerance) {
                if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                return SolverStatus::Success;
            }
            break;  // Simplex too small, can't bisect further
        }

        // Try each edge (longest first) until we find one whose midpoint
        // is genuinely new (not duplicate of an existing vertex).
        bool bisected = false;
        for (const auto& ec : edges) {
            Eigen::VectorXd xmid = 0.5 * (simplex[ec.i].x + simplex[ec.j].x);

            // Skip if midpoint is too close to any existing vertex (duplicate)
            bool duplicate = false;
            for (const auto& sv : simplex) {
                if ((sv.x - xmid).norm() < options.stepTolerance * 1e-3) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate) continue;

            Vertex vmid;
            vmid.x = xmid;
            vmid.F.resize(n);
            if (!safeEval(vmid.x, vmid.F)) {
                // Midpoint eval failed — try perturbed midpoint
                for (int k = 0; k < n; ++k) vmid.x(k) += 0.01 * randn(rng) * scale(k);
                if (!safeEval(vmid.x, vmid.F)) continue;
            }
            vmid.signs = signPattern(vmid.F);

            if (vmid.F.lpNorm<Eigen::Infinity>() < options.tolerance) {
                x = vmid.x;
                if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                return SolverStatus::Success;
            }

            // Replace the vertex with matching sign pattern
            // Prefer replacing the vertex with larger residual
            int replaceIdx = -1;
            for (size_t vi = 0; vi < simplex.size(); ++vi) {
                if (simplex[vi].signs == vmid.signs) {
                    if (replaceIdx < 0 || simplex[vi].F.norm() > simplex[replaceIdx].F.norm())
                        replaceIdx = static_cast<int>(vi);
                }
            }

            if (replaceIdx >= 0) {
                simplex[replaceIdx] = vmid;
            } else {
                // No matching sign pattern — replace the worst vertex
                int worstIdx = 0;
                double worstNorm = 0.0;
                for (size_t vi = 0; vi < simplex.size(); ++vi) {
                    double nr = simplex[vi].F.norm();
                    if (nr > worstNorm) { worstNorm = nr; worstIdx = static_cast<int>(vi); }
                }
                if (vmid.F.norm() < worstNorm)
                    simplex[worstIdx] = vmid;
            }
            bisected = true;
            break;
        }

        // If ALL edges were degenerate (all midpoints duplicate existing vertices),
        // perturb the worst vertex randomly to escape the stagnant simplex.
        if (!bisected) {
            int worstIdx = 0;
            double worstNorm = 0.0;
            for (size_t vi = 0; vi < simplex.size(); ++vi) {
                double nr = simplex[vi].F.norm();
                if (nr > worstNorm) { worstNorm = nr; worstIdx = static_cast<int>(vi); }
            }
            Vertex vperturb;
            vperturb.x = simplex[worstIdx].x;
            for (int k = 0; k < n; ++k)
                vperturb.x(k) += (0.5 + 0.5 * randn(rng)) * std::max(scale(k), 1.0);
            vperturb.F.resize(n);
            if (safeEval(vperturb.x, vperturb.F)) {
                if (vperturb.F.lpNorm<Eigen::Infinity>() < options.tolerance) {
                    x = vperturb.x;
                    if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                    return SolverStatus::Success;
                }
                vperturb.signs = signPattern(vperturb.F);
                simplex[worstIdx] = vperturb;
            }
        }

        if (trace && !trace->iterations.empty())
            trace->iterations.back().stepNorm = maxLen;
    }

    // Return the best vertex found
    int bestVtx = 0;
    double bestResidual = simplex[0].F.lpNorm<Eigen::Infinity>();
    for (size_t vi = 1; vi < simplex.size(); ++vi) {
        double r = simplex[vi].F.lpNorm<Eigen::Infinity>();
        if (r < bestResidual) { bestResidual = r; bestVtx = static_cast<int>(vi); }
    }
    x = simplex[bestVtx].x;

    if (trace) { trace->finalStatus = SolverStatus::MaxIterations; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
    if (detailedError) {
        std::ostringstream ss;
        ss << "BisectionND: max iterations reached (bisMaxIter=" << bisMaxIter
           << " = maxIterations(" << options.maxIterations
           << ") * bisectionNDIterFactor(" << options.bisectionNDIterFactor
           << ")). Best ||F|| = " << bestResidual
           << ". Increase bisectionNDIterFactor to allow more bisection steps.";
        *detailedError = ss.str();
    }
    return SolverStatus::MaxIterations;
}

}  // namespace coolsolve
