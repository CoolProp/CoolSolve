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

/// Maximum dimension for multidimensional bisection (exponential cost above this).
static constexpr int MAX_BISECTION_DIM = 8;

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

    // Guard: only for small systems
    if (n <= 0 || n > MAX_BISECTION_DIM || x.size() != n) {
        if (detailedError)
            *detailedError = "BisectionND: system too large (n=" + std::to_string(n)
                           + ", max=" + std::to_string(MAX_BISECTION_DIM) + ")";
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

    // Probe along each coordinate axis
    std::vector<double> probeFactors = {-2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0};
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

    // Pick best initial vertex (smallest |F|)
    simplex.push_back(*std::min_element(pool.begin(), pool.end(),
        [](const Vertex& a, const Vertex& b) {
            return a.F.lpNorm<Eigen::Infinity>() < b.F.lpNorm<Eigen::Infinity>();
        }));

    while (static_cast<int>(simplex.size()) < n + 1 && !pool.empty()) {
        // For each pool point, compute min sign distance to existing simplex vertices
        int bestIdx = -1;
        int bestMinDist = -1;
        double bestNorm = 1e30;
        for (size_t pi = 0; pi < pool.size(); ++pi) {
            int minDist = n + 1;
            for (const auto& sv : simplex)
                minDist = std::min(minDist, signDist(pool[pi].signs, sv.signs));
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
    // Repeatedly: find the longest edge, bisect it, replace the vertex
    // whose sign pattern matches the midpoint.
    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "BisectionND: timed out";
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

        // Find the longest edge in the simplex
        int edge_i = 0, edge_j = 1;
        double maxLen = 0.0;
        for (size_t i = 0; i < simplex.size(); ++i) {
            for (size_t j = i + 1; j < simplex.size(); ++j) {
                double len = (simplex[i].x - simplex[j].x).norm();
                if (len > maxLen) { maxLen = len; edge_i = static_cast<int>(i); edge_j = static_cast<int>(j); }
            }
        }

        // Check simplex diameter convergence
        if (maxLen < options.stepTolerance) {
            x = simplex[bestVtx].x;
            if (bestResidual < options.lsRelaxedTolerance) {
                if (trace) { trace->finalStatus = SolverStatus::Success; trace->totalTime = std::chrono::high_resolution_clock::now() - startTime; }
                return SolverStatus::Success;
            }
            break;  // Simplex too small, can't bisect further
        }

        // Bisect the longest edge
        Eigen::VectorXd xmid = 0.5 * (simplex[edge_i].x + simplex[edge_j].x);
        Vertex vmid;
        vmid.x = xmid;
        vmid.F.resize(n);
        if (!safeEval(vmid.x, vmid.F)) {
            // Midpoint eval failed — try perturbed midpoint
            vmid.x = 0.5 * (simplex[edge_i].x + simplex[edge_j].x);
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
        ss << "BisectionND: max iterations reached. Best ||F|| = " << bestResidual;
        *detailedError = ss.str();
    }
    return SolverStatus::MaxIterations;
}

}  // namespace coolsolve
