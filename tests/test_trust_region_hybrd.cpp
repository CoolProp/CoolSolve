/**
 * @file test_trust_region_hybrd.cpp
 * @brief Tests for the optional hybrd-style Broyden/QR Jacobian reuse in
 *        TrustRegionSolver (trBroydenRecomputeInterval, trBroydenRestartRejects).
 *
 * These are correctness/regression tests for the *wiring* of the Tier 2.1
 * incremental QR infrastructure (see test_hybrd_qr.cpp) into the
 * trust-region dogleg solver. They are not a substitute for the curated
 * examples-based robustness comparison used to evaluate whether this mode
 * helps on real CoolProp-backed models.
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver.h"
#include <cmath>

using namespace coolsolve;
using Catch::Matchers::WithinAbs;

// ----------------------------------------------------------------------------
// Regression guard: defaults are unchanged (feature is opt-in only).
// ----------------------------------------------------------------------------
TEST_CASE("Defaults leave hybrd Broyden mode disabled", "[trustregion][hybrd][config]") {
    SolverOptions options;
    CHECK(options.trBroydenRecomputeInterval == 0);
    CHECK(options.trBroydenRestartRejects == 2);
}

// ----------------------------------------------------------------------------
// Powell badly-scaled, same problem as the existing default-mode test, but
// with hybrd mode enabled. Confirms enabling Broyden reuse does not break
// convergence on a classic badly-scaled test problem.
// ----------------------------------------------------------------------------
TEST_CASE("TrustRegion hybrd mode - Powell badly-scaled converges", "[trustregion][hybrd]") {
    auto eval = [](const Eigen::VectorXd& x,
                   Eigen::VectorXd& F,
                   Eigen::MatrixXd& J,
                   bool computeJac) {
        F(0) = 1e4 * x(0) * x(1) - 1.0;
        F(1) = std::exp(-x(0)) + std::exp(-x(1)) - 1.0001;
        if (computeJac) {
            J(0, 0) = 1e4 * x(1);          J(0, 1) = 1e4 * x(0);
            J(1, 0) = -std::exp(-x(0));    J(1, 1) = -std::exp(-x(1));
        }
    };

    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = eval;

    Eigen::VectorXd x(2);
    x << 0.0, 1.0;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.maxIterations = 200;
    opts.trBroydenRecomputeInterval = 3;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);

    Eigen::VectorXd F(2);
    Eigen::MatrixXd Jd(2, 2);
    eval(x, F, Jd, false);
    CHECK(F.lpNorm<Eigen::Infinity>() < 1e-6);
}

// ----------------------------------------------------------------------------
// Powell badly-scaled again, but with an aggressive interval (K effectively
// never triggers on its own within maxIterations) so that convergence relies
// on the Powell restart criterion (trBroydenRestartRejects) to recover from
// a stale Broyden approximation. Success here is evidence the restart path
// is functioning as a safety net, not just the periodic K-based refresh.
// ----------------------------------------------------------------------------
TEST_CASE("TrustRegion hybrd mode - restart-on-reject recovers convergence", "[trustregion][hybrd]") {
    auto eval = [](const Eigen::VectorXd& x,
                   Eigen::VectorXd& F,
                   Eigen::MatrixXd& J,
                   bool computeJac) {
        F(0) = 1e4 * x(0) * x(1) - 1.0;
        F(1) = std::exp(-x(0)) + std::exp(-x(1)) - 1.0001;
        if (computeJac) {
            J(0, 0) = 1e4 * x(1);          J(0, 1) = 1e4 * x(0);
            J(1, 0) = -std::exp(-x(0));    J(1, 1) = -std::exp(-x(1));
        }
    };

    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = eval;

    Eigen::VectorXd x(2);
    x << 0.0, 1.0;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.maxIterations = 200;
    opts.trBroydenRecomputeInterval = 1000;
    opts.trBroydenRestartRejects = 2;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);

    Eigen::VectorXd F(2);
    Eigen::MatrixXd Jd(2, 2);
    eval(x, F, Jd, false);
    CHECK(F.lpNorm<Eigen::Infinity>() < 1e-6);
}

// ----------------------------------------------------------------------------
// Mild Rosenbrock sanity check with hybrd mode enabled.
// ----------------------------------------------------------------------------
TEST_CASE("TrustRegion hybrd mode - mild Rosenbrock converges", "[trustregion][hybrd]") {
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool computeJac) {
        F(0) = 1.0 - x(0);
        F(1) = 10.0 * (x(1) - x(0) * x(0));
        if (computeJac) {
            J(0, 0) = -1.0;           J(0, 1) = 0.0;
            J(1, 0) = -20.0 * x(0);   J(1, 1) = 10.0;
        }
    };

    TrustRegionSolver solver;
    Eigen::VectorXd x(2);
    x << -1.0, -1.0;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.trBroydenRecomputeInterval = 2;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinAbs(1.0, 1e-6));
    CHECK_THAT(x(1), WithinAbs(1.0, 1e-6));
}

// ----------------------------------------------------------------------------
// Broyden tridiagonal function (Moré, Garbow, Hillstrom 1981, problem #7):
//   f_i(x) = (3 - 2x_i)x_i - x_{i-1} - 2x_{i+1} + 1,   x_0 = x_{n+1} = 0
// Standard starting point x_i = -1. Chosen because it is the namesake test
// problem for Broyden-type updates and is large enough (n=10) to show a
// meaningful reduction in full-Jacobian evaluations when K > 0.
//
// This counts full-Jacobian evaluations (computeJac == true) via the test's
// own wrapper around the evaluate callback, independent of any internal
// solver instrumentation, and checks that enabling Broyden reuse (K=4)
// requires strictly fewer full-Jacobian evaluations than the legacy
// every-iteration default (K=0) while still converging.
// ----------------------------------------------------------------------------
namespace {
void broydenTridiagonal(const Eigen::VectorXd& x, Eigen::VectorXd& F,
                        Eigen::MatrixXd& J, bool computeJac) {
    int n = static_cast<int>(x.size());
    for (int i = 0; i < n; ++i) {
        double xim1 = (i > 0) ? x(i - 1) : 0.0;
        double xip1 = (i < n - 1) ? x(i + 1) : 0.0;
        F(i) = (3.0 - 2.0 * x(i)) * x(i) - xim1 - 2.0 * xip1 + 1.0;
    }
    if (computeJac) {
        J.setZero();
        for (int i = 0; i < n; ++i) {
            J(i, i) = 3.0 - 4.0 * x(i);
            if (i > 0) J(i, i - 1) = -1.0;
            if (i < n - 1) J(i, i + 1) = -2.0;
        }
    }
}
}  // namespace

TEST_CASE("TrustRegion hybrd mode - Broyden tridiagonal reduces full-Jacobian evaluations", "[trustregion][hybrd]") {
    const int n = 10;

    auto runWithInterval = [&](int K, int& fullJCount) -> SolverStatus {
        fullJCount = 0;
        NonLinearSolver::Problem problem;
        problem.size = n;
        problem.evaluate = [&fullJCount](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                                          Eigen::MatrixXd& J, bool computeJac) {
            if (computeJac) ++fullJCount;
            broydenTridiagonal(x, F, J, computeJac);
        };

        TrustRegionSolver solver;
        Eigen::VectorXd x = Eigen::VectorXd::Constant(n, -1.0);
        SolverOptions opts;
        opts.tolerance = 1e-10;
        opts.maxIterations = 300;
        opts.trBroydenRecomputeInterval = K;
        return solver.solve(problem, x, opts);
    };

    int fullJCount_K0 = 0, fullJCount_K4 = 0;
    auto status0 = runWithInterval(0, fullJCount_K0);
    auto status4 = runWithInterval(4, fullJCount_K4);

    REQUIRE(status0 == SolverStatus::Success);
    REQUIRE(status4 == SolverStatus::Success);
    INFO("Full-Jacobian evaluations: K=0 -> " << fullJCount_K0 << ", K=4 -> " << fullJCount_K4);
    CHECK(fullJCount_K4 < fullJCount_K0);
}
