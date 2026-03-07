/**
 * @file test_nonmonotone.cpp
 * @brief Tests for the non-monotone line search feature (Grippo et al. 1986).
 *
 * Covers:
 *   - NonMonotoneHistory helper class (unit tests)
 *   - Newton solver with non-monotone M > 1 vs monotone M = 1
 *   - TrustRegion solver with non-monotone acceptance
 *   - LM solver with non-monotone acceptance
 *   - Config loading for lsNonMonotoneMemory
 *   - Tracing output includes non-monotone information
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include <cmath>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
using namespace coolsolve;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Test problems
// ============================================================================

// Rosenbrock roots (narrow curved valley — classic non-monotone benchmark)
struct RosenbrockRoots {
    double a = 1.0;
    double b = 100.0;
    void operator()(const Eigen::VectorXd& x,
                    Eigen::VectorXd& F,
                    Eigen::MatrixXd& J,
                    bool computeJacobian) const {
        F(0) = a - x(0);
        F(1) = b * (x(1) - x(0) * x(0));
        if (computeJacobian) {
            J(0, 0) = -1.0;              J(0, 1) = 0.0;
            J(1, 0) = -2.0 * b * x(0);   J(1, 1) = b;
        }
    }
};

// Powell's singular function (4D, near-singular Jacobian at solution)
struct PowellSingular {
    void operator()(const Eigen::VectorXd& x,
                    Eigen::VectorXd& F,
                    Eigen::MatrixXd& J,
                    bool computeJacobian) const {
        F(0) = x(0) + 10.0 * x(1);
        F(1) = std::sqrt(5.0) * (x(2) - x(3));
        F(2) = (x(1) - 2.0 * x(2)) * (x(1) - 2.0 * x(2));
        F(3) = std::sqrt(10.0) * (x(0) - x(3)) * (x(0) - x(3));
        if (computeJacobian) {
            J.setZero();
            J(0, 0) = 1.0;  J(0, 1) = 10.0;
            J(1, 2) = std::sqrt(5.0);  J(1, 3) = -std::sqrt(5.0);
            J(2, 1) = 2.0 * (x(1) - 2.0 * x(2));
            J(2, 2) = -4.0 * (x(1) - 2.0 * x(2));
            J(3, 0) = 2.0 * std::sqrt(10.0) * (x(0) - x(3));
            J(3, 3) = -2.0 * std::sqrt(10.0) * (x(0) - x(3));
        }
    }
};

// Simple 2D transcendental (sanity check)
struct Exp2D {
    void operator()(const Eigen::VectorXd& x,
                    Eigen::VectorXd& F,
                    Eigen::MatrixXd& J,
                    bool computeJacobian) const {
        F(0) = std::exp(x(0)) + x(1) - 2.0;
        F(1) = x(0) + std::exp(x(1)) - 2.0;
        if (computeJacobian) {
            J(0, 0) = std::exp(x(0));  J(0, 1) = 1.0;
            J(1, 0) = 1.0;             J(1, 1) = std::exp(x(1));
        }
    }
};

// ============================================================================
// NonMonotoneHistory unit tests
// ============================================================================

TEST_CASE("NonMonotoneHistory - basic push and maxValue", "[nonmonotone][unit]") {
    NonMonotoneHistory h(5);
    REQUIRE(h.empty());
    REQUIRE(h.maxValue() == 0.0);  // empty → 0

    h.push(3.0);
    CHECK(h.size() == 1);
    CHECK(h.maxValue() == 3.0);

    h.push(1.0);
    h.push(5.0);
    CHECK(h.size() == 3);
    CHECK(h.maxValue() == 5.0);  // max of {3, 1, 5}
}

TEST_CASE("NonMonotoneHistory - memory eviction", "[nonmonotone][unit]") {
    NonMonotoneHistory h(3);  // keep last 3 values
    h.push(10.0);  // {10}
    h.push(2.0);   // {10, 2}
    h.push(5.0);   // {10, 2, 5}
    CHECK(h.maxValue() == 10.0);
    CHECK(h.size() == 3);

    h.push(1.0);   // evict 10 → {2, 5, 1}
    CHECK(h.size() == 3);
    CHECK(h.maxValue() == 5.0);  // 10 was evicted

    h.push(0.5);   // evict 2 → {5, 1, 0.5}
    CHECK(h.maxValue() == 5.0);

    h.push(0.1);   // evict 5 → {1, 0.5, 0.1}
    CHECK(h.maxValue() == 1.0);
}

TEST_CASE("NonMonotoneHistory - memory=1 gives monotone", "[nonmonotone][unit]") {
    NonMonotoneHistory h(1);  // M=1 → monotone
    h.push(5.0);
    CHECK(h.maxValue() == 5.0);

    h.push(3.0);  // evicts 5 immediately
    CHECK(h.maxValue() == 3.0);  // = current value → same as monotone

    h.push(7.0);
    CHECK(h.maxValue() == 7.0);
}

TEST_CASE("NonMonotoneHistory - memory=0 clamped to 1", "[nonmonotone][unit]") {
    NonMonotoneHistory h(0);  // should clamp to 1
    h.push(4.0);
    CHECK(h.size() == 1);
    h.push(2.0);
    CHECK(h.size() == 1);  // old evicted since memory=1
    CHECK(h.maxValue() == 2.0);
}

TEST_CASE("NonMonotoneHistory - negative memory clamped to 1", "[nonmonotone][unit]") {
    NonMonotoneHistory h(-5);  // should clamp to 1
    h.push(4.0);
    h.push(2.0);
    CHECK(h.size() == 1);
    CHECK(h.maxValue() == 2.0);
}

TEST_CASE("NonMonotoneHistory - large memory", "[nonmonotone][unit]") {
    NonMonotoneHistory h(1000);
    for (int i = 0; i < 500; ++i) h.push(static_cast<double>(i));
    CHECK(h.size() == 500);
    CHECK(h.maxValue() == 499.0);

    // Fill beyond 1000
    for (int i = 500; i < 1200; ++i) h.push(static_cast<double>(i));
    CHECK(h.size() == 1000);
    CHECK(h.maxValue() == 1199.0);
}

TEST_CASE("NonMonotoneHistory - boundedRef caps reference", "[nonmonotone][unit]") {
    NonMonotoneHistory h(10);
    // Push a large early value that would normally dominate
    h.push(1000.0);
    h.push(100.0);
    h.push(10.0);
    h.push(1.0);

    // maxValue is the unbounded max
    CHECK(h.maxValue() == 1000.0);

    // boundedRef(currentPhi, maxRatio) caps at currentPhi * maxRatio
    double currentPhi = 1.0;
    CHECK(h.boundedRef(currentPhi, 10.0) == 10.0);    // min(1000, 1*10) = 10
    CHECK(h.boundedRef(currentPhi, 100.0) == 100.0);   // min(1000, 1*100) = 100
    CHECK(h.boundedRef(currentPhi, 2000.0) == 1000.0); // min(1000, 2000) = 1000 (no cap)

    // When current phi equals max, cap is always >= max
    CHECK(h.boundedRef(1000.0, 10.0) == 1000.0); // min(1000, 10000) = 1000

    // Default maxRatio = 10
    CHECK(h.boundedRef(2.0) == 20.0);  // min(1000, 2*10) = 20
}

TEST_CASE("NonMonotoneHistory - boundedRef with M=1", "[nonmonotone][unit]") {
    NonMonotoneHistory h(1);
    h.push(5.0);
    // M=1: maxValue = currentPhi, boundedRef = min(5, 5*10) = 5
    CHECK(h.boundedRef(5.0) == 5.0);
    h.push(3.0);
    CHECK(h.boundedRef(3.0) == 3.0);   // Same as monotone
}

// ============================================================================
// Newton solver – non-monotone vs monotone
// ============================================================================

TEST_CASE("Newton - non-monotone converges on Rosenbrock from poor guess",
          "[nonmonotone][newton]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = RosenbrockRoots();

    Eigen::VectorXd x(2);
    x << -5.0, 10.0;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.maxIterations = 200;
    opts.lsNonMonotoneMemory = 10;  // explicit non-monotone

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-8));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-8));
}

TEST_CASE("Newton - M=1 (monotone) also converges on easy problems",
          "[nonmonotone][newton]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = Exp2D();

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions opts;
    opts.tolerance = 1e-12;
    opts.lsNonMonotoneMemory = 1;  // monotone

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(std::exp(x(0)) + x(1), WithinRel(2.0, 1e-10));
}

TEST_CASE("Newton - non-monotone M=10 vs M=1 on Rosenbrock",
          "[nonmonotone][newton][compare]") {
    // Both should converge, but non-monotone may use fewer iterations
    RosenbrockRoots rosen;

    auto runWith = [&](int memory) -> std::pair<SolverStatus, int> {
        NewtonSolver solver;
        NonLinearSolver::Problem problem;
        problem.size = 2;
        problem.evaluate = rosen;

        Eigen::VectorXd x(2);
        x << -3.0, 5.0;

        SolverOptions opts;
        opts.tolerance = 1e-10;
        opts.maxIterations = 300;
        opts.lsNonMonotoneMemory = memory;

        SolverTrace trace;
        auto status = solver.solve(problem, x, opts, &trace);
        return {status, static_cast<int>(trace.iterations.size())};
    };

    auto [statusM1,  itersM1]  = runWith(1);
    auto [statusM10, itersM10] = runWith(10);

    // Both should converge
    CHECK(statusM1  == SolverStatus::Success);
    CHECK(statusM10 == SolverStatus::Success);

    // Non-monotone version should use no more iterations (typically fewer)
    INFO("Monotone (M=1) iterations: " << itersM1);
    INFO("Non-monotone (M=10) iterations: " << itersM10);
    // We don't require strictly fewer — just that both converge
}

TEST_CASE("Newton - non-monotone on Powell singular",
          "[nonmonotone][newton][powell]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 4;
    problem.evaluate = PowellSingular();

    Eigen::VectorXd x(4);
    x << 3.0, -1.0, 0.0, 1.0;  // standard starting point

    SolverOptions opts;
    opts.tolerance = 1e-8;
    opts.maxIterations = 500;
    opts.lsNonMonotoneMemory = 10;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
    // Solution is the origin
    for (int i = 0; i < 4; ++i)
        CHECK_THAT(x(i), WithinAbs(0.0, 1e-4));
}

// ============================================================================
// TrustRegion solver – non-monotone acceptance
// ============================================================================

TEST_CASE("TrustRegion - non-monotone converges on Rosenbrock",
          "[nonmonotone][trustregion]") {
    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = RosenbrockRoots();

    Eigen::VectorXd x(2);
    x << -3.0, 5.0;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.maxIterations = 500;
    opts.lsNonMonotoneMemory = 10;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-6));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-6));
}

TEST_CASE("TrustRegion - M=1 still works (monotone fallback)",
          "[nonmonotone][trustregion]") {
    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = Exp2D();

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.lsNonMonotoneMemory = 1;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
}

// ============================================================================
// LM solver – non-monotone acceptance
// ============================================================================

TEST_CASE("LM - non-monotone converges on Rosenbrock",
          "[nonmonotone][lm]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = RosenbrockRoots();

    Eigen::VectorXd x(2);
    x << -3.0, 5.0;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.maxIterations = 500;
    opts.lsNonMonotoneMemory = 10;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-6));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-6));
}

TEST_CASE("LM - M=1 still works (monotone fallback)",
          "[nonmonotone][lm]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = Exp2D();

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.lsNonMonotoneMemory = 1;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
}

TEST_CASE("LM - non-monotone on Powell singular",
          "[nonmonotone][lm][powell]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 4;
    problem.evaluate = PowellSingular();

    Eigen::VectorXd x(4);
    x << 3.0, -1.0, 0.0, 1.0;

    SolverOptions opts;
    opts.tolerance = 1e-8;
    opts.maxIterations = 1000;
    opts.lsNonMonotoneMemory = 10;

    auto status = solver.solve(problem, x, opts);
    REQUIRE(status == SolverStatus::Success);
    for (int i = 0; i < 4; ++i)
        CHECK_THAT(x(i), WithinAbs(0.0, 1e-4));
}

// ============================================================================
// Config loading
// ============================================================================

TEST_CASE("Config loads lsNonMonotoneMemory", "[nonmonotone][config]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_nonmono.conf";
    {
        std::ofstream f(configPath);
        REQUIRE(f.is_open());
        f << "lsNonMonotoneMemory = 5\n";
    }
    SolverOptions opts;
    bool loaded = loadSolverOptionsFromFile(configPath.string(), opts);
    fs::remove(configPath);
    REQUIRE(loaded);
    CHECK(opts.lsNonMonotoneMemory == 5);
}

TEST_CASE("Config lsNonMonotoneMemory = 1 gives monotone", "[nonmonotone][config]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_nonmono2.conf";
    {
        std::ofstream f(configPath);
        REQUIRE(f.is_open());
        f << "lsNonMonotoneMemory = 1\n";
    }
    SolverOptions opts;
    bool loaded = loadSolverOptionsFromFile(configPath.string(), opts);
    fs::remove(configPath);
    REQUIRE(loaded);
    CHECK(opts.lsNonMonotoneMemory == 1);
}

// ============================================================================
// Trace output verification
// ============================================================================

TEST_CASE("Newton trace records iterations with non-monotone",
          "[nonmonotone][newton][trace]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = RosenbrockRoots();

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions opts;
    opts.tolerance = 1e-10;
    opts.lsNonMonotoneMemory = 10;

    SolverTrace trace;
    auto status = solver.solve(problem, x, opts, &trace);
    REQUIRE(status == SolverStatus::Success);
    REQUIRE(trace.iterations.size() >= 1);
    REQUIRE(trace.finalStatus == SolverStatus::Success);

    // Residuals should eventually converge to near zero
    CHECK(trace.iterations.back().residualNorm < 1e-9);
}
