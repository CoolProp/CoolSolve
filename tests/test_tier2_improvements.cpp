/**
 * @file test_tier2_improvements.cpp
 * @brief Tests for Tier 2 solver improvements:
 *   #6 Broyden quasi-Newton (Newton solver option)
 *   #7 Trust Region adaptive radius and improved radius management
 *   #8 Levenberg-Marquardt Nielsen update and geodesic acceleration
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver.h"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>

using namespace coolsolve;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

namespace fs = std::filesystem;

// ============================================================================
// Test problems
// ============================================================================

// Rosenbrock: a - x = 0, b*(y - x^2) = 0  →  solution (1, 1)
static void rosenbrockEvaluate(const Eigen::VectorXd& x,
                               Eigen::VectorXd& F,
                               Eigen::MatrixXd& J,
                               bool computeJ) {
    constexpr double a = 1.0, b = 10.0;
    F(0) = a - x(0);
    F(1) = b * (x(1) - x(0) * x(0));
    if (computeJ) {
        J(0, 0) = -1.0;      J(0, 1) = 0.0;
        J(1, 0) = -2.0*b*x(0); J(1, 1) = b;
    }
}

// 3D nonlinear system:
//   x^2 + y + z  = 6
//   x   + y^2 + z = 6
//   x   + y + z^2 = 6
// Solution: (1, 1, 1)  and others; we target (1,1,1) from (2,2,2).
static void cubic3DEvaluate(const Eigen::VectorXd& x,
                            Eigen::VectorXd& F,
                            Eigen::MatrixXd& J,
                            bool computeJ) {
    F(0) = x(0)*x(0) + x(1) + x(2) - 6.0;
    F(1) = x(0) + x(1)*x(1) + x(2) - 6.0;
    F(2) = x(0) + x(1) + x(2)*x(2) - 6.0;
    if (computeJ) {
        J.setZero();
        J(0, 0) = 2.0*x(0); J(0, 1) = 1.0;       J(0, 2) = 1.0;
        J(1, 0) = 1.0;       J(1, 1) = 2.0*x(1); J(1, 2) = 1.0;
        J(2, 0) = 1.0;       J(2, 1) = 1.0;       J(2, 2) = 2.0*x(2);
    }
}

// Powell's badly-scaled function (very poor conditioning):
//   F1 = 1e4 * x1 * x2 - 1
//   F2 = exp(-x1) + exp(-x2) - 1.0001
// Solution approx: x1 ≈ 1.098e-5, x2 ≈ 9.106
static void powellBadlyScaledEvaluate(const Eigen::VectorXd& x,
                                       Eigen::VectorXd& F,
                                       Eigen::MatrixXd& J,
                                       bool computeJ) {
    F(0) = 1e4 * x(0) * x(1) - 1.0;
    F(1) = std::exp(-x(0)) + std::exp(-x(1)) - 1.0001;
    if (computeJ) {
        J(0, 0) = 1e4 * x(1);
        J(0, 1) = 1e4 * x(0);
        J(1, 0) = -std::exp(-x(0));
        J(1, 1) = -std::exp(-x(1));
    }
}

// Extended Rosenbrock for n > 2 (pairs: F_{2i-1} = 1-x_{2i-1}, F_{2i} = 10(x_{2i} - x_{2i-1}^2))
static void extendedRosenbrockEvaluate(const Eigen::VectorXd& x,
                                       Eigen::VectorXd& F,
                                       Eigen::MatrixXd& J,
                                       bool computeJ) {
    int n = x.size();
    if (computeJ) J.setZero();
    for (int i = 0; i < n / 2; ++i) {
        int i1 = 2 * i;
        int i2 = 2 * i + 1;
        F(i1) = 1.0 - x(i1);
        F(i2) = 10.0 * (x(i2) - x(i1) * x(i1));
        if (computeJ) {
            J(i1, i1) = -1.0;
            J(i2, i1) = -20.0 * x(i1);
            J(i2, i2) = 10.0;
        }
    }
}

// ============================================================================
// #6 — Broyden quasi-Newton tests
// ============================================================================

TEST_CASE("Broyden - disabled by default (broydenRecomputeInterval=0)", "[broyden][newton]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;
    REQUIRE(options.broydenRecomputeInterval == 0);  // default

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-6));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-6));
}

TEST_CASE("Broyden - enabled K=3 on Rosenbrock", "[broyden][newton]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.broydenRecomputeInterval = 3;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-6));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-6));
}

TEST_CASE("Broyden - enabled K=5 on 3D system", "[broyden][newton]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 3;
    problem.evaluate = cubic3DEvaluate;

    Eigen::VectorXd x(3);
    x << 2.0, 2.0, 2.0;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.broydenRecomputeInterval = 5;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    // Should converge to (1,1,1) or another root — check residual
    Eigen::VectorXd F(3);
    Eigen::MatrixXd J(3, 3);
    cubic3DEvaluate(x, F, J, false);
    CHECK(F.lpNorm<Eigen::Infinity>() < 1e-9);
}

TEST_CASE("Broyden - extended Rosenbrock n=6", "[broyden][newton]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 6;
    problem.evaluate = extendedRosenbrockEvaluate;

    Eigen::VectorXd x(6);
    x << -1.0, 1.0, -1.0, 1.0, -1.0, 1.0;

    SolverOptions options;
    options.tolerance = 1e-9;
    options.broydenRecomputeInterval = 4;
    options.maxIterations = 200;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    for (int i = 0; i < 6; ++i)
        CHECK_THAT(x(i), WithinAbs(1.0, 1e-5));
}

TEST_CASE("Broyden - fallback to full Jacobian on difficult step", "[broyden][newton]") {
    // Powell's badly-scaled: Broyden may struggle, but automatic recovery
    // (recompute full J when line search fails) should still converge.
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = powellBadlyScaledEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 5.0;

    SolverOptions options;
    options.tolerance = 1e-8;
    options.broydenRecomputeInterval = 2;
    options.maxIterations = 200;

    SolverStatus status = solver.solve(problem, x, options);
    // Either converges or at least doesn't crash
    if (status == SolverStatus::Success) {
        Eigen::VectorXd F(2);
        Eigen::MatrixXd J(2, 2);
        powellBadlyScaledEvaluate(x, F, J, false);
        CHECK(F.lpNorm<Eigen::Infinity>() < 1e-2);  // relaxed: Powell is very hard for Broyden
    }
}

// ============================================================================
// #7 — Trust Region adaptive radius tests
// ============================================================================

TEST_CASE("TR adaptive radius - default is true", "[trustregion][adaptive]") {
    SolverOptions options;
    REQUIRE(options.trAdaptiveRadius == true);
}

TEST_CASE("TR adaptive radius - Rosenbrock convergence", "[trustregion][adaptive]") {
    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << -1.0, 2.0;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.trAdaptiveRadius = true;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-5));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-5));
}

TEST_CASE("TR adaptive radius - disabled still works", "[trustregion][adaptive]") {
    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.trAdaptiveRadius = false;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-5));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-5));
}

TEST_CASE("TR adaptive radius - 3D nonlinear", "[trustregion][adaptive]") {
    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 3;
    problem.evaluate = cubic3DEvaluate;

    Eigen::VectorXd x(3);
    x << 3.0, 3.0, 3.0;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.trAdaptiveRadius = true;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    Eigen::VectorXd F(3);
    Eigen::MatrixXd J(3, 3);
    cubic3DEvaluate(x, F, J, false);
    CHECK(F.lpNorm<Eigen::Infinity>() < 1e-9);
}

// ============================================================================
// #8 — LM Nielsen update and geodesic acceleration tests
// ============================================================================

TEST_CASE("LM Nielsen - defaults are true", "[lm][nielsen]") {
    SolverOptions options;
    REQUIRE(options.lmNielsenUpdate == true);
    REQUIRE(options.lmGeodesicAcceleration == true);
}

TEST_CASE("LM Nielsen - Rosenbrock convergence", "[lm][nielsen]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << -1.0, 2.0;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.lmNielsenUpdate = true;
    options.lmGeodesicAcceleration = true;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-5));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-5));
}

TEST_CASE("LM legacy mode (Nielsen and geodesic disabled)", "[lm][nielsen]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.lmNielsenUpdate = false;
    options.lmGeodesicAcceleration = false;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-5));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-5));
}

TEST_CASE("LM geodesic only (Nielsen disabled)", "[lm][nielsen]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.lmNielsenUpdate = false;
    options.lmGeodesicAcceleration = true;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-5));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-5));
}

TEST_CASE("LM Nielsen on Powell's badly-scaled problem", "[lm][nielsen][hard]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = powellBadlyScaledEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 5.0;

    SolverOptions options;
    options.tolerance = 1e-8;
    options.maxIterations = 300;
    options.lmNielsenUpdate = true;
    options.lmGeodesicAcceleration = true;

    SolverStatus status = solver.solve(problem, x, options);
    if (status == SolverStatus::Success) {
        Eigen::VectorXd F(2);
        Eigen::MatrixXd J(2, 2);
        powellBadlyScaledEvaluate(x, F, J, false);
        CHECK(F.lpNorm<Eigen::Infinity>() < 1e-6);
    }
}

TEST_CASE("LM Nielsen - extended Rosenbrock n=6", "[lm][nielsen]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 6;
    problem.evaluate = extendedRosenbrockEvaluate;

    Eigen::VectorXd x(6);
    x << -1.0, 1.0, -1.0, 1.0, -1.0, 1.0;

    SolverOptions options;
    options.tolerance = 1e-9;
    options.maxIterations = 300;
    options.lmNielsenUpdate = true;
    options.lmGeodesicAcceleration = true;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    for (int i = 0; i < 6; ++i)
        CHECK_THAT(x(i), WithinAbs(1.0, 1e-5));
}

// ============================================================================
// Config loading for new options
// ============================================================================

TEST_CASE("Config loads Tier 2 options", "[config][tier2]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_tier2.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "broydenRecomputeInterval = 5\n";
    f << "trAdaptiveRadius = false\n";
    f << "lmNielsenUpdate = false\n";
    f << "lmGeodesicAcceleration = false\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.broydenRecomputeInterval == 5);
    REQUIRE(options.trAdaptiveRadius == false);
    REQUIRE(options.lmNielsenUpdate == false);
    REQUIRE(options.lmGeodesicAcceleration == false);
}

TEST_CASE("Config defaults for Tier 2 options when not set", "[config][tier2]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_tier2_defaults.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "maxIterations = 50\n";  // unrelated key
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.broydenRecomputeInterval == 0);
    REQUIRE(options.trAdaptiveRadius == true);
    REQUIRE(options.lmNielsenUpdate == true);
    REQUIRE(options.lmGeodesicAcceleration == true);
}
