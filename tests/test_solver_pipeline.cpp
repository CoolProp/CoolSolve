/**
 * @file test_solver_pipeline.cpp
 * @brief Tests for the solver pipeline architecture, Levenberg-Marquardt solver,
 *        configurable fallback chains, and parallel execution mode.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

using namespace coolsolve;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Reusable test problems
// ============================================================================

// Simple 1D: f(x) = x^2 - 2 = 0, solution x = sqrt(2)
static void sqrt2Evaluate(const Eigen::VectorXd& x,
                           Eigen::VectorXd& F,
                           Eigen::MatrixXd& J,
                           bool computeJacobian) {
    F(0) = x(0) * x(0) - 2.0;
    if (computeJacobian) {
        J(0, 0) = 2.0 * x(0);
    }
}

// 2D nonlinear: x^2 + y^2 = 4, x*y = 1
static void circle2DEvaluate(const Eigen::VectorXd& x,
                              Eigen::VectorXd& F,
                              Eigen::MatrixXd& J,
                              bool computeJacobian) {
    F(0) = x(0) * x(0) + x(1) * x(1) - 4.0;
    F(1) = x(0) * x(1) - 1.0;
    if (computeJacobian) {
        J(0, 0) = 2.0 * x(0); J(0, 1) = 2.0 * x(1);
        J(1, 0) = x(1);       J(1, 1) = x(0);
    }
}

// Rosenbrock: a - x = 0, b*(y - x^2) = 0  (a=1, b=100)
static void rosenbrockEvaluate(const Eigen::VectorXd& x,
                                Eigen::VectorXd& F,
                                Eigen::MatrixXd& J,
                                bool computeJacobian) {
    const double a = 1.0, b = 100.0;
    F(0) = a - x(0);
    F(1) = b * (x(1) - x(0) * x(0));
    if (computeJacobian) {
        J(0, 0) = -1.0;            J(0, 1) = 0.0;
        J(1, 0) = -2.0 * b * x(0); J(1, 1) = b;
    }
}

// ============================================================================
// SolverStrategy / pipeline utility tests
// ============================================================================

TEST_CASE("strategyToString round-trips", "[pipeline][utility]") {
    CHECK(strategyToString(SolverStrategy::Newton) == "Newton");
    CHECK(strategyToString(SolverStrategy::TrustRegion) == "TrustRegion");
    CHECK(strategyToString(SolverStrategy::LevenbergMarquardt) == "LevenbergMarquardt");
    CHECK(strategyToString(SolverStrategy::Partitioned) == "Partitioned");
}

TEST_CASE("parseStrategy accepts various names", "[pipeline][utility]") {
    SolverStrategy s;
    REQUIRE(parseStrategy("Newton", s));
    CHECK(s == SolverStrategy::Newton);

    REQUIRE(parseStrategy("trustregion", s));
    CHECK(s == SolverStrategy::TrustRegion);

    REQUIRE(parseStrategy("trust_region", s));
    CHECK(s == SolverStrategy::TrustRegion);

    REQUIRE(parseStrategy("TR", s));
    CHECK(s == SolverStrategy::TrustRegion);

    REQUIRE(parseStrategy("LM", s));
    CHECK(s == SolverStrategy::LevenbergMarquardt);

    REQUIRE(parseStrategy("levenberg_marquardt", s));
    CHECK(s == SolverStrategy::LevenbergMarquardt);

    REQUIRE(parseStrategy("partitioned", s));
    CHECK(s == SolverStrategy::Partitioned);

    REQUIRE_FALSE(parseStrategy("bogus", s));
}

TEST_CASE("pipelineModeToString", "[pipeline][utility]") {
    CHECK(pipelineModeToString(SolverPipelineMode::Sequential) == "Sequential");
    CHECK(pipelineModeToString(SolverPipelineMode::Parallel) == "Parallel");
}

TEST_CASE("createSolver returns correct types", "[pipeline][factory]") {
    auto newton = createSolver(SolverStrategy::Newton);
    REQUIRE(newton != nullptr);

    auto tr = createSolver(SolverStrategy::TrustRegion);
    REQUIRE(tr != nullptr);

    auto lm = createSolver(SolverStrategy::LevenbergMarquardt);
    REQUIRE(lm != nullptr);

    // Partitioned is handled specially; factory returns nullptr
    auto part = createSolver(SolverStrategy::Partitioned);
    CHECK(part == nullptr);
}

// ============================================================================
// Levenberg-Marquardt solver tests
// ============================================================================

TEST_CASE("LM solver - sqrt(2)", "[lm][1d]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = sqrt2Evaluate;

    Eigen::VectorXd x(1);
    x(0) = 1.0;

    SolverOptions options;
    options.tolerance = 1e-12;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(std::sqrt(2.0), 1e-10));
}

TEST_CASE("LM solver - 2D circle", "[lm][2d]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = circle2DEvaluate;

    Eigen::VectorXd x(2);
    x << 1.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0) * x(0) + x(1) * x(1), WithinRel(4.0, 1e-8));
    CHECK_THAT(x(0) * x(1), WithinRel(1.0, 1e-8));
}

TEST_CASE("LM solver - Rosenbrock", "[lm][rosenbrock]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << 0.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-8));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-8));
}

TEST_CASE("LM solver - Rosenbrock from poor guess", "[lm][rosenbrock][hard]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = rosenbrockEvaluate;

    Eigen::VectorXd x(2);
    x << -5.0, 10.0;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.maxIterations = 500;

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-6));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-6));
}

TEST_CASE("LM solver - with tracing", "[lm][trace]") {
    LevenbergMarquardtSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = sqrt2Evaluate;

    Eigen::VectorXd x(1);
    x(0) = 1.0;

    SolverOptions options;
    options.tolerance = 1e-12;

    SolverTrace trace;
    SolverStatus status = solver.solve(problem, x, options, &trace);

    REQUIRE(status == SolverStatus::Success);
    REQUIRE(trace.iterations.size() >= 1);
    CHECK(trace.solverType.find("LevenbergMarquardt") != std::string::npos);
}

// ============================================================================
// Pipeline configuration tests (via config file)
// ============================================================================

TEST_CASE("Config: solverPipeline parsing", "[config][pipeline]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_pipeline_test.conf";
    {
        std::ofstream f(configPath);
        REQUIRE(f.is_open());
        f << "solverPipeline = Newton, LM, TrustRegion\n";
        f << "pipelineMode = parallel\n";
    }

    SolverOptions options;
    bool loaded = loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);

    REQUIRE(loaded);
    REQUIRE(options.solverPipeline.size() == 3);
    CHECK(options.solverPipeline[0] == SolverStrategy::Newton);
    CHECK(options.solverPipeline[1] == SolverStrategy::LevenbergMarquardt);
    CHECK(options.solverPipeline[2] == SolverStrategy::TrustRegion);
    CHECK(options.pipelineMode == SolverPipelineMode::Parallel);
}

TEST_CASE("Config: single solver pipeline", "[config][pipeline]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_single_test.conf";
    {
        std::ofstream f(configPath);
        REQUIRE(f.is_open());
        f << "solverPipeline = LM\n";
        f << "pipelineMode = sequential\n";
    }

    SolverOptions options;
    bool loaded = loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);

    REQUIRE(loaded);
    REQUIRE(options.solverPipeline.size() == 1);
    CHECK(options.solverPipeline[0] == SolverStrategy::LevenbergMarquardt);
    CHECK(options.pipelineMode == SolverPipelineMode::Sequential);
}

TEST_CASE("Config: LM options parsing", "[config][lm]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_lm_test.conf";
    {
        std::ofstream f(configPath);
        REQUIRE(f.is_open());
        f << "lmInitialLambda = 0.01\n";
        f << "lmLambdaIncrease = 5.0\n";
        f << "lmLambdaDecrease = 0.2\n";
    }

    SolverOptions options;
    bool loaded = loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);

    REQUIRE(loaded);
    CHECK(options.lmInitialLambda == 0.01);
    CHECK(options.lmLambdaIncrease == 5.0);
    CHECK(options.lmLambdaDecrease == 0.2);
}

// ============================================================================
// Pipeline execution mode tests (using individual solvers directly)
// ============================================================================

TEST_CASE("Newton-only pipeline solves sqrt(2)", "[pipeline][newton-only]") {
    // Verify that a pipeline with only Newton works
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = sqrt2Evaluate;

    Eigen::VectorXd x(1);
    x(0) = 1.0;

    SolverOptions options;
    options.tolerance = 1e-12;
    options.solverPipeline = {SolverStrategy::Newton};

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(std::sqrt(2.0), 1e-10));
}

TEST_CASE("TrustRegion-only pipeline solves circle", "[pipeline][tr-only]") {
    TrustRegionSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = circle2DEvaluate;

    Eigen::VectorXd x(2);
    x << 1.5, 0.5;

    SolverOptions options;
    options.tolerance = 1e-10;
    options.solverPipeline = {SolverStrategy::TrustRegion};

    SolverStatus status = solver.solve(problem, x, options);
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0) * x(0) + x(1) * x(1), WithinRel(4.0, 1e-8));
}

// ============================================================================
// Default pipeline backward compatibility
// ============================================================================

TEST_CASE("Default pipeline includes all four solvers", "[pipeline][defaults]") {
    SolverOptions defaults;
    REQUIRE(defaults.solverPipeline.size() == 4);
    CHECK(defaults.solverPipeline[0] == SolverStrategy::Newton);
    CHECK(defaults.solverPipeline[1] == SolverStrategy::TrustRegion);
    CHECK(defaults.solverPipeline[2] == SolverStrategy::LevenbergMarquardt);
    CHECK(defaults.solverPipeline[3] == SolverStrategy::Partitioned);
    CHECK(defaults.pipelineMode == SolverPipelineMode::Sequential);
}
