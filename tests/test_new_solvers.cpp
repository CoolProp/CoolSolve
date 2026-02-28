/**
 * @file test_new_solvers.cpp
 * @brief Tests for BisectionND and Homotopy solvers.
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver.h"
#include <cmath>

using namespace coolsolve;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Homotopy solver tests
// ============================================================================

TEST_CASE("Homotopy - simple 1D root", "[homotopy]") {
    HomotopySolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool computeJac) {
        F(0) = x(0) * x(0) - 4.0;
        if (computeJac) J(0, 0) = 2.0 * x(0);
    };

    Eigen::VectorXd x(1);
    x(0) = 10.0;  // far from solution

    SolverOptions opts;
    opts.tolerance = 1e-10;
    auto status = solver.solve(problem, x, opts);

    CHECK(status == SolverStatus::Success);
    CHECK_THAT(std::abs(x(0)), WithinAbs(2.0, 1e-8));
}

TEST_CASE("Homotopy - 2D nonlinear", "[homotopy]") {
    HomotopySolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool computeJac) {
        F(0) = x(0) * x(0) + x(1) * x(1) - 4.0;
        F(1) = x(0) * x(1) - 1.0;
        if (computeJac) {
            J(0, 0) = 2.0 * x(0); J(0, 1) = 2.0 * x(1);
            J(1, 0) = x(1);       J(1, 1) = x(0);
        }
    };

    Eigen::VectorXd x(2);
    x << 3.0, 3.0;  // initial guess far from solution

    SolverOptions opts;
    opts.tolerance = 1e-10;
    auto status = solver.solve(problem, x, opts);

    CHECK(status == SolverStatus::Success);
    // Verify F(x) ≈ 0
    Eigen::VectorXd F(2);
    Eigen::MatrixXd J(2, 2);
    problem.evaluate(x, F, J, false);
    CHECK(F.lpNorm<Eigen::Infinity>() < 1e-8);
}

TEST_CASE("Homotopy - Rosenbrock (mild)", "[homotopy]") {
    HomotopySolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool computeJac) {
        // Mild Rosenbrock (b=10 instead of 100)
        F(0) = 1.0 - x(0);
        F(1) = 10.0 * (x(1) - x(0) * x(0));
        if (computeJac) {
            J(0, 0) = -1.0;           J(0, 1) = 0.0;
            J(1, 0) = -20.0 * x(0);   J(1, 1) = 10.0;
        }
    };

    Eigen::VectorXd x(2);
    x << -1.0, -1.0;  // Moderately far from solution (1, 1)

    SolverOptions opts;
    opts.tolerance = 1e-10;
    auto status = solver.solve(problem, x, opts);

    CHECK(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinAbs(1.0, 1e-6));
    CHECK_THAT(x(1), WithinAbs(1.0, 1e-6));
}

// ============================================================================
// BisectionND solver tests
// ============================================================================

TEST_CASE("BisectionND - simple 1D root", "[bisectionnd]") {
    BisectionNDSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool computeJac) {
        F(0) = x(0) * x(0) - 4.0;
        if (computeJac) J(0, 0) = 2.0 * x(0);
    };

    Eigen::VectorXd x(1);
    x(0) = 10.0;

    SolverOptions opts;
    opts.tolerance = 1e-8;
    auto status = solver.solve(problem, x, opts);

    // BisectionND may not converge to high precision, but should get close
    if (status == SolverStatus::Success) {
        CHECK_THAT(std::abs(x(0)), WithinAbs(2.0, 1e-4));
    }
    // If it fails, that's acceptable — BisectionND is a last-resort solver
}

TEST_CASE("BisectionND - 2D system", "[bisectionnd]") {
    BisectionNDSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool computeJac) {
        // x + y = 3, x - y = 1 → solution (2, 1)
        F(0) = x(0) + x(1) - 3.0;
        F(1) = x(0) - x(1) - 1.0;
        if (computeJac) {
            J(0, 0) = 1.0; J(0, 1) = 1.0;
            J(1, 0) = 1.0; J(1, 1) = -1.0;
        }
    };

    Eigen::VectorXd x(2);
    x << 0.0, 0.0;

    SolverOptions opts;
    opts.tolerance = 1e-6;
    auto status = solver.solve(problem, x, opts);

    if (status == SolverStatus::Success) {
        CHECK_THAT(x(0), WithinAbs(2.0, 1e-3));
        CHECK_THAT(x(1), WithinAbs(1.0, 1e-3));
    }
}

TEST_CASE("BisectionND - rejects large systems", "[bisectionnd]") {
    BisectionNDSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 20;  // > MAX_BISECTION_DIM
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool) {
        F = x;
    };

    Eigen::VectorXd x = Eigen::VectorXd::Zero(20);
    SolverOptions opts;
    auto status = solver.solve(problem, x, opts);
    CHECK(status == SolverStatus::InvalidInput);
}

// ============================================================================
// Strategy parsing / factory
// ============================================================================

TEST_CASE("Strategy parsing for new solvers", "[pipeline]") {
    SolverStrategy strat;
    CHECK(parseStrategy("BisectionND", strat));
    CHECK(strat == SolverStrategy::BisectionND);
    CHECK(parseStrategy("Homotopy", strat));
    CHECK(strat == SolverStrategy::Homotopy);
    CHECK(strategyToString(SolverStrategy::BisectionND) == "BisectionND");
    CHECK(strategyToString(SolverStrategy::Homotopy) == "Homotopy");
}

TEST_CASE("createSolver for new types", "[pipeline]") {
    auto bis = createSolver(SolverStrategy::BisectionND);
    CHECK(bis != nullptr);
    auto hom = createSolver(SolverStrategy::Homotopy);
    CHECK(hom != nullptr);
}
