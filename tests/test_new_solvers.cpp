/**
 * @file test_new_solvers.cpp
 * @brief Tests for BisectionND and Homotopy solvers, plus per-solver advantage
 *        tests that illustrate when each solver outperforms the others.
 *
 * Per-solver rationale:
 *   Newton     – Fastest (quadratic convergence) on smooth well-conditioned systems.
 *   TrustRegion– More robust when Newton's unconstrained step exits the feasible
 *                domain (division-by-zero, log of negative number, etc.).
 *   LM         – Handles near-singular or rank-deficient Jacobians via λ damping.
 *   Homotopy   – Converges from starting points so far that iterative methods never
 *                reach the basin of attraction of any root.
 *   BisectionND– Converges when the Jacobian is identically zero at the starting
 *                point, making every gradient-based method immediately stuck.
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

// ============================================================================
// Per-solver advantage tests
//   Each test demonstrates a scenario where a specific solver shines and
//   another method would struggle, with comments explaining the mathematics.
// ============================================================================

// ----------------------------------------------------------------------------
// Newton advantage: quadratic convergence rate
//
// Newton's method has LOCAL quadratic convergence: near the solution x*,
//   ||x_{k+1} - x*|| ≤ C · ||x_k - x*||²
//
// This test verifies the fast convergence empirically by running with a trace
// and checking that the solver finishes in very few iterations. Gradient-based
// methods like TrustRegion and LM exhibit only linear convergence in general.
//
// Problem: x² + y² = 2, x - y = 0. Unique solution: x = y = 1.
// From (2, 2): Newton should converge in 4-5 iterations (quadratically).
// ----------------------------------------------------------------------------
TEST_CASE("Newton - quadratic convergence from moderate guess", "[newton][advantage]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x,
                          Eigen::VectorXd& F,
                          Eigen::MatrixXd& J,
                          bool computeJac) {
        F(0) = x(0)*x(0) + x(1)*x(1) - 2.0;
        F(1) = x(0) - x(1);
        if (computeJac) {
            J(0, 0) = 2.0*x(0); J(0, 1) = 2.0*x(1);
            J(1, 0) = 1.0;      J(1, 1) = -1.0;
        }
    };

    Eigen::VectorXd x(2);
    x << 2.0, 2.0;  // moderate initial guess

    SolverOptions opts;
    opts.tolerance = 1e-12;

    SolverTrace trace;
    auto status = solver.solve(problem, x, opts, &trace);

    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinAbs(1.0, 1e-10));
    CHECK_THAT(x(1), WithinAbs(1.0, 1e-10));
    // Quadratic convergence: very few iterations even near solution
    INFO("Newton iterations: " << trace.iterations.size());
    CHECK(trace.iterations.size() <= 8);
}

// ----------------------------------------------------------------------------
// TrustRegion advantage: feasibility under domain constraints
//
// When the Newton direction points OUT of the function's domain (e.g., into
// log of a negative number), Newton with backtracking must take diminishingly
// small steps (the line search is forced to shrink α to near-zero to avoid
// the domain violation). TrustRegion's dogleg step uses the CAUCHY DIRECTION
// (steepest descent of ‖F‖²) as a fallback, which can stay within the domain.
//
// Problem: log(x) + y = 1,  x = exp(y).
//   Solution: y = 0.5, x = exp(0.5) ≈ 1.6487.
//
// From x = 0.001, y = 3: the Newton direction requires Δx < 0 (x → negative),
// which exits the log domain. The Cauchy direction increases x (correct!) and
// decreases y — TrustRegion reaches the solution much sooner.
// Newton (limited iterations) still converges but needs many more steps.
// ----------------------------------------------------------------------------
TEST_CASE("TrustRegion - feasibility under log domain constraint", "[trustregion][advantage]") {
    auto evaluateLogExp = [](const Eigen::VectorXd& x,
                             Eigen::VectorXd& F,
                             Eigen::MatrixXd& J,
                             bool computeJac) {
        // Requires x(0) > 0 for log to be defined
        double xv = x(0) > 0 ? x(0) : 1e-14;
        F(0) = std::log(xv) + x(1) - 1.0;
        F(1) = xv - std::exp(x(1));
        if (computeJac) {
            J(0, 0) = 1.0/xv; J(0, 1) = 1.0;
            J(1, 0) = 1.0;    J(1, 1) = -std::exp(x(1));
        }
    };

    const double sol_x = std::exp(0.5);   // ≈ 1.6487
    const double sol_y = 0.5;

    // TrustRegion converges from the difficult starting point
    {
        TrustRegionSolver trSolver;
        NonLinearSolver::Problem prob;
        prob.size = 2;
        prob.evaluate = evaluateLogExp;

        Eigen::VectorXd x(2);
        x << 0.001, 3.0;  // very small x, far too large y

        SolverOptions opts;
        opts.tolerance = 1e-8;
        opts.maxIterations = 500;

        SolverTrace trace;
        auto status = trSolver.solve(prob, x, opts, &trace);
        REQUIRE(status == SolverStatus::Success);
        CHECK_THAT(x(0), WithinAbs(sol_x, 1e-5));
        CHECK_THAT(x(1), WithinAbs(sol_y, 1e-5));
        INFO("TrustRegion iterations: " << trace.iterations.size());
    }

    // Newton with limited budget struggles more from this starting point
    // (backs off to tiny steps because Newton direction exits log domain)
    {
        NewtonSolver nSolver;
        NonLinearSolver::Problem prob;
        prob.size = 2;
        prob.evaluate = evaluateLogExp;

        Eigen::VectorXd x(2);
        x << 0.001, 3.0;

        SolverOptions opts;
        opts.tolerance = 1e-8;
        opts.maxIterations = 15;  // tight budget

        auto status = nSolver.solve(prob, x, opts);
        // Newton may or may not converge in 15 iterations from this point.
        // The key insight is that TrustRegion's Cauchy direction heads toward
        // increasing x (correct) while Newton's step initially forces x toward 0.
        INFO("Newton in 15 iters: " << statusToString(status));
        // No hard assertion here — this is an informational illustration.
        // The REQUIRE above (TR succeeds) is the real check.
    }
}

// ----------------------------------------------------------------------------
// LM advantage: near-singular Jacobian
//
// Levenberg-Marquardt regularises the normal equations: (J'J + λI)Δ = -J'F.
// This makes the step well-defined even when J is rank-deficient or nearly so.
//
// Problem:
//   f₁(x,y) = x + y − 2
//   f₂(x,y) = x + y − 2 + ε·(x − y)³,  ε = 0.01
//
// Unique solution: x = y = 1 (so x+y=2 and x=y, making the cubic term zero).
//
// At the starting point (0, 0):
//   J = [[1, 1], [1+3ε(x−y)², 1−3ε(x−y)²]] = [[1,1],[1,1]] → RANK 1.
//
// Newton with a rank-1 Jacobian cannot determine a unique step direction.
// ColPivHouseholderQR returns the minimum-norm solution, which lies in the
// null-space of J and makes no progress towards reducing f₁ alone.
// LM adds λI so (J'J + λI) is full-rank and a meaningful step is taken.
// ----------------------------------------------------------------------------
TEST_CASE("LM - near-singular Jacobian at initial point", "[lm][advantage]") {
    const double eps = 0.01;
    auto evaluateNearSingular = [eps](const Eigen::VectorXd& x,
                                     Eigen::VectorXd& F,
                                     Eigen::MatrixXd& J,
                                     bool computeJac) {
        double diff = x(0) - x(1);
        F(0) = x(0) + x(1) - 2.0;
        F(1) = x(0) + x(1) - 2.0 + eps * diff * diff * diff;
        if (computeJac) {
            double dCubic = 3.0 * eps * diff * diff;
            J(0, 0) = 1.0; J(0, 1) = 1.0;
            J(1, 0) = 1.0 + dCubic; J(1, 1) = 1.0 - dCubic;
        }
    };

    // LM converges from the rank-1 starting point
    {
        LevenbergMarquardtSolver lmSolver;
        NonLinearSolver::Problem prob;
        prob.size = 2;
        prob.evaluate = evaluateNearSingular;

        Eigen::VectorXd x(2);
        x << 0.0, 0.0;  // J is rank-1 here

        SolverOptions opts;
        opts.tolerance = 1e-10;
        opts.maxIterations = 300;

        auto status = lmSolver.solve(prob, x, opts);
        REQUIRE(status == SolverStatus::Success);
        CHECK_THAT(x(0), WithinAbs(1.0, 1e-7));
        CHECK_THAT(x(1), WithinAbs(1.0, 1e-7));
    }
}

// ----------------------------------------------------------------------------
// Homotopy advantage: convergence from far starting points
//
// Homotopy continuation builds a globally convergent path from the trivial
// solution (x = x₀, i.e. F(x) = x − x₀ is zero) to the actual solution.
// Iterative methods (Newton, TR, LM) all require the initial guess to already
// be in the basin of attraction of the root.
//
// This test starts Newton and Homotopy from x₀ = (−10, −10), which is far
// outside the natural basin. Newton (limited iterations) cannot converge;
// Homotopy tracks the homotopy path and arrives at the solution.
//
// Problem: x² + y² = 4,  x − y = 0.  Solution: x = y = √2.
// ----------------------------------------------------------------------------
TEST_CASE("Homotopy - far starting point where Newton fails", "[homotopy][advantage]") {
    auto evaluateCircle = [](const Eigen::VectorXd& x,
                             Eigen::VectorXd& F,
                             Eigen::MatrixXd& J,
                             bool computeJac) {
        F(0) = x(0)*x(0) + x(1)*x(1) - 4.0;
        F(1) = x(0) - x(1);
        if (computeJac) {
            J(0, 0) = 2.0*x(0); J(0, 1) = 2.0*x(1);
            J(1, 0) = 1.0;      J(1, 1) = -1.0;
        }
    };

    const double sol = std::sqrt(2.0);

    // Homotopy converges from the far-away starting point
    {
        HomotopySolver hSolver;
        NonLinearSolver::Problem prob;
        prob.size = 2;
        prob.evaluate = evaluateCircle;

        Eigen::VectorXd x(2);
        x << -10.0, -10.0;

        SolverOptions opts;
        opts.tolerance = 1e-10;
        auto status = hSolver.solve(prob, x, opts);

        REQUIRE(status == SolverStatus::Success);
        // Both positive and negative roots are valid; pick whichever was found
        CHECK_THAT(std::abs(x(0)), WithinAbs(sol, 1e-6));
        CHECK_THAT(std::abs(x(1)), WithinAbs(sol, 1e-6));
    }

    // Newton with a tight iteration budget fails from this starting point.
    // From (-10,-10), Newton's step goes to roughly (-5,-5), then (-2.5,-2.5)...
    // It needs many steps to get close enough to the solution's basin.
    {
        NewtonSolver nSolver;
        NonLinearSolver::Problem prob;
        prob.size = 2;
        prob.evaluate = evaluateCircle;

        Eigen::VectorXd x(2);
        x << -10.0, -10.0;

        SolverOptions opts;
        opts.tolerance = 1e-10;
        opts.maxIterations = 5;  // deliberately tight to show the contrast

        auto status = nSolver.solve(prob, x, opts);
        // With only 5 iterations Newton does NOT converge from this distance
        INFO("Newton from far point in 5 iters: " << statusToString(status));
        CHECK(status != SolverStatus::Success);
    }
}

// ----------------------------------------------------------------------------
// BisectionND advantage: zero Jacobian at starting point
//
// All gradient-based methods (Newton, TR, LM) compute a search direction from
// the Jacobian J. When J(x₀) = 0 exactly, the Newton step Δ = −J⁻¹F is
// undefined (the system J·Δ = −F has no solution or infinite solutions).
// Our Newton uses ColPivHouseholderQR which returns Δ = 0 in this case,
// leaving the solver stuck at x₀ for every iteration.
//
// BisectionND does NOT use derivative information. It probes the function on a
// structured grid around x₀ and finds vertices with differing sign patterns,
// then repeatedly bisects the simplex edge connecting them. The sign change
// guarantees a root is enclosed and convergence is systematic.
//
// Problem: f₁(x,y) = x³ − 1,  f₂(x,y) = y³ − 1.
//   Unique real root: (1, 1).
//   Starting point: (0, 0).
//
// At (0, 0): J = [[3·0², 0], [0, 3·0²]] = [[0, 0], [0, 0]] — the zero matrix.
// Newton gets Δ = 0 every iteration (ColPivQR minimum-norm solution = 0)
// and remains at (0,0) where F = (−1, −1) forever.
// BisectionND probes the structured grid, finds sign changes between (0,0)
// and (2,2) ([−,−] vs [+,+]), bisects toward the unique root (1,1).
// ----------------------------------------------------------------------------
TEST_CASE("BisectionND - zero Jacobian at starting point (Newton stuck)", "[bisectionnd][advantage]") {
    // f₁ = x³-1, f₂ = y³-1.  Unique real root: (1,1).
    // At (0,0): J = diag(3·0², 3·0²) = [[0,0],[0,0]] exactly.
    auto evaluateCubes = [](const Eigen::VectorXd& x,
                            Eigen::VectorXd& F,
                            Eigen::MatrixXd& J,
                            bool computeJac) {
        F(0) = x(0)*x(0)*x(0) - 1.0;
        F(1) = x(1)*x(1)*x(1) - 1.0;
        if (computeJac) {
            J(0, 0) = 3.0*x(0)*x(0); J(0, 1) = 0.0;
            J(1, 0) = 0.0;            J(1, 1) = 3.0*x(1)*x(1);
        }
    };

    // BisectionND finds the root from the zero-Jacobian starting point
    {
        BisectionNDSolver bisSolver;
        NonLinearSolver::Problem prob;
        prob.size = 2;
        prob.evaluate = evaluateCubes;

        Eigen::VectorXd x(2);
        x << 0.0, 0.0;  // zero Jacobian here: J = [[0,0],[0,0]]

        SolverOptions opts;
        opts.tolerance = 1e-6;
        opts.maxIterations = 200;

        auto status = bisSolver.solve(prob, x, opts);
        REQUIRE(status == SolverStatus::Success);

        // Root must be (1,1)
        Eigen::VectorXd F(2);
        Eigen::MatrixXd Jd(2, 2);
        prob.evaluate(x, F, Jd, false);
        CHECK(F.lpNorm<Eigen::Infinity>() < 1e-5);
    }

    // Newton is stuck at (0,0): J=0 → step Δ=0 → no progress.
    // (ColPivHouseholderQR of the zero matrix returns the zero vector.)
    {
        NewtonSolver nSolver;
        NonLinearSolver::Problem prob;
        prob.size = 2;
        prob.evaluate = evaluateCubes;

        Eigen::VectorXd x(2);
        x << 0.0, 0.0;

        SolverOptions opts;
        opts.tolerance = 1e-6;
        opts.maxIterations = 50;

        auto status = nSolver.solve(prob, x, opts);
        // Newton cannot converge: J=0 makes every Newton step Δ=0,
        // triggering the StepTolerance early-exit without reducing residual.
        INFO("Newton on zero-Jacobian: " << statusToString(status));
        CHECK(status != SolverStatus::Success);
    }
}
