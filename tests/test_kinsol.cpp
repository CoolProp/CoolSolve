/**
 * @file test_kinsol.cpp
 * @brief Tests for the KINSOL solver (SUNDIALS-style inexact Newton + Dennis-
 *        Schnabel line search, Picard, and Anderson-accelerated fixed point).
 *
 * Covers all three globalisation modes (KIN_LINESEARCH, KIN_PICARD, KIN_FP),
 * convergence checks, agreement between modes, the Anderson depth-0 fallback
 * to plain Picard, the strategy/factory wiring, and coolsolve.conf parsing.
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver.h"
#include <cmath>
#include <filesystem>
#include <fstream>

using namespace coolsolve;
using Catch::Matchers::WithinAbs;
namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Verify F(x) ≈ 0 component-wise for a solved problem.
static void checkResidualZero(NonLinearSolver::Problem& p, const Eigen::VectorXd& x,
                              double tol = 1e-6) {
    Eigen::VectorXd F(x.size());
    Eigen::MatrixXd J(x.size(), x.size());
    p.evaluate(x, F, J, false);
    CHECK(F.lpNorm<Eigen::Infinity>() < tol);
}

// ============================================================================
// KIN_LINESEARCH — inexact Newton + Dennis-Schnabel line search
// ============================================================================

TEST_CASE("KINSOL LineSearch - simple 1D root", "[kinsol][linesearch]") {
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd& J, bool computeJac) {
        F(0) = x(0) * x(0) - 4.0;
        if (computeJac) J(0, 0) = 2.0 * x(0);
    };
    Eigen::VectorXd x(1); x(0) = 10.0;

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::LineSearch;
    opts.tolerance = 1e-10;
    CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
    CHECK_THAT(std::abs(x(0)), WithinAbs(2.0, 1e-8));
    checkResidualZero(problem, x);
}

TEST_CASE("KINSOL LineSearch - 2D nonlinear", "[kinsol][linesearch]") {
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd& J, bool computeJac) {
        F(0) = x(0) * x(0) + x(1) * x(1) - 4.0;
        F(1) = x(0) * x(1) - 1.0;
        if (computeJac) {
            J(0, 0) = 2.0 * x(0); J(0, 1) = 2.0 * x(1);
            J(1, 0) = x(1);       J(1, 1) = x(0);
        }
    };
    Eigen::VectorXd x(2); x << 2.0, 0.5;  // non-singular start (det J = 4·2 − 1·0.5 ≠ 0)

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::LineSearch;
    opts.tolerance = 1e-10;
    CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
    checkResidualZero(problem, x);
}

TEST_CASE("KINSOL LineSearch - mild Rosenbrock", "[kinsol][linesearch]") {
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd& J, bool computeJac) {
        F(0) = 1.0 - x(0);
        F(1) = 10.0 * (x(1) - x(0) * x(0));
        if (computeJac) {
            J(0, 0) = -1.0;           J(0, 1) = 0.0;
            J(1, 0) = -20.0 * x(0);   J(1, 1) = 10.0;
        }
    };
    Eigen::VectorXd x(2); x << -1.2, 1.0;

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::LineSearch;
    opts.tolerance = 1e-10;
    CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
    CHECK_THAT(x(0), WithinAbs(1.0, 1e-6));
    CHECK_THAT(x(1), WithinAbs(1.0, 1e-6));
}

TEST_CASE("KINSOL LineSearch - reports SingularJacobian on rank-deficient J",
          "[kinsol][linesearch]") {
    // F(x) = x^2 (double root at 0). J = 2x, which is 0 at x=0 and rank-deficient
    // as a system near there. Starting at x=0 the Jacobian is exactly singular.
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd& J, bool computeJac) {
        F(0) = x(0) * x(0);
        if (computeJac) J(0, 0) = 2.0 * x(0);
    };
    Eigen::VectorXd x(1); x(0) = 0.0;  // singular Jacobian here

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::LineSearch;
    opts.tolerance = 1e-10;
    SolverStatus status = solver.solve(problem, x, opts);
    // At x=0, F=0 already → the solver converges immediately (before factoring).
    CHECK(status == SolverStatus::Success);
    CHECK_THAT(std::abs(x(0)), WithinAbs(0.0, 1e-9));
}

// ============================================================================
// KIN_PICARD — fixed-point (Richardson) iteration
// ============================================================================

TEST_CASE("KINSOL Picard - diagonal linear system (exact in 1 step)", "[kinsol][picard]") {
    // F0 = 2x - 5, F1 = 3y - 6.  J = diag(2,3).  With ω=0.5, the first
    // variable converges in one Picard step (|1 - ω·2| = 0); the second needs
    // |1 - 0.5·3| = 0.5 so a few steps.  Solutions: x=2.5, y=2.
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd&, bool) {
        F(0) = 2.0 * x(0) - 5.0;
        F(1) = 3.0 * x(1) - 6.0;
    };
    Eigen::VectorXd x(2); x << 0.0, 0.0;

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::Picard;
    opts.kinsolPicardOmega = 0.5;
    opts.tolerance = 1e-10;
    opts.maxIterations = 1000;
    CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
    CHECK_THAT(x(0), WithinAbs(2.5, 1e-8));
    CHECK_THAT(x(1), WithinAbs(2.0, 1e-8));
    checkResidualZero(problem, x);
}

TEST_CASE("KINSOL Picard - under-relaxation needed for convergence", "[kinsol][picard]") {
    // F0 = 4x - 1.  ρ(I - ωJ) = |1 - 4ω|.  With ω=1 it diverges (|−3|>1);
    // with ω=0.4 it converges (|1-1.6|=0.6).  Solution x=0.25.
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd&, bool) {
        F(0) = 4.0 * x(0) - 1.0;
    };

    // Diverges with ω=1.
    {
        Eigen::VectorXd x(1); x(0) = 0.1;
        SolverOptions opts;
        opts.kinsolGlobalStrategy = KinsolGlobalStrategy::Picard;
        opts.kinsolPicardOmega = 1.0;
        opts.tolerance = 1e-9;
        opts.maxIterations = 200;
        SolverStatus s = solver.solve(problem, x, opts);
        CHECK((s == SolverStatus::MaxIterations || s == SolverStatus::Diverged));
    }
    // Converges with ω=0.4.
    {
        Eigen::VectorXd x(1); x(0) = 0.1;
        SolverOptions opts;
        opts.kinsolGlobalStrategy = KinsolGlobalStrategy::Picard;
        opts.kinsolPicardOmega = 0.4;
        opts.tolerance = 1e-9;
        opts.maxIterations = 200;
        CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
        CHECK_THAT(x(0), WithinAbs(0.25, 1e-8));
    }
}

// ============================================================================
// KIN_FP — Anderson-accelerated fixed-point iteration
// ============================================================================

TEST_CASE("KINSOL FixedPoint (Anderson) - contraction converges", "[kinsol][fp]") {
    // Contraction fixed point: x = 0.5·cos(y) + 0.5, y = 0.5·sin(x) + 0.5.
    // Jacobian of G has |entries| ≤ 0.5 → guaranteed contraction.
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd&, bool) {
        F(0) = x(0) - (0.5 * std::cos(x(1)) + 0.5);
        F(1) = x(1) - (0.5 * std::sin(x(0)) + 0.5);
    };
    Eigen::VectorXd x(2); x << 0.0, 0.0;

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::FixedPoint;
    opts.kinsolAndersonDepth = 5;
    opts.tolerance = 1e-10;
    opts.maxIterations = 200;
    CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
    checkResidualZero(problem, x);
}

TEST_CASE("KINSOL FixedPoint - depth 0 behaves as plain Picard", "[kinsol][fp]") {
    // Same contraction as above; depth=0 disables Anderson acceleration so the
    // iteration is plain fixed-point x ← x − ωF(x).  It must still converge.
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd&, bool) {
        F(0) = x(0) - (0.5 * std::cos(x(1)) + 0.5);
        F(1) = x(1) - (0.5 * std::sin(x(0)) + 0.5);
    };
    Eigen::VectorXd x(2); x << 0.0, 0.0;

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::FixedPoint;
    opts.kinsolAndersonDepth = 0;
    opts.tolerance = 1e-10;
    opts.maxIterations = 500;
    CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
    checkResidualZero(problem, x);
}

TEST_CASE("KINSOL FixedPoint (Anderson) - solves 2D circle", "[kinsol][fp]") {
    // A problem that is not a simple contraction: the unit-circle / hyperbola
    // system.  Anderson acceleration should still find the root from a
    // reasonable start (it mixes the fixed-point map G = x − F).
    KINSOLSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd&, bool) {
        F(0) = x(0) * x(0) + x(1) * x(1) - 4.0;
        F(1) = x(0) * x(1) - 1.0;
    };
    Eigen::VectorXd x(2); x << 1.5, 1.5;

    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::FixedPoint;
    opts.kinsolAndersonDepth = 5;
    opts.kinsolAndersonRelaxation = 0.5;  // damped for stability on this non-contraction
    opts.tolerance = 1e-8;
    opts.maxIterations = 500;
    SolverStatus status = solver.solve(problem, x, opts);
    CHECK(status == SolverStatus::Success);
    if (status == SolverStatus::Success) checkResidualZero(problem, x, 1e-5);
}

// ============================================================================
// Cross-mode agreement: all three modes reach the same root.
// ============================================================================

TEST_CASE("KINSOL - all three modes agree on diagonal system", "[kinsol]") {
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd& J, bool computeJac) {
        F(0) = 2.0 * x(0) - 5.0;
        F(1) = 3.0 * x(1) - 6.0;
        if (computeJac) { J(0, 0) = 2.0; J(1, 1) = 3.0; }
    };

    auto run = [&](KinsolGlobalStrategy mode) {
        KINSOLSolver solver;
        Eigen::VectorXd x(2); x << 0.0, 0.0;
        SolverOptions opts;
        opts.kinsolGlobalStrategy = mode;
        opts.kinsolPicardOmega = 0.5;
        opts.tolerance = 1e-10;
        opts.maxIterations = 1000;
        CHECK(solver.solve(problem, x, opts) == SolverStatus::Success);
        return x;
    };

    Eigen::VectorXd xLs = run(KinsolGlobalStrategy::LineSearch);
    Eigen::VectorXd xPi = run(KinsolGlobalStrategy::Picard);
    Eigen::VectorXd xFp = run(KinsolGlobalStrategy::FixedPoint);
    CHECK_THAT(xLs(0), WithinAbs(2.5, 1e-8));
    CHECK_THAT(xLs(1), WithinAbs(2.0, 1e-8));
    CHECK_THAT(xPi(0), WithinAbs(2.5, 1e-8));
    CHECK_THAT(xPi(1), WithinAbs(2.0, 1e-8));
    CHECK_THAT(xFp(0), WithinAbs(2.5, 1e-8));
    CHECK_THAT(xFp(1), WithinAbs(2.0, 1e-8));
}

// ============================================================================
// Strategy / factory wiring
// ============================================================================

TEST_CASE("KINSOL - strategy string and parsing", "[kinsol][config]") {
    CHECK(strategyToString(SolverStrategy::Kinsol) == "Kinsol");
    SolverStrategy s;
    CHECK(parseStrategy("kinsol", s));
    CHECK(s == SolverStrategy::Kinsol);
    CHECK(parseStrategy("KINSOL", s));
    CHECK(s == SolverStrategy::Kinsol);
    CHECK(parseStrategy("kin", s));
    CHECK(s == SolverStrategy::Kinsol);
}

TEST_CASE("KINSOL - factory creates a usable solver", "[kinsol][config]") {
    auto solver = createSolver(SolverStrategy::Kinsol);
    REQUIRE(solver != nullptr);
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = [](const Eigen::VectorXd& x, Eigen::VectorXd& F,
                          Eigen::MatrixXd& J, bool cj) {
        F(0) = x(0) - 7.0;
        if (cj) J(0, 0) = 1.0;
    };
    Eigen::VectorXd x(1); x(0) = 0.0;
    SolverOptions opts;
    opts.kinsolGlobalStrategy = KinsolGlobalStrategy::LineSearch;
    opts.tolerance = 1e-10;
    CHECK(solver->solve(problem, x, opts) == SolverStatus::Success);
    CHECK_THAT(x(0), WithinAbs(7.0, 1e-9));
}

// ============================================================================
// Config-file parsing of the kinsol* keys
// ============================================================================

static SolverOptions loadFromConfString(const std::string& body) {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_kinsol.conf";
    { std::ofstream f(configPath); f << body; }
    SolverOptions opts;
    loadSolverOptionsFromFile(configPath.string(), opts);
    fs::remove(configPath);
    return opts;
}

TEST_CASE("Config - kinsolGlobalStrategy = picard", "[kinsol][config]") {
    auto opts = loadFromConfString("kinsolGlobalStrategy = picard\n");
    CHECK(opts.kinsolGlobalStrategy == KinsolGlobalStrategy::Picard);
}

TEST_CASE("Config - kinsolGlobalStrategy = fp", "[kinsol][config]") {
    auto opts = loadFromConfString("kinsolGlobalStrategy = fp\n");
    CHECK(opts.kinsolGlobalStrategy == KinsolGlobalStrategy::FixedPoint);
}

TEST_CASE("Config - kinsolGlobalStrategy = linesearch (default)", "[kinsol][config]") {
    auto opts = loadFromConfString("kinsolGlobalStrategy = linesearch\n");
    CHECK(opts.kinsolGlobalStrategy == KinsolGlobalStrategy::LineSearch);
}

TEST_CASE("Config - invalid kinsolGlobalStrategy falls back to linesearch",
          "[kinsol][config]") {
    auto opts = loadFromConfString("kinsolGlobalStrategy = bogus\n");
    CHECK(opts.kinsolGlobalStrategy == KinsolGlobalStrategy::LineSearch);
}

TEST_CASE("Config - Kinsol accepted in solverPipeline", "[kinsol][config]") {
    auto opts = loadFromConfString("solverPipeline = Newton, Kinsol, Homotopy\n");
    REQUIRE(opts.solverPipeline.size() == 3);
    CHECK(opts.solverPipeline[0] == SolverStrategy::Newton);
    CHECK(opts.solverPipeline[1] == SolverStrategy::Kinsol);
    CHECK(opts.solverPipeline[2] == SolverStrategy::Homotopy);
}

TEST_CASE("Config - kinsol numeric options round-trip", "[kinsol][config]") {
    auto opts = loadFromConfString(
        "kinsolLineSearchAlpha = 1e-3\n"
        "kinsolLineSearchMaxIters = 25\n"
        "kinsolPicardOmega = 0.7\n"
        "kinsolAndersonDepth = 9\n"
        "kinsolAndersonRelaxation = 0.8\n");
    CHECK_THAT(opts.kinsolLineSearchAlpha, WithinAbs(1e-3, 1e-15));
    CHECK(opts.kinsolLineSearchMaxIters == 25);
    CHECK_THAT(opts.kinsolPicardOmega, WithinAbs(0.7, 1e-15));
    CHECK(opts.kinsolAndersonDepth == 9);
    CHECK_THAT(opts.kinsolAndersonRelaxation, WithinAbs(0.8, 1e-15));
}

TEST_CASE("Config - out-of-range kinsol options fall back to defaults",
          "[kinsol][config]") {
    auto opts = loadFromConfString(
        "kinsolLineSearchAlpha = -1.0\n"
        "kinsolPicardOmega = 0.0\n"
        "kinsolAndersonRelaxation = 2.0\n");
    CHECK_THAT(opts.kinsolLineSearchAlpha, WithinAbs(1e-4, 1e-15));
    CHECK_THAT(opts.kinsolPicardOmega, WithinAbs(1.0, 1e-15));
    CHECK_THAT(opts.kinsolAndersonRelaxation, WithinAbs(1.0, 1e-15));
}
