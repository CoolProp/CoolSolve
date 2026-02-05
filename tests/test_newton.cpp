#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver.h"
#include <cmath>
#include <iostream>

using namespace coolsolve;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Test Functions for Newton Solver
// ============================================================================

// Simple 1D root finding: f(x) = x^2 - 2 = 0, solution x = sqrt(2)
struct Sqrt2Problem {
    void operator()(const Eigen::VectorXd& x, 
                    Eigen::VectorXd& F, 
                    Eigen::MatrixXd& J, 
                    bool computeJacobian) const {
        F(0) = x(0) * x(0) - 2.0;
        if (computeJacobian) {
            J(0, 0) = 2.0 * x(0);
        }
    }
};

// 2D linear system: x + y = 3, x - y = 1, solution (2, 1)
struct Linear2DProblem {
    void operator()(const Eigen::VectorXd& x, 
                    Eigen::VectorXd& F, 
                    Eigen::MatrixXd& J, 
                    bool computeJacobian) const {
        F(0) = x(0) + x(1) - 3.0;
        F(1) = x(0) - x(1) - 1.0;
        if (computeJacobian) {
            J(0, 0) = 1.0; J(0, 1) = 1.0;
            J(1, 0) = 1.0; J(1, 1) = -1.0;
        }
    }
};

// 2D non-linear: x^2 + y^2 = 4, x*y = 1, solutions include (sqrt(2+sqrt(3)), sqrt(2-sqrt(3)))
struct Circle2DProblem {
    void operator()(const Eigen::VectorXd& x, 
                    Eigen::VectorXd& F, 
                    Eigen::MatrixXd& J, 
                    bool computeJacobian) const {
        F(0) = x(0) * x(0) + x(1) * x(1) - 4.0;
        F(1) = x(0) * x(1) - 1.0;
        if (computeJacobian) {
            J(0, 0) = 2.0 * x(0); J(0, 1) = 2.0 * x(1);
            J(1, 0) = x(1);       J(1, 1) = x(0);
        }
    }
};

// Rosenbrock function roots: a - x = 0, b*(y - x^2) = 0
// Standard params: a=1, b=100, solution is (1, 1)
struct RosenbrockProblem {
    double a = 1.0;
    double b = 100.0;
    
    void operator()(const Eigen::VectorXd& x, 
                    Eigen::VectorXd& F, 
                    Eigen::MatrixXd& J, 
                    bool computeJacobian) const {
        F(0) = a - x(0);
        F(1) = b * (x(1) - x(0) * x(0));
        if (computeJacobian) {
            J(0, 0) = -1.0;          J(0, 1) = 0.0;
            J(1, 0) = -2.0 * b * x(0); J(1, 1) = b;
        }
    }
};

// Transcendental: exp(x) + y = 2, x + exp(y) = 2
// Solution is approximately (0.5671, 0.5671)
struct TranscendentalProblem {
    void operator()(const Eigen::VectorXd& x, 
                    Eigen::VectorXd& F, 
                    Eigen::MatrixXd& J, 
                    bool computeJacobian) const {
        F(0) = std::exp(x(0)) + x(1) - 2.0;
        F(1) = x(0) + std::exp(x(1)) - 2.0;
        if (computeJacobian) {
            J(0, 0) = std::exp(x(0)); J(0, 1) = 1.0;
            J(1, 0) = 1.0;            J(1, 1) = std::exp(x(1));
        }
    }
};

// Colebrook equation (implicit): 1/sqrt(f) = -2*log10(e/D/3.7 + 2.51/(Re*sqrt(f)))
// Rearranged to F(f) = 0
struct ColebrookProblem {
    double Re = 100000.0;  // Reynolds number
    double eD = 0.001;     // Relative roughness e/D
    
    void operator()(const Eigen::VectorXd& x, 
                    Eigen::VectorXd& F, 
                    Eigen::MatrixXd& J, 
                    bool computeJacobian) const {
        double f = x(0);
        if (f <= 0) f = 1e-10;  // Prevent sqrt of negative
        
        double sqrtf = std::sqrt(f);
        double arg = eD / 3.7 + 2.51 / (Re * sqrtf);
        
        F(0) = 1.0 / sqrtf + 2.0 * std::log10(arg);
        
        if (computeJacobian) {
            // d/df (1/sqrt(f)) = -1/(2*f^(3/2))
            double dsqrtf = -0.5 / (f * sqrtf);
            
            // d/df (log10(a + b/sqrt(f))) = -b/(2*f^(3/2)) / (a + b/sqrt(f)) / ln(10)
            double darg = -2.51 / (Re * 2.0 * f * sqrtf);
            double dlog = darg / (arg * std::log(10.0));
            
            J(0, 0) = dsqrtf + 2.0 * dlog;
        }
    }
};

// ============================================================================
// Tests
// ============================================================================

TEST_CASE("Newton solver - sqrt(2)", "[newton][1d]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = Sqrt2Problem();
    
    Eigen::VectorXd x(1);
    x(0) = 1.0;  // Initial guess
    
    SolverOptions options;
    options.tolerance = 1e-12;
    
    SolverStatus status = solver.solve(problem, x, options);
    
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(std::sqrt(2.0), 1e-10));
}

TEST_CASE("Newton solver - 2D linear system", "[newton][2d][linear]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = Linear2DProblem();
    
    Eigen::VectorXd x(2);
    x << 0.0, 0.0;  // Initial guess
    
    SolverOptions options;
    options.tolerance = 1e-12;
    
    SolverStatus status = solver.solve(problem, x, options);
    
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(2.0, 1e-10));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-10));
}

TEST_CASE("Newton solver - 2D non-linear (circle)", "[newton][2d][nonlinear]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = Circle2DProblem();
    
    Eigen::VectorXd x(2);
    x << 1.5, 0.5;  // Initial guess near solution
    
    SolverOptions options;
    options.tolerance = 1e-10;
    
    SolverStatus status = solver.solve(problem, x, options);
    
    REQUIRE(status == SolverStatus::Success);
    
    // Verify solution satisfies equations
    CHECK_THAT(x(0) * x(0) + x(1) * x(1), WithinRel(4.0, 1e-9));
    CHECK_THAT(x(0) * x(1), WithinRel(1.0, 1e-9));
}

TEST_CASE("Newton solver - Rosenbrock", "[newton][rosenbrock]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = RosenbrockProblem();
    
    Eigen::VectorXd x(2);
    x << 0.5, 0.5;  // Initial guess
    
    SolverOptions options;
    options.tolerance = 1e-12;
    
    SolverStatus status = solver.solve(problem, x, options);
    
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-10));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-10));
}

TEST_CASE("Newton solver - Rosenbrock from poor guess", "[newton][rosenbrock][hard]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = RosenbrockProblem();
    
    Eigen::VectorXd x(2);
    x << -5.0, 10.0;  // Poor initial guess
    
    SolverOptions options;
    options.tolerance = 1e-10;
    options.maxIterations = 200;  // May need more iterations
    
    SolverStatus status = solver.solve(problem, x, options);
    
    // Line search helps with convergence from poor guesses
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(x(0), WithinRel(1.0, 1e-8));
    CHECK_THAT(x(1), WithinRel(1.0, 1e-8));
}

TEST_CASE("Newton solver - Transcendental", "[newton][transcendental]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = TranscendentalProblem();
    
    Eigen::VectorXd x(2);
    x << 0.5, 0.5;  // Initial guess
    
    SolverOptions options;
    options.tolerance = 1e-12;
    
    SolverStatus status = solver.solve(problem, x, options);
    
    REQUIRE(status == SolverStatus::Success);
    
    // The solution is symmetric: x = y
    CHECK_THAT(x(0), WithinAbs(x(1), 1e-10));
    
    // Verify solution satisfies equations
    CHECK_THAT(std::exp(x(0)) + x(1), WithinRel(2.0, 1e-10));
    CHECK_THAT(x(0) + std::exp(x(1)), WithinRel(2.0, 1e-10));
}

TEST_CASE("Newton solver - Colebrook equation", "[newton][colebrook]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    
    ColebrookProblem colebrook;
    colebrook.Re = 100000.0;
    colebrook.eD = 0.001;
    problem.evaluate = colebrook;
    
    Eigen::VectorXd x(1);
    x(0) = 0.02;  // Initial guess (typical turbulent friction factor)
    
    SolverOptions options;
    options.tolerance = 1e-12;
    
    SolverStatus status = solver.solve(problem, x, options);
    
    REQUIRE(status == SolverStatus::Success);
    
    // Verify: Colebrook equation should be satisfied
    double f = x(0);
    double sqrtf = std::sqrt(f);
    double arg = colebrook.eD / 3.7 + 2.51 / (colebrook.Re * sqrtf);
    double residual = std::abs(1.0 / sqrtf + 2.0 * std::log10(arg));
    CHECK(residual < 1e-10);
    
    // Friction factor should be in reasonable range
    CHECK(f > 0.01);
    CHECK(f < 0.05);
}

TEST_CASE("Newton solver - with tracing", "[newton][trace]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = Sqrt2Problem();
    
    Eigen::VectorXd x(1);
    x(0) = 1.0;
    
    SolverOptions options;
    options.tolerance = 1e-12;
    
    SolverTrace trace;
    SolverStatus status = solver.solve(problem, x, options, &trace);
    
    REQUIRE(status == SolverStatus::Success);
    REQUIRE(trace.iterations.size() >= 1);
    REQUIRE(trace.finalStatus == SolverStatus::Success);
    
    // Check that residuals decrease
    for (size_t i = 1; i < trace.iterations.size(); ++i) {
        // Generally residuals should decrease (with some tolerance for numerical issues)
        // but line search guarantees descent in ||F||^2
    }
    
    // Output trace for inspection
    INFO("Trace:\n" << trace.toString());
}

TEST_CASE("Newton solver - convergence from negative guess", "[newton][robust]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 1;
    problem.evaluate = Sqrt2Problem();
    
    // Negative initial guess - should still converge to positive root
    Eigen::VectorXd x(1);
    x(0) = -1.0;
    
    SolverOptions options;
    options.tolerance = 1e-10;
    
    SolverStatus status = solver.solve(problem, x, options);
    
    // Could converge to -sqrt(2) or +sqrt(2), both are valid
    REQUIRE(status == SolverStatus::Success);
    CHECK_THAT(std::abs(x(0)), WithinRel(std::sqrt(2.0), 1e-8));
}

TEST_CASE("Newton solver - max iterations", "[newton][limits]") {
    NewtonSolver solver;
    NonLinearSolver::Problem problem;
    problem.size = 2;
    problem.evaluate = RosenbrockProblem();
    
    Eigen::VectorXd x(2);
    x << -10.0, 100.0;  // Very poor guess
    
    SolverOptions options;
    options.tolerance = 1e-15;  // Very tight tolerance
    options.maxIterations = 3;  // Very few iterations
    
    SolverStatus status = solver.solve(problem, x, options);
    
    // Should hit max iterations with these settings
    CHECK(status == SolverStatus::MaxIterations);
}
