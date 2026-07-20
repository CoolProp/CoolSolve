#include "integrator_test_util.h"
#include "coolsolve/integral/integrator.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace coolsolve;
using namespace coolsolve::test;
using Catch::Matchers::WithinRel;

// ============================================================================
// TDD tests for the classic 4th-order Runge-Kutta integrator.
// ============================================================================

TEST_CASE("RK4: dy/dt = -y recovers e^{-t}", "[integrator][rk4]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK4);
    REQUIRE(std::string(ig->name()) == "RK4");
    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 1.0;

    // 4th-order: halving h cuts error ~16x.
    double err10  = std::abs(marchFixed(*ig, 0.0, 4.0, y0, f, 10, opt)(0)  - std::exp(-4.0));
    double err20  = std::abs(marchFixed(*ig, 0.0, 4.0, y0, f, 20, opt)(0)  - std::exp(-4.0));
    INFO("RK4 err10=" << err10 << " err20=" << err20 << " ratio=" << (err10 / err20));
    CHECK(err20 < err10);
    CHECK_THAT(err10 / err20, WithinRel(16.0, 0.25));  // ~O(h^4)

    // With 100 steps RK4 is essentially exact.
    double yEnd = marchFixed(*ig, 0.0, 4.0, y0, f, 100, opt)(0);
    CHECK_THAT(yEnd, WithinRel(std::exp(-4.0), 1e-6));
}

TEST_CASE("RK4: dy/dt = t^2 integrates polynomial closely", "[integrator][rk4]") {
    RHSFunction f = [](double t, const Eigen::VectorXd& /*y*/) {
        Eigen::VectorXd r(1); r(0) = t * t; return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK4);
    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 0.0;

    // RK4 is exact (to round-off) for polynomials of degree <= 4 in a single
    // step's internal quadrature; the global result is very tight.
    double yEnd = marchFixed(*ig, 0.0, 1.0, y0, f, 10, opt)(0);
    CHECK_THAT(yEnd, WithinRel(1.0 / 3.0, 1e-9));
}

TEST_CASE("RK4: coupled harmonic oscillator", "[integrator][rk4]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(2);
        r(0) = y(1);
        r(1) = -y(0);
        return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK4);
    IntegratorOptions opt;
    Eigen::VectorXd y0(2); y0(0) = 0.0; y0(1) = 1.0;

    Eigen::VectorXd yEnd = marchFixed(*ig, 0.0, 1.0, y0, f, 20, opt);
    CHECK_THAT(yEnd(0), WithinRel(std::sin(1.0), 1e-6));
    CHECK_THAT(yEnd(1), WithinRel(std::cos(1.0), 1e-6));
}

TEST_CASE("RK4: RHS evaluation count is exactly 4 per step", "[integrator][rk4]") {
    // Mirrors the contributing-guide performance discipline: each method must
    // document and respect its eval budget.
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK4);
    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 1.0;
    StepResult r = ig->step(0.0, y0, f, 0.1, opt);
    CHECK(r.rhsEvaluations == 4);
}
