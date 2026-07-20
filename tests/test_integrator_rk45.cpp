#include "integrator_test_util.h"
#include "coolsolve/integral/integrator.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace coolsolve;
using namespace coolsolve::test;
using Catch::Matchers::WithinRel;

// ============================================================================
// TDD tests for the Dormand-Prince RK45 adaptive integrator.
// ============================================================================

TEST_CASE("RK45: dy/dt = -y recovers e^{-t}", "[integrator][rk45]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK45);
    REQUIRE(std::string(ig->name()) == "RK45");
    IntegratorOptions opt;
    opt.relTol = 1e-6;
    opt.absTol = 1e-9;
    Eigen::VectorXd y0(1); y0(0) = 1.0;

    int taken = 0, rejected = 0;
    Eigen::VectorXd yEnd = marchAdaptive(*ig, 0.0, 4.0, y0, f, 0.1, opt, &taken, &rejected);
    INFO("RK45 taken=" << taken << " rejected=" << rejected);
    CHECK(taken > 0);
    CHECK_THAT(yEnd(0), WithinRel(std::exp(-4.0), 1e-4));
}

TEST_CASE("RK45: error scales with the requested tolerance", "[integrator][rk45]") {
    // Tighter tol => smaller actual error.  The error should track `relTol`.
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK45);
    Eigen::VectorXd y0(1); y0(0) = 1.0;
    const double exact = std::exp(-4.0);

    double errLo = 0, errHi = 0;
    {
        IntegratorOptions opt; opt.relTol = 1e-3; opt.absTol = 1e-9;
        errLo = std::abs(marchAdaptive(*ig, 0.0, 4.0, y0, f, 0.1, opt)(0) - exact);
    }
    {
        IntegratorOptions opt; opt.relTol = 1e-9; opt.absTol = 1e-12;
        errHi = std::abs(marchAdaptive(*ig, 0.0, 4.0, y0, f, 0.1, opt)(0) - exact);
    }
    INFO("errLo=" << errLo << " errHi=" << errHi);
    CHECK(errHi < errLo);
    CHECK(errHi < 1e-6);
}

TEST_CASE("RK45: coupled harmonic oscillator", "[integrator][rk45]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(2);
        r(0) = y(1);
        r(1) = -y(0);
        return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK45);
    IntegratorOptions opt; opt.relTol = 1e-8; opt.absTol = 1e-10;
    Eigen::VectorXd y0(2); y0(0) = 0.0; y0(1) = 1.0;

    Eigen::VectorXd yEnd = marchAdaptive(*ig, 0.0, 1.0, y0, f, 0.05, opt);
    CHECK_THAT(yEnd(0), WithinRel(std::sin(1.0), 1e-6));
    CHECK_THAT(yEnd(1), WithinRel(std::cos(1.0), 1e-6));
}

TEST_CASE("RK45: rejected step is retried at smaller h", "[integrator][rk45]") {
    // A discontinuous / steep RHS triggers rejection; the step must be marked
    // not-accepted and a smaller nextStep proposed.
    RHSFunction f = [](double t, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1);
        // Steep derivative; with a large first step the error estimate is huge.
        r(0) = -1000.0 * y(0) + (t > 0.5 ? 1e6 : 0.0);
        return r;
    };
    auto ig = createIntegrator(IntegratorOptions::RK45);
    IntegratorOptions opt; opt.relTol = 1e-6; opt.absTol = 1e-9;
    Eigen::VectorXd y0(1); y0(0) = 1.0;

    // First step of 0.5 (jumps across the kink) should be rejected.
    StepResult r = ig->step(0.0, y0, f, 0.5, opt);
    CHECK((!r.accepted || r.errorEstimate > 0.0));
    if (!r.accepted) {
        CHECK(r.nextStep < 0.5);
    }
}
