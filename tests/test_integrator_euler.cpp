#include "integrator_test_util.h"
#include "coolsolve/integral/integrator.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace coolsolve;
using namespace coolsolve::test;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ============================================================================
// TDD tests for the Euler (explicit + implicit) integrators.
//
// These pin down the numerics BEFORE any parser/IR coupling, against
// analytical solutions.
// ============================================================================

static const double tol5pct = 0.05;   // generous for 1st-order methods
static const double tol1pct = 0.01;

TEST_CASE("Euler explicit: dy/dt = -y recovers e^{-t}", "[integrator][euler]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    auto ig = createIntegrator(IntegratorOptions::EulerExplicit);
    REQUIRE(std::string(ig->name()) == "EulerExplicit");

    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 1.0;

    // Convergence order check: error ~ h.  Halving h halves the error.
    double err100 = std::abs(marchFixed(*ig, 0.0, 4.0, y0, f, 100, opt)(0) - std::exp(-4.0));
    double err200 = std::abs(marchFixed(*ig, 0.0, 4.0, y0, f, 200, opt)(0) - std::exp(-4.0));
    INFO("err100=" << err100 << " err200=" << err200 << " ratio=" << (err100 / err200));
    // ~2x reduction => first-order method.
    CHECK(err200 < err100);
    CHECK_THAT(err100 / err200, WithinRel(2.0, 0.15));

    // Coarse accuracy.
    double yEnd = marchFixed(*ig, 0.0, 4.0, y0, f, 1000, opt)(0);
    CHECK_THAT(yEnd, WithinRel(std::exp(-4.0), tol5pct));
}

TEST_CASE("Euler explicit: dy/dt = t^2 integrates polynomial exactly-ish", "[integrator][euler]") {
    RHSFunction f = [](double t, const Eigen::VectorXd& /*y*/) {
        Eigen::VectorXd r(1); r(0) = t * t; return r;
    };
    auto ig = createIntegrator(IntegratorOptions::EulerExplicit);
    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 0.0;

    // Analytical: y(1) = 1/3.  Euler is first order; with 1000 steps we are close.
    double yEnd = marchFixed(*ig, 0.0, 1.0, y0, f, 1000, opt)(0);
    CHECK_THAT(yEnd, WithinRel(1.0 / 3.0, tol1pct));
}

TEST_CASE("Euler explicit: coupled harmonic oscillator", "[integrator][euler]") {
    // y' = z, z' = -y  =>  y(t) = sin(t), z(t) = cos(t), y(0)=0, z(0)=1.
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(2);
        r(0) = y(1);       // dy/dt = z
        r(1) = -y(0);      // dz/dt = -y
        return r;
    };
    auto ig = createIntegrator(IntegratorOptions::EulerExplicit);
    IntegratorOptions opt;
    Eigen::VectorXd y0(2); y0(0) = 0.0; y0(1) = 1.0;

    // Euler drifts on oscillatory problems; use many steps and a loose band.
    Eigen::VectorXd yEnd = marchFixed(*ig, 0.0, 1.0, y0, f, 10000, opt);
    CHECK_THAT(yEnd(0), WithinRel(std::sin(1.0), 0.05));
    CHECK_THAT(yEnd(1), WithinRel(std::cos(1.0), 0.05));
}

TEST_CASE("Euler implicit: stable on stiff dy/dt = -1000 y", "[integrator][euler]") {
    // Explicit Euler is unstable here unless h < 2/1000 = 0.002.
    // Implicit Euler is A-stable: any step size is stable.
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -1000.0 * y(0); return r;
    };
    Eigen::VectorXd y0(1); y0(0) = 1.0;
    const double tf = 0.01;  // y(tf) = e^{-10} ~ 4.5e-5

    IntegratorOptions opt;

    // Explicit Euler with h = 0.01 (1 step) blows up: |1 - 1000*0.01| = 9.
    auto igE = createIntegrator(IntegratorOptions::EulerExplicit);
    double yExpl = marchFixed(*igE, 0.0, tf, y0, f, 1, opt)(0);
    CHECK(std::abs(yExpl) > 1.0);  // diverged

    // Implicit Euler with the same large step stays bounded and positive.
    auto igI = createIntegrator(IntegratorOptions::EulerImplicit);
    double yImpl = marchFixed(*igI, 0.0, tf, y0, f, 1, opt)(0);
    CHECK(std::isfinite(yImpl));
    CHECK(yImpl >= 0.0);
    CHECK(yImpl < 1.0);  // decaying, not diverging
    // Analytical implicit-Euler step: y1 = y0/(1+1000 h) = 1/11 ~ 0.0909.
    CHECK_THAT(yImpl, WithinRel(1.0 / 11.0, tol1pct));
}

TEST_CASE("Euler implicit: dy/dt = -y recovers e^{-t} (1st order)", "[integrator][euler]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    auto ig = createIntegrator(IntegratorOptions::EulerImplicit);
    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 1.0;

    double err100 = std::abs(marchFixed(*ig, 0.0, 4.0, y0, f, 100, opt)(0) - std::exp(-4.0));
    double err200 = std::abs(marchFixed(*ig, 0.0, 4.0, y0, f, 200, opt)(0) - std::exp(-4.0));
    INFO("implicit err100=" << err100 << " err200=" << err200 << " ratio=" << (err100 / err200));
    CHECK(err200 < err100);
    CHECK_THAT(err100 / err200, WithinRel(2.0, 0.15));  // first order
}
