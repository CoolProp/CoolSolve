#include "integrator_test_util.h"
#include "coolsolve/integral/integrator.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace coolsolve;
using namespace coolsolve::test;
using Catch::Matchers::WithinRel;

// ============================================================================
// Phase 0 — TDD tests for Richardson extrapolation (fixed-step methods).
//
// Richardson cancels the leading O(h^p) error term, raising the effective
// order by one:  I_rich = (2^p · I_{h/2} − I_h) / (2^p − 1).
// We verify the *order gain* (the defining property) by checking the error
// convergence ratio when N doubles, not a hard absolute threshold.
// See docs/integral_table_plan.md §Phase 0 and §2.3.
// ============================================================================

TEST_CASE("Richardson: RK4 gains one order (4 -> 5)", "[richardson][integrator]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    const double exact = std::exp(-4.0);
    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 1.0;

    // Plain RK4: order 4, so halving h cuts error ~16x.
    auto rk4 = createIntegrator(IntegratorOptions::RK4);
    double e10 = std::abs(marchFixed(*rk4, 0.0, 4.0, y0, f, 10, opt)(0) - exact);
    double e20 = std::abs(marchFixed(*rk4, 0.0, 4.0, y0, f, 20, opt)(0) - exact);
    INFO("RK4 plain e10=" << e10 << " e20=" << e20 << " ratio=" << (e10 / e20));
    CHECK_THAT(e10 / e20, WithinRel(16.0, 0.30));

    // Richardson-wrapped RK4: order 5, halving h cuts error ~32x.
    auto rich = wrapRichardson(createIntegrator(IntegratorOptions::RK4), true);
    double r10 = std::abs(marchFixed(*rich, 0.0, 4.0, y0, f, 10, opt)(0) - exact);
    double r20 = std::abs(marchFixed(*rich, 0.0, 4.0, y0, f, 20, opt)(0) - exact);
    INFO("RK4+Richardson r10=" << r10 << " r20=" << r20 << " ratio=" << (r10 / r20));
    CHECK(r10 < e10);                  // strictly better at the same N
    CHECK_THAT(r10 / r20, WithinRel(32.0, 0.40));  // effective order 5
}

TEST_CASE("Richardson: disabled returns the base integrator", "[richardson][integrator]") {
    auto wrapped = wrapRichardson(createIntegrator(IntegratorOptions::RK4), false);
    // When disabled, wrapRichardson returns the base integrator unchanged.
    REQUIRE(std::string(wrapped->name()) == "RK4");
}

TEST_CASE("Richardson: Euler gains one order (1 -> 2)", "[richardson][integrator]") {
    RHSFunction f = [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1); r(0) = -y(0); return r;
    };
    const double exact = std::exp(-4.0);
    IntegratorOptions opt;
    Eigen::VectorXd y0(1); y0(0) = 1.0;

    // Plain Euler: order 1, halving h halves the error (ratio ~2).
    auto euler = createIntegrator(IntegratorOptions::EulerExplicit);
    double e100 = std::abs(marchFixed(*euler, 0.0, 4.0, y0, f, 100, opt)(0) - exact);
    double e200 = std::abs(marchFixed(*euler, 0.0, 4.0, y0, f, 200, opt)(0) - exact);
    INFO("Euler plain e100=" << e100 << " e200=" << e200 << " ratio=" << (e100 / e200));
    CHECK_THAT(e100 / e200, WithinRel(2.0, 0.15));

    // Richardson-wrapped Euler: order 2, halving h cuts error ~4x.
    auto rich = wrapRichardson(createIntegrator(IntegratorOptions::EulerExplicit), true);
    double r100 = std::abs(marchFixed(*rich, 0.0, 4.0, y0, f, 100, opt)(0) - exact);
    double r200 = std::abs(marchFixed(*rich, 0.0, 4.0, y0, f, 200, opt)(0) - exact);
    INFO("Euler+Richardson r100=" << r100 << " r200=" << r200 << " ratio=" << (r100 / r200));
    CHECK(r100 < e100);               // strictly better at the same N
    CHECK_THAT(r100 / r200, WithinRel(4.0, 0.30));  // effective order 2
}
