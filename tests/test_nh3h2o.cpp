#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "coolsolve/evaluator.h"
#include "coolsolve/ir.h"
#include "coolsolve/nh3h2o.h"
#include "coolsolve/parser.h"
#include "coolsolve/taylor_jet.h"
#include "AbstractState.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <iomanip>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using coolsolve::nh3h2o::state;
using coolsolve::nh3h2o::stateJets;
using coolsolve::taylor::Jet;

namespace {

/// Relative difference of two doubles (a / b - 1 in magnitude), for the INFO lines.
double relativeDifference(double a, double b) { return std::abs(a - b) / std::abs(b); }

/// The five-point central finite difference of f(T, rho, x) with respect to the variable `var`
/// (0: T, 1: rho, 2: x), accurate to about h^4; tests may use finite differences, the library may not.
double centralDifference(double (*f)(double, double, double), double T, double rho, double x, int var) {
    const double steps[3] = {1e-3 * T, 1e-3 * rho, 1e-3 * std::min(x, 1.0 - x)};
    const double h = steps[var];
    auto at = [&](double k) {
        double v[3] = {T, rho, x};
        v[var] += k * h;
        return f(v[0], v[1], v[2]);
    };
    return (-at(2.0) + 8.0 * at(1.0) - 8.0 * at(-1.0) + at(-2.0)) / (12.0 * h);
}

Jet<3> jetT2(double T, double rho, double x) {
    const Jet<3> jetT = Jet<3>::variable(T, 0);
    const Jet<3> jetRho = Jet<3>::variable(rho, 1);
    const Jet<3> jetX = Jet<3>::variable(x, 2);
    return jetT * jetT * jetRho + jetRho * jetRho * jetX;   // T^2 rho + rho^2 x
}

double pressure(double T, double rho, double x) { return state(T, rho, x).p; }
double enthalpy(double T, double rho, double x) { return state(T, rho, x).h; }
double internalEnergy(double T, double rho, double x) { return state(T, rho, x).u; }
double entropy(double T, double rho, double x) { return state(T, rho, x).s; }

}  // namespace

// ---------------------------------------------------------------------------
// Jet arithmetic against closed forms
// ---------------------------------------------------------------------------

TEST_CASE("Jet arithmetic of a product with exp and log", "[nh3h2o]") {
    // f(T, rho, x) = exp(T rho / 2) ln(1 + x^2), derivatives by hand
    const double T = 1.3, rho = 0.7, x = 0.4;
    const Jet<3> jetT = Jet<3>::variable(T, 0);
    const Jet<3> jetRho = Jet<3>::variable(rho, 1);
    const Jet<3> jetX = Jet<3>::variable(x, 2);
    const Jet<3> f = exp(jetT * jetRho * 0.5) * log(1.0 + jetX * jetX);

    const double a = T * rho / 2.0, E = std::exp(a);
    const double X = 1.0 + x * x;
    const double L = std::log(X), L1 = 2.0 * x / X, L2 = 2.0 * (1.0 - x * x) / (X * X),
                 L3 = -4.0 * x * (3.0 - x * x) / (X * X * X);
    const double ET = rho / 2.0, ER = T / 2.0;  // a_T, a_rho

    CHECK(f.value() == Approx(E * L).margin(1e-14));
    CHECK(f.gradient(0) == Approx(ET * E * L).margin(1e-14));
    CHECK(f.gradient(1) == Approx(ER * E * L).margin(1e-14));
    CHECK(f.gradient(2) == Approx(E * L1).margin(1e-14));
    // Hessian: f_TT, f_T rho, f_T x, f_rho rho, f_rho x, f_xx
    CHECK(f.hessian(0, 0) == Approx(ET * ET * E * L).margin(1e-14));
    CHECK(f.hessian(0, 1) == Approx(E * L * (0.5 + T * rho / 4.0)).margin(1e-14));
    CHECK(f.hessian(0, 2) == Approx(ET * E * L1).margin(1e-14));
    CHECK(f.hessian(1, 1) == Approx(ER * ER * E * L).margin(1e-14));
    CHECK(f.hessian(1, 2) == Approx(ER * E * L1).margin(1e-14));
    CHECK(f.hessian(2, 2) == Approx(E * L2).margin(1e-14));
    CHECK(f.hessian(2, 0) == Approx(f.hessian(0, 2)).margin(0.0));  // symmetric storage
    // Third derivatives, all ten multi-indices
    CHECK(f.thirdDerivative(0, 0, 0) == Approx(ET * ET * ET * E * L).margin(1e-14));
    CHECK(f.thirdDerivative(0, 0, 1) == Approx(E * L * (rho / 2.0 + T * rho * rho / 8.0)).margin(1e-14));
    CHECK(f.thirdDerivative(0, 0, 2) == Approx(ET * ET * E * L1).margin(1e-14));
    CHECK(f.thirdDerivative(0, 1, 1) == Approx(E * L * (T / 2.0 + T * T * rho / 8.0)).margin(1e-14));
    CHECK(f.thirdDerivative(0, 1, 2) == Approx(E * L1 * (0.5 + T * rho / 4.0)).margin(1e-14));
    CHECK(f.thirdDerivative(0, 2, 2) == Approx(ET * E * L2).margin(1e-14));
    CHECK(f.thirdDerivative(1, 1, 1) == Approx(ER * ER * ER * E * L).margin(1e-14));
    CHECK(f.thirdDerivative(1, 1, 2) == Approx(ER * ER * E * L1).margin(1e-14));
    CHECK(f.thirdDerivative(1, 2, 2) == Approx(ER * E * L2).margin(1e-14));
    CHECK(f.thirdDerivative(2, 2, 2) == Approx(E * L3).margin(1e-14));
}

TEST_CASE("Jet arithmetic of a real power and a quotient", "[nh3h2o]") {
    // g(T, rho, x) = A(u) / B(rho) with u = T x + 2, A = u^(5/3), B = 1 + rho; derivatives by hand
    const double T = 1.1, rho = 2.3, x = 0.6, c = 5.0 / 3.0;
    const Jet<3> jetT = Jet<3>::variable(T, 0);
    const Jet<3> jetRho = Jet<3>::variable(rho, 1);
    const Jet<3> jetX = Jet<3>::variable(x, 2);
    const Jet<3> g = pow(jetT * jetX + 2.0, c) / (1.0 + jetRho);

    const double u = T * x + 2.0, v = 1.0 + rho;
    const double A = std::pow(u, c), A1 = c * std::pow(u, c - 1.0), A2 = c * (c - 1.0) * std::pow(u, c - 2.0),
                 A3 = c * (c - 1.0) * (c - 2.0) * std::pow(u, c - 3.0);
    // T-x derivatives of A: dT = x A', dx = T A', dTT = x^2 A'', dTx = A' + x T A'', dxx = T^2 A'', ...
    CHECK(g.value() == Approx(A / v).margin(1e-14));
    CHECK(g.gradient(0) == Approx(x * A1 / v).margin(1e-14));
    CHECK(g.gradient(1) == Approx(-A / (v * v)).margin(1e-14));
    CHECK(g.gradient(2) == Approx(T * A1 / v).margin(1e-14));
    CHECK(g.hessian(0, 0) == Approx(x * x * A2 / v).margin(1e-14));
    CHECK(g.hessian(0, 1) == Approx(-x * A1 / (v * v)).margin(1e-14));
    CHECK(g.hessian(0, 2) == Approx((A1 + x * T * A2) / v).margin(1e-14));
    CHECK(g.hessian(1, 1) == Approx(2.0 * A / (v * v * v)).margin(1e-14));
    CHECK(g.hessian(1, 2) == Approx(-T * A1 / (v * v)).margin(1e-14));
    CHECK(g.hessian(2, 2) == Approx(T * T * A2 / v).margin(1e-14));
    CHECK(g.thirdDerivative(0, 0, 0) == Approx(x * x * x * A3 / v).margin(1e-14));
    CHECK(g.thirdDerivative(0, 0, 1) == Approx(-x * x * A2 / (v * v)).margin(1e-14));
    CHECK(g.thirdDerivative(0, 0, 2) == Approx((2.0 * x * A2 + x * x * T * A3) / v).margin(1e-14));
    CHECK(g.thirdDerivative(0, 1, 1) == Approx(2.0 * x * A1 / (v * v * v)).margin(1e-14));
    CHECK(g.thirdDerivative(0, 1, 2) == Approx(-(A1 + x * T * A2) / (v * v)).margin(1e-14));
    CHECK(g.thirdDerivative(0, 2, 2) == Approx((2.0 * T * A2 + T * T * x * A3) / v).margin(1e-14));
    CHECK(g.thirdDerivative(1, 1, 1) == Approx(-6.0 * A / (v * v * v * v)).margin(1e-14));
    CHECK(g.thirdDerivative(1, 1, 2) == Approx(2.0 * T * A1 / (v * v * v)).margin(1e-14));
    CHECK(g.thirdDerivative(1, 2, 2) == Approx(-T * T * A2 / (v * v)).margin(1e-14));
    CHECK(g.thirdDerivative(2, 2, 2) == Approx(T * T * T * A3 / v).margin(1e-14));
}

TEST_CASE("Jet identities: quotient, integer powers, sqrt, derivative", "[nh3h2o]") {
    const double T = 1.7, rho = 0.9, x = 0.35;
    const Jet<3> a = Jet<3>::variable(T, 0) * Jet<3>::variable(rho, 1) + 3.0 * Jet<3>::variable(x, 2);
    const Jet<3> b = exp(a) + 2.0;

    const Jet<3> quotient = a / b * b;
    for (int i = 0; i < coolsolve::taylor::jetSize(3); ++i)
        CHECK(quotient.c[i] == Approx(a.c[i]).margin(1e-12));

    CHECK(pow(a, 0) == 1.0);
    const Jet<3> a2 = pow(a, 2), a3 = pow(a, 3);
    for (int i = 0; i < coolsolve::taylor::jetSize(3); ++i) {
        CHECK(a2.c[i] == Approx((a * a).c[i]).margin(1e-13));
        CHECK(a3.c[i] == Approx((a * a * a).c[i]).margin(1e-12));
    }
    const Jet<3> squareRoot = sqrt(a);
    for (int i = 0; i < coolsolve::taylor::jetSize(3); ++i)
        CHECK((squareRoot * squareRoot).c[i] == Approx(a.c[i]).margin(1e-12));
    CHECK(-(-a) == a.value());
    CHECK(a - a == 0.0);
    CHECK((2.0 + a) == (a + 2.0));
    CHECK((2.0 * a).value() == Approx(2.0 * a.value()).margin(0.0));
    CHECK((2.0 / b).value() == Approx(2.0 / b.value()).margin(1e-15));
    CHECK((1.0 / b * b).value() == Approx(1.0).margin(1e-13));

    // derivative() of h = T^2 rho + rho^2 x, checked against hand derivatives of order <= 2
    const Jet<3> h = jetT2(T, rho, x);
    const Jet<2> dhdt = derivative(h, 0);
    CHECK(dhdt.value() == Approx(2.0 * T * rho).margin(1e-14));
    CHECK(dhdt.gradient(0) == Approx(2.0 * rho).margin(1e-14));
    CHECK(dhdt.gradient(1) == Approx(2.0 * T).margin(1e-14));
    CHECK(dhdt.gradient(2) == Approx(0.0).margin(1e-14));
    CHECK(dhdt.hessian(0, 1) == Approx(2.0).margin(1e-14));
    const Jet<2> dhdx = derivative(h, 2);
    CHECK(dhdx.value() == Approx(rho * rho).margin(1e-14));
    CHECK(dhdx.gradient(1) == Approx(2.0 * rho).margin(1e-14));
    CHECK(dhdx.hessian(1, 1) == Approx(2.0).margin(1e-14));
    const Jet<1> d2hdt2 = derivative(dhdt, 0);
    CHECK(d2hdt2.value() == Approx(2.0 * rho).margin(1e-14));
}

// ---------------------------------------------------------------------------
// The formulation against the guideline and CoolProp
// ---------------------------------------------------------------------------

TEST_CASE("IAPWS G4-01 Table 6 single-phase verification values", "[nh3h2o][formulation]") {
    struct Row {
        double x, T, rhoDm3, f, pMPa, cv, w;
    };
    const Row rows[] = {
        {0.1, 600, 35, -13734.1763, 32.1221333, 53.3159544, 883.925596},
        {0.1, 600, 4, -16991.6697, 12.7721090, 52.7644553, 471.762394},
        {0.5, 500, 32, -12109.5369, 21.3208159, 58.0077346, 830.295833},
        {0.5, 500, 1, -18281.3020, 3.6423080, 36.8228098, 510.258362},
        {0.9, 400, 30, -6986.4869, 22.2830797, 51.8072415, 895.748711},
        {0.9, 400, 0.5, -13790.6278, 1.5499708, 32.9703870, 478.608147},
    };
    for (const Row& row : rows) {
        CAPTURE(row.x, row.T, row.rhoDm3);
        const coolsolve::nh3h2o::State s = state(row.T, row.rhoDm3 * 1000.0, row.x);
        INFO("relative differences: f " << relativeDifference(s.f, row.f) << ", p "
                                        << relativeDifference(s.p, row.pMPa * 1e6) << ", cv "
                                        << relativeDifference(s.cv, row.cv) << ", w "
                                        << relativeDifference(s.w, row.w));
        CHECK(s.f == Approx(row.f).epsilon(1e-7));
        CHECK(s.p == Approx(row.pMPa * 1e6).epsilon(1e-7));
        CHECK(s.cv == Approx(row.cv).epsilon(1e-7));
        CHECK(s.w == Approx(row.w).epsilon(1e-7));
    }
}

TEST_CASE("The water residual equals IAPWS-95 (CoolProp Water)", "[nh3h2o][formulation]") {
    std::shared_ptr<CoolProp::AbstractState> water(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double states[][2] = {
        {300.0, 55300.0},    // liquid
        {600.0, 4000.0},     // vapour
        {647.096, 17873.7},  // critical
        {800.0, 30000.0},    // supercritical
        {450.0, 200.0},      // dilute
        {365.0, 41000.0},    // compressed liquid near Ttr
    };
    for (const double (&Trho)[2] : states) {
        const double T = Trho[0], rho = Trho[1];
        water->update(CoolProp::DmolarT_INPUTS, rho, T);
        const double tau = water->tau(), delta = water->delta();
        CAPTURE(T, rho, tau, delta, water->alphar());
        const double mine = coolsolve::nh3h2o::residualWater(tau, delta);
        INFO(std::setprecision(17) << "alphar " << mine << " vs " << water->alphar() << " relative " << relativeDifference(mine, water->alphar()));
        CHECK(mine == Approx(water->alphar()).epsilon(1e-9));
        // first derivatives, as an order-1 jet in (tau, delta)
        const Jet<1> jetTau = Jet<1>::variable(tau, 0);
        const Jet<1> jetDelta = Jet<1>::variable(delta, 1);
        const Jet<1> phi = coolsolve::nh3h2o::residualWater(jetTau, jetDelta);
        CHECK(phi.gradient(0) == Approx(water->dalphar_dTau()).epsilon(1e-9));
        CHECK(phi.gradient(1) == Approx(water->dalphar_dDelta()).epsilon(1e-9));
    }
}

TEST_CASE("The water residual at delta = 1 exactly is finite (non-analytic guard)", "[nh3h2o][formulation]") {
    const Jet<2> jetTau = Jet<2>::variable(1.3, 0);
    const Jet<2> jetDelta = Jet<2>::variable(1.0, 1);  // delta == 1 exactly
    const Jet<2> phi = coolsolve::nh3h2o::residualWater(jetTau, jetDelta);
    CHECK(std::isfinite(phi.value()));
    for (int i = 0; i < 3; ++i) {
        CHECK(std::isfinite(phi.gradient(i)));
        for (int j = 0; j < 3; ++j) CHECK(std::isfinite(phi.hessian(i, j)));
    }
    // the tau = delta = 1 corner: the value is 0 for the DELTA^b factors, no NaN anywhere
    const Jet<2> critical = coolsolve::nh3h2o::residualWater(Jet<2>::variable(1.0, 0), Jet<2>::variable(1.0, 1));
    CHECK(std::isfinite(critical.value()));
}

TEST_CASE("The ammonia residual is the Tillner-Roth surface (CoolProp EOS data)", "[nh3h2o][formulation]") {
    // CoolProp's default "Ammonia" is the newer Gao et al. (2023) equation; the surface of the
    // guideline appendix is the second EOS entry (TillnerRoth-DKV-1993), not reachable through
    // the public API. Its generated state data (in dev/fluids/Ammonia.json) is used instead:
    // at x = 1 the module is exactly that equation (same residual coefficients, same reducing
    // constants, same R), so the pressures must agree to machine precision.
    const double referenceStates[][3] = {
        {195.495, 43035.33929207322, 6091.2231081315085},   // saturated liquid at the triple point
        {195.495, 3.763506027681136, 6091.223108650368},    // saturated vapour at the triple point
        {445.94, 11890.599438881152, 18635643.933141705},   // supercritical anchor state
    };
    for (const double (&TrhoP)[3] : referenceStates) {
        const double p = state(TrhoP[0], TrhoP[1], 1.0).p;
        CAPTURE(TrhoP[0], TrhoP[1], p, TrhoP[2]);
        CHECK(p == Approx(TrhoP[2]).epsilon(1e-9));
    }
}

namespace {

/// A single-phase state of `fluid` at (T, p) taken from CoolProp: the molar density, so that the module is
/// compared at a state where CoolProp's pressure is the equation-of-state one (inside the dome it returns the
/// saturation pressure).
double singlePhaseDensity(CoolProp::AbstractState& fluid, double T, double p) {
    fluid.update(CoolProp::PT_INPUTS, p, T);
    return fluid.rhomolar();
}

}  // namespace

TEST_CASE("Pure-fluid limits of p and cv: water against CoolProp to 1/60000", "[nh3h2o][formulation]") {
    // The guideline (Section 3) expects differences up to 1 part in 60000 against IAPWS-95
    // because of the gas constant (here 8.314471, CoolProp's water uses 8.314371357587) and
    // the rounding of the common ideal-gas reducing constants.
    std::shared_ptr<CoolProp::AbstractState> water(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double states[][2] = {{600.0, 2.0e7}, {600.0, 5.0e6}, {300.0, 1.0e5}, {450.0, 1.0e5}, {700.0, 3.0e7}, {800.0, 1.0e6}};
    for (const double (&Tp)[2] : states) {
        const double T = Tp[0], rho = singlePhaseDensity(*water, T, Tp[1]);
        const coolsolve::nh3h2o::State s = state(T, rho, 0.0);
        CAPTURE(T, rho);
        INFO("p: " << s.p << " vs " << water->p() << " (relative " << relativeDifference(s.p, water->p())
                   << "), cv: " << s.cv << " vs " << water->cvmolar() << " (relative "
                   << relativeDifference(s.cv, water->cvmolar()) << ")");
        CHECK(s.p == Approx(water->p()).epsilon(1.0 / 60000.0));
        CHECK(s.cv == Approx(water->cvmolar()).epsilon(1.0e-4));   // the refit ideal-gas part moves cv a little more
    }
}

TEST_CASE("Pure ammonia p and cv against CoolProp's current Ammonia EOS", "[nh3h2o][formulation]") {
    // A cross-formulation sanity check only: CoolProp's default Ammonia is the Gao et al.
    // (2023) surface, not the guideline appendix surface, so the agreement is that between the
    // two published ammonia equations (a fraction of a percent in the vapour), not the 1/60000 of the guideline.
    std::shared_ptr<CoolProp::AbstractState> ammonia(CoolProp::AbstractState::factory("HEOS", "Ammonia"));
    const double states[][2] = {{350.0, 1.0e5}, {400.0, 1.0e6}, {450.0, 5.0e6}, {500.0, 5.0e5}};   // vapour: no pressure-density amplification
    for (const double (&Tp)[2] : states) {
        const double T = Tp[0], rho = singlePhaseDensity(*ammonia, T, Tp[1]);
        const coolsolve::nh3h2o::State s = state(T, rho, 1.0);
        CAPTURE(T, rho);
        INFO("p: " << s.p << " vs " << ammonia->p() << " (relative " << relativeDifference(s.p, ammonia->p())
                   << "), cv: " << s.cv << " vs " << ammonia->cvmolar() << " (relative "
                   << relativeDifference(s.cv, ammonia->cvmolar()) << ")");
        CHECK(s.p == Approx(ammonia->p()).epsilon(0.01));
        CHECK(s.cv == Approx(ammonia->cvmolar()).epsilon(0.03));
    }
}

TEST_CASE("Jet-derived property derivatives are exact (finite-difference check)", "[nh3h2o][formulation]") {
    const double states[][3] = {{600, 35000, 0.1}, {500, 1000, 0.5}, {400, 30000, 0.9}, {450, 8000, 0.3}};
    for (const double (&TrhoX)[3] : states) {
        const double T = TrhoX[0], rho = TrhoX[1], x = TrhoX[2];
        const coolsolve::nh3h2o::StateJets jets = stateJets(T, rho, x);
        CAPTURE(T, rho, x);
        for (int var = 0; var < 3; ++var) {
            INFO("variable " << var << ": dp " << relativeDifference(jets.p.gradient(var), centralDifference(pressure, T, rho, x, var))
                             << ", dh " << relativeDifference(jets.h.gradient(var), centralDifference(enthalpy, T, rho, x, var))
                             << ", du " << relativeDifference(jets.u.gradient(var), centralDifference(internalEnergy, T, rho, x, var))
                             << ", ds " << relativeDifference(jets.s.gradient(var), centralDifference(entropy, T, rho, x, var)));
            CHECK(jets.p.gradient(var) == Approx(centralDifference(pressure, T, rho, x, var)).epsilon(1e-6));
            CHECK(jets.h.gradient(var) == Approx(centralDifference(enthalpy, T, rho, x, var)).epsilon(1e-6));
            CHECK(jets.u.gradient(var) == Approx(centralDifference(internalEnergy, T, rho, x, var)).epsilon(1e-6));
            CHECK(jets.s.gradient(var) == Approx(centralDifference(entropy, T, rho, x, var)).epsilon(1e-6));
        }
        // one Hessian spot check: d2p/drho2 by a second central difference (looser, round-off)
        const double h = 1e-3 * rho;
        const double d2pdrho2 =
            (pressure(T, rho + h, x) - 2.0 * pressure(T, rho, x) + pressure(T, rho - h, x)) / (h * h);
        CHECK(jets.p.hessian(1, 1) == Approx(d2pdrho2).epsilon(1e-5));
    }
}

// ---------------------------------------------------------------------------
// Domain errors
// ---------------------------------------------------------------------------

TEST_CASE("Domain errors of the ammonia-water core", "[nh3h2o]") {
    REQUIRE_THROWS_AS(state(600.0, 1000.0, -0.1), std::runtime_error);
    REQUIRE_THROWS_WITH(state(600.0, 1000.0, -0.1), ContainsSubstring("mole fraction"));
    REQUIRE_THROWS_WITH(state(600.0, 1000.0, 1.1), ContainsSubstring("mole fraction"));
    REQUIRE_THROWS_WITH(state(0.0, 1000.0, 0.5), ContainsSubstring("temperature"));
    REQUIRE_THROWS_WITH(state(-5.0, 1000.0, 0.5), ContainsSubstring("temperature"));
    REQUIRE_THROWS_WITH(state(600.0, 0.0, 0.5), ContainsSubstring("density"));
    REQUIRE_THROWS_WITH(state(600.0, -100.0, 0.5), ContainsSubstring("density"));
    REQUIRE_THROWS_WITH(state(std::nan(""), 1000.0, 0.5), ContainsSubstring("not finite"));
    REQUIRE_THROWS_AS(stateJets(600.0, 1000.0, 2.0), std::runtime_error);
    REQUIRE_THROWS_AS(coolsolve::nh3h2o::reducedHelmholtz(600.0, 1000.0, -1.0), std::runtime_error);
    // the pure limits themselves are part of the domain
    CHECK(std::isfinite(state(600.0, 1000.0, 0.0).p));
    CHECK(std::isfinite(state(400.0, 10000.0, 1.0).p));
}

// ---------------------------------------------------------------------------
// Phase equilibrium: the verification tables of the guideline (Tables 7 and 8)
// ---------------------------------------------------------------------------

namespace {
using coolsolve::nh3h2o::Saturation;
using coolsolve::nh3h2o::SatSpec;
using coolsolve::nh3h2o::saturation;
}  // namespace

TEST_CASE("IAPWS G4-01 Table 7: bubble points", "[nh3h2o][equilibrium]") {
    struct Row {
        double xL, T, pMPa, xV, rhoL, rhoV;   // densities in mol/dm3
    };
    const Row rows[] = {{0.2, 300, 0.040710, 0.9360, 51.941, 0.01640},
                        {0.4, 400, 2.5545, 0.9363, 43.318, 0.8608},
                        {0.6, 500, 16.698, 0.7844, 25.459, 8.86}};
    for (const Row& row : rows) {
        const Saturation s = saturation(SatSpec::BubbleT, row.T, row.xL);
        CAPTURE(row.xL, row.T);
        INFO("p " << s.p << " (" << row.pMPa * 1e6 << "), xV " << s.xV << " (" << row.xV << "), rhoL " << s.rhoL / 1e3
                  << " (" << row.rhoL << "), rhoV " << s.rhoV / 1e3 << " (" << row.rhoV << ")");
        CHECK(s.p == Approx(row.pMPa * 1e6).epsilon(2e-4));
        CHECK(s.xV == Approx(row.xV).margin(6e-5));
        CHECK(s.rhoL / 1e3 == Approx(row.rhoL).epsilon(2e-4));
        CHECK(s.rhoV / 1e3 == Approx(row.rhoV).epsilon(6e-4));
        if (s.p < 5e6) {   // the bubble point at the pressure gives back the temperature (near the critical locus it need not)
            const Saturation at = saturation(SatSpec::BubbleP, s.p, row.xL);
            CHECK(at.T == Approx(row.T).epsilon(1e-8));
        }
    }
}

TEST_CASE("IAPWS G4-01 Table 8: dew points", "[nh3h2o][equilibrium]") {
    struct Row {
        double xV, T, pMPa, xL, rhoL, rhoV;
    };
    const Row rows[] = {{0.2, 300, 0.00437062, 0.010672, 55.16434, 0.00175506},
                        {0.4, 400, 0.394694, 0.051541, 50.83187, 0.122658},
                        {0.6, 500, 6.52607, 0.22135, 39.93714, 2.00730}};
    for (const Row& row : rows) {
        const Saturation s = saturation(SatSpec::DewT, row.T, row.xV);
        CAPTURE(row.xV, row.T);
        INFO("p " << s.p << " (" << row.pMPa * 1e6 << "), xL " << s.xL << " (" << row.xL << "), rhoL " << s.rhoL / 1e3
                  << " (" << row.rhoL << "), rhoV " << s.rhoV / 1e3 << " (" << row.rhoV << ")");
        CHECK(s.p == Approx(row.pMPa * 1e6).epsilon(2e-4));
        CHECK(s.xL == Approx(row.xL).epsilon(2e-4));
        CHECK(s.rhoL / 1e3 == Approx(row.rhoL).epsilon(2e-4));
        CHECK(s.rhoV / 1e3 == Approx(row.rhoV).epsilon(2e-4));
        const Saturation at = saturation(SatSpec::DewP, s.p, row.xV);
        CHECK(at.T == Approx(row.T).epsilon(1e-8));
    }
}

// ---------------------------------------------------------------------------
// The flash: the three-property states of the EES procedure NH3H2O
// ---------------------------------------------------------------------------

namespace {
using coolsolve::nh3h2o::flash;
using coolsolve::nh3h2o::FlashResult;
using coolsolve::nh3h2o::Prop;
using coolsolve::nh3h2o::singlePhaseExtras;

const char* const kPropNames[8] = {"T", "P", "x", "h", "s", "u", "v", "q"};

/// The 15 combinations of the EES procedure, as property numbers.
const int kCombinations[15][3] = {{1, 2, 3}, {1, 2, 8}, {1, 3, 7}, {1, 3, 8}, {1, 4, 8}, {1, 5, 8}, {1, 6, 8}, {1, 7, 8},
                                  {2, 3, 4}, {2, 3, 5}, {2, 3, 8}, {2, 4, 8}, {2, 5, 8}, {2, 6, 8}, {2, 7, 8}};

FlashResult flashOf(const int combination[3], const FlashResult& state) {
    const double in[3] = {state.value[combination[0] - 1], state.value[combination[1] - 1], state.value[combination[2] - 1]};
    return flash(static_cast<Prop>(combination[0]), static_cast<Prop>(combination[1]), static_cast<Prop>(combination[2]),
                 in[0], in[1], in[2]);
}
}  // namespace

TEST_CASE("NH3H2O flash: the EES reference points", "[nh3h2o][flash]") {
    // EES (Ibrahim & Klein 1993) against this implementation (IAPWS G4-01): the two formulations agree within
    // their own uncertainty, the reference states of h, s and u differ (so differences at fixed x are compared)
    const FlashResult vapour = flash(Prop::T, Prop::P, Prop::X, 450.0, 1.0e6, 0.5);   // EES: v 0.2036, q 1.001
    CHECK(vapour.quality() == 1.001);
    CHECK(vapour.volume() == Approx(0.2036).epsilon(2e-3));
    const FlashResult liquid = flash(Prop::T, Prop::P, Prop::X, 283.15, 1.0e6, 0.5);   // EES: h -194.8 kJ/kg, q -0.001
    CHECK(liquid.quality() == -0.001);
    CHECK_FALSE(liquid.twoPhase);
    CHECK((vapour.enthalpy() - liquid.enthalpy()) / 1e3 == Approx(2224.0 + 194.8).epsilon(1e-2));
    const FlashResult bubble = flash(Prop::T, Prop::P, Prop::Q, 323.15, 1.0e6, 0.0);   // EES: x = 0.5953
    CHECK(bubble.massFraction() == Approx(0.5953).epsilon(2e-2));
    CHECK(bubble.twoPhase);
}

TEST_CASE("NH3H2O flash: every combination returns the state it was built from", "[nh3h2o][flash]") {
    const double references[][3] = {{350.0, 8e5, 0.3}, {300.0, 2e5, 0.7}, {400.0, 2e6, 0.5}, {330.0, 1e6, 0.5}};
    for (const double (&ref)[3] : references) {
        const FlashResult state = flash(Prop::T, Prop::P, Prop::Q, ref[0], ref[1], ref[2]);
        REQUIRE(state.twoPhase);
        for (const int(&combination)[3] : kCombinations) {
            const std::string name = std::string(kPropNames[combination[0] - 1]) + kPropNames[combination[1] - 1] +
                                     kPropNames[combination[2] - 1];
            CAPTURE(ref[0], ref[1], ref[2], name);
            const FlashResult again = flashOf(combination, state);
            CHECK(again.twoPhase);
            // (T, h|u, q) has two solutions on the bubble curve (heat of mixing); the one returned must satisfy the inputs
            const bool ambiguous = combination[0] == 1 && (combination[1] == 4 || combination[1] == 6) && combination[2] == 8;
            for (int k = 0; k < 8; ++k) {
                const bool isInput = k + 1 == combination[0] || k + 1 == combination[1] || k + 1 == combination[2];
                if (ambiguous && !isInput) continue;
                const double scale = std::max(std::abs(state.value[k]), k == 7 ? 1.0 : 1e-3);
                INFO("output " << kPropNames[k] << ": " << again.value[k] << " vs " << state.value[k]);
                CHECK(std::abs(again.value[k] - state.value[k]) / scale < 1e-7);
            }
        }
    }
    // single-phase states: the four combinations without q
    const double singles[][3] = {{300.0, 5e5, 0.4}, {420.0, 5e5, 0.9}, {300.0, 1e7, 0.2}, {600.0, 2e6, 0.5}};
    for (const double (&ref)[3] : singles) {
        const FlashResult state = flash(Prop::T, Prop::P, Prop::X, ref[0], ref[1], ref[2]);
        REQUIRE_FALSE(state.twoPhase);
        for (const int(&combination)[3] : kCombinations) {
            if (combination[2] == 8) continue;
            CAPTURE(ref[0], ref[1], ref[2], combination[0], combination[1], combination[2]);
            const FlashResult again = flashOf(combination, state);
            CHECK_FALSE(again.twoPhase);
            for (int k = 0; k < 8; ++k)
                CHECK(std::abs(again.value[k] - state.value[k]) / std::max(std::abs(state.value[k]), 1e-3) < 1e-7);
        }
    }
}

TEST_CASE("NH3H2O flash: exact derivatives against finite differences", "[nh3h2o][flash]") {
    struct Case {
        int combination[3];
        double inputs[3];
    };
    const Case cases[] = {
        {{1, 2, 3}, {300.0, 5e5, 0.4}},      // sub-cooled liquid
        {{1, 2, 3}, {420.0, 5e5, 0.9}},      // superheated vapour
        {{1, 2, 3}, {350.0, 8e5, 0.5533}},   // two-phase
        {{1, 2, 8}, {350.0, 8e5, 0.3}},
        {{1, 3, 7}, {350.0, 0.5533, 0.06167}},
        {{1, 3, 8}, {330.0, 0.5, 0.0}},
        {{1, 3, 8}, {330.0, 0.5, 0.4}},
        {{1, 4, 8}, {300.0, 1189852.7, 0.7}},
        {{1, 5, 8}, {300.0, 4969.5, 0.7}},
        {{2, 3, 4}, {1e6, 0.4, 1.5e5}},      // sub-cooled liquid
        {{2, 3, 4}, {8e5, 0.5533, 7.0e5}},   // two-phase
        {{2, 3, 5}, {8e5, 0.5533, 2950.0}},
        {{2, 3, 8}, {1e6, 0.5, 0.0}},
        {{2, 3, 8}, {1e6, 0.98, 1.0}},
        {{2, 3, 8}, {8e5, 0.5533, 0.3}},
        {{2, 4, 8}, {8e5, 707713.7, 0.3}},
        {{2, 7, 8}, {8e5, 0.06, 0.3}},
    };
    for (const Case& c : cases) {
        const auto call = [&](const double* in) {
            return flash(static_cast<Prop>(c.combination[0]), static_cast<Prop>(c.combination[1]),
                         static_cast<Prop>(c.combination[2]), in[0], in[1], in[2]);
        };
        FlashResult base;
        try {
            base = call(c.inputs);
        } catch (const std::runtime_error& e) {
            FAIL("flash failed: " << e.what());
        }
        for (int j = 0; j < 3; ++j) {
            const double h = 1e-4 * std::max(std::abs(c.inputs[j]), 1e-2);
            double in[3] = {c.inputs[0], c.inputs[1], c.inputs[2]};
            FlashResult at[4];
            const double k[4] = {2.0, 1.0, -1.0, -2.0};
            bool usable = true;
            for (int m = 0; m < 4; ++m) {
                in[j] = c.inputs[j] + k[m] * h;
                if (c.combination[j] == 8 && (in[j] < 0.0 || in[j] > 1.0)) usable = false;   // q is at its bound
                if (!usable) break;
                try {
                    at[m] = call(in);
                } catch (const std::runtime_error&) {
                    usable = false;
                }
            }
            if (!usable) continue;
            for (int out = 0; out < 7; ++out) {
                if (out + 1 == c.combination[0] || out + 1 == c.combination[1] || out + 1 == c.combination[2]) continue;   // an input: 0 or 1
                const double fd = (-at[0].value[out] + 8.0 * at[1].value[out] - 8.0 * at[2].value[out] + at[3].value[out]) / (12.0 * h);
                const double scale = std::max({std::abs(fd), std::abs(base.derivative[out][j]), 1e-12});
                INFO("combination " << c.combination[0] << c.combination[1] << c.combination[2] << ", d" << kPropNames[out]
                                    << "/d" << kPropNames[c.combination[j] - 1] << ": exact " << base.derivative[out][j]
                                    << ", finite difference " << fd);
                CHECK(std::abs(base.derivative[out][j] - fd) <= 1e-4 * scale + 1e-9 * std::abs(base.value[out]) / std::max(std::abs(c.inputs[j]), 1e-9));
            }
        }
    }
}

TEST_CASE("NH3H2O flash: specific heat and expansion coefficient of single-phase states", "[nh3h2o][flash]") {
    const double T = 300.0, p = 5e5, x = 0.4;
    const auto extras = singlePhaseExtras(Prop::T, Prop::P, Prop::X, T, p, x);
    // cp = (dh/dT)_p and the expansion coefficient (1/v)(dv/dT)_p by finite differences of the flash
    const double h = 1e-3;
    auto at = [&](double TT) { return flash(Prop::T, Prop::P, Prop::X, TT, p, x); };
    const double cpFd = (-at(T + 2 * h).enthalpy() + 8 * at(T + h).enthalpy() - 8 * at(T - h).enthalpy() + at(T - 2 * h).enthalpy()) / (12 * h);
    const double betaFd = (-at(T + 2 * h).volume() + 8 * at(T + h).volume() - 8 * at(T - h).volume() + at(T - 2 * h).volume()) / (12 * h) / at(T).volume();
    CHECK(extras.cp == Approx(cpFd).epsilon(1e-6));
    CHECK(extras.expansion == Approx(betaFd).epsilon(1e-6));
    CHECK(extras.cp > 3000.0);
    // derivatives with respect to the inputs: d cp / d p by finite difference
    auto cpAt = [&](double pp) { return singlePhaseExtras(Prop::T, Prop::P, Prop::X, T, pp, x).cp; };
    const double dp = 10.0;
    const double dcpdp = (-cpAt(p + 2 * dp) + 8 * cpAt(p + dp) - 8 * cpAt(p - dp) + cpAt(p - 2 * dp)) / (12 * dp);
    CHECK(extras.dcp[1] == Approx(dcpdp).epsilon(1e-4));
    REQUIRE_THROWS_WITH(singlePhaseExtras(Prop::T, Prop::P, Prop::Q, 350.0, 8e5, 0.3), ContainsSubstring("single-phase"));
}

TEST_CASE("NH3H2O flash: errors name the state", "[nh3h2o][flash]") {
    REQUIRE_THROWS_WITH(flash(Prop::P, Prop::T, Prop::X, 1e6, 300.0, 0.5), ContainsSubstring("order"));
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::H, Prop::S, 300.0, 1e5, 1e3), ContainsSubstring("unsupported combination"));
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::P, Prop::X, 300.0, 5e5, 1.2), ContainsSubstring("mass fraction"));
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::P, Prop::X, 300.0, -5e5, 0.4), ContainsSubstring("pressure"));
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::P, Prop::X, -300.0, 5e5, 0.4), ContainsSubstring("temperature"));
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::X, Prop::Q, 300.0, 0.4, 1.5), ContainsSubstring("quality"));
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::X, Prop::V, 300.0, 0.4, -1.0), ContainsSubstring("volume"));
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::X, Prop::Q, std::nan(""), 0.4, 0.5), ContainsSubstring("not finite"));
    // a quality input needs a temperature or pressure at which the two phases exist
    REQUIRE_THROWS_WITH(flash(Prop::T, Prop::X, Prop::Q, 700.0, 0.4, 0.5), ContainsSubstring("NH3H2O(T = 700"));
    // an enthalpy below the sub-cooled range
    REQUIRE_THROWS_AS(flash(Prop::P, Prop::X, Prop::H, 1e6, 0.4, -1e8), std::runtime_error);
}

// ---------------------------------------------------------------------------
// The solver is strict: the equilibrium equations hold at what it returns
// ---------------------------------------------------------------------------

namespace {
/// p, mu_1 and mu_2 of a phase at (T, rho, x), recomputed here from the Helmholtz energy by second-order jets, with no
/// code of the solver: p = rho^2 f_rho, mu_i = f + p/rho + (delta_i2 - x) f_x.
struct IndependentPhase {
    double p, mu1, mu2;
};

IndependentPhase independentPhase(double T, double rho, double x) {
    const Jet<2> jT = Jet<2>::variable(T, 0), jR = Jet<2>::variable(rho, 1), jX = Jet<2>::variable(x, 2);
    const Jet<2> f = coolsolve::nh3h2o::kR * jT * coolsolve::nh3h2o::reducedHelmholtz(jT, jR, jX);
    const double fValue = f.value(), fRho = f.gradient(1), fX = f.gradient(2);
    const double p = rho * rho * fRho, g = fValue + p / rho;
    return {p, g - x * fX, g + (1.0 - x) * fX};
}
}  // namespace

TEST_CASE("Saturation: the equilibrium equations hold to 1e-11 of their natural scales", "[nh3h2o][equilibrium]") {
    // pressures are compared on the scale rho R T of a liquid (its pressure is a sum of terms of that size: round-off
    // 3e-13 of it, measured), the chemical potentials on R * 300 K
    using coolsolve::nh3h2o::kR;
    int checked = 0;
    for (double T : {280.0, 320.0, 360.0, 400.0})
        for (double x : {0.05, 0.2, 0.5, 0.8, 0.95}) {
            for (SatSpec spec : {SatSpec::BubbleT, SatSpec::DewT}) {
                Saturation s;
                try {
                    s = saturation(spec, T, x);
                } catch (const std::runtime_error&) {
                    continue;   // above the critical locus: no two-phase state
                }
                const IndependentPhase L = independentPhase(s.T, s.rhoL, s.xL), V = independentPhase(s.T, s.rhoV, s.xV);
                CAPTURE(T, x, static_cast<int>(spec), s.p, s.xL, s.xV);
                CHECK(std::abs(L.p - s.p) <= 1e-11 * s.rhoL * kR * s.T);
                CHECK(std::abs(V.p - s.p) <= 1e-11 * std::max(s.p, s.rhoV * kR * s.T));
                CHECK(std::abs(L.mu1 - V.mu1) <= 1e-11 * kR * 300.0);
                CHECK(std::abs(L.mu2 - V.mu2) <= 1e-11 * kR * 300.0);
                ++checked;
            }
        }
    CHECK(checked >= 30);
}

// ---------------------------------------------------------------------------
// The language: the fluid NH3H2O and CALL NH3H2O in the evaluator
// ---------------------------------------------------------------------------

namespace {

struct AwCall {
    double value;
    std::vector<double> gradient;
};

/// Value and gradient of the property call `call` written with the variables `names` at `values`, straight from the
/// automatic differentiation of the evaluator.
AwCall evalAw(const std::string& call, const std::vector<std::string>& names, const std::vector<double>& values,
              const coolsolve::CoolPropConfig& cfg = coolsolve::CoolPropConfig()) {
    coolsolve::EESParser parser;
    auto parsed = parser.parse("y = " + call + "\n");
    REQUIRE(parsed.success);
    auto ir = coolsolve::IR::fromAST(parsed.program);
    coolsolve::ExpressionEvaluator ev(names.size(), cfg);
    for (size_t i = 0; i < names.size(); ++i)
        ev.setVariable(names[i], coolsolve::ADValue::independent(values[i], i, names.size()));
    const coolsolve::ADValue r = ev.evaluate(ir.getEquations()[0].rhs);
    return {r.value, r.gradient};
}

std::string awError(const std::string& call) {
    try {
        evalAw(call, {}, {});
    } catch (const std::exception& e) {
        return e.what();
    }
    return "";
}

}  // namespace

TEST_CASE("NH3H2O fluid functions: the values of the module in the default units (C, Pa, J/kg)", "[nh3h2o][evaluator]") {
    const FlashResult ref = flash(Prop::T, Prop::P, Prop::X, 283.15, 1.0e6, 0.5);
    const std::vector<std::string> names = {"Tv", "Pv", "xv"};
    const std::vector<double> values = {10.0, 1.0e6, 0.5};
    CHECK(evalAw("Enthalpy(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(ref.enthalpy()).epsilon(1e-12));
    CHECK(evalAw("Entropy(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(ref.entropy()).epsilon(1e-12));
    CHECK(evalAw("IntEnergy(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(ref.internalEnergy()).epsilon(1e-12));
    CHECK(evalAw("Volume(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(ref.volume()).epsilon(1e-12));
    CHECK(evalAw("Density(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(1.0 / ref.volume()).epsilon(1e-12));
    CHECK(evalAw("Quality(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == -0.001);
    CHECK(evalAw("Temperature(NH3H2O, P=Pv, X=xv, H=hv)", {"Pv", "xv", "hv"}, {1.0e6, 0.5, ref.enthalpy()}).value ==
          Approx(10.0).margin(1e-9));
    const auto extras = singlePhaseExtras(Prop::T, Prop::P, Prop::X, 283.15, 1.0e6, 0.5);
    CHECK(evalAw("SpecHeat(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(extras.cp).epsilon(1e-12));
    CHECK(evalAw("Cp(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(extras.cp).epsilon(1e-12));
    CHECK(evalAw("VolExpCoef(NH3H2O, T=Tv, P=Pv, X=xv)", names, values).value == Approx(extras.expansion).epsilon(1e-12));
    // the fluid name is read in any case and as a string
    CHECK(evalAw("Enthalpy('nh3h2o', T=Tv, P=Pv, X=xv)", names, values).value == Approx(ref.enthalpy()).epsilon(1e-12));
    // a two-phase state: the mass fraction of the bubble point at 50 C and 10 bar
    const FlashResult bubble = flash(Prop::T, Prop::P, Prop::Q, 323.15, 1.0e6, 0.0);
    CHECK(evalAw("MassFraction(NH3H2O, T=Tv, P=Pv, Q=qv)", {"Tv", "Pv", "qv"}, {50.0, 1.0e6, 0.0}).value ==
          Approx(bubble.massFraction()).epsilon(1e-12));
}

TEST_CASE("NH3H2O fluid functions: another unit system (K, bar, kJ/kg)", "[nh3h2o][evaluator]") {
    coolsolve::CoolPropConfig cfg;
    cfg.units.temperature = "K";
    cfg.units.pressure = "bar";
    cfg.units.specific_energy = "kJ/kg";
    cfg.units.specific_entropy = "kJ/kg-K";
    cfg.units.specific_heat = "kJ/kg-K";
    const FlashResult ref = flash(Prop::T, Prop::P, Prop::X, 450.0, 1.0e6, 0.5);
    const std::vector<std::string> names = {"Tv", "Pv", "xv"};
    const std::vector<double> values = {450.0, 10.0, 0.5};
    CHECK(evalAw("Enthalpy(NH3H2O, T=Tv, P=Pv, X=xv)", names, values, cfg).value == Approx(ref.enthalpy() / 1e3).epsilon(1e-12));
    CHECK(evalAw("Entropy(NH3H2O, T=Tv, P=Pv, X=xv)", names, values, cfg).value == Approx(ref.entropy() / 1e3).epsilon(1e-12));
    CHECK(evalAw("Pressure(NH3H2O, T=Tv, X=xv, V=vv)", {"Tv", "xv", "vv"}, {450.0, 0.5, ref.volume()}, cfg).value ==
          Approx(10.0).margin(1e-8));
    // the gradient follows the units: d h[kJ/kg] / d P[bar] = dh/dP * 1e5 / 1e3
    const FlashResult plus = flash(Prop::T, Prop::P, Prop::X, 450.0, 1.0e6 + 100.0, 0.5);
    const FlashResult minus = flash(Prop::T, Prop::P, Prop::X, 450.0, 1.0e6 - 100.0, 0.5);
    const double dhdp = (plus.enthalpy() - minus.enthalpy()) / 200.0;
    const AwCall h = evalAw("Enthalpy(NH3H2O, T=Tv, P=Pv, X=xv)", names, values, cfg);
    CHECK(h.gradient[1] == Approx(dhdp * 1e5 / 1e3).epsilon(1e-6));
    // specific heat in kJ/kg-K
    const auto extras = singlePhaseExtras(Prop::T, Prop::P, Prop::X, 450.0, 1.0e6, 0.5);
    CHECK(evalAw("SpecHeat(NH3H2O, T=Tv, P=Pv, X=xv)", names, values, cfg).value == Approx(extras.cp / 1e3).epsilon(1e-12));
}

TEST_CASE("NH3H2O fluid functions: exact derivatives through the evaluator (finite differences of the call)", "[nh3h2o][evaluator]") {
    struct Case {
        std::string call;
        std::vector<std::string> names;
        std::vector<double> values;
    };
    const std::vector<Case> cases = {
        {"Enthalpy(NH3H2O, T=Tv, P=Pv, X=xv)", {"Tv", "Pv", "xv"}, {30.0, 5e5, 0.4}},                 // liquid
        {"Entropy(NH3H2O, T=Tv, P=Pv, X=xv)", {"Tv", "Pv", "xv"}, {150.0, 5e5, 0.9}},                 // vapour
        {"Volume(NH3H2O, T=Tv, P=Pv, X=xv)", {"Tv", "Pv", "xv"}, {77.0, 8e5, 0.5533}},                // two-phase
        {"Density(NH3H2O, T=Tv, P=Pv, X=xv)", {"Tv", "Pv", "xv"}, {30.0, 5e5, 0.4}},
        {"Temperature(NH3H2O, P=Pv, X=xv, Q=qv)", {"Pv", "xv", "qv"}, {8e5, 0.5533, 0.3}},
        {"Pressure(NH3H2O, T=Tv, X=xv, Q=qv)", {"Tv", "xv", "qv"}, {60.0, 0.5, 0.2}},
        {"MassFraction(NH3H2O, T=Tv, P=Pv, Q=qv)", {"Tv", "Pv", "qv"}, {60.0, 8e5, 0.4}},
        {"Temperature(NH3H2O, P=Pv, X=xv, H=hv)", {"Pv", "xv", "hv"}, {1e6, 0.4, 1.5e5}},
        {"SpecHeat(NH3H2O, T=Tv, P=Pv, X=xv)", {"Tv", "Pv", "xv"}, {30.0, 5e5, 0.4}},
        {"VolExpCoef(NH3H2O, T=Tv, P=Pv, X=xv)", {"Tv", "Pv", "xv"}, {30.0, 5e5, 0.4}},
    };
    for (const Case& c : cases) {
        const AwCall base = evalAw(c.call, c.names, c.values);
        for (size_t j = 0; j < c.values.size(); ++j) {
            const double h = 1e-4 * std::max(std::abs(c.values[j]), 1e-2);
            auto at = [&](double k) {
                std::vector<double> v = c.values;
                v[j] += k * h;
                return evalAw(c.call, c.names, v).value;
            };
            const double fd = (-at(2.0) + 8.0 * at(1.0) - 8.0 * at(-1.0) + at(-2.0)) / (12.0 * h);
            INFO(c.call << " d/d" << c.names[j] << ": exact " << base.gradient[j] << ", finite difference " << fd);
            CHECK(std::abs(base.gradient[j] - fd) <= 1e-4 * std::max({std::abs(fd), std::abs(base.gradient[j]), 1e-9}));
        }
    }
}

// ---------------------------------------------------------------------------
// Near-pure ammonia: the refrigerant of an absorption chiller (0.9995 mass fraction)
// ---------------------------------------------------------------------------

TEST_CASE("Near-pure refrigerant: evaporation of 0.9995 ammonia has a large temperature glide and every state solves", "[nh3h2o][flash][near-pure]") {
    // bubble point at -10 C, then the states of increasing quality at that pressure: the temperature rises by ~19 K
    // between the start of the evaporation and the dew point, because the water accumulates in the liquid
    const double x = 0.9995;
    const FlashResult bubble = flash(Prop::T, Prop::X, Prop::Q, 263.15, x, 0.0);
    const double p = bubble.pressure();
    CHECK(p == Approx(2.9057e5).epsilon(2e-3));
    double previousT = bubble.temperature(), previousH = bubble.enthalpy();
    for (double q : {0.1, 0.5, 0.9, 0.95, 0.98, 0.99, 0.999, 1.0}) {
        const FlashResult r = flash(Prop::P, Prop::X, Prop::Q, p, x, q);
        CAPTURE(q);
        REQUIRE(r.twoPhase);
        CHECK(r.temperature() > previousT);
        CHECK(r.enthalpy() > previousH);
        previousT = r.temperature();
        previousH = r.enthalpy();
        // the quality given back by the state: (P, X, H) -> Q
        if (q < 1.0)   // exactly saturated vapour is on the boundary: it is labelled single-phase (1.001), as documented
            CHECK(flash(Prop::P, Prop::X, Prop::H, p, x, r.enthalpy()).quality() == Approx(q).margin(1e-8));
    }
    CHECK(previousT - 273.15 == Approx(8.9).margin(0.3));   // the dew point: +8.9 C
    CHECK(flash(Prop::P, Prop::X, Prop::Q, p, x, 0.99).temperature() - 273.15 == Approx(-8.8).margin(0.1));
}

TEST_CASE("Near-pure compositions: bubble and dew points from the pure-ammonia side", "[nh3h2o][equilibrium][near-pure]") {
    // mole fractions within 1e-2 of pure: the vapour water fraction is below 2e-4 and `1 - x` of a double is known to 1e-10 only
    // there, so the equilibrium residual is checked against that resolution, not against 1e-11
    using coolsolve::nh3h2o::kR;
    for (double T : {263.15, 300.0, 340.0})
        for (double xi : {0.99, 0.998, 0.9995}) {
            const double x = coolsolve::nh3h2o::moleFractionOf(xi);
            Saturation s;
            try {
                s = saturation(SatSpec::BubbleT, T, x);
            } catch (const std::runtime_error& e) {
                FAIL("bubble point failed at T = " << T << ", x = " << xi << ": " << e.what());
            }
            const IndependentPhase L = independentPhase(s.T, s.rhoL, s.xL), V = independentPhase(s.T, s.rhoV, s.xV);
            const double wV = 1.0 - s.xV;   // water mole fraction of the vapour
            CAPTURE(T, xi, s.p, s.xV);
            CHECK(wV < 2e-4);
            CHECK(std::abs(L.p - s.p) <= 1e-9 * s.rhoL * kR * s.T);
            CHECK(std::abs(V.p - s.p) <= 1e-9 * s.p);
            // mu_1 contains ln(1 - x): its resolution is kR T ulp / (1 - x) of the vapour
            const double resolution = kR * s.T * 1.2e-16 / wV;
            CHECK(std::abs(L.mu1 - V.mu1) <= 1e-9 * kR * 300.0 + 4.0 * resolution);
            CHECK(std::abs(L.mu2 - V.mu2) <= 1e-9 * kR * 300.0);
            // the dew point of the vapour of that bubble point gives back its pressure
            const Saturation d = saturation(SatSpec::DewT, T, s.xV);
            CHECK(d.p == Approx(s.p).epsilon(1e-6));
            CHECK(d.xL == Approx(s.xL).epsilon(1e-6));
        }
}

// ---------------------------------------------------------------------------
// Correctness of the saturation and of the flashes: stable roots only, no single phase from a failed saturation, the
// same state whatever the earlier calls of the thread (the warm-start cache of the saturation solver is thread_local)
// ---------------------------------------------------------------------------

namespace {

/// Runs f in a new thread, whose warm-start cache of the saturation solver is empty (a cold call); an exception of f is
/// rethrown here. The Catch2 assertions stay in the calling thread.
void inFreshThread(const std::function<void()>& f) {
    std::exception_ptr error;
    std::thread worker([&] {
        try {
            f();
        } catch (...) {
            error = std::current_exception();
        }
    });
    worker.join();
    if (error) std::rethrow_exception(error);
}

/// dp/drho at fixed T and x: positive for a mechanically stable phase.
double pressureSlope(double T, double rho, double x) { return stateJets(T, rho, x).p.gradient(1); }

/// The result of a flash, or its error message.
struct FlashOutcome {
    bool ok = false;
    FlashResult result;
    std::string error;
};

FlashOutcome tryFlash(const int combination[3], const double in[3]) {
    FlashOutcome out;
    try {
        out.result = flash(static_cast<Prop>(combination[0]), static_cast<Prop>(combination[1]),
                           static_cast<Prop>(combination[2]), in[0], in[1], in[2]);
        out.ok = true;
    } catch (const std::runtime_error& e) {
        out.error = e.what();
    }
    return out;
}

/// The largest difference of the eight values of two states, relative to the value (q and values below 1e-3 absolute).
double largestDifference(const FlashResult& a, const FlashResult& b) {
    double d = 0.0;
    for (int k = 0; k < 8; ++k)
        d = std::max(d, std::abs(a.value[k] - b.value[k]) / std::max(std::abs(b.value[k]), k == 7 ? 1.0 : 1e-3));
    return d;
}

}  // namespace

TEST_CASE("Saturation: only mechanically stable pairs, from a cold or a warm start", "[nh3h2o][equilibrium]") {
    // The T-p equilibrium at the bubble pressure of X = 0.9 at 400 K has an unstable root besides the bubble pair (a
    // "liquid" of 21109 mol/m3 with dp/drho < 0); from a cold cache it used to be returned.
    const double x = coolsolve::nh3h2o::moleFractionOf(0.9);
    Saturation bubble{}, pair{};
    inFreshThread([&] { bubble = saturation(SatSpec::BubbleT, 400.0, x); });
    inFreshThread([&] { pair = saturation(SatSpec::TP, 400.0, bubble.p); });
    CAPTURE(bubble.p, bubble.rhoL, pair.rhoL);
    CHECK(std::abs(pair.rhoL / bubble.rhoL - 1.0) <= 1e-8);
    CHECK(std::abs(pair.xV - bubble.xV) <= 1e-8);
    CHECK(pressureSlope(pair.T, pair.rhoL, pair.xL) > 0.0);
    CHECK(pressureSlope(pair.T, pair.rhoV, pair.xV) > 0.0);

    // a sweep below the bubble pressures, in one thread (each pair starts from the previous one)
    std::vector<Saturation> pairs;
    std::vector<std::string> errors;
    inFreshThread([&] {
        for (double T : {300.0, 350.0, 400.0})
            for (double xi : {0.5, 0.9, 0.99}) {
                const double pb = saturation(SatSpec::BubbleT, T, coolsolve::nh3h2o::moleFractionOf(xi)).p;
                for (double f : {0.9, 0.925, 0.95, 0.975, 1.0}) {
                    try {
                        pairs.push_back(saturation(SatSpec::TP, T, f * pb));
                    } catch (const std::runtime_error& e) {
                        errors.push_back(e.what());
                    }
                }
            }
    });
    std::string messages;
    for (const std::string& e : errors) messages += e + "\n";
    INFO(messages);
    CHECK(errors.empty());
    CHECK(pairs.size() == 45);
    for (const Saturation& s : pairs) {
        CAPTURE(s.T, s.p, s.xL, s.rhoL, s.rhoV);
        CHECK(pressureSlope(s.T, s.rhoL, s.xL) > 0.0);
        CHECK(pressureSlope(s.T, s.rhoV, s.xV) > 0.0);
    }
}

TEST_CASE("NH3H2O flash: (T, x, v) near the critical locus is two-phase, cold and after a warm-up", "[nh3h2o][flash]") {
    // From a cold cache the bracket used the unstable root above and threw "no two-phase state"; after the warm-up below
    // the bubble point came out unstable and the call returned a single-phase liquid at 4.675 MPa without an error.
    FlashResult reference;
    inFreshThread([&] { reference = flash(Prop::T, Prop::X, Prop::Q, 400.0, 0.9, 0.05); });
    REQUIRE(reference.twoPhase);
    FlashResult cold, warm;
    inFreshThread([&] { cold = flash(Prop::T, Prop::X, Prop::V, 400.0, 0.9, reference.volume()); });
    inFreshThread([&] {
        flash(Prop::T, Prop::X, Prop::Q, 400.0, 0.99, 0.0);   // the warm-up: a bubble point of 0.99 at 400 K
        warm = flash(Prop::T, Prop::X, Prop::V, 400.0, 0.9, reference.volume());
    });
    for (const FlashResult* r : {&cold, &warm}) {
        CAPTURE(r == &cold ? "cold" : "warm", r->pressure(), r->quality());
        CHECK(r->twoPhase);
        CHECK(std::abs(r->pressure() / reference.pressure() - 1.0) <= 1e-8);
        CHECK(std::abs(r->quality() - 0.05) <= 1e-8);
    }
}

TEST_CASE("NH3H2O flash: the same state in a fresh thread and after any warm-up", "[nh3h2o][flash]") {
    double volumeB = 0.0, pressureD = 0.0, lowP = 0.0, lowS = 0.0;
    inFreshThread([&] {
        volumeB = flash(Prop::T, Prop::X, Prop::Q, 400.0, 0.9, 0.05).volume();
        const double xD = coolsolve::nh3h2o::moleFractionOf(0.999);
        pressureD = std::sqrt(saturation(SatSpec::BubbleT, 400.0, xD).p * saturation(SatSpec::DewT, 400.0, xD).p);
        const FlashResult low = flash(Prop::T, Prop::X, Prop::Q, 270.0, 0.1, 0.05);
        lowP = low.pressure();
        lowS = low.entropy();
    });
    struct Case {
        int combination[3];
        double inputs[3];
    };
    const std::vector<Case> cases = {
        {{1, 3, 7}, {400.0, 0.9, volumeB}},      // two-phase near the critical locus (the unstable root)
        {{1, 2, 3}, {400.0, pressureD, 0.999}},  // inside the narrow dome of 0.999 at 400 K
        {{2, 3, 5}, {lowP, 0.1, lowS}},          // two-phase below 2 kPa
        {{1, 2, 3}, {330.0, 1e6, 0.5}},
        {{2, 3, 4}, {8e5, 0.5533, 7e5}},
        {{1, 3, 8}, {350.0, 0.3, 0.5}},
        {{2, 3, 8}, {1e6, 0.5, 0.3}},
        {{2, 4, 8}, {8e5, 707713.7, 0.3}},
        {{1, 5, 8}, {300.0, 4969.5, 0.7}},
        {{1, 2, 3}, {300.0, 5e5, 0.4}},          // sub-cooled liquid
        {{2, 3, 4}, {1e6, 0.4, 1.5e5}},          // sub-cooled liquid
        {{1, 2, 3}, {263.15, 2.9e5, 0.9995}},    // near-pure refrigerant
    };
    std::vector<FlashOutcome> cold(cases.size());
    for (size_t i = 0; i < cases.size(); ++i)
        inFreshThread([&] { cold[i] = tryFlash(cases[i].combination, cases[i].inputs); });
    const std::vector<std::function<void()>> warmUps = {
        [] { flash(Prop::T, Prop::X, Prop::Q, 400.0, 0.99, 0.0); },
        [] {
            for (double T : {250.0, 300.0, 350.0, 400.0})
                for (double xi : {0.1, 0.5, 0.9, 0.999})
                    for (double q : {0.0, 0.5, 1.0}) try {
                            flash(Prop::T, Prop::X, Prop::Q, T, xi, q);
                        } catch (const std::runtime_error&) {
                        }
        },
        [] {   // around the bubble pressures near the critical locus
            for (double T : {400.0, 380.0, 360.0})
                for (double xi : {0.99, 0.95, 0.9})
                    for (double f : {0.95, 1.0, 1.05}) try {
                            const double pb = saturation(SatSpec::BubbleT, T, coolsolve::nh3h2o::moleFractionOf(xi)).p;
                            flash(Prop::T, Prop::P, Prop::X, T, f * pb, xi);
                        } catch (const std::runtime_error&) {
                        }
        },
    };
    for (size_t w = 0; w < warmUps.size(); ++w) {
        std::vector<FlashOutcome> warm(cases.size());
        inFreshThread([&] {
            warmUps[w]();
            for (size_t i = 0; i < cases.size(); ++i) warm[i] = tryFlash(cases[i].combination, cases[i].inputs);
        });
        for (size_t i = 0; i < cases.size(); ++i) {
            const Case& c = cases[i];
            CAPTURE(w, c.combination[0], c.combination[1], c.combination[2], c.inputs[0], c.inputs[1], c.inputs[2]);
            INFO("cold: " << (cold[i].ok ? "a state" : cold[i].error) << "; warm: " << (warm[i].ok ? "a state" : warm[i].error));
            REQUIRE(warm[i].ok == cold[i].ok);
            if (!cold[i].ok) {
                CHECK(warm[i].error == cold[i].error);
                continue;
            }
            CHECK(warm[i].result.twoPhase == cold[i].result.twoPhase);
            CHECK(largestDifference(warm[i].result, cold[i].result) <= 1e-9);
        }
    }
    for (size_t i = 0; i < cases.size(); ++i) {   // these states all exist
        INFO(cold[i].error);
        CHECK(cold[i].ok);
    }
}

TEST_CASE("NH3H2O flash: (T, P, x) strictly inside the dome is two-phase", "[nh3h2o][flash]") {
    // a saturation that does not converge must not be taken for the absence of a two-phase region
    int checked = 0;
    for (double T = 250.0; T <= 400.0 + 1e-9; T += 25.0)
        for (double xi : {0.1, 0.5, 0.9, 0.999}) {
            const double x = coolsolve::nh3h2o::moleFractionOf(xi);
            const double pb = saturation(SatSpec::BubbleT, T, x).p, pd = saturation(SatSpec::DewT, T, x).p;
            for (double f : {0.1, 0.5, 0.9}) {
                const double p = pd * std::pow(pb / pd, f);
                const int combination[3] = {1, 2, 3};
                const double in[3] = {T, p, xi};
                const FlashOutcome r = tryFlash(combination, in);
                CAPTURE(T, xi, pd, pb, p);
                INFO(r.error);
                REQUIRE(r.ok);
                CHECK(r.result.twoPhase);
                CHECK(r.result.quality() > 0.0);
                CHECK(r.result.quality() < 1.0);
                ++checked;
            }
        }
    CHECK(checked == 84);
}

TEST_CASE("NH3H2O flash: two-phase and vapour states below 2 kPa", "[nh3h2o][flash]") {
    // the vapour density is below 1 mol/m3 there, the former lower bound of the polish
    const FlashResult state = flash(Prop::T, Prop::X, Prop::Q, 270.0, 0.1, 0.05);
    REQUIRE(state.twoPhase);
    CHECK(state.pressure() == Approx(1985.24).epsilon(1e-5));
    const FlashResult byVolume = flash(Prop::T, Prop::X, Prop::V, 270.0, 0.1, state.volume());
    const FlashResult byEntropy = flash(Prop::P, Prop::X, Prop::S, state.pressure(), 0.1, state.entropy());
    for (const FlashResult* r : {&byVolume, &byEntropy}) {
        CAPTURE(r->pressure(), r->quality());
        CHECK(r->twoPhase);
        CHECK(std::abs(r->pressure() / state.pressure() - 1.0) <= 1e-8);
        CHECK(std::abs(r->quality() - 0.05) <= 1e-8);
    }
    // a vapour of 0.3 mol/m3: the density root is not converged to the polish's criterion, so the polish takes Newton
    // steps from it (they were clamped to 1 mol/m3)
    const FlashResult vapour = flash(Prop::T, Prop::P, Prop::X, 400.0, 1000.0, 0.05);
    CHECK_FALSE(vapour.twoPhase);
    CHECK(vapour.quality() == 1.001);
    CHECK(1.0 / (vapour.volume() * coolsolve::nh3h2o::molarMass(coolsolve::nh3h2o::moleFractionOf(0.05))) < 0.31);
    const FlashResult again = flash(Prop::P, Prop::X, Prop::H, 1000.0, 0.05, vapour.enthalpy());
    CHECK(std::abs(again.temperature() / 400.0 - 1.0) <= 1e-9);
    CHECK(std::abs(again.volume() / vapour.volume() - 1.0) <= 1e-9);
}

TEST_CASE("NH3H2O flash: scan of (T, P, x) inside the dome, cold against warm", "[.nh3h2o-scan]") {
    // the larger version of the grid above, up to the critical locus (about half a minute): every state inside the dome
    // is two-phase or an error that does not say it is outside the two-phase region (on 2026-10-19: 2 errors, 240 K
    // and X <= 0.1 at 44 to 56 Pa, below the freezing line of the water-rich liquid), the same in a fresh thread as in
    // the warm sequence
    int twoPhase = 0, errors = 0, points = 0;
    for (double T = 230.0; T <= 420.0 + 1e-9; T += 10.0)
        for (double xi : {0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.995, 0.999}) {
            const double x = coolsolve::nh3h2o::moleFractionOf(xi);
            double pb = 0.0, pd = 0.0;
            try {
                pb = saturation(SatSpec::BubbleT, T, x).p;
                pd = saturation(SatSpec::DewT, T, x).p;
            } catch (const std::runtime_error&) {
                continue;   // no dome at this temperature (above the critical locus), or no saturation found
            }
            for (double f : {0.02, 0.25, 0.5, 0.75, 0.98}) {
                const double p = pd * std::pow(pb / pd, f);
                const int combination[3] = {1, 2, 3};
                const double in[3] = {T, p, xi};
                const FlashOutcome warm = tryFlash(combination, in);
                FlashOutcome cold;
                inFreshThread([&] { cold = tryFlash(combination, in); });
                CAPTURE(T, xi, pd, pb, p);
                INFO("warm: " << (warm.ok ? "a state" : warm.error) << "; cold: " << (cold.ok ? "a state" : cold.error));
                ++points;
                if (warm.ok) {
                    CHECK(warm.result.twoPhase);
                    ++twoPhase;
                } else {
                    CHECK_THAT(warm.error, !ContainsSubstring("outside"));
                    ++errors;
                }
                REQUIRE(cold.ok == warm.ok);
                if (warm.ok) CHECK(largestDifference(warm.result, cold.result) <= 1e-9);
            }
        }
    INFO(points << " points: " << twoPhase << " two-phase, " << errors << " errors");
    CHECK(points > 900);
    CHECK(errors <= 2);
}
