#include "coolsolve/nh3h2o.h"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace coolsolve::nh3h2o {

using taylor::Jet;

namespace {

// Coefficients of IAPWS G4-01 (Tillner-Roth & Friend 1998), Tables 1 to 3 and A1, and of
// IAPWS-95 for the water residual (published data, cross-checked against CoolProp's
// dev/fluids/Water.json and the TillnerRoth-DKV-1993 entry of dev/fluids/Ammonia.json).

// Table 2: the ideal part, Eq. (2). a[1..8] water, a[9..14] ammonia, theta[4..8] the Planck-
// Einstein frequencies of water, t[12..14] the powers of tau0 of ammonia.
constexpr double kIdealA[15] = {0.0,
                                -7.720435,
                                8.649358,
                                3.006320,
                                0.012436,
                                0.97315,
                                1.279500,
                                0.969560,
                                0.248730,
                                -16.444285,
                                4.036946,
                                -1.0,
                                10.69955,
                                -1.775436,
                                0.82374034};
constexpr double kIdealTheta[9] = {0.0, 0.0, 0.0, 0.0, 1.666, 4.578, 10.018, 11.964, 35.600};
constexpr double kIdealT[15] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, -1.5, -1.75};

// Tables 1 and 3: the reducing functions, Eqs. (4)-(7).
constexpr double kKt = 0.9648407;
constexpr double kKv = 1.2395117;
constexpr double kAlpha = 1.125455;    ///< exponent of the cross term of Tn(x)
constexpr double kBeta = 0.8978069;    ///< exponent of the cross term of rhon(x)
constexpr double kGamma = 0.5248379;   ///< exponent of the composition factor of Eq. (8)
constexpr double kTc12 = kKt * (kTc1 + kTc2) / 2.0;
constexpr double kInvRhoc12 = kKv * (1.0 / kRhoc1 + 1.0 / kRhoc2) / 2.0;  ///< 1/rhoc12 [m3/mol]

/// One term n tau^t delta^d exp(-delta^e), e = 0 meaning no exponential factor.
struct Term {
    double n;
    int d;
    double t;
    int e;
};

// Table 3: the departure function, Eq. (8) (t1..t14, d1..d14, e1..e14).
constexpr Term kDeparture[14] = {
    {-1.855822e-02, 4, 1.5, 0},  {5.258010e-02, 5, 0.5, 1},   {3.552874e-10, 15, 6.5, 1}, {5.451379e-06, 12, 1.75, 1},
    {-5.998546e-13, 12, 15.0, 1}, {-3.687808e-06, 15, 6.0, 2}, {2.586192e-01, 4, -1.0, 1}, {-1.368072e-08, 15, 4.0, 1},
    {1.226146e-02, 4, 3.5, 1},   {-7.181443e-02, 5, 0.0, 1},  {9.970849e-02, 6, -1.0, 2},  {1.0584086e-03, 10, 8.0, 2},
    {-1.963687e-01, 6, 7.5, 2},  {-7.777897e-01, 2, 4.0, 2},
};

// Table A1: the residual part of pure ammonia, Eq. (A1) (the coefficients of the guideline, to
// its 8 digits, identical to CoolProp's TillnerRoth-DKV-1993 entry).
constexpr Term kAmmonia[21] = {
    {-1.858814e+00, 1, 1.5, 0},  {4.5544310e-02, 2, -0.5, 0}, {7.2385480e-01, 1, 0.5, 0},  {1.2294700e-02, 4, 1.0, 0},
    {2.1418820e-11, 15, 3.0, 0}, {-1.4300200e-02, 3, 0.0, 1}, {3.4413240e-01, 3, 3.0, 1},  {-2.8735710e-01, 1, 4.0, 1},
    {2.3525890e-05, 8, 4.0, 1},  {-3.4971110e-02, 2, 5.0, 1}, {1.8311170e-03, 8, 5.0, 2},  {2.3978520e-02, 1, 3.0, 2},
    {-4.0853750e-02, 1, 6.0, 2}, {2.3792750e-01, 2, 8.0, 2},  {-3.5489720e-02, 3, 8.0, 2}, {-1.8237290e-01, 2, 10.0, 2},
    {2.2815560e-02, 4, 10.0, 2}, {-6.6634440e-03, 3, 5.0, 3}, {-8.8474860e-03, 1, 7.5, 3}, {2.2726350e-03, 2, 15.0, 3},
    {-5.5886550e-04, 4, 30.0, 3},
};

/// One power-type term n delta^d tau^t exp(-delta^l) of IAPWS-95 (l = 0: no exponential).
struct WaterPowerTerm {
    double n;
    int d;
    double t;
    int l;
};

/// One Gaussian term n delta^d tau^t exp(-alpha (delta - epsilon)^2 - beta (tau - gamma)^2).
struct WaterGaussianTerm {
    double n;
    int d;
    double t;
    double alpha;
    double epsilon;
    double beta;
    double gamma;
};

/// One non-analytic term n DELTA^b PSI with theta = (1 - tau) + A ((delta-1)^2)^(1/(2 beta)),
/// DELTA = theta^2 + B ((delta-1)^2)^a and PSI = exp(-C (delta-1)^2 - D (tau-1)^2).
struct WaterNonAnalyticTerm {
    double n;
    double a;
    double b;
    double beta;
    double A;
    double B;
    double C;
    double D;
};

// The residual Helmholtz energy of IAPWS-95 for water (Wagner & Pruss 2002): 51 power-type
// terms (7 without and 44 with the exponential), 3 Gaussian and 2 non-analytic terms; published
// coefficients, as held in CoolProp's dev/fluids/Water.json ("alphar", same numbers).
constexpr WaterPowerTerm kWaterPower[51] = {
    {0.012533547935523, 1, -0.5, 0},   {7.8957634722828, 1, 0.875, 0},   {-8.7803203303561, 1, 1, 0},
    {0.31802509345418, 2, 0.5, 0},     {-0.26145533859358, 2, 0.75, 0},  {-0.0078199751687981, 3, 0.375, 0},
    {0.0088089493102134, 4, 1, 0},     {-0.66856572307965, 1, 4, 1},     {0.20433810950965, 1, 6, 1},
    {-6.6212605039687e-05, 1, 12, 1},  {-0.19232721156002, 2, 1, 1},     {-0.25709043003438, 2, 5, 1},
    {0.16074868486251, 3, 4, 1},       {-0.040092828925807, 4, 2, 1},    {3.9343422603254e-07, 4, 13, 1},
    {-7.5941377088144e-06, 5, 9, 1},   {0.00056250979351888, 7, 3, 1},   {-1.5608652257135e-05, 9, 4, 1},
    {1.1537996422951e-09, 10, 11, 1},  {3.6582165144204e-07, 11, 4, 1},  {-1.3251180074668e-12, 13, 13, 1},
    {-6.2639586912454e-10, 15, 1, 1},  {-0.10793600908932, 1, 7, 2},     {0.017611491008752, 2, 1, 2},
    {0.22132295167546, 2, 9, 2},       {-0.40247669763528, 2, 10, 2},    {0.58083399985759, 3, 10, 2},
    {0.0049969146990806, 4, 3, 2},     {-0.031358700712549, 4, 7, 2},    {-0.74315929710341, 4, 10, 2},
    {0.4780732991548, 5, 10, 2},       {0.020527940895948, 6, 6, 2},     {-0.13636435110343, 6, 10, 2},
    {0.014180634400617, 7, 10, 2},     {0.0083326504880713, 9, 1, 2},    {-0.029052336009585, 9, 2, 2},
    {0.038615085574206, 9, 3, 2},      {-0.020393486513704, 9, 4, 2},    {-0.0016554050063734, 9, 8, 2},
    {0.0019955571979541, 10, 6, 2},    {0.00015870308324157, 10, 9, 2},  {-1.638856834253e-05, 12, 8, 2},
    {0.043613615723811, 3, 16, 3},     {0.034994005463765, 4, 22, 3},    {-0.076788197844621, 4, 23, 3},
    {0.022446277332006, 5, 23, 3},     {-6.2689710414685e-05, 14, 10, 4}, {-5.5711118565645e-10, 3, 50, 6},
    {-0.19905718354408, 6, 44, 6},     {0.31777497330738, 6, 46, 6},     {-0.11841182425981, 6, 50, 6},
};

constexpr WaterGaussianTerm kWaterGaussian[3] = {
    {-31.306260323435, 3, 0, 20.0, 1.0, 150.0, 1.21},
    {31.546140237781, 3, 1, 20.0, 1.0, 150.0, 1.21},
    {-2521.3154341695, 3, 4, 20.0, 1.0, 250.0, 1.25},
};

constexpr WaterNonAnalyticTerm kWaterNonAnalytic[2] = {
    {-0.14874640856724, 3.5, 0.85, 0.3, 0.32, 0.2, 28.0, 700.0},
    {0.31806110878444, 3.5, 0.95, 0.3, 0.32, 0.2, 32.0, 800.0},
};

/// The value of a scalar or of a jet (the point the Taylor expansion is about).
double scalarOf(double v) { return v; }
template <int Order>
double scalarOf(const Jet<Order>& jet) {
    return jet.value();
}

/// The zero of the type of `v` (0.0 or the zero jet): the pure-limit convention of the terms
/// whose derivative does not exist at an endpoint of the domain (see the header).
double makeZero(double) { return 0.0; }
template <int Order>
Jet<Order> makeZero(const Jet<Order>&) {
    return Jet<Order>::constant(0.0);
}

/// v ln v with the limit 0 at v = 0 (the mixing entropy terms of Eq. (2) at the pure limits).
template <class S>
S vLogV(const S& v) {
    if (scalarOf(v) == 0.0) return makeZero(v);
    return v * log(v);
}

}  // namespace

template <class S>
S reducingTemperature(const S& x) {
    using std::pow;
    const S cross = x - pow(x, 1.0 + kAlpha);  // x (1 - x^alpha)
    return (1.0 - x) * (1.0 - x) * kTc1 + x * x * kTc2 + 2.0 * cross * kTc12;
}

template <class S>
S reducingDensity(const S& x) {
    using std::pow;
    const S cross = x - pow(x, 1.0 + kBeta);  // x (1 - x^beta)
    const S inverse = (1.0 - x) * (1.0 - x) / kRhoc1 + x * x / kRhoc2 + 2.0 * cross * kInvRhoc12;
    return 1.0 / inverse;
}

template <class S>
S idealHelmholtz(const S& tau0, const S& delta0, const S& x) {
    using std::exp;
    using std::log;
    using std::pow;
    S water = kIdealA[1] + kIdealA[2] * tau0 + kIdealA[3] * log(tau0);
    for (int i = 4; i <= 8; ++i) water = water + kIdealA[i] * log(1.0 - exp(-kIdealTheta[i] * tau0));
    S ammonia = kIdealA[9] + kIdealA[10] * tau0 + kIdealA[11] * log(tau0);
    for (int i = 12; i <= 14; ++i) ammonia = ammonia + kIdealA[i] * pow(tau0, kIdealT[i]);
    return log(delta0) + (1.0 - x) * water + vLogV(1.0 - x) + x * ammonia + vLogV(x);
}

namespace {

/// n tau^t delta^d exp(-delta^e) for one table term (e = 0: no exponential).
template <class S>
S evaluateTerm(const Term& term, const S& tau, const S& delta) {
    using std::exp;
    using std::pow;
    S value = term.n * pow(delta, term.d) * pow(tau, term.t);
    if (term.e > 0) value = value * exp(-pow(delta, term.e));
    return value;
}

}  // namespace

template <class S>
S residualAmmonia(const S& tau, const S& delta) {
    S phi = makeZero(tau);
    for (const Term& term : kAmmonia) phi = phi + evaluateTerm(term, tau, delta);
    return phi;
}

template <class S>
S departureFunction(const S& tau, const S& delta, const S& x) {
    using std::pow;
    S bracket = evaluateTerm(kDeparture[0], tau, delta);
    for (int i = 1; i <= 5; ++i) bracket = bracket + evaluateTerm(kDeparture[i], tau, delta);
    S xWeighted = makeZero(tau);
    for (int i = 6; i <= 12; ++i) xWeighted = xWeighted + evaluateTerm(kDeparture[i], tau, delta);
    bracket = bracket + x * xWeighted + x * x * evaluateTerm(kDeparture[13], tau, delta);
    const S factor = x - pow(x, 1.0 + kGamma);  // x (1 - x^gamma)
    return factor * bracket;
}

template <class S>
S residualWater(const S& tau, const S& delta) {
    using std::exp;
    using std::pow;
    S phi = makeZero(tau);
    for (const WaterPowerTerm& term : kWaterPower) {
        S value = term.n * pow(delta, term.d) * pow(tau, term.t);
        if (term.l > 0) value = value * exp(-pow(delta, term.l));
        phi = phi + value;
    }
    for (const WaterGaussianTerm& term : kWaterGaussian) {
        const S d = delta - term.epsilon;
        const S t = tau - term.gamma;
        phi = phi + term.n * pow(delta, term.d) * pow(tau, term.t) * exp(-term.alpha * (d * d) - term.beta * (t * t));
    }
    // The non-analytic terms: every power of (delta-1)^2 is C^3 in delta at delta = 1 (its
    // derivatives are 0 there), which the zero-value rule of taylor::pow returns; at the exact
    // critical point (tau = delta = 1) the same convention sets the divergent higher tau
    // derivatives of DELTA^b (b < 1) to 0, the value and the first derivatives being 0 as well.
    for (const WaterNonAnalyticTerm& term : kWaterNonAnalytic) {
        const S d = delta - 1.0;
        const S u = d * d;
        const S theta = (1.0 - tau) + term.A * pow(u, 1.0 / (2.0 * term.beta));
        const S mixture = theta * theta + term.B * pow(u, term.a);
        const S psi = exp(-term.C * u - term.D * (tau - 1.0) * (tau - 1.0));
        const S singular = scalarOf(mixture) == 0.0 ? makeZero(mixture) : pow(mixture, term.b);
        phi = phi + term.n * singular * psi;
    }
    return phi;
}

namespace {

[[noreturn]] void domainError(double T, double rho, double x, const char* what) {
    std::ostringstream os;
    os << "nh3h2o: (T = " << T << " K, rho = " << rho << " mol/m3, x = " << x << "): " << what;
    throw std::runtime_error(os.str());
}

void checkDomain(double T, double rho, double x) {
    if (!std::isfinite(T) || !std::isfinite(rho) || !std::isfinite(x)) domainError(T, rho, x, "an argument is not finite");
    if (!(T > 0.0)) domainError(T, rho, x, "the temperature must be positive [K]");
    if (!(rho > 0.0)) domainError(T, rho, x, "the molar density must be positive [mol/m3]");
    if (!(x >= 0.0 && x <= 1.0)) domainError(T, rho, x, "the mole fraction of ammonia must be in [0, 1]");
}

/// Everything `state` and `stateJets` compute, as exact jets in (T, rho, x) from the order-3 jet
/// of Phi (f = R T Phi, p = rho^2 df/drho, u = f - T df/dT, s = -df/dT, h = u + p/rho).
struct Core {
    Jet<2> p;
    Jet<2> u;
    Jet<2> h;
    Jet<2> s;
    double f;
    double cv;
    double cp;
    double w;
};

Core computeCore(double T, double rho, double x) {
    const Jet<3> jetT = Jet<3>::variable(T, 0);
    const Jet<3> jetRho = Jet<3>::variable(rho, 1);
    const Jet<3> jetX = Jet<3>::variable(x, 2);
    const Jet<3> phi = reducedHelmholtz(jetT, jetRho, jetX);

    const Jet<3> f3 = kR * jetT * phi;
    const Jet<2> f2(f3);
    const Jet<2> fT = derivative(f3, 0);
    const Jet<2> fRho = derivative(f3, 1);
    const Jet<2> rho2(jetRho);

    Core core;
    core.p = rho2 * rho2 * fRho;
    core.u = f2 - Jet<2>(jetT) * fT;
    core.s = -fT;
    core.h = core.u + core.p / rho2;
    core.f = f2.value();
    core.cv = derivative(core.u, 0).value();
    const double hT = derivative(core.h, 0).value(), hRho = derivative(core.h, 1).value();
    const double pT = derivative(core.p, 0).value(), pRho = derivative(core.p, 1).value();
    core.cp = hT - hRho * pT / pRho;
    core.w = std::sqrt(core.cp / core.cv * pRho / molarMass(x));
    return core;
}

}  // namespace

template <class S>
S reducedHelmholtz(const S& T, const S& rho, const S& x) {
    checkDomain(scalarOf(T), scalarOf(rho), scalarOf(x));
    const S tau0 = kT0 / T;
    const S delta0 = rho / kRho0;
    const S tau = reducingTemperature(x) / T;
    const S delta = rho / reducingDensity(x);
    return idealHelmholtz(tau0, delta0, x) + (1.0 - x) * residualWater(tau, delta) +
           x * residualAmmonia(tau, delta) + departureFunction(tau, delta, x);
}

State state(double T, double rho, double x) {
    const Core core = computeCore(T, rho, x);
    return {core.p.value(), core.f, core.u.value(), core.h.value(), core.s.value(), core.cv, core.cp, core.w};
}

StateJets stateJets(double T, double rho, double x) {
    const Core core = computeCore(T, rho, x);
    return {core.p, core.u, core.h, core.s};
}

// The formulation is written once as a template; the instantiations the library and the tests
// use are made here so that the coefficients stay in this translation unit.
template double idealHelmholtz(const double&, const double&, const double&);
template double residualWater(const double&, const double&);
template double residualAmmonia(const double&, const double&);
template double reducingTemperature(const double&);
template double reducingDensity(const double&);
template double departureFunction(const double&, const double&, const double&);
template double reducedHelmholtz(const double&, const double&, const double&);
template Jet<1> idealHelmholtz(const Jet<1>&, const Jet<1>&, const Jet<1>&);
template Jet<1> residualWater(const Jet<1>&, const Jet<1>&);
template Jet<1> residualAmmonia(const Jet<1>&, const Jet<1>&);
template Jet<1> reducingTemperature(const Jet<1>&);
template Jet<1> reducingDensity(const Jet<1>&);
template Jet<1> departureFunction(const Jet<1>&, const Jet<1>&, const Jet<1>&);
template Jet<1> reducedHelmholtz(const Jet<1>&, const Jet<1>&, const Jet<1>&);
template Jet<2> idealHelmholtz(const Jet<2>&, const Jet<2>&, const Jet<2>&);
template Jet<2> residualWater(const Jet<2>&, const Jet<2>&);
template Jet<2> residualAmmonia(const Jet<2>&, const Jet<2>&);
template Jet<2> reducingTemperature(const Jet<2>&);
template Jet<2> reducingDensity(const Jet<2>&);
template Jet<2> departureFunction(const Jet<2>&, const Jet<2>&, const Jet<2>&);
template Jet<2> reducedHelmholtz(const Jet<2>&, const Jet<2>&, const Jet<2>&);
template Jet<3> idealHelmholtz(const Jet<3>&, const Jet<3>&, const Jet<3>&);
template Jet<3> residualWater(const Jet<3>&, const Jet<3>&);
template Jet<3> residualAmmonia(const Jet<3>&, const Jet<3>&);
template Jet<3> reducingTemperature(const Jet<3>&);
template Jet<3> reducingDensity(const Jet<3>&);
template Jet<3> departureFunction(const Jet<3>&, const Jet<3>&, const Jet<3>&);
template Jet<3> reducedHelmholtz(const Jet<3>&, const Jet<3>&, const Jet<3>&);

}  // namespace coolsolve::nh3h2o
