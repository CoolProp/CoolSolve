#include "coolsolve/libr_ashrae.h"

#include "coolsolve/ees_string.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace coolsolve::libr_ashrae {

namespace {

// Coefficients of the ASHRAE Handbook (1989) correlations, in increasing powers of X (mass %); see libr_ashrae.h.
// Duhring relation t = B(X) + t' A(X)
constexpr double kDuhringA[4] = {-2.00755, 0.16976, -3.133362e-3, 1.97668e-5};
constexpr double kDuhringB[4] = {124.937, -7.71649, 0.152286, -7.95090e-4};
// log10(p [kPa]) = C + D/T' + E/T'^2 with T' = t' + 273.15 K, the vapour pressure of the refrigerant
constexpr double kVapourC = 7.05, kVapourD = -1596.49, kVapourE = -104095.5;
constexpr double kKelvin = 273.15;
// Enthalpy [kJ/kg] = Ah(X) + t Bh(X) + t^2 Ch(X) (McNeely)
constexpr double kEnthalpyA[5] = {-2024.33, 163.309, -4.88161, 6.302948e-2, -2.913705e-4};
constexpr double kEnthalpyB[5] = {18.2829, -1.1691757, 3.248041e-2, -4.034184e-4, 1.8520569e-6};
constexpr double kEnthalpyC[5] = {-3.7008214e-2, 2.8877666e-3, -8.1313015e-5, 9.9116628e-7, -4.4441207e-9};

constexpr double kLn10 = 2.30258509299404568402;

// What a call of a known but not implemented routine is told (see Routine::notImplemented)
constexpr const char* kAshraeMissing =
    "EES LiBr-H2O routine not implemented in CoolSolve yet (register CS-GAP-FLUIDS-ABS); the ASHRAE routines X_LIBR, "
    "H_LIBR, T_LIBR, P_LIBR are available";
constexpr const char* kPatekKlomfarMissing =
    "EES LiBr-H2O routine not implemented in CoolSolve yet (register CS-GAP-LIBR-PK); the ASHRAE routines X_LIBR, "
    "H_LIBR, T_LIBR, P_LIBR are available";

/// The routines: the four implemented ones (the descriptors of the messages below), then the known EES routines of
/// the absorption pairs that are not implemented: the rest of the ASHRAE library, the Patek-Klomfar library
/// (functions `*_LiBrH2O`, procedure `Q_LiBrH2O`).
const Routine kRoutines[] = {
    {"X_LIBR", "T", "P", &concentration, nullptr},
    {"H_LIBR", "T", "X", &enthalpy, nullptr},
    {"T_LIBR", "P", "X", &temperature, nullptr},
    {"P_LIBR", "T", "X", &pressure, nullptr},
    {"V_LIBR", nullptr, nullptr, nullptr, kAshraeMissing},
    {"Q_LIBR", nullptr, nullptr, nullptr, kAshraeMissing},
    {"VISC_LIBR", nullptr, nullptr, nullptr, kAshraeMissing},
    {"COND_LIBR", nullptr, nullptr, nullptr, kAshraeMissing},
    {"T_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"P_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"x_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"h_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"s_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"cp_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"rho_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"Visc_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"Cond_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"MassFraction_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"MoleFraction_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
    {"Q_LiBrH2O", nullptr, nullptr, nullptr, kPatekKlomfarMissing},
};
const Routine& kConcentrationRoutine = kRoutines[0];
const Routine& kEnthalpyRoutine = kRoutines[1];
const Routine& kTemperatureRoutine = kRoutines[2];
const Routine& kPressureRoutine = kRoutines[3];

bool endsWith(const std::string& s, const char* suffix) {
    const size_t n = std::char_traits<char>::length(suffix);
    return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}

std::string showValue(double v) {
    std::ostringstream os;
    os << v;
    return os.str();
}

/// A call as the messages name it, `X_LIBR('SI', T = 40, P = 0.87)`: the routine and its two numeric arguments.
struct CallText {
    const Routine& routine;
    double first, second;

    [[noreturn]] void refuse(const std::string& reason) const {
        throw std::runtime_error(std::string(routine.name) + "('SI', " + routine.first + " = " + showValue(first) + ", " +
                                 routine.second + " = " + showValue(second) + "): " + reason);
    }

    void requireFinite() const {
        if (!std::isfinite(first) || !std::isfinite(second)) refuse("an argument is not finite");
    }
};

/// Value and derivative of a polynomial with coefficients `c` (increasing powers) at `x`, by Horner's scheme.
struct Polynomial {
    double value, derivative;
};

template <size_t N>
Polynomial polynomial(const double (&c)[N], double x) {
    double value = c[N - 1], derivative = 0.0;
    for (size_t i = N - 1; i-- > 0;) {
        derivative = derivative * x + value;
        value = value * x + c[i];
    }
    return {value, derivative};
}

/// d log10(p) / d T' of the vapour-pressure equation, positive for every T' > 0.
double dLog10pdT(double kelvin) {
    return -kVapourD / (kelvin * kelvin) - 2.0 * kVapourE / (kelvin * kelvin * kelvin);
}

/// Refrigerant temperature t' [degC] of the pressure `p` [kPa] and its derivative dt'/dp: the positive root of
/// (C - log10 p) T'^2 + D T' + E = 0. The pressure is positive and finite (checked by the caller).
struct Refrigerant {
    double celsius, dCelsiusdP;
};

Refrigerant refrigerantAtPressure(const CallText& call, double p) {
    if (!(p > 0.0)) call.refuse("the pressure must be positive, in kPa");
    const double a = kVapourC - std::log10(p);
    if (!(a > 0.0))
        call.refuse("no refrigerant temperature belongs to this pressure: the vapour-pressure equation stops at 10^" +
                    showValue(kVapourC) + " kPa");
    const double kelvin = (-kVapourD + std::sqrt(kVapourD * kVapourD - 4.0 * a * kVapourE)) / (2.0 * a);
    return {kelvin - kKelvin, 1.0 / (p * kLn10 * dLog10pdT(kelvin))};
}

/// Solution temperature of the Duhring relation at the concentration `x` for the refrigerant temperature `tPrime`,
/// with the derivatives dt/dX and dt/dt' = A(X).
struct Duhring {
    double t, dtdX, dtdTPrime;
};

Duhring duhring(double x, double tPrime) {
    const Polynomial a = polynomial(kDuhringA, x), b = polynomial(kDuhringB, x);
    return {b.value + tPrime * a.value, b.derivative + tPrime * a.derivative, a.value};
}

/// The stationary points of the Duhring cubic at `tPrime`, dt/dX = c0 + c1 X + c2 X^2 = 0: roots of the quadratic
/// that are finite (0, 1 or 2 of them), by the form that does not lose digits to cancellation.
struct Roots {
    double value[2];
    int count = 0;
};

Roots stationaryPoints(double tPrime) {
    double c[3];
    for (int n = 0; n < 3; ++n) c[n] = (n + 1) * (kDuhringB[n + 1] + tPrime * kDuhringA[n + 1]);
    Roots roots;
    const double discriminant = c[1] * c[1] - 4.0 * c[2] * c[0];
    if (discriminant < 0.0) return roots;
    const double q = -0.5 * (c[1] + std::copysign(std::sqrt(discriminant), c[1]));
    for (const double root : {c[2] != 0.0 ? q / c[2] : INFINITY, q != 0.0 ? c[0] / q : INFINITY})
        if (std::isfinite(root)) roots.value[roots.count++] = root;
    return roots;
}

}  // namespace

const Routine* findRoutine(const std::string& lowerName) {
    // the names end with _libr or _librh2o: the table is searched only for those
    if (!(endsWith(lowerName, "_libr") || endsWith(lowerName, "_librh2o"))) return nullptr;
    for (const Routine& routine : kRoutines)
        if (eesStringEqual(routine.name, lowerName)) return &routine;
    return nullptr;
}

void requireImplemented(const Routine& routine) {
    if (!routine.evaluate) throw std::runtime_error(std::string(routine.name) + ": " + routine.notImplemented);
}

void checkUnitSystem(const Routine& routine, const std::string& unit) {
    const std::string units = "T in degC, P in kPa, X in mass %, h in kJ/kg";
    if (eesStringEqual(unit, "SI")) return;
    if (eesStringEqual(unit, "Eng"))
        throw std::runtime_error(std::string(routine.name) + ": 'Eng' (English units) is not supported by CoolSolve: use 'SI' (" +
                                 units + ")");
    throw std::runtime_error(std::string(routine.name) + ": the unit system '" + unit +
                             "' is not valid: EES takes 'SI' or 'Eng', and CoolSolve supports 'SI' only (" + units + ")");
}

Result enthalpy(double t, double x) {
    const CallText call{kEnthalpyRoutine, t, x};
    call.requireFinite();
    const Polynomial a = polynomial(kEnthalpyA, x), b = polynomial(kEnthalpyB, x), c = polynomial(kEnthalpyC, x);
    return {a.value + t * b.value + t * t * c.value,   // h
            b.value + 2.0 * t * c.value,               // dh/dT
            a.derivative + t * b.derivative + t * t * c.derivative};   // dh/dX
}

Result temperature(double p, double x) {
    const CallText call{kTemperatureRoutine, p, x};
    call.requireFinite();
    const Refrigerant refrigerant = refrigerantAtPressure(call, p);
    const Duhring solution = duhring(x, refrigerant.celsius);
    return {solution.t, solution.dtdTPrime * refrigerant.dCelsiusdP, solution.dtdX};
}

Result pressure(double t, double x) {
    const CallText call{kPressureRoutine, t, x};
    call.requireFinite();
    // t' from the Duhring relation, then the vapour-pressure equation
    const Polynomial a = polynomial(kDuhringA, x), b = polynomial(kDuhringB, x);
    if (a.value == 0.0) call.refuse("the coefficient A(X) of the Duhring relation is zero: no refrigerant temperature");
    const double tPrime = (t - b.value) / a.value;
    const double kelvin = tPrime + kKelvin;
    if (!(kelvin > 0.0))
        call.refuse("the refrigerant temperature " + showValue(tPrime) + " degC, from the Duhring relation, is below absolute zero");
    const double p = std::pow(10.0, kVapourC + kVapourD / kelvin + kVapourE / (kelvin * kelvin));
    const double dpdTPrime = p * kLn10 * dLog10pdT(kelvin);
    const double dTPrimedX = -(b.derivative + tPrime * a.derivative) / a.value;
    return {p, dpdTPrime / a.value, dpdTPrime * dTPrimedX};
}

Result concentration(double t, double p) {
    const CallText call{kConcentrationRoutine, t, p};
    call.requireFinite();
    const Refrigerant refrigerant = refrigerantAtPressure(call, p);
    const double tPrime = refrigerant.celsius;

    // The Duhring cubic is increasing between its local minimum and its local maximum: take the interval around
    // 60 %, where dt/dX > 0, limited to 0...100 % (it holds 37.6...82.8 % for every t' above -20 degC).
    constexpr double anchor = 60.0;
    if (!(duhring(anchor, tPrime).dtdX > 0.0))
        call.refuse("the Duhring relation does not increase with the concentration around 60 % at the refrigerant temperature " +
                    showValue(tPrime) + " degC");
    double lo = 0.0, hi = 100.0;
    const Roots stationary = stationaryPoints(tPrime);
    for (int i = 0; i < stationary.count; ++i) {
        if (stationary.value[i] < anchor) lo = std::max(lo, stationary.value[i]);
        else hi = std::min(hi, stationary.value[i]);
    }
    const double tLo = duhring(lo, tPrime).t, tHi = duhring(hi, tPrime).t;
    if (!(t >= tLo && t <= tHi))
        call.refuse("no concentration of the Duhring relation gives this temperature: for the refrigerant temperature " +
                    showValue(tPrime) + " degC the solution temperature is between " + showValue(tLo) + " degC (X = " +
                    showValue(lo) + " %) and " + showValue(tHi) + " degC (X = " + showValue(hi) + " %)");

    // Newton's method safeguarded by bisection (the bracket [lo, hi] shrinks with the sign of t(X) - t, increasing)
    double x = anchor;
    bool converged = false;
    Duhring state = duhring(x, tPrime);
    for (int iteration = 0; iteration < 100 && !converged; ++iteration) {
        const double g = state.t - t;
        if (g > 0.0) hi = x;
        else if (g < 0.0) lo = x;
        double next = x - g / state.dtdX;
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        converged = std::abs(next - x) <= 1e-12 * std::abs(next);
        x = next;
        state = duhring(x, tPrime);
    }
    if (!converged) call.refuse("the iterations on the concentration did not converge");

    // implicit function theorem on g(X, t, P) = B(X) + t'(P) A(X) - t = 0
    return {x, 1.0 / state.dtdX, -state.dtdTPrime * refrigerant.dCelsiusdP / state.dtdX};
}

}  // namespace coolsolve::libr_ashrae
