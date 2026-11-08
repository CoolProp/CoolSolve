#include "coolsolve/ideal_gas.h"
#include "AbstractState.h"
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace coolsolve {

namespace {

using CoolProp::parameters;

/// A CoolProp fluid as an ideal gas: its own HEOS state, in the gas phase for good, and its gas constant.
struct Gas {
    std::string name;
    std::unique_ptr<CoolProp::AbstractState> state;
    double R = 0.0;       ///< molar gas constant of the equation of state [J/mol-K]
    double Rs = 0.0;      ///< specific gas constant [J/kg-K]
    double M = 0.0;       ///< molar mass [kg/mol]
    double lastT = 300.0; ///< the temperature of the last solve: where the next one starts
};

/// The state of a fluid, created once per thread. specify_phase(gas) makes every update(DmolarT_INPUTS) a plain
/// evaluation of the Helmholtz terms at (T, rho): no phase determination, no limit of the equation of state.
Gas& gasOf(const std::string& coolPropName) {
    static thread_local std::unordered_map<std::string, Gas> cache;
    auto it = cache.find(coolPropName);
    if (it != cache.end()) return it->second;
    Gas gas;
    gas.name = coolPropName;
    gas.state.reset(CoolProp::AbstractState::factory("HEOS", coolPropName));
    gas.state->specify_phase(CoolProp::iphase_gas);
    gas.R = gas.state->gas_constant();
    gas.M = gas.state->molar_mass();
    gas.Rs = gas.R / gas.M;
    return cache.emplace(coolPropName, std::move(gas)).first->second;
}

std::string numberText(double v) {
    std::ostringstream oss;
    oss.precision(10);
    oss << v;
    return oss.str();
}

std::runtime_error idealError(const Gas& gas, const std::string& what) {
    return std::runtime_error("ideal gas " + gas.name + ": " + what);
}

/// cp0, h0, s0 (and dcp0/dT when asked) at T, from the derivatives of alpha0 with respect to tau = Tr/T at the
/// density of the standard pressure: h = R T (1 + tau a1), s = R (tau a1 - a0), cp = R (1 - tau^2 a2).
IdealGasPoint evaluate(Gas& gas, double T, bool withDcp) {
    if (!(T > 0.0) || !std::isfinite(T))
        throw idealError(gas, "the temperature " + numberText(T) + " K is not positive");
    double tau, a0, a1, a2, a3 = 0.0;
    try {
        gas.state->update(CoolProp::DmolarT_INPUTS, kIdealGasStandardPressureSI / (gas.R * T), T);
        tau = gas.state->tau();
        a0 = gas.state->alpha0();
        a1 = gas.state->dalpha0_dTau();
        a2 = gas.state->d2alpha0_dTau2();
        if (withDcp) a3 = gas.state->d3alpha0_dTau3();
    } catch (const std::exception& e) {
        throw idealError(gas, "the ideal-gas part of its CoolProp equation of state is not defined at T = " +
                                  numberText(T) + " K (" + e.what() + ")");
    }
    IdealGasPoint p;
    p.R = gas.Rs;
    p.M = gas.M;
    p.h = gas.Rs * T * (1.0 + tau * a1);
    p.s = gas.Rs * (tau * a1 - a0);
    p.cp = gas.Rs * (1.0 - tau * tau * a2);
    p.dcp_dT = gas.Rs * (tau / T) * (2.0 * tau * a2 + tau * tau * a3);
    if (!std::isfinite(p.h + p.s + p.cp + p.dcp_dT))
        throw idealError(gas, "the ideal-gas part of its CoolProp equation of state is not finite at T = " +
                                  numberText(T) + " K");
    return p;
}

/// The dilute-gas transport property (viscosity or conductivity) at T: CoolProp's model on the imposed-gas state at
/// 1 Pa, where the density terms vanish (they are below 1e-5 of the dilute term for T >= 200 K).
double transport(Gas& gas, parameters which, double T) {
    if (!(T > 0.0) || !std::isfinite(T))
        throw idealError(gas, "the temperature " + numberText(T) + " K is not positive");
    try {
        gas.state->update(CoolProp::DmolarT_INPUTS, 1.0 / (gas.R * T), T);
        return which == CoolProp::iviscosity ? gas.state->viscosity() : gas.state->conductivity();
    } catch (const std::exception& e) {
        throw idealError(gas, e.what());
    }
}

bool fixesTemperature(parameters p) { return p == CoolProp::iT || p == CoolProp::iHmass || p == CoolProp::iUmass; }

/**
 * Safeguarded Newton iteration for the temperature T of the increasing function f(T) = target, `eval(T, f, df)`.
 * The bracket starts at the whole range; a Newton step that leaves the bracket is replaced by its middle. The
 * iteration starts at the temperature of the last solve. An answer at a bound of the range is an error (the target
 * is outside what the range holds).
 */
template <class Eval>
double solveTemperature(Gas& gas, Eval&& eval, double target, const std::string& quantity) {
    double lo = kIdealGasSolveTMinSI, hi = kIdealGasSolveTMaxSI;
    double T = std::min(std::max(gas.lastT, lo), hi);
    for (int iteration = 0; iteration < 200; ++iteration) {
        double f, df;
        eval(T, f, df);
        const double residual = f - target;
        if (residual > 0.0) hi = std::min(hi, T); else lo = std::max(lo, T);
        double next = T - residual / df;
        if (!(df > 0.0) || !(next > lo && next < hi)) next = 0.5 * (lo + hi);
        if (std::abs(next - T) <= 1e-13 * T || hi - lo <= 1e-13 * hi) {
            T = next;
            if (T - kIdealGasSolveTMinSI < 1e-6 * T || kIdealGasSolveTMaxSI - T < 1e-6 * T) break;
            gas.lastT = T;
            return T;
        }
        T = next;
    }
    throw idealError(gas, "the " + quantity + " " + numberText(target) + " is not reached between " +
                              numberText(kIdealGasSolveTMinSI) + " K and " + numberText(kIdealGasSolveTMaxSI) + " K");
}

/// The temperature T and the pressure P of the state fixed by two inputs.
void solveState(Gas& gas, parameters a, double va, parameters b, double vb, double& T, double& P) {
    const double R = gas.Rs;
    if (fixesTemperature(a) && fixesTemperature(b))
        throw idealError(gas, "two inputs that each fix the temperature (T, h, u) do not fix the state");
    if (!fixesTemperature(a) && fixesTemperature(b)) { std::swap(a, b); std::swap(va, vb); }

    if (fixesTemperature(a)) {
        if (a == CoolProp::iT) T = va;
        else if (a == CoolProp::iHmass)
            T = solveTemperature(gas, [&](double x, double& f, double& d) {
                const IdealGasPoint p = evaluate(gas, x, false); f = p.h; d = p.cp; }, va, "enthalpy");
        else  // internal energy
            T = solveTemperature(gas, [&](double x, double& f, double& d) {
                const IdealGasPoint p = evaluate(gas, x, false); f = p.h - R * x; d = p.cp - R; }, va, "internal energy");
        if (b == CoolProp::iP) P = vb;
        else if (b == CoolProp::iDmass) P = vb * R * T;
        else if (b == CoolProp::iSmass) P = kIdealGasStandardPressureSI * std::exp((evaluate(gas, T, false).s - vb) / R);
        else throw idealError(gas, "an input other than T, P, h, u, d (or v) and s is not a state of an ideal gas");
        return;
    }
    // neither input fixes the temperature: the pairs (P, d), (P, s), (d, s)
    if (a == CoolProp::iSmass || (a == CoolProp::iDmass && b == CoolProp::iP)) { std::swap(a, b); std::swap(va, vb); }
    if (a == CoolProp::iP && b == CoolProp::iDmass) { P = va; T = P / (vb * R); return; }
    if (a == CoolProp::iP && b == CoolProp::iSmass) {
        P = va;
        if (!(P > 0.0)) throw idealError(gas, "the pressure " + numberText(P) + " Pa is not positive");
        const double target = vb + R * std::log(P / kIdealGasStandardPressureSI);
        T = solveTemperature(gas, [&](double x, double& f, double& d) {
            const IdealGasPoint p = evaluate(gas, x, false); f = p.s; d = p.cp / x; }, target, "entropy");
        return;
    }
    if (a == CoolProp::iDmass && b == CoolProp::iSmass) {
        const double rho = va;
        if (!(rho > 0.0)) throw idealError(gas, "the density " + numberText(rho) + " kg/m3 is not positive");
        T = solveTemperature(gas, [&](double x, double& f, double& d) {
            const IdealGasPoint p = evaluate(gas, x, false);
            f = p.s - R * std::log(rho * R * x / kIdealGasStandardPressureSI); d = (p.cp - R) / x; }, vb, "entropy");
        P = rho * R * T;
        return;
    }
    throw idealError(gas, "the two inputs are the same property, or one is not a state variable of an ideal gas");
}

/// An input as a function of the state: its value and its derivatives with respect to T and P.
void inputRow(parameters in, double T, double P, const IdealGasPoint& p, double& v, double& dT, double& dP) {
    switch (in) {
        case CoolProp::iT: v = T; dT = 1.0; dP = 0.0; return;
        case CoolProp::iP: v = P; dT = 0.0; dP = 1.0; return;
        case CoolProp::iHmass: v = p.h; dT = p.cp; dP = 0.0; return;
        case CoolProp::iUmass: v = p.h - p.R * T; dT = p.cp - p.R; dP = 0.0; return;
        case CoolProp::iDmass: v = P / (p.R * T); dT = -v / T; dP = v / P; return;
        case CoolProp::iSmass: v = p.s - p.R * std::log(P / kIdealGasStandardPressureSI); dT = p.cp / T; dP = -p.R / P; return;
        default: break;
    }
    throw std::runtime_error("an input that is not T, P, h, u, d (or v) and s is not a state of an ideal gas");
}

}  // namespace

IdealGasPoint idealGasAt(const std::string& coolPropName, double T) {
    return evaluate(gasOf(coolPropName), T, true);
}

IdealGasValue idealGasProperty(const std::string& coolPropName, parameters output,
                               parameters input1, double value1, parameters input2, double value2, bool derivatives) {
    Gas& gas = gasOf(coolPropName);
    const double R = gas.Rs;
    double T = 0.0, P = 0.0;
    solveState(gas, input1, value1, input2, value2, T, P);
    if (!(P > 0.0) || !std::isfinite(P))
        throw idealError(gas, "the pressure " + numberText(P) + " Pa of the state is not positive");

    // dcp0/dT is always evaluated (it costs a third derivative of alpha0): the validity of the ideal-gas part of the
    // equation of state is the same for the values and for the derivatives, in the solve and in the check
    const IdealGasPoint p = evaluate(gas, T, true);

    // the property and its derivatives with respect to T and P
    double v = 0.0, dT = 0.0, dP = 0.0;
    switch (output) {
        case CoolProp::iT: v = T; dT = 1.0; break;
        case CoolProp::iP: v = P; dP = 1.0; break;
        case CoolProp::iHmass: v = p.h; dT = p.cp; break;
        case CoolProp::iUmass: v = p.h - R * T; dT = p.cp - R; break;
        case CoolProp::iSmass: case CoolProp::iDmass:
            inputRow(output, T, P, p, v, dT, dP); break;
        case CoolProp::iCpmass: v = p.cp; dT = p.dcp_dT; break;
        case CoolProp::iCvmass: v = p.cp - R; dT = p.dcp_dT; break;
        case CoolProp::iisentropic_expansion_coefficient: case CoolProp::ispeed_sound: {
            const double cv = p.cp - R, gamma = p.cp / cv;
            const double dgamma = p.dcp_dT * (1.0 / cv - p.cp / (cv * cv));
            if (output == CoolProp::iisentropic_expansion_coefficient) { v = gamma; dT = dgamma; }
            else { v = std::sqrt(gamma * R * T); dT = 0.5 * v * (dgamma / gamma + 1.0 / T); }
            break;
        }
        case CoolProp::iZ: v = 1.0; break;
        case CoolProp::imolar_mass: v = gas.M; break;
        case CoolProp::iviscosity: case CoolProp::iconductivity: case CoolProp::iPrandtl: {
            // dilute gas: a function of T alone
            auto transportAt = [&](double x, double& mu, double& k) {
                mu = transport(gas, CoolProp::iviscosity, x);
                k = transport(gas, CoolProp::iconductivity, x);
            };
            auto property = [&](double x, double cp) {
                if (output == CoolProp::iviscosity) return transport(gas, CoolProp::iviscosity, x);
                if (output == CoolProp::iconductivity) return transport(gas, CoolProp::iconductivity, x);
                double mu, k; transportAt(x, mu, k); return cp * mu / k;
            };
            v = property(T, p.cp);
            if (derivatives) {   // central difference in T: CoolProp has no derivative of its transport models
                const double h = 1e-5 * T;
                const IdealGasPoint up = evaluate(gas, T + h, false), down = evaluate(gas, T - h, false);
                dT = (property(T + h, up.cp) - property(T - h, down.cp)) / (2.0 * h);
            }
            break;
        }
        default:
            throw idealError(gas, "the property has no ideal-gas value (the output is not a property of the state of an ideal gas)");
    }

    IdealGasValue result;
    result.value = v;
    if (!derivatives) return result;

    // chain rule through the inputs: d(out)/d(in) = grad_(T,P)(out) * J^-1 with J = d(in1, in2)/d(T, P)
    double f1, f1T, f1P, f2, f2T, f2P;
    inputRow(input1, T, P, p, f1, f1T, f1P);
    inputRow(input2, T, P, p, f2, f2T, f2P);
    const double det = f1T * f2P - f1P * f2T;
    if (det == 0.0 || !std::isfinite(det)) throw idealError(gas, "the two inputs do not fix the state (singular)");
    result.d1 = (dT * f2P - dP * f2T) / det;
    result.d2 = (dP * f1T - dT * f1P) / det;
    return result;
}

}  // namespace coolsolve
