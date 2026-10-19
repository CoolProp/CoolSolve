#include "coolsolve/evaluator.h"
#include "coolsolve/nh3h2o.h"
#include "coolsolve/units.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <stdexcept>

// The EES ammonia-water mixture NH3H2O in the evaluator: the fluid functions (`Enthalpy(NH3H2O, T=..., P=..., X=...)`)
// and the procedure `CALL NH3H2O(Code, In1, In2, In3 : T, P, x, h, s, u, v, q)`. The thermodynamics is the module
// `nh3h2o` (IAPWS Guideline G4-01, no CoolProp); this file only converts units, assembles the derivatives and applies the
// same limits to the trial points of the solver as the CoolProp path does. Kept apart from evaluator.cpp so that the
// shared core carries two small dispatch hooks and nothing else.

namespace coolsolve {

namespace {

using nh3h2o::Prop;

/// The property designators of the EES fluid NH3H2O (`T=`, `P=`, `X=` mass fraction, `Q=` vapour mass fraction, `H=`, `S=`,
/// `U=`, `V=`), in any case, short or long.
std::optional<Prop> designatorOf(const std::string& lower) {
    if (lower == "t" || lower == "temperature") return Prop::T;
    if (lower == "p" || lower == "pressure") return Prop::P;
    if (lower == "x" || lower == "massfraction") return Prop::X;
    if (lower == "q" || lower == "quality") return Prop::Q;
    if (lower == "h" || lower == "enthalpy") return Prop::H;
    if (lower == "s" || lower == "entropy") return Prop::S;
    if (lower == "u" || lower == "internalenergy" || lower == "intenergy") return Prop::U;
    if (lower == "v" || lower == "volume") return Prop::V;
    return std::nullopt;
}

const char* propLabel(Prop p) {
    static const char* const labels[8] = {"T", "P", "X", "H", "S", "U", "V", "Q"};
    return labels[static_cast<int>(p) - 1];
}

/// The combinations of three properties of the EES procedure, as its code (the property numbers T=1, P=2, x=3, h=4, s=5,
/// u=6, v=7, q=8 in ascending order).
bool validCode(int code) {
    static const int codes[15] = {123, 128, 137, 138, 148, 158, 168, 178, 234, 235, 238, 248, 258, 268, 278};
    return std::find(std::begin(codes), std::end(codes), code) != std::end(codes);
}

std::string codeList() {
    return "123 (T, P, x), 128 (T, P, q), 137 (T, x, v), 138 (T, x, q), 148 (T, h, q), 158 (T, s, q), 168 (T, u, q), "
           "178 (T, v, q), 234 (P, x, h), 235 (P, x, s), 238 (P, x, q), 248 (P, h, q), 258 (P, s, q), 268 (P, u, q), "
           "278 (P, v, q)";
}

/// One input of the state: its property, its value in SI units (mass basis, x and q dimensionless), its factor
/// d(SI)/d(model unit) and the value of the model with its gradient with respect to the unknowns.
struct StateInput {
    Prop prop = Prop::T;
    double si = 0.0;
    double siPerUnit = 1.0;
    const ADValue* value = nullptr;
};

/// The units of the model for each property; `v` follows the density unit.
struct ModelUnits {
    const UnitSystem& units;
    double densityToSI() const { return UnitConverter::toSI(1.0, UnitType::Density, units.density); }

    /// SI value and factor d(SI)/d(unit) of the property `p` given in the unit system of the model.
    void toSI(Prop p, double v, double& si, double& factor) const {
        auto affine = [&](UnitType type, const std::string& unit) {
            si = UnitConverter::toSI(v, type, unit);
            factor = UnitConverter::scaleToSI(type, unit);
        };
        switch (p) {
            case Prop::T: affine(UnitType::Temperature, units.temperature); break;
            case Prop::P: affine(UnitType::Pressure, units.pressure); break;
            case Prop::H:
            case Prop::U: affine(UnitType::SpecificEnergy, units.specific_energy); break;
            case Prop::S: affine(UnitType::SpecificEntropy, units.specific_entropy); break;
            case Prop::V: factor = 1.0 / densityToSI(); si = v * factor; break;   // v is 1 / (a density in the density unit)
            default: si = v; factor = 1.0; break;
        }
    }
    /// Model-unit value of an SI value of the property `p` and the factor d(unit)/d(SI).
    void fromSI(Prop p, double si, double& v, double& factor) const {
        auto affine = [&](UnitType type, const std::string& unit) {
            v = UnitConverter::fromSI(si, type, unit);
            factor = UnitConverter::scaleFromSI(type, unit);
        };
        switch (p) {
            case Prop::T: affine(UnitType::Temperature, units.temperature); break;
            case Prop::P: affine(UnitType::Pressure, units.pressure); break;
            case Prop::H:
            case Prop::U: affine(UnitType::SpecificEnergy, units.specific_energy); break;
            case Prop::S: affine(UnitType::SpecificEntropy, units.specific_entropy); break;
            case Prop::V: factor = densityToSI(); v = si * factor; break;
            default: v = si; factor = 1.0; break;
        }
    }
};

/// The limits of a trial point of the solver (an input that depends on the unknowns of its block), as for CoolProp:
/// below T = 150 K (the formulation starts at the triple points, ~195 K) and p = 1000 Pa the value is limited, and x and
/// q are kept in [0, 1]. Returns the clamped SI value; `limited` tells whether it was.
double clampTrialPoint(Prop prop, double si, bool& limited) {
    double bound = si;
    switch (prop) {
        case Prop::T: bound = std::max(si, 150.0); break;
        case Prop::P: bound = std::max(si, 1000.0); break;
        case Prop::X:
        case Prop::Q: bound = std::min(std::max(si, 0.0), 1.0); break;
        default: break;
    }
    limited = bound != si;
    return bound;
}

}  // namespace

/// The flash of the three inputs (sorted by property number), with the clamping of the trial points; the inputs carry
/// their gradients.
static nh3h2o::FlashResult solveState(StateInput in[3], const std::string& label, const bool clampingEnabled,
                                      InputClampLog* clampLog) {
    std::sort(in, in + 3, [](const StateInput& a, const StateInput& b) { return static_cast<int>(a.prop) < static_cast<int>(b.prop); });
    if (clampingEnabled) {
        for (int j = 0; j < 3; ++j) {
            const auto& grad = in[j].value->gradient;
            if (std::none_of(grad.begin(), grad.end(), [](double g) { return g != 0.0; })) continue;   // a value of the model
            bool limited = false;
            const double bound = clampTrialPoint(in[j].prop, in[j].si, limited);
            if (!limited) continue;
            if (clampLog) {
                clampLog->fired = true;
                if (clampLog->report) {
                    std::ostringstream oss;
                    oss << label << ": the input " << propLabel(in[j].prop) << " = " << in[j].si << " (SI) is outside its limits ("
                        << bound << "); an input that depends on the unknowns of its block is limited there during the "
                        << "iterations, so the solution would not satisfy this equation";
                    clampLog->report->push_back(oss.str());
                }
            }
            in[j].si = bound;
        }
    }
    try {
        return nh3h2o::flash(in[0].prop, in[1].prop, in[2].prop, in[0].si, in[1].si, in[2].si);
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(label + ": " + e.what());
    }
}

/// d(output)/d(unknowns) = outPerSI * sum_j (d output / d input_j)(SI) * d(input_j SI)/d(unit) * d(input_j)/d(unknowns).
static ADValue assemble(double valueInModelUnits, const double dOutDIn[3], const StateInput in[3], double outPerSI,
                        size_t numVariables) {
    ADValue out(valueInModelUnits, numVariables);
    for (size_t i = 0; i < numVariables; ++i) {
        double g = 0.0;
        for (int j = 0; j < 3; ++j) g += dOutDIn[j] * in[j].siPerUnit * in[j].value->gradient[i];
        out.gradient[i] = outPerSI * g;
    }
    return out;
}

ADValue ExpressionEvaluator::evaluateAmmoniaWaterFunction(const FunctionCall& func, const std::string& funcName,
                                                          std::map<std::string, ADValue>& inputs) {
    const UnitSystem& units = coolpropConfig_.units;
    const ModelUnits model{units};
    const std::string label = func.name + "(NH3H2O)";

    // the output asked
    enum class Output { Prop, Density, SpecificHeat, Expansion };
    Output output = Output::Prop;
    Prop outputProp = Prop::T;
    if (funcName == "enthalpy" || funcName == "h") outputProp = Prop::H;
    else if (funcName == "entropy" || funcName == "s") outputProp = Prop::S;
    else if (funcName == "intenergy" || funcName == "internalenergy" || funcName == "u") outputProp = Prop::U;
    else if (funcName == "volume" || funcName == "v") outputProp = Prop::V;
    else if (funcName == "temperature" || funcName == "t") outputProp = Prop::T;
    else if (funcName == "pressure" || funcName == "p") outputProp = Prop::P;
    else if (funcName == "quality" || funcName == "q") outputProp = Prop::Q;
    else if (funcName == "massfraction" || funcName == "x") outputProp = Prop::X;
    else if (funcName == "density" || funcName == "rho" || funcName == "d") output = Output::Density;
    else if (funcName == "specheat" || funcName == "cp" || funcName == "c_p" || funcName == "c") output = Output::SpecificHeat;
    else if (funcName == "volexpcoef") output = Output::Expansion;
    else if (funcName == "conductivity" || funcName == "viscosity" || funcName == "surfacetension" || funcName == "k" ||
             funcName == "mu")
        throw std::runtime_error(label + ": the transport properties of NH3H2O (conductivity, viscosity, surface tension) are "
                                 "not implemented (register CS-GAP-NH3H2O, remainder); the thermodynamic properties are "
                                 "Enthalpy, Entropy, IntEnergy, Volume, Density, Temperature, Pressure, Quality, "
                                 "MassFraction, SpecHeat and VolExpCoef");
    else
        throw std::runtime_error(label + ": the function " + func.name + " is not available for NH3H2O; the properties are "
                                 "Enthalpy, Entropy, IntEnergy, Volume, Density, Temperature, Pressure, Quality, "
                                 "MassFraction, SpecHeat and VolExpCoef");

    // the three inputs
    if (inputs.size() != 3)
        throw std::runtime_error(label + ": NH3H2O needs exactly three state properties among T, P, X (mass fraction), Q, H, S, U, V "
                                 "(for example T=..., P=..., X=...), got " + std::to_string(inputs.size()));
    StateInput in[3];
    int n = 0;
    for (const auto& [name, ad] : inputs) {
        const auto prop = designatorOf(name);
        if (!prop)
            throw std::runtime_error(label + ": '" + name + "' is not a property designator of NH3H2O (T, P, X mass fraction, Q, H, S, U, V)");
        for (int j = 0; j < n; ++j)
            if (in[j].prop == *prop) throw std::runtime_error(label + ": the property " + propLabel(*prop) + " is given twice");
        in[n].prop = *prop;
        in[n].value = &ad;
        model.toSI(*prop, ad.value, in[n].si, in[n].siPerUnit);
        ++n;
    }
    {   // the combination must be one of the 15 of the EES procedure
        int sorted[3] = {static_cast<int>(in[0].prop), static_cast<int>(in[1].prop), static_cast<int>(in[2].prop)};
        std::sort(sorted, sorted + 3);
        if (!validCode(100 * sorted[0] + 10 * sorted[1] + sorted[2]))
            throw std::runtime_error(label + ": " + propLabel(static_cast<Prop>(sorted[0])) + ", " + propLabel(static_cast<Prop>(sorted[1])) +
                                     " and " + propLabel(static_cast<Prop>(sorted[2])) + " is not one of the combinations that fix an NH3H2O "
                                     "state: " + codeList());
    }
    const size_t numVars = numVariables_;
    const nh3h2o::FlashResult state = solveState(in, label, !disableClamping_, clampLog_);

    auto derivativeRow = [&](int k, double d[3]) { for (int j = 0; j < 3; ++j) d[j] = state.derivative[k][j]; };
    double d[3];
    if (output == Output::Prop) {
        const int k = static_cast<int>(outputProp) - 1;
        derivativeRow(k, d);
        double value, outPerSI;
        model.fromSI(outputProp, state.value[k], value, outPerSI);
        return assemble(value, d, in, outPerSI, numVars);
    }
    if (output == Output::Density) {   // 1 / v
        const double v = state.value[6];
        derivativeRow(6, d);
        for (int j = 0; j < 3; ++j) d[j] *= -1.0 / (v * v);
        double rho, outPerSI;
        const double rhoSI = 1.0 / v;
        rho = UnitConverter::fromSI(rhoSI, UnitType::Density, units.density);
        outPerSI = UnitConverter::scaleFromSI(UnitType::Density, units.density);
        return assemble(rho, d, in, outPerSI, numVars);
    }
    // specific heat and expansion coefficient: single-phase states only
    nh3h2o::SinglePhaseExtras extras;
    try {
        extras = nh3h2o::singlePhaseExtras(in[0].prop, in[1].prop, in[2].prop, in[0].si, in[1].si, in[2].si);
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(label + ": " + e.what());
    }
    if (output == Output::SpecificHeat) {
        const double value = UnitConverter::fromSI(extras.cp, UnitType::SpecificHeat, units.specific_heat);
        const double outPerSI = UnitConverter::scaleFromSI(UnitType::SpecificHeat, units.specific_heat);
        return assemble(value, extras.dcp, in, outPerSI, numVars);
    }
    // (1/v)(dv/dT)_p per unit of the temperature unit of the model
    const double perUnit = 1.0 / UnitConverter::scaleFromSI(UnitType::Temperature, units.temperature);
    return assemble(extras.expansion * perUnit, extras.dexpansion, in, perUnit, numVars);
}

void ExpressionEvaluator::evaluateAmmoniaWaterCall(const ProcedureCall& call) {
    if (call.inputArgs.size() != 4 || call.outputVars.empty() || call.outputVars.size() > 8)
        throw std::runtime_error("CALL NH3H2O(Code, In1, In2, In3 : T, P, x, h, s, u, v, q) takes 4 inputs and up to 8 outputs, got " +
                                 std::to_string(call.inputArgs.size()) + " input(s) and " + std::to_string(call.outputVars.size()) +
                                 " output(s)");
    const ADValue codeValue = evaluate(call.inputArgs[0]);
    const long code = std::lround(codeValue.value);
    if (std::abs(codeValue.value - static_cast<double>(code)) > 1e-9 || !validCode(static_cast<int>(code)))
        throw std::runtime_error("CALL NH3H2O: the code " + std::to_string(codeValue.value) + " is not one of " + codeList());
    const int digits[3] = {static_cast<int>(code / 100), static_cast<int>(code / 10 % 10), static_cast<int>(code % 10)};

    // the procedure has fixed units: T [K], P [bar], x [kg/kg], h, u [kJ/kg], s [kJ/(kg K)], v [m3/kg], q
    auto procedureUnit = [](Prop p) {   // d(SI)/d(unit) of the input of the procedure
        switch (p) {
            case Prop::P: return 1.0e5;
            case Prop::H:
            case Prop::S:
            case Prop::U: return 1.0e3;
            default: return 1.0;
        }
    };
    ADValue values[3];
    StateInput in[3];
    for (int j = 0; j < 3; ++j) {
        values[j] = evaluate(call.inputArgs[1 + j]);
        in[j].prop = static_cast<Prop>(digits[j]);
        in[j].siPerUnit = procedureUnit(in[j].prop);
        in[j].si = values[j].value * in[j].siPerUnit;
        in[j].value = &values[j];
    }
    const std::string label = "CALL NH3H2O(" + std::to_string(code) + ")";
    const size_t numVars = numVariables_;
    const nh3h2o::FlashResult state = solveState(in, label, !disableClamping_, clampLog_);
    for (size_t k = 0; k < call.outputVars.size(); ++k) {
        double d[3];
        for (int j = 0; j < 3; ++j) d[j] = state.derivative[k][j];
        const Prop prop = static_cast<Prop>(k + 1);
        const double perSI = 1.0 / procedureUnit(prop);
        setVariable(resolveVariableName(call.outputVars[k]), assemble(state.value[k] * perSI, d, in, perSI, numVars));
    }
}

}  // namespace coolsolve
