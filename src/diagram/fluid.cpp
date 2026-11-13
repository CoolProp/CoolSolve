// The fluid of a chart (include/coolsolve/diagram/fluid.h).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/native.py
// (fluid_string, port_fluid_string), spec/fluid.py (FluidSpec and HumidAirSpec records), thermo/provider.py (the
// provider interface) and thermo/fluid.py (Fluid).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#include "coolsolve/diagram/fluid.h"

#include "coolsolve/diagram/humid_air.h"
#include "coolsolve/fluids.h"

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/DataStructures.h"

#include <coolplot/core.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace coolsolve {
namespace diagram {

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

/// The CoolProp parameter of a CoolPlot quantity key of a fluid state (thermo.fluid.parameter_index); "v" is
/// handled by the callers, as a density.
CoolProp::parameters parameterIndex(const std::string& key) {
    if (key == "T") return CoolProp::iT;
    if (key == "p") return CoolProp::iP;
    if (key == "h") return CoolProp::iHmass;
    if (key == "s") return CoolProp::iSmass;
    if (key == "rho") return CoolProp::iDmass;
    if (key == "u") return CoolProp::iUmass;
    if (key == "Q") return CoolProp::iQ;
    throw DiagramError("'" + key + "' is not a quantity of a fluid state.");
}

/// A specific volume given as an input becomes the density CoolProp works with (thermo.fluid._as_density).
std::pair<std::string, double> asDensity(const std::string& key, double value) {
    if (key != "v") return {key, value};
    if (!(value > 0.0)) throw std::runtime_error("A specific volume has to be positive, not " + coolplot::general_text(value, 6) + ".");
    return {"rho", 1.0 / value};
}

/// Whether the backend is a tabular one built on an equation of state ("BICUBIC&HEOS", "TTSE&HEOS").
bool isTabular(const std::string& backend) {
    return backend.find('&') != std::string::npos;
}

}  // namespace

int coolpropParameter(const std::string& key) {
    return parameterIndex(key);
}

const std::vector<std::string>& stateKeys() {
    static const std::vector<std::string> keys = {"T", "p", "h", "s", "rho", "v", "u", "Q"};
    return keys;
}

// ---------------------------------------------------------------------------
// ChartFluid (spec/fluid.py, thermo/native.py)
// ---------------------------------------------------------------------------

std::string ChartFluid::name() const {
    if (isHumidAir()) return humidAirName(p_Pa);
    if (isIdealGas()) return idealGas;
    std::string names;
    for (size_t i = 0; i < components.size(); ++i) {
        if (i) names += "&";
        names += components[i];
        if (!fractions.empty()) names += "[" + coolplot::general_text(fractions[i], 12) + "]";
    }
    return names;
}

std::string ChartFluid::fluidString() const {
    if (isHumidAir()) return "HUMIDAIR::" + coolplot::general_text(p_Pa, 17);
    return (isIdealGas() ? std::string("IDEALGAS") : backend) + "::" + name();
}

std::string ChartFluid::portFluidString(const std::string& key) const {
    if (key != "Q" || !isTabular(backend)) return fluidString();
    ChartFluid base = *this;
    base.backend = backend.substr(backend.find('&') + 1);
    return base.fluidString();
}

std::string ChartFluid::kind() const {
    if (isHumidAir()) return "HumidAir";
    if (isIdealGas()) return "IdealGas";
    if (backend.rfind("INCOMP", 0) == 0) return "Incompressible";
    return components.size() > 1 ? "Mixture" : "Real";
}

json ChartFluid::record() const {
    if (isHumidAir()) return {{"model", "humid_air"}, {"p_Pa", p_Pa}};
    json out = {{"backend", backend}, {"components", components}, {"fractions", fractions}, {"basis", basis},
                {"interactions", json::array()}};
    if (isIdealGas()) out["coolsolve"] = {{"model", "ideal_gas"}, {"fluid", idealGas}};
    return out;
}

ChartFluid ChartFluid::fromRecord(const json& record) {
    if (record.is_object() && record.value("model", "") == "humid_air") {  // HumidAirSpec.from_dict
        const json p = record.value("p_Pa", json(DEFAULT_HUMID_AIR_P_PA));
        if (!p.is_number()) throw DiagramError("The total pressure of humid air (\"p_Pa\") is a number of Pa, not " + p.dump() + ".");
        ChartFluid fluid;
        fluid.backend.clear();
        fluid.basis.clear();
        fluid.p_Pa = checkedHumidAirPressure(p.get<double>());
        return fluid;
    }
    if (!record.is_object() || !record.contains("components"))
        throw DiagramError("A fluid record needs its components, as CoolPlot writes it (FluidSpec.to_dict).");
    ChartFluid fluid;
    fluid.backend = record.value("backend", "HEOS");
    if (fluid.backend.empty()) fluid.backend = "HEOS";
    fluid.components = record.at("components").get<std::vector<std::string>>();
    fluid.fractions = record.value("fractions", std::vector<double>());
    fluid.basis = record.value("basis", "mole");
    if (fluid.components.empty()) throw DiagramError("A fluid record needs at least one component.");
    if (!fluid.fractions.empty() && fluid.fractions.size() != fluid.components.size())
        throw DiagramError("A fluid record needs one fraction per component, or none.");
    // The basis the port reads (native.fluid_string): mass fractions for INCOMP, mole fractions otherwise
    const std::string expected = fluid.backend.rfind("INCOMP", 0) == 0 ? "mass" : "mole";
    if (!fluid.fractions.empty() && fluid.basis != expected)
        throw DiagramError("The fractions of a " + fluid.backend + " fluid are drawable as " + expected + " fractions only, not as '" +
                           fluid.basis + "'.");
    if (record.contains("interactions") && !record.at("interactions").empty())
        throw DiagramError("Fluids with changed binary interaction parameters are not drawable.");
    const json coolsolve = record.value("coolsolve", json::object());
    if (coolsolve.is_object() && coolsolve.value("model", "") == "ideal_gas") fluid.idealGas = coolsolve.value("fluid", "");
    return fluid;
}

std::string diagramFluidName(const Fluid& fluid) {
    return fluid.getType() == FluidType::HumidAir ? HUMID_AIR_FLUID : fluid.getName();
}

std::pair<double, double> solutionFractionRange(const std::string& coolpropKey) {
    try {
        std::unique_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory("INCOMP", coolpropKey));
        return {state->trivial_keyed_output(CoolProp::ifraction_min), state->trivial_keyed_output(CoolProp::ifraction_max)};
    } catch (const std::exception& e) {
        throw DiagramError("CoolProp does not know the incompressible fluid '" + coolpropKey + "': " + e.what());
    }
}

ChartFluid chartFluid(const std::string& name, const std::string& backend, const std::vector<double>& fractions,
                      const std::string& basis, std::optional<double> p_Pa) {
    std::string error;
    const std::shared_ptr<Fluid> fluid = FluidRegistry::getFluid(name, &error);
    if (!fluid) throw DiagramError("Unknown fluid '" + name + "'" + (error.empty() ? "" : ": " + error) + ".");
    if (p_Pa && fluid->getType() != FluidType::HumidAir)
        throw DiagramError("'" + name + "' takes no total pressure (\"p_Pa\"): it is the pressure of a chart of humid air.");
    const auto solution = std::dynamic_pointer_cast<IncompressibleFluid>(fluid);
    if ((!fractions.empty() || !basis.empty()) && !(solution && solution->isSolution()))
        throw DiagramError("'" + name + "' takes no fractions: they are the concentration of an incompressible solution (a "
                           "mixture carries its mole fractions in its name).");
    ChartFluid chart;
    chart.backend = backend.empty() ? "HEOS" : backend;
    auto notYet = [&name](const std::string& kind) {
        return DiagramError("'" + name + "' is " + kind + ": its charts are not drawable yet (they come with a later "
                            "version of the Diagrams tab).");
    };
    switch (fluid->getType()) {
        case FluidType::Real:
            chart.components = {fluid->getCoolPropName()};
            return chart;
        case FluidType::Mixture: {
            const auto mixture = std::dynamic_pointer_cast<MixtureFluid>(fluid);
            if (!mixture) break;
            chart.components = mixture->getComponents();
            chart.fractions = mixture->getMoleFractions();
            return chart;
        }
        case FluidType::IdealGas:
            chart.components = {fluid->getCoolPropName()};
            chart.idealGas = fluid->getName();
            return chart;
        case FluidType::Incompressible: {
            if (!solution) break;
            // CoolProp's INCOMP backend whatever the model's backend, at a mass fraction (spec.fluid.FluidSpec.parse)
            chart.backend = "INCOMP";
            chart.components = {solution->getCoolPropKey()};
            if (!solution->isSolution()) return chart;
            chart.basis = "mass";
            if (!basis.empty() && basis != "mass")
                throw DiagramError("The concentration of '" + name + "' is a mass fraction (basis \"mass\"), not a " + basis + " fraction.");
            const auto [lo, hi] = solutionFractionRange(chart.components[0]);
            const std::string range = coolplot::general_text(lo, 6) + " and " + coolplot::general_text(hi, 6);
            if (fractions.size() != 1)
                throw DiagramError("'" + name + "' is a solution: give its concentration as a mass fraction between " + range +
                                   ", e.g. \"fractions\": [" + coolplot::general_text(0.5 * (lo + hi), 6) + "].");
            if (!(fractions[0] >= lo && fractions[0] <= hi))
                throw DiagramError("The concentration of '" + name + "' is a mass fraction between " + range + ", not " +
                                   coolplot::general_text(fractions[0], 6) + ".");
            chart.fractions = fractions;
            return chart;
        }
        case FluidType::HumidAir:  // HumidAirSpec: no components, the total pressure of the chart
            chart.backend.clear();
            chart.basis.clear();
            chart.p_Pa = checkedHumidAirPressure(p_Pa.value_or(DEFAULT_HUMID_AIR_P_PA));
            return chart;
        case FluidType::AmmoniaWater: throw notYet("the ammonia-water mixture of EES");
        default: break;
    }
    throw DiagramError("'" + name + "' cannot be drawn.");
}

// ---------------------------------------------------------------------------
// PropertyProvider (thermo/provider.py)
// ---------------------------------------------------------------------------

State PropertyProvider::statePoint(const std::string& key1, double value1, const std::string& key2, double value2) {
    update(key1, value1, key2, value2);
    State state;
    for (const std::string& key : quantities()) state[key] = get(key);
    return state;
}

std::string PropertyProvider::pathKind(const State& /*a*/, const State& /*b*/, const std::string& kind) const {
    return kind;
}

std::string PropertyProvider::stateProblem(const State& /*state*/) const {
    return "";
}

const State& PropertyProvider::critical() {
    throw std::runtime_error(name() + " has no critical point.");
}

std::pair<double, double> PropertyProvider::saturationRange(const std::string& key, double /*margin*/) {
    const std::string k = coolplot::quantity_key(key);
    if (k != "T" && k != "p") throw DiagramError("Saturation ranges are defined in T or p, not in '" + key + "'.");
    throw DiagramError(name() + " has no saturation curve (CoolProp provides no critical point for it); give absolute "
                                "temperature and pressure limits.");
}

// ---------------------------------------------------------------------------
// FluidModel (thermo/fluid.py)
// ---------------------------------------------------------------------------

FluidModel::FluidModel(const ChartFluid& fluid)
    : name_(fluid.name()), backend_(fluid.backend), components_(fluid.components), fractions_(fluid.fractions), basis_(fluid.basis) {
    state_ = makeState();
    current_ = state_;
}

FluidModel::~FluidModel() = default;

std::shared_ptr<CoolProp::AbstractState> FluidModel::makeState() const {
    std::string names;
    for (size_t i = 0; i < components_.size(); ++i) names += (i ? "&" : "") + components_[i];
    std::shared_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory(backend_, names));
    if (!fractions_.empty()) {
        if (basis_ == "mass") state->set_mass_fractions(fractions_);  // an incompressible solution (thermo.fluid.build_state)
        else state->set_mole_fractions(fractions_);
    }
    return state;
}

std::shared_ptr<CoolProp::AbstractState> FluidModel::guideState() {
    if (!guideTried_) {
        guideTried_ = true;
        try {
            std::shared_ptr<CoolProp::AbstractState> guide = makeState();
            guide->build_phase_envelope("veryfine");
            guide_ = guide;
        } catch (const std::exception&) {  // no guide: callers check for null
        }
    }
    return guide_;
}

void FluidModel::update(const std::string& key1, double value1, const std::string& key2, double value2) {
    const auto [k1, v1] = asDensity(coolplot::quantity_key(key1), value1);
    const auto [k2, v2] = asDensity(coolplot::quantity_key(key2), value2);
    double o1 = 0.0, o2 = 0.0;
    const CoolProp::input_pairs pair = CoolProp::generate_update_pair(parameterIndex(k1), v1, parameterIndex(k2), v2, o1, o2);
    if (pair == CoolProp::INPUT_PAIR_INVALID)
        throw std::runtime_error("CoolProp has no flash for the inputs " + k1 + " and " + k2 + ".");
    current_ = state_;
    try {
        state_->update(pair, o1, o2);
    } catch (const std::exception&) {
        // A mixture state CoolProp cannot find is tried once more with the phase envelope as a guide
        std::shared_ptr<CoolProp::AbstractState> guide = isPure() ? nullptr : guideState();
        if (!guide) throw;
        guide->update(pair, o1, o2);
        current_ = guide;
    }
    const std::string problem = inconsistency(*current_);
    if (!problem.empty()) throw std::runtime_error("CoolProp returned an inconsistent state for " + name_ + ": " + problem);
    // The check used the guide state for its saturation lookups: put the requested state back on it
    if (current_ != state_) current_->update(pair, o1, o2);
}

std::string FluidModel::inconsistency(CoolProp::AbstractState& state) {
    if (isPure()) return "";
    double T = NaN, rho = NaN;
    try {
        const double Q = state.Q();
        if (Q >= 0.0 && Q <= 1.0) return "";  // two-phase states come out of the saturation solver itself
        T = state.T();
        rho = state.rhomolar();
    } catch (const std::exception&) {
        return "";
    }
    if (!hasCriticalPoint() || T >= critical().at("T") * (1.0 - 1e-6)) return "";
    std::shared_ptr<CoolProp::AbstractState> guide = guideState();
    if (!guide) return "";
    double rhoVap = NaN, rhoLiq = NaN;
    try {
        guide->update(CoolProp::QT_INPUTS, 1.0, T);
        rhoVap = guide->rhomolar();
        guide->update(CoolProp::QT_INPUTS, 0.0, T);
        rhoLiq = guide->rhomolar();
    } catch (const std::exception&) {
        return "";  // no saturation state at this temperature: nothing to compare with
    }
    if (rhoVap * (1.0 + 1e-6) < rho && rho < rhoLiq * (1.0 - 1e-6)) {
        std::ostringstream text;
        text << "reported as single phase at T = " << T << " K with rho = " << rho
             << " mol/m3, between the saturated vapour and liquid densities, i.e. inside the two-phase region";
        return text.str();
    }
    return "";
}

double FluidModel::get(const std::string& key) const {
    const std::string k = coolplot::quantity_key(key);
    if (k == "v") return 1.0 / current_->keyed_output(CoolProp::iDmass);
    return current_->keyed_output(parameterIndex(k));
}

double FluidModel::trivial(int key) const {
    try {
        const double value = state_->trivial_keyed_output(static_cast<CoolProp::parameters>(key));
        return std::isfinite(value) ? value : NaN;
    } catch (const std::exception&) {
        return NaN;
    }
}

double FluidModel::lowestTemperature(int key) {
    // For a mixture the highest value of its components: CoolProp reports the mole-weighted average, below
    // where the equation of state of some component is valid
    double value = trivial(key);
    if (!isPure()) {
        for (const std::string& component : components_) {
            try {
                std::unique_ptr<CoolProp::AbstractState> one(CoolProp::AbstractState::factory(backend_, component));
                const double v = one->trivial_keyed_output(static_cast<CoolProp::parameters>(key));
                if (std::isfinite(v) && !(v <= value)) value = v;
            } catch (const std::exception&) {  // a component without this limit
            }
        }
    }
    return value;
}

double FluidModel::tTripleK() { return lowestTemperature(CoolProp::iT_triple); }

double FluidModel::tMinK() {
    const double T_min = lowestTemperature(CoolProp::iT_min);
    // CoolSolve's: an incompressible solution has no state below its freezing point (CoolProp refuses them). A margin
    // of 1e-9 keeps a line at the lowest temperature above it once its value is snapped to 12 digits (coolplot::snap)
    const double T_freeze = backend_.rfind("INCOMP", 0) == 0 ? trivial(CoolProp::iT_freeze) * (1.0 + 1e-9) : NaN;
    return std::isfinite(T_freeze) && !(T_freeze <= T_min) ? T_freeze : T_min;
}

double FluidModel::tMaxK() { return trivial(CoolProp::iT_max); }
double FluidModel::pMinPa() { return trivial(CoolProp::iP_min); }
double FluidModel::pMaxPa() { return trivial(CoolProp::iP_max); }

bool FluidModel::hasCriticalPoint() {
    if (!hasCritical_) {
        try {
            critical();
            hasCritical_ = true;
        } catch (const std::exception&) {
            hasCritical_ = false;
        }
    }
    return *hasCritical_;
}

const State& FluidModel::critical() {
    if (critical_) return *critical_;
    // thermo.fluid._critical_values: CoolProp's values, or the stable point of highest temperature of a mixture
    double T = NaN, p = NaN, rhomolar = NaN;
    try {
        T = state_->T_critical();
        p = state_->p_critical();
        rhomolar = state_->rhomolar_critical();
    } catch (const std::exception&) {
        try {
            for (const CoolProp::CriticalState& point : state_->all_critical_points())
                if (point.stable && !(point.T <= T)) {
                    T = point.T;
                    p = point.p;
                    rhomolar = point.rhomolar;
                }
        } catch (const std::exception&) {  // T stays NaN, reported below
        }
    }
    if (!std::isfinite(T)) throw std::runtime_error("CoolProp could not provide the critical point of " + name_ + ".");
    // thermo.fluid._critical_point: the other properties on a separate state, the phase imposed first
    std::shared_ptr<CoolProp::AbstractState> crit = makeState();
    std::vector<std::pair<CoolProp::input_pairs, std::pair<double, double>>> attempts;
    if (std::isfinite(rhomolar)) attempts.push_back({CoolProp::DmolarT_INPUTS, {rhomolar, T}});
    if (std::isfinite(p)) attempts.push_back({CoolProp::PT_INPUTS, {p, T}});
    for (const auto& attempt : attempts) {
        for (bool impose : {true, false}) {
            try {
                if (impose) crit->specify_phase(CoolProp::iphase_critical_point);
                else crit->unspecify_phase();
                crit->update(attempt.first, attempt.second.first, attempt.second.second);
                State values;
                for (const char* key : {"T", "p", "h", "s", "rho", "u"}) values[key] = crit->keyed_output(parameterIndex(key));
                critical_ = values;
                return *critical_;
            } catch (const std::exception&) {  // the next attempt
            }
        }
    }
    if (!std::isfinite(p)) throw std::runtime_error("Could not evaluate the critical point of " + name_ + ".");
    critical_ = State{{"T", T}, {"p", p}, {"h", NaN}, {"s", NaN}, {"rho", NaN}, {"u", NaN}};
    return *critical_;
}

std::array<double, 4> FluidModel::defaultBox() {
    return {tMinK(), tMaxK(), LIQUID_P_MIN_PA, LIQUID_P_MAX_PA};
}

double FluidModel::defaultPressure() const {
    return backend_.rfind("INCOMP", 0) == 0 ? IncompressibleFluid::kInternalPressureSI : NaN;
}

std::string FluidModel::pseudoPureBlend() {
    if (!pseudoPure_) {
        pseudoPure_ = std::string();
        if (components_.size() == 1 && state_->backend_name().find("Helmholtz") != std::string::npos) {
            try {
                if (state_->fluid_param_string("pure") == "false") {
                    const std::string blend = components_[0] + ".mix";
                    const std::string all = "," + CoolProp::get_global_param_string("predefined_mixtures") + ",";
                    if (all.find("," + blend + ",") != std::string::npos) pseudoPure_ = blend;
                }
            } catch (const std::exception&) {  // a backend without the parameter: not a pseudo-pure fluid
            }
        }
    }
    return *pseudoPure_;
}

std::pair<double, double> FluidModel::saturationRange(const std::string& key, double margin) {
    const std::string k = coolplot::quantity_key(key);
    if (k != "T" && k != "p") throw DiagramError("Saturation ranges are defined in T or p, not in '" + key + "'.");
    if (!hasCriticalPoint())
        throw DiagramError(name_ + " has no saturation curve (CoolProp provides no critical point for it); give absolute "
                                   "temperature and pressure limits.");
    const State& crit = critical();
    double T_lo = 0.0;
    for (double v : {tTripleK(), tMinK()})
        if (std::isfinite(v)) T_lo = std::max(T_lo, v);
    T_lo += crit.at("T") * margin;
    if (k == "T") return {T_lo, crit.at("T") * (1.0 - margin)};
    update("Q", 0.0, "T", T_lo);
    return {get("p") * (1.0 + margin), crit.at("p") * (1.0 - margin)};
}

}  // namespace diagram
}  // namespace coolsolve
