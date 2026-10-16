#include "coolsolve/evaluator.h"
#include "coolsolve/constants.h"
#include "coolsolve/ees_string.h"
#include "coolsolve/lookup_table.h"
#include "coolsolve/units.h"
#include "coolsolve/fluids.h"
#include "coolsolve/user_hints.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "CoolProp.h"
#include "CoolPropLib.h"  // For set_config_bool C-style API
#include "AbstractState.h"
#include "DataStructures.h"
#include "HumidAirProp.h"

#include <memory>
#include <mutex>

namespace coolsolve {

// ============================================================================
// AbstractState Cache (thread-local for thread safety)
// ============================================================================

/// A new AbstractState for a CoolProp fluid name; a mixture `A[x]&B[y]` (see MixtureFluid) is built
/// from its components, with its mole fractions set.
static std::shared_ptr<CoolProp::AbstractState> createState(const std::string& backend,
                                                            const std::string& fluidName) {
    if (fluidName.find('&') == std::string::npos)
        return std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(backend, fluidName));
    std::string reason;
    auto mixture = std::dynamic_pointer_cast<MixtureFluid>(FluidRegistry::getFluid(fluidName, &reason));
    if (!mixture) throw std::runtime_error(reason.empty() ? "Unknown mixture: " + fluidName : reason);
    std::shared_ptr<CoolProp::AbstractState> state(
        CoolProp::AbstractState::factory(backend, mixture->getComponentList()));
    state->set_mole_fractions(mixture->getMoleFractions());
    return state;
}

/// CoolProp's own lower limits of a fluid: its minimum temperature and its minimum pressure (the triple
/// point). 0 when the backend has none (a value CoolProp cannot give: incompressible fluids have no triple point).
struct FluidLimits {
    double Tmin = 0.0;
    double pMin = 0.0;
};

/// Per-thread cache of AbstractState objects keyed by "backend::fluidName".
/// Avoids repeated string parsing, fluid lookup, and EOS initialisation that
/// PropsSI performs on every call.
struct AbstractStateCache {
    std::map<std::string, std::shared_ptr<CoolProp::AbstractState>> states;
    std::map<std::string, FluidLimits> limits;   ///< see fluidLimits
    size_t hits = 0;
    size_t misses = 0;

    std::shared_ptr<CoolProp::AbstractState> getOrCreate(
        const std::string& backend, const std::string& fluidName) {
        std::string key = backend + "::" + fluidName;
        auto it = states.find(key);
        if (it != states.end()) {
            ++hits;
            return it->second;
        }
        ++misses;
        auto state = createState(backend, fluidName);
        states[key] = state;
        return state;
    }

    /// The limits of a fluid, read from CoolProp once and cached; 0 for a limit it does not give.
    FluidLimits fluidLimits(const std::string& backend, const std::string& fluidName) {
        const std::string key = backend + "::" + fluidName;
        auto it = limits.find(key);
        if (it != limits.end()) return it->second;
        FluidLimits found;
        try {
            auto state = getOrCreate(backend, fluidName);
            try { found.Tmin = state->Tmin(); } catch (...) {}
            try { found.pMin = state->p_triple(); } catch (...) {}
            if (!(found.Tmin > 0.0)) found.Tmin = 0.0;   // also NaN
            if (!(found.pMin > 0.0)) found.pMin = 0.0;
        } catch (...) {}   // an unknown fluid: the property call itself reports it
        return limits[key] = found;
    }

    void clear() { states.clear(); limits.clear(); hits = misses = 0; }
};

static thread_local AbstractStateCache g_abstractStateCache;

/// Map two CoolProp parameter enums to the correct input_pairs enum.
static CoolProp::input_pairs getInputPair(CoolProp::parameters p1,
                                           CoolProp::parameters p2) {
    // Use CoolProp's generate_update_pair to determine the canonical ordering
    double out1, out2;
    return CoolProp::generate_update_pair(p1, 0.0, p2, 0.0, out1, out2);
}

/// MOLARMASS is returned in kg/kmol (EES convention); CoolProp's "M" is in kg/mol.
static constexpr double KG_PER_KMOL = 1000.0;

/// Get a property value from an already-updated AbstractState.
static double getOutputValue(CoolProp::AbstractState& state,
                             CoolProp::parameters param) {
    return state.keyed_output(param);
}

static ThermoInputKind thermoKindFromCoolProp(CoolProp::parameters param,
                                              const std::string& inputName) {
    std::string n = inputName;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (param == CoolProp::iP) return ThermoInputKind::Pressure;
    if (param == CoolProp::iT) return ThermoInputKind::Temperature;
    if (param == CoolProp::iHmass) return ThermoInputKind::Enthalpy;
    if (param == CoolProp::iSmass) return ThermoInputKind::Entropy;
    if (param == CoolProp::iQ) return ThermoInputKind::Quality;
    if (param == CoolProp::iDmass) {
        if (n == "v" || n == "volume") return ThermoInputKind::SpecificVolume;
        return ThermoInputKind::Density;
    }
    return ThermoInputKind::Pressure;
}

static ThermoInputKind thermoKindFromHAName(const std::string& name) {
    std::string n = name;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "t") return ThermoInputKind::Temperature;
    if (n == "p") return ThermoInputKind::Pressure;
    if (n == "w") return ThermoInputKind::HumidityRatio;
    if (n == "r") return ThermoInputKind::RelativeHumidity;
    if (n == "d") return ThermoInputKind::Density;
    if (n == "h") return ThermoInputKind::Enthalpy;
    return ThermoInputKind::Temperature;
}

void resetCoolPropUnitHints() {
    resetThermoInputHints();
}

// ============================================================================
// Simple Profiling Instrumentation
// ============================================================================
struct ProfilingStats {
    size_t propsSI_calls = 0;
    size_t hapropsSI_calls = 0;
    size_t abstractState_calls = 0;
    size_t analyticalDeriv_calls = 0;
    size_t block_evaluations = 0;
    size_t expression_evaluations = 0;
    double propsSI_time_ms = 0.0;
    double hapropsSI_time_ms = 0.0;
    double abstractState_time_ms = 0.0;
    double block_eval_time_ms = 0.0;
    double expr_eval_time_ms = 0.0;
    double coolprop_warmup_time_ms = 0.0;
    bool coolprop_warmup_performed = false;
    
    void reset() {
        propsSI_calls = 0;
        hapropsSI_calls = 0;
        abstractState_calls = 0;
        analyticalDeriv_calls = 0;
        block_evaluations = 0;
        expression_evaluations = 0;
        propsSI_time_ms = 0.0;
        hapropsSI_time_ms = 0.0;
        abstractState_time_ms = 0.0;
        block_eval_time_ms = 0.0;
        expr_eval_time_ms = 0.0;
        coolprop_warmup_time_ms = 0.0;
        coolprop_warmup_performed = false;
        g_abstractStateCache.clear();
    }
    
    std::string toString() const {
        std::ostringstream oss;
        oss << "\n=== CoolSolve Profiling Statistics ===\n";
        if (abstractState_calls > 0) {
            oss << "AbstractState calls: " << abstractState_calls
                << " (total: " << abstractState_time_ms << " ms, avg: "
                << (abstractState_time_ms / abstractState_calls) << " ms/call)\n";
            oss << "  Analytical derivative calls: " << analyticalDeriv_calls << "\n";
            oss << "  Cache hits/misses: " << g_abstractStateCache.hits
                << "/" << g_abstractStateCache.misses << "\n";
        }
        oss << "PropsSI calls: " << propsSI_calls 
            << " (total: " << propsSI_time_ms << " ms, avg: " 
            << (propsSI_calls > 0 ? propsSI_time_ms / propsSI_calls : 0.0) << " ms/call)\n";
        oss << "HAPropsSI calls: " << hapropsSI_calls 
            << " (total: " << hapropsSI_time_ms << " ms, avg: "
            << (hapropsSI_calls > 0 ? hapropsSI_time_ms / hapropsSI_calls : 0.0) << " ms/call)\n";
        oss << "Block evaluations: " << block_evaluations
            << " (total: " << block_eval_time_ms << " ms, avg: "
            << (block_evaluations > 0 ? block_eval_time_ms / block_evaluations : 0.0) << " ms/eval)\n";
        oss << "Expression evaluations: " << expression_evaluations << "\n";
        oss << "Total CoolProp time: " << (propsSI_time_ms + hapropsSI_time_ms + abstractState_time_ms) << " ms\n";
        oss << "Non-CoolProp block eval time: " << (block_eval_time_ms - propsSI_time_ms - hapropsSI_time_ms - abstractState_time_ms) << " ms\n";
        if (coolprop_warmup_performed) {
            oss << "CoolProp warmup (first PropsSI) time: " << coolprop_warmup_time_ms << " ms\n";
        }
        oss << "======================================\n\n";
        return oss.str();
    }
};

static ProfilingStats g_profilingStats;

// One-time CoolProp warmup to pay first-call cost up front.
double warmupCoolProp() {
    if (g_profilingStats.coolprop_warmup_performed) {
        return 0.0;
    }

    using clock = std::chrono::high_resolution_clock;
    auto start = clock::now();
    try {
        // A simple, representative PropsSI call that forces CoolProp
        // to load its fluid tables and internal data structures.
        // We intentionally do NOT use timedPropsSI here so that the
        // warmup cost is reported separately from normal PropsSI stats.
        (void)CoolProp::PropsSI("H", "T", 300.0, "P", 101325.0, "Water");
    } catch (...) {
        // If warmup fails, we still mark it as performed to avoid
        // repeatedly paying the cost on every run attempt.
    }
    auto end = clock::now();

    g_profilingStats.coolprop_warmup_performed = true;
    g_profilingStats.coolprop_warmup_time_ms =
        std::chrono::duration<double, std::milli>(end - start).count();

    return g_profilingStats.coolprop_warmup_time_ms;
}

// Wrapper for timed PropsSI calls
static double timedPropsSI(const std::string& output, const std::string& name1, double val1,
                           const std::string& name2, double val2, const std::string& fluid) {
    auto start = std::chrono::high_resolution_clock::now();
    double result = CoolProp::PropsSI(output, name1, val1, name2, val2, fluid);
    auto end = std::chrono::high_resolution_clock::now();
    g_profilingStats.propsSI_calls++;
    g_profilingStats.propsSI_time_ms += std::chrono::duration<double, std::milli>(end - start).count();
    return result;
}

// Wrapper for timed HAPropsSI calls
static double timedHAPropsSI(const std::string& output, const std::string& name1, double val1,
                             const std::string& name2, double val2, const std::string& name3, double val3) {
    auto start = std::chrono::high_resolution_clock::now();
    double result = HumidAir::HAPropsSI(output, name1, val1, name2, val2, name3, val3);
    auto end = std::chrono::high_resolution_clock::now();
    g_profilingStats.hapropsSI_calls++;
    g_profilingStats.hapropsSI_time_ms += std::chrono::duration<double, std::milli>(end - start).count();
    return result;
}

void resetProfilingStats() {
    g_profilingStats.reset();
}

std::string getProfilingStatsString() {
    return g_profilingStats.toString();
}

void applyCoolPropConfig(const CoolPropConfig& config) {
    // Set superancillary configuration using C-style API
    // Note: This must be called before any CoolProp fluid initialization
    ::set_config_bool("ENABLE_SUPERANCILLARIES", config.enableSuperancillaries);
}

// ============================================================================
// CoolProp Helper Functions
// ============================================================================
// Map common EES function names to CoolProp output parameters and unit types
struct CoolPropParamInfo {
    CoolProp::parameters param;
    UnitType unitType;
};

static CoolPropParamInfo getCoolPropOutputParam(const std::string& funcName) {
    std::string name = funcName;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    if (name == "enthalpy" || name == "h") return {CoolProp::iHmass, UnitType::SpecificEnergy};
    if (name == "entropy" || name == "s") return {CoolProp::iSmass, UnitType::SpecificEntropy};
    if (name == "density" || name == "rho" || name == "d") return {CoolProp::iDmass, UnitType::Density};
    if (name == "cp" || name == "specificheat" || name == "c_p") return {CoolProp::iCpmass, UnitType::SpecificHeat};
    if (name == "cv" || name == "c_v") return {CoolProp::iCvmass, UnitType::SpecificHeat};
    if (name == "temperature" || name == "t") return {CoolProp::iT, UnitType::Temperature};
    if (name == "pressure" || name == "p") return {CoolProp::iP, UnitType::Pressure};
    if (name == "quality" || name == "x" || name == "q") return {CoolProp::iQ, UnitType::Dimensionless};
    if (name == "t_sat" || name == "tsat") return {CoolProp::iT, UnitType::Temperature};
    if (name == "p_sat" || name == "psat") return {CoolProp::iP, UnitType::Pressure};
    if (name == "viscosity" || name == "mu" || name == "v") return {CoolProp::iviscosity, UnitType::Viscosity};
    if (name == "conductivity" || name == "k" || name == "l") return {CoolProp::iconductivity, UnitType::Conductivity};
    if (name == "speed_of_sound" || name == "soundspeed" || name == "a") return {CoolProp::ispeed_sound, UnitType::Dimensionless}; 
    if (name == "internalenergy" || name == "u" || name == "intenergy") return {CoolProp::iUmass, UnitType::SpecificEnergy};
    // Additional EES-style function names
    if (name == "volume" || name == "v") return {CoolProp::iDmass, UnitType::Density};  // Inverted to specific volume after PropsSI call
    if (name == "specheat" || name == "c") return {CoolProp::iCpmass, UnitType::SpecificHeat};
    if (name == "molarmass" || name == "mm") return {CoolProp::imolar_mass, UnitType::Dimensionless};
    // Additional EES thermophysical properties (priority 2)
    if (name == "prandtl") return {CoolProp::iPrandtl, UnitType::Dimensionless};
    if (name == "surfacetension") return {CoolProp::isurface_tension, UnitType::Dimensionless}; // N/m (no config unit)
    if (name == "compressibilityfactor") return {CoolProp::iZ, UnitType::Dimensionless};
    if (name == "isentropicexponent") return {CoolProp::iisentropic_expansion_coefficient, UnitType::Dimensionless};
    // Derived / composed properties (handled in evaluateCoolPropFunction before
    // falling into the generic 2-input path).  Return a sentinel so that
    // isThermo-detection still routes them through the CoolProp path.
    if (name == "kinematicviscosity") return {CoolProp::INVALID_PARAMETER, UnitType::Dimensionless};
    if (name == "thermaldiffusivity") return {CoolProp::INVALID_PARAMETER, UnitType::Dimensionless};
    // Pure-fluid constants (dispatched before the 2-input path).
    if (name == "t_crit" || name == "tcrit") return {CoolProp::iT_critical, UnitType::Temperature};
    if (name == "p_crit" || name == "pcrit") return {CoolProp::iP_critical, UnitType::Pressure};
    if (name == "v_crit" || name == "vcrit") return {CoolProp::irhomass_critical, UnitType::Dimensionless};
    if (name == "t_triple" || name == "ttriple") return {CoolProp::iT_triple, UnitType::Temperature};
    if (name == "p_triple" || name == "ptriple") return {CoolProp::INVALID_PARAMETER, UnitType::Pressure};
    if (name == "acentricfactor") return {CoolProp::INVALID_PARAMETER, UnitType::Dimensionless};
    // PHASE$ is string-valued, but we still register it to catch accidental
    // numeric-context calls gracefully.
    if (name == "phase$") return {CoolProp::INVALID_PARAMETER, UnitType::Dimensionless};
    
    // Humid Air specific outputs (mapped to dummy CoolProp params for type info)
    if (name == "humrat" || name == "w") return {CoolProp::INVALID_PARAMETER, UnitType::Dimensionless};
    if (name == "relhum" || name == "r") return {CoolProp::INVALID_PARAMETER, UnitType::Dimensionless};
    if (name == "wetbulb" || name == "b") return {CoolProp::iT, UnitType::Temperature};
    if (name == "dewpoint" || name == "d") return {CoolProp::iT, UnitType::Temperature};
    
    throw std::runtime_error("Unknown CoolProp output: " + funcName);
}

// Map EES input names to CoolProp parameters
static CoolProp::parameters getCoolPropInputParam(const std::string& inputName) {
    std::string name = inputName;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    if (name == "t" || name == "temperature") return CoolProp::iT;
    if (name == "p" || name == "pressure") return CoolProp::iP;
    if (name == "h" || name == "enthalpy") return CoolProp::iHmass;
    if (name == "s" || name == "entropy") return CoolProp::iSmass;
    if (name == "d" || name == "rho" || name == "density") return CoolProp::iDmass;
    if (name == "v" || name == "volume") return CoolProp::iDmass;
    if (name == "q" || name == "x" || name == "quality") return CoolProp::iQ;  // EES uses x for quality
    if (name == "u" || name == "internalenergy" || name == "intenergy") return CoolProp::iUmass;
    
    // For Humid Air, we'll handle these separately
    if (name == "w" || name == "humrat") return CoolProp::INVALID_PARAMETER;
    if (name == "r" || name == "relhum") return CoolProp::INVALID_PARAMETER;
    if (name == "b" || name == "wetbulb") return CoolProp::iT;
    if (name == "d" || name == "dewpoint") return CoolProp::iT;
    
    throw std::runtime_error("Unknown CoolProp input: " + inputName);
}

// Humid Air parameter mapping. CoolProp's HAPropsSI gives H, S, U, V, C and CV per kg of DRY air (as EES does for
// AirH2O); Rho is not a CoolProp output but 1/V (see getHAOutputParam).
enum class HAParam { T, P, W, R, D, B, H, S, V, M, K, C, CV, U, Rho, Unknown };

static HAParam getHAParam(const std::string& name) {
    std::string n = name;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "t" || n == "temperature") return HAParam::T;
    if (n == "p" || n == "pressure") return HAParam::P;
    if (n == "w" || n == "humrat") return HAParam::W;
    if (n == "r" || n == "relhum") return HAParam::R;
    if (n == "d" || n == "dewpoint") return HAParam::D;
    if (n == "b" || n == "wetbulb") return HAParam::B;
    if (n == "h" || n == "enthalpy") return HAParam::H;
    if (n == "s" || n == "entropy") return HAParam::S;
    if (n == "v" || n == "volume") return HAParam::V;
    if (n == "mu" || n == "viscosity") return HAParam::M;
    if (n == "k" || n == "conductivity") return HAParam::K;
    if (n == "cp" || n == "specheat") return HAParam::C;
    if (n == "cv") return HAParam::CV;
    if (n == "u" || n == "internalenergy" || n == "intenergy") return HAParam::U;
    return HAParam::Unknown;
}

/// The humid-air output asked by the EES function `name`: getHAParam, and `density` (EES: 1/volume per kg of dry air,
/// like the other properties of AirH2O; the density of the mixture is density*(1+w)). `density` is an output only:
/// as the name of an input it stays unknown (it must not be read as a volume).
static HAParam getHAOutputParam(const std::string& name) {
    std::string n = name;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "density" || n == "rho") return HAParam::Rho;
    return getHAParam(n);
}

static std::string haParamToString(HAParam p) {
    switch (p) {
        case HAParam::T: return "T";
        case HAParam::P: return "P";
        case HAParam::W: return "W";
        case HAParam::R: return "R";
        case HAParam::D: return "D";
        case HAParam::B: return "B";
        case HAParam::H: return "H";
        case HAParam::S: return "S";
        case HAParam::V: return "V";
        case HAParam::M: return "M";
        case HAParam::K: return "K";
        case HAParam::C: return "C";
        case HAParam::CV: return "CV";
        case HAParam::U: return "U";
        case HAParam::Rho: return "V";   // 1/Vda
        default: return "";
    }
}

static UnitType haParamToUnitType(HAParam p) {
    switch (p) {
        case HAParam::T:
        case HAParam::D:
        case HAParam::B: return UnitType::Temperature;
        case HAParam::P: return UnitType::Pressure;
        case HAParam::H:
        case HAParam::U: return UnitType::SpecificEnergy;
        case HAParam::S: return UnitType::SpecificEntropy;
        case HAParam::V: return UnitType::Density; // m3/kg is inverse density
        case HAParam::Rho: return UnitType::Density;
        case HAParam::M: return UnitType::Viscosity;
        case HAParam::K: return UnitType::Conductivity;
        case HAParam::C:
        case HAParam::CV: return UnitType::SpecificHeat;
        default: return UnitType::Dimensionless;
    }
}

/// Why CoolProp's humid-air model returned no value (NaN/Inf: it does not throw) for the inputs `valsSI` (of the
/// kinds `params`), as the end of an error message in the user's units: "with inputs T=1026.85 C, P=800000 Pa, W=0:
/// <reason>". The reason names the inputs outside the validity range of the model (`check_bounds` of CoolProp's
/// HumidAirProp.cpp: T, Twb and Tdp from -143.15 to 350 C, P from 10 Pa to 10 MPa, W from 0 to 10, R from 0 to 1)
/// and, above 350 C, the way out for dry air; without such an input it is CoolProp's own message.
static std::string humidAirFailureDetail(const std::vector<HAParam>& params, const std::vector<double>& valsSI,
                                         const UnitSystem& units) {
    auto userUnit = [&](UnitType t) -> std::string {
        return t == UnitType::Temperature ? units.temperature : t == UnitType::Pressure ? units.pressure :
               t == UnitType::SpecificEnergy ? units.specific_energy : t == UnitType::SpecificEntropy ? units.specific_entropy :
               t == UnitType::Density ? units.density : "";
    };
    auto inUserUnits = [&](HAParam p, double valueSI) {   // "<number> <unit>", the unit when there is one
        const UnitType t = haParamToUnitType(p);
        std::ostringstream oss;
        oss << std::setprecision(10);
        if (t == UnitType::Dimensionless) oss << valueSI;
        else oss << UnitConverter::fromSI(valueSI, t, userUnit(t)) << " " << userUnit(t);
        return oss.str();
    };

    std::string inputs, reason;
    bool tooHot = false;
    for (size_t i = 0; i < params.size(); ++i) {
        const HAParam p = params[i];
        const double v = valsSI[i];
        const std::string state = haParamToString(p) + "=" + inUserUnits(p, v);
        inputs += (i ? ", " : "") + state;

        double lo = 0.0, hi = 0.0;
        switch (p) {
            case HAParam::T: case HAParam::B: case HAParam::D: lo = 130.0; hi = 623.15; break;
            case HAParam::P: lo = 10.0; hi = 1e7; break;
            case HAParam::W: lo = 0.0; hi = 10.0; break;
            case HAParam::R: lo = 0.0; hi = 1.0; break;
            default: continue;   // the other inputs have no limit
        }
        if (std::isfinite(v) && v >= lo && v <= hi) continue;
        reason += (reason.empty() ? "" : " ") + state;
        reason += !std::isfinite(v) ? " is not a finite number."
            : " is outside the validity range of CoolProp's humid-air model (" + inUserUnits(p, lo) + " to " +
              inUserUnits(p, hi) + ").";
        if (std::isfinite(v) && v > hi && haParamToUnitType(p) == UnitType::Temperature) tooHot = true;
    }
    if (tooHot) reason += " For dry air above " + inUserUnits(HAParam::T, 623.15) + " use the ideal-gas fluid Air.";
    if (reason.empty()) {
        const std::string coolProp = CoolProp::get_global_param_string("errstring");
        reason = coolProp.empty() ? "no further information from CoolProp." : "CoolProp: " + coolProp;
    }
    return "with inputs " + inputs + ": " + reason;
}

// Map CoolProp parameter to PropsSI string
static std::string paramToString(CoolProp::parameters param) {
    switch (param) {
        case CoolProp::iT: return "T";
        case CoolProp::iP: return "P";
        case CoolProp::iHmass: return "H";
        case CoolProp::iSmass: return "S";
        case CoolProp::iDmass: return "D";
        case CoolProp::iQ: return "Q";
        case CoolProp::iCpmass: return "C";
        case CoolProp::iCvmass: return "CVMASS";
        case CoolProp::iviscosity: return "V";
        case CoolProp::iconductivity: return "L";
        case CoolProp::ispeed_sound: return "A";
        case CoolProp::iUmass: return "U";
        case CoolProp::imolar_mass: return "M";
        case CoolProp::iPrandtl: return "Prandtl";
        case CoolProp::isurface_tension: return "I";
        case CoolProp::iZ: return "Z";
        case CoolProp::iisentropic_expansion_coefficient: return "isentropic_expansion_coefficient";
        default: throw std::runtime_error("Cannot convert parameter to string");
    }
}

/// The fluid of a property call. An unknown name, or an invalid mixture string, is an error that says
/// why; a mixture also gets its once-per-run warning (its properties are predictive).
static std::shared_ptr<Fluid> resolveFluid(const std::string& name, DiagnosticCollector* diag) {
    std::string reason;
    auto fluid = FluidRegistry::getFluid(name, &reason);
    if (!fluid) throw std::runtime_error(reason.empty() ? unknownFluidErrorMessage(name) : reason);
    if (fluid->getType() == FluidType::Mixture) warnMixtureFluid(diag, fluid->getName());
    return fluid;
}

// ============================================================================
// ExpressionEvaluator Implementation
// ============================================================================

// Internal procedural statement evaluator (forward declaration)
class ProceduralEvaluator {
public:
    /// Runs the statements of the body of a FUNCTION or PROCEDURE, in order; `routine` ("PROCEDURE 'p'") names
    /// it in the errors. A GOTO jumps to its label, in the block that holds the GOTO or in one that contains it.
    static void run(ExpressionEvaluator& eval, const std::vector<StmtPtr>& body, const std::string& routine);

private:
    /// GOTOs executed by one run of a body beyond which it is taken for a loop that never ends.
    static constexpr long kMaxGotos = 100000;
    /// Iterations of one run of a REPEAT loop beyond which it is taken for a loop that never ends.
    static constexpr long kMaxRepeatIterations = 100000;

    /// The GOTO being carried out: while `pending`, the blocks are left, innermost first, up to the one that
    /// holds `label`.
    struct Jump {
        const std::string& routine;
        bool pending = false;
        int label = 0;
        long executed = 0;
    };

    static void runBlock(ExpressionEvaluator& eval, const std::vector<StmtPtr>& block, Jump& jump);
    static void execute(ExpressionEvaluator& eval, const StmtPtr& stmt, Jump& jump);
};

ExpressionEvaluator::ExpressionEvaluator(size_t numVariables, const CoolPropConfig& config)
    : numVariables_(numVariables), coolpropConfig_(config) {}

void ExpressionEvaluator::setVariable(const std::string& name, const ADValue& value) {
    variables_[name] = value;
}

void ExpressionEvaluator::setStringVariable(const std::string& name, const std::string& value) {
    stringVariables_[name] = value;
}

std::string ExpressionEvaluator::getStringVariable(const std::string& name) const {
    auto it = stringVariables_.find(name);
    if (it != stringVariables_.end()) {
        return it->second;
    }
    throw std::runtime_error("String variable not found: " + name);
}

bool ExpressionEvaluator::hasStringVariable(const std::string& name) const {
    return stringVariables_.find(name) != stringVariables_.end();
}

ADValue ExpressionEvaluator::getVariable(const std::string& name) const {
    auto it = variables_.find(name);
    if (it != variables_.end()) {
        return it->second;
    }
    
    // A built-in constant (pi, g#, R#, ...) is available in every scope: the main program gets it as
    // an IR variable (IR::addVariable), the body of a FUNCTION or PROCEDURE has no IR variables and
    // would otherwise read it as an unset local (CS-BUG-PI-FUNCTION). The lookup above comes first: a
    // body that sets a variable named like a constant (a local `c#`) keeps its own value, as before.
    if (Constants::isConstant(name)) {
        return ADValue::constant(Constants::getValue(name), numVariables_);
    }

    // If not found, use default value of 1.0 and record as missing
    missingVariables_.insert(name);
    return ADValue::constant(1.0, numVariables_);
}

bool ExpressionEvaluator::hasVariable(const std::string& name) const {
    return true; // All variables are now "available" with a default value
}

void ExpressionEvaluator::clear() {
    variables_.clear();
    stringVariables_.clear();
    missingVariables_.clear();
}

void ExpressionEvaluator::registerFunction(const FunctionDefinition& func) {
    userFunctions_[func.name] = func;
}

void ExpressionEvaluator::registerProcedure(const ProcedureDefinition& proc) {
    userProcedures_[proc.name] = proc;
}

std::string ExpressionEvaluator::resolveVariableName(const Variable& var) {
    if (var.indices.empty()) return var.name;
    
    std::string result = var.name + "[";
    for (size_t i = 0; i < var.indices.size(); ++i) {
        if (i > 0) result += ",";
        ADValue idxVal = evaluate(var.indices[i]);
        double val = idxVal.value;
        if (val == static_cast<int>(val)) {
            result += std::to_string(static_cast<int>(val));
        } else {
            result += std::to_string(val);
        }
    }
    result += "]";
    return result;
}

void ExpressionEvaluator::evaluateProcedureCall(const ProcedureCall& call) {
    // Built-in PSYCHPROPS procedure (EES compatibility).
    // Signature:
    //   CALL psychprops(T, P, RH : T, v, h, s, u, W, R, Twb, Tdp)
    // where RH is the relative humidity (0..1).  Inputs are in the configured
    // unit system; outputs are assembled in the configured unit system too.
    {
        std::string lowerCallName = call.name;
        std::transform(lowerCallName.begin(), lowerCallName.end(),
                       lowerCallName.begin(), ::tolower);
        if (lowerCallName == "psychprops") {
            if (call.inputArgs.size() != 3) {
                throw std::runtime_error(
                    "psychprops expects 3 inputs (T, P, R), got " +
                    std::to_string(call.inputArgs.size()));
            }
            if (call.outputVars.size() < 1 || call.outputVars.size() > 9) {
                throw std::runtime_error(
                    "psychprops expects 1 to 9 outputs (T, v, h, s, u, W, R, Twb, Tdp), got " +
                    std::to_string(call.outputVars.size()));
            }
            
            ADValue tIn = evaluate(call.inputArgs[0]);
            ADValue pIn = evaluate(call.inputArgs[1]);
            ADValue rIn = evaluate(call.inputArgs[2]);
            
            const UnitSystem& units = coolpropConfig_.units;
            double tSI = UnitConverter::toSI(tIn.value, UnitType::Temperature, units.temperature);
            double pSI = UnitConverter::toSI(pIn.value, UnitType::Pressure, units.pressure);
            double rVal = rIn.value;
            
            auto callHA = [&](const std::string& out) -> double {
                return timedHAPropsSI(out, "T", tSI, "P", pSI, "R", rVal);
            };
            // Assemble all 9 outputs in SI then convert to configured units.
            struct OutputSpec { std::string name; double valueSI; UnitType ut; std::string unit; };
            std::vector<OutputSpec> outs;
            outs.push_back({"T", tSI, UnitType::Temperature, units.temperature});
            outs.push_back({"v", callHA("Vha"), UnitType::Dimensionless, ""});       // m^3/kg dry air
            outs.push_back({"h", callHA("Hha"), UnitType::SpecificEnergy, units.specific_energy}); // J/kg dry air
            outs.push_back({"s", callHA("Sha"), UnitType::SpecificEntropy, units.specific_entropy});
            // Internal energy u = h - p*v (dry-basis); CoolProp has no direct "Uha"
            double v_si = outs[1].valueSI;
            double h_si = outs[2].valueSI;
            outs.push_back({"u", h_si - pSI * v_si, UnitType::SpecificEnergy, units.specific_energy});
            outs.push_back({"W", callHA("W"), UnitType::Dimensionless, ""});
            outs.push_back({"R", rVal, UnitType::Dimensionless, ""});
            outs.push_back({"Twb", callHA("Twb"), UnitType::Temperature, units.temperature});
            outs.push_back({"Tdp", callHA("Tdp"), UnitType::Temperature, units.temperature});
            
            for (size_t i = 0; i < call.outputVars.size(); ++i) {
                std::string outVarName = resolveVariableName(call.outputVars[i]);
                double raw = outs[i].valueSI;
                if (!std::isfinite(raw)) {
                    throw std::runtime_error("psychprops: CoolProp returned non-finite for '" +
                                             outs[i].name + "' " +
                                             humidAirFailureDetail({HAParam::T, HAParam::P, HAParam::R},
                                                                   {tSI, pSI, rVal}, units));
                }
                double finalVal = (outs[i].ut == UnitType::Dimensionless)
                    ? raw
                    : UnitConverter::fromSI(raw, outs[i].ut, outs[i].unit);
                setVariable(outVarName, ADValue::constant(finalVal, numVariables_));
            }
            return;
        }
    }
    
    auto it = userProcedures_.find(call.name);
    if (it == userProcedures_.end()) {
        std::string lowerName = call.name;
        std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
        // Built-in CALL ERROR('message XXXF1', value): stop the calculation (EES).
        if (lowerName == "error") throw std::runtime_error("CALL ERROR: " + userMessage(call, true));
        // Built-in CALL WARNING('message XXXF1', value): a warning, the calculation goes on (EES). Only the
        // evaluation at the solution has a sink; the iterations skip it (see setCallWarningSink).
        if (lowerName == "warning") {
            if (!callWarningSink_) return;
            Diagnostic warning;
            warning.severity = DiagnosticSeverity::Warning;
            warning.code = "W001";
            warning.message = userMessage(call, false);
            warning.line = call.inputArgs.empty() ? 0 : call.inputArgs[0]->sourceLineNumber;   // of the CALL statement
            warning.source = "evaluator";
            callWarningSink_->push_back(std::move(warning));
            return;
        }
        // Built-in CALL Colebrook(Re, RR : f) (EES library; the parser also takes CALL Colebrook(Re, RR, f)):
        // the Darcy friction factor of the Colebrook-White equation, with its derivatives (AD).
        if (lowerName == "colebrook") {
            if (call.inputArgs.size() != 2 || call.outputVars.size() != 1)   // the parser refuses these calls
                throw std::runtime_error("CALL Colebrook takes 2 inputs (Re, RR) and 1 output (f)");
            setVariable(resolveVariableName(call.outputVars[0]),
                        colebrookFrictionFactor(evaluate(call.inputArgs[0]), evaluate(call.inputArgs[1])));
            return;
        }
        throw std::runtime_error("Unknown procedure: " + call.name);
    }
    
    const auto& proc = it->second;
    if (call.inputArgs.size() != proc.inputs.size()) {
        throw std::runtime_error("Procedure " + call.name + " expected " + std::to_string(proc.inputs.size()) + " inputs, got " + std::to_string(call.inputArgs.size()));
    }
    
    ExpressionEvaluator localEval = makeBodyEvaluator("PROCEDURE", proc.name);
    
    for (size_t i = 0; i < proc.inputs.size(); ++i) {
        const std::string& inputName = proc.inputs[i];
        if (!inputName.empty() && inputName.back() == '$') {
            localEval.setStringVariable(inputName, evaluateString(call.inputArgs[i]));
        } else {
            localEval.setVariable(inputName, evaluate(call.inputArgs[i]));
        }
    }
    
    // Evaluate body
    ProceduralEvaluator::run(localEval, proc.body, localEval.routine_);
    
    // Update output variables in the current scope
    if (call.outputVars.size() != proc.outputs.size()) {
        throw std::runtime_error("Procedure " + call.name + " returns " + std::to_string(proc.outputs.size()) + " outputs, but " + std::to_string(call.outputVars.size()) + " were provided");
    }
    
    for (size_t i = 0; i < proc.outputs.size(); ++i) {
        std::string outParamName = proc.outputs[i];
        std::string callOutVarName = resolveVariableName(call.outputVars[i]);
        if (!outParamName.empty() && outParamName.back() == '$') {
            setStringVariable(callOutVarName, localEval.getStringVariable(outParamName));
        } else {
            setVariable(callOutVarName, localEval.getVariable(outParamName));
        }
    }
}

std::string ExpressionEvaluator::userMessage(const ProcedureCall& call, bool isError) {
    const std::string procedure = isError ? "CALL ERROR" : "CALL WARNING";
    const auto& args = call.inputArgs;
    if (args.empty() || args.size() > 2)   // the parser refuses these calls; a built model cannot hold one
        throw std::runtime_error(procedure + " takes a message and a value at most, and at least one argument");
    const bool hasMessage = isStringExpression(args[0]);
    if (!hasMessage && args.size() > 1)
        throw std::runtime_error("the first argument of " + procedure + " must be the message (a string) when a value follows");
    const ExprPtr valueArg = hasMessage ? (args.size() > 1 ? args[1] : nullptr) : args[0];

    const std::string where = routineName_.empty() ? "the main program" : routineName_;
    const std::string message = hasMessage ? evaluateString(args[0])
        : isError ? "Calculations have been halted because a parameter is out of range.  The value of the parameter is XXX in function " + where + "."
                  : "A warning message was issued due to the value XXX in " + where + ".";
    if (!valueArg) return message;

    // The value: a string for XXX$ (a string variable), else a number in the display format after XXX
    const bool valueIsText = isStringExpression(valueArg);
    const std::string text = valueIsText ? evaluateString(valueArg) : std::string();
    const double number = valueIsText ? 0.0 : evaluate(valueArg).value;
    bool found = false;
    std::string out = eesSubstituteValue(message, [&](const std::string& format) {
        if (valueIsText) return text;
        return eesFormatNumber(number, format == "$" ? "A" : format);
    }, &found);
    if (!found && hasMessage) out += " [" + (valueIsText ? text : eesFormatNumber(number)) + "]";   // keep the value
    return out;
}

ADValue ExpressionEvaluator::evaluate(const ExprPtr& expr) {
    if (!expr) {
        throw std::runtime_error("Null expression pointer");
    }
    
    if (expr->is<NumberLiteral>()) {
        return evaluateNumber(expr->as<NumberLiteral>());
    } else if (expr->is<Variable>()) {
        return evaluateVariable(expr->as<Variable>());
    } else if (expr->is<UnaryOp>()) {
        return evaluateUnaryOp(expr->as<UnaryOp>());
    } else if (expr->is<BinaryOp>()) {
        return evaluateBinaryOp(expr->as<BinaryOp>());
    } else if (expr->is<FunctionCall>()) {
        return evaluateFunctionCall(expr->as<FunctionCall>());
    } else if (expr->is<StringLiteral>()) {
        return ADValue::constant(0.0, numVariables_);
    }
    
    throw std::runtime_error("Unknown expression type");
}

std::string ExpressionEvaluator::evaluateString(const ExprPtr& expr) {
    if (!expr) {
        throw std::runtime_error("Null expression pointer");
    }
    
    if (expr->is<StringLiteral>()) {
        return expr->as<StringLiteral>().value;
    } else if (expr->is<Variable>()) {
        return getStringVariable(resolveVariableName(expr->as<Variable>()));
    } else if (expr->is<FunctionCall>()) {
        const auto& call = expr->as<FunctionCall>();
        std::string lname = call.name;
        std::transform(lname.begin(), lname.end(), lname.begin(), ::tolower);
        
        // STRING$(value) — format a numeric value as a compact decimal string.
        if (lname == "string$" && call.args.size() == 1 && call.namedArgs.empty()) {
            // EES prints integer values without decimal, otherwise up to ~10
            // significant digits, trimming trailing zeros (auto format, shared with
            // the messages of CALL WARNING / CALL ERROR).
            return eesFormatNumber(evaluate(call.args[0]).value);
        }
        
        // PHASE$(Fluid, T=..., P=...) — return the thermodynamic phase as a
        // human-readable string.
        if (lname == "phase$") {
            const std::string eesFluidName = fluidNameOf(call.args.empty() ? nullptr : call.args[0]);
            auto fluid = resolveFluid(eesFluidName, diagnostics_);
            if (fluid->getType() == FluidType::HumidAir) {
                return "gas";
            }
            if (fluid->getType() == FluidType::Unknown) {
                auto unsupported = std::dynamic_pointer_cast<UnsupportedFluid>(fluid);
                throw std::runtime_error("Fluid '" + eesFluidName + "' is not supported: " +
                    (unsupported ? unsupported->getReason() : "Unknown reason"));
            }
            
            // Evaluate named-arg state inputs in SI units
            std::map<std::string, double> stateSI;
            const UnitSystem& units = coolpropConfig_.units;
            for (const auto& [argName, argExpr] : call.namedArgs) {
                std::string lk = argName;
                std::transform(lk.begin(), lk.end(), lk.begin(), ::tolower);
                ADValue v = evaluate(argExpr);
                double raw = v.value;
                double siv = raw;
                if (lk == "t") siv = UnitConverter::toSI(raw, UnitType::Temperature, units.temperature);
                else if (lk == "p") siv = UnitConverter::toSI(raw, UnitType::Pressure, units.pressure);
                else if (lk == "h") siv = UnitConverter::toSI(raw, UnitType::SpecificEnergy, units.specific_energy);
                else if (lk == "s") siv = UnitConverter::toSI(raw, UnitType::SpecificEntropy, units.specific_entropy);
                else if (lk == "d" || lk == "rho") siv = UnitConverter::toSI(raw, UnitType::Density, units.density);
                stateSI[lk] = siv;
            }
            if (stateSI.size() != 2) {
                throw std::runtime_error("phase$() requires exactly 2 state inputs (e.g. T=..., P=...)");
            }
            
            try {
                auto state = g_abstractStateCache.getOrCreate(
                    coolpropConfig_.getBackendString(), fluid->getCoolPropName());
                auto pairStr = [](const std::string& k) {
                    if (k == "t") return CoolProp::iT;
                    if (k == "p") return CoolProp::iP;
                    if (k == "h") return CoolProp::iHmass;
                    if (k == "s") return CoolProp::iSmass;
                    if (k == "d" || k == "rho") return CoolProp::iDmass;
                    if (k == "q" || k == "x") return CoolProp::iQ;
                    return CoolProp::INVALID_PARAMETER;
                };
                auto it1 = stateSI.begin();
                auto it2 = std::next(it1);
                CoolProp::parameters p1 = pairStr(it1->first);
                CoolProp::parameters p2 = pairStr(it2->first);
                double iv1 = 0, iv2 = 0;
                CoolProp::input_pairs ipair = CoolProp::generate_update_pair(
                    p1, it1->second, p2, it2->second, iv1, iv2);
                state->update(ipair, iv1, iv2);
                int ph = static_cast<int>(state->phase());
                switch (ph) {
                    case CoolProp::iphase_liquid: return "liquid";
                    case CoolProp::iphase_supercritical: return "supercritical";
                    case CoolProp::iphase_supercritical_gas: return "supercritical_gas";
                    case CoolProp::iphase_supercritical_liquid: return "supercritical_liquid";
                    case CoolProp::iphase_critical_point: return "critical_point";
                    case CoolProp::iphase_gas: return "gas";
                    case CoolProp::iphase_twophase: return "twophase";
                    case CoolProp::iphase_unknown: return "unknown";
                    case CoolProp::iphase_not_imposed: return "not_imposed";
                    default: return "unknown";
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("phase$() failed for fluid '" + eesFluidName + "': " + e.what());
            }
        }
        
        // LOOKUP$('table', row, col) — the text of a lookup-table cell.
        if (lname == "lookup$") return evaluateLookupString(call);

        // UPPERCASE$(s$) / LOWERCASE$(s$) — the string in upper / lower case.
        if (lname == "uppercase$" || lname == "lowercase$") {
            if (call.args.size() != 1 || !call.namedArgs.empty())
                throw std::runtime_error(call.name + "(): expects 1 argument: a string");
            if (!isStringExpression(call.args[0]))
                throw std::runtime_error(call.name + "(): the argument must be a string ('text', a string variable "
                                         "ending in $, or a string function), not a number; STRING$(x) converts a number");
            const std::string text = evaluateString(call.args[0]);
            return lname == "uppercase$" ? eesUpperCase(text) : eesFoldCase(text);
        }

        throw std::runtime_error("String-valued function not supported: " + call.name);
    }
    
    throw std::runtime_error("Expression is not a string literal or variable");
}

ADValue ExpressionEvaluator::evaluateNumber(const NumberLiteral& num) {
    return ADValue::constant(num.value, numVariables_);
}

ADValue ExpressionEvaluator::evaluateVariable(const Variable& var) {
    std::string name = resolveVariableName(var);
    
    if (hasVariable(name)) {
        return getVariable(name);
    }
    
    if (hasVariable(var.name)) {
        return getVariable(var.name);
    }
    
    throw std::runtime_error("Undefined variable: " + name);
}

ADValue ExpressionEvaluator::evaluateUnaryOp(const UnaryOp& op) {
    ADValue operand = evaluate(op.operand);
    
    if (op.op == "-") {
        return -operand;
    } else if (op.op == "+") {
        return operand;
    }
    
    throw std::runtime_error("Unknown unary operator: " + op.op);
}

std::string ExpressionEvaluator::fluidNameOf(const ExprPtr& arg) {
    if (!arg) return "";
    if (arg->is<StringLiteral>()) return arg->as<StringLiteral>().value;
    if (arg->is<Variable>()) {
        const auto& var = arg->as<Variable>();
        std::string resolvedName = resolveVariableName(var);
        if (hasStringVariable(resolvedName)) return getStringVariable(resolvedName);
        // Fallback: use the name itself (minus $ if present), the unquoted fluid name of EES
        std::string name = var.name;
        if (!name.empty() && name.back() == '$') name.pop_back();
        return name;
    }
    // A string function, such as LOOKUP$('table', i, 1)
    if (arg->is<FunctionCall>() && isStringExpression(arg)) return evaluateString(arg);
    return "";
}

ADValue ExpressionEvaluator::evaluateBinaryOp(const BinaryOp& op) {
    // Conditions of IF statements compare strings with = and <>, without regard to case
    // as EES does (CS-BUG-STRING-CASE; the helper is shared with the $IF directives).
    if ((op.op == "=" || op.op == "<>") && (isStringExpression(op.left) || isStringExpression(op.right))) {
        const bool equal = eesStringEqual(evaluateString(op.left), evaluateString(op.right));
        return ADValue::constant(equal == (op.op == "=") ? 1.0 : -1.0, numVariables_);
    }

    ADValue left = evaluate(op.left);
    ADValue right = evaluate(op.right);
    
    if (op.op == "+") {
        return left + right;
    } else if (op.op == "-") {
        return left - right;
    } else if (op.op == "*") {
        return left * right;
    } else if (op.op == "/") {
        return left / right;
    } else if (op.op == "^") {
        return pow(left, right);
    } else if (op.op == ">") {
        return ADValue::constant(left.value > right.value ? 1.0 : -1.0, left.gradient.size());
    } else if (op.op == "<") {
        return ADValue::constant(left.value < right.value ? 1.0 : -1.0, left.gradient.size());
    } else if (op.op == ">=") {
        return ADValue::constant(left.value >= right.value ? 1.0 : -1.0, left.gradient.size());
    } else if (op.op == "<=") {
        return ADValue::constant(left.value <= right.value ? 1.0 : -1.0, left.gradient.size());
    } else if (op.op == "=") {
        return ADValue::constant(left.value == right.value ? 1.0 : -1.0, left.gradient.size());
    } else if (op.op == "<>") {
        return ADValue::constant(left.value != right.value ? 1.0 : -1.0, left.gradient.size());
    } else if (op.op == "and") {
        return ADValue::constant(left.value > 0.5 && right.value > 0.5 ? 1.0 : -1.0, left.gradient.size());
    } else if (op.op == "or") {
        return ADValue::constant(left.value > 0.5 || right.value > 0.5 ? 1.0 : -1.0, left.gradient.size());
    }
    
    throw std::runtime_error("Unknown binary operator: " + op.op);
}

/// The lookup table / interpolation functions: their arguments mix strings (table, columns) and numbers,
/// so they are evaluated from the call itself (evaluateLookupFunction).
static bool isLookupFunction(const std::string& lowerName) {
    static const std::unordered_set<std::string> names = {
        "lookup", "lookup$row", "lookupcol", "lookupcol1", "lookupcellempty",
        "tablevalue", "tablevalue#", "tablerun#",
        "interpolate", "interpolate1", "interpolate2", "interpolate2dm",
        "nlookuprows", "nlookupcolumns",
        "sumlookup", "avglookup", "maxlookup", "minlookup", "stddevlookup",
    };
    return names.count(lowerName) > 0;
}

    ADValue ExpressionEvaluator::evaluateFunctionCall(const FunctionCall& func) {
    std::string name = func.name;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    // Check for user-defined function first
    auto it = userFunctions_.find(name);
    if (it != userFunctions_.end()) {
        return evaluateUserFunction(it->second, func);
    }
    
    // Check for user-defined procedure called with inline syntax (single output)
    // EES allows: result = proc(a, b, c) when the procedure has exactly one output
    {
        auto procIt = userProcedures_.find(name);
        if (procIt != userProcedures_.end()) {
            const auto& proc = procIt->second;
            if (func.args.size() == proc.inputs.size()) {
                if (proc.outputs.size() != 1) {
                    throw std::runtime_error("Procedure " + func.name + " has " +
                        std::to_string(proc.outputs.size()) +
                        " outputs and cannot be called as a function (inline syntax requires exactly 1 output)");
                }
                // Build a temporary ProcedureCall and evaluate it
                ProcedureCall tempCall;
                tempCall.name = name;
                tempCall.inputArgs = func.args;
                Variable outVar;
                outVar.name = "__proc_inline_" + name;
                tempCall.outputVars = {outVar};
                evaluateProcedureCall(tempCall);
                return getVariable("__proc_inline_" + name);
            }
        }
    }
    
    // ISIDEALGAS(Fluid$) — 1 if the substance is modelled as an ideal gas (the IdealGas fluids of the registry:
    // Air, N2, CO2, H2O, ...), 0 for the real fluids, humid air, incompressibles and mixtures; an unknown fluid is
    // an error, as in a property call. A constant, so its derivative is 0. Before the property calls: it takes
    // the fluid only, and no state (T=..., P=...).
    if (name == "isidealgas") {
        if (func.args.size() != 1 || !func.namedArgs.empty())
            throw std::runtime_error(func.name + "(): expects 1 argument, the fluid (a name, a quoted string or a "
                                     "string variable), no state inputs such as T= or P=");
        auto fluid = resolveFluid(fluidNameOf(func.args[0]), diagnostics_);
        return ADValue::constant(fluid->getType() == FluidType::IdealGas ? 1.0 : 0.0, numVariables_);
    }

    // Check if it's a thermodynamic function handled by CoolProp
    bool isThermo = false;
    try {
        getCoolPropOutputParam(name);
        isThermo = true;
    } catch (...) {
        // Not a thermo function
    }
    
    // Named arguments make a property call (T=20), except in the EES form of INTERPOLATE (Col=value)
    if (isThermo || (!func.namedArgs.empty() && !isLookupFunction(name))) {
        return evaluateCoolPropFunction(func);
    }
    
    // LOOKUP$ gives text (see evaluateString), never a number
    if (name == "lookup$") {
        throw std::runtime_error(func.name + " returns a string and cannot be used as a number; use it where a "
            "string is expected (a fluid argument, a string comparison, a string variable set in a "
            "PROCEDURE or FUNCTION), or STRINGVAL(LOOKUP$(...)) to convert a numeric cell" + lookupWhere());
    }

    // UPPERCASE$ and LOWERCASE$ give text too
    if (name == "uppercase$" || name == "lowercase$") {
        throw std::runtime_error(func.name + " returns a string and cannot be used as a number; use it where a "
            "string is expected (a fluid argument, a string comparison, a string variable set in a "
            "PROCEDURE or FUNCTION)");
    }

    // STRINGVAL(s$) — parse a numeric string to a number.
    if (name == "stringval" && func.args.size() == 1 && func.namedArgs.empty()) {
        std::string s = evaluateString(func.args[0]);
        // Trim whitespace
        size_t a = s.find_first_not_of(" \t\r\n");
        size_t b = s.find_last_not_of(" \t\r\n");
        if (a == std::string::npos) {
            throw std::runtime_error("stringval(): empty string");
        }
        std::string trimmed = s.substr(a, b - a + 1);
        try {
            size_t pos = 0;
            double val = std::stod(trimmed, &pos);
            if (pos != trimmed.size()) {
                throw std::runtime_error("stringval(): trailing characters in '" + s + "'");
            }
            return ADValue::constant(val, numVariables_);
        } catch (const std::exception& e) {
            throw std::runtime_error("stringval() cannot parse '" + s + "': " + e.what());
        }
    }
    
    // CONVERT('from', 'to') — returns conversion factor so that 1 [from] = factor [to]; an unknown unit or
    // units of different dimensions are an error (the parser has already checked literal units)
    if (name == "convert" && func.args.size() == 2) {
        const double factor = unitConversionFactor(evaluateString(func.args[0]), evaluateString(func.args[1]));
        return ADValue::constant(factor, numVariables_);
    }
    
    // CONVERTTEMP('from', 'to', value) — converts a temperature value between scales (C, K, F, R); another unit is
    // an error (the parser has already checked literal units). The scales differ by a factor and an offset, so the
    // derivative with respect to the value is the ratio of the degrees.
    if (name == "converttemp" && func.args.size() == 3) {
        const std::string fromUnit = evaluateString(func.args[0]);
        const std::string toUnit = evaluateString(func.args[1]);
        ADValue result = evaluate(func.args[2]);
        double slope = 1.0;
        result.value = convertTemperature(result.value, fromUnit, toUnit, &slope);
        for (double& g : result.gradient) g *= slope;
        return result;
    }

    // UNITSYSTEM('K') — 1 if that unit setting is in use, 0 otherwise (used by EES
    // library routines to check the unit settings they were written for).
    if (name == "unitsystem" && func.args.size() == 1 && func.namedArgs.empty()) {
        const bool active = unitSettingActive(evaluateString(func.args[0]), coolpropConfig_.units);
        return ADValue::constant(active ? 1.0 : 0.0, numVariables_);
    }
    
    // ----------------------------------------------------------------
    // Lookup table / interpolation functions
    // ----------------------------------------------------------------
    // These functions have mixed arguments: strings (table name, column names, the key of LOOKUP$ROW,
    // possibly string functions) and numbers. They are evaluated from the call itself, so before the
    // arguments are all evaluated as numbers, which a string function cannot be.
    if (isLookupFunction(name)) {
        return evaluateLookupFunction(func);
    }

    std::vector<ADValue> args;
    for (const auto& arg : func.args) {
        args.push_back(evaluate(arg));
    }
    
    // Special case for pi() with 0 arguments
    if (name == "pi" && args.empty()) {
        return ADValue::constant(3.14159265358979323846, numVariables_);
    }

    return evaluateStandardFunction(func.name, args);
}

ExpressionEvaluator ExpressionEvaluator::makeBodyEvaluator(const char* kind, const std::string& name) const {
    // The depth travels with the evaluators (a body evaluator is one level deeper than its caller), so
    // there is no counter to restore when an exception unwinds.
    if (callDepth_ >= kMaxRoutineCallDepth)
        throw std::runtime_error("more than " + std::to_string(kMaxRoutineCallDepth) + " routine calls are nested when entering " +
                                 kind + " '" + name + "': a recursive FUNCTION or PROCEDURE that does not terminate? "
                                 "(recursion is not allowed in EES)");
    ExpressionEvaluator body(numVariables_, coolpropConfig_);
    body.callDepth_ = callDepth_ + 1;
    body.userFunctions_ = userFunctions_;
    body.userProcedures_ = userProcedures_;
    body.lookupTableStore_ = lookupTableStore_;
    body.diagnostics_ = diagnostics_;   // warnings raised in a body (mixtures, CoolProp errors) are reported too
    body.callWarningSink_ = callWarningSink_;   // and so are the CALL WARNING of the bodies
    body.clampLog_ = clampLog_;                 // and the inputs they limit
    body.routine_ = std::string(kind) + " '" + name + "'";
    body.routineName_ = name;
    return body;
}

ADValue ExpressionEvaluator::evaluateUserFunction(const FunctionDefinition& func, const FunctionCall& call) {
    if (call.args.size() != func.parameters.size()) {
        throw std::runtime_error("Function " + func.name + " expected " + std::to_string(func.parameters.size()) + " arguments, got " + std::to_string(call.args.size()));
    }
    
    // Create a local evaluator for the function body
    ExpressionEvaluator localEval = makeBodyEvaluator("FUNCTION", func.name);
    
    for (size_t i = 0; i < call.args.size(); ++i) {
        const std::string& paramName = func.parameters[i];
        if (!paramName.empty() && paramName.back() == '$') {
            localEval.setStringVariable(paramName, evaluateString(call.args[i]));
        } else {
            localEval.setVariable(paramName, evaluate(call.args[i]));
        }
    }
    
    // Evaluate the body procedurally
    ProceduralEvaluator::run(localEval, func.body, localEval.routine_);
    
    // The result is stored in a variable with the function name
    return localEval.getVariable(func.name);
}

/// Error of a single-input call of an ideal gas whose property depends on the pressure (density, entropy, ...).
static std::runtime_error idealGasNeedsPressure(const std::string& funcName) {
    return std::runtime_error("Ideal gas property '" + funcName + "' requires pressure input");
}

/// Error of a single-input call of an ideal gas whose input does not fix the temperature (density, specific volume,
/// entropy, quality, pressure): not even a property of the temperature alone (cp, h, u, ...) is defined by it.
static std::runtime_error idealGasNeedsSecondInput(const std::string& funcName, const std::string& fluidName,
                                                   const std::string& inputName) {
    return std::runtime_error(funcName + "(" + fluidName + ", " + inputName + "=...): an ideal gas needs the "
                              "temperature (or h or u), or two inputs such as T and P / d and P");
}

ADValue ExpressionEvaluator::evaluateCoolPropFunction(const FunctionCall& func) {
    std::string funcName = func.name;
    std::transform(funcName.begin(), funcName.end(), funcName.begin(), ::tolower);
    
    const std::string eesFluidName = fluidNameOf(func.args.empty() ? nullptr : func.args[0]);
    
    auto fluid = resolveFluid(eesFluidName, diagnostics_);
    
    // Derived thermophysical properties — compute by combining base CoolProp calls
    // entirely in AD space so gradients are correct.
    //   KINEMATICVISCOSITY = viscosity / density          (m^2/s)
    //   THERMALDIFFUSIVITY = conductivity / (density*cp)  (m^2/s)
    // For humid air the density and cp are per kg of dry air, so k/(density*cp) is the diffusivity of the mixture, but
    // viscosity/density would be off by (1+w) from its kinematic viscosity: that call is refused below as a humid-air
    // property without output, naming the function.
    if (funcName == "thermaldiffusivity" ||
        (funcName == "kinematicviscosity" && fluid->getType() != FluidType::HumidAir)) {
        // The density of an ideal gas needs its pressure, unlike the viscosity and the conductivity alone: name
        // the call that was written, not the density it divides by
        if (fluid->getType() == FluidType::IdealGas && func.namedArgs.size() == 1)
            throw idealGasNeedsPressure(funcName);
        auto makeSubCall = [&](const std::string& name) {
            FunctionCall sub = func;
            sub.name = name;
            return sub;
        };
        if (funcName == "kinematicviscosity") {
            ADValue mu = evaluateCoolPropFunction(makeSubCall("viscosity"));
            ADValue rho = evaluateCoolPropFunction(makeSubCall("density"));
            return mu / rho;
        } else {
            ADValue k = evaluateCoolPropFunction(makeSubCall("conductivity"));
            ADValue rho = evaluateCoolPropFunction(makeSubCall("density"));
            ADValue cp = evaluateCoolPropFunction(makeSubCall("cp"));
            return k / (rho * cp);
        }
    }
    
    // Extract named argument values
    std::map<std::string, ADValue> inputs;
    for (const auto& [argName, argExpr] : func.namedArgs) {
        std::string lowerArgName = argName;
        std::transform(lowerArgName.begin(), lowerArgName.end(), lowerArgName.begin(), ::tolower);
        inputs[lowerArgName] = evaluate(argExpr);
    }
    
    // Special case for molar mass (constant for pure fluids)
    if (funcName == "molarmass" || funcName == "mm") {
        if (fluid->getType() == FluidType::HumidAir) {
            // Humid air molar mass depends on humidity, but EES molarmass('Air_ha')
            // often refers to the molar mass of dry air if no state is given.
            if (inputs.empty()) {
                return ADValue::constant(28.966, numVariables_);
            }
        } else if (fluid->getType() != FluidType::Unknown) {
            // For pure fluids or incompressible, molar mass is constant
            if (inputs.empty()) {
                std::string cpFluidName = fluid->getCoolPropName();
                try {
                    if (fluid->getType() == FluidType::Mixture) {   // mole-weighted, no state needed
                        auto state = g_abstractStateCache.getOrCreate(coolpropConfig_.getBackendString(), cpFluidName);
                        return ADValue::constant(state->molar_mass() * KG_PER_KMOL, numVariables_);
                    }
                    double mm = timedPropsSI("M", "T", 300, "P", 101325, cpFluidName);
                    return ADValue::constant(mm * KG_PER_KMOL, numVariables_);
                } catch (...) {
                    // Fallback or ignore if PropsSI fails
                }
            }
        }
    }
    
    // Pure-fluid constants that depend only on the fluid itself (no state):
    //   T_CRIT, P_CRIT, V_CRIT, T_TRIPLE, P_TRIPLE, ACENTRICFACTOR
    // These are called in EES with a single positional argument (the fluid)
    // and no named arguments.
    {
        bool isConstant =
            (funcName == "t_crit" || funcName == "tcrit" ||
             funcName == "p_crit" || funcName == "pcrit" ||
             funcName == "v_crit" || funcName == "vcrit" ||
             funcName == "t_triple" || funcName == "ttriple" ||
             funcName == "p_triple" || funcName == "ptriple" ||
             funcName == "acentricfactor");
        if (isConstant && inputs.empty()) {
            if (fluid->getType() == FluidType::HumidAir) {
                throw std::runtime_error("Fluid constant '" + func.name +
                    "' is not defined for humid air (Air_ha / airH2O)");
            }
            if (fluid->getType() == FluidType::Unknown) {
                auto unsupported = std::dynamic_pointer_cast<UnsupportedFluid>(fluid);
                throw std::runtime_error("Fluid '" + eesFluidName + "' is not supported: " +
                    (unsupported ? unsupported->getReason() : "Unknown reason"));
            }
            if (fluid->getType() == FluidType::Incompressible ||
                fluid->getType() == FluidType::Mixture) {
                throw std::runtime_error("Fluid constant '" + func.name +
                    "' is not available for fluid '" + eesFluidName + "'");
            }
            
            std::string cpFluidName = fluid->getCoolPropName();
            try {
                double valueSI = 0.0;
                UnitType outType = UnitType::Dimensionless;
                std::string outUnit;
                
                if (funcName == "t_crit" || funcName == "tcrit") {
                    valueSI = CoolProp::Props1SI(cpFluidName, "Tcrit");
                    outType = UnitType::Temperature;
                    outUnit = coolpropConfig_.units.temperature;
                } else if (funcName == "p_crit" || funcName == "pcrit") {
                    valueSI = CoolProp::Props1SI(cpFluidName, "Pcrit");
                    outType = UnitType::Pressure;
                    outUnit = coolpropConfig_.units.pressure;
                } else if (funcName == "v_crit" || funcName == "vcrit") {
                    // v_crit = 1 / rhomass_critical
                    double rhoCrit = CoolProp::Props1SI(cpFluidName, "rhomass_critical");
                    valueSI = 1.0 / rhoCrit;
                    // m^3/kg — no config unit, report in SI
                    outType = UnitType::Dimensionless;
                } else if (funcName == "t_triple" || funcName == "ttriple") {
                    valueSI = CoolProp::Props1SI(cpFluidName, "Ttriple");
                    outType = UnitType::Temperature;
                    outUnit = coolpropConfig_.units.temperature;
                } else if (funcName == "p_triple" || funcName == "ptriple") {
                    valueSI = CoolProp::Props1SI(cpFluidName, "ptriple");
                    outType = UnitType::Pressure;
                    outUnit = coolpropConfig_.units.pressure;
                } else if (funcName == "acentricfactor") {
                    valueSI = CoolProp::Props1SI(cpFluidName, "acentric");
                    outType = UnitType::Dimensionless;
                }
                
                if (!std::isfinite(valueSI)) {
                    throw std::runtime_error("CoolProp returned non-finite value");
                }
                
                double result = (outType == UnitType::Dimensionless)
                    ? valueSI
                    : UnitConverter::fromSI(valueSI, outType, outUnit);
                return ADValue::constant(result, numVariables_);
            } catch (const std::exception& e) {
                throw std::runtime_error("Cannot compute '" + func.name + "' for fluid '" +
                    eesFluidName + "': " + e.what());
            }
        }
    }
    
    const UnitSystem& units = coolpropConfig_.units;
    
    // --- Humid Air Handling ---
    if (fluid->getType() == FluidType::HumidAir) {
        // Humid Air typically requires 3 inputs (e.g., T, P, W).
        // However, EES sometimes calls it with 2 inputs (e.g., T, P).
        // In such cases, we default the 3rd input to humidity ratio W=0 (dry air).
        if (inputs.size() == 2) {
            inputs["w"] = ADValue::constant(0.0, numVariables_);
        }
        
        if (inputs.size() != 3) {
            throw std::runtime_error("Humid Air properties require 2 or 3 input properties, got " + 
                                    std::to_string(inputs.size()));
        }
        
        // Check for incompatible water-content inputs (CoolProp limitation)
        // R (relative humidity) and W (humidity ratio) cannot both be specified
        bool hasR = inputs.find("r") != inputs.end();
        bool hasW = inputs.find("w") != inputs.end();
        if (hasR && hasW) {
            throw std::runtime_error(
                "CoolProp HumidAir cannot accept both R (relative humidity) and W (humidity ratio) as inputs. "
                "If you need T at R=1 with given W, use dewpoint() instead: T_c = dewpoint(airH2O, P=..., W=..., T=T_ref)");
        }
        
        std::vector<std::string> inputNamesList;
        std::vector<ADValue> inputValues;
        for (const auto& [name, val] : inputs) {
            inputNamesList.push_back(name);
            inputValues.push_back(val);
        }
        
        HAParam outParam = getHAOutputParam(funcName);
        std::vector<HAParam> inParams;
        for (const auto& name : inputNamesList) {
            inParams.push_back(getHAParam(name));
        }
        
        std::string outputStr = haParamToString(outParam);
        std::vector<std::string> inputStrs;
        for (auto p : inParams) {
            inputStrs.push_back(haParamToString(p));
        }
        
        std::vector<double> vals;
        for (size_t i = 0; i < 3; ++i) {
            double v = inputValues[i].value;
            checkThermoInputValueHint(diagnostics_, func.name + "()", units,
                                      thermoKindFromHAName(inputNamesList[i]),
                                      inputNamesList[i], v);
            UnitType type = haParamToUnitType(inParams[i]);
            
            if (type != UnitType::Dimensionless) {
                v = UnitConverter::toSI(v, type, 
                    type == UnitType::Temperature ? units.temperature : 
                    type == UnitType::Pressure ? units.pressure :
                    type == UnitType::SpecificEnergy ? units.specific_energy :
                    type == UnitType::SpecificEntropy ? units.specific_entropy :
                    type == UnitType::Density ? units.density : "");
            }
            vals.push_back(v);
        }
        
        // The output at the inputs `v`; the density is 1/Vda, so its derivatives below are those of 1/Vda as well
        // (a failure of CoolProp stays NaN/Inf: 1/Inf must not become a density of 0).
        auto haOutput = [&](const std::vector<double>& v) {
            const double r = timedHAPropsSI(outputStr, inputStrs[0], v[0], inputStrs[1], v[1], inputStrs[2], v[2]);
            return (outParam == HAParam::Rho && std::isfinite(r)) ? 1.0 / r : r;
        };

        try {
            double result = haOutput(vals);
            
            if (!std::isfinite(result)) {
                // Special handling for relative humidity with supersaturated air:
                // when W > W_sat at the given T and P, HAPropsSI("R",...) returns
                // NaN.  Compute RH = W / W_sat instead (gives RH > 1, which the
                // model can use to detect condensation).
                bool handled = false;
                if (outParam == HAParam::R) {
                    // Find the W and non-W inputs
                    int wIdx = -1;
                    for (size_t k = 0; k < 3; ++k) {
                        if (inParams[k] == HAParam::W) { wIdx = k; break; }
                    }
                    if (wIdx >= 0) {
                        // Build a call with R=1 instead of W to get W_sat
                        std::vector<std::string> satStrs;
                        std::vector<double> satVals;
                        for (size_t k = 0; k < 3; ++k) {
                            if ((int)k == wIdx) {
                                satStrs.push_back("R");
                                satVals.push_back(1.0);
                            } else {
                                satStrs.push_back(inputStrs[k]);
                                satVals.push_back(vals[k]);
                            }
                        }
                        double wSat = timedHAPropsSI("W", satStrs[0], satVals[0],
                                                     satStrs[1], satVals[1],
                                                     satStrs[2], satVals[2]);
                        if (std::isfinite(wSat) && wSat > 0) {
                            result = vals[wIdx] / wSat;
                            handled = true;
                        }
                    }
                }
                if (!handled) {
                    // CoolProp returns NaN/Inf instead of throwing (inputs outside its validity range, output it has
                    // no name for). Reported like every other CoolProp failure: C001 below, with the penalty value
                    // kept for the iterations, and the solution check turns it into an error.
                    if (outParam == HAParam::Unknown)
                        throw std::runtime_error("'" + func.name + "' is not available for humid air (" + eesFluidName +
                                                 "): CoolSolve has no CoolProp humid-air output for it");
                    throw std::runtime_error("CoolProp returned invalid result (NaN or Inf) for " + func.name + "(" +
                                             eesFluidName + ") " + humidAirFailureDetail(inParams, vals, units));
                }
            }
            
            // Convert output from SI
            UnitType outType = haParamToUnitType(outParam);
            if (outType != UnitType::Dimensionless) {
                result = UnitConverter::fromSI(result, outType, 
                    outType == UnitType::Temperature ? units.temperature : 
                    outType == UnitType::Pressure ? units.pressure :
                    outType == UnitType::SpecificEnergy ? units.specific_energy :
                    outType == UnitType::SpecificEntropy ? units.specific_entropy :
                    outType == UnitType::Density ? units.density : "");
            }
            
            // Derivatives
            const double relEps = 1e-6;
            const double absEps = 1e-8;
            std::vector<double> dResult_dInput(3, 0.0);
            
            for (size_t i = 0; i < 3; ++i) {
                double step = std::max(relEps * std::abs(vals[i]), absEps);
                std::vector<double> vPlus = vals;
                std::vector<double> vMinus = vals;
                vPlus[i] += step;
                vMinus[i] -= step;
                
                double rPlus = haOutput(vPlus);
                double rMinus = haOutput(vMinus);
                
                if (std::isfinite(rPlus) && std::isfinite(rMinus)) {
                    dResult_dInput[i] = (rPlus - rMinus) / (2.0 * step);
                }
            }
            
            // Scaling factors for derivatives
            auto getScale = [&](UnitType t, const std::string& u, bool toSI) {
                if (t == UnitType::Dimensionless) return 1.0;
                if (toSI) return UnitConverter::toSI(1.0, t, u) - UnitConverter::toSI(0.0, t, u);
                return UnitConverter::fromSI(1.0, t, u) - UnitConverter::fromSI(0.0, t, u);
            };
            
            double outScale = getScale(outType, 
                outType == UnitType::Temperature ? units.temperature : 
                outType == UnitType::Pressure ? units.pressure :
                outType == UnitType::SpecificEnergy ? units.specific_energy :
                outType == UnitType::SpecificEntropy ? units.specific_entropy :
                outType == UnitType::Density ? units.density : "", false);
                
            ADValue output(result, numVariables_);
            for (size_t j = 0; j < 3; ++j) {
                UnitType inType = haParamToUnitType(inParams[j]);
                double inScale = getScale(inType, 
                    inType == UnitType::Temperature ? units.temperature : 
                    inType == UnitType::Pressure ? units.pressure :
                    inType == UnitType::SpecificEnergy ? units.specific_energy :
                    inType == UnitType::SpecificEntropy ? units.specific_entropy :
                    inType == UnitType::Density ? units.density : "", true);
                    
                double dTotal = dResult_dInput[j] * outScale * inScale;
                for (size_t i = 0; i < numVariables_; ++i) {
                    output.gradient[i] += dTotal * inputValues[j].gradient[i];
                }
            }
            
            return output;
        } catch (const std::exception& e) {
            // Penalty-based HumidAir error handling (consistent with PropsSI)
            if (diagnostics_) {
                diagnostics_->push(DiagnosticSeverity::Warning, "C001",
                    "HumidAir error in " + func.name + "(): " + e.what(),
                    "evaluator");
            }
            constexpr double PENALTY = 1e4;
            ADValue output(PENALTY, numVariables_);
            return output;
        }
    }
    
    // --- Standard CoolProp Handling ---
    if (fluid->getType() == FluidType::Unknown) {
        auto unsupported = std::dynamic_pointer_cast<UnsupportedFluid>(fluid);
        throw std::runtime_error("Fluid '" + eesFluidName + "' is not supported: " + (unsupported ? unsupported->getReason() : "Unknown reason"));
    }
    
    // Incompressible fluid: accept T,P inputs; reject saturation-style inputs.
    // For solutions, resolve concentration from a second positional argument
    // (mass fraction in %, e.g. density('MEG', 30, T=280, P=1e5)) or from a
    // `C` named argument.
    std::string cpFluidName = fluid->getCoolPropName();
    if (fluid->getType() == FluidType::Incompressible) {
        auto incomp = std::dynamic_pointer_cast<IncompressibleFluid>(fluid);
        if (incomp && incomp->isSolution()) {
            double massFraction = -1.0;
            // 2nd positional argument (numeric literal or any expression)
            if (func.args.size() >= 2) {
                try {
                    ADValue v = evaluate(func.args[1]);
                    // Interpret as mass percent 0..100 → fraction 0..1
                    massFraction = v.value / 100.0;
                } catch (...) {
                    // Not an expression — ignore
                }
            }
            // C or concentration named arg (mass percent)
            for (auto it = inputs.begin(); it != inputs.end(); ) {
                if (it->first == "c" || it->first == "conc" || it->first == "concentration") {
                    massFraction = it->second.value / 100.0;
                    it = inputs.erase(it);
                } else {
                    ++it;
                }
            }
            if (massFraction < 0.0 || massFraction > 1.0) {
                throw std::runtime_error("Incompressible solution '" + eesFluidName +
                    "' requires a concentration (mass %, 0..100) as a second positional "
                    "argument or a C=... named argument.");
            }
            cpFluidName = incomp->getCoolPropNameWithFraction(massFraction);
        }
        // Reject quality-based inputs — saturation is not defined for incompressibles.
        for (const auto& [k, v] : inputs) {
            if (k == "q" || k == "x" || k == "quality") {
                throw std::runtime_error("Incompressible fluid '" + eesFluidName +
                    "' does not support saturation (quality) inputs.");
            }
        }
        // Reject request for quality output for incompressibles.
        if (funcName == "quality" || funcName == "x") {
            throw std::runtime_error("Incompressible fluid '" + eesFluidName +
                "' has no vapor quality.");
        }
    }
    
    // True when the pressure was injected below (single-input ideal-gas call): it is
    // an internal constant of the ideal-gas model, not a user input, so it must be the
    // same in every pass — never clamped, never checked against the unit hints.
    bool pressureInjected = false;

    if (fluid->getType() == FluidType::IdealGas && inputs.size() == 1) {
        // The properties of an ideal gas that depend only on temperature {T, h, u,
        // cp, cv} can be computed from any other member of this set without
        // knowing pressure.  When exactly one input is given, it must fix the
        // temperature (T, h or u) and the requested output must belong to this
        // "thermal-only" set: a dummy pressure is then injected so CoolProp gets
        // its required 2-input state.  Any other single input (density, specific
        // volume, entropy, quality, pressure) does not fix the state, whatever the
        // output: a second input is required, the dummy pressure is never a state.
        // The transport properties of the dilute gas (viscosity, conductivity, Prandtl
        // number) depend on the temperature only too: they are outputs of this set,
        // evaluated at the dummy pressure (a few tenths of a percent from the
        // dilute-gas limit); the density-dependent ones (kinematic viscosity, thermal
        // diffusivity) are not, they need the pressure.
        auto isThermalOnly = [](const std::string& prop) {
            std::string p = prop;
            for (auto& c : p) c = std::tolower(c);
            return p == "h" || p == "enthalpy" || p == "u" || p == "internalenergy" ||
                   p == "cp" || p == "cv" || p == "specheat" ||
                   p == "t" || p == "temperature";
        };
        auto isTransport = [](const std::string& prop) {
            return prop == "viscosity" || prop == "mu" || prop == "conductivity" ||
                   prop == "k" || prop == "l" || prop == "prandtl";
        };
        auto fixesTemperature = [](const std::string& prop) {
            std::string p = prop;
            for (auto& c : p) c = std::tolower(c);
            return p == "t" || p == "temperature" || p == "h" || p == "enthalpy" ||
                   p == "u" || p == "internalenergy" || p == "intenergy";
        };
        
        std::string inputName = inputs.begin()->first;
        (void)getCoolPropInputParam(inputName);   // an unknown input name keeps its own error
        if (!fixesTemperature(inputName)) {
            throw idealGasNeedsSecondInput(func.name, eesFluidName, inputName);
        } else if (isThermalOnly(funcName) || isTransport(funcName)) {
            // The input fixes T and the output depends on T only → inject dummy pressure
            // Use fluid-specific dummy pressure (e.g. 100 Pa for H2O to avoid
            // liquid phase at room temperature).
            auto idealGas = std::dynamic_pointer_cast<IdealGasFluid>(fluid);
            double dummyP = idealGas ? idealGas->getDummyPressureSI() : 101325.0;
            inputs["p"] = ADValue::constant(UnitConverter::fromSI(dummyP, UnitType::Pressure, units.pressure), numVariables_);
            pressureInjected = true;
        } else {
            throw idealGasNeedsPressure(funcName);
        }
    }
    
    // Handle saturation functions (T_sat, P_sat, etc.) which take only 1 property input
    // We add quality Q=0.5 to force saturation state in CoolProp
    if (inputs.size() == 1) {
        std::string lowerName = funcName;
        std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
        if (lowerName == "t_sat" || lowerName == "p_sat" || lowerName == "tsat" || lowerName == "psat") {
            if (fluid->getType() == FluidType::Mixture) {
                throw std::runtime_error(func.name + "() is not defined for the mixture '" + fluid->getName() +
                    "': it has a bubble and a dew point, not one saturation point. Use temperature(fluid, P=..., x=0) "
                    "for the bubble point and x=1 for the dew point (pressure(fluid, T=..., x=0 or 1) at given T)");
            }
            inputs["q"] = ADValue::constant(0.5, numVariables_);
        }
    }
    
    if (inputs.size() != 2) {
        throw std::runtime_error(func.name + " requires exactly 2 input properties "
            "(e.g. T=..., P=...), got " + std::to_string(inputs.size()) +
            coolPropUnitsReminder(units));
    }
    
    std::vector<std::string> inputNamesList;
    std::vector<ADValue> inputValues;
    for (const auto& [name, val] : inputs) {
        inputNamesList.push_back(name);
        inputValues.push_back(val);
    }
    
    CoolPropParamInfo outputInfo = getCoolPropOutputParam(funcName);
    CoolProp::parameters input1Param = getCoolPropInputParam(inputNamesList[0]);
    CoolProp::parameters input2Param = getCoolPropInputParam(inputNamesList[1]);
    
    std::string outputStr = paramToString(outputInfo.param);
    std::string input1Str = paramToString(input1Param);
    std::string input2Str = paramToString(input2Param);
    
    // Check for unsupported input pairs
    // CoolProp does not support T,H as an input pair for any output property
    if ((input1Param == CoolProp::iT && input2Param == CoolProp::iHmass) ||
        (input1Param == CoolProp::iHmass && input2Param == CoolProp::iT)) {
        throw std::runtime_error("CoolProp does not support T,H as input pair. Use P,H or T,S instead." +
            coolPropUnitsReminder(units));
    }
    
    double val1 = inputValues[0].value;
    double val2 = inputValues[1].value;

    // The injected ideal-gas pressure is not a user input: no unit hint for it.
    if (!(pressureInjected && input1Param == CoolProp::iP))
        checkThermoInputValueHint(diagnostics_, func.name + "()", units,
                                  thermoKindFromCoolProp(input1Param, inputNamesList[0]),
                                  inputNamesList[0], val1);
    if (!(pressureInjected && input2Param == CoolProp::iP))
        checkThermoInputValueHint(diagnostics_, func.name + "()", units,
                                  thermoKindFromCoolProp(input2Param, inputNamesList[1]),
                                  inputNamesList[1], val2);
    
    // Clamp quality inputs to [0, 1]: CoolProp returns -1 for subcooled and values >1 for
    // superheated states, but when quality is used as *input* to another CoolProp call
    // (e.g. enthalpy(Water, P=P, x=Q)), out-of-range values cause NaN results.
    // Clamping is disabled during final solution verification.
    if (!disableClamping_) {
    if (input1Param == CoolProp::iQ) {
        if (val1 < 0.0) {
            if (diagnostics_) {
                diagnostics_->push(DiagnosticSeverity::Warning, "C002",
                    "Quality input " + inputNamesList[0] + "=" + std::to_string(val1) +
                        " clamped to 0 (subcooled); quality is only valid in [0, 1] for two-phase states.",
                    "evaluator");
            }
            val1 = 0.0;
            for (size_t i = 0; i < numVariables_; ++i) inputValues[0].gradient[i] = 0.0;
        } else if (val1 > 1.0) {
            if (diagnostics_) {
                diagnostics_->push(DiagnosticSeverity::Warning, "C002",
                    "Quality input " + inputNamesList[0] + "=" + std::to_string(val1) +
                        " clamped to 1 (superheated); quality is only valid in [0, 1] for two-phase states.",
                    "evaluator");
            }
            val1 = 1.0;
            for (size_t i = 0; i < numVariables_; ++i) inputValues[0].gradient[i] = 0.0;
        }
    }
    if (input2Param == CoolProp::iQ) {
        if (val2 < 0.0) {
            if (diagnostics_) {
                diagnostics_->push(DiagnosticSeverity::Warning, "C002",
                    "Quality input " + inputNamesList[1] + "=" + std::to_string(val2) +
                        " clamped to 0 (subcooled); quality is only valid in [0, 1] for two-phase states.",
                    "evaluator");
            }
            val2 = 0.0;
            for (size_t i = 0; i < numVariables_; ++i) inputValues[1].gradient[i] = 0.0;
        } else if (val2 > 1.0) {
            if (diagnostics_) {
                diagnostics_->push(DiagnosticSeverity::Warning, "C002",
                    "Quality input " + inputNamesList[1] + "=" + std::to_string(val2) +
                        " clamped to 1 (superheated); quality is only valid in [0, 1] for two-phase states.",
                    "evaluator");
            }
            val2 = 1.0;
            for (size_t i = 0; i < numVariables_; ++i) inputValues[1].gradient[i] = 0.0;
        }
    }
    } // end of quality clamping
    
    UnitType type1 = UnitType::Dimensionless;
    UnitType type2 = UnitType::Dimensionless;
    
    if (input1Param == CoolProp::iT) type1 = UnitType::Temperature;
    else if (input1Param == CoolProp::iP) type1 = UnitType::Pressure;
    else if (input1Param == CoolProp::iHmass) type1 = UnitType::SpecificEnergy;
    else if (input1Param == CoolProp::iSmass) type1 = UnitType::SpecificEntropy;
    else if (input1Param == CoolProp::iDmass) {
        // If the input name was 'v' or 'volume', the value is specific volume (1/rho)
        std::string inputName = inputNamesList[0];
        std::transform(inputName.begin(), inputName.end(), inputName.begin(), ::tolower);
        if (inputName == "v" || inputName == "volume") {
            type1 = UnitType::Density; // We'll treat it as density but invert it
            val1 = 1.0 / val1;
            // Adjust gradient for 1/v: d(1/v)/dx = -1/v^2 * dv/dx
            double invV2 = -1.0 / (inputValues[0].value * inputValues[0].value);
            for (size_t i = 0; i < numVariables_; ++i) {
                inputValues[0].gradient[i] *= invV2;
            }
        } else {
            type1 = UnitType::Density;
        }
    }
    
    if (input2Param == CoolProp::iT) type2 = UnitType::Temperature;
    else if (input2Param == CoolProp::iP) type2 = UnitType::Pressure;
    else if (input2Param == CoolProp::iHmass) type2 = UnitType::SpecificEnergy;
    else if (input2Param == CoolProp::iSmass) type2 = UnitType::SpecificEntropy;
    else if (input2Param == CoolProp::iDmass) {
        // If the input name was 'v' or 'volume', the value is specific volume (1/rho)
        std::string inputName = inputNamesList[1];
        std::transform(inputName.begin(), inputName.end(), inputName.begin(), ::tolower);
        if (inputName == "v" || inputName == "volume") {
            type2 = UnitType::Density;
            val2 = 1.0 / val2;
            // Adjust gradient for 1/v: d(1/v)/dx = -1/v^2 * dv/dx
            double invV2 = -1.0 / (inputValues[1].value * inputValues[1].value);
            for (size_t i = 0; i < numVariables_; ++i) {
                inputValues[1].gradient[i] *= invV2;
            }
        } else {
            type2 = UnitType::Density;
        }
    }
    
    val1 = UnitConverter::toSI(val1, type1, 
        type1 == UnitType::Temperature ? units.temperature : 
        type1 == UnitType::Pressure ? units.pressure :
        type1 == UnitType::SpecificEnergy ? units.specific_energy :
        type1 == UnitType::SpecificEntropy ? units.specific_entropy :
        type1 == UnitType::Density ? units.density : "");
        
    val2 = UnitConverter::toSI(val2, type2, 
        type2 == UnitType::Temperature ? units.temperature : 
        type2 == UnitType::Pressure ? units.pressure :
        type2 == UnitType::SpecificEnergy ? units.specific_energy :
        type2 == UnitType::SpecificEntropy ? units.specific_entropy :
        type2 == UnitType::Density ? units.density : "");
    
    // Fluids with an EES reference state (ideal-gas enthalpy of formation): enthalpy-type
    // inputs are given on the EES basis, CoolProp works on its own — shift them back.
    // The same offset is added to the results (postProcess below).
    const ReferenceState refState = fluid->getReferenceState(coolpropConfig_.getBackendString());
    auto toCoolPropBasis = [&](CoolProp::parameters param, double& v) {
        if (param == CoolProp::iHmass || param == CoolProp::iUmass) v -= refState.h_offset;
        else if (param == CoolProp::iSmass) v -= refState.s_offset;
    };
    toCoolPropBasis(input1Param, val1);
    toCoolPropBasis(input2Param, val2);

    // ── CoolProp input sanitization ──────────────────────────────────
    // During Newton iteration, the solver may evaluate trial points with
    // unphysical thermodynamic inputs (negative pressure, temperature below
    // absolute zero).  Instead of letting CoolProp throw (which causes the
    // line-search / trust-region to shrink the step but may leave no viable
    // step direction), we clamp the inputs to physically valid ranges and
    // let the solver see a finite (but poor) residual, keeping the
    // optimization landscape smooth and allowing it to navigate back.
    //
    // Only an input that depends on the unknowns of the block is a trial point
    // (its AD gradient is not zero). An input that does not (a constant, a
    // variable solved by an earlier block, the injected pressure of an
    // ideal-gas call) is a value of the model: it reaches CoolProp unchanged,
    // here and in the final verification, which evaluates the same call.
    // The lower bound of a trial point is the one below, or the minimum of the
    // fluid in CoolProp when that is lower (Hydrogen, Helium: Tmin < 50 K).
    // A clamp that is still active at the converged point is refused by
    // SystemEvaluator::checkInputClamps.
    // Clamping is disabled during final solution verification so that
    // out-of-range inputs are caught as errors.
    if (!disableClamping_) {
    constexpr double P_MIN_SI  = 1000.0;   // 1000 Pa (0.01 bar) — low but CoolProp-safe
    constexpr double T_MIN_SI  = 50.0;     // 50 K   — below triple point of most fluids
    constexpr double RHO_MIN_SI = 1e-4;    // near-vacuum density floor

    auto clampInput = [&](size_t index, CoolProp::parameters param, double& v) {
        double bound;
        if (param == CoolProp::iP && v < P_MIN_SI) bound = P_MIN_SI;
        else if (param == CoolProp::iT && v < T_MIN_SI) bound = T_MIN_SI;
        else if (param == CoolProp::iDmass && v < RHO_MIN_SI) bound = RHO_MIN_SI;
        else return;   // every evaluation of a value above the bounds ends here
        // Below a bound: the rare case, the one that is worth looking at the gradient and the fluid.
        const auto& grad = inputValues[index].gradient;
        if (std::none_of(grad.begin(), grad.end(), [](double g) { return g != 0.0; })) return;
        const char* quantity = param == CoolProp::iP ? "pressure" : param == CoolProp::iT ? "temperature" : "density";
        const char* unit = param == CoolProp::iP ? " Pa" : param == CoolProp::iT ? " K" : " kg/m³";
        std::string source = std::string("the lowest ") + quantity + " the solver evaluates";
        if (param != CoolProp::iDmass) {
            const FluidLimits limits = g_abstractStateCache.fluidLimits(coolpropConfig_.getBackendString(), cpFluidName);
            const double fluidLimit = (param == CoolProp::iP) ? limits.pMin : limits.Tmin;
            if (fluidLimit > 0.0 && fluidLimit < bound) {
                bound = fluidLimit;
                source = std::string("the minimum ") + quantity + " of " + eesFluidName + " in CoolProp" +
                         (param == CoolProp::iP ? ", its triple point" : "");
            }
            if (!(v < bound)) return;
        }
        if (clampLog_) {
            clampLog_->fired = true;
            if (clampLog_->report) {
                std::ostringstream oss;
                oss << func.name << "(" << eesFluidName;
                for (const auto& named : func.namedArgs) oss << ", " << named.first;
                oss << ")" << (routine_.empty() ? "" : " in " + routine_) << ": " << quantity << " " << v << unit
                    << " is below " << source << " (" << bound << unit
                    << "); an input that depends on the unknowns of its block is limited there during the iterations,"
                       " so the solution would not satisfy this equation";
                clampLog_->report->push_back(oss.str());
            }
        }
        // The gradient is kept: the solver still has a direction back to the valid region
        v = bound;
    };

    clampInput(0, input1Param, val1);
    clampInput(1, input2Param, val2);
    // ─────────────────────────────────────────────────────────────────
    } // end of clamping block
    
    // Helper for unit scale factors (reused by both paths)
    auto getScale = [&](UnitType t, const std::string& u, bool toSI) {
        if (t == UnitType::Dimensionless) return 1.0;
        if (toSI) return UnitConverter::toSI(1.0, t, u) - UnitConverter::toSI(0.0, t, u);
        return UnitConverter::fromSI(1.0, t, u) - UnitConverter::fromSI(0.0, t, u);
    };
    auto unitStringFor = [&](UnitType t) -> std::string {
        if (t == UnitType::Temperature) return units.temperature;
        if (t == UnitType::Pressure) return units.pressure;
        if (t == UnitType::SpecificEnergy) return units.specific_energy;
        if (t == UnitType::SpecificEntropy) return units.specific_entropy;
        if (t == UnitType::Density) return units.density;
        if (t == UnitType::Viscosity) return units.viscosity;
        if (t == UnitType::Conductivity) return units.conductivity;
        if (t == UnitType::SpecificHeat) return units.specific_heat;
        return "";
    };
    
    // SI unit label for CoolProp parameters (used in error messages)
    auto siUnitLabel = [](CoolProp::parameters param) -> std::string {
        switch (param) {
            case CoolProp::iT: return " K";
            case CoolProp::iP: return " Pa";
            case CoolProp::iHmass: return " J/kg";
            case CoolProp::iSmass: return " J/(kg·K)";
            case CoolProp::iDmass: return " kg/m³";
            case CoolProp::iQ: return "";
            default: return "";
        }
    };
    
    // Common post-processing: unit conversion, volume inversion, gradient assembly
    auto postProcess = [&](double result, double dResult_dInput1, double dResult_dInput2) -> ADValue {
        // Clamp quality output
        if (outputInfo.param == CoolProp::iQ) {
            if (result < 0.0) result = 0.0;
            else if (result > 1.0) result = 1.0;
        }
        
        // Unit conversion for result
        if (outputInfo.unitType != UnitType::Dimensionless) {
            result = UnitConverter::fromSI(result, outputInfo.unitType, unitStringFor(outputInfo.unitType));
        }
        // Reference state offsets
        if (outputInfo.unitType == UnitType::SpecificEnergy)
            result += UnitConverter::fromSI(refState.h_offset, outputInfo.unitType, units.specific_energy);
        else if (outputInfo.unitType == UnitType::SpecificEntropy)
            result += UnitConverter::fromSI(refState.s_offset, outputInfo.unitType, units.specific_entropy);
        
        // EES returns the molar mass in kg/kmol, CoolProp in kg/mol
        if (outputInfo.param == CoolProp::imolar_mass) {
            result *= KG_PER_KMOL;
            dResult_dInput1 *= KG_PER_KMOL;
            dResult_dInput2 *= KG_PER_KMOL;
        }

        // Volume inversion (density→specific volume)
        std::string lowerFuncName = func.name;
        std::transform(lowerFuncName.begin(), lowerFuncName.end(), lowerFuncName.begin(), ::tolower);
        if (lowerFuncName == "volume") {
            double densitySI = UnitConverter::toSI(result, UnitType::Density, units.density);
            double invDensity2 = 1.0 / (densitySI * densitySI);
            result = 1.0 / result;
            dResult_dInput1 *= -invDensity2;
            dResult_dInput2 *= -invDensity2;
        }
        
        // Apply unit scale factors to derivatives
        double outScale = getScale(outputInfo.unitType, unitStringFor(outputInfo.unitType), false);
        double in1Scale = getScale(type1, unitStringFor(type1), true);
        double in2Scale = getScale(type2, unitStringFor(type2), true);
        dResult_dInput1 *= outScale * in1Scale;
        dResult_dInput2 *= outScale * in2Scale;
        
        // Assemble ADValue
        ADValue output(result, numVariables_);
        for (size_t i = 0; i < numVariables_; ++i) {
            output.gradient[i] = dResult_dInput1 * inputValues[0].gradient[i] + 
                                 dResult_dInput2 * inputValues[1].gradient[i];
        }
        return output;
    };
    
    try {
        double result = 0.0;
        double dResult_dInput1 = 0.0;
        double dResult_dInput2 = 0.0;
        
        // ── AbstractState path (low-level API) ──────────────────────
        bool usedAbstractState = false;
        if (coolpropConfig_.useAbstractState) {
            try {
                auto startAS = std::chrono::high_resolution_clock::now();
                
                auto state = g_abstractStateCache.getOrCreate(
                    coolpropConfig_.getBackendString(), cpFluidName);
                
                // Determine input pair and canonical ordering
                double iv1 = 0, iv2 = 0;
                CoolProp::input_pairs ipair = CoolProp::generate_update_pair(
                    input1Param, val1, input2Param, val2, iv1, iv2);
                
                state->update(ipair, iv1, iv2);
                result = getOutputValue(*state, outputInfo.param);
                
                auto endAS = std::chrono::high_resolution_clock::now();
                g_profilingStats.abstractState_calls++;
                g_profilingStats.abstractState_time_ms +=
                    std::chrono::duration<double, std::milli>(endAS - startAS).count();
                
                if (!std::isfinite(result)) {
                    std::ostringstream oss;
                    oss << "CoolProp returned invalid result (NaN or Inf) for " << outputStr 
                        << "(" << cpFluidName << ") with inputs: " 
                        << input1Str << "=" << val1 << siUnitLabel(input1Param) << ", "
                        << input2Str << "=" << val2 << siUnitLabel(input2Param)
                        << coolPropUnitsReminder(units);
                    throw std::runtime_error(oss.str());
                }
                
                // ── Derivatives ──────────────────────────────────────
                if (!residualOnly_) {
                    bool derivOk = false;
                    
                    if (coolpropConfig_.enableAnalyticalDerivatives) {
                        // Analytical derivatives via first_partial_deriv() with
                        // forward-FD consistency check.  CoolProp's analytical
                        // derivatives can be wildly wrong near phase boundaries
                        // for pseudo-pure / mixture fluids (e.g. Air).  We
                        // validate by comparing against a quick forward-FD and
                        // fall back to the forward-FD values when they disagree.
                        try {
                            auto dStart = std::chrono::high_resolution_clock::now();
                            double ad1 = state->first_partial_deriv(
                                outputInfo.param, input1Param, input2Param);
                            double ad2 = state->first_partial_deriv(
                                outputInfo.param, input2Param, input1Param);
                            auto dEnd = std::chrono::high_resolution_clock::now();
                            g_profilingStats.analyticalDeriv_calls++;
                            g_profilingStats.abstractState_time_ms +=
                                std::chrono::duration<double, std::milli>(dEnd - dStart).count();
                            
                            if (!std::isfinite(ad1) || !std::isfinite(ad2)) {
                                throw std::runtime_error("Non-finite analytical derivative");
                            }
                            
                            // Forward-FD consistency check (2 extra state evaluations)
                            const double relEps = 1e-6;
                            const double absEps = 1e-8;
                            const double consistencyTol = 0.05; // 5% relative tolerance
                            
                            double step1 = std::max(relEps * std::abs(val1), absEps);
                            double step2 = std::max(relEps * std::abs(val2), absEps);
                            
                            double p1v1, p1v2, p2v1, p2v2;
                            CoolProp::generate_update_pair(
                                input1Param, val1 + step1, input2Param, val2, p1v1, p1v2);
                            state->update(ipair, p1v1, p1v2);
                            double fd1 = (getOutputValue(*state, outputInfo.param) - result) / step1;
                            
                            CoolProp::generate_update_pair(
                                input1Param, val1, input2Param, val2 + step2, p2v1, p2v2);
                            state->update(ipair, p2v1, p2v2);
                            double fd2 = (getOutputValue(*state, outputInfo.param) - result) / step2;
                            
                            // Check consistency for both inputs
                            bool consistent = true;
                            double maxAbs1 = std::max(std::abs(ad1), std::abs(fd1));
                            if (maxAbs1 > 1e-15 &&
                                std::abs(ad1 - fd1) / maxAbs1 > consistencyTol)
                                consistent = false;
                            double maxAbs2 = std::max(std::abs(ad2), std::abs(fd2));
                            if (maxAbs2 > 1e-15 &&
                                std::abs(ad2 - fd2) / maxAbs2 > consistencyTol)
                                consistent = false;
                            
                            if (consistent) {
                                // Analytical derivatives are accurate — use them
                                dResult_dInput1 = ad1;
                                dResult_dInput2 = ad2;
                            } else {
                                // Analytical derivatives inaccurate (likely near
                                // phase boundary) — use forward-FD instead
                                dResult_dInput1 = std::isfinite(fd1) ? fd1 : 0.0;
                                dResult_dInput2 = std::isfinite(fd2) ? fd2 : 0.0;
                            }
                            derivOk = true;
                        } catch (...) {
                            // first_partial_deriv threw — fall through to FD
                        }
                    }
                    
                    if (!derivOk) {
                        // FD derivatives via AbstractState (central differences)
                        const double relEps = 1e-6;
                        const double absEps = 1e-8;
                        double step1 = std::max(relEps * std::abs(val1), absEps);
                        double step2 = std::max(relEps * std::abs(val2), absEps);
                        
                        try {
                            double p1v1, p1v2, m1v1, m1v2;
                            CoolProp::generate_update_pair(input1Param, val1 + step1, input2Param, val2, p1v1, p1v2);
                            state->update(ipair, p1v1, p1v2);
                            double rp1 = getOutputValue(*state, outputInfo.param);
                            CoolProp::generate_update_pair(input1Param, val1 - step1, input2Param, val2, m1v1, m1v2);
                            state->update(ipair, m1v1, m1v2);
                            double rm1 = getOutputValue(*state, outputInfo.param);
                            dResult_dInput1 = (std::isfinite(rp1) && std::isfinite(rm1)) ? (rp1 - rm1) / (2.0 * step1) : 0.0;
                        } catch (...) {
                            dResult_dInput1 = 0.0;
                        }
                        try {
                            double p2v1, p2v2, m2v1, m2v2;
                            CoolProp::generate_update_pair(input1Param, val1, input2Param, val2 + step2, p2v1, p2v2);
                            state->update(ipair, p2v1, p2v2);
                            double rp2 = getOutputValue(*state, outputInfo.param);
                            CoolProp::generate_update_pair(input1Param, val1, input2Param, val2 - step2, m2v1, m2v2);
                            state->update(ipair, m2v1, m2v2);
                            double rm2 = getOutputValue(*state, outputInfo.param);
                            dResult_dInput2 = (std::isfinite(rp2) && std::isfinite(rm2)) ? (rp2 - rm2) / (2.0 * step2) : 0.0;
                        } catch (...) {
                            dResult_dInput2 = 0.0;
                        }
                    }
                }
                // else: residualOnly_ → derivatives stay 0
                
                usedAbstractState = true;
            } catch (const std::exception& e) {
                // A mixture has no other path (PropsSI would only return Inf and lose the reason)
                if (fluid->getType() == FluidType::Mixture) {
                    const std::string what = e.what();
                    if (what.find("with inputs:") != std::string::npos) throw;   // already says which call
                    std::ostringstream oss;
                    oss << what << " [" << outputStr << "(" << cpFluidName << ") with inputs: "
                        << input1Str << "=" << val1 << siUnitLabel(input1Param) << ", "
                        << input2Str << "=" << val2 << siUnitLabel(input2Param) << "]";
                    throw std::runtime_error(oss.str());
                }
                // AbstractState failed — fall through to PropsSI below
                usedAbstractState = false;
            } catch (...) {
                usedAbstractState = false;
            }
        }
        
        // ── PropsSI fallback path ────────────────────────────────────
        if (!usedAbstractState) {
            result = timedPropsSI(outputStr, input1Str, val1, input2Str, val2, cpFluidName);
            
            if (!std::isfinite(result)) {
                std::ostringstream oss;
                oss << "CoolProp returned invalid result (NaN or Inf) for " << outputStr 
                    << "(" << cpFluidName << ") with inputs: " 
                    << input1Str << "=" << val1 << siUnitLabel(input1Param) << ", "
                    << input2Str << "=" << val2 << siUnitLabel(input2Param)
                    << coolPropUnitsReminder(units);
                throw std::runtime_error(oss.str());
            }
            
            if (!residualOnly_) {
                const double relEps = 1e-6;
                const double absEps = 1e-8;
                double step1 = std::max(relEps * std::abs(val1), absEps);
                double step2 = std::max(relEps * std::abs(val2), absEps);
                
                double resultPlus1 = timedPropsSI(outputStr, input1Str, val1 + step1, input2Str, val2, cpFluidName);
                double resultMinus1 = timedPropsSI(outputStr, input1Str, val1 - step1, input2Str, val2, cpFluidName);
                dResult_dInput1 = (std::isfinite(resultPlus1) && std::isfinite(resultMinus1)) ? (resultPlus1 - resultMinus1) / (2.0 * step1) : 0.0;
                
                double resultPlus2 = timedPropsSI(outputStr, input1Str, val1, input2Str, val2 + step2, cpFluidName);
                double resultMinus2 = timedPropsSI(outputStr, input1Str, val1, input2Str, val2 - step2, cpFluidName);
                dResult_dInput2 = (std::isfinite(resultPlus2) && std::isfinite(resultMinus2)) ? (resultPlus2 - resultMinus2) / (2.0 * step2) : 0.0;
            }
        }
        
        return postProcess(result, dResult_dInput1, dResult_dInput2);
    } catch (const std::exception& e) {
        // ── Penalty-based CoolProp error handling ────────────────────
        // Return a large penalty value with zero gradient rather than
        // throwing.  This keeps the residual landscape finite so that
        // trust-region and LM solvers can shrink their step / increase
        // damping smoothly rather than losing the trial evaluation
        // entirely.
        if (diagnostics_) {
            diagnostics_->push(DiagnosticSeverity::Warning, "C001",
                "CoolProp error in " + func.name + "(): " + e.what(),
                "evaluator");
        }
        constexpr double PENALTY = 1e4;
        ADValue output(PENALTY, numVariables_);
        return output;
    }
}

// ============================================================================
// BlockEvaluator Implementation
// ============================================================================

namespace {

using RoutineNames = std::set<std::string, CaseInsensitiveLess>;

bool callsAnyOf(const ExprPtr& expr, const RoutineNames& names);
bool callsAnyOf(const std::vector<StmtPtr>& body, const RoutineNames& names);

bool callsAnyOf(const ProcedureCall& call, const RoutineNames& names) {
    if (names.count(call.name)) return true;
    for (const auto& arg : call.inputArgs)
        if (callsAnyOf(arg, names)) return true;
    for (const auto& out : call.outputExprs)
        if (callsAnyOf(out, names)) return true;
    return false;
}

bool callsAnyOf(const ExprPtr& expr, const RoutineNames& names) {
    if (!expr) return false;
    if (expr->is<FunctionCall>()) {
        const auto& call = expr->as<FunctionCall>();
        if (names.count(call.name)) return true;
        for (const auto& arg : call.args)
            if (callsAnyOf(arg, names)) return true;
        for (const auto& named : call.namedArgs)
            if (callsAnyOf(named.second, names)) return true;
        return false;
    }
    if (expr->is<UnaryOp>()) return callsAnyOf(expr->as<UnaryOp>().operand, names);
    if (expr->is<BinaryOp>())
        return callsAnyOf(expr->as<BinaryOp>().left, names) || callsAnyOf(expr->as<BinaryOp>().right, names);
    if (expr->is<Variable>()) {
        for (const auto& index : expr->as<Variable>().indices)
            if (callsAnyOf(index, names)) return true;
    }
    return false;
}

bool callsAnyOf(const std::vector<StmtPtr>& body, const RoutineNames& names) {
    for (const auto& stmt : body) {
        if (!stmt) continue;
        if (stmt->is<Equation>()) {
            const auto& eq = stmt->as<Equation>();
            if (callsAnyOf(eq.lhs, names) || callsAnyOf(eq.rhs, names)) return true;
        } else if (stmt->is<Assignment>()) {
            const auto& assign = stmt->as<Assignment>();
            if (callsAnyOf(assign.lhs, names) || callsAnyOf(assign.rhs, names)) return true;
        } else if (stmt->is<IfThenElse>()) {
            const auto& branch = stmt->as<IfThenElse>();
            if (callsAnyOf(branch.condition, names) || callsAnyOf(branch.thenBranch, names) ||
                callsAnyOf(branch.elseBranch, names)) return true;
        } else if (stmt->is<ProcedureCall>()) {
            if (callsAnyOf(stmt->as<ProcedureCall>(), names)) return true;
        } else if (stmt->is<Duplicate>()) {
            const auto& loop = stmt->as<Duplicate>();
            if (callsAnyOf(loop.start, names) || callsAnyOf(loop.end, names) || callsAnyOf(loop.body, names)) return true;
        } else if (stmt->is<RepeatUntil>()) {
            const auto& loop = stmt->as<RepeatUntil>();
            if (callsAnyOf(loop.body, names) || callsAnyOf(loop.condition, names)) return true;
        }
    }
    return false;
}

/// The user FUNCTIONs and PROCEDUREs whose body can execute `CALL WARNING`, directly or through
/// the routines it calls (the fixed point of "calls WARNING or a routine that does"). Empty for a
/// model without that call.
RoutineNames routinesThatWarn(const IR& ir) {
    RoutineNames names{"warning"};
    for (bool grown = true; grown;) {
        grown = false;
        for (const auto& func : ir.getFunctions())
            if (!names.count(func.name) && callsAnyOf(func.body, names)) { names.insert(func.name); grown = true; }
        for (const auto& proc : ir.getProcedures())
            if (!names.count(proc.name) && callsAnyOf(proc.body, names)) { names.insert(proc.name); grown = true; }
    }
    names.erase("warning");   // the main program cannot call it
    return names;
}

/// An equation `X$ = ...` that the block evaluator assigns (CS-BUG-STRING-ALIAS): a string defined from another
/// string variable or from a string function. A literal is set by the IR (IR::fromAST) and costs nothing here.
bool assignsString(const EquationInfo& eq) { return isStringAssignment(eq) && !eq.rhs->is<StringLiteral>(); }

}  // namespace

void BlockEvaluator::markCallWarnings(const std::set<std::string, CaseInsensitiveLess>& routines) {
    mayCallWarning_ = false;
    for (const EquationInfo* eq : equations_) {
        if (callsAnyOf(eq->lhs, routines) || callsAnyOf(eq->rhs, routines) ||
            (eq->procedureCall && callsAnyOf(*eq->procedureCall, routines))) {
            mayCallWarning_ = true;
            return;
        }
    }
}

BlockEvaluator::BlockEvaluator(const Block& block, const IR& ir, const CoolPropConfig& config)
    : variables_(block.variables), equationIds_(block.equationIds), config_(config) {
    const auto& allEquations = ir.getEquations();
    for (int eqId : equationIds_) {
        if (eqId >= 0 && eqId < static_cast<int>(allEquations.size())) {
            equations_.push_back(&allEquations[eqId]);
        }
    }
    for (size_t i = 0; i < variables_.size(); ++i) {
        varToIndex_[variables_[i]] = i;
    }
    // CALLs that compute a string (CS-BUG-STRING-CALL-OUT) and string assignments `X$ = ...`
    // (CS-BUG-STRING-ALIAS): the strings they set are read by other equations of the block, so evaluate
    // them first, each after the producers of the strings it reads.
    auto isString = [](const std::string& name) { return !name.empty() && name.back() == '$'; };
    auto stringOutputs = [&](const EquationInfo& eq) {
        std::set<std::string, CaseInsensitiveLess> names;
        if (eq.procedureCall)
            for (const auto& out : eq.procedureCall->outputVars)
                if (isString(out.name)) names.insert(out.name);
        if (assignsString(eq)) names.insert(eq.lhs->as<Variable>().name);
        return names;
    };
    std::vector<size_t> producers;
    for (size_t i = 0; i < equations_.size(); ++i)
        if (!stringOutputs(*equations_[i]).empty()) producers.push_back(i);
    if (producers.empty()) return;
    hasStringCallOutputs_ = true;
    std::vector<bool> placed(equations_.size(), false);
    while (evalOrder_.size() < producers.size()) {
        bool progress = false;
        for (size_t i : producers) {
            if (placed[i]) continue;
            // ready when no still-unplaced producer sets a string this CALL reads
            bool ready = true;
            for (size_t j : producers) {
                if (placed[j] || j == i) continue;
                for (const auto& name : stringOutputs(*equations_[j]))
                    if (equations_[i]->variables.count(name) &&
                        !stringOutputs(*equations_[i]).count(name)) ready = false;
            }
            if (!ready) continue;
            evalOrder_.push_back(i);
            placed[i] = true;
            progress = true;
        }
        if (!progress) {   // strings that depend on each other: keep the order of the block
            for (size_t i : producers)
                if (!placed[i]) { evalOrder_.push_back(i); placed[i] = true; }
        }
    }
    for (size_t i = 0; i < equations_.size(); ++i)
        if (!placed[i]) evalOrder_.push_back(i);
}

EvaluationResult BlockEvaluator::evaluate(const std::vector<double>& x,
                                           const std::map<std::string, double>& externalVars,
                                           const std::map<std::string, std::string>& externalStringVars,
                                           bool computeJacobian,
                                           std::vector<Diagnostic>* callWarnings,
                                           std::vector<std::string>* clampReport) {
    auto eval_start = std::chrono::high_resolution_clock::now();
    g_profilingStats.block_evaluations++;
    
    if (x.size() != variables_.size()) {
        throw std::runtime_error("State vector size mismatch");
    }
    
    EvaluationResult result;
    result.residuals.resize(equations_.size());
    if (computeJacobian) {
        result.jacobian.resize(equations_.size(), std::vector<double>(variables_.size(), 0.0));
    }
    
    ExpressionEvaluator exprEval(variables_.size(), config_);
    exprEval.setResidualOnly(!computeJacobian);
    if (diagnostics_) exprEval.setDiagnostics(diagnostics_);
    if (callWarnings) exprEval.setCallWarningSink(callWarnings);
    InputClampLog clampLog;
    clampLog.report = clampReport;
    exprEval.setInputClampLog(&clampLog);
    if (lookupTableStore_) exprEval.setLookupTableStore(lookupTableStore_);
    for (const auto& [name, func] : functions_) exprEval.registerFunction(func);
    for (const auto& [name, proc] : procedures_) exprEval.registerProcedure(proc);
    
    for (size_t i = 0; i < variables_.size(); ++i) {
        exprEval.setVariable(variables_[i], ADValue::independent(x[i], i, variables_.size()));
    }
    for (const auto& [name, value] : externalVars) {
        if (varToIndex_.find(name) == varToIndex_.end()) {
            exprEval.setVariable(name, ADValue::constant(value, variables_.size()));
        }
    }
    for (const auto& [name, value] : externalStringVars) {
        exprEval.setStringVariable(name, value);
    }
    
    // Case-insensitive string comparison
    auto ciEqual = [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (size_t k = 0; k < a.size(); ++k)
            if (std::tolower(static_cast<unsigned char>(a[k])) !=
                std::tolower(static_cast<unsigned char>(b[k]))) return false;
        return true;
    };

    // Track CALLs already processed (to skip sibling equations)
    std::set<int> processedCallIds;
    
    for (size_t k = 0; k < equations_.size(); ++k) {
        const size_t eq = evalOrder_.empty() ? k : evalOrder_[k];
        const EquationInfo* eqInfo = equations_[eq];
        
        if (eqInfo->procedureCall) {
            // Skip sibling CALL equations already handled by their primary
            if (processedCallIds.count(eqInfo->callId)) {
                continue;
            }
            processedCallIds.insert(eqInfo->callId);

            const auto& call = *eqInfo->procedureCall;
            const size_t nOutputs = call.outputVars.size();

            // Find ALL sibling CALL equations in this block (same callId),
            // sorted by global equation ID.  The IR creates sibling equations
            // with consecutive IDs in output order, so sorting by GID gives the
            // correct output-index-to-equation mapping regardless of the SCC
            // ordering within the block.
            std::vector<std::pair<int, size_t>> siblingGidAndLocal; // (globalId, localIdx)
            for (size_t j = 0; j < equations_.size(); ++j) {
                if (equations_[j]->callId == eqInfo->callId) {
                    siblingGidAndLocal.push_back({equationIds_[j], j});
                }
            }
            std::sort(siblingGidAndLocal.begin(), siblingGidAndLocal.end());

            // Store the initial state of the outputs to compute residuals.
            // These "old" values come from x[] (via the AD-seeded variables
            // set up above), so residuals carry correct AD gradients.
            std::map<std::string, ADValue> oldOutputs;
            for (const auto& var : call.outputVars) {
                std::string name = exprEval.resolveVariableName(var);
                oldOutputs[name] = exprEval.getVariable(name);
            }
            
            // Evaluate procedure call ONCE — this updates exprEval's state
            // for all output variables with the procedure-computed values.
            exprEval.evaluateProcedureCall(call);

            // A string output is no unknown of the numeric system: hand it to the caller
            if (hasStringCallOutputs_) {
                for (const auto& var : call.outputVars) {
                    const std::string name = exprEval.resolveVariableName(var);
                    if (!name.empty() && name.back() == '$' && exprEval.hasStringVariable(name))
                        result.stringOutputs[name] = exprEval.getStringVariable(name);
                }
            }
            
            // Write one residual per output to its corresponding sibling
            // CALL equation's slot.  siblingGidAndLocal[i] maps to the i-th
            // output because the IR creates equations in output order.
            for (size_t i = 0; i < nOutputs; ++i) {
                std::string outName = exprEval.resolveVariableName(call.outputVars[i]);
                ADValue newValue = exprEval.getVariable(outName);
                ADValue residual = oldOutputs[outName] - newValue;

                if (i < siblingGidAndLocal.size()) {
                    size_t targetIdx = siblingGidAndLocal[i].second;
                    result.residuals[targetIdx] = residual.value;
                    if (computeJacobian) {
                        for (size_t var = 0; var < variables_.size(); ++var) {
                            result.jacobian[targetIdx][var] = residual.gradient[var];
                        }
                    }
                }
            }

            // Restore evaluator state for outputs that are NOT matched to
            // one of this CALL's sibling equations.  If a CALL output is
            // matched to a different equation (e.g. constraint "DELTAT_pp =
            // DT_min"), leaving the CALL-computed value in the evaluator
            // would corrupt that equation's residual and Jacobian, because
            // the AD gradient w.r.t. x[varIdx] would be lost.
            for (size_t i = 0; i < nOutputs; ++i) {
                std::string outName = exprEval.resolveVariableName(call.outputVars[i]);

                // Find this output in the block's variable list
                bool isBlockVar = false;
                size_t blockVarIdx = 0;
                for (size_t j = 0; j < variables_.size(); ++j) {
                    if (ciEqual(variables_[j], outName)) {
                        isBlockVar = true;
                        blockVarIdx = j;
                        break;
                    }
                }

                bool restore = false;
                if (!isBlockVar) {
                    // External variable (solved by another block) — restore
                    restore = true;
                } else {
                    // Check whether the equation matched to this variable is
                    // a sibling CALL equation (same CALL).  If not,
                    // the output is "foreign" and must be restored.
                    if (equations_[blockVarIdx]->callId != eqInfo->callId) {
                        restore = true;
                    }
                }

                if (restore) {
                    exprEval.setVariable(outName, oldOutputs[outName]);
                }
            }

            continue;
        }

        if (!eqInfo->lhs || !eqInfo->rhs) continue;

        // `X$ = <string expression>`: X$ gets the string, handed to the caller like a string output of a
        // CALL. The number that stands for X$ in the system stays 0 (as for a string literal).
        if (hasStringCallOutputs_ && assignsString(*eqInfo)) {
            const std::string& name = eqInfo->lhs->as<Variable>().name;
            const std::string value = exprEval.evaluateString(eqInfo->rhs);
            exprEval.setStringVariable(name, value);
            result.stringOutputs[name] = value;
            ADValue standIn = exprEval.evaluate(eqInfo->lhs);
            result.residuals[eq] = standIn.value;
            if (computeJacobian) {
                for (size_t var = 0; var < variables_.size(); ++var) {
                    result.jacobian[eq][var] = standIn.gradient[var];
                }
            }
            continue;
        }

        ADValue lhs = exprEval.evaluate(eqInfo->lhs);
        ADValue rhs = exprEval.evaluate(eqInfo->rhs);
        ADValue residual = lhs - rhs;
        result.residuals[eq] = residual.value;
        if (computeJacobian) {
            for (size_t var = 0; var < variables_.size(); ++var) {
                result.jacobian[eq][var] = residual.gradient[var];
            }
        }
    }
    
    // Collect warnings from missing variables
    for (const auto& var : exprEval.getMissingVariables()) {
        result.warnings.push_back("Undefined variable: " + var + " (using default value 1.0)");
    }
    if (clampLog.fired) inputClamped_->store(true, std::memory_order_relaxed);
    
    auto eval_end = std::chrono::high_resolution_clock::now();
    g_profilingStats.block_eval_time_ms += std::chrono::duration<double, std::milli>(eval_end - eval_start).count();
    
    return result;
}

// ============================================================================
// SystemEvaluator Implementation
// ============================================================================

SystemEvaluator::SystemEvaluator(const IR& ir, const StructuralAnalysisResult& analysisResult, const CoolPropConfig& config)
    : ir_(ir), config_(config) {
    for (const auto& block : analysisResult.blocks) {
        blockEvaluators_.emplace_back(block, ir_, config_);
    }
    
    // Register functions and procedures in all block evaluators
    for (auto& blockEval : blockEvaluators_) {
        for (const auto& func : ir.getFunctions()) {
            blockEval.registerFunction(func);
        }
        for (const auto& proc : ir.getProcedures()) {
            blockEval.registerProcedure(proc);
        }
    }

    // The blocks that can reach a CALL WARNING (none when the model has no such call)
    const RoutineNames warningRoutines = routinesThatWarn(ir);
    if (!warningRoutines.empty())
        for (auto& blockEval : blockEvaluators_) blockEval.markCallWarnings(warningRoutines);
}

void SystemEvaluator::setVariableValue(const std::string& name, double value) {
    currentState_[name] = value;
}

void SystemEvaluator::setStringVariableValue(const std::string& name, const std::string& value) {
    currentStringState_[name] = value;
}

double SystemEvaluator::getVariableValue(const std::string& name) const {
    auto it = currentState_.find(name);
    if (it != currentState_.end()) return it->second;
    
    const auto* varInfo = ir_.getVariable(name);
    if (varInfo && varInfo->solutionValue) return *varInfo->solutionValue;
    if (varInfo && varInfo->guessValue) return *varInfo->guessValue;
    return 1.0;
}

std::string SystemEvaluator::getStringVariableValue(const std::string& name) const {
    auto it = currentStringState_.find(name);
    if (it != currentStringState_.end()) return it->second;
    
    const auto* varInfo = ir_.getVariable(name);
    if (varInfo && varInfo->solutionStringValue) return *varInfo->solutionStringValue;
    return "";
}

void SystemEvaluator::initializeFromGuesses() {
    for (const auto& [name, info] : ir_.getVariables()) {
        if (info.type == VariableType::String) {
            if (info.solutionStringValue) currentStringState_[name] = *info.solutionStringValue;
        } else {
            if (info.solutionValue) currentState_[name] = *info.solutionValue;
            else if (info.guessValue) currentState_[name] = *info.guessValue;
            else currentState_[name] = 1.0;
        }
    }
}

void SystemEvaluator::publishStringOutputs(size_t blockIndex) {
    if (blockIndex >= blockEvaluators_.size() || !blockEvaluators_[blockIndex].hasStringCallOutputs()) return;
    BlockEvaluator& blockEval = blockEvaluators_[blockIndex];
    std::vector<double> x;
    for (const auto& varName : blockEval.getVariables()) x.push_back(getVariableValue(varName));
    std::map<std::string, double> externalVars;
    for (const auto& [name, value] : currentState_) externalVars[name] = value;
    std::map<std::string, std::string> externalStringVars;
    for (const auto& [name, value] : currentStringState_) externalStringVars[name] = value;
    const EvaluationResult result = blockEval.evaluate(x, externalVars, externalStringVars, false);
    for (const auto& [name, value] : result.stringOutputs) currentStringState_[name] = value;
}

void SystemEvaluator::reportCallWarnings(size_t blockIndex, DiagnosticCollector& out) {
    if (blockIndex >= blockEvaluators_.size() || !blockEvaluators_[blockIndex].mayCallWarning()) return;
    std::vector<Diagnostic> raised;
    evaluateBlock(blockIndex, false, &raised);
    for (auto& warning : raised)
        if (reportedCallWarnings_.insert({warning.line, warning.message}).second) out.push(std::move(warning));
}

void SystemEvaluator::checkInputClamps(size_t blockIndex) {
    if (blockIndex >= blockEvaluators_.size() || !blockEvaluators_[blockIndex].takeInputClamped()) return;
    std::vector<std::string> clamps;
    evaluateBlock(blockIndex, false, nullptr, &clamps);
    if (clamps.empty()) return;   // the iterates were limited, the solution is not
    std::string message = clamps.front();
    for (size_t i = 1; i < clamps.size(); ++i)
        if (std::find(clamps.begin(), clamps.begin() + i, clamps[i]) == clamps.begin() + i) message += "; " + clamps[i];
    throw std::runtime_error(message);
}

EvaluationResult SystemEvaluator::evaluateBlock(size_t blockIndex, bool computeJacobian,
                                                std::vector<Diagnostic>* callWarnings,
                                                std::vector<std::string>* clampReport) {
    if (blockIndex >= blockEvaluators_.size()) throw std::out_of_range("Block index out of range");
    BlockEvaluator& blockEval = blockEvaluators_[blockIndex];
    std::vector<double> x;
    for (const auto& varName : blockEval.getVariables()) {
        x.push_back(getVariableValue(varName));
    }
    std::map<std::string, double> externalVars;
    for (const auto& [name, value] : currentState_) externalVars[name] = value;
    std::map<std::string, std::string> externalStringVars;
    for (const auto& [name, value] : currentStringState_) externalStringVars[name] = value;
    return blockEval.evaluate(x, externalVars, externalStringVars, computeJacobian, callWarnings, clampReport);
}

// ============================================================================
// Utility Functions
// ============================================================================

double compareJacobianWithFiniteDifferences(
    BlockEvaluator& evaluator,
    const std::vector<double>& x,
    const std::map<std::string, double>& externalVars,
    const std::map<std::string, std::string>& externalStringVars,
    double epsilon,
    bool verbose) {
    
    EvaluationResult adResult = evaluator.evaluate(x, externalVars, externalStringVars);
    size_t n = x.size();
    size_t m = adResult.residuals.size();
    std::vector<std::vector<double>> numJacobian(m, std::vector<double>(n, 0.0));
    
    for (size_t j = 0; j < n; ++j) {
        std::vector<double> xPlus = x;
        std::vector<double> xMinus = x;
        xPlus[j] += epsilon;
        xMinus[j] -= epsilon;
        EvaluationResult resultPlus = evaluator.evaluate(xPlus, externalVars, externalStringVars);
        EvaluationResult resultMinus = evaluator.evaluate(xMinus, externalVars, externalStringVars);
        for (size_t i = 0; i < m; ++i) {
            numJacobian[i][j] = (resultPlus.residuals[i] - resultMinus.residuals[i]) / (2.0 * epsilon);
        }
    }
    
    double maxDiff = 0.0;
    size_t maxI = 0, maxJ = 0;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double diff = std::abs(adResult.jacobian[i][j] - numJacobian[i][j]);
            if (diff > maxDiff) {
                maxDiff = diff;
                maxI = i;
                maxJ = j;
            }
        }
    }
    
    if (verbose) {
        const auto& vars = evaluator.getVariables();
        
        std::cout << "\n=== Jacobian Comparison ===\n";
        std::cout << "State vector x = [";
        for (size_t j = 0; j < n; ++j) {
            if (j > 0) std::cout << ", ";
            std::cout << vars[j] << "=" << x[j];
        }
        std::cout << "]\n";
        std::cout << "Residuals F(x) = [";
        for (size_t i = 0; i < m; ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << adResult.residuals[i];
        }
        std::cout << "]\n\n";
        
        std::cout << "Jacobian (AD) - computed using automatic differentiation:\n";
        std::cout << std::setw(12) << "";
        for (size_t j = 0; j < n; ++j) {
            std::cout << std::setw(14) << vars[j];
        }
        std::cout << "\n";
        for (size_t i = 0; i < m; ++i) {
            std::cout << std::setw(10) << "Eq" << i << " ";
            for (size_t j = 0; j < n; ++j) {
                std::cout << std::setw(14) << std::scientific << std::setprecision(4) << adResult.jacobian[i][j];
            }
            std::cout << "\n";
        }
        
        std::cout << "\nJacobian (FD) - computed using finite differences (eps=" << epsilon << "):\n";
        std::cout << std::setw(12) << "";
        for (size_t j = 0; j < n; ++j) {
            std::cout << std::setw(14) << vars[j];
        }
        std::cout << "\n";
        for (size_t i = 0; i < m; ++i) {
            std::cout << std::setw(10) << "Eq" << i << " ";
            for (size_t j = 0; j < n; ++j) {
                std::cout << std::setw(14) << std::scientific << std::setprecision(4) << numJacobian[i][j];
            }
            std::cout << "\n";
        }
        
        std::cout << "\nDifference (AD - FD):\n";
        std::cout << std::setw(12) << "";
        for (size_t j = 0; j < n; ++j) {
            std::cout << std::setw(14) << vars[j];
        }
        std::cout << "\n";
        for (size_t i = 0; i < m; ++i) {
            std::cout << std::setw(10) << "Eq" << i << " ";
            for (size_t j = 0; j < n; ++j) {
                double diff = adResult.jacobian[i][j] - numJacobian[i][j];
                std::cout << std::setw(14) << std::scientific << std::setprecision(4) << diff;
            }
            std::cout << "\n";
        }
        
        std::cout << "\nMax absolute difference: " << maxDiff;
        if (m > 0 && n > 0) {
            std::cout << " at [Eq" << maxI << ", " << vars[maxJ] << "]";
        }
        std::cout << "\n=== End Jacobian Comparison ===\n\n";
    }
    
    return maxDiff;
}

void ProceduralEvaluator::run(ExpressionEvaluator& eval, const std::vector<StmtPtr>& body, const std::string& routine) {
    Jump jump{routine};
    runBlock(eval, body, jump);
    if (jump.pending)   // the parser refuses a GOTO whose label is not in the body
        throw std::runtime_error("GOTO " + std::to_string(jump.label) + " in " + routine + ": the body has no such statement label");
}

void ProceduralEvaluator::runBlock(ExpressionEvaluator& eval, const std::vector<StmtPtr>& block, Jump& jump) {
    for (size_t i = 0; i < block.size(); ++i) {
        execute(eval, block[i], jump);
        if (!jump.pending) continue;
        const auto target = std::find_if(block.begin(), block.end(), [&](const StmtPtr& s) {
            return s->is<Label>() && s->as<Label>().number == jump.label;
        });
        if (target == block.end()) return;   // the label is in an enclosing block
        jump.pending = false;
        i = static_cast<size_t>(target - block.begin());   // the loop goes on after the label
    }
}

void ProceduralEvaluator::execute(ExpressionEvaluator& eval, const StmtPtr& stmt, Jump& jump) {
    if (stmt->is<Assignment>()) {
        const auto& assign = stmt->as<Assignment>();
        if (assign.lhs->is<Variable>()) {
            const auto& var = assign.lhs->as<Variable>();
            std::string name = eval.resolveVariableName(var);
            if (!name.empty() && name.back() == '$') {
                eval.setStringVariable(name, eval.evaluateString(assign.rhs));
            } else {
                eval.setVariable(name, eval.evaluate(assign.rhs));
            }
        }
    } else if (stmt->is<Equation>()) {
        // In EES procedural blocks, '=' is also treated as assignment. A left-hand side that is not
        // a variable is refused by the parser's check of the bodies (CS-BUG-IMPLICIT-PROC); a
        // lookup-table cell write `LOOKUP(...) = ...` is reported by the parser and not in the body
        // (CS-BUG-LOOKUP-WRITE).
        const auto& eq = stmt->as<Equation>();
        if (eq.lhs->is<Variable>()) {
            const auto& var = eq.lhs->as<Variable>();
            std::string name = eval.resolveVariableName(var);
            if (!name.empty() && name.back() == '$') {
                eval.setStringVariable(name, eval.evaluateString(eq.rhs));
            } else {
                eval.setVariable(name, eval.evaluate(eq.rhs));
            }
        }
    } else if (stmt->is<IfThenElse>()) {
        const auto& ite = stmt->as<IfThenElse>();
        ADValue cond = eval.evaluate(ite.condition);
        runBlock(eval, cond.value > 0.5 ? ite.thenBranch : ite.elseBranch, jump);
    } else if (stmt->is<ProcedureCall>()) {
        const auto& call = stmt->as<ProcedureCall>();
        eval.evaluateProcedureCall(call);
    } else if (stmt->is<Duplicate>()) {
        const auto& dup = stmt->as<Duplicate>();
        int startVal = static_cast<int>(std::round(eval.evaluate(dup.start).value));
        int endVal = static_cast<int>(std::round(eval.evaluate(dup.end).value));
        for (int i = startVal; i <= endVal && !jump.pending; ++i) {
            eval.setVariable(dup.iteratorVar, ADValue::constant(static_cast<double>(i), eval.getNumVariables()));
            runBlock(eval, dup.body, jump);
        }
    } else if (stmt->is<RepeatUntil>()) {
        const auto& ru = stmt->as<RepeatUntil>();
        for (long iter = 1;; ++iter) {
            runBlock(eval, ru.body, jump);
            if (jump.pending) break;   // a GOTO out of the loop
            ADValue cond = eval.evaluate(ru.condition);
            if (cond.value > 0.5) break;
            if (iter >= kMaxRepeatIterations)
                throw std::runtime_error("REPEAT loop in " + jump.routine + " (line " + std::to_string(stmt->sourceLineNumber) +
                                         ") did not end after " + std::to_string(kMaxRepeatIterations) +
                                         " iterations: an UNTIL condition that never becomes true?");
        }
    } else if (stmt->is<Goto>()) {
        if (++jump.executed > kMaxGotos)
            throw std::runtime_error("more than " + std::to_string(kMaxGotos) + " GOTO statements executed in " + jump.routine +
                                     " (the last one at line " + std::to_string(stmt->sourceLineNumber) +
                                     "): a loop built with GOTO that never ends?");
        jump.pending = true;
        jump.label = stmt->as<Goto>().label;
    }
}

std::string generateEvaluatorReport(const SystemEvaluator& evaluator) {
    std::ostringstream oss;
    oss << "# Evaluator Report\n\n## Blocks Summary\n\n| Block | Size | Type | Variables | Equations |\n|-------|------|------|-----------|-----------|\n";
    for (size_t i = 0; i < evaluator.getNumBlocks(); ++i) {
        const BlockEvaluator& block = evaluator.getBlock(i);
        oss << "| " << i << " | " << block.size() << " | " << (block.isExplicit() ? "Explicit" : "Algebraic loop") << " | ";
        
        const auto& vars = block.getVariables();
        for (size_t j = 0; j < vars.size() && j < 3; ++j) {
            if (j > 0) oss << ", ";
            oss << "`" << vars[j] << "`";
        }
        if (vars.size() > 3) oss << ", ...";
        oss << " | ";

        const auto& eqIds = block.getEquationIds();
        for (size_t j = 0; j < eqIds.size() && j < 3; ++j) {
            if (j > 0) oss << ", ";
            oss << eqIds[j];
        }
        if (eqIds.size() > 3) oss << ", ...";
        oss << " |\n";
    }
    
    // Test evaluability of each block
    oss << "\n## Block Evaluation Test\n\n";
    oss << "Testing residual evaluation at current state:\n\n";
    oss << "| Block | Status | Max Residual |\n|-------|--------|--------------|\n";
    
    for (size_t i = 0; i < evaluator.getNumBlocks(); ++i) {
        oss << "| " << i << " | ";
        try {
            auto result = const_cast<SystemEvaluator&>(evaluator).evaluateBlock(i);
            double maxRes = 0.0;
            for (double r : result.residuals) {
                maxRes = std::max(maxRes, std::abs(r));
            }
            oss << "OK | " << std::scientific << std::setprecision(4) << maxRes << " |\n";
        } catch (const std::exception& e) {
            oss << "Error: " << e.what() << " | - |\n";
        }
    }

    oss << "\n## Current State\n\n";
    const auto& state = evaluator.getAllVariables();
    if (!state.empty()) {
        oss << "| Variable | Value |\n|----------|-------|\n";
        for (const auto& [name, value] : state) oss << "| " << name << " | " << value << " |\n";
    } else oss << "No variables set.\n";
    return oss.str();
}

// ============================================================================
// Lookup table / interpolation function evaluation
// ============================================================================

std::string ExpressionEvaluator::lookupWhere() const {
    // Where an error comes from when the call is inside a FUNCTION / PROCEDURE body.
    return routine_.empty() ? std::string() : " (in " + routine_ + ")";
}

const LookupTable& ExpressionEvaluator::lookupTableFor(const FunctionCall& func, const std::string& tableName) {
    // Resolve a table: emits a diagnostic and throws if not found
    if (!lookupTableStore_) {
        std::string msg = func.name + "(): lookup table '" + tableName + "' not found" +
            (routine_.empty() ? " (" : " (in " + routine_ + "; ") +
            "no table store is available in this context)";
        if (diagnostics_) diagnostics_->push(DiagnosticSeverity::Error, "L002", msg, "lookup_table");
        throw std::runtime_error(msg);
    }
    const LookupTable* tbl = lookupTableStore_->get(tableName);
    if (!tbl) {
        std::string msg = func.name + "(): lookup table '" + tableName + "' not found" + lookupWhere() +
            ". To fix: place '<model>-" + tableName +
            ".csv' next to the model file, or load it via the Lookup Tables panel.";
        if (diagnostics_) diagnostics_->push(DiagnosticSeverity::Error, "L002", msg, "lookup_table");
        throw std::runtime_error(msg);
    }
    return *tbl;
}

const LookupTable& ExpressionEvaluator::defaultLookupTable(const FunctionCall& func) {
    static const std::string defaultName = "lookup_1";
    if (lookupTableStore_ && !lookupTableStore_->has(defaultName)) {
        const std::vector<std::string> names = lookupTableStore_->names();
        if (names.size() == 1) return lookupTableFor(func, names.front());   // unambiguous: EES reads its first table
        if (names.size() > 1) {
            std::string found;
            for (const auto& n : names) found += (found.empty() ? "'" : ", '") + n + "'";
            const std::string msg = func.name + "(): the nameless form " + func.name + "(row, col) reads the default "
                "lookup table '" + defaultName + "', which this model does not have" + lookupWhere() +
                ". To fix: place '<model>-" + defaultName + ".csv' next to the model file, or name the table: " +
                func.name + "('table', row, col). Tables found: " + found + ".";
            if (diagnostics_) diagnostics_->push(DiagnosticSeverity::Error, "L002", msg, "lookup_table");
            throw std::runtime_error(msg);
        }
    }
    return lookupTableFor(func, defaultName);   // not found: the usual error, which names the file to provide
}

size_t ExpressionEvaluator::lookupColumnOf(const FunctionCall& func, const LookupTable& tbl, size_t argIdx) {
    // Resolve a column by name or 1-based number (as double arg)
    // EES allows both "T_K" and a numeric column index.
    if (argIdx >= func.args.size()) {
        throw std::runtime_error(func.name + "(): missing column argument" + lookupWhere());
    }
    const auto& arg = func.args[argIdx];
    if (isStringExpression(arg)) {
        std::string colName = evaluateString(arg);
        size_t col = tbl.columnIndex(colName);
        if (col == 0) {
            throw std::runtime_error(func.name + "(): column '" + colName +
                                     "' not found in table '" + tbl.name() + "'" + lookupWhere());
        }
        return col;
    } else {
        // Numeric column index (0, which no column has, when it is not a positive number)
        const double col = std::round(evaluate(arg).value);
        return col >= 1 && col < 1e18 ? static_cast<size_t>(col) : 0;
    }
}

size_t ExpressionEvaluator::lookupCellColumnOf(const FunctionCall& func, const LookupTable& tbl, size_t argIdx) {
    const size_t col = lookupColumnOf(func, tbl, argIdx);
    if (col < 1 || col > tbl.numCols()) {
        throw std::runtime_error(func.name + "(): column " + std::to_string(col) + " is outside the table '" +
                                 tbl.name() + "' (1 to " + std::to_string(tbl.numCols()) + " columns)" + lookupWhere());
    }
    return col;
}

std::string ExpressionEvaluator::evaluateLookupString(const FunctionCall& func) {
    // LOOKUP$('table', row, col): the text of a cell, the way EES reads strings from a lookup table
    // (the replacement of string arrays: the strings are in a table, the model keeps the row number).
    if (func.args.size() != 3 || !func.namedArgs.empty()) {
        throw std::runtime_error(func.name + "(): expects 3 arguments: (table, row, col)" + lookupWhere());
    }
    const LookupTable& tbl = lookupTableFor(func, evaluateString(func.args[0]));

    const double rowValue = evaluate(func.args[1]).value;
    const double row = std::round(rowValue);
    if (!(row >= 1 && row <= static_cast<double>(tbl.numRows()))) {
        std::ostringstream msg;
        msg << func.name << "(): row " << rowValue << " is outside the table '" << tbl.name() << "' (1 to "
            << tbl.numRows() << " rows)" << lookupWhere();
        throw std::runtime_error(msg.str());
    }
    const size_t col = lookupCellColumnOf(func, tbl, 2);
    return tbl.text(static_cast<size_t>(row), col);
}

ExpressionEvaluator::EesInterpolation ExpressionEvaluator::resolveEesInterpolation(const FunctionCall& func, size_t nCols) {
    if (func.args.size() != nCols + 1 || func.namedArgs.size() != nCols - 1) {
        const bool oneD = nCols == 2;
        throw std::runtime_error(func.name + "(): the EES form takes the table, " +
            (oneD ? "2 columns and 1 named argument" : "3 columns and 2 named arguments") + " <column>=<value> for the input " +
            (oneD ? "column: " : "columns: ") + func.name +
            (oneD ? "('table', 'Col1', 'Col2', Col2=value)" : "('table', 'Col1', 'Col2', 'Col3', Col1=x, Col2=y)") +
            "; found " + std::to_string(func.args.size()) + " positional and " + std::to_string(func.namedArgs.size()) +
            " named argument(s)" + lookupWhere());
    }

    EesInterpolation call;
    call.table = &lookupTableFor(func, evaluateString(func.args[0]));
    std::vector<size_t> cols;
    for (size_t i = 1; i <= nCols; ++i) cols.push_back(lookupColumnOf(func, *call.table, i));

    // Each named argument stands for the first column of the call, not yet named, that has its name
    std::vector<const ExprPtr*> valueOf(nCols, nullptr);
    for (const auto& [argName, valueExpr] : func.namedArgs) {
        const size_t col = call.table->columnIndex(argName);
        size_t pos = nCols;
        bool listed = false;
        for (size_t p = 0; p < nCols && col != 0; ++p) {
            if (cols[p] != col) continue;
            listed = true;
            if (!valueOf[p]) { pos = p; break; }
        }
        if (pos == nCols) {
            std::string columns;
            for (size_t c : cols)
                columns += (columns.empty() ? "'" : ", '") +
                           (c >= 1 && c <= call.table->numCols() ? call.table->columnNames()[c - 1] : std::to_string(c)) + "'";
            throw std::runtime_error(func.name + "(): the named argument '" + argName + "' " +
                (listed ? "names a column of the call a second time" : "is not one of the columns of the call") +
                " (" + columns + ") of table '" + call.table->name() + "'" + lookupWhere());
        }
        valueOf[pos] = &valueExpr;
    }
    for (size_t p = 0; p < nCols; ++p) {
        if (valueOf[p]) {
            call.inputCols.push_back(cols[p]);
            call.inputValues.push_back(evaluate(*valueOf[p]));
        } else {
            call.outputCol = cols[p];
        }
    }
    return call;
}

ADValue ExpressionEvaluator::evaluateLookupFunction(const FunctionCall& func) {
    std::string name = func.name;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);

    // Helper: extract a string literal argument (table name or column name)
    auto getString = [&](size_t idx) -> std::string {
        if (idx >= func.args.size()) {
            throw std::runtime_error(func.name + "(): missing string argument at position " +
                                     std::to_string(idx) + lookupWhere());
        }
        return evaluateString(func.args[idx]);
    };

    // Helper: evaluate a numeric argument
    auto getNum = [&](size_t idx) -> ADValue {
        if (idx >= func.args.size()) {
            throw std::runtime_error(func.name + "(): missing numeric argument at position " +
                                     std::to_string(idx) + lookupWhere());
        }
        return evaluate(func.args[idx]);
    };

    // ------------------------------------------------------------------
    // NLOOKUPROWS('table')
    // ------------------------------------------------------------------
    if (name == "nlookuprows") {
        std::string tableName = getString(0);
        const LookupTable& tbl = lookupTableFor(func, tableName);
        return ADValue::constant(static_cast<double>(tbl.numRows()), numVariables_);
    }

    // ------------------------------------------------------------------
    // NLOOKUPCOLUMNS('table')
    // ------------------------------------------------------------------
    if (name == "nlookupcolumns") {
        std::string tableName = getString(0);
        const LookupTable& tbl = lookupTableFor(func, tableName);
        return ADValue::constant(static_cast<double>(tbl.numCols()), numVariables_);
    }

    // ------------------------------------------------------------------
    // LOOKUP$ROW('table', col, key)
    // The first row whose column `col` (a name or a number) holds the key: a string, matched without regard
    // to case against the text of the cells (what LOOKUP$ reads), or a number. A key that no row holds is an
    // error, never a row 0. The row is an integer: its derivative is zero.
    // ------------------------------------------------------------------
    if (name == "lookup$row") {
        if (func.args.size() != 3 || !func.namedArgs.empty()) {
            throw std::runtime_error(func.name + "(): expects 3 arguments: (table, col, key)" + lookupWhere());
        }
        const LookupTable& tbl = lookupTableFor(func, getString(0));
        const size_t col = lookupCellColumnOf(func, tbl, 1);
        const ExprPtr& keyArg = func.args[2];
        size_t row;
        std::string shownKey;   // the key as the error shows it
        if (isStringExpression(keyArg)) {
            const std::string key = evaluateString(keyArg);
            row = tbl.findRowByText(col, key);
            shownKey = "'" + key + "'";
        } else {
            const double key = evaluate(keyArg).value;
            row = tbl.findRowByValue(col, key);
            shownKey = eesFormatNumber(key);
        }
        if (row == 0) {
            throw std::runtime_error(func.name + "(): no row of table '" + tbl.name() + "' holds " + shownKey +
                                     " in column '" + tbl.columnNames()[col - 1] + "'" + lookupWhere());
        }
        return ADValue::constant(static_cast<double>(row), numVariables_);
    }

    // ------------------------------------------------------------------
    // LOOKUPCOL('table', 'colname')
    // LOOKUPCOL1('table', 'colname')  — identical to LOOKUPCOL in CoolSolve
    // ------------------------------------------------------------------
    if (name == "lookupcol" || name == "lookupcol1") {
        std::string tableName = getString(0);
        std::string colName   = getString(1);
        const LookupTable& tbl = lookupTableFor(func, tableName);
        size_t col = tbl.columnIndex(colName);
        return ADValue::constant(static_cast<double>(col), numVariables_);
    }

    // ------------------------------------------------------------------
    // LOOKUPCELLEMPTY('table', row, col)
    // Returns 1 if the cell is NaN (empty), 0 otherwise.
    // ------------------------------------------------------------------
    if (name == "lookupcellempty") {
        std::string tableName = getString(0);
        const LookupTable& tbl = lookupTableFor(func, tableName);
        ADValue rowV = getNum(1);
        size_t col   = lookupColumnOf(func, tbl, 2);
        size_t row   = static_cast<size_t>(std::round(rowV.value));
        double v     = tbl.value(row, col);
        return ADValue::constant(std::isnan(v) ? 1.0 : 0.0, numVariables_);
    }

    // ------------------------------------------------------------------
    // LOOKUP('table', row, col)
    // LOOKUP(row, col)  — the nameless form: the default table (see defaultLookupTable)
    // TABLEVALUE('table', row, col)
    // TABLEVALUE#('table', row, col)  — same as LOOKUP in CoolSolve
    //
    // EES-compatible behaviour:
    //  - Non-integer row/col: linear (bilinear) interpolation.
    //  - row < 1 or row > nRows: clamp to first/last row (derivative = 0).
    //  - col < 1 or col > nCols: clamp to first/last col (derivative = 0).
    //  - col can be a string (column name → exact integer) or a numeric
    //    value that may itself be non-integer.
    // ------------------------------------------------------------------
    if (name == "lookup" || name == "tablevalue" || name == "tablevalue#" || name == "tablerun#") {
        // Two arguments, the first not a string: LOOKUP(row, col), the table argument is left out
        const bool nameless = name == "lookup" && func.args.size() == 2 && !isStringExpression(func.args[0]);
        if (!nameless && func.args.size() < 3) {
            throw std::runtime_error(func.name + "(): expects 3 arguments: (table, row, col)" +
                                     (name == "lookup" ? " or, for the default table, 2: (row, col)" : "") + lookupWhere());
        }
        const size_t rowIdx = nameless ? 0 : 1;   // the row argument; the column argument follows
        const LookupTable& tbl = nameless ? defaultLookupTable(func) : lookupTableFor(func, getString(0));
        ADValue rowV = getNum(rowIdx);

        // Column: string name → exact integer; numeric → possibly non-integer
        ADValue colV;
        const auto& colArg = func.args[rowIdx + 1];
        if (isStringExpression(colArg)) {
            std::string colName = evaluateString(colArg);
            size_t col = tbl.columnIndex(colName);
            if (col == 0) {
                throw std::runtime_error(func.name + "(): column '" + colName +
                                         "' not found in table '" + tbl.name() + "'" + lookupWhere());
            }
            colV = ADValue::constant(static_cast<double>(col), numVariables_);
        } else {
            colV = evaluate(colArg);
        }

        return tbl.lookup(rowV, colV, diagnostics_);
    }

    // ------------------------------------------------------------------
    // SUMLOOKUP / AVGLOOKUP / MAXLOOKUP / MINLOOKUP / STDDEVLOOKUP
    //   ('table', col)
    // ------------------------------------------------------------------
    if (name == "sumlookup" || name == "avglookup" || name == "maxlookup" ||
        name == "minlookup" || name == "stddevlookup") {
        std::string tableName = getString(0);
        const LookupTable& tbl = lookupTableFor(func, tableName);
        size_t col = lookupColumnOf(func, tbl, 1);

        double result;
        if      (name == "sumlookup")    result = tbl.sumCol(col);
        else if (name == "avglookup")    result = tbl.avgCol(col);
        else if (name == "maxlookup")    result = tbl.maxCol(col);
        else if (name == "minlookup")    result = tbl.minCol(col);
        else /* stddevlookup */          result = tbl.stddevCol(col);

        return ADValue::constant(result, numVariables_);
    }

    // ------------------------------------------------------------------
    // INTERPOLATE('table', 'xcol', 'ycol', xval)
    // INTERPOLATE1 — same as INTERPOLATE
    // EES form: INTERPOLATE('table', 'Col1', 'Col2', Col2=value), Col1 where Col2 = value
    // ------------------------------------------------------------------
    if (name == "interpolate" || name == "interpolate1") {
        if (!func.namedArgs.empty()) {
            const EesInterpolation call = resolveEesInterpolation(func, 2);
            return call.table->interpolate1D(call.inputCols[0], call.outputCol, call.inputValues[0], diagnostics_);
        }
        if (func.args.size() < 4) {
            throw std::runtime_error(func.name + "(): expects 4 arguments: (table, xcol, ycol, xval), or the "
                                     "EES form (table, col1, col2, col2=value)" + lookupWhere());
        }
        std::string tableName = getString(0);
        const LookupTable& tbl = lookupTableFor(func, tableName);
        size_t xcol = lookupColumnOf(func, tbl, 1);
        size_t ycol = lookupColumnOf(func, tbl, 2);
        ADValue xval = getNum(3);
        return tbl.interpolate1D(xcol, ycol, xval, diagnostics_);
    }

    // ------------------------------------------------------------------
    // INTERPOLATE2('table', 'xcol', 'ycol', 'zcol', xval, yval)
    // INTERPOLATE2DM — same as INTERPOLATE2 in CoolSolve
    // EES form: INTERPOLATE2('table', 'Col1', 'Col2', 'Col3', Col1=x, Col2=y), Col3 where Col1 = x and Col2 = y
    // ------------------------------------------------------------------
    if (name == "interpolate2" || name == "interpolate2dm") {
        if (!func.namedArgs.empty()) {
            const EesInterpolation call = resolveEesInterpolation(func, 3);
            return call.table->interpolate2D(call.inputCols[0], call.inputCols[1], call.outputCol,
                                             call.inputValues[0], call.inputValues[1], diagnostics_);
        }
        if (func.args.size() < 6) {
            throw std::runtime_error(func.name + "(): expects 6 arguments: (table, xcol, ycol, zcol, xval, yval), "
                                     "or the EES form (table, col1, col2, col3, col1=x, col2=y)" + lookupWhere());
        }
        std::string tableName = getString(0);
        const LookupTable& tbl = lookupTableFor(func, tableName);
        size_t xcol = lookupColumnOf(func, tbl, 1);
        size_t ycol = lookupColumnOf(func, tbl, 2);
        size_t zcol = lookupColumnOf(func, tbl, 3);
        ADValue xval = getNum(4);
        ADValue yval = getNum(5);
        return tbl.interpolate2D(xcol, ycol, zcol, xval, yval, diagnostics_);
    }

    // Should not reach here given the dispatch set check above
    throw std::runtime_error("evaluateLookupFunction: unhandled function '" + func.name + "'");
}

}  // namespace coolsolve
