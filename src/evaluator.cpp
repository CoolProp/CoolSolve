#include "coolsolve/evaluator.h"
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

/// Per-thread cache of AbstractState objects keyed by "backend::fluidName".
/// Avoids repeated string parsing, fluid lookup, and EOS initialisation that
/// PropsSI performs on every call.
struct AbstractStateCache {
    std::map<std::string, std::shared_ptr<CoolProp::AbstractState>> states;
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
        auto state = std::shared_ptr<CoolProp::AbstractState>(
            CoolProp::AbstractState::factory(backend, fluidName));
        states[key] = state;
        return state;
    }

    void clear() { states.clear(); hits = misses = 0; }
};

static thread_local AbstractStateCache g_abstractStateCache;

/// Map two CoolProp parameter enums to the correct input_pairs enum.
static CoolProp::input_pairs getInputPair(CoolProp::parameters p1,
                                           CoolProp::parameters p2) {
    // Use CoolProp's generate_update_pair to determine the canonical ordering
    double out1, out2;
    return CoolProp::generate_update_pair(p1, 0.0, p2, 0.0, out1, out2);
}

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

// Humid Air parameter mapping
enum class HAParam { T, P, W, R, D, B, H, S, V, M, K, C, Unknown };

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
    if (n == "cp" || n == "cv" || n == "specheat") return HAParam::C;
    return HAParam::Unknown;
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
        default: return "";
    }
}

static UnitType haParamToUnitType(HAParam p) {
    switch (p) {
        case HAParam::T:
        case HAParam::D:
        case HAParam::B: return UnitType::Temperature;
        case HAParam::P: return UnitType::Pressure;
        case HAParam::H: return UnitType::SpecificEnergy;
        case HAParam::S: return UnitType::SpecificEntropy;
        case HAParam::V: return UnitType::Density; // m3/kg is inverse density
        case HAParam::M: return UnitType::Viscosity;
        case HAParam::K: return UnitType::Conductivity;
        case HAParam::C: return UnitType::SpecificHeat;
        default: return UnitType::Dimensionless;
    }
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

// ============================================================================
// ExpressionEvaluator Implementation
// ============================================================================

// Internal procedural statement evaluator (forward declaration)
class ProceduralEvaluator {
public:
    static void evaluate(ExpressionEvaluator& eval, const StmtPtr& stmt);
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
                                             outs[i].name + "'");
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
        throw std::runtime_error("Unknown procedure: " + call.name);
    }
    
    const auto& proc = it->second;
    if (call.inputArgs.size() != proc.inputs.size()) {
        throw std::runtime_error("Procedure " + call.name + " expected " + std::to_string(proc.inputs.size()) + " inputs, got " + std::to_string(call.inputArgs.size()));
    }
    
    ExpressionEvaluator localEval(numVariables_, coolpropConfig_);
    localEval.userFunctions_ = userFunctions_;
    localEval.userProcedures_ = userProcedures_;
    
    for (size_t i = 0; i < proc.inputs.size(); ++i) {
        const std::string& inputName = proc.inputs[i];
        if (!inputName.empty() && inputName.back() == '$') {
            localEval.setStringVariable(inputName, evaluateString(call.inputArgs[i]));
        } else {
            localEval.setVariable(inputName, evaluate(call.inputArgs[i]));
        }
    }
    
    // Evaluate body
    for (const auto& stmt : proc.body) {
        ProceduralEvaluator::evaluate(localEval, stmt);
    }
    
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
            ADValue v = evaluate(call.args[0]);
            std::ostringstream oss;
            // EES prints integer values without decimal, otherwise up to ~10
            // significant digits, trimming trailing zeros.  Keep it simple and
            // compatible.
            if (std::isfinite(v.value) && v.value == std::trunc(v.value) &&
                std::abs(v.value) < 1e16) {
                oss << static_cast<long long>(v.value);
            } else {
                oss.setf(std::ios::fmtflags(0), std::ios::floatfield);
                oss.precision(10);
                oss << v.value;
            }
            return oss.str();
        }
        
        // PHASE$(Fluid, T=..., P=...) — return the thermodynamic phase as a
        // human-readable string.
        if (lname == "phase$") {
            std::string eesFluidName;
            if (!call.args.empty() && call.args[0]->is<StringLiteral>()) {
                eesFluidName = call.args[0]->as<StringLiteral>().value;
            } else if (!call.args.empty() && call.args[0]->is<Variable>()) {
                const auto& var = call.args[0]->as<Variable>();
                std::string resolvedName = resolveVariableName(var);
                if (hasStringVariable(resolvedName)) {
                    eesFluidName = getStringVariable(resolvedName);
                } else {
                    eesFluidName = var.name;
                    if (!eesFluidName.empty() && eesFluidName.back() == '$') {
                        eesFluidName.pop_back();
                    }
                }
            }
            auto fluid = FluidRegistry::getFluid(eesFluidName);
            if (!fluid) {
                throw std::runtime_error(unknownFluidErrorMessage(eesFluidName));
            }
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

ADValue ExpressionEvaluator::evaluateBinaryOp(const BinaryOp& op) {
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
    }
    
    throw std::runtime_error("Unknown binary operator: " + op.op);
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
    
    // Check if it's a thermodynamic function handled by CoolProp
    bool isThermo = false;
    try {
        getCoolPropOutputParam(name);
        isThermo = true;
    } catch (...) {
        // Not a thermo function
    }
    
    if (!func.namedArgs.empty() || isThermo) {
        return evaluateCoolPropFunction(func);
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
    
    // CONVERT('from', 'to') — returns conversion factor so that 1 [from] = factor [to]
    if (name == "convert" && func.args.size() == 2) {
        std::string fromUnit = evaluateString(func.args[0]);
        std::string toUnit = evaluateString(func.args[1]);
        // Find which UnitType these units belong to, then compute factor
        auto tryConvert = [&](UnitType type) -> double {
            double fromSI = UnitConverter::toSI(1.0, type, fromUnit);
            double toVal = UnitConverter::fromSI(fromSI, type, toUnit);
            return toVal;
        };
        // Try all unit types to find a matching pair
        static const UnitType types[] = {
            UnitType::Energy, UnitType::Pressure, UnitType::Mass, UnitType::Length,
            UnitType::Time, UnitType::Power, UnitType::SpecificHeat, UnitType::SpecificEnergy,
            UnitType::SpecificEntropy, UnitType::Conductivity, UnitType::Viscosity,
            UnitType::Density, UnitType::Dimensionless
        };
        for (auto type : types) {
            double factor = tryConvert(type);
            // If the factor differs from 1.0, we found a real conversion
            if (std::abs(factor - 1.0) > 1e-15) {
                return ADValue::constant(factor, numVariables_);
            }
        }
        // If no conversion found, return 1.0 (same units)
        return ADValue::constant(1.0, numVariables_);
    }
    
    // CONVERTTEMP('from', 'to', value) — converts a temperature value between scales
    if (name == "converttemp" && func.args.size() == 3) {
        std::string fromUnit = evaluateString(func.args[0]);
        std::string toUnit = evaluateString(func.args[1]);
        ADValue val = evaluate(func.args[2]);
        double inSI = UnitConverter::toSI(val.value, UnitType::Temperature, fromUnit);
        double result = UnitConverter::fromSI(inSI, UnitType::Temperature, toUnit);
        return ADValue::constant(result, numVariables_);
    }
    
    std::vector<ADValue> args;
    for (const auto& arg : func.args) {
        args.push_back(evaluate(arg));
    }
    
    // Special case for pi() with 0 arguments
    if (name == "pi" && args.empty()) {
        return ADValue::constant(3.14159265358979323846, numVariables_);
    }

    // ----------------------------------------------------------------
    // Lookup table / interpolation functions
    // ----------------------------------------------------------------
    // These functions have mixed arguments: string literals (table name,
    // column names) followed by numeric arguments.  We re-parse the
    // function call directly rather than using the already-evaluated args.
    static const std::unordered_set<std::string> LOOKUP_FUNC_NAMES = {
        "lookup", "lookupcol", "lookupcol1", "lookupcellempty",
        "tablevalue", "tablevalue#", "tablerun#",
        "interpolate", "interpolate1", "interpolate2", "interpolate2dm",
        "nlookuprows", "nlookupcolumns",
        "sumlookup", "avglookup", "maxlookup", "minlookup", "stddevlookup",
    };
    if (LOOKUP_FUNC_NAMES.count(name)) {
        return evaluateLookupFunction(func);
    }

    return evaluateStandardFunction(func.name, args);
}

ADValue ExpressionEvaluator::evaluateUserFunction(const FunctionDefinition& func, const FunctionCall& call) {
    if (call.args.size() != func.parameters.size()) {
        throw std::runtime_error("Function " + func.name + " expected " + std::to_string(func.parameters.size()) + " arguments, got " + std::to_string(call.args.size()));
    }
    
    // Create a local evaluator for the function body
    ExpressionEvaluator localEval(numVariables_, coolpropConfig_);
    // Inherit user functions and procedures
    localEval.userFunctions_ = userFunctions_;
    localEval.userProcedures_ = userProcedures_;
    // Inherit lookup table store
    localEval.lookupTableStore_ = lookupTableStore_;
    
    for (size_t i = 0; i < call.args.size(); ++i) {
        const std::string& paramName = func.parameters[i];
        if (!paramName.empty() && paramName.back() == '$') {
            localEval.setStringVariable(paramName, evaluateString(call.args[i]));
        } else {
            localEval.setVariable(paramName, evaluate(call.args[i]));
        }
    }
    
    // Evaluate the body procedurally
    for (const auto& stmt : func.body) {
        ProceduralEvaluator::evaluate(localEval, stmt);
    }
    
    // The result is stored in a variable with the function name
    return localEval.getVariable(func.name);
}

ADValue ExpressionEvaluator::evaluateCoolPropFunction(const FunctionCall& func) {
    std::string funcName = func.name;
    std::transform(funcName.begin(), funcName.end(), funcName.begin(), ::tolower);
    
    std::string eesFluidName;
    if (!func.args.empty() && func.args[0]->is<StringLiteral>()) {
        eesFluidName = func.args[0]->as<StringLiteral>().value;
    } else if (!func.args.empty() && func.args[0]->is<Variable>()) {
        const auto& var = func.args[0]->as<Variable>();
        std::string resolvedName = resolveVariableName(var);
        if (hasStringVariable(resolvedName)) {
            eesFluidName = getStringVariable(resolvedName);
        } else {
            // Fallback: use the name itself (minus $ if present)
            eesFluidName = var.name;
            if (!eesFluidName.empty() && eesFluidName.back() == '$') {
                eesFluidName.pop_back();
            }
        }
    }
    
    auto fluid = FluidRegistry::getFluid(eesFluidName);
    if (!fluid) {
        throw std::runtime_error(unknownFluidErrorMessage(eesFluidName));
    }
    
    // Derived thermophysical properties — compute by combining base CoolProp calls
    // entirely in AD space so gradients are correct.
    //   KINEMATICVISCOSITY = viscosity / density          (m^2/s)
    //   THERMALDIFFUSIVITY = conductivity / (density*cp)  (m^2/s)
    if (funcName == "kinematicviscosity" || funcName == "thermaldiffusivity") {
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
        } else if (fluid->getType() != FluidType::Unknown && fluid->getType() != FluidType::Mixture) {
            // For pure fluids or incompressible, molar mass is constant
            if (inputs.empty()) {
                std::string cpFluidName = fluid->getCoolPropName();
                try {
                    double mm = timedPropsSI("M", "T", 300, "P", 101325, cpFluidName);
                    return ADValue::constant(mm, numVariables_);
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
        
        HAParam outParam = getHAParam(funcName);
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
        
        try {
            double result = timedHAPropsSI(outputStr, inputStrs[0], vals[0], inputStrs[1], vals[1], inputStrs[2], vals[2]);
            
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
                    // Generic penalty for other NaN/Inf cases
                    constexpr double PENALTY = 1e4;
                    ADValue output(PENALTY, numVariables_);
                    return output;
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
                
                double rPlus = timedHAPropsSI(outputStr, inputStrs[0], vPlus[0], inputStrs[1], vPlus[1], inputStrs[2], vPlus[2]);
                double rMinus = timedHAPropsSI(outputStr, inputStrs[0], vMinus[0], inputStrs[1], vMinus[1], inputStrs[2], vMinus[2]);
                
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
    if (fluid->getType() == FluidType::Mixture) {
        throw std::runtime_error("Fluid '" + eesFluidName + "' (mixture) is not available. "
            "Use REFPROP or an external correlation.");
    } else if (fluid->getType() == FluidType::Unknown) {
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
    
    if (fluid->getType() == FluidType::IdealGas && inputs.size() == 1) {
        // For ideal gases, properties that depend only on temperature {T, h, u,
        // cp, cv} can be computed from any other member of this set without
        // knowing pressure.  When exactly one input is given, check that both
        // the requested output and the provided input belong to this
        // "thermal-only" set.  If so, inject a dummy pressure so CoolProp gets
        // its required 2-input state; otherwise require explicit pressure.
        auto isThermalOnly = [](const std::string& prop) {
            std::string p = prop;
            for (auto& c : p) c = std::tolower(c);
            return p == "h" || p == "enthalpy" || p == "u" || p == "internalenergy" ||
                   p == "cp" || p == "cv" || p == "specheat" ||
                   p == "t" || p == "temperature";
        };
        
        std::string inputName = inputs.begin()->first;
        if (isThermalOnly(funcName) && isThermalOnly(inputName)) {
            // Both output and input are T-dependent-only → inject dummy pressure
            // Use fluid-specific dummy pressure (e.g. 100 Pa for H2O to avoid
            // liquid phase at room temperature).
            auto idealGas = std::dynamic_pointer_cast<IdealGasFluid>(fluid);
            double dummyP = idealGas ? idealGas->getDummyPressureSI() : 101325.0;
            inputs["p"] = ADValue::constant(UnitConverter::fromSI(dummyP, UnitType::Pressure, units.pressure), numVariables_);
        } else if (!fluid->propertyDependsOnPressure(funcName)) {
            // Output doesn't depend on pressure (old path for h(T), cp(T), etc.)
            auto idealGas = std::dynamic_pointer_cast<IdealGasFluid>(fluid);
            double dummyP = idealGas ? idealGas->getDummyPressureSI() : 101325.0;
            inputs["p"] = ADValue::constant(UnitConverter::fromSI(dummyP, UnitType::Pressure, units.pressure), numVariables_);
        } else {
            throw std::runtime_error("Ideal gas property '" + funcName + "' requires pressure input");
        }
    }
    
    // Handle saturation functions (T_sat, P_sat, etc.) which take only 1 property input
    // We add quality Q=0.5 to force saturation state in CoolProp
    if (inputs.size() == 1) {
        std::string lowerName = funcName;
        std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
        if (lowerName == "t_sat" || lowerName == "p_sat" || lowerName == "tsat" || lowerName == "psat") {
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

    checkThermoInputValueHint(diagnostics_, func.name + "()", units,
                              thermoKindFromCoolProp(input1Param, inputNamesList[0]),
                              inputNamesList[0], val1);
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
    
    // ── CoolProp input sanitization ──────────────────────────────────
    // During Newton iteration, the solver may evaluate trial points with
    // unphysical thermodynamic inputs (negative pressure, temperature below
    // absolute zero).  Instead of letting CoolProp throw (which causes the
    // line-search / trust-region to shrink the step but may leave no viable
    // step direction), we clamp the inputs to physically valid ranges and
    // let the solver see a finite (but poor) residual, keeping the
    // optimization landscape smooth and allowing it to navigate back.
    // Clamping is disabled during final solution verification so that
    // out-of-range inputs are caught as errors.
    if (!disableClamping_) {
    constexpr double P_MIN_SI  = 1000.0;   // 1000 Pa (0.01 bar) — low but CoolProp-safe
    constexpr double T_MIN_SI  = 50.0;     // 50 K   — below triple point of most fluids
    constexpr double RHO_MIN_SI = 1e-4;    // near-vacuum density floor
    
    auto clampInput = [&](CoolProp::parameters param, double& v,
                          std::vector<double>& grad) {
        if (param == CoolProp::iP && v < P_MIN_SI) {
            v = P_MIN_SI;
            // Don't zero gradient: keep the derivative information so the solver
            // has a non-singular Jacobian and a direction back to valid region
        } else if (param == CoolProp::iT && v < T_MIN_SI) {
            v = T_MIN_SI;
        } else if (param == CoolProp::iDmass && v < RHO_MIN_SI) {
            v = RHO_MIN_SI;
        }
    };
    
    clampInput(input1Param, val1, inputValues[0].gradient);
    clampInput(input2Param, val2, inputValues[1].gradient);
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
            result += UnitConverter::fromSI(fluid->getReferenceState().h_offset, outputInfo.unitType, units.specific_energy);
        else if (outputInfo.unitType == UnitType::SpecificEntropy)
            result += UnitConverter::fromSI(fluid->getReferenceState().s_offset, outputInfo.unitType, units.specific_entropy);
        
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
            } catch (...) {
                // AbstractState failed — fall through to PropsSI below
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
}

EvaluationResult BlockEvaluator::evaluate(const std::vector<double>& x,
                                           const std::map<std::string, double>& externalVars,
                                           const std::map<std::string, std::string>& externalStringVars,
                                           bool computeJacobian) {
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

    // Track CALL source lines already processed (to skip sibling equations)
    std::set<int> processedCallLines;
    
    for (size_t eq = 0; eq < equations_.size(); ++eq) {
        const EquationInfo* eqInfo = equations_[eq];
        
        if (eqInfo->procedureCall) {
            // Skip sibling CALL equations already handled by their primary
            if (processedCallLines.count(eqInfo->sourceLine)) {
                continue;
            }
            processedCallLines.insert(eqInfo->sourceLine);

            const auto& call = *eqInfo->procedureCall;
            const size_t nOutputs = call.outputVars.size();

            // Find ALL sibling CALL equations in this block (same sourceLine),
            // sorted by global equation ID.  The IR creates sibling equations
            // with consecutive IDs in output order, so sorting by GID gives the
            // correct output-index-to-equation mapping regardless of the SCC
            // ordering within the block.
            std::vector<std::pair<int, size_t>> siblingGidAndLocal; // (globalId, localIdx)
            for (size_t j = 0; j < equations_.size(); ++j) {
                if (equations_[j]->procedureCall &&
                    equations_[j]->sourceLine == eqInfo->sourceLine) {
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
                    // a sibling CALL equation (same source line).  If not,
                    // the output is "foreign" and must be restored.
                    const EquationInfo* matchedEq = equations_[blockVarIdx];
                    if (!(matchedEq->procedureCall &&
                          matchedEq->sourceLine == eqInfo->sourceLine)) {
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

EvaluationResult SystemEvaluator::evaluateBlock(size_t blockIndex) {
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
    return blockEval.evaluate(x, externalVars, externalStringVars);
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

void ProceduralEvaluator::evaluate(ExpressionEvaluator& eval, const StmtPtr& stmt) {
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
        // In EES procedural blocks, '=' is also treated as assignment
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
        if (cond.value > 0.5) {
            for (const auto& s : ite.thenBranch) evaluate(eval, s);
        } else {
            for (const auto& s : ite.elseBranch) evaluate(eval, s);
        }
    } else if (stmt->is<ProcedureCall>()) {
        const auto& call = stmt->as<ProcedureCall>();
        eval.evaluateProcedureCall(call);
    } else if (stmt->is<Duplicate>()) {
        const auto& dup = stmt->as<Duplicate>();
        int startVal = static_cast<int>(std::round(eval.evaluate(dup.start).value));
        int endVal = static_cast<int>(std::round(eval.evaluate(dup.end).value));
        for (int i = startVal; i <= endVal; ++i) {
            eval.setVariable(dup.iteratorVar, ADValue::constant(static_cast<double>(i), eval.getNumVariables()));
            for (const auto& s : dup.body) evaluate(eval, s);
        }
    } else if (stmt->is<RepeatUntil>()) {
        const auto& ru = stmt->as<RepeatUntil>();
        const int maxIter = 10000;
        for (int iter = 0; iter < maxIter; ++iter) {
            for (const auto& s : ru.body) evaluate(eval, s);
            ADValue cond = eval.evaluate(ru.condition);
            if (cond.value > 0.5) break;
        }
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

ADValue ExpressionEvaluator::evaluateLookupFunction(const FunctionCall& func) {
    std::string name = func.name;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);

    // Helper: extract a string literal argument (table name or column name)
    auto getString = [&](size_t idx) -> std::string {
        if (idx >= func.args.size()) {
            throw std::runtime_error(func.name + "(): missing string argument at position " +
                                     std::to_string(idx));
        }
        return evaluateString(func.args[idx]);
    };

    // Helper: evaluate a numeric argument
    auto getNum = [&](size_t idx) -> ADValue {
        if (idx >= func.args.size()) {
            throw std::runtime_error(func.name + "(): missing numeric argument at position " +
                                     std::to_string(idx));
        }
        return evaluate(func.args[idx]);
    };

    // Helper: resolve table — emits a diagnostic and throws if not found
    auto getTable = [&](const std::string& tableName) -> const LookupTable& {
        if (!lookupTableStore_) {
            std::string msg = func.name + "(): lookup table '" + tableName +
                "' not found (no table store is available in this context)";
            if (diagnostics_) diagnostics_->push(DiagnosticSeverity::Error, "L002", msg, "lookup_table");
            throw std::runtime_error(msg);
        }
        const LookupTable* tbl = lookupTableStore_->get(tableName);
        if (!tbl) {
            std::string msg = func.name + "(): lookup table '" + tableName +
                "' not found. To fix: place '" + tableName +
                ".csv' in the model directory, or load it via the Lookup Tables panel.";
            if (diagnostics_) diagnostics_->push(DiagnosticSeverity::Error, "L002", msg, "lookup_table");
            throw std::runtime_error(msg);
        }
        return *tbl;
    };

    // Helper: resolve a column by name or 1-based number (as double arg)
    // EES allows both "T_K" and a numeric column index.
    auto resolveCol = [&](const LookupTable& tbl, size_t argIdx) -> size_t {
        if (argIdx >= func.args.size()) {
            throw std::runtime_error(func.name + "(): missing column argument");
        }
        const auto& arg = func.args[argIdx];
        if (arg->is<StringLiteral>()) {
            std::string colName = arg->as<StringLiteral>().value;
            size_t col = tbl.columnIndex(colName);
            if (col == 0) {
                throw std::runtime_error(func.name + "(): column '" + colName +
                                         "' not found in table '" + tbl.name() + "'");
            }
            return col;
        } else {
            // Numeric column index
            ADValue v = evaluate(arg);
            return static_cast<size_t>(std::round(v.value));
        }
    };

    // ------------------------------------------------------------------
    // NLOOKUPROWS('table')
    // ------------------------------------------------------------------
    if (name == "nlookuprows") {
        std::string tableName = getString(0);
        const LookupTable& tbl = getTable(tableName);
        return ADValue::constant(static_cast<double>(tbl.numRows()), numVariables_);
    }

    // ------------------------------------------------------------------
    // NLOOKUPCOLUMNS('table')
    // ------------------------------------------------------------------
    if (name == "nlookupcolumns") {
        std::string tableName = getString(0);
        const LookupTable& tbl = getTable(tableName);
        return ADValue::constant(static_cast<double>(tbl.numCols()), numVariables_);
    }

    // ------------------------------------------------------------------
    // LOOKUPCOL('table', 'colname')
    // LOOKUPCOL1('table', 'colname')  — identical to LOOKUPCOL in CoolSolve
    // ------------------------------------------------------------------
    if (name == "lookupcol" || name == "lookupcol1") {
        std::string tableName = getString(0);
        std::string colName   = getString(1);
        const LookupTable& tbl = getTable(tableName);
        size_t col = tbl.columnIndex(colName);
        return ADValue::constant(static_cast<double>(col), numVariables_);
    }

    // ------------------------------------------------------------------
    // LOOKUPCELLEMPTY('table', row, col)
    // Returns 1 if the cell is NaN (empty), 0 otherwise.
    // ------------------------------------------------------------------
    if (name == "lookupcellempty") {
        std::string tableName = getString(0);
        const LookupTable& tbl = getTable(tableName);
        ADValue rowV = getNum(1);
        size_t col   = resolveCol(tbl, 2);
        size_t row   = static_cast<size_t>(std::round(rowV.value));
        double v     = tbl.value(row, col);
        return ADValue::constant(std::isnan(v) ? 1.0 : 0.0, numVariables_);
    }

    // ------------------------------------------------------------------
    // LOOKUP('table', row, col)
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
        std::string tableName = getString(0);
        const LookupTable& tbl = getTable(tableName);
        ADValue rowV = getNum(1);

        // Column: string name → exact integer; numeric → possibly non-integer
        ADValue colV;
        const auto& colArg = func.args[2];
        if (colArg->is<StringLiteral>()) {
            std::string colName = colArg->as<StringLiteral>().value;
            size_t col = tbl.columnIndex(colName);
            if (col == 0) {
                throw std::runtime_error(func.name + "(): column '" + colName +
                                         "' not found in table '" + tbl.name() + "'");
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
        const LookupTable& tbl = getTable(tableName);
        size_t col = resolveCol(tbl, 1);

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
    // ------------------------------------------------------------------
    if (name == "interpolate" || name == "interpolate1") {
        if (func.args.size() < 4) {
            throw std::runtime_error(func.name + "(): expects 4 arguments: "
                                     "(table, xcol, ycol, xval)");
        }
        std::string tableName = getString(0);
        const LookupTable& tbl = getTable(tableName);
        size_t xcol = resolveCol(tbl, 1);
        size_t ycol = resolveCol(tbl, 2);
        ADValue xval = getNum(3);
        return tbl.interpolate1D(xcol, ycol, xval, diagnostics_);
    }

    // ------------------------------------------------------------------
    // INTERPOLATE2('table', 'xcol', 'ycol', 'zcol', xval, yval)
    // INTERPOLATE2DM — same as INTERPOLATE2 in CoolSolve
    // ------------------------------------------------------------------
    if (name == "interpolate2" || name == "interpolate2dm") {
        if (func.args.size() < 6) {
            throw std::runtime_error(func.name + "(): expects 6 arguments: "
                                     "(table, xcol, ycol, zcol, xval, yval)");
        }
        std::string tableName = getString(0);
        const LookupTable& tbl = getTable(tableName);
        size_t xcol = resolveCol(tbl, 1);
        size_t ycol = resolveCol(tbl, 2);
        size_t zcol = resolveCol(tbl, 3);
        ADValue xval = getNum(4);
        ADValue yval = getNum(5);
        return tbl.interpolate2D(xcol, ycol, zcol, xval, yval, diagnostics_);
    }

    // Should not reach here given the dispatch set check above
    throw std::runtime_error("evaluateLookupFunction: unhandled function '" + func.name + "'");
}

}  // namespace coolsolve
