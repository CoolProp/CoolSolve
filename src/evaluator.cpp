#include "coolsolve/evaluator.h"
#include "coolsolve/units.h"
#include "coolsolve/fluids.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "CoolProp.h"
#include "CoolPropLib.h"  // For set_config_bool C-style API
#include "AbstractState.h"
#include "DataStructures.h"
#include "HumidAirProp.h"

namespace coolsolve {

// ============================================================================
// Simple Profiling Instrumentation
// ============================================================================
struct ProfilingStats {
    size_t propsSI_calls = 0;
    size_t hapropsSI_calls = 0;
    size_t block_evaluations = 0;
    size_t expression_evaluations = 0;
    double propsSI_time_ms = 0.0;
    double hapropsSI_time_ms = 0.0;
    double block_eval_time_ms = 0.0;
    double expr_eval_time_ms = 0.0;
    double coolprop_warmup_time_ms = 0.0;
    bool coolprop_warmup_performed = false;
    
    void reset() {
        propsSI_calls = 0;
        hapropsSI_calls = 0;
        block_evaluations = 0;
        expression_evaluations = 0;
        propsSI_time_ms = 0.0;
        hapropsSI_time_ms = 0.0;
        block_eval_time_ms = 0.0;
        expr_eval_time_ms = 0.0;
        coolprop_warmup_time_ms = 0.0;
        coolprop_warmup_performed = false;
    }
    
    std::string toString() const {
        std::ostringstream oss;
        oss << "\n=== CoolSolve Profiling Statistics ===\n";
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
        oss << "Total CoolProp time: " << (propsSI_time_ms + hapropsSI_time_ms) << " ms\n";
        oss << "Non-CoolProp block eval time: " << (block_eval_time_ms - propsSI_time_ms - hapropsSI_time_ms) << " ms\n";
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

ADValue ExpressionEvaluator::evaluateUserFunction(const FunctionDefinition& func, const FunctionCall& call) {
    if (call.args.size() != func.parameters.size()) {
        throw std::runtime_error("Function " + func.name + " expected " + std::to_string(func.parameters.size()) + " arguments, got " + std::to_string(call.args.size()));
    }
    
    // Create a local evaluator for the function body
    ExpressionEvaluator localEval(numVariables_, coolpropConfig_);
    // Inherit user functions and procedures
    localEval.userFunctions_ = userFunctions_;
    localEval.userProcedures_ = userProcedures_;
    
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
        throw std::runtime_error("Unknown fluid: " + eesFluidName);
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
                std::ostringstream oss;
                oss << "HumidAir returned invalid result (NaN or Inf) for " << outputStr 
                    << " with inputs: " 
                    << inputStrs[0] << "=" << vals[0] << ", "
                    << inputStrs[1] << "=" << vals[1] << ", "
                    << inputStrs[2] << "=" << vals[2];
                throw std::runtime_error(oss.str());
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
            throw std::runtime_error("HumidAir error: " + std::string(e.what()));
        }
    }
    
    // --- Standard CoolProp Handling ---
    if (fluid->getType() == FluidType::Incompressible) {
        throw std::runtime_error("Incompressible fluids are not yet supported");
    } else if (fluid->getType() == FluidType::Mixture) {
        throw std::runtime_error("Mixtures are not yet supported");
    } else if (fluid->getType() == FluidType::Unknown) {
        auto unsupported = std::dynamic_pointer_cast<UnsupportedFluid>(fluid);
        throw std::runtime_error("Fluid '" + eesFluidName + "' is not supported: " + (unsupported ? unsupported->getReason() : "Unknown reason"));
    }
    
    std::string cpFluidName = fluid->getCoolPropName();
    
    if (fluid->getType() == FluidType::IdealGas && inputs.size() == 1) {
        if (fluid->propertyDependsOnPressure(funcName)) {
            throw std::runtime_error("Ideal gas property '" + funcName + "' requires pressure input");
        }
        inputs["p"] = ADValue::constant(UnitConverter::fromSI(101325.0, UnitType::Pressure, units.pressure), numVariables_);
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
        throw std::runtime_error("CoolProp functions require exactly 2 input properties, got " + 
                                std::to_string(inputs.size()));
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
        throw std::runtime_error("CoolProp does not support T,H as input pair. Use P,H or T,S instead.");
    }
    
    double val1 = inputValues[0].value;
    double val2 = inputValues[1].value;
    
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
    
    try {
        double result = timedPropsSI(outputStr, input1Str, val1, input2Str, val2, cpFluidName);
        
        if (!std::isfinite(result)) {
            std::ostringstream oss;
            oss << "CoolProp returned invalid result (NaN or Inf) for " << outputStr 
                << "(" << cpFluidName << ") with inputs: " 
                << input1Str << "=" << val1 << ", "
                << input2Str << "=" << val2;
            throw std::runtime_error(oss.str());
        }
        
        if (outputInfo.unitType == UnitType::Temperature) result = UnitConverter::fromSI(result, outputInfo.unitType, units.temperature);
        else if (outputInfo.unitType == UnitType::Pressure) result = UnitConverter::fromSI(result, outputInfo.unitType, units.pressure);
        else if (outputInfo.unitType == UnitType::SpecificEnergy) {
            result = UnitConverter::fromSI(result, outputInfo.unitType, units.specific_energy);
            // Apply reference state offset for enthalpy/internal energy
            result += UnitConverter::fromSI(fluid->getReferenceState().h_offset, outputInfo.unitType, units.specific_energy);
        }
        else if (outputInfo.unitType == UnitType::SpecificEntropy) {
            result = UnitConverter::fromSI(result, outputInfo.unitType, units.specific_entropy);
            // Apply reference state offset for entropy
            result += UnitConverter::fromSI(fluid->getReferenceState().s_offset, outputInfo.unitType, units.specific_entropy);
        }
        else if (outputInfo.unitType == UnitType::Density) result = UnitConverter::fromSI(result, outputInfo.unitType, units.density);
        else if (outputInfo.unitType == UnitType::Viscosity) result = UnitConverter::fromSI(result, outputInfo.unitType, units.viscosity);
        else if (outputInfo.unitType == UnitType::Conductivity) result = UnitConverter::fromSI(result, outputInfo.unitType, units.conductivity);
        else if (outputInfo.unitType == UnitType::SpecificHeat) result = UnitConverter::fromSI(result, outputInfo.unitType, units.specific_heat);
        
        const double relEps = 1e-6;
        const double absEps = 1e-8;
        double step1 = std::max(relEps * std::abs(val1), absEps);
        double step2 = std::max(relEps * std::abs(val2), absEps);
        
        double resultPlus1 = timedPropsSI(outputStr, input1Str, val1 + step1, input2Str, val2, cpFluidName);
        double resultMinus1 = timedPropsSI(outputStr, input1Str, val1 - step1, input2Str, val2, cpFluidName);
        double dResult_dInput1 = (std::isfinite(resultPlus1) && std::isfinite(resultMinus1)) ? (resultPlus1 - resultMinus1) / (2.0 * step1) : 0.0;
        
        double resultPlus2 = timedPropsSI(outputStr, input1Str, val1, input2Str, val2 + step2, cpFluidName);
        double resultMinus2 = timedPropsSI(outputStr, input1Str, val1, input2Str, val2 - step2, cpFluidName);
        double dResult_dInput2 = (std::isfinite(resultPlus2) && std::isfinite(resultMinus2)) ? (resultPlus2 - resultMinus2) / (2.0 * step2) : 0.0;
        
        // If the function was volume(), we got density from CoolProp and need to invert
        // to get specific volume: v = 1/ρ, so dv/dX = -1/ρ² * dρ/dX
        // Note: "v" alone maps to viscosity (line 125), so we only check for "volume"
        std::string lowerFuncName = func.name;
        std::transform(lowerFuncName.begin(), lowerFuncName.end(), lowerFuncName.begin(), ::tolower);
        if (lowerFuncName == "volume") {
            // Get the density value in SI units for gradient calculation
            double densitySI = UnitConverter::toSI(result, UnitType::Density, units.density);
            double invDensity2 = 1.0 / (densitySI * densitySI);
            // Invert result: v = 1/ρ
            result = 1.0 / result;
            // Apply chain rule to gradients: dv/dX = -1/ρ² * dρ/dX
            dResult_dInput1 *= -invDensity2;
            dResult_dInput2 *= -invDensity2;
        }
        
        auto getScale = [&](UnitType t, const std::string& u, bool toSI) {
            if (t == UnitType::Dimensionless) return 1.0;
            if (toSI) return UnitConverter::toSI(1.0, t, u) - UnitConverter::toSI(0.0, t, u);
            return UnitConverter::fromSI(1.0, t, u) - UnitConverter::fromSI(0.0, t, u);
        };
        
        double outScale = getScale(outputInfo.unitType, 
            outputInfo.unitType == UnitType::Temperature ? units.temperature : 
            outputInfo.unitType == UnitType::Pressure ? units.pressure :
            outputInfo.unitType == UnitType::SpecificEnergy ? units.specific_energy :
            outputInfo.unitType == UnitType::SpecificEntropy ? units.specific_entropy :
            outputInfo.unitType == UnitType::Density ? units.density :
            outputInfo.unitType == UnitType::Viscosity ? units.viscosity :
            outputInfo.unitType == UnitType::Conductivity ? units.conductivity :
            outputInfo.unitType == UnitType::SpecificHeat ? units.specific_heat : "", false);
            
        double in1Scale = getScale(type1, 
            type1 == UnitType::Temperature ? units.temperature : 
            type1 == UnitType::Pressure ? units.pressure :
            type1 == UnitType::SpecificEnergy ? units.specific_energy :
            type1 == UnitType::SpecificEntropy ? units.specific_entropy :
            type1 == UnitType::Density ? units.density : "", true);
            
        double in2Scale = getScale(type2, 
            type2 == UnitType::Temperature ? units.temperature : 
            type2 == UnitType::Pressure ? units.pressure :
            type2 == UnitType::SpecificEnergy ? units.specific_energy :
            type2 == UnitType::SpecificEntropy ? units.specific_entropy :
            type2 == UnitType::Density ? units.density : "", true);
            
        dResult_dInput1 *= outScale * in1Scale;
        dResult_dInput2 *= outScale * in2Scale;
        
        ADValue output(result, numVariables_);
        for (size_t i = 0; i < numVariables_; ++i) {
            output.gradient[i] = dResult_dInput1 * inputValues[0].gradient[i] + 
                                 dResult_dInput2 * inputValues[1].gradient[i];
        }
        
        return output;
    } catch (const std::exception& e) {
        throw std::runtime_error("CoolProp error in " + func.name + "(" + cpFluidName + "): input " + input1Str + "=" + std::to_string(val1) + ", " + input2Str + "=" + std::to_string(val2) + ": " + e.what());
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
                                           const std::map<std::string, std::string>& externalStringVars) {
    auto eval_start = std::chrono::high_resolution_clock::now();
    g_profilingStats.block_evaluations++;
    
    if (x.size() != variables_.size()) {
        throw std::runtime_error("State vector size mismatch");
    }
    
    EvaluationResult result;
    result.residuals.resize(equations_.size());
    result.jacobian.resize(equations_.size(), std::vector<double>(variables_.size(), 0.0));
    
    ExpressionEvaluator exprEval(variables_.size(), config_);
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
    
    for (size_t eq = 0; eq < equations_.size(); ++eq) {
        const EquationInfo* eqInfo = equations_[eq];
        
        if (eqInfo->procedureCall) {
            // Store the initial state of the outputs to compute residuals
            std::map<std::string, ADValue> oldOutputs;
            for (const auto& var : eqInfo->procedureCall->outputVars) {
                std::string name = exprEval.resolveVariableName(var);
                oldOutputs[name] = exprEval.getVariable(name);
            }
            
            // Evaluate procedure call - this updates exprEval's variable state
            exprEval.evaluateProcedureCall(*eqInfo->procedureCall);
            
            // The residuals are (old_value - new_value) for each output variable
            // We need to fill as many residuals as there are outputs in the CALL.
            // These should map to the equations in the block.
            for (size_t i = 0; i < eqInfo->procedureCall->outputVars.size(); ++i) {
                if (eq + i >= result.residuals.size()) break;
                
                std::string name = exprEval.resolveVariableName(eqInfo->procedureCall->outputVars[i]);
                ADValue newValue = exprEval.getVariable(name);
                ADValue residual = oldOutputs[name] - newValue;
                result.residuals[eq + i] = residual.value;
                for (size_t var = 0; var < variables_.size(); ++var) {
                    result.jacobian[eq + i][var] = residual.gradient[var];
                }
            }
            // Skip the next N-1 residuals in the loop, where N is the number of outputs
            eq += eqInfo->procedureCall->outputVars.size() - 1;
            continue;
        }

        if (!eqInfo->lhs || !eqInfo->rhs) continue;

        ADValue lhs = exprEval.evaluate(eqInfo->lhs);
        ADValue rhs = exprEval.evaluate(eqInfo->rhs);
        ADValue residual = lhs - rhs;
        result.residuals[eq] = residual.value;
        for (size_t var = 0; var < variables_.size(); ++var) {
            result.jacobian[eq][var] = residual.gradient[var];
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

}  // namespace coolsolve
