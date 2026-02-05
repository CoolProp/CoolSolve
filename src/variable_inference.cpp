#include "coolsolve/variable_inference.h"
#include "coolsolve/ast.h"
#include "coolsolve/units.h"
#include "coolsolve/fluids.h"
#include "CoolProp.h"
#include "HumidAirProp.h"
#include <iostream>
#include <algorithm>
#include <chrono>
#include <map>
#include <functional>
#include <cmath>
#include <sstream>

// Profiling for variable inference
static size_t g_inference_propsSI_calls = 0;
static double g_inference_propsSI_time_ms = 0.0;
static double g_inference_fluid_lookup_time_ms = 0.0;

namespace coolsolve {

namespace {

    // Helper to extract string from AST node (StringLiteral or Variable)
    std::string extractString(const ExprPtr& expr, const IR& ir) {
        if (!expr) return "";
        if (expr->is<StringLiteral>()) {
            return expr->as<StringLiteral>().value;
        }
        if (expr->is<Variable>()) {
            // Try to find if this variable has a known string value (e.g. F$ = 'Water')
            const auto& var = expr->as<Variable>();
            const auto* info = ir.getVariable(var.flattenedName());
            if (info && info->solutionStringValue) {
                return *info->solutionStringValue;
            }
            // Fallback: assume the variable name itself is the string (e.g. R744)
            // This handles unquoted fluid names common in EES
            return var.flattenedName();
        }
        return "";
    }

    struct PropMapping {
        std::string code;      // "T", "P", "H", etc.
        UnitType unitType;
    };

    struct FluidGuessState {
        double T_target;
        double P_target;
    };

    std::map<std::string, FluidGuessState> fluidGuessCache;

    std::optional<PropMapping> getPropMapping(const std::string& name) {
        std::string n = name;
        std::transform(n.begin(), n.end(), n.begin(), ::tolower);
        
        if (n == "temperature" || n == "t") return PropMapping{"T", UnitType::Temperature};
        if (n == "pressure" || n == "p") return PropMapping{"P", UnitType::Pressure};
        if (n == "enthalpy" || n == "h") return PropMapping{"H", UnitType::SpecificEnergy};
        if (n == "entropy" || n == "s") return PropMapping{"S", UnitType::SpecificEntropy};
        if (n == "internalenergy" || n == "u") return PropMapping{"U", UnitType::SpecificEnergy};
        if (n == "volume" || n == "v") return PropMapping{"D", UnitType::Density}; // Density or Specific Volume
        if (n == "density" || n == "rho") return PropMapping{"D", UnitType::Density};
        if (n == "quality" || n == "x" || n == "q") return PropMapping{"Q", UnitType::Dimensionless};
        if (n == "viscosity" || n == "mu") return PropMapping{"V", UnitType::Viscosity};
        if (n == "conductivity" || n == "k") return PropMapping{"L", UnitType::Conductivity};
        if (n == "specheat" || n == "cp") return PropMapping{"C", UnitType::SpecificHeat};
        
        return std::nullopt;
    }

    std::string getUnitString(UnitType type, const UnitSystem& units) {
        switch (type) {
            case UnitType::Temperature: return units.temperature;
            case UnitType::Pressure: return units.pressure;
            case UnitType::SpecificEnergy: return units.specific_energy;
            case UnitType::SpecificEntropy: return units.specific_entropy;
            case UnitType::Density: return units.density;
            case UnitType::Viscosity: return units.viscosity;
            case UnitType::Conductivity: return units.conductivity;
            case UnitType::SpecificHeat: return units.specific_heat;
            default: return "";
        }
    }

    // recursive visitor to find function calls
    void findFunctionCalls(const ExprPtr& expr, 
                          std::function<void(const FunctionCall&)> callback) {
        if (!expr) return;

        std::visit([&](const auto& node) {
            using T = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<T, FunctionCall>) {
                callback(node);
                // process args
                for (const auto& arg : node.args) findFunctionCalls(arg, callback);
                for (const auto& [name, arg] : node.namedArgs) findFunctionCalls(arg, callback);
            } else if constexpr (std::is_same_v<T, UnaryOp>) {
                findFunctionCalls(node.operand, callback);
            } else if constexpr (std::is_same_v<T, BinaryOp>) {
                findFunctionCalls(node.left, callback);
                findFunctionCalls(node.right, callback);
            } else if constexpr (std::is_same_v<T, Variable>) {
                for (const auto& idx : node.indices) findFunctionCalls(idx, callback);
            }
        }, expr->node);
    }
}

void inferVariables(IR& ir) {
    UnitSystem defaultUnits; // Use default unit system
    // Note: Ideally we should read $UnitSystem directive if parsed, 
    // but for now we use the defaults in units.h (SI/EES standard)

    // 2. Scan for inputs in all function calls
    auto processCall = [&](const FunctionCall& func) {
        auto funcMapping = getPropMapping(func.name);
        if (funcMapping) {
            // Extract fluid
            std::string fluidName = extractString((!func.args.empty()) ? func.args[0] : nullptr, ir);
            if (fluidName.empty()) return;

            // Check named arguments
            for (const auto& [name, expr] : func.namedArgs) {
                if (expr->is<Variable>()) {
                    const auto& var = expr->as<Variable>();
                    auto argMapping = getPropMapping(name);
                    if (argMapping) {
                        VariableInfo* info = const_cast<VariableInfo*>(ir.getVariable(var.flattenedName()));
                        if (info) {
                            info->inferredProperty = argMapping->code;
                            info->inferredFluid = fluidName;
                            if (info->units.empty()) {
                                info->units = getUnitString(argMapping->unitType, defaultUnits);
                            }
                        }
                    }
                }
            }
        }
    };

    for (const auto& eq : ir.getEquations()) {
        
        // Skip if LHS or RHS is null (e.g. procedure calls)
        if (!eq.lhs || !eq.rhs) {
            // If it is a procedure call, we should still scan input arguments for function calls
            if (eq.procedureCall) {
                for (const auto& arg : eq.procedureCall->inputArgs) {
                    findFunctionCalls(arg, processCall);
                }
            }
            continue;
        }
        
        // 1. Check for Output Assignment: Var = ThermoFunc(...)
        if (eq.lhs->is<Variable>() && eq.rhs->is<FunctionCall>()) {
            const auto& lhsVar = eq.lhs->as<Variable>();
            const auto& func = eq.rhs->as<FunctionCall>();
            
            auto mapping = getPropMapping(func.name);
            if (mapping) {
                // It is a thermo function
                std::string fluidName = extractString((!func.args.empty()) ? func.args[0] : nullptr, ir);
                
                // If fluid found, set property for LHS
                if (!fluidName.empty()) {
                     VariableInfo* info = const_cast<VariableInfo*>(ir.getVariable(lhsVar.flattenedName()));
                     if (info) {
                         info->inferredProperty = mapping->code;
                         info->inferredFluid = fluidName;
                         if (info->units.empty()) {
                             info->units = getUnitString(mapping->unitType, defaultUnits);
                         }
                     }
                }
            }
        }

        findFunctionCalls(eq.rhs, processCall);
        findFunctionCalls(eq.lhs, processCall);
    }
}

// Helper to try computing a property value at given conditions
// Returns true if successful, false otherwise
static bool tryComputeProperty(const std::string& cpFluid, const std::string& prop,
                               double T_K, double P_Pa, double& outValue) {
    try {
        auto t1 = std::chrono::high_resolution_clock::now();
        double valSI = CoolProp::PropsSI(prop, "T", T_K, "P", P_Pa, cpFluid);
        auto t2 = std::chrono::high_resolution_clock::now();
        double call_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
        g_inference_propsSI_calls++;
        g_inference_propsSI_time_ms += call_time;
        
        if (call_time > 100) {
            std::cerr << "[VarInfer] SLOW PropsSI " << prop << "(T=" << T_K << ",P=" << P_Pa << "," << cpFluid << "): " << call_time << " ms\n";
        }
        
        if (std::isfinite(valSI)) {
            outValue = valSI;
            return true;
        }
    } catch (...) {
        // CoolProp threw an exception
    }
    return false;
}

void initializeVariables(IR& ir) {
    UnitSystem units; // Defaults
    
    // Reset profiling counters
    g_inference_propsSI_calls = 0;
    g_inference_propsSI_time_ms = 0.0;
    g_inference_fluid_lookup_time_ms = 0.0;
    
    // Collect errors for variables that couldn't be initialized
    std::vector<std::string> initErrors;

    for (auto& [name, info] : const_cast<std::map<std::string, VariableInfo, CaseInsensitiveLess>&>(ir.getVariables())) {
        // Skip if already has value
        if (info.guessValue || info.solutionValue) continue;
        if (info.type != VariableType::Real) continue; // Skip strings

        double value = 1.0; // Default fallback

        if (!info.inferredFluid.empty() && !info.inferredProperty.empty()) {
            auto t1 = std::chrono::high_resolution_clock::now();
            auto fluid = FluidRegistry::getFluid(info.inferredFluid);
            auto t2 = std::chrono::high_resolution_clock::now();
            g_inference_fluid_lookup_time_ms += std::chrono::duration<double, std::milli>(t2 - t1).count();
            
            if (!fluid) {
                initErrors.push_back("Variable '" + name + "': unknown fluid '" + info.inferredFluid + "'");
                info.guessValue = value;
                continue;
            }
            
            std::string cpFluid = fluid->getCoolPropName();
            bool valueFound = false;
            
            if (fluid->getType() == FluidType::HumidAir) {
                // Humid Air: T=25°C (298.15K), P=1 atm, R=0.5
                double T = 298.15;
                double P = 101325.0;
                double R = 0.5;
                
                std::string outProp = info.inferredProperty;
                if (outProp == "D") outProp = "V"; // HAPropsSI uses V for volume
                
                try {
                    double valSI = HumidAir::HAPropsSI(outProp, "T", T, "P", P, "R", R);
                    
                    if (std::isfinite(valSI)) {
                        auto propMap = getPropMapping(info.inferredProperty);
                        if (propMap) {
                            if (propMap->code == "D" && outProp == "V") {
                                if (valSI != 0) valSI = 1.0 / valSI;
                            }
                            value = UnitConverter::fromSI(valSI, propMap->unitType, info.units);
                        } else {
                            value = valSI;
                        }
                        valueFound = true;
                    }
                } catch (...) {
                    // HAPropsSI failed
                }
                
                if (!valueFound) {
                    initErrors.push_back("Variable '" + name + "': could not compute " + info.inferredProperty + 
                                       " for humid air at reference conditions");
                }
            } else {
                // Real/Ideal Fluid - Simple approach:
                // 1. First try standard conditions: T=20°C (293.15K), P=101325 Pa
                // 2. If that fails, try lower pressure: T=20°C, P=10000 Pa
                // 3. If that fails, try higher temperature: T=100°C (373.15K), P=101325 Pa
                
                double valSI = 0.0;
                
                // Attempt 1: Standard conditions (T=20°C, P=1 atm)
                if (tryComputeProperty(cpFluid, info.inferredProperty, 293.15, 101325.0, valSI)) {
                    valueFound = true;
                }
                // Attempt 2: Lower pressure (useful for gases that might not exist at 1 atm)
                else if (tryComputeProperty(cpFluid, info.inferredProperty, 293.15, 10000.0, valSI)) {
                    valueFound = true;
                }
                // Attempt 3: Higher temperature (useful for fluids that are solid at 20°C)
                else if (tryComputeProperty(cpFluid, info.inferredProperty, 373.15, 101325.0, valSI)) {
                    valueFound = true;
                }
                // Attempt 4: Very low pressure (useful for cryogenic fluids like N2, O2)
                else if (tryComputeProperty(cpFluid, info.inferredProperty, 293.15, 1000.0, valSI)) {
                    valueFound = true;
                }
                
                if (valueFound) {
                    auto propMap = getPropMapping(info.inferredProperty);
                    if (propMap) {
                        if (propMap->unitType == UnitType::SpecificEnergy) {
                            valSI += fluid->getReferenceState().h_offset;
                        }
                        if (propMap->unitType == UnitType::SpecificEntropy) {
                            valSI += fluid->getReferenceState().s_offset;
                        }
                        value = UnitConverter::fromSI(valSI, propMap->unitType, info.units);
                    } else {
                        value = valSI;
                    }
                } else {
                    initErrors.push_back("Variable '" + name + "': could not compute " + info.inferredProperty + 
                                       " for " + cpFluid + " at any reference condition");
                }
            }
            
            // Final check: ensure the computed value is finite
            if (!std::isfinite(value)) {
                initErrors.push_back("Variable '" + name + "': computed value is not finite (inf or NaN) for " + 
                                   info.inferredProperty + " of " + cpFluid);
                value = 1.0; // Reset to default
            }
        }
        
        info.guessValue = value;
    }
    
    // If there were initialization errors, throw an exception to stop the process
    if (!initErrors.empty()) {
        std::ostringstream oss;
        oss << "Variable initialization failed for " << initErrors.size() << " variable(s):\n";
        for (const auto& err : initErrors) {
            oss << "  - " << err << "\n";
        }
        throw std::runtime_error(oss.str());
    }
}

} // namespace coolsolve
