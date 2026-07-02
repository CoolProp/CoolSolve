#pragma once
#include "coolsolve/ir.h"
#include <optional>
#include <string>

namespace coolsolve {

// Infers variable properties (fluid, type) from usage in thermodynamic function calls
void inferVariables(IR& ir);

// Calculates and sets initial values for variables without guesses
// Uses inferred properties to calculate values at T=25C, x=0
// Defaults to 1.0 for others
void initializeVariables(IR& ir);

// Compute a thermodynamic property value at a given reference state, in the
// same units that initializeVariables() would store for a variable with the
// given `units` annotation.  This is the building block used by the solver's
// multi-start fallback (roadmap §4.2) to build internally-consistent starting
// vectors: every thermo variable of a block is re-evaluated at the same
// (T_K, P_Pa), so e.g. h and T stay compatible with h = enthalpy(fluid, T, P).
//
// @param inferredFluid  Fluid name as stored in VariableInfo::inferredFluid
//                       (e.g. "Water", "R744", "Air_ha").
// @param propCode       CoolProp-style property code from
//                       VariableInfo::inferredProperty ("T","P","H","S","D",...).
// @param T_K            Temperature in Kelvin.
// @param P_Pa           Pressure in Pascal.
// @param units          The variable's unit annotation (e.g. "J/kg").
// @return The property value in `units`, or nullopt if CoolProp cannot
//         evaluate it at the requested state.
std::optional<double> computeThermoGuessAt(const std::string& inferredFluid,
                                           const std::string& propCode,
                                           double T_K, double P_Pa,
                                           const std::string& units);

} // namespace coolsolve
