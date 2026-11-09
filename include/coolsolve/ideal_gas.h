#pragma once

#include <string>
#include "DataStructures.h"

namespace coolsolve {

/**
 * @file ideal_gas.h
 * @brief The ideal gases of the fluid registry (`IdealGasFluid`), evaluated from the ideal-gas part (alpha0) of the
 * CoolProp equation of state of the substance: pv = RT, h and u functions of T only, cp - cv = R, and an entropy
 * s(T, P) = s0(T) - R ln(P/P0).
 *
 * Every function works on the CoolProp basis (SI, per kg, the enthalpy, internal energy and entropy references of
 * CoolProp); the reference offsets of `IdealGasFluid::getReferenceState()` are applied by the evaluator, as for every
 * fluid. The thermodynamic derivatives are analytical: dh/dT = cp0, ds/dT = cp0/T, ds/dP = -R/P, v = RT/P, ...; the
 * inputs that are not T and P go through the implicit-function theorem. The only finite difference is the one of
 * the dilute-gas transport properties in T, for which CoolProp has no derivative either.
 *
 * The CoolProp state is a separate per-thread HEOS state with its phase imposed to gas (no flash, no phase
 * determination, no limit of the equation of state: it is evaluated at its density `P0/(R T)`), never shared with the
 * real-fluid path (an imposed phase would leak into it). The configured CoolProp backend does not apply.
 */

/// Pressure of the standard state of the entropy [Pa] (the 1 bar of the NASA/JANAF tables)
constexpr double kIdealGasStandardPressureSI = 1.0e5;

/// Range of the temperature searched when it is solved from h, u or s (at fixed P or d) [K]
constexpr double kIdealGasSolveTMinSI = 50.0;
constexpr double kIdealGasSolveTMaxSI = 20000.0;

/// True for an input that fixes the temperature of an ideal gas by itself (T, h, u): two of them do not fix the state.
inline bool idealGasFixesTemperature(CoolProp::parameters p) {
    return p == CoolProp::iT || p == CoolProp::iHmass || p == CoolProp::iUmass;
}

/// The temperature-only ideal-gas properties of a substance at T, per kg, SI (CoolProp basis).
struct IdealGasPoint {
    double cp;       ///< specific heat cp0 [J/kg-K]
    double dcp_dT;   ///< dcp0/dT [J/kg-K^2]
    double h;        ///< enthalpy h0 [J/kg]
    double s;        ///< entropy at the standard pressure s0 [J/kg-K]
    double R;        ///< specific gas constant of the equation of state R/M [J/kg-K]
    double M;        ///< molar mass of the equation of state [kg/mol]
};

/// The ideal-gas properties of the CoolProp fluid `coolPropName` at T [K]. Throws std::runtime_error naming the fluid
/// and T when the equation of state has no finite ideal-gas part there (T <= 0, Air below 33 K).
IdealGasPoint idealGasAt(const std::string& coolPropName, double T);

/// One property of an ideal gas with its two partial derivatives.
struct IdealGasValue {
    double value = 0.0;
    double d1 = 0.0;   ///< derivative with respect to the first input at constant second input
    double d2 = 0.0;   ///< derivative with respect to the second input at constant first input
};

/**
 * @brief Property `output` of the ideal gas `coolPropName` at the state fixed by two inputs, SI, CoolProp basis.
 *
 * Inputs: iT, iP, iHmass, iUmass, iDmass, iSmass (any two that are independent: not two of T, h, u, which each fix the
 * temperature). Outputs: iT, iP, iHmass, iUmass, iSmass, iDmass, iCpmass, iCvmass, ispeed_sound, iZ,
 * iisentropic_expansion_coefficient, imolar_mass and the dilute-gas iviscosity, iconductivity, iPrandtl (functions of
 * T only: the pressure is ignored).
 * The state comes from the inputs in closed form, or by a safeguarded Newton iteration on T for h, u (monotone: cp0
 * and cp0 - R are positive) and for s with P or d given, searched between kIdealGasSolveTMinSI and
 * kIdealGasSolveTMaxSI. Throws std::runtime_error with a message naming the fluid when the state or the property is
 * not defined.
 * @param derivatives  false: the derivatives are not computed (zero), for a residual-only evaluation
 */
IdealGasValue idealGasProperty(const std::string& coolPropName, CoolProp::parameters output,
                               CoolProp::parameters input1, double value1,
                               CoolProp::parameters input2, double value2, bool derivatives = true);

}  // namespace coolsolve
