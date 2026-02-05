#pragma once

#include <string>
#include <map>

namespace coolsolve {

// Unit types relevant for thermodynamic calculations
enum class UnitType {
    Temperature,    // K, C, F, R
    Pressure,       // Pa, kPa, bar, MPa, psi, atm
    Energy,         // J, kJ
    Mass,           // kg
    Length,         // m
    Time,           // s
    Power,          // W, kW
    Conductivity,   // W/m-K
    Viscosity,      // Pa-s
    SpecificHeat,   // J/kg-K, kJ/kg-K
    SpecificEnergy, // J/kg, kJ/kg (Enthalpy, Internal Energy)
    SpecificEntropy,// J/kg-K, kJ/kg-K (Entropy)
    Density,        // kg/m^3
    Dimensionless
};

// Configuration for the unit system used in the EES model
struct UnitSystem {
    std::string temperature = "C";
    std::string pressure = "Pa";
    std::string energy = "J";
    std::string mass = "kg";
    std::string length = "m";
    std::string time = "s";
    std::string power = "W";
    std::string conductivity = "W/m-K";
    std::string viscosity = "Pa-s";
    std::string specific_heat = "J/kg-K";
    std::string specific_energy = "J/kg";
    std::string specific_entropy = "J/kg-K";
    std::string density = "kg/m^3";
};

class UnitConverter {
public:
    // Convert value from the specified unit to SI (CoolProp standard)
    // SI units: K, Pa, J, kg, m, s, W, W/m-K, Pa-s, J/kg-K, kg/m^3
    static double toSI(double value, UnitType type, const std::string& unit);
    
    // Convert value from SI (CoolProp standard) to the specified unit
    static double fromSI(double value, UnitType type, const std::string& unit);
};

} // namespace coolsolve
