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

/**
 * @brief `UNITSYSTEM('setting')` of EES: whether that unit setting is the one in use.
 *
 * CoolSolve works in SI on a mass basis with angles in degrees (`SI`, `Mass`, `Deg`
 * are always in use, `Eng`, `Molar`, `Rad` never); temperature, pressure and energy
 * follow `units`. Shared by the `UNITSYSTEM()` function (evaluator) and the
 * `$IF UNITSYSTEM('XX')` directive (parser).
 * @throws std::runtime_error for an unknown setting.
 */
bool unitSettingActive(const std::string& setting, const UnitSystem& units);

/**
 * @brief Whether `name` is the symbol of a unit that `CONVERT` knows (K, C, kPa, bar, kJ, kg, m, s, kW, hp, lbf, ...).
 *
 * Case-sensitive, spelled as in EES (`kW`, not `kw`; `J`, not `j`, which is the usual name of an index), and only
 * the units written as one identifier (not `kg/s`, `J/kg-K`). The parser uses it to recognise `a [kW]`, a unit
 * written after a variable, which has the form of the array element `a[kW]` (CS-BUG-VAR-ANNOTATION). It is the
 * table of unit names of unitConversionFactor.
 */
bool isUnitName(const std::string& name);

/**
 * @brief The factor of the EES function `CONVERT('from', 'to')`: 1 [from] = factor [to].
 *
 * A unit string is a product of unit names with `-` or `*` between them, a quotient with `/` (everything after
 * a `/` is in the denominator, `kJ/kg-K` is kJ/(kg*K)), `^n` for an integer power, parentheses and the number 1
 * (`1/s`): `lbm/hr`, `W/m^2-K`, `Btu/hr-ft^2-R`, `ft^3/min`, `kg/(m-s)`. The names are the SI units with their
 * usual prefixes, the English units of EES and others (see isUnitName); a name is matched with its exact
 * capitalization first (`mPa` and `MPa` differ), else case-insensitively when only one unit fits. Temperature
 * units are temperature differences (K = C, R = F = 5/9 K): `CONVERT('C','F')` is 1.8, the conversion of a
 * temperature value with its offset is CONVERTTEMP. The result is cached per pair of strings (a constant
 * CONVERT is evaluated at every iteration of the solver).
 *
 * @throws std::runtime_error "CONVERT('from','to'): <reason>" for an unknown or ambiguous unit name, a unit string
 *         that cannot be read, or units of different dimensions ("incompatible units (energy vs. pressure)").
 */
double unitConversionFactor(const std::string& from, const std::string& to);

class UnitConverter {
public:
    // Convert value from the specified unit to SI (CoolProp standard)
    // SI units: K, Pa, J, kg, m, s, W, W/m-K, Pa-s, J/kg-K, kg/m^3
    static double toSI(double value, UnitType type, const std::string& unit);
    
    // Convert value from SI (CoolProp standard) to the specified unit
    static double fromSI(double value, UnitType type, const std::string& unit);
};

} // namespace coolsolve
