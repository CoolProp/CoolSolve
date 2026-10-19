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
 * @brief What CoolSolve does with a unit setting of EES that is not in use.
 *
 * Used to warn about a `$UnitSystem` directive, which CoolSolve ignores (its units are fixed).
 * @return empty when `setting` is in use (`unitSettingActive`), else a phrase such as
 *         `pressures are read in Pa` (for `kPa`, `bar`, ...), `energies are read in J`,
 *         `temperatures are read in °C`, `angles are read in degrees`,
 *         `CoolSolve works in SI` (`Eng`) or `CoolSolve works on a mass basis` (`Molar`).
 * @throws std::runtime_error for an unknown setting, as unitSettingActive.
 */
std::string unitSettingReplacement(const std::string& setting, const UnitSystem& units);

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

/**
 * @brief Checks that `unit` is a temperature scale of `CONVERTTEMP`: `C`, `K`, `F`, `R` or the names `celsius`,
 *        `kelvin`, `fahrenheit`, `rankine`.
 *
 * The name is matched like those of `CONVERT` (exact capitalization first, else case-insensitively when only one
 * unit fits, which holds for all the scales); spaces around it are ignored. The parser calls it for the literal
 * units of a `CONVERTTEMP`, so that they are an error of their line.
 *
 * @throws std::runtime_error "CONVERTTEMP: 'unit' is not a temperature unit (expected C, K, F or R, ...)".
 */
void checkTemperatureUnit(const std::string& unit);

/**
 * @brief The EES function `CONVERTTEMP('from', 'to', value)`: a temperature `value` expressed in the scale `from`,
 *        in the scale `to` (100 C = 373.15 K = 212 F = 671.67 R).
 *
 * @param slope when not null, receives d(result)/d(value), the ratio of the degrees of the two scales (1, 1.8,
 *        5/9): the automatic differentiation of the evaluator needs it.
 * @throws std::runtime_error as checkTemperatureUnit for a unit that is not a temperature scale.
 */
double convertTemperature(double value, const std::string& from, const std::string& to, double* slope = nullptr);

class UnitConverter {
public:
    // Convert value from the specified unit to SI (CoolProp standard)
    // SI units: K, Pa, J, kg, m, s, W, W/m-K, Pa-s, J/kg-K, kg/m^3
    static double toSI(double value, UnitType type, const std::string& unit);
    
    // Convert value from SI (CoolProp standard) to the specified unit
    static double fromSI(double value, UnitType type, const std::string& unit);

    /**
     * @brief d(SI value) / d(value in `unit`): the factor of the derivatives of a quantity converted to SI by toSI.
     *
     * The conversions are affine (a temperature in C or F has an offset), so the factor is the slope
     * `toSI(1) - toSI(0)` (1000 for kJ, 1 for C and K, 5/9 for F). It is 1 for `UnitType::Dimensionless` and for a
     * unit that toSI does not recognise (toSI passes such a value through unchanged).
     */
    static double scaleToSI(UnitType type, const std::string& unit);

    /**
     * @brief d(value in `unit`) / d(SI value): the factor of the derivatives of a quantity converted from SI by fromSI.
     *
     * The slope `fromSI(1) - fromSI(0)` of the affine conversion (0.001 for kJ, 1 for C and K, 9/5 for F); the inverse
     * of scaleToSI. It is 1 for `UnitType::Dimensionless` and for a unit that fromSI does not recognise.
     */
    static double scaleFromSI(UnitType type, const std::string& unit);
};

} // namespace coolsolve
