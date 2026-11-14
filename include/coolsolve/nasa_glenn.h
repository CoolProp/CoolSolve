#pragma once

#include <string>
#include <vector>

namespace coolsolve::nasa_glenn {

/**
 * @file nasa_glenn.h
 * @brief The thermodynamic properties of ideal-gas species from the NASA-Glenn 9-coefficient polynomials.
 *
 * Source of the data: B. J. McBride, M. J. Zehe, S. Gordon, "NASA Glenn Coefficients for Calculating Thermodynamic
 * Properties of Individual Species", NASA/TP-2002-211556 (2002), a US government publication; the coefficients and
 * the headers are those of `data/thermo.inp` of NASA CEA (2021 revision, github.com/nasa/cea, Apache-2.0), copied
 * verbatim in `nasa_glenn.cpp`.
 *
 * For each species and temperature interval, with a1 ... a7 and b1, b2 the coefficients and R the gas constant,
 * @code
 *   cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
 *   H/R  = -a1 T^-1 + a2 ln T + a3 T + a4 T^2/2 + a5 T^3/3 + a6 T^4/4 + a7 T^5/5 + b1
 *   S/R  = -a1 T^-2/2 - a2 T^-1 + a3 ln T + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2     (at 1 bar)
 * @endcode
 *
 * The enthalpy is **anchored** on the heat of formation of the species header: h(298.15 K) is that value exactly
 * (the elements N2, O2, H2, He, Ar: 0), the polynomials giving only the rise from 298.15 K. (The polynomial's own
 * H(298.15 K) differs from the header by up to 2.7 J/mol (N; 2.2 J/mol for CO2), the rounding of the fit.) The entropy is the polynomial's:
 * the standard entropy of the ideal gas at 1 bar, s(298.15 K) = `Species::standardEntropy`.
 *
 * The module is pure C++ with no dependency on CoolProp or on the evaluator, so that it could one day be proposed
 * upstream as a backend of CoolProp for the species it has no equation of state for (OH, H, O, N, NO, NO2 ...). It
 * is evaluated by `src/ideal_gas.cpp` (the ideal-gas fluids of the registry that have no CoolProp fluid) and gives the
 * formation enthalpies and standard entropies of the reference state of every ideal gas of the registry
 * (`src/fluids.cpp`).
 */

/// Molar gas constant [J/mol-K], CODATA 2018 (the value of the standard entropies of `src/fluids.cpp`; EES applies
/// 8.314, 5.6e-5 below it)
constexpr double kGasConstant = 8.314462618;

/// Temperature of the enthalpy of formation and of the standard entropy [K]
constexpr double kReferenceTemperature = 298.15;

/// Pressure of the standard state of the entropy [Pa]
constexpr double kStandardPressure = 1.0e5;

/// Lowest temperature of the range of a species [K]: the polynomials of the species that start at 300 K are
/// extrapolated down to it
constexpr double kLowestTemperature = 200.0;

/// One temperature interval of a species: the 7 coefficients a1 ... a7 and the 2 integration constants b1, b2.
struct Interval {
    double tLow, tHigh;   ///< [K]
    double a[7];
    double b[2];          ///< b1 of the enthalpy [K], b2 of the entropy
};

/// The number of atoms of each element in a molecule of the species (the elements of the equilibrium of CHEM_EQUIL).
struct Composition {
    int C = 0, H = 0, O = 0, N = 0, Ar = 0, He = 0;
};

/// One species of the table (the name is that of the fluid registry: `OH`, `NO2`, `C8H18`, ...).
struct Species {
    const char* name;
    Composition composition;
    double molarMass;           ///< [kg/mol], the header of the record (IUPAC 2005 atomic weights)
    double formationEnthalpy;   ///< heat of formation of the gas at 298.15 K [J/mol], the header of the record
    const Interval* intervals;  ///< the temperature intervals, in ascending order
    int intervalCount;
    double polynomialH298;      ///< H of the polynomial at 298.15 K [J/mol] (the anchor of the enthalpy, set at load)
    double standardEntropy;     ///< s at 298.15 K and 1 bar [J/mol-K] (set at load)

    /// Lowest temperature of the range [K]: 200 K, or the start of the first interval when that is lower
    double tMin() const { return intervals[0].tLow < kLowestTemperature ? intervals[0].tLow : kLowestTemperature; }
    /// Highest temperature of the range [K]: the end of the last interval
    double tMax() const { return intervals[intervalCount - 1].tHigh; }
};

/// The thermodynamic properties of a species at T, molar, SI.
struct Point {
    double cp;       ///< specific heat at constant pressure [J/mol-K]
    double dcp_dT;   ///< dcp/dT [J/mol-K^2]
    double h;        ///< enthalpy, anchored: the heat of formation at 298.15 K [J/mol]
    double s;        ///< entropy at 1 bar [J/mol-K]
};

/// The species `name` (the name of the fluid registry, without regard to case), or nullptr when the table has none.
const Species* find(const std::string& name);

/// Every species of the table.
const std::vector<Species>& allSpecies();

/**
 * @brief The properties of the species at T [K].
 *
 * The range is [Species::tMin(), Species::tMax()]. The first interval is extrapolated below its start (the species
 * that start at 300 K: down to 200 K); the interval of a temperature that is a boundary is the upper one, so the
 * result is deterministic (the polynomials are continuous there to the rounding of the fit; dcp/dT jumps).
 * No allocation.
 * @throws std::out_of_range naming the species and its range in K when T is outside it (also NaN)
 */
Point evaluate(const Species& species, double T);

}  // namespace coolsolve::nasa_glenn
