#pragma once

#include <string>

namespace coolsolve::libr_ashrae {

/**
 * @file libr_ashrae.h
 * @brief The EES routines `X_LIBR`, `H_LIBR`, `T_LIBR` and `P_LIBR`: lithium bromide - water solution
 *        (absorption machines), by the correlations of the ASHRAE Handbook.
 *
 * The routines are the EES API (`X_LIBR('SI', T, P)`, ...) around the correlations of the 1989 ASHRAE
 * Handbook of Fundamentals (the enthalpy fit is that of McNeely), the equations that EES documents in the
 * help pages of its `libr_ashrae` library. With the solution at the temperature t [degC] and the mass
 * percentage X of LiBr [%], and the refrigerant (water) at the temperature t' [degC]:
 *
 *  - Duhring relation: `t = B(X) + t' A(X)`, A and B cubic polynomials in X (45 <= X <= 75 %);
 *  - vapour pressure of the refrigerant: `log10(p) = C + D/T' + E/T'^2`, `T' = t' + 273.15` K, p in kPa;
 *  - enthalpy: `h = Ah(X) + t Bh(X) + t^2 Ch(X)` in kJ/kg, Ah, Bh, Ch quartic polynomials in X
 *    (40 <= X <= 70 %, 15.6 < t < 165.6 degC).
 *
 * Their units are those of the routines, not those of the unit system of the model: T in degC, P in kPa,
 * X in mass %, h in kJ/kg (the models write `X_LIBR('SI', T, P/1000)` and `H_LIBR(...)*1000`). Only the unit
 * string `'SI'` is accepted. The validity ranges above are not enforced (EES evaluates outside them);
 * an error is thrown only where the equations are undefined, with a message naming the routine and its
 * arguments. The derivatives are analytical (the implicit function theorem for `X_LIBR`).
 *
 * Future work. The property models of the absorption pairs are meant to be one module per correlation, side
 * by side, each giving values and analytical partial derivatives and all reached from the one evaluator entry
 * point (`evaluateAbsorptionFunction`, `src/evaluator.cpp`): `libr_ashrae` now; planned `libr_patek_klomfar`
 * for the functions `*_LiBrH2O` of the EES `LiBrH2O` library (Patek & Klomfar 2006; register item
 * `CS-GAP-LIBR-PK`) and `nh3h2o` for the EES procedure `NH3H2O` (Ibrahim & Klein 1993; register item
 * `CS-GAP-NH3H2O`). Until then the table of this module lists those EES routines, and the remaining ASHRAE
 * ones (`V_LIBR`, `Q_LIBR`, `VISC_LIBR`, `COND_LIBR`; `CS-GAP-FLUIDS-ABS`), as known but not implemented: a
 * call fails with a message naming the routine and the register item.
 */

/// Value of a routine and its partial derivatives with respect to its two numeric arguments, in argument order.
struct Result {
    double value;
    double d1;
    double d2;
};

/// A routine of the EES libraries of the absorption pairs. The four implemented ones are `NAME('SI', first,
/// second)`; the others are known EES routines that CoolSolve does not implement yet (`evaluate` is `nullptr`).
struct Routine {
    const char* name;     ///< EES name, e.g. "X_LIBR"
    const char* first;    ///< implemented: name of the first numeric argument, e.g. "T" (for the messages)
    const char* second;   ///< implemented: name of the second numeric argument, e.g. "P"
    Result (*evaluate)(double, double);   ///< `nullptr` for a routine that is not implemented
    const char* notImplemented;   ///< not implemented: the end of the message (what it is, the register item)
};

/// The routine called `lowerName` (lower case), implemented or not, `nullptr` for any other name. Cheap (a test of
/// the end of the name): meant to be asked of every function call.
const Routine* findRoutine(const std::string& lowerName);

/// Throws for a routine that is known but not implemented: "NAME: EES LiBr-H2O routine not implemented in
/// CoolSolve yet (register CS-GAP-…); …".
/// @throws std::runtime_error
void requireImplemented(const Routine& routine);

/// Checks the unit string of a call (case-insensitive). `'SI'` is accepted; `'Eng'` (English units, accepted by
/// EES) and any other string are an error naming the routine.
/// @throws std::runtime_error
void checkUnitSystem(const Routine& routine, const std::string& unit);

/// `X_LIBR('SI', T, P)`: mass % LiBr of the solution at the temperature `t` [degC] in equilibrium with the
/// refrigerant vapour at the pressure `p` [kPa]. Solves the Duhring relation for X on the increasing branch of
/// the cubic that contains the practical range (the interval around 60 % where dt/dX > 0, within 0...100 %).
/// @throws std::runtime_error for p <= 0 or a pressure beyond the vapour-pressure equation, non-finite
///         arguments, and a temperature that no concentration of the branch gives.
Result concentration(double t, double p);

/// `H_LIBR('SI', T, X)`: specific enthalpy [kJ/kg] of the solution at `t` [degC] and `x` [mass %].
/// @throws std::runtime_error for non-finite arguments.
Result enthalpy(double t, double x);

/// `T_LIBR('SI', P, X)`: temperature [degC] of the solution of `x` [mass %] in equilibrium with the vapour at
/// the pressure `p` [kPa].
/// @throws std::runtime_error for p <= 0 or a pressure beyond the vapour-pressure equation, non-finite arguments.
Result temperature(double p, double x);

/// `P_LIBR('SI', T, X)`: vapour pressure [kPa] over the solution of `x` [mass %] at the temperature `t` [degC].
/// @throws std::runtime_error for non-finite arguments, A(X) = 0 or a refrigerant temperature below absolute zero.
Result pressure(double t, double x);

}  // namespace coolsolve::libr_ashrae
