#pragma once

namespace coolsolve::psychprops {

/**
 * @file psychprops.h
 * @brief The nine outputs of the EES procedure PSYCHPROPS (help `hs786.htm`), shared by the parser and the evaluator.
 *
 * `CALL psychprops(T=T_1, P=P_1, R=RH_1 : T, v, h, s, u, w, R, WetB, DewP)` returns, in this order, the temperature,
 * specific volume, enthalpy, entropy, internal energy, humidity ratio, relative humidity, wet-bulb and dew-point
 * temperatures of the substance AirH2O. Each output is the AirH2O property function of the same quantity at the
 * state given by the three inputs (P and two others), so the procedure returns exactly what `volume(AirH2O, ...)`,
 * `enthalpy(AirH2O, ...)`, ... return: per kg of dry air.
 */

/// One output of the procedure.
struct Output {
    /// The identifier (lower case) of the same quantity as an input (`T=`, `R=`, ...): an output written with the
    /// variable that is that input is the input itself, no equation (`CALL psychprops(T=T_1, ... : T_1, ...)`).
    const char* input;
    /// The AirH2O property function that computes it.
    const char* function;
};

inline constexpr int kOutputCount = 9;
inline constexpr Output kOutputs[kOutputCount] = {
    {"t", "temperature"}, {"v", "volume"},  {"h", "enthalpy"}, {"s", "entropy"}, {"u", "intenergy"},
    {"w", "humrat"},      {"r", "relhum"},  {"b", "wetbulb"},  {"d", "dewpoint"},
};

/// The inputs of the positional form `CALL psychprops(T, P, R : ...)`, in order (identifiers, lower case).
inline constexpr const char* kPositionalInputs[3] = {"t", "p", "r"};

/// The identifiers (lower case) accepted as inputs of the named form: those of the AirH2O property functions that
/// CoolProp's humid-air model takes as inputs. P is required, and so are two of the others.
inline constexpr const char* kInputIdentifiers[] = {"t", "p", "r", "w", "b", "d", "h", "v", "s"};

}  // namespace coolsolve::psychprops
