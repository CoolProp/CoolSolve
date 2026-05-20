#pragma once

#include "coolsolve/diagnostic.h"
#include "coolsolve/units.h"
#include <string>

namespace coolsolve {

/// Thermodynamic input categories used for unit-mismatch heuristics.
enum class ThermoInputKind {
    Pressure,
    Temperature,
    Enthalpy,
    Entropy,
    Density,
    SpecificVolume,
    Quality,
    HumidityRatio,
    RelativeHumidity,
};

/// Clear per-run deduplication state (call at the start of each solve).
void resetThermoInputHints();

/**
 * @brief Emit a C002 warning when a numeric input looks like a common unit mistake.
 *
 * @param context  Call site label, e.g. "enthalpy()" or "enthalpy (line 3)".
 * @param line     Source line (0 if unknown).
 */
void checkThermoInputValueHint(DiagnosticCollector* diag,
                               const std::string& context,
                               const UnitSystem& units,
                               ThermoInputKind kind,
                               const std::string& paramName,
                               double modelValue,
                               int line = 0);

std::string coolPropUnitsReminder(const UnitSystem& units);

std::string unknownFluidErrorMessage(const std::string& eesFluidName);

/// Non-empty when @p name is a known misspelling of a built-in thermo function.
std::string suggestThermoFunctionTypo(const std::string& name);

/// Non-empty when @p name looks like a misspelled or non-CoolProp fluid identifier.
std::string suggestFluidNameHint(const std::string& name);

/// Non-empty when a thermo call uses an unquoted non-string variable as the fluid.
std::string unquotedFluidArgumentHint(const std::string& funcName,
                                      const std::string& fluidIdentifier,
                                      bool endsWithDollar);

}  // namespace coolsolve
