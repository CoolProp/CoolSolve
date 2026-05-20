#include "coolsolve/user_hints.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace coolsolve {

namespace {

thread_local std::unordered_set<std::string> g_thermoHintKeys;

static std::string lowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

static bool isPressurePa(const UnitSystem& units) {
    const std::string p = lowerCopy(units.pressure);
    return p.empty() || p == "pa" || p == "pascal";
}

static bool isTemperatureCelsius(const UnitSystem& units) {
    const std::string t = lowerCopy(units.temperature);
    return t.empty() || t == "c" || t == "celsius";
}

static bool isSpecificEnergyJkg(const UnitSystem& units) {
    const std::string u = lowerCopy(units.specific_energy);
    return u.empty() || u == "j/kg";
}

static bool isSpecificEntropyJkgK(const UnitSystem& units) {
    const std::string u = lowerCopy(units.specific_entropy);
    return u.empty() || u == "j/kg-k" || u == "j/kgk";
}

static bool isDensityKgm3(const UnitSystem& units) {
    const std::string u = lowerCopy(units.density);
    return u.empty() || u == "kg/m^3" || u == "kg/m3";
}

static bool nearValue(double v, double target, double tol) {
    return std::isfinite(v) && std::abs(v - target) <= tol;
}

static void emitHint(DiagnosticCollector* diag,
                     const std::string& key,
                     const std::string& message,
                     int line) {
    if (!diag || g_thermoHintKeys.count(key)) return;
    g_thermoHintKeys.insert(key);
    Diagnostic d;
    d.severity = DiagnosticSeverity::Warning;
    d.code = "C002";
    d.message = message;
    d.source = "user_hints";
    d.line = line;
    diag->push(std::move(d));
}

static std::string formatContext(const std::string& context, int line) {
    if (line > 0) return context + " (line " + std::to_string(line) + ")";
    return context;
}

static std::string paramKeyForDedup(const std::string& paramName) {
    std::string k = paramName;
    std::transform(k.begin(), k.end(), k.begin(), ::tolower);
    return k;
}

static void hintPressure(DiagnosticCollector* diag,
                         const std::string& ctx,
                         const std::string& paramName,
                         double v,
                         int line) {
    const std::string loc = formatContext(ctx, line);
    const std::string keyBase = "P:" + paramKeyForDedup(paramName) + ":" + std::to_string(v);

    if (v >= 12.0 && v <= 16.0) {
        std::ostringstream oss;
        oss << loc << ": " << paramName << "=" << v
            << " looks like psi (" << v << " psi ≈ "
            << std::fixed << std::setprecision(0) << (v * 6894.76)
            << " Pa). CoolSolve expects pressure in Pa.";
        emitHint(diag, keyBase + ":psi", oss.str(), line);
        return;
    }

    if (v >= 50.0 && v <= 200.0 && nearValue(v, 101.325, 5.0)) {
        emitHint(diag, keyBase + ":kpa",
                 loc + ": " + paramName + "=" + std::to_string(v) +
                     " is Pa, not kPa. For ~1 atm use P=101325 or P=1e5.",
                 line);
        return;
    }

    if (v > 0.05 && v < 20.0) {
        std::ostringstream oss;
        oss << loc << ": " << paramName << "=" << v
            << " may be MPa (CoolSolve default is Pa). "
            << "If you meant " << v << " MPa, use P="
            << std::scientific << std::setprecision(4) << (v * 1.0e6) << std::fixed << ".";
        emitHint(diag, keyBase + ":mpa", oss.str(), line);
        return;
    }

    if (v > 0.0 && v < 10000.0) {
        std::ostringstream oss;
        oss << loc << ": " << paramName << "=" << v << " is interpreted as " << v << " Pa. ";
        if (nearValue(v, 1.0, 0.2)) {
            oss << "If you meant 1 atm, use P=101325; ";
        }
        oss << "If you meant " << v << " bar, use P="
            << std::scientific << std::setprecision(4) << (v * 1.0e5) << std::fixed
            << "; if you meant kPa, use P=" << (v * 1.0e3) << ".";
        emitHint(diag, keyBase + ":low", oss.str(), line);
    }
}

static void hintTemperature(DiagnosticCollector* diag,
                            const std::string& ctx,
                            const std::string& paramName,
                            double v,
                            int line) {
    const std::string loc = formatContext(ctx, line);
    const std::string keyBase = "T:" + paramKeyForDedup(paramName) + ":" + std::to_string(v);

    if (v >= 200.0 && v <= 600.0) {
        std::ostringstream oss;
        oss << loc << ": " << paramName << "=" << v
            << " is interpreted as " << v
            << " °C. If you meant " << v << " K, use T="
            << std::fixed << std::setprecision(2) << (v - 273.15) << ".";
        emitHint(diag, keyBase + ":kelvin", oss.str(), line);
    }

    const struct { double f; const char* meaning; } fRefs[] = {
        {32.0, "water freezing (0 °C)"},
        {77.0, "room temperature (~25 °C)"},
        {212.0, "water boiling at 1 atm (100 °C)"},
    };
    for (const auto& ref : fRefs) {
        if (nearValue(v, ref.f, 3.0)) {
            const double c = (v - 32.0) * 5.0 / 9.0;
            std::ostringstream oss;
            oss << loc << ": " << paramName << "=" << v
                << " looks like Fahrenheit (" << ref.meaning
                << "). CoolSolve uses Celsius: T="
                << std::fixed << std::setprecision(2) << c << ".";
            emitHint(diag, keyBase + ":f:" + std::to_string(ref.f), oss.str(), line);
            break;
        }
    }
}

static void hintEnthalpy(DiagnosticCollector* diag,
                         const std::string& ctx,
                         const std::string& paramName,
                         double v,
                         int line) {
    if (v <= 0.0 || v >= 5000.0) return;
    const std::string loc = formatContext(ctx, line);
    std::ostringstream oss;
    oss << loc << ": " << paramName << "=" << v
        << " is very small for enthalpy in J/kg (typical liquid water ≈ 2e5–5e5). "
        << "If you meant " << v << " kJ/kg, use H=" << std::fixed << std::setprecision(0)
        << (v * 1000.0) << ".";
    emitHint(diag, "H:" + paramKeyForDedup(paramName) + ":" + std::to_string(v), oss.str(), line);
}

static void hintEntropy(DiagnosticCollector* diag,
                          const std::string& ctx,
                          const std::string& paramName,
                          double v,
                          int line) {
    if (v <= 0.0 || v >= 30.0) return;
    const std::string loc = formatContext(ctx, line);
    std::ostringstream oss;
    oss << loc << ": " << paramName << "=" << v
        << " is small for entropy in J/(kg·K). "
        << "If you meant " << v << " kJ/(kg·K), use S=" << (v * 1000.0) << ".";
    emitHint(diag, "S:" + paramName + ":" + std::to_string(v), oss.str(), line);
}

static void hintDensity(DiagnosticCollector* diag,
                        const std::string& ctx,
                        const std::string& paramName,
                        double v,
                        int line) {
    if (v <= 0.0 || v >= 10.0) return;
    const std::string loc = formatContext(ctx, line);
    std::ostringstream oss;
    oss << loc << ": " << paramName << "=" << v
        << " is low for density in kg/m³ (liquid water ≈ 1000). "
        << "If you meant " << v << " g/cm³, use D=" << std::fixed << std::setprecision(0)
        << (v * 1000.0) << ".";
    emitHint(diag, "D:" + paramName + ":" + std::to_string(v), oss.str(), line);
}

static void hintSpecificVolume(DiagnosticCollector* diag,
                               const std::string& ctx,
                               const std::string& paramName,
                               double v,
                               int line) {
    if (v <= 0.0 || v >= 100.0) return;
    const std::string loc = formatContext(ctx, line);
    std::ostringstream oss;
    oss << loc << ": " << paramName << "=" << v
        << " m³/kg is unusually large. If this is density in kg/m³, use D="
        << std::fixed << std::setprecision(4) << (1.0 / v) << " instead of V/VOLUME.";
    emitHint(diag, "V:" + paramName + ":" + std::to_string(v), oss.str(), line);
}

static void hintQuality(DiagnosticCollector* diag,
                        const std::string& ctx,
                        const std::string& paramName,
                        double v,
                        int line) {
    const std::string loc = formatContext(ctx, line);
    if (v > 1.5 && v <= 100.0) {
        emitHint(diag, "Q:percent:" + std::to_string(v),
                 loc + ": " + paramName + "=" + std::to_string(v) +
                     " looks like quality in percent; CoolSolve expects 0..1 (use " +
                     std::to_string(v / 100.0) + ").",
                 line);
    } else if (v < -0.01 || v > 1.01) {
        emitHint(diag, "Q:range:" + std::to_string(v),
                 loc + ": " + paramName + "=" + std::to_string(v) +
                     " is outside the two-phase quality range [0, 1]. "
                     "Subcooled/superheated states should not use quality as an input.",
                 line);
    }
}

static void hintHumidityRatio(DiagnosticCollector* diag,
                              const std::string& ctx,
                              const std::string& paramName,
                              double v,
                              int line) {
    if (v <= 0.05) return;
    const std::string loc = formatContext(ctx, line);
    std::ostringstream oss;
    oss << loc << ": " << paramName << "=" << v
        << " is high for humidity ratio in kg/kg (typical values ≈ 0.001–0.03). "
        << "If you meant " << v << " g/kg, use W=" << std::scientific
        << std::setprecision(4) << (v / 1000.0) << ".";
    emitHint(diag, "W:" + std::to_string(v), oss.str(), line);
}

static void hintRelativeHumidity(DiagnosticCollector* diag,
                                 const std::string& ctx,
                                 const std::string& paramName,
                                 double v,
                                 int line) {
    if (v <= 1.0) return;
    const std::string loc = formatContext(ctx, line);
    if (v <= 100.0) {
        emitHint(diag, "R:percent:" + std::to_string(v),
                 loc + ": " + paramName + "=" + std::to_string(v) +
                     " looks like relative humidity in percent; CoolSolve expects 0..1 (use " +
                     std::to_string(v / 100.0) + ").",
                 line);
    } else {
        emitHint(diag, "R:high:" + std::to_string(v),
                 loc + ": " + paramName + "=" + std::to_string(v) +
                     " is outside the valid relative-humidity range [0, 1].",
                 line);
    }
}

static const std::unordered_map<std::string, std::string>& thermoFunctionTypos() {
    static const std::unordered_map<std::string, std::string> m = {
        {"enthalphy", "enthalpy"}, {"enthaly", "enthalpy"}, {"enthalpyy", "enthalpy"},
        {"entropi", "entropy"}, {"entopy", "entropy"},
        {"temperture", "temperature"}, {"temperatur", "temperature"},
        {"presure", "pressure"}, {"pressur", "pressure"},
        {"densitiy", "density"}, {"denisty", "density"},
        {"viscocity", "viscosity"}, {"viscosiy", "viscosity"},
        {"conductivty", "conductivity"},
        {"dewpont", "dewpoint"}, {"dewpoit", "dewpoint"},
        {"relhumidity", "relhum"}, {"relhumid", "relhum"},
        {"humidityratio", "humrat"}, {"humidityrat", "humrat"},
        {"wetbulbtemp", "wetbulb"},
        {"specifcheat", "specheat"}, {"specificheat", "specheat"},
        {"kinematicvisc", "kinematicviscosity"},
        {"psat", "p_sat"}, {"tsat", "t_sat"},
    };
    return m;
}

static const std::unordered_map<std::string, std::string>& fluidNameHints() {
    static const std::unordered_map<std::string, std::string> m = {
        {"watter", "Water"}, {"wate", "Water"}, {"h20", "Water"},
        {"steam", "Water"}, {"vapor", "Water"},
        {"r134", "R134a"}, {"r22", "R22"}, {"r410a", "R410A"},
        {"airh20", "AirH2O"}, {"air-h2o", "AirH2O"},
        {"nh3", "Ammonia"}, {"co2", "CO2"},
    };
    return m;
}

}  // namespace

void resetThermoInputHints() {
    g_thermoHintKeys.clear();
}

void checkThermoInputValueHint(DiagnosticCollector* diag,
                               const std::string& context,
                               const UnitSystem& units,
                               ThermoInputKind kind,
                               const std::string& paramName,
                               double modelValue,
                               int line) {
    if (!diag || !std::isfinite(modelValue)) return;

    switch (kind) {
        case ThermoInputKind::Pressure:
            if (isPressurePa(units)) hintPressure(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::Temperature:
            if (isTemperatureCelsius(units)) hintTemperature(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::Enthalpy:
            if (isSpecificEnergyJkg(units)) hintEnthalpy(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::Entropy:
            if (isSpecificEntropyJkgK(units)) hintEntropy(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::Density:
            if (isDensityKgm3(units)) hintDensity(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::SpecificVolume:
            hintSpecificVolume(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::Quality:
            hintQuality(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::HumidityRatio:
            hintHumidityRatio(diag, context, paramName, modelValue, line);
            break;
        case ThermoInputKind::RelativeHumidity:
            hintRelativeHumidity(diag, context, paramName, modelValue, line);
            break;
    }
}

std::string coolPropUnitsReminder(const UnitSystem& units) {
    return " (CoolSolve default units: T in " + units.temperature +
           ", P in " + units.pressure +
           ", H in " + units.specific_energy +
           ", S in " + units.specific_entropy + ")";
}

std::string unknownFluidErrorMessage(const std::string& eesFluidName) {
    std::string msg = "Unknown fluid: '" + eesFluidName + "'. "
                      "The first argument must be a fluid name (e.g. Water, R134a) or a "
                      "string variable (fluid$).";
    const std::string hint = suggestFluidNameHint(eesFluidName);
    if (!hint.empty()) msg += " " + hint;
    else if (eesFluidName == "fluid" || eesFluidName == "refrigerant" ||
             eesFluidName == "workingfluid") {
        msg += " Did you mean to use a string variable such as fluid$='R134a'?";
    }
    return msg;
}

std::string suggestThermoFunctionTypo(const std::string& name) {
    const auto& m = thermoFunctionTypos();
    const auto it = m.find(lowerCopy(name));
    if (it != m.end()) {
        return "Did you mean '" + it->second + "'?";
    }
    return {};
}

std::string suggestFluidNameHint(const std::string& name) {
    const std::string lower = lowerCopy(name);
    const auto& m = fluidNameHints();
    const auto it = m.find(lower);
    if (it != m.end()) {
        return "Did you mean '" + it->second + "'?";
    }
    if (lower == "air") {
        return "For moist air properties use AirH2O (humrat, wetbulb, …), not 'air'.";
    }
    if (lower == "fluid" || lower == "refrigerant" || lower == "workingfluid" ||
        lower == "coolant" || lower == "working_fluid") {
        return "This looks like a placeholder, not a CoolProp fluid. "
               "Define a string variable (e.g. fluid$='R134a') and pass fluid$.";
    }
    return {};
}

std::string unquotedFluidArgumentHint(const std::string& funcName,
                                      const std::string& fluidIdentifier,
                                      bool endsWithDollar) {
    if (endsWithDollar) return {};
    if (fluidIdentifier.empty()) return {};
    const std::string lower = lowerCopy(fluidIdentifier);
    if (lower == "water" || lower == "r134a" || lower == "airh2o") return {};
    return funcName + "(): '" + fluidIdentifier +
           "' is used as the fluid name without quotes. "
           "Use a quoted literal ('Water') or a string variable (" + fluidIdentifier + "$).";
}

}  // namespace coolsolve
