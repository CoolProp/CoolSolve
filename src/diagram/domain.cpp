// The calculation region of a chart (include/coolsolve/diagram/domain.h).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): core/limits.py (TpLimits and its
// presets), thermo/limits.py (resolve_tp_box, property_range_SI), thermo/provider.py (chart_domain,
// saturation_anchor_values, saturation_pairs, dome_ends, marked_point_states) and engine/planning.py (Domain).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#include "coolsolve/diagram/domain.h"

#include <coolplot/core.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

namespace coolsolve {
namespace diagram {

namespace {

constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

Bound absoluteBound(double value) { return {Bound::Kind::Absolute, value}; }
Bound relativeBound(double factor) { return {Bound::Kind::Relative, factor}; }

json boundToJson(const Bound& bound) {
    if (bound.kind == Bound::Kind::None) return nullptr;
    if (bound.kind == Bound::Kind::Relative) return {{"relative", bound.value}};
    return bound.value;
}

Bound boundFromJson(const json& data, const char* name) {
    if (data.is_null()) return {};
    if (data.is_number()) return absoluteBound(data.get<double>());
    if (data.is_object() && data.contains("relative") && data.at("relative").is_number())
        return relativeBound(data.at("relative").get<double>());
    throw DiagramError(std::string("The limit ") + name + " is a number in SI units, {\"relative\": factor} or null, not " +
                       data.dump() + ".");
}

/// The values of a state that are finite numbers (provider: "if math.isfinite(...)").
State finiteValues(const State& state) {
    State out;
    for (const auto& [key, value] : state)
        if (std::isfinite(value)) out[key] = value;
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// TpLimits (core/limits.py)
// ---------------------------------------------------------------------------

json TpLimits::toJson() const {
    return {{"T_min", boundToJson(T_min)}, {"T_max", boundToJson(T_max)}, {"p_min", boundToJson(p_min)},
            {"p_max", boundToJson(p_max)}};
}

TpLimits TpLimits::fromJson(const json& data) {
    if (data.is_null()) return tpLimitsPreset("DEF");
    if (data.is_string()) return tpLimitsPreset(data.get<std::string>());
    if (!data.is_object()) throw DiagramError("The limits of a chart are a preset name or {T_min, T_max, p_min, p_max}, not " + data.dump() + ".");
    if (data.contains("bounds"))
        throw DiagramError("A chart of a fluid takes TpLimits (bounds on T and p), not a Region.");
    TpLimits limits;
    limits.T_min = boundFromJson(data.value("T_min", json()), "T_min");
    limits.T_max = boundFromJson(data.value("T_max", json()), "T_max");
    limits.p_min = boundFromJson(data.value("p_min", json()), "p_min");
    limits.p_max = boundFromJson(data.value("p_max", json()), "p_max");
    return limits;
}

const std::vector<std::string>& tpLimitsPresetNames() {
    static const std::vector<std::string> names = {"DEF", "ACHP", "ORC", "NONE"};
    return names;
}

TpLimits tpLimitsPreset(const std::string& name) {
    std::string upper = name;
    std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) { return std::toupper(c); });
    TpLimits limits;
    if (upper == "NONE") return limits;
    if (upper == "DEF") {
        limits = {relativeBound(1.01), relativeBound(2.25), relativeBound(1.01), relativeBound(2.25)};
    } else if (upper == "ACHP") {
        limits = {absoluteBound(173.15), absoluteBound(493.15), absoluteBound(0.25e5), relativeBound(2.25)};
    } else if (upper == "ORC") {
        limits = {absoluteBound(273.15), absoluteBound(673.15), absoluteBound(0.25e5), relativeBound(2.25)};
    } else {
        throw DiagramError("Unknown limits preset '" + name + "', expected one of ['ACHP', 'DEF', 'NONE', 'ORC'].");
    }
    return limits;
}

// ---------------------------------------------------------------------------
// Domain (engine/planning.py)
// ---------------------------------------------------------------------------

Range Domain::rangeSI(const std::string& key) const {
    for (const auto& [k, range] : ranges_SI)
        if (k == key) return range;
    throw DiagramError("The range of '" + key + "' could not be determined for this fluid; give the isoline values or "
                       "a range explicitly.");
}

const std::vector<Range>& Domain::anchorPairs(const std::string& key) const {
    for (const auto& [k, pairs] : anchored)
        if (k == key) return pairs;
    throw DiagramError("The calculation region has no values for the lines of constant '" + key +
                       "' placed by saturation states; it was made for another spec.");
}

// ---------------------------------------------------------------------------
// thermo/limits.py
// ---------------------------------------------------------------------------

std::array<double, 4> resolveTpBox(FluidModel& fluid, const TpLimits& limits, double margin) {
    std::optional<Range> saturationT, saturationP;  // looked up only when a bound needs it
    auto absolute = [&](const Bound& bound, const std::string& key, int end, double fallback) {
        double value = fallback;
        if (bound.kind == Bound::Kind::Relative) {
            std::optional<Range>& saturation = key == "T" ? saturationT : saturationP;
            if (!saturation) saturation = fluid.saturationRange(key, margin);
            value = bound.value * (end == 0 ? saturation->first : saturation->second);
        } else if (bound.kind == Bound::Kind::Absolute) {
            value = bound.value;
        }
        if (std::isnan(value))
            throw DiagramError("The " + key + " limit of the plotting region is not a number; give a value, or null for "
                               "the limit of the equation of state.");
        return value;
    };
    double T_lo = absolute(limits.T_min, "T", 0, 0.0);
    double T_hi = absolute(limits.T_max, "T", 1, INF);
    double p_lo = absolute(limits.p_min, "p", 0, 0.0);
    double p_hi = absolute(limits.p_max, "p", 1, INF);
    // Never leave the validity range of the equation of state
    if (std::isfinite(fluid.tMinK())) T_lo = std::max(T_lo, fluid.tMinK());
    if (std::isfinite(fluid.tMaxK())) T_hi = std::min(T_hi, fluid.tMaxK());
    if (std::isfinite(fluid.pMinPa())) p_lo = std::max(p_lo, fluid.pMinPa());
    if (std::isfinite(fluid.pMaxPa())) p_hi = std::min(p_hi, fluid.pMaxPa());
    if (!(T_lo < T_hi && p_lo < p_hi))
        throw DiagramError("Empty plotting region T=[" + coolplot::repr_number(T_lo) + ", " + coolplot::repr_number(T_hi) +
                           "] K, p=[" + coolplot::repr_number(p_lo) + ", " + coolplot::repr_number(p_hi) + "] Pa.");
    return {T_lo, T_hi, p_lo, p_hi};
}

Range propertyRangeSI(FluidModel& fluid, const std::string& key, const std::array<double, 4>& tp_box) {
    if (key == "T") return {tp_box[0], tp_box[1]};
    if (key == "p") return {tp_box[2], tp_box[3]};
    if (key == "Q") return {0.0, 1.0};
    std::vector<double> values;
    for (double T : {tp_box[0], tp_box[1]})
        for (double p : {tp_box[2], tp_box[3]}) {
            try {
                fluid.update("p", p, "T", T);
                const double value = fluid.get(key);
                if (std::isfinite(value)) values.push_back(value);
            } catch (const std::exception&) {  // a corner outside the equation of state does not contribute
            }
        }
    if (values.size() < 2) throw DiagramError("Could not determine the range of '" + key + "' for " + fluid.name() + ".");
    return {*std::min_element(values.begin(), values.end()), *std::max_element(values.begin(), values.end())};
}

// ---------------------------------------------------------------------------
// thermo/provider.py
// ---------------------------------------------------------------------------

Domain chartDomain(FluidModel& fluid, const TpLimits& limits, double margin) {
    Domain domain;
    domain.tp_box = resolveTpBox(fluid, limits, margin);
    for (const std::string& key : stateKeys()) {
        try {
            domain.ranges_SI.push_back({key, propertyRangeSI(fluid, key, domain.tp_box)});
        } catch (const DiagramError&) {  // not defined for this fluid
        }
    }
    domain.saturation = fluid.hasCriticalPoint();
    return domain;
}

std::vector<Range> saturationAnchorValues(FluidModel& fluid, const std::string& key, const std::string& by, double step_SI,
                                          double Q, const std::optional<Range>& range_SI, const std::array<double, 4>& tp_box,
                                          double margin, const coolplot::Unit& unit, std::vector<std::string>& warnings) {
    if (!fluid.hasCriticalPoint())
        throw DiagramError("Lines are placed by saturation states only for fluids with a saturation curve; " + fluid.name() +
                           " has none.");
    auto [lo, hi] = fluid.saturationRange(by, margin);
    const double box_lo = by == "T" ? tp_box[0] : tp_box[2], box_hi = by == "T" ? tp_box[1] : tp_box[3];
    lo = std::max(lo, box_lo);
    hi = std::min(hi, box_hi);
    if (range_SI) {
        lo = std::max(lo, range_SI->first);
        hi = std::min(hi, range_SI->second);
    }
    if (!(hi >= lo)) return {};
    // Round steps in the display unit: a difference converts with the scale only; 10 digits drop the noise of the
    // conversion (20 deg F is 11.1111111111 K, which converts back to 19.99999999998)
    const double step_shown = coolplot::parse_number(coolplot::general_text(std::fabs(step_SI * unit.scale), 10));
    const double lo_shown = std::min(unit.from_SI(lo), unit.from_SI(hi));
    const double hi_shown = std::max(unit.from_SI(lo), unit.from_SI(hi));
    const long long first = static_cast<long long>(std::ceil(lo_shown / step_shown - 1e-9));
    const long long last = static_cast<long long>(std::floor(hi_shown / step_shown + 1e-9));
    if (last - first + 1 > MAX_ANCHORED_LINES)
        throw DiagramError("A step of " + coolplot::general_text(step_shown, 6) + " " +
                           (unit.difference_label.empty() ? unit.label : unit.difference_label) + " would place " +
                           std::to_string(last - first + 1) + " lines; choose a larger step.");
    std::vector<Range> values;
    int failed = 0;
    for (long long k = first; k <= last; ++k) {
        const double shown = coolplot::parse_number(coolplot::general_text(static_cast<double>(k) * step_shown, 10));
        const double saturation_SI = std::min(std::max(unit.to_SI(shown), lo), hi);
        if (key == by) {
            values.push_back({saturation_SI, saturation_SI});
            continue;
        }
        double value = NaN;
        try {
            value = fluid.statePoint(by, saturation_SI, "Q", Q).at(key);
        } catch (const std::exception&) {
        }
        if (std::isfinite(value)) values.push_back({saturation_SI, value});
        else ++failed;  // CoolProp cannot compute this saturation state; the line is left out
    }
    if (failed)
        warnings.push_back(std::to_string(failed) + " of " + std::to_string(last - first + 1) +
                           " saturation states for the lines of constant " + key + " of " + fluid.name() +
                           " could not be computed; they have no line.");
    return values;
}

std::vector<Range> saturationPairs(FluidModel& fluid, const std::string& by, const std::string& key,
                                   const std::vector<double>& values_SI, double Q, std::vector<std::string>& warnings) {
    std::vector<Range> pairs;
    int failed = 0;
    for (double value : values_SI) {
        double result = NaN;
        try {
            result = fluid.statePoint(by, value, "Q", Q).at(key);
        } catch (const std::exception&) {
        }
        if (std::isfinite(result)) pairs.push_back({value, result});
        else ++failed;
    }
    if (failed)
        warnings.push_back(std::to_string(failed) + " of " + std::to_string(values_SI.size()) +
                           " saturation states of the saturation scale of " + fluid.name() +
                           " could not be computed; they have no tick.");
    return pairs;
}

std::vector<DomeEnd> domeEnds(FluidModel& fluid, const std::string& key, const std::vector<double>& values_SI, double margin) {
    if (!fluid.hasCriticalPoint()) return {};
    const auto [lo, hi] = fluid.saturationRange(key, margin);
    std::vector<DomeEnd> ends;
    for (double value : values_SI) {
        if (!(lo <= value && value <= hi)) continue;
        try {
            State bubble = finiteValues(fluid.statePoint(key, value, "Q", 0.0));
            State dew = finiteValues(fluid.statePoint(key, value, "Q", 1.0));
            ends.push_back({value, bubble, dew});
        } catch (const std::exception&) {  // the line only loses an extra label
        }
    }
    return ends;
}

std::vector<std::pair<std::string, State>> markedPointStates(FluidModel& fluid, const std::vector<std::string>& names,
                                                             double margin) {
    if (!fluid.hasCriticalPoint()) return {};
    auto marks = [&names](const char* name) { return std::find(names.begin(), names.end(), name) != names.end(); };
    std::vector<std::pair<std::string, State>> states;
    if (marks("critical")) {
        State critical = fluid.critical();
        const double rho = critical.at("rho");
        if (std::isfinite(rho) && rho > 0.0) critical["v"] = 1.0 / rho;
        states.push_back({"critical", finiteValues(critical)});
    }
    // A pseudo-pure blend (R407C in HEOS) is one component to CoolProp, but its lowest temperature is where the
    // fit ends, not a triple point
    const bool pure = fluid.isPure() && fluid.pseudoPureBlend().empty();
    const double T_triple = pure ? fluid.tTripleK() : NaN;
    if (marks("triple") && std::isfinite(T_triple)) {
        const double T_min = fluid.tMinK();
        // The saturation curve starts at the triple point only when the equation of state reaches down to it; the
        // states are taken where the curve starts, a margin above the triple temperature
        if (!std::isfinite(T_min) || T_triple >= T_min) {
            const double T_lo = fluid.saturationRange("T", margin).first;
            for (const auto& [name, Q] : {std::pair<const char*, double>{"triple_liquid", 0.0}, {"triple_vapour", 1.0}}) {
                try {
                    states.push_back({name, finiteValues(fluid.statePoint("T", T_lo, "Q", Q))});
                } catch (const std::exception&) {
                }
            }
        }
    }
    return states;
}

}  // namespace diagram
}  // namespace coolsolve
