// The calculation region of a chart: the T-p box from the spec's limits, the range of every quantity over it,
// and the saturation states the chart needs (lines placed through saturation states, the marked critical and
// triple points, the ticks of the saturation scale, the ends of the lines labelled inside the dome).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): core/limits.py (TpLimits and its
// presets), thermo/limits.py (resolve_tp_box, property_range_SI), thermo/provider.py (chart_domain,
// saturation_anchor_values, saturation_pairs, dome_ends, marked_point_states) and engine/planning.py (Domain).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#pragma once

#include "coolsolve/diagram/fluid.h"

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace coolplot {
struct Unit;
}

namespace coolsolve {
namespace diagram {

/// (min, max) in SI units.
using Range = std::pair<double, double>;

/// One bound of a TpLimits: none (the limit of the equation of state), an absolute SI value, or a factor of
/// the fluid's saturation range (coolplot.core.limits.Relative).
struct Bound {
    enum class Kind { None, Absolute, Relative };
    Kind kind = Kind::None;
    double value = 0.0;
};

/**
 * @brief Temperature bounds in K and pressure bounds in Pa (coolplot.core.limits.TpLimits).
 */
struct TpLimits {
    Bound T_min, T_max, p_min, p_max;

    /// The record of a spec: {"T_min": number, {"relative": f} or null, ...} (spec.diagram._limits_to_dict).
    json toJson() const;
    /// From a spec's "limits": a preset name or the record above. @throws DiagramError
    static TpLimits fromJson(const json& data);
};

/// A preset by name, case-insensitive: "DEF", "ACHP", "ORC", "NONE" (TP_LIMITS). @throws DiagramError
TpLimits tpLimitsPreset(const std::string& name);
/// The names of the presets, in the order the GUI offers them.
const std::vector<std::string>& tpLimitsPresetNames();

/**
 * @brief Absolute bounds on quantities in SI units, for the charts that no T-p box describes: the psychrometric
 *        chart, bounded in dry-bulb temperature T and humidity ratio W, and the Mollier h,x chart, bounded in its
 *        temperature scale Thx and W (coolplot.core.limits.Region).
 */
struct Region {
    std::vector<std::pair<std::string, Range>> bounds;  ///< (quantity key, (low, high)), sorted by key

    /// Whether the region bounds quantity `key`.
    bool has(const std::string& key) const;
    /// The (low, high) bounds of quantity `key`. @throws DiagramError when the region has none
    Range rangeSI(const std::string& key) const;
    /// The record of a spec: {"bounds": {key: [low, high], ...}} (Region.to_dict).
    json toJson() const;
    /// From a spec's "limits" record (Region.from_dict, with Region's checks: two different finite values per
    /// quantity, each quantity bounded once). @throws DiagramError
    static Region fromJson(const json& data);
};

/// The region of a humid-air chart where the spec gives none: `temperature` ("T", or "Thx" for the Mollier chart)
/// from -10 to 50 °C and humidity ratios from 0 to 30 g/kg (the defaults of CoolPlot's PsychrometricChart).
Region defaultHumidAirRegion(const std::string& temperature);

/**
 * @brief How a family of lines is placed through saturation states (coolplot.spec.SaturationAnchor).
 */
struct SaturationAnchor {
    std::string by = "T";             ///< "T" or "p"
    double step_SI = 10.0;            ///< step as a difference in SI units
    double Q = 1.0;                   ///< vapour quality of the saturation states
    std::optional<Range> range_SI;    ///< saturation values to place lines between
    std::string label = "value";      ///< "value" or "saturation" (read by coolplot-core)
};

/// Where a line labelled inside the dome meets the bubble line and the dew line.
struct DomeEnd {
    double value_SI = 0.0;
    State bubble;
    State dew;
};

/**
 * @brief The calculation region of a chart in SI units (coolplot.engine.planning.Domain).
 *
 * `anchored` holds the (saturation value, line value) pairs of the families placed by saturation states,
 * `points` the states of the marked critical and triple points ("critical", "triple_liquid",
 * "triple_vapour"), `scale` the ticks of the saturation scale (saturation value, value on the vertical axis),
 * `dome` the ends of the lines labelled inside the dome.
 */
struct Domain {
    std::array<double, 4> tp_box{};                              ///< T_min, T_max, p_min, p_max
    std::vector<std::pair<std::string, Range>> ranges_SI;        ///< (key, (min, max)) of every quantity found
    std::vector<std::pair<std::string, std::vector<Range>>> anchored;
    bool saturation = false;                                     ///< the chart has a saturation curve
    std::vector<std::pair<std::string, State>> points;
    std::vector<Range> scale;
    std::vector<DomeEnd> dome;

    /// (min, max) of quantity `key` over the region. @throws DiagramError when it could not be determined
    Range rangeSI(const std::string& key) const;
    /// The (saturation value, line value) pairs of the family `key` placed by saturation states. @throws DiagramError
    const std::vector<Range>& anchorPairs(const std::string& key) const;
};

/// More lines than this from one saturation anchor is a step that is too small (provider.MAX_ANCHORED_LINES).
constexpr int MAX_ANCHORED_LINES = 200;

/**
 * @brief Absolute (T_min_K, T_max_K, p_min_Pa, p_max_Pa) for one fluid (limits.resolve_tp_box): relative bounds
 *        refer to the saturation range (`margin` inside its ends), the box never leaves the equation of state.
 *
 * CoolSolve's addition for a fluid without saturation curve (CoolPlot refuses its relative bounds): a relative bound,
 * and a bound left to a model that has no limit there, take the fluid's default box (PropertyProvider::defaultBox).
 * @throws DiagramError for an empty region
 */
std::array<double, 4> resolveTpBox(PropertyProvider& fluid, const TpLimits& limits, double margin);

/// Smallest and largest value of a quantity on the corners of the T-p box (limits.property_range_SI).
/// @throws DiagramError when fewer than two corners give a value
Range propertyRangeSI(PropertyProvider& fluid, const std::string& key, const std::array<double, 4>& tp_box);

/// The T-p box and the range of every state quantity over it (provider.chart_domain); quantities whose range
/// cannot be found are left out. The saturation flag is set from the fluid's critical point.
Domain chartDomain(PropertyProvider& fluid, const TpLimits& limits, double margin);

/**
 * @brief (saturation value, line value) of lines of constant `key` through saturation states at even steps
 *        (provider.saturation_anchor_values). The steps are round in `unit`, the display unit of `by`.
 * @param warnings receives the warning about saturation states CoolProp could not compute
 * @throws DiagramError for a fluid without saturation curve, or a step that would place too many lines
 */
std::vector<Range> saturationAnchorValues(PropertyProvider& fluid, const std::string& key, const std::string& by, double step_SI,
                                          double Q, const std::optional<Range>& range_SI, const std::array<double, 4>& tp_box,
                                          double margin, const coolplot::Unit& unit, std::vector<std::string>& warnings);

/// (saturation value, value of `key`) at the saturation states `by` = value, quality `Q` (provider.saturation_pairs).
std::vector<Range> saturationPairs(PropertyProvider& fluid, const std::string& by, const std::string& key,
                                   const std::vector<double>& values_SI, double Q, std::vector<std::string>& warnings);

/// Where lines of constant `key` ("T" or "p") meet the bubble and the dew line (provider.dome_ends).
std::vector<DomeEnd> domeEnds(PropertyProvider& fluid, const std::string& key, const std::vector<double>& values_SI, double margin);

/// The states of the critical and triple points named in `names` ("critical", "triple"), as
/// (name, state) pairs with only finite values (provider.marked_point_states).
std::vector<std::pair<std::string, State>> markedPointStates(PropertyProvider& fluid, const std::vector<std::string>& names,
                                                             double margin);

}  // namespace diagram
}  // namespace coolsolve
