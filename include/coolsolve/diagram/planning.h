// Which curves a chart needs, decided without any property calculation: the values of the lines of each
// family (evenly spaced or round), their points, labels, anchors and item ids, the saturation curve and the
// single lines. Also the parts of a CoolPlot spec record that the glue reads (coolplot-core reads the rest).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): engine/planning.py (plan and its
// helpers) and the checks of spec/diagram.py and spec/settings.py (IsolineSet, LineSpec, SaturationScale,
// CalcSettings) for the fields read here.
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#pragma once

#include "coolsolve/diagram/domain.h"

#include <coolplot/core.h>

#include <optional>
#include <string>
#include <vector>

namespace coolsolve {
namespace diagram {

/// How curves are calculated (coolplot.spec.settings.CalcSettings).
struct CalcSettings {
    int points = 250;
    int saturation_points = 250;
    std::string saturation_spacing = "critical";
    std::string solver = "auto";
    double critical_margin = 1e-7;
    double warn_failed_fraction = 0.1;
};

/// One family of lines, values in SI units (coolplot.spec.IsolineSet).
struct IsolineSet {
    std::string key;
    std::optional<std::vector<double>> values_SI;
    int num = 15;
    std::optional<Range> range_SI;
    bool rounding = false;
    int points = 0;  ///< 0: the default of the calculation settings
    bool labels = false;
    std::optional<SaturationAnchor> anchor;
};

/// One single line of constant `key` (coolplot.spec.LineSpec); its value is snapped to 12 digits.
struct LineSpec {
    std::string name;
    std::string key;
    double value_SI = 0.0;
    int points = 0;  ///< 0: the default of the calculation settings

    std::string itemId() const { return "line/" + name; }
};

/// A scale of saturation states next to the vertical axis (coolplot.spec.SaturationScale).
struct SaturationScale {
    std::optional<double> step_SI;
    double Q = 1.0;
};

/// The number of ticks a saturation scale aims for when no step is given (spec.diagram.SCALE_TICKS).
constexpr int SCALE_TICKS = 10;

/**
 * @brief The parts of a CoolPlot chart spec that the domain and the plan read (coolplot.spec.DiagramSpec).
 *
 * `show_saturation` is true, false or unset (a record without the field: the values of the quality family
 * decide, as in CoolPlot before 1.2). The limits of a chart of a fluid are a TpLimits, those of a humid-air chart
 * (kind "humid_air") a Region with bounds on its two axes.
 */
struct PlanSpec {
    coolplot::DiagramType diagram;
    coolplot::UnitSystem units;
    TpLimits limits;                         ///< a chart of a fluid
    Region region;                           ///< a humid-air chart: bounds on its axes (defaultHumidAirRegion without limits)
    CalcSettings calc;
    std::vector<IsolineSet> isolines;
    std::vector<LineSpec> lines;
    std::optional<bool> show_saturation;
    std::vector<std::string> marked_points;  ///< "critical" and/or "triple", in that order
    std::optional<SaturationScale> saturation_scale;
    bool dome_labels = false;                ///< settings.labels.dome

    /**
     * @brief Read a spec record (DiagramSpec.to_dict, fluid and cycles not needed), with CoolPlot's defaults for
     *        the fields it leaves out (DiagramSpec.from_dict) and CoolPlot's checks.
     * @throws DiagramError for a value CoolPlot refuses
     */
    static PlanSpec fromJson(const json& spec);

    /// Whether the chart is a psychrometric or Mollier chart of humid air (the kind of its diagram type).
    bool humidAir() const { return diagram.kind == "humid_air"; }
    /// The limits as a record: the Region of a humid-air chart, else the TpLimits.
    json limitsJson() const { return humidAir() ? region.toJson() : limits.toJson(); }
};

/// The family whose lines are labelled inside the dome: "T", or "p" where T is an axis; empty for p-T (spec.diagram.dome_key).
std::string domeKey(const coolplot::DiagramType& diagram);
/// The quantity of the saturation scale: "T" next to a pressure axis, "p" next to a temperature axis (spec.diagram.scale_key).
/// @throws DiagramError for a chart that cannot have one
std::string scaleKey(const coolplot::DiagramType& diagram);

/// One isoline to compute: constant `key` at `value_SI` with `points` points (planning.CurveRequest).
struct CurveRequest {
    std::string item_id;
    std::string key;
    double value_SI = 0.0;
    int points = 0;
    bool labels = false;
    std::optional<double> anchor_SI;  ///< the saturation T or p the line passes through, for an anchored family

    /// The "request" of a curve of the compose record.
    json toJson() const;
};

/// Stable item id of an isoline: "iso/<key>/<value to 12 significant digits>" (planning.isoline_id).
std::string isolineId(const std::string& key, double value_SI);

/// `num` evenly spaced values from lo to hi, the middle for one (planning.linear_grid, numpy.linspace).
std::vector<double> linearGrid(double lo, double hi, int num);
/// `num` values with a constant ratio from lo to hi, the geometric mean for one (planning.log_grid, numpy.geomspace).
std::vector<double> logGrid(double lo, double hi, int num);
/// Round values with a step of 1, 2, 2.5 or 5 times a power of ten, about `num` of them (planning.nice_linear).
/// @throws DiagramError for a range that is not finite
std::vector<double> niceLinear(double lo, double hi, int num);
/// Round values m 10^k with the mantissa set giving about `num` values (planning.nice_log).
std::vector<double> niceLog(double lo, double hi, int num);
/// About `num` round values between lo and hi (SI), rounded in display units (planning._round_values).
std::vector<double> roundValues(const coolplot::Unit& unit, double lo, double hi, int num, bool log);

/// The constant values of one family in SI units, snapped and without duplicates (planning.isoline_values_SI).
std::vector<double> isolineValuesSI(const PlanSpec& spec, const IsolineSet& iso, const Domain& domain);

/**
 * @brief All the curves a chart needs, in the order of the families of the spec (planning.plan): the lines of
 *        each family, the lines of the saturation curve it did not plan itself (right after the quality family,
 *        or last without one), then the single lines that are not lines of an axis quantity.
 */
std::vector<CurveRequest> plan(const PlanSpec& spec, const Domain& domain);

}  // namespace diagram
}  // namespace coolsolve
