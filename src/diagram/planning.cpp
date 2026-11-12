// Which curves a chart needs (include/coolsolve/diagram/planning.h).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): engine/planning.py (plan and its
// helpers) and the checks of spec/diagram.py and spec/settings.py for the fields read here.
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#include "coolsolve/diagram/planning.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace coolsolve {
namespace diagram {

namespace {

// Mantissas for round values on a logarithmic scale, from coarse to fine (planning._LOG_MANTISSAS)
const std::vector<std::vector<double>>& logMantissas() {
    static const std::vector<std::vector<double>> sets = {
        {1.0},
        {1.0, 2.0, 5.0},
        {1.0, 1.5, 2.0, 3.0, 5.0, 7.0},
        {1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0},
    };
    return sets;
}

// Round steps on a linear scale: 1, 2, 2.5 or 5 times a power of ten (planning._LINEAR_STEPS)
const double LINEAR_STEPS[] = {1.0, 2.0, 2.5, 5.0, 10.0};

/// 10^n as numpy's round computes it (power_of_ten in numpy's calculation.c).
double powerOfTen(int n) {
    static const double p10[] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8};
    if (n < 9) return p10[n];
    double value = 1e9;
    while (n-- > 9) value *= 10.0;
    return value;
}

/// numpy.round(value, decimals) for decimals >= 0: half to even at the last decimal.
double roundDecimals(double value, int decimals) {
    const double f = powerOfTen(decimals);
    return std::nearbyint(value * f) / f;
}

/// m 10^k for the mantissas and every decade from lo to hi, inside the range (planning._log_values).
std::vector<double> logValues(const std::vector<double>& mantissas, double lo, double hi) {
    const int k_lo = static_cast<int>(std::floor(std::log10(lo)));
    const int k_hi = static_cast<int>(std::ceil(std::log10(hi)));
    std::vector<double> values;
    for (int k = k_lo; k <= k_hi; ++k)
        for (double m : mantissas) {
            const double v = coolplot::parse_number(coolplot::general_text(m * std::pow(10.0, k), 6));
            if (lo * (1 - 1e-9) <= v && v <= hi * (1 + 1e-9)) values.push_back(v);
        }
    std::sort(values.begin(), values.end());
    return values;
}

bool sameAsAny(double value, const std::vector<double>& others) {
    return std::any_of(others.begin(), others.end(), [value](double o) { return coolplot::same_number(value, o); });
}

// --- Reading the spec record ------------------------------------------------

double snapped(double value) { return coolplot::snap(value); }

std::optional<Range> rangeFrom(const json& data, const std::string& what, bool snap) {
    if (data.is_null()) return std::nullopt;
    if (!data.is_array() || data.size() != 2) throw DiagramError(what + " is (min, max) or null, not " + data.dump() + ".");
    double a = data.at(0).get<double>(), b = data.at(1).get<double>();
    if (snap) {
        a = snapped(a);
        b = snapped(b);
    }
    return Range{std::min(a, b), std::max(a, b)};
}

std::optional<int> pointsFrom(const json& data, const std::string& what) {
    if (data.is_null()) return std::nullopt;
    const int points = data.get<int>();
    if (points < 2) throw DiagramError(what + " needs at least 2 points, not " + std::to_string(points) + ".");
    return points;
}

CalcSettings calcFrom(const json& data) {
    CalcSettings calc;
    if (data.is_null()) return calc;
    calc.points = data.value("points", calc.points);
    calc.saturation_points = data.value("saturation_points", calc.saturation_points);
    calc.saturation_spacing = data.value("saturation_spacing", calc.saturation_spacing);
    calc.solver = data.value("solver", calc.solver);
    calc.critical_margin = data.value("critical_margin", calc.critical_margin);
    calc.warn_failed_fraction = data.value("warn_failed_fraction", calc.warn_failed_fraction);
    if (calc.points < 2) throw DiagramError("CalcSettings.points has to be at least 2, not " + std::to_string(calc.points) + ".");
    if (calc.saturation_points < 2)
        throw DiagramError("CalcSettings.saturation_points has to be at least 2, not " + std::to_string(calc.saturation_points) + ".");
    if (calc.saturation_spacing != "critical" && calc.saturation_spacing != "linear")
        throw DiagramError("Unknown saturation_spacing '" + calc.saturation_spacing + "', expected one of ('critical', 'linear').");
    if (calc.solver != "auto" && calc.solver != "tracer" && calc.solver != "flash")
        throw DiagramError("Unknown solver '" + calc.solver + "', expected one of ('auto', 'tracer', 'flash').");
    if (!(calc.critical_margin > 0.0 && calc.critical_margin < 0.1))
        throw DiagramError("critical_margin has to be between 0 and 0.1, not " + coolplot::repr_number(calc.critical_margin) + ".");
    if (!(calc.warn_failed_fraction >= 0.0 && calc.warn_failed_fraction <= 1.0))
        throw DiagramError("warn_failed_fraction has to be between 0 and 1.");
    return calc;
}

SaturationAnchor anchorFrom(const json& data) {
    SaturationAnchor anchor;
    anchor.by = coolplot::quantity_key(data.value("by", std::string("T")));
    if (anchor.by != "T" && anchor.by != "p")
        throw DiagramError("Lines can be placed by saturation temperature 'T' or pressure 'p', not '" + anchor.by + "'.");
    // Snapped like isoline values: -40 deg C converted to K is 233.14999999999998
    anchor.step_SI = snapped(data.value("step_SI", 10.0));
    if (!(anchor.step_SI > 0.0))
        throw DiagramError("The step between saturation states has to be positive, not " + coolplot::repr_number(anchor.step_SI) + ".");
    anchor.Q = data.value("Q", 1.0);
    if (!(anchor.Q >= 0.0 && anchor.Q <= 1.0))
        throw DiagramError("The vapour quality of the saturation states has to lie between 0 and 1, not " +
                           coolplot::repr_number(anchor.Q) + ".");
    anchor.label = data.value("label", std::string("value"));
    if (anchor.label != "saturation" && anchor.label != "value")
        throw DiagramError("The label of lines placed by saturation states is one of ('saturation', 'value'), not '" + anchor.label + "'.");
    anchor.range_SI = rangeFrom(data.value("range_SI", json()), "The range of a saturation anchor", true);
    return anchor;
}

IsolineSet isolineSetFrom(const json& data, const coolplot::DiagramType& diagram) {
    IsolineSet iso;
    iso.key = coolplot::quantity_key(data.at("key").get<std::string>());
    if (!diagram.family(iso.key))
        throw DiagramError("Lines of constant '" + iso.key + "' cannot be drawn in a " + diagram.name() + " diagram.");
    const json values = data.value("values_SI", json());
    if (!values.is_null()) iso.values_SI = values.get<std::vector<double>>();
    iso.num = data.value("num", 15);
    if (iso.num < 1) throw DiagramError("Ask for at least one isoline, not " + std::to_string(iso.num) + ".");
    iso.range_SI = rangeFrom(data.value("range_SI", json()), "The range of the lines of constant " + iso.key, false);
    iso.rounding = data.value("rounding", false);
    iso.points = pointsFrom(data.value("points", json()), "An isoline").value_or(0);
    iso.labels = data.value("labels", false);
    const json anchor = data.value("anchor", json());
    if (!anchor.is_null()) {
        if (iso.values_SI || iso.range_SI)
            throw DiagramError("Give either the values or range of the lines or their saturation anchor, not both; the "
                               "anchor has its own range.");
        if (iso.key == "Q") throw DiagramError("Lines of constant vapour quality cannot be placed by saturation states.");
        iso.anchor = anchorFrom(anchor);
    }
    return iso;
}

LineSpec lineSpecFrom(const json& data, const coolplot::DiagramType& diagram) {
    LineSpec line;
    line.name = data.at("name").get<std::string>();
    if (line.name.find_first_not_of(" \t\r\n") == std::string::npos) throw DiagramError("A line needs a name.");
    const std::string suffix = "/label";
    if (line.name.size() >= suffix.size() && line.name.compare(line.name.size() - suffix.size(), suffix.size(), suffix) == 0)
        throw DiagramError("A line name may not end in '/label', not '" + line.name + "'.");
    line.key = coolplot::quantity_key(data.at("key").get<std::string>());
    // Snapped like isoline values, so that the line shares its curve with a family line of the same value
    line.value_SI = snapped(data.at("value_SI").get<double>());
    if (!std::isfinite(line.value_SI))
        throw DiagramError("Line '" + line.name + "': the value has to be a finite number.");
    if (!diagram.family(line.key) && !diagram.is_axis(line.key))
        throw DiagramError("Line '" + line.name + "': lines of constant '" + line.key + "' cannot be drawn in a " +
                           diagram.name() + " diagram.");
    line.points = pointsFrom(data.value("points", json()), "A line").value_or(0);
    return line;
}

}  // namespace

// ---------------------------------------------------------------------------
// The spec (spec/diagram.py, spec/settings.py)
// ---------------------------------------------------------------------------

std::string domeKey(const coolplot::DiagramType& diagram) {
    if (diagram.x != "T" && diagram.y != "T") return "T";
    if (diagram.x != "p" && diagram.y != "p") return "p";
    return "";
}

std::string scaleKey(const coolplot::DiagramType& diagram) {
    const bool y_scaled = diagram.y == "p" || diagram.y == "T", x_scaled = diagram.x == "p" || diagram.x == "T";
    if (!y_scaled || x_scaled)
        throw DiagramError("A saturation scale needs pressure or temperature on the vertical axis and the other quantities on "
                           "the horizontal one; a " + diagram.name() + " diagram has none.");
    return diagram.y == "p" ? "T" : "p";
}

PlanSpec PlanSpec::fromJson(const json& spec) {
    if (!spec.is_object()) throw DiagramError("The chart spec has to be a JSON object (CoolPlot's DiagramSpec record).");
    try {
        PlanSpec out;
        out.diagram = coolplot::DiagramType::parse(spec.value("diagram", std::string("ph")));
        if (out.diagram.kind != "fluid")
            throw DiagramError("Charts of humid air ('" + out.diagram.name() + "') are not drawable yet.");
        out.units = coolplot::UnitSystem::from_json(spec.contains("units") ? spec.at("units") : json("engineering"));
        out.limits = TpLimits::fromJson(spec.value("limits", json("DEF")));
        const json settings = spec.value("settings", json::object());
        out.calc = calcFrom(settings.is_object() ? settings.value("calc", json()) : json());
        if (settings.is_object() && settings.contains("labels") && settings.at("labels").is_object())
            out.dome_labels = settings.at("labels").value("dome", false);
        std::set<std::string> keys;
        for (const json& iso : spec.value("isolines", json::array())) {
            out.isolines.push_back(isolineSetFrom(iso, out.diagram));
            if (!keys.insert(out.isolines.back().key).second)
                throw DiagramError("Each isoline quantity may appear only once, got '" + out.isolines.back().key + "' twice.");
        }
        std::set<std::string> names;
        for (const json& line : spec.value("lines", json::array())) {
            out.lines.push_back(lineSpecFrom(line, out.diagram));
            if (!names.insert(out.lines.back().name).second)
                throw DiagramError("Line names have to be unique, got '" + out.lines.back().name + "' twice.");
        }
        // A record without the field reads as None: the family's values decide (spec.diagram.DiagramSpec.from_dict)
        const json show = spec.value("show_saturation", json());
        if (!show.is_null()) out.show_saturation = show.get<bool>();
        std::vector<std::string> marked;
        for (const json& point : spec.value("marked_points", json::array())) {
            std::string name = point.get<std::string>();
            name.erase(0, name.find_first_not_of(' '));
            name.erase(name.find_last_not_of(' ') + 1);
            std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::tolower(c); });
            if (name != "critical" && name != "triple")
                throw DiagramError("A chart marks the points 'critical', 'triple', not '" + name + "'.");
            marked.push_back(name);
        }
        for (const char* name : {"critical", "triple"})
            if (std::find(marked.begin(), marked.end(), name) != marked.end()) out.marked_points.push_back(name);
        const json scale = spec.value("saturation_scale", json());
        if (!scale.is_null()) {
            SaturationScale s;
            const json step = scale.value("step_SI", json());
            if (!step.is_null()) {
                s.step_SI = snapped(step.get<double>());
                if (!(*s.step_SI > 0.0)) throw DiagramError("The step of the saturation scale has to be positive.");
            }
            s.Q = scale.value("Q", 1.0);
            if (!(s.Q >= 0.0 && s.Q <= 1.0)) throw DiagramError("The vapour quality of the saturation scale has to lie between 0 and 1.");
            scaleKey(out.diagram);  // refuses a chart without a scale
            out.saturation_scale = s;
        }
        return out;
    } catch (const coolplot::Error& e) {
        throw DiagramError(e.what());
    } catch (const json::exception& e) {
        throw DiagramError(std::string("The chart spec is not a valid CoolPlot spec record: ") + e.what());
    }
}

// ---------------------------------------------------------------------------
// Planning (engine/planning.py)
// ---------------------------------------------------------------------------

json CurveRequest::toJson() const {
    return {{"item_id", item_id}, {"key", key},  {"value_SI", value_SI}, {"points", points},
            {"labels", labels},   {"anchor_SI", anchor_SI ? json(*anchor_SI) : json()}};
}

std::string isolineId(const std::string& key, double value_SI) {
    return "iso/" + key + "/" + coolplot::general_text(coolplot::snap(value_SI), 12);
}

std::vector<double> linearGrid(double lo, double hi, int num) {
    if (num <= 1) return {0.5 * (lo + hi)};
    std::vector<double> values(num);
    const double step = (hi - lo) / (num - 1);
    for (int i = 0; i < num; ++i) values[i] = i * step + lo;
    values.back() = hi;
    return values;
}

std::vector<double> logGrid(double lo, double hi, int num) {
    if (num <= 1) return {std::sqrt(lo * hi)};
    const std::vector<double> exponents = linearGrid(std::log10(lo), std::log10(hi), num);
    std::vector<double> values(num);
    for (int i = 0; i < num; ++i) values[i] = std::pow(10.0, exponents[i]);
    values.front() = lo;  // numpy.geomspace keeps the ends exactly
    values.back() = hi;
    return values;
}

std::vector<double> niceLinear(double lo, double hi, int num) {
    if (!(std::isfinite(lo) && std::isfinite(hi)))
        throw DiagramError("Cannot round values over the range (" + coolplot::repr_number(lo) + ", " + coolplot::repr_number(hi) + ").");
    if (!(hi > lo)) return {lo};
    const double raw = (hi - lo) / std::max(num - 1, 1);
    const double power = std::pow(10.0, std::floor(std::log10(raw)));
    // The round step closest to the raw step on a log scale (the first one on a tie, as Python's min)
    double step = LINEAR_STEPS[0] * power;
    for (double m : LINEAR_STEPS)
        if (std::fabs(std::log(m * power / raw)) < std::fabs(std::log(step / raw))) step = m * power;
    const long long first = static_cast<long long>(std::ceil(lo / step - 1e-9));
    const long long last = static_cast<long long>(std::floor(hi / step + 1e-9));
    const int decimals = std::max(0, -static_cast<int>(std::floor(std::log10(step))) + 1);
    std::vector<double> values;
    for (long long k = first; k <= last; ++k) values.push_back(roundDecimals(static_cast<double>(k) * step, decimals));
    return values;
}

std::vector<double> niceLog(double lo, double hi, int num) {
    if (!(hi > lo && lo > 0.0)) return niceLinear(lo, hi, num);
    std::vector<double> best;
    bool found = false;
    for (const std::vector<double>& mantissas : logMantissas()) {
        std::vector<double> values = logValues(mantissas, lo, hi);
        if (values.empty()) continue;
        if (!found || std::abs(static_cast<int>(values.size()) - num) < std::abs(static_cast<int>(best.size()) - num)) {
            best = values;
            found = true;
        }
    }
    if (!found) return niceLinear(lo, hi, num);  // a range narrower than the finest mantissa step
    return best;
}

std::vector<double> roundValues(const coolplot::Unit& unit, double lo, double hi, int num, bool log) {
    // Rounded in display units, so that labels read 25 deg C, not 298.15 K
    const double a = unit.from_SI(lo), b = unit.from_SI(hi);
    const double lo_d = std::min(a, b), hi_d = std::max(a, b);
    std::vector<double> values = log ? niceLog(lo_d, hi_d, num) : niceLinear(lo_d, hi_d, num);
    if (values.empty()) values = {0.5 * (lo_d + hi_d)};  // no round value in a very narrow range
    for (double& v : values) v = unit.to_SI(v);
    return values;
}

std::vector<double> isolineValuesSI(const PlanSpec& spec, const IsolineSet& iso, const Domain& domain) {
    std::vector<double> values;
    if (iso.values_SI) {
        values = *iso.values_SI;
    } else if (iso.anchor) {
        for (const Range& pair : domain.anchorPairs(iso.key)) values.push_back(pair.second);
    } else {
        const auto [lo, hi] = iso.range_SI ? *iso.range_SI : domain.rangeSI(iso.key);
        const bool log = coolplot::get_quantity(iso.key).log_scale && lo > 0.0;
        if (iso.rounding) values = roundValues(spec.units.at(iso.key), lo, hi, iso.num, log);
        else values = log ? logGrid(lo, hi, iso.num) : linearGrid(lo, hi, iso.num);
    }
    std::set<double> unique;
    for (double v : values)
        if (std::isfinite(v)) unique.insert(coolplot::snap(v));
    return {unique.begin(), unique.end()};
}

std::vector<CurveRequest> plan(const PlanSpec& spec, const Domain& domain) {
    const CalcSettings& calc = spec.calc;
    std::optional<coolplot::Family> boundary;
    for (const coolplot::Family& family : spec.diagram.families())
        if (family.on_boundary()) {
            boundary = family;
            break;
        }
    const std::optional<bool>& show = spec.show_saturation;
    bool addSaturation = show == true && domain.saturation && boundary.has_value();
    std::vector<CurveRequest> requests;

    // The lines of the saturation curve the family did not plan itself
    auto saturationRequests = [&](const coolplot::Family& family, const std::vector<double>& planned) {
        const int points = family.key == "Q" ? calc.saturation_points : calc.points;
        for (double b : family.boundary)
            if (!sameAsAny(b, planned)) requests.push_back({isolineId(family.key, b), family.key, b, points, false, std::nullopt});
    };

    for (const IsolineSet& iso : spec.isolines) {
        const coolplot::Family& family = *spec.diagram.family(iso.key);
        const int points = iso.points ? iso.points : (iso.key == "Q" ? calc.saturation_points : calc.points);
        std::map<double, double> saturation;  // line value -> the saturation state that names it (the first one)
        if (iso.anchor)
            for (const Range& pair : domain.anchorPairs(iso.key)) saturation.emplace(coolplot::snap(pair.second), pair.first);
        std::vector<double> planned;
        for (double value : isolineValuesSI(spec, iso, domain)) {
            if (!family.has_line(value)) continue;  // a quality outside 0 to 1 has no line
            if (show == false && sameAsAny(value, family.boundary)) continue;  // the saturation curve is hidden
            planned.push_back(value);
            const auto named = saturation.find(value);
            requests.push_back({isolineId(iso.key, value), iso.key, value, points, iso.labels,
                                named == saturation.end() ? std::nullopt : std::optional<double>(named->second)});
        }
        if (addSaturation && family.key == boundary->key) {  // next to the family's own lines
            saturationRequests(*boundary, planned);
            addSaturation = false;
        }
    }
    if (addSaturation) saturationRequests(*boundary, {});  // the family is not drawn: the saturation curve comes last
    for (const LineSpec& line : spec.lines) {
        if (spec.diagram.is_axis(line.key)) continue;  // straight, drawn by coolplot-core
        const int points = line.points ? line.points : (line.key == "Q" ? calc.saturation_points : calc.points);
        requests.push_back({line.itemId(), line.key, line.value_SI, points, false, std::nullopt});
    }
    return requests;
}

}  // namespace diagram
}  // namespace coolsolve
