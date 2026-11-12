// Tests of the Diagrams module (src/diagram/, docs/coolplot_integration.md): planning against CoolPlot's Python
// planning, the domain, the CoolPlot fixtures of tools/coolplot_fixtures.py (tests/fixtures/coolplot/), charts
// end to end, cycles from solved values, the invariants of §6.2 and the errors of a chart request.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "coolsolve/diagram/engine.h"
#include "coolsolve/fluids.h"

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"

#include <coolplot/compose.h>
#include <coolplot/svg.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <regex>

using namespace coolsolve;
using namespace coolsolve::diagram;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinRel;
namespace fs = std::filesystem;

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

bool sameNumbers(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!coolplot::same_number(a[i], b[i])) return false;
    return true;
}

/// Item ids, points and labels of a plan, as text for a readable comparison.
std::vector<std::string> describe(const std::vector<CurveRequest>& requests) {
    std::vector<std::string> out;
    for (const CurveRequest& r : requests)
        out.push_back(r.item_id + " " + std::to_string(r.points) + (r.labels ? " labels" : ""));
    return out;
}

/// The hand-made domain of the reference cases (no property call).
Domain referenceDomain() {
    Domain domain;
    domain.tp_box = {250.0, 450.0, 1.5e4, 4.2e6};
    domain.ranges_SI = {{"T", {250.0, 450.0}}, {"p", {1.5e4, 4.2e6}},   {"h", {1.6e5, 5.2e5}}, {"s", {850.0, 2100.0}},
                        {"rho", {0.6, 1400.0}}, {"u", {1.5e5, 4.6e5}}, {"Q", {0.0, 1.0}}};
    domain.saturation = true;
    return domain;
}

fs::path fixtureDir() {
    const fs::path here = fs::path(__FILE__).parent_path() / "fixtures" / "coolplot";
    if (fs::exists(here)) return here;
    return fs::path("..") / "tests" / "fixtures" / "coolplot";
}

json readJson(const fs::path& path) {
    std::ifstream file(path);
    REQUIRE(file.good());
    return json::parse(file);
}

/// A position on an axis for comparisons: log10 on a log axis.
double axisValue(double v, bool log) { return log ? std::log10(v) : v; }

/// Largest distance from the valid points of `a` to the polyline `b`, in units of the axis spans (x and y
/// normalised by `xSpan`, `ySpan` after axisValue); a point of `b` alone counts as a point.
double polylineDistance(const std::vector<double>& ax, const std::vector<double>& ay, const std::vector<double>& bx,
                        const std::vector<double>& by, bool xLog, bool yLog, double xSpan, double ySpan) {
    double worst = 0.0;
    for (size_t i = 0; i < ax.size(); ++i) {
        if (!(std::isfinite(ax[i]) && std::isfinite(ay[i]))) continue;
        const double px = axisValue(ax[i], xLog) / xSpan, py = axisValue(ay[i], yLog) / ySpan;
        double best = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < bx.size(); ++j) {
            if (!(std::isfinite(bx[j]) && std::isfinite(by[j]))) continue;
            const double x0 = axisValue(bx[j], xLog) / xSpan, y0 = axisValue(by[j], yLog) / ySpan;
            best = std::min(best, std::hypot(px - x0, py - y0));
            if (j + 1 < bx.size() && std::isfinite(bx[j + 1]) && std::isfinite(by[j + 1])) {
                const double x1 = axisValue(bx[j + 1], xLog) / xSpan, y1 = axisValue(by[j + 1], yLog) / ySpan;
                const double dx = x1 - x0, dy = y1 - y0, len2 = dx * dx + dy * dy;
                if (len2 > 0.0) {
                    const double t = std::clamp(((px - x0) * dx + (py - y0) * dy) / len2, 0.0, 1.0);
                    best = std::min(best, std::hypot(px - (x0 + t * dx), py - (y0 + t * dy)));
                }
            }
        }
        worst = std::max(worst, best);
    }
    return worst;
}

std::vector<double> numbers(const json& values) {
    std::vector<double> out;
    for (const json& v : values) out.push_back(v.is_null() ? NaN : v.get<double>());
    return out;
}

/// Whether two item ids are the same up to the last digit of their values: a family spaced evenly over a range
/// computed by CoolProp ("iso/h/420287.063156") can round its 12th digit the other way with another CoolProp.
bool sameId(const std::string& a, const std::string& b) {
    if (a == b) return true;
    auto parts = [](const std::string& id) {
        std::vector<std::string> out;
        size_t start = 0;
        for (size_t slash; (slash = id.find('/', start)) != std::string::npos; start = slash + 1) out.push_back(id.substr(start, slash - start));
        out.push_back(id.substr(start));
        return out;
    };
    const std::vector<std::string> pa = parts(a), pb = parts(b);
    if (pa.size() != pb.size()) return false;
    for (size_t i = 0; i < pa.size(); ++i) {
        if (pa[i] == pb[i]) continue;
        try {
            const double x = coolplot::parse_number(pa[i]), y = coolplot::parse_number(pb[i]);
            if (std::fabs(x - y) > 1e-11 * std::fabs(y)) return false;
        } catch (const coolplot::Error&) {
            return false;
        }
    }
    return true;
}

bool sameIds(const std::vector<std::string>& a, const std::vector<std::string>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!sameId(a[i], b[i])) return false;
    return true;
}

/// The relative difference of two values in units of a span.
double spanDiff(double a, double b, double span) { return std::fabs(a - b) / span; }

/// The R134a states of a small heat pump cycle in CoolSolve's units (°C, Pa, J/kg, ...).
SolutionValues heatPumpSolution() {
    std::unique_ptr<CoolProp::AbstractState> st(CoolProp::AbstractState::factory("HEOS", "R134a"));
    SolutionValues values;
    const double p_ev = 2.0e5, p_cd = 1.0e6;
    st->update(CoolProp::PQ_INPUTS, p_ev, 1.0);  // 1: saturated vapour
    const double h1 = st->hmass(), s1 = st->smass();
    values["P_1"] = p_ev;
    values["h_1"] = h1;
    values["T_1"] = st->T() - 273.15;
    st->update(CoolProp::PSmass_INPUTS, p_cd, s1);  // 2: isentropic compression
    values["P_2"] = p_cd;
    values["h_2"] = st->hmass();
    st->update(CoolProp::PQ_INPUTS, p_cd, 0.0);  // 3: saturated liquid
    values["P_3"] = p_cd;
    values["h_3"] = st->hmass();
    values["x_3"] = 0.0;
    values["P_4"] = p_ev;  // 4: throttled
    values["h_4"] = st->hmass();
    for (int i = 1; i <= 4; ++i) {
        values["P[" + std::to_string(i) + "]"] = values["P_" + std::to_string(i)];
        values["h[" + std::to_string(i) + "]"] = values["h_" + std::to_string(i)];
    }
    return values;
}

json heatPumpCycle() {
    json states = json::array();
    for (int i = 1; i <= 4; ++i)
        states.push_back({{"number", i}, {"vars", {{"p", "P_" + std::to_string(i)}, {"h", "h_" + std::to_string(i)}}}});
    states[0]["label"] = "su,cp";
    states[0]["vars"]["T"] = "T_1";
    return {{"name", "Heat pump"}, {"states", states}};
}

json defaultSpec(const std::string& diagram) {
    const coolplot::DiagramType type = coolplot::DiagramType::parse(diagram);
    json isolines = json::array();
    for (const coolplot::Family& family : type.families())
        if (family.key != "u") isolines.push_back({{"key", family.key}, {"num", family.default_num}});
    return {{"diagram", diagram},
            {"units", "engineering"},
            {"limits", "DEF"},
            {"isolines", isolines},
            {"show_saturation", true},
            {"settings", {{"branding", {{"references", "none"}, {"versions", false}}}}}};
}

/// The index of the point of a curve closest to (x, y), relative distances.
size_t nearest(const Curve& curve, double x, double y) {
    size_t best = 0;
    double bestDistance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < curve.x.size(); ++i) {
        const double d = std::hypot((curve.x[i] - x) / x, (curve.y[i] - y) / y);
        if (d < bestDistance) {
            bestDistance = d;
            best = i;
        }
    }
    return best;
}

}  // namespace

// ---------------------------------------------------------------------------
// Planning (engine/planning.py)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram planning: round values as CoolPlot's planning", "[diagram]") {
    // Reference values from CoolPlot 7338f48 (Python, coolplot.engine.planning: nice_linear, nice_log,
    // _round_values), computed in a virtual environment with numpy 2.4
    CHECK(sameNumbers(niceLinear(-23.4, 87.6, 10), {-20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0}));
    CHECK(sameNumbers(niceLinear(0.21, 0.93, 7), {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}));
    CHECK(sameNumbers(niceLog(0.37, 45.0, 8), {0.5, 1.0, 2.0, 5.0, 10.0, 20.0}));
    CHECK(sameNumbers(niceLog(1.2, 2.9, 5), {1.2, 1.5, 2.0, 2.5}));
    const coolplot::UnitSystem& eng = coolplot::get_unit_system("engineering");
    CHECK(sameNumbers(roundValues(eng.at("T"), 250.0, 450.0, 10, false),
                      {253.14999999999998, 273.15, 293.15, 313.15, 333.15, 353.15, 373.15, 393.15, 413.15, 433.15}));
    CHECK(sameNumbers(roundValues(eng.at("p"), 1.5e4, 4.2e6, 9, true),
                      {20000.0, 49999.99999999999, 99999.99999999999, 199999.99999999997, 499999.99999999994, 999999.9999999999,
                       1999999.9999999998}));
    CHECK(isolineId("T", 253.14999999999998) == "iso/T/253.15");
    CHECK(isolineId("rho", 7.958114415791234) == "iso/rho/7.95811441579");
    CHECK_THROWS_AS(niceLinear(0.0, std::numeric_limits<double>::infinity(), 5), DiagramError);
}

TEST_CASE("Diagram planning: requests as CoolPlot's plan()", "[diagram]") {
    const Domain domain = referenceDomain();
    SECTION("families, rounding on linear and log axes, explicit values, single lines") {
        // CoolPlot 7338f48: DiagramSpec("R134a", "ph", units="engineering", settings=interactive points (60, 100)),
        // T num 10 rounded with labels, s num 5, rho num 8 rounded, Q values (0, 0.25, 0.5, 1, 1.5) with 40
        // points, lines Ambient (T 298.15), Relief (p 2e6, an axis) and Dry (Q 0.9); plan(spec, domain) gives:
        const json spec = {{"diagram", "ph"},
                           {"units", "engineering"},
                           {"settings", {{"calc", {{"points", 60}, {"saturation_points", 100}}}}},
                           {"isolines",
                            {{{"key", "T"}, {"num", 10}, {"rounding", true}, {"labels", true}},
                             {{"key", "s"}, {"num", 5}},
                             {{"key", "rho"}, {"num", 8}, {"rounding", true}},
                             {{"key", "Q"}, {"values_SI", {0.0, 0.25, 0.5, 1.0, 1.5}}, {"points", 40}}}},
                           {"lines",
                            {{{"name", "Ambient"}, {"key", "T"}, {"value_SI", 298.15}},
                             {{"name", "Relief"}, {"key", "p"}, {"value_SI", 2e6}},
                             {{"name", "Dry"}, {"key", "Q"}, {"value_SI", 0.9}}}},
                           {"show_saturation", true}};
        const std::vector<std::string> expected = {
            "iso/T/253.15 60 labels", "iso/T/273.15 60 labels", "iso/T/293.15 60 labels", "iso/T/313.15 60 labels",
            "iso/T/333.15 60 labels", "iso/T/353.15 60 labels", "iso/T/373.15 60 labels", "iso/T/393.15 60 labels",
            "iso/T/413.15 60 labels", "iso/T/433.15 60 labels", "iso/s/850 60",           "iso/s/1162.5 60",
            "iso/s/1475 60",          "iso/s/1787.5 60",        "iso/s/2100 60",          "iso/rho/1 60",
            "iso/rho/2 60",           "iso/rho/5 60",           "iso/rho/10 60",          "iso/rho/20 60",
            "iso/rho/50 60",          "iso/rho/100 60",         "iso/rho/200 60",         "iso/rho/500 60",
            "iso/rho/1000 60",        "iso/Q/0 40",             "iso/Q/0.25 40",          "iso/Q/0.5 40",
            "iso/Q/1 40",             "line/Ambient 60",        "line/Dry 100"};
        const std::vector<CurveRequest> requests = plan(PlanSpec::fromJson(spec), domain);
        CHECK(describe(requests) == expected);
        CHECK(coolplot::same_number(requests[0].value_SI, 253.15));
        CHECK_FALSE(requests[0].anchor_SI.has_value());
    }
    SECTION("an evenly spaced log family and a range; the saturation curve hidden") {
        // CoolPlot: DiagramSpec("R134a", "Ts", units="SI", isolines=(rho num 4, h num 3 range (2e5, 4e5)),
        // show_saturation=False) -> rho 0.6, 7.95811441579, 105.552641758, 1400 and h 2e5, 3e5, 4e5, 250 points
        const json spec = {{"diagram", "Ts"},
                           {"units", "SI"},
                           {"isolines", {{{"key", "rho"}, {"num", 4}}, {{"key", "h"}, {"num", 3}, {"range_SI", {2e5, 4e5}}}}},
                           {"show_saturation", false}};
        CHECK(describe(plan(PlanSpec::fromJson(spec), domain)) ==
              std::vector<std::string>{"iso/rho/0.6 250", "iso/rho/7.95811441579 250", "iso/rho/105.552641758 250",
                                       "iso/rho/1400 250", "iso/h/200000 250", "iso/h/300000 250", "iso/h/400000 250"});
    }
    SECTION("the saturation curve without its family, and with show_saturation unset") {
        // CoolPlot: T values (273.15,) with show_saturation=True -> the T line, then iso/Q/0 and iso/Q/1; a Q family
        // of 3 values with show_saturation=None -> iso/Q/0, 0.5, 1
        const json shown = {{"diagram", "ph"}, {"isolines", {{{"key", "T"}, {"values_SI", {273.15}}}}}, {"show_saturation", true}};
        CHECK(describe(plan(PlanSpec::fromJson(shown), domain)) ==
              std::vector<std::string>{"iso/T/273.15 250", "iso/Q/0 250", "iso/Q/1 250"});
        const json unset = {{"diagram", "ph"}, {"isolines", {{{"key", "Q"}, {"num", 3}}}}};
        CHECK(describe(plan(PlanSpec::fromJson(unset), domain)) ==
              std::vector<std::string>{"iso/Q/0 250", "iso/Q/0.5 250", "iso/Q/1 250"});
        const json none = {{"diagram", "ph"}, {"isolines", {{{"key", "T"}, {"values_SI", {273.15}}}}}};
        CHECK(describe(plan(PlanSpec::fromJson(none), domain)) == std::vector<std::string>{"iso/T/273.15 250"});
    }
    SECTION("spec values CoolPlot refuses") {
        CHECK_THROWS_WITH(PlanSpec::fromJson({{"diagram", "ph"}, {"isolines", {{{"key", "p"}}}}}),
                          ContainsSubstring("cannot be drawn in a ph diagram"));
        CHECK_THROWS_WITH(PlanSpec::fromJson({{"diagram", "ph"}, {"isolines", {{{"key", "T"}}, {{"key", "T"}}}}}),
                          ContainsSubstring("only once"));
        CHECK_THROWS_WITH(PlanSpec::fromJson({{"settings", {{"calc", {{"solver", "magic"}}}}}}), ContainsSubstring("Unknown solver"));
        CHECK_THROWS_WITH(PlanSpec::fromJson({{"diagram", "xy"}}), ContainsSubstring("Unknown diagram type"));
        CHECK_THROWS_WITH(PlanSpec::fromJson({{"limits", "WIDE"}}), ContainsSubstring("Unknown limits preset"));
        CHECK_THROWS_WITH(PlanSpec::fromJson({{"diagram", "pT"}, {"saturation_scale", {{"Q", 1.0}}}}),
                          ContainsSubstring("saturation scale"));
        CHECK_THROWS_WITH(PlanSpec::fromJson({{"isolines", {{{"key", "Q"}, {"anchor", {{"by", "T"}}}}}}}),
                          ContainsSubstring("vapour quality cannot be placed"));
    }
}

// ---------------------------------------------------------------------------
// Domain (thermo/limits.py, thermo/provider.py)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram domain: region presets and relative bounds", "[diagram]") {
    const TpLimits def = tpLimitsPreset("def");
    CHECK(def.T_min.kind == Bound::Kind::Relative);
    CHECK(def.T_min.value == 1.01);
    CHECK(def.p_max.value == 2.25);
    const TpLimits achp = tpLimitsPreset("ACHP");
    CHECK(achp.T_min.kind == Bound::Kind::Absolute);
    CHECK(achp.T_max.value == 493.15);
    CHECK(achp.p_max.kind == Bound::Kind::Relative);
    CHECK(tpLimitsPreset("NONE").T_min.kind == Bound::Kind::None);
    CHECK_THROWS_WITH(tpLimitsPreset("WIDE"), ContainsSubstring("expected one of ['ACHP', 'DEF', 'NONE', 'ORC']"));
    CHECK(TpLimits::fromJson(def.toJson()).toJson() == def.toJson());

    FluidModel water(chartFluid("Water"));
    const double margin = 1e-7;
    const auto T_sat = water.saturationRange("T", margin);
    const auto p_sat = water.saturationRange("p", margin);
    CHECK_THAT(T_sat.second, WithinRel(647.096 * (1 - margin), 1e-12));
    const std::array<double, 4> box = resolveTpBox(water, def, margin);
    CHECK(box[0] == 1.01 * T_sat.first);
    CHECK(box[1] == std::min(2.25 * T_sat.second, water.tMaxK()));
    CHECK(box[2] == 1.01 * p_sat.first);
    CHECK(box[3] == 2.25 * p_sat.second);
    // Bounds given as {"relative": f} and absolute numbers, clamped to the equation of state
    const TpLimits mixed = TpLimits::fromJson({{"T_min", 280.0}, {"T_max", {{"relative", 1.0}}}, {"p_min", nullptr}, {"p_max", 1e12}});
    const std::array<double, 4> clamped = resolveTpBox(water, mixed, margin);
    CHECK(clamped[0] == 280.0);
    CHECK(clamped[1] == T_sat.second);
    CHECK(clamped[2] == water.pMinPa());
    CHECK(clamped[3] == water.pMaxPa());
    CHECK_THROWS_WITH(resolveTpBox(water, TpLimits::fromJson({{"T_min", 500.0}, {"T_max", 400.0}}), margin),
                      ContainsSubstring("Empty plotting region"));
    CHECK_THROWS_AS(TpLimits::fromJson({{"T_min", "cold"}}), DiagramError);

    const Domain domain = chartDomain(water, def, margin);
    CHECK(domain.saturation);
    CHECK(domain.rangeSI("Q") == Range{0.0, 1.0});
    CHECK(domain.rangeSI("T") == Range{box[0], box[1]});
    CHECK(domain.rangeSI("v").first > 0.0);
    const auto points = markedPointStates(water, {"critical", "triple"}, margin);
    REQUIRE(points.size() == 3);
    CHECK(points[0].first == "critical");
    CHECK_THAT(points[0].second.at("p"), WithinRel(22.064e6, 1e-6));
    CHECK(points[1].first == "triple_liquid");
    CHECK(points[2].second.at("Q") == 1.0);
}

// ---------------------------------------------------------------------------
// The CoolPlot fixtures (tools/coolplot_fixtures.py)
// ---------------------------------------------------------------------------

namespace {

// Comparison tolerances with CoolPlot's Python numerics on the CoolProp 8.0.0 wheel (tests/fixtures/coolplot/
// header.json), set from the differences measured with the pinned CoolProp 8eb3b778 (docs/coolplot_integration.md
// §6.1): at most 7.8e-14 for the domain, 4.6e-8 for the curves (the fixtures keep 8 digits), 1.9e-10 for the
// labels, every curve with the same number of points. They compare two implementations; no tolerance of the glue
// itself depends on them.
constexpr double DOMAIN_TOL = 1e-12;  // ranges, box and saturation states, relative to the span of the quantity
constexpr double CURVE_TOL = 2e-7;    // difference of a point of a curve, in units of the axis spans
constexpr double LABEL_TOL = 1e-9;    // label positions and axis limits, in units of the axis spans

struct FixtureResult {
    double domain = 0.0;      // largest difference of the domain
    double curve = 0.0;       // largest distance between curves
    double label = 0.0;       // largest difference of a text position
    size_t curvesSameSize = 0, curvesCompared = 0;
};

FixtureResult compareFixture(const std::string& name) {
    INFO("fixture " << name);
    const json fixture = readJson(fixtureDir() / (name + ".json"));
    const json& input = fixture.at("input");
    const json& specRecord = input.at("spec");
    const PlanSpec spec = PlanSpec::fromJson(specRecord);
    const Chart chart(ChartFluid::fromRecord(specRecord.at("fluid")), spec);
    FixtureResult result;

    // 1. The plan on CoolPlot's own domain: the same requests, bit for bit
    Domain theirs;
    const json& d = input.at("domain");
    for (int i = 0; i < 4; ++i) theirs.tp_box[i] = d.at("tp_box").at(i).get<double>();
    for (auto it = d.at("ranges_SI").begin(); it != d.at("ranges_SI").end(); ++it)
        theirs.ranges_SI.push_back({it.key(), {it.value().at(0).get<double>(), it.value().at(1).get<double>()}});
    theirs.saturation = true;
    for (const IsolineSet& iso : spec.isolines) {
        if (!iso.anchor) continue;
        std::vector<Range> pairs;
        for (const json& r : fixture.at("requests"))
            if (r.at("key") == iso.key && !r.at("anchor_SI").is_null())
                pairs.push_back({r.at("anchor_SI").get<double>(), r.at("value_SI").get<double>()});
        theirs.anchored.push_back({iso.key, pairs});
    }
    const std::vector<CurveRequest> planned = plan(spec, theirs);
    REQUIRE(planned.size() == fixture.at("requests").size());
    for (size_t i = 0; i < planned.size(); ++i) {
        const json& r = fixture.at("requests")[i];
        INFO("request " << r.dump());
        CHECK(planned[i].item_id == r.at("item_id"));
        CHECK(planned[i].key == r.at("key"));
        CHECK(coolplot::same_number(planned[i].value_SI, r.at("value_SI").get<double>()));
        CHECK(planned[i].points == r.at("points"));
        CHECK(planned[i].labels == r.at("labels"));
        CHECK(planned[i].anchor_SI.has_value() == !r.at("anchor_SI").is_null());
    }

    // 2. The domain computed with the pinned CoolProp
    std::vector<std::string> warnings;
    const Domain ours = Engine::shared().domain(chart, warnings);
    auto span = [&theirs](const std::string& key) {
        const Range r = theirs.rangeSI(key);
        return coolplot::get_quantity(key).log_scale ? std::fabs(r.second) : r.second - r.first;
    };
    for (int i = 0; i < 4; ++i)
        result.domain = std::max(result.domain, std::fabs(ours.tp_box[i] - theirs.tp_box[i]) / std::fabs(theirs.tp_box[i]));
    for (const auto& [key, range] : theirs.ranges_SI) {
        const Range mine = ours.rangeSI(key);
        result.domain = std::max({result.domain, spanDiff(mine.first, range.first, span(key)), spanDiff(mine.second, range.second, span(key))});
    }
    if (d.contains("points")) {
        REQUIRE(ours.points.size() == d.at("points").size());
        for (const auto& [point, state] : ours.points)
            for (const auto& [key, value] : state)
                result.domain = std::max(result.domain, spanDiff(value, d.at("points").at(point).at(key).get<double>(), span(key)));
    }
    if (d.contains("scale")) {
        REQUIRE(ours.scale.size() == d.at("scale").size());
        const std::string by = scaleKey(spec.diagram);
        for (size_t i = 0; i < ours.scale.size(); ++i) {
            result.domain = std::max(result.domain, spanDiff(ours.scale[i].first, d.at("scale")[i].at(0).get<double>(), span(by)));
            result.domain = std::max(result.domain, spanDiff(ours.scale[i].second, d.at("scale")[i].at(1).get<double>(), span(spec.diagram.y)));
        }
    }
    if (d.contains("dome")) {
        REQUIRE(ours.dome.size() == d.at("dome").size());
        for (size_t i = 0; i < ours.dome.size(); ++i)
            for (const auto& [key, value] : ours.dome[i].dew)
                result.domain = std::max(result.domain, spanDiff(value, d.at("dome")[i].at(2).at(key).get<double>(), span(key)));
    }

    // 3. The curves: every request of CoolPlot's plan on our domain, the same lines drawn, the same shape. The
    //    values of a family spaced evenly over a range CoolProp computes may differ in their 12th digit (sameId)
    json diagnostics = json::array();
    const std::vector<CurveRequest> requests = plan(spec, ours);
    const std::vector<PlannedCurve> curves = Engine::shared().curves(chart, ours, requests, diagnostics, warnings);
    std::vector<std::string> ourIds, theirIds;
    for (const json& r : fixture.at("requests")) theirIds.push_back(r.at("item_id"));
    for (const CurveRequest& r : requests) ourIds.push_back(r.item_id);
    CHECK(sameIds(ourIds, theirIds));
    ourIds.clear();
    theirIds.clear();
    for (const auto& [request, _] : curves) ourIds.push_back(request.item_id);
    for (const json& c : input.at("curves")) theirIds.push_back(c.at("request").at("item_id"));
    CHECK(sameIds(ourIds, theirIds));
    const bool xLog = spec.diagram.x_log(), yLog = spec.diagram.y_log();
    auto axisSpan = [&](const std::string& key, bool log) {
        const Range r = theirs.rangeSI(key);
        return log ? std::log10(r.second / r.first) : r.second - r.first;
    };
    const double xSpan = axisSpan(spec.diagram.x, xLog), ySpan = axisSpan(spec.diagram.y, yLog);
    for (const json& c : input.at("curves")) {
        const std::string id = c.at("request").at("item_id");
        const auto mine = std::find_if(curves.begin(), curves.end(), [&id](const PlannedCurve& p) { return sameId(p.first.item_id, id); });
        if (mine == curves.end()) continue;
        // The fixture keeps every `stride`-th point and the last one (tools/coolplot_fixtures.py)
        const std::vector<double> tx = numbers(c.at("curve").at("x")), ty = numbers(c.at("curve").at("y"));
        const size_t stride = c.at("curve").at("stride"), n = c.at("curve").at("n_points");
        const Curve& curve = *mine->second;
        double distance = 0.0;
        if (curve.x.size() == n) {
            // The same points: compared one by one, gaps (NaN) included
            for (size_t k = 0; k < tx.size(); ++k) {
                const size_t i = std::min(k * stride, n - 1);
                const bool theirsValid = std::isfinite(tx[k]) && std::isfinite(ty[k]);
                const bool oursValid = std::isfinite(curve.x[i]) && std::isfinite(curve.y[i]);
                if (theirsValid != oursValid) distance = std::numeric_limits<double>::infinity();
                if (!theirsValid || !oursValid) continue;
                distance = std::max({distance, std::fabs(axisValue(curve.x[i], xLog) - axisValue(tx[k], xLog)) / xSpan,
                                     std::fabs(axisValue(curve.y[i], yLog) - axisValue(ty[k], yLog)) / ySpan});
            }
        } else {
            distance = polylineDistance(tx, ty, curve.x, curve.y, xLog, yLog, xSpan, ySpan);
        }
        if (distance > result.curve) {
            result.curve = distance;
            UNSCOPED_INFO(name << ": largest curve distance " << distance << " on " << id);
        }
        ++result.curvesCompared;
        if (curve.x.size() == n) ++result.curvesSameSize;
    }

    // 4. The scene of our compose record: CoolPlot's items, groups, roles, z order, styles and texts
    json record = composeRecord(specRecord, spec, ours, curves, {}, std::nullopt, json());
    const json scene = coolplot::compose(record).to_json(true);
    const json& expected = fixture.at("scene");
    CHECK(scene.at("axes").at("x_label") == expected.at("axes").at("x_label"));
    CHECK(scene.at("axes").at("y_label") == expected.at("axes").at("y_label"));
    const bool sxLog = expected.at("axes").at("x_log"), syLog = expected.at("axes").at("y_log");
    const json& xl = expected.at("axes").at("x_limits");
    const json& yl = expected.at("axes").at("y_limits");
    const double sxSpan = axisValue(xl[1], sxLog) - axisValue(xl[0], sxLog), sySpan = axisValue(yl[1], syLog) - axisValue(yl[0], syLog);
    for (int i = 0; i < 2; ++i) {
        result.label = std::max(result.label, std::fabs(axisValue(scene["axes"]["x_limits"][i], sxLog) - axisValue(xl[i], sxLog)) / sxSpan);
        result.label = std::max(result.label, std::fabs(axisValue(scene["axes"]["y_limits"][i], syLog) - axisValue(yl[i], syLog)) / sySpan);
    }
    std::vector<std::string> ids, expectedIds;
    for (const json& item : scene.at("items")) ids.push_back(item.at("id"));
    for (const json& item : expected.at("items")) expectedIds.push_back(item.at("id"));
    CHECK(sameIds(ids, expectedIds));  // drawing order included
    for (size_t i = 0; i < expected.at("items").size() && i < scene.at("items").size(); ++i) {
        const json& item = expected.at("items")[i];
        const json* mine = &scene.at("items")[i];
        INFO("item " << item.at("id"));
        for (const char* field : {"type", "group", "role", "z_order", "legend", "visible", "text", "clip"})
            if (item.contains(field)) CHECK(mine->at(field) == item.at(field));
        CHECK(mine->at("style") == expected.at("styles").at(item.at("style").get<size_t>()));  // styles written once each
        if (item.at("type") == "text") {
            const double dx = std::fabs(axisValue(mine->at("x"), sxLog) - axisValue(item.at("x"), sxLog)) / sxSpan;
            const double dy = std::fabs(axisValue(mine->at("y"), syLog) - axisValue(item.at("y"), syLog)) / sySpan;
            result.label = std::max({result.label, dx, dy});
        }
    }
    return result;
}

}  // namespace

TEST_CASE("Diagram fixtures: the compose records and scenes of CoolPlot", "[diagram]") {
    const json header = readJson(fixtureDir() / "header.json");
    CHECK(header.at("coolplot_commit") == COOLSOLVE_COOLPLOT_COMMIT);
    for (const char* name : {"r134a_ph_default", "water_ts", "r134a_th", "r290_ph_achp", "co2_ph_marked"}) {
        const FixtureResult r = compareFixture(name);
        INFO(name << ": domain " << r.domain << ", curves " << r.curve << " (" << r.curvesSameSize << " of " << r.curvesCompared
                  << " with the same points), labels " << r.label);
        CHECK(r.curvesSameSize == r.curvesCompared);
        CHECK(r.domain < DOMAIN_TOL);
        CHECK(r.curve < CURVE_TOL);
        CHECK(r.label < LABEL_TOL);
    }
}

// ---------------------------------------------------------------------------
// Charts end to end (POST /api/v1/diagram/scene)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram scene: the default p-h chart of R134a and a T-s chart of water", "[diagram]") {
    SECTION("R134a p-h") {
        SceneRequest request;
        request.fluid = "R134a";
        request.spec = defaultSpec("ph");
        request.spec["marked_points"] = {"critical"};
        request.figureSizePt = std::make_pair(600.0, 420.0);
        const json response = scene(request);
        const std::string svg = response.at("svg");
        CHECK(svg.rfind("<svg", 0) == 0);
        for (const char* id : {"iso/Q/0", "iso/Q/1", "iso/Q/0.5", "point/critical"})
            CHECK_THAT(svg, ContainsSubstring(std::string("data-item=\"") + id + "\""));
        // No point is drawn where the flash failed: no NaN or infinity anywhere ("dominant-baseline" has "nan")
        CHECK_FALSE(std::regex_search(svg, std::regex("[^A-Za-z](nan|NaN|inf|Infinity)[^A-Za-z]")));
        CHECK(response.at("diagnostics").at("warnings").empty());
        const json& axes = response.at("axes");
        CHECK(axes.at("x_key") == "h");
        CHECK(axes.at("y_key") == "p");
        CHECK(axes.at("y_log") == true);
        CHECK(axes.at("x_unit").at("label") == "kJ/kg");
        CHECK(axes.at("x_unit").at("scale") == 1e-3);
        CHECK(axes.at("y_unit").at("scale") == 1e-5);
        CHECK(axes.at("x_limits").size() == 2);
        CHECK(response.at("spec").at("fluid").at("components") == json::array({"R134a"}));
        CHECK(response.at("spec").at("cycles").empty());
        CHECK(response.at("states").empty());
        for (const char* key : {"curves", "compose", "svg"}) CHECK(response.at("timing_ms").at(key).get<double>() >= 0.0);
        const json& curves = response.at("diagnostics").at("curves");
        CHECK(curves.size() == 15 + 15 + 15 + 11);  // T, s, rho, Q
        // A repeated request is served from the caches and gives the same document
        CHECK(scene(request).at("svg") == svg);
    }
    SECTION("water T-s with labels, the critical and triple points, the model's units") {
        SceneRequest request;
        request.fluid = "Steam";  // a CoolSolve name of Water
        request.spec = defaultSpec("Ts");
        request.spec["isolines"][0]["labels"] = true;
        request.spec["marked_points"] = {"triple", "critical"};
        request.spec["limits"] = "NONE";  // the triple point is drawn when the region reaches down to it
        request.spec["units"] = coolsolveUnitSystem();
        const json response = scene(request);
        const std::string svg = response.at("svg");
        for (const char* id : {"iso/Q/0", "iso/Q/1", "point/critical", "point/triple"})
            CHECK_THAT(svg, ContainsSubstring(std::string("data-item=\"") + id + "\""));
        CHECK(response.at("spec").at("fluid").at("components") == json::array({"Water"}));
        CHECK_FALSE(response.at("labels").empty());
        for (const json& label : response.at("labels")) CHECK(label.at("role") == "isoline_label");
        CHECK(response.at("axes").at("y_unit").at("offset") == -273.15);
        CHECK(response.at("axes").at("x_unit").at("scale") == 1.0);
        // The labels sent back keep their places
        SceneRequest again = request;
        again.previousLabels = response.at("labels");
        CHECK(scene(again).at("labels") == response.at("labels"));
    }
}

TEST_CASE("Diagram invariants: the dome, isolines through it, solved states on their lines", "[diagram]") {
    auto curveOf = [](const std::vector<PlannedCurve>& curves, const std::string& id) -> const Curve& {
        for (const auto& [request, curve] : curves)
            if (request.item_id == id) return *curve;
        FAIL("no curve " << id);
        throw std::logic_error("unreachable");
    };
    SECTION("R134a p-h: the dome closes at the critical point, an isotherm is horizontal inside it") {
        json spec = defaultSpec("ph");
        spec["isolines"] = {{{"key", "T"}, {"values_SI", {273.15, 313.15}}}, {{"key", "Q"}, {"num", 3}}};
        const Chart chart(chartFluid("R134a"), PlanSpec::fromJson(spec));
        std::vector<std::string> warnings;
        json diagnostics;
        const Domain domain = Engine::shared().domain(chart, warnings);
        const std::vector<PlannedCurve> curves = Engine::shared().curves(chart, domain, plan(chart.spec, domain), diagnostics, warnings);
        FluidModel fluid(chartFluid("R134a"));
        const State& crit = fluid.critical();
        const Range h = domain.rangeSI("h");
        // The bubble and dew lines end next to the critical point and meet there
        const Curve& bubble = curveOf(curves, "iso/Q/0");
        const Curve& dew = curveOf(curves, "iso/Q/1");
        auto top = [](const Curve& c) {
            size_t best = 0;
            for (size_t i = 0; i < c.y.size(); ++i)
                if (std::isfinite(c.y[i]) && !(c.y[i] <= c.y[best])) best = i;
            return best;
        };
        const size_t ib = top(bubble), id = top(dew);
        CHECK_THAT(bubble.y[ib], WithinRel(crit.at("p"), 1e-5));
        CHECK_THAT(dew.y[id], WithinRel(crit.at("p"), 1e-5));
        CHECK(std::fabs(bubble.x[ib] - crit.at("h")) < 0.01 * (h.second - h.first));
        CHECK(std::fabs(dew.x[id] - crit.at("h")) < 0.01 * (h.second - h.first));
        // An isotherm holds its saturation states and is horizontal between them
        for (double T : {273.15, 313.15}) {
            const State liquid = fluid.statePoint("T", T, "Q", 0.0), vapour = fluid.statePoint("T", T, "Q", 1.0);
            const Curve& isotherm = curveOf(curves, isolineId("T", T));
            const size_t iL = nearest(isotherm, liquid.at("h"), liquid.at("p")), iV = nearest(isotherm, vapour.at("h"), vapour.at("p"));
            CHECK_THAT(isotherm.x[iL], WithinRel(liquid.at("h"), 1e-9));
            CHECK_THAT(isotherm.x[iV], WithinRel(vapour.at("h"), 1e-9));
            CHECK(std::max(iL, iV) - std::min(iL, iV) >= 3);
            for (size_t i = std::min(iL, iV); i <= std::max(iL, iV); ++i) CHECK_THAT(isotherm.y[i], WithinRel(liquid.at("p"), 1e-9));
        }
    }
    SECTION("water T-s: an isobar is horizontal inside the dome") {
        json spec = defaultSpec("Ts");
        spec["isolines"] = {{{"key", "p"}, {"values_SI", {1e5}}}};
        const Chart chart(chartFluid("Water"), PlanSpec::fromJson(spec));
        std::vector<std::string> warnings;
        json diagnostics;
        const Domain domain = Engine::shared().domain(chart, warnings);
        const std::vector<PlannedCurve> curves = Engine::shared().curves(chart, domain, plan(chart.spec, domain), diagnostics, warnings);
        FluidModel fluid(chartFluid("Water"));
        const State liquid = fluid.statePoint("p", 1e5, "Q", 0.0), vapour = fluid.statePoint("p", 1e5, "Q", 1.0);
        const Curve& isobar = curveOf(curves, "iso/p/100000");
        const size_t iL = nearest(isobar, liquid.at("s"), liquid.at("T")), iV = nearest(isobar, vapour.at("s"), vapour.at("T"));
        CHECK_THAT(isobar.x[iL], WithinRel(liquid.at("s"), 1e-9));
        CHECK_THAT(isobar.x[iV], WithinRel(vapour.at("s"), 1e-9));
        CHECK(std::max(iL, iV) - std::min(iL, iV) >= 3);
        for (size_t i = std::min(iL, iV); i <= std::max(iL, iV); ++i) CHECK_THAT(isobar.y[i], WithinRel(liquid.at("T"), 1e-9));
    }
    SECTION("a solved state lies on the lines of its own properties") {
        const SolutionValues values = heatPumpSolution();
        FluidModel fluid(chartFluid("R134a"));
        const State state2 = fluid.statePoint("p", values.at("P_2"), "h", values.at("h_2"));
        json spec = defaultSpec("ph");
        spec["isolines"] = json::array();
        spec["lines"] = {{{"name", "T2"}, {"key", "T"}, {"value_SI", state2.at("T")}},
                         {{"name", "s2"}, {"key", "s"}, {"value_SI", state2.at("s")}}};
        const Chart chart(chartFluid("R134a"), PlanSpec::fromJson(spec));
        std::vector<std::string> warnings;
        json diagnostics;
        const Domain domain = Engine::shared().domain(chart, warnings);
        const std::vector<PlannedCurve> curves = Engine::shared().curves(chart, domain, plan(chart.spec, domain), diagnostics, warnings);
        const Range h = domain.rangeSI("h"), p = domain.rangeSI("p");
        const double xSpan = h.second - h.first, ySpan = std::log10(p.second / p.first);
        for (const char* id : {"line/T2", "line/s2"}) {
            const Curve& line = curveOf(curves, id);
            INFO(id);
            CHECK(polylineDistance({values.at("h_2")}, {values.at("P_2")}, line.x, line.y, false, true, xSpan, ySpan) < 1e-4);
        }
    }
}

// ---------------------------------------------------------------------------
// Cycles (phase 1: straight chains)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram cycles: states from solved variables, straight chains", "[diagram]") {
    const SolutionValues values = heatPumpSolution();
    const coolplot::DiagramType ph = coolplot::DiagramType::parse("ph");
    FluidModel fluid(chartFluid("R134a"));
    std::vector<std::string> warnings;

    SECTION("explicit states: full states, the solved values kept, one closed loop") {
        const std::vector<ResolvedCycle> cycles = resolveCycles(json::array({heatPumpCycle()}), &values, fluid, ph, warnings);
        CHECK(warnings.empty());
        REQUIRE(cycles.size() == 1);
        const ResolvedCycle& cycle = cycles[0];
        REQUIRE(cycle.states.size() == 4);
        CHECK(cycle.states.at(1).at("h") == values.at("h_1"));
        CHECK(cycle.states.at(1).at("T") == values.at("T_1") + 273.15);  // °C to K
        CHECK(cycle.states.at(2).at("p") == values.at("P_2"));
        CHECK(cycle.states.at(2).at("Q") == -1.0);  // single phase, as CoolProp (and CoolPlot) report it
        CHECK_THAT(cycle.states.at(4).at("Q"), Catch::Matchers::WithinAbs(0.25, 0.2));
        for (const char* key : {"T", "p", "h", "s", "rho", "v", "u", "Q"}) CHECK(cycle.states.at(3).count(key) == 1);
        REQUIRE(cycle.chains.size() == 1);
        REQUIRE(cycle.chains[0].size() == 5);  // straight: the states themselves, 1 2 3 4 1
        CHECK(cycle.chains[0][0] == cycle.states.at(1));
        CHECK(cycle.chains[0][2] == cycle.states.at(3));
        CHECK(cycle.chains[0][4] == cycle.states.at(1));
        CHECK(cycle.spec.at("loops") == json::parse("[[1,2,3,4,1]]"));
        CHECK(cycle.spec.at("paths") == json::parse(R"([[[1,2],"straight"],[[2,3],"straight"],[[3,4],"straight"],[[4,1],"straight"]])"));
        CHECK(cycle.spec.at("labels") == json::parse(R"([[1,"su,cp"]])"));
        CHECK(cycle.spec.at("states")[0][1].at("state").at("h") == values.at("h_1"));
        CHECK(cycle.vars.at(1).at("T") == "T_1");
    }
    SECTION("arrays, given chains, markers alone, paths auto") {
        json cycles = json::array();
        cycles.push_back({{"name", "Arrays"}, {"arrays", {{"P", "P"}, {"h", "h"}}}, {"loops", {{1, 2, 3}, {3, 4}}}});
        cycles.push_back({{"name", "Points"}, {"arrays", {{"p", "p"}, {"h", "H"}}}, {"loops", json::array()}, {"paths", "auto"}});
        cycles.push_back({{"name", "Empty"}, {"states", json::array()}});
        const std::vector<ResolvedCycle> resolved = resolveCycles(cycles, &values, fluid, ph, warnings);
        REQUIRE(resolved.size() == 2);
        CHECK(resolved[0].vars.at(3).at("p") == "P[3]");
        CHECK(resolved[0].chains.size() == 2);
        CHECK(resolved[0].chains[1].size() == 2);
        CHECK(resolved[1].index == 1);
        CHECK(resolved[1].states.size() == 4);
        CHECK(resolved[1].chains.empty());  // markers only
        CHECK(resolved[1].spec.at("loops") == json::parse("[[1,1],[2,2],[3,3],[4,4]]"));
        REQUIRE(warnings.size() == 1);
        CHECK_THAT(warnings[0], ContainsSubstring("phase 2"));
    }
    SECTION("a missing variable, no solve: the cycle is left out with a warning") {
        json cycle = heatPumpCycle();
        cycle["states"][2]["vars"]["p"] = "P_9";
        CHECK(resolveCycles(json::array({cycle}), &values, fluid, ph, warnings).empty());
        REQUIRE(warnings.size() == 1);
        CHECK_THAT(warnings[0], ContainsSubstring("'P_9' has no value"));
        warnings.clear();
        CHECK(resolveCycles(json::array({heatPumpCycle()}), nullptr, fluid, ph, warnings).empty());
        REQUIRE(warnings.size() == 1);
        CHECK_THAT(warnings[0], ContainsSubstring("no solved model"));
        warnings.clear();
        json one = heatPumpCycle();
        one["states"][1]["vars"].erase("h");
        CHECK(resolveCycles(json::array({one}), &values, fluid, ph, warnings).empty());
        CHECK_THAT(warnings.at(0), ContainsSubstring("needs two properties"));
    }
    SECTION("markers alone (loops []): a state that cannot be computed is left out, the others are drawn") {
        json markers = heatPumpCycle();
        markers["loops"] = json::array();
        markers["states"][1]["vars"]["p"] = "P_9";   // no value
        markers["states"][2]["vars"].erase("h");     // one property
        const std::vector<ResolvedCycle> resolved = resolveCycles(json::array({markers}), &values, fluid, ph, warnings);
        REQUIRE(resolved.size() == 1);
        CHECK(resolved[0].states.size() == 2);
        CHECK(resolved[0].states.count(1) == 1);
        CHECK(resolved[0].states.count(4) == 1);
        CHECK(resolved[0].spec.at("loops") == json::parse("[[1,1],[4,4]]"));
        REQUIRE(warnings.size() == 2);
        CHECK_THAT(warnings[0], ContainsSubstring("left out"));
        warnings.clear();
        for (json& state : markers["states"]) state["vars"] = {{"p", "P_9"}};
        CHECK(resolveCycles(json::array({markers}), &values, fluid, ph, warnings).empty());
        CHECK_THAT(warnings.back(), ContainsSubstring("none of its states"));
    }
    SECTION("definitions CoolPlot refuses") {
        CHECK_THROWS_WITH(resolveCycles(json::array({heatPumpCycle(), heatPumpCycle()}), &values, fluid, ph, warnings),
                          ContainsSubstring("unique"));
        CHECK_THROWS_AS(resolveCycles(json::object(), &values, fluid, ph, warnings), DiagramError);
    }
    SECTION("in a scene: the states, the spec CoolPlot reads, the cycle drawn") {
        SceneRequest request;
        request.fluid = "R134a";
        request.spec = defaultSpec("ph");
        request.cycles = json::array({heatPumpCycle()});
        request.solution = &values;
        const json response = scene(request);
        const json& states = response.at("states");
        REQUIRE(states.size() == 4);
        CHECK(states[0].at("cycle") == 0);
        CHECK(states[0].at("number") == 1);
        CHECK(states[0].at("label") == "su,cp");
        CHECK(states[0].at("vars").at("p") == "P_1");
        CHECK(states[1].at("values_SI").at("h") == values.at("h_2"));
        CHECK(response.at("spec").at("cycles")[0].at("name") == "Heat pump");
        const std::string svg = response.at("svg");
        CHECK_THAT(svg, ContainsSubstring("data-item=\"cycle/Heat pump/loop/0\""));
        CHECK_THAT(svg, ContainsSubstring("data-item=\"cycle/Heat pump/states\""));
        request.solution = nullptr;
        const json unsolved = scene(request);
        CHECK(unsolved.at("states").empty());
        CHECK_THAT(unsolved.at("diagnostics").at("warnings")[0].get<std::string>(), ContainsSubstring("is not drawn"));
    }
}

// ---------------------------------------------------------------------------
// Errors, options, definitions
// ---------------------------------------------------------------------------

TEST_CASE("Diagram requests that cannot be drawn", "[diagram]") {
    SceneRequest request;
    request.spec = defaultSpec("ph");
    request.fluid = "Unobtainium";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("Unknown fluid 'Unobtainium'"));
    REQUIRE(FluidRegistry::getFluid("N2")->getType() == FluidType::IdealGas);
    request.fluid = "N2";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("not drawable yet"));
    request.fluid = "AirH2O";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("not drawable yet"));
    request.fluid = "NH3H2O";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("not drawable yet"));
    request.fluid = "R134a";
    request.spec["diagram"] = "xy";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("Unknown diagram type 'xy'"));
    request.spec["diagram"] = "psychrometric";
    request.spec["isolines"] = json::array();
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("not drawable yet"));
    request.spec = defaultSpec("ph");
    request.spec["settings"]["look"] = "neon";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("Unknown theme 'neon'"));
    request.spec = defaultSpec("ph");
    request.spec["units"] = "furlongs";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("Unknown unit system"));
    request.spec = defaultSpec("ph");
    request.figureSizePt = std::make_pair(10.0, 400.0);
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("figure size"));
    request.figureSizePt.reset();
    request.spec = json::array();
    CHECK_THROWS_AS(scene(request), DiagramError);
}

TEST_CASE("Diagram options and definitions", "[diagram]") {
    const json o = options();
    std::vector<std::string> names;
    for (const json& type : o.at("diagram_types")) names.push_back(type.at("name"));
    CHECK(names == std::vector<std::string>{"ph", "Ts", "hs", "ps", "prho", "Trho", "pT", "Th"});
    const json& th = o.at("diagram_types")[7];
    CHECK(th.at("x") == "h");
    CHECK(th.at("y") == "T");
    CHECK(th.at("families").size() == 5);  // p, s, rho, u, Q
    CHECK(o.at("diagram_types")[0].at("families")[0].at("default_num") == 15);
    CHECK(o.at("coolplot").at("core_version") == "1.5.0");
    CHECK(o.at("coolplot").at("commit") == COOLSOLVE_COOLPLOT_COMMIT);
    CHECK(o.at("themes") == json::array({"default", "dark", "print", "coolpack"}));
    CHECK(o.at("label_modes") == json::array({"rails", "edge", "margin", "inline"}));
    CHECK(o.at("regions").size() == 4);
    const json& custom = o.at("unit_systems").back();
    CHECK(custom.at("name") == "coolsolve");
    const coolplot::UnitSystem units = coolplot::UnitSystem::from_json(custom.at("record"));
    for (const std::string& key : stateKeys()) CHECK(units.units.count(key) == 1);
    CHECK(units.at("T").offset == -273.15);
    CHECK(units.at("T").difference_label == "K");
    CHECK(units.at("p").scale == 1.0);

    CHECK_NOTHROW(checkDefinitions(emptyDefinitions()));
    CHECK_NOTHROW(checkDefinitions(json::parse(R"({"diagrams": [{"name": "a", "fluid": "R134a", "cycles": [{"states": []}]}]})")));
    CHECK_THROWS_AS(checkDefinitions(json::array()), DiagramError);
    CHECK_THROWS_AS(checkDefinitions({{"format", "other"}}), DiagramError);
    CHECK_THROWS_AS(checkDefinitions({{"version", "1"}}), DiagramError);
    CHECK_THROWS_AS(checkDefinitions({{"diagrams", json::object()}}), DiagramError);
}

TEST_CASE("Diagram performance: the default charts of R134a and water, cold and warm", "[diagram]") {
    // Budget set from the Release build of 2026-11-12 (Linux, GCC): a cold scene (curves computed) took 28 ms for
    // R134a p-h and 70 ms for water T-s, a warm one (compose and SVG only) 3 ms. The budget is loose, 15 to 30 times
    // that, so that only a real regression trips it on a slower machine (docs/coolplot_integration.md §6.3)
    for (const auto& [fluid, diagram] : {std::pair<const char*, const char*>{"R134a", "ph"}, {"Water", "Ts"}}) {
        SceneRequest request;
        request.fluid = fluid;
        request.spec = defaultSpec(diagram);
        Engine::shared().clear();
        auto start = std::chrono::steady_clock::now();
        const json cold = scene(request);
        const double coldMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
        const json warm = scene(request);
        const double warmMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        INFO(fluid << " " << diagram << ": cold " << coldMs << " ms, warm " << warmMs << " ms, timing " << cold.at("timing_ms").dump());
        CHECK(coldMs < 1000.0);
        CHECK(warmMs < 100.0);
        CHECK(warm.at("svg") == cold.at("svg"));
    }
}
