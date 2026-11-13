// Tests of the Diagrams module (src/diagram/, docs/coolplot_integration.md): planning against CoolPlot's Python
// planning, the domain, the CoolPlot fixtures of tools/coolplot_fixtures.py (tests/fixtures/coolplot/), charts
// end to end, cycles from solved values, the invariants of §6.2 and the errors of a chart request.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "coolsolve/diagram/engine.h"
#include "coolsolve/diagram/humid_air.h"
#include "coolsolve/diagram/ideal_gas.h"
#include "coolsolve/diagram/states.h"
#include "coolsolve/fluids.h"
#include "coolsolve/parser.h"
#include "coolsolve/runner.h"
#include "coolsolve/variable_inference.h"

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include "HumidAirProp.h"

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
constexpr double CYCLE_TOL = 5e-9;    // the states along the paths of a cycle (measured 4.3e-10), relative to the span

struct FixtureResult {
    double domain = 0.0;      // largest difference of the domain
    double curve = 0.0;       // largest distance between curves
    double label = 0.0;       // largest difference of a text position
    double cycle = 0.0;       // largest difference of a state along the paths of a cycle
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
    // The span a difference of a quantity is measured in: its range, the temperature range for the temperatures of
    // humid air that have none (Tdp, and Thx on the psychrometric chart), the value of the one pressure of humid air
    auto span = [&theirs](const std::string& key) {
        const bool ranged = std::any_of(theirs.ranges_SI.begin(), theirs.ranges_SI.end(), [&key](const auto& r) { return r.first == key; });
        const std::string k = ranged ? key : (key == "Tdp" || key == "Thx" ? "T" : key);
        const Range r = theirs.rangeSI(k);
        return coolplot::get_quantity(k).log_scale || !(r.second > r.first) ? std::fabs(r.second) : r.second - r.first;
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

    // 4. The cycles: CoolPlot's states given by two properties become the solved values of a model, resolved as the
    //    Diagrams tab resolves them; the paths between the states must be CoolPlot's
    json definitions = json::array();
    SolutionValues solved;
    for (const json& c : specRecord.value("cycles", json::array())) {
        json states = json::array();
        for (const json& entry : c.at("states")) {
            const int number = entry.at(0);
            json vars = json::object();
            for (const json& prop : entry.at(1).at("props")) {
                const std::string key = prop.at(0), var = key + "_" + std::to_string(number);
                vars[key] = var;
                const bool temperature = key == "T" || key == "Twb" || key == "Tdp";
                solved[var] = temperature ? prop.at(1).get<double>() - 273.15 : prop.at(1).get<double>();  // the model's units
            }
            states.push_back({{"number", number}, {"vars", vars}});
        }
        for (const json& label : c.value("labels", json::array()))
            for (json& state : states)
                if (state.at("number") == label.at(0)) state["label"] = label.at(1);
        definitions.push_back({{"name", c.at("name")}, {"states", states}, {"loops", c.at("loops")}, {"paths", c.at("paths")},
                               {"steps", c.at("steps")}, {"numbers", c.at("numbers")}});
    }
    const std::vector<ResolvedCycle> cycles = Engine::shared().withFluid(chart.fluid, [&](PropertyProvider& fluid) {
        return resolveCycles(definitions, &solved, fluid, spec.diagram, warnings);
    });
    REQUIRE(cycles.size() == input.at("cycles").size());
    for (size_t i = 0; i < cycles.size(); ++i) {
        const json& chains = input.at("cycles")[i].at("chains");
        REQUIRE(cycles[i].chains.size() == chains.size());
        for (size_t c = 0; c < chains.size(); ++c) {
            REQUIRE(cycles[i].chains[c].size() == chains[c].size());
            for (size_t k = 0; k < chains[c].size(); ++k)
                for (auto it = chains[c][k].begin(); it != chains[c][k].end(); ++it) {
                    const double theirsValue = it.value(), oursValue = cycles[i].chains[c][k].at(it.key());
                    result.cycle = std::max(result.cycle, spanDiff(oursValue, theirsValue, it.key() == "Q" ? 1.0 : span(it.key())));
                }
        }
    }

    // 5. The scene of our compose record: CoolPlot's items, groups, roles, z order, styles and texts
    json record = composeRecord(specRecord, spec, ours, curves, cycles, std::nullopt, json());
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
    for (const char* name : {"r134a_ph_default", "water_ts", "r134a_th", "r290_ph_achp", "r134a_ph_cycle", "co2_ph_marked",
                             "psychrometric_ac", "mollier_heating"}) {
        const FixtureResult r = compareFixture(name);
        INFO(name << ": domain " << r.domain << ", curves " << r.curve << " (" << r.curvesSameSize << " of " << r.curvesCompared
                  << " with the same points), labels " << r.label << ", cycle " << r.cycle);
        CHECK(r.curvesSameSize == r.curvesCompared);
        CHECK(r.domain < DOMAIN_TOL);
        CHECK(r.curve < CURVE_TOL);
        CHECK(r.label < LABEL_TOL);
        CHECK(r.cycle < CYCLE_TOL);
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
// Cycles
// ---------------------------------------------------------------------------

TEST_CASE("Diagram cycles: states from solved variables, chains", "[diagram]") {
    const SolutionValues values = heatPumpSolution();
    const coolplot::DiagramType ph = coolplot::DiagramType::parse("ph");
    FluidModel fluid(chartFluid("R134a"));
    std::vector<std::string> warnings;

    SECTION("explicit states: full states, the solved values kept, one closed loop") {
        json definition = heatPumpCycle();
        definition["paths"] = "straight";
        const std::vector<ResolvedCycle> cycles = resolveCycles(json::array({definition}), &values, fluid, ph, warnings);
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
    SECTION("arrays, given chains with auto paths, markers alone") {
        json cycles = json::array();
        cycles.push_back({{"name", "Arrays"}, {"arrays", {{"P", "P"}, {"h", "h"}}}, {"loops", {{1, 2, 3}, {3, 4}}}});
        cycles.push_back({{"name", "Points"}, {"arrays", {{"p", "p"}, {"h", "H"}}}, {"loops", json::array()}});
        cycles.push_back({{"name", "Empty"}, {"states", json::array()}});
        const std::vector<ResolvedCycle> resolved = resolveCycles(cycles, &values, fluid, ph, warnings);
        REQUIRE(resolved.size() == 2);
        CHECK(resolved[0].vars.at(3).at("p") == "P[3]");
        CHECK(resolved[0].chains.size() == 2);
        // 3 to 4 is a throttling: an isenthalp of 20 segments
        REQUIRE(resolved[0].chains[1].size() == 21);
        CHECK(resolved[0].chains[1].front() == resolved[0].states.at(3));
        CHECK(resolved[0].chains[1].back() == resolved[0].states.at(4));
        for (const State& state : resolved[0].chains[1]) CHECK_THAT(state.at("h"), WithinRel(values.at("h_3"), 1e-8));
        CHECK(resolved[0].spec.at("paths") == json::array());  // all "auto": CoolPlot's default
        CHECK(resolved[1].index == 1);
        CHECK(resolved[1].states.size() == 4);
        CHECK(resolved[1].chains.empty());  // markers only
        CHECK(resolved[1].spec.at("loops") == json::parse("[[1,1],[2,2],[3,3],[4,4]]"));
        CHECK(resolved[1].spec.at("paths") == json::parse(R"([[[1,1],"straight"],[[2,2],"straight"],[[3,3],"straight"],[[4,4],"straight"]])"));
        CHECK(warnings.empty());
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
    request.fluid = "AirH2O";  // a fluid's chart of humid air
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("A ph diagram cannot show humid air"));
    request.fluid = "NH3H2O";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("not drawable yet"));
    request.fluid = "R134a";
    request.spec["diagram"] = "xy";
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("Unknown diagram type 'xy'"));
    request.spec["diagram"] = "psychrometric";  // a psychrometric chart of a fluid
    request.spec["isolines"] = json::array();
    request.spec["limits"] = nullptr;
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("A WT diagram cannot show R134a"));
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
    CHECK(names == std::vector<std::string>{"ph", "Ts", "hs", "ps", "prho", "Trho", "pT", "Th", "WT", "ThxW"});
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
    CHECK(o.at("fluid_types") == json::array({"Real", "Mixture", "Incompressible", "IdealGas", "HumidAir"}));
    const json& kinds = o.at("fluid_kinds");
    CHECK(kinds.at("Real").at("diagram_types").size() == 8);
    CHECK(kinds.at("Mixture").at("saturation") == true);
    CHECK(kinds.at("Real").at("default_region") == "DEF");
    CHECK(kinds.at("Mixture").at("default_region") == "ACHP");  // DEF reaches the triple point of the lightest component
    CHECK(kinds.at("IdealGas").at("default_region").is_null());
    CHECK(kinds.at("Incompressible").at("diagram_types") == json::array({"ph", "ps", "prho", "pT"}));
    CHECK(kinds.at("Incompressible").at("families") == json::array({"T", "p", "h", "s"}));
    CHECK(kinds.at("Incompressible").at("fractions") == true);
    CHECK(kinds.at("Incompressible").at("default_box") == json::array({nullptr, nullptr, 1e5, 1e7}));
    CHECK(kinds.at("IdealGas").at("diagram_types") == json::array({"ph", "Ts", "hs", "ps", "prho", "Trho", "pT"}));
    CHECK(kinds.at("IdealGas").at("saturation") == false);
    CHECK(kinds.at("IdealGas").at("default_box") == json::array({200.0, 2000.0, 1e4, 1e7}));
    // Humid air: the psychrometric and Mollier charts, their families, the region and the pressure of a new chart
    CHECK(kinds.at("HumidAir").at("diagram_types") == json::array({"WT", "ThxW"}));
    CHECK(kinds.at("HumidAir").at("families") == json::array({"T", "RH", "Twb", "h", "v"}));
    CHECK(kinds.at("HumidAir").at("chart_kind") == "humid_air");
    CHECK(kinds.at("Real").at("chart_kind") == "fluid");
    CHECK(kinds.at("HumidAir").at("total_pressure") == 101325.0);
    CHECK(kinds.at("HumidAir").at("default_region") == json{{"bounds", {{"T", {263.15, 323.15}}, {"W", {0.0, 0.03}}}}});
    const json& psychro = o.at("diagram_types")[8];
    CHECK(psychro.at("title") == "Psychrometric chart");
    CHECK(psychro.at("kind") == "humid_air");
    CHECK(psychro.at("families").size() == 4);  // RH, Twb, h, v
    const json& mollier = o.at("diagram_types")[9];
    CHECK(mollier.at("families").size() == 5);  // the isotherms too
    CHECK(mollier.at("families")[1].at("symbol") == "phi");  // the Mollier chart writes phi
    CHECK(mollier.at("families")[0].at("label") == "Dry-bulb temperature");
    const json& custom = o.at("unit_systems").back();
    CHECK(custom.at("name") == "coolsolve");
    const coolplot::UnitSystem units = coolplot::UnitSystem::from_json(custom.at("record"));
    for (const std::string& key : stateKeys()) CHECK(units.units.count(key) == 1);
    for (const std::string& key : airStateKeys()) CHECK(units.units.count(key) == 1);
    CHECK(units.at("W").label == "kg/kg");
    CHECK(units.at("RH").scale == 1.0);  // a fraction, as the model's relhum
    CHECK(units.at("Twb").offset == -273.15);
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

// ---------------------------------------------------------------------------
// Paths between the states of a cycle (thermo/cycles/paths.py)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram paths: CoolPlot's kinds between two states", "[diagram]") {
    FluidModel fluid(chartFluid("R134a"));
    // A refrigeration cycle: vapour at 2 bar, compressed to 10 bar, condensed, throttled
    const State a = fluid.statePoint("p", 2e5, "h", 4.0e5), b = fluid.statePoint("p", 1e6, "h", 4.4e5);
    const State c = fluid.statePoint("p", 1e6, "h", 2.4e5), d = fluid.statePoint("p", 2e5, "h", 2.4e5);
    const State e = fluid.statePoint("p", 1e6, "s", a.at("s")), f = fluid.statePoint("T", a.at("T"), "p", 1.5e5);

    SECTION("auto_kind: equal p, then h, s, T, else polytropic, within SAME_RTOL") {
        CHECK(autoKind(b, c) == "isobaric");
        CHECK(autoKind(c, d) == "isenthalpic");
        CHECK(autoKind(a, e) == "isentropic");
        CHECK(autoKind(a, f) == "isothermal");
        CHECK(autoKind(a, b) == "polytropic");
        State close = b;
        close["p"] = c.at("p") * (1.0 + 0.5e-6);
        CHECK(autoKind(close, c) == "isobaric");
        close["p"] = c.at("p") * (1.0 + 2e-6);
        CHECK(autoKind(close, c) == "polytropic");
    }
    SECTION("an isobar through the dome holds both saturation states") {
        const std::vector<State> path = pathStates(fluid, b, c, "auto", DEFAULT_STEPS);
        REQUIRE(path.size() == 2 + 19 + 2);  // the ends, the inner steps, the dew and bubble points
        CHECK(path.front() == b);
        CHECK(path.back() == c);
        const State dew = fluid.statePoint("p", 1e6, "Q", 1.0), bubble = fluid.statePoint("p", 1e6, "Q", 0.0);
        bool hasDew = false, hasBubble = false;
        for (size_t i = 0; i < path.size(); ++i) {
            CHECK_THAT(path[i].at("p"), WithinRel(1e6, 1e-12));
            if (i > 0) CHECK(path[i].at("h") < path[i - 1].at("h"));
            hasDew = hasDew || std::fabs(path[i].at("h") - dew.at("h")) < 1e-6;
            hasBubble = hasBubble || std::fabs(path[i].at("h") - bubble.at("h")) < 1e-6;
        }
        CHECK(hasDew);
        CHECK(hasBubble);
        // An isobar above the critical pressure crosses no dome
        FluidModel co2(chartFluid("CarbonDioxide"));
        const State hot = co2.statePoint("p", 1e7, "h", 5e5), cold = co2.statePoint("p", 1e7, "h", 3e5);
        CHECK(pathStates(co2, hot, cold, "isobaric", DEFAULT_STEPS).size() == 21);
    }
    SECTION("isenthalp, isentrope, isotherm, polytropic, straight") {
        std::vector<State> path = pathStates(fluid, c, d, "isenthalpic", DEFAULT_STEPS);
        REQUIRE(path.size() == 21);
        for (const State& s : path) CHECK_THAT(s.at("h"), WithinRel(c.at("h"), 1e-8));  // CoolProp's p-h flash gives h back to 1e-9
        CHECK_THAT(path[10].at("p"), WithinRel(std::sqrt(1e6 * 2e5), 1e-12));  // pressures spaced logarithmically
        path = pathStates(fluid, a, e, "isentropic", DEFAULT_STEPS);
        REQUIRE(path.size() == 21);
        for (const State& s : path) CHECK_THAT(s.at("s"), WithinRel(a.at("s"), 1e-9));
        path = pathStates(fluid, a, f, "isothermal", DEFAULT_STEPS);
        REQUIRE(path.size() == 21);
        for (const State& s : path) CHECK_THAT(s.at("T"), WithinRel(a.at("T"), 1e-9));
        path = pathStates(fluid, a, b, "polytropic", DEFAULT_STEPS);
        REQUIRE(path.size() == 21);
        CHECK_THAT(path[10].at("p"), WithinRel(std::sqrt(2e5 * 1e6), 1e-12));
        CHECK_THAT(path[10].at("s"), WithinRel(0.5 * (a.at("s") + b.at("s")), 1e-12));  // entropy even with the steps
        CHECK(pathStates(fluid, a, b, "straight", DEFAULT_STEPS) == std::vector<State>{a, b});
        CHECK(pathStates(fluid, a, b, "polytropic", 1) == std::vector<State>{a, b});
        CHECK_THROWS_AS(pathStates(fluid, a, b, "zigzag", DEFAULT_STEPS), DiagramError);
    }
    SECTION("cycle definitions: a kind for every connection, CoolPlot's record, steps") {
        const SolutionValues values = heatPumpSolution();
        const coolplot::DiagramType ph = coolplot::DiagramType::parse("ph");
        std::vector<std::string> warnings;
        json definition = heatPumpCycle();
        definition["paths"] = "isobaric";
        definition["steps"] = 4;
        std::vector<ResolvedCycle> cycles = resolveCycles(json::array({definition}), &values, fluid, ph, warnings);
        REQUIRE(cycles.size() == 1);
        CHECK(cycles[0].spec.at("steps") == 4);
        CHECK(cycles[0].spec.at("paths") == json::parse(R"([[[1,2],"isobaric"],[[2,3],"isobaric"],[[3,4],"isobaric"],[[4,1],"isobaric"]])"));
        definition["paths"] = json::parse(R"([[[1, 2], "isentropic"]])");
        definition["steps"] = 10;
        cycles = resolveCycles(json::array({definition}), &values, fluid, ph, warnings);
        REQUIRE(cycles.size() == 1);
        CHECK(cycles[0].spec.at("paths") == json::parse(R"([[[1,2],"isentropic"]])"));
        // 1-2 isentropic (10 segments), 2-3 an isobar from the vapour to the saturated liquid (the dew point inserted),
        // 3-4 an isenthalp, 4-1 an isobar to the saturated vapour (nothing to insert)
        CHECK(cycles[0].chains[0].size() == 1 + 10 + (9 + 1 + 1) + 10 + (9 + 1));
        CHECK_THAT(cycles[0].chains[0][5].at("s"), WithinRel(cycles[0].states.at(1).at("s"), 1e-9));
        CHECK(warnings.empty());
        definition["paths"] = "zigzag";
        CHECK_THROWS_WITH(resolveCycles(json::array({definition}), &values, fluid, ph, warnings), ContainsSubstring("Unknown path 'zigzag'"));
        definition["paths"] = json::parse(R"([[[1, 2], "spiral"]])");
        CHECK_THROWS_WITH(resolveCycles(json::array({definition}), &values, fluid, ph, warnings), ContainsSubstring("from 1 to 2"));
        definition["paths"] = "auto";
        definition["steps"] = 0;
        CHECK_THROWS_WITH(resolveCycles(json::array({definition}), &values, fluid, ph, warnings), ContainsSubstring("steps"));
    }
}

// ---------------------------------------------------------------------------
// The state at a point of a chart (POST /api/v1/diagram/state)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram state at a point of a chart", "[diagram]") {
    SECTION("inside the dome, in the vapour, on other axes") {
        const json wet = stateAt(chartFluid("R134a"), "ph", 3e5, 2e5);
        REQUIRE(wet.is_object());
        CHECK(wet.at("h") == 3e5);  // the axis values come back
        CHECK(wet.at("p") == 2e5);
        CHECK(wet.at("Q").get<double>() > 0.0);
        CHECK(wet.at("Q").get<double>() < 1.0);
        CHECK_THAT(wet.at("v").get<double>(), WithinRel(1.0 / wet.at("rho").get<double>(), 1e-12));
        for (const std::string& key : stateKeys()) CHECK(wet.at(key).is_number());
        const json vapour = stateAt(chartFluid("R134a"), "ph", 4.5e5, 2e5);
        CHECK(vapour.at("Q") == -1.0);  // single phase, as CoolPlot reports it
        const json ts = stateAt(chartFluid("Water"), "Ts", 6000.0, 500.0);
        REQUIRE(ts.is_object());
        CHECK_THAT(ts.at("T").get<double>(), WithinRel(500.0, 1e-12));
        CHECK_THAT(ts.at("s").get<double>(), WithinRel(6000.0, 1e-9));
    }
    SECTION("no state outside the fluid's range") {
        CHECK(stateAt(chartFluid("R134a"), "ph", 3e5, 1e12).is_null());  // above the highest pressure
        CHECK(stateAt(chartFluid("R134a"), "ph", -5e6, 2e5).is_null());
        CHECK(stateAt(chartFluid("Water"), "Ts", 6000.0, 5000.0).is_null());  // above the highest temperature
    }
    SECTION("requests that cannot be answered") {
        CHECK_THROWS_WITH(stateAt(chartFluid("Unobtainium"), "ph", 3e5, 2e5), ContainsSubstring("Unknown fluid"));
        CHECK_THROWS_WITH(stateAt(chartFluid("R134a"), "xy", 3e5, 2e5), ContainsSubstring("Unknown diagram type"));
        CHECK_THROWS_WITH(stateAt(chartFluid("AirH2O"), "ph", 3e5, 2e5), ContainsSubstring("cannot show humid air"));
        CHECK_THROWS_AS(stateAt(chartFluid("R134a"), "ph", NaN, 2e5), DiagramError);
    }
    SECTION("fast enough for the hover") {
        stateAt(chartFluid("R134a"), "ph", 3e5, 2e5);
        const auto start = std::chrono::steady_clock::now();
        const int n = 2000;
        const ChartFluid r134a = chartFluid("R134a");
        for (int i = 0; i < n; ++i) stateAt(r134a, "ph", 2e5 + 100.0 * i, 1e5 + 1000.0 * i);
        const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start).count() / n;
        INFO("state at a point: " << us << " us per call");
        CHECK(us < 1000.0);  // measured 9 us in a Release build; the GUI debounces its calls
    }
}

// ---------------------------------------------------------------------------
// The states of a model from its property calls (§4.2)
// ---------------------------------------------------------------------------

namespace {

IR irOf(const std::string& source) {
    EESParser parser;
    const ParseResult parsed = parser.parse(source);
    REQUIRE(parsed.success);
    return IR::fromAST(parsed.program);
}

const ModelState* stateWith(const std::vector<ModelState>& states, const std::string& var) {
    for (const ModelState& state : states)
        for (const auto& [_, name] : state.vars)
            if (name == var) return &state;
    return nullptr;
}

}  // namespace

TEST_CASE("Diagram states of a model from its property calls", "[diagram]") {
    SECTION("calls sharing two independent variables are one state; (T, p) alone or one variable is not") {
        const IR ir = irOf(
            "P_2 = 1e6\n"
            "s_2 = 1700\n"
            "h_2 = enthalpy(R134a, P=P_2, s=s_2)\n"
            "T_2 = temperature(R134a, P=P_2, h=h_2)\n"  // shares P_2 and h_2 with the call above: one state
            "P_cond = P_2\n"
            "T_3 = 30\n"
            "h_3 = enthalpy(R134a, P=P_cond, T=T_3)\n"
            "T_sat = temperature(R134a, P=P_cond, x=1)\n"  // shares only P_cond: separate, and T, p fix no state
            "T_1 = 5\n"
            "P_1 = 2e5\n"
            "h_1 = enthalpy(R134a, T=T_1, P=P_1)\n"
            "s_1 = entropy(R134a, T=T_1, P=P_1)\n"     // shares T_1 and P_1 only: not merged
            "h_su_cp = enthalpy(R134a, T=T_su_cp, P=P_su_cp)\n"
            "T_su_cp = 10\n"
            "P_su_cp = 3e5\n"
            "rho_w = density(Water, T=T_w, P=P_w)\n"
            "T_w = 20\n"
            "P_w = 1e5\n");
        const ModelStates states = groupStates(ir, nullptr, false);  // the strict rule: (T, p) never merges
        REQUIRE(states.count("R134a") == 1);
        const std::vector<ModelState>& r134a = states.at("R134a");
        const ModelState* two = stateWith(r134a, "h_2");
        REQUIRE(two != nullptr);
        CHECK(two->vars == std::map<std::string, std::string>{{"p", "P_2"}, {"s", "s_2"}, {"h", "h_2"}, {"T", "T_2"}});
        CHECK(two->number == 2);
        CHECK(two->label.empty());
        CHECK(stateWith(r134a, "T_sat") == nullptr);  // no state of its own, not in the condenser state
        CHECK(stateWith(r134a, "h_3")->vars.count("T") == 1);
        CHECK(stateWith(r134a, "h_1") != stateWith(r134a, "s_1"));  // (T, p) does not merge
        CHECK(stateWith(r134a, "h_1")->number == 1);
        CHECK(stateWith(r134a, "s_1")->label == "1");  // number 1 is taken: the next free number, labelled
        const ModelState* named = stateWith(r134a, "h_su_cp");
        CHECK(named->label == "su_cp");
        CHECK(states.at("Water").size() == 1);
        // Numbers are unique and the JSON has CoolPlot keys
        std::set<int> numbers;
        for (const ModelState& state : r134a) CHECK(numbers.insert(state.number).second);
        const json j = statesToJson(states);
        CHECK(j.at("fluids").at("R134a").size() == r134a.size());
        CHECK(j.at("fluids").at("R134a")[0].contains("label"));
        // The rule as used (the default, decided from the library validation) merges the two (T, p) calls
        // without quality, never a saturation call
        const ModelStates merged = groupStates(ir);
        CHECK(stateWith(merged.at("R134a"), "h_1") == stateWith(merged.at("R134a"), "s_1"));
        CHECK(stateWith(merged.at("R134a"), "T_sat") == nullptr);
    }
    SECTION("arrays give state i; variables without a solved value are left out") {
        const IR ir = irOf(
            "DUPLICATE i=1,2\n"
            "  P[i] = 1e5 * i\n"
            "  T[i] = 300 + 10 * i\n"
            "  h[i] = enthalpy(Water, T=T[i], P=P[i])\n"
            "  s[i] = entropy(Water, P=P[i], h=h[i])\n"
            "END\n");
        const ModelStates states = groupStates(ir);
        REQUIRE(states.at("Water").size() == 2);
        CHECK(states.at("Water")[0].number == 1);
        CHECK(states.at("Water")[1].vars == std::map<std::string, std::string>{{"T", "T[2]"}, {"p", "P[2]"}, {"h", "h[2]"}, {"s", "s[2]"}});
        SolutionValues solved = {{"P[1]", 1e5}, {"T[1]", 26.85}, {"h[1]", 1e5}};
        const ModelStates partial = groupStates(ir, &solved);
        REQUIRE(partial.at("Water").size() == 1);
        CHECK(partial.at("Water")[0].vars.size() == 3);
    }
    SECTION("numbers from the names: the number most variables carry, claimed by the groups agreeing most") {
        using Vars = std::map<std::string, std::string>;
        const std::vector<ModelState> states = numberStates({
            Vars{{"h", "h_2_s"}, {"p", "P_2"}, {"s", "s_1"}},              // no number agreed: next free
            Vars{{"T", "T_2"}, {"h", "h_2"}, {"p", "P_2"}, {"s", "s[2]"}},  // 2, array index included
            Vars{{"T", "T[1]"}, {"h", "h[1]"}, {"p", "p[1]"}, {"Q", "x_1"}},  // 1
            Vars{{"T", "T_su_cp"}, {"p", "P_su_cp"}},                       // labelled su_cp
            Vars{{"T", "T_3"}, {"p", "P_3"}, {"h", "h_x"}},                 // 3 (two of three)
            Vars{{"T", "T_3b"}, {"p", "P_3"}, {"h", "h_3"}},                // 3 too, but agrees less: labelled 3
        });
        auto withVar = [&states](const std::string& var) -> const ModelState* {
            for (const ModelState& s : states)
                for (const auto& [_, v] : s.vars)
                    if (v == var) return &s;
            return nullptr;
        };
        CHECK(withVar("h_2")->number == 2);
        CHECK(withVar("h[1]")->number == 1);
        CHECK(withVar("h_x")->number == 3);
        CHECK(withVar("T_3b")->label == "3");
        CHECK(withVar("T_3b")->number > 3);
        CHECK(withVar("h_2_s")->label.empty());
        CHECK(withVar("h_2_s")->number > 3);
        CHECK(withVar("T_su_cp")->label == "su_cp");
        std::set<int> numbers;
        for (const ModelState& s : states) CHECK(numbers.insert(s.number).second);
    }
    SECTION("the former name heuristic, for comparison") {
        IR ir = irOf("h_1 = enthalpy(R134a, T=T_1, P=P_1)\nT_1 = 5\nP_1 = 2e5\nT_sat_1 = temperature(R134a, P=P_1, x=1)\n");
        inferVariables(ir);
        const ModelStates byName = nameHeuristicStates(ir);
        REQUIRE(byName.at("R134a").size() == 1);  // "1" (T_1, P_1, h_1); "sat_1" has one property only
        CHECK(byName.at("R134a")[0].number == 1);
        CHECK(byName.at("R134a")[0].vars == std::map<std::string, std::string>{{"T", "T_1"}, {"p", "P_1"}, {"h", "h_1"}});
    }
}

// ---------------------------------------------------------------------------
// Phase 3: more fluids (blends on several threads, incompressibles, tabular backends, ideal gases)
// ---------------------------------------------------------------------------

namespace {

/// The curve of an item of a chart, computed by the shared engine (on several threads).
std::vector<PlannedCurve> chartCurves(const Chart& chart, Domain& domain, std::vector<std::string>& warnings) {
    json diagnostics;
    domain = Engine::shared().domain(chart, warnings);
    return Engine::shared().curves(chart, domain, plan(chart.spec, domain), diagnostics, warnings);
}

const Curve* curveWithId(const std::vector<PlannedCurve>& curves, const std::string& id) {
    for (const auto& [request, curve] : curves)
        if (request.item_id == id) return curve.get();
    return nullptr;
}

/// The distance of a point (x, y) to a curve of a chart, in units of the axis spans (polylineDistance).
double distanceToCurve(const Curve& curve, double x, double y, const Domain& domain, const coolplot::DiagramType& diagram) {
    const Range xr = domain.rangeSI(diagram.x), yr = domain.rangeSI(diagram.y);
    const bool xLog = diagram.x_log(), yLog = diagram.y_log();
    const double xSpan = xLog ? std::log10(xr.second / xr.first) : xr.second - xr.first;
    const double ySpan = yLog ? std::log10(yr.second / yr.first) : yr.second - yr.first;
    return polylineDistance({x}, {y}, curve.x, curve.y, xLog, yLog, xSpan, ySpan);
}

/// Solve an EES model given as text (in a scratch folder: the runner writes its outputs next to the model).
SolutionValues solveText(const std::string& name, const std::string& source, std::unique_ptr<CoolSolveRunner>& runner) {
    const fs::path dir = fs::temp_directory_path() / ("coolsolve_diagram_" + name);
    fs::create_directories(dir);
    const fs::path file = dir / (name + ".eescode");
    std::ofstream(file) << source;
    runner = std::make_unique<CoolSolveRunner>(file.string());
    REQUIRE(runner->run());
    REQUIRE(runner->isSolveSuccess());
    const auto& variables = runner->getSolveResult().variables;
    return SolutionValues(variables.begin(), variables.end());
}

/// An air-standard Brayton cycle of the ideal gas `gas` (EES syntax, CoolSolve's units: °C, Pa, J/kg).
std::string braytonModel(const std::string& gas) {
    std::string text = R"("Air-standard Brayton cycle, ideal gas GAS"
T_1 = 25
P_1 = 1e5
r_p = 10
P_2 = r_p*P_1
h_1 = enthalpy(GAS, T=T_1)
s_1 = entropy(GAS, T=T_1, P=P_1)
s_2 = s_1
h_2 = enthalpy(GAS, P=P_2, s=s_2)
T_2 = temperature(GAS, P=P_2, s=s_2)
P_3 = P_2
T_3 = 1100
h_3 = enthalpy(GAS, T=T_3)
s_3 = entropy(GAS, T=T_3, P=P_3)
P_4 = P_1
s_4 = s_3
h_4 = enthalpy(GAS, P=P_4, s=s_4)
T_4 = temperature(GAS, P=P_4, s=s_4)
w_net = (h_3 - h_4) - (h_2 - h_1)
)";
    for (size_t at; (at = text.find("GAS")) != std::string::npos;) text.replace(at, 3, gas);
    return text;
}

json braytonCycle() {
    json states = json::array();
    for (int i = 1; i <= 4; ++i) {
        const std::string n = std::to_string(i);
        states.push_back({{"number", i}, {"vars", {{"T", "T_" + n}, {"p", "P_" + n}, {"h", "h_" + n}, {"s", "s_" + n}}}});
    }
    return {{"name", "Brayton"}, {"states", states}};
}

}  // namespace

TEST_CASE("Diagram blends: the lines of one chart on several threads", "[diagram]") {
    // The lines of a chart are computed on up to MAX_CURVE_THREADS threads, one plot each: they are those of one
    // plot computing them one after the other, bit for bit
    CHECK(MAX_CURVE_THREADS == 4);
    json spec = defaultSpec("ph");
    spec["isolines"] = {{{"key", "T"}, {"num", 6}}, {{"key", "s"}, {"num", 6}}, {{"key", "Q"}, {"num", 5}}};
    const Chart chart(chartFluid("R32[0.7]&R125[0.3]"), PlanSpec::fromJson(spec));
    CHECK(chart.fluid.kind() == "Mixture");
    Engine::shared().clear();
    std::vector<std::string> warnings;
    Domain domain;
    const std::vector<PlannedCurve> parallel = chartCurves(chart, domain, warnings);
    CHECK(warnings.empty());
    REQUIRE(parallel.size() >= 15);
    Engine serial;  // an engine of its own: every line computed again, by one plot
    std::vector<std::string> serialWarnings;
    const Domain serialDomain = serial.domain(chart, serialWarnings);
    for (const auto& [request, curve] : parallel) {
        INFO(request.item_id);
        const std::shared_ptr<const Curve> one = serial.curve(chart, serialDomain, request);
        REQUIRE(one->x.size() == curve->x.size());
        CHECK(sameNumbers(one->x, curve->x));
        CHECK(sameNumbers(one->y, curve->y));
    }
    // The default p-h chart of the blend, cold: 1.55 s on one thread, 0.6 s on four (Release, 2026-11-12); the
    // budget is loose, so that only a real regression trips it on a slower machine
    SceneRequest request;
    request.fluid = "R32[0.7]&R125[0.3]";
    request.spec = defaultSpec("ph");
    Engine::shared().clear();
    const auto start = std::chrono::steady_clock::now();
    const json cold = scene(request);
    const double coldMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    INFO("R32[0.7]&R125[0.3] p-h cold: " << coldMs << " ms");
    CHECK(coldMs < 15000.0);
    Engine::shared().clear();
    CHECK(scene(request).at("svg") == cold.at("svg"));  // the same chart whichever thread computes which line
}

TEST_CASE("Diagram incompressibles: solutions at a mass fraction, absolute limits", "[diagram]") {
    SECTION("the chart fluid: CoolProp's INCOMP backend, CoolPlot's fluid record") {
        const ChartFluid meg = chartFluid("MEG", "TTSE&HEOS", {0.3});
        CHECK(meg.backend == "INCOMP");  // whatever the model's backend
        CHECK(meg.kind() == "Incompressible");
        CHECK(meg.fluidString() == "INCOMP::MEG[0.3]");
        CHECK(meg.name() == "MEG[0.3]");
        // FluidSpec.parse("INCOMP::MEG[0.3]").to_dict()
        CHECK(meg.record() == json::parse(R"({"backend": "INCOMP", "components": ["MEG"], "fractions": [0.3], "basis": "mass",
                                              "interactions": []})"));
        CHECK(ChartFluid::fromRecord(meg.record()).fluidString() == meg.fluidString());
        CHECK(chartFluid("EthyleneGlycol", "HEOS", {0.3}, "mass").fluidString() == "INCOMP::MEG[0.3]");  // an EES name
        CHECK(chartFluid("SeaWater", "HEOS", {0.035}).components == std::vector<std::string>{"MITSW"});
        CHECK(solutionFractionRange("MEG") == std::pair<double, double>{0.0, 0.6});
    }
    SECTION("what cannot be drawn: an unknown solution, a concentration missing or out of range, fractions elsewhere") {
        CHECK_THROWS_WITH(solutionFractionRange("NoSuchBrine"), ContainsSubstring("does not know the incompressible fluid"));
        CHECK_THROWS_WITH(chartFluid("NoSuchBrine", "HEOS", {0.3}), ContainsSubstring("Unknown fluid"));
        CHECK_THROWS_WITH(chartFluid("MEG"), ContainsSubstring("give its concentration as a mass fraction between 0 and 0.6"));
        CHECK_THROWS_WITH(chartFluid("MEG", "HEOS", {0.9}), ContainsSubstring("between 0 and 0.6, not 0.9"));
        CHECK_THROWS_WITH(chartFluid("MEG", "HEOS", {-0.1}), ContainsSubstring("not -0.1"));
        CHECK_THROWS_WITH(chartFluid("MEG", "HEOS", {0.2, 0.1}), ContainsSubstring("give its concentration"));
        CHECK_THROWS_WITH(chartFluid("MEG", "HEOS", {0.3}, "mole"), ContainsSubstring("mass fraction"));
        CHECK_THROWS_WITH(chartFluid("R134a", "HEOS", {0.3}), ContainsSubstring("takes no fractions"));
        CHECK_THROWS_WITH(chartFluid("N2", "HEOS", {0.3}), ContainsSubstring("takes no fractions"));
        SceneRequest request;
        request.fluid = "MEG";
        request.spec = defaultSpec("ph");
        request.fractions = {0.75};
        CHECK_THROWS_WITH(scene(request), ContainsSubstring("not 0.75"));
        CHECK_THROWS_AS(ChartFluid::fromRecord({{"backend", "INCOMP"}, {"components", {"MEG"}}, {"fractions", {0.3}}, {"basis", "mole"}}),
                        DiagramError);
    }
    SECTION("the region: the fluid's temperatures from its freezing point, 1 to 100 bar, absolute in the composed spec") {
        FluidModel meg(chartFluid("MEG", "HEOS", {0.3}));
        CHECK_FALSE(meg.hasCriticalPoint());
        std::unique_ptr<CoolProp::AbstractState> st(CoolProp::AbstractState::factory("INCOMP", "MEG"));
        st->set_mass_fractions({0.3});
        const double T_freeze = st->trivial_keyed_output(CoolProp::iT_freeze);
        CHECK_THAT(meg.tMinK(), WithinRel(T_freeze, 2e-9));
        CHECK(meg.tMinK() > T_freeze);
        const std::array<double, 4> box = resolveTpBox(meg, tpLimitsPreset("DEF"), 1e-7);
        CHECK(box == std::array<double, 4>{meg.tMinK(), meg.tMaxK(), 1e5, 1e7});
        CHECK(resolveTpBox(meg, tpLimitsPreset("NONE"), 1e-7) == box);
        // Absolute bounds are kept inside the fluid's temperatures
        const std::array<double, 4> achp = resolveTpBox(meg, tpLimitsPreset("ACHP"), 1e-7);
        CHECK(achp == std::array<double, 4>{meg.tMinK(), meg.tMaxK(), 0.25e5, 1e7});
        CHECK(resolveTpBox(meg, TpLimits::fromJson({{"T_min", 280.0}, {"T_max", 330.0}, {"p_min", 2e5}, {"p_max", 6e5}}), 1e-7) ==
              std::array<double, 4>{280.0, 330.0, 2e5, 6e5});

        SceneRequest request;
        request.fluid = "MEG";
        request.fractions = {0.3};
        request.spec = defaultSpec("ph");
        request.spec["isolines"] = {{{"key", "T"}, {"num", 8}, {"rounding", true}}, {{"key", "s"}, {"num", 8}, {"rounding", true}}};
        request.spec["marked_points"] = {"critical"};  // none to mark: no point, no error
        const json response = scene(request);
        CHECK(response.at("diagnostics").at("warnings").empty());
        CHECK(response.at("spec").at("limits") == json({{"T_min", box[0]}, {"T_max", box[1]}, {"p_min", 1e5}, {"p_max", 1e7}}));
        CHECK(response.at("spec").at("fluid").at("basis") == "mass");
        const std::string svg = response.at("svg");
        CHECK_THAT(svg, ContainsSubstring("data-item=\"iso/T/"));
        CHECK_THAT(svg, !ContainsSubstring("data-item=\"iso/Q/"));
        CHECK_THAT(svg, !ContainsSubstring("point/critical"));
        // Lines of constant quality cannot be drawn: a diagnostic each, never an error (as for any fluid without dome)
        request.spec = defaultSpec("ph");
        const json withQuality = scene(request);
        bool qualityWarned = false;
        for (const json& w : withQuality.at("diagnostics").at("warnings"))
            qualityWarned |= w.get<std::string>().find("has no two-phase region") != std::string::npos;
        CHECK(qualityWarned);
    }
    SECTION("a solved state lies on its isotherm; a state given by T and h alone is at the pressure of the solver") {
        // A state of a glycol loop as the solver evaluates `h = enthalpy(MEG, T=T, C=30)`: at 1 atm
        std::unique_ptr<CoolProp::AbstractState> st(CoolProp::AbstractState::factory("INCOMP", "MEG"));
        st->set_mass_fractions({0.3});
        st->update(CoolProp::PT_INPUTS, IncompressibleFluid::kInternalPressureSI, 300.0);
        SolutionValues values = {{"T_1", 300.0 - 273.15}, {"h_1", st->hmass()}};
        json spec = defaultSpec("ph");
        spec["isolines"] = json::array();
        spec["lines"] = {{{"name", "T1"}, {"key", "T"}, {"value_SI", 300.0}}};
        SceneRequest request;
        request.fluid = "MEG";
        request.fractions = {0.3};
        request.spec = spec;
        request.solution = &values;
        request.cycles = {{{"name", "Loop"}, {"states", {{{"number", 1}, {"vars", {{"T", "T_1"}, {"h", "h_1"}}}}}}, {"loops", json::array()}}};
        const json response = scene(request);
        REQUIRE(response.at("states").size() == 1);
        const json& state = response.at("states")[0].at("values_SI");
        CHECK(state.at("p") == IncompressibleFluid::kInternalPressureSI);
        CHECK(state.at("h") == values.at("h_1"));
        const Chart chart(chartFluid("MEG", "HEOS", {0.3}), PlanSpec::fromJson(spec));
        std::vector<std::string> warnings;
        Domain domain;
        const std::vector<PlannedCurve> curves = chartCurves(chart, domain, warnings);
        const Curve* isotherm = curveWithId(curves, "line/T1");
        REQUIRE(isotherm);
        CHECK(distanceToCurve(*isotherm, values.at("h_1"), IncompressibleFluid::kInternalPressureSI, domain, chart.spec.diagram) < 1e-4);
        // The state under the cursor of a chart of the solution
        const json point = stateAt(chartFluid("MEG", "HEOS", {0.3}), "ph", values.at("h_1"), IncompressibleFluid::kInternalPressureSI);
        REQUIRE(point.is_object());
        CHECK_THAT(point.at("T").get<double>(), WithinRel(300.0, 1e-9));
    }
    SECTION("the concentration of a model's solutions, from its property calls") {
        const IR ir = irOf("T_1 = 20\nc_glycol = 35\nh_1 = enthalpy(MEG, T=T_1, C=30)\ncp_1 = cp(MPG, T=T_1, C=c_glycol)\n"
                           "rho_1 = density(SeaWater, 3.5, T=T_1)\nh_2 = enthalpy(R134a, T=T_1, x=1)\n");
        std::map<std::string, double> fractions = solutionFractions(ir);
        CHECK(fractions.size() == 2);  // the variable c_glycol has no value without a solve
        CHECK_THAT(fractions.at("MEG"), WithinRel(0.3, 1e-12));
        CHECK_THAT(fractions.at("SeaWater"), WithinRel(0.035, 1e-12));
        const SolutionValues solved = {{"c_glycol", 35.0}};
        fractions = solutionFractions(ir, &solved);
        CHECK_THAT(fractions.at("MPG"), WithinRel(0.35, 1e-12));
        const json document = statesToJson({}, fractions);
        CHECK(document.at("fractions").at("MEG") == json::array({fractions.at("MEG")}));
        CHECK(statesToJson({}).contains("fractions") == false);
    }
}

TEST_CASE("Diagram tabular backends: the lines of constant quality from the equation of state", "[diagram]") {
    // A model whose coolsolve.conf sets coolpropBackend = BICUBIC&HEOS (or TTSE&HEOS): the lines with the tables, the
    // lines of constant quality with the equation of state behind them, as CoolPlot's native.port_fluid_string
    for (const char* backend : {"BICUBIC&HEOS", "TTSE&HEOS"}) {
        INFO(backend);
        const ChartFluid water = chartFluid("Water", backend);
        CHECK(water.kind() == "Real");
        CHECK(water.fluidString() == std::string(backend) + "::Water");
        CHECK(water.portFluidString("p") == water.fluidString());
        CHECK(water.portFluidString("Q") == "HEOS::Water");
        json spec = defaultSpec("Ts");
        spec["isolines"] = {{{"key", "p"}, {"values_SI", {1e5, 1e6}}}, {{"key", "Q"}, {"values_SI", {0.5}}}};
        const Chart chart(water, PlanSpec::fromJson(spec));
        std::vector<std::string> warnings;
        Domain domain;
        const std::vector<PlannedCurve> curves = chartCurves(chart, domain, warnings);
        CHECK(warnings.empty());
        for (const char* id : {"iso/p/100000", "iso/p/1000000", "iso/Q/0.5", "iso/Q/0", "iso/Q/1"}) {
            INFO(id);
            const Curve* curve = curveWithId(curves, id);
            REQUIRE(curve);
            CHECK(curve->nValid() > 20);
        }
        // A state of the model's backend lies on its isobar
        std::unique_ptr<CoolProp::AbstractState> st(CoolProp::AbstractState::factory(backend, "Water"));
        st->update(CoolProp::PT_INPUTS, 1e6, 500.0);
        CHECK(distanceToCurve(*curveWithId(curves, "iso/p/1000000"), st->smass(), 500.0, domain, chart.spec.diagram) < 1e-4);
        // A whole chart, and the hover with the model's backend
        SceneRequest request;
        request.fluid = "Water";
        request.backend = backend;
        request.spec = defaultSpec("Ts");
        const json response = scene(request);
        CHECK_THAT(response.at("svg").get<std::string>(), ContainsSubstring("data-item=\"iso/Q/0.5\""));
        CHECK(response.at("spec").at("fluid").at("backend") == backend);
        const json point = stateAt(water, "Ts", st->smass(), 500.0);
        REQUIRE(point.is_object());
        CHECK_THAT(point.at("p").get<double>(), WithinRel(1e6, 1e-3));  // the (s, T) flash of the tables
    }
}

TEST_CASE("Diagram ideal gases: CoolSolve's provider on the EES/JANAF basis", "[diagram]") {
    SECTION("the provider: the solver's ideal-gas functions with the reference offsets") {
        const ChartFluid co2 = chartFluid("CO2", "TTSE&HEOS");
        CHECK(co2.isIdealGas());
        CHECK(co2.kind() == "IdealGas");
        CHECK(co2.name() == "CO2");
        CHECK(co2.fluidString() == "IDEALGAS::CO2");
        CHECK(co2.record().at("components") == json::array({"CarbonDioxide"}));
        CHECK(co2.record().at("coolsolve") == json({{"model", "ideal_gas"}, {"fluid", "CO2"}}));
        CHECK(ChartFluid::fromRecord(co2.record()).idealGas == "CO2");
        IdealGasModel gas(co2);
        // At 25 °C and 1 bar: the enthalpy of formation and the standard entropy of the NASA tables (fluids.cpp)
        const double M = 44.0095;
        const State standard = gas.statePoint("T", 298.15, "p", 1e5);
        CHECK_THAT(standard.at("h"), WithinRel(-393510e3 / M, 1e-4));
        CHECK_THAT(standard.at("s"), WithinRel(213786.2 / M, 1e-4));
        CHECK_THAT(standard.at("rho") * standard.at("v"), WithinRel(1.0, 1e-14));
        CHECK(std::isnan(standard.at("Q")));
        // Any independent pair gives the same state back
        for (const auto& [a, b] : {std::pair<const char*, const char*>{"p", "h"}, {"p", "s"}, {"T", "s"}, {"rho", "s"}, {"h", "rho"}, {"v", "u"}}) {
            INFO(a << ", " << b);
            const State again = gas.statePoint(a, standard.at(a), b, standard.at(b));
            CHECK_THAT(again.at("T"), WithinRel(298.15, 1e-9));
            CHECK_THAT(again.at("p"), WithinRel(1e5, 1e-9));
        }
        CHECK_THROWS(gas.update("T", 300.0, "h", standard.at("h")));  // both fix only the temperature
        CHECK_THROWS(gas.update("p", 1e5, "Q", 0.5));
        CHECK_FALSE(gas.hasCriticalPoint());
        CHECK_THROWS_AS(gas.saturationRange("T", 1e-7), DiagramError);
        CHECK(resolveTpBox(gas, tpLimitsPreset("DEF"), 1e-7) == IDEAL_GAS_BOX);
        CHECK(resolveTpBox(gas, TpLimits::fromJson({{"T_max", 3000.0}}), 1e-7) == std::array<double, 4>{50.0, 3000.0, 1e4, 1e7});
    }
    SECTION("the charts it has, and those it has not") {
        SceneRequest request;
        request.fluid = "N2";
        request.spec = defaultSpec("Ts");
        request.spec["isolines"] = {{{"key", "p"}, {"num", 8}}, {{"key", "rho"}, {"num", 6}}, {{"key", "h"}, {"num", 6}}};
        const json ts = scene(request);
        CHECK(ts.at("diagnostics").at("warnings").empty());
        CHECK(ts.at("spec").at("limits") == json({{"T_min", 200.0}, {"T_max", 2000.0}, {"p_min", 1e4}, {"p_max", 1e7}}));
        CHECK(ts.at("spec").at("fluid").at("coolsolve").at("fluid") == "N2");
        for (const json& c : ts.at("diagnostics").at("curves")) CHECK(c.at("method") == "ideal_gas");
        const std::string svg = ts.at("svg");
        CHECK_THAT(svg, ContainsSubstring("data-item=\"iso/p/10000\""));
        CHECK_FALSE(std::regex_search(svg, std::regex("[^A-Za-z](nan|NaN|inf|Infinity)[^A-Za-z]")));
        for (const char* diagram : {"ph", "hs", "ps", "prho", "Trho", "pT"}) {
            INFO(diagram);
            request.spec = defaultSpec(diagram);
            json isolines = json::array();
            for (const json& iso : request.spec.at("isolines"))
                if (iso.at("key") != "Q") isolines.push_back(iso);
            request.spec["isolines"] = isolines;
            const json response = scene(request);
            CHECK(response.at("diagnostics").at("warnings").empty());
            size_t drawn = 0;
            for (const json& c : response.at("diagnostics").at("curves")) drawn += c.at("drawn").get<bool>();
            CHECK(drawn >= 20);
        }
        // Quality lines: a diagnostic, as for any fluid without dome
        request.spec = defaultSpec("Ts");
        CHECK_THAT(scene(request).at("diagnostics").at("warnings")[0].get<std::string>(), ContainsSubstring("no two-phase region"));
        // T-h: every state on one line
        request.spec = defaultSpec("Th");
        CHECK_THROWS_WITH(scene(request), ContainsSubstring("An ideal gas has no T-h chart"));
        CHECK_THROWS_WITH(stateAt(chartFluid("N2"), "Th", 3e5, 400.0), ContainsSubstring("no T-h chart"));
        const json point = stateAt(chartFluid("N2"), "Ts", 7000.0, 400.0);
        REQUIRE(point.is_object());
        CHECK(point.at("Q").is_null());
        CHECK_THAT(point.at("T").get<double>(), WithinRel(400.0, 1e-12));
    }
    SECTION("the states of an ideal-gas model: two of T, h, u fix no state") {
        CHECK(independentKeys("T", "h"));
        CHECK_FALSE(independentKeys("T", "h", true));
        CHECK_FALSE(independentKeys("u", "h", true));
        CHECK(independentKeys("p", "h", true));
        const ModelStates states = groupStates(irOf(braytonModel("N2")));
        REQUIRE(states.count("N2") == 1);
        // h_1 = enthalpy(N2, T=T_1) shares T_1 alone with state 1, and (T_1, h_1) is no state of an ideal gas
        const ModelState* first = stateWith(states.at("N2"), "s_1");
        REQUIRE(first);
        CHECK(first->vars == std::map<std::string, std::string>{{"T", "T_1"}, {"p", "P_1"}, {"s", "s_1"}});
        CHECK(stateWith(states.at("N2"), "h_1") == nullptr);
        CHECK(stateWith(states.at("N2"), "h_2") != nullptr);  // (p, s, h) of the compressor outlet
    }
    SECTION("a solved air-standard Brayton cycle lies on its isobars and isentropes") {
        for (const char* gas : {"Air", "N2"}) {
            INFO(gas);
            std::unique_ptr<CoolSolveRunner> runner;
            const SolutionValues values = solveText(std::string("brayton_") + gas, braytonModel(gas), runner);
            CHECK(values.at("w_net") > 0.0);
            const json limits = {{"T_min", 250.0}, {"T_max", 1500.0}, {"p_min", 5e4}, {"p_max", 2e6}};
            // The isobars are curves of the T-s chart, the isentropes curves of the p-h chart
            const json lines = {{"Ts", {{{"name", "p1"}, {"key", "p"}, {"value_SI", values.at("P_1")}},
                                        {{"name", "p2"}, {"key", "p"}, {"value_SI", values.at("P_2")}}}},
                                {"ph", {{{"name", "s1"}, {"key", "s"}, {"value_SI", values.at("s_1")}},
                                        {{"name", "s3"}, {"key", "s"}, {"value_SI", values.at("s_3")}}}}};
            const std::map<std::string, std::vector<std::string>> through = {{"Ts", {"line/p1", "line/p2", "line/p2", "line/p1"}},
                                                                             {"ph", {"line/s1", "line/s1", "line/s3", "line/s3"}}};
            for (const char* diagram : {"Ts", "ph"}) {
                INFO(diagram);
                json spec = defaultSpec(diagram);
                spec["isolines"] = json::array();
                spec["limits"] = limits;
                spec["lines"] = lines.at(diagram);
                SceneRequest request;
                request.fluid = gas;
                request.spec = spec;
                request.cycles = json::array({braytonCycle()});
                request.solution = &values;
                const json response = scene(request);
                CHECK(response.at("diagnostics").at("warnings").empty());
                REQUIRE(response.at("states").size() == 4);
                CHECK_THAT(response.at("svg").get<std::string>(), ContainsSubstring("data-item=\"cycle/Brayton/loop/0\""));
                CHECK(response.at("spec").at("cycles")[0].at("paths") == json::array());  // auto: CoolPlot's isentropes, isobars
                const Chart chart(chartFluid(gas), PlanSpec::fromJson(spec));
                std::vector<std::string> warnings;
                Domain domain;
                const std::vector<PlannedCurve> curves = chartCurves(chart, domain, warnings);
                const coolplot::DiagramType& type = chart.spec.diagram;
                for (int i = 1; i <= 4; ++i) {
                    const std::string n = std::to_string(i), id = through.at(diagram)[i - 1];
                    INFO("state " << i << " on " << id);
                    const json& drawn = response.at("states")[i - 1].at("values_SI");
                    // The marker sits where the solver put it, and the provider agrees with the solver
                    CHECK(drawn.at("s") == values.at("s_" + n));
                    CHECK_THAT(drawn.at("h").get<double>(), WithinRel(values.at("h_" + n), 1e-9));
                    const Curve* line = curveWithId(curves, id);
                    REQUIRE(line);
                    CHECK(distanceToCurve(*line, drawn.at(type.x).get<double>(), drawn.at(type.y).get<double>(), domain, type) < 1e-4);
                }
                // The path 2-3 is drawn along the isobar of state 2 (CoolPlot's auto kind)
                if (std::string(diagram) == "Ts") {
                    const Curve* isobar = curveWithId(curves, "line/p2");
                    std::vector<std::string> unused;
                    const std::vector<ResolvedCycle> resolved = Engine::shared().withFluid(chart.fluid, [&](PropertyProvider& model) {
                        return resolveCycles(json::array({braytonCycle()}), &values, model, type, unused);
                    });
                    REQUIRE(resolved.size() == 1);
                    const std::vector<State>& loop = resolved[0].chains.at(0);
                    REQUIRE(loop.size() == 4 * DEFAULT_STEPS + 1);
                    for (int k = DEFAULT_STEPS; k <= 2 * DEFAULT_STEPS; ++k)
                        CHECK(distanceToCurve(*isobar, loop[k].at("s"), loop[k].at("T"), domain, type) < 1e-4);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Humid air: the psychrometric and Mollier charts (CS-FEAT-PSYCHRO)
// ---------------------------------------------------------------------------

namespace {

/// A chart of humid air as the Diagrams tab writes it: every family, round values, labelled (CoolPlot's defaults for
/// humid air), PsychrometricChart's region (no limits), no badge.
json humidAirSpec(const std::string& diagram) {
    const coolplot::DiagramType type = coolplot::DiagramType::parse(diagram);
    json isolines = json::array();
    for (const coolplot::Family& family : type.families())
        isolines.push_back({{"key", family.key}, {"num", family.default_num}, {"rounding", true}, {"labels", true}});
    return {{"diagram", diagram},
            {"units", "engineering"},
            {"isolines", isolines},
            {"show_saturation", true},
            {"settings", {{"branding", {{"references", "none"}, {"versions", false}}}}}};
}

/// Air conditioning: outdoor air cooled and dehumidified on a coil, then reheated at constant humidity ratio (EES
/// syntax, CoolSolve's units: °C, Pa, J/kg, kg/kg and a fraction for the relative humidity).
const char* const AIR_CONDITIONING = R"("Air conditioning: cooling and dehumidification, then reheating"
P = 101325
T[1] = 32
RH[1] = 0.5
w[1] = humrat(AirH2O, T=T[1], P=P, R=RH[1])
h[1] = enthalpy(AirH2O, T=T[1], P=P, w=w[1])
T_wb[1] = wetbulb(AirH2O, T=T[1], P=P, w=w[1])
T[2] = 13
RH[2] = 0.95
w[2] = humrat(AirH2O, T=T[2], P=P, R=RH[2])
h[2] = enthalpy(AirH2O, T=T[2], P=P, R=RH[2])
w[3] = w[2]
T[3] = 20
h[3] = enthalpy(AirH2O, T=T[3], P=P, w=w[3])
RH[3] = relhum(AirH2O, T=T[3], P=P, w=w[3])
v[3] = volume(AirH2O, T=T[3], P=P, w=w[3])
q_cool = h[1] - h[2]
q_reheat = h[3] - h[2]
)";

}  // namespace

TEST_CASE("Diagram humid air: the provider, the region, the lines", "[diagram]") {
    SECTION("states of humid air from two properties at the chart's total pressure (air_state)") {
        const State s = airState("T", 298.15, "RH", 0.5, 101325.0);
        for (const std::string& key : airStateKeys()) CHECK(std::isfinite(s.at(key)));
        CHECK(s.at("T") == 298.15);  // HAPropsSI hands the inputs back
        CHECK(s.at("RH") == 0.5);
        CHECK(s.at("p") == 101325.0);
        CHECK_THAT(s.at("W"), WithinRel(HumidAir::HAPropsSI("W", "T", 298.15, "P", 101325.0, "R", 0.5), 1e-15));
        CHECK_THAT(s.at("Thx"), WithinRel(273.15 + (s.at("h") - 2500e3 * s.at("W")) / 1006.0, 1e-15));
        // The same state from W and the Mollier temperature scale, and from T and h
        const State m = airState("Thx", s.at("Thx"), "W", s.at("W"), 101325.0);
        CHECK_THAT(m.at("T"), WithinRel(298.15, 1e-9));
        CHECK_THAT(airState("T", 298.15, "h", s.at("h"), 101325.0).at("W"), WithinRel(s.at("W"), 1e-9));
        CHECK_THROWS_WITH(airState("T", 298.15, "T", 300.0, 101325.0), ContainsSubstring("two different properties"));
        CHECK_THROWS_WITH(airState("T", 298.15, "p", 1e5, 101325.0), ContainsSubstring("total pressure of the chart"));
        CHECK_THROWS_WITH(airState("Thx", 298.15, "T", 300.0, 101325.0), ContainsSubstring("together with the humidity ratio"));
        CHECK_THROWS_AS(airState("T", 298.15, "s", 100.0, 101325.0), DiagramError);
        CHECK(std::isnan(saturationW(500.0, 101325.0)));  // above the boiling point: no saturation
    }
    SECTION("the chart fluid: AirH2O at a total pressure, CoolPlot's HumidAirSpec record") {
        const ChartFluid air = chartFluid("AirH2O");
        CHECK(air.isHumidAir());
        CHECK(air.p_Pa == 101325.0);
        CHECK(air.kind() == "HumidAir");
        CHECK(air.name() == "humid air at 101325 Pa");
        CHECK(air.record() == json{{"model", "humid_air"}, {"p_Pa", 101325.0}});
        const ChartFluid low = chartFluid("AirH2O", "HEOS", {}, "", 9e4);
        CHECK(low.fluidString() != air.fluidString());  // separate caches
        CHECK(ChartFluid::fromRecord(low.record()).p_Pa == 9e4);
        CHECK_THROWS_WITH(chartFluid("AirH2O", "HEOS", {}, "", 5e3), ContainsSubstring("between 10000.0 and 10000000.0 Pa"));
        CHECK_THROWS_WITH(chartFluid("R134a", "HEOS", {}, "", 1e5), ContainsSubstring("takes no total pressure"));
        CHECK_THROWS_AS(ChartFluid::fromRecord({{"model", "humid_air"}, {"p_Pa", "sea level"}}), DiagramError);
    }
    SECTION("the region: PsychrometricChart's default, the dry-bulb span of the Mollier chart, CoolPlot's checks") {
        const PlanSpec psychro = PlanSpec::fromJson(humidAirSpec("WT"));
        CHECK(psychro.region.toJson() == json{{"bounds", {{"T", {263.15, 323.15}}, {"W", {0.0, 0.03}}}}});
        const HumidAirModel air(101325.0);
        const Domain d = air.domain(psychro.region);
        CHECK(d.saturation);
        CHECK(d.tp_box == std::array<double, 4>{263.15, 323.15, 101325.0, 101325.0});
        CHECK(d.rangeSI("RH") == Range{0.0, 1.0});
        for (const char* key : {"Twb", "h", "v"}) CHECK(d.rangeSI(key).second > d.rangeSI(key).first);
        // The Mollier chart is bounded by its temperature scale; the dry-bulb temperatures behind it reach further
        const PlanSpec mollier = PlanSpec::fromJson(humidAirSpec("ThxW"));
        CHECK(mollier.region.has("Thx"));
        const Domain dm = air.domain(mollier.region);
        CHECK(dm.rangeSI("Thx") == Range{263.15, 323.15});
        CHECK(dm.tp_box[0] < 263.15);
        CHECK(dm.tp_box[1] > 323.15);
        CHECK_THROWS_WITH(air.domain(Region::fromJson({{"bounds", {{"T", {100.0, 300.0}}, {"W", {0.0, 0.03}}}}})),
                          ContainsSubstring("between 173.15 and 623.15 K"));
        CHECK_THROWS_WITH(air.domain(Region::fromJson({{"bounds", {{"T", {280.0, 300.0}}, {"W", {-0.01, 0.03}}}}})),
                          ContainsSubstring("at least 0"));
        CHECK_THROWS_WITH(air.domain(Region::fromJson({{"bounds", {{"T", {263.15, 273.15}}, {"W", {0.02, 0.03}}}}})),
                          ContainsSubstring("entirely above saturation (fog)"));
        CHECK_THROWS_WITH(Region::fromJson({{"bounds", {{"T", {300.0, 300.0}}, {"W", {0.0, 0.03}}}}}), ContainsSubstring("two different"));
        json spec = humidAirSpec("WT");
        spec["limits"] = "DEF";
        CHECK_THROWS_WITH(PlanSpec::fromJson(spec), ContainsSubstring("needs a Region with bounds on 'T' and 'W'"));
        spec["limits"] = {{"bounds", {{"T", {263.15, 323.15}}}}};
        CHECK_THROWS_WITH(PlanSpec::fromJson(spec), ContainsSubstring("needs a Region"));
        spec = humidAirSpec("WT");
        spec["marked_points"] = {"critical"};
        CHECK_THROWS_WITH(PlanSpec::fromJson(spec), ContainsSubstring("only in charts of a fluid"));
        spec = humidAirSpec("WT");
        spec["isolines"][0]["anchor"] = {{"by", "T"}, {"step_SI", 10.0}};
        CHECK_THROWS_WITH(PlanSpec::fromJson(spec), ContainsSubstring("only in charts of a fluid"));
        spec = humidAirSpec("WT");
        spec["saturation_scale"] = {{"Q", 1.0}};
        CHECK_THROWS_AS(PlanSpec::fromJson(spec), DiagramError);
    }
    SECTION("the lines: on their own value, never above saturation, starting on it; isotherms of the Mollier chart") {
        const Chart chart(chartFluid("AirH2O"), PlanSpec::fromJson(humidAirSpec("WT")));
        std::vector<std::string> warnings;
        Domain domain;
        const std::vector<PlannedCurve> curves = chartCurves(chart, domain, warnings);
        CHECK(warnings.empty());
        REQUIRE(curveWithId(curves, "iso/RH/1"));  // the saturation curve
        CHECK_FALSE(curveWithId(curves, "iso/RH/0"));  // the dry-air axis is no line
        const double p = 101325.0;
        for (const auto& [request, curve] : curves) {
            INFO(request.item_id);
            CHECK(curve->nGaps() == 0);
            for (size_t i = 0; i < curve->x.size(); ++i) {
                const double T = curve->x[i], W = curve->y[i];
                CHECK(W >= 0.0);
                CHECK(W <= saturationW(T, p) * (1.0 + 1e-6));
                // A line of constant RH, h or v holds its value (relative to the range of its family); the saturation
                // line is the humidity ratio of saturated air
                if (request.item_id == "iso/RH/1") {
                    CHECK_THAT(W, WithinRel(saturationW(T, p), 1e-12));
                } else if (W > 1e-6 && request.key != "Twb") {
                    const Range range = domain.rangeSI(request.key);
                    CHECK(std::fabs(humidAir(request.key, "T", T, "W", W, p) - request.value_SI) <= 1e-6 * (range.second - range.first));
                }
            }
        }
        // A wet-bulb line starts on saturation at its own temperature
        const Curve* wetBulb = curveWithId(curves, "iso/Twb/293.15");
        REQUIRE(wetBulb);
        CHECK(wetBulb->x.front() == 293.15);
        CHECK_THAT(wetBulb->y.front(), WithinRel(saturationW(293.15, p), 1e-12));

        const Chart mollier(chartFluid("AirH2O"), PlanSpec::fromJson(humidAirSpec("ThxW")));
        const std::vector<PlannedCurve> mollierCurves = chartCurves(mollier, domain, warnings);
        CHECK(warnings.empty());
        const Curve* isotherm = curveWithId(mollierCurves, "iso/T/293.15");
        REQUIRE(isotherm);
        CHECK(isotherm->x.front() == 0.0);  // from dry air (at Thx = T, but for the real enthalpy of dry air) to saturation
        CHECK_THAT(isotherm->y.front(), Catch::Matchers::WithinAbs(293.15, 0.01));
        CHECK_THAT(isotherm->x.back(), WithinRel(saturationW(293.15, p), 1e-12));
        for (size_t i = 0; i < isotherm->x.size(); ++i)
            CHECK_THAT(humidAir("T", "h", enthalpyFromMollier(isotherm->y[i], isotherm->x[i]), "W", isotherm->x[i], p), WithinRel(293.15, 1e-7));
    }
}

TEST_CASE("Diagram humid air: the states of a model's humid-air calls", "[diagram]") {
    CHECK(independentAirKeys("T", "W"));
    CHECK(independentAirKeys("h", "RH"));
    CHECK_FALSE(independentAirKeys("T", "p"));  // the pressure is the chart's
    CHECK_FALSE(independentAirKeys("W", "Tdp"));  // at one pressure, the dew point is the humidity ratio
    CHECK_FALSE(independentAirKeys("RH", "RH"));
    // Two calls at the same T and P do not describe the same humid air; two that share T and w do
    const ModelStates states = groupStates(irOf(R"(P = 101325
T_1 = 30
h_a = enthalpy(AirH2O, T=T_1, P=P, w=w_a)
h_b = enthalpy(AirH2O, T=T_1, P=P, R=RH_b)
RH_a = relhum(AirH2O, T=T_1, P=P, w=w_a)
T_dp = dewpoint(AirH2O, T=T_d, P=P, w=w_a)
)"));
    REQUIRE(states.count("AirH2O") == 1);
    const ModelState* a = stateWith(states.at("AirH2O"), "h_a");
    REQUIRE(a);
    CHECK(a->vars == std::map<std::string, std::string>{{"T", "T_1"}, {"p", "P"}, {"W", "w_a"}, {"h", "h_a"}, {"RH", "RH_a"}});
    const ModelState* b = stateWith(states.at("AirH2O"), "h_b");
    REQUIRE(b);
    CHECK(b != a);
    CHECK(b->vars.count("W") == 0);
    // The dew-point call shares w_a alone with state a (and P): a state of its own, (T_d, w_a)
    const ModelState* d = stateWith(states.at("AirH2O"), "T_dp");
    REQUIRE(d);
    CHECK(d != a);
}

TEST_CASE("Diagram humid air: charts, a solved process, the state under the cursor", "[diagram]") {
    SECTION("the psychrometric and Mollier charts end to end; the region and the fluid in the spec CoolPlot reads") {
        for (const char* diagram : {"WT", "ThxW"}) {
            INFO(diagram);
            SceneRequest request;
            request.fluid = "AirH2O";
            request.spec = humidAirSpec(diagram);
            request.figureSizePt = std::make_pair(600.0, 420.0);
            const json response = scene(request);
            const std::string svg = response.at("svg");
            for (const char* id : {"iso/RH/1", "iso/RH/0.5", "iso/h/50000"})
                CHECK_THAT(svg, ContainsSubstring(std::string("data-item=\"") + id + "\""));
            CHECK_FALSE(std::regex_search(svg, std::regex("[^A-Za-z](nan|NaN|inf|Infinity)[^A-Za-z]")));
            CHECK(response.at("diagnostics").at("warnings").empty());
            CHECK(response.at("spec").at("fluid") == json{{"model", "humid_air"}, {"p_Pa", 101325.0}});
            CHECK(response.at("spec").at("limits").at("bounds").contains(std::string(diagram) == "WT" ? "T" : "Thx"));
            const json& axes = response.at("axes");
            CHECK(axes.at("x_key") == (std::string(diagram) == "WT" ? "T" : "W"));
            CHECK(axes.at("y_unit").at("scale") == (std::string(diagram) == "WT" ? 1e3 : 1.0));  // g/kg, °C
        }
        // The badge names the humid-air model, as CoolPlot's does
        SceneRequest request;
        request.fluid = "AirH2O";
        request.spec = humidAirSpec("WT");
        request.spec["settings"].erase("branding");
        CHECK_THAT(scene(request).at("svg").get<std::string>(), ContainsSubstring("Herrmann"));
        // The model's units hold the quantities of humid air
        request.spec["units"] = coolsolveUnitSystem();
        CHECK(scene(request).at("axes").at("y_unit").at("label") == "kg/kg");
    }
    SECTION("a solved air-conditioning process: its states grouped, on the lines of their properties, straight paths") {
        std::unique_ptr<CoolSolveRunner> runner;
        const SolutionValues values = solveText("air_conditioning", AIR_CONDITIONING, runner);
        // The states of the model: calls sharing two independent humid-air variables, the pressure aside
        const ModelStates states = groupStates(runner->getIR(), &values);
        REQUIRE(states.count("AirH2O") == 1);
        const std::vector<ModelState>& air = states.at("AirH2O");
        REQUIRE(air.size() == 3);
        CHECK(air[0].vars == std::map<std::string, std::string>{{"T", "T[1]"}, {"p", "P"}, {"RH", "RH[1]"}, {"W", "w[1]"},
                                                                {"h", "h[1]"}, {"Twb", "T_wb[1]"}});
        CHECK(air[1].vars == std::map<std::string, std::string>{{"T", "T[2]"}, {"p", "P"}, {"RH", "RH[2]"}, {"W", "w[2]"}, {"h", "h[2]"}});
        CHECK(air[2].vars == std::map<std::string, std::string>{{"T", "T[3]"}, {"p", "P"}, {"W", "w[3]"}, {"h", "h[3]"},
                                                                {"RH", "RH[3]"}, {"v", "v[3]"}});
        json cycleStates = json::array();
        for (const ModelState& state : air) cycleStates.push_back({{"number", state.number}, {"vars", state.vars}});
        const json cycle = {{"name", "Air conditioning"}, {"states", cycleStates}, {"loops", {{1, 2, 3}}}};
        // The lines of the states' own properties
        const json lines = {{{"name", "RH2"}, {"key", "RH"}, {"value_SI", values.at("RH[2]")}},
                            {{"name", "RH3"}, {"key", "RH"}, {"value_SI", values.at("RH[3]")}},
                            {{"name", "h1"}, {"key", "h"}, {"value_SI", values.at("h[1]")}},
                            {{"name", "Twb1"}, {"key", "Twb"}, {"value_SI", values.at("T_wb[1]") + 273.15}},
                            {{"name", "v3"}, {"key", "v"}, {"value_SI", values.at("v[3]")}}};
        const std::map<int, std::vector<std::string>> through = {{1, {"iso/RH/0.5", "line/h1", "line/Twb1"}}, {2, {"line/RH2"}},
                                                                 {3, {"line/RH3", "line/v3"}}};
        for (const char* diagram : {"WT", "ThxW"}) {
            INFO(diagram);
            json spec = humidAirSpec(diagram);
            spec["lines"] = lines;
            SceneRequest request;
            request.fluid = "AirH2O";
            request.spec = spec;
            request.cycles = json::array({cycle});
            request.solution = &values;
            const json response = scene(request);
            CHECK(response.at("diagnostics").at("warnings").empty());
            REQUIRE(response.at("states").size() == 3);
            CHECK_THAT(response.at("svg").get<std::string>(), ContainsSubstring("data-item=\"cycle/Air conditioning/loop/0\""));
            const Chart chart(chartFluid("AirH2O"), PlanSpec::fromJson(spec));
            std::vector<std::string> warnings;
            Domain domain;
            const std::vector<PlannedCurve> curves = chartCurves(chart, domain, warnings);
            const coolplot::DiagramType& type = chart.spec.diagram;
            for (const auto& [number, ids] : through) {
                const std::string n = "[" + std::to_string(number) + "]";
                const json& drawn = response.at("states")[number - 1].at("values_SI");
                // The marker sits where the solver put it
                CHECK(drawn.at("W") == values.at("w" + n));
                CHECK(drawn.at("T") == values.at("T" + n) + 273.15);
                CHECK_THAT(drawn.at("h").get<double>(), WithinRel(values.at("h" + n), 1e-9));
                CHECK_THAT(drawn.at("Thx").get<double>(), WithinRel(mollierTemperature(drawn.at("h"), drawn.at("W")), 1e-6));
                for (const std::string& id : ids) {
                    INFO("state " << number << " on " << id);
                    const Curve* line = curveWithId(curves, id);
                    REQUIRE(line);
                    CHECK(distanceToCurve(*line, drawn.at(type.x).get<double>(), drawn.at(type.y).get<double>(), domain, type) < 1e-4);
                }
            }
            // Air processes are straight lines between their states (HumidAir.path_kind)
            std::vector<std::string> unused;
            const std::vector<ResolvedCycle> resolved = Engine::shared().withFluid(chart.fluid, [&](PropertyProvider& model) {
                return resolveCycles(json::array({cycle}), &values, model, type, unused);
            });
            REQUIRE(resolved.size() == 1);
            CHECK(resolved[0].chains.at(0).size() == 3);
        }
        // A chart at another total pressure does not draw the states of the model
        SceneRequest request;
        request.fluid = "AirH2O";
        request.p_Pa = 1e5;
        request.spec = humidAirSpec("WT");
        request.cycles = json::array({cycle});
        request.solution = &values;
        const json response = scene(request);
        CHECK(response.at("states").empty());
        REQUIRE(response.at("diagnostics").at("warnings").size() == 1);
        CHECK_THAT(response.at("diagnostics").at("warnings")[0].get<std::string>(),
                   ContainsSubstring("is at 101325 Pa, but the chart shows humid air at 100000 Pa"));
        CHECK(response.at("spec").at("fluid").at("p_Pa") == 1e5);
    }
    SECTION("the state under the cursor") {
        const ChartFluid air = chartFluid("AirH2O");
        const json s = stateAt(air, "WT", 298.15, 0.01);
        REQUIRE(s.is_object());
        for (const std::string& key : airStateKeys()) CHECK(s.at(key).is_number());
        CHECK(s.at("T") == 298.15);
        CHECK(s.at("W") == 0.01);
        CHECK(s.at("RH").get<double>() > 0.0);
        CHECK(s.at("RH").get<double>() < 1.0);
        const json m = stateAt(air, "ThxW", 0.01, 300.0);
        REQUIRE(m.is_object());
        CHECK_THAT(m.at("Thx").get<double>(), WithinRel(300.0, 1e-9));
        CHECK(stateAt(air, "WT", 298.15, 0.05).is_null());  // above saturation: fog
        // The total pressure decides where saturation lies: 21 g/kg at 25 °C is fog at sea level, not at 0.9 bar
        CHECK(stateAt(air, "WT", 298.15, 0.021).is_null());
        CHECK(stateAt(chartFluid("AirH2O", "HEOS", {}, "", 9e4), "WT", 298.15, 0.021).is_object());
        CHECK_THROWS_WITH(stateAt(air, "ph", 3e5, 2e5), ContainsSubstring("cannot show humid air"));
        CHECK_THROWS_WITH(stateAt(chartFluid("R134a"), "WT", 300.0, 0.01), ContainsSubstring("cannot show R134a"));
        const auto start = std::chrono::steady_clock::now();
        const int n = 200;
        for (int i = 0; i < n; ++i) stateAt(air, "WT", 280.0 + 0.1 * i, 0.002 + 1e-5 * i);
        const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start).count() / n;
        INFO("state of humid air at a point: " << us << " us per call");
        CHECK(us < 5000.0);
    }
    SECTION("the default psychrometric and Mollier charts, cold and warm") {
        // Measured in the Release build of 2026-11-13 (Linux, GCC, 250 points per line, the lines on 4 threads): a cold
        // psychrometric chart 140 ms, a cold Mollier chart 130 ms, warm 2.4 and 3.7 ms. The budget is that of the
        // other charts (Diagram performance)
        for (const char* diagram : {"WT", "ThxW"}) {
            SceneRequest request;
            request.fluid = "AirH2O";
            request.spec = humidAirSpec(diagram);
            Engine::shared().clear();
            auto start = std::chrono::steady_clock::now();
            const json cold = scene(request);
            const double coldMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
            start = std::chrono::steady_clock::now();
            const json warm = scene(request);
            const double warmMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
            INFO(diagram << ": cold " << coldMs << " ms, warm " << warmMs << " ms, timing " << cold.at("timing_ms").dump());
            CHECK(coldMs < 1000.0);
            CHECK(warmMs < 100.0);
            CHECK(warm.at("svg") == cold.at("svg"));
        }
    }
}
