// Humid air for psychrometric charts (include/coolsolve/diagram/humid_air.h).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/humid_air.py (HumidAir,
// air_state, psychro_line, _sample_line, _W_by_search, _dry_end, _saturation_crossing, _start_on_saturation,
// _saturation_T, chart_ranges, _dry_bulb_span, saturation_W, _saturation_grid), core/limits.py (within) and
// core/state.py (AIR_STATE_KEYS and the constants of the Mollier chart).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#include "coolsolve/diagram/humid_air.h"

#include "HumidAirProp.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <tuple>

namespace coolsolve {
namespace diagram {

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
constexpr double INF = std::numeric_limits<double>::infinity();

/// Humidity ratio used for "dry air" in calls that need W > 0 (humid_air._W_DRY).
constexpr double W_DRY = 1e-12;
/// Steps a search for W may take before it gives up (humid_air._SEARCH_ITERATIONS).
constexpr int SEARCH_ITERATIONS = 60;
/// The specific heat of water vapour behind the dry-bulb span of a Mollier chart (humid_air._dry_bulb_span).
constexpr double CP_VAPOUR_J_KG_K = 1860.0;

/// The families whose range over the visible region a chart needs for round default values (humid_air._SPAN_KEYS).
const std::vector<std::string>& spanKeys() {
    static const std::vector<std::string> keys = {"RH", "Twb", "h", "v"};
    return keys;
}

/// The name of a quantity for HAPropsSI (core.quantities: humid_air_name). @throws DiagramError for another key
const char* haName(const std::string& key) {
    static const std::map<std::string, const char*> names = {{"T", "T"},   {"p", "P"},   {"W", "W"}, {"RH", "R"},
                                                             {"Twb", "B"}, {"Tdp", "D"}, {"h", "H"}, {"v", "V"}};
    const auto found = names.find(key);
    if (found == names.end())
        throw DiagramError("Unknown humid-air quantity '" + key + "', expected one of ['RH', 'T', 'Tdp', 'Twb', 'W', 'h', 'p', 'v'].");
    return found->second;
}

std::string quantityKey(const std::string& name) {
    try {
        return coolplot::quantity_key(name);
    } catch (const coolplot::Error& e) {
        throw DiagramError(e.what());
    }
}

/// limits.first <= value <= limits.second, allowing for rounding (-100 °C is 173.14999... K) (core.limits.within).
bool within(double value, double lo, double hi) {
    const double slack = 1e-9 * std::max(std::fabs(lo), std::fabs(hi));
    return lo - slack <= value && value <= hi + slack;
}

std::string pairText(const Range& range) {
    return "(" + coolplot::repr_number(range.first) + ", " + coolplot::repr_number(range.second) + ")";
}

/// The humidity ratio on a line of constant `key` at `T_K` (humid_air._W_on_line). @throws std::runtime_error
double wOnLine(const std::string& key, double value, double T_K, double p_Pa) {
    return humidAir("W", "T", T_K, key, value, p_Pa);
}

/// The saturation humidity ratios on the temperature grid of the lines, computed once per chart (humid_air.
/// _saturation_grid, an lru_cache of 32): the lines of one chart are computed on several threads.
std::shared_ptr<const std::vector<double>> saturationGrid(double T_lo, double T_hi, int points, double p_Pa) {
    using Key = std::tuple<double, double, int, double>;
    static std::mutex mutex;
    static std::map<Key, std::shared_ptr<const std::vector<double>>> kept;
    const Key key{T_lo, T_hi, points, p_Pa};
    {
        std::lock_guard<std::mutex> lock(mutex);
        const auto found = kept.find(key);
        if (found != kept.end()) return found->second;
    }
    auto grid = std::make_shared<std::vector<double>>();
    for (double T : linearGrid(T_lo, T_hi, points)) grid->push_back(saturationW(T, p_Pa));
    std::lock_guard<std::mutex> lock(mutex);
    if (kept.size() >= 32) kept.clear();
    return kept.emplace(key, grid).first->second;
}

/// How a search for W on a line ended (humid_air._W_by_search): a point, the line has ended below W = 0 ("dry"), the
/// point is above saturation ("fog"), or no point.
enum class Status { Ok, Dry, Fog, Failed };

/**
 * W on a line of Twb, h or v at T, from HAPropsSI(key, T, W) (humid_air._W_by_search). At fixed T, Twb, h and v all
 * grow with W, so the line value is bracketed between dry and saturated air, nearly linearly: a regula falsi search
 * (the "Illinois" variant, which halves the weight of an end that stays put; Dowell and Jarratt, BIT 11 (1971) 168)
 * needs a handful of fast calls. Used for wet-bulb lines, where asking HAPropsSI for W directly is slow and has gaps,
 * and as a fallback for the others.
 */
std::pair<double, Status> wBySearch(const std::string& key, double value, double T_K, double p_Pa, double W_sat) {
    auto excess = [&](double W) { return humidAir(key, "T", T_K, "W", W, p_Pa) - value; };
    try {
        double a = W_DRY, f_a = excess(W_DRY);
        if (f_a > 0.0) return {NaN, Status::Dry};
        double b = std::isfinite(W_sat) ? W_sat : 10.0;
        double f_b = excess(b);
        if (f_b < 0.0) return {NaN, Status::Fog};
        const double tolerance = 1e-7 * std::max(std::fabs(value), 1.0);
        int side = 0;
        for (int k = 0; k < SEARCH_ITERATIONS; ++k) {
            // A regula falsi step, or the midpoint when the residuals are the same (both zero: f_a <= 0 <= f_b holds)
            const double denominator = f_b - f_a;
            const double c = !coolplot::is_zero(denominator) ? (a * f_b - b * f_a) / denominator : 0.5 * (a + b);
            const double f_c = excess(c);
            if (std::fabs(f_c) <= tolerance || b - a < 1e-9) return {c, Status::Ok};
            if (f_c > 0.0) {
                b = c;
                f_b = f_c;
                if (side == +1) f_a *= 0.5;
                side = +1;
            } else {
                a = c;
                f_a = f_c;
                if (side == -1) f_b *= 0.5;
                side = -1;
            }
        }
        return {NaN, Status::Failed};  // no convergence: not a point to draw
    } catch (const std::exception&) {
        return {NaN, Status::Failed};
    }
}

/// The temperature in [T_lo, T_hi] where a falling line reaches W = 0, by bisection (humid_air._dry_end).
std::optional<double> dryEnd(const std::string& key, double value, double T_lo, double T_hi, double p_Pa) {
    auto excess = [&](double T) { return humidAir(key, "T", T, "W", W_DRY, p_Pa) - value; };
    try {
        if (excess(T_lo) > 0.0 || excess(T_hi) < 0.0) return std::nullopt;
        for (int k = 0; k < 40; ++k) {
            const double T_mid = 0.5 * (T_lo + T_hi);
            if (excess(T_mid) > 0.0) T_hi = T_mid;
            else T_lo = T_mid;
            if (T_hi - T_lo < 1e-3) break;
        }
        return 0.5 * (T_lo + T_hi);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

/// The temperature in [T_lo, T_hi] where the line meets saturation, by bisection (humid_air._saturation_crossing).
std::optional<double> saturationCrossing(const std::string& key, double value, double T_lo, double T_hi, double p_Pa) {
    auto excess = [&](double T) {
        try {
            return wOnLine(key, value, T, p_Pa) - saturationW(T, p_Pa);
        } catch (const std::exception&) {
            return NaN;
        }
    };
    double f_lo = excess(T_lo);
    const double f_hi = excess(T_hi);
    if (!(std::isfinite(f_lo) && std::isfinite(f_hi)) || f_lo * f_hi > 0.0) return std::nullopt;
    for (int k = 0; k < 50; ++k) {
        const double T_mid = 0.5 * (T_lo + T_hi);
        const double f_mid = excess(T_mid);
        if (!std::isfinite(f_mid)) return std::nullopt;
        if (f_mid * f_lo > 0.0) {
            T_lo = T_mid;
            f_lo = f_mid;
        } else {
            T_hi = T_mid;
        }
        if (T_hi - T_lo < 1e-6) break;
    }
    return 0.5 * (T_lo + T_hi);
}

/// A line sampled at the temperatures T (humid_air._sample_line): W (NaN where there is no point), the points
/// HAPropsSI failed for, and the temperature where a falling line reaches W = 0 (the sampling stops there).
struct Sampled {
    std::vector<double> W;
    int n_failed = 0;
    std::optional<double> dry_end;
};

Sampled sampleLine(const std::string& key, double value, const std::vector<double>& T, const std::vector<double>& W_sat, double p_Pa) {
    Sampled out;
    out.W.assign(T.size(), NaN);
    std::vector<double>& W = out.W;
    const bool falling = key != "RH";  // W falls with T along lines of Twb, h and v
    for (size_t i = 0; i < T.size(); ++i) {
        const double t = T[i];
        if (key == "Twb" && t < value) continue;  // a wet-bulb line starts on saturation, at T = Twb
        Status status = Status::Ok;
        if (key == "Twb") {
            std::tie(W[i], status) = wBySearch(key, value, t, p_Pa, W_sat[i]);
            if (status == Status::Failed) ++out.n_failed;
        } else {
            try {
                W[i] = wOnLine(key, value, t, p_Pa);
                status = falling && W[i] < 0.0 ? Status::Dry : Status::Ok;
            } catch (const std::exception&) {
                if (!falling) {
                    if (std::any_of(W.begin(), W.begin() + static_cast<std::ptrdiff_t>(i), [](double w) { return std::isfinite(w); }))
                        break;  // above boiling: the vapour would exceed the total pressure, the line ends
                    ++out.n_failed;
                    continue;
                }
                std::tie(W[i], status) = wBySearch(key, value, t, p_Pa, W_sat[i]);
                if (status == Status::Failed) ++out.n_failed;
            }
        }
        if (status == Status::Dry) {
            // The line has reached W = 0: the exact end is added and the sampling stops, the rest would only fail
            W[i] = NaN;
            if (i > 0 && std::isfinite(W[i - 1])) out.dry_end = dryEnd(key, value, T[i - 1], t, p_Pa);
            break;
        }
    }
    return out;
}

/// The temperature where the line meets saturation, between the last fog point and the first inside one
/// (humid_air._start_on_saturation).
std::optional<double> startOnSaturation(const std::string& key, double value, const std::vector<double>& T,
                                        const std::vector<bool>& inside, double p_Pa) {
    if (key == "Twb" && T.front() <= value && value <= T.back()) return value;  // on saturation at T = Twb by definition
    const size_t first = static_cast<size_t>(std::find(inside.begin(), inside.end(), true) - inside.begin());
    if (first > 0 && first < T.size()) return saturationCrossing(key, value, T[first - 1], T[first], p_Pa);
    return std::nullopt;
}

/// The temperature in [T_lo, T_hi] where saturated air has humidity ratio W, by bisection (humid_air._saturation_T).
std::optional<double> saturationT(double W, double T_lo, double T_hi, double p_Pa) {
    auto excess = [&](double T) {
        const double W_sat = saturationW(T, p_Pa);
        return (std::isfinite(W_sat) ? W_sat : INF) - W;
    };
    if (!(excess(T_lo) < 0.0 && 0.0 < excess(T_hi))) return std::nullopt;
    for (int k = 0; k < 40; ++k) {
        const double T_mid = 0.5 * (T_lo + T_hi);
        if (excess(T_mid) > 0.0) T_hi = T_mid;
        else T_lo = T_mid;
        if (T_hi - T_lo < 1e-3) break;
    }
    return 0.5 * (T_lo + T_hi);
}

/// The dry-bulb temperatures behind a Mollier chart's window, with a margin (humid_air._dry_bulb_span). Along a
/// horizontal of the Mollier chart the dry-bulb temperature falls with W (in °C: t = thx cp_a / (cp_a + cp_v W)), so
/// the corners of the window bound it; 5 K on either side cover the rest, within the limits of HAPropsSI.
Range dryBulbSpan(const Range& Thx_K, const Range& W) {
    double lo = INF, hi = -INF;
    for (double Thx : {Thx_K.first, Thx_K.second})
        for (double w : {W.first, W.second}) {
            const double corner = 273.15 + (Thx - 273.15) * CP_DRY_AIR_J_KG_K / (CP_DRY_AIR_J_KG_K + CP_VAPOUR_J_KG_K * w);
            lo = std::min(lo, corner);
            hi = std::max(hi, corner);
        }
    return {std::max(lo - 5.0, HUMID_AIR_T_MIN_K), std::min(hi + 5.0, HUMID_AIR_T_MAX_K)};
}

}  // namespace

// ---------------------------------------------------------------------------
// States (thermo/humid_air.py, core/state.py)
// ---------------------------------------------------------------------------

const std::vector<std::string>& airStateKeys() {
    static const std::vector<std::string> keys = {"T", "W", "RH", "Twb", "Tdp", "h", "v", "Thx", "p"};
    return keys;
}

const std::vector<std::string>& airInputKeys() {
    static const std::vector<std::string> keys = {"T", "W", "RH", "Twb", "Tdp", "h", "v", "Thx"};
    return keys;
}

double checkedHumidAirPressure(double p_Pa) {
    if (!within(p_Pa, HUMID_AIR_P_MIN_PA, HUMID_AIR_P_MAX_PA))
        throw DiagramError("The total pressure has to lie between " + coolplot::repr_number(HUMID_AIR_P_MIN_PA) + " and " +
                           coolplot::repr_number(HUMID_AIR_P_MAX_PA) + " Pa for HAPropsSI, not " + coolplot::repr_number(p_Pa) + " Pa.");
    return p_Pa;
}

std::string humidAirName(double p_Pa) {
    return "humid air at " + coolplot::general_text(p_Pa, 6) + " Pa";
}

double mollierTemperature(double h_J_kg, double W) {
    return 273.15 + (h_J_kg - R0_J_KG * W) / CP_DRY_AIR_J_KG_K;
}

double enthalpyFromMollier(double Thx_K, double W) {
    return CP_DRY_AIR_J_KG_K * (Thx_K - 273.15) + R0_J_KG * W;
}

double humidAir(const std::string& output, const std::string& key1, double value1, const std::string& key2, double value2,
                double p_Pa) {
    const char* out = haName(output);
    const char* in1 = haName(key1);
    const char* in2 = haName(key2);
    const double value = HumidAir::HAPropsSI(out, in1, value1, "P", p_Pa, in2, value2);
    if (!std::isfinite(value))
        throw std::runtime_error(std::string("HAPropsSI has no ") + out + " at " + in1 + "=" + coolplot::general_text(value1, 6) +
                                 ", " + in2 + "=" + coolplot::general_text(value2, 6) + ", P=" + coolplot::general_text(p_Pa, 6));
    return value;
}

double saturationW(double T_K, double p_Pa) {
    try {
        return humidAir("W", "T", T_K, "RH", 1.0, p_Pa);
    } catch (const std::exception&) {
        return NaN;
    }
}

State airState(const std::string& given1, double value1, const std::string& given2, double value2, double p_Pa) {
    const std::string key1 = quantityKey(given1), key2 = quantityKey(given2);
    if (key1 == key2) throw DiagramError("Give two different properties, not '" + key1 + "' twice.");
    if (key1 == "Thx" || key2 == "Thx") {
        // A point of the Mollier chart: its temperature scale and the humidity ratio give the enthalpy
        const bool first = key1 == "Thx";
        const std::string& other = first ? key2 : key1;
        const double Thx = first ? value1 : value2, W = first ? value2 : value1;
        if (other != "W")
            throw DiagramError("The temperature scale Thx of the Mollier chart is given together with the humidity ratio W, not with '" +
                               other + "'.");
        return airState("h", enthalpyFromMollier(Thx, W), "W", W, p_Pa);
    }
    // The total pressure is the chart's, not one of the two properties: HAPropsSI would return saturated air
    if (key1 == "p" || key2 == "p")
        throw DiagramError("Humid air is at the total pressure of the chart (" + coolplot::general_text(p_Pa, 6) +
                           " Pa); give two of T, W, RH, Twb, Tdp, h and v (or Thx and W).");
    State state;
    for (const std::string& key : airInputKeys())
        if (key != "Thx") state[key] = humidAir(key, key1, value1, key2, value2, p_Pa);
    state["Thx"] = mollierTemperature(state.at("h"), state.at("W"));
    state["p"] = p_Pa;
    return state;
}

// ---------------------------------------------------------------------------
// Ranges and lines (thermo/humid_air.py)
// ---------------------------------------------------------------------------

std::vector<std::pair<std::string, Range>> chartRanges(double p_Pa, const Range& T_K, const Range& W) {
    std::vector<std::pair<double, double>> corners;  // (T, W)
    for (double T : {T_K.first, T_K.second}) {
        const double W_sat = saturationW(T, p_Pa);
        const double W_top = std::isfinite(W_sat) ? std::min(W.second, W_sat) : W.second;  // no saturation above boiling
        if (W_top >= W.first) {  // otherwise this edge of the chart is all fog
            corners.push_back({T, W.first});
            corners.push_back({T, W_top});
        }
    }
    // Where saturation crosses the bottom and top edges: the corners of the visible region when an edge is fog
    for (double w : {W.first, W.second})
        if (const std::optional<double> T_sat = saturationT(w, T_K.first, T_K.second, p_Pa)) corners.push_back({*T_sat, w});
    std::vector<std::pair<std::string, Range>> out = {{"RH", {0.0, 1.0}}};
    for (const char* key : {"Twb", "h", "v"}) {
        double lo = INF, hi = -INF;
        bool any = false;
        for (const auto& [T, w] : corners) {
            try {
                const double value = humidAir(key, "T", T, "W", std::max(w, 1e-9), p_Pa);
                lo = std::min(lo, value);
                hi = std::max(hi, value);
                any = true;
            } catch (const std::exception&) {
            }
        }
        if (any) out.push_back({key, {lo, hi}});
    }
    std::sort(out.begin(), out.end());  // sorted by key, as CoolPlot lists them
    return out;
}

Curve psychroLine(const std::string& key, double value, double p_Pa, const Range& T_K, int points) {
    if (std::find(spanKeys().begin(), spanKeys().end(), key) == spanKeys().end())  // isotherms: the Mollier chart's
        throw DiagramError("Lines of constant '" + key + "' are not offered, expected one of ('RH', 'Twb', 'h', 'v').");
    const std::vector<double> T = linearGrid(T_K.first, T_K.second, points);
    const std::shared_ptr<const std::vector<double>> W_sat = saturationGrid(T_K.first, T_K.second, points, p_Pa);
    const Sampled sampled = sampleLine(key, value, T, *W_sat, p_Pa);
    const std::vector<double>& W = sampled.W;
    std::vector<bool> inside(T.size());
    for (size_t i = 0; i < T.size(); ++i) {
        const double limit = std::isfinite((*W_sat)[i]) ? (*W_sat)[i] * (1.0 + 1e-9) : INF;  // no saturation above boiling
        inside[i] = key == "RH" && value >= 1.0 ? std::isfinite(W[i])  // the saturation line itself
                                                 : std::isfinite(W[i]) && W[i] >= 0.0 && W[i] <= limit;
    }
    Curve curve;
    for (size_t i = 0; i < T.size(); ++i)
        if (inside[i]) {
            curve.x.push_back(T[i]);
            curve.y.push_back(W[i]);
        }
    if (key != "RH" && !curve.x.empty()) {
        const std::optional<double> crossing = startOnSaturation(key, value, T, inside, p_Pa);
        if (crossing && *crossing < curve.x.front()) {
            curve.x.insert(curve.x.begin(), *crossing);
            curve.y.insert(curve.y.begin(), saturationW(*crossing, p_Pa));
        }
    }
    if (sampled.dry_end && !curve.x.empty() && *sampled.dry_end > curve.x.back()) {
        curve.x.push_back(*sampled.dry_end);
        curve.y.push_back(0.0);
    }
    curve.n_failed = sampled.n_failed;
    return curve;
}

// ---------------------------------------------------------------------------
// HumidAirModel (thermo/humid_air.py: HumidAir)
// ---------------------------------------------------------------------------

HumidAirModel::HumidAirModel(double p_Pa) : p_Pa_(checkedHumidAirPressure(p_Pa)), name_(humidAirName(p_Pa)) {}

std::string HumidAirModel::pathKind(const State& /*a*/, const State& /*b*/, const std::string& /*kind*/) const {
    return "straight";
}

std::string HumidAirModel::stateProblem(const State& state) const {
    const auto p = state.find("p");
    // The same tolerance as CoolPlot's chart.check_state
    if (p != state.end() && std::fabs(p->second - p_Pa_) > 1e-6 * p_Pa_)
        return "it is at " + coolplot::general_text(p->second, 6) + " Pa, but the chart shows " + name_ +
               " (set the total pressure of the chart)";
    return "";
}

void HumidAirModel::update(const std::string& key1, double value1, const std::string& key2, double value2) {
    state_.clear();
    state_ = airState(key1, value1, key2, value2, p_Pa_);
}

double HumidAirModel::get(const std::string& key) const {
    const std::string k = quantityKey(key);
    const auto found = state_.find(k);
    if (found != state_.end()) return found->second;
    if (std::find(airStateKeys().begin(), airStateKeys().end(), k) == airStateKeys().end())
        throw DiagramError("'" + key + "' is not a quantity of a state of humid air.");
    throw std::runtime_error(name_ + " has no state yet.");
}

double HumidAirModel::tTripleK() { return NaN; }

Domain HumidAirModel::domain(const Region& region) const {
    const Range W = region.rangeSI("W");
    const bool mollier = region.has("Thx");
    // The bounds as given are checked: the dry-bulb span of a Mollier chart is clamped to the limits
    const Range given = region.rangeSI(mollier ? "Thx" : "T");
    for (double T : {given.first, given.second})
        if (!within(T, HUMID_AIR_T_MIN_K, HUMID_AIR_T_MAX_K))
            throw DiagramError("Temperatures have to lie between " + coolplot::repr_number(HUMID_AIR_T_MIN_K) + " and " +
                               coolplot::repr_number(HUMID_AIR_T_MAX_K) + " K for HAPropsSI, got " + pairText(given) + " K.");
    const Range T = mollier ? dryBulbSpan(given, W) : given;
    if (W.first < 0.0) throw DiagramError("Humidity ratios have to be at least 0, got " + pairText(W) + ".");
    Domain domain;
    domain.tp_box = {T.first, T.second, p_Pa_, p_Pa_};
    domain.ranges_SI = {{"T", T}, {"W", W}, {"p", {p_Pa_, p_Pa_}}};
    if (mollier) domain.ranges_SI.push_back({"Thx", given});
    const std::vector<std::pair<std::string, Range>> spans = chartRanges(p_Pa_, T, W);
    if (spans.size() < spanKeys().size())
        throw DiagramError("The region T = " + pairText(T) + " K, W = " + pairText(W) + " lies entirely above saturation (fog) at " +
                           coolplot::general_text(p_Pa_, 6) + " Pa; there are no humid-air states to draw.");
    domain.ranges_SI.insert(domain.ranges_SI.end(), spans.begin(), spans.end());
    domain.saturation = true;  // the line of 100 % relative humidity (Engine.domain)
    return domain;
}

Curve HumidAirModel::line(const std::string& key, double value_SI, const Domain& domain, int points,
                          const coolplot::DiagramType& diagram) const {
    const Range T = {domain.tp_box[0], domain.tp_box[1]};
    const bool mollier = diagram.kind == "humid_air" && diagram.y == "Thx";
    if (!mollier) return psychroLine(key, value_SI, p_Pa_, T, points);
    if (key == "T") {
        std::optional<Range> W;
        for (const auto& [k, range] : domain.ranges_SI)
            if (k == "W") W = range;
        return isotherm(value_SI, points, W);
    }
    return toMollier(psychroLine(key, value_SI, p_Pa_, T, points));
}

Curve HumidAirModel::isotherm(double T_K, int points, const std::optional<Range>& W) const {
    double W_end = saturationW(T_K, p_Pa_);
    if (W) {
        const double W_cap = 1.2 * W->second;  // no need to compute far beyond the top of the chart
        W_end = std::isfinite(W_end) ? std::min(W_end, W_cap) : W_cap;
    }
    Curve curve;
    if (!(std::isfinite(W_end) && W_end > 0.0)) {
        char text[64];
        std::snprintf(text, sizeof(text), "%.2f", T_K);
        curve.reason = std::string("no humid-air states at ") + text + " K";
        return curve;
    }
    curve.y = linearGrid(0.0, W_end, std::max(points, 2));
    curve.x.assign(curve.y.size(), T_K);
    return toMollier(std::move(curve));
}

Curve HumidAirModel::toMollier(Curve curve) const {
    // The curve in (T, W) becomes (W, Thx); points HAPropsSI cannot compute become gaps (NaN)
    std::vector<double> Thx(curve.y.size(), NaN);
    for (size_t i = 0; i < curve.y.size(); ++i) {
        const double t = curve.x[i], w = curve.y[i];
        if (!(std::isfinite(t) && std::isfinite(w))) continue;
        try {
            Thx[i] = mollierTemperature(humidAir("h", "T", t, "W", std::max(w, 0.0), p_Pa_), w);
        } catch (const std::exception&) {
        }
    }
    curve.x = std::move(curve.y);
    curve.y = std::move(Thx);
    return curve;
}

}  // namespace diagram
}  // namespace coolsolve
