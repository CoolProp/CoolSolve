// Computing what a plan asks for, keeping the results, and turning a chart request into a scene
// (include/coolsolve/diagram/engine.h).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): engine/engine.py (Engine: domain,
// curves, warnings), engine/cache.py (CurveCache), thermo/native.py (the plots of the port), corelib.py
// (compose_record), chart.py (the plot area measured again after drawing, Chart.draw) and render/pdf.py (the PDF
// writer's plot area).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#include "coolsolve/diagram/engine.h"

#include "coolsolve/diagram/humid_air.h"
#include "coolsolve/diagram/ideal_gas.h"
#include "coolsolve/diagram/references.h"

#include "CoolProp/Configuration.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/plotting/CoolPropPlot.h"

#include <coolplot/capi.h>
#include <coolplot/compose.h>
#include <coolplot/pdf.h>
#include <coolplot/svg.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <list>
#include <set>
#include <thread>
#include <unordered_map>

#ifndef COOLSOLVE_COOLPLOT_COMMIT
#define COOLSOLVE_COOLPLOT_COMMIT "unknown"
#endif

namespace coolsolve {
namespace diagram {

namespace {

// ---------------------------------------------------------------------------
// A bounded, thread-safe store (engine/cache.py: CurveCache keeps curves this way)
// ---------------------------------------------------------------------------

/// Values by key; the least recently used one is dropped when full. Every access holds the lock; values are
/// computed outside of it, so two threads may compute the same value: the first one stored is kept.
template <class V>
class LruCache {
public:
    explicit LruCache(size_t capacity) : capacity_(capacity) {}

    std::optional<V> get(const std::string& key) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto found = index_.find(key);
        if (found == index_.end()) return std::nullopt;
        items_.splice(items_.begin(), items_, found->second);
        return found->second->second;
    }

    /// Store `value` unless the key has one; return the value stored.
    V insert(const std::string& key, V value) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto found = index_.find(key);
        if (found != index_.end()) return found->second->second;
        items_.emplace_front(key, std::move(value));
        index_[key] = items_.begin();
        if (items_.size() > capacity_) {
            index_.erase(items_.back().first);
            items_.pop_back();
        }
        return items_.front().second;
    }

    /// The stored value, or `compute()` stored first.
    template <class Compute>
    V getOrCompute(const std::string& key, Compute&& compute) {
        if (std::optional<V> found = get(key)) return *found;
        return insert(key, compute());
    }

    void clear() {
        std::lock_guard<std::mutex> lock(mutex_);
        items_.clear();
        index_.clear();
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return items_.size();
    }

private:
    mutable std::mutex mutex_;
    size_t capacity_;
    std::list<std::pair<std::string, V>> items_;
    std::unordered_map<std::string, typename std::list<std::pair<std::string, V>>::iterator> index_;
};

/// A number in a cache key: every digit, so that two keys differ when the values do.
std::string keyNumber(double value) {
    return coolplot::general_text(value, 17);
}

std::string keyBox(const std::array<double, 4>& box) {
    return keyNumber(box[0]) + "," + keyNumber(box[1]) + "," + keyNumber(box[2]) + "," + keyNumber(box[3]);
}

std::string keyUnit(const coolplot::Unit& unit) {
    return unit.label + "," + keyNumber(unit.scale) + "," + keyNumber(unit.offset);
}

/// What decides CoolProp's answers beyond the fluid: its revision and its configuration (thermo.fluid.Fluid.fingerprint).
std::string coolpropFingerprint() {
    static const std::string revision = CoolProp::get_global_param_string("gitrevision");
    std::string config;
    try {
        config = CoolProp::get_config_as_json_string();
    } catch (const std::exception&) {
    }
    return revision + "/" + std::to_string(std::hash<std::string>{}(config));
}

double elapsedMs(std::chrono::steady_clock::time_point since) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - since).count();
}

/// An optional finite number as JSON (null otherwise), as CoolPlot's records write NaN.
json finiteOrNull(double value) {
    return std::isfinite(value) ? json(value) : json();
}

/// The quantities of the charts: those of a fluid state and those of humid air.
const std::vector<std::string>& chartKeys() {
    static const std::vector<std::string> keys = [] {
        std::vector<std::string> out = stateKeys();
        for (const std::string& key : airStateKeys())
            if (std::find(out.begin(), out.end(), key) == out.end()) out.push_back(key);
        return out;
    }();
    return keys;
}

/// The units of a unit system as {key: [label, scale, offset(, difference)]} for the quantities of the charts.
json unitTable(const coolplot::UnitSystem& units) {
    json table = json::object();
    for (const std::string& key : chartKeys()) {
        const coolplot::Unit& unit = units.at(key);
        json entry = json::array({unit.label, unit.scale, unit.offset});
        if (!unit.difference_label.empty()) entry.push_back(unit.difference_label);
        table[key] = entry;
    }
    return table;
}

}  // namespace

// ---------------------------------------------------------------------------
// Curve (thermo/isolines.py: Curve)
// ---------------------------------------------------------------------------

size_t Curve::nValid() const {
    size_t n = 0;
    for (size_t i = 0; i < x.size(); ++i)
        if (std::isfinite(x[i]) && std::isfinite(y[i])) ++n;
    return n;
}

size_t Curve::nGaps() const {
    size_t first = x.size(), last = 0;
    for (size_t i = 0; i < x.size(); ++i)
        if (std::isfinite(x[i]) && std::isfinite(y[i])) {
            first = std::min(first, i);
            last = i;
        }
    if (first >= x.size()) return 0;
    size_t gaps = 0;
    for (size_t i = first; i <= last; ++i)
        if (!(std::isfinite(x[i]) && std::isfinite(y[i]))) ++gaps;
    return gaps;
}

json Curve::toJson() const {
    json xs = json::array(), ys = json::array();
    for (double v : x) xs.push_back(finiteOrNull(v));
    for (double v : y) ys.push_back(finiteOrNull(v));
    return {{"x", xs},
            {"y", ys},
            {"n_failed", n_failed},
            {"n_fallback", n_fallback},
            {"n_inserted", n_inserted},
            {"method", method},
            {"reason", reason},
            {"n_evaluated", n_evaluated},
            {"n_outside", n_outside}};
}

Chart::Chart(ChartFluid fluid_, PlanSpec spec_) : fluid(std::move(fluid_)), spec(std::move(spec_)), fingerprint(coolpropFingerprint()) {}

// ---------------------------------------------------------------------------
// Engine (engine/engine.py, engine/cache.py, thermo/native.py)
// ---------------------------------------------------------------------------

namespace {

/// Whether a diagram type draws the fluid: a chart of humid air shows humid air, a chart of a fluid a fluid (the
/// check of DiagramSpec). @throws DiagramError otherwise
void checkChartKind(const ChartFluid& fluid, const coolplot::DiagramType& diagram) {
    if ((diagram.kind == "humid_air") != fluid.isHumidAir())
        throw DiagramError("A " + diagram.name() + " diagram cannot show " + fluid.name() + "." +
                           (fluid.isHumidAir() ? " Humid air is drawn in the psychrometric (WT) and Mollier (ThxW) charts."
                                               : " The psychrometric and Mollier charts are those of humid air (AirH2O)."));
}

/// The property provider of a chart fluid (provider.make_provider): HAPropsSI's for humid air, CoolSolve's for an ideal
/// gas, CoolProp's otherwise.
std::unique_ptr<PropertyProvider> makeProvider(const ChartFluid& fluid) {
    if (fluid.isHumidAir()) return std::make_unique<HumidAirModel>(fluid.p_Pa);
    if (fluid.isIdealGas()) return std::make_unique<IdealGasModel>(fluid);
    return std::make_unique<FluidModel>(fluid);
}

}  // namespace

struct Engine::FluidEntry {
    explicit FluidEntry(const ChartFluid& fluid) : model(makeProvider(fluid)) {}
    std::mutex mutex;
    std::unique_ptr<PropertyProvider> model;
};

/// Values computed with a warning that is given again whenever they are used.
struct WarnedPairs {
    std::vector<Range> pairs;
    std::vector<std::string> warnings;
};

struct Engine::Impl {
    LruCache<std::shared_ptr<FluidEntry>> fluids{16};
    LruCache<Domain> domains{64};
    LruCache<std::shared_ptr<const WarnedPairs>> saturationValues{256};  // anchored families and scale ticks
    LruCache<std::vector<std::pair<std::string, State>>> markedPoints{64};
    LruCache<std::vector<DomeEnd>> domes{64};
    LruCache<std::shared_ptr<CoolProp::Plot::PropertyPlot>> plots{32};  // thermo.native._PLOTS_KEPT is 8 per thread
    LruCache<std::shared_ptr<const Curve>> curves{2000};               // CurveCache.max_size
};

Engine& Engine::shared() {
    static Engine engine;
    return engine;
}

Engine::Engine() : impl_(std::make_unique<Impl>()) {}
Engine::~Engine() = default;

std::shared_ptr<Engine::FluidEntry> Engine::fluidEntry(const ChartFluid& fluid) {
    const std::string key = fluid.fluidString() + "|" + coolpropFingerprint();
    return impl_->fluids.getOrCompute(key, [&fluid] {
        try {
            return std::make_shared<FluidEntry>(fluid);
        } catch (const std::exception& e) {
            throw DiagramError("CoolProp cannot load " + fluid.fluidString() + ": " + e.what());
        }
    });
}

PropertyProvider& Engine::model(FluidEntry& entry) { return *entry.model; }
std::mutex& Engine::mutex(FluidEntry& entry) { return entry.mutex; }

void Engine::clear() {
    impl_->fluids.clear();
    impl_->domains.clear();
    impl_->saturationValues.clear();
    impl_->markedPoints.clear();
    impl_->domes.clear();
    impl_->plots.clear();
    impl_->curves.clear();
}

size_t Engine::curveCount() const { return impl_->curves.size(); }

Domain Engine::domain(const Chart& chart, std::vector<std::string>& warnings) {
    const PlanSpec& spec = chart.spec;
    const std::string fluidKey = chart.fluid.fluidString() + "|" + chart.fingerprint;
    const double margin = spec.calc.critical_margin;
    const std::string region = spec.limitsJson().dump() + "|" + keyNumber(margin);
    checkChartKind(chart.fluid, spec.diagram);
    try {
        if (chart.fluid.isHumidAir()) {
            // A region of T (or Thx) and W at one total pressure, without saturation states (HumidAir.domain); the
            // labels inside the dome are those of a fluid's chart (Engine._dome_ends)
            return impl_->domains.getOrCompute(fluidKey + "|domain|" + region,
                                               [&] { return HumidAirModel(chart.fluid.p_Pa).domain(spec.region); });
        }
        Domain domain = impl_->domains.getOrCompute(fluidKey + "|domain|" + region, [&] {
            return withFluid(chart.fluid, [&](PropertyProvider& fluid) { return chartDomain(fluid, spec.limits, margin); });
        });
        // The families placed by saturation states
        for (const IsolineSet& iso : spec.isolines) {
            if (!iso.anchor) continue;
            const SaturationAnchor& anchor = *iso.anchor;
            const coolplot::Unit& unit = spec.units.at(anchor.by);  // the round steps are round in the units shown
            const std::string key = fluidKey + "|anchor|" + region + "|" + iso.key + "|" + anchor.by + "|" +
                                    keyNumber(anchor.step_SI) + "|" + keyNumber(anchor.Q) + "|" +
                                    (anchor.range_SI ? keyNumber(anchor.range_SI->first) + "," + keyNumber(anchor.range_SI->second) : "") +
                                    "|" + keyUnit(unit);
            auto values = impl_->saturationValues.getOrCompute(key, [&] {
                auto computed = std::make_shared<WarnedPairs>();
                withFluid(chart.fluid, [&](PropertyProvider& fluid) {
                    computed->pairs = saturationAnchorValues(fluid, iso.key, anchor.by, anchor.step_SI, anchor.Q, anchor.range_SI,
                                                             domain.tp_box, margin, unit, computed->warnings);
                });
                return std::shared_ptr<const WarnedPairs>(computed);
            });
            domain.anchored.push_back({iso.key, values->pairs});
            warnings.insert(warnings.end(), values->warnings.begin(), values->warnings.end());
        }
        // The marked critical and triple points
        if (!spec.marked_points.empty()) {
            std::string names;
            for (const std::string& name : spec.marked_points) names += name + ",";
            domain.points = impl_->markedPoints.getOrCompute(fluidKey + "|marked|" + keyNumber(margin) + "|" + names, [&] {
                return withFluid(chart.fluid, [&](PropertyProvider& fluid) { return markedPointStates(fluid, spec.marked_points, margin); });
            });
        }
        // The ticks of the saturation scale (Engine._scale_ticks)
        if (spec.saturation_scale) {
            const SaturationScale& scale = *spec.saturation_scale;
            const std::string by = scaleKey(spec.diagram), key = spec.diagram.y;
            const coolplot::Unit& unit = spec.units.at(by);  // the ticks are round in the units shown
            const std::string cacheKey = fluidKey + "|scale|" + region + "|" + (scale.step_SI ? keyNumber(*scale.step_SI) : "") +
                                         "|" + keyNumber(scale.Q) + "|" + by + "|" + key + "|" + keyUnit(unit);
            auto ticks = impl_->saturationValues.getOrCompute(cacheKey, [&] {
                auto computed = std::make_shared<WarnedPairs>();
                withFluid(chart.fluid, [&](PropertyProvider& fluid) {
                    if (!fluid.hasCriticalPoint()) return;  // no saturation curve, no scale
                    if (scale.step_SI) {
                        computed->pairs = saturationAnchorValues(fluid, key, by, *scale.step_SI, scale.Q, std::nullopt, domain.tp_box,
                                                                 margin, unit, computed->warnings);
                        return;
                    }
                    // Round values over the saturation states of the chart: 1-2-5 per decade for pressures
                    auto [lo, hi] = fluid.saturationRange(by, margin);
                    lo = std::max(lo, by == "T" ? domain.tp_box[0] : domain.tp_box[2]);
                    hi = std::min(hi, by == "T" ? domain.tp_box[1] : domain.tp_box[3]);
                    if (!(hi > lo)) return;
                    std::vector<double> values;
                    for (double v : roundValues(unit, lo, hi, SCALE_TICKS, by == "p"))
                        if (lo <= v && v <= hi) values.push_back(v);
                    computed->pairs = saturationPairs(fluid, by, key, values, scale.Q, computed->warnings);
                });
                return std::shared_ptr<const WarnedPairs>(computed);
            });
            domain.scale = ticks->pairs;
            warnings.insert(warnings.end(), ticks->warnings.begin(), ticks->warnings.end());
        }
        // Where the labelled lines of the dome family meet the bubble and dew lines (Engine._dome_ends)
        const std::string dome = domeKey(spec.diagram);
        if (spec.dome_labels && !dome.empty()) {
            std::vector<double> values;
            for (const CurveRequest& request : plan(spec, domain))
                if (request.key == dome && request.labels) values.push_back(request.value_SI);
            if (!values.empty()) {
                std::string key = fluidKey + "|dome|" + keyNumber(margin) + "|" + dome + "|";
                for (double v : values) key += keyNumber(v) + ",";
                domain.dome = impl_->domes.getOrCompute(key, [&] {
                    return withFluid(chart.fluid, [&](PropertyProvider& fluid) { return domeEnds(fluid, dome, values, margin); });
                });
            }
        }
        return domain;
    } catch (const DiagramError&) {
        throw;
    } catch (const std::exception& e) {
        throw DiagramError("The region of the chart of " + chart.fluid.name() + " cannot be computed: " + e.what());
    }
}

namespace {

/// Everything that determines one curve (Engine.cache_key): lines of constant quality do not depend on the T-p box
/// but on the saturation settings; all other lines are swept through the box.
std::string curveKey(const Chart& chart, const Domain& domain, const CurveRequest& request) {
    const CalcSettings& calc = chart.spec.calc;
    std::string key = chart.fluid.portFluidString(request.key) + "|" + chart.fingerprint + "|" + chart.spec.diagram.name() + "|" +
                      request.key + "|" + keyNumber(request.value_SI) + "|" + std::to_string(request.points) + "|";
    if (request.key == "Q") key += keyNumber(calc.critical_margin) + "|" + calc.saturation_spacing;
    else key += keyBox(domain.tp_box) + "|" + calc.solver + "|" + keyNumber(calc.critical_margin);
    // An isotherm of the Mollier chart ends a little beyond the top of the humidity ratio range, which the T-p box
    // does not hold
    if (request.key == "T" && chart.fluid.isHumidAir())
        for (const auto& [k, range] : domain.ranges_SI)
            if (k == "W") key += "|" + keyNumber(range.first) + "," + keyNumber(range.second);
    return key;
}

}  // namespace

std::shared_ptr<const Curve> Engine::curve(const Chart& chart, const Domain& domain, const CurveRequest& request, size_t slot) {
    const PlanSpec& spec = chart.spec;
    const CalcSettings& calc = spec.calc;
    const std::string key = curveKey(chart, domain, request);
    if (std::optional<std::shared_ptr<const Curve>> found = impl_->curves.get(key)) return *found;

    if (chart.fluid.isHumidAir()) {  // HAPropsSI, per thread in CoolProp: a provider of its own (cheap, per call)
        const HumidAirModel air(chart.fluid.p_Pa);
        auto computed = std::make_shared<Curve>(air.line(request.key, request.value_SI, domain, request.points, spec.diagram));
        return impl_->curves.insert(key, std::shared_ptr<const Curve>(computed));
    }
    if (chart.fluid.isIdealGas()) {  // CoolSolve's provider: no plot, a model of its own (cheap, per call)
        IdealGasModel gas(chart.fluid);
        auto computed = std::make_shared<Curve>(idealGasLine(gas, spec.diagram, request.key, request.value_SI, domain, request.points));
        return impl_->curves.insert(key, std::shared_ptr<const Curve>(computed));
    }

    // The plot of the chart for this slot, made once for its lines (thermo.native._plot), always from the absolute box
    const std::string name = chart.fluid.portFluidString(request.key);
    const std::string plotKey = name + "|" + chart.fingerprint + "|" + spec.diagram.y + "|" + spec.diagram.x + "|" + keyBox(domain.tp_box) +
                                "|" + std::to_string(slot);
    std::shared_ptr<CoolProp::Plot::PropertyPlot> plot = impl_->plots.getOrCompute(plotKey, [&] {
        try {
            const std::array<double, 4>& box = domain.tp_box;
            return std::make_shared<CoolProp::Plot::PropertyPlot>(
                name, static_cast<CoolProp::parameters>(coolpropParameter(spec.diagram.y)),
                static_cast<CoolProp::parameters>(coolpropParameter(spec.diagram.x)), box[0], box[1], box[2], box[3]);
        } catch (const std::exception& e) {
            throw DiagramError("CoolProp cannot draw a " + spec.diagram.name() + " chart of " + chart.fluid.name() + ": " + e.what());
        }
    });
    CoolProp::Plot::CalcSettings settings;
    settings.points = request.points;
    settings.saturation_points = request.points;
    settings.solver = calc.solver;
    settings.critical_margin = calc.critical_margin;
    settings.saturation_spacing = calc.saturation_spacing;
    auto computed = std::make_shared<Curve>();
    try {
        CoolProp::Plot::Isoline line =
            plot->calc_isolines(static_cast<CoolProp::parameters>(coolpropParameter(request.key)), {request.value_SI}, settings).at(0);
        computed->x = std::move(line.x);
        computed->y = std::move(line.y);
        computed->n_failed = line.n_failed;
        computed->n_fallback = line.n_fallback;
        computed->n_inserted = line.n_inserted;
        computed->n_evaluated = line.n_evaluated;
        computed->n_outside = line.n_outside;
        computed->method = line.method;
        computed->reason = line.reason;
    } catch (const std::exception& e) {
        // An argument the port refuses for this line (a quality line of a fluid without dome): a diagnostic, not an error
        computed->reason = e.what();
    }
    return impl_->curves.insert(key, std::shared_ptr<const Curve>(computed));
}

std::vector<PlannedCurve> Engine::curves(const Chart& chart, const Domain& domain, const std::vector<CurveRequest>& requests,
                                         json& diagnostics, std::vector<std::string>& warnings) {
    // The curves not kept yet, computed by up to MAX_CURVE_THREADS threads (this one included) that take the next
    // request in turn, each with its own plot. A plot that cannot be made stops them all and is reported.
    std::vector<std::shared_ptr<const Curve>> computed(requests.size());
    std::vector<size_t> missing;
    for (size_t i = 0; i < requests.size(); ++i) {
        if (std::optional<std::shared_ptr<const Curve>> found = impl_->curves.get(curveKey(chart, domain, requests[i])))
            computed[i] = *found;
        else
            missing.push_back(i);
    }
    const size_t cores = std::max(1u, std::thread::hardware_concurrency());
    const size_t threads = std::min({static_cast<size_t>(MAX_CURVE_THREADS), cores, missing.size()});
    std::atomic<size_t> next{0};
    std::vector<std::exception_ptr> failures(threads);
    auto work = [&](size_t slot) {
        try {
            for (size_t k; (k = next++) < missing.size();) computed[missing[k]] = curve(chart, domain, requests[missing[k]], slot);
        } catch (...) {
            failures[slot] = std::current_exception();
            next = missing.size();
        }
    };
    std::vector<std::thread> helpers;
    for (size_t slot = 1; slot < threads; ++slot) helpers.emplace_back(work, slot);
    if (threads > 0) work(0);
    for (std::thread& helper : helpers) helper.join();
    for (const std::exception_ptr& failure : failures)
        if (failure) std::rethrow_exception(failure);

    std::vector<PlannedCurve> result;
    const std::string fluidName = chart.fluid.name();
    const double tolerance = chart.spec.calc.warn_failed_fraction;
    bool pseudoPureWarned = false;
    if (!diagnostics.is_array()) diagnostics = json::array();
    for (size_t i = 0; i < requests.size(); ++i) {
        const CurveRequest& request = requests[i];
        const std::shared_ptr<const Curve>& c = computed[i];
        const size_t total = c->x.size();
        const std::string line = request.key + "=" + coolplot::general_text(request.value_SI, 6) + " in " + fluidName;
        bool drawn = true;
        if (c->nValid() < 2) {
            drawn = false;
            if (c->n_failed == 0 && c->n_outside > 0) {
                // the line lies outside the calculation domain: nothing to draw, nothing wrong
            } else if (request.key == "Q" && !withFluid(chart.fluid, [](PropertyProvider& f) { return f.pseudoPureBlend(); }).empty()) {
                if (!pseudoPureWarned) {  // one warning for the fluid explains every missing quality line
                    pseudoPureWarned = true;
                    warnings.push_back(fluidName + " is CoolProp's pseudo-pure model of a blend: it has no two-phase states, so "
                                       "the lines of constant vapour quality are left out and the dome has no glide.");
                }
            } else {
                const std::string why = c->reason.empty() ? "CoolProp failed for " + std::to_string(c->n_failed) + " of " +
                                                                std::to_string(total) + " points"
                                                          : c->reason;
                warnings.push_back("Skipping the line " + line + ": " + why + ".");
            }
        } else {
            if (!c->reason.empty()) warnings.push_back("The line " + line + " is incomplete: " + c->reason + ".");
            if (total && c->nGaps() > tolerance * total)
                warnings.push_back("The line " + line + " has gaps: CoolProp failed for " + std::to_string(c->nGaps()) + " of " +
                                   std::to_string(total) + " points inside the line.");
            result.push_back({request, c});
        }
        diagnostics.push_back({{"item_id", request.item_id},
                               {"n_failed", c->n_failed},
                               {"n_fallback", c->n_fallback},
                               {"n_inserted", c->n_inserted},
                               {"n_evaluated", c->n_evaluated},
                               {"n_outside", c->n_outside},
                               {"method", c->method},
                               {"reason", c->reason},
                               {"drawn", drawn}});
    }
    return result;
}

// ---------------------------------------------------------------------------
// The compose record (corelib.compose_record)
// ---------------------------------------------------------------------------

json composeRecord(const json& specRecord, const PlanSpec& spec, const Domain& domain, const std::vector<PlannedCurve>& curves,
                   const std::vector<ResolvedCycle>& cycles, const std::optional<std::pair<double, double>>& plotSizePt,
                   const json& previousLabels, const std::vector<std::string>& sources) {
    json ranges = json::object();
    for (const auto& [key, range] : domain.ranges_SI) ranges[key] = json::array({range.first, range.second});
    json domainRecord = {{"tp_box", domain.tp_box}, {"ranges_SI", ranges}};
    // Written only when there are some, as CoolPlot does since 1.5
    if (!spec.marked_points.empty() && !domain.points.empty()) {
        json points = json::object();
        for (const auto& [name, state] : domain.points) points[name] = stateRecord(state);
        domainRecord["points"] = points;
    }
    if (spec.saturation_scale && !domain.scale.empty()) {
        json scale = json::array();
        for (const Range& tick : domain.scale) scale.push_back(json::array({tick.first, tick.second}));
        domainRecord["scale"] = scale;
    }
    if (spec.dome_labels && !domain.dome.empty()) {
        json dome = json::array();
        for (const DomeEnd& end : domain.dome) dome.push_back(json::array({end.value_SI, stateRecord(end.bubble), stateRecord(end.dew)}));
        domainRecord["dome"] = dome;
    }
    json curvesRecord = json::array();
    for (const auto& [request, curve] : curves) curvesRecord.push_back({{"request", request.toJson()}, {"curve", curve->toJson()}});
    json cyclesRecord = json::array();
    for (const ResolvedCycle& cycle : cycles) {
        json states = json::object();
        for (const auto& [number, state] : cycle.states) states[std::to_string(number)] = stateRecord(state);
        json chains = json::array();
        for (const std::vector<State>& chain : cycle.chains) {
            json c = json::array();
            for (const State& state : chain) c.push_back(stateRecord(state));
            chains.push_back(c);
        }
        cyclesRecord.push_back({{"states", states}, {"chains", chains}});
    }
    json record = {{"spec", specRecord},
                   {"domain", domainRecord},
                   {"curves", curvesRecord},
                   {"cycles", cyclesRecord},
                   {"plot_size_pt", plotSizePt ? json::array({plotSizePt->first, plotSizePt->second}) : json()},
                   {"previous_labels", previousLabels.is_array() ? previousLabels : json()}};
    // Only written when there are some, so the record of a chart without badge is the same as before 1.5
    if (!sources.empty()) record["sources"] = sources;
    return record;
}

// ---------------------------------------------------------------------------
// The scene of a chart request, and its files
// ---------------------------------------------------------------------------

namespace {

coolplot::Scene composeScene(const json& record) {
    try {
        return coolplot::compose(record);
    } catch (const coolplot::Error& e) {
        throw DiagramError(std::string("coolplot-core refused the chart: ") + e.what());
    } catch (const json::exception& e) {
        throw DiagramError(std::string("coolplot-core refused the chart: ") + e.what());
    }
}

json isolineLabels(const coolplot::Scene& scene) {
    json labels = json::array();
    for (const coolplot::Item& item : scene.items)
        if (item.role == "isoline_label") labels.push_back(item.to_json());
    return labels;
}

coolplot::SvgOptions svgOptions(const std::pair<double, double>& figure) {
    coolplot::SvgOptions options;
    options.width = static_cast<int>(std::lround(figure.first * coolplot::PX_PER_PT));
    options.height = static_cast<int>(std::lround(figure.second * coolplot::PX_PER_PT));
    return options;
}

coolplot::PdfOptions pdfOptions(const std::pair<double, double>& figure) {
    coolplot::PdfOptions options;
    options.width_pt = figure.first;
    options.height_pt = figure.second;
    return options;
}

/// A chart request composed into a scene, with what scene() reports of it.
struct ComposedChart {
    coolplot::DiagramSpec core;        ///< coolplot-core's reading of the spec (units, symbols)
    std::string xKey, yKey;            ///< the axes
    json specRecord;                   ///< the spec as composed, fluid and cycles written
    std::vector<ResolvedCycle> cycles;
    coolplot::Scene scene;
    std::pair<double, double> figure;  ///< the page in points
    json curveDiagnostics = json::array();
    std::vector<std::string> warnings;
    double curvesMs = 0.0, composeMs = 0.0;
};

/**
 * The scene of a chart request for the plot area of a writer (the SVG writer, or the PDF writer when `pdf`): composed
 * for the plot area of an empty scene, then again if drawing it changes the area by more than 2 % (Chart.draw); only
 * the labels change, nothing is recalculated.
 */
ComposedChart composeChart(const SceneRequest& request, bool pdf) {
    const auto start = std::chrono::steady_clock::now();
    ComposedChart out;
    ChartFluid fluid = chartFluid(request.fluid, request.backend, request.fractions, request.basis, request.p_Pa);
    if (!request.spec.is_object()) throw DiagramError("The chart spec has to be a JSON object (CoolPlot's DiagramSpec record).");
    out.specRecord = request.spec;
    out.specRecord.erase("fluid");  // the server writes them
    out.specRecord.erase("cycles");
    // coolplot-core reads the spec first, so that a value it refuses (a theme, a unit system, a label mode) is
    // reported before any property is computed
    try {
        out.core = coolplot::DiagramSpec::from_json(out.specRecord);
    } catch (const coolplot::Error& e) {
        throw DiagramError(e.what());
    } catch (const json::exception& e) {
        throw DiagramError(std::string("The chart spec is not a valid CoolPlot spec record: ") + e.what());
    }
    const Chart chart(fluid, PlanSpec::fromJson(out.specRecord));
    checkChartKind(fluid, chart.spec.diagram);
    if (fluid.isIdealGas()) checkIdealGasDiagram(chart.spec.diagram);
    out.xKey = chart.spec.diagram.x;
    out.yKey = chart.spec.diagram.y;
    out.figure = {out.core.settings.page.width_pt, out.core.settings.page.height_pt};
    if (request.figureSizePt) out.figure = *request.figureSizePt;
    if (!(out.figure.first >= 50.0 && out.figure.first <= 10000.0 && out.figure.second >= 50.0 && out.figure.second <= 10000.0))
        throw DiagramError("The figure size is two numbers of points between 50 and 10000.");
    out.specRecord["fluid"] = fluid.record();
    const std::vector<std::string> sources = sourceLines(out.specRecord);  // Engine.source_lines

    Engine& engine = Engine::shared();
    const Domain domain = engine.domain(chart, out.warnings);
    if (chart.spec.humidAir())  // the region drawn, PsychrometricChart's default included
        out.specRecord["limits"] = chart.spec.region.toJson();
    else if (!domain.saturation)  // CoolPlot draws a fluid without saturation curve from absolute limits only
        out.specRecord["limits"] = {{"T_min", domain.tp_box[0]}, {"T_max", domain.tp_box[1]}, {"p_min", domain.tp_box[2]},
                                    {"p_max", domain.tp_box[3]}};
    const std::vector<PlannedCurve> curves = engine.curves(chart, domain, plan(chart.spec, domain), out.curveDiagnostics, out.warnings);
    out.cycles = engine.withFluid(fluid, [&](PropertyProvider& model) {
        return resolveCycles(request.cycles, request.solution, model, chart.spec.diagram, out.warnings);
    });
    out.specRecord["cycles"] = json::array();
    for (const ResolvedCycle& cycle : out.cycles) out.specRecord["cycles"].push_back(cycle.spec);
    out.curvesMs = elapsedMs(start);

    const auto composeStart = std::chrono::steady_clock::now();
    auto plotSize = [&](const coolplot::Scene& scene) {
        return pdf ? coolplot::pdf_plot_size_pt(scene, pdfOptions(out.figure)) : coolplot::svg_plot_size_pt(scene, svgOptions(out.figure));
    };
    std::pair<double, double> size = plotSize(coolplot::Scene());
    json record = composeRecord(out.specRecord, chart.spec, domain, curves, out.cycles, size, request.previousLabels, sources);
    out.scene = composeScene(record);
    const std::pair<double, double> drawn = plotSize(out.scene);
    if (std::fabs(drawn.first - size.first) > 0.02 * size.first || std::fabs(drawn.second - size.second) > 0.02 * size.second) {
        size = drawn;
        record["plot_size_pt"] = json::array({size.first, size.second});
        record["previous_labels"] = isolineLabels(out.scene);
        out.scene = composeScene(record);
    }
    out.composeMs = elapsedMs(composeStart);
    return out;
}

}  // namespace

json scene(const SceneRequest& request) {
    ComposedChart chart = composeChart(request, false);
    const coolplot::SvgOptions options = svgOptions(chart.figure);
    const auto svgStart = std::chrono::steady_clock::now();
    const std::string svg = coolplot::scene_to_svg(chart.scene, options);
    const double svgMs = elapsedMs(svgStart);

    // The axes, with the units the GUI needs to turn a zoom window back into SI values
    const coolplot::AxesSpec& axes = chart.scene.axes;
    auto unitRecord = [&](const std::string& key) {
        const coolplot::Unit& unit = chart.core.units.at(key);
        return json{{"label", coolplot::unit_text(unit.label, chart.core.settings.labels.symbols)}, {"scale", unit.scale}, {"offset", unit.offset}};
    };
    const std::pair<double, double> plotSize = coolplot::svg_plot_size_pt(chart.scene, options);
    json axesRecord = {{"x_key", chart.xKey},
                       {"y_key", chart.yKey},
                       {"x_label", axes.x_label},
                       {"y_label", axes.y_label},
                       {"x_limits", json::array({axes.x_limits.lo, axes.x_limits.hi})},
                       {"y_limits", json::array({axes.y_limits.lo, axes.y_limits.hi})},
                       {"x_log", axes.x_log},
                       {"y_log", axes.y_log},
                       {"x_unit", unitRecord(chart.xKey)},
                       {"y_unit", unitRecord(chart.yKey)},
                       {"plot_size_pt", json::array({plotSize.first, plotSize.second})}};
    json states = json::array();
    for (const ResolvedCycle& cycle : chart.cycles)
        for (const auto& [number, state] : cycle.states) {
            const auto label = cycle.labels.find(number);
            states.push_back({{"cycle", cycle.index},
                              {"number", number},
                              {"label", label == cycle.labels.end() ? "" : label->second},
                              {"vars", cycle.vars.at(number)},
                              {"values_SI", stateRecord(state)}});
        }
    return {{"svg", svg},
            {"axes", axesRecord},
            {"labels", isolineLabels(chart.scene)},
            {"states", states},
            {"spec", chart.specRecord},
            {"diagnostics", {{"curves", chart.curveDiagnostics}, {"warnings", chart.warnings}}},
            {"timing_ms", {{"curves", chart.curvesMs}, {"compose", chart.composeMs}, {"svg", svgMs}}}};
}

std::string chartFile(const SceneRequest& request, const std::string& format, std::vector<std::string>* warnings) {
    if (format != "svg" && format != "pdf") throw DiagramError("A chart is written as \"svg\" or \"pdf\", not \"" + format + "\".");
    const ComposedChart chart = composeChart(request, format == "pdf");
    if (warnings) warnings->insert(warnings->end(), chart.warnings.begin(), chart.warnings.end());
    if (format == "svg") return coolplot::scene_to_svg(chart.scene, svgOptions(chart.figure));
    try {
        return coolplot::scene_to_pdf(chart.scene, pdfOptions(chart.figure));
    } catch (const coolplot::Error& e) {
        // No TrueType font of the theme's families on this platform: the server's set-up, not the request
        throw std::runtime_error(std::string("The PDF cannot be written: ") + e.what());
    }
}

std::string chartFormat(const std::string& path) {
    std::string extension = std::filesystem::path(path).extension().string();
    for (char& c : extension) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (extension == ".svg") return "svg";
    if (extension == ".pdf") return "pdf";
    throw DiagramError("A chart file is .svg or .pdf, not '" + path + "'.");
}

const json& findDefinition(const json& document, const std::string& name) {
    checkDefinitions(document);
    static const json none = json::array();
    const json& diagrams = document.contains("diagrams") ? document.at("diagrams") : none;
    if (diagrams.empty()) throw DiagramError("The diagram definitions hold no diagram.");
    if (name.empty()) return diagrams.at(0);
    std::string names;
    for (const json& diagram : diagrams) {
        const std::string own = diagram.value("name", std::string());
        if (own == name) return diagram;
        names += (names.empty() ? "'" : ", '") + own + "'";
    }
    throw DiagramError("The diagram definitions have no diagram named '" + name + "' (they have " + names + ").");
}

void readConcentration(const json& record, std::vector<double>& fractions, std::string& basis) {
    const json given = record.value("fractions", json());
    if (!given.is_null()) {
        if (!given.is_array()) throw DiagramError("\"fractions\" is a list of numbers, e.g. [0.3].");
        for (const json& f : given) {
            if (!f.is_number()) throw DiagramError("\"fractions\" is a list of numbers, e.g. [0.3].");
            fractions.push_back(f.get<double>());
        }
    }
    const json b = record.value("basis", json());
    if (!b.is_null()) {
        if (!b.is_string()) throw DiagramError("\"basis\" is \"mass\" for the concentration of a solution.");
        basis = b.get<std::string>();
    }
}

std::optional<double> readTotalPressure(const json& record) {
    const json given = record.is_object() ? record.value("p_Pa", json()) : json();
    if (given.is_null()) return std::nullopt;
    if (!given.is_number()) throw DiagramError("\"p_Pa\" is the total pressure of a chart of humid air, a number of Pa.");
    return given.get<double>();
}

SceneRequest definitionRequest(const json& definition) {
    if (!definition.is_object() || !definition.contains("fluid") || !definition.at("fluid").is_string())
        throw DiagramError("A diagram definition needs the name of its fluid (\"fluid\").");
    SceneRequest request;
    request.fluid = definition.at("fluid").get<std::string>();
    readConcentration(definition, request.fractions, request.basis);
    request.p_Pa = readTotalPressure(definition);
    request.spec = definition.value("spec", json::object());
    request.cycles = definition.value("cycles", json::array());
    return request;
}

std::vector<std::string> writeDiagramFile(const json& document, const std::string& name, const SolutionValues& solution,
                                          const std::string& backend, const std::string& path) {
    const std::string format = chartFormat(path);
    SceneRequest request = definitionRequest(findDefinition(document, name));
    request.backend = backend.empty() ? "HEOS" : backend;
    request.solution = &solution;
    std::vector<std::string> warnings;
    const std::string file = chartFile(request, format, &warnings);
    std::ofstream out(path, std::ios::binary);
    if (!out.is_open() || !out.write(file.data(), static_cast<std::streamsize>(file.size())))
        throw std::runtime_error("Cannot write the chart to '" + path + "'.");
    return warnings;
}

std::string reportDiagramFile(size_t index) {
    return index == 0 ? "thermo_diagrams.pdf" : "thermo_diagram_" + std::to_string(index + 1) + ".pdf";
}

std::vector<std::pair<std::string, std::string>> reportDiagrams(const json& document, const SolutionValues* solution,
                                                                const std::string& backend, std::vector<std::string>& warnings) {
    std::vector<std::pair<std::string, std::string>> files;
    if (!document.is_object() || !document.contains("diagrams") || !document.at("diagrams").is_array()) return files;
    const json& diagrams = document.at("diagrams");
    for (size_t i = 0; i < diagrams.size() && i < REPORT_DIAGRAMS; ++i) {
        const std::string label = diagrams[i].is_object() ? diagrams[i].value("name", std::string()) : std::string();
        try {
            SceneRequest request = definitionRequest(diagrams[i]);
            request.backend = backend.empty() ? "HEOS" : backend;
            request.solution = solution;
            if (!request.spec.is_object()) request.spec = json::object();
            if (!request.spec["settings"].is_object()) request.spec["settings"] = json::object();
            request.spec["settings"]["look"] = "print";  // the report is printed
            files.emplace_back(reportDiagramFile(i), chartFile(request, "pdf", &warnings));
        } catch (const std::exception& e) {
            warnings.push_back("The diagram '" + label + "' is left out of the report: " + e.what());
        }
    }
    return files;
}

// ---------------------------------------------------------------------------
// The state at a point of a chart
// ---------------------------------------------------------------------------

json stateAt(const ChartFluid& chart, const std::string& diagram, double x_SI, double y_SI) {
    coolplot::DiagramType type;
    try {
        type = coolplot::DiagramType::parse(diagram);
    } catch (const coolplot::Error& e) {
        throw DiagramError(e.what());
    }
    checkChartKind(chart, type);
    if (!std::isfinite(x_SI) || !std::isfinite(y_SI)) throw DiagramError("The point of a chart is two finite numbers in SI units.");
    if (chart.isIdealGas()) checkIdealGasDiagram(type);
    return Engine::shared().withFluid(chart, [&](PropertyProvider& model) -> json {
        try {
            model.update(type.x, x_SI, type.y, y_SI);
        } catch (const std::exception&) {
            return nullptr;  // no state there
        }
        State state;
        for (const std::string& key : model.quantities()) {
            try {
                state[key] = model.get(key);
            } catch (const std::exception&) {
                state[key] = std::numeric_limits<double>::quiet_NaN();
            }
        }
        // Outside the model CoolProp may still answer, by extrapolation: no state there either
        const double T = state.at("T"), p = state.at("p");
        const double T_min = model.tMinK(), T_max = model.tMaxK(), p_min = model.pMinPa(), p_max = model.pMaxPa();
        if (!std::isfinite(T) || !std::isfinite(p) || (std::isfinite(T_min) && T < T_min) || (std::isfinite(T_max) && T > T_max) ||
            (std::isfinite(p_min) && p < p_min) || (std::isfinite(p_max) && p > p_max))
            return nullptr;
        // Humid air is no state above saturation (fog) nor below dry air
        if (chart.isHumidAir() && !(state.at("RH") >= 0.0 && state.at("RH") <= 1.0 + 1e-9 && state.at("W") >= 0.0)) return nullptr;
        json out = json::object();
        for (const std::string& key : model.quantities()) out[key] = finiteOrNull(state.at(key));
        return out;
    });
}

// ---------------------------------------------------------------------------
// Options, units, definitions
// ---------------------------------------------------------------------------

namespace {

/// The diagram types of a fluid the Diagrams tab offers, in its order.
const std::vector<std::string>& fluidDiagramTypes() {
    static const std::vector<std::string> names = {"ph", "Ts", "hs", "ps", "prho", "Trho", "pT", "Th"};
    return names;
}

/// The charts of humid air: the psychrometric chart (W over T) and the Mollier h,x chart (Thx over W).
const std::vector<std::string>& humidAirDiagramTypes() {
    static const std::vector<std::string> names = {"WT", "ThxW"};
    return names;
}

/**
 * What each kind of fluid offers (options().fluid_kinds). Incompressible liquids: the charts with the pressure on the
 * vertical axis, whose axes CoolProp's INCOMP backend has a flash for (in the others, T-s, h-s, T-rho, T-h, every
 * state lies near one line: h, s and rho depend on T alone but for a small pressure term); the lines CoolProp::Plot
 * computes for them (it has no (rho, T) nor (u, p) flash of INCOMP). Ideal gases: every chart but T-h, every family
 * but the quality. Humid air: the psychrometric and Mollier charts at a total pressure (`total_pressure`, the
 * default, within `total_pressure_range`), with no dome, critical point nor saturation scale (`saturation` false: its
 * saturation line, RH = 1, comes with the chart). `default_region`: the region of a new chart; ACHP for mixtures,
 * whose DEF region reaches down to the triple point of the lightest component (0.4 Pa for propane/isobutane), where
 * their lines are slow to compute (CSL-0324: 6.0 s with DEF, 2.8 s with ACHP); null for the fluids without saturation
 * curve; for humid air the Region record of a psychrometric chart (T and W; the Mollier chart bounds Thx instead of
 * T, within `temperature_range`). `chart_kind`: CoolPlot's kind of the charts, "fluid" or "humid_air".
 */
json fluidKinds() {
    const json all = fluidDiagramTypes();
    json idealTypes = json::array();
    for (const std::string& name : fluidDiagramTypes()) {
        try {
            checkIdealGasDiagram(coolplot::DiagramType::parse(name));
            idealTypes.push_back(name);
        } catch (const DiagramError&) {  // a chart an ideal gas does not have
        }
    }
    json saturated = {{"diagram_types", all}, {"families", {"T", "p", "h", "s", "rho", "u", "Q"}}, {"saturation", true},
                      {"fractions", false}, {"default_box", nullptr}, {"default_region", "DEF"}, {"chart_kind", "fluid"}};
    json mixture = saturated;
    mixture["default_region"] = "ACHP";
    return {{"Real", saturated},
            {"Mixture", mixture},
            {"Incompressible", {{"diagram_types", {"ph", "ps", "prho", "pT"}}, {"families", {"T", "p", "h", "s"}}, {"saturation", false},
                                {"fractions", true}, {"default_box", {nullptr, nullptr, LIQUID_P_MIN_PA, LIQUID_P_MAX_PA}},
                                {"default_region", nullptr}, {"chart_kind", "fluid"}}},
            {"IdealGas", {{"diagram_types", idealTypes}, {"families", {"T", "p", "h", "s", "rho", "u"}}, {"saturation", false},
                          {"fractions", false}, {"default_box", IDEAL_GAS_BOX}, {"default_region", nullptr}, {"chart_kind", "fluid"}}},
            {"HumidAir", {{"diagram_types", humidAirDiagramTypes()}, {"families", {"T", "RH", "Twb", "h", "v"}}, {"saturation", false},
                          {"fractions", false}, {"default_box", nullptr}, {"default_region", defaultHumidAirRegion("T").toJson()},
                          {"chart_kind", "humid_air"}, {"total_pressure", DEFAULT_HUMID_AIR_P_PA},
                          {"total_pressure_range", {HUMID_AIR_P_MIN_PA, HUMID_AIR_P_MAX_PA}},
                          {"temperature_range", {HUMID_AIR_T_MIN_K, HUMID_AIR_T_MAX_K}}}}};
}

}  // namespace

json options() {
    json types = json::array();
    std::vector<std::string> offered = fluidDiagramTypes();
    offered.insert(offered.end(), humidAirDiagramTypes().begin(), humidAirDiagramTypes().end());
    for (const std::string& name : offered) {
        const coolplot::DiagramType type = coolplot::DiagramType::parse(name);
        const coolplot::Quantity& y = coolplot::get_quantity(type.y);
        const coolplot::Quantity& x = coolplot::get_quantity(type.x);
        std::string xLabel = x.label;
        xLabel[0] = static_cast<char>(std::tolower(static_cast<unsigned char>(xLabel[0])));
        std::string title = y.label + "-" + xLabel;
        if (type.kind == "humid_air") title = type.y == "Thx" ? "Mollier h,x chart" : "Psychrometric chart";
        json families = json::array();
        for (const coolplot::Family& family : type.families())
            families.push_back({{"key", family.key},
                                {"label", type.quantity_label(family.key)},
                                {"symbol", type.symbol(family.key)},
                                {"default_num", family.default_num}});
        types.push_back({{"name", type.name()},
                         {"title", title},
                         {"kind", type.kind},
                         {"x", type.x},
                         {"y", type.y},
                         {"x_log", type.x_log()},
                         {"y_log", type.y_log()},
                         {"families", families}});
    }
    json units = json::array();
    for (const char* name : {"SI", "KSI", "engineering", "IP"})
        units.push_back({{"name", name}, {"record", name}, {"units", unitTable(coolplot::get_unit_system(name))}});
    const json custom = coolsolveUnitSystem();
    units.push_back({{"name", "coolsolve"}, {"record", custom}, {"units", unitTable(coolplot::UnitSystem::from_json(custom))}});
    json regions = json::array();
    for (const std::string& name : tpLimitsPresetNames()) regions.push_back({{"name", name}, {"limits", tpLimitsPreset(name).toJson()}});
    return {{"diagram_types", types},
            {"fluid_types", {"Real", "Mixture", "Incompressible", "IdealGas", "HumidAir"}},
            {"fluid_kinds", fluidKinds()},
            {"unit_systems", units},
            {"default_unit_system", "engineering"},
            {"regions", regions},
            {"themes", {"default", "dark", "print", "coolpack"}},
            {"label_modes", {"rails", "edge", "margin", "inline"}},
            {"marked_points", {"critical", "triple"}},
            {"coolplot", {{"core_version", coolplot_version()}, {"commit", COOLSOLVE_COOLPLOT_COMMIT}}}};
}

json emptyDefinitions() {
    return {{"format", "coolsolve-diagrams"}, {"version", 1}, {"diagrams", json::array()}};
}

void checkDefinitions(const json& document) {
    if (!document.is_object()) throw DiagramError("The diagram definitions are a JSON object {\"format\", \"version\", \"diagrams\"}.");
    if (document.contains("format") &&
        (!document.at("format").is_string() || document.at("format").get<std::string>() != "coolsolve-diagrams"))
        throw DiagramError("The diagram definitions have the format \"coolsolve-diagrams\", not " + document.at("format").dump() + ".");
    if (document.contains("version") && !document.at("version").is_number())
        throw DiagramError("The version of the diagram definitions is a number, not " + document.at("version").dump() + ".");
    if (document.contains("diagrams")) {
        const json& diagrams = document.at("diagrams");
        if (!diagrams.is_array()) throw DiagramError("\"diagrams\" is a list of diagram definitions.");
        for (const json& d : diagrams)
            if (!d.is_object()) throw DiagramError("A diagram definition is an object, not " + d.dump() + ".");
    }
}

}  // namespace diagram
}  // namespace coolsolve
