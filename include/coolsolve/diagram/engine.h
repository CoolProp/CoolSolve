// Computing what a plan asks for, keeping the results, and turning a chart request into a scene: the plot and
// curve caches, the compose record, the plot size in two passes, the diagnostics, the state at a point of a chart,
// the options of the Diagrams tab and the diagram definitions document. The cycles are in cycles.h.
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): engine/engine.py (Engine: domain,
// curves, warnings), engine/cache.py (CurveCache), thermo/native.py (the plots of the port), corelib.py
// (compose_record) and chart.py (the plot area measured again after drawing, Chart.draw).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
// See docs/coolplot_integration.md §3 for the REST contract this serves.
#pragma once

#include "coolsolve/diagram/cycles.h"
#include "coolsolve/diagram/planning.h"

#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace coolsolve {
namespace diagram {

/**
 * @brief One computed line (the curve record of CoolPlot's contracts): x and y in SI units, NaN where a point
 *        failed, with the diagnostics of CoolProp::Plot::Isoline.
 */
struct Curve {
    std::vector<double> x, y;
    int n_failed = 0;
    int n_fallback = 0;
    int n_inserted = 0;
    int n_evaluated = 0;
    int n_outside = 0;
    std::string method = "flash";
    std::string reason;

    /// Points with finite x and y.
    size_t nValid() const;
    /// Invalid points between the first and the last valid one.
    size_t nGaps() const;
    /// The curve record (non-finite values become null).
    json toJson() const;
};

/// A curve request with its curve.
using PlannedCurve = std::pair<CurveRequest, std::shared_ptr<const Curve>>;

/// What one chart is computed from: the fluid, the spec and the fingerprint of the CoolProp that computes it.
struct Chart {
    ChartFluid fluid;
    PlanSpec spec;
    std::string fingerprint;  ///< CoolProp's git revision and a digest of its configuration (Fluid.fingerprint)

    Chart(ChartFluid fluid, PlanSpec spec);
};

/// The most threads that compute the lines of one chart (the calling one included). The online server is shared by
/// its sessions: a chart request takes at most this many cores, whatever the machine has (and never more than it
/// has). The lines of a blend take seconds on one core (docs/coolplot_integration.md §6.3).
constexpr unsigned MAX_CURVE_THREADS = 4;

/**
 * @brief The process-wide store of everything expensive (coolplot.engine.Engine and CurveCache): fluid models
 *        (property providers), calculation domains, saturation states, CoolProp::Plot::PropertyPlot objects (they
 *        hold the critical point and the phase envelope) and curves, each keyed by what determines it and bounded
 *        (least recently used dropped first).
 *
 * Thread safe: the server's sessions share it. A fluid model is used by one request at a time (it holds a lock);
 * a plot serialises its own calls, so the lines of one chart are computed on up to MAX_CURVE_THREADS threads with
 * one plot each (thermo.native._plot keeps one per thread for the same reason). Two requests may compute the same
 * curve at once; the first one stored wins.
 */
class Engine {
public:
    /// The engine of the process.
    static Engine& shared();

    Engine();
    ~Engine();
    Engine(const Engine&) = delete;
    Engine& operator=(const Engine&) = delete;

    /**
     * @brief The calculation domain of a chart (Engine.domain): the T-p box, the ranges, and the saturation states
     *        the spec asks for (anchored families, marked points, saturation scale, dome ends).
     * @param warnings receives the warnings of saturation states CoolProp could not compute
     * @throws DiagramError for a region CoolPlot refuses or CoolProp cannot compute
     */
    Domain domain(const Chart& chart, std::vector<std::string>& warnings);

    /**
     * @brief One curve, from the cache or computed: by CoolProp::Plot, or by CoolSolve's ideal-gas provider for an
     *        ideal gas (idealGasLine). A line that cannot be computed comes back without points and with its reason.
     * @param slot which of the chart's plots computes it: one per thread, so that threads do not wait for each other
     * @throws DiagramError when the chart's plot cannot be made at all
     */
    std::shared_ptr<const Curve> curve(const Chart& chart, const Domain& domain, const CurveRequest& request, size_t slot = 0);

    /**
     * @brief The curves of all requests that have at least two valid points (Engine.curves), in the order of the
     *        requests. The curves not kept yet are computed on up to MAX_CURVE_THREADS threads. Every request gets an
     *        entry in `diagnostics` (item_id, the counters, method, reason, drawn); a line left out or with too many
     *        gaps gets a warning.
     */
    std::vector<PlannedCurve> curves(const Chart& chart, const Domain& domain, const std::vector<CurveRequest>& requests,
                                     json& diagnostics, std::vector<std::string>& warnings);

    /// Run `work(model)` with the chart fluid's model (its property provider), locked for this call.
    template <class Work>
    auto withFluid(const ChartFluid& fluid, Work&& work) -> decltype(work(std::declval<PropertyProvider&>()));

    /// Forget everything (the next chart is computed cold).
    void clear();
    /// The number of curves kept.
    size_t curveCount() const;

private:
    struct Impl;
    struct FluidEntry;
    std::shared_ptr<FluidEntry> fluidEntry(const ChartFluid& fluid);
    static PropertyProvider& model(FluidEntry& entry);
    static std::mutex& mutex(FluidEntry& entry);
    std::unique_ptr<Impl> impl_;
};

template <class Work>
auto Engine::withFluid(const ChartFluid& fluid, Work&& work) -> decltype(work(std::declval<PropertyProvider&>())) {
    std::shared_ptr<FluidEntry> entry = fluidEntry(fluid);
    std::lock_guard<std::mutex> lock(mutex(*entry));
    return work(model(*entry));
}

/**
 * @brief The compose record of coolplot-core (coolplot/compose.h; corelib.compose_record) for a chart.
 * @param specRecord the spec record with "fluid" and "cycles" written
 */
json composeRecord(const json& specRecord, const PlanSpec& spec, const Domain& domain, const std::vector<PlannedCurve>& curves,
                   const std::vector<ResolvedCycle>& cycles, const std::optional<std::pair<double, double>>& plotSizePt,
                   const json& previousLabels);

/// A chart request of the Diagrams tab (POST /api/v1/diagram/scene).
struct SceneRequest {
    std::string fluid;                       ///< CoolSolve fluid name (FluidRegistry)
    std::string backend = "HEOS";            ///< the model's coolpropBackend
    std::vector<double> fractions;           ///< the mass fraction of an incompressible solution (chartFluid)
    std::string basis;                       ///< what the fractions are: "mass" (the default) for a solution
    json spec = json::object();              ///< CoolPlot spec record, without "fluid" and "cycles"
    json cycles = json::array();             ///< cycle definitions (§4.1)
    std::optional<std::pair<double, double>> figureSizePt;  ///< the whole figure; default: settings.page
    json previousLabels;                     ///< null, or the "labels" of the previous response
    const SolutionValues* solution = nullptr;  ///< the session's last solve, null without one
};

/**
 * @brief The response of POST /api/v1/diagram/scene: svg, axes, labels, states, spec, diagnostics and timing_ms
 *        (docs/coolplot_integration.md §3.1). The plot area is measured on the SVG writer's layout and the chart
 *        composed again when it differs by more than 2 % from the size laid out for (Chart.draw).
 * For a fluid without saturation curve (an incompressible liquid, an ideal gas) the composed spec carries the
 * absolute T-p box of its lines as its limits, which CoolPlot needs for such a fluid.
 * @throws DiagramError for an unknown fluid, a fluid not drawable yet, a solution without its concentration or with
 *         one outside its range, an unknown diagram type, a chart an ideal gas does not have, or a spec value CoolPlot
 *         or coolplot-core refuses
 */
json scene(const SceneRequest& request);

/**
 * @brief The full state at a point of a chart (POST /api/v1/diagram/state, the hover of the Diagrams tab).
 *
 * The fluid's state at the axis values `x_SI`, `y_SI` of `diagram`: every key of stateKeys(), null for a value that
 * cannot be computed, Q = -1 in single phase as CoolProp and CoolPlot report it. One flash with the chart fluid's
 * model (the same as the cycles and the domain use), the caches of the engine kept.
 * @param fluid the chart fluid (chartFluid: the model's backend, the concentration of a solution)
 * @return the state, or null when the point has none: the model cannot compute it, or it lies outside the
 *         temperature and pressure range of the model
 * @throws DiagramError for an unknown diagram type, a value that is not a finite number
 */
json stateAt(const ChartFluid& fluid, const std::string& diagram, double x_SI, double y_SI);

/**
 * @brief The response of GET /api/v1/diagram/options: diagram types with their families, the fluid types drawn and
 *        what each kind of fluid offers (`fluid_kinds`), unit systems (the presets and `coolsolve`, CoolSolve's units
 *        as a custom unit record), region presets, themes, label modes and the versions.
 *
 * `fluid_kinds` maps "Real", "Mixture", "Incompressible" and "IdealGas" to {"diagram_types": names offered,
 * "families": the keys of the families its lines can have, "saturation": whether it has a saturation curve (dome,
 * lines of constant quality, critical point, saturation scale), "fractions": whether it takes a concentration (a
 * mass fraction), "default_box": [T_min_K, T_max_K, p_min_Pa, p_max_Pa] of a fluid without saturation curve where
 * the spec gives no absolute limit (null: the fluid's own temperatures)}.
 */
json options();

/// An empty diagram definitions document (§3.4).
json emptyDefinitions();
/// Check a diagram definitions document loosely (an object; "format" "coolsolve-diagrams" and a numeric
/// "version" when given; "diagrams" a list when given). @throws DiagramError
void checkDefinitions(const json& document);

}  // namespace diagram
}  // namespace coolsolve
