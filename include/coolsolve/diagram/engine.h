// Computing what a plan asks for, keeping the results, and turning a chart request into a scene: the plot and
// curve caches, the compose record, the plot size in two passes, the diagnostics, the cycles of a solved model
// (straight chains in phase 1), the options of the Diagrams tab and the diagram definitions document.
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): engine/engine.py (Engine: domain,
// curves, warnings), engine/cache.py (CurveCache), thermo/native.py (the plots of the port), corelib.py
// (compose_record) and chart.py (the plot area measured again after drawing, Chart.draw).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
// See docs/coolplot_integration.md §3 for the REST contract this serves.
#pragma once

#include "coolsolve/diagram/planning.h"
#include "coolsolve/ir.h"  // CaseInsensitiveLess

#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace coolsolve {
namespace diagram {

/// Values of a solved model by variable name (case-insensitive), in CoolSolve's units (°C, Pa, J/kg, ...).
using SolutionValues = std::map<std::string, double, CaseInsensitiveLess>;

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

/**
 * @brief The process-wide store of everything expensive (coolplot.engine.Engine and CurveCache): CoolProp fluid
 *        models, calculation domains, saturation states, CoolProp::Plot::PropertyPlot objects (they hold the
 *        critical point and the phase envelope) and curves, each keyed by what determines it and bounded
 *        (least recently used dropped first).
 *
 * Thread safe: the server's sessions share it. A fluid model is used by one request at a time (it holds a lock);
 * a plot serialises its own calls. Two requests may compute the same curve at once; the first one stored wins.
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

    /// One curve, from the cache or computed by CoolProp::Plot. A line that cannot be computed comes back
    /// without points and with its reason. @throws DiagramError when the chart's plot cannot be made at all
    std::shared_ptr<const Curve> curve(const Chart& chart, const Domain& domain, const CurveRequest& request);

    /**
     * @brief The curves of all requests that have at least two valid points (Engine.curves), in the order of the
     *        requests. Every request gets an entry in `diagnostics` (item_id, the counters, method, reason, drawn);
     *        a line left out or with too many gaps gets a warning.
     */
    std::vector<PlannedCurve> curves(const Chart& chart, const Domain& domain, const std::vector<CurveRequest>& requests,
                                     json& diagnostics, std::vector<std::string>& warnings);

    /// Run `work(model)` with the chart fluid's CoolProp model, locked for this call.
    template <class Work>
    auto withFluid(const ChartFluid& fluid, Work&& work) -> decltype(work(std::declval<FluidModel&>()));

    /// Forget everything (the next chart is computed cold).
    void clear();
    /// The number of curves kept.
    size_t curveCount() const;

private:
    struct Impl;
    struct FluidEntry;
    std::shared_ptr<FluidEntry> fluidEntry(const ChartFluid& fluid);
    static FluidModel& model(FluidEntry& entry);
    static std::mutex& mutex(FluidEntry& entry);
    std::unique_ptr<Impl> impl_;
};

template <class Work>
auto Engine::withFluid(const ChartFluid& fluid, Work&& work) -> decltype(work(std::declval<FluidModel&>())) {
    std::shared_ptr<FluidEntry> entry = fluidEntry(fluid);
    std::lock_guard<std::mutex> lock(mutex(*entry));
    return work(model(*entry));
}

/**
 * @brief A cycle of a diagram definition resolved against the solved model: the CoolPlot cycle record of the spec
 *        (CycleSpec.to_dict, full states, straight paths), the states by number and the chains through them.
 */
struct ResolvedCycle {
    size_t index = 0;                           ///< position in the request's cycle definitions
    json spec;                                  ///< the entry of spec.cycles
    std::map<int, State> states;                ///< full states in SI units
    std::map<int, std::map<std::string, std::string>> vars;  ///< number -> {quantity key: variable name}
    std::map<int, std::string> labels;          ///< number -> label
    std::vector<std::vector<State>> chains;     ///< the states of each loop, in order (straight paths)
};

/**
 * @brief The cycles of diagram definitions (docs/coolplot_integration.md §4.1) from a solved model.
 *
 * A state is given by `vars` ({quantity key: variable name}) or, for `arrays` ({key: array name}), by the elements
 * `name[i]` of every array. The values (CoolSolve's units) are converted to SI and completed to full states from
 * two independent properties (the chart's axes when the model has both) with the chart fluid; the solved values
 * are kept. `loops` gives the chains of state numbers; without it, one loop through the states in ascending order
 * (closed unless `closed` is false); an empty list draws the state markers alone. Phase 1 joins the states with
 * straight lines ("paths": "auto" gets a warning). A cycle without states is skipped silently; one that cannot be
 * resolved (no solve, a variable without a value, a state CoolProp cannot compute) is left out with a warning.
 * @throws DiagramError for definitions that are not a list of objects, or two cycles of the same name
 */
std::vector<ResolvedCycle> resolveCycles(const json& definitions, const SolutionValues* solution, FluidModel& fluid,
                                         const coolplot::DiagramType& diagram, std::vector<std::string>& warnings);

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
 * @throws DiagramError for an unknown fluid, a fluid not drawable yet, an unknown diagram type or a spec value
 *         CoolPlot or coolplot-core refuses
 */
json scene(const SceneRequest& request);

/// The response of GET /api/v1/diagram/options: diagram types with their families, unit systems (the presets and
/// `coolsolve`, CoolSolve's units as a custom unit record), region presets, themes, label modes and the versions.
json options();

/// CoolSolve's fixed units as a custom unit record of coolplot-core: {"name": "coolsolve", "units": {key: [label,
/// scale, offset(, difference)]}} for every quantity of a fluid chart.
json coolsolveUnitSystem();

/// An empty diagram definitions document (§3.4).
json emptyDefinitions();
/// Check a diagram definitions document loosely (an object; "format" "coolsolve-diagrams" and a numeric
/// "version" when given; "diagrams" a list when given). @throws DiagramError
void checkDefinitions(const json& document);

}  // namespace diagram
}  // namespace coolsolve
