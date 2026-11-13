// Computing what a plan asks for, keeping the results, and turning a chart request into a scene: the plot and
// curve caches, the compose record, the plot size in two passes, the diagnostics, the state at a point of a chart,
// the SVG and PDF files of a chart (downloads, coolsolve --diagram, the figures of the LaTeX report), the options
// of the Diagrams tab and the diagram definitions document. The cycles are in cycles.h, the badge lines in references.h.
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): engine/engine.py (Engine: domain,
// curves, warnings), engine/cache.py (CurveCache), thermo/native.py (the plots of the port), corelib.py
// (compose_record), chart.py (the plot area measured again after drawing, Chart.draw) and render/pdf.py (the PDF
// writer's plot area).
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
 * @param sources the badge lines (references.h: sourceLines), written only when there are some, as CoolPlot does
 */
json composeRecord(const json& specRecord, const PlanSpec& spec, const Domain& domain, const std::vector<PlannedCurve>& curves,
                   const std::vector<ResolvedCycle>& cycles, const std::optional<std::pair<double, double>>& plotSizePt,
                   const json& previousLabels, const std::vector<std::string>& sources = {});

/// A chart request of the Diagrams tab (POST /api/v1/diagram/scene).
struct SceneRequest {
    std::string fluid;                       ///< CoolSolve fluid name (FluidRegistry)
    std::string backend = "HEOS";            ///< the model's coolpropBackend
    std::vector<double> fractions;           ///< the mass fraction of an incompressible solution (chartFluid)
    std::string basis;                       ///< what the fractions are: "mass" (the default) for a solution
    std::optional<double> p_Pa;              ///< humid air: the total pressure of the chart (default 101325 Pa)
    json spec = json::object();              ///< CoolPlot spec record, without "fluid" and "cycles"
    json cycles = json::array();             ///< cycle definitions (§4.1)
    std::optional<std::pair<double, double>> figureSizePt;  ///< the whole figure; default: settings.page
    json previousLabels;                     ///< null, or the "labels" of the previous response
    const SolutionValues* solution = nullptr;  ///< the session's last solve, null without one
};

/**
 * @brief The response of POST /api/v1/diagram/scene: svg, axes, labels, states, spec, diagnostics and timing_ms
 *        (docs/coolplot_integration.md §3.1). The plot area is measured on the SVG writer's layout and the chart
 *        composed again when it differs by more than 2 % from the size laid out for (Chart.draw). The badge holds
 *        the lines settings.branding asks for (references.h: sourceLines).
 * For a fluid without saturation curve (an incompressible liquid, an ideal gas) the composed spec carries the
 * absolute T-p box of its lines as its limits, which CoolPlot needs for such a fluid.
 * A chart of humid air (the psychrometric "WT" and Mollier "ThxW" charts) is at the request's total pressure, over the
 * spec's Region (or PsychrometricChart's default), which the composed spec carries.
 * @throws DiagramError for an unknown fluid, a fluid not drawable yet, a solution without its concentration or with
 *         one outside its range, an unknown diagram type, one that cannot show the fluid (a fluid's chart of humid air,
 *         a psychrometric chart of a fluid), a chart an ideal gas does not have, or a spec value CoolPlot or
 *         coolplot-core refuses
 */
json scene(const SceneRequest& request);

/**
 * @brief The file of a chart (POST /api/v1/diagram/export, `coolsolve --diagram`, the figures of the LaTeX report):
 *        the scene of `request` written by coolplot-core's SVG writer, or by its PDF writer with the TrueType fonts it
 *        finds on the platform embedded as subsets. The page is the request's figure size, else the spec's page
 *        (CoolPlot's settings.page, 576 x 432 pt); the plot area is measured on the writer's own layout, as scene()
 *        does.
 * @param format "svg" or "pdf"
 * @param warnings receives the warnings of the chart (lines or states left out), when given
 * @throws DiagramError as scene(), and for another format; std::runtime_error for a PDF when the platform has no
 *         TrueType font of the theme's families (a server needs a font package, e.g. fonts-dejavu-core)
 */
std::string chartFile(const SceneRequest& request, const std::string& format, std::vector<std::string>* warnings = nullptr);

/// The file format of a chart file name: "svg" or "pdf" from its extension (any case). @throws DiagramError otherwise
std::string chartFormat(const std::string& path);

/**
 * @brief The diagram named `name` of a diagram definitions document (§3.4), or its first one when `name` is empty.
 * @throws DiagramError for a document that is not one, one without diagrams, or a name it does not have
 */
const json& findDefinition(const json& document, const std::string& name);

/**
 * @brief The concentration of an incompressible solution in a chart request or a diagram definition: "fractions"
 *        (mass fractions, as CoolPlot's fluid record) and "basis" ("mass"), both optional; appended to `fractions`
 *        and written to `basis` when given. chartFluid() checks them against the fluid.
 * @throws DiagramError when "fractions" is not a list of numbers or "basis" not a string
 */
void readConcentration(const json& record, std::vector<double>& fractions, std::string& basis);

/**
 * @brief The total pressure of a chart of humid air in a chart request or a diagram definition: "p_Pa" (HumidAirSpec's
 *        field), a number of Pa; none when absent or null. chartFluid() checks it against the fluid and its range.
 * @throws DiagramError when it is not a number
 */
std::optional<double> readTotalPressure(const json& record);

/**
 * @brief The chart request of a diagram definition ({name, fluid, fractions, basis, p_Pa, spec, cycles}); the backend,
 *        the solution and the figure size are left to the caller.
 * @throws DiagramError for a definition without the name of its fluid, or with an invalid concentration
 */
SceneRequest definitionRequest(const json& definition);

/**
 * @brief `coolsolve --diagram`: the diagram `name` (the first one when empty) of a definitions document drawn with
 *        a solved model and written to `path`, as SVG or PDF after its extension, on the spec's page.
 * @param backend the model's coolpropBackend
 * @return the warnings of the chart (lines or states left out)
 * @throws DiagramError for a request the engine refuses (see chartFile, findDefinition); std::runtime_error when the
 *         file cannot be written or no font is found for a PDF
 */
std::vector<std::string> writeDiagramFile(const json& document, const std::string& name, const SolutionValues& solution,
                                          const std::string& backend, const std::string& path);

/// The most diagrams the LaTeX report includes (the loop of latex_report.cpp).
constexpr size_t REPORT_DIAGRAMS = 20;

/// The file of the report's figure of the diagram at `index` of a definitions document, as latex_report.cpp includes
/// it: "thermo_diagrams.pdf" for the first (the name of the raster image the GUI sends too), then
/// "thermo_diagram_2.pdf", "thermo_diagram_3.pdf", ...
std::string reportDiagramFile(size_t index);

/**
 * @brief The figures of the LaTeX report: every diagram of a definitions document (at most REPORT_DIAGRAMS) as a PDF
 *        in the "print" theme on the spec's page, named as latex_report.cpp includes them: "thermo_diagrams.pdf",
 *        then "thermo_diagram_2.pdf", "thermo_diagram_3.pdf", ...
 *
 * A diagram that cannot be drawn is left out with a warning naming it (the report then has no figure for it, or the
 * raster image the GUI sent); so is every one when no TrueType font is found.
 * @return (file name, PDF) pairs
 */
std::vector<std::pair<std::string, std::string>> reportDiagrams(const json& document, const SolutionValues* solution,
                                                                const std::string& backend, std::vector<std::string>& warnings);

/**
 * @brief The full state at a point of a chart (POST /api/v1/diagram/state, the hover of the Diagrams tab).
 *
 * The fluid's state at the axis values `x_SI`, `y_SI` of `diagram`: every key of stateKeys(), null for a value that
 * cannot be computed, Q = -1 in single phase as CoolProp and CoolPlot report it; for humid air every key of
 * airStateKeys() (T, W, RH, Twb, Tdp, h, v, Thx, p). One flash with the chart fluid's model (the same as the cycles
 * and the domain use), the caches of the engine kept.
 * @param fluid the chart fluid (chartFluid: the model's backend, the concentration of a solution, the total pressure
 *        of humid air)
 * @return the state, or null when the point has none: the model cannot compute it, it lies outside the temperature
 *         and pressure range of the model, or (humid air) above saturation or below dry air
 * @throws DiagramError for an unknown diagram type, one that cannot show the fluid, a value that is not a finite number
 */
json stateAt(const ChartFluid& fluid, const std::string& diagram, double x_SI, double y_SI);

/**
 * @brief The response of GET /api/v1/diagram/options: diagram types with their families and short titles ("p-h",
 *        "T-s", "psychrometric"), the fluid types drawn and what each kind of fluid offers (`fluid_kinds`), the units
 *        of the model (`model_units`), region presets, themes, label modes and the versions.
 *
 * `model_units` is the unit system of the model, in which the GUI draws every chart: {"record": the custom unit record
 * to put in a spec's "units" (coolsolveUnitSystem(): CoolSolve's fixed units today), "units": its unit table {key:
 * [label, scale, offset(, difference)]} for every quantity of the charts}.
 *
 * `fluid_kinds` maps "Real", "Mixture", "Incompressible", "IdealGas" and "HumidAir" to {"diagram_types": names
 * offered, "families": the keys of the families its lines can have, "saturation": whether it has a saturation curve
 * (dome, lines of constant quality, critical point, saturation scale), "fractions": whether it takes a concentration
 * (a mass fraction), "default_box": [T_min_K, T_max_K, p_min_Pa, p_max_Pa] of a fluid without saturation curve where
 * the spec gives no absolute limit (null: the fluid's own temperatures), "default_region": the limits of a new chart,
 * "chart_kind": "fluid" or "humid_air"}; humid air adds "total_pressure" (the default, Pa), "total_pressure_range" and
 * "temperature_range" (those of HAPropsSI, Pa and K).
 */
json options();

/// An empty diagram definitions document (§3.4).
json emptyDefinitions();
/// Check a diagram definitions document loosely (an object; "format" "coolsolve-diagrams" and a numeric
/// "version" when given; "diagrams" a list when given). @throws DiagramError
void checkDefinitions(const json& document);

}  // namespace diagram
}  // namespace coolsolve
