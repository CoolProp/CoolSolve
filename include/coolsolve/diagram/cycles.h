// The cycles of a solved model on a chart: full states from the model's variables, and the paths between them
// (isobars through the dome, isenthalps, isentropes, isotherms, polytropic paths) as chains of states.
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/cycles/paths.py (auto_kind,
// path_states, _inner_inputs, _isobar_enthalpies, cycle_paths), thermo/process.py (state_point) and core/cycle.py
// (PATH_KINDS, check_paths, default_loops, CycleSpec.to_dict). The states come from the solve instead of
// CoolPlot's StateInput (resolve_states): in CoolSolve the model is the cycle.
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
// See docs/coolplot_integration.md §4 for the cycle definitions.
#pragma once

#include "coolsolve/diagram/fluid.h"
#include "coolsolve/ir.h"  // CaseInsensitiveLess

#include <coolplot/core.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace coolsolve {
namespace diagram {

/// Values of a solved model by variable name (case-insensitive), in CoolSolve's units (°C, Pa, J/kg, ...).
using SolutionValues = std::map<std::string, double, CaseInsensitiveLess>;

/// How the path between two consecutive states is drawn (core.cycle.PATH_KINDS): "auto", "isobaric",
/// "isenthalpic", "isentropic", "isothermal", "polytropic", "straight".
const std::vector<std::string>& pathKinds();

/// Two states count as having the same property within this relative tolerance when the path kind is chosen
/// automatically (thermo.cycles.paths.SAME_RTOL).
constexpr double SAME_RTOL = 1e-6;

/// Segments per drawn path (CycleSpec.steps).
constexpr int DEFAULT_STEPS = 20;

/// The process between two states, judged from what they have in common: equal p an isobar, equal h an
/// isenthalp, equal s an isentrope, equal T an isotherm, else "polytropic" (paths.auto_kind).
std::string autoKind(const State& a, const State& b);

/**
 * @brief The states along the process from `a` to `b`, both included (paths.path_states): `steps` segments; an
 *        isobar gets the saturation states it crosses. The provider may replace the kind (PropertyProvider::pathKind:
 *        humid air draws every path straight). States CoolProp cannot compute are left out (the path has a longer
 *        straight piece there), as are states without finite T, p, h and s.
 * @throws DiagramError for an unknown kind
 */
std::vector<State> pathStates(PropertyProvider& fluid, const State& a, const State& b, const std::string& kind, int steps);

/// Connections (from, to) and their kinds; a connection that is not listed is "auto".
using PathKinds = std::map<std::pair<int, int>, std::string>;

/// For every chain, the states along all its connections, in order (paths.cycle_paths).
std::vector<std::vector<State>> cyclePaths(PropertyProvider& fluid, const std::map<int, State>& states,
                                           const std::vector<std::vector<int>>& loops, const PathKinds& paths, int steps);

/// The values of a state as a JSON record, without the values that are not finite.
json stateRecord(const State& state);

/// A value of the model (CoolSolve's units, units.h) of quantity `key` in SI units. @throws DiagramError
double solvedToSI(const std::string& key, double value);

/// CoolSolve's fixed units as a custom unit record of coolplot-core: {"name": "coolsolve", "units": {key: [label,
/// scale, offset(, difference)]}} for every quantity of a fluid chart and of a chart of humid air (W in kg/kg, RH as a
/// fraction, Twb and Tdp in °C; Thx is shown as T).
json coolsolveUnitSystem();

/**
 * @brief A cycle of a diagram definition resolved against the solved model: the CoolPlot cycle record of the spec
 *        (CycleSpec.to_dict: full states, loops, paths, steps), the states by number and the chains through them.
 */
struct ResolvedCycle {
    size_t index = 0;                                        ///< position in the request's cycle definitions
    json spec;                                               ///< the entry of spec.cycles
    std::map<int, State> states;                             ///< full states in SI units
    std::map<int, std::map<std::string, std::string>> vars;  ///< number -> {quantity key: variable name}
    std::map<int, std::string> labels;                       ///< number -> label
    std::vector<std::vector<State>> chains;                  ///< the states along each loop, in order
};

/**
 * @brief The cycles of diagram definitions (docs/coolplot_integration.md §4.1) from a solved model.
 *
 * A state is given by `vars` ({quantity key: variable name}) or, for `arrays` ({key: array name}), by the elements
 * `name[i]` of every array. The values (CoolSolve's units) are converted to SI and completed to full states from
 * two independent properties (the chart's axes when the model has both) with the chart fluid; the solved values
 * are kept. A state of humid air is at the chart's total pressure: one whose `p` variable has another value is not
 * drawn. `loops` gives the chains of state numbers; without it, one loop through the states in ascending order
 * (closed unless `closed` is false); an empty list draws the state markers alone (states that cannot be computed
 * are then left out with a warning). `paths` is a kind applied to every connection, or CoolPlot's record
 * [[[a, b], kind], ...] whose unlisted connections are "auto" (the default); `steps` the segments of a path.
 * A cycle without states is skipped silently; one that cannot be resolved (no solve, a variable without a value,
 * a state CoolProp cannot compute) is left out with a warning.
 * @throws DiagramError for definitions that are not a list of objects, two cycles of the same name, an unknown
 *         path kind or a bad `steps`
 */
std::vector<ResolvedCycle> resolveCycles(const json& definitions, const SolutionValues* solution, PropertyProvider& fluid,
                                         const coolplot::DiagramType& diagram, std::vector<std::string>& warnings);

}  // namespace diagram
}  // namespace coolsolve
