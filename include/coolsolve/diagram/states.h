// The thermodynamic states of a model, found from its property calls (docs/coolplot_integration.md §4.2,
// decision D6): each call of a property function on a fluid involves a set of state variables (its inputs and the
// variable it gives), and two calls describe the same state when they share two variables of independent
// properties. T and p are not independent inside the dome: a shared (T, p) pair alone does not merge two calls.
//
// CoolSolve's own: CoolPlot is given its cycles state by state. The numbering follows the former name heuristic
// of the Diagrams tab (T_1, P_1, h_1 -> state 1), which nameHeuristicStates() keeps for comparison.
#pragma once

#include "coolsolve/diagram/cycles.h"  // SolutionValues

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace coolsolve {
namespace diagram {

/// One call of a property function on a fluid in the equations of a model, with its state variables.
struct PropertyCall {
    std::string fluid;                                        ///< CoolSolve fluid name (FluidRegistry::getName)
    std::vector<std::pair<std::string, std::string>> vars;    ///< (variable, quantity key), inputs then output
    bool quality = false;                                     ///< a vapour quality among its inputs or output (a value too)
};

/// The property calls of a model's equations, in the order of the equations: functions of a state (temperature,
/// pressure, enthalpy, entropy, intenergy, volume, density, quality, T_sat, P_sat and their short names) whose
/// first argument is a known fluid; of humid air (AirH2O) the functions and arguments of CoolPlot's humid-air
/// quantities (temperature, pressure, humrat, relhum, wetbulb, dewpoint, enthalpy, volume and their one-letter
/// names T, P, W, R, B, D, H, V; keys T, p, W, RH, Twb, Tdp, h, v). Arguments that are not plain variables are not
/// state variables.
std::vector<PropertyCall> propertyCalls(const IR& ir);

/// Whether two quantity keys fix a state of a pure fluid together: different, and neither (T, p) nor (rho, v); for
/// an ideal gas (`idealGas`) neither two of T, h and u, which each fix the temperature alone.
bool independentKeys(const std::string& a, const std::string& b, bool idealGas = false);

/// Whether two quantity keys fix a state of humid air at a total pressure together: different, neither the pressure,
/// and not the humidity ratio with the dew point (at one pressure each fixes the other).
bool independentAirKeys(const std::string& a, const std::string& b);

/// One state of the model: its number, its label (empty for none) and its variables by quantity key.
struct ModelState {
    int number = 0;
    std::string label;
    std::map<std::string, std::string> vars;  ///< CoolPlot key (T p h s rho u v Q; humid air T p W RH Twb Tdp h v) -> variable name
};

/// The states of a model by fluid.
using ModelStates = std::map<std::string, std::vector<ModelState>>;

/**
 * @brief The states of a model by fluid from its property calls (§4.2): calls on the same fluid that share two
 *        variables of independent properties are one state. Only states with two independent properties are kept
 *        (of an ideal gas: not two of T, h, u). Humid air is at the total pressure of its chart: its calls are one state
 *        when they share two variables of independent properties other than the pressure (independentAirKeys).
 * @param solution when given, variables without a value in it are left out (the solved model)
 * @param mergeOnTp also merge two calls sharing T and p when neither involves a vapour quality (nor is a saturation
 *        function): T and p are independent in a call that takes no quality, since such a call fixes a single-phase
 *        state. On by default (decided from the library validation, docs/coolplot_integration.md §4.2); false gives
 *        the strict rule, for comparison only. Not for humid air, whose T and p fix no state.
 */
ModelStates groupStates(const IR& ir, const SolutionValues* solution = nullptr, bool mergeOnTp = true);

/**
 * @brief The states of the former name heuristic of the Diagrams tab (gui/src/components/diagram/definition.ts,
 *        inferredStates): the inferred variables (inferVariables must have run on `ir`) of a fluid grouped by their
 *        name without its leading letters and underscore, groups of two properties or more. For comparison.
 */
ModelStates nameHeuristicStates(const IR& ir, const SolutionValues* solution = nullptr);

/**
 * @brief Numbers for groups of variables, from their names: the number that more than half of the variables of a
 *        group carry (i of X[i], or the integer ending the name: T_2, h_2, P_2, s[2] -> 2; T[1], x_1 -> 1), the
 *        groups whose names agree most claiming their number first; a group without one, or whose number is taken,
 *        is labelled (by the suffix common to its names, T_su_cp, h_su_cp -> "su_cp", or by its number) and takes
 *        the next free number.
 */
std::vector<ModelState> numberStates(const std::vector<std::map<std::string, std::string>>& groups);

/**
 * @brief The concentration of each incompressible solution of a model as a mass fraction (the default of its charts):
 *        that of the first call on the solution whose concentration (C=..., or the second positional argument, in
 *        mass %, as the evaluator reads it) is a number or, with `solution`, a solved variable.
 */
std::map<std::string, double> solutionFractions(const IR& ir, const SolutionValues* solution = nullptr);

/// GET /api/v1/diagram/states: {"fluids": {fluid: [{"number", "label" (null for none), "vars"}]}}, and
/// "fractions": {solution: [mass fraction]} when the model has incompressible solutions (solutionFractions).
json statesToJson(const ModelStates& states, const std::map<std::string, double>& fractions = {});

}  // namespace diagram
}  // namespace coolsolve
