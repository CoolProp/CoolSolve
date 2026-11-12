// The cycles of a solved model on a chart (include/coolsolve/diagram/cycles.h).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/cycles/paths.py (auto_kind,
// path_states, _inner_inputs, _isobar_enthalpies, cycle_paths), thermo/process.py (state_point) and core/cycle.py
// (PATH_KINDS, check_paths, default_loops, CycleSpec.to_dict).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#include "coolsolve/diagram/cycles.h"

#include "coolsolve/diagram/planning.h"  // linearGrid, logGrid (numpy.linspace, numpy.geomspace)
#include "coolsolve/units.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <set>

namespace coolsolve {
namespace diagram {

namespace {

// ---------------------------------------------------------------------------
// CoolSolve's units (units.h: fixed, °C, Pa, J/kg, J/kg-K, kg/m^3, m^3/kg)
// ---------------------------------------------------------------------------

struct CoolSolveUnit {
    std::string key;         ///< CoolPlot quantity key
    UnitType type;
    std::string unit;        ///< CoolSolve's unit (UnitSystem)
    std::string label;       ///< the same unit spelled as coolplot-core's records write it
    std::string difference;  ///< unit of a difference when it differs from the label
};

const std::vector<CoolSolveUnit>& coolsolveUnits() {
    static const UnitSystem units;  // CoolSolve's units are fixed (units.h); the labels below spell them
    static const std::vector<CoolSolveUnit> table = {
        {"T", UnitType::Temperature, units.temperature, "deg C", "K"},
        {"p", UnitType::Pressure, units.pressure, "Pa", ""},
        {"h", UnitType::SpecificEnergy, units.specific_energy, "J/kg", ""},
        {"s", UnitType::SpecificEntropy, units.specific_entropy, "J/kg/K", ""},
        {"rho", UnitType::Density, units.density, "kg/m3", ""},
        {"v", UnitType::SpecificVolume, units.specific_volume, "m3/kg", ""},
        {"u", UnitType::SpecificEnergy, units.specific_energy, "J/kg", ""},
        {"Q", UnitType::Dimensionless, "", "-", ""},
    };
    return table;
}

/// The full state from the solved values (SI): two independent properties, the chart's axes first, (T, p) last
/// (not independent inside the dome); the solved values are kept. Empty with `error` set when none works.
std::optional<State> completeState(FluidModel& fluid, const State& given, const coolplot::DiagramType& diagram, std::string& error) {
    std::vector<std::pair<std::string, std::string>> pairs;
    auto has = [&given](const std::string& key) { return given.count(key) > 0; };
    auto add = [&pairs](const std::string& a, const std::string& b) {
        for (const auto& [x, y] : pairs)
            if ((x == a && y == b) || (x == b && y == a)) return;
        pairs.push_back({a, b});
    };
    auto isTp = [](const std::string& a, const std::string& b) { return (a == "T" && b == "p") || (a == "p" && b == "T"); };
    if (has(diagram.x) && has(diagram.y) && !isTp(diagram.x, diagram.y)) add(diagram.x, diagram.y);
    const std::vector<std::string>& keys = stateKeys();
    for (size_t i = 0; i < keys.size(); ++i)
        for (size_t j = i + 1; j < keys.size(); ++j)
            if (has(keys[i]) && has(keys[j]) && !isTp(keys[i], keys[j]) && !(keys[i] == "rho" && keys[j] == "v"))
                add(keys[i], keys[j]);
    if (has("T") && has("p")) add("T", "p");
    for (const auto& [a, b] : pairs) {
        try {
            State state = fluid.statePoint(a, given.at(a), b, given.at(b));
            for (const auto& [key, value] : given) state[key] = value;  // the marker sits where the solver put it
            State finite;
            for (const auto& [key, value] : state)
                if (std::isfinite(value)) finite[key] = value;
            return finite;
        } catch (const std::exception& e) {
            error = e.what();
        }
    }
    if (pairs.empty()) error = "no pair of independent properties";
    return std::nullopt;
}

bool samePair(const State& a, const State& b, const char* key) {
    const double x = a.at(key), y = b.at(key);
    return std::fabs(x - y) <= SAME_RTOL * std::max({std::fabs(x), std::fabs(y), 1e-12});
}

/// numpy.linspace(lo, hi, n) without its ends: the n - 2 inner values.
std::vector<double> inner(const std::vector<double>& values) {
    return values.size() <= 2 ? std::vector<double>() : std::vector<double>(values.begin() + 1, values.end() - 1);
}

/// Values from lo to hi, spaced logarithmically when both are positive (paths._log_or_linear).
std::vector<double> logOrLinear(double lo, double hi, int n) {
    return lo > 0.0 && hi > 0.0 ? logGrid(lo, hi, n) : linearGrid(lo, hi, n);
}

/// Enthalpies along an isobar from a to b, with the saturation states added so that the path bends exactly on the
/// dome (paths._isobar_enthalpies).
std::vector<double> isobarEnthalpies(FluidModel& fluid, const State& a, const State& b, const std::vector<double>& t) {
    std::vector<double> h;
    for (double ti : t) h.push_back(a.at("h") + ti * (b.at("h") - a.at("h")));
    for (double q : {0.0, 1.0}) {
        double h_sat = 0.0;
        try {
            h_sat = fluid.statePoint("p", a.at("p"), "Q", q).at("h");
        } catch (const std::exception&) {
            continue;
        }
        if (std::min(a.at("h"), b.at("h")) < h_sat && h_sat < std::max(a.at("h"), b.at("h"))) h.push_back(h_sat);
    }
    if (b.at("h") < a.at("h")) std::sort(h.begin(), h.end(), std::greater<double>());
    else std::sort(h.begin(), h.end());
    return h;
}

struct Inputs {
    std::string k1;
    double v1;
    std::string k2;
    double v2;
};

/// The inputs of the states between a and b along a path kind (paths._inner_inputs).
std::vector<Inputs> innerInputs(FluidModel& fluid, const State& a, const State& b, const std::string& kind, int steps) {
    const std::vector<double> t = inner(linearGrid(0.0, 1.0, steps + 1));
    std::vector<Inputs> out;
    if (kind == "isobaric") {
        for (double h : isobarEnthalpies(fluid, a, b, t)) out.push_back({"p", a.at("p"), "h", h});
        return out;
    }
    if (kind == "isothermal") {
        for (double ti : t) out.push_back({"T", a.at("T"), "s", a.at("s") + ti * (b.at("s") - a.at("s"))});
        return out;
    }
    if (kind != "isenthalpic" && kind != "isentropic" && kind != "polytropic")
        throw DiagramError("Unknown path kind '" + kind + "'.");
    const std::vector<double> p = inner(logOrLinear(a.at("p"), b.at("p"), steps + 1));
    for (size_t i = 0; i < p.size(); ++i) {
        if (kind == "isenthalpic") {
            out.push_back({"p", p[i], "h", a.at("h")});
        } else {
            // Entropy kept (isentropic), or changed evenly with pressure (polytropic)
            const double s = kind == "isentropic" ? a.at("s") : a.at("s") + t[i] * (b.at("s") - a.at("s"));
            out.push_back({"p", p[i], "s", s});
        }
    }
    return out;
}

/// A path kind of a definition, checked (core.cycle.check_paths).
std::string checkedKind(const std::string& kind, const std::string& where) {
    const std::vector<std::string>& kinds = pathKinds();
    if (std::find(kinds.begin(), kinds.end(), kind) != kinds.end()) return kind;
    std::string expected;
    for (const std::string& k : kinds) expected += (expected.empty() ? "'" : ", '") + k + "'";
    throw DiagramError("Unknown path '" + kind + "'" + where + ", expected one of (" + expected + ").");
}

}  // namespace

// ---------------------------------------------------------------------------
// Paths (thermo/cycles/paths.py)
// ---------------------------------------------------------------------------

const std::vector<std::string>& pathKinds() {
    static const std::vector<std::string> kinds = {"auto", "isobaric", "isenthalpic", "isentropic", "isothermal", "polytropic", "straight"};
    return kinds;
}

std::string autoKind(const State& a, const State& b) {
    for (const auto& [key, kind] : {std::pair<const char*, const char*>{"p", "isobaric"}, {"h", "isenthalpic"},
                                    {"s", "isentropic"}, {"T", "isothermal"}})
        if (samePair(a, b, key)) return kind;
    return "polytropic";
}

std::vector<State> pathStates(FluidModel& fluid, const State& a, const State& b, const std::string& kind, int steps) {
    const std::string resolved = kind == "auto" ? autoKind(a, b) : kind;
    if (resolved == "straight" || steps < 2) return {a, b};
    std::vector<State> out = {a};
    for (const Inputs& in : innerInputs(fluid, a, b, resolved, steps)) {
        State state;
        try {
            state = fluid.statePoint(in.k1, in.v1, in.k2, in.v2);
        } catch (const std::exception&) {
            continue;  // a point CoolProp cannot compute is left out
        }
        if (std::isfinite(state.at("T")) && std::isfinite(state.at("p")) && std::isfinite(state.at("h")) && std::isfinite(state.at("s")))
            out.push_back(state);
    }
    out.push_back(b);
    return out;
}

std::vector<std::vector<State>> cyclePaths(FluidModel& fluid, const std::map<int, State>& states,
                                           const std::vector<std::vector<int>>& loops, const PathKinds& paths, int steps) {
    std::vector<std::vector<State>> out;
    for (const std::vector<int>& chain : loops) {
        std::vector<State> line = {states.at(chain.front())};
        for (size_t k = 0; k + 1 < chain.size(); ++k) {
            const auto given = paths.find({chain[k], chain[k + 1]});
            const std::string kind = given == paths.end() ? "auto" : given->second;
            const std::vector<State> path = pathStates(fluid, states.at(chain[k]), states.at(chain[k + 1]), kind, steps);
            line.insert(line.end(), path.begin() + 1, path.end());
        }
        out.push_back(std::move(line));
    }
    return out;
}

// ---------------------------------------------------------------------------
// The model's values and units
// ---------------------------------------------------------------------------

json stateRecord(const State& state) {
    json out = json::object();
    for (const auto& [key, value] : state)
        if (std::isfinite(value)) out[key] = value;
    return out;
}

double solvedToSI(const std::string& key, double value) {
    for (const CoolSolveUnit& u : coolsolveUnits())
        if (u.key == key) return UnitConverter::toSI(value, u.type, u.unit);
    throw DiagramError("'" + key + "' is not a quantity of a fluid state.");
}

json coolsolveUnitSystem() {
    json units = json::object();
    for (const CoolSolveUnit& u : coolsolveUnits()) {
        // display = SI * scale + offset, from CoolSolve's own conversions
        const double offset = u.type == UnitType::Dimensionless ? 0.0 : UnitConverter::fromSI(0.0, u.type, u.unit);
        const double scale = UnitConverter::scaleFromSI(u.type, u.unit);
        json entry = json::array({u.label, scale, offset});
        if (!u.difference.empty()) entry.push_back(u.difference);
        units[u.key] = entry;
    }
    return {{"name", "coolsolve"}, {"units", units}};
}

// ---------------------------------------------------------------------------
// Cycle definitions (§4.1)
// ---------------------------------------------------------------------------

std::vector<ResolvedCycle> resolveCycles(const json& definitions, const SolutionValues* solution, FluidModel& fluid,
                                         const coolplot::DiagramType& diagram, std::vector<std::string>& warnings) {
    std::vector<ResolvedCycle> resolved;
    if (definitions.is_null()) return resolved;
    if (!definitions.is_array()) throw DiagramError("The cycles of a diagram are a list of cycle definitions.");
    std::set<std::string> names;
    for (size_t index = 0; index < definitions.size(); ++index) {
        const json& def = definitions[index];
        if (!def.is_object()) throw DiagramError("A cycle definition is an object, not " + def.dump() + ".");
        const json nameField = def.value("name", json("Cycle " + std::to_string(index + 1)));
        if (!nameField.is_string()) throw DiagramError("The name of a cycle is a text, not " + nameField.dump() + ".");
        const std::string name = nameField.get<std::string>();
        if (!names.insert(name).second) throw DiagramError("Cycle names have to be unique, got '" + name + "' twice.");
        const std::string prefix = "Cycle '" + name + "'";
        auto skip = [&](const std::string& why) { warnings.push_back(prefix + " is not drawn: " + why); };
        try {
            ResolvedCycle cycle;
            cycle.index = index;
            // How the states are joined: a kind for every connection, or CoolPlot's record of the connections that are
            // not "auto"; the segments of a path (CycleSpec.steps)
            const json pathsField = def.value("paths", json("auto"));
            PathKinds paths;
            std::string everyPath;
            if (pathsField.is_string()) {
                everyPath = checkedKind(pathsField.get<std::string>(), "");
            } else if (pathsField.is_array()) {
                for (const json& entry : pathsField) {
                    if (!entry.is_array() || entry.size() != 2 || !entry.at(0).is_array() || entry.at(0).size() != 2)
                        throw DiagramError(prefix + ": a path is [[from, to], kind], not " + entry.dump() + ".");
                    const int a = entry.at(0).at(0).get<int>(), b = entry.at(0).at(1).get<int>();
                    paths[{a, b}] = checkedKind(entry.at(1).get<std::string>(),
                                                " from " + std::to_string(a) + " to " + std::to_string(b));
                }
            } else {
                throw DiagramError(prefix + ": \"paths\" is a path kind or [[[from, to], kind], ...], not " + pathsField.dump() + ".");
            }
            const json stepsField = def.value("steps", json(DEFAULT_STEPS));
            if (!stepsField.is_number_integer() || stepsField.get<int>() < 1)
                throw DiagramError(prefix + ": steps has to be at least 1, not " + stepsField.dump() + ".");
            const int steps = stepsField.get<int>();

            // The variables of each state: explicit states, or the elements of arrays
            const json states = def.value("states", json::array());
            const json arrays = def.value("arrays", json::object());
            if (!states.is_array()) throw DiagramError(prefix + ": \"states\" is a list of states.");
            if (!arrays.is_object()) throw DiagramError(prefix + ": \"arrays\" is {quantity key: array name}.");
            if (!states.empty()) {
                for (size_t j = 0; j < states.size(); ++j) {
                    const json& s = states[j];
                    const int number = s.value("number", static_cast<int>(j) + 1);
                    if (cycle.vars.count(number)) throw DiagramError(prefix + ": state " + std::to_string(number) + " is given twice.");
                    std::map<std::string, std::string>& vars = cycle.vars[number];
                    const json given = s.value("vars", json::object());
                    for (auto it = given.begin(); it != given.end(); ++it)
                        vars[coolplot::quantity_key(it.key())] = it.value().get<std::string>();
                    const std::string label = s.value("label", std::string());
                    if (label.find_first_not_of(" \t") != std::string::npos) cycle.labels[number] = label;
                }
            } else if (!arrays.empty()) {
                if (!solution) {
                    skip("there is no solved model (solve the model first).");
                    continue;
                }
                // State i where every array has an element i
                std::optional<std::set<int>> common;
                std::map<std::string, std::string> keyArray;
                for (auto it = arrays.begin(); it != arrays.end(); ++it) {
                    const std::string key = coolplot::quantity_key(it.key()), array = it.value().get<std::string>();
                    keyArray[key] = array;
                    std::set<int> indices;
                    const std::string start = array + "[";
                    for (auto v = solution->lower_bound(start); v != solution->end(); ++v) {
                        const std::string& var = v->first;
                        if (var.size() <= start.size() || CaseInsensitiveLess{}(start, var.substr(0, start.size())) ||
                            CaseInsensitiveLess{}(var.substr(0, start.size()), start))
                            break;
                        const std::string inside = var.substr(start.size(), var.size() - start.size() - 1);
                        if (var.back() != ']' || inside.empty() || inside.size() > 9 ||
                            inside.find_first_not_of("0123456789") != std::string::npos)
                            continue;  // not an element of a one-dimensional array
                        indices.insert(std::stoi(inside));
                    }
                    if (!common) {
                        common = indices;
                    } else {
                        std::set<int> both;
                        std::set_intersection(common->begin(), common->end(), indices.begin(), indices.end(), std::inserter(both, both.end()));
                        common = both;
                    }
                }
                if (!common || common->empty()) {
                    skip("its arrays have no common elements in the last solve.");
                    continue;
                }
                for (int i : *common)
                    for (const auto& [key, array] : keyArray) cycle.vars[i][key] = array + "[" + std::to_string(i) + "]";
            } else {
                continue;  // a cycle without states yet (the state editor stores it so): nothing to draw
            }
            if (!solution) {
                skip("there is no solved model (solve the model first).");
                continue;
            }
            // The full states. A cycle drawn as markers alone ("loops": []) leaves out a state it cannot compute (a
            // group of inferred variables can be a saturation T and p, which fix no state); a cycle with chains is not drawn
            const bool markersOnly = def.contains("loops") && def.at("loops").is_array() && def.at("loops").empty();
            std::string problem;
            for (const auto& [number, vars] : cycle.vars) {
                std::string stateProblem;
                State given;
                for (const auto& [key, var] : vars) {
                    const auto found = solution->find(var);
                    if (found == solution->end()) {
                        stateProblem = "the variable '" + var + "' has no value in the last solve.";
                        break;
                    }
                    given[key] = solvedToSI(key, found->second);
                }
                if (stateProblem.empty() && given.size() < 2)
                    stateProblem = "state " + std::to_string(number) + " needs two properties, it has " + std::to_string(given.size()) + ".";
                std::optional<State> state;
                if (stateProblem.empty()) {
                    std::string error;
                    state = completeState(fluid, given, diagram, error);
                    if (!state) stateProblem = "state " + std::to_string(number) + " cannot be computed with " + fluid.name() + ": " + error;
                }
                if (stateProblem.empty()) {
                    cycle.states[number] = *state;
                } else if (markersOnly) {
                    warnings.push_back(prefix + ": " + stateProblem + (stateProblem.back() == '.' ? "" : ".") + " It is left out.");
                    cycle.labels.erase(number);
                } else {
                    problem = stateProblem;
                    break;
                }
            }
            if (problem.empty() && cycle.states.empty()) problem = "none of its states can be computed.";
            if (!problem.empty()) {
                skip(problem);
                continue;
            }
            // The chains: given (an empty list: the state markers alone), or without "loops" one loop through the
            // states in ascending order, closed unless "closed" is false (CoolPlot's default_loops)
            std::vector<std::vector<int>> loops;
            if (def.contains("loops") && !def.at("loops").is_null()) {
                loops = def.at("loops").get<std::vector<std::vector<int>>>();
            } else {
                std::vector<int> loop;
                for (const auto& [number, _] : cycle.states) loop.push_back(number);
                if (def.value("closed", true) || loop.size() < 2) loop.push_back(loop.front());
                loops.push_back(loop);
            }
            for (const std::vector<int>& loop : loops) {
                if (loop.size() < 2) problem = "a chain needs at least two states.";
                for (int number : loop)
                    if (!cycle.states.count(number)) problem = "a chain refers to state " + std::to_string(number) + ", which is not defined.";
            }
            if (!problem.empty()) {
                skip(problem);
                continue;
            }
            // One kind for every connection: written out as CoolPlot's record, "auto" being its default
            if (!everyPath.empty() && everyPath != "auto")
                for (const std::vector<int>& loop : loops)
                    for (size_t k = 0; k + 1 < loop.size(); ++k) paths[{loop[k], loop[k + 1]}] = everyPath;
            cycle.chains = cyclePaths(fluid, cycle.states, loops, paths, steps);

            // The cycle as CoolPlot's spec writes it (CycleSpec.to_dict): complete states, the loops, the paths that
            // are not "auto"
            json loopsRecord = json::array(), pathsRecord = json::array();
            for (const std::vector<int>& loop : loops) loopsRecord.push_back(loop);
            if (loops.empty()) {
                // Markers only: coolplot-core draws the states of a cycle without chains, but CoolPlot's CycleSpec reads
                // no loops as its default loop; a straight chain from each state to itself draws nothing there
                paths.clear();
                for (const auto& [number, _] : cycle.states) {
                    loopsRecord.push_back(json::array({number, number}));
                    paths[{number, number}] = "straight";
                }
            }
            for (const auto& [pair, kind] : paths) pathsRecord.push_back(json::array({json::array({pair.first, pair.second}), kind}));
            json statesRecord = json::array();
            for (const auto& [number, state] : cycle.states)
                statesRecord.push_back(json::array({number, json{{"state", stateRecord(state)}}}));
            cycle.spec = {{"name", name},
                          {"states", statesRecord},
                          {"loops", loopsRecord},
                          {"paths", pathsRecord},
                          {"numbers", def.value("numbers", true)},
                          {"steps", steps},
                          {"line_style", def.value("line_style", json())},
                          {"marker_style", def.value("marker_style", json())}};
            if (!cycle.labels.empty()) {
                json labels = json::array();
                for (const auto& [number, text] : cycle.labels) labels.push_back(json::array({number, text}));
                cycle.spec["labels"] = labels;
            }
            resolved.push_back(std::move(cycle));
        } catch (const coolplot::Error& e) {
            throw DiagramError(prefix + ": " + e.what());
        } catch (const json::exception& e) {
            throw DiagramError(prefix + " is not a valid cycle definition: " + e.what());
        }
    }
    return resolved;
}

}  // namespace diagram
}  // namespace coolsolve
