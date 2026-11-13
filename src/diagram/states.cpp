// The thermodynamic states of a model, found from its property calls (include/coolsolve/diagram/states.h,
// docs/coolplot_integration.md §4.2, decision D6). CoolSolve's own: no CoolPlot module mirrored.
#include "coolsolve/diagram/states.h"

#include "coolsolve/ast.h"
#include "coolsolve/fluids.h"
#include "coolsolve/user_hints.h"  // isConcentrationArgName

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <set>

namespace coolsolve {
namespace diagram {

namespace {

std::string lowerCase(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

/// The quantity a property function gives (the names the evaluator reads, findCoolPropOutputParam); `v` is the
/// viscosity there, not a volume, and the transport properties are not state quantities. `saturation`: the
/// function gives a saturation state (T_sat, P_sat).
std::optional<std::string> functionKey(const std::string& name, bool& saturation) {
    const std::string n = lowerCase(name);
    saturation = false;
    if (n == "temperature" || n == "t") return "T";
    if (n == "pressure" || n == "p") return "p";
    if (n == "enthalpy" || n == "h") return "h";
    if (n == "entropy" || n == "s") return "s";
    if (n == "internalenergy" || n == "intenergy" || n == "u") return "u";
    if (n == "volume") return "v";
    if (n == "density" || n == "rho" || n == "d") return "rho";
    if (n == "quality" || n == "x" || n == "q") return "Q";
    saturation = true;
    if (n == "t_sat" || n == "tsat") return "T";
    if (n == "p_sat" || n == "psat") return "p";
    saturation = false;
    return std::nullopt;
}

/// The quantity of a named argument of a property call (variable_inference.cpp, getPropMapping): `v` is a volume.
std::optional<std::string> argumentKey(const std::string& name) {
    const std::string n = lowerCase(name);
    if (n == "temperature" || n == "t") return "T";
    if (n == "pressure" || n == "p") return "p";
    if (n == "enthalpy" || n == "h") return "h";
    if (n == "entropy" || n == "s") return "s";
    if (n == "internalenergy" || n == "u") return "u";
    if (n == "volume" || n == "v") return "v";
    if (n == "density" || n == "rho" || n == "d") return "rho";
    if (n == "quality" || n == "x" || n == "q") return "Q";
    return std::nullopt;
}

/// The CoolSolve name of the fluid of a property call: its first argument, a string or a fluid name (a string
/// variable with a literal value gives that value). Empty when it is not a known fluid.
std::string callFluid(const FunctionCall& call, const IR& ir) {
    if (call.args.empty() || !call.args[0]) return "";
    std::string name;
    const ExprPtr& first = call.args[0];
    if (first->is<StringLiteral>()) {
        name = first->as<StringLiteral>().value;
    } else if (first->is<Variable>()) {
        name = first->as<Variable>().flattenedName();
        const VariableInfo* info = ir.getVariable(name);
        if (info && info->solutionStringValue) name = *info->solutionStringValue;
    } else {
        return "";
    }
    const std::shared_ptr<Fluid> fluid = FluidRegistry::getFluid(name);
    return fluid ? fluid->getName() : "";
}

void visitCalls(const ExprPtr& expr, const std::function<void(const FunctionCall&)>& visit) {
    if (!expr) return;
    std::visit(
        [&](const auto& node) {
            using T = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<T, FunctionCall>) {
                visit(node);
                for (const auto& arg : node.args) visitCalls(arg, visit);
                for (const auto& [_, arg] : node.namedArgs) visitCalls(arg, visit);
            } else if constexpr (std::is_same_v<T, UnaryOp>) {
                visitCalls(node.operand, visit);
            } else if constexpr (std::is_same_v<T, BinaryOp>) {
                visitCalls(node.left, visit);
                visitCalls(node.right, visit);
            }
        },
        expr->node);
}

/// The call as a set of state variables, `output` the variable the equation gives it to (null for none).
std::optional<PropertyCall> stateCall(const FunctionCall& call, const Variable* output, const IR& ir) {
    bool saturation = false;
    const std::optional<std::string> key = functionKey(call.name, saturation);
    if (!key) return std::nullopt;
    PropertyCall out;
    out.fluid = callFluid(call, ir);
    if (out.fluid.empty()) return std::nullopt;  // a FUNCTION of the model with the name of a property, or a fluid unknown
    out.quality = saturation || *key == "Q";
    auto add = [&out](const std::string& variable, const std::string& k) {
        for (const auto& [v, _] : out.vars)
            if (!CaseInsensitiveLess{}(v, variable) && !CaseInsensitiveLess{}(variable, v)) return;
        out.vars.push_back({variable, k});
    };
    for (const auto& [name, arg] : call.namedArgs) {
        const std::optional<std::string> argKey = argumentKey(name);
        if (!argKey || !arg) continue;
        if (*argKey == "Q") out.quality = true;
        if (arg->is<Variable>()) add(arg->as<Variable>().flattenedName(), *argKey);
    }
    if (output) add(output->flattenedName(), *key);
    return out;
}

struct UnionFind {
    std::vector<size_t> parent;
    explicit UnionFind(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    size_t find(size_t i) { return parent[i] == i ? i : parent[i] = find(parent[i]); }
    void join(size_t a, size_t b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[std::max(a, b)] = std::min(a, b);  // the earlier call names the group
    }
};

bool hasIndependentPair(const std::map<std::string, std::string>& vars, bool idealGas) {
    for (auto a = vars.begin(); a != vars.end(); ++a)
        for (auto b = std::next(a); b != vars.end(); ++b)
            if (independentKeys(a->first, b->first, idealGas)) return true;
    return false;
}

/// Whether CoolSolve's fluid of this name is an ideal gas.
bool isIdealGas(const std::string& fluid) {
    const std::shared_ptr<Fluid> found = FluidRegistry::getFluid(fluid);
    return found && found->getType() == FluidType::IdealGas;
}

/// The name of a variable without its leading letters and one underscore, as the name heuristic groups it
/// ("T_su_cp" -> "su_cp", "P1" -> "1"); the name itself when nothing is left.
std::string nameSuffix(const std::string& name) {
    size_t i = 0;
    while (i < name.size() && std::isalpha(static_cast<unsigned char>(name[i]))) ++i;
    if (i == 0) return name;
    if (i < name.size() && name[i] == '_') ++i;
    return i < name.size() ? name.substr(i) : name;
}

/// An integer written without leading zeros, as the name heuristic takes a state number.
bool isStateNumber(const std::string& text) {
    if (text.empty() || text.size() > 9 || text.find_first_not_of("0123456789") != std::string::npos) return false;
    return text == "0" || text[0] != '0';
}

/// The index of an element X[i] of a one-dimensional array, or -1.
int arrayIndex(const std::string& name) {
    const size_t open = name.find('[');
    if (open == std::string::npos || name.back() != ']') return -1;
    const std::string inside = name.substr(open + 1, name.size() - open - 2);
    return isStateNumber(inside) ? std::stoi(inside) : -1;
}

std::vector<std::string> splitUnderscore(const std::string& text) {
    std::vector<std::string> parts;
    size_t start = 0;
    for (size_t pos; (pos = text.find('_', start)) != std::string::npos; start = pos + 1) parts.push_back(text.substr(start, pos - start));
    parts.push_back(text.substr(start));
    return parts;
}

}  // namespace

bool independentKeys(const std::string& a, const std::string& b, bool idealGas) {
    if (a == b) return false;
    auto is = [&](const char* x, const char* y) { return (a == x && b == y) || (a == y && b == x); };
    if (idealGas && (is("T", "h") || is("T", "u") || is("h", "u"))) return false;  // each fixes only the temperature
    return !is("T", "p") && !is("rho", "v");
}

namespace {

/// Every function call of the equations, with the output variable of the call that gives an equation its value
/// (`h_2 = enthalpy(...)`, or the other way round) and null for the others; the inputs of a CALL too.
void visitEquationCalls(const IR& ir, const std::function<void(const FunctionCall&, const Variable*)>& visit) {
    for (const EquationInfo& eq : ir.getEquations()) {
        if (!eq.lhs || !eq.rhs) {
            if (eq.procedureCall)  // the inputs of a CALL can hold property calls, without an output variable
                for (const ExprPtr& arg : eq.procedureCall->inputArgs)
                    visitCalls(arg, [&](const FunctionCall& call) { visit(call, nullptr); });
            continue;
        }
        const FunctionCall* top = nullptr;
        const Variable* output = nullptr;
        if (eq.lhs->is<Variable>() && eq.rhs->is<FunctionCall>()) {
            output = &eq.lhs->as<Variable>();
            top = &eq.rhs->as<FunctionCall>();
        } else if (eq.rhs->is<Variable>() && eq.lhs->is<FunctionCall>()) {
            output = &eq.rhs->as<Variable>();
            top = &eq.lhs->as<FunctionCall>();
        }
        for (const ExprPtr& side : {eq.lhs, eq.rhs})
            visitCalls(side, [&](const FunctionCall& call) { visit(call, &call == top ? output : nullptr); });
    }
}

}  // namespace

std::vector<PropertyCall> propertyCalls(const IR& ir) {
    std::vector<PropertyCall> calls;
    visitEquationCalls(ir, [&](const FunctionCall& call, const Variable* output) {
        if (auto c = stateCall(call, output, ir)) calls.push_back(*c);
    });
    return calls;
}

std::map<std::string, double> solutionFractions(const IR& ir, const SolutionValues* solution) {
    std::map<std::string, double> fractions;
    visitEquationCalls(ir, [&](const FunctionCall& call, const Variable*) {
        const std::string fluid = callFluid(call, ir);
        if (fluid.empty() || fractions.count(fluid)) return;
        const auto incompressible = std::dynamic_pointer_cast<IncompressibleFluid>(FluidRegistry::getFluid(fluid));
        if (!incompressible || !incompressible->isSolution()) return;
        // The concentration in mass %, as the evaluator reads it: C=..., or the second positional argument
        const Expression* given = nullptr;
        for (const auto& [name, arg] : call.namedArgs)
            if (isConcentrationArgName(name) && arg) given = arg.get();
        if (!given && call.args.size() >= 2 && call.args[1]) given = call.args[1].get();
        double percent = std::numeric_limits<double>::quiet_NaN();
        if (given && given->is<NumberLiteral>()) {
            percent = given->as<NumberLiteral>().value;
        } else if (given && given->is<Variable>() && solution) {
            const auto found = solution->find(given->as<Variable>().flattenedName());
            if (found != solution->end()) percent = found->second;
        }
        if (std::isfinite(percent)) fractions[fluid] = percent / 100.0;
    });
    return fractions;
}

ModelStates groupStates(const IR& ir, const SolutionValues* solution, bool mergeOnTp) {
    std::vector<PropertyCall> calls = propertyCalls(ir);
    if (solution)  // the variables of the solved model only
        for (PropertyCall& call : calls)
            call.vars.erase(std::remove_if(call.vars.begin(), call.vars.end(),
                                           [solution](const auto& var) { return !solution->count(var.first); }),
                            call.vars.end());
    // Calls that share a variable, through an index of the variables
    std::map<std::string, std::vector<size_t>, CaseInsensitiveLess> callsOf;
    for (size_t i = 0; i < calls.size(); ++i)
        for (const auto& [var, _] : calls[i].vars) callsOf[var].push_back(i);
    UnionFind groups(calls.size());
    for (size_t i = 0; i < calls.size(); ++i) {
        const bool idealGas = isIdealGas(calls[i].fluid);
        std::set<size_t> others;
        for (const auto& [var, _] : calls[i].vars)
            for (size_t j : callsOf[var])
                if (j > i && calls[j].fluid == calls[i].fluid) others.insert(j);
        for (size_t j : others) {
            std::set<std::string> sharedKeys;  // the quantities of the variables both calls have
            for (const auto& [var, key] : calls[i].vars)
                for (const auto& [other, _] : calls[j].vars)
                    if (!CaseInsensitiveLess{}(var, other) && !CaseInsensitiveLess{}(other, var)) sharedKeys.insert(key);
            bool same = false;
            for (auto a = sharedKeys.begin(); a != sharedKeys.end() && !same; ++a)
                for (auto b = std::next(a); b != sharedKeys.end() && !same; ++b) same = independentKeys(*a, *b, idealGas);
            if (!same && mergeOnTp && sharedKeys.count("T") && sharedKeys.count("p") && !calls[i].quality && !calls[j].quality)
                same = true;  // the refinement: two single-phase states at the same T and p
            if (same) groups.join(i, j);
        }
    }
    // The variables of each group, the first one of each key, in the order of the calls
    std::map<std::string, std::vector<std::map<std::string, std::string>>> byFluid;
    std::map<size_t, std::pair<std::string, size_t>> slot;  // root call -> (fluid, position in byFluid)
    for (size_t i = 0; i < calls.size(); ++i) {
        const size_t root = groups.find(i);
        auto found = slot.find(root);
        if (found == slot.end()) {
            byFluid[calls[i].fluid].push_back({});
            found = slot.emplace(root, std::make_pair(calls[i].fluid, byFluid[calls[i].fluid].size() - 1)).first;
        }
        std::map<std::string, std::string>& vars = byFluid[found->second.first][found->second.second];
        for (const auto& [var, key] : calls[i].vars) vars.emplace(key, var);
    }
    ModelStates states;
    for (auto& [fluid, groupsOfFluid] : byFluid) {
        std::vector<std::map<std::string, std::string>> kept;
        for (auto& vars : groupsOfFluid)
            if (hasIndependentPair(vars, isIdealGas(fluid))) kept.push_back(std::move(vars));
        if (!kept.empty()) states[fluid] = numberStates(kept);
    }
    return states;
}

ModelStates nameHeuristicStates(const IR& ir, const SolutionValues* solution) {
    // The property codes of the variable inference CoolPlot knows (definition.ts, PROPERTY_KEYS)
    static const std::map<std::string, std::string> keys = {{"T", "T"}, {"P", "p"}, {"H", "h"}, {"S", "s"},
                                                            {"D", "rho"}, {"U", "u"}, {"Q", "Q"}};
    std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> groups;  // fluid -> name -> vars
    for (const auto& [name, info] : ir.getVariables()) {
        if (name.find('[') != std::string::npos || info.inferredFluid.empty()) continue;
        if (solution && !solution->count(name)) continue;
        const auto key = keys.find(info.inferredProperty);
        if (key == keys.end()) continue;
        groups[info.inferredFluid][nameSuffix(name)][key->second] = name;
    }
    ModelStates states;
    for (const auto& [fluid, byName] : groups) {
        std::vector<ModelState> numbered, named;
        for (const auto& [name, vars] : byName) {
            if (vars.size() < 2) continue;
            if (isStateNumber(name)) numbered.push_back({std::stoi(name), "", vars});
            else named.push_back({0, name, vars});
        }
        std::set<int> used;
        for (const ModelState& s : numbered) used.insert(s.number);
        int next = 1;
        for (ModelState& s : named) {  // sorted by name: the order of the map
            while (used.count(next)) ++next;
            used.insert(next);
            s.number = next;
            numbered.push_back(s);
        }
        std::sort(numbered.begin(), numbered.end(), [](const ModelState& a, const ModelState& b) { return a.number < b.number; });
        if (!numbered.empty()) states[fluid] = numbered;
    }
    return states;
}

namespace {

/// The state number a variable's name carries: i for an element X[i], else the integer that ends the name after
/// its leading letters (T_2, P2, x_4 -> 2, 2, 4; h_2_s and T_su_cp carry none); 0 for none.
int nameNumber(const std::string& var) {
    const int index = arrayIndex(var);
    if (index >= 0) return index;
    const std::vector<std::string> parts = splitUnderscore(nameSuffix(var));
    return !parts.empty() && isStateNumber(parts.back()) ? std::stoi(parts.back()) : 0;
}

/// The number most variables of a group carry, when more than half of them do (T_2, h_2, P_2, s[2] -> 2;
/// T[1], h[1], x_1 -> 1; h_2_s, P_2, s_1 -> none: 0). `share` receives the fraction of the variables carrying it.
int agreedNumber(const std::map<std::string, std::string>& vars, double& share) {
    std::map<int, size_t> votes;
    for (const auto& [_, var] : vars)
        if (const int n = nameNumber(var); n > 0) ++votes[n];
    int best = 0;
    size_t count = 0;
    for (const auto& [n, c] : votes)
        if (c > count) best = n, count = c;
    share = vars.empty() ? 0.0 : static_cast<double>(count) / static_cast<double>(vars.size());
    return 2 * count > vars.size() ? best : 0;
}

/// The underscore-separated parts that end every name after its leading letters: "su_cp" of T_su_cp, P_su_cp.
std::string commonSuffix(const std::map<std::string, std::string>& vars) {
    std::vector<std::string> common;
    bool first = true;
    for (const auto& [_, var] : vars) {
        const std::vector<std::string> parts = splitUnderscore(nameSuffix(var));
        if (first) {
            common = parts;
            first = false;
            continue;
        }
        size_t n = 0;
        while (n < common.size() && n < parts.size() && common[common.size() - 1 - n] == parts[parts.size() - 1 - n]) ++n;
        common.erase(common.begin(), common.end() - static_cast<std::ptrdiff_t>(n));
    }
    std::string label;
    for (const std::string& part : common) label += (label.empty() ? "" : "_") + part;
    return label;
}

}  // namespace

std::vector<ModelState> numberStates(const std::vector<std::map<std::string, std::string>>& groups) {
    struct Candidate {
        ModelState state;
        double share = 0.0;  ///< the fraction of its variables that carry its number
    };
    std::vector<Candidate> candidates;
    for (const auto& vars : groups) {
        Candidate c;
        c.state.vars = vars;
        c.state.number = agreedNumber(vars, c.share);
        if (c.state.number == 0) c.state.label = commonSuffix(vars);
        candidates.push_back(std::move(c));
    }
    // The numbers go to the groups whose names agree most on them; a group whose number is taken keeps it as label
    std::vector<size_t> order(candidates.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return candidates[a].share > candidates[b].share; });
    std::vector<ModelState> states;
    std::vector<ModelState> named;
    std::set<int> used;
    for (size_t i : order) {
        ModelState& state = candidates[i].state;
        if (state.number > 0 && used.insert(state.number).second) continue;
        if (state.number > 0) state.label = std::to_string(state.number);  // the number is taken
        state.number = 0;
    }
    for (Candidate& c : candidates) (c.state.number > 0 ? states : named).push_back(std::move(c.state));
    // The others take the next free numbers, by label (those without one last, in their order)
    std::stable_sort(named.begin(), named.end(), [](const ModelState& a, const ModelState& b) {
        if (a.label.empty() != b.label.empty()) return b.label.empty();
        return a.label < b.label;
    });
    int next = 1;
    for (ModelState& state : named) {
        while (used.count(next)) ++next;
        used.insert(next);
        state.number = next;
        states.push_back(state);
    }
    std::sort(states.begin(), states.end(), [](const ModelState& a, const ModelState& b) { return a.number < b.number; });
    return states;
}

json statesToJson(const ModelStates& states, const std::map<std::string, double>& fractions) {
    json fluids = json::object();
    for (const auto& [fluid, list] : states) {
        json entries = json::array();
        for (const ModelState& state : list)
            entries.push_back({{"number", state.number}, {"label", state.label.empty() ? json() : json(state.label)}, {"vars", state.vars}});
        fluids[fluid] = entries;
    }
    json out = {{"fluids", fluids}};
    if (!fractions.empty()) {
        json given = json::object();
        for (const auto& [fluid, fraction] : fractions) given[fluid] = json::array({fraction});
        out["fractions"] = given;
    }
    return out;
}

}  // namespace diagram
}  // namespace coolsolve
