#include "coolsolve/ir.h"
#include "coolsolve/parser.h"
#include "coolsolve/constants.h"
#include <nlohmann/json.hpp>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cctype>

namespace coolsolve {

// ============================================================================
// IR Building from AST
// ============================================================================

IR IR::fromAST(const Program& program) {
    IR ir;
    int eqId = 0;
    std::vector<std::string> pendingComments;  // Accumulate "" display comments
    
    for (const auto& stmt : program.statements) {
        // Collect standalone "" display comments and attach to next equation
        if (stmt->is<Comment>()) {
            const auto& cmt = stmt->as<Comment>();
            if (cmt.isDisplay) {
                pendingComments.push_back(cmt.text);
            }
            continue;
        }

        if (stmt->is<Equation>()) {
            const auto& eq = stmt->as<Equation>();
            
            EquationInfo eqInfo;
            eqInfo.id = eqId++;
            eqInfo.sourceLine = stmt->sourceLineNumber;
            eqInfo.lhs = eq.lhs;
            eqInfo.rhs = eq.rhs;
            eqInfo.originalText = astToString(eq.lhs) + " = " + astToString(eq.rhs);
            eqInfo.inlineComment = eq.comment;
            eqInfo.precedingComments = std::move(pendingComments);
            pendingComments.clear();
            
            // Extract variables from both sides
            ir.extractVariables(eq.lhs, eqInfo.variables);
            ir.extractVariables(eq.rhs, eqInfo.variables);
            // RHS only: for explicit-solve detection (if output var not in RHS, we can solve by direct evaluation)
            ir.extractVariables(eq.rhs, eqInfo.rhsVariables);
            
            // Register all variables
            for (const auto& varName : eqInfo.variables) {
                ir.addVariable(varName);
            }

            // Handle string variable assignments (e.g., fluid$ = 'Water')
            if (eq.lhs->is<Variable>() && eq.rhs->is<StringLiteral>()) {
                const auto& var = eq.lhs->as<Variable>();
                if (!var.name.empty() && var.name.back() == '$') {
                    auto it = ir.variables_.find(var.name);
                    if (it != ir.variables_.end()) {
                        it->second.solutionStringValue = eq.rhs->as<StringLiteral>().value;
                        it->second.type = VariableType::String;
                    }
                }
            }
            
            ir.equations_.push_back(eqInfo);
        } else if (stmt->is<FunctionDefinition>()) {
            ir.functions_.push_back(stmt->as<FunctionDefinition>());
        } else if (stmt->is<ProcedureDefinition>()) {
            ir.procedures_.push_back(stmt->as<ProcedureDefinition>());
        } else if (stmt->is<ProcedureCall>()) {
            // Procedure calls involve multiple variables and provide multiple equations
            const auto& call = stmt->as<ProcedureCall>();
            
            // Extract all variables involved in the call
            std::set<std::string, CaseInsensitiveLess> callVars;
            for (const auto& arg : call.inputArgs) {
                ir.extractVariables(arg, callVars);
            }
            for (const auto& var : call.outputVars) {
                ir.addVariable(var.flattenedName());
                callVars.insert(var.flattenedName());
            }
            
            // Register all variables in the IR
            for (const auto& varName : callVars) {
                ir.addVariable(varName);
            }

            // A procedure call with N outputs provides N equations to the system.
            // We create N EquationInfo objects, each representing one output.
            // This ensures the system remains square and matching works correctly.
            if (call.outputVars.empty()) {
                // Special case: CALL with no outputs (unlikely but possible)
                EquationInfo eqInfo;
                eqInfo.id = eqId++;
                eqInfo.sourceLine = stmt->sourceLineNumber;
                eqInfo.originalText = "CALL " + call.name + "(...)";
                eqInfo.procedureCall = call;
                eqInfo.variables = callVars;
                eqInfo.rhsVariables = callVars;  // no LHS/RHS; treat as implicit
                ir.equations_.push_back(eqInfo);
            } else {
                for (size_t i = 0; i < call.outputVars.size(); ++i) {
                    EquationInfo eqInfo;
                    eqInfo.id = eqId++;
                    eqInfo.sourceLine = stmt->sourceLineNumber;
                    eqInfo.originalText = "CALL " + call.name + "(...) [output " + std::to_string(i+1) + "]";
                    eqInfo.procedureCall = call;
                    eqInfo.variables = callVars;
                    eqInfo.rhsVariables = callVars;  // procedure outputs are implicit, not explicit
                    ir.equations_.push_back(eqInfo);
                }
            }
        }
    }
    
    // Build incidence matrix
    ir.incidenceMatrix_.resize(ir.equations_.size());
    for (size_t i = 0; i < ir.equations_.size(); ++i) {
        ir.incidenceMatrix_[i] = ir.equations_[i].variables;
    }
    
    // Assign indices to variables
    int varIdx = 0;
    for (auto& [name, info] : ir.variables_) {
        info.index = varIdx++;
    }
    
    return ir;
}

// Set of EES thermodynamic function names (lowercase for comparison)
// In EES, thermophysical functions have a fluid name as their first positional
// argument, which can be specified without quotes. These functions are skipped
// when extracting variables to avoid treating fluid names as system variables.
static const std::set<std::string> THERMO_FUNCTIONS = {
    // Standard thermodynamic properties
    "enthalpy", "entropy", "temperature", "pressure", "density",
    "quality", "volume", "cp", "cv", "viscosity", "conductivity",
    "t_sat", "p_sat", "t_crit", "p_crit", "v_crit", "t_triple",
    "molarmass", "soundspeed", "specheat", "intenergy",
    // Psychrometric / humid air functions
    "humrat", "wetbulb", "dewpoint", "relhum",
    // Additional EES functions
    "surfacetension", "prandtl", "thermaldiffusivity", "kinematicviscosity",
    "compressibilityfactor", "fugacity", "isentropicexponent", "volexpcoef",
    "isothermalcompress", "joulethomsoncoef", "dipole", "acentricfactor",
    "enthalpy_fusion", "enthalpy_vaporization", "normalboilingpt", "freezingpt",
    "gibbsfreeenergy", "helmholtzfreeenergy", "phase$",
    "higherheatingvalue", "lowerheatingvalue", "enthalpy_formation",
    "massfraction", "molarmass_soln", "saturation%",
    "p_melting", "p_sublimation", "t_melting"
};

// Check if a function name is a thermodynamic function
static bool isThermoFunction(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return THERMO_FUNCTIONS.find(lower) != THERMO_FUNCTIONS.end();
}

void IR::extractVariables(const ExprPtr& expr, std::set<std::string, CaseInsensitiveLess>& vars) {
    if (!expr) return;
    
    std::visit([this, &vars](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        
        if constexpr (std::is_same_v<T, Variable>) {
            // Get the flattened name for arrays
            std::string name = node.flattenedName();
            
            // String variables (ending with $) are now included
            vars.insert(name);
            
            // Also extract variables from array indices
            for (const auto& idx : node.indices) {
                extractVariables(idx, vars);
            }
        } else if constexpr (std::is_same_v<T, UnaryOp>) {
            extractVariables(node.operand, vars);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            extractVariables(node.left, vars);
            extractVariables(node.right, vars);
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            // For thermodynamic functions, skip the first positional arg (fluid name)
            bool isThermo = isThermoFunction(node.name);
            bool firstArg = true;
            
            // Extract from positional args
            for (const auto& arg : node.args) {
                if (isThermo && firstArg) {
                    // Skip first argument (fluid name) for thermo functions
                    firstArg = false;
                    continue;
                }
                extractVariables(arg, vars);
                firstArg = false;
            }
            // Extract from named args
            for (const auto& [name, arg] : node.namedArgs) {
                extractVariables(arg, vars);
            }
        }
    }, expr->node);
}

void IR::addVariable(const std::string& name) {
    if (variables_.find(name) != variables_.end()) return;
    
    VariableInfo info;
    info.originalName = name;
    info.internalId = name;  // Use original name as internal ID for now
    info.latexLabel = variableToLatex(name);
    
    if (!name.empty() && name.back() == '$') {
        info.type = VariableType::String;
    } else {
        info.type = VariableType::Real;
    }
    
    // Check for built-in constants
    if (Constants::isConstant(name)) {
        info.solutionValue = Constants::getValue(name);
        // We could also set units/description if we had fields for them in VariableInfo
        // auto constInfo = Constants::getInfo(name);
        // if (constInfo) {
        //     info.units = constInfo->units;
        //     info.description = constInfo->description;
        // }
    }
    
    variables_[name] = info;
}

const VariableInfo* IR::getVariable(const std::string& name) const {
    auto it = variables_.find(name);
    return (it != variables_.end()) ? &it->second : nullptr;
}

VariableInfo* IR::getVariableMutable(const std::string& name) {
    auto it = variables_.find(name);
    return (it != variables_.end()) ? &it->second : nullptr;
}

int IR::getAlgebraicVariableCount() const {
    int count = 0;
    for (const auto& [name, info] : variables_) {
        if (info.type == VariableType::Real && !Constants::isConstant(name)) {
            count++;
        }
    }
    return count;
}

int IR::getNonConstantVariableCount() const {
    int count = 0;
    for (const auto& [name, info] : variables_) {
        if (!Constants::isConstant(name)) {
            count++;
        }
    }
    return count;
}

std::vector<std::string> IR::getUnmatchedVariables() const {
    std::set<std::string> determined;
    for (const auto& eq : equations_) {
        if (!eq.outputVariable.empty()) {
            determined.insert(eq.outputVariable);
        }
    }
    
    std::vector<std::string> unmatched;
    for (const auto& [name, info] : variables_) {
        // Only count non-constant variables
        if (!Constants::isConstant(name)) {
            if (determined.find(name) == determined.end()) {
                unmatched.push_back(name);
            }
        }
    }
    return unmatched;
}

std::vector<int> IR::getUnmatchedEquations() const {
    std::vector<int> unmatched;
    for (const auto& eq : equations_) {
        if (eq.outputVariable.empty()) {
            unmatched.push_back(eq.id);
        }
    }
    return unmatched;
}

// ============================================================================
// Solution/Initials Loading
// ============================================================================

int IR::loadInitialsFromFile(const std::string& filePath) {
    return loadSolutionFromFile(filePath);
}

int IR::loadSolutionFromFile(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        return 0;
    }
    
    int matchCount = 0;
    std::string line;
    
    while (std::getline(file, line)) {
        // Remove comments
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        
        // Trim whitespace
        while (!line.empty() && std::isspace(line.front())) line.erase(0, 1);
        while (!line.empty() && std::isspace(line.back())) line.pop_back();
        
        if (line.empty()) continue;
        
        // Find equals sign
        size_t eqPos = line.find('=');
        if (eqPos == std::string::npos) continue;
        
        std::string varName = line.substr(0, eqPos);
        std::string valueStr = line.substr(eqPos + 1);
        
        // Trim variable name and value
        while (!varName.empty() && std::isspace(varName.front())) varName.erase(0, 1);
        while (!varName.empty() && std::isspace(varName.back())) varName.pop_back();
        
        while (!valueStr.empty() && std::isspace(valueStr.front())) valueStr.erase(0, 1);
        while (!valueStr.empty() && std::isspace(valueStr.back())) valueStr.pop_back();

        // Check for units in the value string (e.g., 123.4 "kPa")
        std::string units;
        if (valueStr.size() > 2 && valueStr.back() == '"') {
            size_t quoteStart = valueStr.rfind('"', valueStr.size() - 2);
            if (quoteStart != std::string::npos) {
                // Found unit string
                units = valueStr.substr(quoteStart + 1, valueStr.size() - quoteStart - 2);
                valueStr = valueStr.substr(0, quoteStart);
                while (!valueStr.empty() && std::isspace(valueStr.back())) valueStr.pop_back();
            }
        }
        
        // Check if variable exists
        auto it = variables_.find(varName);
        if (it == variables_.end()) continue;
        
        VariableInfo& info = it->second;

        // Set units if found
        if (!units.empty()) {
            info.units = units;
        }
        
        // Check if value is a string
        if (valueStr.size() >= 2 && valueStr.front() == '\'' && valueStr.back() == '\'') {
            info.solutionStringValue = valueStr.substr(1, valueStr.size() - 2);
            info.type = VariableType::String;
            matchCount++;
        } else {
            try {
                double val = std::stod(valueStr);
                info.solutionValue = val;
                info.type = VariableType::Real;
                matchCount++;
            } catch (...) {
                // Ignore parse errors
            }
        }
    }
    
    return matchCount;
}

// ============================================================================
// LaTeX Generation
// ============================================================================

std::string variableToLatex(const std::string& name) {
    // ---------------------------------------------------------------
    // EES-to-LaTeX variable-name translation.
    //
    // Rules (applied in this order):
    //   1. Array brackets:  P[1] → P_{1}  (recursive on base)
    //   2. Split name on underscores into segments.
    //   3. Detect a Greek-letter prefix in the first segment:
    //        - Uppercase (DELTA, OMEGA …) matched as a *prefix* inside
    //          the segment, so "DELTAW" → (\Delta, W).
    //        - Lowercase (eta, alpha …) matched only when the *entire*
    //          first segment equals the Greek word, to avoid false
    //          positives like "pipe" matching "pi".
    //   4. Segments equal to "bar" or "dot" become \bar{} / \dot{}
    //      decorators applied to the core symbol.
    //   5. All remaining segments form the subscript (comma-separated).
    //
    // Examples:
    //   eta_bar_boil_1   → \bar{\eta}_{boil,1}
    //   DELTAW_dot_cp    → \Delta \dot{W}_{cp}
    //   M_dot            → \dot{M}
    //   T_oil_su         → T_{oil,su}
    //   P[1]             → P_{1}
    //   omega            → \omega
    // ---------------------------------------------------------------

    // 0. Strip trailing $ (EES string-variable marker) — not valid in LaTeX math
    std::string cleanName = name;
    if (!cleanName.empty() && cleanName.back() == '$') {
        cleanName.pop_back();
    }
    if (cleanName.empty()) return "\\text{" + name + "}";

    // 1. Array brackets: extract index and fold into the base name
    //    so that T_hf[2] becomes segments ["T","hf","2"] (avoiding
    //    double subscripts like T_{hf}_{2}).
    std::string arrayIndex;
    size_t bk = cleanName.find('[');
    if (bk != std::string::npos) {
        arrayIndex = cleanName.substr(bk + 1);
        if (!arrayIndex.empty() && arrayIndex.back() == ']') arrayIndex.pop_back();
        cleanName = cleanName.substr(0, bk);
    }

    // 2. Split on underscores
    std::vector<std::string> segs;
    {
        std::string seg;
        for (char c : cleanName) {
            if (c == '_') {
                if (!seg.empty()) { segs.push_back(seg); seg.clear(); }
            } else {
                seg += c;
            }
        }
        if (!seg.empty()) segs.push_back(seg);
    }
    if (segs.empty()) return name;

    // 3. Greek-letter prefix detection on the first segment
    //    Uppercase prefixes are sorted longest-first to avoid partial clashes.
    static const std::vector<std::pair<std::string, std::string>> greekMap = {
        // Uppercase (prefix match inside first segment)
        {"LAMBDA", "\\Lambda"}, {"DELTA",  "\\Delta"},
        {"GAMMA",  "\\Gamma"},  {"OMEGA",  "\\Omega"},
        {"SIGMA",  "\\Sigma"},  {"THETA",  "\\Theta"},
        // Lowercase (exact first-segment match only)
        {"epsilon","\\epsilon"},{"lambda", "\\lambda"},
        {"alpha",  "\\alpha"},  {"beta",   "\\beta"},
        {"gamma",  "\\gamma"},  {"delta",  "\\delta"},
        {"theta",  "\\theta"},  {"sigma",  "\\sigma"},
        {"omega",  "\\omega"},  {"eta",    "\\eta"},
        {"mu",     "\\mu"},     {"nu",     "\\nu"},
        {"pi",     "\\pi"},     {"rho",    "\\rho"},
        {"tau",    "\\tau"},    {"phi",    "\\phi"},
    };

    std::string greekLatex;   // e.g. "\\eta"
    std::string baseSymbol;   // e.g. "" or "W"
    bool foundGreek = false;

    const std::string& first = segs[0];
    for (const auto& [gk, gl] : greekMap) {
        bool isUpper = (gk[0] >= 'A' && gk[0] <= 'Z');
        if (isUpper) {
            // Uppercase: prefix match inside first segment
            if (first.size() >= gk.size() &&
                first.substr(0, gk.size()) == gk) {
                greekLatex = gl;
                baseSymbol = first.substr(gk.size());
                foundGreek = true;
                break;
            }
        } else {
            // Lowercase: entire first segment must equal the Greek word
            if (first == gk) {
                greekLatex = gl;
                baseSymbol = "";
                foundGreek = true;
                break;
            }
        }
    }

    if (!foundGreek) {
        baseSymbol = first;
    }

    // 4. Classify remaining segments as modifiers or subscripts
    std::vector<std::string> modifiers;   // "bar", "dot"
    std::vector<std::string> subscripts;
    for (size_t i = 1; i < segs.size(); ++i) {
        if (segs[i] == "bar" || segs[i] == "dot") {
            modifiers.push_back(segs[i]);
        } else {
            subscripts.push_back(segs[i]);
        }
    }
    // Append array index as the last subscript element
    if (!arrayIndex.empty()) {
        subscripts.push_back(arrayIndex);
    }

    // 5. Build the core symbol
    //    Modifiers wrap only the base (or the Greek letter when there is no
    //    additional base), so that  DELTAW_dot → \Delta \dot{W}  rather than
    //    \dot{\Delta W}.
    std::string core = baseSymbol.empty() ? greekLatex : baseSymbol;

    for (const auto& mod : modifiers) {
        if (mod == "bar") core = "\\bar{" + core + "}";
        else if (mod == "dot") core = "\\dot{" + core + "}";
    }

    // Prepend Greek prefix when there is also a separate base symbol
    if (!greekLatex.empty() && !baseSymbol.empty()) {
        core = greekLatex + " " + core;
    }

    // 6. Append subscript
    if (!subscripts.empty()) {
        std::string sub;
        for (size_t i = 0; i < subscripts.size(); ++i) {
            if (i > 0) sub += ",";
            sub += subscripts[i];
        }
        core += "_{" + sub + "}";
    }

    return core;
}

std::string IR::expressionToLatex(const ExprPtr& expr) const {
    if (!expr) return "";
    
    return std::visit([this](const auto& node) -> std::string {
        using T = std::decay_t<decltype(node)>;
        
        if constexpr (std::is_same_v<T, NumberLiteral>) {
            double val = node.value;
            // Format scientific notation nicely
            if (std::abs(val) >= 1e4 || (std::abs(val) < 1e-3 && val != 0)) {
                std::ostringstream oss;
                oss << std::scientific << val;
                std::string s = oss.str();
                // Convert 1e5 to 1 \times 10^{5}
                size_t ePos = s.find('e');
                if (ePos != std::string::npos) {
                    std::string mantissa = s.substr(0, ePos);
                    std::string exp = s.substr(ePos + 1);
                    // Remove leading + or leading zeros in exponent
                    if (exp[0] == '+') exp = exp.substr(1);
                    while (exp.size() > 1 && exp[0] == '0') exp = exp.substr(1);
                    return mantissa + " \\times 10^{" + exp + "}";
                }
            }
            return std::to_string(val);
        } else if constexpr (std::is_same_v<T, StringLiteral>) {
            return "\\mathrm{" + node.value + "}";
        } else if constexpr (std::is_same_v<T, Variable>) {
            return variableToLatex(node.flattenedName());
        } else if constexpr (std::is_same_v<T, UnaryOp>) {
            return node.op + expressionToLatex(node.operand);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            std::string left = expressionToLatex(node.left);
            std::string right = expressionToLatex(node.right);
            
            if (node.op == "/") {
                return "\\frac{" + left + "}{" + right + "}";
            } else if (node.op == "*") {
                return left + " \\cdot " + right;
            } else if (node.op == "^") {
                return left + "^{" + right + "}";
            } else {
                return left + " " + node.op + " " + right;
            }
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            return functionCallToLatex(node);
        }
        return "";
    }, expr->node);
}

std::string functionCallToLatex(const FunctionCall& func) {
    std::string result;
    
    std::string lowerName = func.name;
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
    
    // Use the same THERMO_FUNCTIONS set for consistency
    if (THERMO_FUNCTIONS.find(lowerName) != THERMO_FUNCTIONS.end()) {
        result = "\\mathrm{" + func.name + "}";
    } else if (lowerName == "exp") {
        result = "\\exp";
    } else if (lowerName == "ln" || lowerName == "log") {
        result = "\\ln";
    } else if (lowerName == "sin" || lowerName == "cos" || lowerName == "tan") {
        result = "\\" + lowerName;
    } else if (lowerName == "min" || lowerName == "max") {
        result = "\\mathrm{" + lowerName + "}";
    } else {
        result = "\\mathrm{" + func.name + "}";
    }
    
    result += " \\left(";
    
    bool first = true;
    for (const auto& arg : func.args) {
        if (!first) result += ", ";
        result += expressionToLatex(arg);
        first = false;
    }
    for (const auto& [name, arg] : func.namedArgs) {
        if (!first) result += ", ";
        result += "\\mbox{\\ " + name + "}=" + expressionToLatex(arg);
        first = false;
    }
    
    result += " \\right)";
    return result;
}

std::string expressionToLatex(const ExprPtr& expr) {
    if (!expr) return "";
    
    return std::visit([](const auto& node) -> std::string {
        using T = std::decay_t<decltype(node)>;
        
        if constexpr (std::is_same_v<T, NumberLiteral>) {
            return std::to_string(node.value);
        } else if constexpr (std::is_same_v<T, StringLiteral>) {
            return "\\mathrm{" + node.value + "}";
        } else if constexpr (std::is_same_v<T, Variable>) {
            return variableToLatex(node.flattenedName());
        } else if constexpr (std::is_same_v<T, UnaryOp>) {
            return node.op + expressionToLatex(node.operand);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            std::string left = expressionToLatex(node.left);
            std::string right = expressionToLatex(node.right);
            
            if (node.op == "/") {
                return "\\frac{" + left + "}{" + right + "}";
            } else if (node.op == "*") {
                return left + " \\cdot " + right;
            } else if (node.op == "^") {
                return left + "^{" + right + "}";
            } else {
                return left + " " + node.op + " " + right;
            }
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            return functionCallToLatex(node);
        }
        return "";
    }, expr->node);
}

std::string IR::toLatex() const {
    std::ostringstream oss;
    
    oss << "\\documentclass{article}\n";
    oss << "\\usepackage{amsmath}\n";
    oss << "\\begin{document}\n\n";
    oss << "\\section*{Equations}\n\n";
    
    for (const auto& eq : equations_) {
        oss << "\\begin{equation}\n";
        oss << expressionToLatex(eq.lhs) << " = " << expressionToLatex(eq.rhs) << "\n";
        oss << "\\end{equation}\n\n";
    }
    
    oss << "\\end{document}\n";
    return oss.str();
}

// ============================================================================
// Export Functions
// ============================================================================

std::string IR::variableTableToMarkdown() const {
    std::ostringstream oss;
    
    oss << "| # | Variable | LaTeX | Solution | Units |\n";
    oss << "|---|----------|-------|----------|-------|\n";
    
    int num = 1;
    for (const auto& [name, info] : variables_) {
        oss << "| " << num++;
        oss << " | `" << name << "`";
        oss << " | $" << info.latexLabel << "$";
        if (info.solutionValue) {
            oss << " | " << *info.solutionValue;
        } else if (info.solutionStringValue) {
            oss << " | '" << *info.solutionStringValue << "'";
        } else {
            oss << " | -";
        }
        oss << " | " << info.units;
        oss << " |\n";
    }
    
    return oss.str();
}

std::string IR::variableTableToCSV() const {
    std::ostringstream oss;
    
    oss << "Index,Variable,LaTeX,Solution,Type,Units,Description\n";
    
    for (const auto& [name, info] : variables_) {
        oss << info.index << ",";
        oss << "\"" << name << "\",";
        oss << "\"" << info.latexLabel << "\",";
        if (info.solutionValue) {
            oss << *info.solutionValue;
        } else if (info.solutionStringValue) {
            oss << "'" << *info.solutionStringValue << "'";
        }
        oss << ",";
        oss << "Real,";  // TODO: proper type detection
        oss << "\"" << info.units << "\",";
        oss << "\"" << info.description << "\"\n";
    }
    
    return oss.str();
}

std::string IR::toJSON() const {
    nlohmann::json j;
    
    // Variables
    nlohmann::json varsJson = nlohmann::json::array();
    for (const auto& [name, info] : variables_) {
        nlohmann::json varJson;
        varJson["name"] = name;
        varJson["latex"] = info.latexLabel;
        varJson["index"] = info.index;
        varJson["type"] = "Real";
        varJson["units"] = info.units;
        if (info.solutionValue) {
            varJson["solution"] = *info.solutionValue;
        } else if (info.solutionStringValue) {
            varJson["solution"] = *info.solutionStringValue;
        }
        varsJson.push_back(varJson);
    }
    j["variables"] = varsJson;
    
    // Equations
    nlohmann::json eqsJson = nlohmann::json::array();
    for (const auto& eq : equations_) {
        nlohmann::json eqJson;
        eqJson["id"] = eq.id;
        eqJson["text"] = eq.originalText;
        eqJson["line"] = eq.sourceLine;
        eqJson["variables"] = std::vector<std::string>(eq.variables.begin(), eq.variables.end());
        eqJson["outputVariable"] = eq.outputVariable;
        eqJson["blockId"] = eq.blockId;
        eqsJson.push_back(eqJson);
    }
    j["equations"] = eqsJson;
    
    // Statistics
    j["stats"]["variableCount"] = variables_.size();
    j["stats"]["equationCount"] = equations_.size();
    j["stats"]["isSquare"] = isSquare();
    
    return j.dump(2);
}

}  // namespace coolsolve
