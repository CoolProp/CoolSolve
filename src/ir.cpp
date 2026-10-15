#include "coolsolve/ir.h"
#include "coolsolve/parser.h"
#include "coolsolve/constants.h"
#include <nlohmann/json.hpp>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cctype>
#include <unordered_set>

namespace coolsolve {

// ============================================================================
// Expression/Statement deep-copy with iterator substitution
// ============================================================================

// An index that is an arithmetic expression of numbers (`T[i+1]` once i is replaced by 1) is the integer it
// stands for, so that the copy names the element T[2], not `T[(1 + 1)]`; any other index is returned as is.
static ExprPtr foldIndex(ExprPtr index) {
    if (index && (index->is<BinaryOp>() || index->is<UnaryOp>()))
        if (const auto value = evaluateConstantExpression(index))
            if (const auto integer = asIntegerBound(*value)) return makeNumber(*integer, index->sourceLineNumber);
    return index;
}

// Deep-copy an expression, replacing occurrences of iteratorVar with the given value.
// - In Variable names: if name matches iteratorVar (case-insensitive), replace with NumberLiteral
// - In Variable indices: recursively substitute
// - In all sub-expressions: recursively substitute
static ExprPtr cloneExprWithSubstitution(const ExprPtr& expr, const std::string& iteratorVar, double value) {
    if (!expr) return nullptr;
    auto result = std::make_shared<Expression>();
    result->sourceLineNumber = expr->sourceLineNumber;

    if (expr->is<NumberLiteral>()) {
        result->node = expr->as<NumberLiteral>();
    } else if (expr->is<StringLiteral>()) {
        result->node = expr->as<StringLiteral>();
    } else if (expr->is<Variable>()) {
        const auto& var = expr->as<Variable>();
        // Check if this variable IS the iterator variable (bare reference, no indices)
        if (eesStringEqual(var.name, iteratorVar) && var.indices.empty()) {
            result->node = NumberLiteral{value};
        } else {
            Variable newVar;
            newVar.name = var.name;
            for (const auto& idx : var.indices) {
                newVar.indices.push_back(foldIndex(cloneExprWithSubstitution(idx, iteratorVar, value)));
            }
            result->node = std::move(newVar);
        }
    } else if (expr->is<UnaryOp>()) {
        const auto& op = expr->as<UnaryOp>();
        UnaryOp newOp;
        newOp.op = op.op;
        newOp.operand = cloneExprWithSubstitution(op.operand, iteratorVar, value);
        result->node = std::move(newOp);
    } else if (expr->is<BinaryOp>()) {
        const auto& op = expr->as<BinaryOp>();
        BinaryOp newOp;
        newOp.op = op.op;
        newOp.left = cloneExprWithSubstitution(op.left, iteratorVar, value);
        newOp.right = cloneExprWithSubstitution(op.right, iteratorVar, value);
        result->node = std::move(newOp);
    } else if (expr->is<FunctionCall>()) {
        const auto& fc = expr->as<FunctionCall>();
        FunctionCall newFc;
        newFc.name = fc.name;
        // the index of an indexed aggregation hides an enclosing index of the same name in its expression
        const bool indexed = isIndexedAggregate(fc);
        const bool hidden = indexed && eesStringEqual(fc.namedArgs[0].first, iteratorVar);
        for (size_t k = 0; k < fc.args.size(); ++k) {
            newFc.args.push_back(hidden && k == 0 ? cloneExprWithSubstitution(fc.args[0], std::string(), 0.0)
                                                  : cloneExprWithSubstitution(fc.args[k], iteratorVar, value));
        }
        for (const auto& [key, val] : fc.namedArgs) {
            newFc.namedArgs.push_back({key, cloneExprWithSubstitution(val, iteratorVar, value)});
        }
        result->node = std::move(newFc);
        if (indexed)   // its bounds are numbers once the enclosing index is replaced
            if (auto expanded = expandIndexedAggregate(result)) return expanded;
    }
    return result;
}

bool isIndexedAggregate(const FunctionCall& call) {
    if (call.namedArgs.size() != 1 || call.args.size() != 2) return false;
    std::string lower = call.name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return lower == "sum" || lower == "product";
}

namespace {
constexpr int kMaxAggregateTerms = 100000;

/// The integer value of the `which` bound of `owner` ("DUPLICATE i", "SUM index i") in a copy of the
/// enclosing loop, where the parser left an expression of its index; a ModelError when it is not an integer.
int copyBound(const ExprPtr& expr, const std::string& owner, const char* which, int line) {
    const auto value = evaluateConstantExpression(expr);
    const auto integer = value ? asIntegerBound(*value) : std::nullopt;
    if (!integer) {
        std::ostringstream what;
        what << owner << ": the " << which << " bound is ";
        if (value) what << std::setprecision(15) << *value; else what << astToString(expr);
        what << " in a copy of the enclosing loop, not an integer";
        throw ModelError(line, what.str());
    }
    return *integer;
}
}  // namespace

ExprPtr expandIndexedAggregate(const ExprPtr& callExpr) {
    const auto& call = callExpr->as<FunctionCall>();
    const std::string& index = call.namedArgs[0].first;
    const ExprPtr& firstBound = call.namedArgs[0].second;
    const ExprPtr& lastBound = call.args[1];
    if (!evaluateConstantExpression(firstBound) || !evaluateConstantExpression(lastBound)) return nullptr;

    std::string owner = call.name;
    std::transform(owner.begin(), owner.end(), owner.begin(), ::toupper);
    owner += " index " + index;
    const int line = callExpr->sourceLineNumber;
    const int first = copyBound(firstBound, owner, "first", line);
    const int last = copyBound(lastBound, owner, "last", line);
    if (first > last)
        throw ModelError(line, owner + ": the first bound (" + std::to_string(first) + ") is greater than the last (" +
                               std::to_string(last) + ")");
    if (static_cast<long long>(last) - first >= kMaxAggregateTerms)
        throw ModelError(line, owner + ": more than " + std::to_string(kMaxAggregateTerms) + " terms");

    FunctionCall expanded;
    expanded.name = call.name;
    for (int k = first; k <= last; ++k)
        expanded.args.push_back(cloneExprWithSubstitution(call.args[0], index, static_cast<double>(k)));
    auto result = std::make_shared<Expression>();
    result->node = std::move(expanded);
    result->sourceLineNumber = callExpr->sourceLineNumber;
    return result;
}

// Deep-copy a statement, applying iterator substitution to all contained expressions
static StmtPtr cloneStmtWithSubstitution(const StmtPtr& stmt, const std::string& iteratorVar, double value) {
    if (!stmt) return nullptr;
    auto result = std::make_shared<Statement>();
    result->sourceLineNumber = stmt->sourceLineNumber;

    if (stmt->is<Equation>()) {
        const auto& eq = stmt->as<Equation>();
        Equation newEq;
        newEq.lhs = cloneExprWithSubstitution(eq.lhs, iteratorVar, value);
        newEq.rhs = cloneExprWithSubstitution(eq.rhs, iteratorVar, value);
        newEq.units = eq.units;
        newEq.comment = eq.comment;
        result->node = std::move(newEq);
    } else if (stmt->is<Assignment>()) {
        const auto& asgn = stmt->as<Assignment>();
        Assignment newAsgn;
        newAsgn.lhs = cloneExprWithSubstitution(asgn.lhs, iteratorVar, value);
        newAsgn.rhs = cloneExprWithSubstitution(asgn.rhs, iteratorVar, value);
        result->node = std::move(newAsgn);
    } else if (stmt->is<ProcedureCall>()) {
        const auto& call = stmt->as<ProcedureCall>();
        ProcedureCall newCall;
        newCall.name = call.name;
        for (const auto& arg : call.inputArgs) {
            newCall.inputArgs.push_back(cloneExprWithSubstitution(arg, iteratorVar, value));
        }
        for (const auto& outVar : call.outputVars) {
            Variable newVar;
            newVar.name = outVar.name;
            for (const auto& idx : outVar.indices) {
                newVar.indices.push_back(foldIndex(cloneExprWithSubstitution(idx, iteratorVar, value)));
            }
            newCall.outputVars.push_back(std::move(newVar));
        }
        for (const auto& outExpr : call.outputExprs) {
            newCall.outputExprs.push_back(cloneExprWithSubstitution(outExpr, iteratorVar, value));
        }
        result->node = std::move(newCall);
    } else if (stmt->is<Comment>()) {
        result->node = stmt->as<Comment>();
    } else if (stmt->is<Directive>()) {
        result->node = stmt->as<Directive>();
    } else if (stmt->is<Duplicate>()) {
        // Nested DUPLICATE: clone with substitution in start/end expressions and body
        const auto& dup = stmt->as<Duplicate>();
        Duplicate newDup;
        newDup.iteratorVar = dup.iteratorVar;
        newDup.start = cloneExprWithSubstitution(dup.start, iteratorVar, value);
        newDup.end = cloneExprWithSubstitution(dup.end, iteratorVar, value);
        for (const auto& bodyStmt : dup.body) {
            newDup.body.push_back(cloneStmtWithSubstitution(bodyStmt, iteratorVar, value));
        }
        result->node = std::move(newDup);
    }
    return result;
}

// Recursively expand a single Duplicate statement into a flat list of expanded statements.
// Handles nested DUPLICATE by recursively expanding inner DUPLICATE loops. The parser has resolved the
// bounds (CS-BUG-DUPLICATE-VAR-BOUND): an integer, or for a nested loop an arithmetic expression of the
// indices of the enclosing loops, which are numbers once this copy is made. A bound that is still not an
// integer is an error (the loop was silently dropped).
static void expandDuplicate(const Duplicate& dup, int line, std::vector<StmtPtr>& out) {
    const std::string owner = "DUPLICATE " + dup.iteratorVar;
    const int startVal = copyBound(dup.start, owner, "first", line);
    const int endVal = copyBound(dup.end, owner, "last", line);

    int step = (startVal <= endVal) ? 1 : -1;
    for (int i = startVal; step > 0 ? i <= endVal : i >= endVal; i += step) {
        for (const auto& bodyStmt : dup.body) {
            auto cloned = cloneStmtWithSubstitution(bodyStmt, dup.iteratorVar, static_cast<double>(i));
            if (cloned && cloned->is<Duplicate>()) {
                // Recursively expand nested DUPLICATE
                expandDuplicate(cloned->as<Duplicate>(), cloned->sourceLineNumber, out);
            } else if (cloned) {
                out.push_back(std::move(cloned));
            }
        }
    }
}

// ============================================================================
// IR Building from AST
// ============================================================================

IR IR::fromAST(const Program& program) {
    IR ir;
    int eqId = 0;
    std::vector<std::string> pendingComments;  // Accumulate "" display comments

    // Expand all top-level DUPLICATE statements into flat statement list
    std::vector<StmtPtr> expandedStatements;
    for (const auto& stmt : program.statements) {
        if (stmt->is<Duplicate>()) {
            expandDuplicate(stmt->as<Duplicate>(), stmt->sourceLineNumber, expandedStatements);
        } else {
            expandedStatements.push_back(stmt);
        }
    }
    
    ir.lowerCallExpressionOutputs(expandedStatements);

    for (const auto& stmt : expandedStatements) {
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

namespace {

bool containsVariable(const ExprPtr& expr) {
    if (!expr) return false;
    if (expr->is<Variable>()) return true;
    if (expr->is<UnaryOp>()) return containsVariable(expr->as<UnaryOp>().operand);
    if (expr->is<BinaryOp>())
        return containsVariable(expr->as<BinaryOp>().left) || containsVariable(expr->as<BinaryOp>().right);
    if (expr->is<FunctionCall>()) {
        for (const auto& arg : expr->as<FunctionCall>().args)
            if (containsVariable(arg)) return true;
        for (const auto& named : expr->as<FunctionCall>().namedArgs)
            if (containsVariable(named.second)) return true;
    }
    return false;
}

ExprPtr makeBinary(const std::string& op, ExprPtr left, ExprPtr right) {
    auto expr = std::make_shared<Expression>();
    expr->node = BinaryOp{op, std::move(left), std::move(right)};
    return expr;
}

/// Solve `expr = rhs` for the one variable of `expr` when `expr` is built from that variable with
/// `+`, `-`, `*` and a unary minus, or `/` by something without variables, the other operands
/// having none: `c/1000 = a` gives `c = a*1000`. Returns false (and leaves the outputs alone)
/// for any other shape, which stays an implicit equation.
bool solveForOnlyVariable(const ExprPtr& expr, const ExprPtr& rhs, ExprPtr& variable, ExprPtr& value) {
    if (expr->is<Variable>()) {
        for (const auto& idx : expr->as<Variable>().indices)
            if (containsVariable(idx)) return false;
        variable = expr;
        value = rhs;
        return true;
    }
    if (expr->is<UnaryOp>()) {
        const auto& op = expr->as<UnaryOp>();
        if (op.op == "+") return solveForOnlyVariable(op.operand, rhs, variable, value);
        if (op.op != "-") return false;
        auto negated = std::make_shared<Expression>();
        negated->node = UnaryOp{"-", rhs};
        return solveForOnlyVariable(op.operand, negated, variable, value);
    }
    if (!expr->is<BinaryOp>()) return false;
    const auto& op = expr->as<BinaryOp>();
    const bool inLeft = containsVariable(op.left);
    if (inLeft == containsVariable(op.right)) return false;   // both sides: not a single variable
    const ExprPtr& with = inLeft ? op.left : op.right;
    const ExprPtr& other = inLeft ? op.right : op.left;
    if (op.op == "+") return solveForOnlyVariable(with, makeBinary("-", rhs, other), variable, value);
    if (op.op == "*") return solveForOnlyVariable(with, makeBinary("/", rhs, other), variable, value);
    if (op.op == "-")
        return solveForOnlyVariable(with, inLeft ? makeBinary("+", rhs, other) : makeBinary("-", other, rhs),
                                    variable, value);
    if (op.op == "/" && inLeft) return solveForOnlyVariable(with, makeBinary("*", rhs, other), variable, value);
    return false;
}

}  // namespace

// EES lets an output of a CALL be an expression (`CALL BRINEPROP('Specheat', b$, X, T : c_p/1000)`):
// the output of the routine equals the expression, which the system solves for like any equation
// (CS-GAP-CALL-EXPR-OUT). The output goes to an auxiliary variable `__out<k>_<routine>_L<line>`
// (k = position of the output; `_2`, `_3`... for the copies made by a DUPLICATE loop and when the
// name is already used by the model) and the equation `expression = auxiliary` follows the CALL
// (solved for the variable when the expression is a simple chain of `+ - * /` on one variable).
// A variable output is untouched, and so is a program without expression outputs.
void IR::lowerCallExpressionOutputs(std::vector<StmtPtr>& statements) {
    const auto hasExpressionOutput = [](const StmtPtr& stmt) {
        return stmt->is<ProcedureCall>() && !stmt->as<ProcedureCall>().outputExprs.empty();
    };
    if (std::none_of(statements.begin(), statements.end(), hasExpressionOutput)) return;

    // Every name of the model, so that an auxiliary variable never merges with one of them
    std::set<std::string, CaseInsensitiveLess> used;
    for (const auto& stmt : statements) {
        if (stmt->is<Equation>()) {
            extractVariables(stmt->as<Equation>().lhs, used);
            extractVariables(stmt->as<Equation>().rhs, used);
        } else if (stmt->is<ProcedureCall>()) {
            const auto& call = stmt->as<ProcedureCall>();
            for (const auto& arg : call.inputArgs) extractVariables(arg, used);
            for (const auto& out : call.outputVars) used.insert(out.flattenedName());
            for (const auto& out : call.outputExprs) extractVariables(out, used);
        }
    }

    std::vector<StmtPtr> lowered;
    lowered.reserve(statements.size());
    for (const auto& stmt : statements) {
        if (!hasExpressionOutput(stmt)) {
            lowered.push_back(stmt);
            continue;
        }
        ProcedureCall call = stmt->as<ProcedureCall>();
        std::vector<StmtPtr> ties;
        for (size_t i = 0; i < call.outputExprs.size(); ++i) {
            if (!call.outputExprs[i]) continue;
            const std::string base = "__out" + std::to_string(i + 1) + "_" + call.name + "_L"
                                     + std::to_string(stmt->sourceLineNumber);
            std::string name = base;
            for (int copy = 2; used.count(name); ++copy) name = base + "_" + std::to_string(copy);
            used.insert(name);
            call.outputVars[i] = Variable{name, {}};
            // `c/1000 = aux` is written `c = aux*1000` when it can be: the form of an explicit equation
            // (solved without iterating, and scaled like the equations a user writes by hand)
            ExprPtr auxiliary = makeVariable(name, stmt->sourceLineNumber);
            ExprPtr variable, value;
            if (solveForOnlyVariable(call.outputExprs[i], auxiliary, variable, value))
                ties.push_back(makeEquation(variable, value, stmt->sourceLineNumber));
            else
                ties.push_back(makeEquation(call.outputExprs[i], auxiliary, stmt->sourceLineNumber));
        }
        call.outputExprs.clear();
        auto lowerCall = std::make_shared<Statement>();
        lowerCall->node = std::move(call);
        lowerCall->sourceLineNumber = stmt->sourceLineNumber;
        lowered.push_back(std::move(lowerCall));
        for (auto& tie : ties) lowered.push_back(std::move(tie));
    }
    statements = std::move(lowered);
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
    
    std::visit([&vars](const auto& node) {
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
            if (isIndexedAggregate(node)) {
                // not expanded yet (the check of a body): its index is local, only its bounds and the variables
                // of the expression, read with the index replaced by a number, belong to the equation
                extractVariables(cloneExprWithSubstitution(node.args[0], node.namedArgs[0].first, 1.0), vars);
                extractVariables(node.namedArgs[0].second, vars);
                extractVariables(node.args[1], vars);
                return;
            }
            // For thermodynamic functions, skip the first positional arg (fluid name)
            bool isThermo = isThermoFunction(node.name);

            // For lookup/interpolate functions, all StringLiteral args (table
            // name and column names) must not be treated as solver variables.
            static const std::unordered_set<std::string> LOOKUP_FUNCTIONS = {
                "lookup", "lookup$", "lookupcol", "lookupcol1", "lookupcellempty",
                "tablevalue", "tablevalue#", "tablerun#",
                "interpolate", "interpolate1", "interpolate2", "interpolate2dm",
                "nlookuprows", "nlookupcolumns",
                "sumlookup", "avglookup", "maxlookup", "minlookup", "stddevlookup",
            };
            std::string funcLower = node.name;
            std::transform(funcLower.begin(), funcLower.end(), funcLower.begin(), ::tolower);
            bool isLookup = LOOKUP_FUNCTIONS.count(funcLower) > 0;

            bool firstArg = true;

            // Extract from positional args
            for (const auto& arg : node.args) {
                if (isThermo && firstArg) {
                    // Skip first argument (fluid name) for thermo functions, except a call
                    // (LOOKUP$('table', i, 1)): its row and column are variables of the equation
                    firstArg = false;
                    if (!(arg && arg->template is<FunctionCall>())) continue;
                }
                // Skip string literals in lookup functions (they are table/column names)
                if (isLookup && arg && arg->template is<StringLiteral>()) {
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
        // Remove comments: a '#' that starts the line or follows a blank (a '#' is also a character
        // of names: Nu#_f = 6.2)
        for (size_t k = 0; k < line.size(); ++k) {
            if (line[k] == '#' && (k == 0 || std::isspace(static_cast<unsigned char>(line[k - 1])))) {
                line.erase(k);
                break;
            }
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

std::string escapeLatexText(const std::string& s) {
    std::string out;
    out.reserve(s.size() * 2);
    for (char c : s) {
        switch (c) {
            case '&':  out += "\\&";  break;
            case '%':  out += "\\%";  break;
            case '$':  out += "\\$";  break;
            case '#':  out += "\\#";  break;
            case '_':  out += "\\_";  break;
            case '{':  out += "\\{";  break;
            case '}':  out += "\\}";  break;
            case '~':  out += "\\textasciitilde{}"; break;
            case '^':  out += "\\textasciicircum{}"; break;
            case '\\': out += "\\textbackslash{}"; break;
            default:   out += c;
        }
    }
    return out;
}

std::string escapeLatexMath(const std::string& s) {
    std::string out;
    out.reserve(s.size() * 2);
    for (size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];
        switch (c) {
            case '&':  out += "\\&";  break;
            case '%':  out += "\\%";  break;
            case '$':  out += "\\$";  break;
            case '#':  out += "\\#";  break;
            case '{':  out += "\\{";  break;
            case '}':  out += "\\}";  break;
            case '|':  out += "\\vert{}"; break;
            case '~':  out += "\\sim{}"; break;
            case '^':  out += "\\hat{}"; break;
            case '\\': out += "\\backslash{}"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x80) { out += c; break; }
                // UTF-8: the characters accepted in names (ees_string.h); `\text` needs amsmath
                bool letter = false;
                const size_t len = eesNonAsciiNameCharLength(s, i, &letter);
                if (len == 0) { out += c; break; }
                const std::string ch = s.substr(i, len);
                i += len - 1;
                if (letter) out += ch == "\xC2\xB5" ? "\\mu" : "\\text{" + ch + "}";
                else if (ch == "\xC2\xB0") out += "{}^{\\circ}";
                else out += "{}^{" + std::string(1, ch == "\xC2\xB2" ? '2' : ch == "\xC2\xB3" ? '3' : '1') + "}";
        }
    }
    return out;
}

std::string markdownTableCell(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '|') out += '\\';
        out += c;
    }
    return out;
}

/// A segment of a name in math mode: the `|` of EES starts a superscript (`G|o` is G^o).
static std::string latexSegment(const std::string& seg) {
    const size_t bar = seg.find('|');
    if (bar == std::string::npos) return escapeLatexMath(seg);
    return escapeLatexMath(seg.substr(0, bar)) + "^{" + escapeLatexMath(seg.substr(bar + 1)) + "}";
}

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
    if (cleanName.empty()) return "\\text{" + escapeLatexText(name) + "}";

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
    if (segs.empty()) return "\\text{" + escapeLatexText(name) + "}";

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
                baseSymbol = latexSegment(first.substr(gk.size()));
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
        baseSymbol = latexSegment(first);
    }

    // 4. Classify remaining segments as modifiers or subscripts
    std::vector<std::string> modifiers;   // "bar", "dot"
    std::vector<std::string> subscripts;
    for (size_t i = 1; i < segs.size(); ++i) {
        if (segs[i] == "bar" || segs[i] == "dot") {
            modifiers.push_back(segs[i]);
        } else {
            subscripts.push_back(latexSegment(segs[i]));
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
            return "\\mathrm{" + escapeLatexMath(node.value) + "}";
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
        result = "\\mathrm{" + escapeLatexMath(func.name) + "}";
    } else if (lowerName == "exp") {
        result = "\\exp";
    } else if (lowerName == "ln" || lowerName == "log") {
        result = "\\ln";
    } else if (lowerName == "sin" || lowerName == "cos" || lowerName == "tan") {
        result = "\\" + lowerName;
    } else if (lowerName == "min" || lowerName == "max") {
        result = "\\mathrm{" + lowerName + "}";
    } else {
        result = "\\mathrm{" + escapeLatexMath(func.name) + "}";
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
        result += "\\mbox{\\ " + escapeLatexText(name) + "}=" + expressionToLatex(arg);
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
            return "\\mathrm{" + escapeLatexMath(node.value) + "}";
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
        oss << " | `" << markdownTableCell(name) << "`";
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
