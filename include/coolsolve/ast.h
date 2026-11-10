#pragma once

#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "coolsolve/integral/integral_table.h"  // for IntegralTableSpec (Directive payload)

namespace coolsolve {

// Forward declarations
struct Expression;
struct Statement;

// Type aliases for smart pointers
using ExprPtr = std::shared_ptr<Expression>;
using StmtPtr = std::shared_ptr<Statement>;

// ============================================================================
// Expression Types
// ============================================================================

struct NumberLiteral {
    double value;
};

struct StringLiteral {
    std::string value;
};

struct Variable {
    std::string name;
    std::vector<ExprPtr> indices;  // For array access like P[1], empty for scalars
    
    // Returns flattened name for array variables (e.g., "P_1" for P[1])
    std::string flattenedName() const;
};

/**
 * @brief The operator of a UnaryOp or a BinaryOp as a code, so that an evaluation switches on it instead of
 *        comparing strings. Set from the text of the operator where the node is built (the default member
 *        initializer of `code`); `Other` for a text without a code (the `..` of an array range).
 */
enum class OpCode : unsigned char { Other, Add, Sub, Mul, Div, Pow, Gt, Lt, Ge, Le, Eq, Ne, And, Or };

/// The code of the operator text \p op ("+", "<>", "and", ...); OpCode::Other for any other text.
inline OpCode opCodeOf(const std::string& op) {
    if (op.size() == 1) {
        switch (op[0]) {
            case '+': return OpCode::Add;
            case '-': return OpCode::Sub;
            case '*': return OpCode::Mul;
            case '/': return OpCode::Div;
            case '^': return OpCode::Pow;
            case '>': return OpCode::Gt;
            case '<': return OpCode::Lt;
            case '=': return OpCode::Eq;
            default: return OpCode::Other;
        }
    }
    if (op == ">=") return OpCode::Ge;
    if (op == "<=") return OpCode::Le;
    if (op == "<>") return OpCode::Ne;
    if (op == "and") return OpCode::And;
    if (op == "or") return OpCode::Or;
    return OpCode::Other;
}

struct UnaryOp {
    std::string op;  // "-", "+"
    ExprPtr operand;
    OpCode code = opCodeOf(op);   ///< the code of `op`, set when the node is built: changing `op` later needs a new code
};

struct BinaryOp {
    std::string op;  // "+", "-", "*", "/", "^"
    ExprPtr left;
    ExprPtr right;
    OpCode code = opCodeOf(op);   ///< the code of `op`, set when the node is built: changing `op` later needs a new code
};

struct FunctionCall {
    std::string name;
    std::vector<ExprPtr> args;
    // For thermodynamic functions, we store named arguments like T=T_1, P=P_1
    std::vector<std::pair<std::string, ExprPtr>> namedArgs;
};

// Expression is a variant of all possible expression types
struct Expression {
    std::variant<
        NumberLiteral,
        StringLiteral,
        Variable,
        UnaryOp,
        BinaryOp,
        FunctionCall
    > node;
    
    int sourceLineNumber = 0;  // For error reporting
    
    template<typename T>
    bool is() const { return std::holds_alternative<T>(node); }
    
    template<typename T>
    const T& as() const { return std::get<T>(node); }
    
    template<typename T>
    T& as() { return std::get<T>(node); }
};

// ============================================================================
// Statement Types
// ============================================================================

struct Equation {
    ExprPtr lhs;
    ExprPtr rhs;
    std::string units;  // Units annotation like "[kJ/kg]"
    std::string comment;  // Inline comment
};

struct Comment {
    std::string text;
    bool isBlock;   // true for {} or "" comments, false for // comments
    bool isDisplay; // true for "" display comments, false for {} hidden comments
};

struct Directive {
    std::string name;  // e.g., "$bookmark", "IntegralTable" ($if/$ifnot/$else/$endif are resolved by the parser and never stored)
    std::string content;
    // Optional structured payload, populated by the parser when `name` is
    // "IntegralTable". Other directives leave `hasIntegralTableSpec == false`.
    IntegralTableSpec integralTableSpec;
    bool hasIntegralTableSpec = false;
};

struct IfThenElse {
    ExprPtr condition;
    std::vector<StmtPtr> thenBranch;
    std::vector<StmtPtr> elseBranch;
};

struct Assignment {
    ExprPtr lhs;
    ExprPtr rhs;
};

struct ProcedureCall {
    std::string name;
    std::vector<ExprPtr> inputArgs;
    std::vector<Variable> outputVars;
    /// EES lets an output of a CALL of the main program be an expression (`CALL p(x : c/1000)`):
    /// the routine's output equals the expression. Empty when every output is a variable; otherwise
    /// one entry per output, null for a variable (`outputVars[i]`) and the expression where
    /// `outputVars[i]` is an empty placeholder. `IR::lowerCallExpressionOutputs` lowers them (CS-GAP-CALL-EXPR-OUT).
    std::vector<ExprPtr> outputExprs;
    /// The identifier of each input (`T`, `P`, `R`, ...) of a CALL whose inputs are named
    /// (`CALL psychprops(T=T_1, P=P_1, R=RH_1 : ...)`), parallel to `inputArgs`; empty for the positional form.
    std::vector<std::string> inputNames;
    /// The position, among the outputs of the procedure, of each entry of `outputVars` when some outputs are skipped
    /// (`CALL psychprops(... : , , h_1, , , w_1)`: the slots 2 and 5); empty when `outputVars[i]` is output `i`.
    std::vector<int> outputSlots;
};

struct FunctionDefinition {
    std::string name;
    std::vector<std::string> parameters;
    std::vector<StmtPtr> body;
    /// Empty for a routine of the model; the file it was imported from by `$INCLUDE` otherwise
    /// (`IncludedSource::describe()`: "my_props.eescode", "library model hx_ees_function (CSL-0331)").
    std::string includedFrom;
};

struct ProcedureDefinition {
    std::string name;
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    std::vector<StmtPtr> body;
    std::string includedFrom;   ///< As `FunctionDefinition::includedFrom`.
};

struct Duplicate {
    std::string iteratorVar;     // e.g., "i"
    ExprPtr start;               // start value expression
    ExprPtr end;                 // end value expression
    std::vector<StmtPtr> body;   // loop body statements
};

struct RepeatUntil {
    std::vector<StmtPtr> body;   // loop body statements
    ExprPtr condition;           // condition checked after each iteration
};

/// A numeric statement label (`100:`), the target of a GOTO. EES has labels and GOTO in the body of a
/// FUNCTION or PROCEDURE only; the parser refuses them elsewhere (CS-GAP-GOTO).
struct Label {
    int number;
};

/// `GOTO 100`: the body goes on after the label `100:` of the same body. The label is in the block of the
/// GOTO or in a block that contains it (the parser refuses a jump into a block).
struct Goto {
    int label;
};

struct Statement {
    std::variant<
        Equation,
        Assignment,
        Comment,
        Directive,
        IfThenElse,
        ProcedureCall,
        FunctionDefinition,
        ProcedureDefinition,
        Duplicate,
        RepeatUntil,
        Label,
        Goto
    > node;
    
    int sourceLineNumber = 0;
    
    template<typename T>
    bool is() const { return std::holds_alternative<T>(node); }
    
    template<typename T>
    const T& as() const { return std::get<T>(node); }
    
    template<typename T>
    T& as() { return std::get<T>(node); }
};

// ============================================================================
// Top-level Program
// ============================================================================

struct Program {
    std::vector<StmtPtr> statements;
    std::string sourceFilename;
};

// ============================================================================
// Helper functions for AST construction
// ============================================================================

inline ExprPtr makeNumber(double value, int line = 0) {
    auto expr = std::make_shared<Expression>();
    expr->node = NumberLiteral{value};
    expr->sourceLineNumber = line;
    return expr;
}

inline ExprPtr makeString(const std::string& value, int line = 0) {
    auto expr = std::make_shared<Expression>();
    expr->node = StringLiteral{value};
    expr->sourceLineNumber = line;
    return expr;
}

inline ExprPtr makeVariable(const std::string& name, int line = 0) {
    auto expr = std::make_shared<Expression>();
    expr->node = Variable{name, {}};
    expr->sourceLineNumber = line;
    return expr;
}

inline ExprPtr makeArrayVariable(const std::string& name, std::vector<ExprPtr> indices, int line = 0) {
    auto expr = std::make_shared<Expression>();
    expr->node = Variable{name, std::move(indices)};
    expr->sourceLineNumber = line;
    return expr;
}

inline ExprPtr makeUnaryOp(const std::string& op, ExprPtr operand, int line = 0) {
    auto expr = std::make_shared<Expression>();
    expr->node = UnaryOp{op, std::move(operand)};
    expr->sourceLineNumber = line;
    return expr;
}

inline ExprPtr makeBinaryOp(const std::string& op, ExprPtr left, ExprPtr right, int line = 0) {
    auto expr = std::make_shared<Expression>();
    expr->node = BinaryOp{op, std::move(left), std::move(right)};
    expr->sourceLineNumber = line;
    return expr;
}

inline ExprPtr makeFunctionCall(const std::string& name, std::vector<ExprPtr> args, int line = 0) {
    auto expr = std::make_shared<Expression>();
    expr->node = FunctionCall{name, std::move(args), {}};
    expr->sourceLineNumber = line;
    return expr;
}

inline StmtPtr makeEquation(ExprPtr lhs, ExprPtr rhs, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Equation{std::move(lhs), std::move(rhs), "", ""};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeComment(const std::string& text, bool isBlock = true, bool isDisplay = false, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Comment{text, isBlock, isDisplay};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeDirective(const std::string& name, const std::string& content, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Directive{name, content};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeAssignment(ExprPtr lhs, ExprPtr rhs, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Assignment{std::move(lhs), std::move(rhs)};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeProcedureCall(const std::string& name, std::vector<ExprPtr> inputArgs, std::vector<Variable> outputVars,
                                 int line = 0, std::vector<ExprPtr> outputExprs = {}) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = ProcedureCall{name, std::move(inputArgs), std::move(outputVars), std::move(outputExprs)};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeFunctionDefinition(const std::string& name, std::vector<std::string> parameters, std::vector<StmtPtr> body, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = FunctionDefinition{name, std::move(parameters), std::move(body)};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeProcedureDefinition(const std::string& name, std::vector<std::string> inputs, std::vector<std::string> outputs, std::vector<StmtPtr> body, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = ProcedureDefinition{name, std::move(inputs), std::move(outputs), std::move(body)};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeDuplicate(const std::string& iteratorVar, ExprPtr start, ExprPtr end, std::vector<StmtPtr> body, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Duplicate{iteratorVar, std::move(start), std::move(end), std::move(body)};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeIfThenElse(ExprPtr condition, std::vector<StmtPtr> thenBranch, std::vector<StmtPtr> elseBranch, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = IfThenElse{std::move(condition), std::move(thenBranch), std::move(elseBranch)};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeRepeatUntil(std::vector<StmtPtr> body, ExprPtr condition, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = RepeatUntil{std::move(body), std::move(condition)};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeLabel(int number, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Label{number};
    stmt->sourceLineNumber = line;
    return stmt;
}

inline StmtPtr makeGoto(int label, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Goto{label};
    stmt->sourceLineNumber = line;
    return stmt;
}

/// True for expressions that yield a string: 'literals', variables ending in '$' and
/// string-valued functions such as PHASE$() or LOOKUP$().
inline bool isStringExpression(const ExprPtr& expr) {
    if (expr->is<StringLiteral>()) return true;
    if (expr->is<Variable>()) {
        const std::string& n = expr->as<Variable>().name;
        return !n.empty() && n.back() == '$';
    }
    if (expr->is<FunctionCall>()) {
        const std::string& n = expr->as<FunctionCall>().name;
        return !n.empty() && n.back() == '$';
    }
    return false;
}

}  // namespace coolsolve
