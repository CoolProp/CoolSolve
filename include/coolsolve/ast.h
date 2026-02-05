#pragma once

#include <memory>
#include <string>
#include <variant>
#include <vector>

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

struct UnaryOp {
    std::string op;  // "-", "+"
    ExprPtr operand;
};

struct BinaryOp {
    std::string op;  // "+", "-", "*", "/", "^"
    ExprPtr left;
    ExprPtr right;
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
    bool isBlock;  // true for {} or "" comments, false for // comments
};

struct Directive {
    std::string name;  // e.g., "$ifnot", "$endif", "$bookmark"
    std::string content;
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
};

struct FunctionDefinition {
    std::string name;
    std::vector<std::string> parameters;
    std::vector<StmtPtr> body;
};

struct ProcedureDefinition {
    std::string name;
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    std::vector<StmtPtr> body;
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
        ProcedureDefinition
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

inline StmtPtr makeComment(const std::string& text, bool isBlock = true, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = Comment{text, isBlock};
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

inline StmtPtr makeProcedureCall(const std::string& name, std::vector<ExprPtr> inputArgs, std::vector<Variable> outputVars, int line = 0) {
    auto stmt = std::make_shared<Statement>();
    stmt->node = ProcedureCall{name, std::move(inputArgs), std::move(outputVars)};
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

}  // namespace coolsolve
