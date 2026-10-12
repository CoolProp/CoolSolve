#pragma once

#include "ast.h"
#include "diagnostic.h"
#include <functional>
#include <string>
#include <optional>
#include <vector>

namespace coolsolve {

// ============================================================================
// Parser Error
// ============================================================================

struct ParseError {
    int line;
    int column;
    std::string message;
    std::string context;  // The source line where error occurred
};

// ============================================================================
// Parser Result
// ============================================================================

/**
 * One `$IF` / `$IFNOT` directive as resolved while parsing (CS-GAP-IF-DIRECTIVE;
 * see conditional_directives.h). The statements of the removed branch are not in
 * `ParseResult::program`. Reported in the debug folder (`report.md`).
 */
struct ConditionalDirectiveInfo {
    int line = 0;             ///< Line of the `$IF` / `$IFNOT`.
    std::string text;         ///< The directive as written, comments removed (`$IF fluid$='Water'`).
    bool evaluated = false;   ///< False: inside a removed branch (not evaluated) or not evaluable (an error).
    bool thenKept = false;    ///< Only if `evaluated`: the lines up to `$ELSE`/`$ENDIF` are kept (else: the lines after `$ELSE`).
    int elseLine = 0;         ///< Line of the `$ELSE`, 0 when there is none.
    int endLine = 0;          ///< Line of the `$ENDIF`.
};

struct ParseResult {
    bool success;
    Program program;
    std::vector<ParseError> errors;
    DiagnosticCollector diagnostics;
    std::vector<ConditionalDirectiveInfo> conditionalDirectives;  ///< Every `$IF`/`$IFNOT`, in file order.
    
    // Statistics
    int totalLines = 0;
    int equationCount = 0;
    int commentCount = 0;
    int directiveCount = 0;
};

// ============================================================================
// EES Parser
// ============================================================================

class EESParser {
public:
    EESParser();
    ~EESParser();
    
    // Parse EES source code string
    ParseResult parse(const std::string& source, const std::string& filename = "<input>");
    
    // Parse from file
    ParseResult parseFile(const std::string& filepath);
    
    // Get the last error message for debugging
    std::string getLastError() const;
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

// ============================================================================
// Utility functions
// ============================================================================

// Read a file into a string
std::optional<std::string> readFile(const std::string& filepath);

// Collect all variables from an expression
void collectVariables(const ExprPtr& expr, std::vector<std::string>& vars);

/**
 * Value of an expression known before the solve: numbers combined with `+ - * / ^`, a unary sign and
 * parentheses, a name (a variable without subscript) standing for `lookup(name)`. This is the resolver of
 * the integer bounds of the DUPLICATE loops and of the array ranges `X[a..b]` (the parser resolves the
 * names, the IR evaluates the copies of a nested loop) and of the `$CONSTANT` values.
 *
 * @param lookup  Value of a name; nullopt (or no lookup) when the name has none.
 * @param unknown When given and still empty, receives the first name (or array element) without a value.
 * @return nullopt for a name without value, a function call, a string, a comparison or a result that is
 *         not finite.
 */
std::optional<double> evaluateConstantExpression(
    const ExprPtr& expr, const std::function<std::optional<double>(const std::string&)>& lookup = {},
    std::string* unknown = nullptr);

/// `value` as an integer bound (DUPLICATE, array range) when it is an integer up to rounding errors
/// (`0.1*170` is 17); nullopt otherwise.
std::optional<int> asIntegerBound(double value);

// Convert AST to string representation (for debugging)
std::string astToString(const ExprPtr& expr);
std::string astToString(const StmtPtr& stmt);
std::string astToString(const Program& program);

}  // namespace coolsolve
