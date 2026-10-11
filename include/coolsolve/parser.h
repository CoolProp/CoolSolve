#pragma once

#include "ast.h"
#include "diagnostic.h"
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

// Convert AST to string representation (for debugging)
std::string astToString(const ExprPtr& expr);
std::string astToString(const StmtPtr& stmt);
std::string astToString(const Program& program);

}  // namespace coolsolve
