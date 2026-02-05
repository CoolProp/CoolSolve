#pragma once

#include "ast.h"
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

struct ParseResult {
    bool success;
    Program program;
    std::vector<ParseError> errors;
    
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
