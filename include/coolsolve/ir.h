#pragma once

#include "ast.h"
#include <map>
#include <set>
#include <string>
#include <vector>
#include <optional>
#include <algorithm>

namespace coolsolve {

// Case-insensitive comparator for variable names (EES style)
struct CaseInsensitiveLess {
    bool operator()(const std::string& a, const std::string& b) const {
        std::string la = a;
        std::string lb = b;
        std::transform(la.begin(), la.end(), la.begin(), ::tolower);
        std::transform(lb.begin(), lb.end(), lb.begin(), ::tolower);
        return la < lb;
    }
};

// ============================================================================
// Variable Information
// ============================================================================

enum class VariableType {
    Real,
    String,
    Integer
};

struct VariableInfo {
    std::string originalName;      // Exact name from EES source (e.g., "T_su[1]")
    std::string internalId;        // Unique identifier (e.g., "VAR_001" or "T_su_1")
    std::string latexLabel;        // LaTeX formatted (e.g., "T_{su,1}")
    VariableType type = VariableType::Real;
    
    // Value and bounds
    std::optional<double> guessValue;
    std::optional<double> lowerBound;
    std::optional<double> upperBound;
    std::optional<double> solutionValue;  // Solution value from .solution file (for Real variables)
    std::optional<std::string> solutionStringValue; // Solution value from .solution file (for String variables)
    
    // Metadata
    std::string units;
    std::string description;
    
    // Automatic inference
    std::string inferredFluid;     // E.g., "Water"
    std::string inferredProperty;  // E.g., "T", "P", "H"
    
    // Index in the system (for matrix operations)
    int index = -1;
};

// ============================================================================
// Equation Information
// ============================================================================

struct EquationInfo {
    int id;                        // Unique equation ID
    std::string originalText;      // Original equation text from source
    std::string latexText;         // LaTeX representation
    int sourceLine;                // Line number in source
    
    // Expression trees
    ExprPtr lhs;
    ExprPtr rhs;
    
    // Procedure call (if this is a CALL statement)
    std::optional<ProcedureCall> procedureCall;
    
    // Variables appearing in this equation
    std::set<std::string, CaseInsensitiveLess> variables;
    
    // After matching: the variable this equation determines
    std::string outputVariable;
    
    // Block assignment after SCC analysis
    int blockId = -1;
};

// ============================================================================
// Intermediate Representation
// ============================================================================

class IR {
public:
    // Build IR from parsed AST
    static IR fromAST(const Program& program);
    
    // Access variables
    const std::map<std::string, VariableInfo, CaseInsensitiveLess>& getVariables() const { return variables_; }
    const VariableInfo* getVariable(const std::string& name) const;
    int getVariableCount() const { return static_cast<int>(variables_.size()); }
    
    // Access equations
    const std::vector<EquationInfo>& getEquations() const { return equations_; }
    int getEquationCount() const { return static_cast<int>(equations_.size()); }
    
    // Access functions and procedures
    const std::vector<FunctionDefinition>& getFunctions() const { return functions_; }
    const std::vector<ProcedureDefinition>& getProcedures() const { return procedures_; }
    
    // Get the incidence matrix (which variables appear in which equations)
    // incidenceMatrix[eqId] = set of variable names in that equation
    const std::vector<std::set<std::string, CaseInsensitiveLess>>& getIncidenceMatrix() const { return incidenceMatrix_; }
    
    // Generate LaTeX representation
    std::string toLatex() const;
    
    // Export variable table to various formats
    std::string variableTableToMarkdown() const;
    std::string variableTableToCSV() const;
    std::string toJSON() const;
    
    // Diagnostics
    int getAlgebraicVariableCount() const;
    int getNonConstantVariableCount() const;
    bool isSquare() const { return getNonConstantVariableCount() == getEquationCount(); }
    std::vector<std::string> getUnmatchedVariables() const;
    std::vector<int> getUnmatchedEquations() const;
    
    // Load solution values from .solution or .initials file
    // Format: variable=value # comment
    // Returns number of variables matched
    int loadSolutionFromFile(const std::string& filePath);
    
    // Load initial values from .initials file
    // This is an alias for loadSolutionFromFile but specifically for initials
    int loadInitialsFromFile(const std::string& filePath);
    
private:
    std::map<std::string, VariableInfo, CaseInsensitiveLess> variables_;
    std::vector<EquationInfo> equations_;
    std::vector<FunctionDefinition> functions_;
    std::vector<ProcedureDefinition> procedures_;
    std::vector<std::set<std::string, CaseInsensitiveLess>> incidenceMatrix_;
    
    // Helper functions
    void extractVariables(const ExprPtr& expr, std::set<std::string, CaseInsensitiveLess>& vars);
    std::string expressionToLatex(const ExprPtr& expr) const;
    void addVariable(const std::string& name);
};

// ============================================================================
// LaTeX Generation
// ============================================================================

// Convert a variable name to LaTeX format
// Examples: "T_oil_su" -> "T_{oil,su}", "M_dot" -> "\dot{M}"
std::string variableToLatex(const std::string& name);

// Convert an expression to LaTeX
std::string expressionToLatex(const ExprPtr& expr);

// Convert a function call to LaTeX (handles thermodynamic functions specially)
std::string functionCallToLatex(const FunctionCall& func);

}  // namespace coolsolve
