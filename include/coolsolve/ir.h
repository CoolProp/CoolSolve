#pragma once

#include "ast.h"
#include "ees_string.h"
#include <map>
#include <set>
#include <string>
#include <vector>
#include <optional>
#include <algorithm>
#include <stdexcept>

namespace coolsolve {

/// An error of the model found while the IR is built, at a line of the source: reported like a parse
/// error (CoolSolveRunner::run).
struct ModelError : std::runtime_error {
    int line;
    ModelError(int line, const std::string& message) : std::runtime_error(message), line(line) {}
};

/// True for the EES indexed aggregation `SUM(expr, i=a, b)` or `PRODUCT(expr, i=a, b)` as the parser stores
/// it: the arguments `{expr, b}` and the named argument `{"i", a}`. `i` is local to the aggregation.
bool isIndexedAggregate(const FunctionCall& call);

/// The variadic form of the indexed aggregation `call` (an Expression holding a FunctionCall for which
/// isIndexedAggregate is true): `SUM(expr[i:=a], ..., expr[i:=b])`, the substitution being that of the
/// copies of a DUPLICATE loop (the index is replaced by its value, `x[i+1]` becomes the element `x[2]`).
/// @return nullptr while a bound is not a number, i.e. it uses the index of an enclosing DUPLICATE loop: the
///         copies of that loop expand the aggregation once the index is replaced.
/// @throws ModelError when a bound is not an integer, when the first is greater than the last or when the
///         range holds more than 100000 terms.
ExprPtr expandIndexedAggregate(const ExprPtr& call);

// Case-insensitive comparator for variable names (EES style): ASCII and the accented letters
// of Latin-1 (UTF-8) are folded, see ees_string.h. T_a_réf and T_A_RÉF are one variable.
struct CaseInsensitiveLess {
    bool operator()(const std::string& a, const std::string& b) const {
        return eesCompareFolded(a, b) < 0;
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
    /// Identifies the CALL the equation is an output of (-1 for any other equation). The equations of one
    /// CALL share it; every CALL has its own, also the copies a DUPLICATE loop makes of it (same source line).
    int callId = -1;
    
    // Variables appearing in this equation (LHS and RHS combined)
    std::set<std::string, CaseInsensitiveLess> variables;
    /// Variables appearing only in the RHS (for explicit-solve detection: if output is not here, equation is x = expr(others))
    std::set<std::string, CaseInsensitiveLess> rhsVariables;
    
    // After matching: the variable this equation determines
    std::string outputVariable;
    
    // Block assignment after SCC analysis
    int blockId = -1;

    // Comments
    std::vector<std::string> precedingComments;  // Standalone "" comments before this equation
    std::string inlineComment;                    // Inline "" comment on this equation
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
    VariableInfo* getVariableMutable(const std::string& name);
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
    
    /**
     * Add to `vars` the names of the variables that `expr` reads, the way the equations of the
     * main program count them: an array element by its flattened name (`x[1]`, `x[i]`) plus the
     * variables of its indices; not the function names, the keys of named arguments (`T=`), the
     * fluid name of a property call (`enthalpy(R134a, ...)`) nor string literals. Also used by
     * the parser's check of FUNCTION/PROCEDURE bodies.
     */
    static void extractVariables(const ExprPtr& expr, std::set<std::string, CaseInsensitiveLess>& vars);
    
private:
    std::map<std::string, VariableInfo, CaseInsensitiveLess> variables_;
    std::vector<EquationInfo> equations_;
    std::vector<FunctionDefinition> functions_;
    std::vector<ProcedureDefinition> procedures_;
    std::vector<std::set<std::string, CaseInsensitiveLess>> incidenceMatrix_;
    
    // Helper functions
    void lowerCallExpressionOutputs(std::vector<StmtPtr>& statements);
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

// Escape a text for LaTeX text mode (table cells, comments): % # & _ { } $ \ ^ ~
std::string escapeLatexText(const std::string& s);

// Escape the characters of a name or string for LaTeX math mode: the same set, `|` and `\` as
// symbols, `°` and the accented letters of Latin-1 (\text{é}) so that pdflatex compiles them
std::string escapeLatexMath(const std::string& s);

// Escape a text for a cell of a Markdown table: `|` would end the cell
std::string markdownTableCell(const std::string& s);

// Convert a function call to LaTeX (handles thermodynamic functions specially)
std::string functionCallToLatex(const FunctionCall& func);

}  // namespace coolsolve
