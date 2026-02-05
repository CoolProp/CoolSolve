#pragma once

#include "ast.h"
#include "ir.h"
#include "structural_analysis.h"
#include "autodiff_node.h"
#include "units.h"
#include <map>
#include <string>
#include <vector>
#include <functional>

namespace coolsolve {

// ============================================================================
// Profiling API
// ============================================================================
void resetProfilingStats();
std::string getProfilingStatsString();

// ============================================================================
// CoolProp Configuration API
// ============================================================================
/// Apply CoolProp global configuration settings (call before any CoolProp calls)
void applyCoolPropConfig(const struct CoolPropConfig& config);

// ============================================================================
// Evaluation Result
// ============================================================================

/**
 * @brief Result of evaluating a block of equations.
 * 
 * Contains:
 * - residuals: F(x) values for each equation (should be zero at solution)
 * - jacobian: The Jacobian matrix J[i][j] = dF_i / dx_j
 */
struct EvaluationResult {
    std::vector<double> residuals;           // F(x) for each equation
    std::vector<std::vector<double>> jacobian;  // J[eq][var] = dF_eq/dx_var
    std::vector<std::string> warnings;       // Any warnings during evaluation
    
    size_t numEquations() const { return residuals.size(); }
    size_t numVariables() const { return jacobian.empty() ? 0 : jacobian[0].size(); }
};

// ============================================================================
// CoolProp Backend Configuration
// ============================================================================

/**
 * @brief Configuration for CoolProp backend selection.
 * 
 * Allows switching between:
 * - Full EOS: "HEOS" (default, most accurate)
 * - Interpolation: "BICUBIC&HEOS" or "TTSE&HEOS" (faster, less accurate)
 */
struct CoolPropConfig {
    std::string backend = "HEOS";  // Default to full Helmholtz EOS
    bool cacheEnabled = true;       // Cache AbstractState instances
    double cacheTolerance = 1e-10;  // Re-use cached state if inputs change less than this
    bool enableSuperancillaries = true;  // Enable CoolProp superancillary functions (faster VLE but more init time)
    
    // Factory function string for creating AbstractState
    std::string getBackendString() const {
        return backend;
    }
    
    // Unit system configuration
    UnitSystem units;
};

// ============================================================================
// Expression Evaluator
// ============================================================================

/**
 * @brief Evaluates AST expressions with automatic differentiation.
 * 
 * This class traverses the expression tree and computes both values and 
 * gradients using forward-mode automatic differentiation.
 */
class ExpressionEvaluator {
public:
    ExpressionEvaluator(size_t numVariables, const CoolPropConfig& config = CoolPropConfig());
    
    /**
     * @brief Set the value of a variable for evaluation.
     * @param name Variable name
     * @param value ADValue containing value and gradient
     */
    void setVariable(const std::string& name, const ADValue& value);
    
    /**
     * @brief Set the value of a string variable.
     * @param name Variable name
     * @param value String value
     */
    void setStringVariable(const std::string& name, const std::string& value);
    
    /**
     * @brief Get the current value of a variable.
     * @param name Variable name
     * @return The ADValue for this variable
     */
    ADValue getVariable(const std::string& name) const;
    
    /**
     * @brief Get the list of variables that were not found and used a default value.
     */
    const std::set<std::string, CaseInsensitiveLess>& getMissingVariables() const { return missingVariables_; }
    
    /**
     * @brief Clear the list of missing variables.
     */
    void clearMissingVariables() { missingVariables_.clear(); }
    
    /**
     * @brief Get the current value of a string variable.
     * @param name Variable name
     * @return The string value
     */
    std::string getStringVariable(const std::string& name) const;
    
    /**
     * @brief Check if a variable is set.
     */
    bool hasVariable(const std::string& name) const;
    
    /**
     * @brief Check if a string variable is set.
     */
    bool hasStringVariable(const std::string& name) const;
    
    /**
     * @brief Evaluate an expression tree.
     * @param expr The expression to evaluate
     * @return ADValue with computed value and gradient
     */
    ADValue evaluate(const ExprPtr& expr);
    
    /**
     * @brief Evaluate a string expression (literal or variable).
     */
    std::string evaluateString(const ExprPtr& expr);
    
    /**
     * @brief Clear all variable values (but keep configuration).
     */
    void clear();
    
    /**
     * @brief Register a user-defined function.
     */
    void registerFunction(const FunctionDefinition& func);
    
    /**
     * @brief Register a user-defined procedure.
     */
    void registerProcedure(const ProcedureDefinition& proc);
    
    /**
     * @brief Resolve a variable name, including evaluating indices for array access.
     */
    std::string resolveVariableName(const Variable& var);
    
    /**
     * @brief Evaluate a procedure call.
     */
    void evaluateProcedureCall(const ProcedureCall& call);
    
    /**
     * @brief Get the number of variables this evaluator is configured for.
     */
    size_t getNumVariables() const { return numVariables_; }
    
private:
    size_t numVariables_;
    std::map<std::string, ADValue, CaseInsensitiveLess> variables_;
    std::map<std::string, std::string, CaseInsensitiveLess> stringVariables_;
    std::map<std::string, FunctionDefinition, CaseInsensitiveLess> userFunctions_;
    std::map<std::string, ProcedureDefinition, CaseInsensitiveLess> userProcedures_;
    mutable std::set<std::string, CaseInsensitiveLess> missingVariables_;
    CoolPropConfig coolpropConfig_;
    
    // Helper methods for evaluating different expression types
    ADValue evaluateNumber(const NumberLiteral& num);
    ADValue evaluateVariable(const Variable& var);
    ADValue evaluateUnaryOp(const UnaryOp& op);
    ADValue evaluateBinaryOp(const BinaryOp& op);
    ADValue evaluateFunctionCall(const FunctionCall& func);
    
    // User-defined function evaluation
    ADValue evaluateUserFunction(const FunctionDefinition& func, const FunctionCall& call);
    
    // CoolProp function evaluation
    ADValue evaluateCoolPropFunction(const FunctionCall& func);
};

// ============================================================================
// Block Evaluator
// ============================================================================

/**
 * @brief Evaluates a block of equations and computes residuals and Jacobian.
 * 
 * A Block is a set of equations that must be solved simultaneously (an algebraic loop).
 * The BlockEvaluator:
 * 1. Takes the current state vector x
 * 2. Computes residuals F(x) = LHS - RHS for each equation
 * 3. Computes the Jacobian matrix J = dF/dx using forward-mode AD
 */
class BlockEvaluator {
public:
    /**
     * @brief Construct a block evaluator.
     * @param block The block definition from structural analysis
     * @param ir The IR containing equations and variable info
     * @param config CoolProp configuration
     */
    BlockEvaluator(const Block& block, const IR& ir, 
                   const CoolPropConfig& config = CoolPropConfig());
    
    /**
     * @brief Evaluate the block at the given state.
     * @param x Values for the block's variables (in order of block.variables)
     * @param externalVars Values for variables external to this block
     * @param externalStringVars Values for string variables external to this block
     * @return Residuals and Jacobian
     */
    EvaluationResult evaluate(const std::vector<double>& x,
                              const std::map<std::string, double>& externalVars = {},
                              const std::map<std::string, std::string>& externalStringVars = {});
    
    /**
     * @brief Get the variable names in this block.
     */
    const std::vector<std::string>& getVariables() const { return variables_; }
    
    /**
     * @brief Get the equation IDs in this block.
     */
    const std::vector<int>& getEquationIds() const { return equationIds_; }
    
    /**
     * @brief Get the number of equations/variables in this block.
     */
    size_t size() const { return variables_.size(); }
    
    /**
     * @brief Check if this block is explicit (can be solved directly).
     */
    bool isExplicit() const { return equationIds_.size() == 1; }
    
    /**
     * @brief Register a user-defined function for this block.
     */
    void registerFunction(const FunctionDefinition& func) { functions_[func.name] = func; }
    
    /**
     * @brief Register a user-defined procedure for this block.
     */
    void registerProcedure(const ProcedureDefinition& proc) { procedures_[proc.name] = proc; }
    
private:
    std::vector<std::string> variables_;
    std::vector<int> equationIds_;
    std::vector<const EquationInfo*> equations_;
    std::map<std::string, FunctionDefinition, CaseInsensitiveLess> functions_;
    std::map<std::string, ProcedureDefinition, CaseInsensitiveLess> procedures_;
    CoolPropConfig config_;
    
    // Create mapping from variable name to index in block
    std::map<std::string, size_t> varToIndex_;
};

// ============================================================================
// System Evaluator
// ============================================================================

/**
 * @brief Evaluates the complete equation system using the block structure.
 * 
 * The SystemEvaluator orchestrates the solution of all blocks in topological order.
 * For explicit blocks (size 1), it evaluates directly.
 * For algebraic loops, it provides the interface for Newton iteration.
 */
class SystemEvaluator {
public:
    /**
     * @brief Construct a system evaluator.
     * @param ir The IR containing the full equation system
     * @param analysisResult The structural analysis result with blocks
     * @param config CoolProp configuration
     */
    SystemEvaluator(const IR& ir, 
                    const StructuralAnalysisResult& analysisResult,
                    const CoolPropConfig& config = CoolPropConfig());
    
    /**
     * @brief Get the number of blocks.
     */
    size_t getNumBlocks() const { return blockEvaluators_.size(); }
    
    /**
     * @brief Get a block evaluator by index.
     */
    const BlockEvaluator& getBlock(size_t index) const { return blockEvaluators_[index]; }
    BlockEvaluator& getBlock(size_t index) { return blockEvaluators_[index]; }
    
    /**
     * @brief Set the current solution state for a variable.
     */
    void setVariableValue(const std::string& name, double value);
    
    /**
     * @brief Set the current solution state for a string variable.
     */
    void setStringVariableValue(const std::string& name, const std::string& value);
    
    /**
     * @brief Get the current solution state for a variable.
     */
    double getVariableValue(const std::string& name) const;
    
    /**
     * @brief Get the current solution state for a string variable.
     */
    std::string getStringVariableValue(const std::string& name) const;
    
    /**
     * @brief Get all current variable values.
     */
    const std::map<std::string, double, CaseInsensitiveLess>& getAllVariables() const { return currentState_; }
    
    /**
     * @brief Get all current string variable values.
     */
    const std::map<std::string, std::string, CaseInsensitiveLess>& getAllStringVariables() const { return currentStringState_; }
    
    /**
     * @brief Initialize all variables with guess values from the IR.
     */
    void initializeFromGuesses();
    
    /**
     * @brief Evaluate a specific block with current state.
     * @param blockIndex The block to evaluate
     * @return Evaluation result for that block
     */
    EvaluationResult evaluateBlock(size_t blockIndex);
    
private:
    const IR& ir_;
    std::vector<BlockEvaluator> blockEvaluators_;
    std::map<std::string, double, CaseInsensitiveLess> currentState_;
    std::map<std::string, std::string, CaseInsensitiveLess> currentStringState_;
    CoolPropConfig config_;
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Compare Jacobian computed by AD with numerical finite differences.
 * 
 * Useful for testing that the AD implementation is correct.
 * 
 * @param evaluator The block evaluator
 * @param x Current state
 * @param externalVars External variable values
 * @param externalStringVars External string variable values
 * @param epsilon Finite difference step size
 * @param verbose If true, prints both Jacobians to stdout for comparison
 * @return Maximum absolute difference between AD and numerical Jacobian
 */
double compareJacobianWithFiniteDifferences(
    BlockEvaluator& evaluator,
    const std::vector<double>& x,
    const std::map<std::string, double>& externalVars = {},
    const std::map<std::string, std::string>& externalStringVars = {},
    double epsilon = 1e-7,
    bool verbose = false);

/**
 * @brief Generate a summary report of the evaluation system.
 */
std::string generateEvaluatorReport(const SystemEvaluator& evaluator);

}  // namespace coolsolve
