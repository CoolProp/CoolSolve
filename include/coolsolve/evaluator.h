#pragma once

#include "ast.h"
#include "ir.h"
#include "lookup_table.h"
#include "structural_analysis.h"
#include "autodiff_node.h"
#include "units.h"
#include "diagnostic.h"
#include <map>
#include <string>
#include <vector>
#include <functional>

namespace coolsolve {

// ============================================================================
// Profiling API
// ============================================================================
void resetProfilingStats();
void resetCoolPropUnitHints();
std::string getProfilingStatsString();

/**
 * @brief Perform a one-time CoolProp warmup call and return its duration.
 *
 * This triggers CoolProp's internal initialization (fluid list loading, etc.)
 * so that subsequent property evaluations are not penalized by the first-call
 * overhead. The warmup is safe to call multiple times; only the first call
 * performs work and records timing.
 *
 * @return Warmup duration in milliseconds (0.0 if already initialized or on failure).
 */
double warmupCoolProp();

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
    /// Values of the string variables that a CALL of the block set (name -> value): a string is not
    /// an unknown of the numeric system, so this is the only way it leaves the evaluation.
    std::map<std::string, std::string> stringOutputs;
    
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
    
    // Low-level API: use AbstractState instead of PropsSI
    // When true, property evaluations go through cached AbstractState objects,
    // eliminating string parsing and fluid lookup overhead on every call.
    // When false, falls back to the high-level PropsSI interface (zero overhead).
    bool useAbstractState = true;
    
    // Analytical derivatives: use CoolProp's first_partial_deriv() instead of FD
    // Requires useAbstractState=true. When true, each property evaluation uses
    // exact analytical derivatives (1 call instead of 5). When false, central
    // finite differences are used (the legacy behaviour; zero overhead).
    bool enableAnalyticalDerivatives = true;
    
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
     * @brief Evaluate a string expression (literal, variable or string function: STRING$, PHASE$, LOOKUP$).
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
    
    /**
     * @brief Set residual-only mode (skip CoolProp derivative computations).
     */
    void setResidualOnly(bool residualOnly) { residualOnly_ = residualOnly; }
    bool isResidualOnly() const { return residualOnly_; }
    
    void setDiagnostics(DiagnosticCollector* diag) { diagnostics_ = diag; }

    /**
     * @brief Where the `CALL WARNING` statements executed by this evaluator (and by the bodies of the
     *        routines it calls) put their warnings; null (the default) makes them no-ops.
     *
     * A warning is meaningful only at the solution, not at every Newton iterate, so only the evaluation
     * of a block at its solution gives a sink (`SystemEvaluator::reportCallWarnings`): the evaluations
     * of the iterations skip the statement without even formatting its message.
     */
    void setCallWarningSink(std::vector<Diagnostic>* sink) { callWarningSink_ = sink; }

    /** @brief Provide access to lookup tables for INTERPOLATE / LOOKUP calls. */
    void setLookupTableStore(const LookupTableStore* store) { lookupTableStore_ = store; }
    const LookupTableStore* getLookupTableStore() const { return lookupTableStore_; }

    /// When true, skip input clamping (P >= 1000 Pa, T >= 50 K, etc.).
    /// Used for final solution verification so that out-of-range inputs
    /// are caught as errors rather than silently corrected.
    void setDisableClamping(bool disable) { disableClamping_ = disable; }
    bool isClampingDisabled() const { return disableClamping_; }
    
private:
    size_t numVariables_;
    bool residualOnly_ = false;
    bool disableClamping_ = false;
    DiagnosticCollector* diagnostics_ = nullptr;
    std::vector<Diagnostic>* callWarningSink_ = nullptr;   ///< see setCallWarningSink
    const LookupTableStore* lookupTableStore_ = nullptr;
    /// "FUNCTION 'name'" / "PROCEDURE 'name'" when this evaluator runs the body of a user routine
    /// (empty in the main program); only used to say where a lookup error comes from.
    std::string routine_;
    /// The name of that routine (empty in the main program): the YYY of the generic messages of
    /// CALL WARNING and CALL ERROR.
    std::string routineName_;
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

    /**
     * @brief The text of a `CALL WARNING` / `CALL ERROR` (EES help, WARNING and ERROR procedures).
     *
     * The call is `(message, value)`, `(message)` or `(value)`; the message is a string constant,
     * variable or function, the value a number (or, with `XXX$`, a string). The characters `XXX` of the
     * message, with the display format that may follow (`XXXA1`, `XXXF1`, `XXXE4`, see `eesFormatNumber`),
     * are replaced by the value. Without a message EES writes a generic text naming the routine. A message
     * that has no `XXX` but a value gets the value appended in brackets, so that it is not lost.
     * @param isError  true for CALL ERROR (its generic text), false for CALL WARNING.
     */
    std::string userMessage(const ProcedureCall& call, bool isError);

    /// Evaluator for the body of the user routine `kind` (FUNCTION / PROCEDURE) `name` called from
    /// here. Variables are local to a body (only the arguments cross the boundary); what it
    /// inherits is the context every expression needs: the user routines and the lookup tables.
    ExpressionEvaluator makeBodyEvaluator(const char* kind, const std::string& name) const;
    
    // CoolProp function evaluation
    ADValue evaluateCoolPropFunction(const FunctionCall& func);

    // Lookup table function evaluation (INTERPOLATE, LOOKUP, etc.)
    ADValue evaluateLookupFunction(const FunctionCall& func);

    /// LOOKUP$('table', row, col): the text of a cell (a string-valued function, see evaluateString).
    /// The row is rounded to the nearest integer; the column is a name or a number.
    std::string evaluateLookupString(const FunctionCall& func);

    /// The table `tableName` for the lookup function `func`; an error (and a diagnostic) if there is none.
    const LookupTable& lookupTableFor(const FunctionCall& func, const std::string& tableName);

    /// The 1-based column given by argument `argIdx` of the lookup function `func`: a name or a number
    /// (0 when the number is not positive).
    size_t lookupColumnOf(const FunctionCall& func, const LookupTable& tbl, size_t argIdx);

    /// " (in FUNCTION 'f')" inside a body, empty in the main program: where a lookup error comes from.
    std::string lookupWhere() const;

    /// The fluid named by the first argument of a property function: a string literal, a string variable,
    /// an unquoted name (`Water`) or a string function such as LOOKUP$('table', i, 1). Empty if `arg` is none of these.
    std::string fluidNameOf(const ExprPtr& arg);
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
    /**
     * @param computeJacobian When false, skip derivative computation
     *        (saves CoolProp calls during line search; residuals only).
     * @param callWarnings When given, the `CALL WARNING` statements executed by the evaluation push their
     *        warnings (code `W001`) there; when null they do nothing (the case of every iteration).
     */
    EvaluationResult evaluate(const std::vector<double>& x,
                              const std::map<std::string, double>& externalVars = {},
                              const std::map<std::string, std::string>& externalStringVars = {},
                              bool computeJacobian = true,
                              std::vector<Diagnostic>* callWarnings = nullptr);
    
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
    
    void setDiagnostics(DiagnosticCollector* diag) { diagnostics_ = diag; }

    /** @brief Provide lookup table access for INTERPOLATE / LOOKUP calls. */
    void setLookupTableStore(const LookupTableStore* store) { lookupTableStore_ = store; }

    /** @brief True if a CALL of this block has a string output (`EvaluationResult::stringOutputs`). */
    bool hasStringCallOutputs() const { return hasStringCallOutputs_; }

    /**
     * @brief True if an equation of this block calls a routine that may execute `CALL WARNING`
     *        (see `markCallWarnings`); false for every block of a model without that call.
     */
    bool mayCallWarning() const { return mayCallWarning_; }

    /**
     * @brief Record whether the equations of this block call one of `routines`, the FUNCTIONs and
     *        PROCEDUREs whose body can reach a `CALL WARNING` (directly or through other routines).
     */
    void markCallWarnings(const std::set<std::string, CaseInsensitiveLess>& routines);
    
private:
    bool hasStringCallOutputs_ = false;
    bool mayCallWarning_ = false;
    /// Order in which `evaluate` visits the equations when a CALL of the block computes a string: those
    /// CALLs first, so that every CALL of the block finds the strings it reads. Empty: the order of
    /// `equations_`. The row of an equation in the residual vector and the Jacobian never changes.
    std::vector<size_t> evalOrder_;
    std::vector<std::string> variables_;
    std::vector<int> equationIds_;
    DiagnosticCollector* diagnostics_ = nullptr;
    const LookupTableStore* lookupTableStore_ = nullptr;
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
     * @brief Store the string variables that the CALLs of a solved block compute.
     *
     * A string is not an unknown of the numeric system, so the solver does not carry it from one
     * block to the next: the CALL that sets it is evaluated once more at the solution and the
     * strings it returns become part of the state (CS-BUG-STRING-CALL-OUT). Does nothing, and
     * costs nothing, for a block without a CALL that has a string output.
     */
    void publishStringOutputs(size_t blockIndex);
    
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

    /** @brief Provide lookup table access to all block evaluators. */
    void setLookupTableStore(const LookupTableStore* store) {
        for (auto& be : blockEvaluators_) be.setLookupTableStore(store);
    }
    
    /**
     * @brief Report the `CALL WARNING` that holds at the solution of a block (code `W001`).
     *
     * A user warning describes the solution (an input outside the validity range of a correlation), not
     * the Newton iterates that may cross the range on the way: it is raised by one more evaluation of
     * the block at its solution, made after the block converged, and not by the evaluations of the
     * iterations. A warning identical to one already reported (same line, same text) is not repeated.
     * Does nothing, and costs nothing, for a block that cannot reach a `CALL WARNING`
     * (`BlockEvaluator::mayCallWarning`), that is for every block of a model without that call.
     * @param blockIndex The block, solved
     * @param out Receives the warnings
     */
    void reportCallWarnings(size_t blockIndex, DiagnosticCollector& out);

    /**
     * @brief Evaluate a specific block with current state.
     * @param blockIndex The block to evaluate
     * @param computeJacobian false: residuals only
     * @param callWarnings Receives the `CALL WARNING` raised by the evaluation, none if null
     * @return Evaluation result for that block
     */
    EvaluationResult evaluateBlock(size_t blockIndex, bool computeJacobian = true,
                                   std::vector<Diagnostic>* callWarnings = nullptr);
    
private:
    const IR& ir_;
    std::vector<BlockEvaluator> blockEvaluators_;
    std::map<std::string, double, CaseInsensitiveLess> currentState_;
    std::map<std::string, std::string, CaseInsensitiveLess> currentStringState_;
    std::set<std::pair<int, std::string>> reportedCallWarnings_;   ///< (line, message) of the warnings already reported
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
