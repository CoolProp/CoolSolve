#pragma once
/**
 * @file symbolic_reduction.h
 * @brief Symbolic Block Reduction — pre-processing pass that shrinks algebraic
 *        blocks before the iterative solver is invoked.
 *
 * Three reduction techniques are applied iteratively until a fixed point:
 *
 * 1. **Explicit extraction**: equations where the matched output variable can
 *    be computed directly from external (already-known) values are extracted
 *    from the block and evaluated before solving.
 *
 * 2. **CoolProp call inversion**: thermodynamic equations of the form
 *        output = thermo(fluid, A=a, B=b)
 *    are reformulated so that an unknown input becomes the output, if the
 *    original output and the other input are known.  CoolProp supports many
 *    input pairs (P-T, H-P, P-S, …) so this often succeeds.
 *
 * 3. **Equation substitution**: if a variable determined by one block equation
 *    appears in exactly one other block equation, the first equation's RHS is
 *    symbolically substituted into the second, eliminating one variable and one
 *    equation from the block.
 *
 * The result is a (possibly smaller) block for the iterative solver, plus
 * ordered lists of pre-solve and post-solve direct-evaluation steps.
 *
 * This feature is **optional** and **off by default**; when disabled, zero
 * overhead is added to the solving pipeline.
 */

#include "ir.h"
#include "structural_analysis.h"
#include <string>
#include <vector>
#include <set>
#include <map>

namespace coolsolve {

// ============================================================================
// Reduction result
// ============================================================================

/**
 * @brief One direct-evaluation step produced by the reduction pass.
 *
 * Each step computes exactly one variable either BEFORE the reduced block is
 * solved (pre-solve, its inputs are all external) or AFTER (post-solve, its
 * inputs may include variables solved by the reduced block).
 */
struct ReductionStep {
    std::string variable;          ///< Variable to compute
    int equationId;                ///< Original equation ID used
    bool inverted = false;         ///< True if this is an inverted CoolProp call

    // --- CoolProp inversion details (only when inverted == true) ---
    std::string invertedFuncName;  ///< New CoolProp function name (e.g. "temperature")
    std::string fluidName;         ///< Fluid argument (e.g. "Water")
    /// New named arguments: (param_name, variable_or_literal).
    /// The param_name can be "H", "P", "S", etc.  The value is the original
    /// variable/literal expression text for reconstruction.
    std::vector<std::pair<std::string, ExprPtr>> newNamedArgs;
    ExprPtr fluidExpr;             ///< Original fluid expression from args[0]
};

/**
 * @brief Result of the symbolic reduction analysis for a single block.
 */
struct BlockReductionResult {
    bool reduced = false;              ///< True if the block was reduced at all

    /// Steps to execute BEFORE the iterative solver (all inputs are external).
    std::vector<ReductionStep> preSolveSteps;

    /// Steps to execute AFTER the iterative solver (inputs may include solved vars).
    std::vector<ReductionStep> postSolveSteps;

    /// Remaining equation IDs that need iterative solving.
    std::vector<int> reducedEquationIds;

    /// Remaining variables that need iterative solving.
    std::vector<std::string> reducedVariables;

    /// Number of CoolProp inversions applied
    int inversionsApplied = 0;

    /// Number of explicit extractions applied
    int extractionsApplied = 0;

    /// Number of substitutions applied
    int substitutionsApplied = 0;
};

// ============================================================================
// Analysis entry point
// ============================================================================

/**
 * @brief Analyse a single block and determine possible symbolic reductions.
 *
 * This is a pure analysis function: it does not modify the IR, evaluator, or
 * block data structures.  The caller uses the result to (a) execute pre-solve
 * steps, (b) build a smaller sub-problem, and (c) execute post-solve steps.
 *
 * @param block          The block to analyse.
 * @param ir             The intermediate representation (equations + variables).
 * @param analysis       The structural analysis result (matching info).
 * @return BlockReductionResult  Reduction plan for this block.
 */
BlockReductionResult analyseBlockReduction(
    const Block& block,
    const IR& ir,
    const StructuralAnalysisResult& analysis);

// ============================================================================
// CoolProp inversion helpers (exposed for testing)
// ============================================================================

/**
 * @brief Map a CoolProp output function name to the corresponding input
 *        parameter name used in named arguments.
 *
 * Example: "enthalpy" → "H", "temperature" → "T", "entropy" → "S".
 * Returns empty string if the function is not invertible.
 */
std::string coolpropFuncToInputName(const std::string& funcName);

/**
 * @brief Map a CoolProp input parameter name to the corresponding output
 *        function name.
 *
 * Example: "H" → "enthalpy", "T" → "temperature", "S" → "entropy".
 * Returns empty string if no such function exists.
 */
std::string coolpropInputToFuncName(const std::string& inputName);

/**
 * @brief Check whether (param1, param2) form a valid CoolProp input pair.
 *
 * @param param1Name  First input name, e.g. "T", "H", "P" (case-insensitive).
 * @param param2Name  Second input name.
 * @return true if CoolProp can accept this pair as update inputs.
 */
bool isValidCoolPropInputPair(const std::string& param1Name,
                              const std::string& param2Name);

}  // namespace coolsolve
