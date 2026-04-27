#pragma once

#include "ir.h"
#include "evaluator.h"
#include "lookup_table.h"
#include "diagnostic.h"
#include <map>
#include <string>
#include <vector>

namespace coolsolve {

/**
 * @brief Result of checking one equation against the proposed solution.
 */
struct EquationCheckResult {
    int equationId;
    std::string originalText;
    double lhsValue;
    double rhsValue;
    double residual;       ///< |LHS - RHS|
    double relativeError;  ///< |LHS - RHS| / max(|LHS|, |RHS|, 1)
    bool satisfied;        ///< true if relativeError <= tolerance
    bool isProcedureCall;  ///< true if equation is a CALL statement
};

/**
 * @brief Overall result of the solution verification.
 */
struct SolutionCheckResult {
    bool allSatisfied = true;
    size_t totalEquations = 0;
    size_t satisfiedCount = 0;
    size_t violatedCount = 0;
    size_t skippedCount = 0;       ///< Equations skipped (null LHS/RHS, secondary CALL outputs)
    double maxResidual = 0.0;
    double maxRelativeError = 0.0;
    std::string worstEquationText;
    int worstEquationId = -1;
    std::vector<EquationCheckResult> checks;
    DiagnosticCollector diagnostics;
};

/**
 * @brief Verify that a set of variable values satisfies all equations in the IR.
 *
 * For each equation, evaluates LHS and RHS independently and checks that
 * |LHS - RHS| / max(|LHS|, |RHS|, 1) <= tolerance.
 *
 * Procedure CALLs are evaluated by running the procedure and comparing
 * each output's old value with the value returned by the procedure.
 *
 * The verification tolerance should be significantly looser than the solver's
 * Newton convergence tolerance. The solver tolerance (e.g. 1e-6) controls how
 * tightly each block converges internally, but a full independent re-evaluation
 * of all equations amplifies small numerical differences — especially through
 * procedure interpolation (e.g. therminol66, pinch functions) and accumulated
 * floating-point noise across blocks. A verification tolerance of 1e-3 (0.1%
 * relative error) is appropriate for confirming the solution is physically valid.
 *
 * @param ir            The intermediate representation with all equations
 * @param variables     Numeric variable values (typically from SolveResult or .sol file)
 * @param stringVars    String variable values
 * @param config        CoolProp configuration to use for evaluation
 * @param tolerance     Relative error tolerance for verification (default: 1e-3)
 * @return SolutionCheckResult with per-equation details
 */
SolutionCheckResult checkSolution(
    const IR& ir,
    const std::map<std::string, double, CaseInsensitiveLess>& variables,
    const std::map<std::string, std::string, CaseInsensitiveLess>& stringVars,
    const CoolPropConfig& config = CoolPropConfig(),
    double tolerance = 1e-3,
    const LookupTableStore* lookupStore = nullptr);

/**
 * @brief Print a human-readable summary of the solution check to stdout.
 *
 * Shows violated equations with LHS, RHS, and residual values.
 *
 * @param result  The solution check result
 * @param verbose If true, print all equations (not just violations)
 */
void printSolutionCheckReport(const SolutionCheckResult& result, bool verbose = false);

/**
 * @brief Write a Markdown report of the solution check to a file.
 *
 * Writes a table of all checked equations with LHS, RHS, residual, and
 * relative error.  Violated equations are highlighted.
 *
 * @param path    Output file path (e.g. debugDir / "solution_check.md")
 * @param result  The solution check result
 */
void writeSolutionCheckReport(const std::string& path, const SolutionCheckResult& result);

} // namespace coolsolve
