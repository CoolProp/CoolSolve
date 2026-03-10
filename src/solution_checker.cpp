#include "coolsolve/solution_checker.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace coolsolve {

SolutionCheckResult checkSolution(
    const IR& ir,
    const std::map<std::string, double, CaseInsensitiveLess>& variables,
    const std::map<std::string, std::string, CaseInsensitiveLess>& stringVars,
    const CoolPropConfig& config,
    double tolerance)
{
    SolutionCheckResult result;
    const auto& equations = ir.getEquations();
    result.totalEquations = equations.size();

    // Set up a residual-only expression evaluator (no derivatives needed)
    ExpressionEvaluator exprEval(0, config);
    exprEval.setResidualOnly(true);

    // Register user-defined functions and procedures
    for (const auto& func : ir.getFunctions()) {
        exprEval.registerFunction(func);
    }
    for (const auto& proc : ir.getProcedures()) {
        exprEval.registerProcedure(proc);
    }

    // Load all variable values into the evaluator
    for (const auto& [name, value] : variables) {
        exprEval.setVariable(name, ADValue::constant(value, 0));
    }
    for (const auto& [name, value] : stringVars) {
        exprEval.setStringVariable(name, value);
    }

    // Evaluate each equation
    for (size_t eq = 0; eq < equations.size(); ++eq) {
        const EquationInfo& eqInfo = equations[eq];

        // Handle procedure CALLs
        if (eqInfo.procedureCall) {
            // Store old output values
            std::map<std::string, double> oldOutputs;
            for (const auto& var : eqInfo.procedureCall->outputVars) {
                std::string name = exprEval.resolveVariableName(var);
                oldOutputs[name] = exprEval.getVariable(name).value;
            }

            // Execute the procedure
            try {
                exprEval.evaluateProcedureCall(*eqInfo.procedureCall);
            } catch (const std::exception& e) {
                // If procedure fails, report all outputs as violated
                for (const auto& var : eqInfo.procedureCall->outputVars) {
                    std::string name = exprEval.resolveVariableName(var);
                    EquationCheckResult check;
                    check.equationId = eqInfo.id;
                    check.originalText = eqInfo.originalText + " [output: " + name + ", ERROR: " + e.what() + "]";
                    check.lhsValue = oldOutputs[name];
                    check.rhsValue = std::numeric_limits<double>::quiet_NaN();
                    check.residual = std::numeric_limits<double>::infinity();
                    check.relativeError = std::numeric_limits<double>::infinity();
                    check.satisfied = false;
                    check.isProcedureCall = true;
                    result.checks.push_back(check);
                    result.violatedCount++;
                    result.allSatisfied = false;
                }
                // Skip sibling CALL equations
                size_t nOut = eqInfo.procedureCall->outputVars.size();
                if (nOut > 1) {
                    for (size_t s = 0; s < nOut - 1 && eq + 1 < equations.size(); ++s) {
                        eq++;
                        result.skippedCount++;
                    }
                }
                continue;
            }

            // Compare old values with new procedure outputs
            for (const auto& var : eqInfo.procedureCall->outputVars) {
                std::string name = exprEval.resolveVariableName(var);
                double oldVal = oldOutputs[name];
                double newVal = exprEval.getVariable(name).value;
                double absResidual = std::abs(oldVal - newVal);
                double scale = std::max({std::abs(oldVal), std::abs(newVal), 1.0});
                double relError = absResidual / scale;

                EquationCheckResult check;
                check.equationId = eqInfo.id;
                check.originalText = eqInfo.originalText + " [output: " + name + "]";
                check.lhsValue = oldVal;
                check.rhsValue = newVal;
                check.residual = absResidual;
                check.relativeError = relError;
                check.satisfied = (relError <= tolerance);
                check.isProcedureCall = true;
                result.checks.push_back(check);

                if (check.satisfied) {
                    result.satisfiedCount++;
                } else {
                    result.violatedCount++;
                    result.allSatisfied = false;
                }
                if (absResidual > result.maxResidual) {
                    result.maxResidual = absResidual;
                    result.worstEquationText = check.originalText;
                    result.worstEquationId = eqInfo.id;
                }
                if (relError > result.maxRelativeError) {
                    result.maxRelativeError = relError;
                }

                // Restore old value so other equations still see the .sol value
                exprEval.setVariable(name, ADValue::constant(oldVal, 0));
            }

            // Skip the remaining N-1 sibling equations for this CALL.
            // Each CALL with N outputs creates N consecutive EquationInfo
            // entries, each with procedureCall set. We already checked all
            // outputs above, so skip the duplicates.
            size_t nOutputs = eqInfo.procedureCall->outputVars.size();
            if (nOutputs > 1) {
                size_t toSkip = nOutputs - 1;
                for (size_t s = 0; s < toSkip && eq + 1 < equations.size(); ++s) {
                    eq++;
                    result.skippedCount++;
                }
            }
            continue;
        }

        // Skip placeholder equations (secondary CALL outputs with null LHS/RHS)
        if (!eqInfo.lhs || !eqInfo.rhs) {
            result.skippedCount++;
            continue;
        }

        // Evaluate normal equation: LHS = RHS
        double lhsVal = 0.0, rhsVal = 0.0;
        try {
            lhsVal = exprEval.evaluate(eqInfo.lhs).value;
            rhsVal = exprEval.evaluate(eqInfo.rhs).value;
        } catch (const std::exception& e) {
            EquationCheckResult check;
            check.equationId = eqInfo.id;
            check.originalText = eqInfo.originalText + " [EVAL ERROR: " + e.what() + "]";
            check.lhsValue = lhsVal;
            check.rhsValue = rhsVal;
            check.residual = std::numeric_limits<double>::infinity();
            check.relativeError = std::numeric_limits<double>::infinity();
            check.satisfied = false;
            check.isProcedureCall = false;
            result.checks.push_back(check);
            result.violatedCount++;
            result.allSatisfied = false;
            continue;
        }

        double absResidual = std::abs(lhsVal - rhsVal);
        double scale = std::max({std::abs(lhsVal), std::abs(rhsVal), 1.0});
        double relError = absResidual / scale;

        EquationCheckResult check;
        check.equationId = eqInfo.id;
        check.originalText = eqInfo.originalText;
        check.lhsValue = lhsVal;
        check.rhsValue = rhsVal;
        check.residual = absResidual;
        check.relativeError = relError;
        check.satisfied = (relError <= tolerance);
        check.isProcedureCall = false;
        result.checks.push_back(check);

        if (check.satisfied) {
            result.satisfiedCount++;
        } else {
            result.violatedCount++;
            result.allSatisfied = false;
        }
        if (absResidual > result.maxResidual) {
            result.maxResidual = absResidual;
            result.worstEquationText = eqInfo.originalText;
            result.worstEquationId = eqInfo.id;
        }
        if (relError > result.maxRelativeError) {
            result.maxRelativeError = relError;
        }
    }

    return result;
}

void printSolutionCheckReport(const SolutionCheckResult& result, bool verbose) {
    std::cout << "\n=== Solution Verification ===\n";
    std::cout << "Total equations: " << result.totalEquations << "\n";
    std::cout << "Checked:   " << result.satisfiedCount + result.violatedCount << "\n";
    std::cout << "Satisfied: " << result.satisfiedCount << "\n";
    std::cout << "Violated:  " << result.violatedCount << "\n";
    std::cout << "Skipped:   " << result.skippedCount << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Max |residual|:     " << result.maxResidual << "\n";
    std::cout << "Max relative error: " << result.maxRelativeError << "\n";

    if (result.allSatisfied) {
        std::cout << "Result: ALL EQUATIONS SATISFIED\n";
    } else {
        std::cout << "Result: " << result.violatedCount << " EQUATION(S) VIOLATED\n";
        std::cout << "\nViolated equations:\n";
        std::cout << std::string(80, '-') << "\n";
        for (const auto& c : result.checks) {
            if (!c.satisfied) {
                std::cout << "  Eq " << c.equationId << ": " << c.originalText << "\n";
                std::cout << "    LHS = " << std::setprecision(10) << c.lhsValue
                          << "  RHS = " << c.rhsValue
                          << "  |res| = " << std::setprecision(4) << c.residual
                          << "  rel = " << c.relativeError << "\n";
            }
        }
        std::cout << std::string(80, '-') << "\n";
    }

    if (verbose) {
        std::cout << "\nAll equation checks:\n";
        std::cout << std::string(100, '-') << "\n";
        for (const auto& c : result.checks) {
            std::cout << "  " << (c.satisfied ? "OK  " : "FAIL")
                      << " Eq " << std::setw(3) << c.equationId
                      << ": LHS=" << std::setprecision(8) << std::setw(14) << c.lhsValue
                      << " RHS=" << std::setw(14) << c.rhsValue
                      << " |res|=" << std::setprecision(3) << std::setw(10) << c.residual
                      << " rel=" << std::setw(10) << c.relativeError
                      << "  " << c.originalText << "\n";
        }
        std::cout << std::string(100, '-') << "\n";
    }
}

} // namespace coolsolve
