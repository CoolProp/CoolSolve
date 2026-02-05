#include "coolsolve/runner.h"
#include "coolsolve/variable_inference.h"
#include "coolsolve/evaluator.h" // For getProfilingStatsString
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <chrono>

namespace fs = std::filesystem;

namespace coolsolve {

CoolSolveRunner::CoolSolveRunner(const std::string& inputFile)
    : inputFile_(inputFile) {}

bool CoolSolveRunner::run(const SolverOptions& options, bool enableTracing) {
    auto pipeline_start = std::chrono::high_resolution_clock::now();
    
    // 1. Parse
    auto t1 = std::chrono::high_resolution_clock::now();
    EESParser parser;
    parseResult_ = parser.parseFile(inputFile_);
    auto t2 = std::chrono::high_resolution_clock::now();
    timing_.parse_time_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    if (!parseResult_.success && parseResult_.equationCount == 0) {
        return false;
    }

    // 2. Build IR
    try {
        t1 = std::chrono::high_resolution_clock::now();
        ir_ = std::make_unique<IR>(IR::fromAST(parseResult_.program));
        t2 = std::chrono::high_resolution_clock::now();
        timing_.ir_time_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    } catch (...) {
        return false;
    }

    // 3. Infer & Initialize Variables
    try {
        t1 = std::chrono::high_resolution_clock::now();
        inferVariables(*ir_);
        initializeVariables(*ir_);
        t2 = std::chrono::high_resolution_clock::now();
        timing_.infer_time_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    } catch (...) {
        // Warning?
    }

    // 4. Structural Analysis
    try {
        t1 = std::chrono::high_resolution_clock::now();
        analysisResult_ = StructuralAnalyzer::analyze(*ir_);
        t2 = std::chrono::high_resolution_clock::now();
        timing_.analysis_time_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
        if (!analysisResult_.success) return false;
    } catch (...) {
        return false;
    }

    // 5. Solve
    try {
        // Load initials if present BEFORE constructing solver
        fs::path initialsPath = fs::path(inputFile_);
        initialsPath.replace_extension(".initials");
        if (fs::exists(initialsPath)) {
            ir_->loadInitialsFromFile(initialsPath.string());
        }

        t1 = std::chrono::high_resolution_clock::now();
        Solver solver(*ir_, analysisResult_);
        solveResult_ = solver.solve(options, enableTracing);
        t2 = std::chrono::high_resolution_clock::now();
        timing_.solve_time_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    } catch (...) {
        return false;
    }

    auto pipeline_end = std::chrono::high_resolution_clock::now();
    timing_.total_time_ms = std::chrono::duration<double, std::milli>(pipeline_end - pipeline_start).count();
    
    return solveResult_.success;
}

int CoolSolveRunner::loadInitials(const std::string& initialsPath) {
    if (ir_) {
        return ir_->loadInitialsFromFile(initialsPath);
    }
    return 0;
}

// Helper to write string to file
static void writeFile(const fs::path& path, const std::string& content) {
    std::ofstream file(path);
    if (file.is_open()) {
        file << content;
    }
}

void CoolSolveRunner::generateDebugOutput(const std::string& debugDirStr, const std::string& sourceCode) {
    fs::path debugDir(debugDirStr);
    try {
        fs::create_directories(debugDir);
    } catch (const std::exception& e) {
        std::cerr << "Error: Could not create debug directory '" << debugDirStr << "': " << e.what() << "\n";
        return;
    }
    
    // 1. Copy original source code
    writeFile(debugDir / "original.eescode", sourceCode);
    
    // 2. Generate model statistics
    std::ostringstream stats;
    stats << "# CoolSolve Analysis Report\n\n";
    stats << "**Input file:** " << inputFile_ << "\n\n";
    
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    stats << "**Generated:** " << std::ctime(&time) << "\n";
    
    stats << "## Model Statistics\n\n";
    stats << "| Metric | Value |\n";
    stats << "|--------|-------|\n";
    stats << "| Total lines | " << parseResult_.totalLines << " |\n";
    stats << "| Equations | " << parseResult_.equationCount << " |\n";
    stats << "| Comments | " << parseResult_.commentCount << " |\n";
    stats << "| Variables | " << ir_->getAlgebraicVariableCount() << " |\n";
    stats << "| System square | " << (ir_->isSquare() ? "Yes" : "No") << " |\n";
    stats << "| Total blocks | " << analysisResult_.totalBlocks << " |\n";
    stats << "| Largest block | " << analysisResult_.largestBlockSize << " |\n";
    stats << "| Explicit equations | " << analysisResult_.explicitEquationCount << " |\n";
    
    if (solveResult_.status != SolverStatus::InvalidInput) { // If solver ran
        stats << "| Solver Status | " << (solveResult_.success ? "SUCCESS" : "FAILED") << " |\n";
        stats << "| Solver Method | Newton-Raphson |\n";
        stats << "| Total Iterations | " << solveResult_.totalIterations << " |\n";
    }
    
    if (!analysisResult_.errorMessage.empty()) {
        stats << "\n## Warnings\n\n";
        stats << "- " << analysisResult_.errorMessage << "\n";
    }
    
    // Block summary
    stats << "\n## Block Summary\n\n";
    stats << "| Block | Size | Type | Variables |\n";
    stats << "|-------|------|------|----------|\n";
    for (const auto& block : analysisResult_.blocks) {
        stats << "| " << block.id;
        stats << " | " << block.size();
        stats << " | " << (block.isExplicit() ? "Explicit" : "Algebraic loop");
        stats << " | ";
        for (size_t i = 0; i < block.variables.size() && i < 3; ++i) {
            if (i > 0) stats << ", ";
            stats << block.variables[i];
        }
        if (block.variables.size() > 3) stats << ", ...";
        stats << " |\n";
    }
    
    writeFile(debugDir / "report.md", stats.str());

    // 2b. Solver Errors
    if (!solveResult_.success && solveResult_.status != SolverStatus::InvalidInput) {
        std::ostringstream errs;
        errs << "# Solver Errors\n\n";
        errs << "**Status:** " << coolsolve::statusToString(solveResult_.status) << "\n\n";
        errs << "**Error Message:** " << solveResult_.errorMessage << "\n\n";
        if (!solveResult_.detailedError.empty()) {
            errs << "## Detailed Error\n\n";
            errs << "```\n" << solveResult_.detailedError << "\n```\n";
        }
        writeFile(debugDir / "solver_errors.md", errs.str());
    }
    
    // 3. Variable mapping table
    std::ostringstream varTable;
    varTable << "# Variable Mapping Table\n\n";
    varTable << "| # | Variable Name | LaTeX | Solution | Index | Units | Inferred |\n";
    varTable << "|---|---------------|-------|----------|-------|-------|----------|\n";
    int varNum = 1;
    for (const auto& [name, info] : ir_->getVariables()) {
        varTable << "| " << varNum++;
        varTable << " | `" << name << "`";
        varTable << " | $" << info.latexLabel << "$";
        
        std::string valStr = "-";
        if (solveResult_.variables.find(name) != solveResult_.variables.end()) {
             std::stringstream ss;
             ss << solveResult_.variables.at(name);
             valStr = ss.str();
        } else if (info.solutionValue) {
             std::stringstream ss;
             ss << *info.solutionValue << " (initial)";
             valStr = ss.str();
        } else if (info.guessValue) {
             std::stringstream ss;
             ss << *info.guessValue << " (guess)";
             valStr = ss.str();
        }
        
        varTable << " | " << valStr;
        varTable << " | " << info.index;
        varTable << " | " << info.units;
        varTable << " | " << info.inferredProperty << " (" << info.inferredFluid << ")";
        varTable << " |\n";
    }
    writeFile(debugDir / "variables.md", varTable.str());
    
    // 4. Variables CSV
    writeFile(debugDir / "variables.csv", ir_->variableTableToCSV());
    
    // 5. Equations
    std::ostringstream eqList;
    eqList << "# Equations by Block\n\n";
    std::map<int, std::vector<const coolsolve::EquationInfo*>> blockEquations;
    for (const auto& eq : ir_->getEquations()) {
        blockEquations[eq.blockId].push_back(&eq);
    }
    for (const auto& [blockId, eqs] : blockEquations) {
        eqList << "## Block " << blockId;
        if (eqs.size() == 1) eqList << " (Explicit)\n\n";
        else eqList << " (Algebraic loop, size " << eqs.size() << ")\n\n";
        for (const auto* eq : eqs) {
            eqList << "- **Eq " << eq->id << "** (line " << eq->sourceLine << "): ";
            eqList << "`" << eq->originalText << "`";
            if (!eq->outputVariable.empty()) eqList << " → `" << eq->outputVariable << "`";
            eqList << "\n";
        }
        eqList << "\n";
    }
    writeFile(debugDir / "equations.md", eqList.str());
    
    // 6. JSON
    writeFile(debugDir / "analysis.json", coolsolve::generateAnalysisJSON(*ir_, analysisResult_));
    
    // 7. Residuals (Text)
    std::string eesResiduals = coolsolve::generateResidualsReport(*ir_, analysisResult_);
    writeFile(debugDir / "ees_residuals.txt", 
        "NOTE: This file is formatted for compatibility with EES residuals reports.\n"
        "It does NOT contain actual coolsolve residuals.\n\n" + eesResiduals);
    
    // 8. Residuals (Markdown) & Evaluator
    std::ostringstream csRes;
    csRes << "# CoolSolve Residuals Report\n\n";
    csRes << "| Block | Equation | LHS | RHS | Residual |\n";
    csRes << "|-------|----------|-----|-----|----------|\n";
    
    // Create evaluator for reporting
    coolsolve::SystemEvaluator csEval(*ir_, analysisResult_);
    // Use solve result if available
    if (solveResult_.success) {
        for (const auto& [n, v] : solveResult_.variables) csEval.setVariableValue(n, v);
        for (const auto& [n, v] : solveResult_.stringVariables) csEval.setStringVariableValue(n, v);
    } else {
        csEval.initializeFromGuesses(); // From IR
    }
    
    std::map<int, std::vector<const coolsolve::EquationInfo*>> blockEqs;
    for (const auto& eq : ir_->getEquations()) blockEqs[eq.blockId].push_back(&eq);
    
    for (const auto& [blockId, eqs] : blockEqs) {
        bool evalSuccess = false;
        coolsolve::EvaluationResult evalResult;
        try {
            evalResult = csEval.evaluateBlock(blockId);
            evalSuccess = true;
        } catch (...) { evalSuccess = false; }
        
        for (size_t i = 0; i < eqs.size(); ++i) {
            const auto* eq = eqs[i];
            csRes << "| " << blockId << " | `" << eq->originalText << "`";
            if (evalSuccess && i < evalResult.residuals.size()) {
                try {
                    coolsolve::ExpressionEvaluator exprEval(ir_->getVariableCount());
                    for (const auto& [name, val] : csEval.getAllVariables()) {
                        exprEval.setVariable(name, coolsolve::ADValue::constant(val, ir_->getVariableCount()));
                    }
                    for (const auto& [name, val] : csEval.getAllStringVariables()) {
                        exprEval.setStringVariable(name, val);
                    }
                    auto lhsVal = exprEval.evaluate(eq->lhs);
                    auto rhsVal = exprEval.evaluate(eq->rhs);
                    csRes << " | " << std::scientific << std::setprecision(4) << lhsVal.value;
                    csRes << " | " << rhsVal.value;
                    csRes << " | " << (lhsVal.value - rhsVal.value);
                } catch (...) {
                    csRes << " | - | - | " << evalResult.residuals[i];
                }
            } else {
                csRes << " | - | - | -";
            }
            csRes << " |\n";
        }
    }
    writeFile(debugDir / "coolsolve_residuals.md", csRes.str());
    
    // 9. LaTeX
    writeFile(debugDir / "equations.tex", ir_->toLatex());
    
    // 10. Incidence
    std::ostringstream incidence;
    incidence << "# Incidence Matrix\n\n";
    const auto& incMatrix = ir_->getIncidenceMatrix();
    for (size_t i = 0; i < incMatrix.size(); ++i) {
        incidence << "**Eq " << i << "**: ";
        bool first = true;
        for (const auto& var : incMatrix[i]) {
            if (!first) incidence << ", ";
            incidence << "`" << var << "`";
            first = false;
        }
        incidence << "\n";
    }
    writeFile(debugDir / "incidence.md", incidence.str());
    
    // 11. Evaluator Report
    writeFile(debugDir / "evaluator.md", coolsolve::generateEvaluatorReport(csEval));

    // 12. Profiling Report
    std::ostringstream prof;
    prof << "# Profiling Report\n\n";
    prof << "## Pipeline Timing\n\n";
    prof << "| Stage | Time (ms) |\n|---|---|\n";
    prof << "| Parse | " << timing_.parse_time_ms << " |\n";
    prof << "| IR Build | " << timing_.ir_time_ms << " |\n";
    prof << "| Variable Inference | " << timing_.infer_time_ms << " |\n";
    prof << "| Structural Analysis | " << timing_.analysis_time_ms << " |\n";
    prof << "| Solver | " << timing_.solve_time_ms << " |\n";
    prof << "| **Total** | **" << timing_.total_time_ms << "** |\n";
    
    prof << "\n## CoolProp Statistics\n\n";
    prof << "```\n" << coolsolve::getProfilingStatsString() << "\n```\n";
    writeFile(debugDir / "profiling.md", prof.str());

    // 13. Solver Trace
    if (!solveResult_.blockTraces.empty()) {
        std::ostringstream trace;
        trace << "# Solver Trace\n\n";
        for (size_t i = 0; i < solveResult_.blockTraces.size(); ++i) {
            const auto& tr = solveResult_.blockTraces[i];
            if (tr.iterations.empty()) continue;
            trace << "## Block " << i << "\n\n";
            trace << "```\n" << tr.toString() << "\n```\n\n";
        }
        writeFile(debugDir / "solver_trace.md", trace.str());
    }
    
    // 14. Index
    std::ostringstream index;
    index << "# CoolSolve Debug Output\n\n";
    index << "Analysis of: `" << inputFile_ << "`\n\n";
    index << "## Contents\n\n";
    index << "| File | Description |\n|---|---|\n";
    index << "| [report.md](report.md) | Model statistics and block summary |\n";
    index << "| [variables.md](variables.md) | Variable mapping table |\n";
    index << "| [coolsolve_residuals.md](coolsolve_residuals.md) | Detailed equation values and residuals |\n";
    index << "| [evaluator.md](evaluator.md) | Evaluator structure and test |\n";
    index << "| [equations.md](equations.md) | Equations grouped by block |\n";
    index << "| [incidence.md](incidence.md) | Variable-equation incidence matrix |\n";
    index << "| [equations.tex](equations.tex) | LaTeX formatted equations |\n";
    index << "| [profiling.md](profiling.md) | Performance profiling and stats |\n";
    if (!solveResult_.blockTraces.empty()) {
        index << "| [solver_trace.md](solver_trace.md) | Detailed solver iteration trace |\n";
    }
    writeFile(debugDir / "README.md", index.str());
}

} // namespace coolsolve
