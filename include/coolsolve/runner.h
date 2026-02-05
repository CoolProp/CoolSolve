#pragma once

#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/solver.h"
#include <string>
#include <memory>

namespace coolsolve {

class CoolSolveRunner {
public:
    CoolSolveRunner(const std::string& inputFile);
    
    // Run the full pipeline
    bool run(const SolverOptions& options = SolverOptions(), bool enableTracing = false);
    
    // Generate debug output
    void generateDebugOutput(const std::string& debugDir, const std::string& sourceCode);
    
    struct PipelineTiming {
        double parse_time_ms = 0.0;
        double ir_time_ms = 0.0;
        double infer_time_ms = 0.0;
        double analysis_time_ms = 0.0;
        double solve_time_ms = 0.0;
        double total_time_ms = 0.0;
    };
    
    // Accessors
    const ParseResult& getParseResult() const { return parseResult_; }
    const IR& getIR() const { return *ir_; }
    IR& getIR() { return *ir_; } // Non-const accessor needed for loading initials
    const StructuralAnalysisResult& getAnalysisResult() const { return analysisResult_; }
    const SolveResult& getSolveResult() const { return solveResult_; }
    const PipelineTiming& getTiming() const { return timing_; }
    
    bool isParseSuccess() const { return parseResult_.success; }
    bool isIRSuccess() const { return ir_ != nullptr; }
    bool isAnalysisSuccess() const { return analysisResult_.success; }
    bool isSolveSuccess() const { return solveResult_.success; }

    // Helper to load initials (not solution)
    int loadInitials(const std::string& initialsPath);

private:
    std::string inputFile_;
    ParseResult parseResult_;
    std::unique_ptr<IR> ir_;
    StructuralAnalysisResult analysisResult_;
    SolveResult solveResult_;
    PipelineTiming timing_;
};

} // namespace coolsolve
