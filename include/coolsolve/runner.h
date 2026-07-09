#pragma once

#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/solver.h"
#include "coolsolve/lookup_table.h"
#include "coolsolve/diagnostic.h"
#include "coolsolve/integral/integral_solver.h"
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
        double coolprop_warmup_time_ms = 0.0;
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

    /// True if the last run solved an equation-based dynamic (INTEGRAL) model.
    bool hasIntegralResult() const { return integralResult_.success || integralModel_; }
    /// Result of the last integral solve (trajectory table + diagnostics).
    const IntegralSolveResult& getIntegralResult() const { return integralResult_; }
    
    bool isParseSuccess() const { return parseResult_.success; }
    bool isIRSuccess() const { return ir_ != nullptr; }
    bool isAnalysisSuccess() const { return analysisResult_.success; }
    bool isSolveSuccess() const { return solveResult_.success; }

    // Helper to load initials (not solution)
    int loadInitials(const std::string& initialsPath);

    // LaTeX report generation.
    // Call after a successful solve.  Returns the LaTeX source; also stores
    // it internally so the server can serve it without regenerating.
    std::string generateLatexReportContent(const std::string& modelName) const;

    // Cached LaTeX source from the last generateLatexReportContent() call.
    const std::string& getLatexReportContent() const { return latexReportContent_; }

    // Whether a LaTeX report has been generated for the current solve.
    bool hasLatexReport() const { return !latexReportContent_.empty(); }

    // Aggregated diagnostics from all pipeline phases
    const DiagnosticCollector& getDiagnostics() const { return diagnostics_; }

    // Lookup table store loaded for the current run
    const LookupTableStore& getLookupTableStore() const { return lookupTableStore_; }

    /**
     * @brief Pre-supply a lookup table store before calling run().
     *
     * When set, run() skips disk-based CSV loading and uses this store
     * directly.  This is used by the GUI server, which manages its own
     * table store from the session's in-memory CSVs.  For CLI usage,
     * leave unset so that run() auto-loads the companion CSV file.
     */
    void setLookupTableStore(LookupTableStore store) {
        lookupTableStore_ = std::move(store);
        lookupTableStorePreloaded_ = true;
    }

private:
    std::string inputFile_;
    ParseResult parseResult_;
    std::unique_ptr<IR> ir_;
    StructuralAnalysisResult analysisResult_;
    SolveResult solveResult_;
    IntegralSolveResult integralResult_;
    bool integralModel_ = false;   ///< set when the parsed model contains INTEGRAL()
    PipelineTiming timing_;
    mutable std::string latexReportContent_;  // cached LaTeX source
    DiagnosticCollector diagnostics_;
    LookupTableStore lookupTableStore_;
    bool lookupTableStorePreloaded_ = false;
};

} // namespace coolsolve
