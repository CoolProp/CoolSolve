/**
 * Solver Robustness Diagnosis
 *
 * Tests all example .eescode files with each individual solver strategy and
 * also without initial values, collecting results into a Markdown report.
 *
 * Run with: ./coolsolve_tests "[solver-robustness]"
 *
 * The report is written to examples/solver_robustness_report.md
 */

#include <catch2/catch_test_macros.hpp>
#include "coolsolve/runner.h"
#include "coolsolve/variable_inference.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <csignal>
#include <atomic>
#include <chrono>

namespace fs = std::filesystem;

static std::atomic<bool> g_rob_interrupted{false};

static void robustness_sigint_handler(int) {
    g_rob_interrupted = true;
}

static bool robShouldStop() { return g_rob_interrupted; }

// Default examples directory (relative to build folder)
static const std::string ROBUSTNESS_EXAMPLES_DIR = "../examples/";

struct RobustnessResult {
    std::string filename;
    bool parseOk = false;
    bool irOk = false;
    bool analysisOk = false;
    bool solveOk = false;
    double totalTimeMs = 0.0;
    int totalIterations = 0;
    std::string errorMsg;
    coolsolve::ErrorCategory errorCategory = coolsolve::ErrorCategory::None;

    // Additional diagnostics for better reporting
    int failedBlockIndex = -1;
    int failedBlockSize = -1;
    double initialResidualNorm = -1.0;
    double finalResidualNorm = -1.0;

    /// Short status string for compact table display
    std::string compactStatus() const {
        if (!parseOk) return "PARSE";
        if (!analysisOk) return "ANALYSIS";
        if (solveOk) {
            std::ostringstream ss;
            ss << "OK " << std::fixed << std::setprecision(2) << (totalTimeMs / 1000.0) << "s";
            return ss.str();
        }
        // Build a compact failure tag
        std::string cat = coolsolve::categoryToString(errorCategory);
        std::ostringstream ss;
        ss << cat;
        if (failedBlockSize > 0) ss << " blk" << failedBlockSize;
        if (finalResidualNorm > 0) {
            ss << " |F|=";
            if (finalResidualNorm >= 1e6) ss << std::scientific << std::setprecision(0) << finalResidualNorm;
            else ss << std::fixed << std::setprecision(1) << finalResidualNorm;
        }
        return ss.str();
    }
};

struct SolverConfig {
    std::string label;
    std::vector<coolsolve::SolverStrategy> pipeline;
    bool useInitials;
    bool enableTearing;
};

static bool shouldSkipFile(const fs::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) return true;
    char header[100] = {0};
    file.read(header, 99);
    std::string headerStr(header);
    if (headerStr.substr(0, 5) == "{\\rtf") return true;
    if (headerStr.find("$bookmark") != std::string::npos ||
        headerStr.find("module") != std::string::npos ||
        headerStr.find("MODULE") != std::string::npos) {
        return true;
    }
    return false;
}

static std::vector<fs::path> findFiles(const fs::path& directory) {
    std::vector<fs::path> files;
    if (!fs::exists(directory)) return files;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".eescode") {
            if (!shouldSkipFile(entry.path())) {
                files.push_back(entry.path());
            }
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

static RobustnessResult testFile(const fs::path& filepath, const SolverConfig& cfg) {
    RobustnessResult res;
    res.filename = filepath.filename().string();

    coolsolve::CoolSolveRunner runner(filepath.string());
    coolsolve::SolverOptions opts;
    opts.tolerance = 1e-6;
    opts.timeoutSeconds = 30;
    opts.solverPipeline = cfg.pipeline;
    opts.pipelineMode = coolsolve::SolverPipelineMode::Sequential;
    opts.enableTearing = cfg.enableTearing;

    // If we want to skip initials, we need to run the pipeline stages manually
    // to avoid auto-loading from .initials file
    if (!cfg.useInitials) {
        // We replicate the runner pipeline but skip initials loading

        // 1. Parse
        coolsolve::EESParser parser;
        auto parseResult = parser.parseFile(filepath.string());
        res.parseOk = parseResult.success;
        if (!res.parseOk) { res.errorMsg = "Parse failed"; return res; }

        // 2. Build IR
        std::unique_ptr<coolsolve::IR> ir;
        try {
            ir = std::make_unique<coolsolve::IR>(coolsolve::IR::fromAST(parseResult.program));
            res.irOk = true;
        } catch (...) { res.errorMsg = "IR failed"; return res; }

        // 3. Infer & initialize
        try {
            coolsolve::inferVariables(*ir);
            coolsolve::initializeVariables(*ir);
        } catch (...) {}

        // 4. Structural analysis
        coolsolve::StructuralAnalysisResult analysis;
        try {
            analysis = coolsolve::StructuralAnalyzer::analyze(*ir);
            res.analysisOk = analysis.success;
            if (!analysis.success) { res.errorMsg = "Analysis failed"; return res; }
        } catch (...) { res.errorMsg = "Analysis exception"; return res; }

        // 5. Solve (NO initials loaded)
        try {
            auto t1 = std::chrono::high_resolution_clock::now();
            coolsolve::Solver solver(*ir, analysis);
            auto solveResult = solver.solve(opts, false);
            auto t2 = std::chrono::high_resolution_clock::now();
            res.totalTimeMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
            res.solveOk = solveResult.success;
            res.totalIterations = solveResult.totalIterations;
            if (!solveResult.success) {
                res.errorMsg = solveResult.errorMessage;
                res.errorCategory = coolsolve::categorizeError(res.errorMsg);
                // Find first failed block for diagnostics
                for (const auto& br : solveResult.blockResults) {
                    if (!br.success) {
                        res.failedBlockIndex = static_cast<int>(br.id);
                        res.failedBlockSize = 0; // populated below from analysis
                        res.finalResidualNorm = br.maxResidual;
                        break;
                    }
                }
                // Get block size from analysis
                if (res.failedBlockIndex >= 0) {
                    for (const auto& blk : analysis.blocks) {
                        if (static_cast<int>(blk.id) == res.failedBlockIndex) {
                            res.failedBlockSize = static_cast<int>(blk.variables.size());
                            break;
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            res.errorMsg = std::string("Exception: ") + e.what();
        }
    } else {
        // Normal run with initials
        bool ok = runner.run(opts, false);
        res.totalTimeMs = runner.getTiming().total_time_ms;
        res.parseOk = runner.isParseSuccess();
        res.irOk = runner.isIRSuccess();
        res.analysisOk = runner.isAnalysisSuccess();
        res.solveOk = ok;
        if (runner.isAnalysisSuccess()) {
            res.totalIterations = runner.getSolveResult().totalIterations;
            if (!ok) {
                res.errorMsg = runner.getSolveResult().errorMessage;
                res.errorCategory = coolsolve::categorizeError(res.errorMsg);
                // Find first failed block for diagnostics
                for (const auto& br : runner.getSolveResult().blockResults) {
                    if (!br.success) {
                        res.failedBlockIndex = static_cast<int>(br.id);
                        res.finalResidualNorm = br.maxResidual;
                        break;
                    }
                }
                if (res.failedBlockIndex >= 0) {
                    for (const auto& blk : runner.getAnalysisResult().blocks) {
                        if (static_cast<int>(blk.id) == res.failedBlockIndex) {
                            res.failedBlockSize = static_cast<int>(blk.variables.size());
                            break;
                        }
                    }
                }
            }
        } else if (runner.isParseSuccess()) {
            res.errorMsg = "IR or Analysis failed";
        } else {
            res.errorMsg = "Parse failed";
        }
    }

    return res;
}

TEST_CASE("Solver robustness diagnosis", "[.][solver-robustness]") {
    signal(SIGINT, robustness_sigint_handler);
    signal(SIGTERM, robustness_sigint_handler);
    g_rob_interrupted = false;

    fs::path exDir = fs::path(ROBUSTNESS_EXAMPLES_DIR);
    if (const char* envDir = std::getenv("COOLSOLVE_EXAMPLES_DIR"))
        exDir = fs::path(envDir);

    REQUIRE(fs::exists(exDir));
    auto files = findFiles(exDir);
    REQUIRE(!files.empty());

    // Define solver configurations to test
    std::vector<SolverConfig> configs = {
        // -- With initials --
        {"Default pipeline (with initials)",
         {coolsolve::SolverStrategy::Newton,
          coolsolve::SolverStrategy::TrustRegion,
          coolsolve::SolverStrategy::LevenbergMarquardt,
          coolsolve::SolverStrategy::Homotopy,
          coolsolve::SolverStrategy::Partitioned},
         true, false},

        {"Newton only (with initials)",
         {coolsolve::SolverStrategy::Newton},
         true, false},

        {"TrustRegion only (with initials)",
         {coolsolve::SolverStrategy::TrustRegion},
         true, false},

        {"LevenbergMarquardt only (with initials)",
         {coolsolve::SolverStrategy::LevenbergMarquardt},
         true, false},

        {"Homotopy only (with initials)",
         {coolsolve::SolverStrategy::Homotopy},
         true, false},

        {"Partitioned only (with initials)",
         {coolsolve::SolverStrategy::Partitioned},
         true, false},

        {"Default + Tearing (with initials)",
         {coolsolve::SolverStrategy::Newton,
          coolsolve::SolverStrategy::TrustRegion,
          coolsolve::SolverStrategy::LevenbergMarquardt,
          coolsolve::SolverStrategy::Homotopy,
          coolsolve::SolverStrategy::Partitioned},
         true, true},

        // -- Without initials --
        {"Default pipeline (NO initials)",
         {coolsolve::SolverStrategy::Newton,
          coolsolve::SolverStrategy::TrustRegion,
          coolsolve::SolverStrategy::LevenbergMarquardt,
          coolsolve::SolverStrategy::Homotopy,
          coolsolve::SolverStrategy::Partitioned},
         false, false},

        {"Newton only (NO initials)",
         {coolsolve::SolverStrategy::Newton},
         false, false},

        {"TrustRegion only (NO initials)",
         {coolsolve::SolverStrategy::TrustRegion},
         false, false},

        {"LevenbergMarquardt only (NO initials)",
         {coolsolve::SolverStrategy::LevenbergMarquardt},
         false, false},

        {"Homotopy only (NO initials)",
         {coolsolve::SolverStrategy::Homotopy},
         false, false},

        {"Partitioned only (NO initials)",
         {coolsolve::SolverStrategy::Partitioned},
         false, false},

        {"Default + Tearing (NO initials)",
         {coolsolve::SolverStrategy::Newton,
          coolsolve::SolverStrategy::TrustRegion,
          coolsolve::SolverStrategy::LevenbergMarquardt,
          coolsolve::SolverStrategy::Homotopy,
          coolsolve::SolverStrategy::Partitioned},
         false, true},
    };

    // Collect results: configs x files
    // results[configIdx][fileIdx]
    std::vector<std::vector<RobustnessResult>> allResults(configs.size());

    for (size_t ci = 0; ci < configs.size(); ++ci) {
        if (robShouldStop()) break;
        const auto& cfg = configs[ci];
        std::cout << "\n=== " << cfg.label << " ===\n";
        allResults[ci].reserve(files.size());

        for (const auto& file : files) {
            if (robShouldStop()) break;
            std::cout << "  " << file.filename().string() << "... " << std::flush;
            auto r = testFile(file, cfg);
            allResults[ci].push_back(r);

            if (!r.parseOk)          std::cout << "PARSE FAIL";
            else if (!r.analysisOk)  std::cout << "ANALYSIS FAIL";
            else if (r.solveOk)      std::cout << "OK (" << std::fixed << std::setprecision(2) << (r.totalTimeMs/1000.0) << "s)";
            else                     std::cout << "FAIL";
            std::cout << "\n";
        }
    }

    // ========================================================================
    // Write combined report
    // ========================================================================
    fs::path reportPath = fs::absolute(exDir) / "solver_robustness_report.md";
    std::ofstream report(reportPath);
    REQUIRE(report.is_open());

    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);

    report << "# CoolSolve Solver Robustness Report\n\n";
    report << "**Generated:** " << std::ctime(&time) << "\n";
    report << "Total example files tested: " << files.size() << "\n\n";
    report << "**Legend:** OK = converged, PARSE = parse error, ANALYSIS = structural analysis error.  \n";
    report << "Failure cells show: `ErrorCategory blkN |F|=residual` where N is the failed block size.\n\n";

    // ---- Summary table ----
    report << "## Summary: Solve Success Rate by Configuration\n\n";
    report << "| # | Configuration | Initials | Tearing | Solved | Total | Rate | Avg time (s) |\n";
    report << "|---:|---|:---:|:---:|---:|---:|---:|---:|\n";
    for (size_t ci = 0; ci < configs.size(); ++ci) {
        int solved = 0, total = 0;
        double totalTime = 0;
        for (const auto& r : allResults[ci]) {
            if (r.analysisOk) {
                total++;
                if (r.solveOk) { solved++; totalTime += r.totalTimeMs; }
            }
        }
        double avgTime = solved > 0 ? (totalTime / solved / 1000.0) : 0;
        double rate = total > 0 ? (100.0 * solved / total) : 0;
        report << "| " << (ci + 1)
               << " | " << configs[ci].label
               << " | " << (configs[ci].useInitials ? "Yes" : "No")
               << " | " << (configs[ci].enableTearing ? "Yes" : "No")
               << " | " << solved
               << " | " << total
               << " | " << std::fixed << std::setprecision(1) << rate << "%"
               << " | " << std::setprecision(3) << avgTime
               << " |\n";
    }

    // ---- Helper lambda for compact cell content ----
    auto cellContent = [](const RobustnessResult& r) -> std::string {
        if (!r.parseOk) return "PARSE";
        if (!r.analysisOk) return "ANALYSIS";
        if (r.solveOk) {
            std::ostringstream ss;
            ss << "**OK** (" << std::fixed << std::setprecision(2) << (r.totalTimeMs/1000.0) << "s)";
            return ss.str();
        }
        // Failure: show compact diagnostic
        return r.compactStatus();
    };

    // ---- Per-file detail table (with initials) ----
    // Use short column labels to keep the table readable
    std::vector<std::string> withInitLabels;
    std::vector<size_t> withInitIdx;
    for (size_t ci = 0; ci < configs.size(); ++ci) {
        if (configs[ci].useInitials) {
            withInitIdx.push_back(ci);
            // Build short label: solver name(s) + tearing flag
            std::string shortLabel;
            for (auto s : configs[ci].pipeline) {
                if (!shortLabel.empty()) shortLabel += "+";
                std::string sn = coolsolve::strategyToString(s);
                if (sn == "LevenbergMarquardt") sn = "LM";
                else if (sn == "Partitioned") sn = "Part";
                else if (sn == "TrustRegion") sn = "TR";
                else if (sn == "Newton") sn = "Nwt";
                shortLabel += sn;
            }
            if (configs[ci].enableTearing) shortLabel += "+Tear";
            withInitLabels.push_back(shortLabel);
        }
    }

    report << "\n## Detailed Results: With Initials\n\n";
    report << "| File |";
    for (const auto& l : withInitLabels) report << " " << l << " |";
    report << "\n|---|";
    for (size_t i = 0; i < withInitLabels.size(); ++i) report << "---|";
    report << "\n";

    for (size_t fi = 0; fi < files.size(); ++fi) {
        report << "| " << files[fi].filename().string() << " |";
        for (size_t ci : withInitIdx) {
            if (fi < allResults[ci].size())
                report << " " << cellContent(allResults[ci][fi]) << " |";
            else
                report << " — |";
        }
        report << "\n";
    }

    // ---- Per-file detail table (without initials) ----
    std::vector<std::string> noInitLabels;
    std::vector<size_t> noInitIdx;
    for (size_t ci = 0; ci < configs.size(); ++ci) {
        if (!configs[ci].useInitials) {
            noInitIdx.push_back(ci);
            std::string shortLabel;
            for (auto s : configs[ci].pipeline) {
                if (!shortLabel.empty()) shortLabel += "+";
                std::string sn = coolsolve::strategyToString(s);
                if (sn == "LevenbergMarquardt") sn = "LM";
                else if (sn == "Partitioned") sn = "Part";
                else if (sn == "TrustRegion") sn = "TR";
                else if (sn == "Newton") sn = "Nwt";
                shortLabel += sn;
            }
            if (configs[ci].enableTearing) shortLabel += "+Tear";
            noInitLabels.push_back(shortLabel);
        }
    }

    report << "\n## Detailed Results: Without Initials\n\n";
    report << "| File |";
    for (const auto& l : noInitLabels) report << " " << l << " |";
    report << "\n|---|";
    for (size_t i = 0; i < noInitLabels.size(); ++i) report << "---|";
    report << "\n";

    for (size_t fi = 0; fi < files.size(); ++fi) {
        report << "| " << files[fi].filename().string() << " |";
        for (size_t ci : noInitIdx) {
            if (fi < allResults[ci].size())
                report << " " << cellContent(allResults[ci][fi]) << " |";
            else
                report << " — |";
        }
        report << "\n";
    }

    // ---- Hardest models: which models fail most often? ----
    report << "\n## Model Difficulty Ranking\n\n";
    report << "Models ranked by number of configurations that failed to solve them.\n\n";
    report << "| File | Failures / Configs | Failed Configurations |\n";
    report << "|---|---:|---|\n";

    struct ModelDifficulty { std::string name; int failures; int total; std::string failedCfgs; };
    std::vector<ModelDifficulty> difficulty;
    for (size_t fi = 0; fi < files.size(); ++fi) {
        ModelDifficulty md;
        md.name = files[fi].filename().string();
        md.failures = 0;
        md.total = 0;
        for (size_t ci = 0; ci < configs.size(); ++ci) {
            if (fi >= allResults[ci].size()) continue;
            const auto& r = allResults[ci][fi];
            if (!r.analysisOk) continue; // skip parse/analysis issues
            md.total++;
            if (!r.solveOk) {
                md.failures++;
                if (!md.failedCfgs.empty()) md.failedCfgs += ", ";
                md.failedCfgs += configs[ci].label;
            }
        }
        difficulty.push_back(md);
    }
    std::sort(difficulty.begin(), difficulty.end(),
        [](const ModelDifficulty& a, const ModelDifficulty& b) { return a.failures > b.failures; });
    for (const auto& md : difficulty) {
        if (md.failures == 0 && md.total > 0) continue; // skip fully-solved models
        if (md.total == 0) continue; // skip un-parseable models
        report << "| " << md.name << " | " << md.failures << " / " << md.total << " | "
               << (md.failedCfgs.empty() ? "—" : md.failedCfgs) << " |\n";
    }

    // ---- Error category summary ----
    report << "\n## Error Category Breakdown\n\n";
    report << "Across all configurations and models:\n\n";
    std::map<std::string, int> categoryCount;
    int totalFailures = 0;
    for (size_t ci = 0; ci < configs.size(); ++ci) {
        for (const auto& r : allResults[ci]) {
            if (r.analysisOk && !r.solveOk) {
                categoryCount[coolsolve::categoryToString(r.errorCategory)]++;
                totalFailures++;
            }
        }
    }
    report << "| Error Category | Count | Fraction |\n";
    report << "|---|---:|---:|\n";
    for (const auto& [cat, count] : categoryCount) {
        report << "| " << cat << " | " << count << " | "
               << std::fixed << std::setprecision(1) << (100.0 * count / std::max(1, totalFailures))
               << "% |\n";
    }

    // ---- Detailed error messages for failures ----
    report << "\n## Detailed Error Messages\n\n";
    for (size_t ci = 0; ci < configs.size(); ++ci) {
        bool anyFail = false;
        for (const auto& r : allResults[ci]) {
            if (r.analysisOk && !r.solveOk) { anyFail = true; break; }
        }
        if (!anyFail) continue;

        report << "### " << configs[ci].label << "\n\n";
        report << "| File | Category | Block | Residual | Error (truncated) |\n";
        report << "|---|---|---:|---:|---|\n";
        for (const auto& r : allResults[ci]) {
            if (r.analysisOk && !r.solveOk) {
                std::string errTrunc = r.errorMsg.substr(0, 100);
                for (auto& c : errTrunc) { if (c == '|') c = '/'; if (c == '\n') c = ' '; }
                report << "| " << r.filename
                       << " | " << coolsolve::categoryToString(r.errorCategory)
                       << " | " << (r.failedBlockSize > 0 ? std::to_string(r.failedBlockSize) : "?")
                       << " | ";
                if (r.finalResidualNorm > 0) {
                    if (r.finalResidualNorm >= 1e6) report << std::scientific << std::setprecision(1) << r.finalResidualNorm;
                    else report << std::fixed << std::setprecision(2) << r.finalResidualNorm;
                } else {
                    report << "—";
                }
                report << " | " << errTrunc << " |\n";
            }
        }
        report << "\n";
    }

    report.close();
    std::cout << "\nRobustness report written to: " << reportPath.string() << "\n";
}
