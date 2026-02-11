/**
 * Comprehensive Example Testing
 * 
 * This file recursively tests all .eescode files in a folder.
 * It uses CoolSolveRunner to run the full pipeline for each file.
 * 
 * Run with: ./coolsolve_tests "[examples-comprehensive]"
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include "coolsolve/runner.h"
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

// ============================================================================
// Signal Handling for graceful Ctrl+C
// ============================================================================

static std::atomic<bool> g_interrupted{false};

void sigint_handler(int) {
    g_interrupted = true;
    const char msg[] = "\n[Interrupted by Ctrl+C - finishing current operation...]\n";
#ifdef _WIN32
    std::cout << msg << std::endl;
#else
    ssize_t unused = write(STDOUT_FILENO, msg, sizeof(msg) - 1);
    (void)unused;
#endif
}

void installSignalHandlers() {
    signal(SIGINT, sigint_handler);
    signal(SIGTERM, sigint_handler);
}

bool shouldStop() {
    return g_interrupted;
}

// Default examples directory (relative to build folder)
const std::string DEFAULT_EXAMPLES_DIR = "../examples/";

// Expected solutions for validation (variable name, expected value, tolerance)
struct ExpectedSolution {
    std::string variableName;
    double expectedValue;
    double tolerancePercent = 1.0;  // Default 1% tolerance
};

const std::map<std::string, ExpectedSolution> EXPECTED_SOLUTIONS = {
    {"condenser_3zones.eescode",  {"epsilon_cd_tp", 0.952}},
    {"exchangers1.eescode",       {"Q_dot", 134349}},
    {"exchangers2.eescode",       {"A", 7.378}},
    {"exchangers3.eescode",       {"DELTAT_w", 15.94}},
    {"expander_module.eescode",   {"W_dot_sh_exp", 2202}},
    {"humidair1.eescode",         {"DELTAw", 0.002933}},
    {"humidair2.eescode",         {"M_dot_w_cond", 0.0001046}},
    {"orc_co2.eescode",       {"W_dot_t", 117440}},
    {"orc_complex.eescode",       {"eta_I", 0.108}},
    {"orc_extraction.eescode",    {"eta", 0.1172}},
    {"orc_r245fa.eescode",        {"eta_cycle", 0.1574}},
    {"orc_simple.eescode",        {"eta_cycle", 0.1411}},
    {"orc_solar_complex.eescode", {"eta_overall", 0.05895}},
    {"pressuredrop.eescode",      {"DELTAP", 236849}},
    {"rankine1.eescode",          {"eta", 0.284}},
    {"rankine2.eescode",          {"eta", 0.3033}},
    {"refrigeration1.eescode",    {"COP", 4.818}},
    {"refrigeration2.eescode",    {"COP", 4.472}},
    {"refrigeration3.eescode",    {"COP", 3.481}},
    {"scroll_compressor.eescode", {"epsilon_s_cp", 0.3525}}
};

// Test result structure for reporting
struct ExampleTestResult {
    std::string filename;
    bool parseSuccess = false;
    bool irSuccess = false;
    bool analysisSuccess = false;
    bool solveSuccess = false;
    
    // Solution validation
    bool hasExpectedSolution = false;
    bool solutionValueCorrect = false;
    std::string expectedVarName;
    double expectedValue = 0.0;
    double actualValue = 0.0;
    double percentError = 0.0;
    
    size_t equationCount = 0;
    size_t blockCount = 0;
    
    std::string errorMsg;
    coolsolve::ErrorCategory errorCategory = coolsolve::ErrorCategory::None;
    
    struct BlockStats {
        size_t id;
        bool success;
        int iterations;
        double maxResidual;
        coolsolve::ErrorCategory category;
    };
    std::vector<BlockStats> blockStats;
};

// Check if a file should be skipped (e.g. RTF or unsupported features)
bool shouldSkipFile(const fs::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) return true;
    
    char header[100] = {0};
    file.read(header, 99);
    std::string headerStr(header);
    
    if (headerStr.substr(0, 5) == "{\\rtf") return true;
    
    // Skip files with unsupported advanced features for now
    if (headerStr.find("$bookmark") != std::string::npos ||
        headerStr.find("module") != std::string::npos ||
        headerStr.find("MODULE") != std::string::npos) {
        return true;
    }
    
    return false;
}

std::vector<fs::path> findEescodeFiles(const fs::path& directory) {
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

ExampleTestResult testExampleFile(const fs::path& filepath) {
    ExampleTestResult result;
    result.filename = filepath.filename().string();
    
    if (shouldStop()) return result;
    
    coolsolve::CoolSolveRunner runner(filepath.string());
    coolsolve::SolverOptions options;
    options.tolerance = 1e-6;
    options.timeoutSeconds = 30; // Add 30 second timeout for large models
    
    bool success = runner.run(options);
    
    result.parseSuccess = runner.isParseSuccess();
    if (result.parseSuccess) {
        result.equationCount = runner.getParseResult().equationCount;
        result.irSuccess = runner.isIRSuccess();
        if (result.irSuccess) {
            result.analysisSuccess = runner.isAnalysisSuccess();
            if (result.analysisSuccess) {
                result.blockCount = runner.getAnalysisResult().totalBlocks;
                result.solveSuccess = runner.isSolveSuccess();
                
                // Collect per-block stats
                const auto& solveRes = runner.getSolveResult();
                for (const auto& br : solveRes.blockResults) {
                    ExampleTestResult::BlockStats bs;
                    bs.id = br.id;
                    bs.success = br.success;
                    bs.iterations = br.iterations;
                    bs.maxResidual = br.maxResidual;
                    bs.category = coolsolve::categorizeError(br.errorMessage);
                    result.blockStats.push_back(bs);
                }

                if (!result.solveSuccess) {
                    result.errorMsg = solveRes.errorMessage;
                    result.errorCategory = coolsolve::categorizeError(result.errorMsg);
                } else {
                    // Check solution value if expected solution exists
                    auto it = EXPECTED_SOLUTIONS.find(result.filename);
                    if (it != EXPECTED_SOLUTIONS.end()) {
                        result.hasExpectedSolution = true;
                        result.expectedVarName = it->second.variableName;
                        result.expectedValue = it->second.expectedValue;
                        
                        // Look for the variable in the solution
                        const auto& vars = solveRes.variables;
                        auto varIt = vars.find(result.expectedVarName);
                        if (varIt != vars.end()) {
                            result.actualValue = varIt->second;
                            result.percentError = std::abs(result.actualValue - result.expectedValue) 
                                                  / std::abs(result.expectedValue) * 100.0;
                            result.solutionValueCorrect = (result.percentError <= it->second.tolerancePercent);
                        }
                    }
                }
            } else {
                result.errorMsg = "Structural analysis failed";
                result.errorCategory = coolsolve::ErrorCategory::Other;
            }
        } else {
            result.errorMsg = "IR building failed";
            result.errorCategory = coolsolve::ErrorCategory::Other;
        }
    } else {
        if (!runner.getParseResult().errors.empty()) {
            result.errorMsg = runner.getParseResult().errors[0].message;
            result.errorCategory = coolsolve::categorizeError(result.errorMsg);
        } else {
            result.errorMsg = "Unknown parse error";
            result.errorCategory = coolsolve::ErrorCategory::Other;
        }
    }

    return result;
}

void writeDetailedReport(const fs::path& reportPath, const std::vector<ExampleTestResult>& results) {
    std::ofstream report(reportPath);
    if (!report.is_open()) return;
    
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    
    report << "# CoolSolve Example Files Test Report\n\n";
    report << "**Generated:** " << std::ctime(&time) << "\n";
    
    report << "## Summary\n\n";
    size_t total = results.size();
    size_t parseOk = 0, irOk = 0, analysisOk = 0, solveOk = 0;
    for (const auto& r : results) {
        if (r.parseSuccess) parseOk++;
        if (r.irSuccess) irOk++;
        if (r.analysisSuccess) analysisOk++;
        if (r.solveSuccess) solveOk++;
    }
    
    report << "| Stage | Passed | Total |\n";
    report << "|-------|--------|-------|\n";
    report << "| Parsing | " << parseOk << " | " << total << " |\n";
    report << "| IR Building | " << irOk << " | " << total << " |\n";
    report << "| Analysis | " << analysisOk << " | " << total << " |\n";
    report << "| Solving | " << solveOk << " | " << total << " |\n\n";

    // Error categorization
    std::map<coolsolve::ErrorCategory, size_t> categories;
    size_t totalErrors = 0;
    for (const auto& r : results) {
        if (!r.solveSuccess && r.errorCategory != coolsolve::ErrorCategory::None) {
            categories[r.errorCategory]++;
            totalErrors++;
        }
    }

    if (totalErrors > 0) {
        report << "## Error Categorization\n\n";
        report << "| Category | Count | % |\n";
        report << "|----------|------:|--:|\n";
        for (const auto& [cat, count] : categories) {
            double pct = 100.0 * count / totalErrors;
            report << "| " << coolsolve::categoryToString(cat) << " | " << count 
                   << " | " << std::fixed << std::setprecision(1) << pct << "% |\n";
        }
        report << "\n";
    }
    
    // Solution validation summary
    size_t withExpected = 0, valCorrect = 0;
    for (const auto& r : results) {
        if (r.hasExpectedSolution) {
            withExpected++;
            if (r.solutionValueCorrect) valCorrect++;
        }
    }
    if (withExpected > 0) {
        report << "## Solution Validation\n\n";
        report << "| Metric | Count |\n";
        report << "|--------|------:|\n";
        report << "| Files with expected solution | " << withExpected << " |\n";
        report << "| Correct values (within 1%) | " << valCorrect << " |\n\n";
    }

    report << "## Results by File\n\n";
    report << "| File | Parse | IR | Analysis | Solve | Value Check | Eqs | Blocks |\n";
    report << "|------|-------|----|---------|-------|-------------|----:|-------:|\n";
    for (const auto& r : results) {
        std::string valueCheck = "-";
        if (r.hasExpectedSolution) {
            if (r.solveSuccess) {
                valueCheck = r.solutionValueCorrect ? "OK" : "FAIL";
            } else {
                valueCheck = "N/A";
            }
        }
        report << "| " << r.filename
               << " | " << (r.parseSuccess ? "OK" : "FAIL")
               << " | " << (r.irSuccess ? "OK" : "FAIL")
               << " | " << (r.analysisSuccess ? "OK" : "FAIL")
               << " | " << (r.solveSuccess ? "OK" : "FAIL")
               << " | " << valueCheck
               << " | " << r.equationCount
               << " | " << r.blockCount << " |\n";
    }

    // Detailed errors and value mismatches
    report << "\n## Detailed Errors\n\n";
    for (const auto& r : results) {
        if (!r.solveSuccess) {
            report << "### " << r.filename << "\n\n";
            report << "**Category:** " << coolsolve::categoryToString(r.errorCategory) << "\n\n";
            report << "**Error:** " << r.errorMsg << "\n\n";
            
            if (!r.blockStats.empty()) {
                report << "| Block | Status | Iter | Max Residual | Category |\n";
                report << "|-------|--------|------|--------------|----------|\n";
                for (const auto& bs : r.blockStats) {
                    report << "| " << bs.id 
                           << " | " << (bs.success ? "OK" : "FAIL")
                           << " | " << bs.iterations
                           << " | " << std::scientific << std::setprecision(2) << bs.maxResidual
                           << " | " << coolsolve::categoryToString(bs.category) << " |\n";
                }
                report << "\n";
            }
        } else if (r.hasExpectedSolution && !r.solutionValueCorrect) {
            report << "### " << r.filename << "\n\n";
            report << "**Issue:** Solution value mismatch\n\n";
            report << "| Variable | Expected | Actual | Error |\n";
            report << "|----------|----------|--------|-------|\n";
            report << "| " << r.expectedVarName 
                   << " | " << std::setprecision(6) << r.expectedValue
                   << " | " << std::setprecision(6) << r.actualValue
                   << " | " << std::setprecision(2) << r.percentError << "% |\n\n";
        }
    }
    
    report.close();
}

TEST_CASE("Comprehensive example file testing", "[.][examples-comprehensive]") {
    installSignalHandlers();
    g_interrupted = false;
    
    fs::path examplesDir = fs::path(DEFAULT_EXAMPLES_DIR);
    if (const char* envDir = std::getenv("COOLSOLVE_EXAMPLES_DIR")) {
        examplesDir = fs::path(envDir);
    }
    
    if (!fs::exists(examplesDir)) {
        SKIP("Examples directory not found: " << examplesDir.string());
    }
    
    auto files = findEescodeFiles(examplesDir);
    if (files.empty()) {
        SKIP("No .eescode files found in: " << examplesDir.string());
    }
    
    std::cout << "\nTesting " << files.size() << " example files using CoolSolveRunner...\n";
    std::cout << "(Press Ctrl+C to interrupt gracefully)\n\n";
    
    std::vector<ExampleTestResult> results;
    for (const auto& file : files) {
        if (shouldStop()) break;
        
        std::cout << "  " << file.filename().string() << "... " << std::flush;
        auto result = testExampleFile(file);
        results.push_back(result);
        
        if (result.solveSuccess) {
            std::cout << "OK";
            if (result.hasExpectedSolution) {
                if (result.solutionValueCorrect) {
                    std::cout << " [" << result.expectedVarName << "=" 
                              << std::setprecision(4) << result.actualValue << " OK]";
                } else {
                    std::cout << " [" << result.expectedVarName << "=" 
                              << std::setprecision(4) << result.actualValue 
                              << " EXPECTED " << result.expectedValue 
                              << " ERROR " << std::setprecision(2) << result.percentError << "%]";
                }
            }
        }
        else if (result.analysisSuccess) std::cout << "SOLVE FAIL";
        else if (result.irSuccess) std::cout << "ANALYSIS FAIL";
        else if (result.parseSuccess) std::cout << "IR FAIL";
        else std::cout << "PARSE FAIL";
        std::cout << " (" << result.blockCount << " blocks)\n";
    }
    
    fs::path reportPath = fs::absolute(examplesDir) / "test_examples.md";
    writeDetailedReport(reportPath, results);
    std::cout << "\nDetailed report written to: " << reportPath.string() << "\n";
}

#include <catch2/catch_session.hpp>

int main(int argc, char* argv[]) {
    Catch::Session session;
    std::vector<char*> new_argv;
    for (int i = 0; i < argc; ++i) new_argv.push_back(argv[i]);
    if (argc == 1) {
        static char run_all[] = "~[.], [examples-comprehensive]";
        new_argv.push_back(run_all);
    }
    return session.run(static_cast<int>(new_argv.size()), new_argv.data());
}
