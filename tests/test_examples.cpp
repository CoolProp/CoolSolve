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
#include "coolsolve/solution_checker.h"
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
    // Original models
    {"condenser_3zones.eescode",  {"epsilon_cd_tp", 0.952}},
    {"cooling_coil.eescode",      {"epsilon_c", 0.8111}},
    {"cpbar.eescode",             {"c_bar_p", 1020.3}},
    {"exchangers1.eescode",       {"Q_dot", 134349}},
    {"exchangers2.eescode",       {"A", 7.378}},
    {"exchangers3.eescode",       {"DELTAT_w", 15.94}},
    {"expander_module.eescode",   {"W_dot_sh_exp", 2202}},
    {"humidair1.eescode",         {"DELTAw", 0.002933}},
    {"humidair2.eescode",         {"M_dot_w_cond", 0.0001046}},
    {"orc_co2.eescode",           {"W_dot_t", 117440}},
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
    {"scroll_compressor.eescode", {"epsilon_s_cp", 0.2424}},
    // Models from notsolving/ (fixed)
    {"internal_combustion_engine.eescode",      {"W_dot", 190817}},
    {"internal_combustion_engine_cpbar.eescode", {"M_dot_a", 0.004209}},
    {"piston_compressor.eescode",               {"C", 0.04694}},
    {"turbocompressor.eescode",                 {"M_dot", 0.6660}},
    // Models from solving/
    {"air_screw_compressor.eescode",            {"epsilon_s", 0.4896}},
    {"air_screw_compressor_simple.eescode",     {"epsilon_s", 0.5918}},
    {"boiler_cpbar.eescode",                    {"eta", 0.8853}},
    {"boiler_cpbar2.eescode",                   {"eta", 0.8838}},
    {"compressor_refrigeration_simple.eescode", {"COP_1", 3.934}},
    {"condenser_wet.eescode",                   {"AU_f", 1767.5}},
    {"cooling_tower.eescode",                   {"Q_dot", 1977305}},
    {"cooling_tower2.eescode",                  {"Q_dot", 1977305}},
    {"evaporator.eescode",                      {"epsilon_ev_f", 0.7153}},
    {"heat_pump_MSTh_SB_R10.eescode",           {"COP", 3.311}},
    {"refrigeration_compressor.eescode",        {"epsilon_v_1", 0.4729}},
    {"simple_centrifugal_compressor.eescode",   {"U", 354.8}},
    {"zorlu_heat_pump.eescode",                 {"COP_HP", 6.4167}},
    // Advanced features showcase (exercises all CoolSolve features)
    {"advanced_features.eescode",               {"verification", 2083.5}},
    // Lookup table demo (lookup_demo.csv + lookup_demo_watercp.csv loaded as companion tables)
    {"lookup_demo.eescode",                     {"Cp_water", 4.2675}},
};

// Test result structure for reporting
struct ExampleTestResult {
    std::string filename;
    bool parseSuccess = false;
    bool irSuccess = false;
    bool analysisSuccess = false;
    bool solveSuccess = false;
    
    /** Total pipeline time (ms) from runner's existing timing (parse+IR+infer+analysis+solve). */
    double totalTimeMs = 0.0;
    
    // Solution validation
    bool hasExpectedSolution = false;
    bool solutionValueCorrect = false;
    std::string expectedVarName;
    double expectedValue = 0.0;
    double actualValue = 0.0;
    double percentError = 0.0;
    
    // Full equation verification
    bool equationCheckDone = false;
    bool allEquationsSatisfied = false;
    size_t equationsChecked = 0;
    size_t equationsViolated = 0;
    double maxResidual = 0.0;
    double maxRelativeError = 0.0;
    std::string worstEquationText;
    
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

    // Load coolsolve.conf from the same directory as the model (same as CLI)
    auto configPath = filepath.parent_path() / "coolsolve.conf";
    if (fs::exists(configPath)) {
        coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
        options.tolerance = 1e-6;          // keep our test tolerance
        options.timeoutSeconds = 30;       // keep our test timeout
    }
    
    bool success = runner.run(options);
    
    const auto& timing = runner.getTiming();
    result.totalTimeMs = timing.total_time_ms;
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
                    
                    // Full equation verification: check all equations are satisfied
                    auto checkRes = coolsolve::checkSolution(
                        runner.getIR(), solveRes.variables, solveRes.stringVariables,
                        options.coolpropConfig, 1e-3, &runner.getLookupTableStore());
                    result.equationCheckDone = true;
                    result.allEquationsSatisfied = checkRes.allSatisfied;
                    result.equationsChecked = checkRes.satisfiedCount + checkRes.violatedCount;
                    result.equationsViolated = checkRes.violatedCount;
                    result.maxResidual = checkRes.maxResidual;
                    result.maxRelativeError = checkRes.maxRelativeError;
                    result.worstEquationText = checkRes.worstEquationText;
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
    double totalTimeMs = 0.0;
    for (const auto& r : results) {
        if (r.parseSuccess) parseOk++;
        if (r.irSuccess) irOk++;
        if (r.analysisSuccess) analysisOk++;
        if (r.solveSuccess) solveOk++;
        totalTimeMs += r.totalTimeMs;
    }
    
    report << "| Stage | Passed | Total |\n";
    report << "|-------|--------|-------|\n";
    report << "| Parsing | " << parseOk << " | " << total << " |\n";
    report << "| IR Building | " << irOk << " | " << total << " |\n";
    report << "| Analysis | " << analysisOk << " | " << total << " |\n";
    report << "| Solving | " << solveOk << " | " << total << " |\n";
    report << "| **Total time** | **" << std::fixed << std::setprecision(2) << (totalTimeMs / 1000.0) << " s** | (all models) |\n\n";

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
    report << "| File | Parse | IR | Analysis | Solve | Value Check | Eq Check | Eqs | Blocks | Time (s) |\n";
    report << "|------|-------|----|---------|-------|-------------|----------|----:|-------:|--------:|\n";
    for (const auto& r : results) {
        std::string valueCheck = "-";
        if (r.hasExpectedSolution) {
            if (r.solveSuccess) {
                valueCheck = r.solutionValueCorrect ? "OK" : "FAIL";
            } else {
                valueCheck = "N/A";
            }
        }
        std::string eqCheck = "-";
        if (r.equationCheckDone) {
            if (r.allEquationsSatisfied) {
                eqCheck = "ALL OK";
            } else {
                eqCheck = std::to_string(r.equationsViolated) + " FAIL";
            }
        }
        report << "| " << r.filename
               << " | " << (r.parseSuccess ? "OK" : "FAIL")
               << " | " << (r.irSuccess ? "OK" : "FAIL")
               << " | " << (r.analysisSuccess ? "OK" : "FAIL")
               << " | " << (r.solveSuccess ? "OK" : "FAIL")
               << " | " << valueCheck
               << " | " << eqCheck
               << " | " << r.equationCount
               << " | " << r.blockCount
               << " | " << std::fixed << std::setprecision(2) << (r.totalTimeMs / 1000.0) << " |\n";
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
        // Report equation check failures
        if (r.equationCheckDone && !r.allEquationsSatisfied) {
            report << "### " << r.filename << " - Equation Check\n\n";
            report << "**" << r.equationsViolated << " equation(s) violated** out of "
                   << r.equationsChecked << " checked\n\n";
            report << "- Max |residual|: " << std::scientific << std::setprecision(4) << r.maxResidual << "\n";
            report << "- Max relative error: " << r.maxRelativeError << "\n";
            report << "- Worst equation: " << r.worstEquationText << "\n\n";
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
            std::cout << " (" << std::fixed << std::setprecision(2) << (result.totalTimeMs / 1000.0) << "s)";
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
            if (result.equationCheckDone) {
                if (result.allEquationsSatisfied) {
                    std::cout << " [EqCheck: ALL OK]";
                } else {
                    std::cout << " [EqCheck: " << result.equationsViolated << " VIOLATED"
                              << " maxRel=" << std::scientific << std::setprecision(2)
                              << result.maxRelativeError << "]";
                }
            }
        }
        else if (result.analysisSuccess) {
            auto cat = coolsolve::categorizeError(result.errorMsg);
            if (cat == coolsolve::ErrorCategory::UnsupportedFunction || 
                result.errorMsg.find("Unknown fluid") != std::string::npos)
                std::cout << "UNSUPPORTED";
            else
                std::cout << "SOLVE FAIL";
        }
        else if (result.irSuccess) std::cout << "ANALYSIS FAIL";
        else if (result.parseSuccess) std::cout << "IR FAIL";
        else std::cout << "PARSE FAIL";
        std::cout << " (" << result.blockCount << " blocks";
        if (!result.solveSuccess && result.totalTimeMs > 0.0) std::cout << ", " << std::fixed << std::setprecision(2) << (result.totalTimeMs / 1000.0) << "s";
        std::cout << ")\n";
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
