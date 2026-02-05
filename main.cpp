#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/runner.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"  // For profiling stats
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " [options] <input.eescode>\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -o, --output <file>     Output file (default: stdout)\n";
    std::cerr << "  -f, --format <format>   Output format: json, residuals, latex (default: json)\n";
    std::cerr << "  -c, --compare <file>    Compare with EES .residuals file\n";
    std::cerr << "  -d, --debug [dir]       Debug mode: create output folder with all analysis files\n";
    std::cerr << "                          If dir not specified, creates <input_dir>/<input_name>_coolsolve/\n";
    std::cerr << "  --no-sol                Disable generation of .sol file\n";
    std::cerr << "  -g, --guess             Update .initials file with solution on success\n";
    std::cerr << "  --no-superancillary     Disable CoolProp superancillary functions (faster VLE solving)\n";
    std::cerr << "  -h, --help              Show this help message\n";
}

// Helper to read file content
std::optional<std::string> readFileContent(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) return std::nullopt;
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

int main(int argc, char* argv[]) {
    std::string inputFile;
    std::string outputFile;
    std::string compareFile;
    std::string debugDir;
    std::string format = "json";
    bool debugMode = false;
    bool writeSolFile = true;
    bool enableSuperancillary = true;
    bool updateGuessFile = false;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "--no-sol") {
            writeSolFile = false;
        } else if (arg == "--no-superancillary") {
            enableSuperancillary = false;
        } else if (arg == "-g" || arg == "--guess") {
            updateGuessFile = true;
        } else if (arg == "-d" || arg == "--debug") {
            debugMode = true;
            // Check if next arg is a directory (not starting with - and not ending with .eescode)
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                std::string nextArg = argv[i + 1];
                // If it looks like an input file (.eescode), don't consume it as debug dir
                if (nextArg.size() < 8 || nextArg.substr(nextArg.size() - 8) != ".eescode") {
                    debugDir = argv[++i];
                }
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                outputFile = argv[++i];
            } else {
                std::cerr << "Error: -o requires an argument\n";
                return 1;
            }
        } else if (arg == "-f" || arg == "--format") {
            if (i + 1 < argc) {
                format = argv[++i];
            } else {
                std::cerr << "Error: -f requires an argument\n";
                return 1;
            }
        } else if (arg == "-c" || arg == "--compare") {
            if (i + 1 < argc) {
                compareFile = argv[++i];
            } else {
                std::cerr << "Error: -c requires an argument\n";
                return 1;
            }
        } else if (arg[0] != '-') {
            inputFile = arg;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return 1;
        }
    }
    
    if (inputFile.empty()) {
        std::cerr << "Error: No input file specified\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Apply CoolProp configuration before any CoolProp calls
    coolsolve::CoolPropConfig cpConfig;
    cpConfig.enableSuperancillaries = enableSuperancillary;
    coolsolve::applyCoolPropConfig(cpConfig);
    
    // Use CoolSolveRunner to handle the execution pipeline
    coolsolve::CoolSolveRunner runner(inputFile);
    
    // Configure solver options
    coolsolve::SolverOptions options;
    options.verbose = false; // Disable stdout/stderr verbosity
    
    // Run the pipeline (Parse -> IR -> Infer -> Analyze -> Solve)
    // Note: runner.run() automatically loads .initials if present
    // Pass debugMode as enableTracing
    bool runSuccess = runner.run(options, debugMode);
    
    if (!runner.isParseSuccess()) {
        std::cerr << "Parse failed:\n";
        for (const auto& err : runner.getParseResult().errors) {
            std::cerr << "  Line " << err.line << ": " << err.message << "\n";
        }
        return 1;
    }
    
    // Retrieve results from runner
    const auto& ir = runner.getIR();
    const auto& analysisResult = runner.getAnalysisResult();
    const auto& solveResult = runner.getSolveResult();
    
    // Print statistics
    std::cout << "=== Model Statistics ===\n";
    std::cout << "Equations: " << ir.getEquationCount() << "\n";
    std::cout << "Variables: " << ir.getNonConstantVariableCount() << "\n";
    std::cout << "System square: " << (ir.isSquare() ? "Yes" : "No") << "\n";
    
    if (!ir.isSquare()) {
        auto unmatchedVars = ir.getUnmatchedVariables();
        auto unmatchedEqs = ir.getUnmatchedEquations();
        if (!unmatchedVars.empty()) {
            std::cout << "Unmatched variables (" << unmatchedVars.size() << "): ";
            for (size_t i = 0; i < std::min(unmatchedVars.size(), size_t(10)); ++i) std::cout << unmatchedVars[i] << " ";
            if (unmatchedVars.size() > 10) std::cout << "...";
            std::cout << "\n";
        }
        if (!unmatchedEqs.empty()) {
            std::cout << "Unmatched equations (" << unmatchedEqs.size() << "): ";
            for (size_t i = 0; i < std::min(unmatchedEqs.size(), size_t(10)); ++i) std::cout << unmatchedEqs[i] << " ";
            if (unmatchedEqs.size() > 10) std::cout << "...";
            std::cout << "\n";
        }
    }
    
    // Handle Debug Output generation
    if (debugMode) {
        // Determine debug directory name
        fs::path debugPath;
        if (debugDir.empty()) {
            // Default: <input_basename>_coolsolve/
            fs::path inputPath(inputFile);
            debugPath = inputPath.parent_path() / (inputPath.stem().string() + "_coolsolve");
        } else {
            debugPath = debugDir;
        }
        std::cout << "Debug output: " << fs::weakly_canonical(debugPath) << "\n";
        
        // Read source code for inclusion in debug output
        auto sourceCode = readFileContent(inputFile);
        runner.generateDebugOutput(debugPath.string(), sourceCode.value_or(""));
    }

    // Check for structural analysis errors (e.g., non-square system)
    if (!analysisResult.success) {
        std::cerr << "\n=== Structural Analysis Error ===\n";
        std::cerr << analysisResult.errorMessage << "\n";
        return 1;
    }
    
    // Print structural analysis results
    std::cout << "Total blocks: " << analysisResult.totalBlocks << "\n";
    std::cout << "Largest block: " << analysisResult.largestBlockSize << "\n";
    
    // Report Solver Result
    if (!solveResult.success && solveResult.status != coolsolve::SolverStatus::InvalidInput) {
        std::cout << "\n=== Solver Error ===\n";
        std::cout << "Status: " << coolsolve::statusToString(solveResult.status) << "\n";
        std::cout << "Message: " << solveResult.errorMessage << "\n";
    } else if (solveResult.success) {
        std::cout << "\nSolver: SUCCESS (" << solveResult.totalIterations << " iterations)\n";
    }
    
    // Write .sol file if successful
    if (solveResult.success && writeSolFile) {
        fs::path inputPath(inputFile);
        fs::path solPath = inputPath.parent_path() / (inputPath.stem().string() + ".sol");
        std::ofstream solFile(solPath);
        if (solFile.is_open()) {
            solFile << std::scientific << std::setprecision(12);
            for (const auto& [name, val] : solveResult.variables) {
                solFile << name << " = " << val;
                // Add units if available
                const auto* info = ir.getVariable(name);
                if (info && !info->units.empty()) {
                    solFile << " \"" << info->units << "\"";
                }
                solFile << "\n";
            }
            for (const auto& [name, val] : solveResult.stringVariables) {
                solFile << name << " = '" << val << "'\n";
            }
        } else {
             std::cerr << "Warning: Could not write solution file: " << solPath << "\n";
        }
    }
    
    // Write .initials file if successful and requested
    if (solveResult.success && updateGuessFile) {
        fs::path inputPath(inputFile);
        fs::path initialsPath = inputPath.parent_path() / (inputPath.stem().string() + ".initials");
        std::ofstream initialsFile(initialsPath);
        if (initialsFile.is_open()) {
            initialsFile << std::scientific << std::setprecision(12);
            for (const auto& [name, val] : solveResult.variables) {
                initialsFile << name << " = " << val;
                const auto* info = ir.getVariable(name);
                if (info && !info->units.empty()) {
                    initialsFile << " \"" << info->units << "\"";
                }
                initialsFile << "\n";
            }
            for (const auto& [name, val] : solveResult.stringVariables) {
                initialsFile << name << " = '" << val << "'\n";
            }
        } else {
             std::cerr << "Warning: Could not write initials file: " << initialsPath << "\n";
        }
    }
    
    // Compare with EES reference if requested
    if (!compareFile.empty()) {
        auto comparison = coolsolve::StructuralAnalyzer::compareWithEES(analysisResult, compareFile);
        // Maybe print to stdout if requested, but we don't have verbose mode anymore for stdout.
        // Let's print it always if -c is provided? Or only if fails?
        // User asked to remove verbose output.
        // Assuming comparison result is important if requested.
        std::cout << "EES comparison: " << (comparison.matches ? "PASS" : "FAIL") << "\n";
        if (!comparison.matches) {
            for (const auto& diff : comparison.differences) {
                std::cout << "  - " << diff << "\n";
            }
        }
    }
    
    // Generate output based on format
    std::string output;
    if (format == "json") {
        output = coolsolve::generateAnalysisJSON(ir, analysisResult);
    } else if (format == "residuals") {
        output = coolsolve::generateResidualsReport(ir, analysisResult);
    } else if (format == "latex") {
        output = ir.toLatex();
    } else {
        std::cerr << "Unknown format: " << format << "\n";
        return 1;
    }
    
    // Write standard output
    if (!outputFile.empty()) {
        std::ofstream file(outputFile);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open output file: " << outputFile << "\n";
            return 1;
        }
        file << output;
    } 
    // If no output file specified, do we print to stdout?
    // User complaint was "too much output to stdout".
    // Usually stdout is used for piping. If -o is not specified, piping JSON is standard behavior.
    // But previous code only printed if verbose was on OR output file was empty.
    // "Only print full output to stdout if verbose is on and no output file specified" -> NO, logic was:
    // } else if (verbose) { std::cout << output; }
    // So if NOT verbose and NO output file, it didn't print the huge output?
    // Wait, let's check original code.
    // Line 266: } else if (verbose) { ... print output ... }
    // So previously, by default (no -v, no -o), it printed NOTHING of the JSON/residuals?
    // That seems odd for a CLI tool.
    // But lines 141-161 print Model Statistics to cout.
    // And 187-197 print solver result.
    // So by default it just prints summary.
    // I will preserve this behavior. No verbose -> no huge dump to stdout unless specifically asked via -o (to file) or maybe we should support dumping to stdout if user wants?
    // Usually CLI tools print to stdout if no file.
    // But "options.verbose" controlled it.
    // The user said "-v option has too much output to stdout".
    // So removing -v means we shouldn't have that output.
    // If they want the output, they use -o.
    // Or maybe I should allow dumping if format is specified?
    // Let's stick to: only print summary to stdout. If they want the full JSON/Report, use -o or -d (which creates files).
    
    return solveResult.success ? 0 : 1;
}
