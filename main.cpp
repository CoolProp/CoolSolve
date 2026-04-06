#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/runner.h"
#include "coolsolve/diagnostic.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"  // For profiling stats
#include "coolsolve/solution_checker.h"
#include "coolsolve/latex_report.h"
#ifdef COOLSOLVE_GUI
#include "coolsolve/server.h"
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <iomanip>
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#endif

namespace fs = std::filesystem;

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " [options] <input.eescode>\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -o, --output <file>     Output file (default: stdout)\n";
    std::cerr << "  -f, --format <format>   Output format: json, latex (default: json)\n";
    std::cerr << "  -d, --debug [dir]       Debug mode: create output folder with all analysis files\n";
    std::cerr << "                          If dir not specified, creates <input_dir>/<input_name>_coolsolve/\n";
    std::cerr << "  --no-sol                Disable generation of .sol file\n";
    std::cerr << "  -g, --guess             Update .initials file with solution on success\n";
    std::cerr << "  --no-superancillary     Disable CoolProp superancillary functions (faster VLE solving)\n";
#ifdef COOLSOLVE_GUI
    std::cerr << "  --gui [port]            Start the GUI web server (default port: 8550)\n";
    std::cerr << "  --no-browser            Don't open browser automatically (with --gui)\n";
#endif
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
#ifdef _WIN32
    // Enable UTF-8 output so box-drawing characters render correctly in any terminal
    SetConsoleOutputCP(CP_UTF8);
#endif
    std::string inputFile;
    std::string outputFile;
    std::string debugDir;
    std::string format = "json";
    bool debugMode = false;
    bool writeSolFile = true;
    bool enableSuperancillary = true;
    bool updateGuessFile = false;
#ifdef COOLSOLVE_GUI
    bool guiMode = false;
    int guiPort = 8550;
    bool guiOpenBrowser = true;
#endif
    
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
#ifdef COOLSOLVE_GUI
        } else if (arg == "--gui") {
            guiMode = true;
            // Check if next arg is a port number
            if (i + 1 < argc) {
                try {
                    int port = std::stoi(argv[i + 1]);
                    if (port > 0 && port < 65536) {
                        guiPort = port;
                        ++i;
                    }
                } catch (...) {}
            }
        } else if (arg == "--no-browser") {
            guiOpenBrowser = false;
#endif
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
        } else if (arg[0] != '-') {
            inputFile = arg;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return 1;
        }
    }
    
    // Resolve binary directory for resource discovery (examples, etc.)
    fs::path binaryDir;
#ifdef _WIN32
    {
        wchar_t buf[MAX_PATH];
        GetModuleFileNameW(nullptr, buf, MAX_PATH);
        binaryDir = fs::path(buf).parent_path();
    }
#else
    try { binaryDir = fs::canonical(argv[0]).parent_path(); } catch (...) {}
#endif

    if (inputFile.empty()) {
#ifdef COOLSOLVE_GUI
        // No input file: default to GUI mode (e.g. launched from Start Menu with no args)
        coolsolve::ServerOptions serverOpts;
        serverOpts.port = guiPort;
        serverOpts.openBrowser = guiOpenBrowser;
        serverOpts.examplesDir = (binaryDir / "examples").string();
        return coolsolve::startServer(serverOpts);
#else
        std::cerr << "Error: No input file specified\n";
        printUsage(argv[0]);
        return 1;
#endif
    }

#ifdef COOLSOLVE_GUI
    // GUI mode with input file: open it in the GUI
    if (guiMode) {
        coolsolve::ServerOptions serverOpts;
        serverOpts.port = guiPort;
        serverOpts.openBrowser = guiOpenBrowser;
        serverOpts.initialFile = fs::absolute(inputFile).string();
        serverOpts.examplesDir = (binaryDir / "examples").string();
        return coolsolve::startServer(serverOpts);
    }
#endif
    
    // Configure solver options (defaults from solver.h)
    coolsolve::SolverOptions options;
    options.verbose = false; // Disable stdout/stderr verbosity
    options.coolpropConfig.enableSuperancillaries = enableSuperancillary;

    // Override from coolsolve.conf in the same folder as the input file (not subfolders)
    fs::path configPath = fs::path(inputFile).parent_path() / "coolsolve.conf";
    if (fs::exists(configPath)) {
        coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    }
    
    // Apply CoolProp configuration before any CoolProp calls
    coolsolve::applyCoolPropConfig(options.coolpropConfig);
    
    // Use CoolSolveRunner to handle the execution pipeline
    coolsolve::CoolSolveRunner runner(inputFile);

    // Determine the debug output path early so we can pass it to the solver options
    fs::path debugPath;
    if (debugMode) {
        if (debugDir.empty()) {
            fs::path inputPath(inputFile);
            debugPath = inputPath.parent_path() / (inputPath.stem().string() + "_coolsolve");
        } else {
            debugPath = debugDir;
        }
    }
    
    // Clean up previous solve outputs before starting a new solve
    {
        fs::path inputPath(inputFile);
        // Delete existing .sol file
        fs::path solPath = inputPath.parent_path() / (inputPath.stem().string() + ".sol");
        if (fs::exists(solPath)) {
            fs::remove(solPath);
        }
        // Delete existing debug directory
        if (!debugPath.empty() && fs::exists(debugPath)) {
            std::error_code ec;
            fs::remove_all(debugPath, ec);
        }
    }

    // If debug mode and symbolic reduction enabled, write the reduction debug report
    // into the debug directory (the directory is created by runner.generateDebugOutput later,
    // but writeDebugReductionReport will create it on its own if needed)
    if (debugMode && options.enableSymbolicReduction) {
        std::error_code ec;
        fs::create_directories(debugPath, ec);
        options.debugReductionPath = (debugPath / "symbolic_reduction.md").string();
    }
    
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
    if (options.enableSymbolicReduction && !solveResult.blockResults.empty()) {
        int maxReducedSize = 0;
        int totalBlocksReduced = 0;
        int totalVarsReduced = 0;
        for (const auto& br : solveResult.blockResults) {
            if (br.symbolicReductionApplied) {
                ++totalBlocksReduced;
                totalVarsReduced += (br.originalSize - br.reducedSize);
                maxReducedSize = std::max(maxReducedSize, br.reducedSize);
            } else {
                maxReducedSize = std::max(maxReducedSize, br.originalSize);
            }
        }
        if (totalBlocksReduced > 0) {
            std::cout << "Largest block after reduction: " << maxReducedSize
                      << " (" << totalBlocksReduced << " block(s) reduced, "
                      << totalVarsReduced << " variable(s) eliminated)\n";
        }
    }
    
    // Report Solver Result
    if (!solveResult.success && solveResult.status != coolsolve::SolverStatus::InvalidInput) {
        std::cout << "\n=== Solver Error ===\n";
        std::cout << "Status: " << coolsolve::statusToString(solveResult.status) << "\n";
        std::cout << "Message: " << solveResult.errorMessage << "\n";
    } else if (solveResult.success) {
        std::cout << "\nSolver: SUCCESS (" << solveResult.totalIterations << " iterations)\n";
    }
    
    // Always verify the solution when the solver converges.
    // This catches CoolProp evaluation failures in the final solution
    // (e.g. unphysical inputs that were clamped during iteration).
    bool solutionValid = solveResult.success;
    if (solveResult.success) {
        auto checkResult = coolsolve::checkSolution(
            ir, solveResult.variables, solveResult.stringVariables,
            options.coolpropConfig);
        
        if (!checkResult.allSatisfied) {
            solutionValid = false;
            std::cerr << "\n=== Solution Verification Failed ===\n";
            // Show CoolProp errors from the checker
            for (const auto& d : checkResult.diagnostics.items()) {
                if (d.severity == coolsolve::DiagnosticSeverity::Error) {
                    std::cerr << "error: " << d.message << "\n";
                }
            }
            // Show violated equations
            for (const auto& chk : checkResult.checks) {
                if (!chk.satisfied) {
                    std::cerr << "  equation: " << chk.originalText
                              << " (residual=" << std::scientific << std::setprecision(2) << chk.residual << ")\n";
                }
            }
        }
        
        if (debugMode) {
            coolsolve::printSolutionCheckReport(checkResult);
            coolsolve::writeSolutionCheckReport(
                (debugPath / "solution_check.md").string(), checkResult);
        }
    }
    
    // Generate standalone LaTeX report if enableLatexReport is set (non-debug mode)
    // In debug mode the report is already generated inside generateDebugOutput().
    if (solutionValid && options.enableLatexReport && !debugMode) {
        fs::path inputPath(inputFile);
        std::string stem = inputPath.stem().string();
        fs::path texPath = inputPath.parent_path() / (stem + "_report.tex");
        if (coolsolve::writeLatexReport(
                texPath.string(), ir, analysisResult, solveResult,
                runner.getTiming(), stem, inputFile)) {
            std::cout << "LaTeX report: " << fs::weakly_canonical(texPath) << "\n";
        } else {
            std::cerr << "Warning: Could not write LaTeX report: " << texPath << "\n";
        }
    }
    
    // Write .sol file if successful
    if (solutionValid && writeSolFile) {
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
    if (solutionValid && updateGuessFile) {
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
    
    // Generate output based on format
    std::string output;
    if (format == "json") {
        output = coolsolve::generateAnalysisJSON(ir, analysisResult);
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

    // Print diagnostics (warnings, info messages) if any
    const auto& diagnostics = runner.getDiagnostics();
    if (diagnostics.size() > 0) {
        // Separate C001 CoolProp warnings (can be very numerous) from other diagnostics
        std::map<std::string, int> c001Counts;  // unique message -> count
        std::vector<coolsolve::Diagnostic> otherDiag;
        int totalC001 = 0;
        for (const auto& d : diagnostics.items()) {
            if (d.code == "C001") {
                c001Counts[d.message]++;
                totalC001++;
            } else {
                otherDiag.push_back(d);
            }
        }
        
        // Print non-C001 diagnostics (errors, info, other warnings)
        for (const auto& d : otherDiag) {
            const char* prefix = "";
            switch (d.severity) {
                case coolsolve::DiagnosticSeverity::Error:   prefix = "\033[31merror\033[0m"; break;
                case coolsolve::DiagnosticSeverity::Warning: prefix = "\033[33mwarning\033[0m"; break;
                case coolsolve::DiagnosticSeverity::Info:    prefix = "\033[36minfo\033[0m"; break;
            }
            std::cerr << prefix;
            if (d.line > 0) std::cerr << " (line " << d.line << ")";
            std::cerr << ": " << d.message << "\n";
        }
        
        // Print C001 summary instead of flooding the terminal
        if (totalC001 > 0) {
            std::cerr << "\033[33mwarning\033[0m: " << totalC001 << " CoolProp warning(s) during solving ("
                      << c001Counts.size() << " unique)";
            if (!debugMode) {
                std::cerr << ". Use -d for details";
            }
            std::cerr << "\n";
        }
    }

    return solutionValid ? 0 : 1;
}
