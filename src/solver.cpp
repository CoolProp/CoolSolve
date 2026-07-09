#include "coolsolve/solver.h"
#include "coolsolve/solver_common.h"
#include "coolsolve/symbolic_reduction.h"
#include "coolsolve/variable_inference.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <csignal>
#include <atomic>
#include <thread>
#include <mutex>
#include <future>
#include <memory>
#include <optional>
#include <limits>

#ifdef __unix__
#include <unistd.h>
#endif

namespace coolsolve {

// ============================================================================
// Utility Functions
// ============================================================================

std::string statusToString(SolverStatus status) {
    switch (status) {
        case SolverStatus::Success: return "Success";
        case SolverStatus::MaxIterations: return "MaxIterations";
        case SolverStatus::LineSearchFailed: return "LineSearchFailed";
        case SolverStatus::SingularJacobian: return "SingularJacobian";
        case SolverStatus::InvalidInput: return "InvalidInput";
        case SolverStatus::EvaluationError: return "EvaluationError";
        case SolverStatus::Diverged: return "Diverged";
        case SolverStatus::ParseFailed: return "ParseFailed";
        default: return "Unknown";
    }
}

std::string SolverTrace::toString() const {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6);
    
    // Include solver type in header if available
    if (!solverType.empty()) {
        ss << "Solver: " << solverType << " | ";
    }
    ss << "Iterations: " << iterations.size() << " | Status: "
       << statusToString(finalStatus) << "\n";
    ss << "Total time: " << totalTime.count() << " s\n";
    
    // Show tear variable names if available
    if (!tearVarNames.empty()) {
        ss << "Tear variables: ";
        for (size_t i = 0; i < tearVarNames.size(); ++i) {
            if (i > 0) ss << ", ";
            ss << tearVarNames[i];
        }
        ss << "\n";
    }
    
    ss << std::setw(6) << "Iter"
       << std::setw(15) << "||F||"
       << std::setw(15) << "||dx||"
       << std::setw(12) << "lambda" << "\n";
    ss << std::string(48, '-') << "\n";
    
    for (const auto& it : iterations) {
        ss << std::setw(6) << it.iter
           << std::setw(15) << it.residualNorm
           << std::setw(15) << it.stepNorm
           << std::setw(12) << it.lambda << "\n";
        // Print per-iteration detail (e.g. tear variable values, inner solver progress)
        if (!it.detail.empty()) {
            ss << it.detail;
        }
    }
    return ss.str();
}

// ============================================================================
// Error Categorization
// ============================================================================

ErrorCategory categorizeError(const std::string& errorMsg) {
    if (errorMsg.empty()) return ErrorCategory::None;
    if (errorMsg.find("TIMEOUT") != std::string::npos || errorMsg.find("timeout") != std::string::npos) 
        return ErrorCategory::Timeout;
    if (errorMsg.find("Max iterations") != std::string::npos) 
        return ErrorCategory::MaxIterations;
    if (errorMsg.find("Line search failed") != std::string::npos) 
        return ErrorCategory::LineSearchFailed;
    if (errorMsg.find("Jacobian is rank-deficient") != std::string::npos || errorMsg.find("SingularJacobian") != std::string::npos) 
        return ErrorCategory::SingularJacobian;
    if (errorMsg.find("CoolProp") != std::string::npos) 
        return ErrorCategory::CoolPropError;
    if (errorMsg.find("Undefined variable") != std::string::npos) 
        return ErrorCategory::UndefinedVariable;
    if (errorMsg.find("Unknown or unsupported function") != std::string::npos) 
        return ErrorCategory::UnsupportedFunction;
    if (errorMsg.find("Evaluation failed") != std::string::npos || errorMsg.find("EvaluationError") != std::string::npos) 
        return ErrorCategory::EvaluationError;
    return ErrorCategory::Other;
}

std::string categoryToString(ErrorCategory category) {
    switch (category) {
        case ErrorCategory::None: return "None";
        case ErrorCategory::Timeout: return "Timeout";
        case ErrorCategory::Converged: return "Converged";
        case ErrorCategory::MaxIterations: return "Max iterations";
        case ErrorCategory::LineSearchFailed: return "Line search failed";
        case ErrorCategory::SingularJacobian: return "Singular Jacobian";
        case ErrorCategory::CoolPropError: return "CoolProp error";
        case ErrorCategory::UndefinedVariable: return "Undefined variable";
        case ErrorCategory::UnsupportedFunction: return "Unsupported function";
        case ErrorCategory::EvaluationError: return "Evaluation error";
        case ErrorCategory::Other: return "Other";
        default: return "Unknown";
    }
}

// ============================================================================
// Solver Strategy Utilities
// ============================================================================

std::string strategyToString(SolverStrategy strategy) {
    switch (strategy) {
        case SolverStrategy::Newton:            return "Newton";
        case SolverStrategy::TrustRegion:       return "TrustRegion";
        case SolverStrategy::LevenbergMarquardt: return "LevenbergMarquardt";
        case SolverStrategy::Partitioned:       return "Partitioned";
        case SolverStrategy::BisectionND:       return "BisectionND";
        case SolverStrategy::Homotopy:          return "Homotopy";
        case SolverStrategy::Kinsol:            return "Kinsol";
        default:                                return "Unknown";
    }
}

bool parseStrategy(const std::string& name, SolverStrategy& out) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "newton")                          { out = SolverStrategy::Newton; return true; }
    if (lower == "trustregion" || lower == "trust_region" || lower == "tr")
                                                    { out = SolverStrategy::TrustRegion; return true; }
    if (lower == "levenbergmarquardt" || lower == "levenberg_marquardt" || lower == "lm")
                                                    { out = SolverStrategy::LevenbergMarquardt; return true; }
    if (lower == "partitioned")                     { out = SolverStrategy::Partitioned; return true; }
    if (lower == "bisectionnd" || lower == "bisection_nd" || lower == "bisection")
                                                    { out = SolverStrategy::BisectionND; return true; }
    if (lower == "homotopy" || lower == "continuation")
                                                    { out = SolverStrategy::Homotopy; return true; }
    if (lower == "kinsol" || lower == "kin")        { out = SolverStrategy::Kinsol; return true; }
    return false;
}

std::string pipelineModeToString(SolverPipelineMode mode) {
    switch (mode) {
        case SolverPipelineMode::Sequential: return "Sequential";
        case SolverPipelineMode::Parallel:   return "Parallel";
        default:                             return "Unknown";
    }
}

std::unique_ptr<NonLinearSolver> createSolver(SolverStrategy strategy) {
    switch (strategy) {
        case SolverStrategy::Newton:
            return std::make_unique<NewtonSolver>();
        case SolverStrategy::TrustRegion:
            return std::make_unique<TrustRegionSolver>();
        case SolverStrategy::LevenbergMarquardt:
            return std::make_unique<LevenbergMarquardtSolver>();
        case SolverStrategy::BisectionND:
            return std::make_unique<BisectionNDSolver>();
        case SolverStrategy::Homotopy:
            return std::make_unique<HomotopySolver>();
        case SolverStrategy::Kinsol:
            return std::make_unique<KINSOLSolver>();
        case SolverStrategy::Partitioned:
            // Partitioned is handled specially by the orchestrator (needs structural info).
            return nullptr;
        default:
            return nullptr;
    }
}

// ============================================================================
// Config file loading
// ============================================================================

static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end == std::string::npos ? std::string::npos : end - start + 1);
}

static bool parseBool(const std::string& v) {
    std::string s = v;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (s == "true" || s == "1" || s == "yes") return true;
    if (s == "false" || s == "0" || s == "no") return false;
    return std::stoi(v) != 0;
}

bool loadSolverOptionsFromFile(const std::string& path, SolverOptions& options) {
    std::ifstream f(path);
    if (!f.is_open()) return false;
    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));
        if (key.empty()) continue;
        try {
            if (key == "maxIterations") options.maxIterations = std::stoi(val);
            else if (key == "tolerance") options.tolerance = std::stod(val);
            else if (key == "relativeTolerance") options.relativeTolerance = std::stod(val);
            else if (key == "stepTolerance") options.stepTolerance = std::stod(val);
            else if (key == "verbose") options.verbose = parseBool(val);
            else if (key == "lsAlpha") options.lsAlpha = std::stod(val);
            else if (key == "lsRho") options.lsRho = std::stod(val);
            else if (key == "lsMaxIterations") options.lsMaxIterations = std::stoi(val);
            else if (key == "lsMinStep") options.lsMinStep = std::stod(val);
            else if (key == "lsRelaxedTolerance") options.lsRelaxedTolerance = std::stod(val);
            else if (key == "lsNonMonotoneMemory") options.lsNonMonotoneMemory = std::stoi(val);
            else if (key == "enableScaling") options.enableScaling = parseBool(val);
            else if (key == "broydenRecomputeInterval") options.broydenRecomputeInterval = std::stoi(val);
            else if (key == "trInitialRadius") options.trInitialRadius = std::stod(val);
            else if (key == "trMaxRadius") options.trMaxRadius = std::stod(val);
            else if (key == "trEta") options.trEta = std::stod(val);
            else if (key == "trShrinkFactor") options.trShrinkFactor = std::stod(val);
            else if (key == "trGrowFactor") options.trGrowFactor = std::stod(val);
            else if (key == "trAdaptiveRadius") options.trAdaptiveRadius = parseBool(val);
            else if (key == "trBroydenRecomputeInterval") options.trBroydenRecomputeInterval = std::stoi(val);
            else if (key == "trBroydenRestartRejects") options.trBroydenRestartRejects = std::stoi(val);
            else if (key == "partitionedMaxIterations") options.partitionedMaxIterations = std::stoi(val);
            else if (key == "partitionedRelaxation") options.partitionedRelaxation = std::stod(val);
            else if (key == "partitionedMinDiagonal") options.partitionedMinDiagonal = std::stod(val);
            else if (key == "partitionedMinBlockSize") options.partitionedMinBlockSize = std::stoi(val);
            else if (key == "enableTearing") options.enableTearing = parseBool(val);
            else if (key == "enableSymbolicReduction") options.enableSymbolicReduction = parseBool(val);
            else if (key == "debugReductionPath") options.debugReductionPath = val;
            else if (key == "tearingMaxIterations") options.tearingMaxIterations = std::stoi(val);
            else if (key == "tearingMinBlockSize") options.tearingMinBlockSize = std::stoi(val);
            else if (key == "tearingInnerIterations") options.tearingInnerIterations = std::stoi(val);
            else if (key == "timeoutSeconds") options.timeoutSeconds = std::stoi(val);
            // Levenberg-Marquardt options
            else if (key == "lmInitialLambda") options.lmInitialLambda = std::stod(val);
            else if (key == "lmLambdaIncrease") options.lmLambdaIncrease = std::stod(val);
            else if (key == "lmLambdaDecrease") options.lmLambdaDecrease = std::stod(val);
            else if (key == "lmMinLambda") options.lmMinLambda = std::stod(val);
            else if (key == "lmMaxLambda") options.lmMaxLambda = std::stod(val);
            else if (key == "lmNielsenUpdate") options.lmNielsenUpdate = parseBool(val);
            else if (key == "lmGeodesicAcceleration") options.lmGeodesicAcceleration = parseBool(val);
            // Solver pipeline configuration
            else if (key == "solverPipeline") {
                // Parse comma-separated list of solver names
                // e.g. "Newton, TrustRegion, LM"
                options.solverPipeline.clear();
                std::istringstream ss(val);
                std::string token;
                while (std::getline(ss, token, ',')) {
                    token = trim(token);
                    SolverStrategy strat;
                    if (parseStrategy(token, strat)) {
                        options.solverPipeline.push_back(strat);
                    }
                }
            }
            else if (key == "pipelineMode") {
                std::string lower = val;
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                if (lower == "parallel") {
                    options.pipelineMode = SolverPipelineMode::Parallel;
                } else {
                    options.pipelineMode = SolverPipelineMode::Sequential;
                }
            }
            // BisectionND options
            else if (key == "bisectionNDMaxBlockSize") options.bisectionNDMaxBlockSize = std::stoi(val);
            else if (key == "bisectionNDIterFactor")   options.bisectionNDIterFactor   = std::stod(val);
            // KINSOL options
            else if (key == "kinsolGlobalStrategy") {
                std::string lower = val;
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                if (lower == "linesearch" || lower == "line_search" || lower == "ls" || lower == "newton")
                    options.kinsolGlobalStrategy = KinsolGlobalStrategy::LineSearch;
                else if (lower == "picard" || lower == "richardson")
                    options.kinsolGlobalStrategy = KinsolGlobalStrategy::Picard;
                else if (lower == "fp" || lower == "fixedpoint" || lower == "fixed_point" || lower == "anderson")
                    options.kinsolGlobalStrategy = KinsolGlobalStrategy::FixedPoint;
                else {
                    std::cerr << "[Warning] kinsolGlobalStrategy='" << val
                              << "' is not one of {linesearch, picard, fp}; "
                              << "falling back to default (linesearch).\n";
                    options.kinsolGlobalStrategy = KinsolGlobalStrategy::LineSearch;
                }
            }
            else if (key == "kinsolLineSearchAlpha") {
                double v = std::stod(val);
                if (v <= 0.0 || v >= 1.0) {
                    std::cerr << "[Warning] kinsolLineSearchAlpha=" << v
                              << " out of range (0,1); falling back to default (1e-4).\n";
                    options.kinsolLineSearchAlpha = 1e-4;
                } else {
                    options.kinsolLineSearchAlpha = v;
                }
            }
            else if (key == "kinsolLineSearchMaxIters") options.kinsolLineSearchMaxIters = std::stoi(val);
            else if (key == "kinsolPicardOmega") {
                double v = std::stod(val);
                if (v <= 0.0) {
                    std::cerr << "[Warning] kinsolPicardOmega=" << v
                              << " must be > 0; falling back to default (1.0).\n";
                    options.kinsolPicardOmega = 1.0;
                } else {
                    options.kinsolPicardOmega = v;
                }
            }
            else if (key == "kinsolAndersonDepth")   options.kinsolAndersonDepth   = std::stoi(val);
            else if (key == "kinsolAndersonRelaxation") {
                double v = std::stod(val);
                if (v <= 0.0 || v > 1.0) {
                    std::cerr << "[Warning] kinsolAndersonRelaxation=" << v
                              << " out of range (0,1]; falling back to default (1.0).\n";
                    options.kinsolAndersonRelaxation = 1.0;
                } else {
                    options.kinsolAndersonRelaxation = v;
                }
            }
            // Multi-start fallback (roadmap §4.2)
            else if (key == "multiStartEnabled") {
                options.multiStartEnabled = parseBool(val);
            }
            else if (key == "multiStartMaxRestarts") {
                int v = std::stoi(val);
                if (v < 0) {
                    std::cerr << "[Warning] multiStartMaxRestarts=" << v
                              << " is negative; falling back to default (4).\n";
                    options.multiStartMaxRestarts = 4;
                } else {
                    options.multiStartMaxRestarts = v;
                }
            }
            else if (key == "multiStartNumCores") {
                int v = std::stoi(val);
                if (v < 0) {
                    std::cerr << "[Warning] multiStartNumCores=" << v
                              << " is negative; falling back to default (1 = sequential).\n";
                    options.multiStartNumCores = 1;
                } else {
                    options.multiStartNumCores = v;
                }
            }
            // CoolProp integration options
            else if (key == "coolpropBackend") options.coolpropConfig.backend = val;
            else if (key == "coolpropUseAbstractState") options.coolpropConfig.useAbstractState = parseBool(val);
            else if (key == "coolpropEnableAnalyticalDerivatives") options.coolpropConfig.enableAnalyticalDerivatives = parseBool(val);
            else if (key == "coolpropCacheEnabled") options.coolpropConfig.cacheEnabled = parseBool(val);
            else if (key == "coolpropEnableSuperancillaries") options.coolpropConfig.enableSuperancillaries = parseBool(val);
            // LaTeX report options
            else if (key == "enableLatexReport") options.enableLatexReport = parseBool(val);
            else if (key == "latexCompiler") options.latexCompiler = val;
            // --- Equation-based integration ---
            else if (key == "integralMethod") {
                std::string lower = val;
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                options.integralMethod = lower;  // validated by makeIntegratorOptions
            }
            else if (key == "integralFixedStep")      options.integralFixedStep = std::stod(val);
            else if (key == "integralMaxSteps")       options.integralMaxSteps = std::stoi(val);
            else if (key == "integralRelTol")         options.integralRelTol = std::stod(val);
            else if (key == "integralAbsTol")         options.integralAbsTol = std::stod(val);
            else if (key == "integralMinStep")        options.integralMinStep = std::stod(val);
            else if (key == "integralMaxStep")        options.integralMaxStep = std::stod(val);
            else if (key == "integralRichardson")     options.integralRichardson = parseBool(val);
            else if (key == "integralOutputInterval") options.integralOutputInterval = std::stod(val);
        } catch (...) {
            // Ignore malformed values
        }
    }
    return true;
}

// ============================================================================
// Timeout Handling
// ============================================================================

// Thread-local chrono-based timeout: each thread maintains its own deadline.
// This allows parallel solver invocations (e.g. robustness testing) where
// each thread has an independent timeout.
static thread_local std::chrono::steady_clock::time_point tl_deadline{};
static thread_local bool tl_has_deadline{false};

TimeoutGuard::TimeoutGuard(int seconds) : seconds_(seconds) {
    if (seconds > 0) {
        tl_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(seconds);
        tl_has_deadline = true;
    } else {
        tl_has_deadline = false;
    }
}

TimeoutGuard::~TimeoutGuard() {
    tl_has_deadline = false;
}

bool TimeoutGuard::hasTimedOut() {
    if (tl_has_deadline) {
        return std::chrono::steady_clock::now() >= tl_deadline;
    }
    return false;
}

// Newton, TrustRegion, and LM solver implementations are in separate files:

// Newton, TrustRegion, and LevenbergMarquardt implementations moved to:
//   solver_newton.cpp, solver_trust_region.cpp, solver_lm.cpp

// ============================================================================
// Solver (Orchestrator) Implementation
// ============================================================================

Solver::Solver(const IR& ir, 
               const StructuralAnalysisResult& analysis,
               const CoolPropConfig& config)
    : evaluator_(ir, analysis, config)
    , ir_(ir)
    , analysis_(analysis) {
    // Initialize from guesses in IR
    evaluator_.initializeFromGuesses();
}

void Solver::setGuess(const std::string& name, double value) {
    evaluator_.setVariableValue(name, value);
}

void Solver::setStringVariable(const std::string& name, const std::string& value) {
    evaluator_.setStringVariableValue(name, value);
}

bool Solver::tryExplicitSolve(size_t blockIndex) {
    const Block& block = analysis_.blocks[blockIndex];
    if (block.equationIds.size() != 1 || block.variables.size() != 1)
        return false;

    const int eqId = block.equationIds[0];
    const std::string& outVar = block.variables[0];
    const auto& equations = ir_.getEquations();
    if (eqId < 0 || eqId >= static_cast<int>(equations.size()))
        return false;

    const EquationInfo& eq = equations[eqId];
    // Structurally explicit: equation must be "x = expr(others)" so we can set x = value(RHS).
    // 1) Output variable must not appear in the RHS.
    if (eq.rhsVariables.count(outVar) != 0)
        return false;
    // 2) LHS must be exactly the output variable (not e.g. x^2 or f(x)).
    if (!eq.lhs || !eq.lhs->is<Variable>())
        return false;
    const std::string& lhsName = eq.lhs->as<Variable>().flattenedName();
    auto caseInsensitiveEqual = [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (std::tolower(static_cast<unsigned char>(a[i])) !=
                std::tolower(static_cast<unsigned char>(b[i])))
                return false;
        }
        return true;
    };
    if (!caseInsensitiveEqual(lhsName, outVar))
        return false;

    BlockEvaluator& blockEval = evaluator_.getBlock(blockIndex);

    std::map<std::string, double> externalVars;
    for (const auto& [name, value] : evaluator_.getAllVariables()) {
        if (!caseInsensitiveEqual(name, outVar))
            externalVars[name] = value;
    }
    std::map<std::string, std::string> externalStringVars;
    for (const auto& [name, value] : evaluator_.getAllStringVariables()) {
        externalStringVars[name] = value;
    }

    const double xCurrent = evaluator_.getVariableValue(outVar);
    try {
        EvaluationResult result = blockEval.evaluate({xCurrent}, externalVars, externalStringVars);
        if (result.residuals.empty()) return false;
        // Residual is LHS - RHS; for x = expr we have residual = x - expr, so x_solution = expr = x - residual
        const double xSolution = xCurrent - result.residuals[0];
        evaluator_.setVariableValue(outVar, xSolution);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

SolverStatus Solver::solveBlock(size_t blockIndex, 
                                const SolverOptions& options,
                                SolverTrace* trace,
                                std::string* outErrorMessage) {
    BlockEvaluator& blockEval = evaluator_.getBlock(blockIndex);
    const auto& varNames = blockEval.getVariables();
    const size_t n = varNames.size();
    
    if (n == 0) {
        // Zero-unknown block: this is a "check equation" where all variables are already determined.
        // We still need to verify the residual is near zero - if not, the system is overdetermined.
        
        // Get all current variable values as external vars for evaluation
        std::map<std::string, double> externalVars;
        for (const auto& [name, value] : evaluator_.getAllVariables()) {
            externalVars[name] = value;
        }
        std::map<std::string, std::string> externalStringVars;
        for (const auto& [name, value] : evaluator_.getAllStringVariables()) {
            externalStringVars[name] = value;
        }
        
        // Evaluate the block's equations with all variables as external
        try {
            std::vector<double> emptyX;  // No unknowns
            auto evalResult = blockEval.evaluate(emptyX, externalVars, externalStringVars);
            
            // Check that all residuals are within tolerance
            for (size_t i = 0; i < evalResult.residuals.size(); ++i) {
                double residual = std::abs(evalResult.residuals[i]);
                if (residual > options.tolerance) {
                    if (options.verbose) {
                        std::cerr << "Block " << blockIndex << " (zero unknowns): residual["
                                  << i << "] = " << residual << " exceeds tolerance " 
                                  << options.tolerance << std::endl;
                    }
                    if (outErrorMessage) {
                        std::ostringstream ss;
                        ss << "Check equation failed: residual = " << residual 
                           << " (tolerance = " << options.tolerance << "). "
                           << "This indicates an overdetermined or inconsistent system - "
                           << "the same variable may be defined by multiple conflicting equations.";
                        *outErrorMessage = ss.str();
                    }
                    return SolverStatus::EvaluationError;
                }
            }
            
            if (options.verbose) {
                std::cout << "Block " << blockIndex << " (zero unknowns): check equation satisfied" << std::endl;
            }
        } catch (const std::exception& e) {
            if (outErrorMessage) {
                *outErrorMessage = std::string("Check equation evaluation failed: ") + e.what();
            }
            return SolverStatus::EvaluationError;
        }
        
        return SolverStatus::Success;
    }
    
    // Debug: Print block info
    if (options.verbose) {
        std::cout << "=== Block " << blockIndex << " ===" << std::endl;
        std::cout << "Variables: ";
        for (const auto& v : varNames) std::cout << v << " ";
        std::cout << std::endl;
    }
    
    // Gather initial guess from current state
    Eigen::VectorXd x(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = evaluator_.getVariableValue(varNames[i]);
        if (options.verbose) {
            std::cout << "  Initial " << varNames[i] << " = " << x[i] << std::endl;
        }
    }
    
    // Helper for case-insensitive string comparison
    auto caseInsensitiveEqual = [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (std::tolower(static_cast<unsigned char>(a[i])) != 
                std::tolower(static_cast<unsigned char>(b[i]))) {
                return false;
            }
        }
        return true;
    };
    
    // Convert map to std::map<string, double> for evaluation
    std::map<std::string, double> externalVars;
    for (const auto& [name, value] : evaluator_.getAllVariables()) {
        // Add only if not in this block
        bool inBlock = false;
        for (const auto& bvar : varNames) {
            if (caseInsensitiveEqual(name, bvar)) {
                inBlock = true;
                break;
            }
        }
        if (!inBlock) {
            externalVars[name] = value;
            if (options.verbose) {
                std::cout << "  External " << name << " = " << value << std::endl;
            }
        }
    }
    
    // Get external string variables
    std::map<std::string, std::string> externalStringVars;
    for (const auto& [name, value] : evaluator_.getAllStringVariables()) {
        externalStringVars[name] = value;
    }
    
    // For size-1 blocks, check if it's truly implicit or just needs direct evaluation
    if (n == 1 && tryExplicitSolve(blockIndex)) {
        if (trace) {
            trace->solverType = "Explicit";
            trace->finalStatus = SolverStatus::Success;
        }
        return SolverStatus::Success;
    }

    // For size-1 implicit blocks, use a specialized 1D root-finding solver.
    // Standard Newton with line search can fail when the initial guess is very
    // far from the solution (e.g. enthalpy equations where default guess is 1
    // but solution is ~1e5). This solver uses:
    //   1. Newton steps with adaptive trust-region limiting
    //   2. Sign-change detection → bisection fallback
    //   3. Multiple starting point exploration if Newton diverges
    if (n == 1) {
        auto startTime1D = std::chrono::high_resolution_clock::now();
        if (trace) {
            if (trace->solverType.empty())
                trace->solverType = "Newton1D";
            else
                trace->solverType += " -> Newton1D";
        }
        
        // Lambda for evaluating the 1D residual+derivative
        Newton1DSolver::Eval1D eval1D = [&](double xval) -> std::pair<double, double> {
            std::vector<double> x_std = {xval};
            auto er = blockEval.evaluate(x_std, externalVars, externalStringVars);
            double f = er.residuals[0];
            double j = (er.jacobian.size() > 0 && er.jacobian[0].size() > 0) ? er.jacobian[0][0] : 0.0;
            return {f, j};
        };
        
        double xCur = x[0];
        std::string errorMsg;
        SolverStatus status1D = Newton1DSolver::solve(
            eval1D, xCur, externalVars, options, trace, &errorMsg);
        
        if (status1D == SolverStatus::Success) {
            x[0] = xCur;
            evaluator_.setVariableValue(varNames[0], xCur);
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime1D;
            }
            return SolverStatus::Success;
        }
        
        if (status1D == SolverStatus::EvaluationError) {
            // Fatal evaluation error — propagate immediately
            if (outErrorMessage) *outErrorMessage = errorMsg;
            return SolverStatus::EvaluationError;
        }
        
        // Newton1D failed — reset x and fall through to standard pipeline
        x[0] = evaluator_.getVariableValue(varNames[0]);
        if (trace) {
            trace->iterations.clear(); // Clear Newton1D iterations for clean pipeline trace
        }
    }

    // ========================================================================
    // Symbolic Block Reduction (optional pre-processing)
    // ========================================================================
    // When enabled, attempt to reduce the block size before iterative solving.
    // This can turn e.g. a 3-variable CoolProp block into three size-1 blocks.
    if (options.enableSymbolicReduction && n > 1) {
        const Block& block = analysis_.blocks[blockIndex];
        BlockReductionResult reduction = analyseBlockReduction(block, ir_, analysis_);

        if (reduction.reduced && reduction.reducedVariables.size() < n) {
            if (options.verbose) {
                std::cout << "Block " << blockIndex << ": symbolic reduction "
                          << n << " → " << reduction.reducedVariables.size()
                          << " vars (inversions=" << reduction.inversionsApplied
                          << " extractions=" << reduction.extractionsApplied
                          << " substitutions=" << reduction.substitutionsApplied
                          << ")" << std::endl;
            }

            // Record reduction info in the trace (for reporting)
            if (trace) {
                trace->symbolicReductionApplied = true;
                trace->originalBlockSize = static_cast<int>(n);
                trace->reducedBlockSize = static_cast<int>(reduction.reducedVariables.size());
                trace->symInversions = reduction.inversionsApplied;
                trace->symExtractions = reduction.extractionsApplied;
                trace->symSubstitutions = reduction.substitutionsApplied;

                // Record original equations
                for (int eqId : block.equationIds) {
                    if (eqId >= 0 && eqId < static_cast<int>(ir_.getEquations().size())) {
                        trace->originalEquations.push_back(
                            ir_.getEquations()[eqId].originalText);
                    }
                }

                // Record reduction step descriptions
                for (const auto& step : reduction.preSolveSteps) {
                    std::string desc = "PRE: " + step.variable + " = ";
                    if (step.inverted) {
                        desc += step.invertedFuncName + "(" + step.fluidName + ", ...)";
                        desc += "  [CoolProp inversion]";
                    } else {
                        int eqIdx = step.equationId;
                        if (eqIdx >= 0 && eqIdx < static_cast<int>(ir_.getEquations().size())) {
                            desc += ir_.getEquations()[eqIdx].originalText;
                        }
                        desc += "  [explicit extraction]";
                    }
                    trace->reductionStepDescriptions.push_back(desc);
                }
                for (const auto& step : reduction.postSolveSteps) {
                    std::string desc = "POST: " + step.variable + " = ";
                    if (step.inverted) {
                        desc += step.invertedFuncName + "(" + step.fluidName + ", ...)";
                        desc += "  [CoolProp inversion]";
                    } else {
                        int eqIdx = step.equationId;
                        if (eqIdx >= 0 && eqIdx < static_cast<int>(ir_.getEquations().size())) {
                            desc += ir_.getEquations()[eqIdx].originalText;
                        }
                        desc += "  [explicit extraction]";
                    }
                    trace->reductionStepDescriptions.push_back(desc);
                }

                // Record remaining (reduced) equations
                for (int eqId : reduction.reducedEquationIds) {
                    if (eqId >= 0 && eqId < static_cast<int>(ir_.getEquations().size())) {
                        trace->reducedEquations.push_back(
                            ir_.getEquations()[eqId].originalText);
                    }
                }
            }

            // Helper: evaluate an expression in the current evaluator state
            auto evalExpr = [&](const ExprPtr& expr) -> double {
                ExpressionEvaluator exprEval(0, options.coolpropConfig);
                if (lookupTableStore_) exprEval.setLookupTableStore(lookupTableStore_);
                // Load all current variable values
                for (const auto& [name, val] : evaluator_.getAllVariables()) {
                    exprEval.setVariable(name, ADValue::constant(val, 0));
                }
                for (const auto& [name, val] : evaluator_.getAllStringVariables()) {
                    exprEval.setStringVariable(name, val);
                }
                // Register user functions and procedures
                for (const auto& func : ir_.getFunctions()) {
                    exprEval.registerFunction(func);
                }
                for (const auto& proc : ir_.getProcedures()) {
                    exprEval.registerProcedure(proc);
                }
                return exprEval.evaluate(expr).value;
            };

            // Helper: evaluate a CoolProp inversion step
            auto evalInvertedStep = [&](const ReductionStep& step) -> double {
                // Build a FunctionCall AST node for the inverted call
                FunctionCall invCall;
                invCall.name = step.invertedFuncName;
                invCall.args = {step.fluidExpr};
                invCall.namedArgs = step.newNamedArgs;

                auto expr = std::make_shared<Expression>();
                expr->node = invCall;
                return evalExpr(expr);
            };

            bool reductionFailed = false;

            // Execute pre-solve steps
            for (const auto& step : reduction.preSolveSteps) {
                try {
                    double value;
                    if (step.inverted) {
                        value = evalInvertedStep(step);
                    } else {
                        // Evaluate RHS of the original equation
                        const auto& eq = ir_.getEquations()[step.equationId];
                        value = evalExpr(eq.rhs);
                    }
                    evaluator_.setVariableValue(step.variable, value);
                    externalVars[step.variable] = value;
                    if (options.verbose) {
                        std::cout << "  Pre-solve " << step.variable << " = " << value
                                  << (step.inverted ? " (CoolProp inversion)" : " (explicit)")
                                  << std::endl;
                    }
                } catch (const std::exception& e) {
                    if (options.verbose) {
                        std::cerr << "  Symbolic pre-solve failed for " << step.variable
                                  << ": " << e.what() << std::endl;
                    }
                    reductionFailed = true;
                    break;
                }
            }

            if (!reductionFailed && !reduction.reducedVariables.empty()) {
                // Build reduced block and solve it
                Block reducedBlock;
                reducedBlock.id = block.id;
                reducedBlock.equationIds = reduction.reducedEquationIds;
                reducedBlock.variables = reduction.reducedVariables;

                size_t nReduced = reducedBlock.variables.size();

                // Build reduced external vars (include pre-solved variables)
                std::map<std::string, double> reducedExternalVars;
                for (const auto& [name, value] : evaluator_.getAllVariables()) {
                    bool inReducedBlock = false;
                    for (const auto& rv : reducedBlock.variables) {
                        if (caseInsensitiveEqual(name, rv)) {
                            inReducedBlock = true;
                            break;
                        }
                    }
                    if (!inReducedBlock) {
                        reducedExternalVars[name] = value;
                    }
                }

                if (nReduced == 0) {
                    // All variables were extracted — nothing to solve
                    if (options.verbose) {
                        std::cout << "  Block fully reduced (no iterative solve needed)" << std::endl;
                    }
                } else if (nReduced == 1) {
                    // Reduced to size 1 — try explicit solve or Newton1D
                    BlockEvaluator reducedEval(reducedBlock, ir_, options.coolpropConfig);

                    // Try explicit solve first
                    bool explicitOk = false;
                    const auto& eq = ir_.getEquations()[reducedBlock.equationIds[0]];
                    if (eq.lhs && eq.lhs->is<Variable>()) {
                        const std::string& lhsVar = eq.lhs->as<Variable>().flattenedName();
                        if (caseInsensitiveEqual(lhsVar, reducedBlock.variables[0])) {
                            try {
                                double val = evalExpr(eq.rhs);
                                evaluator_.setVariableValue(reducedBlock.variables[0], val);
                                explicitOk = true;
                                if (options.verbose) {
                                    std::cout << "  Reduced explicit: " << reducedBlock.variables[0]
                                              << " = " << val << std::endl;
                                }
                            } catch (...) {}
                        }
                    }

                    if (!explicitOk) {
                        // 1D Newton solver on the reduced block
                        double xCur = evaluator_.getVariableValue(reducedBlock.variables[0]);
                        Newton1DSolver::Eval1D eval1D = [&](double xval) -> std::pair<double, double> {
                            std::vector<double> x_std = {xval};
                            auto er = reducedEval.evaluate(x_std, reducedExternalVars, externalStringVars);
                            double f = er.residuals[0];
                            double j = (er.jacobian.size() > 0 && er.jacobian[0].size() > 0) ? er.jacobian[0][0] : 0.0;
                            return {f, j};
                        };

                        SolverStatus st = Newton1DSolver::solveSimple(eval1D, xCur, options);

                        if (st == SolverStatus::Success) {
                            evaluator_.setVariableValue(reducedBlock.variables[0], xCur);
                            if (options.verbose) {
                                std::cout << "  Reduced Newton1D: " << reducedBlock.variables[0]
                                          << " = " << xCur << std::endl;
                            }
                        } else {
                            reductionFailed = true;
                        }
                    }
                } else {
                    // Reduced block still size > 1.
                    // --------------------------------------------------------
                    // Re-decomposition: check if the reduced block can be split
                    // into smaller independent sub-SCCs.
                    // --------------------------------------------------------
                    bool solvedViaSubBlocks = false;

                    auto subBlocks = StructuralAnalyzer::redecomposeBlock(
                        reduction.reducedEquationIds,
                        reduction.reducedVariables,
                        ir_, analysis_);

                    if (subBlocks.size() > 1) {
                            // The reduced block splits!
                            if (options.verbose) {
                                std::cout << "  Re-decomposition: " << nReduced
                                          << "-var block → " << subBlocks.size()
                                          << " sub-blocks [";
                                for (size_t si = 0; si < subBlocks.size(); ++si) {
                                    if (si > 0) std::cout << ", ";
                                    std::cout << subBlocks[si].variables.size();
                                }
                                std::cout << "]" << std::endl;
                            }

                            // Record in trace
                            if (trace) {
                                trace->redecompositionApplied = true;
                                trace->numSubBlocks = static_cast<int>(subBlocks.size());
                                for (const auto& sb : subBlocks) {
                                    trace->subBlockSizes.push_back(static_cast<int>(sb.variables.size()));
                                }
                            }

                            // Solve sub-blocks sequentially in topological order
                            bool subBlockFailed = false;
                            for (size_t si = 0; si < subBlocks.size() && !subBlockFailed; ++si) {
                                const auto& sb = subBlocks[si];
                                size_t subN = sb.variables.size();

                                // Build external vars for this sub-block
                                // (everything not in this sub-block)
                                std::map<std::string, double> subExternal;
                                for (const auto& [name, value] : evaluator_.getAllVariables()) {
                                    bool inSub = false;
                                    for (const auto& sv : sb.variables) {
                                        if (caseInsensitiveEqual(name, sv)) { inSub = true; break; }
                                    }
                                    if (!inSub) subExternal[name] = value;
                                }

                                if (subN == 0) {
                                    continue;
                                } else if (subN == 1) {
                                    // Try explicit or Newton1D on sub-block
                                    BlockEvaluator subEval(sb, ir_, options.coolpropConfig);
                                    bool explOk = false;
                                    const auto& eq = ir_.getEquations()[sb.equationIds[0]];
                                    if (eq.lhs && eq.lhs->is<Variable>()) {
                                        const std::string& lhsVar = eq.lhs->as<Variable>().flattenedName();
                                        if (caseInsensitiveEqual(lhsVar, sb.variables[0])) {
                                            try {
                                                double val = evalExpr(eq.rhs);
                                                evaluator_.setVariableValue(sb.variables[0], val);
                                                explOk = true;
                                            } catch (...) {}
                                        }
                                    }
                                    if (!explOk) {
                                        // Newton1D on sub-block
                                        double xCur = evaluator_.getVariableValue(sb.variables[0]);
                                        Newton1DSolver::Eval1D eval1D = [&](double xval) -> std::pair<double, double> {
                                            std::vector<double> x_std = {xval};
                                            auto er = subEval.evaluate(x_std, subExternal, externalStringVars);
                                            double f = er.residuals[0];
                                            double j = (er.jacobian.size() > 0 && er.jacobian[0].size() > 0)
                                                       ? er.jacobian[0][0] : 0.0;
                                            return {f, j};
                                        };

                                        SolverStatus st = Newton1DSolver::solveSimple(eval1D, xCur, options);

                                        if (st == SolverStatus::Success) {
                                            evaluator_.setVariableValue(sb.variables[0], xCur);
                                        } else {
                                            subBlockFailed = true;
                                        }
                                    }
                                } else {
                                    // Multi-variable sub-block — use solver pipeline
                                    BlockEvaluator subEval(sb, ir_, options.coolpropConfig);
                                    Eigen::VectorXd xs(subN);
                                    for (size_t vi = 0; vi < subN; ++vi) {
                                        xs[vi] = evaluator_.getVariableValue(sb.variables[vi]);
                                    }

                                    NonLinearSolver::Problem subProblem;
                                    subProblem.size = static_cast<int>(subN);
                                    subProblem.evaluate = [&subEval, &subExternal, &externalStringVars]
                                        (const Eigen::VectorXd& xv, Eigen::VectorXd& F,
                                         Eigen::MatrixXd& J, bool computeJ) {
                                        std::vector<double> x_std(xv.data(), xv.data() + xv.size());
                                        auto res = subEval.evaluate(x_std, subExternal,
                                                                     externalStringVars, computeJ);
                                        const size_t ne = res.residuals.size();
                                        F.resize(ne);
                                        for (size_t i = 0; i < ne; ++i) F[i] = res.residuals[i];
                                        if (computeJ) {
                                            J.resize(ne, xv.size());
                                            for (size_t i = 0; i < ne; ++i)
                                                for (size_t j = 0; j < res.jacobian[i].size(); ++j)
                                                    J(i, j) = res.jacobian[i][j];
                                        }
                                    };

                                    std::string subErr;
                                    SolverStatus sStatus = solveBlockSequential(
                                        blockIndex, subProblem, subEval,
                                        sb.variables, subExternal, externalStringVars,
                                        xs, options, nullptr, &subErr);

                                    if (sStatus == SolverStatus::Success) {
                                        for (size_t vi = 0; vi < subN; ++vi) {
                                            evaluator_.setVariableValue(sb.variables[vi], xs[vi]);
                                        }
                                    } else {
                                        subBlockFailed = true;
                                    }
                                }

                                // Propagate solved variables by updating evaluator
                                // (already done above via setVariableValue)
                            }

                            if (!subBlockFailed) {
                                solvedViaSubBlocks = true;
                                // Update reducedExternalVars with all solved values
                                for (const auto& rv : reducedBlock.variables) {
                                    reducedExternalVars[rv] = evaluator_.getVariableValue(rv);
                                }
                            }
                        }

                    if (!solvedViaSubBlocks) {
                        // Fallback: solve reduced block as single monolithic system
                        BlockEvaluator reducedEval(reducedBlock, ir_, options.coolpropConfig);

                        Eigen::VectorXd xr(nReduced);
                        for (size_t i = 0; i < nReduced; ++i) {
                            xr[i] = evaluator_.getVariableValue(reducedBlock.variables[i]);
                        }

                        NonLinearSolver::Problem reducedProblem;
                        reducedProblem.size = static_cast<int>(nReduced);
                        reducedProblem.evaluate = [&reducedEval, &reducedExternalVars, &externalStringVars]
                                                  (const Eigen::VectorXd& xv,
                                                   Eigen::VectorXd& F,
                                                   Eigen::MatrixXd& J,
                                                   bool computeJacobian) {
                            std::vector<double> x_std(xv.data(), xv.data() + xv.size());
                            auto result = reducedEval.evaluate(x_std, reducedExternalVars, externalStringVars, computeJacobian);
                            const size_t nEqs = result.residuals.size();
                            F.resize(nEqs);
                            for (size_t i = 0; i < nEqs; ++i) F[i] = result.residuals[i];
                            if (computeJacobian) {
                                J.resize(nEqs, xv.size());
                                for (size_t i = 0; i < nEqs; ++i) {
                                    for (size_t j = 0; j < result.jacobian[i].size(); ++j) {
                                        J(i, j) = result.jacobian[i][j];
                                    }
                                }
                            }
                        };

                        SolverStatus rStatus;
                        if (options.pipelineMode == SolverPipelineMode::Parallel && options.solverPipeline.size() > 1) {
                            rStatus = solveBlockParallel(blockIndex, reducedProblem, reducedEval,
                                                         reducedBlock.variables, reducedExternalVars,
                                                         externalStringVars, xr, options, trace, outErrorMessage);
                        } else {
                            rStatus = solveBlockSequential(blockIndex, reducedProblem, reducedEval,
                                                           reducedBlock.variables, reducedExternalVars,
                                                           externalStringVars, xr, options, trace, outErrorMessage);
                        }

                        if (rStatus == SolverStatus::Success) {
                            for (size_t i = 0; i < nReduced; ++i) {
                                evaluator_.setVariableValue(reducedBlock.variables[i], xr[i]);
                            }
                        } else {
                            reductionFailed = true;
                        }
                    }
                }
            } else if (!reductionFailed && reduction.reducedVariables.empty()) {
                // All variables extracted in pre-solve — nothing to solve iteratively
            }

            // Execute post-solve steps
            if (!reductionFailed) {
                for (const auto& step : reduction.postSolveSteps) {
                    try {
                        double value;
                        if (step.inverted) {
                            value = evalInvertedStep(step);
                        } else {
                            const auto& eq = ir_.getEquations()[step.equationId];
                            value = evalExpr(eq.rhs);
                        }
                        evaluator_.setVariableValue(step.variable, value);
                        if (options.verbose) {
                            std::cout << "  Post-solve " << step.variable << " = " << value
                                      << (step.inverted ? " (CoolProp inversion)" : " (explicit)")
                                      << std::endl;
                        }
                    } catch (const std::exception& e) {
                        if (options.verbose) {
                            std::cerr << "  Symbolic post-solve failed for " << step.variable
                                      << ": " << e.what() << std::endl;
                        }
                        reductionFailed = true;
                        break;
                    }
                }
            }

            if (!reductionFailed) {
                // Update the x vector for the original block variables
                for (size_t i = 0; i < n; ++i) {
                    x[i] = evaluator_.getVariableValue(varNames[i]);
                }
                if (trace) {
                    std::string desc = "SymbolicReduction(" + std::to_string(n) + "→"
                                     + std::to_string(reduction.reducedVariables.size()) + ")";
                    if (!trace->solverType.empty()) trace->solverType += " -> ";
                    trace->solverType += desc;
                    trace->finalStatus = SolverStatus::Success;
                }
                return SolverStatus::Success;
            }

            // Reduction failed: restore original state and fall through
            if (options.verbose) {
                std::cout << "Block " << blockIndex << ": Symbolic reduction failed, "
                          << "falling through to standard solvers" << std::endl;
            }
            // Restore any pre-solve values that may have been set
            for (size_t i = 0; i < n; ++i) {
                evaluator_.setVariableValue(varNames[i], x[i]);
            }
        }
    }

    if (options.enableTearing && n >= static_cast<size_t>(options.tearingMinBlockSize) && n > 1) {
        Eigen::VectorXd x_saved = x;  // Save initial guess for fallback
        SolverStatus tearStatus = solveBlockTearing(blockIndex, blockEval, varNames,
                                                    externalVars, externalStringVars,
                                                    x, options, trace, outErrorMessage);
        if (tearStatus == SolverStatus::Success) {
            for (size_t i = 0; i < n; ++i) {
                evaluator_.setVariableValue(varNames[i], x[i]);
                if (options.verbose) {
                    std::cout << "  Updated " << varNames[i] << " = " << x[i] << std::endl;
                }
            }
            return SolverStatus::Success;
        }
        // Tearing failed: restore initial guess and clear tearing iterations from trace
        x = x_saved;
        if (trace) {
            // Keep the tearing info in the solver type but clear its iterations
            // so the pipeline solver has a clean trace
            trace->iterations.clear();
        }
        if (options.verbose) {
            std::cout << "Block " << blockIndex << ": Tearing failed ("
                      << statusToString(tearStatus) << "), trying pipeline" << std::endl;
        }
    }
    
    // Create problem for solver
    NonLinearSolver::Problem problem;
    problem.size = static_cast<int>(n);
    
    // Lambda to evaluate block
    problem.evaluate = [&blockEval, &varNames, &externalVars, &externalStringVars]
                       (const Eigen::VectorXd& xv,
                        Eigen::VectorXd& F,
                        Eigen::MatrixXd& J,
                        bool computeJacobian) {
        // Convert Eigen vector to std::vector
        std::vector<double> x_std(xv.data(), xv.data() + xv.size());
        
        // Evaluate block (pass computeJacobian to skip CoolProp derivatives)
        auto result = blockEval.evaluate(x_std, externalVars, externalStringVars, computeJacobian);
        
        // Copy residuals
        const size_t nEqs = result.residuals.size();
        F.resize(nEqs);
        for (size_t i = 0; i < nEqs; ++i) {
            F[i] = result.residuals[i];
        }
        
        // Copy Jacobian if requested
        if (computeJacobian) {
            J.resize(nEqs, xv.size());
            for (size_t i = 0; i < nEqs; ++i) {
                for (size_t j = 0; j < result.jacobian[i].size(); ++j) {
                    J(i, j) = result.jacobian[i][j];
                }
            }
        }
    };
    
    // Dispatch to the configured pipeline mode
    SolverStatus status;
    if (options.pipelineMode == SolverPipelineMode::Parallel && options.solverPipeline.size() > 1) {
        status = solveBlockParallel(blockIndex, problem, blockEval, varNames,
                                    externalVars, externalStringVars, x, options,
                                    trace, outErrorMessage);
    } else {
        status = solveBlockSequential(blockIndex, problem, blockEval, varNames,
                                      externalVars, externalStringVars, x, options,
                                      trace, outErrorMessage);
    }

    // Update state with solution
    if (status == SolverStatus::Success ||
        status == SolverStatus::MaxIterations) {
        for (size_t i = 0; i < n; ++i) {
            evaluator_.setVariableValue(varNames[i], x[i]);
            if (options.verbose) {
                std::cout << "  Updated " << varNames[i] << " = " << x[i] << std::endl;
            }
        }
    }
    
    return status;
}

// ----------------------------------------------------------------------------
// runSolverStrategy – run a single solver on a block
// ----------------------------------------------------------------------------

SolverStatus Solver::runSolverStrategy(SolverStrategy strategy,
                                       size_t blockIndex,
                                       NonLinearSolver::Problem& problem,
                                       BlockEvaluator& blockEval,
                                       const std::vector<std::string>& varNames,
                                       const std::map<std::string, double>& externalVars,
                                       const std::map<std::string, std::string>& externalStringVars,
                                       Eigen::VectorXd& x,
                                       const SolverOptions& options,
                                       SolverTrace* trace,
                                       std::string* outErrorMessage) {
    const size_t n = varNames.size();

    // Partitioned solver needs structural info – handle specially
    if (strategy == SolverStrategy::Partitioned) {
        if (n < static_cast<size_t>(options.partitionedMinBlockSize)) {
            // Skip gracefully — don't return InvalidInput which looks like a hard error
            return SolverStatus::MaxIterations;
        }
        SolverOptions partOpts = options;
        if (partOpts.partitionedMaxIterations < options.maxIterations) {
            partOpts.partitionedMaxIterations = options.maxIterations;
        }
        return solveBlockPartitioned(blockIndex, blockEval, varNames,
                                     externalVars, externalStringVars,
                                     x, partOpts, trace, outErrorMessage);
    }

    // All other solvers go through the NonLinearSolver interface
    auto solver = createSolver(strategy);
    if (!solver) {
        if (outErrorMessage) {
            *outErrorMessage = "Unknown solver strategy: " + strategyToString(strategy);
        }
        return SolverStatus::InvalidInput;
    }

    // For trust-region and LM, allow more iterations
    SolverOptions solverOpts = options;
    if (strategy == SolverStrategy::TrustRegion) {
        // Allow more iterations for TR which takes conservative steps
        int minIter = (n > 10) ? 500 : 300;
        if (solverOpts.maxIterations < minIter) solverOpts.maxIterations = minIter;
    }
    if (strategy == SolverStrategy::LevenbergMarquardt) {
        int minIter = (n > 10) ? 500 : 300;
        if (solverOpts.maxIterations < minIter) solverOpts.maxIterations = minIter;
    }
    if (strategy == SolverStrategy::Homotopy) {
        // Homotopy uses internal Newton corrector steps; allow generous budget
        int minIter = 200;
        if (solverOpts.maxIterations < minIter) solverOpts.maxIterations = minIter;
    }
    if (strategy == SolverStrategy::Kinsol) {
        // KINSOL (esp. Picard / fixed-point modes) can be slow to converge;
        // allow a generous budget like TrustRegion/LM.
        int minIter = (n > 10) ? 500 : 300;
        if (solverOpts.maxIterations < minIter) solverOpts.maxIterations = minIter;
    }

    try {
        return solver->solve(problem, x, solverOpts, trace, outErrorMessage);
    } catch (const std::exception& e) {
        if (outErrorMessage) {
            *outErrorMessage = e.what();
        }
        return SolverStatus::EvaluationError;
    }
}

// ----------------------------------------------------------------------------
// solveBlockSequential – fallback chain
// ----------------------------------------------------------------------------

SolverStatus Solver::solveBlockSequential(size_t blockIndex,
                                          NonLinearSolver::Problem& problem,
                                          BlockEvaluator& blockEval,
                                          const std::vector<std::string>& varNames,
                                          const std::map<std::string, double>& externalVars,
                                          const std::map<std::string, std::string>& externalStringVars,
                                          Eigen::VectorXd& x,
                                          const SolverOptions& options,
                                          SolverTrace* trace,
                                          std::string* outErrorMessage,
                                          const Eigen::VectorXd* warmStartGuess) {
    const size_t n = varNames.size();
    SolverStatus status = SolverStatus::InvalidInput;
    std::string lastError;

    // Track the best (lowest-residual) solution found across all solvers so
    // that subsequent solvers can start from a better point instead of always
    // resetting to the original initial guess.
    Eigen::VectorXd bestX = x;
    double bestResidualNorm = std::numeric_limits<double>::max();
    double initialResidualNorm = std::numeric_limits<double>::max();

    // Evaluate residual at the initial guess for comparison
    {
        Eigen::VectorXd F0(n);
        Eigen::MatrixXd J0;
        try {
            problem.evaluate(x, F0, J0, false);
            bestResidualNorm = F0.lpNorm<Eigen::Infinity>();
            initialResidualNorm = bestResidualNorm;
        } catch (...) {
            // If even the initial eval fails, keep max residual
        }
    }

    // ── Multi-round pipeline ────────────────────────────────────────
    // If a full pass through the pipeline reduced the residual
    // substantially (>50%) but did not converge, restart the pipeline
    // from the best solution found so far.  This lets later-stage
    // solvers (such as Newton) finish the job when the initial guess
    // was originally too far away.
    constexpr int    MAX_ROUNDS             = 10;
    constexpr double RESTART_MIN_IMPROVEMENT = 0.05; // restart if residual dropped by at least 5%
    bool skipPartitioned = false; // set if Partitioned worsens the solution

    for (int round = 0; round < MAX_ROUNDS; ++round) {
        double roundStartResidual = bestResidualNorm;

        // Check timeout before each pipeline round
        if (TimeoutGuard::hasTimedOut()) break;

    for (size_t idx = 0; idx < options.solverPipeline.size(); ++idx) {
        SolverStrategy strategy = options.solverPipeline[idx];

        // Check timeout before each solver attempt
        if (TimeoutGuard::hasTimedOut()) break;

        // Skip Partitioned if it previously worsened the solution
        if (strategy == SolverStrategy::Partitioned && skipPartitioned) {
            continue;
        }

        if (options.verbose) {
            std::cout << "Block " << blockIndex << " (size " << n
                      << "): Trying " << strategyToString(strategy)
                      << " [" << (idx + 1) << "/" << options.solverPipeline.size() << "]"
                      << (round > 0 ? " (round " + std::to_string(round + 1) + ")" : "")
                      << std::endl;
        }

        // Start from whichever is better: initial evaluator state or best
        // solution found so far by a previous solver.
        if (idx > 0 || round > 0) {
            Eigen::VectorXd x_init(n);
            if (warmStartGuess) {
                // Thread-local warm-start context (e.g. a multi-start candidate):
                // chain solvers against THIS starting point rather than the shared
                // evaluator state, which may hold an unrelated guess when several
                // candidates run concurrently.
                x_init = *warmStartGuess;
            } else {
                for (size_t i = 0; i < n; ++i) {
                    x_init[i] = evaluator_.getVariableValue(varNames[i]);
                }
            }
            // Evaluate residual at initial guess
            double initResidualNorm = std::numeric_limits<double>::max();
            try {
                Eigen::VectorXd F_init(n);
                Eigen::MatrixXd J_dummy;
                problem.evaluate(x_init, F_init, J_dummy, false);
                initResidualNorm = F_init.lpNorm<Eigen::Infinity>();
            } catch (...) {}

            // Use whichever starting point has the smaller residual
            if (bestResidualNorm < initResidualNorm) {
                x = bestX;
                if (options.verbose) {
                    std::cout << "Block " << blockIndex
                              << ": Warm-starting from best previous solution"
                              << " (||F|| = " << bestResidualNorm
                              << " vs initial " << initResidualNorm << ")"
                              << std::endl;
                }
            } else {
                x = x_init;
            }
        }

        std::string solverError;
        auto attemptStart = std::chrono::high_resolution_clock::now();
        status = runSolverStrategy(strategy, blockIndex, problem, blockEval,
                                   varNames, externalVars, externalStringVars,
                                   x, options, trace, &solverError);
        auto attemptEnd = std::chrono::high_resolution_clock::now();
        double attemptMs = std::chrono::duration<double, std::milli>(attemptEnd - attemptStart).count();

        // Record this attempt in the trace (if tracing is enabled)
        if (trace) {
            SolverAttempt attempt;
            attempt.strategy = strategy;
            attempt.status = status;
            attempt.elapsedMs = attemptMs;
            // Evaluate residual at current x for the attempt record
            try {
                Eigen::VectorXd F_att(n);
                Eigen::MatrixXd J_att;
                problem.evaluate(x, F_att, J_att, false);
                attempt.finalResidual = F_att.lpNorm<Eigen::Infinity>();
            } catch (...) {
                attempt.finalResidual = std::numeric_limits<double>::max();
            }
            if (trace && !trace->iterations.empty()) {
                attempt.iterations = static_cast<int>(trace->iterations.size());
            }
            trace->solverAttempts.push_back(attempt);
        }

        if (status == SolverStatus::Success) {
            return status;
        }

        // Fatal evaluation error (unsupported function/fluid) — stop pipeline immediately
        if (isFatalEvaluationError(solverError)) {
            if (outErrorMessage) {
                *outErrorMessage = solverError;
            }
            return SolverStatus::EvaluationError;
        }

        // Even on failure, check whether this solver found a better point
        {
            Eigen::VectorXd F_curr(n);
            Eigen::MatrixXd J_dummy;
            try {
                problem.evaluate(x, F_curr, J_dummy, false);
                double currResidualNorm = F_curr.lpNorm<Eigen::Infinity>();
                if (currResidualNorm < bestResidualNorm) {
                    bestResidualNorm = currResidualNorm;
                    bestX = x;
                } else if (strategy == SolverStrategy::Partitioned &&
                           currResidualNorm > bestResidualNorm * 10.0) {
                    // Partitioned solver worsened solution dramatically — skip it in future rounds
                    skipPartitioned = true;
                }
            } catch (...) {
                if (strategy == SolverStrategy::Partitioned) {
                    skipPartitioned = true;
                }
                // Current x is infeasible, ignore it
            }
        }

        // Accumulate error messages
        if (!solverError.empty()) {
            if (!lastError.empty()) lastError += "\n";
            lastError += "[" + strategyToString(strategy) + "] " + solverError;
        }

        if (options.verbose) {
            std::cout << "Block " << blockIndex << ": "
                      << strategyToString(strategy) << " failed ("
                      << statusToString(status) << ")" << std::endl;
        }
    } // end of pipeline loop

    // Check if this round made sufficient progress to justify restarting
    if (bestResidualNorm < roundStartResidual * (1.0 - RESTART_MIN_IMPROVEMENT) &&
        bestResidualNorm > options.tolerance) {
        if (options.verbose) {
            std::cout << "Block " << blockIndex
                      << ": Pipeline round " << (round + 1)
                      << " reduced residual from " << roundStartResidual
                      << " to " << bestResidualNorm
                      << " — restarting pipeline" << std::endl;
        }
        // Clear accumulated errors for fresh reporting on retry
        lastError.clear();
        continue; // restart the pipeline
    }
    break; // no significant progress, stop
    } // end of round loop

    // Restore best solution found across all solvers and rounds
    x = bestX;

    // All solvers failed — build informative error message
    if (outErrorMessage) {
        std::ostringstream ss;
        if (!lastError.empty()) {
            ss << lastError;
        }
        // Add summary of progress made
        if (initialResidualNorm < std::numeric_limits<double>::max()) {
            ss << "\nInitial ||F||_inf = " << initialResidualNorm
               << ", best achieved = " << bestResidualNorm;
            if (bestResidualNorm < initialResidualNorm * 0.9) {
                double reduction = (1.0 - bestResidualNorm / initialResidualNorm) * 100.0;
                if (reduction > 99.99) {
                    ss << " (>99.99% reduction)";
                } else {
                    ss << " (" << std::fixed << std::setprecision(1)
                       << reduction << "% reduction)";
                }
            }
        }
        if (!outErrorMessage->empty()) *outErrorMessage += "\n";
        *outErrorMessage += ss.str();
    }
    return status;
}

// ----------------------------------------------------------------------------
// solveBlockParallel – concurrent first-to-solve-wins
// ----------------------------------------------------------------------------

SolverStatus Solver::solveBlockParallel(size_t blockIndex,
                                        NonLinearSolver::Problem& problem,
                                        BlockEvaluator& blockEval,
                                        const std::vector<std::string>& varNames,
                                        const std::map<std::string, double>& externalVars,
                                        const std::map<std::string, std::string>& externalStringVars,
                                        Eigen::VectorXd& x,
                                        const SolverOptions& options,
                                        SolverTrace* trace,
                                        std::string* outErrorMessage) {
    const size_t n = varNames.size();
    const auto& pipeline = options.solverPipeline;

    // Each thread gets its own copy of x and its own Problem lambda.
    // The Problem lambda captures blockEval by reference, but BlockEvaluator::evaluate
    // is const-safe (it doesn't mutate shared state), so concurrent calls are safe
    // as long as each thread uses its own x vector.

    struct ThreadResult {
        SolverStatus status = SolverStatus::InvalidInput;
        Eigen::VectorXd solution;
        SolverTrace trace;
        std::string error;
        SolverStrategy strategy;
    };

    std::vector<std::future<ThreadResult>> futures;
    futures.reserve(pipeline.size());

    // Shared stop flag: fires when any thread wins OR user cancels
    auto parallelStop = std::make_shared<std::atomic<bool>>(false);

    // Monitor thread: propagate external cancel → parallelStop
    std::thread cancelMonitor;
    if (options.cancelToken) {
        cancelMonitor = std::thread([cancelToken = options.cancelToken, parallelStop]() {
            while (!parallelStop->load(std::memory_order_relaxed)) {
                if (cancelToken->load(std::memory_order_relaxed)) {
                    parallelStop->store(true, std::memory_order_release);
                    return;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        });
    }

    for (const auto& strategy : pipeline) {
        // Each thread gets its own copy of x
        Eigen::VectorXd x_copy = x;

        // Each thread needs its own Problem with its own capture of x_copy
        // We create a new problem lambda per thread
        futures.push_back(std::async(std::launch::async,
            [this, strategy, blockIndex, &blockEval, &varNames,
             &externalVars, &externalStringVars, &options,
             x_copy, parallelStop]() mutable -> ThreadResult {
                ThreadResult result;
                result.strategy = strategy;
                result.solution = x_copy;

                // Thread-local options with combined cancel token
                SolverOptions threadOptions = options;
                threadOptions.cancelToken = parallelStop.get();

                // Create a thread-local problem lambda
                NonLinearSolver::Problem localProblem;
                localProblem.size = static_cast<int>(varNames.size());
                localProblem.evaluate = [&blockEval, &varNames, &externalVars, &externalStringVars]
                                        (const Eigen::VectorXd& xv,
                                         Eigen::VectorXd& F,
                                         Eigen::MatrixXd& J,
                                         bool computeJacobian) {
                    std::vector<double> x_std(xv.data(), xv.data() + xv.size());
                    auto evalResult = blockEval.evaluate(x_std, externalVars, externalStringVars, computeJacobian);
                    const size_t nEqs = evalResult.residuals.size();
                    F.resize(nEqs);
                    for (size_t i = 0; i < nEqs; ++i) F[i] = evalResult.residuals[i];
                    if (computeJacobian) {
                        J.resize(nEqs, xv.size());
                        for (size_t i = 0; i < nEqs; ++i)
                            for (size_t j = 0; j < evalResult.jacobian[i].size(); ++j)
                                J(i, j) = evalResult.jacobian[i][j];
                    }
                };

                result.status = runSolverStrategy(strategy, blockIndex, localProblem,
                                                  blockEval, varNames,
                                                  externalVars, externalStringVars,
                                                  result.solution, threadOptions,
                                                  &result.trace, &result.error);

                if (result.status == SolverStatus::Success) {
                    parallelStop->store(true, std::memory_order_release);
                }
                return result;
            }));
    }

    // Collect results – pick the first successful one
    SolverStatus bestStatus = SolverStatus::InvalidInput;
    std::string allErrors;

    for (auto& fut : futures) {
        ThreadResult result = fut.get();

        if (result.status == SolverStatus::Success && bestStatus != SolverStatus::Success) {
            bestStatus = SolverStatus::Success;
            x = result.solution;
            if (trace) *trace = result.trace;
            if (options.verbose) {
                std::cout << "Block " << blockIndex << ": "
                          << strategyToString(result.strategy)
                          << " won (parallel)" << std::endl;
            }
        } else if (bestStatus != SolverStatus::Success) {
            // Keep the "best" failure (prefer MaxIterations over others)
            if (result.status == SolverStatus::MaxIterations ||
                bestStatus == SolverStatus::InvalidInput) {
                bestStatus = result.status;
                x = result.solution;
                if (trace) *trace = result.trace;
            }
        }

        if (!result.error.empty()) {
            if (!allErrors.empty()) allErrors += "\n";
            allErrors += "[" + strategyToString(result.strategy) + "] " + result.error;
        }
    }

    if (bestStatus != SolverStatus::Success && outErrorMessage && !allErrors.empty()) {
        if (!outErrorMessage->empty()) *outErrorMessage += "\n";
        *outErrorMessage += allErrors;
    }

    // Shut down monitor thread
    parallelStop->store(true, std::memory_order_release);
    if (cancelMonitor.joinable()) cancelMonitor.join();

    return bestStatus;
}


SolverStatus Solver::solveBlockPartitioned(size_t blockIndex,
                                           BlockEvaluator& blockEval,
                                           const std::vector<std::string>& varNames,
                                           const std::map<std::string, double>& externalVars,
                                           const std::map<std::string, std::string>& externalStringVars,
                                           Eigen::VectorXd& x,
                                           const SolverOptions& options,
                                           SolverTrace* trace,
                                           std::string* outErrorMessage) {
    auto startTime = std::chrono::high_resolution_clock::now();

    if (trace) {
        if (trace->solverType.empty()) {
            trace->solverType = "Partitioned";
        } else if (trace->solverType.find("Partitioned") == std::string::npos) {
            trace->solverType += " -> Partitioned";
        }
    }

    const size_t n = varNames.size();
    const auto& equationIds = blockEval.getEquationIds();
    if (n == 0 || equationIds.size() != n) {
        if (outErrorMessage) {
            std::ostringstream ss;
            ss << "Partitioned solver requires a square block (vars=" << n
               << ", eqs=" << equationIds.size() << ").";
            *outErrorMessage = ss.str();
        }
        return SolverStatus::InvalidInput;
    }

    // Map variable name -> index (case-insensitive)
    std::map<std::string, size_t, CaseInsensitiveLess> varIndex;
    for (size_t i = 0; i < n; ++i) {
        varIndex[varNames[i]] = i;
    }

    const auto& equations = ir_.getEquations();
    std::vector<int> eqToVarIndex(equationIds.size(), -1);
    std::vector<bool> varUsed(n, false);

    for (size_t eq = 0; eq < equationIds.size(); ++eq) {
        int eqId = equationIds[eq];
        if (eqId < 0 || eqId >= static_cast<int>(equations.size())) {
            continue;
        }
        const auto& eqInfo = equations[eqId];
        if (eqInfo.outputVariable.empty()) {
            continue;
        }
        auto it = varIndex.find(eqInfo.outputVariable);
        if (it == varIndex.end()) {
            continue;
        }
        size_t varIdx = it->second;
        if (varUsed[varIdx]) {
            continue;
        }
        eqToVarIndex[eq] = static_cast<int>(varIdx);
        varUsed[varIdx] = true;
    }

    for (size_t eq = 0; eq < eqToVarIndex.size(); ++eq) {
        if (eqToVarIndex[eq] < 0) {
            if (outErrorMessage) {
                std::ostringstream ss;
                ss << "Partitioned solver missing output-variable mapping for equation "
                   << equationIds[eq] << ".";
                *outErrorMessage = ss.str();
            }
            return SolverStatus::InvalidInput;
        }
    }

    double initialResidualNorm = 0.0;

    for (int iter = 0; iter < options.partitionedMaxIterations; ++iter) {
        if (options.cancelToken && options.cancelToken->load(std::memory_order_relaxed)) {
            if (outErrorMessage) *outErrorMessage = "Partitioned: cancelled";
            return SolverStatus::MaxIterations;
        }
        EvaluationResult evalResult;
        try {
            std::vector<double> x_std(x.data(), x.data() + x.size());
            evalResult = blockEval.evaluate(x_std, externalVars, externalStringVars);
        } catch (const std::exception& e) {
            if (outErrorMessage) {
                *outErrorMessage = std::string("Partitioned solver evaluation failed: ") + e.what();
            }
            if (trace) {
                trace->finalStatus = SolverStatus::EvaluationError;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::EvaluationError;
        }

        Eigen::VectorXd F(evalResult.residuals.size());
        for (size_t i = 0; i < evalResult.residuals.size(); ++i) {
            F[static_cast<int>(i)] = evalResult.residuals[i];
        }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) {
            initialResidualNorm = residualNorm;
        }

        if (options.verbose) {
            std::cout << "Partitioned iter " << iter << ": ||F||_inf = " << residualNorm << std::endl;
        }

        if (trace) {
            SolverTrace::Iteration traceIter;
            traceIter.iter = iter;
            traceIter.residualNorm = residualNorm;
            traceIter.stepNorm = 0.0;
            traceIter.lambda = options.partitionedRelaxation;
            traceIter.x = std::vector<double>(x.data(), x.data() + x.size());
            traceIter.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(traceIter);
        }

        if (residualNorm < options.tolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::Success;
        }

        if (initialResidualNorm > 0 && residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::Success;
        }

        Eigen::VectorXd dx = Eigen::VectorXd::Zero(static_cast<int>(n));
        for (size_t eq = 0; eq < eqToVarIndex.size(); ++eq) {
            int varIdx = eqToVarIndex[eq];
            double diag = evalResult.jacobian[eq][static_cast<size_t>(varIdx)];
            if (std::abs(diag) < options.partitionedMinDiagonal) {
                continue;
            }
            double step = -options.partitionedRelaxation * evalResult.residuals[eq] / diag;
            if (!std::isfinite(step)) {
                continue;
            }
            dx[varIdx] = step;
        }

        double stepNorm = dx.lpNorm<Eigen::Infinity>();
        x += dx;

        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = stepNorm;
        }

        if (stepNorm < options.stepTolerance) {
            if (residualNorm < options.tolerance * 100 || residualNorm < options.lsRelaxedTolerance) {
                if (trace) {
                    trace->finalStatus = SolverStatus::Success;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                return SolverStatus::Success;
            }
        }
    }

    if (trace) {
        trace->finalStatus = SolverStatus::MaxIterations;
        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
    }
    if (outErrorMessage) {
        std::ostringstream ss;
        ss << "Partitioned solver: Max iterations (" << options.partitionedMaxIterations
           << ") reached without convergence.";
        *outErrorMessage = ss.str();
    }
    return SolverStatus::MaxIterations;
}

// ----------------------------------------------------------------------------
// solveBlockTearing – structural tearing (FVS + acyclic solve + Newton on tears)
// ----------------------------------------------------------------------------

SolverStatus Solver::solveBlockTearing(size_t blockIndex,
                                       BlockEvaluator& blockEval,
                                       const std::vector<std::string>& varNames,
                                       const std::map<std::string, double>& externalVars,
                                       const std::map<std::string, std::string>& externalStringVars,
                                       Eigen::VectorXd& x,
                                       const SolverOptions& options,
                                       SolverTrace* trace,
                                       std::string* outErrorMessage) {
    auto startTime = std::chrono::high_resolution_clock::now();
    if (trace) {
        if (trace->solverType.empty()) {
            trace->solverType = "Tearing";
        } else if (trace->solverType.find("Tearing") == std::string::npos) {
            trace->solverType += " -> Tearing";
        }
        trace->varNames = varNames;
    }

    const size_t n = varNames.size();
    if (blockIndex >= analysis_.blocks.size()) {
        if (outErrorMessage) *outErrorMessage = "Invalid block index for tearing";
        return SolverStatus::InvalidInput;
    }
    const Block& block = analysis_.blocks[blockIndex];
    BlockTearSetResult tearResult = computeBlockTearSet(block, ir_);
    if (tearResult.tearVarNames.empty()) {
        if (outErrorMessage) *outErrorMessage = "Tear set empty (block too small or no cycles)";
        return SolverStatus::InvalidInput;
    }

    // Store tear variable names in trace
    if (trace) {
        trace->tearVarNames = tearResult.tearVarNames;
    }

    auto caseInsensitiveEqual = [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (std::tolower(static_cast<unsigned char>(a[i])) !=
                std::tolower(static_cast<unsigned char>(b[i]))) return false;
        }
        return true;
    };

    std::vector<size_t> tearVarIndices;
    for (const auto& tv : tearResult.tearVarNames) {
        for (size_t i = 0; i < n; ++i) {
            if (caseInsensitiveEqual(varNames[i], tv)) {
                tearVarIndices.push_back(i);
                break;
            }
        }
    }
    const auto& equationIds = blockEval.getEquationIds();
    std::vector<size_t> tearEqLocalIndices;
    for (int globalEqId : tearResult.tearEquationIds) {
        for (size_t eq = 0; eq < equationIds.size(); ++eq) {
            if (equationIds[eq] == globalEqId) {
                tearEqLocalIndices.push_back(eq);
                break;
            }
        }
    }
    if (tearEqLocalIndices.size() != tearVarIndices.size()) {
        if (outErrorMessage) *outErrorMessage = "Tear set size mismatch";
        return SolverStatus::InvalidInput;
    }

    std::vector<size_t> nonTearEqLocalIndices;
    for (int globalEqId : tearResult.topoOrderNonTearEqIds) {
        for (size_t eq = 0; eq < equationIds.size(); ++eq) {
            if (equationIds[eq] == globalEqId) {
                nonTearEqLocalIndices.push_back(eq);
                break;
            }
        }
    }

    std::map<std::string, size_t, CaseInsensitiveLess> varNameToIndex;
    for (size_t i = 0; i < n; ++i) varNameToIndex[varNames[i]] = i;

    const size_t nTear = tearVarIndices.size();
    const int maxOuter = options.tearingMaxIterations;
    const int maxInner = options.tearingInnerIterations;

    if (options.verbose) {
        std::cout << "Tearing: block " << blockIndex << " tear vars=" << nTear
                  << " acyclic eqs=" << nonTearEqLocalIndices.size() << std::endl;
        for (size_t i = 0; i < tearResult.tearVarNames.size(); ++i) {
            std::cout << "  tear[" << i << "] = " << tearResult.tearVarNames[i]
                      << " (x[" << tearVarIndices[i] << "])" << std::endl;
        }
    }

    for (int outer = 0; outer < maxOuter; ++outer) {
        // Check timeout at each outer iteration
        if (TimeoutGuard::hasTimedOut()) {
            if (outErrorMessage) *outErrorMessage = "Tearing: TIMEOUT";
            if (trace) {
                trace->finalStatus = SolverStatus::MaxIterations;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::MaxIterations;
        }

        std::vector<double> x_std(x.data(), x.data() + x.size());
        EvaluationResult evalResult;
        try {
            evalResult = blockEval.evaluate(x_std, externalVars, externalStringVars);
        } catch (const std::exception& e) {
            if (outErrorMessage) *outErrorMessage = std::string("Tearing evaluation failed: ") + e.what();
            if (trace) {
                trace->finalStatus = SolverStatus::EvaluationError;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::EvaluationError;
        }

        Eigen::VectorXd F(static_cast<int>(evalResult.residuals.size()));
        Eigen::MatrixXd J(F.size(), static_cast<int>(n));
        for (size_t i = 0; i < evalResult.residuals.size(); ++i) F(static_cast<int>(i)) = evalResult.residuals[i];
        for (size_t i = 0; i < evalResult.jacobian.size(); ++i) {
            for (size_t j = 0; j < evalResult.jacobian[i].size(); ++j) {
                J(static_cast<int>(i), static_cast<int>(j)) = evalResult.jacobian[i][j];
            }
        }

        // Acyclic forward sweep: use the Jacobian from the initial evaluation
        // for a linearized forward substitution.  Each acyclic equation's
        // matched variable is updated with a single Newton step using the
        // diagonal Jacobian entry.  After each update, downstream residuals
        // are corrected via the Jacobian (off-diagonal entries), avoiding
        // re-evaluation of the entire block.  This replaces the previous
        // inner Newton loop that called blockEval.evaluate() for every
        // acyclic equation × up to maxInner iterations — the dominant cost
        // for blocks with expensive CoolProp or procedure-call evaluations.
        for (size_t k = 0; k < nonTearEqLocalIndices.size(); ++k) {
            size_t eq = nonTearEqLocalIndices[k];
            double fEq = F(static_cast<int>(eq));
            double jEq = J(static_cast<int>(eq), static_cast<int>(eq));
            if (std::abs(jEq) < 1e-14) continue;
            double step = -fEq / jEq;
            if (!std::isfinite(step)) continue;
            x[eq] += step;

            // Propagate update to downstream acyclic equations (linearized).
            // For an acyclic (lower-triangular) system this is exact for
            // linear equations and first-order accurate for nonlinear ones.
            for (size_t m = k + 1; m < nonTearEqLocalIndices.size(); ++m) {
                size_t nextEq = nonTearEqLocalIndices[m];
                F(static_cast<int>(nextEq)) += J(static_cast<int>(nextEq), static_cast<int>(eq)) * step;
            }
        }

        x_std.assign(x.data(), x.data() + x.size());
        try {
            evalResult = blockEval.evaluate(x_std, externalVars, externalStringVars);
        } catch (const std::exception& e) {
            if (outErrorMessage) *outErrorMessage = std::string("Tearing evaluation failed: ") + e.what();
            if (trace) {
                trace->finalStatus = SolverStatus::EvaluationError;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::EvaluationError;
        }
        for (size_t i = 0; i < evalResult.residuals.size(); ++i) F(static_cast<int>(i)) = evalResult.residuals[i];
        for (size_t i = 0; i < evalResult.jacobian.size(); ++i) {
            for (size_t j = 0; j < evalResult.jacobian[i].size(); ++j) {
                J(static_cast<int>(i), static_cast<int>(j)) = evalResult.jacobian[i][j];
            }
        }

        Eigen::VectorXd F_tear(static_cast<int>(nTear));
        for (size_t i = 0; i < nTear; ++i) {
            F_tear(static_cast<int>(i)) = F(static_cast<int>(tearEqLocalIndices[i]));
        }

        // Schur complement Jacobian: S = A - B D^{-1} C (total derivative of tear residuals w.r.t. tear vars)
        // A = J_tear,tear, B = J_tear,acyclic, C = J_acyclic,tear, D = J_acyclic,acyclic (lower triangular)
        const size_t nAcyclic = nonTearEqLocalIndices.size();
        Eigen::MatrixXd S(static_cast<int>(nTear), static_cast<int>(nTear));
        if (nAcyclic == 0) {
            for (size_t i = 0; i < nTear; ++i) {
                for (size_t j = 0; j < nTear; ++j) {
                    S(static_cast<int>(i), static_cast<int>(j)) =
                        J(static_cast<int>(tearEqLocalIndices[i]), static_cast<int>(tearVarIndices[j]));
                }
            }
        } else {
            Eigen::MatrixXd A(static_cast<int>(nTear), static_cast<int>(nTear));
            Eigen::MatrixXd B(static_cast<int>(nTear), static_cast<int>(nAcyclic));
            Eigen::MatrixXd C(static_cast<int>(nAcyclic), static_cast<int>(nTear));
            Eigen::MatrixXd D(static_cast<int>(nAcyclic), static_cast<int>(nAcyclic));
            for (size_t i = 0; i < nTear; ++i) {
                for (size_t j = 0; j < nTear; ++j) {
                    A(static_cast<int>(i), static_cast<int>(j)) =
                        J(static_cast<int>(tearEqLocalIndices[i]), static_cast<int>(tearVarIndices[j]));
                }
                for (size_t j = 0; j < nAcyclic; ++j) {
                    size_t acyclicVarIdx = nonTearEqLocalIndices[j];
                    B(static_cast<int>(i), static_cast<int>(j)) =
                        J(static_cast<int>(tearEqLocalIndices[i]), static_cast<int>(acyclicVarIdx));
                }
            }
            for (size_t i = 0; i < nAcyclic; ++i) {
                for (size_t j = 0; j < nTear; ++j) {
                    C(static_cast<int>(i), static_cast<int>(j)) =
                        J(static_cast<int>(nonTearEqLocalIndices[i]), static_cast<int>(tearVarIndices[j]));
                }
                for (size_t j = 0; j < nAcyclic; ++j) {
                    D(static_cast<int>(i), static_cast<int>(j)) =
                        J(static_cast<int>(nonTearEqLocalIndices[i]), static_cast<int>(nonTearEqLocalIndices[j]));
                }
            }
            // V = D^{-1} C by forward substitution (D is lower triangular)
            Eigen::MatrixXd V = D.triangularView<Eigen::Lower>().solve(C);
            S = A - B * V;
        }

        double tearNorm = F_tear.lpNorm<Eigen::Infinity>();
        double fullNorm = F.lpNorm<Eigen::Infinity>();
        if (trace) {
            SolverTrace::Iteration traceIter;
            traceIter.iter = outer;
            traceIter.residualNorm = tearNorm;
            traceIter.stepNorm = 0.0;
            traceIter.lambda = 1.0;
            traceIter.x = std::vector<double>(x.data(), x.data() + x.size());
            traceIter.residuals = std::vector<double>(F_tear.data(), F_tear.data() + F_tear.size());
            // Build detail string with tear variable values and full residual norm
            std::ostringstream detail;
            detail << std::scientific << std::setprecision(6);
            detail << "       Full ||F||_inf = " << fullNorm << "\n";
            detail << "       Tear vars: ";
            for (size_t i = 0; i < nTear; ++i) {
                if (i > 0) detail << ", ";
                detail << tearResult.tearVarNames[i] << "=" << x[tearVarIndices[i]];
            }
            detail << "\n";
            detail << "       Tear residuals: ";
            for (size_t i = 0; i < nTear; ++i) {
                if (i > 0) detail << ", ";
                detail << F_tear(static_cast<int>(i));
            }
            detail << "\n";
            traceIter.detail = detail.str();
            trace->iterations.push_back(traceIter);
        }
        if (options.verbose) {
            std::cout << "Tearing outer " << outer << ": ||F_tear||_inf = " << tearNorm << std::endl;
        }

        // Use the full-system residual norm (including acyclic equations)
        // to declare convergence.  Checking only tear residuals can hide
        // failures in the acyclic forward sweep (e.g. silently skipped
        // equations due to evaluation errors).
        if (fullNorm < options.tolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::Success;
        }

        Eigen::VectorXd dx_tear;
        Eigen::FullPivLU<Eigen::MatrixXd> lu(S);
        if (!lu.isInvertible()) {
            // Schur complement is singular — tear Newton step impossible.
            // However, if the acyclic forward sweep already brought the
            // full residual close to tolerance, accept the current x.
            if (fullNorm < std::max(options.tolerance * 1000.0, 1e-3)) {
                if (trace) {
                    trace->finalStatus = SolverStatus::Success;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                return SolverStatus::Success;
            }
            if (outErrorMessage) *outErrorMessage = "Tearing: singular Schur complement (tear system)";
            if (trace) {
                trace->finalStatus = SolverStatus::SingularJacobian;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::SingularJacobian;
        }
        dx_tear = lu.solve(-F_tear);
        for (size_t i = 0; i < nTear; ++i) {
            x[tearVarIndices[i]] += dx_tear(static_cast<int>(i));
        }
        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = dx_tear.lpNorm<Eigen::Infinity>();
        }
    }

    if (trace) {
        trace->finalStatus = SolverStatus::MaxIterations;
        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
    }
    if (outErrorMessage) {
        std::ostringstream ss;
        ss << "Tearing: max iterations (" << options.tearingMaxIterations << ") reached.";
        *outErrorMessage = ss.str();
    }
    return SolverStatus::MaxIterations;
}

// ============================================================================
// Multi-start fallback (roadmap §4.2)
// ============================================================================
//
// When a multi-variable block fails the entire solver pipeline, the dominant
// cause (for models without a .initials file) is that the default guess lies
// in the wrong convergence basin — typically the wrong *pressure level* for a
// refrigeration / power cycle, since initializeVariables() guesses every
// thermo variable at T=20°C, P=1 atm regardless of the cycle's operating
// point.  The candidates below shift every thermo variable *together* to a
// plausible operating level, which is far more effective than perturbing
// variables independently (the physical states in a cycle are correlated).
//
// References:
//   - Roadmap §4.2: "Multi-Start for the Without-Initials Case".
//   - Allgower & Georg (2003), §1.3 on the role of starting points for
//     continuation / Newton methods.

namespace {

// A physically-coherent operating regime.  Every thermo variable in the block
// is re-evaluated via CoolProp at (T_K, P_Pa) so the candidate is internally
// consistent (h matches enthalpy(fluid, T, P), etc.) — this is what makes the
// candidates land in real convergence basins rather than producing the huge
// residuals of an inconsistent state.  genericScale perturbs the remaining
// (non-thermo) variables to escape scale locks.
struct MultiStartRegime {
    double T_K;          // reference temperature (Kelvin)
    double P_Pa;         // reference pressure (Pascal)
    double quality;      // for quality variables (no CoolProp recompute)
    double genericScale; // multiplier for non-thermo variables
    const char* label;
};

// Ordered by estimated likelihood of rescuing thermodynamic cycles.  The first
// few cover the common low/medium/high pressure bands of refrigeration and
// ORC loops; the last two add mild generic perturbation to escape scale locks.
const std::vector<MultiStartRegime>& multistartRegimes() {
    static const std::vector<MultiStartRegime> regimes = {
        {353.15,  2.0e6, 0.5, 1.0, "medium pressure (2 MPa, 80C)"},
        {423.15,  1.0e7, 0.5, 1.0, "high pressure (10 MPa, 150C)"},
        {298.15,  5.0e5, 0.5, 1.0, "low pressure (0.5 MPa, 25C)"},
        {373.15,  5.0e6, 0.5, 1.0, "mid-high pressure (5 MPa, 100C)"},
        {273.15,  1.0e5, 0.5, 0.7, "atmospheric (0C, downscale)"},
        {523.15,  2.0e7, 0.5, 1.5, "very high (20 MPa, 250C, upscale)"},
    };
    return regimes;
}

// Global scale factors used for purely-algebraic blocks (no thermo variables).
// Many such blocks solve to small values (efficiencies, volumetric ratios,
// clearance C ~ 0.05, swept volume ~ 0.1) while initializeVariables() guesses
// 1.0, so a downscale is the most effective single move.  The set spans three
// orders of magnitude symmetrically around 1.
const std::vector<double>& multistartScales() {
    static const std::vector<double> s = {0.1, 0.3, 3.0, 10.0};
    return s;
}

} // namespace

std::vector<std::pair<Eigen::VectorXd, std::string>>
Solver::generateMultiStartCandidates(size_t blockIndex, int maxRestarts) const {
    std::vector<std::pair<Eigen::VectorXd, std::string>> candidates;
    if (maxRestarts <= 0) return candidates;

    const auto& varNames = evaluator_.getBlock(blockIndex).getVariables();
    const size_t n = varNames.size();
    if (n <= 1) return candidates;

    // Default guess vector (already tried before multi-start; reused as the
    // base for non-thermo variables).
    Eigen::VectorXd defaultGuess(n);
    bool hasThermo = false;
    for (size_t i = 0; i < n; ++i) {
        defaultGuess[i] = evaluator_.getVariableValue(varNames[i]);
        const VariableInfo* info = ir_.getVariable(varNames[i]);
        if (info && !info->inferredFluid.empty() && !info->inferredProperty.empty()) {
            hasThermo = true;
        }
    }

    // Helper to guard against zero guesses (would collapse a Jacobian column).
    auto guarded = [](double v) { return (v == 0.0) ? 1.0 : v; };

    // Build the candidate list.  Two complementary strategies:
    //   - Thermo regime candidates for blocks containing CoolProp-dependent
    //     variables (recompute them consistently at a reference state).
    //   - Global scale candidates for purely-algebraic blocks, where the
    //     failure is a wrong overall magnitude rather than a wrong fluid state.
    // For thermo blocks the regimes dominate and are tried first; for algebraic
    // blocks the scale factors are the only useful move.
    if (hasThermo) {
        const auto& regimes = multistartRegimes();
        int count = std::min<int>(maxRestarts, static_cast<int>(regimes.size()));
        candidates.reserve(count);
        for (int k = 0; k < count; ++k) {
            const auto& r = regimes[k];
            Eigen::VectorXd cand(n);
            for (size_t i = 0; i < n; ++i) {
                const VariableInfo* info = ir_.getVariable(varNames[i]);
                std::string prop = info ? info->inferredProperty : std::string{};
                std::string fluid = info ? info->inferredFluid : std::string{};
                std::string units = info ? info->units : std::string{};

                double val = std::numeric_limits<double>::quiet_NaN();
                if (!fluid.empty() && !prop.empty()) {
                    if (prop == "Q") {
                        val = r.quality;
                    } else {
                        // Recompute the property consistently at the regime state.
                        auto guessed = computeThermoGuessAt(fluid, prop, r.T_K, r.P_Pa, units);
                        if (guessed && std::isfinite(*guessed)) val = *guessed;
                    }
                }
                if (std::isnan(val)) {
                    val = guarded(defaultGuess[i]) * r.genericScale;
                }
                cand[i] = val;
            }
            candidates.emplace_back(std::move(cand), std::string(r.label));
        }
    } else {
        const auto& scales = multistartScales();
        int count = std::min<int>(maxRestarts, static_cast<int>(scales.size()));
        candidates.reserve(count);
        for (int k = 0; k < count; ++k) {
            Eigen::VectorXd cand(n);
            for (size_t i = 0; i < n; ++i) {
                cand[i] = guarded(defaultGuess[i]) * scales[k];
            }
            std::ostringstream lbl;
            lbl << "scale x" << scales[k];
            candidates.emplace_back(std::move(cand), lbl.str());
        }
    }
    return candidates;
}

SolverStatus Solver::solveBlockWithMultiStart(size_t blockIndex,
                                              const SolverOptions& options,
                                              SolverTrace* trace,
                                              std::string* outErrorMessage,
                                              std::string* multistartInfo) {
    if (multistartInfo) multistartInfo->clear();

    // Capture the original starting guess BEFORE the first attempt: solveBlock
    // leaves the evaluator state at its (possibly diverged) iterate after a
    // failure, so reading the "default" guess later would return that garbage.
    // Multi-start candidates must be derived from the genuine initial guess.
    const auto& varNames0 = evaluator_.getBlock(blockIndex).getVariables();
    Eigen::VectorXd originalGuess(varNames0.size());
    for (size_t i = 0; i < varNames0.size(); ++i) {
        originalGuess[i] = evaluator_.getVariableValue(varNames0[i]);
    }

    // 1. Normal attempt with the default starting guess.
    std::string firstError;
    SolverStatus status = solveBlock(blockIndex, options, trace, &firstError);

    // Conditions to engage multi-start:
    //  - option enabled,
    //  - multi-variable block (size-1 has its own multi-probe in Newton1D),
    //  - a non-fatal failure (EvaluationError/InvalidInput/ParseFailed are not
    //    starting-point problems and would just waste time retrying).
    bool fatal = (status == SolverStatus::EvaluationError ||
                  status == SolverStatus::InvalidInput ||
                  status == SolverStatus::ParseFailed);
    if (status == SolverStatus::Success || !options.multiStartEnabled || fatal) {
        if (outErrorMessage) *outErrorMessage = firstError;
        return status;
    }

    const auto& varNames = evaluator_.getBlock(blockIndex).getVariables();
    if (varNames.size() <= 1) {
        if (outErrorMessage) *outErrorMessage = firstError;
        return status;
    }

    // Restore the genuine initial guess so candidate generation reads clean
    // base values rather than the first attempt's diverged endpoint.
    for (size_t i = 0; i < varNames.size(); ++i) {
        evaluator_.setVariableValue(varNames[i], originalGuess[i]);
    }

    // 2. Generate candidate starting vectors and replay the pipeline.
    auto candidates = generateMultiStartCandidates(
        blockIndex, options.multiStartMaxRestarts);
    if (candidates.empty()) {
        if (outErrorMessage) *outErrorMessage = firstError;
        return status;
    }

    auto candidateResidualNorm = [](SolverTrace* t) -> double {
        if (t && !t->iterations.empty()) return t->iterations.back().residualNorm;
        return std::numeric_limits<double>::infinity();
    };

    // Preserve the default-guess attempt as the initial "best" so that, if no
    // candidate improves things, we restore the state solveBlock left behind.
    double bestResidual = candidateResidualNorm(trace);
    int bestCandidate = -1;  // -1 == default-guess attempt
    Eigen::VectorXd bestState(varNames.size());
    for (size_t i = 0; i < varNames.size(); ++i) {
        bestState[i] = evaluator_.getVariableValue(varNames[i]);
    }

    if (options.verbose) {
        std::cerr << "[MultiStart] Block " << blockIndex << " (size "
                  << varNames.size() << ") failed with default guess ("
                  << statusToString(status) << ", |F|=" << bestResidual
                  << "); trying " << candidates.size() << " alternative start(s).\n";
    }

    // Dispatch: run candidates in parallel when more than one core is
    // available.  The parallel path is self-contained (it updates the
    // evaluator, trace, and multistartInfo itself).  The single-core path
    // below stays inline so the default sequential behaviour is unchanged.
    int effCores = options.multiStartNumCores;
    if (effCores == 0) {
        unsigned hw = std::thread::hardware_concurrency();
        effCores = (hw == 0) ? 1 : static_cast<int>(hw);
    }
    if (effCores > 1 && static_cast<int>(candidates.size()) > 1) {
        return solveBlockMultiStartParallel(blockIndex, candidates, originalGuess,
                                            options, effCores, status, firstError,
                                            trace, outErrorMessage, multistartInfo);
    }

    for (size_t k = 0; k < candidates.size(); ++k) {
        if (TimeoutGuard::hasTimedOut()) break;

        const std::string& candLabel = candidates[k].second;

        // Install candidate starting values.
        for (size_t i = 0; i < varNames.size(); ++i) {
            evaluator_.setVariableValue(varNames[i], candidates[k].first[i]);
        }

        SolverTrace candTrace;
        SolverTrace* pt = trace ? &candTrace : nullptr;
        std::string candError;
        auto t0 = std::chrono::high_resolution_clock::now();
        SolverStatus candStatus = solveBlock(blockIndex, options, pt, &candError);
        auto t1 = std::chrono::high_resolution_clock::now();

        if (options.verbose) {
            double initRes = (pt && !pt->iterations.empty())
                ? pt->iterations.front().residualNorm
                : std::numeric_limits<double>::infinity();
            std::cerr << "[MultiStart]   candidate " << (k + 1) << " ("
                      << candLabel << ") init|F|=" << initRes << "\n";
        }

        if (candStatus == SolverStatus::Success) {
            // Winner: reflect the winning attempt in the reported trace.
            if (trace) {
                trace->iterations = candTrace.iterations;
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = t1 - t0;
                trace->solverAttempts = candTrace.solverAttempts;
                trace->symbolicReductionApplied = candTrace.symbolicReductionApplied;
                trace->originalBlockSize = candTrace.originalBlockSize;
                trace->reducedBlockSize = candTrace.reducedBlockSize;
                trace->symInversions = candTrace.symInversions;
                trace->symExtractions = candTrace.symExtractions;
                trace->symSubstitutions = candTrace.symSubstitutions;
                trace->redecompositionApplied = candTrace.redecompositionApplied;
                trace->numSubBlocks = candTrace.numSubBlocks;
                trace->subBlockSizes = candTrace.subBlockSizes;
                trace->solverType = "MultiStart(" + candLabel + ")->" +
                    (candTrace.solverType.empty() ? std::string("pipeline")
                                                  : candTrace.solverType);
            }
            if (multistartInfo) {
                std::ostringstream ss;
                ss << "block " << blockIndex << " rescued by multi-start candidate "
                   << (k + 1) << "/" << candidates.size()
                   << " (" << candLabel << ")";
                *multistartInfo = ss.str();
            }
            if (options.verbose) {
                std::cerr << "[MultiStart] Block " << blockIndex
                          << " converged with candidate " << (k + 1)
                          << " (" << candLabel << ").\n";
            }
            if (outErrorMessage) outErrorMessage->clear();
            return SolverStatus::Success;
        }

        // Candidate failed too — keep the lowest-residual state for restoration.
        double candResidual = candidateResidualNorm(pt);
        if (candStatus != SolverStatus::EvaluationError && candResidual < bestResidual) {
            bestResidual = candResidual;
            bestCandidate = static_cast<int>(k);
            for (size_t i = 0; i < varNames.size(); ++i) {
                bestState[i] = evaluator_.getVariableValue(varNames[i]);
            }
        }

        if (options.verbose) {
            std::cerr << "[MultiStart]   candidate " << (k + 1) << " ("
                      << candLabel << "): "
                      << statusToString(candStatus) << " |F|=" << candResidual << "\n";
        }
    }

    // 3. No candidate converged — restore the best attempt's state and report
    //    the original failure (the default-guess error message is the most
    //    informative).
    for (size_t i = 0; i < varNames.size(); ++i) {
        evaluator_.setVariableValue(varNames[i], bestState[i]);
    }
    if (options.verbose) {
        std::cerr << "[MultiStart] Block " << blockIndex
                  << " still failed after " << candidates.size()
                  << " restart(s); restored best |F|=" << bestResidual
                  << " (candidate " << (bestCandidate + 1) << ").\n";
    }
    if (outErrorMessage) *outErrorMessage = firstError;
    return status;
}

// ----------------------------------------------------------------------------
// Multi-start parallel execution (roadmap §4.2)
// ----------------------------------------------------------------------------
//
// Concurrency model: identical to solveBlockParallel's.  BlockEvaluator::
// evaluate is const-safe (it does not mutate shared state), and CoolProp's
// AbstractState cache is thread-local, so concurrent candidate solves against
// the same block evaluator are safe as long as each thread uses its own x
// vector and trace.  solveBlockSequential additionally reads evaluator_ once
// per solver for its warm-start heuristic — a concurrent read that is safe
// because no thread writes evaluator_ during the parallel section (the winner
// is written back single-threaded afterwards).
//
// Candidates are run in waves of size `numCores` to honour the configured core
// limit.  A shared atomic stop flag gives first-to-converge semantics: as soon
// as one candidate succeeds, the others' solvers see the cancel token and wind
// down quickly.

namespace {
// One candidate attempt's outcome, produced inside a worker thread.
struct MSTaskResult {
    bool converged = false;
    int index = -1;
    std::string label;
    Eigen::VectorXd solution;
    Eigen::VectorXd finalX;   // solver's final iterate (== solution when converged)
    double residualNorm = std::numeric_limits<double>::infinity();
    SolverTrace trace;
};
} // namespace

SolverStatus Solver::solveBlockMultiStartParallel(
    size_t blockIndex,
    const std::vector<std::pair<Eigen::VectorXd, std::string>>& candidates,
    const Eigen::VectorXd& originalGuess,
    const SolverOptions& options,
    int numCores,
    SolverStatus firstAttemptStatus,
    const std::string& firstAttemptError,
    SolverTrace* trace,
    std::string* outErrorMessage,
    std::string* multistartInfo) {

    const auto& varNames = evaluator_.getBlock(blockIndex).getVariables();
    const size_t n = varNames.size();

    // Build the (shared) evaluation problem, mirroring solveBlock().  The
    // evaluate lambda only captures const-safe references (blockEval,
    // externalVars), so it is safe to invoke concurrently from many threads,
    // each passing its own x.
    BlockEvaluator& blockEval = evaluator_.getBlock(blockIndex);

    auto caseInsensitiveEqual = [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i)
            if (std::tolower(static_cast<unsigned char>(a[i])) !=
                std::tolower(static_cast<unsigned char>(b[i]))) return false;
        return true;
    };
    std::map<std::string, double> externalVars;
    for (const auto& [name, value] : evaluator_.getAllVariables()) {
        bool inBlock = false;
        for (const auto& bvar : varNames)
            if (caseInsensitiveEqual(name, bvar)) { inBlock = true; break; }
        if (!inBlock) externalVars[name] = value;
    }
    std::map<std::string, std::string> externalStringVars;
    for (const auto& [name, value] : evaluator_.getAllStringVariables())
        externalStringVars[name] = value;

    NonLinearSolver::Problem problem;
    problem.size = static_cast<int>(n);
    problem.evaluate = [&blockEval, &externalVars, &externalStringVars](
                           const Eigen::VectorXd& xv, Eigen::VectorXd& F,
                           Eigen::MatrixXd& J, bool computeJacobian) {
        std::vector<double> x_std(xv.data(), xv.data() + xv.size());
        auto result = blockEval.evaluate(x_std, externalVars, externalStringVars,
                                         computeJacobian);
        const size_t nEqs = result.residuals.size();
        F.resize(nEqs);
        for (size_t i = 0; i < nEqs; ++i) F[i] = result.residuals[i];
        if (computeJacobian) {
            J.resize(nEqs, xv.size());
            for (size_t i = 0; i < nEqs; ++i)
                for (size_t j = 0; j < result.jacobian[i].size(); ++j)
                    J(i, j) = result.jacobian[i][j];
        }
    };

    // Keep the evaluator on the genuine initial guess (already restored by the
    // caller) so solveBlockSequential's warm-start read is consistent across
    // all worker threads.
    // (No extra action needed; originalGuess is the current evaluator state.)

    auto parallelStop = std::make_shared<std::atomic<bool>>(false);

    // Lowest-residual fallback across all candidates (for the failure case).
    MSTaskResult best;
    best.index = -1;
    best.label = "default guess";
    best.solution = originalGuess;
    best.finalX = originalGuess;

    const int waveSize = std::max(1, numCores);
    size_t launched = 0;

    while (launched < candidates.size()) {
        if (TimeoutGuard::hasTimedOut()) break;

        // Launch one wave of up to waveSize candidates.
        std::vector<std::future<MSTaskResult>> futures;
        size_t waveEnd = std::min(candidates.size(), launched + waveSize);
        for (size_t k = launched; k < waveEnd; ++k) {
            Eigen::VectorXd x0 = candidates[k].first;        // candidate start
            std::string label = candidates[k].second;
            int idx = static_cast<int>(k);
            futures.push_back(std::async(std::launch::async,
                [this, blockIndex, &varNames, &problem, &blockEval,
                 &externalVars, &externalStringVars, &options, parallelStop,
                 x0 = std::move(x0), label = std::move(label), idx]() mutable {
                    MSTaskResult r;
                    r.index = idx;
                    r.label = label;
                    r.finalX = x0;

                    SolverOptions threadOpts = options;
                    threadOpts.cancelToken = parallelStop.get();

                    Eigen::VectorXd x = x0;
                    SolverTrace tlocal;
                    SolverTrace* pt = &tlocal;
                    std::string err;
                    auto tA = std::chrono::high_resolution_clock::now();
                    // Pass the candidate as the warm-start context so solver
                    // chaining stays within this candidate (the shared
                    // evaluator still holds the original guess).
                    SolverStatus st = solveBlockSequential(
                        blockIndex, problem, blockEval, varNames,
                        externalVars, externalStringVars, x, threadOpts, pt, &err,
                        &x0);
                    auto tB = std::chrono::high_resolution_clock::now();
                    tlocal.totalTime = tB - tA;

                    r.converged = (st == SolverStatus::Success);
                    r.finalX = x;
                    if (r.converged) {
                        r.solution = x;
                        r.residualNorm =
                            (pt && !pt->iterations.empty())
                                ? pt->iterations.back().residualNorm
                                : 0.0;
                        r.trace = std::move(tlocal);
                        parallelStop->store(true, std::memory_order_release);
                    } else {
                        r.residualNorm =
                            (pt && !pt->iterations.empty())
                                ? pt->iterations.back().residualNorm
                                : std::numeric_limits<double>::infinity();
                        r.trace = std::move(tlocal);
                    }
                    return r;
                }));
        }
        launched = waveEnd;

        // Collect the wave, polling so the block timeout / cancel can fire.
        for (auto& fut : futures) {
            if (!fut.valid()) continue;
            while (true) {
                auto st = fut.wait_for(std::chrono::milliseconds(50));
                if (st == std::future_status::ready) break;
                if (TimeoutGuard::hasTimedOut())
                    parallelStop->store(true, std::memory_order_release);
            }
            MSTaskResult r = fut.get();

            if (r.converged) {
                // Winner — write solution back, update trace/info, return.
                for (size_t i = 0; i < n; ++i)
                    evaluator_.setVariableValue(varNames[i], r.solution[i]);
                if (trace) {
                    trace->iterations = r.trace.iterations;
                    trace->finalStatus = SolverStatus::Success;
                    trace->totalTime = r.trace.totalTime;
                    trace->solverAttempts = r.trace.solverAttempts;
                    trace->solverType =
                        "MultiStart(" + r.label + ")->" +
                        (r.trace.solverType.empty() ? std::string("pipeline")
                                                    : r.trace.solverType);
                }
                if (multistartInfo) {
                    std::ostringstream ss;
                    ss << "block " << blockIndex
                       << " rescued by multi-start candidate " << (r.index + 1)
                       << "/" << candidates.size() << " (" << r.label
                       << ") [parallel, " << numCores << " cores]";
                    *multistartInfo = ss.str();
                }
                if (options.verbose) {
                    std::cerr << "[MultiStart] Block " << blockIndex
                              << " converged with candidate " << (r.index + 1)
                              << " (" << r.label << ") [parallel].\n";
                }
                if (outErrorMessage) outErrorMessage->clear();
                return SolverStatus::Success;
            }

            if (r.residualNorm < best.residualNorm) best = std::move(r);
        }

        if (parallelStop->load(std::memory_order_relaxed)) break; // cancelled
    }

    // No candidate converged — restore the lowest-residual iterate and report
    // the original failure.
    for (size_t i = 0; i < n; ++i)
        evaluator_.setVariableValue(varNames[i], best.finalX[i]);
    if (options.verbose) {
        std::cerr << "[MultiStart] Block " << blockIndex
                  << " still failed after " << candidates.size()
                  << " restart(s) [parallel]; restored best |F|=" << best.residualNorm
                  << " (candidate " << (best.index + 1) << ").\n";
    }
    if (outErrorMessage) *outErrorMessage = firstAttemptError;
    return firstAttemptStatus;
}

SolveResult Solver::solve(const SolverOptions& options, bool enableTracing) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    SolveResult result;
    result.blocksEvaluated = 0;
    result.totalIterations = 0;
    
    // Always allocate block traces so iteration counts are tracked even
    // without full debug tracing.  The overhead is minimal.
    result.blockTraces.resize(evaluator_.getNumBlocks());
    
    // Solve blocks in topological order
    int totalBlocks = static_cast<int>(evaluator_.getNumBlocks());
    for (size_t blockIdx = 0; blockIdx < evaluator_.getNumBlocks(); ++blockIdx) {
        // Check cancellation before each block
        if (options.cancelToken && options.cancelToken->load()) {
            result.success = false;
            result.status = SolverStatus::MaxIterations;
            result.errorMessage = "Solve cancelled by user";
            result.variables = evaluator_.getAllVariables();
            result.stringVariables = evaluator_.getAllStringVariables();
            result.totalTime = std::chrono::high_resolution_clock::now() - startTime;
            
            if (options.progressCallback) {
                options.progressCallback(static_cast<int>(blockIdx), totalBlocks, "fail", 0, 0.0);
            }
            return result;
        }
        
        SolverTrace* trace = &result.blockTraces[blockIdx];
        
        // Notify progress: block starting
        if (options.progressCallback) {
            options.progressCallback(static_cast<int>(blockIdx), totalBlocks, "start", 0, 0.0);
        }
        
        // Setup timeout protection
        TimeoutGuard timeout(options.timeoutSeconds);
        
        // Set up per-block diagnostic collection for CoolProp errors
        DiagnosticCollector blockDiag;
        evaluator_.getBlock(blockIdx).setDiagnostics(&blockDiag);
        
        std::string blockError;
        std::string multistartInfo;
        // Multi-start fallback (roadmap §4.2): retries the block from
        // alternative starting points if the normal pipeline fails.
        SolverStatus blockStatus = solveBlockWithMultiStart(
            blockIdx, options, trace, &blockError, &multistartInfo);
        
        // Detach diagnostics pointer (blockDiag goes out of scope later)
        evaluator_.getBlock(blockIdx).setDiagnostics(nullptr);
        
        result.blocksEvaluated++;
        
        // Record per-block result
        SolveResult::BlockResult br;
        br.id = blockIdx;
        br.success = (blockStatus == SolverStatus::Success);
        br.status = blockStatus;
        br.iterations = static_cast<int>(trace->iterations.size());
        br.maxResidual = 0.0;
        if (!trace->iterations.empty()) {
            br.maxResidual = trace->iterations.back().residualNorm;
        }
        br.errorMessage = blockError;

        // Block size info (always recorded)
        const auto& blkVars = evaluator_.getBlock(blockIdx).getVariables();
        br.originalSize = static_cast<int>(blkVars.size());
        br.reducedSize = br.originalSize;

        // Copy symbolic reduction info & solver attempts from trace
        if (trace) {
            if (trace->symbolicReductionApplied) {
                br.symbolicReductionApplied = true;
                br.originalSize = trace->originalBlockSize;
                br.reducedSize = trace->reducedBlockSize;
                br.inversionsApplied = trace->symInversions;
                br.extractionsApplied = trace->symExtractions;
                br.substitutionsApplied = trace->symSubstitutions;
            }
            if (trace->redecompositionApplied) {
                br.redecompositionApplied = true;
                br.numSubBlocks = trace->numSubBlocks;
                br.subBlockSizes = trace->subBlockSizes;
            }
            br.solverAttempts = trace->solverAttempts;
        }

        result.blockResults.push_back(br);
        
        // Emit solver diagnostics for this block
        auto& blkDiag = result.blockResults.back().diagnostics;
        blkDiag.merge(blockDiag);  // Merge CoolProp diagnostics collected during evaluation
        if (!multistartInfo.empty()) {
            // V006: multi-start rescued a block that failed with the default guess
            blkDiag.push(DiagnosticSeverity::Info, "V006",
                "Multi-start: " + multistartInfo, "solver");
        }
        if (br.solverAttempts.size() > 1) {
            // V004: pipeline fallthrough
            std::string attemptList;
            for (size_t i = 0; i < br.solverAttempts.size(); ++i) {
                if (i > 0) attemptList += " -> ";
                attemptList += strategyToString(br.solverAttempts[i].strategy);
                attemptList += "(" + statusToString(br.solverAttempts[i].status) + ")";
            }
            blkDiag.push(DiagnosticSeverity::Info, "V004",
                "Block " + std::to_string(br.id) + ": solver pipeline fallthrough: " + attemptList,
                "solver");
        }
        if (!br.success) {
            // V005: per-block convergence failure
            blkDiag.push(DiagnosticSeverity::Error, "V005",
                "Block " + std::to_string(br.id) + " (" + std::to_string(br.originalSize) +
                " vars) failed to converge: " + statusToString(br.status) +
                " (residual=" + std::to_string(br.maxResidual) + ")",
                "solver");
        }
        
        result.totalIterations += static_cast<int>(trace->iterations.size());
        
        if (blockStatus != SolverStatus::Success) {
            // Notify progress: block failed
            if (options.progressCallback) {
                options.progressCallback(static_cast<int>(blockIdx), totalBlocks, "fail",
                    br.iterations, br.maxResidual);
            }
            
            result.success = false;
            result.status = blockStatus;
            result.detailedError = blockError;
            
            // Get block info for error message
            const auto& block = evaluator_.getBlock(blockIdx);
            std::ostringstream ss;
            ss << "Block " << blockIdx << " (size " << block.size() 
               << ", vars: ";
            const auto& vars = block.getVariables();
            for (size_t i = 0; i < std::min(vars.size(), size_t(3)); ++i) {
                if (i > 0) ss << ", ";
                ss << vars[i];
            }
            if (vars.size() > 3) ss << ", ...";
            ss << ") failed: " << statusToString(blockStatus);
            if (!blockError.empty()) {
                ss << " - " << blockError;
            }
            result.errorMessage = ss.str();
            
            // Copy partial solution anyway
            result.variables = evaluator_.getAllVariables();
            result.stringVariables = evaluator_.getAllStringVariables();
            result.totalTime = std::chrono::high_resolution_clock::now() - startTime;

            // Write debug reduction .md file if requested (even on failure)
            if (!options.debugReductionPath.empty() && enableTracing) {
                writeDebugReductionReport(options.debugReductionPath, result);
            }

            return result;
        }
        
        // Notify progress: block done
        if (options.progressCallback) {
            options.progressCallback(static_cast<int>(blockIdx), totalBlocks, "done",
                br.iterations, br.maxResidual);
        }
    }
    
    // Success - copy final solution
    result.success = true;
    result.status = SolverStatus::Success;
    result.variables = evaluator_.getAllVariables();
    result.stringVariables = evaluator_.getAllStringVariables();
    result.totalTime = std::chrono::high_resolution_clock::now() - startTime;

    // Write debug reduction .md file if requested
    if (!options.debugReductionPath.empty() && enableTracing) {
        writeDebugReductionReport(options.debugReductionPath, result);
    }

    return result;
}

// ============================================================================
// Debug Reduction Report
// ============================================================================

void Solver::writeDebugReductionReport(const std::string& path,
                                       const SolveResult& result) {
    std::ofstream out(path);
    if (!out.is_open()) return;

    out << "# CoolSolve — Symbolic Reduction Debug Report\n\n";
    out << "**Status:** " << (result.success ? "SUCCESS" : "FAILED") << "\n\n";

    // Summary table
    out << "## Block Summary\n\n";
    out << "| Block | Original Size | Reduced Size | Sub-Blocks | Inversions | Extractions | Substitutions | Status |\n";
    out << "|---:|---:|---:|---:|---:|---:|---:|---|\n";

    for (size_t bi = 0; bi < result.blockResults.size(); ++bi) {
        const auto& br = result.blockResults[bi];
        std::string subBlockInfo = "—";
        if (br.redecompositionApplied) {
            subBlockInfo = std::to_string(br.numSubBlocks) + " [";
            for (size_t si = 0; si < br.subBlockSizes.size(); ++si) {
                if (si > 0) subBlockInfo += ",";
                subBlockInfo += std::to_string(br.subBlockSizes[si]);
            }
            subBlockInfo += "]";
        }
        out << "| " << bi
            << " | " << br.originalSize
            << " | " << (br.symbolicReductionApplied ? std::to_string(br.reducedSize) : "—")
            << " | " << subBlockInfo
            << " | " << (br.symbolicReductionApplied ? std::to_string(br.inversionsApplied) : "—")
            << " | " << (br.symbolicReductionApplied ? std::to_string(br.extractionsApplied) : "—")
            << " | " << (br.symbolicReductionApplied ? std::to_string(br.substitutionsApplied) : "—")
            << " | " << (br.success ? "OK" : statusToString(br.status))
            << " |\n";
    }

    // Per-block details for blocks where reduction was applied
    out << "\n## Detailed Reduction Steps\n\n";
    for (size_t bi = 0; bi < result.blockTraces.size(); ++bi) {
        const auto& tr = result.blockTraces[bi];
        if (!tr.symbolicReductionApplied) continue;

        out << "### Block " << bi << " (" << tr.originalBlockSize
            << " → " << tr.reducedBlockSize << " variables)\n\n";

        // Original equations
        out << "**Original equations (" << tr.originalEquations.size() << "):**\n\n";
        for (size_t i = 0; i < tr.originalEquations.size(); ++i) {
            out << (i + 1) << ". `" << tr.originalEquations[i] << "`\n";
        }

        // Reduction steps
        if (!tr.reductionStepDescriptions.empty()) {
            out << "\n**Reduction steps:**\n\n";
            for (size_t i = 0; i < tr.reductionStepDescriptions.size(); ++i) {
                out << (i + 1) << ". " << tr.reductionStepDescriptions[i] << "\n";
            }
        }

        // Remaining equations
        if (!tr.reducedEquations.empty()) {
            out << "\n**Remaining equations (" << tr.reducedEquations.size()
                << ", solved iteratively):**\n\n";
            for (size_t i = 0; i < tr.reducedEquations.size(); ++i) {
                out << (i + 1) << ". `" << tr.reducedEquations[i] << "`\n";
            }
        } else if (tr.reducedBlockSize == 0) {
            out << "\n*Block fully reduced — no iterative solve needed.*\n";
        }

        // Re-decomposition info
        if (bi < result.blockResults.size() && result.blockResults[bi].redecompositionApplied) {
            const auto& br = result.blockResults[bi];
            out << "\n**Re-decomposition:** " << br.numSubBlocks << " sub-blocks [";
            for (size_t si = 0; si < br.subBlockSizes.size(); ++si) {
                if (si > 0) out << ", ";
                out << br.subBlockSizes[si];
            }
            out << "]\n";
        }

        out << "\n";
    }

    // Per-block solver pipeline results
    out << "## Solver Pipeline Results\n\n";
    for (size_t bi = 0; bi < result.blockResults.size(); ++bi) {
        const auto& br = result.blockResults[bi];
        if (br.originalSize <= 1 && br.solverAttempts.empty()) continue;

        out << "### Block " << bi << " (size "
            << (br.symbolicReductionApplied ? std::to_string(br.reducedSize) : std::to_string(br.originalSize))
            << ")\n\n";

        if (br.solverAttempts.empty()) {
            out << "- " << (br.success ? "Solved" : "Failed")
                << " (" << statusToString(br.status) << ")\n\n";
            continue;
        }

        out << "| # | Solver | Status | Iterations | Residual | Time (ms) |\n";
        out << "|---:|---|---|---:|---:|---:|\n";
        for (size_t si = 0; si < br.solverAttempts.size(); ++si) {
            const auto& sa = br.solverAttempts[si];
            out << "| " << (si + 1)
                << " | " << strategyToString(sa.strategy)
                << " | " << statusToString(sa.status)
                << " | " << sa.iterations
                << " | " << std::scientific << std::setprecision(2) << sa.finalResidual
                << " | " << std::fixed << std::setprecision(1) << sa.elapsedMs
                << " |\n";
        }
        out << "\n";
    }

    out.close();
}

// ============================================================================
// Report Generation
// ============================================================================

std::string generateSolveReport(const SolveResult& result) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6);
    
    ss << "=== Solve Report ===\n";
    ss << "Status: " << (result.success ? "SUCCESS" : "FAILED") 
       << " (" << statusToString(result.status) << ")\n";
    
    if (!result.errorMessage.empty()) {
        ss << "Error: " << result.errorMessage << "\n";
    }
    
    ss << "Blocks evaluated: " << result.blocksEvaluated << "\n";
    ss << "Total iterations: " << result.totalIterations << "\n";
    ss << "Total time: " << result.totalTime.count() << " s\n";
    
    ss << "\n--- Variables (" << result.variables.size() << ") ---\n";
    for (const auto& [name, value] : result.variables) {
        ss << "  " << std::setw(20) << std::left << name 
           << " = " << std::scientific << std::setprecision(9) << value << "\n";
    }
    
    if (!result.stringVariables.empty()) {
        ss << "\n--- String Variables (" << result.stringVariables.size() << ") ---\n";
        for (const auto& [name, value] : result.stringVariables) {
            ss << "  " << std::setw(20) << std::left << name 
               << " = \"" << value << "\"\n";
        }
    }
    
    return ss.str();
}

}  // namespace coolsolve
