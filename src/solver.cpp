#include "coolsolve/solver.h"
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
            else if (key == "enableScaling") options.enableScaling = parseBool(val);
            else if (key == "trInitialRadius") options.trInitialRadius = std::stod(val);
            else if (key == "trMaxRadius") options.trMaxRadius = std::stod(val);
            else if (key == "trEta") options.trEta = std::stod(val);
            else if (key == "trShrinkFactor") options.trShrinkFactor = std::stod(val);
            else if (key == "trGrowFactor") options.trGrowFactor = std::stod(val);
            else if (key == "partitionedMaxIterations") options.partitionedMaxIterations = std::stoi(val);
            else if (key == "partitionedRelaxation") options.partitionedRelaxation = std::stod(val);
            else if (key == "partitionedMinDiagonal") options.partitionedMinDiagonal = std::stod(val);
            else if (key == "partitionedMinBlockSize") options.partitionedMinBlockSize = std::stoi(val);
            else if (key == "enableTearing") options.enableTearing = parseBool(val);
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
        } catch (...) {
            // Ignore malformed values
        }
    }
    return true;
}

// ============================================================================
// Timeout Handling
// ============================================================================

static std::atomic<bool> g_timed_out{false};

#ifdef __unix__
void handle_sigalrm(int) {
    g_timed_out = true;
}
#endif

TimeoutGuard::TimeoutGuard(int seconds) : seconds_(seconds) {
    g_timed_out = false;
    if (seconds > 0) {
#ifdef __unix__
        signal(SIGALRM, handle_sigalrm);
        alarm(seconds);
#endif
    }
}

TimeoutGuard::~TimeoutGuard() {
    if (seconds_ > 0) {
#ifdef __unix__
        alarm(0);
        signal(SIGALRM, SIG_DFL);
#endif
    }
}

bool TimeoutGuard::hasTimedOut() {
    return g_timed_out;
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
        auto eval1D = [&](double xval) -> std::pair<double, double> {
            std::vector<double> x_std = {xval};
            auto er = blockEval.evaluate(x_std, externalVars, externalStringVars);
            double f = er.residuals[0];
            double j = (er.jacobian.size() > 0 && er.jacobian[0].size() > 0) ? er.jacobian[0][0] : 0.0;
            return {f, j};
        };
        
        // Phase 1: Try Newton with trust-region limiting
        double xCur = x[0];
        double radius = std::max(std::abs(xCur) * 2.0, 100.0);
        bool converged = false;
        
        // Track bracket for bisection fallback
        bool hasBracket = false;
        double xLo = 0, xHi = 0, fLo = 0, fHi = 0;
        
        // Phase 1 uses fewer iterations since Phase 2 probing is more effective
        // for problems where the initial guess is far from the solution
        int phase1MaxIter = std::min(options.maxIterations, 50);
        for (int iter = 0; iter < phase1MaxIter && !converged; ++iter) {
            double f, j;
            try {
                auto [fv, jv] = eval1D(xCur);
                f = fv; j = jv;
            } catch (...) {
                // Evaluation failed — try reducing x toward zero
                xCur *= 0.5;
                continue;
            }
            
            if (trace) {
                SolverTrace::Iteration traceIter;
                traceIter.iter = iter;
                traceIter.residualNorm = std::abs(f);
                traceIter.stepNorm = 0.0;
                traceIter.lambda = 1.0;
                traceIter.x = {xCur};
                traceIter.residuals = {f};
                trace->iterations.push_back(traceIter);
            }
            
            if (std::abs(f) < options.tolerance) {
                converged = true;
                break;
            }
            
            // Update bracket
            if (!hasBracket) {
                if (iter == 0) {
                    xLo = xCur; fLo = f;
                } else {
                    if (f * fLo < 0) {
                        hasBracket = true;
                        xHi = xCur; fHi = f;
                        if (xLo > xHi) { std::swap(xLo, xHi); std::swap(fLo, fHi); }
                    } else {
                        xLo = xCur; fLo = f;
                    }
                }
            } else {
                // Keep bracket tight
                if (f * fLo < 0) { xHi = xCur; fHi = f; }
                else { xLo = xCur; fLo = f; }
                if (xLo > xHi) { std::swap(xLo, xHi); std::swap(fLo, fHi); }
            }
            
            double dx;
            if (hasBracket) {
                // Use bisection-Newton hybrid (Illinois/Dekker-style)
                double xNewton = (std::abs(j) > 1e-30) ? xCur - f / j : (xLo + xHi) / 2.0;
                // If Newton step is within bracket, use it; otherwise bisect
                if (xNewton > xLo && xNewton < xHi) {
                    dx = xNewton - xCur;
                } else {
                    dx = (xLo + xHi) / 2.0 - xCur;
                }
            } else if (std::abs(j) > 1e-30) {
                // Newton step with trust-region clamping
                dx = -f / j;
                if (std::abs(dx) > radius) {
                    dx = (dx > 0 ? 1.0 : -1.0) * radius;
                }
                // Grow trust region on each step
                radius = std::max(radius, std::abs(dx) * 2.0);
            } else {
                // Zero derivative — explore by stepping
                dx = radius;
                radius *= 2.0;
            }
            
            xCur += dx;
            if (trace && !trace->iterations.empty()) {
                trace->iterations.back().stepNorm = std::abs(dx);
            }
        }
        
        // Phase 2: If Newton didn't converge and no bracket found,
        // try multiple starting points to find a sign change.
        // Scan both positive and negative ranges on a log scale, plus
        // intermediate values that may cross sign near poles.
        if (!converged && !hasBracket) {
            // Build comprehensive probe list
            std::vector<double> probes;
            // Dense scan around zero
            for (double v : {0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0}) {
                probes.push_back(v);
                probes.push_back(-v);
            }
            // Log-spaced scan
            for (double scale = 100; scale <= 1e8; scale *= 2.0) {
                probes.push_back(scale);
                probes.push_back(-scale);
            }
            // KEY STRATEGY: Probe near external variable values.
            // Solutions of equations are typically at scales related to their inputs.
            // Without this, narrow sign-change regions near poles are missed.
            for (const auto& [name, val] : externalVars) {
                if (!std::isfinite(val) || val == 0.0) continue;
                double absVal = std::abs(val);
                // Probe at, near, and around the external variable value
                for (double frac : {0.5, 0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 1.5, 2.0}) {
                    probes.push_back(val * frac);
                    if (val > 0) probes.push_back(-val * frac); // Also try negative
                }
            }
            
            // Record all evaluated points with their residuals for bracket detection
            struct ProbeResult { double x; double f; bool valid; };
            std::vector<ProbeResult> results;
            
            double bestF = 1e30;
            double bestX = xCur;

            for (double probe : probes) {
                try {
                    auto [f, j] = eval1D(probe);
                    if (!std::isfinite(f)) continue;
                    results.push_back({probe, f, true});
                    if (std::abs(f) < bestF) {
                        bestF = std::abs(f);
                        bestX = probe;
                    }
                    if (std::abs(f) < options.tolerance) {
                        xCur = probe;
                        converged = true;
                        break;
                    }
                } catch (...) {
                    results.push_back({probe, 0.0, false});
                    continue;
                }
            }
            
            // Look for sign changes between any pair of valid probe results.
            // Collect ALL sign changes, then pick the best one (smallest midpoint |f|)  
            // to avoid locking onto poles instead of roots.
            if (!converged) {
                // Sort by x value
                std::sort(results.begin(), results.end(), 
                          [](const ProbeResult& a, const ProbeResult& b) { return a.x < b.x; });
                
                struct Bracket { double lo, flo, hi, fhi, midF; };
                std::vector<Bracket> brackets;
                
                for (size_t i = 0; i + 1 < results.size(); ++i) {
                    if (!results[i].valid || !results[i+1].valid) continue;
                    if (results[i].f * results[i+1].f < 0) {
                        // Found a sign change — evaluate midpoint to score it
                        double lo = results[i].x, hi = results[i+1].x;
                        double flo = results[i].f, fhi = results[i+1].f;
                        double mid = (lo + hi) / 2.0;
                        double fmid = 1e30;
                        try {
                            auto [fm, jm] = eval1D(mid);
                            fmid = std::abs(fm);
                        } catch (...) {
                            fmid = 1e30; // Evaluation failed at midpoint — likely a pole
                        }
                        brackets.push_back({lo, flo, hi, fhi, fmid});
                    }
                }
                
                if (!brackets.empty()) {
                    // Pick bracket with smallest midpoint |f| (most likely a root, not a pole)
                    auto& best = *std::min_element(brackets.begin(), brackets.end(),
                        [](const Bracket& a, const Bracket& b) { return a.midF < b.midF; });
                    hasBracket = true;
                    xLo = best.lo; fLo = best.flo;
                    xHi = best.hi; fHi = best.fhi;
                    if (xLo > xHi) { std::swap(xLo, xHi); std::swap(fLo, fHi); }
                }
            }
            
            if (!converged && !hasBracket) {
                // Use the best point found as starting point for another Newton attempt
                xCur = bestX;
            }
        }
        
        // Phase 3: If we have a bracket, use bisection + Newton
        if (!converged && hasBracket) {
            for (int iter = 0; iter < options.maxIterations && !converged; ++iter) {
                double xMid = (xLo + xHi) / 2.0;
                try {
                    auto [f, j] = eval1D(xMid);
                    
                    if (trace) {
                        SolverTrace::Iteration traceIter;
                        traceIter.iter = 1000 + iter; // Distinguish bisection iterations
                        traceIter.residualNorm = std::abs(f);
                        traceIter.stepNorm = (xHi - xLo) / 2.0;
                        traceIter.lambda = 1.0;
                        traceIter.x = {xMid};
                        traceIter.residuals = {f};
                        traceIter.detail = "       [bisection] bracket: [" + std::to_string(xLo) + ", " + std::to_string(xHi) + "]\n";
                        trace->iterations.push_back(traceIter);
                    }
                    
                    if (std::abs(f) < options.tolerance) {
                        xCur = xMid;
                        converged = true;
                        break;
                    }
                    
                    // Try Newton step within bracket
                    if (std::abs(j) > 1e-30) {
                        double xNewton = xMid - f / j;
                        if (xNewton > xLo && xNewton < xHi) {
                            try {
                                auto [fn, jn] = eval1D(xNewton);
                                if (std::abs(fn) < std::abs(f)) {
                                    if (fn * fLo < 0) { xHi = xNewton; fHi = fn; }
                                    else { xLo = xNewton; fLo = fn; }
                                    continue;
                                }
                            } catch (...) {}
                        }
                    }
                    
                    // Fall back to bisection
                    if (f * fLo < 0) { xHi = xMid; fHi = f; }
                    else { xLo = xMid; fLo = f; }
                    
                    if (xHi - xLo < options.stepTolerance) {
                        xCur = (xLo + xHi) / 2.0;
                        converged = true;
                    }
                } catch (...) {
                    // Evaluation failed at midpoint — narrow bracket from the other side
                    xHi = xMid;
                }
            }
        }
        
        // Phase 4: If still not converged, try one more round of Newton from best position
        if (!converged) {
            for (int iter = 0; iter < 50 && !converged; ++iter) {
                try {
                    auto [f, j] = eval1D(xCur);
                    if (std::abs(f) < options.tolerance) { converged = true; break; }
                    if (std::abs(f) < options.lsRelaxedTolerance) { converged = true; break; }
                    if (std::abs(j) < 1e-30) break;
                    double dx = -f / j;
                    double maxStep = std::max(std::abs(xCur) * 2.0, 1e6);
                    if (std::abs(dx) > maxStep) dx = (dx > 0 ? 1.0 : -1.0) * maxStep;
                    xCur += dx;
                } catch (...) {
                    break;
                }
            }
        }
        
        if (converged) {
            x[0] = xCur;
            evaluator_.setVariableValue(varNames[0], xCur);
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime1D;
            }
            return SolverStatus::Success;
        }
        
        // Newton1D failed — reset x and fall through to standard pipeline
        x[0] = evaluator_.getVariableValue(varNames[0]);
        if (trace) {
            trace->iterations.clear(); // Clear Newton1D iterations for clean pipeline trace
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
        
        // Evaluate block
        auto result = blockEval.evaluate(x_std, externalVars, externalStringVars);
        
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
                                          std::string* outErrorMessage) {
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

    for (size_t idx = 0; idx < options.solverPipeline.size(); ++idx) {
        SolverStrategy strategy = options.solverPipeline[idx];

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
            for (size_t i = 0; i < n; ++i) {
                x_init[i] = evaluator_.getVariableValue(varNames[i]);
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
        status = runSolverStrategy(strategy, blockIndex, problem, blockEval,
                                   varNames, externalVars, externalStringVars,
                                   x, options, trace, &solverError);

        if (status == SolverStatus::Success) {
            return status;
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
                    auto evalResult = blockEval.evaluate(x_std, externalVars, externalStringVars);
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

        for (size_t k = 0; k < nonTearEqLocalIndices.size(); ++k) {
            size_t eq = nonTearEqLocalIndices[k];
            size_t varIdx = eq;
            for (int inner = 0; inner < maxInner; ++inner) {
                x_std.assign(x.data(), x.data() + x.size());
                EvaluationResult er;
                try {
                    er = blockEval.evaluate(x_std, externalVars, externalStringVars);
                } catch (...) {
                    break;
                }
                double fEq = er.residuals[eq];
                double jEq = (eq < er.jacobian.size() && varIdx < er.jacobian[eq].size())
                             ? er.jacobian[eq][varIdx] : 0.0;
                if (std::abs(jEq) < 1e-14) break;
                double step = -fEq / jEq;
                if (!std::isfinite(step)) break;
                x[varIdx] += step;
                if (std::abs(fEq) < options.tolerance) break;
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

        if (tearNorm < options.tolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::Success;
        }

        Eigen::VectorXd dx_tear;
        Eigen::FullPivLU<Eigen::MatrixXd> lu(S);
        if (!lu.isInvertible()) {
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

SolveResult Solver::solve(const SolverOptions& options, bool enableTracing) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    SolveResult result;
    result.blocksEvaluated = 0;
    result.totalIterations = 0;
    
    if (enableTracing) {
        result.blockTraces.resize(evaluator_.getNumBlocks());
    }
    
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
        
        SolverTrace* trace = enableTracing ? &result.blockTraces[blockIdx] : nullptr;
        
        // Notify progress: block starting
        if (options.progressCallback) {
            options.progressCallback(static_cast<int>(blockIdx), totalBlocks, "start", 0, 0.0);
        }
        
        // Setup timeout protection
        TimeoutGuard timeout(options.timeoutSeconds);
        
        std::string blockError;
        SolverStatus blockStatus = solveBlock(blockIdx, options, trace, &blockError);
        result.blocksEvaluated++;
        
        // Record per-block result
        SolveResult::BlockResult br;
        br.id = blockIdx;
        br.success = (blockStatus == SolverStatus::Success);
        br.status = blockStatus;
        br.iterations = trace ? static_cast<int>(trace->iterations.size()) : 0;
        br.maxResidual = 0.0;
        if (trace && !trace->iterations.empty()) {
            br.maxResidual = trace->iterations.back().residualNorm;
        }
        br.errorMessage = blockError;
        result.blockResults.push_back(br);
        
        if (trace) {
            result.totalIterations += static_cast<int>(trace->iterations.size());
        }
        
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
    
    return result;
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
