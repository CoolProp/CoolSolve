#include "coolsolve/solver.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <csignal>
#include <atomic>

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
        default: return "Unknown";
    }
}

std::string SolverTrace::toString() const {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << "Solver Trace (" << iterations.size() << " iterations, "
       << statusToString(finalStatus) << ")\n";
    ss << "Total time: " << totalTime.count() << " s\n";
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

// ============================================================================
// Newton Solver Implementation
// ============================================================================

Eigen::VectorXd NewtonSolver::computeScalingFactors(const Eigen::VectorXd& x) const {
    const int n = x.size();
    Eigen::VectorXd scale(n);
    
    for (int i = 0; i < n; ++i) {
        double xi = std::abs(x(i));
        if (xi < 1.0) {
            // For small values, use 1.0 as scale (no scaling)
            scale(i) = 1.0;
        } else {
            // Use order of magnitude as scale: scale = 10^floor(log10(x))
            // This gives scale factors like 1, 10, 100, 1000, 1e5, etc.
            double log10x = std::floor(std::log10(xi));
            scale(i) = std::pow(10.0, log10x);
            // Clamp to reasonable range to avoid extreme scaling
            scale(i) = std::max(1e-6, std::min(scale(i), 1e6));
        }
    }
    
    return scale;
}

double NewtonSolver::lineSearch(Problem& problem,
                                const Eigen::VectorXd& x,
                                const Eigen::VectorXd& dx,
                                const Eigen::VectorXd& F,
                                const SolverOptions& options) {
    // Compute initial objective: phi(0) = 0.5 * ||F||^2
    double phi0 = 0.5 * F.squaredNorm();
    
    // Compute directional derivative: dphi = F^T * J * dx = -F^T * F (since J*dx = -F)
    // This is the slope of phi at lambda=0
    double dphi0 = -2.0 * phi0;  // = -||F||^2
    
    double lambda = 1.0;
    Eigen::VectorXd x_new(x.size());
    Eigen::VectorXd F_new(x.size());
    Eigen::MatrixXd J_dummy;  // Not used in line search
    
    for (int lsIter = 0; lsIter < options.lsMaxIterations; ++lsIter) {
        x_new = x + lambda * dx;
        
        try {
            problem.evaluate(x_new, F_new, J_dummy, false);
        } catch (const std::exception& e) {
            // Evaluation failed, try smaller step
            lambda *= options.lsRho;
            continue;
        }
        
        double phi_new = 0.5 * F_new.squaredNorm();
        
        // Armijo condition: phi(lambda) <= phi(0) + alpha * lambda * dphi(0)
        // Since dphi0 is negative, this checks for sufficient decrease
        if (phi_new <= phi0 + options.lsAlpha * lambda * dphi0) {
            return lambda;
        }
        
        // Near minimum: when phi0 is small, numerical noise can prevent Armijo.
        // Accept any step that doesn't significantly increase the residual.
        if (phi0 < 0.05 && phi_new < phi0 * 2.0) {
            return lambda;
        }
        
        // Last resort for very small steps: if we're near convergence (small phi0)
        // and the step doesn't make things worse, accept to avoid line search failure.
        if (lambda < 0.01 && phi0 < 0.1 && phi_new <= phi0 * 1.5) {
            return lambda;
        }
        
        // Reduce step size
        lambda *= options.lsRho;
        
        if (lambda < options.lsMinStep) {
            // Line search failed to find acceptable step
            return 0.0;
        }
    }
    
    return 0.0;  // Line search failed
}

SolverStatus NewtonSolver::solve(Problem& problem,
                                 Eigen::VectorXd& x,
                                 const SolverOptions& options,
                                 SolverTrace* trace,
                                 std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    const int n = problem.size;
    if (n <= 0 || x.size() != n) {
        return SolverStatus::InvalidInput;
    }
    
    // Compute scaling factors (either automatic or unity)
    Eigen::VectorXd scale;
    if (options.enableScaling) {
        scale = computeScalingFactors(x);
        if (options.verbose) {
            std::cout << "Newton: Scaling enabled - factors min=" << scale.minCoeff()
                      << ", max=" << scale.maxCoeff() << std::endl;
        }
    } else {
        scale = Eigen::VectorXd::Ones(n);
    }
    
    // Work in scaled coordinates: y = x ./ scale
    // When scaling is disabled, scale = 1 so y = x (no overhead)
    Eigen::VectorXd y = x.cwiseQuotient(scale);
    
    Eigen::VectorXd F(n);
    Eigen::MatrixXd J(n, n);
    Eigen::MatrixXd J_unscaled(n, n);
    Eigen::VectorXd dy(n);
    Eigen::VectorXd x_unscaled(n);
    
    double initialResidualNorm = 0.0;
    
    for (int iter = 0; iter < options.maxIterations; ++iter) {
        // Check for timeout
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "Solver timed out";
            return SolverStatus::EvaluationError; // Or add a Timeout status
        }
        // Step 1: Evaluate F(x) and J(x) in original coordinates, then scale Jacobian
        try {
            // Transform from scaled to original: x = scale .* y
            x_unscaled = y.cwiseProduct(scale);
            // Evaluate in original coordinates
            problem.evaluate(x_unscaled, F, J_unscaled, true);
            // Scale Jacobian: dF/dy = dF/dx * dx/dy = J_unscaled * diag(scale)
            J = J_unscaled * scale.asDiagonal();
        } catch (const std::exception& e) {
            if (options.verbose) {
                std::cerr << "Newton: Evaluation failed at iter " << iter
                          << ": " << e.what() << std::endl;
            }
            // We can't easily return the error message here without changing the interface
            // but we can rethrow it and catch it in solveBlock
            throw;
        }
        
        // Step 2: Check convergence
        double residualNorm = F.lpNorm<Eigen::Infinity>();
        
        if (iter == 0) {
            initialResidualNorm = residualNorm;
        }
        
        if (options.verbose) {
            std::cout << "Newton iter " << iter << ": ||F||_inf = " << residualNorm << std::endl;
        }
        
        // Record trace (store unscaled x for easier debugging)
        if (trace) {
            SolverTrace::Iteration traceIter;
            traceIter.iter = iter;
            traceIter.residualNorm = residualNorm;
            traceIter.stepNorm = 0.0;
            traceIter.lambda = 1.0;
            Eigen::VectorXd x_unscaled = y.cwiseProduct(scale);
            traceIter.x = std::vector<double>(x_unscaled.data(), x_unscaled.data() + x_unscaled.size());
            traceIter.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(traceIter);
        }
        
        // Check absolute tolerance
        if (residualNorm < options.tolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            // Transform solution back to unscaled coordinates
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }
        
        // Check relative tolerance
        if (initialResidualNorm > 0 &&
            residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            // Transform solution back to unscaled coordinates
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }
        
        // Step 3: Solve J * dy = -F (in scaled coordinates)
        // Use ColPivHouseholderQR for robustness (handles rank-deficient cases)
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(J);
        
        if (qr.rank() < n) {
            if (options.verbose) {
                std::cerr << "Newton: Jacobian is rank-deficient (rank=" << qr.rank()
                          << ", n=" << n << ")" << std::endl;
            }
            // Try pseudo-inverse solution anyway
        }
        
        dy = qr.solve(-F);
        
        // Check for invalid step
        if (!dy.allFinite()) {
            if (options.verbose) {
                std::cerr << "Newton: Invalid Newton step (contains NaN/Inf)" << std::endl;
            }
            if (trace) {
                trace->finalStatus = SolverStatus::SingularJacobian;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            return SolverStatus::SingularJacobian;
        }
        
        // Step 4: Line search (in scaled coordinates)
        // Create a wrapper for line search that handles coordinate transformation
        Problem lsProblem;
        lsProblem.size = n;
        lsProblem.evaluate = [&](const Eigen::VectorXd& y_trial,
                                  Eigen::VectorXd& F_out,
                                  Eigen::MatrixXd& J_out,
                                  bool computeJac) {
            Eigen::VectorXd x_trial = y_trial.cwiseProduct(scale);
            Eigen::MatrixXd J_temp;
            problem.evaluate(x_trial, F_out, J_temp, false);
        };
        double lambda = lineSearch(lsProblem, y, dy, F, options);
        
        if (lambda == 0.0) {
            // Relaxed acceptance: if residual is already small, accept as converged
            if (residualNorm < options.lsRelaxedTolerance) {
                if (trace) {
                    trace->finalStatus = SolverStatus::Success;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                // Transform solution back to unscaled coordinates
                x = y.cwiseProduct(scale);
                return SolverStatus::Success;
            }
            // Fallback: try a minimal step to make progress (in scaled coordinates)
            for (double fallbackLambda : {0.1, 0.01, 0.001}) {
                Eigen::VectorXd y_trial = y + fallbackLambda * dy;
                Eigen::VectorXd F_trial(n);
                Eigen::MatrixXd J_dummy;
                try {
                    Eigen::VectorXd x_trial = y_trial.cwiseProduct(scale);
                    problem.evaluate(x_trial, F_trial, J_dummy, false);
                    double phi_trial = 0.5 * F_trial.squaredNorm();
                    double phi0 = 0.5 * F.squaredNorm();
                    if (phi_trial < phi0) {
                        lambda = fallbackLambda;
                        break;
                    }
                } catch (...) {
                    continue;
                }
            }
            if (lambda == 0.0) {
                if (options.verbose) {
                    std::cerr << "Newton: Line search failed at iter " << iter << std::endl;
                }
                if (trace) {
                    trace->finalStatus = SolverStatus::LineSearchFailed;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                if (detailedError) {
                    std::ostringstream ss;
                    ss << "Line search failed at iteration " << iter << ".\n"
                       << "Residual norm ||F|| = " << residualNorm << "\n"
                       << "Newton step norm ||dy|| = " << dy.norm() << "\n"
                       << "Current variables:\n";
                    *detailedError = ss.str();
                }
                // Transform back to unscaled coordinates before returning
                x = y.cwiseProduct(scale);
                return SolverStatus::LineSearchFailed;
            }
        }
        
        // Step 5: Update y (in scaled coordinates)
        double stepNorm = (lambda * dy).lpNorm<Eigen::Infinity>();
        y += lambda * dy;
        
        // Update trace with actual step info (convert step to unscaled for reporting)
        if (trace && !trace->iterations.empty()) {
            trace->iterations.back().stepNorm = stepNorm;
            trace->iterations.back().lambda = lambda;
        }
        
        // Check for tiny step (convergence stall)
        if (stepNorm < options.stepTolerance) {
            // This might be success if residual is small, or stall otherwise
            if (residualNorm < options.tolerance * 100) {
                if (trace) {
                    trace->finalStatus = SolverStatus::Success;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                // Transform solution back to unscaled coordinates before returning
                x = y.cwiseProduct(scale);
                return SolverStatus::Success;
            }
        }
    }
    
    if (trace) {
        trace->finalStatus = SolverStatus::MaxIterations;
        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
    }
    if (detailedError) {
        std::ostringstream ss;
        ss << "Max iterations (" << options.maxIterations << ") reached without convergence.\n"
           << "Last residual norm ||F|| = " << F.lpNorm<Eigen::Infinity>();
        *detailedError = ss.str();
    }
    // Transform current y back to unscaled x even on failure
    x = y.cwiseProduct(scale);
    return SolverStatus::MaxIterations;
}

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
    // For now, always use Newton for implicit solve
    // In the future, we could try to detect simple explicit patterns like:
    //   x = expr(known_vars)
    // and solve directly.
    
    // For size-1 blocks, we still use Newton, but it will typically 
    // converge in 1-2 iterations for explicit equations.
    return false;
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
        return SolverStatus::Success;
    }
    
    // Create problem for Newton solver
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
    
    // Solve using Newton
    NewtonSolver newton;
    SolverStatus status;
    std::string newtonError;
    try {
        status = newton.solve(problem, x, options, trace, &newtonError);
        if (!newtonError.empty()) {
            if (outErrorMessage) {
                if (!outErrorMessage->empty()) *outErrorMessage += "\n";
                *outErrorMessage += newtonError;
            }
        }
    } catch (const std::exception& e) {
        if (outErrorMessage) *outErrorMessage = e.what();
        return SolverStatus::EvaluationError;
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

SolveResult Solver::solve(const SolverOptions& options, bool enableTracing) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    SolveResult result;
    result.blocksEvaluated = 0;
    result.totalIterations = 0;
    
    if (enableTracing) {
        result.blockTraces.resize(evaluator_.getNumBlocks());
    }
    
    // Solve blocks in topological order
    for (size_t blockIdx = 0; blockIdx < evaluator_.getNumBlocks(); ++blockIdx) {
        SolverTrace* trace = enableTracing ? &result.blockTraces[blockIdx] : nullptr;
        
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
