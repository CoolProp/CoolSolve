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
// Solver Strategy Utilities
// ============================================================================

std::string strategyToString(SolverStrategy strategy) {
    switch (strategy) {
        case SolverStrategy::Newton:            return "Newton";
        case SolverStrategy::TrustRegion:       return "TrustRegion";
        case SolverStrategy::LevenbergMarquardt: return "LevenbergMarquardt";
        case SolverStrategy::Partitioned:       return "Partitioned";
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
        case SolverStrategy::Partitioned:
            // Partitioned is handled specially by the orchestrator (needs structural info).
            // Return nullptr; the caller must handle this case.
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
            else if (key == "useTrustRegion") options.useTrustRegion = parseBool(val);
            else if (key == "trInitialRadius") options.trInitialRadius = std::stod(val);
            else if (key == "trMaxRadius") options.trMaxRadius = std::stod(val);
            else if (key == "trEta") options.trEta = std::stod(val);
            else if (key == "trShrinkFactor") options.trShrinkFactor = std::stod(val);
            else if (key == "trGrowFactor") options.trGrowFactor = std::stod(val);
            else if (key == "usePartitionedSolver") options.usePartitionedSolver = parseBool(val);
            else if (key == "partitionedMaxIterations") options.partitionedMaxIterations = std::stoi(val);
            else if (key == "partitionedRelaxation") options.partitionedRelaxation = std::stod(val);
            else if (key == "partitionedMinDiagonal") options.partitionedMinDiagonal = std::stod(val);
            else if (key == "partitionedMinBlockSize") options.partitionedMinBlockSize = std::stoi(val);
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
    
    // Record solver type in trace for debugging
    if (trace) {
        trace->solverType = "Newton";
    }
    
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
                // Capture J and F for debug diagnostics
                trace->singularJacobianF.assign(F.data(), F.data() + F.size());
                trace->singularJacobianJ.resize(J.rows());
                for (Eigen::Index r = 0; r < J.rows(); ++r)
                    trace->singularJacobianJ[r].assign(J.row(r).data(), J.row(r).data() + J.cols());
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
// Trust Region Dogleg Solver Implementation
// ============================================================================

Eigen::VectorXd TrustRegionSolver::computeScalingFactors(const Eigen::VectorXd& x) const {
    const int n = x.size();
    Eigen::VectorXd scale(n);
    
    for (int i = 0; i < n; ++i) {
        double xi = std::abs(x(i));
        if (xi < 1.0) {
            scale(i) = 1.0;
        } else {
            double log10x = std::floor(std::log10(xi));
            scale(i) = std::pow(10.0, log10x);
            scale(i) = std::max(1e-6, std::min(scale(i), 1e6));
        }
    }
    
    return scale;
}

double TrustRegionSolver::evaluateModel(const Eigen::VectorXd& F,
                                        const Eigen::MatrixXd& J,
                                        const Eigen::VectorXd& p) {
    // m(p) = 0.5 * ||F + J*p||^2
    Eigen::VectorXd F_plus_Jp = F + J * p;
    return 0.5 * F_plus_Jp.squaredNorm();
}

Eigen::VectorXd TrustRegionSolver::doglegStep(const Eigen::VectorXd& dx_n,
                                              const Eigen::VectorXd& dx_c,
                                              double delta) {
    double norm_n = dx_n.norm();
    double norm_c = dx_c.norm();
    
    // Case 1: Newton step is inside trust region
    if (norm_n <= delta) {
        return dx_n;
    }
    
    // Case 2: Cauchy step is outside trust region, scale it
    if (norm_c >= delta) {
        return (delta / norm_c) * dx_c;
    }
    
    // Case 3: Dogleg path - find interpolation between Cauchy and Newton
    // Path: p(tau) = dx_c + tau*(dx_n - dx_c) for tau in [0, 1]
    // Find tau such that ||p(tau)|| = delta
    Eigen::VectorXd diff = dx_n - dx_c;
    
    // Solve: ||dx_c + tau*diff||^2 = delta^2
    // => (diff^T*diff)*tau^2 + 2*(dx_c^T*diff)*tau + (dx_c^T*dx_c - delta^2) = 0
    double a = diff.squaredNorm();
    double b = 2.0 * dx_c.dot(diff);
    double c = dx_c.squaredNorm() - delta * delta;
    
    // Quadratic formula: tau = (-b + sqrt(b^2 - 4ac)) / (2a)
    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) {
        discriminant = 0.0;  // Shouldn't happen, but be safe
    }
    
    // We want the root in [0, 1], which is the positive one
    double tau = (-b + std::sqrt(discriminant)) / (2.0 * a);
    tau = std::max(0.0, std::min(1.0, tau));
    
    return dx_c + tau * diff;
}

SolverStatus TrustRegionSolver::solve(Problem& problem,
                                      Eigen::VectorXd& x,
                                      const SolverOptions& options,
                                      SolverTrace* trace,
                                      std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Record solver type in trace for debugging
    if (trace) {
        trace->solverType = "TrustRegion";
    }
    
    const int n = problem.size;
    if (n <= 0 || x.size() != n) {
        return SolverStatus::InvalidInput;
    }
    
    // Compute scaling factors
    Eigen::VectorXd scale;
    if (options.enableScaling) {
        scale = computeScalingFactors(x);
        if (options.verbose) {
            std::cout << "TrustRegion: Scaling enabled - factors min=" << scale.minCoeff()
                      << ", max=" << scale.maxCoeff() << std::endl;
        }
    } else {
        scale = Eigen::VectorXd::Ones(n);
    }
    
    // Work in scaled coordinates: y = x ./ scale
    Eigen::VectorXd y = x.cwiseQuotient(scale);
    
    Eigen::VectorXd F(n);
    Eigen::MatrixXd J(n, n);
    Eigen::MatrixXd J_unscaled(n, n);
    
    double initialResidualNorm = 0.0;
    double delta = options.trInitialRadius;  // Trust region radius
    
    for (int iter = 0; iter < options.maxIterations; ++iter) {
        // Check for timeout
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "Solver timed out";
            x = y.cwiseProduct(scale);
            return SolverStatus::EvaluationError;
        }
        
        // Evaluate F(x) and J(x)
        try {
            Eigen::VectorXd x_unscaled = y.cwiseProduct(scale);
            problem.evaluate(x_unscaled, F, J_unscaled, true);
            // Scale Jacobian: dF/dy = dF/dx * dx/dy = J_unscaled * diag(scale)
            J = J_unscaled * scale.asDiagonal();
        } catch (const std::exception& e) {
            if (options.verbose) {
                std::cerr << "TrustRegion: Evaluation failed at iter " << iter
                          << ": " << e.what() << std::endl;
            }
            x = y.cwiseProduct(scale);
            throw;
        }
        
        // Check convergence
        double residualNorm = F.lpNorm<Eigen::Infinity>();
        
        if (iter == 0) {
            initialResidualNorm = residualNorm;
        }
        
        if (options.verbose) {
            std::cout << "TrustRegion iter " << iter << ": ||F||_inf = " << residualNorm
                      << ", delta = " << delta << std::endl;
        }
        
        // Record trace
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
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }
        
        // Compute gradient: g = J^T * F
        Eigen::VectorXd g = J.transpose() * F;
        
        // Compute Cauchy step (steepest descent with optimal step length)
        // dx_c = -alpha * g where alpha = (g^T*g) / (g^T*J^T*J*g)
        Eigen::VectorXd dx_c;
        double g_norm_sq = g.squaredNorm();
        if (g_norm_sq < 1e-30) {
            // Gradient is essentially zero, we're at a stationary point
            // but residual is not small enough - might be local minimum
            if (residualNorm < options.lsRelaxedTolerance) {
                if (trace) {
                    trace->finalStatus = SolverStatus::Success;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                x = y.cwiseProduct(scale);
                return SolverStatus::Success;
            }
            if (options.verbose) {
                std::cerr << "TrustRegion: Gradient is zero but residual is large" << std::endl;
            }
            if (trace) {
                trace->finalStatus = SolverStatus::SingularJacobian;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                // Capture J and F for debug diagnostics
                trace->singularJacobianF.assign(F.data(), F.data() + F.size());
                trace->singularJacobianJ.resize(J.rows());
                for (Eigen::Index r = 0; r < J.rows(); ++r)
                    trace->singularJacobianJ[r].assign(J.row(r).data(), J.row(r).data() + J.cols());
            }
            x = y.cwiseProduct(scale);
            return SolverStatus::SingularJacobian;
        }
        
        Eigen::VectorXd Jg = J * g;
        double alpha = g_norm_sq / Jg.squaredNorm();
        dx_c = -alpha * g;
        
        // Compute Newton step by solving J * dx_n = -F
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(J);
        Eigen::VectorXd dx_n = qr.solve(-F);
        
        // Check for invalid Newton step
        if (!dx_n.allFinite()) {
            if (options.verbose) {
                std::cerr << "TrustRegion: Invalid Newton step (contains NaN/Inf)" << std::endl;
            }
            // Fall back to Cauchy step
            dx_n = dx_c;
        }
        
        // Compute dogleg step within trust region
        Eigen::VectorXd dy = doglegStep(dx_n, dx_c, delta);
        
        // Evaluate proposed step
        Eigen::VectorXd y_new = y + dy;
        Eigen::VectorXd F_new(n);
        Eigen::MatrixXd J_dummy;
        
        try {
            Eigen::VectorXd x_new = y_new.cwiseProduct(scale);
            problem.evaluate(x_new, F_new, J_dummy, false);
        } catch (const std::exception& e) {
            // Evaluation failed - reject step and shrink trust region
            if (options.verbose) {
                std::cout << "TrustRegion: Evaluation failed, shrinking delta from " << delta << std::endl;
            }
            delta *= options.trShrinkFactor;
            if (delta < 1e-8) {
                // Trust region too small - reset and try a gradient descent step
                if (options.verbose) {
                    std::cout << "TrustRegion: Delta too small, resetting to initial radius" << std::endl;
                }
                delta = options.trInitialRadius;
                
                // If we've reset multiple times, try a different approach
                if (iter > 50 && residualNorm < initialResidualNorm * 0.5) {
                    // We've made significant progress but stuck - accept current solution
                    if (residualNorm < options.lsRelaxedTolerance) {
                        if (options.verbose) {
                            std::cout << "TrustRegion: Accepting suboptimal solution after significant progress" << std::endl;
                        }
                        if (trace) {
                            trace->finalStatus = SolverStatus::Success;
                            trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                        }
                        x = y.cwiseProduct(scale);
                        return SolverStatus::Success;
                    }
                }
                
                // If trust region keeps collapsing, give up
                if (iter > options.maxIterations * 0.8) {
                    if (trace) {
                        trace->finalStatus = SolverStatus::LineSearchFailed;
                        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                    }
                    if (detailedError) {
                        std::ostringstream ss;
                        ss << "Trust region collapsed at iteration " << iter << ".\n"
                           << "Residual norm ||F|| = " << residualNorm;
                        *detailedError = ss.str();
                    }
                    x = y.cwiseProduct(scale);
                    return SolverStatus::LineSearchFailed;
                }
            }
            continue;
        }
        
        double residualNormNew = F_new.lpNorm<Eigen::Infinity>();
        
        // Compute actual vs predicted reduction
        double phi_old = 0.5 * F.squaredNorm();
        double phi_new = 0.5 * F_new.squaredNorm();
        double model_old = evaluateModel(F, J, Eigen::VectorXd::Zero(n));
        double model_new = evaluateModel(F, J, dy);
        
        double actualReduction = phi_old - phi_new;
        double predictedReduction = model_old - model_new;
        
        double rho = 0.0;
        if (std::abs(predictedReduction) > 1e-30) {
            rho = actualReduction / predictedReduction;
        }
        
        if (options.verbose) {
            std::cout << "TrustRegion: rho = " << rho << ", actual = " << actualReduction
                      << ", predicted = " << predictedReduction << std::endl;
        }
        
        // Update trust region and accept/reject step
        // Accept step if it makes ANY positive progress (actual reduction > 0)
        // This is key for highly nonlinear problems where the model may be poor
        bool acceptStep = (actualReduction > 0);
        
        if (!acceptStep) {
            // Reject step and shrink trust region
            if (options.verbose) {
                std::cout << "TrustRegion: Rejecting step (rho=" << rho
                          << ", actual=" << actualReduction << "), shrinking delta" << std::endl;
            }
            delta *= options.trShrinkFactor;
            
            // If trust region gets too small, try a gradient descent step
            if (delta < 1e-6) {
                if (options.verbose) {
                    std::cout << "TrustRegion: Delta too small, resetting" << std::endl;
                }
                delta = options.trInitialRadius;
            }
        }
        
        if (acceptStep) {
            // Accept step
            y = y_new;
            residualNorm = residualNormNew;  // Update residual for next iteration
            double stepNorm = dy.lpNorm<Eigen::Infinity>();
            
            // Update trace
            if (trace && !trace->iterations.empty()) {
                trace->iterations.back().stepNorm = stepNorm;
            }
            
            // Check for tiny step (convergence stall)
            if (stepNorm < options.stepTolerance) {
                if (residualNormNew < options.tolerance * 100 ||
                    residualNormNew < options.lsRelaxedTolerance) {
                    if (trace) {
                        trace->finalStatus = SolverStatus::Success;
                        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                    }
                    x = y.cwiseProduct(scale);
                    return SolverStatus::Success;
                }
            }
            
            // Grow trust region if step was good
            if (rho > 0.75 && dy.norm() >= 0.9 * delta) {
                delta = std::min(options.trGrowFactor * delta, options.trMaxRadius);
                if (options.verbose) {
                    std::cout << "TrustRegion: Growing delta to " << delta << std::endl;
                }
            }
        }
    }
    
    if (trace) {
        trace->finalStatus = SolverStatus::MaxIterations;
        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
    }
    if (detailedError) {
        std::ostringstream ss;
        ss << "Trust region: Max iterations (" << options.maxIterations
           << ") reached without convergence.";
        *detailedError = ss.str();
    }
    x = y.cwiseProduct(scale);
    return SolverStatus::MaxIterations;
}

// ============================================================================
// Levenberg-Marquardt Solver Implementation
// ============================================================================

Eigen::VectorXd LevenbergMarquardtSolver::computeScalingFactors(const Eigen::VectorXd& x) const {
    const int n = x.size();
    Eigen::VectorXd scale(n);
    for (int i = 0; i < n; ++i) {
        double xi = std::abs(x(i));
        if (xi < 1.0) {
            scale(i) = 1.0;
        } else {
            double log10x = std::floor(std::log10(xi));
            scale(i) = std::pow(10.0, log10x);
            scale(i) = std::max(1e-6, std::min(scale(i), 1e6));
        }
    }
    return scale;
}

SolverStatus LevenbergMarquardtSolver::solve(Problem& problem,
                                              Eigen::VectorXd& x,
                                              const SolverOptions& options,
                                              SolverTrace* trace,
                                              std::string* detailedError) {
    auto startTime = std::chrono::high_resolution_clock::now();

    if (trace) {
        if (trace->solverType.empty())
            trace->solverType = "LevenbergMarquardt";
        else if (trace->solverType.find("LevenbergMarquardt") == std::string::npos)
            trace->solverType += " -> LevenbergMarquardt";
    }

    const int n = problem.size;
    if (n <= 0 || x.size() != n) {
        return SolverStatus::InvalidInput;
    }

    // Compute scaling factors
    Eigen::VectorXd scale;
    if (options.enableScaling) {
        scale = computeScalingFactors(x);
    } else {
        scale = Eigen::VectorXd::Ones(n);
    }

    // Work in scaled coordinates: y = x ./ scale
    Eigen::VectorXd y = x.cwiseQuotient(scale);

    Eigen::VectorXd F(n);
    Eigen::MatrixXd J(n, n);
    Eigen::MatrixXd J_unscaled(n, n);
    Eigen::VectorXd x_unscaled(n);

    double lambda = options.lmInitialLambda;
    double initialResidualNorm = 0.0;

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        if (TimeoutGuard::hasTimedOut()) {
            if (detailedError) *detailedError = "LM solver timed out";
            x = y.cwiseProduct(scale);
            return SolverStatus::EvaluationError;
        }

        // Evaluate F(x) and J(x)
        try {
            x_unscaled = y.cwiseProduct(scale);
            problem.evaluate(x_unscaled, F, J_unscaled, true);
            // Scale Jacobian: dF/dy = J_unscaled * diag(scale)
            J = J_unscaled * scale.asDiagonal();
        } catch (const std::exception& e) {
            if (options.verbose) {
                std::cerr << "LM: Evaluation failed at iter " << iter
                          << ": " << e.what() << std::endl;
            }
            throw;
        }

        double residualNorm = F.lpNorm<Eigen::Infinity>();
        if (iter == 0) initialResidualNorm = residualNorm;

        if (options.verbose) {
            std::cout << "LM iter " << iter << ": ||F||_inf = " << residualNorm
                      << ", lambda = " << lambda << std::endl;
        }

        // Record trace
        if (trace) {
            SolverTrace::Iteration traceIter;
            traceIter.iter = iter;
            traceIter.residualNorm = residualNorm;
            traceIter.stepNorm = 0.0;
            traceIter.lambda = lambda;
            Eigen::VectorXd xu = y.cwiseProduct(scale);
            traceIter.x = std::vector<double>(xu.data(), xu.data() + xu.size());
            traceIter.residuals = std::vector<double>(F.data(), F.data() + F.size());
            trace->iterations.push_back(traceIter);
        }

        // Check convergence
        if (residualNorm < options.tolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }
        if (initialResidualNorm > 0 &&
            residualNorm / initialResidualNorm < options.relativeTolerance) {
            if (trace) {
                trace->finalStatus = SolverStatus::Success;
                trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
            }
            x = y.cwiseProduct(scale);
            return SolverStatus::Success;
        }

        // Solve (J^T J + lambda * I) dy = -J^T F
        Eigen::MatrixXd JtJ = J.transpose() * J;
        Eigen::VectorXd JtF = J.transpose() * F;

        // Add damping: (J^T J + lambda * diag(J^T J)) dy = -J^T F
        // Using the diagonal scaling variant (Marquardt's improvement)
        Eigen::VectorXd diag_JtJ = JtJ.diagonal();
        for (int i = 0; i < n; ++i) {
            JtJ(i, i) += lambda * std::max(diag_JtJ(i), 1e-6);
        }

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(JtJ);
        Eigen::VectorXd dy = qr.solve(-JtF);

        if (!dy.allFinite()) {
            if (options.verbose) {
                std::cerr << "LM: Invalid step (NaN/Inf), increasing lambda" << std::endl;
            }
            lambda *= options.lmLambdaIncrease;
            if (lambda > options.lmMaxLambda) {
                if (trace) {
                    trace->finalStatus = SolverStatus::SingularJacobian;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                x = y.cwiseProduct(scale);
                return SolverStatus::SingularJacobian;
            }
            continue;
        }

        // Evaluate at trial point
        Eigen::VectorXd y_new = y + dy;
        Eigen::VectorXd F_new(n);
        Eigen::MatrixXd J_dummy;
        double phi_old = 0.5 * F.squaredNorm();
        double phi_new;

        try {
            Eigen::VectorXd x_trial = y_new.cwiseProduct(scale);
            problem.evaluate(x_trial, F_new, J_dummy, false);
            phi_new = 0.5 * F_new.squaredNorm();
        } catch (...) {
            // Trial point failed; increase damping and retry
            lambda *= options.lmLambdaIncrease;
            if (lambda > options.lmMaxLambda) {
                if (trace) {
                    trace->finalStatus = SolverStatus::EvaluationError;
                    trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                }
                x = y.cwiseProduct(scale);
                return SolverStatus::EvaluationError;
            }
            continue;
        }

        // Compute gain ratio: actual reduction / predicted reduction
        // For LM with (J^T J + lambda*D) dy = -J^T F, the predicted reduction is:
        //   predicted = 0.5 * dy^T * (lambda * D * dy - J^T * F)
        // where D = diag(max(diag(J^T J), 1e-6))
        Eigen::VectorXd D = diag_JtJ.cwiseMax(1e-6);
        double predicted = 0.5 * dy.dot(lambda * D.cwiseProduct(dy) - JtF);
        double actual = phi_old - phi_new;

        double rho = (std::abs(predicted) > 1e-30) ? actual / predicted : 0.0;

        if (options.verbose) {
            std::cout << "LM: rho = " << rho << ", actual = " << actual
                      << ", predicted = " << predicted << std::endl;
        }

        if (actual > 0) {
            // Accept step
            y = y_new;

            double stepNorm = dy.lpNorm<Eigen::Infinity>();
            if (trace && !trace->iterations.empty()) {
                trace->iterations.back().stepNorm = stepNorm;
            }

            // Decrease lambda (move towards Gauss-Newton)
            if (rho > 0.75) {
                lambda = std::max(lambda * options.lmLambdaDecrease, options.lmMinLambda);
            } else if (rho < 0.25) {
                lambda = std::min(lambda * options.lmLambdaIncrease, options.lmMaxLambda);
            }

            // Check for tiny step
            if (stepNorm < options.stepTolerance) {
                double newNorm = F_new.lpNorm<Eigen::Infinity>();
                if (newNorm < options.tolerance * 100 || newNorm < options.lsRelaxedTolerance) {
                    if (trace) {
                        trace->finalStatus = SolverStatus::Success;
                        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
                    }
                    x = y.cwiseProduct(scale);
                    return SolverStatus::Success;
                }
            }
        } else {
            // Reject step; increase lambda (move towards gradient descent)
            lambda = std::min(lambda * options.lmLambdaIncrease, options.lmMaxLambda);
            if (options.verbose) {
                std::cout << "LM: Rejecting step, lambda -> " << lambda << std::endl;
            }
        }
    }

    if (trace) {
        trace->finalStatus = SolverStatus::MaxIterations;
        trace->totalTime = std::chrono::high_resolution_clock::now() - startTime;
    }
    if (detailedError) {
        std::ostringstream ss;
        ss << "Levenberg-Marquardt: Max iterations (" << options.maxIterations
           << ") reached without convergence.";
        *detailedError = ss.str();
    }
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
            if (outErrorMessage) {
                *outErrorMessage = "Partitioned solver skipped: block too small";
            }
            return SolverStatus::InvalidInput;
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
    if (strategy == SolverStrategy::TrustRegion && solverOpts.maxIterations < 500) {
        solverOpts.maxIterations = 500;
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

    for (size_t idx = 0; idx < options.solverPipeline.size(); ++idx) {
        SolverStrategy strategy = options.solverPipeline[idx];

        if (options.verbose) {
            std::cout << "Block " << blockIndex << " (size " << n
                      << "): Trying " << strategyToString(strategy)
                      << " [" << (idx + 1) << "/" << options.solverPipeline.size() << "]"
                      << std::endl;
        }

        // Reset initial guess from evaluator state before each attempt
        if (idx > 0) {
            for (size_t i = 0; i < n; ++i) {
                x[i] = evaluator_.getVariableValue(varNames[i]);
            }
        }

        std::string solverError;
        status = runSolverStrategy(strategy, blockIndex, problem, blockEval,
                                   varNames, externalVars, externalStringVars,
                                   x, options, trace, &solverError);

        if (status == SolverStatus::Success) {
            return status;
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
    }

    // All solvers failed
    if (outErrorMessage && !lastError.empty()) {
        if (!outErrorMessage->empty()) *outErrorMessage += "\n";
        *outErrorMessage += lastError;
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

    // Shared flag: set to true when any thread succeeds
    auto winnerFound = std::make_shared<std::atomic<bool>>(false);

    for (const auto& strategy : pipeline) {
        // Each thread gets its own copy of x
        Eigen::VectorXd x_copy = x;

        // Each thread needs its own Problem with its own capture of x_copy
        // We create a new problem lambda per thread
        futures.push_back(std::async(std::launch::async,
            [this, strategy, blockIndex, &blockEval, &varNames,
             &externalVars, &externalStringVars, &options,
             x_copy, winnerFound]() mutable -> ThreadResult {
                ThreadResult result;
                result.strategy = strategy;
                result.solution = x_copy;

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
                                                  result.solution, options,
                                                  &result.trace, &result.error);

                if (result.status == SolverStatus::Success) {
                    winnerFound->store(true, std::memory_order_release);
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
