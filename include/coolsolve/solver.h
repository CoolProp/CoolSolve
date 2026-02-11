#pragma once

#include "evaluator.h"
#include "structural_analysis.h"
#include <Eigen/Dense>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <chrono>

namespace coolsolve {

// ============================================================================
// Solver Status & Options
// ============================================================================

/**
 * @brief Status codes for solver operations.
 */
enum class SolverStatus {
    Success,           // Converged to solution
    MaxIterations,     // Reached max iterations without converging
    LineSearchFailed,  // Line search couldn't find a descent direction
    SingularJacobian,  // Jacobian is singular or near-singular
    InvalidInput,      // Invalid input (e.g., wrong dimensions)
    EvaluationError,   // Error during function evaluation
    Diverged           // Solution diverged (residual increasing)
};

/**
 * @brief Convert SolverStatus to string for logging.
 */
std::string statusToString(SolverStatus status);

/**
 * @brief Options for non-linear solvers.
 */
struct SolverOptions {
    int maxIterations = 100;          // Maximum Newton iterations
    double tolerance = 1e-9;          // Convergence tolerance (||F||_inf)
    double relativeTolerance = 1e-9;  // Relative tolerance (||F|| / ||F0||)
    double stepTolerance = 1e-12;     // Minimum step size
    bool verbose = false;             // Print iteration info
    
    // Line search options
    double lsAlpha = 1e-4;            // Armijo condition parameter
    double lsRho = 0.5;               // Step reduction factor
    int lsMaxIterations = 20;         // Max line search iterations
    double lsMinStep = 1e-10;         // Minimum step size in line search
    double lsRelaxedTolerance = 1e-2; // Accept as converged when ||F|| < this (line search fail or max iter)

    // Variable scaling
    bool enableScaling = true;        // Enable automatic variable scaling for improved conditioning

    // Performance and safety
    int timeoutSeconds = 0;           // Timeout in seconds (0 = none)
};

/**
 * @brief Categories for solver/evaluator errors.
 */
enum class ErrorCategory {
    None,
    Timeout,
    Converged,
    MaxIterations,
    LineSearchFailed,
    SingularJacobian,
    CoolPropError,
    UndefinedVariable,
    UnsupportedFunction,
    EvaluationError,
    Other
};

/**
 * @brief Categorize an error message into a high-level category.
 */
ErrorCategory categorizeError(const std::string& errorMsg);

/**
 * @brief Convert ErrorCategory to string.
 */
std::string categoryToString(ErrorCategory category);

// ============================================================================
// Solver Trace (Debug Information)
// ============================================================================

/**
 * @brief Records the iteration history for debugging.
 */
struct SolverTrace {
    struct Iteration {
        int iter;
        double residualNorm;
        double stepNorm;
        double lambda;  // Line search step size
        std::vector<double> x;
        std::vector<double> residuals;
    };
    
    std::vector<Iteration> iterations;
    SolverStatus finalStatus;
    std::chrono::duration<double> totalTime;
    
    std::string toString() const;
};

// ============================================================================
// Non-Linear Solver Interface (Strategy Pattern)
// ============================================================================

/**
 * @brief Abstract base class for non-linear solvers.
 * 
 * Defines the interface for solving F(x) = 0 where F: R^n -> R^n.
 */
class NonLinearSolver {
public:
    /**
     * @brief Problem definition for the solver.
     */
    struct Problem {
        /**
         * @brief Callback to evaluate F(x) and optionally J(x).
         * 
         * @param x Current state vector
         * @param F Output: Residual vector F(x)
         * @param J Output: Jacobian matrix (if not null)
         * @param computeJacobian If true, compute and store Jacobian in J
         */
        std::function<void(const Eigen::VectorXd& x, 
                          Eigen::VectorXd& F, 
                          Eigen::MatrixXd& J, 
                          bool computeJacobian)> evaluate;
        int size;  // Number of equations/variables
    };

    virtual ~NonLinearSolver() = default;

    /**
     * @brief Solve the non-linear system.
     * 
     * @param problem Problem definition
     * @param x_guess Initial guess (modified in place to contain solution)
     * @param options Solver options
     * @param trace Optional trace for debugging (nullptr to disable)
     * @return SolverStatus indicating success or failure mode
     */
    virtual SolverStatus solve(Problem& problem, 
                               Eigen::VectorXd& x_guess, 
                               const SolverOptions& options = SolverOptions(),
                               SolverTrace* trace = nullptr,
                               std::string* detailedError = nullptr) = 0;
};

// ============================================================================
// Newton Solver (Damped Newton-Raphson with Line Search)
// ============================================================================

/**
 * @brief Damped Newton-Raphson solver with backtracking line search.
 * 
 * Algorithm:
 * 1. Compute F(x), J(x)
 * 2. Check convergence: ||F||_inf < tol
 * 3. Solve J * dx = -F for Newton step dx
 * 4. Line search: find lambda in (0, 1] such that ||F(x + lambda*dx)|| < ||F(x)||
 * 5. Update x <- x + lambda * dx
 * 6. Repeat until convergence or max iterations
 */
class NewtonSolver : public NonLinearSolver {
public:
    SolverStatus solve(Problem& problem, 
                       Eigen::VectorXd& x_guess, 
                       const SolverOptions& options = SolverOptions(),
                       SolverTrace* trace = nullptr,
                       std::string* detailedError = nullptr) override;
    
private:
    /**
     * @brief Compute automatic scaling factors for variables.
     *
     * Variables with different magnitudes (T~300, P~1e7, h~3e5) create
     * ill-conditioned Jacobians. Scaling improves convergence.
     *
     * @param x Initial guess vector
     * @return Vector of scaling factors (one per variable)
     */
    Eigen::VectorXd computeScalingFactors(const Eigen::VectorXd& x) const;
    
    /**
     * @brief Perform backtracking line search.
     *
     * Finds lambda such that phi(lambda) = ||F(x + lambda*dx)||^2 < phi(0)
     * using the Armijo condition.
     *
     * @param problem Problem definition
     * @param x Current point
     * @param dx Newton direction
     * @param F Current residuals F(x)
     * @param options Solver options
     * @return Step size lambda, or 0 if line search failed
     */
    double lineSearch(Problem& problem,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& dx,
                      const Eigen::VectorXd& F,
                      const SolverOptions& options);
};

// ============================================================================
// Main Solver (Orchestrator)
// ============================================================================

/**
 * @brief Result of solving the complete system.
 */
struct SolveResult {
    bool success = false;
    SolverStatus status = SolverStatus::InvalidInput;
    std::string errorMessage;
    
    // Solution values
    std::map<std::string, double, CaseInsensitiveLess> variables;
    std::map<std::string, std::string, CaseInsensitiveLess> stringVariables;
    
    // Statistics
    int totalIterations = 0;
    int blocksEvaluated = 0;
    std::chrono::duration<double> totalTime{0};
    
    // Per-block results
    struct BlockResult {
        size_t id;
        bool success;
        SolverStatus status;
        int iterations;
        double maxResidual;
        std::string errorMessage;
    };
    std::vector<BlockResult> blockResults;
    
    // Per-block traces (for debugging)
    std::vector<SolverTrace> blockTraces;
    
    // Detailed error from evaluation or solver
    std::string detailedError;
};

/**
 * @brief Main solver class that orchestrates the solution of the entire system.
 * 
 * Iterates through blocks in topological order:
 * - For size-1 blocks: Attempts direct evaluation, falls back to Newton if implicit
 * - For larger blocks: Uses Newton solver
 */
class Solver {
public:
    /**
     * @brief Construct a solver with the given system.
     * 
     * @param ir The intermediate representation of the equation system
     * @param analysis The structural analysis result
     * @param config CoolProp configuration
     */
    Solver(const IR& ir, 
           const StructuralAnalysisResult& analysis,
           const CoolPropConfig& config = CoolPropConfig());
    
    /**
     * @brief Solve the complete system.
     * 
     * @param options Solver options
     * @param enableTracing If true, record iteration traces for each block
     * @return SolveResult containing solution and status
     */
    SolveResult solve(const SolverOptions& options = SolverOptions(),
                      bool enableTracing = false);
    
    /**
     * @brief Set initial guess for a variable.
     */
    void setGuess(const std::string& name, double value);
    
    /**
     * @brief Set a string variable value.
     */
    void setStringVariable(const std::string& name, const std::string& value);
    
    /**
     * @brief Get the system evaluator (for inspection).
     */
    const SystemEvaluator& getEvaluator() const { return evaluator_; }
    SystemEvaluator& getEvaluator() { return evaluator_; }
    
private:
    SystemEvaluator evaluator_;
    const IR& ir_;
    const StructuralAnalysisResult& analysis_;
    
    /**
     * @brief Solve a single block.
     * 
     * @param blockIndex Index of the block to solve
     * @param options Solver options
     * @param trace Optional trace for debugging
     * @param outErrorMessage Optional output for detailed error message
     * @return Status of solving this block
     */
    SolverStatus solveBlock(size_t blockIndex, 
                           const SolverOptions& options,
                           SolverTrace* trace,
                           std::string* outErrorMessage = nullptr);
    
    /**
     * @brief Try to solve an explicit block directly.
     * 
     * For size-1 blocks where the equation can be rearranged to x = expr(other_vars).
     * 
     * @return true if solved directly, false if Newton iteration needed
     */
    bool tryExplicitSolve(size_t blockIndex);
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Generate a summary report of the solve process.
 */
std::string generateSolveReport(const SolveResult& result);

/**
 * @brief RAII class to handle timeouts using SIGALRM.
 */
class TimeoutGuard {
public:
    TimeoutGuard(int seconds);
    ~TimeoutGuard();
    static bool hasTimedOut();
private:
    int seconds_;
};

}  // namespace coolsolve
