#pragma once

#include "evaluator.h"
#include "structural_analysis.h"
#include <Eigen/Dense>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <future>

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
    Diverged,          // Solution diverged (residual increasing)
    ParseFailed        // Input failed to parse
};

/**
 * @brief Convert SolverStatus to string for logging.
 */
std::string statusToString(SolverStatus status);

/**
 * @brief Identifies a solver algorithm that can be used in the solver pipeline.
 *
 * Each value corresponds to a concrete NonLinearSolver subclass.
 * Solvers can be composed into fallback chains or run in parallel.
 */
enum class SolverStrategy {
    Newton,            ///< Damped Newton-Raphson with backtracking line search
    TrustRegion,       ///< Trust-region dogleg method
    LevenbergMarquardt,///< Levenberg-Marquardt (damped least-squares)
    Partitioned,       ///< Per-variable diagonal updates (tear-based)
    BisectionND,       ///< Multi-dimensional bisection (small systems only, n≤8)
    Homotopy           ///< Homotopy continuation (global convergence)
};

/**
 * @brief Convert SolverStrategy to a human-readable string.
 */
std::string strategyToString(SolverStrategy strategy);

/**
 * @brief Parse a solver strategy name (case-insensitive) to enum.
 * @return true if parsing succeeded, false if the name is unknown.
 */
bool parseStrategy(const std::string& name, SolverStrategy& out);

/**
 * @brief Execution mode for the solver pipeline.
 *
 * - Sequential: solvers are tried one after another (fallback chain).
 * - Parallel:   solvers are launched simultaneously in separate threads;
 *               the first one to converge wins and the others are cancelled.
 */
enum class SolverPipelineMode {
    Sequential,  ///< Try solvers in order; stop at first success (default)
    Parallel     ///< Run solvers concurrently; first success wins
};

/**
 * @brief Convert SolverPipelineMode to string.
 */
std::string pipelineModeToString(SolverPipelineMode mode);

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
    double lsRho = 0.5;              // Step reduction factor
    int lsMaxIterations = 100;        // Max line search iterations
    double lsMinStep = 1e-10;         // Minimum step size in line search
    double lsRelaxedTolerance = 1e-2; // Accept as converged when ||F|| < this (line search fail or max iter)

    // Non-monotone line search (Grippo, Lampariello, Lucidi 1986).
    // Instead of requiring strict decrease at every step (monotone Armijo),
    // compare against the maximum merit value over the last M iterations.
    // This helps escape narrow curved valleys and saddle points that trap
    // monotone methods.  M=1 gives standard monotone behavior.
    // Applied to Newton (Armijo), TrustRegion (acceptance), and LM (acceptance).
    int lsNonMonotoneMemory = 10;     // Non-monotone memory M (default 10; 1 = monotone)

    // Variable scaling
    bool enableScaling = true;        // Enable automatic variable scaling for improved conditioning

    // --- Broyden quasi-Newton (Jacobian reuse) ---
    // When > 0, the Newton solver computes a full Jacobian every K iterations
    // and uses Broyden rank-1 updates in between, saving O(n) CoolProp calls
    // per intermediate iteration.  If a Broyden step fails the line search,
    // a full Jacobian is automatically recomputed.  0 = disabled (default).
    int broydenRecomputeInterval = 0;

    // Trust region options
    double trInitialRadius = 10.0;    // Initial trust region radius (larger for aggressive steps)
    double trMaxRadius = 1000.0;      // Maximum trust region radius
    double trEta = 0.05;              // Threshold for accepting step (rho >= eta) - more lenient
    double trShrinkFactor = 0.5;      // Factor to shrink trust region on rejection (less aggressive)
    double trGrowFactor = 2.0;        // Factor to grow trust region on good steps
    bool trAdaptiveRadius = true;     // Adaptive initial radius based on Cauchy step norm

    // Levenberg-Marquardt options
    double lmInitialLambda = 1e-3;    // Initial damping parameter
    double lmLambdaIncrease = 10.0;   // Factor to increase lambda on bad step (fallback; see Nielsen)
    double lmLambdaDecrease = 0.1;    // Factor to decrease lambda on good step (fallback; see Nielsen)
    double lmMinLambda = 1e-12;       // Minimum damping parameter
    double lmMaxLambda = 1e8;         // Maximum damping parameter
    bool lmNielsenUpdate = true;      // Use Nielsen's smooth lambda adaptation (default: true)
    bool lmGeodesicAcceleration = true; // Add geodesic acceleration correction (1 extra eval/iter)

    // Partitioned block solver options (variable-wise updates)
    int partitionedMaxIterations = 300;   // Max iterations for partitioned solver
    double partitionedRelaxation = 0.6;   // Relaxation factor for updates (0 < w <= 1)
    double partitionedMinDiagonal = 1e-12; // Minimum |dF_i/dx_i| to update variable
    int partitionedMinBlockSize = 4;      // Only use partitioned solver for blocks >= this size

    // --- Tearing (structural decomposition of algebraic loops) ---
    // When enabled, blocks above tearingMinBlockSize are solved by selecting a tear set
    // (feedback vertex set), solving the acyclic part sequentially, then Newton on tear residuals.
    bool enableTearing = false;           // Use tearing for large blocks when true
    int tearingMaxIterations = 100;      // Max outer (tear) Newton iterations
    int tearingMinBlockSize = 3;         // Only apply tearing to blocks of this size or larger
    int tearingInnerIterations = 5;      // Max 1D iterations per equation in acyclic solve

    // --- Symbolic Block Reduction ---
    // When enabled, blocks are pre-processed to reduce their size before the
    // iterative solver runs.  Techniques: explicit extraction, CoolProp call
    // inversion, and equation substitution.  Off by default (zero overhead).
    bool enableSymbolicReduction = false;

    // --- BisectionND options ---
    // BisectionND is a derivative-free sign-change bisection solver.
    // It is only feasible for small blocks because it requires 2^n function evaluations
    // per probe phase, and the iteration cost also grows exponentially with n.
    // Blocks larger than bisectionNDMaxBlockSize are automatically skipped (InvalidInput).
    int bisectionNDMaxBlockSize = 8;   // Skip BisectionND for blocks with more than this many unknowns.
                                       // Default: 8. Increase with caution: cost is exponential in n.
    // Multiplier applied to maxIterations for the BisectionND bisection loop.
    // Because bisection converges slowly (linear), more iterations than Newton-type
    // solvers are often needed.  bisectionNDMaxBlockSize is unaffected by this factor.
    // Example: bisectionNDIterFactor = 5 gives 500 bisection steps when maxIterations = 100.
    double bisectionNDIterFactor = 1.0; // Multiplier for max bisection iterations (default: 1.0).

    // --- Solver pipeline configuration ---
    // The pipeline defines which solvers to try and in what order.
    // Default: Newton -> TrustRegion -> LM -> BisectionND (small blocks only)
    //       -> Homotopy -> Partitioned.
    // BisectionND automatically returns InvalidInput for blocks > bisectionNDMaxBlockSize and is skipped.
    std::vector<SolverStrategy> solverPipeline = {
        SolverStrategy::Newton,
        SolverStrategy::TrustRegion,
        SolverStrategy::LevenbergMarquardt,
        SolverStrategy::BisectionND,
        SolverStrategy::Homotopy,
        SolverStrategy::Partitioned
    };

    /// Execution mode: sequential fallback or parallel (first-to-solve wins)
    SolverPipelineMode pipelineMode = SolverPipelineMode::Sequential;

    // Performance and safety
    int timeoutSeconds = 0;           // Timeout in seconds (0 = none)

    // Progress callback (optional, for GUI/SSE reporting)
    // Called at the start/end of each block solve.
    // Parameters: blockIndex, totalBlocks, event ("start"/"done"/"fail"), iterations, residualNorm
    using ProgressCallback = std::function<void(
        int blockIndex, int totalBlocks,
        const std::string& event,
        int iterations, double residualNorm
    )>;
    ProgressCallback progressCallback = nullptr;

    // Cancellation token (optional, for GUI stop button)
    // When set, the solver checks this between blocks and aborts if cancelled.
    std::atomic<bool>* cancelToken = nullptr;
    
    // CoolProp integration options
    CoolPropConfig coolpropConfig;

    // --- LaTeX Report ---
    // When true, a comprehensive LaTeX report is generated after a successful
    // solve.  The report includes model equations, variable solutions, block
    // structure, and solver statistics.  In debug mode (-d) this is always
    // enabled.  In non-debug mode the .tex file is written next to the model.
    // Configurable via coolsolve.conf: enableLatexReport = true/false
    bool enableLatexReport = true;

    // LaTeX compiler command (default: pdflatex).  Used by the GUI to compile
    // the generated .tex file into a PDF report.
    std::string latexCompiler = "pdflatex";

    // --- Debug output ---
    // When non-empty, the solver writes a Markdown file (one per solve) showing
    // original vs reduced equations for every block where symbolic reduction
    // was applied.  Useful for inspecting what the reduction pass did.
    std::string debugReductionPath;
};

/**
 * @brief Load solver options from a coolsolve.conf file.
 * Only keys present in the file are applied; others keep their current values.
 * @param path Path to the config file (e.g. dir/coolsolve.conf)
 * @param options Options to override (modified in place)
 * @return true if the file was found and read, false otherwise
 */
bool loadSolverOptionsFromFile(const std::string& path, SolverOptions& options);

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
// Solver Attempt (Per-Solver Result in Pipeline)
// ============================================================================

/**
 * @brief Result of a single solver strategy attempt on one block.
 */
struct SolverAttempt {
    SolverStrategy strategy;
    SolverStatus status;
    int iterations = 0;
    double finalResidual = 0.0;
    double elapsedMs = 0.0;
};

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
        std::string detail;  // Optional per-iteration detail (e.g. inner solver info for tearing)
    };
    
    std::vector<Iteration> iterations;
    SolverStatus finalStatus;
    std::chrono::duration<double> totalTime;
    std::string solverType;  // "Newton", "TrustRegion", "Tearing", etc.

    // Variable names for enriched trace output (especially for tearing)
    std::vector<std::string> varNames;
    // Tear variable names (only populated for tearing solver)
    std::vector<std::string> tearVarNames;

    // Singular Jacobian diagnostics (only populated when finalStatus == SingularJacobian)
    std::vector<double> singularJacobianF;           // Residual vector F at failure
    std::vector<std::vector<double>> singularJacobianJ;  // Jacobian matrix J at failure

    // Per-solver attempt results recorded during sequential pipeline execution
    std::vector<SolverAttempt> solverAttempts;

    // Symbolic reduction info recorded by solveBlock when reduction was applied
    bool symbolicReductionApplied = false;
    int originalBlockSize = 0;
    int reducedBlockSize = 0;
    int symInversions = 0;
    int symExtractions = 0;
    int symSubstitutions = 0;
    /// Original equation texts (before reduction)
    std::vector<std::string> originalEquations;
    /// Descriptions of reduction steps applied
    std::vector<std::string> reductionStepDescriptions;
    /// Remaining equation texts (after reduction)
    std::vector<std::string> reducedEquations;

    // Re-decomposition info (when reduced block splits into sub-SCCs)
    bool redecompositionApplied = false;
    int numSubBlocks = 0;
    std::vector<int> subBlockSizes;   // size of each sub-block after re-decomposition
    
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
     * @brief Perform backtracking line search with non-monotone Armijo condition.
     *
     * Finds lambda such that φ(x + λ dx) ≤ refPhi + α λ ∇φ·dx,
     * where refPhi is the non-monotone reference (max of recent merit values).
     * When refPhi == φ(x), this is the standard monotone Armijo condition.
     *
     * @param problem Problem definition
     * @param x Current point
     * @param dx Newton direction
     * @param F Current residuals F(x)
     * @param options Solver options
     * @param refPhi Non-monotone reference merit value (max of recent history)
     * @return Step size lambda, or 0 if line search failed
     */
    double lineSearch(Problem& problem,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& dx,
                      const Eigen::VectorXd& F,
                      const SolverOptions& options,
                      double refPhi);
};

// ============================================================================
// Trust Region Dogleg Solver
// ============================================================================

/**
 * @brief Trust region solver using the dogleg method.
 *
 * More robust than line search for highly nonlinear problems.
 * Adapts step size based on model accuracy and switches between
 * Newton and gradient descent directions as needed.
 *
 * Algorithm:
 * 1. Compute Newton step dx_n by solving J*dx_n = -F
 * 2. Compute Cauchy (steepest descent) step dx_c
 * 3. If ||dx_n|| <= Delta, use Newton step
 * 4. If ||dx_c|| >= Delta, use scaled Cauchy step
 * 5. Otherwise, use dogleg path interpolation
 * 6. Adapt trust radius Delta based on actual vs predicted reduction
 */
class TrustRegionSolver : public NonLinearSolver {
public:
    SolverStatus solve(Problem& problem,
                       Eigen::VectorXd& x_guess,
                       const SolverOptions& options = SolverOptions(),
                       SolverTrace* trace = nullptr,
                       std::string* detailedError = nullptr) override;
    
private:
    /**
     * @brief Compute automatic scaling factors for variables.
     */
    Eigen::VectorXd computeScalingFactors(const Eigen::VectorXd& x) const;
    
    /**
     * @brief Compute the dogleg step within the trust region.
     *
     * @param dx_n Newton step
     * @param dx_c Cauchy step (steepest descent)
     * @param delta Trust region radius
     * @return The dogleg step
     */
    Eigen::VectorXd doglegStep(const Eigen::VectorXd& dx_n,
                               const Eigen::VectorXd& dx_c,
                               double delta);
    
    /**
     * @brief Evaluate the model at a proposed step.
     *
     * m(p) = 0.5 * ||F + J*p||^2
     *
     * @param F Current residual
     * @param J Current Jacobian
     * @param p Proposed step
     * @return Model value
     */
    double evaluateModel(const Eigen::VectorXd& F,
                         const Eigen::MatrixXd& J,
                         const Eigen::VectorXd& p);
};

// ============================================================================
// Levenberg-Marquardt Solver
// ============================================================================

/**
 * @brief Levenberg-Marquardt solver for nonlinear least-squares problems.
 *
 * Blends Gauss-Newton and gradient descent via an adaptive damping parameter λ.
 * Particularly effective when the initial guess is far from the solution,
 * because the damping prevents oversized steps that would cause divergence.
 *
 * Solves: min 0.5*||F(x)||^2  by iterating
 *   (J^T J + λ I) dx = -J^T F
 *
 * - When λ is large → gradient descent (safe, slow)
 * - When λ is small → Gauss-Newton (fast, quadratic convergence near solution)
 * - λ is adapted based on actual vs predicted reduction (gain ratio).
 */
class LevenbergMarquardtSolver : public NonLinearSolver {
public:
    SolverStatus solve(Problem& problem,
                       Eigen::VectorXd& x_guess,
                       const SolverOptions& options = SolverOptions(),
                       SolverTrace* trace = nullptr,
                       std::string* detailedError = nullptr) override;

private:
    Eigen::VectorXd computeScalingFactors(const Eigen::VectorXd& x) const;
};

// ============================================================================
// Multi-Dimensional Bisection Solver (small systems only)
// ============================================================================

/**
 * @brief Multi-dimensional bisection for small nonlinear systems (n ≤ 8).
 *
 * Maintains a simplex of vertices with diverse residual sign patterns and
 * iteratively bisects the longest edge.  Guaranteed to converge if a simplex
 * with the right sign structure is found.  Only used as a last resort for
 * small blocks that resist Newton-type methods.
 */
class BisectionNDSolver : public NonLinearSolver {
public:
    SolverStatus solve(Problem& problem,
                       Eigen::VectorXd& x_guess,
                       const SolverOptions& options = SolverOptions(),
                       SolverTrace* trace = nullptr,
                       std::string* detailedError = nullptr) override;
};

// ============================================================================
// Homotopy Continuation Solver
// ============================================================================

/**
 * @brief Homotopy continuation solver for global convergence.
 *
 * Constructs H(x,t) = t·F(x) + (1-t)·(x-x0) and tracks the solution
 * from t=0 (trivial) to t=1 (target) using an existing solver as corrector.
 * Adaptive step-size control: increase dt on success, decrease on failure.
 */
class HomotopySolver : public NonLinearSolver {
public:
    SolverStatus solve(Problem& problem,
                       Eigen::VectorXd& x_guess,
                       const SolverOptions& options = SolverOptions(),
                       SolverTrace* trace = nullptr,
                       std::string* detailedError = nullptr) override;
};

// ============================================================================
// Specialized 1D Root-Finding Solver
// ============================================================================

/**
 * @brief Multi-phase 1D root-finding solver for size-1 implicit blocks.
 *
 * Unlike the n-dimensional NonLinearSolver hierarchy, this solver works
 * directly on a scalar evaluation function (residual + derivative).
 * It uses a four-phase approach:
 *   - Phase 1: Trust-region Newton with bracket detection
 *   - Phase 2: Multi-probe sign-change exploration (~900 probes)
 *   - Phase 3: Bisection + Newton hybrid within bracket
 *   - Phase 4: Final Newton polish with relaxed tolerance
 *
 * A simplified solver (solveSimple) is also provided for use in symbolic
 * reduction where the problem is already well-conditioned.
 */
class Newton1DSolver {
public:
    /// 1D evaluation function: given x, returns (residual, derivative).
    using Eval1D = std::function<std::pair<double, double>(double)>;

    /**
     * @brief Full multi-phase solve with probing and bracket detection.
     *
     * @param eval            Evaluation function f(x) → (residual, df/dx)
     * @param x               [in/out] Initial guess → converged solution
     * @param externalVars    External variable values (used for probe generation)
     * @param options         Solver options (tolerance, maxIterations, etc.)
     * @param trace           Optional trace for iteration recording
     * @param outErrorMessage Optional detailed error message on failure
     * @return SolverStatus::Success if converged, MaxIterations or EvaluationError otherwise
     */
    static SolverStatus solve(
        Eval1D& eval,
        double& x,
        const std::map<std::string, double>& externalVars,
        const SolverOptions& options,
        SolverTrace* trace = nullptr,
        std::string* outErrorMessage = nullptr);

    /**
     * @brief Simplified Newton solver for well-conditioned 1D problems.
     *
     * Basic Newton iteration without probing or bracket detection.
     * Used inside symbolic reduction where the reduced block is already
     * close to the solution.
     *
     * @param eval    Evaluation function f(x) → (residual, df/dx)
     * @param x       [in/out] Initial guess → converged solution
     * @param options Solver options (tolerance, maxIterations)
     * @return SolverStatus::Success if converged
     */
    static SolverStatus solveSimple(
        Eval1D& eval,
        double& x,
        const SolverOptions& options);
};

// ============================================================================
// Solver Factory
// ============================================================================

/**
 * @brief Create a NonLinearSolver instance for the given strategy.
 *
 * @param strategy The solver algorithm to instantiate
 * @return A unique_ptr to the solver (never null)
 */
std::unique_ptr<NonLinearSolver> createSolver(SolverStrategy strategy);

// ============================================================================
// Cancellation Token (for parallel solver execution)
// ============================================================================

/**
 * @brief Thread-safe cancellation token used to stop parallel solvers
 *        once one of them has converged.
 */
class CancellationToken {
public:
    CancellationToken() : cancelled_(false) {}

    /// Request cancellation (thread-safe).
    void cancel() { cancelled_.store(true, std::memory_order_release); }

    /// Check whether cancellation has been requested (thread-safe).
    bool isCancelled() const { return cancelled_.load(std::memory_order_acquire); }

private:
    std::atomic<bool> cancelled_;
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

        // Symbolic reduction info
        int originalSize = 0;             ///< Original block size (number of variables)
        int reducedSize = 0;              ///< Size after symbolic reduction (== originalSize if not applied)
        bool symbolicReductionApplied = false;
        int inversionsApplied = 0;
        int extractionsApplied = 0;
        int substitutionsApplied = 0;

        // Per-solver attempt results (only populated when tracing is enabled)
        std::vector<SolverAttempt> solverAttempts;

        // Re-decomposition info
        bool redecompositionApplied = false;
        int numSubBlocks = 0;
        std::vector<int> subBlockSizes;
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
 * - For larger blocks: Uses the configured solver pipeline (fallback chain or parallel)
 *
 * The solver pipeline is configured via SolverOptions::solverPipeline and
 * SolverOptions::pipelineMode.  In **sequential** mode the solvers are tried
 * one after another; in **parallel** mode they are launched concurrently and
 * the first one to converge wins.
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
     * @brief Solve a single block using the configured solver pipeline.
     *
     * Dispatches to solveBlockSequential() or solveBlockParallel() depending
     * on the pipeline mode.
     *
     * @param blockIndex Index of the block to solve
     * @param options Solver options (includes pipeline config)
     * @param trace Optional trace for debugging
     * @param outErrorMessage Optional output for detailed error message
     * @return Status of solving this block
     */
    SolverStatus solveBlock(size_t blockIndex,
                           const SolverOptions& options,
                           SolverTrace* trace,
                           std::string* outErrorMessage = nullptr);

    /**
     * @brief Run a single NonLinearSolver strategy on a block.
     *
     * Handles the Partitioned strategy specially (it needs the structural
     * matching information).  All other strategies go through the
     * NonLinearSolver interface.
     *
     * @return SolverStatus from the chosen strategy.
     */
    SolverStatus runSolverStrategy(SolverStrategy strategy,
                                   size_t blockIndex,
                                   NonLinearSolver::Problem& problem,
                                   BlockEvaluator& blockEval,
                                   const std::vector<std::string>& varNames,
                                   const std::map<std::string, double>& externalVars,
                                   const std::map<std::string, std::string>& externalStringVars,
                                   Eigen::VectorXd& x,
                                   const SolverOptions& options,
                                   SolverTrace* trace,
                                   std::string* outErrorMessage);

    /**
     * @brief Sequential fallback: try each solver in the pipeline in order.
     */
    SolverStatus solveBlockSequential(size_t blockIndex,
                                      NonLinearSolver::Problem& problem,
                                      BlockEvaluator& blockEval,
                                      const std::vector<std::string>& varNames,
                                      const std::map<std::string, double>& externalVars,
                                      const std::map<std::string, std::string>& externalStringVars,
                                      Eigen::VectorXd& x,
                                      const SolverOptions& options,
                                      SolverTrace* trace,
                                      std::string* outErrorMessage);

    /**
     * @brief Parallel execution: launch all pipeline solvers concurrently.
     *
     * Each solver runs in its own thread with a copy of the initial guess.
     * The first solver to converge sets the solution; the others are
     * (logically) abandoned.
     */
    SolverStatus solveBlockParallel(size_t blockIndex,
                                    NonLinearSolver::Problem& problem,
                                    BlockEvaluator& blockEval,
                                    const std::vector<std::string>& varNames,
                                    const std::map<std::string, double>& externalVars,
                                    const std::map<std::string, std::string>& externalStringVars,
                                    Eigen::VectorXd& x,
                                    const SolverOptions& options,
                                    SolverTrace* trace,
                                    std::string* outErrorMessage);

    /**
     * @brief Partitioned block solve using per-equation variable updates.
     *
     * Uses the equation-to-output-variable mapping from structural matching
     * to apply diagonal Newton-like updates per variable. This is a robust
     * fallback for highly nonlinear or ill-conditioned blocks.
     */
    SolverStatus solveBlockPartitioned(size_t blockIndex,
                                       BlockEvaluator& blockEval,
                                       const std::vector<std::string>& varNames,
                                       const std::map<std::string, double>& externalVars,
                                       const std::map<std::string, std::string>& externalStringVars,
                                       Eigen::VectorXd& x,
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

    /**
     * @brief Solve a block using structural tearing (feedback vertex set + acyclic solve + Newton on tears).
     *
     * When enableTearing is true and the block is large enough, this can be tried first.
     */
    SolverStatus solveBlockTearing(size_t blockIndex,
                                   BlockEvaluator& blockEval,
                                   const std::vector<std::string>& varNames,
                                   const std::map<std::string, double>& externalVars,
                                   const std::map<std::string, std::string>& externalStringVars,
                                   Eigen::VectorXd& x,
                                   const SolverOptions& options,
                                   SolverTrace* trace,
                                   std::string* outErrorMessage = nullptr);

    /**
     * @brief Write a Markdown debug report showing original vs reduced equations
     *        for every block where symbolic reduction was applied.
     */
    static void writeDebugReductionReport(const std::string& path,
                                          const SolveResult& result);
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Generate a summary report of the solve process.
 */
std::string generateSolveReport(const SolveResult& result);

/**
 * @brief RAII class to handle timeouts.
 *
 * Uses a thread-local chrono deadline so that multiple threads can each
 * have independent timeouts (needed for parallel robustness testing).
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
