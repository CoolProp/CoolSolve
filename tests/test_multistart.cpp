/**
 * Multi-start fallback tests (roadmap §4.2).
 *
 * Covers:
 *   - A synthetic algebraic block with a singular Jacobian at the default
 *     guess (1.0) whose solution lies at a much smaller scale (0.1): multi-start
 *     must rescue it by scaling the starting point.
 *   - The option gating (disabled -> fail, enabled -> success).
 *   - The candidate-count bound.
 *   - Targeted real models (piston_compressor, refrigeration_compressor) that
 *     are known to fail without initials and to be rescued by multi-start.
 *
 * Run with: ./coolsolve_tests "[multistart]"
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/parser.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/solver.h"
#include "coolsolve/variable_inference.h"

#include <filesystem>
#include <memory>

using namespace coolsolve;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace fs = std::filesystem;

namespace {
// Solve a source snippet from scratch (no .initials) with the given options.
// Returns the SolveResult; the IR is kept alive for inspection.
struct SolveOutcome {
    bool success = false;
    SolveResult result;
    std::unique_ptr<IR> ir;
};

SolveOutcome solveFromSource(const std::string& code, SolverOptions opts) {
    SolveOutcome out;
    EESParser parser;
    auto parseResult = parser.parse(code);
    if (!parseResult.success) return out;

    out.ir = std::make_unique<IR>(IR::fromAST(parseResult.program));
    try {
        inferVariables(*out.ir);
        initializeVariables(*out.ir);
    } catch (...) {
        return out;
    }

    StructuralAnalysisResult analysis = StructuralAnalyzer::analyze(*out.ir);
    if (!analysis.success) return out;

    Solver solver(*out.ir, analysis);
    out.result = solver.solve(opts, /*enableTracing=*/true);
    out.success = out.result.success;
    return out;
}

// Default options: Newton-only pipeline so the singular-Jacobian failure is
// not masked by BisectionND/Homotopy fallbacks — this isolates multi-start.
// Multi-start is set to Always so legacy `multiStartEnabled = true/false`
// overrides below continue to work as a pure on/off switch for these tests.
SolverOptions newtonOnlyOpts() {
    SolverOptions opts;
    opts.tolerance = 1e-8;
    opts.maxIterations = 50;
    opts.solverPipeline = {SolverStrategy::Newton};
    opts.multiStartMode = MultiStartMode::Always;
    return opts;
}
} // namespace

// ----------------------------------------------------------------------------
// Synthetic block: 1/x + 1/y = 20  and  x*y = 0.01.
//   - Solution: x = y = 0.1.
//   - At the default guess (1, 1) the Jacobian [ -1/x^2, -1/y^2 ; y, x ] =
//     [ -1, -1 ; 1, 1 ] is rank-deficient, so Newton stalls.
//   - At (0.1, 0.1) both residuals are already 0.
// This is the minimal analogue of piston_compressor's volumetric-efficiency
// block, where the solution scale (~0.05-0.1) is far below the default 1.0.
// ----------------------------------------------------------------------------
static const char* kSingularScaleBlock =
    "1/x + 1/y = 20\n"
    "x*y = 0.01\n";

TEST_CASE("Multi-start rescues a block with singular Jacobian at default guess",
          "[solver][multistart]") {
    // Without multi-start: Newton-only fails (singular Jacobian / stall).
    SolverOptions optsOff = newtonOnlyOpts();
    optsOff.multiStartEnabled = false;
    auto off = solveFromSource(kSingularScaleBlock, optsOff);
    REQUIRE(off.ir != nullptr);
    REQUIRE_FALSE(off.success);

    // With multi-start: the scale x0.1 candidate lands exactly on the solution.
    SolverOptions optsOn = newtonOnlyOpts();
    optsOn.multiStartEnabled = true;
    optsOn.multiStartMaxRestarts = 4;
    auto on = solveFromSource(kSingularScaleBlock, optsOn);
    REQUIRE(on.success);

    // Verify the solution values are recovered (stored in result.variables).
    REQUIRE(on.result.variables.count("x"));
    REQUIRE(on.result.variables.count("y"));
    REQUIRE_THAT(on.result.variables.at("x"), WithinAbs(0.1, 1e-4));
    REQUIRE_THAT(on.result.variables.at("y"), WithinAbs(0.1, 1e-4));

    // The multi-start rescue must be recorded as an info diagnostic.
    bool foundDiag = false;
    for (const auto& br : on.result.blockResults) {
        for (const auto& d : br.diagnostics.items()) {
            if (d.code == "V006") foundDiag = true;
        }
    }
    REQUIRE(foundDiag);
}

TEST_CASE("Multi-start disabled leaves the singular block unsolved",
          "[solver][multistart]") {
    SolverOptions opts = newtonOnlyOpts();
    opts.multiStartEnabled = false;
    auto out = solveFromSource(kSingularScaleBlock, opts);
    REQUIRE_FALSE(out.success);
}

TEST_CASE("Multi-start respects multiStartMaxRestarts=0 (no candidates)",
          "[solver][multistart]") {
    // maxRestarts=0 produces no candidates, so the block is not retried.
    SolverOptions opts = newtonOnlyOpts();
    opts.multiStartEnabled = true;
    opts.multiStartMaxRestarts = 0;
    auto out = solveFromSource(kSingularScaleBlock, opts);
    REQUIRE_FALSE(out.success);
}

TEST_CASE("Multi-start is zero-overhead when the block converges first time",
          "[solver][multistart]") {
    // A trivially-solvable block: multi-start must not engage.
    const char* code =
        "a = 2\n"
        "b = 3\n"
        "c = a + b\n";
    SolverOptions opts = newtonOnlyOpts();
    opts.multiStartEnabled = true;
    auto out = solveFromSource(code, opts);
    REQUIRE(out.success);
    // No V006 (multi-start) diagnostic should be emitted.
    for (const auto& br : out.result.blockResults) {
        for (const auto& d : br.diagnostics.items()) {
            REQUIRE(d.code != "V006");
        }
    }
}

TEST_CASE("Multi-start parallel mode rescues the singular block",
          "[solver][multistart][parallel]") {
    // Same singular-Jacobian block as above, but candidates run concurrently.
    // The parallel path must recover the same solution (x = y = 0.1).
    SolverOptions opts = newtonOnlyOpts();
    opts.multiStartEnabled = true;
    opts.multiStartMaxRestarts = 4;
    opts.multiStartNumCores = 2;  // enable parallel multi-start
    auto out = solveFromSource(kSingularScaleBlock, opts);
    REQUIRE(out.success);
    REQUIRE(out.result.variables.count("x"));
    REQUIRE(out.result.variables.count("y"));
    REQUIRE_THAT(out.result.variables.at("x"), WithinAbs(0.1, 1e-4));
    REQUIRE_THAT(out.result.variables.at("y"), WithinAbs(0.1, 1e-4));

    bool foundDiag = false;
    for (const auto& br : out.result.blockResults) {
        for (const auto& d : br.diagnostics.items()) {
            if (d.code == "V006") foundDiag = true;
        }
    }
    REQUIRE(foundDiag);
}

TEST_CASE("Multi-start auto cores (0) rescues the singular block",
          "[solver][multistart][parallel]") {
    SolverOptions opts = newtonOnlyOpts();
    opts.multiStartEnabled = true;
    opts.multiStartMaxRestarts = 4;
    opts.multiStartNumCores = 0;  // auto = use hardware_concurrency
    auto out = solveFromSource(kSingularScaleBlock, opts);
    REQUIRE(out.success);
    REQUIRE_THAT(out.result.variables.at("x"), WithinAbs(0.1, 1e-4));
}

// ----------------------------------------------------------------------------
// Targeted real models (run only when the examples directory is available).
// These are the models the roadmap §4.2 identifies as the multi-start targets:
// they solve WITH initials but FAIL without, and the failure is a wrong-scale
// start on a small algebraic block.
// ----------------------------------------------------------------------------
TEST_CASE("Multi-start rescues piston_compressor without initials",
          "[.][solver][multistart][examples]") {
    fs::path exDir = "../examples/";
    if (const char* e = std::getenv("COOLSOLVE_EXAMPLES_DIR")) exDir = fs::path(e);
    fs::path model = exDir / "piston_compressor.eescode";
    if (!fs::exists(model)) { SUCCEED("model not found, skipping"); return; }

    EESParser parser;
    auto pr = parser.parseFile(model.string());
    REQUIRE(pr.success);
    auto ir = std::make_unique<IR>(IR::fromAST(pr.program));
    try { inferVariables(*ir); initializeVariables(*ir); } catch (...) {}
    auto an = StructuralAnalyzer::analyze(*ir);

    SolverOptions opts;
    opts.tolerance = 1e-6;
    opts.timeoutSeconds = 30;
    opts.solverPipeline = {
        SolverStrategy::Newton, SolverStrategy::TrustRegion,
        SolverStrategy::LevenbergMarquardt, SolverStrategy::BisectionND,
        SolverStrategy::Homotopy, SolverStrategy::Partitioned,
    };
    opts.multiStartMode = MultiStartMode::Always;

    Solver solver(*ir, an);
    auto res = solver.solve(opts, true);
    REQUIRE(res.success);
}

TEST_CASE("Multi-start rescues refrigeration_compressor without initials",
          "[.][solver][multistart][examples]") {
    fs::path exDir = "../examples/";
    if (const char* e = std::getenv("COOLSOLVE_EXAMPLES_DIR")) exDir = fs::path(e);
    fs::path model = exDir / "refrigeration_compressor.eescode";
    if (!fs::exists(model)) { SUCCEED("model not found, skipping"); return; }

    EESParser parser;
    auto pr = parser.parseFile(model.string());
    REQUIRE(pr.success);
    auto ir = std::make_unique<IR>(IR::fromAST(pr.program));
    try { inferVariables(*ir); initializeVariables(*ir); } catch (...) {}
    auto an = StructuralAnalyzer::analyze(*ir);

    SolverOptions opts;
    opts.tolerance = 1e-6;
    opts.timeoutSeconds = 30;
    opts.solverPipeline = {
        SolverStrategy::Newton, SolverStrategy::TrustRegion,
        SolverStrategy::LevenbergMarquardt, SolverStrategy::BisectionND,
        SolverStrategy::Homotopy, SolverStrategy::Partitioned,
    };
    opts.multiStartMode = MultiStartMode::Always;

    Solver solver(*ir, an);
    auto res = solver.solve(opts, true);
    REQUIRE(res.success);
}

// ----------------------------------------------------------------------------
// "Try Harder" / Deep Search integration.
//
// Verifies the two essential properties of the feature:
//   1. With default options (Newton-only pipeline, multiStartMode = InDeepSearch),
//      the singular-scale block FAILS in a normal solve because multi-start
//      does not engage (deepSearch = false).
//   2. Setting `deepSearch = true` engages the full deep-search pipeline AND
//      flips multi-start on (because mode == InDeepSearch); the block now
//      converges.
// ----------------------------------------------------------------------------
TEST_CASE("Deep search (Try Harder) rescues the singular block",
          "[solver][multistart][deepsearch]") {
    // 1. Baseline: defaults fail.
    {
        SolverOptions opts;
        opts.tolerance = 1e-8;
        opts.maxIterations = 50;
        // Default: solverPipeline = {Newton}, multiStartMode = InDeepSearch,
        // deepSearch = false → multi-start NOT engaged.
        auto out = solveFromSource(kSingularScaleBlock, opts);
        REQUIRE_FALSE(out.success);
    }

    // 2. Try Harder: same options but with deepSearch = true.
    //    The deep-search pipeline (default = all solvers, sequential) now
    //    runs, tearing + symbolic reduction are forced on, AND multi-start
    //    engages (mode == InDeepSearch, deepSearch == true).  The block now
    //    converges via one of these mechanisms.
    {
        SolverOptions opts;
        opts.tolerance = 1e-8;
        opts.maxIterations = 50;
        opts.deepSearch = true;
        auto out = solveFromSource(kSingularScaleBlock, opts);
        REQUIRE(out.success);
        REQUIRE_THAT(out.result.variables.at("x"), WithinAbs(0.1, 1e-4));
        REQUIRE_THAT(out.result.variables.at("y"), WithinAbs(0.1, 1e-4));
    }
}

TEST_CASE("isMultiStartActive follows the mode + deepSearch combination",
          "[solver][multistart][deepsearch]") {
    SolverOptions opts;
    // Default: InDeepSearch, deepSearch off → not active.
    REQUIRE_FALSE(opts.isMultiStartActive());

    // Flip deepSearch on → active.
    opts.deepSearch = true;
    REQUIRE(opts.isMultiStartActive());

    // Always → active regardless of deepSearch.
    opts.multiStartMode = MultiStartMode::Always;
    opts.deepSearch = false;
    REQUIRE(opts.isMultiStartActive());

    // Never → never active.
    opts.multiStartMode = MultiStartMode::Never;
    opts.deepSearch = true;
    REQUIRE_FALSE(opts.isMultiStartActive());

    // Legacy bool override: false wins regardless of mode.
    opts.multiStartMode = MultiStartMode::Always;
    opts.deepSearch = true;
    opts.multiStartEnabled = false;
    REQUIRE_FALSE(opts.isMultiStartActive());
}
