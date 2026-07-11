#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/variable_inference.h"
#include "coolsolve/integral/integral_solver.h"

#include <catch2/catch_test_macros.hpp>
#include <atomic>
#include <chrono>
#include <cmath>
#include <thread>

using namespace coolsolve;

namespace {
// Parse + build IR + analyse + construct an IntegralSolver on the model.
// `opts` defaults to a fresh SolverOptions; tests that exercise the
// progressCallback / cancelToken wiring pass a custom one.
IntegralSolver makeSolver(const std::string& src, std::unique_ptr<IR>& irOut,
                          SolverOptions opts = SolverOptions{}) {
    EESParser parser;
    auto pr = parser.parse(src);
    REQUIRE(pr.success);
    irOut = std::make_unique<IR>(IR::fromAST(pr.program));
    inferVariables(*irOut);
    initializeVariables(*irOut);
    auto analysis = StructuralAnalyzer::analyze(*irOut);
    // The full IR is intentionally non-square (integration var t is free), so
    // `analysis.success` may be false. The IntegralSolver builds its own
    // reduced (square) algebraic analysis internally; the full analysis is only
    // consulted for the problem extraction and tolerates failure.
    static thread_local std::vector<std::pair<std::string, Program>> keepalive;
    keepalive.push_back({src, std::move(pr.program)});
    return IntegralSolver(keepalive.back().second, *irOut, analysis, opts);
}
}  // namespace

// ============================================================================
// Phase 5 — end-to-end tests for the IntegralSolver time-march loop.
// See docs/integral_table_plan.md §Phase 5.
// ============================================================================

TEST_CASE("Integral e2e: exponential decay", "[integral][e2e]") {
    std::unique_ptr<IR> ir;
    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 4)\n"
        "dydt = -y\n", ir);

    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    opt.maxSteps = 2000;
    IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);
    REQUIRE(r.totalSteps > 0);

    const double exact = std::exp(-4.0);
    double yEnd = r.table.interpolate("y", 4.0);
    INFO("y(4)=" << yEnd << " exact=" << exact);
    CHECK(std::abs(yEnd - exact) < 1e-4);
}

TEST_CASE("Integral e2e: harmonic oscillator (two coupled states)", "[integral][e2e]") {
    std::unique_ptr<IR> ir;
    auto solver = makeSolver(
        "y = 0 + integral(dydt, t, 0, 1)\n"
        "z = 1 + integral(dzdt, t, 0, 1)\n"
        "dydt = z\n"
        "dzdt = -y\n", ir);

    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    opt.maxSteps = 2000;
    IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);

    double yEnd = r.table.interpolate("y", 1.0);
    double zEnd = r.table.interpolate("z", 1.0);
    INFO("y(1)=" << yEnd << " exact=" << std::sin(1.0));
    CHECK(std::abs(yEnd - std::sin(1.0)) < 1e-4);
    CHECK(std::abs(zEnd - std::cos(1.0)) < 1e-4);
}

TEST_CASE("Integral e2e: RK45 adaptive on decay", "[integral][e2e]") {
    std::unique_ptr<IR> ir;
    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 4)\n"
        "dydt = -y\n", ir);

    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK45;
    opt.relTol = 1e-7;
    opt.absTol = 1e-9;
    IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);
    CHECK(r.rejectedSteps >= 0);

    const double exact = std::exp(-4.0);
    double yEnd = r.table.interpolate("y", 4.0);
    CHECK(std::abs(yEnd - exact) < 1e-5);
}

TEST_CASE("Integral e2e: algebraic variable coupled to a state", "[integral][e2e]") {
    // y' = -y, with an auxiliary algebraic variable w = 2*y that the table records.
    std::unique_ptr<IR> ir;
    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 2)\n"
        "dydt = -y\n"
        "w = 2*y\n"
        "$IntegralTable t y w\n", ir);

    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    opt.maxSteps = 2000;
    IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);

    double yEnd = r.table.interpolate("y", 2.0);
    double wEnd = r.table.interpolate("w", 2.0);
    double yExact = std::exp(-2.0);
    CHECK(std::abs(yEnd - yExact) < 1e-4);
    CHECK(std::abs(wEnd - 2.0 * yExact) < 2e-4);
}

TEST_CASE("Integral e2e: fixed step from INTEGRAL 5th argument", "[integral][e2e]") {
    std::unique_ptr<IR> ir;
    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 1, 0.1)\n"
        "dydt = -y\n", ir);

    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);
    // 10 steps of 0.1 over [0,1].
    CHECK(r.totalSteps == 10);
    double yEnd = r.table.interpolate("y", 1.0);
    CHECK(std::abs(yEnd - std::exp(-1.0)) < 1e-4);
}

TEST_CASE("Integral e2e: $IntegralTable directive columns honoured", "[integral][e2e]") {
    std::unique_ptr<IR> ir;
    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 1)\n"
        "dydt = -y\n"
        "w = 3*y\n"
        "$IntegralTable t:0.5 y w\n", ir);

    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    opt.maxSteps = 1000;
    IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);
    const auto& cols = r.table.columns();
    REQUIRE(cols.size() == 3);
    CHECK(cols[0] == "t");
    CHECK(cols[1] == "y");
    CHECK(cols[2] == "w");
    // Output interval 0.5 => at least t=0, 0.5, 1.0 rows.
    CHECK(r.table.numRows() >= 3);
}

// ============================================================================
// Regression tests: IntegralSolver ↔ GUI/server integration contract.
//
// These guard against a class of bugs that the numerical e2e tests above cannot
// catch because they construct the solver with default SolverOptions (no
// progressCallback, no cancelToken). In server mode the runner propagates both
// into the IntegralSolver; the time-march loop must (a) NOT relay the per-block
// progressCallback to its thousands of internal algebraic solves, and (b) still
// honour the cancelToken so the Stop button works. See docs/integral_table_plan.md
// §"2026-07-11 — Bugfix: GUI freeze on integral solve".
// ============================================================================

TEST_CASE("Integral e2e: internal solves do not relay the progress callback",
          "[integral][e2e][progress]") {
    // Regression for the GUI freeze: the time-march loop reuses the algebraic
    // Solver thousands of times. The per-block progressCallback is meant for a
    // single top-level solve; if it fires inside every step it floods the GUI's
    // SSE stream (24k+ events for a 1000-step RK4 solve) and freezes the
    // browser. The IntegralSolver must suppress it for internal solves.
    std::unique_ptr<IR> ir;
    SolverOptions opts;
    int callbackCalls = 0;
    opts.progressCallback = [&callbackCalls](int, int, const std::string&,
                                             int, double) {
        ++callbackCalls;
    };

    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 4)\n"
        "dydt = -y\n", ir, opts);

    IntegratorOptions iopt;
    iopt.method = IntegratorOptions::RK4;
    iopt.maxSteps = 1000;  // ~24 000 per-block events if the callback leaked
    IntegralSolveResult r = solver.solve(iopt);
    REQUIRE(r.success);

    // Not a single per-block event should leak from the internal solves.
    // (Pre-fix this asserted value was ~24 000.)
    INFO("progressCallback fired " << callbackCalls << " times (expected 0)");
    CHECK(callbackCalls == 0);
    // Sanity: suppressing the callback must not break the numerics.
    CHECK(std::abs(r.table.interpolate("y", 4.0) - std::exp(-4.0)) < 1e-4);
}

TEST_CASE("Integral e2e: pre-set cancel token aborts before integration",
          "[integral][e2e][cancel]") {
    // The IntegralSolver must carry SolverOptions::cancelToken through to its
    // algebraic solves (the GUI Stop button sets this flag). With the token
    // already raised, the solve aborts without completing the march.
    std::unique_ptr<IR> ir;
    SolverOptions opts;
    std::atomic<bool> cancelled{true};
    opts.cancelToken = &cancelled;

    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 4)\n"
        "dydt = -y\n", ir, opts);

    IntegratorOptions iopt;
    iopt.method = IntegratorOptions::RK4;
    iopt.maxSteps = 1000;
    IntegralSolveResult r = solver.solve(iopt);

    REQUIRE(!r.success);
    CHECK(r.totalSteps == 0);  // never reached the march loop
}

TEST_CASE("Integral e2e: march loop stops promptly when cancelled mid-flight",
          "[integral][e2e][cancel]") {
    // Guards the per-step cancel poll added to the march loop. A huge step
    // budget makes the full march take well over a second; raising the token
    // partway must interrupt it long before it completes.
    std::unique_ptr<IR> ir;
    SolverOptions opts;
    std::atomic<bool> cancelled{false};
    opts.cancelToken = &cancelled;

    auto solver = makeSolver(
        "y = 1 + integral(dydt, t, 0, 4)\n"
        "dydt = -y\n", ir, opts);

    IntegratorOptions iopt;
    iopt.method = IntegratorOptions::RK4;
    iopt.maxSteps = 1000000;  // h = 4e-6 → ~10⁶ steps, ~1 s of work
    IntegralSolveResult r;
    std::thread worker([&] { r = solver.solve(iopt); });

    // Let the march get underway, then request cancellation.
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    cancelled.store(true);
    worker.join();

    REQUIRE(!r.success);
    // The march did not run to completion.
    CHECK(r.totalSteps < iopt.maxSteps);
}
