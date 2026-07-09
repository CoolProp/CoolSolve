#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/variable_inference.h"
#include "coolsolve/integral/integral_solver.h"

#include <catch2/catch_test_macros.hpp>
#include <cmath>

using namespace coolsolve;

namespace {
// Parse + build IR + analyse + construct an IntegralSolver on the model.
IntegralSolver makeSolver(const std::string& src, std::unique_ptr<IR>& irOut) {
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
    return IntegralSolver(keepalive.back().second, *irOut, analysis, SolverOptions{});
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
