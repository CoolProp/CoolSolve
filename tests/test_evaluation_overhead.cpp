/**
 * @file test_evaluation_overhead.cpp
 * @brief The evaluation without copies of the state and the time march of the blocks that change
 *        (docs/solver_roadmap.md §6.1, steps 1 to 4).
 *
 * - BlockEvaluator::evaluate reads the variables a block does not own from the state, by reference;
 *   an explicit assignment evaluated for its value carries no gradient.
 * - The routine tables are referenced, not copied, by the expression evaluators.
 * - Newton1D on a block of one unknown probes around the other variables of the state only.
 * - Solver::solveBlocks solves part of the system; the INTEGRAL march solves only the blocks that
 *   depend on the time or on a state.
 * - callsPropertyFunctions decides the CoolProp warm-up.
 * - EulerImplicit fails on a step whose iteration does not converge; integral.md names the method used;
 *   a fixed step too large for an explicit method says so.
 */
#include "coolsolve/evaluator.h"
#include "coolsolve/integral/integral_solver.h"
#include "coolsolve/integral/integrator.h"
#include "coolsolve/ir.h"
#include "coolsolve/parser.h"
#include "coolsolve/solver.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/variable_inference.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <memory>
#include <string>

using namespace coolsolve;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

/// A model parsed, with its IR (inferred and initialised) and its structural analysis.
struct Model {
    Program program;
    std::unique_ptr<IR> ir;
    StructuralAnalysisResult analysis;
};

std::unique_ptr<Model> buildModel(const std::string& src) {
    auto m = std::make_unique<Model>();
    EESParser parser;
    auto pr = parser.parse(src);
    REQUIRE(pr.success);
    m->program = std::move(pr.program);
    m->ir = std::make_unique<IR>(IR::fromAST(m->program));
    inferVariables(*m->ir);
    initializeVariables(*m->ir);
    m->analysis = StructuralAnalyzer::analyze(*m->ir);
    return m;
}

/// The index of the block that determines `name`.
size_t blockOf(const StructuralAnalysisResult& analysis, const std::string& name) {
    for (size_t b = 0; b < analysis.blocks.size(); ++b)
        for (const auto& v : analysis.blocks[b].variables)
            if (eesCompareFolded(v, name) == 0) return b;
    FAIL("no block determines " << name);
    return 0;
}

}  // namespace

// ============================================================================
// BlockEvaluator: the external values are read by reference
// ============================================================================

TEST_CASE("BlockEvaluator reads the other variables from the state it is given, its own from x",
          "[eval-overhead]") {
    auto m = buildModel("a = 3\nb = 4\ny = 2*A + b\n");
    REQUIRE(m->analysis.success);
    BlockEvaluator be(m->analysis.blocks[blockOf(m->analysis, "y")], *m->ir);
    REQUIRE(be.getVariables().size() == 1);

    // The state holds the block variable too (a stale value): x wins. The lookup ignores the case.
    VariableValues state{{"A", 3.0}, {"b", 4.0}, {"y", 1000.0}};
    const auto withJac = be.evaluate({5.0}, state);
    REQUIRE(withJac.residuals.size() == 1);
    CHECK(withJac.residuals[0] == 5.0 - 10.0);
    REQUIRE(withJac.jacobian.size() == 1);
    CHECK(withJac.jacobian[0][0] == 1.0);
    CHECK(withJac.warnings.empty());

    // The state is read, not copied: a later change is seen by the next evaluation
    state["b"] = 40.0;
    CHECK(be.evaluate({5.0}, state).residuals[0] == 5.0 - 46.0);

    // A variable of the equation that the state does not hold keeps the former default, 1, with a warning
    const VariableValues partial{{"a", 3.0}};
    const auto missing = be.evaluate({0.0}, partial);
    CHECK(missing.residuals[0] == -(6.0 + 1.0));
    REQUIRE(missing.warnings.size() == 1);
    CHECK_THAT(missing.warnings[0], ContainsSubstring("Undefined variable: b"));
}

TEST_CASE("An explicit assignment is evaluated without gradient, to the same value", "[eval-overhead]") {
    auto m = buildModel(
        "a = 0.3\n"
        "y = sin(a)^2 + exp(a)/3 + max(a, 0.1)\n"     // explicit
        "z^2 = a\n"                                    // not an assignment of z
        "w = 0.5*w + a\n");                            // w on both sides
    REQUIRE(m->analysis.success);
    BlockEvaluator y(m->analysis.blocks[blockOf(m->analysis, "y")], *m->ir);
    BlockEvaluator z(m->analysis.blocks[blockOf(m->analysis, "z")], *m->ir);
    BlockEvaluator w(m->analysis.blocks[blockOf(m->analysis, "w")], *m->ir);
    CHECK(y.isExplicitAssignment());
    CHECK_FALSE(z.isExplicitAssignment());
    CHECK_FALSE(w.isExplicitAssignment());

    const VariableValues state{{"a", 0.3}};
    const auto full = y.evaluate({0.0}, state, {}, true);
    const auto valueOnly = y.evaluate({0.0}, state, {}, false);
    CHECK(valueOnly.jacobian.empty());
    CHECK(valueOnly.scale.empty());
    CHECK(valueOnly.residuals[0] == full.residuals[0]);   // bit for bit

    // The other blocks keep their gradients in a residual-only evaluation (the clamping of property inputs
    // looks at them); the Jacobian is still not computed
    CHECK(w.evaluate({1.0}, state, {}, false).jacobian.empty());
    CHECK(w.evaluate({1.0}, state, {}, true).jacobian[0][0] == 0.5);
}

TEST_CASE("Solver: an explicit block is solved from the state, other blocks see its value", "[eval-overhead]") {
    auto m = buildModel(
        "a = 2\n"
        "b = a^2 + 1\n"
        "c^3 + c = b\n");     // 1D implicit after the explicit ones
    REQUIRE(m->analysis.success);
    Solver solver(*m->ir, m->analysis);
    const SolveResult r = solver.solve();
    REQUIRE(r.success);
    CHECK(r.variables.at("b") == 5.0);
    CHECK_THAT(r.variables.at("c"), WithinRel(1.5159802, 1e-6));
    CHECK(r.blockTraces[blockOf(m->analysis, "b")].solverType == "Explicit");
}

// ============================================================================
// Routine tables and external values of ExpressionEvaluator
// ============================================================================

TEST_CASE("ExpressionEvaluator: shared routine tables, external values", "[eval-overhead]") {
    auto m = buildModel(
        "FUNCTION f(x)\n  f := 2*x\nEND\n"
        "FUNCTION g(x)\n  g := x + 1\nEND\n"
        "q = 4\n"
        "y = f(3) + g(1) + q\n");
    REQUIRE(m->ir->getFunctions().size() == 2);
    const EquationInfo* eq = nullptr;
    for (const auto& e : m->ir->getEquations())
        if (e.lhs && e.lhs->is<Variable>() && e.lhs->as<Variable>().name == "y") eq = &e;
    REQUIRE(eq);

    FunctionTable onlyF;
    for (const auto& func : m->ir->getFunctions())
        if (func.name == "f") onlyF[func.name] = func;
    const ProcedureTable none;
    const VariableValues state{{"q", 4.0}};

    ExpressionEvaluator ev(0);
    ev.setRoutines(&onlyF, &none);
    ev.setExternalValues(&state, nullptr);
    // g is not in the shared table: registering it copies the table first and leaves it unchanged
    for (const auto& func : m->ir->getFunctions())
        if (func.name == "g") ev.registerFunction(func);
    CHECK(onlyF.size() == 1);

    // f(3) = 6, g(1) = 2, q = 4 from the external values
    const ADValue value = ev.evaluate(eq->rhs);
    CHECK(value.value == 6.0 + 2.0 + 4.0);

    // A value set in the evaluator hides the external one
    ev.setVariable("q", ADValue::constant(10.0, 0));
    CHECK(ev.evaluate(eq->rhs).value == 6.0 + 2.0 + 10.0);
}

TEST_CASE("ExpressionEvaluator: external string values are read when not set locally", "[eval-overhead]") {
    const StringValues strings{{"F$", "Water"}};
    ExpressionEvaluator ev(0);
    ev.setExternalValues(nullptr, &strings);
    CHECK(ev.hasStringVariable("f$"));
    CHECK(ev.getStringVariable("f$") == "Water");
    ev.setStringVariable("f$", "R134a");
    CHECK(ev.getStringVariable("F$") == "R134a");
    CHECK_FALSE(ev.hasStringVariable("g$"));
    CHECK_THROWS(ev.getStringVariable("g$"));
}

// ============================================================================
// Newton1D: the probes of a block of one unknown
// ============================================================================

TEST_CASE("Newton1D on a block reads its probes from the state, without its own variable", "[eval-overhead]") {
    // F = 1 except in a narrow dip around 777.7 (two roots at 777.7 +- 0.42): only a probe at the value of
    // a variable near the dip finds the sign change (Phase 2)
    const double c = 777.7;
    Newton1DSolver::Eval1D eval = [&](double x) {
        const double u = (x - c) / 0.5;
        const double f = 1.0 - 2.0 * std::exp(-u * u);
        const double df = 2.0 * std::exp(-u * u) * 2.0 * u / 0.5;
        return Newton1DSolver::Eval1DResult(f, df);
    };
    SolverOptions options;
    options.tolerance = 1e-9;

    // The dip is at the value of another variable of the state: found, as with the map of the other overload
    double x = 5000.0;
    const VariableValues state{{"c", c}, {"x", 5000.0}};
    REQUIRE(Newton1DSolver::solve(eval, x, state, "x", options) == SolverStatus::Success);
    double xMap = 5000.0;
    REQUIRE(Newton1DSolver::solve(eval, xMap, std::map<std::string, double>{{"c", c}}, options) == SolverStatus::Success);
    CHECK(x == xMap);
    CHECK_THAT(std::abs(x - c), WithinAbs(0.5 * std::sqrt(std::log(2.0)), 1e-6));

    // The value of the block's own variable (its current guess) is not a probe: the dip is not found
    double own = 5000.0;
    const VariableValues onlyOwn{{"X", c}};
    CHECK(Newton1DSolver::solve(eval, own, onlyOwn, "x", options) != SolverStatus::Success);
}

// ============================================================================
// Solver::solveBlocks
// ============================================================================

TEST_CASE("Solver::solveBlocks solves the given blocks only, with the message of solve()", "[eval-overhead]") {
    auto m = buildModel("a = 2\nb = 3*a\nc = b + 1\n");
    REQUIRE(m->analysis.success);
    Solver solver(*m->ir, m->analysis);
    REQUIRE(solver.solve().success);
    const auto& state = solver.getEvaluator().getAllVariables();
    CHECK(state.at("c") == 7.0);

    // b is not solved again: c follows the value the state holds
    solver.setGuess("b", 100.0);
    std::string error;
    REQUIRE(solver.solveBlocks({blockOf(m->analysis, "c")}, SolverOptions{}, &error) == SolverStatus::Success);
    CHECK(error.empty());
    CHECK(state.at("b") == 100.0);
    CHECK(state.at("c") == 101.0);
    REQUIRE(solver.solveBlocks({blockOf(m->analysis, "b"), blockOf(m->analysis, "c")}, SolverOptions{}) ==
            SolverStatus::Success);
    CHECK(state.at("c") == 7.0);

    // A block that fails: the status and the message are those of solve()
    auto bad = buildModel("x^2 + 1 = 0\n");
    REQUIRE(bad->analysis.success);
    Solver failing(*bad->ir, bad->analysis);
    const SolveResult r = failing.solve();
    REQUIRE_FALSE(r.success);
    Solver failingAgain(*bad->ir, bad->analysis);
    std::string message;
    CHECK(failingAgain.solveBlocks({0}, SolverOptions{}, &message) == r.status);
    CHECK(message == r.errorMessage);
    CHECK_THAT(message, ContainsSubstring("Block 0 (size 1, vars: x) failed"));
}

// ============================================================================
// The march solves the blocks that depend on t or on a state
// ============================================================================

TEST_CASE("Integral march: the blocks independent of t and of the states keep their value", "[eval-overhead]") {
    // k (an implicit block, k = 1) and g are solved once at the initial point; z depends on the state
    auto m = buildModel(
        "y = 1 + integral(dydt, t, 0, 2)\n"
        "k^3 + k = 2\n"
        "g = 2*k\n"
        "z^2 = y + g\n"
        "dydt = -k*y + 0*z\n");
    IntegralSolver solver(m->program, *m->ir, m->analysis, SolverOptions{});
    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    opt.maxSteps = 200;
    const IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);
    CHECK(r.method == "RK4");
    CHECK_THAT(r.table.interpolate("y", 2.0), WithinRel(std::exp(-2.0), 1e-8));
    // The final point holds the static variables and those of the last solve
    const auto& vars = r.algebraicResult.variables;
    CHECK_THAT(vars.at("k"), WithinAbs(1.0, 1e-12));
    CHECK_THAT(vars.at("g"), WithinAbs(2.0, 1e-12));
    CHECK_THAT(vars.at("z"), WithinRel(std::sqrt(std::exp(-2.0) + 2.0), 1e-8));
    CHECK(r.algebraicResult.success);
}

TEST_CASE("Integral march: a check equation that reads t or a state is verified at every evaluation", "[eval-overhead]") {
    // z is defined twice: the second equation is a check equation, true at t = 0 only
    auto m = buildModel(
        "y = 1 + integral(dydt, t, 0, 1)\n"
        "dydt = -y\n"
        "z = 2*y\n"
        "z = 2*y*(1 + t)\n"
        "w = u\n");
    IntegralSolver solver(m->program, *m->ir, m->analysis, SolverOptions{});
    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    opt.maxSteps = 20;
    const IntegralSolveResult r = solver.solve(opt);
    REQUIRE_FALSE(r.success);
    CHECK_THAT(r.errorMessage, ContainsSubstring("Check equation failed"));

    // The same check equation, true at every t, does not stop the march
    auto holds = buildModel(
        "y = 1 + integral(dydt, t, 0, 1)\n"
        "dydt = -y\n"
        "z = 2*y\n"
        "z = 2*y*(1 + 0*t)\n"
        "w = u\n");
    IntegralSolver fine(holds->program, *holds->ir, holds->analysis, SolverOptions{});
    CHECK(fine.solve(opt).success);
}

TEST_CASE("IntegralSolveResult names the method used, Richardson included", "[eval-overhead]") {
    auto m = buildModel("y = 1 + integral(dydt, t, 0, 1)\ndydt = -y\n");
    IntegralSolver solver(m->program, *m->ir, m->analysis, SolverOptions{});
    IntegratorOptions opt;
    opt.maxSteps = 20;
    opt.method = IntegratorOptions::EulerImplicit;
    CHECK(solver.solve(opt).method == "EulerImplicit");
    opt.method = IntegratorOptions::RK4;
    opt.richardson = true;
    CHECK(solver.solve(opt).method == "RK4 + Richardson");
    opt.method = IntegratorOptions::RK45;   // Richardson applies to the fixed-step methods only
    CHECK(solver.solve(opt).method == "RK45");
}

TEST_CASE("A fixed step too large for an explicit method says that it may be too large", "[eval-overhead]") {
    // w^2 = 10 - y has no solution once y > 10: the algebraic solve fails when the states diverge
    const std::string model =
        "y = 1 + integral(dydt, t, 0, 10, %STEP%)\n"
        "dydt = -50*y\n"
        "w^2 = 10 - y\n";
    auto withStep = [&](const std::string& h) {
        std::string src = model;
        src.replace(src.find("%STEP%"), 6, h);
        return src;
    };

    auto unstable = buildModel(withStep("1"));
    IntegralSolver diverging(unstable->program, *unstable->ir, unstable->analysis, SolverOptions{});
    IntegratorOptions opt;
    opt.method = IntegratorOptions::RK4;
    const IntegralSolveResult r = diverging.solve(opt);
    REQUIRE_FALSE(r.success);
    CHECK_THAT(r.errorMessage, ContainsSubstring("Algebraic solve failed at t="));
    CHECK_THAT(r.errorMessage, ContainsSubstring("The step h = 1 may be too large for RK4"));

    // A stable step solves the same model; a failure of RK45 carries no such hint (it controls its step)
    auto stable = buildModel(withStep("0.01"));
    IntegralSolver fine(stable->program, *stable->ir, stable->analysis, SolverOptions{});
    CHECK(fine.solve(opt).success);
}

// ============================================================================
// EulerImplicit: a step whose iteration does not converge is an error
// ============================================================================

TEST_CASE("EulerImplicit fails on a step that does not converge, with the time and the residual", "[eval-overhead]") {
    // f = -10 sign(y): the Newton iteration of a step of 1 jumps between -9 and 11 and never converges
    RHSFunction f = [](double, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1);
        r(0) = y(0) > 0.0 ? -10.0 : 10.0;
        return r;
    };
    auto ig = createIntegrator(IntegratorOptions::EulerImplicit);
    Eigen::VectorXd y0(1);
    y0(0) = 1.0;
    IntegratorOptions opt;
    try {
        ig->step(2.0, y0, f, 1.0, opt);
        FAIL("the step was accepted");
    } catch (const std::runtime_error& e) {
        const std::string what = e.what();
        CHECK_THAT(what, ContainsSubstring("EulerImplicit"));
        CHECK_THAT(what, ContainsSubstring("from t = 2 to 3"));
        CHECK_THAT(what, ContainsSubstring("did not converge after 30 iterations"));
        CHECK_THAT(what, ContainsSubstring("largest residual"));
    }
}

TEST_CASE("EulerImplicit accepts a step converged to the rounding of a large state", "[eval-overhead]") {
    // y ~ 1E8: a correction cannot go below ulp(y) ~ 1.5E-8, far above the absolute test 1E-12
    RHSFunction f = [](double, const Eigen::VectorXd& y) { return Eigen::VectorXd(-y); };
    auto ig = createIntegrator(IntegratorOptions::EulerImplicit);
    Eigen::VectorXd y0(1);
    y0(0) = 1.0e8;
    IntegratorOptions opt;
    const StepResult sr = ig->step(0.0, y0, f, 0.1, opt);
    CHECK(sr.accepted);
    CHECK_THAT(sr.yNew(0), WithinRel(1.0e8 / 1.1, 1e-14));
}

// ============================================================================
// CoolProp warm-up only for a model that calls CoolProp
// ============================================================================

TEST_CASE("callsPropertyFunctions: property functions, routines that call them, CALL psychprops", "[eval-overhead]") {
    CHECK_FALSE(callsPropertyFunctions(*buildModel("x = 2\ny = exp(x) + max(x, 1)\n")->ir));
    CHECK(callsPropertyFunctions(*buildModel("h = enthalpy(Water, T=20, P=101.325)\n")->ir));
    CHECK(callsPropertyFunctions(*buildModel(
        "FUNCTION hw(T)\n  hw := enthalpy(Water, T=T, P=101.325)\nEND\n"
        "FUNCTION twice(T)\n  twice := 2*hw(T)\nEND\n"
        "y = twice(20)\n")->ir));
    CHECK_FALSE(callsPropertyFunctions(*buildModel(
        "FUNCTION unused(T)\n  unused := enthalpy(Water, T=T, P=101.325)\nEND\n"
        "y = 3\n")->ir));   // a routine that no equation calls
    CHECK(callsPropertyFunctions(*buildModel(
        "CALL psychprops(25, 101.325, 0.5 : Tout, v, h, s, u, W, R, Twb, Tdp)\n")->ir));
}
