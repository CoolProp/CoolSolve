/**
 * @file test_integrator_bdf.cpp
 * @brief The variable-order BDF integrator and the Newton iteration it shares with EulerImplicit
 *        (docs/solver_roadmap.md §6.1, step 6).
 *
 * - BDF against analytical solutions, its error following the tolerance, a stiff problem in few steps, its
 *   interpolant, a new history at another point, the rejection of a step whose iteration does not converge.
 * - The march of IntegralSolver with BDF: rows at an output interval solved at the interpolated states, the count
 *   of the evaluations of the right-hand side.
 * - EulerImplicit: a weighted test, a Jacobian kept from step to step, the full Newton iteration across a kink.
 */
#include "integrator_test_util.h"
#include "coolsolve/integral/integral_solver.h"
#include "coolsolve/integral/integrator.h"
#include "coolsolve/ir.h"
#include "coolsolve/parser.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/variable_inference.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <stdexcept>

using namespace coolsolve;
using namespace coolsolve::test;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

/// dy/dt = -y, counting the evaluations
RHSFunction decay(int* evaluations = nullptr) {
    return [evaluations](double, const Eigen::VectorXd& y) {
        if (evaluations) ++*evaluations;
        return Eigen::VectorXd(-y);
    };
}

Eigen::VectorXd one(double value) {
    Eigen::VectorXd y(1);
    y(0) = value;
    return y;
}

/// A model parsed, with its IR and analysis, kept alive with the solver built on it.
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

}  // namespace

// ============================================================================
// BDF on analytical solutions
// ============================================================================

TEST_CASE("BDF: names, factory and parsing", "[integrator][bdf]") {
    auto ig = createIntegrator(IntegratorOptions::BDF);
    CHECK(std::string(ig->name()) == "BDF");
    CHECK(methodToString(IntegratorOptions::BDF) == "BDF");
    IntegratorOptions::Method m = IntegratorOptions::RK4;
    CHECK(parseIntegralMethod("bdf", m));
    CHECK(m == IntegratorOptions::BDF);
    m = IntegratorOptions::RK4;
    CHECK(parseIntegralMethod("NDF", m));
    CHECK(m == IntegratorOptions::BDF);
    CHECK_FALSE(parseIntegralMethod("bdf2", m));
}

TEST_CASE("BDF: dy/dt = -y, the error follows the tolerance", "[integrator][bdf]") {
    const double exact = std::exp(-4.0);
    double previousError = 1.0;
    int previousEvaluations = 0;
    for (double relTol : {1e-3, 1e-5, 1e-7, 1e-9}) {
        auto ig = createIntegrator(IntegratorOptions::BDF);
        IntegratorOptions opt;
        opt.relTol = relTol;
        opt.absTol = relTol * 1e-3;
        int evaluations = 0, taken = 0, rejected = 0;
        const double yEnd = marchAdaptive(*ig, 0.0, 4.0, one(1.0), decay(&evaluations), 0.5, opt, &taken, &rejected)(0);
        const double error = std::abs(yEnd - exact) / exact;
        INFO("relTol " << relTol << ": error " << error << ", " << taken << " steps, " << rejected << " rejected, "
                       << evaluations << " evaluations");
        CHECK(error < 50.0 * relTol);
        CHECK(error < previousError);
        CHECK(evaluations > previousEvaluations);
        previousError = error;
        previousEvaluations = evaluations;
    }
}

TEST_CASE("BDF: coupled harmonic oscillator, higher orders", "[integrator][bdf]") {
    RHSFunction f = [](double, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(2);
        r(0) = y(1);
        r(1) = -y(0);
        return r;
    };
    auto ig = createIntegrator(IntegratorOptions::BDF);
    IntegratorOptions opt;
    opt.relTol = 1e-8;
    opt.absTol = 1e-11;
    Eigen::VectorXd y0(2);
    y0 << 0.0, 1.0;
    int taken = 0;
    const Eigen::VectorXd yEnd = marchAdaptive(*ig, 0.0, 2.0, y0, f, 0.1, opt, &taken);
    CHECK_THAT(yEnd(0), WithinAbs(std::sin(2.0), 1e-6));
    CHECK_THAT(yEnd(1), WithinAbs(std::cos(2.0), 1e-6));
    // An order-1 method needs about 1/sqrt(1e-8) = 1e4 steps for this accuracy: the order rises
    INFO(taken << " steps");
    CHECK(taken < 1000);
}

TEST_CASE("BDF: a stiff problem in few steps, where RK45 is held by stability", "[integrator][bdf]") {
    // y' = -1000 (y - cos t): the solution follows cos t after a transient of 1/1000; RK45 is limited to h ~ 3e-3
    RHSFunction f = [](double t, const Eigen::VectorXd& y) { return Eigen::VectorXd(-1000.0 * (y.array() - std::cos(t))); };
    IntegratorOptions opt;
    opt.relTol = 1e-6;
    opt.absTol = 1e-9;
    opt.maxSteps = 100000;
    int bdfSteps = 0, rkSteps = 0;
    auto bdf = createIntegrator(IntegratorOptions::BDF);
    const double yBdf = marchAdaptive(*bdf, 0.0, 5.0, one(0.0), f, 0.1, opt, &bdfSteps)(0);
    auto rk = createIntegrator(IntegratorOptions::RK45);
    const double yRk = marchAdaptive(*rk, 0.0, 5.0, one(0.0), f, 0.1, opt, &rkSteps)(0);
    // y(t) = cos t + 1e-3 sin t - (1 + ...) e^{-1000 t}, to 1e-6: the steady part
    const double exact = (1e6 * std::cos(5.0) + 1e3 * std::sin(5.0)) / (1e6 + 1.0);
    INFO("BDF " << bdfSteps << " steps, RK45 " << rkSteps << " steps");
    CHECK_THAT(yBdf, WithinAbs(exact, 1e-5));
    CHECK_THAT(yRk, WithinAbs(exact, 1e-5));
    CHECK(bdfSteps * 5 < rkSteps);
}

TEST_CASE("BDF: the interpolant of the last step", "[integrator][bdf]") {
    auto ig = createIntegrator(IntegratorOptions::BDF);
    IntegratorOptions opt;
    opt.relTol = 1e-8;
    opt.absTol = 1e-11;
    Eigen::VectorXd yq;
    CHECK_FALSE(ig->denseOutput(0.0, yq));   // no step yet

    const RHSFunction f = decay();
    double t = 0.0, h = 0.1;
    Eigen::VectorXd y = one(1.0);
    int checked = 0;
    while (t < 3.0) {
        const StepResult sr = ig->step(t, y, f, std::min(h, 3.0 - t), opt);
        h = sr.nextStep;
        if (!sr.accepted) continue;
        const double tPrev = t;
        t += sr.stepTaken;
        y = sr.yNew;
        REQUIRE(ig->denseOutput(t, yq));
        CHECK_THAT(yq(0), WithinRel(y(0), 1e-12));   // the interpolant passes through the new point
        const double tMid = 0.5 * (tPrev + t);
        REQUIRE(ig->denseOutput(tMid, yq));
        CHECK_THAT(yq(0), WithinRel(std::exp(-tMid), 1e-6));
        ++checked;
    }
    CHECK(checked > 5);
}

TEST_CASE("BDF: a point other than the end of the last step starts a new history", "[integrator][bdf]") {
    IntegratorOptions opt;
    auto used = createIntegrator(IntegratorOptions::BDF);
    StepResult sr;
    do sr = used->step(0.0, one(1.0), decay(), 0.1, opt);
    while (!sr.accepted);
    // Another point: the same step as a fresh integrator, bit for bit
    auto fresh = createIntegrator(IntegratorOptions::BDF);
    const StepResult a = used->step(5.0, one(2.0), decay(), 0.1, opt);
    const StepResult b = fresh->step(5.0, one(2.0), decay(), 0.1, opt);
    CHECK(a.accepted == b.accepted);
    CHECK(a.stepTaken == b.stepTaken);
    CHECK(a.yNew(0) == b.yNew(0));
    CHECK(a.rhsEvaluations == b.rhsEvaluations);
}

TEST_CASE("BDF: the first step is chosen from f, within the step proposed", "[integrator][bdf]") {
    IntegratorOptions opt;
    opt.relTol = 1e-6;
    opt.absTol = 1e-9;
    auto ig = createIntegrator(IntegratorOptions::BDF);
    const StepResult sr = ig->step(0.0, one(1.0), decay(), 10.0, opt);
    CHECK(sr.accepted);
    CHECK(sr.stepTaken < 0.1);   // not the 10 proposed: an error of 1e-6 needs a short first step at order 1
    CHECK(sr.stepTaken > 0.0);
    auto small = createIntegrator(IntegratorOptions::BDF);
    CHECK(small->step(0.0, one(1.0), decay(), 1e-5, opt).stepTaken == 1e-5);   // never more than proposed
}

TEST_CASE("BDF: a step whose iteration fails is rejected and halved, with the reason", "[integrator][bdf]") {
    // The right-hand side cannot be evaluated above y = 1.5: a step that tries it is rejected
    RHSFunction f = [](double, const Eigen::VectorXd& y) {
        if (y(0) > 1.5) throw std::runtime_error("Algebraic solve failed at y > 1.5");
        return Eigen::VectorXd(Eigen::VectorXd::Constant(1, 1.0));
    };
    IntegratorOptions opt;
    auto ig = createIntegrator(IntegratorOptions::BDF);
    StepResult sr = ig->step(0.0, one(1.0), f, 1e-3, opt);
    REQUIRE(sr.accepted);
    sr = ig->step(sr.stepTaken, sr.yNew, f, 1.0, opt);   // to y = 2, beyond the limit
    CHECK_FALSE(sr.accepted);
    CHECK(sr.nextStep == 0.5);
    CHECK_THAT(sr.rejectReason, ContainsSubstring("Algebraic solve failed at y > 1.5"));
}

TEST_CASE("BDF and EulerImplicit: the Newton test is weighted, a large state converges at once", "[integrator][bdf]") {
    // y ~ 1e6: an absolute test would take every iteration of every step
    IntegratorOptions opt;
    int evaluations = 0, steps = 0, rejected = 0;
    auto bdf = createIntegrator(IntegratorOptions::BDF);
    const double yBdf = marchAdaptive(*bdf, 0.0, 1.0, one(1.0e6), decay(&evaluations), 0.01, opt, &steps, &rejected)(0);
    CHECK_THAT(yBdf, WithinRel(1.0e6 * std::exp(-1.0), 1e-5));
    INFO("BDF: " << evaluations << " evaluations for " << steps << " steps, " << rejected << " rejected");
    CHECK(evaluations <= 3 * (steps + rejected) + 5);

    evaluations = 0;
    auto euler = createIntegrator(IntegratorOptions::EulerImplicit);
    const double yEuler = marchFixed(*euler, 0.0, 1.0, one(1.0e6), decay(&evaluations), 100, opt)(0);
    CHECK_THAT(yEuler, WithinRel(1.0e6 * std::pow(1.0 / 1.01, 100), 1e-9));
    INFO("EulerImplicit: " << evaluations << " evaluations for 100 steps");
    CHECK(evaluations <= 3 * 100 + 1);   // one Jacobian, then a predictor and two iterations per step
}

// ============================================================================
// EulerImplicit on the Newton iteration of BDF
// ============================================================================

TEST_CASE("EulerImplicit: the Jacobian is kept from step to step", "[integrator][euler]") {
    int evaluations = 0;
    auto ig = createIntegrator(IntegratorOptions::EulerImplicit);
    IntegratorOptions opt;
    const Eigen::VectorXd y = marchFixed(*ig, 0.0, 1.0, one(1.0), decay(&evaluations), 100, opt);
    CHECK_THAT(y(0), WithinRel(std::pow(1.0 / 1.01, 100), 1e-9));
    // One Jacobian (1 column) at the first step, then a predictor and two iterations per step
    CHECK(evaluations <= 3 * 100 + 1);
}

TEST_CASE("EulerImplicit: across a kink of f, the full Newton iteration converges", "[integrator][euler]") {
    // f = g(t) - 50 max(y - 1, 0), g = 10 from t = 1: the predictor (y = 0.5, where f = 0) is below the kink, the
    // solution above it. With the slope of the predictor (0), the simplified iteration oscillates; at every iterate,
    // it converges.
    RHSFunction f = [](double t, const Eigen::VectorXd& y) {
        Eigen::VectorXd r(1);
        r(0) = (t >= 1.0 ? 10.0 : 0.0) - 50.0 * std::max(y(0) - 1.0, 0.0);
        return r;
    };
    auto ig = createIntegrator(IntegratorOptions::EulerImplicit);
    IntegratorOptions opt;
    const StepResult sr = ig->step(0.0, one(0.5), f, 1.0, opt);
    REQUIRE(sr.accepted);
    CHECK_THAT(sr.yNew(0), WithinRel(60.5 / 51.0, 1e-9));   // y = 0.5 + 10 - 50 (y - 1)
}

// ============================================================================
// The march of IntegralSolver with BDF
// ============================================================================

TEST_CASE("Integral march with BDF: rows solved at the interpolated states, evaluations counted", "[integral][bdf]") {
    // dy/dt = -k y with k = 0.5 from an algebraic block; w is an algebraic variable of the state
    auto m = buildModel(
        "y = 1 + integral(dydt, t, 0, 10)\n"
        "k^3 + k = 0.625\n"
        "dydt = -k*y\n"
        "w = 2*y + t\n"
        "$IntegralTable t:0.5 y dydt w\n");
    IntegralSolver solver(m->program, *m->ir, m->analysis, SolverOptions{});
    IntegratorOptions opt;
    opt.method = IntegratorOptions::BDF;
    opt.relTol = 1e-8;
    opt.absTol = 1e-11;
    const IntegralSolveResult r = solver.solve(opt);
    REQUIRE(r.success);
    CHECK(r.method == "BDF");
    CHECK(r.rhsEvaluations > r.totalSteps);
    REQUIRE(r.table.numRows() == 21);
    for (size_t i = 0; i < r.table.numRows(); ++i) {
        const double t = 0.5 * static_cast<double>(i);
        INFO("t = " << t);
        const double y = r.table.interpolate("y", t);
        CHECK_THAT(y, WithinRel(std::exp(-0.5 * t), 1e-6));
        // Every column of a row is that of the model solved at the interpolated state, not interpolated
        CHECK_THAT(r.table.interpolate("dydt", t), WithinRel(-0.5 * y, 1e-12));
        CHECK_THAT(r.table.interpolate("w", t), WithinRel(2.0 * y + t, 1e-12));
    }
    // Far fewer steps than rows of a fixed method at the same accuracy
    CHECK(r.totalSteps < 200);
    CHECK_THAT(r.algebraicResult.variables.at("y"), WithinRel(std::exp(-5.0), 1e-6));
}

TEST_CASE("Integral march with BDF: an exhausted budget is an error, not a truncated success", "[integral][bdf]") {
    auto m = buildModel("y = 1 + integral(dydt, t, 0, 10)\ndydt = -50*y\n");
    IntegralSolver solver(m->program, *m->ir, m->analysis, SolverOptions{});
    IntegratorOptions opt;
    opt.method = IntegratorOptions::BDF;
    opt.maxSteps = 2;
    const IntegralSolveResult r = solver.solve(opt);
    CHECK_FALSE(r.success);
    CHECK_THAT(r.errorMessage, ContainsSubstring("before reaching the upper limit"));

    // With a step in the model the budget is at least that of its steps (1000 here), whatever integralMaxSteps
    auto stepped = buildModel("y = 1 + integral(dydt, t, 0, 10, 0.01)\ndydt = -50*y\n");
    IntegralSolver steppedSolver(stepped->program, *stepped->ir, stepped->analysis, SolverOptions{});
    CHECK(steppedSolver.solve(opt).success);
}

TEST_CASE("Integral march with BDF and RK45: a state that starts at 0 on a jump of the right-hand side",
          "[integral][bdf]") {
    // dy/dt jumps by 0.2 % as y leaves 0 (the start of examples/ice_storage_tank): the error control lets the first
    // step pass only where absTol dominates the weight of y, many orders below 1e-9 of the interval
    auto m = buildModel("y = integral(dydt, t, 0, 50000)\ndydt = IF(y, 1.002E6, 1E6)\n");
    for (auto method : {IntegratorOptions::BDF, IntegratorOptions::RK45}) {
        INFO(methodToString(method));
        IntegralSolver solver(m->program, *m->ir, m->analysis, SolverOptions{});
        IntegratorOptions opt;
        opt.method = method;
        const IntegralSolveResult r = solver.solve(opt);
        REQUIRE(r.success);
        CHECK_THAT(r.algebraicResult.variables.at("y"), WithinRel(1.002e6 * 50000.0, 1e-5));
        CHECK(r.totalSteps + r.rejectedSteps < 300);
    }
}
