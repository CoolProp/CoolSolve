#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/integral/integral_problem.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cmath>

using namespace coolsolve;

namespace {
// Parse, build IR, analyse, and extract the IntegralProblem.
IntegralProblem extractFrom(const std::string& src) {
    EESParser parser;
    auto pr = parser.parse(src);
    REQUIRE(pr.success);
    auto ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    return extractIntegralProblem(ir, analysis);
}

bool diagContains(const IntegralProblem& p, const std::string& needle) {
    for (const auto& d : p.diagnostics)
        if (d.find(needle) != std::string::npos) return true;
    return false;
}
}  // namespace

// ============================================================================
// TDD tests for extractIntegralProblem.
// ============================================================================

TEST_CASE("Extraction: decay model", "[integral][extraction]") {
    auto p = extractFrom(
        "y = integral(dydt, t, 0, 4)\n"
        "dydt = -y\n");
    REQUIRE(p.hasIntegral);
    REQUIRE(p.valid);
    REQUIRE(p.states.size() == 1);
    CHECK(p.states[0].name == "y");
    CHECK(p.states[0].integrandVar == "dydt");
    CHECK(p.integrationVar == "t");
    CHECK(p.limitsConstant);
    CHECK(std::abs(p.lowerLimit - 0.0) < 1e-12);
    CHECK(std::abs(p.upperLimit - 4.0) < 1e-12);
    CHECK(p.fixedStep == 0.0);  // auto step
    CHECK(p.integralEquationIds.size() == 1);
}

TEST_CASE("Extraction: base expression for initial value", "[integral][extraction]") {
    // y = 1 + integral(...)  =>  baseExpr evaluates to 1 at t=t0.
    auto p = extractFrom(
        "y = 1 + integral(dydt, t, 0, 4)\n"
        "dydt = -y\n");
    REQUIRE(p.valid);
    REQUIRE(p.states.size() == 1);
    REQUIRE(p.states[0].baseExpr != nullptr);
    // baseExpr should fold to 1.0 (the integral term blanked to 0).
    CHECK((p.states[0].baseExpr->is<BinaryOp>() || p.states[0].baseExpr->is<NumberLiteral>()));
}

TEST_CASE("Extraction: two coupled states (oscillator)", "[integral][extraction]") {
    auto p = extractFrom(
        "y = 0 + integral(dydt, t, 0, 1)\n"
        "z = 1 + integral(dzdt, t, 0, 1)\n"
        "dydt = z\n"
        "dzdt = -y\n");
    REQUIRE(p.valid);
    REQUIRE(p.states.size() == 2);
    CHECK(p.stateNames == std::vector<std::string>{"y", "z"});
    CHECK(p.derivativeNames == std::vector<std::string>{"dydt", "dzdt"});
    CHECK(p.integrationVar == "t");
}

TEST_CASE("Extraction: fixed step from 5th argument", "[integral][extraction]") {
    auto p = extractFrom(
        "y = integral(dydt, t, 0, 4, 0.1)\n"
        "dydt = -y\n");
    REQUIRE(p.valid);
    CHECK(std::abs(p.fixedStep - 0.1) < 1e-12);
}

TEST_CASE("Extraction: rejects inconsistent integration variable", "[integral][extraction]") {
    auto p = extractFrom(
        "y = integral(dydt, t, 0, 4)\n"
        "z = integral(dzdt, x, 0, 4)\n"
        "dydt = -y\n"
        "dzdt = -z\n");
    CHECK_FALSE(p.valid);
    CHECK(p.errorMessage.find("integration variable") != std::string::npos);
}

TEST_CASE("Extraction: rejects inconsistent limits", "[integral][extraction]") {
    auto p = extractFrom(
        "y = integral(dydt, t, 0, 4)\n"
        "z = integral(dzdt, t, 0, 2)\n"
        "dydt = -y\n"
        "dzdt = -z\n");
    CHECK_FALSE(p.valid);
    CHECK(p.errorMessage.find("interval") != std::string::npos);
}

TEST_CASE("Extraction: rejects nested integrals", "[integral][extraction]") {
    auto p = extractFrom(
        "y = integral(integral(a, t, 0, 1), t, 0, 1)\n"
        "a = 1\n");
    CHECK_FALSE(p.valid);
    CHECK(p.errorMessage.find("Nested INTEGRAL") != std::string::npos);
}

TEST_CASE("Extraction: rejects wrong argument count", "[integral][extraction]") {
    auto p = extractFrom("y = integral(dydt, t)\n");
    CHECK_FALSE(p.valid);
    CHECK(p.errorMessage.find("4 or 5 arguments") != std::string::npos);
}

TEST_CASE("hasIntegral: false on a plain model", "[integral][extraction]") {
    EESParser parser;
    auto pr = parser.parse("x = a + b\ny = x * 2\n");
    REQUIRE(pr.success);
    auto ir = IR::fromAST(pr.program);
    CHECK_FALSE(hasIntegral(ir));
}

TEST_CASE("hasIntegral: true when an INTEGRAL call is present", "[integral][extraction]") {
    EESParser parser;
    auto pr = parser.parse("y = integral(dydt, t, 0, 1)\ndydt = -y\n");
    REQUIRE(pr.success);
    auto ir = IR::fromAST(pr.program);
    CHECK(hasIntegral(ir));
}

TEST_CASE("Extraction: classifies algebraic vs state variables", "[integral][extraction]") {
    auto p = extractFrom(
        "y = integral(dydt, t, 0, 1)\n"   // state y
        "dydt = -k * y\n"                  // algebraic: dydt
        "k = 2\n");                        // algebraic: k (parameter)
    REQUIRE(p.valid);
    // State variables excluded from algebraic list; integration var 't' too.
    CHECK(p.algebraicVars.size() >= 2);    // at least {dydt, k}
    // 'y' must not appear as algebraic; it is a state.
    bool yAlgebraic = false;
    for (const auto& v : p.algebraicVars) if (v == "y") yAlgebraic = true;
    CHECK_FALSE(yAlgebraic);
}

TEST_CASE("Extraction: non-square algebraic subsystem is flagged", "[integral][extraction]") {
    // Over-constrained: an extra algebraic equation pinning y to 5.
    auto p = extractFrom(
        "y = integral(dydt, t, 0, 1)\n"
        "dydt = -y\n"
        "y = 5\n");
    CHECK(diagContains(p, "not square"));
}
