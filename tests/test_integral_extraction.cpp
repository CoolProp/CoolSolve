#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/integral/integral_problem.h"
#include "coolsolve/parser.h"

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

// ============================================================================
// INTEGRAL calls inside larger expressions: auxiliary states (CS-BUG-INTEGRAL-FACTOR),
// symbolic limits (CS-GAP-INTEGRAL-LIMITS)
// ============================================================================

namespace {
bool mentionsIntegral(const ExprPtr& e) {
    return astToString(e).find("integral(") != std::string::npos;
}
}  // namespace

TEST_CASE("Extraction: the plain additive form is a user state, nothing internal", "[integral][extraction][integral-factor]") {
    auto p = extractFrom(
        "y0 = 3\n"
        "y = y0 + integral(dydt, t, 0, 4, 0.1)\n"
        "dydt = -y\n");
    REQUIRE(p.valid);
    REQUIRE(p.states.size() == 1);
    CHECK_FALSE(p.states[0].auxiliary);
    CHECK(p.states[0].name == "y");
    CHECK(p.states[0].tableName == "y");
    CHECK(p.states[0].integrandVar == "dydt");
    CHECK(p.rewrittenEquations.empty());
    CHECK(p.internalNames.empty());
    CHECK(p.integrandEquations.empty());
    CHECK(p.auxEquationIds.empty());
    CHECK(p.integralEquationIds.size() == 1);
    // constant limits and step: folded, no internal variable for them
    CHECK(p.limitsConstant);
    CHECK(p.stepConstant);
    REQUIRE(p.calls.size() == 1);
    CHECK(p.calls[0].loVar.empty());
    CHECK(p.calls[0].hiVar.empty());
}

TEST_CASE("Extraction: additive forms with a sign are tracked", "[integral][extraction][integral-factor]") {
    auto userState = [](const std::string& rhs) {
        auto p = extractFrom("w = 1\nv = 2\ny = " + rhs + "\n");
        REQUIRE(p.valid);
        REQUIRE(p.states.size() == 1);
        return !p.states[0].auxiliary;
    };
    CHECK(userState("integral(v, t, 0, 1, 0.1)"));
    CHECK(userState("w + integral(v, t, 0, 1, 0.1)"));
    CHECK(userState("integral(v, t, 0, 1, 0.1) + w"));
    CHECK(userState("integral(v, t, 0, 1, 0.1) - w"));
    CHECK(userState("w - 2 + integral(v, t, 0, 1, 0.1) + 5"));
    CHECK(userState("w - (w - integral(v, t, 0, 1, 0.1))"));
    CHECK_FALSE(userState("w - integral(v, t, 0, 1, 0.1)"));
    CHECK_FALSE(userState("-integral(v, t, 0, 1, 0.1)"));
    CHECK_FALSE(userState("w*integral(v, t, 0, 1, 0.1)"));
    CHECK_FALSE(userState("integral(v, t, 0, 1, 0.1)/3.6E6"));
    CHECK_FALSE(userState("w + 2*integral(v, t, 0, 1, 0.1)"));
    CHECK_FALSE(userState("sqrt(integral(v, t, 0, 1, 0.1))"));
    CHECK_FALSE(userState("w + integral(v, t, 0, 1, 0.1)^2"));
}

TEST_CASE("Extraction: an integral inside a factor gets an auxiliary state", "[integral][extraction][integral-factor]") {
    auto p = extractFrom(
        "m_dot = 0.15\n"
        "y = m_dot*integral(v, P, 0, 10, 1)\n"
        "v = 2\n");
    REQUIRE(p.valid);
    REQUIRE(p.states.size() == 1);
    const auto& s = p.states[0];
    CHECK(s.auxiliary);
    CHECK(s.name != "y");
    CHECK(s.integrandVar == "v");
    CHECK(s.tableName == "y");              // the default table shows the left-hand side
    CHECK(p.isInternal(s.name));
    CHECK(p.internalNames == std::vector<std::string>{s.name});
    CHECK(p.integralEquationIds.empty());
    CHECK(p.auxEquationIds.size() == 1);
    REQUIRE(p.rewrittenEquations.size() == 1);
    CHECK_FALSE(mentionsIntegral(p.rewrittenEquations[0].rhs));
    CHECK(astToString(p.rewrittenEquations[0].rhs).find(s.name) != std::string::npos);
    CHECK(mentionsIntegral(p.rewrittenEquations[0].originalRhs));   // the original is untouched
    CHECK(p.integrationVar == "P");
    // y is an ordinary algebraic variable, the state is internal
    bool yAlgebraic = false;
    for (const auto& v : p.algebraicVars) if (v == "y") yAlgebraic = true;
    CHECK(yAlgebraic);
}

TEST_CASE("Extraction: several calls in one equation, an integral on the left, an integrand expression",
          "[integral][extraction][integral-factor]") {
    SECTION("two calls") {
        auto p = extractFrom("a = 1\nb = 2\ny = integral(a, t, 0, 1) - 2*integral(b, t, 0, 1)\n");
        REQUIRE(p.valid);
        REQUIRE(p.states.size() == 2);
        CHECK(p.states[0].auxiliary);
        CHECK(p.states[1].auxiliary);
        CHECK(p.states[0].name != p.states[1].name);
        CHECK(p.states[0].integrandVar == "a");
        CHECK(p.states[1].integrandVar == "b");
        CHECK(p.states[0].tableName.empty());   // two calls: no single variable to show
        REQUIRE(p.rewrittenEquations.size() == 1);
        CHECK_FALSE(mentionsIntegral(p.rewrittenEquations[0].rhs));
    }
    SECTION("call on the left-hand side") {
        auto p = extractFrom("f = 2\nintegral(f, t, 0, 1) = y\n");
        REQUIRE(p.valid);
        REQUIRE(p.states.size() == 1);
        CHECK(p.states[0].auxiliary);
        REQUIRE(p.rewrittenEquations.size() == 1);
        CHECK_FALSE(mentionsIntegral(p.rewrittenEquations[0].lhs));
        CHECK(hasIntegral(IR::fromAST(([&] { EESParser q; return q.parse("f = 2\nintegral(f, t, 0, 1) = y\n").program; })())));
    }
    SECTION("integrand that is not a variable") {
        auto p = extractFrom("y = 2*integral(x^2 + 1, x, 0, 3)\n");
        REQUIRE(p.valid);
        REQUIRE(p.states.size() == 1);
        REQUIRE(p.integrandEquations.size() == 1);
        CHECK(p.states[0].integrandVar == p.integrandEquations[0].first);
        CHECK(p.isInternal(p.states[0].integrandVar));
        CHECK(p.internalNames.size() == 2);
        // the plain form keeps its integrand expression too: an equation of its own
        auto q = extractFrom("y = 0 + integral(x^2 + 1, x, 0, 3)\n");
        REQUIRE(q.valid);
        CHECK_FALSE(q.states[0].auxiliary);
        CHECK(q.integrandEquations.size() == 1);
    }
    SECTION("array element as the integrand is a variable") {
        auto p = extractFrom("x[2] = 1\ny = 2*integral(x[2], t, 0, 3)\n");
        REQUIRE(p.valid);
        CHECK(p.states[0].integrandVar == "x[2]");
        CHECK(p.integrandEquations.empty());
    }
}

TEST_CASE("Extraction: a base that follows the integration variable or a state forces an auxiliary state",
          "[integral][extraction][integral-factor]") {
    auto aux = [](const std::string& model) {
        auto p = extractFrom(model);
        REQUIRE(p.valid);
        return p.states[0].auxiliary;
    };
    CHECK(aux("y = 2*t + integral(d, t, 0, 1)\nd = 1\n"));
    CHECK(aux("y = y/2 + integral(d, t, 0, 1)\nd = 1\n"));
    CHECK_FALSE(aux("a = 3\ny = a + integral(d, t, 0, 1)\nd = 1\n"));
    auto chain = extractFrom("y1 = 1 + integral(d1, t, 0, 1)\ny2 = y1 + integral(d2, t, 0, 1)\nd1 = 1\nd2 = 1\n");
    REQUIRE(chain.valid);
    REQUIRE(chain.states.size() == 2);
    CHECK_FALSE(chain.states[0].auxiliary);
    CHECK(chain.states[1].auxiliary);          // its base is the state y1
    CHECK(chain.states[1].tableName == "y2");
    // a state that is the integration variable (case-insensitive) is an error
    auto clash = extractFrom("T = 1 + integral(d, t, 0, 1)\nd = 1\n");
    CHECK_FALSE(clash.valid);
    CHECK(clash.errorMessage.find("integration variable") != std::string::npos);
}

TEST_CASE("Extraction: internal names never clash with the model's variables", "[integral][extraction][integral-factor]") {
    auto p = extractFrom("__integral_1 = 5\n__integrand_1 = 2\ny = __integral_1*integral(x^2, x, 0, 3)\n");
    REQUIRE(p.valid);
    REQUIRE(p.states.size() == 1);
    CHECK(p.states[0].name != "__integral_1");
    CHECK(p.states[0].integrandVar != "__integrand_1");
    CHECK(p.states[0].name.find("__integral_1") == 0);
}

TEST_CASE("Extraction: symbolic limits and steps", "[integral][extraction][integral-limits]") {
    SECTION("variables") {
        auto p = extractFrom("t_1 = 0\nt_2 = 4\ny = 1 + integral(dydt, t, t_1, t_2)\ndydt = -y\n");
        REQUIRE(p.valid);
        CHECK_FALSE(p.limitsConstant);
        CHECK(p.stepConstant);
        REQUIRE(p.calls.size() == 1);
        CHECK(!p.calls[0].loVar.empty());
        CHECK(!p.calls[0].hiVar.empty());
        CHECK(p.calls[0].stepVar.empty());
        CHECK(p.isInternal(p.calls[0].loVar));
        CHECK(p.isInternal(p.calls[0].hiVar));
    }
    SECTION("only the step is symbolic") {
        auto p = extractFrom("h = 0.1\ny = 1 + integral(dydt, t, 0, 4, h)\ndydt = -y\n");
        REQUIRE(p.valid);
        CHECK(p.limitsConstant);
        CHECK_FALSE(p.stepConstant);
        CHECK(!p.calls[0].stepVar.empty());
    }
    SECTION("constant expressions are folded as before") {
        auto p = extractFrom("y = 1 + integral(dydt, t, 2*0.5 - 1, 8/2, 0.5*0.2)\ndydt = -y\n");
        REQUIRE(p.valid);
        CHECK(p.limitsConstant);
        CHECK(p.stepConstant);
        CHECK(std::abs(p.upperLimit - 4.0) < 1e-12);
        CHECK(std::abs(p.fixedStep - 0.1) < 1e-12);
    }
    SECTION("a negative constant step is an error, 0 means automatic") {
        auto bad = extractFrom("y = 1 + integral(dydt, t, 0, 4, -1)\ndydt = -y\n");
        CHECK_FALSE(bad.valid);
        CHECK(bad.errorMessage.find("step must be positive") != std::string::npos);
        auto autoStep = extractFrom("y = 1 + integral(dydt, t, 0, 4, 0)\ndydt = -y\n");
        REQUIRE(autoStep.valid);
        CHECK(autoStep.fixedStep == 0.0);
    }
    SECTION("different explicit steps of constants: a warning, the first one is used") {
        auto p = extractFrom("y = 1 + integral(a, t, 0, 4, 0.1)\nz = 1 + integral(a, t, 0, 4, 0.2)\na = 1\n");
        REQUIRE(p.valid);
        CHECK(std::abs(p.fixedStep - 0.1) < 1e-12);
        REQUIRE(p.warnings.size() == 1);
        CHECK(p.warnings[0].find("different steps") != std::string::npos);
    }
}
