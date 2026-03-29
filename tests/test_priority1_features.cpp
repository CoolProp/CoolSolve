/**
 * Tests for Priority 1 features:
 * - Inverse hyperbolic functions (ARCSINH, ARCCOSH, ARCTANH)
 * - Rounding functions (CEIL, FLOOR, ROUND, TRUNC, SIGN)
 * - MOD function
 * - Inline IF function
 * - Aggregation functions (SUM, SUM2D, AVERAGE, PRODUCT, STDDEV)
 * - CONVERT / CONVERTTEMP
 * - DUPLICATE loop
 * - REPEAT-UNTIL loop
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "coolsolve/autodiff_node.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// Helper: parse a single equation, evaluate the RHS with given variable values
static coolsolve::ADValue evalRHS(const std::string& source,
                                  const std::vector<std::pair<std::string, double>>& vars) {
    coolsolve::EESParser parser;
    auto parseResult = parser.parse(source);
    REQUIRE(parseResult.success);
    auto ir = coolsolve::IR::fromAST(parseResult.program);
    const auto& equations = ir.getEquations();
    REQUIRE(!equations.empty());

    size_t n = vars.size();
    coolsolve::ExpressionEvaluator eval(n);
    for (size_t i = 0; i < n; ++i) {
        eval.setVariable(vars[i].first,
                         coolsolve::ADValue::independent(vars[i].second, i, n));
    }
    return eval.evaluate(equations[0].rhs);
}

// ============================================================================
// Inverse Hyperbolic Functions
// ============================================================================

TEST_CASE("Inverse hyperbolic functions", "[priority1][math]") {
    using namespace coolsolve;

    SECTION("ARCSINH") {
        std::vector<ADValue> args = { ADValue::independent(1.0, 0, 1) };
        auto result = evaluateStandardFunction("arcsinh", args);
        REQUIRE_THAT(result.value, WithinAbs(std::asinh(1.0), 1e-12));
        // d(asinh(x))/dx = 1/sqrt(x^2+1)
        REQUIRE_THAT(result.gradient[0], WithinAbs(1.0 / std::sqrt(2.0), 1e-12));
    }

    SECTION("ARCCOSH") {
        std::vector<ADValue> args = { ADValue::independent(2.0, 0, 1) };
        auto result = evaluateStandardFunction("arccosh", args);
        REQUIRE_THAT(result.value, WithinAbs(std::acosh(2.0), 1e-12));
        // d(acosh(x))/dx = 1/sqrt(x^2-1)
        REQUIRE_THAT(result.gradient[0], WithinAbs(1.0 / std::sqrt(3.0), 1e-12));
    }

    SECTION("ARCTANH") {
        std::vector<ADValue> args = { ADValue::independent(0.5, 0, 1) };
        auto result = evaluateStandardFunction("arctanh", args);
        REQUIRE_THAT(result.value, WithinAbs(std::atanh(0.5), 1e-12));
        // d(atanh(x))/dx = 1/(1-x^2)
        REQUIRE_THAT(result.gradient[0], WithinAbs(1.0 / (1.0 - 0.25), 1e-12));
    }
}

// ============================================================================
// Rounding Functions
// ============================================================================

TEST_CASE("Rounding functions", "[priority1][math]") {
    using namespace coolsolve;

    SECTION("CEIL") {
        std::vector<ADValue> args = { ADValue::constant(2.3, 1) };
        auto result = evaluateStandardFunction("ceil", args);
        REQUIRE(result.value == 3.0);
    }

    SECTION("FLOOR") {
        std::vector<ADValue> args = { ADValue::constant(2.7, 1) };
        auto result = evaluateStandardFunction("floor", args);
        REQUIRE(result.value == 2.0);
    }

    SECTION("ROUND") {
        std::vector<ADValue> args = { ADValue::constant(2.5, 1) };
        auto result = evaluateStandardFunction("round", args);
        // std::round(2.5) rounds away from zero = 3.0
        REQUIRE(result.value == 3.0);

        args = { ADValue::constant(2.4, 1) };
        result = evaluateStandardFunction("round", args);
        REQUIRE(result.value == 2.0);
    }

    SECTION("TRUNC") {
        std::vector<ADValue> args = { ADValue::constant(2.9, 1) };
        auto result = evaluateStandardFunction("trunc", args);
        REQUIRE(result.value == 2.0);

        args = { ADValue::constant(-2.9, 1) };
        result = evaluateStandardFunction("trunc", args);
        REQUIRE(result.value == -2.0);
    }

    SECTION("SIGN") {
        std::vector<ADValue> args = { ADValue::constant(5.0, 1) };
        REQUIRE(evaluateStandardFunction("sign", args).value == 1.0);

        args = { ADValue::constant(-3.0, 1) };
        REQUIRE(evaluateStandardFunction("sign", args).value == -1.0);

        args = { ADValue::constant(0.0, 1) };
        REQUIRE(evaluateStandardFunction("sign", args).value == 0.0);
    }
}

// ============================================================================
// MOD Function
// ============================================================================

TEST_CASE("MOD function", "[priority1][math]") {
    using namespace coolsolve;
    std::vector<ADValue> args = {
        ADValue::constant(7.0, 1),
        ADValue::constant(3.0, 1)
    };
    auto result = evaluateStandardFunction("mod", args);
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-12));

    args = { ADValue::constant(10.0, 1), ADValue::constant(4.0, 1) };
    result = evaluateStandardFunction("mod", args);
    REQUIRE_THAT(result.value, WithinAbs(2.0, 1e-12));
}

// ============================================================================
// Inline IF Function
// ============================================================================

TEST_CASE("Inline IF function", "[priority1][math]") {
    using namespace coolsolve;

    SECTION("Condition true (positive)") {
        std::vector<ADValue> args = {
            ADValue::constant(1.0, 1),   // condition > 0 → true
            ADValue::constant(10.0, 1),  // true_value
            ADValue::constant(20.0, 1)   // false_value
        };
        auto result = evaluateStandardFunction("if", args);
        REQUIRE(result.value == 10.0);
    }

    SECTION("Condition false (negative)") {
        std::vector<ADValue> args = {
            ADValue::constant(-1.0, 1),
            ADValue::constant(10.0, 1),
            ADValue::constant(20.0, 1)
        };
        auto result = evaluateStandardFunction("if", args);
        REQUIRE(result.value == 20.0);
    }

    SECTION("Condition false (zero)") {
        std::vector<ADValue> args = {
            ADValue::constant(0.0, 1),
            ADValue::constant(10.0, 1),
            ADValue::constant(20.0, 1)
        };
        auto result = evaluateStandardFunction("if", args);
        REQUIRE(result.value == 20.0);
    }
}

// ============================================================================
// Aggregation Functions (SUM, AVERAGE, PRODUCT, STDDEV)
// ============================================================================

TEST_CASE("Aggregation functions", "[priority1][math]") {
    using namespace coolsolve;

    SECTION("SUM") {
        std::vector<ADValue> args = {
            ADValue::constant(1.0, 1),
            ADValue::constant(2.0, 1),
            ADValue::constant(3.0, 1),
            ADValue::constant(4.0, 1)
        };
        auto result = evaluateStandardFunction("sum", args);
        REQUIRE_THAT(result.value, WithinAbs(10.0, 1e-12));
    }

    SECTION("SUM2D") {
        std::vector<ADValue> args = {
            ADValue::constant(5.0, 1),
            ADValue::constant(10.0, 1),
            ADValue::constant(15.0, 1)
        };
        auto result = evaluateStandardFunction("sum2d", args);
        REQUIRE_THAT(result.value, WithinAbs(30.0, 1e-12));
    }

    SECTION("AVERAGE") {
        std::vector<ADValue> args = {
            ADValue::constant(2.0, 1),
            ADValue::constant(4.0, 1),
            ADValue::constant(6.0, 1)
        };
        auto result = evaluateStandardFunction("average", args);
        REQUIRE_THAT(result.value, WithinAbs(4.0, 1e-12));
    }

    SECTION("PRODUCT") {
        std::vector<ADValue> args = {
            ADValue::constant(2.0, 1),
            ADValue::constant(3.0, 1),
            ADValue::constant(4.0, 1)
        };
        auto result = evaluateStandardFunction("product", args);
        REQUIRE_THAT(result.value, WithinAbs(24.0, 1e-12));
    }

    SECTION("STDDEV") {
        // stddev of {2, 4, 4, 4, 5, 5, 7, 9} = sqrt(32/8) = 2.0
        std::vector<ADValue> args = {
            ADValue::constant(2.0, 1),
            ADValue::constant(4.0, 1),
            ADValue::constant(4.0, 1),
            ADValue::constant(4.0, 1),
            ADValue::constant(5.0, 1),
            ADValue::constant(5.0, 1),
            ADValue::constant(7.0, 1),
            ADValue::constant(9.0, 1)
        };
        auto result = evaluateStandardFunction("stddev", args);
        // Mean = 40/8 = 5, sum of (x-5)^2 = 9+1+1+1+0+0+4+16 = 32
        // stddev = sqrt(32/8) = sqrt(4) = 2
        REQUIRE_THAT(result.value, WithinAbs(2.0, 1e-10));
    }
}

// ============================================================================
// CONVERT / CONVERTTEMP (via parser + evaluator integration)
// ============================================================================

TEST_CASE("CONVERT function", "[priority1][convert]") {
    // CONVERT('kJ', 'J') should return 1000
    auto result = evalRHS("y = CONVERT('kJ', 'J')", {});
    REQUIRE_THAT(result.value, WithinAbs(1000.0, 1e-6));
}

TEST_CASE("CONVERTTEMP function", "[priority1][convert]") {
    // Convert 100 °C to K → 373.15
    auto result = evalRHS("y = CONVERTTEMP('C', 'K', 100)", {});
    REQUIRE_THAT(result.value, WithinAbs(373.15, 0.01));

    // Convert 32 °F to °C → 0
    result = evalRHS("y = CONVERTTEMP('F', 'C', 32)", {});
    REQUIRE_THAT(result.value, WithinAbs(0.0, 0.01));
}

// ============================================================================
// DUPLICATE loop (parser + AST)
// ============================================================================

TEST_CASE("DUPLICATE parsing", "[priority1][duplicate]") {
    coolsolve::EESParser parser;
    std::string source = R"(
DUPLICATE i = 1, 3
    x[i] = i * 2
END
)";
    auto parseResult = parser.parse(source);
    REQUIRE(parseResult.success);

    // Find the Duplicate statement
    bool foundDuplicate = false;
    for (const auto& stmt : parseResult.program.statements) {
        if (stmt->is<coolsolve::Duplicate>()) {
            const auto& dup = stmt->as<coolsolve::Duplicate>();
            REQUIRE(dup.iteratorVar == "i");
            REQUIRE(dup.body.size() == 1);
            foundDuplicate = true;
        }
    }
    REQUIRE(foundDuplicate);
}

// ============================================================================
// REPEAT-UNTIL loop (parser + AST)
// ============================================================================

TEST_CASE("REPEAT-UNTIL parsing", "[priority1][repeat]") {
    coolsolve::EESParser parser;
    std::string source = R"(
REPEAT
    x = x + 1
UNTIL(x > 10)
)";
    auto parseResult = parser.parse(source);
    REQUIRE(parseResult.success);

    bool foundRepeat = false;
    for (const auto& stmt : parseResult.program.statements) {
        if (stmt->is<coolsolve::RepeatUntil>()) {
            const auto& ru = stmt->as<coolsolve::RepeatUntil>();
            REQUIRE(ru.body.size() == 1);
            REQUIRE(ru.condition != nullptr);
            foundRepeat = true;
        }
    }
    REQUIRE(foundRepeat);
}

// ============================================================================
// Integration: functions used in parsed expressions
// ============================================================================

TEST_CASE("Math functions via parser", "[priority1][integration]") {
    SECTION("ceil via parser") {
        auto result = evalRHS("y = ceil(2.3)", {});
        REQUIRE(result.value == 3.0);
    }

    SECTION("floor via parser") {
        auto result = evalRHS("y = floor(2.7)", {});
        REQUIRE(result.value == 2.0);
    }

    SECTION("round via parser") {
        auto result = evalRHS("y = round(2.6)", {});
        REQUIRE(result.value == 3.0);
    }

    SECTION("trunc via parser") {
        auto result = evalRHS("y = trunc(-2.9)", {});
        REQUIRE(result.value == -2.0);
    }

    SECTION("sign via parser") {
        auto result = evalRHS("y = sign(-5)", {});
        REQUIRE(result.value == -1.0);
    }

    SECTION("mod via parser") {
        auto result = evalRHS("y = mod(7, 3)", {});
        REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-12));
    }

    SECTION("IF via parser") {
        auto result = evalRHS("y = IF(1, 42, 99)", {});
        REQUIRE(result.value == 42.0);

        result = evalRHS("y = IF(-1, 42, 99)", {});
        REQUIRE(result.value == 99.0);
    }

    SECTION("SUM via parser") {
        auto result = evalRHS("y = SUM(1, 2, 3, 4)", {});
        REQUIRE_THAT(result.value, WithinAbs(10.0, 1e-12));
    }

    SECTION("AVERAGE via parser") {
        auto result = evalRHS("y = AVERAGE(2, 4, 6)", {});
        REQUIRE_THAT(result.value, WithinAbs(4.0, 1e-12));
    }

    SECTION("PRODUCT via parser") {
        auto result = evalRHS("y = PRODUCT(2, 3, 4)", {});
        REQUIRE_THAT(result.value, WithinAbs(24.0, 1e-12));
    }

    SECTION("arcsinh via parser") {
        auto result = evalRHS("y = arcsinh(1)", {});
        REQUIRE_THAT(result.value, WithinAbs(std::asinh(1.0), 1e-10));
    }

    SECTION("arccosh via parser") {
        auto result = evalRHS("y = arccosh(2)", {});
        REQUIRE_THAT(result.value, WithinAbs(std::acosh(2.0), 1e-10));
    }

    SECTION("arctanh via parser") {
        auto result = evalRHS("y = arctanh(0.5)", {});
        REQUIRE_THAT(result.value, WithinAbs(std::atanh(0.5), 1e-10));
    }
}

// ============================================================================
// DUPLICATE expansion in IR (top-level DUPLICATE generates equations)
// ============================================================================

TEST_CASE("DUPLICATE expansion in IR", "[priority1][duplicate][ir]") {
    using namespace coolsolve;
    EESParser parser;

    SECTION("Simple DUPLICATE expands to individual equations") {
        std::string source = R"(
DUPLICATE i = 1, 3
    x[i] = i * 2
END
)";
        auto parseResult = parser.parse(source);
        REQUIRE(parseResult.success);

        auto ir = IR::fromAST(parseResult.program);
        // Should produce 3 equations: x[1]=1*2, x[2]=2*2, x[3]=3*2
        REQUIRE(ir.getEquationCount() == 3);
        // Variables x[1], x[2], x[3] should be registered
        REQUIRE(ir.getVariables().count("x[1]") == 1);
        REQUIRE(ir.getVariables().count("x[2]") == 1);
        REQUIRE(ir.getVariables().count("x[3]") == 1);
    }

    SECTION("Nested DUPLICATE expands correctly") {
        std::string source = R"(
DUPLICATE i = 1, 2
    DUPLICATE j = 1, 2
        A[i,j] = i * 10 + j
    END
END
)";
        auto parseResult = parser.parse(source);
        REQUIRE(parseResult.success);

        auto ir = IR::fromAST(parseResult.program);
        // 2*2 = 4 equations
        REQUIRE(ir.getEquationCount() == 4);
    }

    SECTION("DUPLICATE inside FUNCTION body parsed correctly") {
        std::string source = R"(
FUNCTION MySum(N)
    MySum := 0
    DUPLICATE i = 1, N
        MySum := MySum + i
    END
END

result = MySum(5)
)";
        auto parseResult = parser.parse(source);
        REQUIRE(parseResult.success);

        // The DUPLICATE should be inside the function body as a Duplicate AST node
        bool foundFunc = false;
        for (const auto& stmt : parseResult.program.statements) {
            if (stmt->is<FunctionDefinition>()) {
                const auto& func = stmt->as<FunctionDefinition>();
                if (func.name == "MySum") {
                    foundFunc = true;
                    // Body should contain: Assignment + Duplicate
                    bool foundDup = false;
                    for (const auto& bodyStmt : func.body) {
                        if (bodyStmt->is<Duplicate>()) {
                            foundDup = true;
                            const auto& dup = bodyStmt->as<Duplicate>();
                            REQUIRE(dup.iteratorVar == "i");
                        }
                    }
                    REQUIRE(foundDup);
                }
            }
        }
        REQUIRE(foundFunc);
    }
}
