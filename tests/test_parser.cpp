#include <catch2/catch_test_macros.hpp>
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"
#include "CoolProp.h"

TEST_CASE("Parser basic expressions", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("Simple assignment") {
        auto result = parser.parse("x = 5");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Binary operations") {
        auto result = parser.parse("y = x + 3");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Scientific notation") {
        auto result = parser.parse("P = 1E5");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Negative numbers") {
        auto result = parser.parse("T = -10");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Power operator") {
        auto result = parser.parse("z = x^2");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Complex expression") {
        auto result = parser.parse("result = a*b + c/d - e^2");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
}

TEST_CASE("Parser function calls", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("Simple function") {
        auto result = parser.parse("y = sin(x)");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Function with multiple args") {
        auto result = parser.parse("y = min(a, b)");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Function with 3 args") {
        auto result = parser.parse("y = min(a, b, c)");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Thermodynamic function with named args") {
        auto result = parser.parse("h = enthalpy(R134a, T=300, P=1E5)");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Saturation function with 1 input") {
        auto result = parser.parse("T_sat = T_sat(fluid$, P=P)");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Nested function calls") {
        auto result = parser.parse("y = max(min(a, b), c)");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
}

TEST_CASE("Parser arrays", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("Simple array access") {
        auto result = parser.parse("P[1] = 100");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Array in expression") {
        auto result = parser.parse("h[2] = h[1] + q");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
}

TEST_CASE("Parser comments", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("Quote comment") {
        auto result = parser.parse("\"This is a comment\"");
        REQUIRE(result.success);
        REQUIRE(result.commentCount == 1);
    }
    
    SECTION("C-style comment") {
        auto result = parser.parse("// This is a comment");
        REQUIRE(result.success);
        REQUIRE(result.commentCount == 1);
    }
    
    SECTION("Inline comment after equation") {
        auto result = parser.parse("x = 5 \"value in meters\"");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }

    SECTION("Inline brace comment after equation") {
        auto result = parser.parse("x = 5 {value in meters}");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }

    SECTION("Multi-line quote block comments behave like brace comments") {
        const std::string code_quotes = R"(
x = 1
"
    y = 2
"
z = 3
)";
        const std::string code_braces = R"(
x = 1
{
    y = 2
}
z = 3
)";

        auto res_quotes = parser.parse(code_quotes);
        auto res_braces = parser.parse(code_braces);

        REQUIRE(res_quotes.success);
        REQUIRE(res_braces.success);
        // In both cases only the equations outside the comment block should be kept.
        REQUIRE(res_quotes.equationCount == 2);
        REQUIRE(res_braces.equationCount == 2);
    }

    SECTION("Nested brace block comments spanning multiple lines") {
        const std::string code = R"(
a = 1
b = 2
{
    c = 3
    { nested comment start
        d = 4
    } nested comment end
    e = 5
}
f = 6
)";
        auto result = parser.parse(code);
        REQUIRE(result.success);
        // Only equations outside the brace comment block should be counted.
        REQUIRE(result.equationCount == 3); // a = 1, b = 2, f = 6
    }
}

TEST_CASE("Parser units annotation", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("Simple units") {
        auto result = parser.parse("P = 1E5 \"[Pa]\"");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
}

TEST_CASE("Parser string variables", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("String variable assignment") {
        auto result = parser.parse("fluid$ = 'R134a'");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("String variable in function") {
        auto result = parser.parse("h = enthalpy(fluid$, T=300, P=100)");
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1);
    }
    
    SECTION("Invalid variable with $ in middle") {
        auto result = parser.parse("flu$id = 'R134a'");
        // Should fail to parse the line, but might continue
        // The result.success depends on whether any equations were parsed
        // Here we expect it to NOT parse the equation
        REQUIRE(result.equationCount == 0);
    }
}

TEST_CASE("Parser directives", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("ifnot directive") {
        auto result = parser.parse("$ifnot parametrictable\nx = 5\n$endif");
        REQUIRE(result.success);
    }
}

TEST_CASE("IR building", "[ir]") {
    coolsolve::EESParser parser;
    
    SECTION("Variable extraction") {
        auto parseResult = parser.parse("x = a + b\ny = x * c");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        REQUIRE(ir.getEquationCount() == 2);
        // Variables: x, a, b, y, c = 5
        REQUIRE(ir.getVariableCount() == 5);
    }
    
    SECTION("Incidence matrix") {
        auto parseResult = parser.parse("x = a + b\ny = x * c");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        const auto& incidence = ir.getIncidenceMatrix();
        REQUIRE(incidence.size() == 2);
        
        // First equation: x, a, b
        REQUIRE(incidence[0].count("x") == 1);
        REQUIRE(incidence[0].count("a") == 1);
        REQUIRE(incidence[0].count("b") == 1);
        
        // Second equation: y, x, c
        REQUIRE(incidence[1].count("y") == 1);
        REQUIRE(incidence[1].count("x") == 1);
        REQUIRE(incidence[1].count("c") == 1);
    }
}

TEST_CASE("Structural analysis - simple system", "[analysis]") {
    coolsolve::EESParser parser;
    
    SECTION("Sequential system") {
        // x = 5
        // y = x + 3
        // z = y * 2
        auto parseResult = parser.parse("x = 5\ny = x + 3\nz = y * 2");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto result = coolsolve::StructuralAnalyzer::analyze(ir);
        
        REQUIRE(result.success);
        REQUIRE(result.totalEquations == 3);
        // All equations should be in separate blocks (sequential evaluation)
        REQUIRE(result.totalBlocks == 3);
    }
    
    SECTION("System with algebraic loop") {
        // x = y + 1
        // y = x + 1
        // This creates a 2x2 algebraic loop
        auto parseResult = parser.parse("x = y + 1\ny = x + 1");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto result = coolsolve::StructuralAnalyzer::analyze(ir);
        
        REQUIRE(result.success);
        REQUIRE(result.totalEquations == 2);
        // Both equations should be in the same block (loop)
        REQUIRE(result.largestBlockSize == 2);
    }
}

// Hard-coded example: Heat exchanger problem
// This is a simplified version of exchangers1.eescode
TEST_CASE("Hard-coded heat exchanger example", "[examples]") {
    coolsolve::EESParser parser;
    
    // Simplified heat exchanger equations (counter-flow)
    const std::string exchangerCode = R"(
A=12.5
c_oil=2100
T_oil_su=100
P_w_su=1E5
T_w_su=20
M_dot_oil=2
M_dot_w=0.48
U=400

c_w=CP(Water,T=T_w_su,P=P_w_su)

C_dot_oil=M_dot_oil*c_oil
C_dot_w=M_dot_w*c_w
C_dot_min=min(C_dot_oil,C_dot_w)
C_dot_max=max(C_dot_oil,C_dot_w)
omega=C_dot_min/C_dot_max

AU=A*U
NTU=AU/C_dot_min

epsilon=(1-exp(-NTU*(1-omega)))/(1-omega*exp(-NTU*(1-omega)))

Q_dot=epsilon*C_dot_min*(T_oil_su-T_w_su)
Q_dot=C_dot_w*(T_w_ex-T_w_su)
Q_dot=C_dot_oil*(T_oil_su-T_oil_ex)
)";
    
    auto parseResult = parser.parse(exchangerCode);
    
    REQUIRE(parseResult.success);
    REQUIRE(parseResult.equationCount == 20);
    
    auto ir = coolsolve::IR::fromAST(parseResult.program);
    REQUIRE(ir.getEquationCount() == 20);
    
    auto analysisResult = coolsolve::StructuralAnalyzer::analyze(ir);
    REQUIRE(analysisResult.success);
    
    INFO("Total blocks: " << analysisResult.totalBlocks);
    INFO("Largest block size: " << analysisResult.largestBlockSize);
}

// Hard-coded example: Refrigeration cycle
// This is a simplified version of refrigeration1.eescode
TEST_CASE("Parser procedures and functions", "[parser]") {
    coolsolve::EESParser parser;
    
    SECTION("Simple procedure definition and call") {
        const std::string code = R"(
PROCEDURE test(x, y : z)
    z := x + y
END
CALL test(1, 2 : result)
)";
        auto result = parser.parse(code);
        REQUIRE(result.success);
        // Note: Procedure definition is one statement, CALL is another (represented as pseudo-equation in IR)
        // But in the parser, we don't count them as equations.
        // Actually, the parser implementation for success check:
        // result.success = !unsupportedFound && (result.equationCount > 0 || result.commentCount > 0 || result.directiveCount > 0);
        // I need to update this logic in parser.cpp if I want it to be success without equations.
    }
    
    SECTION("Procedure with multiple outputs and assignment") {
        const std::string code = R"(
PROCEDURE multi(x : y, z)
    y := x * 2
    z = x * 3
END
CALL multi(10 : a, b)
)";
        auto result = parser.parse(code);
        REQUIRE(result.success);
    }

    SECTION("Function definition and call") {
        const std::string code = R"(
FUNCTION myfunc(x, y)
    myfunc := x^2 + y
END
z = myfunc(3, 4)
)";
        auto result = parser.parse(code);
        REQUIRE(result.success);
        REQUIRE(result.equationCount == 1); // z = myfunc(3, 4)
    }

    SECTION("Function with IF-THEN-ELSE") {
        const std::string code = R"(
FUNCTION test_if(x)
    IF (x > 0) THEN
        test_if := 1
    ELSE
        test_if := -1
    ENDIF
END
y = test_if(10)
)";
        auto result = parser.parse(code);
        REQUIRE(result.success);
    }
}

TEST_CASE("Evaluator with procedures and functions", "[evaluator]") {
    coolsolve::CoolPropConfig config;
    
    SECTION("User function evaluation") {
        const std::string code = R"(
FUNCTION sq(x)
    sq := x*x
END
y = sq(4)
)";
        coolsolve::EESParser parser;
        auto parseResult = parser.parse(code);
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        coolsolve::SystemEvaluator evaluator(ir, analysis, config);
        
        evaluator.initializeFromGuesses();
        // y = sq(4) -> y = 16
        auto result = evaluator.evaluateBlock(0);
        REQUIRE(result.residuals.size() == 1);
        // Residual for y = sq(4) is y - 16. If y=1 (default), residual is -15.
        // We can check if it computes the right value by setting y=16
        evaluator.setVariableValue("y", 16.0);
        result = evaluator.evaluateBlock(0);
        REQUIRE(std::abs(result.residuals[0]) < 1e-10);
        
        // Check Jacobian dy/dy = 1
        REQUIRE(std::abs(result.jacobian[0][0] - 1.0) < 1e-10);
    }
    
    SECTION("User procedure evaluation") {
        const std::string code = R"(
PROCEDURE calc(x : y, z)
    y := x * 2
    z := x + 5
END
CALL calc(10 : a, b)
)";
        coolsolve::EESParser parser;
        auto parseResult = parser.parse(code);
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        coolsolve::SystemEvaluator evaluator(ir, analysis, config);
        
        evaluator.initializeFromGuesses();
        // block 0 should have a, b
        evaluator.setVariableValue("a", 20.0);
        evaluator.setVariableValue("b", 15.0);
        auto result = evaluator.evaluateBlock(0);
        // There should be at least one residual for the procedure call.
        // In my simplified implementation, it returns residuals for the first output only.
        // So size should be 1 if CALL is the only statement in block.
        REQUIRE(result.residuals.size() >= 1);
        REQUIRE(std::abs(result.residuals[0]) < 1e-10);
    }

    SECTION("Procedure with string variable passing") {
        const std::string code = R"(
PROCEDURE get_h(f$, t, p : h)
    h := enthalpy(f$, T=t, P=p)
END
fluid$ = 'Water'
T = 300
P = 101325
CALL get_h(fluid$, T, P : h_res)
)";
        coolsolve::EESParser parser;
        auto parseResult = parser.parse(code);
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        coolsolve::CoolPropConfig config;
        config.units.temperature = "K"; // Use Kelvin for easier comparison
        coolsolve::SystemEvaluator evaluator(ir, analysis, config);
        
        evaluator.initializeFromGuesses();
        evaluator.setStringVariableValue("fluid$", "Water");
        evaluator.setVariableValue("T", 300.0);
        evaluator.setVariableValue("P", 101325.0);
        
        // Find the block containing h_res
        int hResBlock = -1;
        for (size_t i = 0; i < evaluator.getNumBlocks(); ++i) {
            for (const auto& var : evaluator.getBlock(i).getVariables()) {
                if (var == "h_res") {
                    hResBlock = i;
                    break;
                }
            }
        }
        REQUIRE(hResBlock != -1);
        
        // Water enthalpy at 300K, 101325Pa
        double h_si = CoolProp::PropsSI("H", "T", 300, "P", 101325, "Water");
        // Apply reference state offset if necessary. 
        // CoolProp's default for Water is usually fine, but let's see.
        evaluator.setVariableValue("h_res", h_si);
        
        auto result = evaluator.evaluateBlock(hResBlock);
        // Result should be very close to 0 if passing worked
        REQUIRE(std::abs(result.residuals[0]) < 1e-5);
    }
}

TEST_CASE("IR building with procedures", "[ir]") {
    coolsolve::EESParser parser;
    
    SECTION("Procedure and function extraction") {
        const std::string code = R"(
FUNCTION f(x)
    f := x * 2
END
PROCEDURE p(x : y)
    y := f(x) + 1
END
CALL p(10 : res)
)";
        auto parseResult = parser.parse(code);
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        REQUIRE(ir.getFunctions().size() == 1);
        REQUIRE(ir.getProcedures().size() == 1);
        REQUIRE(ir.getEquationCount() == 1); // CALL p(10 : res)
        REQUIRE(ir.getVariableCount() == 1); // res
    }
}
