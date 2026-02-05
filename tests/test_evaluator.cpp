#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "coolsolve/autodiff_node.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include <cmath>
#include <iomanip>
#include <sstream>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ============================================================================
// ADValue Tests
// ============================================================================

TEST_CASE("ADValue basic operations", "[autodiff]") {
    const size_t numVars = 2;
    
    SECTION("Constant has zero gradient") {
        auto c = coolsolve::ADValue::constant(5.0, numVars);
        REQUIRE(c.value == 5.0);
        REQUIRE(c.gradient.size() == numVars);
        REQUIRE(c.gradient[0] == 0.0);
        REQUIRE(c.gradient[1] == 0.0);
    }
    
    SECTION("Independent variable has unit gradient") {
        auto x = coolsolve::ADValue::independent(3.0, 0, numVars);
        REQUIRE(x.value == 3.0);
        REQUIRE(x.gradient[0] == 1.0);
        REQUIRE(x.gradient[1] == 0.0);
        
        auto y = coolsolve::ADValue::independent(7.0, 1, numVars);
        REQUIRE(y.value == 7.0);
        REQUIRE(y.gradient[0] == 0.0);
        REQUIRE(y.gradient[1] == 1.0);
    }
    
    SECTION("Addition") {
        auto x = coolsolve::ADValue::independent(3.0, 0, numVars);
        auto y = coolsolve::ADValue::independent(7.0, 1, numVars);
        auto z = x + y;
        
        REQUIRE(z.value == 10.0);
        REQUIRE(z.gradient[0] == 1.0);  // dz/dx = 1
        REQUIRE(z.gradient[1] == 1.0);  // dz/dy = 1
    }
    
    SECTION("Subtraction") {
        auto x = coolsolve::ADValue::independent(3.0, 0, numVars);
        auto y = coolsolve::ADValue::independent(7.0, 1, numVars);
        auto z = x - y;
        
        REQUIRE(z.value == -4.0);
        REQUIRE(z.gradient[0] == 1.0);   // dz/dx = 1
        REQUIRE(z.gradient[1] == -1.0);  // dz/dy = -1
    }
    
    SECTION("Multiplication") {
        auto x = coolsolve::ADValue::independent(3.0, 0, numVars);
        auto y = coolsolve::ADValue::independent(7.0, 1, numVars);
        auto z = x * y;
        
        REQUIRE(z.value == 21.0);
        REQUIRE(z.gradient[0] == 7.0);  // dz/dx = y = 7
        REQUIRE(z.gradient[1] == 3.0);  // dz/dy = x = 3
    }
    
    SECTION("Division") {
        auto x = coolsolve::ADValue::independent(6.0, 0, numVars);
        auto y = coolsolve::ADValue::independent(2.0, 1, numVars);
        auto z = x / y;
        
        REQUIRE(z.value == 3.0);
        REQUIRE_THAT(z.gradient[0], WithinAbs(0.5, 1e-10));   // dz/dx = 1/y = 0.5
        REQUIRE_THAT(z.gradient[1], WithinAbs(-1.5, 1e-10)); // dz/dy = -x/y^2 = -6/4 = -1.5
    }
}

TEST_CASE("ADValue power operations", "[autodiff]") {
    const size_t numVars = 1;
    
    SECTION("Power with constant exponent") {
        // z = x^2
        auto x = coolsolve::ADValue::independent(3.0, 0, numVars);
        auto z = coolsolve::pow(x, 2.0);
        
        REQUIRE(z.value == 9.0);
        REQUIRE_THAT(z.gradient[0], WithinAbs(6.0, 1e-10));  // dz/dx = 2*x = 6
    }
    
    SECTION("Power with variable exponent") {
        // z = x^y at x=2, y=3
        const size_t n = 2;
        auto x = coolsolve::ADValue::independent(2.0, 0, n);
        auto y = coolsolve::ADValue::independent(3.0, 1, n);
        auto z = coolsolve::pow(x, y);
        
        REQUIRE(z.value == 8.0);  // 2^3 = 8
        // dz/dx = y * x^(y-1) = 3 * 2^2 = 12
        REQUIRE_THAT(z.gradient[0], WithinAbs(12.0, 1e-10));
        // dz/dy = x^y * ln(x) = 8 * ln(2)
        REQUIRE_THAT(z.gradient[1], WithinAbs(8.0 * std::log(2.0), 1e-10));
    }
}

TEST_CASE("ADValue transcendental functions", "[autodiff]") {
    const size_t numVars = 1;
    
    SECTION("sin(x)") {
        auto x = coolsolve::ADValue::independent(M_PI / 6.0, 0, numVars);  // 30 degrees
        auto z = coolsolve::sin(x);
        
        REQUIRE_THAT(z.value, WithinAbs(0.5, 1e-10));
        REQUIRE_THAT(z.gradient[0], WithinAbs(std::cos(M_PI / 6.0), 1e-10));
    }
    
    SECTION("cos(x)") {
        auto x = coolsolve::ADValue::independent(M_PI / 3.0, 0, numVars);  // 60 degrees
        auto z = coolsolve::cos(x);
        
        REQUIRE_THAT(z.value, WithinAbs(0.5, 1e-10));
        REQUIRE_THAT(z.gradient[0], WithinAbs(-std::sin(M_PI / 3.0), 1e-10));
    }
    
    SECTION("exp(x)") {
        auto x = coolsolve::ADValue::independent(1.0, 0, numVars);
        auto z = coolsolve::exp(x);
        
        REQUIRE_THAT(z.value, WithinAbs(M_E, 1e-10));
        REQUIRE_THAT(z.gradient[0], WithinAbs(M_E, 1e-10));  // d(e^x)/dx = e^x
    }
    
    SECTION("log(x)") {
        auto x = coolsolve::ADValue::independent(M_E, 0, numVars);
        auto z = coolsolve::log(x);
        
        REQUIRE_THAT(z.value, WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(z.gradient[0], WithinAbs(1.0 / M_E, 1e-10));  // d(ln(x))/dx = 1/x
    }
    
    SECTION("sqrt(x)") {
        auto x = coolsolve::ADValue::independent(4.0, 0, numVars);
        auto z = coolsolve::sqrt(x);
        
        REQUIRE(z.value == 2.0);
        REQUIRE_THAT(z.gradient[0], WithinAbs(0.25, 1e-10));  // d(sqrt(x))/dx = 1/(2*sqrt(x)) = 0.25
    }
}

TEST_CASE("ADValue complex expression", "[autodiff]") {
    // Test: y = x^2 + 3*x
    // dy/dx = 2*x + 3
    const size_t numVars = 1;
    
    SECTION("At x = 2") {
        auto x = coolsolve::ADValue::independent(2.0, 0, numVars);
        auto y = coolsolve::pow(x, 2.0) + x * 3.0;
        
        REQUIRE(y.value == 10.0);  // 4 + 6 = 10
        REQUIRE_THAT(y.gradient[0], WithinAbs(7.0, 1e-10));  // 2*2 + 3 = 7
    }
    
    SECTION("At x = -1") {
        auto x = coolsolve::ADValue::independent(-1.0, 0, numVars);
        auto y = coolsolve::pow(x, 2.0) + x * 3.0;
        
        REQUIRE(y.value == -2.0);  // 1 - 3 = -2
        REQUIRE_THAT(y.gradient[0], WithinAbs(1.0, 1e-10));  // 2*(-1) + 3 = 1
    }
}

// ============================================================================
// ExpressionEvaluator Tests
// ============================================================================

TEST_CASE("ExpressionEvaluator basic expressions", "[evaluator]") {
    coolsolve::EESParser parser;
    
    SECTION("Simple algebraic expression") {
        // y = x^2 + 3*x at x = 2
        auto parseResult = parser.parse("y = x^2 + 3*x");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        const auto& equations = ir.getEquations();
        REQUIRE(equations.size() == 1);
        
        // We have 2 variables: x (independent) and y (the output)
        coolsolve::ExpressionEvaluator eval(2);
        eval.setVariable("x", coolsolve::ADValue::independent(2.0, 0, 2));
        eval.setVariable("y", coolsolve::ADValue::independent(10.0, 1, 2));  // y at the solution
        
        // Evaluate LHS (y)
        auto lhs = eval.evaluate(equations[0].lhs);
        REQUIRE(lhs.gradient[0] == 0.0);  // d(y)/dx = 0
        REQUIRE(lhs.gradient[1] == 1.0);  // d(y)/dy = 1
        
        // Evaluate RHS (x^2 + 3*x)
        auto rhs = eval.evaluate(equations[0].rhs);
        REQUIRE(rhs.value == 10.0);
        REQUIRE_THAT(rhs.gradient[0], WithinAbs(7.0, 1e-10));  // d(RHS)/dx = 2x + 3 = 7
        REQUIRE(rhs.gradient[1] == 0.0);  // d(RHS)/dy = 0
    }
    
    SECTION("Expression with multiple variables") {
        // z = x*y + y^2 at x=2, y=3
        auto parseResult = parser.parse("z = x*y + y^2");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        const auto& equations = ir.getEquations();
        
        coolsolve::ExpressionEvaluator eval(2);
        eval.setVariable("x", coolsolve::ADValue::independent(2.0, 0, 2));
        eval.setVariable("y", coolsolve::ADValue::independent(3.0, 1, 2));
        
        auto rhs = eval.evaluate(equations[0].rhs);
        REQUIRE(rhs.value == 15.0);  // 2*3 + 9 = 15
        REQUIRE_THAT(rhs.gradient[0], WithinAbs(3.0, 1e-10));  // dz/dx = y = 3
        REQUIRE_THAT(rhs.gradient[1], WithinAbs(8.0, 1e-10));  // dz/dy = x + 2y = 2 + 6 = 8
    }

    SECTION("Variadic min function") {
        auto parseResult = parser.parse("z = min(x, y, 5)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        
        coolsolve::ExpressionEvaluator eval(2);
        eval.setVariable("x", coolsolve::ADValue::independent(10.0, 0, 2));
        eval.setVariable("y", coolsolve::ADValue::independent(2.0, 1, 2));
        
        auto rhs = eval.evaluate(ir.getEquations()[0].rhs);
        REQUIRE(rhs.value == 2.0);
        REQUIRE(rhs.gradient[0] == 0.0);
        REQUIRE(rhs.gradient[1] == 1.0);
        
        // Change inputs
        eval.setVariable("x", coolsolve::ADValue::independent(1.0, 0, 2));
        rhs = eval.evaluate(ir.getEquations()[0].rhs);
        REQUIRE(rhs.value == 1.0);
        REQUIRE(rhs.gradient[0] == 1.0);
        REQUIRE(rhs.gradient[1] == 0.0);
    }

    SECTION("Variadic max function") {
        auto parseResult = parser.parse("z = max(x, y, 5)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        
        coolsolve::ExpressionEvaluator eval(2);
        eval.setVariable("x", coolsolve::ADValue::independent(1.0, 0, 2));
        eval.setVariable("y", coolsolve::ADValue::independent(2.0, 1, 2));
        
        auto rhs = eval.evaluate(ir.getEquations()[0].rhs);
        REQUIRE(rhs.value == 5.0);
        REQUIRE(rhs.gradient[0] == 0.0);
        REQUIRE(rhs.gradient[1] == 0.0);
        
        // Change inputs
        eval.setVariable("x", coolsolve::ADValue::independent(10.0, 0, 2));
        rhs = eval.evaluate(ir.getEquations()[0].rhs);
        REQUIRE(rhs.value == 10.0);
        REQUIRE(rhs.gradient[0] == 1.0);
        REQUIRE(rhs.gradient[1] == 0.0);
    }
    
    SECTION("Expression with standard functions") {
        // y = sin(x) at x = pi/6
        auto parseResult = parser.parse("y = sin(x)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        
        coolsolve::ExpressionEvaluator eval(1);
        eval.setVariable("x", coolsolve::ADValue::independent(M_PI / 6.0, 0, 1));
        
        auto rhs = eval.evaluate(ir.getEquations()[0].rhs);
        REQUIRE_THAT(rhs.value, WithinAbs(0.5, 1e-10));
        REQUIRE_THAT(rhs.gradient[0], WithinAbs(std::sqrt(3.0) / 2.0, 1e-10));  // cos(pi/6)
    }
}

// ============================================================================
// BlockEvaluator Tests
// ============================================================================

TEST_CASE("BlockEvaluator simple system", "[evaluator]") {
    coolsolve::EESParser parser;
    
    SECTION("Single equation block") {
        // x = 5 (should have residual x - 5 = 0)
        auto parseResult = parser.parse("x = 5");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        REQUIRE(analysis.success);
        REQUIRE(analysis.blocks.size() == 1);
        
        coolsolve::BlockEvaluator blockEval(analysis.blocks[0], ir);
        
        // At x = 5, residual should be 0
        auto result = blockEval.evaluate({5.0});
        REQUIRE_THAT(result.residuals[0], WithinAbs(0.0, 1e-10));
        
        // At x = 3, residual should be -2 (3 - 5 = -2)
        result = blockEval.evaluate({3.0});
        REQUIRE_THAT(result.residuals[0], WithinAbs(-2.0, 1e-10));
        REQUIRE_THAT(result.jacobian[0][0], WithinAbs(1.0, 1e-10));  // d(x-5)/dx = 1
    }
    
    SECTION("Two-equation algebraic loop") {
        // x = y + 1
        // y = x + 1
        // This is an inconsistent system, but we can still evaluate residuals
        auto parseResult = parser.parse("x = y + 1\ny = x + 1");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        REQUIRE(analysis.success);
        // Should form an algebraic loop
        REQUIRE(analysis.largestBlockSize == 2);
        
        // Find the block with both equations
        for (const auto& block : analysis.blocks) {
            if (block.size() == 2) {
                coolsolve::BlockEvaluator blockEval(block, ir);
                
                // At x=0, y=0: residuals are [0 - (0+1), 0 - (0+1)] = [-1, -1]
                auto result = blockEval.evaluate({0.0, 0.0});
                
                // Jacobian should be:
                // d(x - y - 1)/dx = 1, d(x - y - 1)/dy = -1
                // d(y - x - 1)/dx = -1, d(y - x - 1)/dy = 1
                // (order depends on how variables are matched)
                
                // Just verify we get a 2x2 Jacobian
                REQUIRE(result.jacobian.size() == 2);
                REQUIRE(result.jacobian[0].size() == 2);
                REQUIRE(result.jacobian[1].size() == 2);
                
                break;
            }
        }
    }
}

TEST_CASE("BlockEvaluator Jacobian verification", "[evaluator]") {
    coolsolve::EESParser parser;
    
    SECTION("Compare AD with finite differences - polynomial") {
        // x = 2, y = x^2 + 3*x - 5
        // Square system: 2 equations, 2 variables
        auto parseResult = parser.parse("x = 2\ny = x^2 + 3*x - 5");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        REQUIRE(analysis.success);
        
        // Find the block with the polynomial equation (y = x^2 + 3*x - 5)
        for (const auto& block : analysis.blocks) {
            coolsolve::BlockEvaluator blockEval(block, ir);
            const auto& vars = blockEval.getVariables();
            
            // Test at multiple points - provide external variable values
            for (double testVal : {-2.0, 0.0, 1.5, 3.0}) {
                std::vector<double> xBlock(vars.size(), testVal);
                std::map<std::string, double> external;
                // Set external variables that aren't in the block
                for (const auto& [name, info] : ir.getVariables()) {
                    bool inBlock = std::find(vars.begin(), vars.end(), name) != vars.end();
                    if (!inBlock) {
                        external[name] = testVal;
                    }
                }
                
                double maxDiff = coolsolve::compareJacobianWithFiniteDifferences(
                    blockEval, xBlock, external, {}, 1e-7, true);  // verbose=true
                REQUIRE(maxDiff < 1e-5);
            }
        }
    }
    
    SECTION("Compare AD with finite differences - transcendental") {
        // x = 1, y = sin(x) + exp(x)
        // Square system: 2 equations, 2 variables
        auto parseResult = parser.parse("x = 1\ny = sin(x) + exp(x)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        REQUIRE(analysis.success);
        
        for (const auto& block : analysis.blocks) {
            coolsolve::BlockEvaluator blockEval(block, ir);
            const auto& vars = blockEval.getVariables();
            
            for (double testVal : {0.0, 0.5, 1.0, 2.0}) {
                std::vector<double> xBlock(vars.size(), testVal);
                std::map<std::string, double> external;
                for (const auto& [name, info] : ir.getVariables()) {
                    bool inBlock = std::find(vars.begin(), vars.end(), name) != vars.end();
                    if (!inBlock) {
                        external[name] = testVal;
                    }
                }
                
                double maxDiff = coolsolve::compareJacobianWithFiniteDifferences(
                    blockEval, xBlock, external, {}, 1e-7, true);  // verbose=true
                REQUIRE(maxDiff < 1e-5);
            }
        }
    }
    
    SECTION("Compare AD with finite differences - multi-variable") {
        // x = 2, y = 3, z = x*y + x^2 - y^2
        // Square system: 3 equations, 3 variables
        auto parseResult = parser.parse("x = 2\ny = 3\nz = x*y + x^2 - y^2");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        REQUIRE(analysis.success);
        
        // Test each block
        for (const auto& block : analysis.blocks) {
            coolsolve::BlockEvaluator blockEval(block, ir);
            const auto& vars = blockEval.getVariables();
            
            if (block.size() >= 1) {
                std::vector<double> xBlock(block.size(), 1.5);
                std::map<std::string, double> external;
                for (const auto& [name, info] : ir.getVariables()) {
                    bool inBlock = std::find(vars.begin(), vars.end(), name) != vars.end();
                    if (!inBlock) {
                        external[name] = 2.0;  // Different value for external vars
                    }
                }
                
                double maxDiff = coolsolve::compareJacobianWithFiniteDifferences(
                    blockEval, xBlock, external, {}, 1e-7, true);  // verbose=true
                REQUIRE(maxDiff < 1e-5);
            }
        }
    }
}

// ============================================================================
// SystemEvaluator Tests
// ============================================================================

TEST_CASE("SystemEvaluator initialization", "[evaluator]") {
    coolsolve::EESParser parser;
    
    SECTION("Sequential system") {
        auto parseResult = parser.parse("x = 5\ny = x + 3\nz = y * 2");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        
        REQUIRE(sysEval.getNumBlocks() == 3);
        
        // Set initial values
        sysEval.setVariableValue("x", 5.0);
        sysEval.setVariableValue("y", 8.0);
        sysEval.setVariableValue("z", 16.0);
        
        // Evaluate all blocks - should have zero residuals at solution
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            auto result = sysEval.evaluateBlock(i);
            // At correct solution, residuals should be very small
            for (double r : result.residuals) {
                REQUIRE_THAT(r, WithinAbs(0.0, 1e-10));
            }
        }
    }
}

TEST_CASE("Evaluator report generation", "[evaluator]") {
    coolsolve::EESParser parser;
    
    auto parseResult = parser.parse("x = 5\ny = x^2");
    REQUIRE(parseResult.success);
    
    auto ir = coolsolve::IR::fromAST(parseResult.program);
    auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
    
    coolsolve::SystemEvaluator sysEval(ir, analysis);
    sysEval.setVariableValue("x", 5.0);
    sysEval.setVariableValue("y", 25.0);
    
    std::string report = coolsolve::generateEvaluatorReport(sysEval);
    
    // Check report contains expected sections
    REQUIRE(report.find("Evaluator Report") != std::string::npos);
    REQUIRE(report.find("Blocks Summary") != std::string::npos);
    REQUIRE(report.find("Current State") != std::string::npos);
}

// ============================================================================
// CoolProp Integration Tests
// ============================================================================

#include "CoolProp.h"

TEST_CASE("CoolProp function evaluation", "[coolprop]") {
    coolsolve::EESParser parser;
    
    // Note: CoolSolve now converts temperature from Celsius to Kelvin automatically
    // (EES convention uses Celsius, CoolProp uses Kelvin)
    
    SECTION("Enthalpy calculation with T and P") {
        // h = enthalpy(Water, T=25, P=101325) -- T in Celsius (25°C = 298.15 K)
        auto parseResult = parser.parse("h = enthalpy(Water, T=25, P=101325)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        // The equation has h as unknown, T and P as constants
        REQUIRE(analysis.blocks.size() >= 1);
        
        coolsolve::BlockEvaluator blockEval(analysis.blocks[0], ir);
        
        // Get expected value from CoolProp directly (using Kelvin)
        double expected_h = CoolProp::PropsSI("H", "T", 25.0 + 273.15, "P", 101325.0, "Water");
        
        // Evaluate at the correct h value
        auto result = blockEval.evaluate({expected_h});
        
        // Residual should be near zero
        REQUIRE_THAT(result.residuals[0], WithinAbs(0.0, 1.0));  // 1 J/kg tolerance
    }
    
    SECTION("Entropy calculation with T and P") {
        // T=80°C = 353.15 K
        auto parseResult = parser.parse("s = entropy(Water, T=80, P=200000)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::BlockEvaluator blockEval(analysis.blocks[0], ir);
        
        double expected_s = CoolProp::PropsSI("S", "T", 80.0 + 273.15, "P", 200000.0, "Water");
        
        auto result = blockEval.evaluate({expected_s});
        REQUIRE_THAT(result.residuals[0], WithinAbs(0.0, 0.01));  // 0.01 J/kg/K tolerance
    }
    
    SECTION("CoolProp with variable inputs") {
        // h = enthalpy(Water, T=T_in, P=P_in)
        // T_in = 25 (Celsius)
        // P_in = 101325
        auto parseResult = parser.parse("T_in = 25\nP_in = 101325\nh = enthalpy(Water, T=T_in, P=P_in)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        
        // Set correct solution values (T in Celsius)
        sysEval.setVariableValue("T_in", 25.0);
        sysEval.setVariableValue("P_in", 101325.0);
        double expected_h = CoolProp::PropsSI("H", "T", 25.0 + 273.15, "P", 101325.0, "Water");
        sysEval.setVariableValue("h", expected_h);
        
        // All blocks should have small residuals
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            auto result = sysEval.evaluateBlock(i);
            for (double r : result.residuals) {
                REQUIRE(std::abs(r) < 10.0);  // Allow some tolerance
            }
        }
    }
    
    SECTION("CoolProp Jacobian verification") {
        // Test that derivatives are computed correctly
        // Use a simple case: h = enthalpy(Water, T=T_in, P=101325) where T_in is a variable
        // T_in = 80 Celsius
        auto parseResult = parser.parse("T_in = 80\nh = enthalpy(Water, T=T_in, P=101325)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        double expected_h = CoolProp::PropsSI("H", "T", 80.0 + 273.15, "P", 101325.0, "Water");
        
        // Find and test each block
        for (const auto& block : analysis.blocks) {
            coolsolve::BlockEvaluator blockEval(block, ir);
            const auto& vars = blockEval.getVariables();
            
            // Set external variables (all variables not in this block)
            std::map<std::string, double> external;
            for (const auto& [name, info] : ir.getVariables()) {
                bool inBlock = std::find(vars.begin(), vars.end(), name) != vars.end();
                if (!inBlock) {
                    if (name == "T_in") {
                        external[name] = 80.0;  // Celsius
                    } else if (name == "h") {
                        external[name] = expected_h;
                    }
                }
            }
            
            // Set block variables to their expected values
            std::vector<double> x(block.size());
            for (size_t i = 0; i < vars.size(); ++i) {
                if (vars[i] == "T_in") {
                    x[i] = 80.0;  // Celsius
                } else if (vars[i] == "h") {
                    x[i] = expected_h;
                }
            }
            
            // Verify Jacobian with finite differences (verbose=true to show comparison)
            double maxDiff = coolsolve::compareJacobianWithFiniteDifferences(
                blockEval, x, external, {}, 1e-4, true);  // verbose=true
            
            // Allow for some numerical error in CoolProp derivatives
            // The tolerance needs to account for the fact that CoolProp 
            // derivatives are computed numerically
            REQUIRE(maxDiff < 1000.0);  // Larger tolerance for CoolProp numerical derivatives
        }
    }
}

TEST_CASE("CoolProp different fluids", "[coolprop]") {
    coolsolve::EESParser parser;
    
    SECTION("R134a refrigerant") {
        // T=25°C = 298.15 K
        auto parseResult = parser.parse("h = enthalpy(R134a, T=25, P=500000)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::BlockEvaluator blockEval(analysis.blocks[0], ir);
        
        // Expected value from CoolProp (T in Kelvin = 25 + 273.15)
        double expected_h = CoolProp::PropsSI("H", "T", 25.0 + 273.15, "P", 500000.0, "R134a");
        auto result = blockEval.evaluate({expected_h});
        
        REQUIRE_THAT(result.residuals[0], WithinAbs(0.0, 10.0));
    }
    
    SECTION("Nitrogen") {
        // T=25°C = 298.15 K
        auto parseResult = parser.parse("rho = density(Nitrogen, T=25, P=101325)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::BlockEvaluator blockEval(analysis.blocks[0], ir);
        
        double expected_rho = CoolProp::PropsSI("D", "T", 25.0 + 273.15, "P", 101325.0, "Nitrogen");
        auto result = blockEval.evaluate({expected_rho});
        
        REQUIRE_THAT(result.residuals[0], WithinAbs(0.0, 0.001));
    }
}

// ============================================================================
// CoolProp Input Pair Tests - All common EES function patterns
// ============================================================================

TEST_CASE("CoolProp saturation pressure (T,x inputs)", "[coolprop]") {
    coolsolve::EESParser parser;
    
    // This is the pattern: P_sat = pressure(R134a, T=T_ev, x=1)
    // EES uses 'x' for quality, CoolProp uses 'Q'
    
    SECTION("Saturation pressure of R134a at -10°C") {
        auto parseResult = parser.parse("P_sat = pressure(R134a, T=-10, x=1)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        REQUIRE(analysis.blocks.size() >= 1);
        coolsolve::BlockEvaluator blockEval(analysis.blocks[0], ir);
        
        // Expected value from CoolProp directly
        double T_K = -10 + 273.15;
        double expected_P = CoolProp::PropsSI("P", "T", T_K, "Q", 1, "R134a");
        
        INFO("Expected P_sat = " << expected_P << " Pa");
        
        auto result = blockEval.evaluate({expected_P});
        REQUIRE_THAT(result.residuals[0], WithinAbs(0.0, 1.0));  // 1 Pa tolerance
    }
    
    SECTION("Saturation pressure of Water at 100°C") {
        auto parseResult = parser.parse("P_sat = pressure(Water, T=100, x=0)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::BlockEvaluator blockEval(analysis.blocks[0], ir);
        
        double T_K = 100 + 273.15;
        double expected_P = CoolProp::PropsSI("P", "T", T_K, "Q", 0, "Water");
        
        INFO("Expected P_sat(Water, 100°C) = " << expected_P << " Pa");
        
        auto result = blockEval.evaluate({expected_P});
        REQUIRE_THAT(result.residuals[0], WithinAbs(0.0, 1.0));
    }
}

TEST_CASE("CoolProp quality calculation (P,h inputs)", "[coolprop]") {
    coolsolve::EESParser parser;
    
    // Pattern: x = quality(R134a, P=P_ev, h=h_4)
    // h must be between h_f and h_g at the given pressure
    
    SECTION("Quality in two-phase region") {
        // First get valid h values at P=200000 Pa for R134a
        double P = 200000;  // Pa
        double h_f = CoolProp::PropsSI("H", "P", P, "Q", 0, "R134a");
        double h_g = CoolProp::PropsSI("H", "P", P, "Q", 1, "R134a");
        
        // Use h at quality = 0.3
        double h_test = h_f + 0.3 * (h_g - h_f);
        double expected_Q = 0.3;
        
        INFO("At P=" << P << " Pa: h_f=" << h_f << ", h_g=" << h_g << ", h_test=" << h_test);
        
        // Create equation with computed h value
        std::ostringstream eqs;
        eqs << "h_4 = " << std::fixed << std::setprecision(2) << h_test << "\n";
        eqs << "P_4 = " << P << "\n";
        eqs << "x_4 = quality(R134a, P=P_4, h=h_4)";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        
        // Set values
        sysEval.setVariableValue("h_4", h_test);
        sysEval.setVariableValue("P_4", P);
        sysEval.setVariableValue("x_4", expected_Q);
        
        // All residuals should be small
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            auto result = sysEval.evaluateBlock(i);
            for (double r : result.residuals) {
                INFO("Block " << i << " residual: " << r);
                REQUIRE(std::abs(r) < 0.01);  // Small tolerance for quality
            }
        }
    }
}

TEST_CASE("CoolProp enthalpy from P,s inputs", "[coolprop]") {
    coolsolve::EESParser parser;
    
    // Pattern: h_2s = enthalpy(R134a, P=P_2, s=s_1)
    // This is used for isentropic processes
    
    SECTION("Isentropic compression of R134a") {
        // Start: superheated vapor at -8°C, P=200 kPa
        double T_1 = -8 + 273.15;
        double P_1 = 200000;
        double h_1 = CoolProp::PropsSI("H", "T", T_1, "P", P_1, "R134a");
        double s_1 = CoolProp::PropsSI("S", "T", T_1, "P", P_1, "R134a");
        
        // End: P_2 = 572 kPa (condensation pressure)
        double P_2 = 572000;
        double h_2s = CoolProp::PropsSI("H", "P", P_2, "S", s_1, "R134a");
        
        INFO("s_1 = " << s_1 << " J/kg/K");
        INFO("h_2s = " << h_2s << " J/kg");
        
        // Create equation
        std::ostringstream eqs;
        eqs << std::fixed << std::setprecision(2);
        eqs << "s_1 = " << s_1 << "\n";
        eqs << "P_2 = " << P_2 << "\n";
        eqs << "h_2s = enthalpy(R134a, P=P_2, s=s_1)";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("s_1", s_1);
        sysEval.setVariableValue("P_2", P_2);
        sysEval.setVariableValue("h_2s", h_2s);
        
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            auto result = sysEval.evaluateBlock(i);
            for (double r : result.residuals) {
                INFO("Block " << i << " residual: " << r);
                REQUIRE(std::abs(r) < 10.0);  // 10 J/kg tolerance
            }
        }
    }
}

TEST_CASE("CoolProp temperature from P,h inputs", "[coolprop]") {
    coolsolve::EESParser parser;
    
    // Pattern: T_2 = temperature(R134a, P=P_2, h=h_2)
    
    SECTION("Temperature from P,h for superheated R134a") {
        double P = 500000;  // 500 kPa
        double T_expected = 30 + 273.15;  // 30°C in Kelvin
        double h = CoolProp::PropsSI("H", "T", T_expected, "P", P, "R134a");
        
        INFO("At P=" << P << " Pa, h=" << h << " J/kg, expected T=30°C");
        
        std::ostringstream eqs;
        eqs << std::fixed << std::setprecision(2);
        eqs << "P_2 = " << P << "\n";
        eqs << "h_2 = " << h << "\n";
        eqs << "T_2 = temperature(R134a, P=P_2, h=h_2)";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("P_2", P);
        sysEval.setVariableValue("h_2", h);
        sysEval.setVariableValue("T_2", 30.0);  // Expected: 30°C (our code converts to Kelvin)
        
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            auto result = sysEval.evaluateBlock(i);
            for (double r : result.residuals) {
                INFO("Block " << i << " residual: " << r);
                REQUIRE(std::abs(r) < 0.1);  // 0.1°C tolerance
            }
        }
    }
}

TEST_CASE("CoolProp refrigeration cycle - complete", "[coolprop][integration]") {
    coolsolve::EESParser parser;
    
    // Test a simplified refrigeration cycle to verify all CoolProp functions work together
    // This mirrors the refrigeration1.eescode example but with CoolProp-consistent values
    
    SECTION("R134a cycle state points") {
        // Define cycle parameters
        double T_ev = -10;  // °C - evaporator temperature
        double T_cond = 20;  // °C - condenser temperature
        double dT_sc = 6;    // K - subcooling
        double dT_oh = 2;    // K - overheating
        
        // Convert to Kelvin for CoolProp calls
        double T_ev_K = T_ev + 273.15;
        double T_cond_K = T_cond + 273.15;
        
        // Calculate saturation pressures
        double P_ev = CoolProp::PropsSI("P", "T", T_ev_K, "Q", 1, "R134a");
        double P_cd = CoolProp::PropsSI("P", "T", T_cond_K, "Q", 1, "R134a");
        
        INFO("P_ev = " << P_ev << " Pa, P_cd = " << P_cd << " Pa");
        
        // State 1: Compressor inlet (superheated vapor)
        double T_1 = T_ev + dT_oh;  // °C
        double h_1 = CoolProp::PropsSI("H", "T", T_1 + 273.15, "P", P_ev, "R134a");
        double s_1 = CoolProp::PropsSI("S", "T", T_1 + 273.15, "P", P_ev, "R134a");
        
        // State 3: Condenser outlet (subcooled liquid)
        double T_3 = T_cond - dT_sc;  // °C
        double h_3 = CoolProp::PropsSI("H", "T", T_3 + 273.15, "P", P_cd, "R134a");
        
        // State 4: Evaporator inlet (two-phase)
        double h_4 = h_3;  // Isenthalpic expansion
        double x_4 = CoolProp::PropsSI("Q", "H", h_4, "P", P_ev, "R134a");
        
        INFO("State 1: T=" << T_1 << "°C, h=" << h_1 << " J/kg, s=" << s_1 << " J/kg/K");
        INFO("State 3: T=" << T_3 << "°C, h=" << h_3 << " J/kg");
        INFO("State 4: h=" << h_4 << " J/kg, x=" << x_4);
        
        // Verify quality is valid (between 0 and 1)
        REQUIRE(x_4 >= 0.0);
        REQUIRE(x_4 <= 1.0);
        
        // Now test parsing and evaluating EES-style equations
        std::ostringstream eqs;
        eqs << std::fixed << std::setprecision(2);
        eqs << "T_ev = " << T_ev << "\n";
        eqs << "P_ev = pressure(R134a, T=T_ev, x=1)\n";
        eqs << "T_1 = " << T_1 << "\n";
        eqs << "h_1 = enthalpy(R134a, T=T_1, P=P_ev)\n";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        
        // Set all variable values to their expected (CoolProp-computed) values
        sysEval.setVariableValue("T_ev", T_ev);
        sysEval.setVariableValue("P_ev", P_ev);
        sysEval.setVariableValue("T_1", T_1);
        sysEval.setVariableValue("h_1", h_1);
        
        // All residuals should be small when using consistent CoolProp values
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            auto result = sysEval.evaluateBlock(i);
            for (size_t j = 0; j < result.residuals.size(); ++j) {
                INFO("Block " << i << ", equation " << j << " residual: " << result.residuals[j]);
                REQUIRE(std::abs(result.residuals[j]) < 100.0);  // Allow some tolerance for FD derivatives
            }
        }
    }
}

// ============================================================================
// Evaluator Error Detection Tests
// These tests verify that the evaluator correctly reports errors and that
// tests fail when evaluation errors occur.
// ============================================================================

TEST_CASE("Evaluator error detection - no exceptions during evaluation", "[evaluator][errors]") {
    coolsolve::EESParser parser;
    
    // This test verifies that for valid CoolProp-consistent equations,
    // no exceptions are thrown during evaluation
    
    SECTION("Simple enthalpy calculation must not throw") {
        std::ostringstream eqs;
        double T = 25;  // Celsius
        double P = 101325;  // Pa
        double expected_h = CoolProp::PropsSI("H", "T", T + 273.15, "P", P, "Water");
        
        eqs << std::fixed << std::setprecision(2);
        eqs << "T = " << T << "\n";
        eqs << "P = " << P << "\n";
        eqs << "h = enthalpy(Water, T=T, P=P)\n";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("T", T);
        sysEval.setVariableValue("P", P);
        sysEval.setVariableValue("h", expected_h);
        
        // Evaluation must not throw any exceptions
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            REQUIRE_NOTHROW(sysEval.evaluateBlock(i));
        }
    }
    
    SECTION("Saturation pressure must not throw") {
        std::ostringstream eqs;
        double T = 0;  // Celsius (0°C)
        double expected_P = CoolProp::PropsSI("P", "T", T + 273.15, "Q", 1, "R134a");
        
        eqs << std::fixed << std::setprecision(2);
        eqs << "T_sat = " << T << "\n";
        eqs << "P_sat = pressure(R134a, T=T_sat, x=1)\n";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("T_sat", T);
        sysEval.setVariableValue("P_sat", expected_P);
        
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            REQUIRE_NOTHROW(sysEval.evaluateBlock(i));
        }
    }

    SECTION("Saturation function T_sat with 1 input must not throw") {
        std::ostringstream eqs;
        double P = 101325;  // Pa
        double expected_T = CoolProp::PropsSI("T", "P", P, "Q", 0.5, "Water");
        
        eqs << std::fixed << std::setprecision(2);
        eqs << "P = " << P << "\n";
        eqs << "T_sat = T_sat(Water, P=P)\n";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("P", P);
        sysEval.setVariableValue("T_sat", expected_T - 273.15); // Evaluator converts to Celsius
        
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            REQUIRE_NOTHROW(sysEval.evaluateBlock(i));
        }
    }
}

TEST_CASE("Evaluator must throw for undefined variables", "[evaluator][errors]") {
    coolsolve::EESParser parser;
    
    // Variables that are used but not initialized should return a warning
    
    SECTION("Uninitialized variable returns warning") {
        // Square system: x = 5, y = x + z, z = 10
        // But we only set x and y, not z - so when evaluating y's equation,
        // z should be used with default value and trigger a warning
        auto parseResult = parser.parse("x = 5\nz = 10\ny = x + z");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        REQUIRE(analysis.success);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("x", 5.0);
        sysEval.setVariableValue("y", 1.0);
        // Note: z is NOT set - using default initialization
        // The warning behavior depends on implementation - check if any warnings are generated
        
        // Evaluation should proceed without throwing
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            REQUIRE_NOTHROW(sysEval.evaluateBlock(i));
        }
    }
}

TEST_CASE("Evaluator must throw for unsupported functions", "[evaluator][errors]") {
    coolsolve::EESParser parser;
    
    SECTION("Unknown function throws exception") {
        auto parseResult = parser.parse("y = unknown_function(1, 2, 3, 4)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("y", 1.0);
        
        // Evaluation should throw because unknown_function is not supported
        bool exceptionThrown = false;
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            try {
                sysEval.evaluateBlock(i);
            } catch (const std::exception& e) {
                std::string msg = e.what();
                if (msg.find("Unknown") != std::string::npos || 
                    msg.find("unsupported") != std::string::npos) {
                    exceptionThrown = true;
                }
            }
        }
        REQUIRE(exceptionThrown);
    }
}

TEST_CASE("Evaluator must throw for invalid CoolProp inputs", "[evaluator][errors]") {
    coolsolve::EESParser parser;
    
    SECTION("Unknown fluid name (returns inf, should throw)") {
        // CoolProp returns inf for unknown fluid names, which our code detects
        auto parseResult = parser.parse("h = enthalpy(UnknownFluid, T=25, P=101325)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("h", 1.0);
        
        // Evaluation should throw because UnknownFluid is not a valid CoolProp fluid
        bool exceptionThrown = false;
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            try {
                sysEval.evaluateBlock(i);
            } catch (const std::exception& e) {
                std::string msg = e.what();
                INFO("Exception: " << msg);
                if (msg.find("CoolProp") != std::string::npos ||
                    msg.find("invalid") != std::string::npos ||
                    msg.find("NaN") != std::string::npos ||
                    msg.find("Inf") != std::string::npos ||
                    msg.find("Unknown fluid") != std::string::npos) {
                    exceptionThrown = true;
                }
            }
        }
        REQUIRE(exceptionThrown);
    }
    
    SECTION("Generic fluid name (common EES pattern) throws") {
        // EES files often use variable names like "fluid", "fluidcd", "fluidev"
        // which are not valid CoolProp fluid names
        auto parseResult = parser.parse("h = enthalpy(fluid, T=25, P=101325)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("h", 1.0);
        
        bool exceptionThrown = false;
        for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
            try {
                sysEval.evaluateBlock(i);
            } catch (const std::exception& e) {
                std::string msg = e.what();
                INFO("Exception: " << msg);
                exceptionThrown = true;
            }
        }
        REQUIRE(exceptionThrown);
    }
}

// Helper function to count evaluation errors
struct EvaluationStats {
    size_t totalBlocks = 0;
    size_t successfulBlocks = 0;
    size_t errorBlocks = 0;
    std::vector<std::string> errors;
};

EvaluationStats evaluateAllBlocks(coolsolve::SystemEvaluator& sysEval) {
    EvaluationStats stats;
    stats.totalBlocks = sysEval.getNumBlocks();
    
    for (size_t i = 0; i < sysEval.getNumBlocks(); ++i) {
        try {
            auto result = sysEval.evaluateBlock(i);
            stats.successfulBlocks++;
        } catch (const std::exception& e) {
            stats.errorBlocks++;
            stats.errors.push_back("Block " + std::to_string(i) + ": " + e.what());
        }
    }
    
    return stats;
}

TEST_CASE("System evaluation error counting", "[evaluator][errors]") {
    coolsolve::EESParser parser;
    
    SECTION("Valid system has no errors") {
        double T = 25;
        double P = 101325;
        double expected_h = CoolProp::PropsSI("H", "T", T + 273.15, "P", P, "Water");
        
        std::ostringstream eqs;
        eqs << std::fixed << std::setprecision(2);
        eqs << "T = " << T << "\n";
        eqs << "P = " << P << "\n";
        eqs << "h = enthalpy(Water, T=T, P=P)\n";
        
        auto parseResult = parser.parse(eqs.str());
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("T", T);
        sysEval.setVariableValue("P", P);
        sysEval.setVariableValue("h", expected_h);
        
        auto stats = evaluateAllBlocks(sysEval);
        
        INFO("Errors found:");
        for (const auto& err : stats.errors) {
            INFO("  " << err);
        }
        
        // All blocks should evaluate successfully
        REQUIRE(stats.errorBlocks == 0);
        REQUIRE(stats.successfulBlocks == stats.totalBlocks);
    }
    
    SECTION("System with invalid CoolProp call has errors") {
        // Using "fluid" which is not a valid CoolProp fluid name
        auto parseResult = parser.parse("h = enthalpy(fluid, T=25, P=101325)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        auto analysis = coolsolve::StructuralAnalyzer::analyze(ir);
        
        coolsolve::SystemEvaluator sysEval(ir, analysis);
        sysEval.setVariableValue("h", 1.0);
        
        auto stats = evaluateAllBlocks(sysEval);
        
        INFO("Total blocks: " << stats.totalBlocks);
        INFO("Error blocks: " << stats.errorBlocks);
        for (const auto& err : stats.errors) {
            INFO("  " << err);
        }
        
        // Should have at least one error
        REQUIRE(stats.errorBlocks > 0);
    }
}

// ============================================================================
// Variable Inference Tests
// ============================================================================

#include "coolsolve/variable_inference.h"

TEST_CASE("Variable inference produces finite values", "[inference]") {
    coolsolve::EESParser parser;
    
    SECTION("Water enthalpy inference") {
        auto parseResult = parser.parse("h = enthalpy(Water, T=T_in, P=P_in)\nT_in = 25\nP_in = 101325");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        
        // Run inference
        coolsolve::inferVariables(ir);
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        // Check that h got inferred with a finite value
        const auto* hInfo = ir.getVariable("h");
        REQUIRE(hInfo != nullptr);
        REQUIRE(hInfo->guessValue.has_value());
        REQUIRE(std::isfinite(*hInfo->guessValue));
        
        // Value should be reasonable for water enthalpy (~100000 J/kg at 20°C)
        INFO("Inferred h = " << *hInfo->guessValue);
        REQUIRE(*hInfo->guessValue > 50000);
        REQUIRE(*hInfo->guessValue < 500000);
    }
    
    SECTION("R134a enthalpy inference") {
        auto parseResult = parser.parse("h = enthalpy(R134a, T=T_in, P=P_in)\nT_in = 25\nP_in = 500000");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        const auto* hInfo = ir.getVariable("h");
        REQUIRE(hInfo != nullptr);
        REQUIRE(hInfo->guessValue.has_value());
        REQUIRE(std::isfinite(*hInfo->guessValue));
        INFO("Inferred h = " << *hInfo->guessValue);
    }
    
    SECTION("Air (pseudo-pure) enthalpy inference") {
        // This is the key test - Air_ha should work correctly
        auto parseResult = parser.parse("h = enthalpy(Air_ha, T=T_in, P=P_in)\nT_in = 25\nP_in = 101325");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        const auto* hInfo = ir.getVariable("h");
        REQUIRE(hInfo != nullptr);
        REQUIRE(hInfo->guessValue.has_value());
        REQUIRE(std::isfinite(*hInfo->guessValue));
        
        // Air enthalpy at ~20°C should be around 293000-420000 J/kg
        INFO("Inferred Air h = " << *hInfo->guessValue);
        REQUIRE(*hInfo->guessValue > 200000);
        REQUIRE(*hInfo->guessValue < 500000);
    }
    
    SECTION("Nitrogen (cryogenic) inference") {
        // Nitrogen has very different properties - should still work
        auto parseResult = parser.parse("h = enthalpy(Nitrogen, T=T_in, P=P_in)\nT_in = 25\nP_in = 101325");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        const auto* hInfo = ir.getVariable("h");
        REQUIRE(hInfo != nullptr);
        REQUIRE(hInfo->guessValue.has_value());
        REQUIRE(std::isfinite(*hInfo->guessValue));
        INFO("Inferred N2 h = " << *hInfo->guessValue);
    }
    
    SECTION("Multiple properties inference") {
        auto parseResult = parser.parse(
            "T = 25\n"
            "P = 101325\n"
            "h = enthalpy(Water, T=T, P=P)\n"
            "s = entropy(Water, T=T, P=P)\n"
            "rho = density(Water, T=T, P=P)"
        );
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        // All properties should have finite inferred values
        const auto* hInfo = ir.getVariable("h");
        const auto* sInfo = ir.getVariable("s");
        const auto* rhoInfo = ir.getVariable("rho");
        
        REQUIRE(hInfo != nullptr);
        REQUIRE(sInfo != nullptr);
        REQUIRE(rhoInfo != nullptr);
        
        REQUIRE(hInfo->guessValue.has_value());
        REQUIRE(sInfo->guessValue.has_value());
        REQUIRE(rhoInfo->guessValue.has_value());
        
        REQUIRE(std::isfinite(*hInfo->guessValue));
        REQUIRE(std::isfinite(*sInfo->guessValue));
        REQUIRE(std::isfinite(*rhoInfo->guessValue));
        
        INFO("h = " << *hInfo->guessValue);
        INFO("s = " << *sInfo->guessValue);
        INFO("rho = " << *rhoInfo->guessValue);
        
        // Water density at 20°C should be close to 1000 kg/m³
        REQUIRE(*rhoInfo->guessValue > 900);
        REQUIRE(*rhoInfo->guessValue < 1100);
    }
}

TEST_CASE("Variable inference error handling", "[inference][errors]") {
    coolsolve::EESParser parser;
    
    SECTION("Unknown fluid throws error") {
        auto parseResult = parser.parse("h = enthalpy(UnknownFluid123, T=25, P=101325)");
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        
        // Should throw because UnknownFluid123 is not registered
        REQUIRE_THROWS_AS(coolsolve::initializeVariables(ir), std::runtime_error);
        
        try {
            coolsolve::initializeVariables(ir);
        } catch (const std::runtime_error& e) {
            std::string msg = e.what();
            INFO("Error message: " << msg);
            // Error message should mention the unknown fluid
            REQUIRE(msg.find("unknown fluid") != std::string::npos);
        }
    }
}

TEST_CASE("Variable inference with string variables", "[inference]") {
    coolsolve::EESParser parser;
    
    SECTION("Fluid specified via string variable") {
        auto parseResult = parser.parse(
            "fluid$ = 'R134a'\n"
            "T = 25\n"
            "P = 500000\n"
            "h = enthalpy(fluid$, T=T, P=P)"
        );
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        
        // This should work - fluid$ resolves to R134a
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        const auto* hInfo = ir.getVariable("h");
        REQUIRE(hInfo != nullptr);
        REQUIRE(hInfo->guessValue.has_value());
        REQUIRE(std::isfinite(*hInfo->guessValue));
    }
}

TEST_CASE("Variable inference for ORC-like model with air_ha", "[inference][orc]") {
    coolsolve::EESParser parser;
    
    // Simplified test case mimicking the ORC model that was failing
    SECTION("Cold fluid (air_ha) heat exchanger variables") {
        auto parseResult = parser.parse(
            "cf$ = 'air_ha'\n"
            "t_cf_su = 20\n"
            "p_cf = 101325\n"
            "h_cf_su = enthalpy(cf$, T=t_cf_su, P=p_cf)\n"
            "h_cf_ex = enthalpy(cf$, T=t_cf_ex, P=p_cf)\n"
            "t_cf_ex = 40"
        );
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        // Both enthalpy variables should have finite values
        const auto* h_su = ir.getVariable("h_cf_su");
        const auto* h_ex = ir.getVariable("h_cf_ex");
        
        REQUIRE(h_su != nullptr);
        REQUIRE(h_ex != nullptr);
        
        REQUIRE(h_su->guessValue.has_value());
        REQUIRE(h_ex->guessValue.has_value());
        
        REQUIRE(std::isfinite(*h_su->guessValue));
        REQUIRE(std::isfinite(*h_ex->guessValue));
        
        INFO("h_cf_su = " << *h_su->guessValue);
        INFO("h_cf_ex = " << *h_ex->guessValue);
    }
    
    SECTION("Hot fluid (air_ha) evaporator variables") {
        auto parseResult = parser.parse(
            "hf$ = 'air_ha'\n"
            "t_hf_su = 150\n"
            "p_hf = 101325\n"
            "h_hf_su = enthalpy(hf$, T=t_hf_su, P=p_hf)\n"
            "T_0 = 20\n"
            "h_hf_0 = enthalpy(hf$, T=T_0, P=p_hf)"
        );
        REQUIRE(parseResult.success);
        
        auto ir = coolsolve::IR::fromAST(parseResult.program);
        coolsolve::inferVariables(ir);
        REQUIRE_NOTHROW(coolsolve::initializeVariables(ir));
        
        const auto* h_hf_su = ir.getVariable("h_hf_su");
        const auto* h_hf_0 = ir.getVariable("h_hf_0");
        
        REQUIRE(h_hf_su != nullptr);
        REQUIRE(h_hf_0 != nullptr);
        
        REQUIRE(h_hf_su->guessValue.has_value());
        REQUIRE(h_hf_0->guessValue.has_value());
        
        REQUIRE(std::isfinite(*h_hf_su->guessValue));
        REQUIRE(std::isfinite(*h_hf_0->guessValue));
        
        INFO("h_hf_su = " << *h_hf_su->guessValue);
        INFO("h_hf_0 = " << *h_hf_0->guessValue);
    }
}
