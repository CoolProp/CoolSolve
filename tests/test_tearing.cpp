/**
 * Tests for structural tearing (feedback vertex set) and solver option enableTearing.
 */

#include <catch2/catch_test_macros.hpp>
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/solver.h"

using namespace coolsolve;

// ============================================================================
// Tear set computation (computeBlockTearSet)
// ============================================================================

TEST_CASE("computeBlockTearSet returns empty for size-1 block", "[tearing][structural]") {
    std::string code = R"(x = 1)";
    EESParser parser;
    auto parseResult = parser.parse(code);
    REQUIRE(parseResult.success);
    IR ir = IR::fromAST(parseResult.program);
    StructuralAnalysisResult analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);
    REQUIRE(analysis.blocks.size() >= 1);
    const Block& block = analysis.blocks[0];
    REQUIRE(block.size() == 1);

    BlockTearSetResult result = computeBlockTearSet(block, ir);
    REQUIRE(result.tearVarNames.empty());
    REQUIRE(result.topoOrderNonTearEqIds.empty());
    REQUIRE(result.tearEquationIds.empty());
}

TEST_CASE("computeBlockTearSet on 2-equation cycle yields one tear", "[tearing][structural]") {
    // x = y + 1, y = x - 1  -> one block of size 2, one tear variable
    std::string code = R"(
        x = y + 1
        y = x - 1
    )";
    EESParser parser;
    auto parseResult = parser.parse(code);
    REQUIRE(parseResult.success);
    IR ir = IR::fromAST(parseResult.program);
    StructuralAnalysisResult analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);
    REQUIRE(analysis.blocks.size() == 1);
    const Block& block = analysis.blocks[0];
    REQUIRE(block.size() == 2);

    BlockTearSetResult result = computeBlockTearSet(block, ir);
    REQUIRE(result.tearVarNames.size() == 1);
    REQUIRE(result.tearEquationIds.size() == 1);
    REQUIRE(result.topoOrderNonTearEqIds.size() == 1);
}

TEST_CASE("computeBlockTearSet on 3-equation cycle yields tear set and acyclic order", "[tearing][structural]") {
    // a = b, b = c, c = a  -> one block of size 3
    std::string code = R"(
        a = b
        b = c
        c = a
    )";
    EESParser parser;
    auto parseResult = parser.parse(code);
    REQUIRE(parseResult.success);
    IR ir = IR::fromAST(parseResult.program);
    StructuralAnalysisResult analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);
    REQUIRE(analysis.blocks.size() == 1);
    const Block& block = analysis.blocks[0];
    REQUIRE(block.size() == 3);

    BlockTearSetResult result = computeBlockTearSet(block, ir);
    REQUIRE_FALSE(result.tearVarNames.empty());
    REQUIRE(result.tearVarNames.size() >= 1);
    REQUIRE(result.tearEquationIds.size() == result.tearVarNames.size());
    REQUIRE(result.topoOrderNonTearEqIds.size() + result.tearEquationIds.size() == 3);
}

// ============================================================================
// Solver with enableTearing
// ============================================================================

TEST_CASE("Solver with enableTearing solves 2-equation cycle", "[tearing][solver][integration]") {
    std::string code = R"(
        x = y + 1
        y = x - 1
    )";
    EESParser parser;
    auto parseResult = parser.parse(code);
    REQUIRE(parseResult.success);
    IR ir = IR::fromAST(parseResult.program);
    StructuralAnalysisResult analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    CoolPropConfig config;
    Solver solver(ir, analysis, config);
    solver.setGuess("x", 1.0);
    solver.setGuess("y", 0.0);

    SolverOptions options;
    options.enableTearing = true;
    options.tearingMinBlockSize = 2;
    options.verbose = false;

    SolveResult result = solver.solve(options, false);
    REQUIRE(result.success);
    REQUIRE(result.variables.count("x") > 0);
    REQUIRE(result.variables.count("y") > 0);
    double x = result.variables.at("x");
    double y = result.variables.at("y");
    REQUIRE(std::abs((x - (y + 1))) < 1e-6);
    REQUIRE(std::abs((y - (x - 1))) < 1e-6);
}

TEST_CASE("Solver without tearing still solves 2-equation cycle", "[tearing][solver][integration]") {
    std::string code = R"(
        x = y + 1
        y = x - 1
    )";
    EESParser parser;
    auto parseResult = parser.parse(code);
    REQUIRE(parseResult.success);
    IR ir = IR::fromAST(parseResult.program);
    StructuralAnalysisResult analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    CoolPropConfig config;
    Solver solver(ir, analysis, config);
    solver.setGuess("x", 1.0);
    solver.setGuess("y", 0.0);

    SolverOptions options;
    options.enableTearing = false;

    SolveResult result = solver.solve(options, false);
    REQUIRE(result.success);
}

TEST_CASE("Tearing with Schur complement solves strongly-coupled nonlinear system", "[tearing][solver][integration]") {
    // x + y = 10, y = x^2 - 6  ->  strongly coupled nonlinear block.
    // The partial-derivative Jacobian (dF_tear/dx_tear only) would miss the
    // coupling through y, making convergence slow or impossible.
    // The Schur complement captures the total derivative and converges fast.
    // Solutions: x ≈ 3.3166... (positive root of x^2 + x - 16 = 0)
    std::string code = R"(
        x + y = 10
        y = x^2 - 6
    )";
    EESParser parser;
    auto parseResult = parser.parse(code);
    REQUIRE(parseResult.success);
    IR ir = IR::fromAST(parseResult.program);
    StructuralAnalysisResult analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    CoolPropConfig config;
    Solver solver(ir, analysis, config);
    solver.setGuess("x", 5.0);
    solver.setGuess("y", 5.0);

    SolverOptions options;
    options.enableTearing = true;
    options.tearingMinBlockSize = 2;
    options.verbose = false;

    SolveResult result = solver.solve(options, false);
    REQUIRE(result.success);
    double x = result.variables.at("x");
    double y = result.variables.at("y");
    // Check both equations are satisfied
    REQUIRE(std::abs(x + y - 10) < 1e-6);
    REQUIRE(std::abs(y - (x * x - 6)) < 1e-6);
}
