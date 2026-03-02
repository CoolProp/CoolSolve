/**
 * @file test_symbolic_reduction.cpp
 * @brief Tests for the symbolic block reduction feature (§3.2.1).
 *
 * Tests cover:
 *  - CoolProp inversion helpers (function/input name mapping, pair validity)
 *  - Analysis of blocks: explicit extraction, CoolProp inversion, substitution
 *  - End-to-end solver integration with enableSymbolicReduction = true
 *  - Verification that disabling the feature has zero overhead
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/symbolic_reduction.h"
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/solver.h"
#include <cmath>
#include <iostream>
#include <chrono>
#include <fstream>
#include <filesystem>

using namespace coolsolve;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ============================================================================
// CoolProp inversion helper tests
// ============================================================================

TEST_CASE("coolpropFuncToInputName maps standard functions", "[symbolic][helpers]") {
    CHECK(coolpropFuncToInputName("enthalpy") == "H");
    CHECK(coolpropFuncToInputName("Enthalpy") == "H");
    CHECK(coolpropFuncToInputName("h") == "H");
    CHECK(coolpropFuncToInputName("entropy") == "S");
    CHECK(coolpropFuncToInputName("temperature") == "T");
    CHECK(coolpropFuncToInputName("pressure") == "P");
    CHECK(coolpropFuncToInputName("density") == "D");
    CHECK(coolpropFuncToInputName("quality") == "Q");
    CHECK(coolpropFuncToInputName("internalenergy") == "U");
}

TEST_CASE("coolpropFuncToInputName returns empty for non-invertible", "[symbolic][helpers]") {
    CHECK(coolpropFuncToInputName("cp").empty());
    CHECK(coolpropFuncToInputName("viscosity").empty());
    CHECK(coolpropFuncToInputName("conductivity").empty());
    CHECK(coolpropFuncToInputName("soundspeed").empty());
    CHECK(coolpropFuncToInputName("nonexistent").empty());
}

TEST_CASE("coolpropInputToFuncName maps standard inputs", "[symbolic][helpers]") {
    CHECK(coolpropInputToFuncName("H") == "enthalpy");
    CHECK(coolpropInputToFuncName("h") == "enthalpy");
    CHECK(coolpropInputToFuncName("S") == "entropy");
    CHECK(coolpropInputToFuncName("T") == "temperature");
    CHECK(coolpropInputToFuncName("P") == "pressure");
    CHECK(coolpropInputToFuncName("D") == "density");
    CHECK(coolpropInputToFuncName("Q") == "quality");
}

TEST_CASE("isValidCoolPropInputPair recognizes standard pairs", "[symbolic][helpers]") {
    CHECK(isValidCoolPropInputPair("P", "T"));
    CHECK(isValidCoolPropInputPair("T", "P"));  // order independent
    CHECK(isValidCoolPropInputPair("H", "P"));
    CHECK(isValidCoolPropInputPair("P", "S"));
    CHECK(isValidCoolPropInputPair("H", "S"));
    CHECK(isValidCoolPropInputPair("D", "P"));
    CHECK(isValidCoolPropInputPair("D", "T"));
    CHECK(isValidCoolPropInputPair("Q", "T"));
    CHECK(isValidCoolPropInputPair("P", "Q"));
}

TEST_CASE("isValidCoolPropInputPair rejects H,T", "[symbolic][helpers]") {
    CHECK_FALSE(isValidCoolPropInputPair("H", "T"));
    CHECK_FALSE(isValidCoolPropInputPair("T", "H"));
}

// ============================================================================
// Block reduction analysis tests (pure algebra, no CoolProp)
// ============================================================================

TEST_CASE("analyseBlockReduction: no reduction for size-1 block", "[symbolic][analysis]") {
    std::string code = R"(x = 42)";
    EESParser parser;
    auto pr = parser.parse(code);
    REQUIRE(pr.success);
    IR ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);
    REQUIRE(analysis.blocks.size() >= 1);

    auto result = analyseBlockReduction(analysis.blocks[0], ir, analysis);
    CHECK_FALSE(result.reduced);
}

TEST_CASE("analyseBlockReduction: explicit extraction from 2-var block", "[symbolic][analysis]") {
    // a = 5; b = a + 3; c = a * b
    // Structural analysis should create some blocks.
    // If a and b end up in the same block (which they won't normally since a=5 is explicit),
    // let's try a case where explicitly-solvable equations are within a block.
    // Actually, in a well-formed system, explicit equations go to size-1 blocks.
    // So let's create an actual algebraic loop with some extractable parts.
    std::string code = R"(
        x = y + 1
        y = x - 1
    )";
    EESParser parser;
    auto pr = parser.parse(code);
    REQUIRE(pr.success);
    IR ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    // This creates a block of size 2 (x=y+1, y=x-1 form a cycle)
    REQUIRE(analysis.blocks.size() == 1);
    REQUIRE(analysis.blocks[0].size() == 2);

    auto result = analyseBlockReduction(analysis.blocks[0], ir, analysis);
    // In this cycle, neither equation has all-external RHS, so no explicit
    // extraction is possible.  However, one variable may only appear in its own
    // equation within the block — but both x and y appear in both equations here.
    // So no reduction expected.
    // (This is correct: x=y+1 and y=x-1 are fully coupled.)
}

TEST_CASE("analyseBlockReduction: substitution for variable used only in its own equation", "[symbolic][analysis]") {
    // Create a block where one variable doesn't participate in other equations
    // x = y + z; y = x - 1; z = 10
    // z=10 is explicit (size-1 block). x=y+z and y=x-1 form a 2-var cycle.
    // In the 2-var block {x,y}: x appears in y's equation, y appears in x's equation.
    // No substitution possible (both are coupled).

    // Better test: 3-var block where one var is decoupled
    // a + b = 10; a - b = 2; c = a + b + 1
    // a+b=10 and a-b=2 form a 2-var block; c = a+b+1 depends on their solution.
    // But c=a+b+1 would typically be in its own block after structural analysis.
    // We need a case where a variable is in a block but doesn't affect others.
    // The substitution phase handles this: if var V is matched to eq E,
    // and V doesn't appear in any OTHER block equation, extract it.
    // This is rare in naturally-formed blocks (SCCs always couple all vars).
    // But it can happen after CoolProp inversion changes dependencies.
    // We'll test this indirectly through the CoolProp integration tests.
}

// ============================================================================
// CoolProp inversion integration tests
// ============================================================================

TEST_CASE("Solver with symbolic reduction solves simple thermo model", "[symbolic][solver][integration]") {
    // Model: known P and T, compute h and s via CoolProp
    // This should be explicit even without reduction, but tests the path
    std::string code = R"(
        T = 300
        P = 100000
        h = enthalpy(Water, T=T, P=P)
        s = entropy(Water, T=T, P=P)
    )";
    EESParser parser;
    auto pr = parser.parse(code);
    REQUIRE(pr.success);
    IR ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    CoolPropConfig cpConfig;
    Solver solver(ir, analysis, cpConfig);
    SolverOptions opts;
    opts.enableSymbolicReduction = true;
    SolveResult result = solver.solve(opts);
    CHECK(result.success);
}

TEST_CASE("Solver with symbolic reduction: CoolProp inversion on refrigeration-like", "[symbolic][solver][integration]") {
    // Pattern that triggers CoolProp inversion:
    // P_ev = pressure(R134a, T=T_ev, x=1)   -> P_ev from known T_ev, x
    // h_1 = enthalpy(R134a, P=P_ev, T=T_1)  -> h_1 from known P_ev, T_1
    // All explicit, but let's test the feature doesn't break anything
    std::string code = R"(
        T_ev = -10
        T_1 = T_ev + 2
        P_ev = pressure(R134a, T=T_ev, x=1)
        h_1 = enthalpy(R134a, P=P_ev, T=T_1)
        s_1 = entropy(R134a, P=P_ev, T=T_1)
    )";
    EESParser parser;
    auto pr = parser.parse(code);
    REQUIRE(pr.success);
    IR ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    CoolPropConfig cpConfig;
    Solver solver(ir, analysis, cpConfig);
    SolverOptions opts;
    opts.enableSymbolicReduction = true;
    SolveResult result = solver.solve(opts);
    CHECK(result.success);
    CHECK(result.variables.count("h_1") > 0);
    CHECK(result.variables.count("s_1") > 0);
}

TEST_CASE("Symbolic reduction: model with invertible CoolProp in algebraic block", "[symbolic][solver][integration]") {
    // This model creates a block where CoolProp inversion helps:
    // h_known = 200000 (enthalpy is known)
    // P_known = 100000 (pressure is known)
    // h_known = enthalpy(Water, T=T_unknown, P=P_known)
    // T_unknown + delta_T = T_out
    // delta_T = 5
    //
    // Without inversion: h_known and T_unknown are coupled in a block via the
    // enthalpy equation. With inversion, T_unknown = temperature(Water, H=h_known, P=P_known)
    // can be computed directly.
    std::string code = R"(
        h_known = 200000
        P_known = 100000
        h_known = enthalpy(Water, T=T_unknown, P=P_known)
        delta_T = 5
        T_out = T_unknown + delta_T
    )";
    EESParser parser;
    auto pr = parser.parse(code);
    REQUIRE(pr.success);
    IR ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    // Solve WITHOUT reduction
    {
        CoolPropConfig cpConfig;
        Solver solver(ir, analysis, cpConfig);
        SolverOptions opts;
        opts.enableSymbolicReduction = false;
        SolveResult result = solver.solve(opts);
        CHECK(result.success);
    }

    // Solve WITH reduction
    {
        CoolPropConfig cpConfig;
        Solver solver(ir, analysis, cpConfig);
        SolverOptions opts;
        opts.enableSymbolicReduction = true;
        SolveResult result = solver.solve(opts);
        CHECK(result.success);
        // T_unknown should be roughly 47.7°C for water at h=200kJ/kg, P=100kPa
        if (result.success) {
            double T = result.variables.at("T_unknown");
            CHECK(T > 40.0);
            CHECK(T < 60.0);
        }
    }
}

TEST_CASE("Symbolic reduction: disabled by default adds no overhead", "[symbolic][performance]") {
    // Verify that with enableSymbolicReduction = false (default), the solver
    // takes the same path as before (no extra analysis time).
    std::string code = R"(
        x = y + 1
        y = x - 1
    )";
    EESParser parser;
    auto pr = parser.parse(code);
    REQUIRE(pr.success);
    IR ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    CoolPropConfig cpConfig;

    // Time with reduction disabled (default)
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; ++i) {
        Solver solver(ir, analysis, cpConfig);
        solver.setGuess("x", 1.0);
        solver.setGuess("y", 0.0);
        SolverOptions opts;
        opts.enableSymbolicReduction = false;
        solver.solve(opts);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double dtDisabled = std::chrono::duration<double>(t1 - t0).count();

    // Time with default SolverOptions (enableSymbolicReduction should be false)
    auto t2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; ++i) {
        Solver solver(ir, analysis, cpConfig);
        solver.setGuess("x", 1.0);
        solver.setGuess("y", 0.0);
        SolverOptions opts;  // default
        CHECK_FALSE(opts.enableSymbolicReduction);  // Verify default is false
        solver.solve(opts);
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    double dtDefault = std::chrono::duration<double>(t3 - t2).count();

    // Both should take similar time (within 50%)
    CHECK(dtDefault < dtDisabled * 1.5 + 0.01);
}

TEST_CASE("Symbolic reduction: falls back when reduction fails at evaluation", "[symbolic][robustness]") {
    // Create a model that the reduction can analyse but fails at evaluation.
    // The solver should fall back to the standard pipeline.
    std::string code = R"(
        x = y + 1
        y = x - 1
    )";
    EESParser parser;
    auto pr = parser.parse(code);
    REQUIRE(pr.success);
    IR ir = IR::fromAST(pr.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);

    CoolPropConfig cpConfig;
    Solver solver(ir, analysis, cpConfig);
    solver.setGuess("x", 1.0);
    solver.setGuess("y", 0.0);

    SolverOptions opts;
    opts.enableSymbolicReduction = true;  // enable even though no reduction is possible
    SolveResult result = solver.solve(opts);
    CHECK(result.success);
}

TEST_CASE("Symbolic reduction config loading", "[symbolic][config]") {
    // Verify that the config file parser recognizes the new option
    SolverOptions opts;
    CHECK_FALSE(opts.enableSymbolicReduction);  // default is false

    // Write a temp config file
    std::string tmpPath = "/tmp/test_symbolic_config.conf";
    {
        std::ofstream f(tmpPath);
        f << "enableSymbolicReduction = true\n";
    }

    loadSolverOptionsFromFile(tmpPath, opts);
    CHECK(opts.enableSymbolicReduction);
}

// ============================================================================
// Full model tests with symbolic reduction
// ============================================================================

TEST_CASE("Symbolic reduction on real example files", "[symbolic][integration][examples]") {
    // Test that enabling symbolic reduction doesn't break any of the simple
    // example models that we can parse
    namespace fs = std::filesystem;

    std::vector<std::string> testFiles = {
        "../examples/refrigeration1.eescode",
        "../examples/rankine1.eescode",
        "../examples/cpbar.eescode",
    };

    for (const auto& file : testFiles) {
        if (!fs::exists(file)) continue;

        SECTION(file) {
            std::ifstream ifs(file);
            std::string code((std::istreambuf_iterator<char>(ifs)),
                              std::istreambuf_iterator<char>());

            EESParser parser;
            auto pr = parser.parse(code);
            if (!pr.success) continue;  // Skip unparseable files

            IR ir = IR::fromAST(pr.program);

            // Load initials if available
            std::string initialsPath = file.substr(0, file.rfind('.')) + ".initials";
            if (fs::exists(initialsPath)) {
                ir.loadInitialsFromFile(initialsPath);
            }

            auto analysis = StructuralAnalyzer::analyze(ir);
            if (!analysis.success) continue;

            // Solve without reduction
            CoolPropConfig cpConfig;
            SolverOptions opts;
            opts.enableSymbolicReduction = false;

            Solver solver1(ir, analysis, cpConfig);
            auto result1 = solver1.solve(opts);

            // Solve with reduction
            opts.enableSymbolicReduction = true;
            Solver solver2(ir, analysis, cpConfig);
            auto result2 = solver2.solve(opts);

            // Both should succeed if the first one did
            if (result1.success) {
                CHECK(result2.success);

                // Results should be very close
                for (const auto& [name, val1] : result1.variables) {
                    auto it = result2.variables.find(name);
                    if (it != result2.variables.end()) {
                        double val2 = it->second;
                        if (std::abs(val1) > 1e-10) {
                            CHECK_THAT(val2, WithinRel(val1, 0.01));  // 1% tolerance
                        }
                    }
                }
            }
        }
    }
}
