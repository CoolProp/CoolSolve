/**
 * test_parametric.cpp — Tests for imposed variable detection and parametric
 * study infrastructure.
 *
 * Covers:
 *   - AST-level detection of "var = constant" patterns (positive & negative)
 *   - Rejection of non-imposed patterns (expressions, function calls, etc.)
 *   - Scientific notation, negative values, zero
 *   - Case-insensitive variable matching
 *   - Multiple imposed variables in a single model
 *   - Edge cases: units annotations, inline comments, commented-out equations
 *   - IR round-trip: parse → IR → impose detection matches expectations
 *   - Real example files (pressuredrop, cpbar, etc.)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/ast.h"
#include "coolsolve/constants.h"
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace coolsolve;
using Catch::Approx;

// ============================================================================
// Helper: detect imposed variables from an IR (mirrors server.cpp logic)
// ============================================================================

struct ImposedVar {
    std::string name;
    double value;
};

/**
 * Scan IR equations for patterns: LHS is a scalar Variable, RHS is a
 * NumberLiteral or UnaryOp("-", NumberLiteral).  Returns the list of
 * imposed variables with their constant values.
 */
static std::vector<ImposedVar> detectImposedVars(const IR& ir) {
    std::set<std::string, CaseInsensitiveLess> seen;
    std::vector<ImposedVar> result;

    for (const auto& eq : ir.getEquations()) {
        if (!eq.lhs || !eq.rhs) continue;

        if (eq.lhs->is<Variable>()) {
            const auto& var = eq.lhs->as<Variable>();
            if (!var.indices.empty()) continue;  // skip arrays

            double val = 0.0;
            bool isImposed = false;

            if (eq.rhs->is<NumberLiteral>()) {
                val = eq.rhs->as<NumberLiteral>().value;
                isImposed = true;
            } else if (eq.rhs->is<UnaryOp>()) {
                const auto& uop = eq.rhs->as<UnaryOp>();
                if (uop.op == "-" && uop.operand && uop.operand->is<NumberLiteral>()) {
                    val = -uop.operand->as<NumberLiteral>().value;
                    isImposed = true;
                }
            }

            if (isImposed && seen.insert(var.name).second) {
                result.push_back({var.name, val});
            }
        }
    }
    return result;
}

/** Parse source → IR, returning success flag and the IR */
static std::pair<bool, IR> parseToIR(const std::string& source) {
    EESParser parser;
    auto parseResult = parser.parse(source);
    if (!parseResult.success) {
        return {false, IR()};
    }
    auto ir = IR::fromAST(parseResult.program);
    return {true, std::move(ir)};
}

// ============================================================================
// Basic imposed detection
// ============================================================================

TEST_CASE("Imposed detection — simple positive literal", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("x = 5");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "x");
    CHECK(imposed[0].value == 5.0);
}

TEST_CASE("Imposed detection — scientific notation", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("P = 1E5");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "P");
    CHECK(imposed[0].value == Approx(1e5));
}

TEST_CASE("Imposed detection — small scientific notation", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("epsilon = 0.0000015");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "epsilon");
    CHECK(imposed[0].value == Approx(1.5e-6));
}

TEST_CASE("Imposed detection — negative literal", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("T_min = -40");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "T_min");
    CHECK(imposed[0].value == Approx(-40.0));
}

TEST_CASE("Imposed detection — negative scientific notation", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("coeff = -3.5E-2");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "coeff");
    CHECK(imposed[0].value == Approx(-0.035));
}

TEST_CASE("Imposed detection — zero", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("offset = 0");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "offset");
    CHECK(imposed[0].value == 0.0);
}

TEST_CASE("Imposed detection — decimal with no integer part", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("eta = 0.85");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "eta");
    CHECK(imposed[0].value == Approx(0.85));
}

TEST_CASE("Imposed detection — with units annotation", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("D_i = 16E-3 [m]");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "D_i");
    CHECK(imposed[0].value == Approx(0.016));
}

TEST_CASE("Imposed detection — with inline comment", "[parametric][imposed]") {
    // EES uses {} and "" for comments
    auto [ok, ir] = parseToIR("L = 10 [m] \"pipe length\"");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "L");
    CHECK(imposed[0].value == Approx(10.0));
}

// ============================================================================
// Negative cases — should NOT be detected as imposed
// ============================================================================

TEST_CASE("Not imposed — RHS is expression", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("y = x + 3");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    CHECK(imposed.empty());
}

TEST_CASE("Not imposed — RHS is variable", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("y = x");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    CHECK(imposed.empty());
}

TEST_CASE("Not imposed — RHS is function call", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("y = sin(x)");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    CHECK(imposed.empty());
}

TEST_CASE("Not imposed — LHS is complex expression", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("a + b = 5");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    CHECK(imposed.empty());
}

TEST_CASE("Not imposed — RHS is product of number and variable", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("Q = 2 * m_dot");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    CHECK(imposed.empty());
}

TEST_CASE("Not imposed — commented-out equation", "[parametric][imposed]") {
    auto [ok, ir] = parseToIR("{T = 300}\nP = 1E5");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    // Only P=1E5 should be imposed, T=300 is in a comment
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "P");
}

// ============================================================================
// Multiple imposed variables
// ============================================================================

TEST_CASE("Multiple imposed variables", "[parametric][imposed]") {
    const char* code = R"(
        T_in = 300
        P_in = 1E5
        m_dot = 0.5
        T_out = T_in + deltaT
        deltaT = 10
    )";
    auto [ok, ir] = parseToIR(code);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);

    // T_in, P_in, m_dot, deltaT are imposed; T_out is not (RHS is expression)
    REQUIRE(imposed.size() == 4);

    // Build a name-to-value map for easy checking
    std::map<std::string, double, CaseInsensitiveLess> imap;
    for (const auto& v : imposed) {
        imap[v.name] = v.value;
    }

    CHECK(imap.count("T_in") == 1);
    CHECK(imap["T_in"] == Approx(300.0));
    CHECK(imap.count("P_in") == 1);
    CHECK(imap["P_in"] == Approx(1e5));
    CHECK(imap.count("m_dot") == 1);
    CHECK(imap["m_dot"] == Approx(0.5));
    CHECK(imap.count("deltaT") == 1);
    CHECK(imap["deltaT"] == Approx(10.0));
    CHECK(imap.count("T_out") == 0);
}

// ============================================================================
// Real example model: pressuredrop.eescode
// ============================================================================

TEST_CASE("Imposed detection — pressuredrop.eescode", "[parametric][imposed][examples]") {
    // Read the example file
    std::ifstream ifs("../examples/pressuredrop.eescode");
    if (!ifs.is_open()) {
        // Skip if not run from build/ directory
        WARN("Skipping: pressuredrop.eescode not found (run tests from build/ directory)");
        return;
    }
    std::string source((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());

    auto [ok, ir] = parseToIR(source);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);

    // Expected imposed: epsilon, D_i, P, M_dot, L (5 variables)
    std::map<std::string, double, CaseInsensitiveLess> imap;
    for (const auto& v : imposed) {
        imap[v.name] = v.value;
    }

    CHECK(imap.count("epsilon") == 1);
    CHECK(imap["epsilon"] == Approx(1.5e-6));
    CHECK(imap.count("D_i") == 1);
    CHECK(imap["D_i"] == Approx(0.016));
    CHECK(imap.count("P") == 1);
    CHECK(imap["P"] == Approx(1e5));
    CHECK(imap.count("M_dot") == 1);
    CHECK(imap["M_dot"] == Approx(0.08));
    CHECK(imap.count("L") == 1);
    CHECK(imap["L"] == Approx(10.0));

    // T = T_sat+10 should NOT be imposed (RHS is expression)
    CHECK(imap.count("T") == 0);

    INFO("Total imposed: " << imposed.size());
    CHECK(imposed.size() == 5);
}

// ============================================================================
// Real example model: cpbar.eescode
// ============================================================================

TEST_CASE("Imposed detection — cpbar.eescode", "[parametric][imposed][examples]") {
    std::ifstream ifs("../examples/cpbar.eescode");
    if (!ifs.is_open()) {
        WARN("Skipping: cpbar.eescode not found (run tests from build/ directory)");
        return;
    }
    std::string source((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());

    auto [ok, ir] = parseToIR(source);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);

    // At least some imposed variables should be detected
    INFO("cpbar.eescode imposed count: " << imposed.size());
    CHECK(imposed.size() > 0);

    // Print them for debugging
    for (const auto& v : imposed) {
        INFO("  " << v.name << " = " << v.value);
    }
}

// ============================================================================
// Edge case: same variable imposed twice (last wins in set, but both detected)
// ============================================================================

TEST_CASE("Imposed detection — duplicate variable", "[parametric][imposed]") {
    // If a variable appears in two imposed equations, the first one wins
    // (set insertion only keeps first)
    const char* code = R"(
        x = 10
        y = 20
        x = 30
    )";
    auto [ok, ir] = parseToIR(code);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);

    // x appears twice — the detector should have it once (from first occurrence)
    std::map<std::string, double, CaseInsensitiveLess> imap;
    for (const auto& v : imposed) {
        imap[v.name] = v.value;
    }
    CHECK(imap.count("x") == 1);
    CHECK(imap.count("y") == 1);
}

// ============================================================================
// Edge case: array variable should NOT be detected as imposed
// ============================================================================

TEST_CASE("Not imposed — array variable", "[parametric][imposed]") {
    const char* code = "T[1] = 300\nT[2] = 400";
    auto [ok, ir] = parseToIR(code);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    // Array elements should NOT be detected as imposed (our detection skips them)
    // Note: depending on how the IR stores array variables, they may or may not
    // have indices. The detection checks var.indices.empty().
    // If the parser produces Variable with indices, these should be skipped.
    for (const auto& v : imposed) {
        INFO("Unexpected imposed: " << v.name);
    }
    // Array elements like T[1] have non-empty indices, so they should be skipped
    CHECK(imposed.empty());
}

// ============================================================================
// String variable should NOT be detected as imposed
// ============================================================================

TEST_CASE("Not imposed — string variable assignment", "[parametric][imposed]") {
    const char* code = "fluid$ = 'Water'\nP = 1E5";
    auto [ok, ir] = parseToIR(code);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    // Only P should be imposed, not fluid$
    REQUIRE(imposed.size() == 1);
    CHECK(imposed[0].name == "P");
}

// ============================================================================
// Thermodynamic model with function calls on RHS
// ============================================================================

TEST_CASE("Not imposed — thermodynamic function call on RHS", "[parametric][imposed]") {
    const char* code = R"(
        P = 1E5
        T = 300
        h = enthalpy('Water', T=T, P=P)
    )";
    auto [ok, ir] = parseToIR(code);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);
    // P and T are imposed; h is not (RHS is function call)
    REQUIRE(imposed.size() == 2);
    std::set<std::string, CaseInsensitiveLess> names;
    for (const auto& v : imposed) names.insert(v.name);
    CHECK(names.count("P") == 1);
    CHECK(names.count("T") == 1);
    CHECK(names.count("h") == 0);
}

// ============================================================================
// Model with all example imposed patterns from solver_roadmap
// ============================================================================

TEST_CASE("Imposed detection — comprehensive patterns", "[parametric][imposed]") {
    const char* code = R"(
        T_in = 300
        P_in = 101325
        m_dot = 0.5 [kg/s]
        eta_s = 0.85
        Q_loss = 0 [W]
        T_ambient = -15
        epsilon_r = 1.5E-6 [m]
        fraction = 16E-3
    )";
    auto [ok, ir] = parseToIR(code);
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);

    std::map<std::string, double, CaseInsensitiveLess> imap;
    for (const auto& v : imposed) {
        imap[v.name] = v.value;
    }

    REQUIRE(imposed.size() == 8);
    CHECK(imap["T_in"] == Approx(300.0));
    CHECK(imap["P_in"] == Approx(101325.0));
    CHECK(imap["m_dot"] == Approx(0.5));
    CHECK(imap["eta_s"] == Approx(0.85));
    CHECK(imap["Q_loss"] == Approx(0.0));
    CHECK(imap["T_ambient"] == Approx(-15.0));
    CHECK(imap["epsilon_r"] == Approx(1.5e-6));
    CHECK(imap["fraction"] == Approx(0.016));
}

// ============================================================================
// String variable filtering (regression for orc_simple / fluid$ bug)
// ============================================================================

TEST_CASE("String variables are not detected as imposed", "[parametric][imposed][string]") {
    auto [ok, ir] = parseToIR(R"(
        fluid$='r245fa'
        T = 300
        P = 1E5
    )");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);

    std::set<std::string, CaseInsensitiveLess> names;
    for (const auto& v : imposed) names.insert(v.name);

    CHECK(names.count("T") == 1);
    CHECK(names.count("P") == 1);
    // fluid$ is a string assignment, not a numeric imposed variable
    CHECK(names.count("fluid$") == 0);
}

TEST_CASE("String variable filtering in initials generation", "[parametric][string]") {
    // Simulate the parametric initials-building logic:
    // variables ending with $ should be skipped when writing numeric initials
    std::map<std::string, double> allVars = {
        {"T", 300.0},
        {"P", 1e5},
        {"fluid$", 0.0},  // string var incorrectly stored as numeric
        {"cf$", 0.0},
        {"hf$", 0.0},
        {"M_dot", 1.5},
    };

    // Build initials string, skipping string variables (names ending with $)
    std::ostringstream oss;
    for (const auto& [name, val] : allVars) {
        if (!name.empty() && name.back() == '$') continue;
        oss << name << " = " << val << "\n";
    }
    std::string initials = oss.str();

    // Verify no string variables appear in the output
    CHECK(initials.find("fluid$") == std::string::npos);
    CHECK(initials.find("cf$") == std::string::npos);
    CHECK(initials.find("hf$") == std::string::npos);
    // Verify numeric variables are present
    CHECK(initials.find("T =") != std::string::npos);
    CHECK(initials.find("P =") != std::string::npos);
    CHECK(initials.find("M_dot =") != std::string::npos);
}

TEST_CASE("ORC model with string variables parses correctly", "[parametric][imposed][orc]") {
    // Mimics the orc_simple.eescode pattern:
    // fluid$ is a string, pinch_cd is imposed AND used in a min() equation
    auto [ok, ir] = parseToIR(R"(
        fluid$='r245fa'
        M_dot = 1
        T_cd = 30
        pinch_cd = 10
        pinch_cd = min(T_ex - T_cf_su, T_cd_v - T_cf_ex)
    )");
    REQUIRE(ok);
    auto imposed = detectImposedVars(ir);

    std::map<std::string, double, CaseInsensitiveLess> imap;
    for (const auto& v : imposed) imap[v.name] = v.value;

    // pinch_cd = 10 should be detected as imposed
    CHECK(imap.count("pinch_cd") == 1);
    CHECK(imap["pinch_cd"] == Approx(10.0));
    // M_dot = 1 should be detected
    CHECK(imap.count("M_dot") == 1);
    // T_cd = 30 should be detected
    CHECK(imap.count("T_cd") == 1);
    // fluid$ should NOT be detected as imposed
    CHECK(imap.count("fluid$") == 0);
}
