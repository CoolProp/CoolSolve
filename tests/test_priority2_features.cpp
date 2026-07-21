/**
 * Tests for Priority 2 EES-compatibility features:
 *
 *   Math:
 *     - ERF, ERFC
 *     - LMTD_F (log-mean temperature difference)
 *   String:
 *     - STRING$(x)  number  →  string
 *     - STRINGVAL(s$) string → number
 *   Thermophysical properties:
 *     - PRANDTL, SURFACETENSION
 *     - KINEMATICVISCOSITY (derived)
 *     - THERMALDIFFUSIVITY (derived)
 *     - COMPRESSIBILITYFACTOR, ISENTROPICEXPONENT
 *     - T_CRIT, P_CRIT, V_CRIT, T_TRIPLE (pure-fluid constants)
 *     - PHASE$
 *   Fluids:
 *     - Incompressible solutions (MEG, MPG, CaCl2, ...)
 *     - Unsupported-fluid error messages (NH3H2O, LiBrH2O)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "coolsolve/autodiff_node.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/fluids.h"
#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// Helper: parse a single equation, evaluate the RHS with given variable values.
// Uses a SI/Kelvin unit system so T=300 really means 300 K, matching natural
// thermodynamic reference values in this test file.
static coolsolve::ADValue evalRHS(const std::string& source,
                                  const std::vector<std::pair<std::string, double>>& vars = {}) {
    coolsolve::EESParser parser;
    auto parseResult = parser.parse(source);
    REQUIRE(parseResult.success);
    auto ir = coolsolve::IR::fromAST(parseResult.program);
    const auto& equations = ir.getEquations();
    REQUIRE(!equations.empty());

    size_t n = std::max<size_t>(1, vars.size());
    coolsolve::CoolPropConfig cfg;
    cfg.units.temperature = "K";
    coolsolve::ExpressionEvaluator eval(n, cfg);
    for (size_t i = 0; i < vars.size(); ++i) {
        eval.setVariable(vars[i].first,
                         coolsolve::ADValue::independent(vars[i].second, i, n));
    }
    return eval.evaluate(equations[0].rhs);
}

// Variant for string-valued expressions using K unit system
static std::string evalRHSString(const std::string& source) {
    coolsolve::EESParser parser;
    auto parseResult = parser.parse(source);
    REQUIRE(parseResult.success);
    auto ir = coolsolve::IR::fromAST(parseResult.program);
    const auto& equations = ir.getEquations();
    REQUIRE(!equations.empty());

    coolsolve::CoolPropConfig cfg;
    cfg.units.temperature = "K";
    coolsolve::ExpressionEvaluator eval(1, cfg);
    return eval.evaluateString(equations[0].rhs);
}

// ============================================================================
// Math: ERF / ERFC
// ============================================================================

TEST_CASE("ERF and ERFC math functions", "[priority2][math]") {
    using namespace coolsolve;

    SECTION("erf(0) = 0, erf(inf) -> 1") {
        std::vector<ADValue> args = { ADValue::independent(0.0, 0, 1) };
        auto r = evaluateStandardFunction("erf", args);
        REQUIRE_THAT(r.value, WithinAbs(0.0, 1e-12));
        // d/dx erf(0) = 2/sqrt(pi)
        REQUIRE_THAT(r.gradient[0],
                     WithinAbs(2.0 / std::sqrt(M_PI), 1e-12));
    }

    SECTION("erf(1) matches std::erf") {
        std::vector<ADValue> args = { ADValue::independent(1.0, 0, 1) };
        auto r = evaluateStandardFunction("erf", args);
        REQUIRE_THAT(r.value, WithinAbs(std::erf(1.0), 1e-12));
        REQUIRE_THAT(r.gradient[0],
                     WithinAbs(2.0 / std::sqrt(M_PI) * std::exp(-1.0), 1e-12));
    }

    SECTION("erfc(0) = 1, erfc(x) = 1 - erf(x)") {
        std::vector<ADValue> args = { ADValue::independent(0.0, 0, 1) };
        auto r = evaluateStandardFunction("erfc", args);
        REQUIRE_THAT(r.value, WithinAbs(1.0, 1e-12));
        args = { ADValue::independent(1.2, 0, 1) };
        r = evaluateStandardFunction("erfc", args);
        REQUIRE_THAT(r.value, WithinAbs(std::erfc(1.2), 1e-12));
    }

    SECTION("erf via parser") {
        auto r = evalRHS("y = erf(1)", {});
        REQUIRE_THAT(r.value, WithinAbs(std::erf(1.0), 1e-12));
    }
}

// ============================================================================
// Math: LMTD_F
// ============================================================================

TEST_CASE("LMTD_F function", "[priority2][math]") {
    using namespace coolsolve;

    SECTION("Standard case DT1=10, DT2=2") {
        std::vector<ADValue> args = {
            ADValue::independent(10.0, 0, 2),
            ADValue::independent(2.0, 1, 2)
        };
        auto r = evaluateStandardFunction("lmtd_f", args);
        double expected = (10.0 - 2.0) / std::log(10.0 / 2.0);
        REQUIRE_THAT(r.value, WithinAbs(expected, 1e-10));
    }

    SECTION("Degenerate case DT1=DT2 handled smoothly") {
        std::vector<ADValue> args = {
            ADValue::independent(5.0, 0, 2),
            ADValue::independent(5.0, 1, 2)
        };
        auto r = evaluateStandardFunction("lmtd_f", args);
        REQUIRE_THAT(r.value, WithinAbs(5.0, 1e-8));
    }

    SECTION("LMTD_F via parser") {
        auto r = evalRHS("y = LMTD_F(8, 2)", {});
        double expected = (8.0 - 2.0) / std::log(4.0);
        REQUIRE_THAT(r.value, WithinAbs(expected, 1e-10));
    }
}

// ============================================================================
// String: STRING$ / STRINGVAL
// ============================================================================

TEST_CASE("STRINGVAL parses numeric strings", "[priority2][string]") {
    using namespace coolsolve;
    coolsolve::EESParser parser;
    // Build an IR with a string variable and call STRINGVAL
    auto parseResult = parser.parse("y = STRINGVAL('3.14')");
    REQUIRE(parseResult.success);
    auto ir = IR::fromAST(parseResult.program);
    const auto& equations = ir.getEquations();
    REQUIRE(!equations.empty());

    ExpressionEvaluator eval(1);
    auto r = eval.evaluate(equations[0].rhs);
    REQUIRE_THAT(r.value, WithinAbs(3.14, 1e-12));
}

TEST_CASE("STRING$ formats numbers", "[priority2][string]") {
    using namespace coolsolve;
    coolsolve::EESParser parser;
    auto parseResult = parser.parse("s$ = STRING$(42)");
    REQUIRE(parseResult.success);
    auto ir = IR::fromAST(parseResult.program);
    const auto& equations = ir.getEquations();
    REQUIRE(!equations.empty());

    ExpressionEvaluator eval(1);
    std::string s = eval.evaluateString(equations[0].rhs);
    REQUIRE(s == "42");

    parseResult = parser.parse("s$ = STRING$(3.14)");
    REQUIRE(parseResult.success);
    ir = IR::fromAST(parseResult.program);
    const auto& eqs2 = ir.getEquations();
    REQUIRE(!eqs2.empty());
    s = eval.evaluateString(eqs2[0].rhs);
    REQUIRE(s.substr(0, 4) == "3.14");
}

TEST_CASE("STRINGVAL(STRING$(x)) round-trip", "[priority2][string]") {
    using namespace coolsolve;
    auto r = evalRHS("y = STRINGVAL(STRING$(2.718))", {});
    REQUIRE_THAT(r.value, WithinAbs(2.718, 1e-8));
}

// ============================================================================
// Thermophysical properties — PRANDTL, SURFACETENSION
// ============================================================================

TEST_CASE("PRANDTL number for liquid water", "[priority2][thermo]") {
    // Prandtl of liquid water near 300 K ≈ 5.83, but varies slightly by
    // CoolProp version.  Accept a generous range.
    auto r = evalRHS("y = PRANDTL('Water', T=300, P=101325)", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 5.0);
    REQUIRE(r.value < 7.0);
}

TEST_CASE("SURFACETENSION of water at saturation", "[priority2][thermo]") {
    // Surface tension of water at T=300 K, x=0 ~ 0.07 N/m
    auto r = evalRHS("y = SURFACETENSION('Water', T=300, x=0)", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 0.05);
    REQUIRE(r.value < 0.08);
}

// ============================================================================
// Derived properties — KINEMATICVISCOSITY, THERMALDIFFUSIVITY
// ============================================================================

TEST_CASE("KINEMATICVISCOSITY is mu/rho", "[priority2][thermo]") {
    // Compare KINEMATICVISCOSITY('Water', T=300, P=1e5) with
    // viscosity / density.
    auto nu = evalRHS("y = KINEMATICVISCOSITY('Water', T=300, P=101325)", {});
    auto mu = evalRHS("y = VISCOSITY('Water', T=300, P=101325)", {});
    auto rho = evalRHS("y = DENSITY('Water', T=300, P=101325)", {});
    REQUIRE(std::isfinite(nu.value));
    REQUIRE(std::isfinite(mu.value));
    REQUIRE(std::isfinite(rho.value));
    REQUIRE_THAT(nu.value, WithinRel(mu.value / rho.value, 1e-6));
}

TEST_CASE("THERMALDIFFUSIVITY is k/(rho*cp)", "[priority2][thermo]") {
    auto alpha = evalRHS("y = THERMALDIFFUSIVITY('Water', T=300, P=101325)", {});
    auto k = evalRHS("y = CONDUCTIVITY('Water', T=300, P=101325)", {});
    auto rho = evalRHS("y = DENSITY('Water', T=300, P=101325)", {});
    auto cp = evalRHS("y = CP('Water', T=300, P=101325)", {});
    REQUIRE(std::isfinite(alpha.value));
    REQUIRE_THAT(alpha.value, WithinRel(k.value / (rho.value * cp.value), 1e-6));
}

// ============================================================================
// Pure-fluid constants — T_CRIT, P_CRIT, V_CRIT, T_TRIPLE
// ============================================================================

TEST_CASE("T_CRIT of water ~ 647 K", "[priority2][thermo]") {
    auto r = evalRHS("y = T_CRIT('Water')", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 640.0);
    REQUIRE(r.value < 650.0);
}

TEST_CASE("P_CRIT of water ~ 22.06 MPa", "[priority2][thermo]") {
    auto r = evalRHS("y = P_CRIT('Water')", {});
    REQUIRE(std::isfinite(r.value));
    // SI (Pa): ≈ 2.206e7; kPa: ≈ 22060.  Default is Pa in config.
    REQUIRE(r.value > 2.0e7);
    REQUIRE(r.value < 2.3e7);
}

TEST_CASE("V_CRIT of water > 0", "[priority2][thermo]") {
    auto r = evalRHS("y = V_CRIT('Water')", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 0.0);
    REQUIRE(r.value < 0.01); // m^3/kg
}

TEST_CASE("T_TRIPLE of water ~ 273.16 K", "[priority2][thermo]") {
    auto r = evalRHS("y = T_TRIPLE('Water')", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE_THAT(r.value, WithinAbs(273.16, 0.2));
}

// ============================================================================
// Ratio properties — COMPRESSIBILITYFACTOR, ISENTROPICEXPONENT
// ============================================================================

TEST_CASE("COMPRESSIBILITYFACTOR of water vapor < 1", "[priority2][thermo]") {
    auto r = evalRHS("y = COMPRESSIBILITYFACTOR('Water', T=400, P=101325)", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 0.9);
    REQUIRE(r.value <= 1.0);
}

TEST_CASE("ISENTROPICEXPONENT of water vapor > 1", "[priority2][thermo]") {
    auto r = evalRHS("y = ISENTROPICEXPONENT('Water', T=400, P=101325)", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 1.0);
    REQUIRE(r.value < 2.0);
}

// ============================================================================
// PHASE$
// ============================================================================

TEST_CASE("PHASE$ returns human-readable phase", "[priority2][thermo][string]") {
    // At T=400 K, P=101325 Pa, water is superheated vapor
    std::string s = evalRHSString("ph$ = PHASE$('Water', T=400, P=101325)");
    REQUIRE((s == "gas" || s == "supercritical_gas"));
    s = evalRHSString("ph$ = PHASE$('Water', T=300, P=101325)");
    REQUIRE(s == "liquid");
}

// ============================================================================
// Incompressible solutions
// ============================================================================

TEST_CASE("Incompressible solution: MEG density", "[priority2][incomp]") {
    // Ethylene glycol, 30% mass fraction, T=280 K, P=101325 Pa
    // Expected density ~ 1040-1050 kg/m^3 per CoolProp
    auto r = evalRHS("y = DENSITY('MEG', 30, T=280, P=101325)", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 1000.0);
    REQUIRE(r.value < 1100.0);
}

TEST_CASE("Incompressible solution requires concentration", "[priority2][incomp]") {
    // Calling density('MEG', T=..., P=...) without a concentration should fail
    coolsolve::EESParser parser;
    auto parseResult = parser.parse("y = DENSITY('MEG', T=280, P=101325)");
    REQUIRE(parseResult.success);
    auto ir = coolsolve::IR::fromAST(parseResult.program);
    const auto& equations = ir.getEquations();
    REQUIRE(!equations.empty());

    coolsolve::ExpressionEvaluator eval(1);
    REQUIRE_THROWS_AS(eval.evaluate(equations[0].rhs), std::runtime_error);
}

TEST_CASE("Incompressible concentration via named arg C", "[priority2][incomp]") {
    auto r = evalRHS("y = CP('EthyleneGlycol', T=280, P=101325, C=30)", {});
    REQUIRE(std::isfinite(r.value));
    REQUIRE(r.value > 2000.0);
    REQUIRE(r.value < 4500.0);
}

// ============================================================================
// Unsupported-fluid error messages
// ============================================================================

TEST_CASE("NH3H2O is marked unsupported with a clear error", "[priority2][fluids]") {
    coolsolve::EESParser parser;
    auto parseResult = parser.parse("y = DENSITY('NH3H2O', T=300, P=101325)");
    REQUIRE(parseResult.success);
    auto ir = coolsolve::IR::fromAST(parseResult.program);
    const auto& equations = ir.getEquations();
    coolsolve::ExpressionEvaluator eval(1);
    try {
        eval.evaluate(equations[0].rhs);
        FAIL("Expected runtime_error for NH3H2O");
    } catch (const std::exception& e) {
        std::string msg = e.what();
        REQUIRE(msg.find("NH3H2O") != std::string::npos);
    }
}

// ============================================================================
// PSYCHPROPS built-in procedure
// ============================================================================

TEST_CASE("PSYCHPROPS procedure returns 9 psychrometric properties",
          "[priority2][psychprops]") {
    using namespace coolsolve;
    coolsolve::EESParser parser;
    // Moist air at 25°C, 101325 Pa, 50% RH
    std::string src = R"(
T_1 = 298.15
P_1 = 101325
R_1 = 0.5
CALL psychprops(T_1, P_1, R_1: Tout, v, h, s, u, W, R, Twb, Tdp)
)";
    auto parseResult = parser.parse(src);
    REQUIRE(parseResult.success);
    auto ir = IR::fromAST(parseResult.program);

    CoolPropConfig cfg;
    cfg.units.temperature = "K";
    ExpressionEvaluator eval(0, cfg);

    // Evaluate all equations/calls in order at top level
    const auto& eqs = ir.getEquations();
    bool foundProcCall = false;
    for (const auto& eq : eqs) {
        if (eq.procedureCall.has_value()) {
            eval.evaluateProcedureCall(*eq.procedureCall);
            foundProcCall = true;
        } else if (eq.lhs && eq.lhs->is<Variable>()) {
            auto& v = eq.lhs->as<Variable>();
            ADValue val = eval.evaluate(eq.rhs);
            eval.setVariable(v.name, val);
        }
    }
    REQUIRE(foundProcCall);

    // Validate outputs against expected psychrometrics at ~25°C, 50% RH
    REQUIRE_THAT(eval.getVariable("Tout").value, WithinAbs(298.15, 1e-6));
    REQUIRE_THAT(eval.getVariable("R").value, WithinAbs(0.5, 1e-6));
    // Humidity ratio at 25°C, 50% RH ≈ 0.00988 kg/kg dry air
    REQUIRE_THAT(eval.getVariable("W").value, WithinAbs(0.00988, 0.001));
    // Wet-bulb temperature at 25°C, 50% RH ≈ 290.75 K (17.6°C)
    REQUIRE(eval.getVariable("Twb").value > 285.0);
    REQUIRE(eval.getVariable("Twb").value < 295.0);
    // Dew point at 25°C, 50% RH ≈ 287.0 K (~13.9°C)
    REQUIRE(eval.getVariable("Tdp").value > 283.0);
    REQUIRE(eval.getVariable("Tdp").value < 291.0);
    // Enthalpy per kg dry air at 25°C, 50% RH ~ 50 kJ/kg (50000 J/kg SI)
    REQUIRE(eval.getVariable("h").value > 40000.0);
    REQUIRE(eval.getVariable("h").value < 60000.0);
}

// ============================================================================
// FluidRegistry sanity check
// ============================================================================

TEST_CASE("FluidRegistry exposes new priority-2 fluids", "[priority2][fluids]") {
    using namespace coolsolve;
    REQUIRE(FluidRegistry::getFluid("MEG") != nullptr);
    REQUIRE(FluidRegistry::getFluid("EthyleneGlycol") != nullptr);
    REQUIRE(FluidRegistry::getFluid("MPG") != nullptr);
    REQUIRE(FluidRegistry::getFluid("PropyleneGlycol") != nullptr);
    REQUIRE(FluidRegistry::getFluid("CaCl2") != nullptr);
    REQUIRE(FluidRegistry::getFluid("NaCl") != nullptr);
    REQUIRE(FluidRegistry::getFluid("LiBr") != nullptr);
    REQUIRE(FluidRegistry::getFluid("Glycerol") != nullptr);
    
    auto f = FluidRegistry::getFluid("NH3H2O");
    REQUIRE(f != nullptr);
    REQUIRE(f->getType() == FluidType::Unknown);
}
