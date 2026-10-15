/**
 * Regression tests for the gaps and bugs found while importing models into the
 * CoolSolve Library (see docs/model_library_support.md, §4 and §5). Each test
 * case is named after the register entry it covers and is built from the
 * minimal reproducer of that entry.
 *
 * Run with: ./coolsolve_tests "[library-gaps]"
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/autodiff_node.h"
#include "coolsolve/constants.h"
#include "coolsolve/ees_string.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/ir.h"
#include "coolsolve/lookup_table.h"
#include "coolsolve/parser.h"
#include "coolsolve/runner.h"
#include "coolsolve/solution_checker.h"
#include "CoolProp.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace coolsolve;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

namespace fs = std::filesystem;

/// Outcome of running a model through the whole pipeline, like the CLI does.
struct ModelRun {
    bool parseOk = false;
    bool solveOk = false;
    bool verified = false;   ///< every equation re-evaluated with the solution (rel. 1e-3)
    std::string message;     ///< first parse error, or the solver error
    std::vector<std::string> warnings;   ///< messages of all warning diagnostics of the run
    std::vector<Diagnostic> warningItems;   ///< the warning diagnostics themselves (code, line, message)
    std::map<std::string, double, CaseInsensitiveLess> vars;
    std::map<std::string, std::string, CaseInsensitiveLess> strings;   ///< string variables of the solution

    double operator[](const std::string& name) const { return vars.at(name); }
};

/// Write `code` to a temporary model file, solve it and verify the solution.
/// `tables` maps a lookup-table name to its CSV text (written as `<model>-<name>.csv`);
/// `initials`, when not empty, is the text of the model's `.initials` file (guess values).
ModelRun runModel(const std::string& code,
                  const std::map<std::string, std::string>& tables = {},
                  const std::string& initials = "") {
    static int counter = 0;
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_library_gaps";
    fs::create_directories(dir);
    const std::string stem = "model_" + std::to_string(counter++);
    const fs::path file = dir / (stem + ".eescode");
    {
        std::ofstream f(file);
        f << code;
    }
    std::vector<fs::path> tableFiles;
    for (const auto& [name, csv] : tables) {
        tableFiles.push_back(dir / (stem + "-" + name + ".csv"));
        std::ofstream f(tableFiles.back());
        f << csv;
    }

    if (!initials.empty()) {
        tableFiles.push_back(dir / (stem + ".initials"));
        std::ofstream f(tableFiles.back());
        f << initials;
    }

    ModelRun run;
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    options.tolerance = 1e-9;
    run.solveOk = runner.run(options);
    run.parseOk = runner.isParseSuccess();
    for (const auto& d : runner.getDiagnostics().items())
        if (d.severity == DiagnosticSeverity::Warning) {
            run.warnings.push_back(d.message);
            run.warningItems.push_back(d);
        }
    if (!run.parseOk && !runner.getParseResult().errors.empty())
        run.message = runner.getParseResult().errors[0].message;
    if (run.solveOk) {
        const auto& res = runner.getSolveResult();
        run.vars = res.variables;
        run.strings = res.stringVariables;
        run.verified = checkSolution(runner.getIR(), res.variables, res.stringVariables,
                                     options.coolpropConfig, 1e-3,
                                     &runner.getLookupTableStore()).allSatisfied;
    } else {
        run.message = runner.getSolveResult().errorMessage;
    }
    fs::remove(file);
    for (const auto& t : tableFiles) fs::remove(t);
    return run;
}

/// Parse `code` only and return the messages of all parse errors.
std::vector<std::string> parseErrors(const std::string& code) {
    EESParser parser;
    auto result = parser.parse(code);
    std::vector<std::string> messages;
    for (const auto& e : result.errors) messages.push_back(e.message);
    return messages;
}

double molarMassKgPerKmol(const std::string& coolPropFluid) {
    return 1000.0 * CoolProp::PropsSI("M", "T", 300.0, "P", 101325.0, coolPropFluid);
}

bool anyContains(const std::vector<std::string>& messages, const std::string& text) {
    for (const auto& m : messages)
        if (m.find(text) != std::string::npos) return true;
    return false;
}

}  // namespace

// ============================================================================
// CS-BUG-IF-IGNORED: IF conditions inside FUNCTION/PROCEDURE bodies
// ============================================================================

TEST_CASE("CS-BUG-IF-IGNORED: reproducer of the register", "[library-gaps][if-then-else]") {
    SECTION("single-line IF in a FUNCTION with a literal-true condition") {
        auto run = runModel(R"(
FUNCTION f(v)
  f = -1
  IF (5>0) THEN f = 1
END
a = f(0)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
    }

    SECTION("single-line IF ... THEN ... ELSE in a PROCEDURE, both branches") {
        auto run = runModel(R"(
PROCEDURE p(v : r)
  IF (v>0) THEN r = 1 ELSE r = -1
END
CALL p(5 : b)
CALL p(-5 : c)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["b"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(-1.0, 1e-12));
    }
}

TEST_CASE("CS-BUG-IF-IGNORED: block form with ELSE and ENDIF", "[library-gaps][if-then-else]") {
    // Layout of the ULiege cpbar library: THEN followed by a comment, ELSE and
    // ENDIF on their own lines, statements separated by semicolons, brace comment after ELSE.
    const std::string code = R"(
FUNCTION g(f)
  x = 0 ; y = 0
  IF (f>0) Then "products of combustion"
    x = 10
    y = x + 1
  ELSE {pure air}
    x = 20
    y = x + 2
  ENDIF
  g = x + y
END
a = g(1)
b = g(0)
)";
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(21.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(42.0, 1e-12));
}

TEST_CASE("CS-BUG-IF-IGNORED: nested IF and ENDIF in any case", "[library-gaps][if-then-else]") {
    auto run = runModel(R"(
FUNCTION cls(T)
  IF T < 0 THEN
    cls = 1
  ELSE
    IF (T < 100) then
      cls = 2
    else
      cls = 3
    Endif
  ENDIF
END
a = cls(-5)
b = cls(50)
c = cls(150)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(3.0, 1e-12));
}

TEST_CASE("CS-BUG-IF-IGNORED: ELSE followed by a nested IF", "[library-gaps][if-then-else]") {
    auto run = runModel(R"(
FUNCTION sgn(v)
  IF (v>0) THEN sgn = 1 ELSE IF (v<0) THEN sgn = -1 ELSE sgn = 0
END
a = sgn(3)
b = sgn(-3)
c = sgn(0)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(-1.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(0.0, 1e-12));
}

TEST_CASE("CS-BUG-IF-IGNORED: conditions with =, <>, AND, OR and strings", "[library-gaps][if-then-else]") {
    auto run = runModel(R"(
FUNCTION cond(a, b, t$)
  cond = 0
  IF (a = b) THEN cond = cond + 1
  IF (a <> b) THEN cond = cond + 2
  IF (a >= 1) AND (b <= 5) THEN cond = cond + 4
  IF (a > 100) OR (b = 2) THEN cond = cond + 8
  IF (t$ = 'N') THEN cond = cond + 16
  IF (t$ <> 'N') THEN cond = cond + 32
  IF ((a < 0) OR (b < 0)) AND (t$ = 'K') THEN cond = cond + 64
END
r1 = cond(2, 2, 'N')
r2 = cond(1, 3, 'K')
r3 = cond(-1, 7, 'K')
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["r1"], WithinAbs(1 + 4 + 8 + 16, 1e-12));
    CHECK_THAT(run["r2"], WithinAbs(2 + 4 + 32, 1e-12));
    CHECK_THAT(run["r3"], WithinAbs(2 + 32 + 64, 1e-12));
}

TEST_CASE("CS-BUG-IF-IGNORED: the branch not taken is not executed", "[library-gaps][if-then-else]") {
    // Pure-air call of the cpbar library: the untaken branch divides by the zero
    // parameter (EES never executes it). All branches executing gave NaN here.
    auto run = runModel(R"(
FUNCTION inv(v)
  IF (v > 0) THEN
    inv = 1/v
  ELSE
    inv = -1
  ENDIF
END
a = inv(0)
b = inv(4)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinAbs(-1.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(0.25, 1e-12));
}

TEST_CASE("CS-BUG-IF-IGNORED: derivatives flow through the branch taken", "[library-gaps][if-then-else]") {
    // The unknown enters the condition and the taken branch: Newton needs d f/d x.
    auto solve = runModel(R"(
FUNCTION f(v)
  IF (v > 0) THEN
    f = v^2
  ELSE
    f = -v
  ENDIF
END
f(x) = 9
)");
    REQUIRE(solve.solveOk);
    CHECK(solve.verified);
    CHECK_THAT(std::abs(solve["x"]), WithinAbs(3.0, 1e-6));
}

TEST_CASE("CS-BUG-IF-IGNORED: IF statements in DUPLICATE loops of a procedure", "[library-gaps][if-then-else]") {
    auto run = runModel(R"(
FUNCTION countpos(n)
  countpos = 0
  DUPLICATE i = 1, 4
    IF (i > 2) THEN countpos = countpos + 1
  END
END
a = countpos(4)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-BUG-IF-IGNORED: IF statements inside multi-line comments are not executed", "[library-gaps][if-then-else]") {
    // Layout of the commented-out unit-system guard of the cpbar library: the IF on the
    // second line of a multi-line brace comment must not run (it calls unitsystem()/error()).
    auto run = runModel(R"EES(
FUNCTION g(f)
  {Check of the units setting
   IF (f>0)  THEN CALL error('never executed')
   If (f>0)  then CALL error('never executed either')}
  "Another one in a quote comment
   IF (f>0) THEN CALL error('not executed')"
  g = f + 1
END
a = g(1)
)EES");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-BUG-IF-IGNORED: malformed IF statements are reported", "[library-gaps][if-then-else]") {
    SECTION("IF in the main program") {
        auto errs = parseErrors("a = 3\nIF (a>0) THEN b = 1\n");
        REQUIRE(anyContains(errs, "only allowed inside FUNCTION and PROCEDURE"));
    }
    SECTION("block IF in the main program does not leave stray ENDIF errors") {
        auto errs = parseErrors("a = 3\nIF (a>0) THEN\n b = 1\nELSE\n b = 2\nENDIF\n");
        REQUIRE(errs.size() == 1);
        REQUIRE(anyContains(errs, "only allowed inside FUNCTION and PROCEDURE"));
    }
    SECTION("IF in a DUPLICATE of the main program") {
        auto errs = parseErrors("DUPLICATE i = 1, 2\n IF (i>1) THEN x[i] = 1\nEND\n");
        REQUIRE(anyContains(errs, "only allowed inside FUNCTION and PROCEDURE"));
    }
    SECTION("missing ENDIF") {
        auto errs = parseErrors("FUNCTION f(v)\n IF (v>0) THEN\n  f = 1\nEND\na = f(1)\n");
        REQUIRE(anyContains(errs, "without a matching ENDIF"));
    }
    SECTION("ENDIF without IF: ignored with a warning in a body (EES accepts it), an error outside") {
        // CS-GAP-ELSEIF-CHAIN: the ladders of the EES libraries have ENDIFs that close nothing
        EESParser parser;
        const auto inBody = parser.parse("FUNCTION f(v)\n f = 1\n ENDIF\nEND\na = f(1)\n");
        REQUIRE(inBody.errors.empty());
        bool warned = false;
        for (const auto& d : inBody.diagnostics.items())
            if (d.code == "P008" && d.line == 3) warned = true;
        REQUIRE(warned);
        auto errs = parseErrors("DUPLICATE i = 1, 2\n x[i] = 1\n ENDIF\nEND\n");
        REQUIRE(anyContains(errs, "ENDIF without a matching IF"));
    }
    SECTION("ELSE without IF") {
        auto errs = parseErrors("FUNCTION f(v)\n f = 1\n ELSE\nEND\na = f(1)\n");
        REQUIRE(anyContains(errs, "ELSE without a matching IF"));
    }
    SECTION("two ELSE parts") {
        auto errs = parseErrors("FUNCTION f(v)\n IF (v>0) THEN\n f=1\n ELSE\n f=2\n ELSE\n f=3\n ENDIF\nEND\na = f(1)\n");
        REQUIRE(anyContains(errs, "two ELSE parts"));
    }
    SECTION("IF without THEN") {
        auto errs = parseErrors("FUNCTION f(v)\n IF (v>0) f = 1\nEND\na = f(1)\n");
        REQUIRE(anyContains(errs, "IF statement without THEN"));
    }
    SECTION("unparsable condition") {
        auto errs = parseErrors("FUNCTION f(v)\n IF (v >) THEN f = 1\nEND\na = f(1)\n");
        REQUIRE(anyContains(errs, "Could not parse IF condition"));
    }
}

// ============================================================================
// CS-BUG-MOLARMASS: MOLARMASS in kg/kmol, like EES
// ============================================================================

TEST_CASE("CS-BUG-MOLARMASS: molar masses are in kg/kmol", "[library-gaps][molarmass]") {
    auto run = runModel(R"(
mm_CO2 = molarmass(CO2)
mm_N2 = molarmass(N2)
mm_H2O = molarmass(H2O)
mm_air = molarmass(Air_ha)
fluid$ = 'R134a'
mm_R134a = molarmass(fluid$)
R_N2 = 8314/molarmass(N2)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["mm_CO2"], WithinRel(44.01, 1e-3));
    CHECK_THAT(run["mm_N2"], WithinRel(28.013, 1e-3));
    CHECK_THAT(run["mm_H2O"], WithinRel(18.015, 1e-3));
    CHECK_THAT(run["mm_air"], WithinRel(28.965, 1e-3));
    CHECK_THAT(run["mm_R134a"], WithinRel(102.03, 1e-3));
    // gas constant in J/kg-K from the universal constant in J/kmol-K, the EES pattern
    CHECK_THAT(run["R_N2"], WithinRel(296.8, 1e-3));
}

TEST_CASE("CS-BUG-MOLARMASS: unsupported fluids still report an error", "[library-gaps][molarmass]") {
    auto run = runModel("mm = molarmass(NoSuchFluid)\n");
    REQUIRE_FALSE(run.solveOk);
}

// ============================================================================
// CS-BUG-SINGLE-INPUT-PAIR: single-input (ideal-gas) property calls
// ============================================================================

TEST_CASE("CS-BUG-SINGLE-INPUT-PAIR: reproducer of the register", "[library-gaps][single-input]") {
    // EES pattern for mean specific heats; the unit-system-free form of the register
    // reproducer (MOLARMASS is in kg/kmol now). Before the fix the solve pass used an
    // internal pressure of 1000 Pa (clamped) and the verification pass 100 Pa: the
    // equation failed the verification by 6.3e3, and no .sol was written.
    auto pair = runModel(R"(
MM_H2O = molarmass(H2O)
dh_H2O = (enthalpy(H2O,T=200)-enthalpy(H2O,T=25))*MM_H2O
)");
    REQUIRE(pair.solveOk);
    CHECK(pair.verified);

    // The same with separate equations (the old workaround) gives the same number.
    auto split = runModel(R"(
MM_H2O = molarmass(H2O)
h200 = enthalpy(H2O,T=200)
h25 = enthalpy(H2O,T=25)
dh_H2O = (h200-h25)*MM_H2O
)");
    REQUIRE(split.solveOk);
    CHECK(split.verified);
    CHECK_THAT(pair["dh_H2O"], WithinRel(split["dh_H2O"], 1e-12));
    // ~ cp(H2O vapour) * 175 K * 18.015 kg/kmol, in J/kmol
    CHECK_THAT(pair["dh_H2O"], WithinRel(5.98e6, 5e-3));
}

TEST_CASE("CS-BUG-SINGLE-INPUT-PAIR: the same pair with an unknown temperature", "[library-gaps][single-input]") {
    // The temperature is solved by Newton iteration: every residual evaluation of the
    // pair must agree with the verification pass.
    auto run = runModel(R"(
MM_H2O = molarmass(H2O)
dh = 3e6
dh = (enthalpy(H2O,T=T_x)-enthalpy(H2O,T=25))*MM_H2O
MM_N2 = molarmass(N2)
cbar_N2 = (enthalpy(N2,T=T_x)-enthalpy(N2,T=25))/(T_x-25)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["cbar_N2"], WithinRel(1040.0, 2e-2));
}

TEST_CASE("CS-BUG-SINGLE-INPUT-PAIR: the internal pressure raises no unit warning", "[library-gaps][single-input]") {
    // The pressure injected for single-input ideal-gas calls (100 Pa for H2O) is not a
    // user input: "p=100 is Pa, not kPa" was a spurious warning.
    auto run = runModel("h = enthalpy(H2O,T=100)\ncp_N2 = cp(N2,T=50)\n");
    REQUIRE(run.solveOk);
    for (const auto& w : run.warnings)
        CHECK(w.find("not kPa") == std::string::npos);
}

TEST_CASE("CS-BUG-SINGLE-INPUT-PAIR: user-supplied pressures are still checked", "[library-gaps][single-input]") {
    // Negative test: an explicit P=100 (a typical kPa/Pa mix-up) keeps its unit hint.
    auto run = runModel("h = enthalpy(N2,T=100,P=100)\n");
    REQUIRE(run.solveOk);
    bool hinted = false;
    for (const auto& w : run.warnings)
        if (w.find("not kPa") != std::string::npos) hinted = true;
    CHECK(hinted);
}

// ============================================================================
// CS-GAP-UNITSYSTEM-FUNC: UNITSYSTEM('...') and CALL ERROR('...')
// ============================================================================

TEST_CASE("CS-GAP-UNITSYSTEM-FUNC: unit-system guards of EES library routines", "[library-gaps][unitsystem]") {
    // The guards of the ULiege combustion library (CombCmHn_SI_PNG2003_V2.LIB), uncommented:
    // CoolSolve works in degrees Celsius on a mass basis, so neither guard fires.
    auto run = runModel(R"(
FUNCTION g(x)
  IF (unitsystem('K')=1)  THEN CALL error('Please, set C for temperature units!')
  If (unitsystem('Molar')=1)  then CALL error('Please, set Mass basis!')
  g = 2*x
END
a = g(3)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(6.0, 1e-12));
    for (const auto& w : run.warnings) {
        CHECK(w.find("Unknown function 'unitsystem'") == std::string::npos);
        CHECK(w.find("Unknown function 'error'") == std::string::npos);
    }
}

TEST_CASE("CS-GAP-UNITSYSTEM-FUNC: UNITSYSTEM reports the settings in use", "[library-gaps][unitsystem]") {
    auto run = runModel(R"(
u_SI = unitsystem('SI')
u_Eng = unitsystem('Eng')
u_Mass = unitsystem('Mass')
u_Molar = unitsystem('Molar')
u_Deg = unitsystem('Deg')
u_Rad = unitsystem('Rad')
u_C = unitsystem('C')
u_K = unitsystem('K')
u_F = unitsystem('F')
u_Pa = unitsystem('Pa')
u_kPa = unitsystem('kPa')
u_bar = unitsystem('bar')
u_J = unitsystem('J')
u_kJ = unitsystem('kJ')
)");
    REQUIRE(run.solveOk);
    CHECK(run["u_SI"] == 1.0);   CHECK(run["u_Eng"] == 0.0);
    CHECK(run["u_Mass"] == 1.0); CHECK(run["u_Molar"] == 0.0);
    CHECK(run["u_Deg"] == 1.0);  CHECK(run["u_Rad"] == 0.0);
    CHECK(run["u_C"] == 1.0);    CHECK(run["u_K"] == 0.0);   CHECK(run["u_F"] == 0.0);
    CHECK(run["u_Pa"] == 1.0);   CHECK(run["u_kPa"] == 0.0); CHECK(run["u_bar"] == 0.0);
    CHECK(run["u_J"] == 1.0);    CHECK(run["u_kJ"] == 0.0);
}

TEST_CASE("CS-GAP-UNITSYSTEM-FUNC: unknown unit setting is an error", "[library-gaps][unitsystem]") {
    auto run = runModel("a = unitsystem('furlong')\n");
    REQUIRE_FALSE(run.solveOk);
    CHECK(run.message.find("unknown unit setting 'furlong'") != std::string::npos);
}

TEST_CASE("CS-GAP-UNITSYSTEM-FUNC: CALL ERROR stops the calculation with its message", "[library-gaps][unitsystem]") {
    auto bad = runModel(R"(
FUNCTION g(x)
  IF (x < 0) THEN CALL error('x must be positive', x)
  g = sqrt(x)
END
a = g(-4)
)");
    REQUIRE_FALSE(bad.solveOk);
    CHECK(bad.message.find("CALL ERROR: x must be positive") != std::string::npos);
    CHECK(bad.message.find("-4") != std::string::npos);   // the extra argument is reported

    auto good = runModel(R"(
FUNCTION g(x)
  IF (x < 0) THEN CALL error('x must be positive', x)
  g = sqrt(x)
END
a = g(4)
)");
    REQUIRE(good.solveOk);
    CHECK_THAT(good["a"], WithinAbs(2.0, 1e-12));
}

// ============================================================================
// CS-GAP-FORMATION-ENTHALPY: ideal-gas enthalpies include the heats of formation
// ============================================================================

TEST_CASE("CS-GAP-FORMATION-ENTHALPY: heats of combustion from species enthalpies", "[library-gaps][formation-enthalpy]") {
    // The register reproducer: lower heating values built from the enthalpies of the species at
    // the reference temperature. CoolProp's own reference states gave -5.53 MJ/kmol for CO.
    auto run = runModel(R"(
MM_CO = molarmass(CO) ; MM_O2 = molarmass(O2) ; MM_CO2 = molarmass(CO2)
MM_H2 = molarmass(H2) ; MM_H2O = molarmass(H2O)
T_ref = 25
LHV_CO_mol = enthalpy(CO,T=T_ref)*MM_CO + 0.5*enthalpy(O2,T=T_ref)*MM_O2 - enthalpy(CO2,T=T_ref)*MM_CO2
LHV_H2_mol = enthalpy(H2,T=T_ref)*MM_H2 + 0.5*enthalpy(O2,T=T_ref)*MM_O2 - enthalpy(H2O,T=T_ref)*MM_H2O
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // NASA-Glenn heats of formation (thermo.inp): CO2 -393 510, CO -110 535.196, H2O(g) -241 826 kJ/kmol,
    // hence 282 974.8 and 241 826 (EES 10, same data with R = 8.314: about 282 957 and 241 811). The Cengel & Boles
    // table A-26 values used before gave 282 990 and 241 820 (EES 7.966, stored Q_4 of the cpbar test file of the
    // ULiege library: 282 989.9).
    CHECK_THAT(run["LHV_CO_mol"], WithinRel(282.974804e6, 1e-6));
    CHECK_THAT(run["LHV_H2_mol"], WithinRel(241.826e6, 1e-6));
}

TEST_CASE("CS-GAP-FORMATION-ENTHALPY: h(25 C) is the enthalpy of formation", "[library-gaps][formation-enthalpy]") {
    auto run = runModel(R"(
hf_CO2 = enthalpy(CO2,T=25)*molarmass(CO2)
hf_CO = enthalpy(CO,T=25)*molarmass(CO)
hf_H2O = enthalpy(H2O,T=25)*molarmass(H2O)
hf_CH4 = enthalpy(CH4,T=25)*molarmass(CH4)
hf_C2H6 = enthalpy(C2H6,T=25)*molarmass(C2H6)
hf_C3H8 = enthalpy(C3H8,T=25)*molarmass(C3H8)
h_N2 = enthalpy(N2,T=25)
h_O2 = enthalpy(O2,T=25)
h_H2 = enthalpy(H2,T=25)
h_He = enthalpy(He,T=25)
h_Ar = enthalpy(Ar,T=25)
u_CO2 = intenergy(CO2,T=25,P=101325)
)");
    REQUIRE(run.solveOk);
    // J/kmol: standard enthalpies of formation at 298.15 K of the NASA-Glenn thermodynamic data (McBride, Zehe,
    // Gordon, NASA/TP-2002-211556; thermo.inp of NASA CEA, field "heat of formation" of the species header, J/mol),
    // the data family of EES (EES 10 help: CO2 -393 486, H2O -241 811, CH4 -74 595 kJ/kmol, i.e. the NASA
    // polynomials with R = 8.314). The Cengel & Boles values used before were -393 520, -110 530, -241 820,
    // -74 850, -84 680, -103 850 (JANAF-derived, rounded; ethane and propane 0.83 MJ/kmol away from NASA).
    CHECK_THAT(run["hf_CO2"], WithinRel(-393.510e6, 1e-9));
    CHECK_THAT(run["hf_CO"], WithinRel(-110.535196e6, 1e-9));
    CHECK_THAT(run["hf_H2O"], WithinRel(-241.826e6, 1e-9));
    CHECK_THAT(run["hf_CH4"], WithinRel(-74.600e6, 1e-9));
    CHECK_THAT(run["hf_C2H6"], WithinRel(-83.851544e6, 1e-9));
    CHECK_THAT(run["hf_C3H8"], WithinRel(-104.680e6, 1e-9));
    // elements in their standard state (He and Ar too: EES help, "reference state for specific enthalpy is 0 at 25 C")
    CHECK_THAT(run["h_N2"], WithinAbs(0.0, 1e-3));
    CHECK_THAT(run["h_O2"], WithinAbs(0.0, 1e-3));
    CHECK_THAT(run["h_H2"], WithinAbs(0.0, 1e-3));
    CHECK_THAT(run["h_He"], WithinAbs(0.0, 1e-3));
    CHECK_THAT(run["h_Ar"], WithinAbs(0.0, 1e-3));
    // the internal energy follows the same reference (u = h - R T)
    CHECK(run["u_CO2"] < run["hf_CO2"] / molarMassKgPerKmol("CarbonDioxide"));
}

TEST_CASE("CS-GAP-FORMATION-ENTHALPY: enthalpy differences and inversion are untouched", "[library-gaps][formation-enthalpy]") {
    // Only the reference of the enthalpy changes: enthalpy differences are CoolProp's and an enthalpy
    // given as input is shifted back. (The entropy has its own absolute reference, see CS-GAP-IDEAL-ENTROPY.)
    auto run = runModel(R"(
dh_CO2 = enthalpy(CO2,T=200) - enthalpy(CO2,T=25)
h_200 = enthalpy(CO2,T=200,P=101325)
T_back = temperature(CO2,P=101325,h=h_200)
s_CO2 = entropy(CO2,T=25,P=101325)
h_air = enthalpy(Air,T=25,P=101325)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    const double P = 100000.0 + 1325.0;
    const double dhRef = CoolProp::PropsSI("H", "T", 473.15, "P", P, "CarbonDioxide") -
                         CoolProp::PropsSI("H", "T", 298.15, "P", P, "CarbonDioxide");
    CHECK_THAT(run["dh_CO2"], WithinRel(dhRef, 5e-4));   // 101325 Pa internal pressure vs the same
    CHECK_THAT(run["T_back"], WithinAbs(200.0, 1e-6));
    // the entropy is the third-law one: s(25 C, 1 atm) = s0(1 bar) - R ln(1.01325) = 213.786 - 0.109 kJ/kmol-K,
    // minus the real-gas residual of CO2 at 1 atm (0.1 kJ/kmol-K, CoolProp's equation of state)
    CHECK_THAT(run["s_CO2"] * molarMassKgPerKmol("CarbonDioxide") / 1000.0, WithinAbs(213.786 - 0.109, 0.2));
    // Air has no formation enthalpy (EES references it to 0 K): it keeps CoolProp's reference state
    CHECK_THAT(run["h_air"], WithinRel(CoolProp::PropsSI("H", "T", 298.15, "P", P, "Air"), 1e-9));
}

TEST_CASE("CS-GAP-FORMATION-ENTHALPY: the cpbar library reproduces the EES test file", "[library-gaps][formation-enthalpy]") {
    // examples/cpbar.eescode is the ULiege combustion library with the test program of the original
    // EES file (CombCmHn_SI_PNG2003_V2_test.EES, EES 7.966). Stored EES solution, f = 0.07:
    // c_bar_p 1090.598, Q_4 1.517249e8, x 0.9251884, e_min -0.3359375, e -0.02513201;
    // pure air: gamma 1.394538, MM_prod 28.85006.
    const fs::path example = fs::path("..") / "examples" / "cpbar.eescode";
    if (!fs::exists(example)) SKIP("examples/cpbar.eescode not found (run from the build folder)");
    std::ifstream in(example);
    const std::string code((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["c_bar_p"], WithinRel(1090.598, 3e-3));   // ideal-gas tables: EES vs CoolProp
    CHECK_THAT(run["Q_4"], WithinRel(1.517249e8, 2e-3));     // needs the heats of formation
    CHECK_THAT(run["x"], WithinRel(0.9251884, 1e-3));
    CHECK_THAT(run["e_min"], WithinRel(-0.3359375, 1e-9));
    CHECK_THAT(run["e"], WithinRel(-0.02513201, 2e-3));
    CHECK_THAT(run["gamma"], WithinRel(1.394538, 2e-3));
    CHECK_THAT(run["MM_prod"], WithinRel(28.85006, 1e-4));
}

TEST_CASE("CS-GAP-FORMATION-ENTHALPY: He and Ar are elements, h = 0 at 25 C, as in EES", "[library-gaps][formation-enthalpy]") {
    // EES help (He, Ar): "The reference state for specific enthalpy is 0 at 25 C". They kept CoolProp's own enthalpy
    // reference (+6.2 MJ/kmol at 25 C and 1 atm, both) before the elements were given a zero heat of formation.
    // Monatomic ideal gases: h(100 C) - h(25 C) = 5/2 R 75 K, to the real-gas residual of CoolProp at 1 atm
    // (Ar 0.18 %, He less than 0.1 %).
    auto run = runModel(R"(
h_He = molarmass(He)*enthalpy(He,T=25)
h_Ar = molarmass(Ar)*enthalpy(Ar,T=25)
dh_He = molarmass(He)*(enthalpy(He,T=100) - enthalpy(He,T=25))
dh_Ar = molarmass(Ar)*(enthalpy(Ar,T=100) - enthalpy(Ar,T=25))
T_He = temperature(He,h=enthalpy(He,T=100))
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["h_He"], WithinAbs(0.0, 1e-3));
    CHECK_THAT(run["h_Ar"], WithinAbs(0.0, 1e-3));
    const double dhIdeal = 2.5 * 8314.462618 * 75.0;   // J/kmol between 25 C and 100 C
    CHECK_THAT(run["dh_He"], WithinRel(dhIdeal, 1e-3));
    CHECK_THAT(run["dh_Ar"], WithinRel(dhIdeal, 3e-3));
    CHECK_THAT(run["T_He"], WithinAbs(100.0, 1e-6));   // the shift is undone when the enthalpy is an input
}

TEST_CASE("CS-GAP-FORMATION-ENTHALPY: Gibbs energies of reaction share one source with the NASA-Glenn tables", "[library-gaps][formation-enthalpy][ideal-entropy]") {
    // Standard Gibbs energy of the species g = M*h(T) - T*(M*s(T,P) + R ln(P/1 bar)) [J/kmol], from the enthalpy
    // (heat of formation + CoolProp's sensible part) and the absolute entropy of CoolSolve, at P = 1500 Pa where
    // every species is an ideal-gas vapour. Reaction dG0 = sum(nu g). Reference: the same quantities from the
    // NASA-Glenn 9-coefficient polynomials of thermo.inp (NASA/TP-2002-211556, R = 8.314462618 J/mol-K):
    //   g(T) = dHf(298.15 K) + [H(T) - H(298.15 K)] - T S(T)       (dHf: header of the species, S: 1 bar)
    // computed in Python (H, S: integrals of the polynomial of the interval containing T; extrapolated to 298.15 K
    // for C2H6, C3H8 whose first interval starts at 300 K). dG0 in MJ/kmol:
    //                          298.15 K     2000 K      3000 K
    //   CO2 -> CO + 1/2 O2      257.2007    110.2552     27.5656      (ln K = -6.6303, -1.1051)
    //   H2O -> H2 + 1/2 O2      228.5803    135.5622     77.4940      (ln K = -8.1522, -3.1068)
    //   CH4 + 2 O2 -> CO2 + 2 H2O    -801.0044
    //   C2H6 + 7/2 O2 -> 2 CO2 + 3 H2O  -1442.6051
    //   C3H8 + 5 O2 -> 3 CO2 + 4 H2O    -2073.1526
    // At 298.15 K the specific heats do not enter: dG0 = dH0 - T dS0 is fixed by the heats of formation and the
    // standard entropies, so a heat of formation from another source shows up directly (Cengel & Boles instead of
    // NASA: CO2 -> CO + 1/2 O2 +15 kJ/kmol, H2O -> H2 + 1/2 O2 -6, CH4 +252, C2H6 +829, C3H8 -829; CoolProp's own
    // entropy references: 2 to 120 MJ/kmol). At 2000 and 3000 K the test also contains CoolProp's ideal-gas specific
    // heats against NASA's (about 0.1 % of the enthalpy rise) and cannot tell the formation enthalpies of the CO2 and
    // H2O reactions apart (shift of ln K: 0.0009 and 0.0004): it guards the common basis of h and s.
    const double P = 1500.0;
    auto gibbs = [&](double TK, const std::vector<std::string>& species) {
        std::ostringstream code;
        code << std::setprecision(12);
        const double TC = TK - 273.15;
        for (const auto& s : species)
            code << "g_" << s << " = molarmass(" << s << ")*enthalpy(" << s << ",T=" << TC << ") - " << TK
                 << "*(molarmass(" << s << ")*entropy(" << s << ",T=" << TC << ",P=" << P << ") + 8314.462618*ln(" << P
                 << "/1E5))\n";
        return code.str();
    };
    auto solveFor = [&](double TK, const std::vector<std::string>& species, const std::string& reactions) {
        auto run = runModel(gibbs(TK, species) + reactions);
        INFO("T = " << TK << " K: " << run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        return run;
    };

    SECTION("298.15 K: heats of formation and standard entropies") {
        auto run = solveFor(298.15, {"CO2", "CO", "H2O", "H2", "O2", "CH4", "C2H6", "C3H8"},
                            "dG_CO2 = g_CO + 0.5*g_O2 - g_CO2\n"
                            "dG_H2O = g_H2 + 0.5*g_O2 - g_H2O\n"
                            "dG_CH4 = g_CO2 + 2*g_H2O - g_CH4 - 2*g_O2\n"
                            "dG_C2H6 = 2*g_CO2 + 3*g_H2O - g_C2H6 - 3.5*g_O2\n"
                            "dG_C3H8 = 3*g_CO2 + 4*g_H2O - g_C3H8 - 5*g_O2\n");
        // J/kmol. 3 kJ/kmol of rounding and residuals, plus 10 kJ/kmol per H2O of the reaction: CoolProp's Water at
        // 1500 Pa and 25 C (47 % of saturation) is a real vapour, its residual entropy -0.03 kJ/kmol-K is -T s
        // = +9.4 kJ/kmol in g (measured: dG0 of H2O -> H2 + 1/2 O2 is 9.4 off, CH4 + 2 O2 18.8, C2H6 28.2,
        // C3H8 37.6; CO2 -> CO + 1/2 O2: 0.1). The other sources are 15 kJ/kmol (CO2) to 829 kJ/kmol away.
        CHECK_THAT(run["dG_CO2"], WithinAbs(257.2007e6, 3e3));
        CHECK_THAT(run["dG_H2O"], WithinAbs(228.5803e6, 13e3));
        CHECK_THAT(run["dG_CH4"], WithinAbs(-801.0044e6, 23e3));
        CHECK_THAT(run["dG_C2H6"], WithinAbs(-1442.6051e6, 33e3));
        CHECK_THAT(run["dG_C3H8"], WithinAbs(-2073.1526e6, 43e3));
    }
    SECTION("2000 K and 3000 K: equilibrium constants of CO2 and H2O dissociation") {
        // measured deviations from the NASA reference: ln K(CO2) 0.0003 and 0.0009, ln K(H2O) 0.004 and 0.017
        const double R = 8314.462618;
        struct Point { double TK, lnK_CO2, lnK_H2O, tol_CO2, tol_H2O; };
        for (const Point& pt : {Point{2000.0, -6.6303, -8.1522, 2e-3, 1e-2}, Point{3000.0, -1.1051, -3.1068, 3e-3, 3e-2}}) {
            auto run = solveFor(pt.TK, {"CO2", "CO", "H2O", "H2", "O2"},
                                "dG_CO2 = g_CO + 0.5*g_O2 - g_CO2\n"
                                "dG_H2O = g_H2 + 0.5*g_O2 - g_H2O\n");
            INFO("T = " << pt.TK << " K");
            CHECK_THAT(-run["dG_CO2"] / (R * pt.TK), WithinAbs(pt.lnK_CO2, pt.tol_CO2));
            CHECK_THAT(-run["dG_H2O"] / (R * pt.TK), WithinAbs(pt.lnK_H2O, pt.tol_H2O));
        }
    }
}

// ============================================================================
// CS-GAP-IDEAL-ENTROPY: the ideal gases have absolute (third-law) entropies
// ============================================================================

namespace {

/// Standard molar entropy at 298.15 K and 1 bar [kJ/kmol-K], NASA-Glenn thermodynamic data
/// (McBride, Zehe, Gordon, NASA/TP-2002-211556; github.com/nasa/cea data/thermo.inp, R = 8.314462618).
/// NIST-JANAF differs by at most 0.12 (CH4); EES is built on the same NASA data.
struct StandardEntropy {
    const char* ees;        ///< EES name of the ideal gas
    const char* coolProp;   ///< CoolProp fluid behind it
    double s0;              ///< [kJ/kmol-K]
};
const StandardEntropy kStandardEntropies[] = {
    {"N2", "Nitrogen", 191.6086},    {"O2", "Oxygen", 205.1483},     {"H2", "Hydrogen", 130.6803},
    {"He", "Helium", 126.1530},      {"Ar", "Argon", 154.8458},      {"CO2", "CarbonDioxide", 213.7862},
    {"CO", "CarbonMonoxide", 197.6587}, {"CH4", "Methane", 186.3700}, {"C2H6", "Ethane", 229.2198},
    {"C3H8", "Propane", 270.3133},   {"H2O", "Water", 188.8280},
};
constexpr double kGasConstant = 8.314462618;   // [kJ/kmol-K]

const StandardEntropy& standardEntropy(const std::string& ees) {
    for (const auto& e : kStandardEntropies)
        if (ees == e.ees) return e;
    FAIL("no standard entropy for " << ees);
    return kStandardEntropies[0];
}

}  // namespace

TEST_CASE("CS-GAP-IDEAL-ENTROPY: s(25 C) of the 11 species is the third-law entropy minus R ln(P/1 bar)", "[library-gaps][ideal-entropy]") {
    // Low pressure (1500 Pa, above the 1 kPa input floor, H2O is vapour at 25 C below 3.17 kPa): the real-gas
    // residual of CoolProp's equation of state vanishes (<= 0.03 kJ/kmol-K, H2O), what is left is the reference.
    // The molar entropy is M*s [kJ/kmol-K]. The CoolProp reference states were 14 to 145 kJ/kmol-K away
    // (H2 -23, He -14, CO -85.5, CO2 -93.2, CH4 -79.1, C2H6 -126.9, C3H8 -144.8, H2O -63.3).
    std::string code;
    for (const auto& e : kStandardEntropies) {
        const std::string n = e.ees;
        code += "s_" + n + " = molarmass(" + n + ")*entropy(" + n + ",T=25,P=1500)/1000\n";
    }
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (const auto& e : kStandardEntropies) {
        INFO("species " << e.ees);
        const double expected = e.s0 - kGasConstant * std::log(1500.0 / 1.0e5);
        CHECK_THAT(run[std::string("s_") + e.ees], WithinAbs(expected, 0.05));
    }
}

TEST_CASE("CS-GAP-IDEAL-ENTROPY: at 1 bar the residual of the real-gas equation of state is the only difference", "[library-gaps][ideal-entropy]") {
    // N2, O2, Ar, CO2, ... evaluated with CoolProp's real-gas equation of state at 1 bar: s differs from the
    // ideal-gas standard entropy by the (negative) residual entropy, 0.02 (N2) to 0.29 (C3H8) kJ/kmol-K.
    // H2O is liquid at 25 C and 1 bar (the ideal vapour is called at low pressure, previous test).
    std::string code;
    for (const auto& e : kStandardEntropies) {
        if (std::string(e.ees) == "H2O") continue;
        const std::string n = e.ees;
        code += "s_" + n + " = molarmass(" + n + ")*entropy(" + n + ",T=25,P=1E5)/1000\n";
    }
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    for (const auto& e : kStandardEntropies) {
        if (std::string(e.ees) == "H2O") continue;
        INFO("species " << e.ees);
        const double s = run[std::string("s_") + e.ees];
        CHECK(s <= e.s0 + 0.01);
        CHECK(s >= e.s0 - 0.35);
    }
}

TEST_CASE("CS-GAP-IDEAL-ENTROPY: entropy of the ideal gases depends on the pressure as -R ln(P)", "[library-gaps][ideal-entropy]") {
    auto run = runModel(R"(
M_He = molarmass(He)
M_CO2 = molarmass(CO2)
M_CH4 = molarmass(CH4)
ds_He = M_He*(entropy(He,T=100,P=1E6) - entropy(He,T=100,P=1E5))/1000
ds_CO2 = M_CO2*(entropy(CO2,T=25,P=6000) - entropy(CO2,T=25,P=1500))/1000
ds_CH4 = M_CH4*(entropy(CH4,T=25,P=6000) - entropy(CH4,T=25,P=1500))/1000
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["ds_He"], WithinAbs(-kGasConstant * std::log(10.0), 0.01));   // real-gas residual of He at 10 bar: 0.004
    CHECK_THAT(run["ds_CO2"], WithinAbs(-kGasConstant * std::log(4.0), 0.01));
    CHECK_THAT(run["ds_CH4"], WithinAbs(-kGasConstant * std::log(4.0), 0.01));
}

TEST_CASE("CS-GAP-IDEAL-ENTROPY: isentropic round trip with the entropy as an input", "[library-gaps][ideal-entropy]") {
    // temperature(X, s=, P=) and enthalpy(X, s=, P=): the entropy given as input is on the absolute basis, and
    // the temperature of the isentropic state is CoolProp's, whatever the reference.
    auto run = runModel(R"(
s_CO2 = entropy(CO2,T=25,P=1E5)
T_CO2 = temperature(CO2,s=s_CO2,P=3E5)
h_CO2 = enthalpy(CO2,s=s_CO2,P=3E5)
h_CO2_T = enthalpy(CO2,T=T_CO2,P=3E5)
s_back_CO2 = entropy(CO2,T=T_CO2,P=3E5)
s_H2O = entropy(H2O,T=150,P=5000)
T_H2O = temperature(H2O,s=s_H2O,P=10000)
h_H2O = enthalpy(H2O,s=s_H2O,P=10000)
h_H2O_T = enthalpy(H2O,T=T_H2O,P=10000)
s_back_H2O = entropy(H2O,T=T_H2O,P=10000)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // the isentropic temperatures, from CoolProp on its own basis
    const double sCO2 = CoolProp::PropsSI("S", "T", 298.15, "P", 1e5, "CarbonDioxide");
    CHECK_THAT(run["T_CO2"] + 273.15, WithinAbs(CoolProp::PropsSI("T", "S", sCO2, "P", 3e5, "CarbonDioxide"), 1e-6));
    const double sH2O = CoolProp::PropsSI("S", "T", 423.15, "P", 5000.0, "Water");
    CHECK_THAT(run["T_H2O"] + 273.15, WithinAbs(CoolProp::PropsSI("T", "S", sH2O, "P", 10000.0, "Water"), 1e-6));
    // enthalpy and entropy of the recovered state are the ones of the input state
    CHECK_THAT(run["h_CO2"], WithinRel(run["h_CO2_T"], 1e-9));
    CHECK_THAT(run["s_back_CO2"], WithinRel(run["s_CO2"], 1e-9));
    CHECK_THAT(run["h_H2O"], WithinRel(run["h_H2O_T"], 1e-9));
    CHECK_THAT(run["s_back_H2O"], WithinRel(run["s_H2O"], 1e-9));
}

TEST_CASE("CS-GAP-IDEAL-ENTROPY: a literal absolute entropy as input", "[library-gaps][ideal-entropy]") {
    // He is almost ideal at 1 bar: the third-law entropy of the tables (126.153 kJ/kmol-K at 298.15 K, 1 bar)
    // gives 25 C back, and a state at 100 C the entropy of the tables there.
    const double M = molarMassKgPerKmol("Helium");
    std::ostringstream code;
    code << std::setprecision(12)
         << "T_25 = temperature(He, s=" << 126.1530e3 / M << ", P=1E5)\n"
         << "s_100 = molarmass(He)*entropy(He,T=100,P=1E5)/1000\n";
    auto run = runModel(code.str());
    REQUIRE(run.solveOk);
    CHECK_THAT(run["T_25"], WithinAbs(25.0, 0.02));   // real-gas residual of He at 1 bar: 0.0004 kJ/kmol-K = 0.006 K
    // cp = 5/2 R: s(373.15 K) = s(298.15 K) + 5/2 R ln(373.15/298.15)
    CHECK_THAT(run["s_100"], WithinAbs(126.1530 + 2.5 * kGasConstant * std::log(373.15 / 298.15), 5e-3));
}

TEST_CASE("CS-GAP-IDEAL-ENTROPY: entropy and Gibbs energy of the reaction CH4 + 2 O2 -> CO2 + 2 H2O", "[library-gaps][ideal-entropy]") {
    // Species entropies are combined: 25 C, every species at the same pressure (2000 Pa, where H2O is a vapour),
    // 3 mol of gas on each side so that the entropy of reaction does not depend on the pressure.
    // NASA: dS0 = -5.224 kJ/kmol-K (NIST-JANAF -5.08); the CoolProp references gave about -146.
    // Gibbs energy dG0 = dH0 - T dS0 with the heats of formation of the species enthalpies (NASA-Glenn thermo.inp:
    // CO2 -393 510, H2O(g) -241 826, CH4 -74 600 kJ/kmol): dH0 = -802.562 MJ/kmol, dG0 = -802.562 + 298.15*5.2244e-3
    // = -801.004 MJ/kmol. The tabulated dG0f of the textbooks (CO2 -394.359, H2O(g) -228.572, CH4 -50.72 kJ/mol)
    // give -800.78, 0.03 % higher: CH4 has dHf = -74 873 in NIST-JANAF (Cengel -74 850, NASA 2002 -74 600).
    auto run = runModel(R"(
P = 2000
T = 25
MM_CH4 = molarmass(CH4) ; MM_O2 = molarmass(O2) ; MM_CO2 = molarmass(CO2) ; MM_H2O = molarmass(H2O)
dS = (MM_CO2*entropy(CO2,T=T,P=P) + 2*MM_H2O*entropy(H2O,T=T,P=P) - MM_CH4*entropy(CH4,T=T,P=P) - 2*MM_O2*entropy(O2,T=T,P=P))/1000
dH = MM_CO2*enthalpy(CO2,T=T) + 2*MM_H2O*enthalpy(H2O,T=T) - MM_CH4*enthalpy(CH4,T=T) - 2*MM_O2*enthalpy(O2,T=T)
dG = dH - (T+273.15)*dS*1000
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["dS"], WithinAbs(-5.224, 0.2));      // kJ/kmol-K
    CHECK_THAT(run["dH"], WithinRel(-802.562e6, 1e-6));  // J/kmol, NASA-Glenn heats of formation (Cengel: -802.31)
    CHECK_THAT(run["dG"], WithinRel(-801.004e6, 1e-4));  // J/kmol
}

TEST_CASE("CS-GAP-IDEAL-ENTROPY: Air, real fluids and humid air keep CoolProp's reference state", "[library-gaps][ideal-entropy]") {
    // Air is not a third-law species in EES either; Air_ha, the real fluids and AirH2O are not ideal-gas species.
    const double P = 101325.0;
    auto run = runModel(R"(
s_air = entropy(Air,T=25,P=101325)
s_air_ha = entropy(Air_ha,T=25,P=101325)
s_N2_real = entropy(Nitrogen,T=25,P=101325)
s_CO2_real = entropy(CarbonDioxide,T=25,P=101325)
s_water = entropy(Water,T=25,P=101325)
h_air = enthalpy(Air,T=25,P=101325)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["s_air"], WithinRel(CoolProp::PropsSI("S", "T", 298.15, "P", P, "Air"), 1e-9));
    CHECK_THAT(run["s_air_ha"], WithinRel(CoolProp::PropsSI("S", "T", 298.15, "P", P, "Air"), 1e-9));
    CHECK_THAT(run["s_N2_real"], WithinRel(CoolProp::PropsSI("S", "T", 298.15, "P", P, "Nitrogen"), 1e-9));
    CHECK_THAT(run["s_CO2_real"], WithinRel(CoolProp::PropsSI("S", "T", 298.15, "P", P, "CarbonDioxide"), 1e-9));
    CHECK_THAT(run["s_water"], WithinRel(CoolProp::PropsSI("S", "T", 298.15, "P", P, "Water"), 1e-9));
    CHECK_THAT(run["h_air"], WithinRel(CoolProp::PropsSI("H", "T", 298.15, "P", P, "Air"), 1e-9));
}

// ============================================================================
// CS-GAP-MULTILINE-COMMENT: "..." comment trailing an equation, over several lines
// ============================================================================

TEST_CASE("CS-GAP-MULTILINE-COMMENT: reproducer of the register", "[library-gaps][multiline-comment]") {
    auto run = runModel("x = 1 \"a comment\nspanning two lines\"\ny = 2*x\n");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-GAP-MULTILINE-COMMENT: the text of the comment is never parsed", "[library-gaps][multiline-comment]") {
    // three lines, the middle one looks like an equation
    auto run = runModel(R"EES(
a = 5 "first line
  b = 99
  third line"
c = a + 1
T = 25 "[C]"
d = 2*T "units annotations are not comments that continue
  on the next line"
e = 3 // it's a "quoted word: a C-style comment does not open anything
f = e + 1
)EES");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["d"], WithinAbs(50.0, 1e-12));
    CHECK_THAT(run["f"], WithinAbs(4.0, 1e-12));
    CHECK(run.vars.count("b") == 0);
}

TEST_CASE("CS-GAP-MULTILINE-COMMENT: inside function and procedure bodies", "[library-gaps][multiline-comment]") {
    auto run = runModel(R"EES(
FUNCTION g(x)
  y = x + 1 "comment on the
     next line: z = 99"
  g = 2*y
END
PROCEDURE p(x : r)
  r = x * 3     "wrapped comment of
     a procedure"
END
a = g(1)
CALL p(2 : b)
)EES");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(6.0, 1e-12));
}

TEST_CASE("CS-GAP-MULTILINE-COMMENT: a comment that is never closed is reported", "[library-gaps][multiline-comment]") {
    EESParser parser;
    auto result = parser.parse("a = 1\nb = 2 \"never closed\nc = 3\nd = 4\n");
    // the equations before the comment are kept, the rest is ignored...
    CHECK(result.equationCount == 2);
    // ...and the user is told
    bool warned = false;
    for (const auto& d : result.diagnostics.items())
        if (d.severity == DiagnosticSeverity::Warning &&
            d.message.find("opened by a '\"' at line 2 is never closed") != std::string::npos) warned = true;
    CHECK(warned);

    // a comment that is closed raises no warning
    auto ok = parser.parse("a = 1 \"closed\non the next line\"\nb = 2\n");
    CHECK(ok.equationCount == 2);
    for (const auto& d : ok.diagnostics.items())
        CHECK(d.message.find("never closed") == std::string::npos);
}

// ============================================================================
// CS-BUG-INTERP-DESC: INTERPOLATE with a descending x column
// ============================================================================

TEST_CASE("CS-BUG-INTERP-DESC: reproducer of the register", "[library-gaps][interpolate]") {
    // Compressor-map table of the two-shaft gas turbine: rows in descending order of N_rN.
    // INTERPOLATE('t','N_rN','M',0.95) returned 370 instead of 420.
    const std::map<std::string, std::string> tables = {{"t", "N_rN,M\n1,454\n0.95,420\n0.9,370\n"}};
    auto run = runModel(R"(
M_node = INTERPOLATE('t', 'N_rN', 'M', 0.95)
M_mid = INTERPOLATE('t', 'N_rN', 'M', 0.925)
M_top = INTERPOLATE('t', 'N_rN', 'M', 0.975)
M_below = INTERPOLATE('t', 'N_rN', 'M', 0.5)
M_above = INTERPOLATE('t', 'N_rN', 'M', 1.5)
)", tables);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["M_node"], WithinAbs(420.0, 1e-9));
    CHECK_THAT(run["M_mid"], WithinAbs(395.0, 1e-9));
    CHECK_THAT(run["M_top"], WithinAbs(437.0, 1e-9));
    // outside the table: linear extrapolation through the two end nodes, whatever the order
    // (CS-GAP-INTERP-EXTRAP; flat before): 370 + (0.5 - 0.9) * (420 - 370)/0.05, 454 + (1.5 - 1) * (454 - 420)/0.05
    CHECK_THAT(run["M_below"], WithinAbs(-30.0, 1e-9));
    CHECK_THAT(run["M_above"], WithinAbs(794.0, 1e-9));
}

TEST_CASE("CS-BUG-INTERP-DESC: ascending and descending tables give the same result", "[library-gaps][interpolate]") {
    const std::map<std::string, std::string> asc = {{"t", "N_rN,M\n0.9,370\n0.95,420\n1,454\n"}};
    const std::map<std::string, std::string> desc = {{"t", "N_rN,M\n1,454\n0.95,420\n0.9,370\n"}};
    const std::string model = "y1 = INTERPOLATE('t', 'N_rN', 'M', 0.93)\ny2 = INTERPOLATE('t', 'N_rN', 'M', 0.9)\n"
                              "y3 = INTERPOLATE('t', 'N_rN', 'M', 1)\n";
    auto a = runModel(model, asc);
    auto d = runModel(model, desc);
    REQUIRE(a.solveOk);
    REQUIRE(d.solveOk);
    for (const char* v : {"y1", "y2", "y3"}) CHECK_THAT(d[v], WithinAbs(a[v], 1e-9));
}

TEST_CASE("CS-BUG-INTERP-DESC: the derivative of a descending table is right (Newton on the x column)", "[library-gaps][interpolate]") {
    // Inverse lookup: find N_rN such that M = 395 -> 0.925. Needs d M / d N_rN from the table.
    const std::map<std::string, std::string> tables = {{"t", "N_rN,M\n1,454\n0.95,420\n0.9,370\n"}};
    auto run = runModel("INTERPOLATE('t', 'N_rN', 'M', N) = 395\n", tables);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["N"], WithinAbs(0.925, 1e-6));
}

// ============================================================================
// CS-BUG-LOOKUP-PATH: companion tables with a bare model file name
// ============================================================================

TEST_CASE("CS-BUG-LOOKUP-PATH: companion tables are loaded for a bare file name", "[library-gaps][lookup-path]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_lookup_path";
    fs::create_directories(dir);
    {
        std::ofstream f(dir / "m-data.csv");
        f << "x,y\n1,10\n2,20\n";
    }
    const fs::path previous = fs::current_path();
    fs::current_path(dir);
    // as in `cd examples && coolsolve lookup_demo.eescode`: no directory part
    auto store = loadLookupTableForModel("m.eescode");
    fs::current_path(previous);
    fs::remove_all(dir);
    REQUIRE(store.has("data"));
}

// ============================================================================
// CS-GAP-IF5: the 5-argument IF(A, B, X, Y, Z) of EES
// ============================================================================

TEST_CASE("CS-GAP-IF5: reproducer of the register and the three branches", "[library-gaps][if5]") {
    // EES: IF(A,B,X,Y,Z) = X if A<B, Y if A=B, Z if A>B. Register reproducer: 18000 s is past
    // the 3600 s switch, so the schedule gives its last value (0 W).
    auto run = runModel(R"(
tau = 18000
Q = if(tau, 3600, 10000, 10000, 0)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["Q"], WithinAbs(0.0, 1e-9));

    // distinct values for X, Y and Z tell the branches apart
    auto below = runModel("a = 1000\ny = if(a, 3600, 11, 22, 33)\n");
    auto equal = runModel("a = 3600\ny = if(a, 3600, 11, 22, 33)\n");
    auto above = runModel("a = 7200\ny = if(a, 3600, 11, 22, 33)\n");
    REQUIRE(below.solveOk);
    REQUIRE(equal.solveOk);
    REQUIRE(above.solveOk);
    CHECK_THAT(below["y"], WithinAbs(11.0, 1e-12));   // A < B -> X
    CHECK_THAT(equal["y"], WithinAbs(22.0, 1e-12));   // A = B -> Y
    CHECK_THAT(above["y"], WithinAbs(33.0, 1e-12));   // A > B -> Z
    CHECK(below.verified);
    CHECK(equal.verified);
    CHECK(above.verified);
}

TEST_CASE("CS-GAP-IF5: case-insensitive name, expressions as arguments", "[library-gaps][if5]") {
    auto run = runModel(R"(
a = 2
b = 3
y1 = IF(a, b, a + b, 0, -1)
y2 = If(b*2, a*3, 100, 200, 300*a)
y3 = iF(a - b, -1, 1, 2, 3) + 1
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["y1"], WithinAbs(5.0, 1e-12));     // 2 < 3
    CHECK_THAT(run["y2"], WithinAbs(200.0, 1e-12));   // 6 = 6
    CHECK_THAT(run["y3"], WithinAbs(3.0, 1e-12));     // -1 = -1 -> 2, plus 1
}

TEST_CASE("CS-GAP-IF5: nested IF5 and products of IF5 (vertical_ghes_refsim)", "[library-gaps][if5]") {
    // the penalty function of CSL-0081: IF(n_b,4,0,0,IF(n_b,144,1,0,0)) * ...
    const std::string model = R"(
n_b = %N%
lnt = %L%
pen_n = IF(n_b,4,0,0,IF(n_b,144,1,0,0))
pen_t = IF(lnt,-2,0,0,IF(lnt,3,1,0,0))
f_penalty = IF(n_b,4,0,0,IF(n_b,144,1,0,0))*IF(lnt,-2,0,0,IF(lnt,3,1,0,0))
)";
    auto make = [&](const char* n, const char* l) {
        std::string m = model;
        m.replace(m.find("%N%"), 3, n);
        m.replace(m.find("%L%"), 3, l);
        return runModel(m);
    };
    auto inside = make("36", "1");
    auto belowN = make("2", "1");
    auto aboveN = make("200", "1");
    auto belowT = make("36", "-5");
    REQUIRE(inside.solveOk);
    REQUIRE(belowN.solveOk);
    REQUIRE(aboveN.solveOk);
    REQUIRE(belowT.solveOk);
    CHECK(inside.verified);
    CHECK_THAT(inside["f_penalty"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(inside["pen_n"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(belowN["f_penalty"], WithinAbs(0.0, 1e-12));
    CHECK_THAT(aboveN["pen_n"], WithinAbs(0.0, 1e-12));
    CHECK_THAT(belowT["pen_t"], WithinAbs(0.0, 1e-12));
    CHECK_THAT(belowT["f_penalty"], WithinAbs(0.0, 1e-12));

    // nesting in any of the five positions (here in B and in X)
    auto deep = runModel("a = 5\ny = if(a, if(a, 3, 10, 10, 20), if(a, 9, 1, 2, 3), 7, 8)\n");
    REQUIRE(deep.solveOk);
    // inner B = if(5,3,..) = 20 (5 > 3 -> Z); a=5 < 20 -> X = if(5,9,1,2,3) = 1
    CHECK_THAT(deep["y"], WithinAbs(1.0, 1e-12));
}

TEST_CASE("CS-GAP-IF5: inside FUNCTION and PROCEDURE bodies", "[library-gaps][if5]") {
    auto run = runModel(R"(
FUNCTION sched(t)
  sched = if(t, 3600, 10000, 5000, 0)
END
PROCEDURE pwr(t : P, flag)
  P = IF(t, 3600, 10000, 5000, 0)
  flag = IF(t, 3600, IF(t, 1800, -2, -1, 0), 1, 2)
END
Q1 = sched(1000)
Q2 = sched(3600)
Q3 = sched(7200)
CALL pwr(7200 : P3, f3)
CALL pwr(1000 : P1, f1)
CALL pwr(1800 : P2, f2)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["Q1"], WithinAbs(10000.0, 1e-9));
    CHECK_THAT(run["Q2"], WithinAbs(5000.0, 1e-9));
    CHECK_THAT(run["Q3"], WithinAbs(0.0, 1e-9));
    CHECK_THAT(run["P3"], WithinAbs(0.0, 1e-9));
    CHECK_THAT(run["f3"], WithinAbs(2.0, 1e-9));     // t > 3600 -> 2
    CHECK_THAT(run["P1"], WithinAbs(10000.0, 1e-9));
    CHECK_THAT(run["f1"], WithinAbs(-2.0, 1e-9));    // t < 1800 -> -2
    CHECK_THAT(run["f2"], WithinAbs(-1.0, 1e-9));    // t = 1800 -> -1
}

TEST_CASE("CS-GAP-IF5: the derivative is that of the selected branch (AD)", "[library-gaps][if5]") {
    // x0 = 2, x1 = 3 (independent variables), n = 3 variables: a, b, c
    auto a = ADValue::independent(2.0, 0, 3);
    auto b = ADValue::independent(3.0, 1, 3);
    auto c = ADValue::independent(5.0, 2, 3);
    auto x = a * c;      // 10, d/da = 5, d/dc = 2
    auto y = b + c;      // 8,  d/db = 1, d/dc = 1
    auto z = c * c;      // 25, d/dc = 10

    SECTION("A < B selects X; A and B contribute nothing") {
        auto r = evaluateStandardFunction("if", {a, b, x, y, z});
        CHECK(r.value == 10.0);
        REQUIRE(r.gradient.size() == 3);
        CHECK(r.gradient[0] == 5.0);
        CHECK(r.gradient[1] == 0.0);
        CHECK(r.gradient[2] == 2.0);
    }
    SECTION("A = B selects Y") {
        auto r = evaluateStandardFunction("IF", {a, a, x, y, z});
        CHECK(r.value == 8.0);
        CHECK(r.gradient[0] == 0.0);
        CHECK(r.gradient[1] == 1.0);
        CHECK(r.gradient[2] == 1.0);
    }
    SECTION("A > B selects Z") {
        auto r = evaluateStandardFunction("if", {b, a, x, y, z});
        CHECK(r.value == 25.0);
        CHECK(r.gradient[0] == 0.0);
        CHECK(r.gradient[1] == 0.0);
        CHECK(r.gradient[2] == 10.0);
    }
    SECTION("the equality test is exact, as in EES") {
        auto s = ADValue::constant(0.1 + 0.2, 0);    // 0.30000000000000004
        auto t = ADValue::constant(0.3, 0);
        auto k = [](double v) { return ADValue::constant(v, 0); };
        CHECK(evaluateStandardFunction("if", {s, t, k(1), k(2), k(3)}).value == 3.0);
        CHECK(evaluateStandardFunction("if", {t, t, k(1), k(2), k(3)}).value == 2.0);
    }
    SECTION("the 3-argument form is unchanged (condition > 0)") {
        auto k = [](double v) { return ADValue::constant(v, 3); };
        CHECK(evaluateStandardFunction("if", {k(1), x, y}).value == 10.0);
        CHECK(evaluateStandardFunction("if", {k(0), x, y}).value == 8.0);
        CHECK(evaluateStandardFunction("if", {k(-1), x, y}).value == 8.0);
        CHECK(evaluateStandardFunction("if", {k(1), x, y}).gradient[0] == 5.0);
    }
}

TEST_CASE("CS-GAP-IF5: unknowns selected inside a coupled block (Jacobian)", "[library-gaps][if5]") {
    // Q and Q_wet form a 2x2 block: Q = IF(...) picks Q_wet (X) or Q_dry (Z), Q_wet depends on Q.
    // Newton needs dQ/dQ_wet = 1 on the X branch; a zero derivative would make the block singular.
    const std::string model = R"(
Q_dry = %QD%
Q = if(Q_dry, Q_wet, Q_wet, Q_dry, Q_dry)
Q_wet = 30 + 0.5*Q
)";
    auto make = [&](const char* qd) {
        std::string m = model;
        m.replace(m.find("%QD%"), 4, qd);
        return runModel(m);
    };
    // Q_dry = 40: X branch, Q = 30 + 0.5 Q -> 60 (Q_wet = 60 > 40, consistent)
    auto wet = make("40");
    REQUIRE(wet.solveOk);
    CHECK(wet.verified);
    CHECK_THAT(wet["Q"], WithinAbs(60.0, 1e-6));
    CHECK_THAT(wet["Q_wet"], WithinAbs(60.0, 1e-6));
    // Q_dry = 80: Z branch, Q = 80 (Q_wet = 70 < 80, consistent)
    auto dry = make("80");
    REQUIRE(dry.solveOk);
    CHECK(dry.verified);
    CHECK_THAT(dry["Q"], WithinAbs(80.0, 1e-6));
    CHECK_THAT(dry["Q_wet"], WithinAbs(70.0, 1e-6));
}

TEST_CASE("CS-GAP-IF5: the arguments A and B are structural incidences", "[library-gaps][if5]") {
    EESParser parser;
    auto parsed = parser.parse("y = if(a, b, x, z, w)\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    REQUIRE(ir.getEquations().size() == 1);
    const auto& vars = ir.getEquations()[0].variables;
    for (const char* v : {"y", "a", "b", "x", "z", "w"}) CHECK(vars.count(v) == 1);

    // B is the only place where the unknown Q_wet appears: the block must still use it
    auto run = runModel(R"(
Q_dry = 40
Q_wet = 30 + 0.25*Q
Q = if(Q_dry, Q_wet, 1, 2, 3)*10
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // Q_wet = 37.5 < 40 -> A > B -> Z = 3 -> Q = 30 -> Q_wet = 37.5
    CHECK_THAT(run["Q"], WithinAbs(30.0, 1e-6));
    CHECK_THAT(run["Q_wet"], WithinAbs(37.5, 1e-6));
}

TEST_CASE("CS-GAP-IF5: unit annotations on the literal arguments", "[library-gaps][if5]") {
    // EES style from the cooling-coil / DHW models: the literals carry their unit in brackets
    auto run = runModel(R"(
tau = 18000 [s]
Q = if(tau, 3600[s], 10000[W], 10000[W], 0[W])
Q2 = if(1800[s], 3600[s], 10000[W], 10000[W], 0[W])
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["Q"], WithinAbs(0.0, 1e-9));
    CHECK_THAT(run["Q2"], WithinAbs(10000.0, 1e-9));
}

TEST_CASE("CS-GAP-IF5: the 3-argument IF still works", "[library-gaps][if5]") {
    auto run = runModel(R"(
x = 3
y1 = if(x - 2, 10, 20)
y2 = IF(x - 5, 10, 20)
y3 = if(x - 3, 10, 20)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["y1"], WithinAbs(10.0, 1e-12));   // condition > 0
    CHECK_THAT(run["y2"], WithinAbs(20.0, 1e-12));
    CHECK_THAT(run["y3"], WithinAbs(20.0, 1e-12));   // 0 is not > 0
}

TEST_CASE("CS-GAP-IF5: IF with another number of arguments is rejected with the accepted forms", "[library-gaps][if5]") {
    for (const char* call : {"if(1, 2)", "if(1, 2, 3, 4)", "if(1, 2, 3, 4, 5, 6)", "IF(1)"}) {
        INFO(call);
        auto run = runModel(std::string("y = ") + call + "\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "IF takes 3 or 5 arguments"));
        CHECK(anyContains({run.message}, "IF(cond, a, b)"));
        CHECK(anyContains({run.message}, "IF(A, B, X, Y, Z)"));
    }
    // the number of arguments received is reported
    auto four = runModel("y = if(1, 2, 3, 4)\n");
    CHECK(anyContains({four.message}, "got 4"));
}

TEST_CASE("CS-GAP-IF5: LaTeX rendering is that of the 3-argument form", "[library-gaps][if5]") {
    EESParser parser;
    auto parsed = parser.parse("y = if(tau, 3600, Q_a, Q_b, 0)\nz = if(tau, Q_a, Q_b)\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    const std::string five = coolsolve::expressionToLatex(ir.getEquations()[0].rhs);
    const std::string three = coolsolve::expressionToLatex(ir.getEquations()[1].rhs);
    CHECK(five.find("\\mathrm{if}") == 0);
    CHECK(three.find("\\mathrm{if}") == 0);
    // all five arguments are rendered, in order
    CHECK(five.find("\\tau") != std::string::npos);
    CHECK(five.find("3600") != std::string::npos);
    CHECK(five.find("Q_{a}") != std::string::npos);
    CHECK(five.find("Q_{b}") != std::string::npos);
    size_t commas = 0;
    for (char ch : five) if (ch == ',') ++commas;
    CHECK(commas == 4);
}

TEST_CASE("CS-GAP-IF5: in a FUNCTION called with an unknown (Jacobian through a user function)", "[library-gaps][if5]") {
    // f(x) = x^2 below 2, 3 x above: the derivative of the selected branch drives Newton
    const std::string model = R"(
FUNCTION f(x)
  f = if(x, 2, x^2, 4, 3*x)
END
f(y) = %RHS%
)";
    auto make = [&](const char* rhs) {
        std::string m = model;
        m.replace(m.find("%RHS%"), 5, rhs);
        return runModel(m);
    };
    auto lowBranch = make("2.25");
    REQUIRE(lowBranch.solveOk);
    CHECK_THAT(lowBranch["y"], WithinAbs(1.5, 1e-6));
    auto highBranch = make("9");
    REQUIRE(highBranch.solveOk);
    CHECK_THAT(highBranch["y"], WithinAbs(3.0, 1e-6));
}

TEST_CASE("CS-GAP-IF5: in an INTEGRAL (dynamic) model", "[library-gaps][if5]") {
    // schedule dy/dt = 3 until t = 2.2, then 0: y(4) = 6.6 (fixed step 0.01 -> error < 0.05)
    SECTION("through a variable, evaluated by the algebraic solver at each step") {
        auto run = runModel(R"(
y = 0 + INTEGRAL(dydt, t, 0, 4, 0.01)
dydt = if(t, 2.2, 3, 3, 0)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(6.6, 0.05));
    }
    SECTION("directly as the integrand expression") {
        auto run = runModel(R"(
y = 0 + INTEGRAL(if(t, 2.2, 3, 3, 0), t, 0, 4, 0.01)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(6.6, 0.05));
    }
    SECTION("other built-in functions in the integrand expression are evaluated as well (were 0)") {
        auto run = runModel(R"(
y = 0 + INTEGRAL(sqrt(t) + max(t, 1), t, 0, 4, 0.01)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(16.0 / 3.0 + 8.5, 0.05));
    }
    SECTION("in the initial value of the state") {
        auto run = runModel(R"(
y = if(1, 2, 5, 6, 7) + INTEGRAL(dydt, t, 0, 4, 0.01)
dydt = 1
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(5.0 + 4.0, 0.05));   // 1 < 2 -> 5
    }
}


// ============================================================================
// CS-BUG-ARG-DROP: an argument that cannot be parsed is an error, not dropped
// ============================================================================

TEST_CASE("CS-BUG-ARG-DROP: reproducer of the register is a parse error", "[library-gaps][arg-drop]") {
    // was SUCCESS with y = 2: the argument '3 4' was skipped and max(2) evaluated
    auto run = runModel("y = max(3 4, 2)\n");
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(anyContains({run.message}, "Could not parse argument 1 of 'max': '3 4'"));
    CHECK(parseErrors("y = max(3 4, 2)\n") == std::vector<std::string>{"Could not parse argument 1 of 'max': '3 4'"});
}

TEST_CASE("CS-BUG-ARG-DROP: the message gives the position and the text of the argument", "[library-gaps][arg-drop]") {
    using Messages = std::vector<std::string>;
    CHECK(parseErrors("y = max(2, 3 4)\n") == Messages{"Could not parse argument 2 of 'max': '3 4'"});
    CHECK(parseErrors("y = max(1, 2, 3 4, 5)\n") == Messages{"Could not parse argument 3 of 'max': '3 4'"});
    // a user function; a thermophysical function (the position counts the named arguments too)
    CHECK(parseErrors("y = f(1 2, 3)\n") == Messages{"Could not parse argument 1 of 'f': '1 2'"});
    CHECK(parseErrors("h = enthalpy(Water, T=20, 3 4)\n") ==
          Messages{"Could not parse argument 3 of 'enthalpy': '3 4'"});

    // one error, reported on the line of the call; the other lines are not concerned
    EESParser parser;
    auto result = parser.parse("a = 1\nb = 2\ny = max(3 4, 2)\nc = 3\n");
    CHECK_FALSE(result.success);
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].line == 3);
}

TEST_CASE("CS-BUG-ARG-DROP: nested calls report the innermost faulty call", "[library-gaps][arg-drop]") {
    using Messages = std::vector<std::string>;
    CHECK(parseErrors("y = 1 + sqrt(max(3 4, 2))\n") == Messages{"Could not parse argument 1 of 'max': '3 4'"});
    CHECK(parseErrors("y = sqrt(2 * max(1, 2 3)) + 1\n") == Messages{"Could not parse argument 2 of 'max': '2 3'"});
    CHECK(parseErrors("y = x * (2 + abs(3 4))\n") == Messages{"Could not parse argument 1 of 'abs': '3 4'"});
    // on the left-hand side
    CHECK(parseErrors("max(3 4, 2) = y\n") == Messages{"Could not parse argument 1 of 'max': '3 4'"});
    // several statements on one line: only the faulty one is refused
    CHECK(parseErrors("a = 1 ; y = max(3 4, 2) ; b = 2\n") ==
          Messages{"Could not parse argument 1 of 'max': '3 4'"});
}

TEST_CASE("CS-BUG-ARG-DROP: an empty argument is an error", "[library-gaps][arg-drop]") {
    using Messages = std::vector<std::string>;
    CHECK(parseErrors("y = max(3, )\n") == Messages{"Could not parse argument 2 of 'max': ''"});
    CHECK(parseErrors("y = max(, 3)\n") == Messages{"Could not parse argument 1 of 'max': ''"});
    CHECK(parseErrors("y = max(1,, 3)\n") == Messages{"Could not parse argument 2 of 'max': ''"});
}

TEST_CASE("CS-BUG-ARG-DROP: a call without argument is still valid", "[library-gaps][arg-drop]") {
    // `FUNCTION g( )` and `PROCEDURE z( : r)` used to declare one parameter with an empty name
    auto run = runModel(R"(
FUNCTION g( )
  g = 7
END
FUNCTION h()
  h = 8
END
PROCEDURE z( : r)
  r = 5
END
CALL z( : q)
x = pi()
y = pi( )
z = g()
w = g( ) + h()
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(3.14159265358979, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(3.14159265358979, 1e-12));
    CHECK_THAT(run["z"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["w"], WithinAbs(15.0, 1e-12));
    CHECK_THAT(run["q"], WithinAbs(5.0, 1e-12));
}

TEST_CASE("CS-BUG-ARG-DROP: the same defect in an array index", "[library-gaps][arg-drop]") {
    using Messages = std::vector<std::string>;
    // was x[3 4] -> the plain variable x (an index dropped without any message)
    CHECK(parseErrors("y = x[3 4]\n") == Messages{"Could not parse index 1 of 'x': '3 4'"});
    CHECK(parseErrors("y = A[1, 2 3]\n") == Messages{"Could not parse index 2 of 'A': '2 3'"});
    CHECK(parseErrors("x[1 2] = 5\n") == Messages{"Could not parse index 1 of 'x': '1 2'"});
    CHECK(parseErrors("y = x[1, ]\n") == Messages{"Could not parse index 2 of 'x': ''"});
    // an index range is an argument list (CS-GAP-ARRAY-ARG), its bounds are parsed like the other indices
    CHECK(parseErrors("y = sum(x[1..3])\n").empty());
    CHECK(parseErrors("y = sum(x[1..3 4])\n") == Messages{"Could not parse index 1 of 'x': '1..3 4'"});
    // an index inside a call argument: one message, the innermost
    CHECK(parseErrors("y = max(x[3 4], 2)\n") == Messages{"Could not parse index 1 of 'x': '3 4'"});
}

TEST_CASE("CS-BUG-ARG-DROP: well-formed arguments and indices of every kind still work", "[library-gaps][arg-drop]") {
    auto run = runModel(R"(
DUPLICATE i=1,3
  x[i] = 10*i
END
A[1,2] = 5
fl$ = 'Water'
y1 = x[2] + A[1, 2] + max(x[1], x[3])
h1 = enthalpy(Water, T=20, P=1E5)
h2 = enthalpy('Water', T=20, P=1E5)
h3 = enthalpy(fl$, T=20, P=1E5)
n1 = max(min(3, 4), abs(-5))
n2 = max(1, 2) {a comment}
n3 = max(1, 3) "another comment"
n4 = max((1+2)*3, 2^3-1)
n5 = max(-1, -2)
n6 = max(1E-3, 1.5e+2)
n7 = max(3600[s], 2)
n8 = max( 1 , 2 )
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["y1"], WithinAbs(20.0 + 5.0 + 30.0, 1e-9));
    CHECK_THAT(run["h1"], WithinRel(run["h2"], 1e-12));
    CHECK_THAT(run["h3"], WithinRel(run["h2"], 1e-12));
    CHECK_THAT(run["n1"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["n2"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["n3"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["n4"], WithinAbs(9.0, 1e-12));
    CHECK_THAT(run["n5"], WithinAbs(-1.0, 1e-12));
    CHECK_THAT(run["n6"], WithinAbs(150.0, 1e-12));
    CHECK_THAT(run["n7"], WithinAbs(3600.0, 1e-12));
    CHECK_THAT(run["n8"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-BUG-ARG-DROP: statements in FUNCTION/PROCEDURE bodies are checked too", "[library-gaps][arg-drop]") {
    using Messages = std::vector<std::string>;
    // was: the statement was dropped without any message (y = 1 !)
    CHECK(parseErrors("FUNCTION f(a)\n  f = max(3 4, a)\nEND\ny = f(1)\n") ==
          Messages{"Could not parse argument 1 of 'max': '3 4'"});
    CHECK(parseErrors("PROCEDURE p(a : r)\n  IF (a > 0) THEN r = max(3 4, a) ELSE r = 0\nEND\nCALL p(1 : y)\n") ==
          Messages{"Could not parse argument 1 of 'max': '3 4'"});
    CHECK(parseErrors("PROCEDURE p(a : r)\n  r = x[3 4] ; r = a\nEND\nCALL p(1 : y)\n") ==
          Messages{"Could not parse index 1 of 'x': '3 4'"});
    // an IF condition names the faulty argument
    CHECK(parseErrors("PROCEDURE p(a : r)\n  IF (max(3 4, a) > 0) THEN r = 1 ELSE r = 0\nEND\nCALL p(1 : y)\n") ==
          Messages{"Could not parse IF condition '(max(3 4, a) > 0)': Could not parse argument 1 of 'max': '3 4'"});
}

TEST_CASE("CS-BUG-ARG-DROP: the same defect in the arguments of a CALL", "[library-gaps][arg-drop]") {
    using Messages = std::vector<std::string>;
    const std::string proc = "PROCEDURE p(a, b : c)\n  c = a + b\nEND\n";
    // was: the argument was dropped and the call refused at solve time (or worse: not at all)
    CHECK(parseErrors(proc + "CALL p(3 4, 2 : r)\n") == Messages{"Could not parse argument 1 of 'p': '3 4'"});
    CHECK(parseErrors(proc + "CALL p(1, 2 3 : r)\n") == Messages{"Could not parse argument 2 of 'p': '2 3'"});
    CHECK(parseErrors(proc + "CALL p(max(3 4, 2), 1 : r)\n") ==
          Messages{"Could not parse argument 1 of 'max': '3 4'"});
    // an output may be an expression in the main program (CS-GAP-CALL-EXPR-OUT: EES equates the output to
    // it; it was dropped silently, then refused with "must be a variable"); one that cannot be parsed is not
    CHECK(parseErrors(proc + "CALL p(1, 2 : a+b)\n").empty());
    CHECK(parseErrors(proc + "CALL p(1, 2 : r, z/1000)\n").empty());
    CHECK(parseErrors(proc + "CALL p(1, 2 : r, 3 4)\n") == Messages{"Could not parse output 2 of 'p': '3 4'"});
    // in a procedure body
    CHECK(parseErrors(proc + "PROCEDURE q(a : r)\n  CALL p(3 4, a : r)\nEND\nCALL q(1 : y)\n") ==
          Messages{"Could not parse argument 1 of 'p': '3 4'"});

    // well-formed calls: expressions, nested calls, array elements as inputs and outputs, strings
    auto run = runModel(proc + R"(
x[1] = 4
CALL p(max(1, 2), x[1] + 1 : r)
CALL p(1, 2 : y[2])
CALL p(unitsystem('SI'), 2 : s)
CALL p( 1 , 2 : t )
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["r"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["y[2]"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["s"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["t"], WithinAbs(3.0, 1e-12));
}

// ============================================================================
// CS-BUG-ARG-DROP (named-argument part, register section 8, C-147): `name=value`
// is a named argument only for the functions that take them
// ============================================================================

TEST_CASE("CS-BUG-ARG-DROP: IF(a=5, 1, 2) is a clear error, not 'Unknown fluid'", "[library-gaps][named-arg]") {
    // was: parsed as a property call with a named argument 'a', solved with
    // 'Unknown fluid: ""' (C-147). EES has no comparison inside the IF function.
    auto run = runModel("a = 5\nt = IF(a=5, 1, 2)\n");
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(anyContains({run.message}, "IF does not take named arguments ('a=5')"));
    CHECK(anyContains({run.message}, "IF(A, B, X, Y, Z)"));
    CHECK_FALSE(anyContains({run.message}, "Unknown fluid"));

    // the text of the argument is reported as written, with the name of the function as written
    CHECK(anyContains(parseErrors("t = if(a = 5, 1, 2)\n"), "if does not take named arguments ('a = 5')"));
    // the 5-argument form too
    CHECK(anyContains(parseErrors("t = IF(a=5, 1, 2, 3, 4)\n"), "IF does not take named arguments ('a=5')"));
    // inside a FUNCTION body (the statement was dropped without any message)
    CHECK(anyContains(parseErrors("FUNCTION f(a)\n  f = IF(a=5, 1, 2)\nEND\ny = f(5)\n"),
                      "IF does not take named arguments ('a=5')"));
}

TEST_CASE("CS-BUG-ARG-DROP: the forms that do the test are accepted", "[library-gaps][named-arg]") {
    // what the error message recommends: the 5-argument form tests A = B exactly
    auto run = runModel(R"(
a = 5
t_eq = IF(a, 5, 2, 1, 2)
t_lt = IF(a - 6, 1, 2)
t_gt = IF(a - 4, 1, 2)
t_cmp = IF(a < 6, 1, 2)
t_cmp2 = IF(a >= 6, 1, 2)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["t_eq"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["t_lt"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["t_gt"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["t_cmp"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["t_cmp2"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-BUG-ARG-DROP: the other built-in functions take no named argument either", "[library-gaps][named-arg]") {
    CHECK(anyContains(parseErrors("y = max(a=1, 2)\n"), "max does not take named arguments ('a=1')"));
    CHECK(anyContains(parseErrors("y = sin(x=30)\n"), "sin does not take named arguments ('x=30')"));
    CHECK(anyContains(parseErrors("y = 2 + sqrt(abs(a=1))\n"), "abs does not take named arguments ('a=1')"));
    // the message says which functions do take them
    CHECK(anyContains(parseErrors("y = max(a=1, 2)\n"), "thermophysical"));
    // the error is the only one of the line
    CHECK(parseErrors("y = max(a=1, 2)\n").size() == 1);
}

TEST_CASE("CS-BUG-ARG-DROP: the functions that take named arguments are untouched", "[library-gaps][named-arg]") {
    // thermophysical property functions, also humid air, string-valued fluid, upper case
    auto run = runModel(R"(
fl$ = 'R134a'
h1 = enthalpy(Water, T=20, P=1E5)
h2 = ENTHALPY(fl$, T=20, x=1)
w = humrat(AirH2O, T=25, P=101325, R=0.5)
u = intenergy(Water, T=20, P=1E5)
Ts = T_sat(Water, P=1E5)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["h1"], WithinAbs(CoolProp::PropsSI("H", "T", 293.15, "P", 1e5, "Water"), 1e-3));
    CHECK_THAT(run["Ts"], WithinAbs(99.6, 0.1));
    CHECK(run["w"] > 0.0);

    // the indexed aggregation and the EES INTERPOLATE form carry `name=value` arguments: not
    // refused by the parser (the evaluation of the first is CS-GAP-SUM-INDEXED, that of the second
    // CS-GAP-INTERP-EES)
    CHECK(parseErrors("s = sum(x[i], i=1, 4)\n").empty());
    CHECK(parseErrors("y = interpolate('t', 'x', 'y', x=3)\n").empty());
}

// ============================================================================
// CS-GAP-UNIT-SUBEXPR / CS-GAP-UNIT-NAMEDARG: a number followed by [unit] inside an expression
// ============================================================================
// The annotation is documentation only: the number is used as written (no unit conversion).

TEST_CASE("CS-GAP-UNIT-SUBEXPR: reproducer of the register", "[library-gaps][unit-literal]") {
    // ULiege style: each term of an absolute temperature carries its unit
    auto run = runModel(R"(
T_A = 25.0[K]+273.15 [K]
v_l = 1
x = T_A - 273.15
y = 2*x
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["T_A"], WithinAbs(298.15, 1e-12));
    CHECK_THAT(run["x"], WithinAbs(25.0, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(50.0, 1e-12));
}

TEST_CASE("CS-GAP-UNIT-SUBEXPR: the annotation at the end of the right-hand side, as before", "[library-gaps][unit-literal]") {
    auto run = runModel(R"(
a = 25.0+273.15 [K]
b = (25.0+273.15) [K]
c = 25.0 [K] + 273.15
P_in = 1E5 [Pa]
d = 25 [C] + 3 [K] "[K]"
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(298.15, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(298.15, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(298.15, 1e-12));
    CHECK_THAT(run["P_in"], WithinAbs(1e5, 1e-9));
    CHECK_THAT(run["d"], WithinAbs(28.0, 1e-12));
}

TEST_CASE("CS-GAP-UNIT-SUBEXPR: annotation on each term, in parentheses, on products and powers", "[library-gaps][unit-literal]") {
    auto run = runModel(R"(
a = 1[m] + 2 [m] + 3[m]
b = (25[C] + 273.15 [K]) * 2
c = 2 * (3 [m] - 1[m])
d = 3[kg/m^3] * 2[m^3]
e = 4 [J/kg-K] + 1 [J/kg-K]
f = 2[m]^2
g = 10^2[m]
h = 12[m] / 4 [s]
i = max(3[m], 5 [m])
j = 7[m] + a
k = (25 + 273.15) [K] + 2 [K]
l = 2 * (3 + 4) [m]
m1 = ((1 + 2) [m] + 3) [m]
n = (25[C] + 273.15 [K]) [K]
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(596.3, 1e-9));
    CHECK_THAT(run["c"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["d"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["e"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["f"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["g"], WithinAbs(100.0, 1e-12));
    CHECK_THAT(run["h"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["i"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["j"], WithinAbs(13.0, 1e-12));
    CHECK_THAT(run["k"], WithinAbs(300.15, 1e-9));
    CHECK_THAT(run["l"], WithinAbs(14.0, 1e-12));
    CHECK_THAT(run["m1"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["n"], WithinAbs(298.15, 1e-9));
}

TEST_CASE("CS-GAP-UNIT-SUBEXPR: negative numbers and exponent notation with a space before the unit", "[library-gaps][unit-literal]") {
    auto run = runModel(R"(
a = -5[C]
b = 3 * -5[C]
c = 2 - 5[C]
d = -5 [C] + 10[C]
e = (-5[C])
p = 1E5 [Pa]
q = 1.5e-3 [m] * 2
r = 2E+3 [W] + 1[W]
s = 1E-3[m] + 2.5E+1 [m]
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(-5.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(-15.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(-3.0, 1e-12));
    CHECK_THAT(run["d"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["e"], WithinAbs(-5.0, 1e-12));
    CHECK_THAT(run["p"], WithinAbs(1e5, 1e-9));
    CHECK_THAT(run["q"], WithinAbs(3e-3, 1e-15));
    CHECK_THAT(run["r"], WithinAbs(2001.0, 1e-9));
    CHECK_THAT(run["s"], WithinAbs(25.001, 1e-12));
}

TEST_CASE("CS-GAP-UNIT-SUBEXPR: the unit is documentation only, the number is not converted", "[library-gaps][unit-literal]") {
    auto run = runModel(R"(
T_k = 25[K]
T_c = 25[C]
P_b = 1[bar]
m_g = 5[g]
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["T_k"], WithinAbs(25.0, 1e-12));
    CHECK_THAT(run["T_c"], WithinAbs(25.0, 1e-12));
    CHECK_THAT(run["P_b"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["m_g"], WithinAbs(5.0, 1e-12));
}

TEST_CASE("CS-GAP-UNIT-SUBEXPR: array elements are still indexed and malformed annotations are refused", "[library-gaps][unit-literal]") {
    auto run = runModel(R"(
x[1] = 4
x[2] = 6
y = x[2] + 3[m] + x[1]
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(13.0, 1e-12));

    CHECK_FALSE(parseErrors("y = 5[m\n").empty());        // unclosed
    CHECK_FALSE(parseErrors("y = 5[m]]\n").empty());       // extra closing bracket
    CHECK_FALSE(parseErrors("y = 5[m] 3\n").empty());      // text after the annotation
    CHECK_FALSE(parseErrors("y = 5[m] + 3 [\n").empty());
    // a parenthesis followed by something else than an annotation is not a call of a function without a name
    CHECK_FALSE(parseErrors("y = (1 + 2) 3\n").empty());
    CHECK_FALSE(parseErrors("y = (1 + 2)(3 + 4)\n").empty());
}

TEST_CASE("CS-GAP-UNIT-NAMEDARG: reproducer of the register", "[library-gaps][unit-literal]") {
    auto run = runModel("h = enthalpy(Water,P=1E5[Pa],T=20)\n");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    auto plain = runModel("h = enthalpy(Water,P=1E5,T=20)\n");
    REQUIRE(plain.solveOk);
    CHECK(run["h"] == plain["h"]);
    CHECK_THAT(run["h"], WithinAbs(CoolProp::PropsSI("H", "T", 293.15, "P", 1e5, "Water"), 1e-3));
}

TEST_CASE("CS-GAP-UNIT-NAMEDARG: variants (space before the unit, both arguments, negative number, expression)", "[library-gaps][unit-literal]") {
    auto run = runModel(R"(
h_ref = enthalpy(Water,P=1E5,T=20)
h1 = enthalpy(Water,P=1E5 [Pa],T=20 [C])
h2 = enthalpy(Water,T=20[C],P=1E5[Pa])
h3 = enthalpy(Water,P=1.0E+5 [Pa],T=15[C]+5 [K])
h4 = enthalpy('Water',P=(1E5[Pa]),T=20[C])
h5 = 2*enthalpy(Water,P=1E5[Pa],T=20 [C]) + 1[J/kg]
hn_ref = enthalpy(R134a,T=-5,P=2E5)
hn = enthalpy(R134a,T=-5[C],P=2E5[Pa])
w_ref = humrat(AirH2O,T=25,P=101325,R=0.5)
w = humrat(AirH2O,T=25[C],P=101325[Pa],R=0.5)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK(run["h1"] == run["h_ref"]);
    CHECK(run["h2"] == run["h_ref"]);
    CHECK(run["h3"] == run["h_ref"]);
    CHECK(run["h4"] == run["h_ref"]);
    CHECK_THAT(run["h5"], WithinAbs(2.0 * run["h_ref"] + 1.0, 1e-9));
    CHECK(run["hn"] == run["hn_ref"]);
    CHECK(run["w"] == run["w_ref"]);
}

TEST_CASE("CS-GAP-UNIT-NAMEDARG: a literal with its unit as a named argument of a coupled block (AD)", "[library-gaps][unit-literal]") {
    // the unknown temperature is found by Newton through the property call whose other argument carries a unit
    auto run = runModel(R"(
h_t = enthalpy(Water, P=1E5[Pa], T=T_x)
h_t = 1E5[J/kg]
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["T_x"], WithinAbs(CoolProp::PropsSI("T", "H", 1e5, "P", 1e5, "Water") - 273.15, 1e-6));
}

TEST_CASE("CS-GAP-UNIT-NAMEDARG: as in the ULiege files, the unit is not converted (1e2[kPa] is 100 Pa)", "[library-gaps][unit-literal]") {
    auto run = runModel(R"(
T_su_ref = 20
h_a = enthalpy(Water,P=1e2[kPa],T=T_su_ref)
h_pa = enthalpy(Water,P=100,T=T_su_ref)
h_kpa = enthalpy(Water,P=1e5,T=T_su_ref)
)");
    REQUIRE(run.solveOk);
    CHECK(run["h_a"] == run["h_pa"]);
    CHECK(run["h_a"] != run["h_kpa"]);
}

// ============================================================================
// CS-BUG-BODY-LINE-DROP: a statement that cannot be parsed in a FUNCTION/PROCEDURE
// body is an error (it was dropped without any message)
// ============================================================================

TEST_CASE("CS-BUG-BODY-LINE-DROP: reproducer of the register is a parse error", "[library-gaps][body-line-drop]") {
    // was SUCCESS with y = 1: the statement `f = 3 4` was gone
    const std::string code = "FUNCTION f(a)\nf = 3 4\nEND\ny = f(1)\n";
    auto run = runModel(code);
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(anyContains({run.message}, "Line 2: Could not parse line"));

    EESParser parser;
    auto result = parser.parse(code);
    CHECK_FALSE(result.success);
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].line == 2);   // the line number of the faulty statement
    CHECK(result.errors[0].context == "f = 3 4");
}

TEST_CASE("CS-BUG-BODY-LINE-DROP: same message as in the main program, for every kind of unparsable line", "[library-gaps][body-line-drop]") {
    using Messages = std::vector<std::string>;
    const Messages main = parseErrors("y = 3 4\n");
    REQUIRE(main == Messages{"Could not parse line"});

    // in a FUNCTION body and in a PROCEDURE body
    for (const std::string line : {"f = 3 4", "f", "f = ", "= 3", "3 4", "f 3", "f = a b + 1"}) {
        INFO("body line: " << line);
        CHECK(parseErrors("FUNCTION f(a)\n  " + line + "\nEND\ny = f(1)\n") == main);
        CHECK(parseErrors("PROCEDURE p(a : r)\n  " + line + "\nEND\nCALL p(1 : y)\n") == main);
    }

    // the line number is that of the statement, whatever its position in the body
    EESParser parser;
    auto result = parser.parse("FUNCTION f(a)\n  s = 1\n  t = 2\n  f = s + t 5\n  u = 3\nEND\ny = f(1)\nz = 2\n");
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].line == 4);
}

TEST_CASE("CS-BUG-BODY-LINE-DROP: faulty statements inside the blocks of a body", "[library-gaps][body-line-drop]") {
    using Messages = std::vector<std::string>;
    const Messages oneLine{"Could not parse line"};
    const std::string tail = "END\ny = f(1)\n";
    // single-line IF, both branches
    CHECK(parseErrors("FUNCTION f(a)\n  IF (a > 0) THEN f = 3 4\n" + tail) == oneLine);
    CHECK(parseErrors("FUNCTION f(a)\n  IF (a > 0) THEN f = 1 ELSE f = 3 4\n" + tail) == oneLine);
    // block IF with ELSE
    CHECK(parseErrors("FUNCTION f(a)\n  IF (a > 0) THEN\n    f = 3 4\n  ENDIF\n" + tail) == oneLine);
    CHECK(parseErrors("FUNCTION f(a)\n  IF (a > 0) THEN\n    f = 1\n  ELSE\n    f = 3 4\n  ENDIF\n" + tail) == oneLine);
    CHECK(parseErrors("FUNCTION f(a)\n  IF (a > 0) THEN\n    f = 1\n  ELSE f = 3 4\n  ENDIF\n" + tail) == oneLine);
    // DUPLICATE and REPEAT ... UNTIL
    CHECK(parseErrors("FUNCTION f(a)\n  DUPLICATE i = 1, 2\n    x[i] = 3 4\n  END\n  f = 1\n" + tail) == oneLine);
    CHECK(parseErrors("FUNCTION f(a)\n  k := 0\n  REPEAT\n    k := k 1\n  UNTIL (k > 3)\n  f = k\n" + tail) == oneLine);
    // several statements on one line: the faulty segment is named
    CHECK(parseErrors("FUNCTION f(a)\n  s = 1 ; t = 3 4 ; f = s\n" + tail) ==
          Messages{"Could not parse segment: t = 3 4"});
    CHECK(parseErrors("FUNCTION f(a)\n  IF (a > 0) THEN s = 1 ; t = 3 4\n  f = s\n" + tail) ==
          Messages{"Could not parse segment: t = 3 4"});
    // the other lines are not concerned: one error only
    EESParser parser;
    auto result = parser.parse("FUNCTION f(a)\n  f = a\n  IF (a > 0) THEN\n    f = 3 4\n  ENDIF\nEND\ny = f(1)\n");
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].line == 4);
}

TEST_CASE("CS-BUG-BODY-LINE-DROP: a specific message still wins over the generic one", "[library-gaps][body-line-drop]") {
    using Messages = std::vector<std::string>;
    CHECK(parseErrors("FUNCTION f(a)\n  f = max(3 4, a)\nEND\ny = f(1)\n") ==
          Messages{"Could not parse argument 1 of 'max': '3 4'"});
    CHECK(parseErrors("FUNCTION f(a)\n  f == 3\nEND\ny = f(1)\n") ==
          Messages{"Invalid operator '=='. EES uses '=' for equations, not '=='"});
}

TEST_CASE("CS-BUG-BODY-LINE-DROP: every legitimate form of a body line is still accepted", "[library-gaps][body-line-drop]") {
    // comments of all the EES forms (single-line and multi-line, braces, quotes, //), blank lines,
    // directives, statements separated by semicolons (with a trailing separator or comment),
    // assignments, CALL, single-line and block IF ... ELSE ... ENDIF, REPEAT ... UNTIL, DUPLICATE ... END
    const std::string code = R"(
PROCEDURE pp(a : r)
  r = a + 100
END

FUNCTION g(a)
{single-line brace comment}
{multi-line brace
 comment with x = 1 and a nested {brace}
 3 4 and f = (
}
"a quoted comment"
"a multi-line quoted comment
 that goes on, with 3 4 in it"
"
standalone-quote block comment
3 4
"
// a C-style comment 3 4
$Bookmark a directive in a body

  s = 0 ; // a comment after a separator
  t := 2 ; {brace after a separator} ;
  u = 3 "trailing quote comment"
  v = 4 {trailing brace comment}
  w = 5 // trailing C comment
  w2 = 6 ;
  ;
  z1 = 1 ; z2 = 2 ; z3 = 3 {last}
  CALL pp(1 : r)
  CALL pp(2 : r2) ; // procedure call followed by a comment
  IF (a > 0) THEN s = 1 ELSE s = -1
  IF (a > 5) THEN
    s = s + 10 // comment
  ELSE // comment
    s = s + 20 {comment}
  ENDIF "end of the test"
  IF (a > 0) THEN
    q = 1
  ELSE
    IF (a < 0) THEN
      q = 2
    ELSE
      q = 3
    ENDIF
  ENDIF
  IF (a > 0) THEN q = q + 1 ; q = q + 1
  DUPLICATE i = 1, 3
    x[i] = i
  END
  k := 0
  REPEAT
    k := k + 1
  UNTIL (k > 3) // loop end
  g = s + t + u + v + w + w2 + z1 + z2 + z3 + r + r2 + q + x[3] + k
END // end of the function
y = g(1)
)";
    // s = 1 then 21; t..w2 = 2+3+4+5+6; z = 6; r = 101, r2 = 102; q = 1 + 2; x[3] = 3; k = 4
    auto run = runModel(code);
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(21 + (2 + 3 + 4 + 5 + 6) + 6 + 101 + 102 + 3 + 3 + 4, 1e-12));
    CHECK(parseErrors(code).empty());
}

TEST_CASE("CS-BUG-BODY-LINE-DROP: legitimate bodies of PROCEDURE with outputs and single-line forms", "[library-gaps][body-line-drop]") {
    auto run = runModel(R"(
PROCEDURE p(a : r, s)
  "outputs"
  r = a ; s = 2*a ; // two statements
  IF (a > 100) THEN r = 0 ; s = 0
END
CALL p(3 : b, c)
FUNCTION h(a) // header comment
  h = a {the value}
END // end
d = h(4)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["b"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["d"], WithinAbs(4.0, 1e-12));
}

// ============================================================================
// CS-BUG-TRAILING-TEXT: text after the closing ) of a call or the ] of an array element
// is an error (it was ignored); a unit annotation [unit] is the only thing allowed there
// ============================================================================

TEST_CASE("CS-BUG-TRAILING-TEXT: reproducers of the register are parse errors", "[library-gaps][trailing-text]") {
    using Messages = std::vector<std::string>;
    // was SUCCESS with y = 2
    auto run = runModel("y = max(1,2) 3\n");
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(anyContains({run.message}, "Unexpected text '3' after the call 'max(1,2)'"));
    CHECK(parseErrors("y = max(1,2) x\n") == Messages{"Unexpected text 'x' after the call 'max(1,2)'"});
    // was y = x[1]
    CHECK(parseErrors("x[1]=5\ny = x[1] 2\n") == Messages{"Unexpected text '2' after the array element 'x[1]'"});
    auto run2 = runModel("x[1]=5\ny = x[1] 2\n");
    CHECK_FALSE(run2.parseOk);
    CHECK_FALSE(run2.solveOk);

    // one error, on the line of the statement
    EESParser parser;
    auto result = parser.parse("a = 1\nb = 2\ny = max(1,2) 3\nc = 3\n");
    CHECK_FALSE(result.success);
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].line == 3);
}

TEST_CASE("CS-BUG-TRAILING-TEXT: calls of every kind, on either side, nested, in arguments", "[library-gaps][trailing-text]") {
    using Messages = std::vector<std::string>;
    // user function, thermophysical function with named arguments, no argument, several texts
    CHECK(parseErrors("y = f(x) 3\n") == Messages{"Unexpected text '3' after the call 'f(x)'"});
    CHECK(parseErrors("h = enthalpy(Water, T=20, P=1E5) 3\n") ==
          Messages{"Unexpected text '3' after the call 'enthalpy(Water, T=20, P=1E5)'"});
    CHECK(parseErrors("y = pi() 3\n") == Messages{"Unexpected text '3' after the call 'pi()'"});
    CHECK(parseErrors("y = max(1,2) 3 z\n") == Messages{"Unexpected text '3 z' after the call 'max(1,2)'"});
    // two calls, the matching parenthesis is the one of the first call, not the last in the text
    CHECK(parseErrors("y = f(x)(z)\n") == Messages{"Unexpected text '(z)' after the call 'f(x)'"});
    CHECK(parseErrors("y = f(x) g(z)\n") == Messages{"Unexpected text 'g(z)' after the call 'f(x)'"});
    // an index written after a call
    CHECK(parseErrors("y = f(x)[1]\n") == Messages{"Unexpected text '[1]' after the call 'f(x)'"});
    // left-hand side
    CHECK(parseErrors("max(1,2) 3 = y\n") == Messages{"Unexpected text '3' after the call 'max(1,2)'"});
    // inside an expression, an argument, an index, a named argument
    CHECK(parseErrors("y = 1 + max(1,2) 3\n") == Messages{"Unexpected text '3' after the call 'max(1,2)'"});
    CHECK(parseErrors("y = max(1,2) 3 * 2\n") == Messages{"Unexpected text '3' after the call 'max(1,2)'"});
    CHECK(parseErrors("y = sqrt(max(1,2) 3)\n") == Messages{"Unexpected text '3' after the call 'max(1,2)'"});
    CHECK(parseErrors("y = x[max(1,2) 3]\n") == Messages{"Unexpected text '3' after the call 'max(1,2)'"});
    CHECK(parseErrors("h = enthalpy(Water, T=max(1,2) 3, P=1E5)\n") ==
          Messages{"Unexpected text '3' after the call 'max(1,2)'"});
    // several statements on a line: only the faulty one is refused
    CHECK(parseErrors("a = 1 ; y = max(1,2) 3 ; b = 2\n") == Messages{"Unexpected text '3' after the call 'max(1,2)'"});
    // the innermost faulty call is named
    CHECK(parseErrors("y = max(min(1,2) 3, 4)\n") == Messages{"Unexpected text '3' after the call 'min(1,2)'"});
    // strings with parentheses inside do not confuse the matching
    CHECK(parseErrors("y = f('a)', 'b(') 3\n") == Messages{"Unexpected text '3' after the call 'f('a)', 'b(')'"});
}

TEST_CASE("CS-BUG-TRAILING-TEXT: array elements and strings", "[library-gaps][trailing-text]") {
    using Messages = std::vector<std::string>;
    CHECK(parseErrors("y = x[1] 2\n") == Messages{"Unexpected text '2' after the array element 'x[1]'"});
    CHECK(parseErrors("y = A[1,2] z\n") == Messages{"Unexpected text 'z' after the array element 'A[1,2]'"});
    CHECK(parseErrors("y = x[i+1] 2\n") == Messages{"Unexpected text '2' after the array element 'x[i+1]'"});
    CHECK(parseErrors("y = x[y[1]] 2\n") == Messages{"Unexpected text '2' after the array element 'x[y[1]]'"});
    // an element indexed twice is not C: A[i,j] is the syntax (and [1] cannot be a unit)
    CHECK(parseErrors("y = A[1][2]\n") == Messages{"Unexpected text '[2]' after the array element 'A[1]'"});
    CHECK(parseErrors("y = x[1](2)\n") == Messages{"Unexpected text '(2)' after the array element 'x[1]'"});
    CHECK(parseErrors("x[1] 2 = 5\n") == Messages{"Unexpected text '2' after the array element 'x[1]'"});
    CHECK(parseErrors("y = 1 + x[1] 2\n") == Messages{"Unexpected text '2' after the array element 'x[1]'"});
    CHECK(parseErrors("y = max(x[1] 2, 3)\n") == Messages{"Unexpected text '2' after the array element 'x[1]'"});
    // a string literal (the closing quote ends it)
    CHECK(parseErrors("fl$ = 'Water' 3\n") == Messages{"Unexpected text '3' after the string 'Water'"});
    CHECK(parseErrors("fl$ = 'Water' x\n") == Messages{"Unexpected text 'x' after the string 'Water'"});
    CHECK(parseErrors("h = enthalpy('Water' 3, T=20, P=1E5)\n") == Messages{"Unexpected text '3' after the string 'Water'"});
}

TEST_CASE("CS-BUG-TRAILING-TEXT: the rule applies in FUNCTION/PROCEDURE bodies, IF conditions and CALL", "[library-gaps][trailing-text]") {
    using Messages = std::vector<std::string>;
    const Messages call{"Unexpected text '3' after the call 'max(1,2)'"};
    // was: the statement was kept with the call (or dropped)
    CHECK(parseErrors("FUNCTION f(a)\n  f = max(1,2) 3\nEND\ny = f(1)\n") == call);
    CHECK(parseErrors("PROCEDURE p(a : r)\n  IF (a > 0) THEN r = max(1,2) 3 ELSE r = 0\nEND\nCALL p(1 : y)\n") == call);
    CHECK(parseErrors("PROCEDURE p(a : r)\n  r = x[1] 2 ; r = a\nEND\nCALL p(1 : y)\n") ==
          Messages{"Unexpected text '2' after the array element 'x[1]'"});
    CHECK(parseErrors("PROCEDURE p(a : r)\n  IF (max(1,2) 3 > 0) THEN r = 1 ELSE r = 0\nEND\nCALL p(1 : y)\n") ==
          Messages{"Could not parse IF condition '(max(1,2) 3 > 0)': Unexpected text '3' after the call 'max(1,2)'"});

    // the arguments of a CALL
    const std::string proc = "PROCEDURE p(a, b : c)\n  c = a + b\nEND\n";
    CHECK(parseErrors(proc + "CALL p(max(1,2) 3, 1 : r)\n") == call);
    CHECK(parseErrors(proc + "CALL p(1, x[1] 2 : r)\n") == Messages{"Unexpected text '2' after the array element 'x[1]'"});
    CHECK(parseErrors(proc + "CALL p(1, 2 : y[1] 3)\n") == Messages{"Unexpected text '3' after the array element 'y[1]'"});
    // the text after the CALL itself (was ignored)
    CHECK(parseErrors(proc + "CALL p(1, 2 : r) 3\n") == Messages{"Unexpected text '3' after the arguments of 'p'"});
    CHECK(parseErrors(proc + "CALL p(1, 2 : r) x = 2\n") == Messages{"Unexpected text 'x = 2' after the arguments of 'p'"});
    CHECK(parseErrors(proc + "PROCEDURE q(a : r)\n  CALL p(a, 1 : r) 3\nEND\nCALL q(1 : y)\n") ==
          Messages{"Unexpected text '3' after the arguments of 'p'"});
    // ... but a separator, comments and blanks are not text
    auto run = runModel(proc + "CALL p(1, 2 : r);\nCALL p(3, 4 : s)  // a comment\nCALL p(5, 6 : t) \"a comment\" ;\n"
                        "PROCEDURE q(a : r)\n  CALL p(a, 1 : r) ; r = r + 1\nEND\nCALL q(1 : u)\n"
                        "CALL p(7, 8 : v) \"[kg/s]\"\nCALL p(1, 1 : w) \"[a]\" \"[b]\" ;\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["r"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["s"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["t"], WithinAbs(11.0, 1e-12));
    CHECK_THAT(run["u"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["v"], WithinAbs(15.0, 1e-12));
    CHECK_THAT(run["w"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-BUG-TRAILING-TEXT: a unit annotation after a call or an array element is accepted, at the end of the right-hand side", "[library-gaps][trailing-text]") {
    // this is the layout of the ULiege and library files: `Cmin = min(C_h, C_c)     [W/K]`
    auto run = runModel(R"(
x[1] = 5
x[2] = 7
A[1,2] = 11
FUNCTION f(a)
  f = 2*a
END
y1 = max(1,2) [kW]
y2 = x[1] [K]
y3 = A[1,2] [kg/m^3]
y4 = f(3) [J/kg-K]
y5 = enthalpy(Water, T=20, P=1E5) [J/kg]
y6 = max(1,2)[kW]
y7 = x[2]   [K]
y8 = f(x[1]) [W]
y9 = max(x[1], x[2]) [-]
y10 = x[ 2 ] [K]
y11 = -max(1, 2) [K]
y12 = sqrt(16) [m] "[m]"
y13 = sqrt(16) "[m]"
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["y1"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["y2"], WithinAbs(5.0, 1e-12));    // was x (without the index) before CS-BUG-ARG-DROP
    CHECK_THAT(run["y3"], WithinAbs(11.0, 1e-12));
    CHECK_THAT(run["y4"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["y5"], WithinAbs(CoolProp::PropsSI("H", "T", 293.15, "P", 1e5, "Water"), 1e-3));
    CHECK_THAT(run["y6"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["y7"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["y8"], WithinAbs(10.0, 1e-12));
    CHECK_THAT(run["y9"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["y10"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["y11"], WithinAbs(-2.0, 1e-12));
    CHECK_THAT(run["y12"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["y13"], WithinAbs(4.0, 1e-12));
}

TEST_CASE("CS-BUG-TRAILING-TEXT: a unit annotation after a call or an array element inside an expression", "[library-gaps][trailing-text]") {
    // the annotation ends the term it follows, as 3600[s] and (...) [unit] do
    auto run = runModel(R"(
x[1] = 5
x[2] = 7
FUNCTION f(a)
  f = 2*a
END
y1 = max(1,2) [kW] + 2
y2 = 2 + x[1] [K] + x[2] [K]
y3 = f(3) [kW] * 2 [-]
y4 = 10 - f(1) [W] / 2
y5 = f(2)[W]^2
y6 = max(x[1] [K], 2 [K]) + f(x[2] [K]) [W]
y7 = (max(1,2) [kW] + 1) [kW]
y8 = x[1] [K] > 3
y9 = if(max(1,2) [kW] - 1, 5, 6)
y10 = enthalpy(Water, T=f(10) [K], P=1E5 [Pa]) [J/kg] - enthalpy(Water, T=20, P=1E5)
y11 = -f(1) [W] + x[2] [K]
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["y1"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["y2"], WithinAbs(14.0, 1e-12));
    CHECK_THAT(run["y3"], WithinAbs(12.0, 1e-12));
    CHECK_THAT(run["y4"], WithinAbs(9.0, 1e-12));
    CHECK_THAT(run["y5"], WithinAbs(16.0, 1e-12));
    CHECK_THAT(run["y6"], WithinAbs(5.0 + 14.0, 1e-12));
    CHECK_THAT(run["y7"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["y8"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["y9"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["y10"], WithinAbs(0.0, 1e-6));
    CHECK_THAT(run["y11"], WithinAbs(5.0, 1e-12));
}

TEST_CASE("CS-BUG-TRAILING-TEXT: the annotation in FUNCTION/PROCEDURE bodies, CALL arguments and IF conditions", "[library-gaps][trailing-text]") {
    auto run = runModel(R"(
PROCEDURE p(a, b : c)
  c = a + min(b, 1E9) [-]
END
FUNCTION g(a)
  x[1] = 3
  s = max(a, 1) [K] + x[1] [K]
  IF (max(a, 1) [K] > 2) THEN s = s + 100 [K] ELSE s = s + 200
  g = s
END
x[1] = 5
CALL p(max(1,2) [kW], x[1] [K] : r)
y = g(4)
z = g(0)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["r"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(4.0 + 3.0 + 100.0, 1e-12));
    CHECK_THAT(run["z"], WithinAbs(1.0 + 3.0 + 200.0, 1e-12));
}

TEST_CASE("CS-BUG-TRAILING-TEXT: what follows the annotation, an empty or numeric annotation, an unclosed one", "[library-gaps][trailing-text]") {
    using Messages = std::vector<std::string>;
    // text after the annotation
    CHECK(parseErrors("y = max(1,2) [kW] 3\n") == Messages{"Unexpected text '[kW] 3' after the call 'max(1,2)'"});
    CHECK_FALSE(parseErrors("x[1] = 5\ny = x[1] [K] 3\n").empty());
    CHECK_FALSE(parseErrors("y = max(1,2) [kW] z\n").empty());
    // a number between brackets is an index, not a unit
    CHECK(parseErrors("y = max(1,2) [3]\n") == Messages{"Unexpected text '[3]' after the call 'max(1,2)'"});
    CHECK(parseErrors("y = A[1] [2]\n") == Messages{"Unexpected text '[2]' after the array element 'A[1]'"});
    // empty or unclosed
    CHECK_FALSE(parseErrors("y = max(1,2) []\n").empty());
    CHECK_FALSE(parseErrors("y = max(1,2) [kW\n").empty());
    CHECK_FALSE(parseErrors("x[1] = 5\ny = x[1] [K\n").empty());
    // an annotation does not make an unparsable argument or index valid
    CHECK(parseErrors("y = max(3 4, 2) [kW]\n") == Messages{"Could not parse argument 1 of 'max': '3 4'"});
    CHECK(parseErrors("y = x[3 4] [K]\n") == Messages{"Could not parse index 1 of 'x': '3 4'"});
    CHECK(parseErrors("y = max(1,2) [kW] + x[3 4]\n") == Messages{"Could not parse index 1 of 'x': '3 4'"});
}

TEST_CASE("CS-BUG-TRAILING-TEXT: calls and elements without trailing text are untouched", "[library-gaps][trailing-text]") {
    auto run = runModel(R"(
x[1] = 5
A[1,2] = 11
fl$ = 'Water'
n1 = max(1,2)
n2 = x[1]
n3 = A[1, 2] + max(x[1], 3)
n4 = max( 1 , 2 )
n5 = ( max(1,2) ) + ( x[1] )
n6 = max(1,2) {a comment}
n7 = max(1,2) "a comment"
n8 = max(1,2) // a comment
n9 = x[1] "[K]"
h = enthalpy(fl$, T=20, P=1E5)
w = max(max(1,2),x[1])
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["n1"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["n2"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["n3"], WithinAbs(16.0, 1e-12));
    CHECK_THAT(run["n4"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["n5"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["n6"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["n7"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["n8"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["n9"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["w"], WithinAbs(5.0, 1e-12));
}

// ============================================================================
// CS-GAP-IF-DIRECTIVE: compile-time $IF / $IFNOT / $ELSE / $ENDIF
// ============================================================================

namespace {

/// Parse only and keep the whole result (errors with their line, resolved directives).
ParseResult parseFull(const std::string& code) {
    EESParser parser;
    return parser.parse(code);
}

/// "line: message" of every parse error.
std::vector<std::string> parseErrorsWithLines(const std::string& code) {
    std::vector<std::string> out;
    for (const auto& e : parseFull(code).errors) out.push_back(std::to_string(e.line) + ": " + e.message);
    return out;
}

bool hasVar(const ModelRun& run, const std::string& name) { return run.vars.count(name) > 0; }

}  // namespace

TEST_CASE("CS-GAP-IF-DIRECTIVE: register reproducer, $ELSE keeps the true branch and the system is square", "[library-gaps][if-directive]") {
    // was: all branches kept, "3 equations and 2 unknowns"
    auto run = runModel("act$ = 'a'\n$if act$='a'\nx = 5\n$else\nx = 7\n$endif\n");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["x"], WithinAbs(5.0, 1e-12));

    // the other value of the string selects the ELSE branch
    auto other = runModel("act$ = 'b'\n$if act$='a'\nx = 5\n$else\nx = 7\n$endif\n");
    REQUIRE(other.solveOk);
    CHECK_THAT(other["x"], WithinAbs(7.0, 1e-12));

    // the string may be set after the branch only if it is set before the directive: here, before
    auto later = runModel("$if act$='a'\nx = 5\n$endif\nact$ = 'a'\n");
    CHECK_FALSE(later.parseOk);
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: register reproducer, without $ELSE a false branch is removed", "[library-gaps][if-directive]") {
    // was: x = 5 solved with SUCCESS although the condition is false
    auto run = runModel("act$ = 'a'\n$if act$='b'\nx = 5\n$endif\ny = 1\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_FALSE(hasVar(run, "x"));
    CHECK_THAT(run["y"], WithinAbs(1.0, 1e-12));

    // the true condition keeps the equation
    auto kept = runModel("act$ = 'a'\n$if act$='a'\nx = 5\n$endif\ny = 1\n");
    REQUIRE(kept.solveOk);
    CHECK_THAT(kept["x"], WithinAbs(5.0, 1e-12));

    // the removed equation is not an equation of the model at all
    auto pr = parseFull("act$ = 'a'\n$if act$='b'\nx = 5\n$endif\ny = 1\n");
    REQUIRE(pr.success);
    CHECK(pr.equationCount == 2);   // act$ = 'a' and y = 1
    CHECK(pr.program.statements.size() == 2);
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: $IFNOT, and the code in the branches of both", "[library-gaps][if-directive]") {
    auto run = runModel(R"(
act$ = 'a'
$ifnot act$='b'
x = 5
$else
x = 7
$endif
$IFNOT act$='a'
y = 5
$ELSE
y = 7
$ENDIF
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(7.0, 1e-12));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: a string variable and the constant are compared without regard to case", "[library-gaps][if-directive]") {
    // the library files rely on it: expandertype_exp$='hermetic' ... $else expanderType_exp$='Hermetic'
    auto run = runModel(R"(
expandertype_exp$ = 'hermetic'
$if expanderType_exp$='Hermetic'
x = 1
$else
x = 2
$endif
FLUID$ = 'ABC'
$IF fluid$ = 'abc'
y = 1
$ELSE
y = 2
$ENDIF
$if Fluid$='aBc'
z = 1
$endif
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["z"], WithinAbs(1.0, 1e-12));

    // accented letters (UTF-8), as in the French models of the ULiege collection
    auto accents = runModel(R"(
activite$ = 'Définie par le métabolisme'
$if activite$='définie par le MÉTABOLISME'
x = 1
$else
x = 2
$endif
energie$ = 'énergie'
$if energie$='ÉNERGIE'
y = 1
$endif
$if energie$='energie'
z = 1
$endif
)");
    REQUIRE(accents.solveOk);
    CHECK_THAT(accents["x"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(accents["y"], WithinAbs(1.0, 1e-12));
    CHECK_FALSE(hasVar(accents, "z"));   // e and é are different letters
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: spellings of the directive that EES files use", "[library-gaps][if-directive]") {
    CHECK(parseErrorsWithLines("a$ = 'x'\n$If a$='x'\nb = 1\n$EndIf\n$ELSE\n") ==
          std::vector<std::string>{"5: $ELSE without a matching $IF or $IFNOT"});

    auto ok = runModel(R"(
a$ = 'x'
$If a$='x'
b1 = 1
$Else
b1 = 2
$EndIf
$iFnOt a$='x'
b2 = 1
$ELSE
b2 = 2
$ENDIF
$IF a$ ='x'
b3 = 1
$endif
$IF a$= 'x'
b4 = 1
$endif
$IF   a$   =   'x'
b5 = 1
$endif
$if a$='x' then
b6 = 1
$else
b6 = 2
$endif
$if a$='y' THEN
b7 = 1
$else a$='x'
b7 = 2
$endif
)");
    REQUIRE(ok.solveOk);
    for (const char* v : {"b1", "b3", "b4", "b5", "b6"}) CHECK_THAT(ok[v], WithinAbs(1.0, 1e-12));
    for (const char* v : {"b2", "b7"}) CHECK_THAT(ok[v], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: comments after the condition, after $ELSE and $ENDIF, directives inside comments", "[library-gaps][if-directive]") {
    auto run = runModel(R"(
a$ = 'x'
$IF a$='x' "the first case"
b1 = 1
$ELSE " a$ = 'y' "
b1 = 2
$ENDIF "end of the case"
$IF a$='x' {a brace comment}
b2 = 1
$ENDIF {another}
$IF a$='x' // a line comment
b3 = 1
$ENDIF // end
$IF a$='y' "it's a comment with an apostrophe"
b4 = 1
$ELSE
b4 = 2
$ENDIF
{ $IF a$='y'  commented out }
b5 = 1
"$ENDIF in a comment"
$IF a$='y'
b6 = 1
{ $ELSE }
$ENDIF
b7 = 1
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["b1"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["b2"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["b3"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["b4"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["b5"], WithinAbs(1.0, 1e-12));
    CHECK_FALSE(hasVar(run, "b6"));
    CHECK_THAT(run["b7"], WithinAbs(1.0, 1e-12));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: the single quotes around the constant are optional, other characters are kept", "[library-gaps][if-directive]") {
    auto run = runModel(R"(
a$ = 'x'
$IF a$=x
b1 = 1
$ENDIF
$IF a$=X
b2 = 1
$ENDIF
$IF a$=y
b3 = 1
$ENDIF
regime$ = 'on\off'
$IF regime$='ON\OFF'
b4 = 1
$ENDIF
$IF regime$='on/off'
b5 = 1
$ENDIF
w$ = 'with spaces and (parentheses)'
$IF w$='With Spaces and (parentheses)'
b6 = 1
$ENDIF
$IF (w$='no') OR (w$='with spaces and (parentheses)')
b7 = 1
$ENDIF
)");
    REQUIRE(run.solveOk);
    for (const char* v : {"b1", "b2", "b4", "b6", "b7"}) CHECK(hasVar(run, v));
    for (const char* v : {"b3", "b5"}) CHECK_FALSE(hasVar(run, v));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: string variables copy constants, a computed string has no compile-time value", "[library-gaps][if-directive]") {
    // B$ = C$ (EES manual: string variables set equal to other string variables)
    auto run = runModel("c$ = 'abc'\nb$ = c$\n$IF b$='ABC'\nx = 1\n$ENDIF\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));

    // set by a CALL output (library model CSL-0166: regime$ comes from a procedure)
    const std::string proc = "PROCEDURE p(a : r$)\n  r$ = 'ON'\nEND\n";
    auto errors = parseErrorsWithLines(proc + "CALL p(1 : regime$)\n$if regime$='ON'\nx = 1\n$endif\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].rfind("5: ", 0) == 0);
    CHECK(errors[0].find("regime$") != std::string::npos);
    CHECK(errors[0].find("CALL") != std::string::npos);
    // set from another computed string
    errors = parseErrorsWithLines(proc + "CALL p(1 : regime$)\nm$ = regime$\n$if m$='ON'\nx = 1\n$endif\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].find("m$") != std::string::npos);
    // never set (EES would take it from the Diagram window)
    errors = parseErrorsWithLines("$if never$='a'\nx = 1\n$endif\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].rfind("1: ", 0) == 0);
    CHECK(errors[0].find("never$") != std::string::npos);
    // set later in the file: EES needs the definition before the directive
    errors = parseErrorsWithLines("$if later$='a'\nx = 1\n$endif\nlater$ = 'a'\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].find("later$") != std::string::npos);
    CHECK(errors[0].find("before") != std::string::npos);
    // a numerical variable is not a string variable
    CHECK_FALSE(parseErrorsWithLines("a = 1\n$if a='1'\nx = 1\n$endif\n").empty());
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: the string set in a removed branch does not count", "[library-gaps][if-directive]") {
    auto run = runModel(R"(
mode$ = 'a'
$if mode$='b'
fluid$ = 'Water'
$else
fluid$ = 'R134a'
$endif
$if fluid$='R134a'
x = 1
$else
x = 2
$endif
$if fluid$='Water'
y = 1
$endif
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));
    CHECK_FALSE(hasVar(run, "y"));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: nested directives", "[library-gaps][if-directive]") {
    const std::string head = "a$ = 'x'\nb$ = 'y'\n";
    // the four combinations
    const std::string nested =
        "$if a$='x'\n"
        "  $if b$='y'\n    r = 11\n  $else\n    r = 10\n  $endif\n"
        "$else\n"
        "  $if b$='y'\n    r = 1\n  $else\n    r = 0\n  $endif\n"
        "$endif\n";
    for (const char* a : {"x", "z"}) {
        for (const char* b : {"y", "z"}) {
            auto run = runModel(std::string("a$ = '") + a + "'\nb$ = '" + b + "'\n" + nested);
            REQUIRE(run.solveOk);
            const double expected = (std::string(a) == "x" ? 10.0 : 0.0) + (std::string(b) == "y" ? 1.0 : 0.0);
            CHECK_THAT(run["r"], WithinAbs(expected, 1e-12));
        }
    }
    // the same directive form for three levels, with equations between the directives
    auto deep = runModel(head + R"(
$if a$='x'
p1 = 1
$if b$='y'
p2 = 2
$ifnot a$='x'
p3 = 3
$else
p4 = 4
$endif
p5 = 5
$endif
p6 = 6
$endif
)");
    REQUIRE(deep.solveOk);
    for (const char* v : {"p1", "p2", "p4", "p5", "p6"}) CHECK(hasVar(deep, v));
    CHECK_FALSE(hasVar(deep, "p3"));

    // a nested directive inside a removed branch is not evaluated: no error for a string with no value there
    auto removed = runModel(head + "$if a$='z'\n$if never$='q'\nx = 1\n$else\nx = 2\n$endif\n$endif\ny = 1\n");
    REQUIRE(removed.parseOk);
    REQUIRE(removed.solveOk);
    CHECK_FALSE(hasVar(removed, "x"));
    // ... nor is its $ELSE branch taken
    CHECK_THAT(removed["y"], WithinAbs(1.0, 1e-12));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: documented flag conditions, in a CoolSolve solve", "[library-gaps][if-directive]") {
    // EES help ($IF and $IFNOT): PARAMETRICTABLE, MIN/MAX (MINMAX), MIN/MAXTABLE (MINMAXTABLE),
    // UNCERTAINTY, UNCERTAINTYTABLE, DIAGRAMWINDOW, MACRO, PROFESSIONAL, TABLERUN#, UNITSYSTEM('XX')
    // A solve is not a table run: only PROFESSIONAL is true (CoolSolve has no feature restriction).
    auto run = runModel(R"(
$IFNOT PARAMETRICTABLE
k1 = 1
$ENDIF
$IF PARAMETRICTABLE
f1 = 1
$ENDIF
$IF ParametricTable='Table 1'
f2 = 1
$ENDIF
$IFNOT parametrictable='Table 1'
k2 = 1
$ENDIF
$IF MINMAX
f3 = 1
$ENDIF
$IF MIN/MAX
f4 = 1
$ENDIF
$IF MINMAXTABLE
f5 = 1
$ENDIF
$IF MIN/MAXTABLE
f6 = 1
$ENDIF
$IF MINMAXTABLE='Table 1'
f7 = 1
$ENDIF
$IF UNCERTAINTY
f8 = 1
$ENDIF
$IF UNCERTAINTYTABLE
f9 = 1
$ENDIF
$IF UNCERTAINTYTABLE='Table 1'
f10 = 1
$ENDIF
$IF DIAGRAMWINDOW
f11 = 1
$ENDIF
$IF MACRO
f12 = 1
$ENDIF
$IFNOT DIAGRAMWINDOW
k3 = 1
$ENDIF
$IF PROFESSIONAL
k4 = 1
$ENDIF
$IF TABLERUN#=1
f13 = 1
$ENDIF
$IF Tablerun#<13
f14 = 1
$ENDIF
$IF TABLERUN# > 0
f15 = 1
$ENDIF
$IFNOT TABLERUN#=1
k5 = 1
$ENDIF
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    for (const char* v : {"k1", "k2", "k3", "k4", "k5"}) CHECK(hasVar(run, v));
    for (int i = 1; i <= 15; ++i) CHECK_FALSE(hasVar(run, "f" + std::to_string(i)));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: UNITSYSTEM and INTEGRALTABLE conditions", "[library-gaps][if-directive]") {
    // CoolSolve: SI, mass basis, degrees, C, Pa, J (as the UNITSYSTEM() function)
    auto run = runModel(R"(
$IF UNITSYSTEM('SI')
u1 = 1
$ENDIF
$IF UNITSYSTEM('C')=1
u2 = 1
$ENDIF
$IF UNITSYSTEM('F')=1
f1 = (100-32)/1.8
$ELSE
f1 = 100
$ENDIF
$IF UNITSYSTEM('K')
f2 = 1
$ENDIF
$IF UNITSYSTEM('kPa')
f3 = 1
$ENDIF
$IF unitsystem('Pa') AND UNITSYSTEM('J') AND UNITSYSTEM('Mass') AND UNITSYSTEM('Deg')
u3 = 1
$ENDIF
$IF UNITSYSTEM('Eng')=0
u4 = 1
$ENDIF
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    for (const char* v : {"u1", "u2", "u3", "u4"}) CHECK(hasVar(run, v));
    CHECK_THAT(run["f1"], WithinAbs(100.0, 1e-12));
    CHECK_FALSE(hasVar(run, "f2"));
    CHECK_FALSE(hasVar(run, "f3"));
    // unknown setting
    auto errors = parseErrorsWithLines("$IF UNITSYSTEM('zz')\nx = 1\n$ENDIF\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].find("zz") != std::string::npos);

    // INTEGRALTABLE: true when a $IntegralTable directive appears in the file, above or below the test
    // (EES compiles in two passes for this condition)
    const std::string test = "$IF INTEGRALTABLE\na = 1\n$ELSE\na = 2\n$ENDIF\n";
    auto without = runModel(test);
    REQUIRE(without.solveOk);
    CHECK_THAT(without["a"], WithinAbs(2.0, 1e-12));
    auto below = parseFull(test + "$IntegralTable t y\n");
    auto above = parseFull("$IntegralTable t y\n" + test);
    for (const auto* pr : {&below, &above}) {
        REQUIRE(pr->success);
        REQUIRE(pr->conditionalDirectives.size() == 1);
        CHECK(pr->conditionalDirectives[0].thenKept);
    }
    // the directive of a removed branch does not make the test true
    auto removed = parseFull("m$ = 'a'\n$IF m$='b'\n$IntegralTable t y\n$ENDIF\n" + test);
    REQUIRE(removed.success);
    REQUIRE(removed.conditionalDirectives.size() == 2);
    CHECK_FALSE(removed.conditionalDirectives[1].thenKept);
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: AND, OR and parentheses in the condition", "[library-gaps][if-directive]") {
    // the example of the EES help
    auto run = runModel(R"(
A$='A'; B$='D'
$if (A$='A') and ((B$='B') or (B$='C'))
z=1
$else
z=2
$endif
$IF ((A$='A') OR (A$='K')) AND (B$='d')
w1 = 1
$ENDIF
$IF A$='A' AND B$='B' OR B$='D'
w2 = 1
$ENDIF
$IF A$='A' OR B$='B' AND B$='C'
w3 = 1
$ENDIF
$IF (A$='X') or (B$='X')
w4 = 1
$ENDIF
$IFNOT (A$='X') or (B$='X')
w5 = 1
$ENDIF
$IF (A$='A') AND PARAMETRICTABLE
w6 = 1
$ENDIF
$IF (A$='A') AND NOT_A_CONDITION
w7 = 1
$ENDIF
)");
    // the last test is not a condition: that is an error, whatever the first term says
    CHECK_FALSE(run.parseOk);
    auto ok = runModel(R"(
A$='A'; B$='D'
$if (A$='A') and ((B$='B') or (B$='C'))
z=1
$else
z=2
$endif
$IF ((A$='A') OR (A$='K')) AND (B$='d')
w1 = 1
$ENDIF
$IF A$='A' AND B$='B' OR B$='D'
w2 = 1
$ENDIF
$IF A$='A' OR B$='B' AND B$='C'
w3 = 1
$ENDIF
$IF (A$='X') or (B$='X')
w4 = 1
$ENDIF
$IFNOT (A$='X') or (B$='X')
w5 = 1
$ENDIF
$IF (A$='A') AND PARAMETRICTABLE
w6 = 1
$ENDIF
)");
    REQUIRE(ok.solveOk);
    CHECK_THAT(ok["z"], WithinAbs(2.0, 1e-12));
    // AND binds tighter than OR; w1: true; w2 = (A and B='B') or B='D': true; w3 = A or (B='B' and B='C'): true
    for (const char* v : {"w1", "w2", "w3", "w5"}) CHECK(hasVar(ok, v));
    for (const char* v : {"w4", "w6"}) CHECK_FALSE(hasVar(ok, v));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: everything in a removed branch is removed, parsed or not", "[library-gaps][if-directive]") {
    // CALL, DUPLICATE, FUNCTION, PROCEDURE, other directives, comments, text CoolSolve cannot parse
    auto run = runModel(R"(
mode$ = 'a'
PROCEDURE p(a : b)
  b = 2*a
END
$if mode$='b'
$Bookmark Never
{ a comment }
"another comment"
PROCEDURE q(a : b)
  b = 3*a
END
FUNCTION g(a)
  g = a
END
CALL p(1 : y1)
CALL q(1 : y2)
DUPLICATE i=1,3
  z[i] = i
END
y3 = g(3)
y4 = (this is not valid EES
MODULE m(a : b)
  b = a
END
$integralautostep 0.1
$else
y5 = 5
$endif
CALL p(4 : v)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y5"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["v"], WithinAbs(8.0, 1e-12));
    for (const char* v : {"y1", "y2", "y3", "y4", "z[1]", "z[2]", "z[3]"}) CHECK_FALSE(hasVar(run, v));
    // no warning (unknown directive, unknown function...) comes from the removed lines either
    for (const auto& w : run.warnings) CHECK(w.find("Never") == std::string::npos);

    // the same lines in a kept branch are parsed and executed
    auto kept = runModel(R"(
mode$ = 'b'
PROCEDURE q(a : b)
  b = 3*a
END
$if mode$='b'
DUPLICATE i=1,3
  z[i] = i
END
CALL q(1 : y2)
$endif
)");
    REQUIRE(kept.solveOk);
    CHECK_THAT(kept["y2"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(kept["z[3]"], WithinAbs(3.0, 1e-12));

    // an error outside a removed branch is still reported, one inside is not
    auto errors = parseErrorsWithLines("mode$ = 'a'\n$if mode$='b'\ny = (1 +\n$endif\nx = (2 +\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].rfind("5: ", 0) == 0);
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: the definition of a function or procedure chosen by a directive", "[library-gaps][if-directive]") {
    auto model = [](const std::string& choice) {
        return "FUN$ = '" + choice + "'\n"
               "$if FUN$='double'\n"
               "FUNCTION f(a)\n  f = 2*a\nEND\n"
               "$else\n"
               "FUNCTION f(a)\n  f = 3*a\nEND\n"
               "$endif\n"
               "y = f(5)\n";
    };
    auto run = runModel(model("double"));
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(10.0, 1e-12));
    auto other = runModel(model("triple"));
    REQUIRE(other.solveOk);
    CHECK_THAT(other["y"], WithinAbs(15.0, 1e-12));
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: directives in FUNCTION and PROCEDURE bodies", "[library-gaps][if-directive]") {
    auto run = runModel(R"(
mode$ = 'b'
FUNCTION f(a)
$if mode$='a'
  f = a
$else
  f = 2*a
$endif
END
PROCEDURE p(a : r)
  r = a
  $ifnot PARAMETRICTABLE
  r = r + 100
  $endif
  $if PARAMETRICTABLE
  r = -1
  $endif
END
y = f(3)
CALL p(1 : z)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["z"], WithinAbs(101.0, 1e-12));

    // a branch that ends in the middle of a body keeps the END of that body: not a block boundary
    auto errors = parseErrorsWithLines("mode$ = 'a'\nFUNCTION f(a)\n$if mode$='b'\n  f = a\nEND\n$else\n  f = 2*a\nEND\n$endif\ny = f(1)\n");
    CHECK(errors.empty());
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: unbalanced directives and unusable conditions are errors with a line", "[library-gaps][if-directive]") {
    using Messages = std::vector<std::string>;
    // $ELSE or $ENDIF without $IF
    CHECK(parseErrorsWithLines("x = 1\n$else\ny = 2\n") == Messages{"2: $ELSE without a matching $IF or $IFNOT"});
    CHECK(parseErrorsWithLines("x = 1\n$endif\n") == Messages{"2: $ENDIF without a matching $IF or $IFNOT"});
    CHECK(parseErrorsWithLines("a$ = 'a'\n$if a$='a'\nx = 1\n$endif\n$endif\n") ==
          Messages{"5: $ENDIF without a matching $IF or $IFNOT"});
    // $IF without $ENDIF: the line of the $IF
    CHECK(parseErrorsWithLines("a$ = 'a'\n$if a$='a'\nx = 1\n") == Messages{"2: $IF without a matching $ENDIF"});
    CHECK(parseErrorsWithLines("a$ = 'a'\n$ifnot a$='b'\nx = 1\n$else\ny = 1\n") == Messages{"2: $IFNOT without a matching $ENDIF"});
    // the outer and the inner one are both open
    auto open = parseErrorsWithLines("a$ = 'a'\n$if a$='a'\n$if a$='a'\nx = 1\n");
    REQUIRE(open.size() == 2);
    CHECK(open[0].rfind("3: ", 0) == 0);
    CHECK(open[1].rfind("2: ", 0) == 0);
    // two $ELSE
    auto twice = parseErrorsWithLines("a$ = 'a'\n$if a$='a'\nx = 1\n$else\ny = 1\n$else\nz = 1\n$endif\n");
    REQUIRE(twice.size() == 1);
    CHECK(twice[0].rfind("6: ", 0) == 0);
    CHECK(twice[0].find("line 2") != std::string::npos);
    // no condition
    auto empty = parseErrorsWithLines("$if\nx = 1\n$endif\n");
    REQUIRE(empty.size() == 1);
    CHECK(empty[0].rfind("1: ", 0) == 0);
    // unknown condition: named
    auto unknown = parseErrorsWithLines("$if SOMETHING\nx = 1\n$endif\n");
    REQUIRE(unknown.size() == 1);
    CHECK(unknown[0].rfind("1: ", 0) == 0);
    CHECK(unknown[0].find("SOMETHING") != std::string::npos);
    CHECK_FALSE(parseErrorsWithLines("x = 1\n$ifnot y > 3\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("a$ = 'a'\n$if a$='a' garbage\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("a$ = 'a'\n$if (a$='a'\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("a$ = 'a'\n$if a$='a')\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("a$ = 'a'\n$if a$='a\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("a$ = 'a'\n$if a$='a' AND\nz = 1\n$endif\n").empty());
    // != and <> are not EES conditions: $IFNOT is
    auto ne = parseErrorsWithLines("a$ = 'a'\n$if a$<>'b'\nz = 1\n$endif\n");
    REQUIRE(ne.size() == 1);
    CHECK(ne[0].find("$IFNOT") != std::string::npos);
    // a table name is a string constant: a string variable is not accepted
    CHECK_FALSE(parseErrorsWithLines("t$ = 'T'\n$if PARAMETRICTABLE=t$\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("$if PROFESSIONAL='x'\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("$if TABLERUN#='x'\nz = 1\n$endif\n").empty());
    CHECK_FALSE(parseErrorsWithLines("$if TABLERUN# = 1.5\nz = 1\n$endif\n").empty());
    // VERSION: CoolSolve has no EES version number
    auto version = parseErrorsWithLines("$if VERSION>10.5\nz = 1\n$endif\n");
    REQUIRE(version.size() == 1);
    CHECK(version[0].find("VERSION") != std::string::npos);
    // $ELSEIF is $ELSE in the EES help; with a condition it is not a documented directive
    CHECK(parseErrorsWithLines("a$ = 'a'\n$if a$='b'\nz = 1\n$elseif\nz = 2\n$endif\n").empty());
    auto elseif = parseErrorsWithLines("a$ = 'a'\n$if a$='b'\nz = 1\n$elseif a$='a'\nz = 2\n$endif\n");
    REQUIRE(elseif.size() == 1);
    CHECK(elseif[0].rfind("4: ", 0) == 0);
    // an error in the condition of the first directive does not stop the report on the next
    CHECK(parseErrorsWithLines("$if A\nx = 1\n$endif\n$if B\ny = 1\n$endif\n").size() == 2);
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: the directives are listed with the branch kept", "[library-gaps][if-directive]") {
    auto pr = parseFull(
        "a$ = 'x'\n"                  // 1
        "$if a$='x'\n"                // 2
        "p = 1\n"                     // 3
        "$else\n"                     // 4
        "p = 2\n"                     // 5
        "$endif\n"                    // 6
        "$ifnot a$='x' \"comment\"\n" // 7
        "q = 1\n"                     // 8
        "$if never$='y'\n"            // 9   (inside a removed branch)
        "q2 = 1\n"
        "$endif\n"
        "$endif\n"                    // 12
        "$IF PARAMETRICTABLE\n"       // 13
        "r = 1\n"
        "$ELSE\n"
        "r = 2\n"
        "$ENDIF\n");                  // 17
    REQUIRE(pr.success);
    REQUIRE(pr.conditionalDirectives.size() == 4);
    const auto& d = pr.conditionalDirectives;
    CHECK(d[0].line == 2);
    CHECK(d[0].text == "$if a$='x'");
    CHECK(d[0].evaluated);
    CHECK(d[0].thenKept);
    CHECK(d[0].elseLine == 4);
    CHECK(d[0].endLine == 6);
    CHECK(d[1].line == 7);
    CHECK(d[1].text == "$ifnot a$='x'");   // the comment is not part of it
    CHECK(d[1].evaluated);
    CHECK_FALSE(d[1].thenKept);              // $IFNOT of a true test
    CHECK(d[1].elseLine == 0);
    CHECK(d[1].endLine == 12);
    CHECK(d[2].line == 9);
    CHECK_FALSE(d[2].evaluated);             // inside the removed branch of the previous one
    CHECK(d[2].endLine == 11);
    CHECK(d[3].line == 13);
    CHECK(d[3].evaluated);
    CHECK_FALSE(d[3].thenKept);
    CHECK(d[3].elseLine == 15);
    CHECK(d[3].endLine == 17);
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: debug folder says which branches were kept", "[library-gaps][if-directive]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_library_gaps_debug";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string code =
        "a$ = 'x'\n"
        "$if a$='x'\n"
        "p = 1\n"
        "$else\n"
        "p = 2\n"
        "$endif\n"
        "$if a$='y'\n"
        "q = 1\n"
        "$endif\n"
        "r = 1\n";
    const fs::path file = dir / "model.eescode";
    { std::ofstream f(file); f << code; }
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    REQUIRE(runner.run(options));
    runner.generateDebugOutput((dir / "debug").string(), code);
    std::ifstream in(dir / "debug" / "report.md");
    REQUIRE(in.good());
    std::string report((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    CHECK(report.find("## Conditional directives") != std::string::npos);
    // one line per directive: its line, its text, the branch kept and the lines removed
    CHECK(report.find("| 2 | `$if a$='x'` | THEN branch kept (line 3), ELSE branch removed (line 5) |") != std::string::npos);
    CHECK(report.find("| 7 | `$if a$='y'` | THEN branch removed (line 8), no $ELSE |") != std::string::npos);
    fs::remove_all(dir);

    // a model without directive has no such section
    const fs::path dir2 = fs::temp_directory_path() / "coolsolve_test_library_gaps_debug2";
    fs::remove_all(dir2);
    fs::create_directories(dir2);
    const fs::path file2 = dir2 / "model.eescode";
    { std::ofstream f(file2); f << "x = 1\n"; }
    CoolSolveRunner runner2(file2.string());
    REQUIRE(runner2.run(options));
    runner2.generateDebugOutput((dir2 / "debug").string(), "x = 1\n");
    std::ifstream in2(dir2 / "debug" / "report.md");
    std::string report2((std::istreambuf_iterator<char>(in2)), std::istreambuf_iterator<char>());
    CHECK(report2.find("Conditional directives") == std::string::npos);
    fs::remove_all(dir2);
}

TEST_CASE("CS-GAP-IF-DIRECTIVE: every entry point goes through the same parser (file, string, INTEGRAL model)", "[library-gaps][if-directive]") {
    const std::string code = "act$ = 'a'\n$if act$='a'\nx = 5\n$else\nx = 7\n$endif\n";
    // EESParser::parse (REST /parse, parametric study) and parseFile (CoolSolveRunner, CLI)
    CHECK(parseFull(code).equationCount == 2);
    const fs::path file = fs::temp_directory_path() / "coolsolve_test_library_gaps_if.eescode";
    { std::ofstream f(file); f << code; }
    EESParser parser;
    auto fromFile = parser.parseFile(file.string());
    fs::remove(file);
    CHECK(fromFile.equationCount == 2);
    REQUIRE(fromFile.conditionalDirectives.size() == 1);

    // an equation-based dynamic model: the branch is chosen before the INTEGRAL extraction
    const std::string integral = R"(
mode$ = 'decay'
$if mode$='decay'
dydt = -y
$else
dydt = y
$endif
y = 1 + integral(dydt, t, 0, 2)
$IntegralTable t:0.5 y
)";
    const fs::path ifile = fs::temp_directory_path() / "coolsolve_test_library_gaps_if_integral.eescode";
    { std::ofstream f(ifile); f << integral; }
    CoolSolveRunner runner(ifile.string());
    SolverOptions options;
    options.tolerance = 1e-9;
    const bool ok = runner.run(options);
    fs::remove(ifile);
    REQUIRE(ok);
    REQUIRE(runner.getIntegralResult().success);
    CHECK_THAT(runner.getIntegralResult().table.interpolate("y", 2.0), WithinAbs(std::exp(-2.0), 1e-3));
}

// ============================================================================
// CS-BUG-STRING-CASE: strings are compared without regard to case
// ============================================================================

TEST_CASE("CS-BUG-STRING-CASE: register reproducer", "[library-gaps][string-case]") {
    // was: z = 2 (EES: 1)
    auto run = runModel("PROCEDURE pc(a$:y)\nIF (a$ = 'ABC') THEN y = 1 ELSE y = 2\nEND\nCALL pc('abc':z)\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["z"], WithinAbs(1.0, 1e-12));
    // a different string is still different
    auto other = runModel("PROCEDURE pc(a$:y)\nIF (a$ = 'ABC') THEN y = 1 ELSE y = 2\nEND\nCALL pc('abd':z)\n");
    REQUIRE(other.solveOk);
    CHECK_THAT(other["z"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-BUG-STRING-CASE: = and <> in FUNCTION and PROCEDURE bodies, every form", "[library-gaps][string-case]") {
    auto run = runModel(R"(
PROCEDURE pc(a$, b$ : y1, y2, y3, y4, y5)
  IF (a$ = 'ABC') THEN y1 = 1 ELSE y1 = 2
  IF (a$ <> 'abc') THEN y2 = 1 ELSE y2 = 2
  IF ('Abc' = a$) THEN y3 = 1 ELSE y3 = 2
  IF (a$ = b$) THEN y4 = 1 ELSE y4 = 2
  IF (a$ <> b$) THEN
    y5 = 1
  ELSE
    y5 = 2
  ENDIF
END
FUNCTION f(a$)
  IF (a$ = 'WATER' OR a$ = 'R134A') THEN
    f = 1
  ELSE
    f = 0
  ENDIF
END
CALL pc('aBc', 'AbC' : r1, r2, r3, r4, r5)
CALL pc('aBc', 'AbD' : s1, s2, s3, s4, s5)
w1 = f('water')
w2 = f('R134a')
w3 = f('Air')
fl$ = 'Water'
w4 = f(fl$)
w5 = f('WaTeR')
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["r1"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["r2"], WithinAbs(2.0, 1e-12));   // <> is false for strings equal but for the case
    CHECK_THAT(run["r3"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["r4"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["r5"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["s4"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["s5"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["w1"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["w2"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["w3"], WithinAbs(0.0, 1e-12));
    CHECK_THAT(run["w4"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["w5"], WithinAbs(1.0, 1e-12));
}

TEST_CASE("CS-BUG-STRING-CASE: accented letters, empty strings, strings that only look alike", "[library-gaps][string-case]") {
    auto run = runModel(R"(
PROCEDURE pc(a$, b$ : y)
  IF (a$ = b$) THEN y = 1 ELSE y = 0
END
CALL pc('Énergie', 'énergie' : r1)
CALL pc('ÉNERGIE', 'énergie' : r2)
CALL pc('Energie', 'énergie' : r3)
CALL pc('abc', 'abcd' : r4)
CALL pc('abc', 'ab' : r5)
CALL pc('', '' : r6)
CALL pc('a b', 'A B' : r7)
CALL pc('a b', 'ab' : r8)
CALL pc('a_b', 'A_B' : r9)
CALL pc('à', 'À' : r10)
)");
    REQUIRE(run.solveOk);
    const std::map<std::string, double> expected = {{"r1", 1}, {"r2", 1}, {"r3", 0}, {"r4", 0}, {"r5", 0},
                                                    {"r6", 1}, {"r7", 1}, {"r8", 0}, {"r9", 1}, {"r10", 1}};
    for (const auto& [name, value] : expected) CHECK_THAT(run[name], WithinAbs(value, 1e-12));
}

TEST_CASE("CS-BUG-STRING-CASE: eesStringEqual, the helper shared by IF conditions and $IF directives", "[library-gaps][string-case]") {
    CHECK(eesStringEqual("abc", "ABC"));
    CHECK(eesStringEqual("AirH2O", "airh2o"));
    CHECK(eesStringEqual("", ""));
    CHECK(eesStringEqual("R134a", "r134A"));
    CHECK(eesStringEqual("Définie", "DÉFINIE"));
    CHECK_FALSE(eesStringEqual("abc", "abd"));
    CHECK_FALSE(eesStringEqual("abc", "abcd"));
    CHECK_FALSE(eesStringEqual("abc", ""));
    CHECK_FALSE(eesStringEqual("é", "e"));
    CHECK_FALSE(eesStringEqual("a b", "ab"));
    // the same equality decides a $IF directive and an IF statement
    auto run = runModel(R"(
PROCEDURE pc(a$:y)
  IF (a$ = 'AIRH2O') THEN y = 1 ELSE y = 2
END
fluid$ = 'airh2o'
CALL pc(fluid$ : z)
$if fluid$ = 'AirH2O'
w = 1
$else
w = 2
$endif
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["z"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["w"], WithinAbs(1.0, 1e-12));
}


// ============================================================================
// Equation-based dynamic (INTEGRAL) models: CS-BUG-INTEGRAL-FACTOR,
// CS-GAP-INTEGRAL-LIMITS, CS-BUG-INTEGRAL-MAXSTEPS, CS-BUG-INTEGRAL-TABLE-SEP,
// CS-BUG-INTEGRAL-TABLE-CASE
// ============================================================================

namespace {

/// Outcome of a run of an INTEGRAL model: result variables (at the upper limit), the trajectory table.
struct IntegralRun {
    bool parseOk = false;
    bool solveOk = false;
    std::string message;                 ///< solver (or parse) error
    std::vector<std::string> warnings;   ///< warning diagnostics of the run
    std::map<std::string, double, CaseInsensitiveLess> vars;
    IntegralTable table;
    int steps = 0;
    int rejected = 0;

    double operator[](const std::string& name) const { return vars.at(name); }
    bool has(const std::string& name) const { return vars.count(name) > 0; }
    /// Tabulated value of column `col` at row `row`.
    double cell(const std::string& col, int row) const { return table.column(col).at(row); }
    bool hasInternalVariable() const {
        for (const auto& [name, value] : vars) {
            (void)value;
            if (name.rfind("__int", 0) == 0) return true;
        }
        return false;
    }
};

/// `tables` maps a lookup-table name to its CSV text (written as `<model>-<name>.csv`).
IntegralRun runIntegral(const std::string& code,
                        const std::function<void(SolverOptions&)>& configure = {},
                        const std::map<std::string, std::string>& tables = {}) {
    static int counter = 0;
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_library_gaps_integral";
    fs::create_directories(dir);
    const std::string stem = "model_" + std::to_string(counter++);
    const fs::path file = dir / (stem + ".eescode");
    { std::ofstream f(file); f << code; }
    std::vector<fs::path> tableFiles;
    for (const auto& [name, csv] : tables) {
        tableFiles.push_back(dir / (stem + "-" + name + ".csv"));
        std::ofstream f(tableFiles.back());
        f << csv;
    }
    IntegralRun run;
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    options.tolerance = 1e-9;
    if (configure) configure(options);
    run.solveOk = runner.run(options);
    run.parseOk = runner.isParseSuccess();
    for (const auto& d : runner.getDiagnostics().items())
        if (d.severity == DiagnosticSeverity::Warning) run.warnings.push_back(d.message);
    if (!run.parseOk && !runner.getParseResult().errors.empty())
        run.message = runner.getParseResult().errors[0].message;
    run.message = run.message.empty() ? runner.getSolveResult().errorMessage : run.message;
    run.vars = runner.getSolveResult().variables;
    const auto& ires = runner.getIntegralResult();
    run.table = ires.table;
    run.steps = ires.totalSteps;
    run.rejected = ires.rejectedSteps;
    fs::remove(file);
    for (const auto& t : tableFiles) fs::remove(t);
    return run;
}

/// The columns and the output interval of the $IntegralTable directive of `code` (parse only).
IntegralTableSpec parseTableSpec(const std::string& body) {
    auto pr = parseFull("$IntegralTable " + body + "\n");
    for (const auto& st : pr.program.statements)
        if (st->is<Directive>() && st->as<Directive>().hasIntegralTableSpec)
            return st->as<Directive>().integralTableSpec;
    return IntegralTableSpec{};
}

using Columns = std::vector<std::string>;

}  // namespace

// ---- CS-BUG-INTEGRAL-FACTOR -------------------------------------------------

TEST_CASE("CS-BUG-INTEGRAL-FACTOR: register reproducer (m_dot*integral) returns m_dot times the integral",
          "[library-gaps][integral-factor]") {
    // was: y = 20 (the factor was dropped); EES: 0.15 * 20 = 3
    auto run = runIntegral("m_dot = 0.15\ny = m_dot*integral(v, P, 0, 10, 1)\nv = 2\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(3.0, 1e-9));
    CHECK_THAT(run["P"], WithinAbs(10.0, 1e-12));   // every variable holds its value at the upper limit
    // the split form (the workaround of the variants) gives the same
    auto split = runIntegral("m_dot = 0.15\nw = 0 + integral(v, P, 0, 10, 1)\ny = m_dot*w\nv = 2\n");
    REQUIRE(split.solveOk);
    CHECK_THAT(split["y"], WithinAbs(3.0, 1e-9));
    // the trajectory: columns of the default table, y(P) = 0.3*P
    CHECK(run.table.columns() == Columns{"P", "y", "v"});
    REQUIRE(run.table.numRows() == 11);
    CHECK_THAT(run.cell("y", 5), WithinAbs(1.5, 1e-9));
    CHECK_THAT(run.cell("P", 5), WithinAbs(5.0, 1e-12));
    // the auxiliary state of the call is internal: not in the results
    CHECK_FALSE(run.hasInternalVariable());
}

TEST_CASE("CS-BUG-INTEGRAL-FACTOR: every form of an INTEGRAL call inside an equation", "[library-gaps][integral-factor]") {
    struct Form { const char* name; const char* code; double expected; };
    const Form forms[] = {
        {"division", "v = 7.2E6\ny = integral(v, t, 0, 10, 1)/3.6E6\n", 20.0},
        {"negative factor", "N = 0.2\nv = 3\ny = -N*integral(v, P, 1, 5, 0.5)\n", -2.4},
        {"unary minus", "v = 3\ny = -integral(v, P, 1, 5, 0.5)\n", -12.0},
        {"a + b*integral", "a = 10\nb = 2\nv = 3\ny = a + b*integral(v, P, 1, 5, 0.5)\n", 34.0},
        {"a - integral", "a = 5\nv = 2\ny = a - integral(v, t, 0, 1, 0.1)\n", 3.0},
        {"integral - a (still additive)", "w = 1\nv = 3\ny = integral(v, t, 0, 2, 0.5) - w\n", 5.0},
        {"two calls in one equation", "f1 = 2\nf2 = t\ny = integral(f1, t, 0, 4, 0.5) - 2*integral(f2, t, 0, 4, 0.5)\n", -8.0},
        {"product of two integrals", "a = 1\nb = 3\ny = integral(a, t, 0, 2, 0.5)*integral(b, t, 0, 2, 0.5)\n", 12.0},
        {"function of the integral", "v = 3\ny = sqrt(integral(v, t, 0, 4, 0.5))\n", std::sqrt(12.0)},
        {"scaled left-hand side", "f1 = 2\n2*y = integral(f1, t, 0, 4, 0.5)\n", 4.0},
        {"call on the left-hand side", "f1 = 2\nintegral(f1, t, 0, 4, 0.5) = y\n", 8.0},
        {"array element on the left", "f1 = 2\ny[1] = 3*integral(f1, t, 0, 1, 0.1)\ny = y[1]\n", 6.0},
        {"integrand that is not a variable", "y = 5 - integral(2, t, 0, 1, 0.1)\n", 3.0},
        {"integrand with an array element", "x[2] = 3\ny = 2*integral(x[2], t, 0, 4, 0.5)\n", 24.0},
        {"integrand expression, factor", "y = 0.5*integral(x^2, x, 0, 3, 0.1)\n", 4.5},
    };
    for (const auto& f : forms) {
        INFO(f.name);
        auto run = runIntegral(f.code);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(f.expected, 1e-6));
        CHECK_FALSE(run.hasInternalVariable());
    }
}

TEST_CASE("CS-BUG-INTEGRAL-FACTOR: the integral feeds back into its own integrand", "[library-gaps][integral-factor]") {
    // y = 2*I, dI/dt = 1 - y  =>  y(t) = 1 - exp(-2 t)
    auto run = runIntegral("y = 2*integral(dydt, t, 0, 1, 0.01)\ndydt = 1 - y\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(1.0 - std::exp(-2.0), 1e-6));
    CHECK_THAT(run.table.interpolate("y", 0.5), WithinAbs(1.0 - std::exp(-1.0), 1e-6));
    // the integrand is a constant expression of the state: dI/dt = -I, I(0) = 0 (an equation of its own)
    auto dec = runIntegral("y = 3 + 2*integral(-z, t, 0, 1, 0.01)\nz = y - 3\n");
    REQUIRE(dec.solveOk);
    CHECK_THAT(dec["y"], WithinAbs(3.0, 1e-9));   // I stays 0
}

TEST_CASE("CS-BUG-INTEGRAL-FACTOR: a base that varies along the trajectory is not an initial value", "[library-gaps][integral-factor]") {
    // EES: y(t) = 2 t + t = 3 at t = 1. The base was evaluated once at t0 (y = 1).
    SECTION("the base uses the integration variable") {
        auto run = runIntegral("y = 2*t + integral(dydt, t, 0, 1, 0.1)\ndydt = 1\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(3.0, 1e-9));
    }
    SECTION("the base depends on t through an algebraic equation") {
        auto run = runIntegral("a = 2*t\ny = a + integral(dydt, t, 0, 1, 0.1)\ndydt = 1\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(3.0, 1e-9));
        CHECK_THAT(run["a"], WithinAbs(2.0, 1e-9));
    }
    SECTION("the base is another state: y2(1) = y1(1) + 1") {
        auto run = runIntegral("y1 = 1 + integral(d1, t, 0, 1, 0.1)\nd1 = 1\ny2 = y1 + integral(d2, t, 0, 1, 0.1)\nd2 = 1\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y1"], WithinAbs(2.0, 1e-9));
        CHECK_THAT(run["y2"], WithinAbs(3.0, 1e-9));
    }
    SECTION("a constant base keeps the plain state form") {
        auto run = runIntegral("y0 = 4\ny = y0 + integral(dydt, t, 0, 1, 0.1)\ndydt = -y\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(4.0 * std::exp(-1.0), 1e-5));
    }
}

TEST_CASE("CS-BUG-INTEGRAL-FACTOR: negative tests", "[library-gaps][integral-factor]") {
    SECTION("wrong arity") {
        auto run = runIntegral("y = 2*integral(dydt, t)\ndydt = 1\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "4 or 5 arguments"));
    }
    SECTION("nested call") {
        auto run = runIntegral("y = 2*integral(integral(a, t, 0, 1), t, 0, 1)\na = 1\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "Nested INTEGRAL"));
    }
    SECTION("a call in a limit") {
        auto run = runIntegral("w = integral(a, t, 0, 1)\ny = integral(a, t, 0, 2*integral(a, t, 0, 1))\na = 1\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "Nested INTEGRAL"));
    }
    SECTION("different integration variables are still rejected (CS-GAP-INTEGRAL-MULTIVAR)") {
        auto run = runIntegral("p = 1\ny = 2*integral(p, v1, 0, 1, 0.1)\nz = 3*integral(p, v2, 0, 1, 0.1)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "same integration variable"));
    }
    SECTION("a negative step") {
        auto run = runIntegral("y = 2*integral(1, t, 0, 1, -0.1)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "step must be positive"));
    }
    SECTION("an equation of a DUPLICATE loop is a clear error, not a wrong value") {
        auto run = runIntegral("DUPLICATE i=1,2\n y[i] = 1 + integral(dydt[i], t, 0, 1, 0.1)\n dydt[i] = -i*y[i]\nEND\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "DUPLICATE"));
    }
}

TEST_CASE("CS-BUG-INTEGRAL-FACTOR: default table of an equation with a factor", "[library-gaps][integral-factor]") {
    auto run = runIntegral("y = 2*integral(dydt, t, 0, 1, 0.25)\ndydt = 4*t\n");
    REQUIRE(run.solveOk);
    CHECK(run.table.columns() == Columns{"t", "y", "dydt"});
    CHECK_THAT(run["y"], WithinAbs(4.0, 1e-9));
    CHECK_THAT(run.table.interpolate("y", 0.5), WithinAbs(1.0, 1e-9));
    CHECK(run.steps == 4);
}

// ---- CS-GAP-INTEGRAL-LIMITS ---------------------------------------------------

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: register reproducer (limits that are variables)", "[library-gaps][integral-limits]") {
    auto run = runIntegral("t_1 = 0\nt_2 = 4\ny = 1 + integral(dydt, t, t_1, t_2)\ndydt = -y\n",
                           [](SolverOptions& o) { o.integralMaxSteps = 400; });
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(std::exp(-4.0), 1e-5));
    CHECK_THAT(run["t"], WithinAbs(4.0, 1e-12));
    CHECK(run.steps == 400);
    CHECK_FALSE(run.hasInternalVariable());
    CHECK(run.table.columns() == Columns{"t", "y", "dydt"});
}

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: constant expressions and the library model CSL-0047", "[library-gaps][integral-limits]") {
    // the native file of CSL-0047: limits are expressions of constants, factor -N in front of the call
    auto run = runIntegral(R"(
N = 0.2
V1 = 2
V2 = V1*2
T_K = 350
R = 8314
v1_bar = V1/N
v2_bar = V2/N
k = 1E4
v_bar*(P + k/v_bar^2) = R*T_K
W_b_EES = -N*integral(P, v_bar, v1_bar, v2_bar, 0.01)
W_b_hand = -N*(R*T_K*ln(v2_bar/v1_bar) + k*(1/v2_bar - 1/v1_bar))
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["W_b_hand"], WithinAbs(-403297.8, 0.1));          // EES value of the README
    CHECK_THAT(run["W_b_EES"], WithinRel(run["W_b_hand"], 1e-8));    // was +2 016 489 J (no -N, limits rejected)
    CHECK_THAT(run["v_bar"], WithinAbs(20.0, 1e-9));                  // at the upper limit
    CHECK_THAT(run["P"], WithinRel(8314.0*350.0/20.0 - 1e4/400.0, 1e-9));
}

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: limits from expressions, a lower limit that is not zero", "[library-gaps][integral-limits]") {
    SECTION("lower limit and upper limit are computed") {
        auto run = runIntegral("t_lo = 1\nt_hi = t_lo + 2\ny = integral(2*t, t, t_lo, t_hi, 0.1)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(8.0, 1e-9));    // t^2 from 1 to 3
        CHECK_THAT(run["t"], WithinAbs(3.0, 1e-12));
        REQUIRE(run.table.numRows() > 1);
        CHECK_THAT(run.cell("t", 0), WithinAbs(1.0, 1e-12));
        CHECK_THAT(run.cell("y", 0), WithinAbs(0.0, 1e-12));
    }
    SECTION("array elements as limits (the Stirling cycle model)") {
        auto run = runIntegral("v[1] = 1\nv[2] = 3*v[1]\ny = integral(1/v12, v12, v[1], v[2], 0.01)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(std::log(3.0), 1e-9));
    }
    SECTION("limits from a user function") {
        auto run = runIntegral("FUNCTION twice(x)\n twice = 2*x\nEND\nhi = twice(2)\ny = integral(t, t, 0, hi, 0.1)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(8.0, 1e-9));
    }
    SECTION("the limit is the solution of a nonlinear equation") {
        auto run = runIntegral("a*exp(a/10) = 3*exp(0.3)\ny = integral(1, x, 0, a, 0.1)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(3.0, 1e-7));
        CHECK_THAT(run["y"], WithinAbs(3.0, 1e-7));
    }
    SECTION("the limit is read from a coupled system, the integration variable has no other equation") {
        auto run = runIntegral("x + y = 7\nx - y = 1\nz = integral(2, s, 0, x, 0.5)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["z"], WithinAbs(8.0, 1e-9));    // x = 4
    }
    SECTION("symbolic limits of a state equation with an additive base and a factor") {
        auto run = runIntegral("tf = 2\nTw = 20 + integral(dTdt, t, 0, tf, 0.01)\ndTdt = -0.5*(Tw - 10)\nQ = 3*integral(dTdt, t, 0, tf, 0.01)\n");
        REQUIRE(run.solveOk);
        const double Tend = 10.0 + 10.0*std::exp(-1.0);
        CHECK_THAT(run["Tw"], WithinAbs(Tend, 1e-6));
        CHECK_THAT(run["Q"], WithinAbs(3.0*(Tend - 20.0), 1e-6));
    }
}

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: steady quadrature, the integration variable is a pressure (CSL-0010)", "[library-gaps][integral-limits]") {
    // isothermal compression of an ideal gas: w = R*T*ln(p_out/p_in), W = m_dot*w
    auto run = runIntegral(R"(
m_dot = 0.5
R = 287
T = 300
p_in = 1E5
p_out = 6.25E5
v_isoT = R*T/P
W_dot_isoT = m_dot*integral(v_isoT, P, p_in, p_out, 1000)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["W_dot_isoT"], WithinRel(0.5*287.0*300.0*std::log(6.25), 1e-8));
    CHECK_THAT(run["P"], WithinAbs(6.25E5, 1e-6));
    CHECK(run.steps == 525);   // (6.25E5 - 1E5)/1000
}

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: the step may be a variable or an expression", "[library-gaps][integral-limits]") {
    SECTION("a variable") {
        auto run = runIntegral("h = 0.25\ny = integral(2*t, t, 0, 1, h)\n");
        REQUIRE(run.solveOk);
        CHECK(run.steps == 4);
        CHECK_THAT(run["y"], WithinAbs(1.0, 1e-12));
    }
    SECTION("an expression") {
        auto run = runIntegral("dt = 0.5\ny = integral(2*t, t, 0, 1, dt/2)\n");
        REQUIRE(run.solveOk);
        CHECK(run.steps == 4);
    }
    SECTION("symbolic limits and a constant step") {
        auto run = runIntegral("hi = 1\ny = integral(2*t, t, 0, hi, 0.1)\n");
        REQUIRE(run.solveOk);
        CHECK(run.steps == 10);
    }
}

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: several calls share the interval, defined by variables", "[library-gaps][integral-limits]") {
    auto run = runIntegral("t_1 = 0\nt_2 = 2\ny = 2*integral(a, t, t_1, t_2, 0.5)\nz = integral(b, t, 0, 2, 0.5)\nw = integral(a, t, 0, t_2, 0.5)\na = 1\nb = 2\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(4.0, 1e-9));
    CHECK_THAT(run["z"], WithinAbs(4.0, 1e-9));
    CHECK_THAT(run["w"], WithinAbs(2.0, 1e-9));
    // different intervals are an error, whether the limits are constants or variables
    auto differ = runIntegral("a = 1\nb = 2\ny = integral(1, t, 0, a, 0.1)\nz = integral(1, t, 0, b, 0.1)\n");
    CHECK_FALSE(differ.solveOk);
    CHECK(anyContains({differ.message}, "same [lo, hi] interval"));
    auto constants = runIntegral("y = integral(1, t, 0, 1, 0.1)\nz = integral(1, t, 0, 2, 0.1)\n");
    CHECK_FALSE(constants.solveOk);
    CHECK(anyContains({constants.message}, "same [lo, hi] interval"));
}

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: a limit that changes along the trajectory is an error", "[library-gaps][integral-limits]") {
    SECTION("the upper limit depends on the integration variable") {
        auto run = runIntegral("tf = 2 + t\ny = 1 + integral(dydt, t, 0, tf)\ndydt = -y\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "upper limit 'tf'"));
        CHECK(anyContains({run.message}, "depends on the integration variable"));
    }
    SECTION("the upper limit depends on a state") {
        auto run = runIntegral("tf = 2 + y\ny = 1 + integral(dydt, t, 0, tf)\ndydt = -y\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "depends on the integration variable"));
    }
    SECTION("the lower limit depends on a state") {
        auto run = runIntegral("t0 = 0.1*y\ny = 1 + integral(dydt, t, t0, 2)\ndydt = -y\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "lower limit"));
    }
    SECTION("the step depends on a state") {
        auto run = runIntegral("h = 0.1*y\ny = 1 + integral(dydt, t, 0, 2, h)\ndydt = -y\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "step"));
    }
    SECTION("through an algebraic loop containing the integration variable") {
        auto run = runIntegral("hi = 1 + 0.5*v\nv = 2 - 0.5*t + 0*hi\ny = integral(1, t, 0, hi, 0.1)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "depends on the integration variable"));
    }
}

TEST_CASE("CS-GAP-INTEGRAL-LIMITS: interval errors", "[library-gaps][integral-limits]") {
    SECTION("upper limit below the lower limit") {
        auto run = runIntegral("a = 2\nb = 1\ny = integral(1, t, a, b, 0.1)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "Invalid integration interval"));
    }
    SECTION("a limit that nothing defines") {
        auto run = runIntegral("y = integral(1, t, 0, hi_undefined, 0.1)\nz = 2\n");
        CHECK_FALSE(run.solveOk);
    }
    SECTION("the equations that define the limit cannot be solved") {
        auto run = runIntegral("a = sqrt(0 - 1)\ny = integral(1, t, 0, a, 0.1)\n");
        CHECK_FALSE(run.solveOk);
    }
}

// ---- CS-BUG-INTEGRAL-MAXSTEPS ----------------------------------------------

TEST_CASE("CS-BUG-INTEGRAL-MAXSTEPS: register reproducer, the explicit step of the model wins", "[library-gaps][integral-maxsteps]") {
    // was: SUCCESS after 4*integralMaxSteps = 4000 steps, trajectory and .sol at t = 4000
    auto run = runIntegral("y = 1 + integral(dydt, t, 0, 5000, 1)\ndydt = -y\n");
    REQUIRE(run.solveOk);
    CHECK(run.steps == 5000);
    CHECK_THAT(run["t"], WithinAbs(5000.0, 1e-9));
    REQUIRE(run.table.numRows() == 5001);
    CHECK_THAT(run.cell("t", 5000), WithinAbs(5000.0, 1e-9));
    CHECK_THAT(run.cell("t", 4000), WithinAbs(4000.0, 1e-9));
}

TEST_CASE("CS-BUG-INTEGRAL-MAXSTEPS: the number of steps follows from (hi - lo)/step", "[library-gaps][integral-maxsteps]") {
    const char* model = "y = integral(2, t, 0, 1, 0.1)\n";
    for (int maxSteps : {3, 10, 1000}) {
        INFO("integralMaxSteps = " << maxSteps);
        auto run = runIntegral(model, [&](SolverOptions& o) { o.integralMaxSteps = maxSteps; });
        REQUIRE(run.solveOk);
        CHECK(run.steps == 10);
        CHECK_THAT(run["y"], WithinAbs(2.0, 1e-12));
    }
    SECTION("a step that does not divide the interval: the last step is shorter") {
        auto run = runIntegral("y = integral(2, t, 0, 1, 0.3)\n");
        REQUIRE(run.solveOk);
        CHECK(run.steps == 4);
        CHECK_THAT(run["t"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["y"], WithinAbs(2.0, 1e-12));
        REQUIRE(run.table.numRows() == 5);
        CHECK_THAT(run.cell("t", 3), WithinAbs(0.9, 1e-12));
    }
    SECTION("a step larger than the interval takes one step") {
        auto run = runIntegral("y = integral(2, t, 0, 1, 5)\n");
        REQUIRE(run.solveOk);
        CHECK(run.steps == 1);
        CHECK_THAT(run["y"], WithinAbs(2.0, 1e-12));
    }
    SECTION("accumulated rounding does not add a spurious step") {
        auto run = runIntegral("y = integral(2, t, 0, 1, 0.01)\n");
        REQUIRE(run.solveOk);
        CHECK(run.steps == 100);
        auto big = runIntegral("y = integral(2, P, 1E5, 6.25E5, 1000)\n");
        REQUIRE(big.solveOk);
        CHECK(big.steps == 525);
    }
    SECTION("the fixed step of the configuration is explicit as well") {
        auto run = runIntegral("y = integral(2, t, 0, 1)\n",
                               [](SolverOptions& o) { o.integralFixedStep = 0.2; o.integralMaxSteps = 3; });
        REQUIRE(run.solveOk);
        CHECK(run.steps == 5);
    }
    SECTION("the step of the model wins over the configuration") {
        auto run = runIntegral("y = integral(2, t, 0, 1, 0.5)\n",
                               [](SolverOptions& o) { o.integralFixedStep = 0.2; });
        REQUIRE(run.solveOk);
        CHECK(run.steps == 2);
    }
}

TEST_CASE("CS-BUG-INTEGRAL-MAXSTEPS: without a step, integralMaxSteps is the number of steps", "[library-gaps][integral-maxsteps]") {
    for (const char* method : {"RK4", "EulerExplicit", "EulerImplicit"}) {
        INFO(method);
        auto run = runIntegral("y = integral(2, t, 0, 1)\n",
                               [&](SolverOptions& o) { o.integralMethod = method; o.integralMaxSteps = 50; });
        REQUIRE(run.solveOk);
        CHECK(run.steps == 50);
        CHECK_THAT(run["t"], WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("CS-BUG-INTEGRAL-MAXSTEPS: an adaptive run that exhausts its budget is an error, not a success", "[library-gaps][integral-maxsteps]") {
    // dy/dt = -50 y needs ~150 RK45 steps over [0, 10]; the budget is 4*integralMaxSteps attempts
    const char* model = "y = 1 + integral(dydt, t, 0, 10)\ndydt = -50*y\n";
    auto starved = runIntegral(model, [](SolverOptions& o) {
        o.integralMethod = "RK45";
        o.integralMaxSteps = 3;
    });
    CHECK_FALSE(starved.solveOk);
    CHECK(anyContains({starved.message}, "integralMaxSteps"));
    CHECK(anyContains({starved.message}, "before reaching the upper limit"));
    // with a sufficient budget the same model succeeds
    auto fed = runIntegral(model, [](SolverOptions& o) {
        o.integralMethod = "RK45";
        o.integralMaxSteps = 1000;
    });
    REQUIRE(fed.solveOk);
    CHECK_THAT(fed["t"], WithinAbs(10.0, 1e-9));
    CHECK(std::abs(fed["y"]) < 1e-6);
}

TEST_CASE("CS-BUG-INTEGRAL-MAXSTEPS: an absurd number of steps is refused up front", "[library-gaps][integral-maxsteps]") {
    auto run = runIntegral("y = integral(2, t, 0, 1, 1E-12)\n");
    CHECK_FALSE(run.solveOk);
    CHECK(anyContains({run.message}, "more than 2e9 steps"));
}

// ---- CS-BUG-INTEGRAL-TABLE-SEP -------------------------------------------------

TEST_CASE("CS-BUG-INTEGRAL-TABLE-SEP: register reproducer, columns separated by ;", "[library-gaps][integral-table]") {
    // was: header "t,y;,dydt" and an empty y column
    auto run = runIntegral("y = 1 + integral(dydt, t, 0, 4)\ndydt = -y\n$IntegralTable t:0.5; y; dydt\n",
                           [](SolverOptions& o) { o.integralMaxSteps = 400; });
    REQUIRE(run.solveOk);
    CHECK(run.table.columns() == Columns{"t", "y", "dydt"});
    CHECK_THAT(run.table.interpolate("y", 0.5), WithinAbs(std::exp(-0.5), 1e-6));
    CHECK_THAT(run.table.interpolate("dydt", 0.5), WithinAbs(-std::exp(-0.5), 1e-6));
    CHECK(run.warnings.empty());
}

TEST_CASE("CS-BUG-INTEGRAL-TABLE-SEP: spaces, commas and semicolons, in any mix", "[library-gaps][integral-table]") {
    const Columns expected{"t", "y", "dydt", "z"};
    for (const char* body : {"t:0.5 y dydt z", "t:0.5; y; dydt; z", "t:0.5, y, dydt, z", "t:0.5,y,dydt,z",
                             "t:0.5;y;dydt;z", "t:0.5, y; dydt z", "t:0.5 ,y ;dydt, z", "  t:0.5\ty,\tdydt;z  ",
                             "t:0.5, y, dydt, z,", "t:0.5; y; dydt; z;", "t:0.5,, y;; dydt, z"}) {
        INFO("$IntegralTable " << body);
        auto spec = parseTableSpec(body);
        CHECK(spec.valid);
        CHECK(spec.columns == expected);
        CHECK_THAT(spec.outputInterval, WithinAbs(0.5, 1e-12));
    }
    // no output interval
    CHECK(parseTableSpec("t, y; dydt").columns == Columns{"t", "y", "dydt"});
    CHECK(parseTableSpec("t, y; dydt").outputInterval == 0.0);
    // array ranges and subscripts: a delimiter inside brackets belongs to the subscript
    CHECK(parseTableSpec("t; X[1..3], z").columns == Columns{"t", "X[1]", "X[2]", "X[3]", "z"});
    CHECK(parseTableSpec("t, A[1,2], B[ 3, 4 ]; C").columns == Columns{"t", "A[1,2]", "B[3,4]", "C"});
    // names of the EES models of the collection
    CHECK(parseTableSpec("tau:1; dU\\dtau; Q_dot_tap; t_tank").columns ==
          Columns{"tau", "dU\\dtau", "Q_dot_tap", "t_tank"});
    CHECK(parseTableSpec("tau:10,  x_ice, u, t_w, DELTAu").columns ==
          Columns{"tau", "x_ice", "u", "t_w", "DELTAu"});
    // a comment after the list is not a column
    CHECK(parseTableSpec("t:0.5; y {the state}").columns == Columns{"t", "y"});
}

TEST_CASE("CS-BUG-INTEGRAL-TABLE-SEP: the same trajectory whatever the separator", "[library-gaps][integral-table]") {
    const std::string model = "y = 1 + integral(dydt, t, 0, 2)\ndydt = -y\nz = 3*y\n";
    auto space = runIntegral(model + "$IntegralTable t:0.5 y dydt z\n", [](SolverOptions& o) { o.integralMaxSteps = 200; });
    auto semi = runIntegral(model + "$IntegralTable t:0.5; y; dydt; z\n", [](SolverOptions& o) { o.integralMaxSteps = 200; });
    auto comma = runIntegral(model + "$IntegralTable t:0.5,y,dydt,z\n", [](SolverOptions& o) { o.integralMaxSteps = 200; });
    REQUIRE(space.solveOk);
    REQUIRE(semi.solveOk);
    REQUIRE(comma.solveOk);
    CHECK(semi.table.toCSV() == space.table.toCSV());
    CHECK(comma.table.toCSV() == space.table.toCSV());
    CHECK(space.table.numRows() == 5);
}

TEST_CASE("CS-BUG-INTEGRAL-TABLE-SEP: an unknown column is reported, not silently empty", "[library-gaps][integral-table]") {
    auto run = runIntegral("y = 1 + integral(dydt, t, 0, 1)\ndydt = -y\n$IntegralTable t:0.5; y; nothere; dydt\n",
                           [](SolverOptions& o) { o.integralMaxSteps = 100; });
    REQUIRE(run.solveOk);
    CHECK(run.warnings.size() == 1);
    CHECK(anyContains(run.warnings, "'nothere' is not a variable of the model"));
    // the known columns are filled, the unknown one stays empty and keeps its place
    CHECK(run.table.columns() == Columns{"t", "y", "nothere", "dydt"});
    CHECK_THAT(run.cell("y", 0), WithinAbs(1.0, 1e-12));
    CHECK(std::isnan(run.cell("nothere", 0)));
    // the separator typo of the register is reported the same way (a name that is not a variable)
    auto typo = runIntegral("y = 1 + integral(dydt, t, 0, 1)\ndydt = -y\n$IntegralTable t:0.5 y: dydt\n",
                            [](SolverOptions& o) { o.integralMaxSteps = 100; });
    REQUIRE(typo.solveOk);
    CHECK(anyContains(typo.warnings, "'y:' is not a variable of the model"));
}

TEST_CASE("CS-BUG-INTEGRAL-TABLE-SEP: a column listed twice, the output interval as a variable", "[library-gaps][integral-table]") {
    SECTION("a column listed twice is kept once") {
        auto run = runIntegral("y = 1 + integral(dydt, t, 0, 1)\ndydt = -y\n$IntegralTable t:0.5 y Y dydt y\n",
                               [](SolverOptions& o) { o.integralMaxSteps = 100; });
        REQUIRE(run.solveOk);
        CHECK(run.table.columns() == Columns{"t", "y", "dydt"});
        CHECK(anyContains(run.warnings, "listed twice"));
        CHECK(run.table.numRows() == 3);
    }
    SECTION("the output interval may be a variable holding a constant (EES)") {
        auto run = runIntegral("dt_out = 0.5\ny = 1 + integral(dydt, t, 0, 2)\ndydt = -y\n$IntegralTable t:dt_out, y\n",
                               [](SolverOptions& o) { o.integralMaxSteps = 200; });
        REQUIRE(run.solveOk);
        CHECK(run.table.numRows() == 5);
        CHECK(run.warnings.empty());
        auto spec = parseTableSpec("tau:DELTAtau_table,tau_summer,time_day");
        CHECK(spec.valid);
        CHECK(spec.columns == Columns{"tau", "tau_summer", "time_day"});
    }
    SECTION("an output interval that is not a variable is reported") {
        auto run = runIntegral("y = 1 + integral(dydt, t, 0, 1)\ndydt = -y\n$IntegralTable t:nothere y\n",
                               [](SolverOptions& o) { o.integralMaxSteps = 10; });
        REQUIRE(run.solveOk);
        CHECK(anyContains(run.warnings, "output interval 'nothere'"));
        CHECK(run.table.numRows() == 11);   // every step
    }
    SECTION("an invalid output interval is a warning of the parser") {
        auto spec = parseTableSpec("t:%$ y");
        CHECK_FALSE(spec.valid);
        CHECK(anyContains({spec.errorMessage}, "invalid output interval"));
    }
}

TEST_CASE("CS-BUG-INTEGRAL-TABLE-SEP: the first name must be the integration variable", "[library-gaps][integral-table]") {
    auto run = runIntegral("y = 1 + integral(dydt, t, 0, 1)\ndydt = -y\n$IntegralTable y:0.5 t\n",
                           [](SolverOptions& o) { o.integralMaxSteps = 100; });
    REQUIRE(run.solveOk);
    CHECK(anyContains(run.warnings, "is not the integration variable"));
}

// ---- CS-BUG-INTEGRAL-TABLE-CASE ---------------------------------------------

TEST_CASE("CS-BUG-INTEGRAL-TABLE-CASE: register reproducer (t_w in the table, T_w in the equations)", "[library-gaps][integral-table]") {
    // was: the t_w column empty
    auto run = runIntegral("u = 0 + INTEGRAL(u_dot, tau, 0, 10)\nu_dot = 1\nq = 5*(5 - T_w)\nt_w = 20\n$IntegralTable tau:10 t_w u\n",
                           [](SolverOptions& o) { o.integralMaxSteps = 10; });
    REQUIRE(run.solveOk);
    REQUIRE(run.table.numRows() == 2);
    CHECK(run.table.columns() == Columns{"tau", "T_w", "u"});   // header: spelling of the model
    CHECK_THAT(run.cell("T_w", 0), WithinAbs(20.0, 1e-12));
    CHECK_THAT(run.cell("T_w", 1), WithinAbs(20.0, 1e-12));
    CHECK_THAT(run.cell("u", 1), WithinAbs(10.0, 1e-9));
    CHECK(run.warnings.empty());
    // the solution file holds the variable once, with the spelling of the model
    CHECK(run.has("T_w"));
    CHECK_THAT(run["t_w"], WithinAbs(20.0, 1e-12));
}

TEST_CASE("CS-BUG-INTEGRAL-TABLE-CASE: state, derivative, integration variable and algebraic columns in any case", "[library-gaps][integral-table]") {
    const std::string model = "Y = 1 + integral(Dydt, TAU, 0, 2)\ndydt = -y\nZ = 3*y\n";
    auto run = runIntegral(model + "$IntegralTable tau:0.5 y DYDT z\n", [](SolverOptions& o) { o.integralMaxSteps = 200; });
    REQUIRE(run.solveOk);
    REQUIRE(run.table.numRows() == 5);
    CHECK(run.table.columns() == Columns{"TAU", "Y", "Dydt", "Z"});
    CHECK_THAT(run.cell("TAU", 2), WithinAbs(1.0, 1e-9));
    CHECK_THAT(run.cell("Y", 4), WithinAbs(std::exp(-2.0), 1e-6));
    CHECK_THAT(run.cell("Dydt", 4), WithinAbs(-std::exp(-2.0), 1e-6));
    CHECK_THAT(run.cell("Z", 4), WithinAbs(3.0*std::exp(-2.0), 1e-6));
    CHECK(run.warnings.empty());
    // the default table (no directive) is unchanged
    auto plain = runIntegral(model, [](SolverOptions& o) { o.integralMaxSteps = 200; });
    REQUIRE(plain.solveOk);
    CHECK(plain.table.columns() == Columns{"TAU", "Y", "Dydt"});
}

TEST_CASE("CS-BUG-INTEGRAL-TABLE-CASE: the integration variable is the same whatever its case in the calls", "[library-gaps][integral-table]") {
    // EES names are case-insensitive: INTEGRAL(..., T, ...) and INTEGRAL(..., t, ...) share one variable
    auto run = runIntegral("y = integral(a, T, 0, 2, 0.5)\nz = integral(b, t, 0, 2, 0.5)\na = 1\nb = 2\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(2.0, 1e-9));
    CHECK_THAT(run["z"], WithinAbs(4.0, 1e-9));
    // genuinely different variables are still rejected
    auto other = runIntegral("y = integral(a, T, 0, 2, 0.5)\nz = integral(b, s, 0, 2, 0.5)\na = 1\nb = 2\n");
    CHECK_FALSE(other.solveOk);
    CHECK(anyContains({other.message}, "same integration variable"));
}

// ============================================================================
// CS-GAP-LOOKUP-PROC: lookup functions inside FUNCTION / PROCEDURE bodies
// ============================================================================

namespace {

/// `Brine2` of the register reproducer (BrineProp mean concentration and temperature of EG).
const std::map<std::string, std::string> kBrine2Table = {{"Brine2", "Column1,Column2\n38.1615,6.3333\n"}};

/// Table `t`: x = 1, 2, 3; y = 10 x; z = 100 x^2 (in rows 1 to 3).
const std::string kXyzCsv = "x,y,z\n1,10,100\n2,20,400\n3,30,900\n";
const std::map<std::string, std::string> kXyzTable = {{"t", kXyzCsv}};

struct ValueAndSlope { double value; double slope; };

/// Evaluate the first equation of `code` (`name = expression` or a CALL whose output is `result`)
/// with the table `t` loaded and the variable `r` independent (value `rv`): the value and the
/// derivative with respect to `r`, straight from the automatic differentiation.
ValueAndSlope evalWithSlope(const std::string& code, const std::string& csv, double rv,
                            const std::string& result = "") {
    EESParser parser;
    auto parsed = parser.parse(code);
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    LookupTableStore store;
    store.add(LookupTable::fromCSV("t", csv));
    ExpressionEvaluator ev(1);
    ev.setLookupTableStore(&store);
    for (const auto& f : ir.getFunctions()) ev.registerFunction(f);
    for (const auto& p : ir.getProcedures()) ev.registerProcedure(p);
    ev.setVariable("r", ADValue::independent(rv, 0, 1));
    const auto& eq = ir.getEquations().at(0);
    ADValue v;
    if (eq.procedureCall) {
        ev.evaluateProcedureCall(*eq.procedureCall);
        v = ev.getVariable(result);
    } else {
        v = ev.evaluate(eq.rhs);
    }
    return {v.value, v.gradient[0]};
}

}  // namespace

TEST_CASE("CS-GAP-LOOKUP-PROC: register reproducer, LOOKUP inside a PROCEDURE", "[library-gaps][lookup-proc]") {
    // was: "lookup(): lookup table 'Brine2' not found (no table store is available in this context)"
    auto run = runModel(R"(
PROCEDURE pk(r:y)
y=lookup('Brine2',r,2)
END
CALL pk(1:yv)
z=yv
w=lookup('Brine2',1,2)
)", kBrine2Table);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["z"], WithinAbs(6.3333, 1e-12));
    CHECK_THAT(run["w"], WithinAbs(6.3333, 1e-12));   // the main program always worked
    CHECK(run.warnings.empty());
}

TEST_CASE("CS-GAP-LOOKUP-PROC: FUNCTION, nested routines, table name given as an argument", "[library-gaps][lookup-proc]") {
    auto run = runModel(R"(
FUNCTION g(r)
  g = lookup('t', r, 2)
END
FUNCTION h(r)
  h = 2*g(r) + 1
END
PROCEDURE inner(r : y)
  y = lookup('t', r, 'z')
END
PROCEDURE outer(r : y1, y2, y3)
  CALL inner(r : y1)
  y2 = g(r)
  y3 = h(r)
END
FUNCTION tv(tab$, r)
  tv = lookup(tab$, r, 2)
END
FUNCTION local_name(r)
  tab$ = 't'
  local_name = lookup(tab$, r, 3)
END
a = g(2.5)
b = h(2.5)
CALL outer(2.5 : c1, c2, c3)
d = tv('t', 2.5)
e = local_name(2)
)", kXyzTable);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinAbs(25.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(51.0, 1e-12));
    CHECK_THAT(run["c1"], WithinAbs(650.0, 1e-9));   // procedure called from a procedure
    CHECK_THAT(run["c2"], WithinAbs(25.0, 1e-12));   // function called from a procedure
    CHECK_THAT(run["c3"], WithinAbs(51.0, 1e-12));   // function called from a function called from a procedure
    CHECK_THAT(run["d"], WithinAbs(25.0, 1e-12));
    CHECK_THAT(run["e"], WithinAbs(400.0, 1e-12));
}

TEST_CASE("CS-GAP-LOOKUP-PROC: REPEAT and DUPLICATE loops and IF conditions of a body (the BrineProp pattern)",
          "[library-gaps][lookup-proc]") {
    // the structure of PROCEDURE BRINEPROP of the ULiege Brineprop.lib: the coefficient rows are read
    // in a REPEAT loop (c[k] = lookup('Brine1', Row[k], Pr)) and the mean values from 'Brine2'
    const std::map<std::string, std::string> tables = {
        {"Brine1", "p1,p2\n1,10\n2,20\n3,30\n4,40\n"},    // 2 fluids x 2 coefficients, 2 properties
        {"Brine2", "xm,ym\n5,6\n7,8\n"}};
    auto run = runModel(R"(
PROCEDURE brinemini(Fl, Pr, Temp : output)
  k = 0
  REPEAT
    k = k + 1
    Row[k] = k + (Fl-1)*2
    c[k] = lookup('Brine1', Row[k], Pr)
  UNTIL (k >= 2)
  xm = lookup('Brine2', Fl, 1)
  ym = lookup('Brine2', Fl, 2)
  y = Temp - ym
  output = c[1] + c[2]*y + xm
END
PROCEDURE sumcol(n : s)
  s = 0
  DUPLICATE i = 1, 4
    v[i] = lookup('Brine1', i, 2)
  END
  s = v[1] + v[2] + v[3] + v[4]
END
PROCEDURE pick(r : y)
  IF (lookup('Brine1', r, 1) > 2.5) THEN
    y = 1
  ELSE
    y = 2
  ENDIF
END
FUNCTION pick1(r)
  pick1 = 2
  IF (lookup('Brine2', r, 1) > 6) THEN pick1 = 1
END
CALL brinemini(2, 2, 10 : out)
CALL sumcol(0 : total)
CALL pick(1 : lo)
CALL pick(4 : hi)
f1 = pick1(1)
f2 = pick1(2)
)", tables);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // Fl = 2, Pr = 2: rows 3 and 4 of column 2 (30, 40), xm = 7, ym = 8, y = 2: 30 + 40*2 + 7
    CHECK_THAT(run["out"], WithinAbs(117.0, 1e-9));
    CHECK_THAT(run["total"], WithinAbs(100.0, 1e-9));
    CHECK_THAT(run["lo"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["hi"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["f1"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["f2"], WithinAbs(1.0, 1e-12));
}

TEST_CASE("CS-GAP-LOOKUP-PROC: every table function works inside a body, as in the main program", "[library-gaps][lookup-proc]") {
    const std::map<std::string, std::string> tables = {
        {"t", kXyzCsv},
        {"g", "x,y,z\n1,1,11\n2,1,12\n1,2,21\n2,2,22\n"}};
    auto run = runModel(R"(
PROCEDURE all(r : a, b, c, d, e, f, g1, h, i, j, k, l)
  a = lookup('t', r, 2)
  b = tablevalue('t', r, 'z')
  c = interpolate('t', 'x', 'z', r)
  d = nlookuprows('t')
  e = nlookupcolumns('t')
  f = lookupcol('t', 'z')
  g1 = sumlookup('t', 'y')
  h = maxlookup('t', 'z')
  i = minlookup('t', 2)
  j = avglookup('t', 'y')
  k = interpolate2('g', 'x', 'y', 'z', 1.5, 1.5)
  l = lookupcellempty('t', 1, 2)
END
CALL all(2.5 : a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1)
a2 = lookup('t', 2.5, 2)
b2 = tablevalue('t', 2.5, 'z')
c2 = interpolate('t', 'x', 'z', 2.5)
d2 = nlookuprows('t')
e2 = nlookupcolumns('t')
f2 = lookupcol('t', 'z')
g2 = sumlookup('t', 'y')
h2 = maxlookup('t', 'z')
i2 = minlookup('t', 2)
j2 = avglookup('t', 'y')
k2 = interpolate2('g', 'x', 'y', 'z', 1.5, 1.5)
l2 = lookupcellempty('t', 1, 2)
)", tables);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    const double expected[] = {25.0, 650.0, 650.0, 3.0, 3.0, 3.0, 60.0, 900.0, 10.0, 20.0, 16.5, 0.0};
    const char* suffix[] = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"};
    for (int n = 0; n < 12; ++n) {
        INFO("function #" << n + 1);
        CHECK_THAT(run[std::string(suffix[n]) + "1"], WithinAbs(expected[n], 1e-9));
        CHECK_THAT(run[std::string(suffix[n]) + "2"], WithinAbs(expected[n], 1e-9));
    }
}

TEST_CASE("CS-GAP-LOOKUP-PROC: the derivative through a routine is analytic", "[library-gaps][lookup-proc]") {
    // LOOKUP at a fractional row interpolates: the slope of column y is 10, of column z 500 in [2, 3]
    const std::string proc = "PROCEDURE pk(r : y)\n  y = lookup('t', r, 'z')\nEND\nCALL pk(r : yv)\n";
    auto viaProc = evalWithSlope(proc, kXyzCsv, 2.5, "yv");
    CHECK_THAT(viaProc.value, WithinAbs(650.0, 1e-9));
    CHECK_THAT(viaProc.slope, WithinAbs(500.0, 1e-9));
    const std::string func = "FUNCTION g(r)\n  g = interpolate('t', 'x', 'y', 2*r)\nEND\nz = g(r)\n";
    auto viaFunc = evalWithSlope(func, kXyzCsv, 1.25, "");
    CHECK_THAT(viaFunc.value, WithinAbs(25.0, 1e-9));
    CHECK_THAT(viaFunc.slope, WithinAbs(20.0, 1e-9));
    // nested: the chain rule crosses a procedure and a function
    const std::string nested =
        "FUNCTION g(r)\n  g = lookup('t', r, 'y')\nEND\n"
        "PROCEDURE pk(r : y)\n  y = 3*g(r)\nEND\nCALL pk(r : yv)\n";
    auto viaBoth = evalWithSlope(nested, kXyzCsv, 1.5, "yv");
    CHECK_THAT(viaBoth.value, WithinAbs(45.0, 1e-9));
    CHECK_THAT(viaBoth.slope, WithinAbs(30.0, 1e-9));
}

TEST_CASE("CS-GAP-LOOKUP-PROC: Newton converges when the routine output is an unknown of a block", "[library-gaps][lookup-proc]") {
    // r is found from the output of a routine that interpolates a table
    auto viaProc = runModel(R"(
PROCEDURE pk(r : y)
  y = lookup('t', r, 'y')
END
CALL pk(r : yv)
yv = 25
)", kXyzTable);
    REQUIRE(viaProc.solveOk);
    CHECK(viaProc.verified);
    CHECK_THAT(viaProc["r"], WithinAbs(2.5, 1e-9));

    auto viaFunc = runModel(R"(
FUNCTION g(x)
  g = interpolate('t', 'x', 'z', x)
END
v = g(x)
v = 650
)", kXyzTable);
    REQUIRE(viaFunc.solveOk);
    CHECK(viaFunc.verified);
    CHECK_THAT(viaFunc["x"], WithinAbs(2.5, 1e-9));

    // a coupled block: the table value and a linear equation, two unknowns solved together
    auto coupled = runModel(R"(
PROCEDURE pk(r : y)
  y = lookup('t', r, 'y')
END
CALL pk(r : yv)
yv + r = 27
)", kXyzTable);
    REQUIRE(coupled.solveOk);
    CHECK(coupled.verified);
    CHECK_THAT(coupled["r"], WithinAbs(27.0 / 11.0, 1e-9));   // 10 r + r = 27
    CHECK_THAT(coupled["yv"], WithinAbs(270.0 / 11.0, 1e-9));
}

TEST_CASE("CS-GAP-LOOKUP-PROC: the table store given to a Solver reaches the routines", "[library-gaps][lookup-proc]") {
    // what the parametric study does: its own Solver with a store built from the session tables
    EESParser parser;
    auto parsed = parser.parse("PROCEDURE pk(r : y)\n  y = lookup('t', r, 'y')\nEND\nCALL pk(2.5 : yv)\nw = 2*yv\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);
    LookupTableStore store;
    store.add(LookupTable::fromCSV("t", kXyzCsv));
    SolverOptions options;
    Solver solver(ir, analysis, options.coolpropConfig);
    solver.setLookupTableStore(&store);
    auto result = solver.solve(options, false);
    REQUIRE(result.success);
    CHECK_THAT(result.variables.at("w"), WithinAbs(50.0, 1e-9));
}

TEST_CASE("CS-GAP-LOOKUP-PROC: negative tests, errors name the table and the routine", "[library-gaps][lookup-proc]") {
    SECTION("unknown table in a PROCEDURE") {
        auto run = runModel("PROCEDURE pk(r : y)\n  y = lookup('nope', r, 2)\nEND\nCALL pk(1 : yv)\n", kBrine2Table);
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "lookup table 'nope' not found"));
        CHECK(anyContains({run.message}, "PROCEDURE 'pk'"));
    }
    SECTION("unknown table in a FUNCTION") {
        auto run = runModel("FUNCTION gf(r)\n  gf = interpolate('nope', 'x', 'y', r)\nEND\nz = gf(1)\n", kXyzTable);
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "lookup table 'nope' not found"));
        CHECK(anyContains({run.message}, "FUNCTION 'gf'"));
    }
    SECTION("unknown table in a nested routine: the innermost one is named") {
        auto run = runModel(
            "FUNCTION inner_f(r)\n  inner_f = lookup('nope', r, 1)\nEND\n"
            "PROCEDURE outer_p(r : y)\n  y = inner_f(r)\nEND\nCALL outer_p(1 : yv)\n", kXyzTable);
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "lookup table 'nope' not found"));
        CHECK(anyContains({run.message}, "FUNCTION 'inner_f'"));
    }
    SECTION("the file name to provide is the one the loader expects") {
        auto run = runModel("PROCEDURE pk(r : y)\n  y = lookup('nope', r, 2)\nEND\nCALL pk(1 : yv)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "-nope.csv"));
    }
    SECTION("unknown column in a body") {
        auto run = runModel("PROCEDURE pk(r : y)\n  y = lookup('t', r, 'nocol')\nEND\nCALL pk(1 : yv)\n", kXyzTable);
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "column 'nocol' not found in table 't'"));
        CHECK(anyContains({run.message}, "PROCEDURE 'pk'"));
    }
    SECTION("wrong number of arguments, in a body and in the main program") {
        auto inBody = runModel("PROCEDURE pk(r : y)\n  y = lookup('t', r)\nEND\nCALL pk(1 : yv)\n", kXyzTable);
        CHECK_FALSE(inBody.solveOk);
        CHECK(anyContains({inBody.message}, "expects 3 arguments"));
        CHECK(anyContains({inBody.message}, "PROCEDURE 'pk'"));
        auto inMain = runModel("z = lookup('t', 1)\n", kXyzTable);   // was: "Unknown expression type"
        CHECK_FALSE(inMain.solveOk);
        CHECK(anyContains({inMain.message}, "expects 3 arguments"));
        auto interp = runModel("FUNCTION gf(r)\n  gf = interpolate('t', 'x', r)\nEND\nz = gf(1)\n", kXyzTable);
        CHECK_FALSE(interp.solveOk);
        CHECK(anyContains({interp.message}, "expects 4 arguments"));
        CHECK(anyContains({interp.message}, "FUNCTION 'gf'"));
    }
    SECTION("an evaluator without table store says so and names the routine") {
        EESParser parser;
        auto parsed = parser.parse("PROCEDURE pk(r : y)\n  y = lookup('t', r, 2)\nEND\nCALL pk(1 : yv)\n");
        REQUIRE(parsed.success);
        IR ir = IR::fromAST(parsed.program);
        ExpressionEvaluator ev(0);
        for (const auto& p : ir.getProcedures()) ev.registerProcedure(p);
        try {
            ev.evaluateProcedureCall(*ir.getEquations().at(0).procedureCall);
            FAIL("expected an error");
        } catch (const std::runtime_error& e) {
            const std::string msg = e.what();
            CHECK(msg.find("no table store is available") != std::string::npos);
            CHECK(msg.find("PROCEDURE 'pk'") != std::string::npos);
        }
    }
}

TEST_CASE("CS-GAP-LOOKUP-PROC: debug folder evaluates the lookups of the residual report", "[library-gaps][lookup-proc]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_library_gaps_lookup_debug";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string code =
        "PROCEDURE pk(r : y)\n  y = lookup('t', r, 'y')\nEND\n"
        "CALL pk(2.5 : yv)\n"
        "a = lookup('t', 2.5, 'z')\n";
    const fs::path file = dir / "model.eescode";
    { std::ofstream f(file); f << code; }
    { std::ofstream f(dir / "model-t.csv"); f << kXyzCsv; }
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    REQUIRE(runner.run(options));
    runner.generateDebugOutput((dir / "debug").string(), code);
    std::ifstream in(dir / "debug" / "coolsolve_residuals.md");
    REQUIRE(in.good());
    std::string line, aRow, callRow;
    while (std::getline(in, line)) {
        if (line.find("a = lookup") != std::string::npos) aRow = line;
        if (line.find("CALL pk") != std::string::npos) callRow = line;
    }
    // was: "| - | - | -" on both rows (the evaluator of the report had no table store)
    REQUIRE_FALSE(aRow.empty());
    CHECK(aRow.find("| - | - | -") == std::string::npos);
    CHECK(aRow.find("6.5000e+02") != std::string::npos);
    REQUIRE_FALSE(callRow.empty());
    CHECK(callRow.find("| - | - | -") == std::string::npos);
    fs::remove_all(dir);
}

// ============================================================================
// CS-BUG-LOOKUP-COL-ARG: a column named by a string variable (found when the routines could read tables)
// ============================================================================

TEST_CASE("CS-BUG-LOOKUP-COL-ARG: register reproducer, the column name comes through an argument", "[library-gaps][lookup-col-arg]") {
    // was: x1 = 38.16 (the first column, silently) instead of 6.33
    const std::map<std::string, std::string> tab = {{"tab", "W,M\n38.16,6.33\n"}};
    auto run = runModel(R"(
FUNCTION g1(tab$, col$, r)
  g1 = lookup(tab$, r, col$)
END
x1 = g1('tab', 'M', 1)
x2 = lookup('tab', 1, 'M')
)", tab);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["x1"], WithinAbs(6.33, 1e-12));
    CHECK_THAT(run["x2"], WithinAbs(6.33, 1e-12));
}

TEST_CASE("CS-BUG-LOOKUP-COL-ARG: string variables as column names, in the main program and in bodies", "[library-gaps][lookup-col-arg]") {
    auto run = runModel(R"(
PROCEDURE rd(pipe$, r : v)
  v = lookup('t', r, pipe$)
END
FUNCTION local_col(r)
  c$ = 'z'
  local_col = lookup('t', r, c$)
END
FUNCTION interp(xc$, yc$, x)
  interp = interpolate('t', xc$, yc$, x)
END
FUNCTION agg(c$)
  agg = sumlookup('t', c$)
END
col$ = 'Y'
a = lookup('t', 2, col$)
CALL rd('Z', 3 : b)
c = local_col(2)
d = interp('x', 'z', 2.5)
e = agg('y')
f = lookup('t', 2, 'y')
g = lookup('t', 2, 2)
)", kXyzTable);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinAbs(20.0, 1e-12));    // main program, string variable (was the first column: 2)
    CHECK_THAT(run["b"], WithinAbs(900.0, 1e-12));   // table literal, column in a string argument (the Pipedata pattern)
    CHECK_THAT(run["c"], WithinAbs(400.0, 1e-12));   // string variable local to the body
    CHECK_THAT(run["d"], WithinAbs(650.0, 1e-12));   // INTERPOLATE with both columns by name
    CHECK_THAT(run["e"], WithinAbs(60.0, 1e-12));    // aggregate
    CHECK_THAT(run["f"], WithinAbs(20.0, 1e-12));    // literal name, any case
    CHECK_THAT(run["g"], WithinAbs(20.0, 1e-12));    // numeric column unchanged
    // part (b) of CS-BUG-LOOKUP-COLNAME: the Pipedata layout, accented column names, name in a string variable
    auto pipe = runModel(R"(
pipe$ = 'ø88.9 Steel Seamless'
x = lookup('pipedata', 1, pipe$)
y = lookup('pipedata', 1, 'ø88.9 Steel Seamless')
)", {{"pipedata", "ø101.6 Steel Seamless,ø88.9 Steel Seamless\n76.1,88.9\n3.2,3.2\n"}});
    REQUIRE(pipe.solveOk);
    CHECK_THAT(pipe["x"], WithinAbs(88.9, 1e-12));   // was the first column, 76.1
    CHECK_THAT(pipe["y"], WithinAbs(88.9, 1e-12));
    // a name that is not a column is an error naming it, not a silent first column
    auto bad = runModel("FUNCTION g(c$)\n  g = lookup('t', 1, c$)\nEND\nx = g('nocol')\n", kXyzTable);
    CHECK_FALSE(bad.solveOk);
    CHECK(anyContains({bad.message}, "column 'nocol' not found in table 't'"));
    CHECK(anyContains({bad.message}, "FUNCTION 'g'"));
}

// ============================================================================
// CS-GAP-INTEGRAL-LOOKUP: lookup functions inside an INTEGRAL model
// ============================================================================

TEST_CASE("CS-GAP-INTEGRAL-LOOKUP: INTERPOLATE as the forcing of a state", "[library-gaps][integral-lookup]") {
    // was: "Initial algebraic solve failed: ... INTERPOLATE(): lookup table 'f' not found (no table
    // store is available in this context)". dy/dt = v(t) = 2 t, y(0) = 0: y(10) = 100 (RK4 integrates it exactly).
    const std::map<std::string, std::string> tables = {{"f", "t,v\n0,0\n5,10\n10,20\n"}};
    auto run = runIntegral(
        "y = 0 + integral(dydt, t, 0, 10, 0.5)\n"
        "dydt = interpolate('f', 't', 'v', t)\n"
        "$IntegralTable t:2.5 y dydt\n",
        [](SolverOptions& o) { o.integralMethod = "RK4"; }, tables);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(100.0, 1e-9));
    CHECK_THAT(run["t"], WithinAbs(10.0, 1e-12));
    CHECK(run.warnings.empty());
    // the trajectory carries the table-driven derivative: y(5) = 25 and dy/dt(5) = 10
    REQUIRE(run.table.numRows() == 5);
    CHECK_THAT(run.cell("t", 2), WithinAbs(5.0, 1e-9));
    CHECK_THAT(run.cell("y", 2), WithinAbs(25.0, 1e-9));
    CHECK_THAT(run.cell("dydt", 2), WithinAbs(10.0, 1e-9));
}

TEST_CASE("CS-GAP-INTEGRAL-LOOKUP: a weather table drives a thermal state (the pattern of building_rc_network)",
          "[library-gaps][integral-lookup]") {
    // dT/dt = (T_out(t) - T)/10 with T_out(t) = 30 for the whole table: T(t) = 30 - 10 exp(-t/10)
    const std::map<std::string, std::string> tables = {{"week", "tau,T_out,f\n0,30,0\n20,30,1\n40,30,0\n"}};
    auto run = runIntegral(
        "T = 20 + integral(dTdt, tau, 0, 40, 0.1)\n"
        "dTdt = (T_out - T)/10\n"
        "T_out = INTERPOLATE('week', 'tau', 'T_out', tau)\n"
        "f_occ = INTERPOLATE('week', 'tau', 'f', tau)\n"
        "$IntegralTable tau:20 T T_out f_occ\n",
        [](SolverOptions& o) { o.integralMethod = "RK4"; }, tables);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["T"], WithinAbs(30.0 - 10.0 * std::exp(-4.0), 1e-5));
    CHECK_THAT(run["T_out"], WithinAbs(30.0, 1e-12));
    REQUIRE(run.table.numRows() == 3);
    CHECK(run.table.columns() == Columns{"tau", "T", "T_out", "f_occ"});
    CHECK_THAT(run.cell("T", 1), WithinAbs(30.0 - 10.0 * std::exp(-2.0), 1e-5));
    CHECK_THAT(run.cell("f_occ", 1), WithinAbs(1.0, 1e-12));
    CHECK_THAT(run.cell("f_occ", 2), WithinAbs(0.0, 1e-12));
}

TEST_CASE("CS-GAP-INTEGRAL-LOOKUP: lookups in routines of an INTEGRAL model, and in its limits", "[library-gaps][integral-lookup]") {
    const std::map<std::string, std::string> tables = {
        {"f", "t,v\n0,0\n5,10\n10,20\n"},
        {"cfg", "tf,step\n8,0.5\n"}};
    // the forcing is read inside a FUNCTION and inside a PROCEDURE (the two kinds of body)
    auto viaFunc = runIntegral(
        "FUNCTION force(x)\n  force = interpolate('f', 't', 'v', x)\nEND\n"
        "y = integral(dydt, t, 0, 10, 0.5)\n"
        "dydt = force(t)\n", {}, tables);
    REQUIRE(viaFunc.solveOk);
    CHECK_THAT(viaFunc["y"], WithinAbs(100.0, 1e-6));
    auto viaProc = runIntegral(
        "PROCEDURE force(x : v)\n  v = lookup('f', 1 + x/5, 'v')\nEND\n"
        "y = integral(dydt, t, 0, 10, 0.5)\n"
        "CALL force(t : dydt)\n", {}, tables);
    REQUIRE(viaProc.solveOk);
    CHECK_THAT(viaProc["y"], WithinAbs(100.0, 1e-6));
    // the limit and the step are read from a table: they are solved on a second algebraic system
    auto limits = runIntegral(
        "tf = lookup('cfg', 1, 'tf')\n"
        "dt = lookup('cfg', 1, 'step')\n"
        "y = integral(v, t, 0, tf, dt)\n"
        "v = 2\n", {}, tables);
    REQUIRE(limits.solveOk);
    CHECK_THAT(limits["y"], WithinAbs(16.0, 1e-9));
    CHECK_THAT(limits["t"], WithinAbs(8.0, 1e-12));
    CHECK(limits.steps == 16);
}

TEST_CASE("CS-GAP-INTEGRAL-LOOKUP: negative tests", "[library-gaps][integral-lookup]") {
    // an unknown table is reported with its name (not silently ignored), in the model and in a routine
    auto direct = runIntegral("y = integral(dydt, t, 0, 1, 0.5)\ndydt = interpolate('nope', 't', 'v', t)\n");
    CHECK_FALSE(direct.solveOk);
    CHECK(anyContains({direct.message}, "lookup table 'nope' not found"));
    auto inFunc = runIntegral(
        "FUNCTION force(x)\n  force = lookup('nope', x, 1)\nEND\n"
        "y = integral(dydt, t, 0, 1, 0.5)\ndydt = force(t)\n");
    CHECK_FALSE(inFunc.solveOk);
    CHECK(anyContains({inFunc.message}, "lookup table 'nope' not found"));
    CHECK(anyContains({inFunc.message}, "FUNCTION 'force'"));
    // a column that does not exist
    const std::map<std::string, std::string> tables = {{"f", "t,v\n0,0\n10,20\n"}};
    auto column = runIntegral("y = integral(dydt, t, 0, 1, 0.5)\ndydt = interpolate('f', 't', 'nocol', t)\n", {}, tables);
    CHECK_FALSE(column.solveOk);
    CHECK(anyContains({column.message}, "column 'nocol' not found in table 'f'"));
    // a lookup written inline in the initial value of a state is refused (it is not evaluated by the
    // algebraic solver); the value written as a variable works
    const std::map<std::string, std::string> init = {{"c", "a,b\n5,7\n"}};
    auto inlineBase = runIntegral("y = lookup('c', 1, 1) + integral(dydt, t, 0, 2, 0.5)\ndydt = 3\n", {}, init);
    CHECK_FALSE(inlineBase.solveOk);
    CHECK(anyContains({inlineBase.message}, "Initial value of a state"));
    auto viaVariable = runIntegral("y0 = lookup('c', 1, 1)\ny = y0 + integral(dydt, t, 0, 2, 0.5)\ndydt = 3\n", {}, init);
    REQUIRE(viaVariable.solveOk);
    CHECK_THAT(viaVariable["y"], WithinAbs(11.0, 1e-9));
    // the table of an INTEGRAL model without any lookup is untouched
    auto plain = runIntegral("y = integral(dydt, t, 0, 1, 0.5)\ndydt = 2\n");
    REQUIRE(plain.solveOk);
    CHECK_THAT(plain["y"], WithinAbs(2.0, 1e-9));
}

// ============================================================================
// CS-GAP-CALL-EXPR-OUT: an expression as an output of a CALL (EES: the output is equated to it)
// ============================================================================

namespace {

/// A procedure with one output, 2*a, shared by the cases below.
const std::string kDoubler = "PROCEDURE dbl(a : y)\n  y = 2*a\nEND\n";

/// Names of the variables of a run that start with the prefix of the auxiliary outputs.
std::vector<std::string> auxiliaryNames(const ModelRun& run) {
    std::vector<std::string> names;
    for (const auto& [name, value] : run.vars)
        if (name.rfind("__", 0) == 0) names.push_back(name);
    return names;
}

/// Residual vector and Jacobian of the block that holds `blockVar`, evaluated at `point`.
struct BlockSnapshot {
    std::vector<std::string> vars;                 ///< columns
    std::vector<std::string> equations;            ///< original text of the rows
    std::vector<std::vector<double>> jacobian;     ///< [row][column]

    int column(const std::string& prefix) const {
        for (size_t j = 0; j < vars.size(); ++j)
            if (vars[j].rfind(prefix, 0) == 0) return static_cast<int>(j);
        return -1;
    }
    int row(const std::string& text) const {
        for (size_t i = 0; i < equations.size(); ++i)
            if (equations[i].find(text) != std::string::npos) return static_cast<int>(i);
        return -1;
    }
};

BlockSnapshot snapshotOfBlock(const std::string& code, const std::string& blockVar,
                              const std::map<std::string, double>& point) {
    EESParser parser;
    auto parsed = parser.parse(code);
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);
    SystemEvaluator system(ir, analysis);
    system.initializeFromGuesses();
    for (const auto& [name, value] : point) system.setVariableValue(name, value);
    BlockSnapshot snap;
    for (size_t b = 0; b < system.getNumBlocks(); ++b) {
        const auto& vars = system.getBlock(b).getVariables();
        if (std::find(vars.begin(), vars.end(), blockVar) == vars.end()) continue;
        snap.vars = vars;
        for (int id : system.getBlock(b).getEquationIds())
            snap.equations.push_back(ir.getEquations()[id].originalText);
        snap.jacobian = system.evaluateBlock(b).jacobian;
        return snap;
    }
    FAIL("no block holds " << blockVar);
    return snap;
}

}  // namespace

TEST_CASE("CS-GAP-CALL-EXPR-OUT: register reproducer, the output is equated to the expression",
          "[library-gaps][call-expr-out]") {
    // valid EES: z/1000 = 2 * 21, so z = 42000 (was: parse error 'Output 1 of 'pg' must be a variable')
    auto run = runModel(R"(
PROCEDURE pg(a:y)
y=a*2
END
CALL pg(21:z/1000)
w=z
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["z"], WithinAbs(42000.0, 1e-6));
    CHECK_THAT(run["w"], WithinAbs(42000.0, 1e-6));

    // the two directions of a unit conversion
    auto scaled = runModel(kDoubler + "CALL dbl(4 : mu*1000)\nCALL dbl(4 : c/1000)\n");
    REQUIRE(scaled.solveOk);
    CHECK_THAT(scaled["mu"], WithinAbs(0.008, 1e-12));
    CHECK_THAT(scaled["c"], WithinAbs(8000.0, 1e-9));
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: forms of expression", "[library-gaps][call-expr-out]") {
    SECTION("sign, offset, array element, power, function, parentheses") {
        auto run = runModel(kDoubler + R"(
CALL dbl(5 : -x)
CALL dbl(300 : T[2] + 273.15)
CALL dbl(8 : q^3)
CALL dbl(0.5 : exp(v))
CALL dbl(3 : (u + 1)/2)
CALL dbl(2 : 4*s - 1)
)");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinAbs(-10.0, 1e-9));
        CHECK_THAT(run["T[2]"], WithinAbs(600.0 - 273.15, 1e-9));
        CHECK_THAT(run["q"], WithinAbs(std::cbrt(16.0), 1e-9));
        CHECK_THAT(run["v"], WithinAbs(0.0, 1e-9));         // exp(v) = 1
        CHECK_THAT(run["u"], WithinAbs(11.0, 1e-9));        // (u+1)/2 = 6
        CHECK_THAT(run["s"], WithinAbs(5.0 / 4.0, 1e-9));   // 4 s - 1 = 4
    }

    SECTION("a plain variable and an array element are unchanged") {
        auto run = runModel(kDoubler + "CALL dbl(2 : r)\nCALL dbl(3 : y[4])\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["r"], WithinAbs(4.0, 1e-12));
        CHECK_THAT(run["y[4]"], WithinAbs(6.0, 1e-12));
        CHECK(auxiliaryNames(run).empty());
        CHECK(run.vars.size() == 2);
    }

    SECTION("several outputs, expressions and plain variables mixed, in order") {
        auto run = runModel(R"(
PROCEDURE three(a : p, q, r)
  p = 10*a
  q = 20*a
  r = 30*a
END
CALL three(1 : x/2, y, z*3)
CALL three(2 : u, v + 1, w)
)");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinAbs(20.0, 1e-9));
        CHECK_THAT(run["y"], WithinAbs(20.0, 1e-9));
        CHECK_THAT(run["z"], WithinAbs(10.0, 1e-9));
        CHECK_THAT(run["u"], WithinAbs(20.0, 1e-9));
        CHECK_THAT(run["v"], WithinAbs(39.0, 1e-9));
        CHECK_THAT(run["w"], WithinAbs(60.0, 1e-9));
    }

    SECTION("the library pattern: a property routine with a selector and a unit conversion") {
        auto run = runModel(R"(
PROCEDURE BrineProp(prop$, fluid$, X, T : v)
  IF (prop$ = 'Specheat') THEN
    v = 3.8 + 0.002*T + 0.5*X
  ELSE
    v = 1.5e-3 - 1e-5*T
  ENDIF
END
brine$ = 'MEG'
X_brine = 0.2
t_su = 5
CALL BRINEPROP('Specheat', brine$, X_brine, t_su : c_p/1000)
CALL BRINEPROP('Dynvisc', brine$, X_brine, t_su : mu*1000)
Q = 0.3 * c_p * 4
)");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["c_p"], WithinAbs(1000.0 * (3.8 + 0.01 + 0.1), 1e-9));
        CHECK_THAT(run["mu"], WithinAbs((1.5e-3 - 5e-5) / 1000.0, 1e-15));
        CHECK_THAT(run["Q"], WithinAbs(0.3 * 3910.0 * 4.0, 1e-6));
    }
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: the unknown of the expression may be an input of the routine",
          "[library-gaps][call-expr-out]") {
    // known output (through the expression), unknown input: the routine is inverted by Newton
    auto run = runModel(kDoubler + "CALL dbl(u : w/1000)\nw = 8000\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["u"], WithinAbs(4.0, 1e-9));

    // nothing determines the expression's variable but the call: it is the unknown of the output
    auto both = runModel(kDoubler + "CALL dbl(u : w/1000)\nu = 4\n");
    REQUIRE(both.solveOk);
    CHECK_THAT(both["w"], WithinAbs(8000.0, 1e-9));
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: coupled block, the Jacobian is analytic", "[library-gaps][call-expr-out]") {
    // c_p depends on the mean temperature, which depends on the outlet temperature, which depends on c_p
    const std::string code = R"(
PROCEDURE cpw(T : c)
  c = 4.0 + 0.5*exp(-T/100)
END
T_in = 10
Q = 5000
m = 0.5
CALL cpw(T_m : c_p/1000)
Q = m * c_p * (T_out - T_in)
T_m = (T_in + T_out)/2
)";
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // fixed point of T_out = T_in + Q/(m c_p(T_m))
    double tOut = 10.0;
    for (int it = 0; it < 200; ++it)
        tOut = 10.0 + 5000.0 / (0.5 * 1000.0 * (4.0 + 0.5 * std::exp(-(10.0 + tOut) / 200.0)));
    const double tMean = (10.0 + tOut) / 2.0;
    CHECK_THAT(run["T_out"], WithinAbs(tOut, 1e-8));
    CHECK_THAT(run["T_m"], WithinAbs(tMean, 1e-8));
    CHECK_THAT(run["c_p"], WithinAbs(1000.0 * (4.0 + 0.5 * std::exp(-tMean / 100.0)), 1e-5));

    // one block of four: c_p, T_m, T_out and the output of the routine
    auto snap = snapshotOfBlock(code, "c_p", {{"c_p", 4000.0}, {"T_m", 20.0}, {"T_out", 30.0},
                                              {"T_in", 10.0}, {"m", 0.5}, {"Q", 5000.0}});
    REQUIRE(snap.vars.size() == 4);
    const int aux = snap.column("__");
    const int cp = snap.column("c_p");
    const int tm = snap.column("T_m");
    const int to = snap.column("T_out");
    REQUIRE(aux >= 0);
    REQUIRE(cp >= 0);
    const int callRow = snap.row("CALL cpw");
    const int exprRow = snap.row("c_p = ");   // c_p/1000 = output is written c_p = output*1000
    const int energyRow = snap.row("Q");
    REQUIRE(callRow >= 0);
    REQUIRE(exprRow >= 0);
    REQUIRE(energyRow >= 0);
    // CALL: residual = (previous output) - c(T_m): d/d(output) = 1, d/dT_m = 0.005 exp(-T_m/100)
    CHECK_THAT(snap.jacobian[callRow][aux], WithinAbs(1.0, 1e-12));
    CHECK_THAT(snap.jacobian[callRow][tm], WithinAbs(0.005 * std::exp(-0.2), 1e-12));
    // c_p = output*1000
    CHECK_THAT(snap.jacobian[exprRow][cp], WithinAbs(1.0, 1e-12));
    CHECK_THAT(snap.jacobian[exprRow][aux], WithinAbs(-1000.0, 1e-12));
    CHECK_THAT(snap.jacobian[exprRow][to], WithinAbs(0.0, 1e-12));
    // Q = m c_p (T_out - T_in): d/dc_p = -m (T_out - T_in) = -10, d/dT_out = -m c_p = -2000
    CHECK_THAT(snap.jacobian[energyRow][cp], WithinAbs(-10.0, 1e-9));
    CHECK_THAT(snap.jacobian[energyRow][to], WithinAbs(-2000.0, 1e-9));
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: an expression that is not solved for the variable stays implicit",
          "[library-gaps][call-expr-out]") {
    // sqrt(c_p) = output: the equation is kept as written, its derivative comes from the AD
    const std::string code = R"(
PROCEDURE cpw(T : c)
  c = 60 + 0.5*exp(-T/100)
END
T_in = 10
Q = 5000
m = 0.5
CALL cpw(T_m : sqrt(c_p))
Q = m * c_p * (T_out - T_in)
T_m = (T_in + T_out)/2
)";
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // c_p = c(T_m)^2
    const double cp = run["c_p"];
    CHECK_THAT(cp, WithinRel(std::pow(60.0 + 0.5 * std::exp(-run["T_m"] / 100.0), 2), 1e-8));
    CHECK_THAT(run["Q"], WithinAbs(0.5 * cp * (run["T_out"] - 10.0), 1e-5));

    auto snap = snapshotOfBlock(code, "c_p", {{"c_p", 3600.0}, {"T_m", 20.0}, {"T_out", 30.0},
                                              {"T_in", 10.0}, {"m", 0.5}, {"Q", 5000.0}});
    REQUIRE(snap.vars.size() == 4);
    const int aux = snap.column("__");
    const int cpCol = snap.column("c_p");
    const int exprRow = snap.row("sqrt(c_p)");
    REQUIRE(aux >= 0);
    REQUIRE(exprRow >= 0);
    CHECK_THAT(snap.jacobian[exprRow][cpCol], WithinAbs(0.5 / 60.0, 1e-12));   // d sqrt(c_p) = 0.5 / sqrt(c_p)
    CHECK_THAT(snap.jacobian[exprRow][aux], WithinAbs(-1.0, 1e-12));
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: forms solved for the variable and forms kept implicit",
          "[library-gaps][call-expr-out]") {
    // every shape gives the value that EES's equation expression = output gives
    auto run = runModel(kDoubler + R"(
CALL dbl(5 : 100 - a)
CALL dbl(5 : b - 100)
CALL dbl(5 : 3*(c + 2) - 1)
CALL dbl(5 : 40/d)
CALL dbl(5 : e*e)
CALL dbl(5 : f/g)
g = 2
CALL dbl(5 : -(h*4 + 1)/2)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinAbs(90.0, 1e-9));
    CHECK_THAT(run["b"], WithinAbs(110.0, 1e-9));
    CHECK_THAT(run["c"], WithinAbs(11.0 / 3.0 - 2.0, 1e-9));
    CHECK_THAT(run["d"], WithinAbs(4.0, 1e-9));
    CHECK_THAT(run["e"] * run["e"], WithinAbs(10.0, 1e-8));
    CHECK_THAT(run["f"], WithinAbs(20.0, 1e-9));
    CHECK_THAT(run["h"], WithinAbs((-20.0 - 1.0) / 4.0, 1e-9));
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: CALL inside a DUPLICATE loop of the main program", "[library-gaps][call-expr-out]") {
    // every iteration has its own output of the routine
    auto run = runModel(kDoubler + R"(
DUPLICATE i = 1, 4
  CALL dbl(i : z[i]/10)
END
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (int i = 1; i <= 4; ++i) CHECK_THAT(run["z[" + std::to_string(i) + "]"], WithinAbs(20.0 * i, 1e-9));
    CHECK(auxiliaryNames(run).size() == 4);

    // the iterator inside the expression, nested loops, a neighbouring plain output
    auto nested = runModel(kDoubler + R"(
DUPLICATE i = 1, 2
  DUPLICATE j = 1, 3
    CALL dbl(10*i + j : a[i,j] + 100*i)
  END
END
)");
    REQUIRE(nested.solveOk);
    CHECK(nested.verified);
    for (int i = 1; i <= 2; ++i)
        for (int j = 1; j <= 3; ++j)
            CHECK_THAT(nested["a[" + std::to_string(i) + "," + std::to_string(j) + "]"],
                       WithinAbs(2.0 * (10 * i + j) - 100.0 * i, 1e-9));
    CHECK(auxiliaryNames(nested).size() == 6);

    // two expression outputs in the same loop body
    auto two = runModel(R"(
PROCEDURE pr(a : s, t)
  s = a + 1
  t = a * 3
END
DUPLICATE k = 1, 3
  CALL pr(k : p[k]*2, q[k] - 1)
END
)");
    REQUIRE(two.solveOk);
    for (int k = 1; k <= 3; ++k) {
        CHECK_THAT(two["p[" + std::to_string(k) + "]"], WithinAbs((k + 1) / 2.0, 1e-9));
        CHECK_THAT(two["q[" + std::to_string(k) + "]"], WithinAbs(3.0 * k + 1.0, 1e-9));
    }
    CHECK(auxiliaryNames(two).size() == 6);
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: in an INTEGRAL model and in the parametric solver", "[library-gaps][call-expr-out]") {
    // the time-march solves the algebraic part at every step: the expression output is re-lowered there
    auto run = runIntegral(kDoubler + "CALL dbl(t : c/1000)\ny = integral(dydt, t, 0, 1, 0.1)\ndydt = c\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run.vars.at("y"), WithinAbs(1000.0, 1e-6));   // dy/dt = 2000 t

    // a Solver built by hand from the parsed program (what the parametric study does)
    EESParser parser;
    auto parsed = parser.parse(kDoubler + "CALL dbl(a : w/1000)\nz = 3*w\na = 4\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    auto analysis = StructuralAnalyzer::analyze(ir);
    REQUIRE(analysis.success);
    SolverOptions options;
    Solver solver(ir, analysis, options.coolpropConfig);
    auto result = solver.solve(options, false);
    REQUIRE(result.success);
    CHECK_THAT(result.variables.at("w"), WithinAbs(8000.0, 1e-9));
    CHECK_THAT(result.variables.at("z"), WithinAbs(24000.0, 1e-9));
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: the auxiliary variable is documented and unique", "[library-gaps][call-expr-out]") {
    // __out<position>_<procedure>_L<line>: the raw output of the routine, in the units of the routine
    auto run = runModel("PROCEDURE pg(a:y)\ny=a*2\nEND\nCALL pg(21:z/1000)\n");
    REQUIRE(run.solveOk);
    REQUIRE(auxiliaryNames(run).size() == 1);
    CHECK(auxiliaryNames(run)[0] == "__out1_pg_L4");
    CHECK_THAT(run["__out1_pg_L4"], WithinAbs(42.0, 1e-9));

    // models without an expression output have no such variable
    auto plain = runModel(kDoubler + "CALL dbl(4 : r)\n");
    REQUIRE(plain.solveOk);
    CHECK(auxiliaryNames(plain).empty());
    CHECK(plain.vars.size() == 1);

    // an equal name written by the user is never merged with the auxiliary variable
    auto clash = runModel("PROCEDURE pg(a:y)\ny=a*2\nEND\n__out1_pg_L5 = 7\nCALL pg(21:z/1000)\n");
    REQUIRE(clash.solveOk);
    CHECK_THAT(clash["__out1_pg_L5"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(clash["z"], WithinAbs(42000.0, 1e-6));
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: the debug folder and the solution file", "[library-gaps][call-expr-out]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_library_gaps_exprout_debug";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string code = "PROCEDURE pg(a:y)\ny=a*2\nEND\nCALL pg(21:z/1000)\nw=z\n";
    const fs::path file = dir / "model.eescode";
    { std::ofstream f(file); f << code; }
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    REQUIRE(runner.run(options));
    runner.generateDebugOutput((dir / "debug").string(), code);
    auto slurp = [&](const char* name) {
        std::ifstream in(dir / "debug" / name);
        REQUIRE(in.good());
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };
    // the equation that ties the expression to the output is listed (z/1000 = output, solved for z)
    const std::string equations = slurp("equations.md");
    CHECK(equations.find("z = (__out1_pg_L4 * 1000") != std::string::npos);
    const std::string variables = slurp("variables.md");
    CHECK(variables.find("__out1_pg_L4") != std::string::npos);
    fs::remove_all(dir);
}

TEST_CASE("CS-GAP-CALL-EXPR-OUT: negative tests", "[library-gaps][call-expr-out]") {
    SECTION("an unparsable output is still a parse error") {
        CHECK(parseErrors(kDoubler + "CALL dbl(1 : r, 3 4)\n") ==
              std::vector<std::string>{"Could not parse output 2 of 'dbl': '3 4'"});
    }

    SECTION("a CALL inside a FUNCTION or PROCEDURE body: statements run in order, nothing to solve") {
        // EES: outputs of a CALL in a routine are variables; the expression cannot be inverted
        const auto inProc = parseErrors(kDoubler +
            "PROCEDURE q(a : r)\n  CALL dbl(a : r/1000)\nEND\nCALL q(1 : y)\n");
        REQUIRE(inProc.size() == 1);
        CHECK(anyContains(inProc, "Output 1 of 'dbl' is the expression 'r/1000'"));
        CHECK(anyContains(inProc, "FUNCTION or PROCEDURE"));
        const auto inFunc = parseErrors(kDoubler +
            "FUNCTION f(a)\n  CALL dbl(a : -r)\n  f = r\nEND\ny = f(1)\n");
        CHECK(anyContains(inFunc, "Output 1 of 'dbl' is the expression '-r'"));
        // inside an IF of a body too
        const auto inIf = parseErrors(kDoubler +
            "PROCEDURE q(a : r)\n  IF (a > 0) THEN\n    CALL dbl(a : r + 1)\n  ENDIF\nEND\nCALL q(1 : y)\n");
        CHECK(anyContains(inIf, "Output 1 of 'dbl' is the expression 'r + 1'"));
        // a variable output in a body is fine
        auto fine = runModel(kDoubler + "PROCEDURE q(a : r)\n  CALL dbl(a : t)\n  r = t + 1\nEND\nCALL q(1 : y)\n");
        REQUIRE(fine.solveOk);
        CHECK_THAT(fine["y"], WithinAbs(3.0, 1e-12));
    }

    SECTION("too many outputs are still counted, with the expressions") {
        auto run = runModel(kDoubler + "CALL dbl(1 : z/1000, w)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "returns 1 outputs, but 2 were provided"));
    }

    SECTION("an expression output that contradicts a given value is reported") {
        // 2*5 = 10 but the output is required to be 11
        auto run = runModel(kDoubler + "CALL dbl(5 : 11)\n");
        CHECK_FALSE(run.solveOk);
        auto over = runModel(kDoubler + "CALL dbl(5 : z/1000)\nz = 1\n");
        CHECK_FALSE(over.solveOk);
    }
}

// ============================================================================
// CS-BUG-STRING-CALL-OUT: a string produced only by a CALL, used by later statements
// ============================================================================

namespace {

/// `gen` computes a string from a number (no literal assignment of the string in the main program).
const std::string kStringGen =
    "PROCEDURE gen(x : r$)\n"
    "  IF (x > 1) THEN\n"
    "    r$ := 'BIG'\n"
    "  ENDIF\n"
    "  IF (x <= 1) THEN\n"
    "    r$ := 'SMALL'\n"
    "  ENDIF\n"
    "END\n"
    "PROCEDURE use(s$ : y)\n"
    "  IF (s$ = 'BIG') THEN\n"
    "    y = 1\n"
    "  ELSE\n"
    "    y = 2\n"
    "  ENDIF\n"
    "END\n";

}  // namespace

TEST_CASE("CS-BUG-STRING-CALL-OUT: a string computed by a CALL reaches the next CALL", "[library-gaps][string-call-out]") {
    // was: "Block 1 (size 1, vars: z) failed: EvaluationError - [Newton] String variable not found: r$"
    auto big = runModel(kStringGen + "CALL gen(2 : r$)\nCALL use(r$ : z)\n");
    REQUIRE(big.solveOk);
    CHECK(big.verified);
    CHECK_THAT(big["z"], WithinAbs(1.0, 1e-12));
    CHECK(big.strings.at("r$") == "BIG");

    auto small = runModel(kStringGen + "CALL gen(1 : r$)\nCALL use(r$ : z)\n");
    REQUIRE(small.solveOk);
    CHECK_THAT(small["z"], WithinAbs(2.0, 1e-12));
    CHECK(small.strings.at("r$") == "SMALL");

    // the order of the lines does not matter: the structural analysis orders the blocks
    auto reversed = runModel(kStringGen + "CALL use(r$ : z)\nCALL gen(2 : r$)\n");
    REQUIRE(reversed.solveOk);
    CHECK_THAT(reversed["z"], WithinAbs(1.0, 1e-12));

    // the number given to the first routine is itself computed
    auto chained = runModel(kStringGen + "n = 3 - 1\nCALL gen(n : r$)\nCALL use(r$ : z)\nw = 10*z\n");
    REQUIRE(chained.solveOk);
    CHECK(chained.verified);
    CHECK_THAT(chained["w"], WithinAbs(10.0, 1e-12));
}

TEST_CASE("CS-BUG-STRING-CALL-OUT: the register reproducer", "[library-gaps][string-call-out]") {
    // The reproducer of the register writes into a lookup table (CS-BUG-LOOKUP-WRITE, a separate bug).
    // Whatever becomes of those statements, the string must be found.
    auto run = runModel(R"(
PROCEDURE gen(x:r$)
IF (x>1) THEN
r$:='BIG'
ENDIF
IF (x<=1) THEN
r$:='SMALL'
ENDIF
END
PROCEDURE wr(i,regime$:res$)
lookup('noTable',i,1)=i
lookup$('noTable',i,2)=regime$
res$='OK'
END
i=1
CALL gen(2:r$)
CALL wr(i,r$:res$)
)");
    CHECK_FALSE(anyContains({run.message}, "String variable not found"));
    // with the string given by a literal the file solves, as before
    auto literal = runModel(R"(
PROCEDURE wr(i,regime$:res$)
res$='OK'
END
i=1
r$='ON'
CALL wr(i,r$:res$)
)");
    CHECK(literal.solveOk);
}

TEST_CASE("CS-BUG-STRING-CALL-OUT: other uses of the string", "[library-gaps][string-call-out]") {
    SECTION("several string outputs, mixed with numeric ones") {
        auto run = runModel(R"(
PROCEDURE props(k : name$, v, kind$)
  IF (k = 1) THEN
    name$ := 'Water'
    kind$ := 'liquid'
  ELSE
    name$ := 'R134a'
    kind$ := 'refrigerant'
  ENDIF
  v := 10*k
END
PROCEDURE isliq(kind$ : f)
  IF (kind$ = 'liquid') THEN
    f = 1
  ELSE
    f = 0
  ENDIF
END
CALL props(1 : n$, a, t$)
CALL props(2 : n2$, b, t2$)
CALL isliq(t$ : f1)
CALL isliq(t2$ : f2)
)");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK(run.strings.at("n$") == "Water");
        CHECK(run.strings.at("t$") == "liquid");
        CHECK(run.strings.at("n2$") == "R134a");
        CHECK(run.strings.at("t2$") == "refrigerant");
        CHECK_THAT(run["a"], WithinAbs(10.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(20.0, 1e-12));
        CHECK_THAT(run["f1"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["f2"], WithinAbs(0.0, 1e-12));
    }

    SECTION("the string is the fluid of a property call of the main program") {
        auto run = runModel(R"(
PROCEDURE pick(k : fl$)
  IF (k = 1) THEN
    fl$ := 'Water'
  ELSE
    fl$ := 'R134a'
  ENDIF
END
CALL pick(1 : fluid$)
h = enthalpy(fluid$, T=100, P=1e5)
)");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["h"], WithinRel(CoolProp::PropsSI("H", "T", 373.15, "P", 1e5, "Water"), 1e-6));
    }

    SECTION("the string is the fluid of a property call inside another routine") {
        auto run = runModel(R"(
PROCEDURE pick(k : fl$)
  fl$ := 'Water'
END
PROCEDURE hh(fl$, T : h)
  h := enthalpy(fl$, T=T, P=1e5)
END
CALL pick(1 : fluid$)
CALL hh(fluid$, 100 : h)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["h"], WithinRel(CoolProp::PropsSI("H", "T", 373.15, "P", 1e5, "Water"), 1e-6));
    }

    SECTION("a coupled numeric block reads the string") {
        auto run = runModel(kStringGen + R"(
PROCEDURE scaled(s$, a : y)
  IF (s$ = 'BIG') THEN
    y = 2*a
  ELSE
    y = a/2
  ENDIF
END
CALL gen(2 : r$)
CALL scaled(r$, a : y)
a + y = 9
)");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        // 'BIG': y = 2 a, so a = 3 and y = 6 (with 'SMALL' it would be a = 6)
        CHECK_THAT(run["a"], WithinAbs(3.0, 1e-9));
        CHECK_THAT(run["y"], WithinAbs(6.0, 1e-9));
    }

    SECTION("a string given as a literal is untouched") {
        auto run = runModel(kStringGen + "r$ = 'BIG'\nCALL use(r$ : z)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["z"], WithinAbs(1.0, 1e-12));
        CHECK(run.strings.at("r$") == "BIG");
    }

    SECTION("the string is computed inside the coupled block that reads it") {
        // a is the input of the routine that sets the string, and the string changes the equation that
        // fixes a: one block of three unknowns (a, y, r$), whatever the order of the lines
        const std::string scaled =
            "PROCEDURE scaled(s$, a : y)\n  IF (s$ = 'BIG') THEN\n    y = 2*a\n  ELSE\n    y = a/2\n  ENDIF\nEND\n";
        for (const char* order : {"CALL gen(a : r$)\nCALL scaled(r$, a : y)\na + y = 9\n",
                                  "a + y = 9\nCALL scaled(r$, a : y)\nCALL gen(a : r$)\n"}) {
            auto run = runModel(kStringGen + scaled + order);
            REQUIRE(run.solveOk);
            CHECK(run.verified);
            CHECK_THAT(run["a"], WithinAbs(3.0, 1e-9));   // 'BIG' (a > 1): y = 2 a, a + 2 a = 9
            CHECK_THAT(run["y"], WithinAbs(6.0, 1e-9));
            CHECK(run.strings.at("r$") == "BIG");
        }
        // a chain of two string producers in the block
        auto chain = runModel(kStringGen + scaled +
            "PROCEDURE rename(s$ : t$)\n  t$ := s$\nEND\n"
            "CALL scaled(q$, a : y)\nCALL rename(r$ : q$)\na + y = 9\nCALL gen(a : r$)\n");
        REQUIRE(chain.solveOk);
        CHECK(chain.verified);
        CHECK_THAT(chain["a"], WithinAbs(3.0, 1e-9));
        CHECK(chain.strings.at("q$") == "BIG");
    }
}

TEST_CASE("CS-BUG-STRING-CALL-OUT: the string is in the solution and the debug folder", "[library-gaps][string-call-out]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_library_gaps_string_debug";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string code = kStringGen + "CALL gen(2 : r$)\nCALL use(r$ : z)\n";
    const fs::path file = dir / "model.eescode";
    { std::ofstream f(file); f << code; }
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    REQUIRE(runner.run(options));
    const auto& res = runner.getSolveResult();
    CHECK(res.stringVariables.at("r$") == "BIG");
    // the post-solve check evaluates the second CALL with the string of the first
    auto check = checkSolution(runner.getIR(), res.variables, res.stringVariables, options.coolpropConfig, 1e-3);
    CHECK(check.allSatisfied);
    runner.generateDebugOutput((dir / "debug").string(), code);
    std::ifstream in(dir / "debug" / "variables.md");
    REQUIRE(in.good());
    const std::string variables((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    CHECK(variables.find("`r$`") != std::string::npos);
    fs::remove_all(dir);
}

TEST_CASE("CS-BUG-STRING-CALL-OUT: negative tests", "[library-gaps][string-call-out]") {
    SECTION("a string that nothing computes is reported with its name") {
        // never$ is an unknown without equation: the system is not square
        auto run = runModel(kStringGen + "CALL use(never$ : z)\n");
        CHECK_FALSE(run.solveOk);
    }

    SECTION("a routine that does not set its string output") {
        auto run = runModel("PROCEDURE lazy(x : r$)\n  y = x\nEND\n" + std::string(
            "PROCEDURE use(s$ : y)\n  y = 1\nEND\nCALL lazy(1 : r$)\nCALL use(r$ : z)\n"));
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "r$"));
    }
}

// ============================================================================
// Characters of variable names: CS-GAP-NAME-SYMBOL (%), CS-GAP-NAME-PIPE (|),
// CS-GAP-ACCENTED-IDENT (accented letters, °) and CS-BUG-HASH-IDENT (# inside a name)
//
// EES help, "Variable names": a name starts with a letter and holds any keyboard character
// except ( ) * / + - ^ { } ; (older manual: also ' | : "). `$` ends a string variable,
// `[ ]` index an array. The valid-EES evidence of each character is in the register.
// ============================================================================

namespace {

/// The variables of `run` whose name is not a CoolSolve auxiliary (`__...`).
std::vector<std::string> userVariables(const ModelRun& run) {
    std::vector<std::string> names;
    for (const auto& [name, value] : run.vars) {
        (void)value;
        if (name.rfind("__", 0) != 0) names.push_back(name);
    }
    return names;
}

/// Run the report generator on a solved model and return the LaTeX source.
std::string latexOf(const std::string& code) {
    static int counter = 0;
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_name_chars";
    fs::create_directories(dir);
    const fs::path file = dir / ("latex_" + std::to_string(counter++) + ".eescode");
    { std::ofstream f(file); f << code; }
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    options.tolerance = 1e-9;
    const bool ok = runner.run(options, false);
    fs::remove(file);
    REQUIRE(ok);
    return runner.generateLatexReportContent("name_chars");
}

/// True when pdflatex is installed.
bool havePdflatex() { return std::system("pdflatex --version > /dev/null 2>&1") == 0; }

/// Compile `tex` with pdflatex; true when the run ends without error.
bool pdflatexCompiles(const std::string& tex, std::string* log = nullptr) {
    static int counter = 0;
    const fs::path dir = fs::temp_directory_path() / ("coolsolve_test_name_chars_tex_" + std::to_string(counter++));
    fs::remove_all(dir);
    fs::create_directories(dir);
    { std::ofstream f(dir / "r.tex"); f << tex; }
    const std::string cmd = "cd '" + dir.string() + "' && pdflatex -interaction=nonstopmode -halt-on-error r.tex > out.txt 2>&1";
    const bool ok = std::system(cmd.c_str()) == 0;
    if (log) {
        std::ifstream in(dir / "r.log");
        *log = std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    fs::remove_all(dir);
    return ok;
}

}  // namespace

// ---------------------------------------------------------------------------
// CS-GAP-NAME-SYMBOL: `%` in a name
// ---------------------------------------------------------------------------

TEST_CASE("CS-GAP-NAME-SYMBOL: register reproducer, the fuel mass fractions C%, H%, O%", "[library-gaps][name-symbol]") {
    // valid EES (orc_complex.EES, EES 9.920): C% = 0.5, H% = 0.06, O% = 1-C%-H%, m = C%/12 stored as 0.04167
    // was: "Parse failed: Line 1 … Line 5: Could not parse line"
    auto run = runModel("C% = 0.5\nH% = 0.06\nm_boil = C%/12\nn_boil = H%/1\nO% = 1-C%-H%\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["m_boil"], WithinAbs(0.5 / 12.0, 1e-15));
    CHECK_THAT(run["n_boil"], WithinAbs(0.06, 1e-15));
    CHECK_THAT(run["O%"], WithinAbs(0.44, 1e-15));
    CHECK_THAT(run["C%"], WithinAbs(0.5, 1e-15));
    CHECK_THAT(run["m_boil"], WithinAbs(0.04167, 1e-5));   // the value EES stores
}

TEST_CASE("CS-GAP-NAME-SYMBOL: the % in every position of the model", "[library-gaps][name-symbol]") {
    SECTION("start, middle, end of a name, in expressions and in function arguments") {
        auto run = runModel(
            "C% = 86.9\nH% = 100 - C%\nO% = 0\nS% = 0\n"
            "HHV_kJ = 4.186*(81*C%+300*(H%-O%/8)+25*S%)\n"      // the correlation of CSL-0150
            "Q_dot_2_% = 100*3/12\n"                             // trailing % (CSL-0152)
            "m = (C%/100)/12.011\n"
            "y = max(C%, H%) + min(C%/2, 7)\n"
            "p%q = sqrt(C%)\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["HHV_kJ"], WithinRel(4.186 * (81 * 86.9 + 300 * 13.1), 1e-12));
        CHECK_THAT(run["Q_dot_2_%"], WithinAbs(25.0, 1e-12));
        CHECK_THAT(run["m"], WithinRel(0.869 / 12.011, 1e-12));
        CHECK_THAT(run["y"], WithinRel(86.9 + 7.0, 1e-12));
        CHECK_THAT(run["p%q"], WithinRel(std::sqrt(86.9), 1e-12));
    }

    SECTION("FUNCTION and PROCEDURE: parameters, locals and CALL outputs") {
        auto run = runModel(
            "FUNCTION hhv(C%, H%)\n  S% = 0\n  hhv = 4.186*(81*C%+300*H%+25*S%)\nEND\n"
            "PROCEDURE comp(C%, H% : m%, n%)\n  m% = C%/12\n  n% = H%/1\nEND\n"
            "C% = 86.9\nH% = 100 - C%\n"
            "x = hhv(C%, H%)\n"
            "CALL comp(C%, H% : mm%, nn%)\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinRel(4.186 * (81 * 86.9 + 300 * 13.1), 1e-12));
        CHECK_THAT(run["mm%"], WithinRel(86.9 / 12.0, 1e-12));
        CHECK_THAT(run["nn%"], WithinRel(13.1, 1e-12));
        // the locals of a routine are not variables of the model
        CHECK(run.vars.count("S%") == 0);
        CHECK(run.vars.count("m%") == 0);
    }

    SECTION("IF conditions and loops inside a routine") {
        auto run = runModel(
            "FUNCTION g(a%)\n"
            "  IF (a% > 50) THEN g = 1 ELSE g = 2\n"
            "END\n"
            "PROCEDURE cnt(n% : s%)\n"
            "  s% = 0\n"
            "  i% = 0\n"
            "  REPEAT\n"
            "    i% = i% + 1\n"
            "    s% = s% + i%\n"
            "  UNTIL (i% >= n%)\n"
            "END\n"
            "a = g(86.9)\nb = g(10)\nCALL cnt(4 : c%)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["c%"], WithinAbs(10.0, 1e-12));
    }

    SECTION("array elements X%[i], DUPLICATE and sums") {
        auto run = runModel(
            "C% = 0.5\n"
            "DUPLICATE i = 1, 3\n  X%[i] = i*C%\nEND\n"
            "S = X%[1] + X%[2] + X%[3]\n"
            "Y%[2] = 5\nz = Y%[2] * 2\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["X%[2]"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["S"], WithinAbs(3.0, 1e-12));
        CHECK_THAT(run["z"], WithinAbs(10.0, 1e-12));
    }

    SECTION("a string variable and a $IF directive") {
        auto run = runModel(
            "fuel%$ = 'CH4'\n"
            "$IF fuel%$='ch4'\nx = 1\n$ELSE\nx = 2\n$ENDIF\n"
            "$IFNOT fuel%$='CH4'\ny = 1\n$ELSE\ny = 2\n$ENDIF\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["y"], WithinAbs(2.0, 1e-12));
        CHECK(run.strings.at("fuel%$") == "CH4");
    }

    SECTION("thermophysical function with a name=value argument") {
        auto run = runModel("T% = 100\nP% = pressure(Water, T=T%, x=1)\n");   // °C, the default unit system
        REQUIRE(run.solveOk);
        CHECK_THAT(run["P%"], WithinRel(CoolProp::PropsSI("P", "T", 373.15, "Q", 1, "Water"), 1e-6));
    }
}

TEST_CASE("CS-GAP-NAME-SYMBOL: state of an INTEGRAL model, the trajectory and its CSV", "[library-gaps][name-symbol]") {
    auto run = runIntegral(
        "y% = 0 + integral(dydt%, t, 0, 2, 0.01)\n"
        "dydt% = 3*t\n"
        "$IntegralTable t:1 y% dydt%\n",
        [](SolverOptions& o) { o.integralMethod = "RK4"; });
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y%"], WithinAbs(6.0, 1e-9));
    REQUIRE(run.table.numRows() == 3);
    CHECK(run.table.columns() == std::vector<std::string>{"t", "y%", "dydt%"});
    CHECK_THAT(run.cell("y%", 1), WithinAbs(1.5, 1e-9));
    CHECK(run.table.toCSV().rfind("t,y%,dydt%", 0) == 0);
}

TEST_CASE("CS-GAP-NAME-SYMBOL: names are case-insensitive and the debug folder lists them", "[library-gaps][name-symbol]") {
    auto run = runModel("C% = 0.5\ny = c%*2\nc_x% = 3\nz = C_X%\n");
    REQUIRE(run.solveOk);
    CHECK(userVariables(run).size() == 4);      // C%, y, c_x%, z: c% and C% are one variable
    CHECK_THAT(run["y"], WithinAbs(1.0, 1e-15));
    CHECK_THAT(run["z"], WithinAbs(3.0, 1e-15));

    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_name_symbol_debug";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string code = "C% = 0.5\nm_boil = C%/12\n";
    { std::ofstream f(dir / "model.eescode"); f << code; }
    CoolSolveRunner runner((dir / "model.eescode").string());
    SolverOptions options;
    REQUIRE(runner.run(options));
    runner.generateDebugOutput((dir / "debug").string(), code);
    auto slurp = [&](const char* name) {
        std::ifstream in(dir / "debug" / name);
        REQUIRE(in.good());
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };
    CHECK(slurp("variables.md").find("`C%`") != std::string::npos);
    CHECK(slurp("equations.md").find("m_boil = (C% / 12") != std::string::npos);
    CHECK(slurp("incidence.md").find("`C%`") != std::string::npos);
    CHECK(slurp("variables.csv").find("\"C%\"") != std::string::npos);
    CHECK(slurp("analysis.json").find("\"C%\"") != std::string::npos);
    fs::remove_all(dir);
}

TEST_CASE("CS-GAP-NAME-SYMBOL: .initials and .sol files with a % in the name", "[library-gaps][name-symbol]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_name_symbol_files";
    fs::remove_all(dir);
    fs::create_directories(dir);
    EESParser parser;
    auto parsed = parser.parse("C% = 0.5\nx% = 2*C%\nc2% = x%^2\n");
    REQUIRE(parsed.success);
    SECTION("initials: values, units and comments") {
        { std::ofstream f(dir / "m.initials"); f << "C% = 0.4 \"kg/kg\"\n# a comment line\nx% = 3 # trailing comment\n"; }
        IR ir = IR::fromAST(parsed.program);
        CHECK(ir.loadInitialsFromFile((dir / "m.initials").string()) == 2);
        REQUIRE(ir.getVariable("C%") != nullptr);
        CHECK_THAT(*ir.getVariable("C%")->solutionValue, WithinAbs(0.4, 1e-15));
        CHECK(ir.getVariable("C%")->units == "kg/kg");
        CHECK_THAT(*ir.getVariable("x%")->solutionValue, WithinAbs(3.0, 1e-15));
        CHECK_FALSE(ir.getVariable("c2%")->solutionValue.has_value());
    }
    SECTION("the .sol written by a solve is read back") {
        auto run = runModel("C% = 0.5\nx% = 2*C%\nc2% = x%^2\n");
        REQUIRE(run.solveOk);
        std::ofstream f(dir / "m.sol");
        f << std::scientific << std::setprecision(12);
        for (const auto& [name, value] : run.vars) f << name << " = " << value << "\n";
        f << "\n# IntegralTable\nt,x%\n0,1\n";
        f.close();
        IR ir = IR::fromAST(parsed.program);
        CHECK(ir.loadSolutionFromFile((dir / "m.sol").string()) == 3);
        CHECK_THAT(*ir.getVariable("c2%")->solutionValue, WithinAbs(1.0, 1e-12));
    }
    fs::remove_all(dir);
}

TEST_CASE("CS-GAP-NAME-SYMBOL: LaTeX output escapes the % (a comment character of LaTeX)", "[library-gaps][name-symbol][latex]") {
    CHECK(variableToLatex("C%") == "C\\%");
    CHECK(variableToLatex("Q_dot_2_%") == "\\dot{Q}_{2,\\%}");
    CHECK(variableToLatex("X%[2]") == "X\\%_{2}");
    const std::string tex = latexOf("C% = 0.5\nm_boil = C%/12\nQ_dot_2_% = 100*3/12\n");
    // no raw % starts a comment in the middle of a math formula or of a table cell
    CHECK(tex.find("C\\%") != std::string::npos);
    CHECK(tex.find("$C% ") == std::string::npos);
    CHECK(tex.find("\\frac{C\\%}{12") != std::string::npos);
    if (havePdflatex()) {
        std::string log;
        INFO(log);
        CHECK(pdflatexCompiles(tex, &log));
    }
}

TEST_CASE("CS-GAP-NAME-SYMBOL: negative tests, characters that stay invalid give a clear error", "[library-gaps][name-symbol]") {
    SECTION("a name starts with a letter: not with %, # or |") {
        const auto e1 = parseErrors("%x = 3\n");
        REQUIRE_FALSE(e1.empty());
        CHECK(anyContains(e1, "Invalid variable name '%x'"));
        CHECK(anyContains(parseErrors("x = |a + 1\n"), "Invalid variable name '|a'"));
        CHECK(anyContains(parseErrors("x = 2 + #a\n"), "Invalid variable name '#a'"));
        CHECK_FALSE(parseErrors("1x = 3\n").empty());   // a name does not start with a digit
    }
    SECTION("characters EES does not take in names, or that CoolSolve needs as syntax") {
        CHECK(anyContains(parseErrors("x@y = 3\n"), "Invalid variable name 'x@y': the character '@'"));
        CHECK(anyContains(parseErrors("a&b = 3\n"), "the character '&'"));
        CHECK(anyContains(parseErrors("a~b = 3\n"), "the character '~'"));
        CHECK(anyContains(parseErrors("a?b = 3\n"), "the character '?'"));
        CHECK(anyContains(parseErrors("x = a@b + 1\n"), "Invalid variable name 'a@b'"));
        // inside a call, a CALL argument and a CALL output the same message
        CHECK(anyContains(parseErrors("x = max(a@b, 1)\n"), "Invalid variable name 'a@b'"));
        CHECK(anyContains(parseErrors("PROCEDURE p(a : b)\nb = a\nEND\nCALL p(1 : r@s)\n"), "Invalid variable name 'r@s'"));
        CHECK(anyContains(parseErrors("PROCEDURE p(a : b)\nb = a\nEND\nCALL p(r@s : y)\n"), "Invalid variable name 'r@s'"));
        // a % inside a function body is as good as outside, a bad one is as loud as outside
        CHECK(anyContains(parseErrors("FUNCTION f(a)\n  f = a@b\nEND\ny = f(1)\n"), "Invalid variable name 'a@b'"));
    }
    SECTION("the quote opens a string and the dot is a decimal point: neither can be part of a name") {
        CHECK_FALSE(parseErrors("a'b = 3\n").empty());
        CHECK(anyContains(parseErrors("a.b = 3\n"), "Invalid variable name 'a.b': the character '.'"));
    }
    SECTION("operators are still operators, a name stops where one starts") {
        auto run = runModel("a% = 6\nb% = 3\nx = a%+b%\ny = a%-b%\nz = a%*b%\nw = a%/b%\nv = a%^2\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(9.0, 1e-12));
        CHECK_THAT(run["y"], WithinAbs(3.0, 1e-12));
        CHECK_THAT(run["z"], WithinAbs(18.0, 1e-12));
        CHECK_THAT(run["w"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["v"], WithinAbs(36.0, 1e-12));
    }
}

// ---------------------------------------------------------------------------
// CS-GAP-NAME-PIPE: `|` in a name
// ---------------------------------------------------------------------------

TEST_CASE("CS-GAP-NAME-PIPE: register reproducer", "[library-gaps][name-pipe]") {
    // valid EES (the K* coefficient of the cooling-coil bank): was "Line 1/2: Could not parse line"
    auto run = runModel("K|star_a = 106\nK|star_b = 2*K|star_a\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["K|star_b"], WithinAbs(212.0, 1e-12));
    CHECK(run.warnings.empty());   // was: "Unknown function 'K|star_r_n/K|star_r*'"
}

TEST_CASE("CS-GAP-NAME-PIPE: the cooling-coil pattern (CSL-0074, CSL-0163)", "[library-gaps][name-pipe]") {
    // R_r_coil = R_r_coil_n*K|star_r_n/K|star_r*(M_dot_r_coil_n / M_dot_r_coil)^0.8, K|star_r from the properties
    auto run = runModel(
        "K|star_r_n = 106\n"
        "mu_r_coil = 2.1e-4\nk_r_coil = 0.085\nc_p_r_coil = 1300\nn = 0.8\nm = 0.2\n"
        "K|star_r=mu_r_coil^(n-m)*k_r_coil^(1-n)*c_p_r_coil^n\n"
        "R_r_coil_n = 0.0125\nM_dot_r_coil_n = 0.04\nM_dot_r_coil = 0.03\n"
        "R_r_coil = R_r_coil_n*K|star_r_n/K|star_r*(M_dot_r_coil_n / M_dot_r_coil)^0.8  \"turbulent flow in circular tube\"\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK(run.warnings.empty());
    const double kstar = std::pow(2.1e-4, 0.6) * std::pow(0.085, 0.2) * std::pow(1300.0, 0.8);
    CHECK_THAT(run["K|star_r"], WithinRel(kstar, 1e-12));
    CHECK_THAT(run["R_r_coil"], WithinRel(0.0125 * 106.0 / kstar * std::pow(0.04 / 0.03, 0.8), 1e-12));
}

TEST_CASE("CS-GAP-NAME-PIPE: the compressor pattern, | in the values of name=value arguments (CSL-0116)", "[library-gaps][name-pipe]") {
    // PRESSURE(fluid$, v=v_in_cp|star, u=u_in_cp|star), enthalpy(...) and implicit equations between starred variables
    auto run = runModel(
        "fluid$ = 'R134a'\n"
        "T_in = 27\nx_in = 1\n"
        "v_in_cp|star = volume(fluid$, T=T_in, x=x_in)\n"
        "u_in_cp|star = intenergy(fluid$, T=T_in, x=x_in)\n"
        "p_in_cp|star = PRESSURE(fluid$, v=v_in_cp|star, u=u_in_cp|star)\n"
        "h_in_cp|star = enthalpy(fluid$, v=v_in_cp|star, u=u_in_cp|star)\n"
        "w_in2_cp_v|star = v_in_cp|star * (p_in_cp|star - 100000)\n"
        "{w_in2_cp=if(p_in_cp;p_ex2_cp;w_in2_cp|star;w_in2_cp|star;w_in2_cp_v)}\n");   // a comment with | is a comment
    REQUIRE(run.solveOk);
    CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
    const double p = CoolProp::PropsSI("P", "T", 300.15, "Q", 1, "R134a");
    CHECK_THAT(run["p_in_cp|star"], WithinRel(p, 1e-4));
    CHECK_THAT(run["w_in2_cp_v|star"], WithinRel(run["v_in_cp|star"] * (run["p_in_cp|star"] - 100000.0), 1e-12));
}

TEST_CASE("CS-GAP-NAME-PIPE: the | in every position of the model", "[library-gaps][name-pipe]") {
    SECTION("routines: parameters, locals, outputs, IF conditions") {
        auto run = runModel(
            "FUNCTION f(a|o, b|star)\n  c|loc = a|o*2\n  IF (b|star > 1) THEN f = c|loc + b|star ELSE f = c|loc\nEND\n"
            "PROCEDURE p(a|o : r|o, s|star)\n  r|o = a|o + 1\n  s|star = a|o*3\nEND\n"
            "x = f(2, 3)\ny = f(2, 0)\nCALL p(4 : q|o, w|star)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(7.0, 1e-12));
        CHECK_THAT(run["y"], WithinAbs(4.0, 1e-12));
        CHECK_THAT(run["q|o"], WithinAbs(5.0, 1e-12));
        CHECK_THAT(run["w|star"], WithinAbs(12.0, 1e-12));
        CHECK(run.vars.count("c|loc") == 0);
    }
    SECTION("arrays and DUPLICATE") {
        auto run = runModel("DUPLICATE i = 1, 3\n  G|o[i] = i^2\nEND\ns = G|o[1] + G|o[2] + G|o[3]\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["s"], WithinAbs(14.0, 1e-12));
        CHECK_THAT(run["G|o[3]"], WithinAbs(9.0, 1e-12));
    }
    SECTION("string variables and $IF") {
        auto run = runModel("mode|a$ = 'ON'\n$IF mode|a$='on'\nx = 1\n$ELSE\nx = 2\n$ENDIF\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));
    }
    SECTION("the state of an INTEGRAL model") {
        auto run = runIntegral(
            "T|star = 1 + integral(dTdt|star, t, 0, 1, 0.01)\n"
            "dTdt|star = -T|star\n"
            "$IntegralTable t:0.5 T|star\n",
            [](SolverOptions& o) { o.integralMethod = "RK4"; });
        REQUIRE(run.solveOk);
        CHECK_THAT(run["T|star"], WithinAbs(std::exp(-1.0), 1e-8));
        CHECK(run.table.columns() == std::vector<std::string>{"t", "T|star"});
    }
    SECTION("an implicit equation between starred variables is solved") {
        auto run = runModel("x|star^3 + x|star = 10\ny|star = x|star*2\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x|star"], WithinAbs(2.0, 1e-9));
    }
    SECTION("names are case-insensitive") {
        auto run = runModel("K|STAR = 4\ny = k|Star*2\n");
        REQUIRE(run.solveOk);
        CHECK(userVariables(run).size() == 2);
        CHECK_THAT(run["y"], WithinAbs(8.0, 1e-12));
    }
}

TEST_CASE("CS-GAP-NAME-PIPE: LaTeX output, the | starts a superscript as in EES", "[library-gaps][name-pipe][latex]") {
    CHECK(variableToLatex("G|o") == "G^{o}");
    CHECK(variableToLatex("K|star_r_n") == "K^{star}_{r,n}");
    CHECK(variableToLatex("p_in_cp|star") == "p_{in,cp^{star}}");
    CHECK(variableToLatex("C|o_m") == "C^{o}_{m}");
    CHECK(variableToLatex("a||b") == "a^{\\vert{}b}");
    CHECK(variableToLatex("a|") == "a^{}");
    const std::string tex = latexOf("K|star_a = 106\nK|star_b = 2*K|star_a\n");
    CHECK(tex.find("K^{star}_{a}") != std::string::npos);
    if (havePdflatex()) {
        std::string log;
        INFO(log);
        CHECK(pdflatexCompiles(tex, &log));
    }
}

TEST_CASE("CS-GAP-NAME-PIPE: the debug folder keeps the markdown tables intact", "[library-gaps][name-pipe]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_name_pipe_debug";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string code = "K|star_a = 106\nK|star_b = 2*K|star_a\n";
    { std::ofstream f(dir / "model.eescode"); f << code; }
    CoolSolveRunner runner((dir / "model.eescode").string());
    SolverOptions options;
    REQUIRE(runner.run(options));
    runner.generateDebugOutput((dir / "debug").string(), code);
    auto slurp = [&](const char* name) {
        std::ifstream in(dir / "debug" / name);
        REQUIRE(in.good());
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };
    // a raw | in a table cell would add a column: it is written \| (also inside the code span)
    const std::string variables = slurp("variables.md");
    CHECK(variables.find("`K\\|star_a`") != std::string::npos);
    CHECK(variables.find("`K|star_a`") == std::string::npos);
    const std::string residuals = slurp("coolsolve_residuals.md");
    CHECK(residuals.find("K\\|star_b") != std::string::npos);
    // every row of the variable table has the same number of cells as its header
    std::istringstream rows(variables);
    std::string row;
    int checked = 0;
    while (std::getline(rows, row)) {
        if (row.rfind("| ", 0) != 0) continue;
        std::string unescaped;
        for (size_t i = 0; i < row.size(); ++i) {
            if (row[i] == '\\' && i + 1 < row.size() && row[i + 1] == '|') { ++i; continue; }
            unescaped += row[i];
        }
        CHECK(std::count(unescaped.begin(), unescaped.end(), '|') == 8);
        ++checked;
    }
    CHECK(checked >= 3);
    CHECK(slurp("variables.csv").find("\"K|star_a\"") != std::string::npos);
    fs::remove_all(dir);
}

TEST_CASE("CS-GAP-NAME-PIPE: negative tests", "[library-gaps][name-pipe]") {
    // an operator is not part of a name: a pipe next to one is still a pipe of a name
    CHECK_FALSE(parseErrors("x = a +| b\n").empty());
    CHECK_FALSE(parseErrors("x = a || b\n").empty());
    // a comment is a comment, whatever it holds
    auto run = runModel("a|b = 1 \"a|b = 5\"\nc = a|b\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["c"], WithinAbs(1.0, 1e-12));
}

// ---------------------------------------------------------------------------
// CS-GAP-ACCENTED-IDENT: accented letters and ° in names
// ---------------------------------------------------------------------------

TEST_CASE("CS-GAP-ACCENTED-IDENT: register reproducer", "[library-gaps][accented-ident]") {
    // valid EES (the ISO 3046 file of CSL-0151): was "Parse failed: Line 1: Could not parse line"
    auto run = runModel("T_a_réf = 298\nx = T_a_réf/2\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["x"], WithinAbs(149.0, 1e-12));
    CHECK(run.warnings.empty());   // was a spurious "Unknown function 'p_atm_réf - a*'"
}

TEST_CASE("CS-GAP-ACCENTED-IDENT: the ISO 3046 pattern with names=value arguments and the reference values (CSL-0151)",
          "[library-gaps][accented-ident]") {
    auto run = runModel(
        "p_atm_réf = 100E3\nT_a_réf = 24.85\nphi_a_réf = 30\n"
        "p_w_s_réf = pressure(STEAM, T=T_a_réf, x=1)\n"
        "a = 0.7\nTest_réf = (phi_a_réf/100)*p_w_s_réf\n"
        "num = p_atm_réf - a*(phi_a_réf/100)*p_w_s_réf\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
    const double ps = CoolProp::PropsSI("P", "T", 24.85 + 273.15, "Q", 1, "Water");
    CHECK_THAT(run["p_w_s_réf"], WithinRel(ps, 1e-6));
    CHECK_THAT(run["num"], WithinRel(100e3 - 0.7 * 0.3 * run["p_w_s_réf"], 1e-12));
}

TEST_CASE("CS-GAP-ACCENTED-IDENT: the degree sign in names (CSL-0158)", "[library-gaps][accented-ident]") {
    // PROCEDURE Fittings(..., n°_90°elbow_reg, ...) and DELTAp_90°elbow_reg = n°_90°elbow_reg*K_90°elbow_reg*rho_i*(U^2/2)
    auto run = runModel(
        "PROCEDURE Fittings(f, U : n°_90°elbow_reg, K_90°elbow_reg)\n"
        "  n°_90°elbow_reg = f*2\n"
        "  K_90°elbow_reg = 0.3\n"
        "END\n"
        "rho_i = 1.2\nU = 10\n"
        "CALL Fittings(1, U : n°_90°elbow_reg, K_90°elbow_reg)\n"
        "DELTAp_90°elbow_reg=n°_90°elbow_reg*K_90°elbow_reg*rho_i*(U^2/2)\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["n°_90°elbow_reg"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["DELTAp_90°elbow_reg"], WithinRel(2.0 * 0.3 * 1.2 * 50.0, 1e-12));
}

TEST_CASE("CS-GAP-ACCENTED-IDENT: accented names in every position", "[library-gaps][accented-ident]") {
    SECTION("a name that starts with a non-ASCII letter") {
        auto run = runModel("é_x = 3\nÉnergie = é_x*2\nàb = É_X + 1\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["Énergie"], WithinAbs(6.0, 1e-12));
        CHECK_THAT(run["àb"], WithinAbs(4.0, 1e-12));
        CHECK(userVariables(run).size() == 3);   // é_x and É_X are one variable
    }
    SECTION("routines, locals, CALL outputs, arrays") {
        auto run = runModel(
            "FUNCTION fé(a_réf, b)\n  c_loc_é = a_réf*b\n  fé = c_loc_é + 1\nEND\n"
            "PROCEDURE pé(a_réf : r_éq)\n  r_éq = a_réf + 1\nEND\n"
            "x = fé(2, 3)\nCALL pé(4 : q_éq)\n"
            "DUPLICATE i = 1, 3\n  T_réf[i] = i*10\nEND\ns = T_réf[1] + T_réf[3]\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(7.0, 1e-12));
        CHECK_THAT(run["q_éq"], WithinAbs(5.0, 1e-12));
        CHECK_THAT(run["s"], WithinAbs(40.0, 1e-12));
        CHECK_THAT(run["T_réf[2]"], WithinAbs(20.0, 1e-12));
    }
    SECTION("string variable and $IF, accents in the value too") {
        auto run = runModel(
            "régime$ = 'Été'\n"
            "$IF régime$='ÉTÉ'\nx = 1\n$ELSE\nx = 2\n$ENDIF\n"
            "$IF RÉGIME$='hiver'\ny = 1\n$ELSE\ny = 2\n$ENDIF\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["y"], WithinAbs(2.0, 1e-12));
    }
    SECTION("the state of an INTEGRAL model and its CSV") {
        auto run = runIntegral(
            "T_réf = 1 + integral(dTdt_réf, t, 0, 1, 0.01)\n"
            "dTdt_réf = -T_réf\n"
            "$IntegralTable t:0.5 T_réf\n",
            [](SolverOptions& o) { o.integralMethod = "RK4"; });
        REQUIRE(run.solveOk);
        CHECK_THAT(run["T_réf"], WithinAbs(std::exp(-1.0), 1e-8));
        CHECK(run.table.toCSV().rfind("t,T_réf", 0) == 0);
        // the state is the same variable in any case (the equations, the integration variable)
        auto upper = runIntegral(
            "T_RÉF = 1 + integral(DTDT_RÉF, T, 0, 1, 0.01)\n"
            "dTdt_réf = -t_réf\n"
            "$IntegralTable t:0.5 T_réf\n",
            [](SolverOptions& o) { o.integralMethod = "RK4"; });
        REQUIRE(upper.solveOk);
        CHECK_THAT(upper["T_réf"], WithinAbs(std::exp(-1.0), 1e-8));
    }
    SECTION("the iterator of DUPLICATE is compared like any other name") {
        auto run = runModel("DUPLICATE ï = 1, 2\n  Y[Ï] = ï*3\nEND\ns = Y[1] + Y[2]\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["s"], WithinAbs(9.0, 1e-12));
    }
}

TEST_CASE("CS-GAP-ACCENTED-IDENT: names stay case-insensitive with accented letters", "[library-gaps][accented-ident]") {
    auto run = runModel("T_a_réf = 298\nx = T_A_RÉF + 1\ny = t_a_RÉf * 2\n");
    REQUIRE(run.solveOk);
    CHECK(userVariables(run).size() == 3);   // T_a_réf, x, y: the three spellings are one variable
    CHECK_THAT(run["x"], WithinAbs(299.0, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(596.0, 1e-12));
    // the comparator of every name map, and the helper shared with the string comparisons
    CaseInsensitiveLess less;
    CHECK_FALSE(less("T_a_réf", "T_A_RÉF"));
    CHECK_FALSE(less("T_A_RÉF", "T_a_réf"));
    CHECK(less("e", "é"));          // distinct letters stay distinct, and the order is a strict one
    CHECK_FALSE(less("é", "e"));
    CHECK(less("f", "é"));        // by UTF-8 bytes: after every ASCII letter
    CHECK(less("abc", "ABCD"));
    CHECK_FALSE(less("ABCD", "abc"));
    CHECK(less("a", "ä"));
    CHECK_FALSE(less("ä", "Ä"));
    CHECK_FALSE(less("Ä", "ä"));
    CHECK(eesStringEqual("n°_90°elbow", "N°_90°ELBOW"));
    CHECK_FALSE(eesStringEqual("n°", "n"));
    CHECK(eesFoldCase("T_A_RÉF°Ç") == "t_a_réf°ç");
    CHECK(eesFoldCase("ÀÖÞß×ÿ") == "àöþß×ÿ");   // × is not a letter, ß and ÿ have no capital in Latin-1
    // a std::map keyed with it finds any spelling
    std::map<std::string, int, CaseInsensitiveLess> m;
    m["Détarage"] = 1;
    CHECK(m.count("DÉTARAGE") == 1);
    CHECK(m.count("détarage") == 1);
    CHECK(m.count("Detarage") == 0);
}

TEST_CASE("CS-GAP-ACCENTED-IDENT: .initials, solution file and JSON in UTF-8", "[library-gaps][accented-ident]") {
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_accented_files";
    fs::remove_all(dir);
    fs::create_directories(dir);
    EESParser parser;
    auto parsed = parser.parse("T_a_réf = 298\nn°_90°elbow = T_a_réf/2\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    { std::ofstream f(dir / "m.initials"); f << "t_A_RÉF = 300 \"K\"\nn°_90°elbow = 150\n"; }   // any case
    CHECK(ir.loadInitialsFromFile((dir / "m.initials").string()) == 2);
    CHECK_THAT(*ir.getVariable("T_a_réf")->solutionValue, WithinAbs(300.0, 1e-12));
    CHECK(ir.getVariable("T_a_réf")->units == "K");
    // valid UTF-8 JSON
    const std::string json = ir.toJSON();
    CHECK(json.find("\"T_a_réf\"") != std::string::npos);
    CHECK(json.find("n°_90°elbow") != std::string::npos);
    fs::remove_all(dir);
}

TEST_CASE("CS-GAP-ACCENTED-IDENT: LaTeX output, accents and ° are typeset in math mode", "[library-gaps][accented-ident][latex]") {
    CHECK(variableToLatex("T_a_réf") == "T_{a,r\\text{é}f}");
    CHECK(variableToLatex("n°_90°elbow_reg") == "n{}^{\\circ}_{90{}^{\\circ}elbow,reg}");
    CHECK(variableToLatex("Énergie") == "\\text{É}nergie");
    CHECK(variableToLatex("µ_x") == "\\mu_{x}");
    const std::string tex = latexOf(
        "T_a_réf = 298\nx_réf = T_a_réf/2\nn°_90°elbow_reg = 2\nK_90°elbow = 0.5\nD = n°_90°elbow_reg*K_90°elbow\n");
    CHECK(tex.find("T_{a,r\\text{é}f}") != std::string::npos);
    CHECK(tex.find("n{}^{\\circ}") != std::string::npos);
    if (havePdflatex()) {
        std::string log;
        INFO(log);
        CHECK(pdflatexCompiles(tex, &log));
    }
}

TEST_CASE("CS-GAP-ACCENTED-IDENT: negative tests, characters that stay invalid give a clear error", "[library-gaps][accented-ident]") {
    SECTION("letters outside Latin-1 and symbols that are not part of a name") {
        CHECK(anyContains(parseErrors("Δp = 3\n"), "Invalid variable name 'Δp': the character 'Δ'"));
        CHECK(anyContains(parseErrors("x = a×b\n"), "Invalid variable name 'a×b': the character '×'"));
        CHECK(anyContains(parseErrors("x = a−b\n"), "Invalid variable name 'a−b': the character '−'"));   // U+2212 minus sign
        CHECK(anyContains(parseErrors("x = a\xC2\xA0+ b\n"), "Invalid variable name"));                  // a non-breaking space
        CHECK(anyContains(parseErrors("x_€ = 3\n"), "the character '€'"));
    }
    SECTION("a file saved in Latin-1 is not UTF-8: the byte is named") {
        const auto errors = parseErrors("T_a_r\xE9" "f = 298\n");   // é in Latin-1
        REQUIRE_FALSE(errors.empty());
        CHECK(anyContains(errors, "not valid UTF-8"));
        CHECK(anyContains(errors, "0xE9"));
    }
    SECTION("accented text in comments, strings and units is not a name") {
        auto run = runModel("x = 3 \"température [°C]\"\n{commentaire é}\ny$ = 'Été'\nz = 3 [°C]\n");
        REQUIRE(run.solveOk);
        CHECK(run.strings.at("y$") == "Été");
    }
    SECTION("a keyword followed by accents is a name") {
        auto run = runModel("ifé = 3\ncallé = ifé + 1\nendé = callé*2\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["endé"], WithinAbs(8.0, 1e-12));
    }
}

// ---------------------------------------------------------------------------
// CS-BUG-HASH-IDENT: `#` inside a name
// ---------------------------------------------------------------------------

TEST_CASE("CS-BUG-HASH-IDENT: register reproducer, the local Nusselt number of a procedure", "[library-gaps][hash-ident]") {
    // EES: yv = 6. Was a silent drop of the assignment inside the body (yv = 1), a parse error outside.
    auto run = runModel(
        "procedure th(x:y, z)\n"
        "Nu#_f = 2*x\n"
        "y = Nu#_f\n"
        "z = 1\n"
        "end\n"
        "x = 3\n"
        "call th(x: yv, zv)\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["yv"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["zv"], WithinAbs(1.0, 1e-12));
    // the same two lines in the main program
    auto main = runModel("Nu#_f = 6\ny = Nu#_f\n");
    REQUIRE(main.solveOk);
    CHECK_THAT(main["y"], WithinAbs(6.0, 1e-12));
}

TEST_CASE("CS-BUG-HASH-IDENT: the martin procedure pattern (CSL-0111, CSL-0161)", "[library-gaps][hash-ident]") {
    // h = Nu#_f*k/d_h with Nu#_f = Pr^(1/3)*(mu/mu_w)^(1/6)*0.122*(xi*Re^2*sin(2*beta))^0.374
    auto run = runModel(
        "PROCEDURE martin(Re_f, Pr_f, mu_f, mu_w, xi, beta, k_f, d_h : h_f)\n"
        "  Nu#_f = Pr_f^(1/3)*(mu_f/mu_w)^(1/6)*0.122*(xi*Re_f^2*sin(2*beta))^(0.374)\n"
        "  h_f = Nu#_f*k_f/d_h\n"
        "END\n"
        "beta = 60\nCALL martin(500, 3.5, 8e-4, 7e-4, 0.3, beta, 0.6, 0.004 : h)\n"
        "Nu#_f = h*0.004/0.6\n");   // and a variable of the main program with the same name
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    const double nu = std::pow(3.5, 1.0 / 3.0) * std::pow(8e-4 / 7e-4, 1.0 / 6.0) * 0.122 *
                      std::pow(0.3 * 500.0 * 500.0 * std::sin(2.0 * 60.0 * M_PI / 180.0), 0.374);
    CHECK_THAT(run["h"], WithinRel(nu * 0.6 / 0.004, 1e-12));
    CHECK_THAT(run["Nu#_f"], WithinRel(nu, 1e-12));   // the local of the procedure and the global do not mix
}

TEST_CASE("CS-BUG-HASH-IDENT: the built-in constants are unchanged", "[library-gaps][hash-ident]") {
    auto run = runModel("x = g#*2\ny = pi\ng#a = 3\nz = g#a + g#\nTableRun#_x = 5\nw = tablerun#_x\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(2.0 * 9.80665, 1e-12));
    CHECK_THAT(run["y"], WithinAbs(M_PI, 1e-15));
    CHECK_THAT(run["z"], WithinAbs(3.0 + 9.80665, 1e-12));
    CHECK_THAT(run["w"], WithinAbs(5.0, 1e-15));
    // g# is the constant, g#a is a variable
    CHECK(run.vars.count("g#a") == 1);
}

TEST_CASE("CS-BUG-HASH-IDENT: # in a name does not start a comment in .initials and .sol", "[library-gaps][hash-ident]") {
    // the readers cut a line at '#': the name was truncated and the value lost
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_hash_files";
    fs::remove_all(dir);
    fs::create_directories(dir);
    EESParser parser;
    auto parsed = parser.parse("Nu#_f = 6\ny = Nu#_f*2\nR#2 = y\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    {
        std::ofstream f(dir / "m.initials");
        f << "# initial values\n"
          << "Nu#_f = 4.5e+00 \"-\"\n"
          << "  y = 9 # comment after a space\n"
          << "R#2 = 3#glued comment\n"
          << "# R#2 = 100\n";
    }
    CHECK(ir.loadInitialsFromFile((dir / "m.initials").string()) == 3);
    CHECK_THAT(*ir.getVariable("Nu#_f")->solutionValue, WithinAbs(4.5, 1e-15));
    CHECK_THAT(*ir.getVariable("y")->solutionValue, WithinAbs(9.0, 1e-15));
    CHECK_THAT(*ir.getVariable("R#2")->solutionValue, WithinAbs(3.0, 1e-15));
    fs::remove_all(dir);
}

TEST_CASE("CS-BUG-HASH-IDENT: # in a name, in routines, arrays, $IF and INTEGRAL states", "[library-gaps][hash-ident]") {
    SECTION("function with # parameters and local, IF condition") {
        auto run = runModel(
            "FUNCTION f(Re#, n#)\n  c# = Re#*n#\n  IF (c# > 5) THEN f = c# ELSE f = -c#\nEND\n"
            "x = f(2, 3)\ny = f(1, 3)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(6.0, 1e-12));
        CHECK_THAT(run["y"], WithinAbs(-3.0, 1e-12));
    }
    SECTION("arrays and string variables") {
        auto run = runModel(
            "DUPLICATE i = 1, 2\n  Nu#[i] = i*4\nEND\ns = Nu#[1] + Nu#[2]\n"
            "fluid#$ = 'Water'\n$IF fluid#$='water'\nx = 1\n$ELSE\nx = 2\n$ENDIF\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["s"], WithinAbs(12.0, 1e-12));
        CHECK_THAT(run["x"], WithinAbs(1.0, 1e-12));
    }
    SECTION("INTEGRAL state") {
        auto run = runIntegral(
            "T#_f = 1 + integral(dT#_f, t, 0, 1, 0.01)\ndT#_f = -T#_f\n$IntegralTable t:0.5 T#_f\n",
            [](SolverOptions& o) { o.integralMethod = "RK4"; });
        REQUIRE(run.solveOk);
        CHECK_THAT(run["T#_f"], WithinAbs(std::exp(-1.0), 1e-8));
    }
    SECTION("names are case-insensitive") {
        auto run = runModel("NU#_F = 4\ny = nu#_f*2\n");
        REQUIRE(run.solveOk);
        CHECK(userVariables(run).size() == 2);
        CHECK_THAT(run["y"], WithinAbs(8.0, 1e-12));
    }
}

TEST_CASE("CS-BUG-HASH-IDENT: LaTeX output escapes the #", "[library-gaps][hash-ident][latex]") {
    CHECK(variableToLatex("Nu#_f") == "Nu\\#_{f}");
    CHECK(variableToLatex("Nu#") == "Nu\\#");
    const std::string tex = latexOf("Nu#_f = 6\ny = Nu#_f*2\n");
    CHECK(tex.find("Nu\\#_{f}") != std::string::npos);
    if (havePdflatex()) {
        std::string log;
        INFO(log);
        CHECK(pdflatexCompiles(tex, &log));
    }
}

TEST_CASE("CS-BUG-HASH-IDENT: negative tests", "[library-gaps][hash-ident]") {
    // a name starts with a letter
    CHECK(anyContains(parseErrors("#x = 3\n"), "Invalid variable name '#x'"));
    // a body that cannot be parsed is still loud (CS-BUG-BODY-LINE-DROP)
    CHECK(anyContains(parseErrors("FUNCTION f(a)\n  f = a@#\nEND\ny = f(1)\n"), "Invalid variable name"));
}

// ---------------------------------------------------------------------------
// The rest of the family: \ in names, keywords next to the new characters, the quote
// ---------------------------------------------------------------------------

TEST_CASE("CS-GAP-NAME-SYMBOL: a backslash inside a name (EES namespace separator) works and is typeset", "[library-gaps][name-symbol][latex]") {
    // dU\dtau and time_h\day: names of the EES files of the collection
    auto run = runModel("dU\\dtau = 3\ntime_h\\day = 2\nq = dU\\dtau*time_h\\day\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["q"], WithinAbs(6.0, 1e-12));
    CHECK(variableToLatex("dU\\dtau") == "dU\\backslash{}dtau");
    CHECK(variableToLatex("time_h\\day") == "time_{h\\backslash{}day}");
    const std::string tex = latexOf("dU\\dtau = 3\ntime_h\\day = 2\nq = dU\\dtau*time_h\\day\n");
    if (havePdflatex()) {
        std::string log;
        INFO(log);
        CHECK(pdflatexCompiles(tex, &log));
    }
    // a backslash at the end of a name, or before a digit, is not a namespace separator
    CHECK_FALSE(parseErrors("a\\ = 3\n").empty());
}

TEST_CASE("CS-GAP-NAME-SYMBOL: a keyword followed by one of the new characters is a name", "[library-gaps][name-symbol]") {
    // CALL%, IF|, FUNCTION# ... start with a keyword but are variables
    auto run = runModel(
        "call% = 1\nif|x = 2\nfunction#a = 3\nduplicate% = 4\nrepeat|n = 5\nmodule° = 6\nsubprogram#1 = 7\n"
        "s = call%+if|x+function#a+duplicate%+repeat|n+module°+subprogram#1\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["s"], WithinAbs(28.0, 1e-12));
    // inside a routine too
    auto body = runModel(
        "FUNCTION f(a)\n  if% = a*2\n  end|x = if% + 1\n  f = end|x\nEND\ny = f(3)\n");
    REQUIRE(body.solveOk);
    CHECK_THAT(body["y"], WithinAbs(7.0, 1e-12));
    // the keywords themselves are unchanged
    auto kw = runModel("FUNCTION f(a)\n  IF (a > 1) THEN f = 1 ELSE f = 2\nEND\ny = f(3)\n");
    REQUIRE(kw.solveOk);
    CHECK_THAT(kw["y"], WithinAbs(1.0, 1e-12));
}

TEST_CASE("CS-GAP-NAME-SYMBOL: LaTeX of a whole report with every kind of name compiles", "[library-gaps][name-symbol][latex]") {
    if (!havePdflatex()) SKIP("pdflatex is not installed");
    const std::string tex = latexOf(
        "C% = 0.5\nH% = 0.06\nO% = 1-C%-H%\nm_boil = C%/12\n"
        "K|star_a = 106\nK|star_b = 2*K|star_a\nG|o_m = 3\n"
        "T_a_réf = 298\nx_réf = T_a_réf/2\nn°_90°elbow_reg = 2\nK_90°elbow = 0.5\nD_90°elbow = n°_90°elbow_reg*K_90°elbow\n"
        "Nu#_f = 6\ny = Nu#_f*2\n"
        "dU\\dtau = 3\nq = dU\\dtau*2\n"
        "fluid%$ = 'R134a'\nh%q = enthalpy(fluid%$, T=300, x=0)\nÉnergie|2 = 3\nµ_x = 2\n"
        "X%[1] = 1\nX%[2] = X%[1]*2\n");
    std::string log;
    INFO(log);
    CHECK(pdflatexCompiles(tex, &log));
}

// ============================================================================
// CS-GAP-FLUIDS-MIXTURE: CoolProp mixture fluids `A[x]&B[y]` in property calls
// ============================================================================

namespace {

constexpr const char* kMix = "R245fa[0.5]&R134a[0.5]";

/// Number as model text, with all its digits.
std::string numText(double v) {
    std::ostringstream o;
    o << std::setprecision(17) << v;
    return o.str();
}

/// Value and derivatives of a property call `call` (written with the variables Tv, Pv) at (T [C], P [Pa]).
struct MixtureCall { double value, dT, dP; };

MixtureCall evalMixtureCall(const std::string& call, double T_C, double P_Pa, const CoolPropConfig& cfg) {
    EESParser parser;
    auto parsed = parser.parse("y = " + call + "\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    ExpressionEvaluator ev(2, cfg);
    ev.setVariable("Tv", ADValue::independent(T_C, 0, 2));
    ev.setVariable("Pv", ADValue::independent(P_Pa, 1, 2));
    ADValue r = ev.evaluate(ir.getEquations()[0].rhs);
    return {r.value, r.gradient[0], r.gradient[1]};
}

/// Warnings of a run that talk about a mixture.
size_t mixtureWarnings(const ModelRun& run) {
    size_t n = 0;
    for (const auto& w : run.warnings)
        if (w.find("mixture") != std::string::npos) ++n;
    return n;
}

}  // namespace

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: register reproducer, enthalpy of a R245fa/R134a mixture", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
fl$ = 'R245fa[0.5]&R134a[0.5]'
h = enthalpy(fl$, T=50, P=5E5)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["h"], WithinRel(439716.94640022865, 1e-9));   // CoolProp 7.2, as PropsSI
    CHECK_THAT(run["h"], WithinRel(CoolProp::PropsSI("H", "T", 323.15, "P", 5e5, kMix), 1e-9));
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: writing the fluid as a literal, with HEOS:: prefix, spaces, EES name case", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
h1 = enthalpy('R245fa[0.5]&R134a[0.5]', T=50, P=5E5)
h2 = enthalpy('HEOS::R245fa[0.5]&R134a[0.5]', T=50, P=5E5)
h3 = enthalpy(' r245fa[ 0.5 ] & R134a[0.5] ', T=50, P=5E5)
fl$ = 'heos::R245fa[0.5]&R134a[0.5]'
h4 = enthalpy(fl$, T=50, P=5E5)
h5 = enthalpy('R245fa[5E-1]&R134a[.5]', T=50, P=5E5)
hw1 = enthalpy('Water[0.5]&Ethanol[0.5]', T=50, P=1E5)
hw2 = enthalpy('Steam[0.5]&Ethanol[0.5]', T=50, P=1E5)
hw3 = enthalpy('R718[0.5]&Ethanol[0.5]', T=50, P=1E5)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (const char* v : {"h2", "h3", "h4", "h5"})
        CHECK_THAT(run[v], WithinRel(run["h1"], 1e-12));
    CHECK_THAT(run["h1"], WithinRel(439716.94640022865, 1e-9));
    CHECK_THAT(run["hw2"], WithinRel(run["hw1"], 1e-12));   // EES aliases of the registry name the same component
    CHECK_THAT(run["hw3"], WithinRel(run["hw1"], 1e-12));
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: the pure end points equal the pure fluids", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
mixA$ = 'R245fa[1]&R134a[0]'
mixB$ = 'R245fa[0]&R134a[1]'
hA = enthalpy(mixA$, T=50, P=5E5)
hB = enthalpy(mixB$, T=50, P=5E5)
hA_pure = enthalpy('R245fa', T=50, P=5E5)
hB_pure = enthalpy('R134a', T=50, P=5E5)
TA = temperature(mixA$, P=5E5, x=0)
TB = temperature(mixB$, P=5E5, x=0)
TA_pure = temperature('R245fa', P=5E5, x=0)
TB_pure = temperature('R134a', P=5E5, x=0)
TdA = temperature(mixA$, P=5E5, x=1)
rhoA = density(mixA$, T=20, P=5E5)
rhoA_pure = density('R245fa', T=20, P=5E5)
h_nopair = enthalpy('R245fa[1]&Water[0]', T=50, P=5E5)
Tsat_pure = T_sat('R245fa[1]&R134a[0]', P=5E5)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["hA"], WithinRel(266776.2, 1e-6));        // values given by CoolProp 7.2
    CHECK_THAT(run["hB"], WithinRel(440111.7, 1e-6));
    CHECK_THAT(run["TA"] + 273.15, WithinAbs(335.913, 1e-3));
    CHECK_THAT(run["TB"] + 273.15, WithinAbs(288.885, 1e-3));
    CHECK_THAT(run["hA"], WithinRel(run["hA_pure"], 1e-6));
    CHECK_THAT(run["hB"], WithinRel(run["hB_pure"], 1e-6));
    CHECK_THAT(run["TA"], WithinAbs(run["TA_pure"], 1e-4));
    CHECK_THAT(run["TB"], WithinAbs(run["TB_pure"], 1e-4));
    CHECK_THAT(run["TdA"], WithinAbs(run["TA"], 1e-4));      // no glide for a pure composition
    CHECK_THAT(run["rhoA"], WithinRel(run["rhoA_pure"], 1e-6));
    // a zero fraction drops the component, so a pair without binary parameters is no problem either,
    // and the result is the pure fluid in every respect: its functions (T_sat), and no mixture warning
    CHECK_THAT(run["h_nopair"], WithinRel(run["hA_pure"], 1e-12));
    CHECK_THAT(run["Tsat_pure"], WithinAbs(run["TA_pure"], 1e-9));
    CHECK(mixtureWarnings(run) == 0);
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: bubble and dew temperatures at 5 bar show the glide", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
fl$ = 'R245fa[0.5]&R134a[0.5]'
Tb = temperature(fl$, P=5E5, x=0)
Td = temperature(fl$, P=5E5, x=1)
Tmid = temperature(fl$, P=5E5, x=0.5)
rho_b = density(fl$, P=5E5, x=0)
rho_d = density(fl$, P=5E5, x=1)
rho_pureA = density('R245fa', T=Tb, x=0)
rho_pureB = density('R134a', T=Tb, x=0)
Pb = pressure(fl$, T=20, x=0)
Pd = pressure(fl$, T=20, x=1)
q = quality(fl$, T=40, P=5E5)
h_b = enthalpy(fl$, P=5E5, x=0)
h_d = enthalpy(fl$, P=5E5, x=1)
h_m = enthalpy(fl$, P=5E5, x=0.5)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["Tb"] + 273.15, WithinRel(306.0, 1e-3));      // values given by CoolProp 7.2
    CHECK_THAT(run["Td"] + 273.15, WithinRel(319.4, 1e-3));
    CHECK_THAT(run["Tb"] + 273.15, WithinRel(306.0143512292315, 1e-9));
    CHECK_THAT(run["Td"] + 273.15, WithinRel(319.39495958183375, 1e-9));
    CHECK(run["Td"] - run["Tb"] > 10.0);                          // glide of 13 K
    CHECK(run["Tmid"] > run["Tb"]);
    CHECK(run["Tmid"] < run["Td"]);
    // liquid density between the pure-component values at the same temperature
    CHECK(run["rho_b"] > std::min(run["rho_pureA"], run["rho_pureB"]));
    CHECK(run["rho_b"] < std::max(run["rho_pureA"], run["rho_pureB"]));
    CHECK(run["rho_d"] < run["rho_b"] / 10.0);
    CHECK(run["Pb"] > run["Pd"]);                                 // bubble pressure above dew pressure
    CHECK_THAT(run["Pb"], WithinRel(CoolProp::PropsSI("P", "T", 293.15, "Q", 0, kMix), 1e-9));
    CHECK_THAT(run["Pd"], WithinRel(CoolProp::PropsSI("P", "T", 293.15, "Q", 1, kMix), 1e-9));
    CHECK_THAT(run["q"], WithinRel(CoolProp::PropsSI("Q", "T", 313.15, "P", 5e5, kMix), 1e-9));   // two-phase in the glide
    CHECK(run["q"] > 0.0);
    CHECK(run["q"] < 1.0);
    CHECK(run["h_b"] < run["h_m"]);
    CHECK(run["h_m"] < run["h_d"]);
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: the property functions of a mixture", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
fl$ = 'R245fa[0.5]&R134a[0.5]'
h = enthalpy(fl$, T=50, P=5E5)
s = entropy(fl$, T=50, P=5E5)
rho = density(fl$, T=50, P=5E5)
v = volume(fl$, T=50, P=5E5)
c_p = cp(fl$, T=50, P=5E5)
c_v = cv(fl$, T=50, P=5E5)
u = intenergy(fl$, T=50, P=5E5)
mu = viscosity(fl$, T=50, P=5E5)
lambda = conductivity(fl$, T=50, P=5E5)
Pr = prandtl(fl$, T=50, P=5E5)
a = soundspeed(fl$, T=50, P=5E5)
Z = compressibilityfactor(fl$, T=50, P=5E5)
P_back = pressure(fl$, T=50, rho=rho)
rho_l = density(fl$, T=20, P=5E5)
mu_l = viscosity(fl$, T=20, P=5E5)
lambda_l = conductivity(fl$, T=20, P=5E5)
mu_sat = viscosity(fl$, P=5E5, x=1)
nu = kinematicviscosity(fl$, T=50, P=5E5)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    const double T = 323.15, P = 5e5;
    auto ps = [&](const char* out, double T_K = 323.15) { return CoolProp::PropsSI(out, "T", T_K, "P", P, kMix); };
    CHECK_THAT(run["h"], WithinRel(ps("H"), 1e-9));
    CHECK_THAT(run["s"], WithinRel(ps("S"), 1e-9));
    CHECK_THAT(run["rho"], WithinRel(ps("D"), 1e-9));
    CHECK_THAT(run["v"], WithinRel(1.0 / ps("D"), 1e-9));
    CHECK_THAT(run["c_p"], WithinRel(ps("C"), 1e-9));
    CHECK_THAT(run["c_v"], WithinRel(ps("CVMASS"), 1e-9));
    CHECK_THAT(run["u"], WithinRel(ps("U"), 1e-9));
    CHECK_THAT(run["mu"], WithinRel(ps("V"), 1e-9));
    CHECK_THAT(run["lambda"], WithinRel(ps("L"), 1e-9));
    CHECK_THAT(run["Pr"], WithinRel(ps("Prandtl"), 1e-9));
    CHECK_THAT(run["a"], WithinRel(ps("A"), 1e-9));
    CHECK_THAT(run["Z"], WithinRel(ps("Z"), 1e-9));
    CHECK_THAT(run["P_back"], WithinRel(P, 1e-6));
    CHECK_THAT(run["rho_l"], WithinRel(ps("D", 293.15), 1e-9));
    CHECK_THAT(run["mu_l"], WithinRel(ps("V", 293.15), 1e-9));
    CHECK_THAT(run["lambda_l"], WithinRel(ps("L", 293.15), 1e-9));
    CHECK_THAT(run["mu_sat"], WithinRel(CoolProp::PropsSI("V", "P", P, "Q", 1, kMix), 1e-9));
    CHECK_THAT(run["nu"], WithinRel(ps("V") / ps("D"), 1e-9));
    (void)T;
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: a mixture of three components", "[library-gaps][fluids-mixture]") {
    const std::string ng = "Methane[0.7]&Ethane[0.2]&Propane[0.1]";
    auto run = runModel("fl$ = '" + ng + "'\nrho = density(fl$, T=27, P=1E6)\nh = enthalpy(fl$, T=27, P=1E6)\nM = molarmass(fl$)\n");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["rho"], WithinRel(CoolProp::PropsSI("D", "T", 300.15, "P", 1e6, ng), 1e-9));
    CHECK_THAT(run["h"], WithinRel(CoolProp::PropsSI("H", "T", 300.15, "P", 1e6, ng), 1e-9));
    CHECK_THAT(run["M"], WithinRel(0.7 * molarMassKgPerKmol("Methane") + 0.2 * molarMassKgPerKmol("Ethane") +
                                   0.1 * molarMassKgPerKmol("Propane"), 1e-9));
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: the molar mass is mole-weighted, in kg/kmol", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
fl$ = 'R245fa[0.3]&R134a[0.7]'
M = molarmass(fl$)
M_lit = molarmass('R245fa[0.3]&R134a[0.7]')
M_state = molarmass(fl$, T=50, P=5E5)
M_A = molarmass('R245fa')
M_B = molarmass('R134a')
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    const double expected = 0.3 * molarMassKgPerKmol("R245fa") + 0.7 * molarMassKgPerKmol("R134a");
    CHECK_THAT(run["M"], WithinRel(expected, 1e-9));
    CHECK_THAT(run["M_lit"], WithinRel(expected, 1e-9));
    CHECK_THAT(run["M_state"], WithinRel(expected, 1e-9));
    CHECK_THAT(run["M"], WithinRel(0.3 * run["M_A"] + 0.7 * run["M_B"], 1e-9));
    CHECK(run["M"] > 100.0);   // kg/kmol, not kg/mol
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: a mixture in FUNCTION and PROCEDURE bodies", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
FUNCTION hmix(fl$, T, P)
  hmix = enthalpy(fl$, T=T, P=P)
END
FUNCTION hlit(T, P)
  hlit = enthalpy('R245fa[0.5]&R134a[0.5]', T=T, P=P)
END
PROCEDURE props(fl$, T, P : rho, mu)
  rho = density(fl$, T=T, P=P)
  mu = viscosity('R245fa[0.5]&R134a[0.5]', T=T, P=P)
END
fl$ = 'R245fa[0.5]&R134a[0.5]'
h1 = hmix(fl$, 50, 5E5)
h2 = hlit(50, 5E5)
h3 = hmix('R245fa[0.5]&R134a[0.5]', 50, 5E5)
CALL props(fl$, 20, 5E5 : rho, mu)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["h1"], WithinRel(439716.94640022865, 1e-9));
    CHECK_THAT(run["h2"], WithinRel(run["h1"], 1e-12));
    CHECK_THAT(run["h3"], WithinRel(run["h1"], 1e-12));
    CHECK_THAT(run["rho"], WithinRel(CoolProp::PropsSI("D", "T", 293.15, "P", 5e5, kMix), 1e-9));
    CHECK_THAT(run["mu"], WithinRel(CoolProp::PropsSI("V", "T", 293.15, "P", 5e5, kMix), 1e-9));
    CHECK(mixtureWarnings(run) == 1);   // the same mixture in the equations, a FUNCTION and a PROCEDURE: one warning
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: a mixture used only inside a FUNCTION body still warns", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
FUNCTION hlit(T, P)
  hlit = enthalpy('R245fa[0.5]&R134a[0.5]', T=T, P=P)
END
h = hlit(50, 5E5)
)");
    REQUIRE(run.solveOk);
    CHECK(mixtureWarnings(run) == 1);
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: the warning is emitted once per mixture and run", "[library-gaps][fluids-mixture]") {
    SECTION("one mixture used many times: one warning, saying that the properties are predictive and unreliable") {
        auto run = runModel(R"(
fl$ = 'R245fa[0.5]&R134a[0.5]'
h1 = enthalpy(fl$, T=50, P=5E5)
h2 = enthalpy(fl$, T=40, P=5E5)
s1 = entropy('R245fa[0.5]&R134a[0.5]', T=50, P=5E5)
T1 = temperature(fl$, P=5E5, x=0)
rho = density(fl$, T=T1, x=1)
)");
        REQUIRE(run.solveOk);
        REQUIRE(mixtureWarnings(run) == 1);
        for (const auto& w : run.warnings) {
            if (w.find("mixture") == std::string::npos) continue;
            CHECK(w.find("R245fa[0.5]&R134a[0.5]") != std::string::npos);
            CHECK(w.find("CoolProp") != std::string::npos);
            CHECK(w.find("unreliable") != std::string::npos);
            CHECK(w.find("two-phase") != std::string::npos);
            CHECK(w.find("transport") != std::string::npos);
        }
    }
    SECTION("two different mixtures: two warnings; a pure fluid: none") {
        auto run = runModel(R"(
h1 = enthalpy('R245fa[0.5]&R134a[0.5]', T=50, P=5E5)
h2 = enthalpy('R245fa[0.2]&R134a[0.8]', T=50, P=5E5)
h3 = enthalpy('R245fa', T=50, P=5E5)
)");
        REQUIRE(run.solveOk);
        CHECK(mixtureWarnings(run) == 2);
    }
    SECTION("pure fluids only: no warning") {
        auto run = runModel("h = enthalpy('R245fa', T=50, P=5E5)\nrho = density('Water', T=50, P=1E5)\n");
        REQUIRE(run.solveOk);
        CHECK(mixtureWarnings(run) == 0);
    }
    SECTION("the warning comes again in the next run") {
        const std::string code = "h = enthalpy('R245fa[0.5]&R134a[0.5]', T=50, P=5E5)\n";
        CHECK(mixtureWarnings(runModel(code)) == 1);
        CHECK(mixtureWarnings(runModel(code)) == 1);
    }
    SECTION("no unit or fluid-name hint fires on the mixture string") {
        auto run = runModel("fl$ = 'R245fa[0.5]&R134a[0.5]'\nh = enthalpy(fl$, T=50, P=5E5)\n"
                            "h2 = enthalpy('R134a[0.5]&R245fa[0.5]', T=50, P=5E5)\n"
                            "h3 = enthalpy('Water[0.5]&Ethanol[0.5]', T=50, P=1E5)\n");
        REQUIRE(run.solveOk);
        CHECK(run.warnings.size() == mixtureWarnings(run));   // the mixture warnings are the only ones
        CHECK(run.warnings.size() == 3);
    }
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: Newton solves an unknown temperature through a mixture property", "[library-gaps][fluids-mixture]") {
    SECTION("single-phase vapor") {
        auto run = runModel(R"(
fl$ = 'R245fa[0.5]&R134a[0.5]'
h1 = 439716.94640022865
h1 = enthalpy(fl$, T=T1, P=5E5)
)");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["T1"], WithinAbs(50.0, 1e-6));
    }
    SECTION("inside the glide (two-phase), and subcooled liquid") {
        const double hGlide = CoolProp::PropsSI("H", "T", 313.15, "P", 5e5, kMix);
        const double hLiquid = CoolProp::PropsSI("H", "T", 293.15, "P", 5e5, kMix);
        auto run = runModel(
            "fl$ = 'R245fa[0.5]&R134a[0.5]'\n"
            "hg = " + numText(hGlide) + "\nhg = enthalpy(fl$, T=Tg, P=5E5)\n"
            "hl = " + numText(hLiquid) + "\nhl = enthalpy(fl$, T=Tl, P=5E5)\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["Tg"], WithinAbs(40.0, 1e-5));
        CHECK_THAT(run["Tl"], WithinAbs(20.0, 1e-5));
    }
    SECTION("two unknowns: pressure and temperature of a state given by density and enthalpy-type equations") {
        const double rho = CoolProp::PropsSI("D", "T", 323.15, "P", 5e5, kMix);
        const double h = CoolProp::PropsSI("H", "T", 323.15, "P", 5e5, kMix);
        auto run = runModel(
            "fl$ = 'R245fa[0.5]&R134a[0.5]'\n"
            "rho0 = " + numText(rho) + "\nh0 = " + numText(h) + "\n"
            "rho0 = density(fl$, T=Tx, P=Px)\nh0 = enthalpy(fl$, T=Tx, P=Px)\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["Tx"], WithinAbs(50.0, 1e-4));
        CHECK_THAT(run["Px"], WithinRel(5e5, 1e-6));
    }
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: derivatives of mixture properties are analytic and agree with finite differences", "[library-gaps][fluids-mixture]") {
    CoolPropConfig analytic;                       // AbstractState + CoolProp's first_partial_deriv()
    analytic.useAbstractState = true;
    analytic.enableAnalyticalDerivatives = true;
    auto central = [&](const std::string& call, double T, double P, double& dT, double& dP) {
        const double hT = 1e-3, hP = 1.0;           // K, Pa
        dT = (evalMixtureCall(call, T + hT, P, analytic).value - evalMixtureCall(call, T - hT, P, analytic).value) / (2 * hT);
        dP = (evalMixtureCall(call, T, P + hP, analytic).value - evalMixtureCall(call, T, P - hP, analytic).value) / (2 * hP);
    };
    struct State { const char* label; double T_C, P; };
    // vapor, subcooled liquid, and inside the glide (two-phase)
    for (const State& st : {State{"vapor", 50.0, 5e5}, State{"liquid", 20.0, 5e5}, State{"glide", 40.0, 5e5}}) {
        for (const char* fn : {"enthalpy", "density", "entropy", "cp"}) {
            const std::string call = std::string(fn) + "('" + kMix + "', T=Tv, P=Pv)";
            INFO(st.label << " " << call);
            MixtureCall a = evalMixtureCall(call, st.T_C, st.P, analytic);
            double dT, dP;
            central(call, st.T_C, st.P, dT, dP);
            CHECK_THAT(a.dT, WithinRel(dT, 2e-3));
            CHECK_THAT(a.dP, WithinRel(dP, 2e-3));
            CHECK(std::abs(a.dT) > 0.0);
        }
    }
    // the PropsSI path (no AbstractState) gives the same values
    CoolPropConfig propsSI;
    propsSI.useAbstractState = false;
    propsSI.enableAnalyticalDerivatives = false;
    MixtureCall viaAS = evalMixtureCall(std::string("enthalpy('") + kMix + "', T=Tv, P=Pv)", 50.0, 5e5, analytic);
    MixtureCall viaPS = evalMixtureCall(std::string("enthalpy('") + kMix + "', T=Tv, P=Pv)", 50.0, 5e5, propsSI);
    CHECK_THAT(viaPS.value, WithinRel(viaAS.value, 1e-9));
    CHECK_THAT(viaPS.dT, WithinRel(viaAS.dT, 1e-3));
    CHECK_THAT(viaPS.dP, WithinRel(viaAS.dP, 1e-3));
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: PHASE$ of a mixture", "[library-gaps][fluids-mixture]") {
    auto phase = [](double T_C) {
        EESParser parser;
        auto parsed = parser.parse("ph$ = PHASE$('R245fa[0.5]&R134a[0.5]', T=" + numText(T_C) + ", P=5E5)\n");
        REQUIRE(parsed.success);
        IR ir = IR::fromAST(parsed.program);
        ExpressionEvaluator ev(1);
        return ev.evaluateString(ir.getEquations()[0].rhs);
    };
    CHECK(phase(50.0) == "gas");
    CHECK(phase(20.0) == "liquid");
    CHECK(phase(40.0) == "twophase");   // inside the glide
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: pure fluids and other fluid types are unaffected", "[library-gaps][fluids-mixture]") {
    auto run = runModel(R"(
h = enthalpy('R245fa', T=50, P=5E5)
rho = density(Water, T=20, P=1E5)
hN2 = enthalpy('N2', T=25)
rho_air = density('AirH2O', T=25, P=101325, R=0.5)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["h"], WithinRel(266776.2, 1e-6));
    CHECK_THAT(run["rho"], WithinRel(CoolProp::PropsSI("D", "T", 293.15, "P", 1e5, "Water"), 1e-9));
    CHECK(mixtureWarnings(run) == 0);
    auto unknown = runModel("h = enthalpy('NoSuchFluid', T=50, P=5E5)\n");
    CHECK_FALSE(unknown.solveOk);
    CHECK(anyContains({unknown.message}, "Unknown fluid"));
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: negative tests, the mixture string is checked with a clear error", "[library-gaps][fluids-mixture]") {
    auto fails = [](const std::string& fluid, const std::string& needle, const std::string& call = "enthalpy") {
        auto run = runModel("h = " + call + "('" + fluid + "', T=50, P=5E5)\n");
        INFO(fluid << " -> " << run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, needle));
        CHECK_FALSE(anyContains({run.message}, "Unknown fluid"));
    };
    SECTION("fractions that do not sum to 1") {
        fails("R245fa[0.5]&R134a[0.3]", "sum");
        fails("R245fa[0.5]&R134a[0.3]", "0.8");
        fails("R245fa[0.6]&R134a[0.6]", "sum");
        fails("R245fa[0.5]&R134a[0.5000015]", "sum");   // beyond the 1e-6 tolerance
    }
    SECTION("a fraction that is negative, larger than 1 or not a number") {
        fails("R245fa[1.5]&R134a[-0.5]", "mole fraction");
        fails("R245fa[-0.2]&R134a[1.2]", "mole fraction");
        fails("R245fa[x]&R134a[0.5]", "not a number");
        fails("R245fa[]&R134a[0.5]", "not a number");
    }
    SECTION("an unknown component names the component") {
        fails("R245fa[0.5]&Foo[0.5]", "'Foo'");
        fails("Bar[0.5]&R134a[0.5]", "'Bar'");
        fails("R245fa[0.5]&R134a[0.3]&Baz[0.2]", "'Baz'");
    }
    SECTION("a pair without binary parameters is reported by CoolProp") {
        fails("R245fa[0.5]&Water[0.5]", "binary pair");
        fails("R245fa[0.5]&R245fa[0.5]", "twice");
    }
    SECTION("malformed strings") {
        fails("R245fa[0.5]&", "empty");
        fails("&R134a[0.5]", "empty");
        fails("R245fa[0.5]&&R134a[0.5]", "empty");
        fails("R245fa&R134a", "mole fraction");              // fractions are required
        fails("R245fa[0.5]&R134a", "mole fraction");
        fails("R245fa[0.5&R134a[0.5]", "mole fraction");     // unbalanced brackets
        fails("R245fa[0.5]]&R134a[0.5]", "mole fraction");
        fails("[0.5]&R134a[0.5]", "component");
        fails("REFPROP::R245fa[0.5]&R134a[0.5]", "HEOS");   // other backends are not accepted
    }
    SECTION("the same errors reach every property function and a string variable") {
        fails("R245fa[0.5]&R134a[0.3]", "sum", "density");
        fails("R245fa[0.5]&Foo[0.5]", "'Foo'", "molarmass");
        auto run = runModel("fl$ = 'R245fa[0.5]&Foo[0.5]'\nh = enthalpy(fl$, T=50, P=5E5)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "'Foo'"));
    }
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: what CoolProp cannot do for a mixture gives a clear error", "[library-gaps][fluids-mixture]") {
    SECTION("T_sat and P_sat: a mixture has bubble and dew points, not one saturation point") {
        auto run = runModel("T = T_sat('R245fa[0.5]&R134a[0.5]', P=5E5)\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "mixture"));
        CHECK(anyContains({run.message}, "x=0"));
        auto run2 = runModel("P = P_sat('R245fa[0.5]&R134a[0.5]', T=20)\n");
        CHECK_FALSE(run2.solveOk);
        CHECK(anyContains({run2.message}, "mixture"));
    }
    SECTION("fluid constants such as T_crit") {
        auto run = runModel("T = T_crit('R245fa[0.5]&R134a[0.5]')\n");
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "'R245fa[0.5]&R134a[0.5]'"));
        CHECK(anyContains({run.message}, "not available"));
    }
    SECTION("surface tension is not implemented for mixtures: CoolProp's message is reported") {
        auto run = runModel("sigma = surfacetension('R245fa[0.5]&R134a[0.5]', T=20, x=0)\n");
        CHECK_FALSE(run.verified);
        bool reported = false;
        for (const auto& w : run.warnings)
            if (w.find("CoolProp error in surfacetension()") != std::string::npos &&
                w.find("not implemented for mixtures") != std::string::npos) reported = true;
        CHECK(reported);
    }
}

TEST_CASE("CS-GAP-FLUIDS-MIXTURE: LaTeX report of a model with a mixture compiles", "[library-gaps][fluids-mixture][latex]") {
    if (!havePdflatex()) SKIP("pdflatex is not installed");
    const std::string tex = latexOf("fl$ = 'R245fa[0.5]&R134a[0.5]'\nh = enthalpy(fl$, T=50, P=5E5)\n"
                                    "rho = density('R245fa[0.5]&R134a[0.5]', T=T1, P=5E5)\nT1 = 30\n");
    std::string log;
    INFO(log);
    CHECK(pdflatexCompiles(tex, &log));
}

// ============================================================================
// CS-BUG-COMMON-PROC: $COMMON is not supported (maintainer decision): the body of a
// FUNCTION/PROCEDURE saw 0 for the variables of its $COMMON list while the run reported
// SUCCESS. It is now an error that names the directive, the variables and the routine and
// says how to rewrite the model: pass the shared variables as input arguments.
// ============================================================================

namespace {

/// The $COMMON error of `code`, or an empty string.
std::string commonError(const std::string& code) {
    for (const auto& m : parseErrors(code))
        if (m.find("$COMMON") != std::string::npos) return m;
    return "";
}

}  // namespace

TEST_CASE("CS-BUG-COMMON-PROC: the register reproducer is rejected, never SUCCESS", "[library-gaps][common-proc]") {
    const std::string code =
        "PROCEDURE pq(a:y)\n"
        "$common b\n"
        "y = a + b\n"
        "END\n"
        "b = 21\n"
        "a = 1\n"
        "CALL pq(a:y)\n"
        "z = y*2\n"
        "w = z + b\n";
    auto run = runModel(code);
    INFO(run.message);
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    // the message names the directive, the variable and the routine, and says how to rewrite
    CHECK(anyContains({run.message}, "$COMMON is not supported"));
    CHECK(anyContains({run.message}, "pass b as an input argument of PROCEDURE pq"));
    CHECK(anyContains({run.message}, "in each CALL"));
}

TEST_CASE("CS-BUG-COMMON-PROC: the error points at the $COMMON line", "[library-gaps][common-proc]") {
    EESParser parser;
    auto result = parser.parse(
        "x = 1\n"
        "PROCEDURE pq(a:y)\n"
        "$common b\n"
        "y = a + b\n"
        "END\n"
        "CALL pq(x:y)\n");
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].line == 3);
    CHECK(result.errors[0].context.find("$common b") != std::string::npos);
}

TEST_CASE("CS-BUG-COMMON-PROC: a FUNCTION with $COMMON is rejected too", "[library-gaps][common-proc]") {
    auto run = runModel(
        "FUNCTION fq(a)\n"
        "$common b\n"
        "fq = a + b\n"
        "END\n"
        "b = 21\n"
        "a = 1\n"
        "y = fq(a)\n");
    INFO(run.message);
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(anyContains({run.message}, "$COMMON is not supported"));
    CHECK(anyContains({run.message}, "pass b as an input argument of FUNCTION fq"));
    CHECK(anyContains({run.message}, "in each call"));
}

TEST_CASE("CS-BUG-COMMON-PROC: the message lists every variable of the directive", "[library-gaps][common-proc]") {
    SECTION("upper case $COMMON, comma separated") {
        const std::string msg = commonError(
            "PROCEDURE pq(a:y)\n$COMMON b, c, d\ny = a + b + c + d\nEND\n");
        INFO(msg);
        CHECK(msg.find("pass b, c, d as input arguments of PROCEDURE pq") != std::string::npos);
    }
    SECTION("semicolon separators") {
        const std::string msg = commonError(
            "PROCEDURE pq(a:y)\n$COMMON b; c ; d\ny = a + b + c + d\nEND\n");
        INFO(msg);
        CHECK(msg.find("pass b, c, d as input arguments of PROCEDURE pq") != std::string::npos);
    }
    SECTION("mixed case, a string variable in the list, a trailing comment") {
        const std::string msg = commonError(
            "FUNCTION fq(a)\n$Common b, fl$, c   {shared with the main program}\nfq = a + b + c\nEND\n");
        INFO(msg);
        CHECK(msg.find("pass b, fl$, c as input arguments of FUNCTION fq") != std::string::npos);
    }
    SECTION("only a string variable") {
        const std::string msg = commonError(
            "PROCEDURE pq(a:y)\n$common fl$\ny = a\nEND\n");
        INFO(msg);
        CHECK(msg.find("pass fl$ as an input argument of PROCEDURE pq") != std::string::npos);
    }
    SECTION("an array element keeps its subscripts") {
        const std::string msg = commonError(
            "PROCEDURE pq(a:y)\n$common x[1,2], b\ny = a + b\nEND\n");
        INFO(msg);
        CHECK(msg.find("pass x[1,2], b as input arguments of PROCEDURE pq") != std::string::npos);
    }
    SECTION("an empty list still names the directive and the routine") {
        const std::string msg = commonError("PROCEDURE pq(a:y)\n$common\ny = a\nEND\n");
        INFO(msg);
        CHECK(msg.find("$COMMON is not supported") != std::string::npos);
        CHECK(msg.find("PROCEDURE pq") != std::string::npos);
    }
}

TEST_CASE("CS-BUG-COMMON-PROC: every routine with $COMMON is reported", "[library-gaps][common-proc]") {
    EESParser parser;
    auto result = parser.parse(
        "PROCEDURE p1(a:y)\n$common b\ny = a + b\nEND\n"
        "FUNCTION f2(a)\n$common c\nf2 = a + c\nEND\n"
        "PROCEDURE p3(a:y)\ny = a\nEND\n"
        "b = 1\nc = 2\n"
        "CALL p1(b:y1)\nz = f2(c)\nCALL p3(b:y3)\n");
    REQUIRE(result.errors.size() == 2);
    CHECK(result.errors[0].message.find("PROCEDURE p1") != std::string::npos);
    CHECK(result.errors[1].message.find("FUNCTION f2") != std::string::npos);
}

TEST_CASE("CS-BUG-COMMON-PROC: $COMMON inside nested blocks of a body is rejected", "[library-gaps][common-proc]") {
    const std::string msg = commonError(
        "PROCEDURE pq(a:y)\n"
        "IF (a > 0) THEN\n"
        "$common b\n"
        "y = a + b\n"
        "ELSE\n"
        "y = a\n"
        "ENDIF\n"
        "END\n");
    INFO(msg);
    CHECK(msg.find("pass b as an input argument of PROCEDURE pq") != std::string::npos);
}

TEST_CASE("CS-BUG-COMMON-PROC: passing the shared variables as arguments gives the EES values", "[library-gaps][common-proc]") {
    SECTION("PROCEDURE: the register reproducer rewritten") {
        auto run = runModel(
            "PROCEDURE pq(a, b:y)\n"
            "y = a + b\n"
            "END\n"
            "b = 21\n"
            "a = 1\n"
            "CALL pq(a, b:y)\n"
            "z = y*2\n"
            "w = z + b\n");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["y"], WithinAbs(22.0, 1e-9));
        CHECK_THAT(run["z"], WithinAbs(44.0, 1e-9));
        CHECK_THAT(run["w"], WithinAbs(65.0, 1e-9));
    }
    SECTION("PROCEDURE: the product, which the zero-valued $COMMON made 0") {
        auto run = runModel(
            "PROCEDURE pq(a, b:y)\n"
            "y = a*b\n"
            "END\n"
            "b = 21\n"
            "a = 2\n"
            "CALL pq(a, b:y)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(42.0, 1e-9));
    }
    SECTION("FUNCTION") {
        auto run = runModel(
            "FUNCTION fq(a, b)\n"
            "fq = a + b\n"
            "END\n"
            "b = 21\n"
            "a = 1\n"
            "y = fq(a, b)\n"
            "z = y*2\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(22.0, 1e-9));
        CHECK_THAT(run["z"], WithinAbs(44.0, 1e-9));
    }
}

TEST_CASE("CS-BUG-COMMON-PROC: no false positive for $common in a comment or a string", "[library-gaps][common-proc]") {
    SECTION("brace, quote and C-style comments in a body") {
        const std::string code =
            "PROCEDURE pq(a, b:y)\n"
            "{ $common b }\n"
            "\"$common b\"\n"
            "// $common b\n"
            "y = a + b {$common c}\n"
            "END\n"
            "b = 21\n"
            "a = 1\n"
            "CALL pq(a, b:y)\n";
        CHECK(commonError(code).empty());
        auto run = runModel(code);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(22.0, 1e-9));
    }
    SECTION("a line of a multi-line brace comment that starts with $common") {
        const std::string code =
            "FUNCTION fq(a, b)\n"
            "{ the former version used\n"
            "$common b\n"
            "which is not supported }\n"
            "fq = a + b\n"
            "END\n"
            "y = fq(1, 21)\n";
        CHECK(commonError(code).empty());
        auto run = runModel(code);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(22.0, 1e-9));
    }
    SECTION("a line of a multi-line quote comment that starts with $common") {
        const std::string code =
            "FUNCTION fq(a, b)\n"
            "\"\n"
            "$common b\n"
            "\"\n"
            "fq = a + b\n"
            "END\n"
            "y = fq(1, 21)\n";
        CHECK(commonError(code).empty());
    }
    SECTION("a string that contains the word") {
        auto run = runModel(
            "FUNCTION fq(a, b)\n"
            "fq = a + b\n"
            "END\n"
            "note$ = '$common b'\n"
            "y = fq(1, 21)\n");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(22.0, 1e-9));
        CHECK(run.strings.at("note$") == "$common b");
    }
    SECTION("a $IF branch that is switched off") {
        const std::string code =
            "PROCEDURE pq(a, b:y)\n"
            "$IF PARAMETRICTABLE\n"
            "$common b\n"
            "$ENDIF\n"
            "y = a + b\n"
            "END\n"
            "CALL pq(1, 21:y)\n";
        CHECK(commonError(code).empty());
    }
    SECTION("a similar directive name is not $COMMON") {
        const std::string code =
            "PROCEDURE pq(a, b:y)\n"
            "$commons b\n"
            "y = a + b\n"
            "END\n"
            "CALL pq(1, 21:y)\n";
        CHECK(commonError(code).empty());
    }
}

TEST_CASE("CS-BUG-COMMON-PROC: a $COMMON line outside a FUNCTION or PROCEDURE is still accepted", "[library-gaps][common-proc]") {
    // it shares nothing with any routine: only the bodies are concerned
    auto run = runModel(
        "$common b\n"
        "b = 21\n"
        "a = 1\n"
        "y = a + b\n");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(22.0, 1e-9));
}

// ============================================================================
// CS-GAP-NASA: the EES `CALL NASA(X$, T : cp, h, s)` is not a CoolSolve procedure; models use the
// ideal-gas species as fluids (h = MolarMass(X$)*enthalpy(X$,T=T), J/kmol). These tests pin the
// values of that rewrite against the JANAF tables.
//
// Reference: NIST-JANAF Thermochemical Tables, 4th ed. (M. W. Chase, J. Phys. Chem. Ref. Data,
// Monograph 9, 1998), as reproduced by the Shomate fits of the NIST Chemistry WebBook:
// H(2400 K) - H(298.15 K) in kJ/mol N2 70.640, O2 74.453, CO2 115.779, H2O 93.744, CO 71.327,
// H2 66.928; cp(2400 K) in J/mol-K N2 36.510, O2 38.632, CO2 61.277, H2O 53.438, CO 36.702,
// H2 35.559. 2400 K = 2126.85 degC.
// ============================================================================

TEST_CASE("CS-GAP-NASA: molar enthalpy rise of the combustion species against JANAF", "[library-gaps][nasa]") {
    // h = MolarMass(X)*enthalpy(X,T=T) [J/kmol] replaces CALL NASA(X$,T:cp,h,s) (h in kJ/kmol)
    auto run = runModel(R"(
T_hot = 2126.85
T_ref = 25
dH_N2 = MolarMass(N2)*(enthalpy(N2,T=T_hot)-enthalpy(N2,T=T_ref))
dH_O2 = MolarMass(O2)*(enthalpy(O2,T=T_hot)-enthalpy(O2,T=T_ref))
dH_CO2 = MolarMass(CO2)*(enthalpy(CO2,T=T_hot)-enthalpy(CO2,T=T_ref))
dH_H2O = MolarMass(H2O)*(enthalpy(H2O,T=T_hot)-enthalpy(H2O,T=T_ref))
dH_CO = MolarMass(CO)*(enthalpy(CO,T=T_hot)-enthalpy(CO,T=T_ref))
dH_H2 = MolarMass(H2)*(enthalpy(H2,T=T_hot)-enthalpy(H2,T=T_ref))
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // J/kmol, tolerance 0.2 % (CoolProp's ideal-gas parts deviate by 0.01 to 0.10 %)
    CHECK_THAT(run["dH_N2"], WithinRel(70.640e6, 2e-3));
    CHECK_THAT(run["dH_O2"], WithinRel(74.453e6, 2e-3));
    CHECK_THAT(run["dH_CO2"], WithinRel(115.779e6, 2e-3));
    CHECK_THAT(run["dH_H2O"], WithinRel(93.744e6, 2e-3));
    CHECK_THAT(run["dH_CO"], WithinRel(71.327e6, 2e-3));
    CHECK_THAT(run["dH_H2"], WithinRel(66.928e6, 2e-3));
}

TEST_CASE("CS-GAP-NASA: the molar enthalpy at 25 C is the enthalpy of formation, in J/kmol", "[library-gaps][nasa]") {
    // NASA returns h of formation at 298.15 K (elements 0), the reference of the heats of reaction
    auto run = runModel(R"(
hf_CO2 = MolarMass(CO2)*enthalpy(CO2,T=25)
hf_H2O = MolarMass(H2O)*enthalpy(H2O,T=25)
hf_CO = MolarMass(CO)*enthalpy(CO,T=25)
hf_CH4 = MolarMass(CH4)*enthalpy(CH4,T=25)
hf_N2 = MolarMass(N2)*enthalpy(N2,T=25)
hf_O2 = MolarMass(O2)*enthalpy(O2,T=25)
hf_H2 = MolarMass(H2)*enthalpy(H2,T=25)
)");
    REQUIRE(run.solveOk);
    // NASA-Glenn thermo.inp heats of formation (the data CALL NASA reads); Cengel & Boles gave -393.520, -241.820,
    // -110.530, -74.850
    CHECK_THAT(run["hf_CO2"], WithinRel(-393.510e6, 1e-9));
    CHECK_THAT(run["hf_H2O"], WithinRel(-241.826e6, 1e-9));
    CHECK_THAT(run["hf_CO"], WithinRel(-110.535196e6, 1e-9));
    CHECK_THAT(run["hf_CH4"], WithinRel(-74.600e6, 1e-9));
    CHECK_THAT(run["hf_N2"], WithinAbs(0.0, 1.0));
    CHECK_THAT(run["hf_O2"], WithinAbs(0.0, 1.0));
    CHECK_THAT(run["hf_H2"], WithinAbs(0.0, 1.0));
}

TEST_CASE("CS-GAP-NASA: molar specific heat against JANAF", "[library-gaps][nasa]") {
    // cp = MolarMass(X)*cp(X,T=T) [J/kmol-K] replaces the cp output of CALL NASA (kJ/kmol-K)
    auto run = runModel(R"(
T_hot = 2126.85
cp_N2 = MolarMass(N2)*cp(N2,T=T_hot)
cp_O2 = MolarMass(O2)*cp(O2,T=T_hot)
cp_CO2 = MolarMass(CO2)*cp(CO2,T=T_hot)
cp_H2O = MolarMass(H2O)*cp(H2O,T=T_hot)
cp_CO = MolarMass(CO)*cp(CO,T=T_hot)
cp_H2 = MolarMass(H2)*cp(H2,T=T_hot)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["cp_N2"], WithinRel(36.510e3, 5e-3));
    CHECK_THAT(run["cp_O2"], WithinRel(38.632e3, 5e-3));
    CHECK_THAT(run["cp_CO2"], WithinRel(61.277e3, 5e-3));
    CHECK_THAT(run["cp_H2O"], WithinRel(53.438e3, 5e-3));
    CHECK_THAT(run["cp_CO"], WithinRel(36.702e3, 5e-3));
    CHECK_THAT(run["cp_H2"], WithinRel(35.559e3, 5e-3));
}

TEST_CASE("CS-GAP-NASA: the species name from a string variable", "[library-gaps][nasa]") {
    // the CSL-0150 pattern with the species held in string variables; Prod$[i] (a string-array
    // element as fluid argument) is not accepted by CoolSolve, see CS-GAP-STRING-ARRAY
    auto run = runModel(R"(
Sp1$ = 'CO2'
Sp2$ = 'H2O'
T_hot = 2126.85
T_ref = 25
H_ref1 = MolarMass(Sp1$)*enthalpy(Sp1$,T=T_ref)
H_hot1 = MolarMass(Sp1$)*enthalpy(Sp1$,T=T_hot)
H_ref2 = MolarMass(Sp2$)*enthalpy(Sp2$,T=T_ref)
H_hot2 = MolarMass(Sp2$)*enthalpy(Sp2$,T=T_hot)
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["H_ref1"], WithinRel(-393.510e6, 1e-9));   // NASA-Glenn heats of formation
    CHECK_THAT(run["H_hot1"] - run["H_ref1"], WithinRel(115.779e6, 2e-3));
    CHECK_THAT(run["H_ref2"], WithinRel(-241.826e6, 1e-9));
    CHECK_THAT(run["H_hot2"] - run["H_ref2"], WithinRel(93.744e6, 2e-3));
}

TEST_CASE("CS-GAP-NASA: adiabatic flame temperature of CH4 with the rewritten species enthalpies", "[library-gaps][nasa]") {
    // CSL-0150 without dissociation (CH4, 20 % excess air, reactants at 288 K), species enthalpies from
    // property calls. Hand balance with the NASA-Glenn data (thermo.inp: polynomials and heats of formation,
    // computed in Python): H_react = -78.344 MJ/kmol fuel, T_adiab = 2061.0 K = 1787.85 degC. EES stores
    // H_react = -78 333.02 kJ/kmol for the same reactants (Ex_3-4_CH4-Fuel Oil_Tadiab...EES; its CH4 is -74 595 against
    // NASA's -74 600). The Cengel & Boles heats of formation gave -78.601 MJ/kmol and 1788.2 degC (hand balance
    // with the JANAF enthalpies: -78.593 and 1788.5).
    auto run = runModel(R"(
m = 1 ; n = 4 ; e = 0.2
T_fuel = 14.85 ; T_air = 14.85
H_react = molarmass(CH4)*enthalpy(CH4,T=T_fuel)+(m+n/4)*(1+e)*(molarmass(O2)*enthalpy(O2,T=T_air)+(79/21)*molarmass(N2)*enthalpy(N2,T=T_air))
nkmol[1] = m ; nkmol[2] = n/2 ; nkmol[3] = (m+n/4)*e ; nkmol[4] = (m+n/4)*(79/21)*(1+e)
H_prod = nkmol[1]*molarmass(CO2)*enthalpy(CO2,T=T_adiab)+nkmol[2]*molarmass(H2O)*enthalpy(H2O,T=T_adiab)+nkmol[3]*molarmass(O2)*enthalpy(O2,T=T_adiab)+nkmol[4]*molarmass(N2)*enthalpy(N2,T=T_adiab)
H_react = H_prod
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["H_react"], WithinRel(-78.333e6, 5e-4));
    CHECK_THAT(run["T_adiab"], WithinAbs(1787.85, 2.0));   // CoolProp's specific heats: +0.8 K
}

TEST_CASE("CS-GAP-NASA: the species of the NASA database that CoolProp lacks are unknown fluids", "[library-gaps][nasa]") {
    // CS-GAP-FLUIDS-NASA-SPECIES: EES evaluates enthalpy(OH, T=2000) from its NASA database
    for (const std::string species : {"OH", "H", "O", "N", "NO", "NO2"}) {
        auto run = runModel("h = enthalpy(" + species + ", T=2000)\n");
        INFO(species);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("Unknown fluid: '" + species + "'") != std::string::npos);
    }
}

TEST_CASE("CS-GAP-NASA: CALL NASA is not a CoolSolve procedure", "[library-gaps][nasa]") {
    // the model is rewritten with property calls; the call must not be silently ignored
    auto run = runModel(R"(
T = 2000
CALL NASA('CO2',T:CP,H,S)
)");
    CHECK_FALSE(run.solveOk);
    CHECK(run.message.find("Unknown procedure: NASA") != std::string::npos);
}

// ============================================================================
// CS-GAP-STRING-ARRAY: string arrays are rejected with an error (maintainer decision,
// 2026-10-10: not supported). The strings go into a lookup table, the array holds the
// row number, and LOOKUP$('table', row, col) reads the string.
// ============================================================================

namespace {

/// Every error of the parse of `code`: (line, message).
std::vector<std::pair<int, std::string>> strArrayParseErrors(const std::string& code) {
    EESParser parser;
    auto result = parser.parse(code);
    std::vector<std::pair<int, std::string>> errors;
    for (const auto& e : result.errors) errors.emplace_back(e.line, e.message);
    return errors;
}

/// The error says what is not supported and how to rewrite the model.
void checkStringArrayMessage(const std::string& message) {
    INFO(message);
    CHECK(message.find("String arrays are not supported") != std::string::npos);
    CHECK(message.find("lookup table") != std::string::npos);
    CHECK(message.find("LOOKUP$('table', row, col)") != std::string::npos);
}

/// The species table of the CSL-0150 pattern: names in column 1, molar masses in column 2.
const std::map<std::string, std::string> kSpeciesTable = {
    {"Species", "Name,M\nN2,28.0134\nO2,31.9988\nCO2,44.0095\n"}};

}  // namespace

TEST_CASE("CS-GAP-STRING-ARRAY: the register reproducer of the main program is an error, never a silent 0",
          "[library-gaps][string-array]") {
    // was: Prod$[1]='N2' stored as the number 0, then "Unknown fluid: 'Prod'"
    const std::string code =
        "Prod$[1]='N2'\n"
        "M = MolarMass(Prod$[1])\n";
    auto errors = strArrayParseErrors(code);
    REQUIRE(errors.size() == 2);
    CHECK(errors[0].first == 1);   // the assignment
    CHECK(errors[1].first == 2);   // the read as fluid argument
    for (const auto& e : errors) checkStringArrayMessage(e.second);
    CHECK(errors[0].second.find("'Prod$[1]'") != std::string::npos);

    auto run = runModel(code);
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(run.vars.empty());   // nothing is solved, in particular no Prod$[1] = 0
}

TEST_CASE("CS-GAP-STRING-ARRAY: the register reproducer of a PROCEDURE is an error on every element",
          "[library-gaps][string-array]") {
    // was: assignment accepted, then "String variable not found: U$[2]" at the read
    auto errors = strArrayParseErrors(
        "PROCEDURE pi(n:y)\n"          // 1
        "U$[1]='a'\n"                  // 2
        "U$[2]='b'\n"                  // 3
        "w$=U$[n]\n"                   // 4
        "y=n\n"                        // 5
        "END\n"                        // 6
        "CALL pi(2:yv)\n"
        "z=yv\n");
    REQUIRE(errors.size() == 3);
    CHECK(errors[0].first == 2);
    CHECK(errors[1].first == 3);
    CHECK(errors[2].first == 4);
    for (const auto& e : errors) checkStringArrayMessage(e.second);
    CHECK(errors[2].second.find("'U$[n]'") != std::string::npos);
}

TEST_CASE("CS-GAP-STRING-ARRAY: FUNCTION bodies, DUPLICATE, CALL arguments and 2-D elements", "[library-gaps][string-array]") {
    SECTION("FUNCTION body") {
        auto errors = strArrayParseErrors(
            "FUNCTION f(n)\n"
            "L$[1] = 'a'\n"
            "f = n\n"
            "END\n"
            "x = f(1)\n");
        REQUIRE(errors.size() == 1);
        CHECK(errors[0].first == 2);
        checkStringArrayMessage(errors[0].second);
    }
    SECTION("DUPLICATE: the species names of CSL-0150") {
        auto errors = strArrayParseErrors(
            "T = 300\n"                                  // 1
            "Prod$[1] = 'N2' ; Prod$[2] = 'O2'\n"        // 2 (two statements on one line)
            "DUPLICATE i=1,2\n"                          // 3
            "  h[i] = enthalpy(Prod$[i],T=T)\n"          // 4
            "END\n");
        REQUIRE(errors.size() == 3);
        CHECK(errors[0].first == 2);
        CHECK(errors[1].first == 2);
        CHECK(errors[2].first == 4);
        for (const auto& e : errors) checkStringArrayMessage(e.second);
    }
    SECTION("an element as argument and as output of a CALL") {
        auto errors = strArrayParseErrors(
            "PROCEDURE p(a$:b$)\n"   // 1
            "b$ = a$\n"              // 2
            "END\n"                  // 3
            "CALL p(U$[2] : w$)\n"   // 4
            "CALL p('x' : V$[1])\n"); // 5
        REQUIRE(errors.size() == 2);
        CHECK(errors[0].first == 4);
        CHECK(errors[1].first == 5);
        for (const auto& e : errors) checkStringArrayMessage(e.second);
    }
    SECTION("2-D element and an expression as index") {
        auto errors = strArrayParseErrors("S$[1,2] = 'a'\nk = 2\nw$ = S$[k+1]\n");
        REQUIRE(errors.size() == 2);
        CHECK(errors[0].first == 1);
        CHECK(errors[1].first == 3);
        for (const auto& e : errors) checkStringArrayMessage(e.second);
    }
    SECTION("an element inside a condition") {
        auto errors = strArrayParseErrors(
            "PROCEDURE p(n:y)\n"
            "IF (U$[n] = 'a') THEN y = 1 ELSE y = 2\n"
            "END\n"
            "CALL p(1:yv)\n");
        REQUIRE_FALSE(errors.empty());
        CHECK(errors[0].first == 2);
        checkStringArrayMessage(errors[0].second);
    }
}

TEST_CASE("CS-GAP-STRING-ARRAY: scalar string variables and numeric arrays are unaffected", "[library-gaps][string-array]") {
    auto run = runModel(R"(
F$ = 'N2'
Note$ = 'see U$[1] in the manual'
P[1] = 2 ; P[2] = 3
Q = P[2] + P[1]
M = MolarMass(F$)
w = 5 [m]
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["Q"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["M"], WithinRel(molarMassKgPerKmol("N2"), 1e-9));
    CHECK(run.strings.at("F$") == "N2");
    CHECK(run.strings.at("Note$") == "see U$[1] in the manual");   // a '$[' inside a string literal is text
}

// ----------------------------------------------------------------------------
// LOOKUP$
// ----------------------------------------------------------------------------

TEST_CASE("CS-GAP-LOOKUPROW: LookupTable keeps the text of the cells", "[library-gaps][lookup-string]") {
    SECTION("a table with text cells keeps every cell as written in the CSV") {
        auto t = LookupTable::fromCSV("t", "Name,x,y\nN2, 1.50 ,1e3\nO2,2,\n\"a,b\",3,4\n");
        REQUIRE(t.numRows() == 3);
        CHECK(t.text(1, 1) == "N2");
        CHECK(t.text(1, 2) == "1.50");    // a numeric cell: its text as written, whitespace trimmed
        CHECK(t.text(1, 3) == "1e3");
        CHECK(t.text(2, 3) == "");        // an empty cell
        CHECK(t.text(3, 1) == "a,b");     // a quoted cell
        // the numeric view is unchanged
        CHECK_THAT(t.value(1, 2), WithinAbs(1.5, 1e-15));
        CHECK_THAT(t.value(1, 3), WithinAbs(1000.0, 1e-12));
        CHECK(std::isnan(t.value(1, 1)));
        CHECK(std::isnan(t.value(2, 3)));
        // toCSV keeps the text cells
        auto back = LookupTable::fromCSV("t", t.toCSV());
        CHECK(back.text(1, 1) == "N2");
        CHECK(back.text(3, 1) == "a,b");
    }
    SECTION("a purely numeric table gives the number as text") {
        auto t = LookupTable::fromCSV("n", "x,y\n1,0.5\n2,\n");
        CHECK(t.text(1, 1) == "1");
        CHECK(t.text(1, 2) == "0.5");
        CHECK(t.text(2, 2) == "");
    }
    SECTION("out of range") {
        auto t = LookupTable::fromCSV("t", "Name,x\nN2,1\n");
        CHECK_THROWS_AS(t.text(2, 1), std::out_of_range);
        CHECK_THROWS_AS(t.text(1, 3), std::out_of_range);
        CHECK_THROWS_AS(t.text(0, 1), std::out_of_range);
    }
}

TEST_CASE("CS-GAP-LOOKUPROW: LOOKUP$ returns the text of a cell, the column as a number or a name",
          "[library-gaps][lookup-string]") {
    auto run = runModel(R"(
PROCEDURE cell(r, c : s$)
  s$ = LOOKUP$('Species', r, c)
END
PROCEDURE cellname(r : s$)
  s$ = LOOKUP$('Species', r, 'Name')
END
CALL cell(1, 1 : A$)
CALL cell(3, 1 : B$)
CALL cell(2, 2 : C$)
CALL cellname(2 : D$)
)", kSpeciesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.strings.at("A$") == "N2");
    CHECK(run.strings.at("B$") == "CO2");
    CHECK(run.strings.at("C$") == "31.9988");   // a numeric cell: its text as written
    CHECK(run.strings.at("D$") == "O2");        // the column by name, case of the table
}

TEST_CASE("CS-GAP-LOOKUPROW: LOOKUP$ of a numeric cell keeps the text as written in the CSV", "[library-gaps][lookup-string]") {
    auto run = runModel(R"(
PROCEDURE cell(r, c : s$)
  s$ = LOOKUP$('t', r, c)
END
CALL cell(1, 2 : A$)
CALL cell(2, 2 : B$)
CALL cell(3, 2 : C$)
x = stringval(A$) + stringval(B$)
)", {{"t", "name,v\n a ,1.50\nb,1e3\nc,\n"}});
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.strings.at("A$") == "1.50");
    CHECK(run.strings.at("B$") == "1e3");
    CHECK(run.strings.at("C$") == "");          // an empty cell is the empty string
    CHECK_THAT(run["x"], WithinAbs(1001.5, 1e-9));
}

TEST_CASE("CS-GAP-LOOKUPROW: LOOKUP$ as the fluid argument, the CSL-0150 pattern in a DUPLICATE",
          "[library-gaps][lookup-string]") {
    // Prod$[i] of the EES source becomes LOOKUP$('Species', i, 1)
    auto run = runModel(R"(
T = 300
DUPLICATE i=1,3
  H[i] = MolarMass(LOOKUP$('Species',i,1))*enthalpy(LOOKUP$('Species',i,1),T=T)
END
Href[1] = MolarMass(N2)*enthalpy(N2,T=T)
Href[2] = MolarMass(O2)*enthalpy(O2,T=T)
Href[3] = MolarMass(CO2)*enthalpy(CO2,T=T)
)", kSpeciesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (int i = 1; i <= 3; ++i) {
        INFO("species " << i);
        CHECK_THAT(run["H[" + std::to_string(i) + "]"], WithinRel(run["Href[" + std::to_string(i) + "]"], 1e-12));
    }
}

TEST_CASE("CS-GAP-LOOKUPROW: LOOKUP$ as the fluid of the other property functions", "[library-gaps][lookup-string]") {
    auto run = runModel(R"(
T = 300
P = 1E5
a = cp(LOOKUP$('Species',2,1),T=T)
a_ref = cp(O2,T=T)
b = MolarMass(LOOKUP$('Species',3,1))
b_ref = MolarMass(CO2)
c = entropy(LOOKUP$('Species',1,1),T=T,P=P)
c_ref = entropy(N2,T=T,P=P)
d = density(LOOKUP$('Species',1,1),T=T,P=P)
d_ref = density(N2,T=T,P=P)
tc = T_crit(LOOKUP$('Species',3,1))
tc_ref = T_crit(CO2)
)", kSpeciesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinRel(run["a_ref"], 1e-12));
    CHECK_THAT(run["b"], WithinRel(run["b_ref"], 1e-12));
    CHECK_THAT(run["c"], WithinRel(run["c_ref"], 1e-12));
    CHECK_THAT(run["d"], WithinRel(run["d_ref"], 1e-12));
    CHECK_THAT(run["tc"], WithinRel(run["tc_ref"], 1e-12));
}

TEST_CASE("CS-GAP-LOOKUPROW: the row of LOOKUP$ is a variable computed by another equation", "[library-gaps][lookup-string]") {
    // the equation that uses the fluid must wait for idx (a dependency hidden in the fluid argument
    // would evaluate it with the default value 1, i.e. the wrong species, silently)
    auto run = runModel(R"(
H = enthalpy(LOOKUP$('Species', idx, 1), T=300)
idx = k + 1
k = 1
H_ref = enthalpy(O2, T=300)
)", kSpeciesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["idx"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["H"], WithinRel(run["H_ref"], 1e-12));
}

TEST_CASE("CS-GAP-LOOKUPROW: an unknown solved through a property call of the LOOKUP$ fluid", "[library-gaps][lookup-string]") {
    // the temperature at which O2 (row 2 of the table) has the molar enthalpy of 800 degC
    auto run = runModel(R"(
Href = MolarMass(O2)*enthalpy(O2,T=800)
Href = MolarMass(LOOKUP$('Species',2,1))*enthalpy(LOOKUP$('Species',2,1),T=Tx)
)", kSpeciesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["Tx"], WithinAbs(800.0, 1e-6));
}

TEST_CASE("CS-GAP-LOOKUPROW: LOOKUP$ in FUNCTION and PROCEDURE bodies", "[library-gaps][lookup-string]") {
    SECTION("FUNCTION: fluid argument and string comparison") {
        auto run = runModel(R"(
FUNCTION hm(i, T)
  hm = MolarMass(LOOKUP$('Species',i,1))*enthalpy(LOOKUP$('Species',i,1),T=T)
END
FUNCTION isO2(r)
  isO2 = 0
  IF (LOOKUP$('Species', r, 1) = 'o2') THEN isO2 = 1
END
FUNCTION isgas(i)
  isgas = 0
  IF (PHASE$(LOOKUP$('Species', i, 1), T=25, P=1E5) = 'gas') THEN isgas = 1
END
a = hm(3, 400)
a_ref = MolarMass(CO2)*enthalpy(CO2,T=400)
b = isO2(2)
c = isO2(1)
d = isgas(3)
)", kSpeciesTable);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinRel(run["a_ref"], 1e-12));
        CHECK_THAT(run["b"], WithinAbs(1.0, 1e-12));   // the comparison ignores the case, as EES does
        CHECK_THAT(run["c"], WithinAbs(0.0, 1e-12));
        CHECK_THAT(run["d"], WithinAbs(1.0, 1e-12));   // PHASE$ takes the LOOKUP$ fluid too
    }
    SECTION("PROCEDURE: a string variable of the body, the BrineProp pattern") {
        // UO$ = U$[Pro] of the BrineProp library becomes UO$ = LOOKUP$('Units', Pro, 1)
        auto run = runModel(R"(
PROCEDURE unitlabel(Pro : UO$, n)
  UO$ = LOOKUP$('Units', Pro, 1)
  n = Pro
END
CALL unitlabel(3 : A$, n3)
CALL unitlabel(5 : B$, n5)
)", {{"Units", "Unit\nC\nkg/m^3\nkJ/kg-K\nW/m-K\nmilliPa-s\n"}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.strings.at("A$") == "kJ/kg-K");
        CHECK(run.strings.at("B$") == "milliPa-s");
    }
    SECTION("table name from a string variable") {
        auto run = runModel(R"(
PROCEDURE pick(tab$, r : s$)
  s$ = LOOKUP$(tab$, r, 1)
END
CALL pick('Species', 2 : A$)
)", kSpeciesTable);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.strings.at("A$") == "O2");
    }
}

TEST_CASE("CS-GAP-LOOKUPROW: LOOKUP$ negative tests give clear errors", "[library-gaps][lookup-string]") {
    auto fail = [&](const std::string& call, const std::string& expected) {
        auto run = runModel("PROCEDURE p(r : s$)\n  s$ = " + call + "\nEND\nCALL p(1 : A$)\n", kSpeciesTable);
        INFO(call);
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find(expected) != std::string::npos);
    };
    SECTION("row out of range") {
        fail("LOOKUP$('Species', 4, 1)", "LOOKUP$(): row 4 is outside the table 'Species' (1 to 3");
        fail("LOOKUP$('Species', 0, 1)", "LOOKUP$(): row 0 is outside the table 'Species' (1 to 3");
        fail("LOOKUP$('Species', -2, 1)", "row -2 is outside");
    }
    SECTION("column out of range or unknown") {
        fail("LOOKUP$('Species', 1, 3)", "LOOKUP$(): column 3 is outside the table 'Species' (1 to 2");
        fail("LOOKUP$('Species', 1, 0)", "column 0 is outside");
        fail("LOOKUP$('Species', 1, 'Nope')", "LOOKUP$(): column 'Nope' not found in table 'Species'");
    }
    SECTION("unknown table") {
        fail("LOOKUP$('Nope', 1, 1)", "LOOKUP$(): lookup table 'Nope' not found");
    }
    SECTION("arity") {
        fail("LOOKUP$('Species', 1)", "LOOKUP$(): expects 3 arguments: (table, row, col)");
    }
    SECTION("used as a number") {
        auto run = runModel("x = LOOKUP$('Species', 1, 2)\n", kSpeciesTable);
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("LOOKUP$ returns a string") != std::string::npos);
    }
    SECTION("an unknown fluid text from the table is reported as such") {
        auto run = runModel("h = enthalpy(LOOKUP$('Bad', 1, 1), T=300)\n", {{"Bad", "Name\nUnobtainium\n"}});
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("Unknown fluid: 'Unobtainium'") != std::string::npos);
    }
}

TEST_CASE("CS-GAP-LOOKUPROW: a numeric LOOKUP next to a text column keeps its value and its derivative",
          "[library-gaps][lookup-string]") {
    // 0*NaN: the text column must not spoil the cell of the numeric column (was NaN)
    const std::string csv = "name,v\nA,10\nB,20\nC,40\n";
    auto exact = evalWithSlope("y = lookup('t', r, 2)\n", csv, 2.0);
    CHECK_THAT(exact.value, WithinAbs(20.0, 1e-12));
    CHECK_THAT(exact.slope, WithinAbs(20.0, 1e-12));   // the slope of the segment 2..3 (40 - 20)
    auto between = evalWithSlope("y = lookup('t', r, 'v')\n", csv, 2.5);
    CHECK_THAT(between.value, WithinAbs(30.0, 1e-12));
    CHECK_THAT(between.slope, WithinAbs(20.0, 1e-12));
    auto first = evalWithSlope("y = lookup('t', r, 2)\n", csv, 1.0);
    CHECK_THAT(first.value, WithinAbs(10.0, 1e-12));
    CHECK_THAT(first.slope, WithinAbs(10.0, 1e-12));
}

TEST_CASE("CS-GAP-LOOKUPROW: the numeric lookup functions are unchanged on a table that has text cells",
          "[library-gaps][lookup-string]") {
    auto run = runModel(R"(
a = lookup('Species', 2, 2)
b = lookup('Species', 2.5, 'M')
c = nlookuprows('Species')
d = nlookupcolumns('Species')
e = sumlookup('Species', 'M')
f = interpolate('Species', 'M', 'M', 30)
)", kSpeciesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(31.9988, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(0.5 * (31.9988 + 44.0095), 1e-12));
    CHECK_THAT(run["c"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["d"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["e"], WithinAbs(28.0134 + 31.9988 + 44.0095, 1e-9));
    CHECK_THAT(run["f"], WithinAbs(30.0, 1e-9));
}

// ============================================================================
// CS-GAP-CALL-WARNING: CALL WARNING('message XXXA1', value) in functions and procedures
// ============================================================================

namespace {

/// The warnings raised by CALL WARNING (code W001), in the order reported.
std::vector<Diagnostic> userWarnings(const ModelRun& run) {
    std::vector<Diagnostic> out;
    for (const auto& d : run.warningItems)
        if (d.code == "W001") out.push_back(d);
    return out;
}

/// The messages of the CALL WARNING warnings of `run`.
std::vector<std::string> userWarningMessages(const ModelRun& run) {
    std::vector<std::string> out;
    for (const auto& d : userWarnings(run)) out.push_back(d.message);
    return out;
}

}  // namespace

TEST_CASE("CS-GAP-CALL-WARNING: the EES help example warns once at the solution and the calculation continues",
          "[library-gaps][call-warning]") {
    // warning_procedure.htm of the EES help: the sum is -1, so the branch is taken
    const std::string model = R"(
Function abc(X,Y)
  G:=X+Y
  If (G<=0) Then Call warning('The sum of X and Y should be greater than zero.  The sum was XXXA1 ',G)
  abc:=G
End
g=abc(3,-4)
)";
    auto run = runModel(model);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["g"], WithinAbs(-1.0, 1e-12));
    CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
    const auto warnings = userWarnings(run);
    REQUIRE(warnings.size() == 1);
    CHECK(warnings[0].message == "The sum of X and Y should be greater than zero.  The sum was -1 ");
    CHECK(warnings[0].line == 4);   // the line of the CALL WARNING statement

    // the branch is not taken: nothing is reported
    auto quiet = runModel(model.substr(0, model.rfind("g=abc(3,-4)")) + "g=abc(3,4)\n");
    REQUIRE(quiet.solveOk);
    CHECK_THAT(quiet["g"], WithinAbs(7.0, 1e-12));
    CHECK(userWarnings(quiet).empty());
}

TEST_CASE("CS-GAP-CALL-WARNING: the value placeholders XXX, XXXA1, XXXF1, XXXE3, XXXN4, XXXG3 and XXX$",
          "[library-gaps][call-warning]") {
    auto run = runModel(R"(
FUNCTION fmt(v, s$)
  CALL WARNING('plain=XXX|', v)
  CALL WARNING('auto=XXXA3|', v)
  CALL WARNING('fixed1=XXXF1|', v)
  CALL WARNING('fixed0=XXXF0|', v)
  CALL WARNING('exp=XXXE2|', v)
  CALL WARNING('sig=XXXN4|', v)
  CALL WARNING('eng=XXXG3|', v)
  CALL WARNING('string=XXX$|', s$)
  CALL WARNING('no value, text only')
  fmt = v
END
r$ = 'R134a'
x = fmt(123.4567, r$)
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    const auto messages = userWarningMessages(run);
    REQUIRE(messages.size() == 9);
    CHECK(messages[0] == "plain=123.4567|");
    CHECK(messages[1] == "auto=123.4567|");
    CHECK(messages[2] == "fixed1=123.5|");
    CHECK(messages[3] == "fixed0=123|");
    CHECK(messages[4] == "exp=1.23E2|");
    CHECK(messages[5] == "sig=123.5|");
    CHECK(messages[6] == "eng=123E0|");
    CHECK(messages[7] == "string=R134a|");
    CHECK(messages[8] == "no value, text only");
}

TEST_CASE("CS-GAP-CALL-WARNING: a message without placeholder keeps the value, a numeric argument alone gets the generic text",
          "[library-gaps][call-warning]") {
    auto run = runModel(R"(
Function abc(X)
  If (X>0) Then Call warning('X is positive', X)
  If (X>1) Then Call warning(X)
  abc = X
End
g = abc(3)
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    const auto messages = userWarningMessages(run);
    REQUIRE(messages.size() == 2);
    CHECK(messages[0] == "X is positive [3]");   // no XXX in the text: the value is not lost
    CHECK(messages[1] == "A warning message was issued due to the value 3 in abc.");   // EES generic text
}

TEST_CASE("CS-GAP-CALL-WARNING: the range checker of the ULiege air_nozzle_ejector, as written there (lower-case xxx)",
          "[library-gaps][call-warning]") {
    // verbatim from the library model CSL-0118: `if(...) or (...) then call warning(' ... xxxE4 ...', var)`
    auto run = runModel(R"(
procedure warning_error(var1_in,var2_in : var1_out,var2_out)
if(var1_in < 0.18) or (var1_in  > 0.46) then call warning(' FF should be between 0.18 and 0.46.  A value of  xxxE4 was supplied',var1_in)
if(var2_in< 126) or (var2_in  > 265) then call warning(' FAD_tot  should be between 126 and 265.  A value of  xxxE4 was supplied',var2_in)
var1_out := var1_in
var2_out := var2_in
end
FF = 0.1
FAD_tot = 260
call warning_error(FF,FAD_tot:var1_out,var2_out)
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    // only FF is outside its range
    CHECK(userWarningMessages(run) ==
          std::vector<std::string>{" FF should be between 0.18 and 0.46.  A value of  1.0000E-1 was supplied"});
}

TEST_CASE("CS-GAP-CALL-WARNING: CALL WARNING in a PROCEDURE and in routines called by a routine",
          "[library-gaps][call-warning]") {
    auto run = runModel(R"(
FUNCTION inner(v)
  IF (v > 5) THEN CALL WARNING('inner: v = XXXF1', v)
  inner = 2*v
END
PROCEDURE outer(a : b, c)
  b = inner(a)
  c = a + 1
  IF (c > 3) THEN CALL WARNING('outer: c = XXXF2', c)
END
a = 4
CALL outer(a : b, c)
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["b"], WithinAbs(8.0, 1e-12));
    const auto messages = userWarningMessages(run);
    REQUIRE(messages.size() == 1);   // inner does not warn (4 <= 5), outer does (5 > 3)
    CHECK(messages[0] == "outer: c = 5.00");

    auto both = runModel(R"(
FUNCTION inner(v)
  IF (v > 5) THEN CALL WARNING('inner: v = XXXF1', v)
  inner = 2*v
END
PROCEDURE outer(a : b)
  b = inner(a)
END
a = 6
CALL outer(a : b)
)");
    REQUIRE(both.solveOk);
    const auto m2 = userWarningMessages(both);
    REQUIRE(m2.size() == 1);   // a warning raised two calls deep
    CHECK(m2[0] == "inner: v = 6.0");
}

TEST_CASE("CS-GAP-CALL-WARNING: reported once, at the solution, not once per iteration",
          "[library-gaps][call-warning]") {
    // x + f(x) = 12 is solved by iteration (x = 3 from the guess 1); x > 2.5 holds at the solution and at
    // several of the iterates: one message, with the value of the solution
    auto run = runModel(R"(
FUNCTION f(x)
  IF (x > 2.5) THEN CALL WARNING('x is above 2.5: XXXF4', x)
  f = x^2
END
x + f(x) = 12
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(3.0, 1e-6));
    const auto messages = userWarningMessages(run);
    REQUIRE(messages.size() == 1);
    CHECK(messages[0] == "x is above 2.5: 3.0000");
}

TEST_CASE("CS-GAP-CALL-WARNING: a branch taken by an iterate but not at the solution is not reported",
          "[library-gaps][call-warning]") {
    // The solver starts from the guess x = 1, where x < 2 holds; the solution x = 3 is outside the branch.
    auto run = runModel(R"(
FUNCTION f(x)
  IF (x < 2) THEN CALL WARNING('x is below 2: XXXF2', x)
  f = x^2
END
x + f(x) = 12
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["x"], WithinAbs(3.0, 1e-6));
    CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
    CHECK(userWarnings(run).empty());
}

TEST_CASE("CS-GAP-CALL-WARNING: identical warnings are reported once, distinct ones each", "[library-gaps][call-warning]") {
    const std::string routine = R"(
FUNCTION g(v)
  IF (v > 1) THEN CALL WARNING('v is large: XXXF1', v)
  g = v
END
)";
    auto same = runModel(routine + R"(
DUPLICATE i = 1, 3
  y[i] = g(5)
END
)");
    INFO(same.message);
    REQUIRE(same.solveOk);
    CHECK(userWarningMessages(same) == std::vector<std::string>{"v is large: 5.0"});

    auto distinct = runModel(routine + R"(
DUPLICATE i = 1, 3
  y[i] = g(5*i)
END
)");
    REQUIRE(distinct.solveOk);
    CHECK(userWarningMessages(distinct) ==
          std::vector<std::string>{"v is large: 5.0", "v is large: 10.0", "v is large: 15.0"});
}

TEST_CASE("CS-GAP-CALL-WARNING: a warning is reported whatever the $Warnings directive says (EES help)",
          "[library-gaps][call-warning]") {
    // warning_procedure.htm: "always displayed no matter what the setting ... or whether or not a
    // $Warnings On directive is provided"; CoolSolve recognises $Warnings but does not interpret it
    auto run = runModel(R"(
$Warnings OFF
FUNCTION f(x)
  IF (x > 0) THEN CALL WARNING('x = XXXA1', x)
  f = x
END
y = f(2)
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(userWarningMessages(run) == std::vector<std::string>{"x = 2"});
}

TEST_CASE("CS-GAP-CALL-WARNING: CALL ERROR uses the same message text as CALL WARNING", "[library-gaps][call-warning]") {
    const std::string routine = R"(
FUNCTION g(x)
  IF (x < 0) THEN CALL error('x must be positive: XXXF1 was given', x)
  g = sqrt(x)
END
)";
    auto placeholder = runModel(routine + "a = g(-4)\n");
    REQUIRE_FALSE(placeholder.solveOk);
    CHECK(placeholder.message.find("CALL ERROR: x must be positive: -4.0 was given") != std::string::npos);
    CHECK(placeholder.message.find("was given [") == std::string::npos);   // substituted, not appended

    auto generic = runModel(R"(
FUNCTION g(x)
  IF (x < 0) THEN CALL error(x)
  g = sqrt(x)
END
a = g(-4)
)");
    REQUIRE_FALSE(generic.solveOk);
    CHECK(generic.message.find("Calculations have been halted because a parameter is out of range.  "
                               "The value of the parameter is -4 in function g.") != std::string::npos);

    auto exponential = runModel(R"(
FUNCTION g(x)
  IF (x > 0) THEN CALL error('value XXXE3!', x)
  g = x
END
a = g(512.3)
)");
    REQUIRE_FALSE(exponential.solveOk);
    CHECK(exponential.message.find("CALL ERROR: value 5.123E2!") != std::string::npos);   // 5.123E2 = 512.3 (EES help, Display Formats)
}

TEST_CASE("CS-GAP-CALL-WARNING: malformed CALL WARNING and CALL ERROR are parse errors", "[library-gaps][call-warning]") {
    const auto inBody = [](const std::string& call) {
        return parseErrors("FUNCTION f(x)\n  " + call + "\n  f = x\nEND\ny = f(1)\n");
    };
    auto expectError = [](const std::vector<std::string>& errors, const std::string& text) {
        INFO("expected: " << text);
        bool found = false;
        for (const auto& e : errors) {
            INFO(e);
            if (e.find(text) != std::string::npos) found = true;
        }
        CHECK(found);
    };
    SECTION("no argument") {
        expectError(inBody("CALL WARNING()"), "CALL WARNING needs an argument");
        expectError(inBody("CALL ERROR()"), "CALL ERROR needs an argument");
    }
    SECTION("a message that is not a string, followed by a value") {
        expectError(inBody("CALL WARNING(x, 2)"), "the first argument of CALL WARNING must be the message");
        expectError(inBody("CALL ERROR(x + 1, 2)"), "the first argument of CALL ERROR must be the message");
    }
    SECTION("too many arguments") {
        expectError(inBody("CALL WARNING('a', x, 2)"), "CALL WARNING takes a message and a value at most");
    }
    SECTION("outputs") {
        expectError(inBody("CALL WARNING('a', x : y)"), "CALL WARNING has no outputs");
    }
    SECTION("the main program") {
        expectError(parseErrors("x = 1\nCALL WARNING('a', x)\n"),
                    "CALL WARNING can only be used inside a FUNCTION or PROCEDURE");
        expectError(parseErrors("x = 1\nCALL ERROR('a', x)\n"),
                    "CALL ERROR can only be used inside a FUNCTION or PROCEDURE");
    }
    SECTION("the valid forms parse") {
        CHECK(inBody("CALL WARNING('a')").empty());
        CHECK(inBody("CALL WARNING('a XXXA1', x)").empty());
        CHECK(inBody("CALL WARNING(x)").empty());
        // s$ is assigned first: a body reads no local it has not assigned (CS-GAP-STRING-LOCAL-UNSET)
        CHECK(inBody("s$ = 'b'\n  CALL warning(s$)").empty());
        CHECK(inBody("s$ = 'b'\n  CALL warning('a XXX$', s$)").empty());
        CHECK(inBody("CALL ERROR('')").empty());
    }
}

TEST_CASE("CS-GAP-CALL-WARNING: display formats of the value (EES help, Display Formats and STRING$)",
          "[library-gaps][call-warning]") {
    // the examples of the EES help
    CHECK(eesFormatNumber(9.1234565, "F1") == "9.1");        // STRING$ help
    CHECK(eesFormatNumber(1234567.0, "F0") == "1234567");    // STRING$ help
    CHECK(eesFormatNumber(123.4567, "N3") == "123");         // Display Formats help
    CHECK(eesFormatNumber(123.4567, "N4") == "123.5");
    CHECK(eesFormatNumber(512.3, "E3") == "5.123E2");
    // the other formats
    CHECK(eesFormatNumber(1234567.0, "N3") == "1.23E6");     // outside 1E-6..1E6: E format
    CHECK(eesFormatNumber(1.23456e-7, "N3") == "1.23E-7");
    CHECK(eesFormatNumber(1.23456e-5, "N3") == "0.0000123");
    CHECK(eesFormatNumber(0.012345, "N3") == "0.0123");
    CHECK(eesFormatNumber(99.96, "N3") == "100");            // the rounding carries into the next digit
    CHECK(eesFormatNumber(-3.0, "E1") == "-3.0E0");
    CHECK(eesFormatNumber(12345.678, "G3") == "12.3E3");
    CHECK(eesFormatNumber(0.00123, "G2") == "1.2E-3");
    CHECK(eesFormatNumber(0.0, "F2") == "0.00");
    CHECK(eesFormatNumber(123.456, "FM") == "123.456");      // M: the most digits
    CHECK(eesFormatNumber(123.456, "f1") == "123.5");        // lower case
    // auto: whole numbers without decimals, else up to 10 significant digits; also the unknown codes
    CHECK(eesFormatNumber(3.0) == "3");
    CHECK(eesFormatNumber(-4.0, "A1") == "-4");
    CHECK(eesFormatNumber(1.5, "A3") == "1.5");
    CHECK(eesFormatNumber(0.1 + 0.2) == "0.3");
    CHECK(eesFormatNumber(1.0 / 3.0) == "0.3333333333");
    CHECK(eesFormatNumber(2.5, "T0") == "2.5");
    CHECK(eesFormatNumber(2.5, "") == "2.5");
}

TEST_CASE("CS-GAP-CALL-WARNING: the placeholder XXX and the format code that may follow", "[library-gaps][call-warning]") {
    const auto tag = [](const std::string& format) { return "<" + format + ">"; };
    CHECK(eesSubstituteValue("a XXXA1 b", tag) == "a <A1> b");
    CHECK(eesSubstituteValue("a XXX b", tag) == "a <> b");
    CHECK(eesSubstituteValue("v=XXXF3.", tag) == "v=<F3>.");
    CHECK(eesSubstituteValue("v=XXXe4", tag) == "v=<e4>");
    CHECK(eesSubstituteValue("A value of  xxxE4 was supplied", tag) == "A value of  <E4> was supplied");   // library files write xxx
    CHECK(eesSubstituteValue("a Xxx b", tag) == "a <> b");
    CHECK(eesSubstituteValue("s=XXX$.", tag) == "s=<$>.");
    CHECK(eesSubstituteValue("XXXE4 and XXXE4", tag) == "<E4> and <E4>");
    CHECK(eesSubstituteValue("XXXAverage", tag) == "<>Average");   // a letter that is not followed by a digit is text
    CHECK(eesSubstituteValue("XXXFM.", tag) == "<FM>.");
    CHECK(eesSubstituteValue("XXXFMAX", tag) == "<>FMAX");
    bool found = true;
    CHECK(eesSubstituteValue("no placeholder, XX only", tag, &found) == "no placeholder, XX only");
    CHECK_FALSE(found);
    CHECK(eesSubstituteValue("XXXA1", tag, &found) == "<A1>");   // at the very end of the text
    CHECK(found);
}

TEST_CASE("CS-GAP-CALL-WARNING: only the evaluation at the solution reports, only blocks that can warn evaluate again",
          "[library-gaps][call-warning]") {
    const std::string routines = R"(
FUNCTION f(x)
  IF (x < 2) THEN CALL WARNING('x is below 2: XXXF2', x)
  f = x^2
END
FUNCTION viaF(x)
  viaF = 3*f(x)
END
)";
    // build the evaluator of a model: its blocks, found by the unknown they solve
    struct Built {
        IR ir;
        StructuralAnalysisResult analysis;
        size_t blockOf(const std::string& variable) const {
            for (size_t i = 0; i < analysis.blocks.size(); ++i)
                for (const auto& v : analysis.blocks[i].variables)
                    if (v == variable) return i;
            FAIL("no block solves " << variable);
            return 0;
        }
    };
    const auto build = [](const std::string& code) {
        EESParser parser;
        auto parsed = parser.parse(code);
        REQUIRE(parsed.success);
        Built built{IR::fromAST(parsed.program), {}};
        built.analysis = StructuralAnalyzer::analyze(built.ir);
        return built;
    };

    SECTION("a block that calls a routine that warns, directly or through another one") {
        const auto m = build(routines + "x + f(x) = 12\nz = 4\nw = viaF(z) + 1\nv = 2*z\n");
        SystemEvaluator ev(m.ir, m.analysis);
        CHECK(ev.getBlock(m.blockOf("x")).mayCallWarning());
        CHECK(ev.getBlock(m.blockOf("w")).mayCallWarning());   // viaF calls f, which warns
        CHECK_FALSE(ev.getBlock(m.blockOf("z")).mayCallWarning());
        CHECK_FALSE(ev.getBlock(m.blockOf("v")).mayCallWarning());
    }

    SECTION("an iterate takes the branch and nothing is reported; the solution does not and nothing is") {
        const auto m = build(routines + "x + f(x) = 12\n");
        SystemEvaluator ev(m.ir, m.analysis);
        const size_t block = m.blockOf("x");

        // the guess x = 1 of the first iteration is inside the branch: the evaluation of an iterate
        // (no sink) leaves no trace, while a sink collects the warning it would raise
        DiagnosticCollector diagnostics;
        ev.getBlock(block).setDiagnostics(&diagnostics);
        ev.setVariableValue("x", 1.0);
        ev.evaluateBlock(block);
        CHECK(diagnostics.empty());
        std::vector<Diagnostic> raised;
        ev.evaluateBlock(block, false, &raised);
        REQUIRE(raised.size() == 1);
        CHECK(raised[0].code == "W001");
        CHECK(raised[0].severity == DiagnosticSeverity::Warning);
        CHECK(raised[0].message == "x is below 2: 1.00");

        // reported at a point of the branch, once however often it is asked ...
        DiagnosticCollector out;
        ev.reportCallWarnings(block, out);
        ev.reportCallWarnings(block, out);
        REQUIRE(out.size() == 1);
        CHECK(out.items()[0].message == "x is below 2: 1.00");
        CHECK(diagnostics.empty());   // the diagnostics of the iterations are not used

        // ... and nothing at the solution x = 3
        DiagnosticCollector atSolution;
        ev.setVariableValue("x", 3.0);
        ev.reportCallWarnings(block, atSolution);
        CHECK(atSolution.empty());
    }

    SECTION("a model without CALL WARNING marks no block, whatever routines it has") {
        const auto m = build(R"(
FUNCTION f(x)
  IF (x < 2) THEN CALL ERROR('x is below 2: XXXF2', x)
  f = x^2
END
x + f(x) = 12
y = 2*x
)");
        SystemEvaluator ev(m.ir, m.analysis);
        for (size_t i = 0; i < ev.getNumBlocks(); ++i) CHECK_FALSE(ev.getBlock(i).mayCallWarning());
        DiagnosticCollector out;
        for (size_t i = 0; i < ev.getNumBlocks(); ++i) ev.reportCallWarnings(i, out);
        CHECK(out.empty());
    }
}

// ============================================================================
// CS-GAP-NDAY: the EES date function nDay_(month, day), the day of the year
// ============================================================================

TEST_CASE("CS-GAP-NDAY: reproducer of the register", "[library-gaps][nday]") {
    // From the EES stored solution of library model CSL-0119: nDay_(12, 4) = 338
    auto run = runModel(R"(
month = 12
day_ = 4
n = nDay_(month, day_)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["n"], WithinAbs(338.0, 1e-12));
    CHECK_FALSE(anyContains(run.warnings, "Unknown function 'nDay_'"));
}

TEST_CASE("CS-GAP-NDAY: day of a non-leap year, month by month", "[library-gaps][nday]") {
    auto run = runModel(R"(
a = nDay_(1, 1)
b = nDay_(2, 28)
c = nDay_(3, 1)
d = nDay_(12, 31)
n01 = nDay_(1, 31)
n02 = nDay_(2, 1)
n04 = nDay_(4, 1)
n05 = nDay_(5, 1)
n06 = nDay_(6, 1)
n07 = nDay_(7, 1)
n08 = nDay_(8, 1)
n09 = nDay_(9, 1)
n10 = nDay_(10, 1)
n11 = nDay_(11, 1)
n12 = nDay_(12, 1)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(59.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(60.0, 1e-12));      // 29 February does not exist
    CHECK_THAT(run["d"], WithinAbs(365.0, 1e-12));
    // first day of each month = cumulative length of the months before + 1
    CHECK_THAT(run["n01"], WithinAbs(31.0, 1e-12));
    CHECK_THAT(run["n02"], WithinAbs(32.0, 1e-12));
    CHECK_THAT(run["n04"], WithinAbs(91.0, 1e-12));
    CHECK_THAT(run["n05"], WithinAbs(121.0, 1e-12));
    CHECK_THAT(run["n06"], WithinAbs(152.0, 1e-12));
    CHECK_THAT(run["n07"], WithinAbs(182.0, 1e-12));
    CHECK_THAT(run["n08"], WithinAbs(213.0, 1e-12));
    CHECK_THAT(run["n09"], WithinAbs(244.0, 1e-12));
    CHECK_THAT(run["n10"], WithinAbs(274.0, 1e-12));
    CHECK_THAT(run["n11"], WithinAbs(305.0, 1e-12));
    CHECK_THAT(run["n12"], WithinAbs(335.0, 1e-12));
}

TEST_CASE("CS-GAP-NDAY: the name is case-insensitive", "[library-gaps][nday]") {
    auto run = runModel(R"(
a = NDAY_(3, 1)
b = nday_(3, 1)
c = Nday_(3, 1)
d = nDay_(3, 1)
)");
    REQUIRE(run.solveOk);
    for (const char* name : {"a", "b", "c", "d"}) CHECK_THAT(run[name], WithinAbs(60.0, 1e-12));
    CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
}

TEST_CASE("CS-GAP-NDAY: the day is used as given, the month is rounded to the nearest integer", "[library-gaps][nday]") {
    // SETP.LIB: monthx = round(month), so any month value is rounded, whatever the distance to an integer
    auto run = runModel(R"(
a = nDay_(12, 4.5)
m = 36 / 3 + 1e-12
b = nDay_(m, 4)
c = nDay_(12.0001, 4)
d = nDay_(2.5, 1)
e = nDay_(2.4, 1)
f = nDay_(11.6, 1)
g = nDay_(2, 28.0005)
h = nDay_(7, 31.0005)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(338.5, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(338.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(338.0, 1e-12));
    CHECK_THAT(run["d"], WithinAbs(60.0, 1e-12));      // 2.5 -> month 3
    CHECK_THAT(run["e"], WithinAbs(32.0, 1e-12));      // 2.4 -> month 2
    CHECK_THAT(run["f"], WithinAbs(335.0, 1e-12));     // 11.6 -> month 12
    // the day may exceed the length of its month by up to 0.001
    CHECK_THAT(run["g"], WithinAbs(59.0005, 1e-12));
    CHECK_THAT(run["h"], WithinAbs(212.0005, 1e-12));
}

TEST_CASE("CS-GAP-NDAY: automatic differentiation, 1 with respect to the day and 0 to the month", "[library-gaps][nday]") {
    // two independent variables: x0 = month, x1 = day
    std::vector<ADValue> args = {ADValue::independent(3.0, 0, 2), ADValue::independent(4.0, 1, 2)};
    const ADValue r = evaluateStandardFunction("nDay_", args);
    CHECK_THAT(r.value, WithinAbs(63.0, 1e-12));
    REQUIRE(r.gradient.size() == 2);
    CHECK(r.gradient[0] == 0.0);
    CHECK(r.gradient[1] == 1.0);

    // chain rule through the day argument: d nDay_(3, 2 x) / dx = 2
    const ADValue x = ADValue::independent(10.5, 1, 2);
    std::vector<ADValue> chained = {ADValue::constant(3.0, 2), x * ADValue::constant(2.0, 2)};
    const ADValue c = evaluateStandardFunction("NDAY_", chained);
    CHECK_THAT(c.value, WithinAbs(80.0, 1e-12));
    CHECK(c.gradient[0] == 0.0);
    CHECK_THAT(c.gradient[1], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-GAP-NDAY: the solver goes through nDay_ (Jacobian by automatic differentiation)", "[library-gaps][nday]") {
    SECTION("the day is the unknown") {
        auto run = runModel("nDay_(12, d) = 340\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["d"], WithinAbs(6.0, 1e-8));
    }
    SECTION("the day is an expression of the unknown") {
        auto run = runModel("nDay_(3, 2*x) = 80\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinAbs(10.5, 1e-8));
    }
    SECTION("coupled with another equation") {
        auto run = runModel("nDay_(7, d) + y = 200\ny = 2 * d\n");   // 181 + d + 2 d = 200
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["d"], WithinAbs(19.0 / 3.0, 1e-8));
        CHECK_THAT(run["y"], WithinAbs(38.0 / 3.0, 1e-8));
    }
}

TEST_CASE("CS-GAP-NDAY: usable in a FUNCTION and in the argument of other functions", "[library-gaps][nday]") {
    auto run = runModel(R"(
FUNCTION decl(n)
  decl = 23.45 * sin(360 * (284 + n) / 365)
END
n = nDay_(6, 21)
delta = decl(nDay_(6, 21))
m = max(nDay_(1, 1), nDay_(2, 1))
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["n"], WithinAbs(172.0, 1e-12));
    CHECK_THAT(run["delta"], WithinRel(23.45 * std::sin(360.0 * (284.0 + 172.0) / 365.0 * M_PI / 180.0), 1e-12));
    CHECK_THAT(run["m"], WithinAbs(32.0, 1e-12));
}

TEST_CASE("CS-GAP-NDAY: a month that does not round to 1..12 is an error, not a value", "[library-gaps][nday]") {
    for (const std::string month : {"13", "0", "-1", "12.5", "0.4", "13.4"}) {
        INFO("month = " << month);
        auto run = runModel("m = " + month + "\nn = nDay_(m, 1)\n");
        REQUIRE_FALSE(run.solveOk);
        CHECK(run.message.find("nDay_") != std::string::npos);
        CHECK(run.message.find("month must be between 1 and 12 after rounding") != std::string::npos);
        CHECK(run.message.find("got " + month) != std::string::npos);
        CHECK(run.vars.count("n") == 0);
    }
    // literal arguments, no variable involved
    auto literal = runModel("n = nDay_(13, 4)\n");
    REQUIRE_FALSE(literal.solveOk);
    CHECK(literal.message.find("nDay_: month must be between 1 and 12 after rounding, got 13") != std::string::npos);
}

TEST_CASE("CS-GAP-NDAY: a day below 1 or beyond the length of its month is an error, not a value", "[library-gaps][nday]") {
    SECTION("day below 1") {
        for (const std::string day : {"0", "-3", "0.5"}) {
            INFO("day = " << day);
            auto run = runModel("d = " + day + "\nn = nDay_(7, d)\n");
            REQUIRE_FALSE(run.solveOk);
            CHECK(run.message.find("nDay_: day must be at least 1, got " + day) != std::string::npos);
        }
    }
    SECTION("day beyond the length of the month, in a non-leap year") {
        const int lengths[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        for (int month = 1; month <= 12; ++month) {
            const int len = lengths[month - 1];
            INFO("month = " << month);
            // the last day of the month is accepted, the next one is not
            auto ok = runModel("n = nDay_(" + std::to_string(month) + ", " + std::to_string(len) + ")\n");
            CHECK(ok.solveOk);
            auto bad = runModel("n = nDay_(" + std::to_string(month) + ", " + std::to_string(len + 1) + ")\n");
            REQUIRE_FALSE(bad.solveOk);
            CHECK(bad.message.find("nDay_: day " + std::to_string(len + 1) + " exceeds the " + std::to_string(len) +
                                   " days of month " + std::to_string(month)) != std::string::npos);
        }
        // the cases of the register: no 29 February, no 31 April
        CHECK_FALSE(runModel("n = nDay_(2, 29)\n").solveOk);
        CHECK_FALSE(runModel("n = nDay_(4, 31)\n").solveOk);
        // a fraction of a day beyond the last day is an error from 0.001 on; the month is the rounded one
        auto frac = runModel("d = 28.002\nn = nDay_(2, d)\n");
        REQUIRE_FALSE(frac.solveOk);
        CHECK(frac.message.find("nDay_: day 28.002 exceeds the 28 days of month 2") != std::string::npos);
        auto rounded = runModel("m = 2.5\nn = nDay_(m, 32)\n");
        REQUIRE_FALSE(rounded.solveOk);
        CHECK(rounded.message.find("nDay_: day 32 exceeds the 31 days of month 3") != std::string::npos);
    }
}

TEST_CASE("CS-GAP-NDAY: the number of arguments is checked", "[library-gaps][nday]") {
    for (const std::string call : {"nDay_(3)", "nDay_(3, 1, 2)", "nDay_()"}) {
        INFO("call = " << call);
        auto run = runModel("n = " + call + "\n");
        REQUIRE_FALSE(run.solveOk);
        CHECK(run.message.find("nDay_ takes 2 arguments (month, day)") != std::string::npos);
    }
}

// ============================================================================
// CS-BUG-IMPLICIT-PROC, CS-BUG-LOCAL-UNASSIGNED, CS-GAP-STRING-LOCAL-UNSET: the statements of a
// FUNCTION/PROCEDURE body are assignments executed in order, as in EES (maintainer decision of
// 2026-10-12). An implicit statement in a body, and a read of a variable that the body has not
// assigned before (its locals are not initialised), are errors of the check of the bodies, reported
// with their line like the parse errors; they were skipped / read as 1, SUCCESS with a wrong value.
// A routine may still be called implicitly from the main program.
// ============================================================================

namespace {

using Messages = std::vector<std::string>;

const std::string kImplicitBodyText =
    "the statements of a FUNCTION or PROCEDURE are assignments executed in order, so the left-hand side "
    "must be a variable. Solve the equation in the main program, or iterate explicitly in the body (REPEAT ... UNTIL)";
const std::string kUnassignedLocalText =
    "before it is assigned: the variables of a FUNCTION or PROCEDURE other than its inputs are local and "
    "not initialised (pass main-program values as input arguments)";

/// "line: message" (parseErrorsWithLines) of an implicit statement in `routine` ("PROCEDURE 'pq'").
std::string implicitBodyError(int line, const std::string& routine) {
    return std::to_string(line) + ": implicit equation in " + routine + ": " + kImplicitBodyText;
}

/// "line: message" (parseErrorsWithLines) of a read of `var` before its assignment in `routine`.
std::string unassignedLocalError(int line, const std::string& var, const std::string& routine) {
    return std::to_string(line) + ": variable '" + var + "' is read in " + routine + " " + kUnassignedLocalText;
}

const std::string kColebrookImplicit = "1/sqrt(f) = -2*log10((r/1000)/3.7+2.51/(1E5*sqrt(f)))";

}  // namespace

TEST_CASE("CS-BUG-IMPLICIT-PROC: reproducer (a), an implicit statement after an assignment is an error, not skipped",
          "[library-gaps][implicit-proc]") {
    // was SUCCESS with fv = 1: the statement was skipped and f, never assigned, read as 1
    const std::string code = "PROCEDURE pq(r : f, g)\ng = r*2\n" + kColebrookImplicit + "\nend\ncall pq(1 : fv, gv)\n";
    CHECK(parseErrorsWithLines(code) == Messages{implicitBodyError(3, "PROCEDURE 'pq'")});

    const auto result = parseFull(code);
    CHECK_FALSE(result.success);
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].context == kColebrookImplicit);   // the source line, as for the other parse errors

    auto run = runModel(code);   // like any parse error: the model is not solved
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(run.message.find("Line 3: implicit equation in PROCEDURE 'pq'") != std::string::npos);
}

TEST_CASE("CS-BUG-IMPLICIT-PROC: reproducer (b), a body made of the implicit statement only", "[library-gaps][implicit-proc]") {
    const std::string code = "PROCEDURE pq(r : f)\n" + kColebrookImplicit + "\nend\ncall pq(1 : fv)\n";
    CHECK(parseErrorsWithLines(code) == Messages{implicitBodyError(2, "PROCEDURE 'pq'")});
    auto run = runModel(code);
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
}

TEST_CASE("CS-BUG-IMPLICIT-PROC: reproducer (c), an implicit statement in a FUNCTION", "[library-gaps][implicit-proc]") {
    const std::string code =
        "FUNCTION ff(r)\n1/sqrt(ff) = -2*log10((r/1000)/3.7+2.51/(1E5*sqrt(ff)))\nend\nfv = ff(1)\n";
    CHECK(parseErrorsWithLines(code) == Messages{implicitBodyError(2, "FUNCTION 'ff'")});
    auto run = runModel(code);
    CHECK_FALSE(run.parseOk);
    CHECK(run.message.find("Line 2: implicit equation in FUNCTION 'ff'") != std::string::npos);
}

TEST_CASE("CS-BUG-IMPLICIT-PROC: implicit statements in IF, REPEAT and DUPLICATE blocks, with := and semicolons",
          "[library-gaps][implicit-proc]") {
    const std::string fn = "FUNCTION 'f'";
    const std::string tail = "END\ny = f(1)\n";
    // block IF and single-line IF ... ELSE
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  IF (a > 0) THEN\n    2*f = a\n  ELSE\n    f = a\n  ENDIF\n" + tail) ==
          Messages{implicitBodyError(3, fn)});
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  IF (a > 0) THEN f = a ELSE f^2 = a\n" + tail) ==
          Messages{implicitBodyError(2, fn)});
    // REPEAT ... UNTIL and DUPLICATE (a variable of an implicit statement is not reported again later)
    CHECK(parseErrorsWithLines("PROCEDURE p(a : y)\n  k = 0\n  REPEAT\n    k + 1 = 2*k\n  UNTIL (k > 3)\n  y = k\nEND\n"
                               "CALL p(1 : z)\n") == Messages{implicitBodyError(4, "PROCEDURE 'p'")});
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  DUPLICATE i = 1, 3\n    x[i] + 1 = a*i\n  END\n  f = x[3]\n" + tail) ==
          Messages{implicitBodyError(3, fn)});
    // `:=` is an assignment as well, and a statement after a semicolon has the line of its line
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  f := a\n  2*f := a\n" + tail) == Messages{implicitBodyError(3, fn)});
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  s = a ; 2*f = s\n" + tail) == Messages{implicitBodyError(2, fn)});
    // every implicit statement is reported
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  2*s = a\n  2*f = s\n" + tail) ==
          Messages{implicitBodyError(2, fn), implicitBodyError(3, fn)});
}

TEST_CASE("CS-BUG-IMPLICIT-PROC: a lookup-table cell write in a body is unchanged (CS-BUG-LOOKUP-WRITE)",
          "[library-gaps][implicit-proc]") {
    // an open item of its own: no new error, the statement is still skipped when the body runs
    auto run = runModel("PROCEDURE w(M, h : y)\n  lookup('ev', 1, 1) = M*h\n  Lookup('ev', 2, 1) = M\n  y = M + h\nEND\n"
                        "CALL w(2, 3 : z)\n");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["z"], WithinAbs(5.0, 1e-12));
}

TEST_CASE("CS-BUG-IMPLICIT-PROC: an explicit routine may still be called implicitly from the main program",
          "[library-gaps][implicit-proc]") {
    SECTION("ff(x) = 10 solves x through the explicit FUNCTION ff") {
        auto run = runModel("FUNCTION ff(x)\n  ff = x^2 + 1\nEND\nff(x) = 10\n");
        REQUIRE(run.parseOk);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(3.0, 1e-9));
        CHECK(run.verified);
    }
    SECTION("a CALL inside an algebraic loop of the main program") {
        auto run = runModel("PROCEDURE sq(x : y)\n  y = x^2\nEND\nCALL sq(x : y)\ny + x = 6\n");
        REQUIRE(run.parseOk);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x"], WithinAbs(2.0, 1e-9));
        CHECK_THAT(run["y"], WithinAbs(4.0, 1e-9));
    }
}

TEST_CASE("CS-BUG-IMPLICIT-PROC: the explicit rewrite of reproducer (a), a fixed-point REPEAT ... UNTIL on Colebrook",
          "[library-gaps][implicit-proc]") {
    // Re = 1E5, RR = 0.001: f = 0.0221745359... (the value EES gives for the implicit form)
    auto run = runModel(R"(
PROCEDURE pq(r : f, g)
  g = r*2
  f = 0.02
  REPEAT
    f_old = f
    f = (-2*log10((r/1000)/3.7 + 2.51/(1E5*sqrt(f_old))))^(-2)
  UNTIL (abs(f - f_old) < 1E-12)
END
call pq(1 : fv, gv)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["fv"], WithinRel(0.022174535944515, 1e-9));
    CHECK_THAT(run["gv"], WithinAbs(2.0, 1e-12));
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: reproducer, a local read before its assignment is an error", "[library-gaps][local-unassigned]") {
    // was SUCCESS with yv = 2: b read as 1
    const std::string code = "PROCEDURE p(a : y)\ny = 2 * b\nb = a + 1\nEND\nCALL p(5 : yv)\n";
    CHECK(parseErrorsWithLines(code) == Messages{unassignedLocalError(2, "b", "PROCEDURE 'p'")});
    const auto result = parseFull(code);
    REQUIRE(result.errors.size() == 1);
    CHECK(result.errors[0].context == "y = 2 * b");
    auto run = runModel(code);
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
    CHECK(run.message.find("Line 2: variable 'b' is read in PROCEDURE 'p' before it is assigned") != std::string::npos);

    // the two statements in the order of execution: yv = 12
    auto ordered = runModel("PROCEDURE p(a : y)\nb = a + 1\ny = 2 * b\nEND\nCALL p(5 : yv)\n");
    REQUIRE(ordered.solveOk);
    CHECK_THAT(ordered["yv"], WithinAbs(12.0, 1e-12));
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: a name never assigned in the body, even when the main program defines it",
          "[library-gaps][local-unassigned]") {
    // a body does not see the main program: bb must be passed as an argument
    CHECK(parseErrorsWithLines("PROCEDURE p(a : y)\ny = 2 * bb\nEND\nbb = 10\nCALL p(5 : yv)\n") ==
          Messages{unassignedLocalError(2, "bb", "PROCEDURE 'p'")});
    // a misspelt name
    CHECK(parseErrorsWithLines("FUNCTION area(D)\n  r = D/2\n  area = 3.14159*rr^2\nEND\nA = area(1)\n") ==
          Messages{unassignedLocalError(3, "rr", "FUNCTION 'area'")});
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: in a FUNCTION", "[library-gaps][local-unassigned]") {
    // was SUCCESS with q = 2
    const std::string code = "FUNCTION f(a)\nf = 2*b\nb = a+1\nEND\nq = f(5)\n";
    CHECK(parseErrorsWithLines(code) == Messages{unassignedLocalError(2, "b", "FUNCTION 'f'")});
    auto run = runModel(code);
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: the outputs are not pre-assigned", "[library-gaps][local-unassigned]") {
    // the FUNCTION's own name read before it is assigned
    CHECK(parseErrorsWithLines("FUNCTION ff(x)\n  ff = ff + x\nEND\ny = ff(1)\n") ==
          Messages{unassignedLocalError(2, "ff", "FUNCTION 'ff'")});
    // a PROCEDURE output read before it is assigned
    CHECK(parseErrorsWithLines("PROCEDURE p(a : y, z)\n  z = y + a\n  y = a\nEND\nCALL p(1 : u, v)\n") ==
          Messages{unassignedLocalError(2, "y", "PROCEDURE 'p'")});
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: every output must be assigned by the body", "[library-gaps][local-unassigned]") {
    // was SUCCESS with the unassigned result read as 1
    const std::string never = "never assigns its output";
    auto fn = parseErrorsWithLines("FUNCTION ff(x)\n  y = 2*x\nEND\nq = ff(3)\n");
    REQUIRE(fn.size() == 1);
    CHECK(fn[0].rfind("1: FUNCTION 'ff' never assigns its output 'ff'", 0) == 0);
    auto proc = parseErrorsWithLines("x = 1\nPROCEDURE p(a : y, z)\n  y = a\nEND\nCALL p(1 : u, v)\n");
    REQUIRE(proc.size() == 1);
    CHECK(proc[0].rfind("2: PROCEDURE 'p' never assigns its output 'z'", 0) == 0);
    auto run = runModel("FUNCTION ff(x)\n  y = 2*x\nEND\nq = ff(3)\n");
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);

    // assigned in one branch only, or through a CALL of the body: accepted (flow-insensitive)
    CHECK(parseErrorsWithLines("FUNCTION sgn(x)\n  IF x > 0 THEN\n    sgn = 1\n  ELSE\n    sgn = -1\n  ENDIF\nEND\nq = sgn(3)\n").empty());
    CHECK(parseErrorsWithLines("PROCEDURE inner(a : b)\n  b = 2*a\nEND\nPROCEDURE outer(a : y)\n  CALL inner(a : y)\nEND\n"
                               "CALL outer(2 : yv)\n").empty());
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: every kind of read is checked", "[library-gaps][local-unassigned]") {
    const std::string fn = "FUNCTION 'f'";
    const std::string tail = "END\ny = f(1)\n";
    const std::string q = "PROCEDURE q(a : b)\n  b = a\nEND\n";   // lines 1-3
    // IF and UNTIL conditions
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  IF (c > 0) THEN f = 1 ELSE f = 2\n" + tail) ==
          Messages{unassignedLocalError(2, "c", fn)});
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  k = 0\n  REPEAT\n    k = k + 1\n  UNTIL (k > n)\n  f = k\n" + tail) ==
          Messages{unassignedLocalError(5, "n", fn)});
    // input arguments of a CALL, and the index of one of its outputs
    CHECK(parseErrorsWithLines(q + "FUNCTION f(a)\n  CALL q(w : v)\n  f = v\n" + tail) ==
          Messages{unassignedLocalError(5, "w", fn)});
    CHECK(parseErrorsWithLines(q + "FUNCTION f(a)\n  CALL q(a : x[k])\n  f = a\n" + tail) ==
          Messages{unassignedLocalError(5, "k", fn)});
    // array indices, on the right-hand side and on the left-hand side
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  x[1] = a\n  f = x[j]\n" + tail) ==
          Messages{unassignedLocalError(3, "j", fn)});
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  x[j] = a\n  f = a\n" + tail) == Messages{unassignedLocalError(2, "j", fn)});
    // arguments of a function call, values of the named arguments of a property call
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  f = max(a, m)\n" + tail) == Messages{unassignedLocalError(2, "m", fn)});
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  f = enthalpy(Water, T=t1, P=a)\n" + tail) ==
          Messages{unassignedLocalError(2, "t1", fn)});
    // the bounds of a DUPLICATE
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  DUPLICATE i = 1, n\n    x[i] = i\n  END\n  f = a\n" + tail) ==
          Messages{unassignedLocalError(2, "n", fn)});
    // the right-hand side read before the left-hand side is assigned: x = x + 1
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  x = x + 1\n  f = x\n" + tail) == Messages{unassignedLocalError(2, "x", fn)});
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: inside a REPEAT loop the text order is the rule", "[library-gaps][local-unassigned]") {
    // t is assigned later in the loop: its first read, in the first pass, finds no value
    const std::string loop = "  REPEAT\n    s = t + 1\n    t = k\n    k = k + 1\n  UNTIL (k > 3)\n  f = s\nEND\ny = f(1)\n";
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  k = 0\n" + loop) == Messages{unassignedLocalError(4, "t", "FUNCTION 'f'")});
    // assigned before the loop: k 0 -> 4, s = 1, 1, 2, 3
    auto run = runModel("FUNCTION f(a)\n  k = 0\n  t = 0\n" + loop);
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(3.0, 1e-12));
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: one error per variable, at its first read; a syntax error is not followed by others",
          "[library-gaps][local-unassigned]") {
    CHECK(parseErrorsWithLines("PROCEDURE p(a : y)\n  u = b + 1\n  v = 2*b\n  y = u + v + c\nEND\nCALL p(1 : z)\n") ==
          Messages{unassignedLocalError(2, "b", "PROCEDURE 'p'"), unassignedLocalError(4, "c", "PROCEDURE 'p'")});
    // a body that does not parse is not checked: the statement that failed would assign t
    CHECK(parseErrorsWithLines("FUNCTION f(a)\n  t = 3 4\n  f = t\nEND\ny = f(1)\n") == Messages{"2: Could not parse line"});
}

TEST_CASE("CS-BUG-LOCAL-UNASSIGNED: what a body may read keeps solving", "[library-gaps][local-unassigned]") {
    SECTION("a variable assigned in the branches of an IF, read after the IF") {
        auto run = runModel("FUNCTION f(a)\n  IF (a > 0) THEN\n    s = 10\n  ELSE\n    s = 20\n  ENDIF\n  f = s + a\nEND\n"
                            "y1 = f(1)\ny2 = f(-1)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y1"], WithinAbs(11.0, 1e-12));
        CHECK_THAT(run["y2"], WithinAbs(19.0, 1e-12));
        // the check is flow-insensitive: an assignment in one branch counts
        auto oneBranch = runModel("FUNCTION g(a)\n  IF (a > 0) THEN t = 5\n  g = t\nEND\ny = g(1)\n");
        REQUIRE(oneBranch.solveOk);
        CHECK_THAT(oneBranch["y"], WithinAbs(5.0, 1e-12));
    }
    SECTION("a REPEAT loop with a counter initialised before the loop") {
        auto run = runModel("FUNCTION nsteps(a)\n  k = 0\n  x = a\n  REPEAT\n    x = x/2\n    k = k + 1\n  UNTIL (x < 1)\n"
                            "  nsteps = k\nEND\nn = nsteps(10)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["n"], WithinAbs(4.0, 1e-12));
    }
    SECTION("a DUPLICATE loop filling an array, then summed") {
        auto run = runModel("FUNCTION sq(n)\n  DUPLICATE i = 1, 4\n    x[i] = i^2\n  END\n  s = 0\n  DUPLICATE i = 1, 4\n"
                            "    s = s + x[i]\n  END\n  sq = s + n\nEND\ntotal = sq(0)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["total"], WithinAbs(30.0, 1e-12));
    }
    SECTION("a CALL inside a body whose output is read afterwards") {
        auto run = runModel("PROCEDURE dbl(a : b)\n  b = 2*a\nEND\nFUNCTION f(a)\n  CALL dbl(a : d)\n  f = d + 1\nEND\ny = f(3)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(7.0, 1e-12));
    }
    SECTION("a FUNCTION reads its own name after assigning it") {
        auto run = runModel("FUNCTION ff(x)\n  ff = 1\n  ff = ff + x\nEND\ny = ff(4)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["y"], WithinAbs(5.0, 1e-12));
    }
    SECTION("pi and a # constant are no locals") {
        // their value in a body is CS-BUG-PI-FUNCTION; this test is about the check only
        auto run = runModel("FUNCTION fa(D)\n  fa = pi/4*D^2 + 0*g#\nEND\nx = fa(1)\n");
        REQUIRE(run.parseOk);
        CHECK(run.solveOk);
    }
    SECTION("an unquoted fluid name, named arguments and a units annotation in a property call") {
        auto run = runModel("FUNCTION hR(T, P)\n  hR = enthalpy(R134a, T=T, P=P) \"[J/kg]\"\nEND\nh1 = hR(300, 1E5)\n"
                            "h2 = enthalpy(R134a, T=300, P=1E5)\n");
        REQUIRE(run.parseOk);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["h1"], WithinRel(run["h2"], 1e-12));
    }
    SECTION("an array input written SP[1..5] in the header is an input (its form is CS-GAP-ARRAY-ARG)") {
        CHECK(parseErrorsWithLines("PROCEDURE p(SP[1..2] : y)\n  y = SP[1] + SP[2]\nEND\n").empty());
    }
    SECTION("string inputs, string locals assigned before they are read, string literals") {
        auto run = runModel("PROCEDURE p(f$, a : y, z)\n  IF f$ = 'water' THEN y = a ELSE y = 0\n  g$ = 'glycol'\n"
                            "  IF g$ = 'glycol' THEN z = 1 ELSE z = a\nEND\nCALL p('Water', 3 : yv, zv)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["yv"], WithinAbs(3.0, 1e-12));
        CHECK_THAT(run["zv"], WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("CS-GAP-STRING-LOCAL-UNSET: reproducer, a string local read before any assignment is an error",
          "[library-gaps][local-unassigned]") {
    // was the run-time error "String variable not found: fl$"; now the same check as the numeric locals
    const std::string body = "PROCEDURE p(a : y)\nIF fl$ = 'glycol' THEN\ny = 1\nELSE\ny = a\nENDIF\nEND\nCALL p(2 : yv)\n";
    const Messages expected{unassignedLocalError(2, "fl$", "PROCEDURE 'p'")};
    CHECK(parseErrorsWithLines(body) == expected);
    CHECK(parseErrorsWithLines(body + "fl$ = 'glycol'\n") == expected);   // a body does not see the main program
    auto run = runModel(body);
    CHECK_FALSE(run.parseOk);
    CHECK(run.message.find("Line 2: variable 'fl$' is read in PROCEDURE 'p' before it is assigned") != std::string::npos);

    // the CSL-0153 form: the string is an input argument of a CALL in the body
    CHECK(parseErrorsWithLines("PROCEDURE prop(f$, T : c)\n  c = T\nEND\nPROCEDURE hx(T : c)\n  CALL prop(fluidev$, T : c)\nEND\n"
                               "CALL hx(300 : cv)\n") == Messages{unassignedLocalError(5, "fluidev$", "PROCEDURE 'hx'")});
}

// ============================================================================
// CS-GAP-NDAY (continued): the other time functions of the EES library SETP.LIB
// (Duffie & Beckman): NumDay_(month), AveDay_(month), EqnTime_(n)
// ============================================================================

namespace {

constexpr int kMonthLength[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
constexpr int kAverageDay[12] = {17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344};

/// Equation of time in minutes, Duffie & Beckman eq. 1.5.3, written independently of the implementation.
double equationOfTimeMinutes(double n) {
    const double b = (n - 1.0) * 2.0 * M_PI / 365.0;
    return 229.2 * (0.000075 + 0.001868 * std::cos(b) - 0.032077 * std::sin(b) - 0.014615 * std::cos(2.0 * b) -
                    0.04089 * std::sin(2.0 * b));
}

/// d EqnTime_ / d n, by hand.
double equationOfTimeSlope(double n) {
    const double k = 2.0 * M_PI / 365.0, b = (n - 1.0) * k;
    return 229.2 * k * (-0.001868 * std::sin(b) - 0.032077 * std::cos(b) + 2.0 * 0.014615 * std::sin(2.0 * b) -
                        2.0 * 0.04089 * std::cos(2.0 * b));
}

/// EqnTime_ and its derivative with respect to n through the evaluator.
ADValue eqnTimeAD(double n) {
    return evaluateStandardFunction("EqnTime_", {ADValue::independent(n, 0, 1)});
}

}  // namespace

TEST_CASE("CS-GAP-NDAY: NumDay_, the number of days of each month (non-leap year)", "[library-gaps][nday]") {
    std::string code;
    for (int m = 1; m <= 12; ++m) code += "nd" + std::to_string(m) + " = NumDay_(" + std::to_string(m) + ")\n";
    // the sum over the year, written as one expression
    code += "total = NumDay_(1)";
    for (int m = 2; m <= 12; ++m) code += " + NumDay_(" + std::to_string(m) + ")";
    code += "\n";
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (int m = 1; m <= 12; ++m) {
        INFO("month = " << m);
        CHECK_THAT(run["nd" + std::to_string(m)], WithinAbs(kMonthLength[m - 1], 1e-12));
    }
    CHECK_THAT(run["total"], WithinAbs(365.0, 1e-12));
}

TEST_CASE("CS-GAP-NDAY: AveDay_, the day number of the average day of each month (Table 1.6.1)", "[library-gaps][nday]") {
    std::string code;
    for (int m = 1; m <= 12; ++m) code += "ad" + std::to_string(m) + " = AveDay_(" + std::to_string(m) + ")\n";
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (int m = 1; m <= 12; ++m) {
        INFO("month = " << m);
        CHECK_THAT(run["ad" + std::to_string(m)], WithinAbs(kAverageDay[m - 1], 1e-12));
    }
    // the same days as dates: 17 Jan, 16 Feb, 16 Mar, 15 Apr, 15 May, 11 Jun, 17 Jul, 16 Aug, 15 Sep, 15 Oct, 14 Nov, 10 Dec
    const int dayOfMonth[12] = {17, 16, 16, 15, 15, 11, 17, 16, 15, 15, 14, 10};
    std::string dates;
    for (int m = 1; m <= 12; ++m)
        dates += "a" + std::to_string(m) + " = AveDay_(" + std::to_string(m) + ") - nDay_(" + std::to_string(m) + ", " +
                 std::to_string(dayOfMonth[m - 1]) + ")\n";
    auto diff = runModel(dates);
    REQUIRE(diff.solveOk);
    for (int m = 1; m <= 12; ++m) {
        INFO("month = " << m);
        CHECK(diff["a" + std::to_string(m)] == 0.0);
    }
}

TEST_CASE("CS-GAP-NDAY: NumDay_, AveDay_ and nDay_ agree with each other", "[library-gaps][nday]") {
    std::string code;
    for (int m = 1; m <= 12; ++m) {
        const std::string s = std::to_string(m);
        // the last day of the month is the cumulative number of days of the year up to it
        code += "last" + s + " = nDay_(" + s + ", NumDay_(" + s + "))\n";
        code += "ave" + s + " = AveDay_(" + s + ")\n";
        code += "first" + s + " = nDay_(" + s + ", 1)\n";
    }
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    double cumulative = 0.0;
    for (int m = 1; m <= 12; ++m) {
        const std::string s = std::to_string(m);
        cumulative += kMonthLength[m - 1];
        INFO("month = " << m);
        CHECK_THAT(run["last" + s], WithinAbs(cumulative, 1e-12));
        CHECK_THAT(run["first" + s], WithinAbs(cumulative - kMonthLength[m - 1] + 1.0, 1e-12));
        CHECK(run["ave" + s] >= run["first" + s]);   // the average day is a day of its month
        CHECK(run["ave" + s] <= run["last" + s]);
    }
    CHECK_THAT(run["last12"], WithinAbs(365.0, 1e-12));
}

TEST_CASE("CS-GAP-NDAY: NumDay_ and AveDay_ round the month, an error when the rounded month is not 1..12", "[library-gaps][nday]") {
    auto run = runModel(R"(
a = NumDay_(2.4)
b = NumDay_(2.5)
c = NumDay_(12.0001)
d = AveDay_(11.6)
e = AveDay_(0.9)
m = 36 / 3 + 1e-12
f = AveDay_(m)
)");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(28.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(31.0, 1e-12));     // 2.5 -> month 3
    CHECK_THAT(run["c"], WithinAbs(31.0, 1e-12));
    CHECK_THAT(run["d"], WithinAbs(344.0, 1e-12));    // 11.6 -> month 12
    CHECK_THAT(run["e"], WithinAbs(17.0, 1e-12));     // 0.9 -> month 1
    CHECK_THAT(run["f"], WithinAbs(344.0, 1e-12));

    for (const std::string fn : {"NumDay_", "AveDay_"}) {
        for (const std::string month : {"13", "0", "-1", "12.5", "0.4", "13.4"}) {
            INFO(fn << " month = " << month);
            auto bad = runModel("m = " + month + "\nn = " + fn + "(m)\n");
            REQUIRE_FALSE(bad.solveOk);
            CHECK(bad.message.find(fn + ": month must be between 1 and 12 after rounding, got " + month) != std::string::npos);
            CHECK(bad.vars.count("n") == 0);
        }
    }
    auto literal = runModel("n = NumDay_(13)\n");
    REQUIRE_FALSE(literal.solveOk);
    CHECK(literal.message.find("NumDay_: month must be between 1 and 12 after rounding, got 13") != std::string::npos);
}

TEST_CASE("CS-GAP-NDAY: EqnTime_, the equation of time in minutes (eq. 1.5.3)", "[library-gaps][nday]") {
    // values of the formula evaluated in Python, B = (n-1)*360/365 degrees
    auto run = runModel(R"(
e1 = EqnTime_(1)
e32 = EqnTime_(32)
e40 = EqnTime_(40)
e172 = EqnTime_(172)
e338 = EqnTime_(338)
e365 = EqnTime_(365)
e0 = EqnTime_(0)
e1005 = EqnTime_(100.5)
n = nDay_(12, 4)
edec4 = EqnTime_(n)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["e1"], WithinAbs(-2.9044223999999996, 1e-10));
    CHECK_THAT(run["e32"], WithinAbs(-13.17909744608791, 1e-10));
    CHECK_THAT(run["e40"], WithinAbs(-14.10758891845034, 1e-10));
    CHECK_THAT(run["e172"], WithinAbs(-1.3247285561653666, 1e-10));
    CHECK_THAT(run["e338"], WithinAbs(9.5932959833589191, 1e-10));
    CHECK_THAT(run["edec4"], WithinAbs(9.5932959833589191, 1e-10));
    CHECK_THAT(run["e365"], WithinAbs(-2.4533484116245572, 1e-10));
    CHECK_THAT(run["e1005"], WithinAbs(-1.4795807393866627, 1e-10));
    // no range check (the library has none): the formula is periodic in n with period 365
    CHECK_THAT(run["e0"], WithinAbs(run["e365"], 1e-10));
    for (double n : {1.0, 32.0, 100.5, 338.0}) CHECK_THAT(eqnTimeAD(n).value, WithinAbs(equationOfTimeMinutes(n), 1e-12));

    // the known extremes of the equation of time: about -14.3 min in mid-February, +16.4 min in early November
    double lo = 1e9, hi = -1e9;
    int nLo = 0, nHi = 0;
    for (int d = 1; d <= 365; ++d) {
        const double e = eqnTimeAD(d).value;
        if (e < lo) { lo = e; nLo = d; }
        if (e > hi) { hi = e; nHi = d; }
    }
    CHECK_THAT(lo, WithinAbs(-14.27, 0.01));
    CHECK(nLo >= 40);
    CHECK(nLo <= 50);
    CHECK_THAT(hi, WithinAbs(16.39, 0.01));
    CHECK(nHi >= 300);
    CHECK(nHi <= 312);
}

TEST_CASE("CS-GAP-NDAY: EqnTime_ derivative with respect to n is exact", "[library-gaps][nday]") {
    // by hand, at the points of the value test
    for (double n : {1.0, 32.0, 40.0, 172.0, 338.0, 365.0}) {
        INFO("n = " << n);
        const ADValue r = eqnTimeAD(n);
        REQUIRE(r.gradient.size() == 1);
        CHECK_THAT(r.gradient[0], WithinAbs(equationOfTimeSlope(n), 1e-12));
    }
    // against a central finite difference over the whole year
    const double h = 1e-5;
    for (int d = 1; d <= 365; d += 4) {
        INFO("n = " << d);
        const double fd = (eqnTimeAD(d + h).value - eqnTimeAD(d - h).value) / (2.0 * h);
        CHECK_THAT(eqnTimeAD(d).gradient[0], WithinAbs(fd, 1e-8));
    }
    // two independent variables: x0 = something else, x1 = n; and the chain rule through the argument
    const ADValue second = evaluateStandardFunction("eqntime_", {ADValue::independent(40.0, 1, 2)});
    REQUIRE(second.gradient.size() == 2);
    CHECK(second.gradient[0] == 0.0);
    CHECK_THAT(second.gradient[1], WithinAbs(-0.064297463611859645, 1e-12));
    const ADValue x = ADValue::independent(20.0, 0, 1);
    const ADValue chained = evaluateStandardFunction("EQNTIME_", {x * ADValue::constant(2.0, 1)});
    CHECK_THAT(chained.value, WithinAbs(-14.10758891845034, 1e-10));
    CHECK_THAT(chained.gradient[0], WithinAbs(-0.12859492722371929, 1e-12));
}

TEST_CASE("CS-GAP-NDAY: the solver goes through EqnTime_ (Jacobian by automatic differentiation)", "[library-gaps][nday]") {
    SECTION("the day number is the unknown, from a guess a few days away") {
        auto run = runModel("EqnTime_(n) = 9.5932959833589191\n", {}, "n = 330\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["n"], WithinAbs(338.0, 1e-6));
    }
    SECTION("near the minimum of the equation of time, where the slope is small") {
        auto run = runModel("EqnTime_(n) = -14.10758891845034\n", {}, "n = 35\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["n"], WithinAbs(40.0, 1e-5));
    }
    SECTION("the day of the month is the unknown, through nDay_") {
        auto run = runModel("EqnTime_(nDay_(12, d)) = 9.5932959833589191\n", {}, "d = 3\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["d"], WithinAbs(4.0, 1e-6));
    }
    SECTION("coupled with another equation") {
        auto run = runModel("t = EqnTime_(n) + 2 * n\nn = 338 + y\ny = 0.5\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["t"], WithinAbs(equationOfTimeMinutes(338.5) + 677.0, 1e-9));
    }
}

TEST_CASE("CS-GAP-NDAY: the time functions are case-insensitive", "[library-gaps][nday]") {
    auto run = runModel(R"(
a1 = NUMDAY_(2)
a2 = numday_(2)
a3 = Numday_(2)
b1 = AVEDAY_(2)
b2 = aveday_(2)
b3 = Aveday_(2)
c1 = EQNTIME_(32)
c2 = eqntime_(32)
c3 = Eqntime_(32)
c4 = EqnTime_(32)
)");
    REQUIRE(run.solveOk);
    for (const char* name : {"a1", "a2", "a3"}) CHECK_THAT(run[name], WithinAbs(28.0, 1e-12));
    for (const char* name : {"b1", "b2", "b3"}) CHECK_THAT(run[name], WithinAbs(47.0, 1e-12));
    for (const char* name : {"c1", "c2", "c3", "c4"}) CHECK_THAT(run[name], WithinAbs(-13.17909744608791, 1e-10));
    CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
}

TEST_CASE("CS-GAP-NDAY: the time functions are usable in a FUNCTION", "[library-gaps][nday]") {
    auto run = runModel(R"(
FUNCTION solarterm(m, d)
  n = nDay_(m, d)
  solarterm = EqnTime_(n) + NumDay_(m) + AveDay_(m)
END
a = solarterm(12, 4)
b = solarterm(2, 28)
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinAbs(9.5932959833589191 + 31.0 + 344.0, 1e-9));
    CHECK_THAT(run["b"], WithinAbs(equationOfTimeMinutes(59.0) + 28.0 + 47.0, 1e-9));
}

TEST_CASE("CS-GAP-NDAY: the number of arguments of the time functions is checked", "[library-gaps][nday]") {
    for (const std::string call : {"NumDay_()", "NumDay_(1, 2)", "AveDay_()", "AveDay_(1, 2)", "EqnTime_()", "EqnTime_(1, 2)"}) {
        INFO("call = " << call);
        const std::string fn = call.substr(0, call.find('('));
        const std::string got = (call.find(',') != std::string::npos) ? "2" : "0";
        const std::string params = (fn == "EqnTime_") ? "(n)" : "(month)";
        auto run = runModel("n = " + call + "\n");
        REQUIRE_FALSE(run.solveOk);
        CHECK(run.message.find(fn + " takes 1 argument " + params + ", got " + got) != std::string::npos);
    }
}

// ============================================================================
// CS-BUG-NO-SOL-DELETE: `--no-sol` deleted an existing <model>.sol (the clean-up before the solve ran
// whatever the option), so a check run removed the regression baselines of the model library.
// These tests run the CLI binary built next to coolsolve_tests (the tests run from the build folder).
// ============================================================================

namespace {

/// Run `coolsolve <args> m.eescode` in `dir`; true when the command succeeds.
bool runCli(const fs::path& binary, const fs::path& dir, const std::string& args) {
    const std::string cmd = "cd '" + dir.string() + "' && '" + binary.string() + "' " + args + " m.eescode > out.txt 2>&1";
    return std::system(cmd.c_str()) == 0;
}

std::string fileText(const fs::path& path) {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

}  // namespace

TEST_CASE("CS-BUG-NO-SOL-DELETE: --no-sol leaves an existing .sol untouched and writes none", "[library-gaps][no-sol]") {
    const fs::path binary = fs::absolute("coolsolve");
    if (!fs::exists(binary)) SKIP("coolsolve binary not found next to the tests: " << binary.string());
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_no_sol";
    fs::remove_all(dir);
    fs::create_directories(dir);
    { std::ofstream f(dir / "m.eescode"); f << "a = 1\nb = a + 1\n"; }

    // a normal run writes the .sol
    REQUIRE(runCli(binary, dir, ""));
    REQUIRE(fs::exists(dir / "m.sol"));
    const std::string baseline = fileText(dir / "m.sol");
    CHECK(baseline.find("b = 2.") != std::string::npos);

    // --no-sol neither deletes nor rewrites it, even when the model changed (was: m.sol removed)
    { std::ofstream f(dir / "m.eescode"); f << "a = 5\nb = a + 1\n"; }
    REQUIRE(runCli(binary, dir, "--no-sol"));
    REQUIRE(fs::exists(dir / "m.sol"));
    CHECK(fileText(dir / "m.sol") == baseline);

    // --no-sol on a model without .sol does not create one
    fs::remove(dir / "m.sol");
    REQUIRE(runCli(binary, dir, "--no-sol"));
    CHECK_FALSE(fs::exists(dir / "m.sol"));

    // a normal run still replaces the solution of the previous run
    REQUIRE(runCli(binary, dir, ""));
    CHECK(fileText(dir / "m.sol").find("b = 6.") != std::string::npos);
    fs::remove_all(dir);
}

// ============================================================================
// CS-GAP-COLEBROOK: the built-in procedure Colebrook(Re, RR : f) of the EES library, the Darcy friction
// factor of the implicit Colebrook-White equation 1/sqrt(f) = -2*log10(RR/3.7 + 2.51/(Re*sqrt(f))),
// called with or without the colon (`Call Colebrook(Re_D, RR, f)` in the LABORELEC pipe files, CSL-0158).
// The reference values are those of the equation solved to 40 digits (mpmath).
// ============================================================================

namespace {

/// f of `CALL Colebrook(Re, RR : f)` run by the evaluator, with Re and RR the independent variables
/// 0 and 1: the gradient of f holds df/dRe and df/dRR.
ADValue colebrookCall(double re, double rr) {
    EESParser parser;
    const auto parsed = parser.parse("CALL Colebrook(Re, RR : f)\n");
    REQUIRE(parsed.errors.empty());
    REQUIRE(parsed.program.statements.size() == 1);
    REQUIRE(parsed.program.statements[0]->is<ProcedureCall>());
    ExpressionEvaluator eval(2);
    eval.setVariable("Re", ADValue::independent(re, 0, 2));
    eval.setVariable("RR", ADValue::independent(rr, 1, 2));
    eval.evaluateProcedureCall(parsed.program.statements[0]->as<ProcedureCall>());
    return eval.getVariable("f");
}

/// The message of the error thrown by colebrookCall; empty when it throws none.
std::string colebrookCallError(double re, double rr) {
    try {
        colebrookCall(re, rr);
    } catch (const std::exception& e) {
        return e.what();
    }
    return "";
}

/// Residual of the Colebrook-White equation for f.
double colebrookResidual(double f, double re, double rr) {
    return 1.0 / std::sqrt(f) + 2.0 * std::log10(rr / 3.7 + 2.51 / (re * std::sqrt(f)));
}

}  // namespace

TEST_CASE("CS-GAP-COLEBROOK: reproducer of the register, with and without the colon", "[library-gaps][colebrook]") {
    const double reference = 0.018513866077471643;   // Re = 1E5, RR = 1E-4
    for (const std::string call : {"Call Colebrook(Re, RR, f)", "CALL Colebrook(Re, RR : f)"}) {
        INFO(call);
        auto run = runModel("Re = 1E5\nRR = 0.0001\n" + call + "\nz = f\n");
        REQUIRE(run.parseOk);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["f"], WithinRel(reference, 1e-12));
        CHECK_THAT(run["z"], WithinRel(reference, 1e-12));
        CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
    }
}

TEST_CASE("CS-GAP-COLEBROOK: the name is case-insensitive and the arguments are expressions", "[library-gaps][colebrook]") {
    auto run = runModel(R"(
D = 0.05
eps = 0.000005
Re = 1E5
CALL COLEBROOK(Re, eps/D : f1)
call colebrook(2*Re/2, 1E-4, f2)
Call CoLeBrOoK(Re, (eps/D) : f3) "comment"
)");
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (const char* name : {"f1", "f2", "f3"}) CHECK_THAT(run[name], WithinRel(0.018513866077471643, 1e-12));
}

TEST_CASE("CS-GAP-COLEBROOK: reference values over the turbulent range", "[library-gaps][colebrook]") {
    struct Point { double re, rr, f; };
    const Point points[] = {
        {4e3, 0.0, 0.039907014055634898},   {1e5, 1e-4, 0.018513866077471643},
        {1e6, 1e-3, 0.019943465840476866},  {1e8, 0.0, 0.0059404663516367614},
        {1e8, 0.05, 0.071550904091083255},  {4e3, 0.05, 0.076986834889224867},
        {2.5e4, 2e-5, 0.024570135667965701}, {1e5, 1e-3, 0.022174535944515075},   // the last one: as EES
    };
    for (const auto& p : points) {
        INFO("Re = " << p.re << ", RR = " << p.rr);
        CHECK_THAT(colebrookCall(p.re, p.rr).value, WithinRel(p.f, 1e-12));
    }
}

TEST_CASE("CS-GAP-COLEBROOK: f satisfies the Colebrook-White equation to 1e-12 for Re 4E3...1E8, RR 0...0.05",
          "[library-gaps][colebrook]") {
    for (double re : {4e3, 6e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8}) {
        for (double rr : {0.0, 1e-7, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 0.03, 0.05}) {
            INFO("Re = " << re << ", RR = " << rr);
            const double f = colebrookCall(re, rr).value;
            CHECK(f > 0.005);
            CHECK(f < 0.08);
            CHECK(std::abs(colebrookResidual(f, re, rr)) < 1e-12);
        }
    }
}

TEST_CASE("CS-GAP-COLEBROOK: the CSL-0158 point, inside the PROCEDURE body as in the original (IF ... ELSE CALL)",
          "[library-gaps][colebrook]") {
    // Re_D is local to PIPE_GAS and not stored: 1164166.4513 is rebuilt from the stored values of the EES file
    // (V_dot, rho_i, Pr, Q, T_int: cp, then mu and lambda); RR = epsilon/di with di = 82.5 mm. EES: f = 0.01752150787
    // (the classical equation gives it to 5e-9). The body passes the check of the routine bodies (CS-BUG-IMPLICIT-PROC).
    const std::string code = R"(
PROCEDURE PIPE_GAS(Re_D, di_mm : f, DpDx)
  di = di_mm/1000
  epsilon=0.000046	"roughness of the pipe"
  RR=epsilon/di
  if (Re_D<=2000) then
    f = 64/Re_D
  else
   Call  Colebrook(Re_D, RR , f) "turbulent flow"
  endif
  DpDx = 0.5*f*7.776430299515753*24.05591817048639^2/di
END
CALL PIPE_GAS(1164166.4513, 82.5 : f, DpDx)
CALL PIPE_GAS(1500, 82.5 : f_lam, DpDx_lam)
)";
    CHECK(parseErrors(code).empty());
    auto run = runModel(code);
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["f"], WithinRel(0.01752150787, 1e-6));
    CHECK_THAT(run["DpDx"], WithinRel(0.5 * 0.01752150787 * 7.776430299515753 * 24.05591817048639 * 24.05591817048639 / 0.0825, 1e-6));
    CHECK_THAT(run["f_lam"], WithinRel(64.0 / 1500.0, 1e-9));   // the laminar branch, to the solver tolerance

    // the Re_D and RR of the CoolSolve variant (CoolProp viscosity): its f, from the implicit equation
    CHECK_THAT(colebrookCall(1.184661580325e+06, 5.575757575756e-04).value, WithinRel(1.751477912659e-02, 1e-11));
}

TEST_CASE("CS-GAP-COLEBROOK: usable in a FUNCTION, the output being the function itself or a local",
          "[library-gaps][colebrook]") {
    auto run = runModel(R"(
FUNCTION fric(Re, RR)
  IF (Re < 2300) THEN
    fric = 64/Re
  ELSE
    CALL Colebrook(Re, RR, fric)
  ENDIF
END
FUNCTION fric2(Re, RR)
  CALL Colebrook(Re, RR : ff)
  fric2 = 2*ff
END
a = fric(1E5, 1E-4)
b = fric(1000, 1E-4)
c = fric2(1E5, 1E-4)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a"], WithinRel(0.018513866077471643, 1e-12));
    CHECK_THAT(run["b"], WithinRel(0.064, 1e-14));
    CHECK_THAT(run["c"], WithinRel(2 * 0.018513866077471643, 1e-12));
}

TEST_CASE("CS-GAP-COLEBROOK: derivatives by the implicit function theorem (AD)", "[library-gaps][colebrook]") {
    // at Re = 1E5, RR = 1E-4 (derivatives of the 40-digit solution)
    const ADValue f = colebrookCall(1e5, 1e-4);
    REQUIRE(f.gradient.size() == 2);
    CHECK_THAT(f.gradient[0], WithinRel(-3.4602179397150064e-8, 1e-10));
    CHECK_THAT(f.gradient[1], WithinRel(5.0696335336773551, 1e-10));

    // against central finite differences over the range
    for (double re : {4e3, 1e5, 1e7}) {
        for (double rr : {1e-5, 1e-3, 0.05}) {
            INFO("Re = " << re << ", RR = " << rr);
            const ADValue g = colebrookCall(re, rr);
            const double hRe = re * 1e-4, hRR = rr * 1e-4;   // truncation error ~1e-8, round-off negligible
            const double dRe = (colebrookCall(re + hRe, rr).value - colebrookCall(re - hRe, rr).value) / (2 * hRe);
            const double dRR = (colebrookCall(re, rr + hRR).value - colebrookCall(re, rr - hRR).value) / (2 * hRR);
            CHECK_THAT(g.gradient[0], WithinRel(dRe, 1e-6));
            CHECK_THAT(g.gradient[1], WithinRel(dRR, 1e-6));
        }
    }
    // RR = 0 (smooth pipe): df/dRR is still defined
    const ADValue smooth = colebrookCall(1e5, 0.0);
    const double h = 1e-9;
    CHECK_THAT(smooth.gradient[1], WithinRel((colebrookCall(1e5, h).value - smooth.value) / h, 1e-4));
}

TEST_CASE("CS-GAP-COLEBROOK: the solver goes through Colebrook (Re unknown, f given)", "[library-gaps][colebrook]") {
    const double re = 67137.863981363882;   // f(Re, 1E-4) = 0.02
    SECTION("CALL in the main program") {
        auto run = runModel("CALL Colebrook(Re, 1E-4 : f)\nf = 0.02\n", {}, "Re = 50000\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["Re"], WithinRel(re, 1e-7));
    }
    SECTION("colon-less CALL coupled with other equations") {
        auto run = runModel("Re = rho*U*D/mu\nrho = 1000\nD = 0.05\nmu = 0.001\nCall Colebrook(Re, 1E-4, f)\nf = 0.02\n",
                            {}, "U = 1\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["U"], WithinRel(re * 0.001 / (1000 * 0.05), 1e-7));
    }
    SECTION("through a FUNCTION called with the unknown") {
        auto run = runModel("FUNCTION fric(Re, RR)\n  CALL Colebrook(Re, RR : fric)\nEND\nfric(Re, 1E-4) = 0.02\n",
                            {}, "Re = 50000\n");
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["Re"], WithinRel(re, 1e-7));
    }
}

TEST_CASE("CS-GAP-COLEBROOK: Re <= 0, RR < 0, NaN and RR >= 3.7 are errors naming Colebrook", "[library-gaps][colebrook]") {
    SECTION("in a model") {
        const std::pair<std::string, std::string> cases[] = {
            {"Re = 0\nRR = 1E-4\n", "Colebrook: the Reynolds number Re must be positive"},
            {"Re = -1E5\nRR = 1E-4\n", "Colebrook: the Reynolds number Re must be positive"},
            {"Re = 1E5\nRR = -1E-4\n", "Colebrook: the relative roughness RR"},
        };
        for (const auto& [inputs, message] : cases) {
            INFO(inputs);
            auto run = runModel(inputs + "CALL Colebrook(Re, RR : f)\n");
            REQUIRE(run.parseOk);
            REQUIRE_FALSE(run.solveOk);
            CHECK(run.message.find(message) != std::string::npos);
        }
    }
    SECTION("the messages") {
        const double nan = std::nan("");
        CHECK(colebrookCallError(0.0, 1e-4).find("Colebrook: the Reynolds number Re must be positive") != std::string::npos);
        CHECK(colebrookCallError(-1e5, 1e-4).find("Colebrook: the Reynolds number Re must be positive") != std::string::npos);
        CHECK(colebrookCallError(nan, 1e-4).find("Colebrook: the Reynolds number Re must be positive") != std::string::npos);
        CHECK(colebrookCallError(1e5, -1e-4).find("Colebrook: the relative roughness RR") != std::string::npos);
        CHECK(colebrookCallError(1e5, nan).find("Colebrook: the relative roughness RR") != std::string::npos);
        CHECK(colebrookCallError(1e5, 3.7).find("Colebrook: the relative roughness RR") != std::string::npos);
        CHECK(colebrookCallError(1e5, 0.0).empty());
    }
}

TEST_CASE("CS-GAP-COLEBROOK: the number of arguments and the output are checked", "[library-gaps][colebrook]") {
    const std::string arity = "CALL Colebrook takes 2 inputs (Re, RR) and 1 output (f)";
    for (const std::string call : {"CALL Colebrook(1E5 : f)", "CALL Colebrook(1E5, 1E-4, 2 : f)",
                                   "CALL Colebrook(1E5, 1E-4 : f, g)", "CALL Colebrook(1E5, 1E-4)",
                                   "CALL Colebrook(1E5, 1E-4, f, g)", "CALL Colebrook(1E5)"}) {
        INFO(call);
        const auto errors = parseErrors(call + "\n");
        REQUIRE(errors.size() == 1);
        CHECK(errors[0].find(arity) != std::string::npos);
    }
    // the same in a routine body
    CHECK(anyContains(parseErrors("PROCEDURE p(x : y)\n  CALL Colebrook(x : y)\nEND\nCALL p(1E5 : z)\n"), arity));

    // the third argument of the colon-less form is the output: a variable
    const auto constant = parseErrors("CALL Colebrook(1E5, 1E-4, 0.02)\n");
    REQUIRE(constant.size() == 1);
    CHECK(constant[0].find("the third argument of CALL Colebrook(Re, RR, f) is its output f and must be a variable") !=
          std::string::npos);
    CHECK(anyContains(parseErrors("CALL Colebrook(1E5, 1E-4, 2*f)\n"), "must be a variable"));
    // inside a routine, an expression output is refused as for any CALL
    CHECK(anyContains(parseErrors("PROCEDURE p(x : y)\n  CALL Colebrook(x, 1E-4 : 2*y)\nEND\nCALL p(1E5 : z)\n"),
                      "Output 1 of 'Colebrook' is the expression"));
    // an array element is a variable
    auto run = runModel("Re[1] = 1E5\nCALL Colebrook(Re[1], 1E-4, f[1])\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["f[1]"], WithinRel(0.018513866077471643, 1e-12));
}

TEST_CASE("CS-GAP-COLEBROOK: a PROCEDURE Colebrook of the model replaces the built-in", "[library-gaps][colebrook]") {
    SECTION("same signature") {
        auto run = runModel(R"(
PROCEDURE Colebrook(Re, RR : f)
  f = 0.316*Re^(-0.25)   "Blasius"
END
CALL Colebrook(1E5, 0 : f1)
CALL COLEBROOK(1E5, 0 : f2)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["f1"], WithinRel(0.316 * std::pow(1e5, -0.25), 1e-14));
        CHECK_THAT(run["f2"], WithinRel(0.316 * std::pow(1e5, -0.25), 1e-14));
    }
    SECTION("another signature: the built-in rules do not apply") {
        auto run = runModel(R"(
PROCEDURE colebrook(Re, eps, D : f, RR)
  RR = eps/D
  f = 64/Re + RR
END
CALL Colebrook(1000, 0.001, 0.1 : f, RR)
)");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["f"], WithinRel(0.074, 1e-14));
        CHECK_THAT(run["RR"], WithinRel(0.01, 1e-14));
    }
}

// ============================================================================
// CS-BUG-PI-FUNCTION: the built-in constants (pi, g#, R#, ...) inside FUNCTION/PROCEDURE bodies.
// The body evaluator started with no variables, so `pi` was read as an unset local (= 1) and every
// `pi/4*D^2` of a body was silently wrong.
// ============================================================================

TEST_CASE("CS-BUG-PI-FUNCTION: pi in the body of a FUNCTION", "[library-gaps][pi-function]") {
    const double pi = 3.14159265358979323846;
    // the reproducer of the register: the area of a circle
    auto run = runModel(R"(
FUNCTION fa(D)
  a = PI/4*D^2
  fa = a
END
x1 = fa(1)
y1 = PI/4
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.warnings.empty());
    CHECK_THAT(run["x1"], WithinRel(pi / 4, 1e-14));
    CHECK_THAT(run["y1"], WithinRel(pi / 4, 1e-14));

    // the three spellings of the constant, and the function form
    auto spellings = runModel(R"(
FUNCTION f(x)
  f = PI*x + 10*pi*x + 100*Pi*x + 1000*pi()*x
END
y = f(2)
)");
    REQUIRE(spellings.solveOk);
    CHECK_THAT(spellings["y"], WithinRel(1111 * pi * 2, 1e-14));
}

TEST_CASE("CS-BUG-PI-FUNCTION: pi in the body of a PROCEDURE, and in a function called by a body", "[library-gaps][pi-function]") {
    const double pi = 3.14159265358979323846;
    auto run = runModel(R"(
PROCEDURE circle(D : A, P)
  A = pi/4*D^2
  P = PI*D
END
CALL circle(2 : A1, P1)
FUNCTION inner(x)
  inner = Pi*x
END
FUNCTION outer(x)
  outer = inner(x) + pi
END
PROCEDURE proc2(x : y)
  y = outer(x)
END
z = outer(3)
CALL proc2(3 : w)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["A1"], WithinRel(pi, 1e-14));
    CHECK_THAT(run["P1"], WithinRel(2 * pi, 1e-14));
    CHECK_THAT(run["z"], WithinRel(4 * pi, 1e-14));
    CHECK_THAT(run["w"], WithinRel(4 * pi, 1e-14));
}

TEST_CASE("CS-BUG-PI-FUNCTION: every built-in constant, in every spelling, inside a body", "[library-gaps][pi-function]") {
    // the value is the one of the main program: src/constants.cpp is the single table. The results are
    // scaled by 1E30 so that the 1E-34 constants are not lost in the Newton step from the guess value 1
    // (x - (x - 1E-19) is 0 in double precision; unrelated to the constants).
    int checked = 0;
    for (const auto& [name, info] : Constants::getAll()) {
        INFO(name);
        auto run = runModel("FUNCTION fc(x)\n  fc = " + name + "*x*1E30\nEND\nPROCEDURE pc(x : y)\n  y = " + name +
                            "*x*1E30\nEND\nr = fc(1)\nCALL pc(1 : s)\nm = " + name + "*1E30\n");
        REQUIRE(run.parseOk);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["r"], WithinRel(info.value * 1e30, 1e-13));
        CHECK_THAT(run["s"], WithinRel(info.value * 1e30, 1e-13));
        CHECK_THAT(run["m"], WithinRel(info.value * 1e30, 1e-13));
        ++checked;
    }
    CHECK(checked == static_cast<int>(Constants::getAll().size()));
    CHECK(checked >= 17);   // PI, pi, Pi and the # constants
    // spot checks of the values themselves
    CHECK(Constants::getValue("g#") == 9.80665);
    auto g = runModel("FUNCTION w(m)\n  w = m*g#\nEND\nF = w(2)\n");
    REQUIRE(g.solveOk);
    CHECK_THAT(g["F"], WithinRel(2 * 9.80665, 1e-14));
}

TEST_CASE("CS-BUG-PI-FUNCTION: constants in the IF branches and the loops of a body", "[library-gaps][pi-function]") {
    const double pi = 3.14159265358979323846;
    auto run = runModel(R"(
FUNCTION f(x)
  IF (x > 0) THEN
    r = pi*x
  ELSE
    r = -pi*x
  ENDIF
  IF (pi > 3) THEN s = 1 ELSE s = 0
  IF (x > 100) THEN t = 0 ELSE t = g#
  k := 0
  REPEAT
    k := k + 1
    u := pi*k
  UNTIL (u > 10)
  DUPLICATE i = 1, 3
    a[i] = pi*i
  END
  f = r + 10*s + 100*t + 1000*k + 10000*u + 100000*(a[1] + a[2] + a[3])
END
y1 = f(2)
y2 = f(-2)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    const double common = 10 + 100 * 9.80665 + 1000 * 4 + 10000 * (4 * pi) + 100000 * 6 * pi;
    CHECK_THAT(run["y1"], WithinRel(2 * pi + common, 1e-14));
    CHECK_THAT(run["y2"], WithinRel(2 * pi + common, 1e-14));
}

TEST_CASE("CS-BUG-PI-FUNCTION: derivatives through a body that uses pi", "[library-gaps][pi-function]") {
    const double pi = 3.14159265358979323846;
    // an implicit main-program equation solved through a FUNCTION: the Jacobian comes from the body
    auto run = runModel(R"(
FUNCTION fa(D)
  a = PI/4*D^2
  fa = a
END
fa(D) = 1
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["D"], WithinRel(std::sqrt(4 / pi), 1e-9));

    // and through a PROCEDURE that is part of a coupled block, with g# in the body
    auto coupled = runModel(R"(
PROCEDURE fall(h : v)
  v = sqrt(2*g#*h)
END
CALL fall(H : v)
v = 5
Area = pi/4*Dia^2
Area = 2
)");
    REQUIRE(coupled.solveOk);
    CHECK_THAT(coupled["H"], WithinRel(25 / (2 * 9.80665), 1e-9));
    CHECK_THAT(coupled["Dia"], WithinRel(std::sqrt(8 / pi), 1e-9));
}

TEST_CASE("CS-BUG-PI-FUNCTION: a local variable of a body named like a constant keeps its own value", "[library-gaps][pi-function]") {
    // `c#`, `h#`, `k#`... are also plausible local names (the speed of light is a constant of the main program only
    // when nothing sets it): the local assignment wins in the body, the main program keeps the constant
    auto run = runModel(R"(
FUNCTION f(x)
  c# = 2*x
  f = c#
END
FUNCTION g(x)
  pi = 3
  g = pi*x
END
FUNCTION h(x)
  h = pi*x + c#*1E-8
END
y = f(4)
z = g(2)
w = h(1)
m = c#*1E-8
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    const double pi = 3.14159265358979323846;
    CHECK_THAT(run["y"], WithinRel(8.0, 1e-14));
    CHECK_THAT(run["z"], WithinRel(6.0, 1e-14));
    // the unset local of another body is the constant, not 1
    CHECK_THAT(run["w"], WithinRel(pi + 299792458.0 * 1e-8, 1e-14));
    CHECK_THAT(run["m"], WithinRel(299792458.0 * 1e-8, 1e-12));
}

// ============================================================================
// CS-GAP-ARRAY-ARG: an array range X[a..b] stands for its elements X[a], ..., X[b] (EES help, Array Range
// Notation): in the arguments of a function or of a CALL (outputs included), in the header of a FUNCTION or
// PROCEDURE and in the aggregate built-ins (SUM, AVERAGE, MAX, ...). Its bounds are integers known before
// the solve, resolved like the bounds of a DUPLICATE loop (CS-BUG-DUPLICATE-VAR-BOUND): numbers, $CONSTANTs,
// arithmetic of those, and variables of the main program set to a constant expression (`n = 3`).
// ============================================================================

namespace {

/// Number of equations of the IR built from `code`, which must parse without error.
int equationCount(const std::string& code) {
    const auto parsed = parseFull(code);
    REQUIRE(parsed.errors.empty());
    return IR::fromAST(parsed.program).getEquationCount();
}

}  // namespace

TEST_CASE("CS-GAP-ARRAY-ARG: reproducer of the register", "[library-gaps][array-arg]") {
    // was: "Line 7: Could not parse index 1 of 'SP': '1..3'" (EES: z = 6)
    auto run = runModel("SP[1]=1\nSP[2]=2\nSP[3]=3\nPROCEDURE pp(a, b, c : y)\ny = a+b+c\nend\ncall pp(SP[1..3] : yv)\n"
                        "z = yv\n");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["z"], WithinAbs(6.0, 1e-12));
}

TEST_CASE("CS-GAP-ARRAY-ARG: a range in the header of a PROCEDURE and of a FUNCTION", "[library-gaps][array-arg]") {
    // the header p(a, X[1..3], b) declares 5 inputs, X being the array of the body (constant subscripts, a
    // REPEAT index, a DUPLICATE index). Any 5 values may be passed: a matching range, single elements,
    // expressions, a shorter range completed by other values.
    auto run = runModel(R"(
PROCEDURE p(a, X[1..3], b : s, m)
  s = a*1000 + X[1]*100 + X[2]*10 + X[3] + b/10
  m := 0
  i := 0
  REPEAT
    i := i + 1
    m := m + X[i]^2
  UNTIL (i >= 3)
END
FUNCTION f(Y[1..2], c)
  t := 0
  DUPLICATE k = 1, 2
    t := t + k*Y[k]
  END
  f = t + c
END
V[1] = 4
V[2] = 5
V[3] = 6
CALL p(1, V[1..3], 2 : s1, m1)
CALL p(1, V[1], V[2], V[3], 2 : s2, m2)
CALL p(1, 7, 8, 2+7, 2 : s3, m3)
CALL p(1, V[2..3], 3, 2 : s4, m4)
y1 = f(V[2..3], 1)
y2 = f(5, 6, 1)
y3 = f(V[1..2], V[3])
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["s1"], WithinAbs(1456.2, 1e-9));
    CHECK_THAT(run["m1"], WithinAbs(16.0 + 25.0 + 36.0, 1e-12));
    CHECK_THAT(run["s2"], WithinAbs(1456.2, 1e-9));
    CHECK_THAT(run["m2"], WithinAbs(16.0 + 25.0 + 36.0, 1e-12));
    CHECK_THAT(run["s3"], WithinAbs(1789.2, 1e-9));
    CHECK_THAT(run["m3"], WithinAbs(49.0 + 64.0 + 81.0, 1e-12));
    CHECK_THAT(run["s4"], WithinAbs(1563.2, 1e-9));
    CHECK_THAT(run["y1"], WithinAbs(5.0 + 2 * 6.0 + 1.0, 1e-12));
    CHECK_THAT(run["y2"], WithinAbs(5.0 + 2 * 6.0 + 1.0, 1e-12));
    CHECK_THAT(run["y3"], WithinAbs(4.0 + 2 * 5.0 + 6.0, 1e-12));
}

TEST_CASE("CS-GAP-ARRAY-ARG: ranges as outputs of a CALL and of a PROCEDURE header", "[library-gaps][array-arg]") {
    auto run = runModel(R"(
PROCEDURE sq(X[1..3] : Y[1..3], tot)
  DUPLICATE k = 1, 3
    Y[k] = X[k]^2
  END
  tot = Y[1] + Y[2] + Y[3]
END
DUPLICATE i = 1, 3
  x[i] = i + 1
END
CALL sq(x[1..3] : y[1..3], t)
CALL sq(x[1..3] : w[1], w[2], w[3], u)
CALL sq(1, 2, 3 : z[2..4], v)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (int i = 1; i <= 3; ++i) {
        CHECK_THAT(run["y[" + std::to_string(i) + "]"], WithinAbs((i + 1.0) * (i + 1.0), 1e-12));
        CHECK_THAT(run["w[" + std::to_string(i) + "]"], WithinAbs((i + 1.0) * (i + 1.0), 1e-12));
        CHECK_THAT(run["z[" + std::to_string(i + 1) + "]"], WithinAbs(1.0 * i * i, 1e-12));
    }
    CHECK_THAT(run["t"], WithinAbs(29.0, 1e-12));
    CHECK_THAT(run["u"], WithinAbs(29.0, 1e-12));
    CHECK_THAT(run["v"], WithinAbs(14.0, 1e-12));
}

TEST_CASE("CS-GAP-ARRAY-ARG: bounds from a $CONSTANT and from variables set to a constant, in the aggregate built-ins",
          "[library-gaps][array-arg]") {
    const std::string code = R"(
$CONSTANT N#=4
n = 3
m = n - 1
x[1] = 1
x[2] = 2
x[3] = 4
x[4] = 8
s_c = SUM(x[1..N#])
s_v = SUM(x[1..n])
a = AVERAGE(x[1..N#])
mx = MAX(x[m..n])
mn = MIN(x[2..N#-1])
p = PRODUCT(x[1..3])
sd = STDDEV(x[1..2])
mixed = SUM(x[1..2], 10, x[N#])
)";
    CHECK(equationCount(code) == 14);
    auto run = runModel(code);
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["s_c"], WithinAbs(15.0, 1e-12));
    CHECK_THAT(run["s_v"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["a"], WithinAbs(3.75, 1e-12));
    CHECK_THAT(run["mx"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["mn"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["p"], WithinAbs(8.0, 1e-12));
    CHECK_THAT(run["sd"], WithinAbs(0.5, 1e-12));
    CHECK_THAT(run["mixed"], WithinAbs(21.0, 1e-12));
}

TEST_CASE("CS-GAP-ARRAY-ARG: $CONSTANT sizes as in orc_solar_complex, the constant defined below its uses",
          "[library-gaps][array-arg]") {
    // header ranges, a CALL with input and output ranges, a function call with a range inside a DUPLICATE,
    // AVERAGE of a range; the PROCEDURE has an input named like the constant (it is its local)
    auto run = runModel(R"(
FUNCTION q(T, i, EPS[1..imax#])
  q = EPS[i]*T
END
PROCEDURE props(imax#, T[1..imax#] : EPS[1..imax#], tot)
  j := 0
  tot := 0
  REPEAT
    j := j + 1
    EPS[j] := T[j]/100
    tot := tot + EPS[j]
  UNTIL (j >= imax#)
END
DUPLICATE i = 1, imax#
  T_3[i] = 10*i
END
CALL props(imax#, T_3[1..imax#] : eps_3[1..imax#], total)
DUPLICATE i = 1, imax#
  q_3[i] = q(T_3[i], i, eps_3[1..imax#])
END
T_bar = AVERAGE(T_3[1..imax#])
$CONSTANT imax#=5
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    for (int i = 1; i <= 5; ++i) {
        CHECK_THAT(run["eps_3[" + std::to_string(i) + "]"], WithinAbs(0.1 * i, 1e-12));
        CHECK_THAT(run["q_3[" + std::to_string(i) + "]"], WithinAbs(1.0 * i * i, 1e-12));
    }
    CHECK_THAT(run["total"], WithinAbs(1.5, 1e-12));
    CHECK_THAT(run["T_bar"], WithinAbs(30.0, 1e-12));
}

TEST_CASE("CS-GAP-ARRAY-ARG: ranges in the CALLs and the function calls of a body", "[library-gaps][array-arg]") {
    auto run = runModel(R"(
$CONSTANT K#=3
PROCEDURE inner(a, b, c : s)
  s = a + 10*b + 100*c
END
FUNCTION total(Z[1..K#])
  total = SUM(Z[1..K#]) + MAX(Z[2..3])
END
PROCEDURE outer(X[1..3] : r, t, u)
  CALL inner(X[1..3] : r)
  CALL inner(X[2..K#], 0 : t)
  u = total(X[1..K#])
END
CALL outer(1, 2, 3 : r, t, u)
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["r"], WithinAbs(321.0, 1e-12));
    CHECK_THAT(run["t"], WithinAbs(32.0, 1e-12));
    CHECK_THAT(run["u"], WithinAbs(9.0, 1e-12));
}

TEST_CASE("CS-GAP-ARRAY-ARG: a range on one index of a 2-D array, and on both", "[library-gaps][array-arg]") {
    auto run = runModel(R"(
FUNCTION corner(M[1..2,1..2])
  corner = M[1,2] - M[2,1]
END
A[1,1] = 1
A[1,2] = 2
A[2,1] = 3
A[2,2] = 4
s_row = SUM(A[2,1..2])
s_col = SUM(A[1..2,1])
s_all = SUM(A[1..2,1..2])
c = corner(A[1..2,1..2])
)");
    REQUIRE(run.parseOk);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["s_row"], WithinAbs(7.0, 1e-12));
    CHECK_THAT(run["s_col"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["s_all"], WithinAbs(10.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(-1.0, 1e-12));
}

TEST_CASE("CS-GAP-ARRAY-ARG: a number of values that does not match the header is an error", "[library-gaps][array-arg]") {
    // the register reproducer with a range one element short: counted after the ranges are expanded
    auto errors = parseErrorsWithLines("SP[1]=1\nSP[2]=2\nPROCEDURE pp(a, b, c : y)\ny = a+b+c\nend\ncall pp(SP[1..2] : yv)\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].rfind("6: ", 0) == 0);
    CHECK(anyContains(errors, "CALL pp gives 2 inputs"));
    CHECK(anyContains(errors, "PROCEDURE pp (line 3) takes 3"));
    // a header range: f(X[1..3]) takes 3 arguments
    errors = parseErrorsWithLines("FUNCTION f(X[1..3])\n  f = X[1]\nEND\ny = f(1, 2)\nz = f(1, 2, 3, 4)\n");
    REQUIRE(errors.size() == 2);
    CHECK(errors[0].rfind("4: ", 0) == 0);
    CHECK(anyContains(errors, "f(...) gives 2 arguments"));
    CHECK(anyContains(errors, "f(...) gives 4 arguments"));
    CHECK(anyContains(errors, "FUNCTION f (line 1) takes 3"));
    // outputs
    errors = parseErrorsWithLines("PROCEDURE q(a : Y[1..2])\n  Y[1] = a\n  Y[2] = a\nEND\nCALL q(1 : z[1..3])\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].rfind("5: ", 0) == 0);
    CHECK(anyContains(errors, "CALL q gives 3 outputs"));
    CHECK(anyContains(errors, "PROCEDURE q (line 1) returns 2"));
    // the same in a body
    errors = parseErrorsWithLines("PROCEDURE pp(a, b : y)\n  y = a+b\nEND\nPROCEDURE o(X[1..3] : r)\n  CALL pp(X[1..3] : r)\nEND\n");
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].rfind("5: ", 0) == 0);
    CHECK(anyContains(errors, "CALL pp gives 3 inputs"));
    // the solve is not attempted
    auto run = runModel("SP[1]=1\nSP[2]=2\nPROCEDURE pp(a, b, c : y)\ny = a+b+c\nend\ncall pp(SP[1..2] : yv)\n");
    CHECK_FALSE(run.parseOk);
    CHECK_FALSE(run.solveOk);
}

TEST_CASE("CS-GAP-ARRAY-ARG: a reversed, non-integer, misplaced or unknown range is an error", "[library-gaps][array-arg]") {
    const std::string x = "x[1]=1\nx[2]=2\nx[3]=3\n";
    SECTION("reversed and non-integer bounds") {
        CHECK(parseErrorsWithLines(x + "y = SUM(x[3..1])\n") ==
              std::vector<std::string>{"4: array range 'x[3..1]': the first index (3) is greater than the last (1)"});
        CHECK(parseErrorsWithLines(x + "y = SUM(x[1..2.5])\n") ==
              std::vector<std::string>{"4: array range 'x[1..2.5]': the last index is 2.5, not an integer"});
        CHECK(parseErrorsWithLines(x + "n = 1.5\ny = SUM(x[n..2])\n") ==
              std::vector<std::string>{"5: array range 'x[n..2]': the first index is 1.5, not an integer"});
        CHECK_FALSE(parseErrors(x + "y = SUM(x[1...3])\n").empty());
        CHECK_FALSE(parseErrors(x + "y = SUM(x[1..])\n").empty());
    }
    SECTION("a bound that is not known before the solve") {
        auto errors = parseErrorsWithLines(x + "n = y + 1\nz = SUM(x[1..n])\n");
        REQUIRE(errors.size() == 1);
        CHECK(errors[0].rfind("5: array range 'x[1..n]'", 0) == 0);
        CHECK(anyContains(errors, "'n = y + 1' (line 4)"));
        errors = parseErrorsWithLines(x + "z = SUM(x[1..nn])\n");
        REQUIRE(errors.size() == 1);
        CHECK(anyContains(errors, "'nn' is not set"));
    }
    SECTION("a range that is not an argument") {
        for (const std::string line : {"y = x[1..3]\n", "y = 2*SUM(x[1..2]) + x[1..2]\n", "x[4..5] = 5\n",
                                       "y = SUM(x[1..2]^2)\n", "y = x[x[1..2]]\n"}) {
            INFO(line);
            const auto errors = parseErrorsWithLines(x + line);
            REQUIRE(errors.size() == 1);
            CHECK(errors[0].rfind("4: array range", 0) == 0);
            CHECK(anyContains(errors, "is only allowed as an argument"));
        }
    }
    SECTION("bounds that depend on the index of a DUPLICATE") {
        const auto errors = parseErrorsWithLines(x + "DUPLICATE i=1,3\n  y[i] = SUM(x[1..i])\nEND\n");
        REQUIRE(errors.size() == 1);
        CHECK(errors[0].rfind("5: array range 'x[1..i]'", 0) == 0);
        CHECK(anyContains(errors, "index 'i' of the enclosing DUPLICATE"));
    }
    SECTION("in a body and in a header, only numbers and $CONSTANTs") {
        auto errors = parseErrorsWithLines("FUNCTION g(n, X[1..3])\n  g = SUM(X[1..n])\nEND\nz = g(2, 1, 2, 3)\n");
        REQUIRE(errors.size() == 1);
        CHECK(errors[0].rfind("2: array range 'X[1..n]' in FUNCTION 'g'", 0) == 0);
        CHECK(anyContains(errors, "numbers or $CONSTANTs"));
        errors = parseErrorsWithLines("FUNCTION h(n, A[1..n])\n  h = A[1]\nEND\nz = h(1, 5)\n");
        REQUIRE(errors.size() == 1);
        CHECK(errors[0].rfind("1: array range 'A[1..n]' in the header of FUNCTION 'h'", 0) == 0);
        CHECK(anyContains(errors, "numbers or $CONSTANTs"));
    }
}

TEST_CASE("CS-GAP-ARRAY-ARG: the $CONSTANT directive", "[library-gaps][array-arg]") {
    SECTION("a value usable in the main program and in the bodies, defined above or below") {
        auto run = runModel(R"(
$CONSTANT N#=4
$Constant M# = N#-1 [m]
$CONSTANT L#=2 [m] 6.56 [ft]   "the English value is ignored"
FUNCTION f(x)
  f = x*N# + P#
END
y = 2*N#
z = f(3)
w = M#*L#
$constant P#=0.5
)");
        REQUIRE(run.parseOk);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK(run.warnings.empty());   // no "Unknown directive"
        CHECK_THAT(run["y"], WithinAbs(8.0, 1e-12));
        CHECK_THAT(run["z"], WithinAbs(12.5, 1e-12));
        CHECK_THAT(run["w"], WithinAbs(6.0, 1e-12));
        CHECK_FALSE(hasVar(run, "N#"));   // a constant, not a variable
    }
    SECTION("errors") {
        CHECK(anyContains(parseErrorsWithLines("$CONSTANT N#=3\nN# = 5\n"), "2: 'N#' is a constant ($CONSTANT, line 1)"));
        CHECK(anyContains(parseErrorsWithLines("$CONSTANT N#=3\nFUNCTION f(x)\n  N# = x\n  f = N#\nEND\ny = f(1)\n"),
                          "3: 'N#' is a constant ($CONSTANT, line 1)"));
        CHECK(anyContains(parseErrorsWithLines("$CONSTANT N#=3\nPROCEDURE p(x : y)\n  y = x\nEND\nCALL p(1 : N#)\n"),
                          "5: 'N#' is a constant ($CONSTANT, line 1)"));
        CHECK(anyContains(parseErrorsWithLines("$CONSTANT N# = 2 * 3\n"), "1: $CONSTANT"));
        CHECK(anyContains(parseErrorsWithLines("x = 2\n$CONSTANT N#=x+1\n"), "2: $CONSTANT N#: 'x' is not a constant"));
        CHECK(anyContains(parseErrorsWithLines("$CONSTANT N=3\n"), "1: $CONSTANT: the name of a constant ends with #"));
        CHECK(anyContains(parseErrorsWithLines("$CONSTANT R#$='R134a'\n"), "1: $CONSTANT R#$: string constants are not supported"));
        CHECK(anyContains(parseErrorsWithLines("$CONSTANT N#=3\n$CONSTANT N#=4\n"), "2: $CONSTANT N# is already defined (line 1)"));
    }
}

// ============================================================================
// CS-BUG-DUPLICATE-VAR-BOUND: the bounds of a DUPLICATE loop of the main program may be variables set to a
// constant expression (`n = 17`, `n = 2*m`), $CONSTANTs and arithmetic of those (EES: "previously set"); a
// bound that cannot be known before the solve is an error with its line. The loop was silently dropped.
// ============================================================================

TEST_CASE("CS-BUG-DUPLICATE-VAR-BOUND: reproducer of the register", "[library-gaps][duplicate-var-bound]") {
    // was: Equations: 1, SUCCESS (the loop dropped); EES: 3 equations, x[1] = 2, x[2] = 4
    const std::string code = "n=2\nDUPLICATE i=1,n\nx[i]=2*i\nEND\n";
    CHECK(equationCount(code) == 3);
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["n"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["x[1]"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["x[2]"], WithinAbs(4.0, 1e-12));
}

TEST_CASE("CS-BUG-DUPLICATE-VAR-BOUND: n = 17, a chain of variables, a $CONSTANT, expressions", "[library-gaps][duplicate-var-bound]") {
    SECTION("n = 17 as in CSL-0157, several equations per row") {
        const std::string code = "n=17\nduplicate i=1,n\n  a[i] = i\n  b[i] = a[i]^2\nend\ntot = sum(b[1..n])\n";
        CHECK(equationCount(code) == 1 + 2 * 17 + 1);
        auto run = runModel(code);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["b[17]"], WithinAbs(289.0, 1e-12));
        CHECK_THAT(run["tot"], WithinAbs(1785.0, 1e-9));
    }
    SECTION("n = 2*m, m = 3, defined after the loop") {
        const std::string code = "DUPLICATE i=1,n\n  x[i] = i\nEND\nn = 2*m\nm = 3\n";
        CHECK(equationCount(code) == 6 + 2);
        auto run = runModel(code);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x[6]"], WithinAbs(6.0, 1e-12));
    }
    SECTION("a $CONSTANT, alone and in an expression") {
        // (no index arithmetic such as x[i+1] in the loop: that is CS-BUG-DUPLICATE-INDEX-EXPR)
        const std::string code = "$CONSTANT N#=4\nDUPLICATE i=1,N#\n  x[i] = i\nEND\nDUPLICATE i = 1, (N#-1)\n"
                                 "  d[i] = 2*x[i]\nEND\n";
        CHECK(equationCount(code) == 4 + 3);
        auto run = runModel(code);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["d[3]"], WithinAbs(6.0, 1e-12));
        CHECK_FALSE(hasVar(run, "x[5]"));
    }
    SECTION("a start bound that is a variable, an end bound that is an expression") {
        const std::string code = "k = 2\nDUPLICATE i=k,2*k+1\n  x[i] = i\nEND\n";
        CHECK(equationCount(code) == 1 + 4);
        auto run = runModel(code);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["x[2]"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["x[5]"], WithinAbs(5.0, 1e-12));
    }
}

TEST_CASE("CS-BUG-DUPLICATE-VAR-BOUND: nested loops whose bounds use the outer index", "[library-gaps][duplicate-var-bound]") {
    // `j = i+1` was dropped too (an expression, not a number, once i is substituted)
    const std::string code = "n = 4\nDUPLICATE i=1,n-1\n  DUPLICATE j=i+1,n\n    c[i,j] = 10*i + j\n  END\nEND\n";
    CHECK(equationCount(code) == 1 + 6);
    auto run = runModel(code);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["c[1,2]"], WithinAbs(12.0, 1e-12));
    CHECK_THAT(run["c[3,4]"], WithinAbs(34.0, 1e-12));
}

TEST_CASE("CS-BUG-DUPLICATE-VAR-BOUND: a DUPLICATE in a body takes its bounds from the inputs when it runs",
          "[library-gaps][duplicate-var-bound]") {
    auto run = runModel("$CONSTANT M#=3\nFUNCTION tri(n)\n  s := 0\n  DUPLICATE k=1,n\n    s := s + k\n  END\n"
                        "  DUPLICATE k=1,M#\n    s := s + 100\n  END\n  tri = s\nEND\ny = tri(4)\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(310.0, 1e-12));
}

TEST_CASE("CS-BUG-DUPLICATE-VAR-BOUND: a bound that is not known before the solve is an error, not a dropped loop",
          "[library-gaps][duplicate-var-bound]") {
    const std::string loop = "DUPLICATE i=1,n\n  x[i] = i\nEND\n";
    SECTION("n set by an equation that depends on an unknown") {
        const auto errors = parseErrorsWithLines("n = y + 1\n" + loop);
        REQUIRE(errors.size() == 1);
        CHECK(errors[0].rfind("2: DUPLICATE i: the last bound 'n'", 0) == 0);
        CHECK(anyContains(errors, "'n = y + 1' (line 1)"));
        auto run = runModel("n = y + 1\n" + loop);
        CHECK_FALSE(run.parseOk);
        CHECK_FALSE(run.solveOk);
    }
    SECTION("n set by nothing, by a CALL, by an equation that is not of the form n = ...") {
        CHECK(anyContains(parseErrorsWithLines(loop), "1: DUPLICATE i: the last bound 'n'"));
        CHECK(anyContains(parseErrorsWithLines(loop), "'n' is not set"));
        CHECK(anyContains(parseErrorsWithLines("PROCEDURE p(a : b)\n  b = a\nEND\nCALL p(3 : n)\n" + loop),
                          "5: DUPLICATE i: the last bound 'n'"));
        CHECK(anyContains(parseErrorsWithLines("2*n = 6\n" + loop), "2: DUPLICATE i: the last bound 'n'"));
    }
    SECTION("a non-integer bound, a cycle, an array element") {
        CHECK(parseErrorsWithLines("n = 2.5\n" + loop) ==
              std::vector<std::string>{"2: DUPLICATE i: the last bound 'n' is 2.5, not an integer"});
        CHECK(anyContains(parseErrorsWithLines("n = m\nm = n + 0\n" + loop), "3: DUPLICATE i: the last bound 'n'"));
        CHECK(anyContains(parseErrorsWithLines("N[1] = 3\nDUPLICATE i=1,N[1]\n  x[i] = i\nEND\n"),
                          "2: DUPLICATE i: the last bound 'N[1]'"));
    }
    SECTION("the first bound too") {
        CHECK(anyContains(parseErrorsWithLines("DUPLICATE i=k,3\n  x[i] = i\nEND\n"), "1: DUPLICATE i: the first bound 'k'"));
    }
    SECTION("a nested bound that is not an integer in one of the copies") {
        // known only once the outer index is substituted (IR): reported like a parse error, with its line
        auto run = runModel("n = 3\nDUPLICATE i=1,n\n  DUPLICATE j=1,i/2\n    x[i,j] = 1\n  END\nEND\n");
        CHECK_FALSE(run.parseOk);
        CHECK_FALSE(run.solveOk);
        CHECK(anyContains({run.message}, "Line 3: DUPLICATE j: the last bound is 0.5 in a copy of the enclosing loop, "
                                         "not an integer"));
    }
}

// ============================================================================
// CS-GAP-ELSEIF-CHAIN: IF ladders closed by `;`-joined ENDIFs, and single-line IFs
// ending with a bare ELSE whose ELSE part runs on the following lines up to its ENDIF
// (Brineprop.lib of the ULiege collection, library model CSL-0079)
// ============================================================================

TEST_CASE("CS-GAP-ELSEIF-CHAIN: register reproducer, nested block IFs closed by ENDIF;ENDIF",
          "[library-gaps][elseif-chain]") {
    // was: "Line 2: IF ... THEN without a matching ENDIF" (the whole `ENDIF;ENDIF` line closed one IF)
    const std::string proc =
        "PROCEDURE p(x : y)\n"
        "  IF (x>0) THEN\n"
        "    IF (x>1) THEN\n"
        "      y=20\n"
        "    ELSE\n"
        "      y=10\n"
        "    ENDIF;ENDIF\n"
        "END\n";
    REQUIRE(parseErrorsWithLines(proc).empty());
    auto run = runModel(proc + "CALL p(2 : a)\nCALL p(1 : b)\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(20.0, 1e-12));   // EES: x>1
    CHECK_THAT(run["b"], WithinAbs(10.0, 1e-12));   // EES: 0<x<=1
}

TEST_CASE("CS-GAP-ELSEIF-CHAIN: ladders of block IFs closed by ENDIF;ENDIF;ENDIF in a FUNCTION",
          "[library-gaps][elseif-chain]") {
    SECTION("ELSE and IF on separate lines, one ENDIF per IF on one line") {
        auto run = runModel(
            "FUNCTION lad(x)\n"
            "  IF (x<1) THEN\n"
            "    lad = 1\n"
            "  ELSE\n"
            "  IF (x<2) THEN\n"
            "    lad = 2\n"
            "  ELSE\n"
            "  IF (x<3) THEN\n"
            "    lad = 3\n"
            "  ELSE\n"
            "    lad = 4\n"
            "  ENDIF;ENDIF;ENDIF\n"
            "END\n"
            "a = lad(0.5)\nb = lad(1.5)\nc = lad(2.5)\nd = lad(3.5)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(3.0, 1e-12));
        CHECK_THAT(run["d"], WithinAbs(4.0, 1e-12));
    }
    SECTION("ELSE IF (c) THEN on one line, spaces and a comment around the ENDIFs, any letter case") {
        auto run = runModel(
            "FUNCTION lad(x)\n"
            "  If (x<1) Then\n"
            "    lad = 1\n"
            "  Else If (x<2) Then\n"
            "    lad = 2\n"
            "  ELSE IF (x<3) THEN\n"
            "    lad = 3\n"
            "  ELSE IF (x<4) THEN\n"
            "    lad = 4\n"
            "  ELSE\n"
            "    lad = 5\n"
            "  endif ; ENDIF;  Endif ;ENDIF  \"end of the ladder\"\n"
            "END\n"
            "a = lad(0.5)\nb = lad(1.5)\nc = lad(2.5)\nd = lad(3.5)\ne = lad(4.5)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(3.0, 1e-12));
        CHECK_THAT(run["d"], WithinAbs(4.0, 1e-12));
        CHECK_THAT(run["e"], WithinAbs(5.0, 1e-12));
    }
}

TEST_CASE("CS-GAP-ELSEIF-CHAIN: ladders closed by ENDIF;ENDIF in a PROCEDURE, with a CALL in the branches",
          "[library-gaps][elseif-chain]") {
    // layout of the first ladder of BRINEPROP: THEN part on the next line, ELSE on its own line,
    // one ENDIF per IF on the last line, the last IF without ELSE
    const std::string code =
        "PROCEDURE chk(t$, v : r)\n"
        "  r = 0\n"
        "IF (t$='A') AND (v<0) THEN\n"
        "\tCall ERROR('A: v must be >= 0')\n"
        "\tELSE\n"
        "\tIF (t$='A') AND (v>10) THEN\n"
        "\tCall ERROR('A: v must be <= 10')\n"
        "\tELSE\n"
        "IF (t$='B') AND (v<5) THEN\n"
        "\tCall ERROR('B: v must be >= 5')\n"
        "ENDIF;ENDIF;ENDIF\n"
        "  r = v\n"
        "END\n";
    REQUIRE(parseErrorsWithLines(code).empty());
    auto run = runModel(code + "CALL chk('A', 3 : a)\nCALL chk('B', 6 : b)\nCALL chk('C', -1 : c)\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(6.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(-1.0, 1e-12));
    CHECK_FALSE(runModel(code + "CALL chk('A', -3 : a)\n").solveOk);   // the CALL ERROR branch is taken
    CHECK_FALSE(runModel(code + "CALL chk('A', 11 : a)\n").solveOk);
    CHECK_FALSE(runModel(code + "CALL chk('B', 4 : b)\n").solveOk);
}

TEST_CASE("CS-GAP-ELSEIF-CHAIN: ENDIF; followed by statements, ELSE or END", "[library-gaps][elseif-chain]") {
    SECTION("statements after ENDIF;ENDIF; run unconditionally") {
        auto run = runModel(
            "FUNCTION f(x)\n"
            "  s = 0\n"
            "  IF (x>0) THEN\n"
            "    IF (x>1) THEN\n"
            "      s = 10\n"
            "    ENDIF;ENDIF; s = s + 1; f = s\n"
            "END\n"
            "a = f(2)\nb = f(1)\nc = f(0)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(11.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(1.0, 1e-12));
    }
    SECTION("ENDIF;ELSE: the ELSE belongs to the outer IF") {
        auto run = runModel(
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN\n"
            "    f = 1\n"
            "    IF (x>1) THEN\n"
            "      f = 10\n"
            "    ENDIF;ELSE\n"
            "    f = -1\n"
            "  ENDIF\n"
            "END\n"
            "a = f(2)\nb = f(1)\nc = f(-3)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(10.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(-1.0, 1e-12));
    }
    SECTION("ENDIF;END ends the function, as if END were on the next line") {
        auto run = runModel(
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN f = 1 ELSE\n"
            "    f = 2\n"
            "  ENDIF;END\n"
            "a = f(2)\nb = f(-1)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(2.0, 1e-12));
    }
    SECTION("a multi-line comment opened after ENDIF;ENDIF does not swallow the second ENDIF") {
        for (const char* open : {"{note on the ladder\n  second line of the note}", "\"note on the ladder\n  second line of the note\""}) {
            auto run = runModel(
                std::string("FUNCTION f(x)\n"
                            "  f = 0\n"
                            "  IF (x>0) THEN\n"
                            "    IF (x>1) THEN f = 2 ELSE\n"
                            "      f = 1\n"
                            "    ENDIF;ENDIF ") + open + "\n"
                            "  f = f + 10\n"
                            "END\n"
                            "a = f(5)\nb = f(0.5)\nc = f(-1)\n");
            REQUIRE(run.solveOk);
            CHECK_THAT(run["a"], WithinAbs(12.0, 1e-12));
            CHECK_THAT(run["b"], WithinAbs(11.0, 1e-12));
            CHECK_THAT(run["c"], WithinAbs(10.0, 1e-12));
        }
    }
    SECTION("an empty segment between the ENDIFs is allowed") {
        auto run = runModel(
            "FUNCTION f(x)\n"
            "  f = 0\n"
            "  IF (x>0) THEN\n"
            "    IF (x>1) THEN f = 2 ELSE\n"
            "      f = 1\n"
            "    ENDIF;;ENDIF;\n"
            "END\n"
            "a = f(5)\nb = f(0.5)\nc = f(-1)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("CS-GAP-ELSEIF-CHAIN: single-line IF ending with a bare ELSE, string conditions as in BrineProp",
          "[library-gaps][elseif-chain]") {
    // was: the empty ELSE part was taken as empty, the next lines ran unconditionally and the
    // ENDIFs were orphans
    const std::string proc =
        "PROCEDURE sel(Pr$ : Pro)\n"
        "Pro=0\n"
        "Freeze=1;Density=2;SpecHeat=3;ThermalC=4;DynVisc=5\n"
        "IF (Pr$='FREEZE') THEN Pro:=Freeze ELSE\n"
        "\tIF (Pr$='DENSITY') THEN Pro:=Density ELSE\n"
        "\t\tIF (Pr$='SPECHEAT') THEN Pro:=SpecHeat ELSE\n"
        "\t\t\tIF (Pr$='THERMALC') THEN Pro:=ThermalC ELSE\n"
        "\t\t\t\tIF (Pr$='DYNVISC') THEN Pro:=DynVisc\n"
        "ENDIF;ENDIF;ENDIF;ENDIF\n"
        "END\n";
    REQUIRE(parseErrorsWithLines(proc).empty());
    auto run = runModel(proc +
        "CALL sel('FREEZE' : a1)\nCALL sel('DENSITY' : a2)\nCALL sel('SPECHEAT' : a3)\n"
        "CALL sel('THERMALC' : a4)\nCALL sel('DYNVISC' : a5)\nCALL sel('NONE' : a6)\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a1"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["a2"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["a3"], WithinAbs(3.0, 1e-12));
    CHECK_THAT(run["a4"], WithinAbs(4.0, 1e-12));
    CHECK_THAT(run["a5"], WithinAbs(5.0, 1e-12));
    CHECK_THAT(run["a6"], WithinAbs(0.0, 1e-12));   // no branch taken: Pro keeps its initial value
}

TEST_CASE("CS-GAP-ELSEIF-CHAIN: the ELSE part of a bare ELSE runs to its ENDIF", "[library-gaps][elseif-chain]") {
    SECTION("several statements in the ELSE part, comments after ELSE and after ENDIF, statements after") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN r = 1; s = 1 ELSE {else part below}\n"
            "    r = 2\n"
            "    s = 20\n"
            "  ENDIF \"end\"\n"
            "  f = r + s\n"
            "END\n"
            "a = f(1)\nb = f(-1)\n";
        REQUIRE(parseErrorsWithLines(code).empty());
        auto run = runModel(code);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(22.0, 1e-12));
    }
    SECTION("a block IF with its own ELSE and ENDIF inside the ELSE part") {
        auto run = runModel(
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN f = 1 ELSE\n"
            "    IF (x<-5) THEN\n"
            "      f = -2\n"
            "    ELSE\n"
            "      f = -1\n"
            "    ENDIF\n"
            "  ENDIF\n"
            "END\n"
            "a = f(3)\nb = f(-1)\nc = f(-9)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(-1.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(-2.0, 1e-12));
    }
    SECTION("single-line IFs with a trailing ELSE and block IFs mixed, closed by one ENDIF line") {
        auto run = runModel(
            "FUNCTION f(x)\n"
            "  IF (x=1) THEN f = 10 ELSE\n"
            "  IF (x=2) THEN\n"
            "    f = 20\n"
            "  ELSE\n"
            "    IF (x=3) THEN f = 30 ELSE\n"
            "      f = 99\n"
            "  ENDIF;ENDIF;ENDIF\n"
            "END\n"
            "a = f(1)\nb = f(2)\nc = f(3)\nd = f(4)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(10.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(20.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(30.0, 1e-12));
        CHECK_THAT(run["d"], WithinAbs(99.0, 1e-12));
    }
    SECTION("the single-line forms without a trailing ELSE are unchanged") {
        auto run = runModel(
            "FUNCTION f(x)\n"
            "  f = 0\n"
            "  IF (x>0) THEN f = 1\n"
            "  IF (x>1) THEN f = 2 ELSE f = f + 10\n"
            "END\n"
            "a = f(3)\nb = f(0.5)\nc = f(-1)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(2.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(11.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(10.0, 1e-12));
    }
}

TEST_CASE("CS-GAP-ELSEIF-CHAIN: ladders inside DUPLICATE and REPEAT loops of a body", "[library-gaps][elseif-chain]") {
    auto run = runModel(
        "FUNCTION g(n)\n"
        "  s = 0\n"
        "  DUPLICATE i = 1, n\n"
        "    IF (i > 2) THEN\n"
        "      IF (i > 3) THEN\n"
        "        s = s + 10\n"
        "      ELSE\n"
        "        s = s + 1\n"
        "      ENDIF;ENDIF\n"
        "  END\n"
        "  g = s\n"
        "END\n"
        "FUNCTION h(n)\n"
        "  s = 0\n"
        "  k = 0\n"
        "  REPEAT\n"
        "    k = k + 1\n"
        "    IF (k > 2) THEN s = s + 1 ELSE\n"
        "      IF (k > 1) THEN s = s + 100\n"
        "    ENDIF;\n"
        "  UNTIL (k >= n)\n"
        "  h = s\n"
        "END\n"
        "a = g(5)\nb = h(5)\n");
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(21.0, 1e-12));    // i = 3: +1, i = 4, 5: +10 each
    CHECK_THAT(run["b"], WithinAbs(103.0, 1e-12));   // k = 2: +100, k = 3, 4, 5: +1 each
}

namespace {

/// "line: message" of the P008 warnings (an ENDIF that closes no IF) of the diagnostics `items`.
std::vector<std::string> ignoredEndifWarnings(const std::vector<Diagnostic>& items) {
    std::vector<std::string> out;
    for (const auto& d : items)
        if (d.code == "P008") {
            CHECK(d.severity == DiagnosticSeverity::Warning);
            out.push_back(std::to_string(d.line) + ": " + d.message);
        }
    return out;
}

}  // namespace

TEST_CASE("CS-GAP-ELSEIF-CHAIN: an ENDIF that closes no IF in a body is a warning, as EES accepts it",
          "[library-gaps][elseif-chain]") {
    const std::string orphan = "ENDIF without a matching IF ... THEN: ignored (EES accepts it)";
    SECTION("one ENDIF too many on a joined line: one warning at its line, the model solves") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN\n"
            "    f = 1\n"
            "  ELSE\n"
            "    f = 2\n"
            "  ENDIF;ENDIF\n"
            "END\n"
            "a = f(1)\nb = f(-1)\n";
        const auto parsed = parseFull(code);
        CHECK(parsed.errors.empty());
        CHECK(ignoredEndifWarnings(parsed.diagnostics.items()) == std::vector<std::string>{"6: " + orphan});
        auto run = runModel(code);
        REQUIRE(run.parseOk);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(2.0, 1e-12));
        CHECK(ignoredEndifWarnings(run.warningItems) == std::vector<std::string>{"6: " + orphan});
    }
    SECTION("ENDIF;ENDIF with no open IF: one warning for the line, not one per ENDIF; one per line") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  f = 1\n"
            "  ENDIF;ENDIF\n"
            "  f = f + x\n"
            "  ENDIF\n"
            "END\n"
            "a = f(1)\n";
        const auto parsed = parseFull(code);
        CHECK(parsed.errors.empty());
        CHECK(ignoredEndifWarnings(parsed.diagnostics.items()) ==
              std::vector<std::string>{"3: " + orphan, "5: " + orphan});
        auto run = runModel(code);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(2.0, 1e-12));
    }
    SECTION("the surplus ENDIF of a loop in a body is a warning too, and statements after it run") {
        auto run = runModel(
            "FUNCTION g(n)\n"
            "  s = 0\n"
            "  DUPLICATE i = 1, n\n"
            "    s = s + i\n"
            "    ENDIF; s = s + 100\n"
            "  END\n"
            "  g = s\n"
            "END\n"
            "a = g(3)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(306.0, 1e-12));
        CHECK(ignoredEndifWarnings(run.warningItems) == std::vector<std::string>{"5: " + orphan});
    }
    SECTION("the two-THEN single line of the register: its orphan ENDIFs are ignored with one warning") {
        // The help of EES documents one THEN on a single-line IF; this line is read as nested single-line
        // IFs, and what closes nothing is reported (it used to be two errors).
        const std::string code =
            "PROCEDURE p(x : y)\n"
            "IF (x=1) THEN y=10 ELSE IF (x=2) THEN y=20\n"
            "ENDIF;ENDIF\n"
            "END\n"
            "CALL p(2 : yv)\nz=yv\n";
        auto run = runModel(code);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["z"], WithinAbs(20.0, 1e-12));
        CHECK(ignoredEndifWarnings(run.warningItems) == std::vector<std::string>{"3: " + orphan});
    }
    SECTION("reduced copies of the BrineProp ladders, which have more ENDIFs than IFs") {
        // 3 block IFs closed by 4 ENDIFs (23 for 22 in the library), 3 single-line IFs of which two end
        // with ELSE closed by 3 ENDIFs (11 for 10 in the library)
        const std::string code =
            "PROCEDURE brine(Fl$, Conc : output)\n"
            "IF (Fl$='EG') AND (Conc<0) THEN\n"
            "\tCall ERROR('EG: Conc must be >= 0')\n"
            "\tELSE\n"
            "\tIF (Fl$='EG') AND (Conc>56.1) THEN\n"
            "\tCall ERROR('EG: Conc must be <= 56.1')\n"
            "\tELSE\n"
            "IF (Fl$='PG') AND (Conc<15.2) THEN\n"
            "\tCall ERROR('PG: Conc must be >= 15.2')\n"
            "ENDIF;ENDIF;ENDIF;ENDIF\n"
            "IF (Fl$='EG') THEN Flu:=1 ELSE\n"
            "\tIF (Fl$='PG') THEN Flu:=2 ELSE\n"
            "\t\tIF (Fl$='EA') THEN Flu:=3\n"
            "ENDIF;ENDIF;ENDIF\n"
            "output = Flu*100 + Conc\n"
            "END\n";
        const auto parsed = parseFull(code);
        CHECK(parsed.errors.empty());
        CHECK(ignoredEndifWarnings(parsed.diagnostics.items()) ==
              std::vector<std::string>{"10: " + orphan, "14: " + orphan});
        auto run = runModel(code + "CALL brine('EG', 25 : a)\nCALL brine('PG', 30 : b)\nCALL brine('EA', 5 : c)\n");
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(125.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(230.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(305.0, 1e-12));
        CHECK(ignoredEndifWarnings(run.warningItems) == std::vector<std::string>{"10: " + orphan, "14: " + orphan});
        CHECK_FALSE(runModel(code + "CALL brine('EG', -1 : a)\n").solveOk);        // the CALL ERROR branches
        CHECK_FALSE(runModel(code + "CALL brine('EG', 60 : a)\n").solveOk);
        CHECK_FALSE(runModel(code + "CALL brine('PG', 10 : a)\n").solveOk);
    }
}

TEST_CASE("CS-GAP-ELSEIF-CHAIN: unmatched ENDIFs and ELSEs: what stays an error", "[library-gaps][elseif-chain]") {
    SECTION("one ENDIF too few: the outer IF is reported at its own line") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  f = 0\n"
            "  IF (x>0) THEN\n"
            "    IF (x>1) THEN\n"
            "      f = 1\n"
            "    ENDIF\n"
            "END\n"
            "a = f(1)\n";
        CHECK(parseErrorsWithLines(code) ==
              std::vector<std::string>{"3: IF ... THEN without a matching ENDIF"});
        CHECK_FALSE(runModel(code).parseOk);
    }
    SECTION("one ENDIF too few on a joined line") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  f = 0\n"
            "  IF (x>0) THEN\n"
            "    IF (x>1) THEN\n"
            "      IF (x>2) THEN\n"
            "        f = 1\n"
            "      ENDIF;ENDIF\n"
            "END\n"
            "a = f(1)\n";
        CHECK(parseErrorsWithLines(code) ==
              std::vector<std::string>{"3: IF ... THEN without a matching ENDIF"});
    }
    SECTION("an ENDIF with no IF outside a body (DUPLICATE of the main program): one error per ENDIF") {
        const std::string code =
            "DUPLICATE i = 1, 2\n"
            "  x[i] = i\n"
            "  ENDIF;ENDIF\n"
            "END\n";
        const auto parsed = parseFull(code);
        CHECK(parseErrorsWithLines(code) ==
              std::vector<std::string>{"3: ENDIF without a matching IF ... THEN",
                                       "3: ENDIF without a matching IF ... THEN"});
        CHECK(ignoredEndifWarnings(parsed.diagnostics.items()).empty());
    }
    SECTION("an ELSE with no IF is an error, in a body too") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  f = 1\n"
            "  ELSE\n"
            "  f = 2\n"
            "END\n"
            "a = f(1)\n";
        CHECK(parseErrorsWithLines(code) ==
              std::vector<std::string>{"3: ELSE without a matching IF ... THEN"});
    }
    SECTION("text after ENDIF that is not separated by ';' is an error, never dropped") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  f = 0\n"
            "  IF (x>0) THEN\n"
            "    f = 1\n"
            "  ENDIF f = 5\n"
            "END\n"
            "a = f(1)\n";
        CHECK(parseErrorsWithLines(code) ==
              std::vector<std::string>{"5: Unexpected text 'f = 5' after ENDIF (separate statements with ';')"});
        CHECK_FALSE(runModel(code).parseOk);
        // the same text after a ';', or a comment after the ENDIF, is fine
        CHECK(parseErrorsWithLines("FUNCTION f(x)\n  f = 0\n  IF (x>0) THEN\n    f = 1\n  ENDIF; f = f + 5\nEND\na = f(1)\n").empty());
        CHECK(parseErrorsWithLines("FUNCTION f(x)\n  f = 0\n  IF (x>0) THEN\n    f = 1\n  ENDIF // end\n  f = f + 5\nEND\na = f(1)\n").empty());
        CHECK(parseErrorsWithLines("FUNCTION f(x)\n  f = 0\n  IF (x>0) THEN\n    f = 1\n  ENDIF {end} \"done\"\n  f = f + 5\nEND\na = f(1)\n").empty());
        // `ENDIF ELSE` is no ENDIF followed by an ELSE either
        CHECK(anyContains(parseErrorsWithLines("FUNCTION f(x)\n  IF (x>0) THEN\n    IF (x>1) THEN f = 2 ELSE\n      f = 1\n    ENDIF ELSE\n    f = 0\n  ENDIF\nEND\n"),
                          "Unexpected text 'ELSE' after ENDIF"));
    }
    SECTION("a bare ELSE without ENDIF: the IF is reported at its line (it used to run the rest unconditionally)") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN f = 1 ELSE\n"
            "    f = 2\n"
            "END\n"
            "a = f(1)\n";
        CHECK(parseErrorsWithLines(code) ==
              std::vector<std::string>{"2: IF ... THEN without a matching ENDIF"});
    }
    SECTION("two ELSE parts") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN f = 1 ELSE\n"
            "    f = 2\n"
            "  ELSE\n"
            "    f = 3\n"
            "  ENDIF\n"
            "END\n"
            "a = f(1)\n";
        CHECK(parseErrorsWithLines(code) ==
              std::vector<std::string>{"4: IF statement with two ELSE parts"});
    }
    SECTION("the one-word ELSEIF is not EES and stays a parse error") {
        const std::string code =
            "FUNCTION f(x)\n"
            "  IF (x>0) THEN\n"
            "    f = 1\n"
            "  ELSEIF (x<0) THEN\n"
            "    f = -1\n"
            "  ENDIF\n"
            "END\n"
            "a = f(1)\n";
        const auto errors = parseErrorsWithLines(code);
        REQUIRE_FALSE(errors.empty());
        CHECK(errors[0] == "4: Could not parse line");
        CHECK_FALSE(runModel(code).parseOk);
    }
}

// ============================================================================
// CS-GAP-UPPERCASE: the string functions Uppercase$(S$) and Lowercase$(S$)
// ============================================================================

TEST_CASE("CS-GAP-UPPERCASE: reproducer of the register", "[library-gaps][uppercase]") {
    // was: warning "Unknown function 'Uppercase$'", then "String-valued function not supported: Uppercase$"
    auto run = runModel(
        "PROCEDURE pe(Pr$:y)\n"
        "w$=Uppercase$('ab')\n"
        "y=1\n"
        "END\n"
        "CALL pe('ab':yv)\n"
        "z=yv\n");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["yv"], WithinAbs(1.0, 1e-12));
    CHECK_THAT(run["z"], WithinAbs(1.0, 1e-12));
    CHECK_FALSE(anyContains(run.warnings, "Unknown function"));
    CHECK(parseErrors("w$ = Uppercase$('ab')\n").empty());
    // the case of the function name does not matter
    CHECK_FALSE(anyContains(runModel("PROCEDURE p(s$ : n)\nw$ = UPPERCASE$(s$)\nv$ = lowerCASE$(w$)\nn = 1\nEND\nCALL p('a' : x)\n").warnings,
                            "Unknown function"));
}

TEST_CASE("CS-GAP-UPPERCASE: the text that comes out of the two functions", "[library-gaps][uppercase]") {
    // the strings come out of the PROCEDURE as outputs: the comparison of EES ignores the case, so only the
    // text of a string variable shows what the function did
    auto run = runModel(R"(
PROCEDURE conv(s$ : U$, L$, n)
  U$ = Uppercase$(s$)
  L$ = Lowercase$(s$)
  n = 1
END
CALL conv('Brine-1 ab_c, 9%' : U1$, L1$, n1)
CALL conv('' : U2$, L2$, n2)
CALL conv('ÉcoLe à ß × ÷ ÿ ç' : U3$, L3$, n3)
CALL conv('already UPPER' : U4$, L4$, n4)
)");
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.strings.at("U1$") == "BRINE-1 AB_C, 9%");
    CHECK(run.strings.at("L1$") == "brine-1 ab_c, 9%");
    CHECK(run.strings.at("U2$") == "");                 // the empty string
    CHECK(run.strings.at("L2$") == "");
    // the accented letters of Latin-1 follow, as in the folding of the comparisons (ees_string.h); the
    // signs, the sharp s and y with diaeresis (whose capital is not in Latin-1) are kept
    CHECK(run.strings.at("U3$") == "ÉCOLE À ß × ÷ ÿ Ç");
    CHECK(run.strings.at("L3$") == "école à ß × ÷ ÿ ç");
    CHECK(run.strings.at("U4$") == "ALREADY UPPER");
    CHECK(run.strings.at("L4$") == "already upper");
}

TEST_CASE("CS-GAP-UPPERCASE: BrineProp pattern, an input of the PROCEDURE is reassigned", "[library-gaps][uppercase]") {
    // PROCEDURE BRINEPROP(Pr$, ... : ...): `Pr$=Uppercase$(Pr$)` then `IF (Pr$='FREEZE') THEN ...`
    const std::string code = R"(
PROCEDURE brine(Pr$, Conc : out, Name$)
  Pr$ = Uppercase$(Pr$)
  out = -1
  IF (Pr$ = 'FREEZE') THEN out = 10*Conc
  IF (Pr$ = 'DENSITY') THEN out = 1000 + Conc
  IF (Pr$ <> 'FREEZE') AND (Pr$ <> 'DENSITY') THEN CALL ERROR('unknown property XXX$', Pr$)
  Name$ = Pr$
END
f$ = 'freeze'
CALL brine(f$, 3 : a, A$)
CALL brine('Density', 3 : b, B$)
CALL brine('FREEZE', 3 : c, C$)
)";
    auto run = runModel(code);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a"], WithinAbs(30.0, 1e-12));
    CHECK_THAT(run["b"], WithinAbs(1003.0, 1e-12));
    CHECK_THAT(run["c"], WithinAbs(30.0, 1e-12));
    CHECK(run.strings.at("A$") == "FREEZE");
    CHECK(run.strings.at("B$") == "DENSITY");
    CHECK(run.strings.at("C$") == "FREEZE");
    CHECK(run.strings.at("f$") == "freeze");   // the input is a copy: the caller's string is not changed

    auto bad = runModel(code + "CALL brine('nonsense', 3 : d, D$)\n");
    INFO(bad.message);
    CHECK_FALSE(bad.solveOk);
    CHECK(bad.message.find("unknown property NONSENSE") != std::string::npos);   // the converted text
}

TEST_CASE("CS-GAP-UPPERCASE: in FUNCTION and PROCEDURE bodies: locals, IF operands, nesting, CALL arguments",
          "[library-gaps][uppercase]") {
    SECTION("a local string of a FUNCTION") {
        auto run = runModel(R"(
FUNCTION isEG(s$)
  w$ = Lowercase$(s$)
  isEG = 0
  IF (w$ = 'eg') THEN isEG = 1
END
a = isEG('EG')
b = isEG('Eg')
c = isEG('PG')
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(0.0, 1e-12));
    }
    SECTION("an operand of an IF comparison, either side, = and <>") {
        auto run = runModel(R"(
FUNCTION f(s$)
  f = 0
  IF (Uppercase$(s$) = 'EG') THEN f = f + 1
  IF ('eg' = Lowercase$(s$)) THEN f = f + 10
  IF (Uppercase$(s$) <> Lowercase$('XX')) THEN f = f + 100
  IF (Uppercase$(s$) = Lowercase$(s$)) THEN f = f + 1000
END
a = f('Eg')
b = f('xx')
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1111.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(1000.0, 1e-12));
    }
    SECTION("nested, and applied to the other string functions") {
        auto run = runModel(R"(
PROCEDURE nest(s$ : A$, B$, C$, D$, E$, n)
  A$ = Uppercase$(Lowercase$(s$))
  B$ = Lowercase$(Uppercase$(Lowercase$(s$)))
  C$ = Uppercase$(String$(2.5))
  D$ = Uppercase$(LOOKUP$('Names', 2, 1))
  E$ = Lowercase$(Uppercase$('aB'))
  n = 1
END
CALL nest('aBcD' : A$, B$, C$, D$, E$, n)
)", {{"Names", "Name\nnitrogen\nwater\n"}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.strings.at("A$") == "ABCD");
        CHECK(run.strings.at("B$") == "abcd");
        CHECK(run.strings.at("C$") == "2.5");
        CHECK(run.strings.at("D$") == "WATER");
        CHECK(run.strings.at("E$") == "ab");
    }
    SECTION("as the fluid of a property call, in a FUNCTION and a PROCEDURE") {
        auto run = runModel(R"(
FUNCTION h(f$, T, P)
  h = enthalpy(Uppercase$(f$), T=T, P=P)
END
PROCEDURE props(f$, T, P : rho, cpv)
  rho = density(Lowercase$(f$), T=T, P=P)
  cpv = cp(Uppercase$(Lowercase$(f$)), T=T, P=P)
END
h1 = h('water', 80, 2E5)
h2 = h('WaTeR', 80, 2E5)
h_ref = enthalpy(Water, T=80, P=2E5)
CALL props('WATER', 80, 2E5 : rho, cpv)
rho_ref = density(Water, T=80, P=2E5)
cp_ref = cp(Water, T=80, P=2E5)
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["h1"], WithinRel(run["h_ref"], 1e-12));
        CHECK_THAT(run["h2"], WithinRel(run["h_ref"], 1e-12));
        CHECK_THAT(run["rho"], WithinRel(run["rho_ref"], 1e-12));
        CHECK_THAT(run["cpv"], WithinRel(run["cp_ref"], 1e-12));
    }
    SECTION("as the fluid argument of PHASE$ and as the argument of a CALL, inside a body") {
        auto run = runModel(R"(
PROCEDURE inner(s$ : T$, n)
  T$ = s$
  n = 1
END
PROCEDURE outer(s$ : O$, g)
  CALL inner(Uppercase$(s$) : O$, n)
  g = 0
  IF (PHASE$(Uppercase$(s$), T=120, P=1E5) = 'gas') THEN g = 1
END
CALL outer('water' : A$, g)
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.strings.at("A$") == "WATER");
        CHECK_THAT(run["g"], WithinAbs(1.0, 1e-12));
    }
    SECTION("as the argument of a user FUNCTION and of STRINGVAL") {
        auto run = runModel(R"(
FUNCTION isEG(s$)
  isEG = 0
  IF (s$ = 'EG') THEN isEG = 1
END
f$ = 'eg'
a = isEG(Uppercase$(f$))
b = isEG(Lowercase$('XX'))
c = StringVal(Uppercase$('2.5e1'))
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(1.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(0.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(25.0, 1e-12));
    }
    SECTION("as the message of CALL WARNING and of CALL ERROR") {
        auto run = runModel(R"(
FUNCTION w(s$)
  CALL WARNING(Uppercase$('low text XXX$'), s$)
  CALL WARNING(Lowercase$(s$))
  w = 1
END
x = w('Value')
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        const auto messages = userWarningMessages(run);
        REQUIRE(messages.size() == 2);
        CHECK(messages[0] == "LOW TEXT Value");   // the converted text: the placeholder is the value of s$
        CHECK(messages[1] == "value");
        auto err = runModel("FUNCTION w(s$)\n  CALL ERROR(Uppercase$(s$))\n  w = 1\nEND\nx = w('stop here')\n");
        INFO(err.message);
        CHECK_FALSE(err.solveOk);
        CHECK(err.message.find("STOP HERE") != std::string::npos);
    }
}

TEST_CASE("CS-GAP-UPPERCASE: in the main program", "[library-gaps][uppercase]") {
    SECTION("as the fluid of a property call, literal or string variable") {
        auto run = runModel(R"(
f$ = 'water'
h1 = enthalpy(Uppercase$('water'), T=80, P=2E5)
h2 = enthalpy(Uppercase$(f$), T=80, P=2E5)
h3 = enthalpy(Lowercase$('WATER'), T=80, P=2E5)
h_ref = enthalpy(Water, T=80, P=2E5)
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["h1"], WithinRel(run["h_ref"], 1e-12));
        CHECK_THAT(run["h2"], WithinRel(run["h_ref"], 1e-12));
        CHECK_THAT(run["h3"], WithinRel(run["h_ref"], 1e-12));
    }
    SECTION("as an input of a CALL, the result a string output") {
        auto run = runModel(R"(
PROCEDURE copy(s$ : t$, n)
  t$ = s$
  n = 1
END
f$ = 'Water'
CALL copy(Uppercase$(f$) : A$, n1)
CALL copy(Lowercase$('AbC') : B$, n2)
)");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.strings.at("A$") == "WATER");
        CHECK(run.strings.at("B$") == "abc");
    }
    SECTION("a string assigned from the function is refused as for LOOKUP$: an error, never a silent 0") {
        // same rule as `g$ = LOOKUP$(...)` (CS-BUG-STRING-ALIAS): use it where the string is needed or in a body
        for (const std::string call : {"Uppercase$('ab')", "Lowercase$('AB')", "Uppercase$(f$)"}) {
            auto run = runModel("f$ = 'ab'\ng$ = " + call + "\nx = 1\n");
            INFO(call);
            INFO(run.message);
            CHECK_FALSE(run.solveOk);
            CHECK(run.message.find("returns a string and cannot be used as a number") != std::string::npos);
        }
        auto lookup = runModel("g$ = LOOKUP$('Species', 1, 1)\n", kSpeciesTable);
        CHECK_FALSE(lookup.solveOk);
        CHECK(lookup.message.find("returns a string and cannot be used as a number") != std::string::npos);
    }
}

TEST_CASE("CS-GAP-UPPERCASE: negative tests give clear errors", "[library-gaps][uppercase]") {
    auto fail = [&](const std::string& expr, const std::string& expected) {
        auto run = runModel("PROCEDURE p(s$, v : t$, n)\n  t$ = " + expr + "\n  n = v\nEND\nCALL p('a', 1 : A$, n)\n");
        INFO(expr);
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find(expected) != std::string::npos);
    };
    SECTION("wrong number of arguments") {
        fail("Uppercase$('a', 'b')", "Uppercase$(): expects 1 argument: a string");
        fail("Lowercase$(s$, s$)", "Lowercase$(): expects 1 argument: a string");
        fail("Uppercase$()", "Uppercase$(): expects 1 argument: a string");
        fail("Lowercase$(s$, T=1)", "Lowercase$ does not take named arguments ('T=1')");   // a parse error
    }
    SECTION("a number is not a string") {
        fail("Uppercase$(3)", "Uppercase$(): the argument must be a string");
        fail("Lowercase$(v)", "Lowercase$(): the argument must be a string");
        fail("Uppercase$(v + 1)", "Uppercase$(): the argument must be a string");
        fail("Uppercase$(sqrt(v))", "Uppercase$(): the argument must be a string");
    }
    SECTION("used as a number") {
        auto run = runModel("PROCEDURE p(s$ : n)\n  n = Uppercase$(s$)\nEND\nCALL p('a' : x)\n");
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("Uppercase$ returns a string and cannot be used as a number") != std::string::npos);
        auto top = runModel("x = Lowercase$('A') + 1\n");
        CHECK_FALSE(top.solveOk);
        CHECK(top.message.find("Lowercase$ returns a string and cannot be used as a number") != std::string::npos);
    }
    SECTION("an unknown fluid text is reported with the converted name") {
        auto run = runModel("h = enthalpy(Uppercase$('unobtainium'), T=300, P=1E5)\n");
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("Unknown fluid: 'UNOBTAINIUM'") != std::string::npos);
    }
    SECTION("a misspelt name still gives the unknown function warning") {
        auto run = runModel("PROCEDURE p(s$ : n)\n  w$ = Uppercas$(s$)\n  n = 1\nEND\nCALL p('a' : x)\n");
        CHECK(anyContains(run.warnings, "Unknown function 'Uppercas$'"));
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("String-valued function not supported: Uppercas$") != std::string::npos);
    }
}

// ============================================================================
// CS-GAP-INTERP-EES: the EES call form INTERPOLATE('Table','Col1','Col2',Col2=value)
// ============================================================================
// The named argument names the column that is the input; the other column of the call is returned
// (EES: the value of Col1 where Col2 = value). The positional form INTERPOLATE('t','xcol','ycol',x)
// of CoolSolve is unchanged. Was: "Unknown fluid: 'T'" (the call was evaluated as a property call).

namespace {
// not equally spaced, increasing in both columns: (1,10) (2,30) (4,70)
const std::map<std::string, std::string> kEesTable = {{"T", "A,B\n1,10\n2,30\n4,70\n"}};
// a 2-D grid (x, y in {1, 2}): z = 11, 12 (y = 1) and 21, 22 (y = 2)
const std::map<std::string, std::string> kEesGrid = {{"g", "x,y,z\n1,1,11\n2,1,12\n1,2,21\n2,2,22\n"}};
}  // namespace

TEST_CASE("CS-GAP-INTERP-EES: reproducer of the register", "[library-gaps][interp-ees]") {
    // Col1 = A, Col2 = B: the value of A where B = v (linear: B = 10, 30, 70 at A = 1, 2, 4)
    auto run = runModel(R"(
v = 50
a_v = INTERPOLATE('T','A','B',B=v)
a_lit = INTERPOLATE('T','A','B',B=20)
a_node = INTERPOLATE('T','A','B',B=30)
a_expr = INTERPOLATE('T','A','B',B=v/5+20)
)", kEesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["a_v"], WithinAbs(3.0, 1e-12));      // between (2,30) and (4,70)
    CHECK_THAT(run["a_lit"], WithinAbs(1.5, 1e-12));    // between (1,10) and (2,30)
    CHECK_THAT(run["a_node"], WithinAbs(2.0, 1e-12));
    CHECK_THAT(run["a_expr"], WithinAbs(2.0, 1e-12));   // B = 30: any expression is a value
}

TEST_CASE("CS-GAP-INTERP-EES: INTERPOLATE1 and the column order", "[library-gaps][interp-ees]") {
    auto run = runModel(R"(
a1 = INTERPOLATE1('T','A','B',B=20)
b1 = INTERPOLATE('T','A','B',A=3)
b2 = INTERPOLATE1('T','A','B',A=3)
b3 = INTERPOLATE('T','B','A',A=3)
a3 = INTERPOLATE('T','B','A',B=20)
)", kEesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["a1"], WithinAbs(1.5, 1e-12));   // INTERPOLATE1 takes the same form
    // the named column is the input whichever position it has in the call: A = 3 gives B = 50
    CHECK_THAT(run["b1"], WithinAbs(50.0, 1e-12));
    CHECK_THAT(run["b2"], WithinAbs(50.0, 1e-12));
    CHECK_THAT(run["b3"], WithinAbs(50.0, 1e-12));
    CHECK_THAT(run["a3"], WithinAbs(1.5, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EES: the EES form gives what the positional form gives", "[library-gaps][interp-ees]") {
    // the positional form is (x column, y column, x value): the EES form with the roles of the columns reversed
    auto run = runModel(R"(
a1 = INTERPOLATE('T','A','B',B=15)
p1 = INTERPOLATE('T','B','A',15)
a2 = INTERPOLATE('T','A','B',B=64)
p2 = INTERPOLATE('T','B','A',64)
a3 = INTERPOLATE('T','A','B',B=5)
p3 = INTERPOLATE('T','B','A',5)
a4 = INTERPOLATE('T','A','B',B=500)
p4 = INTERPOLATE('T','B','A',500)
)", kEesTable);
    REQUIRE(run.solveOk);
    for (const char* n : {"1", "2", "3", "4"}) {
        INFO("point " << n);
        CHECK(run[std::string("a") + n] == run[std::string("p") + n]);   // the same interpolation, bit for bit
    }
    // outside the table both forms extrapolate linearly (CS-GAP-INTERP-EXTRAP; flat before)
    CHECK_THAT(run["a3"], WithinAbs(1.0 + (5.0 - 10.0) / 20.0, 1e-12));
    CHECK_THAT(run["a4"], WithinAbs(4.0 + (500.0 - 70.0) / 20.0, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EES: names are case-insensitive and columns may be numbers", "[library-gaps][interp-ees]") {
    auto run = runModel(R"(
a1 = INTERPOLATE('T','a','B',b=20)
a2 = INTERPOLATE('T','A','b',B=20)
a3 = INTERPOLATE('T',1,2,B=20)
a4 = INTERPOLATE('T',2,1,A=3)
)", kEesTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    for (const char* n : {"a1", "a2", "a3"}) CHECK_THAT(run[n], WithinAbs(1.5, 1e-12));
    CHECK_THAT(run["a4"], WithinAbs(50.0, 1e-12));   // the second column of the call is named: A = 3 gives B
}

TEST_CASE("CS-GAP-INTERP-EES: the three call shapes of the library models", "[library-gaps][interp-ees]") {
    SECTION("CSL-0035: the variable has the name of the column and of the argument") {
        const std::map<std::string, std::string> tables = {
            {"lookup_1", "r_p,M_r,epsilon_s\n4.5,34.3,0.77\n4.7,33.8,0.79\n5.0,32.9,0.86\n"}};
        auto run = runModel(R"(
p_su_bar = 0.9
p_ex_bar = 4.6*p_su_bar
r_p = p_ex_bar/p_su_bar
epsilon_s = INTERPOLATE('lookup_1','epsilon_s','r_p',r_p=r_p)
M_r = INTERPOLATE('lookup_1','M_r','r_p',r_p=r_p)
)", tables);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["epsilon_s"], WithinAbs(0.77 + (0.79 - 0.77) * 0.5, 1e-9));
        CHECK_THAT(run["M_r"], WithinAbs(34.3 + (33.8 - 34.3) * 0.5, 1e-9));
    }
    SECTION("CSL-0011: a descending table, a literal value for the column") {
        const std::map<std::string, std::string> tables = {
            {"lookup_1", "N_rN,M_r_cp,r_p,epsilon_s_cp\n1,454,4.6,0.86\n0.95,420,4,0.865\n0.9,370,3.6,0.86\n"}};
        auto run = runModel(R"(
M_r_cp = INTERPOLATE('lookup_1','M_r_cp','N_rN',N_rN=0.95)
r_p = INTERPOLATE('lookup_1','r_p','N_rN',N_rN=0.95)
epsilon_s_cp = INTERPOLATE('lookup_1','epsilon_s_cp','N_rN',N_rN=0.925)
m2 = INTERPOLATE('lookup_1','N_rN','M_r_cp',0.95)
)", tables);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["M_r_cp"], WithinAbs(420.0, 1e-9));
        CHECK_THAT(run["r_p"], WithinAbs(4.0, 1e-9));
        CHECK_THAT(run["epsilon_s_cp"], WithinAbs(0.8625, 1e-9));
        CHECK_THAT(run["m2"], WithinAbs(420.0, 1e-9));   // the positional form of the `_ambient` variant: the same value
    }
    SECTION("gas engine: the argument name is the first-listed input column and the value is 0") {
        const std::map<std::string, std::string> tables = {
            {"lookup_1", "W_dot_sh,DELTAP_throttling\n18189.2927,5000\n17738.76861,6939\n17288.24354,8878\n"}};
        auto run = runModel("DELTAP_throttling_stb=INTERPOLATE('lookup_1','DELTAP_throttling','W_dot_sh',W_dot_sh=0)\n"
                            "dp_at = INTERPOLATE('lookup_1','DELTAP_throttling','W_dot_sh',W_dot_sh=17500)\n", tables);
        INFO(run.message);
        REQUIRE(run.solveOk);
        // W = 0 is below the table (descending column): linear extrapolation through the two last rows
        // (CS-GAP-INTERP-EXTRAP; the end value 8878 before)
        const double slopeLast = (8878.0 - 6939.0) / (17288.24354 - 17738.76861);
        CHECK_THAT(run["DELTAP_throttling_stb"], WithinAbs(8878.0 + slopeLast * (0.0 - 17288.24354), 1e-6));
        CHECK_THAT(run["dp_at"], WithinAbs(6939.0 + (8878.0 - 6939.0) * (17738.76861 - 17500.0) / (17738.76861 - 17288.24354), 1e-6));
    }
}

TEST_CASE("CS-GAP-INTERP-EES: the named value is an unknown solved by the block", "[library-gaps][interp-ees]") {
    // Newton needs the exact slope of the table: the value of the named argument counts as a
    // variable of the equation, as the positional x does
    SECTION("the value is the input, a second equation closes the block") {
        auto run = runModel(R"(
y = INTERPOLATE('T','A','B',B=x)
x + y = 53
)", kEesTable);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinAbs(50.0, 1e-9));   // A = 2 + (B - 30)/20
        CHECK_THAT(run["y"], WithinAbs(3.0, 1e-9));
    }
    SECTION("the named column is the first of the call") {
        // guess values inside the table: from the table end (x = 1, slope 0 there) Newton does not move, as with the positional form
        auto run = runModel(R"(
y = INTERPOLATE('T','B','A',A=x)
x + y = 53
)", kEesTable, "x = 3\ny = 40\n");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinAbs(3.0, 1e-9));   // B = 20 A - 10
        CHECK_THAT(run["y"], WithinAbs(50.0, 1e-9));
    }
    SECTION("inverse lookup: the call is a side of the equation") {
        auto run = runModel("INTERPOLATE('T','A','B',B=N) = 2.5\nINTERPOLATE('T','B','A',A=M) = 50\n", kEesTable);
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["N"], WithinAbs(40.0, 1e-9));   // B at A = 2.5: 30 + 40 * 0.25
        CHECK_THAT(run["M"], WithinAbs(3.0, 1e-9));
    }
    SECTION("the incidence: only the value counts, never the name of the argument") {
        EESParser parser;
        auto parsed = parser.parse("a = INTERPOLATE('T','A','B',B=v)\n");
        REQUIRE(parsed.success);
        IR ir = IR::fromAST(parsed.program);
        CHECK(ir.getAlgebraicVariableCount() == 2);
        CHECK(ir.getVariable("v") != nullptr);
        CHECK(ir.getVariable("a") != nullptr);
        CHECK(ir.getVariable("B") == nullptr);
    }
}

TEST_CASE("CS-GAP-INTERP-EES: the derivative is the slope of the table", "[library-gaps][interp-ees]") {
    // x column of the lookup t: 1 2 3; y: 10 20 30; z: 100 400 900 (r independent)
    auto second = evalWithSlope("v = INTERPOLATE('t','z','x',x=r)\n", kXyzCsv, 2.5);
    CHECK_THAT(second.value, WithinAbs(650.0, 1e-9));
    CHECK_THAT(second.slope, WithinAbs(500.0, 1e-9));
    auto first = evalWithSlope("v = INTERPOLATE('t','x','z',x=r)\n", kXyzCsv, 2.5);   // same input column, named first
    CHECK_THAT(first.value, WithinAbs(650.0, 1e-9));
    CHECK_THAT(first.slope, WithinAbs(500.0, 1e-9));
    // chain rule through the value expression: y = 10 x with x = 2 r
    auto chain = evalWithSlope("v = INTERPOLATE('t','y','x',x=2*r)\n", kXyzCsv, 1.25);
    CHECK_THAT(chain.value, WithinAbs(25.0, 1e-9));
    CHECK_THAT(chain.slope, WithinAbs(20.0, 1e-9));
    // outside the table the slope is the end slope, as for the positional form (CS-GAP-INTERP-EXTRAP; 0 before)
    auto outside = evalWithSlope("v = INTERPOLATE('t','y','x',x=r)\n", kXyzCsv, 5.0);
    CHECK_THAT(outside.value, WithinAbs(50.0, 1e-9));
    CHECK_THAT(outside.slope, WithinAbs(10.0, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EES: inside FUNCTION and PROCEDURE bodies", "[library-gaps][interp-ees]") {
    auto run = runModel(R"(
FUNCTION g(r)
  g = INTERPOLATE('t','y','x',x=2*r)
END
PROCEDURE pk(r : y, z)
  y = INTERPOLATE('t','z','x',x=r)
  z = INTERPOLATE1('t','x','y',y=r*10)
END
fv = g(1.25)
CALL pk(2.5 : yv, zv)
fp = INTERPOLATE('t','y','x',x=2.5)
)", kXyzTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["fv"], WithinAbs(25.0, 1e-9));
    CHECK_THAT(run["yv"], WithinAbs(650.0, 1e-9));
    CHECK_THAT(run["zv"], WithinAbs(2.5, 1e-9));    // x at y = 25
    CHECK_THAT(run["fp"], WithinAbs(25.0, 1e-9));

    // the output of a routine that reads the table is an unknown of the block (Newton, analytic slope)
    auto inverse = runModel(R"(
FUNCTION g(x)
  g = INTERPOLATE('t','z','x',x=x)
END
v = g(x)
v = 650
)", kXyzTable);
    INFO(inverse.message);
    REQUIRE(inverse.solveOk);
    CHECK(inverse.verified);
    CHECK_THAT(inverse["x"], WithinAbs(2.5, 1e-9));
}

TEST_CASE("CS-GAP-INTERP-EES: the EES form as the forcing of an INTEGRAL model", "[library-gaps][interp-ees]") {
    // dy/dt = v(t) = 2 t read in the table (as CS-GAP-INTEGRAL-LOOKUP, with the EES call form): y(10) = 100
    const std::map<std::string, std::string> tables = {{"f", "t,v\n0,0\n5,10\n10,20\n"}};
    auto run = runIntegral(
        "y = 0 + integral(dydt, t, 0, 10, 0.5)\n"
        "dydt = interpolate('f', 'v', 't', t=t)\n",
        [](SolverOptions& o) { o.integralMethod = "RK4"; }, tables);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(100.0, 1e-9));
}

TEST_CASE("CS-GAP-INTERP-EES: 2-D tables, INTERPOLATE2 and INTERPOLATE2DM", "[library-gaps][interp-ees]") {
    // the two named arguments are the input columns, the third column of the call is returned
    auto run = runModel(R"(
xv = 1.25
yv = 1.75
z1 = INTERPOLATE2('g','x','y','z',x=xv,y=yv)
z2 = INTERPOLATE2('g','x','y','z',y=yv,x=xv)
z3 = INTERPOLATE2DM('g','x','y','z',x=xv,y=yv)
z4 = INTERPOLATE2('g','z','x','y',x=xv,y=yv)
z5 = INTERPOLATE2('g','x','y','z',xv,yv)
)", kEesGrid);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["z1"], WithinAbs(18.75, 1e-12));   // 11 + 1*0.25 + 10*0.75
    for (const char* n : {"z2", "z3", "z4", "z5"}) {
        INFO(n);
        CHECK(run[n] == run["z1"]);   // the argument order, the alias, the position of the output column: the same value
    }

    // the exact derivative (r is x, then y): z = 11 + (x - 1) + 10 (y - 1)
    const std::string grid = kEesGrid.at("g");
    auto dx = evalWithSlope("v = INTERPOLATE2('t','x','y','z',x=r,y=1.5)\n", grid, 1.25);
    CHECK_THAT(dx.value, WithinAbs(16.25, 1e-12));
    CHECK_THAT(dx.slope, WithinAbs(1.0, 1e-12));
    auto dy = evalWithSlope("v = INTERPOLATE2('t','x','y','z',x=1.5,y=r)\n", grid, 1.75);
    CHECK_THAT(dy.value, WithinAbs(11.5 + 7.5, 1e-12));
    CHECK_THAT(dy.slope, WithinAbs(10.0, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EES: the EES form renders in LaTeX with its named argument", "[library-gaps][interp-ees]") {
    EESParser parser;
    auto parsed = parser.parse("a = INTERPOLATE('T','A','B',B=v_in)\n");
    REQUIRE(parsed.success);
    IR ir = IR::fromAST(parsed.program);
    const std::string latex = coolsolve::expressionToLatex(ir.getEquations()[0].rhs);
    INFO(latex);
    CHECK(latex.find("INTERPOLATE") != std::string::npos);
    CHECK(latex.find("\\mathrm{A}") != std::string::npos);
    CHECK(latex.find("\\mbox{\\ B}=") != std::string::npos);   // the named argument, then its value
    CHECK(latex.find("v_{in}") != std::string::npos);
}

TEST_CASE("CS-GAP-INTERP-EES: the positional form and the other functions are unchanged", "[library-gaps][interp-ees]") {
    auto run = runModel(R"(
p1 = INTERPOLATE('T','A','B',3)
p2 = INTERPOLATE1('T','A','B',1.5)
h = enthalpy(Water, T=20, P=1E5)
)", kEesTable);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["p1"], WithinAbs(50.0, 1e-12));
    CHECK_THAT(run["p2"], WithinAbs(20.0 * 0.5 + 10.0, 1e-12));
    CHECK_THAT(run["h"], WithinAbs(CoolProp::PropsSI("H", "T", 293.15, "P", 1e5, "Water"), 1e-3));
    // the other table functions still refuse a named argument
    for (const char* call : {"LOOKUP('T', 1, 2, B=3)", "NLOOKUPROWS('T', B=3)", "SUMLOOKUP('T', 'A', A=1)"}) {
        auto errors = parseErrors(std::string("y = ") + call + "\n");
        INFO(call);
        REQUIRE_FALSE(errors.empty());
        CHECK(errors[0].find("does not take named arguments") != std::string::npos);
    }
}

TEST_CASE("CS-GAP-INTERP-EES: negative tests give clear errors", "[library-gaps][interp-ees]") {
    auto fail = [&](const std::string& expr, const std::vector<std::string>& expected) {
        auto run = runModel("a = " + expr + "\n", kEesTable);
        INFO(expr);
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        for (const auto& text : expected) CHECK(run.message.find(text) != std::string::npos);
    };
    SECTION("the named argument is not a column of the call") {
        fail("INTERPOLATE('T','A','B',C=3)",
             {"INTERPOLATE(): the named argument 'C' is not one of the columns of the call ('A', 'B') of table 'T'"});
        fail("INTERPOLATE1('T','A','B',a2=3)", {"INTERPOLATE1(): the named argument 'a2' is not one of the columns"});
    }
    SECTION("a column of the table that the call does not list") {
        auto run = runModel("a = INTERPOLATE('T','A','A',B=3)\n", kEesTable);   // B is a column of T, not of this call
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("the named argument 'B' is not one of the columns of the call ('A', 'A')") != std::string::npos);
    }
    SECTION("two named arguments for a 1-D call") {
        fail("INTERPOLATE('T','A','B',A=1,B=2)",
             {"INTERPOLATE(): the EES form takes the table, 2 columns and 1 named argument", "found 3 positional and 2 named"});
    }
    SECTION("a named argument and a fourth positional value") {
        fail("INTERPOLATE('T','A','B',3,B=3)", {"INTERPOLATE(): the EES form takes the table, 2 columns and 1 named argument",
                                               "found 4 positional and 1 named"});
        fail("INTERPOLATE('T','A','B',B=3,3)", {"found 4 positional and 1 named"});
    }
    SECTION("a missing column or value") {
        fail("INTERPOLATE('T','A',B=3)", {"the EES form takes the table, 2 columns", "found 2 positional and 1 named"});
        fail("INTERPOLATE('T',B=3)", {"found 1 positional and 1 named"});
        fail("INTERPOLATE('T','A','B')", {"INTERPOLATE(): expects 4 arguments: (table, xcol, ycol, xval)", "EES form"});
        auto noValue = parseErrors("a = INTERPOLATE('T','A','B',B=)\n");   // `B=` without a value is a parse error
        REQUIRE_FALSE(noValue.empty());
        CHECK(noValue[0].find("Could not parse value for named argument 'B'") != std::string::npos);
    }
    SECTION("unknown table and unknown column keep their messages") {
        auto table = runModel("a = INTERPOLATE('nope','A','B',B=3)\n", kEesTable);
        CHECK_FALSE(table.solveOk);
        CHECK(table.message.find("lookup table 'nope' not found") != std::string::npos);
        fail("INTERPOLATE('T','A','Q',Q=3)", {"INTERPOLATE(): column 'Q' not found in table 'T'"});
        fail("INTERPOLATE('T','Q','B',B=3)", {"column 'Q' not found in table 'T'"});
    }
    SECTION("2-D: the number of named arguments and their names") {
        auto fail2 = [&](const std::string& expr, const std::vector<std::string>& expected) {
            auto run = runModel("a = " + expr + "\n", kEesGrid);
            INFO(expr);
            INFO(run.message);
            CHECK_FALSE(run.solveOk);
            for (const auto& text : expected) CHECK(run.message.find(text) != std::string::npos);
        };
        fail2("INTERPOLATE2('g','x','y','z',x=1)",
              {"INTERPOLATE2(): the EES form takes the table, 3 columns and 2 named arguments", "found 4 positional and 1 named"});
        fail2("INTERPOLATE2('g','x','y','z',x=1,y=1,z=1)", {"found 4 positional and 3 named"});
        fail2("INTERPOLATE2DM('g','x','y','z',1,x=1,y=1)", {"INTERPOLATE2DM(): the EES form", "found 5 positional and 2 named"});
        fail2("INTERPOLATE2('g','x','y','z',x=1,w=1)",
              {"INTERPOLATE2(): the named argument 'w' is not one of the columns of the call ('x', 'y', 'z')"});
        fail2("INTERPOLATE2('g','x','y','z',x=1,X=2)",
              {"INTERPOLATE2(): the named argument 'X' names a column of the call a second time"});
        fail2("INTERPOLATE2('g','x','y','z')", {"INTERPOLATE2(): expects 6 arguments", "EES form"});
    }
    SECTION("in a body the routine is named") {
        auto run = runModel("FUNCTION gf(r)\n  gf = INTERPOLATE('T','A','B',C=r)\nEND\nz = gf(1)\n", kEesTable);
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("the named argument 'C' is not one of the columns of the call") != std::string::npos);
        CHECK(run.message.find("FUNCTION 'gf'") != std::string::npos);
    }
}

// ============================================================================
// CS-GAP-INTERP-EXTRAP: INTERPOLATE extrapolates linearly beyond the table ends
// ============================================================================
// EES continues the line through the two end nodes (evidence: the stored solution of TM-0061, whose
// operating point r_p = 6 lies beyond the last node 5.0 of its 3-node map: epsilon_s = 1.093333 and
// M_r = 29.9, where clamping gives 0.86 and 32.9; the lower end is assumed symmetric). Was: clamped
// to the end values with a zero derivative, silently. 2-D interpolation and LOOKUP are not changed.

namespace {
// x = 4.5, 4.7, 5.0: slope 0.1 over the first interval, 7/30 over the last one
const std::string kExtrapCsv = "x,y\n4.5,0.77\n4.7,0.79\n5.0,0.86\n";
const std::string kExtrapDescCsv = "x,y\n5.0,0.86\n4.7,0.79\n4.5,0.77\n";
constexpr double kExtrapSlopeUp = 0.07 / 0.3;     // last interval
constexpr double kExtrapSlopeLow = 0.02 / 0.2;    // first interval
}  // namespace

TEST_CASE("CS-GAP-INTERP-EXTRAP: reproducer of the register", "[library-gaps][interp-extrap]") {
    auto run = runModel(R"(
x_hi = 6
x_lo = 3
y_hi = INTERPOLATE('t','x','y',x_hi)
y_lo = INTERPOLATE('t','x','y',x_lo)
)", {{"t", kExtrapCsv}});
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["y_hi"], WithinAbs(0.86 + kExtrapSlopeUp, 1e-12));   // 1.093333 (was 0.86)
    CHECK_THAT(run["y_hi"], WithinAbs(1.0933333333, 1e-9));
    CHECK_THAT(run["y_lo"], WithinAbs(0.62, 1e-12));                    // 0.77 + (3 - 4.5) * 0.1 (was 0.77)
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: the EES point of CSL-0035 (r_p = 6 beyond a 3-node map)",
          "[library-gaps][interp-extrap]") {
    const std::map<std::string, std::string> tables = {
        {"lookup_1", "r_p,M_r,epsilon_s\n4.5,34.3,0.77\n4.7,33.8,0.79\n5.0,32.9,0.86\n"}};
    auto run = runModel(R"(
p_su_bar = 1
p_ex_bar = 6*p_su_bar
r_p = p_ex_bar/p_su_bar
epsilon_s = INTERPOLATE('lookup_1','epsilon_s','r_p',r_p=r_p)
M_r = INTERPOLATE('lookup_1','M_r','r_p',r_p=r_p)
r_low = 4
epsilon_low = INTERPOLATE('lookup_1','epsilon_s','r_p',r_p=r_low)
M_low = INTERPOLATE('lookup_1','M_r','r_p',r_p=r_low)
)", tables);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["epsilon_s"], WithinAbs(1.093333, 1e-6));   // EES: 0.86 + (6 - 5) (0.86 - 0.79)/(5.0 - 4.7)
    CHECK_THAT(run["M_r"], WithinAbs(29.9, 1e-9));             // EES: 32.9 - 3.0 (6 - 5)
    CHECK_THAT(run["epsilon_low"], WithinAbs(0.77 - 0.5 * 0.1, 1e-12));   // symmetric below the first node
    CHECK_THAT(run["M_low"], WithinAbs(34.3 + 0.5 * 2.5, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: both call forms, INTERPOLATE and INTERPOLATE1, give the same line",
          "[library-gaps][interp-extrap]") {
    const std::map<std::string, std::string> tables = {{"t", kExtrapCsv}};
    auto run = runModel(R"(
p1 = INTERPOLATE('t','x','y',6)
p2 = INTERPOLATE1('t','x','y',6)
e1 = INTERPOLATE('t','y','x',x=6)
e2 = INTERPOLATE1('t','y','x',x=6)
e3 = INTERPOLATE('t','x','y',x=6)
q1 = INTERPOLATE('t','x','y',3)
q2 = INTERPOLATE1('t','y','x',x=3)
)", tables);
    INFO(run.message);
    REQUIRE(run.solveOk);
    for (const char* n : {"p2", "e1", "e2", "e3"}) {
        INFO(n);
        CHECK(run[n] == run["p1"]);   // one interpolation behind every form, bit for bit
    }
    CHECK_THAT(run["p1"], WithinAbs(0.86 + kExtrapSlopeUp, 1e-12));
    CHECK(run["q2"] == run["q1"]);
    CHECK_THAT(run["q1"], WithinAbs(0.62, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: ascending and descending x columns give the same line",
          "[library-gaps][interp-extrap]") {
    const std::string model = R"(
a = INTERPOLATE('t','x','y',6)
b = INTERPOLATE('t','x','y',3)
c = INTERPOLATE('t','x','y',5)
d = INTERPOLATE('t','x','y',4.5)
e = INTERPOLATE('t','x','y',4.6)
f = INTERPOLATE('t','y','x',y=1.5)
g = INTERPOLATE('t','y','x',y=0.5)
)";
    auto asc = runModel(model, {{"t", kExtrapCsv}});
    auto desc = runModel(model, {{"t", kExtrapDescCsv}});
    REQUIRE(asc.solveOk);
    REQUIRE(desc.solveOk);
    CHECK_THAT(desc["a"], WithinAbs(0.86 + kExtrapSlopeUp, 1e-12));   // above the first row of a descending column
    CHECK_THAT(desc["b"], WithinAbs(0.62, 1e-12));                    // below its last row
    CHECK_THAT(desc["c"], WithinAbs(0.86, 1e-12));                    // the end nodes themselves
    CHECK_THAT(desc["d"], WithinAbs(0.77, 1e-12));
    for (const char* n : {"a", "b", "c", "d", "e", "f", "g"}) {
        INFO(n);
        CHECK_THAT(desc[n], WithinAbs(asc[n], 1e-12));
    }
    // y is increasing with x: x as a function of y outside the table, in the same way
    CHECK_THAT(asc["f"], WithinAbs(5.0 + (1.5 - 0.86) / kExtrapSlopeUp, 1e-9));
    CHECK_THAT(asc["g"], WithinAbs(4.5 + (0.5 - 0.77) / kExtrapSlopeLow, 1e-9));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: values inside the table and at its nodes are unchanged",
          "[library-gaps][interp-extrap]") {
    auto run = runModel(R"(
n1 = INTERPOLATE('t','x','y',4.5)
n2 = INTERPOLATE('t','x','y',4.7)
n3 = INTERPOLATE('t','x','y',5)
m1 = INTERPOLATE('t','x','y',4.6)
m2 = INTERPOLATE('t','x','y',4.85)
)", {{"t", kExtrapCsv}});
    REQUIRE(run.solveOk);
    CHECK(run["n1"] == 0.77);   // the end nodes come out exactly
    CHECK(run["n2"] == 0.79);
    CHECK(run["n3"] == 0.86);
    CHECK_THAT(run["m1"], WithinAbs(0.78, 1e-12));
    CHECK_THAT(run["m2"], WithinAbs(0.825, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: the derivative is the end slope outside the table (exact AD)",
          "[library-gaps][interp-extrap]") {
    // x column of the lookup t: 1 2 3; y: 10 20 30 (slope 10); z: 100 400 900 (slopes 300, 500)
    auto yAbove = evalWithSlope("v = INTERPOLATE('t','x','y',r)\n", kXyzCsv, 5.0);
    CHECK_THAT(yAbove.value, WithinAbs(50.0, 1e-12));
    CHECK_THAT(yAbove.slope, WithinAbs(10.0, 1e-12));
    auto zAbove = evalWithSlope("v = INTERPOLATE('t','x','z',r)\n", kXyzCsv, 5.0);
    CHECK_THAT(zAbove.value, WithinAbs(900.0 + 2.0 * 500.0, 1e-9));   // the last interval, not the whole table
    CHECK_THAT(zAbove.slope, WithinAbs(500.0, 1e-12));
    auto zBelow = evalWithSlope("v = INTERPOLATE('t','x','z',r)\n", kXyzCsv, 0.0);
    CHECK_THAT(zBelow.value, WithinAbs(100.0 - 300.0, 1e-9));
    CHECK_THAT(zBelow.slope, WithinAbs(300.0, 1e-12));
    // the named form, and the chain rule through the value: z at x = 2 r
    auto named = evalWithSlope("v = INTERPOLATE('t','z','x',x=2*r)\n", kXyzCsv, 2.5);
    CHECK_THAT(named.value, WithinAbs(900.0 + 2.0 * 500.0, 1e-9));
    CHECK_THAT(named.slope, WithinAbs(1000.0, 1e-9));
    // at the end nodes: the end slope (it was 0 there)
    auto first = evalWithSlope("v = INTERPOLATE('t','x','z',r)\n", kXyzCsv, 1.0);
    CHECK_THAT(first.value, WithinAbs(100.0, 1e-12));
    CHECK_THAT(first.slope, WithinAbs(300.0, 1e-12));
    auto last = evalWithSlope("v = INTERPOLATE('t','x','z',r)\n", kXyzCsv, 3.0);
    CHECK_THAT(last.value, WithinAbs(900.0, 1e-12));
    CHECK_THAT(last.slope, WithinAbs(500.0, 1e-12));
    // a descending column
    const std::string desc = "x,z\n3,900\n2,400\n1,100\n";
    auto dHigh = evalWithSlope("v = INTERPOLATE('t','x','z',r)\n", desc, 5.0);
    CHECK_THAT(dHigh.value, WithinAbs(1900.0, 1e-9));
    CHECK_THAT(dHigh.slope, WithinAbs(500.0, 1e-12));
    auto dLow = evalWithSlope("v = INTERPOLATE('t','x','z',r)\n", desc, 0.0);
    CHECK_THAT(dLow.value, WithinAbs(-200.0, 1e-9));
    CHECK_THAT(dLow.slope, WithinAbs(300.0, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: an unknown solved through an extrapolated table converges to the linear root",
          "[library-gaps][interp-extrap]") {
    // y = f(x), y + x = c: the root lies beyond the table, so Newton needs the end slope as derivative
    SECTION("above the last node, positional form") {
        auto run = runModel("y = INTERPOLATE('t','x','y',x)\ny + x = 8\n", {{"t", kExtrapCsv}}, "x = 4.9\ny = 0.8\n");
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        // 0.86 + (x - 5) s + x = 8
        const double xRoot = (8.0 - 0.86 + 5.0 * kExtrapSlopeUp) / (1.0 + kExtrapSlopeUp);
        CHECK_THAT(run["x"], WithinAbs(xRoot, 1e-9));
        CHECK_THAT(run["y"], WithinAbs(8.0 - xRoot, 1e-9));
        CHECK(run["x"] > 5.0);   // (clamping gives 7.14)
    }
    SECTION("above the last node, from the default guess") {
        auto run = runModel("y = INTERPOLATE('t','x','y',x)\ny + x = 8\n", {{"t", kExtrapCsv}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinAbs((8.0 - 0.86 + 5.0 * kExtrapSlopeUp) / (1.0 + kExtrapSlopeUp), 1e-9));
    }
    SECTION("below the first node, EES form") {
        auto run = runModel("y = INTERPOLATE('t','y','x',x=x)\ny + x = 2\n", {{"t", kExtrapCsv}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        // 0.77 + (x - 4.5) s + x = 2
        CHECK_THAT(run["x"], WithinAbs((2.0 - 0.77 + 4.5 * kExtrapSlopeLow) / (1.0 + kExtrapSlopeLow), 1e-9));
        CHECK(run["x"] < 4.5);   // (clamping gives 1.23)
    }
    SECTION("descending column") {
        auto run = runModel("y = INTERPOLATE('t','x','y',x)\ny + x = 8\n", {{"t", kExtrapDescCsv}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["x"], WithinAbs((8.0 - 0.86 + 5.0 * kExtrapSlopeUp) / (1.0 + kExtrapSlopeUp), 1e-9));
    }
    SECTION("inverse lookup: the call is a side of an equation") {
        // x where y = 1.5 (beyond the last node: the column x is the output)
        auto run = runModel("INTERPOLATE('t','x','y',y=N) = 5.5\nINTERPOLATE('t','y','x',x=M) = 1.5\n", {{"t", kExtrapCsv}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK(run.verified);
        CHECK_THAT(run["N"], WithinAbs(0.86 + 0.5 * kExtrapSlopeUp, 1e-9));   // y at x = 5.5
        CHECK_THAT(run["M"], WithinAbs(5.0 + (1.5 - 0.86) / kExtrapSlopeUp, 1e-9));
    }
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: tables of one and two rows", "[library-gaps][interp-extrap]") {
    SECTION("two rows: the line through them, on both sides") {
        auto run = runModel(R"(
a = INTERPOLATE('t','x','y',5)
b = INTERPOLATE('t','x','y',-1)
c = INTERPOLATE('t','x','y',2)
d = INTERPOLATE('t','y','x',y=50)
)", {{"t", "x,y\n1,10\n3,30\n"}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(50.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(-10.0, 1e-12));
        CHECK_THAT(run["c"], WithinAbs(20.0, 1e-12));
        CHECK_THAT(run["d"], WithinAbs(5.0, 1e-12));
        auto slope = evalWithSlope("v = INTERPOLATE('t','x','y',r)\n", "x,y\n1,10\n3,30\n", 7.0);
        CHECK_THAT(slope.value, WithinAbs(70.0, 1e-12));
        CHECK_THAT(slope.slope, WithinAbs(10.0, 1e-12));
    }
    SECTION("one row: the single value, zero derivative (unchanged)") {
        auto run = runModel("a = INTERPOLATE('t','x','y',5)\nb = INTERPOLATE('t','x','y',-1)\n", {{"t", "x,y\n1,10\n"}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(10.0, 1e-12));
        CHECK_THAT(run["b"], WithinAbs(10.0, 1e-12));
        auto slope = evalWithSlope("v = INTERPOLATE('t','x','y',r)\n", "x,y\n1,10\n", 7.0);
        CHECK_THAT(slope.value, WithinAbs(10.0, 1e-12));
        CHECK_THAT(slope.slope, WithinAbs(0.0, 1e-12));
    }
    SECTION("a repeated end node: no division by zero, the end value is held") {
        auto run = runModel("a = INTERPOLATE('t','x','y',5)\nb = INTERPOLATE('t','x','y',0)\n",
                            {{"t", "x,y\n1,10\n3,30\n3,30\n"}});
        INFO(run.message);
        REQUIRE(run.solveOk);
        CHECK_THAT(run["a"], WithinAbs(30.0, 1e-12));   // beyond the repeated node: its value, slope 0
        CHECK_THAT(run["b"], WithinAbs(0.0, 1e-12));    // below the first node: the first interval is regular
    }
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: inside FUNCTION and PROCEDURE bodies", "[library-gaps][interp-extrap]") {
    auto run = runModel(R"(
FUNCTION g(r)
  g = INTERPOLATE('t','x','z',r)
END
PROCEDURE pk(r : y, z)
  y = INTERPOLATE1('t','z','x',x=r)
  z = INTERPOLATE('t','x','y',2*r)
END
fv = g(5)
fl = g(0)
CALL pk(5 : yv, zv)
)", kXyzTable);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK(run.verified);
    CHECK_THAT(run["fv"], WithinAbs(1900.0, 1e-9));
    CHECK_THAT(run["fl"], WithinAbs(-200.0, 1e-9));
    CHECK_THAT(run["yv"], WithinAbs(1900.0, 1e-9));
    CHECK_THAT(run["zv"], WithinAbs(100.0, 1e-9));   // y = 10 x at x = 10

    // the output of a routine that reads the table beyond its end is an unknown of the block: Newton with the end slope
    auto inverse = runModel(R"(
FUNCTION g(x)
  g = INTERPOLATE('t','x','z',x)
END
v = g(x)
v = 1900
)", kXyzTable, "x = 2.5\nv = 650\n");
    INFO(inverse.message);
    REQUIRE(inverse.solveOk);
    CHECK(inverse.verified);
    CHECK_THAT(inverse["x"], WithinAbs(5.0, 1e-9));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: the forcing of an INTEGRAL model beyond the table", "[library-gaps][interp-extrap]") {
    // dy/dt = v(t) with the table v = 2 t on t = 0 .. 5, read up to t = 10: the line continues, y(10) = 100
    const std::map<std::string, std::string> tables = {{"f", "t,v\n0,0\n5,10\n"}};
    auto run = runIntegral(
        "y = 0 + integral(dydt, t, 0, 10, 0.5)\n"
        "dydt = interpolate('f', 't', 'v', t)\n",
        [](SolverOptions& o) { o.integralMethod = "RK4"; }, tables);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["y"], WithinAbs(100.0, 1e-9));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: the 2-D functions and LOOKUP are not changed (they clamp)",
          "[library-gaps][interp-extrap]") {
    // The behaviour of EES 2-D extrapolation is not established: INTERPOLATE2 holds the edge values with a zero
    // derivative outside the grid, as before; LOOKUP clamps the row and the column index.
    auto run = runModel(R"(
z_hi = INTERPOLATE2('g','x','y','z',5,5)
z_lo = INTERPOLATE2('g','x','y','z',-3,-3)
z_x = INTERPOLATE2('g','x','y','z',5,1.5)
z_e = INTERPOLATE2('g','x','y','z',x=5,y=1.5)
l_hi = LOOKUP('g', 99, 3)
)", kEesGrid);
    INFO(run.message);
    REQUIRE(run.solveOk);
    CHECK_THAT(run["z_hi"], WithinAbs(22.0, 1e-12));
    CHECK_THAT(run["z_lo"], WithinAbs(11.0, 1e-12));
    CHECK_THAT(run["z_x"], WithinAbs(17.0, 1e-12));   // x held at 2, y = 1.5: (12 + 22)/2
    CHECK(run["z_e"] == run["z_x"]);
    CHECK_THAT(run["l_hi"], WithinAbs(22.0, 1e-12));
    auto slope = evalWithSlope("v = INTERPOLATE2('t','x','y','z',r,1.5)\n", kEesGrid.at("g"), 5.0);
    CHECK_THAT(slope.slope, WithinAbs(0.0, 1e-12));
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: extrapolation is silent (no warning, as in EES)", "[library-gaps][interp-extrap]") {
    auto run = runModel("y = INTERPOLATE('t','x','y',6)\nz = INTERPOLATE('t','x','y',3)\n", {{"t", kExtrapCsv}});
    REQUIRE(run.solveOk);
    CHECK(run.warnings.empty());
}

TEST_CASE("CS-GAP-INTERP-EXTRAP: negative tests, the errors of INTERPOLATE are kept", "[library-gaps][interp-extrap]") {
    auto fail = [&](const std::string& expr, const std::vector<std::string>& expected) {
        auto run = runModel("a = " + expr + "\n", {{"t", kExtrapCsv}});
        INFO(expr);
        INFO(run.message);
        CHECK_FALSE(run.solveOk);
        for (const auto& text : expected) CHECK(run.message.find(text) != std::string::npos);
    };
    SECTION("a table or a column that does not exist, beyond the table ends too") {
        fail("INTERPOLATE('nope','x','y',6)", {"lookup table 'nope' not found"});
        fail("INTERPOLATE('t','x','q',6)", {"INTERPOLATE(): column 'q' not found in table 't'"});
        fail("INTERPOLATE('t','x','y',q=6)", {"INTERPOLATE(): the named argument 'q' is not one of the columns of the call ('x', 'y')"});
    }
    SECTION("a call without a value is still an error") {
        fail("INTERPOLATE('t','x','y')", {"INTERPOLATE(): expects 4 arguments: (table, xcol, ycol, xval)"});
    }
    SECTION("an empty table (header only) is an error, not an extrapolation") {
        auto run = runModel("a = INTERPOLATE('t','x','y',6)\n", {{"t", "x,y\n"}});
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("table is empty") != std::string::npos);
    }
    SECTION("a column index out of the table") {
        auto run = runModel("a = INTERPOLATE('t',1,7,6)\n", {{"t", kExtrapCsv}});
        CHECK_FALSE(run.solveOk);
        CHECK(run.message.find("LookupTable 't': interpolate1D column index out of range") != std::string::npos);
    }
}
