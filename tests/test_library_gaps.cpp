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
#include "coolsolve/ir.h"
#include "coolsolve/lookup_table.h"
#include "coolsolve/parser.h"
#include "coolsolve/runner.h"
#include "coolsolve/solution_checker.h"
#include "CoolProp.h"

#include <filesystem>
#include <fstream>
#include <map>
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
    std::map<std::string, double, CaseInsensitiveLess> vars;

    double operator[](const std::string& name) const { return vars.at(name); }
};

/// Write `code` to a temporary model file, solve it and verify the solution.
/// `tables` maps a lookup-table name to its CSV text (written as `<model>-<name>.csv`).
ModelRun runModel(const std::string& code,
                  const std::map<std::string, std::string>& tables = {}) {
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

    ModelRun run;
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    options.tolerance = 1e-9;
    run.solveOk = runner.run(options);
    run.parseOk = runner.isParseSuccess();
    for (const auto& d : runner.getDiagnostics().items())
        if (d.severity == DiagnosticSeverity::Warning) run.warnings.push_back(d.message);
    if (!run.parseOk && !runner.getParseResult().errors.empty())
        run.message = runner.getParseResult().errors[0].message;
    if (run.solveOk) {
        const auto& res = runner.getSolveResult();
        run.vars = res.variables;
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
    SECTION("ENDIF without IF") {
        auto errs = parseErrors("FUNCTION f(v)\n f = 1\n ENDIF\nEND\na = f(1)\n");
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
    // EES: 282 989.9 kJ/kmol (stored Q_4 of the cpbar test file of the ULiege library)
    CHECK_THAT(run["LHV_CO_mol"], WithinRel(282.990e6, 1e-5));
    CHECK_THAT(run["LHV_H2_mol"], WithinRel(241.820e6, 1e-5));
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
u_CO2 = intenergy(CO2,T=25,P=101325)
)");
    REQUIRE(run.solveOk);
    // J/kmol
    CHECK_THAT(run["hf_CO2"], WithinRel(-393.52e6, 1e-9));
    CHECK_THAT(run["hf_CO"], WithinRel(-110.53e6, 1e-9));
    CHECK_THAT(run["hf_H2O"], WithinRel(-241.82e6, 1e-9));
    CHECK_THAT(run["hf_CH4"], WithinRel(-74.85e6, 1e-9));
    CHECK_THAT(run["hf_C2H6"], WithinRel(-84.68e6, 1e-9));
    CHECK_THAT(run["hf_C3H8"], WithinRel(-103.85e6, 1e-9));
    // elements in their standard state
    CHECK_THAT(run["h_N2"], WithinAbs(0.0, 1e-3));
    CHECK_THAT(run["h_O2"], WithinAbs(0.0, 1e-3));
    CHECK_THAT(run["h_H2"], WithinAbs(0.0, 1e-3));
    // the internal energy follows the same reference (u = h - R T)
    CHECK(run["u_CO2"] < run["hf_CO2"] / molarMassKgPerKmol("CarbonDioxide"));
}

TEST_CASE("CS-GAP-FORMATION-ENTHALPY: enthalpy differences, inversion and entropy are untouched", "[library-gaps][formation-enthalpy]") {
    // Only the reference of the enthalpy changes: enthalpy differences are CoolProp's, an enthalpy
    // given as input is shifted back, and entropies keep CoolProp's reference state.
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
    CHECK_THAT(run["s_CO2"], WithinRel(CoolProp::PropsSI("S", "T", 298.15, "P", P, "CarbonDioxide"), 1e-9));
    // species without a formation enthalpy keep CoolProp's reference state
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
    // outside the table: flat extrapolation at the nearest end, whatever the order
    CHECK_THAT(run["M_below"], WithinAbs(370.0, 1e-9));
    CHECK_THAT(run["M_above"], WithinAbs(454.0, 1e-9));
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
