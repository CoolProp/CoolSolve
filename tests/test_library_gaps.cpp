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
    // an index range is not an expression (only the $IntegralTable directive reads x[1..3])
    CHECK(parseErrors("y = sum(x[1..3])\n") == Messages{"Could not parse index 1 of 'x': '1..3'"});
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
    // an output must be a variable
    CHECK(anyContains(parseErrors(proc + "CALL p(1, 2 : 3)\n"), "Output 1 of 'p' must be a variable, not the expression '3'"));
    CHECK(anyContains(parseErrors(proc + "CALL p(1, 2 : a+b)\n"), "Output 1 of 'p' must be a variable, not the expression 'a+b'"));
    // (register CS-GAP-CALL-EXPR-OUT: EES accepts CALL BRINEPROP(... : c_p/1000); it was dropped silently)
    CHECK(anyContains(parseErrors(proc + "CALL p(1, 2 : r, z/1000)\n"), "Output 2 of 'p' must be a variable"));
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
    // refused by the parser (their evaluation is CS-GAP-SUM-INDEXED / CS-GAP-INTERP-EES)
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
$COMMON a

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
