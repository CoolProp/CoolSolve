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

#include "coolsolve/parser.h"
#include "coolsolve/runner.h"
#include "coolsolve/solution_checker.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <string>

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
    std::map<std::string, double, CaseInsensitiveLess> vars;

    double operator[](const std::string& name) const { return vars.at(name); }
};

/// Write `code` to a temporary model file, solve it and verify the solution.
ModelRun runModel(const std::string& code) {
    static int counter = 0;
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_library_gaps";
    fs::create_directories(dir);
    const fs::path file = dir / ("model_" + std::to_string(counter++) + ".eescode");
    {
        std::ofstream f(file);
        f << code;
    }

    ModelRun run;
    CoolSolveRunner runner(file.string());
    SolverOptions options;
    options.tolerance = 1e-9;
    run.solveOk = runner.run(options);
    run.parseOk = runner.isParseSuccess();
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
