/**
 * Tests for `$INCLUDE` (CS-FEAT-IMPORT, CS-GAP-INCLUDE): the resolution of a name (a file of the model
 * folder first, then a library model by name or ID), the import of the FUNCTIONs and PROCEDUREs only,
 * nested includes, the precedence rules, and the errors (unknown name, path, EES file, cycle, error in
 * the included file). Uses an in-memory library, so it does not depend on the embedded snapshot.
 *
 * Run with: ./coolsolve_tests "[include]"
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "coolsolve/include_resolver.h"
#include "coolsolve/latex_report.h"
#include "coolsolve/library_store.h"
#include "coolsolve/parser.h"
#include "coolsolve/runner.h"

#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;
using namespace coolsolve;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinRel;

namespace {

const char* kIndex = R"({
  "snapshot": {"content_hash": "sha256:abc", "models": 3},
  "models": [
    {"id": "CSL-9001", "name": "lib_funcs", "title": "Library functions", "category": "fundamentals",
     "status": "verified", "level": 1, "kind": "function", "path": "models/fundamentals/lib_funcs",
     "main_file": "lib_funcs.eescode"},
    {"id": "CSL-9002", "name": "lib_nested", "title": "Nested functions", "category": "fundamentals",
     "status": "verified", "level": 1, "kind": "function", "path": "models/fundamentals/lib_nested",
     "main_file": "lib_nested.eescode"},
    {"id": "CSL-9003", "name": "lib_other_twice", "title": "Another twice", "category": "fundamentals",
     "status": "verified", "level": 1, "kind": "function", "path": "models/fundamentals/lib_other_twice",
     "main_file": "lib_other_twice.eescode"}
  ],
  "functions": [
    {"function": "twice", "type": "FUNCTION", "model_id": "CSL-9001"},
    {"function": "split", "type": "PROCEDURE", "model_id": "CSL-9001"},
    {"function": "quad", "type": "FUNCTION", "model_id": "CSL-9002"}
  ]
})";

// A function model: definitions, then a demonstration program that an $INCLUDE ignores
const char* kLibFuncs = R"(FUNCTION twice(x)
  twice := 2*x
END

PROCEDURE split(x : a, b)
  a := x/4
  b := 3*x/4
END

"Demonstration program"
y_demo = twice(3)
)";

const char* kLibNested = R"($INCLUDE lib_funcs
FUNCTION quad(x)
  quad := twice(twice(x))
END
z_demo = quad(1)
)";

const char* kLibOtherTwice = R"(FUNCTION twice(x)
  twice := 20*x
END
)";

const LibraryStore& testLibrary() {
    static const LibraryStore store(LibraryFileSource::fromMap({
        {"library.json", kIndex},
        {"models/fundamentals/lib_funcs/lib_funcs.eescode", kLibFuncs},
        {"models/fundamentals/lib_nested/lib_nested.eescode", kLibNested},
        {"models/fundamentals/lib_other_twice/lib_other_twice.eescode", kLibOtherTwice},
    }));
    return store;
}

ParseResult parseWith(const std::string& source, std::map<std::string, std::string> projectFiles = {},
                      const LibraryStore* library = &testLibrary()) {
    EESParser parser;
    parser.setIncludeResolver(IncludeResolver::forFiles(std::move(projectFiles), "main.eescode", library));
    return parser.parse(source);
}

std::string allErrors(const ParseResult& r) {
    std::string out;
    for (const auto& e : r.errors) out += "line " + std::to_string(e.line) + ": " + e.message + "\n";
    return out;
}

bool hasDiagnostic(const ParseResult& r, const std::string& code, const std::string& text) {
    for (const auto& d : r.diagnostics.items())
        if (d.code == code && d.message.find(text) != std::string::npos) return true;
    return false;
}

const FunctionDefinition* findFunction(const ParseResult& r, const std::string& name) {
    for (const auto& s : r.program.statements)
        if (s->is<FunctionDefinition>() && s->as<FunctionDefinition>().name == name) return &s->as<FunctionDefinition>();
    return nullptr;
}

int countFunctions(const ParseResult& r, const std::string& name) {
    int n = 0;
    for (const auto& s : r.program.statements)
        if (s->is<FunctionDefinition>() && s->as<FunctionDefinition>().name == name) ++n;
    return n;
}

/// A fresh folder with the given files, removed at the end of the test
struct TempProject {
    fs::path dir;
    explicit TempProject(const std::map<std::string, std::string>& files) {
        dir = fs::temp_directory_path() / ("coolsolve_include_test_" + std::to_string(std::rand()));
        fs::remove_all(dir);
        fs::create_directories(dir);
        for (const auto& [name, content] : files) std::ofstream(dir / name) << content;
    }
    ~TempProject() {
        std::error_code ec;
        fs::remove_all(dir, ec);
    }
};

double solveValue(const fs::path& model, const std::string& variable, const LibraryStore* library = &testLibrary()) {
    CoolSolveRunner runner(model.string());
    runner.setIncludeResolver(IncludeResolver::forDirectory(model.parent_path().string(), model.filename().string(),
                                                           library));
    SolverOptions options;
    options.verbose = false;
    REQUIRE(runner.run(options));
    return runner.getSolveResult().variables.at(variable);
}

}  // namespace

TEST_CASE("Include name: comments, quotes and the .eescode extension are removed", "[include]") {
    std::string name, error;
    REQUIRE(parseIncludeName("hx_ees_function", name, error));
    CHECK(name == "hx_ees_function");
    REQUIRE(parseIncludeName("  'my_props.eescode'  \"the properties\"", name, error));
    CHECK(name == "my_props");
    REQUIRE(parseIncludeName("CSL-0331 {by ID} // comment", name, error));
    CHECK(name == "CSL-0331");
}

TEST_CASE("Include name: paths, EES files and an empty name are refused", "[include]") {
    std::string name, error;
    CHECK_FALSE(parseIncludeName("", name, error));
    CHECK_THAT(error, ContainsSubstring("needs the name"));
    CHECK_FALSE(parseIncludeName("sub/props", name, error));
    CHECK_THAT(error, ContainsSubstring("without path"));
    CHECK_FALSE(parseIncludeName("C:\\EES32\\USERLIB\\props.lib", name, error));
    CHECK_FALSE(parseIncludeName("../props", name, error));
    CHECK_FALSE(parseIncludeName("library:props", name, error));
    CHECK_FALSE(parseIncludeName("props.lib", name, error));
    CHECK_THAT(error, ContainsSubstring("convert the EES file into props.eescode"));
    CHECK_FALSE(parseIncludeName("props.txt /R", name, error));   // the /R flag of EES: a path to the parser
}

TEST_CASE("Include a library model: its functions and procedures, not its main program", "[include]") {
    auto r = parseWith("$INCLUDE lib_funcs\ny = twice(4)\nCALL split(8 : a, b)\n");
    REQUIRE(r.success);
    const auto* twice = findFunction(r, "twice");
    REQUIRE(twice != nullptr);
    CHECK(twice->includedFrom == "library model lib_funcs (CSL-9001)");
    // The demonstration program of the library model is not imported
    for (const auto& s : r.program.statements) {
        if (!s->is<Equation>()) continue;
        CHECK(s->as<Equation>().lhs->as<Variable>().name != "y_demo");
    }
    REQUIRE(r.includes.size() == 1);
    CHECK(r.includes[0].line == 1);
    CHECK(r.includes[0].includedBy.empty());
    CHECK(r.includes[0].source.fromLibrary);
    CHECK(r.includes[0].source.libraryId == "CSL-9001");
    CHECK(r.includes[0].routines == std::vector<std::string>{"FUNCTION twice", "PROCEDURE split"});
    // The imported routines are known: no "Unknown function" warning
    CHECK_FALSE(hasDiagnostic(r, "P003", "twice"));
}

TEST_CASE("Include: a library model by ID, in any case", "[include]") {
    REQUIRE(parseWith("$INCLUDE CSL-9001\ny = twice(4)\n").success);
    auto r = parseWith("$include LIB_FUNCS\ny = twice(4)\n");
    REQUIRE(r.success);
    REQUIRE(findFunction(r, "twice") != nullptr);
}

TEST_CASE("Include: a file of the model folder takes precedence over the library model", "[include]") {
    auto r = parseWith("$INCLUDE lib_funcs\ny = twice(4)\n",
                       {{"lib_funcs.eescode", "FUNCTION twice(x)\n  twice := 3*x\nEND\n"}});
    REQUIRE(r.success);
    REQUIRE(r.includes.size() == 1);
    CHECK_FALSE(r.includes[0].source.fromLibrary);
    CHECK(findFunction(r, "twice")->includedFrom == "lib_funcs.eescode");
}

TEST_CASE("Include: a routine of the model replaces an included one (warning)", "[include]") {
    auto r = parseWith("$INCLUDE lib_funcs\nFUNCTION twice(x)\n  twice := 5*x\nEND\ny = twice(4)\n");
    REQUIRE(r.success);
    CHECK(countFunctions(r, "twice") == 1);
    CHECK(findFunction(r, "twice")->includedFrom.empty());
    CHECK(hasDiagnostic(r, "P011", "FUNCTION twice of the model is used"));
}

TEST_CASE("Include: nested includes are imported, a library model resolves in the library only", "[include]") {
    // The project has its own lib_funcs.eescode, but the library model lib_nested includes the library's
    auto r = parseWith("$INCLUDE lib_nested\nz = quad(1)\n",
                       {{"lib_funcs.eescode", "FUNCTION twice(x)\n  twice := 3*x\nEND\n"}});
    REQUIRE(r.success);
    REQUIRE(findFunction(r, "quad") != nullptr);
    REQUIRE(findFunction(r, "twice") != nullptr);
    CHECK(findFunction(r, "twice")->includedFrom == "library model lib_funcs (CSL-9001)");
    REQUIRE(r.includes.size() == 2);
    CHECK(r.includes[0].source.name == "lib_nested");
    CHECK(r.includes[1].source.name == "lib_funcs");
    CHECK(r.includes[1].includedBy == "library model lib_nested (CSL-9002)");
}

TEST_CASE("Include: a routine of a project file wins over the same routine of a library model", "[include]") {
    auto r = parseWith("$INCLUDE lib_nested\n$INCLUDE my_twice\nz = quad(1)\n",
                       {{"my_twice.eescode", "FUNCTION twice(x)\n  twice := 3*x\nEND\n"}});
    REQUIRE(r.success);
    CHECK(countFunctions(r, "twice") == 1);
    CHECK(findFunction(r, "twice")->includedFrom == "my_twice.eescode");
    CHECK(hasDiagnostic(r, "P011", "takes precedence over the library"));
}

TEST_CASE("Include: two library models with the same routine, the first $INCLUDE provides it", "[include]") {
    auto r = parseWith("$INCLUDE lib_funcs\n$INCLUDE lib_other_twice\ny = twice(1)\n");
    REQUIRE(r.success);
    CHECK(countFunctions(r, "twice") == 1);
    CHECK(findFunction(r, "twice")->includedFrom == "library model lib_funcs (CSL-9001)");
    CHECK(hasDiagnostic(r, "P011", "the first $INCLUDE that defines a name provides it"));
}

TEST_CASE("Include: the same file reached twice is imported once", "[include]") {
    // Directly twice: a warning; through a nested include (lib_nested includes lib_funcs): silently
    auto twiceDirect = parseWith("$INCLUDE lib_funcs\n$INCLUDE lib_funcs\ny = twice(1)\n");
    REQUIRE(twiceDirect.success);
    CHECK(countFunctions(twiceDirect, "twice") == 1);
    CHECK(hasDiagnostic(twiceDirect, "P011", "already included above"));

    auto diamond = parseWith("$INCLUDE lib_funcs\n$INCLUDE lib_nested\ny = quad(1)\n");
    REQUIRE(diamond.success);
    CHECK(countFunctions(diamond, "twice") == 1);
    CHECK_FALSE(hasDiagnostic(diamond, "P011", "twice"));
}

TEST_CASE("Include: unknown name, with the closest library names", "[include]") {
    auto r = parseWith("$INCLUDE lib_func\ny = 1\n");
    REQUIRE_FALSE(r.success);
    REQUIRE(r.errors.size() == 1);
    CHECK(r.errors[0].line == 1);
    CHECK_THAT(r.errors[0].message, ContainsSubstring("no file lib_func.eescode in the folder of the model"));
    CHECK_THAT(r.errors[0].message, ContainsSubstring("did you mean lib_funcs"));
}

TEST_CASE("Include: a build without library says so", "[include]") {
    auto r = parseWith("$INCLUDE hx_ees_function\ny = 1\n", {}, nullptr);
    REQUIRE_FALSE(r.success);
    CHECK_THAT(r.errors[0].message, ContainsSubstring("this build does not embed the CoolSolve Library"));
}

TEST_CASE("Include: a cycle and the model itself are errors", "[include]") {
    auto cycle = parseWith("$INCLUDE a\ny = 1\n", {{"a.eescode", "$INCLUDE b\n"}, {"b.eescode", "$INCLUDE a\n"}});
    REQUIRE_FALSE(cycle.success);
    CHECK_THAT(allErrors(cycle), ContainsSubstring("Circular $INCLUDE: a.eescode → b.eescode → a.eescode"));

    auto self = parseWith("$INCLUDE main\ny = 1\n", {{"main.eescode", "y = 1"}});
    REQUIRE_FALSE(self.success);
    CHECK_THAT(allErrors(self), ContainsSubstring("main.eescode is the model itself"));
}

TEST_CASE("Include: an error of the included file is reported at the directive, with its file and line",
          "[include]") {
    auto r = parseWith("y = 1\n$INCLUDE broken\n", {{"broken.eescode", "FUNCTION f(x)\n  f := x\n"}});
    REQUIRE_FALSE(r.success);
    REQUIRE(r.errors.size() == 1);
    CHECK(r.errors[0].line == 2);
    CHECK_THAT(r.errors[0].message, ContainsSubstring("In broken.eescode, line 1: Function 'f' missing END"));
}

TEST_CASE("Include: only in the main program", "[include]") {
    auto inBody = parseWith("FUNCTION f(x)\n$INCLUDE lib_funcs\n  f := x\nEND\ny = f(1)\n");
    REQUIRE_FALSE(inBody.success);
    CHECK_THAT(allErrors(inBody), ContainsSubstring("$INCLUDE belongs to the main program"));
}

TEST_CASE("Include: in the removed branch of an $IF, the directive is ignored", "[include]") {
    auto r = parseWith("act$ = 'a'\n$IF act$ = 'b'\n$INCLUDE does_not_exist\n$ENDIF\ny = 1\n");
    REQUIRE(r.success);
    CHECK(r.includes.empty());
}

TEST_CASE("Include: the $CONSTANTs and ranges of the included file are resolved in it", "[include]") {
    auto r = parseWith("$INCLUDE consts\ny = scaled(2)\n",
                       {{"consts.eescode", "$CONSTANT k# = 7\nFUNCTION scaled(x)\n  scaled := k#*x\nEND\n"}});
    REQUIRE(r.success);
    REQUIRE(findFunction(r, "scaled") != nullptr);
}

TEST_CASE("Include: solve a model with a project file and a library model (CLI folder)", "[include]") {
    TempProject project({
        {"main.eescode", "$INCLUDE helpers\n$INCLUDE lib_funcs\ny = triple(twice(2))\nCALL split(8 : a, b)\n"},
        {"helpers.eescode", "FUNCTION triple(x)\n  triple := 3*x\nEND\nunused = 1\n"},
    });
    const fs::path model = project.dir / "main.eescode";
    CHECK_THAT(solveValue(model, "y"), WithinRel(12.0, 1e-12));
    CHECK_THAT(solveValue(model, "b"), WithinRel(6.0, 1e-12));

    // Without a resolver, parseFile() looks in the folder of the file (whatever the working directory)
    EESParser parser;
    auto r = parser.parseFile(model.string());
    REQUIRE(r.includes.size() >= 1);
    CHECK(r.includes[0].source.fileName == "helpers.eescode");
}

TEST_CASE("Include: the debug report lists the included files; the LaTeX report names, not lists, the routines",
          "[include][latex]") {
    TempProject project({
        {"main.eescode", "$INCLUDE helpers\ny = triple(2)\n"},
        {"helpers.eescode", "FUNCTION triple(x)\n  triple := 3*x + 0.123456789*0\nEND\n"},
    });
    const fs::path model = project.dir / "main.eescode";
    CoolSolveRunner runner(model.string());
    SolverOptions options;
    options.verbose = false;
    REQUIRE(runner.run(options));

    const std::string tex = runner.generateLatexReportContent("main");
    CHECK_THAT(tex, ContainsSubstring("Included routines"));
    CHECK_THAT(tex, ContainsSubstring("helpers.eescode"));
    CHECK_THAT(tex, ContainsSubstring("\\texttt{triple}"));
    CHECK_THAT(tex, !ContainsSubstring("0.123456789"));   // the body of the routine is not in the report

    const fs::path debugDir = project.dir / "debug";
    runner.generateDebugOutput(debugDir.string(), "$INCLUDE helpers\ny = triple(2)\n");
    std::ifstream in(debugDir / "report.md");
    std::stringstream report;
    report << in.rdbuf();
    CHECK_THAT(report.str(), ContainsSubstring("## Included files"));
    CHECK_THAT(report.str(), ContainsSubstring("FUNCTION triple"));
}

TEST_CASE("Include: a model without $INCLUDE is parsed as before", "[include]") {
    auto r = parseWith("FUNCTION f(x)\n  f := x\nEND\ny = f(1)\n");
    REQUIRE(r.success);
    CHECK(r.includes.empty());
    CHECK(findFunction(r, "f")->includedFrom.empty());
}
