/**
 * @file test_latex_report.cpp
 * @brief Tests for the LaTeX report generator.
 *
 * Validates:
 *   - generateLatexReport() produces valid LaTeX for a simple model
 *   - writeLatexReport() creates a .tex file on disk
 *   - Config loading of enableLatexReport / latexCompiler
 *   - Runner integration (generateLatexReportContent + cache)
 *
 * Run with: ./coolsolve_tests "[latex]"
 */

#include <catch2/catch_test_macros.hpp>
#include "coolsolve/latex_report.h"
#include "coolsolve/ir.h"
#include "coolsolve/runner.h"
#include "coolsolve/solver.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;
using namespace coolsolve;

// ============================================================================
// Helpers
// ============================================================================

static fs::path getExamplesDir() {
    const char* env = std::getenv("COOLSOLVE_EXAMPLES_DIR");
    if (env && env[0] != '\0') return fs::path(env);
    return fs::path("../examples");
}

static std::string readFile(const fs::path& p) {
    std::ifstream f(p);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// ============================================================================
// Config tests
// ============================================================================

TEST_CASE("Config: enableLatexReport defaults to true", "[latex][config]") {
    SolverOptions opts;
    REQUIRE(opts.enableLatexReport == true);
    REQUIRE(opts.latexCompiler == "pdflatex");
}

TEST_CASE("Config: enableLatexReport and latexCompiler are loaded", "[latex][config]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_latex.conf";
    {
        std::ofstream f(configPath);
        REQUIRE(f.is_open());
        f << "enableLatexReport = true\n";
        f << "latexCompiler = xelatex\n";
    }
    SolverOptions opts;
    bool loaded = loadSolverOptionsFromFile(configPath.string(), opts);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(opts.enableLatexReport == true);
    REQUIRE(opts.latexCompiler == "xelatex");
}

// ============================================================================
// Report generation with a real model (cpbar — small and fast)
// ============================================================================

TEST_CASE("LaTeX report from cpbar model", "[latex]") {
    fs::path exDir = getExamplesDir();
    fs::path eesPath = exDir / "cpbar.eescode";
    if (!fs::exists(eesPath)) {
        SKIP("cpbar.eescode not found in examples directory");
    }

    CoolSolveRunner runner(eesPath.string());
    SolverOptions opts;
    opts.verbose = false;
    bool ok = runner.run(opts, /*enableTracing=*/false);
    REQUIRE(ok);

    SECTION("generateLatexReport produces a complete document") {
        std::string tex = generateLatexReport(
            runner.getIR(),
            runner.getAnalysisResult(),
            runner.getSolveResult(),
            runner.getTiming(),
            "cpbar",
            eesPath.string()
        );
        // Non-empty
        REQUIRE(!tex.empty());
        // Starts with documentclass
        REQUIRE(tex.find("\\documentclass") != std::string::npos);
        // Contains expected sections
        REQUIRE(tex.find("\\section{Model Overview}") != std::string::npos);
        REQUIRE(tex.find("\\section{Equations}") != std::string::npos);
        REQUIRE(tex.find("\\section{Variable Solutions}") != std::string::npos);
        REQUIRE(tex.find("\\section{Plots}") != std::string::npos);
        // Removed sections should NOT be present
        REQUIRE(tex.find("\\section{Block Structure}") == std::string::npos);
        REQUIRE(tex.find("\\section{Solver Performance}") == std::string::npos);
        REQUIRE(tex.find("\\section{Timing Profile}") == std::string::npos);
        // Equations section should not contain block subsections
        REQUIRE(tex.find("\\subsection{Block") == std::string::npos);
        // Variable table should not contain a # column header
        REQUIRE(tex.find("\\textbf{\\#}") == std::string::npos);
        // Contains end document
        REQUIRE(tex.find("\\end{document}") != std::string::npos);
        // Contains the model name
        REQUIRE(tex.find("cpbar") != std::string::npos);
    }

    SECTION("writeLatexReport writes to disk") {
        fs::path tmpTex = fs::temp_directory_path() / "coolsolve_test_report.tex";
        bool written = writeLatexReport(
            tmpTex.string(),
            runner.getIR(),
            runner.getAnalysisResult(),
            runner.getSolveResult(),
            runner.getTiming(),
            "cpbar",
            eesPath.string()
        );
        REQUIRE(written);
        REQUIRE(fs::exists(tmpTex));
        std::string content = readFile(tmpTex);
        REQUIRE(content.find("\\documentclass") != std::string::npos);
        fs::remove(tmpTex);
    }
}

// ============================================================================
// Runner integration
// ============================================================================

TEST_CASE("Runner: generateLatexReportContent caches report", "[latex]") {
    fs::path exDir = getExamplesDir();
    fs::path eesPath = exDir / "cpbar.eescode";
    if (!fs::exists(eesPath)) {
        SKIP("cpbar.eescode not found in examples directory");
    }

    CoolSolveRunner runner(eesPath.string());
    SolverOptions opts;
    opts.verbose = false;
    bool ok = runner.run(opts, false);
    REQUIRE(ok);

    // Initially no report cached
    REQUIRE_FALSE(runner.hasLatexReport());

    // Generate and cache
    std::string tex = runner.generateLatexReportContent("cpbar");
    REQUIRE(!tex.empty());
    REQUIRE(runner.hasLatexReport());
    REQUIRE(runner.getLatexReportContent() == tex);
}

// ============================================================================
// EES-to-LaTeX variable name translation
// ============================================================================

TEST_CASE("variableToLatex: EES name translation", "[latex]") {
    SECTION("Simple name") {
        CHECK(variableToLatex("T") == "T");
        CHECK(variableToLatex("P") == "P");
    }
    SECTION("Underscore subscripts") {
        CHECK(variableToLatex("T_oil_su") == "T_{oil,su}");
        CHECK(variableToLatex("h_liq") == "h_{liq}");
    }
    SECTION("Array brackets") {
        CHECK(variableToLatex("P[1]") == "P_{1}");
        CHECK(variableToLatex("T_su[2]") == "T_{su,2}");
        CHECK(variableToLatex("T_hf[3]") == "T_{hf,3}");
        CHECK(variableToLatex("s[10]") == "s_{10}");
    }
    SECTION("Greek letters (whole name)") {
        CHECK(variableToLatex("eta") == "\\eta");
        CHECK(variableToLatex("omega") == "\\omega");
        CHECK(variableToLatex("alpha") == "\\alpha");
    }
    SECTION("Greek letters as first segment with subscript") {
        CHECK(variableToLatex("eta_s") == "\\eta_{s}");
        CHECK(variableToLatex("eta_boil_1") == "\\eta_{boil,1}");
    }
    SECTION("_dot modifier") {
        CHECK(variableToLatex("M_dot") == "\\dot{M}");
        CHECK(variableToLatex("Q_dot_evap") == "\\dot{Q}_{evap}");
    }
    SECTION("_bar modifier") {
        CHECK(variableToLatex("h_bar") == "\\bar{h}");
        CHECK(variableToLatex("eta_bar_boil_1") == "\\bar{\\eta}_{boil,1}");
    }
    SECTION("Uppercase Greek prefix (DELTA)") {
        CHECK(variableToLatex("DELTAP") == "\\Delta P");
        CHECK(variableToLatex("DELTAW_dot_cp") == "\\Delta \\dot{W}_{cp}");
        CHECK(variableToLatex("DELTAT_lm") == "\\Delta T_{lm}");
    }
    SECTION("No false positive on partial lowercase Greek") {
        // "pipe" should NOT match the Greek "pi"
        CHECK(variableToLatex("pipe_dia") == "pipe_{dia}");
    }
    SECTION("String variable $ suffix stripped") {
        CHECK(variableToLatex("fluid$") == "fluid");
        CHECK(variableToLatex("R$") == "R");
    }
}

// ============================================================================
// Pressuredrop model (string vars, $ in names, "" comments)
// ============================================================================

TEST_CASE("LaTeX report from pressuredrop model", "[latex]") {
    fs::path exDir = getExamplesDir();
    fs::path eesPath = exDir / "pressuredrop.eescode";
    if (!fs::exists(eesPath)) {
        SKIP("pressuredrop.eescode not found in examples directory");
    }

    CoolSolveRunner runner(eesPath.string());
    SolverOptions opts;
    opts.verbose = false;
    bool ok = runner.run(opts, false);
    REQUIRE(ok);

    std::string tex = runner.generateLatexReportContent("pressuredrop");
    REQUIRE(!tex.empty());

    // No bare $ inside equation environments (would break pdflatex)
    // Find each \begin{equation} ... \end{equation} and ensure no bare $
    size_t pos = 0;
    while ((pos = tex.find("\\begin{equation}", pos)) != std::string::npos) {
        size_t eqEnd = tex.find("\\end{equation}", pos);
        REQUIRE(eqEnd != std::string::npos);
        std::string eqBody = tex.substr(pos + 16, eqEnd - pos - 16);
        // A bare $ would break math mode; escaped \$ is fine
        for (size_t i = 0; i < eqBody.size(); ++i) {
            if (eqBody[i] == '$') {
                // Check it's not escaped
                bool escaped = (i > 0 && eqBody[i-1] == '\\');
                REQUIRE(escaped);
            }
        }
        pos = eqEnd + 14;
    }

    // String assignment equations (fluid$='...') should NOT appear in equation blocks
    // (but 'fluid' may still appear in the variable solutions table)
    {
        size_t eqPos = 0;
        while ((eqPos = tex.find("\\begin{equation}", eqPos)) != std::string::npos) {
            size_t eqEnd = tex.find("\\end{equation}", eqPos);
            REQUIRE(eqEnd != std::string::npos);
            std::string eqBody = tex.substr(eqPos, eqEnd - eqPos);
            // The string assignment should not be an equation
            REQUIRE(eqBody.find("n-pentane") == std::string::npos);
            eqPos = eqEnd + 14;
        }
    }

    // "" display comments should appear as bold text
    REQUIRE(tex.find("\\textbf{") != std::string::npos);

    // Document should be structurally complete
    REQUIRE(tex.find("\\begin{document}") != std::string::npos);
    REQUIRE(tex.find("\\end{document}") != std::string::npos);
}

// ============================================================================
// Comments in LaTeX report
// ============================================================================

TEST_CASE("LaTeX report includes display comments", "[latex]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path tmpEes = tmpDir / "latex_comments_test.eescode";
    {
        std::ofstream f(tmpEes);
        REQUIRE(f.is_open());
        f << "\"Section A\"\n";
        f << "x + y = 10\n";
        f << "x - y = 2 \"difference equation\"\n";
        f << "{this is a hidden comment}\n";
        f << "\"Section B should not appear after hidden\"\n";
    }

    CoolSolveRunner runner(tmpEes.string());
    SolverOptions opts;
    opts.verbose = false;
    bool ok = runner.run(opts, false);
    fs::remove(tmpEes);
    REQUIRE(ok);

    std::string tex = runner.generateLatexReportContent("latex_comments_test");

    // "" standalone comment "Section A" should appear as bold text
    REQUIRE(tex.find("\\textbf{Section A}") != std::string::npos);

    // "" inline comment should appear as italic
    REQUIRE(tex.find("\\textit{difference equation}") != std::string::npos);

    // {} hidden comment should NOT appear in the report text
    REQUIRE(tex.find("hidden comment") == std::string::npos);
}

// ============================================================================
// LaTeX escaping edge cases
// ============================================================================

TEST_CASE("LaTeX report handles model with special characters", "[latex]") {
    // Create a temporary .eescode file with underscores and other chars in
    // variable names — these must be escaped correctly in the LaTeX output.
    fs::path tmpDir = fs::temp_directory_path();
    fs::path tmpEes = tmpDir / "latex_escape_test.eescode";
    {
        std::ofstream f(tmpEes);
        REQUIRE(f.is_open());
        // Simple 2-equation system: x_val + y_val = 3, x_val - y_val = 1
        f << "x_val + y_val = 3\n";
        f << "x_val - y_val = 1\n";
    }

    CoolSolveRunner runner(tmpEes.string());
    SolverOptions opts;
    opts.verbose = false;
    bool ok = runner.run(opts, false);
    fs::remove(tmpEes);
    REQUIRE(ok);

    std::string tex = runner.generateLatexReportContent("latex_escape_test");
    REQUIRE(!tex.empty());

    // The underscore in variable names should be escaped or rendered in math mode
    // (not raw _ which would break LaTeX compilation)
    // Check that the document is structurally complete
    REQUIRE(tex.find("\\begin{document}") != std::string::npos);
    REQUIRE(tex.find("\\end{document}") != std::string::npos);
}
