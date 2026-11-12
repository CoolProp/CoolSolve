// Tests of the outputs of the Diagrams module (docs/coolplot_integration.md, phase 4): the badge lines
// (src/diagram/references.cpp) against CoolPlot's text, the SVG and PDF files of a chart, the diagram of a
// definitions document written by `coolsolve --diagram` (writeDiagramFile) and the figures of the LaTeX report.
//
// The expected references are CoolPlot's (thermo/references.py at 7338f48): its golden compose records
// (cpp/tests/data/compose_badge_r290_ph*.json, compose_rails_water_ts_*.json, compose_rails_blend_ph_*.json,
// compose_rails_psychrometric_single_lines.json of the CoolPlot repository) and, for R134a, its rules applied to
// the entry TillnerRoth-JPCRD-1994 of CoolProp's bibliography.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "coolsolve/diagram/engine.h"
#include "coolsolve/diagram/references.h"
#include "coolsolve/runner.h"

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"

#include <coolplot/capi.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>

using namespace coolsolve;
using namespace coolsolve::diagram;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::StartsWith;
namespace fs = std::filesystem;

namespace {

json pureFluid(const std::string& name) {
    return {{"backend", "HEOS"}, {"components", {name}}, {"fractions", json::array()}, {"basis", "mole"}, {"interactions", json::array()}};
}

json mixture(const std::vector<std::string>& components, const std::vector<double>& fractions) {
    return {{"backend", "HEOS"}, {"components", components}, {"fractions", fractions}, {"basis", "mole"}, {"interactions", json::array()}};
}

/// The p-h chart of R134a with every family but the internal energy; `branding` null keeps CoolPlot's defaults.
json chartSpec(const json& branding = json()) {
    json spec = {{"diagram", "ph"},
                 {"units", "engineering"},
                 {"limits", "DEF"},
                 {"isolines", {{{"key", "T"}, {"num", 8}}, {{"key", "s"}, {"num", 8}}, {{"key", "Q"}, {"num", 5}}}},
                 {"settings", {{"labels", {{"mode", "rails"}, {"symbols", "unicode"}}}}}};
    if (!branding.is_null()) spec["settings"]["branding"] = branding;
    return spec;
}

/// The R134a states of a heat pump in CoolSolve's units (°C, Pa, J/kg) and its cycle definition.
SolutionValues heatPump() {
    std::unique_ptr<CoolProp::AbstractState> st(CoolProp::AbstractState::factory("HEOS", "R134a"));
    SolutionValues values;
    st->update(CoolProp::PQ_INPUTS, 2.0e5, 1.0);
    const double s1 = st->smass();
    values["P_1"] = 2.0e5;
    values["h_1"] = st->hmass();
    st->update(CoolProp::PSmass_INPUTS, 1.0e6, s1);
    values["P_2"] = 1.0e6;
    values["h_2"] = st->hmass();
    st->update(CoolProp::PQ_INPUTS, 1.0e6, 0.0);
    values["P_3"] = 1.0e6;
    values["h_3"] = st->hmass();
    values["P_4"] = 2.0e5;
    values["h_4"] = st->hmass();
    return values;
}

json heatPumpCycle() {
    json states = json::array();
    for (int i = 1; i <= 4; ++i)
        states.push_back({{"number", i}, {"vars", {{"p", "P_" + std::to_string(i)}, {"h", "h_" + std::to_string(i)}}}});
    return {{"name", "Heat pump"}, {"states", states}};
}

json definitions(const json& diagrams) {
    return {{"format", "coolsolve-diagrams"}, {"version", 1}, {"diagrams", diagrams}};
}

fs::path examplesDir() {
    const char* env = std::getenv("COOLSOLVE_EXAMPLES_DIR");
    if (env && env[0] != '\0') return fs::path(env);
    const fs::path here = fs::path(__FILE__).parent_path().parent_path() / "examples";
    return fs::exists(here) ? here : fs::path("..") / "examples";
}

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

/// A PDF of coolplot-core: a header, an embedded TrueType subset, the end of file.
void checkPdf(const std::string& pdf) {
    CHECK_THAT(pdf, StartsWith("%PDF-"));
    CHECK_THAT(pdf, ContainsSubstring("/FontFile2"));
    CHECK_THAT(pdf, ContainsSubstring("/BaseFont /COOLPL+"));
    CHECK(pdf.find("%%EOF") != std::string::npos);
}

}  // namespace

// ---------------------------------------------------------------------------
// The badge lines (thermo/references.py)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram references: LaTeX of BibTeX fields as text", "[diagram]") {
    CHECK(latexToText("Wagner, W. and Pru{\\ss}, A.") == "Wagner, W. and Pruß, A.");
    CHECK(latexToText("M{\\\"u}ller") == "Müller");
    CHECK(latexToText("\\\"{U}ber \\'e \\` a") == "Über é à");
    CHECK(latexToText("{\\v s}t\\v{c} \\r{A}ngstr\\\"om") == "štč Ångström");
    CHECK(latexToText("\\vs") == "");  // a command, not an accent: \v needs a space or braces before its letter
    CHECK(latexToText("CO$_2$ and $^{13}$C") == "CO2 and 13C");
    CHECK(latexToText("Heat \\& Mass  Transfer\n 3141--3180") == "Heat & Mass Transfer 3141-3180");
    CHECK(latexToText("{\\textit{Ann. Phys.}}") == "Ann. Phys.");
    CHECK(latexToText("Ga{\\l}ecki, {\\O}stergaard, \\aa{}rhus") == "Gałecki, Østergaard, århus");
}

TEST_CASE("Diagram references: the equation of state as CoolPlot cites it", "[diagram]") {
    REQUIRE(hasBibliography());

    SECTION("pure fluids, short (the default since CoolPlot 1.5) and full") {
        CHECK(fluidReferences(pureFluid("R290")) == std::vector<std::string>{"R290: Lemmon et al. (2009)"});
        CHECK(fluidReferences(pureFluid("R290"), "full") ==
              std::vector<std::string>{"R290: Lemmon et al., J. Chem. Eng. Data 54 (2009) 3141-3180"});
        CHECK(fluidReferences(pureFluid("R134a")) == std::vector<std::string>{"R134a: Tillner-Roth and Baehr (1994)"});
        CHECK(fluidReferences(pureFluid("R134a"), "full") ==
              std::vector<std::string>{"R134a: Tillner-Roth and Baehr, J. Phys. Chem. Ref. Data 23 (1994) 657-729"});
        CHECK(fluidReferences(pureFluid("Water")) == std::vector<std::string>{"Water: Wagner and Pruß (2002)"});
    }
    SECTION("the fluids of CoolSolve's names: the components of the chart's fluid record") {
        CHECK(fluidReferences(chartFluid("R134a").record()) == std::vector<std::string>{"R134a: Tillner-Roth and Baehr (1994)"});
        CHECK(fluidReferences(chartFluid("Steam").record()) == std::vector<std::string>{"Water: Wagner and Pruß (2002)"});
    }
    SECTION("a blend given by its mole fractions is named after CoolProp's predefined blend of that composition") {
        CHECK(fluidReferences(mixture({"R32", "R125"}, {0.7, 0.3})) ==
              std::vector<std::string>{"R410A: R32 50.3 %, R125 49.7 % by mass", "R32: Tillner-Roth and Yokozeki (1997)",
                                       "R125: Lemmon and Jacobsen (2005)", "Mixture model: Lemmon and Jacobsen (2004)"});
        CHECK(fluidReferences(chartFluid("R32[0.7]&R125[0.3]").record()) == fluidReferences(mixture({"R32", "R125"}, {0.7, 0.3})));
        CHECK(blendName({"R32", "R125"}, {0.5, 0.5}) == "R410A");
        CHECK(blendName({"R32", "R125"}, {0.6, 0.4}).empty());
    }
    SECTION("a mixture without a predefined blend: its composition by mole") {
        CHECK(fluidReferences(mixture({"R32", "R1234yf"}, {0.5, 0.5})) ==
              std::vector<std::string>{"Mixture: R32 50 %, R1234yf 50 % by mole", "R32: Tillner-Roth and Yokozeki (1997)",
                                       "R1234yf: Lemmon and Akasaka (2022)", "Mixture model: Bell (2023)"});
    }
    SECTION("a predefined blend, humid air, other backends") {
        const std::vector<std::string> r407c = fluidReferences(pureFluid("R407C.mix"));
        REQUIRE(r407c.size() == 5);
        CHECK_THAT(r407c[0], StartsWith("R407C: R32 23 %, R125 25 %, R134a 52 % by mass"));
        CHECK_THAT(r407c[4], StartsWith("Mixture model: "));
        CHECK(fluidReferences({{"model", "humid_air"}, {"p_Pa", 101325.0}}) == std::vector<std::string>{"Humid air: Herrmann et al. (2009)"});
        json incomp = pureFluid("MEG");
        incomp["backend"] = "INCOMP";
        incomp["fractions"] = {0.3};
        incomp["basis"] = "mass";
        CHECK(fluidReferences(incomp) == std::vector<std::string>{"MEG 30 % by mass: CoolProp incompressible fluid"});
        json tabular = pureFluid("R134a");
        tabular["backend"] = "BICUBIC&HEOS";
        CHECK(fluidReferences(tabular) == fluidReferences(pureFluid("R134a")));
    }
    SECTION("a key the bibliography does not have reads as the key; an unknown style is refused") {
        CHECK(cite("No-Such-Key-2042") == "No-Such-Key-2042");
        CHECK(cite("Lemmon-JCED-2009") == "Lemmon et al. (2009)");
        CHECK_THROWS_AS(fluidReferences(pureFluid("R290"), "long"), DiagramError);
    }
}

TEST_CASE("Diagram references: the versions line and the badge lines of a spec", "[diagram]") {
    // CoolPlot names its own version; CoolSolve draws with coolplot-core, so the core's version stands for it
    std::string coolprop = CoolProp::get_global_param_string("version");
    if (coolprop.find("dev") != std::string::npos) coolprop += " (" + CoolProp::get_global_param_string("gitrevision").substr(0, 7) + ")";
    const std::string versions = "CoolProp " + coolprop + ", CoolPlot " + coolplot_version();
    CHECK(versionsLine() == versions);
    CHECK_THAT(versions, StartsWith("CoolProp 8."));
    CHECK_THAT(versions, ContainsSubstring(", CoolPlot 1.5.0"));

    json spec = chartSpec();
    spec["fluid"] = pureFluid("R290");
    CHECK(sourceLines(spec) == std::vector<std::string>{"R290: Lemmon et al. (2009)", versions});  // CoolPlot's defaults
    spec["settings"]["branding"] = {{"references", "full"}, {"versions", false}};
    CHECK(sourceLines(spec) == std::vector<std::string>{"R290: Lemmon et al., J. Chem. Eng. Data 54 (2009) 3141-3180"});
    spec["settings"]["branding"] = {{"references", "none"}, {"versions", false}};
    CHECK(sourceLines(spec).empty());
    spec["settings"]["branding"] = {{"references", "none"}};
    CHECK(sourceLines(spec) == std::vector<std::string>{versions});
    // A record before CoolPlot 1.2 (its other settings, no branding) keeps the chart without badge it was drawn as
    spec["settings"] = {{"calc", {{"points", 100}}}};
    CHECK(sourceLines(spec).empty());
    spec["settings"] = {{"branding", {{"references", "long"}}}};
    CHECK_THROWS_WITH(sourceLines(spec), ContainsSubstring("Unknown references style 'long'"));
}

TEST_CASE("Diagram scene: the badge names the equation of state and the versions", "[diagram]") {
    SceneRequest request;
    request.fluid = "R134a";
    request.spec = chartSpec();
    request.figureSizePt = std::make_pair(600.0, 420.0);
    const json response = scene(request);
    const std::string svg = response.at("svg");
    CHECK_THAT(svg, ContainsSubstring("R134a: Tillner-Roth and Baehr (1994)"));
    CHECK_THAT(svg, ContainsSubstring(versionsLine()));
    CHECK(response.at("diagnostics").at("warnings").empty());  // no warning about the badge any more
    // The spec as composed keeps CoolPlot's defaults, so that CoolPlot draws the same badge
    CHECK_FALSE(response.at("spec").at("settings").contains("branding"));

    request.spec = chartSpec({{"references", "none"}, {"versions", false}});
    const std::string plain = scene(request).at("svg");
    CHECK(plain.find("Tillner-Roth") == std::string::npos);
    CHECK(plain.find("CoolPlot 1.5.0") == std::string::npos);
    request.spec = chartSpec({{"references", "long"}});
    CHECK_THROWS_WITH(scene(request), ContainsSubstring("Unknown references style"));
}

// ---------------------------------------------------------------------------
// Chart files (POST /api/v1/diagram/export)
// ---------------------------------------------------------------------------

TEST_CASE("Diagram files: SVG and PDF of a solved cycle", "[diagram]") {
    const SolutionValues values = heatPump();
    SceneRequest request;
    request.fluid = "R134a";
    request.spec = chartSpec();
    request.cycles = json::array({heatPumpCycle()});
    request.solution = &values;

    SECTION("SVG: the document of the scene, at the figure size given") {
        request.figureSizePt = std::make_pair(600.0, 420.0);
        std::vector<std::string> warnings;
        const std::string svg = chartFile(request, "svg", &warnings);
        CHECK(svg == scene(request).at("svg").get<std::string>());
        CHECK_THAT(svg, ContainsSubstring("data-item=\"cycle/Heat pump/loop/0\""));
        CHECK_THAT(svg, ContainsSubstring("width=\"800\""));
        CHECK(warnings.empty());
    }
    SECTION("PDF: coolplot-core's writer with a TrueType font of the platform, CoolPlot's page by default") {
        const std::string pdf = chartFile(request, "pdf");
        checkPdf(pdf);
        CHECK_THAT(pdf, ContainsSubstring("/MediaBox [0 0 576 432]"));
        request.figureSizePt = std::make_pair(400.0, 300.0);
        const std::string small = chartFile(request, "pdf");
        CHECK_THAT(small, ContainsSubstring("/MediaBox [0 0 400 300]"));
        // The same chart gives the same bytes (no date, no file identifier)
        CHECK(chartFile(request, "pdf") == small);
    }
    SECTION("the cycle left out without a solve, with a warning") {
        request.solution = nullptr;
        std::vector<std::string> warnings;
        const std::string svg = chartFile(request, "svg", &warnings);
        CHECK(svg.find("cycle/Heat pump") == std::string::npos);
        REQUIRE(warnings.size() == 1);
        CHECK_THAT(warnings[0], ContainsSubstring("is not drawn"));
    }
    SECTION("formats and file names") {
        CHECK_THROWS_WITH(chartFile(request, "png"), ContainsSubstring("\"svg\" or \"pdf\""));
        CHECK(chartFormat("chart.svg") == "svg");
        CHECK(chartFormat("/tmp/Chart.PDF") == "pdf");
        CHECK_THROWS_AS(chartFormat("chart.png"), DiagramError);
        CHECK_THROWS_AS(chartFormat("chart"), DiagramError);
    }
}

// ---------------------------------------------------------------------------
// coolsolve --diagram (writeDiagramFile) and the figures of the report
// ---------------------------------------------------------------------------

TEST_CASE("Diagram definitions: the diagram of a document by name", "[diagram]") {
    const json first = {{"name", "Cycle on log(p)-h"}, {"fluid", "R134a"}, {"spec", chartSpec()}, {"cycles", json::array()}};
    const json second = {{"name", "T-s"}, {"fluid", "Water"}, {"spec", {{"diagram", "Ts"}}}};
    const json document = definitions({first, second});
    CHECK(&findDefinition(document, "") == &document.at("diagrams")[0]);
    CHECK(&findDefinition(document, "T-s") == &document.at("diagrams")[1]);
    CHECK_THROWS_WITH(findDefinition(document, "p-h"), ContainsSubstring("no diagram named 'p-h' (they have 'Cycle on log(p)-h', 'T-s')"));
    CHECK_THROWS_WITH(findDefinition(definitions(json::array()), ""), ContainsSubstring("no diagram"));
    CHECK_THROWS_AS(findDefinition(json::array(), ""), DiagramError);
    const SceneRequest request = definitionRequest(second);
    CHECK(request.fluid == "Water");
    CHECK(request.spec.at("diagram") == "Ts");
    CHECK(request.cycles == json::array());
    CHECK_THROWS_WITH(definitionRequest(json{{"spec", chartSpec()}}), ContainsSubstring("\"fluid\""));
}

TEST_CASE("coolsolve --diagram: a solved model drawn as SVG and PDF", "[diagram][cli]") {
    // The model is solved as the command line solves it, on a copy (a run writes next to its model)
    const fs::path dir = fs::temp_directory_path() / "coolsolve_test_diagram_cli";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const fs::path model = dir / "refrigeration1.eescode";
    fs::copy_file(examplesDir() / "refrigeration1.eescode", model);
    CoolSolveRunner runner(model.string());
    SolverOptions options;
    options.verbose = false;
    REQUIRE(runner.run(options, false));
    const SolutionValues solution(runner.getSolveResult().variables.begin(), runner.getSolveResult().variables.end());

    json states = json::array();
    for (int i = 1; i <= 4; ++i)
        states.push_back({{"number", i}, {"vars", {{"p", "P_" + std::to_string(i)}, {"h", "h_" + std::to_string(i)}}}});
    const json document = definitions({{{"name", "Water"}, {"fluid", "Water"}, {"spec", {{"diagram", "Ts"}}}},
                                       {{"name", "Refrigeration"}, {"fluid", "R134a"}, {"spec", chartSpec()},
                                        {"cycles", {{{"name", "Cycle"}, {"states", states}}}}}});

    const fs::path svg = dir / "chart.svg";
    CHECK(writeDiagramFile(document, "Refrigeration", solution, "HEOS", svg.string()).empty());
    const std::string text = readFile(svg);
    CHECK_THAT(text, StartsWith("<svg"));
    CHECK_THAT(text, ContainsSubstring("data-item=\"cycle/Cycle/loop/0\""));
    CHECK_THAT(text, ContainsSubstring("R134a: Tillner-Roth and Baehr (1994)"));

    const fs::path pdf = dir / "chart.pdf";
    writeDiagramFile(document, "", solution, "HEOS", pdf.string());  // the first diagram
    checkPdf(readFile(pdf));

    CHECK_THROWS_AS(writeDiagramFile(document, "Nope", solution, "HEOS", svg.string()), DiagramError);
    CHECK_THROWS_AS(writeDiagramFile(document, "", solution, "HEOS", (dir / "chart.png").string()), DiagramError);
    CHECK_THROWS_AS(writeDiagramFile(document, "", solution, "HEOS", (dir / "no" / "such" / "dir.svg").string()), std::runtime_error);
    fs::remove_all(dir);
}

TEST_CASE("Diagram figures of the LaTeX report: every diagram as a PDF in the print theme", "[diagram][latex]") {
    const SolutionValues values = heatPump();
    const json document = definitions({{{"name", "Heat pump"}, {"fluid", "R134a"}, {"spec", chartSpec()}, {"cycles", {heatPumpCycle()}}},
                                       {{"name", "Unknown"}, {"fluid", "Unobtainium"}, {"spec", chartSpec()}},
                                       {{"name", "Water"}, {"fluid", "Water"}, {"spec", {{"diagram", "Ts"}}}}});
    std::vector<std::string> warnings;
    const auto files = reportDiagrams(document, &values, "HEOS", warnings);
    REQUIRE(files.size() == 2);
    CHECK(files[0].first == "thermo_diagrams.pdf");
    CHECK(files[1].first == "thermo_diagram_3.pdf");  // the place of the diagram in the document
    for (const auto& [name, pdf] : files) checkPdf(pdf);
    REQUIRE(warnings.size() == 1);
    CHECK_THAT(warnings[0], ContainsSubstring("'Unknown' is left out of the report"));
    // The print theme: the same chart as the PDF of the spec with look "print"
    SceneRequest request = definitionRequest(document.at("diagrams")[0]);
    request.solution = &values;
    request.spec["settings"]["look"] = "print";
    CHECK(chartFile(request, "pdf") == files[0].second);
    request.spec["settings"]["look"] = "default";
    CHECK(chartFile(request, "pdf") != files[0].second);

    warnings.clear();
    CHECK(reportDiagrams(definitions(json::array()), &values, "HEOS", warnings).empty());
    CHECK(reportDiagrams(json::object(), nullptr, "HEOS", warnings).empty());
    CHECK(warnings.empty());
}
