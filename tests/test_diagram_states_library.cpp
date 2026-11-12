/**
 * The grouping of a model's variables into states (src/diagram/states.cpp, docs/coolplot_integration.md §4.2,
 * decision D6) on the models of the CoolSolve Library, against the former name heuristic of the Diagrams tab.
 *
 * The plan validates the rule on the library models before it replaces the heuristic. For every runnable model
 * (the main file of a model whose status is verified or modified) the model is parsed and its IR built (no solve:
 * the rule reads the property calls of the equations), then three groupings are compared, fluid by fluid:
 *   - the strict rule: calls sharing two variables of independent properties are one state ((T, p) alone does not
 *     merge; groupStates(.., false));
 *   - the rule as used, with "(T, p) merges two calls without vapour quality" (groupStates(..), the default);
 *   - the name heuristic (T_1, P_1, h_1 -> state 1), on the variables inferVariables() types.
 * The comparison is written to ./diagram_states_library.md ($COOLSOLVE_DIAGRAM_STATES_REPORT). It measures, it does
 * not judge: the test fails only when no model can be read.
 *
 * Hidden tag; run with
 *     ./coolsolve_tests "[diagram-states-library]"
 * The library is $COOLSOLVE_LIBRARY_DIR (a checkout with library.json and models/), else the snapshot embedded in
 * the binary (the checkout CMake embedded, as for [library-regression]).
 */

#include <catch2/catch_test_macros.hpp>

#include "coolsolve/diagram/states.h"
#include "coolsolve/include_resolver.h"
#include "coolsolve/ir.h"
#include "coolsolve/library_store.h"
#include "coolsolve/parser.h"
#include "coolsolve/variable_inference.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace coolsolve;
using namespace coolsolve::diagram;
namespace fs = std::filesystem;

namespace {

/// A library checkout on disk as a LibraryFileSource (the files of the root and of models/).
LibraryFileSource folderSource(const fs::path& root) {
    LibraryFileSource src;
    src.read = [root](const std::string& path, std::string& content) {
        std::ifstream in(root / path, std::ios::binary);
        if (!in) return false;
        std::ostringstream buffer;
        buffer << in.rdbuf();
        content = buffer.str();
        return true;
    };
    src.list = [root](const std::string& prefix) {
        std::vector<std::string> files;
        std::error_code ec;
        for (fs::directory_iterator it(root, ec), end; !ec && it != end; it.increment(ec))
            if (it->is_regular_file(ec)) files.push_back(it->path().filename().generic_string());
        for (fs::recursive_directory_iterator it(root / "models", ec), end; !ec && it != end; it.increment(ec))
            if (it->is_regular_file(ec)) files.push_back(fs::relative(it->path(), root, ec).generic_string());
        std::vector<std::string> out;
        for (const auto& f : files)
            if (f.compare(0, prefix.size(), prefix) == 0) out.push_back(f);
        std::sort(out.begin(), out.end());
        return out;
    };
    return src;
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

using VarSet = std::set<std::string>;

VarSet varsOf(const ModelState& state) {
    VarSet out;
    for (const auto& [_, var] : state.vars) out.insert(lower(var));
    return out;
}

std::string describe(const ModelState& state) {
    std::string out = "{";
    for (const auto& [key, var] : state.vars) out += (out.size() > 1 ? ", " : "") + key + ": " + var;
    return out + "}";
}

struct ModelRow {
    std::string id, name;
    size_t rule = 0, ruleTp = 0, names = 0, identical = 0;
    size_t identicalTp = 0;   ///< heuristic states identical to a state of the rule + (T,p)
    size_t namesSubset = 0;   ///< heuristic states inside one rule state (the rule finds more variables)
    size_t namesSpan = 0;     ///< heuristic states over several rule states (the names mix states)
    size_t namesOutside = 0;  ///< heuristic states with a variable in no rule state
    size_t ruleOnly = 0;      ///< rule states that no heuristic state matches exactly
    double groupMs = 0.0;
    std::vector<std::string> examples;
};

}  // namespace

TEST_CASE("Diagram states: the property-call rule against the name heuristic on the library models",
          "[.][diagram-states-library]") {
    std::unique_ptr<LibraryStore> own;
    const LibraryStore* store = &embeddedLibrary();
    if (const char* dir = std::getenv("COOLSOLVE_LIBRARY_DIR"); dir && *dir) {
        own = std::make_unique<LibraryStore>(folderSource(dir));
        store = own.get();
    }
    if (!store->available()) SKIP("No CoolSolve Library: set COOLSOLVE_LIBRARY_DIR or build with the library embedded");
    std::string indexText;
    REQUIRE(store->readFile("library.json", indexText));
    const nlohmann::json index = nlohmann::json::parse(indexText);

    std::vector<ModelRow> rows;
    std::vector<std::string> problems;
    size_t runnable = 0;
    for (const nlohmann::json& entry : index.at("models")) {
        const std::string status = entry.value("status", "");
        if (status != "verified" && status != "modified") continue;
        ++runnable;
        ModelRow row;
        row.id = entry.value("id", "");
        row.name = entry.value("name", "");
        LibraryModel model;
        std::string error;
        if (!store->loadModel(row.id, "", model, error)) {
            problems.push_back(row.id + ": " + error);
            continue;
        }
        EESParser parser;
        parser.setIncludeResolver(IncludeResolver::forLibraryModel(row.id, model.file, store));
        const ParseResult parsed = parser.parse(model.eescode, model.file);
        if (!parsed.success) {
            problems.push_back(row.id + ": does not parse");
            continue;
        }
        ModelStates rule, ruleTp, names;
        try {
            IR ir = IR::fromAST(parsed.program);
            const auto start = std::chrono::steady_clock::now();
            rule = groupStates(ir, nullptr, false);
            row.groupMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
            ruleTp = groupStates(ir);
            inferVariables(ir);
            names = nameHeuristicStates(ir);
        } catch (const std::exception& e) {
            problems.push_back(row.id + ": " + e.what());
            continue;
        }
        std::set<std::string> fluids;
        for (const auto& [fluid, _] : rule) fluids.insert(fluid);
        for (const auto& [fluid, _] : names) fluids.insert(fluid);  // both name fluids by FluidRegistry::getName
        for (const auto& [fluid, states] : ruleTp) row.ruleTp += states.size();
        for (const std::string& fluid : fluids) {
            const std::vector<ModelState> none;
            const auto r = rule.find(fluid);
            const auto n = names.find(fluid);
            const std::vector<ModelState>& ruleStates = r == rule.end() ? none : r->second;
            const std::vector<ModelState>& nameStates = n == names.end() ? none : n->second;
            row.rule += ruleStates.size();
            row.names += nameStates.size();
            std::vector<VarSet> ruleSets, tpSets;
            for (const ModelState& s : ruleStates) ruleSets.push_back(varsOf(s));
            if (const auto t = ruleTp.find(fluid); t != ruleTp.end())
                for (const ModelState& s : t->second) tpSets.push_back(varsOf(s));
            std::set<size_t> matched;
            for (const ModelState& s : nameStates) {
                const VarSet vars = varsOf(s);
                if (std::find(tpSets.begin(), tpSets.end(), vars) != tpSets.end()) ++row.identicalTp;
                const auto same = std::find(ruleSets.begin(), ruleSets.end(), vars);
                if (same != ruleSets.end()) {
                    ++row.identical;
                    matched.insert(static_cast<size_t>(same - ruleSets.begin()));
                    continue;
                }
                std::set<size_t> holders;
                bool outside = false;
                for (const std::string& var : vars) {
                    bool found = false;
                    for (size_t i = 0; i < ruleSets.size(); ++i)
                        if (ruleSets[i].count(var)) {
                            holders.insert(i);
                            found = true;
                        }
                    outside = outside || !found;
                }
                std::string why;
                if (outside) {
                    ++row.namesOutside;
                    why = "has variables in no state of the rule";
                } else if (holders.size() > 1) {
                    ++row.namesSpan;
                    why = "mixes " + std::to_string(holders.size()) + " states of the rule";
                } else {
                    ++row.namesSubset;
                    why = "is part of the rule's " + describe(ruleStates[*holders.begin()]);
                }
                if (row.examples.size() < 3)
                    row.examples.push_back(fluid + ": name group " + std::to_string(s.number) + (s.label.empty() ? "" : " (" + s.label + ")") +
                                           " " + describe(s) + " " + why);
            }
            row.ruleOnly += ruleStates.size() - matched.size();
        }
        rows.push_back(row);
    }

    // The report
    size_t identicalTp = 0;
    size_t total = 0, totalTp = 0, totalNames = 0, identical = 0, subset = 0, span = 0, outside = 0, ruleOnly = 0, withStates = 0;
    double maxMs = 0.0, sumMs = 0.0;
    std::vector<double> times;
    for (const ModelRow& row : rows) {
        total += row.rule;
        totalTp += row.ruleTp;
        totalNames += row.names;
        identical += row.identical;
        identicalTp += row.identicalTp;
        subset += row.namesSubset;
        span += row.namesSpan;
        outside += row.namesOutside;
        ruleOnly += row.ruleOnly;
        withStates += (row.rule + row.names) > 0;
        maxMs = std::max(maxMs, row.groupMs);
        sumMs += row.groupMs;
        times.push_back(row.groupMs);
    }
    std::sort(times.begin(), times.end());
    auto percentile = [&times](double q) { return times.empty() ? 0.0 : times[static_cast<size_t>(q * (times.size() - 1))]; };
    std::string slowest;
    for (const ModelRow& row : rows)
        if (row.groupMs == maxMs) slowest = row.id;
    const char* reportEnv = std::getenv("COOLSOLVE_DIAGRAM_STATES_REPORT");
    const std::string reportPath = reportEnv && *reportEnv ? reportEnv : "diagram_states_library.md";
    std::ofstream md(reportPath);
    md << std::fixed << std::setprecision(3);
    md << "# States of the library models: property calls against the name heuristic\n\n"
       << "Written by `./coolsolve_tests \"[diagram-states-library]\"` (tests/test_diagram_states_library.cpp). The rule of "
          "docs/coolplot_integration.md §4.2 (calls sharing two variables of independent properties are one state; T and p "
          "alone do not merge) is compared, fluid by fluid, with the name heuristic of the former Diagrams tab (T_1, P_1, "
          "h_1 -> state 1). \"Rule + (T,p)\" also merges two calls sharing T and p when neither has a vapour quality. Only "
          "states of two independent properties (rule) or two properties (heuristic) count.\n\n"
       << "Library: " << (own ? std::string(std::getenv("COOLSOLVE_LIBRARY_DIR")) : std::string("the snapshot embedded in the binary"))
       << ", snapshot " << store->snapshot().value("content_hash", std::string("?")) << ".\n\n"
       << "## Summary\n\n| | Count |\n|---|---:|\n"
       << "| Runnable models (verified, modified) | " << runnable << " |\n"
       << "| Models read (parsed, IR built) | " << rows.size() << " |\n"
       << "| Models with states (either way) | " << withStates << " |\n"
       << "| States, rule | " << total << " |\n"
       << "| States, rule + (T,p) | " << totalTp << " |\n"
       << "| States, name heuristic | " << totalNames << " |\n"
       << "| Heuristic states identical to a rule state (same variables) | " << identical << " |\n"
       << "| Heuristic states identical to a rule + (T,p) state | " << identicalTp << " |\n"
       << "| Heuristic state inside one rule state (the rule finds more variables) | " << subset << " |\n"
       << "| Heuristic state over several rule states (the names mix states) | " << span << " |\n"
       << "| Heuristic state with variables in no rule state | " << outside << " |\n"
       << "| Rule states without an identical heuristic state | " << ruleOnly << " |\n"
       << "| Grouping time per model, median / 95 % / max (ms) | " << percentile(0.5) << " / " << percentile(0.95) << " / " << maxMs
       << " (" << slowest << ") |\n\n"
       << "The first models pay the loading of CoolProp's fluids by the fluid registry (the server has done it before a "
          "solve).\n\n";
    if (!problems.empty()) {
        md << "## Models not read\n\n";
        for (const std::string& p : problems) md << "- " << p << "\n";
        md << "\n";
    }
    md << "## Models\n\n| Model | Rule | Rule + (T,p) | Names | Identical | Inside | Mixed | Outside | Rule only | Examples |\n"
       << "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n";
    for (const ModelRow& row : rows) {
        if (row.rule + row.names == 0) continue;
        std::string examples;
        for (const std::string& e : row.examples) examples += (examples.empty() ? "" : "<br>") + e;
        md << "| " << row.id << " " << row.name << " | " << row.rule << " | " << row.ruleTp << " | " << row.names << " | " << row.identical
           << " | " << row.namesSubset << " | " << row.namesSpan << " | " << row.namesOutside << " | " << row.ruleOnly << " | " << examples
           << " |\n";
    }
    md.close();
    WARN("States of " << rows.size() << " library models: rule " << total << ", rule + (T,p) " << totalTp << ", names " << totalNames
                      << ", identical " << identical << " (" << identicalTp << " with (T,p)) (inside " << subset << ", mixed " << span << ", outside " << outside
                      << ", rule only " << ruleOnly << "); grouping median " << percentile(0.5) << " ms, 95 % " << percentile(0.95)
                      << " ms, max " << maxMs << " ms (" << slowest << "); report " << reportPath);
    REQUIRE(!rows.empty());
}
