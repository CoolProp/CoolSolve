/**
 * Library regression test (CS-FEAT-LIB-TEST)
 *
 * Solves every runnable model of the CoolSolve Library (https://github.com/IntSusEnergySystems/CoolSolve_Library)
 * and compares the fresh solution with the `.sol` stored next to the model, the library's regression baseline.
 * It is the in-process twin of the library's `tools/test_models.py`, which drives the CLI: both select the same
 * targets and apply the same comparison, so they give the same verdicts for the same CoolSolve build.
 *
 * Targets (as in tools/test_models.py):
 *   - the main file (`main_file` of `model.json`) of every model whose status is verified or modified, `CSL-xxxx`;
 *   - every other `*.eescode` file of a model folder that has a sibling `.sol` of the same stem (a runnable
 *     variant such as `<name>_coolsolve.eescode`), whatever the status of the model, `CSL-xxxx:variant`.
 * Each target is solved from a copy of its model folder (so its `.initials`, `coolsolve.conf`, companion CSV
 * tables and `$INCLUDE`d project files are used as the CLI uses them) with the CLI's pipeline: solver options
 * from `coolsolve.conf`, solver, then the solution verification that gates the `.sol` of the CLI.
 *
 * Comparison: for every numeric `name = value` line of the baseline (names case-insensitive, strings and the
 * IntegralTable rows are not values), the fresh value must satisfy
 *     |new - baseline| <= 1e-9 + 1e-6 * max(|new|, |baseline|)
 * and exist. A value that became undefined or infinite is a deviation too (tools/test_models.py lets it pass).
 * The test fails when a target does not solve, deviates, misses a variable, or has a baseline without any
 * numeric value (an empty or malformed `.sol` would compare nothing: tools/test_models.py lets it pass).
 * A main file without `.sol` is only solved.
 *
 * Slow (about five minutes), hence the hidden tag. Run with
 *     ./coolsolve_tests "[library-regression]"
 * The library is, in this order: $COOLSOLVE_LIBRARY_DIR; the checkout CMake embedded in the build
 * (the FetchContent copy in .fetchcontent_cache/, or FETCHCONTENT_SOURCE_DIR_COOLSOLVELIBRARY); else the test is
 * skipped. The report is written to ./library_regression.md ($COOLSOLVE_LIBRARY_REPORT).
 *
 * Environment (all optional):
 *   COOLSOLVE_LIBRARY_DIR            library checkout to test (a folder holding library.json and models/)
 *   COOLSOLVE_LIBRARY_TEST_FILTER    comma-separated targets: "CSL-0001" (the model and its variants),
 *                                    "CSL-0009:coolsolve" (one variant); default: every target
 *   COOLSOLVE_LIBRARY_TEST_EXCLUDE   comma-separated targets to leave out, e.g. the slow "CSL-0009:coolsolve"
 *   COOLSOLVE_LIBRARY_REPORT         report file (default library_regression.md in the working directory)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>

#include "coolsolve/evaluator.h"
#include "coolsolve/include_resolver.h"
#include "coolsolve/library_store.h"
#include "coolsolve/runner.h"
#include "coolsolve/solution_checker.h"
#include "coolsolve/solver.h"
#include "coolsolve/version.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifndef COOLSOLVE_LIBRARY_SOURCE_DIR
#define COOLSOLVE_LIBRARY_SOURCE_DIR ""   // set by CMake: the library checkout embedded in the build
#endif

namespace fs = std::filesystem;

namespace {

// Defaults of tools/test_models.py (--rtol, --atol, --timeout)
constexpr double kRelTol = 1e-6;
constexpr double kAbsTol = 1e-9;
constexpr int kTimeoutSeconds = 600;   // per target, unless its coolsolve.conf sets one

const std::set<std::string> kRunnableStatuses = {"verified", "modified"};

// ----------------------------------------------------------------------------
// Small helpers
// ----------------------------------------------------------------------------

bool startsWith(const std::string& s, const std::string& prefix) {
    return s.compare(0, prefix.size(), prefix) == 0;
}

bool endsWith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string trim(const std::string& s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::set<std::string> splitList(const char* text) {
    std::set<std::string> out;
    if (!text) return out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item = trim(item);
        if (!item.empty()) out.insert(item);
    }
    return out;
}

std::string num(double v, int digits = 6) {
    std::ostringstream os;
    os << std::setprecision(digits) << v;
    return os.str();
}

/// A message on one line, for a table cell: line breaks and `|` removed, cut to `maxLen` characters.
std::string oneLine(const std::string& text, size_t maxLen = 300) {
    std::string out;
    bool space = false;
    for (char c : text) {
        if (c == '\n' || c == '\r' || c == '\t' || c == ' ') {
            space = !out.empty();
            continue;
        }
        if (space) out += ' ';
        space = false;
        out += (c == '|') ? '/' : c;
    }
    if (out.size() > maxLen) out = out.substr(0, maxLen - 3) + "...";
    return out;
}

// SIGINT: finish the target being solved, write the report of what was run
std::atomic<bool> g_interrupted{false};
void libraryRegressionSigint(int) { g_interrupted = true; }

// ----------------------------------------------------------------------------
// The library on disk
// ----------------------------------------------------------------------------

/// The files of a library checkout for a LibraryStore (the `$INCLUDE` of library models), so that the
/// models are solved against the library under test and not against the snapshot embedded in the binary.
/// Listing: the files of the root and of `models/` (all that the store reads).
coolsolve::LibraryFileSource fileSystemSource(const fs::path& root) {
    coolsolve::LibraryFileSource src;
    src.read = [root](const std::string& path, std::string& content) {
        const fs::path file = root / path;
        std::error_code ec;
        if (!fs::is_regular_file(file, ec)) return false;
        std::ifstream in(file, std::ios::binary);
        if (!in) return false;
        std::ostringstream buffer;
        buffer << in.rdbuf();
        content = buffer.str();
        return true;
    };
    auto files = std::make_shared<std::vector<std::string>>();
    auto listed = std::make_shared<bool>(false);
    src.list = [root, files, listed](const std::string& prefix) {
        if (!*listed) {
            *listed = true;
            std::error_code ec;
            for (fs::directory_iterator it(root, ec), end; !ec && it != end; it.increment(ec))
                if (it->is_regular_file(ec)) files->push_back(it->path().filename().generic_string());
            for (fs::recursive_directory_iterator it(root / "models", ec), end; !ec && it != end; it.increment(ec))
                if (it->is_regular_file(ec))
                    files->push_back(fs::relative(it->path(), root, ec).generic_string());
            std::sort(files->begin(), files->end());
        }
        std::vector<std::string> out;
        for (const auto& f : *files)
            if (startsWith(f, prefix)) out.push_back(f);
        return out;
    };
    return src;
}

/// The `COOLPROP_BASELINE` constant of the library's tools/test_models.py: the CoolProp git revision (or
/// prefix) that produced the baselines. Empty when the library has no such file.
std::string libraryCoolPropBaseline(const fs::path& root) {
    std::ifstream in(root / "tools" / "test_models.py");
    std::string line;
    while (std::getline(in, line)) {
        if (!startsWith(line, "COOLPROP_BASELINE")) continue;
        const auto open = line.find('"');
        const auto close = open == std::string::npos ? open : line.find('"', open + 1);
        if (close != std::string::npos) return lower(line.substr(open + 1, close - open - 1));
    }
    return "";
}

/// Copy of a model folder without what the library's own regression leaves out too: `original/`,
/// `figures/` and the `*_coolsolve` debug folders.
void copyModelFolder(const fs::path& from, const fs::path& to) {
    fs::create_directories(to);
    for (const auto& entry : fs::directory_iterator(from)) {
        const std::string name = entry.path().filename().string();
        if (name == "original" || name == "figures" || endsWith(name, "_coolsolve")) continue;
        if (entry.is_directory())
            copyModelFolder(entry.path(), to / name);
        else if (entry.is_regular_file())
            fs::copy_file(entry.path(), to / name, fs::copy_options::overwrite_existing);
    }
}

// ----------------------------------------------------------------------------
// Targets
// ----------------------------------------------------------------------------

struct Target {
    std::string id;         ///< "CSL-0012"
    std::string name;
    std::string status;
    std::string variant;    ///< empty: the main file of the model
    std::string file;       ///< "x.eescode"
    fs::path folder;
    fs::path baseline;      ///< the `.sol` of the file; empty when there is none

    std::string label() const { return variant.empty() ? id : id + ":" + variant; }
};

struct Selection {
    std::vector<Target> targets;                  ///< to solve
    size_t mainFiles = 0, variants = 0;           ///< of `targets`
    size_t skippedNotRunnable = 0;                ///< main files of models whose status is not runnable
    std::vector<std::string> variantsWithoutSol;  ///< variants without baseline (not tested)
    std::vector<std::string> problems;            ///< unreadable model.json, unknown filter entries
    std::string restriction;                      ///< the filter and exclusion in force, empty when every target is selected
};

Selection selectTargets(const fs::path& root, const std::set<std::string>& filter,
                        const std::set<std::string>& exclude) {
    Selection sel;
    std::vector<fs::path> modelJsons;
    for (const auto& entry : fs::recursive_directory_iterator(root / "models"))
        if (entry.is_regular_file() && entry.path().filename() == "model.json") modelJsons.push_back(entry.path());
    std::sort(modelJsons.begin(), modelJsons.end());

    std::set<std::string> matched;   // filter entries that name an existing target (runnable or not)
    auto wanted = [&](const Target& t) {
        if (filter.empty()) return true;
        bool hit = false;
        for (const std::string& key : {t.id, t.label()})
            if (filter.count(key)) {
                matched.insert(key);
                hit = true;
            }
        return hit;
    };

    for (const auto& mj : modelJsons) {
        nlohmann::json meta;
        try {
            std::ifstream in(mj);
            meta = nlohmann::json::parse(in);
        } catch (const std::exception& e) {
            sel.problems.push_back(mj.string() + ": cannot read model.json: " + e.what());
            continue;
        }
        const fs::path folder = mj.parent_path();
        const std::string mainFile = meta.value("main_file", "");
        const std::string mainStem = fs::path(mainFile).stem().string();

        // The main file first, then the other *.eescode files by name (the order of tools/test_models.py)
        std::vector<fs::path> files = {folder / mainFile};
        std::vector<fs::path> others;
        for (const auto& f : fs::directory_iterator(folder))
            if (f.is_regular_file() && f.path().extension() == ".eescode" && f.path().filename() != mainFile)
                others.push_back(f.path());
        std::sort(others.begin(), others.end());
        files.insert(files.end(), others.begin(), others.end());

        for (const auto& f : files) {
            Target t;
            t.id = meta.value("id", "");
            t.name = meta.value("name", "");
            t.status = meta.value("status", "");
            t.file = f.filename().string();
            t.folder = folder;
            const bool isMain = (t.file == mainFile);
            if (!isMain) {
                const std::string stem = f.stem().string();
                t.variant = startsWith(stem, mainStem + "_") ? stem.substr(mainStem.size() + 1) : stem;
            }
            fs::path sol = fs::path(f).replace_extension(".sol");
            if (fs::exists(sol)) t.baseline = sol;

            if (!wanted(t) || exclude.count(t.label())) continue;
            if (isMain && !kRunnableStatuses.count(t.status)) {
                ++sel.skippedNotRunnable;
                continue;
            }
            if (!isMain && t.baseline.empty()) {
                sel.variantsWithoutSol.push_back(t.label());
                continue;
            }
            (isMain ? sel.mainFiles : sel.variants)++;
            sel.targets.push_back(std::move(t));
        }
    }
    for (const auto& entry : filter)
        if (!matched.count(entry)) sel.problems.push_back("COOLSOLVE_LIBRARY_TEST_FILTER: no such model or variant: " + entry);
    auto joined = [](const std::set<std::string>& items) {
        std::string out;
        for (const auto& item : items) out += (out.empty() ? "" : ",") + item;
        return out;
    };
    if (!filter.empty()) sel.restriction = "only " + joined(filter);
    if (!exclude.empty()) sel.restriction += std::string(sel.restriction.empty() ? "" : ", ") + "without " + joined(exclude);
    return sel;
}

// ----------------------------------------------------------------------------
// Solving: the pipeline of main.cpp (the CLI), in this process
// ----------------------------------------------------------------------------

/// A value as the CLI writes it to the `.sol` (12 digits after the point, scientific), which is also the
/// precision of the baselines.
double asWrittenToSol(double v) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(12) << v;
    return std::strtod(os.str().c_str(), nullptr);
}

struct SolveOutcome {
    bool ok = false;
    std::string message;                    ///< why it did not (a solver status and message)
    std::map<std::string, double> values;   ///< lower-cased name -> value of every numeric variable
};

SolveOutcome solveLikeCli(const fs::path& file, const coolsolve::LibraryStore& library) {
    SolveOutcome out;
    try {
        coolsolve::SolverOptions options;
        options.verbose = false;
        const fs::path conf = file.parent_path() / "coolsolve.conf";
        if (fs::exists(conf)) coolsolve::loadSolverOptionsFromFile(conf.string(), options);
        if (options.timeoutSeconds == 0) options.timeoutSeconds = kTimeoutSeconds;
        coolsolve::applyCoolPropConfig(options.coolpropConfig);

        coolsolve::CoolSolveRunner runner(file.string());
        runner.setIncludeResolver(coolsolve::IncludeResolver::forDirectory(
            file.parent_path().string(), file.filename().string(), &library));
        runner.run(options, false);

        if (!runner.isParseSuccess()) {
            const auto& errors = runner.getParseResult().errors;
            out.message = "Parse failed: " + (errors.empty() ? std::string("unknown error")
                          : "line " + std::to_string(errors[0].line) + ": " + errors[0].message);
            return out;
        }
        if (!runner.isIRSuccess()) {
            out.message = "IR construction failed";
            return out;
        }
        const auto& analysis = runner.getAnalysisResult();
        if (!analysis.success && !runner.hasIntegralResult()) {
            out.message = "Structural analysis failed: " + analysis.errorMessage;
            return out;
        }
        const auto& solve = runner.getSolveResult();
        if (!solve.success) {
            out.message = "Solver " + coolsolve::statusToString(solve.status) + ": " +
                          (solve.errorMessage.empty() ? std::string("no message") : solve.errorMessage);
            return out;
        }
        // The CLI writes its .sol only when this verification passes (tolerance 1e-3)
        const auto check = coolsolve::checkSolution(runner.getIR(), solve.variables, solve.stringVariables,
                                                    options.coolpropConfig, 1e-3, &runner.getLookupTableStore());
        if (!check.allSatisfied) {
            out.message = "Solution verification failed: " + std::to_string(check.violatedCount) +
                          " equation(s) violated";
            for (const auto& chk : check.checks)
                if (!chk.satisfied) {
                    out.message += ", e.g. " + chk.originalText + " (relative error " + num(chk.relativeError, 3) + ")";
                    break;
                }
            for (const auto& d : check.diagnostics.items())
                if (d.severity == coolsolve::DiagnosticSeverity::Error) {
                    out.message += "; " + d.message;
                    break;
                }
            return out;
        }
        for (const auto& [name, value] : solve.variables) out.values[lower(name)] = asWrittenToSol(value);
        out.ok = true;
    } catch (const std::exception& e) {
        out.message = std::string("Exception: ") + e.what();
    } catch (...) {
        out.message = "Unknown exception";
    }
    return out;
}

// ----------------------------------------------------------------------------
// Baseline and comparison (tools/test_models.py: read_sol, compare)
// ----------------------------------------------------------------------------

/// The numeric `name = value` lines of a `.sol`, names lower-cased. A line is `name = number`, with an
/// optional `"unit"` after it; strings (`name = 'text'`), comments and the rows of the IntegralTable are not.
std::map<std::string, double> readSol(const fs::path& path) {
    static const std::string kNumberChars = "-+0123456789.eEinfINFnaN";
    std::map<std::string, double> sol;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = trim(line.substr(0, eq));
        if (key.empty() || key.find_first_of(" \t\r\n\f\v") != std::string::npos) continue;
        const std::string rest = line.substr(eq + 1);
        size_t i = 0;
        while (i < rest.size() && std::isspace(static_cast<unsigned char>(rest[i]))) ++i;
        size_t j = i;
        while (j < rest.size() && kNumberChars.find(rest[j]) != std::string::npos) ++j;
        if (j == i) continue;
        const std::string token = rest.substr(i, j - i);
        std::string tail = trim(rest.substr(j));
        if (!tail.empty()) {
            if (tail.front() != '"' || tail.find('"', 1) != tail.size() - 1) continue;   // only a "unit" may follow
        }
        char* end = nullptr;
        const double value = std::strtod(token.c_str(), &end);
        if (end != token.c_str() + token.size()) continue;
        sol[lower(key)] = value;
    }
    return sol;
}

struct Comparison {
    size_t compared = 0;       ///< variables of the baseline
    size_t outside = 0;        ///< beyond the tolerance
    size_t missing = 0;        ///< in the baseline, not in the solution
    std::string worstName;     ///< the variable closest to (or farthest beyond) the tolerance
    double worstBaseline = 0, worstFresh = 0;
    double worstRatio = -1;    ///< |difference| / tolerance of the worst variable (> 1: outside)
    double worstRelative = 0;  ///< |difference| / max(|new|, |baseline|) of the worst variable
    std::vector<std::string> details;   ///< the first deviations, for the report
};

Comparison compareWithBaseline(const std::map<std::string, double>& fresh, const std::map<std::string, double>& baseline) {
    Comparison c;
    for (const auto& [name, a] : baseline) {
        ++c.compared;
        const auto it = fresh.find(name);
        if (it == fresh.end()) {
            ++c.missing;
            if (c.details.size() < 10) c.details.push_back(name + " missing from the solution");
            continue;
        }
        const double b = it->second;
        double ratio = 0, relative = 0;
        bool outside = false;
        if (a == b || (std::isnan(a) && std::isnan(b))) {
            // equal (the same infinity, or both undefined)
        } else if (!std::isfinite(a) || !std::isfinite(b)) {
            outside = true;   // finite against infinite or undefined
            ratio = relative = INFINITY;
        } else {
            const double scale = std::max(std::fabs(a), std::fabs(b));
            const double limit = kAbsTol + kRelTol * scale;
            outside = std::fabs(b - a) > limit;   // the comparison of tools/test_models.py
            ratio = std::fabs(b - a) / limit;
            relative = std::fabs(b - a) / scale;
        }
        if (outside) {
            ++c.outside;
            if (c.details.size() < 10) c.details.push_back(name + ": " + num(a) + " -> " + num(b));
        }
        if (ratio > c.worstRatio) {
            c.worstRatio = ratio;
            c.worstName = name;
            c.worstBaseline = a;
            c.worstFresh = b;
            c.worstRelative = relative;
        }
    }
    return c;
}

// ----------------------------------------------------------------------------
// Result of one target and the report
// ----------------------------------------------------------------------------

enum class Verdict { Pass, Deviation, Fail, NoBaseline };

const char* verdictName(Verdict v) {
    switch (v) {
        case Verdict::Pass: return "PASS";
        case Verdict::Deviation: return "DEVIATION";
        case Verdict::Fail: return "FAIL";
        case Verdict::NoBaseline: return "NO BASELINE";
    }
    return "?";
}

struct TargetResult {
    Target target;
    Verdict verdict = Verdict::Fail;
    std::string message;       ///< FAIL: why
    Comparison comparison;     ///< PASS, DEVIATION
    double seconds = 0;
};

/// "max dev 3.2e-09 (x)" / "3 of 120 beyond tolerance, worst x ..." / the solver message: one line.
std::string summaryOf(const TargetResult& r) {
    const Comparison& c = r.comparison;
    switch (r.verdict) {
        case Verdict::Fail: return oneLine(r.message);
        case Verdict::NoBaseline: return "solved, no `.sol` baseline";
        case Verdict::Pass:
            return c.worstRatio < 0 ? "no numeric value in the baseline"
                                    : "max deviation " + num(c.worstRelative, 3) + " (" + c.worstName + ")";
        case Verdict::Deviation: {
            std::string s;
            if (c.outside) s += std::to_string(c.outside) + " of " + std::to_string(c.compared) + " beyond tolerance";
            if (c.missing) s += std::string(s.empty() ? "" : ", ") + std::to_string(c.missing) + " missing";
            if (c.outside) s += "; worst " + c.worstName + ": " + num(c.worstBaseline) + " -> " + num(c.worstFresh) +
                                " (" + num(c.worstRelative, 3) + ")";
            return s;
        }
    }
    return "";
}

/// Per verdict (indexed by `Verdict`), the number of results.
std::array<size_t, 4> countVerdicts(const std::vector<TargetResult>& results) {
    std::array<size_t, 4> counts = {0, 0, 0, 0};
    for (const auto& r : results) ++counts[static_cast<size_t>(r.verdict)];
    return counts;
}

/// Solve one target from a copy of its folder in `scratch` and compare it with its baseline.
TargetResult runTarget(const Target& target, const fs::path& scratch, const coolsolve::LibraryStore& library) {
    TargetResult r;
    r.target = target;
    const auto start = std::chrono::steady_clock::now();
    std::string dirName = target.label();
    std::replace(dirName.begin(), dirName.end(), ':', '_');
    const fs::path work = scratch / dirName;
    SolveOutcome outcome;
    try {
        copyModelFolder(target.folder, work);
        outcome = solveLikeCli(work / target.file, library);
    } catch (const std::exception& e) {
        outcome.message = std::string("Exception: ") + e.what();
    }
    { std::error_code ec; fs::remove_all(work, ec); }
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (!outcome.ok) {
        r.verdict = Verdict::Fail;
        r.message = outcome.message;
    } else if (target.baseline.empty()) {
        r.verdict = Verdict::NoBaseline;
    } else {
        const auto baseline = readSol(target.baseline);
        if (baseline.empty()) {   // a baseline that compares nothing would pass whatever the solution is
            r.verdict = Verdict::Fail;
            r.message = "The baseline " + target.baseline.filename().string() +
                        " holds no numeric value (unreadable or malformed)";
        } else {
            r.comparison = compareWithBaseline(outcome.values, baseline);
            r.verdict = (r.comparison.outside + r.comparison.missing > 0) ? Verdict::Deviation : Verdict::Pass;
        }
    }
    return r;
}

void writeReport(const fs::path& path, const std::vector<TargetResult>& results, const Selection& sel,
                 const fs::path& libraryRoot, const nlohmann::json& snapshot, const std::string& coolpropWarning,
                 double totalSeconds, size_t notRun) {
    std::ofstream md(path);
    if (!md.is_open()) return;

    const auto counts = countVerdicts(results);
    const auto info = coolsolve::getVersionInfo();
    const std::time_t now = std::time(nullptr);
    char stamp[64] = "";
    std::strftime(stamp, sizeof(stamp), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    md << "# CoolSolve Library Regression Report\n\n";
    md << "- Generated: " << stamp << " (" << std::fixed << std::setprecision(1) << totalSeconds << " s)\n";
    md << "- CoolSolve " << info.coolsolve << ", CoolProp " << info.coolprop << " (git revision "
       << info.coolpropRevision.substr(0, 8) << ")\n";
    md << "- Library: `" << libraryRoot.string() << "`, " << snapshot.value("models", 0) << " models, snapshot "
       << snapshot.value("content_hash", std::string("unknown")) << "\n";
    md << "- Rule: |new - baseline| <= " << std::defaultfloat << kAbsTol << " + " << kRelTol
       << " * max(|new|, |baseline|) for every numeric variable of the `.sol` baseline, as `tools/test_models.py`\n";
    if (!sel.restriction.empty()) md << "- **Partial run**: " << sel.restriction << "\n";
    if (!coolpropWarning.empty()) md << "\n**" << coolpropWarning << "**\n";
    if (notRun) md << "\n**Interrupted: " << notRun << " target(s) were not run.**\n";

    md << "\n## Summary\n\n";
    md << "| Verdict | Targets |\n|---|---|\n";
    md << "| PASS | " << counts[0] << " |\n";
    md << "| DEVIATION | " << counts[1] << " |\n";
    md << "| FAIL | " << counts[2] << " |\n";
    md << "| NO BASELINE | " << counts[3] << " |\n";
    md << "| **Tested** | **" << results.size() << "** (" << sel.mainFiles << " model main files, " << sel.variants
       << " variants) |\n\n";
    md << "Not tested: " << sel.skippedNotRunnable << " main files of models whose status is not verified or modified";
    md << ", " << sel.variantsWithoutSol.size() << " variants without `.sol` baseline";
    if (!sel.variantsWithoutSol.empty()) {
        md << " (";
        for (size_t i = 0; i < sel.variantsWithoutSol.size(); ++i) md << (i ? ", " : "") << sel.variantsWithoutSol[i];
        md << ")";
    }
    md << ".\n";

    const bool anyBad = counts[1] + counts[2] > 0;
    if (anyBad) {
        md << "\n## Deviations and failures\n\n";
        md << "| Target | Model | File | Status | Verdict | Detail | Time (s) |\n|---|---|---|---|---|---|---|\n";
        for (const auto& r : results) {
            if (r.verdict != Verdict::Deviation && r.verdict != Verdict::Fail) continue;
            md << "| " << r.target.label() << " | " << r.target.name << " | " << r.target.file << " | "
               << r.target.status << " | " << verdictName(r.verdict) << " | " << summaryOf(r) << " | "
               << std::fixed << std::setprecision(2) << r.seconds << " |\n";
        }
        for (const auto& r : results) {
            if (r.verdict != Verdict::Deviation || r.comparison.details.empty()) continue;
            md << "\n### " << r.target.label() << " " << r.target.name << "\n\n";
            for (const auto& d : r.comparison.details) md << "- `" << d << "`\n";
            const size_t shown = r.comparison.details.size();
            if (r.comparison.outside + r.comparison.missing > shown)
                md << "- ... " << (r.comparison.outside + r.comparison.missing - shown) << " more\n";
        }
    }

    md << "\n## All targets\n\n";
    md << "| Target | Model | File | Status | Verdict | Detail | Time (s) |\n|---|---|---|---|---|---|---|\n";
    for (const auto& r : results) {
        md << "| " << r.target.label() << " | " << r.target.name << " | " << r.target.file << " | " << r.target.status
           << " | " << verdictName(r.verdict) << " | " << summaryOf(r) << " | " << std::fixed << std::setprecision(2)
           << r.seconds << " |\n";
    }
}

/// A fresh folder in the temporary directory, removed with the object.
struct TempDir {
    fs::path path;
    explicit TempDir(const std::string& prefix) {
        path = fs::temp_directory_path() / (prefix + std::to_string(std::time(nullptr)) + "_" +
                                            std::to_string(std::random_device{}() % 100000));
        fs::create_directories(path);
    }
    ~TempDir() { std::error_code ec; fs::remove_all(path, ec); }
};

struct LibraryLocation {
    fs::path root;
    std::string problem;      ///< why there is no library to test; empty when `root` is one
    bool requested = false;   ///< `root` comes from COOLSOLVE_LIBRARY_DIR: a bad one is an error, not a skip
};

bool isLibraryCheckout(const fs::path& dir) {
    return !dir.empty() && fs::exists(dir / "library.json") && fs::is_directory(dir / "models");
}

/// The library to test: the folder asked for (`requestedDir`, COOLSOLVE_LIBRARY_DIR), else the checkout
/// embedded in the build (`embeddedDir`).
LibraryLocation locateLibrary(const char* requestedDir, const std::string& embeddedDir) {
    LibraryLocation location;
    if (requestedDir && *requestedDir) {
        location.root = requestedDir;
        location.requested = true;
        if (!isLibraryCheckout(location.root))
            location.problem = "COOLSOLVE_LIBRARY_DIR=" + location.root.string() +
                               " is not a CoolSolve Library checkout (it needs library.json and models/)";
    } else if (isLibraryCheckout(embeddedDir)) {
        location.root = embeddedDir;
    } else {
        location.problem = "No CoolSolve Library to test: set COOLSOLVE_LIBRARY_DIR to a checkout of "
                           "https://github.com/IntSusEnergySystems/CoolSolve_Library, or build with "
                           "COOLSOLVE_EMBED_LIBRARY=ON, which tests the checkout embedded in the build (" +
                           (embeddedDir.empty() ? std::string("none in this build") : "'" + embeddedDir + "' has none") +
                           ")";
    }
    return location;
}

}  // namespace

TEST_CASE("Library regression: every runnable library model against its .sol baseline",
          "[.][library-regression]") {
    // --- Where the library is ---
    const LibraryLocation location = locateLibrary(std::getenv("COOLSOLVE_LIBRARY_DIR"), COOLSOLVE_LIBRARY_SOURCE_DIR);
    if (!location.problem.empty()) {
        if (location.requested) FAIL(location.problem);
        SKIP(location.problem);
    }
    const fs::path root = fs::weakly_canonical(location.root);

    const coolsolve::LibraryStore library{fileSystemSource(root)};
    REQUIRE(library.available());

    const Selection sel = selectTargets(root, splitList(std::getenv("COOLSOLVE_LIBRARY_TEST_FILTER")),
                                        splitList(std::getenv("COOLSOLVE_LIBRARY_TEST_EXCLUDE")));
    for (const auto& p : sel.problems) FAIL_CHECK(p);
    REQUIRE(!sel.targets.empty());

    const auto info = coolsolve::getVersionInfo();
    const std::string baselineRevision = libraryCoolPropBaseline(root);
    std::string coolpropWarning;
    if (!baselineRevision.empty() && !startsWith(info.coolpropRevision, baselineRevision) &&
        !startsWith(baselineRevision, info.coolpropRevision))
        coolpropWarning = "WARNING: CoolProp revision " + info.coolpropRevision.substr(0, 8) +
                          ", the library baselines were produced with CoolProp " + baselineRevision +
                          ": deviations may be property-library differences, not model regressions";

    std::cout << "\nLibrary regression: " << sel.targets.size() << " targets (" << sel.mainFiles
              << " model main files, " << sel.variants << " variants) of " << root.string() << "\n";
    std::cout << "CoolSolve " << info.coolsolve << ", CoolProp " << info.coolprop << " (" << info.coolpropRevision.substr(0, 8)
              << ")\n";
    if (!coolpropWarning.empty()) std::cout << coolpropWarning << "\n";
    if (!sel.restriction.empty()) std::cout << "Partial run: " << sel.restriction << "\n";
    std::cout << "(Press Ctrl+C to stop after the current target and write the report)\n\n";

    // --- Solve every target from a copy of its folder ---
    const TempDir scratch("coolsolve_library_regression_");
    const auto savedFlags = std::cout.flags();
    const auto savedPrecision = std::cout.precision();
    g_interrupted = false;
    auto previousHandler = std::signal(SIGINT, libraryRegressionSigint);

    std::vector<TargetResult> results;
    const auto runStart = std::chrono::steady_clock::now();
    for (const auto& target : sel.targets) {
        if (g_interrupted) break;
        std::cout << "[" << std::setw(3) << results.size() + 1 << "/" << sel.targets.size() << "] " << std::left
                  << std::setw(24) << target.label() << " " << std::setw(47) << target.name.substr(0, 46) << std::right << " "
                  << std::flush;   // a crash leaves the model that caused it on screen
        results.push_back(runTarget(target, scratch.path, library));
        const TargetResult& r = results.back();
        std::cout << verdictName(r.verdict) << " (" << std::fixed << std::setprecision(2) << r.seconds << " s)";
        if (r.verdict != Verdict::Pass) std::cout << " " << summaryOf(r);
        std::cout << std::defaultfloat << "\n" << std::flush;
    }
    const double totalSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - runStart).count();
    std::signal(SIGINT, previousHandler);
    const size_t notRun = sel.targets.size() - results.size();

    // --- Report ---
    fs::path reportPath = "library_regression.md";
    if (const char* env = std::getenv("COOLSOLVE_LIBRARY_REPORT"); env && *env) reportPath = env;
    reportPath = fs::absolute(reportPath);
    writeReport(reportPath, results, sel, root, library.snapshot(), coolpropWarning, totalSeconds, notRun);

    const auto counts = countVerdicts(results);
    std::cout << "\n" << results.size() << " target(s) tested: " << counts[0] << " PASS, " << counts[1] << " DEVIATION, "
              << counts[2] << " FAIL, " << counts[3] << " without baseline, in " << std::fixed << std::setprecision(0)
              << totalSeconds << " s\n";
    std::cout << "Report written to: " << reportPath.string() << "\n";
    std::cout.flags(savedFlags);
    std::cout.precision(savedPrecision);

    // --- Verdict: every target must solve and match its baseline ---
    for (const auto& r : results)
        if (r.verdict == Verdict::Deviation || r.verdict == Verdict::Fail)
            FAIL_CHECK(r.target.label() << " " << r.target.name << " (" << r.target.file << "): " << verdictName(r.verdict)
                       << ": " << summaryOf(r));
    if (notRun) FAIL_CHECK("Interrupted: " << notRun << " of " << sel.targets.size() << " target(s) were not run");
}

// ============================================================================
// Self-tests of the machinery, on a small library written to a temporary folder (fast: they run with
// the unit tests). The slow test above is what tests the real library.
// ============================================================================

namespace {

/// A library checkout of a few models: library.json (the index `$INCLUDE` looks models up in),
/// models/<name>/model.json and the files given.
struct MiniLibrary {
    TempDir dir{"coolsolve_library_regression_test_"};
    nlohmann::json index = {{"snapshot", {{"models", 0}}}, {"models", nlohmann::json::array()}};

    MiniLibrary() {
        write("library.json", index.dump());
        fs::create_directories(root() / "models");
    }
    const fs::path& root() const { return dir.path; }

    void write(const std::string& relativePath, const std::string& text) const {
        const fs::path file = root() / relativePath;
        fs::create_directories(file.parent_path());
        std::ofstream(file) << text;
    }
    void addModel(const std::string& id, const std::string& name, const std::string& status,
                  const std::map<std::string, std::string>& files) {
        const nlohmann::json meta = {{"id", id}, {"name", name}, {"status", status},
                                     {"main_file", name + ".eescode"}, {"path", "models/" + name}};
        write("models/" + name + "/model.json", meta.dump());
        for (const auto& [file, text] : files) write("models/" + name + "/" + file, text);
        index["models"].push_back(meta);
        write("library.json", index.dump());
    }
};

const std::string kSolvable = "a = 2\nb = 3 * a\n";
const std::string kSolutionOk = "a = 2.000000000000e+00 \"m\"\nB = 6.000000000000e+00\n# a comment\nlabel = 'text'\n";

}  // namespace

TEST_CASE("Library regression: locating the library, a clean skip message or an error when there is none",
          "[library-regression]") {
    const MiniLibrary lib;
    const std::string missing = (lib.root() / "does_not_exist").string();

    SECTION("nothing asked, nothing embedded: a skip message that says what to do") {
        for (const std::string& embedded : {std::string(), missing}) {
            const LibraryLocation location = locateLibrary(nullptr, embedded);
            CHECK_FALSE(location.requested);   // the test skips, it does not fail
            CHECK_THAT(location.problem, Catch::Matchers::ContainsSubstring("COOLSOLVE_LIBRARY_DIR"));
            CHECK_THAT(location.problem, Catch::Matchers::ContainsSubstring("COOLSOLVE_EMBED_LIBRARY"));
        }
        CHECK_THAT(locateLibrary("", missing).problem, Catch::Matchers::ContainsSubstring(missing));
    }
    SECTION("COOLSOLVE_LIBRARY_DIR that is not a library: an error naming the folder") {
        const LibraryLocation location = locateLibrary(missing.c_str(), lib.root().string());   // the embedded one is not a fallback
        CHECK(location.requested);
        CHECK_THAT(location.problem, Catch::Matchers::ContainsSubstring(missing));
        const TempDir empty("coolsolve_library_regression_test_");   // exists, but has no library.json
        CHECK(locateLibrary(empty.path.string().c_str(), "").requested);
        CHECK_FALSE(locateLibrary(empty.path.string().c_str(), "").problem.empty());
    }
    SECTION("a library: the folder asked for first, else the embedded one") {
        const LibraryLocation embedded = locateLibrary(nullptr, lib.root().string());
        CHECK(embedded.problem.empty());
        CHECK(embedded.root == lib.root());
        CHECK_FALSE(embedded.requested);
        const LibraryLocation asked = locateLibrary(lib.root().string().c_str(), missing);
        CHECK(asked.problem.empty());
        CHECK(asked.root == lib.root());
        CHECK(asked.requested);
    }
}

TEST_CASE("Library regression: reading a .sol baseline", "[library-regression]") {
    const TempDir dir("coolsolve_library_regression_test_");
    const fs::path sol = dir.path / "x.sol";
    std::ofstream(sol) <<
        "a = 1.5e+00 \"kW\"\n"
        "B_1[2] = -3\n"
        "s = 'text'\n"
        "# IntegralTable\n"
        "t,y\n"
        "0,1\n"
        " c = nan\n"
        "d = inf\n"
        "g = +.5\n"
        "bad = 1e\n"          // not a number
        "e = 2.0 \"unit\" x\n"   // something after the unit
        "f==3\n"
        "=4\n"
        "h = 7 \"unit\n";     // unclosed unit
    const auto values = readSol(sol);
    CHECK(values.size() == 5);
    CHECK(values.at("a") == 1.5);
    CHECK(values.at("b_1[2]") == -3.0);   // names are case-insensitive
    CHECK(std::isnan(values.at("c")));
    CHECK(std::isinf(values.at("d")));
    CHECK(values.at("g") == 0.5);

    CHECK(readSol(dir.path / "missing.sol").empty());   // no file: nothing to compare, no crash
}

TEST_CASE("Library regression: comparison rule of tools/test_models.py", "[library-regression]") {
    auto compare = [](double baseline, double fresh) {
        return compareWithBaseline({{"x", fresh}}, {{"x", baseline}});
    };
    CHECK(compare(100, 100.00005).outside == 0);   // 5e-7 relative, within 1e-6
    CHECK(compare(100, 100.0002).outside == 1);    // 2e-6 relative
    CHECK(compare(0, 5e-10).outside == 0);         // within the absolute 1e-9
    CHECK(compare(0, 5e-9).outside == 1);
    CHECK(compare(-1e5, -1e5).worstRatio == 0);
    CHECK(compare(1, std::nan("")).outside == 1);  // a value that became undefined is a deviation...
    CHECK(compare(std::nan(""), std::nan("")).outside == 0);   // ...unless it was undefined in the baseline
    CHECK(compare(INFINITY, INFINITY).outside == 0);
    CHECK(compare(INFINITY, 1).outside == 1);

    const Comparison missing = compareWithBaseline({}, {{"x", 1.0}});
    CHECK(missing.missing == 1);
    CHECK(missing.outside == 0);
    CHECK_THAT(missing.details.at(0), Catch::Matchers::ContainsSubstring("x missing"));

    // The worst variable is the one nearest to (or farthest beyond) its tolerance, the rest are not reported
    const Comparison c = compareWithBaseline({{"a", 1.0}, {"b", 2.000004}, {"c", 3.0}}, {{"a", 1.0}, {"b", 2.0}, {"c", 3.0}});
    CHECK(c.compared == 3);
    CHECK(c.outside == 1);
    CHECK(c.worstName == "b");
    CHECK_THAT(c.worstRelative, Catch::Matchers::WithinRel(2e-6, 0.01));
    // Variables of the solution that are not in the baseline are not compared
    CHECK(compareWithBaseline({{"a", 1.0}, {"extra", 9.0}}, {{"a", 1.0}}).compared == 1);
}

TEST_CASE("Library regression: targets of a library and their verdicts", "[library-regression]") {
    MiniLibrary lib;
    lib.addModel("CSL-9001", "pass_model", "verified", {{"pass_model.eescode", kSolvable}, {"pass_model.sol", kSolutionOk}});
    lib.addModel("CSL-9002", "deviation_model", "modified",
                 {{"deviation_model.eescode", kSolvable}, {"deviation_model.sol", "b = 6.1\n"}});
    lib.addModel("CSL-9003", "missing_variable_model", "verified",
                 {{"missing_variable_model.eescode", kSolvable}, {"missing_variable_model.sol", "b = 6\nzzz = 1\n"}});
    lib.addModel("CSL-9004", "parse_error_model", "verified",
                 {{"parse_error_model.eescode", "a = = 2\n"}, {"parse_error_model.sol", "a = 2\n"}});
    lib.addModel("CSL-9005", "solver_error_model", "verified",
                 {{"solver_error_model.eescode", "x^2 = -1\n"}, {"solver_error_model.sol", "x = 1\n"}});
    lib.addModel("CSL-9006", "malformed_sol_model", "verified",
                 {{"malformed_sol_model.eescode", kSolvable}, {"malformed_sol_model.sol", "not a solution\nb = six\n"}});
    lib.addModel("CSL-9007", "no_sol_model", "verified", {{"no_sol_model.eescode", kSolvable}});
    // A blocked model is not tested, but its runnable variant is, when it has a baseline
    lib.addModel("CSL-9008", "blocked_model", "blocked",
                 {{"blocked_model.eescode", "a = = 2\n"},
                  {"blocked_model_coolsolve.eescode", kSolvable}, {"blocked_model_coolsolve.sol", kSolutionOk},
                  {"blocked_model_simplified.eescode", kSolvable}});
    // $INCLUDE of a library model resolves in the library under test, not in the one embedded in the binary
    lib.addModel("CSL-9009", "twice_lib", "verified",
                 {{"twice_lib.eescode", "FUNCTION twice(x)\n  twice := 2*x\nEND\ny = twice(3)\n"}, {"twice_lib.sol", "y = 6\n"}});
    lib.addModel("CSL-9010", "include_model", "verified",
                 {{"include_model.eescode", "$INCLUDE twice_lib\nz = twice(4)\n"}, {"include_model.sol", "z = 8\n"}});
    // A model.json that cannot be read is reported, and does not stop the others
    lib.write("models/broken/model.json", "{ not json");

    SECTION("selection") {
        const Selection sel = selectTargets(lib.root(), {}, {});
        CHECK(sel.mainFiles == 9);
        CHECK(sel.variants == 1);
        CHECK(sel.skippedNotRunnable == 1);
        CHECK(sel.variantsWithoutSol == std::vector<std::string>{"CSL-9008:simplified"});
        REQUIRE(sel.problems.size() == 1);
        CHECK_THAT(sel.problems[0], Catch::Matchers::ContainsSubstring("model.json"));

        const Selection filtered = selectTargets(lib.root(), {"CSL-9001", "CSL-9008:coolsolve", "CSL-0000"}, {});
        CHECK(filtered.targets.size() == 2);
        REQUIRE(filtered.problems.size() == 2);   // the broken model.json, the unknown target
        CHECK_THAT(filtered.problems.back(), Catch::Matchers::ContainsSubstring("CSL-0000"));

        const Selection excluded = selectTargets(lib.root(), {}, {"CSL-9001", "CSL-9008:coolsolve"});
        CHECK(excluded.targets.size() == sel.targets.size() - 2);

        CHECK(sel.restriction.empty());   // a partial run says so in its report
        CHECK_THAT(filtered.restriction, Catch::Matchers::ContainsSubstring("only CSL-0000,CSL-9001,CSL-9008:coolsolve"));
        CHECK_THAT(excluded.restriction, Catch::Matchers::ContainsSubstring("without CSL-9001"));
    }

    SECTION("verdicts") {
        const Selection sel = selectTargets(lib.root(), {}, {});
        const coolsolve::LibraryStore library{fileSystemSource(lib.root())};
        REQUIRE(library.available());
        const TempDir scratch("coolsolve_library_regression_test_");
        std::map<std::string, TargetResult> results;
        for (const auto& target : sel.targets) results[target.label()] = runTarget(target, scratch.path, library);
        REQUIRE(results.size() == 10);

        CHECK(results["CSL-9001"].verdict == Verdict::Pass);   // case-insensitive names, string and comment lines ignored
        CHECK(results["CSL-9001"].comparison.compared == 2);

        CHECK(results["CSL-9002"].verdict == Verdict::Deviation);
        CHECK(results["CSL-9002"].comparison.outside == 1);
        CHECK(results["CSL-9002"].comparison.worstName == "b");
        CHECK_THAT(summaryOf(results["CSL-9002"]), Catch::Matchers::ContainsSubstring("worst b: 6.1 -> 6"));

        CHECK(results["CSL-9003"].verdict == Verdict::Deviation);
        CHECK(results["CSL-9003"].comparison.missing == 1);
        CHECK_THAT(results["CSL-9003"].comparison.details.at(0), Catch::Matchers::ContainsSubstring("zzz missing"));

        CHECK(results["CSL-9004"].verdict == Verdict::Fail);
        CHECK_THAT(results["CSL-9004"].message, Catch::Matchers::ContainsSubstring("Parse failed: line 1"));

        CHECK(results["CSL-9005"].verdict == Verdict::Fail);
        CHECK_THAT(results["CSL-9005"].message, Catch::Matchers::ContainsSubstring("Solver SingularJacobian"));

        CHECK(results["CSL-9006"].verdict == Verdict::Fail);   // a malformed baseline is reported, it does not pass
        CHECK_THAT(results["CSL-9006"].message, Catch::Matchers::ContainsSubstring("malformed_sol_model.sol"));
        CHECK_THAT(results["CSL-9006"].message, Catch::Matchers::ContainsSubstring("no numeric value"));

        CHECK(results["CSL-9007"].verdict == Verdict::NoBaseline);   // a main file without .sol is only solved

        CHECK(results["CSL-9008:coolsolve"].verdict == Verdict::Pass);

        CHECK(results["CSL-9009"].verdict == Verdict::Pass);
        CHECK(results["CSL-9010"].verdict == Verdict::Pass);

        // The scratch copies are removed
        CHECK(fs::is_empty(scratch.path));

        // The report names every target and its verdict
        std::vector<TargetResult> all;
        for (const auto& [label, r] : results) all.push_back(r);
        const fs::path report = scratch.path / "report.md";
        writeReport(report, all, sel, lib.root(), library.snapshot(), "", 1.0, 0);
        std::ifstream in(report);
        const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        CHECK_THAT(text, Catch::Matchers::ContainsSubstring("| DEVIATION | 2 |"));
        CHECK_THAT(text, Catch::Matchers::ContainsSubstring("| FAIL | 3 |"));
        CHECK_THAT(text, Catch::Matchers::ContainsSubstring("| CSL-9005 | solver_error_model |"));
    }
}
