/**
 * test_usage_log.cpp — Unit tests for the GUI solve-attempt usage log.
 *
 * Covers appending entries (CSV round-trip, escaping, header creation) and
 * the aggregated statistics (totals, daily counts, histogram, rankings,
 * percentile sanity, malformed-line handling and cache invalidation).
 *
 * Run with: ./coolsolve_tests "[usage-log]"
 */

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "coolsolve/usage_log.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace coolsolve;

// ============================================================================
// Helpers
// ============================================================================

/// Unique temporary path for one test run
static std::string tempLogPath(const std::string& tag) {
    static std::uint64_t counter = 0;
    static const auto seed = std::random_device{}();
    std::ostringstream name;
    name << "coolsolve_test_usage_" << tag << "_" << seed << "_" << counter++ << ".log";
    auto path = fs::temp_directory_path() / name.str();
    fs::remove(path);
    return path.string();
}

static UsageLogEntry makeEntry(const std::string& outcome = "success",
                               double durationMs = 100.0,
                               const std::string& model = "model_a",
                               const std::string& ip = "1.2.3.4",
                               const std::string& kind = "solve") {
    UsageLogEntry e;
    e.kind = kind;
    e.outcome = outcome;
    e.durationMs = durationMs;
    e.modelBytes = 500;
    e.equations = 10;
    e.maxBlockSize = 4;
    e.modelName = model;
    e.clientIp = ip;
    e.version = "0.3";
    return e;
}

/// Read the raw file content (for header/format checks)
static std::string readFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// ============================================================================
// Writing
// ============================================================================

TEST_CASE("Usage log creates file with single header", "[usage-log]") {
    auto path = tempLogPath("header");
    appendUsageLog(path, makeEntry());
    appendUsageLog(path, makeEntry("failed"));

    const std::string content = readFile(path);
    REQUIRE(content.rfind("# timestamp,", 0) == 0);

    // Exactly one header line
    size_t headers = 0;
    std::istringstream ss(content);
    std::string line;
    while (std::getline(ss, line)) {
        if (!line.empty() && line[0] == '#') headers++;
    }
    REQUIRE(headers == 1);

    // Two data lines
    size_t dataLines = 0;
    ss.clear();
    ss.str(content);
    while (std::getline(ss, line)) {
        if (!line.empty() && line[0] != '#') dataLines++;
    }
    REQUIRE(dataLines == 2);
    std::remove(path.c_str());
}

TEST_CASE("Usage log escapes commas and quotes in names", "[usage-log]") {
    auto path = tempLogPath("escape");
    appendUsageLog(path, makeEntry("success", 100.0, "heat, \"pump\" v2"));

    const std::string content = readFile(path);
    REQUIRE(content.find("\"heat, \"\"pump\"\" v2\"") != std::string::npos);

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.topModels.size() == 1);
    REQUIRE(stats.topModels[0].name == "heat, \"pump\" v2");
    std::remove(path.c_str());
}

TEST_CASE("Quoted model name does not shift the remaining fields", "[usage-log]") {
    // Regression: the field separator following a quoted model name must be
    // consumed, otherwise the ip/version columns are silently lost.
    auto path = tempLogPath("quoted-shift");
    appendUsageLog(path, makeEntry("success", 100.0, "heat, pump", "203.0.113.7"));
    appendUsageLog(path, makeEntry("success", 100.0, "plain", "203.0.113.7"));

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.malformedLines == 0);
    REQUIRE(stats.uniqueIps == 1);
    REQUIRE(stats.topIps.size() == 1);
    REQUIRE(stats.topIps[0].name == "203.0.113.7");
    REQUIRE(stats.topIps[0].count == 2);   // both lines, quoted one included
    std::remove(path.c_str());
}

TEST_CASE("Model names with newlines stay on one line", "[usage-log]") {
    // The log is read back line by line, so control characters must never
    // reach the file — quoting alone would split the record in two.
    auto path = tempLogPath("newline");
    appendUsageLog(path, makeEntry("success", 100.0, "line1\nline2\ttab"));

    const std::string content = readFile(path);
    size_t lines = 0;
    std::istringstream ss(content);
    std::string line;
    while (std::getline(ss, line)) lines++;
    REQUIRE(lines == 2);   // header + exactly one record

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.totalAttempts == 1);
    REQUIRE(stats.malformedLines == 0);
    REQUIRE(stats.topModels.size() == 1);
    REQUIRE(stats.topModels[0].name == "line1 line2 tab");
    std::remove(path.c_str());
}

TEST_CASE("Long durations are written in fixed notation", "[usage-log]") {
    // The default float format switches to scientific past 6 significant
    // digits, which is unreadable in a CSV and loses sub-ms resolution.
    auto path = tempLogPath("bignum");
    appendUsageLog(path, makeEntry("success", 1234567.89));

    const std::string content = readFile(path);
    REQUIRE(content.find("e+") == std::string::npos);   // no scientific notation
    REQUIRE(content.find("1234567.890") != std::string::npos);

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.totalDurationMs == Catch::Approx(1234567.89).epsilon(1e-9));
    std::remove(path.c_str());
}

TEST_CASE("Empty model name becomes (unnamed)", "[usage-log]") {
    auto path = tempLogPath("unnamed");
    UsageLogEntry e = makeEntry();
    e.modelName.clear();
    appendUsageLog(path, e);

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.totalAttempts == 1);
    REQUIRE(stats.topModels.size() == 1);
    REQUIRE(stats.topModels[0].name == "(unnamed)");
    std::remove(path.c_str());
}

// ============================================================================
// Statistics
// ============================================================================

TEST_CASE("Statistics aggregate totals, kinds and outcomes", "[usage-log]") {
    auto path = tempLogPath("stats");

    appendUsageLog(path, makeEntry("success", 100.0));
    appendUsageLog(path, makeEntry("success", 300.0));
    appendUsageLog(path, makeEntry("failed", 200.0));
    appendUsageLog(path, makeEntry("parse_error", 5.0));
    appendUsageLog(path, makeEntry("success", 1000.0, "m", "9.9.9.9", "tryharder"));

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.totalAttempts == 5);
    REQUIRE(stats.successes == 3);
    REQUIRE(stats.failures == 1);
    REQUIRE(stats.parseErrors == 1);
    REQUIRE(stats.tryHarderAttempts == 1);
    REQUIRE(stats.uniqueIps == 2);
    REQUIRE(stats.daily.size() >= 1);

    // Today's bucket contains all attempts
    bool foundToday = false;
    for (const auto& [date, count] : stats.daily) {
        if (count.attempts > 0 && date.size() == 10) {
            foundToday = true;
            break;
        }
    }
    REQUIRE(foundToday);
    std::remove(path.c_str());
}

TEST_CASE("Duration histogram clamps into fixed bins", "[usage-log]") {
    auto path = tempLogPath("hist");
    appendUsageLog(path, makeEntry("success", 0.01));      // below first edge
    appendUsageLog(path, makeEntry("success", 100.0));     // middle bin
    appendUsageLog(path, makeEntry("success", 1e12));      // above last edge

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);

    long total = 0;
    for (long c : stats.histogramCounts) total += c;
    REQUIRE(total == 3);
    REQUIRE(stats.histogramEdges.front() == 1.0);   // first edge = 1 ms
    REQUIRE(stats.histogramCounts.front() >= 1);    // tiny duration clamped low
    REQUIRE(stats.histogramCounts.back() >= 1);     // huge duration clamped high

    // Percentiles must lie within the observed range of magnitudes
    REQUIRE(stats.medianMs >= 1.0);
    REQUIRE(stats.p95Ms <= stats.histogramEdges.back() * 2.0);
    std::remove(path.c_str());
}

TEST_CASE("Rankings are sorted descending and capped", "[usage-log]") {
    auto path = tempLogPath("rank");
    for (int i = 0; i < 5; i++) appendUsageLog(path, makeEntry("success", 50.0, "popular"));
    for (int i = 0; i < 3; i++) appendUsageLog(path, makeEntry("failed", 60.0, "medium"));
    appendUsageLog(path, makeEntry("success", 70.0, "rare"));

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.topModels.size() == 3);
    REQUIRE(stats.topModels[0].name == "popular");
    REQUIRE(stats.topModels[0].count == 5);
    REQUIRE(stats.topModels[1].name == "medium");
    REQUIRE(stats.topModels[2].name == "rare");
    std::remove(path.c_str());
}

TEST_CASE("Malformed lines are skipped and counted", "[usage-log]") {
    auto path = tempLogPath("malformed");
    appendUsageLog(path, makeEntry());

    // Append garbage directly (bypassing the writer API)
    {
        std::ofstream f(path, std::ios::app);
        f << "not,a,valid,line\n";
        f << "\n";                       // empty lines ignored, not counted
        f << "# a comment line\n";       // comments ignored, not counted
        f << "2026-08-26T10:00:00Z,solve,success,abc,500,10,4,m,ip,v\n";   // bad number
        f << "2026-08-26T10:00:00Z,solve,success\n";                       // too few fields
    }
    appendUsageLog(path, makeEntry());

    auto stats = computeUsageLogStats(path);
    REQUIRE(stats.valid);
    REQUIRE(stats.totalAttempts == 2);       // only the two valid entries
    REQUIRE(stats.malformedLines == 3);      // the three structurally bad ones
    std::remove(path.c_str());
}

TEST_CASE("Missing file yields invalid statistics", "[usage-log]") {
    auto stats = computeUsageLogStats("/nonexistent/coolsolve/no_such_log.log");
    REQUIRE_FALSE(stats.valid);
    REQUIRE(stats.totalAttempts == 0);
}

TEST_CASE("Cache invalidates when the file grows", "[usage-log]") {
    auto path = tempLogPath("cache");
    appendUsageLog(path, makeEntry());

    auto first = computeUsageLogStats(path);
    REQUIRE(first.totalAttempts == 1);

    // Cached second call must reflect the same state
    auto cached = computeUsageLogStats(path);
    REQUIRE(cached.totalAttempts == 1);

    // Growing the file must invalidate the cache (size change)
    appendUsageLog(path, makeEntry("failed"));
    auto updated = computeUsageLogStats(path);
    REQUIRE(updated.totalAttempts == 2);
    REQUIRE(updated.failures == 1);
    std::remove(path.c_str());
}
