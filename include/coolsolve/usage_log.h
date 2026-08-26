#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace coolsolve {

// ============================================================================
// GUI usage log
//
// The embedded HTTP server records one line per solve attempt ("Solve" or
// "Try Harder") in a small CSV file. The log lives next to the server's
// working directory (default `coolsolve_gui.log`, see ServerOptions /
// COOLSOLVE_GUI_LOG), is git-ignored, and is written *after* the solve has
// completed, off the solver's critical path, so it adds no measurable
// overhead. CLI-only runs are never logged.
//
// File format (one header line starting with '#', then one CSV line per
// attempt):
//
//   # timestamp,kind,outcome,duration_ms,model_bytes,equations,max_block,model,ip,version
//   2026-08-26T14:03:11Z,solve,success,1523.4,412,8,4,"heat_pump.eescode",203.0.113.7,0.3
// ============================================================================

/// One solve attempt as recorded in the GUI usage log.
struct UsageLogEntry {
    std::string kind;           ///< "solve" or "tryharder"
    std::string outcome;        ///< "success", "failed" or "parse_error"
    double durationMs = 0.0;    ///< Wall-clock duration of the attempt [ms]
    std::uint64_t modelBytes = 0;   ///< Size of the .eescode source [bytes]
    int equations = -1;         ///< Number of equations (-1 = unknown, e.g. parse error)
    int maxBlockSize = -1;      ///< Largest solution block (-1 = unknown)
    std::string modelName;      ///< User-facing model name ("" -> "(unnamed)")
    std::string clientIp;       ///< Requesting client IP (X-Forwarded-For aware)
    std::string version;        ///< CoolSolve version string
};

/// Append one CSV line to the log file at `path`. Creates the file and its
/// header line on first use. Best-effort: I/O errors are silently ignored so
/// that logging can never break solving.
void appendUsageLog(const std::string& path, const UsageLogEntry& entry);

/// Per-day attempt counters (key: "YYYY-MM-DD", UTC).
struct UsageDailyCount {
    long attempts = 0;
    long successes = 0;
};

/// A "name → count" pair used for top-model / top-IP rankings.
struct UsageNameCount {
    std::string name;
    long count = 0;
};

/// Aggregated usage statistics. Computed in a single streaming pass over the
/// log file with bounded memory: per-day counts grow with time (one entry per
/// day), while model/IP maps are pruned to a fixed cap and durations are
/// aggregated into a fixed-size logarithmic histogram (used for the median /
/// p95 estimates).
struct UsageLogStats {
    bool valid = false;             ///< false when the file is missing/unreadable
    long totalAttempts = 0;
    long successes = 0;
    long failures = 0;
    long parseErrors = 0;
    long tryHarderAttempts = 0;
    double totalDurationMs = 0.0;

    /// Log-spaced duration histogram (fixed size): edges[i]..edges[i+1] ms.
    std::vector<double> histogramEdges;
    std::vector<long> histogramCounts;

    std::map<std::string, UsageDailyCount> daily;   ///< Sorted by date

    std::vector<UsageNameCount> topModels;  ///< Descending, capped (10)
    std::vector<UsageNameCount> topIps;     ///< Descending, capped (10)
    /// Distinct client IPs seen. Exact while the log holds no more than the
    /// internal name cap (20 000) of them; beyond that the map is pruned and
    /// this becomes a lower bound.
    long uniqueIps = 0;

    long malformedLines = 0;        ///< Lines that could not be parsed

    /// Median / 95th-percentile duration estimated from the histogram [ms].
    double medianMs = 0.0;
    double p95Ms = 0.0;
};

/// Compute aggregated statistics for the log file at `path`. Results are
/// cached and only recomputed when the file's size or modification time
/// changes, so repeated requests stay O(1) even as the log grows.
UsageLogStats computeUsageLogStats(const std::string& path);

} // namespace coolsolve
