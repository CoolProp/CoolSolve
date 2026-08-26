#include "coolsolve/usage_log.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <unordered_map>

namespace fs = std::filesystem;

namespace coolsolve {

// ============================================================================
// Constants
// ============================================================================

/// Header line written once when the log file is created.
static const char* kLogHeader =
    "# timestamp,kind,outcome,duration_ms,model_bytes,equations,max_block,model,ip,version";

/// Number of logarithmic duration bins (edges: 1 ms .. 2^25 ms).
static constexpr int kHistogramBins = 25;
static constexpr double kHistogramBaseMs = 1.0;

/// Maximum number of distinct names tracked before pruning (memory bound).
static constexpr std::size_t kNameMapCap = 20000;
static constexpr std::size_t kNameKeepAfterPrune = 5000;

/// Number of entries kept in the top rankings.
static constexpr std::size_t kTopEntries = 10;

// ============================================================================
// Helpers
// ============================================================================

/// Current UTC time formatted as ISO 8601 ("2026-08-26T14:03:11Z").
static std::string iso8601UtcNow() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tmUtc{};
#if defined(_MSC_VER)
    gmtime_s(&tmUtc, &t);
#else
    gmtime_r(&t, &tmUtc);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tmUtc);
    return buf;
}

/// Quote a CSV field if it contains characters that would break the format.
/// Newlines and other control characters are replaced by a space first: the
/// log is read back line by line, so a record must never span two lines
/// (quoting alone would not help a line-oriented reader).
static std::string csvField(const std::string& s) {
    bool needsQuote = false;
    bool needsSanitize = false;
    for (unsigned char c : s) {
        if (c == ',' || c == '"') needsQuote = true;
        else if (c < 0x20 || c == 0x7f) needsSanitize = true;
    }
    if (!needsQuote && !needsSanitize) return s;

    std::string cleaned;
    cleaned.reserve(s.size());
    for (unsigned char c : s) cleaned += (c < 0x20 || c == 0x7f) ? ' ' : static_cast<char>(c);
    if (!needsQuote) return cleaned;

    std::string out = "\"";
    for (char c : cleaned) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += '"';
    return out;
}

// ============================================================================
// Writing
// ============================================================================

void appendUsageLog(const std::string& path, const UsageLogEntry& entry) {
    // Serialize outside the lock: formatting is pure work on local data.
    // Fixed notation with 3 decimals: the default float format switches to
    // scientific past 6 significant digits, which is unreadable in a CSV and
    // loses sub-ms resolution for solves longer than ~17 minutes.
    std::ostringstream line;
    line << iso8601UtcNow() << ','
         << entry.kind << ','
         << entry.outcome << ','
         << std::fixed << std::setprecision(3) << entry.durationMs << std::defaultfloat << ','
         << entry.modelBytes << ','
         << entry.equations << ','
         << entry.maxBlockSize << ','
         << csvField(entry.modelName.empty() ? std::string("(unnamed)") : entry.modelName) << ','
         << csvField(entry.clientIp) << ','
         << csvField(entry.version)
         << '\n';

    // One global mutex: appends are rare (once per solve attempt) and cheap,
    // so a simple lock is preferable to per-file handles or an async queue.
    static std::mutex appendMutex;
    std::lock_guard<std::mutex> lock(appendMutex);

    // Create the file with its header on first use, then always append.
    // Errors are intentionally swallowed: logging must never break solving.
    std::error_code ec;
    if (!fs::exists(path, ec)) {
        std::ofstream create(path, std::ios::binary | std::ios::trunc);
        if (create.is_open()) create << kLogHeader << '\n';
    }
    std::ofstream file(path, std::ios::binary | std::ios::app);
    if (file.is_open()) file << line.str();
}

// ============================================================================
// Reading / aggregation
// ============================================================================

namespace {

using NameCountMap = std::unordered_map<std::string, long>;

/// Keep the map bounded: when it grows past the cap, retain only the most
/// frequent names. Rankings remain exact unless a single file holds more
/// distinct names than the cap.
void pruneNameMap(NameCountMap& map) {
    if (map.size() <= kNameMapCap) return;
    std::vector<std::pair<long, const std::string*>> items;
    items.reserve(map.size());
    for (const auto& [name, count] : map) items.emplace_back(count, &name);
    std::nth_element(items.begin(), items.begin() + kNameKeepAfterPrune, items.end(),
                     [](const auto& a, const auto& b) { return a.first > b.first; });
    NameCountMap kept;
    kept.reserve(kNameKeepAfterPrune);
    for (std::size_t i = 0; i < kNameKeepAfterPrune && i < items.size(); ++i)
        kept[*items[i].second] = items[i].first;
    map = std::move(kept);
}

/// Convert a NameCountMap to a descending-sorted vector of at most kTopEntries.
std::vector<UsageNameCount> topOf(const NameCountMap& map) {
    std::vector<UsageNameCount> all;
    all.reserve(map.size());
    for (const auto& [name, count] : map) all.push_back({name, count});
    std::partial_sort(all.begin(),
                      all.begin() + std::min(all.size(), kTopEntries),
                      all.end(),
                      [](const UsageNameCount& a, const UsageNameCount& b) {
                          if (a.count != b.count) return a.count > b.count;
                          return a.name < b.name;
                      });
    if (all.size() > kTopEntries) all.resize(kTopEntries);
    return all;
}

/// Split the leading fixed-format fields of a CSV log line
/// (timestamp..max_block). Returns false when the line is structurally invalid.
bool parseFixedFields(const std::string& line, std::string fields[8]) {
    int idx = 0;
    std::size_t start = 0;
    for (; idx < 7; ++idx) {
        std::size_t comma = line.find(',', start);
        if (comma == std::string::npos) return false;
        fields[idx] = line.substr(start, comma - start);
        start = comma + 1;
    }
    fields[7] = line.substr(start);   // tail: model,ip,version (parsed later)
    return true;
}

/// Parse the trailing "model,ip,version" section, honouring CSV quoting.
bool parseTailFields(const std::string& tail, std::string& model,
                     std::string& ip, std::string& version) {
    std::size_t pos = 0;
    if (pos < tail.size() && tail[pos] == '"') {
        std::string unescaped;
        bool closed = false;
        ++pos;
        while (pos < tail.size()) {
            if (tail[pos] == '"') {
                if (pos + 1 < tail.size() && tail[pos + 1] == '"') {
                    unescaped += '"';
                    pos += 2;
                } else { ++pos; closed = true; break; }
            } else {
                unescaped += tail[pos++];
            }
        }
        // An unterminated quote, or anything other than the field separator
        // right after the closing quote, means the record is corrupt.
        if (!closed || pos >= tail.size() || tail[pos] != ',') return false;
        ++pos;   // consume the separator
        model = unescaped;
    } else {
        std::size_t comma = tail.find(',', pos);
        if (comma == std::string::npos) return false;
        model = tail.substr(pos, comma - pos);
        pos = comma + 1;
    }
    std::size_t comma = tail.find(',', pos);
    if (comma == std::string::npos) return false;
    ip = tail.substr(pos, comma - pos);
    version = tail.substr(comma + 1);
    return true;
}

} // anonymous namespace

UsageLogStats computeUsageLogStats(const std::string& path) {
    // --- Cache: recompute only when the file changed ------------------------
    struct CacheEntry {
        std::mutex mutex;
        bool valid = false;
        std::string path;
        std::uintmax_t mtime = 0;
        std::uintmax_t size = 0;
        UsageLogStats stats;
    };
    static CacheEntry cache;

    std::error_code ec;
    fs::path filePath(path);
    if (!fs::exists(filePath, ec)) {
        UsageLogStats missing;
        missing.valid = false;
        return missing;
    }
    auto times = fs::last_write_time(filePath, ec);
    std::uintmax_t mtime = ec ? 0 : times.time_since_epoch().count();
    std::uintmax_t size = fs::file_size(filePath, ec);
    if (ec) size = 0;

    {
        std::lock_guard<std::mutex> lock(cache.mutex);
        if (cache.valid && cache.path == path && cache.mtime == mtime && cache.size == size) {
            return cache.stats;
        }
    }

    // --- Single streaming pass, bounded memory ------------------------------
    UsageLogStats stats;
    stats.histogramEdges.resize(kHistogramBins + 1);
    for (int i = 0; i <= kHistogramBins; ++i)
        stats.histogramEdges[i] = kHistogramBaseMs * static_cast<double>(1ull << i);
    stats.histogramCounts.assign(kHistogramBins, 0);

    NameCountMap models;
    NameCountMap ips;

    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return stats;   // valid stays false

    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;

        std::string fixed[8];
        if (!parseFixedFields(line, fixed)) { stats.malformedLines++; continue; }

        const std::string& timestamp = fixed[0];
        const std::string& kind = fixed[1];
        const std::string& outcome = fixed[2];

        char* end = nullptr;
        double durationMs = std::strtod(fixed[3].c_str(), &end);
        if (end == fixed[3].c_str()) { stats.malformedLines++; continue; }

        std::string model, ip, version;
        if (!parseTailFields(fixed[7], model, ip, version)) { stats.malformedLines++; continue; }

        // ---- Aggregate ----
        stats.totalAttempts++;
        stats.totalDurationMs += durationMs;
        if (kind == "tryharder") stats.tryHarderAttempts++;

        bool success = (outcome == "success");
        if (success) stats.successes++;
        else if (outcome == "parse_error") stats.parseErrors++;
        else stats.failures++;

        if (timestamp.size() >= 10) {
            UsageDailyCount& day = stats.daily[timestamp.substr(0, 10)];
            day.attempts++;
            if (success) day.successes++;
        }

        // Duration histogram (clamped at both ends)
        double d = std::max(durationMs, 0.0);
        int bin = kHistogramBins - 1;
        if (d < stats.histogramEdges[kHistogramBins]) {
            bin = static_cast<int>(std::floor(std::log2(std::max(d, 1e-9) / kHistogramBaseMs)));
            if (bin < 0) bin = 0;
            if (bin >= kHistogramBins) bin = kHistogramBins - 1;
        }
        stats.histogramCounts[bin]++;

        if (!model.empty()) { models[model]++; pruneNameMap(models); }
        if (!ip.empty()) { ips[ip]++; pruneNameMap(ips); }
    }

    stats.topModels = topOf(models);
    stats.topIps = topOf(ips);
    stats.uniqueIps = static_cast<long>(ips.size());

    // Median / p95 from the histogram (geometric interpolation inside bins).
    auto percentile = [&](double q) -> double {
        long total = 0;
        for (long c : stats.histogramCounts) total += c;
        if (total == 0) return 0.0;
        double target = q * static_cast<double>(total);
        long cumulative = 0;
        for (int i = 0; i < kHistogramBins; ++i) {
            long after = cumulative + stats.histogramCounts[i];
            if (static_cast<double>(after) >= target && stats.histogramCounts[i] > 0) {
                double lo = stats.histogramEdges[i];
                double hi = stats.histogramEdges[i + 1];
                double frac = (target - cumulative) / static_cast<double>(stats.histogramCounts[i]);
                return lo * std::pow(hi / lo, frac);
            }
            cumulative = after;
        }
        return stats.histogramEdges[kHistogramBins];
    };
    stats.medianMs = percentile(0.5);
    stats.p95Ms = percentile(0.95);
    stats.valid = true;

    {
        std::lock_guard<std::mutex> lock(cache.mutex);
        cache.valid = true;
        cache.path = path;
        cache.mtime = mtime;
        cache.size = size;
        cache.stats = stats;
    }
    return stats;
}

} // namespace coolsolve
