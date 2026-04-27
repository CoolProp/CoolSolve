#include "coolsolve/lookup_table.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace coolsolve {

namespace fs = std::filesystem;

// ============================================================================
// Internal CSV helpers
// ============================================================================

static std::string trimWhitespace(const std::string& s) {
    const char* ws = " \t\r\n";
    size_t a = s.find_first_not_of(ws);
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

// Parse one CSV row into fields.  Handles double-quoted fields (including
// embedded commas and doubled quotes).
static std::vector<std::string> parseCSVRow(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    bool inQuote = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (inQuote) {
            if (c == '"') {
                // Doubled quote → literal "
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    cur += '"';
                    ++i;
                } else {
                    inQuote = false;
                }
            } else {
                cur += c;
            }
        } else {
            if (c == '"') {
                inQuote = true;
            } else if (c == ',') {
                fields.push_back(cur);
                cur.clear();
            } else {
                cur += c;
            }
        }
    }
    fields.push_back(cur);
    return fields;
}

// ============================================================================
// LookupTable constructors
// ============================================================================

LookupTable::LookupTable(std::string name,
                         std::vector<std::string> columnNames,
                         std::vector<std::vector<double>> data)
    : name_(std::move(name))
    , columnNames_(std::move(columnNames))
    , data_(std::move(data))
{}

// ============================================================================
// LookupTable::fromCSV
// ============================================================================

LookupTable LookupTable::fromCSV(const std::string& name, const std::string& csv) {
    std::istringstream stream(csv);
    std::string line;

    // Skip empty lines at the start; first non-empty line is the header
    std::vector<std::string> header;
    while (std::getline(stream, line)) {
        std::string trimmed = trimWhitespace(line);
        if (trimmed.empty()) continue;
        auto fields = parseCSVRow(trimmed);
        for (auto& f : fields) f = trimWhitespace(f);
        header = std::move(fields);
        break;
    }
    if (header.empty()) {
        throw std::runtime_error("LookupTable '" + name + "': CSV has no header row");
    }

    std::vector<std::vector<double>> data;
    size_t numCols = header.size();

    while (std::getline(stream, line)) {
        std::string trimmed = trimWhitespace(line);
        if (trimmed.empty()) continue;

        auto fields = parseCSVRow(trimmed);
        std::vector<double> row(numCols, std::numeric_limits<double>::quiet_NaN());
        for (size_t c = 0; c < numCols && c < fields.size(); ++c) {
            std::string s = trimWhitespace(fields[c]);
            if (s.empty()) {
                row[c] = std::numeric_limits<double>::quiet_NaN();
            } else {
                try {
                    size_t pos = 0;
                    row[c] = std::stod(s, &pos);
                    // Accept trailing whitespace but not other characters
                    if (pos < s.size() && !std::all_of(s.begin() + pos, s.end(),
                            [](char ch){ return ch == ' ' || ch == '\t'; })) {
                        row[c] = std::numeric_limits<double>::quiet_NaN();
                    }
                } catch (...) {
                    row[c] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
        data.push_back(std::move(row));
    }

    return LookupTable(name, std::move(header), std::move(data));
}

// ============================================================================
// LookupTable::toCSV
// ============================================================================

std::string LookupTable::toCSV() const {
    std::ostringstream out;
    // Header row
    for (size_t c = 0; c < columnNames_.size(); ++c) {
        if (c > 0) out << ',';
        out << columnNames_[c];
    }
    out << '\n';
    // Data rows
    for (const auto& row : data_) {
        for (size_t c = 0; c < row.size(); ++c) {
            if (c > 0) out << ',';
            if (std::isnan(row[c])) {
                // Empty cell
            } else {
                out << row[c];
            }
        }
        out << '\n';
    }
    return out.str();
}

// ============================================================================
// LookupTable metadata helpers
// ============================================================================

size_t LookupTable::columnIndex(const std::string& colName) const {
    std::string lower = colName;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    for (size_t i = 0; i < columnNames_.size(); ++i) {
        std::string cn = columnNames_[i];
        std::transform(cn.begin(), cn.end(), cn.begin(), ::tolower);
        if (cn == lower) return i + 1; // 1-based
    }
    return 0;
}

double LookupTable::value(size_t row1, size_t col1) const {
    if (row1 < 1 || row1 > data_.size()) {
        throw std::out_of_range("LookupTable '" + name_ + "': row " +
                                std::to_string(row1) + " out of range [1," +
                                std::to_string(data_.size()) + "]");
    }
    if (col1 < 1 || col1 > columnNames_.size()) {
        throw std::out_of_range("LookupTable '" + name_ + "': column " +
                                std::to_string(col1) + " out of range [1," +
                                std::to_string(columnNames_.size()) + "]");
    }
    return data_[row1 - 1][col1 - 1];
}

// ============================================================================
// LookupTable::findInterval
// ============================================================================

size_t LookupTable::findInterval(size_t col0, double x) const {
    // Returns lo such that data_[lo][col0] <= x < data_[lo+1][col0].
    // Clamps to [0, n-2] so interpolation is always possible when n >= 2.
    size_t n = data_.size();
    if (n < 2) return 0;

    // Binary search
    size_t lo = 0, hi = n - 1;
    while (lo + 1 < hi) {
        size_t mid = (lo + hi) / 2;
        if (data_[mid][col0] <= x) lo = mid;
        else hi = mid;
    }
    return lo;
}

// ============================================================================
// LookupTable::interpolate1D
// ============================================================================

ADValue LookupTable::interpolate1D(size_t xcol1, size_t ycol1,
                                    const ADValue& xval,
                                    DiagnosticCollector* diags) const {
    if (xcol1 < 1 || xcol1 > columnNames_.size() ||
        ycol1 < 1 || ycol1 > columnNames_.size()) {
        throw std::runtime_error("LookupTable '" + name_ +
                                 "': interpolate1D column index out of range");
    }
    size_t xc = xcol1 - 1;
    size_t yc = ycol1 - 1;
    size_t n  = data_.size();

    if (n == 0) {
        throw std::runtime_error("LookupTable '" + name_ + "': table is empty");
    }
    if (n == 1) {
        // Single row: return the single value with zero derivative
        double yv = data_[0][yc];
        return ADValue::constant(yv, xval.gradient.size());
    }

    double x = xval.value;
    size_t lo = findInterval(xc, x);
    size_t hi = lo + 1;

    double x0 = data_[lo][xc];
    double x1 = data_[hi][xc];
    double y0 = data_[lo][yc];
    double y1 = data_[hi][yc];

    double dx = x1 - x0;
    double slope = (std::abs(dx) > 1e-300) ? (y1 - y0) / dx : 0.0;

    // Clamp x to [x0, x1] for extrapolation: slope becomes 0 outside
    double t;
    double effectiveSlope;
    if (x <= x0) {
        t = 0.0;
        effectiveSlope = 0.0; // flat extrapolation: derivative = 0
    } else if (x >= x1) {
        t = 1.0;
        effectiveSlope = 0.0;
    } else {
        t = (x - x0) / dx;
        effectiveSlope = slope;
    }

    double yv = y0 + slope * (x - x0);
    // Clamp output for flat extrapolation
    if (x <= x0) yv = y0;
    if (x >= x1) yv = y1;
    (void)t; // used conceptually above

    // Propagate gradient: dy/dx_input = effectiveSlope
    size_t ng = xval.gradient.size();
    ADValue result(yv, ng);
    for (size_t i = 0; i < ng; ++i) {
        result.gradient[i] = effectiveSlope * xval.gradient[i];
    }
    return result;
}

// ============================================================================
// LookupTable::interpolate2D
// ============================================================================

ADValue LookupTable::interpolate2D(size_t xcol1, size_t ycol1, size_t zcol1,
                                    const ADValue& xval, const ADValue& yval,
                                    DiagnosticCollector* diags) const {
    if (xcol1 < 1 || xcol1 > columnNames_.size() ||
        ycol1 < 1 || ycol1 > columnNames_.size() ||
        zcol1 < 1 || zcol1 > columnNames_.size()) {
        throw std::runtime_error("LookupTable '" + name_ +
                                 "': interpolate2D column index out of range");
    }
    size_t xc = xcol1 - 1;
    size_t yc = ycol1 - 1;
    size_t zc = zcol1 - 1;
    size_t n  = data_.size();

    if (n < 4) {
        throw std::runtime_error("LookupTable '" + name_ +
                                 "': bilinear interpolation requires at least 4 rows");
    }

    // Collect unique sorted x values
    std::vector<double> xs;
    for (size_t r = 0; r < n; ++r) {
        double v = data_[r][xc];
        if (!std::isnan(v)) {
            if (xs.empty() || v > xs.back() + 1e-14) xs.push_back(v);
        }
    }
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end(),
        [](double a, double b){ return std::abs(b - a) < 1e-14; }), xs.end());

    // Collect unique sorted y values
    std::vector<double> ys;
    for (size_t r = 0; r < n; ++r) {
        double v = data_[r][yc];
        if (!std::isnan(v)) {
            if (ys.empty() || v > ys.back() + 1e-14) ys.push_back(v);
        }
    }
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end(),
        [](double a, double b){ return std::abs(b - a) < 1e-14; }), ys.end());

    size_t nx = xs.size();
    size_t ny = ys.size();

    if (nx < 2 || ny < 2) {
        if (diags) {
            diags->push(DiagnosticSeverity::Warning, "L002",
                "LookupTable '" + name_ + "': not enough distinct x or y values for bilinear interpolation",
                "lookup_table");
        }
        // Fall back to 1D on x
        return interpolate1D(xcol1, zcol1, xval, diags);
    }

    // Build z grid indexed as z[ix][iy]
    std::vector<std::vector<double>> zgrid(nx, std::vector<double>(ny,
        std::numeric_limits<double>::quiet_NaN()));

    for (size_t r = 0; r < n; ++r) {
        double xv = data_[r][xc];
        double yv = data_[r][yc];
        double zv = data_[r][zc];
        if (std::isnan(xv) || std::isnan(yv) || std::isnan(zv)) continue;

        // Find position in xs and ys
        auto itx = std::lower_bound(xs.begin(), xs.end(), xv - 1e-14);
        auto ity = std::lower_bound(ys.begin(), ys.end(), yv - 1e-14);
        if (itx == xs.end() || ity == ys.end()) continue;
        size_t ix = static_cast<size_t>(itx - xs.begin());
        size_t iy = static_cast<size_t>(ity - ys.begin());
        if (ix < nx && iy < ny) zgrid[ix][iy] = zv;
    }

    // Clamp x and y to grid extents
    double xq = xval.value;
    double yq = yval.value;
    double xClamped = std::max(xs.front(), std::min(xs.back(), xq));
    double yClamped = std::max(ys.front(), std::min(ys.back(), yq));

    // Find surrounding intervals
    auto xlo_it = std::upper_bound(xs.begin(), xs.end(), xClamped);
    if (xlo_it != xs.begin()) --xlo_it;
    size_t xlo = static_cast<size_t>(xlo_it - xs.begin());
    if (xlo + 1 >= nx) xlo = nx - 2;

    auto ylo_it = std::upper_bound(ys.begin(), ys.end(), yClamped);
    if (ylo_it != ys.begin()) --ylo_it;
    size_t ylo = static_cast<size_t>(ylo_it - ys.begin());
    if (ylo + 1 >= ny) ylo = ny - 2;

    size_t xhi = xlo + 1;
    size_t yhi = ylo + 1;

    double x0 = xs[xlo], x1 = xs[xhi];
    double y0 = ys[ylo], y1 = ys[yhi];
    double z00 = zgrid[xlo][ylo];
    double z10 = zgrid[xhi][ylo];
    double z01 = zgrid[xlo][yhi];
    double z11 = zgrid[xhi][yhi];

    double dx = x1 - x0;
    double dy = y1 - y0;

    double tx = (std::abs(dx) > 1e-300) ? (xClamped - x0) / dx : 0.0;
    double ty = (std::abs(dy) > 1e-300) ? (yClamped - y0) / dy : 0.0;

    // Bilinear interpolation
    double z = z00 * (1 - tx) * (1 - ty)
             + z10 * tx       * (1 - ty)
             + z01 * (1 - tx) * ty
             + z11 * tx       * ty;

    // Partial derivatives (only non-zero inside the grid)
    double dzdt_x = 0.0, dzdt_y = 0.0;
    bool insideX = (xq > xs.front() && xq < xs.back());
    bool insideY = (yq > ys.front() && yq < ys.back());

    if (insideX && std::abs(dx) > 1e-300) {
        // dz/dx = (1/dx) * (dz/dtx)
        double dz_dtx = (z10 - z00) * (1 - ty) + (z11 - z01) * ty;
        dzdt_x = dz_dtx / dx;
    }
    if (insideY && std::abs(dy) > 1e-300) {
        double dz_dty = (z01 - z00) * (1 - tx) + (z11 - z10) * tx;
        dzdt_y = dz_dty / dy;
    }

    size_t ng = xval.gradient.size();
    ADValue result(z, ng);
    for (size_t i = 0; i < ng; ++i) {
        result.gradient[i] = dzdt_x * xval.gradient[i] + dzdt_y * yval.gradient[i];
    }
    return result;
}

// ============================================================================
// LookupTable::lookup
// ============================================================================

ADValue LookupTable::lookup(const ADValue& rowVal, const ADValue& colVal,
                             DiagnosticCollector* /*diags*/) const {
    size_t nRows = data_.size();
    size_t nCols = columnNames_.size();

    if (nRows == 0 || nCols == 0) {
        throw std::runtime_error("LookupTable '" + name_ + "': table is empty");
    }

    double r = rowVal.value;
    double c = colVal.value;

    // Derivatives are zero when a coordinate is outside the valid range
    // (flat extrapolation / clamping matches EES behaviour).
    bool rInRange = (r >= 1.0 && r <= static_cast<double>(nRows));
    bool cInRange = (c >= 1.0 && c <= static_cast<double>(nCols));

    // Clamp to [1, N] in each direction
    double rC = std::max(1.0, std::min(static_cast<double>(nRows),  r));
    double cC = std::max(1.0, std::min(static_cast<double>(nCols), c));

    // Row interval (0-based).  For rC in [1, nRows]:
    //   rlo0 = floor(rC) - 1, clamped so that rhi0 = rlo0+1 stays in range.
    size_t rlo0 = 0, rhi0 = 0;
    double tr = 0.0;
    if (nRows > 1) {
        rlo0 = static_cast<size_t>(std::floor(rC)) - 1;
        if (rlo0 >= nRows - 1) rlo0 = nRows - 2; // keep rhi0 in bounds
        rhi0 = rlo0 + 1;
        tr = rC - static_cast<double>(rlo0 + 1); // fractional part in row space
        tr = std::max(0.0, std::min(1.0, tr));
    }

    // Column interval (0-based), same logic.
    size_t clo0 = 0, chi0 = 0;
    double tc = 0.0;
    if (nCols > 1) {
        clo0 = static_cast<size_t>(std::floor(cC)) - 1;
        if (clo0 >= nCols - 1) clo0 = nCols - 2;
        chi0 = clo0 + 1;
        tc = cC - static_cast<double>(clo0 + 1);
        tc = std::max(0.0, std::min(1.0, tc));
    }

    double v00 = data_[rlo0][clo0];
    double v10 = data_[rhi0][clo0];
    double v01 = data_[rlo0][chi0];
    double v11 = data_[rhi0][chi0];

    double result = (1.0 - tr) * (1.0 - tc) * v00
                  + tr         * (1.0 - tc) * v10
                  + (1.0 - tr) * tc          * v01
                  + tr         * tc          * v11;

    // Analytical partial derivatives w.r.t. the row and column indices.
    // Zero when clamped outside the valid range (flat extrapolation).
    double dvdr = 0.0;
    double dvdc = 0.0;
    if (rInRange && nRows > 1) {
        dvdr = (1.0 - tc) * (v10 - v00) + tc * (v11 - v01);
    }
    if (cInRange && nCols > 1) {
        dvdc = (1.0 - tr) * (v01 - v00) + tr * (v11 - v10);
    }

    size_t ng = rowVal.gradient.size();
    ADValue res(result, ng);
    for (size_t i = 0; i < ng; ++i) {
        res.gradient[i] = dvdr * rowVal.gradient[i] + dvdc * colVal.gradient[i];
    }
    return res;
}

// ============================================================================
// Aggregate functions
// ============================================================================

static std::vector<double> finiteValuesInCol(const std::vector<std::vector<double>>& data,
                                              size_t col0) {
    std::vector<double> vals;
    vals.reserve(data.size());
    for (const auto& row : data) {
        if (col0 < row.size() && !std::isnan(row[col0])) {
            vals.push_back(row[col0]);
        }
    }
    return vals;
}

double LookupTable::sumCol(size_t col1) const {
    if (col1 < 1 || col1 > columnNames_.size()) return std::numeric_limits<double>::quiet_NaN();
    auto vals = finiteValuesInCol(data_, col1 - 1);
    return std::accumulate(vals.begin(), vals.end(), 0.0);
}

double LookupTable::avgCol(size_t col1) const {
    if (col1 < 1 || col1 > columnNames_.size()) return std::numeric_limits<double>::quiet_NaN();
    auto vals = finiteValuesInCol(data_, col1 - 1);
    if (vals.empty()) return std::numeric_limits<double>::quiet_NaN();
    return std::accumulate(vals.begin(), vals.end(), 0.0) / static_cast<double>(vals.size());
}

double LookupTable::maxCol(size_t col1) const {
    if (col1 < 1 || col1 > columnNames_.size()) return std::numeric_limits<double>::quiet_NaN();
    auto vals = finiteValuesInCol(data_, col1 - 1);
    if (vals.empty()) return std::numeric_limits<double>::quiet_NaN();
    return *std::max_element(vals.begin(), vals.end());
}

double LookupTable::minCol(size_t col1) const {
    if (col1 < 1 || col1 > columnNames_.size()) return std::numeric_limits<double>::quiet_NaN();
    auto vals = finiteValuesInCol(data_, col1 - 1);
    if (vals.empty()) return std::numeric_limits<double>::quiet_NaN();
    return *std::min_element(vals.begin(), vals.end());
}

double LookupTable::stddevCol(size_t col1) const {
    if (col1 < 1 || col1 > columnNames_.size()) return std::numeric_limits<double>::quiet_NaN();
    auto vals = finiteValuesInCol(data_, col1 - 1);
    if (vals.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double variance = 0.0;
    for (double v : vals) {
        double diff = v - mean;
        variance += diff * diff;
    }
    variance /= static_cast<double>(vals.size()); // population stddev
    return std::sqrt(variance);
}

// ============================================================================
// LookupTableStore
// ============================================================================

void LookupTableStore::add(LookupTable table) {
    tables_[table.name()] = std::move(table);
}

const LookupTable* LookupTableStore::get(const std::string& name) const {
    auto it = tables_.find(name);
    return (it != tables_.end()) ? &it->second : nullptr;
}

bool LookupTableStore::has(const std::string& name) const {
    return tables_.find(name) != tables_.end();
}

void LookupTableStore::remove(const std::string& name) {
    tables_.erase(name);
}

std::vector<std::string> LookupTableStore::names() const {
    std::vector<std::string> result;
    result.reserve(tables_.size());
    for (const auto& [k, _] : tables_) result.push_back(k);
    return result;
}

// ============================================================================
// loadLookupTablesFromDirectory
// ============================================================================

LookupTableStore loadLookupTablesFromDirectory(const std::string& dir,
                                                DiagnosticCollector* diags) {
    LookupTableStore store;
    if (!fs::exists(dir) || !fs::is_directory(dir)) return store;

    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        auto path = entry.path();
        if (path.extension() != ".csv") continue;

        std::string tableName = path.stem().string();
        try {
            std::ifstream file(path);
            if (!file.is_open()) {
                if (diags) {
                    diags->push(DiagnosticSeverity::Warning, "L001",
                        "Could not open lookup table file: " + path.string(),
                        "lookup_table");
                }
                continue;
            }
            std::string content((std::istreambuf_iterator<char>(file)),
                                 std::istreambuf_iterator<char>());
            store.add(LookupTable::fromCSV(tableName, content));
        } catch (const std::exception& e) {
            if (diags) {
                diags->push(DiagnosticSeverity::Warning, "L001",
                    "Failed to parse lookup table '" + tableName + "': " + e.what(),
                    "lookup_table");
            }
        }
    }
    return store;
}

// ============================================================================
// loadLookupTableForModel
// ============================================================================

LookupTableStore loadLookupTableForModel(const std::string& modelFilePath,
                                          DiagnosticCollector* diags) {
    LookupTableStore store;
    fs::path model(modelFilePath);
    fs::path dir = model.parent_path();
    std::string stem = model.stem().string();

    // Scan the model's directory for every *.csv file whose stem starts with
    // "<stem>-".  The part after the hyphen is the table name, so that
    // "mymodel-watercp.csv" → table "watercp" and is callable as
    // LOOKUP('watercp', ...).  Files that do not match this pattern are
    // ignored, including a bare "<stem>.csv" (the table name would be empty).
    if (!fs::exists(dir) || !fs::is_directory(dir)) return store;

    const std::string prefix = stem + "-";

    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        auto p = entry.path();
        if (p.extension() != ".csv") continue;
        std::string fileStem = p.stem().string();
        // Require "<stem>-<tableName>.csv" with a non-empty table name
        if (fileStem.size() <= prefix.size() ||
            fileStem.compare(0, prefix.size(), prefix) != 0) {
            continue;
        }
        std::string tableName = fileStem.substr(prefix.size());
        try {
            std::ifstream file(p);
            if (!file.is_open()) {
                if (diags) {
                    diags->push(DiagnosticSeverity::Warning, "L001",
                        "Could not open lookup table file: " + p.string(),
                        "lookup_table");
                }
                continue;
            }
            std::string content((std::istreambuf_iterator<char>(file)),
                                 std::istreambuf_iterator<char>());
            store.add(LookupTable::fromCSV(tableName, content));
        } catch (const std::exception& e) {
            if (diags) {
                diags->push(DiagnosticSeverity::Warning, "L001",
                    "Failed to parse lookup table '" + tableName + "': " + e.what(),
                    "lookup_table");
            }
        }
    }

    return store;
}

} // namespace coolsolve
