/**
 * @file integral_table.cpp
 * @brief Implementation of `IntegralTable` and helpers.
 *
 * Columnar storage with deterministic column order; linear interpolation by
 * binary search over the (monotonically increasing) integration-variable
 * column.  CSV output mirrors the parametric-study convention.
 */
#include "coolsolve/integral/integral_table.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

namespace coolsolve {

IntegralTable::IntegralTable(std::string integrationVar)
    : integrationVar_(std::move(integrationVar)) {
    if (!integrationVar_.empty()) {
        columnOrder_.push_back(integrationVar_);
        data_[integrationVar_];  // create empty column
    }
}

void IntegralTable::setColumns(const std::vector<std::string>& cols) {
    columnOrder_.clear();
    data_.clear();
    if (!cols.empty()) {
        integrationVar_ = cols[0];
    }
    for (const auto& c : cols) {
        columnOrder_.push_back(c);
        data_[c];  // ensure the column exists
    }
}

void IntegralTable::appendRow(double t, const std::map<std::string, double>& values) {
    // Ensure the integration-variable column exists and receives `t`.
    if (!integrationVar_.empty()) {
        data_[integrationVar_].push_back(t);
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto& colName : columnOrder_) {
        if (colName == integrationVar_) continue;
        auto it = values.find(colName);
        data_[colName].push_back(it != values.end() ? it->second : nan);
    }
}

void IntegralTable::appendRow(const std::vector<double>& row) {
    const int n = static_cast<int>(columnOrder_.size());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < n; ++i) {
        double v = (i < static_cast<int>(row.size())) ? row[i] : nan;
        data_[columnOrder_[i]].push_back(v);
    }
}

int IntegralTable::numRows() const {
    if (integrationVar_.empty()) return 0;
    auto it = data_.find(integrationVar_);
    if (it == data_.end()) return 0;
    return static_cast<int>(it->second.size());
}

double IntegralTable::value(int row, int colIndex) const {
    if (row < 0 || colIndex < 0 || colIndex >= static_cast<int>(columnOrder_.size()))
        return std::numeric_limits<double>::quiet_NaN();
    auto it = data_.find(columnOrder_[colIndex]);
    if (it == data_.end() || row >= static_cast<int>(it->second.size()))
        return std::numeric_limits<double>::quiet_NaN();
    return it->second[row];
}

const std::vector<double>& IntegralTable::column(const std::string& name) const {
    auto it = data_.find(name);
    if (it == data_.end()) return empty_;
    return it->second;
}

double IntegralTable::interpolate(const std::string& name, double t) const {
    auto cit = data_.find(name);
    if (cit == data_.end() || cit->second.empty()) return std::numeric_limits<double>::quiet_NaN();
    const auto& ys = cit->second;

    // Integration-variable column (the abscissa).
    auto tit = data_.find(integrationVar_);
    if (tit == data_.end() || tit->second.empty()) return std::numeric_limits<double>::quiet_NaN();
    const auto& xs = tit->second;
    const int n = static_cast<int>(xs.size());
    if (static_cast<int>(ys.size()) < n) return std::numeric_limits<double>::quiet_NaN();

    // Single point: return it.
    if (n == 1) return ys[0];

    // Clamp below.
    if (t <= xs.front()) return ys.front();
    // Clamp above.
    if (t >= xs.back()) return ys.back();

    // Binary search for the bracket [lo, hi] with xs[lo] <= t <= xs[hi].
    auto lower = std::lower_bound(xs.begin(), xs.end(), t);
    int hi = static_cast<int>(lower - xs.begin());
    if (hi <= 0) return ys.front();
    int lo = hi - 1;
    // Guard against a zero-width interval (duplicate timestamps).
    double dx = xs[hi] - xs[lo];
    if (dx == 0.0) return ys[lo];
    double w = (t - xs[lo]) / dx;
    return ys[lo] + w * (ys[hi] - ys[lo]);
}

void IntegralTable::clear() {
    for (auto& kv : data_) kv.second.clear();
}

std::string IntegralTable::toCSV(char sep) const {
    std::ostringstream os;
    // Header.
    for (size_t i = 0; i < columnOrder_.size(); ++i) {
        if (i) os << sep;
        os << columnOrder_[i];
    }
    os << '\n';
    // Rows.
    const int rows = numRows();
    for (int r = 0; r < rows; ++r) {
        for (size_t c = 0; c < columnOrder_.size(); ++c) {
            if (c) os << sep;
            double v = value(r, static_cast<int>(c));
            if (std::isnan(v))
                os << "";  // empty cell for missing values
            else
                os << v;
        }
        os << '\n';
    }
    return os.str();
}

bool IntegralTable::writeCSV(const std::string& path) const {
    std::ofstream f(path);
    if (!f.is_open()) return false;
    f << toCSV();
    return f.good();
}

}  // namespace coolsolve
