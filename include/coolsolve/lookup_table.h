#pragma once

#include "coolsolve/autodiff_node.h"
#include "coolsolve/diagnostic.h"
#include "coolsolve/ir.h"
#include <map>
#include <string>
#include <vector>

namespace coolsolve {

// ============================================================================
// LookupTable
// ============================================================================

/**
 * @brief An in-memory lookup table backed by a CSV file.
 *
 * A table is a named rectangular matrix of doubles with a header row.
 * Row and column indices are 1-based throughout the public API, matching
 * EES conventions.  Empty CSV cells are stored as NaN.
 *
 * Linear interpolation is provided for 1-D and 2-D (bilinear) cases with
 * full forward-mode AD derivative support.  Extrapolation clamps to the
 * nearest boundary value (flat extrapolation), consistent with EES.
 */
class LookupTable {
public:
    LookupTable() = default;

    /** @brief Construct from explicit column names and row data. */
    LookupTable(std::string name,
                std::vector<std::string> columnNames,
                std::vector<std::vector<double>> data);

    // ------------------------------------------------------------------
    // Factories
    // ------------------------------------------------------------------

    /**
     * @brief Parse a CSV string into a LookupTable.
     *
     * The first non-empty line is treated as the header row (column names).
     * Subsequent lines are numeric data rows.  Empty cells become NaN.
     * Whitespace around fields is stripped.
     *
     * @param name  Table name (used for error messages and as the key in
     *              LookupTableStore).
     * @param csv   Raw CSV content.
     * @return The parsed table.  Throws std::runtime_error on malformed
     *         input.
     */
    static LookupTable fromCSV(const std::string& name, const std::string& csv);

    /** @brief Serialise the table back to CSV (header row + data rows). */
    std::string toCSV() const;

    // ------------------------------------------------------------------
    // Metadata
    // ------------------------------------------------------------------

    const std::string& name()    const { return name_; }
    size_t             numRows() const { return data_.size(); }
    size_t             numCols() const { return columnNames_.size(); }

    const std::vector<std::string>&          columnNames() const { return columnNames_; }
    const std::vector<std::vector<double>>&  data()        const { return data_; }

    /**
     * @brief Return the 1-based column index for a column name.
     *
     * Matching is case-insensitive.  Returns 0 if not found.
     */
    size_t columnIndex(const std::string& name) const;

    /**
     * @brief Return the value at 1-based (row, col).
     *
     * Throws std::out_of_range if indices are out of bounds.
     */
    double value(size_t row1, size_t col1) const;

    // ------------------------------------------------------------------
    // Interpolation
    // ------------------------------------------------------------------

    /**
     * @brief 1-D linear interpolation (EES INTERPOLATE / INTERPOLATE1).
     *
     * Given x from column xcol1 and y from column ycol1 (both 1-based),
     * compute y at xval by linear interpolation.  The x column is assumed
     * to be monotonically increasing.  Values outside the range are clamped
     * (flat extrapolation).
     *
     * The AD gradient is propagated through the piecewise-linear slope.
     *
     * @param xcol1  1-based index of the independent variable column.
     * @param ycol1  1-based index of the dependent variable column.
     * @param xval   The query point (with AD gradient).
     * @param diags  Optional diagnostic collector for warnings.
     */
    ADValue interpolate1D(size_t xcol1, size_t ycol1,
                          const ADValue& xval,
                          DiagnosticCollector* diags = nullptr) const;

    /**
     * @brief 2-D bilinear interpolation (EES INTERPOLATE2 / INTERPOLATE2DM).
     *
     * The table must contain a structured grid: xcol1 lists unique x values
     * (repeated for each y level) and ycol1 lists unique y values (in a
     * regular pattern).  The result column zcol1 provides z at each (x,y)
     * grid node.
     *
     * Partial derivatives w.r.t. both x and y are computed analytically and
     * propagated through the bilinear formula.
     *
     * @param xcol1  1-based index of the first independent variable column.
     * @param ycol1  1-based index of the second independent variable column.
     * @param zcol1  1-based index of the dependent variable column.
     * @param xval   Query value for x (with AD gradient).
     * @param yval   Query value for y (with AD gradient).
     * @param diags  Optional diagnostic collector for warnings/errors.
     */
    ADValue interpolate2D(size_t xcol1, size_t ycol1, size_t zcol1,
                          const ADValue& xval, const ADValue& yval,
                          DiagnosticCollector* diags = nullptr) const;

    /**
     * @brief EES-compatible LOOKUP: retrieve a value at a (possibly
     * non-integer) row and column index, with linear interpolation and
     * boundary clamping.
     *
     * Behaviour matches EES LOOKUP:
     * - Row and column indices are 1-based.
     * - Non-integer indices trigger linear interpolation between adjacent
     *   rows/columns (bilinear when both are non-integer).
     * - Indices below 1 clamp to the first row/column.
     * - Indices above the last row/column clamp to the last row/column.
     * - Analytical AD derivatives are propagated; the derivative is zero
     *   when the index is outside [1, N] (flat extrapolation).
     *
     * @param rowVal  1-based row index (with AD gradient).
     * @param colVal  1-based column index (with AD gradient).
     * @param diags   Optional diagnostic collector for warnings.
     */
    ADValue lookup(const ADValue& rowVal, const ADValue& colVal,
                   DiagnosticCollector* diags = nullptr) const;

    // ------------------------------------------------------------------
    // Aggregate functions (EES SUMLOOKUP etc.)
    // ------------------------------------------------------------------

    /** @brief Sum of all finite values in col1 (1-based). */
    double sumCol(size_t col1) const;

    /** @brief Mean of all finite values in col1 (1-based).  NaN if empty. */
    double avgCol(size_t col1) const;

    /** @brief Maximum of all finite values in col1 (1-based).  NaN if empty. */
    double maxCol(size_t col1) const;

    /** @brief Minimum of all finite values in col1 (1-based).  NaN if empty. */
    double minCol(size_t col1) const;

    /** @brief Population standard deviation of finite values in col1.  NaN if < 2. */
    double stddevCol(size_t col1) const;

private:
    std::string                      name_;
    std::vector<std::string>         columnNames_;
    std::vector<std::vector<double>> data_;   // data_[row][col]

    // Internal helper: find the interval [lo, hi) in a sorted column such
    // that data_[lo][col0] <= x < data_[hi][col0].  Returns lo == hi == 0
    // when the column has fewer than 2 rows.
    //
    // col0 is 0-based internally.
    size_t findInterval(size_t col0, double x) const;
};

// ============================================================================
// LookupTableStore
// ============================================================================

/**
 * @brief Registry of LookupTable objects for a single solve session.
 *
 * Tables are keyed by their name (case-insensitive).  The store is passed
 * to ExpressionEvaluator so that lookup/interpolate calls can resolve tables
 * at evaluation time.
 */
class LookupTableStore {
public:
    /** @brief Add or replace a table (key is table.name(), case-insensitive). */
    void add(LookupTable table);

    /** @brief Return a pointer to the named table, or nullptr if not found. */
    const LookupTable* get(const std::string& name) const;

    /** @brief True if a table with that name exists. */
    bool has(const std::string& name) const;

    /** @brief Remove a table by name.  No-op if not present. */
    void remove(const std::string& name);

    /** @brief Return all table names in sorted order. */
    std::vector<std::string> names() const;

    /** @brief Return the number of tables. */
    size_t size() const { return tables_.size(); }

    /** @brief True if the store contains no tables. */
    bool empty() const { return tables_.empty(); }

private:
    std::map<std::string, LookupTable, CaseInsensitiveLess> tables_;
};

/**
 * @brief Load all *.csv files in a directory into a LookupTableStore.
 *
 * Each file name (without extension) becomes the table name.  Files that
 * fail to parse emit a warning to @p diags (if provided) and are skipped.
 */
LookupTableStore loadLookupTablesFromDirectory(const std::string& dir,
                                                DiagnosticCollector* diags = nullptr);

/**
 * @brief Load the companion lookup tables for a model file.
 *
 * Scans the model's directory for CSV files following the convention
 * **`<modelStem>-<tableName>.csv`**.  The part after the first hyphen is
 * used as the table name, so `mymodel-watercp.csv` becomes the table
 * `watercp` and is referenced from equations as `LOOKUP('watercp', ...)`,
 * `INTERPOLATE('watercp', ...)`, etc.
 *
 * Files that do not match this pattern are ignored, including a bare
 * `<modelStem>.csv` (since it would imply an empty table name).  Files
 * that fail to parse emit a warning to @p diags (if provided).
 */
LookupTableStore loadLookupTableForModel(const std::string& modelFilePath,
                                          DiagnosticCollector* diags = nullptr);

} // namespace coolsolve
