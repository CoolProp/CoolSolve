#pragma once
/**
 * @file integral_table.h
 * @brief IntegralTable — the time-series store for the equation-based
 *        dynamic solver, plus the `IntegralTableSpec` captured from the
 *        `$IntegralTable` directive.
 *
 * The table holds one column per tabulated variable, with the integration
 * variable always in column 0.  It supports append, direct cell access, and
 * linear interpolation (the backing store for the `INTEGRALVALUE(t,'X')`
 * built-in).  CSV serialisation mirrors the lookup-table / parametric-study
 * workflow (`<modelname>-integral.csv`).
 *
 * See `docs/integral_table_plan.md` §Phase 2.
 */
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace coolsolve {

/**
 * @brief Specification parsed from a `$IntegralTable` directive.
 *
 * Example EES syntax:  `$IntegralTable t:0.1  y  dydt  X[1..5]`
 *   - `integrationVar` = "t"
 *   - `outputInterval` = 0.1  (record a row roughly every 0.1 in t; 0 = every step)
 *   - `columns`        = ["t", "y", "dydt", "X[1]", "X[2]", ...]
 *
 * The parser (Phase 3) fills this in; `IntegralTable` and `IntegralSolver`
 * consume it.  Array-range expansion (`X[1..5]`) happens at parse time.
 */
struct IntegralTableSpec {
    std::string integrationVar;            ///< The independent variable (e.g. "t").
    double outputInterval = 0.0;            ///< Output spacing in t; 0 = every step.
    std::vector<std::string> columns;       ///< Ordered column names, integration var first.
    bool valid = false;                     ///< Set once the directive parsed cleanly.
    std::string errorMessage;               ///< Diagnostic if `valid` is false.

    /// True if the directive was present and parsed (even partially).
    bool isPresent() const { return !integrationVar.empty(); }
};

/**
 * @brief Append-only table of integration-trajectory values.
 *
 * Storage is columnar (`std::map<string, vector<double>>`) keyed by column
 * name, with the column order preserved separately so CSV output and GUI
 * rendering stay deterministic.  The integration variable is always column 0.
 */
class IntegralTable {
public:
    /// Construct an empty table whose column 0 is `integrationVar`.
    explicit IntegralTable(std::string integrationVar = "t");

    /// Define the full column set. The first entry MUST be the integration variable.
    void setColumns(const std::vector<std::string>& cols);

    /// Ordered column names (integration variable first).
    const std::vector<std::string>& columns() const { return columnOrder_; }

    /// Name of the integration variable (column 0).
    const std::string& integrationVar() const { return integrationVar_; }

    /// Append a row from a name→value map. Missing columns are left as NaN.
    /// Uses the default `std::less` comparator so brace-enclosed initializer
    /// lists (`{{"y", 1.0}, ...}`) bind directly. The caller is responsible
    /// for case-normalisation if needed.
    void appendRow(double t, const std::map<std::string, double>& values);

    /// Append a full ordered row (size must equal column count, col 0 = t).
    void appendRow(const std::vector<double>& row);

    /// Number of rows currently stored.
    int numRows() const;

    /// Access the value at (row, colIndex). Returns NaN if out of range.
    double value(int row, int colIndex) const;

    /// Read-only access to a named column. Empty if the name is unknown.
    const std::vector<double>& column(const std::string& name) const;

    /**
     * @brief Linear interpolation of column `name` at integration-variable
     *        value `t`.
     *
     * Clamps to the first/last value when `t` is outside the stored range,
     * and returns NaN if the table has no rows or the column is unknown.
     * This is the numerical core behind `INTEGRALVALUE(t,'X')`.
     */
    double interpolate(const std::string& name, double t) const;

    /// Remove every row (columns preserved).
    void clear();

    /// Render the table as CSV (header row + data rows). `sep` defaults to ','.
    std::string toCSV(char sep = ',') const;

    /**
     * @brief Write the table to `path` as CSV.
     * @return true on success.
     */
    bool writeCSV(const std::string& path) const;

private:
    std::string integrationVar_;
    std::vector<std::string> columnOrder_;
    std::map<std::string, std::vector<double>> data_;
    // Sentinel returned by column() for unknown names (mutable so the class
    // remains move-assignable).
    std::vector<double> empty_;
};

}  // namespace coolsolve
