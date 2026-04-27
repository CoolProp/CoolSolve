/**
 * Unit tests for LookupTable and LookupTableStore
 *
 * Covers: CSV round-trip, 1D interpolation (exact, midpoint, clamping),
 * AD gradient propagation, 2D bilinear interpolation with AD, aggregate
 * functions, metadata helpers, and LookupTableStore operations.
 *
 * Run with: ./coolsolve_tests "[lookup-table]"
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "coolsolve/lookup_table.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace coolsolve;

// Helper: build an ADValue scalar with no gradient variables
static ADValue constant(double v) {
    return ADValue::constant(v, 0);
}

// Helper: build an ADValue with one gradient variable (itself)
static ADValue independent(double v) {
    return ADValue::independent(v, 0, 1);
}

// Helper: build an ADValue with two gradient variables
static ADValue indep2(double v, size_t idx) {
    return ADValue::independent(v, idx, 2);
}

// Helper: constant with explicit gradient size (for mixing with independent variables)
static ADValue constant_n(double v, size_t n) {
    return ADValue::constant(v, n);
}

// ============================================================================
// CSV round-trip
// ============================================================================

TEST_CASE("LookupTable CSV round-trip", "[lookup-table]") {
    const std::string csv =
        "T,h,s\n"
        "100,2675.6,7.354\n"
        "150,2745.9,7.528\n"
        "200,2826.8,7.834\n";

    SECTION("fromCSV parses header and data") {
        auto tbl = LookupTable::fromCSV("steam", csv);
        REQUIRE(tbl.name() == "steam");
        REQUIRE(tbl.numCols() == 3);
        REQUIRE(tbl.numRows() == 3);
        REQUIRE_THAT(tbl.value(1, 1), WithinAbs(100.0, 1e-9));
        REQUIRE_THAT(tbl.value(2, 2), WithinAbs(2745.9, 1e-6));
        REQUIRE_THAT(tbl.value(3, 3), WithinAbs(7.834, 1e-6));
    }

    SECTION("toCSV round-trip preserves values") {
        auto tbl = LookupTable::fromCSV("steam", csv);
        std::string csv2 = tbl.toCSV();
        auto tbl2 = LookupTable::fromCSV("steam", csv2);
        REQUIRE(tbl2.numRows() == tbl.numRows());
        REQUIRE(tbl2.numCols() == tbl.numCols());
        for (size_t r = 1; r <= tbl.numRows(); ++r) {
            for (size_t c = 1; c <= tbl.numCols(); ++c) {
                REQUIRE_THAT(tbl2.value(r, c), WithinAbs(tbl.value(r, c), 1e-9));
            }
        }
    }

    SECTION("columnIndex is case-insensitive and 1-based") {
        auto tbl = LookupTable::fromCSV("steam", csv);
        REQUIRE(tbl.columnIndex("T") == 1);
        REQUIRE(tbl.columnIndex("t") == 1);
        REQUIRE(tbl.columnIndex("H") == 2);
        REQUIRE(tbl.columnIndex("S") == 3);
        REQUIRE(tbl.columnIndex("nonexistent") == 0);
    }

    SECTION("empty CSV throws") {
        REQUIRE_THROWS(LookupTable::fromCSV("empty", ""));
    }

    SECTION("quoted fields parsed correctly") {
        const std::string qcsv = "\"Name\",\"Value\"\n\"a,b\",1.5\n";
        auto tbl = LookupTable::fromCSV("q", qcsv);
        REQUIRE(tbl.numRows() == 1);
        REQUIRE_THAT(tbl.value(1, 2), WithinAbs(1.5, 1e-9));
    }
}

// ============================================================================
// 1D Interpolation
// ============================================================================

TEST_CASE("LookupTable interpolate1D", "[lookup-table]") {
    // Simple table: x = [0, 1, 2, 4], y = [0, 10, 20, 40]  (linear: y = 10x)
    const std::string csv = "x,y\n0,0\n1,10\n2,20\n4,40\n";
    auto tbl = LookupTable::fromCSV("lin", csv);
    const size_t xcol = 1, ycol = 2;

    SECTION("exact node values") {
        for (double xv : {0.0, 1.0, 2.0, 4.0}) {
            auto res = tbl.interpolate1D(xcol, ycol, constant(xv), nullptr);
            REQUIRE_THAT(res.value, WithinAbs(10.0 * xv, 1e-9));
        }
    }

    SECTION("midpoint interpolation") {
        auto res = tbl.interpolate1D(xcol, ycol, constant(0.5), nullptr);
        REQUIRE_THAT(res.value, WithinAbs(5.0, 1e-9));

        auto res2 = tbl.interpolate1D(xcol, ycol, constant(3.0), nullptr);
        REQUIRE_THAT(res2.value, WithinAbs(30.0, 1e-9));
    }

    SECTION("below-range clamping to first value") {
        auto res = tbl.interpolate1D(xcol, ycol, constant(-1.0), nullptr);
        REQUIRE_THAT(res.value, WithinAbs(0.0, 1e-9));
    }

    SECTION("above-range clamping to last value") {
        auto res = tbl.interpolate1D(xcol, ycol, constant(10.0), nullptr);
        REQUIRE_THAT(res.value, WithinAbs(40.0, 1e-9));
    }

    SECTION("AD gradient at interior point") {
        // x is an independent variable; y = 10x so dy/dx = 10 (slope between x=1 and x=2)
        auto xad = independent(1.5);
        auto res = tbl.interpolate1D(xcol, ycol, xad, nullptr);
        REQUIRE_THAT(res.value, WithinAbs(15.0, 1e-9));
        REQUIRE_THAT(res.gradient[0], WithinAbs(10.0, 1e-9));
    }

    SECTION("AD gradient at clamped point is zero") {
        auto xad = independent(-5.0); // below range
        auto res = tbl.interpolate1D(xcol, ycol, xad, nullptr);
        REQUIRE_THAT(res.gradient[0], WithinAbs(0.0, 1e-9));

        auto xad2 = independent(100.0); // above range
        auto res2 = tbl.interpolate1D(xcol, ycol, xad2, nullptr);
        REQUIRE_THAT(res2.gradient[0], WithinAbs(0.0, 1e-9));
    }

    SECTION("slope changes across intervals") {
        // Between x=2 and x=4: slope = (40-20)/(4-2) = 10
        // Between x=0 and x=1: slope = 10
        // Both should give slope 10 for this uniform-slope table
        auto xad3 = independent(3.0);
        auto res3 = tbl.interpolate1D(xcol, ycol, xad3, nullptr);
        REQUIRE_THAT(res3.gradient[0], WithinAbs(10.0, 1e-9));
    }
}

// ============================================================================
// 2D Bilinear Interpolation
// ============================================================================

TEST_CASE("LookupTable interpolate2D bilinear", "[lookup-table]") {
    // 2x2 grid: x in {0,1}, y in {0,1}, z = x + y
    // Row layout: (x=0,y=0,z=0), (x=0,y=1,z=1), (x=1,y=0,z=1), (x=1,y=1,z=2)
    const std::string csv =
        "x,y,z\n"
        "0,0,0\n"
        "0,1,1\n"
        "1,0,1\n"
        "1,1,2\n";
    auto tbl = LookupTable::fromCSV("grid", csv);

    SECTION("exact corners") {
        REQUIRE_THAT(tbl.interpolate2D(1, 2, 3, constant(0.0), constant(0.0), nullptr).value, WithinAbs(0.0, 1e-9));
        REQUIRE_THAT(tbl.interpolate2D(1, 2, 3, constant(1.0), constant(0.0), nullptr).value, WithinAbs(1.0, 1e-9));
        REQUIRE_THAT(tbl.interpolate2D(1, 2, 3, constant(0.0), constant(1.0), nullptr).value, WithinAbs(1.0, 1e-9));
        REQUIRE_THAT(tbl.interpolate2D(1, 2, 3, constant(1.0), constant(1.0), nullptr).value, WithinAbs(2.0, 1e-9));
    }

    SECTION("center of grid") {
        auto res = tbl.interpolate2D(1, 2, 3, constant(0.5), constant(0.5), nullptr);
        REQUIRE_THAT(res.value, WithinAbs(1.0, 1e-9)); // 0.5 + 0.5
    }

    SECTION("AD gradient for z = x + y => dz/dx = 1, dz/dy = 1") {
        auto xad = indep2(0.5, 0);
        auto yad = indep2(0.5, 1);
        auto res = tbl.interpolate2D(1, 2, 3, xad, yad, nullptr);
        REQUIRE_THAT(res.value, WithinAbs(1.0, 1e-9));
        REQUIRE_THAT(res.gradient[0], WithinAbs(1.0, 1e-6)); // dz/dx
        REQUIRE_THAT(res.gradient[1], WithinAbs(1.0, 1e-6)); // dz/dy
    }
}

// ============================================================================
// Aggregate functions
// ============================================================================

TEST_CASE("LookupTable aggregate functions", "[lookup-table]") {
    // Column 1: [1, 2, 3, 4, 5]
    const std::string csv = "v\n1\n2\n3\n4\n5\n";
    auto tbl = LookupTable::fromCSV("agg", csv);

    SECTION("sumCol") {
        REQUIRE_THAT(tbl.sumCol(1), WithinAbs(15.0, 1e-9));
    }
    SECTION("avgCol") {
        REQUIRE_THAT(tbl.avgCol(1), WithinAbs(3.0, 1e-9));
    }
    SECTION("maxCol") {
        REQUIRE_THAT(tbl.maxCol(1), WithinAbs(5.0, 1e-9));
    }
    SECTION("minCol") {
        REQUIRE_THAT(tbl.minCol(1), WithinAbs(1.0, 1e-9));
    }
    SECTION("stddevCol population stddev of [1,2,3,4,5]") {
        // Population stddev = sqrt(2) ≈ 1.4142135...
        REQUIRE_THAT(tbl.stddevCol(1), WithinAbs(std::sqrt(2.0), 1e-9));
    }
}

// ============================================================================
// LookupTableStore
// ============================================================================

TEST_CASE("LookupTableStore operations", "[lookup-table]") {
    LookupTableStore store;
    const std::string csv = "a,b\n1,2\n3,4\n";

    SECTION("initially empty") {
        REQUIRE(store.names().empty());
        REQUIRE(!store.has("any"));
    }

    SECTION("add and get") {
        store.add(LookupTable::fromCSV("tbl1", csv));
        REQUIRE(store.has("tbl1"));
        REQUIRE(store.has("TBL1")); // case-insensitive
        REQUIRE(store.get("tbl1")->numRows() == 2);
        REQUIRE(store.get("TBL1")->numRows() == 2);
    }

    SECTION("names lists all tables") {
        store.add(LookupTable::fromCSV("alpha", csv));
        store.add(LookupTable::fromCSV("beta", csv));
        auto names = store.names();
        REQUIRE(names.size() == 2);
    }

    SECTION("remove") {
        store.add(LookupTable::fromCSV("tmp", csv));
        REQUIRE(store.has("tmp"));
        store.remove("tmp");
        REQUIRE(!store.has("tmp"));
    }

    SECTION("get nonexistent returns nullptr") {
        REQUIRE(store.get("nope") == nullptr);
    }
}

// ============================================================================
// LOOKUP: EES-compatible row/col interpolation and clamping
// ============================================================================

TEST_CASE("LookupTable::lookup EES-compatible cell access", "[lookup-table]") {
    // 3-row, 2-column table:  col1 = [10, 20, 30],  col2 = [100, 200, 300]
    const std::string csv = "a,b\n10,100\n20,200\n30,300\n";
    auto tbl = LookupTable::fromCSV("lk", csv);

    SECTION("integer row and col returns exact cell value") {
        REQUIRE_THAT(tbl.lookup(constant(1), constant(1), nullptr).value, WithinAbs(10.0,  1e-9));
        REQUIRE_THAT(tbl.lookup(constant(2), constant(2), nullptr).value, WithinAbs(200.0, 1e-9));
        REQUIRE_THAT(tbl.lookup(constant(3), constant(1), nullptr).value, WithinAbs(30.0,  1e-9));
    }

    SECTION("non-integer row interpolates between rows") {
        // row=1.5, col=1: midpoint of rows 1 and 2 in col 1 = (10+20)/2 = 15
        REQUIRE_THAT(tbl.lookup(constant(1.5), constant(1), nullptr).value, WithinAbs(15.0,  1e-9));
        // row=2.5, col=2: midpoint of rows 2 and 3 in col 2 = (200+300)/2 = 250
        REQUIRE_THAT(tbl.lookup(constant(2.5), constant(2), nullptr).value, WithinAbs(250.0, 1e-9));
    }

    SECTION("non-integer col interpolates between columns") {
        // row=2, col=1.5: midpoint of cols 1 and 2 in row 2 = (20+200)/2 = 110
        REQUIRE_THAT(tbl.lookup(constant(2), constant(1.5), nullptr).value, WithinAbs(110.0, 1e-9));
    }

    SECTION("non-integer row and col does bilinear interpolation") {
        // row=1.5, col=1.5: bilinear of (10,100,20,200) corners at t=0.5,0.5
        // = 0.25*10 + 0.25*100 + 0.25*20 + 0.25*200 = (10+100+20+200)/4 = 82.5
        REQUIRE_THAT(tbl.lookup(constant(1.5), constant(1.5), nullptr).value, WithinAbs(82.5, 1e-9));
    }

    SECTION("row < 1 clamps to row 1") {
        REQUIRE_THAT(tbl.lookup(constant(0),   constant(1), nullptr).value, WithinAbs(10.0, 1e-9));
        REQUIRE_THAT(tbl.lookup(constant(-5),  constant(2), nullptr).value, WithinAbs(100.0, 1e-9));
    }

    SECTION("row > nRows clamps to last row") {
        REQUIRE_THAT(tbl.lookup(constant(4),   constant(1), nullptr).value, WithinAbs(30.0,  1e-9));
        REQUIRE_THAT(tbl.lookup(constant(100), constant(2), nullptr).value, WithinAbs(300.0, 1e-9));
    }

    SECTION("col < 1 clamps to col 1") {
        REQUIRE_THAT(tbl.lookup(constant(2), constant(0),  nullptr).value, WithinAbs(20.0, 1e-9));
        REQUIRE_THAT(tbl.lookup(constant(1), constant(-3), nullptr).value, WithinAbs(10.0, 1e-9));
    }

    SECTION("col > nCols clamps to last col") {
        REQUIRE_THAT(tbl.lookup(constant(2), constant(3),  nullptr).value, WithinAbs(200.0, 1e-9));
        REQUIRE_THAT(tbl.lookup(constant(3), constant(99), nullptr).value, WithinAbs(300.0, 1e-9));
    }

    SECTION("AD derivative w.r.t. row (non-integer, inside range)") {
        // At row=1.5, col=1: slope w.r.t. row = v[row2][col1] - v[row1][col1] = 20 - 10 = 10
        auto rowAD = independent(1.5);
        auto colC  = constant_n(1.0, 1); // col fixed, same gradient size
        auto res = tbl.lookup(rowAD, colC, nullptr);
        REQUIRE_THAT(res.value,       WithinAbs(15.0, 1e-9));
        REQUIRE_THAT(res.gradient[0], WithinAbs(10.0, 1e-9));
    }

    SECTION("AD derivative w.r.t. col (non-integer, inside range)") {
        // At row=2, col=1.5: slope w.r.t. col = v[row2][col2] - v[row2][col1] = 200 - 20 = 180
        auto rowC  = constant_n(2.0, 1); // row fixed, same gradient size
        auto colAD = independent(1.5);
        auto res = tbl.lookup(rowC, colAD, nullptr);
        REQUIRE_THAT(res.value,       WithinAbs(110.0, 1e-9));
        REQUIRE_THAT(res.gradient[0], WithinAbs(180.0, 1e-9));
    }

    SECTION("AD derivative is zero when clamped") {
        auto rowAD = independent(0.0);  // below range
        auto res = tbl.lookup(rowAD, constant_n(1.0, 1), nullptr);
        REQUIRE_THAT(res.gradient[0], WithinAbs(0.0, 1e-9));

        auto rowAD2 = independent(10.0);  // above range
        auto res2 = tbl.lookup(rowAD2, constant_n(1.0, 1), nullptr);
        REQUIRE_THAT(res2.gradient[0], WithinAbs(0.0, 1e-9));

        auto colAD = independent(0.0);  // below range
        auto res3 = tbl.lookup(constant_n(2.0, 1), colAD, nullptr);
        REQUIRE_THAT(res3.gradient[0], WithinAbs(0.0, 1e-9));
    }

    SECTION("single-row table returns that row's value") {
        auto tbl1 = LookupTable::fromCSV("one", "v\n42\n");
        REQUIRE_THAT(tbl1.lookup(constant(1),  constant(1), nullptr).value, WithinAbs(42.0, 1e-9));
        REQUIRE_THAT(tbl1.lookup(constant(0),  constant(1), nullptr).value, WithinAbs(42.0, 1e-9));
        REQUIRE_THAT(tbl1.lookup(constant(99), constant(1), nullptr).value, WithinAbs(42.0, 1e-9));
    }
}

// ============================================================================
// loadLookupTableForModel: <modelStem>-<tableName>.csv naming convention
// ============================================================================

TEST_CASE("loadLookupTableForModel uses <modelStem>-<tableName>.csv convention",
          "[lookup-table][lookup-naming]") {
    namespace fs = std::filesystem;
    auto tmpDir = fs::temp_directory_path() / "coolsolve_lookup_naming_test";
    fs::remove_all(tmpDir);
    fs::create_directories(tmpDir);

    // Helper: write a tiny CSV file with a single header + value
    auto writeCsv = [&](const std::string& filename, const std::string& content) {
        std::ofstream f(tmpDir / filename);
        f << content;
    };

    // The model file itself need not exist on disk for the loader, but its
    // path is parsed for the parent directory and stem.  We still create it
    // so the test mirrors a real model layout.
    const std::string modelStem = "mymodel";
    fs::path modelPath = tmpDir / (modelStem + ".eescode");
    { std::ofstream f(modelPath); f << "x = 1\n"; }

    SECTION("file <stem>-<tableName>.csv loads as table 'tableName'") {
        writeCsv("mymodel-data.csv",    "v\n1\n2\n3\n");
        writeCsv("mymodel-watercp.csv", "v\n10\n20\n");
        auto store = loadLookupTableForModel(modelPath.string());
        REQUIRE(store.size() == 2);
        REQUIRE(store.has("data"));
        REQUIRE(store.has("watercp"));
        REQUIRE(store.get("data")->numRows()    == 3);
        REQUIRE(store.get("watercp")->numRows() == 2);
        // The old "full file stem" naming MUST NOT be used
        REQUIRE(!store.has("mymodel-data"));
        REQUIRE(!store.has("mymodel_data"));
    }

    SECTION("bare <stem>.csv is ignored (no table name after the hyphen)") {
        writeCsv("mymodel.csv", "v\n1\n");
        auto store = loadLookupTableForModel(modelPath.string());
        REQUIRE(store.empty());
    }

    SECTION("legacy <stem>_<suffix>.csv naming is no longer accepted") {
        writeCsv("mymodel_watercp.csv", "v\n1\n2\n");
        auto store = loadLookupTableForModel(modelPath.string());
        REQUIRE(store.empty());
    }

    SECTION("CSVs not starting with <stem>- are ignored") {
        writeCsv("other-data.csv",  "v\n1\n");
        writeCsv("variables.csv",    "v\n1\n");
        writeCsv("mymodel-keep.csv", "v\n1\n2\n");
        auto store = loadLookupTableForModel(modelPath.string());
        REQUIRE(store.size() == 1);
        REQUIRE(store.has("keep"));
    }

    SECTION("table name preserves additional hyphens after the first one") {
        writeCsv("mymodel-air-cp.csv", "v\n1\n");
        auto store = loadLookupTableForModel(modelPath.string());
        REQUIRE(store.has("air-cp"));
    }

    SECTION("malformed CSV emits a diagnostic and is skipped") {
        // Empty file → fromCSV throws; the loader must not throw out
        writeCsv("mymodel-bad.csv", "");
        DiagnosticCollector diags;
        auto store = loadLookupTableForModel(modelPath.string(), &diags);
        REQUIRE(!store.has("bad"));
        // At least one warning about the failed parse
        bool sawWarning = false;
        for (const auto& d : diags.items()) {
            if (d.code == "L001") { sawWarning = true; break; }
        }
        REQUIRE(sawWarning);
    }

    fs::remove_all(tmpDir);
}
