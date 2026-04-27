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
