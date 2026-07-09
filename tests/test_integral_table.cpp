#include "coolsolve/integral/integral_table.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <fstream>
#include <cstdio>

using namespace coolsolve;
using Catch::Matchers::WithinRel;

// ============================================================================
// Phase 2 — TDD tests for IntegralTable storage + interpolation + CSV.
// See docs/integral_table_plan.md §Phase 2.
// ============================================================================

TEST_CASE("IntegralTable: columns and append-row", "[integral][integral-table]") {
    IntegralTable t("time");
    t.setColumns({"time", "y", "dydt"});
    REQUIRE(t.columns() == std::vector<std::string>{"time", "y", "dydt"});
    REQUIRE(t.integrationVar() == "time");

    t.appendRow(0.0, {{"y", 1.0}, {"dydt", -1.0}});
    t.appendRow(1.0, {{"y", 0.5}, {"dydt", -0.5}});
    REQUIRE(t.numRows() == 2);
    CHECK_THAT(t.value(0, 0), WithinRel(0.0));   // time,row0
    CHECK_THAT(t.value(1, 1), WithinRel(0.5));   // y,row1
    CHECK_THAT(t.value(0, 2), WithinRel(-1.0));  // dydt,row0
}

TEST_CASE("IntegralTable: missing column becomes NaN", "[integral][integral-table]") {
    IntegralTable t("t");
    t.setColumns({"t", "y", "z"});
    t.appendRow(0.0, {{"y", 2.0}});  // z omitted
    REQUIRE(t.numRows() == 1);
    CHECK(std::isnan(t.value(0, 2)));
    CHECK_THAT(t.value(0, 1), WithinRel(2.0));
}

TEST_CASE("IntegralTable: linear interpolation midpoints and clamping", "[integral][integral-table]") {
    IntegralTable t("t");
    t.setColumns({"t", "y"});
    // y = 2t on [0, 4]; mid-point at t=1 is y=2.
    for (int i = 0; i <= 4; ++i) {
        double tt = static_cast<double>(i);
        t.appendRow(tt, {{"y", 2.0 * tt}});
    }
    CHECK_THAT(t.interpolate("y", 1.0), WithinRel(2.0));   // exact node
    CHECK_THAT(t.interpolate("y", 1.5), WithinRel(3.0));   // midpoint
    CHECK_THAT(t.interpolate("y", 2.25), WithinRel(4.5));  // quarter
    // Clamping outside range.
    CHECK_THAT(t.interpolate("y", -5.0), WithinRel(0.0));  // below -> first
    CHECK_THAT(t.interpolate("y", 99.0), WithinRel(8.0));  // above -> last
}

TEST_CASE("IntegralTable: interpolate on empty / single-row / unknown column", "[integral][integral-table]") {
    IntegralTable t("t");
    t.setColumns({"t", "y"});
    CHECK(std::isnan(t.interpolate("y", 0.0)));  // no rows
    t.appendRow(0.0, {{"y", 7.0}});
    CHECK_THAT(t.interpolate("y", 0.0), WithinRel(7.0));  // single row
    CHECK(std::isnan(t.interpolate("nope", 0.0)));        // unknown column
}

TEST_CASE("IntegralTable: CSV round-trip", "[integral][integral-table]") {
    IntegralTable t("t");
    t.setColumns({"t", "y"});
    t.appendRow(0.0, {{"y", 0.0}});
    t.appendRow(0.5, {{"y", 0.25}});
    t.appendRow(1.0, {{"y", 1.0}});

    std::string csv = t.toCSV();
    CHECK(csv.find("t,y") != std::string::npos);
    CHECK(csv.find("0.5") != std::string::npos);

    const std::string path = "integral_table_test_out.csv";
    REQUIRE(t.writeCSV(path));
    std::ifstream f(path);
    std::string contents((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());
    CHECK(contents.find("t,y") != std::string::npos);
    std::remove(path.c_str());
}

TEST_CASE("IntegralTableSpec: presence flag", "[integral][integral-table]") {
    IntegralTableSpec spec;
    CHECK(!spec.isPresent());
    spec.integrationVar = "t";
    spec.columns = {"t", "y"};
    spec.valid = true;
    CHECK(spec.isPresent());
}
