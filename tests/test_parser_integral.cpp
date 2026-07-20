#include "coolsolve/parser.h"
#include "coolsolve/ast.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <algorithm>
#include <cctype>
#include <cmath>

using namespace coolsolve;

namespace {
// Find the first Directive statement with the given (lower-cased) name.
const Directive* findDirective(const Program& p, const std::string& lowerName) {
    for (const auto& s : p.statements) {
        if (s->is<Directive>()) {
            const auto& d = s->as<Directive>();
            std::string ln = d.name;
            std::transform(ln.begin(), ln.end(), ln.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (ln == lowerName) return &d;
        }
    }
    return nullptr;
}

bool diagnosticContains(const ParseResult& r, const std::string& needle) {
    for (const auto& d : r.diagnostics.items())
        if (d.message.find(needle) != std::string::npos) return true;
    return false;
}
}  // namespace

// ============================================================================
// Parser & AST hooks for INTEGRAL / $IntegralTable.
// ============================================================================

TEST_CASE("Parser: $IntegralTable captures spec", "[parser][integral]") {
    EESParser parser;
    auto r = parser.parse("$IntegralTable t:0.1 y dydt\n");
    REQUIRE(r.success);
    const Directive* d = findDirective(r.program, "integraltable");
    REQUIRE(d != nullptr);
    REQUIRE(d->hasIntegralTableSpec);
    REQUIRE(d->integralTableSpec.valid);
    CHECK(d->integralTableSpec.integrationVar == "t");
    CHECK(std::abs(d->integralTableSpec.outputInterval - 0.1) < 1e-9);
    CHECK(d->integralTableSpec.columns == std::vector<std::string>{"t", "y", "dydt"});
}

TEST_CASE("Parser: $IntegralTable expands X[1..5] ranges", "[parser][integral]") {
    EESParser parser;
    auto r = parser.parse("$IntegralTable t y X[1..3]\n");
    REQUIRE(r.success);
    const Directive* d = findDirective(r.program, "integraltable");
    REQUIRE(d != nullptr);
    REQUIRE(d->integralTableSpec.valid);
    const auto& cols = d->integralTableSpec.columns;
    REQUIRE(cols.size() == 5);
    CHECK(cols[0] == "t");
    CHECK(cols[1] == "y");
    CHECK(cols[2] == "X[1]");
    CHECK(cols[3] == "X[2]");
    CHECK(cols[4] == "X[3]");
}

TEST_CASE("Parser: $IntegralTable without interval defaults to 0", "[parser][integral]") {
    EESParser parser;
    auto r = parser.parse("$IntegralTable t y\n");
    REQUIRE(r.success);
    const Directive* d = findDirective(r.program, "integraltable");
    REQUIRE(d != nullptr);
    REQUIRE(d->integralTableSpec.valid);
    CHECK(d->integralTableSpec.outputInterval == 0.0);
    CHECK(d->integralTableSpec.columns.size() == 2);
}

TEST_CASE("Parser: INTEGRAL() function call parses without warning", "[parser][integral]") {
    EESParser parser;
    auto r = parser.parse("y = integral(dydt, t, 0, 4)\n");
    REQUIRE(r.success);
    // Must NOT warn that 'integral' is an unknown function.
    CHECK_FALSE(diagnosticContains(r, "Unknown function 'integral'"));
}

TEST_CASE("Parser: integralvalue() function call parses without warning", "[parser][integral]") {
    EESParser parser;
    auto r = parser.parse("z = integralvalue(t, 'y')\n");
    REQUIRE(r.success);
    CHECK_FALSE(diagnosticContains(r, "Unknown function 'integralvalue'"));
}

TEST_CASE("Parser: $IntegralAutoStep / $IntegralStop warned and ignored", "[parser][integral]") {
    EESParser parser;
    auto r = parser.parse("$IntegralAutoStep on\n$IntegralStop y>10\n");
    REQUIRE(r.success);
    CHECK(diagnosticContains(r, "not interpreted"));
    CHECK(diagnosticContains(r, "coolsolve.conf"));
    // Not flagged as an unknown directive.
    CHECK_FALSE(diagnosticContains(r, "Unknown directive"));
}

TEST_CASE("Parser: empty $IntegralTable warns and is invalid", "[parser][integral]") {
    EESParser parser;
    auto r = parser.parse("$IntegralTable\n");
    REQUIRE(r.success);
    const Directive* d = findDirective(r.program, "integraltable");
    REQUIRE(d != nullptr);
    REQUIRE(d->hasIntegralTableSpec);
    CHECK_FALSE(d->integralTableSpec.valid);
    CHECK_FALSE(d->integralTableSpec.errorMessage.empty());
}
