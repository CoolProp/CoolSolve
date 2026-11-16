// The substitute of the floating-point std::to_chars / std::from_chars that coolplot-core is compiled with on macOS
// (cmake/coolplot_charconv_compat.h, option COOLSOLVE_COOLPLOT_CHARCONV_COMPAT): the same text and the same values
// as the standard functions, which this test compares where they are available (not on macOS before 26.0).
#include <catch2/catch_test_macros.hpp>
#include "../cmake/coolplot_charconv_compat.h"

#include <clocale>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

std::string compatText(double value, std::chars_format format, int precision) {
    char buffer[400];
    auto result = coolsolve_charconv::to_chars(buffer, buffer + sizeof(buffer), value, format, precision);
    REQUIRE(result.ec == std::errc());
    return std::string(buffer, result.ptr);
}

std::vector<double> sampleValues() {
    std::vector<double> values = {0.0, -0.0, 1.0, -1.0, 0.1, 1.0 / 3.0, 2.5, 1e-3, 9.999e-4, 1e6, 123456.789,
                                  -273.15, 101325.0, 1e-300, 1e300, 5e-324, std::numeric_limits<double>::max(),
                                  std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> mantissa(-10.0, 10.0);
    std::uniform_int_distribution<int> exponent(-30, 30);
    for (int i = 0; i < 2000; ++i) values.push_back(mantissa(rng) * std::pow(10.0, exponent(rng)));
    return values;
}

}  // namespace

TEST_CASE("charconv substitute writes numbers in the C locale", "[charconv]") {
    CHECK(compatText(1.5, std::chars_format::fixed, 3) == "1.500");
    CHECK(compatText(1234.5, std::chars_format::scientific, 2) == "1.23e+03");
    CHECK(compatText(0.000123, std::chars_format::general, 12) == "0.000123");
    CHECK(compatText(1e20, std::chars_format::general, 6) == "1e+20");

    // A locale with a decimal comma does not change the text
    const char* german = std::setlocale(LC_NUMERIC, "de_DE.UTF-8");
    if (german) {
        CHECK(compatText(1.5, std::chars_format::fixed, 1) == "1.5");
        double value = 0.0;
        const std::string text = "2.25";
        auto result = coolsolve_charconv::from_chars(text.data(), text.data() + text.size(), value);
        CHECK(result.ec == std::errc());
        CHECK(value == 2.25);
        std::setlocale(LC_NUMERIC, "C");
    }

    char small[16];
    volatile size_t room = 4;  // not a constant: the compiler would warn of the truncation it is meant to test
    auto tooLong = coolsolve_charconv::to_chars(small, small + room, 123.456, std::chars_format::fixed, 3);
    CHECK(tooLong.ec == std::errc::value_too_large);
}

TEST_CASE("charconv substitute reads numbers as std::from_chars", "[charconv]") {
    auto parse = [](const std::string& text, double& value) {
        return coolsolve_charconv::from_chars(text.data(), text.data() + text.size(), value);
    };
    double value = 7.0;
    auto r = parse("-1.25e3", value);
    CHECK((r.ec == std::errc() && value == -1250.0));
    const std::string partial = "1e";  // the number is "1", the "e" is left
    r = coolsolve_charconv::from_chars(partial.data(), partial.data() + partial.size(), value);
    CHECK((r.ec == std::errc() && value == 1.0 && r.ptr == partial.data() + 1));
    value = 7.0;
    CHECK(parse(" 1", value).ec == std::errc::invalid_argument);
    CHECK(parse("+1", value).ec == std::errc::invalid_argument);
    CHECK(parse("abc", value).ec == std::errc::invalid_argument);
    CHECK(parse("", value).ec == std::errc::invalid_argument);
    CHECK(value == 7.0);
    CHECK(parse("1e999", value).ec == std::errc::result_out_of_range);
    CHECK(parse("-1e-999", value).ec == std::errc::result_out_of_range);
    r = parse("0x1A", value);
    CHECK((r.ec == std::errc() && value == 0.0));
    CHECK((parse("inf", value).ec == std::errc() && std::isinf(value)));
}

#if !defined(__APPLE__)
TEST_CASE("charconv substitute gives the text and values of the standard functions", "[charconv]") {
    const std::vector<double> values = sampleValues();
    const std::chars_format formats[] = {std::chars_format::fixed, std::chars_format::scientific,
                                         std::chars_format::general};
    for (double v : values) {
        for (std::chars_format format : formats) {
            for (int precision : {0, 1, 3, 6, 12, 17}) {
                if (format == std::chars_format::fixed && std::abs(v) > 1e200) continue;  // 300-digit texts
                char buffer[400];
                auto std_result = std::to_chars(buffer, buffer + sizeof(buffer), v, format, precision);
                REQUIRE(std_result.ec == std::errc());
                const std::string expected(buffer, std_result.ptr);
                const std::string actual = compatText(v, format, precision);
                if (actual != expected) FAIL_CHECK(expected << " != " << actual);

                double a = 0.0, b = 0.0;
                auto ra = coolsolve_charconv::from_chars(expected.data(), expected.data() + expected.size(), a);
                auto rb = std::from_chars(expected.data(), expected.data() + expected.size(), b);
                REQUIRE(ra.ec == rb.ec);
                REQUIRE(ra.ptr == rb.ptr);
                if (ra.ec == std::errc()) {
                    if (!(a == b || (std::isnan(a) && std::isnan(b)))) FAIL_CHECK(expected << ": " << a << " != " << b);
                }
            }
        }
    }
}
#endif
