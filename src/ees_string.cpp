#include "coolsolve/ees_string.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>

namespace coolsolve {

std::string eesFoldCase(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) out += static_cast<char>(eesFoldByteAt(s, i));
    return out;
}

bool eesStringEqual(const std::string& a, const std::string& b) {
    return a.size() == b.size() && eesCompareFolded(a, b) == 0;
}

namespace {

/// The digits of a display-format code: 0 to 9, or 15 for 'M' (the most the display holds).
bool formatDigits(char c, int& digits) {
    if (c >= '0' && c <= '9') { digits = c - '0'; return true; }
    if (c == 'M') { digits = 15; return true; }
    return false;
}

/// `value` rounded to `significant` figures, as "<sign><digits>" and the decimal exponent
/// (the rounding can raise it: 9.99 to 2 figures is 10, "10", exponent 1).
void significantDigits(double value, int significant, std::string& digits, int& exponent, bool& negative) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*e", significant - 1, std::fabs(value));
    digits.clear();
    const char* p = buf;
    for (; *p && *p != 'e'; ++p)
        if (*p != '.') digits += *p;
    exponent = *p == 'e' ? std::atoi(p + 1) : 0;
    negative = std::signbit(value);
}

/// `E<d>`: one digit left of the point, `d` right of it, the exponent without sign padding.
std::string exponentialFormat(double value, int decimals) {
    std::string digits;
    int exponent = 0;
    bool negative = false;
    significantDigits(value, decimals + 1, digits, exponent, negative);
    std::string out = negative ? "-" : "";
    out += digits.substr(0, 1);
    if (decimals > 0) out += "." + digits.substr(1);
    return out + "E" + std::to_string(exponent);
}

/// `G<s>`: `s` significant figures, the exponent a multiple of 3.
std::string engineeringFormat(double value, int significant) {
    if (value == 0.0) return "0";
    std::string digits;
    int exponent = 0;
    bool negative = false;
    significantDigits(value, significant, digits, exponent, negative);
    const int power = static_cast<int>(std::floor(exponent / 3.0)) * 3;
    const std::size_t before = static_cast<std::size_t>(exponent - power) + 1;   // 1 to 3 digits left of the point
    if (digits.size() < before) digits.append(before - digits.size(), '0');
    std::string out = negative ? "-" : "";
    out += digits.substr(0, before);
    if (digits.size() > before) out += "." + digits.substr(before);
    return out + "E" + std::to_string(power);
}

/// `N<s>`: `s` significant figures, in E format outside 1E-6..1E6.
std::string significantFormat(double value, int significant) {
    if (significant < 1) significant = 1;
    const double magnitude = std::fabs(value);
    if (magnitude != 0.0 && (magnitude >= 1e6 || magnitude < 1e-6)) return exponentialFormat(value, significant - 1);
    std::string digits;
    int exponent = 0;
    bool negative = false;
    significantDigits(value, significant, digits, exponent, negative);
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", std::max(0, significant - 1 - exponent), value);
    return buf;
}

std::string autoFormat(double value) {
    if (std::isfinite(value) && value == std::trunc(value) && std::fabs(value) < 1e16)
        return std::to_string(static_cast<long long>(value));
    std::ostringstream oss;
    oss.precision(10);
    oss << value;
    return oss.str();
}

bool isFormatLetter(char c) {
    return c == 'A' || c == 'N' || c == 'F' || c == 'E' || c == 'G' ||
           c == 'a' || c == 'n' || c == 'f' || c == 'e' || c == 'g';
}

}  // namespace

std::string eesFormatNumber(double value, const std::string& format) {
    int digits = 0;
    if (!std::isfinite(value) || format.size() < 2 || !formatDigits(format[1], digits)) return autoFormat(value);
    char buf[512];
    switch (std::toupper(static_cast<unsigned char>(format[0]))) {
        case 'F':
            if (format[1] == 'M') std::snprintf(buf, sizeof buf, "%.15g", value);
            else std::snprintf(buf, sizeof buf, "%.*f", digits, value);
            return buf;
        case 'E': return exponentialFormat(value, digits);
        case 'G': return engineeringFormat(value, std::max(digits, 1));
        case 'N': return significantFormat(value, digits);
        default:  return autoFormat(value);
    }
}

std::string eesSubstituteValue(const std::string& message,
                               const std::function<std::string(const std::string& format)>& render,
                               bool* found) {
    std::string out;
    bool any = false;
    for (std::size_t i = 0; i < message.size();) {
        // EES is not case sensitive: the library files write xxxE4
        if (i + 3 > message.size() || !std::equal(message.begin() + i, message.begin() + i + 3, "xxx",
                [](char a, char b) { return std::tolower(static_cast<unsigned char>(a)) == b; })) {
            out += message[i++];
            continue;
        }
        any = true;
        std::size_t next = i + 3;
        std::string format;
        int digits = 0;
        if (next < message.size() && message[next] == '$') {
            format = "$";
            next += 1;
        } else if (next + 1 < message.size() && isFormatLetter(message[next]) && formatDigits(message[next + 1], digits) &&
                   // 'M' (most digits) only ends a word: XXXFMAX is the value followed by "FMAX"
                   !(message[next + 1] == 'M' && next + 2 < message.size() &&
                     std::isalnum(static_cast<unsigned char>(message[next + 2])))) {
            format = message.substr(next, 2);
            next += 2;
        }
        out += render(format);
        i = next;
    }
    if (found) *found = any;
    return out;
}

}  // namespace coolsolve
