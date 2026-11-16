// coolplot_charconv_compat.h -- the floating-point std::to_chars and std::from_chars of coolplot-core, for macOS.
//
// coolplot-core (src/core.cpp, src/compose.cpp) writes numbers with std::to_chars(first, last, value, format,
// precision) and reads them with std::from_chars(first, last, value). The libc++ of macOS declares these
// floating-point overloads unavailable before macOS 13.3 (to_chars) and 26.0 (from_chars), so a build for older
// systems does not compile. On Apple, CMakeLists.txt compiles copies of the two files whose calls go to the
// functions below instead (option COOLSOLVE_COOLPLOT_CHARCONV_COMPAT).
//
// The text is the same: to_chars with a format and a precision is defined as printf with %.*f, %.*e or %.*g in the
// "C" locale, which snprintf gives here with the "C" locale made current for the thread; from_chars reads what
// strtod reads in the "C" locale, less the leading spaces, the plus sign and the hexadecimal prefix it refuses.

#pragma once

#include <cctype>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <locale.h>
#include <string>
#include <system_error>
#if defined(__APPLE__)
#include <xlocale.h>
#endif

namespace coolsolve_charconv {

/// The "C" locale, current for this thread for the lifetime of the object (the global locale is untouched)
class CLocaleScope {
public:
    CLocaleScope() : previous_(uselocale(c_locale())) {}
    ~CLocaleScope() { uselocale(previous_); }
    CLocaleScope(const CLocaleScope&) = delete;
    CLocaleScope& operator=(const CLocaleScope&) = delete;

private:
    static locale_t c_locale() {
        static const locale_t c = newlocale(LC_ALL_MASK, "C", static_cast<locale_t>(0));
        return c;
    }
    locale_t previous_;
};

inline std::to_chars_result to_chars(char* first, char* last, double value, std::chars_format format,
                                     int precision) {
    const char* spec = format == std::chars_format::fixed        ? "%.*f"
                       : format == std::chars_format::scientific ? "%.*e"
                                                                 : "%.*g";
    if (format == std::chars_format::hex) return {last, std::errc::invalid_argument};  // not used by coolplot-core
    int n = 0;
    {
        CLocaleScope c;
        n = std::snprintf(first, static_cast<size_t>(last - first), spec, precision, value);
    }
    if (n < 0 || n >= last - first) return {last, std::errc::value_too_large};
    return {first + n, std::errc()};
}

inline std::from_chars_result from_chars(const char* first, const char* last, double& value) {
    const char* p = first;
    if (p != last && *p == '-') ++p;
    // No leading space or plus sign: the first character after the sign starts the number
    if (p == last || !(std::isdigit(static_cast<unsigned char>(*p)) || *p == '.' || *p == 'i' || *p == 'I' ||
                       *p == 'n' || *p == 'N')) {
        return {first, std::errc::invalid_argument};
    }
    // "0x1A": std::from_chars reads the 0 and stops at the x
    if (last - p > 1 && p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) {
        value = p != first ? -0.0 : 0.0;
        return {p + 1, std::errc()};
    }
    const std::string text(first, last);
    char* end = nullptr;
    double v = 0.0;
    int error = 0;
    {
        CLocaleScope c;
        errno = 0;
        v = std::strtod(text.c_str(), &end);
        error = errno;
    }
    if (end == text.c_str()) return {first, std::errc::invalid_argument};
    const char* ptr = first + (end - text.c_str());
    if (error == ERANGE && (v == 0.0 || std::isinf(v))) return {ptr, std::errc::result_out_of_range};
    value = v;
    return {ptr, std::errc()};
}

}  // namespace coolsolve_charconv
