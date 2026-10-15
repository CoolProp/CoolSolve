#pragma once

#include <cstddef>
#include <functional>
#include <string>

namespace coolsolve {

/**
 * @file ees_string.h
 * @brief Case-insensitive comparison of names and user strings, the non-ASCII characters
 *        accepted in names, and the text of numbers in user messages.
 *
 * EES does not distinguish upper and lower case in variable names or in string values
 * (`'AirH2O'` and `'airh2o'` are the same fluid, `IF (a$ = 'ABC')` matches `'abc'`, the
 * EES help of `$IF` says so for string constants). The folding done here covers the ASCII
 * letters and the accented letters of Latin-1 written in UTF-8 (`é`/`É`, `à`/`À`, ...);
 * every other byte is compared as it is. It is locale independent, never allocates in a
 * comparison and is the single definition behind `CaseInsensitiveLess` (every map of
 * variables), the string comparisons of the evaluator and of `$IF`, and the iterator of
 * `DUPLICATE`.
 */

/// The byte of `s[i]` after folding: an ASCII capital becomes lower case, and so does the
/// second byte of a capital of the Latin-1 supplement (UTF-8 `C3 80`..`C3 9E` except the
/// multiplication sign `C3 97`; the lower case is the capital + 0x20). A continuation byte
/// never equals 0xC3, so a preceding 0xC3 is always the lead byte.
inline unsigned char eesFoldByteAt(const std::string& s, std::size_t i) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c >= 'A' && c <= 'Z') return static_cast<unsigned char>(c + ('a' - 'A'));
    if (c >= 0x80 && c <= 0x9E && c != 0x97 && i > 0 && static_cast<unsigned char>(s[i - 1]) == 0xC3)
        return static_cast<unsigned char>(c + 0x20);
    return c;
}

/// Three-way comparison of the folded forms of `a` and `b` (< 0, 0, > 0): the order of
/// `CaseInsensitiveLess`, which makes a name and its accented or capitalised spelling one key.
inline int eesCompareFolded(const std::string& a, const std::string& b) {
    const std::size_t n = a.size() < b.size() ? a.size() : b.size();
    for (std::size_t i = 0; i < n; ++i) {
        // Equal bytes fold alike: the folded bytes before i are equal, and the only byte that
        // changes the folding of the next one is 0xC3, which folds to itself and to nothing else.
        if (a[i] == b[i]) continue;
        const unsigned char x = eesFoldByteAt(a, i);
        const unsigned char y = eesFoldByteAt(b, i);
        if (x != y) return x < y ? -1 : 1;
    }
    if (a.size() == b.size()) return 0;
    return a.size() < b.size() ? -1 : 1;
}

/// The lower-case form of `s` (see above); the same length as `s`.
std::string eesFoldCase(const std::string& s);

/// The upper-case form of `s`: the inverse of the folding above (ASCII letters and the Latin-1
/// letters in UTF-8: `é` becomes `É`); `ß`, `ÿ` and every other byte are kept. The same length as `s`.
std::string eesUpperCase(const std::string& s);

/**
 * @brief Equality of two user strings as EES defines it: without regard to case.
 *
 * Used for every comparison of user strings: `=` and `<>` in the `IF` statements
 * of FUNCTION/PROCEDURE bodies (the evaluator) and the `$IF` directives (the
 * parser). Fluid names are out of its scope.
 */
bool eesStringEqual(const std::string& a, const std::string& b);

/**
 * @brief Length in bytes of the non-ASCII character of a name that starts at `s[i]`, 0 when
 *        that character is not accepted in a name.
 *
 * EES (Windows, ANSI) takes any keyboard character in a name; CoolSolve reads UTF-8 and
 * accepts the characters of Latin-1 that a name can use: the letters (`é`, `ç`, `ß`, `µ`, ... not
 * the signs `×` and `÷`), the degree sign `°` and the superscript digits `²`, `³`, `¹`.
 * Other letters (Greek, Cyrillic, ...), typographic signs, spaces (U+00A0) and bytes that are
 * not UTF-8 are not accepted: they stay an error rather than a silently different variable.
 * `isLetter`, when given, tells whether the character may start a name.
 */
inline std::size_t eesNonAsciiNameCharLength(const std::string& s, std::size_t i, bool* isLetter = nullptr) {
    if (i + 1 >= s.size()) return 0;
    const unsigned char lead = static_cast<unsigned char>(s[i]);
    const unsigned char next = static_cast<unsigned char>(s[i + 1]);
    if (lead == 0xC3 && next >= 0x80 && next <= 0xBF && next != 0x97 && next != 0xB7) {   // À..ÿ without × and ÷
        if (isLetter) *isLetter = true;
        return 2;
    }
    if (lead == 0xC2 && next == 0xB5) {                                                    // µ
        if (isLetter) *isLetter = true;
        return 2;
    }
    if (lead == 0xC2 && (next == 0xB0 || next == 0xB2 || next == 0xB3 || next == 0xB9)) {  // ° ² ³ ¹
        if (isLetter) *isLetter = false;
        return 2;
    }
    return 0;
}

/**
 * @brief A number as text in an EES display format (EES help, "Display Formats").
 *
 * `format` is the two-character code of EES: `A` (auto, the second character is ignored; also the
 * result for an empty or unknown code), `F<d>` (fixed: `d` digits right of the decimal point),
 * `E<d>` (exponential: one digit left of the point, `d` right of it, and the exponent as a plain
 * integer: `5.123E2`), `G<d>` (engineering: exponent a multiple of 3, `d` significant figures),
 * `N<d>` (`d` significant figures, shown in `E` format outside 1E-6..1E6: 123.4567 is `123` in `N3`,
 * `123.5` in `N4`), `d` being a digit 0 to 9 or `M` (the most digits: 15). The auto format writes an
 * integral value without decimals and any other value with up to 10 significant digits; it is the
 * text of `STRING$(x)`. The time and date formats of EES (`T`, `D`) are not supported.
 */
std::string eesFormatNumber(double value, const std::string& format = "A");

/**
 * @brief Replace the value placeholder of a CALL WARNING / CALL ERROR message (EES help, WARNING
 *        and ERROR procedures).
 *
 * Every `XXX` of `message` (also `xxx`: EES is not case sensitive) is replaced by `render(format)`, `format` being what follows the three
 * letters: a display format code such as `A3`, `F1` or `E4` (a letter among A, N, F, E, G and a
 * digit or `M`, see `eesFormatNumber`), `$` for a string value (`XXX$`) or the empty string when
 * nothing follows (automatic format). The text after the placeholder is kept as it is.
 * `found`, when given, tells whether `message` held a placeholder.
 */
std::string eesSubstituteValue(const std::string& message,
                               const std::function<std::string(const std::string& format)>& render,
                               bool* found = nullptr);

}  // namespace coolsolve
