#pragma once

#include <cstddef>
#include <string>

namespace coolsolve {

/**
 * @file ees_string.h
 * @brief Case-insensitive comparison of names and user strings, and the non-ASCII characters
 *        accepted in names.
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

}  // namespace coolsolve
