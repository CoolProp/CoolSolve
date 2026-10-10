#pragma once

#include <string>

namespace coolsolve {

/**
 * @brief Equality of two user strings as EES defines it: without regard to case.
 *
 * EES is case-insensitive for names and for string values alike (`'AirH2O'` and
 * `'airh2o'` are the same fluid, `IF (a$ = 'ABC')` matches `'abc'`; the EES help
 * of `$IF` says so for string constants). ASCII letters and the accented letters
 * of Latin-1 written in UTF-8 (`é`/`É`, `à`/`À`, ...) are folded; every other byte
 * is compared as it is.
 *
 * Used for every comparison of user strings: `=` and `<>` in the `IF` statements
 * of FUNCTION/PROCEDURE bodies (the evaluator) and the `$IF` directives (the
 * parser). Fluid names are out of its scope.
 */
bool eesStringEqual(const std::string& a, const std::string& b);

}  // namespace coolsolve
