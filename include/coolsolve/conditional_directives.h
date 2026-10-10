#pragma once

#include "coolsolve/ir.h"      // CaseInsensitiveLess
#include "coolsolve/units.h"

#include <map>
#include <set>
#include <string>

namespace coolsolve {

/**
 * @file conditional_directives.h
 * @brief Conditions of the EES compile-time directives `$IF` / `$IFNOT` (CS-GAP-IF-DIRECTIVE).
 *
 * EES compiles the Equations window top to bottom; the lines between `$IF cond`
 * (or `$IFNOT cond`), an optional `$ELSE` and `$ENDIF` are compiled only if the
 * condition holds (not for `$IFNOT`), the others are ignored. The parser
 * (parser.cpp) keeps the stack of open directives and skips the lines of the
 * removed branch; this file decides whether a condition is true.
 *
 * Conditions (EES help, "$IF and $IFNOT"): `StringVariable$='text'` (quotes
 * optional, case-insensitive, the variable must be set to a string constant
 * *before* the directive), the flags below, `TABLERUN# = < > n`,
 * `UNITSYSTEM('XX')`, and any combination with `AND`, `OR` and parentheses
 * (`AND` binds tighter than `OR`). Values in a CoolSolve solve, which is neither
 * a parametric table run nor an optimisation:
 *
 * | Condition                                              | Value |
 * |--------------------------------------------------------|-------|
 * | `PARAMETRICTABLE`, `PARAMETRICTABLE='name'`            | false (the GUI parametric study replaces the imposed value in the text of each point and solves an ordinary model) |
 * | `MINMAX` / `MIN/MAX`, `MINMAXTABLE` / `MIN/MAXTABLE` [`='name'`] | false |
 * | `UNCERTAINTY`, `UNCERTAINTYTABLE` [`='name'`]          | false |
 * | `TABLERUN# = n`, `< n`, `> n`                          | false (no table run) |
 * | `DIAGRAMWINDOW`, `MACRO`                               | false (no Diagram window, no macro) |
 * | `PROFESSIONAL`                                         | true (no license restriction: a feature CoolSolve lacks fails loudly) |
 * | `INTEGRALTABLE`                                        | true when a (kept) `$IntegralTable` directive is in the file, above or below (EES compiles in two passes for it) |
 * | `UNITSYSTEM('XX')` [`= n`]                             | as the `UNITSYSTEM()` function (SI, Mass, Deg; C, Pa, J) |
 * | `VERSION < n`, `> n`                                   | error (CoolSolve has no EES version number) |
 *
 * Anything else is an error, as is a test on a string variable that has no
 * constant value at that point.
 */

/** What the conditions can query: the compile-time knowledge of the parser at the directive. */
struct CompileTimeState {
    /// String variables set to a string constant (`X$ = 'text'`, `X$ = Y$` of a constant) before the directive.
    std::map<std::string, std::string, CaseInsensitiveLess> strings;
    /// String variables set by something else (a CALL output, an expression): value unknown at compile time.
    std::set<std::string, CaseInsensitiveLess> computedStrings;
    /// A kept `$IntegralTable` directive exists in the file (`INTEGRALTABLE` condition).
    bool integralTable = false;
    UnitSystem units;
};

/** Outcome of `evaluateDirectiveCondition`. */
struct ConditionResult {
    bool ok = false;                 ///< False: `error` says why the condition cannot be evaluated.
    bool value = false;              ///< Value of the condition (not negated for `$IFNOT`).
    bool usedIntegralTable = false;  ///< The `INTEGRALTABLE` test was evaluated (it needs `CompileTimeState::integralTable`).
    std::string error;
};

/**
 * @brief Evaluate the text that follows `$IF` / `$IFNOT` (comments already removed).
 *
 * Every test of the condition is parsed and checked, even when an earlier one
 * already decides the result: an unknown or unusable test is never ignored.
 * A trailing `then` is accepted (seen in valid EES files).
 */
ConditionResult evaluateDirectiveCondition(const std::string& condition, const CompileTimeState& state);

}  // namespace coolsolve
