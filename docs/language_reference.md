# CoolSolve Language Reference

This page describes the CoolSolve language — an EES-compatible subset with
thermodynamic property functions backed by CoolProp. Only features that are
currently implemented and covered by tests are documented here.

> **⚠ Unit system (fixed).** In its current state, **the unit system of
> CoolSolve cannot be changed.** Temperature is always in **°C**, pressure in
> **Pa**, energy in **J**, etc. (see `include/coolsolve/units.h`). CoolProp
> works internally in SI (K), so values are converted `°C → K` on the way in
> and `K → °C` on the way out — this is transparent to the model. **Never use
> Kelvin values for `T=...` arguments** or you will be ~273 °C off.

---

## 1. Syntax basics

### Equations and assignments

Every non-empty line that is not a directive, comment, or loop construct is
an **equation**:

```ees
x = 3 * y + 2
h_1 = enthalpy('Water', T=T_1, P=P_1)
```

Equations are not executed top-to-bottom — they are solved simultaneously.

A `;` separates several statements on one line (equations, assignments and
`CALL`s), as in the body of a FUNCTION or PROCEDURE. A comment may follow a
separator, and a `;` inside a string (`'a;b'`) or a comment does not separate:

```ees
a = 1 ; b = 2 ; "two constants"
z = 2 ; CALL p(1 : y) ; // an equation and a CALL
```

`IF ... THEN`, `REPEAT`, `GOTO` and statement labels are not statements of the
main program (see section 6): they are errors there, after a `;` too.

### Comments

Three comment styles are supported:

| Syntax            | Kind                                                 |
|:------------------|:-----------------------------------------------------|
| `"..."` (dquote)  | Block or inline comment; may span several lines      |
| `{...}`           | Brace comment (may span several lines, may nest)     |
| `// ...`          | Single-line comment                                  |

A `"` comment that is not closed on its line — after an equation or alone on
a line — continues up to the next `"`, as in EES:

```ees
Q_dot = m_dot*(h_1 - h_2)   "heat flow rate of the evaporator,
                             including the superheat"
```

The text of the comment is never parsed. A comment that is never closed hides
the rest of the file: the parser warns (`P007`) with the line that opened it.

### Numbers and identifiers

- Numbers: `1`, `3.14`, `1.5e-3`, `-6E10`.
- Identifiers (variable names) start with a letter or `_` and may contain letters, digits,
  `_` and the symbols `%`, `|` and `#`, anywhere in the name: `C%`, `Q_dot_2_%`, `K|star_r_n`,
  `Nu#_f`. A `\` may join two parts of a name (`dU\dtau`, `time_h\day`) and must be followed
  by a letter or `_`. The built-in constants `pi`, `g#`, `R#`... keep their meaning, in the
  main program and in the bodies of `FUNCTION`s and `PROCEDURE`s alike (a body does not see the
  variables of the main program, but it does see the constants). A variable of a body that has the
  name of a constant (a local `c#`) is a local one: it keeps the value the body gives it.
- Letters are the ASCII letters and the accented letters of Latin-1 (`é`, `ç`, `ß`, `µ`...),
  and the characters `°`, `²`, `³`, `¹` may follow them: `T_a_réf`, `n°_90°elbow_reg`. The
  model file must be UTF-8. Other letters (Greek, Cyrillic...), signs such as `×` or `−`,
  spaces that are not ASCII and bytes that are not UTF-8 (a file saved in Latin-1) are errors
  that name the character: *Invalid variable name 'Δp': the character 'Δ' (U+0394) is not
  allowed in a name*.
- Names are **case-insensitive**, as in EES, accented letters included: `T_a_réf` and
  `T_A_RÉF` are one variable (the spelling of the first occurrence is kept in the results).
- Identifiers ending with `$` are **string variables** (e.g. `fluid$`, `fuel%$`); `$` is
  only allowed at the end of a name. Array elements are written `X%[i]`, `T_réf[3]`; a string
  variable has no elements (string arrays are not supported, see *String arrays* below).
- EES takes any keyboard character in a name except `( ) * / + - ^ { } ;` (and, in the older
  manual, `' | : "`); CoolSolve also refuses the characters that are part of its own syntax
  (`'`, `"`, `.`, `,`, `:`, `=`, `<`, `>`, `[`, `]`) and `@ & ~ ? !`. A name that starts with
  a digit or with `%`, `|`, `#` is an error. A keyword followed by one of these characters
  (`call%`, `if|x`) is a name, not the keyword.
- In the LaTeX report the `%`, `#`, `&`, `\` and `°` of a name are escaped, an accented
  letter is typeset with `\text{}` and the `|` starts a superscript as in the EES formatted
  equations (`K|star_r_n` is typeset K<sup>star</sup><sub>r,n</sub>).
- `.initials` and `.sol` files list the names as written (`Nu#_f = 6.0`); in them a `#` starts a
  comment only at the start of a line or after a blank.

### Operators and precedence

From tightest to loosest:

1. `^` (power, right-associative)
2. Unary `+`, `-`
3. `*`, `/`
4. `+`, `-`
5. `<`, `<=`, `>`, `>=` (return 1/-1 in condition contexts). As in EES, `=>` is also read as `>=` and `=<` as `<=`
   (written without a blank between the two characters: `= >` and `==` are errors)

### Units annotations

An optional **unit annotation** `"[unit]"` may follow an equation:

```ees
T_1 = 25  "[°C]"
```

The annotation is purely informational — it is **not** used to convert the
value. The unit system is fixed (temperature in °C, pressure in Pa, energy in
J, …); see the warning at the top of this page.

A **numeric literal**, a **parenthesised group**, a **function call** or an
**array element** inside an expression, a function call or the value of a named
argument may carry its own annotation in brackets: `if(tau, 3600[s], 10000[W],
10000[W], 0[W])`, `T_A = 25.0[K] + 273.15 [K]`, `T = -5[C]`, `p = 1E5 [Pa]`,
`(25 + 273.15) [K]`, `h = enthalpy(Water, P=1E5[Pa], T=20)`, `Cmin = min(C_h,
C_c) [W/K]`, `y = x[2] [K] + max(a, b) [kW]`. The annotation ends the term it
follows, at the end of a right-hand side as well as inside an expression. The
number (or the term) is used as written and the unit is dropped: **nothing is
converted** (`P=1e2[kPa]` is 100 Pa, not 100 kPa). One annotation per term;
after a call or an array element it must hold a unit, not a number (`f(x)[2]`
and `A[1][2]` are not annotations, see below).

A **plain variable** does not carry an annotation: `y = a [kW]` has the form of
the array element `a[kW]`, so the unit would be read as an index. When the
index is a single identifier that names a unit (`K`, `C`, `Pa`, `kPa`, `bar`,
`J`, `kJ`, `kg`, `m`, `s`, `W`, `kW`, ... spelled as in EES: `[kw]` is not
recognised, nor is a lower-case `j` or `k`, usual names of indices) and is not
set in the scope, it is a parse error with the line: the main program is the scope of
its equations, a `FUNCTION` / `PROCEDURE` the scope of its inputs and of the
variables its body assigns. A loop index (`DUPLICATE`, the index of an indexed
`SUM`) or a variable of that name that the scope sets stays a valid index
(`DUPLICATE s = 1, 3` ⏎ `x[s] = 2*s`; `K = 2` ⏎ `y = x[K]`). Compound units
(`a [kg/s]`) are not covered. Put the unit in the comment of the equation
(`y = a "[kW]"`) or after a number:

```ees
y = a [kW]           "'a [kW]': a unit annotation after a variable is not supported ('kW' would be read as an array index, and nothing sets 'kW' here); ..."
y = a "[kW]"         "accepted: the unit in the comment of the equation"
y = 3600 [s] * a     "accepted: the unit after a number"
```

**Nothing else may follow** the closing parenthesis of a call or the closing
bracket of an array element (nor the closing quote of a string): the text used
to be ignored (`y = max(1,2) 3` gave 2, `y = x[1] 2` gave `x[1]`). It is a
parse error that names the text and the call:

```ees
y = max(1,2) 3       "Unexpected text '3' after the call 'max(1,2)'"
y = x[1] 2           "Unexpected text '2' after the array element 'x[1]'"
y = A[1][2]          "Unexpected text '[2]' after the array element 'A[1]'"  (write A[1,2])
CALL p(1 : r) 3      "Unexpected text '3' after the arguments of 'p'"
```

### Function arguments and indices

An argument or an index that cannot be parsed is a parse error that names the
function (or array), the position and the text; it is never skipped:

```ees
y = max(3 4, 2)      "Could not parse argument 1 of 'max': '3 4'"
y = x[3 4]           "Could not parse index 1 of 'x': '3 4'"
y = max(3, )         "Could not parse argument 2 of 'max': ''"
```

The same holds for the arguments of `CALL name(inputs : outputs)`. An output
is a variable or, in the main program, an expression (`CALL p(a : b/1000)`, see
*Outputs of a CALL* in §7). A call without argument is written `f()`. An array
range `x[1..3]` stands for the arguments `x[1], x[2], x[3]` (see *Array ranges*
in §8); it is not an expression of its own.

`name=value` arguments (`T=`, `P=`, `x=`, …) belong to the thermophysical
property functions (§3). A built-in function that takes positional arguments
only (`IF`, `max`, `sin`, …) refuses them with a parse error — in particular a
comparison inside `IF(a=5, 1, 2)` is not valid (see *Inline selection* below).

### Directives

Parser-level directives begin with `$` and span a single line:

- `$if <cond>` / `$ifnot <cond>` / `$else` / `$endif` — conditional compilation,
  see below.
- `$UnitSystem SI Mass Deg Pa C J` — recognised for EES compatibility but **ignored**; the unit
  system is fixed (see the warning at the top of this page). A directive that asks for
  other units gets a **warning**, see below.
- `$warnings ON` / `OFF` — recognised and **ignored**; the warnings of `CALL WARNING` are
  always reported, as in EES (see *CALL ERROR and CALL WARNING* in §7).
- `$CONSTANT Name# = value [units]` — a constant of the model, see below.

#### `$UnitSystem`

CoolSolve's units are fixed: **SI, mass basis, degrees, Pa, °C, J**, whatever the
directive says. A model written for another setting (kPa, bar, MPa, kJ, K, molar, radians,
English units) would therefore be read in the fixed units, with wrong results and no other
sign, so the parser warns (`P010`, at the line of the directive) when the settings of a
`$UnitSystem` directive differ from them, and names what CoolSolve uses instead:

```ees
$UnitSystem SI Mass Deg kPa C kJ
P = 100                          "read as 100 Pa"
```

```
warning (line 1): $UnitSystem SI Mass Deg kPa C kJ: CoolSolve's units are fixed (SI, mass basis,
degrees, Pa, °C, J) and the directive does not change them: this model asks for kPa (pressures
are read in Pa), kJ (energies are read in J). Convert the model to SI-C-Pa-J (docs/ees_import.md §6)
```

- The settings are those of the `UNITSYSTEM()` function, in any case: `SI`/`Eng`, `Mass`/`Molar`,
  `Deg`/`Rad`, `C`/`K`/`F`/`R`, `Pa`/`kPa`/`bar`/`MPa`/`psia`/`atm`, `J`/`kJ`/`Btu`/`kcal`
  (separated by blanks or commas). `SI Mass Deg Pa C J`, or any part of it, gives no warning.
- A word that is none of these settings is named in the same warning
  (`'furlong' is not a unit setting CoolSolve knows`); it is not an error.
- There is one warning per `$UnitSystem` line that differs. A directive in a comment or in a
  removed `$IF` branch is not read.
- The warning does not change the solve: the model is solved in the fixed units. `UNITSYSTEM('kPa')`
  and `$IF UNITSYSTEM('kPa')` stay false, `UNITSYSTEM('Pa')` stays true, whatever the directive says.

To use such a model, convert it by hand
([EES import §6](ees_import.md#6-unit-system--manual-conversion)).

#### `$CONSTANT`

```ees
$CONSTANT imax# = 20                "number of segments"
$CONSTANT L_seg# = 2/imax# [m]      "an expression of numbers and constants defined above"
```

As in EES, a constant is a value that cannot change and that **every part of the
model sees**: the main program and the bodies of all `FUNCTION`s and `PROCEDURE`s,
whether the directive is above or below them. Its name ends with `#`. The value
is a number or an expression **without spaces** of numbers and of the constants
defined above (`+ - * / ^`, parentheses, the built-in constants such as `pi`); the
`[units]` that may follow, and the English-unit value EES allows after them, are
ignored. Every use of the constant is replaced by its value when the model is
compiled, so a constant is not a variable of the solution; it is the usual way to
size arrays (`DUPLICATE i = 1, imax#`, `T[1..imax#]`, §6 and §8). A routine that
declares the name in its header (`PROCEDURE p(imax#, …)`) keeps it as its own
variable. Errors, with the line: assigning a constant (*"'imax#' is a constant
($CONSTANT, line 3): its value cannot be changed"*), a value that uses a variable
or a later constant, a name without `#`, a second definition of the name, and a
string constant (`$CONSTANT R#$ = 'R134a'`), which is not supported.

#### Conditional compilation: `$IF`, `$IFNOT`, `$ELSE`, `$ENDIF`

As in EES, the lines between `$IF <cond>` and `$ELSE` / `$ENDIF` are compiled
only if the condition is true, the lines between `$ELSE` and `$ENDIF` only if it
is false; `$IFNOT` inverts the condition. The lines of the other branch are
**removed before anything else looks at them** — they are not parsed (they may
hold text CoolSolve does not support), define no variable and no equation, and
appear in no report. Everything may be in a branch: equations, comments, `CALL`s,
`DUPLICATE` blocks, `FUNCTION`/`PROCEDURE` definitions, other directives.
Directives may be nested, and may also be used inside the body of a
`FUNCTION`/`PROCEDURE`. The keywords are case-insensitive.

```ees
act$ = 'a'
$if act$='a'
x = 5
$else
x = 7                      "removed: x = 5"
$endif
```

The condition is either a test or several tests combined with `AND`, `OR` and
parentheses (`AND` binds tighter than `OR`):

| Condition | Value in a CoolSolve solve |
|:----------|:---------------------------|
| `X$='text'` | true if the string variable `X$` equals `text` **without regard to case** (`X$` and `'text'` as in EES: the quotes are optional). `X$` must be set to a **string constant** (`X$ = 'text'`, or `X$ = Y$` with `Y$` a constant) **before** the directive, in the main program, outside a removed branch; the value of a `CALL` output or of an expression is not known when the equations are compiled and is an error naming the variable, as is a variable that has no value (EES also accepts a value from its Diagram window, which CoolSolve does not have). `<>` is not an EES condition: use `$IFNOT`. |
| `PARAMETRICTABLE`, `PARAMETRICTABLE='name'` | **false**: a CoolSolve solve is not an EES table run, and the GUI parametric study replaces the imposed value in the text of each point and solves an ordinary model. So `$IFNOT PARAMETRICTABLE` blocks are always kept |
| `MINMAX` (`MIN/MAX`), `MINMAXTABLE` (`MIN/MAXTABLE`) [`='name'`], `UNCERTAINTY`, `UNCERTAINTYTABLE` [`='name'`], `TABLERUN# = n` / `< n` / `> n`, `DIAGRAMWINDOW`, `MACRO` | false |
| `PROFESSIONAL` | true (no license restriction: what CoolSolve does not implement fails loudly) |
| `INTEGRALTABLE` | true if the file contains a (kept) `$IntegralTable` directive, above or below the test (EES compiles in two passes for it) |
| `UNITSYSTEM('XX')` [`= n`] | as the `UNITSYSTEM()` function: `SI`, `Mass`, `Deg`, `C`, `Pa`, `J` are in use, `Eng`, `Molar`, `Rad`, `K`, `kPa`, … are not |

`VERSION < n` / `> n` and any other condition are errors: CoolSolve has no EES
version number. A `$IF` takes only `=` (and `<` / `>` after `TABLERUN#`): `>=`, `<=` and
their EES spellings `=>` and `=<` are errors there. Tolerated EES spellings: a trailing `then` after the condition
(`$if type$='ari' then`), any text after `$else` / `$endif` (`$else fuel$='ch4'`),
a quoted or `{ }` comment after the condition, spaces around `=`, and a bare
`$ELSEIF` (it is a `$ELSE` in the EES help; `$ELSEIF <cond>` is an error).
Unbalanced directives are errors with their line: `$ELSE` or `$ENDIF` without
`$IF`, a second `$ELSE`, a `$IF` without `$ENDIF`.

With `-d`, the `report.md` of the debug folder lists every `$IF`/`$IFNOT` with
the branch that was kept. See also [`$IF` in the EES help](https://fchartsoftware.com/ees/eeshelp/_if_and__ifnot.htm).

---

## 2. Built-in mathematical functions

### Elementary functions

| Function                 | Meaning                               |
|:-------------------------|:--------------------------------------|
| `sin`, `cos`, `tan`      | Trigonometric, argument in **degrees** (`sin(30) = 0.5`) |
| `sinh`, `cosh`, `tanh`   | Hyperbolic                            |
| `arcsin` / `asin`        | Inverse sine, result in degrees       |
| `arccos` / `acos`        | Inverse cosine, result in degrees     |
| `arctan` / `atan`        | Inverse tangent, result in degrees    |
| `arcsinh` / `asinh`      | Inverse hyperbolic sine               |
| `arccosh` / `acosh`      | Inverse hyperbolic cosine             |
| `arctanh` / `atanh`      | Inverse hyperbolic tangent            |
| `atan2(y, x)`            | Two-argument arctangent, result in degrees |
| `exp(x)`, `ln(x)`        | Natural exponential / logarithm       |
| `log(x)` / `log10(x)`    | Base-10 logarithm                     |
| `sqrt(x)`                | Square root                           |
| `abs(x)`                 | Absolute value                        |
| `pi` / `pi()`            | Mathematical constant π               |

The trigonometric functions use **degrees**, like the default EES setting
`$UnitSystem … DEG`; the hyperbolic functions have no angle. A model written
for the EES `RAD` setting must convert its angles (`sin(x*180/pi)`,
`arctan(y)*pi/180`), see [EES import §6](ees_import.md#6-unit-system--manual-conversion).

### Special functions

| Function                 | Meaning                                            |
|:-------------------------|:---------------------------------------------------|
| `erf(x)`                 | Gauss error function                               |
| `erfc(x)`                | Complementary error function                       |
| `lmtd_f(DT1, DT2)` / `lmtd` | Log-mean temperature difference `(ΔT₁−ΔT₂)/ln(ΔT₁/ΔT₂)`, smooth at `ΔT₁≈ΔT₂` |

### Rounding & arithmetic utilities

| Function                 | Meaning                         |
|:-------------------------|:--------------------------------|
| `ceil`, `floor`          | Nearest integer up / down       |
| `round`                  | Round half away from zero       |
| `trunc`                  | Truncate toward zero            |
| `sign(x)`                | −1, 0, or +1                    |
| `mod(x, y)`              | `x` modulo `y`                  |

### Date and time functions

The time functions of the EES library `SETP.LIB` (Duffie & Beckman, *Solar
Engineering of Thermal Processes*) are built in. Names are case-insensitive.

| Function                 | Meaning                         |
|:-------------------------|:--------------------------------|
| `nDay_(month, day)`      | Day number `n` of the year, e.g. `nDay_(12, 4) = 338` |
| `NumDay_(month)`         | Number of days of the month: 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 |
| `AveDay_(month)`         | Day number `n` of the "average day" of the month (Table 1.6.1): 17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344 |
| `EqnTime_(n)`            | Equation of time in minutes (eq. 1.5.3) on the day number `n` |

All of them count in a **non-leap year** (365 days, the `n` of Duffie &
Beckman): `nDay_(1, 1) = 1`, `nDay_(2, 28) = 59`, `nDay_(3, 1) = 60`,
`nDay_(12, 31) = 365`, and `NumDay_` summed over the 12 months is 365. The
average days of `AveDay_` are the days whose extraterrestrial radiation is about
the monthly mean (17 Jan, 16 Feb, 16 Mar, 15 Apr, 15 May, 11 Jun, 17 Jul,
16 Aug, 15 Sep, 15 Oct, 14 Nov, 10 Dec); they do not hold above the Arctic
circle or below the Antarctic circle.

**Month.** As in `SETP.LIB`, `month` is rounded to the nearest integer (half
away from zero) before use: `nDay_(2.5, 1)` is the first of month 3 and
`12.0001` is month 12. Only a month that rounds outside 1 to 12 is an error.

**Day.** `nDay_` uses `day` as given, so it may have a fraction
(`nDay_(12, 4.5) = 338.5`). It must be at least 1 and at most the length of its
month in a non-leap year plus `0.001`: `nDay_(2, 29)` and `nDay_(4, 31)` are
errors, `nDay_(2, 28.0005)` is accepted.

**Equation of time.** `EqnTime_(n)` is `229.2 (0.000075 + 0.001868 cos B -
0.032077 sin B - 0.014615 cos 2B - 0.04089 sin 2B)` with
`B = 360 (n - 1) / 365` degrees, e.g. `EqnTime_(338) = 9.593` (4 December) and
a minimum of about -14.3 minutes in mid-February. Like the library function it
has no range check: any `n` is accepted and the result is periodic with period
365.

**Errors.** An invalid argument stops the calculation with a message such as
`NumDay_: month must be between 1 and 12 after rounding, got 13`,
`nDay_: day 29 exceeds the 28 days of month 2` or
`nDay_: day must be at least 1, got 0`, as does a call with the wrong number of
arguments (`nDay_` takes 2, the others 1). An iterate of the solver that leaves
the valid range is an evaluation error.

**Derivatives.** `nDay_` has derivative 1 with respect to `day` and 0 with
respect to `month`; `NumDay_` and `AveDay_` are piecewise constant (derivative
0); `EqnTime_` is smooth and its derivative with respect to `n` is exact (forward-mode
automatic differentiation). All four can therefore be used with an unknown
argument: `nDay_(12, d) = 340` solves for `d = 6`, and `EqnTime_(n) = -10`
solves for a day number `n` (the equation of time is not monotonic, so give a
guess value in the `.initials` file).

### Aggregation functions

| Function                 | Meaning                         |
|:-------------------------|:--------------------------------|
| `sum(a, b, …)`           | Sum of arguments                |
| `sum2d(a, b, …)`         | Sum of arguments (same as `sum`) |
| `average(a, b, …)`       | Arithmetic mean                 |
| `product(a, b, …)`       | Product of arguments            |
| `stddev(a, b, …)`        | Population standard deviation   |

`min(a, b, …)` and `max(a, b, …)` take any number of arguments too. An array range
is a list of arguments: `average(T[1..imax#])`, `sum(x[1..n], 10, y)` (§8).

**Indexed form** (EES). `sum(expr, i=a, b)` is the sum of `expr` for the index `i`
from `a` to `b`, and `product(expr, i=a, b)` the product:
`mm = sum(x[i]*MM[i], i=1, 4)`, `h = sum(x[i]*MM[i]*h[i], i=1, 4)/mm`. When the file
is read, CoolSolve replaces it by the plain `sum` (`product`) of the copies of `expr`
with `i` replaced by each value, as a `DUPLICATE` loop does (the equation, and its
LaTeX form, show the terms).

- The index is local to the aggregation: it is never a variable of the model, and a
  model variable (or the index of an enclosing `DUPLICATE`) of the same name is not
  affected. The expression need not use it (`sum(2*y, k=1, 5)` is `10*y`).
- The bounds are those of a `DUPLICATE` (§6): numbers, `$CONSTANT`s, variables set to
  a constant expression (`Nprod = 4`, then `sum(…, i=1, Nprod)`), arithmetic of those,
  or an expression of the index of an enclosing `DUPLICATE` (`u[k] = sum(x[i], i=1, k)`,
  expanded in each copy). A bound that is not known before the solve, or not an
  integer, is an error with the line, as for a loop. An empty range (the first bound
  greater than the last) is an error too, not a sum of no term.
- An index with arithmetic is evaluated: `sum(x[i+1] - x[i], i=1, 4)` uses `x[2] - x[1]`,
  …; sums nest (`sum(sum(a[i,j], j=1, i), i=1, n)`).
- In a `FUNCTION` or `PROCEDURE` the bounds are numbers or `$CONSTANT`s (the body is
  expanded when the file is read); a bound that is an input of the routine is an error.
- The form is `sum(expression, i=first, last)` and nothing else: `sum(x[i], i=1)`,
  `sum(x[i], 4, i=1)` are errors. Only `sum` and `product` have it: `average`, `stddev`
  and `sum2d` refuse an `i=…` argument with an error naming the line (*"AVERAGE: the
  indexed form AVERAGE(expression, i=first, last) is implemented for SUM and PRODUCT
  only"*); for them use a range (`average(x[1..n])`).

### Inline selection

Two forms, selected by the number of arguments (case-insensitive, usable in the
main program and inside `FUNCTION` / `PROCEDURE` bodies, nestable, also in
`INTEGRAL` models):

| Form | Result |
|:-----|:-------|
| `if(cond, a, b)` | `a` if `cond > 0`, otherwise `b` |
| `IF(A, B, X, Y, Z)` | the EES 5-argument form: `X` if `A < B`, `Y` if `A = B`, `Z` if `A > B` |

```ees
Q = if(tau, 3600, 10000, 10000, 0)      "10000 W until tau = 3600 s, then 0"
```

- The test `A = B` is **exact** floating-point equality, as in EES (no
  tolerance): `0.1 + 0.2` is greater than `0.3`.
- The value and the derivatives are those of the selected branch (`X`, `Y` or
  `Z`, or `a` / `b`); the tested arguments (`cond`, or `A` and `B`) are
  piecewise-constant selectors and contribute no derivative. A branch may
  therefore hold unknowns of a coupled block (`Q = if(Q_dry, Q_wet, Q_wet, Q_dry, Q_dry)`).
- All arguments are evaluated, as for any function: the branch that is not
  selected must be computable too.
- Any other number of arguments is an error that names the two accepted forms.
  The string form `IF$(...)` is not implemented.
- There is no comparison inside the function: `IF(a=5, 1, 2)` is the parse
  error *"IF does not take named arguments ('a=5'); …"* (the `=` would be read
  as a named argument). Test an equality with the 5-argument form,
  `IF(a, 5, 2, 1, 2)` (2 if `a < 5`, 1 if `a = 5`, 2 if `a > 5`), and a sign
  with `IF(a - 5, X, Y)` (`X` if `a - 5 > 0`).

### String / number conversions

| Function                 | Meaning                              |
|:-------------------------|:-------------------------------------|
| `STRING$(x)`             | Format a numeric value as a string   |
| `STRINGVAL(s$)`          | Parse a numeric string to a number   |
| `UPPERCASE$(s$)`         | The string in upper case             |
| `LOWERCASE$(s$)`         | The string in lower case             |

`UPPERCASE$` and `LOWERCASE$` take one string (a `'literal'`, a string variable or another
string function, nesting included); a number is an error (*Uppercase$(): the argument must be a
string*, `STRING$(x)` converts a number), and so is any other number of arguments. The ASCII
letters and the accented letters of Latin-1 are converted (`é` becomes `É`); every other
character (digits, signs, `ß`, ...) is kept. The result is a string, used where a string is
expected: as the value of a string variable in a `FUNCTION` or `PROCEDURE` body, also to
reassign an input (`Pr$ = Uppercase$(Pr$)` changes the copy of the procedure, not the string
of the caller), as an operand of `=` and `<>` in an `IF`, as the fluid of a property function
(`enthalpy(Uppercase$(f$), T=T, P=P)`), as the string argument of a `CALL` or of a `FUNCTION`,
and as the message of `CALL WARNING` / `CALL ERROR`. String comparisons already ignore the case
(`IF (Pr$ = 'FREEZE')` matches `'freeze'`), so the conversion matters when the text itself is
kept or shown: a string output, a message. In the main program `s$ = Uppercase$(f$)` gives the
string variable `s$` the converted text (see *String variables of the main program* below);
used as a number (`x = Uppercase$('a') + 1`) it is an error (*Uppercase$ returns a string and
cannot be used as a number*).

### String variables of the main program

A string variable of the main program (`g$`) gets its string from one of these equations, in any
order of the lines (the structural analysis orders them, so a chain `h$ = g$`, `g$ = f$`,
`f$ = 'Water'` works whichever line comes first):

| Equation | The string of `g$` |
|:---------|:-------------------|
| `g$ = 'Water'` | the literal |
| `g$ = f$` | the string of the string variable `f$`, whatever gives `f$` its value (a literal, another alias, a function, a `CALL` output) |
| `g$ = Uppercase$(f$)`, `Lowercase$(...)`, `STRING$(x)`, `PHASE$(fluid, T=..., P=...)`, `LOOKUP$('table', row, col)` | the text of the function; its arguments (`f$`, `x`, the state, `row`, `col`) are computed first, as for any equation |
| `CALL proc(... : g$)` | the string output of the procedure |

The string is set before the equations that read it are solved: `h = enthalpy(g$, T=100, P=1e5)`
uses it as the fluid wherever the line is, and `CALL p(g$ : y)` passes it to a procedure whose body
compares it (`IF (s$ = 'Water')`) or hands it on. Such an equation is never solved as a number; the
solution check compares the strings (case aside). The strings of the solution are in the `.sol`
file (`g$ = 'Water'`), the variable tables and the JSON result.

Errors, which name the variables and are reported before the solve:

- strings defined in terms of each other, directly or through a function (`f$ = g$` with
  `g$ = f$`, `f$ = f$`): *String variables defined from each other: f$ -> g$ -> f$ (a string
  needs a value: give one of them a string constant, e.g. f$ = 'text')*;
- a string variable that nothing sets (`g$ = f$` with no equation for `f$`, or `enthalpy(fluid$, ...)`
  with no `fluid$`): the model is not square and the message ends with *String variable(s) with no
  value: f$ (a string gets its value from X$ = 'text', ...)*;
- a string defined from a variable that holds a number (`f$ = 3`, `g$ = f$`): *String variable not
  found: f$*.

A string variable assigned a number or an arithmetic expression (`g$ = 3`) is not a string: it is
solved as a number, as before (it is not in the string variables of the result).
`$IF` directives still need a string constant (see §1, *Conditional compilation*), whatever the
alias.

### String arrays are not supported

EES accepts arrays of strings (`Prod$[1]='CO2'`, `UO$=U$[Pro]`). CoolSolve does not, by
decision: an element of a string variable (`X$[...]`, in an assignment or in a read, in the
main program, in a `FUNCTION` or `PROCEDURE` body, in a `DUPLICATE`, in a `CALL`) is an
error that names the line and says how to rewrite the model, instead of a string that is
silently lost (stored as the number 0 in the main program, unreadable in a body):

```
Line 3: String arrays are not supported: 'Prod$[i]' is an element of the string array Prod$.
Put the strings in a lookup table (one string per row, file <model>-<table>.csv), keep the row
number in the array instead (a numeric array) and read the string with LOOKUP$('table', row, col) ...
```

The rewrite: the strings go into a lookup table with one string per row (§11.8), the model
keeps the row number, and `LOOKUP$('table', row, col)` gives the string wherever a string is
expected, in particular as the fluid of a property call:

```ees
{EES}
Prod$[1] = 'CO2' ; Prod$[2] = 'H2O'
DUPLICATE i=1,2
  H[i] = MolarMass(Prod$[i])*enthalpy(Prod$[i], T=T)
END

{CoolSolve: the file model-Species.csv has the header  Name  and the rows  CO2  and  H2O}
DUPLICATE i=1,2
  H[i] = MolarMass(LOOKUP$('Species', i, 1))*enthalpy(LOOKUP$('Species', i, 1), T=T)
END
```

### Unit conversions

| Function                           | Meaning                                     |
|:-----------------------------------|:--------------------------------------------|
| `CONVERT('from', 'to')`            | Conversion factor so `1 [from] = factor [to]` |
| `CONVERTTEMP('from', 'to', T)`     | Convert a temperature value between the scales `C`, `K`, `F`, `R` |
| `UNITSYSTEM('setting')`            | 1 if the unit setting is in use, 0 otherwise (`'SI'`, `'Eng'`, `'Mass'`, `'Molar'`, `'Deg'`, `'Rad'`, `'C'`, `'K'`, `'F'`, `'R'`, `'Pa'`, `'kPa'`, `'bar'`, `'MPa'`, `'psia'`, `'atm'`, `'J'`, `'kJ'`, `'Btu'`, `'kcal'`) — CoolSolve always answers 1 for `SI`, `Mass`, `Deg`, `C`, `Pa` and `J` |

`UNITSYSTEM` lets EES library routines check the unit settings they were
written for; in CoolSolve such a guard never fires (the unit system is fixed).
It reports the fixed units whatever a `$UnitSystem` directive of the file says (that
directive only gets a warning, see *Directives*). An unknown setting name is an error.

**`CONVERT` unit strings.** The units are strings (`CONVERT('lbm/hr', 'kg/s')`) or, as in EES, unit
names written without quotes (`CONVERT(lbm/hr, kg/s)`, see below). A unit string is a product of unit
names with `-` or `*` between them, a quotient with `/` (everything after a `/` is in the denominator:
`kJ/kg-K` is kJ/(kg·K)), `^n` for an integer power (`m^2`, `K^-1`), parentheses and the number
1 (`1/s`): `lbm/hr`, `W/m^2-K`, `Btu/hr-ft^2-R`, `ft^3/min`, `kg/(m-s)`. The result is the
product of the factors, checked for the dimension:

```
m1 = CONVERT('lbm/hr', 'kg/s')          {1.259978805556e-4}
h  = CONVERT('W/m^2-K', 'Btu/hr-ft^2-R') {0.17611}
cp = 4.18*CONVERT('kJ/kg-K', 'Btu/lbm-R')
```

- **Unit names**: the SI units with their usual prefixes (`m`, `mm`, `cm`, `km`, `g`, `kg`, `mg`,
  `s`, `ms`, `Pa`, `kPa`, `MPa`, `mPa`, `bar`, `mbar`, `J`, `kJ`, `MJ`, `W`, `kW`, `MW`, `N`, `kN`, `L`,
  `mol`, `kmol`, `A`, `V`, `Hz`, ...), the English and other units of EES (`lbm`, `lb`, `lbf`, `oz`,
  `slug`, `ft`, `in`, `yd`, `mi`, `psi`, `psia`, `atm`, `torr`, `mmHg`, `inHg`, `Btu`, `cal`, `kcal`,
  `kWh`, `hp`, `gal`, `qt`, `acre`, `lbmol`, `mph`, `knot`, `min`, `hr`, `day`, `rpm`, `rev`, `rad`,
  `deg`, `cP`, `poise`, ...) and the long names `meter`, `inch`, `second`, `hour`, `joule`, `watt`, `pascal`,
  `kilogram`, `gram`, `liter`. A name is read with its exact capitalization first (`mPa` and `MPa` differ,
  `mW` and `MW`); else without regard to case when only one unit fits (`KPA`, `BTU`, `Lbm/Hr`); `mpa` is
  ambiguous and is an error.
- **Temperatures** inside a unit are temperature *differences* (K = C, R = F = 5/9 K), so `kJ/kg-K` and
  `kJ/kg-C` are the same unit; a bare `CONVERT('C', 'F')` is the ratio of the degrees, 1.8 (assumption: the
  EES help is not available to check it). The conversion of a temperature *value*, with its offset, is
  `CONVERTTEMP`.
- **Errors** (never a silent factor): an unknown unit name (`CONVERT('lbm/hrr','kg/s'): unknown unit 'hrr'
  in 'lbm/hrr' (did you mean 'hr'?)`; `kg/m3` is one too, a power is written `m^3`), a unit string that
  cannot be read (`kg//s`, `m^x`, `kg s`, an empty string), units of different dimensions (`CONVERT('kJ','Pa'): incompatible units
  (energy vs. pressure)`). When both units are string literals the error is a parse error with the
  line; for a unit computed at run time (a string variable of a FUNCTION) it is an evaluation error. A
  number of arguments other than two is an evaluation error. The factor of a pair of units is computed
  once.
- **Quotes are optional.** The first two arguments of `CONVERT` and `CONVERTTEMP` may be written
  without quotes, as in EES: `CONVERT(lbm/hr, kg/s)`, `CONVERT(kPa, Pa)`, `CONVERT(kJ/kg-K, Btu/lbm-R)`,
  `ConvertTemp(C, K, T)`. An argument made of unit characters only (letters, digits, `-`, `*`, `/`, `^`, `.`,
  parentheses and spaces, with at least one letter) is read as the unit text *as written*, before it becomes an
  expression (`lbm/hr` is not a division of two variables, `kJ/kg-K` is not a subtraction), and is checked like
  the quoted unit: `CONVERT(foo, Pa)` is a parse error with the line (`CONVERT('foo','Pa'): unknown unit 'foo'`).
  Quoted and unquoted arguments can be mixed (`CONVERT(kPa, 'Pa')`). When the unquoted text is also the name
  of a variable of the model (`C = 5` then `ConvertTemp(C, K, T)`), it is the unit that is meant: the
  variable is not read. A string variable (`u$`) or a string function is a run-time unit as before. The
  unit names are not variables of the model (no unknown, nothing in the LaTeX report). An argument with
  other characters (`a+b`, `x[1]`) is an expression, not a unit.

**`CONVERTTEMP` scales.** `CONVERTTEMP('from', 'to', T)` converts the temperature value `T` (with the offset
of the scales, unlike `CONVERT`) and accepts the four temperature scales of EES, `C`, `K`, `F` and `R`, or
their long names `celsius`, `kelvin`, `fahrenheit` and `rankine`; the name is read with its exact
capitalization first, else without regard to case (`c`, `Kelvin`, `FAHRENHEIT`), and spaces around it are
ignored: `CONVERTTEMP('C', 'F', 100)` is 212, `CONVERTTEMP('F', 'R', 212)` is 671.67
(100 C = 373.15 K = 212 F = 671.67 R). `T` can be an unknown of the system: the derivative with respect
to it is the ratio of the degrees of the two scales (1, 1.8 or 5/9), so a block can be solved through
`CONVERTTEMP`. Any other unit is an error that names it (`CONVERTTEMP: 'Foo' is not a temperature unit
(expected C, K, F or R, or celsius, kelvin, fahrenheit or rankine)`), a unit of another quantity (`'kPa'`) or
a compound unit (`'C/s'`) included. When a unit is a string literal the error is a parse error with the line;
for a unit computed at run time (a string variable of a FUNCTION) it is an evaluation error. A number of
arguments other than three is an evaluation error. As for `CONVERT`, the quotes around the units are
optional: `ConvertTemp(C, K, T)`.

---

## 3. Thermophysical property functions (CoolProp-backed)

All thermophysical functions take a fluid as the first positional argument
and state specifications as **named arguments** (typically two state
variables); they are the functions that accept `name=value` arguments (see
*Function arguments and indices* in §1). Example:

```ees
h_1 = enthalpy('Water', T=T_1, P=P_1)
rho = density('R134a', T=10, P=500e3)
```

Inputs and outputs follow the **fixed unit system** of CoolSolve:
temperature is in **°C** (never K), pressure in Pa, energy in J/kg, entropy in
J/(kg·K). See the warning at the top of this page.

### Standard properties

| Function                             | Result (CoolSolve units)   |
|:-------------------------------------|:---------------------------|
| `temperature` / `T`                  | °C                         |
| `pressure` / `P`                     | Pa                         |
| `enthalpy` / `H`                     | J/kg                       |
| `entropy` / `S`                      | J/(kg·K)                   |
| `internalenergy` / `U` / `intenergy` | J/kg                       |
| `density` / `rho` / `D`              | kg/m³                      |
| `volume` / `V`                       | m³/kg (= 1/ρ)              |
| `quality` / `x` / `Q`                | –                          |
| `cp` / `c_p` / `specheat`            | J/(kg·K)                   |
| `cv` / `c_v`                         | J/(kg·K)                   |
| `viscosity` / `mu`                   | Pa·s                       |
| `conductivity` / `k`                 | W/(m·K)                    |
| `speed_of_sound` / `soundspeed`      | m/s                        |
| `molarmass` / `MM`                   | kg/kmol, as in EES (no state needed); also the atomic weight of an element symbol (below) |

### Molar mass of a chemical element

`MOLARMASS` also takes the **symbol of a chemical element**, written like a fluid (`molarmass(C)`, `molarmass('C')`,
a string variable): it returns the standard atomic weight in kg/kmol, from the IUPAC 2005 table
(M. E. Wieser, *Pure Appl. Chem.* **78** (2006) 2051, table 2), which are the values EES returns (C 12.0107,
H 1.00794). It is the usual way to build a fuel C<sub>m</sub>H<sub>n</sub> from its mass composition:

```ees
C% = 86.9                                  "mass percentage of carbon"
m = (C%/100)/MolarMass(C)                  "kmol of C per kg of fuel: 0.0723522"
n = ((100-C%)/100)/MolarMass(H)            "kmol of H per kg of fuel: 0.129968"
M_CH4 = molarmass(C) + 4*molarmass(H)      "16.0425 kg/kmol"
```

The table holds every element that has a standard atomic weight (hydrogen to bismuth except the radioactive
elements that have none, Tc, Pm, Po, At, Rn, Fr, Ra and Ac, and Th, Pa and U), the symbol is read in any case, and
the result is a constant. A name that is a
**registered fluid** is the fluid, not the element: `molarmass(N2)`, `H2`, `O2`, `CO` (so the symbol `Co` of
cobalt is `CO`, carbon monoxide), `He` and `Ar` (the same value as the element) keep the molar mass of the molecule
or of the fluid. A symbol that is neither a fluid nor an element is the *Unknown fluid* error, and so is an element
in any other property function (`enthalpy(C, T=25)`) or with state inputs (`molarmass(C, T=25)`): CoolSolve has no
thermodynamic data for the elements as substances (`C`, `H`, `O`, `N` atoms).

### Saturation properties

| Function                             | Inputs                     |
|:-------------------------------------|:---------------------------|
| `t_sat` / `tsat`                     | `P=...` (returns saturation temperature) |
| `p_sat` / `psat`                     | `T=...` (returns saturation pressure)    |

### Inputs below the usual range (vacuum, cryogenic fluids)

Pressures below 1 kPa (`entropy(Water, T=25, P=100)`, `T_sat(Water, P=500)`) and cryogenic fluids below 50 K
(`enthalpy(Hydrogen, T=-250, P=1E5)`, `enthalpy(Helium, T=-268, P=1E5)`) are valid inputs. An input that does
not depend on the unknowns of its block (a number, a variable computed by an earlier block, the internal
pressure of an ideal-gas call) reaches CoolProp unchanged, in the solve and in the final verification. If
CoolProp cannot evaluate it (liquid water at -100 °C, Hydrogen at 10 K), the call is an error that names it
(*CoolProp evaluation failed in final solution*).

An input that depends on the unknowns of its block is a Newton trial point: during the iterations it is brought
back to a lower bound, so that a trial point with a negative pressure gives a finite residual instead of an
error. The bound is 1000 Pa for a pressure, 50 K for a temperature and 1e-4 kg/m³ for a density, or the minimum
of the fluid in CoolProp when that is lower (triple-point pressure: 611.655 Pa for `Water`; minimum temperature:
13.957 K for `Hydrogen`, 2.1768 K for `Helium`). If the solution of the block sits at the bound, the call was not
evaluated at the solution and the solve ends with an error that names the call, the input and the bound:

```
enthalpy(Hydrogen, T, P): temperature 10 K is below the minimum temperature of Hydrogen in CoolProp
(13.957 K); an input that depends on the unknowns of its block is limited there during the iterations, so
the solution would not satisfy this equation
```

### Additional EES-compatible thermophysical properties (priority 2)

| Function                               | Result                      |
|:---------------------------------------|:----------------------------|
| `PRANDTL`                              | Prandtl number (-)          |
| `SURFACETENSION`                       | N/m                         |
| `KINEMATICVISCOSITY`                   | μ/ρ, m²/s                   |
| `THERMALDIFFUSIVITY`                   | k/(ρ·cp), m²/s              |
| `COMPRESSIBILITYFACTOR`                | Z = pv/(RT) (-)             |
| `ISENTROPICEXPONENT`                   | isentropic exponent (-)     |

### Pure-fluid constants

These take *only* the fluid (no state):

| Function                 | Result                                     |
|:-------------------------|:-------------------------------------------|
| `T_CRIT`                 | Critical temperature                       |
| `P_CRIT`                 | Critical pressure                          |
| `V_CRIT`                 | Critical specific volume (m³/kg)           |
| `T_TRIPLE`               | Triple-point temperature                   |
| `P_TRIPLE`               | Triple-point pressure                      |
| `ACENTRICFACTOR`         | Acentric factor (-)                        |

### Phase inspection (string-valued)

`PHASE$(Fluid, T=..., P=...)` returns one of:

`liquid`, `gas`, `twophase`, `supercritical`, `supercritical_gas`,
`supercritical_liquid`, `critical_point`, `unknown`, `not_imposed`.

### Ideal-gas test

`ISIDEALGAS(Fluid)` returns **1** if CoolSolve models the substance as an ideal gas, **0** otherwise. It
selects the arguments of property calls, since a real fluid needs two state inputs and an ideal gas needs only
the temperature for `enthalpy`, `cp`, `cv`, `viscosity`, `conductivity` and `prandtl`:

```ees
PROCEDURE gam(Fluid$, Temp, Pres : Gamma)
  IF isidealgas(Fluid$) = 1 THEN
    Gamma = cp(Fluid$, T=Temp)/cv(Fluid$, T=Temp)
  ELSE
    Gamma = cp(Fluid$, T=Temp, P=Pres)/cv(Fluid$, T=Temp, P=Pres)
  ENDIF
END
CALL gam('N2', 20, 1E5 : g_ideal)            "N2 is an ideal gas"
CALL gam('Nitrogen', 20, 1E5 : g_real)       "Nitrogen is the real fluid"
```

- The ideal gases are the substances of the ideal-gas group of §5 (`Air`, `N2`, `O2`, `H2`, `He`, `Ar`, `CO2`,
  `CO`, `CH4`, `C2H6`, `C3H8`, `C4H10`, `C5H12`, `C6H14`, `C7H16`, `C8H18`, `H2O`); the answer is read from the
  fluid registry. Every other known fluid gives
  0: the real fluids (`Water`, `R134a`, `Nitrogen`, `CarbonDioxide`, `Air_ha`, …), the incompressible substances, the
  CoolProp mixtures and the humid air `AirH2O` (not an ideal-gas substance; the answer 0 is to be confirmed against EES).
- The fluid is written as in a property call: a quoted name, an unquoted name (`isidealgas(Air)`), a string variable
  or a string function such as `LOOKUP$('Species', i, 1)`, in any case. An unknown fluid is an error (*Unknown fluid*),
  as in a property call.
- The function takes the fluid and nothing else (no `T=` or `P=`). It is a constant of the fluid, with derivative 0,
  usable in the main program and in the bodies of functions and procedures.
- It only tells the substances apart: what a property call accepts is unchanged (the single-input call of an ideal
  gas is limited to the temperature-only properties of §5: `enthalpy`, `cp`, `cv`, `viscosity`, `conductivity`,
  `prandtl`; `density(Air, T=T)` still needs `P=`).

---

## 4. Humid-air / psychrometric functions

Humid air is accessed through the pseudo-fluid `airH2O` (or equivalently
`Air_ha` when treated as real dry air). Three state inputs are required
(typically `T`, `P`, and one of `W`, `R`, `Tdp`, or `Twb`).

| Function                | Output                                 |
|:------------------------|:---------------------------------------|
| `humrat` / `W`          | Humidity ratio kg_w/kg_da              |
| `relhum` / `R`          | Relative humidity (0..1)               |
| `wetbulb` / `B`         | Wet-bulb temperature                   |
| `dewpoint` / `D`        | Dew-point temperature                  |
| `enthalpy`              | Specific enthalpy per kg dry air       |
| `entropy`               | Specific entropy per kg dry air        |
| `intenergy`             | Specific internal energy per kg dry air (= h − P·v) |
| `volume`                | Specific volume m³/kg dry air          |
| `density`               | 1/`volume`, kg dry air per m³ (EES; the density of the mixture is `density`·(1+`humrat`)) |
| `cp`, `cv`              | Specific heats per kg dry air (CoolProp's `C` and `CV`; `cv` < `cp`) |
| `viscosity`, `conductivity` | Dynamic viscosity and thermal conductivity of the mixture |
| `temperature`           | The dry-bulb temperature given the other three inputs (e.g. `H`, `P`, `W`) |
| `thermaldiffusivity`    | k/(`density`·`cp`) of the mixture, m²/s |

Example:

```ees
W_1   = humrat  ( 'airH2O', T=25, P=101325, R=0.5 )
Twb_1 = wetbulb ( 'airH2O', T=25, P=101325, R=0.5 )
```

**Validity range.** CoolProp's humid-air model is limited to temperatures (also `Twb`, `Tdp`) from
−143.15 to 350 °C, pressures from 10 Pa to 10 MPa, `W` from 0 to 10 and `R` from 0 to 1; it gives no
value outside, and none for the properties that CoolSolve has no humid-air output for (`soundspeed`,
`prandtl`, `compressibilityfactor`, `kinematicviscosity`, the saturation and critical properties…).
Such a call is an error: during the iterations it has the placeholder value 10 000, so that
a trial point outside the range is a bad point and not a failure of the solver, but a solution that
contains it is refused (exit status 1, no `.sol` file) with *CoolProp evaluation failed in final
solution*, which names the call, its inputs in your units and the limit that is exceeded:

```
HumidAir error in enthalpy(): CoolProp returned invalid result (NaN or Inf) for enthalpy(AirH2O)
with inputs P=800000 Pa, T=1026.85 C, W=0: T=1026.85 C is outside the validity range of CoolProp's
humid-air model (-143.15 C to 350 C). For dry air above 350 C use the ideal-gas fluid Air.
```

For dry air above 350 °C (gas turbines) use the ideal-gas substance `Air`. The relative humidity of
supersaturated air (`W` above the saturation value) is not an error: `relhum` returns `W/W_sat`, which
is greater than 1.

### PSYCHPROPS built-in procedure

```ees
CALL psychprops(T_1, P_1, R_1 : T, v, h, s, u, W, R, Twb, Tdp)
```

Returns all nine psychrometric properties at the specified state
(`T` in °C, `P` in Pa, `R` = relative humidity in [0..1]). Outputs are
returned in CoolSolve's fixed unit system (temperature in °C, see the warning
at the top of this page). Inputs outside the validity range above are an error that names the input
and the limit.

---

## 5. Fluids

Fluids are grouped into the following categories. Names are
case-insensitive. Use the fluid name as the first positional argument of a
property function: `density('Water', T=25, P=101325)`. The fluid is a string literal, an
unquoted name (`Water`), a string variable (`fluid$`) or any **string expression**, such as
`LOOKUP$('Species', i, 1)` (§11.8), which is evaluated at each call.

### Real fluids (equation-of-state based)

Water / Steam / R718; Ammonia / R717; Propane / R290; Methane; Ethane;
Nitrogen / N2 (also as ideal gas); Oxygen / O2; Argon; CarbonDioxide / R744;
Hydrogen / H2; Helium; Toluene, Acetone, Benzene; refrigerants R11 – R41,
R113 – R161, R218, R22 – R32, R125 – R152a, R227ea, R236ea/fa,
R245ca/R245fa, R365mfc, R404A, R407C, R410A, R507A, R134a, R1233zd(E),
R1234yf, R1234ze(E/Z), R1243zf, RC318, SES36; alkanes n-Butane / IsoButane /
n-Pentane / IsoPentane / n-Hexane / IsoHexane / n-Heptane / n-Octane /
IsoOctane; cyclics CycloHexane / CycloPentane; butenes, HeavyWater, Deuterium,
DimethylEther / DiethylEther, EthylBenzene, Ethylene, Fluorine, HydrogenChloride,
HydrogenSulfide, Krypton, MM / MDM / MD2M / MD3M / MD4M, Neon, Neopentane,
NitrousOxide, Novec649, OrthoHydrogen / ParaHydrogen, Propylene,
SulfurDioxide, SulfurHexafluoride, Xenon, Air_ha (real dry air).

### Ideal gases

`Air`, `N2`, `O2`, `H2`, `He`, `Ar`, `CO2`, `CO`, `CH4`, `C2H6`, `C3H8`,
`C4H10`, `C5H12`, `C6H14`, `C7H16`, `C8H18` (the normal alkanes n-butane to n-octane, as vapor,
at a very low pressure only), `H2O` (as vapor, low pressure only).

The properties that depend on temperature only (`enthalpy`, `intenergy`,
`cp`, `cv`) can be called with `T` alone, as in EES — typically for mean
specific heats. Every spelling of these functions takes `T` alone, as it takes `T` and `P`
(`intenergy`, `internalenergy`, `u`; `h`; `specheat`, `specificheat`, `c_p`, `c`; `c_v`; `t`):

```ees
c_bar_p = (enthalpy(CO2,T=T_p1) - enthalpy(CO2,T=T_p2)) / (T_p1 - T_p2)
```

The single input must fix the temperature: `T`, or `h` or `u` (`cp(N2, h=h_1)`,
`temperature(N2, h=h_1)`). Any other single input (`d`, `v`, `s`, `x`, `P`) does not fix the
state, for any property: the call is an error that names it (*cp(Air, d=...): an ideal gas needs the
temperature (or h or u), or two inputs such as T and P / d and P*). Give two inputs, such as
`cp(Air, d=1.2, P=101325)`.

CoolSolve evaluates them at an internal reference pressure (101325 Pa; 100 Pa
for `H2O`, which keeps the vapor phase down to the triple point, 0.01 °C (CoolProp has no state of water
below it: an error); 0.01 Pa for the alkanes `C4H10` to `C8H18`,
see below). This pressure is part
of the model: it is never clamped by the solver and is identical in the solve
and in the final verification. `intenergy(CO2, T=25)` is the call with `T=25, P=101325`, and
`h − P/ρ` of that state: the internal energy has the reference of the enthalpy (the enthalpy of
formation for the species of the table below), and since CoolSolve evaluates the gases with CoolProp's
real-gas equations of state, `u = h − Z·R·T` with `R = R̄/M` and the compressibility factor `Z` at the
internal pressure (at 25 °C: `Z` = 0.9950 for `CO2`, 0.9998 for `N2`, 0.9997 for `Air`, 0.99994 for `H2O` vapor
at 100 Pa): the ideal-gas value `h − R·T` to that factor, not to the last digit.

The transport properties of the dilute gas depend on the temperature only too:
`viscosity`, `conductivity` and `prandtl` (alias `mu`, `k`, `l`) can be called with `T`
alone (or with an `h` or `u` that fixes it), at the same internal pressure, with the
derivative with respect to `T`:

```ees
mu_air = viscosity(Air, T=25)        "= 1.845e-5 Pa-s"
k_n2 = conductivity(N2, T=25)        "= 0.0258 W/m-K"
Pr = prandtl(Air, T=25)              "= 0.707"
```

The value is CoolProp's at the internal pressure, a few tenths of a percent from its
dilute-gas limit (`viscosity(Air, T=25)`: 0.08 % above CoolProp at 100 Pa, 0.001 % from its value at 1 bar;
`conductivity(CO2, T=25)`: 0.3 % above the limit). A model that needs the property at another
pressure gives it (`viscosity(Air, T=T, P=P)`: with `T` and `P` the call is unchanged, for ideal gases
too). What depends on the density needs `P` and stays an error with `T` alone (*Ideal gas property
'kinematicviscosity' requires pressure input*): `density`, `entropy`, `kinematicviscosity`,
`thermaldiffusivity`, `soundspeed`. The limits are CoolProp's: `CO` has no transport model, and the transport
properties of `H2O` vapor are defined from its triple point (0.01 °C) up; below, the call is an error
(*CoolProp returned invalid result*).

**Reference states of the enthalpy and of the entropy.** As in EES (JANAF convention), the
ideal gases `N2`, `O2`, `H2`, `He`, `Ar`, `CO2`, `CO`, `CH4`, `C2H6`, `C3H8`, `C4H10`, `C5H12`, `C6H14`,
`C7H16`, `C8H18` and `H2O` (vapor) have absolute properties from **one source**, the NASA-Glenn thermodynamic
data (McBride, Zehe and Gordon, NASA/TP-2002-211556; `thermo.inp` of NASA CEA), the
data family of the EES ideal gases: `h(25 °C)` is the standard **enthalpy of formation**
(the elements `N2`, `O2`, `H2`, `He`, `Ar` are 0) and `s(25 °C, 1 bar)` the standard
**third-law entropy** `s°`.

| Species | `h·M` at 25 °C [kJ/kmol] | `s°` (25 °C, 1 bar) [kJ/kmol-K] | Species | `h·M` at 25 °C [kJ/kmol] | `s°` [kJ/kmol-K] |
|:--|--:|--:|:--|--:|--:|
| `N2` | 0 | 191.609 | `CO` | −110 535.2 | 197.659 |
| `O2` | 0 | 205.148 | `CH4` | −74 600 | 186.370 |
| `H2` | 0 | 130.680 | `C2H6` | −83 851.5 | 229.220 |
| `He` | 0 | 126.153 | `C3H8` | −104 680 | 270.313 |
| `Ar` | 0 | 154.846 | `H2O` (vapor) | −241 826 | 188.828 |
| `CO2` | −393 510 | 213.786 | `C4H10` | −125 790 | 309.880 |
| `C5H12` | −146 760 | 349.558 | `C6H14` | −166 920 | 388.848 |
| `C7H16` | −187 780 | 428.087 | `C8H18` | −208 750 | 467.347 |

The enthalpies of formation are the values of the species headers of `thermo.inp`; `s°` is
the 9-coefficient polynomial of the same file at 298.15 K (R = 8.314462618 J/mol-K). EES uses the
same data with R = 8.314, which puts its values 6e-5 lower (EES help: `CO2` −393 486,
`H2O` −241 811, `CH4` −74 595 kJ/kmol); NIST-JANAF entropies differ from NASA's by at most
0.12 kJ/kmol-K (`CH4`). `Air` is not in the table: it keeps CoolProp's reference states (see below).
Earlier CoolSolve versions used the rounded JANAF values of Çengel & Boles (Table A-26),
5 to 830 kJ/kmol away (`C2H6` and `C3H8` the most), and `He` and `Ar` kept CoolProp's
enthalpy reference.

**The alkanes `C4H10` to `C8H18`** are the normal alkanes n-butane, n-pentane, n-hexane,
n-heptane and n-octane (the NASA-Glenn species `C4H10,n-butane` … `C8H18,n-octane`, evaluated with
the CoolProp fluids `n-Butane` … `n-Octane`; the first interval of their polynomials starts at 300 K, so
`s°` is extrapolated by 1.85 K, as for `C2H6` and `C3H8`). The real fluids `n-Butane`, `n-Octane`, …
stay real fluids with CoolProp's reference states. These vapors condense at room conditions (n-octane
at 25 °C: saturation pressure 1.87 kPa; n-butane: 243 kPa), so the calls with `T` alone (`enthalpy`,
`intenergy`, `cp`, `cv`, `viscosity`, `conductivity`, `prandtl`) are evaluated at an internal pressure of
**0.01 Pa**, below the triple-point pressure of the five fluids (0.08 to 2.1 Pa): the state is a gas at every
temperature CoolProp can evaluate, from the triple point up (−138.3 °C `C4H10`, −129.7 °C `C5H12`,
−95.3 °C `C6H14`, −90.6 °C `C7H16`, −56.8 °C `C8H18`), and the real-gas residual is below 1e-6 of the
ideal-gas value (`h − u = R·T` to 0.05 J/kg). Below that temperature the call is an error
(*CoolProp returned invalid result*), never the enthalpy of the liquid. Their enthalpy rise from 25 °C
follows the specific heat of CoolProp's ideal gas, which is within 0.7 % of NASA's (EES) up to 1500 °C
(specific heat at 25 °C: within 0.7 %, `C8H18`), except `C5H12` above 500 °C (1.2 % at 1000 °C, 3.2 % at
1500 °C).

With `T` **and** `P` given, the call is CoolProp's real-fluid state at that pressure, as for `H2O`
(EES gives the ideal gas at any pressure): `C8H18` at 25 °C is a liquid above 1.87 kPa, and
`enthalpy(C8H18,T=25,P=1E5)`, `entropy(C8H18,T=25,P=1E5)` are the properties of the liquid. Evaluate the
entropy of the vapor at a pressure below its saturation pressure and add `−R ln(P/P_low)`
(`entropy(C8H18,T=25,P=100)` is `s° + R ln(1000)`).

**Enthalpy.** The enthalpy and internal energy of the species include the enthalpy of
formation, so heats of reaction can be written with species enthalpies:

```ees
LHV_CO = enthalpy(CO,T=25)*molarmass(CO) + 0.5*enthalpy(O2,T=25)*molarmass(O2) - enthalpy(CO2,T=25)*molarmass(CO2)
"= 282 974.8 kJ/kmol (in J/kmol: 2.829748e8)"
```

The reference applies to all calls, including `h=` given as an input
(`temperature(CO2, P=P, h=h)`); enthalpy differences are unchanged. The
enthalpies of `Air` and of the real fluids (`Water`, `CarbonDioxide`…) keep
CoolProp's reference states: EES's `Air` is not a JANAF species either.

**Entropy.** The ideal-gas entropy at 25 °C and **1 bar** is `s°`, and it depends on the
pressure as `−R ln(P/1 bar)`. Entropies of different species can therefore be combined:
entropy of reaction, Gibbs energy, equilibrium constant, exergy of combustion. Together with the
enthalpies of formation, the Gibbs energy `h − T·s` of the species is on the NASA-Glenn basis, e.g.
the equilibrium constants of `CO2 ⇌ CO + ½ O2` and `H2O ⇌ H2 + ½ O2` at 2000 K and 3000 K are those of the
NASA polynomials, to the difference between CoolProp's and NASA's ideal-gas specific heats
(0.0003 to 0.017 on `ln K`). For `He` EES documents a different entropy convention, zero at 25 °C and
101.3 kPa: CoolSolve gives the third-law value.

```ees
"entropy of reaction CH4 + 2 O2 -> CO2 + 2 H2O(g) at 25 C (3 mol of gas on each side: any pressure, 2000 Pa keeps H2O a vapor)"
dS = (molarmass(CO2)*entropy(CO2,T=25,P=2000) + 2*molarmass(H2O)*entropy(H2O,T=25,P=2000) - molarmass(CH4)*entropy(CH4,T=25,P=2000) - 2*molarmass(O2)*entropy(O2,T=25,P=2000))/1000
"= -5.2 kJ/kmol-K"
```

Entropies given as input (`temperature(CO2, s=s, P=P)`, `enthalpy(CO2, s=s, P=P)`) are on the
same absolute basis, and isentropic calculations are unchanged. CoolSolve evaluates the
species with CoolProp's real-gas equations of state, so `entropy(X, T, P)` is the
third-law ideal-gas value plus CoolProp's residual entropy: it is 0.02 (`N2`) to 0.3 (`C3H8`)
kJ/kmol-K below `s°` at 1 bar and tends to the ideal-gas value at low pressure (EES gives
the ideal-gas value at any pressure). `H2O` is a vapor for CoolProp only below its saturation
pressure: at 25 °C and 1 bar `entropy(H2O,T=25,P=1E5)` is the entropy of liquid water, so evaluate
the vapor at a lower pressure (at least 1 kPa) and add `−R ln(P/P_low)`.

**`Air` keeps CoolProp's reference state for the entropy** (as for the enthalpy). EES's `Air` is not
a third-law species either ("`Air` does NOT use the JANAF table reference states", EES help): its
entropy differs from CoolSolve's by a constant (+1815 J/kg-K for EES), so absolute entropies of
`Air` cannot be compared with EES or combined with those of the species above; entropy differences
and isentropic relations are right. For combustion, use a mixture of `O2` and `N2` instead
of `Air`, as the EES help recommends. `Air_ha`, `AirH2O` and the real fluids (`Water`,
`CarbonDioxide`, `Nitrogen`…) also keep their own CoolProp reference states.

### Humid air

`airH2O` — requires three state inputs.

### Incompressible substances

CoolSolve's incompressible substances are the liquid **solutions** of CoolProp's incompressible library (`INCOMP::…`):
no pure incompressible is registered. The EES solids `Aluminum` and `Copper` are not available (CoolProp has no solid
properties, see *Unsupported EES fluids* below).

**Solutions** (require a concentration as the *second positional argument*
in mass %, or the named `C=...` argument):

| Alias / primary name          | Description              |
|:------------------------------|:-------------------------|
| `EthyleneGlycol` / `EG` / `MEG` | Ethylene glycol/water   |
| `PropyleneGlycol` / `PG` / `MPG`| Propylene glycol/water  |
| `CaCl2` / `MCA`               | Calcium chloride brine   |
| `NaCl` / `MNA`                | Sodium chloride brine    |
| `LiCl` / `MLI`                | Lithium chloride brine   |
| `Methanol` / `MMA`            | Methanol/water           |
| `EthanolSolution` / `MEA`     | Ethanol/water            |
| `Glycerol` / `MGL`            | Glycerol/water           |
| `LiBr`                        | Lithium bromide/water solution (CoolProp, Pátek-Klomfar fits, 0–75 %, 0–227 °C; density, cp, h, s; no viscosity/conductivity) |
| `SeaWater`                    | Sea water, salinity in mass % (0–12 %, 3.5 % = 35 g/kg), Sharqawy et al. 2010 via CoolProp, 0–120 °C, liquid only; density, cp, viscosity, conductivity, h, s |

Example:

```ees
rho = density('MEG', 30, T=10, P=101325)   "30 % MEG at 10 °C"
cp  = cp('EthyleneGlycol', T=10, P=101325, C=30)
cp_w = cp(EG, T=50, C=30)                  "the EES form: T and the concentration alone, 3802.5 J/kg-K"
rho_sw = density(SeaWater, 3.5, T=20, P=101325)   "sea water, salinity 3.5 % (35 g/kg): 1024.9 kg/m3"
```

**The concentration can be an unknown of the model.** `rho = density(MEG, C, T=20, P=1E5)` ⏎ `rho = 1040` solves for
`C` (31.39764 % in 4 iterations), with the second positional argument as well as with `C=`, `conc=` or `concentration=`,
and for every solution of the table. The concentration reaches CoolProp **at full precision** (not rounded), and the
derivative of the property with respect to it is part of the Jacobian. The derivatives with respect to the other inputs
(`T`, `P`, …) are CoolProp's analytical ones; CoolProp has **none with respect to the concentration** of an
incompressible, so that one is a **central finite difference** in the mass fraction (step 1E-6, one-sided at the ends of
the fluid's range: 0–60 % for `MEG`, 0–75 % for `LiBr`, 0–12 % for `SeaWater`), taken only when the concentration depends
on the unknowns of the block: a constant concentration costs nothing more. During the iterations a concentration
outside 0–100 % is the usual error asking for a concentration, and one inside that but above the range of the fluid is a
CoolProp failure (a penalty value for the iterations, an error in the final check of the solution).

**The EES form: the temperature and the concentration alone.** EES evaluates an incompressible substance with the
temperature as its only state input and the concentration `C` in mass %: `cp('EG',T=T,C=30)`,
`density(EG,T=T,C=30)`, `viscosity(...)`, `conductivity(...)`. CoolSolve accepts it for every substance of this section
(the concentration as the second positional argument, `C=`, `conc=` or `concentration=`; any fluid spelling, a string
variable too). The call is evaluated at an **internal pressure of 101325 Pa**: its value is exactly that of the same call with
`P=101325`, with the same derivative with respect to `T` (the pressure is a constant of the model, never clamped by
the solver). The pressure matters for few outputs: CoolProp's `density`, `cp`, `cv`, `viscosity`, `conductivity` (and
what is built on them: `prandtl`, `kinematicviscosity`, `thermaldiffusivity`, `volume`) of an incompressible do **not**
depend on it, so the value is the EES one whatever the pressure; `enthalpy`, `entropy` and `intenergy` do, through a
`p/ρ` term (for 30 % ethylene glycol at 50 °C the enthalpy rises by 0.8 J/kg per kPa): the T-only value is the one at
101325 Pa, and a model that needs another pressure writes `P=`. Any other single input (`P`, `h`, `s`, `d`, …) does not
fix the state and stays the error naming the call (*cp requires exactly 2 input properties … got 1*); `pressure` and the
saturation functions of an incompressible are not evaluated at the internal pressure (no pressure is ever returned as a
result); without a concentration the call is the usual error of the solutions. A real fluid still needs two state
inputs. The named `C=`, `conc=` and `concentration=` are never read as a state input by the unit hints (no warning
for `C=30`).

**`LiBr` is the aqueous solution, not a pure fluid**: it needs its concentration like the other solutions
(`density(LiBr, 50, T=40, P=1E5)` or `C=50`; without it the call is the usual error asking for one). CoolProp has the
density and the specific heat of the solution, and its enthalpy, entropy and internal energy integrated from the specific
heat; it has no viscosity and no conductivity (`viscosity(LiBr, 50, …)` is an error naming the call). The concentration
is limited to 0–75 % (an error naming the call above it). The solution has **no enthalpy of mixing**: its enthalpy and
entropy are integrated from `cp` at each concentration, with a reference state that is not that of `H_LIBR`, so they
cannot be used in the energy balance between streams of different concentration (the absorber and the generator of an
absorption machine): use `H_LIBR` and `X_LIBR` for that (next section). CoolSolve warns about it (code `C004`, in the CLI
and in the GUI console) when the enthalpy, the entropy or the internal energy of `LiBr` is an output or an input of a
call, once per run; the density and the specific heat give no warning.

**`SeaWater` is a salinity solution, not a pure fluid**: it needs its salinity like the other solutions, in mass %
(`density(SeaWater, 3.5, T=20, P=1E5)` or `C=3.5`; 3.5 % is 35 g/kg, a mean sea water; without it the call is the
usual error asking for a concentration). The properties are those of CoolProp's `INCOMP::MITSW` (the correlations of
Sharqawy, Lienhard and Zubair 2010, MIT seawater library): density (1024.9 kg/m³ at 3.5 % and 20 °C), specific heat,
viscosity and conductivity, and the enthalpy, entropy and internal energy integrated from the specific heat. The
validity range is a salinity of 0–12 % (an error naming the call above it), a temperature of 0–120 °C and the **liquid
phase only**: CoolProp refuses a state whose pressure is below the saturation pressure of the solution, so that the
limit at a given pressure is the boiling point of the solution (3.5 %: 100.3 °C at 1 bar, 100.6 °C at 101325 Pa, the
pressure of the temperature-only form above); a higher temperature needs a higher pressure (`P=5E5`). The vapour
pressure of sea water is not available (`p_sat` and `t_sat` of an incompressible are errors). As for every
incompressible solution, the enthalpy and the entropy are zero at 20 °C and 101325 Pa for **each** salinity (no
enthalpy of mixing): differences of enthalpy between streams of different salinity (desalination, brine mixing) are
not energy balances. The salinity is CoolSolve's convention for the solutions (mass % as second positional argument
or `C=`); the EES syntax for sea water was not checked, no compatibility with it is claimed.

### CoolProp mixtures

A **mixture** is written in CoolProp's syntax `A[x]&B[y]`: two or more components, each followed by
its **mole fraction** in brackets. It is used as any fluid, as a literal or in a string variable,
and the properties are computed by CoolProp's mixture model:

```ees
fl$ = 'R245fa[0.5]&R134a[0.5]'
h_1 = enthalpy(fl$, T=50, P=5E5)
T_bub = temperature(fl$, P=5E5, x=0)    "bubble point, 32.9 °C"
T_dew = temperature(fl$, P=5E5, x=1)    "dew point, 46.2 °C: a glide of 13 K"
M = molarmass(fl$)                      "mole-weighted, in kg/kmol: 118.04"
```

- The mole fractions are plain numbers (no expressions, no variables) and must sum to 1 (within
  1e-6); the composition is part of the string. An optional `HEOS::` prefix is accepted. Other
  backends (`REFPROP::…`) are not.
- The components are CoolProp fluid names (the names of the registry above, such as `R718` or `Steam`,
  are accepted too, in any case). A malformed string, fractions not summing to 1, an unknown component
  (the error names it) or a pair for which CoolProp has no binary interaction parameters is an
  evaluation error that says why.
- A component with a zero fraction is dropped, as CoolProp does: `'R245fa[1]&R134a[0]'` **is** pure
  R245fa (same values, no warning). This makes a composition parameter usable down to its end points.
- All the functions of this section work on mixtures with the inputs that CoolProp supports for them
  (`T,P`, `T,x`, `P,x`, …), with analytical derivatives where CoolProp provides them and finite
  differences inside the two-phase glide, as for the pure fluids. A property or input pair that CoolProp
  cannot compute for a mixture (e.g. `SURFACETENSION`) is reported as for any property error.
  `T_SAT` / `P_SAT` and the fluid constants (`T_CRIT`, …) are not defined for a mixture (use the
  bubble and dew points above).
- **Warning.** Mixture properties are *predictive*: they come from CoolProp's mixture model, with
  binary interaction parameters, not from data of the mixture itself, and may be unreliable,
  especially in two-phase (bubble/dew) states and for transport properties (viscosity,
  conductivity). CoolSolve says so with a warning (code `C003`, in the CLI and in the GUI console),
  once per mixture and run.
- The saturation dome of the GUI diagrams is not drawn for mixtures.

### LiBr-H2O solutions (EES ASHRAE routines)

The lithium bromide - water solution of the absorption machines (vapour pressure, enthalpy of mixing) is not a CoolProp
fluid; only its density and specific heat are (the solution `LiBr` of [Incompressible substances](#incompressible-substances)).
CoolSolve has the four routines of the EES *ASHRAE* library, which take a unit string as first argument (`'SI'`) and two
numbers:

| Routine                  | Returns                                                          |
|:-------------------------|:-----------------------------------------------------------------|
| `X_LIBR('SI', T, P)`     | mass % of LiBr of the solution at `T` in equilibrium with the vapour at `P` |
| `H_LIBR('SI', T, X)`     | specific enthalpy of the solution at `T` and `X`                 |
| `T_LIBR('SI', P, X)`     | temperature of the solution of concentration `X` in equilibrium with the vapour at `P` |
| `P_LIBR('SI', T, X)`     | vapour pressure over the solution of concentration `X` at `T`    |

```ees
x_abs = X_LIBR('SI', T_abs, P_ev/1000)      "58.31 mass % at 40 °C and 0.8726 kPa"
h_abs = H_LIBR('SI', T_abs, x_abs)*1000     "107.2 kJ/kg, in J/kg"
```

- **Units are fixed by the routines**, not by the unit system of the model: `T` in °C, `P` in **kPa**, `X` in
  **mass %**, `h` in **kJ/kg** (the models write `X_LIBR('SI', T, P/1000)` and `H_LIBR(...)*1000` for a model in Pa
  and J/kg). Only the unit string `'SI'` is accepted, in any case and as a constant, a string variable or a string
  function; `'Eng'` (English units, accepted by EES) and any other string are an error.
- **Correlations.** The 1989 ASHRAE Handbook of Fundamentals (the enthalpy fit is that of McNeely): the Dühring
  relation `t = B(X) + t'·A(X)` (A and B cubic polynomials of `X`) between the temperature `t` of the solution and
  that of the refrigerant `t'`, the vapour-pressure equation `log10(p) = C + D/T' + E/T'^2` of the refrigerant
  (`T'` in K, `p` in kPa) and the enthalpy polynomial `h = Ah(X) + t·Bh(X) + t^2·Ch(X)`. `T_LIBR` and `P_LIBR`
  evaluate them directly (the inverse of the vapour-pressure equation is a quadratic in `T'`), `X_LIBR` solves the
  Dühring relation for `X` on its increasing branch (the interval around 60 % where dt/dX > 0, within 0...100 %).
  They are implemented in `src/libr_ashrae.cpp` (no CoolProp call).
- **Validity ranges** quoted by EES: `X` 45...75 % for the Dühring relation, 40...70 % and 15.6...165.6 °C for the
  enthalpy. They are not enforced: as in EES, the routines evaluate outside them (the stored EES solution of the
  `water_libr` model uses `H_LIBR` at 74.8 %), the result being an extrapolation of the polynomials.
- **Agreement with EES.** `H_LIBR` reproduces the EES values to 1e-9. `X_LIBR` is the exact root of the ASHRAE
  equations: the EES values of the stored solutions (58.30348 and 74.77572 %) are 6.7e-5 and 1.3e-5 below
  it (58.30739 and 74.77668 %), as if EES used a refrigerant temperature 0.007 K higher (origin not identified).
  The `water_libr` example gives `COP` = 0.74665 against 0.74670 for EES (the water properties come from CoolProp).
  The enthalpies of `H_LIBR` are those of the ASHRAE fit; the models add them to those of `ENTHALPY(Water, …)`
  in the energy balances of the cycle as EES does.
- **Errors**, naming the routine and its arguments (*X_LIBR('SI', T = 300, P = 0.8726): no concentration of the
  Duhring relation gives this temperature …*): a pressure that is not positive (or beyond 10^7.05 kPa, where the
  vapour-pressure equation stops), a refrigerant temperature below absolute zero or `A(X) = 0` in `P_LIBR`, a
  temperature for which `X_LIBR` finds no concentration on the branch, non-finite arguments, a wrong number of
  arguments, a number instead of the unit string and named arguments (`T=…`).
- **Derivatives** are analytical (the implicit function theorem for `X_LIBR`), so `T`, `P` or `X` may be unknowns;
  an implicit solve that leaves the validity of `X_LIBR` stops with the error above, an initial guess near the
  solution avoids it.
- A `FUNCTION` of the model with the name of a routine replaces it.
- **Density and specific heat** of the solution (which the routines above do not give, `V_LIBR` and the transport
  routines being missing) can be taken from the CoolProp solution `LiBr` of
  [Incompressible substances](#incompressible-substances): `density(LiBr, X, T=T, P=P)`, `cp(LiBr, X, T=T, P=P)`, with
  `X` in mass %. Its enthalpy and entropy have no enthalpy of mixing and another reference state than `H_LIBR`
  (warning `C004`): keep `H_LIBR` for the energy balances.
- **Not implemented yet** (a call is an error that names the routine and its register item, see
  [Model library support](model_library_support.md) §4): `V_LIBR`, `Q_LIBR`, `VISC_LIBR` and `COND_LIBR`
  (`CS-GAP-FLUIDS-ABS`); the `*_LiBrH2O` functions of the Pátek-Klomfar library and `Q_LiBrH2O` (`CS-GAP-LIBR-PK`).
  The ammonia-water procedure and fluid `NH3H2O` are implemented (next section). The `*_LiBrSSC` library of EES (Yuan & Herold 2005) is
  not implemented either (`CS-GAP-LIBR-PK`; its names are unknown functions). The fluid `LiBrH2O` of the property
  calls (`ENTHALPY(LiBrH2O, …)`) is an error that points to the routines above.

### Ammonia-water (`NH3H2O`)

`NH3H2O` is the EES fluid of the ammonia-water mixture, with its property functions and the procedure
`CALL NH3H2O(Code, In1, In2, In3 : T, P, x, h, s, u, v, q)`. CoolSolve evaluates it with the formulation of
Tillner-Roth & Friend (1998), adopted by IAPWS as Guideline G4-01 (2001), written in `src/nh3h2o.cpp` and
`src/nh3h2o_flash.cpp` (no CoolProp call: CoolProp has no ammonia-water pair). **EES uses another model** (Ibrahim &
Klein 1993, a Gibbs-energy formulation whose coefficients are not published), so the results agree within the
difference between two models (see *Agreement with EES*), not to the digit.

**Fluid form.** The state is fixed by **three** properties among `T`, `P`, `X` (the **mass fraction of ammonia**: `X` is not
the quality here), `Q` (the vapour **mass** fraction), `H`, `S`, `U`, `V` (specific volume):

```ees
h = Enthalpy(NH3H2O, T=T1, P=P1, X=x1)          "sub-cooled, saturated or superheated: Q tells which"
x = MassFraction(NH3H2O, T=50, P=1e6, Q=0)      "ammonia mass fraction of the saturated liquid"
T = Temperature(NH3H2O, P=P2, X=x2, H=h2)
```

Functions: `Enthalpy`, `Entropy`, `IntEnergy`, `Volume`, `Density`, `Temperature`, `Pressure`, `Quality`, `MassFraction`,
`SpecHeat` (or `Cp`) and `VolExpCoef`. Inputs and outputs are in the units of the model (as for any fluid); `V` follows
the density unit. The fifteen combinations of the EES procedure are accepted, written here as the inputs in the order
of the table: `T P X`, `T P Q`, `T X V`, `T X Q`, `T H Q`, `T S Q`, `T U Q`, `T V Q`, `P X H`, `P X S`, `P X Q`, `P H Q`,
`P S Q`, `P U Q`, `P V Q` (the order of the arguments in the call is free). Any other combination, or a number of inputs
other than three, is an error that lists them. `Quality` is **-0.001** for a sub-cooled liquid, **1.001** for a
superheated vapour (EES convention) and between 0 and 1 inside the dome; a `Q=` input must be in [0, 1].

**`CALL NH3H2O`.** `Code` is the three-digit code of the inputs, the numbers of their properties (T = 1, P = 2, x = 3,
h = 4, s = 5, u = 6, v = 7, q = 8) in ascending order: `123`, `128`, `137`, `138`, `148`, `158`, `168`, `178`, `234`, `235`,
`238`, `248`, `258`, `268`, `278`. `In1`, `In2`, `In3` are the three properties of the code, in this order. The units of the
procedure are fixed whatever the unit system of the model: **T in K, P in bar, x mass fraction, h and u in kJ/kg, s in
kJ/(kg K), v in m3/kg, q** as above. All eight outputs are returned (fewer may be asked, in the order given).

**Formulation and reference state.** The ideal-gas parts of the pure fluids are fitted so that the internal energy and
the entropy of each pure saturated liquid vanish at its triple point (water 273.16 K, ammonia 195.495 K); EES uses
other reference states. `h`, `s` and `u` therefore differ from EES by a quantity that depends on the composition only:
**compare differences at a given composition** (or between states of the same composition), not absolute values.

**Agreement with EES** (the two formulations are independent, so this is the order of magnitude of their
difference): the EES help example `T = 450 K, P = 10 bar, x = 0.5` gives `v` = 0.2037 m3/kg (EES 0.2036) and quality
1.001; the enthalpy difference between that state and `T = 10 °C, P = 1000 kPa, x = 0.5` is 2432 kJ/kg (EES 2419, 0.6 %); the
saturated liquid at 50 °C and 1000 kPa has the mass fraction 0.5895 (EES 0.5953, 1 %; the guideline states an
uncertainty of 0.01 in the compositions). The implementation reproduces the verification tables of the guideline: the
single-phase values (Table 6) to 3e-9 for the Helmholtz energy and 6e-8 for p, c_v and the speed of sound, the bubble and
dew points (Tables 7 and 8) to the digits of the tables, and pure water to the 1 part in 60 000 the guideline quotes.

**Numerics.** The derivatives are exact (the Helmholtz energy is evaluated with truncated Taylor arithmetic in `T`, `rho` and
`x`; the flashes use the implicit function theorem on the converged equations, no finite difference), so any of the inputs
may be an unknown of the model. The equations of every flash are solved to a scaled residual below 1e-12, followed by a
full Newton polish; an iteration that cannot go lower is accepted only when it is within two units of the noise floor
measured at that point (mainly a composition within 1e-5 of pure: `1 - x` of a double is known to 1e-16 / (1 - x) only).
A call costs about 0.1 to 2 ms once warm (the previous solution is the start of the next call; a first call, or a jump,
costs up to some tens of ms).

**Limits.** The formulation covers the range from the triple-point line to the critical locus and pressures up to 40 MPa
(EES: 230 to 600 K, 0.2 to 110 bar). There is no two-phase state above the critical locus (an ammonia-rich mixture above
about 405 to 440 K): `Q` inputs there, and `P`/`T`/`X` inputs that need a saturation, are errors that name the state; a
single-phase state is still given. `T`, `P`, `X` and `Q` of a trial point of the solver are limited (T below 150 K, P below
1000 Pa, X and Q outside [0, 1]) as for CoolProp inputs, and a limit still active at the solution is refused.
`SpecHeat` and `VolExpCoef` are defined for single-phase states only. For the combinations `T H Q`, `T S Q`, `T U Q`, `P H Q`,
`P S Q`, `P U Q` the property along the saturation curve is not monotone in the composition (heat of mixing), so two
solutions can exist; the one of **lowest ammonia content** is returned.

**Not implemented.** The transport properties of EES (`Conductivity`, `Viscosity`, `SurfaceTension`: Cuenca 2013, Conde 2004)
are errors that say so (register `CS-GAP-NH3H2O`, remainder); the properties are those of the list above.
A `PROCEDURE NH3H2O` of the model replaces the built-in.

### Unsupported EES fluids

Some EES fluids are not modelled by CoolProp. Attempting to use them
triggers a clear runtime error with the fluid name and a remediation
suggestion:

`LiBrH2O`, `Acetylene`, `Chlorine`, `Chloroethene`,
`EthyleneOxide`, `Helium3`, `Potassium`, `Sodium`, `Aluminum`, `Copper`.

`Aluminum` and `Copper` are EES solids of its property library; CoolProp has no solid properties, so their
properties are not available (the call is an error naming the fluid, never a NaN): give them as constants or as a
correlation in the model (`CS-BUG-INCOMP-SOLIDS`).

The ammonia-water mixture `NH3H2O` is supported (previous section). The
lithium-bromide-water solution is available through the routines `X_LIBR`,
`H_LIBR`, `T_LIBR` and `P_LIBR` above, not as a fluid.

---

## 6. Control flow

### DUPLICATE loop

```ees
DUPLICATE i = 1, 3
    x[i] = i * 2
END
```

The loop is **expanded at parse time** — `DUPLICATE i = 1, N` produces N
copies of the body with `i` substituted.

The bounds of a loop of the main program must be **integers known before the
solve** (EES: "an integer number or a variable set to an integer constant"):

- a number or an arithmetic expression (`+ - * / ^`, parentheses) of numbers;
- a `$CONSTANT` (`DUPLICATE i = 1, imax#`, `DUPLICATE i = 1, (imax#-1)`);
- a variable of the main program set to a constant expression by an equation
  `n = <expression>`, above or below the loop, possibly through other such
  variables (`n = 2*m`, `m = 3`); the equation stays in the model;
- in a nested loop, an expression of the indices of the enclosing loops
  (`DUPLICATE j = i+1, n` inside `DUPLICATE i = 1, n-1`), evaluated for each copy.

A loop whose first bound is greater than the last runs downward
(`DUPLICATE i = 3, 1` gives `i` = 3, 2, 1).

A comment may follow the bounds on the header line (`DUPLICATE i = 1, 3 "loop over the
nodes"`, `// ...` or `{...}`), in the main program and in a body: its text is never read, a comma
or an `=` in it included (a `"..."` or `//` comment there used to be read as part of the last bound:
*DUPLICATE could not parse range expressions*).

An index with arithmetic is evaluated in each copy, in equations, function
arguments and `CALL` inputs and outputs, with one or several indices:
`T[i+1] = T[i] + 1` gives `T[2] = T[1] + 1`, `T[3] = T[2] + 1`, …, and
`a[i, j+1]` in nested loops names the element `a[1,2]`, `a[1,3]`, ….

Any other bound is an error that names it, with the line of the loop — a
variable set by an equation that depends on an unknown, by a `CALL`, or not set
at all, an array element, a value that is not an integer:

```
Line 2: DUPLICATE i: the last bound 'n' is not known before the solve: 'n' is set
by 'n = y + 1' (line 1), which is not a constant expression. A bound is a number,
a $CONSTANT, or a variable set to a constant expression by an equation of the main
program (n = 17), combined with + - * / ^
Line 2: DUPLICATE i: the last bound 'n' is 2.5, not an integer
```

A `DUPLICATE` inside a `FUNCTION` or `PROCEDURE` body is a loop executed when the
routine runs: its bounds may be the inputs of the routine or its local variables.

### REPEAT-UNTIL loop

```ees
REPEAT
    x = x + 1
UNTIL(x > 10)
```

`REPEAT ... UNTIL` is a loop of the body of a `FUNCTION` or `PROCEDURE` (§7): the
body runs, then the condition is tested, so the body runs at least once. The condition
is the whole text after `UNTIL`, written like the condition of an `IF` (§ IF-THEN-ELSE):
the parentheses around it are optional (`UNTIL x > 10` is `UNTIL(x > 10)`), it may use
`=`, `<>`, `<`, `>`, `<=`, `>=` (and the spellings `=>`, `=<`), compare strings
(`UNTIL s$ = 'done'`, without regard to case) and combine relations with `AND` and `OR`
(`UNTIL (a > 1) AND (b < 2)`); a value without a relation is true when it is greater than 0.5.
An `UNTIL` with no condition (nothing, a comment or `()` after it), or with one that
cannot be read, is an error (*"REPEAT missing UNTIL(condition)"*, with the line of the `REPEAT`).
A run of the
loop that has not ended after **100 000** iterations (the `UNTIL` condition never becomes true)
is an error that names the routine and the line of the `REPEAT`, and the model is not
solved (before this limit existed, the loop stopped silently at 10 000 iterations and the routine
returned the value it held then):

```
REPEAT loop in FUNCTION 'cnt' (line 3) did not end after 100000 iterations: an UNTIL condition that never becomes true?
```

The limit is per run of a loop: a loop nested in another one has its own count each time it
starts. A loop left by a `GOTO` (§7) is not concerned.

A comment may follow `REPEAT` and the `UNTIL` condition (`UNTIL x > 10 // done (or not)`): a `)` in it is
not the closing parenthesis of the condition.

As for `IF ... THEN` and `GOTO`, the loop does not exist in the main program, where the
equations are solved simultaneously and nothing is executed in order: a `REPEAT` there (also in the
body of a `DUPLICATE` of the main program) is a parse error with the line of the `REPEAT`, and the model
is not solved (such a loop used to be parsed and silently never executed). Write the iteration as
equations, which the solver iterates on, or put the loop in a `FUNCTION` or `PROCEDURE`:

```
Line 2: REPEAT ... UNTIL loops are only allowed inside FUNCTION and PROCEDURE bodies (in the main program an iteration is solved by the solver as equations; to run a loop, put it in a FUNCTION or PROCEDURE)
```

### IF-THEN-ELSE (inside functions and procedures)

As in EES, `IF` statements exist only in the body of a `FUNCTION` or
`PROCEDURE` (§7), where statements are executed in order. In the main program
the equations are solved simultaneously: use the inline `if(cond, a, b)` or
`IF(A, B, X, Y, Z)` function there (§ Inline selection), an `IF ... THEN`
statement is a parse error.

```ees
IF (T < 0) THEN f = 1 ELSE f = 2        "single line: no ENDIF"

IF (Re < 2300) THEN                     "block form: closed by ENDIF"
    f = 64/Re
ELSE
    f = 0.316*Re^(-0.25)
ENDIF

IF (Pr$ = 'DENSITY') THEN Pro = 2 ELSE  "single line ending with ELSE: the ELSE part follows"
    IF (Pr$ = 'SPECHEAT') THEN Pro = 3 ELSE
        IF (Pr$ = 'DYNVISC') THEN Pro = 5
ENDIF;ENDIF                             "one ENDIF per IF ending with ELSE"
```

- Only the branch selected by the condition is executed; the other one is
  never evaluated (a division by zero in an untaken branch is harmless).
- The condition compares numbers or strings with `=`, `<>`, `<`, `>`, `<=`,
  `>=` (EES also writes the last two `=>` and `=<`, which are the same operators:
  `IF (Z => 0.001) AND (W > 0) THEN`; strings are compared without regard to case, as in EES:
  `fluid$ = 'R134a'` also matches `'r134A'`), and combines
  relations with `AND`, `OR` and parentheses. The parentheses around the whole
  condition are optional.
- A branch holds any statements of a body, including other `IF` statements and
  `DUPLICATE` / `REPEAT` loops; `ELSE IF (c) THEN ...` is an `IF` statement
  inside the `ELSE` branch (a block-form `IF` always has its own `ENDIF`).
- The `;` that separates statements also separates `ENDIF`s: `ENDIF;ENDIF` closes
  two open `IF`s (one `ENDIF` per `IF`), and the text after an `ENDIF;` is read as
  the next line (`ENDIF; y = 1`, `ENDIF;ELSE`, `ENDIF;END`). Text after an `ENDIF`
  that is not separated by `;` is an error (*"Unexpected text 'y = 1' after ENDIF"*).
  An `IF` left open at the `END` of the body is an error (*"IF ... THEN without a
  matching ENDIF"*). An `ENDIF` that closes no `IF` is ignored in a body, with a
  warning (`P008`, one per line), because EES accepts it (the `IF` ladders of EES
  libraries often have one `ENDIF` too many); outside a body, and for an `ELSE` with
  no `IF`, it is an error (*"ENDIF without a matching IF ... THEN"*).
- A single-line `IF (c) THEN a ELSE` whose line ends with `ELSE` (only a comment
  may follow) has its `ELSE` part on the following lines, up to the `ENDIF` that
  closes it: this is how EES libraries write `IF` ladders (see the example; the
  last `IF` of the ladder has no `ELSE` and needs no `ENDIF`).
- The one-word `ELSEIF` is not an EES keyword: it is a parse error.
- Derivatives are propagated through the branch taken, so the function can be
  used on unknowns.

### GOTO and statement labels (inside functions and procedures)

A `GOTO` and the numeric **statement labels** it jumps to are supported in the body of
a `FUNCTION` or `PROCEDURE`, where statements are executed in order (like `IF`, §7); in
the main program (equations solved in no particular order) they are errors. Legacy EES
libraries use them for ladders and for loops:

```ees
PROCEDURE sumto(n : s)
    s = 0
    j = 0
90: j = j + 1                      "a label, followed by a statement on the same line"
    s = s + j
    IF (j < n) THEN GOTO 90        "backward jump: a loop"
END

PROCEDURE regime(t, d, w : q)
    q = 1
    IF (t <= d) THEN GOTO 10 ELSE GOTO 20
10:
    q = 2
    GOTO 1000
20:
    q = 3
1000:
END
```

- `GOTO n` is a statement: on a line of its own, after a `;`, in the branches of a
  single-line `IF` (`IF (c) THEN GOTO 10 ELSE GOTO 20`) and inside block `IF`s,
  `DUPLICATE` and `REPEAT` loops. The body goes on with the statement that follows the
  label `n:`.
- A label is an integer followed by a colon (`100:`, `10 :`), alone on its line or
  followed by a statement (`90 : j = j - 1`); comments may follow. A label at the end
  of a body (`1000:` before `END`) is valid. The numbers need not be ordered.
- The label must be in the same body, in the block of the `GOTO` or in a block that
  contains it: a jump **within** a block, **backward** (a loop) or **out of** nested `IF`
  branches and loops is valid. A jump **into** a block (a label inside an `IF` branch, a
  `DUPLICATE` or `REPEAT` body that does not contain the `GOTO`) is refused, as is a label
  that two statements define.
- A body that would execute more than 100 000 `GOTO`s in one call is stopped with an error
  that names the routine and the line (a loop built with `GOTO` that never ends); a
  `REPEAT ... UNTIL` loop has the same limit of 100 000 iterations per run (§ REPEAT-UNTIL loop).
- The read-before-assignment check of the bodies (§7) goes in the order of the text, a
  `GOTO` does not change it; derivatives are propagated through the path taken, as for
  an `IF`.

Errors, all with their line and not solved:

```
Line 3: GOTO 100 has no target: PROCEDURE 'p' has no statement label '100:'
Line 5: label 100 is defined twice in PROCEDURE 'p' (lines 3 and 5): a GOTO could not tell which one it means
Line 3: GOTO 10 would jump into a block: the label '10:' (line 5) is inside the THEN branch of the IF at line 4; a GOTO can only jump to a label of its own block or of a block that contains it
Line 2: GOTO is only allowed inside FUNCTION and PROCEDURE bodies (the statements of the main program are equations, solved in no particular order)
Line 4: statement labels are only allowed inside FUNCTION and PROCEDURE bodies (...)
```

---

## 7. User-defined functions and procedures

### Functions

```ees
FUNCTION area(width, height)
    area = width * height
END
```

Functions return a single numeric value through the variable named after
the function.

In the body of a function or a procedure every statement is parsed as in the
main program: a line that cannot be parsed is an error (*"Could not parse
line"*, with its line number), it is never skipped. Comments of every form,
blank lines and directives are not statements. The `END` that closes the
body may be followed by a comment (`END // area`). `GOTO` and numeric statement
labels (`10:`) are supported in a body, see §6 (*GOTO and statement labels*).

The `END` of a `FUNCTION` or `PROCEDURE` may also name the routine, as EES allows: the kind
(`END FUNCTION`, `End procedure`) or its own name (`END area`), in any case, with an optional
comment after it. A different name or kind (`END foo` in `PROCEDURE p`, `END FUNCTION` closing
a `PROCEDURE`) is the error *"END foo does not close PROCEDURE p"*, with its line. Inside the
body a `DUPLICATE` is still closed by a bare `END` (and a `REPEAT` by `UNTIL`): a named `END`
while a loop is open is an error too. `END MODULE` and `END SUBPROGRAM` stay unsupported
(MODULE and SUBPROGRAM are not handled).

```ees
PROCEDURE divmod(x, y : quotient, remainder)
    quotient  = FLOOR(x / y)
    remainder = x - quotient * y
END divmod

FUNCTION area(width, height)
    area = width * height
End Function
```

### Procedures

```ees
PROCEDURE divmod(x, y : quotient, remainder)
    quotient  = FLOOR(x / y)
    remainder = x - quotient * y
END

CALL divmod(17, 5 : q, r)
```

The `:` separates input arguments (left) from output variables (right).

#### Outputs of a CALL

An output is a variable (`CALL p(a : y)`, `CALL p(a : T[2])`) or, in the main
program, **an expression**, as in EES: the output of the procedure is equated to
the expression, and the system solves that equation like any other, for the
unknown it contains (or for an input of the procedure, when the expression is
known).

```ees
CALL BRINEPROP('Specheat', brine$, X, t : c_p/1000)   "c_p = 1000 x the output of the procedure"
CALL BRINEPROP('Dynvisc',  brine$, X, t : mu*1000)
CALL props(x : T+273.15, P_out/1000, rho/MM)   "three outputs, each equated to its expression"
```

- Every output position is independent: expressions and variables may be mixed
  in one `CALL`, and several outputs may be expressions.
- The expression may contain any unknown, also one that depends on the
  procedure's input (a coupled block): the Jacobian comes from the automatic
  differentiation of the procedure like for a variable output.
- An output that is a constant or has no unknown is an equation that must
  hold: with nothing else unknown the system is over-determined and is
  reported as such.
- Inside a `DUPLICATE` loop of the main program the iterator is substituted in
  the expression as in any equation (`CALL p(i : z[i]/10)`).
- **Inside a `FUNCTION` or `PROCEDURE` body an output must be a variable**: the
  statements run in order and nothing solves for the expression. An expression is the parse error *"Output 1 of 'p' is the
  expression 'r/1000': inside a FUNCTION or PROCEDURE …"*; write the result into
  a variable and compute the expression from it.
- Internally the output goes to an **auxiliary variable**
  `__out<k>_<procedure>_L<line>` (`k` = position of the output, `line` = source
  line of the `CALL`; `_2`, `_3`, … for the copies made by a `DUPLICATE` loop or
  when the model already uses that name) and the equation `expression = auxiliary`
  is added (`c_p/1000 = __out1_BRINEPROP_L12`). When the expression is one
  variable under `+`, `-`, `*`, `/` (by something without variables) and a unary
  minus, the equation is written solved for the variable (`c_p = __out1_… * 1000`,
  an explicit equation, like the one a user writes by hand); any other
  expression stays implicit (`sqrt(c_p) = __out1_…`). The auxiliary variable is a
  variable of the model like the others — it is in the solution (`.sol`, the
  variable tables) with the raw value of the procedure output, in the units of the
  procedure, and in the debug folder — and it exists only for expression outputs:
  a model with variable outputs only has none. Names that begin with `__out` are
  best left unused.

**String outputs.** A string variable (`name$`) set by a procedure is a result
like the numeric ones: it is available to the CALLs and the property functions
that follow (as the fluid name, `enthalpy(fluid$, …)`), whether or not the main
program also assigns it a literal, in the order the structural analysis finds
(also inside a coupled block).

### Bodies are assignments executed in order

As in EES, the body of a `FUNCTION` or `PROCEDURE` is **100 % explicit**: a
sequence of assignment statements executed from top to bottom (`=` and `:=` both
mean assignment there). Nothing in a body is solved simultaneously. The routine
itself may still be **called implicitly** from the main program: `ff(x) = 10`
solves `x` through the explicit function `ff`, and a `CALL` inside an algebraic
loop is iterated by the solver like any equation.

Three rules follow, checked once when the model is parsed (statements inside `IF`
branches, `REPEAT ... UNTIL` and `DUPLICATE` loops included). A violation is an
error with its line, like a syntax error, and the model is not solved:

- **The left-hand side is a variable or an array element.** Anything else is an
  implicit equation, which a body cannot solve:

  ```
  Line 3: implicit equation in PROCEDURE 'pq': the statements of a FUNCTION or PROCEDURE are
  assignments executed in order, so the left-hand side must be a variable. Solve the equation
  in the main program, or iterate explicitly in the body (REPEAT ... UNTIL)
  ```

  (A lookup-table cell write `LOOKUP('t', 1, 1) = x` is not concerned: it is no
  statement of the body, see §11.10.)
- **A variable is assigned before it is read.** The variables of a body other
  than its inputs are local and *"not initialised"* (EES help, *Internal
  Functions*): the outputs (the function's own name, the outputs of a procedure)
  and every other name have no value until a statement assigns them, and a body
  does not see the variables of the main program.

  ```
  Line 2: variable 'b' is read in PROCEDURE 'p' before it is assigned: the variables of a
  FUNCTION or PROCEDURE other than its inputs are local and not initialised (pass
  main-program values as input arguments)
  ```

  A variable counts as assigned from the first statement that assigns it, in
  the order of the text: an assignment in either branch of an `IF`, the output
  of an earlier `CALL` and one element of an array all count. The iterator of a
  `DUPLICATE` loop, the inputs, the built-in constants (`pi`, `g#`, `R#`,
  `sigma#`, …) and the `$CONSTANT`s are always available. Inside a `REPEAT` loop the text order is
  the rule too: a variable read before its assignment in the loop must also be
  assigned before the loop. Every read counts — right-hand sides, indices,
  `IF` and `UNTIL` conditions, the input arguments of a `CALL`, the arguments
  of functions — and string variables (`fluid$`) follow the same rule; the
  names that are not variables (function names, the `T=` keys of property
  calls, an unquoted fluid name such as `R134a`) are ignored, as in the main
  program. Each variable is reported once, at its first read.
- **Every output is assigned.** An output that no statement of the body assigns
  (not even in one branch of an `IF`, or as the output of a `CALL`) would have
  no value; the error is reported on the header line of the routine:

  ```
  Line 1: FUNCTION 'ff' never assigns its output 'ff': every output of a FUNCTION or
  PROCEDURE must be assigned by its body
  ```

The explicit rewrites, valid EES with the same results:

```ees
"A value of the main program: pass it as an input argument"
PROCEDURE p(a, b : y)       "was PROCEDURE p(a : y), reading the main-program b"
    y = 2*b + a
END

"An implicit equation: iterate explicitly in the body..."
PROCEDURE friction(Re, RR : f)
    f = 0.02
    REPEAT
        f_old = f
        f = (-2*log10(RR/3.7 + 2.51/(Re*sqrt(f_old))))^(-2)     "Colebrook, fixed point"
    UNTIL (abs(f - f_old) < 1E-12)
END
CALL friction(1E5, 0.001 : f)                                   "f = 0.022175"

"...or solve it in the main program and pass the result as an input argument"
1/sqrt(f) = -2*log10(RR/3.7 + 2.51/(Re*sqrt(f)))
CALL pipe(f, L, D : DELTAp)
```

A body that holds a syntax error is checked once that error is fixed.

### Recursion is not allowed

As in EES (*Internal Functions* of the EES help, which refuses a function that calls itself), a `FUNCTION`
or `PROCEDURE` cannot call itself, nor call a routine that calls it back: a **call cycle between the
routines of the model is a parse error**, with the routines along the cycle and the line of the call that
closes it, and the model is not solved:

```
Line 2: Recursive call: FUNCTION f calls itself; recursion is not allowed in EES
Line 5: Recursive calls: FUNCTION a → FUNCTION b → FUNCTION a; recursion is not allowed in EES
```

Every call counts, wherever it is written in a body: statements, `IF` branches and conditions, `DUPLICATE`
and `REPEAT` loops, the arguments of other calls (built-ins included), the input arguments and the output
indices of a `CALL`, the inline form of a procedure (`y = p(x)`). A recursion that would terminate (a guard
in an `IF`) is an error too. Calls to built-ins, to property functions and to names that the model does not
define are not part of the cycle search; a `FUNCTION` and a `PROCEDURE` of the same name are two routines;
a `PROCEDURE Colebrook` of the model replaces the built-in (it is then a routine like the others, so one that
calls itself is an error, whereas a routine that calls the built-in `Colebrook` is not).

Rewrite a recursion as a loop in the body (`DUPLICATE`, `REPEAT ... UNTIL`), or as a chain of routines that
call one another in one direction only. As a last line of defence for the callers of the evaluator that do
not go through the parser, the evaluator refuses more than 200 nested routine calls with an evaluation
error naming the routine, instead of overflowing the stack.

### `$COMMON` is not supported

In EES, `$COMMON a, b$, c` on the first line of a `FUNCTION` or `PROCEDURE` body
lets the body read variables of the main program without receiving them as
arguments. CoolSolve does not support it: a `$COMMON` line in a body is an
error (the run stops, in the CLI and in the GUI console), naming the variables,
the routine and the rewrite:

```
Line 2: $COMMON is not supported: pass b as an input argument of PROCEDURE pq (and in each CALL)
```

Pass the shared variables as input arguments, after the existing ones, and in every call;
the body does not change, and the model is valid EES with the same results:

```ees
PROCEDURE pq(a, b : y)      "was PROCEDURE pq(a : y) with $COMMON b"
    y = a + b
END

b = 21
a = 1
CALL pq(a, b : y)           "y = 22"
```

String variables (`b$`) are passed the same way. A `$COMMON` line in the main
program shares nothing and is ignored; one inside a comment or a string is text.

### CALL ERROR and CALL WARNING

`CALL ERROR` and `CALL WARNING` are the EES built-in procedures with which a `FUNCTION` or
`PROCEDURE` checks its inputs, usually in an `IF`. `CALL ERROR` stops the calculation and reports the
message; `CALL WARNING` reports a warning and the calculation continues. As in EES they are valid only
inside a `FUNCTION` or `PROCEDURE` body (in the main program they are a parse error) and have no outputs:

```ees
FUNCTION abc(X, Y)
  G := X + Y
  IF (X <= 0) THEN CALL ERROR('X must be greater than 0.  A value of XXXF1 was supplied', X)
  IF (G <= 0) THEN CALL WARNING('The sum of X and Y should be greater than zero.  The sum was XXXA1 ', G)
  abc := G
END
g = abc(3, -4)       "g = -1, with the warning above"
```

| Form | Message |
|:-----|:--------|
| `CALL WARNING('message XXXF1', value)` | the message with the value in place of `XXX` |
| `CALL WARNING('message')` | the message |
| `CALL WARNING(value)` | the generic text *A warning message was issued due to the value XXX in YYY.* (`YYY` is the routine) |
| `CALL ERROR(value)` | the generic text *Calculations have been halted because a parameter is out of range.  The value of the parameter is XXX in function YYY.* |

The message is a string constant, a string variable (`s$`) or a string function; the value is a number
or, written `XXX$`, a string. Where `XXX` is followed by a display-format code of EES the value takes
that format: `A<d>` automatic (the digit is ignored, as is the case without a code), `F<d>` fixed (`d` digits
after the decimal point), `E<d>` exponential (`d` digits after the point, exponent without padding:
`5.123E2`), `N<d>` `d` significant figures (`123` for 123.4567 with `N3`, `123.5` with `N4`; `E` format beyond
1E6 and below 1E-6), `G<d>` engineering (`d` significant figures, exponent a multiple of 3), `d` being a digit
0 to 9 or `M` (the most digits). The automatic format writes a whole number without decimals and any other number
with up to 10 significant digits. A message that has no `XXX` but a value gets the value appended in brackets
(`'X is positive [3]'`), so that it is not lost; the `||` line break of EES is not interpreted, and the time and
date formats are not supported.

**When a warning is reported.** A warning describes the solution, not the Newton iterates: an intermediate
iterate may cross the validity range of a correlation that the solution respects. CoolSolve therefore ignores
`CALL WARNING` while it iterates, and when a block has converged it evaluates the block once more at the
solution and reports the warnings of that evaluation, as a warning diagnostic with code `W001` and the line
of the `CALL WARNING` statement (`warning (line N): …` on the standard error of the CLI, and in the GUI
console). Two warnings with the same line and the same text are reported once, whatever the number of calls
that raised them. Only a block that calls a routine which can reach a `CALL WARNING`, directly or through
other routines, is evaluated again: a model without the statement runs exactly as before. This is the
`$Warnings ON /LAST` behaviour of EES; EES can also list the warnings of the earlier iterations (in gray, option
`/ALL`), which CoolSolve does not. The `$Warnings` directive does not change anything: the help of EES states
that a warning of the `WARNING` procedure is always displayed.

### CALL Colebrook (Darcy friction factor)

`Colebrook` is the procedure of an EES add-on library that EES models call for the Darcy friction factor `f`
of turbulent pipe flow. It is built in: it solves the implicit Colebrook-White equation

```
1/sqrt(f) = -2*log10(RR/3.7 + 2.51/(Re*sqrt(f)))
```

with `Re` the Reynolds number and `RR = e/D` the relative roughness. It has two inputs and one output, and
the EES files call it in either form (the name is case-insensitive):

```ees
CALL Colebrook(Re, RR : f)
Call Colebrook(Re_D, RR, f)      "the EES form without colon: the third argument is the output"
```

In the colon-less form the third argument must be a variable (or an array element); with the colon the output
follows the rules of any `CALL` (an expression is accepted in the main program, see *Outputs of a CALL*).
The procedure is usable in the main program and inside a `FUNCTION` or `PROCEDURE`: it iterates internally,
so for the body of a routine it is an explicit statement that assigns its output.

- **Turbulent flow only.** There is no laminar branch: the models switch themselves, e.g.
  `IF (Re_D <= 2000) THEN f = 64/Re_D ELSE CALL Colebrook(Re_D, RR, f)`.
- **Numerics.** Newton iterations on `x = 1/sqrt(f)` from the Swamee-Jain approximation, to full double
  precision (relative step below 1E-14); for instance `f = 0.0185139` at `Re = 1E5`, `RR = 1E-4`.
- **Errors.** `Re <= 0`, `RR < 0`, `RR >= 3.7` (the equation has no solution) and NaN stop the calculation with
  a message naming `Colebrook` (*Colebrook: the Reynolds number Re must be positive (and finite), got 0*); a
  call with another number of inputs or outputs is a parse error.
- **Derivatives.** `f` is an implicit function of `Re` and `RR`: its derivatives follow from the implicit
  function theorem at the solution (automatic differentiation), so `Re` or `RR` may be unknowns
  (`CALL Colebrook(Re, 1E-4 : f)` with `f = 0.02` solves for `Re`).
- A `PROCEDURE Colebrook` defined by the model (before its calls, as EES requires) replaces the built-in.

---

## 8. Arrays

Arrays use `[index1, index2, ...]` indexing and are resolved at solve
time:

```ees
DUPLICATE i = 1, 5
    T[i] = T_in + i * DeltaT
END
```

### Array ranges

`X[a..b]` stands for the elements `X[a], X[a+1], …, X[b]` (EES *array range
notation*), wherever a list of values is expected:

```ees
$CONSTANT imax# = 5
FUNCTION q(T, i, EPS[1..imax#])          "a header: 2 + 5 inputs, EPS is the array of the body"
  q = EPS[i]*T
END
PROCEDURE props(T[1..imax#] : EPS[1..imax#], tot)
  ...
END
CALL props(T_3[1..imax#] : eps_3[1..imax#], total)   "inputs and outputs"
q_3 = q(T_3[2], 2, eps_3[1..imax#])                  "a function argument"
T_bar = AVERAGE(T_3[1..imax#])                       "an aggregate built-in"
```

- **Where**: in the arguments of a function, of a user `FUNCTION` and of a
  `CALL` (inputs and outputs), in the main program and in the bodies, and in the
  header of a `FUNCTION` or `PROCEDURE`. Anywhere else (`y = x[1..3]`, an
  operand, an index, the left-hand side of an equation) it is an error; the EES
  assignment `A[1..3] = [1, 2, 3]` is not supported.
- **Arguments**: the range expands argument by argument, in order: the header
  `p(a, X[1..3], b)` declares 5 inputs, `a, X[1], X[2], X[3], b`, and in the body
  `X` is an array (`X[2]`, `X[i]` with a loop index). A call supplies the 5 values
  in any form — a matching range, single elements, expressions, a shorter range
  completed by other values (`CALL p(1, V[2..3], 7, 2 : …)`). A call that uses a
  range, or calls a routine whose header has one, must give as many values as the
  header declares: *"CALL pp gives 2 inputs (an array range counts its elements),
  but PROCEDURE pp (line 3) takes 3"* (outputs and function arguments alike).
- **Two dimensions**: a range on one index (`A[2,1..3]`) or on both
  (`A[1..2,1..2]`, the last index varying fastest: `A[1,1], A[1,2], A[2,1], A[2,2]`).
- **Bounds**: integers known before the solve, as for `DUPLICATE` (§6): numbers,
  `$CONSTANT`s, arithmetic of those and, in the main program, variables set to a
  constant expression (`n = 3`, then `x[1..n]`). The first bound must not exceed
  the last. In a header or a body the bounds are numbers or `$CONSTANT`s: a body
  does not see the variables of the main program and its own are known only when
  it runs (EES's variable-size header arrays, `FUNCTION f(n, A[1..n])`, are not
  supported). A bound that uses the index of an enclosing `DUPLICATE`
  (`sum(x[1..i])` inside the loop) is not supported. Each of these cases is an
  error with the line: *"array range 'x[3..1]': the first index (3) is greater
  than the last (1)"*, *"array range 'x[1..2.5]': the last index is 2.5, not an
  integer"*, *"array range 'X[1..n]' in FUNCTION 'g': the bounds of a range in a
  FUNCTION or PROCEDURE are numbers or $CONSTANTs …"*.

---

## 9. Initial guesses

Guess values live in a `.initials` file next to the model, one
`name = value` per line. The GUI and CLI auto-load it when present.

---

## 10. Configuration

See `examples/coolsolve.conf` for the full list of options. Key sections:

- `[coolprop]` — backend (`HEOS`, `BICUBIC&HEOS`, `TTSE&HEOS`), analytical
  derivatives, AbstractState caching.
- `[solver]` — pipeline order (Newton, TrustRegion, LM, BisectionND,
  Homotopy, Partitioned, Kinsol), tolerances, iteration limits.
- `[deepsearch]` — `deepSearchPipeline` / `deepSearchPipelineMode` used by
  the GUI **Try Harder** button (defaults to the full sequential chain).
- `[multistart]` — `multiStartMode` (`always` / `deepsearch` / `never`),
  `multiStartMaxRestarts`, `multiStartNumCores` (default 4). The legacy
  `multiStartEnabled` boolean is still accepted.
- `[tearing]` — SCC-based block reduction options.

> **Note on units.** The unit system itself is **not** configurable through
> `coolsolve.conf`: temperature is always in °C, pressure in Pa, energy in J.
> See the warning at the top of this page.

---

## 11. Lookup tables

Lookup tables let you interpolate from external CSV data files directly inside
equations.  CoolSolve loads companion tables for a model automatically by
scanning the same directory for CSV files whose name follows the convention:

**`<modelname>-<tablename>.csv`**

The part after the first hyphen is the **table name** that you reference in
equations.  For example, alongside `steam_cycle.eescode`:

- `steam_cycle-data.csv`    → callable as `LOOKUP('data', …)`
- `steam_cycle-watercp.csv` → callable as `INTERPOLATE('watercp', …)`

Files that do not match this pattern are ignored (including a bare
`<modelname>.csv`, since it would imply an empty table name).

Additional tables can be loaded via:

- **GUI Lookup Tables panel** — create, import, or edit tables directly in the
  browser.
- **ZIP bundle upload** — include any number of `.csv` files whose stems
  become the table names.
- **REST API** — `PUT /api/v1/tables/{name}` with CSV body.

If a lookup function references a table name that has not been loaded,
CoolSolve reports an error (`L002`) and the solve fails with a descriptive
message that names the table (and the `FUNCTION`/`PROCEDURE` the call is in,
see §11.7).

Each CSV must have a header row with column names (case-insensitive).

### 11.1 1D interpolation

```ees
h = INTERPOLATE('data', 'T_C', 'h_kJ_per_kg', T)
```

Performs linear interpolation of column `h_kJ_per_kg` as a function of
column `T_C` at value `T`.  The rows of the `T_C` column must be sorted, in
ascending **or descending** order (EES accepts both).  Outside the range the line
through the **two end nodes** continues (linear extrapolation, as EES does, with no
error and no warning): beyond the last node the value is the last one plus the slope
of the last interval times the distance, and the same below the first node. A table
with a single row gives its value everywhere.

Equivalent aliases: `INTERPOLATE1('table', 'xcol', 'ycol', x)`.

**EES call form.** EES writes the same interpolation with the input column named in
a last argument, `Column=value`:

```ees
M_r = INTERPOLATE('lookup_1', 'M_r', 'r_p', r_p=r_p)    "M_r where the column r_p equals the variable r_p"
```

The call lists the table and two columns; the named argument gives the **input
column** (it must be one of the two, whatever its position, case-insensitive) and
the value to look up (any expression, an unknown of a block included), and the
**other** column is returned. `INTERPOLATE('t', 'A', 'B', B=v)` therefore gives the
column `A` where `B` = `v` and is the positional `INTERPOLATE('t', 'B', 'A', v)`;
`INTERPOLATE('t', 'A', 'B', A=v)` gives `B` where `A` = `v`. `INTERPOLATE1` takes the
same form. Both forms use the same linear interpolation and derivative. The EES form
takes exactly the table, two columns and one named argument (not a fourth positional
value besides it), and the name must be a column of the call: otherwise the error
names the function, the table and the columns of the call.

> CoolSolve interpolates **linearly**; EES `INTERPOLATE` is cubic, so the results can
> differ slightly inside the table. The linear extrapolation through the two end nodes
> is what an EES solution shows beyond the last node (one library model, an operating
> point beyond its 3-node map); the same rule is assumed below the first node, which
> no EES solution has shown. Operating points far beyond the table are
> extrapolations: check them against the physics of the model.

Analytical derivative (for Newton/AD): slope of the active interval; outside the
table the slope of the end interval (not zero), so that an unknown whose value lies
beyond the table is solved with the exact derivative. At the end nodes it is the slope
of the end interval too; a single-row table has a zero derivative.

### 11.2 2D bilinear interpolation

```ees
z = INTERPOLATE2('mytable', 'x_col', 'y_col', 'z_col', x, y)
```

Performs bilinear interpolation of `z_col` on a grid defined by the unique
values of `x_col` and `y_col`.  The CSV rows may appear in any order;
CoolSolve constructs the grid automatically.  Requires at least 2 distinct
values in each axis.

Alias: `INTERPOLATE2DM(...)`.

EES form: the call lists the table and **three** columns and names the two input
columns, `INTERPOLATE2('mytable', 'x_col', 'y_col', 'z_col', x_col=x, y_col=y)`;
the column that is not named is returned, and the two named arguments may come in
either order (the two inputs are taken in the order of the columns of the call).

Analytical partial derivatives for the Newton/AD engine.

Outside the grid, `INTERPOLATE2` holds the values of the nearest edge (flat, zero
derivative): unlike the 1-D `INTERPOLATE` it does **not** extrapolate linearly (how EES
extrapolates a 2-D table has not been established).

### 11.3 Table metadata

| Function | Description |
|---|---|
| `NLOOKUPROWS('table')` | Number of data rows (excluding header) |
| `NLOOKUPCOLUMNS('table')` | Number of columns |
| `LOOKUPCOL('table', 'colname')` | 1-based column index of the named column |
| `LOOKUPCOL1('table', 'colname')` | Same as `LOOKUPCOL` |
| `LOOKUPCELLEMPTY('table', row, col)` | 1 if the cell is NaN/empty, 0 otherwise |

### 11.4 Direct cell access

```ees
P_at_row3 = LOOKUP('data', 3, 2)
```

Returns the value at the given row and column (both 1-based), matching EES
behaviour:

- **Non-integer row or column** — linear interpolation between the two
  adjacent rows/columns.  When both are non-integer, bilinear interpolation
  is applied across the four surrounding cells.
- **Out-of-range clamping** — if the row is below 1 the first row is
  returned; if it exceeds the number of rows the last row is returned.
  The same rule applies to the column.  The analytical derivative is zero
  in the clamped region (flat extrapolation).
- **Column by name** — the column argument may be a quoted string (column
  header name, case-insensitive) or a numeric index.

```ees
P_mid  = LOOKUP('data', 2.5, 2)             {midpoint between rows 2 and 3}
P_last = LOOKUP('data', 99, 2)              {clamped to last row}
h_r2   = LOOKUP('data', 2, 'h_kJ_per_kg')   {column by name}
```

Aliases: `TABLEVALUE('table', row, col)`, `TABLEVALUE#(...)`,
`TABLERUN#(...)`.

**Nameless form `LOOKUP(row, col)`.** Without a table name, `LOOKUP` reads the **default
table** of the model, as in EES where it is the first Lookup table (*Lookup 1*). The default
table is the companion table named `lookup_1` (`<modelname>-lookup_1.csv`, name matched
without regard to case: this is how an imported EES file stores its *Lookup 1*, see
[EES import](ees_import.md)); a model with no `lookup_1` but exactly one table reads that
table. With several tables and none named `lookup_1` the call is an error that names the file
to provide (`<model>-lookup_1.csv`) and lists the tables found. The form is the named one
with the table left out: the row (rounded or interpolated), the column (a number or a name),
the clamping, the derivatives and the use of a variable as the row or the column are the
same, and the form works in the main program, in `DUPLICATE`, and in `FUNCTION` / `PROCEDURE`
bodies and `INTEGRAL` models (§11.7). It exists for `LOOKUP` only: `LOOKUP$`, the table
metadata and the aggregate functions always take the table name.

```ees
DUPLICATE i = 1, 17
  T_cd[i] = LOOKUP(i, 1)          {same as LOOKUP('lookup_1', i, 1)}
  W[i]    = LOOKUP(i, 'W_dot')
END
```

### 11.5 Aggregate functions

Each function takes a table name and a column index (1-based integer or
column name string):

| Function | Description |
|---|---|
| `SUMLOOKUP('table', col)` | Sum of all non-NaN values in the column |
| `AVGLOOKUP('table', col)` | Arithmetic mean |
| `MAXLOOKUP('table', col)` | Maximum value |
| `MINLOOKUP('table', col)` | Minimum value |
| `STDDEVLOOKUP('table', col)` | Population standard deviation |

### 11.6 Example

For a model file `steam_cycle.eescode` with companion
`steam_cycle-data.csv`:

```
T_C,P_Pa,h_kJ_per_kg
100,101325,2675.6
150,476100,2745.9
200,1554900,2826.8
```

```ees
T = 150  "°C"
h = INTERPOLATE('data', 'T_C', 'h_kJ_per_kg', T)  "h = 2745.9 kJ/kg"
n = NLOOKUPROWS('data')                            "n = 3"
P_row2 = LOOKUP('data', 2, 2)                     "P_row2 = 476100 Pa (integer row)"
P_mid  = LOOKUP('data', 1.5, 2)                   "P_mid  = 288712 Pa (interpolated)"
P_last = LOOKUP('data', 99, 2)                    "P_last = 1554900 Pa (clamped)"
```

The **Lookup Tables** tab in the GUI lets you create, view, and edit tables
directly in the browser — no need to manage CSV files manually.

### 11.7 Lookups inside routines and `INTEGRAL` models

Every lookup function above can be called from the body of a `FUNCTION` or a
`PROCEDURE` (also from a routine called by another routine, in `REPEAT` /
`DUPLICATE` loops and in `IF` conditions), and from the equations of an
`INTEGRAL` model (§12), where the table is read at each integration step. This is
how EES libraries read their coefficient tables, e.g. the brine properties of the
ULiège `Brineprop.lib`:

```ees
PROCEDURE brinemini(Fl, Pr, Temp : output)
  k = 0
  REPEAT
    k = k + 1
    Row[k] = k + (Fl-1)*2
    c[k] = LOOKUP('Brine1', Row[k], Pr)     {coefficient rows of table Brine1}
  UNTIL (k >= 2)
  xm = LOOKUP('Brine2', Fl, 1)
  output = c[1] + c[2]*(Temp - LOOKUP('Brine2', Fl, 2)) + xm
END
```

The derivative of a lookup with respect to its row / interpolation argument is
analytic, so a routine output that is an unknown of a block still converges with
Newton. A wrong call is an error that names the table and the routine, e.g.
`lookup(): lookup table 'Brine1' not found (in PROCEDURE 'brinemini')`; `LOOKUP`
needs 3 arguments (or 2, `LOOKUP(row, col)` on the default table, §11.4), `INTERPOLATE` 4 (or 3 and a named argument, §11.1), `INTERPOLATE2` 6 (or 4 and two
named arguments, §11.2).

### 11.8 Text cells: `LOOKUP$`

`LOOKUP$('table', row, col)` returns the **text** of a cell (EES `LOOKUP$`). A table may
mix text and numbers: a cell that is not a number is text (`N2`, `kJ/kg-K`), and a CSV field
in double quotes may contain commas. This is how a model replaces a string array (see
*String arrays are not supported*, §2): the strings are rows of a table and the model keeps
the row number.

`species.eescode` with `species-Species.csv`:

```
Name,M
N2,28.0134
O2,31.9988
CO2,44.0095
```

```ees
DUPLICATE i=1,3
  h[i] = MolarMass(LOOKUP$('Species', i, 'Name'))*enthalpy(LOOKUP$('Species', i, 1), T=T)
END

PROCEDURE unitlabel(Pro : UO$)
  UO$ = LOOKUP$('Units', Pro, 1)      "the label of row Pro of the table Units"
END
```

- `row` is rounded to the nearest integer; `col` is a number or a column name. A row or a
  column outside the table, an unknown table or column, and a wrong number of arguments are
  errors that name the table (and the `FUNCTION`/`PROCEDURE` the call is in), e.g.
  *LOOKUP$(): row 4 is outside the table 'Species' (1 to 3 rows)*.
- A numeric cell gives its text as written in the CSV (`1.50` stays `1.50`, without the
  surrounding blanks) when the table has at least one text cell; an empty cell gives the
  empty string. In a table without any text cell the number is formatted with up to 15
  significant digits. `STRINGVAL(LOOKUP$(...))` converts a numeric cell to a number.
- The result is a string. It can be used as the fluid of a property function (`enthalpy`,
  `cp`, `MolarMass`, `T_CRIT`, `PHASE$`, ...), as an operand of `=` and `<>` (without regard
  to case), as the table or column argument of another lookup function and as the value of a
  string variable inside a `FUNCTION` or `PROCEDURE` body (`UO$ = LOOKUP$('Units', Pro, 1)`).
  A variable in `row` or `col` is a variable of the equation: it is computed first. In the
  main program `s$ = LOOKUP$('Species', 2, 1)` gives the string variable `s$` the text of the cell
  (*String variables of the main program*, §2). Used as a number it is an error.
- Numeric functions work as before on a table with text cells: `LOOKUP('Species', 2, 2)` is
  31.9988, `LOOKUP('Species', 2.5, 'M')` interpolates between the rows, and a text or empty
  cell (read as NaN) spoils the result only where it is really read. The row of a key is
  found with `LOOKUP$ROW` (§11.9).
- The GUI Lookup Tables panel shows and edits the text cells like the other cells (the table
  is kept as CSV text); no GUI setting is needed.

### 11.9 Rows found by key: `LOOKUP$ROW`

`LOOKUP$ROW('table', col, key)` returns the **row number** of the first row whose column `col`
holds `key` (EES `LOOKUP$ROW`). It is how a table with a key column is read: the row of the key is
found once, and the other columns are read at that row with `LOOKUP` / `LOOKUP$`.

`thermal_comfort.eescode` with `thermal_comfort-veture.csv`:

```
veture,I_cl
Tropicale,0.35
Ete,0.5
Interieur (legere),0.7
```

```ees
Veture$ = 'Ete'
row = LOOKUP$ROW('veture', 'veture', Veture$)      "row = 2"
I_cl = LOOKUP('veture', row, 'I_cl')              "I_cl = 0.5"
label$ = LOOKUP$('veture', row, 1)                 "label$ = 'Ete'"

PROCEDURE ReadPV(PVName$ : A_c)
  Row = LOOKUP$ROW('PVModules', 'PV Name', PVName$)   "a column name may contain spaces"
  A_c = LOOKUP('PVModules', Row, 'A_c')
END ReadPV
```

- `col` is a column name (case-insensitive, blanks allowed) or a number. `key` is a string (a
  literal, a string variable, or a string function such as `LOOKUP$` or `Uppercase$`) or a number.
- A **string** key is compared without regard to case (accented Latin-1 letters included, as `=`
  on strings) with the text of the cells, the text that `LOOKUP$` reads: `'HIVER'` finds `Hiver`,
  and `'2'` finds a numeric cell `2`. A **number** key is compared with the numeric value of the
  cells, equal up to a relative 1e-9 (a computed key such as `0.1*3` still finds `0.3`); a text or
  empty cell never matches a number. EES documents the string key (the help of EES was not
  available to check the numeric one, which CoolSolve also accepts).
- When several rows hold the key, the **first** is returned.
- A key that no row holds is an **error** that names the table, the column and the key, e.g.
  *LOOKUP$ROW(): no row of table 'veture' holds 'Nope' in column 'veture'* (and the
  `FUNCTION`/`PROCEDURE` the call is in), never a row 0. An unknown table or column, a column
  number outside the table, or a call without 3 arguments are errors too.
- The result is a number (an integer row, at least 1), usable wherever a row is: `LOOKUP`,
  `LOOKUP$`, `LOOKUPCELLEMPTY`. It is piecewise constant, so its derivative is zero: an unknown
  cannot be solved *through* `LOOKUP$ROW`, the key must be known when the call is evaluated.
  A string variable used as the key is a variable of the equation: the equation that gives it
  its string (a literal, an alias, a function, a `CALL` output) is solved first. In a
  `PROCEDURE` or `FUNCTION` body the key may be a string input.

### 11.10 Writing a cell: `LOOKUP('table', row, col) = value`

In EES a `LOOKUP` or `LOOKUP$` on the **left-hand side** writes a cell of a lookup table (a way to
export results). **CoolSolve does not support writing a cell and never does it silently**: the
results are in the solution, the GUI and the debug folder, not in a written table. The statement
(`=` or `:=`, in the main program, in a `FUNCTION` or `PROCEDURE` body, in an `IF` or a loop, in a
`DUPLICATE`) is

- **not an equation** and not a statement of the body: the model has the equations it has without
  it (the main program of `x = 5`, `y = 2*x`, `lookup('t', 1, 1) = x` has 2 equations in 2
  unknowns, it used to be refused as *"There are 3 equations and 2 unknowns"*), and a body does the
  same as without it. Its value is **not evaluated** (no side effect: an expression that fails, or
  a variable that nothing sets, does no harm);
- reported with a **warning** `P009` at its line, **one per source line** (not one per copy of a
  `DUPLICATE`, not one per statement of a line with `;`):

  ```
  warning (line 3): lookup('t', 1, 1) = ...: writing a lookup-table cell is not supported;
  the value is not stored (the results are in the solution)
  ```

- an **error** if the model **also reads** a table it writes (any lookup function on the same
  table name, anywhere: `LOOKUP`, `LOOKUP$`, `INTERPOLATE`, `NLOOKUPROWS`, the aggregates, ..., in
  the main program or in a body, before or after the write, names compared without regard to case;
  the nameless `LOOKUP(row, col)` is the table `lookup_1`), because the results would depend on a
  value that CoolSolve does not store. One error per table, at its first write, with the line of
  its first read:

  ```
  Line 2: lookup table 't' is written at line 2 and read at line 3: a lookup-table cell write is
  not supported (the value is not stored), so the read would not see it; remove the write or the
  read (the results are in the solution)
  ```

A table written and never read is only a warning, whether or not the model has a file for it: a
table that the model reads and does not write is read from its file as always. When the table is
named by a string variable (`LOOKUP(table$, row, 1) = x`) the name is not known when the file is
read, so only the warning is given (no check against the reads of the model). Only `LOOKUP` and
`LOOKUP$` on the left-hand side are cell writes: any other implicit statement keeps its
behaviour (an equation in the main program, an error in a body, *"implicit equation in ..."*).

---

## 12. Equation-based integration (dynamic/DAE solving)

CoolSolve can solve **initial-value differential–algebraic equation (DAE)
systems** written in the EES integral form, and tabulate the trajectory of
selected variables. This covers the vast majority of transient thermal models
(tank dynamics, heat exchangers, control loops, etc.).

The EES integral convention declares an ODE as an algebraic equation:

```ees
y = y0 + integral(dydt, t, t0, tf)     "declares a state variable y"
dydt = ...                              "its derivative (an algebraic expression)"
y0 = 1                                  "initial condition"
```

is equivalent to  `dy/dt = f(t, y, …)`, `y(t0) = y0`. When a model contains any
top-level `integral(...)` call, CoolSolve automatically routes it to the
dynamic solver — no flag is needed — and **non-integral models are completely
unaffected** (zero overhead).

### 12.1 `INTEGRAL` — a value, and the state variables it declares

```ees
y = <base> + INTEGRAL(integrand, t, lo, hi)
y = <base> + INTEGRAL(integrand, t, lo, hi, step)   "fixed-step form"
```

`INTEGRAL(integrand, t, lo, hi)` is a **value**: the integral of `integrand` with
respect to `t` from `lo` to the current value of `t`. It may appear anywhere in an
equation, as in EES — `y = m_dot*INTEGRAL(v, P, p_in, p_out, 1000)`,
`Q = INTEGRAL(Q_dot, tau, 0, 18000, 1)/3.6E6`, `W = -N*INTEGRAL(P, v, v1, v2)`,
`z = INTEGRAL(a, t, 0, 1) - 2*INTEGRAL(b, t, 0, 1)`. At the end of the run every
variable holds its value at `t = hi`.

| Argument | Meaning |
|:---------|:--------|
| `integrand` | The derivative (the `f` in `dy/dt = f`): a named variable defined by its own equation, or any expression (it is evaluated by the algebraic solver at each step, so array elements, property calls, lookup tables (§11.7) and user functions are allowed). |
| `t` | The integration variable (time, a pressure, a volume...). All `INTEGRAL` calls in a model must share the **same** `t` (names are not case-sensitive). |
| `lo`, `hi` | Lower and upper limits: constants, variables or expressions (§12.7). All calls must share the **same** `[lo, hi]`; `hi` must be greater than `lo`. |
| `step` | Optional 5th argument: the step size (a constant, a variable or an expression). `0`, or omitting it, lets the method choose (§12.5). |

A call that is the whole right-hand side, or is added to / subtracted from the rest
(`y = <base> + INTEGRAL(...)`), declares `y` a **state variable** owned by the
integrator: `<base>` is the non-integral part, and at `t = lo` the integral term is
zero, so `y(lo)` is determined by the algebraic solve of the remaining equations (no
separate initial-value extraction is needed). This is the fast path, and it applies
only when `<base>` is constant along the trajectory — a base that uses the integration
variable, another state, or any variable that depends on them is not an initial value,
and the equation is then treated like the forms below.

Any other form — a factor, a division, a function, two calls in one equation, a call on
the left-hand side, a subtracted call — is an ordinary algebraic equation in which each
call stands for its own **auxiliary state** (`I(lo) = 0`, `dI/dt = integrand`), created
internally. Auxiliary states are an implementation detail: they never appear in the
solution (`.sol`), in the trajectory table (`-integral.csv`, GUI) or in the variable
list. The variable on the left-hand side is then an algebraic variable like any other;
its column in the default table is its value along the trajectory.

Every unknown that is not a state is *algebraic* — it is solved from the remaining
equations at each time step by the same algebraic solver used for static models.

### 12.2 Coupled ODEs and algebraic variables

Any number of state variables may share the same integration variable and
interval (coupled ODEs). Algebraic variables may depend on the states and on
each other; together this is a **semi-explicit index-1 DAE**
`y' = f(t,y,z)`, `0 = g(t,y,z)`:

```ees
"Harmonic oscillator:  y' = z,  z' = -y"
y = 0 + integral(dydt, t, 0, 10)
z = 1 + integral(dzdt, t, 0, 10)
dydt = z
dzdt = -y
```

```ees
"A state coupled to an algebraic variable (heat-transfer cell)"
T = T0 + integral(dTdt, t, 0, 100)
dTdt = (T_amb - T) / tau        "derivative"
Q = h * (T_amb - T)             "algebraic variable, solved each step"
```

### 12.3 `$IntegralTable` — tabulating the trajectory

```ees
$IntegralTable t:0.1  y  z  dydt
```

The directive lists the variables to record, with the integration variable
first. Names are separated by spaces or, as in EES, by commas or semicolons
(`$IntegralTable t:0.1; y; dydt` and `$IntegralTable t:0.1, y, dydt` are the same
table). Names are not case-sensitive and the column header uses the spelling of the
model. An optional output interval follows the integration variable as
`name:interval` (`t:0.1` → one row every 0.1 units); the interval may also be a
variable that holds a constant (`t:dt_out`). If omitted, a row is
written at every step. Array ranges expand with the `..` notation:

```ees
$IntegralTable t  X[1..5]    "records t and X[1]..X[5]"
```

A name that is not a variable of the model is reported as a **warning** and its
column stays empty; a name listed twice is kept once (warning). Without a
`$IntegralTable` directive the table holds the integration variable, the variables
defined by the `INTEGRAL` equations and their integrands.

After a successful solve, CoolSolve writes the trajectory to a CSV named after
the model — **`<modelname>-integral.csv`** — next to the `.eescode` file. The
first column is the integration variable; subsequent columns are the
`$IntegralTable` variables in order. This file is regenerated on each run and
also appears in the solve JSON and the GUI's **Integral** tab (see
[GUI §6.11](gui.md#611-integral-table-tab--integraltabletsx)).

### 12.4 `INTEGRALVALUE` — retrieving a tabulated value

```ees
y_prev = INTEGRALVALUE(t-0.5, 'y')   "interpolate the trajectory of y"
```

`INTEGRALVALUE(t, 'X')` returns the value of variable `X` at integration
variable value `t` by **linear interpolation** of the trajectory built so far
(with flat clamping at the endpoints). It is meaningful only *during* an
integration step. The parser recognises the function; full evaluator dispatch is
a deferred follow-up (see [Dynamic Solving §7.2](integral_table.md#72-planned-improvements)).

### 12.5 Integration methods and configuration

Integration is controlled by `coolsolve.conf` keys (all inert by default —
uncomment only what you need; see `examples/coolsolve.conf`):

| Key | Default | Meaning |
|:----|:--------|:--------|
| `integralMethod` | `RK4` | `RK4`, `RK45`, `EulerExplicit`, or `EulerImplicit` |
| `integralFixedStep` | `0.0` | Fixed step size (used when the `INTEGRAL` call gives none). `0` ⇒ `integralMaxSteps` steps (fixed methods) or adapt (RK45) |
| `integralMaxSteps` | `1000` | Number of steps of a fixed-step method that has no explicit step; budget of RK45 (§12.5) |
| `integralRelTol` | `1e-6` | RK45 local relative error control |
| `integralAbsTol` | `1e-9` | RK45 absolute error floor |
| `integralMinStep` | `0.0` | Minimum step (0 = auto) |
| `integralMaxStep` | `0.0` | Maximum step (0 = auto) |
| `integralRichardson` | `false` | Richardson extrapolation (fixed-step methods only) |
| `integralOutputInterval` | `0.0` | Default row interval when `$IntegralTable` omits `:n`. `0` = every step |

**Step size and number of steps.** A step given in the model — the 5th argument of
`INTEGRAL`, as in EES — is used as written and the number of steps follows from
`(hi - lo)/step` (the last step is shortened when the step does not divide the
interval); `integralMaxSteps` does not apply. Without it, `integralFixedStep` is used
the same way. With neither, a fixed-step method takes exactly `integralMaxSteps`
steps (`(hi - lo)/integralMaxSteps` each), and RK45 chooses its own steps within a
budget of `4·integralMaxSteps` attempted steps (accepted plus rejected; the explicit
step of the model is then its first step). A run that stops before `hi` is an
**error** (*"The integration stopped at … before reaching the upper limit …"*), never
a truncated success.

Method guidance:

- **`RK4`** (default) — classic 4th-order Runge–Kutta, fixed step. Good accuracy
  for smooth, non-stiff systems.
- **`RK45`** — adaptive Dormand–Prince (embedded 4th/5th order). Adjusts the step
  from the error estimate; the right choice when the dynamics have fast and slow
  phases. Watch the rejected-step count on stiff systems.
- **`EulerExplicit`** — 1st-order, cheapest; useful for quick checks.
- **`EulerImplicit`** — 1st-order, A-stable; the quick choice for stiff systems.

The `$IntegralAutoStep` and `$IntegralStop` directives are recognised for EES
compatibility but **not interpreted** — they emit a diagnostic pointing at the
`integral*` config keys instead.

### 12.6 Worked example

`integral_decay.eescode` — exponential decay `dy/dt = -y`, `y(0) = 1`, whose
analytical solution is `y(t) = e^{-t}`:

```ees
"Exponential decay:  dy/dt = -y,  y(0) = 1  =>  y(t) = e^{-t}"

y = 1 + integral(dydt, t, 0, 4)
dydt = -y

$IntegralTable t:0.5  y  dydt
```

Solving gives `y(4) = 0.0183156 = e^{-4}`, recorded in
`integral_decay-integral.csv`:

```
t,y,dydt
0,1,-1
0.5,0.606531,-0.606531
…
4,0.0183156,-0.0183156
```

### 12.7 Limits, and limitations

**Limits and step.** `lo`, `hi` and `step` may be constants, variables or expressions of
the model (`t_1 = 0`, `t_2 = 18000`, `v1_bar = V1/N`, `v[2]` computed by a property
call, `p_out`, ...). CoolSolve solves the model once at the initial point (`t = lo`,
every integral zero) and reads them there; constant limits are folded without that
solve and keep their cost. A limit may not depend on the integration variable or on a
state — it must stay the same along the integration: such a model is rejected with
*"The upper limit '…' of INTEGRAL() depends on the integration variable … or on a state
of the model"*. All calls must resolve to the same interval.

- **One integration variable and one interval per model.** All `INTEGRAL` calls
  must share the same `t` and `[lo, hi]`. Nested integrals and calls with different
  integration variables are rejected with a clear message.
- **Increasing integration only.** `hi` must be greater than `lo`.
- **`INTEGRAL` inside a `DUPLICATE` loop** is rejected (write the equations one by one).
- **Semi-explicit index-1 DAE only.** Higher-index systems (where an algebraic
  constraint couples to a derivative that must be differentiated to solve) are
  detected and reported; index reduction is not yet implemented. A *warning* that
  a state variable appears in an algebraic equation is benign for ordinary
  index-1 thermal models.
- **2-arg table-based `INTEGRAL(integrand, var)`** (integration over a parametric
  table) is not supported — it requires an EES-style Parametric table CoolSolve
  does not have.
- **Stiff systems.** An explicit method (RK4/RK45/EulerExplicit) may need many
  steps on a stiff system; use `EulerImplicit` as a workaround until a BDF/stiff
  integrator is added (see `docs/solver_roadmap.md` §9.4).
