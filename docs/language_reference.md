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
5. `<`, `<=`, `>`, `>=` (return 1/-1 in condition contexts)

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
and `A[1][2]` are not annotations, see below). A plain variable does not carry
an annotation: `x [m]` is the array element `x[m]`.

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
- `$unitSystem SI` — recognised for EES compatibility but currently **ignored**;
  the unit system is fixed (see the warning at the top of this page).
- `$warnings ON` / `OFF` — recognised and **ignored**; the warnings of `CALL WARNING` are
  always reported, as in EES (see *CALL ERROR and CALL WARNING* in §7).
- `$CONSTANT Name# = value [units]` — a constant of the model, see below.

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
version number. Tolerated EES spellings: a trailing `then` after the condition
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
| `sum2d(a, b, …)`         | Sum of squares                  |
| `average(a, b, …)`       | Arithmetic mean                 |
| `product(a, b, …)`       | Product of arguments            |
| `stddev(a, b, …)`        | Population standard deviation   |

`min(a, b, …)` and `max(a, b, …)` take any number of arguments too. An array range
is a list of arguments: `average(T[1..imax#])`, `sum(x[1..n], 10, y)` (§8). The
indexed form `sum(x[i], i=1, n)` is not supported.

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
kept or shown: a string output, a message. In the main program a string variable is only set
by a literal or by a `CALL` output, as for `LOOKUP$` (§11.8): `s$ = Uppercase$('ab')` at the
top level is an error (*Uppercase$ returns a string and cannot be used as a number*), not a
silent 0; compute the string in a procedure, or use the call where the string is needed.

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
| `CONVERTTEMP('from', 'to', T)`     | Convert a temperature value between scales  |
| `UNITSYSTEM('setting')`            | 1 if the unit setting is in use, 0 otherwise (`'SI'`, `'Eng'`, `'Mass'`, `'Molar'`, `'Deg'`, `'Rad'`, `'C'`, `'K'`, `'F'`, `'R'`, `'Pa'`, `'kPa'`, `'bar'`, `'MPa'`, `'psia'`, `'atm'`, `'J'`, `'kJ'`, `'Btu'`, `'kcal'`) — CoolSolve always answers 1 for `SI`, `Mass`, `Deg`, `C`, `Pa` and `J` |

`UNITSYSTEM` lets EES library routines check the unit settings they were
written for; in CoolSolve such a guard never fires (the unit system is fixed).
An unknown setting name is an error.

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
| `molarmass` / `MM`                   | kg/kmol, as in EES (no state needed) |

### Saturation properties

| Function                             | Inputs                     |
|:-------------------------------------|:---------------------------|
| `t_sat` / `tsat`                     | `P=...` (returns saturation temperature) |
| `p_sat` / `psat`                     | `T=...` (returns saturation pressure)    |

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
| `enthalpy`              | Specific enthalpy of the mixture       |
| `entropy`               | Specific entropy                       |
| `volume`                | Specific volume m³/kg dry-air          |

Example:

```ees
W_1   = humrat  ( 'airH2O', T=25, P=101325, R=0.5 )
Twb_1 = wetbulb ( 'airH2O', T=25, P=101325, R=0.5 )
```

### PSYCHPROPS built-in procedure

```ees
CALL psychprops(T_1, P_1, R_1 : T, v, h, s, u, W, R, Twb, Tdp)
```

Returns all nine psychrometric properties at the specified state
(`T` in °C, `P` in Pa, `R` = relative humidity in [0..1]). Outputs are
returned in CoolSolve's fixed unit system (temperature in °C, see the warning
at the top of this page).

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
`H2O` (as vapor, low pressure only).

The properties that depend on temperature only (`enthalpy`, `intenergy`,
`cp`, `cv`) can be called with `T` alone, as in EES — typically for mean
specific heats:

```ees
c_bar_p = (enthalpy(CO2,T=T_p1) - enthalpy(CO2,T=T_p2)) / (T_p1 - T_p2)
```

CoolSolve evaluates them at an internal reference pressure (101325 Pa; 100 Pa
for `H2O`, which keeps the vapor phase down to −20 °C). This pressure is part
of the model: it is never clamped by the solver and is identical in the solve
and in the final verification.

**Reference states of the enthalpy and of the entropy.** As in EES (JANAF convention), the
ideal gases `N2`, `O2`, `H2`, `He`, `Ar`, `CO2`, `CO`, `CH4`, `C2H6`, `C3H8` and `H2O`
(vapor) have absolute properties from **one source**, the NASA-Glenn thermodynamic
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
| `CO2` | −393 510 | 213.786 | | | |

The enthalpies of formation are the values of the species headers of `thermo.inp`; `s°` is
the 9-coefficient polynomial of the same file at 298.15 K (R = 8.314462618 J/mol-K). EES uses the
same data with R = 8.314, which puts its values 6e-5 lower (EES help: `CO2` −393 486,
`H2O` −241 811, `CH4` −74 595 kJ/kmol); NIST-JANAF entropies differ from NASA's by at most
0.12 kJ/kmol-K (`CH4`). `Air` is not in the table: it keeps CoolProp's reference states (see below).
Earlier CoolSolve versions used the rounded JANAF values of Çengel & Boles (Table A-26),
5 to 830 kJ/kmol away (`C2H6` and `C3H8` the most), and `He` and `Ar` kept CoolProp's
enthalpy reference.

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

Pure: `Aluminum`, `Copper`, `LiBr`, `SeaWater`.

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

Example:

```ees
rho = density('MEG', 30, T=10, P=101325)   "30 % MEG at 10 °C"
cp  = cp('EthyleneGlycol', T=10, P=101325, C=30)
```

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

### Unsupported EES fluids

Some EES fluids are not modelled by CoolProp. Attempting to use them
triggers a clear runtime error with the fluid name and a remediation
suggestion:

`NH3H2O`, `LiBrH2O`, `Acetylene`, `Chlorine`, `Chloroethene`,
`EthyleneOxide`, `Helium3`, `Potassium`, `Sodium`.

For the ammonia-water and lithium-bromide-water systems, use external
correlations written as `FUNCTION`s (the REFPROP library is not available in
CoolSolve).

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
  `>=` (strings are compared without regard to case, as in EES:
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
body may be followed by a comment (`END // area`). `GOTO` and statement labels
(`10:`) are not supported and are refused the same way.

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

  (A lookup-table cell write `LOOKUP('t', 1, 1) = x` is not concerned: writing a
  cell is not supported yet and the statement has no effect.)
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
ascending **or descending** order (EES accepts both).  Values outside the range
are clamped to the nearest endpoint (flat extrapolation, zero derivative).

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

> CoolSolve interpolates **linearly** and holds the end values outside the table;
> EES `INTERPOLATE` is cubic and extrapolates, so the results can differ slightly
> inside the table and widely outside it.

Analytical derivative (for Newton/AD): slope of the active interval;
zero when clamped at an endpoint.

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
needs 3 arguments, `INTERPOLATE` 4 (or 3 and a named argument, §11.1), `INTERPOLATE2` 6 (or 4 and two
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
  main program a string variable is only set by a literal or by a `CALL` output, so
  `s$ = LOOKUP$(...)` at the top level is an error (*LOOKUP$ returns a string and cannot be
  used as a number*): use the call where the string is needed, or compute the string in a
  procedure. Used as a number it is an error too.
- Numeric functions work as before on a table with text cells: `LOOKUP('Species', 2, 2)` is
  31.9988, `LOOKUP('Species', 2.5, 'M')` interpolates between the rows, and a text or empty
  cell (read as NaN) spoils the result only where it is really read. `LOOKUP$ROW` is not
  implemented.
- The GUI Lookup Tables panel shows and edits the text cells like the other cells (the table
  is kept as CSV text); no GUI setting is needed.

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
