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
- Identifiers: start with letter or `_` / `$`; may contain `[a-zA-Z0-9_$#]`.
- Identifiers ending with `$` are **string variables** (e.g. `fluid$`).

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
*Outputs of a CALL* in §7). A call without argument is written `f()`. Array
ranges (`x[1..3]`) are not expression syntax, in `CALL` arguments neither.

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

### Aggregation functions

| Function                 | Meaning                         |
|:-------------------------|:--------------------------------|
| `sum(a, b, …)`           | Sum of arguments                |
| `sum2d(a, b, …)`         | Sum of squares                  |
| `average(a, b, …)`       | Arithmetic mean                 |
| `product(a, b, …)`       | Product of arguments            |
| `stddev(a, b, …)`        | Population standard deviation   |

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
property function: `density('Water', T=25, P=101325)`.

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

**Enthalpy reference state.** As in EES (JANAF convention), the enthalpy and
internal energy of the species below include the standard enthalpy of
formation: `h(25 °C)` equals the heat of formation (elements: 0), so heats of
reaction can be written with species enthalpies:

| Species | `h·M` at 25 °C [kJ/kmol] | Species | `h·M` at 25 °C [kJ/kmol] |
|:--|--:|:--|--:|
| `N2`, `O2`, `H2` | 0 | `CH4` | −74 850 |
| `CO2` | −393 520 | `C2H6` | −84 680 |
| `CO` | −110 530 | `C3H8` | −103 850 |
| `H2O` (vapor) | −241 820 | | |

```ees
LHV_CO = enthalpy(CO,T=25)*molarmass(CO) + 0.5*enthalpy(O2,T=25)*molarmass(O2) - enthalpy(CO2,T=25)*molarmass(CO2)
"= 282 990 kJ/kmol (in J/kmol: 2.8299e8)"
```

The reference applies to all calls, including `h=` given as an input
(`temperature(CO2, P=P, h=h)`); enthalpy differences are unchanged. Entropies,
and the enthalpies of `Air`, `He`, `Ar` and of the real fluids (`Water`,
`CarbonDioxide`…), keep CoolProp's reference states. Values: Çengel & Boles,
Table A-26, which reproduce the heats of combustion of EES.

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

### Unsupported EES fluids

Some EES fluids are not modelled by CoolProp. Attempting to use them
triggers a clear runtime error with the fluid name and a remediation
suggestion:

`NH3H2O`, `LiBrH2O`, `Acetylene`, `Chlorine`, `Chloroethene`,
`EthyleneOxide`, `Helium3`, `Potassium`, `Sodium`.

Use `REFPROP::Fluid1&Fluid2` backends or external correlations for the
ammonia-water / lithium-bromide-water systems.

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

### CALL ERROR

`CALL ERROR('message')` stops the calculation and reports the message — the
EES way to reject inputs in a function or procedure, usually in an `IF`:

```ees
IF (T_p1 = T_p2) THEN CALL ERROR('T_p1 and T_p2 must differ', T_p1)
```

Additional numeric arguments are appended to the message in brackets. (`CALL
WARNING` is not available yet.)

---

## 8. Arrays

Arrays use `[index1, index2, ...]` indexing and are resolved at solve
time:

```ees
DUPLICATE i = 1, 5
    T[i] = T_in + i * DeltaT
END
```

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
needs 3 arguments, `INTERPOLATE` 4, `INTERPOLATE2` 6.

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
