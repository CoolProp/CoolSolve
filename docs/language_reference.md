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
| `"..."` (dquote)  | Block or inline comment (may span multiple tokens)   |
| `{...}`           | Brace comment (inline)                               |
| `// ...`          | Single-line comment                                  |

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

### Directives

Parser-level directives begin with `$` and span a single line:

- `$if <cond>` / `$ifnot <cond>` / `$endif`
- `$unitSystem SI` — recognised for EES compatibility but currently **ignored**;
  the unit system is fixed (see the warning at the top of this page).

---

## 2. Built-in mathematical functions

### Elementary functions

| Function                 | Meaning                               |
|:-------------------------|:--------------------------------------|
| `sin`, `cos`, `tan`      | Trigonometric (radians)               |
| `sinh`, `cosh`, `tanh`   | Hyperbolic                            |
| `arcsin` / `asin`        | Inverse sine                          |
| `arccos` / `acos`        | Inverse cosine                        |
| `arctan` / `atan`        | Inverse tangent                       |
| `arcsinh` / `asinh`      | Inverse hyperbolic sine               |
| `arccosh` / `acosh`      | Inverse hyperbolic cosine             |
| `arctanh` / `atanh`      | Inverse hyperbolic tangent            |
| `atan2(y, x)`            | Two-argument arctangent               |
| `exp(x)`, `ln(x)`        | Natural exponential / logarithm       |
| `log(x)` / `log10(x)`    | Base-10 logarithm                     |
| `sqrt(x)`                | Square root                           |
| `abs(x)`                 | Absolute value                        |
| `pi` / `pi()`            | Mathematical constant π               |

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

`if(cond, val_if_true, val_if_false)` — `cond > 0` selects the true branch.

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

---

## 3. Thermophysical property functions (CoolProp-backed)

All thermophysical functions take a fluid as the first positional argument
and state specifications as **named arguments** (typically two state
variables). Example:

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
| `molarmass` / `MM`                   | kg/mol (no state needed)   |

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

### Procedures

```ees
PROCEDURE divmod(x, y : quotient, remainder)
    quotient  = FLOOR(x / y)
    remainder = x - quotient * y
END

CALL divmod(17, 5 : q, r)
```

The `:` separates input arguments (left) from output variables (right).

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
message.

Each CSV must have a header row with column names (case-insensitive).

### 11.1 1D interpolation

```ees
h = INTERPOLATE('data', 'T_C', 'h_kJ_per_kg', T)
```

Performs linear interpolation of column `h_kJ_per_kg` as a function of
column `T_C` at value `T`.  Values outside the range are clamped to the
nearest endpoint (flat extrapolation, zero derivative).

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

### 12.1 `INTEGRAL` — declaring a state variable

```ees
y = <base> + INTEGRAL(integrand, t, t0, tf)
y = <base> + INTEGRAL(integrand, t, t0, tf, step)   "fixed-step form"
```

| Argument | Meaning |
|:---------|:--------|
| `integrand` | The derivative expression (the `f` in `dy/dt = f`). Usually a named variable defined by its own equation. |
| `t` | The integration variable (e.g. time). All `INTEGRAL` calls in a model must share the **same** `t`. |
| `t0`, `tf` | Lower and upper limits. Must be constants (or constant-foldable). All calls must share the **same** `[t0, tf]`. |
| `step` | Optional 5th argument: a **fixed** step size. Omit it to let the method choose (fixed default derived from `integralMaxSteps`, or adaptive for RK45). |

`<base>` is the non-integral part of the right-hand side; at `t = t0` the
integral term is zero, so `y(t0)` is determined by the algebraic solve of the
remaining equations (no separate initial-value extraction is needed).

The variable on the left-hand side of an integral equation is a **state
variable** (owned by the integrator). Every other unknown is *algebraic* — it is
solved from the remaining equations at each time step by the same algebraic
solver used for static models.

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
first. An optional output interval follows the integration variable as
`name:interval` (`t:0.1` → one row every 0.1 units); if omitted, a row is
written at every step (or every `integralOutputInterval`, see §12.5). Array
ranges expand with the `..` notation:

```ees
$IntegralTable t  X[1..5]    "records t and X[1]..X[5]"
```

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
| `integralFixedStep` | `0.0` | Fixed step size. `0` ⇒ derive from `integralMaxSteps` (fixed methods) or adapt (RK45) |
| `integralMaxSteps` | `1000` | Upper bound on the number of steps |
| `integralRelTol` | `1e-6` | RK45 local relative error control |
| `integralAbsTol` | `1e-9` | RK45 absolute error floor |
| `integralMinStep` | `0.0` | Minimum step (0 = auto) |
| `integralMaxStep` | `0.0` | Maximum step (0 = auto) |
| `integralRichardson` | `false` | Richardson extrapolation (fixed-step methods only) |
| `integralOutputInterval` | `0.0` | Default row interval when `$IntegralTable` omits `:n`. `0` = every step |

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

### 12.7 Limitations

- **One integration variable and one interval per model.** All `INTEGRAL` calls
  must share the same `t` and `[t0, tf]`. Nested (multi-variable) integration is
  detected and rejected with a clear message.
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


