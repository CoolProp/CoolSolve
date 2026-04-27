# CoolSolve Language Reference

[← Back to overview](../README.md)

This page describes the CoolSolve language — an EES-compatible subset with
thermodynamic property functions backed by CoolProp. Only features that are
currently implemented and covered by tests are documented here.

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
T_1 = 300  "[K]"
```

Units annotations are informational; the solver uses the unit system
configured in `coolsolve.conf`.

### Directives

Parser-level directives begin with `$` and span a single line:

- `$if <cond>` / `$ifnot <cond>` / `$endif`
- `$unitSystem SI` / `$unitSystem K` etc.

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
rho = density('R134a', T=280, P=500e3)
```

Inputs and outputs are in the unit system configured in `coolsolve.conf`
(default temperature is Celsius).

### Standard properties

| Function                             | Result (SI)                |
|:-------------------------------------|:---------------------------|
| `temperature` / `T`                  | K                          |
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
(`T` [K|°C], `P` [Pa], `R` = relative humidity in [0..1]). Outputs are
returned in the configured unit system.

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
rho = density('MEG', 30, T=280, P=101325)   "30 % MEG"
cp  = cp('EthyleneGlycol', T=280, P=101325, C=30)
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

- `[units]` — temperature, pressure, energy, …
- `[coolprop]` — backend (`HEOS`, `BICUBIC&HEOS`, `TTSE&HEOS`), analytical
  derivatives, AbstractState caching.
- `[solver]` — pipeline order (Newton, TrustRegion, LM, BisectionND,
  Homotopy, Partitioned), tolerances, iteration limits.
- `[tearing]` — SCC-based block reduction options.

---

## 11. Lookup tables

Lookup tables let you interpolate from external CSV data files directly inside
equations.  CoolSolve loads companion tables for a model automatically by
scanning the same directory for CSV files whose name matches the model stem:

1. **`<modelname>.csv`** — main companion table; table name = `<modelname>`.
2. **`<modelname>_<suffix>.csv`** — auxiliary tables; table name =
   `<modelname>_<suffix>` (the full file stem).  Use this convention to
   bundle multiple lookup tables alongside a model, e.g.
   `steam_cycle_watercp.csv`, `steam_cycle_refprops.csv`.

No other CSV files in the directory are loaded automatically.

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
h = INTERPOLATE('mymodel', 'T_C', 'h_kJ_per_kg', T)
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
P_at_row3 = LOOKUP('mymodel', 3, 2)
```

Returns the value in the given row and column (both 1-based).  No
interpolation — direct cell lookup.

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

For a model file `steam_cycle.eescode` with companion `steam_cycle.csv`:

```
T_C,P_kPa,h_kJ_per_kg
100,101.3,2675.6
150,476.1,2745.9
200,1554.9,2826.8
```

```ees
T = 150  "°C"
h = INTERPOLATE('steam_cycle', 'T_C', 'h_kJ_per_kg', T)  "h = 2745.9 kJ/kg"
n = NLOOKUPROWS('steam_cycle')                             "n = 3"
P_row2 = LOOKUP('steam_cycle', 2, 2)                      "P_row2 = 476.1 kPa"
```

The **Lookup Tables** tab in the GUI lets you create, view, and edit tables
directly in the browser — no need to manage CSV files manually.

