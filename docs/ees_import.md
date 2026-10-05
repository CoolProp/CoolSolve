# Importing EES Models into CoolSolve

This guide explains how to read binary EES files (`.ees`), extract everything
they contain (equations, guess values, stored solution, lookup and parametric
tables, unit settings), turn them into a CoolSolve model folder, solve them and
**verify the result against the solution stored by EES**. It is the reference
procedure used to populate the
[CoolSolve Library](https://github.com/IntSusEnergySystems/CoolSolve_Library).

The procedure was tested on 2026-10-04 on three EES files (see
[§12 Worked examples](#12-worked-examples)) and the extraction tool was run on
the 457 EES files (EES 6.x to 10.x) of the ULiège teaching/research collection.

| Tool | Purpose |
|------|---------|
| [`tools/ees_extract.py`](../tools/ees_extract.py) | Decode a binary `.ees` file into a CoolSolve folder (`.eescode`, `.initials`, lookup `.csv`, reference solution, report) |
| [`tools/compare_solution.py`](../tools/compare_solution.py) | Compare a CoolSolve `.sol` file with the EES reference solution or with one run of an EES parametric table |

Both scripts are plain Python 3 (no dependency).

---

## 1. Principles

1. **Keep the native EES language.** CoolSolve aims to read native EES code.
   When an EES construct is not supported, or behaves differently in
   CoolSolve, **do not rewrite the model around it**: keep the EES syntax and
   report the gap in [model_library_support.md](model_library_support.md)
   (§4, *gap register*). The model is then marked *blocked* by that gap.
2. **Edits that are allowed** (and must be logged in the model's README):
   - translating comments to English and adding the standard header;
   - converting the model **by hand** to SI / °C / Pa / J (§6);
   - fixing genuine errors of the original model;
   - giving default values to inputs that EES took from a parametric table (§8);
   - copying library functions that EES loaded implicitly (§9), until
     CoolSolve can import them;
   - simplifications needed to make a model converge (documented, and kept as a
     separate variant when they change the physics).
3. **Never trust a conversion that was not verified.** EES stores the values of
   the last calculation in the file: every import ends with a numerical
   comparison (§11).
4. **Reproducibility without copies.** Use the scripts. The extraction is done
   in a temporary work folder; the library keeps only the CoolSolve files of
   the model and references the source file by its path (with `~` for the home
   directory), its URL, or its name in an existing library. The work folder
   is deleted once the model is included.

---

## 2. What an EES file contains

| Content | Where EES keeps it | Extracted by `ees_extract.py` | CoolSolve counterpart |
|---------|-------------------|:---:|---|
| Equations window | Binary file, plain Windows-1252 text or RTF ("formatted equations") | ✅ all versions | `<model>.eescode` |
| Unit system (SI/Eng, mass/molar, °C/K, kPa/bar/Pa, deg/rad, J/kJ) | Bytes following the equations | ✅ EES 7–10 (older files: partially) | Fixed: SI, mass, °C, Pa, J, degrees (§6) |
| Variable information: stored value, guess, lower/upper bounds, units | Fixed-size records | ✅ EES 7–10 (EES 6: best effort) | `<model>.initials` (guesses only; bounds unsupported) |
| Embedded Lookup tables | Grid of columns (name, units, values) | ✅ | `<model>-<table>.csv` |
| Parametric tables | Same grid structure, string columns allowed | ✅ | No file equivalent (GUI parametric studies) — kept as reference results |
| External lookup files `.lkt`, `.txt`, `.csv` | Separate files next to the `.ees` | `.lkt`: ❌ (not decoded yet) | `<model>-<table>.csv` |
| Function libraries `.lib` | Separate plain-text files (EES `USERLIB` folder) | n/a (text) | Function model (§9) |
| Compiled external functions `.dlf`, `.dlp`, `.fdl` | DLLs | ❌ | None (re-implement as EES `FUNCTION` if the algorithm is known) |
| Min/Max (optimisation) settings, uncertainty propagation | GUI settings | ❌ | Not available (optimisation is a CoolSolve gap) |
| Diagram window, plots, formatted text | GUI objects | ❌ | Re-create figures for the README if useful |
| Integral table settings | `$IntegralTable` directive in the text | ✅ (text) | `$IntegralTable` (space-separated columns) |

Things EES does **behind the scenes** and that are therefore *not* in the
equations: functions from the user library (`USERLIB` `.lib` files loaded at
start-up, e.g. `cpbar`), values given by a parametric table, and the unit
system. These are the three most common reasons why an extracted model is not
complete.

---

## 3. Quick start

```bash
# 1. Extract (creates the folder, never modifies the .ees file)
python3 tools/ees_extract.py "path/to/Model.EES" -o work/my_model --name my_model

# 2. Read the report: version, unit system, tables, functions called, warnings
cat work/my_model/reference/ees_extract_report.md

# 3. Complete / convert the model (§5–§9), then solve (use ./ or an absolute path, see §7)
cd work/my_model && /path/to/coolsolve -d ./my_model.eescode

# 4. Verify against the solution stored by EES (or a parametric-table run)
python3 /path/to/CoolSolve/tools/compare_solution.py my_model.sol reference/ees_variables.csv
python3 /path/to/CoolSolve/tools/compare_solution.py my_model.sol reference/ees_parametric-table_1.csv --row 1
```

`ees_extract.py --info file.ees` prints a one-screen summary without writing
anything (useful to triage many files). The output folder is a **temporary
work folder** (e.g. `work/<name>/` in the library, ignored by git); only the
model files are moved to the library, the rest is deleted after the import:

```
my_model/
├── my_model.eescode              equations (UTF-8), with a $UnitSystem line reflecting the EES settings
├── my_model.initials             guess values = values stored by EES (last solution)
├── my_model-<table>.csv          one file per embedded Lookup table (CoolSolve naming)
└── reference/
    ├── ees_extract_report.md     provenance, unit system, features, functions called, warnings
    ├── ees_variables.csv         name, value, guess, lower, upper, units  (EES stored solution)
    └── ees_parametric-<t>.csv    parametric tables (reference results or inputs)
```

---

## 4. Text encoding and EES-only tags

- EES writes the equations in **Windows-1252** (`é`, `à`, `°`, `²`, `µ`…) with
  CRLF line endings. The tool converts to UTF-8 / LF. Never open an `.ees` file
  as UTF-8 text: the accents are lost (`d�bit`).
- About 5 % of the files store **RTF** (`{\rtf1\ansi…`): the tool converts it to
  plain text (`\par` → newline, `\'e9` → `é`, colour/font groups dropped).
  Formatting (bold, colours) is lost, the text is not.
- Tags such as `{$ID$ #1206: Jean Lebrun, Laboratoire de Thermodynamique…}`
  (EES licence registration), `{$PX$96}` or `{$ST$OFF}` (display settings)
  are removed; they are listed in the report. The `$ID$` tag identifies the
  **EES licence** used to save the file, *not necessarily the author* of the
  model.
- **Decimal comma.** With the European number format, EES stores the text as
  typed: `A_t = 0,5` and `volume(methane;P=p_t;x=0)` (`;` separates
  arguments). CoolSolve cannot parse it (`CS-GAP-DECIMAL-COMMA`). The tool
  detects this format (a `;` inside parentheses, or a decimal comma after `=`)
  and converts the code — not the comments — to the dot convention; the report
  says so. This is a display convention of EES, not a change of the model.
  It concerns 183 of the 457 files of the ULiège collection (40 %); a random
  sample of 40 converted files all parse in CoolSolve.
- Variable names are case-insensitive in both EES and CoolSolve (`t_1` and
  `T_1` are the same variable): no renaming is needed.
- Comments (`"..."`, `{...}`) in French or Dutch are translated to English in
  the library version; the path of the source file is recorded in the model's
  metadata.

---

## 5. Variables, guesses and the stored solution

- `ees_variables.csv` holds, for each variable of the main program, the
  **value stored by EES** (the last solution), the guess value, the bounds and
  the units entered in the Variable Information dialog.
- `<model>.initials` is written from the stored values (falling back on
  non-default guesses). These are excellent initial guesses for CoolSolve, as
  long as the unit system is unchanged (§6).
- EES writes **−9999** for values cleared before a new calculation. When most
  values are −9999, the file was last run from a **parametric table** (or never
  solved): the tool warns, and the parametric table is the reference (§8).
- **Bounds** (lower/upper limits) are not supported by CoolSolve's `.initials`
  file. When EES relied on them to stay in a physical region (e.g. a quality
  between 0 and 1), note it in the model README; this is gap `CS-GAP-BOUNDS`.
- Variables of FUNCTIONs and PROCEDUREs are local and have no stored value.

---

## 6. Unit system — manual conversion

CoolSolve's unit system is fixed: **SI, mass basis, °C, Pa, J, degrees for
trigonometric functions** — the EES setting `$UnitSystem SI MASS DEG PA C J`.
The extractor decodes the unit settings stored in the EES file and writes them
as a `$UnitSystem` line at the top of the extracted `.eescode`, so that the
units of the original are explicit. In the ULiège collection, 288 of the 540
EES files already use SI-C-Pa-J; the others use kPa/kJ, K, bar or molar units,
and some the RAD trigonometric setting.

> **Rule: when the EES model uses another unit system, the whole model is
> converted by hand, equation by equation, by the person doing the import.**
> Do not rely on the `$UnitSystem` directive (CoolSolve parses but ignores it:
> a kPa/kJ model would silently give wrong results), do not scale results
> automatically, and do not wrap property calls in conversion factors: the
> converted model must read naturally in SI-°C-Pa-J, like any other library
> model. Every converted value and every modified equation is checked by hand
> and listed in the README conversion log.

### 6.1 Prepare a conversion table

1. Read the unit system in `reference/ees_extract_report.md` (work folder).
2. Read the units EES stored for each variable (`reference/ees_variables.csv`,
   column `units`): they tell which variables are pressures, energies,
   absolute temperatures… Variables without stored units must be identified
   from the equations and comments.
3. Write down, for the model, the conversion of every quantity whose unit
   changes:

| Quantity | Typical EES unit | CoolSolve unit | Conversion of a value |
|---|---|---|---|
| Pressure | kPa · bar · MPa · atm | Pa | ×1 000 · ×100 000 · ×10⁶ · ×101 325 |
| Specific energy (h, u, w, q, Δh) | kJ/kg | J/kg | ×1 000 |
| Specific entropy, heat capacity, gas constant | kJ/kg-K | J/kg-K | ×1 000 |
| Power, heat rate | kW · MW | W | ×1 000 · ×10⁶ |
| Energy | kJ · MJ · kWh | J | ×1 000 · ×10⁶ · ×3.6·10⁶ |
| Conductance UA, heat-transfer coefficient | kW/K · kW/m²-K | W/K · W/m²-K | ×1 000 |
| Absolute temperature | K | °C | −273.15 |
| Temperature difference | K | K (= °C difference) | unchanged |
| Molar quantities | kJ/kmol, kmol/s | J/kg, kg/s | through the molar mass |
| Angle in trigonometric functions | rad (RAD setting) | degree | argument ×180/π |

### 6.2 Convert the model, equation by equation

Go through the whole Equations window from top to bottom, including the
FUNCTIONs and PROCEDUREs:

1. **Input values.** Rewrite each value in the new unit and keep the original
   value in the comment, e.g.
   `p_t = 500E3 [Pa]  "500 kPa in the original"`,
   `T_r = 1.85 [C]  "275 K in the original"`, `cp = 1005 "[J/kg-K] (1.005 kJ/kg-K)"`.
2. **Property calls.** Arguments and results of `enthalpy`, `pressure`,
   `entropy`… are now in Pa, J and °C: remove any factor the original applied
   around them (`h = enthalpy(...)*1000`, `P = p_bar*100`, `T = temperature(...)-273.15`).
3. **Balances and other homogeneous equations** usually stay unchanged when all
   their terms change consistently (kW and kJ/kg → W and J/kg). Check every
   equation that mixes quantities converted with different factors, or that
   contains an explicit unit factor (`/1000`, `*1000`, `/3600`, `*100`).
4. **Hidden unit-dependent constants**: gas constants (`R = 0.287` → `287`),
   `101.325` (atmospheric pressure in kPa), `9.81` combined with kJ,
   `4.18`/`4.186` (water cp in kJ/kg-K), `2501` (latent heat in kJ/kg). The
   numeric coefficients of **empirical correlations** were fitted in the
   original units: keep the correlation as it is and feed it with a local
   variable in the original unit (`P_kPa = P/1000`), with a comment.
5. **Absolute temperatures.** In a model written in K, every relation that
   needs an absolute temperature must use `(T + 273.15)` once `T` is in °C:
   ideal-gas law, isentropic relations `T2/T1 = (P2/P1)^((k-1)/k)`, Carnot
   factors `1 - T0/T`, radiation `σ·T⁴`, exergy terms `(T0 + 273.15)*(s - s0)`,
   `ln(T2/T1)`, Arrhenius laws. Alternatively keep explicit absolute-temperature
   variables (`T_K = T + 273.15`) for these relations. Temperature
   *differences* are unchanged.
6. **Molar basis** → mass basis (divide molar properties by the molar mass,
   `MOLARMASS(fluid$)`), **RAD** → degrees in trigonometric arguments
   (`sin(x)` → `sin(x*180/pi)`) or angles stored in degrees.
7. **Units annotations and comments**: update `[kPa]`, `"[kJ/kg]"`… to the new
   units.
8. **Tables**: convert the columns of the lookup tables (`<model>-<table>.csv`)
   and, for the verification, the parametric-table values that changed unit.
9. **Guess values**: convert `<model>.initials` (pressures ×1000, enthalpies
   ×1000, temperatures −273.15), or delete the converted entries and obtain new
   guesses with a simplified model
   ([Debugging Models §5](debugging_models.md#step-5-break-the-loop-simplified-model)).
10. Replace the `$UnitSystem` line by `$UnitSystem SI MASS DEG PA C J`.
11. Log the conversion in the README: the converted inputs, the modified
    equations and the constants found.

`CONVERT('kPa','Pa')` and `CONVERTTEMP('K','C',T)` exist (with quoted unit
names — `ConvertTemp(K, C, T)` without quotes, valid in EES, fails in CoolSolve:
`CS-GAP-CONVERT-UNQUOTED`). In the library, prefer converted numbers with the
original value in a comment; keep `CONVERT` only where a model deliberately
reports a result in another unit (e.g. an energy in kWh).

### 6.3 Verify the conversion

Solve, then compare with the solution stored by EES, converting the reference
values with the units EES stored:

```bash
python3 tools/compare_solution.py model.sol reference/ees_variables.csv --ees-units
```

`--ees-units` converts kPa/bar/MPa → Pa, kJ → J, kW → W and K → °C (for
variables whose name starts with `T`). Every variable must agree within the
tolerances of §11. A remaining deviation by a factor 1 000 (or 100 000), or by
273.15, points to a missed conversion; a deviation on a variable EES stored
without units must be checked by hand.

| Symptom after conversion | Likely cause |
|---|---|
| CoolProp error: pressure or enthalpy out of range | a pressure left in kPa / an enthalpy left in kJ/kg |
| Results off by a factor 1 000 | constant or correlation coefficient in kJ or kPa (step 4) |
| Efficiency or ideal-gas results wrong, balances right | absolute temperature missing in a relation (step 5) |
| Trigonometric results wrong | RAD setting (step 6) |
| Convergence failure from the start | `.initials` not converted (step 9) |

Worked example: §12, Example 3.

## 7. Lookup tables

- **Embedded tables** are written as `<model>-<table>.csv` (the CoolSolve
  companion-table convention, see
  [Language Reference §11](language_reference.md#11-lookup-tables)). Table
  names that are not safe file names (EES default `Lookup 1`) are renamed
  (`lookup_1`) and the quoted references in the equations are updated
  accordingly (logged in the report).
- **External tables**: EES functions can name a file instead of a table
  (`INTERPOLATE('Pipedata.lkt', …)`). Convert `.txt`/`.csv` files directly; the
  binary `.lkt` format is not decoded yet (open it in EES and *Save As* `.csv`,
  or extend the tool — library task T-TOOL).
- Known CoolSolve gaps (found while testing, see
  [model_library_support.md](model_library_support.md)):
  - the native EES call `INTERPOLATE('Table','Col1','Col2',Col2=value)`
    (returns `Col1` at `Col2=value`) is **not supported**: it is parsed as a
    property call and fails with "Unknown fluid" (`CS-GAP-INTERP-EES`);
    CoolSolve only accepts its own positional form
    `INTERPOLATE('table','xcol','ycol',x)`, with the column roles reversed;
  - EES `INTERPOLATE` is cubic (`INTERPOLATE1` is linear); CoolSolve is linear;
  - CoolSolve returns a **wrong value without error** when the x column is in
    descending order (`CS-BUG-INTERP-DESC`);
  - companion tables were not loaded when the model was given as a bare file
    name on the command line (`coolsolve model.eescode`): fixed after v0.3.0
    (`CS-BUG-LOOKUP-PATH`); with older builds use `coolsolve ./model.eescode`.

Following principle 1, models using the native EES `INTERPOLATE` form keep it
and are marked *blocked* until the gap is closed.

---

## 8. Parametric tables

EES parametric tables are extracted to `reference/ees_parametric-<t>.csv`
(numeric and string columns). They play two roles:

- **Reference results**: each row is a run; compare one CoolSolve run with
  `compare_solution.py … --row N`.
- **Inputs**: a variable may be defined *only* in the table (its equation
  commented out in the Equations window, e.g. `"fluid$='R22'"`). The extracted
  model is then not square or fails ("Unknown fluid: 'fluid'"). Give such
  variables the value of the first run in the equations (the "default run") and
  describe the table in the README.

CoolSolve parametric studies (GUI) can reproduce numeric sweeps, but there is
no parametric-table file in a model folder and string variables (fluid names)
cannot be swept (`CS-GAP-PARAMETRIC`).

---

## 9. Functions, procedures and libraries

- `FUNCTION` and `PROCEDURE` are supported; `MODULE`/`SUBPROGRAM` are not
  (`CS-GAP-MODULE`).
- The report lists the **functions called but not defined** in the file. Each
  must be either a CoolSolve built-in (see the
  [Language Reference](language_reference.md)) or a user-library function that
  EES loaded automatically (e.g. `cpbar`, `gamma` from the ULiège combustion
  library). The target is `$INCLUDE library:<name>`, which imports all the
  functions and procedures of a library model (`CS-FEAT-IMPORT`, see
  [model_library_support.md §2.3](model_library_support.md#23-importing-functions-and-procedures-from-eescode-files-cs-feat-import)).
  Until CoolSolve supports it, copy the required definitions at the top of the
  model, in a block marked
  `{--- Library functions copied from <library model ID> ---}`.
- **`.lib` files** (plain text, often starting with `{$DS.}`) contain only
  functions/procedures. They become a *function model* — a regular `.eescode`
  file, no `.lib` is kept: the definitions at the top, followed by a short main
  program that calls them with typical values (the pattern of
  [`examples/cpbar.eescode`](../examples/cpbar.eescode)). Other models then
  include it with `$INCLUDE library:<name>`.
- Compiled external functions (`.dlf`, `.dlp`, `.fdl`) cannot be imported; a
  model depending on them is *blocked* unless the algorithm can be
  re-implemented as an EES `FUNCTION`.

---

## 10. Dynamic and optimisation models

- **`INTEGRAL`** models are supported with the limitations of
  [Language Reference §12.7](language_reference.md#127-limitations) (one
  integration variable, constant limits, no table-based `INTEGRAL`).
  `$IntegralTable` columns must be space-separated.
- **Min/Max optimisation** is configured in the EES GUI (not in the text) and
  CoolSolve has no optimiser (`CS-GAP-OPTIM`). Import the model at the
  **optimum found by EES** (decision variables fixed to their stored values):
  it then runs as a steady model and the optimisation is documented as blocked.

---

## 11. Verification policy

Run `compare_solution.py` and paste its table (or a summary) in the model
README. Tolerances:

| Situation | Expected agreement |
|---|---|
| Fixed inputs, explicit algebra, no property call | exact (≤ 1e-9 relative) |
| Real-fluid properties, same equation of state | ≤ 0.1 % |
| Real-fluid properties, different equation of state (EES vs CoolProp) | ≤ 0.5 % typical |
| Older EES fluid models (e.g. R22 in EES 7) or near the critical point | up to a few %, **explain it** |
| Absolute enthalpy/entropy | may differ by a reference-state offset: compare differences and derived results (powers, efficiencies, COP, temperatures) |

A model is *verified* when all derived results agree within these tolerances;
the comparison and any explained deviation are recorded in the README.

---

## 12. Worked examples

All tests were run on 2026-10-04 with CoolSolve v0.3.0.

### Example 1 — refrigeration cycle with a simple compressor model (*verified*)

Source: `MSTh-SB-R1-Ex3.EES` (EES 7.458, ULiège course *Machines et systèmes
thermiques*). Now library model `CSL-0001`
(`models/cycles/refrigeration_heat_pumps/refrigeration_cycle_simple_compressor`).

1. `ees_extract.py` → 48-line model in French, unit system SI-C-Pa-J, 31
   variables, one parametric table (3 runs × 4 columns, with a string column
   `fluid$` = R22 / R134a / Propane), no lookup table, licence tag removed.
2. The report warned that most stored values are −9999: the file was last run
   from the parametric table. Solving the raw extraction failed with
   *"Unknown fluid: 'fluid'"*: `fluid$` was defined only by the parametric
   table (§8). Fix: `fluid$='R22'` as default run.
3. Comments translated, standard header added (no change to the equations).
4. The three runs were solved and compared with the parametric table
   (`compare_solution.py --row N`):

| Run | Fluid | COP (EES / CoolSolve) | Q̇_ev [W] (EES / CoolSolve) | max. rel. diff. |
|---|---|---|---|---|
| 1 | R22 | 5.218 / 5.264 | 163 575 / 165 952 | 1.4 % |
| 2 | R134a | 4.981 / 4.980 | 105 664 / 105 709 | 0.06 % |
| 3 | Propane | 5.102 / 5.107 | 139 709 / 140 117 | 0.3 % |

R134a and propane agree within the property tolerance. The 1.4 % deviation for
R22 is most likely due to differences between the R22 property formulations of
EES 7.4 and CoolProp (Kamei et al., 1995); it was not investigated further and
is documented in the model README.

### Example 2 — two-shaft gas turbine with a compressor map (*blocked*)

Source: `Revision JL051223-01.EES` (EES 7.458, revision exercise).

- Extraction: 45 variables with their stored solution, one embedded lookup
  table (`Lookup 1`, 3 rows × 4 columns: reduced speed, reduced flow, pressure
  ratio, isentropic efficiency) written to `<model>-lookup_1.csv`; references
  renamed to `'lookup_1'`.
- The report lists `cpbar` and `gamma` as **called but not defined**: they come
  from the ULiège combustion library (`CombCmHn_SI_PNG2003_V2.LIB`, also in
  `examples/cpbar.eescode`) that EES loads automatically.
- Gaps found: native `INTERPOLATE(…, N_rN=0.95)` not supported
  (`CS-GAP-INTERP-EES`); positional `INTERPOLATE` returns a wrong value on this
  table because its x column is descending (`CS-BUG-INTERP-DESC`: 370 instead of
  420); lookup tables ignored with a bare file name (`CS-BUG-LOOKUP-PATH`);
  library functions must be copied (`CS-FEAT-IMPORT`).
- Status: *blocked* by `CS-GAP-INTERP-EES` (kept in native syntax).

### Example 3 — unit conversion of a K / kPa / kJ exercise (*verified*)

Source: `R1_E3_2022.EES` (EES 10, course *Thermodynamique appliquée*
2022-2023): liquid methane evaporating in a tank, volume flow of the heated
vapour.

- Extraction: unit system `SI MASS DEG KPA K KJ` (warning issued), **decimal
  comma** detected and converted (`0,5` → `0.5`, `volume(methane;P=p_t;x=0)` →
  `volume(methane,P=p_t,x=0)`), `{$ST$OFF}` tag removed, 12 variables with
  their stored solution.
- Manual conversion (§6): `$UnitSystem SI MASS DEG PA C J`;
  `p_t = 500E3 [Pa] "500 kPa in the original"`, same for `p_r`;
  `T_r = 1.85 [C] "275 K in the original"`. Review of the equations: no energy
  term, no relation using an absolute temperature, no unit factor — nothing
  else to change. (A first attempt with `ConvertTemp(K, C, 275)` failed:
  `CS-GAP-CONVERT-UNQUOTED`.)
- Verification: `compare_solution.py methane_tank.sol reference/ees_variables.csv
  --ees-units` → all 12 variables agree (largest deviation 3·10⁻⁵ on the
  vapour specific volume; both tools use the reference equation of state of
  methane).

---

## 13. EES binary format notes (for tool developers)

Reverse-engineered on files written by EES 6.137 to 10.x (Delphi program:
little-endian integers, 80-bit *Extended* floats, length-prefixed *short
strings*).

| Offset / structure | Content |
|---|---|
| byte 0, bytes 1… | length + version string (`X7.458`, `X9.920`, `X10.589`; `W6.585` for `.lkt` files) |
| 16–19 | uint32: length *n* of the Equations-window text |
| 20 … 20+n | equations, Windows-1252 text or RTF |
| 20+n, 6 bytes | unit settings `[system, basis, T, P, trig, energy]`: basis 0 mass / 1 molar; T 0 °C / 1 K; P 0 kPa / 1 bar / 2 Pa (3 MPa assumed); trig 0 deg / 1 rad; energy 0 J / 1 kJ. Files older than EES 7 have no energy byte (kJ). |
| next 10 bytes | trigonometric factor (`0.017453293` for degrees, `1.0` for radians) — used to detect the 5/6-byte layout |
| variable record | short string *name*; relative to its length byte: value @+31, guess @+41, lower @+51, upper @+61 (Extended); format bytes @+71–76 (byte @+74 = 3); units short string @+77. Records are 122 bytes (EES 7) or 338 bytes (EES 9–10, with an upper-case copy of the name @+122). ±∞ bounds are stored as ±1E+3999; −9999 marks a cleared value. |
| table column | short string `name\nunits` in a 31-byte field, then *nrows* Extended values, then 6 trailing bytes (Lookup) or 3 (Parametric); string columns are followed by *nrows* × (uint32 len+1, short string). Columns of one table follow each other; the table name (`Lookup 1`, `Table 1`, or a user name) precedes the grid. |

Known limitations of `ees_extract.py`: `.lkt` files (EES 6 layout) not
decoded; EES 6 variable records decoded on a best-effort basis; the meaning of
some unit-setting codes is assumed (English units, MPa); Diagram windows,
plots, Min/Max and uncertainty settings are ignored; string-variable values are
not exported to `.initials`.

---

## 14. Checklist

- [ ] `ees_extract.py` run; report read (version, unit system, tables, functions called, warnings)
- [ ] Missing inputs restored (parametric-table variables, external lookup files, library functions)
- [ ] Unit system converted **by hand** to `SI MASS DEG PA C J` (§6), conversion logged
- [ ] Comments in English, standard header added, EES syntax otherwise unchanged
- [ ] Solved with CoolSolve (`-d` folder inspected if it fails, see [Debugging Models](debugging_models.md))
- [ ] Compared with the EES reference (`compare_solution.py`), deviations explained
- [ ] Unsupported features reported as gaps in [model_library_support.md](model_library_support.md)
- [ ] Only the model files moved to the library; source referenced by path/URL; work folder deleted
