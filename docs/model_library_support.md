# CoolSolve Support for the Model Library

The [CoolSolve Library](https://github.com/IntSusEnergySystems/CoolSolve_Library)
gathers hundreds of thermodynamic models written in the EES language. This page
specifies **how CoolSolve integrates that library** and collects what CoolSolve
still needs to run its models:

1. the integration design (§2): a **snapshot of the library embedded in the
   CoolSolve binary at compile time**, an **HTML model explorer** in the GUI,
   and the **import of functions and procedures from any `.eescode` file**
   (no `.lib` files);
2. the list of features with their implementation plan (§3);
3. the **gap register**: EES features that are missing or behave differently in
   CoolSolve, found while converting library models (§4), with the bugs found
   on the way (§5).

The register is filled by the people converting models (see
[ees_import.md](ees_import.md) and the library's `docs/model_workflow.md`):
each library model that cannot run lists the IDs of the gaps blocking it in
its `model.json` (`missing_features`), so the number of blocked models gives the
development priority. [`ees_vs_coolsolve.csv`](ees_vs_coolsolve.csv) remains the
exhaustive feature inventory; this page is the prioritised, library-driven view.

---

## 1. Design principles

- CoolSolve should read **native EES code**. Library models are kept in EES
  syntax; when CoolSolve differs from EES, the fix belongs in CoolSolve.
- A library model folder is a CoolSolve project folder ("the folder is the
  project", [GUI §7](gui.md)): `<model>.eescode`, `<model>.initials`,
  `coolsolve.conf`, companion tables `<model>-<table>.csv`, `<model>.sol`, plus
  `README.md`, `model.json` and `figures/`.
- **Each CoolSolve build carries the library as it was when it was compiled.**
  The desktop application works offline, the online demo serves exactly the
  same content, and the library version is part of the CoolSolve version.
- **Functions are shared through `.eescode` files.** Any `.eescode` file — a
  library model or a file of the user — can provide its `FUNCTION`s and
  `PROCEDURE`s to another model; no `.lib` files are needed.

---

## 2. Integration design

### 2.1 Build-time library snapshot (`CS-FEAT-LIB-EMBED`)

**Where the library comes from** (CMake cache variables, resolved in this
order):

| Variable | Use |
|---|---|
| `COOLSOLVE_LIBRARY_DIR` | Path of a local `CoolSolve_Library` checkout (development). |
| *(automatic)* `../CoolSolve_Library` | Sibling checkout next to the CoolSolve repository, if it exists. |
| `COOLSOLVE_LIBRARY_TAG` | Otherwise the library is fetched with `FetchContent` (same mechanism and `.fetchcontent_cache/` as CoolProp) at this tag or commit — default `main`; **releases pin a library tag**. |
| `COOLSOLVE_EMBED_LIBRARY` | `ON` by default when the GUI is built; `OFF` builds without library (the explorer then says so). |

**What is embedded** (allow-list, under the URL prefix `/library/`):
`library.json`, `functions.csv`, `taxonomy.json`, `docs/taxonomy.md`, and in
each model folder `README.md`, `model.json`, `*.eescode`, `*.initials`,
`*-*.csv`, `coolsolve.conf`, `*.sol` and `figures/*.{png,svg,jpg,webp}`.
Never `sources/`, `tools/`, `templates/` or `work/`.

**Validation and identity.** The build runs the library's
`tools/build_index.py --check` (Python is already required at configure time)
and stops on errors (warning only for a development build). `library.json`
carries the library commit or tag, its generation date, the number of models
and a **content hash** of all embedded files; CoolSolve exposes them
(`/api/v1/library`, `coolsolve --version`) and the explorer displays them
("Library snapshot 2027-03-01 · 312 models · abc1234"). `docs/versions.md`
records the library tag of each CoolSolve release, as it does for CoolProp.

**Mechanism.** Reuse `cmake/embed_assets.cmake`, which already embeds the GUI,
`docs/*.md` and `README.md` into a generated `embedded_assets.cpp` served by
`getEmbeddedAsset()`. Two extensions are needed: recursive extra directories
(today `EXTRA_DIRn` are scanned flat) with an extension allow-list, and a
**build step** instead of the configure-time `execute_process`
(`add_custom_command` depending on the library's `library.json`, whose content
hash changes whenever an embedded file changes), so that `cmake --build`
picks up a new library state. Size budget: text ≈ 15–20 kB per model and one
figure ≤ 100 kB per model (rule of the library workflow), i.e. ≈ 40 MB for 300
models; the build prints the snapshot size. If hex-encoding in CMake becomes
too slow at that size, generate the table with a Python script, or embed one
packed archive (C++26 `#embed`, `.incbin`, Windows resource) behind the same
lookup function.

**Run-time override.** `--library-dir <path>` (or `COOLSOLVE_LIBRARY_DIR`)
serves a library checkout from disk instead of the snapshot — for library
authors testing changes without rebuilding, or to use a newer library with an
older binary.

All accesses go through one C++ `LibraryStore` (index, taxonomy, function
index, file lookup by model ID, name or path) used by the server, the runner
(`$INCLUDE`) and the CLI.

### 2.2 HTML model explorer (`CS-FEAT-EXPLORER`)

- **Entry points**: a *Library* button in the toolbar next to *Examples*, and
  addressable views `#/library` and `#/library/CSL-0012` (shareable links on
  the online demo).
- **Explorer view**: taxonomy tree with model counts (`taxonomy.json`); search
  box over titles, summaries, tags, fluids and function names; facet filters —
  level 🟢🔵🟠🔴, kind ⚙️⏱️🎯🧩, status ✅☑️⚠️⛔❌📄, fluid, origin type; results
  as cards (title, badges, one-line summary, thumbnail of the README figure).
  Non-runnable models can be hidden by default.
- **Model page**: the README rendered with `marked.js` (already used by
  `docs.html`), relative image links rewritten to `/library/raw/<model path>/…`;
  a metadata box (ID, category, source and authors, CoolSolve version of the
  verification, gaps with links to §4); actions:
  - **Open** — loads `.eescode`, `.initials`, `coolsolve.conf`, companion tables
    and `.sol` into the session, like the *Examples* menu (with the *Back*
    snapshot), from the embedded files;
  - **Include functions** (models that define `FUNCTION`/`PROCEDURE`) — inserts
    `$INCLUDE library:<name>` at the top of the current model (§2.3);
  - **Copy code**, **Download ZIP** (bundle format of [GUI §7.2](gui.md)).
- **Server routes**: `GET /api/v1/library` (index and snapshot information),
  `GET /api/v1/library/taxonomy`, `GET /api/v1/library/functions`,
  `GET /library/raw/{path}` (files, as `/docs/raw/` for the documentation),
  `POST /api/v1/library/open` `{"id": "CSL-0012"}`.
- **CLI**: `coolsolve --library list [filter]`, `--library show <id>`,
  `--library get <id> [dir]` (writes the model folder), all from the snapshot.

### 2.3 Importing functions and procedures from `.eescode` files (`CS-FEAT-IMPORT`)

One directive, on its own line as in EES:

```ees
"Functions of a library model, by name or by permanent ID"
$INCLUDE library:cpbar_combustion_products
$INCLUDE library:CSL-0012
"Functions of a file of the user, relative to the model folder"
$INCLUDE my_correlations.eescode
```

- **What is imported**: from a `.eescode` target, **all its `FUNCTION`s and
  `PROCEDURE`s** (later `MODULE`s and `SUBPROGRAM`s); its main program — e.g. the
  demonstration equations of a library *function model* — is ignored. A
  library model can therefore be run on its own *and* serve as a function
  library; no `.lib` file is involved.
- **Resolution**: `library:<name>` and `library:CSL-NNNN` go through the
  `LibraryStore` (snapshot or `--library-dir`); IDs never change and names are
  unique, so includes survive any re-organisation of the library folders. Paths
  are resolved from the folder of the model, then from the `--library-dir`
  root.
- **Native EES files**, for compatibility only: `$INCLUDE file.txt` inserts
  the equations of the text file (EES behaviour) and `$INCLUDE file.lib` imports
  the functions and procedures of an EES library file; the optional `/R` flag of
  EES is accepted and ignored.
- **Nesting**: `$INCLUDE` lines inside an included `.eescode` file are followed
  (with cycle detection).
- **Conflicts**: a definition in the model wins over an included one (warning);
  two includes defining the same name are an error unless the definitions are
  identical.
- **Diagnostics**: unknown library name (with the closest names), file not
  found, errors reported with the included file name and line.
- **Reproducibility**: the `.sol` file and the debug report list the imported
  definitions with their origin (library model ID and snapshot hash, or file
  path).
- **Discovery in the editor**: auto-completion and hover for library functions
  (from `functions.csv`: name, signature, model ID) and a quick-fix on
  *unknown function `cpbar`*: "Include `library:cpbar_combustion_products`".
- **Portability** (`CS-FEAT-EXPORT-EES`): *Export for EES* writes a
  self-contained file in which every `$INCLUDE` is replaced by the imported
  definitions, so that the model also runs in EES; the ZIP bundle gets the same
  "flatten includes" option.
- **Implementation**: `src/parser.cpp` already has a `preprocess()` step and a
  generic `Directive` rule. Parse the main file, then parse each include target
  with the same parser, keep only its `FUNCTION`/`PROCEDURE` definitions (AST)
  and merge them before IR building (`src/ir.cpp`); the `LibraryStore`
  supplies the text of library targets.
- **Library conventions** that make this work: function models are plain
  `.eescode` files (definitions + demonstration program); function names are
  unique across the library (`tools/build_index.py` warns about duplicates);
  all library functions use SI-°C-Pa-J.

---

## 3. Features and implementation plan

| ID | Feature | Section | Priority |
|---|---|---|---|
| `CS-FEAT-LIB-EMBED` | Library snapshot embedded at build time, with version identity and run-time override | §2.1 | P1 |
| `CS-FEAT-EXPLORER` | HTML model explorer: taxonomy, search, facets, model pages, *Open* and *Include functions* | §2.2 | P1 |
| `CS-FEAT-IMPORT` | `$INCLUDE` of library models and `.eescode` files (all their functions and procedures); native `.txt`/`.lib` for compatibility | §2.3 | P1 |
| `CS-FEAT-LIB-TEST` | Library regression test: solve every library model whose status is *verified* or *runs* and compare with its `.sol` (logic of `test_examples.cpp`, driven by `library.json`), in the CoolSolve CI and in the library CI | – | P1 |
| `CS-FEAT-LIB-CLI` | `coolsolve --library list / show / get` | §2.2 | P2 |
| `CS-FEAT-EXPORT-EES` | *Export for EES*: flatten `$INCLUDE`s into a self-contained file; ZIP option | §2.3 | P2 |
| `CS-FEAT-PARAM-FILE` | Parametric-table file in the model folder (e.g. `<model>-parametric.csv`: one column per variable, one row per run), with **string columns** (fluid names) — needed to import EES parametric tables natively (`CS-GAP-PARAMETRIC`) | – | P2 |
| `CS-FEAT-DIAGRAM-IDEAL` | T-s / P-v / h-s diagrams without dome for ideal gases (`Air`, `N2`, `CO2`…) with the solved states — README figures of air-standard cycles, gas turbines and engines (meanwhile the library uses parametric sweep plots) | – | P2 |
| `CS-FEAT-PSYCHRO` | Psychrometric chart with the solved states of `AirH2O` models and PNG/SVG export — figures of the HVAC models (meanwhile parametric sweep plots) | – | P2 |
| `CS-FEAT-OVERLAY-STATES` | Let the user pick the variables of each state point in the Diagram tab (today: arrays or name suffixes), so that names such as `T_su_cp`, `h_ex_ev` can be drawn without adding equations | – | P2 |
| `CS-FEAT-METADATA` | "Save as library model": create the folder skeleton (README and `model.json` templates) from the current session | – | P3 |
| `CS-FEAT-PLOT-CLI` | Generate a diagram or a parametric plot from the CLI, to regenerate README figures when a model changes | – | P3 |

Suggested order of implementation:

1. **Snapshot and minimal explorer** — recursive allow-listed embedding,
   `library.json` identity, `/api/v1/library` and `/library/raw/`, explorer list
   + model page + *Open*. The library becomes visible in every build.
2. **Imports** — `$INCLUDE` for `.eescode` paths and `library:` targets,
   `/api/v1/library/functions`, *Include functions* button, editor quick-fix.
   Library models then use `$INCLUDE` instead of copied functions.
3. **Tooling** — library regression in CI, CLI `--library`, *Export for EES*,
   run-time `--library-dir`, release procedure pinning the library tag.

---

## 4. Gap register (EES features missing or different)

Priority rule: **P1** = silent wrong result, or blocks ≥ 5 library models;
**P2** = blocks 2–4 models or forces a model rewrite; **P3** = 1 model or
convenience. Update the *Blocked models* column when you mark a model blocked.

| ID | EES feature | CoolSolve behaviour | Workaround in the library | Priority | Blocked models |
|---|---|---|---|---|---|
| `CS-GAP-INTERP-EES` | `INTERPOLATE('T','Col1','Col2',Col2=v)` (and `INTERPOLATE1`, named argument selecting the input column, returns `Col1`) | Parsed as a property call: *"Unknown fluid: 'T'"*. Only CoolSolve's positional form `INTERPOLATE('t','xcol','ycol',x)` exists — column roles reversed w.r.t. EES. | None (keep EES syntax, mark blocked) | P1 | Potentially every model using `INTERPOLATE` (29 embedded lookup tables detected in the ULiège collection, before de-duplication); test model *two-shaft gas turbine* |
| `CS-GAP-INTERP-CUBIC` | `INTERPOLATE` is cubic in EES (`INTERPOLATE1` linear, `INTERPOLATE2D…`) | Linear only | Accept deviation if small, document | P3 | – |
| `CS-GAP-UNITSYSTEM` | `$UnitSystem` and the EES unit settings (kPa, bar, MPa, kJ, K, molar, radians) | Directive parsed but ignored: a model in kPa/kJ/K gives wrong results **silently**. Quick win: warn when the directive differs from `SI MASS DEG PA C J`. | Library models are converted **by hand** to SI-C-Pa-J ([ees_import.md §6](ees_import.md#6-unit-system--manual-conversion)); native support is not needed by the library | P1 (warning) / P3 (native support) | 44 % of the ULiège EES files (before de-duplication) are not in SI-C-Pa-J |
| `CS-GAP-PARAMETRIC` | Parametric tables as model input (variables defined only by table columns, string columns, `TABLERUN#`, `TABLEVALUE`) | No parametric-table file; GUI sweeps are numeric only | Default run in the equations, table documented in the README | P2 | – |
| `CS-GAP-BOUNDS` | Variable lower/upper bounds (Variable Information, `$VarInfo`) | Not supported (`.initials` holds guesses only) | Better guesses, simplified-model bootstrap | P2 | – |
| `CS-GAP-MODULE` | `MODULE` / `SUBPROGRAM` | Parse error ("not yet handled") | Rewrite as PROCEDURE only if trivial; otherwise blocked | P2 | – |
| `CS-GAP-INCLUDE` | `$Include` / `$Load`, implicit `USERLIB` functions | Not supported | Copy the definitions until `CS-FEAT-IMPORT` (§2.3) is available; then `$INCLUDE library:<name>` | P1 | test model *two-shaft gas turbine* (`cpbar`, `gamma`) |
| `CS-GAP-OPTIM` | Min/Max optimisation (single and multi-variable) | Not available | Import at the EES optimum (decision variables fixed) | P2 | – |
| `CS-GAP-LKT` | External binary lookup files (`.lkt`) and file-based tables in `INTERPOLATE('file.txt',…)` | Not readable | Convert to `<model>-<table>.csv` (by EES or a future decoder) | P3 | – |
| `CS-GAP-FLUIDS-ABS` | `NH3H2O`, `LiBrH2O`, `X_LIBR`/`H_LIBR` (absorption machines) | Unsupported fluids/functions (CoolProp has no such mixtures) | External correlations as EES FUNCTIONs | P2 | CoolSolve example `water_libr` |
| `CS-GAP-CONVERT-UNQUOTED` | Unit names without quotes in `CONVERT(kPa, Pa)` and `CONVERTTEMP(K, C, T)` (quotes are optional in EES) | The unit names are taken as variables: system not square | Add the quotes (still valid EES) | P2 | – |
| `CS-GAP-DECIMAL-COMMA` | Files written with the European number format: `0,5` and `;` as argument separator | Parse errors / "Unknown fluid: ''" | `ees_extract.py` converts to the dot convention automatically | P3 | – |
| `CS-GAP-INTEGRAL-LOOKUP` | Lookup functions (`INTERPOLATE`…) inside an `INTEGRAL` model | Lookup store not wired into the `IntegralSolver` ([integral_table.md §7.1](integral_table.md#71-current-limitations)) | Analytic expressions instead of tables | P2 | CoolSolve example `building_rc_network` |
| `CS-GAP-PROC-COMMENT` | Comment strings inside procedure argument lists, e.g. `PROCEDURE p(a, b "<= inputs" : c "outputs =>")` (accepted by EES) | Parse error | Remove the comments (allowed edit: comments only) | P2 | CoolSolve example `orc_complex` |
| `CS-GAP-CALL-WARNING` | `CALL WARNING('message', value)` (and `CALL NOTE`) in functions and procedures, used by the ULiège correlation toolkits to flag out-of-range inputs, e.g. `IF (Re<1E4) THEN CALL WARNING('Re must be above 1E4. A value of xxxE4 was supplied', Re)` | `CALL ERROR` exists (`CS-GAP-UNITSYSTEM-FUNC`, closed) but `CALL WARNING` raises *"Unknown procedure: warning"* when its branch is taken | Keep the call (valid EES); a branch that is never taken is harmless | P3 | – (found while surveying the collection: ~12 calls in 2 files) |
| `CS-GAP-IDEAL-ENTROPY` | Absolute entropies of the ideal-gas species (JANAF `s°` at 25 °C, 1 atm: N2 191.61, O2 205.04, CO2 213.80, CO 197.65, H2O(g) 188.83 kJ/kmol·K), needed for Gibbs-energy, chemical-equilibrium and exergy calculations with species | `entropy(CO2,T=…,P=…)` keeps CoolProp's per-fluid reference state (only the enthalpy got the EES reference, `CS-GAP-FORMATION-ENTHALPY`) | Differences of entropy are right; avoid absolute entropies of species | P3 | – |
| `CS-BUG-IDEAL-INTENERGY` | The alias `intenergy` (and `internalenergy`) is not in the list of temperature-only properties of ideal gases: `intenergy(CO2,T=25)` fails with *"Ideal gas property 'intenergy' requires pressure input"* while `enthalpy`, `cp`, `cv` accept a single input (found while testing `CS-GAP-FORMATION-ENTHALPY`) | – | Add `P=` to the call | P3 | – |

## 5. Bugs and documentation issues found while importing

| ID | Description | Reproducer | Priority |
|---|---|---|---|
| `CS-BUG-LOOKUP-HINT` | The error hint says *"place 'data.csv' in the model directory"* whereas the convention is `<model>-data.csv`. | same as above | P3 |
| `CS-BUG-UNIT-VOLUME` | Unit inference labels specific volumes from `VOLUME(...)` as `kg/m^3` in the `.sol` file (should be `m^3/kg`). | `v = VOLUME('R22', T=8, P=5E5)` → `v = 0.049 "kg/m^3"` | P3 |
| `CS-BUG-SOL-STRING` | String variables appear twice in the `.sol` file: once as a number (`fluid$ = 0.0`) and once as a string (`fluid$ = 'R22'`). | any model with `fluid$ = 'R22'` | P3 |
| `CS-DOC-EESCSV` | `ees_vs_coolsolve.csv` still lists `INTERPOLATE`, `LOOKUP`, `NLOOKUPROWS`… as *No* although they are implemented (Language Reference §11). | – | P3 |

---

## 6. Reporting a new gap

1. Check this register and [`ees_vs_coolsolve.csv`](ees_vs_coolsolve.csv).
2. Write a **minimal reproducer** (≤ 10 lines of EES code) and run it with the
   current CoolSolve build.
3. Add a row to §4 (missing/different feature) or §5 (bug), with a new ID
   `CS-GAP-<SHORT>` / `CS-BUG-<SHORT>`, the workaround used in the library and
   the IDs of the blocked models.
4. Reference the ID in the model's `model.json` (`missing_features`) and README.
5. When a gap is closed: update `ees_vs_coolsolve.csv`, move the row to a
   *Closed* table with the CoolSolve version, and re-test the blocked models
   (the library's roadmap has a task for that).

---

## 7. Closed gaps and bugs

Entries moved here when fixed (§6, step 5). *Fixed in* gives the CoolSolve
version, the branch and the commit; the regression test lives in
`tests/test_library_gaps.cpp` unless stated otherwise (run with
`./coolsolve_tests "[library-gaps]"`).

| ID | What was wrong | Fixed in | Regression test | Notes |
|---|---|---|---|---|
| `CS-BUG-IF-IGNORED` | Conditions of `IF … THEN … ELSE` inside FUNCTION/PROCEDURE bodies were ignored (every branch executed, later assignments overriding earlier ones). The parser never built `IfThenElse` nodes; multi-line `{…}` comments of bodies were not skipped either, so commented-out `IF` lines would now run. | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-BUG-IF-IGNORED: …` (10 cases); example `if_then_else_function` | Single-line and block forms with `ENDIF`, nested IF, `=` `<>` `<` `>` `<=` `>=`, `AND`/`OR`, string conditions; an `IF … THEN` in the main program is now a parse error (EES allows it only in functions/procedures). Existing examples changed: `zorlu_heat_pump` (pinch procedures now work; COP unchanged), `engine_weibe_cycle` (`x_b`). See [language_reference §6](language_reference.md#6-control-flow). |
| `CS-BUG-MOLARMASS` | `MOLARMASS(fluid)` returned kg/mol (CoolProp's `M`) where EES returns kg/kmol: silent 1000× error in native EES code (gas constants in J/kmol·K, stoichiometric ratios). | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-BUG-MOLARMASS: …` (2 cases) | Now kg/kmol for every fluid. Existing examples: `air_screw_compressor*` (`r = 8314/MM_air`; the cross-check `C_thr_leak_bis` now agrees with `C_thr_leak`), `cpbar` (`c_bar_p` 1093 J/kg·K, EES 1090.6), `orc_*` (`mm` unused). |
| `CS-BUG-SINGLE-INPUT-PAIR` | Two single-input property calls of the same ideal gas in one expression (`enthalpy(X,T=T1) − enthalpy(X,T=T2)`) evaluated inconsistently between the solve and the verification passes: the solver clamps every pressure to ≥ 1000 Pa, which also hit the internal pressure injected for single-input calls of `H2O` (100 Pa). Residual 6.3e3 on the register reproducer (exactly the 100 Pa vs 1000 Pa difference of Δh), no `.sol` written. | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-BUG-SINGLE-INPUT-PAIR: …` (4 cases) | The injected ideal-gas pressure is no longer clamped nor checked against the kPa/Pa unit hints (it was reported as a user input: *"p=100 is Pa, not kPa"*). Existing examples: only `cpbar` (`c_bar_p` +1.5e-4). The register reproducer uses `molarmass(H2O)` without the `1000*` (kg/kmol since `CS-BUG-MOLARMASS`). |
| `CS-GAP-UNITSYSTEM-FUNC` | `UNITSYSTEM('K')` was an unknown function and `IF … THEN CALL error('msg')` was not supported: the unit-system guards of EES `.lib` files had to be commented out. | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-GAP-UNITSYSTEM-FUNC: …` (4 cases) | `UNITSYSTEM('name')` follows `coolpropConfig.units` (SI, Mass, Deg always 1; C, Pa, J by default; 0 for Eng, Molar, Rad, K, kPa, kJ…); unknown names are an error. `CALL ERROR('msg' [, values])` throws `CALL ERROR: msg`. A user procedure named `error` still wins. `CALL WARNING` (also in the ULiège toolkits) is not implemented: see `CS-GAP-CALL-WARNING`. |
| `CS-GAP-FORMATION-ENTHALPY` | Ideal-gas species enthalpies lacked the heats of formation (EES/JANAF convention): CoolProp's per-fluid reference states are arbitrary, so `LHV_CO = h_CO·MM_CO + 0.5·h_O2·MM_O2 − h_CO2·MM_CO2` gave −5.53 MJ/kmol instead of 282.99 MJ/kmol, silently (blocked `Q_4` of CSL-0005). | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-GAP-FORMATION-ENTHALPY: …` (4 cases, incl. `examples/cpbar.eescode` against the stored EES solution) | `IdealGasFluid` carries the heat of formation (Çengel Table A-26: CO2 −393 520, CO −110 530, H2O(g) −241 820, CH4 −74 850, C2H6 −84 680, C3H8 −103 850, N2/O2/H2 0 kJ/kmol — they reproduce EES's `Q_4` to 4·10⁻⁵ of the LHV of CO); the enthalpy offset is computed once per CoolProp backend and applied to results *and* to `h=`/`u=` inputs. `Air`, `He`, `Ar` and real fluids are unchanged; entropies keep CoolProp's reference (see `CS-GAP-IDEAL-ENTROPY`). Existing examples: only `cpbar` (`Q_4` = 1.518e8, EES 1.517e8). |
| `CS-GAP-MULTILINE-COMMENT` | A `"…"` comment trailing an equation and continuing on the next line was a parse error on each line (*"Could not parse line"*); only standalone comment lines that wrapped were skipped. Found with library model `rankine_cycle_60mw`. | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-GAP-MULTILINE-COMMENT: …` (4 cases); `test_parser.cpp` *quote comment trailing an equation* | Any `"` not closed on its line opens a comment that runs to the next `"` (code before it is kept); same rule in the main program and in function/procedure bodies (shared `skipCommentLine`). A comment never closed raises warning `P007` with its opening line. The library model `rankine_cycle_60mw` can have its two joined comments un-joined (checked on a re-wrapped copy: identical solution). |
| `CS-BUG-LOOKUP-PATH` | Companion tables `<model>-<table>.csv` were ignored when the CLI got a bare file name (`coolsolve model.eescode`): `loadLookupTableForModel()` returned early because `parent_path()` is empty. (Works with `./model.eescode` or an absolute path.) | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-BUG-LOOKUP-PATH: …` | The current directory is used when the path has no directory part. The library workflow's *note the ./* remark is no longer needed with a fixed build. |
| `CS-BUG-INTERP-DESC` | `INTERPOLATE` returned a wrong value **without error** when the x column was in descending order (EES accepts both orders): the interval search assumed ascending rows. Register reproducer: `N_rN = 1, 0.95, 0.9` / `M = 454, 420, 370`, `INTERPOLATE('t','N_rN','M',0.95)` gave 370 instead of 420. | post-v0.3.0, branch `fix/library-gaps`, commit pending | `CS-BUG-INTERP-DESC: …` (3 cases) | The direction of the x column is detected (first vs last row); nodes, interior points, flat extrapolation at both ends and the derivative (Newton on the x column) are right in both orders. Ascending tables are bit-identical. The native EES call form is still `CS-GAP-INTERP-EES`. |
| `CS-DOC-TRIG` | [Language Reference §2](language_reference.md) said `sin/cos/tan` take radians; the implementation uses degrees (`sin(30) = 0.5`), consistent with the EES `DEG` setting. | post-v0.3.0, branch `fix/library-gaps`, commit pending | documentation only (checked: `sin(30)`, `cos(60)`, `arcsin(0.5)`, `atan2(1,2)` = 26.57°) | The inverse functions return degrees too; the page now says how to convert models written for the EES `RAD` setting. |
