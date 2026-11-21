# CoolSolve — Version history

[← Back to overview](../README.md)

This page lists all released versions of CoolSolve, their installer
download links (Windows; macOS and Linux since v0.5), the CoolProp version bundled in each release, and the main
changes introduced.

---

## Next release (unreleased)

Notes of the changes merged since v0.5, gathered for the next changelog.

### Security

- **LaTeX compilation**: `POST /api/v1/latex/compile` ran the compiler named in the request body, or in
  `latexCompiler` of the session's `coolsolve.conf`, through the shell, so a client could run any command on the
  server. Only `pdflatex`, `xelatex` and `lualatex` are accepted now (anything else: `400`), plot image names must
  be plain `.png` names, and on Linux and macOS the compiler runs with `-no-shell-escape`. A `latexCompiler`
  set to a full path or with extra arguments is now refused: use the bare name and put the TeX binaries on `PATH`.
- **Cross-site requests**: the `Access-Control-Allow-Origin: *` header is gone, and a `POST`, `PUT` or `DELETE`
  whose `Origin` names another site is refused (`403`). Scripts that send no `Origin` are not affected.
- **DNS rebinding**: the server answers only to requests addressed to `localhost`, an IP address, or a name given
  with the new option `--allowed-host <name>` (others: `403`). **Action for server deployments**: a reverse proxy
  that forwards the browser's host name (Apache `ProxyPreserveHost On`, as in the
  [Deployment Guide](deployment_ubuntu_apache.md)) needs `--allowed-host <your domain>` added to the service
  command line, or every page returns `403`.
- **LaTeX report**: the index of an array variable (`X[...]`) is now escaped in the generated `.tex` like the rest
  of the name.
- **Usage statistics**: `GET /api/v1/stats/log` and the `/stats` dashboard show client networks (`topNetworks`,
  IPv4 /24 and IPv6 /48) instead of full client addresses (`topIps` is removed). The log file is unchanged.

### New

- **Conditions in `IF(cond, a, b)`**: the condition is written as the one of an `IF ... THEN` statement: `=` and
  `<>` beside the orderings, strings compared with `=` and `<>` (`IF(layout$ = 'staggered', a, b)`), `AND`, `OR`
  and parentheses (GitHub issues #3 and #4). `IF(a = 5, 1, 2)`, until now a parse error, tests `a = 5`.

### Fixed

- **String variables written twice**: the `.sol` and `.initials` files and the variable table of the GUI listed
  each string variable also as the number 0 (`fluid$ = 0.0` beside `fluid$ = 'R22'`); it is now written once, as
  its string (GitHub issue #4). A tool that compares `.sol` files should ignore the `name$ = 0.0` lines of older
  ones, as the library regression of CoolSolve now does.
- **Strings in `IF`**: `IF(s$, 'banana', 10, 20, 30)` read both strings as 0 and always took the "equal" branch,
  without any message. `IF` compares and returns numbers, as in EES: a string where `IF` expects a number is now a
  parse error, which names the condition that compares the strings, `IF(s$ = 'banana', X, Y)` (GitHub issue #9).

---

## v0.5 — November 2026 (current)

| Item | Detail |
|------|--------|
| **Windows installer** | **[CoolSolve v0.5 for Windows](https://github.com/CoolProp/CoolSolve/releases/download/v0.5/CoolSolve_v0.5_Installer.exe)** (Windows 10/11, x64) |
| **macOS disk image** | **[CoolSolve v0.5 for macOS](https://github.com/CoolProp/CoolSolve/releases/download/v0.5/CoolSolve_v0.5_macOS_universal.dmg)** (macOS 12+, Apple Silicon and Intel; not signed by Apple, see *Installing*) |
| **Linux package** | **[CoolSolve v0.5 for Debian/Ubuntu](https://github.com/CoolProp/CoolSolve/releases/download/v0.5/coolsolve_0.5.0_amd64.deb)** (Ubuntu 22.04+, Debian 12+, x64) |
| **Release notes** | [GitHub release v0.5](https://github.com/CoolProp/CoolSolve/releases/tag/v0.5) |
| **CoolProp** | master branch, commit `8eb3b778` (8.1.0dev, 31 October 2026) |
| **Model library** | [CoolSolve Library](https://github.com/IntSusEnergySystems/CoolSolve_Library) tag `v0.5` (336 models, 307 verified), embedded in the GUI binary |
| **Online demo** | [https://coolsolve.squoilin.eu/](https://coolsolve.squoilin.eu/) |

### Installing

- **Windows**: run the installer. It is not signed: SmartScreen shows *Windows protected your PC*, click *More info*
  → *Run anyway*. CoolSolve is in the Start Menu, and a double click on a `.eescode` file opens it.
- **macOS** (new): open the `.dmg` and drag CoolSolve to Applications. The app is not signed by Apple: macOS refuses
  the first opening. Go to *System Settings → Privacy & Security* and click *Open Anyway* (repeat if macOS asks
  again), or run once in Terminal `xattr -dr com.apple.quarantine /Applications/CoolSolve.app`. CoolSolve opens a
  Terminal window (the server; close it to stop CoolSolve) and the GUI in the browser. The command-line tool is
  `/Applications/CoolSolve.app/Contents/MacOS/coolsolve`.
- **Linux** (new): `sudo apt install ./coolsolve_0.5.0_amd64.deb`. *CoolSolve* is in the applications menu (Science)
  and opens `.eescode` files; `coolsolve` is the command-line tool, `coolsolve --gui` the GUI.

See [Desktop Application](desktop_app_options.md) for how the packages are built.

### Behaviour changes

Model results can move: re-check models that use the items below.

- **Ideal gases are exactly ideal** (`Air`, `N2`, `O2`, `H2`, `He`, `Ar`, `CO2`, `CO`, `CH4` … `C8H18`, `H2O`): the
  ideal-gas part of their CoolProp equation of state, `Pv = RT`, `h(T)`, `u = h − RT`, at any T and P, as in EES.
  They were real-gas CoolProp states: results move by the real-gas residual (typically 0.1–1 %); `C8H18` or `H2O`
  at room conditions are the gas, never the liquid; the transport properties are those of the dilute gas (a given P
  is ignored); `quality`, `T_sat`, `P_sat` of an ideal gas are errors naming the real fluid
  ([Language Reference §5](language_reference.md)).
- **`QUALITY` outside the saturation dome is 100 (vapour side) or -100 (liquid side), as EES**, instead of 0, for
  every input pair; supercritical states are 100 above the critical temperature and -100 below it. A model that
  compares the quality with 0 or 1 outside the dome must account for it (*Quality outside the saturation dome*).
- **`AND`/`OR` in IF and UNTIL conditions are read from left to right**, as EES (`AND` had priority).
- **A `CALL ERROR` reached by any iterate stops the calculation**, as EES; a guarded model may need guesses that
  keep the iterates out of the guarded range ([Debugging Models](debugging_models.md)).
- **`CALL psychprops` returns v, h, s, u per kg of dry air** (as the AirH2O functions and EES) and accepts named
  inputs (`T=`, `P=`, `R=`, …).
- **Convergence test relative to the size of each equation's terms**: tiny values (1e-20) are solved instead of
  being set to 0, and large energy balances converge at the rounding level of their terms; new warning `W002` for
  an implicitly solved value below the tolerance (*Convergence tolerance*, Language Reference §10).
- **The GUI server listens on `127.0.0.1` only** (it listened on every interface, `0.0.0.0`): the GUI is no
  longer reachable from other machines of the network. `coolsolve --gui --host 0.0.0.0` restores it; a reverse
  proxy on the same machine (the online demo) needs nothing ([GUI §12.1](gui.md#121-local-mode)).
- **Try Harder also turns on equation scaling** (`enableEquationScaling`, below), with tearing and symbolic
  reduction: a model that needed *Try Harder* may converge to another solution of a block with several.
- **`VOLUME` is labelled `m^3/kg`**, and `$IntegralTable` rows are at multiples of the output interval
  (`integralOutputInterval` is now applied).

### New

- **Desktop packages for macOS and Linux**, next to the Windows installer: a `.dmg` with `CoolSolve.app` (universal
  binary) and a `.deb` package with a menu entry, the `.eescode` file type and the command-line tool. The three
  packages are built from the release commit by the GitHub Actions workflow `packages.yml`
  ([Desktop Application](desktop_app_options.md)).
- **Thermodynamic diagrams on CoolPlot**: the Diagrams tab is rebuilt on the C++ core of
  [CoolPlot](https://github.com/CoolProp/CoolPlot) (coolplot-core 1.5.0, CoolPlot `7338f48`) and on
  `CoolProp::Plot` of the pinned CoolProp (`8eb3b778`): isolines in the model's units, short
  diagram names (p-h, T-s, …), themes, zoom; the solved states and any number of arrays overlays
  with paths along the processes, and the hover; pure fluids, mixtures,
  incompressible solutions, ideal gases and psychrometric and Mollier charts of humid air; light
  isolines by default; SVG/PDF downloads, vector figures in the LaTeX report,
  `coolsolve --diagram`. The charts are saved in the ZIP bundle (`diagrams.json`); the former
  `POST /api/v1/coolprop/saturation` is removed ([CoolPlot integration](coolplot_integration.md)).
- **`$INCLUDE <name>`**: a model imports the functions and procedures of another `.eescode` file of its folder or
  of a library model (e.g. `$INCLUDE hx_ees_function`), in the CLI and the GUI alike; the Library menu copies the
  line, an include bar opens each included file, the ZIP bundle carries them
  ([Language Reference §7](language_reference.md#include-functions-and-procedures-of-another-file)).
- **Combustion**: `OH`, `H`, `O`, `N`, `NO` and `NO2` are ideal gases from the NASA-Glenn polynomials (200 K to
  20000 K); `CALL NASA(Species$, T : cp, h, s)` on the same table of 22 species; `LowerHeatingValue(Fuel)` and
  `HigherHeatingValue(Fuel)` of the ideal-gas fuels of C, H and O; `CALL CHEM_EQUIL(P, T, A\O, C\O, H\O, N\O : …)`,
  the chemical equilibrium of 14 combustion species by minimum Gibbs energy, with exact derivatives: adiabatic flame
  temperatures with dissociation are solved directly ([Language Reference §5, §7](language_reference.md)).
- **Option `enableEquationScaling`** (default off; *Try Harder* turns it on): Newton steps on the equilibrated
  system (rows and columns scaled by powers of two), for blocks whose equations have terms of very different sizes,
  such as the adiabatic flame with dissociation (*Equation scaling*, Language Reference §10).
- **Dynamic models**: a variable-order BDF integrator for stiff models (`integralMethod = BDF`); faster evaluation
  (variables by slot instead of by name); the check equations that read `t` or a state are verified during the
  march; the step budget of RK45 and BDF follows the step of the model ([Dynamic Solving](integral_table.md)).
- **Open reads the ZIP files of any zip tool**: Deflate-compressed archives (the default of `zip`,
  7-Zip, Windows and macOS), data descriptors, ZIP64, the ZIP of a folder and Windows' code
  page 437 names; before, only uncompressed archives opened, and a Deflate one ended in a bare
  *Internal Server Error*. A ZIP that cannot be opened (password, BZIP2/LZMA, damaged, a file that
  is not UTF-8 text, …) is refused with the reason, and an uncaught server error reaches the
  Console with its message ([GUI §7.2](gui.md#72-zip-bundle-format)).
- **Debug folder**: `solver_final_residuals.md`, the residual of every equation of the block that did not converge
  at its solver's last iterate, sorted worst first; `coolsolve_residuals.md` after a failure is evaluated where the
  solve stopped instead of at the guesses ([Debugging Models](debugging_models.md)).
- **GUI**: `F2` solves the model. The server listens on `127.0.0.1` only, `--host <address>` otherwise (see above).
- **EES import tools** (`tools/ees_extract.py`, `tools/compare_solution.py`) no longer drop or invent data silently:
  names with `\ % | # °` and accents, the range comma of `SUM`/`PRODUCT`/`DUPLICATE`, tables after the first block
  of records ([Importing EES Models](ees_import.md)).
- **Library regression test** `[library-regression]`: every runnable library model against its `.sol` baseline.

### Fixed

- **`x^y` with a negative base** (`q^2`, `q^3` with `q < 0`) had a NaN derivative, which sent blocks to the fallback
  solvers.
- **Tearing**: a bounded backtracking line search on the outer step (the full step could overshoot to a non-physical
  root).
- **Newton1D** no longer accepts a bracket collapsed onto a jump of F (a failed CoolProp call at the edge of its
  domain was taken for a root).
- **Dynamic models**: a Jacobian whose evaluation fails leaves the previous one intact.

### Dependencies

- CoolProp pinned to `8eb3b778`; coolplot-core 1.5.0 pinned like CoolProp (`COOLSOLVE_COOLPLOT_PINNED_TAG`);
  nlohmann/json 3.12.0; CMake 3.16. A server drawing PDF diagrams needs a TrueType font (`fonts-dejavu-core`).

---

## v0.4 — October 2026

| Item | Detail |
|------|--------|
| **Windows installer** | **[CoolSolve v0.4 for Windows](https://github.com/CoolProp/CoolSolve/releases/download/v0.4/CoolSolve_v0.4_Installer.exe)** ([release notes](https://github.com/CoolProp/CoolSolve/releases/tag/v0.4)) |
| **CoolProp** | master branch, commit `75af7816` (8.1.0dev, October 2026) |
| **Model library** | [CoolSolve Library](https://github.com/IntSusEnergySystems/CoolSolve_Library) tag `v0.4` (248 models), embedded in the GUI binary |
| **Online demo** | [https://coolsolve.squoilin.eu/](https://coolsolve.squoilin.eu/) |

### What's new in v0.4

- **Model library in the GUI**: the CoolSolve Library — documented, verified EES models — is embedded in the binary.
  The toolbar **Library** button opens a whole model folder (code, initial values, lookup tables, configuration) in
  one click; the dashboard at `/library/` lists every model with its status, reference and figure.
- **EES compatibility, driven by the library** (about 80 register items closed, see
  [model_library_support.md](model_library_support.md) §7): `IF` ladders and `GOTO` in routine bodies, the EES
  forms of `INTERPOLATE`, `LOOKUP`, `LOOKUP$` and `LOOKUP$ROW`, indexed `SUM`/`PRODUCT`, `END PROCEDURE`, `=>` and
  `=<`, `UNTIL` without parentheses, `Uppercase$`/`Lowercase$`, `ISIDEALGAS`, `MOLARMASS` of an element,
  `CALL Colebrook`, the time functions of `SETP.LIB`, `CALL WARNING`, `$IF` directives, identifier characters
  (`%`, `|`, `#`, accents), `CONVERT` with compound units, array ranges `X[a..b]`, `DUPLICATE` with variable bounds
  and index arithmetic.
- **No silently wrong answers**: unparsable arguments, trailing text, dropped body lines, lookup-table cell writes,
  never-ending `REPEAT` loops, call cycles, unit hints on constants and clamped trial points are now reported with
  an error or a warning instead of being ignored or clamped silently; explicit equations store the exact
  right-hand side; `--no-sol` no longer deletes an existing `.sol`.
- **More fluids**: the normal alkanes `C4H10` … `C8H18` as ideal-gas substances; ideal gases with absolute
  third-law entropies and NASA-Glenn formation enthalpies; the ASHRAE lithium bromide-water routines `X_LIBR`,
  `H_LIBR`, `T_LIBR`, `P_LIBR`; the incompressible solutions in the EES form `cp('EG', T=T, C=30)`; humid-air `cv`
  and density as EES; CoolProp mixtures as fluids (with a warning).
- **Ammonia-water (`NH3H2O`)**: the IAPWS Guideline G4-01 Helmholtz formulation, with
  15 flash input combinations, `CALL NH3H2O` and the usual property functions; see
  [language_reference.md](language_reference.md) and [nh3h2o_review.md](nh3h2o_review.md).
- **Solver**: step control in the one-dimensional Newton solver (no more silent divergence); the input clamp only
  acts on Newton trial points and a clamp at the solution is an error.
- **Build and release**: CoolProp pinned to a tested commit (`coolsolve --version` prints it and the library
  commit); the Windows installer script checks `git` and `cmake`, names the installer after the version and takes the
  library tag from `LIBRARY_TAG`.

---

## v0.3 — July 2026

| Item | Detail |
|------|--------|
| **Windows installer** | **[CoolSolve v0.3 for Windows](https://github.com/CoolProp/CoolSolve/releases/download/v0.3/CoolSolve_v0.3_Installer.exe)** ([release notes](https://github.com/CoolProp/CoolSolve/releases/tag/v0.3)) |
| **CoolProp** | master branch, commit `a540d91` (8.0.0-dev, June 2026) |
| **Online demo** | [https://coolsolve.squoilin.eu/](https://coolsolve.squoilin.eu/) |

### What's new in v0.3

- **Equation-based dynamic solving (`INTEGRAL`)**: solve initial-value
  differential–algebraic equation (DAE) models written in the EES integral form
  `y = y0 + INTEGRAL(dydt, t, t0, tf)`. Fixed-step solvers (Euler explicit,
  Euler implicit, **RK4** — the default) and a variable-step adaptive
  **Dormand–Prince RK45** are provided, with optional Richardson
  extrapolation on fixed steps. Coupled ODEs plus algebraic variables
  (semi-explicit index-1 DAE) are supported; the algebraic subsystem is
  solved at every step by the existing `Solver`, unmodified.
- **`$IntegralTable` directive**: declare which variables to tabulate and the
  output interval (`$IntegralTable t:0.1  y  dydt`), including `X[1..5]`
  range expansion. The trajectory is written to an auto-generated
  `<modelname>-integral.csv`, embedded in the solve JSON, and shown in a new
  **Integral** tab in the bottom panel (table + Plotly line plot + CSV export).
- **`INTEGRALVALUE`**: parser recognises the function; the tabulated trajectory
  is interpolated via `IntegralTable::interpolate`.
- **Configuration**: nine `integral*` keys in `coolsolve.conf` (method, fixed
  step, max steps, rel/abs tol, min/max step, Richardson, output interval),
  all inert by default — zero overhead on non-integral models.
- **Round-trip**: the integral CSV travels in the ZIP bundle
  (`<modelname>-integral.csv` exported and restored on upload).
- **Debug**: new `integral.md` report + `integral_table.csv` copy in the debug
  folder.
- **Multi-start solver**: retry a failed block from alternative starting points
  (CoolProp-consistent seeds for thermo blocks, scale-based seeds for algebraic
  blocks). Configurable policy via `multiStartMode` (`always` / `deepsearch` /
  `never`), with parallel candidate execution (`multiStartNumCores`,
  first-to-converge wins).
- **KINSOL-style solver** (`solver_kinsol.cpp`): three opt-in modes —
  Dennis–Schnabel line search, Picard fixed-point iteration, and
  Anderson-accelerated fixed-point. Selectable via `Kinsol` in
  `solverPipeline`.
- **Trust-Region + hybrd Broyden reuse**: opt-in quasi-Newton Jacobian reuse
  inside the trust-region dogleg solver (`trustRegionUseHybrdBroyden`,
  `trustRegionBroydenMaxReuse`), mirroring Powell's hybrd algorithm.
- **Try Harder button (GUI)**: after a failed solve, the *Solve* button morphs
  into *Try Harder*; clicking it re-runs the model with the full Deep Search
  pipeline (`deepSearchPipeline`), tearing and symbolic reduction forced on,
  and the configured multi-start policy. Editing the model, initials, or
  configuration restores the normal *Solve* button.
- **User hints module** (`user_hints`): targeted, actionable warnings for the
  most common modelling mistakes (e.g. inconsistent units, missing initials,
  CoolProp input conflicts), surfaced in both the CLI and the GUI Console.
- **Better error reporting for malformed expressions** in the parser.
- **Persistence of `coolsolve.conf` on the server**: the configuration file is
  now part of the ZIP bundle round-trip, so the same model + config travel
  together between the CLI and the GUI.
- **Progress / cancellation bug fix**: internal solves (integral steps,
  tearing sub-solves, multi-start candidates) no longer emit progress events
  that froze the GUI; cancellation is now responsive during long integrations.
- **Frontend lint pass**: all `eslint` errors resolved; Monaco editor helpers
  refactored into `gui/src/components/editorUtils.ts`.

---

## v0.2 — April 2026

| Item | Detail |
|------|--------|
| **Windows installer** | **[CoolSolve v0.2 for Windows](https://dox.uliege.be/index.php/s/xCcllAEs4Y8hsnx/download)** |
| **CoolProp** | master branch, commit `bc9341f` (7.2.1-dev, April 2026) |
| **Online demo** | [https://coolsolve.squoilin.eu/](https://coolsolve.squoilin.eu/) |

### What's new in v0.2

- **Lookup tables**: External CSV data files loadable from equations via
  `INTERPOLATE`, `INTERPOLATE2`, `LOOKUP`, `TABLEVALUE`, and related
  aggregation functions (`SUMLOOKUP`, `AVGLOOKUP`, …); GUI **Lookup Tables**
  panel for in-browser creation and editing.
- **Symbolic block reduction**: Optional pre-processing pass that shrinks
  algebraic blocks via explicit extraction, CoolProp call inversion, and
  equation substitution, with automatic re-decomposition of the reduced block
  into independent sub-blocks (`enableSymbolicReduction`).
- **Multi-dimensional bisection (BisectionND)**: Derivative-free sign-change
  solver for small blocks, robust when the Jacobian is singular or zero
  (`bisectionNDMaxBlockSize`, default n ≤ 8).
- **Homotopy continuation**: New solver strategy for difficult starting points
  where gradient methods fail.
- **Non-monotone line search**: Configurable memory window for the Newton line
  search (`lsNonMonotoneMemory`).
- **Solution verification**: Post-solve independent re-evaluation of every
  equation, checking LHS ≈ RHS within a configurable tolerance; enabled
  automatically in debug mode and in the comprehensive test suite.
- **Analytical CoolProp derivatives**: `first_partial_deriv()` integration
  with forward-FD consistency check and automatic fallback near phase
  boundaries.
- **GUI**: Parametric study panel, Thermo diagram, Plotly charts, debug
  viewer, lookup table editor, and ZIP bundle import/export.
- **Version string**: `coolsolve --help` now reports the version on the first
  line.

---

## v0.1 — initial release

| Item | Detail |
|------|--------|
| **Windows installer** | **[CoolSolve v0.1 for Windows](https://dox.uliege.be/index.php/s/Iy6Tl0iDvKPCCve/download)** |
| **CoolProp** | master branch (8fa873a, March 2026) |
| **Online demo** | [https://coolsolve.squoilin.eu/](https://coolsolve.squoilin.eu/) |

### What was in v0.1

- Parser for the EES-compatible `.eescode` language (scalars, arrays,
  operators, built-in functions, `FUNCTION`/`PROCEDURE` blocks, `CALL`,
  `IF-THEN-ELSE`, `DUPLICATE`, `REPEAT-UNTIL`, comments, directives, units
  annotations, CoolProp property calls).
- Structural analysis: Hopcroft-Karp matching, Tarjan SCCs, block
  decomposition.
- Forward-mode automatic differentiation with exact analytical derivatives.
- Solver pipeline: explicit size-1 solve, Newton1D, Newton + line search,
  Trust-Region Dogleg, Levenberg-Marquardt, partitioned block updates.
- Parallel solver execution (first-to-converge).
- CoolProp integration via the low-level `AbstractState` API with thread-local
  caching.
- JSON and LaTeX output formats; `.sol` and `.initials` file support.
- Debug mode (`-d`): comprehensive output folder with Markdown analysis files.
- Embedded web GUI (`--gui`) with Monaco editor, variable table, and
  configuration editor.
- REST API for programmatic access.
