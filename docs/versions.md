# CoolSolve — Version history

[← Back to overview](../README.md)

This page lists all released versions of CoolSolve, their Windows installer
download links, the CoolProp version bundled in each release, and the main
changes introduced.

---

## Unreleased

### Security

- **LaTeX compilation**: `POST /api/v1/latex/compile` ran the compiler named in the request body, or in
  `latexCompiler` of the session's `coolsolve.conf`, through the shell, so a client could run any command on the
  server. Only `pdflatex`, `xelatex` and `lualatex` are accepted now (anything else: `400`), plot image names must
  be plain `.png` names, and on Linux and macOS the compiler runs with `-no-shell-escape`. A `latexCompiler`
  set to a full path or with extra arguments is now refused: use the bare name and put the TeX binaries on `PATH`.
- **Listening address**: the GUI server listens on `127.0.0.1` instead of every network interface. The new option
  `--host <address>` restores network access (`--host 0.0.0.0`) and prints a warning, since the GUI has no
  authentication. Reverse-proxy deployments that proxy to `127.0.0.1` need no change.
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

---

## v0.4 — October 2026 (current)

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
- **Ammonia-water (`NH3H2O`)** *(if included in the release)*: the IAPWS Guideline G4-01 Helmholtz formulation, with
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
