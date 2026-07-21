# CoolSolve — Version history

[← Back to overview](../README.md)

This page lists all released versions of CoolSolve, their Windows installer
download links, the CoolProp version bundled in each release, and the main
changes introduced.

---

## v0.3 — July 2026 (current)

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
