# CoolSolve — Version history

[← Back to overview](../README.md)

This page lists all released versions of CoolSolve, their Windows installer
download links, the CoolProp version bundled in each release, and the main
changes introduced.

---

## v0.2 — April 2026 (current)

| Item | Detail |
|------|--------|
| **Windows installer** | **[CoolSolve v0.2 for Windows](https://dox.uliege.be/index.php/s/XXXXXXXXXXXXXX/download)** |
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
| **CoolProp** | master branch (version at time of release) |
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
