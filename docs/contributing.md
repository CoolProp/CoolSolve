# Contributing to CoolSolve

[← Back to overview](../README.md)

Thank you for considering a contribution to **CoolSolve**! This guide describes
the workflow that every new feature, enhancement, or bug fix should follow so
the codebase stays clean, consistent, and easy to maintain.

CoolSolve is a tightly integrated stack — a single new feature usually touches
the parser, the solver, the GUI, the configuration file, the tests, the
documentation, and the LaTeX/debug reports. Skipping any of these layers
leads to silent regressions or to a feature that "works on the CLI but is
invisible in the GUI". The checklists below are designed to prevent that.

---

## 1. Guiding Principles

Before writing any code, make sure you accept the following constraints. They
are non-negotiable and apply to every contribution.

1. **Prefer simple solutions.** Avoid clever code; readability beats brevity.
2. **No code duplication.** Before adding new logic, search for similar
   functionality elsewhere in the codebase and reuse or extend it.
3. **Only make changes that are requested.** Do not refactor unrelated code
   along the way; submit refactors as separate, focused contributions.
4. **Keep the codebase clean and organised.** Files have well-defined
   responsibilities (see §3) — respect them.
5. **Document everything that is non-obvious.** Use Doxygen-style comments
   on public C++ APIs and JSDoc-style comments on TypeScript exports.
   Avoid redundant comments that just restate the code.
6. **Never decrease computational efficiency.** CoolSolve is performance
   critical: CoolProp calls dominate runtime, and any new feature must be
   either zero-overhead by default or hidden behind an opt-in flag in
   `SolverOptions`. Profile before merging if you have any doubt.
7. **Backward compatibility.** New configuration keys must have sensible
   defaults so that existing models and `coolsolve.conf` files keep working
   unchanged.

---

## 2. Recommended Workflow

Every contribution should follow these phases in order. Do not skip phases.

```
┌──────────────┐   ┌─────────┐   ┌────────────────┐   ┌──────────────┐   ┌─────────────────┐
│ 1. Plan      │ → │ 2. Code │ → │ 3. Targeted    │ → │ 4. Full test │ → │ 5. Docs &       │
│  (checklist) │   │         │   │    tests       │   │    suite     │   │    review       │
└──────────────┘   └─────────┘   └────────────────┘   └──────────────┘   └─────────────────┘
```

### 2.1 Plan

Walk through the **integration checklist in §3** and write down, for each
applicable item, *what* needs to change and *why*. A typical contribution
touches between four and eight items in that list. The output of this phase
is a short Markdown plan you can paste into the pull request description.

### 2.2 Implement

Make the code changes. Follow §1 strictly. When introducing new public
identifiers (functions, options, REST endpoints, GUI components), keep names
consistent with the existing convention in the affected module.

### 2.3 Run targeted tests

Add or update Catch2 tests covering the new behaviour. Also add negative tests, 
where the feature should not work, and proper error reporting messages should
be provided. Then run only those tests first:

```bash
cd build && cmake --build . -j$(nproc)
./coolsolve_tests "[my-new-feature]"
```

Iterate on the implementation until the targeted tests pass. This loop is
much faster than the full suite.

### 2.4 Run the full test suite

Once the targeted tests pass, run the full suite to catch regressions:

```bash
# Unit tests (parser, evaluator, solvers, config, …)
./coolsolve_tests

# Comprehensive example-file tests (solves all .eescode examples and
# verifies expected values + LHS≈RHS for every equation)
./coolsolve_tests "[examples-comprehensive]"

# Solver robustness suite (every example × every pipeline configuration)
./coolsolve_tests "[solver-robustness]"
```

A contribution is **not ready for review** until all three of these pass on
your machine. Pay special attention to `examples/test_examples.md` and
`examples/solver_robustness_report.md` — large changes in iteration count or
runtime are red flags worth investigating before merging.

### 2.5 Update documentation

Even purely internal changes usually require a documentation pass — at the
very least the `README.md` features list and the language reference if a new
keyword/function/option is added. See §3 for the full list.

---

## 3. Integration Checklist

For **every** contribution, walk through this list and explicitly note which
items apply. Mark each item as **needed**, **not applicable**, or **already
covered**, and address every "needed" item before submitting the PR.

> The order below roughly follows the data flow: source code →
> parser → solver → outputs → user-facing surfaces. Following this order
> when implementing a feature minimises rework.

### 3.1 Language and Parser

If the feature introduces new syntax, keywords, operators, built-in
functions, fluids, or thermophysical properties:

- [ ] **`src/parser.cpp`** — extend the PEG grammar / lexer; add the new
      production and any keyword recognition (case-insensitive).
- [ ] **`include/coolsolve/ast.h`** — add new AST node kinds if needed.
- [ ] **`src/ir.cpp`** — handle the new AST node when building the IR
      (variable extraction, incidence matrix, LaTeX rendering).
- [ ] **`src/autodiff_node.cpp`** + **`include/coolsolve/autodiff_node.h`** —
      every new mathematical function MUST also implement its analytical
      derivative (forward-mode AD). Finite-difference fallbacks are
      considered a regression — see the principle in §1.7.
- [ ] **`src/variable_inference.cpp`** — if the feature implies a unit or
      a CoolProp property, extend the inference table.
- [ ] **`src/solution_checker.cpp`** — make sure new equation types are
      re-evaluated correctly during post-solve verification.
- [ ] **`tests/test_parser.cpp`** and **`tests/test_evaluator.cpp`** —
      add unit tests for parsing and AD propagation.

### 3.2 Solvers and Numerical Core

If the feature changes the solver pipeline, adds a new solver, or modifies
how blocks are decomposed:

- [ ] **`include/coolsolve/solver.h`** — add new fields to `SolverOptions`
      with sensible defaults (zero-overhead when disabled).
- [ ] **`src/solver.cpp`** — extend `loadSolverOptionsFromFile()` to parse
      the new keys from `coolsolve.conf`. **Forgetting this step makes the
      option silently ignored.**
- [ ] **`src/solver_*.cpp`** — new solvers go in their own translation
      unit (`solver_<name>.cpp`) following the pattern of
      `solver_newton.cpp`, `solver_lm.cpp`, …
- [ ] **`src/structural_analysis.cpp`** / **`src/solver_symbolic.cpp`** —
      update if block decomposition or symbolic reduction is affected.
- [ ] **`tests/test_solver_pipeline.cpp`**, **`tests/test_newton.cpp`**,
      **`tests/test_solver_robustness.cpp`** — add coverage. New solvers
      must be exercised by the robustness suite.

### 3.3 Static Configuration File

The canonical example file is `examples/coolsolve.conf`. Whenever
`SolverOptions` gains a new field:

- [ ] Add the new key to `examples/coolsolve.conf` **commented out** (so it
      shows the default), with a clear comment block explaining:
      what it does, when to enable it, and any trade-offs.
- [ ] Make sure `loadSolverOptionsFromFile()` parses it (§3.2).
- [ ] Add a regression test in `tests/test_config.cpp` to verify the key
      round-trips correctly.

### 3.4 GUI — Configuration Editor

The HTML configuration editor is the user-facing mirror of `coolsolve.conf`:

- [ ] **`gui/src/components/ConfigEditor.tsx`** — add a `ConfigField`
      entry for the new option, inside the appropriate `ConfigGroup`. Use
      a clear label, the correct type (`number`/`boolean`/`string`), the
      same default value as `SolverOptions`, and a description matching
      the comment in `examples/coolsolve.conf`.
- [ ] If the option belongs to the solver pipeline, update
      `PIPELINE_PRESETS` so the relevant presets remain consistent.
- [ ] Verify the change in dev mode (`npm run dev`) and in the built
      binary (`coolsolve --gui`).

### 3.5 GUI — Frontend (other components)

Beyond the config editor, several other GUI surfaces may need updates:

- [ ] **`gui/src/api/types.ts`** — extend the TypeScript types if the
      REST API response shape changes.
- [ ] **`gui/src/api/client.ts`** — add the corresponding API call.
- [ ] **`gui/src/stores/modelStore.ts`** / **`uiStore.ts`** — store new
      state fields (Zustand).
- [ ] **`gui/src/components/CodeEditor.tsx`** — only changes if editor
      behaviour itself changes.
- [ ] **`gui/src/languages/ees.ts`** — add new keywords, built-in
      functions, or fluid names so Monaco highlights them properly.
- [ ] **`gui/src/components/VariableTable.tsx`**,
      **`ArrayTable.tsx`**, **`ParametricStudy.tsx`**,
      **`ThermoDiagram.tsx`**, **`PlotlyChart.tsx`**,
      **`DebugViewer.tsx`**, **`Toolbar.tsx`**, **`Tooltip.tsx`** —
      update only the components that actually expose the feature.

### 3.6 Server / REST API

If the feature has to be reachable from the GUI:

- [ ] **`src/server.cpp`** — add or extend an HTTP route. Keep the
      handler thin: it should call into `CoolSolveRunner` exactly like
      `main.cpp` does, never reimplement core logic.
- [ ] **`include/coolsolve/server.h`** — update only if the public
      `ServerOptions` surface changes.
- [ ] Update SSE progress events if the feature changes solver progress
      semantics.
- [ ] Add at least one targeted test (Python or C++) that hits the new
      route end-to-end. Existing example: `tests/test_parametric_api.py`.

### 3.7 File Loading and Saving (ZIP Bundle)

The GUI exchanges complete model state as a ZIP bundle through
`/api/v1/files/upload` and `/api/v1/files/bundle`:

- [ ] **`src/server.cpp`** — if the feature persists user-editable data,
      add it to both `createZipBundle` (export) and `extractZipBundle`
      (import). Keep file extensions stable; add new ones rather than
      reusing existing ones.
- [ ] Confirm round-trip: download a bundle, upload it back, verify the
      session state is identical.
- [ ] If you change the bundle layout, document the new convention here
      and in `docs/gui.md`.

### 3.8 Debug Mode (`-d` flag)

CoolSolve's debug mode is the contract by which advanced internals become
inspectable. Anything that produces analytical output **must integrate
into the debug folder** rather than printing to stdout.

- [ ] **`src/runner.cpp` → `CoolSolveRunner::generateDebugOutput()`** —
      write a new `<feature>.md` file in the debug folder. Follow the
      existing pattern: a `# Title`, a brief explanation paragraph, then
      tables/sections with the analytical content.
- [ ] Update the **debug index** (the table written into `README.md` of
      the debug folder at the bottom of `generateDebugOutput()`) so the
      new file appears with a one-line description.
- [ ] Also update the **public debug-folder table** in the main
      `README.md` ("The debug folder contains:" section) so users can
      discover the new file from the project landing page.
- [ ] Markdown output in debug mode must be **deterministic** when
      possible (sorted keys, fixed precision) so debug folders can be
      diffed across runs.

### 3.9 stdout / stderr and GUI Console

CoolSolve uses three message channels with distinct rules:

| Channel               | Use for                                              |
|-----------------------|------------------------------------------------------|
| `stdout`              | The actual program output (JSON / LaTeX / text).     |
| `stderr`              | Progress, warnings, and human-readable errors.       |
| Debug-mode `.md` files| Detailed analytical artefacts (see §3.8).            |

Add user-facing messages with the right channel:

- [ ] User-visible errors must be emitted as **diagnostics**
      (`include/coolsolve/diagnostic.h` → `DiagnosticCollector`) so they
      appear in both the CLI stderr stream and the GUI Console with the
      correct severity.
- [ ] In the GUI, the `Console` component
      (`gui/src/components/Console.tsx`) classifies lines based on
      keywords (`ERROR`, `[Warning]`, `SUCCESS`, …). Use these markers
      consistently so colour-coding works.
- [ ] Never print verbose internals to stdout — that would corrupt the
      JSON/LaTeX output. Verbose internals belong in `-d` debug files.

### 3.10 LaTeX Report

The comprehensive LaTeX report is generated after a successful solve.

- [ ] **`src/latex_report.cpp`** + **`include/coolsolve/latex_report.h`** —
      if the feature contributes new analytical content (new equations,
      new variable categories, new diagnostic plots), include it in the
      report. The report must compile with plain `pdflatex` and degrade
      gracefully when optional plot files are missing (use
      `\IfFileExists`).
- [ ] **`tests/test_latex_report.cpp`** — add a regression test that the
      generated `.tex` contains the new section/macro and is well-formed.
- [ ] When working on the report locally, compile with:
      `pdflatex -synctex=1 -interaction=nonstopmode <file>.tex`.

### 3.11 Examples and Regression Models

CoolSolve's strongest regression net is the suite of `.eescode` examples.

- [ ] If the feature unlocks a new modelling pattern, add a minimal
      `examples/<feature_name>.eescode` that exercises it.
- [ ] Add an entry to `EXPECTED_SOLUTIONS` in
      `tests/test_examples.cpp` with a known target value and a 1 %
      tolerance. **Without this entry the example is parsed but its
      result is not validated.**
- [ ] If the example needs initial guesses, add a
      `examples/<feature_name>.initials` file.

### 3.12 README and Documentation

The `README.md` and the `docs/` folder are served verbatim by the GUI
binary (see `docs/documentation_strategy.md`), so any change there is
immediately visible to users.

- [ ] **`README.md`** — update the *Features* bullet list, the *Command
      Line Options* table, the *Project Structure* tree (if files were
      added), the *Documentation* table (if a new doc page was added),
      and the *File Formats* / *Debug folder* tables when applicable.
- [ ] **`docs/language_reference.md`** — required for any new keyword,
      operator, built-in function, fluid, or thermophysical property.
      Only document features that are implemented and tested.
- [ ] **`docs/debugging_models.md`** — extend if the feature introduces a
      new failure mode or a new diagnostic file.
- [ ] **`docs/solver_roadmap.md`** — update if a roadmap item is now
      delivered, or if the new feature replaces a planned approach.
- [ ] **`docs/gui.md`** — required for any GUI-visible change.
- [ ] **`docs/symbolic_redecomposition.md`** — update only if the
      symbolic-reduction algorithm itself changes.
- [ ] **`docs/docs.html`** — if you add a brand-new `.md` page to
      `docs/`, append it to the sidebar nav so the in-app docs viewer
      can find it.

### 3.13 Build System

- [ ] **`CMakeLists.txt`** — `src/*.cpp` and `tests/*.cpp` are picked up
      automatically by `GLOB_RECURSE`, so most contributions need no
      build-system change. Only update CMake when adding a **new
      dependency**, a **new build option**, or **non-glob sources**
      (e.g. resource files).
- [ ] If you add a dependency, prefer `FetchContent` and pin a tag.
      Update the *Dependencies* table in `README.md`.
- [ ] Verify the **Windows build path** still works: `build_installer.bat`
      and `coolsolve.nsi` should not need changes for normal feature
      additions, but should be regenerated and tested for installer-level
      changes.

---

## 4. Style and Coding Conventions

### C++ (C++17)

- Follow the existing file structure: declarations in
  `include/coolsolve/*.h`, definitions in `src/*.cpp`.
- Public APIs are commented with Doxygen blocks (`/** … */`).
- Prefer `const` and `noexcept` where they apply.
- Avoid raw pointers for ownership; use `std::unique_ptr` /
  `std::shared_ptr`.
- Use the existing `coolsolve::` namespace for all new symbols.
- New numerical algorithms must include a one-paragraph comment naming
  the source paper or textbook reference (see `solver_lm.cpp` for the
  pattern).

### TypeScript / React

- Functional components with hooks; no class components.
- State lives in Zustand stores (`gui/src/stores/`), never in component
  internal state if it crosses component boundaries.
- API calls go through `gui/src/api/client.ts`; do not call `fetch`
  directly from components.
- Prefer existing CSS variables (light/dark theme support) over
  hard-coded colours.

### Markdown / Documentation

- Use sentence case in headings.
- Wrap lines at ~80 characters where reasonable.
- Use code fences with language tags.
- Reference files with backticks (`` `src/solver.cpp` ``).

---

## 5. Performance Discipline

CoolSolve must remain a fast solver. The following rules apply to every
contribution:

1. **No new allocations in inner loops.** The Newton, LM, and
   TrustRegion solvers are called millions of times — any new code path
   reachable from `BlockEvaluator::evaluateBlock()` or the residual
   evaluation routines must be benchmarked.
2. **Default off for any new analysis pass.** If your feature adds
   pre-processing or post-processing work, gate it behind a
   `SolverOptions` flag that defaults to `false` (see
   `enableSymbolicReduction` for the canonical pattern). Document
   explicitly that "when disabled, zero overhead is added".
3. **CoolProp calls are expensive.** Cache derived quantities; use the
   thread-local `AbstractState` cache (`coolpropCacheEnabled`) rather
   than going through `PropsSI`.
4. **Compare runtimes before and after.** Run
   `./coolsolve_tests "[solver-robustness]"` and inspect
   `examples/solver_robustness_report.md` for per-example iteration and
   timing changes. Anything beyond ±10 % on the established models needs
   justification in the PR description.
5. **Release builds only for benchmarks.** Debug builds are 10–50× slower
   (see *Build Type: Release vs Debug* in the README); never benchmark
   in Debug mode.

---

## 6. Pull Request Checklist

Copy the block below into your PR description and tick each box.

```markdown
### Integration checklist (see docs/contributing.md §3)

- [ ] §3.1  Language / parser / AST / IR / AD / inference / solution checker
- [ ] §3.2  Solver pipeline / SolverOptions / loadSolverOptionsFromFile
- [ ] §3.3  examples/coolsolve.conf entry + test_config coverage
- [ ] §3.4  GUI ConfigEditor entry + presets
- [ ] §3.5  GUI components, API types & client, stores, EES Monaco language
- [ ] §3.6  REST endpoints in src/server.cpp
- [ ] §3.7  ZIP bundle round-trip (createZipBundle / extractZipBundle)
- [ ] §3.8  Debug-folder Markdown file + index entry
- [ ] §3.9  stdout / stderr / GUI Console wiring via Diagnostics
- [ ] §3.10 LaTeX report contribution + test_latex_report
- [ ] §3.11 New example .eescode + EXPECTED_SOLUTIONS entry
- [ ] §3.12 README + language_reference + relevant docs/*.md
- [ ] §3.13 CMakeLists.txt / dependency table (only if needed)

### Tests run locally

- [ ] `./coolsolve_tests` passes
- [ ] `./coolsolve_tests "[examples-comprehensive]"` passes
- [ ] (if GUI changed) Manual smoke test in `coolsolve --gui`

### Performance discipline (see §5)

- [ ] No new allocations in inner solver loops
- [ ] Any new pass is opt-in with default off
- [ ] Robustness report iteration counts within ±10 % of baseline
```

---

## 7. Quick Reference — File Map

When in doubt, the table below tells you where the canonical owner of a
concern lives.

| Concern                                  | Canonical file(s)                                    |
|------------------------------------------|------------------------------------------------------|
| Parsing CoolSolve syntax                 | `src/parser.cpp`, `include/coolsolve/ast.h`          |
| Equation graph / blocks                  | `src/ir.cpp`, `src/structural_analysis.cpp`          |
| Automatic differentiation                | `src/autodiff_node.cpp`                              |
| Block evaluation                         | `src/evaluator.cpp`                                  |
| Solver pipeline                          | `src/solver.cpp`                                     |
| Per-strategy solvers                     | `src/solver_<name>.cpp`                              |
| `coolsolve.conf` parsing                 | `loadSolverOptionsFromFile()` in `src/solver.cpp`    |
| Static config example                    | `examples/coolsolve.conf`                            |
| Solver options struct                    | `include/coolsolve/solver.h`                         |
| Post-solve verification                  | `src/solution_checker.cpp`                           |
| Debug-folder generation                  | `CoolSolveRunner::generateDebugOutput()` in `src/runner.cpp` |
| LaTeX report                             | `src/latex_report.cpp`                               |
| HTTP server / REST endpoints             | `src/server.cpp`                                     |
| ZIP bundle import/export                 | `createZipBundle` / `extractZipBundle` in `src/server.cpp` |
| GUI config editor                        | `gui/src/components/ConfigEditor.tsx`                |
| GUI Monaco syntax highlighting           | `gui/src/languages/ees.ts`                           |
| GUI Zustand stores                       | `gui/src/stores/modelStore.ts`, `uiStore.ts`         |
| GUI REST client                          | `gui/src/api/client.ts`, `types.ts`                  |
| GUI Console                              | `gui/src/components/Console.tsx`                     |
| Comprehensive example tests              | `tests/test_examples.cpp` (`EXPECTED_SOLUTIONS`)     |
| Solver robustness tests                  | `tests/test_solver_robustness.cpp`                   |
| Configuration parsing tests              | `tests/test_config.cpp`                              |
| Documentation strategy                   | `docs/documentation_strategy.md`                     |
| In-app docs sidebar                      | `docs/docs.html`                                     |

---

## 8. Reporting Issues

If you find a bug but do not have time to fix it, please file an issue
that includes:

1. The CoolSolve version (`./coolsolve --help` shows it on the first line).
2. The minimal `.eescode` reproducer (and `.initials` if relevant).
3. The output of `./coolsolve -d <model>.eescode` — the debug folder is
   the single most useful artefact for triage. Attach it as a ZIP.
4. The expected vs. observed behaviour.

For solver-convergence bugs, also attach the relevant section of
`solver_robustness_report.md` if available, and read
[Debugging Models](debugging_models.md) before opening the issue.

---

Thanks again for contributing! Following this guide keeps CoolSolve fast,
correct, and pleasant to maintain.
