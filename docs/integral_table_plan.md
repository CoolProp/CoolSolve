# Implementation Plan — `INTEGRAL` / `$IntegralTable` (Equation-Based Dynamic Solver)

[← Back to overview](../README.md)

This document is the **planning blueprint** for adding the EES
equation-based `INTEGRAL` function and the `$IntegralTable` directive to
CoolSolve. It is the result of reading the EES help pages
(`integration_parameters.htm`, `2_ygicd.htm`, `_integraltable.htm`,
`_integralautostep.htm`, `richardson_extrapolation.htm`,
`integralvalue.htm`, `hs830.htm`, `multi_variable_integration.htm`,
`finding_a_limit_of_integration.htm`) and of a thorough exploration of
the existing CoolSolve codebase.

**Status:** Phases 0–10 delivered (2026-07-09). See §12 *Progress Log*.
Phase 11 (documentation) pending.

Reference for conventions: [`docs/contributing.md`](contributing.md).

---

## 1. Goals and scope

### 1.1 Goal

Enable CoolSolve to solve **initial-value differential–algebraic
equation (DAE) systems** written in the EES integral form:

```ees
"State variables are declared via the equation-based INTEGRAL function."
y = y_0 + integral(dydt, t, t_low, t_high)     "declares a state variable y"
dydt = ...                                      "derivative (algebraic in y, t, ...)"
y_0 = 1                                         "initial condition"

$IntegralTable t:0.1  y  dydt                   "tabulate trajectory"
```

and to **tabulate the trajectory** of selected variables in an Integral
Table, mirroring EES.

### 1.2 In scope (this plan)

| Feature | EES source | Notes |
|---------|-----------|-------|
| Equation-based `INTEGRAL(integrand, var, lo, hi)` | `2_ygicd.htm` | Core. Auto step. |
| Equation-based `INTEGRAL(..., lo, hi, step)` | `2_ygicd.htm` | Core. Fixed step. |
| `$IntegralTable` directive | `_integraltable.htm` | Lists variables to tabulate + optional output interval. |
| `INTEGRALVALUE(t,'X')` function | `integralvalue.htm` | Retrieve a past tabulated value by interpolation. |
| Richardson extrapolation (fixed step) | `richardson_extrapolation.htm` | Opt-in option. |
| Fixed-step solvers: **Euler explicit**, **Euler implicit**, **RK4** | this plan | RK4 = recommended default. |
| Variable-step solver: **Adaptive RK45 (Dormand–Prince)** | this plan | Non-stiff ODEs; error-estimation step control. |
| Integral-table output (CLI JSON, `.sol`, debug folder, **auto CSV file**) | this plan | CSV filename based on the model name, mirroring the parametric-table workflow. |
| **Integral Table GUI tab** (bottom panel) + CSV export + ZIP bundle round-trip | this plan | Same UX as the parametric study: a tab at the bottom of the screen showing the trajectory, persisted with the model. |
| Configuration via `coolsolve.conf` (`SolverOptions`) | this plan | **No `$IntegralAutoStep` / `$IntegralStop` directives** — parameters live in the config file (per user preference). |

### 1.3 Out of scope (deferred, with placeholders in code)

| Feature | Reason | Placeholder |
|---------|--------|-------------|
| `$IntegralAutoStep` / `$IntegralStop` directives | User prefers `coolsolve.conf`. | Parsed-and-warned only; not interpreted. |
| Multi-variable (nested) integration | EES warns it is "much slower"; rare in practice. | Detection + clear error message. |
| Table-based `INTEGRAL(integrand,var)` (2-arg) | Needs an EES-style Parametric table that CoolSolve does not have. | Detection + clear error message. |
| Index reduction (Pantelides / dummy derivatives) | Only needed for high-index DAEs; index-1 semi-explicit systems work without it. | Commented placeholders in the integrator; structural detection raises a clear "high-index system not supported in this version" error. |
| BDF / DASSL-style stiff variable-step solver | User selected adaptive RK45. Stiff systems may need many steps. | `Integrator` base class designed so a future `IntegratorBDF` slots in. |

### 1.4 Guiding design constraints (from `contributing.md` §1)

1. **Leave existing code untouched where possible.** The integral
   subsystem is a new, self-contained module. The only edits to existing
   files are: (a) minimal parser hooks to recognise the new syntax, (b)
   a dispatch point in the runner, (c) additive `SolverOptions` fields,
   (d) additive documentation rows.
2. **Distinct files**, in a dedicated folder (`src/integral/`,
   `include/coolsolve/integral/`) so the algebraic solver and the
   dynamic solver are clearly separated.
3. **Reuse existing methods** for initialisation, block decomposition,
   and algebraic solving — the integrator calls the existing `Solver`
   at each time step rather than reimplementing Newton/LM/etc.
4. **Tests first (TDD).** Every new numerical method gets a dedicated
   Catch2 test against a known analytical solution *before* the
   integrator is wired into the pipeline.
5. **Zero overhead by default.** If a model has no `INTEGRAL` call, no
   integral code path is executed; performance is identical to today.

---

## 2. Mathematical model

### 2.1 The EES integral convention

EES expresses a first-order ODE as an integral equation:

```
y = y0 + INTEGRAL(f, t, t0, tf)      ⟺  dy/dt = f,  y(t0) = y0
```

A model may contain several such equations sharing the **same**
integration variable `t` and the **same** interval `[t0, tf]`
(coupled ODEs). Any variable that does **not** appear as the LHS of an
integral equation is *algebraic* (it is solved, at each time step, from
the remaining equations).

### 2.2 Semi-explicit index-1 DAE

The general form CoolSolve will handle in this first version is:

```
y' = f(t, y, z)        ← differential equations  (state variables y)
0  = g(t, y, z)        ← algebraic equations     (algebraic variables z)
```

where `z` can be solved uniquely from `g` given `(t, y)` (index-1
condition). This covers the vast majority of transient thermal models
(heat exchangers, tank dynamics, control loops) without requiring index
reduction.

**Higher-index systems** (where an algebraic constraint couples to a
derivative that must be differentiated to solve) are **detected and
rejected** with a clear message in this version; placeholders are left
for future Pantelides-style index reduction (§9).

### 2.3 How each method advances one step `h` from `(t_n, y_n, z_n)`

For all methods, evaluating `f` or `g` at a candidate point requires
solving the algebraic subsystem `g(t̂, ŷ, ẑ) = 0` for `ẑ` with the
existing algebraic `Solver`. The integrator never re-implements Newton.

| Method | Step formula | Algebraic solves per step | Notes |
|--------|-------------|---------------------------|-------|
| Euler explicit | `y_{n+1} = y_n + h·f(t_n, y_n, z_n)` | 1 (for `z_n`) | Cheapest; 1st order. |
| Euler implicit | `y_{n+1} = y_n + h·f(t_{n+1}, y_{n+1}, z_{n+1})`, solve `(y_{n+1}, z_{n+1})` together | 1 coupled nonlinear solve | Stable for stiff systems; 1st order. Reuses block solver with `y_{n+1}` as unknowns. |
| RK4 | 4 stages, classic Butcher tableau | 4 | 4th order; recommended default. |
| RK45 (Dormand–Prince) | 6 stages embedded pair (4th + 5th) | 6 | Adaptive; `h` adjusted from error estimate. |
| Richardson (option, fixed step only) | Combine one `h` step + two `h/2` steps: `I ≈ (4·I(h/2) − I(h))/3` | triples the cost | Reduces truncation error by 1–2 orders (EES docs). |

---

## 3. Architecture

### 3.1 New file layout

```
include/coolsolve/integral/
├── integral_problem.h      IntegralProblem struct (state vars, derivatives, ICs, interval)
├── integral_table.h        IntegralTable (time-series columns + interpolation)
├── integrator.h            Integrator abstract base + IntegratorOptions
├── integral_solver.h       IntegralSolver (orchestrates the time-march loop)
└── integral_extraction.h   extractIntegralProblem(IR) -> IntegralProblem

src/integral/
├── integral_problem.cpp
├── integral_table.cpp
├── integrator_euler_explicit.cpp    Fixed step, explicit Euler (+Richardson)
├── integrator_euler_implicit.cpp    Fixed step, implicit Euler (+Richardson)
├── integrator_rk4.cpp               Fixed step, classic RK4 (+Richardson)
├── integrator_rk45.cpp              Variable step, Dormand-Prince RK45
├── richardson.cpp                   Richardson extrapolation helper (shared)
├── integral_solver.cpp              Time-march loop, reuses algebraic Solver
└── integral_extraction.cpp          Walk IR, classify state vs algebraic vars

tests/
├── test_integrator_euler.cpp        TDD: Euler explicit/implicit vs analytical
├── test_integrator_rk4.cpp          TDD: RK4 vs analytical
├── test_integrator_rk45.cpp         TDD: RK45 + adaptive step control
├── test_richardson.cpp              TDD: Richardson error reduction
├── test_integral_extraction.cpp     TDD: IR -> IntegralProblem classification
├── test_integral_table.cpp          TDD: IntegralTable store/interpolate/toCSV
└── test_integral_e2e.cpp            End-to-end: .eescode -> trajectory -> CSV

gui/
├── src/components/IntegralTable.tsx     Bottom-panel tab (mirrors ParametricStudy.tsx)
├── src/api/types.ts                     + IntegralTableData type
└── src/languages/ees.ts                 highlight integral / $IntegralTable / integralvalue
```

**Why a dedicated folder?** The contributing guide (`§7`) notes that
per-strategy algebraic solvers live in `src/solver_<name>.cpp`. The
integral subsystem is a *different kind* of solver (time-marching, not
single-shot root finding). Keeping it under `src/integral/` with a
matching `include/coolsolve/integral/` namespace makes the separation
explicit and avoids confusion. Both folders are picked up automatically
by the existing `GLOB_RECURSE` in `CMakeLists.txt` (`CMakeLists.txt:141`,
`CMakeLists.txt:265`), so **no build-system change is needed**.

### 3.2 How it connects to the existing pipeline

```
                 ┌──────────────────────────────────────────────────────────┐
 main.cpp /      │  CoolSolveRunner::run()                                  │
 runner.cpp      │   parse → IR → infer → analyze                           │
                 └───────────────┬──────────────────────────────────────────┘
                                 │
                  ┌──────────────▼───────────────┐
                  │  NEW: hasIntegral(IR)?        │   ← cheap check: any
                  │  yes → IntegralSolver path    │     INTEGRAL() call?
                  │  no  → existing Solver path   │     (zero overhead)
                  └──────┬────────────────┬───────┘
                         │ no             │ yes
                         ▼                ▼
              ┌────────────────┐  ┌─────────────────────────────────────┐
              │  Solver::solve │  │  NEW: IntegralSolver::solve         │
              │  (unchanged)   │  │   extract IntegralProblem           │
              └────────────────┘  │   for each time step:               │
                                  │     fix (t, y_n) as externals       │
                                  │     call Solver::solve on algebraic │
                                  │       subset (REUSED, unmodified)   │
                                  │     read back f (= derivatives)     │
                                  │     advance y by Integrator.step()  │
                                  │     record into IntegralTable       │
                                  │   return SolveResult + table        │
                                  └─────────────────────────────────────┘
```

The key invariant: **`Solver` is called unmodified at each time step**.
The integral layer only decides *when* to call it, *what* to fix as
external, and *how* to advance the state.

### 3.3 Key new types (signatures, illustrative)

```cpp
// include/coolsolve/integral/integral_problem.h
struct StateVariable {
    std::string name;                 // e.g. "y"
    std::string integrandVar;         // e.g. "dydt" (the f in dy/dt = f)
    double initialValue;              // y(t0)
    ExprPtr integrandExpr;            // optional: direct expression for f
};

struct IntegralProblem {
    std::string integrationVar;       // "t"
    double lowerLimit, upperLimit;    // [t0, tf]
    double fixedStep = 0.0;           // 0 => use auto/variable step
    std::vector<StateVariable> states;
    std::vector<std::string> algebraicVars;   // solved at each step
    std::vector<size_t>      algebraicBlocks; // block ids of g equations
    IntegralTableSpec tableSpec;      // from $IntegralTable
    bool valid = false;
    std::string errorMessage;
};

IntegralProblem extractIntegralProblem(const IR& ir,
                                       const AnalysisResult& analysis);
```

```cpp
// include/coolsolve/integral/integrator.h
struct IntegratorOptions {
    enum Method { EulerExplicit, EulerImplicit, RK4, RK45 };
    Method method = RK4;                  // default
    double fixedStep = 0.0;               // 0 => derive from maxSteps
    int    maxSteps = 1000;
    double relTol = 1e-6;                 // RK45 error control
    double absTol = 1e-9;
    double minStep = 0.0, maxStep = 0.0;  // bounds (0 = auto)
    bool   richardson = false;            // fixed-step only
};

class Integrator {
public:
    // One step. rhs(t, y) -> f. Returns y_{n+1} and the step actually taken.
    virtual StepResult step(double t, const Eigen::VectorXd& y,
                            const RHSFunction& rhs,
                            const IntegratorOptions& opt) = 0;
    virtual ~Integrator() = default;
};

// Factory: returns the concrete integrator for opt.method.
std::unique_ptr<Integrator> createIntegrator(IntegratorOptions::Method);
```

```cpp
// include/coolsolve/integral/integral_solver.h
struct IntegralSolveResult {
    SolveResult algebraicResult;                 // final-step algebraic solve
    IntegralTable table;                         // trajectory
    int totalSteps = 0;
    std::vector<double> acceptedStepSizes;       // for debug/diagnostics
    bool success = false;
    std::string errorMessage;
};

class IntegralSolver {
public:
    IntegralSolver(const IR& ir, const AnalysisResult& analysis,
                   const SolverOptions& opt);
    IntegralSolveResult solve(const IntegratorOptions& intOpt);
};
```

### 3.4 RHS evaluation strategy

At time `t` with state `y`, the RHS `f(t, y)` is evaluated by:

1. Write `t` and the components of `y` into the `SystemEvaluator`'s
   current state as **fixed/external** values (the same mechanism the
   parametric sweep uses — `server.cpp:1848` caches analysis and re-solves
   with overrides).
2. Call `Solver::solve()` restricted to the algebraic blocks only.
3. Read the integrand variables (`dydt`, ...) from the solution — those
   are `f`.

For **implicit Euler**, step (1–3) is wrapped inside the nonlinear solve
for `y_{n+1}`: the block solver receives `y_{n+1}` as additional
unknowns appended to the algebraic block. This reuses the existing
Newton/LM/TrustRegion machinery.

---

## 4. Phase-by-phase implementation plan

> **Rule:** each phase is independently testable. Do not start a phase
> until the previous phase's targeted tests pass. See `contributing.md §2`.

### Phase 0 — TDD tests for the numerical integrators (no EES yet)

**Files:** `tests/test_integrator_euler.cpp`,
`tests/test_integrator_rk4.cpp`, `tests/test_integrator_rk45.cpp`,
`tests/test_richardson.cpp`.

Write these **first**, against the `Integrator` interface in §3.3, using
pure-lambda RHS functions (no IR, no evaluator). Reference problems with
known analytical solutions:

| Problem | Analytical | Used to test |
|---------|-----------|--------------|
| `dy/dt = -y, y(0)=1` on `[0,4]` | `y = e^{-t}` | convergence order of every method |
| `dy/dt = t² , y(0)=0` on `[0,1]` | `y = t³/3` | polynomial exactness |
| Coupled: `y'=z, z'=-y` (harmonic oscillator) | `sin/cos` | multi-state coupling |
| Stiff-ish: `dy/dt = -1000·y` | `e^{-1000 t}` | shows why implicit Euler + RK45 step shrinkage matter |
| RK45 tolerance ladder: same problem at `relTol = 1e-3, 1e-6, 1e-9` | — | error scales as `tol` |

These tests pin down the numerics **before** any parser/IR coupling.

### Phase 1 — `Integrator` implementations

**Files:** `src/integral/integrator_euler_explicit.cpp`,
`integrator_euler_implicit.cpp`, `integrator_rk4.cpp`,
`integrator_rk45.cpp`, `richardson.cpp`, plus headers.

- Each file implements `Integrator::step()`.
- A `createIntegrator(Method)` factory mirrors the existing
  `createSolver(SolverStrategy)` pattern (`solver.cpp:170`).
- Richardson extrapolation is a wrapper: it runs the base integrator
  once at `h` and twice at `h/2`, then combines `I ≈ (4·I_{h/2} − I_h)/3`.
  Applicable only when the step is fixed.
- One-paragraph comment at the top of each file citing the reference
  (Hairer & Wanner; Press et al. *Numerical Recipes*; Dormand & Prince
  1980) — per `contributing.md §4`.

Phase-0 tests now compile and pass.

### Phase 2 — `IntegralTable` + `INTEGRALVALUE`

**Files:** `src/integral/integral_table.cpp` + header,
`tests/test_integral_table.cpp`.

- `IntegralTable` stores named columns `std::map<string, vector<double>>`
  keyed by the integration variable first.
- Append-row, resize, and **linear interpolation** `interpolate(var, t)`.
- `INTEGRALVALUE(t, 'X')` in the evaluator becomes a thin call to
  `IntegralTable::interpolate` — but only meaningful *during* a step
  (the table is being filled). Outside an integral context it raises a
  diagnostic.

### Phase 3 — Parser & AST hooks (minimal, additive)

**Goal:** let an `.eescode` file containing `INTEGRAL(...)` and
`$IntegralTable` parse without warnings, and capture the directive
content for the IR layer.

| File | Change | Type |
|------|--------|------|
| `src/parser.cpp:400-418` (`knownBuiltinFunctions`) | Add `integral` and `integralvalue`. | additive |
| `src/parser.cpp:908-945` (`tryParseDirective`) | When `name == "integraltable"`, additionally store the parsed content on a small new `IntegralTableSpec` (integration var, optional output interval, variable list incl. `X[1..5]` ranges). Keep storing the raw `Directive` so existing behaviour is unchanged. | additive |
| `src/parser.cpp` (`parseArrayAccess`) | Recognise `..` range notation so `$IntegralTable X[1..5]` lists expand. Only active inside the directive mini-parser; existing equation parsing untouched. | additive |
| `include/coolsolve/ast.h` | Add `IntegralTableSpec` struct carried alongside the `Directive` (no new Statement variant needed — keep it as metadata on the existing `Directive`). | additive |

**Tests:** add `tests/test_parser.cpp` cases that the directive and the
function call parse without warnings and that the spec is captured.

> **Note on `$IntegralAutoStep` / `$IntegralStop`:** these are *detected
> and warned* ("use coolsolve.conf keys `integral*` instead") but **not
> interpreted**. No code path depends on them. This matches the user's
> preference and keeps a future opt-in possible.

### Phase 4 — IR extraction of the `IntegralProblem`

**Files:** `src/integral/integral_extraction.cpp` + header,
`tests/test_integral_extraction.cpp`.

`extractIntegralProblem(ir, analysis)`:

1. Scan all equations for a top-level `INTEGRAL(...)` call on the RHS.
2. For each such equation `lhs = <integral>`, classify `lhs` as a
   **state variable**; record its integrand (1st arg), integration var
   (2nd arg), and limits.
3. Validate: all integral equations share the **same** integration var
   and the **same** `[lo, hi]`. If not → diagnostic + `valid=false`.
   (Multi-interval / multi-variable integration is out of scope; raise a
   clear message.)
4. Partition all variables into **state** (LHS of an integral eqn),
   **algebraic** (everything else that is unknown), and **parameters**
   (set by the user or by initial conditions).
5. Partition blocks: blocks containing a state-variable's derivative are
   "derivative blocks"; the rest are "algebraic blocks".
6. **High-index detection (placeholder):** if a purely algebraic equation
   constrains a state variable without going through its derivative,
   flag it as potential high-index and emit a diagnostic
   `"High-index DAE detected; index reduction is not yet supported"`.
   This is a conservative rejection; do not attempt reduction (§9
   placeholder).

Tests: feed small hand-built IRs (reusing the `IR::fromAST` path from
existing parser tests) and assert the classification.

### Phase 5 — `IntegralSolver` orchestration

**Files:** `src/integral/integral_solver.cpp` + header,
`tests/test_integral_e2e.cpp`.

The time-march loop:

```
extract IntegralProblem
initialise IntegralTable columns
t := t0;  y := y0
record row (t, y, z0)
while t < tf:
    h := stepSize (fixed) or adaptive (RK45)
    f := evaluateRHS(t, y)         // sets externals, solves algebraic blocks
    (y_new, err) := integrator.step(t, y, f, h)
    if RK45 and err > tol:  reject, shrink h, retry   (no row written)
    y := y_new; t := t + h
    z := last algebraic solution
    record row (t, y, z)   honouring the output interval from $IntegralTable
evaluateRHS / algebraic solve for the final (t, y) so that algebraicResult
is consistent and the post-solve solution_check has a complete state.
```

Reuses: `Solver` (unmodified), `SystemEvaluator` (unmodified),
`DiagnosticCollector` (unmodified).

**End-to-end tests** (`tests/test_integral_e2e.cpp`): small `.eescode`
models solved against analytical solutions:

- Exponential decay `dy/dt = -y` ⇒ `y(4) = e^{-4}`.
- Free fall with drag (the EES example in `_integraltable.htm`).
- Coupled linear system (harmonic oscillator).
- A model with an algebraic variable coupled to a state (heat exchanger
  cell) to exercise the algebraic-resolve path.

Each end-to-end test asserts both the final value and selected
`IntegralTable` rows.

### Phase 6 — Wire into `CoolSolveRunner` (single dispatch point)

**File:** `src/runner.cpp` (one new branch in `run()`).

```cpp
if (integral_extraction::hasIntegral(*ir_)) {
    IntegralSolver isolver(*ir_, *analysisResult_, options);
    integralResult_ = isolver.solve(integratorOptionsFrom(options));
    // map IntegralSolveResult -> solveResult_ for downstream code
} else {
    solveResult_ = solver.solve(options, enableTracing);   // unchanged
}
```

`integratorOptionsFrom(SolverOptions)` maps the new `integral*` fields
(§5) into `IntegratorOptions`. This is the **only** edit to `runner.cpp`
and it is additive.

### Phase 7 — Output surfaces (CLI, files, debug)

| Surface | Change | File |
|---------|--------|------|
| **Auto CSV file** | After a successful integral solve, write `<modelname>-integral.csv` next to the `.eescode` (same folder convention as companion lookup-table CSVs `<modelname>-<table>.csv`). First column = integration variable; then one column per `$IntegralTable` variable; header row = variable names. Overwritten on each solve. This is the primary on-disk artefact and the source the GUI re-reads. | new helper `IntegralTable::writeCSV()` in `src/integral/integral_table.cpp`, called from `main.cpp` and from the server solve handler. |
| `.sol` file | After the scalar variables, append an `# IntegralTable` block (one column-header line + CSV rows). Backward compatible: readers that stop at the first unknown line are unaffected. | `main.cpp` |
| CLI JSON | Add `"integralTable": { "t": [...], "y": [...], ... }` and `"integralCsv": "<modelname>-integral.csv"` to the analysis/solution JSON only when an integral was solved. | `src/structural_analysis.cpp` (`generateAnalysisJSON`) or a new helper. |
| Debug folder (`-d`) | New `integral.md` file: problem summary, method, step count, min/max/avg step, a preview of the table; plus add it to the debug `README.md` index and to the README table. A full copy of `<modelname>-integral.csv` is also dropped into the debug folder. | `src/runner.cpp` (`generateDebugOutput`) |
| Solution verification (`solution_checker.cpp`) | Extend to re-check the integral equations at a few sampled table rows (not just the final state). | `src/solution_checker.cpp` |
| LaTeX report | Optional: include the integral equations and a small plot of the first tabulated variable. Degraded gracefully if the plot is absent. | `src/latex_report.cpp` |

> The auto-written CSV mirrors the **parametric-table** UX: a file named
> after the model, regenerated on each successful run, that the GUI
> loads and that travels in the ZIP bundle (Phase 9).

### Phase 8 — `coolsolve.conf` and `SolverOptions` (config layer)

See §5 for the full key list. Edits:

- `include/coolsolve/solver.h` — add the `integral*` fields to
  `SolverOptions`.
- `src/solver.cpp` (`loadSolverOptionsFromFile`) — add `else if` branches
  for each new key (pattern already used for every other option).
- `examples/coolsolve.conf` — add the new keys, commented out, with the
  usual explanatory block.
- `tests/test_config.cpp` — round-trip regression test for each key.

### Phase 9 — Server, REST API, and ZIP bundle round-trip

The GUI exchanges complete model state as a ZIP bundle
(`/api/v1/files/bundle`, `/api/v1/files/upload`) and follows the same
pattern as the parametric study (`session.lastParametricResult`,
`server.cpp:2463-2469`). The integral table is wired in identically.

| Change | File | Detail |
|--------|------|--------|
| Persist the table on the session | `src/server.cpp` | New `Session` field `lastIntegralCSV` (the CSV text) + `lastIntegralResult` (JSON), guarded by a mutex, mirroring `lastParametricResult` (`server.cpp:196`, `2033`, `2466`). |
| Embed in the solve response | `src/server.cpp` (`solveResultToJSON`, ~line 489) | Add `integralTable` (columnar JSON) and `integralCsvName` to the `{type:"done"}` event so the frontend receives the table with every solve — no new endpoint. |
| Auto-write the CSV in server mode | `src/server.cpp` (solve handler) | Write `<modelname>-integral.csv` into the session working dir, exactly as the CLI does in Phase 7. |
| **Export in the ZIP bundle** | `src/server.cpp` (`createZipBundle`, ~line 2463) | Push `<modelname>-integral.csv` (and optionally `integral_table.json`) into the `files` list, right after the `parametric_studies.json` entry. |
| **Import from the ZIP bundle** | `src/server.cpp` (`extractZipBundle`, ~line 2506) | Recognise `*-integral.csv` in the uploaded ZIP and restore it into `session.lastIntegralCSV` so the Integral Table tab repopulates on load — round-trip identical to the parametric study. |
| Clear on new model | `src/server.cpp` (new/open/reset paths, ~line 2526) | Reset `lastIntegralCSV`/`lastIntegralResult` alongside the other session fields. |
| New GET endpoint (optional convenience) | `src/server.cpp` | `GET /api/v1/integral/result` returning the last integral JSON, mirroring `GET /api/v1/parametric/result` (`server.cpp:2033`). |
| Tests | `tests/test_parametric_api.py` pattern | Add a Python end-to-end test: solve an integral model via `/api/v1/solve`, assert `integralTable` is present, download the bundle, re-upload it, assert the table is restored. |

**Round-trip contract:** download bundle → upload bundle → the Integral
Table tab shows the same columns and rows. This is the same guarantee
already enforced for the parametric study (`contributing.md §3.7`).

### Phase 10 — GUI Integral Table tab (bottom panel)

The bottom panel already hosts the parametric study as the `sensitivity`
tab (`gui/src/App.tsx:141-166`). The integral table becomes a sibling
tab. This is a **full feature phase, not future work** — it ships with
the rest.

| Change | File | Detail |
|--------|------|--------|
| New tab in the bottom panel | `gui/src/App.tsx` | Add an `integral` button to the bottom tab bar (alongside `console` / `sensitivity` / `lookuptables`) and render `<IntegralTable/>` when active. Follow the exact toggle pattern of the existing tabs. |
| New component | `gui/src/components/IntegralTable.tsx` | Mirror `ParametricStudy.tsx`: a scrollable columnar table (integration variable first), a "Export CSV" button, a row count, and an empty-state message when no integral has been solved. Read-only (EES Integral Tables are not editable). Optional: a `PlotlyChart` line plot of the first tabulated variable vs the integration variable. |
| Types | `gui/src/api/types.ts` | Add `IntegralTableData { columns: string[]; rows: number[][]; integrationVar: string; csvName: string }` and extend `SolveResponse` with an optional `integralTable?: IntegralTableData`. |
| Store | `gui/src/stores/modelStore.ts` | Add `integralTable: IntegralTableData \| null`, set it from the solve response, clear it on new/open/reset (same lifecycle as `parametricStudies`, `modelStore.ts:111-114`). Persist/restore it as part of the bundle state. |
| Bundle round-trip (frontend) | `gui/src/api/client.ts` | The CSV is part of the ZIP; no separate call needed. On bundle load, populate the store from the restored session. |
| Syntax highlighting | `gui/src/languages/ees.ts` | Add `integral`, `integralvalue` to built-in functions and `$IntegralTable` to directives so Monaco highlights them. |
| Config editor | `gui/src/components/ConfigEditor.tsx` | Add a new `ConfigGroup` "Integration" exposing the `integral*` keys from §5 (method dropdown, step, tolerances, Richardson toggle), following the existing `ConfigField` pattern. Update `PIPELINE_PRESETS` if relevant. |
| Smoke test | manual | `npm run dev`, solve an integral example, verify the tab populates, export the bundle, re-import, verify repopulation. |

### Phase 11 — Documentation

| Document | Change |
|----------|--------|
| `README.md` | New bullet in *Features*; new row(s) in *Command Line Options* if any; update *Project Structure* tree with `src/integral/`; add `integral.md` to the *Debug folder* table; add this plan (or a user-facing spin-off) to the *Documentation* table. |
| `docs/language_reference.md` | New section "Equation-based integration" covering `INTEGRAL`, `$IntegralTable`, `INTEGRALVALUE`, and the `coolsolve.conf` keys. |
| `docs/debugging_models.md` | New subsection on diagnosing integration failures (step rejected, high-index). |
| `docs/solver_roadmap.md` | Mark dynamic-solving as delivered; note BDF/stiff and index reduction as future. |
| `docs/gui.md` | Document the new Integral Table bottom-panel tab and the bundle round-trip. |
| `docs/docs.html` | Append the new doc page to the sidebar nav. |
| **`docs/ees_vs_coolsolve.csv`** | **Update** the rows for `INTEGRAL`, `INTEGRALVALUE`, `$IntegralTable`, `$IntegralAutoStep`, `$IntegralStop`, and add a "Solver Features → Dynamic/DAE solving" row, flipping `Implemented` to `Yes`/`Partial` where appropriate and filling the `Notes` column with the file references. |
| `docs/versions.md` | Changelog entry at release time. |

---

## 5. New `SolverOptions` fields (`coolsolve.conf` keys)

All default to "off / inert" so that a non-integral model behaves
exactly as today (`contributing.md §1.7`, `§5.2`).

| `coolsolve.conf` key | C++ field | Type | Default | Meaning |
|----------------------|-----------|------|---------|---------|
| `integralMethod` | `integralMethod` | enum `{EulerExplicit, EulerImplicit, RK4, RK45}` | `RK4` | Fixed-step default; RK45 = variable. |
| `integralFixedStep` | `integralFixedStep` | double | `0.0` | `0` ⇒ derive from `integralMaxSteps`. |
| `integralMaxSteps` | `integralMaxSteps` | int | `1000` | Upper bound on number of steps. |
| `integralRelTol` | `integralRelTol` | double | `1e-6` | RK45 local error control. |
| `integralAbsTol` | `integralAbsTol` | double | `1e-9` | RK45 absolute error floor. |
| `integralMinStep` | `integralMinStep` | double | `0.0` | Min step (0 = auto). |
| `integralMaxStep` | `integralMaxStep` | double | `0.0` | Max step (0 = auto). |
| `integralRichardson` | `integralRichardson` | bool | `false` | Richardson extrapolation (fixed step only). |
| `integralOutputInterval` | `integralOutputInterval` | double | `0.0` | Default output interval when `$IntegralTable` omits the `:n`. `0` = every step. |

---

## 6. Testing strategy (TDD)

Following `contributing.md §2.3`:

1. **Phase 0–1:** unit tests for each `Integrator` against analytical
   solutions. Run with:
   ```bash
   ./coolsolve_tests "[integrator]"
   ```
2. **Phase 2–4:** unit tests for `IntegralTable` and
   `extractIntegralProblem`.
3. **Phase 5:** end-to-end tests with `.eescode` snippets.
4. **Phase 6–8:** regression: run the full suite to prove zero impact on
   existing models.
   ```bash
   ./coolsolve_tests
   ./coolsolve_tests "[examples-comprehensive]"
   ./coolsolve_tests "[solver-robustness]"
   ```
   Inspect `examples/test_examples.md` and
   `examples/solver_robustness_report.md` — iteration counts must be
   unchanged on non-integral models (`contributing.md §5.4`).
5. **Negative tests:** high-index system rejected with the right message;
   nested `INTEGRAL` rejected; table-based 2-arg `INTEGRAL` rejected;
   `$IntegralAutoStep` warned-and-ignored.

Add example `.eescode` files to `examples/` (e.g.
`examples/integral_decay.eescode`,
`examples/integral_free_fall.eescode`) with entries in
`EXPECTED_SOLUTIONS` (`tests/test_examples.cpp`) per `contributing.md
§3.11`.

---

## 7. Performance discipline (`contributing.md §5`)

- The `hasIntegral(IR)` check is a single linear scan; it short-circuits
  to the existing solver for every non-integral model — **zero overhead
  by default**.
- The time-march loop reuses the algebraic `Solver`'s warm-start and
  caching (CoolProp `AbstractState` cache is thread-local and persists
  across steps). No new allocations in the inner algebraic path.
- RK45 step rejection is the only genuinely new cost; it is bounded by
  `integralMaxSteps` and reported in the debug folder.
- All new options default to inert; no existing model pays anything.

---

## 8. Mapping to the contributing checklist (`contributing.md §3`)

| Item | Applies? | Where in this plan |
|------|----------|--------------------|
| §3.1 Language / parser / AST / IR / AD | Yes | Phase 3 (parser), Phase 4 (IR) |
| §3.2 Solver pipeline / SolverOptions | Yes (additive) | Phase 5, Phase 8 |
| §3.3 `examples/coolsolve.conf` + `test_config` | Yes | Phase 8 |
| §3.4 GUI ConfigEditor | Yes | Phase 10 ("Integration" config group) |
| §3.5 GUI components, API types & client, stores, EES Monaco language | Yes | Phase 10 (`IntegralTable.tsx`, store, types, `ees.ts`) |
| §3.6 REST endpoints in src/server.cpp | Yes | Phase 9 (table rides on `/solve`; optional `/integral/result`) |
| §3.7 ZIP bundle round-trip (createZipBundle / extractZipBundle) | Yes | Phase 9 (`<modelname>-integral.csv` exported + restored) |
| §3.8 Debug folder | Yes | Phase 7 (`integral.md` + CSV copy) |
| §3.9 stdout/stderr/diagnostics | Yes | Phase 4–5 (diagnostics) |
| §3.10 LaTeX report | Optional | Phase 7 |
| §3.11 Example `.eescode` + `EXPECTED_SOLUTIONS` | Yes | Phase 5–6 |
| §3.12 README + language_reference + docs | Yes | Phase 11 |
| §3.13 CMake | **No change** (auto-globbed) | — |

---

## 9. Future work (placeholders left in code)

1. **Index reduction** (Pantelides 1988; dummy-derivative index
   reduction). Commented hook in `extractIntegralProblem` where the
   high-index check currently rejects. This is the single biggest
   enabler for fully-implicit DAE models.
2. **BDF / DASSL-style stiff variable-step solver.** The
   `Integrator` base class and the `createIntegrator` factory make this
   a drop-in addition; it would reuse the same algebraic-solve
   infrastructure.
3. **Multi-variable (nested) integration.** Detection is already
   specified (reject with a clear message); a future implementation would
   nest two `IntegralSolver` instances.
4. **`$IntegralStop`** as a runtime stop condition (could be re-expressed
   as a `coolsolve.conf` key `integralStopCondition = "var<value"`).
5. **Event detection** (zero-crossing) inside the RK45 step for precise
   stop locations.
6. **Plotting in the GUI tab** beyond the single-variable quick plot
   (multi-axis, P-h/T-s overlay, etc.) — the basic table + CSV export
   ships in Phase 10.

---

## 10. Open items to confirm before coding starts

- Final naming for the `coolsolve.conf` keys (§5) — consistent with the
  existing `bisectionND*`, `tearing*`, `multiStart*` prefix convention.
- Whether the `.sol` IntegralTable block should be opt-in (only written
  with a flag) to keep `.sol` files small for non-integral users
  (default: write only when an integral was actually solved, so it is
  naturally opt-in).
- Default `integralMethod`: this plan proposes `RK4`. Confirm.

---

## 12. Progress Log

A running record of what has been implemented and verified. Each entry is
added after the phase's targeted tests pass.

### 2026-07-05 — Phase 0 + Phase 1 (integrators) ✅

**Delivered:**

- `include/coolsolve/integral/integrator.h` — public interface:
  `IntegratorOptions` (method enum + tolerances), `RHSFunction`,
  `StepResult` (carries `yNew`, `stepTaken`, `nextStep`, `errorEstimate`,
  `accepted`, `rhsEvaluations`), abstract `Integrator` (with `order()` for
  Richardson), `createIntegrator()` factory, `wrapRichardson()`,
  `methodToString()`, `parseIntegralMethod()`.
- `src/integral/integrator.cpp` — factory + helpers.
- `src/integral/integrator_euler_explicit.cpp` — forward Euler (order 1).
- `src/integral/integrator_euler_implicit.cpp` — backward Euler with an
  internal Newton + finite-difference Jacobian (A-stable, order 1).
- `src/integral/integrator_rk4.cpp` — classic RK4 (order 4), 4 evals/step.
- `src/integral/integrator_rk45.cpp` — Dormand-Prince DOPRI5 embedded pair
  with the Hairer §II.4 step-size controller (order 5, 7 stages, adaptive).
- `src/integral/richardson.cpp` — Richardson wrapper using the **general**
  combination `(2^p·I_{h/2} − I_h)/(2^p − 1)` driven by the base method's
  reported order. (The EES-doc formula `(4·I_{h/2} − I_h)/3` is the p=2
  special case; the original plan text used it generically, which is only
  correct for order-2 methods — corrected here.)
- `src/integral/integrators_internal.h` — internal `make*()` declarations.
- `tests/integrator_test_util.h` — shared `marchFixed` / `marchAdaptive`
  helpers (the outer time-march loop that becomes `IntegralSolver` in §5).
- `tests/test_integrator_euler.cpp`, `test_integrator_rk4.cpp`,
  `test_integrator_rk45.cpp`, `test_richardson.cpp` — TDD tests against
  analytical solutions (exp decay, t² polynomial, harmonic oscillator,
  stiff `dy/dt = -1000 y`, RK45 tolerance ladder, convergence-order
  verification for every method and for Richardson's order gain).

**Build:** no `CMakeLists.txt` change — the new files are auto-globbed by
the existing `GLOB_RECURSE` (`CMakeLists.txt:142`, `:266`). One `cmake ..`
reconfigure picks them up.

**Test results:** `./coolsolve_tests "[integrator],[richardson]"` →
16 test cases, 38 assertions, all pass. Existing suites (config, newton)
unaffected — the integral module is fully isolated under
`coolsolve::` / `src/integral/`.

**Decisions / deviations from the plan:**

- The `Integrator::step()` signature takes the proposed step `h` as an
  explicit argument (not bundled inside `IntegratorOptions`), and
  `StepResult::nextStep` carries the adaptive controller's recommendation
  back to the caller. This keeps the step-size logic inside RK45 (where it
  belongs) and leaves the outer march loop trivial — it will map cleanly
  onto `IntegralSolver` in Phase 5.
- `Integrator::order()` was added so Richardson picks the correct 2^p factor
  per method (plan §2.3 only quoted the p=2 formula).
- DOPRI5 is implemented stateless (7 evals/step) rather than FSAL-cached;
  one extra RHS eval/step is negligible next to CoolProp algebraic solves.

**Next:** Phase 2 (`IntegralTable` + `INTEGRALVALUE`).

### 2026-07-05 — Phase 2 (IntegralTable + interpolation) ✅

**Delivered:**

- `include/coolsolve/integral/integral_table.h` — `IntegralTableSpec`
  (integration var, output interval, expanded column list) and the
  `IntegralTable` class (columnar `std::map<string, vector<double>>`
  storage, deterministic column order, integration var always column 0).
- `src/integral/integral_table.cpp` — `appendRow` (name-map and ordered-row
  overloads), `value`, `column`, **`interpolate`** (binary-search linear
  interpolation with clamping at the endpoints), `clear`, `toCSV`,
  `writeCSV`.
- `tests/test_integral_table.cpp` — TDD coverage: column setup, append,
  missing-column NaN, linear-interpolation midpoints + clamping,
  empty/single-row/unknown-column edge cases, CSV round-trip,
  `IntegralTableSpec::isPresent`.

**Test results:** `./coolsolve_tests "[integral-table]"` → 6 test cases,
23 assertions, all pass. Existing suites unaffected.

**Decision / deviation from the plan:** the `INTEGRALVALUE(t,'X')` evaluator
dispatch is deferred to Phase 5. It needs both (a) the active-table context
that only `IntegralSolver` can supply during a step, and (b) parser
recognition of `integralvalue` (Phase 3). The interpolation numerical core
it calls — `IntegralTable::interpolate` — is fully implemented and tested
here, so Phase 5 just wires the evaluator name to it.

**Next:** Phase 3 (parser & AST hooks).

### 2026-07-05 — Phase 3 (parser & AST hooks) ✅

**Delivered (all changes additive):**

- `include/coolsolve/ast.h` — `Directive` carries an optional
  `IntegralTableSpec` payload (`hasIntegralTableSpec` flag). Reuses the
  struct defined in `integral_table.h`; no new Statement variant.
- `src/parser.cpp`
  - `knownBuiltinFunctions()` — added `integral` and `integralvalue` so the
    calls parse without "Unknown function" warnings.
  - `knownDirectives` — added `integralautostep` and `integralstop`
    (recognised for compatibility).
  - new `parseIntegralTableContent()` — parses `$IntegralTable t:0.1 y X[1..5]`
    into an `IntegralTableSpec`, expanding `X[lo..hi]` ranges into a flat
    column list at parse time.
  - `tryParseDirective()` — attaches the spec to the `IntegralTable` directive,
    warns on empty/invalid specs (P006), and emits a dedicated diagnostic
    (P005) for `$IntegralAutoStep`/`$IntegralStop` directing users to the
    `integral*` `coolsolve.conf` keys.
- `tests/test_parser_integral.cpp` — 7 TDD cases: spec capture, interval
  parsing, `X[1..5]` expansion, default-zero interval, `INTEGRAL`/`integralvalue`
  parse without warnings, `$IntegralAutoStep`/`$IntegralStop` warned-and-ignored,
  empty `$IntegralTable` negative case.

**Test results:** `./coolsolve_tests "[parser],[integral]"` → all pass;
`./coolsolve_tests "[parser]"` → 22 test cases, 154 assertions, no regressions.

**Decisions:**

- `IntegralTableSpec` lives canonically in `integral_table.h`; `ast.h`
  includes that header to carry the optional payload. The dependency is
  one-way (integral headers never include `ast.h`), so no layering cycle.
- `makeDirective(name, content)` still works unchanged: the new aggregate
  members are value-initialised (empty spec, `hasIntegralTableSpec=false`).
- `$IntegralTable` columns keep the bracket form `X[1]` as written; the
  IR variable-name reconciliation (e.g. flattened `X_1`) happens in Phase 4.

**Next:** Phase 4 (IR extraction of the `IntegralProblem`).

### 2026-07-05 — Phase 4 (IR extraction of IntegralProblem) ✅

**Delivered:**

- `include/coolsolve/integral/integral_problem.h` — `StateVariable` (name,
  integrand var/expr, base expression, equation id) and `IntegralProblem`
  (integration var, limit expressions + constant-folded values, optional
  fixed step, state/algebraic variable classification, equation-id
  partition, `IntegralTableSpec`, validity + diagnostics). Declares
  `hasIntegral()` and `extractIntegralProblem()`.
- `src/integral/integral_extraction.cpp` — walks every equation for
  top-level `INTEGRAL(...)` calls (4- or 5-arg form), validates a single
  shared integration variable + `[lo, hi]` interval across all states,
  rejects nested integrals (multi-variable integration), classifies state
  vs algebraic variables (excluding the integration var, which the
  integrator owns), partitions equations, and runs a conservative
  high-index *warning* plus a structural squareness check on the algebraic
  subsystem.
- `tests/test_integral_extraction.cpp` — 12 TDD cases built through the
  real parser→IR→analyze path: decay, base expression, two-state oscillator,
  fixed step, inconsistent var/limits (rejected), nested integrals
  (rejected), wrong arg count, `hasIntegral` true/false, variable
  classification, non-square subsystem flagging.

**Key design decision — initial values:** the initial state is *not*
extracted here. At `t = t0` the integral term has a zero-width interval and
evaluates to 0, so `y(t0)` falls out of the algebraic solve at the first
step (handled by `IntegralSolver` in Phase 5). The non-integral part of the
RHS (`baseExpr`, computed by blanking the INTEGRAL call to 0) is what
recovers `y(t0)`.

**Deviation from the plan — high-index handling:** the plan described
high-index detection as a *rejection*. That would false-positive on the
feature's primary use case: legitimate index-1 thermo models routinely
have state variables (e.g. temperatures) appearing in algebraic equations
(heat-transfer laws). It is therefore emitted as a **warning diagnostic**
(`prob.diagnostics`), not `valid=false`. Hard structural rejection is left
to the algebraic-subsystem squareness check. True Pantelides-style index
reduction remains a §9 placeholder.

**Test results:** `./coolsolve_tests "[extraction]"` → 12 cases, 48
assertions, all pass. All 41 integral-module tests pass; parser/IR/analysis
suites (29 cases, 197 assertions) show no regressions.

**Next:** Phase 5 (`IntegralSolver` orchestration + end-to-end tests).

### 2026-07-05 — Phase 5 (IntegralSolver orchestration) ✅

**Delivered:**

- `include/coolsolve/integral/integral_solver.h` — `IntegralSolveResult`
  (success, problem, table, step counts, accepted step sizes, final algebraic
  result) and the `IntegralSolver` class.
- `src/integral/integral_solver.cpp` — the time-march loop. Construction:
  harvests the `$IntegralTable` spec from the AST, builds a **reduced
  algebraic IR** (integral equations removed + driver equations added),
  re-analyses it, and constructs a long-lived algebraic `Solver` reused at
  every step. `solve()`: computes the initial state at `t0` (via each
  state's `baseExpr`), then marches with the chosen `Integrator`, recording
  rows into an `IntegralTable` (honouring the output interval).
- `tests/test_integral_e2e.cpp` — 6 end-to-end cases: exponential decay,
  coupled harmonic oscillator, RK45 adaptive, an algebraic variable coupled
  to a state, fixed step from the 5th `INTEGRAL` argument, and
  `$IntegralTable` column/interval honking.

**The key engineering problem and its solution.** Two non-obvious obstacles
had to be overcome to reuse the algebraic `Solver` per step:

1. **The full IR is non-square.** An integral model is *deliberately*
   under-determined algebraically: the integration variable `t` is a free
   parameter and the states are owned by the integrator. But
   `StructuralAnalyzer::analyze()` rejects non-square systems up front
   (`structural_analysis.cpp:532`). Solution: the `IntegralSolver` never asks
   the analyzer to solve the full model — it builds its own **reduced IR**.
2. **The reduced IR is *also* non-square** (states + `t` have no defining
   equation there). Solution: for every integrator-owned variable the
   constructor adds a **driver equation** `v = <NumberLiteral>`. Each step
   the integrator mutates that literal's `value` to the current `y`/`t`
   (exploiting that `IR::fromAST` *shares* the AST `ExprPtr`s, `ir.cpp:203`),
   so the explicit driver block resolves to the integrator's value and every
   other block sees it as an external. This keeps the reduced IR square and
   pins the integrator-owned variables without touching the algebraic solver.

**Decisions / deviations:**

- `IntegralSolver` stores `SolverOptions` **by value** (the plan showed a
  reference); the caller (runner/tests) may pass a temporary, so a copy is
  the safe choice.
- The `IntegralTableSpec` is harvested from the AST inside the constructor
  (the plan had the runner do it); this keeps Phase 6's runner change minimal
  and makes the solver self-contained.
- Implicit-Euler coupling of state + algebraic unknowns into one Newton
  solve (plan §3.4) is **deferred** — the current `EulerImplicit` integrator
  solves the implicit step with its own internal Newton on the RHS, which is
  correct but re-solves the algebraic subsystem inside each Newton iteration.
  This is a performance item, not a correctness gap; all methods pass their
  analytical end-to-end tests.
- `INTEGRALVALUE(t,'X')` evaluator dispatch remains deferred (it needs the
  active-table context, which now exists; wiring is a small follow-up).

**Test results:** `./coolsolve_tests "[e2e]"` → 6 cases, 27 assertions, all
pass. All 62 integral-module tests pass; no regressions in parser/IR/analysis.

**Next:** Phase 6 (wire the `IntegralSolver` into `CoolSolveRunner`).

### 2026-07-05 — Phase 6 (runner dispatch) ✅

**Delivered:**

- `include/coolsolve/runner.h` — `IntegralSolveResult integralResult_` member,
  an `integralModel_` flag, and accessors `hasIntegralResult()` /
  `getIntegralResult()`. Includes `integral_solver.h`.
- `src/runner.cpp` — `run()` now branches on `hasIntegral(*ir_)`:
  - the structural-analysis step no longer hard-aborts on `success=false`
    (integral models are intentionally non-square);
  - integral models dispatch to `IntegralSolver` (program + IR + analysis +
    options), and the `IntegralSolveResult` is mapped onto `solveResult_`
    (variables, status, error) so all downstream debug/JSON/.sol code keeps
    working unchanged;
  - the algebraic path is byte-for-byte unchanged.
- A file-local `makeIntegratorOptions(SolverOptions)` helper returns sensible
  defaults (RK4, 1000-step budget) — Phase 8 replaces it with the config-driven
  mapping.
- `main.cpp` — the "Structural Analysis Error" early-exit is skipped for
  integral models (`!analysisResult.success && !runner.hasIntegralResult()`).
- `src/solution_checker.cpp` — integral-declaring equations are *skipped*
  during post-solve verification (they are ODE declarations whose correctness
  is established by the trajectory, and the evaluator has no standalone
  `integral` built-in). Phase 7 adds trajectory-based re-checking.
- `examples/integral_decay.eescode` + an `EXPECTED_SOLUTIONS` entry
  (`y = 0.01832 = e⁻⁴`).

**Verification (CLI):**
```
$ ./coolsolve examples/integral_decay.eescode
... Solver: SUCCESS ...
```
The debug folder's `solution_check.md` shows `y(4) = 1.83156e-02 = e⁻⁴` ✓.

**Test results:** all 351 unit tests pass (3096 assertions);
`[examples-comprehensive]` passes (including the new integral example).
No regressions on the algebraic path (the dispatch is guarded by a single
`hasIntegral()` scan that is false for every existing model — zero overhead).

**Next:** Phase 7 (output surfaces: auto-CSV, `.sol`, JSON, debug, solution check).

### 2026-07-05 — Phase 7 (output surfaces) ✅

**Delivered:**

- **Auto CSV** (`main.cpp`): after a successful integral solve, writes
  `<modelname>-integral.csv` next to the `.eescode` (first column = integration
  variable, then the `$IntegralTable` columns). Added `*-integral.csv` to
  `.gitignore` so the artefact is never committed.
- **`.sol` block** (`main.cpp`): appends a `# IntegralTable` section (the CSV)
  after the scalar variables — backward compatible.
- **CLI JSON** (`main.cpp`): injects `"integralTable"` (columnar arrays) and
  `"integralCsv"` into the analysis JSON for dynamic models. Done by re-parsing
  the dumped JSON locally to avoid a circular include between `structural_analysis`
  and the integral module.
- **Debug folder** (`src/runner.cpp::generateDebugOutput`): new `integral.md`
  (problem summary, method, step count, min/avg/max step size, trajectory
  preview, diagnostics) plus a full `integral_table.csv` copy; both registered
  in the debug `README.md` index.
- **README** public debug-folder table: `integral.md` and `integral_table.csv`
  rows added.

**Verification:**
```
$ ./coolsolve examples/integral_decay.eescode -o out.json
# out.json: integralTable.y[-1] = 0.01832 = e^-4, integralCsv = "integral_decay-integral.csv"
# integral_decay-integral.csv: t,y,dydt with y(4)=0.0183156
# -d folder: integral.md + integral_table.csv
# solution_check: ALL EQUATIONS SATISFIED (1 skipped = the ODE declaration)
```

**Deferred items:** the trajectory-based re-checking of integral equations at
sampled rows (plan §7) is deferred — the integral equations are *skipped* in
the post-solve checker (Phase 6) since the trajectory itself is the evidence
of correctness, and the solve is independently validated against analytical
solutions in the Phase 5 e2e tests. A finite-difference slope-vs-derivative
check can be added later as an extra safety net.

**Test results:** all 351 unit tests pass; `[examples-comprehensive]` passes.

**Next:** Phase 8 (`coolsolve.conf` + `SolverOptions` + `test_config`).

### 2026-07-05 — Phase 8 (configuration layer) ✅

**Delivered:**

- `include/coolsolve/solver.h` — nine new `SolverOptions` fields
  (`integralMethod`, `integralFixedStep`, `integralMaxSteps`, `integralRelTol`,
  `integralAbsTol`, `integralMinStep`, `integralMaxStep`, `integralRichardson`,
  `integralOutputInterval`), all inert by default. `integralMethod` is stored as
  a (lower-cased) string so `solver.h` stays decoupled from the integrator
  module's enum.
- `src/solver.cpp` (`loadSolverOptionsFromFile`) — parses all nine keys
  (`std::stod`/`std::stoi`/`parseBool`, case-insensitive method string).
- `src/runner.cpp` — `makeIntegratorOptions()` now maps the `SolverOptions`
  `integral*` fields into `IntegratorOptions` (via `parseIntegralMethod`),
  replacing the placeholder defaults.
- `examples/coolsolve.conf` — all nine keys added, commented out, with an
  explanatory block (mirroring the existing option-documentation style).
- `tests/test_config.cpp` — a round-trip regression test ("Integral options are
  loaded from config") writes every key, reloads it, and asserts each field,
  plus that the method string parses back to `IntegratorOptions::RK45`.

**Verification (config-driven CLI):** placing
```
integralMethod = RK45
integralRelTol = 1e-7
```
in `coolsolve.conf` next to `integral_decay.eescode` solves the model and
produces `y(4) = 0.0183156 = e⁻⁴`.

**Test results:** all 352 unit tests pass (3109 assertions); the new config
test passes; `[examples-comprehensive]` unaffected.

**Next:** Phase 9 (server / REST API / ZIP bundle round-trip).

### 2026-07-05 — Phase 9 (server / REST API / ZIP round-trip) ✅

**Delivered (all in `src/server.cpp`, additive):**

- **Session state** — new fields `lastIntegralResult` (columnar JSON),
  `lastIntegralCSV` (CSV text), `lastIntegralCsvName`, guarded by a new
  `integralMutex` (mirrors `parametricMutex`).
- **Solve handler** — after `runner.run()`, if `runner.hasIntegralResult()`,
  builds the integral JSON (`integralResultToJSON` helper), stores it on the
  session, auto-writes `<model>-integral.csv` into the session temp dir, and
  embeds `integralTable` + `integralCsvName` into both the SSE `"done"` event
  and the `GET /api/v1/solve/result` payload (so the GUI receives the table
  with every solve — no separate fetch needed).
- **ZIP export** (`GET /api/v1/files/bundle`) — pushes `<stem>-integral.csv`
  right after the parametric artefact.
- **ZIP import** (`POST /api/v1/files/upload`) — a new `endsWith("-integral.csv")`
  branch restores the CSV into `lastIntegralCSV` **before** the generic
  lookup-table CSV branch (so it is not mis-filed as a lookup table).
- **Reset points** — integral state cleared on new-model, open-file,
  per-solve, and upload (all four `session.hasResult = false` sites), so a
  stale trajectory never lingers.
- **New endpoint** — `GET /api/v1/integral/result` returns the last integral
  JSON (404 when none), mirroring `/api/v1/parametric/result`.
- **Python test** — `tests/test_integral_api.py` (12 checks, `--auto`
  self-starting): solve → table in result → `y(4)=e⁻⁴` → bundle contains
  `*-integral.csv` → `/new` clears → re-upload restores the CSV (round-trip).

**Verification:**
```
$ python3 tests/test_integral_api.py --auto
Results: 12/12 passed, 0 failed
```

**Design note:** the JSON `lastIntegralResult` is produced only by a live
solve; on ZIP re-upload only the CSV is restored (parsing CSV→columnar JSON on
import was judged not worth the complexity). The GUI tab (Phase 10) will parse
the CSV directly. The round-trip contract — *download → upload → same columns
and rows* — holds via the CSV.

**Test results:** all 352 C++ unit tests pass (3109 assertions); the new
Python API test passes 12/12.

**Next:** Phase 10 (GUI Integral Table tab).

---

## 13. Remaining work (Phase 11)

The **full stack — backend, REST API, ZIP round-trip, and the GUI tab — is
feature-complete and tested.** An integral model can be parsed, solved,
tabulated, exported to CSV, embedded in the solve response, round-tripped
through a ZIP bundle, and visualised in the bottom-panel Integral tab. What
remains is the **user-facing documentation** (Phase 11).

### 13.1 Phase 10 — delivered (see §12 progress log)

Phase 10 (GUI Integral Table tab) shipped on 2026-07-09. See the
**2026-07-09 — Phase 10 (GUI Integral Table tab)** entry in §12 for the full
delivery record. The frontend can now solve, display, plot, export, and
round-trip an integral trajectory through the bottom-panel Integral tab.

### 13.2 Phase 11 — Documentation

| Document | Change |
|----------|--------|
| `README.md` | New *Features* bullet for equation-based integration; add `src/integral/` to the *Project Structure* tree; add the `integral*` keys to the *coolsolve.conf* options list (already in `examples/coolsolve.conf`). |
| `docs/language_reference.md` | New section "Equation-based integration" covering `INTEGRAL`, `$IntegralTable`, `INTEGRALVALUE`, and the `coolsolve.conf` keys. |
| `docs/debugging_models.md` | New subsection on diagnosing integration failures (step rejected, high-index warning, non-constant limits). |
| `docs/solver_roadmap.md` | Mark dynamic solving as delivered; note BDF/stiff and index reduction as future. |
| `docs/gui.md` | Document the Integral Table bottom-panel tab and the CSV bundle round-trip. |
| `docs/docs.html` | Append any new doc page to the sidebar nav. |
| `docs/ees_vs_coolsolve.csv` | Flip `INTEGRAL` / `INTEGRALVALUE` / `$IntegralTable` to `Yes`; add a "Dynamic/DAE solving" row. |
| `docs/versions.md` | Changelog entry at release time. |

> Note: `INTEGRALVALUE(t,'X')` evaluator dispatch is the one functional item
> still deferred from the original plan (it needs the active-table context,
> which now exists in `IntegralSolver`; wiring is a small follow-up in
> `src/evaluator.cpp` + a thread-local pointer to the current `IntegralTable`).
> It is not exercised by any current example or test and does not block the
> core dynamic-solving workflow.

### 2026-07-09 — Phase 10 (GUI Integral Table tab) ✅

**Delivered** — the bottom-panel Integral tab ships with the same UX as the
parametric study: a scrollable columnar table (integration variable first), a
row/step count, an "Export CSV" button, an optional `PlotlyChart` line plot
of one tabulated variable versus the integration variable, and an empty-state
message. Read-only (EES Integral Tables are not user-editable).

**Frontend changes (all additive):**

| File | Change |
|------|--------|
| `gui/src/components/IntegralTable.tsx` (**new**) | The tab component. Resolves the active table from two sources in priority order: (1) the live columnar JSON (`integralTable`), (2) the restored CSV text (`integralCSV`) parsed client-side via a small CSV splitter that honours quoting. Renders the columnar table (integration var highlighted with the existing `.sweep-col` style), a Y-variable selector + Plotly line plot, and an "Export CSV" button that re-serialises the table to a downloadable `<csvName>`. |
| `gui/src/api/types.ts` | New `IntegralTableData` interface (mirrors the backend `integralResultToJSON` payload: `integrationVar`, `columns`, `data` keyed by column, `numRows`, `csvName`, optional `totalSteps`/`rejectedSteps`). `SolveResponse` extended with optional `integralTable?` and `integralCsvName?`. |
| `gui/src/api/client.ts` | `getIntegralResult()` (columnar JSON via `GET /integral/result`) and `getIntegralCSV()` (raw CSV text via `GET /integral/csv`, returning `null` on 404 so the round-trip degrades gracefully when there is no table). |
| `gui/src/stores/modelStore.ts` | New `integralTable: IntegralTableData \| null` + `integralCSV: string` state, with `setIntegralTable`/`setIntegralCSV` actions. Both are reset on `loadFile` and `clearModel` (same lifecycle as `parametricStudies`). |
| `gui/src/stores/uiStore.ts` | `'integral'` added to the `BottomTab` union. |
| `gui/src/App.tsx` | New `Integral` button in the bottom-panel tab bar (alongside Console/Parametric/Lookup Tables) following the exact toggle-on-reclick pattern; renders `<IntegralTable/>` when active. |
| `gui/src/components/Toolbar.tsx` | Solve `done` handler now: (a) stores `result.integralTable` when present and clears `integralCSV`, otherwise clears the table and fetches the CSV as a fallback; (b) ensures the bottom panel is open when an integral was produced. Bundle `handleOpen` and `handleBack` both restore the trajectory by fetching `/integral/result` (columnar) and `/integral/csv` (CSV) after the file load, so the tab repopulates on bundle round-trip. |
| `gui/src/components/ConfigEditor.tsx` | New `Integration` config group exposing all nine `integral*` keys from §5 (method, fixed step, max steps, rel/abs tol, min/max step, Richardson, output interval) with the standard `ConfigField` renderer and full tooltips. |
| `gui/src/languages/ees.ts` | `integral` and `integralvalue` added to `builtinFunctions` so Monaco highlights them as predefined. (`$IntegralTable` / `$IntegralAutoStep` / `$IntegralStop` are already matched by the generic `\$\w+` directive rule.) |
| `gui/src/App.css` | New `.integral-table-panel` / `.integral-header` / `.integral-meta` styles (the table reuses `.parametric-table` / `.sweep-col`). |

**Backend change (minimal additive):**

| File | Change |
|------|--------|
| `src/server.cpp` | New `GET /api/v1/integral/csv` endpoint returning the trajectory CSV text (`text/csv`, 404 when none). This closes the bundle round-trip gap: on ZIP re-upload only the CSV is restored server-side (`session.lastIntegralCSV`), so the GUI tab needs a way to fetch it. Mirrors the existing `GET /tables/{name}` lookup-table CSV pattern. |

**Why the extra backend endpoint?** The plan's Phase 10 table said "no
separate call needed" for the bundle round-trip, but that assumed the solve
response's columnar JSON is always available. On a pure bundle load (no live
solve since server restart), `session.lastIntegralResult` (the JSON) is empty
— only `session.lastIntegralCSV` (the CSV) survives the round-trip (Phase 9
design note, line ~1058). A 21-line `GET /integral/csv` endpoint is the
cleanest way to expose that restored CSV to the frontend, keeping the tab
fully functional on bundle load. The frontend parses the CSV client-side
(`parseCSVToTable` in `IntegralTable.tsx`).

**Verification (live smoke test):**

Started the server in `--gui` mode serving the freshly built `gui/dist`,
opened `examples/integral_decay.eescode`, solved, and confirmed:

```
solve.integralTable present: True | cols: ['t', 'y', 'dydt'] | csvName: integral_decay-integral.csv
  y(4) = 0.018315638888891046  (e^-4 = 0.01831563888873418)
GET /integral/csv -> first line: t,y,dydt
GET /integral/csv -> last line: 4,0.0183156,-0.0183156
```

The bundle round-trip (`*-integral.csv` exported → re-uploaded → CSV restored
via `/integral/csv`) is covered by the extended `tests/test_integral_api.py`.

**Test results:**

- `./coolsolve_tests` → **352 test cases, 3109 assertions, all pass**
  (no regressions; the `server.cpp` change is purely additive).
- `python3 tests/test_integral_api.py --auto` → **14/14 passed**
  (was 12/12; the 2 new checks verify `GET /integral/csv` returns CSV text
  with the correct header and the `y(4)=e⁻⁴` row after a bundle round-trip).
- Frontend: `tsc -b` clean; `npm run build` succeeds (1758 modules, 4.95 s);
  `npm run lint` introduces **zero new errors** (22 problems, identical to
  the pre-Phase-10 baseline — all pre-existing `any`/`_removed` warnings).

**Decisions / deviations from the plan:**

- The plan suggested a "method dropdown" for `integralMethod` in the config
  editor. The standard `ConfigField` renderer only dropdowns booleans, so
  `integralMethod` uses a text input with a tooltip listing the four valid
  values (`EulerExplicit`, `EulerImplicit`, `RK4`, `RK45`) — consistent with
  the existing `kinsolGlobalStrategy` field. A dedicated enum-dropdown
  renderer would be a general ConfigEditor improvement, out of scope here.
- The Integral tab does **not** auto-switch on solve (only ensures the bottom
  panel is open), so the user keeps their console focus during a solve. The
  tab badge/population is immediate once the user clicks it.
- CSV parsing on the frontend is a small bespoke splitter (handles quoting)
  rather than pulling in PapaParse/JSZip — keeps the bundle size unchanged
  and is sufficient for the numeric, headered CSV the backend emits.

#### Pre-existing GUI lint baseline (out of scope — flagged for future cleanup)

Phase 10's `npm run lint` reports **22 problems (19 errors, 3 warnings)**,
all of which pre-date this phase (verified by `git stash` against the
Phase-9 baseline). They are **general GUI health issues, not integral-feature
bugs**: the production build (`tsc -b && vite build`) succeeds and runtime
behaviour is unaffected. Line numbers below are as of the 2026-07-09 delivery
and will drift with future edits — re-run `npm run lint` for current
locations.

| # | Rule | Count | Files (line) | Why it occurs | Severity / impact | Suggested fix |
|---|------|-------|--------------|---------------|-------------------|---------------|
| 1 | `@typescript-eslint/no-explicit-any` | 11 | `Toolbar.tsx` (195, 212, 246, 289, 317, 335, 345, 357, 394, 410), `DebugViewer.tsx` (35) | All are `catch (err: any)` blocks that read `err.message`. Predates TS 4.4's `useUnknownInCatchVariables` default; `any` was the historical convenience. | **Low.** eslint-level only; no runtime/correctness impact. | Change to `catch (err: unknown)` + `err instanceof Error ? err.message : String(err)` (or a small `toMsg(err)` helper). Mechanical, ~15-min sweep. |
| 2 | `@typescript-eslint/no-explicit-any` | 6 | `CodeEditor.tsx` (11, 16, 32, 71, 104, 166) | Monaco editor integration typed loosely: the exported `editorInstance` ref, the `toggleBrace/QuoteComment(editor)` params, an `edits: any[]`, a `useRef<any>`, and `(window as any).monaco` (Monaco is loaded via a `<script>` tag, not the typed npm import). | **Low.** Works because Monaco's runtime shape matches; loses editor/IDE autocomplete on these locals. | `import type * as Monaco from 'monaco-editor'` and type the editor params/refs; replace `(window as any).monaco` with a `declare global` augmentation or the `@monaco-editor/react` loader. |
| 3 | `@typescript-eslint/no-explicit-any` | 1 | `exportPlots.ts` (40) | `(Plotly as any).toImage(...)` — the `plotly.js-cartesian-dist-min` type definitions omit the `toImage` imperative method. | **Low.** Single cast; the call is correct at runtime. | Define a small local interface (`interface PlotlyExport { toImage(gd, opts): Promise<string> }`) and cast through it, or add a scoped `// eslint-disable-next-line` with a comment. |
| 4 | `@typescript-eslint/no-unused-vars` | 1 | `modelStore.ts` (146) | The destructure-to-omit pattern `const { [name]: _removed, ...rest } = state.lookupTableCSVs` in `removeLookupTable`. The `_removed` binding is intentionally unused (underscore prefix is the conventional "I know" signal), but the eslint config has no `varsIgnorePattern: "^_"`. | **Trivial.** | Add `"varsIgnorePattern": "^_"` (and `argsIgnorePattern`) to the `no-unused-vars` rule in `eslint.config.js`, or extract an `omitKey(obj, key)` helper. |
| 5 | `react-refresh/only-export-components` | 2 (warn) | `CodeEditor.tsx` (16, 71) | The file exports both the `CodeEditor` component **and** non-component values (`editorInstance`, `toggleBraceComment`, `toggleQuoteComment`). React Fast Refresh requires component-only files for hot swapping. | **Dev-experience only.** Editing `CodeEditor.tsx` triggers a full page reload instead of an HMR swap; no production impact. | Move `editorInstance` + the two comment-toggle helpers into a sibling non-component module (e.g. `editorUtils.ts`) and import them where needed. |
| 6 | `react-hooks/exhaustive-deps` | 1 (warn) | `ParametricStudy.tsx` (464) | The `plotData` `useMemo` reads `sortedResults` but lists only `[activeResult, plotYVar, plotType]` in its deps. | **Latent.** Currently masked because `activeResult` changes whenever `sortedResults` would, but it is a potential stale-plot bug if the derivation ever decouples. | Add `sortedResults` to the dependency array. |

**Severity summary:** none of these affect correctness or the shipped bundle.
The `no-explicit-any` cluster (items 1–3, 18 of the 19 errors) is the bulk of
the noise and is a straightforward typing sweep; item 4 is a one-line eslint
config change; items 5–6 are dev-experience/latent-hygiene warnings. Tackling
them would bring `npm run lint` to a clean exit and is independent of the
integral feature.

**Next:** Phase 11 (user-facing documentation).
