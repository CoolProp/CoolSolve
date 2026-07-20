# Equation-Based Dynamic Solving

CoolSolve can solve **initial-value differential–algebraic equation (DAE)
systems** written in the EES integral form and tabulate the trajectory of
selected variables. This covers the vast majority of transient thermal models
(tank dynamics, heat exchangers, control loops, etc.).

---

## 1. Mathematical model

### 1.1 The EES integral convention

EES expresses a first-order ODE as an algebraic equation:

```ees
y = y0 + INTEGRAL(dydt, t, t0, tf)     "declares a state variable y"
dydt = ...                              "derivative (algebraic in y, t, …)"
y0 = 1                                  "initial condition"
```

This is equivalent to `dy/dt = f(t, y, …)` with `y(t0) = y0`. A model may
contain several such equations — they all share the **same** integration
variable `t` and the **same** interval `[t0, tf]` (coupled ODEs). Any
variable that does **not** appear as the LHS of an integral equation is
*algebraic* — it is solved at every time step by the existing algebraic
`Solver`.

### 1.2 Semi-explicit index-1 DAE

The general form supported in this version is:

```
y' = f(t, y, z)        ← differential equations  (state variables y)
0  = g(t, y, z)        ← algebraic equations     (algebraic variables z)
```

where `z` can be solved uniquely from `g` given `(t, y)` (index-1
condition). This covers most transient thermal models without requiring
index reduction.

**Higher-index systems** are detected and rejected with a clear message;
placeholders are left for future Pantelides-style index reduction.

---

## 2. Algorithm overview

### 2.1 Time-march loop

When CoolSolve detects any top-level `INTEGRAL(...)` call in the model, it
automatically routes the solve to the **IntegralSolver** path — no flag is
needed. Non-integral models take the standard algebraic path (zero overhead).

The integral solve works as follows:

1. **Extract** the `IntegralProblem` from the IR: walk all equations,
   classify state vs. algebraic variables, partition blocks.
2. **Build a reduced IR** (integral equations removed + driver equations
   added) and re-analyse it. This makes the algebraic subsystem square.
3. **Compute the initial state** at `t = t0`: the integral term has
   zero-width interval, so `y(t0)` falls out of the algebraic solve of the
   remaining equations.
4. **March in time**: at each step, fix `(t, y)` as externals, call the
   algebraic `Solver` on the reduced system to evaluate `f`, then advance
   `y` using the chosen `Integrator`.
5. **Record rows** into an `IntegralTable` honouring the output interval
   from `$IntegralTable`.

The key invariant: **`Solver` is called unmodified at each time step**.
The integral layer only decides *when* to call it, *what* to fix as
external, and *how* to advance the state.

### 2.2 Why a reduced IR?

Two non-obvious obstacles had to be overcome to reuse the algebraic
`Solver` per step:

1. **The full IR is non-square.** An integral model is *deliberately*
   under-determined algebraically: the integration variable `t` is a free
   parameter and the states are owned by the integrator. But
   `StructuralAnalyzer::analyze()` rejects non-square systems up front.
   Solution: `IntegralSolver` builds its own **reduced IR**.

2. **The reduced IR is *also* non-square** (states + `t` have no defining
   equation there). Solution: for every integrator-owned variable the
   constructor adds a **driver equation** `v = <NumberLiteral>`. Each step
   the integrator mutates that literal's `value` to the current `y`/`t`
   (exploiting that `IR::fromAST` *shares* the AST `ExprPtr`s), so the
   explicit driver block resolves to the integrator's value and every
   other block sees it as an external. This keeps the reduced IR square
   and pins the integrator-owned variables without touching the algebraic
   solver.

### 2.3 Integrators

| Method | Step formula | Algebraic solves/step | Notes |
|--------|-------------|-----------------------|-------|
| **Euler explicit** | `y_{n+1} = y_n + h·f(t_n, y_n, z_n)` | 1 | Cheapest; 1st order |
| **Euler implicit** | Implicit solve for `(y_{n+1}, z_{n+1})` | 1 coupled Newton | A-stable; 1st order |
| **RK4** | 4 stages, classic Butcher tableau | 4 | 4th order; **default** |
| **RK45** | 6 stages, Dormand–Prince embedded pair (4th + 5th) | 6 | Adaptive; error-estimation step control |
| **Richardson** | Combine one `h` step + two `h/2` steps | triples | Reduces truncation error; fixed-step only |

The general Richardson combination is `(2^p·I_{h/2} − I_h)/(2^p − 1)`
where `p` is the base method's order (the EES formula `(4·I_{h/2} − I_h)/3`
is the `p = 2` special case).

---

## 3. Architecture

### 3.1 File layout

```
include/coolsolve/integral/
├── integral_problem.h      IntegralProblem struct (state vars, derivatives, ICs, interval)
├── integral_table.h        IntegralTable (time-series columns + interpolation)
├── integrator.h            Integrator abstract base + IntegratorOptions
├── integral_solver.h       IntegralSolver (orchestrates the time-march loop)
└── integral_extraction.h   extractIntegralProblem(IR) → IntegralProblem

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
├── test_integral_extraction.cpp     TDD: IR → IntegralProblem classification
├── test_integral_table.cpp          TDD: IntegralTable store/interpolate/toCSV
└── test_integral_e2e.cpp            End-to-end: .eescode → trajectory → CSV

gui/
├── src/components/IntegralTable.tsx     Bottom-panel tab
├── src/api/types.ts                     + IntegralTableData type
└── src/languages/ees.ts                 highlight integral / $IntegralTable
```

### 3.2 Pipeline integration

```
CoolSolveRunner::run()
  parse → IR → infer → analyze
               │
        hasIntegral(IR)?
          yes → IntegralSolver::solve()
          no  → Solver::solve()       (unchanged)
```

The `hasIntegral(IR)` check is a single linear scan that short-circuits
to the existing solver for non-integral models.

### 3.3 RHS evaluation strategy

At time `t` with state `y`, the RHS `f(t, y)` is evaluated by:

1. Write `t` and the components of `y` into the `SystemEvaluator`'s
   current state as **fixed/external** values (same mechanism the parametric
   sweep uses).
2. Call `Solver::solve()` restricted to the algebraic blocks only.
3. Read the integrand variables (`dydt`, ...) from the solution — those are
   `f`.

For **implicit Euler**, steps 1–3 are wrapped inside the nonlinear solve
for `y_{n+1}`: the block solver receives `y_{n+1}` as additional unknowns
appended to the algebraic block, reusing the existing Newton/LM/TrustRegion
machinery.

---

## 4. IntegralTable

The `IntegralTable` stores named columns `std::map<string, vector<double>>`
keyed by the integration variable first. It supports:

- **Append-row** and **resize** operations.
- **Linear interpolation** via binary search with clamping at endpoints.
- **CSV export** (`writeCSV` / `toCSV`) for file output and GUI display.

`INTEGRALVALUE(t, 'X')` in the evaluator is a thin call to
`IntegralTable::interpolate` — meaningful only *during* an integration step.
The parser recognises the function; full evaluator dispatch is a deferred
follow-up.

---

## 5. Configuration

All nine `integral*` configuration keys live in `coolsolve.conf` and
default to "off / inert" so that non-integral models behave exactly as
before:

| Key | Default | Meaning |
|-----|---------|---------|
| `integralMethod` | `RK4` | `RK4`, `RK45`, `EulerExplicit`, or `EulerImplicit` |
| `integralFixedStep` | `0.0` | Fixed step size. `0` ⇒ derive from `integralMaxSteps` |
| `integralMaxSteps` | `1000` | Upper bound on number of steps |
| `integralRelTol` | `1e-6` | RK45 local relative error control |
| `integralAbsTol` | `1e-9` | RK45 absolute error floor |
| `integralMinStep` | `0.0` | Minimum step (0 = auto) |
| `integralMaxStep` | `0.0` | Maximum step (0 = auto) |
| `integralRichardson` | `false` | Richardson extrapolation (fixed-step only) |
| `integralOutputInterval` | `0.0` | Default `$IntegralTable` row interval (`0` = every step) |

The `$IntegralAutoStep` and `$IntegralStop` directives are recognised for
EES compatibility but **not interpreted** — they emit a diagnostic pointing
at the `integral*` config keys instead.

---

## 6. Output

After a successful integral solve, CoolSolve writes the trajectory to:

| Surface | Detail |
|---------|--------|
| **Auto CSV** | `<modelname>-integral.csv` next to the `.eescode` file (first column = integration var, then `$IntegralTable` columns in order). Regenerated on each run. |
| **`.sol` file** | Appended `# IntegralTable` block (backward compatible). |
| **CLI JSON** | `"integralTable"` and `"integralCsv"` fields in the analysis/solution JSON. |
| **Debug folder** | `integral.md` (problem summary, method, step count, step statistics, trajectory preview) + copy of the CSV. |
| **GUI** | **Integral** tab in the bottom panel: scrollable columnar table, optional Plotly line plot, CSV export button. |
| **ZIP bundle** | The CSV is exported in the ZIP bundle and restored on upload (round-trip identical to the parametric study). |

---

## 7. Limitations and future work

### 7.1 Current limitations

- **One integration variable and one interval per model.** All `INTEGRAL`
  calls must share the same `t` and `[t0, tf]`. Nested (multi-variable)
  integration is detected and rejected.
- **Semi-explicit index-1 DAE only.** Higher-index systems are detected
  and reported; index reduction is not yet implemented.
- **2-arg table-based `INTEGRAL(integrand, var)`** is not supported — it
  requires an EES-style Parametric table that CoolSolve does not have.
- **Stiff systems.** An explicit method may need many steps on a stiff
  system; use `EulerImplicit` as a workaround until a BDF/stiff integrator
  is added.

### 7.2 Planned improvements

| # | Improvement | Priority |
|---|-------------|----------|
| 1 | **BDF / DASSL-style stiff variable-step integrator** — slots in via the existing `Integrator` factory | High |
| 2 | **Pantelides-style index reduction** — enable fully-implicit (high-index) DAE models | Medium |
| 3 | **Event detection** (zero-crossing) inside the RK45 step for precise stop locations | Low |
| 4 | **`INTEGRALVALUE` full evaluator dispatch** — wire the active-table context during a step | Low |

---

## 8. Reference implementations

Bundled examples:

- `examples/integral_decay.eescode` — exponential decay `dy/dt = -y`,
  `y(0) = 1`, analytical solution `y(t) = e^{-t}`.
- `examples/integral_free_fall.eescode` — free fall with drag (EES
  example from the EES documentation).
- `examples/integral_harmonic_oscillator.eescode` — coupled `y' = z`,
  `z' = -y` (harmonic oscillator, tests multi-state coupling).
- `examples/integral_coupled_algebraic.eescode` — state coupled to an
  algebraic variable (heat-transfer cell, exercises the algebraic-resolve
  path).
