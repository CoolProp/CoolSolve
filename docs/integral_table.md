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

More generally, as in EES, `INTEGRAL(f, t, lo, hi[, step])` is a **value** — the
integral of `f` from `lo` to the current `t` — that may appear anywhere in an
equation (`y = m_dot*INTEGRAL(…)`, `INTEGRAL(…)/3.6E6`, two calls in one equation…);
at the end of the run every variable holds its value at `t = hi`. The plain additive
form above makes `y` a state of the integrator; every other form gets an
**auxiliary state** per call (§2.4).

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
   classify state vs. algebraic variables (one state per `INTEGRAL` call, user
   or auxiliary — §2.4), partition blocks.
2. **Build a reduced IR** (integral equations removed + driver equations
   added) and re-analyse it. This makes the algebraic subsystem square.
3. **Resolve the interval**: constant limits and step are folded at extraction;
   symbolic ones are read from the model solved at the initial point (§2.5).
4. **Compute the initial state** at `t = t0`: the integral term has
   zero-width interval, so `y(t0)` falls out of the algebraic solve of the
   remaining equations.
5. **March in time**: at each step, fix `(t, y)` as externals, call the
   algebraic `Solver` on the blocks of the reduced system that depend on `t` or
   on a state to evaluate `f` (§3.3), then advance `y` using the chosen
   `Integrator`. The number of steps of a fixed-step method is set by the model
   step (or `integralMaxSteps` without one); RK45 and BDF choose their steps
   within a budget. A march that ends before `tf` is an error, never a success.
6. **Record rows** into an `IntegralTable` (columns matched to variables without
   regard to case): at every step without an output interval, else at
   `t0 + k·n` and at `tf`, independent of the steps. A row between two steps is
   interpolated linearly (`IntegralTable::interpolate`) between the rows of the
   two steps around it, evaluated once each (the row of a step is reused by the
   next one); a row on a step is that step's row. With BDF, whose steps are long,
   the states at a row between two steps come from the method's own interpolant
   (`Integrator::denseOutput`) and the model is solved there: every column of the
   row is that of the model at that time, an input that switches inside the step
   included.

The key invariant: **the algebraic `Solver` is reused unmodified at each time step**.
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
| **Euler implicit** | Newton on `y_{n+1} − y_n − h·f(t_{n+1}, y_{n+1}) = 0` | 1 + 2 or 3 iterations (a Jacobian, n evaluations, when renewed) | A-stable; 1st order; a step that does not converge is an error |
| **RK4** | 4 stages, classic Butcher tableau | 4 | 4th order; **default** |
| **RK45** | 6 stages, Dormand–Prince embedded pair (4th + 5th) | 6 | Adaptive; error-estimation step control |
| **BDF** | Variable-order (1–5) NDF, `y_{n+1}` from the backward differences | 2 or 3 iterations (a Jacobian, n + 1 evaluations, when renewed) | Adaptive step and order; for stiff models |
| **Richardson** | Combine one `h` step + two `h/2` steps | triples | Reduces truncation error; fixed-step only |

The general Richardson combination is `(2^p·I_{h/2} − I_h)/(2^p − 1)`
where `p` is the base method's order (the EES formula `(4·I_{h/2} − I_h)/3`
is the `p = 2` special case).

#### BDF and the implicit Newton iteration

`BDF` (`integralMethod = BDF`) is the variable-order, variable-step method of
Shampine and Reichelt (the numerical differentiation formulas NDF1–NDF5 of
MATLAB's `ode15s` and scipy's `BDF`). It is meant for stiff models: those whose
fastest time constant, not the accuracy asked, sets the step of an explicit
method. It keeps from one step to the next the backward differences of the
solution, a finite-difference Jacobian of `f` and the factorisation of the
Newton matrix `I − c·J`. Each step:

- predicts the state from the differences and solves the implicit formula with a
  simplified Newton iteration (`src/integral/implicit_newton.h`), tested with the
  weights `integralAbsTol + integralRelTol·|y|`; when it does not converge, the
  Jacobian is renewed and the iteration tried again; if it still does not
  converge, the step is rejected and halved;
- estimates its local error from the Newton correction and is rejected and
  shortened above the tolerance;
- after `order + 1` steps of the same size, compares the errors of the orders
  `k − 1`, `k`, `k + 1` and takes the order and step that promise the longest step.

The first step is chosen from `f` at the initial point, within the step of the
model (the 5th argument of `INTEGRAL`, else a tenth of the interval).
`integral.md` gives the number of evaluations of the right-hand side, those of the
Jacobians and of the rejected steps included.

Its gain is the accuracy given up: on the greenhouse models of the library (10
and 15 states, 6 days), `integralRelTol = 1e-4` (`integralAbsTol = 1e-6`) keeps
the hourly rows within 0.03 K of RK4 at 15 s with 2 900 and 4 800 evaluations,
against 17 280 and 34 560 for RK4 at the steps of the models (0.13 s and 0.18 s
per run instead of 0.6 s and 0.9 s); the default `1e-6` gives 0.0005 K with
10 700 and 13 400 evaluations, where RK4 at the steps of the models is within
0.003 K and 0.02 K.

`EulerImplicit` uses the same iteration with its fixed step: the Jacobian kept
from step to step, renewed at the predicted point when the iteration does not
converge, and, since a fixed step cannot be shortened, the full Newton iteration
(a Jacobian at every iterate) when even the renewed Jacobian fails, as where a
control starts to act between the predictor and the solution. A step that still
does not converge is an error.

**Where an implicit method fails.** It needs a right-hand side continuous in the
states. A jump that the trajectory crosses is passed with short steps. A jump on
which the trajectory stays is not: the intake of `examples/engine_weibe_cycle`,
`IF(p - p_out, 0, m_dot)`, switches the flow off as soon as the pressure reaches
the outside pressure and on again below it, so the pressure follows the outside
pressure with a flow that switches at every step (a sliding mode). The implicit
step then has no solution however short it is: BDF and `EulerImplicit` stop
there with a message that says so. RK4 steps across the switch, with a flow
that alternates from step to step and an error of the order of the step; RK45
follows the switch with very short steps and needs `integralMaxSteps` of about
50 000 on that example.

**The smallest step.** An error-controlled method (BDF, RK45) passes a jump of `f`
at the initial point, with a state that starts at 0, only with a first step short
enough for `integralAbsTol` to dominate the weight of that state (below 1e-9 s
for the 50 000 s of `examples/ice_storage_tank`). The automatic minimum step is
therefore the resolution of the integration variable (ten times the spacing of
the floating-point numbers at `t`), not a fraction of the interval;
`integralMinStep` sets another one.

### 2.4 Auxiliary states (an `INTEGRAL` inside a larger expression)

An equation is a *user state* when it has the plain additive form
`y = base + INTEGRAL(...)` (one call, reached by additions and subtractions with a net
plus sign) **and** its `base` is constant along the trajectory: the base does not use
the integration variable, another state, or any variable that depends on them
(checked on the reduced system's block structure: the dependents of the integration
variable and of the states are the closure of those seeds through the blocks). The
equation is then removed and `y` is driven by the integrator, exactly as above — this is
the historical path, with unchanged cost and results.

Every other form — `y = m_dot*INTEGRAL(...)`, `-INTEGRAL(...)`, `INTEGRAL(...)/3.6E6`,
two calls in one equation, a call on the left-hand side, a base that varies — is handled
by one **auxiliary state** `I_k` per call: an internal variable (`__integral_<k>`, named
so that it cannot clash with the model's variables) with `I_k(lo) = 0` and
`dI_k/dt = integrand_k`. The call is replaced by `I_k` in a copy of the equation, which
stays in the algebraic subsystem as an ordinary equation (`y = m_dot*I_k`); the
integrator owns `I_k` the way it owns a user state. An integrand that is not a plain
variable gets an equation of its own (`__integrand_<k> = expr`), so it is evaluated by
the algebraic solver like any other expression (arrays, property calls, user
functions). Auxiliary states and these internal variables are removed from the final
result: they are not in the `.sol`, the `-integral.csv`, the JSON or the GUI table; the
default table shows the left-hand side `y` instead.

### 2.5 Symbolic limits and step

When a limit or the step of an `INTEGRAL` call is not a constant, the solver builds a
second reduced system for the initial point: the states are pinned at their initial
values (auxiliary states at 0), the integration variable is defined by the equation
`t = <lower limit>`, and each limit / step is defined by an equation of its own
(`__int_lo_k = <expr>`, `__int_hi_k`, `__int_step_k`), so that the algebraic solver
evaluates them like any other variable. One solve gives `t0`, `tf` and the step.
Because the limits must be the same along the whole integration, their structural
dependence on the integration variable or on a state is checked on this system (same
closure as §2.4) and rejected with a clear message. All calls must resolve to the same
interval. Models with constant limits and step never build this system.

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
├── integrator_bdf.cpp               Variable order and step, BDF (NDF) for stiff models
├── implicit_newton.h / .cpp         Newton iteration and Jacobian of BDF and EulerImplicit
├── richardson.cpp                   Richardson extrapolation helper (shared)
├── integral_solver.cpp              Time-march loop, reuses algebraic Solver
└── integral_extraction.cpp          Walk IR, classify state vs algebraic vars

tests/
├── test_integrator_euler.cpp        TDD: Euler explicit/implicit vs analytical
├── test_integrator_rk4.cpp          TDD: RK4 vs analytical
├── test_integrator_rk45.cpp         TDD: RK45 + adaptive step control
├── test_integrator_bdf.cpp          BDF, its interpolant and the march with it; EulerImplicit's iteration
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
2. Call `Solver::solveBlocks()` on the blocks that depend on `t` or on a state
   (the closure of `t` and the states through the blocks, `dependentVariables`).
   The other blocks (parameters, constants of the model and everything computed
   from them only) were solved by the complete solve at `t0` and keep their
   values: they are not solved again at every stage.
3. Read the integrand variables (`dydt`, ...) from the state of the solver — those
   are `f`. No `SolveResult` is built for a stage: a table row reads its columns
   from the state too, and the result at `tf` is the state after the last row.

Each block reads the variables it does not own from the state by reference
(`BlockEvaluator::evaluate`), and an explicit block `x = expr` is evaluated without
derivatives, so the cost of a stage is that of the evaluation of its equations: for
the greenhouse models of the library (about 200 blocks, no property call) a stage
takes 0.1 ms (docs/solver_roadmap.md §6.1).

The lookup tables of the model (companion `<model>-<table>.csv` files, or the tables of
the GUI session) are owned by the runner and given to `IntegralSolver`, which hands them
to every algebraic `Solver` it builds (the march and the system that resolves symbolic
limits). `INTERPOLATE` / `LOOKUP` therefore work in the equations of an `INTEGRAL` model
and in its `FUNCTION`s / `PROCEDURE`s, with the weather (or any forcing) table read at
each step (`examples/building_rc_network.eescode`).

For **implicit Euler**, steps 1–3 are wrapped inside the nonlinear solve
for `y_{n+1}`: a Newton iteration on `G(y_{n+1}) = y_{n+1} − y_n − h·f(t_{n+1}, y_{n+1})`
with a finite-difference Jacobian of `f` (one RHS evaluation per state and per
iteration), at most 30 iterations. It converges when `|G| < 1E-12` or the correction is
below 1E-12 (or, after the 30 iterations, at a few ulps of `y`, the rounding of a large
state). A step that does not converge is an error with the time and the residual: its
state would be off the trajectory.

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
| `integralMethod` | `RK4` | `RK4`, `RK45`, `BDF`, `EulerExplicit`, or `EulerImplicit` |
| `integralFixedStep` | `0.0` | Fixed step size when the `INTEGRAL` call has none. `0` ⇒ `integralMaxSteps` steps (fixed methods) |
| `integralMaxSteps` | `1000` | Number of steps of a fixed-step method without an explicit step; RK45 and BDF budget (4 attempted steps for each, see below) |
| `integralRelTol` | `1e-6` | RK45 and BDF local relative error control; weights of the Newton test of BDF and EulerImplicit |
| `integralAbsTol` | `1e-9` | RK45 and BDF absolute error floor (same uses) |
| `integralMinStep` | `0.0` | Minimum step of RK45 and BDF (0 = auto: the resolution of the integration variable, §2.3) |
| `integralMaxStep` | `0.0` | Maximum step of RK45 and BDF (0 = auto: the interval) |
| `integralRichardson` | `false` | Richardson extrapolation (fixed-step only) |
| `integralOutputInterval` | `0.0` | Default `$IntegralTable` row interval when the directive has no `:n` (`0` = every step) |

**Which step wins.** With a fixed-step method, the step given in the model (5th
`INTEGRAL` argument) is used as written, whatever `integralMaxSteps` says: the number of
steps is `⌈(hi − lo)/step⌉` (the last one is shortened). Next comes `integralFixedStep`.
Without either, the method takes exactly `integralMaxSteps` steps. RK45 and BDF adapt
(the step of the model is their first trial step) within a budget of 4 attempted steps,
accepted or rejected, for each of `integralMaxSteps`, or for each of the
`⌈(hi − lo)/step⌉` steps of the model when there are more. Stopping before `hi` is an error.

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
| **Debug folder** | `integral.md` (problem summary, method, step count, step statistics, evaluations of the right-hand side, trajectory preview) + copy of the CSV. |
| **GUI** | **Integral** tab in the bottom panel: scrollable columnar table, optional Plotly line plot, CSV export button. |
| **ZIP bundle** | The CSV is exported in the ZIP bundle and restored on upload (round-trip identical to the parametric study). |

The auxiliary states and internal variables of §2.4 appear in none of these outputs.
`$IntegralTable` names are matched to the model's variables without regard to case (the
header uses the model's spelling); a name that is not a variable is a warning and gives an
empty column.

---

## 7. Limitations and future work

### 7.1 Current limitations

- **One integration variable and one interval per model.** All `INTEGRAL`
  calls must share the same `t` and `[t0, tf]`. Nested (multi-variable)
  integration and calls with different integration variables are detected and
  rejected (`CS-GAP-INTEGRAL-MULTIVAR`).
- **Increasing integration only**, and **no `INTEGRAL` inside a `DUPLICATE` loop**
  (clear error).
- **Limits and step must stay constant along the integration**: constants, variables
  or expressions of the model that do not depend on `t` or on a state (§2.5).
- **Semi-explicit index-1 DAE only.** Higher-index systems are detected
  and reported; index reduction is not yet implemented.
- **A state whose initial value is an inline lookup** (`y = LOOKUP('t', 1, 2) + INTEGRAL(...)`,
  likewise an inline property call) is refused (*Initial value of a state: Unknown or
  unsupported function*): write the value as a variable (`y0 = LOOKUP('t', 1, 2)`,
  `y = y0 + INTEGRAL(...)`).
- **2-arg table-based `INTEGRAL(integrand, var)`** is not supported — it
  requires an EES-style Parametric table that CoolSolve does not have.
- **Stiff systems.** An explicit method may need many steps on a stiff
  system (a fixed step above its stability limit ends in a failed algebraic solve,
  whose message says that the step may be too large); use `BDF` with a tolerance
  matched to the accuracy needed (§2.3). An implicit method (BDF, `EulerImplicit`)
  needs a right-hand side continuous in the states (§2.3).

### 7.2 Planned improvements

| # | Improvement | Priority |
|---|-------------|----------|
| 1 | **Pantelides-style index reduction** — enable fully-implicit (high-index) DAE models | Medium |
| 2 | **Event detection** (zero-crossing) inside the RK45 and BDF steps for precise stop locations, and to restart BDF at a switch of the right-hand side | Low |
| 3 | **`INTEGRALVALUE` full evaluator dispatch** — wire the active-table context during a step | Low |

---

## 8. Reference implementations

Bundled examples:

- `examples/integral_decay.eescode` — exponential decay `dy/dt = -y`,
  `y(0) = 1`, analytical solution `y(t) = e^{-t}`.
- `examples/integral_free_fall.eescode` — free fall with drag (EES
  example from the EES documentation).
- `examples/building_rc_network.eescode` — 3R2C building over one week, the
  weather read from the companion table `building_rc_network-week.csv` with
  `INTERPOLATE` at every step.
- `examples/integral_harmonic_oscillator.eescode` — coupled `y' = z`,
  `z' = -y` (harmonic oscillator, tests multi-state coupling).
- `examples/integral_coupled_algebraic.eescode` — state coupled to an
  algebraic variable (heat-transfer cell, exercises the algebraic-resolve
  path).
