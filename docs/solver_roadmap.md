# CoolSolve — Solver Architecture & Roadmap

*Last updated: November 2026.*

This document describes the architecture, algorithm portfolio, and
forward-looking roadmap of CoolSolve's solver subsystem.

---

## 1. Architecture

### 1.1 Solver pipeline

CoolSolve uses a configurable multi-solver pipeline. The default
(sequential mode) tries each strategy in order until one succeeds:

**Newton → TrustRegion → LevenbergMarquardt → BisectionND → Homotopy →
Partitioned**

A parallel mode (first-to-solve-wins via `std::async` threads) is also
available. BisectionND is automatically skipped for blocks larger than
`bisectionNDMaxBlockSize` (default 8).

When `enableSymbolicReduction` is enabled, each multi-variable block is
first symbolically reduced (CoolProp call inversion + explicit
extraction + equation substitutions), then automatically re-decomposed
into independent sub-blocks before entering the solver pipeline. See
[Symbolic Block Reduction](symbolic_redecomposition.md).

### 1.2 File structure

| File                       | Content                                                               |
|----------------------------|-----------------------------------------------------------------------|
| `solver.h`                 | All solver class declarations, `SolverOptions`, enums, `Newton1DSolver` |
| `solver_common.h`          | Shared utilities: `computeScalingFactors()`, `NonMonotoneHistory`, `isFatalEvaluationError()` |
| `solver.cpp`               | Orchestrator, tearing, partitioned, pipeline, symbolic-reduction integration |
| `solver_newton.cpp`        | Newton + non-monotone Armijo backtracking line search + Broyden quasi-Newton |
| `solver_newton1d.cpp`      | Specialized 1D root-finding: 4-phase Newton + probing + bisection hybrid, Newton steps that increase \|F\| are rejected |
| `solver_trust_region.cpp`  | Trust Region Dogleg + optional hybrd-style Broyden/QR reuse |
| `solver_lm.cpp`            | Levenberg-Marquardt |
| `solver_bisection_nd.cpp`  | N-dimensional bisection (sign-pattern simplex) |
| `solver_homotopy.cpp`      | Homotopy continuation (predictor-corrector) |
| `solver_kinsol.cpp`        | KINSOL (SUNDIALS-style): inexact Newton + Dennis-Schnabel line search, Picard, Anderson FP |
| `solver_symbolic.cpp`      | Symbolic block reduction (inversion + extraction + substitution) |
| `solver_hybrd_qr.h`/`.cpp` | Incremental QR factorization (`computeInitialQR()`, `rank1QRUpdate()`) |
| `symbolic_reduction.h`     | Symbolic reduction data structures and API |
| `structural_analysis.h`    | Graph decomposition API (Tarjan SCC + re-decomposition) |
| `structural_analysis.cpp`  | Tarjan SCC, Dulmage-Mendelsohn, block re-decomposition |

### 1.3 CoolProp integration

| Concern | Implementation |
|---------|----------------|
| Backend | Low-level `AbstractState` API (2–5× faster than `PropsSI`); `PropsSI` only as a fallback |
| Caching | `thread_local AbstractStateCache` keyed by `(backend, fluid)`. Each parallel-solver thread gets its own cache automatically. Toggle via `coolpropCacheEnabled` (default `true`). |
| Derivatives | `first_partial_deriv()` provides exact analytical derivatives. A forward-FD consistency check (5 % tolerance, 2 extra `update()` calls) catches bad analytical derivatives near phase boundaries for pseudo-pure fluids and falls back to finite differences. Toggle via `coolpropEnableAnalyticalDerivatives` (default `true`). |
| Superancillaries | Enabled by default (`coolpropEnableSuperancillaries`); disable with `--no-superancillary`. |

**Supported CoolProp input pairs** (used by symbolic CoolProp call
inversion):

| Input pair            | Inputs                | Use case                              |
|-----------------------|-----------------------|---------------------------------------|
| `PT_INPUTS`           | Pressure, Temperature | Standard forward evaluation           |
| `HmassP_INPUTS`       | Enthalpy, Pressure    | When `h` is known, solve for `T`      |
| `PSmass_INPUTS`       | Pressure, Entropy     | Isentropic processes                  |
| `HmassSmass_INPUTS`   | Enthalpy, Entropy     | When both `h` and `s` are known       |
| `DmassP_INPUTS`       | Density, Pressure     | When `ρ` is known                     |
| `QT_INPUTS`           | Quality, Temperature  | Saturation properties                 |

---

## 2. Solver algorithms

All algorithms below are implemented, tested, and integrated into the
pipeline (or available as opt-in pipeline entries). "Default" means the
algorithm runs as part of the default `solverPipeline`; "opt-in" means
the user must add it to `solverPipeline` explicitly.

| Algorithm | Mode | Notes |
|-----------|------|-------|
| **Newton + non-monotone line search** | Default workhorse | Grippo-Lampariello-Lucidi (1986) non-monotone Armijo (`lsNonMonotoneMemory`, default 10, bounded by `R × φ_current`); Broyden quasi-Newton rank-1 updates (`broydenRecomputeInterval`, default 0 = always full J) |
| **Trust Region Dogleg** | Default | Adaptive initial radius (Cauchy-step-based), graduated ρ-threshold radius adaptation, gradient-direction recovery after >15 consecutive rejections. Optional hybrd-style Broyden/QR reuse via `trBroydenRecomputeInterval` (default 0 = disabled) — see §4.6 |
| **Levenberg-Marquardt** | Default | Nielsen's λ adaptation; cumulative Marquardt diagonal scaling; geodesic acceleration (Transtrum & Sethna 2012) — one extra residual eval/iter, often halves iteration count on curved problems |
| **Newton1D** | Automatic for size-1 blocks | 4-phase Newton + probing + bisection hybrid. A Newton step that increases \|F\| (e.g. the jump over the pole of `phi = V/(U*A)`) is rejected and retried with half the step; a step that changes the sign of F is kept as a bracket. A bracket that shrinks below `stepTolerance` without a residual below `lsRelaxedTolerance` encloses a jump of F (a pole, the penalty value of a failed CoolProp call), not a root: it is given up and the next bracket tried (Phase 2's, best first) |
| **BisectionND** | Default (skipped when block > `bisectionNDMaxBlockSize`) | Sign-pattern simplex, derivative-free — robust when Jacobian is singular or zero |
| **Homotopy continuation** | Default | Predictor-corrector with adaptive step, Newton + LM corrector fallback, final polish |
| **KINSOL** | Opt-in (`Kinsol`) | In-tree port of the three SUNDIALS KINSOL globalisation modes: inexact Newton + Dennis-Schnabel line search (`linesearch`, default), Picard fixed-point (`picard`), Anderson-accelerated fixed-point (`fp`). Not in the default pipeline. |
| **Partitioned** | Default (last resort) | Per-variable diagonal updates; returns `MaxIterations` for small blocks |
| **Tearing** | Optional (`enableTearing`) | Greedy FVS + acyclic sequential solve + Newton on tear residuals (Schur complement; the acyclic variables follow the tear step to first order). Each outer step (sweep + tear step) is damped by a bounded backtracking line search (non-monotone Armijo shared with Newton, at most four halvings, then the shortest step is kept): a full step no longer throws the iterate onto another root, e.g. a negative leakage velocity in `air_screw_compressor` |
| **Symbolic block reduction** | Optional (`enableSymbolicReduction`) | CoolProp call inversion + explicit extraction + equation substitution; automatic re-decomposition into sub-blocks. See [Symbolic Block Reduction](symbolic_redecomposition.md). |
| **Multi-start fallback** | Default (`multiStartMode = deepsearch`) | Per-block retry from CoolProp-consistent pressure regimes (thermo blocks) or scale factors (algebraic blocks); parallel candidates (`multiStartNumCores`, default 4). Engages on every solve (`always`), only during Deep Search (`deepsearch`, default), or never. See §3. |
| **Explicit solve** | Automatic | `tryExplicitSolve()` detects `x = expr(known)` pattern, evaluates RHS directly |
| **Cancel token propagation** | All solvers | Checked inside iteration loops of every solver + parallel stop signal |

### 2.1 Performance & profiling infrastructure

| Feature | Notes |
|---------|-------|
| `PipelineTiming` | Per-solve timing breakdown, displayed in GUI |
| CoolProp call profiling | Per-call timing for `PropsSI`, `AbstractState`, `HAPropsSI`; warmup tracked separately |
| Variable scaling | Shared `computeScalingFactors()` in `solver_common.h`, used by Newton/TR/LM/BisectionND |
| `--no-superancillary` flag | Disables the CoolProp superancillary functions |

---

## 3. Robustness results

From `examples/solver_robustness_report.md` (42 example files × multiple
configurations, 30 s per-solve timeout):

| Configuration                        | With initials | Without initials |
|--------------------------------------|--------------:|-----------------:|
| Default pipeline (no multi-start)    | 92.7%         | 67.5%            |
| Default pipeline + sequential multi-start (default) | 92.7% | 75.0% |
| Default pipeline + **parallel** multi-start (auto cores) | 92.7% | **80.0%** |

### 3.1 Multi-start contribution

Multi-start is **strictly ≥ baseline**: it never regresses a
with-initials solve (every block converges first try, so multi-start
never engages) and rescues several without-initials solves. Parallel
multi-start adds a small constant-factor benefit because each candidate
gets ~the full wall-clock budget within the 30 s timeout.

The remaining without-initials failures (`orc_co2`, `cooling_tower`,
`zorlu_heat_pump`, …) have initial residuals up to |F|=1e12 from deeply
coupled 5th-order polynomial fits and CoolProp calls whose convergence
basin is too narrow for a handful of starting points — these need good
`.initials` or KINSOL.

### 3.2 Hardest models

The largest and most coupled models in the suite (used for stress
testing new solver work):

`zorlu_heat_pump` (63 vars), `condenser_3zones` (62 vars),
`heat_pump_MSTh_SB_R10` (39 vars), `orc_r245fa` (38 vars),
`orc_co2` (28 vars), `cooling_tower` (22 vars),
`air_screw_compressor` (13-var near-singular compressor-leak block).

To explore: since the CoolProp pin `75af7816` (p,T flash rework, CoolProp #3360, bisected to
`3f06ccf7`), `air_screw_compressor` no longer converges from the default guesses at the default
`tolerance = 1e-9` (SingularJacobian in block 13, tearing at its iteration limit); it converges from
1e-8, and at 1e-9 from its `.initials` (shipped since). Mechanism not analysed (flash precision vs the
residual floor of the leakage block?).

---

## 4. Implementation notes

The sections below describe non-obvious implementation choices that are
still relevant as background reading when modifying these solvers.

### 4.1 Symbolic block reduction + re-decomposition

Controlled by `enableSymbolicReduction` (default `false`). Three
techniques applied iteratively until a fixed point:

1. **CoolProp call inversion**: if the model has
   `h = enthalpy(water, T=T, P=P)` and `h` is the matched output but
   one named-arg input is an unknown block variable while the other is
   external, reformulate as `T = temperature(water, H=h, P=P)`.
2. **Explicit extraction**: equations where all RHS variables are
   external or already-reduced are directly evaluated, removing their
   output variable from the block.
3. **Equation substitution**: variables that appear in zero other
   block equations after extraction are removed along with their
   defining equation.

After reduction, **post-reduction re-decomposition** (Tarjan SCC on the
reduced dependency graph) splits the remaining monolithic block into
independent sub-blocks solved sequentially. Example reduction ratios on
real models:

| Model | Original | Reduced | Sub-blocks |
|-------|---------:|--------:|------------|
| `condenser_3zones` | 62 | 56 | 13 (largest 30) |
| `heat_pump_MSTh_SB_R10` | 39 | 34 | 19 (largest 8) |
| `orc_co2` | 28 | 25 | 13 (largest 8) |
| `air_screw_compressor` | 13 | 6 | 6 (fully decomposed) |

Full algorithm details in [Symbolic Block Reduction](symbolic_redecomposition.md).

### 4.2 Non-monotone line search

Implemented in all three gradient-based solvers (Newton, TrustRegion,
LM). Controlled by `lsNonMonotoneMemory` (default `10`; set to `1` for
classic monotone behaviour).

The monotone Armijo condition is replaced with the Grippo-Lampariello-
Lucidi (1986) non-monotone variant: instead of requiring
`φ(x_{k+1}) < φ(x_k)`, the solver compares against
`max(φ(x_{k-M+1}), ..., φ(x_k))` where `M` is the memory parameter.
This allows the solver to temporarily accept larger merit values to
escape narrow curved valleys and saddle points, with the same global
convergence guarantees.

**Bounded non-monotone reference**: the raw `max(history)` can be
arbitrarily larger than the current merit when a single early bad
iterate inflates the window. To prevent this, all solvers use
`boundedRef(phi, R = 10)` which caps the reference at
`min(max(history), R × currentPhi)`.

When `M = 1`, all three solvers degenerate to exact monotone behaviour.

### 4.3 Trust Region (current state)

- **Adaptive initial radius** (`trAdaptiveRadius`, default `true`):
  `delta = min(trInitialRadius, max(cauchyNorm, 1.0))`.
- **Smoother rho-based radius adaptation**: graduated ρ thresholds
  (grow at ρ > 0.75, hold at ρ > 0.5, shrink proportionally otherwise
  using `max(0.1, 1 − ρ) × delta`).
- **Gradient-based recovery**: after >15 consecutive rejections, the
  solver resets to a gradient-direction step with a small radius.
- **Optional hybrd-style Broyden/QR reuse** via
  `trBroydenRecomputeInterval` (default 0 = disabled) — see §4.6.

### 4.4 Levenberg-Marquardt (current state)

- **Nielsen's λ adaptation** (`lmNielsenUpdate`, default `true`):
  `λ = λ × max(1/3, 1 − (2ρ − 1)³)` on acceptance and exponential
  increase (`λ × ν` with ν doubling) on consecutive rejections.
- **Cumulative Marquardt diagonal scaling**:
  `D_k = max(D_{k-1}, diag(J^T J))`.
- **Geodesic acceleration** (`lmGeodesicAcceleration`, default `true`):
  adds a second-order correction to the LM step by evaluating the
  directional second derivative of `F` along the velocity step. Costs
  one extra residual evaluation per iteration but can halve iteration
  count on curved problems. Includes acceleration/velocity ratio guard.

### 4.5 Broyden quasi-Newton (Newton solver)

Implemented as a hybrid Newton-Broyden mode in `src/solver_newton.cpp`.

Broyden's method avoids recomputing the full Jacobian at every iteration
by maintaining rank-1 updates:

```
B_{k+1} = B_k + (ΔF - B_k Δx) Δx^T / (Δx^T Δx)     [Broyden Type I]
```

- **`broydenRecomputeInterval`** option (default `0` = always full
  Jacobian; `K > 0` = recompute every `K` iterations).
- Computes true Jacobian on iteration 0, then uses Broyden Type-I
  rank-1 updates for intermediate iterations.
- **Automatic fallback**: forces full Jacobian recomputation when
  (a) Broyden step fails line search, (b) residual stalls for
  ≥ 3 iterations, or (c) scheduled recompute interval reached.
- Primarily an **efficiency** improvement (fewer CoolProp calls on
  large blocks). For a 30-var block, full Newton costs ~31 CoolProp
  calls/iteration (30 for Jacobian + 1 for residual). With Broyden
  (`K = 5`), on average only ~7 calls/iteration.

### 4.6 Incremental QR factorization & hybrd-style TrustRegion

The `solver_hybrd_qr.h`/`.cpp` module implements incremental QR
factorization via rank-1 updates (Golub & Van Loan §12.5; the same
numerically-stable building block as MINPACK's `r1updt`/`r1mpyq`).

When `trBroydenRecomputeInterval > 0`, the TrustRegion solver uses the
QR factorization as the Broyden state representation rather than a dense
matrix. Powell's restart criterion (force a full-Jacobian recompute
after `trBroydenRestartRejects` consecutive trust-region rejections,
default 2) is implemented.

**Empirical findings.** Validated non-regressively against the full unit
suite and on a curated subset of the hardest models. The feature causes
zero regressions and the underlying QR module is independently valuable,
but on real CoolProp models (n ≤ 63) the wall-clock benefit is
**neutral** — the saved CoolProp calls do not currently dominate
per-iteration cost, because each Broyden iteration still pays the full
`O(n³)` cost of reconstructing `J = Q·R` and re-running
`ColPivHouseholderQR`. Realizing the speedup would require a genuinely
"hybrd-native" solve against the maintained `R` directly (see §6,
future work).

---

## 5. Dynamic / DAE solving (`INTEGRAL`)

CoolSolve solves equation-based initial-value differential–algebraic
equation (DAE) systems written in the EES integral form. The feature is
delivered (RK4 default + RK45 adaptive + Euler explicit/implicit +
optional Richardson extrapolation), lives in its own module
(`src/integral/`, `include/coolsolve/integral/`), and reuses the
algebraic `Solver` unmodified at each time step.

See [Dynamic Solving](integral_table.md) for the full documentation —
mathematical model, algorithm overview, architecture, configuration,
output, and limitations.

---

## 6. Future work

### 6.1 Evaluation overhead and time marching (November 2026)

**Steps 1 to 4 implemented** (branch `perf/eval-overhead`, under review;
results below). Step 5 is to be decided on those results, step 6 is item 1
of §6.2. Before this work one algebraic solve cost about 0.21 µs × blocks ×
variables whatever the equations contained, and a dynamic model paid it four
times per RK4 step.

#### Diagnosis

Measured on the three greenhouse models of the CoolSolve Library
(`models/buildings/greenhouses`, library commit `eedb4ee`), CoolSolve
`bc7c92c`, Release build, one thread, with `perf` and temporary timers
around the phases of a block solve. The models make no property call.

| Model | CoolSolve | Python reference of the model |
|-------|-----------|-------------------------------|
| `greenhouse_steady_state` (147 equations) | 0.49 s for the process, 4.6 ms for the solve | 0.3 ms for the solve (scipy `root`, 5 unknowns) |
| `greenhouse_climate_dynamic` (10 states, RK4, 4320 steps of 120 s) | 141 s, 8.1 ms per RHS evaluation | 6.4 s for 17280 steps of 30 s, 0.09 ms per RHS evaluation |
| `greenhouse_heating_pipe_dynamic` (15 states, RK4, 8640 steps of 60 s) | 459 s, 13.2 ms per RHS evaluation | 1.9 s for 17280 steps of 30 s, 0.03 ms per RHS evaluation |

The READMEs of the models report 477 s and 1503 s against 20 s and
6.8 s: the same ratios, measured on a run about 3.3 times slower.

Once the states and the time are pinned, the structure is already the
right one: the pipe model is 243 explicit blocks of one equation (Newton
never runs), the climate model 192 explicit blocks and one block of 5
variables (the pipe balance, 2.7 Newton iterations per stage). The time
goes into what is done around each block:

| Share of one RHS evaluation | climate | pipe |
|-----------------------------|--------:|-----:|
| Three copies of all the variables, for every block: `externalVars` in `Solver::solveBlock`, the same map again in `Solver::tryExplicitSolve`, then one `setVariable` per variable into a new `ExpressionEvaluator` in `BlockEvaluator::evaluate` | 83% | 85% |
| Evaluating the expressions of the model | ≈ 2% | 2% |
| The implicit block, its evaluations included | 2% | – |
| Rest: freeing the copies, per-block trace and result | ≈ 12% | 13% |

- The cost is 0.21 µs per block and per variable in the three models
  (40 µs per block with 197 variables, 53 µs with 243): it grows with
  the square of the model size and is paid by every solve, also in a
  parametric table or a multi-start, not only by `INTEGRAL`.
- The copies exist because `SystemEvaluator` holds its state in a
  `std::map` with `CaseInsensitiveLess` and `BlockEvaluator::evaluate`
  takes a `std::map<std::string, double>`: each caller converts.
- Every explicit evaluation computes a gradient and a residual scale
  that are not read (`computeJacobian` defaults to `true`).
- 93 of the 193 march blocks of the climate model and 97 of the 243 of
  the pipe model depend neither on the time nor on a state, and are
  solved again at every stage.
- The steady-state solve has the same profile: 3.6 ms for its 111
  explicit blocks, 1.0 ms for the 36-variable Newton block (7
  iterations). The other 0.47 s of the process is `warmupCoolProp()`,
  which loads the CoolProp fluid library for a model that never uses
  it. `--no-superancillary` does not shorten it: the load builds the
  superancillaries of every fluid (about 80% of its time).
- Not the cause: Newton, the Jacobian and the linear algebra (2% at
  most), RK4 itself, the length of the run (the time per stage is
  constant). `enableTearing` and `enableSymbolicReduction` leave the
  time of the 5-variable block and of the 36-variable block unchanged.
- Integrators: the system is mildly stiff. RK45 keeps an average step
  of 177 s at `integralRelTol` 1e-6 and 198 s at 1e-4 (33 and 36
  rejections over 6 h), so stability sets its step and it needs more
  RHS evaluations than RK4 at 120 s. `EulerImplicit` needs 28 RHS
  evaluations per step at 120 s (a finite-difference Jacobian at every
  inner iteration).

A throwaway prototype of steps 1 and 2 below (explicit blocks only,
kept outside the repository) gave, with trajectory tables identical to
the stock ones to the last digit:

| | Stock | Step 1 | Steps 1 + 2 |
|---|------:|-------:|------------:|
| Climate, full run | 141 s | 8.7 s | 6.6 s |
| Pipe, full run | 459 s | 14.5 s | 10.6 s |
| Per RHS evaluation, climate / pipe | 8.1 / 13.2 ms | 0.46 / 0.40 ms | 0.35 / 0.29 ms |
| Steady-state: solve / process (with step 3) | 4.6 ms / 0.49 s | 1.3 ms / 0.01 s | – |

#### Plan

Steps 1 to 4 are independent and small. Step 5 is to be decided on the
measurements made after steps 1 and 2.

| # | Step | Change | Gain | Effort |
|---|------|--------|------|--------|
| 1 | **Evaluate a block without copying the state** | `BlockEvaluator::evaluate` reads the external variables from the state of `SystemEvaluator` by reference and binds only the variables its equations refer to (`EquationInfo::variables`, listed once in the constructor); any other name is looked up in the state when it is read, instead of defaulting to 1. `Solver::solveBlock` tries the explicit evaluation before anything else. The five copy loops (`solveBlock`, `tryExplicitSolve`, the zero-unknown check, `SystemEvaluator::evaluateBlock`, `publishStringOutputs`) become that one path. Explicit blocks are evaluated without the Jacobian. The Newton path gets the same binding. The function and procedure tables are referenced, not copied into each new `ExpressionEvaluator` (seen in the code, not exercised by these models) | Measured: 16× (climate), 32× (pipe), 3.5× on the steady-state solve | 2–3 days |
| 2 | **March only what changes** | The blocks that depend neither on the integration variable nor on a state are solved once at the initial point: `dependentVariables()` of `integral_solver.cpp` already gives the set, and `Solver::solve` takes the list of blocks to solve. `IntegralSolver::evaluateRHS` reads the integrands from the evaluator state; a `SolveResult` is built only for a table row and for the final point (today one per stage, copied into `lastAlgebraic_`: 15–20 µs, about 5% of a stage after step 1) | Measured for the first part: a further 24% (climate) and 27% (pipe) | 1–2 days |
| 3 | **No CoolProp warm-up without a property call** | `CoolSolveRunner::run` calls `warmupCoolProp()` after the IR is built, and only when the model calls a property function, directly or in a routine (the test of `SystemEvaluator` for the unit hints, `withCallingRoutines`). The GUI server keeps its background warm-up | Measured: 0.49 s → 0.01 s for the steady-state model. A model with property calls keeps the 0.45 s: building every superancillary at load is CoolProp's, to raise upstream | 0.5 day |
| 4 | **Defects found on the way** | (a) `EulerImplicit` accepts a step whose inner Newton iteration has not converged: at a step of 600 s on the climate model all 36 steps end with a relative residual of order 1 after 30 iterations, the air temperature after 6 h is 20.88 °C instead of 18.93 °C and the status is success. It must fail with the time and the residual. (b) `integral.md` always prints `Method: RK4` (`runner.cpp` builds the label from a default `SolverOptions`). (c) A fixed step too large for an explicit method ends as `Block … LineSearchFailed` in the algebraic solve: the message should say that the step may be too large | Correctness and diagnosis | 0.5 day |
| 5 | **Compiled block evaluation** | After steps 1 and 2 a stage is still 4 to 10 times the Python one. What is left is the evaluation itself: the AST is walked with a new `ExpressionEvaluator` per call, each variable reference is a case-insensitive string lookup (for an array element, of a name built from its indices), each intermediate value is an `ADValue` with a heap-allocated gradient. That is 0.6–1.1 µs per explicit equation, 60 µs per evaluation of the 36-variable block and 135 µs per solve of the 5-variable block (36% of a climate stage after step 2). Compile each block once into a flat list of operations on integer slots of one value array, without gradients when no Jacobian is asked; a node the compiler does not handle (CALL, user function, string) keeps the present evaluator | Estimated, not measured: 10× or more on a stage | 1–2 weeks |
| 6 | **Stiff integrator** (item 1 of §6.2) | What these models ask of it: a Jacobian kept from step to step, a convergence test weighted with `integralRelTol` / `integralAbsTol` (`EulerImplicit` uses an absolute 1e-12), a step rejected and shortened when the iteration does not converge, an error-controlled step. It reduces the number of steps (the weather is hourly, the explicit stability limit is below 300 s), not the cost of a stage: it comes after steps 1 and 2 | Not measured | as §6.2 |

#### Results of steps 1 to 4

Same models, same machine (Release, one thread), full runs; every
`-integral.csv` and `.sol` file is identical to the stock one, byte for byte.

| | Stock | Steps 1 + 3 | + step 2 | + two hot spots found on the way |
|---|------:|-------:|------:|------:|
| Climate, full run | 141 s | 3.90 s | 2.59 s | **1.77 s** (80×) |
| Pipe, full run | 459 s | 8.57 s | 6.09 s | **4.01 s** (114×) |
| Per RHS evaluation, climate / pipe | 8.1 / 13.2 ms | 0.22 / 0.25 ms | 0.15 / 0.18 ms | **0.10 / 0.12 ms** |
| Steady state: solve / process | 4.96 ms / 0.53 s | – | – | **1.65 ms / 0.01 s** |

The targets of the acceptance (climate under 8 s, pipe under 12 s, the
steady-state process under 0.05 s) are met after steps 1 and 3. Against the
Python references: the climate run is now 3.6 times shorter than the Python
one (which takes four times more steps) and costs about the same per RHS
evaluation (0.10 ms against 0.09 ms); the pipe costs 4 times the Python per
evaluation (0.12 ms against 0.03 ms) and the run takes 2.1 times longer.

What was done, beyond the plan:

- Step 1 binds the external variables of an equation once per evaluation and
  reads any other name from the state; an explicit assignment evaluated for
  its value carries no gradient at all (`BlockEvaluator::isExplicitAssignment`).
  Newton1D, which probes around the values of the other variables when its
  Newton phase fails, reads them from the state only then, in the former order
  (`Newton1DSolver::solve` with the state and the block variable).
- Step 2: `Solver::solveBlocks` solves a list of blocks without building a
  `SolveResult`; `solve()` shares its per-block work with it
  (`solveSystemBlock`, `blockFailureMessage`). The march blocks are computed
  once in the `IntegralSolver` constructor; the result at `tf` is the state
  after the last row.
- Hot spots found while measuring after step 2: every call of a function that
  is not a property function (`exp`, `max`, `IF`, ...) threw and caught an
  exception to find out (`getCoolPropOutputParam`; now a lookup returning an
  optional): 30% of a stage. The check of the clamped inputs after each block
  made a locked exchange on an atomic flag that is nearly always clear: 9% of
  a stage of the pipe model.
- Step 4: (a) EulerImplicit keeps its convergence tests (|G| < 1e-12, or a
  correction below 1e-12) and accepts, after its 30 iterations, a last
  correction within 4 ulps of y (a state above 1e3 cannot get a correction
  below 1e-12); anything else is an error with the time, the step and the
  largest residual. On the climate model the steps of 120 s converge as
  before (identical table) and a step of 600 s now fails at t = 0 with a
  residual of 232. (b) `integral.md` names the method used
  (`IntegralSolveResult::method`). (c) A failed algebraic solve during the march
  of `RK4` or `EulerExplicit` adds that the step may be too large.

What remains of a stage (`perf` on the pipe model, self time): the
evaluation itself. About a third is the case-insensitive lookups of variable
names in `std::map`s (the evaluator's local map, the state, `setVariableValue`,
`getVariableValue`), an eighth the strings, allocations and `ExpressionEvaluator`
built for each evaluation, most of the rest the walk of the AST with an
`ADValue` per node. That is what step 5 would remove. Estimate, not measured:
it would bring the pipe model to the Python cost per evaluation, a factor of
about 4 on a stage, and less on the climate model.

#### Acceptance

- The three models are the benchmark, each step timed against the
  tables above. Target after steps 1 to 3: climate under 8 s, pipe
  under 12 s, the steady-state process under 0.05 s.
- Steps 1 to 3 change no number: the `-integral.csv` and `.sol` files
  of the three models are reproduced to the last digit, the unit tests
  pass and `examples/solver_robustness_report.md` is unchanged. Step 1
  changes how a block sees the variables it does not own: the library
  regression is the check for CALLs, user functions, string variables
  and array elements with a computed index.
- To measure again: `perf record --call-graph dwarf` on the short
  models (the stacks need DWARF unwinding in a Release build), and the
  block evaluation count and time of `profiling.md` (`-d`).

### 6.2 Other items

Items are listed by priority. Each is independent.

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|-----------------|------------------|
| 1 | **BDF / DASSL-style stiff variable-step integrator** | A genuinely stiff ODE/DAE integrator for the dynamic solver. RK45 takes many steps on stiff systems; BDF reuses the same per-step algebraic-solve infrastructure and slots in via the existing `Integrator` factory. Requirements from the greenhouse models: §6.1, step 6 | 1–2 weeks |
| 2 | **Pantelides-style index reduction** | Enable fully-implicit (high-index) DAE models in the dynamic solver. Currently high-index systems are detected and rejected; dummy-derivative index reduction would make them solvable | 2–3 weeks |
| 3 | **Hybrd-native `O(n²)` triangular solve for TrustRegion Broyden mode** | Realize the speedup §4.6 found missing: solve directly against the maintained `R` instead of reconstructing dense `J` and re-running `ColPivHouseholderQR` every iteration | 3–5 days, plus a careful robustness re-validation (loses rank-revealing pivoting, reopening the numerical-instability concern) |
| 4 | **Pseudo-arclength continuation in Homotopy** | Pass turning points (where the homotopy path folds back in `t`), a known failure mode for hard thermodynamic models with non-monotonic `t`-paths. Keller (1977) augmented system with the Tier-1 hybrd-style TrustRegion as corrector | 1–2 weeks |
| 5 | **Dynamic column scaling in Newton and LM** | Extend the `diag ← max(diag, ‖J[:,j]‖)` rule (already in TR) to Newton and LM via a shared `updateScalingFactors()` helper in `solver_common.h`. Gated by a renamed `dynamicScaling` option | 3–5 days |

### Not recommended

| Improvement | Why skip |
|-------------|----------|
| **Anderson acceleration as a pipeline solver** | Does not use the Jacobian → slower than Newton when Newton works, and the pipeline already has 6 fallback strategies. The "no Jacobian needed" advantage is addressed better by BisectionND. Broyden is the better quasi-Newton approach since it starts from a true Jacobian. |
| **Code complexity monitor** | The codebase is well-structured after the file split. LOC tracking is noise, not signal. |
| **Nonlinear preconditioning (residual weighting)** | Variable scaling (already implemented) addresses the same problem. Residual weighting adds complexity for marginal benefit given that analytical derivatives also improve conditioning. |
| **Better initial guess via simplified model** | Requires model-specific knowledge and a meta-solving infrastructure. The pipeline + homotopy + multi-start already handle poor initial guesses well. Document manual workflows instead. |

---

## 7. References

**Trust-region / hybrd / Broyden:**

- Powell, M. J. D. (1970): *A Hybrid Method for Nonlinear Equations*.
  In P. Rabinowitz (ed.), *Numerical Methods for Nonlinear Algebraic
  Equations*, Gordon and Breach. (The `scipy.optimize.fsolve` algorithm.)
- Moré, J. J., Garbow, B. S., Hillstrom, K. E. (1980): *User Guide for
  MINPACK-1* (ANL-80-74), Argonne National Laboratory.
- Broyden, C. G. (1965): *A class of methods for solving nonlinear
  simultaneous equations*. Math. Comp. 19.
- Dennis, J. E., Moré, J. J. (1977): *Quasi-Newton methods, motivation
  and theory*. SIAM Review 19.
- Kelley, C. T. (2003): *Solving Nonlinear Equations with Newton's
  Method* (SIAM), Chapter 4.

**Linear algebra:**

- Golub, G. H., Van Loan, C. F. (2013): *Matrix Computations*, §12.5
  for rank-1 QR updates.

**Continuation:**

- Keller, H. B. (1977): *Numerical Solution of Bifurcation and
  Nonlinear Eigenvalue Problems*. In *Applications of Bifurcation
  Theory*, Academic Press.
- Allgower, E. L., Georg, K. (2003): *Introduction to Numerical
  Continuation Methods*, SIAM.

**Line search / least squares:**

- Grippo, L., Lampariello, F., Lucidi, S. (1986): *A nonmonotone line
  search technique for Newton's method*. SIAM J. Numer. Anal. 23.
- Zhang, H., Hager, W. W. (2004): *A nonmonotone line search technique
  and its application to unconstrained optimization*. SIAM J. Optim. 14.
- Madsen, K., Nielsen, H. B., Tingleff, O. (2004): *Methods for
  Non-Linear Least Squares Problems* (2nd ed.), DTU. Source of
  Nielsen's λ update.
- Transtrum, M. K., Sethna, J. P. (2012): *Geodesic acceleration and
  the small-data limit of nonlinear least squares*. PNAS 109.

**DAE / integration:**

- Hairer, E., Wanner, G. (1996): *Solving Ordinary Differential
  Equations II: Stiff and Differential-Algebraic Problems* (Springer).
- Dormand, J. R., Prince, P. J. (1980): *A family of embedded Runge-
  Kutta formulae*. J. Comp. Appl. Math. 6.
- Pantelides, C. C. (1988): *The consistent initialization of
  differential-algebraic systems*. SIAM J. Sci. Stat. Comput. 9.
