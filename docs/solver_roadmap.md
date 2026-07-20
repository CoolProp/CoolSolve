# CoolSolve — Solver Architecture & Roadmap

*Last updated: July 2026.*

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
| `solver_newton1d.cpp`      | Specialized 1D root-finding: 4-phase Newton + probing + bisection hybrid |
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
| Superancillaries | Enabled by default (`coolpropEnableSuperancillaries`); disable with `--no-superancillary` for faster VLE solving when critical-point accuracy is not needed. |

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
| **Newton1D** | Automatic for size-1 blocks | 4-phase Newton + probing + bisection hybrid |
| **BisectionND** | Default (skipped when block > `bisectionNDMaxBlockSize`) | Sign-pattern simplex, derivative-free — robust when Jacobian is singular or zero |
| **Homotopy continuation** | Default | Predictor-corrector with adaptive step, Newton + LM corrector fallback, final polish |
| **KINSOL** | Opt-in (`Kinsol`) | In-tree port of the three SUNDIALS KINSOL globalisation modes: inexact Newton + Dennis-Schnabel line search (`linesearch`, default), Picard fixed-point (`picard`), Anderson-accelerated fixed-point (`fp`). Not in the default pipeline. |
| **Partitioned** | Default (last resort) | Per-variable diagonal updates; returns `MaxIterations` for small blocks |
| **Tearing** | Optional (`enableTearing`) | Greedy FVS + acyclic sequential solve + Newton on tear residuals |
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
| `--no-superancillary` flag | Avoids loading the cached superancillary data |

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

Items are listed by priority. Each is independent.

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|-----------------|------------------|
| 1 | **BDF / DASSL-style stiff variable-step integrator** | A genuinely stiff ODE/DAE integrator for the dynamic solver. RK45 takes many steps on stiff systems; BDF reuses the same per-step algebraic-solve infrastructure and slots in via the existing `Integrator` factory | 1–2 weeks |
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
