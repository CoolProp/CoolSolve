# CoolSolve — Solver Roadmap & Status

*Last updated: June 2026 — Tier 2.1 (incremental QR) and Tier 1 (hybrd-style
TrustRegion + Broyden) both implemented, on the corrected order (QR first).
Validated non-regressive; empirically neutral on the tested models — see
§3.7.*

This document consolidates the performance, robustness, and architecture
plans for CoolSolve's solver subsystem. It replaces the earlier
`performance_plan.md`, `solver_fix_strategy.md`, and
`solver_robustness_diagnosis.md`.

---

## 1. Diagnosis: Why "TR + Broyden" Currently Does Not Improve Robustness

The algorithm behind `scipy.optimize.fsolve` is **MINPACK `hybrd`
(Powell 1970)**. Its well-known robustness comes from running **trust
region and Broyden at every iteration**, not as alternative fallback
strategies. CoolSolve currently keeps them in **separate, never-combined
solvers**, which is why their individual contributions look
underwhelming.

| CoolSolve solver today           | Jacobian                       | Step acceptance            |
|----------------------------------|--------------------------------|----------------------------|
| `Newton` (`broydenRecomputeInterval > 0`) | Broyden rank-1         | Line search                |
| `TrustRegion`                    | Full J **every iteration**     | Dogleg trust region        |
| `LM`                             | Full J every iteration         | LM damping                 |
| **`hybrd` (the missing piece)**  | **Broyden, full J only after 2 consecutive TR rejections** | **Dogleg trust region** |

Concrete flaws that explain the symptom:

1. **Broyden is disabled by default** (`broydenRecomputeInterval = 0`,
   `solver.h:114`). Most users never benefit from it.
2. **The TrustRegion solver never reuses a Jacobian**: it recomputes `J`
   unconditionally on every iteration
   (`solver_trust_region.cpp:108-110`), discarding the rank-1
   information that makes hybrd cheap and stable on large blocks.
3. **Restart criterion is wrong**: TR's gradient-recovery kicks in after
   **15 consecutive rejections** (`solver_trust_region.cpp:269-278`).
   Powell's criterion is **2 consecutive rejections** (`if (ncfail .eq.
   2) go to 290` in `hybrd.f`). CoolSolve keeps a stale Jacobian roughly
   7× too long.
4. **No dynamic scaling**: hybrd updates the column scaling every full-J
   iteration (`diag(j) = max(diag(j), wa2(j))` in `hybrd.f:160`).
   CoolSolve computes scaling once at iteration 0
   (`solver_common.h:26-40`).

Tier 1 below fixes all four flaws inside the existing
`TrustRegionSolver`, without adding any new solver slot, dependency, or
pipeline entry.

---

## 2. Architecture Overview

### 2.1 Solver Pipeline

CoolSolve uses a configurable multi-solver pipeline. The default
(sequential mode) tries each strategy in order until one succeeds:

**Newton → TrustRegion → LevenbergMarquardt → BisectionND → Homotopy →
Partitioned**

A parallel mode (first-to-solve-wins via `std::async` threads) is also
available. BisectionND is automatically skipped for blocks larger than
`bisectionNDMaxBlockSize` (default 8).

When `enableSymbolicReduction` is enabled, each multi-variable block is
first symbolically reduced (CoolProp call inversion + explicit
extraction + equation substitution), then automatically re-decomposed
into independent sub-blocks before entering the solver pipeline.

### 2.2 File Structure

| File                       | Lines | Content                                                               |
|----------------------------|------:|-----------------------------------------------------------------------|
| `solver.h`                 | ~906  | All solver class declarations, `SolverOptions`, enums, `Newton1DSolver` |
| `solver_common.h`          | ~140  | Shared utilities: `computeScalingFactors()`, `NonMonotoneHistory`, `isFatalEvaluationError()` |
| `solver.cpp`               | ~2414 | Orchestrator, tearing, partitioned, pipeline, symbolic reduction integration |
| `solver_newton.cpp`        | ~308  | Newton + non-monotone Armijo backtracking line search + Broyden quasi-Newton |
| `solver_newton1d.cpp`      | 355   | Specialized 1D root-finding: 4-phase Newton + probing + bisection hybrid |
| `solver_trust_region.cpp`  | 297   | Trust Region Dogleg (target of the Tier 1 upgrade, §3)                |
| `solver_lm.cpp`            | 268   | Levenberg-Marquardt                                                    |
| `solver_bisection_nd.cpp`  | 450   | N-dimensional bisection (sign-pattern simplex)                         |
| `solver_homotopy.cpp`      | 209   | Homotopy continuation (predictor-corrector)                            |
| `solver_symbolic.cpp`      | 611   | Symbolic block reduction (inversion + extraction + substitution)       |
| `symbolic_reduction.h`     | 147   | Symbolic reduction data structures and API                             |
| `structural_analysis.h`    | 155   | Graph decomposition API (Tarjan SCC + redecomposition)                 |
| `structural_analysis.cpp`  | 763   | Tarjan SCC, Dulmage-Mendelsohn, block re-decomposition                 |

### 2.3 Test Results (current baseline)

`./coolsolve_tests` runs 274 unit tests (2832 assertions). All
passing. The comprehensive example suite solves **36 / 42** `.eescode`
files (≈ 2.07 s total).

### 2.4 Robustness Results Summary

From `examples/solver_robustness_report.md` (20 configurations × 42
files, with a 30 s per-solve timeout):

| Configuration                        | With initials | Without initials |
|--------------------------------------|--------------:|-----------------:|
| Default pipeline (6 solvers)         | 92.7%         | 67.5%            |
| Newton only                          | 82.9%         | 50.0%            |
| TrustRegion only                     | 87.8%         | 47.5%            |
| LevenbergMarquardt only              | 75.6%         | 40.0%            |
| Homotopy only                        | 82.9%         | 57.5%            |
| BisectionND only                     | 41.5%         | 35.0%            |
| Partitioned only                     | 41.5%         | 35.0%            |
| Default + Tearing                    | 92.7%         | 67.5%            |
| Default + SymbolicReduction          | 92.7%         | 67.5%            |
| Default + Tearing + SymbolicReduction | 92.7%        | 67.5%            |

The hardest models (largest algebraic blocks, near phase-boundary
CoolProp calls, near-singular Jacobians) are: `zorlu_heat_pump`
(63 vars), `condenser_3zones` (62 vars), `heat_pump_MSTh_SB_R10`
(39 vars), `orc_r245fa` (38 vars), `orc_co2` (28 vars),
`cooling_tower` (22 vars).

The **without-initials** column (67.5%) is the real robustness gap and
is addressed by Tier 2.2 (multi-start).

This table is the full 42-file × 20-configuration baseline and has not
been regenerated wholesale (the full suite takes too long to re-run for
every change). It remains accurate because Tier 1 (§3) defaults to off.
For a targeted, curated-subset comparison of legacy `TrustRegion` vs.
the new hybrd mode specifically, see §3.7.

---

## 3. Tier 1 — Hybrd-Style TrustRegion + Broyden (the `scipy.fsolve` Algorithm)

**Status (June 2026): implemented (take 2), on top of the QR infrastructure
from §4.1, following the corrected order from §3.6.** Validated
non-regressive against the full unit suite and examples-comprehensive.
Empirical robustness/performance on a curated subset of the hardest models
is **neutral**: no new failures and no measurable net benefit either — see
§3.7 for the full, honest results and what they mean for whether this
feature is worth keeping enabled. It remains **off by default**
(`trBroydenRecomputeInterval = 0`), so no existing user is affected either
way. §3.1–§3.5 below are the *original* plan as written before
implementation; §3.7 documents what was actually built, since it differs
from the original plan in some specifics (notably: QR-maintained Broyden
state rather than a dense matrix, `trDynamicScaling` deferred, and a
default of `0` rather than `5`).

**Goal**: convert the existing `TrustRegionSolver` into a faithful
implementation of MINPACK `hybrd` (Powell 1970). Single, well-scoped
change. No new solver slot, no new dependency.

### 3.1 New Options (in `solver.h`, after the existing TR options ~line 127)

```cpp
/** Full-Jacobian recompute period for the TrustRegion solver.
 *  0  = disabled (legacy behaviour: full J every iteration)
 *  K>0 = hybrd mode: full J every K accepted iterations; Broyden
 *        rank-1 updates in between. Default = 5.
 *  See Tier 1 of docs/solver_roadmap.md. */
int  trBroydenRecomputeInterval = 5;

/** Powell's restart criterion: force a full-Jacobian recompute after N
 *  consecutive trust-region rejections. hybrd uses 2. Default = 2. */
int  trBroydenRestartRejects    = 2;

/** Dynamic column scaling: on every full-J iter, update
 *  scale_j <- max(scale_j, ||J_unscaled[:,j]||). Matches hybrd's
 *  diag(j) = dmax1(diag(j), wa2(j)). Default = true. */
bool trDynamicScaling           = true;
```

Setting `trBroydenRecomputeInterval = 0` recovers the exact current TR
behaviour, so the change is backward-compatible.

> **What was actually implemented (§3.7): `trBroydenRecomputeInterval`
> and `trBroydenRestartRejects` defaulted to `0`/`2` respectively (not
> `5`), the Broyden state is a QR factorization rather than a dense
> matrix `B` (see §4.1), and `trDynamicScaling` was not implemented.
> The two fields that were shipped match the names/semantics above.**

### 3.2 Implementation Steps in `src/solver_trust_region.cpp`

Mirror the Broyden state pattern already proven in
`src/solver_newton.cpp:110-163`:

1. **Add Broyden state to `TrustRegionSolver::solve`**: members `B`
   (Broyden Jacobian approximation in scaled coords), `F_prev`,
   `y_prev`, `itersSinceFullJ`, `consecutiveRejects`,
   `forceFullJacobian`. Use the same names as Newton for symmetry.

2. **Replace the unconditional Jacobian evaluation** at
   `solver_trust_region.cpp:108-110` with a dual branch. The flag
   `fullJacobianThisIter` is true when *any* of:
   - `iter == 0`,
   - `itersSinceFullJ >= options.trBroydenRecomputeInterval`,
   - `forceFullJacobian`,
   - `consecutiveRejects >= options.trBroydenRestartRejects`.

   When false: call
   `problem.evaluate(xu, F, J_dummy, /*computeJacobian=*/false)`, then
   apply the Broyden Type-I update
   `B += ((dF - B·dy) * dy.transpose()) / dy.squaredNorm()` (with
   `dy = y - y_prev`, `dF = F - F_prev`), and alias `J = B`.

   > **As shipped**: this dense update was replaced with
   > `rank1QRUpdate()` against the maintained `(Q, R)` factors (§4.1),
   > per the §3.6 finding that the dense version is numerically
   > unstable. `J` is then reconstructed as `Q * R` for the existing
   > dogleg code below. See §3.7.

3. **On trust-region rejection** (`solver_trust_region.cpp:260-268`):
   increment `consecutiveRejects`; if it reaches
   `trBroydenRestartRejects`, set `forceFullJacobian = true` for the
   next iteration (Powell's restart). Reset to 0 on accept.

4. **On non-finite dogleg step** (`solver_trust_region.cpp:197`) when
   the current iteration used Broyden: set `forceFullJacobian = true`
   and retry without advancing the iterate (mirror
   `solver_newton.cpp:221-226`).

5. **On successful step**: update `F_prev = F_new`, `y_prev = y`,
   `itersSinceFullJ++`, `consecutiveRejects = 0`.

6. **Dynamic scaling** (`trDynamicScaling = true`): on every full-J
   iteration, update `scale_j <- max(scale_j, ||J_unscaled[:,j]||)`
   before constructing the scaled Jacobian
   `J = J_unscaled * scale.asDiagonal()`. Initial scaling via
   `computeScalingFactors()` is unchanged; dynamic scaling only widens
   it, never shrinks.

7. **Header comment**: add a one-paragraph reference to Powell (1970),
   *A Hybrid Method for Nonlinear Equations*, and to MINPACK `hybrd`,
   following the citation pattern in `solver_lm.cpp:1-21` (required by
   `docs/contributing.md` §4).

### 3.3 Integration Checklist (per `docs/contributing.md` §3)

Walk through these explicitly when implementing Tier 1:

- **§3.2 Solvers**: `include/coolsolve/solver.h` (3 new fields with
  defaults — zero-overhead when `trBroydenRecomputeInterval = 0`),
  `src/solver_trust_region.cpp` (the dual-branch Jacobian + Powell
  restart + dynamic scaling). No new translation unit: Tier 1 modifies
  the existing `TrustRegionSolver` in place, per the user's choice to
  keep the slot name.
- **§3.2 Config parsing**: `src/solver.cpp`
  `loadSolverOptionsFromFile()` — parse the three new keys. Add a
  console/trace log when hybrd mode is active (analogous to Newton's
  Broyden trace at `solver_newton.cpp` verbose branch).
- **§3.3 `examples/coolsolve.conf`**: add the three new keys
  **commented out** (showing the defaults), with comment blocks
  explaining what they do, when to enable hybrd mode, and the
  trade-off (K = 0 → exact legacy TR; K = 5 → 4-5× fewer CoolProp
  calls/iter on large blocks).
- **§3.4 GUI**: add three `ConfigField` entries to the existing
  "Trust Region" `ConfigGroup` in `gui/src/components/ConfigEditor.tsx`
  (~line 213). No new solver name to register, so `ALL_SOLVERS` and
  `PIPELINE_PRESETS` are unchanged.
- **§3.8 Debug mode**: when verbose/hybrd mode is on, emit a per-iter
  trace line showing `fullJacobianThisIter` and `consecutiveRejects`
  (similar to Newton's Broyden trace).
- **§3.12 README**: update §"Available Solver Algorithms" →
  "Trust-Region Dogleg" subsection (around line 559) to describe the
  hybrd mode, the Powell restart criterion, and the new default value.
- **§3.12 `docs/solver_roadmap.md`**: this very document; mark Tier 1
  as delivered once the implementation lands.

### 3.4 Tests for Tier 1

All tests follow the existing pattern in `tests/test_new_solvers.cpp`
(`NonLinearSolver::Problem` lambda + `solver.solve(problem, x, opts)`
+ `CHECK_THAT(... WithinAbs(...))`).

> **As shipped (§3.7)**: the actual tests landed in the new
> `tests/test_trust_region_hybrd.cpp` (6 cases: config defaults, Powell
> badly-scaled, restart-on-reject, mild Rosenbrock, and a
> Broyden-tridiagonal full-Jacobian-eval-count reduction check) plus 2
> new cases in `tests/test_config.cpp`, rather than exactly the T1–T8 /
> A1–A2 list below. The list below is retained as the original test
> *design* — most of its scenarios (Powell badly-scaled, restart
> behaviour, full-J eval counting, config round-trip) are covered by
> what shipped; the ones tied to `trDynamicScaling` (T8) were not
> applicable since that option was not implemented.

#### 3.4.1 Unit Tests (Catch2 tag `[trustregion][broyden]`)

Add to `tests/test_new_solvers.cpp`:

- **T1 — Powell badly-scaled** (the canonical hybrd test):
  `F1 = 1e4 * x1 * x2 - 1`, `F2 = exp(-x1) + exp(-x2) - 1.0001`.
  Start from `(0, 1)`. Assert `Success` and
  `‖x - (1.098e-5, 9.106)‖∞ < 1e-6` with `tolerance = 1e-10`.
- **T2 — Freudenstein-Roth**: start from `(11, -9.5)` (hard far-start
  where Newton + Broyden line-search diverges). Assert `Success`.
- **T3 — Rosenbrock, hard variant** (`b = 100`): start from
  `(-1.2, 1)`. Assert `Success` in ≤ 30 iterations.
- **T4 — Broyden tridiagonal** (`n = 20`): assert `Success`. Count
  `evaluate(..., computeJacobian = true)` calls; with
  `trBroydenRecomputeInterval = 5`, expect ≤ 25 full-J evals vs. ≥ 100
  with `= 0`.
- **T5 — Powell singular** (rank-deficient Jacobian at start):
  assert `Success`; verify via the verbose trace that
  `forceFullJacobian` fires after exactly
  `trBroydenRestartRejects` consecutive rejections.
- **T6 — Regression guard**: with
  `trBroydenRecomputeInterval = 0`, `trDynamicScaling = false`, the
  solver reproduces the current TR behaviour bit-for-bit on Rosenbrock
  (same iteration count ±0).
- **T7 — Restart instrumentation**: on a synthetic 4-var block whose
  dogleg step lands in the CoolProp penalty region twice in a row,
  verify via the verbose trace that the full-J recomputation happens
  on the third iteration.
- **T8 — Dynamic scaling never shrinks**: on Powell badly-scaled,
  capture the `scale` vector at every iteration and assert
  `scale_j >= scale_j_prev` for every `j` and every iteration.

#### 3.4.2 Advantage Tests (tag `[trustregion][broyden][advantage]`)

Following the pattern at `tests/test_new_solvers.cpp:210-557` with a
header comment naming the mathematical reason:

- **A1 — Powell badly-scaled**: `TrustRegion` (hybrd mode, default
  options) succeeds; `Newton` (any Broyden setting) and `LM` fail
  within `maxIterations = 50`.
- **A2 — CoolProp-penalty block**: synthetic 4-var block whose
  Newton step lands in the penalty region (mimics `cooling_tower`'s
  blk11). `TrustRegion` with `trBroydenRecomputeInterval = 5`
  converges; `TrustRegion` with `= 0` hits `MaxIterations`.

#### 3.4.3 Config Tests (in `tests/test_config.cpp`)

Mirror the `broydenRecomputeInterval` parsing test in
`tests/test_tier2_improvements.cpp:402-428`:

- Parse `trBroydenRecomputeInterval = 5` from a `coolsolve.conf`
  string and verify the field round-trips.
- Parse `trBroydenRestartRejects = 2` likewise.
- Parse `trDynamicScaling = false` likewise.
- Negative test: an out-of-range value
  (`trBroydenRecomputeInterval = -1`) emits a diagnostic via
  `DiagnosticCollector` (per `docs/contributing.md` §3.9) and falls
  back to the default.

#### 3.4.4 Robustness Coverage

The robustness suite (`tests/test_solver_robustness.cpp`) already
includes a `TrustRegion only` configuration
(`test_solver_robustness.cpp:332-454`). No change needed: the
behaviour change is automatic once `trBroydenRecomputeInterval`
defaults to 5. The regenerated
`examples/solver_robustness_report.md` will reflect the new numbers.

### 3.5 Validation Order (Strict)

Per `docs/contributing.md` §2.4, run in this order and stop on the
first failure:

1. `cmake --build build -j$(nproc)` — must compile cleanly.
2. `./build/coolsolve_tests` (default unit tests) — must pass.
3. `./build/coolsolve_tests "[trustregion]"` — new Tier 1 tests pass.
4. `./build/coolsolve_tests "[config]"` — new config-parsing tests
   pass.
5. `./build/coolsolve_tests "[examples-comprehensive]"` — must still
   reach 36/42 (no regression).
6. **Only last** (slow, ~10 min):
   `./build/coolsolve_tests "[solver-robustness]"` regenerates
   `examples/solver_robustness_report.md`.

   > **As run (§3.7)**: per instructions, the full 42-file suite was
   > not re-run for this change (too slow to iterate on). Instead,
   > `COOLSOLVE_EXAMPLES_DIR` was pointed at a curated 9-file subset of
   > the hardest models plus two easy ones, with two new configs added
   > to `test_solver_robustness.cpp` (`TrustRegion only` vs.
   > `TrustRegion + Hybrd K=5`). See §3.7 for results.

**Targets after Tier 1 lands:**

- With initials: ≥ 95% (currently 92.7%).
- Without initials: ≥ 75% (currently 67.5%; the bigger jump to ≥ 90%
  requires Tier 2.2 multi-start).
- `scroll_compressor` solve time: ≤ 0.10 s (currently 0.11 s).
- For 30-var blocks: full-J evals/iter drop from ~31 (1 + n) to ~7
  (1 every 5 iters + 4 residual-only).

> **Not validated against these exact targets (§3.7)**: the curated
> 9-file subset used for the actual validation run is a different,
> smaller sample than the full 42-file population these percentages are
> computed over, and on that subset legacy `TrustRegion only` already
> solves 100% of the with-initials cases (it happens not to include
> any of the with-initials failures from the full-suite baseline), so
> the ≥95%/≥75% targets could not be directly confirmed or refuted.
> What *was* confirmed: zero new failures on this subset, and no
> measurable net speedup (see §3.7) — so these targets remain open
> questions rather than met goals.

### 3.6 Investigation Findings (June 2026) — Why Tier 1 (Take 1) Was Reverted

An implementation attempt followed the plan above: `SolverOptions` gained
`trBroydenRecomputeInterval` (default 5), `trBroydenRestartRejects`
(default 2), `trDynamicScaling` (default true), and
`TrustRegionSolver` was given the dual-branch Jacobian (full J vs.
Broyden rank-1 update) described in §3.2.

**Symptom**: scoped robustness testing against real thermodynamic
models showed the change *increased* iteration counts and, in some
cases, turned a previously-solving model into a failure:

| Model                    | Effect of enabling Broyden (K=5) |
|---------------------------|-----------------------------------|
| `condenser_3zones`        | +40% iterations                   |
| `zorlu_heat_pump`         | +75% iterations                   |
| `air_screw_compressor`    | TrustRegion-only: `Success` → `MaxIterations` |

**Root cause**: the dense Broyden Type-I rank-1 update
(`B += ((dF - B·dy)·dyᵀ) / dyᵀdy`) is numerically unstable for the
ill-conditioned Jacobians typical of thermodynamic property blocks,
where variables span many orders of magnitude (`T ~ 300`, `P ~ 1e7`,
`h ~ 1e5`) even after scaling. Each rank-1 update accumulates error in
directions the update doesn't observe, and this error compounds over
the K iterations between full-Jacobian refreshes. MINPACK's actual
`hybrd` never forms this dense matrix at all: it maintains the **QR
factorization** of the Jacobian incrementally, via the numerically
stable Givens-rotation-based `r1updt`/`r1mpyq` routines. That
incremental-QR foundation is what §4.1 (Tier 2.1) describes — and it
had not been implemented yet, so the dense-matrix version in this
attempt was the numerically fragile approximation, not the real
algorithm.

**Verification performed**: rather than re-running the full 42-file
robustness suite, the diagnosis was confirmed by:

1. `git diff` against every solver `.cpp` file — confirmed the revert
   left `solver_trust_region.cpp`, `solver_newton.cpp`, `solver_lm.cpp`,
   and `solver_common.h` **byte-identical** to the pre-Tier-1 baseline.
2. A full default unit test run (264 test cases / 1904 assertions) —
   all passing, unaffected.
3. A scoped robustness run (`COOLSOLVE_EXAMPLES_DIR` pointed at a
   10-model subset including every model flagged above) — reproduced
   the pre-Tier-1 baseline almost exactly (identical pass/fail category
   and residual for 9/10 models across all 10 solver configurations,
   including `TrustRegion only` now solving `air_screw_compressor`
   again). The one exception was the combined *Default pipeline*
   column on `air_screw_compressor`'s block 13 (a 13-variable,
   near-singular compressor-leak block that is already the hardest
   case in the suite): it alternated between success and failure
   across separate runs. This is a pre-existing run-to-run sensitivity
   in `solveBlockSequential`'s warm-start logic (later pipeline stages
   reuse the best residual found so far, so Newton's exact — and
   CoolProp-cache-sensitive — failure trajectory changes what starting
   point TrustRegion receives), not something introduced by this
   attempt: the warm-start code itself is untouched.

**Conclusion**: at the point this investigation concluded, there was
**no active code regression** — the working tree's solver behaviour is
functionally identical to the pre-Tier-1 baseline. The regression
previously visible in `examples/solver_robustness_report.md` reflected
an intermediate, since-reverted state of the dense-Broyden experiment
and has been refreshed. The speculative `trBroydenRecomputeInterval` /
`trBroydenRestartRejects` / `trDynamicScaling` options were removed
(rather than left in as inert/dead config surface — setting them would
have silently done nothing) so the codebase reflects exactly what is
implemented.

**Corrected path forward**: implement §4.1 (incremental QR via
`r1updt`/`r1mpyq`) *first*, as a standalone, independently testable
change. Only once TrustRegion maintains a numerically stable
incremental factorization should Broyden updates be layered on top —
at that point they operate on the QR factors, not a dense matrix, and
should no longer suffer the instability documented above. In other
words, Tier 2.1 is a **prerequisite** for Tier 1, not a follow-on
optimisation; the tiers should be re-ordered in any future planning
pass.

### 3.7 Take 2 (June 2026) — QR-Based Implementation & Empirical Findings

The corrected path forward above was followed exactly: §4.1
(`solver_hybrd_qr.h`/`.cpp`, incremental QR via Golub & Van Loan §12.5
rank-1 updates — the same numerically-stable building block as MINPACK's
`r1updt`/`r1mpyq`) was implemented and unit-tested **first** (916
assertions across random-matrix agreement tests, tag `[hybrd][qr]`).
Only then was `TrustRegionSolver` given the dual-branch Jacobian.

**What was actually built, and how it differs from §3.1/§3.2's original
plan**:

- The Broyden approximation is maintained as a QR factorization (`Q`,
  `R`), updated via `rank1QRUpdate()`, **not** as a dense matrix `B`
  updated via the raw Type-I formula. This is the fix for the §3.6
  instability.
- Each iteration where Broyden mode is active still reconstructs a
  dense Jacobian `J = Q * R` and re-factors it with
  `ColPivHouseholderQR` for the actual dogleg solve, rather than doing
  a fully "hybrd-native" O(n²) triangular solve directly against `R`.
  This deliberately trades away part of the theoretical speedup (see
  below) to keep `ColPivHouseholderQR`'s rank-revealing pivoting — the
  incrementally-updated `R` has no pivoting, so a direct triangular
  solve on a near-singular Broyden-updated Jacobian could silently
  produce a large-but-finite garbage step that would not trip the
  existing `!dx_n.allFinite()` safety net.
- `trBroydenRecomputeInterval` defaults to **0** (disabled), not 5.
  Given the empirical findings below, shipping it enabled-by-default
  would not be justified.
- `trDynamicScaling` was **not implemented**. Dynamically growing the
  `scale` vector mid-solve requires rescaling the live iterate
  (`y = x/scale`) to preserve the physical-position invariant — a
  nontrivial correctness-sensitive change. Deferred; not needed to
  evaluate the core hybrd-vs-legacy question.
- Powell's restart criterion (force full-J after `trBroydenRestartRejects`
  consecutive rejections, default 2) **was** implemented as planned.

**Validation performed** (per §3.5, substituting a curated subset for
the full robustness suite as instructed — the full 42-file ×
22-configuration matrix is too slow to run on every iteration of this
work):

1. Full default unit suite: **274 test cases / 2832 assertions, all
   passing** (baseline was 264/1904 at the time of §3.6; growth is the
   new QR + hybrd tests, not a change in existing behaviour).
2. Targeted tests: `[hybrd]` tag — **936 assertions / 11 test cases**
   (the 4 QR-correctness cases from §4.1 plus 7 new TrustRegion-hybrd
   cases in `tests/test_trust_region_hybrd.cpp` and `test_config.cpp`,
   covering config defaults, Powell-badly-scaled convergence, the
   restart-on-reject safety net, mild Rosenbrock, and a
   full-Jacobian-eval-count reduction check on a Broyden-tridiagonal
   problem). `[trustregion]` tag — **38 assertions / 15 test cases**.
3. `[examples-comprehensive]`: still **36/42**, identical to the
   pre-existing baseline (expected: default pipeline uses K=0).
4. **Curated-subset robustness comparison**: a 9-file subset chosen as
   "hardest known models" (`zorlu_heat_pump` 63 vars, `condenser_3zones`
   62 vars, `heat_pump_MSTh_SB_R10` 39 vars, `orc_r245fa` 38 vars,
   `orc_co2` 28 vars, `cooling_tower` 22 vars, `air_screw_compressor`)
   plus two easy models for contrast (`cpbar`, `rankine1`), run through
   `tests/test_solver_robustness.cpp` with two new configurations
   added specifically for this comparison, `TrustRegion only` vs.
   `TrustRegion + Hybrd K=5`, both with-initials and without-initials
   (22 configs × 9 files = 198 combinations total, ≈524s wall time).

**Empirical results (the honest verdict)**:

| Configuration                              | With initials  | Without initials |
|---------------------------------------------|---------------:|------------------:|
| TrustRegion only                             | 9/9 (100.0%), avg 0.079 s | 2/9 (22.2%), avg 0.054 s |
| TrustRegion + Hybrd K=5                      | 9/9 (100.0%), avg 0.104 s | 2/9 (22.2%), avg 0.060 s |

- **No regressions**: every file that legacy `TrustRegion only` solves,
  hybrd mode also solves — both with and without initials. This is the
  key difference from the §3.6 attempt (which turned
  `air_screw_compressor` from `Success` into `MaxIterations`). Without
  initials, both configurations fail on exactly the same 7/9 files
  (`condenser_3zones`, `cooling_tower`, `cpbar`, `heat_pump_MSTh_SB_R10`,
  `orc_co2`, `orc_r245fa`, `zorlu_heat_pump`), mostly via
  `MaxIterations` in both cases.
- **No improvement either**: identical success counts in both
  scenarios. Hybrd mode does not rescue any without-initials failure
  and does not add any with-initials failure.
- **Slightly slower wall-clock, not faster**, on the two largest blocks
  in the with-initials set — the case this feature is meant to help:
  `condenser_3zones` 0.18s → 0.35s (+94%), `zorlu_heat_pump` 0.09s →
  0.13s (+44%). Smaller blocks (`cpbar`, `orc_co2`, `rankine1`, etc.)
  show no measurable difference either way.
- A dedicated unit test (`"Broyden tridiagonal reduces full-Jacobian
  evaluations"`) confirms the *mechanism* works as designed — K=4
  genuinely performs fewer full-Jacobian evaluations than K=0 on that
  synthetic problem. The wall-clock results above show that this
  reduction in Jacobian *evaluations* does not currently translate
  into a reduction in overall *solve time* on real CoolProp models.

**Root cause of the missing speedup**: per the design trade-off above,
the Broyden branch still pays the full `O(n³)` cost of reconstructing
`J = Q·R` and re-running `ColPivHouseholderQR` on every iteration —
identical linear-algebra cost to the legacy full-Jacobian path. The
only thing actually saved is the CoolProp Jacobian-evaluation call
itself. At this codebase's problem sizes (n ≤ 63) that call is
apparently not the dominant per-iteration cost, so skipping it does
not show up as a net win — and any extra iterations the Broyden
approximation needs to reach convergence (compared to a fresh exact
Jacobian every step) tip the balance further against it.

**Conclusion / recommendation**: keep the feature in the codebase —
it is fully opt-in, thoroughly tested, causes zero regressions, and
the underlying QR module (§4.1) is independently valuable — but do
**not** claim it improves robustness or speed, and do **not** enable it
by default. Realizing the originally-hoped-for speedup would require
a genuinely "hybrd-native" solve that uses the maintained `R` directly
(an `O(n²)` triangular solve, forgoing `ColPivHouseholderQR`'s
pivoting) rather than reconstructing and re-factoring a dense matrix
every iteration — deferred as future work in §5, since it reopens the
robustness-vs-speed trade-off discussed in §3.6 and would need its own
careful validation.

---

## 4. Tier 2 — Additional Improvements (after Tier 1 lands)

The four items below are independent and may be implemented in any
order. Each follows the full
`docs/contributing.md` §2 workflow (plan → code → targeted tests →
full test suite → docs).

### 4.1 Incremental QR Factorization (Port MINPACK `r1updt` + `r1mpyq`)

**Status (June 2026): DELIVERED.** Implemented *before* Tier 1, per the
corrected ordering in §3.6/§3.7 (this item is a prerequisite for Tier 1,
not a follow-on — the section numbering here is left as originally
planned to avoid renumbering the whole document, but the actual
delivery order was 4.1 then 3). Lives in
`include/coolsolve/solver_hybrd_qr.h` / `src/solver_hybrd_qr.cpp`,
exposing `computeInitialQR()` and `rank1QRUpdate()`. Validated with 916
assertions (tag `[hybrd][qr]`) checking agreement against a from-scratch
`ColPivHouseholderQR` after sequential rank-1 updates. Consumed by
`TrustRegionSolver`'s hybrd mode (§3.7); not yet consumed by
`solver_newton.cpp`'s pre-existing dense-Broyden path (that would be a
separate, independent follow-on, not required by anything above).

**Why**: When Broyden updates `B`, the Tier 1 solver recomputes
`B` and then calls `ColPivHouseholderQR` from scratch
(`solver_trust_region.cpp:195-197`) — `O(n³)` every iteration. hybrd
instead maintains the QR of `B` incrementally via the 80-line `r1updt`
routine, giving `O(n²)` per update. On 60-var blocks
(`zorlu_heat_pump`, `condenser_3zones`) this is roughly a 60× speedup
on the linear-algebra portion of each iteration.

**Note added in §3.7**: this `O(n²)` speedup is only realized if the
consumer uses the maintained `R` directly for the solve. The current
`TrustRegionSolver` integration does not do this (it reconstructs dense
`J = Q·R` and re-runs `ColPivHouseholderQR`, for pivoting/robustness
reasons) — so the linear-algebra speedup described here is not yet
observed end-to-end. The QR module itself works exactly as specified
and is validated in isolation.

**Implementation**:

- Copy `r1updt.c` and `r1mpyq.c` from cminpack (MINPACK license,
  BSD-compatible) into a new translation unit
  `src/solver_hybrd_qr.cpp` with header
  `include/coolsolve/solver_hybrd_qr.h`.
- Expose `updateQR(MatrixXd& R, MatrixXd& Q, const VectorXd& u, const
  VectorXd& v)` updating the thin QR factorization in place by the
  rank-1 update `R + u·vᵀ`.
- Inside `TrustRegionSolver`, when Broyden mode is active, replace
  the `ColPivHouseholderQR` recompute with cached `R`, `Q` updated
  via `updateQR`. On `forceFullJacobian`, refresh `R, Q` from the new
  `B`.

**Tests (tag `[hybrd][qr]`)**:

- Random-matrix agreement: for 50 random `(n × n)` matrices with `n
  ∈ {2, 5, 10, 30}`, perform 10 sequential rank-1 updates with random
  `u, v`. Assert `‖R_updated − chol(B_updated)‖∞ < 1e-12` against a
  from-scratch QR.
- Speed benchmark on `condenser_3zones` (62-var block): with
  incremental QR enabled, the linear-algebra time per iteration must
  drop by ≥ 3× at parity of convergence (same iteration count, same
  solution).

### 4.2 Multi-Start for the Without-Initials Case

**Status (July 2026): DELIVERED.** Implemented as a block-level multi-start
fallback in `Solver::solveBlockWithMultiStart()` (`src/solver.cpp`), gated by
`multiStartEnabled` (default `true`) and `multiStartMaxRestarts` (default 4).
Both `SolverOptions` keys are parsed from `coolsolve.conf`, exposed in the GUI
ConfigEditor, and covered by `tests/test_multistart.cpp` + `test_config.cpp`.

**Empirical results (curated 44-file subset, without initials, full default
pipeline, 30 s per-solve timeout):**

| Configuration                | Without initials |
|------------------------------|-----------------:|
| Default pipeline (baseline)  | 27/44            |
| Default pipeline + multi-start | 30/44          |

- **3 models rescued** from FAIL → OK: `piston_compressor` (blk 4, scale
  candidate ×0.1), `refrigeration_compressor` (blk 4, scale candidate ×0.1),
  and `orc_r245fa` (blk 12, a thermo block rescued by the medium-pressure
  regime candidate).
- **Zero regressions**: every model the baseline solves, multi-start also
  solves (with initials, multi-start never engages because all blocks converge
  first try — verified by `[examples-comprehensive]` still reaching 36/42 with
  identical timings).
- The roadmap's original target of "+4 to +6 models" was only partially met
  (+3). The remaining without-initials failures (`orc_co2`, `cooling_tower`,
  `heat_pump_MSTh_SB_R10`, `zorlu_heat_pump`, …) have huge initial residuals
  (|F| up to 1e12) from deeply coupled 5th-order polynomial fits and CoolProp
  calls whose convergence basin is too narrow for a handful of starting points
  — these need good `.initials` or KINSOL (Tier 3).

**What was actually built, and how it differs from the original plan:**

- The fallback triggers per-block on **any non-fatal pipeline failure**
  (MaxIterations / SingularJacobian / LineSearchFailed / Diverged), not only
  when no `.initials` is loaded. This gives the intended "rescue the
  without-initials case" behaviour automatically (those are exactly the blocks
  that fail) while keeping **zero overhead when every block converges** — the
  search only runs after a failure. Size-1 blocks are skipped (Newton1D already
  does its own multi-probe).
- Two complementary candidate strategies are generated by
  `Solver::generateMultiStartCandidates()`, picked automatically from the
  block's variables:
  1. **Thermo blocks** (containing ≥1 variable with `inferredProperty` /
     `inferredFluid`): each thermo variable is re-evaluated via the new
     `computeThermoGuessAt()` helper (`src/variable_inference.cpp`) at a
     coherent reference `(T, P)` state, so the candidate is internally
     consistent (`h` matches `enthalpy(fluid, T, P)`). Six pressure regimes
     cover the low/medium/high bands of refrigeration and ORC loops. This is
     the principled fix for the "wrong pressure level" failure mode — an
     earlier, inconsistent version that independently snapped variables
     produced residuals of ~1e12 because `h` and `(T, P)` did not match.
  2. **Purely-algebraic blocks** (no thermo variables): the default guess is
     scaled by `{0.1, 0.3, 3, 10}`. This is the high-value move for blocks
     whose solution is at a different *magnitude* than the default 1.0
     (e.g. `piston_compressor`'s volumetric block: clearance `C ≈ 0.05` and
     swept volume `V_dot_s ≈ 0.1` guessed as 1.0).
- Candidates are deterministic (no randomness) so debug output is diffable, in
  line with `docs/contributing.md` §3.8. The winning candidate is reflected in
  the block trace (`solverType` prefixed `MultiStart(...)->`) and an info
  diagnostic `V006` is emitted.

**Key implementation note (bug found and fixed during development):** the
original guess must be captured **before** the first `solveBlock` attempt.
`solveBlock` leaves the evaluator state at its (possibly diverged) iterate
after a failure, so reading the "default" guess afterwards returns that
garbage — scaling garbage by ×0.1 still gives garbage. The fix captures the
initial vector up-front and restores it before candidate generation.

**Why**: The biggest robustness gap is the **without-initials** case
(67.5%). When no `.initials` file is provided, the pipeline sees a
single bad starting point and gives up. CoolSolve already knows each
variable's physical kind (temperature, pressure, dimensionless, ...) from
unit annotations — this information is currently unused for
initialisation.

**Implementation** (original plan, retained as the design record):

- In `src/solver.cpp` `Solver::solve()` (around line 2079), if no
  `.initials` file was loaded, generate a candidate set of starting
  vectors:
  1. All-ones (default for dimensionless variables).
  2. Physical defaults from unit tags: `T = 300 K`, `P = 1e5 Pa`,
     `h = 1e5 J/kg`, `m_dot = 1.0 kg/s`, ...
  3. ±50% perturbations of (1) and (2) for the first 2-4 hardest
     blocks (largest SCCs).
- Try each candidate on the largest block only; keep the best
  residual. Smaller blocks keep their default starting point.
- Cap at `multiStartMaxRestarts = 4` candidates (configurable).

**New options** (in `solver.h`):

```cpp
bool multiStartEnabled       = true;  // active only when no .initials loaded
int  multiStartMaxRestarts   = 4;     // candidate count for the hardest block
```

**Tests (tag `[multistart]`)**:

- Targeted unit tests: synthesize a 4-var block with a known solution
  far from `(1, 1, 1, 1)`. Without multi-start, the default pipeline
  fails. With multi-start, it converges.
- Without-initials robustness on `orc_co2`, `heat_pump_MSTh_SB_R10`,
  `cooling_tower`, `piston_compressor`: each must transition from
  FAIL to OK with multi-start enabled.
- Target: +4 to +6 models solved without initials (67.5% → 80-85%).

### 4.3 Dynamic Column Scaling in Newton and LM

Tier 1 already adds hybrd-style `diag ← max(diag, ‖J[:,j]‖)` inside
`TrustRegion`. To extend the same rule to the other two gradient
solvers:

**Implementation**:

- Factor the scaling update out of `solver_trust_region.cpp` into
  `solver_common.h` as
  `void updateScalingFactors(VectorXd& scale, const MatrixXd&
  J_unscaled)`.
- Call from `solver_newton.cpp` and `solver_lm.cpp` on every full-J
  iteration (gated by the existing `trDynamicScaling` option, renamed
  to `dynamicScaling` since it now applies to all three solvers).

**Tests (tag `[scaling][dynamic]`)**:

- Unit test: verify `scale_j` is monotonically non-decreasing across
  iterations on Powell badly-scaled, for each of Newton/TR/LM.
- Regression test: with `dynamicScaling = false`, the three solvers
  reproduce their current iteration counts bit-for-bit.

### 4.4 Pseudo-Arclength Continuation in Homotopy

**Why**: The current Homotopy solver uses predictor-corrector with
adaptive step in the homotopy parameter `t`. It cannot pass turning
points (where the path folds back in `t`), which is a known failure
mode for hard thermodynamic models with non-monotonic `t`-paths.

**Implementation**:

- In `src/solver_homotopy.cpp`, replace the `t`-stepping with the
  arclength parameterization `(s, x(s), t(s))` of Keller (1977). Solve
  the augmented system
  `[H(x, t); ẋ·Δx + ṫ·Δt − Δs] = 0` using the Tier 1 hybrd-style
  TrustRegion as the corrector (Tier 1 must land first).
- Add an option `homotopyArclength` (default `true`) so the current
  `t`-stepping behaviour can be recovered for debugging.

**Tests (tag `[homotopy][arclength]`)**:

- Synthetic fold test: continue `x² + y² + t·x·y = 1` from
  `(1, 0, 0)` to `(0, 1, 1)` (path has a fold in `t`). Current
  Homotopy fails; pseudo-arclength passes.
- Regression: on all currently-solved Homotopy-only cases, the new
  arclength mode converges in ≤ the same number of iterations.

---

## 5. Tier 3 — Future Work (Not in This Round)

| #  | Improvement                       | Primary benefit                                              | Estimated effort |
|----|-----------------------------------|--------------------------------------------------------------|------------------|
| 13 | **KINSOL (SUNDIALS) integration** | Robustness for very large blocks (> 30 vars) and preconditioning | 1-2 weeks        |
| 14 | **Hybrd-native `O(n²)` triangular solve for TrustRegion Broyden mode** | Realize the speedup §3.7 found missing: solve directly against the maintained `R` instead of reconstructing dense `J` and re-running `ColPivHouseholderQR` every iteration | 3-5 days, plus a careful robustness re-validation (loses rank-revealing pivoting, reopening §3.6's instability concern) |

KINSOL remains the recommended path when even Tier 1 + Tier 2 cannot
crack the largest industrial models (e.g. `orc_solar_complex`).

### Not Recommended

| Improvement                              | Why skip                                                                                                                                                                                       |
|------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Anderson acceleration**                | Does not use the Jacobian → slower than Newton when Newton works, and the pipeline already has 6 fallback strategies. The "no Jacobian needed" advantage is addressed better by BisectionND. Broyden (§3.2.5 of the old roadmap, now §3 of this document) is the better quasi-Newton approach since it starts from a true Jacobian. |
| **Code complexity monitor**              | The codebase is well-structured after the file split. LOC tracking is noise, not signal.                                                                                                       |
| **Nonlinear preconditioning (residual weighting)** | Variable scaling (already implemented and extended in Tier 2.3) addresses the same problem. Residual weighting adds complexity for marginal benefit given that analytical derivatives also improve conditioning. |
| **Better initial guess via simplified model** | Requires model-specific knowledge and a meta-solving infrastructure. The pipeline + homotopy + Tier 2.2 multi-start already handle poor initial guesses well. Document manual workflows instead. |

---

## 6. What Has Been Implemented (Status Reference)

This section documents the layers that the Tier 1 / Tier 2 work builds
on. It is preserved from earlier revisions of this document as a
reference; nothing here is open work.

### 6.1 Solver Algorithms

| Feature                                    | Status   | Notes                                                                                                                  |
|--------------------------------------------|----------|------------------------------------------------------------------------------------------------------------------------|
| Newton + non-monotone line search          | **Done** | Default workhorse; Grippo et al. 1986 non-monotone Armijo (M = 10); Broyden quasi-Newton rank-1 updates                |
| Trust Region Dogleg                        | **Done** | Adaptive initial radius (Cauchy step norm), smoother rho-based radius adaptation, gradient-based rejection recovery; optional hybrd-style Broyden/QR reuse (`trBroydenRecomputeInterval`, off by default — see §3.7) |
| Levenberg-Marquardt                        | **Done** | Nielsen's λ adaptation, cumulative Marquardt diagonal scaling, geodesic acceleration (Transtrum & Sethna 2012)         |
| BisectionND (simplex bisection)            | **Done** | Sign-pattern simplex; `bisectionNDMaxBlockSize` (default 8), `bisectionNDIterFactor`                                   |
| Homotopy continuation                      | **Done** | Predictor-corrector with adaptive step, Newton + LM corrector fallback, final polish                                   |
| Partitioned solver                         | **Done** | Per-variable diagonal updates; gracefully returns `MaxIterations` for small blocks                                     |
| Newton1D                                   | **Done** | Specialized for size-1 blocks: multi-probe + bisection + Newton hybrid (extracted to `solver_newton1d.cpp`)            |
| Tearing                                    | **Done** | Greedy FVS + acyclic sequential solve + Newton on tear residuals                                                       |
| Explicit solve (size-1 blocks)             | **Done** | `tryExplicitSolve()` detects `x = expr(known)` pattern, evaluates RHS directly                                        |
| Symbolic block reduction                   | **Done** | CoolProp call inversion + explicit extraction + equation substitution; auto re-decomposition into sub-blocks           |
| Sequential pipeline                        | **Done** | Fallback chain with multi-round restart (up to 10 rounds)                                                              |
| Parallel pipeline                          | **Done** | First-to-solve-wins via `std::async`, shared stop flag                                                                 |
| Multi-start fallback                       | **Done** | Per-block retry from CoolProp-consistent pressure regimes (thermo blocks) or scale factors (algebraic blocks); `multiStartEnabled`, `multiStartMaxRestarts` (roadmap §4.2) |
| Cancel token propagation                   | **Done** | Checked inside iteration loops of all 6 solvers + parallel stop signal                                                 |

### 6.2 Performance and Profiling

| Feature                  | Status   | Notes                                                                                            |
|--------------------------|----------|--------------------------------------------------------------------------------------------------|
| Solve time tracker       | **Done** | `PipelineTiming` in `runner.h`; displayed in GUI                                                 |
| CoolProp call profiling  | **Done** | Per-call timing for `PropsSI`, `AbstractState`, and `HAPropsSI`; warmup time tracked separately |
| `--no-superancillary` flag | **Done** | Avoids loading the cached superancillary data                                                  |
| Detailed profiling       | **Done** | First CoolProp call triggers a warmup for fluid init; analytical derivative call counts tracked |
| Variable scaling         | **Done** | Shared `computeScalingFactors()` in `solver_common.h`, used by Newton/TR/LM/BisectionND         |

### 6.3 Architecture

| Feature                            | Status   | Notes                                                                                  |
|------------------------------------|----------|----------------------------------------------------------------------------------------|
| Solver file splitting              | **Done** | 6 separate solver files + orchestrator                                                 |
| Shared scaling utility             | **Done** | `solver_common.h` (each solver class has a thin wrapper)                               |
| Partitioned InvalidInput fix       | **Done** | Returns `MaxIterations` for blocks below min size                                      |
| GUI config editor                  | **Done** | Pipeline dropdown (8 presets), editable solver list, BisectionND params                |
| Configurable BisectionND params    | **Done** | `bisectionNDMaxBlockSize`, `bisectionNDIterFactor` in config and GUI                   |

### 6.4 Detailed Implementation Notes

The sections below describe implementation choices that are still
relevant as background. They are preserved verbatim from earlier
revisions of this document.

#### 6.4.1 Symbolic Block Reduction + Re-Decomposition

Fully implemented in `src/solver_symbolic.cpp` +
`include/coolsolve/symbolic_reduction.h` +
`src/structural_analysis.cpp`. Controlled by `enableSymbolicReduction`
(default `false`). Re-decomposition is automatic when symbolic
reduction is active.

Three techniques are applied iteratively until a fixed point:

1. **CoolProp call inversion**: if the model has
   `h = enthalpy(water, T=T, P=P)` and `h` is the matched output but
   one named-arg input is an unknown block variable while the other is
   external, reformulate as `T = temperature(water, H=h, P=P)`.
   CoolProp supports all standard input pairs (`HmassP_INPUTS`,
   `PSmass_INPUTS`, ...).
2. **Explicit extraction**: equations where all RHS variables are
   external or already-reduced are directly evaluated, removing their
   output variable from the block.
3. **Equation substitution**: variables that appear in zero other
   block equations after extraction are removed along with their
   defining equation.

After reduction, **post-reduction re-decomposition** (Tarjan SCC on the
reduced dependency graph) splits the remaining monolithic block into
independent sub-blocks that are solved sequentially. Examples:

- `condenser_3zones`: 62-var block → reduced to 56 → re-decomposed into
  13 sub-blocks (largest: 30, 15).
- `heat_pump_MSTh_SB_R10`: 39-var block → reduced to 34 → 19 sub-blocks
  (largest: 8, 7).
- `orc_co2`: 28-var block → reduced to 25 → 13 sub-blocks (largest:
  8, 5).
- `air_screw_compressor`: 13-var block → reduced to 6 → fully
  decomposed into 6 scalar sub-blocks.

#### 6.4.2 Non-Monotone Line Search

Fully implemented in all three gradient-based solvers (Newton,
TrustRegion, LM). Controlled by `lsNonMonotoneMemory` (default `10`;
set to `1` for classic monotone behaviour).

Replaced the monotone Armijo condition with the Grippo-Lampariello-
Lucidi (1986) non-monotone variant. Instead of requiring
`φ(x_{k+1}) < φ(x_k)`, the solver compares against
`max(φ(x_{k-M+1}), ..., φ(x_k))` where `M` is the memory parameter.
This allows the solver to temporarily accept larger merit values to
escape narrow curved valleys and saddle points, while maintaining the
same global convergence guarantees.

**Bounded non-monotone reference**: the raw `max(history)` can be
arbitrarily larger than the current merit when a single early bad
iterate inflates the window. To prevent this, all solvers use
`boundedRef(phi, R = 10)` which caps the reference at
`min(max(history), R × currentPhi)`.

Implementation details:

- `NonMonotoneHistory` helper struct in `solver_common.h`: circular
  buffer tracking last `M` merit values, with `boundedRef()` for safe
  non-monotone reference.
- **Newton**: `lineSearch()` takes a `refPhi` parameter (bounded max of
  history); Armijo condition uses `refPhi` instead of current φ. The
  directional derivative `∇φ·d` is still computed from the current
  iterate.
- **TrustRegion**: step acceptance uses `phi_new < refPhi` instead of
  `actual > 0`. Gain ratio ρ for radius management still uses current
  φ.
- **LM**: step acceptance uses `phi_new < refPhi` instead of
  `actual > 0`. Gain ratio ρ for lambda management still uses current
  φ.
- When `M = 1`, all three solvers degenerate to exact monotone
  behaviour (no change from previous implementation).

#### 6.4.3 Trust Region Improvements (Pre-Tier 1)

- **Adaptive initial radius** (`trAdaptiveRadius`, default `true`):
  `delta = min(trInitialRadius, max(cauchyNorm, 1.0))`.
- **Smoother rho-based radius adaptation**: replaces binary
  grow/no-grow with graduated rho thresholds (grow at ρ > 0.75, hold
  at ρ > 0.5, shrink proportionally otherwise using
  `max(0.1, 1 − ρ) × delta`).
- **Gradient-based recovery**: after many consecutive rejections
  (> 15), the solver resets to a gradient-direction step with a small
  radius, escaping degenerate trust regions. **Tier 1 lowers this
  threshold effectively to 2 by forcing a full-Jacobian recompute on
  Powell's criterion**, which usually resolves the stall before the
  gradient-direction fallback fires.

#### 6.4.4 Levenberg-Marquardt Improvements

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

#### 6.4.5 Broyden Quasi-Newton Updates (Newton Solver Only, Pre-Tier 1)

Implemented as a hybrid Newton-Broyden mode in `src/solver_newton.cpp`.
**Note**: this section refers to the *Newton* solver's Broyden mode
only. The Tier 1 work in §3 of this document extends Broyden to the
*TrustRegion* solver in a fundamentally different way (combined with
trust-region acceptance, the `scipy.fsolve` algorithm).

Broyden's method is a quasi-Newton approach that avoids recomputing the
full Jacobian at every iteration by maintaining rank-1 updates:

```
B_{k+1} = B_k + (ΔF - B_k Δx) Δx^T / (Δx^T Δx)     [Broyden Type I]
```

- **`broydenRecomputeInterval`** option (default `0` = always full
  Jacobian; `K > 0` = recompute every `K` iterations). Recommended
  starting point: `K = 5`.
- Computes true Jacobian on iteration 0, then uses Broyden Type-I
  rank-1 updates for intermediate iterations.
- **Automatic fallback**: forces full Jacobian recomputation when
  (a) Broyden step fails line search, (b) residual stalls for
  ≥ 3 iterations, or (c) scheduled recompute interval reached.
- Primarily an **efficiency** improvement (fewer CoolProp calls on
  large blocks) rather than robustness — see §3 of this document for
  the *combined* TR + Broyden mode that *does* improve robustness.

**Cost savings example**: for a 30-var block, full Newton costs ~31
CoolProp calls/iteration (30 for Jacobian + 1 for residual). With
Broyden (`K = 5`), on average only ~7 calls/iteration (1 full Jacobian
every 5 iters + 4 residual-only iters).

References:

- Broyden (1965): "A class of methods for solving nonlinear
  simultaneous equations".
- Dennis & Moré (1977): "Quasi-Newton methods, motivation and theory".
- Kelley (2003): *Solving Nonlinear Equations with Newton's Method*
  (SIAM), Chapter 4.
- Powell (1970): *A Hybrid Method for Nonlinear Equations* — the
  algorithm that combines Broyden with trust region, implemented in
  Tier 1.
- Moré, Garbow, Hillstrom (1980): MINPACK `hybrd` source code
  (`http://www.netlib.org/minpack/hybrd.f`).

---

## 7. Key Metrics to Track

| Metric                                  | Previous    | Current      | Tier 1 target        | Tier 2 target                  |
|-----------------------------------------|-------------|--------------|----------------------|--------------------------------|
| Default pipeline, with initials         | 17/18 (94%) | 36/42 (86%)  | ≥ 40/42 (95%)        | ≥ 41/42 (98%)                  |
| Default pipeline, without initials      | 14/17 (82%) | 28/42 (67%)  | ≥ 31/42 (75%)        | ≥ 38/42 (90%) with multi-start |
| Unit tests                              | 117 (805)   | 262 (1896)   | +10 (Tier 1.4)       | +15 (Tier 2)                   |
| CoolProp evals per property             | 5 (1+4 FD)  | 1-3 (AbsState + check) | 1-3 (unchanged) | 1-3 (unchanged)               |
| CoolProp evals per iter, 30-var block   | ~31         | ~31          | ~7 (hybrd mode)      | ~7 with `O(n²)` QR             |
| `scroll_compressor` solve time          | 0.42 s      | 0.11 s       | ≤ 0.10 s             | ≤ 0.08 s                       |

> **Note**: the model count grew from 18 to 42 across versions, so
> percentages are not directly comparable across the *Previous* and
> *Current* columns. The Tier 1 and Tier 2 targets are anchored to the
> current 42-model baseline. Tier 1 (this table's target column) has
> **not** been reached yet — see §3.6 for why the first attempt was
> reverted and what is required before re-attempting it.

---

## 8. Implementation Notes

### 8.1 CoolProp Input Pair Reference

For CoolProp call inversion (§6.4.1), these are the supported
`AbstractState::update()` input pairs:

| Input pair            | Inputs                | Use case                              |
|-----------------------|-----------------------|---------------------------------------|
| `PT_INPUTS`           | Pressure, Temperature | Standard forward evaluation           |
| `HmassP_INPUTS`       | Enthalpy, Pressure    | When `h` is known, solve for `T`      |
| `PSmass_INPUTS`       | Pressure, Entropy     | Isentropic processes                  |
| `HmassSmass_INPUTS`   | Enthalpy, Entropy     | When both `h` and `s` are known       |
| `DmassP_INPUTS`       | Density, Pressure     | When `ρ` is known                     |
| `QT_INPUTS`           | Quality, Temperature  | Saturation properties                 |

By detecting which variables are known vs unknown in a block, the
evaluator can choose the best input pair to minimise the number of
unknowns that need iterative solving.

### 8.2 AbstractState Caching Strategy

```
Thread-local cache: thread_local AbstractStateCache g_abstractStateCache
  - Key: (backend, canonical_fluid_name), e.g. ("HEOS", "Water")
  - Create on first use: AbstractState::factory(backend, fluid)
  - Reuse for all subsequent calls on the same thread
  - Thread safety: thread_local storage gives each thread its own cache
  - Automatic fallback: if AbstractState throws, falls back to PropsSI
```

The thread-local approach (rather than mutex-protected shared cache) was
chosen because `AbstractState::update()` is not thread-safe and cloning
`AbstractState` is expensive. Each solver thread in parallel mode gets
its own cache automatically.

### 8.3 Analytical Derivative Consistency Check

CoolProp's `first_partial_deriv()` can return incorrect derivatives
near phase boundaries for pseudo-pure fluids (mixtures treated as
pure, like Air = N₂/O₂/Ar). The issue was identified during debugging
of the `air_screw_compressor` model regression:

| State                                | Property | Analytical | FD       | Relative Error |
|--------------------------------------|----------|------------|----------|----------------|
| Air, P = 1 bar, T ≈ 81 K (near sat)  | dT/dS    | 0.0747     | 0.0011   | **6693 %**     |
| Air, P = 1 bar, T ≈ 81 K (near sat)  | dT/dP    | 2.29e-4    | 8.67e-5  | **164 %**      |
| Air, P = 1 bar, T ≈ 81 K (near sat)  | dD/dS    | -4.74e-3   | -1.81e-3 | **162 %**      |
| Air, P = 1 bar, T = 293 K (normal)   | dT/dP    | 8.32e-4    | 8.32e-4  | < 1e-9         |

The consistency check (5 % tolerance on forward-FD comparison) adds 2
extra `state->update()` calls per property but catches these errors
before they corrupt the Jacobian. For pure fluids away from phase
boundaries, the check passes and analytical derivatives are used
(accurate to better than `1e-9` relative error).

### 8.4 References

**Tier 1 (hybrd-style TR + Broyden)**:

- Powell, M. J. D. (1970): *A Hybrid Method for Nonlinear Equations*.
  In P. Rabinowitz (ed.), *Numerical Methods for Nonlinear Algebraic
  Equations*, Gordon and Breach.
- Moré, J. J., Garbow, B. S., Hillstrom, K. E. (1980): *User Guide for
  MINPACK-1* (ANL-80-74), Argonne National Laboratory. The `hybrd`
  source is at `http://www.netlib.org/minpack/hybrd.f`.
- Devernay, F.: *C/C++ Minpack* (`http://devernay.github.io/cminpack`),
  a BSD-licensed C/C++ translation of MINPACK used as the reference
  for the Tier 1 algorithm (not linked in; consulted for the
  algorithm).

**Tier 2 (incremental QR, multi-start, scaling, arclength)**:

- Golub, G. H., Van Loan, C. F. (2013): *Matrix Computations*, §12.5
  for rank-1 QR updates.
- Keller, H. B. (1977): *Numerical Solution of Bifurcation and
  Nonlinear Eigenvalue Problems*, in *Applications of Bifurcation
  Theory*, Academic Press. Origin of pseudo-arclength continuation.
- Allgower, E. L., Georg, K. (2003): *Introduction to Numerical
  Continuation Methods*, SIAM.

**General (preserved from earlier revisions)**:

- Grippo, Lampariello, Lucidi (1986): original non-monotone line
  search.
- Zhang, Hager (2004): improved variant requiring *average* decrease
  rather than *max* reference.
- Madsen, Nielsen, Tingleff (2004): *Methods for Non-Linear Least
  Squares Problems* (2nd ed.), Informatics and Mathematical
  Modelling, DTU. Source of Nielsen's λ update.
- Transtrum, Sethna (2012): *Geodesic acceleration and the small-data
  limit of nonlinear least squares*, PNAS.
- Broyden (1965): *A class of methods for solving nonlinear
  simultaneous equations*.
- Dennis, Moré (1977): *Quasi-Newton methods, motivation and theory*.
- Kelley (2003): *Solving Nonlinear Equations with Newton's Method*
  (SIAM).
