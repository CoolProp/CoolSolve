# CoolSolve — Solver Roadmap & Status

*Last updated: March 2026*

This document consolidates the performance, robustness, and architecture plans for CoolSolve's solver subsystem. It replaces the earlier `performance_plan.md`, `solver_fix_strategy.md`, and `solver_robustness_diagnosis.md`.

---

## 1. Architecture Overview

### 1.1 Solver Pipeline

CoolSolve uses a configurable multi-solver pipeline. The default (sequential mode) tries each strategy in order until one succeeds:

**Newton → TrustRegion → LevenbergMarquardt → BisectionND → Homotopy → Partitioned**

A parallel mode (first-to-solve-wins via `std::async` threads) is also available. BisectionND is automatically skipped for blocks larger than `bisectionNDMaxBlockSize` (default 8).

When `enableSymbolicReduction` is enabled, each multi-variable block is first symbolically reduced (CoolProp call inversion + explicit extraction + equation substitution), then automatically re-decomposed into independent sub-blocks before entering the solver pipeline.

### 1.2 File Structure

| File | Lines | Content |
|------|------:|---------|
| `solver.h` | 822 | All solver class declarations, `SolverOptions`, enums |
| `solver_common.h` | ~100 | Shared `computeScalingFactors()` utility + `NonMonotoneHistory` helper |
| `solver.cpp` | ~2640 | Orchestrator, Newton1D, tearing, partitioned, pipeline, symbolic reduction integration |
| `solver_newton.cpp` | ~217 | Newton + non-monotone Armijo backtracking line search |
| `solver_trust_region.cpp` | 238 | Trust Region Dogleg |
| `solver_lm.cpp` | 177 | Levenberg-Marquardt |
| `solver_bisection_nd.cpp` | 450 | N-dimensional bisection (sign-pattern simplex) |
| `solver_homotopy.cpp` | 209 | Homotopy continuation (predictor-corrector) |
| `solver_symbolic.cpp` | 611 | Symbolic block reduction (inversion + extraction + substitution) |
| `symbolic_reduction.h` | 147 | Symbolic reduction data structures and API |
| `structural_analysis.h` | 155 | Graph decomposition API (Tarjan SCC + redecomposition) |
| `structural_analysis.cpp` | 763 | Tarjan SCC, Dulmage-Mendelsohn, block re-decomposition |

### 1.3 Test Results (current baseline)

149 unit tests (1111 assertions), all passing. Models include:

| Model | Eqs | Blocks | Solve time | Status |
|-------|----:|-------:|-----------:|--------|
| condenser_3zones | 99 | 50 | 0.01 s | OK |
| cooling_tower2 | — | — | 4.86 s | OK (hard) |
| heat_pump_MSTh_SB_R10 | — | — | 0.01 s | OK |
| orc_co2 | 139 | 112 | 0.02 s | OK |
| orc_complex | — | — | — | FAIL (parse: unsupported `MODULE`) |
| orc_extraction | 133 | 113 | 0.02 s | OK |
| orc_r245fa | 173 | 151 | 0.04 s | OK |
| orc_simple | 172 | 150 | 0.03 s | OK |
| scroll_compressor | 99 | 66 | 0.11 s | OK |
| turbocompressor_interpolate | — | — | — | FAIL (evaluation error) |
| water_libr | — | — | — | FAIL (unsupported function) |
| storage_integraltable | — | — | — | FAIL (structural analysis) |

35/38 parseable models solve with the default pipeline + initials.

### 1.4 Robustness Results Summary

| Configuration | With initials | Without initials |
|---|---:|---:|
| Default pipeline (6 solvers) | 35/38 (92.1%) | 27/37 (73.0%) |
| Newton only | 30/38 (78.9%) | 20/37 (54.1%) |
| TrustRegion only | 33/38 (86.8%) | 19/37 (51.4%) |
| LevenbergMarquardt only | 30/38 (78.9%) | 17/37 (45.9%) |
| Homotopy only | 31/38 (81.6%) | 23/37 (62.2%) |
| BisectionND only | 15/38 (39.5%) | 13/37 (35.1%) |
| Partitioned only | 17/38 (44.7%) | 14/37 (37.8%) |
| Default + Tearing | 35/38 (92.1%) | 28/37 (75.7%) |
| Default + SymbolicReduction | 35/38 (92.1%) | 27/37 (73.0%) |
| Default + Tearing + SymbolicReduction | 35/38 (92.1%) | 28/37 (75.7%) |

39 example models tested (38 parseable, 37 without-initials testable). The hardest models (scroll_compressor, cooling_tower2, orc_co2) share: large algebraic blocks (28–62 vars), CoolProp evaluations near phase boundaries, and near-singular Jacobians.

**Models that fail without initials** (default pipeline): cooling_coil (blk12), cooling_tower (blk11), heat_pump (blk39), internal_combustion_engine_cpbar (blk5), orc_co2 (blk28), piston_compressor (blk4), and ~4 others. Most failures are "Max iterations" — the solver runs but gets stuck, suggesting line search or step acceptance improvements would help.

---

## 2. What Has Been Implemented

### 2.1 Solver Algorithms

| Feature | Status | Notes |
|---------|--------|-------|
| Newton + non-monotone line search | **Done** | Default workhorse; Grippo et al. 1986 non-monotone Armijo (M=10) |
| Trust Region Dogleg | **Done** | Full implementation with configurable radius, eta, shrink/grow factors |
| Levenberg-Marquardt | **Done** | Damped least-squares with diagonal scaling |
| BisectionND (simplex bisection) | **Done** | Sign-pattern simplex; `bisectionNDMaxBlockSize` (default 8), `bisectionNDIterFactor` |
| Homotopy continuation | **Done** | Predictor-corrector with adaptive step, Newton+LM corrector fallback, final polish |
| Partitioned solver | **Done** | Per-variable diagonal updates; gracefully returns `MaxIterations` for small blocks |
| Newton1D | **Done** | Specialized for size-1 blocks: multi-probe + bisection + Newton hybrid |
| Tearing | **Done** | Greedy FVS + acyclic sequential solve + Newton on tear residuals |
| Explicit solve (size-1 blocks) | **Done** | `tryExplicitSolve()` detects `x = expr(known)` pattern, evaluates RHS directly |
| Symbolic block reduction | **Done** | CoolProp call inversion + explicit extraction + equation substitution; auto re-decomposition into sub-blocks |
| Sequential pipeline | **Done** | Fallback chain with multi-round restart (up to 10 rounds) |
| Parallel pipeline | **Done** | First-to-solve-wins via `std::async`, shared stop flag |
| Cancel token propagation | **Done** | Checked inside iteration loops of all 6 solvers + parallel stop signal |

### 2.2 Performance & Profiling

| Feature | Status | Notes |
|---------|--------|-------|
| Solve time tracker | **Done** | `PipelineTiming` in runner.h; displayed in GUI |
| CoolProp call profiling | **Done** | Per-call timing for PropsSI, AbstractState, and HAPropsSI; warmup time tracked separately |
| `--no-superancillary` flag | **Done** | Avoids loading the cached superancillary data |
| Detailed profiling | **Done** | First CoolProp call triggers a warmup for fluid init; analytical derivative call counts tracked |
| Variable scaling | **Done** | Shared `computeScalingFactors()` in `solver_common.h`, used by Newton/TR/LM/BisectionND |

### 2.3 Architecture

| Feature | Status | Notes |
|---------|--------|-------|
| Solver file splitting | **Done** | 6 separate solver files + orchestrator |
| Shared scaling utility | **Done** | `solver_common.h` (each solver class has a thin wrapper) |
| Partitioned InvalidInput fix | **Done** | Returns `MaxIterations` for blocks below min size |
| GUI config editor | **Done** | Pipeline dropdown (8 presets), editable solver list, BisectionND params |
| Configurable BisectionND params | **Done** | `bisectionNDMaxBlockSize`, `bisectionNDIterFactor` in config and GUI |

---

## 3. Implemented & Remaining Improvements

### 3.1 High Impact — CoolProp Integration ✅

All four CoolProp integration items have been implemented. CoolProp calls previously dominated solve time, with each property evaluation requiring **5 PropsSI calls** (1 for the value + 4 for finite-difference Jacobian). The new implementation reduces this to **1–3 AbstractState evaluations** per property (1 for value + 2 for forward-FD consistency check when analytical derivatives are enabled, or 1 for value only in residual-only mode).

#### 3.1.1 Low-Level CoolProp API (`AbstractState`) — **Done**

Replaced `PropsSI` calls with the low-level `AbstractState` API:
- Thread-local `AbstractStateCache` maps `(backend, fluid)` → cached `AbstractState` objects
- Property evaluation via `state->update(input_pair, val1, val2)` + `state->hmass()` etc.
- Automatic fallback to `PropsSI` if `AbstractState` throws (e.g., unsupported fluid)
- Configurable via `coolpropUseAbstractState` (default: `true`) and `coolpropBackend` (default: `HEOS`)

**Files changed**: `evaluator.h` (`CoolPropConfig` struct), `evaluator.cpp` (thread-local cache, `getOutputValue()` helper), `solver.cpp` (config loading), `runner.cpp`, `main.cpp`, `server.cpp`.

#### 3.1.2 Superancillary Equations — **Configurable**

The `--no-superancillary` flag and `coolpropEnableSuperancillaries` config key control whether CoolProp loads superancillary data. This is now integrated into the unified `CoolPropConfig` struct and exposed in the GUI config editor. Direct use of superancillaries as a fast BisectionND backend remains a future optimization (see §4).

#### 3.1.3 Analytical Derivatives with Consistency Check — **Done**

Implemented `AbstractState::first_partial_deriv()` for exact gradients, with a **forward-FD consistency check** to handle inaccurate derivatives near phase boundaries:

1. Compute analytical derivatives via `first_partial_deriv()` (fast, 0 extra state evals)
2. Compute forward-FD check (2 extra state evaluations)
3. If analytical and forward-FD agree within 5%: use analytical derivatives
4. If they disagree: use forward-FD values (near phase boundary)
5. If `first_partial_deriv()` throws: fall through to central FD

This consistency check was essential for robustness: **CoolProp's `first_partial_deriv()` returns wildly incorrect derivatives near phase boundaries for pseudo-pure/mixture fluids** (e.g., Air ≈ N2/O2/Ar). At temperatures near Air's saturation (~79K), temperature derivatives can be off by 68× and density derivatives by 1.6×. During Newton iteration, intermediate states can drift near these regions, producing corrupted Jacobians. The consistency check detects this and falls back to FD automatically.

Configurable via `coolpropEnableAnalyticalDerivatives` (default: `true`).

**Files changed**: `evaluator.h`, `evaluator.cpp` (derivative section rewritten), `solver.cpp` (config key), GUI `ConfigEditor.tsx`.

#### 3.1.4 Residual-Only Evaluation Mode — **Done**

Added `bool computeJacobian` parameter to `BlockEvaluator::evaluate()`. When `false`:
- Derivatives are set to zero (no FD perturbations, no `first_partial_deriv()` calls)
- Jacobian matrix is not allocated
- `ExpressionEvaluator::residualOnly_` flag skips gradient propagation

Used during line search backtracking in all gradient-based solvers (Newton, TrustRegion, LM) and passed through the parallel solver lambda.

**Files changed**: `evaluator.h` (`residualOnly_`, `setResidualOnly()`, `isResidualOnly()`), `evaluator.cpp`, `solver.cpp` (lambdas pass `computeJacobian`).

### 3.2 Medium Impact — Solver Improvements

#### 3.2.1 Symbolic Block Reduction + Re-Decomposition ✅ Done

**Status**: Fully implemented in `src/solver_symbolic.cpp` + `include/coolsolve/symbolic_reduction.h` + `src/structural_analysis.cpp`.  
Controlled by `enableSymbolicReduction` (default: `false`). Re-decomposition is automatic when symbolic reduction is active.

**What was done**: Before solving, the solver attempts to symbolically reduce block size via three techniques applied iteratively until a fixed point:

1. **CoolProp call inversion**: If the model has `h = enthalpy(water, T=T, P=P)` and `h` is the matched output but one named-arg input is an unknown block variable while the other is external, reformulate as `T = temperature(water, H=h, P=P)`. CoolProp supports all standard input pairs (`HmassP_INPUTS`, `PSmass_INPUTS`, etc.).

2. **Explicit extraction**: Equations where all RHS variables are external or already-reduced are directly evaluated, removing their output variable from the block.

3. **Equation substitution**: Variables that appear in zero other block equations after extraction are removed along with their defining equation.

After reduction, **post-reduction re-decomposition** (Tarjan SCC on the reduced dependency graph) splits the remaining monolithic block into independent sub-blocks that are solved sequentially. This is particularly effective:

- `condenser_3zones`: 62-var block → reduced to 56 → re-decomposed into 13 sub-blocks (largest: 30, 15)
- `heat_pump_MSTh_SB_R10`: 39-var block → reduced to 34 → 19 sub-blocks (largest: 8, 7)
- `orc_co2`: 28-var block → reduced to 25 → 13 sub-blocks (largest: 8, 5)
- `air_screw_compressor`: 13-var block → reduced to 6 → fully decomposed into 6 scalar sub-blocks

**Robustness results** (from `examples/solver_robustness_report.md`):
- With initials, default pipeline: 35/38 (92.1%) — same success rate, with faster solve on reduced models
- Without initials, default pipeline + SymbolicReduction: 27/37 (73.0%) — same success rate as baseline
- Feature has zero overhead when disabled (verified by unit tests)

#### 3.2.2 Non-Monotone Line Search ✅ Done

**Status**: Fully implemented in all three gradient-based solvers (Newton, TrustRegion, LM).
Controlled by `lsNonMonotoneMemory` (default: `10`; set to `1` for classic monotone behavior).

**What was done**: Replaced the monotone Armijo condition with the Grippo-Lampariello-Lucidi (1986) non-monotone variant. Instead of requiring φ(x_{k+1}) < φ(x_k), the solver compares against max(φ(x_{k-M+1}), ..., φ(x_k)) where M is the memory parameter. This allows the solver to temporarily accept larger merit values to escape narrow curved valleys and saddle points, while maintaining the same global convergence guarantees.

**Bounded non-monotone reference**: the raw `max(history)` can be arbitrarily larger than the current merit when a single early bad iterate inflates the window. This makes the Armijo condition trivially satisfied, allowing steps that destroy recent progress. To prevent this, all solvers use `boundedRef(phi, R=10)` which caps the reference at `min(max(history), R × currentPhi)`. This ensures the solver can accept modest temporary increases (up to 10×) while preventing catastrophic acceptance of bad steps.

**Implementation details**:
- `NonMonotoneHistory` helper struct in `solver_common.h`: circular buffer tracking last M merit values, with `boundedRef()` for safe non-monotone reference
- **Newton**: `lineSearch()` takes a `refPhi` parameter (bounded max of history); Armijo condition uses `refPhi` instead of current φ. The directional derivative `∇φ·d` is still computed from the current iterate.
- **TrustRegion**: Step acceptance uses `phi_new < refPhi` instead of `actual > 0`. Gain ratio ρ for radius management still uses current φ.
- **LM**: Step acceptance uses `phi_new < refPhi` instead of `actual > 0`. Gain ratio ρ for lambda management still uses current φ.
- When M=1, all three solvers degenerate to exact monotone behavior (no change from previous implementation).
- Verbose output shows `refPhi` and memory parameter when non-monotone is active (M > 1).

**Robustness results** (from `examples/solver_robustness_report.md`):
- Newton (with initials): 31/38 (81.6%) — unchanged from baseline
- LM (with initials): 30/38 (+1 vs baseline, scroll_compressor now converges)
- Homotopy (with initials): 32/38 (+1 vs baseline)
- Homotopy (NO initials): 25/37 (+4 vs baseline, 67.6%)
- Newton (NO initials): 20/37 (-1 vs baseline, scroll_compressor marginal case)
- Default pipeline unchanged: 35/38 (92.1%) with initials, 27/37 (73.0%) without
- Net: +5 improvements, -1 marginal regression across all configurations

**Files changed**: `solver_common.h` (NonMonotoneHistory + boundedRef), `solver.h` (lsNonMonotoneMemory option, lineSearch signature), `solver_newton.cpp`, `solver_trust_region.cpp`, `solver_lm.cpp`, `solver.cpp` (config loader), `coolsolve.conf`, `ConfigEditor.tsx`.

**Tests**: 20 tests in `test_nonmonotone.cpp` covering NonMonotoneHistory helper (8 unit tests incl. boundedRef), Newton/TR/LM with non-monotone and monotone fallback, config loading, and trace output.

**References**:
- Grippo, Lampariello, Lucidi (1986): original non-monotone line search
- Zhang, Hager (2004): improved variant requiring *average* decrease rather than *max* reference
- NNES (Rod Bain): Fortran implementation at netlib.org/opt/nnes

#### 3.2.3 Trust Region Improvements

**Current state**: The TR solver underperforms Newton in practice (49% without initials vs 57% for Newton). Known issues:
- Overly aggressive delta shrinking (oscillation between too-small and reset)
- No smooth Cauchy-Newton transition (hard boundary)
- Always starts with `trInitialRadius=10` regardless of prior solver attempts

**What to do**: 
- Implement quadratic interpolation between Cauchy and Newton steps
- Warm-start delta from previous solver in the pipeline
- Fix the shrink/reset oscillation with a more gradual radius adaptation
- Target: TR should match Newton's success rate with initials

#### 3.2.4 Levenberg-Marquardt Improvements

**Current state**: LM solves 76% with initials (vs 82% for Newton). It fails on cases Newton handles.

**What to do**:
- Use Marquardt diagonal scaling more aggressively
- Add geodesic acceleration (Nielsen improvement) — small code change (~30 lines) that improves convergence speed significantly
- Implement delayed gratification: accept temporarily larger residuals if the step direction improves conditioning

#### 3.2.5 Broyden Quasi-Newton Updates (Jacobian Reuse)

**Background**: Broyden's method is a quasi-Newton method for solving F(x) = 0 that avoids recomputing the full Jacobian at every iteration. Instead, it starts with an initial Jacobian J₀ (e.g., from the first Newton iteration) and maintains it via rank-1 updates:

```
B_{k+1} = B_k + (ΔF - B_k Δx) Δx^T / (Δx^T Δx)     [Broyden "good" / Type I]
```

There are two main variants:
- **Broyden Type I ("good")**: Updates the Jacobian approximation B. More robust; solves B·dx = -F each iteration (cheap for moderate block sizes using pre-factored LU with rank-1 update via Sherman-Morrison).
- **Broyden Type II ("bad")**: Updates the *inverse* Jacobian approximation H directly, avoiding the linear solve entirely. Less stable but faster per iteration.

The broader quasi-Newton family includes BFGS (for optimization) and Anderson mixing/acceleration (for fixed-point iterations). Anderson acceleration is already listed as "not recommended" in this roadmap because it doesn't leverage the Jacobian structure. Broyden is fundamentally different: it *starts* from a true Jacobian and makes small corrections, preserving first-order accuracy near the solution.

**Relevance to CoolSolve**: In CoolSolve, each Jacobian column requires O(1) CoolProp evaluations (via forward-FD or analytical derivatives + consistency check). For a block of size n, the full Jacobian costs O(n) CoolProp calls. Broyden updates cost O(n²) arithmetic operations but **zero CoolProp calls** — only one residual evaluation per iteration is needed.

**Cost savings example**: For `condenser_3zones` (30-var sub-block after re-decomposition), a full Newton iteration costs ~30 CoolProp evaluations for the Jacobian + 1 for the residual = 31 calls. A Broyden iteration costs 1 call. If Newton takes 8 iterations and Broyden takes 12 (typical 1.5× factor), the total is Newton: 8 × 31 = 248 calls vs. Broyden: 1 × 31 (initial) + 11 × 1 = 42 calls — a **6× reduction** in CoolProp calls.

**Practical approach — Hybrid Newton-Broyden**:
1. Compute true Jacobian on iteration 0 (and optionally every K iterations)
2. Use Broyden rank-1 updates for intermediate iterations
3. Recompute true Jacobian when: (a) residual stalls for 3+ iterations, (b) Broyden step is rejected by line search, or (c) `B` becomes ill-conditioned
4. Implement as an enhancement to `NewtonSolver` rather than a separate solver class

**What to do**:
- Add a `broydenRecomputeInterval` option (0 = always full Jacobian = current behavior, default; K > 0 = recompute every K iterations)
- Store `B` (or `H` for Type II) between iterations; apply Sherman-Morrison for efficient rank-1 updates
- Add a condition-number check to trigger Jacobian refresh when approximation degrades
- Estimated effort: 1–2 days

**Expected impact**: Primarily **efficiency** (fewer CoolProp calls, faster convergence time) rather than robustness. Won't help models that fail due to singular Jacobians or phase boundaries, but will significantly speed up convergence on already-solvable large blocks. Best combined with non-monotone line search.

**References**:
- Broyden (1965): "A class of methods for solving nonlinear simultaneous equations"
- Dennis & Moré (1977): "Quasi-Newton methods, motivation and theory" (comprehensive survey)
- Kelley (2003): *Solving Nonlinear Equations with Newton's Method* (SIAM), Chapter 4 — practical Broyden implementation guidance

### 3.3 Lower Impact — Architecture & Testing (Not Yet Implemented)

#### 3.3.1 Newton1D Extraction

**Current state**: Newton1D is ~490 lines inlined in `Solver::solveBlock()` in solver.cpp. It's a self-contained algorithm (multi-probe + bisection + Newton for size-1 blocks).

**What to do**: Extract into `solver_newton1d.cpp` as a standalone class. This is purely a code quality improvement — it doesn't change behavior.

#### 3.3.2 Dedicated Regression Test Baseline

**Current state**: `test_solver_robustness.cpp` runs all models with various configs (tagged `[.][solver-robustness]`, slow). `test_solver_pipeline.cpp` and `test_new_solvers.cpp` cover specific solver behaviors. But there's no fast regression test that codifies "these models must solve with the default pipeline."

**What to do**: Add a lightweight `[solver-regression]` test that asserts the default pipeline solves all expected models. Update the pass-list when a solver improvement expands coverage. This gives instant CI feedback on regressions.

#### 3.3.3 Code Complexity Monitor

**What it was**: A script to track total LOC over time.

**Assessment**: Low value. The codebase is manageable (~4000 lines of solver code) and well-structured after the file split. LOC tracking doesn't provide actionable information. **Skip this.**

---

## 4. Prioritized Future Improvements

Ranked by impact on the two main objectives: **solver robustness** and **computational efficiency**.

### Tier 1 — High Impact, Bounded Effort (All Done)

| # | Improvement | Primary benefit | Status |
|---|-------------|----------------|--------|
| 1 | ~~**AbstractState caching** (§3.1.1)~~ | ~~Efficiency: 2–5× faster property evals~~ | **Done** |
| 2 | ~~**Analytical derivatives** (§3.1.3)~~ | ~~Efficiency + robustness: 5× fewer CoolProp calls, exact gradients~~ | **Done** |
| 3 | ~~**Residual-only mode** (§3.1.4)~~ | ~~Efficiency: halve CoolProp calls in line search~~ | **Done** |
| 4 | ~~**Symbolic block reduction + re-decomposition** (§3.2.1)~~ | ~~Robustness: reduce block sizes by reformulation + sub-block splitting~~ | **Done** |

All Tier 1 items are complete. The biggest remaining gap is **without-initials robustness** (73% → target ≥92%).

### Tier 2 — Medium Impact (Current Focus)

| # | Improvement | Primary benefit | Estimated effort | Priority |
|---|-------------|----------------|-----------------|----------|
| 5 | ~~**Non-monotone line search** (§3.2.2)~~ | ~~Robustness: escape narrow valleys~~ | ~~0.5–1 day~~ | **Done** |
| 6 | ⭐ **Broyden quasi-Newton updates** (§3.2.5) | Efficiency: 3–6× fewer CoolProp calls on large blocks | 1–2 days | **Next** |
| 7 | **Trust Region fixes** (§3.2.3) | Robustness: raise TR from 49% to ~80%+ without initials | 2–3 days | |
| 8 | **LM improvements** (§3.2.4) | Robustness: raise LM from 76% to ~87%+ with initials | 1–2 days | |
| 9 | **Superancillary fast evaluation** for BisectionND (§3.1.2) | Efficiency: orders of magnitude faster bisection | 2–3 days | |

### Tier 3 — Nice to Have

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|----------------|-----------------|
| 10 | **Newton1D extraction** (§3.3.1) | Code quality | 0.5 day |
| 11 | **Regression test baseline** (§3.3.2) | CI quality | 0.5 day |
| 12 | **Pseudo-arclength continuation** | Robustness: handle turning points in homotopy | 2–3 days |
| 13 | **KINSOL (SUNDIALS) integration** | Robustness: for very large blocks (>30 vars) | 1–2 weeks |

### Not Recommended

| Improvement | Why skip |
|-------------|---------|
| **Anderson acceleration** | Doesn't use Jacobian → slower than Newton when Newton works, and the pipeline already has 6 fallback strategies. The "no Jacobian needed" advantage is addressed better by BisectionND. Broyden (§3.2.5) is the better quasi-Newton approach since it starts from a true Jacobian. |
| **Code complexity monitor** | The codebase is well-structured after the file split. LOC tracking is noise, not signal. |
| **Nonlinear preconditioning (residual weighting)** | Variable scaling (already implemented) addresses the same problem. Residual weighting adds complexity for marginal benefit given that analytical derivatives (item 2) will also improve conditioning. |
| **Better initial guess via simplified model** | Requires model-specific knowledge and a meta-solving infrastructure. The pipeline + homotopy already handle poor initial guesses well. Document manual workflows instead. |

---

## 5. Key Metrics to Track

| Metric | Previous | Current | Target |
|--------|----------|---------|--------|
| Default pipeline, with initials | 17/18 (94%) | 35/38 (92.1%) | 38/38 (fix MODULE parse + remaining) |
| Default pipeline, without initials | 14/17 (82%) | 27/37 (73.0%) | ≥34/37 (92%) |
| Unit tests | 117 (805 asserts) | 149 (1111 asserts) | Expanding |
| CoolProp evals per property | 5 (1 + 4 FD) | 1–3 (AbstractState + consistency check) | 1 (analytical only, for pure fluids) |
| scroll_compressor solve time | 0.42 s | 0.11 s | ≤0.10 s |
| Cold-start warmup for R22 | 32 s | ~same | <1 s (via superancillary fast eval) |

> **Note**: The model count increased from 18 to 39 (new test examples added), so the percentages are not directly comparable. The absolute numbers show improvement.

---

## 6. Implementation Notes

### CoolProp Input Pair Reference

For CoolProp call inversion (§3.2.1), these are the supported `AbstractState::update()` input pairs:

| Input pair | Inputs | Use case |
|-----------|--------|----------|
| `PT_INPUTS` | Pressure, Temperature | Standard forward evaluation |
| `HmassP_INPUTS` | Enthalpy, Pressure | When h is known, solve for T |
| `PSmass_INPUTS` | Pressure, Entropy | Isentropic processes |
| `HmassSmass_INPUTS` | Enthalpy, Entropy | When both h and s are known |
| `DmassP_INPUTS` | Density, Pressure | When ρ is known |
| `QT_INPUTS` | Quality, Temperature | Saturation properties |

By detecting which variables are known vs unknown in a block, the evaluator can choose the best input pair to minimize the number of unknowns that need iterative solving.

### AbstractState Caching Strategy (Implemented)

```
Thread-local cache: thread_local AbstractStateCache g_abstractStateCache
  - Key: (backend, canonical_fluid_name), e.g. ("HEOS", "Water")
  - Create on first use: AbstractState::factory(backend, fluid)
  - Reuse for all subsequent calls on the same thread
  - Thread safety: thread_local storage gives each thread its own cache
  - Automatic fallback: if AbstractState throws, falls back to PropsSI
```

The thread-local approach (rather than mutex-protected shared cache) was chosen because `AbstractState::update()` is not thread-safe and cloning `AbstractState` is expensive. Each solver thread in parallel mode gets its own cache automatically.

### Analytical Derivative Consistency Check

CoolProp's `first_partial_deriv()` can return incorrect derivatives near phase boundaries for pseudo-pure fluids (mixtures treated as pure, like Air = N₂/O₂/Ar). The issue was identified during debugging of the `air_screw_compressor` model regression:

| State | Property | Analytical | FD | Relative Error |
|-------|----------|-----------|-----|---------------|
| Air, P=1bar, T≈81K (near sat) | dT/dS | 0.0747 | 0.0011 | **6693%** |
| Air, P=1bar, T≈81K (near sat) | dT/dP | 2.29e-4 | 8.67e-5 | **164%** |
| Air, P=1bar, T≈81K (near sat) | dD/dS | -4.74e-3 | -1.81e-3 | **162%** |
| Air, P=1bar, T=293K (normal) | dT/dP | 8.32e-4 | 8.32e-4 | <1e-9 |

The consistency check (5% tolerance on forward-FD comparison) adds 2 extra `state->update()` calls per property but catches these errors before they corrupt the Jacobian. For pure fluids away from phase boundaries, the check passes and analytical derivatives are used (accurate to better than 1e-9 relative error).
