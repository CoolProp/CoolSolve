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
| `solver.h` | ~906 | All solver class declarations, `SolverOptions`, enums, `Newton1DSolver` |
| `solver_common.h` | ~140 | Shared utilities: `computeScalingFactors()`, `NonMonotoneHistory`, `isFatalEvaluationError()` |
| `solver.cpp` | ~2378 | Orchestrator, tearing, partitioned, pipeline, symbolic reduction integration |
| `solver_newton.cpp` | ~270 | Newton + non-monotone Armijo backtracking line search + Broyden quasi-Newton |
| `solver_newton1d.cpp` | 355 | Specialized 1D root-finding: 4-phase Newton + probing + bisection hybrid |
| `solver_trust_region.cpp` | 238 | Trust Region Dogleg |
| `solver_lm.cpp` | 177 | Levenberg-Marquardt |
| `solver_bisection_nd.cpp` | 450 | N-dimensional bisection (sign-pattern simplex) |
| `solver_homotopy.cpp` | 209 | Homotopy continuation (predictor-corrector) |
| `solver_symbolic.cpp` | 611 | Symbolic block reduction (inversion + extraction + substitution) |
| `symbolic_reduction.h` | 147 | Symbolic reduction data structures and API |
| `structural_analysis.h` | 155 | Graph decomposition API (Tarjan SCC + redecomposition) |
| `structural_analysis.cpp` | 763 | Tarjan SCC, Dulmage-Mendelsohn, block re-decomposition |

### 1.3 Test Results (current baseline)

195 unit tests (1291 assertions), all passing. Models include:

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
| Newton + non-monotone line search | **Done** | Default workhorse; Grippo et al. 1986 non-monotone Armijo (M=10); Broyden quasi-Newton rank-1 updates (§3.2.5) |
| Trust Region Dogleg | **Done** | Adaptive initial radius (Cauchy step norm), smoother rho-based radius adaptation, gradient-based rejection recovery (§3.2.3) |
| Levenberg-Marquardt | **Done** | Nielsen's λ adaptation, cumulative Marquardt diagonal scaling, geodesic acceleration (Transtrum & Sethna 2012) (§3.2.4) |
| BisectionND (simplex bisection) | **Done** | Sign-pattern simplex; `bisectionNDMaxBlockSize` (default 8), `bisectionNDIterFactor` |
| Homotopy continuation | **Done** | Predictor-corrector with adaptive step, Newton+LM corrector fallback, final polish |
| Partitioned solver | **Done** | Per-variable diagonal updates; gracefully returns `MaxIterations` for small blocks |
| Newton1D | **Done** | Specialized for size-1 blocks: multi-probe + bisection + Newton hybrid (extracted to `solver_newton1d.cpp`) |
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

#### 3.2.3 Trust Region Improvements ✅ Done

**Status**: Fully implemented in `src/solver_trust_region.cpp`.

**What was done**:
- **Adaptive initial radius** (`trAdaptiveRadius`, default `true`): Sets the initial trust radius from the Cauchy step norm on the first iteration (`delta = min(trInitialRadius, max(cauchyNorm, 1.0))`), automatically scaling to the problem geometry.
- **Smoother rho-based radius adaptation**: Replaces binary grow/no-grow with graduated rho thresholds (grow at ρ > 0.75, hold at ρ > 0.5, shrink proportionally otherwise using `max(0.1, 1−ρ) × delta`).
- **Gradient-based recovery**: After many consecutive rejections (>15), the solver resets to a gradient-direction step with a small radius, escaping degenerate trust regions.

**Files changed**: `solver_trust_region.cpp`, `solver.h` (`trAdaptiveRadius` option).

#### 3.2.4 Levenberg-Marquardt Improvements ✅ Done

**Status**: Fully implemented in `src/solver_lm.cpp`.

**What was done**:
- **Nielsen's λ adaptation** (`lmNielsenUpdate`, default `true`): Uses `λ = λ × max(1/3, 1 − (2ρ−1)³)` on acceptance and exponential increase (`λ × ν` with `ν` doubling) on consecutive rejections (Madsen et al. 2004). Provides smoother, faster-converging λ transitions.
- **Cumulative Marquardt diagonal scaling**: `D_k = max(D_{k−1}, diag(J^T J))`, preventing scale collapse when the Jacobian changes dramatically between iterations.
- **Geodesic acceleration** (`lmGeodesicAcceleration`, default `true`): Adds a second-order correction to the LM step by evaluating the directional second derivative of F along the velocity step (Transtrum & Sethna 2012). Costs 1 extra residual evaluation per iteration but can halve iteration count on curved problems. Includes acceleration/velocity ratio guard to avoid oversized corrections.
- **"Delayed gratification"** is effectively provided by the non-monotone acceptance criterion (§3.2.2).

**Files changed**: `solver_lm.cpp`, `solver.h` (`lmNielsenUpdate`, `lmGeodesicAcceleration` options).

#### 3.2.5 Broyden Quasi-Newton Updates (Jacobian Reuse) ✅ Done

**Status**: Fully implemented as a hybrid Newton-Broyden mode in `src/solver_newton.cpp`.

**What was done**: Broyden's method is a quasi-Newton approach that avoids recomputing the full Jacobian at every iteration by maintaining rank-1 updates:

```
B_{k+1} = B_k + (ΔF - B_k Δx) Δx^T / (Δx^T Δx)     [Broyden Type I]
```

- **`broydenRecomputeInterval`** option (default `0` = always full Jacobian; `K > 0` = recompute every K iterations). Recommended starting point: K = 5.
- Computes true Jacobian on iteration 0, then uses Broyden Type-I rank-1 updates for intermediate iterations.
- **Automatic fallback**: Forces full Jacobian recomputation when (a) Broyden step fails line search, (b) residual stalls for ≥3 iterations, or (c) scheduled recompute interval reached.
- Primarily an **efficiency** improvement (fewer CoolProp calls on large blocks) rather than robustness.

**Cost savings example**: For a 30-var block, full Newton costs ~31 CoolProp calls/iteration (30 for Jacobian + 1 for residual). With Broyden (K=5), on average only ~7 calls/iteration (1 full Jacobian every 5 iters + 4 residual-only iters).

**Files changed**: `solver_newton.cpp` (Broyden state, rank-1 updates, fallback logic), `solver.h` (`broydenRecomputeInterval`), `solver.cpp` (config loading).

**References**:
- Broyden (1965): "A class of methods for solving nonlinear simultaneous equations"
- Dennis & Moré (1977): "Quasi-Newton methods, motivation and theory"
- Kelley (2003): *Solving Nonlinear Equations with Newton's Method* (SIAM), Chapter 4

---

## 4. Prioritized Future Improvements

Ranked by impact on the two main objectives: **solver robustness** and **computational efficiency**.


### Tier 3 — Nice to Have

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|----------------|-----------------|
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
