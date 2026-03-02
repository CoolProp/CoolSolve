# CoolSolve — Solver Roadmap & Status

*Last updated: March 2026*

This document consolidates the performance, robustness, and architecture plans for CoolSolve's solver subsystem. It replaces the earlier `performance_plan.md`, `solver_fix_strategy.md`, and `solver_robustness_diagnosis.md`.

---

## 1. Architecture Overview

### 1.1 Solver Pipeline

CoolSolve uses a configurable multi-solver pipeline. The default (sequential mode) tries each strategy in order until one succeeds:

**Newton → TrustRegion → LevenbergMarquardt → BisectionND → Homotopy → Partitioned**

A parallel mode (first-to-solve-wins via `std::async` threads) is also available. BisectionND is automatically skipped for blocks larger than `bisectionNDMaxBlockSize` (default 8).

### 1.2 File Structure

| File | Lines | Content |
|------|------:|---------|
| `solver.h` | 751 | All solver class declarations, `SolverOptions`, enums |
| `solver_common.h` | 49 | Shared `computeScalingFactors()` utility |
| `solver.cpp` | ~1940 | Orchestrator, Newton1D, tearing, partitioned, pipeline |
| `solver_newton.cpp` | 218 | Newton + line search |
| `solver_trust_region.cpp` | ~240 | Trust Region Dogleg |
| `solver_lm.cpp` | ~178 | Levenberg-Marquardt |
| `solver_bisection_nd.cpp` | 451 | N-dimensional bisection (sign-pattern simplex) |
| `solver_homotopy.cpp` | 210 | Homotopy continuation (predictor-corrector) |

### 1.3 Test Results (current baseline)

| Model | Eqs | Blocks | Solve time | Status |
|-------|----:|-------:|-----------:|--------|
| condenser_3zones | 99 | 50 | 0.03 s | OK |
| exchangers1 | 20 | 20 | 0.01 s | OK |
| exchangers2 | 29 | 26 | 0.00 s | OK |
| exchangers3 | 35 | 33 | 0.03 s | OK |
| humidair1 | 25 | 25 | 0.00 s | OK |
| humidair2 | 37 | 33 | 0.02 s | OK |
| orc_co2 | 139 | 112 | 0.05 s | OK |
| orc_complex | — | — | — | FAIL (parse: unsupported `MODULE`) |
| orc_extraction | 133 | 113 | 0.08 s | OK |
| orc_r245fa | 173 | 151 | 0.30 s | OK |
| orc_simple | 172 | 150 | 0.10 s | OK |
| pressuredrop | 26 | 26 | 0.01 s | OK |
| rankine1 | 30 | 30 | 0.04 s | OK |
| rankine2 | 45 | 42 | 0.09 s | OK |
| refrigeration1 | 39 | 39 | 0.05 s | OK |
| refrigeration2 | 39 | 39 | 0.04 s | OK |
| refrigeration3 | 38 | 38 | 0.01 s | OK |
| scroll_compressor | 99 | 66 | 0.42 s | OK |

117 unit tests (805 assertions), all passing.

### 1.4 Robustness Results Summary

| Configuration | With initials | Without initials |
|---|---:|---:|
| Default pipeline (6 solvers) | 35/38 (92%) | 27/37 (73%) |
| Newton only | 31/38 (82%) | 21/37 (57%) |
| TrustRegion only | 33/38 (87%) | 18/37 (49%) |
| LevenbergMarquardt only | 29/38 (76%) | 17/37 (46%) |
| Partitioned only | 17/38 (45%) | 14/37 (38%) |
| Default + Tearing | 35/38 (92%) | 28/37 (76%) |

39 example models tested (38 parseable). The hardest models (scroll_compressor, cooling_tower2, orc_co2) share: large algebraic blocks (28–62 vars), CoolProp evaluations near phase boundaries, and near-singular Jacobians.

---

## 2. What Has Been Implemented

### 2.1 Solver Algorithms

| Feature | Status | Notes |
|---------|--------|-------|
| Newton + adaptive line search | **Done** | Default workhorse; strongest single solver |
| Trust Region Dogleg | **Done** | Full implementation with configurable radius, eta, shrink/grow factors |
| Levenberg-Marquardt | **Done** | Damped least-squares with diagonal scaling |
| BisectionND (simplex bisection) | **Done** | Sign-pattern simplex; `bisectionNDMaxBlockSize` (default 8), `bisectionNDIterFactor` |
| Homotopy continuation | **Done** | Predictor-corrector with adaptive step, Newton+LM corrector fallback, final polish |
| Partitioned solver | **Done** | Per-variable diagonal updates; gracefully returns `MaxIterations` for small blocks |
| Newton1D | **Done** | Specialized for size-1 blocks: multi-probe + bisection + Newton hybrid |
| Tearing | **Done** | Greedy FVS + acyclic sequential solve + Newton on tear residuals |
| Explicit solve (size-1 blocks) | **Done** | `tryExplicitSolve()` detects `x = expr(known)` pattern, evaluates RHS directly |
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

### 3.2 Medium Impact — Solver Improvements (Not Yet Implemented)

#### 3.2.1 Symbolic Block Reduction

**Current state**: Blocks are determined by structural analysis (finding strongly connected components in the dependency graph). Once formed, a block is solved as-is.

**What to do**: Before solving, attempt to symbolically reduce block size:

1. **Equation substitution**: If a block equation is `x = f(known_vars, y)` and `x` only appears in one other equation, substitute `f(known_vars, y)` for `x` in that equation and remove `x` from the block. This directly reduces block size.

2. **CoolProp call inversion**: If the model has `h = enthalpy(water, T=T, P=P)` and `h` is known but `T` is unknown, reformulate as `T = temperature(water, H=h, P=P)`. CoolProp supports all standard input pairs (`HmassP_INPUTS`, `PSmass_INPUTS`, etc.), so many thermodynamic equations can be inverted to solve directly for the unknown. This can turn a 3-variable block (h, T, P) into three size-1 blocks.

3. **Algebraic simplification**: For purely algebraic equations (no CoolProp), recognize patterns like `a*x + b = 0` and solve analytically.

**Expected impact**: Could reduce the hardest blocks (28–62 variables) significantly. Fewer variables means better conditioning, faster convergence, and potentially eliminating the need for heavyweight solvers on those blocks.

**Caveat**: This requires parsing and manipulating the equation IR, which adds complexity. Start with CoolProp inversion (high value, bounded effort), then algebraic substitution.

#### 3.2.2 Lazy Fluid Initialization

**Current state**: The first `CoolProp::PropsSI` call for any fluid triggers expensive one-time initialization (~32 s for R22 in debug mode and a few hundreds of ms in release mode). This happens during variable inference, before solving starts.

**What to do**: Defer the first PropsSI call by using heuristic-based default values in `initializeVariables()`. Only call CoolProp when the solver actually needs property evaluations. Alternatively, start fluid initialization in a background thread during parsing so it completes by the time solving begins.

**Expected impact**: Eliminates the warmup for cold starts. With `--no-superancillary`, the warmup is significantly decreased, so this becomes less critical — but it still matters for interactive use (GUI) where the user expects fast response.

#### 3.2.3 Non-Monotone Line Search

**Background**: The SNEES solver (based on Rod Bain's NNES, available at netlib.org/opt/nnes) uses a *non-monotone* line search strategy, originally proposed by Grippo, Lampariello, and Lucidi (1986). This was chosen for its superior convergence in "difficult solution landscapes" and reportedly achieves performance comparable to EES's internal solvers on nonlinear equation systems. The technique is coupled with a sparse equation solver in SNEES for additional speed on larger systems.

**Current state**: CoolSolve's Newton solver uses a standard **monotone Armijo** backtracking line search: the merit function `φ(x) = ½‖F(x)‖²` must decrease at every step. This can cause the solver to stall in narrow curved valleys or near saddle points, where every descent direction leads to an increase before reaching the next valley floor.

**What to do**: Replace the monotone Armijo condition with a non-monotone variant. Instead of requiring:

```
φ(x + λd) ≤ φ(x) - c₁ λ |∇φ·d|
```

require:

```
φ(x + λd) ≤ max(φ(x_{k-M}), ..., φ(x_k)) - c₁ λ |∇φ·d|
```

where M is a memory parameter (typically 5–15 recent function values). The implementation change is small: maintain a circular buffer of the last M merit function values and compare against the maximum rather than the current value.

This can be applied to all three gradient-based solvers (Newton, Trust Region, Levenberg-Marquardt). For TR, the non-monotone condition governs the trust region acceptance ratio; for LM, it governs the step acceptance criterion.

**Expected impact**: Better convergence on the hardest models (scroll_compressor, humidair2, orc_co2) which have near-singular Jacobians and highly nonlinear CoolProp evaluations near phase boundaries. Non-monotone methods:
- Navigate narrow curved valleys that trap monotone line search
- Accept unit step lengths more often, reducing backtracking overhead
- Achieve global convergence under weaker conditions than monotone methods
- Maintain the same worst-case complexity guarantees

**References**:
- Grippo, Lampariello, Lucidi (1986): original non-monotone line search
- Zhang, Hager (2004): improved variant requiring *average* decrease rather than *max* reference
- NNES (Rod Bain): Fortran implementation at netlib.org/opt/nnes
- SNEES: engineering equation solver using NNES + sparse linear algebra

#### 3.2.4 Trust Region Improvements

**Current state**: The TR solver underperforms Newton in practice (53% without initials vs 82% for Newton). Known issues:
- Overly aggressive delta shrinking (oscillation between too-small and reset)
- No smooth Cauchy-Newton transition (hard boundary)
- Always starts with `trInitialRadius=10` regardless of prior solver attempts

**What to do**: 
- Implement quadratic interpolation between Cauchy and Newton steps
- Warm-start delta from previous solver in the pipeline
- Fix the shrink/reset oscillation with a more gradual radius adaptation
- Target: TR should match Newton's success rate with initials

#### 3.2.5 Levenberg-Marquardt Improvements

**Current state**: LM solves 83% with initials (vs 94% for Newton). It fails on cases Newton handles.

**What to do**:
- Use Marquardt diagonal scaling more aggressively
- Add geodesic acceleration (Nielsen improvement) — small code change (~30 lines) that improves convergence speed significantly
- Implement delayed gratification: accept temporarily larger residuals if the step direction improves conditioning

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

### Tier 1 — High Impact, Bounded Effort

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|----------------|-----------------|
| 1 | ~~**AbstractState caching** (§3.1.1)~~ | ~~Efficiency: 2–5× faster property evals~~ | **Done** |
| 2 | ~~**Analytical derivatives** (§3.1.3)~~ | ~~Efficiency + robustness: 5× fewer CoolProp calls, exact gradients~~ | **Done** |
| 3 | ~~**Residual-only mode** (§3.1.4)~~ | ~~Efficiency: halve CoolProp calls in line search~~ | **Done** |
| 4 | **CoolProp call inversion** (§3.2.1) | Robustness: reduce block sizes by reformulating property calls | 2–3 days |

Items 1–3 are independent and attack the core bottleneck: CoolProp is called 5× more than necessary. Combined, they could yield a **5–10× speedup** on CoolProp-heavy models. Item 4 directly reduces block sizes, which is the single most effective way to improve solver convergence.

### Tier 2 — Medium Impact

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|----------------|-----------------|
| 5 | **Non-monotone line search** (§3.2.3) | Robustness: better convergence on difficult landscapes | 0.5–1 day |
| 6 | **Lazy fluid initialization** (§3.2.2) | Efficiency: eliminate 32 s cold-start warmup | 1 day |
| 7 | **Trust Region fixes** (§3.2.4) | Robustness: raise TR from 53% to ~80%+ without initials | 2–3 days |
| 8 | **LM improvements** (§3.2.5) | Robustness: raise LM from 83% to ~94%+ with initials | 1–2 days |
| 9 | **Superancillary fast evaluation** for BisectionND (§3.1.2) | Efficiency: orders of magnitude faster bisection | 2–3 days (superancillary config done, fast eval backend remaining) |
| 10 | **Symbolic equation substitution** (§3.2.1) | Robustness: further block size reduction | 3–5 days |

### Tier 3 — Nice to Have

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|----------------|-----------------|
| 11 | **Newton1D extraction** (§3.3.1) | Code quality | 0.5 day |
| 12 | **Regression test baseline** (§3.3.2) | CI quality | 0.5 day |
| 13 | **Pseudo-arclength continuation** | Robustness: handle turning points in homotopy | 2–3 days |
| 14 | **KINSOL (SUNDIALS) integration** | Robustness: for very large blocks (>30 vars) | 1–2 weeks |

### Not Recommended

| Improvement | Why skip |
|-------------|---------|
| **Anderson acceleration** | Doesn't use Jacobian → slower than Newton when Newton works, and the pipeline already has 6 fallback strategies. The "no Jacobian needed" advantage is addressed better by BisectionND. |
| **Code complexity monitor** | The codebase is well-structured after the file split. LOC tracking is noise, not signal. |
| **Nonlinear preconditioning (residual weighting)** | Variable scaling (already implemented) addresses the same problem. Residual weighting adds complexity for marginal benefit given that analytical derivatives (item 2) will also improve conditioning. |
| **Better initial guess via simplified model** | Requires model-specific knowledge and a meta-solving infrastructure. The pipeline + homotopy already handle poor initial guesses well. Document manual workflows instead. |

---

## 5. Key Metrics to Track

| Metric | Previous | Current | Target |
|--------|----------|---------|--------|
| Default pipeline, with initials | 17/18 (94%) | 35/38 (92%) | 38/38 (fix MODULE parse + remaining) |
| Default pipeline, without initials | 14/17 (82%) | 27/37 (73%) | ≥34/37 (92%) |
| CoolProp evals per property | 5 (1 + 4 FD) | 1–3 (AbstractState + consistency check) | 1 (analytical only, for pure fluids) |
| Avg CoolProp calls per model solve | ~500–5000 | ~100–1500 (via §3.1) | ~100–500 (via §3.2.1 block reduction) |
| scroll_compressor solve time | 0.42 s | 0.11 s | ≤0.10 s |
| Cold-start warmup for R22 | 32 s | ~same | <1 s (via item 6) |

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
