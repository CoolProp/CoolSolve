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

110 unit tests, all passing.

### 1.4 Robustness Results Summary

| Configuration | With initials | Without initials |
|---|---:|---:|
| Default pipeline (6 solvers) | 17/18 (94%) | 14/17 (82%) |
| Newton only | 17/18 (94%) | 14/17 (82%) |
| TrustRegion only | 17/18 (94%) | 9/17 (53%) |
| LevenbergMarquardt only | 15/18 (83%) | 10/17 (59%) |
| Partitioned only | 10/18 (56%) | 10/17 (59%) |

The hardest models (scroll_compressor, humidair2, orc_co2) share: large algebraic blocks (28–62 vars), CoolProp evaluations near phase boundaries, and near-singular Jacobians.

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
| CoolProp call profiling | **Done** | Per-call timing for PropsSI and HAPropsSI; warmup time tracked separately |
| `--no-superancillary` flag | **Done** | Avoids loading the cached superancillary data |
| Detailed profiling | **Done** | First CoolProp call triggers a warmup for fluid init |
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

## 3. What Has NOT Been Implemented

### 3.1 High Impact — CoolProp Integration

These are the highest-leverage changes. CoolProp calls dominate solve time, and each property evaluation currently requires **5 PropsSI calls** (1 for the value + 4 for finite-difference Jacobian).

#### 3.1.1 Low-Level CoolProp API (`AbstractState`)

**Current state**: All CoolProp calls go through the high-level `PropsSI()` function, which performs string parsing, fluid lookup, and EOS initialization on every call. The `AbstractState` infrastructure is declared in `evaluator.h` (`cacheEnabled`, `cacheTolerance`) and the header is included, but the actual evaluation code in `evaluator.cpp` still calls `PropsSI()`.

**What to do**: Replace `PropsSI` calls with:
1. Create and cache `AbstractState` objects per fluid (one-time cost)
2. Call `state->update(input_pair, val1, val2)` + `state->hmass()` etc. for property evaluation
3. This eliminates string parsing and fluid lookup overhead on every call

**Expected impact**: 2–5× speedup on property evaluations (the dominant cost in most models). The `update()` path skips string matching, unit conversion, and backend selection that `PropsSI` performs every time.

#### 3.1.2 Superancillary Equations for Fast Saturation

**Current state**: `--no-superancillary` exists and allows significant speedup with the high level interface by using simpler saturation routines. However this disables superancillaries globally.

**What to do**: Instead of disabling superancillaries, *use them directly* as fast function evaluators during bisection-type solvers. Superancillary equations are polynomial fits for saturation properties; they can be evaluated in microseconds (vs milliseconds for full EOS). For BisectionND, which needs many evaluations but doesn't need high accuracy in intermediate steps, superancillaries could be the evaluation backend.

**Practical approach**: Use full-accuracy `AbstractState` for Newton/TR/LM solvers (they need few but accurate evaluations), and fast superancillary evaluations for BisectionND (it needs many cheap evaluations to narrow down the search space). Switch to full accuracy for the final polish.

#### 3.1.3 Analytical Derivatives

**Current state**: Jacobians for CoolProp functions use central finite differences (4 extra PropsSI calls per property). `AbstractState::first_partial_deriv()` exists in CoolProp but is not used.

**What to do**: After switching to `AbstractState` (3.1.1), use `first_partial_deriv()` for exact gradients. This eliminates the 4 finite-difference calls entirely.

**Expected impact**: 5× reduction in CoolProp calls (1 evaluation instead of 5 per property). Exact gradients also improve solver robustness: finite-difference gradients can be noisy near phase boundaries, causing line search failures and slow convergence.

**Caveat**: `first_partial_deriv()` works well for single-phase states but may have issues at two-phase boundaries. A fallback to finite differences for problematic states would be prudent.

#### 3.1.4 Residual-Only Evaluation Mode

**Current state**: `BlockEvaluator::evaluate()` always computes both residuals and the Jacobian. During line search, only the residual norm is needed.

**What to do**: Add a `bool computeJacobian` flag to `evaluate()`. Skip Jacobian computation during line search steps.

**Expected impact**: During line search backtracking (which can take 5–15 steps per Newton iteration), each step only needs 1 CoolProp call instead of 5. For models where line search is active, this could halve total CoolProp calls.

### 3.2 Medium Impact — Solver Improvements

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

#### 3.2.3 Trust Region Improvements

**Current state**: The TR solver underperforms Newton in practice (53% without initials vs 82% for Newton). Known issues:
- Overly aggressive delta shrinking (oscillation between too-small and reset)
- No smooth Cauchy-Newton transition (hard boundary)
- Always starts with `trInitialRadius=10` regardless of prior solver attempts

**What to do**: 
- Implement quadratic interpolation between Cauchy and Newton steps
- Warm-start delta from previous solver in the pipeline
- Fix the shrink/reset oscillation with a more gradual radius adaptation
- Target: TR should match Newton's success rate with initials

#### 3.2.4 Levenberg-Marquardt Improvements

**Current state**: LM solves 83% with initials (vs 94% for Newton). It fails on cases Newton handles.

**What to do**:
- Use Marquardt diagonal scaling more aggressively
- Add geodesic acceleration (Nielsen improvement) — small code change (~30 lines) that improves convergence speed significantly
- Implement delayed gratification: accept temporarily larger residuals if the step direction improves conditioning

### 3.3 Lower Impact — Architecture & Testing

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
| 1 | **AbstractState caching** (§3.1.1) | Efficiency: 2–5× faster property evals | 1–2 days |
| 2 | **Analytical derivatives** (§3.1.3) | Efficiency + robustness: 5× fewer CoolProp calls, exact gradients | 1–2 days (after #1) |
| 3 | **Residual-only mode** (§3.1.4) | Efficiency: halve CoolProp calls in line search | 0.5 day |
| 4 | **CoolProp call inversion** (§3.2.1) | Robustness: reduce block sizes by reformulating property calls | 2–3 days |

Items 1–3 are independent and attack the core bottleneck: CoolProp is called 5× more than necessary. Combined, they could yield a **5–10× speedup** on CoolProp-heavy models. Item 4 directly reduces block sizes, which is the single most effective way to improve solver convergence.

### Tier 2 — Medium Impact

| # | Improvement | Primary benefit | Estimated effort |
|---|-------------|----------------|-----------------|
| 5 | **Lazy fluid initialization** (§3.2.2) | Efficiency: eliminate 32 s cold-start warmup | 1 day |
| 6 | **Trust Region fixes** (§3.2.3) | Robustness: raise TR from 53% to ~80%+ without initials | 2–3 days |
| 7 | **LM improvements** (§3.2.4) | Robustness: raise LM from 83% to ~94%+ with initials | 1–2 days |
| 8 | **Superancillary fast evaluation** for BisectionND (§3.1.2) | Efficiency: orders of magnitude faster bisection | 2–3 days |
| 9 | **Symbolic equation substitution** (§3.2.1) | Robustness: further block size reduction | 3–5 days |

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
| **Anderson acceleration** | Doesn't use Jacobian → slower than Newton when Newton works, and the pipeline already has 6 fallback strategies. The "no Jacobian needed" advantage is addressed better by BisectionND. |
| **Code complexity monitor** | The codebase is well-structured after the file split. LOC tracking is noise, not signal. |
| **Nonlinear preconditioning (residual weighting)** | Variable scaling (already implemented) addresses the same problem. Residual weighting adds complexity for marginal benefit given that analytical derivatives (item 2) will also improve conditioning. |
| **Better initial guess via simplified model** | Requires model-specific knowledge and a meta-solving infrastructure. The pipeline + homotopy already handle poor initial guesses well. Document manual workflows instead. |

---

## 5. Key Metrics to Track

| Metric | Current | Target |
|--------|---------|--------|
| Default pipeline, with initials | 17/18 (94%) | 18/18 (fix MODULE parse) |
| Default pipeline, without initials | 14/17 (82%) | ≥16/17 (94%) |
| Avg CoolProp calls per model solve | ~500–5000 | ~100–1000 (via items 1–3) |
| orc_r245fa solve time | 0.30 s | ≤0.10 s (via items 1–3) |
| scroll_compressor solve time | 0.42 s | ≤0.20 s |
| Cold-start warmup for R22 | 32 s | <1 s (via item 5) |

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

### AbstractState Caching Strategy

```
Per-fluid cache: map<string, shared_ptr<AbstractState>>
  - Key: canonical fluid name (e.g., "Water", "R245fa")
  - Create on first use: AbstractState::factory("HEOS", fluid)
  - Reuse for all subsequent calls
  - Thread-safety: one AbstractState per thread (or mutex-protected)
```

When switching to parallel pipeline mode, each solver thread needs its own `AbstractState` copy (or a thread-local cache) since `AbstractState::update()` is not thread-safe.
