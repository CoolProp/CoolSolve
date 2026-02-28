# CoolSolve Solver Robustness Diagnosis & Roadmap

## 1. Current State: Robustness Test Results

### 1.1 Summary Table (with initial values)

| Configuration | Solved/Total | Rate |
|---|---:|---:|
| Default pipeline (Nwt+TR+LM+Part) | 17/18 | 94.4% |
| Newton only | 17/18 | 94.4% |
| TrustRegion only | 17/18 | 94.4% |
| LevenbergMarquardt only | 15/18 | 83.3% |
| Partitioned only | 10/18 | 55.6% |
| Default + Tearing | 17/18 | 94.4% |

### 1.2 Summary Table (without initial values)

| Configuration | Solved/Total | Rate |
|---|---:|---:|
| Default pipeline (Nwt+TR+LM+Part) | 14/17 | 82.4% |
| Newton only | 14/17 | 82.4% |
| TrustRegion only | 9/17 | 52.9% |
| LevenbergMarquardt only | 10/17 | 58.8% |
| Partitioned only | 10/17 | 58.8% |
| Default + Tearing | 14/17 | 82.4% |

> Note: `orc_complex.eescode` fails at parse (uses unsupported `MODULE` construct) and is excluded from "without initials" totals.

### 1.3 Model Difficulty Ranking

| Model | Failure count (across 12 configs) | Root cause |
|---|---:|---|
| scroll_compressor.eescode | 7/12 | Large block (34 vars), stiff nonlinear coupling |
| humidair2.eescode | 7/12 | HumidAir evaluation errors (invalid inputs near boundaries) |
| orc_co2.eescode | 7/12 | Large block (28 vars), singular Jacobian near phase boundaries |
| orc_complex.eescode | 6/6 | Parse failure (MODULE construct not supported) |
| orc_extraction.eescode | 5/12 | 21-var block, singular Jacobian without good initials |
| exchangers3.eescode | 4/12 | Small block (3 vars) but TR and LM fail without initials |
| condenser_3zones.eescode | 3/12 | Huge block (62 vars), only Newton+pipeline succeeds without initials |

### 1.4 Error Category Breakdown

| Error Category | Count | Fraction |
|---|---:|---:|
| Max iterations | 17 | 37.0% |
| Evaluation error | 12 | 26.1% |
| Other (incl. Partitioned min block size) | 12 | 26.1% |
| Singular Jacobian | 5 | 10.9% |

### 1.5 Key Observations

1. **Newton is the strongest single solver.** With initials, Newton alone matches the full pipeline (94.4%). Without initials, it still leads (82.4%) — it solves `scroll_compressor` in 15s where even the full pipeline fails.

2. **The pipeline adds no value over Newton alone (with initials).** When initial values are provided, Newton never fails on cases that the pipeline solves. This means TR, LM, and Partitioned are not contributing to the with-initials success rate.

3. **TrustRegion is surprisingly weak without initials (52.9%).** It reaches max iterations on many models (exchangers3, orc_co2, orc_extraction, orc_r245fa, rankine2, scroll_compressor). The trust region tends to shrink too aggressively and stall.

4. **LevenbergMarquardt is also weak without initials (58.8%).** It struggles on the same blocks as TR and additionally fails on orc_simple.

5. **Partitioned is a niche solver** that only works on well-structured problems. It fails on blocks smaller than `partitionedMinBlockSize=4` (returning InvalidInput), limiting its utility. Its 55.6% with-initials rate makes it the weakest solver.

6. **Tearing helps on large blocks** (condenser_3zones: 62 vars, scroll_compressor: 34 vars) but adds significant overhead when it doesn't help. With tearing, most solve times increase 2-5×.

7. **The three hardest models** (scroll_compressor, humidair2, orc_co2) all involve:
   - Large algebraic blocks (28-62 variables)
   - CoolProp evaluations near phase boundaries
   - Singular or near-singular Jacobians

---

## 2. Solver Architecture Analysis

### 2.1 Current Architecture: Single File Problem

All solver implementations live in `src/solver.cpp` (~2850 lines):

| Section | Lines | Content |
|---|---:|---|
| Utilities, config loading, timeout | 1-390 | statusToString, categorizeError, config parser |
| NewtonSolver::solve + lineSearch | 391-710 | ~320 lines |
| TrustRegionSolver::solve + dogleg | 712-1040 | ~330 lines |
| LevenbergMarquardtSolver::solve | 1044-1280 | ~240 lines |
| Newton1D (inline in solveBlock) | 1300-1790 | ~490 lines |
| Solver orchestrator + pipeline | 1790-2240 | ~450 lines |
| Partitioned solver | 2240-2420 | ~180 lines |
| Tearing solver | 2422-2700 | ~280 lines |
| Solver::solve (main) + report | 2700-2850 | ~150 lines |

**Recommendation: Split into separate files.** Each solver is already a self-contained class with the `NonLinearSolver` interface. Splitting would:
- Make each solver independently testable and reviewable
- Reduce merge conflicts when multiple solvers change
- Make the codebase easier to navigate

Proposed structure:
```
src/solver/
  solver_common.cpp       # SolverOptions, config loading, timeout, error utils
  solver_newton.cpp       # NewtonSolver + lineSearch
  solver_trust_region.cpp # TrustRegionSolver + dogleg
  solver_lm.cpp           # LevenbergMarquardtSolver
  solver_newton1d.cpp     # Newton1D root finder (extract from solveBlock)
  solver_partitioned.cpp  # Partitioned solver
  solver_tearing.cpp      # Tearing solver
  solver_pipeline.cpp     # Pipeline orchestrator (sequential, parallel)
  solver.cpp              # Solver::solve (main orchestrator)
```

### 2.2 Code Duplication Issues

Several patterns are duplicated across all solvers:

1. **`computeScalingFactors()`** — identical code in Newton, TR, and LM (~15 lines each). Should be a free function or base class method.

2. **Convergence checking** — each solver has its own tolerance/relative-tolerance/step-tolerance checks. Should be a shared utility.

3. **Trace recording** — every solver has near-identical trace-recording boilerplate. A helper like `recordTrace(trace, iter, residualNorm, ...)` would reduce ~30 lines per solver.

4. **Jacobian factorization** — Newton, TR, and LM all use `ColPivHouseholderQR`. A shared "solve linear system with fallback" utility would be cleaner.

### 2.3 Are All Solvers Useful?

| Solver | Useful? | Assessment |
|---|---|---|
| **Newton + Line Search** | **Essential** | The workhorse. Best single solver. Must keep. |
| **Trust Region Dogleg** | **Questionable** | Theoretically more robust than Newton for badly-conditioned problems, but in practice it performs *worse* than Newton on all test cases. The TR implementation may need debugging — its shrink/reset logic is complex and may be overly conservative. |
| **Levenberg-Marquardt** | **Keep but improve** | Fills a useful niche (good when initial guess is far from solution) but its implementation could be improved. Currently fails on cases where Newton succeeds. |
| **Partitioned** | **Niche, keep as last resort** | Designed for specific block structures. The `partitionedMinBlockSize=4` threshold causes it to fail with InvalidInput on small blocks (2-3 vars) that the other solvers handle fine. Useful as a last-resort stabilizer for specific block types. |
| **Tearing** | **Valuable for large blocks** | Demonstrably helps on large blocks (62-var condenser, 34-var scroll_compressor). Worth keeping but should not be first choice due to overhead. |
| **Newton1D** | **Essential** | Specialized code for 1D blocks is very effective. The multi-probe + bisection approach is well-designed. Must keep. |

### 2.4 Is There Redundancy?

**Newton vs TrustRegion:** There is partial redundancy. Both use the same Jacobian factorization. The key difference is globalization strategy (line search vs trust region). In theory, TR should be more robust, but the current implementation isn't. Rather than removing TR, it should be improved — a well-implemented TR solver should outperform Newton on ill-conditioned problems.

**Newton vs LM:** LM is a different algorithm (damped least-squares). It's complementary to Newton. LM should converge from farther initial guesses but is slower near the solution. This is the right solver to keep as a fallback.

**Partitioned:** Not redundant — it's a fundamentally different approach (per-variable updates). However, its current failure on small blocks (partitionedMinBlockSize check) is a bug-like behavior since it returns InvalidInput instead of simply skipping gracefully in the pipeline.

---

## 3. Trust Region Implementation Issues

The current TR solver has several issues that explain its poor performance:

1. **Overly aggressive shrinking:** When evaluation fails, `delta *= trShrinkFactor(0.5)`. After a few failures, delta can reach 1e-8, triggering a full reset to `trInitialRadius`. This oscillation between too-small and too-large is pathological.

2. **No damping relaxation:** Unlike LM which smoothly transitions between gradient descent and Newton, TR has a hard boundary between Cauchy and Newton steps. For highly nonlinear thermodynamic systems, this can cause oscillation.

3. **Missing warm-restart of delta:** When TR is called as a fallback after Newton failure, it always starts with `trInitialRadius=10`. It should inherit a smaller radius if Newton already found a partial solution.

4. **No bisection or bracket exploitation:** Unlike Newton1D, the multidimensional TR doesn't maintain any bracket information. A sign-change along any coordinate direction is discarded.

---

## 4. Multidimensional Bisection: Feasibility Analysis

### 4.1 Classical Multidimensional Bisection (Eiger-Sikorski-Stenger, Vrahatis)

These methods generalize 1D bisection to $\mathbb{R}^n$ by maintaining a simplex/polyhedron with vertices that have specific sign patterns for each component of $F(x)$. Each bisection step evaluates $F$ at one new point and replaces a vertex.

**Key constraint:** The initial simplex requires $O(2^n)$ (Eiger) or $O(2n)$ (Vrahatis) vertices with specific sign patterns. For CoolSolve's block sizes (28-62 variables), finding such a simplex is the fundamental bottleneck.

| Method | Initial vertices | Cost per step | Practical for n=28-62? |
|---|---|---|---|
| Eiger-Sikorski-Stenger | $2^n$ | 1 F-eval | **No** ($2^{28} \approx 10^8$) |
| Vrahatis characteristic bisection | $2n$ | 1 F-eval | **Marginal** (finding characteristic polyhedron is hard) |
| Interval Newton (Kearfott/Hansen) | 1 box | $O(n^3)$ per step | **Only for verification** |

**Conclusion:** Pure multidimensional bisection is **not practical** for large problem sizes. The exponential scaling with dimension and the difficulty of constructing initial simplices with proper sign patterns make it infeasible for blocks with more than ~10 variables. It should be implemented, but with a check function that makes sure that it is not applied to blocks that are too large.

### 4.2 What Works most likely: Homotopy Continuation

Homotopy/continuation methods are the practical alternative that provides global convergence guarantees without the curse of dimensionality. They work for any $n$ with $O(n^3)$ cost per step (same as Newton).

#### Probability-One Homotopy

**Core idea:** Construct a parameterized family of problems:

$$H(x, t) = t \cdot F(x) + (1-t) \cdot G(x) = 0$$

where $G(x) = x - x_0$ has a known solution at $x = x_0$. At $t=0$, the solution is trivially $x_0$. At $t=1$, we recover our target $F(x)=0$. By tracking the solution curve from $t=0$ to $t=1$, we reach a root of $F$.

**Theoretical guarantee (Chow-Mallet-Paret-Yorke, 1978):** For almost all starting points $x_0$, the homotopy path is a smooth curve that reaches a solution of $F(x) = 0$. This is a **global convergence guarantee**.

**Implementation plan for CoolSolve:**

1. **Phase 1 — Simple continuation (natural parameter):**
   - Step $t$ from 0 to 1 in increments $\Delta t$
   - At each step, use Newton (existing solver) with the previous solution as initial guess
   - Adaptive step-size: reduce $\Delta t$ if Newton fails, increase if it converges quickly
   - Estimated effort: ~200 lines of C++
   
2. **Phase 2 — Pseudo-arclength continuation:**
   - Parameterize by arclength $s$ instead of $t$
   - Solves an $(n+1) \times (n+1)$ bordered system per step
   - Handles turning points (fold bifurcations)
   - Estimated effort: ~300 additional lines

3. **Phase 3 — Engineering-specific homotopies:**
   - Instead of the generic $G(x) = x - x_0$, use physically meaningful transformations:
     - Ideal gas → real gas (gradually introduce CoolProp)
     - Simplified geometry → full geometry
     - Low pressure → operating pressure
   - These "physics-aware" homotopies tend to have smoother paths

#### Available Libraries

| Library | Language | Integration effort | Notes |
|---|---|---|---|
| HOMPACK90 | Fortran 90 | Medium (extern "C") | Classic, well-tested, probability-one guarantee |
| LOCA (Trilinos/NOX) | C++ | High | Production-quality but heavyweight dependency |
| Custom implementation | C++ | Low | Uses existing Newton solver as corrector, ~200-500 LOC |

**Recommendation:** Start with a custom implementation (Phase 1) since it reuses existing Newton/LM solvers as the corrector. This gives immediate value with minimal code. If more robustness is needed, integrate HOMPACK90 or implement Phase 2.

### 4.3 Implementation Plan: Homotopy Solver for CoolSolve

#### New class: `HomotopySolver`

```cpp
class HomotopySolver : public NonLinearSolver {
    // Wraps any existing solver as the "corrector"
    SolverStatus solve(Problem& problem, VectorXd& x, const SolverOptions& opts, ...) {
        // Phase 1: Natural parameter continuation
        VectorXd x0 = x;  // Starting point
        double t = 0.0, dt = 0.1;
        
        while (t < 1.0) {
            t = min(t + dt, 1.0);
            
            // Build H(x, t) = t*F(x) + (1-t)*(x - x0)
            Problem homotopyProblem = wrapHomotopy(problem, x0, t);
            
            // Solve with existing Newton/LM solver
            auto status = corrector_->solve(homotopyProblem, x, opts);
            
            if (status == Success) {
                dt = min(dt * 2.0, 0.2);  // Increase step on success
            } else {
                t -= dt;
                dt *= 0.5;  // Reduce step on failure
                if (dt < 1e-6) return MaxIterations;
            }
        }
        return Success;
    }
};
```

#### Integration points:
1. Add `SolverStrategy::Homotopy` enum value
2. Add it to the default pipeline: `Newton, TrustRegion, LM, Homotopy, Partitioned`
3. Homotopy would be tried after the standard solvers fail, before Partitioned

---

## 5. Unit Test Improvements

### 5.1 Current Testing Gaps

- **No per-solver regression tests:** If a solver implementation changes, there's no way to detect if models that previously solved now fail.
- **No performance benchmarks:** Solve times aren't tracked.
- **No "without initials" CI tests:** Only the robustness test (manual, slow) covers this case.

### 5.2 Proposed Testing Strategy

#### Tier 1: Fast CI tests (run on every commit)
- **Existing `[examples-comprehensive]` tests** — keep as-is, validate solutions
- **Add `[solver-regression]` test:** For each model that currently solves, assert it still solves with the default pipeline. This catches regressions quickly.

#### Tier 2: Solver-specific regression tests (run weekly or on solver changes)
- For each solver strategy, maintain a list of models it should solve:
  ```
  Newton:    {all 17 models with initials}
  TR:        {all 17 models with initials}  
  LM:        {15 models with initials}  (excluding orc_extraction, scroll_compressor)
  Partitioned: {10 models with initials}
  ```
- When a solver is improved, update the expected-pass list.

#### Tier 3: Full robustness test (run manually or nightly)
- The improved `[solver-robustness]` test with the Markdown report
- Compare reports across versions to track progress

#### Tier 4: Without-initials regression
- Track which models solve without initials per solver
- This is the **key metric for robustness improvement** — the gap between "with initials" (94%) and "without initials" (82%) is what needs to close.

### 5.3 Concrete Changes to Test Files

1. **`test_solver_robustness.cpp`** — Already improved (shorter column labels, error categories in cells, model difficulty ranking, error category breakdown). The test now runs with a 30s timeout per solve.

2. **New `test_solver_regression.cpp`** — A fast test that asserts the default pipeline solves all expected models. Runs in CI:
   ```cpp
   TEST_CASE("Solver regression: all models with default pipeline", "[solver-regression]") {
       // For each model in examples/, run with default pipeline + initials
       // REQUIRE(result.success) for the expected-pass list
   }
   ```

3. **Solver-specific expected-pass lists** — Codify the current results as expected Pass/Fail baselines. When a solver is improved, update the list and the CI will enforce the new baseline.

---

## 6. Prioritized Roadmap

### Phase 1: Quick Wins (1-2 days)

1. **Extract `computeScalingFactors()` as a shared utility** — eliminate code duplication across Newton, TR, LM.

2. **Fix Partitioned solver's InvalidInput on small blocks** — instead of returning InvalidInput when block < partitionedMinBlockSize, return MaxIterations or skip gracefully in the pipeline. This eliminates the misleading "Other" error category.

3. **Add the `[solver-regression]` test** — codify current pass/fail baselines.

### Phase 2: Solver Improvements (1-2 weeks)

4. **Improve TrustRegion implementation:**
   - Fix the delta shrink/reset oscillation pattern
   - Add a smooth transition between Cauchy and Newton (quadratic interpolation)
   - Inherit delta from previous solver in the pipeline
   - Target: TR should match Newton's success rate (94.4%) with initials

5. **Improve LM implementation:**
   - Use the Marquardt diagonal scaling variant more aggressively
   - Add a geodesic acceleration term (Nielsen improvement)
   - Target: LM should solve at least 17/18 with initials

6. **Implement simple Homotopy continuation solver:**
   - Natural parameter continuation using Newton as corrector
   - Adaptive step-size control
   - Add to pipeline as fallback before Partitioned
   - Target: solve orc_co2 and humidair2 without initials (currently 82.4% → ≥90%)

### Phase 3: Architecture Refactoring (1 week)

7. **Split `solver.cpp` into per-solver files** (as described in §2.1)

8. **Extract Newton1D as a standalone class** — it's currently inlined in `Solver::solveBlock()` at ~490 lines. It should be its own class conforming to the NonLinearSolver interface (or a specialized 1D interface).

9. **Create shared convergence-checking and trace-recording utilities**

### Phase 4: Advanced Solvers (2-4 weeks)

10. **Pseudo-arclength continuation** — extend the homotopy solver to handle turning points.

11. **Physics-aware homotopy** — create homotopy paths that follow physically meaningful parameterizations (ideal→real gas, low→high pressure).

12. **SUNDIALS KINSOL integration** — for large blocks (>30 variables), KINSOL's preconditioned Newton with Anderson acceleration may outperform our custom solvers.

### Phase 5: Validation & Hardening

13. **Add more challenging test models** — models with phase transitions, near-critical conditions, and large algebraic loops.

14. **Automated nightly regression** — run full robustness test, diff against baseline, flag regressions.

15. **Solution verification** — use interval Newton to verify that converged solutions are locally unique.

---

## 7. Summary

| Aspect | Current | Target |
|---|---|---|
| With initials, default pipeline | 94.4% (17/18) | 100% (fix orc_complex parse) |
| Without initials, default pipeline | 82.4% (14/17) | ≥94% (with homotopy) |
| Hardest model (scroll_compressor) | 5/12 configs succeed | ≥10/12 |
| Solver file size | 2850 lines in one file | Split into ~9 focused files |
| Regression testing | Manual only | Automated per-solver baselines |
| Global convergence method | None | Homotopy continuation |

The **single highest-impact change** is implementing a homotopy continuation solver. It directly addresses the 82.4%→94%+ gap for without-initials solving, and the theory guarantees it will converge from (almost) any starting point. This is more impactful than fixing TR or LM, which are incremental improvements.
