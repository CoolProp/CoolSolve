# §3.2.2 — Post-Reduction Block Re-Decomposition

**Status:** Design document  
**Depends on:** §3.2.1 Symbolic Block Reduction (`enableSymbolicReduction`)  
**Authors:** Auto-generated design analysis  
**Date:** 2026-03-03

---

## 1. Problem Statement

The current symbolic reduction pipeline (§3.2.1) removes equations and
variables from strongly connected blocks via CoolProp inversions, explicit
extractions, and dead-variable elimination.  After removal, the *remaining*
equations are solved as a **single monolithic block**.

However, removing a node from a strongly connected component (SCC) can
**break** it into multiple smaller SCCs.  If that happens, the current code
misses an opportunity: instead of solving one `n`-variable nonlinear system,
we could solve two or more independent systems of sizes `n₁ + n₂ + … = n`,
which is dramatically easier (Newton cost scales as O(n³) per iteration).

**The question is: does this actually happen in practice, and what is the
best strategy to exploit it?**

---

## 2. Theoretical Analysis

### 2.1 Can Removing Nodes from an SCC Create Multiple SCCs?

**Yes.** This is a fundamental property of directed graphs.

A strongly connected component (SCC) is a maximal set of nodes where every
node can reach every other node through directed paths.  If we remove a
node `v` that is an *articulation point* — meaning it lies on the only
directed path between some pair of other nodes — then the remaining graph
may have two or more SCCs.

**Graph-theory argument:**

Consider a block with equations `{E1, …, E5}` coupled through variables
`{a, b, c, d, e}`:

```
E1(a, b, c)  →  depends on E3 via c
E2(a, b)     →  depends on E1 via a
E3(c, a)     →  depends on E1 via a;  matched output = c
E4(d, e, c)  →  depends on E3 via c
E5(d, e)     →  depends on E4 via d
```

The dependency graph (eq → eq edges) is:
```
E1 → E3 → E1  (cycle through c, a)
E3 → E4 → E5 → E4  (cycle through d, e)
E1 → E2 → E1  (cycle through a, b)
```
All nodes can reach each other through E3 (which bridges via variable `c`),
so `{E1..E5}` is one SCC.

Now suppose symbolic reduction inverts E3 away (CoolProp inversion: we can
compute `c` directly).  Remove E3 and `c` from the block.  The remaining
graph is:

```
E1(a, b) ↔ E2(a, b)     — SCC₁ = {E1, E2}
E4(d, e) ↔ E5(d, e)     — SCC₂ = {E4, E5}
```

No edge connects the two groups because `c` was the only bridging
variable.  **The single SCC has split into two independent SCCs.**

### 2.2 Is This Likely in Thermodynamic Models?

**Very likely.**  CoolProp property calls are often the "glue" that
couples sub-systems.  Typical patterns:

| Physical model             | Coupling variable      | Bridge equation                           |
|----------------------------|------------------------|-------------------------------------------|
| Compressor ↔ condenser     | `h_ex` (enthalpy)     | `h_ex = enthalpy(R$, T=T_ex, P=P_ex)`    |
| Evaporator ↔ expansion     | `T_su` (temperature)  | `T_su = temperature(R$, h=h_su, P=P_su)` |
| Sub-cooler ↔ two-phase     | `P_r_ex` (pressure)   | `t_out = temperature(R$, x=0, P=P_r_ex)` |

When such an equation is inverted (e.g., `T_ex = temperature(R$, H=h_ex,
P=P_ex)`), the variable that bridged two physical sub-systems disappears
from the remaining equations, and the sub-systems become independent.

### 2.3 Concrete Example: `condenser_3zones.eescode`

Block 38 has 62 equations.  After symbolic reduction, 6 variables are
eliminated (5 CoolProp inversions + 1 extraction), leaving 56 equations
solved as one monolithic system.

The removed variables include:
- `h_r_ex_cd_sh` — enthalpy at superheated zone exit
- `P_r_ex_cd_sh` — pressure at superheated zone exit (CoolProp inversion)
- `h_cf_ex_cd_sh`, `h_cf_ex_cd_tp`, `h_cf_ex_cd_sc` — coolant enthalpies
- `P_r_cd` — average refrigerant pressure

These variables are shared between the three heat exchange zones
(superheated, two-phase, sub-cooled).  Their removal may decouple these
zones into independent sub-blocks, reducing a 56-variable problem into
several smaller ones.

**Whether this block actually splits can only be determined by running
Tarjan's SCC algorithm on the 56-equation reduced graph.** This is exactly
the proposed improvement.

---

## 3. Proposed Strategy: Lightweight In-Block Re-Decomposition

### 3.1 High-Level Design

After symbolic reduction produces a set of reduced equations and variables,
**run a local BLT decomposition on that subset** before handing it to the
nonlinear solver.  If the reduced graph has multiple SCCs (or acyclic
components), solve them sequentially in topological order.

```
Original pipeline:

  StructuralAnalysis → blocks → for each block:
    SymbolicReduction → reducedBlock(nReduced)
    Solve(reducedBlock)

Proposed pipeline:

  StructuralAnalysis → blocks → for each block:
    SymbolicReduction → reducedBlock(nReduced)
    LocalBLT(reducedBlock) → sub-blocks [b₁, b₂, …, bₖ]
    for each sub-block bᵢ (topological order):
      Solve(bᵢ)
      Propagate solved variables to next sub-blocks
```

### 3.2 Algorithm: `redecomposeReducedBlock()`

Input:
- `reducedEquationIds` (vector of global equation IDs)
- `reducedVariables` (vector of variable names)
- Reference to the IR (for `EquationInfo.variables`)
- The global `matching` (for output-variable assignments)

Steps:

1. **Build local matching.**  
   For each reduced equation `eqId`, its output variable is
   `analysis.matching[eqId]`.  Build `localOutput[i] = matching[reducedEqIds[i]]`
   for local index `i`.

2. **Build local adjacency list.**  
   For each local equation index `i`, scan its variable set
   (`EquationInfo.variables`).  For each variable `v` in that set:
   - Skip if `v` is the output of equation `i`
   - Skip if `v` is not in `reducedVariables`
   - Find local equation index `j` such that `localOutput[j] = v`
   - Add edge `i → j`

   This is O(k²) in the worst case for k = |reducedEquationIds|, but in
   practice each equation references only a few variables, so it is O(k·m)
   where m is the average number of variables per equation (~3–5).

3. **Run Tarjan's SCC.**  
   Call `StructuralAnalyzer::tarjanSCC(k, localAdj)`.
   Returns sub-SCCs in topological order (dependencies first — this is
   already the convention in the codebase, see `structural_analysis.cpp`
   line 152).

4. **Build sub-blocks.**  
   For each SCC, create a `Block` with:
   - `equationIds` = global equation IDs of the SCC members
   - `variables` = output variables of those equations

5. **Return** `vector<Block>` in topological solve order.

**Complexity:** O(k + E) for Tarjan, where k ≤ 60 and E ≤ 300 typically.
Negligible compared to even a single Newton iteration.

### 3.3 Integration into `solveBlock()`

The change is localized to the `else` branch in `solveBlock()` (the
`nReduced > 1` case at ~line 1096 in `solver.cpp`).

**Current code (simplified):**
```cpp
} else {
    // Reduced block still size > 1 — use normal pipeline
    BlockEvaluator reducedEval(reducedBlock, ir_, options.coolpropConfig);
    // ... build problem, solve monolithically
}
```

**Proposed replacement:**
```cpp
} else {
    // Re-decompose reduced block into sub-SCCs
    auto subBlocks = redecomposeReducedBlock(
        reduction.reducedEquationIds,
        reduction.reducedVariables,
        ir_, analysis_);

    if (subBlocks.size() <= 1) {
        // Still a single SCC — solve monolithically (existing code)
        BlockEvaluator reducedEval(reducedBlock, ir_, options.coolpropConfig);
        // ...
    } else {
        // Multiple sub-SCCs — solve sequentially
        for (auto& subBlock : subBlocks) {
            size_t subN = subBlock.variables.size();
            // Build external vars (everything not in subBlock)
            auto subExternal = buildExternalVars(subBlock);
            if (subN == 0) continue;
            if (subN == 1) {
                // Explicit or Newton1D — same logic as existing size-1 path
            } else {
                // Full pipeline solve for this sub-block
            }
            // Propagate solved variables for subsequent sub-blocks
            for (const auto& var : subBlock.variables) {
                reducedExternalVars[var] = evaluator_.getVariableValue(var);
            }
        }
    }
}
```

### 3.4 Trace / Debug Output

When re-decomposition splits a block, record it in the `SolverTrace`:

```
Block 38: symbolic reduction 62 → 56 vars
  Re-decomposition: 56-var block → 3 sub-blocks [22, 18, 16]
  Sub-block 38a (22 vars): Newton SUCCESS (2 iter)
  Sub-block 38b (18 vars): Newton SUCCESS (3 iter)
  Sub-block 38c (16 vars): Newton SUCCESS (1 iter)
```

The `symbolic_reduction.md` debug file should include a section like:

```markdown
### Block 38: Re-Decomposition

After symbolic reduction (62 → 56 vars), the reduced block was decomposed
into 3 independent sub-blocks:

| Sub-block | Size | Variables (preview) | Solver | Iterations |
|-----------|------|---------------------|--------|------------|
| 38.1      | 22   | T_r_su_cd_sh, ...   | Newton | 2          |
| 38.2      | 18   | P_r_su_cd_tp, ...   | Newton | 3          |
| 38.3      | 16   | t_r_ex_cd_sc, ...   | Newton | 1          |
```

### 3.5 Activation Control

Re-decomposition is automatically enabled when `enableSymbolicReduction`
is true.  It is a natural extension of symbolic reduction — it only makes
sense when variables have been removed, and has zero cost when no
reduction occurred.  No separate option is needed.

**When symbolic reduction is disabled:** zero overhead.  The
re-decomposition code path is never reached because it sits inside the
symbolic reduction branch.

---

## 4. Implementation Plan

### 4.1 Files to Modify

| File | Change |
|------|--------|
| `include/coolsolve/structural_analysis.h` | Add static method: `redecomposeBlock(eqIds, vars, ir, analysis) → vector<Block>` |
| `src/structural_analysis.cpp` | Implement `redecomposeBlock()` using local adjacency + `tarjanSCC()` |
| `include/coolsolve/solver.h` | Extend `SolverTrace` and `BlockResult` with sub-block info |
| `src/solver.cpp` | In `solveBlock()`: call `redecomposeBlock()` when `nReduced > 1`, solve sub-blocks sequentially |
| `src/solver.cpp` | In `writeDebugReductionReport()`: add re-decomposition section |
| `src/solver.cpp` | In `writeDebugReductionReport()`: add re-decomposition section |

### 4.2 New Function Signature

```cpp
// In structural_analysis.h:
class StructuralAnalyzer {
public:
    /// Re-decompose a reduced block into sub-SCCs after symbolic reduction.
    /// Returns sub-blocks in topological solve order.
    /// If the reduced block is still a single SCC, returns a vector of size 1.
    static std::vector<Block> redecomposeBlock(
        const std::vector<int>& reducedEquationIds,
        const std::vector<std::string>& reducedVariables,
        const IR& ir,
        const StructuralAnalysisResult& analysis);
    // ...
};
```

### 4.3 Estimated Complexity

- **`redecomposeBlock()`:** ~50 lines (adjacency build + Tarjan call + Block construction)
- **`solveBlock()` changes:** ~40 lines (sub-block solve loop, replacing the current monolithic path)
- **Trace/debug extensions:** ~30 lines
- **Tests:** ~60 lines (2–3 test cases with synthetic blocks)

Total: ~180 lines of new/modified code.

---

## 5. Architecture Compatibility

### 5.1 Interaction with Tearing

#### 5.1.1 Current Execution Order

In `solveBlock()`, the three optimisations are applied as **sequential
fallbacks**, not complementary passes:

```
Step 1  →  Symbolic Reduction    (if enableSymbolicReduction && n > 1)
           On SUCCESS → return
           On FAILURE → fall through

Step 2  →  Tearing               (if enableTearing && n >= tearingMinBlockSize)
           On SUCCESS → return    operates on ORIGINAL block (not reduced)
           On FAILURE → fall through

Step 3  →  Standard Pipeline     (Newton → TR → LM → BisectionND → Homotopy → Partitioned)
```

**Critical observations:**

1. Tearing never sees the reduced block.  When symbolic reduction fails
   and falls through, tearing operates on the **full original block**
   (variables `varNames` and block evaluator `blockEval` are unchanged).
   
2. When symbolic reduction succeeds, tearing is skipped entirely.  It
   only activates as a fallback when reduction failed.

3. When *both* are disabled, only the standard pipeline runs.

#### 5.1.2 How Tearing Works

The tearing algorithm implements a **Greedy Feedback Vertex Set (FVS)**:

1. Build the block's internal dependency graph
   (`buildBlockDependencyGraph`, `structural_analysis.cpp` L218)
2. Find any cycle via DFS
3. Remove the node with maximum out-degree from the cycle
4. Repeat until the graph is acyclic
5. Removed nodes become **tear variables** (outer Newton loop)
6. Remaining nodes are solved in **topological order** (inner sequential
   evaluation) at each outer iteration

The torn system is solved via a Schur-complement Newton method: the inner
(acyclic) equations are solved by sequential forward-substitution, then a
Newton update is computed for the tear variables only.

#### 5.1.3 Impact of Re-Decomposition on Tearing

| Scenario | Current behaviour | With re-decomposition |
|----------|------------------|-----------------------|
| Reduction succeeds, block doesn't split | Tearing skipped | Tearing skipped (same) |
| Reduction succeeds, block splits into sub-blocks | N/A (not implemented) | Each sub-block solved independently; tearing only if sub-block ≥ `tearingMinBlockSize` |
| Reduction fails → tearing runs | Tearing on full original block | Same (tearing is a fallback) |
| Reduction succeeds, sub-block solve fails | N/A | Sub-block falls through to tearing on *that sub-block* |

**Key insight:** Re-decomposition generally **reduces or eliminates the
need for tearing**, because:

- Smaller sub-blocks have simpler cycle structures, requiring fewer (or
  zero) tear variables.
- Size-1 sub-blocks need no tearing at all.
- The total tear variables across all sub-blocks is typically ≤ the tear
  set of the original block, because some cycles in the original block
  spanned equations that now live in different sub-blocks.

**When tearing is still needed:** If a sub-block is large enough (≥
`tearingMinBlockSize`, default 4), it goes through the standard pipeline
which includes the tearing fallback.  Re-decomposition does not disable
tearing — it just makes it less necessary.

#### 5.1.4 Could Tearing Guide Which Equations to Reduce?

The current symbolic reduction picks CoolProp inversions in scan order,
without considering which removals have the most structural impact.  An
improvement would be to **prioritise inversions that are articulation
points** (bridge nodes) in the dependency graph — i.e., inversions that
would break the SCC into multiple sub-blocks.

This is a future enhancement (not in scope for this initial
implementation), but the design accommodates it: `redecomposeBlock()` uses
the same `tarjanSCC()` that could be called speculatively during the
reduction analysis to evaluate candidate inversions.

**Potential future algorithm:**

```
For each candidate CoolProp inversion (eq_i, var_i):
    Simulate removing eq_i from the dependency graph
    Run tarjanSCC on the remaining graph
    Score = number of resulting SCCs (more = better)
Sort candidates by score (descending)
Apply inversions in that order
```

This adds O(m · (V + E)) to the reduction analysis (m candidates, V+E
for each Tarjan run), which is negligible for practical block sizes.

### 5.2 Why Re-Decomposition Integrates Cleanly

The codebase already supports all the required primitives:

1. **`Block` is a plain value struct** (`equationIds` + `variables`).
   No special construction needed — we can build sub-blocks from any
   subset of equations and variables.

2. **`BlockEvaluator` accepts any `Block`**.  The solver already does
   this for reduced blocks (`reducedBlock` at solver.cpp ~line 1001).
   Sub-blocks work identically.

3. **`tarjanSCC()` is a static method** taking a generic adjacency list.
   No coupling to the global analysis — it can be called on any sub-graph.

4. **External variable propagation** is already handled: each solved
   block's variable values are added to `externalVars` before solving
   the next block.  The same pattern applies to sub-blocks.

5. **The matching is preserved.**  Symbolic reduction does not change the
   equation-to-variable matching from the original structural analysis.
   The output variable of each reduced equation is still
   `analysis.matching[eqId]`.

### 5.3 Zero-Overhead When Disabled

When `enableSymbolicReduction = false`:
- No code in `redecomposeBlock()` executes
- The re-decomposition logic is entirely inside the symbolic reduction
  branch, so there is zero overhead

### 5.4 Interaction with Tearing (Summary)

See §5.1 for the full analysis.  In brief:
- **Tearing** is a fallback applied after symbolic reduction fails.
- **Re-decomposition** is applied *within* symbolic reduction, before the
  solver pipeline.
- After re-decomposition, sub-blocks that are still large enough may
  still fall through to tearing if the standard pipeline fails.
- Re-decomposition generally reduces or eliminates the need for tearing.

---

## 6. Expected Impact

### 6.1 Performance

For a block of original size `n` that splits into `k` sub-blocks of sizes
`n₁, n₂, …, nₖ`:

- **Jacobian cost:** O(n₁³ + n₂³ + … + nₖ³) vs O(n³)
- Example: 56-var block → 3 sub-blocks of [22, 18, 16]
  - Before: 56³ = 175,616 (units)
  - After: 22³ + 18³ + 16³ = 10,648 + 5,832 + 4,096 = 20,576 (units)
  - **8.5× faster** per iteration
- Smaller blocks also converge in fewer iterations, so the real speedup
  may be **10–20×** for the block solve.

### 6.2 Robustness

Smaller blocks are:
- Less likely to have singular Jacobians
- Less likely to diverge from bad initial guesses
- Faster to fall through the solver pipeline if one strategy fails
- More amenable to BisectionND (which struggles above ~10 variables)

### 6.3 Risk

**Low.** The re-decomposition is a pure refinement — if the block does not
split, the behavior is identical to the current code.  If it does split,
we get smaller independent problems that can only be easier, never harder,
than solving them together.

The only risk is a bug in the local adjacency construction.  This is
mitigated by:
1. The function is simple (~50 lines, reusing existing `tarjanSCC`)
2. Easy to validate: sum of sub-block sizes must equal the reduced block
   size
3. Sits inside the symbolic reduction branch — disabling symbolic
   reduction also disables re-decomposition

---

## 7. Testing Strategy

### 7.1 Unit Tests

| Test | Description |
|------|-------------|
| `Redecomposition: block that splits into 2` | Synthetic block with bridge variable; after removing bridge, verify 2 sub-blocks |
| `Redecomposition: block that stays monolithic` | Synthetic fully-connected block; verify `redecomposeBlock()` returns 1 block |
| `Redecomposition: chain splits into 3` | Synthetic block with two bridge variables; after removing both, verify 3 sub-blocks |
| `Redecomposition: topological order` | Verify sub-blocks are returned in dependency order (causal blocks first) |
| `End-to-end: condenser_3zones with redecomposition` | Solve condenser_3zones.eescode with both flags on; verify success + smaller blocks in trace |

### 7.2 Robustness Tests

Run the existing robustness test suite with a new configuration:
```
"Default + SymbolicReduction + Redecomposition (with initials)"
```
and compare success rates against the existing `SymbolicReduction` config
to verify no regressions.

---

## 8. Summary

| Aspect | Details |
|--------|---------|
| **Correctness** | Mathematically sound: decomposing independent sub-systems does not change the solution |
| **Efficiency** | O(k + E) re-decomposition cost, negligible vs solver; potential 10× speedup for blocks that split |
| **Integration** | Uses existing `Block`, `BlockEvaluator`, `tarjanSCC` — no architectural changes |
| **Activation** | Automatic when `enableSymbolicReduction = true` (no separate option) |
| **Overhead when off** | Single `if` check per block — effectively zero |
| **Debug output** | Recorded in `SolverTrace`, shown in `symbolic_reduction.md` and terminal output |
| **Implementation size** | ~180 lines across 4 files |
| **Risk** | Low — pure refinement, identical behavior when block doesn't split |
