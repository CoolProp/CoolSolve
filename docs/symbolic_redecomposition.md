# Post-Reduction Block Re-Decomposition

**Feature area:** Symbolic reduction pipeline  
**Depends on:** `enableSymbolicReduction = true`  
**Automatic:** Yes — enabled whenever symbolic reduction is active  

---

## 1. Overview

When CoolSolve's symbolic reduction eliminates equations from a strongly
connected component (SCC), the remaining equations may no longer be mutually
dependent.  In that case the single large block can be **split into several
independent sub-blocks**, each solved separately and in dependency order.

This is called *post-reduction block re-decomposition*.  It is implemented
in `StructuralAnalyzer::redecomposeBlock()` and integrated directly into the
`solveBlock()` path in `solver.cpp`.

### Why it matters

Newton's method cost scales as O(n³) per iteration in block size *n*.
Splitting a 56-variable block into sub-blocks of sizes 42, 9, and 5
reduces the per-iteration cost from 175,616 units to 74,088 + 729 + 125 ≈
**75,000 units — a 2× improvement** in that step alone.  With fewer
iterations needed in smaller blocks, real-world speedups are typically
larger.  Smaller blocks also make convergence more reliable, especially for
bad initial guesses.

### Typical trigger

CoolProp property calls are often the "glue" that couples physically
distinct sub-systems inside one SCC.  When symbolic inversion removes such
a call (e.g., inferring enthalpy from a temperature-pressure pair), the
variable that bridged two sub-systems disappears, and the sub-systems
become independent.

---

## 2. How to Enable

Re-decomposition is **automatic** whenever symbolic reduction is enabled:

```conf
# coolsolve.conf
enableSymbolicReduction = true
```

or programmatically:

```cpp
SolverOptions opts;
opts.enableSymbolicReduction = true;   // re-decomposition is included
```

There is **no separate option**.  When `enableSymbolicReduction = false`
(the default), re-decomposition code is never reached, so there is zero
overhead.

---

## 3. What Happens at Runtime

For each multi-variable block, `solveBlock()` runs:

```
1. Symbolic reduction
   → eliminates equations via CoolProp inversion, explicit extraction,
     or dead-variable substitution
   → if the reduced block is empty (all eliminated): done, no Newton needed
   → if the reduced block has size 1: scalar Newton or explicit solve

2. Re-decomposition (only when reduced block size > 1)
   → calls StructuralAnalyzer::redecomposeBlock()
   → builds a local dependency graph on the remaining equations
   → runs Tarjan's SCC algorithm
   → returns sub-blocks in topological solve order

3a. If the block splits (sub-blocks > 1):
    → solve each sub-block independently in order:
        size 0: skip
        size 1: explicit assignment or Newton1D
        size > 1: full solver pipeline (Newton → TR → LM → BisectionND
                                        → Homotopy → Partitioned)
    → propagate solved values before each subsequent sub-block
    → on any sub-block failure: fall through to monolithic solve (step 3b)

3b. If the block does not split (single SCC):
    → solve as a single system with the full solver pipeline
    → this is the same path as without re-decomposition
```

---

## 4. Terminal Output

When `verbose = true` and a block splits, CoolSolve prints a summary line:

```
  Re-decomposition: 56-var block → 3 sub-blocks [42, 9, 5]
```

No output is produced when the block does not split, or when symbolic
reduction is disabled.

---

## 5. Debug Report (`symbolic_reduction.md`)

When `debugReductionPath` is set to a file path, CoolSolve writes a full
Markdown report.  The **Block Summary** table includes a `Sub-Blocks`
column:

```markdown
| Block | Original Size | Reduced Size | Sub-Blocks     | Inversions | ... |
|------:|-------------:|-------------:|---------------:|-----------:|-----|
|    38 |           62 |           56 | 3 [42,9,5]     |          5 | ... |
|    12 |           10 |            8 | —              |          2 | ... |
```

`—` means re-decomposition ran but the reduced block remained a single SCC.

The **Detailed Reduction Steps** section for each block also includes:

```markdown
**Re-decomposition:** 3 sub-blocks [42, 9, 5]
```

### Setting the debug path

```conf
# coolsolve.conf
debugReductionPath = symbolic_reduction.md
```

or:

```cpp
opts.debugReductionPath = "symbolic_reduction.md";
```

---

## 6. Programmatic Access

### `SolverTrace` (per-block tracing, requires `enableTracing = true`)

```cpp
struct SolverTrace {
    bool redecompositionApplied = false;
    int  numSubBlocks = 0;
    std::vector<int> subBlockSizes;   // size of each sub-block
    // ... (also: symbolicReductionApplied, originalBlockSize, reducedBlockSize, ...)
};
```

### `BlockResult` (always populated)

```cpp
struct BlockResult {
    bool redecompositionApplied = false;
    int  numSubBlocks = 0;
    std::vector<int> subBlockSizes;
    // ...
};
```

Example check:

```cpp
auto result = solver.solve(opts, /*enableTracing=*/true);
for (const auto& br : result.blockResults) {
    if (br.redecompositionApplied) {
        std::cout << "Block split into " << br.numSubBlocks << " sub-blocks: ";
        for (int sz : br.subBlockSizes) std::cout << sz << " ";
        std::cout << "\n";
    }
}
```

---

## 7. Interaction with Other Features

### Tearing (`enableTearing`)

Tearing is a separate fallback that runs **only when symbolic reduction
fails**.  Re-decomposition runs **within** symbolic reduction, before the
solver pipeline.

Execution order in `solveBlock()`:

```
Symbolic reduction (with re-decomposition)  →  if fails, fall through
Tearing                                      →  if fails, fall through
Standard pipeline (Newton → TR → … → Partitioned)
```

Re-decomposition and tearing are therefore independent and complementary:
- Re-decomposition converts one large system into several smaller ones.
- If a sub-block is still large and the standard pipeline fails, tearing
  may be applied to that sub-block as a last resort.
- Re-decomposition generally reduces or eliminates the need for tearing,
  because smaller blocks have simpler cycle structures.

### Initial guesses (`.initials` files)

Initial guesses are loaded before structural analysis and are available to
sub-block solvers in the same way as for full blocks.  There is no
interaction with initial guess loading.

### Solver pipeline mode (sequential vs parallel)

Each sub-block uses the configured solver pipeline
(`SolverOptions::solverPipeline`) in sequential mode.  Sub-blocks
themselves are always solved sequentially, in topological order.

---

## 8. Algorithm Details

`StructuralAnalyzer::redecomposeBlock()` takes the reduced equation IDs and
variable names and proceeds as follows:

1. **Build local matching.**  
   For each reduced equation `eqId`, the output variable is
   `analysis.matching[eqId]` (unchanged from the original structural
   analysis — symbolic reduction preserves the matching).

2. **Build local adjacency list.**  
   For each reduced equation *i*, scan its variable set from the IR
   (`EquationInfo::variables`).  For each variable *v*:
   - Skip if *v* is the output of equation *i* (self-loop)
   - Skip if *v* is not in the reduced variable set (external dependency)
   - Add an edge *i → j* where *j* is the equation whose output is *v*

3. **Run Tarjan's SCC.**  
   `tarjanSCC(k, localAdj)` is called on the local *k*-node graph.
   SCCs are returned in topological order (independent equations first,
   equations that depend on them later).

4. **Build output blocks.**  
   Each SCC becomes a `Block` with the corresponding global equation IDs
   and output variable names.

Complexity: O(k + E) where *k* ≤ reduced block size and *E* ≤ k × (average
variables per equation).  For realistic models this runs in microseconds.

---

## 9. Tests

The following tests cover this feature:

| Test | File | Tag |
|------|------|-----|
| `redecomposeBlock: trivial size-0` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp]` |
| `redecomposeBlock: trivial size-1` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp]` |
| `redecomposeBlock: two independent equations split into 2 sub-blocks` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp]` |
| `redecomposeBlock: mutually coupled equations stay monolithic` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp]` |
| `redecomposeBlock: chain of 3 splits into 3 sub-blocks in topological order` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp]` |
| `redecomposeBlock: mixed — 1 cycle + 1 independent forms 2 sub-blocks` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp]` |
| `redecomposeBlock: variables outside reduced set are ignored` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp]` |
| `Re-decomposition: SolverTrace fields populated when split occurs` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp][integration]` |
| `Re-decomposition: results match non-reduction baseline` | `test_symbolic_reduction.cpp` | `[symbolic][redecomp][integration]` |

Run them with:

```bash
./coolsolve_tests "[redecomp]"
```

The robustness test suite (`[solver-robustness]`) also tracks
re-decomposition statistics via the `redecompositionCount` and
`redecompDetail` fields in `RobustnessResult`.  These are reported in the
`Re-decomposition` column of the Symbolic Reduction Impact table.

---

## 10. Known Limitations

- **No articulation-point prioritisation.**  The current symbolic reduction
  picks inversion candidates in scan order, without preferring candidates
  whose removal would split the SCC.  An optimised selection strategy could
  achieve more splits, but is not yet implemented.

- **Sub-block failures fall back to monolithic.**  If any sub-block solve
  fails, the solver retries the entire reduced block as a single system.
  This is safe but loses the performance benefit.

- **Tearing is not applied at sub-block level** as a first-attempt
  strategy; tearing is only a last resort in the full pipeline.

---

## 11. Future Work

- **Articulation-point-guided inversion ordering:** during the symbolic
  reduction analysis, speculatively evaluate which inversion candidates are
  articulation points (bridge nodes) in the dependency graph.  Prioritise
  those candidates to maximise the number of sub-blocks produced.

- **Per-sub-block tearing as a first strategy:** for sub-blocks above a
  threshold size, apply tearing before Newton as a first attempt rather
  than as a last resort.
