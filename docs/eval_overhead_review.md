# Review notes: evaluation overhead and time marching

[← Back to overview](../README.md)

Branch `perf/eval-overhead`, from `a215c3c` (the plan of
[Solver Roadmap §6.1](solver_roadmap.md#61-evaluation-overhead-and-time-marching-november-2026),
where the results are). Not committed, not merged: the changes touch the
shared solver core and are for the maintainer to review.

## 1. Files touched

| File | Core or module | Change |
|------|----------------|--------|
| `include/coolsolve/evaluator.h`, `src/evaluator.cpp` | **core** | `VariableValues` / `StringValues` (case-insensitive maps); `BlockEvaluator::evaluate` takes them by reference, binds the external variables of its equations once and reads any other name from the map (`ExpressionEvaluator::setExternalValues`); explicit assignments without gradient (`isExplicitAssignment`); routine tables referenced (`setRoutines`, copy on a later `registerFunction`); `findCoolPropOutputParam` (no exception); `callsPropertyFunctions`; `takeInputClamped` loads before it exchanges |
| `include/coolsolve/solver.h`, `src/solver.cpp` | **core** | `solveBlock` tries the explicit evaluation first and uses the state by reference (the five copy loops are gone, also those of the symbolic reduction); signatures retyped; `solve()` split into `effectiveOptions`, `solveSystemBlock`, `blockFailureMessage`, shared with the new `solveBlocks` |
| `src/solver_newton1d.cpp` | **core** | overload taking the state and the block variable: the probes of Phase 2 are built from it only when Phase 2 is reached, with the former content and order |
| `src/runner.cpp` | core (runner) | CoolProp warm-up after the IR, only when `callsPropertyFunctions`; `integral.md` prints the method used |
| `include/coolsolve/integral/integral_solver.h`, `src/integral/integral_solver.cpp` | INTEGRAL module | march blocks (`marchBlocks_`), `solveBlocks` per stage, rows and final result read from the state, `IntegralSolveResult::method`, hint for a step too large |
| `src/integral/integrator_euler_implicit.cpp` | INTEGRAL module | a step that does not converge is an error |
| `tests/test_evaluation_overhead.cpp` (new), `tests/test_evaluator.cpp` | tests | 13 targeted test cases `[eval-overhead]`; the external maps of 4 tests retyped |
| `docs/solver_roadmap.md`, `docs/integral_table.md`, `docs/debugging_models.md`, `docs/language_reference.md` | docs | results, march strategy, new messages |

No `SolverOptions` field, configuration key, GUI or REST change.

## 2. Why other models are unaffected

- **Same values seen by every equation.** Before, a block evaluation copied
  the whole state (minus the block's variables) into its evaluator; now it
  binds the external variables its equations list and reads any other name
  from the state. The lookup order is unchanged: local value (block variables
  first), then the state, then the built-in constants, then the default 1 with
  its warning. The bodies of FUNCTIONs and PROCEDUREs still see only their
  arguments.
- **No change of snapshot semantics.** Passing the state by reference instead
  of a snapshot would differ only if the state changed during a block solve.
  Within `solveBlock` it is written only where the snapshot was updated too
  (symbolic pre-solve) or for the block's own variables, which `x` shadows.
  The parallel multi-start keeps its copy, because a winner writes the state
  while the other candidates of its wave may still be evaluating.
- **Explicit blocks.** They are evaluated without the Jacobian: the value is
  computed by the same code (`evaluate` instead of `evaluateTerms`, documented
  as bit for bit), and with no gradient the clamping of a property input sees
  "does not depend on the unknowns", as it did with a zero gradient. The
  explicit test is the former one, now in the `BlockEvaluator` constructor
  (names compared with the EES folding, which also folds accented letters).
- **Newton1D probes.** The same map, in the same order, excluding the block
  variable with the same ASCII comparison, built only when Phase 2 runs.
- **`solve()`** runs the same per-block steps in the same order, with the same
  messages; only its structure changed.
- **INTEGRAL march.** The blocks that depend neither on `t` nor on a state are
  solved by the complete solve at `t0` and keep their values; before, they
  were solved again at every stage from their converged value. Only the
  variables and strings of the final `algebraicResult` are used downstream
  (`runner.cpp`), and they are now the state after the last row.
- **Behaviour changes, intended:** a non-converged `EulerImplicit` step is an
  error (it was accepted); verbose mode (`verbose = true`) no longer prints the
  externals of the explicit blocks (now solved before that printout); a model
  without a property call no longer loads CoolProp, so `profiling.md` shows no
  warm-up line for it.

## 3. Baselines

| Check | Stock (`a215c3c`) | This branch |
|-------|-------------------|-------------|
| Three greenhouse models (full runs) | – | `-integral.csv` and `.sol` identical byte for byte |
| `EulerImplicit` at 120 s on the climate model (6 h) | – | table identical; 600 s now fails at t = 0 (residual 232) |
| Library model `air_nozzle_ejector` (`EulerImplicit`) | – | `.sol` identical |
| `./coolsolve_tests` | all pass | all pass (1334 + 13 new test cases) |
| `[examples-comprehensive]` | 46 OK, `cooling_tower2` and `turbocompressor_interpolate` fail, `orc_complex` does not parse | identical report except the times (5.09 s → 0.90 s in all) |
| `[solver-robustness]` | not run (too slow, at the maintainer's request) | – |
| `[library-regression]` | not run | – |

The assertion count of `./coolsolve_tests` differs from stock by a few:
`tests/test_fluids.cpp` loops over `FluidRegistry::getAllFluids()`, a set
ordered by pointer, inside `SECTION`s, so its count depends on the heap
history of the run. Both binaries give the same count when the same test
files run (537 assertions for `test_config`, `test_evaluator`,
`test_examples`, `test_fluids`).

## 4. Open points

- The library regression is the check named by the plan for CALLs, user
  functions, string variables and array elements with a computed index; it
  has not been run on this branch.
- Step 5 (compiled block evaluation) is to be decided on the results in
  Roadmap §6.1.
