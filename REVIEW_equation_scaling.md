# Review: the option `enableEquationScaling` (branch `feat/equation-scaling`)

Core solver change (register `CS-BUG-BLOCK-SCALE-SPREAD`): for review by the maintainer before any merge. Branch from
`main` @ `3c150bd`, committed on the branch after the supervisor's review and robustness run (§7).

## 1. What the option does, and what it does not

- `SolverOptions::enableEquationScaling` (bool, **default false**), coolsolve.conf key `enableEquationScaling`, the row
  counterpart of `enableScaling` (columns).
- When true, **every** Newton step of Newton, TrustRegion (its Newton step, before the dogleg) and KINSOL (line-search
  mode) is `solveEquilibrated(J, -F)`: rows of J scaled so that their largest entry is in [0.5, 1), then columns likewise,
  by powers of two (exact; LAPACK xGEEQU), then the same `ColPivHouseholderQR`. No trigger, no comparison of two
  steps (the rejected `a84f0e4` design). KINSOL's singularity test then uses the rank of the equilibrated matrix.
- When false, the step is the plain `ColPivHouseholderQR<MatrixXd> qr(J); qr.solve(-F)` of `main`.
- Unchanged in both cases: line search, merit function, convergence tests, variable scaling, Broyden updates.
  Not covered: LM (Marquardt-scaled normal equations), BisectionND, Partitioned, Newton1D (scalar), the outer step of
  the tearing solver (`FullPivLU` of the Schur complement). Homotopy's Newton corrector gets the option (its options
  are built field by field: `correctorOpts.enableEquationScaling` copied like `enableScaling`); every other inner solve
  copies the whole `SolverOptions` (multi-start, parallel pipeline, partitioned, reduced blocks).
- Deep Search (Try Harder) forces it on in `Solver::effectiveOptions`, as it forces tearing and symbolic reduction.

## 2. Files touched

Shared core:
- `include/coolsolve/solver_common.h`: `inversePowerOfTwo`, `solveEquilibrated` (restored from `a84f0e4`) and
  `solveNewtonStep(J, rhs, equationScaling, rank*)`, the one-branch selector.
- `include/coolsolve/solver.h`: the field and its comment; Deep Search comments.
- `src/solver.cpp`: conf key; `effectiveOptions` sets `enableEquationScaling = true` under `deepSearch`.
- `src/solver_newton.cpp`, `src/solver_trust_region.cpp`, `src/solver_kinsol.cpp`: the QR lines replaced by
  `solveNewtonStep(J, -F, options.enableEquationScaling[, &rank])`.
- `src/solver_homotopy.cpp`: one line, the corrector copies the option.

Other:
- Tests: `tests/test_library_gaps.cpp` (5 cases `[library-gaps][block-scale-spread]`, plus a comment of the existing
  `kAdiabaticFlameCh4`), `tests/test_config.cpp` (1 case `[config][solver][equation-scaling]`),
  `tests/test_solver_robustness.cpp` (field `enableEquationScaling` in `SolverConfig`, configurations
  "Default + EquationScaling (with initials)" and "(NO initials)", short label `+EqScale`).
- Config surfaces: `examples/coolsolve.conf` (commented key with when to use it), `gui/src/components/ConfigEditor.tsx`
  (entry beside `enableScaling`, group renamed "Scaling"; GUI rebuilt and re-embedded), comments in
  `gui/src/components/Toolbar.tsx` and `gui/src/api/client.ts`.
- Docs: `docs/language_reference.md` (§10 groups, new "Equation scaling" subsection, CHEM_EQUIL flame paragraph),
  `docs/debugging_models.md` (Try Harder list, CHEM_EQUIL flame remedy first), `docs/gui.md` (Try Harder row, Scaling
  group row), `docs/solver_roadmap.md` (Equation scaling row), `docs/model_library_support.md` (row moved §5 → §7,
  removed from the §9 queue, CSL-0150 follow-up), `README.md` (Try Harder, Convergence line).

## 3. The default path is unchanged

Guarantee: with the option off, `solveNewtonStep` factorises the same J with the same `ColPivHouseholderQR` and
solves the same right-hand side (−F is exact whether negated lazily or into a temporary), so every component of the
step is identical; the extra work is one branch per step. Evidence (Release build):

| Check | Result |
|---|---|
| Unit test "the Newton step is the plain QR step, bit for bit, when enableEquationScaling is off" | exact equality of the step and the rank on a well-conditioned, a singular consistent and an ill-scaled system |
| `./coolsolve_tests` | 1467 test cases, all passed (main: 1461; +6 new) |
| `[examples-comprehensive]` | 47/50 solving (fail: `cooling_tower2`, `orc_complex`, `turbocompressor_interpolate`), 47/47 correct, as the baseline |
| CLI on the 52 top-level examples, main build vs branch build (scratch copies) | same statuses, same Newton iterations (3866 in total), the 47 `.sol` files byte-identical |
| `[library-regression]` on library `a4ea304` vs `library_regression_final.md` | 219 PASS / 101 DEVIATION / 0 FAIL, 0 verdict changes, 0 detail changes; whole report identical except the times (sum 69.1 → 67.6 s, CSL-0326 58.9 → 57.3 s) |

## 4. Evaluation with the option on everywhere (default temporarily true, then reverted)

- **Unit suite**: 1465 / 1467; the 2 failures are the two assertions that the default is false.
- **`[examples-comprehensive]`**: report identical to the default run (verdicts, checked values, failed-block
  iterations); total time 0.82 → 0.91 s.
- **CLI on the 52 examples**: same statuses; Newton iterations 3866 → 3865 (`heat_pump_MSTh_SB_R10` 4 → 3).
  `.sol` changes: rounding level in 7 models (≤ 1.1e-9 relative), `exchangers2` at the tolerance level (A 7.350854 →
  7.350850, 5.6e-7), noise variables in `condenser_3zones` (ΔP ~1e-22), and `ice_storage_tank` (see §6).
- **`[library-regression]`**: 218 PASS / 102 DEVIATION / 0 FAIL. One verdict change: `CSL-0026` PASS → DEVIATION
  (`t_wall` −311.942 → −312.265, 0.1 %; `a`, `h_w`, `u` 2e-4–6e-4): a singular but consistent block with a family of
  solutions, of which the equilibrated QR picks another member (already seen with the always-equilibrated first
  version of `a84f0e4`). 33 detail changes with the same verdict: all at rounding level (≤ 4.6e-7) or on noise
  variables (`CSL-0008` pressure drops of ~0, `CSL-0142` balance residual `q_bal`). Sum of target times 67.6 →
  122.2 s, almost all `CSL-0326`: 57.3 → 111.7 s (+95 %; CLI: 11 → 15 iterations over all blocks, 60.0 → 107.5 s
  wall; its 3599-unknown block makes each extra Newton iteration cost several seconds). No other target changed by
  more than 0.25 s.
- **Native `CSL-0150`** (scratch copy of the export, `T = T_adiab_K` → `T = T_adiab`): off, *Block 17 (45 vars)
  failed to converge: LineSearchFailed* from the library's guesses (residual 6.3e6) and from `T_adiab = 1775`,
  `ntot = 12.5` (residual 0.32); on (its `coolsolve.conf`), SUCCESS in 10 iterations from the library's guesses,
  `T_adiab` = 1776.0848 °C (2049.2348 K), `ntot` = 12.4332, and in 30 iterations from 1775 / 12.5 (same T within
  5e-7 K). The W002 warnings list the mole fractions below the tolerance (`x_CH4` = −1.4e-20): in a block with terms
  of 1E9 they are only known to the absolute tolerance (existing behaviour of the convergence test).
- **Tests**: the flame with the 14 unknowns `H[i]` solves with the option (CALL NASA: 2048.2268 K; property functions:
  1776.0848 °C, equal to the formulation without `H[i]` within 1e-4 K) and with `deepSearch = true`.

## 5. Try Harder

Try Harder now also forces equation scaling: in its pipeline, Newton, TrustRegion, Homotopy's corrector and KINSOL
take equilibrated steps; LM, BisectionND and Partitioned are unchanged. It only runs after a failed solve, so the
default path is untouched; a model that failed for a reason unrelated to scale may take another path in Try Harder
than before (better or worse; the singular family of `CSL-0026` is an example of "another solution"). The
robustness suite has no Deep Search configuration; the two new "Default + EquationScaling" configurations show the
option on the default pipeline.

## 6. Not done, open questions

- Default later? On everywhere: no verdict change on the examples, one on the library (`CSL-0026`, singular family),
  and CSL-0326 +95 % time. Not a default candidate as it stands.
- LM, the tearing outer step and BisectionND are not covered; equilibrating the rows of J before LM's normal equations
  would be a separate change.
- KINSOL: with the option, a matrix that the plain QR finds rank-deficient only through scale (the 3×3 test system)
  is full rank, so KINSOL proceeds where it returned *SingularJacobian*; intended, but a behaviour change of that test.
- Cost when on: one n×n copy and O(n²) scaling per step beside the O(n³) QR; when off, one branch.
- `CSL-0150`: its library `coolsolve.conf` can set `enableEquationScaling = true` (maintainer's file).
- `examples/ice_storage_tank.eescode` was ill-posed (`U_dot_sto = M * u_dot` and `u_dot = U_dot_sto / M`, no energy
  balance: a family of solutions, `U_dot_sto` at its guess 1 with the option off, 65537 with it on). Fixed on this branch
  (maintainer request, 2026-11-16): the duplicate replaced by the first law `Q_dot_env_sto + H_dot_br_sto = U_dot_sto`, as
  in library model `CSL-0041`. The tank now melts in 30 419 s and warms to 5 °C: `E_discharged` = 4.247211e10 J, the closed
  form of the two phases to 7 digits, equal to `DELTAU_sto` (the former 6.56833e10 was the initial power times 50 000 s);
  identical with the option on; BDF 32 rejected of 169 steps (was 235 of 521, `docs/integral_table.md`).
- Library: `CSL-0150` updated in the library working tree (uncommitted): `T = T_adiab`, `coolsolve.conf` with
  `enableEquationScaling = true`, native `.sol`, verified (EES stored run: fuel oil, e = 0, 1 atm, T_adiab −0.21 %). It
  needs this branch: with `main` the key is unknown and the flame block fails, so the library change goes in with the merge.

## 7. `[solver-robustness]` (run by the supervisor, Release, report committed on the branch)

- Default path: the 30 configurations of `main` give the same result in every cell (1470 example × configuration cells,
  no OK → failure, no failure → OK; compared by configuration name, the two new columns excluded).
- The two new configurations, the default pipeline with `enableEquationScaling = true`:

| Configuration | Default pipeline | + EquationScaling |
|---|---:|---:|
| with initials | 44 / 46 | 44 / 46 (same examples) |
| NO initials | 36 / 45 | 34 / 45 |

  Without initials, `cooling_coil` (block 12, MaxIterations, residual 96) and `heat_pump_MSTh_SB_R10` (block 39,
  MaxIterations, residual 1e8) solve with the default pipeline and fail with the option; no example gains. The option is
  a remedy for scale-spread blocks, not a free improvement: this supports the default off.
- Try Harder: maintainer decision (2026-11-16): Try Harder keeps equation scaling forced on, to maximise the chance that a
  failing model starts working. Try Harder only runs after a failed solve; a model it solved before for another reason
  than scale may take another path (not measured: the robustness suite has no Deep Search configuration).