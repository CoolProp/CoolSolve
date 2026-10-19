# Review guide: the ammonia-water mixture (`NH3H2O`) on branch `feat/nh3h2o`

**Status: not merged, not pushed.** The work is on the local branch `feat/nh3h2o` (from `main` @ `a77efaa`). It is **not** to be merged into
`main` before the maintainer has reviewed the points of §3, because it touches the shared evaluator core (a few lines, listed
below). This page lists what changed, what is shared core and what is private to the new module, what guarantees that the other
models are unaffected, and how to check it.

## 1. What the branch adds

The EES fluid `NH3H2O` (three of `T`, `P`, `X`, `Q`, `H`, `S`, `U`, `V`; functions `Enthalpy`, `Entropy`, `IntEnergy`, `Volume`,
`Density`, `Temperature`, `Pressure`, `Quality`, `MassFraction`, `SpecHeat`, `VolExpCoef`) and the procedure
`CALL NH3H2O(Code, In1, In2, In3 : T, P, x, h, s, u, v, q)` (register `CS-GAP-NH3H2O`). User documentation:
[language_reference.md](language_reference.md#ammonia-water-nh3h2o). EES uses the model of Ibrahim & Klein (1993), whose
coefficients are not published; CoolSolve implements the **IAPWS Guideline G4-01 (Tillner-Roth & Friend 1998)** and says so: results
agree with EES within the difference of two formulations (about 1 %), with another reference state for `h`, `s`, `u`.

Commits (oldest first): `27c6813` the thermodynamic core and flashes, `db0246a` strict convergence (see §4), `6ecd20b` the
language integration, `b3f7c35` documentation, `3d2ca05` Windows installer script (unrelated to NH3H2O, own commit).

## 2. Files, by how much review they need

**Core solver files modified: none.** `git diff --name-only main` contains no `solver*.cpp`, no `solver_*.h`, no
`autodiff_node.*`, no `structural_analysis.cpp`, no `symbolic_reduction.cpp`, no `variable_inference.cpp`, no `ir.cpp`: the Newton
solver, the line searches, the block decomposition and the automatic differentiation of CoolSolve are untouched.

### A. New files, self-contained (no code path of any other model reaches them)

| File | Lines | Content |
|---|---:|---|
| `include/coolsolve/taylor_jet.h` | 331 | truncated Taylor polynomial in 3 variables (order 1 to 3), value and derivatives |
| `include/coolsolve/nh3h2o.h`, `src/nh3h2o.cpp` | 235, 401 | the Helmholtz formulation (Eqs. 1 to 8 and Tables 1 to 3, A1 of the guideline; IAPWS-95 for water) written once over a scalar type; the properties of a state |
| `src/nh3h2o_flash.cpp` | 1694 | the saturation solver, the 15 flashes, derivatives by the implicit function theorem. **Contains its own small Newton (`newtonSolve`), Brent and linear solver**, private to the file (anonymous namespace) |
| `src/evaluator_nh3h2o.cpp` | 303 | the glue: units of the model, designators, the trial-point limits, assembly of the `ADValue` gradient, the two entry points |
| `tests/test_nh3h2o.cpp` | 1139 | 29 test cases (tags `[nh3h2o]`, `[formulation]`, `[equilibrium]`, `[flash]`, `[evaluator]`) and two hidden scans (`[.nh3h2o-scan]`, about 35 s) |

### B. Shared files modified (the part to review): 11 files, about 100 lines (75 added, 24 removed)

| File | Change | Why other models are not affected |
|---|---|---|
| `src/evaluator.cpp` (+11/-7) | (1) in `evaluateProcedureCall`: `if (lowerName == "nh3h2o") { evaluateAmmoniaWaterCall(call); return; }` placed after the user-procedure lookup and next to the Colebrook hook; (2) in `evaluateCoolPropFunction`, right after the named inputs are evaluated: `if (fluid->getType() == FluidType::AmmoniaWater) return evaluateAmmoniaWaterFunction(...)`; (3) the two `getScale` lambdas of `evaluateCoolPropFunction` (the derivatives of the humid-air path and of the CoolProp path) now call `UnitConverter::scaleToSI` / `scaleFromSI` instead of writing `toSI(1) - toSI(0)` / `fromSI(1) - fromSI(0)` and the early return of 1 for `Dimensionless` themselves | (1) only reached for a `CALL` named `NH3H2O` when the model defines no such procedure; (2) only for the new fluid type, which no other fluid has. Both before any CoolProp code of the call, so nothing after them changes. (3) Behaviour identical: the helper evaluates the same expression with the same operands, so the factors are bit-identical; for `Dimensionless` it gives 1 - 0 = 1 exactly, as the early return did |
| `include/coolsolve/evaluator.h` (+6) | declares the two member functions | declarations only |
| `include/coolsolve/fluids.h` (+13) | new enum value `FluidType::AmmoniaWater` (appended before `Unknown`) and the class `AmmoniaWaterFluid` | an added enumerator; the existing ones keep their values. Every `switch` over `FluidType` was checked (`server.cpp`: two, edited) |
| `src/fluids.cpp` (+7/-4) | `NH3H2O` registered as `AmmoniaWaterFluid` instead of `UnsupportedFluid` | the one intended behaviour change: `NH3H2O` was an error, it now evaluates |
| `src/parser.cpp` (+1) | `"nh3h2o"` added to `knownBuiltinFunctions()` (as `"colebrook"`) | silences the "Unknown function" warning of `CALL NH3H2O`; the effect on the parser is that a user-defined `FUNCTION nh3h2o` would be taken for a built-in name in `isPositionalOnlyBuiltin` (named arguments refused); no known model |
| `include/coolsolve/units.h` (+17) | declares two static helpers of `UnitConverter`: `scaleToSI(type, unit)` = d(SI value)/d(value in unit) and `scaleFromSI(type, unit)` = d(value in unit)/d(SI value), with a Doxygen comment | declarations of two new functions; `toSI` and `fromSI` are unchanged |
| `src/units.cpp` (+8) | defines them: `toSI(1) - toSI(0)` and `fromSI(1) - fromSI(0)` (the slope of an affine unit conversion), a formula that was written seven times (twice in `src/evaluator.cpp`, five times in `src/evaluator_nh3h2o.cpp`) | new functions; called by the `getScale` lambdas of `src/evaluator.cpp` and by `src/evaluator_nh3h2o.cpp` only |
| `src/server.cpp` (+2) | the two `switch` over `FluidType` get a case | labels of the fluid list and of a message |
| `src/libr_ashrae.cpp`, `include/coolsolve/libr_ashrae.h` (-8) | the row `NH3H2O` of the table of known-but-not-implemented EES routines is removed (it would otherwise answer "not implemented" to a function call `NH3H2O(...)`) | a function call `NH3H2O(…)` is now an unknown function instead of a "not implemented" error; `CALL NH3H2O` is handled first |
| `tests/test_library_gaps.cpp`, `tests/test_priority2_features.cpp` | five assertions or titles in four existing tests that said "NH3H2O is unsupported" now say the new behaviour (two properties are refused with a message naming `NH3H2O`; the `CALL` tests of the other missing routines stay); three new `[library-gaps][nh3h2o]` test cases (94 lines) | the assertions of every other test are unchanged |

Hooks to read first: `git diff main -- src/evaluator.cpp src/fluids.cpp src/parser.cpp include/coolsolve/fluids.h include/coolsolve/evaluator.h`.

### C. Not NH3H2O, but on the branch

`build_installer.bat` and `docs/contributing.md` (Windows installer script, commit `3d2ca05`; not run on Windows). Documentation:
`docs/language_reference.md`, `docs/model_library_support.md` (register row moved from §4 to §7), `docs/ees_vs_coolsolve.csv`,
`README.md`, `gui/src/languages/ees.ts` (keywords).

## 3. Questions for the reviewer

1. **Trial-point limits** (`src/evaluator_nh3h2o.cpp`, `clampTrialPoint`): during the iterations an input that depends on the unknowns is
   limited (T below 150 K, P below 1000 Pa, X and Q outside [0, 1]), logged in the same `InputClampLog` as the CoolProp limits so that a
   limit active at the solution is refused. Same policy as CoolProp inputs; the bounds are new. Is this the wanted behaviour for X and Q?
2. **`FluidType::AmmoniaWater`** as a type of its own (rather than `Mixture`, whose code paths assume CoolProp mixtures). Alternative:
   fold it into `Mixture` with a flag. The separate type keeps all existing mixture code untouched.
3. **Variable inference / initial guesses**: `variable_inference.cpp` computes guesses with CoolProp from the fluid name; for
   `NH3H2O` (empty CoolProp name) the guess falls back to its default without an error in the models tried (see the tests); a
   guess from the module itself could be added.
4. The **two-root ambiguity** of the flashes `T H Q`, `T S Q`, `T U Q`, `P H Q`, `P S Q`, `P U Q` (the property along the saturation curve is
   not monotone in the composition): the root of lowest ammonia content is returned and documented. EES behaviour unknown.
5. **Thread safety**: the saturation solver keeps one warm-start state per kind in `thread_local` storage (`g_satCache`); no locking. The
   result does not depend on it: a call gives the same state (or the same error) in a fresh thread and after other calls (§4).
6. **Reference state**: IAPWS (u = s = 0 for the pure saturated liquids at their triple points), not EES's: documented as a difference,
   not corrected (differences at fixed composition agree).

## 4. The solver inside the module (the private Newton)

`newtonSolve` (`src/nh3h2o_flash.cpp`) is the only solver written for this work; it is not shared with the rest of CoolSolve.
**Convergence is strict**: a scaled residual below 1e-12 (pressure of a liquid on the scale `rho R T`, the scale of its round-off, measured 3e-13
of it; chemical potentials on `R * 300 K`). An earlier version had a loose fallback (step below 1e-9 and residual below 1e-7); it was
**removed** (`db0246a`) after it was found to accept under-converged states in a third of the iterations. The one remaining exit
for an iteration that has stopped lowering the residual requires every residual to be within two units of a **measured** noise floor
(`noiseFloor`: the resolution of the representation of the variables, `sum |J_ij| ulp(v_j)`, plus the jitter of the evaluation: the largest
deviation of F from its linearisation over six perturbations of the non-composition variables by 1 to 3 ulp, taken once, not the sum of the six samples
(review item 4; the measured exits did not change)). The test `Saturation: the equilibrium equations hold to 1e-11 of their natural scales`
recomputes the equilibrium residuals from the Helmholtz energy with none of the solver's code, from 230 K and up to mole fractions of
0.999999 (a hidden scan, `[.nh3h2o-scan]`, covers 230 to 560 K and the pressure-specified points).

**The water mole fraction without rounding** (review item 5). A mole fraction near 1 is a double of spacing 1.1e-16, so `1 - x` computed from
it is known to `1.1e-16 / (1 - x)` only, and the chemical potential of water inherits that error through the ideal mixing term
`(1 - x) ln(1 - x)`: 5.4e-7 of `R * 300 K` at the near-pure saturations, 9.5e-11 already for liquids up to 0.98 (their vapour is purer).
The formulation takes the complement as an argument, `reducedHelmholtz(T, rho, x, xc)` and the same for `idealHelmholtz`,
`reducingTemperature`, `reducingDensity` (the three-argument overloads pass `1 - x`: public API and Table 6 unchanged; the departure function
has no `1 - x`); for a jet the complement is `complementJet(x, xc)`, value `xc` exactly and the derivatives of `1 - x`. The solvers iterate a
composition as its logit `z` and evaluate `x = 1 / (1 + exp(-z))` with `xc = 1 / (1 + exp(z))`, each to its own round-off; `Sat`, the
warm-start cache and the two-phase `Solution` carry the complements, a given mass fraction `xi` gives its complement from `1 - xi`, and the
public `Saturation` returns them (`wL`, `wV`). Measured with the review's `exit_stats` (scans of the handover, §7.3):

| Section | Solved / threw | Non-strict exits before (HEAD `7d96d6d`) | After items 4 and 5 |
|---|---|---|---|
| Regular saturations (230 to 560 K, x 0.001 to 0.98) | 549 / 95 | 19 noise floor (worst 95 × tol) | 1 noise floor (1.8 × tol) |
| Near-pure saturations (x 0.99 to 0.999999) | 168 / 154 | 58 noise floor (3.9e5 × tol), 24 frozen (5.4e5 × tol) | 12 frozen, every row below tol |
| Pressure-specified saturations | 196 / 0 | 17 noise floor (2.8e3 × tol) | none |
| Round trips of the 15 flashes from two-phase states | 1523 / 7 | 5 183 noise floor (6.2e5 × tol), 935 frozen (3.4e5 × tol) | 121 frozen (15 above tol, worst 31 × tol) |
| (T, P, X) inside the dome | 935: all two-phase | 201 noise floor (1.1e3 × tol), 5 frozen (118 × tol) | 2 frozen (3.3 × tol) |

The chemical potentials recomputed independently with the complement (`exit_stats`, every saturation of the grids): `mu1` at most 1.0e-12 of
`R * 300 K` (was 5.4e-7), `mu2` 1.8e-12, p 1.3e-12. Recomputed with `1 - x` of the returned double instead, `mu1` still shows 5.4e-7: that
is the error of the subtraction, not of the state.

- **Frozen-variable polish** (`frozenPolish`, inside `newtonSolve`, with the composition-quantisation machinery `NewtonSetup::resolution` /
  `composition` and `logitResolution`): **still in the code**, because the measured exits above do not meet the condition the review set for
  its removal (no frozen-polish exit). It was written for the quantised composition: at a stall the most quantised variable and the row it
  dominates are frozen, the other rows are solved strictly, the frozen variable is moved by Newton on the dropped row (the Schur complement of
  the Jacobian), and the dropped row is accepted within two units of its quantisation floor `|S| u`. With the exact complement the logit is
  no longer quantised as `ulp(x) / (x (1 - x))` (`logitResolution`, up to 1e5 times the true resolution `ulp(z)`), so the polish now runs
  only as a rescue: it is called whenever one Newton step does not halve the residual, also far from the solution (seen at 1e9 × tol on T-p
  pairs of X = 0.99 at 270 K and on bubble points at p below 200 K), and its dropped-row criterion still uses that model. Measured without the
  polish and the machinery (the default resolution `ulp(z)`, `atResolutionFloor` as the only non-strict exit): the same solves succeed in
  every section of the table, every exit strict except the one noise-floor exit at 1.8 × tol, independent `mu1` at most 1.0e-12, the review
  cases A to E OK; the change is prepared, not applied.
- **Initial guesses** (`guessBubbleT`, `guessDewT`, `guessTP`, `guessAtPressure`): the liquid is evaluated from its fugacities against a vapour
  that is first an ideal gas, then corrected by the fugacity coefficients of the real vapour (they matter within a few percent of the pressure of
  pure ammonia, where the ideal-vapour guess has no root). `satSolve` tries the real-vapour guess first and the ideal-vapour one if Newton
  does not converge from it (at 16 MPa the correction over-corrects). Guesses never affect the accuracy: Newton polishes every result.

**Correctness of the saturation and of the flashes** (review of 2026-10-19, item 3; tests in the last section of `tests/test_nh3h2o.cpp`):

- **Stable roots only.** The equilibrium equations also have roots with a phase between its spinodals: the T-p equilibrium at the bubble pressure
  of X = 0.9 at 400 K gave a "liquid" of 21 109 mol/m3 with dp/drho = -8364 (the bubble pair: 27 144 mol/m3, +3205), and the bubble points at
  400 K of X = 0.55 to 0.9 came out unstable from a warm cache. `satNewton` now returns false unless dp/drho > 0 in both phases (the phases of
  the last Newton evaluation when it is the root, so no extra cost), and `satAttempt` tries its next start; after the polish of `finishTwoPhase`
  the same check is an error naming the state. No tangent-plane test. The trivial solution is rejected when |xL - xV| < 1e-4 and
  |ln(rhoL/rhoV)| < 1e-3 (was 1e-9 and 1e-6): near a critical point Newton stopped on two copies of one phase (2.8e-5 and 3.7e-8, a dew point at
  395 K and X = 0.999 placed 4 % above the true one); a real pair is that close only within about 1e-9 of the critical point.
- **A failed saturation is never a single phase.** `satAttempt` returns a status: `Found`, `NoGuess` (every guess function failed: the evidence
  of a state outside the two-phase region or above the critical locus) or `NotConverged` (a start existed, Newton failed from all of them or
  its root was rejected). The codes 123, 137 and 234/235 decide a single phase only on the side of a saturation point that was found, or from
  `NoGuess`; a `NotConverged` saturation that the decision needs is an error naming the state ("the dew point at p and x (...) did not
  converge to a pair of stable phases [state: NH3H2O(...)]"). Only the `NoGuess` message says "probably outside the two-phase region". The
  scans of the unknown-composition flashes keep the boolean `trySat` (a failed node is skipped).
- **Saturations inside a known bracket** (the T-p pairs of 123, 137, 138, 234/235, 238 and of the refinement of the unknown-composition flashes)
  are started from the bracket (the bubble and dew pairs interpolated in T or p) and a failure is an error that says the state is between its
  bubble and dew points (case D: 400 K, X = 0.999 inside the dome failed from the guesses).
- **(T, P, x) with one saturation point only** (above the critical temperature of x, near the critical locus, near a pure end): the T-p pair
  started from the point that was found decides two-phase if its compositions bracket x. The boundary margin `kBoundaryTol` (1e-7) now applies to
  123 as to 137 and 234/235, so that a state saturated to round-off has the same label whatever the history. 234/235 decide on the side of a
  found point whatever the status of the other one.
- **Pure components** (mole fraction at the floor 1e-10, X = 0 or 1): the mixture cannot represent their saturation pair from both sides (the
  bubble point of pure ammonia needs a vapour purer than the floor), so no saturation is solved and the rules of a state without a two-phase
  region apply, as before (lowest Gibbs energy for 123).
- **Starts.** In its near-pure-ammonia fallback, `guessTP` evaluates the liquid at p rather than at its ideal-vapour bubble pressure (the
  difference is of the size of the margin to the pure-ammonia pressure; cold failures of the T-p pair in a scan 230 to 405 K: 6 of 385, were
  102). When every guess of a bubble or dew point fails, the T-p pair at the T and p of the guess is the last start.
- **Bounds.** `finishTwoPhase` bounds rhoV below by 1e-6 mol/m3 like `satNewton` (was 1), `finishSingle` rho by 1e-6 (was 1): two-phase states
  below about 2 kPa (270 K, X = 0.1: rhoV = 0.88) and vapours below about 1 kPa whose start needs Newton steps now solve.
- **History independence**, measured: 2 311 flashes (the 15 codes from 212 reference states: two-phase 250 to 400 K, single-phase 280 to
  500 K) give the same result in a fresh thread and in one warm sequence (largest difference 2.2e-10 relative, identical errors; the previous
  code threw in 164 of them warm and 181 cold, now 47); so do 945 (T, P, x) states inside the dome (hidden scan) and 345 near the critical
  locus (380 to 450 K, X 0.7 to 0.999; the previous code: 58 differences, 31 single-phase answers inside the dome from a cold cache). The
  acceptance of the warm-start cache (`close()`, within 25 % of both specifications) is unchanged.

Reviewing `newtonSolve` (with `noiseFloor` and `frozenPolish`), `satNewton`, `satAttempt`, the phase decisions of `solveFlash`, `finishTwoPhase`
and `finishSingle` covers the numerics.

Comparison tolerances of the tests (against the guideline tables, EES, CoolProp) are separate from the solver tolerance and are stated in the tests:
Table 6 to 1e-7 relative (the table has 9 to 10 digits and the coefficients of the guideline 7 to 8), pure water to 1/60 000 as the guideline says,
the EES points to 0.2 % (volume) to 2 % (composition).

## 5. Verification and baselines

| Check | Before (`main` @ `a77efaa`) | Branch |
|---|---|---|
| Unit suite (`cd build && ./coolsolve_tests`) | 1182 test cases | 1214 test cases (+32), all pass, 21 041 assertions |
| Examples (`[examples-comprehensive]`) | 43/46 solve (`cooling_tower2`, `orc_complex`, `turbocompressor_interpolate` do not) | 43/46, the same three |
| Library regression (`tools/test_models.py`, scratch copies) | 249 targets, 5 expected deviations (`CSL-0008` ×3, `CSL-0146`, `CSL-0160:coolsolve`) | 249 targets, **the same 5** deviations, none new (run with the binary of `6ecd20b`, the evaluator integration; the later commits change only files of the module that no other model reaches) |

To reproduce: `cmake -S . -B build && cmake --build build -j && cd build && ./coolsolve_tests` and, in the library,
`python3 tools/test_models.py --coolsolve <build>/coolsolve --exclude CSL-0009:coolsolve`.

Verification of the physics: Table 6 (single phase), Tables 7 and 8 (bubble and dew points, compositions and densities), the pure-fluid limits
against CoolProp, the EES help examples, a round trip over all fifteen flashes on liquid, vapour, supercritical and two-phase states, and the
exactness of the derivatives against finite differences of the flash and of the language call (tests, finite differences allowed there).

## 6. Known limits

Transport properties (`Conductivity`, `Viscosity`, `SurfaceTension`) are not implemented (errors that say so). No two-phase state above the critical
locus (error naming the state). Saturation solves meet the strict criterion up to near-pure compositions, the water mole fraction being
evaluated without the rounding of `1 - x` (§4; a T-p pair near pure ammonia can still leave the frozen polish with the chemical potentials of water of
its two phases equal to a few 1e-11 of `R * 300 K` only). Within about 1e-8 of pure ammonia (wider at low temperature) the bubble point cannot be represented (its vapour would be purer than the
composition floor): liquid states there are an error, not a guess. A pure component (X = 0 or 1) given by `(T, x, v)` or `(P, x, h|s)` inside
its saturation dome comes out single-phase, and `(P, X = 0, H|S)` (pure water) fails (no density at the 200 K start of its search), both as
before. Near the critical locus a first call can still fail where a later one succeeds (T-p pairs within about 1e-5 in ln p of the pure-ammonia
pressure at 400 to 405 K). A first call costs up to some tens of ms (guesses), a warm call 0.1 to 2 ms.
