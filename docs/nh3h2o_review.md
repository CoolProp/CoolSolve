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
| `include/coolsolve/nh3h2o.h`, `src/nh3h2o.cpp` | 206, 363 | the Helmholtz formulation (Eqs. 1 to 8 and Tables 1 to 3, A1 of the guideline; IAPWS-95 for water) written once over a scalar type; the properties of a state |
| `src/nh3h2o_flash.cpp` | 1290 | the saturation solver, the 15 flashes, derivatives by the implicit function theorem. **Contains its own small Newton (`newtonSolve`), Brent and linear solver**, private to the file (anonymous namespace) |
| `src/evaluator_nh3h2o.cpp` | 306 | the glue: units of the model, designators, the trial-point limits, assembly of the `ADValue` gradient, the two entry points |
| `tests/test_nh3h2o.cpp` | ~850 | 24 test cases (tags `[nh3h2o]`, `[formulation]`, `[equilibrium]`, `[flash]`, `[evaluator]`) |

### B. Shared files modified (the part to review): 9 files, about 50 lines

| File | Change | Why other models are not affected |
|---|---|---|
| `src/evaluator.cpp` (+10/-2) | (1) in `evaluateProcedureCall`: `if (lowerName == "nh3h2o") { evaluateAmmoniaWaterCall(call); return; }` placed after the user-procedure lookup and next to the Colebrook hook; (2) in `evaluateCoolPropFunction`, right after the named inputs are evaluated: `if (fluid->getType() == FluidType::AmmoniaWater) return evaluateAmmoniaWaterFunction(...)` | (1) only reached for a `CALL` named `NH3H2O` when the model defines no such procedure; (2) only for the new fluid type, which no other fluid has. Both before any CoolProp code of the call, so nothing after them changes |
| `include/coolsolve/evaluator.h` (+6) | declares the two member functions | declarations only |
| `include/coolsolve/fluids.h` (+13) | new enum value `FluidType::AmmoniaWater` (appended before `Unknown`) and the class `AmmoniaWaterFluid` | an added enumerator; the existing ones keep their values. Every `switch` over `FluidType` was checked (`server.cpp`: two, edited) |
| `src/fluids.cpp` (+7/-4) | `NH3H2O` registered as `AmmoniaWaterFluid` instead of `UnsupportedFluid` | the one intended behaviour change: `NH3H2O` was an error, it now evaluates |
| `src/parser.cpp` (+1) | `"nh3h2o"` added to `knownBuiltinFunctions()` (as `"colebrook"`) | silences the "Unknown function" warning of `CALL NH3H2O`; the effect on the parser is that a user-defined `FUNCTION nh3h2o` would be taken for a built-in name in `isPositionalOnlyBuiltin` (named arguments refused); no known model |
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
5. **Thread safety**: the saturation solver keeps one warm-start state per kind in `thread_local` storage (`g_satCache`); no locking.
6. **Reference state**: IAPWS (u = s = 0 for the pure saturated liquids at their triple points), not EES's: documented as a difference,
   not corrected (differences at fixed composition agree).

## 4. The solver inside the module (the private Newton)

`newtonSolve` (`src/nh3h2o_flash.cpp`) is the only solver written for this work; it is not shared with the rest of CoolSolve.
**Convergence is strict**: a scaled residual below 1e-12 (pressure of a liquid on the scale `rho R T`, the scale of its round-off, measured 3e-13
of it; chemical potentials on `R * 300 K`). An earlier version had a loose fallback (step below 1e-9 and residual below 1e-7); it was
**removed** (`db0246a`) after it was found to accept under-converged states in a third of the iterations. The one remaining exit
for an iteration that has stopped lowering the residual requires every residual to be within two units of a **measured** noise floor
(`noiseFloor`: the resolution of the representation of the variables, `sum |J_ij| ulp(v_j)`, which matters for a mole fraction within 1e-5 of 1
because `1 - x` is then known to `1e-16 / (1 - x)` only, plus the jitter of the evaluation measured at the point). The test `Saturation: the
equilibrium equations hold to 1e-11 of their natural scales` recomputes the equilibrium residuals from the Helmholtz energy with none of the
solver's code. 
Two further pieces matter for compositions near pure ammonia (the refrigerant of a chiller, 0.9995 mass fraction):

- **Frozen-variable polish** (`frozenPolish`, inside `newtonSolve`). When a mole fraction within about 1e-5 of 1 cannot move by less than one ulp, the
  full Newton step asks for a change the double cannot represent and the other variables overshoot. At a stall the quantised variable and the
  row it dominates are frozen, the other rows are solved **strictly** (scaled residual below 1e-12) for the other variables, and the frozen variable
  is moved by Newton on the dropped row, whose sensitivity is the Schur complement of the Jacobian. The dropped row is accepted within two units of
  its quantisation floor (`|S| ulp`) plus the measured jitter. So at a near-pure state all rows but one meet the strict tolerance and the last one is
  as small as the representation allows.
- **Initial guesses** (`guessBubbleT`, `guessDewT`, `guessTP`, `guessAtPressure`): the liquid is evaluated from its fugacities against a vapour
  that is first an ideal gas, then corrected by the fugacity coefficients of the real vapour (they matter within a few percent of the pressure of
  pure ammonia, where the ideal-vapour guess has no root). `satSolve` tries the real-vapour guess first and the ideal-vapour one if Newton
  does not converge from it (at 16 MPa the correction over-corrects). Guesses never affect the accuracy: Newton polishes every result.

Reviewing `newtonSolve` (with `noiseFloor` and `frozenPolish`), `satNewton`, `finishTwoPhase` and `finishSingle` covers the numerics.

Comparison tolerances of the tests (against the guideline tables, EES, CoolProp) are separate from the solver tolerance and are stated in the tests:
Table 6 to 1e-7 relative (the table has 9 to 10 digits and the coefficients of the guideline 7 to 8), pure water to 1/60 000 as the guideline says,
the EES points to 0.2 % (volume) to 2 % (composition).

## 5. Verification and baselines

| Check | Before (`main` @ `a77efaa`) | Branch |
|---|---|---|
| Unit suite (`cd build && ./coolsolve_tests`) | 1182 test cases | 1209 test cases (+27), all pass, 19 843 assertions |
| Examples (`[examples-comprehensive]`) | 43/46 solve (`cooling_tower2`, `orc_complex`, `turbocompressor_interpolate` do not) | 43/46, the same three |
| Library regression (`tools/test_models.py`, scratch copies) | 249 targets, 5 expected deviations (`CSL-0008` ×3, `CSL-0146`, `CSL-0160:coolsolve`) | 249 targets, **the same 5** deviations, none new (run with the binary of `6ecd20b`, the evaluator integration; the later commits change only files of the module that no other model reaches) |

To reproduce: `cmake -S . -B build && cmake --build build -j && cd build && ./coolsolve_tests` and, in the library,
`python3 tools/test_models.py --coolsolve <build>/coolsolve --exclude CSL-0009:coolsolve`.

Verification of the physics: Table 6 (single phase), Tables 7 and 8 (bubble and dew points, compositions and densities), the pure-fluid limits
against CoolProp, the EES help examples, a round trip over all fifteen flashes on liquid, vapour, supercritical and two-phase states, and the
exactness of the derivatives against finite differences of the flash and of the language call (tests, finite differences allowed there).

## 6. Known limits

Transport properties (`Conductivity`, `Viscosity`, `SurfaceTension`) are not implemented (errors that say so). No two-phase state above the critical
locus (error naming the state). Saturation solves with a composition within about 1e-5 of pure are accurate to the resolution of a double
(see §4). A first call costs up to some tens of ms (guesses), a warm call 0.1 to 2 ms.
