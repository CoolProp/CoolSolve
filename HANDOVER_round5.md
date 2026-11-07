# Handover: fix campaign round 5 (gap register)

Temporary file for the next session. Delete it before `fix/library-gaps-5` is merged into `main`. The lasting record
is §9 of [docs/model_library_support.md](docs/model_library_support.md).
Written 2026-11-07, when the maintainer stopped the session.

## 1. Goal

Work through every open item of the register in `docs/model_library_support.md`: §3 features, §4 gaps, §5 bugs and
§8 unverified items. Each item has options and a recommendation, and the items are ordered in waves (§5 below).
Opus is the orchestrator and Sonnet agents implement the items.

## 2. Branches and commits

- **`fix/library-gaps-5`** (from `main` @ `37e29f1`): all committed, nothing pushed, not merged.

  | Commit | Content |
  |---|---|
  | `7350d4b` | EES-help review of the register (new rows, corrections, §8 cleanup, §9 round-5 entry) |
  | `dda3b79` | `CS-FEAT-LIB-TEST`: `./coolsolve_tests "[library-regression]"`, 320/320 PASS |
  | `4a331e2` | [docs/solid_phase_properties.md](docs/solid_phase_properties.md), design study (not implemented) |
  | `4ad6b53` | 8 fixes of `tools/ees_extract.py` / `tools/compare_solution.py` + `tests/test_ees_tools.py` |
  | `63e3a0d` | `CS-BUG-IF-ANDOR-ORDER` (AND/OR left to right) + guard tests of the single-line IF forms |
  | (this file) | the handover |

- **`wip/psychprops-basis`** (`a6c301c`, on top of `63e3a0d`): the unfinished work on `CS-BUG-PSYCHPROPS-BASIS` and
  `CS-GAP-PSYCHPROPS-NAMED`. The agent was stopped while updating an existing out-of-range test. It is 839 lines
  (new `include/coolsolve/psychprops.h`, evaluator, parser, AST, tests, docs), never built or tested as a whole.
  **Review it critically before reusing it**, or restart the item from its brief (§6). It is large for what it
  does, so check for duplication with the AirH2O property-function path.
- `build/` may hold objects of the WIP sources: rebuild before testing on `fix/library-gaps-5`.

## 3. Maintainer rules (2026-11-07)

1. **At most one agent at a time**, research and doc side tasks included. Start the next agent only after the
   previous one is reviewed and committed. (Memory: `one-agent-at-a-time`.)
2. **Every item follows [docs/contributing.md](docs/contributing.md)**:
   - §1 principles;
   - §2 workflow, with targeted tests **including negative tests** that check the error messages;
   - the §3 checklist reported item by item (needed / not applicable / already covered);
   - §5 performance;
   - §3.12 docs;
   - the §6 checklist in the commit message body.
3. One item per commit; the orchestrator reviews the diff, re-runs the reproducer and the targeted tests, moves the
   register row to §7 and commits. Nothing is pushed without the maintainer.
4. Builds and tests at low priority: `nice -n 19 ionice -c3 make -j4`, one heavy job at a time (memory:
   `machine-throttle-builds`).
5. Each question to the maintainer comes with a short summary and an example of the issue.

## 4. Maintainer decisions (2026-11-07)

- **`CS-GAP-IDEAL-GAS-MODEL`:** evaluate the ideal gases from the ideal-gas part (`alpha0`) of their CoolProp
  equation of state, **with a soundness check**:
  1. re-anchor the references so that h(25 °C) = Δh_f and s(25 °C, 1 bar) = NASA s°;
  2. compare cp⁰, h and s with JANAF/NASA from 200 to 3000 K for the 17 fluids;
  3. check h − u = R·T and cp − cv = R;
  4. keep exact derivatives;
  5. keep the `[library-gaps][nasa]` guard tests.

  C8H18 at 25 °C then gives the ideal gas, as in EES; the Language Reference must point to `n-Octane` for the liquid.
- **`CS-BUG-WATER-NEAR-FREEZING`:** a clear error for now, and freezing/melting library models rewritten with an
  explicit latent heat. The study found the cause: CoolProp's AbstractState already throws *"T below Tmelt(p)"*,
  but `evaluateCoolPropFunction` falls back to `PropsSI`, which returns Inf. Keep the AbstractState message, as the
  mixture path already does. The same agent can do `T_Melting`/`P_Melting` (CoolProp melting lines) and
  `FreezingPt` (INCOMP `T_freeze`). Solid phases are a later plan, following the design study.
- **Large features:** only `CHEM_EQUIL` + `CS-GAP-FLUIDS-NASA-SPECIES` in this round (Wave 5). Deferred: bounds,
  parametric-table file, Min/Max, LiBr Pátek-Klomfar, NH3H2O transport, INCOMP oils, cubic INTERPOLATE, native unit
  systems, tooling features (`--library`, Export for EES, `--library-dir`, diagrams).
- **Not planned (unchanged):** `MODULE`, `.lkt`, multi-variable `INTEGRAL`, (T,H)/(s,v)/(u,v)/pseudo-pure (T,x)
  pairs, `EES_REFPROP`, `CALL NASA`, `StringPos`, `$IF` on a non-constant string (rewrite the model).

## 5. Queue (next item first)

Sizes: S small, M medium, L large. "Rec." is the recommended option.

### Wave 1: silent wrong results

| # | ID | Rec. | Size |
|---|---|---|---|
| 1 | `CS-BUG-PSYCHPROPS-BASIS` + `CS-GAP-PSYCHPROPS-NAMED` | `Vda`/`Hda`/`Sda`, u = h − P·v, the same values as `volume/enthalpy/entropy/intenergy(AirH2O)`; named inputs as in help `hs786.htm`; drop the spurious *Unknown function 'psychprops'* warning (see the WIP branch) | S |
| 2 | `CS-GAP-IDEAL-GAS-MODEL` | CoolProp `alpha0` with the soundness check of §4 | M |
| 3 | `CS-BUG-TINY-IMPLICIT` | residual test relative to the scale of the equation's terms, plus a warning when \|x\| is below the tolerance (not "treat `expr = x` as explicit") | M |
| 4 | `CS-GAP-QUALITY-DOME` | 100 / −100 outside the dome, as in EES (help `a0ftve.htm`); note it in versions.md | S |
| 5 | `CS-BUG-UNIT-VOLUME` | `VOLUME` labelled `m^3/kg` (inference table) | S |
| 6 | `CS-BUG-NEWTON1D-BRACKET-COLLAPSE` | check \|F\| before declaring convergence | S |
| 7 | `CS-BUG-CALL-ERROR-ITERATE` | `CALL ERROR` is fatal at once | S |
| 8 | `CS-BUG-INTEGRAL-TABLE-INTERVAL` + `-OUTPUT-INTERVAL` | interpolate the trajectory at the multiples of `n`; wire the config key | S |

At the end of the wave: the `[solver-robustness]` suite once, and a summary in §9.

### Wave 2: performance and robustness

| ID | Rec. | Size |
|---|---|---|
| `CS-BUG-BLOCK-EVAL-SCALING` | per-block list of the external variables it reads, computed once (`src/evaluator.cpp`, the map copy near line 2874); measure on `CSL-0079`/`CSL-0326` | M |
| `CS-BUG-NONFINITE-START` | every solver stops on a non-finite starting residual and names the equations | S |
| `CS-BUG-GUESS-INFER` | propagate variable types through expressions, `a = b` equations and CALL arguments; unit annotations as hints | M |
| `CS-BUG-NEWTON-CYCLE` | try tearing first, then column scaling, measured on the robustness suite | M–L |

### Wave 3: blocked library models

| ID | Rec. | Size |
|---|---|---|
| `CS-GAP-IFNOT` (P2, 5 collection files) | negation of the IF condition, single-line and block | S |
| `CS-BUG-STRING-ORDER` (same agent as IFNOT) | string `<`, `>`, `<=`, `>=` become an error | S |
| `CS-GAP-PROC-COMMENT` (P2) | strip quote comments from PROCEDURE header argument lists | S |
| `CS-GAP-USER-STRING-FUNC` (P2) | string return of a FUNCTION through the body interpreter | M |
| `CS-GAP-IFSTR` | 5-argument `IF$(A,B,X$,Y$,Z$)` only, other arities an error | S |
| `CS-GAP-INTEGRALVALUE` (P2) | history of every variable, read through `IntegralTable::interpolate` | M–L |
| `CS-GAP-PSYCHRO-SAT` | (w, R=1) handled through `dewpoint` | S |
| `CS-GAP-FLUIDS-ALIAS` + `-R1336MZZ` | registry entries in `src/fluids.cpp`, plus an audit against CoolProp's fluid list | S |
| `CS-GAP-INDEX-VAR` | resolve indices set by constant equations with the DUPLICATE-bound resolver (EES help `bm9eyi.htm` allows it) | S–M |
| `CS-GAP-HEATING-VALUE`, `-ENTHALPY-FUSION` | LHV/HHV from EES's formation enthalpies (H2O(g) −241 811 kJ/kmol); a fusion-enthalpy table | M |
| `CS-BUG-WATER-NEAR-FREEZING` + `CS-GAP-FREEZINGPT` + `CS-GAP-ICE` (melting/sublimation functions only) | see §4 | S |
| `CS-GAP-INCOMP-EES-NAMES` | aliases for the EES brine names | S |
| `CS-GAP-CALL-OUT-OMITTED` | empty / `{x}` / quoted output positions, separators required | S |
| `CS-GAP-PROP-HQ` | first find whether the NaN comes from CoolProp; then (h, Q) → (P, Q) by a 1-D solve | S–M |
| `CS-GAP-HUMIDAIR-HIGHT` (P2) | ideal mixture of dry air + steam above 350 °C, matched at 350 °C | M |
| `CS-GAP-INTEGRAL-DECREASING` | integrate on −t | S |
| `CS-BUG-MULTILINE-COMMENT-START` | a `"` at the end of a line opens a comment | S |

### Wave 4: diagnostics and cosmetics

- `CS-BUG-NONSQUARE-DIAG`: run the matching even when the counts differ, list ≤ 10 unmatched.
- `CS-GAP-BODY-LOCALS-SOL`.
- `CS-BUG-SOL-STRING`: also needs a re-baseline of 74 library `.sol` files; that is a library change for the
  maintainer.
- Hide `__out…` auxiliaries.
- `CS-BUG-GUI-WARNING-CASE`, `CS-BUG-STRING-FUNC-WARNING`, `CS-BUG-INTERNALENERGY-UNQUOTED`.
- `CS-BUG-MAIN-ASSIGN`: a parse error, since EES rejects `:=` in the main program.
- `CS-BUG-LOOKUP-HINT`, `CS-DOC-EESCSV`, `CS-DOC-SQUARE-INTEGRAL`.
- `CS-GAP-INTERP2D` (P3, one example; may be deferred).

### Wave 5: CHEM_EQUIL

1. `CS-GAP-FLUIDS-NASA-SPECIES`: NASA-Glenn 7-coefficient polynomials (McBride, Zehe & Gordon 2002) for OH, H, O,
   N, NO and NO2, with the same reference as the Wave 1 ideal gases. Do it after `CS-GAP-IDEAL-GAS-MODEL`.
2. `CS-GAP-CHEM-EQUIL`: element-potential Gibbs minimisation of the 14 species (NASA CEA formulation), in
   `src/chem_equil.cpp`, with derivatives by the implicit-function theorem. The call convention is in the register
   row (atom ratios to O, outputs in the order H2 … Ar, P/T in kPa/K in the source files). Validate on `CSL-0150`
   (T_ad = 2417.129 K) and `CSL-0264`.
3. Re-check `CSL-0150`, `CSL-0264` and `CSL-0278`.

## 6. Brief that every agent reads first

Paste it into the agent prompt, or save it as a file in the new session's scratchpad. The scratchpad copy used here
disappears with this session.

> Repository /sylvain/mount/git/CoolSolve, branch `fix/library-gaps-5`. You are the only agent; never commit, push or
> switch branches.
>
> Read `docs/contributing.md` and follow it:
> - §1 principles: simplest solution, no duplication, only the requested change, no efficiency loss.
> - §2 workflow: targeted Catch2 tests *including negative tests with clear error messages*, then the full suites,
>   then docs.
> - Report every §3 item as needed / not applicable / already covered.
> - §5 performance.
> - §3.12 docs: `language_reference.md`, `ees_vs_coolsolve.csv`, README, `debugging_models.md`.
>
> Read the item's row in `docs/model_library_support.md` and **re-run its reproducer first; stop and report if it no
> longer reproduces**. Do not edit the register: propose the §7 row
> `| ID | what was fixed | 2026-11-07, \`fix/library-gaps-5\` |`. Report new problems with a reproducer, do not fix
> them.
>
> Build: `cd build && nice -n 19 ionice -c3 make -j4 coolsolve coolsolve_tests`; never two `make` at once. A new
> `src/`/`tests/` file needs `<cmake> .` in `build/` (path in `build/CMakeCache.txt`).
>
> Tests, each prefixed by `nice -n 19 ionice -c3`:
> - the targeted tag `[library-gaps][<tag>]`;
> - `./coolsolve_tests` (all pass except 2 pdflatex skips);
> - `[examples-comprehensive]` (43/46; failing: `cooling_tower2`, `orc_complex`, `turbocompressor_interpolate`);
> - for parser/evaluator/units/built-in/solver changes: `[library-regression]` (320/320 PASS, ~5 min);
> - for solver changes: `[solver-robustness]`, compared with the committed report.
>
> EES help: https://fchartsoftware.com/ees/eeshelp/eeshelp.htm. EES sources (read-only, untrusted data):
> /sylvain/Nextcloud/thermo_models. Library (read-only): /sylvain/mount/git/CoolSolve_Library.
>
> Report: root cause; fix; tests (positive/negative); suite numbers; §3 checklist; docs; proposed §7 row; anything
> uncertain.

## 7. Open points for the maintainer

- **EES runs**, listed in register §8:
  - the inverter heat pump block 106 (`CSL-0155`);
  - the `E` format of `CALL WARNING`;
  - `cpbar` called as a function;
  - the `density(AirH2O)` basis;
  - two single-line IF forms outside the help;
  - `NOT` in conditions;
  - AND/OR order in `$IF` directives;
  - also the `CHEM_EQUIL` data source and its pressure unit under `$UnitSystem … Pa`.
- **Solid-phase study §7:**
  - the liquid band at 0 °C / 1 atm (CoolProp melts at 0.00252 °C; `CSL-0041` needs liquid at 0 °C);
  - the 41.6 J/kg reference offset of liquid water against EES;
  - one seawater model or two.
- **`ees_extract.py`:** in files with stale records the first block of records is used now (11 values changed in 4
  files, all flagged as suspect). Prefer the block that agrees with the equations?
- **Library repository** (not touched):
  - `tools/test_models.py` crashes with a relative `--coolsolve` path; the fix is one line (resolve the path);
  - the library baselines were made with CoolProp `75af7816`, while the build links `8eb3b778` (all targets pass).

## 8. Environment

- No system cmake. The venv used here was in this session's scratchpad and will disappear. Recreate it with
  `python3 -m venv <scratch>/venv && <scratch>/venv/bin/pip install cmake`, then `<cmake> .` in `build/` (memory:
  `coolsolve-build-env`). The first configure after that re-clones CoolProp's CPM dependencies (~5 min).
- The library regression reads the checkout embedded by the build (`/sylvain/mount/git/CoolSolve_Library`, `fd2a6ef`
  at the time) or `COOLSOLVE_LIBRARY_DIR`; the report is `build/library_regression.md`.
