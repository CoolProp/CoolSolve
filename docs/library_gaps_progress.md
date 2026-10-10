# Library-gap fix campaign — progress tracker

Working file of the campaign that fixes the CoolSolve gaps and bugs found while
building the [CoolSolve Library](https://github.com/IntSusEnergySystems/CoolSolve_Library)
(register: [model_library_support.md](model_library_support.md) §4–§5). It is
kept up to date after every step so that the work can be **restarted after an
interruption**: read §1–§3, then continue at the first item of §4 that is not ✅.

Started 2026-10-10. Supervisor: Claude Code (Opus) session; implementation by
Claude Sonnet agents.

---

## 1. Decisions (maintainer, 2026-10-10)

- **Git:** all CoolSolve work on branch `fix/library-gaps-2` (from `main` at
  `9039b7a`), **one commit per fix**, made by the supervisor after review; no
  push. The library repository (`/home/sylvain/svn/CoolSolve_Library`) is
  **never committed** (library rule D9): changes stay in its working tree.
- **Agents:** Sonnet agents, **one at a time, strictly sequential** (no two
  agents running together, even fix + re-check).
- **Order:** by number of library models unblocked, silent P1 bugs early (§4).
  `CS-GAP-INCLUDE` / `CS-FEAT-IMPORT` is a later, separate feature.
- **Process per item** (maintainer):
  1. follow [contributing.md](contributing.md);
  2. write the test of the problem first (`tests/test_library_gaps.cpp`, tag
     `[library-gaps][<tag>]`, names `"<ID>: …"`), check that it fails;
  3. implement until the test passes;
  4. run `./coolsolve_tests` and `./coolsolve_tests "[examples-comprehensive]"`
     (**not** `[solver-robustness]`, too slow);
  5. re-run the library models blocked by the item, check them, clean the model
     folder and document the improvement in the library;
  6. move to the next item.

## 2. Baselines (before any fix, `main` @ `9039b7a`, Release build)

- `./coolsolve_tests`: all tests passed (3252 assertions in 379 test cases).
- `./coolsolve_tests "[examples-comprehensive]"`: passes; 46 files, 39 solved,
  39/39 expected values correct. Not solving at baseline: `building_rc_network`
  (analysis), `cooling_tower2`, `internal_combustion_engine_cpbar`,
  `orc_complex` (parse), `simple_centrifugal_compressor`,
  `turbocompressor_interpolate`, `water_libr`; `advanced_features` has
  `1 FAIL` in Eq Check.
- Library regression (`python3 tools/test_models.py --coolsolve <baseline binary>`,
  148 models: 98 verified / 7 runs / 43 blocked): 151 targets (105 main files,
  46 variants), **7 failures already at baseline** (not caused by the campaign;
  the library's `.sol` files were produced with another build, probably another
  CoolProp version): `CSL-0114` (a_ev, deltap_ev), `CSL-0038` and `CSL-0039`
  (`balance` ≈ 1e-9, absolute-tolerance noise), `CSL-0019` (R245fa `mu`, `k`
  +15–20 %: transport-property difference), `CSL-0153:coolsolve` (97 diffs),
  `CSL-0158:coolsolve` (≈ 1e-5), `CSL-0124:coolsolve` (deltap, alpha_tp ≈ 2–6 %).
  Variants without `.sol` (not tested): `CSL-0155:coolsolve`,
  `CSL-0162:coolsolve`. Run time 3 min 38 s.
- **CoolProp is pinned** for the whole campaign to the baseline commit
  `5b9c32ac` (2026-08-22): `cmake -DCOOLSOLVE_COOLPROP_TAG=5b9c32ac ..` in
  `build/` (local cache setting only). The default tag `master` made a
  `cmake .` re-run pull `75af7816` (2026-10-09), which changed the failing
  example `internal_combustion_engine_cpbar` (block 39 residual) — a CoolProp
  drift, not a regression. Unpin (`-DCOOLSOLVE_COOLPROP_TAG=master`) at the end.
- The `[library-gaps]` tests were not compiled by the baseline build (the build
  tree predated `tests/test_library_gaps.cpp`; a `cmake .` fixed it): with them
  the unit suite has 424 test cases after fix #1 (379 + 45).

## 3. Procedures

### 3.1 Fix task (one Sonnet agent)

The agent prompt points to the *fix protocol* below (copy it to the session
scratchpad as `fix_protocol.md` if needed) and gives the item ID, the expected
EES behaviour and the native library files that use the feature (read-only).

> **Fix protocol.** Read `docs/contributing.md` (§1–§3), the register row(s) of
> the item in `docs/model_library_support.md` and §7 (format of closed
> entries), and `tests/test_library_gaps.cpp` (style). Then: (1) reproduce the
> register reproducer with `build/coolsolve` in a scratch folder; (2) write the
> tests first — register reproducer with the EES value, realistic usage from
> the blocked library models, derivatives if numeric (analytic AD only),
> negative tests with clear error messages — and check they fail; (3) implement
> the simplest fix reusing existing code, zero overhead when unused; (4) build
> and iterate until `[library-gaps]` passes; (5) run `./coolsolve_tests` and
> `./coolsolve_tests "[examples-comprehensive]"` and compare with §2 (no example
> may regress); (6) docs: `language_reference.md`, `ees_vs_coolsolve.csv`,
> `gui/src/languages/ees.ts` for new built-ins, move the register row from
> §4/§5 to §7 (*Fixed in* = `post-v0.3.0, branch fix/library-gaps-2`; the
> supervisor adds the commit hash), README only if a feature bullet changes;
> (7) do not commit; report root cause, design, files, tests (fail before /
> pass after), suite results, docs, open points.

### 3.2 Supervisor review and commit

Read the diff, re-run `./coolsolve_tests "[library-gaps]"`, `./coolsolve_tests`
and `[examples-comprehensive]`, commit with message
`Fix <ID>: <summary>` (+ co-author line), then record the hash in §4 here and in
the register §7 *Fixed in* column (in the next commit).

### 3.3 Library re-check task (one Sonnet agent, after the commit)

For every model that lists the ID in `model.json` `missing_features`
(`python3 - <<EOF` scan of `models/**/model.json`):

- remove the ID from `missing_features`;
- if no gap is left: solve the native file with the new build in a scratch copy
  of the folder, compare with the EES reference (the variant's `.sol`, and
  `tools/compare_solution.py` against the EES source when available); if it
  agrees: status `verified` (or `runs`), native `.sol` written, `verification`
  block updated (`coolsolve_version` `0.3.0@<commit>`, date, max_rel_diff,
  notes), README updated (status row, *How to run*, *Verification*,
  *Conversion log*, *Limitations and CoolSolve gaps*), the now useless
  `_coolsolve` variant files deleted (clean folder), `history` entry
  `recheck`;
- if other gaps remain: README and `model.json` updated (the ID is closed),
  variant kept;
- then `python3 tools/build_index.py`, `python3 tools/test_models.py <ids>`,
  one line in the library `roadmap.md` §8 progress log. No library commit.

## 4. Queue and status

Legend: ☐ to do · 🔄 in progress · ✅ done (commit) · ⏸ deferred.
*Models* = library models listing the ID at the start of the campaign
(models freed by the item alone in bold).

| # | Item(s) | Models | Fix | Commit | Library re-check |
|---|---|---|---|---|---|
| 1 | `CS-GAP-IF5` | 11: 0009, 0041, **0073**, 0074, 0075, 0080, 0081, 0153, **0155**, 0162, 0166 | ✅ | `59b2862` | ✅ |
| 1b | `CS-BUG-ARG-DROP` (new, P1 silent: unparsable positional argument dropped) + close `CS-GAP-UNIT-SUBEXPR` / `CS-GAP-UNIT-NAMEDARG` (already parse since `59b2862`; tests needed) + named-argument misparse `IF(a=5,1,2)` → *Unknown fluid ''* (§8 C-147) | – | ✅ | `6c36a30` | – (no model listed the IDs; see log) |
| 1c | `CS-BUG-BODY-LINE-DROP` (unparsable body line dropped) + `CS-BUG-TRAILING-TEXT` (text after `)`/`]` ignored) — both new, P1 silent, found in #1b | – | ✅ | `8238a94` | `CSL-0103` comment fix (see log) |
| 2 (scope) | `CS-GAP-IF-DIRECTIVE` + `CS-BUG-STRING-CASE` (one shared case-insensitive string equality) | see #2 | ✅ | `9423934` | ✅ |
| 2 | `CS-GAP-IF-DIRECTIVE` (`$if/$ifnot/$else/$endif`) | 10: 0012, 0040, 0081, 0116, 0119, 0150, 0151, 0161, 0162, 0166 | ☐ | | ☐ |
| 3 | INTEGRAL family: `CS-GAP-INTEGRAL-LIMITS`, `CS-BUG-INTEGRAL-FACTOR`, `CS-BUG-INTEGRAL-TABLE-SEP`, `CS-BUG-INTEGRAL-MAXSTEPS`, `CS-BUG-INTEGRAL-TABLE-CASE` (then `CS-GAP-INTEGRAL-MULTIVAR`) | ✅ | `4cd0ca5` | ✅ |
| 4 | `CS-GAP-LOOKUP-PROC` + `CS-GAP-INTEGRAL-LOOKUP` (same root cause: table store not wired; frees CSL-0042, example `building_rc_network`) | ✅ | `d5d6b37` | ✅ |
| 5 | `CS-GAP-CALL-EXPR-OUT` + `CS-BUG-STRING-CALL-OUT` (CSL-0166) | ✅ | `57b22d3` | 🔄 |
| 6 | Identifier characters: `CS-GAP-NAME-SYMBOL` (`C%`), `CS-GAP-NAME-PIPE` (`K\|star`), `CS-GAP-ACCENTED-IDENT`, `CS-BUG-HASH-IDENT` (`Nu#_f`, moved from item 9) | 0040, 0074, 0080, 0116, 0150, 0151, 0152, 0158, 0163 | ☐ | | ☐ |
| 7 | BrineProp family: `CS-GAP-ELSEIF-CHAIN`, `CS-GAP-UPPERCASE`, `CS-GAP-STRING-ARRAY` | 0072, 0078, 0079 | ☐ | | ☐ |
| 8 | `CS-GAP-INTERP-EES` + `CS-GAP-INTERP-EXTRAP` (silent) | 0011, 0035, 0149 | ☐ | | ☐ |
| 9 | Silent P1 bugs: `CS-BUG-VAR-ANNOTATION` (new, silent in bodies; needs EES evidence of `a [kW]`), `CS-BUG-MAIN-SEMICOLON` + `CS-BUG-DUPLICATE-HEADER-COMMENT` (new, P3, parser hygiene, same family), `CS-BUG-PI-FUNCTION`, `CS-BUG-CONVERT-UNIT`, `CS-BUG-HASH-IDENT`, `CS-BUG-COMMON-PROC`, `CS-BUG-IMPLICIT-PROC`, `CS-BUG-LOOKUP-COL-ARG`, `CS-BUG-LOOKUP-COLNAME`, `CS-BUG-DUPLICATE-VAR-BOUND`, `CS-BUG-DUPLICATE-INDEX-EXPR`, `CS-GAP-GOTO`, `CS-BUG-HUMIDAIR-PROPS`, `CS-GAP-AIR-ENTHALPY-REF`, `CS-BUG-STRING-CASE` | 0075, 0080, 0081, 0118, 0123, 0124, 0127, 0129, 0131, 0157, 0158, 0159, 0161 | ☐ | | ☐ |
| 10 | Cosmetic: one convention for internal variables (`__…` names hidden from `.sol`/JSON/GUI/debug, as the INTEGRAL auxiliaries already are; `__out…` of #5 visible). Single-model items (by register priority): `CS-GAP-PSYCHRO-SAT`, `CS-GAP-FLUIDS-C8H18`, `CS-GAP-FLUIDS-ALIAS`, `CS-GAP-INTEGRAL-LOOKUP`, `CS-GAP-REPEAT-UNTIL-BARE`, `CS-BUG-LOOKUP-WRITE`, `CS-GAP-CALL-WARNING`, `CS-GAP-NDAY`, `CS-GAP-END-PROCEDURE`, `CS-GAP-CONVERT-UNQUOTED`, `CS-GAP-LOOKUPROW`, `CS-BUG-LOOKUP-STRING`, `CS-GAP-IFSTR`, `CS-GAP-ISIDEALGAS`, `CS-GAP-SUM-INDEXED`, `CS-GAP-GE-ARROW`, … | – | ☐ | | ☐ |
| 11 | `CS-GAP-INCLUDE` / `CS-FEAT-IMPORT` (separate feature) | 0011, 0074, 0116, 0149, 0162, 0163 | ⏸ | | |

## 5. Log

| Date | Step | Result |
|---|---|---|
| 2026-10-10 | Set-up | Branch `fix/library-gaps-2` created; baselines §2 recorded; frozen baseline binary in the session scratchpad (`coolsolve_baseline`). |
| 2026-10-10 | #1 `CS-GAP-IF5` | Fixed (`59b2862`): 3- and 5-argument `if` in `evaluateStandardFunction`; also unit-annotated literal arguments (`3600[s]`) no longer dropped, and `INTEGRAL` base expressions evaluate built-in functions instead of 0. 13 tests. Suites: 424/424 unit, examples identical to baseline. Native CSL-0073 solves, identical to its variant `.sol`; CSL-0155 parses but fails on block 106 (151 vars, NaN at start, not IF-related). New §5 row `CS-BUG-ARG-DROP`. |
| 2026-10-10 | #1 re-check | `CSL-0073` blocked → **verified** (native identical to the former variant, 34/34; EES max 2.71e-2 humid-air, unchanged), variant removed, native `.sol` added. `CSL-0155` blocked → **failing** (`missing_features` []): all IFs evaluate, block 106 (151 vars) does not converge even from CoolProp-consistent guesses; `.initials` rebuilt, non-runnable variant removed. 9 other models: IF5 removed from `missing_features`, README/history updated, none stops on IF5 any more. `build_index` 0 errors; `test_models` on the 11 ids OK except baseline `CSL-0153:coolsolve`. Library roadmap §1/§8 updated. Register: new `CS-GAP-REF-STATE` (EES R410A ASHRAE reference vs CoolProp IIR, offsets 142 244.6 J/kg / 773.6 J/kg·K confirmed), §8 C-147 resolved. |
| 2026-10-10 | #5 | Fixed (`57b22d3`): expression CALL outputs lowered to an auxiliary output + equation (solved-for form for one variable under + - * /), error inside bodies; strings set by CALLs published after each block (zero cost otherwise), string-producing CALLs first in a block. 16 tests; 552/552 unit; examples identical; library regression 151 / same 7. Natives: CSL-0074/0163 CALL lines parse (next: `K\|star`), CSL-0162 parses and is square (next: INCLUDE BRINEPROP, then NEWTON-CYCLE), CSL-0166 with `$if` resolved by hand = variant to 8.2e-7, CSL-0072 with BRINEPROP rewritten = variant to 2.1e-10. New rows `CS-GAP-CALL-OUT-OMITTED`, `CS-BUG-STRING-ALIAS`, `CS-BUG-DUPLICATE-CALL-BLOCK`. Auxiliary `__out…` variables are visible in outputs (item 3 hides its `__integral_k`): queued a cosmetic unification in item 10. |
| 2026-10-10 | #4 re-check | **Freed: CSL-0042 blocked → verified** (native = variant on 56 vars, max 5e-12; independent RK4 9.8e-5; variant + `clamp01` tables removed) and **CSL-0123 blocked → runs** (native = variant on 13 vars; no independent reference: the stored EES run used a lost table version, as the README says). LOOKUP-PROC removed from CSL-0072/0078/0079/0116/0158/0160 (BRINEPROP lookups verified in scratch: ≤ 5e-10). Register: `CS-BUG-LOOKUP-COLNAME` closed (part a not reproducible even at baseline — checked by the supervisor; part b fixed in `d5d6b37`), removed from CSL-0158; §8 C-83 (CSL-0078 "redundant equations") resolved — see §6. Library roadmap §1: 103 verified / 8 runs / 36 blocked / 1 failing. |
| 2026-10-10 | #4 | Fixed (`d5d6b37`): table store wired into PROCEDURE bodies (`makeBodyEvaluator`), IntegralSolver, debug report, server parametric fast path; + `CS-BUG-LOOKUP-COL-ARG` (string column argument read as a number — would have become silent in bodies; CSL-0123 native now equals its variant). 15 tests; 536/536 unit; example `building_rc_network` now solves (T_in 30.6825, = CSL-0042 variant; EXPECTED_SOLUTIONS added; examples 40/46, baseline copy updated); integral API 14/14; library regression 151 / same 7. Library natives: CSL-0042 now solves; CSL-0079 and CSL-0072 lookups verified in scratch copies (≤ 5e-10 vs variants) but still blocked by ELSEIF-CHAIN/UPPERCASE/STRING-ARRAY/CALL-EXPR-OUT; CSL-0078 non-square; CSL-0116/0158/0160 parse errors first. New row `CS-GAP-INTEGRAL-BASE-CALL`. |
| 2026-10-10 | #3 re-check | **Freed: CSL-0047, CSL-0010, CSL-0009 blocked → verified** (natives identical to their variants; vs EES: 0047 max 9e-8, 0010 max 4.5e-6 = EES trapezoid vs CoolSolve quadrature, 0009 ≤ 6.1e-5 on curves decoded from the EES plots); variants (and 0009's `coolsolve.conf`) removed, native `.sol` added; CSL-0010 keeps `CS-GAP-OPTIM` for information. CSL-0041 now blocked only by `CS-BUG-WATER-NEAR-FREEZING` (solves with the 3 Water calls replaced). CSL-0059: MULTIVAR + DECREASING. CSL-0081: GOTO + INTEGRALVALUE. Header catch-up cleanup done (IF5/IF-DIRECTIVE mentions). Library roadmap §1: 102 verified / 7 runs / 38 blocked / 1 failing. Register §8: `$IntegralTable` interval rounding note (unverified). |
| 2026-10-10 | #3 | Fixed (`4cd0ca5`): FACTOR (fast path for `y = base + INTEGRAL`, auxiliary states `__integral_k` otherwise, expression integrands get their own equation, internal names hidden), symbolic LIMITS (solved at the initial point), MAXSTEPS (explicit step honoured, early stop = error), TABLE-SEP/CASE (`, ;` separators, case-insensitive, I001 warning). 35 tests; 521/521 unit; examples identical (ice_storage_tank `-integral.csv` T_w column now filled); `test_integral_api.py --auto` 14/14; library regression 150 / same 7; CSL-0009:coolsolve 1m20 before/after. Natives: CSL-0047, CSL-0010, CSL-0009 solve, identical to their variants; CSL-0041 stops on WATER-NEAR-FREEZING; CSL-0059 on MULTIVAR (+ new DECREASING); CSL-0081 on GOTO. New rows `CS-GAP-INTEGRAL-DECREASING`, `CS-GAP-INTEGRAL-DUPLICATE`, `CS-BUG-INTEGRAL-OUTPUT-INTERVAL`. |
| 2026-10-10 | #3 (earlier) | Fix agent dispatched; interrupted by an API usage limit mid-implementation (partial diff ≈1800 lines left in the working tree, uncommitted); resumed with its context after the reset. |
| 2026-10-10 | #2 re-check | 11 models (+ CSL-0074, CSL-0153 for CALL-EXPR-OUT): none freed, all stay blocked; closed IDs removed, `CS-GAP-IF-NO-CONSTANT` added to CSL-0040/0166, `CS-GAP-CALL-EXPR-OUT` added to CSL-0074/0153; READMEs say where a variant's `$IF` workaround is no longer needed. Evidence: CSL-0075 native solves identically to the variant once its 2 humid-air rewrites are applied; CSL-0151 solves once `é`→`e`; CSL-0161 once `Nu#_f` renamed; CSL-0119 identical with `nDay_`→338. Gaps reached only on alternate branches kept (CSL-0150 NAME-SYMBOL on `'Oil'`, CSL-0161 CALL-EXPR-OUT on `'EG'`). build_index 0 errors; test_models OK (only baseline CSL-0153:coolsolve). **Pending cleanup** (give to the next re-check): native-file header comments still cite closed gaps — `CS-GAP-IF5` in 6 natives (boiler_modulating, glycol_runaround, iso5167, dhw_tank, ice_storage, orc_whr_refprop_cost), `CS-GAP-IF-DIRECTIVE` in 8 natives (thermal_comfort, boiler_modulating, glycol_runaround, plate_hx_thermal_resistances, adiabatic_flame, orc_biomass, engine_derating, solar_geometry); protocol B now requires the header update. |
| 2026-10-10 | #2 | Fixed (`9423934`): `$IF/$IFNOT/$ELSE/$ENDIF` resolved while parsing (removed branches skipped unparsed, main program and bodies; EES help semantics: string constant set before the directive, case-insensitive, AND/OR, EES flags; errors for non-constant strings, unknown conditions, unbalanced directives; `report.md` table); `eesStringEqual` (case-insensitive incl. Latin-1 accents) for IF statements and `$IF`. New files `src/conditional_directives.cpp`, `src/ees_string.cpp` (+ headers). 24 tests; 486/486 unit (one old test asserted case-sensitivity); examples identical; library regression 151 / same 7. No model freed by #2 alone: CSL-0012 → `LOOKUP$ROW`, 0119 → `nDay_`, 0151 → accented names, 0161 → `#` names, 0081 → GOTO, 0116 → `=>`, 0150 → not square (CHEM-EQUIL/NASA), 0162 → CALL-EXPR-OUT, 0040 → Diagram-window strings (new `CS-GAP-IF-NO-CONSTANT`) + `C%`, 0166 → `regime$` from a CALL (`CS-GAP-IF-NO-CONSTANT`), 0075 → now AirH2O branch, SingularJacobian block 16 (ISIDEALGAS/HUMIDAIR-PROPS open). PARAMETRICTABLE stays false in GUI parametric studies (documented). |
| 2026-10-10 | #1c | Fixed (`8238a94`): unparsable body lines reported ("Could not parse line"), `END // comment` closes a body, text after `)`/`]`/string is an error, unit annotation after a call/element accepted explicitly (also inside expressions). 15 tests; 462/462 unit; examples identical. AST comparison: 6 files gain errors on previously dropped lines — `orc_solar_complex` example (`Nu#_` names), blocked natives CSL-0081 (GOTO), CSL-0153/0160 (`END procedure`/`END name`), CSL-0158 (`°` names) — and **CSL-0103 (verified)**: six `!` comment lines in `FUNCTION dP_Kern` (not an EES comment form, silently dropped before) → fixed in the library as `//` comments, regression OK, history entry added. Register: GOTO and HASH-IDENT now loud (P1→P3); new rows `CS-BUG-VAR-ANNOTATION` (silent in bodies), `CS-BUG-MAIN-SEMICOLON`, `CS-BUG-DUPLICATE-HEADER-COMMENT` → item 9. |
| 2026-10-10 | #1b | Fixed (`6c36a30`): unparsable positional arguments, indices, CALL inputs/outputs are parse errors (also in bodies); positional-only built-ins refuse `name=value` (`IF(a=5,…)`); `(…) [unit]` accepted; `CS-GAP-UNIT-SUBEXPR`/`-NAMEDARG` closed with tests. 23 tests; 447/447 unit; examples identical; library regression 151 targets, same 7 baseline failures. AST comparison over 267 files: 14 differ, all blocked natives/failing examples now erroring on known gaps (`CALL-EXPR-OUT`, `NAME-SYMBOL`, `ACCENTED-IDENT`, `ARRAY-ARG`). **Library follow-up for item 5:** `CSL-0153` and `CSL-0074` use expression CALL outputs but do not list `CS-GAP-CALL-EXPR-OUT`. New rows `CS-BUG-BODY-LINE-DROP`, `CS-BUG-TRAILING-TEXT` → item 1c. |

## 6. Follow-ups for the maintainer (library decisions, not CoolSolve fixes)

- **CSL-0078** `brine_to_water_heat_pump_refsim`: the *native* file assigns `P_ev`, `P_cd`, `M_dot_r` (lines 248–250) to the stored EES values; the EES source never does (they are unknowns closed by the evaporator-effectiveness, condenser and swept-volume relations). That is why the file had 164 equations for 161 unknowns (the old §8 "redundant check equations" note). Making the native file faithful (removing the three lines) gives a square system whose 76-variable block does not converge yet (solver robustness). Decide whether to restore faithfulness now (status would become `failing`) or keep it as is until the solver handles it.
- **CSL-0158** `gas_pipe_insulated_pressure_drop`: its shipped companion CSV lost the commas and inch marks of the EES `Pipedata.LKT` column names (`ø88.9 Steel, Seamless`, `3/16" Copper (ø4.76)`), so `lookup('pipedata',1,pipe$)` fails with *column not found*; quoting the headers in the CSV fixes it (CoolSolve reads quoted headers correctly).
- Models that **worked around `CS-GAP-LOOKUP-PROC` without listing it** could now be simplified: CSL-0162 (BRINEPROP flattened) and the `ht` libraries `tube_bank_nusselt`, `tube_bank_dp_bell_delaware`, `free_conv_enclosed_and_jackets`, `shell_and_tube_sizing` (table functions skipped or transcribed as IF ladders).
- **CSL-0155** is `failing` (block 106, 151 variables, does not converge) and **CSL-0042** was verified against an independent RK4 only (no EES file): both choices of the re-check agents, to confirm.

---

## Appendix A. Full fix protocol (given to the fix agents)

### Protocol for fixing a CoolSolve library gap (read fully before starting)

Repositories:
- CoolSolve: `/home/sylvain/svn/CoolSolve` (C++17, branch `fix/library-gaps-2`, Release build in `build/`).
- CoolSolve Library: `/home/sylvain/svn/CoolSolve_Library` (do NOT touch it in a fix task; a separate re-check task handles it).

Mandatory reading before coding:
- `docs/contributing.md` (§1 principles, §2 workflow, §3 integration checklist). Follow it.
- The register row(s) of your item in `docs/model_library_support.md` (§4 gaps, §5 bugs; grep the ID), and §7 (format of
  closed entries).
- `tests/test_library_gaps.cpp` (style of the existing regression tests: helpers, tags, names `"<ID>: <what>"`).

Steps, in this order:
1. **Reproduce** the register reproducer with `build/coolsolve` on a small file in your scratch folder
   (`/tmp/claude-1000/-home-sylvain-svn-CoolSolve/9a4edc7b-2ae3-4dde-b675-6b24688bae79/scratchpad/<item>/`).
2. **Write the tests first** in `tests/test_library_gaps.cpp`, tag `[library-gaps][<short-tag>]`, names starting with the
   register ID. Cover: the register reproducer with the EES value, the realistic usage found in the blocked library
   models (look at the native `.eescode` files listed in the row — read-only), derivatives/Jacobian if the feature is
   numeric (AD must be analytic, no finite differences), and **negative tests** (wrong arity, malformed input → a clear
   error message). Build and check that the new tests FAIL before the fix.
3. **Implement** the fix with the simplest design that reuses existing code (no duplication; look for similar
   functionality first). Zero overhead for models not using the feature. Keep changes minimal and on-topic.
4. Build (`cd build && cmake --build . -j$(nproc) 2>&1 | grep -E 'error|warning: unused' ; echo built`) and iterate until
   `./coolsolve_tests "[library-gaps]"` passes.
5. **Full suites** (both must pass; do NOT run `[solver-robustness]`, too slow):
   `./coolsolve_tests` (baseline: all 379 test cases pass) and `./coolsolve_tests "[examples-comprehensive]"` (baseline:
   the test case passes; compare the per-example lines with the baseline in
   `/tmp/claude-1000/-home-sylvain-svn-CoolSolve/9a4edc7b-2ae3-4dde-b675-6b24688bae79/scratchpad/baseline_examples.txt`
   and `examples/test_examples.md` — no example may regress; report any change).
6. **Documentation** (contributing §3.12 and the checklist items that apply):
   - `docs/language_reference.md` for any new/changed syntax or built-in (only what is implemented and tested);
   - `docs/ees_vs_coolsolve.csv` rows of the feature (support column);
   - `gui/src/languages/ees.ts` for new built-in names/keywords (Monaco highlighting) — no GUI build needed;
   - `docs/model_library_support.md`: **move** the row from §4/§5 to the §7 *Closed gaps and bugs* table (same
     columns as the existing §7 rows: ID | What was wrong | Fixed in | Regression test | Notes). *Fixed in* =
     `post-v0.3.0, branch \`fix/library-gaps-2\`` (the supervisor adds the commit hash later). In *Notes*, list what
     changed for existing examples, and the library models that listed the ID (they are re-checked separately).
     Delete the §4/§5 row. Keep any still-open part of the gap as a new, narrower §4/§5 row if the fix is partial.
   - `README.md` only if a feature bullet / table changes.
7. **Do not commit** (the supervisor reviews and commits). Do not touch `examples/test_examples.md` or
   `examples/solver_robustness_report.md` beyond what the test run regenerates. Leave no stray files in the repo
   (scratch files go to the scratch folder).

Final report (concise, ≤ 40 lines): root cause; design of the fix and files changed; tests added (names) and evidence
that they failed before and pass after; full-suite results (counts) and any example whose result changed; docs updated;
anything left open, doubtful or out of scope (with a minimal reproducer).

Rules: valid EES behaviour is the reference (the EES manual/register evidence). Never weaken an existing test to make it
pass. If the fix turns out much larger than expected or needs a design decision, stop and report instead of improvising.

## Appendix B. Full re-check protocol (given to the re-check agents)

### Protocol for a library re-check after a CoolSolve fix (T-RECHECK) — read fully before starting

Repositories:
- CoolSolve Library: `/home/sylvain/svn/CoolSolve_Library` — you edit it. **Never commit** (rule D9), never delete
  anything outside the model folders you are re-checking, never copy student names.
- CoolSolve: `/home/sylvain/svn/CoolSolve` — **read-only** for you (docs and `tools/compare_solution.py`,
  `tools/ees_extract.py`). Do NOT build CoolSolve. Use the frozen binary given in your task (absolute path).

Mandatory reading: library `docs/model_workflow.md` §6 (gaps, blocked models, variants, re-checks) and §8 (definition of
done); `README.md` and `model.json` of each model you re-check; an earlier re-check for the style:
`models/fundamentals/combustion/cpbar_combustion_products/` (history action *"re-checked with CoolSolve fix/library-gaps
(T-RECHECK): blocked -> verified …"*) and the roadmap §8 line *"T-RECHECK `CSL-0005`, `CSL-0006`"*.

For every model listed in your task (all models whose `model.json` `missing_features` contains the closed ID):
1. Remove the closed ID from `missing_features`.
2. **If no gap is left**: copy the model folder to your scratch folder and solve the native main file there with the
   frozen binary (`cd <copy> && <binary> ./<main>.eescode`; add `-d` only to diagnose). If it solves, compare it with
   the reference the model was verified against (the `_coolsolve` variant's `.sol` on the common variables, and/or
   `python3 /home/sylvain/svn/CoolSolve/tools/compare_solution.py` against the EES source when the README documents
   it). If it agrees as well as the variant did:
   - status `verified` (or `runs` if the model has no EES reference), `verification` block updated: `date`
     2026-10-10, `coolsolve_version` = `0.3.0@<commit>` (strict format, the branch goes in `notes`), `max_rel_diff`, `reference`, `notes`
     (state that the native file is now the verified file);
   - write the native `.sol` into the model folder (the regression baseline; same format as produced by CoolSolve);
     copy a `.initials` only if the native file needs one to converge;
   - **clean the folder**: delete the `_coolsolve` variant files (`<name>_coolsolve.eescode/.sol/.initials`,
     `<name>_coolsolve-*.csv`) when the native file reproduces what the variant gave; keep a variant only if it still
     serves a purpose stated in the README (say which, in your report);
   - README: status row (CoolSolve line), *How to run*, *Results* (if it described the variant), *Verification*,
     *Conversion log* (add a line: re-checked with CoolSolve `fix/library-gaps-2` @ <commit>, gap closed, variant
     removed), *Limitations and CoolSolve gaps* (the gap is closed). Keep the style and length of the README;
   - `history`: append `{"date": "2026-10-10", "action": "re-checked with CoolSolve fix/library-gaps-2 @<commit>
     (T-RECHECK, <ID> closed): blocked -> verified; native file verified, _coolsolve variant removed", "by": "Claude
     (Sonnet) for S. Quoilin"}` (adapt the text to what you did).
   Native-file header: update its status/gap comment lines (comments only, never code) so that they no longer cite
   the closed ID as blocking; variant headers keep their history (they explain why the variant was made) but
   say the closed gap no longer applies when the variant is kept.
   If it does not solve or does not agree: keep `blocked`, keep the variant, and find out why (`-d` debug folder,
   minimal reproducer in a scratch file). If the cause is a CoolSolve problem not yet in the register
   (`/home/sylvain/svn/CoolSolve/docs/model_library_support.md` §4–§5), do NOT edit the register: describe it in your
   final report with a minimal reproducer (≤ 10 lines) and the evidence; list the model's remaining problem in the
   README *Limitations* section.
3. **If other gaps remain**: keep `blocked` and the variant; update the README sentences that cite the closed ID (it is
   now closed in CoolSolve `fix/library-gaps-2`) and the `history` (action *"re-checked … (<ID> closed); still blocked by
   …"*). Optionally confirm by running the native file that the closed gap no longer triggers (the first error should
   now come from a remaining gap) — report what you see.
4. After all models: run `python3 tools/build_index.py` (must pass; it regenerates `library.csv`, `library.json`,
   `CATALOG.md`, …) and `python3 tools/test_models.py --coolsolve <frozen binary> <the CSL ids you touched>` (every
   target OK, except the pre-existing failures listed in your task).
5. Add one line to the library `roadmap.md` §8 *Progress log* (end of the table): date 2026-10-10, step
   `T-RECHECK <ID>`, result (models re-checked, status changes, variants removed, open problems). Update the counts of
   roadmap §1 *Status at a glance* (verified / runs / blocked; blocked with a verified variant / not runnable) and the
   *CoolSolve gaps and bugs open/closed* lines if your change affects them.

Final report (≤ 40 lines): per model — before/after status, what was run, comparison result, files deleted/added;
build_index/test_models results; any new CoolSolve problem with its reproducer.

## Appendix C. Prepared prompts of the open items

Each prompt is removed when its item is done. The fix protocol it refers to is Appendix A (save it as `fix_protocol.md` in the session scratchpad and adapt the paths).

### C.6 — `prompt_6.md`

````text
You are fixing CoolSolve items. First read and follow the protocol file:
/tmp/claude-1000/-home-sylvain-svn-CoolSolve/9a4edc7b-2ae3-4dde-b675-6b24688bae79/scratchpad/fix_protocol.md
(The build tree is configured with CoolProp pinned (`COOLSOLVE_COOLPROP_TAG=5b9c32ac`); never re-run cmake with another
tag; if you add .cpp files, `cmake .` in `build/` keeps the cached tag. Take the current `./coolsolve_tests` pass count as
your baseline; the `[examples-comprehensive]` per-example table must stay identical to
/tmp/claude-1000/-home-sylvain-svn-CoolSolve/9a4edc7b-2ae3-4dde-b675-6b24688bae79/scratchpad/baseline_test_examples.md,
except for intended improvements that you report.)

Items (register /home/sylvain/svn/CoolSolve/docs/model_library_support.md §4–§5) — one family: **characters that EES
accepts in variable names and CoolSolve does not**. Each row has valid-EES evidence (files that EES solved):
1. `CS-GAP-NAME-SYMBOL` — `%` in names: `C%`, `H%`, `O%` (fuel mass fractions; CSL-0040, CSL-0080, CSL-0150,
   CSL-0152, example `orc_complex`).
2. `CS-GAP-NAME-PIPE` — `|` in names: `K|star_r_n`, `K|star_r` (CSL-0074, CSL-0116, CSL-0163).
3. `CS-GAP-ACCENTED-IDENT` — non-ASCII letters (UTF-8) and `°` in names: `T_a_réf`, `p_atm_réf`, `n°_90°elbow_reg`
   (CSL-0151, CSL-0158).
4. `CS-BUG-HASH-IDENT` — `#` inside a name (`Nu#_f`, beyond the built-in constants `g#`, `pi#`…) (CSL-0161, CSL-0111).
Check also, and report (fix if trivial and safe): `\` in names (`dU\dtau`, `time_h\day` — probably already supported),
a name starting with a non-ASCII letter (`é_x`), and `'` in names if EES allows it (probably not: it opens a string).
Find out the EES rule (EES help "Variable names" / "Variable information" — WebFetch
`https://fchartsoftware.com/ees/eeshelp/` pages; fchart.com may return 403) and implement the characters it allows
that do not conflict with CoolSolve's syntax; document the rule in `docs/language_reference.md` (identifiers section).
Names stay **case-insensitive** (EES), including accented letters if the helper `eesStringEqual` / the existing
case-insensitive comparator can do it consistently — check how variable names are lower-cased (IR, `CaseInsensitiveLess`)
and make sure `T_a_réf` and `T_A_RÉF` are the same variable or explain why not.

Every place that handles a variable name must accept them: the parser (`isIdentChar`, `isIdentifier`, the PEG grammar
if used, keyword detection at word boundaries, `$IF` conditions, CALL arguments, DUPLICATE, array names `X%[i]`,
string variables `name$`), the IR/structural analysis, the evaluator, `.initials` reading, `.sol` writing, the
`-integral.csv`/JSON/CSV outputs, the GUI/server JSON (UTF-8), the debug folder (`variables.md`, `incidence.md`…), and
**LaTeX output** (`src/ir.cpp` LaTeX rendering and `src/latex_report.cpp`: `%`, `#`, `|`, `\`, `°` and accents must be
escaped or rendered so that `pdflatex` compiles — `tests/test_latex_report.cpp` + compile one generated report with
`pdflatex -interaction=nonstopmode` if pdflatex is installed). Diagram/plot code that splits names (state suffixes) must
not break.

Blocked library models to look at (read-only; do not edit the library): CSL-0040, CSL-0074, CSL-0080, CSL-0111 (verified;
its README says a local `Nu#_f` was renamed — check the native/original form if documented), CSL-0116, CSL-0150, CSL-0151,
CSL-0152, CSL-0158, CSL-0161, CSL-0163. Find them with `grep -l '<ID>' /home/sylvain/svn/CoolSolve_Library/models/*/*/*/model.json`.
Run each native file (scratch copy of its folder) before/after and report the first remaining error; where a model was
blocked only by these gaps (the previous re-check found CSL-0151 solves once `é` is replaced and CSL-0161 once `Nu#_f` is
renamed), check that the native file now gives the variant's `.sol` values.

Tests: tags `[library-gaps][name-symbol]`, `[name-pipe]`, `[accented-ident]`, `[hash-ident]`: register reproducers with
EES values, names in every position (LHS, RHS, function/procedure arguments and locals, CALL outputs, arrays, `$IF`,
INTEGRAL state, `.initials`), case-insensitivity, `.sol` round trip, LaTeX escaping, and negative tests (characters that
remain invalid give a clear error).

**Extra regression check (parser change):** library regression `cd /home/sylvain/svn/CoolSolve_Library && python3
tools/test_models.py --coolsolve /home/sylvain/svn/CoolSolve/build/coolsolve` (pre-existing failures to ignore are given
in your task message) and a parse comparison (statement list / first error) of every `.eescode` of the library
(`models/**`), `examples/` and `tests/` between the "before" binary and your build — every change must be explained
(the previous agents built an AST dumper in the scratchpad `item1b/` and `item1c/pd/` — reuse it if convenient).
````
