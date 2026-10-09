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
| 1 | `CS-GAP-IF5` | 11: 0009, 0041, **0073**, 0074, 0075, 0080, 0081, 0153, **0155**, 0162, 0166 | ✅ | `59b2862` | 🔄 |
| 1b | `CS-BUG-ARG-DROP` (new, P1 silent: unparsable positional argument dropped) + close `CS-GAP-UNIT-SUBEXPR` / `CS-GAP-UNIT-NAMEDARG` (already parse since `59b2862`; tests needed) | – | ☐ | | – |
| 2 | `CS-GAP-IF-DIRECTIVE` (`$if/$ifnot/$else/$endif`) | 10: 0012, 0040, 0081, 0116, 0119, 0150, 0151, 0161, 0162, 0166 | ☐ | | ☐ |
| 3 | INTEGRAL family: `CS-GAP-INTEGRAL-LIMITS`, `CS-BUG-INTEGRAL-FACTOR`, `CS-BUG-INTEGRAL-TABLE-SEP`, `CS-BUG-INTEGRAL-MAXSTEPS`, `CS-BUG-INTEGRAL-TABLE-CASE` (then `CS-GAP-INTEGRAL-MULTIVAR`) | 0009, 0010, 0041, 0047, 0059, 0081 | ☐ | | ☐ |
| 4 | `CS-GAP-LOOKUP-PROC` | 6: 0072, 0078, 0079, 0116, 0158, 0160 | ☐ | | ☐ |
| 5 | `CS-GAP-CALL-EXPR-OUT` | 5: 0072, 0078, 0161, 0162, 0163 | ☐ | | ☐ |
| 6 | Identifier characters: `CS-GAP-NAME-SYMBOL` (`C%`), `CS-GAP-NAME-PIPE` (`K\|star`), `CS-GAP-ACCENTED-IDENT` | 0040, 0074, 0080, 0116, 0150, 0151, 0152, 0158, 0163 | ☐ | | ☐ |
| 7 | BrineProp family: `CS-GAP-ELSEIF-CHAIN`, `CS-GAP-UPPERCASE`, `CS-GAP-STRING-ARRAY` | 0072, 0078, 0079 | ☐ | | ☐ |
| 8 | `CS-GAP-INTERP-EES` + `CS-GAP-INTERP-EXTRAP` (silent) | 0011, 0035, 0149 | ☐ | | ☐ |
| 9 | Silent P1 bugs: `CS-BUG-PI-FUNCTION`, `CS-BUG-CONVERT-UNIT`, `CS-BUG-HASH-IDENT`, `CS-BUG-COMMON-PROC`, `CS-BUG-IMPLICIT-PROC`, `CS-BUG-LOOKUP-COL-ARG`, `CS-BUG-LOOKUP-COLNAME`, `CS-BUG-DUPLICATE-VAR-BOUND`, `CS-BUG-DUPLICATE-INDEX-EXPR`, `CS-GAP-GOTO`, `CS-BUG-HUMIDAIR-PROPS`, `CS-GAP-AIR-ENTHALPY-REF`, `CS-BUG-STRING-CASE` | 0075, 0080, 0081, 0118, 0123, 0124, 0127, 0129, 0131, 0157, 0158, 0159, 0161 | ☐ | | ☐ |
| 10 | Single-model items (by register priority): `CS-GAP-PSYCHRO-SAT`, `CS-GAP-FLUIDS-C8H18`, `CS-GAP-FLUIDS-ALIAS`, `CS-GAP-INTEGRAL-LOOKUP`, `CS-GAP-REPEAT-UNTIL-BARE`, `CS-BUG-LOOKUP-WRITE`, `CS-GAP-CALL-WARNING`, `CS-GAP-NDAY`, `CS-GAP-END-PROCEDURE`, `CS-GAP-CONVERT-UNQUOTED`, `CS-GAP-LOOKUPROW`, `CS-BUG-LOOKUP-STRING`, `CS-GAP-IFSTR`, `CS-GAP-ISIDEALGAS`, `CS-GAP-SUM-INDEXED`, `CS-GAP-GE-ARROW`, … | – | ☐ | | ☐ |
| 11 | `CS-GAP-INCLUDE` / `CS-FEAT-IMPORT` (separate feature) | 0011, 0074, 0116, 0149, 0162, 0163 | ⏸ | | |

## 5. Log

| Date | Step | Result |
|---|---|---|
| 2026-10-10 | Set-up | Branch `fix/library-gaps-2` created; baselines §2 recorded; frozen baseline binary in the session scratchpad (`coolsolve_baseline`). |
| 2026-10-10 | #1 `CS-GAP-IF5` | Fixed (`59b2862`): 3- and 5-argument `if` in `evaluateStandardFunction`; also unit-annotated literal arguments (`3600[s]`) no longer dropped, and `INTEGRAL` base expressions evaluate built-in functions instead of 0. 13 tests. Suites: 424/424 unit, examples identical to baseline. Native CSL-0073 solves, identical to its variant `.sol`; CSL-0155 parses but fails on block 106 (151 vars, NaN at start, not IF-related). New §5 row `CS-BUG-ARG-DROP`. |
