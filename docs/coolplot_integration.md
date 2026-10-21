# CoolPlot Integration Plan

**Status:** proposal (2026-10-20), not started.
**Reference:** [CoolPlot](https://github.com/CoolProp/CoolPlot) 0.4.0.dev0, `master` at `22937c8`.

This document plans how the **Diagrams** tab of the GUI adopts CoolPlot, the
property-diagram package of the CoolProp project, instead of growing its own plotting
code. It records the decision (§1), the principles (§2), the architecture (§3), how a
solved model becomes a cycle on a chart (§4), the phases (§5), the tests (§6) and the
decisions still open (§7).

---

## 1. Context and decision

### 1.1 The Diagrams tab today

`gui/src/components/ThermoDiagram.tsx` (≈450 lines, Plotly.js) draws:

- the saturation dome of one `Real` fluid in T-s, P-h, h-s or T-h, from
  `POST /api/v1/coolprop/saturation` (`src/server.cpp`: N temperatures between the
  triple and the critical point, liquid and vapour arrays in SI);
- the solved state points, grouped by a name heuristic (`T_1`, `P_1`, `h_1` → state `1`)
  from `GET /api/v1/variables/inferred`;
- one array path (`X[i]`, `Y[i]`), optionally closed.

There are no isolines, no unit choice (°C, Pa, J/kg only), no labels, and only `Real`
fluids are offered. The wanted extensions are listed in [gui.md §14.4](gui.md) (isolines,
units, hover, persistence) and in the feature register of
[model_library_support.md §3](model_library_support.md): `CS-FEAT-DIAGRAM-IDEAL`,
`CS-FEAT-PSYCHRO`, `CS-FEAT-OVERLAY-STATES`, `CS-FEAT-PLOT-CLI`.

### 1.2 CoolPlot

CoolPlot (Python, MIT, ≈13,500 lines in `src/`) draws log(p)-h, T-s, h-s, p-s, p-ρ, T-ρ
and p-T diagrams of pure fluids and blends, with cycles on top, and psychrometric
charts. What matters for CoolSolve:

- **Numerics already solved**: isolines are swept in T or p with fast flash pairs (its
  architecture review measured 35 s with 37 % failed points for the legacy approach on a
  T-s chart of water, 0.54 s with none for the new one); saturation states are inserted
  where an isoline crosses the dome; failed points are gaps, never interpolated; blends
  use a tracer ported from CoolProp 8 and the dome comes from the phase envelope.
- **Presentation already designed**: unit systems (`SI`, `KSI`, `EUR`), region presets
  (`DEF`, `ACHP`, `ORC`), round isoline values, themes, label placement on the frame
  edges with keep-out zones around cycles, cycle paths that follow the physics
  (isobars bend on the dome, isenthalps, isentropes).
- **A pipeline of plain data**: `spec → plan → curves → compose → layout → scene →
  renderer`. Three records cross the layers and are JSON contracts
  (`docs/contracts.md`): the chart spec, the curve record and the scene. That document
  expects ports to other languages: "a C++ numerics library produces curve records, a
  C++ presentation library turns a spec plus curve records into a scene, and a renderer
  in any language draws a scene".
- **A plotly backend** (`render/plotly.py`): a Python plotly figure is the JSON that
  Plotly.js reads, so the GUI renderer can follow it almost line by line.

### 1.3 Options considered

| Option | Verdict |
|---|---|
| **A. Run CoolPlot's Python code** (embedded CPython, a Python process next to the server, or Pyodide in the browser) | Rejected. The Windows build is one self-contained `.exe`; Python, numpy and the CoolProp wheel would have to ship with it (or a large Pyodide download, with no certainty that the CoolProp package runs there). The curves would come from a second CoolProp, not the pinned one, and without CoolSolve's reference-state offsets (ideal gases are on the EES/JANAF basis, `IdealGasFluid::getReferenceState`): solved states would sit beside the lines. CoolSolve's own fluids (`NH3H2O`, LiBr) could never be drawn. |
| **B. Keep developing our own plots** | Rejected. Every item of §1.1 has been designed, measured and reviewed in CoolPlot; doing it again means meeting the same problems again. |
| **C. Port CoolPlot's design and contracts to C++ and TypeScript** | **Chosen.** Numerics and presentation in the C++ server, through CoolSolve's fluid layer; the scene drawn by Plotly.js in the GUI. CoolPlot stays the reference implementation and the source of the test fixtures. |

---

## 2. Principles

1. **CoolPlot is the reference.** Contracts, item ids, groups, roles, styles, texts and
   numerical methods follow CoolPlot. CoolSolve follows one recorded CoolPlot commit
   (the header of this document), as it follows one pinned CoolProp commit; moving to a
   newer one is a deliberate change (§6.4).
2. **Properties come from CoolSolve's fluid layer**, never from a second property
   library: the pinned CoolProp, the same reference states as the solver
   (`Fluid::getReferenceState`), the same fluid names (`FluidRegistry`). A state point of
   a solved model lies exactly on the lines of its own properties.
3. **The contracts are CoolPlot's, unchanged.** The GUI sends a chart spec in CoolPlot's
   format and receives a scene in CoolPlot's format. A feature CoolPlot lacks (the T-h
   diagram, model units) is proposed upstream first; a CoolSolve-only extension is a
   documented, namespaced addition, never a silent change of meaning.
4. **Port by layer, with the same boundaries and names.** Each C++ module mirrors one
   CoolPlot module, so that a CoolPlot change can be found and ported. Each ported file
   names its source file and the CoolPlot commit, and keeps CoolPlot's MIT notice
   ("Copyright 2014-2019 Jorrit Wronski and the other developers").
5. **No Python at run time or at test time.** Python is needed only to regenerate the
   test fixtures with a script in `tools/` (§6.1); the fixtures are committed.
6. **Not solver core.** The work touches the server, a new `diagram` module and the GUI,
   not the parser, evaluator or solvers. It is still developed on a branch
   (`feat/coolplot-diagrams`) and reviewed before merging.

---

## 3. Architecture

### 3.1 Data flow

```
GUI (Diagrams tab)                          coolsolve server (C++)
------------------                          -----------------------------------------------
chart spec (CoolPlot JSON)    POST          spec  --plan-->  curve requests
+ plot size in pt             -------->           --engine-> curves (SI, process-wide cache)
+ view (zoom window)          /api/v1/                         ^  thermo: providers over
                              diagram/scene                    |  FluidRegistry / AbstractState
                                                  --compose-> scene (display units, styles)
Plotly.js  <-- TS renderer <-- scene JSON <------ --layout--> scene (labels for the plot size)
(port of render/plotly.py)
```

- One endpoint, `POST /api/v1/diagram/scene`, takes the chart spec, the plot size and
  the view, and returns the scene and the curve diagnostics (failed points, fallbacks,
  lines that could not be computed and why).
- A second one, `GET /api/v1/diagram/options`, lists what the GUI offers: diagram types
  per fluid kind, unit systems, region presets, isoline families per diagram.
- `GET /api/v1/coolprop/fluids` and `GET /api/v1/variables/inferred` stay.
  `POST /api/v1/coolprop/saturation` is removed once the new tab replaces the old one.
- Curves are cached process-wide (the server has several sessions), keyed as in CoolPlot
  by fluid, CoolProp fingerprint and reference state, diagram, family, value, point
  count and calculation settings; the cache is bounded and thread-safe. Changing units,
  theme, labels or the view recomputes no property: only compose and layout run again.

### 3.2 Where compose and layout run (decision D1)

The plan puts compose and layout in C++ with the numerics, and keeps the browser side to
a renderer:

- one port language for everything but the renderer;
- the scene is the single boundary, and CoolPlot can produce the expected scene of any
  spec, so the whole C++ side is tested against fixtures (§6);
- it is the split CoolPlot's contracts anticipate, so the C++ code can later serve
  CoolPlot itself (D2).

The cost is a request on every unit switch, zoom or resize. On the desktop the server is
local; on the online demo the request carries a few kB and no property calculation (the
curves are cached). The alternative, compose and layout in TypeScript, avoids the
requests but splits the port across two languages and leaves the TypeScript layers
without a reference to test against.

### 3.3 Module map

New C++ module `include/coolsolve/diagram/`, `src/diagram/`; new GUI folder
`gui/src/components/diagram/`.

| CoolPlot module | Lines | CoolSolve | Phase |
|---|---|---|---|
| `core/quantities.py`, `units.py`, `diagram_type.py`, `limits.py`, `numbers.py`, `state.py` | 724 | `diagram/core.{h,cpp}` | 1 |
| `spec/fluid.py`, `diagram.py`, `settings.py` (+ JSON round trip) | 829 | `diagram/spec.{h,cpp}` | 1 |
| `thermo/provider.py`, `fluid.py`, `isolines.py`, `limits.py` | 1,245 | `diagram/thermo.{h,cpp}` (providers over `FluidRegistry`) | 1 |
| `engine/plan.py`, `engine/cache.py` | 240 | `diagram/engine.{h,cpp}` | 1 |
| `style/values.py`, `theme.py` | 747 | `diagram/style.{h,cpp}` | 2 |
| `compose/builder.py`, `scene.py` | 307 + 445 | `diagram/compose.{h,cpp}`, `diagram/scene.{h,cpp}` | 2 |
| `render/plotly.py`, `render/ticks.py` | 695 | `gui/src/components/diagram/plotlyRenderer.ts` | 2 |
| `core/cycle.py`, `thermo/cycles/paths.py`, `thermo/process.py` | 402 | `diagram/cycles.{h,cpp}` | 3 |
| `compose/layout.py` | 552 | `diagram/layout.{h,cpp}` | 4 |
| `thermo/tracer.py` (+ envelope dome in `fluid.py`) | 799 | `diagram/tracer.{h,cpp}` | 5 |
| `thermo/humid_air.py` | 386 | `diagram/humid_air.{h,cpp}` | 5 |

Not ported:

| CoolPlot part | Why |
|---|---|
| `engine/background.py` | Works around Python's interpreter lock with threads and processes. The server computes a scene in the request; coarse-then-fine can be added later if a chart is measured to be slow. |
| `core/events.py`, `interaction.py` | Dragging cycle states edits a CoolPlot cycle; in CoolSolve the states come from the solve. Zoom and resize are handled by the renderer (§5, phase 4). |
| `render/mpl.py`, `svg.py`, `web.py`, `pdf_reader.py`, `export.py` | Plotly.js is the only renderer. Print-ready PDF goes through CoolPlot itself (§5, phase 6). |
| `thermo/cycles/models.py` | CoolPlot's cycle models compute a cycle; in CoolSolve the model is the cycle. |
| `diagram.py`, `chart.py`, `psychrometric.py` | Python facades; the GUI builds the spec directly. |

---

## 4. From a solved model to a cycle on the chart

CoolPlot draws a cycle from numbered states (each fixed by two properties), chains of
state numbers (`loops`) and paths between them (`"auto"`: isobar if the pressures are
equal, isenthalp, isentrope, isotherm, otherwise a polytropic path). The question is how
the GUI finds the states of a solved model.

### 4.1 Sources of states

1. **Arrays** (today's array overlay): `P[i]`, `h[i]` give state `i`; the loop is the
   index order. CoolPlot accepts cycles given as arrays.
2. **Inferred states**: variables used as arguments or results of the same property
   calls (`inferVariables`, `src/variable_inference.cpp`) are grouped into states (§4.2).
3. **Picked by the user** (`CS-FEAT-OVERLAY-STATES`): a state editor in the Diagrams
   tab where each state is a set of variables (`T_su_cp`, `P_su_cp`, …) and the loops are
   chains of states. This is the fallback when the automatic grouping is wrong, and the
   way to name states.

The cycle definition (states as variable names, loops, paths) is saved in the ZIP bundle
with the chart spec (`diagram.json`), not in the `.eescode`: EES keeps its plots in the
GUI too, and the model file stays EES code.

### 4.2 Grouping inferred variables into states (decision D6)

Today's name heuristic merges everything that shares a suffix. A rule based on the
equations: each property call on a fluid involves a set of state variables (inputs and
output); two calls describe the same state when they share **two variables of
independent properties** (state postulate for a pure fluid). `h_2 = enthalpy(R134a,
P=P_2, s=s_2)` and `T_2 = temperature(R134a, P=P_2, h=h_2)` share `P_2` and `h_2`: one
state. `T_sat = temperature(R134a, P=P_cond, x=1)` shares only `P_cond` with the
condenser states: a separate state. T and P are not independent inside the dome, so a
pair (T, P) alone does not merge two calls. This needs a test on the library models
before it replaces the heuristic.

### 4.3 Values

States go into the spec in SI (CoolPlot's `CycleSpec` is SI); the GUI converts the
model units (°C, Pa, J/kg) on the way. Each state is drawn from its solved values, never
recomputed. One chart belongs to one fluid; a model with several fluids gets one chart
per fluid, chosen in the fluid list as today.

---

## 5. Phases

Each phase ends with the Diagrams tab usable and the old code it replaces removed.

### Phase 0 — Preparation

- Settle D2 with the CoolPlot maintainers (where the C++ code lives) and record the
  reference commit.
- `tools/coolplot_fixtures.py`: from a list of chart specs, run CoolPlot at the
  reference commit and write the curve records, the scenes and the plotly figures to
  `tests/fixtures/coolplot/` (§6.1).
- Choose the fixture specs: R134a, water, CO2, R290 in every diagram type; `ACHP` and
  `DEF` regions; all units; one cycle per diagram type.

### Phase 1 — Curves for pure fluids

- `diagram/core`, `spec`, `thermo`, `engine` (§3.3): quantities and units, diagram types
  and their isoline families, regions, round values, isolines by flash with saturation
  states inserted, the dome, the plan and the cache.
- Providers: `Real` fluids through CoolSolve's cached `AbstractState`, with the fluid's
  reference state.
- Tests: curve fixtures and the invariants of §6.2.

About 3,000 lines of Python to port; no change in the GUI yet.

### Phase 2 — Scene, renderer, new Diagrams tab

- `diagram/style`, `scene`, `compose`; `POST /api/v1/diagram/scene`,
  `GET /api/v1/diagram/options`.
- `plotlyRenderer.ts`: scene → traces and annotations, points to pixels (4/3), CoolPlot's
  tick placement and number format, log axes.
- Diagrams tab: fluid, diagram type (the seven CoolPlot types, T-h depending on D3),
  units, region preset, isoline families (count, rounding, on/off), theme. Today's state
  and array overlays stay, drawn as CoolPlot markers and lines.
- The chart spec is saved in the ZIP bundle and restored (gui.md §14.4 "Diagram state
  persistence").
- `exportPlots.ts` keeps exporting `thermo_diagrams.png` for the LaTeX report (the new
  tab keeps the `thermo-diagram` class it looks for).
- Remove `POST /api/v1/coolprop/saturation` and the old dome code.

### Phase 3 — Cycles from the model

- `diagram/cycles`: states, loops, automatic paths, state numbers, one legend entry,
  keep-out zones.
- The three state sources of §4.1, the grouping rule of §4.2, the state editor; hover
  shows the variable names and values of each state, as today.
- Closes `CS-FEAT-OVERLAY-STATES`.

### Phase 4 — Labels

- `diagram/layout`: edge labels rotated along the lines, no overlaps, keep-out zones
  around cycles, stable positions between updates; margin and inline modes.
- Zoom and resize: the renderer reports the new view (`plotly_relayout`) and plot size,
  and the GUI asks for a new scene; only compose and layout run.

### Phase 5 — More fluids

- Blends and CoolProp mixtures: the tracer and the dome from the phase envelope. The
  pinned CoolProp is 8.1, which the tracer requires.
- Ideal gases (`CS-FEAT-DIAGRAM-IDEAL`): T-s, h-s and p-ρ (a p-v type proposed upstream
  if wanted) over the gas region, with absolute limits, through `IdealGasFluid` and its
  EES/JANAF reference state.
- Humid air (`CS-FEAT-PSYCHRO`): the psychrometric chart of `AirH2O` models, with air
  processes drawn from the solved states.
- Incompressibles, with absolute limits.
- Later, through the same provider interface: `NH3H2O` and LiBr-H2O charts (Dühring,
  h-x). CoolPlot has no such charts; they would be proposed upstream as new chart kinds.

### Phase 6 — Outputs

- *Export chart spec*: the chart spec with the model's cycle, as JSON that CoolPlot
  reads (`DiagramSpec.from_dict`), for print-ready PDF and matplotlib figures in Python.
- `CS-FEAT-PLOT-CLI` (`coolsolve --diagram spec.json model.eescode`, for README figures):
  needs a non-browser renderer, i.e. a port of `render/svg.py` (830 lines) to C++. Only
  if wanted.

---

## 6. Tests and keeping in step with CoolPlot

### 6.1 Fixtures

`tools/coolplot_fixtures.py` writes, for each spec, CoolPlot's curve records, scene and
plotly figure (`PlotlyRenderer`) to `tests/fixtures/coolplot/`, with the CoolPlot commit,
the Python CoolProp version and the date in a header file. The fixtures are committed;
the C++ and TypeScript tests read them without Python.

- **Curves** (C++, Catch2 `[diagram]`): same lines, same point counts and gaps;
  coordinates within a relative tolerance, because the Python CoolProp wheel is not the
  pinned C++ commit. This is a comparison with another implementation, so the tolerance
  may be set from measured differences; no numerical tolerance of the port itself is
  relaxed to make a fixture pass.
- **Scenes** (C++): ids, groups, roles, z order, styles and texts exactly; coordinates
  as for curves.
- **Plotly figures** (TypeScript): for a fixture scene, `plotlyRenderer.ts` must produce
  the keys CoolPlot's plotly backend sets, with the same values. The GUI has no test
  runner yet (D5).

### 6.2 Invariants

Independent of CoolPlot, on every fixture fluid:

- the dome closes at the critical point; liquid and vapour lines meet there;
- an isoline crossing the dome contains the saturation states; inside it, an isobar of
  a T-s chart and an isotherm of a p-h chart are horizontal;
- no point is drawn where the flash failed;
- consistency with the solver: a state solved by CoolSolve lies on the isolines of its
  own properties (checked for a `Real` fluid and for an ideal gas with a JANAF reference
  state).

### 6.3 Performance

Measure, in a Release build, the cold time of the default chart of each fixture fluid
and the warm time of compose plus layout; set the budget from the measurements and add
it to the test suite ([contributing.md §5](contributing.md)).

### 6.4 Moving to a newer CoolPlot

Bump the reference commit, regenerate the fixtures, port what the failing tests point
to, update the header of this document and note it in [versions.md](versions.md).

---

## 7. Open decisions

| ID | Question | Proposal |
|---|---|---|
| D1 | Compose and layout in C++ (server) or in TypeScript (browser)? | C++ (§3.2). |
| D2 | Does the C++ code live in CoolSolve (`src/diagram/`) or in a library shared with CoolPlot? | Start in CoolSolve with CoolPlot's module boundaries, so it can be extracted. To agree with the CoolPlot maintainers: CoolPlot's `docs/contracts.md` points to `docs/v1_review.md` §6 for a C++ port, a file not in the current CoolPlot checkout. |
| D3 | T-h diagram: CoolPlot has none, today's tab has one (useful for heat exchangers). | Propose `Th` upstream; keep it CoolSolve-only meanwhile. |
| D4 | Default units of the Diagrams tab. | `EUR` (°C, bar, kJ/kg), as printed charts; a "model units" system (°C, Pa, J/kg) proposed upstream as a named system. |
| D5 | A TypeScript test runner (e.g. Vitest, dev dependency only) for the renderer fixtures. | Add it in phase 2. |
| D6 | Automatic grouping of states (§4.2) or user picking only. | Both: automatic by default, the editor to correct and name. Validate the rule on the library models first. |

---

## 8. Documentation to update when the work lands

- [gui.md](gui.md): §4.7 (endpoints), §6.8 (Diagrams tab), §8.3 (variable inference for
  diagrams), §14.4 (items done).
- [model_library_support.md §3](model_library_support.md): `CS-FEAT-OVERLAY-STATES`,
  `CS-FEAT-DIAGRAM-IDEAL`, `CS-FEAT-PSYCHRO`, `CS-FEAT-PLOT-CLI`.
- [contributing.md](contributing.md): file map and the CoolPlot reference commit next to
  the pinned CoolProp.
- README: features, documentation table, future work.
