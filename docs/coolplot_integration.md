# CoolPlot Integration Plan

**Status:** proposal, revised 2026-10-31 (first version 2026-10-20), not started.
**References:** [CoolPlot](https://github.com/CoolProp/CoolPlot) 1.5.0, `master` at `16895a9`
(the `v1.5.0` tag is not pushed yet); [CoolProp](https://github.com/CoolProp/CoolProp) `master`
from `fc91a27f` on (`CoolProp::Plot`, CoolProp PRs #3512 and #3519, in no CoolProp release yet).

This document plans how the **Diagrams** tab of the GUI adopts CoolPlot, the
property-diagram package of the CoolProp project, instead of growing its own plotting
code. It records the decision (§1), the principles (§2), the architecture (§3), how a
solved model becomes a cycle on a chart (§4), the phases (§5), the tests (§6) and the
decisions still open (§7).

**What changed in this revision.** The first version planned to port CoolPlot 0.4 to
C++ and TypeScript inside CoolSolve: about 6,000 lines of Python (numerics, compose,
layout, plotly renderer). CoolPlot has since done both C++ ports itself (§1.2): the
numerics are in CoolProp as `CoolProp::Plot`, the presentation is a C++17 library,
coolplot-core, with SVG and PDF writers. This version uses the two libraries as they
are, and plans only the code between them (about 1,000 lines, §1.3) and the GUI. Four of
the six open decisions of the first version are settled (§7.2); the renderer is open
again, because coolplot-core writes SVG itself (§3.2).

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

### 1.2 CoolPlot 1.5

CoolPlot went from 0.4.0.dev0 (`22937c8`) to 1.5.0 in 194 commits. It draws log(p)-h,
T-s, h-s, p-s, p-ρ, T-ρ, p-T and any other pair of fluid axes (`Th`), of pure fluids,
blends and incompressibles, with cycles on top, and psychrometric charts (ASHRAE and
Mollier h,x) with air processes. The pipeline is still `spec → plan → curves → compose
→ layout → scene → renderer`, with its records written down in CoolPlot's
`docs/reference/contracts.md`. What matters for CoolSolve is that two parts of it now
exist in C++, as CoolPlot's review planned (`docs/contributing/v1_review.md` §6).

**Numerics in CoolProp: `CoolProp::Plot`.** `include/CoolProp/plotting/CoolPropPlot.h`
and `src/CoolPropPlot.cpp` (2,780 lines with its Catch2 tests) are a port of CoolPlot's
`thermo/isolines.py`, `tracer.py` and `stepping.py`, merged into CoolProp `master` on
2026-10-21 (#3512; #3519 on 2026-10-24 interpolates the envelope crossings with a cubic).
It replaces the dormant 820-line `CoolPropPlot.cpp` of 2025, which is what CoolSolve's
pinned commit (`75af7816`, 2026-10-09) still has.

- Lines are stepped with steps chosen from their shape on the page; nothing is computed
  outside the axes; isotherms and isobars are followed through the dome in vapour
  quality; other lines get the saturation state inserted where they cross the dome;
  mixtures are traced (Newton in T and ρ from the previous point) and their dome comes
  from the phase envelope; a failed point is a NaN gap, never interpolated.
- Output: one `Isoline` per value, `x`, `y` in SI with the diagnostics of CoolPlot's
  curve record (`n_failed`, `n_fallback`, `n_inserted`, `n_evaluated`, `n_outside`,
  `method`, `reason`).
- API: `PropertyPlot(fluid_string, ykey, xkey, TPLimits)` with the presets `Def`,
  `Achp`, `Orc`, or with an absolute T-p box (an incompressible needs one);
  `calc_isolines(key, values, CalcSettings)`; `isoline_range(key)`; `tp_box()`;
  `supported_isoline_keys()`; `value_at(key, x, y)` (the state under the cursor).
- A plot builds its own `AbstractState` from the fluid string, so it takes the backend
  (`HEOS`, `TTSE&HEOS`, `INCOMP`, …) and the global configuration of the CoolProp it is
  linked with. Its calls hold a lock; lines are computed in parallel with one plot per
  thread.
- CoolPlot's Python sends its fluid lines there when the installed CoolProp has the port,
  and checks it point for point against its Python numerics (`tests/test_native.py`).
  Humid-air lines and the states of cycles are not part of the port (CoolPlot's
  `docs/explanation/numerics.md`).

**Presentation in C++: coolplot-core.** `cpp/` of the CoolPlot repository, version
1.5.0, ≈7,200 lines of C++17 plus nlohmann/json 3.12 (vendored); it uses neither
CoolProp nor a plotting library.

- Ports `coolplot.core` (number format, quantities, unit systems `SI`, `KSI`,
  `engineering`, `IP`, diagram types and families), `coolplot.style` (themes `default`,
  `dark`, `print`, `coolpack`), `coolplot.scene`, `coolplot.compose` (builder, layout,
  rails, the automatic legend and badge corners) and the SVG writer, and adds a PDF writer
  with TrueType subsets.
- Its input is a fourth record, the **compose record** (`cpp/include/coolplot/compose.h`):
  the chart spec, the calculation domain (`tp_box`, `ranges_SI`), the curves with the
  request each was computed for (`item_id`, `key`, `value_SI`, `points`, `labels`,
  `anchor_SI`), the cycles as resolved states and chains of states (the paths already
  computed), the plot size, the previous labels (stable positions between updates) and
  the badge lines (`sources`).
- Interfaces: C++ (`coolplot::compose`, `scene_to_svg`, `scene_to_pdf`), a C ABI of
  JSON text in and JSON, SVG or PDF out (`coolplot_compose`, `coolplot_scene_to_svg`,
  `coolplot_scene_to_pdf`), and a nanobind module for Python.
- Scenes and SVG documents are byte for byte those of the Python code (golden files in
  `cpp/tests/data`; CI runs CoolPlot's whole Python suite through the core). Measured by
  CoolPlot: 12 ms to compose a p-h chart with labels, JSON reading and writing included.
  CI builds it on Linux with GCC only.
- `src/coolplot/render/static/scene_svg.js` (657 lines) is the same SVG writer in
  JavaScript. The plotly backend (`render/plotly.py`, 833 lines) has no port.

Features added since 0.4 that the GUI gets with the core: labels on rails and a legend
placed inside the plot area by itself (the defaults since 1.2), lines placed through
saturation states (isentropes every 10 K of saturation temperature, `anchor`), state
labels on cycles, symbols written with their own characters (ρ, °C), the `IP` units and
`engineering` (the old `EUR`, still read from specs), one colour per quantity checked for
colour vision deficiencies, and a badge naming the equation of state ("R290: Lemmon et
al. (2009)", on by default since 1.5).

**Not in either library** (the "glue"):

| Step | CoolPlot (Python) | Lines | Needs CoolProp |
|---|---|---|---|
| T-p box from the spec's limits, ranges of the line families | `thermo/limits.py`, `thermo/provider.py` (`chart_domain`) | ≈130 | yes; `PropertyPlot` does the presets and absolute boxes, the relative bounds of a spec (`{"relative": 1.01}`) are resolved by the caller |
| Values of the lines placed through saturation states | `thermo/provider.py` (`saturation_anchor_values`), `engine/engine.py` | ≈80 | yes |
| Plan: which lines at which values (round values), item ids | `engine/planning.py` | 233 | no |
| Curve cache, compose record | `engine/cache.py`, `engine/engine.py`, `corelib.compose_record` | ≈250 | no |
| Cycle states and the paths between them | `thermo/cycles/paths.py`, `thermo/process.py`, `core/cycle.py` | ≈460 | yes |
| Humid-air lines | `thermo/humid_air.py` | 569 | yes (`HAPropsSI`) |
| Badge lines: references of the equation of state, versions | `thermo/references.py`, `engine.source_lines` | ≈440 | yes (BibTeX keys and bibliography) |

### 1.3 Options considered

| Option | Verdict |
|---|---|
| **A. Run CoolPlot's Python code** (embedded CPython, a Python process next to the server, or Pyodide in the browser) | Rejected. The Windows build is one self-contained `.exe`; Python, numpy and the CoolProp wheel would have to ship with it (or a large Pyodide download, with no certainty that the CoolProp package runs there). The curves would come from a second CoolProp, not the pinned one, and without CoolSolve's reference-state offsets (ideal gases are on the EES/JANAF basis, `IdealGasFluid::getReferenceState`): solved states would sit beside the lines. CoolSolve's own fluids (`NH3H2O`, LiBr) could never be drawn. |
| **B. Keep developing our own plots** | Rejected. Every item of §1.1 has been designed, measured and reviewed in CoolPlot. |
| **C. Port CoolPlot to C++ and TypeScript inside CoolSolve** (first version of this plan) | Superseded: CoolPlot made both C++ ports itself, tested against its Python code. A second port would drift from them. |
| **D. Link `CoolProp::Plot` and coolplot-core, write the glue** | **Chosen.** The numerics come with the pinned CoolProp, the same library and configuration the solver uses; coolplot-core is built from source like the other dependencies; CoolSolve writes the glue of §1.2 (about 1,000 lines of C++, less if parts of it go upstream, D2) and the Diagrams tab. |

---

## 2. Principles

1. **CoolPlot is the reference.** Contracts, item ids, groups, roles, styles, texts and
   numerical methods are CoolPlot's. CoolSolve now follows two pinned versions: CoolProp
   (the numerics, already pinned for the solver) and CoolPlot (coolplot-core and the
   reference for the glue, the header of this document). Moving either is a deliberate
   change (§6.4).
2. **Properties come from CoolSolve's fluid layer**, never from a second property
   library. `CoolProp::Plot` is the pinned CoolProp linked into `coolsolve`, with the
   process's configuration (superancillaries) and CoolProp's reference states, which are
   those of the solver for `Real` fluids (`Fluid::getReferenceState` gives no offset).
   The fluid string and the backend come from `FluidRegistry` and the model's
   `coolpropBackend`. Fluids that CoolProp alone does not describe (ideal gases on the
   EES/JANAF basis, `NH3H2O`, LiBr-H2O) get CoolSolve providers that write curve records
   themselves. A state point of a solved model lies exactly on the lines of its own
   properties.
3. **The contracts are CoolPlot's, unchanged.** The GUI sends a chart spec in CoolPlot's
   format; the server builds a compose record and coolplot-core returns the scene. A
   feature CoolPlot lacks is proposed upstream first; a CoolSolve-only extension is a
   documented, namespaced addition, never a silent change of meaning.
4. **Use the libraries, do not copy them.** `CoolProp::Plot` and coolplot-core are linked
   as they are, never patched in CoolSolve: a defect or a missing feature is fixed
   upstream and reaches CoolSolve by moving a pin. The glue mirrors CoolPlot's Python
   modules (names, boundaries, order of the steps), names its source file and the
   CoolPlot commit, and keeps CoolPlot's MIT notice ("Copyright 2014-2019 Jorrit Wronski
   and the other developers"), so that it can move upstream (D2).
5. **No Python at run time or at test time.** Python is needed only to regenerate the
   test fixtures with a script in `tools/` (§6.1); the fixtures are committed.
6. **Not solver core, except the CoolProp pin.** The work touches the server, a new
   `diagram` module and the GUI, not the parser, evaluator or solvers; it is developed on
   a branch (`feat/coolplot-diagrams`) and reviewed before merging. Moving the CoolProp
   pin to a commit with `CoolProp::Plot` changes property values for the solver too, so
   it is a change of its own (§5, phase 0).

---

## 3. Architecture

### 3.1 Data flow

```
GUI (Diagrams tab)                     coolsolve server (C++)
------------------                     ------------------------------------------------------------
chart spec (CoolPlot JSON)   POST      spec --> domain, plan --> curve requests
+ plot size in pt            ------->          curves <-- CoolProp::Plot (pinned CoolProp), cached
+ view (zoom window)         /api/v1/          cycles: states, paths <-- CoolSolve fluid layer
+ labels of the shown scene  diagram/          badge lines
                             scene     = compose record --coolplot-core--> scene --> SVG
SVG + interaction layer  <-- scene JSON, SVG, diagnostics <-----------------
(zoom, hover, state names)
```

- `POST /api/v1/diagram/scene` takes the chart spec, the plot size, the view and the
  label items of the scene on screen, and returns the scene, its SVG (D1) and the curve
  diagnostics (failed points, fallbacks, lines that could not be computed and why).
- `GET /api/v1/diagram/options` lists what the GUI offers: diagram types per fluid kind,
  unit systems, region presets, themes, label modes, line families per diagram.
- `POST /api/v1/diagram/state` returns the full state at a point of a chart
  (`PropertyPlot::value_at`), for the hover of phase 2.
- `GET /api/v1/coolprop/fluids` and `GET /api/v1/variables/inferred` stay.
  `POST /api/v1/coolprop/saturation` is removed once the new tab replaces the old one.
- Two caches, process-wide (the server has several sessions), bounded and thread-safe:
  `PropertyPlot` objects, which hold the critical point and the phase envelope of a
  fluid, keyed by fluid string, axes and T-p box (a few per key, so that sessions do not
  wait on one lock); and curves, keyed as in CoolPlot by fluid string, CoolProp git
  revision and configuration (backend, superancillaries), diagram, family, value and
  calculation settings. Changing units, theme, labels, view or plot size computes no
  property: only the compose step runs again.

### 3.2 How the GUI draws the scene (decision D1)

coolplot-core composes the scene and lays out its labels for the plot size and the
view, so every zoom and resize asks the server for a new scene in any case. What remains
is how the browser draws it:

| Option | For | Against |
|---|---|---|
| **a. The SVG of coolplot-core**, written by the server next to the scene | No renderer to write or keep in step: the document is CoolPlot's, byte for byte, tested upstream. The same writer gives the CLI and the report their figures, with the PDF writer next to it. Each item is a `<g data-item="<id>">` with its role and group as classes, and the root carries the axis limits and log flags, so a thin TypeScript layer maps the pointer to data coordinates and finds the state markers. | Zoom (a rubber band that sends a new view), reset, hover and the tooltip are written by us. |
| b. The scene drawn in the browser by CoolPlot's `scene_svg.js` | The same document, without SVG in the response. | Nothing gained over (a): the round trip stays for the labels; one more vendored file to keep in step. |
| c. A TypeScript port of `render/plotly.py` (833 lines) | Plotly's zoom, pan, hover and PNG export, as in the other plots of the GUI. | A port with no reference implementation in TypeScript to test against; plotly's limitations (annotations always above traces, labels not clipped, no line caps); after a client-side zoom the labels are wrong until the new scene arrives. |

Proposal: (a).

### 3.3 Module map

New C++ module `include/coolsolve/diagram/`, `src/diagram/`; new GUI folder
`gui/src/components/diagram/`. Each glue module mirrors the CoolPlot module in the
third column.

| CoolSolve | Does | Mirrors (CoolPlot) | Phase |
|---|---|---|---|
| `diagram/fluid.{h,cpp}` | the CoolProp fluid string of a model's fluid (`FluidRegistry` name, backend, mixture fractions); which fluids the port takes, which need a CoolSolve provider | `thermo/provider.py`, `thermo/native.py` (`fluid_string`) | 1 |
| `diagram/domain.{h,cpp}` | the T-p box (presets, absolute and relative bounds), `ranges_SI` through `PropertyPlot::isoline_range`, the values of lines placed through saturation states | `thermo/limits.py`, `thermo/provider.py` | 1 |
| `diagram/planning.{h,cpp}` | the curve requests of a spec: values (round values), points, labels, anchors, item ids | `engine/planning.py` | 1 |
| `diagram/engine.{h,cpp}` | the plot and curve caches, the compose record, the diagnostics | `engine/engine.py`, `engine/cache.py`, `corelib.compose_record` | 1 |
| `diagram/cycles.{h,cpp}` | the states of a solved model as full states, the paths between them as chains | `thermo/cycles/paths.py`, `thermo/process.py`, `core/cycle.py` | 2 |
| `diagram/humid_air.{h,cpp}` | the lines of the psychrometric and Mollier charts | `thermo/humid_air.py` | 3 |
| `diagram/references.{h,cpp}` | the badge lines: references of the equation of state, CoolProp and CoolSolve versions | `thermo/references.py` | 4 |
| `server.cpp` | the endpoints of §3.1 | | 1 |
| `gui/src/components/diagram/` | the tab, the SVG and its interaction layer, the state editor | | 1, 2 |

Used, not written: `CoolProp::Plot` (numerics) and coolplot-core (spec, units, themes,
compose, layout, scene, SVG, PDF). That replaces the whole porting table of the first
version: `core/`, `spec/`, `style/`, `compose/`, `scene.py`, `thermo/isolines.py`,
`thermo/tracer.py`, `render/plotly.py`, `render/ticks.py`.

Not used:

| CoolPlot part | Why |
|---|---|
| `engine/background.py`, `engine/store.py` | Threads and processes around Python's interpreter lock, and a disk cache. The server computes a scene in the request; the lines of a blend chart can be computed on several threads (one `PropertyPlot` each) if they are measured to be slow (phase 3). |
| `core/events.py`, `interaction.py` | Dragging cycle states edits a CoolPlot cycle; in CoolSolve the states come from the solve. Zoom, resize and hover are the GUI's (§3.2). |
| `render/` (matplotlib, plotly, SVG, web, PDF) | coolplot-core writes SVG and PDF (D1). |
| `thermo/cycles/models.py` | CoolPlot's cycle models compute a cycle; in CoolSolve the model is the cycle. |
| `chart.py`, `diagram.py`, `psychrometric.py` | Python facades; the GUI builds the spec directly. |

---

## 4. From a solved model to a cycle on the chart

CoolPlot draws a cycle from numbered states, chains of state numbers (`loops`) and paths
between them (`"auto"`: isobar if the pressures are equal, isenthalp, isentrope,
isotherm, otherwise a polytropic path). coolplot-core composes a cycle from states and
chains already computed (§1.2), so the paths are part of the glue (`diagram/cycles`).
The question here is how the GUI finds the states of a solved model.

### 4.1 Sources of states

1. **Arrays** (today's array overlay): `P[i]`, `h[i]` give state `i`; the loop is the
   index order. CoolPlot accepts cycles given as arrays.
2. **Inferred states**: variables used as arguments or results of the same property
   calls (`inferVariables`, `src/variable_inference.cpp`) are grouped into states (§4.2).
3. **Picked by the user** (`CS-FEAT-OVERLAY-STATES`): a state editor in the Diagrams
   tab where each state is a set of variables (`T_su_cp`, `P_su_cp`, …) and the loops are
   chains of states. This is the fallback when the automatic grouping is wrong, and the
   way to name states (CoolPlot's state labels, since 1.2, write the name next to the
   number).

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

States go into the spec in SI (CoolPlot's records are SI); the GUI converts the model
units (°C, Pa, J/kg) on the way. The compose record needs full states (every quantity of
the chart: T, p, h, s, ρ, v, u, Q): `diagram/cycles` keeps the solved values and
completes the others from two of them (the axes of the chart when the model has both)
with CoolSolve's fluid layer, so the marker sits where the solver put it. One chart
belongs to one fluid; a model with several fluids gets one chart per fluid, chosen in the
fluid list as today.

---

## 5. Phases

Each phase ends with the Diagrams tab usable and the old code it replaces removed. The
first version had a phase for labels; coolplot-core lays them out, so it is gone.

### Phase 0 — Preparation

- **CoolProp pin** (D3): move `COOLSOLVE_COOLPROP_PINNED_TAG` from `75af7816` to a
  `master` commit with `CoolProp::Plot` (`fc91a27f` or later), on `main` before the
  feature branch. It brings every CoolProp change since 2026-10-09 to the solver as well
  (for instance #3521, pure SRK/PR flashes through the superancillary, and #3528, critical
  points of mixtures), so the examples, the library regression and the robustness report
  are run again ([contributing.md](contributing.md#pinned-coolprop-version)).
- **coolplot-core as a dependency**: `FetchContent` of the CoolPlot repository at the
  reference commit with `SOURCE_SUBDIR cpp` and `COOLPLOT_CORE_TESTS=OFF` (no second
  Catch2), linked statically into `coolsolve`. Three points to settle on the way:
  - nlohmann/json (D4): coolplot-core vendors 3.12 and exports its `third_party`
    include directory publicly; CoolSolve fetches 3.11.3.
  - The Windows installer is built with MSVC; coolplot-core's CI builds it with GCC
    only. Build it on Windows first; a fix goes upstream.
  - zlib compresses the PDF streams when CMake finds it; without it the PDF is larger,
    not wrong.
- Settle D1 and D2.
- `tools/coolplot_fixtures.py`: from a list of chart specs, run CoolPlot at the
  reference commit and write, for each, the compose record (`coolplot.corelib.
  compose_record`: domain, requests, curves, cycle chains) and the scene to
  `tests/fixtures/coolplot/` (§6.1). Specs: R134a, water, CO2, R290 in every diagram
  type, T-h included; `DEF` and `ACHP` regions; every unit system; a family placed
  through saturation states; one cycle per diagram type.

### Phase 1 — Charts of pure fluids

- `diagram/fluid`, `domain`, `planning`, `engine` (§3.3); `POST /api/v1/diagram/scene`,
  `GET /api/v1/diagram/options`.
- Diagrams tab: fluid, diagram type (the seven CoolPlot types and T-h), units (`SI`,
  `KSI`, `engineering`, `IP`, and the model's units as a custom unit record, D7), region
  preset, line families (count, rounding, on/off, through saturation states), theme
  (`default` or `dark` following the GUI, `print`, `coolpack`), label mode (rails, edge,
  margin, inline).
- The SVG and its interaction layer (D1): rubber-band zoom and reset send a new view,
  a resize sends the new plot size, and the GUI sends back the label items of the scene
  on screen so that labels keep their places.
- Today's state and array overlays stay, as cycles with straight chains between the
  states, until phase 2 draws the paths.
- The chart spec is saved in the ZIP bundle and restored (gui.md §14.4 "Diagram state
  persistence"). Specs written by the GUI set `branding.references` to `"none"` until
  phase 4 computes the badge lines, so that an exported spec draws in CoolPlot as in the
  GUI.
- `exportPlots.ts` exports the chart for the LaTeX report from the SVG (rasterised in
  the browser) instead of `Plotly.toImage`; the new tab keeps the `thermo-diagram` class
  it looks for.
- Remove `POST /api/v1/coolprop/saturation` and `ThermoDiagram.tsx`.

### Phase 2 — Cycles from the model

- `diagram/cycles`: full states from the solved values (§4.3), paths as chains (isobars
  through the dome, isenthalps, isentropes, isotherms, polytropic paths), state numbers
  and labels, one legend entry; coolplot-core keeps the labels of the lines off the
  cycles.
- The three state sources of §4.1, the grouping rule of §4.2, the state editor.
- Hover: the variable names and values of a state marker (the markers of a cycle are in
  ascending state number); elsewhere, the full state under the cursor from
  `POST /api/v1/diagram/state`, debounced.
- Closes `CS-FEAT-OVERLAY-STATES`.

### Phase 3 — More fluids

- **Blends and CoolProp mixtures**: the port traces them and draws the dome from the
  phase envelope; the fluid string is enough. Measure the time of a blend chart (CoolPlot
  measured 5 s for the default p-h chart of R410A in Python, of which the port saves a
  tenth) and, if needed, compute the lines of one chart on several threads.
- **Incompressibles**: absolute limits (`PropertyPlot`'s absolute constructor).
- **Tabular backends** (`TTSE&HEOS`, `BICUBIC&HEOS`): the port draws the lines of
  constant quality with the equation of state behind the tables, as CoolPlot does.
- **Ideal gases** (`CS-FEAT-DIAGRAM-IDEAL`): a CoolSolve provider over `IdealGasFluid`
  (EES/JANAF offsets, enthalpy at the dummy pressure) writes the curve records of T-s,
  h-s and p-ρ charts over the gas region, with absolute limits; coolplot-core composes
  them as for any fluid.
- **Humid air** (`CS-FEAT-PSYCHRO`): coolplot-core composes the psychrometric (`WT`) and
  Mollier (`ThxW`) charts; their lines come from `diagram/humid_air` (or from
  `CoolProp::Plot`, D2); air processes from the solved states, with straight and
  `sensible` paths.
- Later, through the same provider interface: `NH3H2O` and LiBr-H2O charts (Dühring,
  h-x). Neither library knows such charts; they would be proposed upstream as new chart
  kinds, since coolplot-core has to know a diagram type to compose it.

### Phase 4 — Outputs

- The LaTeX report takes `thermo_diagrams.pdf` from coolplot-core's PDF writer instead
  of a PNG. It embeds TrueType fonts found on the platform; the online server needs a
  font package (CoolPlot's CI uses `fonts-dejavu-core`).
- SVG and PDF downloads in the tab, and *Export chart spec*: the chart spec with the
  model's cycle, as JSON that CoolPlot reads (`DiagramSpec.from_dict`), for matplotlib
  figures in Python.
- `CS-FEAT-PLOT-CLI` (`coolsolve --diagram spec.json model.eescode -o chart.svg`, `.pdf`,
  for README figures): a small task now, since the writers are in coolplot-core (the
  first version needed a port of `render/svg.py`).
- The badge: `diagram/references` computes the references of the equation of state from
  CoolProp's BibTeX keys (`get_fluid_param_string(fluid, "BibTeX-EOS")`) and its
  bibliography (`CoolPropBibTeXLibrary.bib` in the CoolProp sources, embedded at build
  time), and a versions line (CoolProp, CoolSolve); the GUI then writes CoolPlot's default
  (`"short"`).

---

## 6. Tests and keeping in step with CoolPlot

### 6.1 Fixtures

`tools/coolplot_fixtures.py` writes, for each spec, CoolPlot's compose record and scene to
`tests/fixtures/coolplot/`, with the CoolPlot commit, the Python CoolProp version (and
whether CoolPlot used the port or its Python numerics) and the date in a header file.
The fixtures are committed; the C++ tests (Catch2 `[diagram]`) read them without Python.

- **Compose records**: the curve requests (ids, keys, values, points, labels, anchors)
  exactly; the domain, the curves and the cycle chains within a relative tolerance,
  because the Python CoolProp wheel is not the pinned C++ commit; the same gaps. This is
  a comparison with another implementation, so the tolerance may be set from measured
  differences; no numerical tolerance of the glue itself is relaxed to make a fixture
  pass.
- **Scenes**: ids, groups, roles, z order, styles and texts as CoolPlot's. Label
  positions follow the curves, so they are compared within the tolerance of the curves.
  That coolplot-core composes what CoolPlot's Python composes, and that `CoolProp::Plot`
  computes what CoolPlot's Python numerics compute, is tested upstream; CoolSolve tests
  that it calls them right.
- No renderer fixtures with D1 (a): the GUI shows coolplot-core's SVG. The interaction
  layer's mapping from pointer to data coordinates gets unit tests if it grows (D5).

### 6.2 Invariants

Independent of CoolPlot, on every fixture fluid:

- the dome closes at the critical point; liquid and vapour lines meet there;
- an isoline crossing the dome contains the saturation states; inside it, an isobar of
  a T-s chart and an isotherm of a p-h chart are horizontal;
- no point is drawn where the flash failed;
- consistency with the solver: a state solved by CoolSolve lies on the isolines of its
  own properties, with the model's backend (checked for a `Real` fluid, a tabular
  backend and an ideal gas with a JANAF reference state).

### 6.3 Performance

Measure, in a Release build, the cold time of the default chart of each fixture fluid
and the warm time of a new scene (compose only; CoolPlot measured 12 ms for a p-h chart
with labels); set the budget from the measurements and add it to the test suite
([contributing.md §5](contributing.md)).

### 6.4 Moving to newer versions

- **CoolPlot**: bump the reference commit, regenerate the fixtures, adapt the glue to
  what the failing tests point to, update the header of this document and note it in
  [versions.md](versions.md). coolplot-core refuses spec fields it does not know, so the
  GUI offers a new CoolPlot feature only once the pin has it.
- **CoolProp**: the pin moves for the solver ([contributing.md](contributing.md#pinned-coolprop-version))
  and takes the numerics with it; the `[diagram]` tests run with the other reports.

---

## 7. Decisions

### 7.1 Open

| ID | Question | Proposal |
|---|---|---|
| D1 | How the GUI draws the scene: coolplot-core's SVG, CoolPlot's `scene_svg.js`, or a TypeScript port of the plotly backend? | coolplot-core's SVG with a TypeScript interaction layer (§3.2). |
| D2 | Where the glue lives. The plan (`engine/planning.py`) needs no CoolProp; cycle paths, humid-air lines and badge references do. | Upstream: `plan()` in coolplot-core (CoolPlot's review lists "line values with rounding" in it) and cycle paths and humid-air lines in `CoolProp::Plot` (CoolPlot's C0 plan had them). Write them in `src/diagram/` with CoolPlot's module boundaries when upstream is not ready at the phase that needs them, and drop the copy once it is. To agree with the CoolPlot maintainers. |
| D3 | CoolProp pin: a `master` commit with `CoolProp::Plot` now, or wait for a CoolProp release with it? | A `master` commit (the pin is one already: `v8.0.0` + 159 commits), moved on `main` ahead of the feature branch with the reports run again. |
| D4 | nlohmann/json: CoolSolve fetches 3.11.3, coolplot-core vendors 3.12 and puts it on the public include path. | Move CoolSolve to 3.12.0; ask upstream for an option to use the caller's `nlohmann_json` target. |
| D5 | A TypeScript test runner (e.g. Vitest, dev dependency only). | Only if the interaction layer grows beyond coordinate mapping. |
| D6 | Automatic grouping of states (§4.2) or user picking only. | Both: automatic by default, the editor to correct and name. Validate the rule on the library models first. |
| D7 | Default units of the Diagrams tab. | `engineering` (°C, bar, kJ/kg), as printed charts; the model's units (°C, Pa, J/kg) as a custom unit record of the spec (`{"name": ..., "units": {...}}`), which needs no upstream change. |

### 7.2 Settled since the first version

| First version | Now |
|---|---|
| D1: compose and layout in C++ or in TypeScript | In C++, by coolplot-core. |
| D2: C++ code in CoolSolve or in a library shared with CoolPlot | Shared: the numerics in CoolProp (`CoolProp::Plot`), the presentation in coolplot-core (`cpp/` of the CoolPlot repository; CoolPlot's plan, decision D7, gives it its own repository later). What is left is the glue (D2 above). |
| D3: T-h diagram | Upstream since CoolPlot 0.9: any pair of fluid axes (`Th`), in Python and in coolplot-core; `CoolProp::Plot` takes any pair CoolProp can flash (`HmassT`). Not tested upstream against the port, so it is a fixture of phase 0. |
| D4: units | `EUR` was renamed `engineering`; custom unit records exist (D7 above). |

---

## 8. Documentation to update when the work lands

- [gui.md](gui.md): §4.7 (endpoints), §6.8 (Diagrams tab), §8.3 (variable inference for
  diagrams), §14.4 (items done).
- [model_library_support.md §3](model_library_support.md): `CS-FEAT-OVERLAY-STATES`,
  `CS-FEAT-DIAGRAM-IDEAL`, `CS-FEAT-PSYCHRO`, `CS-FEAT-PLOT-CLI`.
- [contributing.md](contributing.md): file map, and the CoolPlot reference commit next to
  the pinned CoolProp.
- [versions.md](versions.md): the CoolPlot commit of each release.
- README: features, dependencies table (CoolPlot / coolplot-core), documentation table,
  future work.
