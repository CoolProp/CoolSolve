# CoolPlot Integration Plan

**Status:** implemented on branch `feat/coolplot-diagrams` (revised 2026-11-13; earlier
versions 2026-10-20 and 2026-10-31), to be reviewed before merging. Phases 0 to 4 are done
(§5); what remains is in §9.
**References:** [CoolPlot](https://github.com/CoolProp/CoolPlot) `v1.5.0` = `7338f48`
(2026-11-03), coolplot-core 1.5.0 from its `cpp/` folder; [CoolProp](https://github.com/CoolProp/CoolProp)
`8eb3b778` (2026-10-31), the commit CoolSolve is pinned to, which has `CoolProp::Plot`
(CoolProp PRs #3512 and #3519; unchanged on `master` since).

This document plans how the **Diagrams** tab of the GUI adopts CoolPlot, the
property-diagram package of the CoolProp project, instead of growing its own plotting
code. It records the decision (§1), the principles (§2), the architecture and the REST
contract (§3), how a solved model becomes a cycle on a chart (§4), the phases (§5), the
tests (§6) and the decisions (§7).

**What changed in this revision.** CoolSolve's CoolProp pin already has
`CoolProp::Plot` (`d9545e6`), so phase 0 no longer moves it. CoolPlot released 1.5.0 at
`7338f48`, 32 commits after the `16895a9` of the previous revision: coolplot-core 1.5
composes single lines, marked critical and triple points, a saturation scale and labels
inside the dome, whose states the caller computes (§1.2, glue table). coolplot-core is
built and linked (`5993f9d`). The REST contract and the diagram definitions are written
down (§3.1, §3.4), with an eye on HTML pages kept in model folders, which may later show
the diagrams of their model. The open decisions are settled (§7).

---

## 1. Context and decision

### 1.1 The Diagrams tab before this work

`gui/src/components/ThermoDiagram.tsx` (452 lines, Plotly.js) draws:

- the saturation dome of one `Real` fluid in T-s, P-h, h-s or T-h, from
  `POST /api/v1/coolprop/saturation` (`src/server.cpp`: `PropsSI` with the fluid's
  CoolProp name, so the default backend whatever the model's `coolpropBackend`; 200
  temperatures log-spaced in T_crit − T between T_triple + 1 K and T_crit − 0.1 K; arrays
  in SI with T in °C);
- the solved state points, grouped by a name heuristic (`T_1`, `P_1`, `h_1` → state `1`)
  from `GET /api/v1/variables/inferred`;
- one array path (`X[i]`, `Y[i]`), optionally closed.

There are no isolines, no unit choice (°C, Pa, J/kg only), no labels, and only `Real`
fluids are offered. The wanted extensions are listed in [gui.md §14.4](gui.md) (isolines,
units, hover, persistence) and in the feature register of
[model_library_support.md §3](model_library_support.md): `CS-FEAT-DIAGRAM-IDEAL` (T-s,
P-v, h-s of ideal gases), `CS-FEAT-PSYCHRO` (psychrometric chart, with PNG/SVG export),
`CS-FEAT-OVERLAY-STATES`, `CS-FEAT-PLOT-CLI` (diagrams and parametric plots from the
command line).

### 1.2 CoolPlot 1.5

CoolPlot went from 0.4.0.dev0 (`22937c8`) to 1.5.0 (`7338f48`) in 226 commits. It draws
log(p)-h, T-s, h-s, p-s, p-ρ, T-ρ, p-T and any other pair of fluid axes (`Th`), of pure
fluids, blends and incompressibles, with cycles on top, and psychrometric charts (ASHRAE
and Mollier h,x) with air processes. The pipeline is `spec → plan → curves → compose →
layout → scene → renderer`, with its records in CoolPlot's `docs/reference/contracts.md`
(spec, curve and scene records; the compose record is documented in
`cpp/include/coolplot/compose.h` only). Two parts of it exist in C++.

**Numerics in CoolProp: `CoolProp::Plot`** (`include/CoolProp/plotting/CoolPropPlot.h`,
`src/CoolPropPlot.cpp`, 2,780 lines with their Catch2 tests), a port of CoolPlot's
`thermo/isolines.py`, `tracer.py` and `stepping.py`, compiled into CoolProp's static
library by its `src/*.cpp` glob.

- Lines are stepped with steps chosen from their shape on the page; nothing is computed
  outside the axes; isotherms and isobars are followed through the dome in vapour
  quality; other lines get the saturation state inserted where they cross the dome;
  mixtures are traced (Newton in T and ρ from the previous point) and their dome comes
  from the phase envelope; a failed point is a NaN gap, never interpolated.
- `PropertyPlot(fluid_string, ykey, xkey, TPLimits)` with the presets `None`, `Def`,
  `Achp`, `Orc`, or `PropertyPlot(fluid_string, ykey, xkey, T_min, T_max, p_min, p_max)`
  (absolute box, ±inf for the limit of the equation of state, clamped to it; an
  incompressible needs one). The plot is move-only.
- `calc_isolines(key, values, CalcSettings)` returns one `Isoline` per value: `x`, `y` in SI
  (NaN where a point failed), `value`, and the diagnostics of CoolPlot's curve record
  (`n_failed`, `n_fallback`, `n_inserted`, `n_evaluated`, `n_outside`, `method`,
  `reason`). `CalcSettings`: `points`, `saturation_points`, `solver` (`auto`, `tracer`,
  `flash`), `critical_margin`, `saturation_spacing`. A line that cannot be computed comes
  back with its `reason`; an axis key, a Q outside [0, 1] or a Q line without critical
  point throws.
- Also `isoline_range(key)` (from the four corners of the box), `tp_box()`,
  `supported_isoline_keys()`, `sweep_key(key)`, `is_mixture()` and
  `value_at(key, x, y)` (the state under the cursor, NaN on failure).
- A plot builds its own `AbstractState` from the fluid string, so it takes the backend
  (`HEOS`, `INCOMP`, …) and the global configuration of the CoolProp it is linked with.
  Its calls hold a per-plot lock: lines are computed in parallel with one plot per thread.
  It has no tabular handling: CoolPlot's Python draws the lines of constant quality of a
  `TTSE&HEOS` or `BICUBIC&HEOS` fluid with the base backend itself.
- CoolPlot's Python (`thermo/native.py`) sends its fluid lines there when the installed
  CoolProp has the port, always through the absolute constructor with a box it resolved
  itself, and checks it point for point against its Python numerics
  (`tests/test_native.py`; these tests skip in CoolPlot's CI, which installs CoolProp
  8.0.0, without the port). Humid-air lines and cycle states are not part of the port.

**Presentation in C++: coolplot-core 1.5.0** (`cpp/` of the CoolPlot repository, 7,735
lines of C++17 plus nlohmann/json 3.12.0, vendored and on its public include path; it uses
neither CoolProp nor a plotting library).

- Ports `coolplot.core` (number format, quantities, unit systems `SI`, `KSI`,
  `engineering`, `IP` and custom unit records, diagram types and their families,
  `snap`), `coolplot.style` (themes `default`, `dark`, `print`, `coolpack`),
  `coolplot.scene`, `coolplot.compose` (builder, layout, rails, the automatic legend and
  badge corners) and the SVG writer, and adds a PDF writer with TrueType subsets.
- Interfaces: C++ (`coolplot::compose`, `compose_json`, `scene_to_svg`,
  `svg_plot_size_pt`, `scene_to_pdf`, `pdf_plot_size_pt`), a C ABI of JSON text in and
  JSON, SVG or PDF out (`coolplot_compose`, `coolplot_scene_to_svg`,
  `coolplot_scene_to_pdf`, `coolplot_svg_plot_size_pt`, `coolplot_version`,
  `coolplot_free`), and a nanobind module for Python. `coolplot::json` is
  `nlohmann::json`: with the same version on both sides, CoolSolve hands its records over
  in memory.
- Its input is the **compose record**: `spec` (the chart spec record); `domain`
  (`ranges_SI` of the quantities, of which the composer uses the two axes; `tp_box`, read
  with `points` only; since 1.5 `points`, the states of the critical and triple points the
  spec marks, `scale`, the ticks of the saturation scale, and `dome`, the bubble and dew
  states of the lines labelled inside the dome); `curves` (each with its `request`:
  `item_id`, `key`, `value_SI`, `points`, `labels`, `anchor_SI`, and the curve record, of
  which only `x` and `y` are read; single lines come as `line/<name>`); `cycles` (resolved
  states and chains of states, paired by index with the spec's cycles);
  `plot_size_pt`; `previous_labels` (the isoline labels of the scene on screen, so that
  they keep their places); `sources` (the badge lines).
- Unknown fields of a spec are ignored, unknown values (a theme, a unit system, a label
  mode, a quantity) refused, as are logos and watermarks. A caller therefore offers a
  CoolPlot feature only once its pinned core has it (CoolPlot's Python bridge gates them
  in `corelib._NEEDS_A_NEWER_CORE`).
- Scenes and SVG documents are byte for byte those of the Python code (golden files in
  `cpp/tests/data`, written by `make_data.py`). Measured by CoolPlot: 12 ms to compose a
  p-h chart with labels, JSON included. CoolPlot's CI builds it on Linux with GCC only.
  Compose and the SVG writer keep no mutable state, so requests may run them
  concurrently; the PDF font cache holds a lock.

Features since 0.4 that the GUI gets with the core: labels on rails and a legend placed
inside the plot area by itself, lines placed through saturation states (isentropes every
10 K of saturation temperature, `anchor`), state labels on cycles, symbols written with
their own characters (ρ, °C), the `IP` and `engineering` units, one colour per quantity
checked for colour vision deficiencies, a badge naming the equation of state and the
versions (on by default in CoolPlot since 1.5), single lines (`lines`: a highlighted
isotherm, the pressure of a burst disk), the critical and triple points marked with a dot
(`marked_points`), a saturation scale next to the vertical axis (`saturation_scale`), and
labels inside the dome (`settings.labels.dome`).

**Not in either library** (the "glue"; Python line counts):

| Step | CoolPlot (Python) | Lines | Needs CoolProp |
|---|---|---|---|
| T-p box from the spec's limits (presets, absolute and relative bounds), ranges of the quantities | `thermo/limits.py` (`resolve_tp_box`, `property_range_SI`), `thermo/provider.py` (`chart_domain`) | ≈80 | yes; mostly `PropertyPlot` (box, `isoline_range`) |
| Values of the lines placed through saturation states | `thermo/provider.py` (`saturation_anchor_values`), `engine/engine.py` | ≈70 | yes |
| Plan: which lines at which values (rounding), points, item ids, the saturation curve, single lines | `engine/planning.py` | 254 | no; coolplot-core exports `snap`, `UnitSystem`, `DiagramType`, `Family` |
| Curves: one request each, cache, lines with fewer than two valid points dropped with a warning | `thermo/native.py`, `engine/engine.py`, `engine/cache.py` | ≈190 | yes (the port) |
| Marked points, saturation scale, dome ends (1.5) | `provider.marked_point_states`, `saturation_pairs`, `dome_ends`; `Engine._scale_ticks`, `_dome_ends` | ≈150 | yes |
| Compose record, plot size in two passes (recompose when the plot area of the SVG differs by more than 2 %) | `corelib.compose_record`, `chart.py` | ≈40 | no |
| Cycle paths between full states | `thermo/cycles/paths.py` (`auto_kind`, `path_states`, `cycle_paths`) | ≈70 | yes |
| Humid-air lines | `thermo/humid_air.py` | 569 | yes (`HAPropsSI`) |
| Badge lines: references of the equation of state, versions | `thermo/references.py`, `Engine.source_lines` | ≈440 | yes (`get_BibTeXKey(fluid, "EOS")` and `CoolPropBibTeXLibrary.bib`) |

### 1.3 Options considered

| Option | Verdict |
|---|---|
| **A. Run CoolPlot's Python code** (embedded CPython, a Python process next to the server, or Pyodide in the browser) | Rejected. The Windows build is one self-contained `.exe`; Python, numpy and the CoolProp wheel would have to ship with it. The curves would come from a second CoolProp, not the pinned one, and without CoolSolve's reference-state offsets (ideal gases are on the EES/JANAF basis, `IdealGasFluid::getReferenceState`): solved states would sit beside the lines. CoolSolve's own fluids (`NH3H2O`, LiBr) could never be drawn. |
| **B. Keep developing our own plots** | Rejected. Every item of §1.1 has been designed, measured and reviewed in CoolPlot. |
| **C. Port CoolPlot to C++ and TypeScript inside CoolSolve** (first version of this plan) | Superseded: CoolPlot made both C++ ports itself, tested against its Python code. |
| **D. Link `CoolProp::Plot` and coolplot-core, write the glue** | **Chosen.** The numerics come with the pinned CoolProp, the same library and configuration the solver uses; coolplot-core is built from source like the other dependencies; CoolSolve writes the glue of §1.2 and the Diagrams tab. |

---

## 2. Principles

1. **CoolPlot is the reference.** Contracts, item ids, groups, roles, styles, texts and
   numerical methods are CoolPlot's. CoolSolve follows two pinned versions: CoolProp
   (the numerics, pinned for the solver) and CoolPlot (coolplot-core and the reference
   for the glue, `COOLSOLVE_COOLPLOT_PINNED_TAG`). Moving either is a deliberate change
   (§6.4).
2. **Properties come from CoolSolve's fluid layer**, never from a second property
   library. `CoolProp::Plot` is the pinned CoolProp linked into `coolsolve`, with the
   process's configuration (superancillaries) and CoolProp's reference states, which are
   those of the solver for `Real` fluids (`Fluid::getReferenceState` gives no offset).
   The fluid string comes from `FluidRegistry`, the backend from the session's solver
   option `coolpropBackend` (`coolpropConfig.backend`). Fluids that CoolProp alone does
   not describe (ideal gases on the EES/JANAF basis, `NH3H2O`, LiBr-H2O) get CoolSolve
   providers that write curve records themselves. A state point of a solved model lies
   exactly on the lines of its own properties.
3. **The contracts are CoolPlot's, unchanged.** The GUI sends a chart spec in CoolPlot's
   format; the server builds a compose record and coolplot-core returns the scene. A
   feature CoolPlot lacks is proposed upstream first; a CoolSolve-only extension is a
   documented, namespaced addition next to the spec (the cycle definitions of §4), never
   a silent change of meaning.
4. **Use the libraries, do not copy them.** `CoolProp::Plot` and coolplot-core are linked
   as they are, never patched in CoolSolve: a defect or a missing feature is fixed
   upstream and reaches CoolSolve by moving a pin. The glue mirrors CoolPlot's Python
   modules (names, boundaries, order of the steps), names its source file and the
   CoolPlot commit, and keeps CoolPlot's MIT notice ("Copyright 2014-2019 Jorrit Wronski
   and the other developers"), so that it can move upstream (D2).
5. **No Python at run time or at test time.** Python is needed only to regenerate the
   test fixtures with a script in `tools/` (§6.1); the fixtures are committed.
6. **Not solver core.** The work touches the server, a new `diagram` module and the GUI,
   not the parser, evaluator or solvers; it is developed on `feat/coolplot-diagrams` and
   reviewed before merging.
7. **A diagram is a document of the model, drawn by the server.** A diagram definition
   (the CoolPlot spec, the fluid and the cycles as variable names, §3.4) is plain JSON
   kept with the model, and the server turns it into SVG with the session's solution in
   one request. The Diagrams tab is one client of that request; an HTML page of a model
   folder can later be another, without new code on the server side.

---

## 3. Architecture

### 3.1 Data flow and REST contract

```
GUI (Diagrams tab)                     coolsolve server (C++)
------------------                     ------------------------------------------------------------
diagram definition (§3.4)    POST      spec --> domain, plan --> curve requests
+ figure size in pt          ------->          curves <-- CoolProp::Plot (pinned CoolProp), cached
+ labels of the shown scene  /api/v1/          cycles: variables --> full states, paths <-- fluid layer
                             diagram/          marked points, scale, dome ends, badge lines
                             scene     = compose record --coolplot-core--> scene --> SVG
SVG + interaction layer  <-- SVG, axes, labels, states, diagnostics <-------------
(zoom, hover, state names)
```

**`POST /api/v1/diagram/scene`** — a chart of the session's model.

```jsonc
// request
{
  "fluid": "R134a",              // a CoolSolve fluid name (FluidRegistry), as in the model
  "spec": { ... },               // CoolPlot spec record (DiagramSpec.to_dict), without "fluid"
                                 // and "cycles": the server writes them
  "cycles": [ ... ],             // cycle definitions, states as variable names (§4.1);
                                 // resolved with the session's last solve, skipped without one
  "figure_size_pt": [w, h],      // the whole figure; the server finds the plot area (two passes)
  "previous_labels": [ ... ]     // optional: "labels" of the previous response
}
// response
{
  "svg": "<svg ...>",            // coolplot-core's document
  "axes": { ... },               // x_key, y_key, x_label, y_label, x_limits, y_limits (display
                                 // units), x_log, y_log, x_unit, y_unit ({label, scale, offset}:
                                 // display = SI * scale + offset), plot_size_pt
  "labels": [ ... ],             // isoline label items, to send back as previous_labels
  "states": [ ... ],             // per cycle state: number, label, variables, values in SI
  "spec": { ... },               // the spec as composed (fluid and cycles filled in): CoolPlot
                                 // draws it (Export chart spec)
  "diagnostics": { "curves": [ ... ], "warnings": [ ... ] },   // per curve: item_id, drawn, the
                                 // counts of the curve record, method, reason
  "timing_ms": { "curves": ..., "compose": ..., "svg": ... }
}
```

An unknown fluid, a fluid kind the server cannot draw yet, an unknown diagram type or a
spec value coolplot-core refuses is a 400 with `{"error": ...}`; a line that cannot be
computed is not an error but a diagnostic. An incompressible solution takes its
concentration as CoolPlot's fluid record does: `"fractions": [0.3]` (mass fractions) and
`"basis": "mass"` next to `fluid` (400 when it is missing, out of the fluid's range, or given
to a fluid that takes none). For a fluid without saturation curve (incompressibles, ideal
gases) the composed `spec.limits` is the absolute T-p box of its lines, as CoolPlot needs.
Humid air (`AirH2O`, charts `WT` and `ThxW`) takes the total pressure of the chart as
`"p_Pa"` (default 101325 Pa) and a CoolPlot `Region` as limits (`{"bounds": {"T": [...],
"W": [...]}}`, by default PsychrometricChart's −10 to 50 °C and 0 to 30 g/kg).

**`POST /api/v1/diagram/export`** (phase 4) — the chart as a file: the body of
`/diagram/scene` and `"format": "svg"` or `"pdf"`; the file with
`Content-Disposition: attachment; filename="<model>_<diagram>.<format>"`. The figure size is
the request's, else the spec's page (576 × 432 pt). The PDF embeds the TrueType fonts the
server finds (§5, phase 4); without any, the route answers 500 and SVG still works.

**`GET /api/v1/diagram/options`** — what the GUI offers: diagram types (name, short title
from the axis symbols as the chart writes them, vertical first, `p-h`, `T-s`, `p-ρ`, and
`psychrometric`, `Mollier h-x`; axes, families with their default count), the units of the
model (`model_units`: `record`, the custom unit record of CoolSolve's units that goes into a
spec's `units`, and `units`, its unit table; every chart of the GUI is in these units, D7),
region presets (`DEF`, `ACHP`, `ORC`), themes, label modes,
the coolplot-core version and CoolPlot commit, the fluid types drawn and, per type,
`fluid_kinds`: the diagram types and families it offers, whether it has a saturation curve,
whether it takes a concentration, its default T-p box, and the region of a new chart
(`DEF`; `ACHP` for mixtures, whose `DEF` reaches the triple point of their lightest
component, where their lines are slow to compute; the default `Region` of humid air, with
its total pressure and the temperature range of `HAPropsSI`).

**`GET /api/v1/diagram/definitions`**, **`PUT /api/v1/diagram/definitions`** — the
diagram definitions of the session (§3.4): stored with the model, written to and read
from the ZIP bundle.

**`POST /api/v1/diagram/state`** (phase 2) — the full state at a point of a chart, for the
hover: `{fluid, diagram, x_SI, y_SI}` (and `p_Pa` for humid air, whose state is
`{T, W, RH, Twb, Tdp, h, v, Thx, p}`) → `{"state": {T, p, h, s, rho, u, v, Q}}` in SI (Q = −1 in
single phase, a key that cannot be computed is null) or `{"state": null}` where the point has no
state. One flash of the two axis values on the cached fluid model (9 µs in process, 0.5 ms over
HTTP), with the model's backend.

**`GET /api/v1/diagram/states`** (phase 2) — the states of the solved model grouped by the
property calls they share (§4.2): `{"fluids": {"R134a": [{"number": 1, "label": null, "vars":
{"T": "T_1", "p": "P_1", "h": "h_1", "s": "s_1"}}, ...]}}`, empty without a successful solve.
Computed at solve time next to the inferred variables (median 0.013 ms per library model), kept in
the snapshot of *Back*.

`GET /api/v1/coolprop/fluids` and `GET /api/v1/variables/inferred` stay.
`POST /api/v1/coolprop/saturation` is removed with `ThermoDiagram.tsx`.

**Caches**, process-wide (the server has several sessions), bounded and thread-safe:
`PropertyPlot` objects, which hold the critical point and the phase envelope of a fluid,
keyed by fluid string, axes and T-p box; and curves, keyed as in CoolPlot by fluid
string, CoolProp git revision and configuration (backend, superancillaries), diagram,
family, value and calculation settings. Changing units, theme, labels, view or figure
size computes no property: only the compose step runs again.

### 3.2 How the GUI draws the scene (decision D1)

coolplot-core composes the scene and lays out its labels for the plot size and the
view, so every zoom and resize asks the server for a new scene in any case. The browser
shows **coolplot-core's SVG** as it comes: each item is a `<g data-item="<id>">` with its
role and group as classes, the plot area is the `rect.cp-background`, and the root
carries the axis limits and log flags. A thin TypeScript layer maps the pointer to data
coordinates, draws the rubber band of the zoom (a new `view_SI`), resets it, and finds the
state markers for the hover; CoolPlot's `render/web.py` (pointer client) and
`render/static/scene_svg.js` are its references. The same writer gives the CLI and the
report their figures, with the PDF writer next to it. Rejected: CoolPlot's `scene_svg.js`
in the browser (the round trip stays for the labels) and a TypeScript port of
`render/plotly.py` (no reference to test against; labels wrong after a client-side zoom).

### 3.3 Module map

New C++ module `include/coolsolve/diagram/`, `src/diagram/`; new GUI folder
`gui/src/components/diagram/`. Each glue module mirrors the CoolPlot module in the
third column.

| CoolSolve | Does | Mirrors (CoolPlot) | Phase |
|---|---|---|---|
| `diagram/fluid.{h,cpp}` | the CoolProp fluid string of a model's fluid (`FluidRegistry` name, backend, mixture fractions); which fluids the port takes, which need a CoolSolve provider | `thermo/native.py` (`fluid_string`), `spec/fluid.py` | 1 |
| `diagram/domain.{h,cpp}` | the T-p box (presets, absolute and relative bounds), `ranges_SI`, the values of lines placed through saturation states, marked points, saturation scale, dome ends | `thermo/limits.py`, `thermo/provider.py` | 1 |
| `diagram/planning.{h,cpp}` | the curve requests of a spec: values (rounding), points, labels, anchors, item ids, the saturation curve, single lines | `engine/planning.py` | 1 |
| `diagram/engine.{h,cpp}` | the plot and curve caches, the compose record, the two-pass plot size, the diagnostics, the definitions document | `engine/engine.py`, `engine/cache.py`, `corelib.compose_record` | 1 |
| `diagram/cycles.{h,cpp}` | the states of a solved model as full states, the paths between them as chains | `thermo/cycles/paths.py`, `thermo/process.py`, `core/cycle.py` | 2 |
| `diagram/states.{h,cpp}` | the grouping of the variables of property calls into states (§4.2) | (CoolSolve only) | 2 |
| `diagram/ideal_gas.{h,cpp}` | ideal gases as a `PropertyProvider` (`fluid.h`, the interface of CoolPlot's providers) on the solver's ideal-gas functions with the EES/JANAF offsets, and their lines | `thermo/provider.py` (interface); the lines are CoolSolve only | 3 |
| `diagram/humid_air.{h,cpp}` | the lines of the psychrometric and Mollier charts | `thermo/humid_air.py` | 3 |
| `diagram/references.{h,cpp}` | the badge lines: references of the equation of state (CoolProp's bibliography, embedded by `cmake/embed_bibtex.cmake`) and the versions line, as CoolPlot writes them; the SVG and PDF files of a chart (export, report, CLI) are in `engine` | `thermo/references.py`, `Engine.source_lines` | 4 |
| `server.cpp` | the endpoints of §3.1 (thin handlers), `diagrams.json` in the ZIP bundle | | 1 |
| `gui/src/components/diagram/` | the tab, the SVG and its interaction layer, the overlays | | 1, 2 |

Used, not written: `CoolProp::Plot` (numerics) and coolplot-core (spec, units, themes,
compose, layout, scene, SVG, PDF).

Not used:

| CoolPlot part | Why |
|---|---|
| `engine/background.py`, `engine/store.py` | Threads and processes around Python's interpreter lock, and a disk cache. The server computes a scene in the request; the lines not yet cached are computed on up to 4 threads (`MAX_CURVE_THREADS`, never more than the machine has), one cached `PropertyPlot` each, which gives the same lines bit for bit. |
| `core/events.py`, `interaction.py` | Dragging cycle states edits a CoolPlot cycle; in CoolSolve the states come from the solve. Zoom, resize and hover are the GUI's (§3.2). |
| `render/` (matplotlib, plotly, SVG, web, PDF) | coolplot-core writes SVG and PDF. |
| `thermo/cycles/models.py`, `resolve_states` | CoolPlot's cycle models compute a cycle; in CoolSolve the model is the cycle and its states are solved. |
| `chart.py`, `diagram.py`, `psychrometric.py` | Python facades; the GUI builds the spec directly. |

### 3.4 Diagram definitions (`diagrams.json`)

The diagrams of a model are one JSON document, `diagrams.json`:

```jsonc
{
  "format": "coolsolve-diagrams",
  "version": 1,
  "diagrams": [
    {
      "name": "Cycle on log(p)-h",    // unique in the document: how a page refers to it
      "fluid": "R134a",               // with "fractions" and "basis" for a solution, "p_Pa"
                                      // for humid air (§3.1)
      "spec": { ... },                // CoolPlot spec record (no "fluid", no "cycles")
      "cycles": [ ... ]               // §4.1
    }
  ]
}
```

- It travels in the ZIP bundle next to `parametric_studies.json`, and the GUI reads it back
  through `showSessionModel()` (opening a bundle, *Back*, reloading a library tab). The
  Diagrams tab edits one diagram of the list (phase 1: the first; a selector later).
- The `.eescode` stays EES code: EES keeps its plots in the GUI too.
- **Later, in model folders.** The same file can sit in a model folder of the library
  (next to the `.eescode`, its `.initials` and its tables), so that opening a model shows
  the chart its author made; and an HTML page of the folder can ask the server for the
  SVG of a named diagram of the solved model and show it among its overlaid variables.
  That needs the library's allow-list of embedded files (`cmake/embed_library.cmake` and
  `EMBED_FILES` of the library's `tools/build_index.py`) and is not part of this plan;
  the document and the scene request are made so that it needs no other change.

---

## 4. From a solved model to a cycle on the chart

CoolPlot draws a cycle from numbered states, chains of state numbers (`loops`) and paths
between them (`"auto"`: isobar if the pressures are equal, isenthalp, isentrope,
isotherm, otherwise a polytropic path). coolplot-core composes a cycle from states and
chains already computed (§1.2), so the paths are part of the glue (`diagram/cycles`).

### 4.1 Cycle definitions

A cycle of a diagram definition names the model's variables; the server reads their
values from the session's last solve:

```jsonc
{
  "name": "Heat pump",
  "numbers": true,                    // CoolPlot's CycleSpec fields pass through
  "states": [                         // explicit states...
    { "number": 1, "label": "su,cp", "vars": { "p": "P_su_cp", "h": "h_su_cp", "T": "T_su_cp" } },
    { "number": 2, "vars": { "p": "P_ex_cp", "h": "h_ex_cp" } }
  ],
  "arrays": { "p": "P", "h": "h" },   // ...or arrays: state i from P[i], h[i]
  "loops": [[1, 2, 3, 4, 1]],         // chains of state numbers; default: the states in order
  "closed": true,                     // with no "loops": close the loop of the states (default)
  "paths": "auto",                    // one kind for every connection, or CoolPlot's record
                                      // [[[1, 2], "isentropic"], ...] (unlisted connections: auto)
  "steps": 20                         // points computed along each path
}
```

- `paths`: `auto` (the default: isobaric if the pressures are equal, else isenthalpic,
  isentropic, isothermal, otherwise polytropic, as CoolPlot's `auto_kind` with its relative
  tolerance 1e-6), `straight`, `isobaric` (through the dome, with both saturation states
  inserted), `isenthalpic`, `isentropic`, `isothermal`, `polytropic`. An unknown kind is a 400.
  Points that cannot be computed are left out, as CoolPlot does. The exported spec carries
  CoolPlot's record (the connections that are not `auto`) and `steps`.

- No `loops`: CoolPlot's default, the states in order, closed (unless `"closed": false`).
- `"loops": []`: markers and numbers only, no lines (the inferred states). In the exported
  spec this becomes chains of one state each, because an empty `loops` means the default loop
  to CoolPlot.
- `"states": []` (an overlay switched on before any solve) is stored and skipped when drawing.
- A missing variable, a state that cannot be completed or no successful solve leaves the
  state or the cycle out with a warning; it is never an error. A cycle drawn with lines
  needs all its states; a cycle of markers only (`"loops": []`) draws the states that can be
  computed and names the others in a warning (the name heuristic groups a saturation T and p,
  `T_ev` and `P_ev`, which fix no state).
- State keys are CoolPlot's (`T`, `p`, `h`, `s`, `rho`, `u`, `v`, `Q`); the GUI maps the codes
  of the variable inference (`P` → `p`, `D` → `rho`, `VOL` → `v`, …).
- The Diagrams tab writes two kinds of cycles: the inferred states, `{"name": "Inferred states",
  "numbers": true, "states": [...], "loops": []}` (markers only), and any number of arrays
  overlays, `{"name": "T, s", "numbers": true, "arrays": {"s": "s", "T": "T"}, "closed": false,
  "paths": "auto"}`: the arrays of the chart's horizontal and vertical quantities (on the Mollier
  chart a temperature array is keyed `T`), joined in index order by the server, named after their
  arrays (unique in the definition). Other cycles of a stored definition (those of the former state
  editor) are not drawn by the tab and are dropped at its next edit; the server, the report and
  `coolsolve --diagram` still draw any cycle of this section.

Sources of states:

1. **Arrays** (the former array overlay): `P[i]`, `h[i]` give state `i`; the loop is the
   index order, closed on request. The tab's overlays pick the arrays of the chart's axis
   quantities.
2. **Inferred states**: variables used as arguments or results of the same property
   calls (`inferVariables`, `src/variable_inference.cpp`) are grouped into states (§4.2)
   and drawn as markers.
3. **Picked by the user** (`CS-FEAT-OVERLAY-STATES`): explicit states of any variables and
   chains of states, as in the example above. The Diagrams tab had a state editor for them
   (phase 2); it was removed on 2026-11-13 at the maintainer's request, to keep the tab simple:
   the states come from the grouping and from arrays. Such cycles can still be written in
   `diagrams.json` for the report and `coolsolve --diagram`.

### 4.2 Grouping inferred variables into states (decision D6)

The former name heuristic merges everything that shares a suffix. A rule based on the
equations: each property call on a fluid involves a set of state variables (inputs and
output); two calls describe the same state when they share **two variables of
independent properties** (state postulate for a pure fluid). `h_2 = enthalpy(R134a,
P=P_2, s=s_2)` and `T_2 = temperature(R134a, P=P_2, h=h_2)` share `P_2` and `h_2`: one
state. `T_sat = temperature(R134a, P=P_cond, x=1)` shares only `P_cond` with the
condenser states: a separate state. T and P are not independent inside the dome, so a
pair (T, P) alone does not merge two calls — **except** when neither call involves a vapour
quality nor is a saturation function (`T_sat`, `P_sat`): a call that takes no quality and has T
and p among its variables fixes a single-phase state, where T and p are independent.
`h_1 = enthalpy(R134a, T=T_1, P=P_1)` and `s_1 = entropy(R134a, T=T_1, P=P_1)` are one state.

As implemented (`src/diagram/states.cpp`): the property functions `temperature`, `pressure`,
`enthalpy`, `entropy`, `intenergy`, `volume`, `density`, `quality`, `T_sat`, `P_sat` (and their
short names; `v(...)` is the viscosity) on a known fluid, with their named-argument variables
and the variable of `var = call`; a fluid given as a string variable is resolved when the IR
knows its value. A state needs two independent properties; for an ideal gas, two of T, h
and u are not independent. Numbering from the names: the number that more than half of the
variables of a state carry (`X[i]`, or the integer ending the name: `T_2`, `h_2`, `s[2]` → 2),
the states whose names agree most claiming their number first; the others are labelled by
their common suffix (`su_cp`) or by their number when it is taken, and take the next free
numbers.

Humid air: calls on `AirH2O` are one state when they share two independent humid-air
variables other than the total pressure (the humidity ratio and the dew point are not
independent); two calls at the same T and p never merge.

**Validation** on the 287 runnable models of the library (hidden test
`[diagram-states-library]`, table in `build/diagram_states_library.md`): the name heuristic
finds 861 states, the strict rule 2683, the rule with the (T, p) exception 2315 (2836 and 2468
before ideal gases counted T, h and u as dependent). The 368 states the exception removes are
states evaluated property by property at (T, P), split by the strict rule. Where they differ, the rule completes states through shared variables the names miss
(`P_atm`, `p_steam`, `p[2]`) and finds the states of arrays, which the heuristic ignored; the
heuristic mixes states (a dead state `T_0`, `p_0` reused by several calls) and keeps pairs
(T_sat, P_cd) that fix no state. The rule replaces the heuristic in the Diagrams tab.

### 4.3 Values

The variables hold the model's units (°C, Pa, J/kg, J/(kg·K), kg/m³); the server
converts them to SI with CoolSolve's `UnitConverter`. The compose record needs full
states (T, p, h, s, ρ, v, u, Q): `diagram/cycles` keeps the solved values and completes
the others from two of them (the axes of the chart when the model has both, else an
independent pair) with CoolSolve's fluid layer (CoolProp with the model's backend, the
ideal-gas provider for ideal gases), so the marker sits where the solver put it. A state of
an incompressible given without a pressure is placed at 101325 Pa, the pressure the solver
uses for its T-only calls. One chart belongs to one fluid; a model with several fluids gets
one diagram per fluid.

---

## 5. Phases

Each phase ends with the Diagrams tab usable and the old code it replaces removed.

### Phase 0 — Preparation (done)

- CoolProp pin with `CoolProp::Plot`: `8eb3b778` since `d9545e6` (2026-11-03), on `main`.
- coolplot-core as a dependency (`5993f9d`): `FetchContent_Populate` of the CoolPlot
  repository at `COOLSOLVE_COOLPLOT_PINNED_TAG` (the `v1.5.0` commit) and
  `add_subdirectory(cpp EXCLUDE_FROM_ALL)`, tests and Python module off, linked into
  `coolsolve_lib`; `-DFETCHCONTENT_SOURCE_DIR_COOLPLOT=<checkout>` builds a local
  checkout. nlohmann/json moved to 3.12.0, the version coolplot-core vendors (D4). The
  zlib of the platform compresses the PDF streams when CMake finds it (Linux); without
  it the PDF is larger, not wrong. coolplot-core needs CMake 3.16.
- Not verified: the MSVC build of coolplot-core (CoolPlot's CI builds it with GCC only).
  A fix goes upstream.

### Phase 1 — Charts of pure fluids (done)

Done on `feat/coolplot-diagrams`: `src/diagram/` (`fluid`, `domain`, `planning`, `engine`),
the endpoints but the state, the new tab (`gui/src/components/diagram/`), `diagrams.json` in
the bundle, the report figure from the SVG, `ThermoDiagram.tsx` and the saturation route
removed, `coolsolve --version` naming coolplot-core and its CoolPlot commit. CoolProp mixtures
are drawn too (the fluid string is enough; pseudo-pure blends get a warning). Lines of constant
quality of a tabular backend use the base backend (wired, not tested yet). What the list below
planned:

- `diagram/fluid`, `domain`, `planning`, `engine` (§3.3), with the marked points, the
  saturation scale and the dome ends of coolplot-core 1.5; the endpoints of §3.1 but the
  state.
- Diagrams tab: fluid, diagram type (the fluid types of coolplot-core, T-h included),
  units (`SI`, `KSI`, `engineering`, `IP`, the model's units as `coolsolve`), region
  preset, line families (on/off, count, rounding), marked critical point, saturation
  scale, theme (`default` or `dark` following the GUI, `print`, `coolpack`), label mode
  (rails, edge, margin, inline).
- The SVG and its interaction layer (§3.2): rubber-band zoom and reset send a new view,
  a resize sends the new figure size, and the GUI sends back the labels of the scene on
  screen.
- The former overlays become cycle definitions with straight paths: the inferred states
  (name heuristic) and the arrays, until phase 2.
- `diagrams.json` in the ZIP bundle and `showSessionModel()` (gui.md §14.4 "Diagram state
  persistence"). Specs written by the GUI set `branding` to `{"references": "none",
  "versions": false}` until phase 4 computes the badge lines, so that an exported spec
  draws in CoolPlot as in the GUI.
- The LaTeX report: `exportPlots.ts` looks only for Plotly charts; it gets a branch that
  rasterises the SVG of the diagram (`.thermo-diagram svg`) to `thermo_diagrams.png`.
- Remove `POST /api/v1/coolprop/saturation` and `ThermoDiagram.tsx`.

### Phase 2 — Cycles from the model (done)

Done on `feat/coolplot-diagrams`: `diagram/cycles` (the paths of CoolPlot's `paths.py`),
`diagram/states` (§4.2, validated on the library), `POST /api/v1/diagram/state`,
`GET /api/v1/diagram/states`, the state editor (`CycleEditor.tsx`, removed on 2026-11-13: D6),
the hover (`useDiagramHover.ts`). Follow-ups: two calls of one state that share a single variable
(`h_1 = enthalpy(P=P_1, x=x_1)`, `s_1 = entropy(T=T_1, x=x_1)`) stay two states, which the
structure cannot resolve; models with many property calls give many small states (17 for an
ORC of the library), which makes *From inferred states* noisy; CoolPlot's legend, placed
automatically, can cover a cycle on a short chart. What the list below planned:

- `diagram/cycles`: full states from the solved values (§4.3), paths as chains (isobars
  through the dome, isenthalps, isentropes, isotherms, polytropic paths), state numbers
  and labels, one legend entry; coolplot-core keeps the labels of the lines off the
  cycles.
- The grouping rule of §4.2 (`diagram/states`, `GET /api/v1/diagram/states`) and the state
  editor.
- Hover: the variable names and values of a state marker; elsewhere, the full state under
  the cursor from `POST /api/v1/diagram/state`, debounced.
- Closes `CS-FEAT-OVERLAY-STATES`.

### Phase 3 — More fluids (done)

Done on `feat/coolplot-diagrams`:

- **Blends**: the lines not yet cached are computed on up to 4 threads; CSL-0324
  (propane/isobutane, `DEF`) went from 19.7 s to 6.0 s, `R32[0.7]&R125[0.3]` p-h from 1.54 s
  to 0.60 s; new charts of mixtures use `ACHP` (2.8 s for CSL-0324). The lines of constant
  quality of that blend stall near 150 K and 19 Pa in `CoolProp::Plot`'s mixture stepper
  ("evaluation budget used up"): an upstream issue, reproduced by
  `PropertyPlot("HEOS::Propane[0.5]&IsoButane[0.5]", iP, iHmass, 114.87, 612.5, 0.4019, 9.429e6)
  .calc_isolines(iQ, {0.5}, {})`.
- **Incompressibles**: p-h, p-s, p-ρ and p-T charts (in T-s, h-s, T-ρ and T-h every state lies
  near one line), with the concentration in the definition; CoolPlot refuses relative bounds
  for them, so CoolSolve gives a default box: the fluid's temperatures from its freezing point,
  1 to 100 bar, overridden by the spec's absolute limits. The concentration of a new chart is
  that of the model's calls, else the middle of the fluid's range.
- **Tabular backends**: tested; lines of constant quality with the base backend.
- **Ideal gases** (`CS-FEAT-DIAGRAM-IDEAL`): `diagram/ideal_gas`, p-h, T-s, h-s, p-s, p-ρ, T-ρ
  and p-T (T-h is refused: an isotherm is an isenthalp) over 200-2000 K and 0.1-100 bar; solved
  air-standard cycles lie on their isobars and isentropes. p-ρ stands for the register's P-v.
  The exported spec marks the fluid as CoolSolve's ideal gas (`fluid.coolsolve`): CoolPlot
  would draw the real fluid from it.
- **Humid air** (`CS-FEAT-PSYCHRO`): `diagram/humid_air`, a port of CoolPlot's
  `thermo/humid_air.py` (`HAPropsSI`, lines on up to 4 threads), the psychrometric (`WT`) and
  Mollier h,x (`ThxW`) charts at the total pressure of the model; air processes from the solved
  states, drawn straight as CoolPlot does; a state solved at another total pressure is left
  out with a warning. `/api/v1/coolprop/fluids` names humid air `AirH2O` (it listed `Air_ha`,
  the name of the real fluid, twice).

What the list below planned:

- **Blends and CoolProp mixtures**: the port traces them and draws the dome from the
  phase envelope; the fluid string is enough. Measure the time of a blend chart (CoolPlot
  measured 5 s for the default p-h chart of R410A in Python) and, if needed, compute the
  lines of one chart on several threads.
- **Incompressibles**: absolute limits (`PropertyPlot`'s absolute constructor).
- **Tabular backends** (`TTSE&HEOS`, `BICUBIC&HEOS`): the lines of constant quality are
  computed with the base backend, as CoolPlot's Python does (`native.py`); the others with
  the model's backend.
- **Ideal gases** (`CS-FEAT-DIAGRAM-IDEAL`): a CoolSolve provider over `IdealGasFluid`
  (EES/JANAF offsets) writes the curve records of T-s, h-s and p-ρ charts over the gas
  region, with absolute limits. The register asks for P-v: CoolPlot has no `v` axis, so
  either it is proposed upstream or p-ρ stands for it.
- **Humid air** (`CS-FEAT-PSYCHRO`): coolplot-core composes the psychrometric (`WT`) and
  Mollier (`ThxW`) charts; their lines come from `diagram/humid_air` (or from
  `CoolProp::Plot`, D2); air processes from the solved states.
- Later, through the same provider interface: `NH3H2O` and LiBr-H2O charts (Dühring,
  h-x), proposed upstream as new chart kinds first.

### Phase 4 — Outputs (done)

Done on `feat/coolplot-diagrams`: `diagram/references` (checked against CoolPlot's Python on
its 335 BibTeX keys and 137 fluids, in both styles), `POST /api/v1/diagram/export`, the *SVG*,
*PDF* and *Spec* buttons of the tab, the report's figures (`/api/v1/latex/compile` draws every
stored definition with the last solve, in the print theme, as `thermo_diagrams.pdf`,
`thermo_diagram_2.pdf`, …; the `.tex` takes them with `\IfFileExists` and falls back to the PNG
of the browser, kept because the tab's default chart is not stored until it is edited and a
server without fonts writes no PDF), and `coolsolve --diagram`. The versions line is
CoolPlot's (`CoolProp 8.1.0dev (8eb3b77), CoolPlot 1.5.0`, the core's version): it does not
name CoolSolve, so that an exported spec draws the same badge in CoolPlot. New definitions
write no `branding`, i.e. CoolPlot's defaults; definitions saved in phases 1 to 3 keep
`{"references": "none", "versions": false}`.

**Fonts.** The PDF writer embeds TrueType fonts it finds on the platform (the themes ask for
Inter, Helvetica Neue, Arial, DejaVu Sans). A Linux server needs `fonts-dejavu-core`. Where
`fonts-dejavu-extra` is installed too, coolplot-core picks `DejaVuSans-ExtraLight.ttf` for the
regular face (the first file whose bold and italic flags match): thin text, a coolplot-core
issue to fix upstream.

What the list below planned:

- The LaTeX report takes `thermo_diagrams.pdf` from coolplot-core's PDF writer instead
  of a PNG. It embeds TrueType fonts found on the platform; the online server needs a
  font package (CoolPlot's CI uses `fonts-dejavu-core`).
- SVG and PDF downloads in the tab, and *Export chart spec*: the spec with the model's
  cycle (the `spec` of the scene response), which CoolPlot reads
  (`DiagramSpec.from_dict`), for matplotlib figures in Python.
- `CS-FEAT-PLOT-CLI` (`coolsolve --diagram diagrams.json model.eescode -o chart.svg`,
  `.pdf`, for README figures).
- The badge: `diagram/references` computes the references of the equation of state from
  CoolProp's BibTeX keys (`get_BibTeXKey(fluid, "EOS")`) and its bibliography
  (`CoolPropBibTeXLibrary.bib` in the CoolProp sources, embedded at build time), and the
  versions line (CoolProp, CoolPlot from `coolplot_version()`, CoolSolve); the GUI then
  writes CoolPlot's defaults.

---

## 6. Tests and keeping in step with CoolPlot

### 6.1 Fixtures

`tools/coolplot_fixtures.py`, modelled on CoolPlot's `cpp/tests/make_data.py`
(`compose_case`), runs CoolPlot (Python, at the pinned commit, in a virtual environment)
on a list of chart specs and writes for each the compose record and the scene to
`tests/fixtures/coolplot/`, with the CoolPlot commit, the Python CoolProp version
(whether CoolPlot used the port or its Python numerics) and the date in a header file.
The fixtures are committed; the C++ tests (Catch2 `[diagram]`) read them without Python.
Specs: R134a (no upstream golden), water, CO2 and R290 in the fluid diagram types, T-h
included (no upstream golden); `DEF` and `ACHP` regions; every unit system; rounding on
and off; a family placed through saturation states; marked points and a saturation
scale; one cycle (phase 2).

- **Compose records**: the curve requests (ids, keys, values, points, labels, anchors)
  exactly; the domain, the curves and the cycle chains within a relative tolerance,
  because the Python CoolProp wheel is not the pinned C++ commit; the same gaps. This is
  a comparison with another implementation, so the tolerance may be set from measured
  differences; no numerical tolerance of the glue itself is relaxed to make a fixture
  pass.
- **Scenes**: ids, groups, roles, z order, styles and texts as CoolPlot's. Label
  positions follow the curves, so they are compared within the tolerance of the curves.
  That coolplot-core composes what CoolPlot's Python composes is tested upstream;
  CoolSolve tests that it calls it right.
- No renderer fixtures: the GUI shows coolplot-core's SVG.

Measured with the first fixtures (CoolPlot 1.5.0 on the CoolProp 8.0.0 wheel, which has no
port, so CoolPlot used its Python numerics; five cases: R134a p-h, water T-s, R134a T-h,
R290 p-h with `ACHP`, `engineering` and rounding, CO2 p-h with marked points and a saturation
scale): the plans on CoolPlot's domain match bit for bit and every curve has the same number of
points; the largest differences are 7.8e-14 of the span for the domain, 4.6e-8 of the axis span
for the curves (the fixtures keep 8 digits and every third point) and 1.9e-10 for the labels.
The comparison tolerances are set from these: 1e-12, 2e-7 and 1e-9. Two humid-air fixtures
(phase 3: `psychrometric_ac`, a cooling and reheating process at sea level; `mollier_heating`
at 0.9 bar) match to 9.1e-8 of the axis span, with the same number of points per curve. A sixth fixture (phase 2),
an R134a refrigeration cycle given by (p, h) with every path `auto`, checks the paths: its 84
states match CoolPlot's to 4.3e-10 of the span (tolerance 5e-9). Scenes match in ids,
drawing order, group, role, z order, legend, style and text. The ids of families spaced evenly
over a range computed by CoolProp can differ in the twelfth digit, so ids on CoolSolve's own
domain are compared to 1e-11.

### 6.2 Invariants

Independent of CoolPlot, on every fixture fluid:

- the dome closes at the critical point; liquid and vapour lines meet there;
- an isoline crossing the dome contains the saturation states; inside it, an isobar of
  a T-s chart and an isotherm of a p-h chart are horizontal;
- no point is drawn where the flash failed;
- consistency with the solver: a state solved by CoolSolve lies on the isolines of its
  own properties, with the model's backend.

### 6.3 Performance

Measured in a Release build, in process, 250 points per line:

| Chart | Cold (lines computed) | Warm (compose and SVG) |
|---|---|---|
| R134a p-h | 28 ms | 3 ms |
| Water T-s | 70 ms | 2.8 ms |
| CO2 p-h | 145 ms | |
| R410A p-h | 129 ms | |
| R32[0.7]&R125[0.3] p-h | 1.55 s | |

Over HTTP, R134a p-h takes 47 ms cold and 10 ms warm. The state under the cursor takes 9 µs
in process (0.5 ms over HTTP); the grouping of states at solve time a median 0.013 ms per
library model (95 %: 0.12 ms). The `[diagram]` tests check a loose
budget (cold under 1 s, warm under 100 ms).

With the lines on up to 4 threads (phase 3; cold, in process):

| Chart | One thread | Four threads |
|---|---|---|
| R134a p-h | 31.5 ms | 17.2 ms |
| Water T-s | 72 ms | 31 ms |
| R410A p-h | 87 ms | 36 ms |
| R32[0.7]&R125[0.3] p-h | 1.54 s | 0.60 s |
| Propane[0.5]&IsoButane[0.5] (CSL-0324) p-h, `DEF` | 19.7 s | 6.0 s |
| Propane[0.5]&IsoButane[0.5] (CSL-0324) p-h, `ACHP` | 10.3 s | 2.8 s |

The psychrometric and Mollier charts take 138 and 133 ms cold, 2.4 and 3.7 ms warm; the state
of humid air under the cursor 0.7 to 0.85 ms over HTTP.

### 6.4 Moving to newer versions

- **CoolPlot**: bump `COOLSOLVE_COOLPLOT_PINNED_TAG`, regenerate the fixtures, adapt the
  glue to what the failing tests point to, update the header of this document and note
  it in [versions.md](versions.md). coolplot-core ignores spec fields it does not know,
  so the GUI offers a new CoolPlot feature only once the pin has it.
- **CoolProp**: the pin moves for the solver ([contributing.md](contributing.md#pinned-coolprop-version))
  and takes the numerics with it; the `[diagram]` tests run with the other reports.

---

## 7. Decisions

| ID | Question | Decision |
|---|---|---|
| D1 | How the GUI draws the scene | coolplot-core's SVG with a TypeScript interaction layer (§3.2). It also serves the CLI, the report and pages of model folders. |
| D2 | Where the glue lives | In `src/diagram/`, with CoolPlot's module boundaries. Nothing has moved upstream yet: coolplot-core has no `plan()` (CoolPlot's review lists "line values with rounding" for it), `CoolProp::Plot` has no cycle paths or humid-air lines (CoolPlot's C0 plan had them). A part that moves upstream later is dropped here then. To agree with the CoolPlot maintainers, with a request to document the compose record and the 1.5 scene ids in `contracts.md`. |
| D3 | CoolProp pin | Settled: `8eb3b778` (`d9545e6`). |
| D4 | nlohmann/json | Settled: CoolSolve uses 3.12.0, the vendored version of coolplot-core (`5993f9d`). An upstream option to use the caller's `nlohmann_json` target is still worth asking for. |
| D5 | A TypeScript test runner | Not now: the interaction layer maps coordinates only; the GUI is tested in headless Chromium as `tests/test_gui_shortcuts.py` does. |
| D6 | Grouping of states | Automatic (§4.2), and arrays. The state editor that corrected and named the states (phase 2) was removed on 2026-11-13 at the maintainer's request (the tab was too complex). |
| D7 | Units of the Diagrams tab | The unit system of the model, always (2026-11-13, at the maintainer's request; before: a selector, `engineering` by default): today CoolSolve's fixed units (°C, Pa, J/kg), the custom unit record `model_units.record` of the options; a later unit system per model takes its place. A definition saved in other units is drawn in the model's. |
| D8 | Default look of the lines (2026-11-13, after the first charts: the isolines of CoolPlot's themes, 0.75 pt and opaque, competed with the saturation curve, 1.5 pt, and the cycles, 2 pt) | New definitions carry CoolPlot `style_overrides` drawing every family of lines at 0.5 pt and 55 % opacity (*Line weight: light*; *as the theme* removes them). Compared on p-h and T-s charts of a solved R134a cycle, light and dark themes, against 0.4 pt at 45 % (the lines of constant quality fade out). A spec field of CoolPlot, so an exported spec draws the same; the report's print theme takes it too. |

Settled in the first revision: compose and layout are in C++ (coolplot-core); the
numerics are shared through CoolProp; T-h charts are upstream (`Th`); `EUR` was renamed
`engineering` and custom unit records exist.

---

## 8. Documentation

Updated with each phase: [gui.md](gui.md) (§4.7, §6.8, §7.2, §8.4, §10, §11, §14.1, §14.4),
[model_library_support.md §3](model_library_support.md) (`CS-FEAT-OVERLAY-STATES`,
`CS-FEAT-DIAGRAM-IDEAL`, `CS-FEAT-PSYCHRO`, `CS-FEAT-PLOT-CLI`; `CS-BUG-HUMIDAIR-NAME` in §5),
[contributing.md](contributing.md) (file map, the CoolPlot pin, fonts),
[versions.md](versions.md) (next release: CoolProp `8eb3b778`, CoolPlot `7338f48`),
[deployment_ubuntu_apache.md](deployment_ubuntu_apache.md) (`fonts-dejavu-core`), README.

## 9. Open points

- **CoolSolve**: `CS-BUG-HUMIDAIR-NAME` ([model_library_support.md](model_library_support.md)
  §5): the variable inference names humid air `Air_ha`, the real fluid, so humid-air variables
  get the initial guesses of real air; the Diagrams module names it `AirH2O` itself. A solver
  change, left for review. The total pressure of a humid-air chart is not found when the model
  writes it as a literal inside the call; `CALL psychprops` outputs are not grouped into states.
- **Upstream (CoolPlot, CoolProp)**: the stalled quality lines of the propane/isobutane blend
  (§5, phase 3); the regular face chosen among the DejaVu files (§5, phase 4); documenting the
  compose record and the 1.5 scene ids in `contracts.md`; an option of coolplot-core to use the
  caller's `nlohmann_json` (D4); `plan()`, cycle paths and humid-air lines upstream (D2).
- **Not verified**: the MSVC build of coolplot-core and of `src/diagram/` (CoolPlot's CI builds
  the core with GCC only; `references.cpp` declares CoolProp's `get_mixture_binary_pair_data`
  itself instead of including an internal header).
- **GUI**: no control of the badge (an old definition keeps none); one diagram per model in the
  tab (the document holds a list); the legend placed by CoolPlot can cover a cycle on a short
  chart, and its entries are cut ("Wet-bulb temperature T_…"); models with many property calls
  give many small states.
- **Later**: `NH3H2O` and LiBr-H2O charts (proposed upstream first); diagrams kept in library
  model folders and shown by their HTML pages (§3.4).
