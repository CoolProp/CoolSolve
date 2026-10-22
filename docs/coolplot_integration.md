# CoolPlot Integration Plan

**Status:** proposal, revision 2 (2026-10-22), not started.
**References:**

| Project | Version | Commit | Holds |
|---|---|---|---|
| [CoolPlot](https://github.com/CoolProp/CoolPlot) | 1.2.0, [PyPI](https://pypi.org/project/coolplot/1.2.0/) (2026-10-22) and git tag `v1.2.0` | `03c71b5` | the Python reference, `coolplot-core` 1.2.0 (C++17, `cpp/`), the JavaScript SVG writer |
| [CoolProp](https://github.com/CoolProp/CoolProp) | 8.1.0 development, `master` | `2d51b57e` | `CoolProp::Plot`, the C++ port of CoolPlot's numerics (PR #3512) |

One CoolPlot version is used in two forms: the Python package from PyPI
(`coolplot==1.2.0`) for the test fixtures, and the C++ library `coolplot-core` from the
git tag of the same version (`v1.2.0`), since the PyPI distribution holds only the Python
package.

This document plans how the **Diagrams** tab of the GUI adopts CoolPlot, the
property-diagram package of the CoolProp project, instead of growing its own plotting
code. It records the decision (section 1), the principles and the rule that decides in
which project a piece of work is done (section 2), the architecture (section 3), how a
solved model becomes a cycle on a chart (section 4), the phases (section 5), the tests
(section 6), the work to do in CoolPlot and CoolProp (section 7) and the open decisions
(section 8).

**What changed since revision 1** (2026-10-20, written against CoolPlot 0.4): revision 1
planned to port CoolPlot's Python code to C++ and TypeScript inside CoolSolve, about
8,000 lines. Since then the numerics have been ported into CoolProp (`CoolProp::Plot`)
and the presentation into CoolPlot itself (`coolplot-core`), both tested against the
Python reference. What is left is the glue, a few pieces CoolPlot still has in Python
only, and what is specific to CoolSolve. Revision 1 is in the git history of this file.

---

## 1. Context and decision

### 1.1 The Diagrams tab today

`gui/src/components/ThermoDiagram.tsx` (about 450 lines, Plotly.js) draws:

- the saturation dome of one `Real` fluid in T-s, P-h, h-s or T-h, from
  `POST /api/v1/coolprop/saturation` (`src/server.cpp`: N temperatures between the
  triple and the critical point, liquid and vapour arrays in SI);
- the solved state points, grouped by a name heuristic (`T_1`, `P_1`, `h_1` give state
  `1`) from `GET /api/v1/variables/inferred`;
- one array path (`X[i]`, `Y[i]`), optionally closed.

There are no isolines, no unit choice (deg C, Pa, J/kg only), no labels, and only `Real`
fluids are offered. The wanted extensions are listed in [gui.md section 14.4](gui.md)
(isolines, units, hover, persistence) and in the feature register of
[model_library_support.md section 3](model_library_support.md): `CS-FEAT-DIAGRAM-IDEAL`,
`CS-FEAT-PSYCHRO`, `CS-FEAT-OVERLAY-STATES`, `CS-FEAT-PLOT-CLI`.

### 1.2 What exists upstream

CoolPlot's pipeline is `spec -> plan -> curves -> compose -> layout -> scene ->
renderer`, with three JSON contracts between the layers (CoolPlot
`docs/reference/contracts.md`): the chart spec, the curve record and the scene. Its
Python code is the reference; the C++ and JavaScript code is tested against it, number
for number and byte for byte.

| Layer | Python reference (CoolPlot 1.2) | C++ or JavaScript | Status |
|---|---|---|---|
| Isolines, dome, quality lines, mixture tracer, adaptive stepping, T-p presets, `value_at` (hover) | `thermo/fluid.py`, `tracer.py`, `stepping.py`, `isolines.py` | `CoolProp::Plot` in CoolProp (`include/CoolProp/plotting/CoolPropPlot.h`, `src/CoolPropPlot.cpp`, 2,700 lines, Catch2 `[Plot]`) | done; the Python numerics have not changed since the port |
| Quantities, units, number format, diagram types | `core/` | `coolplot-core`: `core.h` | done |
| Themes, styles | `style/` | `coolplot-core`: `style.h` | done |
| Scene and its JSON | `scene.py` | `coolplot-core`: `scene.h` | done |
| Compose: display units, texts, labels on rails, edge and margin labels, automatic legend corner, state numbers and state labels, lines placed by saturation states, the Mollier chart | `compose/builder.py`, `layout.py`, `rails.py`, `corners.py` | `coolplot-core`: `compose.h`, `layout.h` | done |
| SVG and PDF writers, tick placement | `render/svg.py`, `ticks.py` | `coolplot-core`: `svg.h`, `pdf.h` (PDF with TrueType subsets); `render/static/scene_svg.js` (SVG in the browser) | done |
| C interface | | `coolplot-core`: `capi.h` (JSON text in, scene, SVG or PDF out) | done |
| Plan: which lines, round values | `engine/planning.py` (232 lines) | none | **Python only** |
| Values of lines placed by saturation states | `thermo/provider.py: saturation_anchor_values` | none | **Python only** |
| Engine and curve cache | `engine/engine.py`, `cache.py`, `store.py` | none | **Python only** |
| Cycle states and the paths between them | `thermo/cycles/paths.py`, `thermo/process.py`, `core/cycle.py` | none | **Python only** |
| Humid-air lines and states (psychrometric and Mollier charts) | `thermo/humid_air.py` (569 lines) | none | **Python only** |
| Badge with the references of the equation of state | `compose/branding.py`, `thermo/references.py` | none | **Python only** (`coolplot-core` refuses a branded chart, Python composes it) |
| Plotly renderer | `render/plotly.py` (790 lines) | none | **Python only** |

The two C++ parts are split by what they depend on: `CoolProp::Plot` needs an
`AbstractState` and knows nothing about units, styles or drawing; `coolplot-core` needs
no CoolProp and does no property calculation. CoolPlot's plan for the C++ packages
(`docs/contributing/v1_plan.md` section 3, `v1_review.md` section 6.2) puts the cycle
state resolution and the humid-air lines in `CoolProp::Plot` and the line values with
rounding in `coolplot-core`; these are the "Python only" rows above that are not done yet.

### 1.3 Options considered

| Option | Verdict |
|---|---|
| **A. Run CoolPlot's Python code** (embedded CPython, a Python process next to the server, or Pyodide in the browser) | Rejected. The Windows build is one self-contained `.exe`; Python, numpy and the CoolProp wheel would have to ship with it. The curves would come from a second CoolProp, not the pinned one, without CoolSolve's reference-state offsets, and CoolSolve's own fluids (`NH3H2O`, LiBr) could never be drawn. |
| **B. Keep developing our own plots** | Rejected. Everything in section 1.1 has been designed, measured and reviewed in CoolPlot. |
| **C. Port CoolPlot to C++ and TypeScript inside CoolSolve** (revision 1) | Superseded: most of the port now exists upstream. A second copy would drift from it. |
| **D. Link `CoolProp::Plot` and `coolplot-core`, add the missing pieces where they belong, write the CoolSolve glue** | **Chosen.** |

---

## 2. Principles

1. **CoolPlot is the reference.** Contracts, item ids, groups, roles, styles, texts and
   numerical methods follow CoolPlot. CoolSolve follows one CoolPlot release from PyPI
   (the header of this document) as it follows one pinned CoolProp commit; moving to a
   newer one is a deliberate change (section 6.4).
2. **Properties come from CoolSolve's pinned CoolProp**, never from a second property
   library: the same build of CoolProp computes the solved states and the lines, with the
   same reference states (`Fluid::getReferenceState`). A state point of a solved model
   lies exactly on the lines of its own properties.
3. **The contracts are CoolPlot's, unchanged.** The GUI sends a chart spec in CoolPlot's
   format and receives a scene in CoolPlot's format. A chart spec exported from CoolSolve
   draws in CoolPlot's Python package.
4. **Do the work where it belongs** (section 2.1). A fix to a defect found in an upstream
   part is made upstream, never in a patched copy in CoolSolve.
5. **No Python at run time or at test time in CoolSolve.** Upstream tests its own parts;
   CoolSolve commits the few end-to-end fixtures it needs (section 6.1).
6. **Not solver core.** The work touches the build, the server, a new `diagram` module
   and the GUI, not the parser, evaluator or solvers, with one exception: the CoolProp
   pin (section 5, phase 0). It is developed on a branch (`chp/coolplot-diagrams`) and
   reviewed before merging.

### 2.1 Which project does the work

| The piece of work | Goes to | Because |
|---|---|---|
| Needs an `AbstractState` or `HAPropsSI` and works for any CoolProp fluid: lines, saturation states, the states along a process, humid-air lines, `value_at` | **CoolProp** (`CoolProp::Plot`), with the Python reference in CoolPlot | Every CoolProp wrapper gets it; CoolPlot's Python already calls `CoolProp::Plot` when present (`thermo/native.py`). |
| Does no property calculation: planning and round values, compose, layout, styles, texts, SVG and PDF, the C interface | **CoolPlot** (`coolplot-core`), with the Python reference | One implementation of the presentation for every language, tested byte for byte against Python. |
| Draws a scene in a browser | **CoolPlot** (JavaScript next to `render/static/scene_svg.js`) | The reference (`render/plotly.py`) and the Python to generate expected outputs are there; CoolPlot's own browser backend can use it. |
| Build-system and packaging problems of the upstream code met by CoolSolve | **CoolPlot** or **CoolProp** | Any project that embeds them meets the same problems. |
| Knows CoolSolve: models, the solve, variable names, state grouping, fluids outside CoolProp (`NH3H2O`, LiBr ASHRAE), the EES/JANAF reference states of ideal gases, sessions, the server, the GUI, the ZIP bundle, the LaTeX report | **CoolSolve** | Nothing outside CoolSolve can use it. |
| A new chart kind (h-x of ammonia-water, Duehring chart of LiBr-water) | numerics in **CoolSolve**, chart kind (axes, families, texts) in **CoolPlot** | The presentation of a chart kind is CoolPlot's vocabulary; the properties are CoolSolve's. |

When CoolSolve needs an upstream piece before it is there, CoolSolve may carry an
**interim implementation** in `src/diagram/interim/`: a port of the Python reference with
the upstream names and the same JSON, the CoolPlot commit it was ported from in its
header, and a line in section 7 naming the upstream item that replaces it. It is deleted
when CoolSolve moves to the upstream version. An interim implementation is never
extended beyond the Python reference; a change goes upstream first.

---

## 3. Architecture

### 3.1 Data flow

```
GUI (Diagrams tab)                     coolsolve server (C++)
------------------                     ------------------------------------------------
chart spec (CoolPlot JSON)   POST      spec --provider--> domain (tp_box, ranges_SI)
+ cycles from the solve      ------->       --plan------> curve requests     [coolplot-core *]
+ plot size in pt, view      /api/v1/       --curves----> curve records (SI, cached)
+ previous labels            diagram/          CoolProp::Plot, or a CoolSolve provider
                             scene          --paths-----> cycle chains       [CoolProp::Plot *]
                                            --compose---> scene              [coolplot-core]
Plotly.js <- scene_plotly.js <- scene JSON <-------------
            [CoolPlot JS *]
                                       * = not upstream yet: interim implementation
```

- `POST /api/v1/diagram/scene` takes the chart spec, the cycles, the plot size, the
  view and the previous label items (so that labels stay where they were), and returns
  the scene and the curve diagnostics (failed points, fallbacks, the `reason` of lines
  that could not be computed).
- `POST /api/v1/diagram/state_at` returns the full state under the cursor
  (`PropertyPlot::value_at` for each quantity) for hover.
- `GET /api/v1/diagram/options` lists what the GUI offers: diagram types per fluid kind,
  unit systems, region presets, isoline families per diagram.
- `GET /api/v1/diagram/svg` and `.../pdf` (phase 5) return the same chart through the
  writers of `coolplot-core`.
- `GET /api/v1/coolprop/fluids` and `GET /api/v1/variables/inferred` stay;
  `POST /api/v1/coolprop/saturation` is removed with the old tab.

Changing units, theme, labels, the view or the plot size recomputes no property: the
curves come from the cache and only compose runs again (12 ms for a p-h chart with labels
in CoolPlot's measurement, including its JSON reading and writing).

### 3.2 Fluids

| CoolSolve fluid (`FluidType`) | Lines from | Notes |
|---|---|---|
| `Real` | `CoolProp::Plot` with the CoolProp name | reference offsets are zero |
| `Mixture` (`A[x]&B[y]`) | `CoolProp::Plot` (tracer, dome from the phase envelope) | predefined blends as `R407C.mix`, not the pseudo-pure `R407C` |
| `Incompressible` | `CoolProp::Plot`, absolute limits | no dome |
| `IdealGas` | `CoolProp::Plot` over the gas region, absolute limits, then the h, u and s axes and isoline values shifted by `IdealGasFluid::getReferenceState` | CoolSolve's EES/JANAF basis; CoolSolve only. To verify: temperature-only calls of CoolSolve use the dummy pressure, the lines use the real equation of state at every pressure, so states and lines agree only if the solve used the same pressure |
| `HumidAir` (`AirH2O`) | humid-air lines (psychrometric `WT` and Mollier `ThxW`, both composed by `coolplot-core` 1.2) | lines are Python only today (section 7, U5) |
| `AmmoniaWater` (`NH3H2O`), LiBr ASHRAE | CoolSolve's own curve code (`nh3h2o.h`, `libr_ashrae.h`), curve records in CoolPlot's format | needs chart kinds CoolPlot does not have (section 7, U9) |

A provider interface in `src/diagram/provider.h` hides the difference: domain, lines,
states along a path, state at a point. The engine never looks inside a provider, as in
CoolPlot's `thermo/provider.py`.

### 3.3 Module map

New module `include/coolsolve/diagram/`, `src/diagram/`; new GUI folder
`gui/src/components/diagram/`.

| CoolSolve file | Does | Mirrors |
|---|---|---|
| `diagram/provider.{h,cpp}` | provider interface; providers over `CoolProp::Plot` (one `PropertyPlot` per fluid, chart and thread), ideal-gas offset shift | `thermo/provider.py`, `thermo/native.py` |
| `diagram/engine.{h,cpp}` | domain, plan, curves, cycles, compose input; bounded thread-safe curve cache shared by the sessions | `engine/engine.py`, `engine/cache.py` (not the disk store) |
| `diagram/cycles_from_model.{h,cpp}` | states of a solved model, grouping of inferred variables (section 4) | CoolSolve only |
| `diagram/interim/planning.{h,cpp}` | interim until U2 | `engine/planning.py` |
| `diagram/interim/anchors.{h,cpp}` | interim until U3 | `thermo/provider.py: saturation_anchor_values` |
| `diagram/interim/paths.{h,cpp}` | interim until U4 | `thermo/cycles/paths.py`, `thermo/process.py` |
| `diagram/interim/humid_air.{h,cpp}` | interim until U5 | `thermo/humid_air.py` |
| `diagram/nh3h2o_lines.{h,cpp}` | `NH3H2O` and LiBr curves | CoolSolve only |
| `src/server.cpp` | the endpoints of section 3.1 | |
| `gui/src/components/diagram/` | the Diagrams tab; `scene_plotly.js` from CoolPlot (interim copy in TypeScript until U6) | `render/plotly.py` |

Not used from CoolPlot: the background engine (`engine/background.py`, works around
Python's interpreter lock), the disk store of curves (`engine/store.py`; the server's
cache in memory is enough until a chart is measured to be slow), dragging cycle states
(`core/events.py`, `interaction.py`: in CoolSolve the states come from the solve), the
cycle models (`thermo/cycles/models.py`: in CoolSolve the model is the cycle), the
matplotlib backend and the Python facades.

### 3.4 Threads

A `PropertyPlot` computes all its lines with one `AbstractState` and holds a lock while
it does; two threads that need lines of the same fluid at the same time need two plots.
The engine keeps a small pool of plots per (fluid, chart, region), bounded like the
curve cache. Building a plot is the expensive step for a mixture (critical point and
phase envelope), so plots are reused, never built per request.

### 3.5 Build

- `CoolProp::Plot` is part of the CoolProp library CoolSolve already links; it needs a
  CoolProp pin at or after `2d51b57e` (phase 0).
- `coolplot-core` is fetched with `FetchContent` (`SOURCE_SUBDIR cpp`), static, with
  `COOLPLOT_CORE_TESTS=OFF`, from the git tag `v<version>` of `COOLSOLVE_COOLPLOT_VERSION`
  (the version also pinned on PyPI for the fixtures). zlib is optional (uncompressed PDF
  streams without it).
- `coolplot-core` carries nlohmann/json 3.12 in its public include path; CoolSolve
  fetches 3.11.3. A source file that includes both headers sees one version and the
  library was built with the other: different ABI namespaces, so link errors at best.
  CoolSolve moves to nlohmann/json 3.12.0 (phase 0); upstream gets an option to use the
  host's copy (U1).
- `coolsolve --version` prints the CoolPlot version (`coolplot_version()`) next to the
  CoolProp version.

---

## 4. From a solved model to a cycle on the chart

CoolPlot draws a cycle from numbered states, chains of state numbers (`loops`) and paths
between them (`"auto"`: isobar if the pressures are equal, isenthalp, isentrope,
isotherm, otherwise a straight line). A state may be given as a full state
(`{"state": {...}}`); CoolSolve always gives full states from the solve, so CoolPlot's
`via` resolution (compression with an efficiency, and so on) is not needed. The paths
between states still are (`paths.py`): a condenser follows its isobar through the dome.

### 4.1 Sources of states

1. **Arrays** (today's array overlay): `P[i]`, `h[i]` give state `i`; the loop is the
   index order.
2. **Inferred states**: variables used as arguments or results of the same property
   calls (`inferVariables`, `src/variable_inference.cpp`) are grouped into states
   (section 4.2).
3. **Picked by the user** (`CS-FEAT-OVERLAY-STATES`): a state editor in the Diagrams tab
   where each state is a set of variables (`T_su_cp`, `P_su_cp`, ...) and the loops are
   chains of states. This corrects the automatic grouping and names the states.

Each state carries the names of its variables; CoolPlot 1.2's state labels
(`labels: [[number, text]]`) can show a name next to the state number, and hover shows
the variable names and values, as today. The cycle definition (states as variable names,
loops, paths) is saved in the ZIP bundle with the chart spec (`diagram.json`), not in the
`.eescode`: EES keeps its plots in the GUI too.

### 4.2 Grouping inferred variables into states (decision D6)

Today's name heuristic merges everything that shares a suffix. A rule based on the
equations: each property call on a fluid involves a set of state variables (inputs and
output); two calls describe the same state when they share **two variables of
independent properties** (state postulate for a pure fluid). `h_2 = enthalpy(R134a,
P=P_2, s=s_2)` and `T_2 = temperature(R134a, P=P_2, h=h_2)` share `P_2` and `h_2`: one
state. `T_sat = temperature(R134a, P=P_cond, x=1)` shares only `P_cond` with the
condenser states: a separate state. T and P are not independent inside the dome, so a
pair (T, P) alone does not merge two calls. This needs a test on the library models
before it replaces the heuristic. It is CoolSolve's own work.

### 4.3 Values and units

States go into the spec in SI; the GUI converts the model units on the way. Each state
is drawn from its solved values, never recomputed. One chart belongs to one fluid; a
model with several fluids gets one chart per fluid, chosen in the fluid list as today.

The default unit system is CoolPlot's `engineering` (deg C, bar, kJ/kg; the name was
"EUR" before CoolPlot 1.0). A "model units" choice needs no upstream change: the chart
spec takes a custom unit record (`{"name": "model", "units": {"p": ["Pa", 1.0, 0.0],
...}}`) built from the units of the model.

---

## 5. Phases

Each phase ends with the Diagrams tab usable and the old code it replaces removed. Items
are marked with the project that does them; an upstream item that is not merged when the
phase starts is replaced by its interim implementation (section 2.1).

### Phase 0: foundations

- [CoolSolve] Move the CoolProp pin to `2d51b57e` or later, in its own commit, with the
  procedure of [contributing.md](contributing.md) ("Pinned CoolProp version"). The 17
  commits after the current pin `75af7816` include changes that move solver results:
  DQ/HQ/QS flashes no longer leave the phase imposed (`ca2de185`), the mixture PT-flash
  rework (`a089befb`), an IF97 input fix (`83f3ca35`), PC-SAFT `kij` read as doubles
  (`2eb5c05c`). Every changed test, example and library baseline is explained.
- [CoolSolve] nlohmann/json 3.11.3 to 3.12.0.
- [CoolSolve] Pin CoolPlot 1.2.0: `coolplot==1.2.0` from PyPI for the fixtures, the git
  tag `v1.2.0` for `coolplot-core` (section 3.5).
- [CoolPlot] `coolplot-core` as a sub-project (U1).
- [CoolSolve] Fetch `coolplot-core`, print its version, an empty `src/diagram/` that
  composes a fixed scene in a Catch2 test (`[diagram]`) to prove the build on Linux and
  with `build_installer.bat` on Windows.

### Phase 1: a chart of a pure fluid on the server

- [CoolPlot] `plan()` in `coolplot-core` (U2); interim `planning.cpp` meanwhile.
- [CoolSolve] Providers for `Real`, `Mixture` and `Incompressible` over `CoolProp::Plot`;
  the domain record (`tp_box()`, `isoline_range()` per quantity); the engine and the
  cache; `POST /api/v1/diagram/scene` and `GET /api/v1/diagram/options`.
- [CoolSolve] T-h (decision D3): `coolplot-core` reads `"Th"` as a generic axis pair and
  CoolProp has `HmassT_INPUTS`; check that the lines through the dome are right before
  the tab offers it.
- Tests: section 6.

### Phase 2: the new Diagrams tab

- [CoolPlot] `scene_plotly.js` (U6); an interim TypeScript port in
  `gui/src/components/diagram/` meanwhile, written to be replaced by the upstream file.
- [CoolSolve] The tab: fluid, diagram type, units (including model units), region
  preset, isoline families (count, rounding, lines placed by saturation states, on/off),
  label mode, theme. Today's state and array overlays drawn as CoolPlot cycles with
  straight paths until phase 3.
- [CoolSolve] Zoom and resize: the renderer reports the new view (`plotly_relayout`) and
  plot size, the GUI asks for a new scene with the previous labels; hover through
  `POST /api/v1/diagram/state_at`.
- [CoolSolve] The chart spec is saved in the ZIP bundle and restored (gui.md section
  14.4 "Diagram state persistence"); `exportPlots.ts` keeps exporting
  `thermo_diagrams.png` (the new tab keeps the `thermo-diagram` class).
- [CoolSolve] Remove `POST /api/v1/coolprop/saturation` and `ThermoDiagram.tsx`.

### Phase 3: cycles from the model

- [CoolProp] States along a process in `CoolProp::Plot` (U4); interim `paths.cpp`
  meanwhile.
- [CoolSolve] The three state sources of section 4.1, the grouping rule of section 4.2,
  the state editor, state labels with variable names. Closes `CS-FEAT-OVERLAY-STATES`.

### Phase 4: more fluids

- [CoolSolve] Ideal gases (`CS-FEAT-DIAGRAM-IDEAL`): T-s, h-s and p-rho over the gas
  region with the offset shift of section 3.2, and a test that solved states of a JANAF
  species lie on their isolines. A p-v diagram would be a CoolPlot change (v is not a
  line family there); propose it only if p-rho is not accepted by users.
- [CoolProp] Humid-air lines (U5); interim `humid_air.cpp` meanwhile.
- [CoolSolve] Psychrometric and Mollier charts of `AirH2O` models (`CS-FEAT-PSYCHRO`),
  with the air processes of the model as straight paths.
- [CoolSolve] Lines placed by saturation states need U3 (or its interim).
- [CoolSolve] + [CoolPlot] `NH3H2O` and LiBr-water charts: curves in CoolSolve, chart
  kinds proposed upstream (U9). Only after the other fluids are done.

### Phase 5: outputs

- [CoolSolve] SVG and PDF through `coolplot-core`: a vector PDF figure in the LaTeX
  report instead of `thermo_diagrams.png`; `coolsolve --diagram spec.json
  model.eescode` (`CS-FEAT-PLOT-CLI`) for README figures. Neither needs a renderer
  port any more.
- [CoolSolve] *Export chart spec*: the chart spec with the model's cycle, as JSON that
  CoolPlot's Python reads (`DiagramSpec.from_dict`).
- [CoolPlot] The badge with the references of the equation of state in `coolplot-core`
  (U7), if the reports want it.

---

## 6. Tests and keeping in step

### 6.1 What CoolSolve tests

The numerics are tested in CoolProp (`[Plot]`, mirroring CoolPlot's numerics tests) and
the presentation in CoolPlot (`cpp/tests`, golden files made by the Python composer and
SVG writer). CoolSolve does not repeat those tests. It tests its own parts (Catch2
`[diagram]`):

- **Glue**: fluid names to `CoolProp::Plot` (blends, fractions, incompressibles), the
  domain record, cache keys (a change of CoolProp configuration or reference state never
  serves an old curve), thread safety of the plot pool.
- **End to end**, a handful of charts (R134a p-h with a cycle, water T-s, an R32/R125
  blend, humid air): the scene of the server against the scene CoolPlot's Python makes of
  the same spec, written once by `tools/coolplot_fixtures.py` into
  `tests/fixtures/coolplot/` with the CoolPlot and Python CoolProp versions and the date.
  The tool runs in a virtual environment built from
  `tools/coolplot_fixtures_requirements.txt`: `coolplot==1.2.0` and `CoolProp==8.0.0` from
  PyPI, pinned with their hashes (`pip install --require-hashes`). CoolProp 8.0.0 has no
  `CoolProp::Plot`, so CoolPlot computes the lines with its Python numerics, which are the
  reference. Ids, groups, roles, z order, styles and texts exactly; coordinates within a
  tolerance set from measured differences (the Python CoolProp wheel is not the pinned
  C++ commit). The tolerance of a comparison with another implementation may be set from
  measurements; no numerical tolerance of CoolSolve's own code is relaxed to make a
  fixture pass.
- **Interim implementations**: the golden files of the Python reference for their piece
  (planning: the requests of a spec; paths: the states along each path kind), deleted
  with the interim code.
- **CoolSolve providers**: ideal-gas shift, `NH3H2O` and LiBr curves (closure of the
  dome, saturation states on the lines).

### 6.2 Invariants

On every fixture fluid:

- the dome closes at the critical point; liquid and vapour lines meet there;
- an isoline crossing the dome contains the saturation states; inside it, an isobar of a
  T-s chart and an isotherm of a p-h chart are horizontal;
- no point is drawn where a property call failed;
- consistency with the solver: a state solved by CoolSolve lies on the isolines of its
  own properties (for a `Real` fluid, a mixture and an ideal gas with a JANAF reference
  state).

### 6.3 Performance

Measure, in a Release build, the cold time of the default chart of each fixture fluid
and the warm time of compose; set the budget from the measurements and add it to the
test suite ([contributing.md section 5](contributing.md)).

### 6.4 Moving to a newer CoolPlot or CoolProp

Move to a CoolPlot version published on PyPI (the fixture requirements and
`COOLSOLVE_COOLPLOT_VERSION` together), run the `[diagram]` tests and regenerate the end-to-end fixtures, delete
the interim implementations that upstream now replaces, update the header of this
document and note it in [versions.md](versions.md). A CoolProp bump follows the pin
procedure of contributing.md as well.

---

## 7. Work in CoolPlot and CoolProp

Items CoolSolve needs from upstream, each to be filed as an issue or pull request in the
project named. A change to `coolplot-core` or `CoolProp::Plot` lands together with the
Python reference, as CoolPlot's own rule requires.

| ID | Project | Item | Why it belongs there | Needed in | Interim in CoolSolve |
|---|---|---|---|---|---|
| U1 | CoolPlot | `cpp/CMakeLists.txt` for embedding: tests default to `PROJECT_IS_TOP_LEVEL`; an alias target `coolplot::core`; an option to use the host's `nlohmann_json::nlohmann_json` instead of the vendored copy | every program that embeds `coolplot-core` meets this | phase 0 | set `COOLPLOT_CORE_TESTS=OFF`, align nlohmann/json to 3.12.0 |
| U2 | CoolPlot | `plan(spec, domain)` in `coolplot-core`: curve requests and round values, port of `engine/planning.py`, golden files from Python | no property calculation; CoolPlot's C1 plan lists line values with rounding in the core | phase 1 | `interim/planning.cpp` |
| U3 | CoolProp | values of lines placed by saturation states (`saturation_anchor_values`) as a `PropertyPlot` method | needs saturation states of the fluid | phase 4 | `interim/anchors.cpp` |
| U4 | CoolProp | states along a process (isobaric through the dome, isothermal, isenthalpic, isentropic, isochoric) as a `PropertyPlot` method, port of `paths.py` | needs an `AbstractState`; in CoolPlot's C0 plan ("cycle state resolution") | phase 3 | `interim/paths.cpp` |
| U5 | CoolProp | humid-air lines and states (`humid_air.py`) as a `CoolProp::Plot` class over `HAPropsSI`, including the Mollier temperature scale | needs `HAPropsSI`; in CoolPlot's C0 plan | phase 4 | `interim/humid_air.cpp` |
| U6 | CoolPlot | `scene_plotly.js`: a Plotly.js renderer of the scene, port of `render/plotly.py`, next to `scene_svg.js`, with a node test against the Python backend's figures as `tests/js/check_scene_svg.js` does for SVG | the reference and the Python to make expected outputs are there; CoolPlot's browser backend can use it | phase 2 | TypeScript copy in the GUI |
| U7 | CoolPlot | the badge (`branding.py`, `references.py`) in `coolplot-core`; references from CoolProp's BibTeX keys | presentation; today a branded chart falls back to Python | phase 5, optional | none (no badge) |
| U8 | CoolPlot | merge `chp/1.2`, tag `v1.2.0` | CoolSolve pins tags where it can | phase 0 | **done** (2026-10-22, `03c71b5`) |
| U9 | CoolPlot | chart kinds for absorption pairs: h-x (enthalpy over ammonia mass fraction) and Duehring (log p over -1/T with lines of constant concentration), as diagram types with their families, in Python and `coolplot-core` | the vocabulary of charts is CoolPlot's; the curves come from CoolSolve as curve records | phase 4, last | none |
| U10 | CoolProp | a `PropertyPlot` constructor from an existing `AbstractState` | CoolPlot's `native.py` falls back to Python for states that a fluid string cannot name (changed interaction parameters, mass-fraction bases); CoolSolve could hand over its configured state | later, optional | fluid strings |

Things CoolSolve does not ask upstream for: the reference-state offsets of ideal gases
(CoolSolve's EES/JANAF basis), the state grouping, the server, the curve cache of a
multi-session server.

---

## 8. Decisions

| ID | Question | Decision or proposal |
|---|---|---|
| D1 | Compose and layout in C++ (server) or in TypeScript (browser)? | **Settled:** C++, `coolplot-core`. |
| D2 | Where does the C++ code live? | **Settled:** numerics in CoolProp (`CoolProp::Plot`), presentation in CoolPlot (`coolplot-core`), the glue in CoolSolve; section 2.1 for new pieces. |
| D3 | T-h diagram | `coolplot-core` and CoolPlot's Python read `"Th"` as a generic axis pair; `CoolProp::Plot` accepts any pair CoolProp has a flash for. Verify the lines in phase 1; no upstream change expected. |
| D4 | Default units | `engineering` (deg C, bar, kJ/kg); model units as a custom unit record, no upstream change. |
| D5 | A TypeScript test runner for the renderer | Not needed in CoolSolve if U6 is upstream and tested there; Vitest (dev dependency only) for the interim TypeScript renderer if U6 is late. |
| D6 | Automatic grouping of states or user picking only | Both: automatic by default, the editor to correct and name. Validate the rule on the library models first. |
| D7 | The engine (domain, plan, curves, cycles, cache) in CoolSolve or as an optional CoolPlot target that links CoolProp | Start in CoolSolve with the names of `engine.py` and a provider interface; offer it to CoolPlot once it is stable, if a second C++ program wants it. `coolplot-core` stays free of CoolProp either way. |
| D8 | Plotly.js or SVG in the tab | Plotly.js for the tab (hover, zoom); the SVG and PDF of `coolplot-core` for files and the report. |

---

## 9. Documentation to update when the work lands

- [gui.md](gui.md): section 4.7 (endpoints), 6.8 (Diagrams tab), 8.3 (variable
  inference for diagrams), 14.4 (items done).
- [model_library_support.md section 3](model_library_support.md):
  `CS-FEAT-OVERLAY-STATES`, `CS-FEAT-DIAGRAM-IDEAL`, `CS-FEAT-PSYCHRO`,
  `CS-FEAT-PLOT-CLI`.
- [contributing.md](contributing.md): file map, the CoolPlot pin next to the CoolProp pin.
- [versions.md](versions.md): the CoolPlot version of each release.
- README: features, documentation table, future work, third-party licences (CoolPlot,
  MIT).
