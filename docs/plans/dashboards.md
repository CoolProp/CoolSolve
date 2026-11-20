# Dashboards plan

**Status:** proposal, 2026-11-20, not started (Phase 0 done: this plan, a prototype runtime and
two example dashboards in [`dashboards/`](dashboards/)).
**Scope:** GUI, server, CLI, ZIP bundle, model library. No change to the language or the solver.
**Not served:** this file lives in `docs/plans/`, a subfolder of `docs/`. `cmake/embed_assets.cmake`
embeds the files of `docs/` only (no subfolders) and the server has no filesystem fallback in a
release build (`ServerOptions::docsDir` is never set), so nothing under `docs/plans/` is reachable
from `/docs/` on a CoolSolve server.

A **dashboard** is an HTML page that shows the results of a CoolSolve model the way an engineer
wants to see them: values written on the schematic of the plant, the inputs of the model as
sliders, the diagrams of the session, the *Solve* button, and, for `INTEGRAL` models, an animation
of the results (plots, schematics, a moving engine). It plays the role of the EES *Diagram window*,
but it is an ordinary HTML file: written by an LLM from a generated prompt, or by hand from a
template, opened live in the GUI or on its own as a snapshot, and saved with the model.

| | |
|---|---|
| [§1](#1-requirements) | Requirements, and where this plan answers each of them |
| [§2](#2-what-exists-today) | What CoolSolve has today |
| [§3](#3-concepts) | Concepts: saved plot, dashboard, runtime, host, snapshot |
| [§4](#4-user-experience) | User experience |
| [§5](#5-authoring-contract-v1) | The authoring contract (elements, bindings, expressions, JavaScript API, layout) |
| [§6](#6-data-contract) | The data contract (what a dashboard receives) |
| [§7](#7-inputs) | Inputs: the parametric-table rule, how a change reaches the model |
| [§8](#8-animation) | Animation of `INTEGRAL` models |
| [§9](#9-architecture) | Architecture: runtime, host, security, server, bundle, CLI, plots |
| [§10](#10-llm-workflow) | The LLM workflow (prompt, checker, agents with the binary) |
| [§11](#11-test-models) | Library models selected to test the feature |
| [§12](#12-phases) | Phases |
| [§13](#13-integration-checklist) | Integration checklist (`docs/contributing.md` §3) |
| [§14](#14-ideas) | Ideas that make the feature attractive |
| [§15](#15-decisions) | Decisions taken while hesitating |
| [§16](#16-risks-and-open-questions) | Risks and open questions |
| [§17](#17-documentation-to-update-when-it-lands) | Documentation to update when it lands |
| [§18](#18-the-example-files) | The example files of this folder |

---

## 1. Requirements

| # | Requirement (from the request) | Answer |
|---|---|---|
| R1 | Dashboards are HTML, served by the CoolSolve GUI, and can also be opened on their own to explore results | §4.2, §4.4, §9.1, §9.8 |
| R2 | They can show the diagrams generated in the session | §4.1, §5.2 (`cs-plot`), §9.10 |
| R3 | Several diagrams: keep the current Diagrams window, *save* it, and it becomes a new tab at the bottom | §4.1, D2 |
| R4 | Easy for an LLM to write from good instructions (given when the user clicks *Add dashboard*); doable by hand by editing a simple HTML template, not inside the GUI | §4.3, §4.5, §10 |
| R5 | Variables placed on a schematic, at the right place whatever the window size | §5.6, D4 |
| R6 | Variables are outputs (read from the results) or inputs; inputs follow the rule of the parametric table: only isolated, directly defined variables | §7 |
| R7 | `INTEGRAL` models: animate plots, diagrams and schematics, live, from the integral table, with proper animation buttons; animated drawings (e.g. an internal combustion engine) | §8, example `dashboards/engine.html` |
| R8 | Users are encouraged to go through an LLM with clear instructions; better still if the LLM can run the `coolsolve` binary | §10.3 |
| R9 | The important GUI buttons (Solve, Update guesses, …) can be placed on a dashboard; they do nothing when the page is not served by the GUI | §5.2 (`cs-button`), §9.3 |
| R10 | The plan is in `docs/` but not served online | header of this file, D1 |
| R11 | Dashboards are saved and loaded with the other files of the model (ZIP export/import) | §9.7 |
| R12 | A few library models to test the feature | §11 |
| R13 | Nice features that make CoolSolve attractive; more customisable and practical than the EES Diagram window | §14 |
| R14 | No questions: decisions taken while hesitating are documented | §15 |

---

## 2. What exists today

- **Diagrams tab** — `gui/src/components/ThermoDiagram.tsx` (Plotly): the saturation dome of one
  fluid (`POST /api/v1/coolprop/saturation`, computed in `src/server.cpp`), the solved states
  grouped by name, one array path (`X[i]`, `Y[i]`). Its configuration is not saved
  ([gui.md §14.4](../gui.md), *Diagram state persistence*). A redesign on CoolPlot is planned
  ([coolplot_integration.md](../coolplot_integration.md)); this plan only needs a diagram to be
  describable by a *spec* and drawable from it (§9.10), whatever draws it.
- **Parametric tab** — `ParametricStudy.tsx`: 1-D line plots and 2-D contours of a sweep. The
  sweep variables are the variables flagged `isImposed` by `POST /api/v1/parse` (lines
  `name = number`, scalars only, detected in the parse handler of `src/server.cpp`). Each grid point
  is solved on a copy of the code where the line is rewritten by the lambda `replaceImposedValue`
  of the parametric handler. Each point keeps **all** its variables (`ParametricGridPoint.variables`).
  The last study is saved in the bundle as `parametric_studies.json`.
- **Integral tab** — `IntegralTable.tsx`: the trajectory table of an `INTEGRAL` model and a line plot
  of one column. The trajectory is `GET /api/v1/integral/result` (columnar JSON) and, after a bundle
  round trip, the CSV `<model>-integral.csv`. The `.sol` file also carries it (`# IntegralTable`).
- **Exports of plots** — `gui/src/utils/exportPlots.ts` captures every visible Plotly chart as PNG
  (LaTeX report).
- **ZIP bundle** — `createZipBundle` / `extractZipBundle` in `src/server.cpp` (uncompressed ZIP,
  subfolders allowed: `debug_output/`). The upload ignores the file types it does not know.
- **Sessions** — per-tab sessions (`X-CoolSolve-Session`, `?sid=`), the `/?library=` and
  `/?include=` links that open a model in a window of its own (`gui/src/api/session.ts`,
  [gui.md §5, §7.3, §7.5](../gui.md)).
- **Library** — model folders embedded through an allow-list that mirrors `EMBED_FILES` of the
  library's `tools/build_index.py` (`cmake/embed_library.cmake`: `README.md`, `model.json`,
  `coolsolve.conf`, `*.eescode`, `*.initials`, `*.sol`, `<name>-<table>.csv`, `figures/*`), loaded
  by `LibraryStore::loadModel`.
- **EES Diagram window** — not imported: register row `CS-GAP-DIAGRAM-WINDOW` of
  [model_library_support.md §4](../model_library_support.md) (animation sliders and shapes,
  edit boxes, radio buttons, `$IFNOT DiagramWindow`) with three verified library cases, `CSL-0339`,
  `CSL-0342`, `CSL-0345`; and case (a) of `CS-GAP-IF-NO-CONSTANT` (a `$IF X$='…'` whose string was
  given in the Diagram window). The inputs were copied into the equations as default values.

What EES offers, and where this plan goes further: the EES Diagram window is one window per file,
drawn in EES's own editor, with text items, input boxes, buttons, plots and animations positioned
in pixels. Dashboards are several per model, plain HTML (any drawing, any layout, any library of
the web), written by an LLM in minutes, resizable, live in a second window or on a second screen,
and exportable as a single file that anyone opens in a browser.

---

## 3. Concepts

| Term | Meaning |
|---|---|
| **Saved plot** | A named chart of the session — a thermodynamic diagram, a parametric plot or an integral plot — stored as a **spec** (what to draw, not the picture) and redrawn from the current results. Shown as a tab of the bottom panel; usable in dashboards by its id. |
| **Dashboard** | An HTML file of the model folder, `dashboards/<name>.html`, with its images in `dashboards/assets/`. Written against the contract of §5. |
| **Runtime** | `coolsolve-dashboard.js`: the script a dashboard includes. It finds the data, binds the elements, draws the plots, runs the player and talks to the host. |
| **Host** | The CoolSolve GUI page that shows a dashboard (a bottom-panel tab, or a window of its own). It feeds the data and carries out the actions. |
| **Data** | One JSON document, `coolsolve-dashboard-data/1` (§6): results, inputs, integral table, parametric study, saved plots, solve status. |
| **Hosted mode** | Inside the GUI: live data, working inputs and buttons. |
| **Standalone mode** | Without the GUI: a **snapshot** of the data (inlined in an exported file, or `results.js` next to the dashboard). Inputs and buttons are shown but disabled, with a tooltip. Animation, plots, zoom and the layout helper work. |

---

## 4. User experience

### 4.1 Saved plots: several diagrams, as tabs

- The **Diagrams** tab stays what it is: the scratch diagram of the session.
- A **Save** button (pin icon) on the Diagrams tab, on the plot of the Parametric tab and on the
  plot of the Integral tab asks for a title and stores the current configuration as a saved plot.
- Each saved plot becomes a **tab of the bottom panel**, after *Console · Parametric · Lookup
  Tables · Integral*, with the plot icon and its title. The tab redraws the plot from the current
  results after every solve (a thermodynamic diagram of a cycle follows the cycle). Tab menu:
  rename, duplicate, delete, *pop out* (its own window), *copy id* (to paste in a dashboard).
- The bottom tab bar scrolls horizontally when it is full; the panel can be maximised
  (double-click on its header) so that a tab fills the window.
- Saved plots are saved in the bundle (`plots.json`, §9.7).

### 4.2 Dashboards in the GUI

- A **Dashboards** drop-down in the toolbar: the dashboards of the model, *Add dashboard…* (§4.3),
  *Import HTML…*, and for the dashboard shown: *Export standalone file*, *Download source*,
  *Check*, *Delete*.
- A dashboard opens as a **tab of the bottom panel** (next to the saved plots), with three buttons
  in the tab header: **pop out** (a window of its own, which can go to a second screen),
  **presentation** (full screen, no CoolSolve chrome: lectures, demos) and **reload** (re-read
  the file).
- **Live**: after each solve the dashboard receives the new results; while a solve runs it shows
  the status (`solving`, block *i*/*n*); a failed solve leaves the last results, marked stale.
  A popped-out window stays in sync with the main window (§9.4).
- Inputs changed on a dashboard rewrite the line of the model in the editor, as an undoable edit
  (§7.3). *Solve* on a dashboard is the toolbar *Solve* (it becomes *Try Harder* after a failure).

### 4.3 *Add dashboard*

A dialog with three routes. The first is the recommended one and is shown first.

1. **With an AI assistant.** A text box *"What should the dashboard show?"* (e.g. *"the cycle on
   a schematic with T, P and h at each state, the COP in big, sliders for the evaporation and
   condensation temperatures, and the P-h diagram"*), options (*include the model code*, *include
   the values*, *dashboard name*), and a generated **prompt** (§10.2) with a *Copy* button and a
   short guide: paste it into an assistant (ChatGPT, Claude, Gemini, Le Chat, a local model), save
   the HTML it returns, then **paste or drop the HTML here**. The dialog runs the checker (§10.4):
   problems are listed with a *Copy a fix-up prompt* button that writes the follow-up message for
   the assistant. *Save* adds the dashboard to the model and opens it.
2. **From a template.** Four templates (blank page with header, status and buttons; values on a
   schematic; animation of an integral model; cycle with a saved diagram), filled with the names
   of the model (its inputs, a few outputs, its integral columns). *Download* saves the HTML to edit
   in any text editor; *Import* brings it back (route 3).
3. **Import an HTML file** (and its images): checked, then added.

The GUI never calls an LLM itself (D16): no API key, no cost, works offline and on the public demo.

### 4.4 Standalone use

- **Export standalone file** (GUI) or `coolsolve --dashboard-export` (CLI, §9.9) writes **one HTML
  file** with everything inside: the runtime, the plotting library, the data snapshot and the images
  (data URIs). It opens offline in any browser, can be e-mailed, attached to a report, put on a web
  page or on a course platform. Inputs and buttons are visible but disabled ("Open this dashboard in
  CoolSolve to change inputs and solve").
- **The unzipped bundle** works too: the bundle carries `dashboards/coolsolve-dashboard.js` and a
  snapshot `dashboards/results.js` next to the dashboards, so a double-click on
  `dashboards/<name>.html` shows the last results, from `file://`. Plots need the plotting library:
  it is loaded from a pinned CDN URL when the folder copy is opened online; offline, use the
  standalone export (D13).
- **Library models**: the model page of `/library/` shows *View dashboard*, a standalone snapshot
  built on the fly from the model's `.sol` baseline (values and integral table) without solving,
  next to *Open in CoolSolve* (§9.11).

### 4.5 Writing or editing a dashboard by hand

- The source is a plain HTML file. The user downloads it (or a template), edits it in any editor
  (VS Code shows the HTML, CSS, SVG and JavaScript with colours and completion), and imports it
  again. A model opened from disk (command line, file association) reads its `dashboards/`
  folder at each *reload*: edit, save, click *reload*.
- The **layout helper** of the runtime — press <kbd>L</kbd> on a dashboard, or open it with
  `#layout` — shows the design coordinates under the pointer and lets the user drag the placed
  elements; each drop shows and copies the updated tag (`<cs-value x="586" y="323" …>`) to paste
  into the file. It is a helper of the page, not an editor of the GUI (R4).

---

## 5. Authoring contract (v1)

The contract is what an author — a person or an LLM — writes. It is small on purpose: a few custom
elements, a few `data-cs-*` attributes, a safe expression language and one JavaScript hook; the
rest is ordinary HTML, CSS and SVG. The prototype `dashboards/coolsolve-dashboard.js` implements it
(Plotly aside) and the two examples use every part of it.

### 5.1 File

```html
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="coolsolve-dashboard" content="1">          <!-- contract version -->
  <title>Heat pump dashboard</title>
  <style> /* the author's CSS; the runtime only sets defaults (--cs-* variables) */ </style>
</head>
<body>
  <cs-dashboard autosolve="false"></cs-dashboard>        <!-- optional settings -->
  … HTML, inline SVG, cs-* elements …
  <script src="coolsolve-dashboard.js"></script>
  <script> /* optional: CoolSolve.onFrame(...) */ </script>
</body>
</html>
```

`<cs-dashboard>` attributes: `autosolve` (`false` by default, §7.4), `results` (file of the
snapshot in standalone mode, default `results.js`), `theme` (`auto`, `light`, `dark`; `auto`
follows the GUI).

### 5.2 Elements

| Element | Attributes | Renders |
|---|---|---|
| `<cs-stage>` | `width`, `height` (design units), `fit` (`width` default, `contain`) | A canvas in design coordinates, scaled as one piece (§5.6). Children with `x`, `y` (and `anchor`) are placed on it. |
| `<cs-value>` | `name` or `expr`, `unit` (display unit, §5.7), `scale`, `decimals`, `label`, `hide-unit`, `warn-above`, `warn-below`, `units` (unit of an `expr`) | Label, number, unit. Tooltip: the variable name and whether it is an input. Red when outside the warning limits. |
| `<cs-input>` | `name`, `label`, `min`, `max`, `step`, `unit`; `widget` (`slider` default, `number`, `select`, `radio`), `options` (string inputs) | An input of the model (§7): slider and number box, or a choice for a string variable. |
| `<cs-button>` | `action`: `solve`, `update-guesses`, `stop`, `reset-inputs`; text content = label | A GUI action (§9.3). Disabled in standalone mode, with a tooltip. `solve` becomes *Try Harder* as in the toolbar. |
| `<cs-status>` | — | `solved` / `failed` / `solving` / `stale` pill, message, time; "snapshot of …" in standalone mode. |
| `<cs-plot>` | `ref` (a saved plot, §4.1) **or** `x`, `y` (columns of the integral table or of the parametric study, `y` may list several), `source` (`integral` default, `parametric`), `x-unit`, `y-unit`, `log-x`, `log-y`, `x-range`, `y-range`, `cursor`, `band-x` + `band-label`, `title` | An interactive chart (Plotly: zoom, hover). With `cursor`, a marker follows the player. |
| `<cs-player>` | `duration` (s for the whole range at ×1, default 10), `realtime`, `autoplay`, `no-loop`, `source` (`integral` default, `parametric`) | Animation controls (§8): start, step back, play/pause, step forward, scrubber, position, speed, loop. |
| `<cs-table>` | `names` (list), `units`, `decimals` | A small value table (name, value, unit), e.g. the state points. *(v1.1)* |
| `<cs-log>` | `lines` | The last lines of the GUI console (hosted mode). *(v1.1)* |

Elements marked v1.1 come after v1. The prototype has the v1 elements except `widget` and
`options` of `cs-input`, `source="parametric"` and `realtime`.

### 5.3 Attribute bindings (any element, SVG included)

Custom elements do not work inside SVG, so the drawing is bound with attributes. Each takes a name
or an expression (§5.4).

| Attribute | Effect |
|---|---|
| `data-cs-value="T_2"` (+ `data-cs-unit`, `data-cs-decimals`, `data-cs-hide-unit`) | Text content = formatted value (an SVG `<text>`, a table cell) |
| `data-cs-fill`, `data-cs-stroke` + `data-cs-colormap="0:#2563eb,60:#f97316,90:#dc2626"` | Colour from a value (pipes coloured by temperature, a tank by its temperature) |
| `data-cs-show="x_b > 0 && x_b < 1"` | Visible while the expression is true (a valve open, a heater on, a flame) |
| `data-cs-rotate="theta"` + `data-cs-origin="62 540"` | Rotation in degrees (a crank, a needle, a fan) |
| `data-cs-translate="0, -L*2000"` | Translation in drawing units *(not in the prototype)* |
| `data-cs-scale-y="level/H"` + `data-cs-origin` | Scale (a liquid level, a bar) *(not in the prototype)* |
| `data-cs-flow="M_dot_r"` + `data-cs-flow-ref="0.16"` | Dashes moving along a path, period 1 s when the value equals the reference: flows come alive, also for steady models |

### 5.4 Expressions

A small language evaluated without `eval` (safe in a sandbox and under a content security policy):
numbers, variable names, `+ - * / ^`, comparisons, `&& || !`, parentheses and the functions `abs`,
`sqrt`, `exp`, `ln`, `log10`, `min`, `max`, `floor`, `round`, `clamp`, `sin`, `cos`, `tan`
(**degrees**, as CoolSolve), and `pi`. Names are **case-insensitive** (EES semantics), array
elements are written `h[3]`. An expression mixes integral columns (current frame) and solved
values: `100 * W_tot / Q_comb`.

### 5.5 JavaScript API

For what the declarative part cannot draw (a slider-crank mechanism, a custom gauge), the author
registers a function:

```js
CoolSolve.onFrame((ctx) => {
  ctx.get('theta')      // value at the current frame (integral column interpolated, or solved value)
  ctx.t                 // position of the player (integration variable), NaN for steady models
  ctx.playing           // true while the animation runs
  ctx.format(v, 2)      // number formatting of the runtime
  ctx.data              // the whole data document (§6)
});
```

It is called when data arrive (each solve) and at each animation frame. Also planned in v1:
`CoolSolve.value(name)`, `CoolSolve.onData(fn)`, `CoolSolve.player` (`play()`, `pause()`,
`seek(t)`), `CoolSolve.setInput(name, value)` and `CoolSolve.action(name)` (both no-ops, with a
problem notice, in standalone mode). Nothing else of the page is reachable from the host (§9.2).

### 5.6 Placement and resizing

A `<cs-stage width="1000" height="560">` is drawn in **design coordinates**: the author draws the
schematic in that coordinate system (inline SVG with the same `viewBox`, or an image of that size)
and gives each label its `x`, `y` and `anchor` (`top-left` default, `left`, `center`, `right`,
`top`, `bottom`, …). The runtime scales the whole stage — drawing, labels and input boxes — by one
factor computed from the width of its container (or from width and height with `fit="contain"`),
like a slide. A label therefore stays on the same point of the schematic at any window size, and
the text grows and shrinks with the drawing (D4). The rest of the page (cards, plots, side panels)
uses normal responsive CSS.

Checked with the prototype at 1280 px and at 375 px wide: the labels stay on their state points,
no horizontal scroll (`dashboards/heat_pump.html`).

### 5.7 Units

Values arrive in the units of the model (CoolSolve's fixed SI: Pa, J, W, kg, m, °C). A display unit
converts by a factor: `Pa → kPa, bar, MPa`; `J/kg → kJ/kg`; `J/kg-K → kJ/kg-K`; `J → kJ, MJ`;
`W → kW, MW`; `kg/s → g/s`; `kg → g, mg`; `m³ → cm³, L`; `m → mm, cm`; dimensionless `→ %`.
Anything else takes an explicit `scale`. A unit without conversion is a binding problem. Units
of the data come from the model (`[unit]` of the code, inferred units of the variable table).

### 5.8 Problems

A name that is not a variable, an input that is not an input, a saved plot that does not exist, an
expression that cannot be read, an exception in `onFrame`: each is listed in a **problem badge**
at the bottom right of the page (the element and the reason). A dashboard written for another
version of the model, or by an LLM that guessed a name, shows at once what to fix. The same rules
are checked statically by the checker (§10.4).

### 5.9 Versions

`<meta name="coolsolve-dashboard" content="1">` names the contract version. The runtime keeps
reading every published version; a dashboard without the tag is read as version 1 with a notice.
The data document has its own version (`schema`).

---

## 6. Data contract

One JSON document, built by one C++ function (`buildDashboardData()`, `src/dashboard.cpp`) used by
the server and the CLI (D7). The prototype snapshots `dashboards/*.results.js` follow it.

```json
{
  "schema": "coolsolve-dashboard-data/1",
  "model": { "name": "heat_pump_cycle_r22_30kw", "title": "R22 heat pump cycle delivering 30 kW",
             "libraryId": "CSL-0025" },
  "status": { "state": "solved", "message": "Solved in 41 ms", "solvedAt": "2026-11-20T12:00:00",
              "coolsolve": "CoolSolve 0.4.0", "canTryHarder": false },
  "variables": {
    "T_ev": { "value": 8, "units": "°C", "input": true, "line": 27,
              "description": "evaporation temperature [°C]" },
    "COP":  { "value": 4.4745, "units": "-", "input": false },
    "R$":   { "value": "R134a", "units": "", "input": true, "type": "string" }
  },
  "integral": { "var": "theta", "units": { "theta": "deg", "p": "Pa" },
                "columns": { "theta": [-360, -358], "p": [100000, 99810] } },
  "parametric": { "sweep": ["T_ev"], "runs": [ { "T_ev": 0, "COP": 3.9 } ] },
  "plots": { "ph": { "title": "R22 P-h diagram", "kind": "thermo", "fluid": "R22", "diagram": "Ph",
                     "x": { "name": "h", "units": "J/kg" }, "y": { "name": "P", "units": "Pa", "log": true },
                     "dome": { "x": [], "y": [] }, "path": { "x": "h", "y": "P", "close": true } } }
}
```

- `variables`: every variable of the result (arrays as their elements `h[1]`, …, strings with
  `type: "string"`), its units, whether it is an input (§7.1), the line of an input, and the
  comment of the defining line when there is one (useful to LLMs, §10.2).
- `integral`: present for `INTEGRAL` models; the columns of the trajectory (`$IntegralTable`, or
  what the solver writes without it). A standalone export decimates a table larger than a limit
  (default 5 000 rows, option `--dashboard-max-rows`) by keeping evenly spaced rows plus the first
  and the last; the GUI sends it whole.
- `parametric`: the last parametric study, every run with all its variables (`ParametricGridPoint`).
- `plots`: the saved plots with what they need to be drawn without CoolSolve (a thermodynamic
  diagram carries its dome, §9.10).
- Sources: the last solve of the session, or a `.sol` file (values, units, integral table) — which
  is how a library model shows its dashboard without solving (§4.4, §9.11).

---

## 7. Inputs

### 7.1 Which variables are inputs

The rule of the parametric table (R6): **a variable defined alone, by a constant, in the main
program** — `T_ev = 8`, `P_atm = 101151.7 [Pa]`, `Q_dot_23 = -30E3`. Extended in two ways, both
needed by the test models (D10):

- **array elements with a constant index**: `T[1] = 50`, `P[1] = 1000E3` (`CSL-0345`);
- **string constants**: `R$ = 'R134a'` (`CSL-0345`'s refrigerant, the radio button of the EES
  Diagram window), which also drive `$IF R$='…'` directives.

Not inputs: a variable computed by an expression (`N = 3000/60` is not an input; the checker says
so and suggests writing `N = 50`), a line inside a `FUNCTION`/`PROCEDURE`, a line of a `$IF` branch
that is removed, a variable defined twice. The detection moves out of the parse handler of
`src/server.cpp` into one function (`findModelInputs(const IR&)`, in a new
`src/model_inputs.cpp`) used by `/api/v1/parse`, the parametric study (which keeps sweeping numbers
only), the data builder, the checker and the prompt builder.

### 7.2 Widgets

Numbers: slider + number box (`min`, `max`, `step`; without them the runtime proposes ±50 % around
the value). Strings: `widget="select"` or `widget="radio"` with `options="R134a,R1234yf"`; without
`options` the checker suggests the values the model compares the string with (`$IF` directives,
`IF(X$ = '…')`).

### 7.3 How a change reaches the model

An input change **rewrites the line of the model** (D9): `T_ev = 8` becomes `T_ev = 5`, the comment
and the `[unit]` of the line are kept. The code stays the single source of truth: the editor shows
the new value, the bundle saves it, the CLI solves the same model.

- The rewrite reuses the parametric substitution (`replaceImposedValue`, moved from the lambda of
  the parametric handler to `src/model_inputs.cpp`, extended to strings).
- **The window that owns the editor applies it.** The editor's text is sent to the server only at
  solve and save, so the host sends the change to the main window (directly, or through the
  `BroadcastChannel` of §9.4 from a popped-out window); the main window rewrites its editor text
  with the same rule (`POST /api/v1/inputs` on the current text, which returns the new text and
  the changed lines) and applies it as an **undoable edit** (Monaco `pushEditOperations`): Ctrl+Z
  in the editor undoes a slider move. Without a main window (a popped-out dashboard whose main
  window was closed) the host rewrites the session's code on the server.
- **Reset inputs** returns to the values the inputs had when the dashboard was opened.

### 7.4 When to solve

Manual by default: the user moves the sliders and clicks *Solve* (D11). With
`<cs-dashboard autosolve="true">` a change solves after 400 ms without another change; a solve still
running is cancelled first. Recommended for small models only; the template comment says so.

---

## 8. Animation

### 8.1 What is played

The **integral table of the last solve**, in the browser: nothing is re-solved during the animation
(D12). "Live" means that the animation runs in real time on the page and that a new solve replaces
the table in place (the player keeps its position when the integration variable is the same).

### 8.2 The player

`<cs-player>`: ⏮ start, ◀◀ previous row, ▶/❚❚ play-pause, ▶▶ next row, a scrubber, the position
(`theta = -251.2 deg`), the speed (×0.25 to ×4) and *loop*. Keys: space, ←/→, Home. The whole
range lasts `duration` seconds at ×1 (default 10 s); `realtime` plays one second of a `time` variable
in one second. Between two rows the values are **interpolated linearly**, so the animation is
smooth at 60 frames per second whatever the number of rows (the DHW tank of §11 has 18 000; an
engine at 2° per row has 271).

### 8.3 What moves

Everything bound to a name follows the player: `cs-value`, the `data-cs-*` bindings (colours,
visibility, rotations, translations, flows), the cursors of the plots (`cursor`: a marker on a p–V
loop, a vertical line on a time plot), and the `onFrame` functions (mechanisms). Steady values
stay as solved. See `dashboards/engine.html`: crank, rod, piston and intake valve follow the crank
angle, the gas colour follows its temperature, the spark and the flame show at the right angles,
the intake flow moves while the valve is open, the p–V loop and the pressure and temperature
curves carry the cursor.

### 8.4 Later

- **Parametric player** (`source="parametric"`): step through the runs of a 1-D sweep (each run
  keeps all its variables), so that the schematic and the cycle on the P-h diagram move with the
  swept evaporation temperature. Cheap, and very telling in lectures.
- **Follow the solve**: stream the rows of the trajectory while the integration runs (SSE) and
  play them as they come. Needs a solver callback per accepted step; worth it for long runs only.

---

## 9. Architecture

### 9.1 Components

```
gui/src/dashboard/            runtime (TypeScript) → gui/dist/dashboard/coolsolve-dashboard.js
  runtime.ts, elements/*.ts,    (second Vite entry; embedded in the binary with the GUI)
  expr.ts, units.ts, player.ts
  AUTHORING (from docs/dashboards.md, §10.2) and templates/*.html
gui/src/plots/figures.ts      Plotly figures of thermo / parametric / integral plots, shared by
                              ThermoDiagram, the saved-plot tabs and the runtime
gui/src/components/
  DashboardTab.tsx            host: sandboxed iframe + message bridge (§9.2, §9.3)
  DashboardMenu.tsx           toolbar drop-down, Add dashboard dialog (§4.3)
  SavedPlotTab.tsx            a saved plot as a bottom tab (§4.1)
include/coolsolve/dashboard.h, src/dashboard.cpp
                              data builder (§6), checker (§10.4), prompt builder (§10.2),
                              standalone exporter (§9.8)
src/model_inputs.cpp          inputs of a model and the rewrite of their line (§7)
src/server.cpp                routes (§9.5), session state, bundle (§9.7)
main.cpp                      --dashboard-prompt / --dashboard-check / --dashboard-export (§9.9)
```

The runtime is written in TypeScript and built by Vite so that it shares the figure code with the
GUI (D5). The prototype of this folder is plain JavaScript only so that the examples run without
a build.

### 9.2 Security

A dashboard is code: written by an LLM, downloaded with someone's bundle, or opened from the
library on the public server. It never runs with the rights of the GUI (D8):

- The host shows it in an `<iframe sandbox="allow-scripts">` (no `allow-same-origin`): the page has
  an opaque origin, cannot read the GUI's cookies or storage, cannot call `/api/` with the user's
  session, cannot touch the editor.
- It is served with a content security policy: `default-src 'none'; script-src 'self'
  'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; font-src 'self'
  data:; connect-src 'none'; form-action 'none'`, plus the CSP `sandbox allow-scripts` directive so
  that the page is sandboxed even when its URL is opened directly. No network: a dashboard cannot
  send the model anywhere.
- It talks to the host with `postMessage` only (§9.3); the host checks the sender (`event.source`
  is its iframe), accepts the listed messages only, and checks every input name against the
  inputs of the model before rewriting a line.
- A standalone export is a file on the user's disk; the export dialog says that it contains the
  model's results and, when chosen, its code.

### 9.3 Message protocol

| From | Message | Meaning |
|---|---|---|
| dashboard | `{source: 'coolsolve-dashboard', version: 1, type: 'ready'}` | Loaded; send the data |
| host | `{source: 'coolsolve-host', type: 'data', data, theme}` | Data (§6) after load and after each solve |
| host | `{source: 'coolsolve-host', type: 'status', status}` | Solve progress, failure, *Try Harder* available |
| dashboard | `{… type: 'action', action}` | `solve`, `update-guesses`, `stop`, `reset-inputs` |
| dashboard | `{… type: 'set-input', name, value}` | An input changed (§7.3) |

Without an answer to `ready` within 400 ms the runtime takes the standalone path. Checked with the
prototype: a test page that answers `ready` enables the inputs and buttons, applies the theme, and
receives `action` and `set-input` messages.

### 9.4 Windows

- *Pop out* opens `/?dashboard=<name>` in a new window that uses the session of the main window:
  the cookie session as it is, or a tab session handed over through `localStorage` as the include
  links do (`coolsolve_include_opener`), then adopted (`adoptTabSession`); the parameter is
  consumed and removed from the address. The window shows the dashboard alone, full size.
- The windows of one session keep in step through a `BroadcastChannel('coolsolve-<session>')`:
  `solved` (re-read the data), `code-changed`, `set-input` (applied by the window that owns the
  editor). Same browser only; a second device would need server events (`GET
  /api/v1/session/events`, SSE), listed in §14.

### 9.5 Server routes

| Method | Route | Purpose |
|---|---|---|
| `GET` | `/api/v1/dashboards` | `[{name, title, bytes, modified}]` of the session |
| `GET` / `PUT` / `DELETE` | `/api/v1/dashboards/{name}` | Source of a dashboard (`text/html`); `PUT` returns the checker report |
| `PUT` / `DELETE` | `/api/v1/dashboards/assets/{file}` | An image of the dashboards (`svg`, `png`, `jpg`, `webp`) |
| `POST` | `/api/v1/dashboards/check` | Checker report of an HTML text against the current model (§10.4) |
| `POST` | `/api/v1/dashboards/prompt` | The prompt for an assistant (§10.2); body: request, options |
| `GET` | `/api/v1/dashboards/{name}/standalone` | The single-file export (§9.8) |
| `GET` | `/api/v1/dashboard-data` | The data document of the session (§6) |
| `POST` | `/api/v1/inputs` | `{eescode, values: {T_ev: 5}}` → `{eescode, changed: [{name, line, old, new}]}` (§7.3) |
| `GET` / `PUT` / `DELETE` | `/api/v1/plots`, `/api/v1/plots/{id}` | Saved plots (§4.1) |
| `GET` | `/d/{session}/{name}.html`, `/d/{session}/assets/{file}`, `/d/{session}/coolsolve-dashboard.js`, `/d/{session}/plotly.js` | The dashboard document and its files, with the policy of §9.2 |

The session ID is in the path so that the relative URLs of a dashboard (`assets/x.svg`,
`coolsolve-dashboard.js`) resolve as in the bundle folder; it is the same capability as `?sid=`.
New session state: `dashboards` (name → HTML), `dashboardAssets` (file → bytes), `savedPlots`
(JSON). *New*, *Open* and the library loader reset them as they reset the lookup tables;
`showSessionModel()` (`gui/src/utils/sessionModel.ts`) reads them back.

### 9.6 Names

Dashboard names: `[A-Za-z0-9_-]+`, as table names. Saved-plot ids: the same, generated from the
title, unique in the session.

### 9.7 Bundle and model folder

```
my_model/
├── my_model.eescode, .initials, .sol, coolsolve.conf, my_model-<table>.csv, …
├── plots.json                    saved plots (§4.1)
└── dashboards/
    ├── main.html                 a dashboard (source)
    ├── animation.html
    ├── assets/schematic.svg      images used by the dashboards
    ├── coolsolve-dashboard.js    written at export, ignored at import
    └── results.js                snapshot of the last results, written at export, ignored at import
```

- `createZipBundle` writes `plots.json`, `dashboards/*.html`, `dashboards/assets/*`, and, when the
  model has dashboards, the runtime and `results.js` (`window.COOLSOLVE_DATA = {…}`, a script and
  not JSON so that `file://` pages can read it); `extractZipBundle` reads the first three and
  ignores the generated two. Round trip: download → upload → same dashboards, same saved plots.
- A model opened from disk (`POST /api/v1/files/open`) reads `plots.json` and `dashboards/` of its
  folder, like its lookup tables.
- Library models may carry the same layout (§9.11).

### 9.8 Standalone export

One HTML file: the dashboard, with its `<script src="coolsolve-dashboard.js">` replaced by the
inlined runtime, Plotly inlined (the cartesian bundle the GUI already uses, about 1.4 MB, only when
the dashboard has a `cs-plot`), the data inlined as `<script type="application/json"
id="coolsolve-data">`, and the images as data URIs. Deterministic output (same model, same file),
so that exports can be committed or compared. Written by the C++ exporter (server and CLI); the
runtime and Plotly come from the embedded GUI assets, so the export needs a GUI build (D5); a
CLI-only build says so.

### 9.9 Command line

```bash
coolsolve --dashboard-prompt model.eescode [--request "…"] [--no-code] > prompt.md
coolsolve --dashboard-check dashboards/main.html model.eescode [--format json]
coolsolve --dashboard-export dashboards/main.html model.eescode [-o main_standalone.html] [--from-sol]
```

- `--dashboard-prompt` prints the prompt of §10.2 for the model (values from the `.sol` file when
  there is one; the model is not solved).
- `--dashboard-check` parses the model (no solve) and checks the dashboard (§10.4): Markdown or
  JSON report, exit code 1 when there are errors.
- `--dashboard-export` solves the model with the usual pipeline (or reads its `.sol` with
  `--from-sol`) and writes the standalone file (§9.8).

### 9.10 Plots

- A saved plot is a spec: `{id, title, kind: 'thermo' | 'parametric' | 'integral', …}` with the
  settings of the tab it was saved from (fluid, diagram, overlays, array path; sweep axes; columns).
- Its Plotly figure is built by one module, `gui/src/plots/figures.ts`, used by the Diagrams tab,
  the saved-plot tabs and the runtime: the three draw the same picture.
- What needs CoolProp is computed by the server and travels in the data (§6): the saturation dome
  (the computation of `POST /api/v1/coolprop/saturation` moved to a function the CLI exporter also
  calls). When the CoolPlot integration lands ([coolplot_integration.md](../coolplot_integration.md))
  a thermo plot may carry the SVG scene of coolplot-core instead, which a dashboard shows as is.

### 9.11 Library

- Allow-list: `dashboards/*.html` and `dashboards/assets/*.{svg,png,jpg,jpeg,webp}` of a model
  folder, and `plots.json`, in `cmake/embed_library.cmake` **and** `EMBED_FILES` of the library's
  `tools/build_index.py` (they mirror each other); `LibraryStore::loadModel` loads them into the
  session; `model.json` lists them (`"dashboards": ["main"]`).
- The library dashboard (`library.html`) shows a dashboard badge, and *View dashboard* on the model
  page: `GET /library/<model path>/dashboards/<name>.html` served as a standalone snapshot built
  from the embedded `.sol` (values and `# IntegralTable`), with no solve. The public site gains a
  gallery of animated models that anyone can explore.
- The library workflow (`docs/model_workflow.md` of the library) gains a step: a dashboard is
  checked with `--dashboard-check` by the library tests, and its screenshot may be the README figure.

---

## 10. LLM workflow

### 10.1 Why

Laying out a schematic with labels, writing SVG, an animation of a mechanism: hard for most users,
easy for today's LLMs, provided they know the contract and the names of the model. The contract
(§5) is built for that: standard HTML, few concepts, names checked, errors visible on the page.
The GUI makes the LLM route the obvious one (first route of *Add dashboard*), without depending on
any provider.

### 10.2 The prompt

Generated by the C++ prompt builder (GUI dialog and `--dashboard-prompt` print the same text), in
Markdown:

1. **Role and goal**: "Write one HTML file, a CoolSolve dashboard for the model below, that …" +
   the user's request.
2. **The contract**: the authoring reference, taken from the user documentation
   `docs/dashboards.md` between `<!-- prompt:begin -->` and `<!-- prompt:end -->` markers, embedded
   in every build (one source for the users and for the LLMs, D16).
3. **Rules**: one file; inline SVG for drawings (no external images unless provided); no network
   (no CDN, no fonts from the web — blocked in the GUI); only names listed below; inputs only from
   the input list; place labels on a `cs-stage` in design coordinates; keep the dashboard usable at
   375 px wide; comment the parts a person may want to adjust.
4. **The model**: name, title, kind (steady / dynamic), the variables table (name, units, value,
   input or not, description from the line comment), the arrays, the integral columns and the
   integration range, the saved plots (id, title, kind), the parametric study if any, and — when
   chosen — the model code.
5. **Answer format**: the HTML file only, in one fenced block.

A sketch of the generated text for `CSL-0025` is in [§18](#18-the-example-files); the two example
dashboards are what such a prompt should produce.

### 10.3 Agents with the binary

An agent that can run commands (Claude Code, Codex, Aider, …) gets a loop with feedback, which
gives better dashboards than a one-shot answer:

```bash
coolsolve --dashboard-prompt model.eescode --request "…" > prompt.md     # instructions + names
# … the agent writes dashboards/main.html …
coolsolve --dashboard-check dashboards/main.html model.eescode           # fix until no error
coolsolve --dashboard-export dashboards/main.html model.eescode -o out.html
# … the agent opens out.html in a headless browser, looks at the screenshot, adjusts …
```

The prompt ends with these commands when it is generated by the CLI, so an agent knows them.
A Claude Code skill (`coolsolve-dashboard`) packaging the loop is listed in §14.

### 10.4 The checker

Static checks of a dashboard against the model (parse only), shared by the GUI (*Add dashboard*,
`PUT`) and the CLI:

| Level | Check |
|---|---|
| error | a name of `cs-value`, `cs-input`, `data-cs-*`, an expression, `ctx.get('…')` / `CoolSolve.value('…')` in a script that is not a variable of the model |
| error | a `cs-input` on a variable that is not an input (§7.1), with the reason and the fix (`N = 3000/60` → write a constant) |
| error | a `cs-plot ref` that is not a saved plot; `x`/`y` that are not integral columns (with the `$IntegralTable` line to add) |
| warning | a `cs-player` or a `cursor` for a model without `INTEGRAL`; an unknown `cs-*` element or attribute (typo) |
| warning | a URL to the network (blocked in the GUI); a missing asset; a file larger than 2 MB |
| info | stage children without `x`/`y`; a missing `<meta name="coolsolve-dashboard">` |

The report has a *fix-up prompt* form: the list of problems with the valid names close to each
wrong one, ready to paste back to the assistant.

### 10.5 Evaluation

Phase 6 tries the prompt on the test models (§11) with two or three assistants, records the
prompts, the first answers and the number of fix-up rounds, and improves the contract text where
they stumble. The prompts and outcomes are kept with the plan.

---

## 11. Test models

Chosen in the library to cover the features, from small to large, with the three models of the
register row `CS-GAP-DIAGRAM-WINDOW` (the EES originals had a Diagram window or an animation, so
the dashboards restore what the import dropped).

| Model | Kind | Exercises | Dashboard |
|---|---|---|---|
| `CSL-0025` *heat_pump_cycle_r22_30kw* (55 eq.) | steady, R22 | values on a schematic, numeric inputs, *Solve* / *Update guesses*, a saved P-h diagram from the arrays `P[i]`, `h[i]`, pipes coloured by temperature, flows moving with the mass flow, a warning limit | `dashboards/heat_pump.html` (prototype); to become the library's `dashboards/main.html` |
| `CSL-0345` *r134a_isenthalpic_throttle* (16 eq., variant R1234yf) | steady | the input-rule extension: **string input** (radio R134a / R1234yf) and **array-element inputs** (`T[1]`, `P[1]`); parity with the radio button and edit boxes of the EES Diagram window | throttle sketch, inlet/outlet states, T-h saved plot |
| `CSL-0339` *gas_spring_oscillator* (28 eq.) | dynamic (time, 5 s) | the EES animation of the original (`Ball.top`): `data-cs-translate`, colour by gas temperature, `realtime` player, L(t) and P(t) with cursors. The file has no `$IntegralTable`: the checker must tell which line to add for P and T | ball in a tube over a bottle |
| `CSL-0024` *single_cylinder_engine_weibe* (72 eq.) | dynamic (crank angle) | a mechanism drawn by `onFrame`, declarative show / rotate / fill / flow, p–V loop on log axes with a cursor, a band for the combustion window | `dashboards/engine.html` (prototype), on a variant integrating from −360° to +180° (below) |
| `CSL-0009` *dhw_tank_dynamic* (20 eq.) | dynamic, 18 000 rows | a large integral table (decimation of exports, interpolation speed), schedules shown with `data-cs-show` (heater on, tap draw), tank colour by temperature, 5 h played in 15 s | tank, heater, tap, T(t), heat flows |
| `CSL-0019` *orc_simple_r245fa* (182 eq.) | steady, large | about forty labels on a P&ID-like schematic, speed of the data builder and of the page, a T-s saved plot | stretch goal |

**The engine variant.** The library model integrates the intake stroke only (−360° to −180°), which
animates half a crank turn without combustion. The prototype data come from the same file with the
five `INTEGRAL` calls and `theta_2` set to +180° (intake, compression, combustion, expansion): it
solves with RK4 at 1° in about 10 s, peak 90 bar at θ = 10°, work 1.33 kJ for 2.67 kJ released.
Proposed as a library variant `single_cylinder_engine_weibe_cycle_coolsolve.eescode` with its
`.sol`. While building the dashboard, the animation showed something to look at before: during the
intake the charge cools to 265 K although it mixes air at 300 K with residual gas at 1400 K — the
energy brought by the incoming air is written `M_dot·C_v·(T_out − T)`, without the flow work
`M_dot·r·T_out` (the enthalpy). To be examined by the library workflow when the variant is made;
it is a fair advertisement for dashboards that a picture shows it in seconds.

---

## 12. Phases

Each phase ends with its tests green and its documentation, following `docs/contributing.md` §2.
Sizes are rough.

| Phase | Content | Done when |
|---|---|---|
| **0 — Plan and prototype** *(done)* | This plan; prototype runtime and two examples (§18) | — |
| **1 — Contract, runtime, standalone export** | `gui/src/dashboard/` (port of the prototype, Plotly charts, ~1 500 lines TS); `src/dashboard.cpp`: data builder (live result or `.sol`), checker, exporter (~1 200 lines C++); `--dashboard-check`, `--dashboard-export`; `docs/dashboards.md` (contract reference) | An agent makes a standalone dashboard of a library model with the CLI alone |
| **2 — GUI host and bundle** | Session state, routes (§9.5), bundle and disk folder (§9.7), toolbar menu, bottom tabs, sandboxed host and message bridge, buttons, live updates, pop-out windows and `BroadcastChannel`, presentation mode | A dashboard of a bundle opens in the GUI, solves from its button, follows a solve made in the editor, and survives a download/upload round trip |
| **3 — Inputs** | `src/model_inputs.cpp` (detection with array elements and strings, rewrite), parametric study on the shared detection, `/api/v1/inputs`, undoable editor edits, `autosolve`, *Reset inputs*, string widgets | `CSL-0345`: the radio button switches the refrigerant, the editor shows it, Ctrl+Z undoes it |
| **4 — Saved plots** | *Save* on the three plot tabs, `plots.json`, saved-plot tabs, `gui/src/plots/figures.ts`, dome in the data, `cs-plot ref` | The P-h diagram of `CSL-0025` is a tab and is in its dashboard, live and in the standalone file |
| **5 — Animation** | Player complete (`realtime`, events, `CoolSolve.player`), `data-cs-translate` / `data-cs-scale-y`, cursors on all plot kinds, parametric player | `CSL-0339`, `CSL-0024` (variant) and `CSL-0009` animate smoothly |
| **6 — LLM workflow** | Prompt builder, `--dashboard-prompt`, *Add dashboard* dialog (three routes), templates, fix-up prompt, evaluation (§10.5) | A user without HTML knowledge makes the dashboard of a test model through an assistant in a few minutes |
| **7 — Library** | Allow-list on both sides, loader, `model.json`, *View dashboard* from the `.sol`, dashboards of the test models (library workflow) | The test models have their dashboards on the public library pages |
| **8 — Ideas** | A selection from §14 | — |

Phase 6 needs only Phase 1 for the CLI loop and Phase 2 for the dialog: it can run in parallel with
Phases 3–5. Phase 7 needs Phases 1–2.

---

## 13. Integration checklist

`docs/contributing.md` §3, item by item.

| Item | Status | What |
|---|---|---|
| 3.1 Language and parser | not applicable | No new syntax. The input detection (§7.1) works on the IR. |
| 3.2 Solvers | not applicable | — (a callback per accepted integration step only if *follow the solve* is done, §8.4) |
| 3.3 `coolsolve.conf` | not applicable | Dashboard settings live in the dashboard (`<cs-dashboard>`), not in the solver configuration |
| 3.4 Config editor | not applicable | — |
| 3.5 GUI frontend | needed | `types.ts` (dashboards, saved plots, data), `client.ts` (routes through `request()`), `modelStore.ts` / `uiStore.ts` (dashboards, saved plots, bottom tabs), `sessionModel.ts` (read them back), `App.tsx` (bottom tabs), `Toolbar.tsx` (menu), `ThermoDiagram.tsx` / `ParametricStudy.tsx` / `IntegralTable.tsx` (*Save*), new components of §9.1, `CodeEditor.tsx` (undoable input edits) |
| 3.6 Server / REST | needed | Routes of §9.5 through `getSession()`; session IDs of `/d/{session}/…` checked with `isValidSessionId()`; `tests/test_dashboard_api.py` (`--auto`) |
| 3.7 ZIP bundle | needed | §9.7 in `createZipBundle` / `extractZipBundle`; round-trip test; layout documented in `contributing.md` §3.7 and `gui.md` §7.2 |
| 3.8 Debug mode | not applicable | — |
| 3.9 Console | needed | Checker problems and input rewrites reported in the console with the usual markers (`[Warning]`, `ERROR`) |
| 3.10 LaTeX report | optional | Saved plots in the report (they are Plotly charts; `exportPlots.ts` already captures visible charts) |
| 3.11 Examples | needed | `examples/dashboards/` with a small model and a dashboard for the CLI tests (`--dashboard-check`, `--dashboard-export`) |
| 3.12 README and docs | needed | §17 |
| 3.13 Build system | needed | Second Vite entry for the runtime; the runtime and `docs/dashboards.md` embedded (the latter in every build, for the prompt); library allow-list on both sides |
| Tests | needed | Catch2: data builder (live and `.sol`), input detection and rewrite (numbers, array elements, strings, comments and units kept), checker (each rule, negative cases), exporter (deterministic, decimation). Headless Chromium (pattern of `tests/test_gui_shortcuts.py`): the examples render without problem badge; a label keeps its place on the schematic from 1280 to 375 px; the player advances and a cursor moves; a hosted dashboard's button solves and an input rewrites the editor; the sandbox cannot `fetch('/api/v1/files/eescode')` |

---

## 14. Ideas

Ranked by value for effort; the first ones are cheap once the core exists.

1. **What-if comparison.** *Pin as baseline*: the dashboard keeps the current results as a
   reference; values show their difference (`COP 4.47 ▲ +0.21`), plots draw the baseline in grey
   under the new curve. The most useful thing an engineer does with a dashboard, rarely offered.
2. **Flows that move** (`data-cs-flow`, in the prototype): pipes with dashes moving at a speed
   proportional to the mass flow; with colours by temperature (`data-cs-stroke`), a steady cycle
   looks alive.
3. **Gallery of the library** (§9.11): animated, explorable models on the public site, no install,
   no solve. Each one is an invitation to open CoolSolve.
4. **Presentation mode and second screen** (§4.2, §9.4): the dashboard full screen on the projector,
   the code on the laptop; a lecturer moves a slider and solves in front of the students.
5. **Parametric player** (§8.4): watch the cycle move on the P-h diagram as the evaporation
   temperature is swept.
6. **More widgets** for LLMs to use: `cs-gauge` (dial with zones), `cs-bar` (horizontal bar against a
   limit), `cs-check` (a constraint shown green or red, e.g. "discharge temperature < 120 °C"),
   `cs-sankey` (energy flows of a plant from a list of `Q_dot`/`W_dot`: a Sankey diagram of the cycle,
   drawn in SVG by the runtime), `cs-table`.
7. **Exports from the page**: PNG of the dashboard; a video (WebM) of an animation, recorded frame by
   frame from the SVG; a dashboard snapshot in the LaTeX report.
8. **Claude Code skill / agent recipe** `coolsolve-dashboard` (§10.3) published with the library:
   "make a dashboard for this model" in one sentence for users of coding agents.
9. **Starter dashboard generated without an LLM**: from the model alone, a page with the inputs as
   sliders, the outputs as a table, the saved plots and, for an `INTEGRAL` model, a player and one
   plot per column. A useful default that also serves as a template.
10. **Import of the EES Diagram window**: decode the text items, edit boxes and their positions in
    the EES file (`ees_extract.py`) into a dashboard skeleton. Research: the binary layout of these
    objects is not known yet ([ees_import.md](../ees_import.md)).
11. **Server events** (§9.4) so that a dashboard on a tablet follows a solve made on the laptop.
12. **Sharing**: on the public server, *Share* publishes a standalone snapshot at a short URL
    (read-only, expiring), for a forum post or a course page.
13. **Display unit switch** on the page (SI / bar-kJ / °C-K) applied to every `cs-value` without unit.

---

## 15. Decisions

Taken while hesitating, without asking (R14). Each can be revisited; the alternatives are kept.

| # | Decision | Alternatives | Why |
|---|---|---|---|
| D1 | The plan and its examples are in `docs/plans/` | A flat `docs/` file (served); a folder outside `docs/` | Subfolders of `docs/` are not embedded (`embed_assets.cmake` globs the files of `docs/` only) and the release server has no filesystem fallback: not served, yet in `docs/` as asked. `docs/plans/README.md` says so. |
| D2 | "Several diagrams" = **saved plots**: the Diagrams tab stays the scratch diagram; *Save* adds a bottom tab. Extended to the parametric and integral plots | Several Diagrams tabs on the right; saved plots only inside dashboards | The request says "keep the current window, save it, new tab at the bottom"; the same gesture on the three plot kinds is consistent; dashboards reuse saved plots by id |
| D3 | A dashboard is **a plain HTML file** with custom elements, `data-cs-*` attributes and optional JavaScript | A JSON layout drawn by the GUI; a GUI editor (EES-like) | LLMs write HTML/SVG/JS very well; users edit it in any editor; anything the web can draw is possible; a standalone file is natural. The request excludes editing in the GUI |
| D4 | Labels are placed on a **stage in design coordinates scaled as one piece** | Percent positions; SVG `foreignObject`; CSS container units | Exact placement at any size with text that scales with the drawing; works the same for inline SVG and images; `foreignObject` has browser quirks for inputs; percent positions keep the font size and drift with the aspect ratio. Checked at 1280 and 375 px |
| D5 | The runtime is **TypeScript in `gui/src/dashboard/`, built by Vite**; the standalone export needs a GUI build | Hand-written JavaScript embedded in every build | Shares the figure code with the GUI (no duplication, §1 of the contributing guide); the distributed binaries (installer, server) are GUI builds. The prototype is plain JS only to run without a build |
| D6 | Plots use **Plotly**, as the GUI | A small SVG chart library (as in the prototype) | Same pictures as the GUI tabs, zoom and hover for free, log axes, contours; 1.4 MB inlined only in exports that have a plot |
| D7 | The data document is built by **one C++ function** for the server and the CLI | Assembled in the GUI from its store | One builder for the GUI, the CLI export and the library snapshots (from `.sol`) |
| D8 | Dashboards run in a **sandboxed iframe** and speak `postMessage` only — library dashboards too | Same-origin iframe calling the REST API | Dashboards are code from LLMs and from other people's bundles, on a public server too; an opaque origin and `connect-src 'none'` make a malicious page harmless |
| D9 | An input change **rewrites the line of the model**, as an undoable editor edit | An override layer next to the code (EES keeps Diagram-window values apart from the equations) | The code stays the single source of truth for the editor, the bundle and the CLI; reuses the parametric substitution; no new state to save |
| D10 | Inputs include **constant-index array elements and string constants** | Scalars and numbers only, as the parametric table today | Needed by `CSL-0345` (`T[1]`, `P[1]`, `R$`) and by EES models whose Diagram window set strings for `$IF` directives; the parametric table gets array elements from the same detection |
| D11 | **Manual solve** by default; `autosolve` per dashboard | Always solve on change | Many library models take seconds; a slider that solves at each move would freeze the page |
| D12 | Animation is a **playback of the integral table** in the browser, interpolated; *follow the solve* later | Re-solving per frame; streaming during the solve | "Live, but using the results of the integral table" (the request); smooth, instant, works in standalone files |
| D13 | The bundle carries the runtime and `results.js` (generated, ignored at import), not Plotly | Nothing generated in the bundle; Plotly in every bundle | The unzipped folder opens offline from `file://` for a few kB; 1.4 MB in every bundle is too much, and the standalone export covers offline plots |
| D14 | Dashboards and saved plots are **bottom-panel tabs**, with pop-out and presentation | A right-panel tab; a separate page only | The request places saved diagrams at the bottom; one place for both; the panel is resizable and maximisable, and pop-out gives a full window |
| D15 | Buttons: `solve` (with *Try Harder*), `update-guesses`, `stop`, `reset-inputs` | Also *Save*, *Open*, *Debug solve* | The actions that make sense while exploring; file actions stay in the toolbar (a dashboard must not trigger downloads) |
| D16 | The GUI **does not call an LLM**: it generates the prompt, the user pastes it; the contract text has one source, `docs/dashboards.md` (markers), embedded in every build | Built-in LLM calls with an API key; a separate authoring guide | No key, no cost, offline, provider-neutral; one text for users and LLMs cannot drift |
| D17 | CLI options `--dashboard-prompt`, `--dashboard-check`, `--dashboard-export` | Sub-commands (`coolsolve dashboard check`) | Consistent with the existing long options of `main.cpp` |
| D18 | Display units by **scale factors only** in v1 (and `scale`) | Offsets (°C ↔ K, °F) | CoolSolve's units are fixed SI with °C; factors cover pressures, energies, powers, flows; offsets wait for a real need |
| D19 | Names are **case-insensitive**, trigonometric functions in **degrees** | JavaScript conventions | EES/CoolSolve semantics: `T_ev` and `t_ev` are one variable; `sin(30) = 0.5` in the model and on the dashboard |
| D20 | The test models are `CSL-0025`, `CSL-0345`, `CSL-0339`, `CSL-0024` (+ variant), `CSL-0009`, and `CSL-0019` as a stretch | Others (building `CSL-0042`, absorption chillers) | Small to large, steady and dynamic, the three Diagram-window cases of the register, a mechanism, a long table |
| D21 | The prototype's examples use **real CoolSolve results** (library `.sol`, and a run of the engine variant), copied into `*.results.js` | Synthetic data | Shows what a real dashboard looks like and already revealed a modelling question (§11) |

---

## 16. Risks and open questions

- **CSP and opaque origins.** `'self'` in a policy of a sandboxed document must still match the
  server's URLs for the runtime and assets (it is the URL's origin in CSP 3); to verify on Chromium,
  Firefox and Safari in Phase 2. Fallback: inline the runtime into the served document.
- **Plotly in a sandbox.** Plotly needs no network or storage; its hover and zoom work in sandboxed
  iframes (to confirm in Phase 1).
- **Editor ownership** (§7.3). Several windows of one session could edit at once; the rule "the
  window with the editor applies input changes" must be enforced (one owner, announced on the
  `BroadcastChannel`).
- **Large tables.** 18 000 rows × 10 columns is ~2 MB of JSON: fine for `postMessage`, decimated for
  exports. Models with 100 000 rows would need decimation in the GUI too.
- **Contract stability.** Once library dashboards exist, the v1 contract must stay readable; new
  features are additions. The prototype fixes the names now so that LLM prompts tested in Phase 6
  stay valid.
- **LLM quality varies.** The checker, the problem badge and the fix-up prompt are the safety net;
  the evaluation of §10.5 tells whether the contract text needs examples inside the prompt (it may:
  a short example per element is cheap in tokens).
- **Online server load.** *View dashboard* of the library builds snapshots from `.sol` (no solve);
  hosted dashboards solve like the GUI does. No new load pattern.
- **Open:** should a dashboard be able to show *another* model's results (comparison of two cycles)?
  Not in v1: one dashboard, one model.

---

## 17. Documentation to update when it lands

- `README.md`: *Features* (dashboards and saved plots), *Command Line Options* (three options),
  *File Formats* (`dashboards/`, `plots.json`), *Project Structure* (new files), *Documentation*
  table (`docs/dashboards.md`).
- New **`docs/dashboards.md`**: user guide and contract reference (the prompt is taken from it,
  §10.2), with screenshots; in the sidebar of `docs/docs.html`.
- `docs/gui.md`: routes (§4), session state (§5), layout and tabs (§6), bundle (§7.2), library
  (§7.3); remove *Diagram state persistence* from §14.4 once saved plots exist.
- `docs/contributing.md`: §3.7 (bundle layout), §3.5 (new components).
- `docs/model_library_support.md`: a feature row `CS-FEAT-DASHBOARD` in §3; the resolution path of
  `CS-GAP-DIAGRAM-WINDOW` and of case (a) of `CS-GAP-IF-NO-CONSTANT` (string inputs on a dashboard).
- Library: `docs/model_workflow.md` (dashboards of a model, checker in the tests), `EMBED_FILES`.
- `docs/versions.md` at the release.

---

## 18. The example files

In [`dashboards/`](dashboards/), next to this plan (not served either):

| File | What |
|---|---|
| `coolsolve-dashboard.js` | Prototype runtime (plain JavaScript, no dependency): stage scaling, `cs-value`, `cs-input`, `cs-button`, `cs-status`, `cs-plot` (small SVG charts instead of Plotly), `cs-player`, the `data-cs-*` bindings, the expression language, units, problem badge, layout helper (<kbd>L</kbd>), and the host protocol of §9.3 |
| `engine.html` + `engine.results.js` | `CSL-0024` animated: slider-crank mechanism, valves, spark, flame, gas coloured by temperature, p–V loop and p(θ), T(θ) with cursors, player (the "animated internal combustion engine" of the request) |
| `heat_pump.html` + `heat_pump.results.js` | `CSL-0025`: values on the schematic at the four states, COP, pipes coloured by temperature with moving flows, inputs and buttons (disabled standalone), saved P-h diagram |

**Opening them.** Double-click the `.html` file (it reads its `*.results.js` from the same folder,
`file://` is fine), or serve the folder: `python3 -m http.server -d docs/plans/dashboards`. Keys:
space (play/pause), ←/→ (step), <kbd>L</kbd> (layout helper; drag a label, its new tag is copied).

**Where the data come from.** `heat_pump.results.js`: the `.sol` baseline of `CSL-0025` and the R22
dome of `POST /api/v1/coolprop/saturation`. `engine.results.js`: a CLI run of `CSL-0024` with the
integration window extended to +180° (§11). Inputs were detected with the rule of §7.1, units taken
from the comments of the model. They follow the data contract of §6.

**Sketch of a generated prompt** (`CSL-0025`, abbreviated):

````markdown
You are writing a CoolSolve dashboard: one HTML file that shows the results of the model below.
Request: the cycle on a schematic with T, P and h at each state point, the COP in big, sliders for
T_ev, T_cond, DELTA_T_sc and epsilon_s_cp, a Solve button, and the saved P-h diagram "ph".

## Dashboard contract (version 1)
… (from docs/dashboards.md: file skeleton, elements, data-cs-* bindings, expressions,
CoolSolve.onFrame, cs-stage placement, units) …

## Rules
- Answer with the HTML file only, in one ```html block. Inline SVG for drawings; no network.
- Use only the variable names listed below (case does not matter). Inputs only from the inputs list.
- Put the schematic and its labels on one <cs-stage> in design coordinates (e.g. 1000 × 560).
- The page must work from 375 px wide; comment what a person may want to adjust.

## Model: heat_pump_cycle_r22_30kw — R22 heat pump cycle delivering 30 kW (steady)
| name | units | value | input | description |
|---|---|---|---|---|
| T_ev | °C | 8 | yes | evaporation temperature [°C] |
| T_cond | °C | 60 | yes | condensation temperature [°C] |
| COP | - | 4.4745 | | coefficient of performance [-] |
| h[1..4], P[1..4], T[1..4], s[1..4] | J/kg, Pa, °C, J/kg-K | … | | state points of the cycle |
| … 55 variables … |
Saved plots: ph (thermo, R22, P-h, path h[i], P[i], closed)
Integral table: none.

## Check your file
coolsolve --dashboard-check dashboards/main.html heat_pump_cycle_r22_30kw.eescode
````
