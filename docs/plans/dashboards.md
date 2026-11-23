# Dashboards plan

**Status:** proposal, 2026-11-20, not started (Phase 0 done: this plan, a prototype runtime and
two example dashboards in [`dashboards/`](dashboards/)). Amended 2026-11-21 (R15): actions on any
element and preset buttons, actions from pop-up windows, dashboards that open in a pop-up window;
then, after the merge into `main`, diagrams built on the CoolPlot Diagrams tab (`diagrams.json`,
server-drawn SVG) instead of Plotly.
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
| [§3](#3-concepts) | Concepts: diagram, saved plot, dashboard, runtime, host, snapshot |
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
| R9 | The important GUI buttons (Solve, Update guesses, …) can be placed on a dashboard; they do nothing when the page is not served by the GUI | §5.2 (`cs-button`), §5.3 (`data-cs-action`), §9.3, §9.4 |
| R10 | The plan is in `docs/` but not served online | header of this file, D1 |
| R11 | Dashboards are saved and loaded with the other files of the model (ZIP export/import) | §9.7 |
| R12 | A few library models to test the feature | §11 |
| R13 | Nice features that make CoolSolve attractive; more customisable and practical than the EES Diagram window | §14 |
| R14 | No questions: decisions taken while hesitating are documented | §15 |
| R15 | *(follow-up, 2026-11-21)* A dashboard acts on the current session through predefined actions and buttons written in its HTML (e.g. *Solve*), and can appear as a pop-up window; only when served by the GUI | §4.2, §5.3, §9.3, §9.4, D22–D24 |

---

## 2. What exists today

- **Diagrams tab** — `gui/src/components/diagram/` ([gui.md §6.8](../gui.md),
  [coolplot_integration.md](../coolplot_integration.md), done in four phases): property charts of
  a fluid (p-h, T-s, h-s, …, psychrometric) **drawn by the server** with CoolPlot's C++ core.
  `POST /api/v1/diagram/scene` takes a **diagram definition** (fluid, CoolPlot spec, cycles naming
  the model's variables: inferred states or arrays) and returns the SVG, the axes (limits in
  display units, log flags, `x_unit`/`y_unit` as SI × scale + offset), the states of the cycles
  and the labels. The SVG itself names its axis limits and log flags (`data-x-limits`,
  `data-x-log`, …) and its plot area (`rect.cp-background`), and tags its items (`data-item=
  "cycle/Cycle/states"`, `"iso/T/…"`): the tab maps the pointer onto it this way
  (`svgCoords.ts`). The definitions are a list, `diagrams.json` (`{"format":
  "coolsolve-diagrams", "diagrams": [{name, fluid, spec, cycles}]}`), stored in the session
  (`GET`/`PUT /api/v1/diagram/definitions`) and in the ZIP bundle; the tab edits the first one
  (*Multiple diagrams* is still to do, gui.md §14.4). Downloads SVG, PDF and the chart spec; the
  LaTeX report draws every stored definition as a vector figure; `coolsolve --diagram
  diagrams.json -o chart.svg model.eescode` draws one from the command line. The CoolPlot plan
  foresaw this one: "an HTML page of the folder can ask the server for the SVG of a named diagram
  of the solved model and show it among its overlaid variables" (coolplot_integration.md §3.4).
- **Parametric tab** — `ParametricStudy.tsx`: 1-D line plots and 2-D contours of a sweep. The
  sweep variables are the variables flagged `isImposed` by `POST /api/v1/parse` (lines
  `name = number`, scalars only, detected in the parse handler of `src/server.cpp`). Each grid point
  is solved on a copy of the code where the line is rewritten by the lambda `replaceImposedValue`
  of the parametric handler. Each point keeps **all** its variables (`ParametricGridPoint.variables`).
  The last study is saved in the bundle as `parametric_studies.json`.
- **Integral tab** — `IntegralTable.tsx`: the trajectory table of an `INTEGRAL` model and a line plot
  of one column. The trajectory is `GET /api/v1/integral/result` (columnar JSON) and, after a bundle
  round trip, the CSV `<model>-integral.csv`. The `.sol` file also carries it (`# IntegralTable`).
- **Exports of plots** — `gui/src/utils/exportPlots.ts` captures the visible Plotly charts
  (parametric, integral) as PNG for the LaTeX report; the diagrams come from the server as PDF.
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
| **Diagram** | A named definition of `diagrams.json` (§2): a CoolPlot chart of a fluid with the model's cycles, drawn by the server as SVG from the current results. |
| **Saved plot** | A named parametric or integral plot, stored as a **spec** (what to draw, not the picture) in `plots.json` and drawn with Plotly from the current results. |
| **Ref** | The name of a diagram or of a saved plot (one namespace): a bottom tab shows it, a dashboard uses it (`<cs-plot ref>`). |
| **Dashboard** | An HTML file of the model folder, `dashboards/<name>.html`, with its images in `dashboards/assets/`. Written against the contract of §5. |
| **Runtime** | `coolsolve-dashboard.js`: the script a dashboard includes. It finds the data, binds the elements, draws the plots, runs the player and talks to the host. |
| **Host** | The CoolSolve GUI page that shows a dashboard (a bottom-panel tab, or a pop-up window of its own). It feeds the data and passes the actions on. |
| **Editor window** | The GUI window that owns the editor of the session (the main window). It carries out the actions of every dashboard of the session, popped-out ones included (§9.4). |
| **Data** | One JSON document, `coolsolve-dashboard-data/1` (§6): results, inputs, integral table, parametric study, the diagrams (SVG) and saved plots a dashboard uses, solve status. |
| **Hosted mode** | Inside the GUI: live data, working inputs and buttons. |
| **Standalone mode** | Without the GUI: a **snapshot** of the data (inlined in an exported file, or `results.js` next to the dashboard). Inputs and buttons are shown but disabled, with a tooltip. Animation, plots, zoom and the layout helper work. |

---

## 4. User experience

### 4.1 Several diagrams and plots, as tabs

- The **Diagrams** tab (right panel) stays what it is: the chart being worked on.
- A **Save** button (pin icon) in its status line asks for a name and adds a copy of the current
  definition to `diagrams.json` under that name. This is the *Multiple diagrams* item of
  gui.md §14.4, done the way the request puts it ("keep the current window, save it, new tab").
  The same button on the plot of the Parametric tab and on the plot of the Integral tab stores
  that plot's settings as a saved plot (`plots.json`).
- Each saved diagram or plot becomes a **tab of the bottom panel**, after *Console · Parametric ·
  Lookup Tables · Integral*, with the plot icon and its name. A diagram tab is a `DiagramView`
  (which "draws any definition", gui.md §6.8: zoom, hover and badge included); a plot tab draws
  its Plotly figure. Both are redrawn after every solve (a cycle on a p-h chart follows the
  cycle). Tab menu: rename, duplicate, delete, *edit in the Diagrams tab* (diagrams: the tab then
  edits that definition, through the selector coolplot_integration.md §3.4 announced), *pop out*
  (its own window), *copy name* (to paste in a dashboard).
- The bottom tab bar scrolls horizontally when it is full; the panel can be maximised
  (double-click on its header) so that a tab fills the window.
- Both files are in the bundle (§9.7), and the LaTeX report already draws every definition of
  `diagrams.json` as a vector figure: saved diagrams are in the report for free.

### 4.2 Dashboards in the GUI

- A **Dashboards** drop-down in the toolbar: the dashboards of the model, *Add dashboard…* (§4.3),
  *Import HTML…*, and for the dashboard shown: *Export standalone file*, *Download source*,
  *Check*, *Delete*.
- A dashboard opens as a **tab of the bottom panel** (next to the saved plots), with three buttons
  in the tab header: **pop out** (a pop-up window of its own, without browser toolbars, which can go
  to a second screen), **presentation** (full screen, no CoolSolve chrome: lectures, demos) and
  **reload** (re-read the file).
- **Open in a window by default.** A dashboard can ask to open directly in its pop-up window:
  `<cs-dashboard open="window">` (§5.1), set by its author or by *Always open in a window* in the
  tab menu, which changes that attribute in the file (D24). Its tab then stays as a small
  placeholder with *Bring to front* and *Show here*. Windows open on a click only: loading a model
  adds the tabs of its dashboards and opens no window.
- **Live**: after each solve the dashboard receives the new results; while a solve runs it shows
  the status (`solving`, block *i*/*n*); a failed solve leaves the last results, marked stale.
  A popped-out window stays in sync with the main window (§9.4).
- **Actions**: the buttons and clickable elements of a dashboard (§5.2, §5.3) act on the session:
  *Solve* on a dashboard is the toolbar *Solve* (it becomes *Try Harder* after a failure), also from
  a pop-up window, where the main window carries it out (§9.4); the console says which dashboard
  asked. Inputs changed on a dashboard rewrite the line of the model in the editor, as an undoable
  edit (§7.3). All of this exists only in the GUI: a standalone dashboard shows the same controls,
  disabled (§4.4).

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
follows the GUI), `open` (`tab` by default, or `window`: the GUI opens the dashboard in a pop-up
window, §4.2) and `window-size` (`1100x750`; by default the size of the stage plus the margins).
`open` and `window-size` are read by the GUI (from the source, without running the page) and
ignored elsewhere.

### 5.2 Elements

| Element | Attributes | Renders |
|---|---|---|
| `<cs-stage>` | `width`, `height` (design units), `fit` (`width` default, `contain`) | A canvas in design coordinates, scaled as one piece (§5.6). Children with `x`, `y` (and `anchor`) are placed on it. |
| `<cs-value>` | `name` or `expr`, `unit` (display unit, §5.7), `scale`, `decimals`, `label`, `hide-unit`, `warn-above`, `warn-below`, `units` (unit of an `expr`) | Label, number, unit. Tooltip: the variable name and whether it is an input. Red when outside the warning limits. |
| `<cs-input>` | `name`, `label`, `min`, `max`, `step`, `unit`; `widget` (`slider` default, `number`, `select`, `radio`), `options` (string inputs) | An input of the model (§7): slider and number box, or a choice for a string variable. |
| `<cs-button>` | `action`: `solve`, `update-guesses`, `stop`, `reset-inputs`; text content = label | A GUI action (§9.3, §9.4), carried out by the editor window wherever the dashboard is. Disabled in standalone mode, with a tooltip, and while the action is not available (as the toolbar: *Stop* only while solving). `solve` becomes *Try Harder* as in the toolbar. The styled shorthand of `data-cs-action` (§5.3). |
| `<cs-status>` | — | `solved` / `failed` / `solving` / `stale` pill, message, time; "snapshot of …" in standalone mode. |
| `<cs-plot>` | `ref` (a diagram or a saved plot, §4.1) **or** `x`, `y` (columns of the integral table or of the parametric study, `y` may list several), `source` (`integral` default, `parametric`), `x-unit`, `y-unit`, `log-x`, `log-y`, `x-range`, `y-range`, `cursor`, `band-x` + `band-label`, `title`; for a diagram, `width` and `height` (the figure the server draws, in points, default 640 × 440; it then scales with the page) | A **diagram**: the CoolPlot SVG drawn by CoolSolve (the picture of the Diagrams tab, the report and `--diagram`), with the states' values as tooltips; `cursor="h_x, P_x"` places a marker that follows the player (§8.3). A saved plot or columns: an interactive chart (Plotly: zoom, hover); with `cursor`, a marker follows the player. |
| `<cs-player>` | `duration` (s for the whole range at ×1, default 10), `realtime`, `autoplay`, `no-loop`, `source` (`integral` default, `parametric`) | Animation controls (§8): start, step back, play/pause, step forward, scrubber, position, speed, loop. |
| `<cs-table>` | `names` (list), `units`, `decimals` | A small value table (name, value, unit), e.g. the state points. *(v1.1)* |
| `<cs-log>` | `lines` | The last lines of the GUI console (hosted mode). *(v1.1)* |

Elements marked v1.1 come after v1. The prototype has the v1 elements except `widget` and
`options` of `cs-input`, `source="parametric"`, `realtime`, saved plots of `plots.json`, and
`width`, `height` and `cursor` of a diagram (it shows the SVG of its data as is).

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

**Actions on any element.** `<cs-button>` is a ready-styled button; any other element, HTML or
SVG, can carry the same actions, so that a schematic is clickable like a real control panel and a
scenario is one line of HTML (D23):

| Attribute | Effect |
|---|---|
| `data-cs-action="solve"` | On click (or <kbd>Enter</kbd> / <kbd>Space</kbd>), carries out a GUI action: `solve`, `update-guesses`, `stop`, `reset-inputs` (§9.4). The runtime makes the element focusable (`role="button"`, `tabindex="0"`, pointer cursor) and disables it as a `cs-button` (standalone mode, action not available): class `cs-disabled`, tooltip |
| `data-cs-set="T_ev = -5; T_cond = 45"` | On click, sets inputs (§7.1) to constants — numbers, or `'R1234yf'` for a string input — as one undoable edit: a **preset**. With a `data-cs-action` on the same element, the action follows: `data-cs-action="solve"` sets and solves |

```html
<button data-cs-set="T_ev = -5; T_cond = 45" data-cs-action="solve">Cold day</button>
<button data-cs-set="R$ = 'R1234yf'" data-cs-action="solve">Switch to R1234yf</button>
<g class="compressor" data-cs-action="solve"> … </g>   <!-- click the compressor to solve -->
```

The runtime puts the class `cs-hosted` or `cs-standalone` on `<html>`, so that the author can hide
what only makes sense in the GUI: `:root:not(.cs-hosted) .gui-only { display: none }` (hidden
also during the 400 ms the runtime waits for a host).

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
`seek(t)`), `CoolSolve.setInput(name, value)` (or `setInputs({T_ev: 5, T_cond: 45})`, one edit)
and `CoolSolve.action(name)` (no-ops in standalone mode, with a console notice),
`CoolSolve.hosted` (true in the GUI) and `CoolSolve.onStatus(fn)` (the `status` of §9.3, for a
custom control that follows the solve). Nothing else of the page is reachable from the host
(§9.2).

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

A name that is not a variable, an input that is not an input, a diagram or saved plot that does not exist, an
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
  "plots": {
    "ph": { "kind": "diagram", "title": "R22 p-h diagram", "fluid": "R22", "diagram": "ph",
            "svg": "<svg … data-x-limits=\"90.7 577.6\" data-y-log=\"1\" …>…</svg>",
            "axes": { "x": { "key": "h", "units": "kJ/kg", "scale": 0.001, "offset": 0 },
                      "y": { "key": "p", "units": "bar", "scale": 1e-05, "offset": 0 } },
            "states": [ { "cycle": "Cycle", "number": 1,
                          "vars": { "p": "P[1]", "h": "h[1]", "T": "T[1]" } } ] },
    "cop_vs_tev": { "kind": "parametric", "title": "COP vs T_ev", "x": "T_ev", "y": ["COP"] }
  }
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
- `plots`: the diagrams and saved plots **the dashboard uses** (its `cs-plot ref`), with what they
  need to be drawn without CoolSolve (§9.10). A diagram carries the SVG of the server, drawn for the
  last solve at the size the dashboard asks; `axes` say how a model value maps onto it (display =
  value × `scale` + `offset`, from the scene's `x_unit` and the model's units, so °C works too);
  `states` are the cycles' states of the scene, with their variables. A saved plot carries its
  spec; its data are the parametric runs or the integral columns above.
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
  in the editor undoes a slider move. Several values sent at once (a preset button, §5.3) are one
  edit, undone by one Ctrl+Z. Without a main window (a popped-out dashboard whose main window was
  closed) the host rewrites the session's code on the server (§9.4).
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
loop, a vertical line on a time plot, a point on a CoolPlot diagram placed through its `axes`,
§9.10), and the `onFrame` functions (mechanisms). Steady values stay as solved. See
`dashboards/engine.html`: crank, rod, piston and intake valve follow the crank angle, the gas
colour follows its temperature, the spark and the flame show at the right angles, the intake flow
moves while the valve is open, the p–V loop and the pressure and temperature curves carry the
cursor.

### 8.4 Later

- **Parametric player** (`source="parametric"`): step through the runs of a 1-D sweep (each run
  keeps all its variables), so that the schematic and the cycle on the p-h diagram move with the
  swept evaporation temperature. Cheap, and very telling in lectures. On a diagram, the server's
  cycle is the one of the solved run: while the player runs, the runtime hides it (its items are
  tagged `data-item="cycle/…"`) and draws the current run's states over the fluid lines, joined by
  straight segments, through the `axes` mapping. Exact process paths per run would need a scene per
  run (`steps` points along each path): possible for a few runs, left for later.
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
gui/src/plots/figures.ts      Plotly figures of parametric and integral plots, shared by
                              ParametricStudy, IntegralTable, the saved-plot tabs and the runtime
gui/src/components/
  DashboardTab.tsx            host: sandboxed iframe + message bridge (§9.2, §9.3)
  DashboardWindow.tsx         the page of a popped-out dashboard: slim header + host (§9.4)
  DashboardMenu.tsx           toolbar drop-down, Add dashboard dialog (§4.3)
  SavedPlotTab.tsx            a saved diagram (diagram/DiagramView) or plot as a bottom tab (§4.1)
gui/src/utils/
  sessionActions.ts           solve, stop, update guesses: shared by toolbar, shortcuts, hosts
  sessionChannel.ts           the BroadcastChannel of a session, editor ownership (§9.4)
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
| dashboard | `{source: 'coolsolve-dashboard', version: 1, type: 'ready', plots: [{ref, width, height}]}` | Loaded; send the data, with these diagrams drawn at these sizes (§9.10) |
| host | `{source: 'coolsolve-host', type: 'data', data, theme}` | Data (§6) after load and after each solve |
| host | `{source: 'coolsolve-host', type: 'status', status}` | Solve progress, failure, *Try Harder* available |
| dashboard | `{… type: 'action', action}` | `solve`, `update-guesses`, `stop`, `reset-inputs` |
| dashboard | `{… type: 'set-input', values}` | One or several inputs changed, `{T_ev: 5}` (§7.3), applied as one undoable edit |

`status` carries `state`, `message`, `canTryHarder`, `actions` (which actions are available now,
as the toolbar: `stop` only while solving) and `editor` (`open`, or `closed` in a pop-up without
an editor window, §9.4). The host carries out the messages of one dashboard **in order**, so that
a `solve` sent after a `set-input` (a preset button, §5.3) solves the rewritten code.

Without an answer to `ready` within 400 ms the runtime takes the standalone path. Checked with the
prototype, the host's messages being posted to the page: `data` enables the inputs and buttons
and applies the theme; `cs-button`, `cs-input` (on release of the slider), the clickable
compressor (click and <kbd>Enter</kbd>) and a preset (one `set-input` with two values, then
`solve`) post the messages above; a `status` *solving* disables them all and nothing is posted;
*failed* with `canTryHarder` turns *Solve* into *Try Harder*. An unknown action, a `data-cs-set`
on an output, a non-constant or a string given to a number are listed in the problem badge and
their element is disabled.

### 9.4 Windows and actions

**Opening a window.** *Pop out* (tab header), or opening a dashboard whose file says
`open="window"` (§4.2), calls `window.open('/?dashboard=<name>', 'coolsolve-<session>-<name>',
'popup,width=…,height=…')`:

- a window without browser toolbars, sized by `window-size` or else by the stage of the dashboard;
  the window name makes a second click bring the same window to the front instead of opening
  another one;
- it uses the session of the main window: the cookie session as it is, or a tab session handed over
  through `localStorage` as the include links do (`coolsolve_include_opener`), then adopted
  (`adoptTabSession`); the parameter is consumed and removed from the address;
- it shows the dashboard alone, full size, under a slim header: title, status, *presentation*,
  *show in the main window* (closes the pop-up, the dashboard goes back to its tab);
- it is opened from a click only (menu, tab, button), never at load, since browsers block pop-ups
  without a user gesture. When the browser blocks it anyway, the dashboard opens in its tab with a
  notice. Safari ignores the `popup` feature and opens a tab: the dashboard works the same.

**Who carries out an action.** Every action and every input change is carried out by the
**editor window** (the main window, which owns the editor of the session), wherever the dashboard
is (D22):

- `solve` sends the editor's text and initials (`POST /api/v1/solve` takes them from the client:
  edits not yet solved exist only in the editor) and runs the progress stream; `update-guesses`
  refreshes the initials of its store; input changes and `reset-inputs` edit its editor, undoably
  (§7.3); `stop` cancels the solve it started.
- The toolbar, the keyboard shortcuts and the hosts call the same functions
  (`gui/src/utils/sessionActions.ts`, moved out of `Toolbar.tsx`), so *Solve* on a dashboard is the
  toolbar's *Solve*, *Try Harder* included.
- A dashboard in a bottom tab of the editor window calls them directly. A popped-out host forwards
  them on the channel below.
- The console of the editor window says which dashboard asked (`>>> Solve requested by dashboard
  main`): a solve started from a second screen is not a mystery.

**The channel.** The windows of one session keep in step through a
`BroadcastChannel('coolsolve-<session>')`. Only the hosts (GUI pages) are on it; a dashboard, with
its opaque origin, cannot join it (§9.2).

| Message | From | Meaning |
|---|---|---|
| `who-owns` | pop-up | At load: which window owns the editor? |
| `owner` | editor window | I own the editor (at load, on focus, in answer to `who-owns`) |
| `owner-gone` | editor window | Sent on `pagehide` (closed or reloaded) |
| `action {id, dashboard, action}` | pop-up | Carry out this action |
| `set-input {id, dashboard, values}` | pop-up | Apply these input changes (§7.3) |
| `ack {id}` | editor window | Received; the outcome follows as `status` and `solved` |
| `status {status}` | editor window | Solve progress, relayed to the dashboards of every window |
| `solved` | editor window | Re-read the data (`GET /api/v1/dashboard-data`) |
| `code-changed` | editor window | The code changed (the list of inputs may have changed) |

If two windows with an editor share a session (a cookie session open in two tabs), the one focused
last claims ownership with `owner` and the other yields. Same browser only; a second device would
need server events (`GET /api/v1/session/events`, SSE), listed in §14.

**Without an editor window.** When a pop-up has no owner (no answer to `who-owns`, an
`owner-gone`, or no `ack` within 500 ms), it carries out the actions itself on the session's code as
the server last received it (last solve or save): it solves that code, `stop` cancels its own solve,
`update-guesses` updates the session's initials, input changes rewrite the session's code on the
server (§7.3). The status of its dashboards says *editor closed: solving the code of the last
solve*, and the header offers *Open the editor*, which opens `/?sid=<session>` in a new window that
shows the session's model and becomes the owner. A second-screen dashboard thus stays usable when
the laptop window is closed, and says what it solves.

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
| `GET` / `PUT` / `DELETE` | `/api/v1/plots`, `/api/v1/plots/{id}` | Saved plots (§4.1); diagrams keep the existing `GET`/`PUT /api/v1/diagram/definitions` and `POST /api/v1/diagram/scene` |
| `GET` | `/d/{session}/{name}.html`, `/d/{session}/assets/{file}`, `/d/{session}/coolsolve-dashboard.js`, `/d/{session}/plotly.js` | The dashboard document and its files, with the policy of §9.2 |

The session ID is in the path so that the relative URLs of a dashboard (`assets/x.svg`,
`coolsolve-dashboard.js`) resolve as in the bundle folder; it is the same capability as `?sid=`.
New session state: `dashboards` (name → HTML), `dashboardAssets` (file → bytes), `savedPlots`
(JSON). *New*, *Open* and the library loader reset them as they reset the lookup tables;
`showSessionModel()` (`gui/src/utils/sessionModel.ts`) reads them back.

### 9.6 Names

Dashboard names: `[A-Za-z0-9_-]+`, as table names. Refs (§3): the `name` of a diagram of
`diagrams.json` (free text today, "unique in the document") or the id of a saved plot (generated
from its title); the GUI keeps them unique across both files and proposes `[A-Za-z0-9_-]` names
when saving, so that a ref is easy to type in HTML.

### 9.7 Bundle and model folder

```
my_model/
├── my_model.eescode, .initials, .sol, coolsolve.conf, my_model-<table>.csv, …
├── diagrams.json                 diagrams (exists today; several of them with §4.1)
├── plots.json                    saved parametric and integral plots (§4.1)
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
  `diagrams.json` travels already.
- A model opened from disk (`POST /api/v1/files/open`) reads `diagrams.json`, `plots.json` and
  `dashboards/` of its folder, like its lookup tables.
- Library models may carry the same layout (§9.11).

### 9.8 Standalone export

One HTML file: the dashboard, with its `<script src="coolsolve-dashboard.js">` replaced by the
inlined runtime, Plotly inlined (the cartesian bundle the GUI already uses, about 1.4 MB, only when
the dashboard has a saved plot or a plot of columns: **diagrams need no Plotly**, they are SVG in
the data), the data inlined as `<script type="application/json" id="coolsolve-data">`, and the
images as data URIs. A diagram is drawn by the same engine as `coolsolve --diagram`, at the
dashboard's size, in its theme (`theme="auto"` embeds the light and the dark variant and shows one
by `prefers-color-scheme`; about 80 kB each for the p-h chart of `CSL-0025`). Deterministic output
(same model, same file), so that exports can be committed or compared. Written by the C++ exporter
(server and CLI); the runtime and Plotly come from the embedded GUI assets, so the export needs a
GUI build (D5); a CLI-only build says so.

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
  `--from-sol`) and writes the standalone file (§9.8); the diagrams come from the `diagrams.json`
  next to the model (or `--diagram <file>`, the option `--diagram` already has), drawn as
  `--diagram` draws them.

### 9.10 Plots

**Diagrams: CoolPlot SVG drawn by the server** (D6).

- A diagram of a dashboard is the definition of `diagrams.json` its `ref` names, drawn by the
  engine of the Diagrams tab (`src/diagram/`, `POST /api/v1/diagram/scene`): the same picture as the
  tab, the report and `coolsolve --diagram` — fluid lines, cycles of the model, labels laid out by
  CoolPlot, theme, badge. No CoolProp in the browser.
- **Only what the dashboard uses, at its size.** The runtime lists its diagrams and their sizes in
  `ready` (`cs-plot ref`, `width`, `height` in points, default 640 × 440); the host asks one scene
  per diagram after each solve and when the theme changes, and sends the SVGs in `data`. The
  checker and the exporter read the same list from the HTML. CSL-0025's p-h chart takes about a
  second with its solve from the command line (`--diagram`, 1.3 s in all).
- **Shown as is, scaled with the page** (`width: 100%`, aspect kept). The runtime prefixes the ids
  of each SVG (`coolplot-plot-area` clip path) so that several charts can share a page.
- **Values onto the chart.** The SVG names its axis limits and log flags (`data-x-limits`,
  `data-y-log`, …) and its plot area (`rect.cp-background`); with the `axes` of the data
  (display = value × scale + offset) the runtime maps any model value onto the chart, as
  `svgCoords.ts` maps the pointer the other way. Used for the tooltips of the states (v1), the
  `cursor` of a player (§8.3) and the cycle of the parametric player (§8.4). Checked with the
  prototype: the four state positions it computes for `CSL-0025` fall on CoolPlot's own markers
  to 0.01 viewBox unit.
- Not in v1: zoom (hosted mode could ask a scene with a `view_SI`, as `DiagramView` does) and the
  state under the pointer (`POST /api/v1/diagram/state` through the host).

**Saved plots: Plotly** (parametric and integral plots).

- A saved plot is a spec: `{id, title, kind: 'parametric' | 'integral', …}` with the settings of
  the tab it was saved from (sweep axes, contour or lines; columns).
- Its Plotly figure is built by one module, `gui/src/plots/figures.ts`, used by the Parametric and
  Integral tabs, the saved-plot tabs and the runtime: the three draw the same picture. The plots of
  integral columns written in a dashboard (`cs-plot x y`) use it too.

### 9.11 Library

- Allow-list: `dashboards/*.html` and `dashboards/assets/*.{svg,png,jpg,jpeg,webp}` of a model
  folder, `diagrams.json` (foreseen by coolplot_integration.md §3.4: "opening a model shows the
  chart its author made") and `plots.json`, in `cmake/embed_library.cmake` **and** `EMBED_FILES`
  of the library's `tools/build_index.py` (they mirror each other); `LibraryStore::loadModel`
  loads them into the session; `model.json` lists them (`"dashboards": ["main"]`).
- The library dashboard (`library.html`) shows a dashboard badge, and *View dashboard* on the model
  page: `GET /library/<model path>/dashboards/<name>.html` served as a standalone snapshot built
  from the embedded `.sol` (values and `# IntegralTable`), with no solve. Its diagrams are drawn
  from the `.sol` values (the diagram engine takes a map of values, as its tests do) and cached by
  model, diagram, size and theme: the fluid lines cost CPU, the public server draws each once. The
  public site gains a gallery of animated models that anyone can explore.
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
   integration range, the diagrams (name, fluid, chart, the variables of their states) and saved
   plots (id, title, kind), the parametric study if any, and — when
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
| error | a `cs-plot ref` that is neither a diagram nor a saved plot (with the names there are); `x`/`y` that are not integral columns (with the `$IntegralTable` line to add) |
| error | an action of `cs-button` / `data-cs-action` that is not one of §9.4; a `data-cs-set` on a variable that is not an input, or whose value is not a constant (a number, or a quoted string for a string input) |
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
| `CSL-0025` *heat_pump_cycle_r22_30kw* (55 eq.) | steady, R22 | values on a schematic, numeric inputs, *Solve* / *Update guesses*, its p-h diagram (CoolPlot, the cycle from the arrays `P[i]`, `h[i]`) in `diagrams.json`, pipes coloured by temperature, flows moving with the mass flow, a warning limit | `dashboards/heat_pump.html` (prototype); to become the library's `dashboards/main.html` |
| `CSL-0345` *r134a_isenthalpic_throttle* (16 eq., variant R1234yf) | steady | the input-rule extension: **string input** (radio R134a / R1234yf) and **array-element inputs** (`T[1]`, `P[1]`); parity with the radio button and edit boxes of the EES Diagram window | throttle sketch, inlet/outlet states, T-h diagram |
| `CSL-0339` *gas_spring_oscillator* (28 eq.) | dynamic (time, 5 s) | the EES animation of the original (`Ball.top`): `data-cs-translate`, colour by gas temperature, `realtime` player, L(t) and P(t) with cursors. The file has no `$IntegralTable`: the checker must tell which line to add for P and T | ball in a tube over a bottle |
| `CSL-0024` *single_cylinder_engine_weibe* (72 eq.) | dynamic (crank angle) | a mechanism drawn by `onFrame`, declarative show / rotate / fill / flow, p–V loop on log axes with a cursor, a band for the combustion window | `dashboards/engine.html` (prototype), on a variant integrating from −360° to +180° (below) |
| `CSL-0009` *dhw_tank_dynamic* (20 eq.) | dynamic, 18 000 rows | a large integral table (decimation of exports, interpolation speed), schedules shown with `data-cs-show` (heater on, tap draw), tank colour by temperature, 5 h played in 15 s | tank, heater, tap, T(t), heat flows |
| `CSL-0019` *orc_simple_r245fa* (182 eq.) | steady, large | about forty labels on a P&ID-like schematic, speed of the data builder and of the page, a T-s diagram | stretch goal |

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
| **2 — GUI host and bundle** | Session state, routes (§9.5), bundle and disk folder (§9.7), toolbar menu, bottom tabs, sandboxed host and message bridge, `sessionActions.ts` shared with the toolbar, buttons and `data-cs-action`, live updates, pop-up windows (`open="window"`, *Always open in a window*), `BroadcastChannel`, editor ownership and the pop-up fallback (§9.4), presentation mode | A dashboard of a bundle opens in the GUI, solves from its button — also from its pop-up window, and still after the main window was closed — follows a solve made in the editor, and survives a download/upload round trip |
| **3 — Inputs** | `src/model_inputs.cpp` (detection with array elements and strings, rewrite), parametric study on the shared detection, `/api/v1/inputs`, undoable editor edits, `autosolve`, *Reset inputs*, string widgets, `data-cs-set` presets | `CSL-0345`: the radio button switches the refrigerant, the editor shows it, Ctrl+Z undoes it; a preset button of `CSL-0025` sets two inputs and solves once |
| **4 — Diagrams and saved plots** | *Save* on the Diagrams tab (several definitions in `diagrams.json`, the selector of the tab: *Multiple diagrams* of gui.md §14.4) and on the Parametric and Integral plots (`plots.json`); bottom tabs (`DiagramView`, `gui/src/plots/figures.ts`); scenes in the data at the size and theme a dashboard asks (§9.10), ids prefixed, state tooltips; diagrams in the standalone export and in `--dashboard-export`; `cs-plot ref` | The p-h diagram of `CSL-0025` is a bottom tab and is in its dashboard, live and in the standalone file |
| **5 — Animation** | Player complete (`realtime`, events, `CoolSolve.player`), `data-cs-translate` / `data-cs-scale-y`, cursors on all plot kinds (on diagrams through the `axes` mapping), parametric player (with the run's cycle drawn over a diagram, §8.4) | `CSL-0339`, `CSL-0024` (variant) and `CSL-0009` animate smoothly |
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
| 3.5 GUI frontend | needed | `types.ts` (dashboards, saved plots, data), `client.ts` (routes through `request()`), `modelStore.ts` / `uiStore.ts` (dashboards, saved plots, bottom tabs), `sessionModel.ts` (read them back), `App.tsx` (bottom tabs, the `?dashboard=` window), `Toolbar.tsx` (menu; its solve, stop and update-guesses handlers move to `sessionActions.ts`), `diagram/DiagramTab.tsx` (*Save*, selector of the definitions) / `ParametricStudy.tsx` / `IntegralTable.tsx` (*Save*), `diagram/DiagramView.tsx` reused by the bottom tabs, new components of §9.1, `CodeEditor.tsx` (undoable input edits) |
| 3.6 Server / REST | needed | Routes of §9.5 through `getSession()`; session IDs of `/d/{session}/…` checked with `isValidSessionId()`; `tests/test_dashboard_api.py` (`--auto`) |
| 3.7 ZIP bundle | needed | §9.7 in `createZipBundle` / `extractZipBundle`; round-trip test; layout documented in `contributing.md` §3.7 and `gui.md` §7.2 |
| 3.8 Debug mode | not applicable | — |
| 3.9 Console | needed | Checker problems and input rewrites reported in the console with the usual markers (`[Warning]`, `ERROR`) |
| 3.10 LaTeX report | mostly done | Every definition of `diagrams.json` is already a vector figure of the report (saved diagrams included); saved plots are Plotly charts, which `exportPlots.ts` captures when visible |
| 3.11 Examples | needed | `examples/dashboards/` with a small model and a dashboard for the CLI tests (`--dashboard-check`, `--dashboard-export`) |
| 3.12 README and docs | needed | §17 |
| 3.13 Build system | needed | Second Vite entry for the runtime; the runtime and `docs/dashboards.md` embedded (the latter in every build, for the prompt); library allow-list on both sides |
| Tests | needed | Catch2: data builder (live and `.sol`), input detection and rewrite (numbers, array elements, strings, comments and units kept), checker (each rule, negative cases), exporter (deterministic, decimation). Headless Chromium (pattern of `tests/test_gui_shortcuts.py`): the examples render without problem badge; a label keeps its place on the schematic from 1280 to 375 px; the player advances and a cursor moves; a hosted dashboard's button solves and an input rewrites the editor; a `data-cs-action` element of an SVG solves; a preset sets two inputs as one undoable edit, then solves; a button of a popped-out dashboard solves through the main window, then on its own once the main window is closed; the sandbox cannot `fetch('/api/v1/files/eescode')` |

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
   sliders, the outputs as a table, the diagrams and saved plots and, for an `INTEGRAL` model, a player and one
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
| D2 | "Several diagrams" = **several named definitions of `diagrams.json`**: the Diagrams tab stays the chart being worked on; *Save* adds a copy under a name, shown as a bottom tab. The same gesture saves the parametric and integral plots (`plots.json`) | Several Diagrams tabs on the right; saved plots only inside dashboards; one new file for every kind of plot | The request says "keep the current window, save it, new tab at the bottom"; `diagrams.json` is already a list made for it (coolplot_integration.md §3.4), read by the report and `--diagram`; dashboards use both by name |
| D3 | A dashboard is **a plain HTML file** with custom elements, `data-cs-*` attributes and optional JavaScript | A JSON layout drawn by the GUI; a GUI editor (EES-like) | LLMs write HTML/SVG/JS very well; users edit it in any editor; anything the web can draw is possible; a standalone file is natural. The request excludes editing in the GUI |
| D4 | Labels are placed on a **stage in design coordinates scaled as one piece** | Percent positions; SVG `foreignObject`; CSS container units | Exact placement at any size with text that scales with the drawing; works the same for inline SVG and images; `foreignObject` has browser quirks for inputs; percent positions keep the font size and drift with the aspect ratio. Checked at 1280 and 375 px |
| D5 | The runtime is **TypeScript in `gui/src/dashboard/`, built by Vite**; the standalone export needs a GUI build | Hand-written JavaScript embedded in every build | Shares the figure code with the GUI (no duplication, §1 of the contributing guide); the distributed binaries (installer, server) are GUI builds. The prototype is plain JS only to run without a build |
| D6 | **Diagrams are the server's CoolPlot SVG**, sent in the data at the size the dashboard asks; parametric and integral plots use **Plotly**, as their tabs | Plotly for every plot (the first version of this plan, written before the Diagrams tab moved to CoolPlot); the scene's curves sent as data and drawn in the browser | One picture of a diagram everywhere (tab, report, `--diagram`, dashboards) with CoolPlot's labels, themes and badge, no CoolProp in the browser; the SVG names its axes and plot area, so values can be placed on it; tens of kB instead of 1.4 MB of Plotly in an export with diagrams only. Plotly stays where the GUI uses it, for zoom, hover and contours |
| D7 | The data document is built by **one C++ function** for the server and the CLI | Assembled in the GUI from its store | One builder for the GUI, the CLI export and the library snapshots (from `.sol`) |
| D8 | Dashboards run in a **sandboxed iframe** and speak `postMessage` only — library dashboards too | Same-origin iframe calling the REST API | Dashboards are code from LLMs and from other people's bundles, on a public server too; an opaque origin and `connect-src 'none'` make a malicious page harmless |
| D9 | An input change **rewrites the line of the model**, as an undoable editor edit | An override layer next to the code (EES keeps Diagram-window values apart from the equations) | The code stays the single source of truth for the editor, the bundle and the CLI; reuses the parametric substitution; no new state to save |
| D10 | Inputs include **constant-index array elements and string constants** | Scalars and numbers only, as the parametric table today | Needed by `CSL-0345` (`T[1]`, `P[1]`, `R$`) and by EES models whose Diagram window set strings for `$IF` directives; the parametric table gets array elements from the same detection |
| D11 | **Manual solve** by default; `autosolve` per dashboard | Always solve on change | Many library models take seconds; a slider that solves at each move would freeze the page |
| D12 | Animation is a **playback of the integral table** in the browser, interpolated; *follow the solve* later | Re-solving per frame; streaming during the solve | "Live, but using the results of the integral table" (the request); smooth, instant, works in standalone files |
| D13 | The bundle carries the runtime and `results.js` (generated, ignored at import), not Plotly | Nothing generated in the bundle; Plotly in every bundle | The unzipped folder opens offline from `file://` for a few kB; 1.4 MB in every bundle is too much, and the standalone export covers offline plots |
| D14 | Dashboards and saved plots are **bottom-panel tabs**, with pop-out and presentation; a dashboard may ask to open in its pop-up window (D24) | A right-panel tab; a separate page only | The request places saved diagrams at the bottom; one place for both; the panel is resizable and maximisable, and pop-out gives a full window |
| D15 | Buttons: `solve` (with *Try Harder*), `update-guesses`, `stop`, `reset-inputs` | Also *Save*, *Open*, *Debug solve* | The actions that make sense while exploring; file actions stay in the toolbar (a dashboard must not trigger downloads) |
| D16 | The GUI **does not call an LLM**: it generates the prompt, the user pastes it; the contract text has one source, `docs/dashboards.md` (markers), embedded in every build | Built-in LLM calls with an API key; a separate authoring guide | No key, no cost, offline, provider-neutral; one text for users and LLMs cannot drift |
| D17 | CLI options `--dashboard-prompt`, `--dashboard-check`, `--dashboard-export` | Sub-commands (`coolsolve dashboard check`) | Consistent with the existing long options of `main.cpp` |
| D18 | Display units by **scale factors only** in v1 (and `scale`) | Offsets (°C ↔ K, °F) | CoolSolve's units are fixed SI with °C; factors cover pressures, energies, powers, flows; offsets wait for a real need |
| D19 | Names are **case-insensitive**, trigonometric functions in **degrees** | JavaScript conventions | EES/CoolSolve semantics: `T_ev` and `t_ev` are one variable; `sin(30) = 0.5` in the model and on the dashboard |
| D20 | The test models are `CSL-0025`, `CSL-0345`, `CSL-0339`, `CSL-0024` (+ variant), `CSL-0009`, and `CSL-0019` as a stretch | Others (building `CSL-0042`, absorption chillers) | Small to large, steady and dynamic, the three Diagram-window cases of the register, a mechanism, a long table |
| D21 | The prototype's examples use **real CoolSolve results** (library `.sol`, and a run of the engine variant), copied into `*.results.js` | Synthetic data | Shows what a real dashboard looks like and already revealed a modelling question (§11) |
| D22 | Actions are carried out by the **editor window**; a pop-up forwards them on the `BroadcastChannel`; without an editor window, the pop-up acts on the session's code as the server last received it, and says so | Each window solves on its own from the server's copy; actions disabled in a pop-up once the main window is closed | Edits not yet solved exist only in the editor, and so do the guesses shown and the undo history; one place carrying out actions means two windows never solve two versions of the model at once; the fallback keeps a second-screen dashboard usable when the laptop window is closed |
| D23 | Actions on **any element** (`data-cs-action`) and **preset buttons** (`data-cs-set`), `cs-button` being the styled shorthand | `cs-button` only; actions from JavaScript only (`CoolSolve.action`) | A schematic becomes clickable like a real control panel; a scenario ("cold day", another refrigerant) is one line of HTML for an LLM, without JavaScript; the checker validates attributes statically, which it cannot do for a script |
| D24 | Opening in a window is a **setting of the dashboard file** (`<cs-dashboard open="window">`), toggled from the tab menu; windows open on a click only | A preference in the browser (`localStorage`); a manifest `dashboards/dashboards.json`; a floating panel inside the GUI page; windows opened at load | The setting travels with the dashboard (bundle, library) and the author or an LLM can write it; the GUI changes one attribute of a tag, as it changes one line for an input; a floating panel cannot go to a second screen; browsers block pop-ups opened without a click |

---

## 16. Risks and open questions

- **CSP and opaque origins.** `'self'` in a policy of a sandboxed document must still match the
  server's URLs for the runtime and assets (it is the URL's origin in CSP 3); to verify on Chromium,
  Firefox and Safari in Phase 2. Fallback: inline the runtime into the served document.
- **CoolPlot SVG in a page.** Labels are laid out by the server with its font metrics (`Inter,
  'Helvetica Neue', Arial, 'DejaVu Sans'`); a viewer without these fonts may see labels shift a
  little. Several charts on one page share ids unless prefixed (done by the prototype). The cost
  of the fluid lines is paid per solve and per diagram used: cache scenes until the next solve,
  and per model on the public server (§9.11).
- **Plotly in a sandbox.** Plotly needs no network or storage; its hover and zoom work in sandboxed
  iframes (to confirm in Phase 1).
- **Editor ownership** (§7.3, §9.4). Several windows of one session could edit at once. Designed:
  one owner announced on the `BroadcastChannel`, the last focused editor window, and a pop-up that
  acts on the server's copy when no owner answers. To verify with two tabs of one cookie session
  and with a reload of the main window during a solve requested from a pop-up.
- **Pop-up windows.** Blocked without a user gesture, so never opened at load; `popup` window
  features differ between browsers (Safari opens a tab). The tab stays the fallback.
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

- `README.md`: *Features* (dashboards, saved diagrams and plots), *Command Line Options* (three options),
  *File Formats* (`dashboards/`, `plots.json`), *Project Structure* (new files), *Documentation*
  table (`docs/dashboards.md`).
- New **`docs/dashboards.md`**: user guide and contract reference (the prompt is taken from it,
  §10.2), with screenshots; in the sidebar of `docs/docs.html`.
- `docs/gui.md`: routes (§4), session state (§5), layout and tabs (§6), bundle (§7.2), library
  (§7.3), §6.8 (*Save*, the selector of the definitions); remove *Multiple diagrams* from §14.4.
- `docs/contributing.md`: §3.7 (bundle layout), §3.5 (new components).
- `docs/model_library_support.md`: a feature row `CS-FEAT-DASHBOARD` in §3; the resolution path of
  `CS-GAP-DIAGRAM-WINDOW` and of case (a) of `CS-GAP-IF-NO-CONSTANT` (string inputs on a dashboard).
- `docs/coolplot_integration.md`: §3.4 (*Later, in model folders*) points to this plan once the
  library carries `diagrams.json`.
- Library: `docs/model_workflow.md` (dashboards of a model, checker in the tests), `EMBED_FILES`.
- `docs/versions.md` at the release.

---

## 18. The example files

In [`dashboards/`](dashboards/), next to this plan (not served either):

| File | What |
|---|---|
| `coolsolve-dashboard.js` | Prototype runtime (plain JavaScript, no dependency): stage scaling, `cs-value`, `cs-input`, `cs-button`, `cs-status`, `cs-plot` (diagrams: the CoolPlot SVG with state tooltips placed through its axes; columns: small SVG charts instead of Plotly), `cs-player`, the `data-cs-*` bindings, `data-cs-action` and `data-cs-set`, `CoolSolve.action` / `setInput(s)` / `onStatus` / `hosted`, the `cs-hosted` / `cs-standalone` classes, the expression language, units, problem badge, layout helper (<kbd>L</kbd>), and the host protocol of §9.3 (`data`, `status` with available actions and *Try Harder*) |
| `engine.html` + `engine.results.js` | `CSL-0024` animated: slider-crank mechanism, valves, spark, flame, gas coloured by temperature, p–V loop and p(θ), T(θ) with cursors, player (the "animated internal combustion engine" of the request) |
| `heat_pump.html` + `heat_pump.results.js` | `CSL-0025`: values on the schematic at the four states, COP, pipes coloured by temperature with moving flows, inputs and buttons, two preset buttons (`data-cs-set` + `data-cs-action="solve"`), a compressor that solves when clicked with a *click to solve* hint shown in the GUI only (all disabled standalone), saved P-h diagram |

**Opening them.** Double-click the `.html` file (it reads its `*.results.js` from the same folder,
`file://` is fine), or serve the folder: `python3 -m http.server -d docs/plans/dashboards`. Keys:
space (play/pause), ←/→ (step), <kbd>L</kbd> (layout helper; drag a label, its new tag is copied).

**Where the data come from.** `heat_pump.results.js`: the `.sol` baseline of `CSL-0025`, and its
p-h diagram drawn by CoolSolve 0.5.0 (`coolsolve --no-sol --diagram diagrams.json -o ph.svg
heat_pump_cycle_r22_30kw.eescode`) from this definition:

```json
{"format": "coolsolve-diagrams", "version": 1, "diagrams": [
 {"name": "ph", "fluid": "R22",
  "spec": {"diagram": "ph", "units": "engineering", "limits": "ACHP",
           "isolines": [{"key": "T", "num": 8}, {"key": "s", "num": 8}, {"key": "Q", "num": 5}],
           "settings": {"labels": {"mode": "rails", "symbols": "unicode"}}},
  "cycles": [{"name": "Cycle", "numbers": true, "arrays": {"p": "P", "h": "h"}, "closed": true}]}]}
```

(`"units": "engineering"` draws kJ/kg and bar; the Diagrams tab stores the model's units, which
the server, the report and dashboards draw as stored.) `engine.results.js`: a CLI run of `CSL-0024`
with the integration window extended to +180° (§11). Inputs were detected with the rule of §7.1,
units taken from the comments of the model. They follow the data contract of §6.

**Sketch of a generated prompt** (`CSL-0025`, abbreviated):

````markdown
You are writing a CoolSolve dashboard: one HTML file that shows the results of the model below.
Request: the cycle on a schematic with T, P and h at each state point, the COP in big, sliders for
T_ev, T_cond, DELTA_T_sc and epsilon_s_cp, a Solve button, and the p-h diagram "ph".

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
Diagrams: ph (R22, p-h, engineering units; cycle from the arrays P[i], h[i], states 1-4)
Integral table: none.

## Check your file
coolsolve --dashboard-check dashboards/main.html heat_pump_cycle_r22_30kw.eescode
````
