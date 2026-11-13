# CoolSolve GUI — Implementation Reference

## 1. Overview

CoolSolve ships with a **single-page web application (SPA)** served by an
embedded HTTP server (cpp-httplib) inside the CoolSolve binary.  Running
`coolsolve --gui` starts the server on port 8550 (by default) and opens the
user's default browser — no external runtime dependencies are required.

The GUI works identically in two deployment scenarios:

| Scenario | Backend | Frontend |
|----------|---------|----------|
| **Local** | CoolSolve binary starts an embedded HTTP server on `localhost` and opens the user's default browser | Same SPA served by the embedded server |
| **Online** | CoolSolve binary behind a reverse-proxy (nginx / Caddy) with session-based sandboxing | Same SPA served by the proxy or CDN |

**Design principle**: the GUI is a pure frontend layer.  All core functionality
stays in the existing libraries (`parser.cpp`, `solver.cpp`, `evaluator.cpp`,
`runner.h`).  The HTTP layer is a thin adapter that calls `CoolSolveRunner`
exactly as `main.cpp` does in CLI mode.

---

## 2. Technology Stack

### 2.1 Backend (C++)

| Component | Library | Notes |
|-----------|---------|-------|
| HTTP server | [cpp-httplib](https://github.com/yhirose/cpp-httplib) (header-only) | Zero-dependency, MIT-licensed, supports SSE and multipart uploads |
| JSON serialization | nlohmann/json (existing dependency) | REST request/response bodies |
| Real-time progress | Server-Sent Events (SSE) via cpp-httplib chunked content provider | Streams per-block solve progress without adding a WebSocket library |
| Static asset embedding | `cmake/embed_assets.cmake` custom build step | Embeds the compiled SPA into the binary as C++ byte arrays |

### 2.2 Frontend

| Component | Library | Version |
|-----------|---------|---------|
| Framework | React + TypeScript | 18.3 / 5.6 |
| Build tool | Vite | 5.4 |
| Code editor | Monaco Editor (`@monaco-editor/react`) | 4.7 |
| Plotting | Plotly.js (`react-plotly.js` + `plotly.js-basic-dist-min`) | 3.3 |
| State management | Zustand | 5.0 |
| Icons | Lucide React | 0.564 |

AG Grid (planned for the variable table) and Radix UI / Tailwind CSS (planned
for UI chrome) were considered but not yet adopted — the current UI uses plain
React and custom CSS with light/dark theme variables.

### 2.3 Rationale for Embedded-Server Approach

An Electron or Tauri-based desktop app was evaluated and rejected:

- Electron ships a full Chromium (~150 MB) and Node.js runtime.  CoolSolve
  already carries CoolProp (~30 MB); adding Electron triples the binary size.
- Tauri requires a Rust toolchain and introduces a second systems language.
- The embedded-server approach delivers the same UX with **zero additional
  runtime dependencies**, and the same frontend works online without
  modification.

---

## 3. Core C++ Changes

### 3.1 New Files

| File | Purpose |
|------|---------|
| `src/server.cpp` | HTTP server — all route handlers, session management, asset serving |
| `include/coolsolve/server.h` | Public API: `ServerOptions` struct + `startServer()` declaration |
| `cmake/embed_assets.cmake` | CMake script: reads `gui/dist/*`, generates `build/embedded_assets.cpp` and `build/embedded_assets.h` |
| `cmake/embed_library.cmake` | CMake script run at every build: embeds the model library snapshot (`build/embedded_library.cpp`, §9.6) |
| `cmake/embed_common.cmake` | Helpers shared by the two embedding scripts (MIME types, file → C array) |
| `src/library_store.cpp` / `include/coolsolve/library_store.h` | `LibraryStore`: index of the library snapshot, menu of verified models, loading of a model folder (no HTTP, no embedding: testable) |

### 3.2 Changes to Existing Files

| File | Change |
|------|--------|
| `main.cpp` | Added `--gui [port]` and `--no-browser` flags; calls `coolsolve::startServer()` |
| `CMakeLists.txt` | Added `COOLSOLVE_BUILD_GUI` option (default `ON`); FetchContent for cpp-httplib; asset embedding build step |
| `solver.h` / `solver.cpp` | Added `ProgressCallback` to `SolverOptions` for real-time SSE progress |

### 3.3 No Overhead on CLI Path

- `server.cpp` is only linked when `COOLSOLVE_BUILD_GUI=ON`.
- cpp-httplib is guarded behind `#ifdef COOLSOLVE_GUI`.
- Frontend assets are embedded as `const unsigned char[]` arrays; never
  touched unless the server starts.
- When `progressCallback` is `nullptr` (CLI mode), the compiler eliminates
  the null check — zero overhead on every Newton iteration.

---

## 4. REST API Reference

All endpoints are under `/api/v1/`.  Request and response bodies are JSON
unless noted otherwise.

### 4.1 Health & Lifecycle

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/health` | Health check; returns `{"status":"ok","coolpropReady":true/false}` |
| `POST` | `/api/v1/warmup` | Trigger CoolProp warmup (called automatically on startup) |
| `POST` | `/api/v1/keepalive` | Marks the caller's session as in use (204); every open tab calls it every 15 min, so that its session is not deleted as idle (§5) |

On server start, CoolProp warmup runs in a background thread.  The frontend
polls `/health` and shows a loading indicator until `coolpropReady` is `true`
(typically 1–3 s on first launch).

### 4.2 Session & Model State

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/api/v1/new` | Create a new empty session (clears code, initials, sol, debug output) |
| `POST` | `/api/v1/back` | Restore the previous session snapshot (one-level undo) |
| `GET` | `/api/v1/model-name` | Get the current model name |
| `PUT` | `/api/v1/model-name` | Set the model name |

Sessions are managed via a `coolsolve_session` cookie (32-hex-char ID).
Each session gets isolated in-memory state.  A `SessionSnapshot` is saved
before each destructive operation to support the "Back" button.

A browser tab can also be an **independent instance**: it sends its own
session ID in an `X-CoolSolve-Session` header (or a `?sid=` query parameter
for `EventSource` and downloads, which cannot set headers) and the server
uses it instead of the cookie. This is how the *Open in CoolSolve* links of
the model library open a model without touching the model of another tab
(§5, §7.3).

### 4.3 File Operations

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/files/eescode` | Get current source code |
| `PUT` | `/api/v1/files/eescode` | Set source code; body: `{"content":"..."}` |
| `GET` | `/api/v1/files/initials` | Get initial guesses text |
| `PUT` | `/api/v1/files/initials` | Set initial guesses text |
| `GET` | `/api/v1/files/sol` | Get last solution file content |
| `GET` | `/api/v1/files/conf` | Get solver configuration (as text) |
| `PUT` | `/api/v1/files/conf` | Set solver configuration |
| `POST` | `/api/v1/files/open` | Open the `.eescode` file passed on the command line (file association); auto-discovers companion `.initials`, `.sol`, `coolsolve.conf` and lookup tables. Any other path is refused (403): a client never makes the server read a file of its choice |
| `GET` | `/api/v1/files/bundle` | Download a ZIP bundle (eescode + initials + sol + conf + the files of its `$INCLUDE`s + debug files, §7.2) |
| `POST` | `/api/v1/files/upload` | Upload a multipart `.zip` bundle or individual files; returns `mainFile`, the `.eescode` taken as the model when there are several (§7.2) |
| `POST` | `/api/v1/includes/open` | Open a file that the model includes in a **new** session, for a new window; body: `{"name":"<as in $INCLUDE>"}`, resolved as for the solve (§7.5). Returns `{sessionId, modelName, description}`; 404 with the resolution error |

The file workflow is **ZIP-centric**:
- Open loads files from disk (local mode) or from an uploaded ZIP.
- Save produces a ZIP bundle containing all companion files plus any debug output.
- Models shipped with CoolSolve are opened from the **model library**
  (§4.12), not from a disk folder: the `examples/` folder of the repository is
  only used by the CLI and by the regression tests.

### 4.4 Solver Operations

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/api/v1/solve` | Start an async solve; returns immediately with `{"status":"started"}` |
| `GET` | `/api/v1/solve/result` | Get the last solve result (variables, residuals, timing, blocks) |
| `GET` | `/api/v1/solve/stream` | SSE endpoint: streams `blockStart`, `blockDone`, and `result` events in real time |
| `POST` | `/api/v1/solve/cancel` | Cancel a running solve; returns 409 if no solve is in progress |
| `POST` | `/api/v1/update-guesses` | Copy last `.sol` → `.initials` (update guesses from solution) |

The solve runs in a background thread.  A `ProgressCallback` in `SolverOptions`
pushes events onto the SSE channel between blocks.  The frontend uses the
browser's native `EventSource` API to consume these events.

### 4.5 Parse & Variables

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/api/v1/parse` | Parse source code (no solve); returns equation/variable counts, parse errors and `includes`, the `$INCLUDE` directives of the model as resolved (`{line, name, file, origin: "library"\|"project", libraryId, description, routines, nested}`) |
| `GET` | `/api/v1/variables` | Get solved variable values (name, value, isArray) |
| `GET` | `/api/v1/variables/inferred` | Get solved variables enriched with inferred property and fluid metadata |

`/api/v1/parse` is called with a short debounce as the user types in the
editor, enabling real-time syntax error markers.

### 4.6 Debug Output

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/debug/files` | List debug output files for the current session |
| `GET` | `/api/v1/debug/file?name=...` | Get content of a specific debug file |

A debug solve is triggered by calling `POST /api/v1/solve` with
`"debug": true` in the request body.  This runs the solver with
`enableTracing=true` and calls `runner.generateDebugOutput()`, producing
~15 files (report.md, equations.md, variables.md, analysis.json,
solver_trace.md, incidence.md, tearing.md, evaluator.md, residuals.txt,
equations.tex, etc.) in a session-isolated temp directory.

### 4.7 Thermodynamic Diagrams

The charts are drawn by [CoolPlot](https://github.com/CoolProp/CoolPlot)'s C++ core (coolplot-core) from lines
computed by `CoolProp::Plot`, through the module `src/diagram/` (design and status:
[CoolPlot integration](coolplot_integration.md)).

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/coolprop/fluids` | List all registered fluids with `name`, `type` (`Real`, `IdealGas`, `HumidAir`, `Incompressible`, `Mixture`, `AmmoniaWater`), `hasDome`, `fractionRange` for the solutions (`[0, 0.6]` for `MEG`), and `modelFluids` from the current session |
| `GET` | `/api/v1/diagram/options` | What a chart can be: `diagram_types` (name such as `ph`, `Ts`, `Th`; title; axes `x`, `y` and their log flags; `families` with `key`, `label`, `symbol`, `default_num`), `fluid_types` drawn (`Real`, `Mixture`, `Incompressible`, `IdealGas`, `HumidAir`) and `fluid_kinds` (per type: `diagram_types`, `families`, `saturation`, `fractions`, `default_box`, `default_region`), `unit_systems` (`SI`, `KSI`, `engineering`, `IP` and `coolsolve`, the model's units, each with the `record` to put in a spec), `default_unit_system`, `regions` (`DEF`, `ACHP`, `ORC`, `NONE` with their limits), `themes`, `label_modes`, `marked_points`, `coolplot` (`core_version`, `commit`) |
| `POST` | `/api/v1/diagram/scene` | A chart of the session's model. Body: `fluid` (CoolSolve fluid name), `fractions` and `basis` (the mass fraction of an incompressible solution, `[0.3]`, `"mass"`), `p_Pa` (the total pressure of a humid-air chart, default 101325), `spec` (CoolPlot spec record without `fluid` and `cycles`), `cycles` (cycle definitions naming the model's variables), `figure_size_pt` `[w, h]`, `previous_labels` (optional). Returns `svg` (coolplot-core's document), `axes` (keys, labels, limits in display units, log flags, `x_unit`/`y_unit` with display = SI × `scale` + `offset`, `plot_size_pt`), `labels` (to send back as `previous_labels`), `states` (each cycle state: `cycle`, `number`, `label`, `vars`, `values_SI`), `spec` (as composed, readable by CoolPlot), `diagnostics` (`curves`, `warnings`), `timing_ms` (`curves`, `compose`, `svg`). The CoolProp backend is the model's `coolpropBackend`; the cycles use the last successful solve and are skipped without one. 400 for an unknown fluid, a fluid kind not drawable yet (ideal gases, humid air, incompressibles, `NH3H2O`), an unknown diagram or a spec value coolplot-core refuses; a line that cannot be computed is a diagnostic, not an error |
| `POST` | `/api/v1/diagram/export` | The chart as a file: the body of `/diagram/scene` and `"format"` (`"svg"` or `"pdf"`); returns the file with `Content-Disposition: attachment; filename="<model>_<diagram>.<format>"`. The PDF embeds the TrueType fonts the server finds (500 when there is none: install `fonts-dejavu-core`) |
| `POST` | `/api/v1/diagram/state` | The state under the cursor: `{fluid, fractions, basis, diagram, x_SI, y_SI}` → `{"state": {T, p, h, s, rho, u, v, Q}}` in SI (Q = −1 in single phase) or `{"state": null}` |
| `GET` | `/api/v1/diagram/states` | The states of the last solved model by fluid, grouped by the property calls they share: `{"fluids": {"R134a": [{"number", "label", "vars": {"T": "T_1", ...}}]}}`; computed at solve time, kept by *Back*, cleared by *New* |
| `GET` | `/api/v1/diagram/definitions` | The diagram definitions of the session (`diagrams.json`: `{"format": "coolsolve-diagrams", "version": 1, "diagrams": [...]}`, empty list when none) |
| `PUT` | `/api/v1/diagram/definitions` | Store the diagram definitions of the session; they go into the ZIP bundle, into the snapshot of *Back*, and are cleared by *New*, *Open*, a library model and an upload |

### 4.8 Parametric Study

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/api/v1/parametric` | Start an async parametric sweep; body: `{"sweepVariables":[{"name":"T_in","values":[20,25,30]}]}`; supports 1D (1 variable) and 2D (2 variables, Cartesian product) |
| `GET` | `/api/v1/parametric/result` | Get the result of the last parametric study |

The parametric study runs in a background thread with **warm-starting**:
each grid point uses the nearest previously-converged solution as initial
guesses. For the first point, the current `.initials` file is used. Progress
events are streamed via the existing SSE channel (`/api/v1/solve/stream`).

The sweep works by modifying the source code at each grid point: imposed
variable assignments (e.g. `T_in = 25`) are replaced with the sweep value
using regex substitution. Each point runs a full `CoolSolveRunner::solve()`.

**Parse endpoint enhancement**: `POST /api/v1/parse` now returns additional
fields per variable:
- `isImposed` (boolean): `true` if the variable appears as `var = constant`
- `imposedValue` (number): the constant value (for imposed variables)
- `unitSource` (string): `"code"`, `"inferred"`, or `""` indicating unit origin

**Bundle persistence**: ZIP bundles now include `parametric_studies.json`
containing saved parametric study results for load/restore round-trips.

### 4.9 Lookup Tables

| Method | Endpoint | Description |
|--------|----------|--------------|
| `GET` | `/api/v1/tables` | List all lookup tables in the current session; returns `[{name, columns, rows}]` |
| `GET` | `/api/v1/tables/:name` | Get the raw CSV content of a named table |
| `PUT` | `/api/v1/tables/:name` | Create or replace a table; body is plain `text/csv` |
| `DELETE` | `/api/v1/tables/:name` | Remove a table from the session |

Table names are restricted to `[A-Za-z0-9_-]+`.  Tables are stored in the
session and written to the solver's working directory as `<name>.csv` before
each solve or parametric run.  They are included in the ZIP bundle and
restored on upload.

### 4.10 Integral Table

When a model contains an equation-based `INTEGRAL(...)` call, the solve produces
a time-series trajectory.  It is delivered to the frontend in two forms and rides
on the standard solve response — no separate solve call is needed:

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/integral/result` | Last integral solve as columnar JSON (`integrationVar`, `columns`, `data`, `numRows`, `csvName`, `totalSteps`, `rejectedSteps`); 404 when no integral has been solved |
| `GET` | `/api/v1/integral/csv` | Raw trajectory CSV text (`text/csv`); 404 when none. Exposes the CSV restored from a bundle round-trip where the columnar JSON is not available |

In addition, the `POST /api/v1/solve` `"done"` SSE event and
`GET /api/v1/solve/result` both carry optional `integralTable` (columnar) and
`integralCsvName` fields whenever an integral was solved.  The auto-written
`<model>-integral.csv` is exported in the ZIP bundle (see §7.2) and restored to
the session on upload, so the Integral tab repopulates on a bundle round-trip
(download → upload → same columns and rows).

### 4.11 Usage Log & Statistics

Every GUI solve attempt — a *Solve* click or a *Try Harder* run — appends one
CSV line to a small log file. CLI runs are never logged, parametric studies are
not logged (they are automated batches), and dynamic (`INTEGRAL`) models appear
as ordinary solves.

**Log file location** — resolved once at server start:

1. `ServerOptions::usageLogFile` if set programmatically;
2. else the `COOLSOLVE_GUI_LOG` environment variable;
3. else `coolsolve_gui.log` in the server's working directory.

The resolved path is printed at startup (`Usage log: /path/to/file`). The
directory must be writable by the server process: if it is not, a single
`[Warning] Cannot write the usage log at ...` is emitted, solving continues
unaffected, and the `/stats` dashboard stays empty. Hardened deployments need
an explicit writable path — see
[Deployment Guide §4](deployment_ubuntu_apache.md#4-systemd-service-auto-start),
where `ProtectSystem=strict` makes the working directory read-only and the log
lives in `/var/lib/coolsolve` instead.

The file is git-ignored. It is created with a header line on first write:

```
# timestamp,kind,outcome,duration_ms,model_bytes,equations,max_block,model,ip,version
2026-08-26T12:12:34Z,solve,success,34.891,90,4,1,usage_api_test,127.0.0.1,0.3.0
2026-08-26T12:12:39Z,tryharder,success,1.866,90,4,1,usage_api_test,127.0.0.1,0.3.0
```

| Field | Meaning |
|-------|---------|
| `timestamp` | UTC, ISO 8601 |
| `kind` | `solve` or `tryharder` (Deep Search) |
| `outcome` | `success`, `failed`, `parse_error` |
| `duration_ms` | Wall-clock duration of the whole attempt |
| `model_bytes` | Size of the `.eescode` source ("model length") |
| `equations` / `max_block` | Equation count and largest block (-1 when unknown) |
| `model` | Model name as of the start of the attempt; `(unnamed)` when empty; CSV-quoted when it contains `,` or `"`, with control characters replaced by spaces so a record is always one line |
| `ip` | Client IP — honours `X-Forwarded-For` behind the Apache reverse proxy |
| `version` | CoolSolve version |

**Overhead**: the entry is written once, after the solve has finished, in the
background solve thread (open → append → close, ≈2 µs measured). Nothing is
added to the solver's hot paths and the CLI binary never executes this code.

**Statistics endpoint**

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/stats/log` | Aggregated statistics as JSON: totals per outcome/kind, per-day counts, logarithmic duration histogram (with `medianMs`/`p95Ms`), top models, top clients, unique client count |

The aggregation streams the file in a single pass with bounded memory
(model/IP maps are pruned to a fixed cap; durations live in a fixed-size
logarithmic histogram). Results are cached and only recomputed when the file's
size or modification time changes, so dashboard requests stay O(1) even for
multi-megabyte logs (≈175 ms cold compute for 1 M lines, ≈10 µs cached —
measured on a Linux desktop with a `-O2` build).

**Hidden dashboard**: the SPA serves a statistics page at `/stats`
(`gui/src/components/UsageStats.tsx`). It renders summary cards plus Plotly
charts (attempts per day, outcome split, duration distribution, top models,
top clients). No link points to it anywhere in the UI — it is reachable only
by typing the URL.

### 4.12 Model Library

The CoolSolve Library is embedded in the binary (§9.6) and served by these routes. A model is a
**folder**: code, initial values, solver configuration, lookup tables, solution baseline and a README.

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/library/menu` | The toolbar *Library* menu: the taxonomy (category → sub-category) with the **verified** models only (`{id, name, title, level, kind, functions}`, `functions` being the number of FUNCTIONs and PROCEDUREs that `$INCLUDE <name>` imports), the number of models and the snapshot information (`commit`, `content_hash`, `models`). `available: false` with a `reason` when the build has no library |
| `POST` | `/api/v1/library/open` | Load a model folder into the session; body: `{"id":"CSL-0012"}` (ID or name; optional `"file":"x_coolsolve.eescode"` for another `.eescode` of the folder). **Any status can be opened**, blocked models included. Returns `{libraryId, title, status, modelName, file, variants, hasInitials, hasSol, hasConf, tables}`; 404 for an unknown model or file |
| `GET` | `/library/` | The library dashboard (`library.html`): **all** the models, blocked ones included, with their gaps, functions and a link to open each in CoolSolve |
| `GET` | `/library/<path>` | Any embedded file of the snapshot (README figures, `model.json`, …); only the allow-listed files exist (§9.6): `sources/`, `tools/` and `templates/` are not embedded |
| `GET` | `/docs/raw/library/<path>/README.md` | A model README, rendered by the documentation viewer (`/docs/#library%2Fmodels%2F…%2FREADME.md`) with a banner that opens the model in CoolSolve |

`GET /api/v1/files/eescode` returns `libraryId` (empty when the model does not come from the library), so the
GUI shows a link to the model page of the dashboard next to the model name. *New* and ZIP *Upload* clear it.

---

## 5. Session Management

Session state is stored in a `std::map<string, shared_ptr<Session>>`, keyed
by a 32-hex-char ID.  The session ID is issued as a `coolsolve_session` cookie
on the first request and persists across calls.

Each `Session` holds:
- `eescode`, `initials`, `sol`, `conf` — current file content strings
- `lastResult` — last `SolveResult` from the solver
- `inferredVariables` — enriched variable list (property, fluid, units) extracted after solve
- `modelFluids` — unique fluids detected via variable inference
- `sessionSnapshot` — one-level undo snapshot
- `debugDir` — temporary directory path for debug output files
- `solveInProgress`, `cancelToken` — concurrency control for async solve
- `lookupTableCSVs` — map from table name to raw CSV text (in-memory store for lookup tables)
- `lastIntegralResult` / `lastIntegralCSV` / `lastIntegralCsvName` — last integral-solve trajectory (columnar JSON + raw CSV), guarded by a dedicated mutex; cleared on new/open/reset

A request is attributed to a session in this order: `X-CoolSolve-Session` header, `?sid=` query parameter
(both only if the ID is 32 lowercase hex digits — it names a temp folder — anything else is ignored), then
the cookie; without any, a new cookie session is created. The GUI gives a tab its own ID with
`enableTabSession()` (`gui/src/api/session.ts`, kept in `sessionStorage`, so it survives a reload but not a new
tab); only tabs opened by a library link do it, all the others share the cookie session as before.

A session that has been idle for 24 h and has no solve running is deleted, with its temp folder (checked at most
every 10 minutes), so that tabs opened from the library do not accumulate on a long-running server. A session
costs a few kilobytes (the model texts and the last result) and no CPU until a solve is started. An open tab is
never idle: it calls `POST /api/v1/keepalive` every 15 minutes (`App.tsx`), so only the sessions of closed tabs
expire; and a request never deletes the session it is made for, even when it comes after a long pause.

On a reload, a tab with its own session shows again the model that this session holds (`showSessionModel()` in
`gui/src/utils/sessionModel.ts`: files, lookup tables, last result), with the edits sent to the server by the last
solve or save. A tab of the cookie session starts empty, as before.

Session-scoped temp directories live at `/tmp/coolsolve_sessions/{id}/`.
The "New" action clears all in-memory state while keeping the session alive.

---

## 6. UI Layout

```
┌──────────────────────────────────────────────────────────────────────────┐
│  TOOLBAR                                                                 │
│  [← Back] [New] [📂 Open] [⬇ Save] [▶ Solve] [🐛 Debug] [⏹ Stop]     │
│  [↻ Update Guesses] [{ } Comment] [" " Comment] [☀/☾ Theme]            │
├────────────────────────────┬─────────────────────────────────────────────┤
│                            │                                             │
│  CODE EDITOR (Monaco)      │  RIGHT PANEL (Tabs)                         │
│                            │  ┌──────────────────────────────────────┐   │
│  T_in = 25                 │  │Variables│Arrays│Config│Debug│Diagrams│   │
│                            │  │ (right panel unchanged)               │   │
│  P = 101325                │  ├──────────────────────────────────────┤   │
│  h = enthalpy('Water',     │  │  Variable table (filtered, sortable) │   │
│        T=T_in, P=P)        │  │  ┌──────┬──────┬───────┬──────┐      │   │
│  ...                       │  │  │ Name │Value │Initial│Units │      │   │
│                            │  │  ├──────┼──────┼───────┼──────┤      │   │
│                            │  │  │ T_in │ 25   │  25   │      │      │   │
│                            │  │  │ P    │101325│101325 │      │      │   │
│                            │  │  │ h    │104929│  —    │      │      │   │
│                            │  └──────────────────────────────────────┘   │
├────────────────────────────┴─────────────────────────────────────────────┤
│  BOTTOM PANEL (collapsible, tabbed)                                      │
│  Console │ Parametric │ Lookup Tables │ Integral                         │
└──────────────────────────────────────────────────────────────────────────┘
```

### 6.1 Split Panes

Layout is managed by `SplitPane.tsx`, a reusable resizable split-pane
component (horizontal/vertical, min/max constraints, drag divider).  The
editor/right-panel and main/console splits are both independently resizable.

### 6.2 Toolbar

Implemented in `Toolbar.tsx`.  All actions include keyboard shortcuts:

| Action | Shortcut | Description |
|--------|----------|-------------|
| Back | — | Restore previous session snapshot |
| New | — | Clear editor and start fresh |
| Open | `Ctrl+O` | Load a ZIP bundle (`.eescode` + companion files) |
| Library | — | Drop-down that follows the taxonomy of the model library (category → sub-category → model) with a search box; lists the **verified** models only. Picking one opens the whole model folder (§7.3). The last entry opens the dashboard `/library/` with all the models |
| Save | `Ctrl+S` | Download ZIP bundle |
| Solve | `F2` or `Ctrl+Enter` | Run async solve, stream progress via SSE. `F2` is the EES key; it also works while the editor has the focus (Monaco's own `F2`, *Rename Symbol*, is overridden) |
| Try Harder | `F2` or `Ctrl+Enter` (after a failure) | Re-run with `deepSearch: true` — substitutes `deepSearchPipeline`, forces tearing + symbolic reduction on, and consults `multiStartMode`.  The Solve button morphs into Try Harder automatically when the previous run failed; editing the model, initials, or config reverts it to Solve. |
| Debug Solve | `Ctrl+Shift+Enter` | Solve with tracing enabled, generate debug output |
| Stop | `Escape` | Cancel running solve |
| Update Guesses | `Ctrl+G` | Copy solved values into initial guesses |
| Toggle `{ }` Comment | `Ctrl+/` | Wrap/unwrap selection in brace comments |
| Toggle `" "` Comment | `Ctrl+Shift+/` | Wrap/unwrap selection in double-quote comments |
| Theme | — | Toggle light/dark theme (persisted in `localStorage`) |

### 6.3 Code Editor (Left Pane) — `CodeEditor.tsx`

Monaco Editor with a custom language definition (`languages/ees.ts`):
- **Keywords**: `FUNCTION`, `PROCEDURE`, `END`, `CALL`, `IF`, `THEN`, `ELSE`,
  `ENDIF`, `REPEAT`, `UNTIL`, `GOTO`, `$ifnot`, `$endif`, `$common`, `$include`
- **Built-in functions** (highlighted as `support.function`): `enthalpy`,
  `entropy`, `density`, `pressure`, `temperature`, `quality`, `cp`, `cv`,
  `specheat`, `conductivity`, `viscosity`, `volume`, `soundspeed`,
  `surfacetension`, `sin`, `cos`, `tan`, `exp`, `log`, `log10`, `sqrt`,
  `abs`, `min`, `max`, `ln`
- **Comment styles**: `"..."` (visible/documentation), `{...}` (invisible),
  `//...` (line comment)
- **String literals**: `'R134a'`, `'Water'` (fluid names)
- **Number literals**: integer and floating-point with optional exponent
- Bracket matching, auto-indentation, find/replace, minimap

Real-time parse feedback: a debounced `POST /api/v1/parse` call on every
edit returns parse errors, which are shown as red squiggly underlines via
Monaco `setModelMarkers`.

**Include bar** (`IncludeBar.tsx`): only when the model has an `$INCLUDE`, a
thin bar above the editor lists the included files (the library model name and
ID, or the file name of the folder), from the `includes` of the last parse.
The tooltip gives the routines imported and the nested includes; a click opens
the file in a **new CoolSolve window** (§7.5). A model without `$INCLUDE` shows
no bar.

**Comment toggle** (toolbar buttons + keyboard shortcuts):
- `{ }` invisible comment: wraps selected lines in `{ ... }`.  If already
  wrapped, removes the braces.  Uses Monaco `executeEdits`.
- `" "` visible comment: wraps selection in `"..."`.  Idempotent toggle.

### 6.4 Variables Tab — `VariableTable.tsx`

A filterable HTML table (search box at top) with columns:

| Column | Editable | Description |
|--------|----------|-------------|
| Name | No | Variable name; array vars (`T[1]`, `T[2]`) grouped visually |
| Value | No | Solved value (blank before solving) |
| Initial | **Yes** | Editable initial guess; parses `.initials` text, commits on Enter/blur via `PUT /api/v1/files/initials` |
| Units | No | Unit annotation from inferred metadata (where available) |

### 6.5 Arrays Tab — `ArrayTable.tsx`

A spreadsheet grid for variables of the form `var[i]`:
- Rows = integer indices (sorted numerically)
- Columns = distinct base names (sorted alphabetically)
- Array detection regex: `^[A-Za-z_]\w*\[\d+\]$`
- Filter input at top
- Shows "N arrays × M indices" count
- "Copy to Clipboard" exports as TSV

If no array variables exist, a placeholder message is shown.

### 6.6 Config Tab — `ConfigEditor.tsx`

A form-based editor for `coolsolve.conf` values, organised into collapsible
groups covering all `SolverOptions` keys (~40+ fields):

| Group | Key fields |
|-------|-----------|
| Main Iteration | `maxIterations`, `tolerance`, `relativeTolerance`, `verbose` |
| Line Search | `lsAlpha`, `lsRho`, `lsMaxIterations` |
| Trust Region | `trInitialRadius`, `trMaxRadius`, `trEta` |
| Levenberg-Marquardt | `lmInitialLambda`, `lmLambdaFactor` |
| Partitioned Solver | `partitionedMaxIterations`, `partitionedTolerance` |
| Tearing | `enableTearing`, `tearingMinBlockSize`, `tearingMethod` |
| Pipeline | `solverPipeline` (ordered list), `pipelineMode` (radio) |
| Deep Search Pipeline | `deepSearchPipeline`, `deepSearchPipelineMode` — used by the *Try Harder* button |
| Multi-Start | `multiStartMode` (`always` / `deepsearch` / `never`), `multiStartMaxRestarts`, `multiStartNumCores` |
| Safety | `timeoutSeconds` |
| Output | `writeSolFile`, `writeResiduals` |
| Integration | `integralMethod`, `integralFixedStep`, `integralMaxSteps`, `integralRelTol`, `integralAbsTol`, `integralMinStep`, `integralMaxStep`, `integralRichardson`, `integralOutputInterval` |

Changed values are highlighted.  A "Reset All" button restores defaults.
Changes sync to the session via `PUT /api/v1/files/conf` on "Apply" or "Solve".

### 6.7 Debug Tab — `DebugViewer.tsx`

Displayed after a Debug Solve.  Shows the contents of the debug output folder:
- File list panel (priority-sorted: report.md first, then equations.md,
  variables.md, etc.)
- Content viewer for the selected file
- Refresh button

Files typically generated (~15):
`report.md`, `equations.md`, `variables.md`, `analysis.json`,
`solver_trace.md`, `incidence.md`, `tearing.md`, `evaluator.md`,
`residuals.txt`, `equations.tex`, and block-level debug files.

### 6.8 Diagrams Tab — `components/diagram/`

Property diagrams of the model's fluid with its states, drawn by the server: the tab sends a
diagram definition, the server computes the lines with `CoolProp::Plot` (the pinned CoolProp, the
model's backend) and returns the SVG that CoolPlot's C++ core composes (§4.7,
[CoolPlot integration](coolplot_integration.md)). The tab is located in the right panel alongside
Variables, Arrays, Config, and Debug; it stays mounted (hidden) so that the chart is not redrawn
when the tab comes back.

#### Controls (`DiagramControls.tsx`)

- **Fluid**: from `GET /api/v1/coolprop/fluids`, the model fluids first (variable inference).
  Pure fluids, CoolProp mixtures, incompressibles (with a **Mass fraction** field for the
  solutions, from the model's calls by default) and ideal gases (CoolSolve's, on the EES/JANAF
  basis) and humid air (`AirH2O`: the psychrometric and Mollier h,x charts, with a **total
  pressure** field, the model's by default, and T and W range fields) are drawn; `NH3H2O` is
  not yet (the server answers "not drawable yet").
- **Diagram**: the fluid charts of coolplot-core: log(p)-h, T-s, h-s, p-s, p-ρ, T-ρ, p-T and T-h;
  incompressibles: p-h, p-s, p-ρ, p-T; ideal gases: all but T-h. The families offered follow the
  fluid's kind (no lines of constant quality without a saturation curve).
- **Units**: `engineering` (°C, bar, kJ/kg; the default), `SI`, `KSI`, `IP`, and `coolsolve`, the
  model's units (°C, Pa, J/kg).
- **Region**: `DEF` (relative to the saturation range; the default), `ACHP` (heat pumps; the
  default of mixtures), `ORC`. A fluid without saturation curve shows **T** and **p** range
  fields instead (its default box: 200-2000 K and 0.1-100 bar for an ideal gas, the fluid's
  temperatures and 1-100 bar for an incompressible).
- **Lines**: each family of the chart (isotherms, isentropes, isochores, lines of constant
  quality, …) on or off, with its count and rounded values.
- **Critical point** marked with a dot; **saturation scale** next to the vertical axis (charts
  with p or T vertical); both only for fluids with a saturation curve.
- **Theme**: as the GUI (`default` or `dark`), `print` (black and white, families told apart by
  dash patterns; the report's figures), `coolpack` (the look of CoolPack's log(p)-h charts).
- **Line weight**: *light* (the default: every family of lines at 0.5 pt and 55 % opacity, so
  that the saturation curve and the cycles stand out) or *as the theme* (CoolPlot's 0.75 pt,
  opaque). Stored as CoolPlot's `style_overrides` in the spec.
- **Labels**: `rails`, `edge`, `margin`, `inline`.

The chart is redrawn automatically (debounced) when a control changes, when the panel is resized
and after a solve; there is no *Generate* button. The status line shows the warnings of the server
(lines that could not be computed, states left out) and the time of the lines, the composition and
the SVG.

#### Overlays (`DiagramOverlays.tsx`) and cycles (`CycleEditor.tsx`)

Each overlay or cycle is a cycle definition that names the model's variables; the server reads
their values from the last solve, converts them to SI, completes each state from two of its
properties and draws the paths between them.

- **Inferred states**: the states of the solved model from `GET /api/v1/diagram/states`. The server
  groups the variables of property calls that share two variables of independent properties
  (`h_2 = enthalpy(R134a, P=P_2, s=s_2)` and `T_2 = temperature(R134a, P=P_2, h=h_2)`: one state;
  two calls at the same T and p without quality: one state; a saturation temperature from
  `P_cond` alone: no state) and numbers them from the names (the number most variables carry,
  `T_2`, `h_2`, `s[2]` → 2; otherwise a label such as `su_cp` and the next free number). Drawn as
  numbered markers without lines.
- **Arrays**: two array variables (`h[i]`, `P[i]`) give state i; drawn as a chain in index order,
  closed on request, with paths `auto` (the default) or straight.
- **Cycles** (the collapsible editor): cycles of any variables. Each has a name, states (number,
  optional label, property/variable pairs, the variables of the fluid suggested first), loops
  typed as `1-2-3-4-1; 2-5` (empty: markers only), a path kind for the whole cycle (`auto`:
  isobaric if the pressures are equal, else isenthalpic, isentropic, isothermal, otherwise
  polytropic; or one of these, or straight) and optional per-connection overrides. *From inferred
  states* copies the server's states into an editable cycle with one closed loop. Adding the first
  cycle switches the inferred states off (their markers would sit on the same states); the
  checkbox brings them back.

#### Hover (`useDiagramHover.ts`)

Near a state marker, a tooltip names the cycle and the state and lists its variables with their
values in the chart's units. Elsewhere in the plot area, the state under the cursor
(`POST /api/v1/diagram/state`, debounced, stale requests cancelled): T, p, h, s, ρ, and the
quality inside the dome. Nothing is asked while a zoom is dragged.

#### Zoom and size (`DiagramView.tsx`)

`DiagramView` shows coolplot-core's SVG inline and asks for a figure the size of its container.
Drag a rectangle on the plot area to zoom (a new `view_SI` in the spec: the server lays out the
labels again for the window); double-click or the reset button returns to the whole chart. When only
the size or the zoom changes, the labels of the scene on screen are sent back so that they keep their
places, and no property is computed again. `DiagramView` does not depend on the tab's controls: it
draws any definition.

#### Badge, persistence and export

The badge in a corner of the chart names the equation of state of the fluid
(`R134a: Tillner-Roth and Baehr (1994)`) and the versions (`CoolProp 8.1.0dev (8eb3b77),
CoolPlot 1.5.0`), as CoolPlot writes them; definitions saved before it keep no badge.

The definition is stored in the session (`PUT /api/v1/diagram/definitions`, debounced) and travels
in the ZIP bundle as `diagrams.json` (§7.2); opening a bundle, *Back* and the reload of a library
tab restore it (`showSessionModel()`). The status line has three downloads: **SVG** and **PDF**
(`POST /api/v1/diagram/export`) and **Spec**, the chart spec as composed (with the fluid and the
cycles), which CoolPlot reads (`DiagramSpec.from_dict`) to draw the same chart in Python.

The LaTeX report takes the charts as vector figures: `/api/v1/latex/compile` draws every stored
definition with the last solve in the print theme (`thermo_diagrams.pdf`,
`thermo_diagram_2.pdf`, …). The PNG rasterised in the browser from the SVG (`exportPlots.ts`,
`thermo_diagrams.png`) is the fallback, for the default chart (not stored until it is edited) and
a server without fonts. The same files come from the command line: `coolsolve --diagram
diagrams.json -o chart.pdf model.eescode` (README, *Command Line Options*).

### 6.9 Console (Bottom Panel) — `Console.tsx`

Auto-scrolling log panel with color-coded lines:
- Parse results (equation/variable count, errors)
- Solve progress via SSE (block N/M, iteration, residual norm)
- Final result: SUCCESS / FAIL with timing breakdown
- CoolProp warnings and errors

### 6.10 Lookup Tables Tab — `LookupTableEditor.tsx`

A spreadsheet-style manager for CSV lookup tables in the bottom panel.

#### List view

Displays a table of all tables in the current session (name, row count,
column names).  Populated by `GET /api/v1/tables` automatically after each
successful solve and when the tab is opened.

- **Create**: enter a table name and press Enter or click *Create*; an empty
  2-column template is written to the session via `PUT /api/v1/tables/:name`.
- **Delete**: click the trash icon; prompts for confirmation.
- **Open**: click any row to open the grid editor for that table.

#### Grid editor

An inline editable HTML table:
- Header row cells are editable (rename columns).
- Data cells are plain text inputs (a number, empty for NaN, or text read by `LOOKUP$`; the table is kept as CSV text, so a text cell is shown and saved as typed).
- **+ Row / + Col** buttons append rows and columns.
- **× buttons** (per row/column) delete that row or column (header row
  is protected; at least 1 data row and 1 column are always kept).
- **Import CSV**: loads a `.csv` file from disk, replacing the current grid.
- **Export CSV**: downloads the current grid as a `.csv` file.
- **Save**: calls `PUT /api/v1/tables/:name` with the serialised CSV and
  updates the in-memory store (`modelStore.lookupTableCSVs`).

Tables saved here are automatically used in the next solve — no need to
restart or re-upload files.

### 6.11 Integral Table Tab — `IntegralTable.tsx`

A read-only viewer for the trajectory of an equation-based dynamic (`INTEGRAL`)
model, shown as a sibling **Integral** tab in the bottom panel (alongside
Console / Parametric / Lookup Tables).  It mirrors the Parametric study UX.

The tab resolves its data from two sources, in priority order:

1. The live columnar JSON from the solve response (`integralTable`).
2. The restored CSV text (`integralCSV`), parsed client-side with a small
   CSV splitter that honours quoting — this is the path taken after a pure
   bundle load (no live solve since server restart), where only the CSV
   survives the round-trip.

Contents:

- A scrollable columnar table with the integration variable as the first
  (highlighted) column, followed by the `$IntegralTable` variables.
- A row/step count and, when available, the accepted/rejected step counts
  (RK45, BDF).
- A Y-variable selector + `PlotlyChart` line plot of one tabulated variable
  versus the integration variable.
- An **Export CSV** button that re-serialises the table to a downloadable
  file named after the model (`<csvName>`).

The tab is read-only (EES Integral Tables are not user-editable).  An empty-state
message is shown when no integral has been solved.  It does not auto-switch on
solve — it only ensures the bottom panel is open — so the user keeps console
focus during a solve; the tab populates immediately once clicked.

---

## 7. File Management

### 7.1 Local Mode

CoolSolve works directly with files on disk:

```
my_model/
├── my_model.eescode        # Source code
├── my_model.initials        # Initial guesses
├── my_model.sol             # Solution (generated)
├── coolsolve.conf           # Solver configuration
└── my_model_coolsolve/      # Debug output (generated)
```

- **Open** (`POST /api/v1/files/open`, for the file passed on the command line):
  reads `.eescode` and auto-discovers `.initials`, `.sol`, `coolsolve.conf` and
  the lookup tables in the same directory.
- **Save** (`GET /api/v1/files/bundle`): downloads a ZIP bundle named
  `<model_name>_coolsolve.zip`, containing eescode + initials + sol + conf +
  debug files.  The browser's File System Access API (`showSaveFilePicker`)
  is used when available, with a download-prompt fallback.

The folder *is* the project — no proprietary format is needed.

### 7.2 ZIP Bundle Format

The ZIP uses a minimal uncompressed implementation (custom CRC32-based ZIP
writer/reader in `server.cpp` — no external library):

| Entry | Contents |
|-------|----------|
| `<name>.eescode` | Current editor content |
| `<name>.initials` | Current initial guesses |
| `<name>.sol` | Last solution (if solve was successful) |
| `coolsolve.conf` | Solver configuration |
| `<tablename>.csv` | One entry per lookup table in the session |
| `<name>-integral.csv` | Integral trajectory (when an `INTEGRAL` model was solved) |
| `parametric_studies.json` | Parametric study results (§4.8) |
| `diagrams.json` | Diagram definitions of the Diagrams tab (§6.8): fluid, CoolPlot spec, cycles as variable names |
| `<other>.eescode`, `<other>-<table>.csv` | Every file the model includes (`$INCLUDE`, nested includes too; library models as their `.eescode` file) with its lookup tables, and the other `.eescode` files of an uploaded project |
| `debug/*` | Debug output files (if present) |

**Included files.** Because a file of the model folder takes precedence over a library model of the same
name, the bundle runs on its own: unzipped and run with the CLI, or uploaded again, the model includes the
copies of the bundle, i.e. the library version it was saved with (also in a build without library). On
upload, the model is the `.eescode` file named like the ZIP, else the one that no other file of the ZIP
includes (preferring one that has its `.initials` or `.sol`); the other `.eescode` files become the
session's project files, which the `$INCLUDE`s of the model find.

### 7.3 Library Models (a Model Is a Folder)

A model of the [CoolSolve Library](https://github.com/IntSusEnergySystems/CoolSolve_Library) is a project folder
like §7.1: `<name>.eescode`, `<name>.initials`, `coolsolve.conf`, lookup tables `<name>-<table>.csv`, the
`<name>.sol` baseline, plus `README.md`, `model.json` and `figures/`. The binary carries those files; opening a
model (`POST /api/v1/library/open`, `LibraryStore::loadModel`) loads **everything the session understands** exactly
as *Open* does for a folder on disk: the table name is what follows `<name>-` in the CSV file name, and the
`.initials`/`.sol` of the model are used by its variants (`x_coolsolve.eescode`). The README and the figures stay
on the server side and are shown by the documentation viewer. Because everything lives in the session, the usual
*Save* (ZIP bundle) round-trips a library model like any other.

*Open in CoolSolve* links (dashboard, README banner) point to `/?library=CSL-0012[&file=…]`. The GUI then
gives the tab its own session (§5), opens the model, and removes the parameter from the address: a reload shows
the tab's session, edits included, instead of opening the model again. The tab shows an info icon next to the
model name that links back to the model page of the dashboard.

### 7.4 Session Snapshot (Back Button)

Before each destructive operation (New, Open, Upload), the current session
state is saved to a `SessionSnapshot`.  The `POST /api/v1/back` endpoint
restores it.  This provides one-level undo.

### 7.5 Included Files (`$INCLUDE`)

The `$INCLUDE`s of the model ([Language Reference §7](language_reference.md)) are resolved by the server
exactly as by the CLI (`includeResolverFor()` in `server.cpp`, `include/coolsolve/include_resolver.h`): a file
`<name>.eescode` of the **project folder** first, then the library model `<name>`. The project folder of a
session is:

| The model was opened… | Project folder |
|---|---|
| from disk (command line, file association) | the folder of that file, on disk: only `<name>.eescode` files of that folder are read |
| from the library | the folder of the library model in the snapshot |
| from a ZIP, or *New* | the other `.eescode` files of the ZIP (none for a new model) |

The parse, the solve and the parametric study use the same resolution, so the editor markers and the
solve agree.

**Including a library model.** In the *Library* menu, the models that define functions or procedures have a
copy button next to their title: it copies `$INCLUDE <name>` to the clipboard, to paste in the editor. The model
page of the dashboard (`/library/#CSL-…`) shows the same line with a *Copy* button.

**Opening an included file.** Each file of the include bar is a plain link to `/?include=<name>` with
`target="_blank"` (no script pop-up for a blocker to stop). The click notes the session of the window in
`localStorage` (`coolsolve_include_opener`, valid one minute); the new window calls `POST /api/v1/includes/open`
**as that session** (its `X-CoolSolve-Session`, or the cookie), which creates a new session holding the file (a
library model with its folder, a file of the project with its companion files and the same project folder);
the window adopts it (`adoptTabSession`, §5), shows its model and removes the parameter from the address. Both
windows are independent: edits of the included file in its window do not change the other window (a file on
disk is read again at each solve; files of a ZIP project are copies).

---

## 8. CoolProp Integration

### 8.1 Warmup

On server start, a background thread calls `warmupCoolProp()`:

```
Server start
  └── warmupCoolProp()          // ~1–3 s
      └── PropsSI("T", "P", 101325, "Q", 0, "Water")
```

The frontend polls `GET /api/v1/health` and shows "Initialising thermodynamic
engine…" until `coolpropReady` is `true`.

### 8.2 During Solving

CoolProp `AbstractState` instances are cached per fluid (via
`CoolPropConfig::cacheEnabled`).  Each solve task runs in its own background
thread with its own `CoolSolveRunner` instance — no global state is shared.

### 8.3 Variable Inference for Diagrams

`inferVariables(IR&)` is called as step 3 of every `runner.run()` (regardless
of whether debug mode is on).  It populates `VariableInfo.inferredProperty`
(T/P/H/S/D/Q/V/L/C, and `VOL` for the specific volume, which CoolProp has no
code for: `V` is the viscosity) and `VariableInfo.inferredFluid` for every variable that
appears as an argument to a thermodynamic function.

After the solve thread completes, the server extracts this data into the
session's `inferredVariables` and `modelFluids` fields (~0.1 ms overhead).
This extraction is the only new per-solve cost — negligible.

### 8.4 Fluid Types

The `FluidRegistry` (`src/fluids.cpp`) registers 90+ real fluids, the ideal
gases, humid air, the incompressible solutions (glycols, brines, `LiBr`,
`SeaWater`) and, built from their string, CoolProp mixtures.  The Diagrams tab
draws fluids of type `FluidType::Real`, CoolProp mixtures (`FluidType::Mixture`),
incompressibles, ideal gases and humid air (`AirH2O`; the variable inference still names it
`Air_ha`, `CS-BUG-HUMIDAIR-NAME`); `NH3H2O` is not drawn yet
([CoolPlot integration](coolplot_integration.md) §5).

---

## 9. Build & Packaging

### 9.1 Asset Embedding

`cmake/embed_assets.cmake` runs at CMake configure time:

1. Reads all files under `gui/dist/`
2. For each file, generates a `const unsigned char[]` array in
   `build/embedded_assets.cpp` with the correct MIME type
3. Generates a `lookupEmbeddedAsset(path)` function in
   `build/embedded_assets.h`

The server's catch-all route (`GET /.*`) serves embedded assets when built
with `COOLSOLVE_EMBEDDED_ASSETS` defined, or reads from the filesystem at
`COOLSOLVE_GUI_DEV_DIR` during development.

### 9.2 Development Mode (Hot-Reload)

```bash
# Terminal 1: start backend
cd build && cmake .. && make -j$(nproc) coolsolve
./coolsolve --gui --no-browser
# Backend at http://localhost:8550

# Terminal 2: start Vite dev server
cd gui && npm install && npm run dev
# Frontend at http://localhost:5173
# /api, /docs and /library are proxied to port 8550 via vite.config.ts
```

### 9.3 Production Build (Single Binary)

```bash
# Build frontend
cd gui && npm run build
# Produces gui/dist/ (~5 MB JS + ~10 kB CSS)

# Rebuild backend — CMake embeds gui/dist/ into the binary
cd build && cmake .. && make -j$(nproc) coolsolve

# Run
./coolsolve --gui
# Opens browser at http://localhost:8550
```

### 9.4 CMake Options

| Option | Default | Effect |
|--------|---------|--------|
| `COOLSOLVE_BUILD_GUI` | `ON` | Compile `server.cpp`, link cpp-httplib, embed frontend assets |
| `COOLSOLVE_EMBED_LIBRARY` | `ON` | With the GUI: fetch and embed the CoolSolve Library (§9.6) |
| `COOLSOLVE_LIBRARY_REPO` / `COOLSOLVE_LIBRARY_TAG` | GitHub repository / `main` | Where the library is fetched from, and which branch, tag or commit |
| `FETCHCONTENT_SOURCE_DIR_COOLSOLVELIBRARY` | — | Use a local library checkout instead of fetching |
| `FETCHCONTENT_FULLY_DISCONNECTED` | `OFF` | Never touch the network (uses what `.fetchcontent_cache/` already holds) |

When `OFF`, the binary is identical to the pre-GUI version in all respects.

### 9.5 SPA Fallback Routing

The server's catch-all route (`GET /.*`) returns `index.html` (HTTP 200) for
any non-API path not found in embedded assets.  This supports client-side
routing if it is ever added.

### 9.6 Library Embedding

The model library is embedded by `cmake/embed_library.cmake` (shared helpers in `cmake/embed_common.cmake`):

- **Where it comes from**: fetched with `FetchContent` from `COOLSOLVE_LIBRARY_REPO`
  (`IntSusEnergySystems/CoolSolve_Library`), exactly like CoolProp, into `.fetchcontent_cache/`. The default
  `COOLSOLVE_LIBRARY_TAG` is `main`, so **every `cmake` configure pulls the latest pushed library** (`build_gui.sh`
  and `build_installer.bat` always configure). If GitHub is unreachable the copy already fetched is used, with a
  warning. A release pins a tag (`-DCOOLSOLVE_LIBRARY_TAG=<tag>`).
- **Working on the library at the same time**: `-DFETCHCONTENT_SOURCE_DIR_COOLSOLVELIBRARY=/path/to/CoolSolve_Library`
  embeds that checkout, uncommitted changes included (`coolsolve --version` then prints `<commit>+local-changes`).
- **Re-embedding at every build**: the `coolsolve_library_snapshot` target runs the script at each
  `cmake --build`; it hashes the files to embed and rewrites `embedded_library.cpp` only when the hash changed,
  so an unchanged library costs 0.1 s and no recompilation.
- **What is embedded** (allow-list, mirror of `EMBED_FILES` in the library's `tools/build_index.py`):
  `library.html`, `library.json`, `taxonomy.json`, `docs/taxonomy.md` and, in each model folder, `README.md`,
  `model.json`, `coolsolve.conf`, `*.eescode`, `*.initials`, `*.sol`, `<name>-<table>.csv`,
  `figures/*.{png,svg,jpg,jpeg,webp}`. The script warns when the number of files differs from the one in
  `library.json` (stale index or diverging allow-lists).
- **Cost**: ≈ 770 files, 7 MB → a 38 MB generated source that compiles in ≈ 6 s with ≈ 1 GB of memory, and
  ≈ 8 MB more in the binary. `-DCOOLSOLVE_EMBED_LIBRARY=OFF` builds without library (the Library menu then says so).
- **Identity**: `coolsolve --version` prints `CoolSolve Library <commit>, snapshot <hash>`; the menu response
  carries the same information.

---

## 10. Frontend Source Structure

```
gui/
├── package.json                   # Dependencies: react, @monaco-editor/react, zustand, lucide-react, plotly.js-basic-dist-min
├── tsconfig.json / tsconfig.app.json / tsconfig.node.json
├── vite.config.ts                 # Vite config; /api proxy → port 8550
├── eslint.config.js
├── index.html                     # HTML entry point
├── public/
│   └── favicon.svg
└── src/
    ├── main.tsx                   # React entry point
    ├── App.tsx                    # Root layout: toolbar + resizable split panes + tabs
    ├── App.css                    # Full CSS with light/dark theme variables
    ├── plotly.d.ts                # Plotly type declarations
    ├── vite-env.d.ts
    ├── api/
    │   ├── client.ts              # Typed fetch wrappers for all API endpoints
    │   ├── session.ts             # Per-tab session ID (X-CoolSolve-Session header / ?sid=)
    │   └── types.ts               # TypeScript interfaces for all API responses (incl. parametric types)
    ├── stores/
    │   ├── modelStore.ts          # Zustand: eescode, initials, sol, conf, parse errors, solve result, parametric studies, lookup tables, integral table, user unit overrides
    │   └── uiStore.ts             # Zustand: theme, active tabs (incl. 'lookuptables', 'integral'), panel visibility (persisted to localStorage)
    ├── languages/
    │   └── ees.ts                 # Monaco Monarch language definition for CoolSolve syntax
    ├── utils/
    │   ├── diagrams.ts            # Diagram definitions: empty document, debounced PUT to the session
    │   ├── download.ts            # downloadBlob(): save a file of the server (diagram downloads)
    │   ├── library.ts             # openLibraryModel(): load a library model folder into the session; dashboard URLs
    │   └── sessionModel.ts        # showSessionModel(): editor state read back from the session (open, Back, upload, reload)
    └── components/
        ├── Toolbar.tsx            # Top toolbar with all actions + keyboard shortcuts
        ├── LibraryMenu.tsx        # "Library" drop-down: taxonomy tree of the verified library models, search, link to /library/
        ├── CodeEditor.tsx         # Monaco editor wrapper with syntax highlighting + debounced parse
        ├── VariableTable.tsx      # Variable table with units column (color-coded), imposed badges, filter modes
        ├── ArrayTable.tsx         # Spreadsheet view for array variables (var[i])
        ├── ConfigEditor.tsx       # Form-based coolsolve.conf editor (9 groups, 30+ fields)
        ├── DebugViewer.tsx        # Debug output file list + content viewer
        ├── Console.tsx            # Auto-scrolling console with color-coded log lines
        ├── diagram/               # Diagrams tab (§6.8): coolplot-core's SVG drawn by the server
        │   ├── DiagramTab.tsx     # Options, the definition edited, the request
        │   ├── DiagramControls.tsx # Fluid, diagram, units, region, lines, theme, labels
        │   ├── DiagramOverlays.tsx # Inferred states and arrays as cycle definitions
        │   ├── CycleEditor.tsx    # The user's cycles: states, loops, paths
        │   ├── useDiagramHover.ts # Tooltips: state markers, the state under the cursor
        │   ├── options.ts         # Options loaded once, display units, value format
        │   ├── DiagramView.tsx    # Reusable view: fetch, inline SVG, zoom, resize
        │   ├── definition.ts      # Spec and cycle builders
        │   └── svgCoords.ts       # Pointer → display units → SI (log axes)
        ├── ParametricStudy.tsx    # Parametric sweep: variable selector, range inputs, 1D/2D plots, results table
        ├── LookupTableEditor.tsx  # Lookup table manager: list view + inline editable grid + CSV import/export
        ├── IntegralTable.tsx      # Read-only trajectory viewer for INTEGRAL models: table + Plotly plot + CSV export
        ├── UsageStats.tsx         # Solve usage statistics dashboard (hidden /stats page)
        ├── SplitPane.tsx          # Reusable resizable split-pane (horizontal/vertical)
        └── PlotlyChart.tsx        # Thin Plotly wrapper component
```

---

## 11. Verified Test Results

All the following tests pass against the current codebase:

### 11.1 Backend / API (integration tests in `build/check.sh`)

| Test | Result |
|------|--------|
| Backend compilation (C++17) | ✅ Pass |
| All C++ unit tests (708 assertions, 94 cases) | ✅ Pass |
| `GET /api/v1/health` | ✅ `{"status":"ok","coolpropReady":true}` |
| `POST /api/v1/parse` with source code | ✅ Returns equation/variable counts |
| `POST /api/v1/files/open` with a `.eescode` file | ✅ Loads file + discovers companion files |
| Library menu / open / per-tab sessions / bundle round-trip / `/library/` pages / `/files/open` restricted to the command-line file | ✅ `tests/test_library_api.py --auto` (67 checks) |
| `F2` in headless Chromium: Solve and Try Harder reach `POST /api/v1/solve`, in and outside the editor; no solve on `Shift`/`Ctrl`/`Alt+F2`, an empty model or a running solve | ✅ `tests/test_gui_shortcuts.py --auto` (16 checks) |
| `POST /api/v1/solve` with rankine1 example | ✅ Full solve result with block details |
| `GET /api/v1/variables` | ✅ Returns solved variable values |
| `GET /api/v1/solve/stream` (SSE, curl) | ✅ 30/30 blocks streamed; final result JSON embedded |
| `POST /api/v1/solve` async return | ✅ Returns `{"status":"started"}` immediately |
| `POST /api/v1/files/save-as` | ✅ Creates file at new path with companion files |
| `POST /api/v1/solve/cancel` (no solve) | ✅ Returns 409 "No solve is in progress" |
| `POST /api/v1/solve/cancel` (during solve) | ✅ Cancels at block 79/112; SSE shows "Solve cancelled by user" |
| Session cookie | ✅ `Set-Cookie: coolsolve_session=...` on first request |
| Session isolation | ✅ New session gets empty state; original retains files |
| ZIP bundle download | ✅ Valid ZIP, 4 files, 6972 bytes |
| ZIP bundle upload | ✅ Re-uploaded ZIP; all companion files detected |
| Debug files (before solve) | ✅ `{"files":[]}` |
| Debug files (after debug solve) | ✅ 15+ files generated |
| Debug file content (`report.md`) | ✅ Full analysis content |
| `GET /api/v1/coolprop/fluids` | ✅ 128 fluids; Water `hasDome:true`; Air `hasDome:false` |
| `GET /api/v1/coolprop/fluids` (modelFluids, before parse) | ✅ `"modelFluids":[]` |
| Diagram options, scenes (R134a p-h, errors), definitions round trip, `diagrams.json` in the ZIP bundle, New/Back, a solved model drawn with its states, paths, the state under the cursor and the states of the model, incompressibles, ideal gases, a tabular backend, blends, the badge, the export route and the report's PDF figures | ✅ `tests/test_diagram_api.py --auto` (120 checks) |
| Diagrams tab in headless Chromium: the chart drawn on load at the container size, the default spec, the states of a solved R134a cycle as markers, zoom (`view_SI`, labels sent back) and reset, the report figure `thermo_diagrams.png`, a cycle from the editor drawn along its paths, the hover on a marker and elsewhere, an ideal gas with its cycle, a solution at two concentrations, a blend, the badge and the SVG/PDF/Spec downloads | ✅ `tests/test_gui_diagrams.py --auto` (52 checks) |
| `GET /api/v1/variables/inferred` (before solve) | ✅ 0 variables |
| `GET /api/v1/variables/inferred` (after solve) | ✅ 39 vars with T/P/H/S/Q for R134a |
| `GET /api/v1/coolprop/fluids` (modelFluids after solve) | ✅ `"modelFluids":["R134a"]` |
| modelFluids cleared on new model | ✅ Reset to `[]` after `POST /api/v1/new` |

### 11.2 Frontend

| Test | Result |
|------|--------|
| TypeScript type-check (zero errors) | ✅ Pass |
| Production build (`npm run build`) | ✅ ~5 MB JS (includes plotly.js) + 9.4 kB CSS in `gui/dist/` |
| Embedded assets (index.html, JS, CSS, SVG) | ✅ Served with correct MIME types |
| SPA fallback routing | ✅ Unknown paths return `index.html` (HTTP 200) |
| React app load from embedded binary | ✅ Pass (end-to-end) |

---

## 12. Deployment Architecture

### 12.1 Local Mode

```
coolsolve --gui [port] [--no-browser]
```

Opens the default browser at `http://localhost:<port>` after a short warmup.
`--no-browser` suppresses automatic browser launch (useful for headless
environments or when tunnelling).

### 12.2 Online Mode (Reverse Proxy)

```
Browser ──HTTPS──► nginx/Caddy (reverse proxy) ──HTTP──► coolsolve --gui
```

A Docker image packages the binary with the embedded frontend:

```dockerfile
FROM ubuntu:24.04 AS runtime
COPY coolsolve /usr/local/bin/coolsolve
EXPOSE 8550
ENTRYPOINT ["coolsolve", "--gui", "8550"]
```

### 12.3 Cross-Platform Notes

| Concern | Solution |
|---------|----------|
| Browser opening | `xdg-open` (Linux), `open` (macOS), `start` (Windows) |
| File paths | Normalised with `std::filesystem`; REST API always uses forward slashes |
| Port conflicts | Default 8550; if busy, prints error and exits — auto-increment is a future improvement |
| Binary distribution | Static linking of CoolProp + cpp-httplib + embedded frontend; single binary for Linux x86_64, macOS arm64/x86_64, Windows x86_64 |

### 12.4 Security (Online Deployment)

| Threat | Mitigation |
|--------|------------|
| Arbitrary code execution | CoolSolve only parses `.eescode` — the PEG grammar is the sandbox; no `eval` or shell execution |
| Resource exhaustion | Per-session `timeoutSeconds`, session expiry on idle |
| File system access | Online sessions use isolated temp directories named after validated session IDs; `/files/open` only opens the file passed on the command line (none on a server) |
| Input size | Max upload size enforced by cpp-httplib |

---

## 13. Coding-Agent Workflow

End-to-end verification without a browser:

```bash
# 1. Build backend
cd build && cmake -DCOOLSOLVE_BUILD_GUI=ON .. && make -j$(nproc)

# 2. Run C++ unit tests
./coolsolve_tests

# 3. Build frontend
cd ../gui && npm ci && npm run build

# 4. Re-embed assets and rebuild binary
cd ../build && make -j$(nproc)

# 5. Start server in background
./coolsolve --gui 8551 &
SERVER_PID=$!
sleep 4  # CoolProp warmup

# 6. Run API integration tests
bash check.sh

# 7. Stop server
kill $SERVER_PID
```

Every step has a binary pass/fail exit code.  Steps 1–4 are decoupled and
can be run independently.

---

## 14. Future Work

The following items are not yet implemented.  They are listed roughly in
priority order.

### 14.1 LaTeX Report Generation

`POST /api/v1/report` endpoint to:
- Generate a `.tex` report from the current model (variables, equations,
  solve results) using the existing `runner.generateDebugOutput()` data
- Download as `.tex` source or, if `pdflatex` is available on the server,
  as a compiled PDF
- Frontend: "LaTeX Report" button in the toolbar with format selector

**Open question**: server-side `pdflatex` vs. client-side WASM LaTeX engine
(e.g. SwiftLaTeX) for PDF compilation.

### 14.2 Online Deployment Hardening

| Item | Description |
|------|-------------|
| Session expiry | Automatically destroy sessions after N minutes of inactivity; clean up temp dirs |
| Docker image | Publish a minimal Docker image (`ubuntu:24.04` base with static binary) |
| nginx config | Reference reverse-proxy config with TLS, rate limiting, max body size |
| Authentication | OAuth (GitHub/Google) for a multi-user hosted service; cookie sessions are sufficient for a public demo |

### 14.3 Sensitivity Analysis — Remaining Work

The core parametric sensitivity analysis feature is now implemented (see
sections below). Remaining enhancements:

| Item | Description |
|------|-------------|
| CSV export | Download sweep results as CSV |
| Multi-output plots | Multiple output variables on the same 1D chart |
| 3D surface plots | `plotly.js-basic-dist-min` → `plotly.js-dist-min` for 3D surfaces (currently using heatmap/contour) |
| Progress bar | Visual progress bar during parametric study (currently shows text status) |
| Cancel button | Abort a running parametric study via `POST /api/v1/parametric/cancel` |

#### Implemented Features

- **Imposed variable detection**: Parse endpoint detects `var = constant`
  patterns and reports `isImposed`, `imposedValue` per variable
- **Units column**: Variable table shows units with color coding — blue
  (code-specified), green (inferred from CoolProp), orange (user-specified)
- **Parametric sweep API**: `POST /api/v1/parametric` accepts 1 or 2 sweep
  variables, runs grid with warm-starting in background thread, sends progress
  via SSE, results via `GET /api/v1/parametric/result`
- **Parametric UI**: Bottom panel "Parametric" tab with sweep variable
  selectors, range inputs, 1D line plots and 2D contour/heatmap plots (Plotly),
  results table with swept variables highlighted
- **Bundle persistence**: Parametric study results are saved/loaded in ZIP
  bundles (`parametric_studies.json`)

### 14.4 Thermodynamic Diagram Enhancements

The Diagrams tab is built on CoolPlot (§6.8); its phases are in
[CoolPlot integration](coolplot_integration.md) §5. Done in phase 1: units, isolines, region
presets, themes, zoom, persistence in the ZIP bundle; in phase 2: cycle paths, the state
editor, states grouped from the property calls, the hover; in phase 3: blends on several
threads, incompressibles, tabular backends, ideal gases, humid air; in phase 4: the badge, SVG/PDF
downloads, the report's vector figures, *Export chart spec*, `coolsolve --diagram`. Still to do:

| Item | Description |
|------|-------------|
| Badge control | Choose the references (none, short, full) and the versions line in the tab |
| Multiple diagrams | Several definitions per model (the document holds a list; the tab edits the first) |

### 14.5 AG Grid Variable Table

Replace the current plain HTML table in `VariableTable.tsx` with AG Grid
Community for:
- Sorting and multi-column filtering
- Column pinning and reordering
- Faster rendering for models with hundreds of variables
- In-cell editing with validation

### 14.6 Advanced GUI Features

| Item | Description |
|------|-------------|
| Auto-save | In local mode, auto-save editor content every 30 s to `.coolsolve.autosave`; offer recovery on next startup |
| Variable quick-jump | Click variable in table → highlight in editor; `Ctrl+Click` variable in editor → popup showing value/units/determining equation |
| Equation dependency graph | D3.js or Mermaid graph of block structure: which blocks depend on which, block sizes, converged/failed status |
| Smart initial guesses | Highlight variables with no initial guess in orange; suggest ranges for thermodynamic properties based on fluid and context |
| Diff view | After solving, side-by-side diff of initial guesses vs. converged values; large jumps highlighted |
| Unit conversion helper | Small utility panel for common unit conversions (°C↔K, bar↔Pa, kJ↔J) |
| Export to Python/Julia | Generate a Python or Julia script replicating the model using CoolProp bindings |
| Responsive layout | On narrow screens (tablets), right panel moves below editor; bottom panel collapses to a status bar |
| Port auto-increment | If default port 8550 is busy, try 8551, 8552, … and print the actual URL |

### 14.7 Testing Infrastructure

| Item | Description |
|------|-------------|
| Vitest unit tests | `VariableTable.test.tsx`, `ArrayTable.test.tsx`, `CodeEditor.test.tsx`, `ConfigEditor.test.tsx`, `stores/modelStore.test.ts`, `utils/arrayVars.test.ts` |
| Playwright E2E tests | `solve-flow.spec.ts`, `file-roundtrip.spec.ts`, `comment-toggle.spec.ts`, `config-edit.spec.ts` |
| Example sweep script | `tests/test_all_examples_api.sh` — solve every example in `examples/` via the REST API |
| CI pipeline | GitHub Actions workflow: backend build + unit tests, frontend build + unit tests, API integration sweep |

### 14.8 Notebook-Style Project Format

For future sensitivity and plotting features, a richer `.coolsolve` project
format (JSON) could bundle:
- CoolSolve source code
- Initial guesses and solver configuration
- A list of "analysis cells" (sensitivity sweeps, diagram configs)
- Cached results (last solution, sweep data)

This is deferred until the sensitivity analysis feature is designed in detail.

### 14.9 Long-Term / Speculative

| Item | Description |
|------|-------------|
| WASM build | Compile CoolSolve + CoolProp to WebAssembly — run entirely in the browser with no server |
| Multi-model comparison | Solve two `.eescode` files side by side |
| Collaborative editing | Online shared sessions via CRDT or OT |
| Mobile support | Responsive breakpoints for tablet use (code editor experience on phones is inherently limited) |
