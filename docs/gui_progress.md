# CoolSolve GUI — Implementation Progress

This document tracks the current state of the GUI implementation.
Refer to `docs/gui_plan.md` for the full architecture and design specification.

---

## Architecture Overview

The GUI is a **single-page web application** (React + TypeScript) served by an
embedded HTTP server (cpp-httplib) inside the CoolSolve binary. Running
`coolsolve --gui` starts the server and opens the browser — no external
dependencies are needed.

For development, the Vite dev server (port 5173) proxies API calls to the
CoolSolve backend (port 8550), enabling hot-reload. For production, the
frontend is compiled and embedded directly into the binary via
`cmake/embed_assets.cmake`.

---

## What Is Done

### Backend (C++)

| Component | Status | Notes |
|-----------|--------|-------|
| cpp-httplib integration | ✅ Done | FetchContent in CMake, guarded by `COOLSOLVE_BUILD_GUI` option |
| `include/coolsolve/server.h` | ✅ Done | `ServerOptions` struct + `startServer()` declaration |
| `src/server.cpp` — REST API | ✅ Done | All Phase 1 endpoints implemented (see table below) |
| `--gui [port]` flag in `main.cpp` | ✅ Done | Default port 8550, `--no-browser` flag available |
| CoolProp warmup on startup | ✅ Done | Background thread, `/api/v1/health` reports readiness |
| CORS middleware | ✅ Done | For dev mode (Vite on different port) |
| Cross-platform browser opening | ✅ Done | `xdg-open` / `open` / `start` |
| Asset embedding | ✅ Done | `cmake/embed_assets.cmake` generates `embedded_assets.cpp/h` at configure time |
| SPA fallback | ✅ Done | Non-API 404s serve `index.html` for client-side routing |

#### Implemented API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/health` | Health check + CoolProp readiness |
| `POST` | `/api/v1/warmup` | Trigger CoolProp warmup |
| `GET` | `/api/v1/files/eescode` | Get current EES source code |
| `PUT` | `/api/v1/files/eescode` | Set EES source code |
| `GET` | `/api/v1/files/initials` | Get initial guesses |
| `PUT` | `/api/v1/files/initials` | Set initial guesses |
| `GET` | `/api/v1/files/sol` | Get solution file content |
| `GET` | `/api/v1/files/conf` | Get config |
| `PUT` | `/api/v1/files/conf` | Set config |
| `POST` | `/api/v1/files/open` | Open file from disk (discovers companion files) |
| `POST` | `/api/v1/files/save` | Save files to disk |
| `POST` | `/api/v1/parse` | Parse EES code, return errors + variable info |
| `POST` | `/api/v1/solve` | Async solve — returns immediately, runs in background thread |
| `GET` | `/api/v1/solve/result` | Get last solve result |
| `GET` | `/api/v1/solve/stream` | Real-time SSE progress events (chunked content provider, per-block callbacks) |
| `GET` | `/api/v1/variables` | Get solved variable values |
| `POST` | `/api/v1/files/save-as` | Save files to a new path (creates parent dirs, updates session) |
| `POST` | `/api/v1/update-guesses` | Copy .sol → .initials |
| `POST` | `/api/v1/solve/cancel` | Cancel a running solve |
| `GET` | `/api/v1/examples` | List example files |
| `GET` | `/api/v1/files/bundle` | Download ZIP bundle of eescode+initials+sol+conf |
| `POST` | `/api/v1/files/upload` | Upload multipart .zip or individual files |
| `GET` | `/api/v1/debug/files` | List debug output files for current session |
| `GET` | `/api/v1/debug/file` | Get content of a specific debug file (`?name=...`) |
| `GET` | `/api/v1/variables/inferred` | **Phase 4** — Solved variables with inferred property/fluid/units |
| `GET` | `/api/v1/coolprop/fluids` | **Phase 4** — List all fluids with type, hasDome, + model fluids |
| `POST` | `/api/v1/coolprop/saturation` | **Phase 4** — Compute saturation dome (liquid+vapor T,P,H,S,D + critical point) |

### Frontend (React + TypeScript + Vite)

| Component | Status | Notes |
|-----------|--------|-------|
| Vite project scaffold | ✅ Done | React 18 + TypeScript, Vite 5 |
| Vite API proxy | ✅ Done | `/api` → backend on port 8550 |
| API client (`api/client.ts`) | ✅ Done | Typed fetch wrappers for all endpoints |
| API types (`api/types.ts`) | ✅ Done | TypeScript interfaces for all API responses |
| Zustand model store | ✅ Done | EES code, initials, sol, conf, parse errors, solve result, console |
| Zustand UI store | ✅ Done | Theme, tabs, panel visibility, localStorage persistence |
| EES Monaco language (`languages/ees.ts`) | ✅ Done | Keywords, thermodynamic functions, comment styles, bracket matching, folding |
| `CodeEditor.tsx` | ✅ Done | Monaco editor with EES language, debounced parse-on-type, error markers |
| `VariableTable.tsx` | ✅ Done | Filterable table showing solved variable values + editable Initial column |
| `Console.tsx` | ✅ Done | Auto-scrolling console with color-coded log lines |
| `Toolbar.tsx` | ✅ Done | Open, Save, Examples, Solve, Debug, Update Guesses, Theme toggle; keyboard shortcuts (Ctrl+Enter, Ctrl+S, Ctrl+O, Ctrl+G) |
| `ThermoDiagram.tsx` | ✅ Done | **Phase 4** — Plotly diagrams (T-s, P-h, h-s, T-h), fluid selector with model fluids prioritized, saturation dome, state point overlay, array path overlay with close-cycle, PNG/SVG export |
| `SplitPane.tsx` | ✅ Done | Reusable resizable split pane component (horizontal/vertical, min/max constraints, drag divider) |
| `App.tsx` + `App.css` | ✅ Done | Root layout with resizable split panes, toolbar, editor pane, right panel (tabbed: Variables, Arrays, Config, Debug, Diagrams), bottom panel, status bar, light/dark theme |
| Production build | ✅ Done | ~5 MB JS (includes plotly.js) + 9.4 kB CSS in `gui/dist/` |

### Build & Packaging

| Component | Status | Notes |
|-----------|--------|-------|
| `cmake/embed_assets.cmake` | ✅ Done | Reads `gui/dist/*`, generates C++ byte arrays with MIME types, builds lookup function |
| Embedded asset serving | ✅ Done | `COOLSOLVE_EMBEDDED_ASSETS` code path in `server.cpp` |
| Filesystem fallback serving | ✅ Done | `COOLSOLVE_GUI_DEV_DIR` code path for dev mode |
| Single-binary deployment | ✅ Done | `coolsolve --gui` serves everything from one executable |

### Testing

| Test | Status | Result |
|------|--------|--------|
| Backend compilation | ✅ Pass | C++17 compatible, no errors |
| `/api/v1/health` | ✅ Pass | `{"status":"ok","coolpropReady":true}` |
| `/api/v1/parse` with EES code | ✅ Pass | Returns equation/variable counts |
| `/api/v1/files/open` with example file | ✅ Pass | Loads file + discovers companion files |
| `/api/v1/solve` with rankine1 example | ✅ Pass | Returns full solve result with block details |
| `/api/v1/variables` | ✅ Pass | Returns solved variable values |
| Frontend TypeScript type-check | ✅ Pass | Zero errors |
| Frontend production build | ✅ Pass | Vite build succeeds |
| Embedded assets (all paths) | ✅ Pass | index.html, JS, CSS, SVG all served with correct MIME types |
| SPA fallback routing | ✅ Pass | Unknown paths return index.html (200) |
| End-to-end: server + GUI in browser | ✅ Pass | React app loads from embedded binary |
| SSE solve streaming (curl) | ✅ Pass | 30/30 blocks streamed with start/done events, final result JSON embedded |
| Async `POST /solve` | ✅ Pass | Returns `{"status":"started"}` immediately, solve runs in background |
| All C++ unit tests | ✅ Pass | 708 assertions in 94 test cases |
| `/api/v1/files/save-as` | ✅ Pass | Creates file at new path with companion files |
| `/api/v1/solve/cancel` (no solve) | ✅ Pass | Returns 409 "No solve is in progress" |
| `/api/v1/solve/cancel` (during solve) | ✅ Pass | Cancels at block 79/112, SSE stream shows "Solve cancelled by user" |
| Session cookie | ✅ Pass | `Set-Cookie: coolsolve_session=...` on first request, persists across calls |
| Session isolation | ✅ Pass | New session (no cookie) gets empty state; original session retains loaded files |
| ZIP bundle download | ✅ Pass | Valid ZIP with 4 files (eescode, initials, sol, conf), 6972 bytes |
| ZIP bundle upload | ✅ Pass | Re-uploaded ZIP, all companion files detected and loaded |
| Debug files (before solve) | ✅ Pass | Returns `{"files":[]}` |
| Debug files (after debug solve) | ✅ Pass | 15+ files generated (report.md, equations.md, variables.md, analysis.json, etc.) |
| Debug file content (report.md) | ✅ Pass | Full analysis content with session-isolated temp path |
| Phase 4: CoolProp fluids list | ✅ Pass | 128 fluids returned with correct type/hasDome properties |
| Phase 4: Water fluid properties | ✅ Pass | type=Real, hasDome=true |
| Phase 4: Model fluids initially empty | ✅ Pass | `modelFluids: []` before parse |
| Phase 4: Water saturation dome | ✅ Pass | 50 points, critical T=647.1 K, P=22064000 Pa, ~67 ms |
| Phase 4: R134a saturation dome | ✅ Pass | 100 points, critical T=374.2 K |
| Phase 4: Invalid fluid error | ✅ Pass | Returns HTTP 400 |
| Phase 4: Ideal gas error | ✅ Pass | Air returns HTTP 400 (no dome) |
| Phase 4: Inferred vars empty | ✅ Pass | 0 vars before solve |
| Phase 4: Inferred vars after solve | ✅ Pass | 39 vars with properties {T, P, H, S, Q} for R134a |
| Phase 4: Model fluids after solve | ✅ Pass | `modelFluids: ["R134a"]` |
| Phase 4: Cleared on new model | ✅ Pass | Both inferredVars and modelFluids reset to empty |
| Phase 4: Frontend TypeScript | ✅ Pass | Zero errors (excluding pre-existing Toolbar.tsx unused var) |
| Phase 4: Frontend production build | ✅ Pass | Vite build succeeds |

---

## What Remains To Be Done

### Phase 1 — MVP Completion

These items are needed to make the GUI genuinely usable for interactive work:

| Task | Priority | Description |
|------|----------|-------------|
| ~~SSE solve progress~~ | ~~High~~ | ✅ Done — `ProgressCallback` in `SolverOptions`, async solve in background thread, chunked SSE via `set_chunked_content_provider`, frontend uses `EventSource` |
| ~~Resizable split panes~~ | ~~Medium~~ | ✅ Done — `SplitPane.tsx` component with drag resize, horizontal/vertical modes, integrated into `App.tsx` for editor/panel and main/console splits |
| ~~Editable initial guesses~~ | ~~Medium~~ | ✅ Done — `VariableTable.tsx` parses initials text, shows editable Initial column, inline editing with commit on Enter/blur, syncs via `PUT /files/initials` |
| ~~File Open dialog (browser)~~ | ~~Medium~~ | ✅ Done — File System Access API (`showOpenFilePicker`) with fallback to `<input type="file" multiple>`, companion file discovery for .initials and .conf |
| ~~Save As~~ | ~~Medium~~ | ✅ Done — `POST /api/v1/files/save-as` endpoint + File System Access API (`showSaveFilePicker`) with prompt fallback, creates parent dirs, ensures `.eescode` extension, saves companion files |
| ~~Comment toggle~~ | ~~Low~~ | ✅ Done — Monaco editor actions for `{ }` (invisible) and `" "` (visible) comments with `Ctrl+/` and `Ctrl+Shift+/` keybindings, toolbar buttons |
| ~~Stop/cancel solve~~ | ~~Low~~ | ✅ Done — `POST /api/v1/solve/cancel` endpoint, `cancelToken` in `SolverOptions`, checked between blocks in solver loop, Stop button (red, shown only while solving), `Escape` shortcut |

### Phase 2 — Online Deployment & Polish

| Task | Description |
|------|-------------|
| ~~Session management~~ | ✅ Done — Cookie-based `SessionManager` class, `coolsolve_session` cookie with 32-hex-char IDs, `std::map<string, shared_ptr<Session>>`, session-isolated temp dirs at `/tmp/coolsolve_sessions/{id}/`, all 17+ handlers updated |
| ~~ZIP bundle upload/download~~ | ✅ Done — Custom minimal uncompressed ZIP creation/extraction with CRC32, `GET /files/bundle` creates ZIP of eescode+initials+sol+conf, `POST /files/upload` accepts multipart .zip or individual files |
| ~~Config editor tab~~ | ✅ Done — `ConfigEditor.tsx` with full schema (9 groups × 30+ fields matching all `SolverOptions` keys), collapsible groups, inline editing, boolean dropdowns, Reset All button, syncs via `PUT /files/conf` |
| ~~Debug output viewer tab~~ | ✅ Done — `DebugViewer.tsx` with file list panel + content viewer, priority-sorted files, refresh button, integrated with `runner.generateDebugOutput()` (~15 files: report.md, equations.md, variables.md, analysis.json, solver_trace.md, etc.) |
| ~~Array variables tab~~ | ✅ Done — `ArrayTable.tsx` spreadsheet grid parsing `var[i]` from solve results, columns = base names (sorted alpha), rows = indices (sorted numeric), filter input, "N arrays × M indices" count |
| LaTeX report generation | `/api/v1/report` endpoint + download as `.tex` or PDF |
| AG Grid for variable table | Replace the simple HTML table with AG Grid (sorting, filtering, column pinning, in-cell editing) |
| Deployment config | Docker image, nginx reverse proxy, session expiry |

### Phase 3 — Sensitivity Analysis

| Task | Description |
|------|-------------|
| Sweep API | Identify imposed variables, define sweep ranges, batch solve |
| Sweep UI | Config form, result table, Plotly line charts |
| CSV export | Export sweep results |
| Multi-variable sweeps | 2D heatmaps |

### Phase 4 — Thermodynamic Diagrams ✅ Done

**Implementation complete.** All backend endpoints, frontend component, and
integration tests are implemented and passing.

**Summary of implementation:**
- 3 new REST endpoints (`/coolprop/fluids`, `/coolprop/saturation`, `/variables/inferred`)
- `ThermoDiagram.tsx` component with Plotly charts (T-s, P-h, h-s, T-h)
- Fluid selector with model fluids (from variable inference) prioritized
- Saturation dome generation: logarithmic T spacing, ~200 points, ~0.1-0.3s
- State point overlay: solved scalar variables plotted on the diagram
- Array path overlay: array variables plotted as cycle path with close-cycle option
- 12 integration tests passing, 708 unit tests passing
- Frontend builds cleanly with TypeScript

#### Phase 4 — Detailed Design (reference)

Phase 4 adds interactive thermodynamic property diagrams to the GUI.
The design is **user-driven** — diagrams are never computed automatically.
The user explicitly selects a fluid and diagram type, then clicks "Generate".
This ensures **zero overhead** on standard parse/solve operations.

---

#### 4.0  Design Principles

1. **On-demand only** — no CoolProp calls for diagrams unless the user
   requests them.  Standard `POST /solve` is completely unaffected.
2. **Saturation dome first** — isolines (isotherms, isobars, iso-quality)
   are deferred to a later iteration.
3. **Leverage existing inference** — CoolSolve already runs
   `inferVariables(IR&)` as step 3 of every `runner.run()` call
   (see `src/runner.cpp` lines 50–51).  This populates
   `VariableInfo.inferredProperty` (T/P/H/S/D/Q/V/L/C) and
   `VariableInfo.inferredFluid` (e.g. "Water") for every variable that
   appears as an argument to a thermodynamic function.  This data is
   available **on every run**, not just debug runs.
4. **Minimal new dependencies** — `react-plotly.js` + `plotly.js` on the
   frontend.  No new C++ dependencies.

---

#### 4.1  CoolProp Saturation Curve — How It Works

CoolProp does **not** provide a pre-computed saturation curve.  Each point
must be obtained via `PropsSI()`, but this is fast (~0.05–0.2 ms per call)
and is done on-demand only.

**Algorithm to build the saturation dome for a fluid:**

```
1.  Retrieve critical and triple-point properties:
      T_crit  = PropsSI("Tcrit",    "", 0, "", 0, fluid)
      P_crit  = PropsSI("Pcrit",    "", 0, "", 0, fluid)
      T_triple = PropsSI("T_triple", "", 0, "", 0, fluid)

2.  Generate N temperature points (e.g. N=200) logarithmically spaced
    from T_triple + ε  to  T_crit − ε  (avoid exact endpoints where
    CoolProp may throw).

3.  For each T_i, compute liquid (Q=0) and vapor (Q=1) properties:
      P_liq[i]  = PropsSI("P", "T", T_i, "Q", 0, fluid)
      H_liq[i]  = PropsSI("H", "T", T_i, "Q", 0, fluid)
      S_liq[i]  = PropsSI("S", "T", T_i, "Q", 0, fluid)
      D_liq[i]  = PropsSI("D", "T", T_i, "Q", 0, fluid)
      P_vap[i]  = PropsSI("P", "T", T_i, "Q", 1, fluid)
      H_vap[i]  = PropsSI("H", "T", T_i, "Q", 1, fluid)
      S_vap[i]  = PropsSI("S", "T", T_i, "Q", 1, fluid)
      D_vap[i]  = PropsSI("D", "T", T_i, "Q", 1, fluid)

4.  For pure fluids, P_liq == P_vap and T values are identical by
    definition.  The dome is traced by concatenating:
      liquid branch (low-S to crit) + critical point + vapor branch (crit to low-S)

5.  Return all arrays as JSON.  The frontend picks the appropriate pair
    for the selected diagram type (T vs S, P vs H, H vs S, etc.).
```

Each dome generation requires ~200 × 8 = 1,600 `PropsSI` calls, taking
roughly **0.1–0.3 s** total — acceptable for a user-initiated action.

**Important**: only fluids of type `Real` (i.e. `FluidType::Real` from
`FluidRegistry`) have saturation curves.  Ideal gases, humid air,
incompressibles, and mixtures do not.  The endpoint must filter for this.

---

#### 4.2  New Backend Endpoints

Three new endpoints under `/api/v1/coolprop/`:

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/coolprop/fluids` | List all registered fluids with type, CoolProp name, and whether they have a saturation dome (`hasDome = type == Real`) |
| `POST` | `/api/v1/coolprop/saturation` | Compute the saturation dome for a given fluid.  Body: `{ "fluid": "R134a", "nPoints": 200 }`.  Returns liquid + vapor arrays for T, P, H, S, D.  Also returns `Tcrit`, `Pcrit`, `Hcrit`, `Scrit`. |
| `GET` | `/api/v1/variables/inferred` | Return solved variables enriched with inference data: `{ name, value, inferredProperty, inferredFluid, units, isArray }[]`.  Only available after a successful solve. |

**`GET /api/v1/coolprop/fluids` response:**

```json
{
  "fluids": [
    { "name": "Water",   "coolpropName": "Water",   "type": "Real",      "hasDome": true  },
    { "name": "R134a",   "coolpropName": "R134a",   "type": "Real",      "hasDome": true  },
    { "name": "Air",     "coolpropName": "Air",     "type": "IdealGas",  "hasDome": false },
    { "name": "AirH2O",  "coolpropName": "AirH2O",  "type": "HumidAir",  "hasDome": false }
  ],
  "modelFluids": ["Water", "R134a"]
}
```

The `modelFluids` array contains fluids detected in the current model via
variable inference (unique values from `VariableInfo.inferredFluid`).
These are highlighted in the frontend dropdown for easy selection.
If no model has been parsed yet, this array is empty.

**`POST /api/v1/coolprop/saturation` response:**

```json
{
  "fluid": "R134a",
  "critical": { "T": 374.21, "P": 4059280, "H": 389860, "S": 1562.3, "D": 511.9 },
  "triplePoint": { "T": 169.85 },
  "liquid": {
    "T": [170.0, 172.5, ...],
    "P": [389.6, 467.2, ...],
    "H": [71560, 73980, ...],
    "S": [428.1, 442.7, ...],
    "D": [1529.0, 1524.1, ...]
  },
  "vapor": {
    "T": [170.0, 172.5, ...],
    "P": [389.6, 467.2, ...],
    "H": [381200, 382400, ...],
    "S": [1889.0, 1882.1, ...],
    "D": [2.37, 2.82, ...]
  },
  "nPoints": 200,
  "computeTime_ms": 185.3
}
```

All values are in **SI units** (K, Pa, J/kg, J/(kg·K), kg/m³).  The
frontend will convert to display units if needed (future work).

**`GET /api/v1/variables/inferred` response:**

```json
{
  "variables": [
    { "name": "T_su[1]", "value": 353.15, "inferredProperty": "T", "inferredFluid": "Water", "units": "[K]", "isArray": true },
    { "name": "P_su[1]", "value": 101325, "inferredProperty": "P", "inferredFluid": "Water", "units": "[Pa]", "isArray": true },
    { "name": "h_su[1]", "value": 334940, "inferredProperty": "H", "inferredFluid": "Water", "units": "[J/kg]", "isArray": true },
    { "name": "T_ex",    "value": 293.15, "inferredProperty": "T", "inferredFluid": "R134a", "units": "[K]", "isArray": false }
  ]
}
```

---

#### 4.3  Backend Implementation Details

##### 4.3.1  Storing Inference Data in the Session

Currently, the `runner` object is created inside the solve thread and
destroyed when the thread ends (see `src/server.cpp` lines 850–970).
Only `session.lastResult` (a `SolveResult`) and `session.lastTiming`
are preserved.  The IR (with its inference data) is **not** persisted.

**Change required**: add a new field to `Session`:

```cpp
struct InferredVariable {
    std::string name;
    double value;
    std::string inferredProperty;  // "T", "P", "H", "S", "D", "Q", etc.
    std::string inferredFluid;     // "Water", "R134a", etc.
    std::string units;             // "[K]", "[Pa]", etc.
    bool isArray;
};
std::vector<InferredVariable> inferredVariables;
std::vector<std::string> modelFluids;  // unique fluids in the model
```

At the end of the solve thread (after `session.lastResult = runner.getSolveResult()`
and before the runner goes out of scope), extract inference data:

```cpp
// Extract inference data from IR (runner is still alive here)
const auto& ir = runner.getIR();
session.inferredVariables.clear();
std::set<std::string> fluidSet;
for (const auto& [name, info] : ir.getVariables()) {
    if (Constants::isConstant(name)) continue;
    auto it = runner.getSolveResult().variables.find(name);
    if (it == runner.getSolveResult().variables.end()) continue;
    InferredVariable iv;
    iv.name = name;
    iv.value = it->second;
    iv.inferredProperty = info.inferredProperty;
    iv.inferredFluid = info.inferredFluid;
    iv.units = info.units;
    iv.isArray = (name.find('[') != std::string::npos);
    session.inferredVariables.push_back(iv);
    if (!info.inferredFluid.empty()) {
        fluidSet.insert(info.inferredFluid);
    }
}
session.modelFluids.assign(fluidSet.begin(), fluidSet.end());
```

This extraction is a lightweight map traversal with no CoolProp calls.

##### 4.3.2  Saturation Endpoint Implementation

The endpoint creates a temporary CoolProp evaluation context:

```cpp
svr.Post("/api/v1/coolprop/saturation", [](const Request& req, Response& res) {
    auto body = json::parse(req.body);
    std::string fluidName = body["fluid"];
    int nPoints = body.value("nPoints", 200);

    // Validate fluid
    auto fluid = FluidRegistry::getFluid(fluidName);
    if (!fluid || fluid->getType() != FluidType::Real) {
        // return 400 error
    }
    std::string cpFluid = fluid->getCoolPropName();

    // Get critical and triple point
    double Tcrit = CoolProp::PropsSI("Tcrit", "", 0, "", 0, cpFluid);
    double Pcrit = CoolProp::PropsSI("Pcrit", "", 0, "", 0, cpFluid);
    double Ttriple = CoolProp::PropsSI("T_triple", "", 0, "", 0, cpFluid);

    // Sweep T from Ttriple+1 to Tcrit-0.01
    // For each T: compute Q=0 and Q=1 properties (P, H, S, D)
    // Collect into arrays, skip points where CoolProp throws
    // Return JSON with liquid{}, vapor{}, critical{}, triplePoint{}
});
```

The endpoint is **stateless** — it does not depend on any session or
solved model.  It could even be cached (same fluid always gives the
same dome), though caching is not required for the initial implementation.

##### 4.3.3  Fluids Endpoint Implementation

```cpp
svr.Get("/api/v1/coolprop/fluids", [&](const Request& req, Response& res) {
    auto sessionPtr = getSession(req, res);
    json j;
    json fluidsArr = json::array();
    for (const auto& fluid : FluidRegistry::getAllFluids()) {
        fluidsArr.push_back({
            {"name", fluid->getName()},
            {"coolpropName", fluid->getCoolPropName()},
            {"type", fluidTypeToString(fluid->getType())},
            {"hasDome", fluid->getType() == FluidType::Real}
        });
    }
    j["fluids"] = fluidsArr;
    j["modelFluids"] = sessionPtr->modelFluids;
    res.set_content(j.dump(), "application/json");
});
```

---

#### 4.4  Frontend Component — `ThermoDiagram.tsx`

A new tab in the **right panel** (alongside Variables, Arrays, Config,
Debug).  The right panel is preferred over the bottom panel because:
- Diagrams need significant 2D space; the right panel can accommodate
  a square-ish plot with controls above it.
- The right panel already hosts result-viewing tabs (Variables, Arrays),
  and diagrams are another way to view results.
- The bottom panel (Console) is narrow and text-oriented.

**Component structure:**

```
┌─────────────────────────────────────────────┐
│  Diagrams                                    │
├─────────────────────────────────────────────┤
│  Fluid: [Water      ▼]  Type: [T-s  ▼]     │
│  [Generate Diagram]                          │
│                                              │
│  ☑ Overlay solved states (Water: 12 vars)   │
│  ☐ Overlay array path   X: [s ▼]  Y: [T ▼] │
│                                              │
│  ┌─────────────────────────────────────────┐ │
│  │                                         │ │
│  │          Plotly chart area               │ │
│  │                                         │ │
│  │     ╭───────────╮                       │ │
│  │    ╱              ╲     ● T_su[1]       │ │
│  │   │                │    ● T_ex[1]       │ │
│  │   │   Saturation   │    ──── cycle      │ │
│  │   │     Dome       │                    │ │
│  │    ╲              ╱                     │ │
│  │     ╰───────────╯                       │ │
│  │                                         │ │
│  └─────────────────────────────────────────┘ │
│  [Export PNG]  [Export SVG]                   │
└─────────────────────────────────────────────┘
```

##### 4.4.1  Controls

- **Fluid dropdown**: populated from `GET /coolprop/fluids`.  Model fluids
  (from `modelFluids`) appear first in a "Model Fluids" group, separated
  from the full "All Fluids" list.  Only fluids with `hasDome: true` are
  selectable (others are shown grayed out with a tooltip "No saturation
  curve for ideal gases / humid air").

- **Diagram type dropdown**: T-s, P-h, h-s, T-h.  Each type maps an
  x-axis and y-axis property:

  | Type | X-axis | Y-axis | X-label | Y-label |
  |------|--------|--------|---------|---------|
  | T-s  | S      | T      | s [J/(kg·K)] | T [K] |
  | P-h  | H      | P      | h [J/kg] | P [Pa] |
  | h-s  | S      | H      | s [J/(kg·K)] | h [J/kg] |
  | T-h  | H      | T      | h [J/kg] | T [K] |

  The y-axis for P-h uses **logarithmic scale** (standard convention).

- **Generate Diagram** button: calls `POST /coolprop/saturation` with the
  selected fluid and renders the dome.  Shows a spinner during computation.
  The button is disabled while the request is in flight.

##### 4.4.2  State Overlay (Checkbox: "Overlay solved states")

When checked, the component calls `GET /api/v1/variables/inferred` and:

1. Filters variables matching the selected fluid
   (`inferredFluid == selectedFluid`).
2. Groups scalar variables by their `inferredProperty`.  For the selected
   diagram type (e.g. T-s), it needs both an S value and a T value to
   plot a point.
3. **Grouping heuristic for state points**: variables that share the same
   index suffix (e.g. `T_su[1]` + `s_su[1]` + `h_su[1]`) or the same
   base naming pattern (e.g. `T_1` + `s_1` + `h_1`) are treated as
   belonging to the same thermodynamic state.  The grouping algorithm:
   - For array variables (`var[i]`): group by index `i`.
   - For scalar variables: attempt name-pattern matching (strip known
     prefixes like T_, P_, h_, s_ and match the remainder).
   - Each group becomes a scatter point on the diagram if it has the
     two properties needed for the selected axes.
   - Points are labeled with their index or matched suffix.
4. If a group is missing one axis property but has two other known
   properties (e.g. has T and P but needs S), optionally compute the
   missing property via a new endpoint or a frontend CoolProp.js call.
   **For the initial implementation, only plot points where both axis
   properties are directly available from inference.**

**Visual treatment**: state points are rendered as Plotly markers
(circles, size 10, colored distinctly from the dome curve, with text
labels showing the index/name).

##### 4.4.3  Array Path Overlay (Checkbox: "Overlay array path")

When checked, the user selects X and Y base variable names from
dropdowns populated by the array table's base names (parsed from
`lastResult.variables`).  For example:

- X = `s` (maps to `s[1]`, `s[2]`, ..., `s[N]`)
- Y = `T` (maps to `T[1]`, `T[2]`, ..., `T[N]`)

The component:
1. Extracts paired values `(X[i], Y[i])` for all matching indices.
2. Sorts by index order.
3. Plots as a Plotly scatter+line trace (markers + connected line).
4. **Close cycle** option (checkbox): if checked, connects the last
   point back to the first (common for thermodynamic cycles like
   Rankine, ORC, refrigeration).
5. Points are labeled with their index (`1`, `2`, ..., `N`).

This is independent of the inference system and uses raw solved values,
so it works even for variables without inferred fluid/property metadata.

##### 4.4.4  Export

- **Export PNG**: uses Plotly's built-in `Plotly.downloadImage()`.
- **Export SVG**: uses `Plotly.downloadImage({ format: 'svg' })`.
- Both include the dome, overlays, labels, and legend in the export.

---

#### 4.5  Data Flow Summary

```
  User clicks "Generate Diagram"
       │
       ▼
  POST /coolprop/saturation { fluid, nPoints }
       │     (backend: N × 8 PropsSI calls, ~0.2s)
       ▼
  Response: liquid{T,P,H,S,D}, vapor{T,P,H,S,D}, critical{}
       │
       ▼
  Frontend: pick X/Y arrays for selected diagram type
       │     (e.g. T-s → X=S, Y=T from liquid+vapor)
       ▼
  Plotly renders saturation dome
       │
       ▼
  If "Overlay states" is checked:
       │
       ▼
  GET /variables/inferred
       │     (returns {name, value, property, fluid, units}[])
       ▼
  Filter by fluid → group into state points → plot markers
       │
       ▼
  If "Overlay array path" is checked:
       │     (no API call — uses lastResult.variables from store)
       ▼
  Extract X[i], Y[i] pairs → plot connected line + markers
```

---

#### 4.6  Computational Overhead Analysis

| Operation | When it runs | Cost | Impact on standard solve |
|-----------|-------------|------|------------------------|
| `inferVariables(IR&)` | Every `runner.run()` call (step 3) | ~1–5 ms (map traversal, no CoolProp calls) | **Already runs** — no new cost |
| `initializeVariables(IR&)` | Every `runner.run()` call (step 3) | ~10–50 ms (some CoolProp calls for guesses) | **Already runs** — no new cost |
| Inference data extraction to session | End of solve thread | ~0.1 ms (map copy) | Negligible |
| `POST /coolprop/saturation` | User clicks "Generate" only | ~100–300 ms (1600 PropsSI calls) | **None** — on demand |
| `GET /variables/inferred` | User checks "Overlay" only | ~0.1 ms (read stored data) | **None** — on demand |
| Plotly rendering | Frontend only | ~50–200 ms | **None** — client-side |

**Conclusion**: The only new cost per standard solve is copying inference
data to the session (~0.1 ms), which is negligible.  All diagram-related
CoolProp calls happen on-demand in separate endpoints.

---

#### 4.7  Implementation Steps

| Step | Description | Files Modified/Created |
|------|-------------|----------------------|
| 4.1 | Add `InferredVariable` struct and session fields | `src/server.cpp` |
| 4.2 | Extract inference data at end of solve thread | `src/server.cpp` (solve handler) |
| 4.3 | Implement `GET /api/v1/variables/inferred` endpoint | `src/server.cpp` |
| 4.4 | Implement `GET /api/v1/coolprop/fluids` endpoint | `src/server.cpp` |
| 4.5 | Implement `POST /api/v1/coolprop/saturation` endpoint | `src/server.cpp` |
| 4.6 | Add Plotly.js dependency | `gui/package.json` |
| 4.7 | Create `ThermoDiagram.tsx` component with controls | `gui/src/components/ThermoDiagram.tsx` |
| 4.8 | Add API client functions for new endpoints | `gui/src/api/client.ts`, `gui/src/api/types.ts` |
| 4.9 | Add Diagrams tab to right panel | `gui/src/App.tsx`, `gui/src/stores/uiStore.ts` |
| 4.10 | Implement saturation dome rendering | `ThermoDiagram.tsx` |
| 4.11 | Implement state point overlay | `ThermoDiagram.tsx` |
| 4.12 | Implement array path overlay | `ThermoDiagram.tsx` |
| 4.13 | Implement PNG/SVG export | `ThermoDiagram.tsx` |
| 4.14 | Integration tests | `build/check.sh` |

---

#### 4.8  Required Tests

##### 4.8.1  Backend Integration Tests (curl-based, added to `build/check.sh`)

| # | Test | Method | What it validates |
|---|------|--------|-------------------|
| 1 | **Fluids list** | `GET /coolprop/fluids` | Returns JSON array with ≥90 fluids; each has `name`, `type`, `hasDome`; Water has `hasDome: true`; Air has `hasDome: false` |
| 2 | **Model fluids after parse** | `PUT /files/eescode` + `POST /parse` + `GET /coolprop/fluids` | After parsing a model that uses Water and R134a, `modelFluids` contains both |
| 3 | **Saturation dome — Water** | `POST /coolprop/saturation {"fluid":"Water","nPoints":50}` | Returns 200 OK; `liquid.T` and `vapor.T` arrays have length 50; `critical.T` ≈ 647.1 K (±0.1); `liquid.S[0] < vapor.S[0]`; all values are finite |
| 4 | **Saturation dome — R134a** | `POST /coolprop/saturation {"fluid":"R134a","nPoints":10}` | Returns 200 OK; critical T ≈ 374.2 K; arrays have length 10 |
| 5 | **Saturation dome — invalid fluid** | `POST /coolprop/saturation {"fluid":"FakeFluid"}` | Returns 400 with error message |
| 6 | **Saturation dome — ideal gas** | `POST /coolprop/saturation {"fluid":"Air"}` | Returns 400 with error "No saturation curve for ideal gas" |
| 7 | **Inferred variables — before solve** | `GET /variables/inferred` | Returns empty array or 404 (no solve has been performed) |
| 8 | **Inferred variables — after solve** | Solve rankine1 example, then `GET /variables/inferred` | Returns non-empty array; at least one variable has `inferredProperty: "T"` and `inferredFluid: "Water"`; all entries have `name` and `value` fields |
| 9 | **Inferred variables — property codes** | After solving rankine1, check that at least T, P, H, S properties appear | Validates that `inferVariables()` correctly maps thermodynamic function outputs |
| 10 | **Saturation dome values are SI** | `POST /coolprop/saturation {"fluid":"Water","nPoints":5}` | All `liquid.T` values are > 273 (Kelvin, not Celsius); all `liquid.P` values are > 100 (Pascal, not kPa) |
| 11 | **Model fluids cleared on new model** | `POST /new` + `GET /coolprop/fluids` | `modelFluids` is empty after creating a new model |

##### 4.8.2  Frontend Manual Test Cases

These can be verified visually during development and later converted to
automated tests (Playwright or Cypress):

| # | Test | Steps | Expected result |
|---|------|-------|-----------------|
| F1 | Dome renders | Select Water + T-s → Generate | Saturation dome appears as a closed curve (bell shape) |
| F2 | Diagram type switch | Switch from T-s to P-h → Generate | Plot axes change; P-h has log y-axis |
| F3 | Model fluids priority | Load a model using R134a | R134a appears at top of fluid dropdown under "Model Fluids" |
| F4 | State overlay | Solve rankine1, check "Overlay states" | State points appear as labeled markers on the T-s dome |
| F5 | Array path overlay | Solve orc_simple, check "Overlay array path", X=s, Y=T | Connected line with markers overlaid on dome |
| F6 | Close cycle | Check "Close cycle" checkbox | Last point connects back to first |
| F7 | Export PNG | Click "Export PNG" | Browser downloads a PNG image of the diagram |
| F8 | Ideal gas blocked | Select Air in dropdown | "Generate" button is disabled or Air is grayed out |
| F9 | No dome without Generate | Switch to Diagrams tab | Empty area with instructions, no automatic computation |
| F10 | Spinner during generation | Click Generate for a fluid | Spinner shown while API call is in flight |

##### 4.8.3  Unit Tests (C++ — optional but recommended)

| # | Test | What it validates |
|---|------|-------------------|
| U1 | `FluidRegistry::getAllFluids()` includes Water, R134a, Ammonia | Fluid registration is correct |
| U2 | `FluidRegistry::getFluid("Water")->getType() == FluidType::Real` | Type classification is correct |
| U3 | `FluidRegistry::getFluid("Air")->getType() == FluidType::IdealGas` | Ideal gas excluded from dome |
| U4 | `PropsSI("Tcrit", "", 0, "", 0, "Water")` returns ~647.1 | CoolProp critical point access works |
| U5 | `PropsSI("S", "T", 373.15, "Q", 0, "Water")` returns finite value | Saturation property call works |

---

#### 4.9  Key Code References

These are the existing code locations relevant to Phase 4:

| Location | Relevance |
|----------|-----------|
| `src/variable_inference.cpp` lines 113–200 (`inferVariables()`) | Populates `inferredProperty` and `inferredFluid` for each variable by scanning thermodynamic function calls |
| `src/variable_inference.cpp` lines 57–76 (`getPropMapping()`) | Maps function names to property codes: temperature→T, pressure→P, enthalpy→H, entropy→S, density→D, quality→Q, etc. |
| `include/coolsolve/ir.h` lines 36–58 (`VariableInfo` struct) | Defines `inferredFluid`, `inferredProperty`, `units` fields |
| `src/runner.cpp` lines 50–51 | `inferVariables()` and `initializeVariables()` called as step 3 of every `runner.run()` — confirms inference runs on every solve, not just debug |
| `src/runner.cpp` lines 200–231 | `generateDebugOutput()` writes `variables.md` with the Inferred column — shows format of inference data |
| `src/evaluator.cpp` line 101 (`timedPropsSI()`) | CoolProp wrapper used for all real fluid property evaluations |
| `src/fluids.cpp` lines 11–242 (`FluidRegistry`) | 90+ real fluids, 11 ideal gases, 2 incompressibles, 2 mixtures; `getCoolPropName()` for CoolProp API |
| `include/coolsolve/fluids.h` lines 129–138 (`FluidRegistry` class) | `getFluid()`, `getAllFluids()`, `FluidType` enum |
| `src/server.cpp` lines 850–970 (solve thread) | Where `runner` is created and destroyed — inference data must be extracted here |
| `src/server.cpp` lines 1069–1098 (`/variables` endpoint) | Current endpoint returns name/value/isArray only — new `/variables/inferred` adds property/fluid/units |
| `gui/src/components/ArrayTable.tsx` | Array variable parsing (`var[i]` pattern) — same logic reused for array path overlay |

---

#### 4.10  Future Extensions (Not Part of This Phase)

- **Isolines**: isotherms, isobars, iso-quality lines on the dome.
  Each isoline requires a sweep along one property at a fixed other
  (e.g. T=const, sweep Q from 0 to 1 for isotherms inside the dome,
  or sweep P from P_min to P_max for isotherms outside).
- **Unit conversion**: display properties in user-friendly units
  (°C, kPa, kJ/kg, kJ/(kg·K)) instead of raw SI.
- **CoolProp.js**: run CoolProp in the browser via WebAssembly to
  compute missing state properties without backend calls.
- **Multiple fluids on one diagram**: overlay domes for comparison.
- **Diagram state persistence**: save/restore diagram configuration
  in the `.coolsolve` project bundle.
- **Interactive hover**: show property values on hover over the dome.

### Phase 5 — Advanced Features

| Task | Description |
|------|-------------|
| `.coolsolve` project format | Notebook-style bundling of code, sweeps, plots, cached results |
| Multi-model comparison | Solve two `.eescode` files side by side |
| Undo/redo | Full state snapshots |
| Collaborative editing | Online CRDT/OT (far future) |

---

## File Inventory

### Backend files (created/modified)

| File | Type |
|------|------|
| `include/coolsolve/server.h` | Created |
| `src/server.cpp` | Created |
| `cmake/embed_assets.cmake` | Created |
| `CMakeLists.txt` | Modified (added GUI option, cpp-httplib, asset embedding) |
| `main.cpp` | Modified (added `--gui` and `--no-browser` flags) |

### Frontend files (`gui/`)

| File | Purpose |
|------|---------|
| `gui/package.json` | Dependencies: react, @monaco-editor/react, zustand, lucide-react |
| `gui/vite.config.ts` | Vite config with `/api` proxy to backend |
| `gui/tsconfig.json` | TypeScript config |
| `gui/index.html` | HTML entry point |
| `gui/src/main.tsx` | React entry point |
| `gui/src/App.tsx` | Root layout component |
| `gui/src/App.css` | Full CSS with light/dark theme variables |
| `gui/src/api/types.ts` | TypeScript interfaces for API responses |
| `gui/src/api/client.ts` | Fetch-based API client |
| `gui/src/stores/modelStore.ts` | Zustand store: model state (code, results, console) |
| `gui/src/stores/uiStore.ts` | Zustand store: UI state (theme, tabs, panels) |
| `gui/src/languages/ees.ts` | Monaco language definition for EES |
| `gui/src/components/CodeEditor.tsx` | Monaco editor wrapper |
| `gui/src/components/VariableTable.tsx` | Variable table with filtering |
| `gui/src/components/Console.tsx` | Console log panel |
| `gui/src/components/Toolbar.tsx` | Toolbar with all actions + keyboard shortcuts |
| `gui/src/components/ConfigEditor.tsx` | Form-based config editor with 9 groups, 30+ fields |
| `gui/src/components/DebugViewer.tsx` | Debug output file list + content viewer |
| `gui/src/components/ArrayTable.tsx` | Spreadsheet view for array variables |

---

## How to Build & Run

### Development mode (hot-reload)

```bash
# Terminal 1: start backend
cd build && cmake .. && make -j$(nproc) coolsolve
./coolsolve --gui --no-browser

# Terminal 2: start frontend dev server
cd gui && npm install && npm run dev
# Open http://localhost:5173 — API calls are proxied to port 8550
```

### Production mode (single binary)

```bash
# Build frontend
cd gui && npm install && npm run build

# Rebuild backend (embeds gui/dist/ into binary)
cd build && cmake .. && make -j$(nproc) coolsolve

# Run
./coolsolve --gui
# Opens browser at http://localhost:8550
```
