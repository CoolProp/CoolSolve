# CoolSolve GUI — Architecture & Implementation Plan

## 1. Executive Summary

This document describes the design for a graphical user interface for CoolSolve.
The GUI must work **identically** in two deployment scenarios:

| Scenario | Backend | Frontend |
|----------|---------|----------|
| **Local** | CoolSolve binary starts an embedded HTTP server on `localhost` and opens the user's default browser | Same SPA served by the embedded server |
| **Online** | CoolSolve binary behind a reverse-proxy (nginx / Caddy) with session-based sandboxing | Same SPA served by the proxy or a CDN |

A **single web-based interface** is the recommended approach.  Building a
separate native GUI (Qt, GTK, Electron) would double the maintenance burden
for no meaningful usability gain — modern browsers provide file-system access,
drag-and-drop, offline support via Service Workers, and GPU-accelerated
rendering.  A single HTML/JS/CSS frontend also ensures perfect
cross-platform parity (Windows, macOS, Linux).

---

## 2. Technology Stack

### 2.1 Backend (C++ additions to CoolSolve)

| Component | Library | Rationale |
|-----------|---------|-----------|
| HTTP server | [cpp-httplib](https://github.com/yhirose/cpp-httplib) (header-only) | Zero-dependency, MIT-licensed, supports TLS, serves static files, handles multipart uploads. Already header-only like cpp-peglib. |
| JSON serialization | nlohmann/json (already a dependency) | REST request/response bodies |
| WebSocket (optional) | Built-in SSE (Server-Sent Events) via cpp-httplib | Real-time solve progress without adding a WebSocket library |
| Static asset embedding | CMake `file(READ ... HEX)` or [cmrc](https://github.com/vector-of-bool/cmrc) | Embed the compiled SPA into the binary so `coolsolve --gui` is a single executable with no external files |

### 2.2 Frontend (Single-Page Application)

| Component | Library | Rationale |
|-----------|---------|-----------|
| Framework | **React 18 + TypeScript** | Mature ecosystem, excellent component libraries, strong type safety |
| Code editor | **Monaco Editor** (`@monaco-editor/react`) | VS Code's editor engine — syntax highlighting, bracket matching, find/replace, minimap. A custom EES language definition provides keyword and function highlighting |
| Variable table | **AG Grid Community** | High-performance virtualised grid with in-cell editing, sorting, filtering, column pinning. Free tier is sufficient |
| Plotting | **Plotly.js** (`react-plotly.js`) | Publication-quality interactive plots (T-s, P-h, sensitivity sweeps). WebGL renderer for large datasets |
| UI chrome | **Radix UI + Tailwind CSS** | Accessible unstyled primitives + utility-first CSS. Light/dark theme with minimal bundle size |
| Build tool | **Vite** | Fast HMR during development; optimised production bundle |
| State management | **Zustand** | Lightweight, no boilerplate, works well with React concurrent features |
| Icons | **Lucide React** | Clean, consistent SVG icon set (MIT) |

### 2.3 Why Not Electron / Tauri?

- Electron ships a full Chromium (~150 MB) and a Node.js runtime.  CoolSolve
  already carries CoolProp (~30 MB); adding Electron triples the binary size
  and complicates updates.
- Tauri is lighter but still requires a Rust toolchain to build and introduces
  a second systems language into the project.
- The embedded-server approach delivers the same UX (the browser is always
  available) with **zero additional runtime dependencies**, and the same
  frontend works online without modification.

---

## 3. Modifications to CoolSolve Core

> **Design principle**: the GUI is an optional frontend.  All core
> functionality stays in the existing libraries (`parser.cpp`, `solver.cpp`,
> `evaluator.cpp`, `runner.h`).  The HTTP layer is a thin adapter that calls
> `CoolSolveRunner` exactly as `main.cpp` does today.

### 3.1 New Files

| File | Purpose |
|------|---------|
| `src/server.cpp` | HTTP server implementation: routes, SSE, static file serving |
| `include/coolsolve/server.h` | Public API: `startServer(port, options)` |
| `gui/` | Frontend source tree (React/Vite project, see §6) |

### 3.2 Changes to Existing Files

| File | Change | Overhead when `--gui` is not used |
|------|--------|----------------------------------|
| `main.cpp` | Add `--gui [port]` flag that calls `startServer()` instead of the CLI pipeline | **None** — the flag is simply not triggered |
| `CMakeLists.txt` | Add optional `COOLSOLVE_BUILD_GUI` option (default `ON`). When ON, builds `server.cpp`, links cpp-httplib, and embeds frontend assets. When OFF, the binary is identical to today | **None** when OFF |
| `include/coolsolve/runner.h` | Add `nlohmann::json toJSON() const` to `SolveResult` and `StructuralAnalysisResult` so the server can serialise results without custom code. These are thin wrappers around the existing `generateAnalysisJSON()` function | Negligible (unused code is stripped by the linker) |

### 3.3 No Overhead Guarantee

- `server.cpp` is compiled into a separate object file and only linked when
  `COOLSOLVE_BUILD_GUI=ON`.
- cpp-httplib is header-only but guarded behind `#ifdef COOLSOLVE_GUI`.
- Frontend assets are embedded as `const unsigned char[]` arrays; they occupy
  read-only data and are never touched unless the server starts.
- The CLI path (`main.cpp` without `--gui`) follows exactly the same code
  path as today with zero additional branches.

---

## 4. REST API Design

All endpoints are prefixed with `/api/v1/`.  Request and response bodies are JSON
unless noted otherwise.

### 4.1 Session Management

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/sessions` | Create a new session (returns `sessionId`). Online mode: required. Local mode: a default session is auto-created |
| `DELETE` | `/sessions/:id` | Destroy session and clean up temp files |

In **local mode**, a single implicit session is used and file paths map
directly to the local filesystem.  In **online mode**, each session gets an
isolated temp directory.

### 4.2 File Operations

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/files/eescode` | Return current .eescode content |
| `PUT` | `/files/eescode` | Save .eescode content (body: `{ "content": "..." }`) |
| `GET` | `/files/initials` | Return current .initials content |
| `PUT` | `/files/initials` | Save .initials content |
| `GET` | `/files/conf` | Return current coolsolve.conf as structured JSON |
| `PUT` | `/files/conf` | Save coolsolve.conf from structured JSON |
| `GET` | `/files/sol` | Return .sol content (after successful solve) |
| `POST` | `/files/upload` | Multipart upload of .eescode + .initials + .conf bundle |
| `GET` | `/files/bundle` | Download ZIP bundle (.eescode, .initials, .sol, .conf) |
| `POST` | `/files/open` | **Local only**: open a native file dialog (returns file path) |
| `POST` | `/files/save` | **Local only**: save to the original file path |
| `POST` | `/files/save-as` | **Local only**: save to a new file path |

### 4.3 Solver Operations

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/solve` | Start solving. Body: `{ "eescode": "...", "initials": "...", "conf": {...} }`. Returns immediately with a `taskId` |
| `GET` | `/solve/:taskId/status` | Poll solve status (pending / running / done / failed) |
| `GET` | `/solve/:taskId/result` | Get full solve result (variables, residuals, timing, blocks) |
| `GET` | `/solve/:taskId/stream` | **SSE endpoint**: streams progress events (block N/M, iteration, residual norm) in real time |
| `POST` | `/solve-debug` | Same as `/solve` but runs with `enableTracing=true` and generates debug output |
| `GET` | `/solve/:taskId/debug` | Retrieve debug output files (as JSON index + individual file contents) |
| `POST` | `/update-guesses` | Copy current .sol → .initials |

### 4.4 Analysis (read-only)

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/parse` | Parse .eescode and return variables, equations, errors (no solve). Useful for real-time editor feedback |
| `GET` | `/variables` | Return variable list with metadata (name, units, initial value, solved value, array index) |
| `GET` | `/latex` | Return LaTeX-formatted equations |
| `GET` | `/report` | Download compiled LaTeX report as PDF (requires `pdflatex` on server) or as `.tex` source |

### 4.5 Thermodynamic Helpers (future)

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/coolprop/props` | Evaluate CoolProp properties (for the diagram overlay). Body: `{ "fluid": "...", "input1": ..., "input2": ... }` |
| `POST` | `/coolprop/saturation-curve` | Return saturation dome points for a fluid |
| `GET` | `/coolprop/fluids` | List available CoolProp fluids |

### 4.6 CoolProp Warmup

On server start, `/api/v1/warmup` is called automatically to trigger
`warmupCoolProp()`.  The frontend shows a brief loading indicator until this
completes (~1–3 s on first launch).  Subsequent property evaluations are fast.

---

## 5. UI Layout

### 5.1 Overall Structure

```
┌──────────────────────────────────────────────────────────────────────────┐
│  TOOLBAR                                                                 │
│  [📂 Open] [💾 Save] [▶ Solve] [🐛 Debug] [↻ Update Guesses]           │
│  [{ } Comment] [" " Comment] [📥 Download Bundle] [📄 LaTeX Report]     │
├────────────────────────────┬─────────────────────────────────────────────┤
│                            │                                             │
│  CODE EDITOR (Monaco)      │  RIGHT PANEL (Tabs)                         │
│                            │  ┌─────────────────────────────────────┐    │
│  T_in = 25                 │  │ Variables │ Arrays │ Config │ Debug │    │
│  P = 101325                │  ├─────────────────────────────────────┤    │
│  h = enthalpy(Water,       │  │                                     │    │
│      T=T_in, P=P)          │  │  Variable table (AG Grid)           │    │
│  s = entropy(Water,        │  │  ┌──────┬──────┬───────┬──────┐     │    │
│      T=T_in, P=P)          │  │  │ Name │ Value│Initial│ Unit │     │    │
│                            │  │  ├──────┼──────┼───────┼──────┤     │    │
│                            │  │  │ T_in │ 25   │ 25    │ [°C] │     │    │
│                            │  │  │ P    │101325│101325 │ [Pa] │     │    │
│                            │  │  │ h    │104929│  —    │[J/kg]│     │    │
│                            │  │  │ s    │ 367  │  —    │[J/…] │     │    │
│                            │  │  └──────┴──────┴───────┴──────┘     │    │
│                            │  │                                     │    │
│                            │  └─────────────────────────────────────┘    │
├────────────────────────────┴─────────────────────────────────────────────┤
│  BOTTOM PANEL (collapsible)                                              │
│  ┌──────────────────────────────────────────────────────────────────────┐│
│  │ Console │ Sensitivity │ Thermo Diagrams │                            ││
│  ├──────────────────────────────────────────────────────────────────────┤│
│  │ Solve log, errors, timing, block progress                            ││
│  └──────────────────────────────────────────────────────────────────────┘│
└──────────────────────────────────────────────────────────────────────────┘
```

### 5.2 Detailed Panel Descriptions

#### A. Toolbar

A horizontal bar across the top with icon+label buttons, grouped logically:

| Group | Buttons | Behaviour |
|-------|---------|-----------|
| **File** | Open, Save, Save As, Download Bundle | Local: native file dialog via backend. Online: browser upload/download |
| **Solve** | Solve, Debug Solve, Stop | Solve triggers POST `/solve`, shows progress in console. Debug Solve runs with tracing. Stop cancels via task API |
| **Guesses** | Update Guesses | Copies solved values into the initials column (and optionally saves to .initials) |
| **Edit** | Toggle `{ }` Comment, Toggle `" "` Comment | Wraps/unwraps selected lines in the editor. Implemented as Monaco editor actions (no backend call) |
| **Export** | LaTeX Report, Download `.tex` | Calls `/report` endpoint. If `pdflatex` is available on the server, offers PDF; otherwise downloads `.tex` |
| **Settings** | Theme toggle (light/dark) | Stored in `localStorage` |

#### B. Code Editor (Left Pane)

- **Monaco Editor** with a custom `ees` language definition:
  - Keywords: `FUNCTION`, `PROCEDURE`, `END`, `CALL`, `IF`, `THEN`, `ELSE`,
    `$ifnot`, `$endif`, etc.
  - Built-in functions highlighted: `enthalpy`, `entropy`, `pressure`,
    `temperature`, `density`, `quality`, `cp`, `cv`, `specheat`, etc.
  - String/fluid literals: `'R134a'`, `'Water'`
  - Comment styles: `"..."`, `{...}`, `//...`
  - Bracket matching and auto-indentation
- Real-time parse feedback: as the user types, a debounced call to `/parse`
  returns errors, which are shown as red squiggly underlines (Monaco
  `setModelMarkers`).
- The editor occupies ~50% of the horizontal space by default; the divider
  is draggable.

#### C. Variables Tab (Right Pane — Default)

An **AG Grid** table with the following columns:

| Column | Editable | Description |
|--------|----------|-------------|
| **Name** | No | Variable name (sorted alphabetically, array vars grouped) |
| **Value** | No | Solved value (blank before solving, green highlight on convergence) |
| **Initial** | **Yes** | Editable initial guess. Synced with the in-memory .initials |
| **Units** | **Yes** | Editable unit annotation (e.g. `[kJ/kg]`) |
| **Block** | No | The SCC block index this variable belongs to |
| **Residual** | No | Per-equation residual at convergence (hidden by default, toggleable) |

- **Filtering**: a search box at the top filters by variable name.
- **Scalar vs Array**: scalar variables (`T_in`, `P`) and array variables
  (`h[1]`, `h[2]`) are shown together, with array vars visually grouped
  (indented or collapsible rows).
- **Update Guesses**: when the user clicks "Update Guesses" in the toolbar,
  the "Value" column is copied into the "Initial" column for all variables.

#### D. Arrays Tab (Right Pane)

For variables of the form `var[i]` (where `i` is an integer), a **spreadsheet
view** is displayed:

- Rows are the integer indices (1, 2, 3, …).
- Columns are the distinct base names (`T`, `P`, `h`, `s`, …).
- Each cell shows the solved value.
- A "Copy to Clipboard" button exports the table as TSV.
- A "Download CSV" button exports as CSV.
- If no array variables exist, the tab shows a placeholder message.

The array detection heuristic: any variable whose name matches the regex
`^[A-Za-z_]\w*\[\d+\]$` is treated as an array variable.

#### E. Config Tab (Right Pane)

A **form-based editor** for `coolsolve.conf` values, grouped into sections
matching the config file's structure:

- Main Iteration (maxIterations, tolerance, relativeTolerance, …)
- Line Search (lsAlpha, lsRho, …)
- Trust Region (trInitialRadius, trMaxRadius, …)
- Levenberg-Marquardt (lmInitialLambda, …)
- Partitioned Solver (partitionedMaxIterations, …)
- Tearing (enableTearing, tearingMinBlockSize, …)
- Pipeline (solverPipeline as a drag-and-drop reorderable list, pipelineMode
  as a radio button)
- Safety (timeoutSeconds)

Each field has its default value shown as placeholder text.  Changed values
are highlighted.  A "Reset to Defaults" button is available.  Changes are
sent to the backend on "Apply" or "Solve".

#### F. Debug Tab (Right Pane)

Displayed after a Debug Solve.  Shows the contents of the debug output folder
as a navigable tree:

- `report.md` rendered as formatted text
- `variables.md` rendered as a table
- `equations.md` rendered as formatted text with LaTeX math
- `incidence.md` rendered as a matrix view
- `evaluator.md` rendered as formatted text
- `tearing.md` rendered (if tearing was enabled)
- `analysis.json` viewable in a collapsible JSON tree
- `residuals.txt` in a code block
- `equations.tex` with LaTeX rendering (KaTeX)

A "Download All" button zips the debug folder for download.

#### G. Console (Bottom Panel — Default)

A scrollable log area that shows:

- Parse results (equation/variable count, errors)
- Solve progress (via SSE): block N of M, current iteration, residual norm
- Final result: SUCCESS / FAIL with timing breakdown
- CoolProp warnings
- Any errors with clickable line references that jump to the editor

#### H. Sensitivity Tab (Bottom Panel — Future)

**Phase 1 (skeleton):**
- A form to select a variable that is "imposed" (appears alone on one side of
  an equation, e.g. `T_in = 25`), define a range (`min`, `max`, `steps`), and
  choose output variables to track.
- A "Run Sweep" button that solves the system for each value in the range
  (calling `/solve` in batch) and collects results.
- A Plotly chart showing the output variables vs. the swept input variable.

**Phase 2 (future):**
- Multi-variable sweeps (2D heatmaps).
- Export sweep data as CSV.
- Save/load sweep configurations.

#### I. Thermodynamic Diagrams Tab (Bottom Panel — Future)

**Phase 1 (skeleton):**
- Fluid selector (dropdown with all CoolProp fluids).
- Diagram type selector: T-s, P-h, h-s, T-h.
- "Generate Diagram" button that calls `/coolprop/saturation-curve` and plots
  the dome using Plotly.

**Phase 2 (future):**
- **State-point overlay**: After solving, automatically plot the state points
  from the solution on the diagram.  Array variables `T[i]`, `P[i]`, `h[i]`,
  `s[i]` are connected as a cycle.
- **Iso-lines**: Add isotherms, isobars, iso-quality lines.
- **CoolProp Online integration**: Optionally link to [CoolProp
  Online](http://www.coolprop.org/fluid_properties/fluids.html) for
  interactive exploration.
- Export diagram as SVG/PNG.

---

## 6. Frontend Project Structure

```
gui/
├── package.json
├── tsconfig.json
├── vite.config.ts
├── index.html
├── public/
│   └── favicon.svg
├── src/
│   ├── main.tsx                    # App entry point
│   ├── App.tsx                     # Root layout (toolbar + panels)
│   ├── api/
│   │   ├── client.ts              # Axios/fetch wrapper for REST API
│   │   ├── types.ts               # TypeScript types matching API responses
│   │   └── sse.ts                 # SSE helper for solve progress
│   ├── stores/
│   │   ├── modelStore.ts          # Zustand: eescode, initials, conf, solve state
│   │   └── uiStore.ts            # Zustand: panel sizes, active tabs, theme
│   ├── components/
│   │   ├── Toolbar.tsx            # Top toolbar
│   │   ├── CodeEditor.tsx         # Monaco wrapper with EES language
│   │   ├── VariableTable.tsx      # AG Grid variable table
│   │   ├── ArrayTable.tsx         # Array variable spreadsheet
│   │   ├── ConfigEditor.tsx       # coolsolve.conf form editor
│   │   ├── DebugViewer.tsx        # Debug output viewer
│   │   ├── Console.tsx            # Bottom log panel
│   │   ├── SensitivityPanel.tsx   # Sensitivity sweep (future skeleton)
│   │   ├── ThermoDiagram.tsx      # Thermo diagram (future skeleton)
│   │   └── SplitPane.tsx          # Resizable split pane wrapper
│   ├── languages/
│   │   └── ees.ts                 # Monaco language definition for EES
│   └── utils/
│       ├── arrayVars.ts           # Detect and group array variables
│       └── formatting.ts          # Number formatting, unit display
└── tests/
    └── ...                         # Vitest unit tests
```

---

## 7. Project File Strategy

### 7.1 Local Mode — Direct Filesystem

When running locally, CoolSolve works directly with files on disk,
just as it does in CLI mode:

```
my_model/
├── my_model.eescode        # Source code
├── my_model.initials        # Initial guesses
├── my_model.sol             # Solution (generated)
├── coolsolve.conf           # Solver configuration
└── my_model_coolsolve/      # Debug output (generated, if debug solve used)
```

- **Open File** reads `.eescode` and auto-discovers `.initials`, `.sol`, and
  `coolsolve.conf` in the same directory (exact same logic as CLI).
- **Save** writes back to the same paths.
- **Save As** allows choosing a new directory/filename.

This is the simplest and most robust approach.  No proprietary project format
is needed — **the folder *is* the project**, exactly like the CLI workflow.

### 7.2 Online Mode — ZIP Bundles

- **Upload**: user selects one or more files (`.eescode` required, `.initials`
  and `coolsolve.conf` optional).  These are stored in a session-scoped temp
  directory.
- **Download**: the "Download Bundle" button generates a ZIP containing:
  - `.eescode` (current editor content)
  - `.initials` (current initials, including any user edits)
  - `.sol` (if a solve was successful)
  - `coolsolve.conf` (if modified from defaults)
- The ZIP file is named `<model_name>_coolsolve.zip`.

### 7.3 Future: Notebook-Style Projects

For the sensitivity analysis and plotting features (Phase 2+), a richer
project format may be desirable.  Inspired by Jupyter notebooks:

- A `.coolsolve` file (JSON) that bundles:
  - The `.eescode` source
  - Initial guesses
  - Solver configuration
  - A list of "analysis cells" (sensitivity sweeps, diagram configs, plot
    settings)
  - Cached results (last solution, sweep data)
- This is a **future extension** and is not needed for the initial
  implementation.  The decision on format should be deferred until the
  sensitivity and plotting features are designed in detail.

---

## 8. CoolProp Lifecycle

### 8.1 Warmup on Startup

When the server starts (or when the SPA loads), a single warmup call
initialises CoolProp:

```
Server start
  └── warmupCoolProp()          // ~1-3 seconds
      └── PropsSI("T", "P", 101325, "Q", 0, "Water")
```

The frontend shows a loading splash ("Initialising thermodynamic engine…")
until a `GET /api/v1/health` returns `{ "coolpropReady": true }`.

### 8.2 During Solving

CoolProp `AbstractState` instances are cached per fluid (already implemented
in CoolSolve via `CoolPropConfig::cacheEnabled`).  No additional changes are
needed.

### 8.3 Online Multi-User

In online mode, CoolProp's internal state is thread-safe when using separate
`AbstractState` instances per thread.  Each solve task runs in its own thread
with its own `CoolSolveRunner` instance.  No global state is shared.

---

## 9. Cross-Platform Considerations

| Concern | Solution |
|---------|----------|
| **Binary distribution** | Static linking of CoolProp + cpp-httplib + embedded frontend. Single binary (`coolsolve` / `coolsolve.exe`), no runtime dependencies. Ship pre-built binaries for Linux x86_64, macOS arm64/x86_64, Windows x86_64 |
| **Browser opening** | On `--gui`, use `xdg-open` (Linux), `open` (macOS), `start` (Windows) to launch the default browser |
| **File paths** | The backend normalises paths using `std::filesystem`. The REST API uses forward slashes. The frontend never manipulates raw paths |
| **File dialogs (local)** | Two options: (a) use the browser's `<input type="file">` and `showSaveFilePicker()` APIs (File System Access API), or (b) on platforms without browser FS API support, the backend spawns a native file dialog via a tiny helper (e.g., [tinyfiledialogs](https://sourceforge.net/projects/tinyfiledialogs/), header-only C) |
| **Port conflicts** | Default port 8550 (unlikely to conflict). Auto-increment if busy. Print actual URL to terminal |
| **Firewall (online)** | Online deployment behind reverse proxy handles TLS, rate limiting, and session isolation |

---

## 10. Security (Online Deployment)

| Threat | Mitigation |
|--------|------------|
| Arbitrary code execution | CoolSolve only parses `.eescode` — it does not execute arbitrary code. The PEG grammar is the sandbox |
| Resource exhaustion | Per-session timeout (`timeoutSeconds`), max concurrent solves, session expiry (30 min idle) |
| File system access | Online sessions use isolated temp directories. The `/files/open` and `/files/save` endpoints are **disabled** in online mode |
| Input size | Max upload size (e.g. 1 MB for `.eescode`) enforced by the HTTP server |

---

## 11. Implementation Phases

### Phase 1 — Core GUI (MVP)

**Goal**: A usable local GUI that replaces the CLI for interactive work.

| Step | Task | Effort |
|------|------|--------|
| 1.1 | Add cpp-httplib to CMake (FetchContent or vendored header) | 1 day |
| 1.2 | Implement REST API in `src/server.cpp`: session, file, solve, parse endpoints | 3–4 days |
| 1.3 | Implement SSE progress streaming from `CoolSolveRunner` (requires a callback hook in the solver loop) | 1–2 days |
| 1.4 | Add `--gui [port]` flag to `main.cpp` | 0.5 day |
| 1.5 | Scaffold React frontend (Vite + TypeScript) | 0.5 day |
| 1.6 | Implement Monaco editor with EES language definition | 2 days |
| 1.7 | Implement variable table (AG Grid) with editable initials/units | 2 days |
| 1.8 | Implement toolbar (Solve, Open, Save, Update Guesses, Comment toggles) | 1–2 days |
| 1.9 | Implement console panel with SSE progress | 1 day |
| 1.10 | Embed frontend assets into binary (CMake build step) | 1 day |
| 1.11 | Cross-platform testing (Linux, macOS, Windows) | 2 days |

**Estimated total: ~2–3 weeks**

### Phase 2 — Online Deployment & Polish

| Step | Task | Effort |
|------|------|--------|
| 2.1 | Session management and temp-directory sandboxing | 2 days |
| 2.2 | ZIP bundle upload/download | 1 day |
| 2.3 | Config editor tab (form-based coolsolve.conf editing) | 1–2 days |
| 2.4 | Debug output viewer tab | 1–2 days |
| 2.5 | Array variables tab (spreadsheet view) | 1 day |
| 2.6 | LaTeX report generation and download | 1–2 days |
| 2.7 | Deployment configuration (Docker image, nginx reverse proxy) | 1–2 days |
| 2.8 | Online demo deployment | 1 day |

**Estimated total: ~2 weeks**

### Phase 3 — Sensitivity Analysis

| Step | Task |
|------|------|
| 3.1 | Design sweep API: identify imposed variables, define sweep ranges |
| 3.2 | Batch solve endpoint (solves N times, returns table of results) |
| 3.3 | Frontend: sweep config form, result table, Plotly line charts |
| 3.4 | Export sweep data as CSV |

### Phase 4 — Thermodynamic Diagrams

| Step | Task |
|------|------|
| 4.1 | CoolProp saturation-curve / isoline endpoint |
| 4.2 | Frontend: Plotly T-s, P-h, h-s diagram generation |
| 4.3 | State-point overlay from solved array variables |
| 4.4 | Export diagrams as SVG/PNG |

### Phase 5 — Advanced Features

| Step | Task |
|------|------|
| 5.1 | `.coolsolve` notebook-style project format (bundling sweeps, plots, code) |
| 5.2 | Multi-model comparison (solve two `.eescode` files side by side) |
| 5.3 | Undo/redo with full state snapshots |
| 5.4 | Collaborative editing (online, via CRDT or OT — far future) |

---

## 12. Solver Callback Hook (Core Change)

To enable real-time progress reporting via SSE, the solver needs a callback
mechanism.  This is the **only meaningful change** to core solving code:

```cpp
// In solver.h — add to SolverOptions:
using ProgressCallback = std::function<void(
    int blockIndex,       // Current block being solved
    int totalBlocks,      // Total number of blocks
    int iteration,        // Current iteration within block
    double residualNorm   // Current ||F||
)>;

std::function<ProgressCallback> progressCallback = nullptr;
```

The callback is invoked at the start of each Newton iteration.  When
`progressCallback` is `nullptr` (CLI mode), there is zero overhead — the
compiler eliminates the null check.  The server sets this callback to
push events onto an SSE channel.

---

## 13. Monaco EES Language Definition

A custom Monaco language contribution provides:

```typescript
// gui/src/languages/ees.ts

export const eesLanguage: monaco.languages.IMonarchLanguage = {
  keywords: [
    'FUNCTION', 'PROCEDURE', 'END', 'CALL',
    'IF', 'THEN', 'ELSE', 'ENDIF',
    'REPEAT', 'UNTIL', 'GOTO',
    '$ifnot', '$endif', '$common', '$include',
  ],
  
  builtinFunctions: [
    'enthalpy', 'entropy', 'density', 'pressure', 'temperature',
    'quality', 'cp', 'cv', 'specheat', 'conductivity', 'viscosity',
    'volume', 'soundspeed', 'surfacetension',
    'sin', 'cos', 'tan', 'exp', 'log', 'log10', 'sqrt', 'abs',
    'min', 'max', 'ln',
  ],
  
  tokenizer: {
    root: [
      [/"[^"]*"/, 'comment'],           // "double-quote comments"
      [/\{[^}]*\}/, 'comment'],         // {brace comments}
      [/\/\/.*$/, 'comment'],           // // line comments
      [/'[^']*'/, 'string'],            // 'string literals' (fluid names)
      [/\b\d+(\.\d+)?([eE][+-]?\d+)?\b/, 'number'],
      [/[A-Za-z_]\w*/, {
        cases: {
          '@keywords': 'keyword',
          '@builtinFunctions': 'support.function',
          '@default': 'variable',
        }
      }],
    ],
  },
};
```

This gives users syntax highlighting comparable to the EES IDE.

---

## 14. Comment Toggle Implementation

Two toolbar buttons for commenting, implemented entirely in the Monaco editor
(no backend call):

### Invisible Comment (`{ }`)

- **Toggle on**: wrap selected lines in `{ ... }` (the parser ignores this).
  If the selection is already wrapped, remove the braces.
- Implementation: Monaco `executeEdits` API.

### Visible Comment (`" "`)

- **Toggle on**: wrap selected text in `"..."` (displayed as a remark in EES).
  If already wrapped, remove the quotes.
- This is EES's "documentation comment" — it's shown in the output but
  the equation inside is ignored by the solver.

Both operations register as Monaco `editor.action` contributions so they
also get keyboard shortcuts (e.g., `Ctrl+/` for invisible, `Ctrl+Shift+/`
for visible).

---

## 15. Deployment Architecture (Online)

```
                    ┌─────────────┐
                    │   Browser   │
                    │   (SPA)     │
                    └──────┬──────┘
                           │ HTTPS
                    ┌──────▼──────┐
                    │   nginx /   │
                    │   Caddy     │
                    │  (reverse   │
                    │   proxy)    │
                    └──────┬──────┘
                           │ HTTP (localhost:8550)
                    ┌──────▼──────┐
                    │  CoolSolve  │
                    │  --gui      │
                    │  (embedded  │
                    │   server)   │
                    └─────────────┘
```

A **Docker image** packages the CoolSolve binary with the embedded frontend:

```dockerfile
FROM ubuntu:24.04 AS runtime
COPY coolsolve /usr/local/bin/coolsolve
EXPOSE 8550
ENTRYPOINT ["coolsolve", "--gui", "--online", "8550"]
```

The `--online` flag enables session management, disables filesystem
endpoints, and sets appropriate resource limits.

---

## 16. Additional Ideas

### 16.1 Auto-Save and Recovery

- In local mode, auto-save the editor content every 30 seconds to a
  `.coolsolve.autosave` file.  On next startup, if this file is newer than
  the `.eescode`, offer to recover.
- In online mode, use `localStorage` to persist the editor content across
  browser refreshes.

### 16.2 Variable Quick-Jump

- Clicking a variable name in the table highlights it in the editor (and
  vice versa).
- `Ctrl+Click` on a variable in the editor opens a popup showing its
  current value, units, and the equation that determines it.

### 16.3 Equation Dependency Graph

- A visual graph (using D3.js or Mermaid) showing the block structure:
  which blocks depend on which, with block sizes and solve status
  (green = converged, red = failed).
- Accessible from the Debug tab or as a separate panel.

### 16.4 Smart Initial Guesses

- When a variable has no initial guess, highlight it in orange in the
  variable table.
- For thermodynamic properties, suggest physically reasonable ranges
  based on the fluid and typical operating conditions.

### 16.5 Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| `Ctrl+Enter` | Solve |
| `Ctrl+Shift+Enter` | Debug Solve |
| `Ctrl+S` | Save |
| `Ctrl+O` | Open |
| `Ctrl+/` | Toggle `{ }` comment |
| `Ctrl+Shift+/` | Toggle `" "` comment |
| `Ctrl+G` | Update guesses |
| `Escape` | Stop solve |

### 16.6 Responsive Design

- On narrow screens (tablets), the right panel moves below the editor
  (vertical split instead of horizontal).
- The bottom panel collapses to a thin status bar showing only the last
  log line, expandable on click.

### 16.7 Export to Python / Julia

- Generate a Python or Julia script that reproduces the model using
  CoolProp's native Python/Julia bindings.  This would allow users to
  extend their models in a general-purpose language.

### 16.8 Model Library / Examples

- The GUI ships with the `examples/` folder accessible from a "File →
  Open Example" menu.  Each example has a short description.
- Online mode: examples are pre-loaded and selectable from a dropdown.

### 16.9 Diff View for Guesses vs Solution

- After solving, a side-by-side diff showing how each variable moved
  from its initial guess to the converged value.  Large jumps are
  highlighted, helping users identify poor guesses.

### 16.10 Unit Conversion Helper

- A small utility panel where users can convert between common
  engineering units (°C↔K, bar↔Pa, kJ↔J, etc.).
- Could also flag variables whose solved values seem inconsistent with
  their declared units.

---

## 17. Summary of Core Code Changes

| Change | Files Affected | Impact on CLI |
|--------|---------------|---------------|
| Add `progressCallback` to `SolverOptions` | `solver.h`, `solver.cpp` | None (null check eliminated by compiler) |
| Add `toJSON()` to result structs | `runner.h`, `runner.cpp` | None (unused methods stripped by linker) |
| New `server.cpp` + `server.h` | New files only | None (not linked when `COOLSOLVE_BUILD_GUI=OFF`) |
| Add `--gui` flag | `main.cpp` | None (flag not triggered) |
| Add `COOLSOLVE_BUILD_GUI` CMake option | `CMakeLists.txt` | None when `OFF` (default can be `ON` since overhead is only build-time) |
| Add cpp-httplib dependency | `CMakeLists.txt` | None (header-only, only `#include`d in `server.cpp`) |

**Total lines of core C++ changes: ~50 lines across `solver.h` and `runner.h`.**  
Everything else is additive (new files).

---

## 18. Testing Strategy

The GUI must be testable **end-to-end by a coding agent** (e.g. GitHub
Copilot, Cursor, Aider) that can compile, launch the server, hit the API,
and verify results — all without human interaction.  The strategy is layered:

### 18.1 Test Pyramid

```
                    ┌─────────────────────┐
                    │  E2E / Playwright   │  ← browser automation (few)
                    ├─────────────────────┤
                    │  API Integration    │  ← curl / httpie / Python (many)
                    ├─────────────────────┤
                    │  C++ Unit Tests     │  ← Catch2 (existing + new)
                    ├─────────────────────┤
                    │  Frontend Unit      │  ← Vitest + React Testing Lib
                    └─────────────────────┘
```

### 18.2 C++ Backend Tests (Catch2)

Extend the existing Catch2 suite with server-specific tests:

| Test file | Tag | What it covers |
|-----------|-----|----------------|
| `test_server.cpp` | `[server]` | Starts the embedded HTTP server on a random port, makes HTTP requests using cpp-httplib's client, verifies JSON responses |
| `test_server_solve.cpp` | `[server-solve]` | Uploads an `.eescode` via `POST /api/v1/solve`, polls `/status`, retrieves `/result`, and asserts convergence + variable values |
| `test_server_files.cpp` | `[server-files]` | Tests file upload/download, bundle ZIP creation, initials round-trip |
| `test_server_sse.cpp` | `[server-sse]` | Connects to the SSE stream during a solve and asserts that progress events are emitted |

**Key design: random port allocation.**  Each test creates a server on port 0
(OS-assigned), avoiding conflicts when tests run in parallel or when another
server is already running.

```cpp
// Example: test_server_solve.cpp
TEST_CASE("Solve via REST API", "[server-solve]") {
    // Start server on random port
    int port = startTestServer();  // returns actual port
    httplib::Client cli("localhost", port);

    // Upload and solve orc_simple
    auto eescode = readFile("../examples/orc_simple.eescode");
    auto initials = readFile("../examples/orc_simple.initialsOK");
    json body = {{"eescode", eescode}, {"initials", initials}};

    auto res = cli.Post("/api/v1/solve", body.dump(), "application/json");
    REQUIRE(res->status == 202);
    auto taskId = json::parse(res->body)["taskId"].get<std::string>();

    // Poll until done (max 30 s)
    json result;
    for (int i = 0; i < 300; ++i) {
        auto status = cli.Get("/api/v1/solve/" + taskId + "/status");
        auto j = json::parse(status->body);
        if (j["status"] == "done" || j["status"] == "failed") {
            result = j;
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    REQUIRE(result["status"] == "done");

    // Retrieve variables
    auto vars = cli.Get("/api/v1/solve/" + taskId + "/result");
    auto jvars = json::parse(vars->body);
    REQUIRE(jvars["success"] == true);

    stopTestServer(port);
}
```

### 18.3 API Integration Tests (Shell / Python)

A script that a **coding agent can execute in a terminal** to verify the full
stack without needing a browser.  Lives in `tests/test_api.sh` (zero
dependencies beyond `curl` and `jq`):

```bash
#!/usr/bin/env bash
# tests/test_api.sh — Headless API integration tests
# Usage: ./test_api.sh [port]
# If port is omitted, starts coolsolve --gui on a random port.

set -euo pipefail

PORT=${1:-0}
PASS=0; FAIL=0; TOTAL=0

# ── Helpers ──────────────────────────────────────────────────────────
assert_eq()  { ((TOTAL++)); if [[ "$1" == "$2" ]]; then ((PASS++)); else ((FAIL++)); echo "FAIL: $3 (expected '$2', got '$1')"; fi }
assert_ok()  { ((TOTAL++)); if [[ "$1" -ge 200 && "$1" -lt 300 ]]; then ((PASS++)); else ((FAIL++)); echo "FAIL: $3 (HTTP $1)"; fi }
cleanup()    { [[ -n "${SERVER_PID:-}" ]] && kill "$SERVER_PID" 2>/dev/null; }
trap cleanup EXIT

# ── Start server ─────────────────────────────────────────────────────
if [[ "$PORT" == "0" ]]; then
    PORT=8551  # fixed fallback; could parse from stdout
    ../build/coolsolve --gui "$PORT" &
    SERVER_PID=$!
    sleep 3  # wait for server + CoolProp warmup
fi
BASE="http://localhost:${PORT}/api/v1"

# ── Health check ─────────────────────────────────────────────────────
HTTP=$(curl -s -o /dev/null -w '%{http_code}' "$BASE/health")
assert_ok "$HTTP" 200 "GET /health"

# ── Parse-only endpoint ──────────────────────────────────────────────
HTTP=$(curl -s -o /dev/null -w '%{http_code}' -X POST "$BASE/parse" \
    -H 'Content-Type: application/json' \
    -d '{"eescode":"T=25\nP=101325"}')
assert_ok "$HTTP" 200 "POST /parse"

# ── Full solve (orc_simple) ──────────────────────────────────────────
EESCODE=$(cat ../examples/orc_simple.eescode)
INITIALS=$(cat ../examples/orc_simple.initialsOK)
RESPONSE=$(curl -s -X POST "$BASE/solve" \
    -H 'Content-Type: application/json' \
    -d "$(jq -n --arg e "$EESCODE" --arg i "$INITIALS" '{eescode:$e, initials:$i}')")
TASK_ID=$(echo "$RESPONSE" | jq -r '.taskId')

# Poll until done
for i in $(seq 1 60); do
    STATUS=$(curl -s "$BASE/solve/${TASK_ID}/status" | jq -r '.status')
    [[ "$STATUS" == "done" || "$STATUS" == "failed" ]] && break
    sleep 0.5
done
assert_eq "$STATUS" "done" "Solve orc_simple"

# Verify variables
RESULT=$(curl -s "$BASE/solve/${TASK_ID}/result")
SUCCESS=$(echo "$RESULT" | jq -r '.success')
assert_eq "$SUCCESS" "true" "orc_simple converged"

# ── Download bundle ──────────────────────────────────────────────────
HTTP=$(curl -s -o /dev/null -w '%{http_code}' "$BASE/files/bundle")
assert_ok "$HTTP" 200 "GET /files/bundle"

# ── Report ───────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════"
echo "  Results: $PASS / $TOTAL passed"
[[ "$FAIL" -gt 0 ]] && echo "  FAILURES: $FAIL" && exit 1
echo "  All tests passed."
echo "═══════════════════════════════════"
```

A **Python variant** (`tests/test_api.py`) using `requests` + `pytest` is
also provided for richer assertions and easier extension:

```python
# tests/test_api.py — pytest-based API tests
# Usage:  pytest tests/test_api.py --port 8550
#    or:  pytest tests/test_api.py              (auto-starts server)

import subprocess, time, requests, pytest, os, json

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")

@pytest.fixture(scope="session")
def server(request):
    port = request.config.getoption("--port", default=None)
    if port:
        yield f"http://localhost:{port}/api/v1"
        return
    port = 8552
    proc = subprocess.Popen(
        ["../build/coolsolve", "--gui", str(port)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    time.sleep(4)  # CoolProp warmup
    yield f"http://localhost:{port}/api/v1"
    proc.terminate()

def solve_example(server, name):
    """Solve an example and return the result JSON."""
    eescode = open(os.path.join(EXAMPLES_DIR, f"{name}.eescode")).read()
    initials_candidates = [f"{name}.initials", f"{name}.initialsOK", f"{name}.initialsl"]
    initials = ""
    for ic in initials_candidates:
        path = os.path.join(EXAMPLES_DIR, ic)
        if os.path.exists(path):
            initials = open(path).read()
            break
    conf_path = os.path.join(EXAMPLES_DIR, "coolsolve.conf")
    conf = open(conf_path).read() if os.path.exists(conf_path) else ""
    r = requests.post(f"{server}/solve",
                      json={"eescode": eescode, "initials": initials, "conf": conf})
    assert r.status_code == 202
    task_id = r.json()["taskId"]
    for _ in range(120):
        s = requests.get(f"{server}/solve/{task_id}/status").json()
        if s["status"] in ("done", "failed"):
            break
        time.sleep(0.5)
    result = requests.get(f"{server}/solve/{task_id}/result").json()
    return result

class TestHealth:
    def test_health(self, server):
        r = requests.get(f"{server}/health")
        assert r.status_code == 200
        assert r.json()["coolpropReady"] is True

class TestParse:
    def test_simple(self, server):
        r = requests.post(f"{server}/parse",
                          json={"eescode": "T=25\nP=101325"})
        assert r.status_code == 200
        data = r.json()
        assert data["variableCount"] == 2
        assert data["equationCount"] == 2
        assert len(data.get("errors", [])) == 0

class TestSolve:
    @pytest.mark.parametrize("example", [
        "orc_simple", "rankine1", "refrigeration1",
        "exchangers1", "pressuredrop",
    ])
    def test_example_solves(self, server, example):
        result = solve_example(server, example)
        assert result["success"], f"{example} failed: {result.get('errorMessage')}"

    def test_variables_returned(self, server):
        result = solve_example(server, "rankine1")
        assert "variables" in result
        assert len(result["variables"]) > 0

class TestFiles:
    def test_bundle_download(self, server):
        # First solve something
        solve_example(server, "orc_simple")
        r = requests.get(f"{server}/files/bundle")
        assert r.status_code == 200
        assert r.headers["Content-Type"] == "application/zip"

    def test_update_guesses(self, server):
        solve_example(server, "orc_simple")
        r = requests.post(f"{server}/update-guesses")
        assert r.status_code == 200
```

### 18.4 Example-Based Regression Sweep

A dedicated script iterates over **all** `.eescode` files in `examples/` via
the REST API, replicating what `test_examples.cpp` does for the CLI but
exercising the entire GUI stack:

```bash
#!/usr/bin/env bash
# tests/test_all_examples_api.sh — Solve every example via the REST API
set -euo pipefail
PORT=${1:-8551}
BASE="http://localhost:${PORT}/api/v1"
PASS=0; FAIL=0

for f in ../examples/*.eescode; do
    NAME=$(basename "$f" .eescode)
    EESCODE=$(cat "$f")

    # Find initials
    INITIALS=""
    for ext in .initials .initialsOK .initialsl; do
        IPATH="../examples/${NAME}${ext}"
        [[ -f "$IPATH" ]] && INITIALS=$(cat "$IPATH") && break
    done

    # Solve
    RESPONSE=$(curl -s -X POST "$BASE/solve" \
        -H 'Content-Type: application/json' \
        -d "$(jq -n --arg e "$EESCODE" --arg i "$INITIALS" '{eescode:$e,initials:$i}')")
    TASK_ID=$(echo "$RESPONSE" | jq -r '.taskId')

    # Poll (60 s timeout)
    STATUS="pending"
    for i in $(seq 1 120); do
        STATUS=$(curl -s "$BASE/solve/${TASK_ID}/status" | jq -r '.status')
        [[ "$STATUS" == "done" || "$STATUS" == "failed" ]] && break
        sleep 0.5
    done

    if [[ "$STATUS" == "done" ]]; then
        echo "  OK   $NAME"
        ((PASS++))
    else
        echo "  FAIL $NAME (status: $STATUS)"
        ((FAIL++))
    fi
done

echo ""
echo "Results: $PASS passed, $FAIL failed ($(( PASS + FAIL )) total)"
[[ "$FAIL" -gt 0 ]] && exit 1
exit 0
```

### 18.5 Frontend Tests (Vitest)

Unit tests for React components and utility functions:

| Test | What it covers |
|------|----------------|
| `VariableTable.test.tsx` | Renders variable grid, editable initials update store |
| `ArrayTable.test.tsx` | Array variable grouping heuristic (regex matching) |
| `CodeEditor.test.tsx` | Monaco loads, EES language registered, comment toggle |
| `ConfigEditor.test.tsx` | Form renders all config keys, validates numeric inputs |
| `api/client.test.ts` | Mock fetch, request/response serialization |
| `stores/modelStore.test.ts` | State transitions: idle → solving → solved → error |
| `utils/arrayVars.test.ts` | `T[1]` → grouped, `T_in` → scalar, `T[abc]` → scalar |

Run with:
```bash
cd gui && npx vitest run
```

### 18.6 E2E Browser Tests (Playwright)

For critical user workflows that must be tested in a real browser:

| Test | Scenario |
|------|----------|
| `solve-flow.spec.ts` | Open editor → type code → click Solve → verify variables appear in table |
| `file-roundtrip.spec.ts` | Upload `.eescode` → solve → download bundle → re-upload → verify same result |
| `comment-toggle.spec.ts` | Select lines → click `{ }` comment → verify braces added → click again → braces removed |
| `config-edit.spec.ts` | Open Config tab → change `maxIterations` → Solve → verify config was used |

Run with:
```bash
cd gui && npx playwright test
```

Playwright tests are **optional** for a coding agent (they require a headed
or headless browser).   The API tests (§18.3) cover the same backend logic
without a browser dependency and are the primary agent-driven test suite.

### 18.7 Coding-Agent Workflow (Step by Step)

A coding agent implementing a change to the GUI should follow this exact
sequence.  Every step is a terminal command with a deterministic exit code:

```bash
# 1. Build the C++ backend (including server)
cd build && cmake -DCOOLSOLVE_BUILD_GUI=ON .. && make -j$(nproc)
# Exit 0 = build OK, non-zero = compile error → fix and retry

# 2. Run C++ unit tests (fast, no server needed)
./coolsolve_tests "[server]"
# Exit 0 = all pass

# 3. Build the frontend
cd ../gui && npm ci && npm run build
# Exit 0 = build OK

# 4. Run frontend unit tests
npx vitest run
# Exit 0 = all pass

# 5. Re-embed frontend assets and rebuild binary
cd ../build && make -j$(nproc)
# (The CMake build step copies gui/dist/ into the binary)

# 6. Start the server in the background
./coolsolve --gui 8551 &
SERVER_PID=$!
sleep 4  # CoolProp warmup

# 7. Run API integration tests
cd ../tests && bash test_api.sh 8551
# Exit 0 = all pass

# 8. (Optional) Run full example sweep
bash test_all_examples_api.sh 8551
# Exit 0 = all examples solve

# 9. Stop the server
kill $SERVER_PID
```

**Key properties for agent autonomy:**
- Every step has a binary pass/fail exit code
- No interactive prompts or browser windows needed for steps 1–8
- The `examples/` folder provides ready-made test inputs with known-good solutions
- Errors are reported to stderr with enough context to diagnose the issue
- Steps 1–4 can be run independently (frontend and backend are decoupled)
- Steps 6–8 test the integrated stack end-to-end without a browser

### 18.8 CI Pipeline

The same steps integrate into GitHub Actions / GitLab CI:

```yaml
# .github/workflows/gui.yml
name: GUI Tests
on: [push, pull_request]
jobs:
  backend:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Build
        run: |
          mkdir -p build && cd build
          cmake -DCOOLSOLVE_BUILD_GUI=ON ..
          make -j$(nproc)
      - name: Unit tests
        run: cd build && ./coolsolve_tests "[server]"
      - name: API integration
        run: |
          cd build && ./coolsolve --gui 8551 &
          sleep 4
          cd ../tests && bash test_api.sh 8551
          kill %1

  frontend:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with: { node-version: '20' }
      - name: Install & test
        run: cd gui && npm ci && npx vitest run
      - name: Build
        run: cd gui && npm run build

  examples:
    runs-on: ubuntu-latest
    needs: backend
    steps:
      - uses: actions/checkout@v4
      - name: Build
        run: mkdir -p build && cd build && cmake -DCOOLSOLVE_BUILD_GUI=ON .. && make -j$(nproc)
      - name: Sweep all examples
        timeout-minutes: 15
        run: |
          cd build && ./coolsolve --gui 8551 &
          sleep 4
          cd ../tests && bash test_all_examples_api.sh 8551
          kill %1
```

### 18.9 Test Data & Fixtures

The `examples/` directory already contains all the test data needed:

| File pattern | Count | Purpose |
|-------------|-------|---------|
| `*.eescode` | ~25 | Model source code |
| `*.initials` / `*.initialsOK` | ~20 | Initial guesses |
| `*.sol` | ~8 | Known-good solutions (for regression checks) |
| `*.residuals` | ~20 | EES reference outputs (for structural comparison) |
| `coolsolve.conf` | 1 | Default solver configuration |

The `.sol` files serve as **golden references**: after solving via the API,
the test can compare each variable against the `.sol` value with a 1%
tolerance (matching the existing `test_examples.cpp` logic).

No synthetic test data needs to be created — the examples provide
comprehensive coverage of:
- Simple systems (2–4 equations)
- Medium systems (10–30 equations with CoolProp calls)
- Complex systems (50+ equations, multiple algebraic loops, tearing)
- Edge cases (array variables, procedures, string variables)

---

## 19. Open Questions

1. **PDF generation**: Should LaTeX report compilation (`pdflatex`) be a
   server-side feature, or should the frontend generate PDFs client-side
   (e.g., via a WASM LaTeX engine like SwiftLaTeX)?  Server-side is simpler
   but requires `pdflatex` to be installed.

2. **Authentication (online)**: Is the online deployment public, or does it
   require login?  For a public demo, sessions with expiry are sufficient.
   For a multi-user service, OAuth (GitHub/Google) could be added later.

3. **WASM fallback**: Compiling CoolSolve + CoolProp to WebAssembly would
   allow the online version to run entirely in the browser (no server
   needed).  This is a significant undertaking but would eliminate server
   costs.  Worth exploring as a long-term option.

4. **Mobile support**: Is there a need for the GUI to work on phones?  The
   current design targets desktop browsers (≥1024px width).  Minimal
   mobile support could be added via responsive breakpoints, but the code
   editor experience on phones is inherently limited.

5. **Test coverage threshold**: Should a minimum coverage percentage be enforced
   in CI (e.g., 80% for frontend, 90% for API routes)?  This can be added
   once the test suites stabilize.
