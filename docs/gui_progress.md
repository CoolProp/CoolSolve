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
| `SplitPane.tsx` | ✅ Done | Reusable resizable split pane component (horizontal/vertical, min/max constraints, drag divider) |
| `App.tsx` + `App.css` | ✅ Done | Root layout with resizable split panes, toolbar, editor pane, right panel (tabbed), bottom panel, status bar, light/dark theme |
| Production build | ✅ Done | 185 kB JS + 6.3 kB CSS in `gui/dist/` |

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

### Phase 4 — Thermodynamic Diagrams

| Task | Description |
|------|-------------|
| CoolProp endpoints | `/coolprop/saturation-curve`, `/coolprop/props`, `/coolprop/fluids` |
| Diagram generation | Plotly T-s, P-h, h-s, T-h diagrams |
| State-point overlay | Plot solved array variables on diagrams |
| Iso-lines | Isotherms, isobars, iso-quality lines |
| Export | SVG/PNG export |

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
