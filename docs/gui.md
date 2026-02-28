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

### 4.3 File Operations

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/files/eescode` | Get current EES source code |
| `PUT` | `/api/v1/files/eescode` | Set EES source code; body: `{"content":"..."}` |
| `GET` | `/api/v1/files/initials` | Get initial guesses text |
| `PUT` | `/api/v1/files/initials` | Set initial guesses text |
| `GET` | `/api/v1/files/sol` | Get last solution file content |
| `GET` | `/api/v1/files/conf` | Get solver configuration (as text) |
| `PUT` | `/api/v1/files/conf` | Set solver configuration |
| `POST` | `/api/v1/files/open` | Open a `.eescode` file from disk; auto-discovers companion `.initials` and `coolsolve.conf` |
| `GET` | `/api/v1/files/bundle` | Download a ZIP bundle (eescode + initials + sol + conf + debug files) |
| `POST` | `/api/v1/files/upload` | Upload a multipart `.zip` bundle or individual files |
| `GET` | `/api/v1/examples` | List built-in example files |

The file workflow is **ZIP-centric**:
- Open loads files from disk (local mode) or from an uploaded ZIP.
- Save produces a ZIP bundle containing all companion files plus any debug output.
- Examples are read-only on the server and opened via `POST /api/v1/files/open`.

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
| `POST` | `/api/v1/parse` | Parse EES code (no solve); returns equation/variable counts and parse errors |
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

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/coolprop/fluids` | List all registered fluids with `name`, `type`, `hasDome`, and `modelFluids` from the current session |
| `POST` | `/api/v1/coolprop/saturation` | Compute saturation dome for a fluid; body: `{"fluid":"R134a","nPoints":200}`; returns liquid+vapor arrays for T, P, H, S, D and critical/triple-point data |

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
│  STATUS BAR / CONSOLE (collapsible bottom panel)                         │
│  Solve log, errors, block progress, timing                               │
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
| Open | `Ctrl+O` | Load `.eescode` from disk + companion files |
| Save | `Ctrl+S` | Download ZIP bundle |
| Solve | `Ctrl+Enter` | Run async solve, stream progress via SSE |
| Debug Solve | `Ctrl+Shift+Enter` | Solve with tracing enabled, generate debug output |
| Stop | `Escape` | Cancel running solve |
| Update Guesses | `Ctrl+G` | Copy solved values into initial guesses |
| Toggle `{ }` Comment | `Ctrl+/` | Wrap/unwrap selection in brace comments |
| Toggle `" "` Comment | `Ctrl+Shift+/` | Wrap/unwrap selection in double-quote comments |
| Theme | — | Toggle light/dark theme (persisted in `localStorage`) |

### 6.3 Code Editor (Left Pane) — `CodeEditor.tsx`

Monaco Editor with a custom `ees` language definition (`languages/ees.ts`):
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

A form-based editor for `coolsolve.conf` values, organised into 9 collapsible
groups covering all `SolverOptions` keys (~30+ fields):

| Group | Key fields |
|-------|-----------|
| Main Iteration | `maxIterations`, `tolerance`, `relativeTolerance`, `verbose` |
| Line Search | `lsAlpha`, `lsRho`, `lsMaxIterations` |
| Trust Region | `trInitialRadius`, `trMaxRadius`, `trEta` |
| Levenberg-Marquardt | `lmInitialLambda`, `lmLambdaFactor` |
| Partitioned Solver | `partitionedMaxIterations`, `partitionedTolerance` |
| Tearing | `enableTearing`, `tearingMinBlockSize`, `tearingMethod` |
| Pipeline | `solverPipeline` (ordered list), `pipelineMode` (radio) |
| Safety | `timeoutSeconds` |
| Output | `writeSolFile`, `writeResiduals` |

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

### 6.8 Diagrams Tab — `ThermoDiagram.tsx`

Interactive thermodynamic property diagrams powered by Plotly.  The tab is
located in the right panel alongside Variables, Arrays, Config, and Debug.

#### Controls

- **Fluid dropdown**: populated from `GET /api/v1/coolprop/fluids`.  Model
  fluids (detected via variable inference) appear first.  Only fluids with
  `hasDome: true` (i.e. CoolProp `Real` type) are selectable; ideal gases and
  humid air are shown greyed out.
- **Diagram type**: T-s, P-h, h-s, T-h (P-h uses logarithmic y-axis)
- **Generate button**: calls `POST /api/v1/coolprop/saturation`; shows spinner
  during computation (~100–300 ms for 200 points)

#### Saturation Dome

The backend generates the dome using N logarithmically spaced temperature
points between the triple point and critical point:
1. For each T: compute Q=0 (liquid) and Q=1 (vapor) values of P, H, S, D
   via `CoolProp::PropsSI()`
2. Concatenate liquid branch + critical point + vapor branch
3. Return arrays for T, P, H, S, D plus critical and triple-point data (all SI units)

Each dome generation requires ~200 × 8 = 1,600 `PropsSI` calls, taking
approximately **0.1–0.3 s** in practice — acceptable for a user-initiated action.

#### State Point Overlay

After solving, checking "Overlay solved states" calls
`GET /api/v1/variables/inferred` and:
1. Filters variables matching the selected fluid (`inferredFluid == selectedFluid`)
2. Groups by index `i` for array variables, or by name-pattern matching for scalars
3. Plots each group as a labeled scatter marker on the diagram (where both
   axis properties are directly available from inference)

#### Array Path Overlay

Checking "Overlay array path" lets the user select X and Y base variable names.
The component:
1. Extracts paired values `(X[i], Y[i])` for all matching indices
2. Sorts by index order
3. Plots as a connected scatter+line trace with index labels
4. Optional "Close cycle" checkbox connects the last point back to the first

#### Export

PNG and SVG export via Plotly's built-in `Plotly.downloadImage()`.

### 6.9 Console (Bottom Panel) — `Console.tsx`

Auto-scrolling log panel with color-coded lines:
- Parse results (equation/variable count, errors)
- Solve progress via SSE (block N/M, iteration, residual norm)
- Final result: SUCCESS / FAIL with timing breakdown
- CoolProp warnings and errors

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

- **Open** (`POST /api/v1/files/open`): reads `.eescode` and auto-discovers
  `.initials`, `.sol`, and `coolsolve.conf` in the same directory.
- **Save** (`GET /api/v1/files/bundle`): downloads a ZIP bundle named
  `<model_name>_coolsolve.zip`, containing eescode + initials + sol + conf +
  debug files.  The browser's File System Access API (`showSaveFilePicker`)
  is used when available, with a download-prompt fallback.

The folder *is* the project — no proprietary format is needed.

### 7.2 ZIP Bundle Format

The ZIP uses a minimal uncompressed implementation (custom CRC32-based ZIP
writer/reader in `server.cpp` — no external library):

| Entry | Contents |
|-------|---------|
| `<name>.eescode` | Current editor content |
| `<name>.initials` | Current initial guesses |
| `<name>.sol` | Last solution (if solve was successful) |
| `coolsolve.conf` | Solver configuration |
| `debug/*` | Debug output files (if present) |

### 7.3 Session Snapshot (Back Button)

Before each destructive operation (New, Open, Upload), the current session
state is saved to a `SessionSnapshot`.  The `POST /api/v1/back` endpoint
restores it.  This provides one-level undo.

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
(T/P/H/S/D/Q/V/L/C) and `VariableInfo.inferredFluid` for every variable that
appears as an argument to a thermodynamic function.

After the solve thread completes, the server extracts this data into the
session's `inferredVariables` and `modelFluids` fields (~0.1 ms overhead).
This extraction is the only new per-solve cost — negligible.

### 8.4 Fluid Types

The `FluidRegistry` (`src/fluids.cpp`) registers 90+ real fluids, 11 ideal
gases, 2 incompressibles, and 2 mixtures.  Only fluids of type `FluidType::Real`
have saturation curves and are enabled in the diagram fluid selector.

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
# API calls proxied to port 8550 via vite.config.ts
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

When `OFF`, the binary is identical to the pre-GUI version in all respects.

### 9.5 SPA Fallback Routing

The server's catch-all route (`GET /.*`) returns `index.html` (HTTP 200) for
any non-API path not found in embedded assets.  This supports client-side
routing if it is ever added.

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
    │   └── types.ts               # TypeScript interfaces for all API responses
    ├── stores/
    │   ├── modelStore.ts          # Zustand: eescode, initials, sol, conf, parse errors, solve result, console log
    │   └── uiStore.ts             # Zustand: theme, active tabs, panel visibility (persisted to localStorage)
    ├── languages/
    │   └── ees.ts                 # Monaco Monarch language definition for EES
    └── components/
        ├── Toolbar.tsx            # Top toolbar with all actions + keyboard shortcuts
        ├── CodeEditor.tsx         # Monaco editor wrapper with EES language + debounced parse
        ├── VariableTable.tsx      # Variable table with search filter + editable Initial column
        ├── ArrayTable.tsx         # Spreadsheet view for array variables (var[i])
        ├── ConfigEditor.tsx       # Form-based coolsolve.conf editor (9 groups, 30+ fields)
        ├── DebugViewer.tsx        # Debug output file list + content viewer
        ├── Console.tsx            # Auto-scrolling console with color-coded log lines
        ├── ThermoDiagram.tsx      # Plotly thermo diagrams (T-s, P-h, h-s, T-h) with overlays + export
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
| `POST /api/v1/parse` with EES code | ✅ Returns equation/variable counts |
| `POST /api/v1/files/open` with example file | ✅ Loads file + discovers companion files |
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
| `POST /api/v1/coolprop/saturation` — Water, 50 pts | ✅ `critical.T≈647.1 K`, `computeTime_ms≈67` |
| `POST /api/v1/coolprop/saturation` — R134a, 100 pts | ✅ `critical.T≈374.2 K` |
| `POST /api/v1/coolprop/saturation` — invalid fluid | ✅ HTTP 400 |
| `POST /api/v1/coolprop/saturation` — Air (ideal gas) | ✅ HTTP 400 "no dome" |
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
| File system access | Online sessions use isolated temp directories; `/files/open` works with uploaded content only |
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

### 14.3 Sensitivity Analysis (Phase 3)

A sweep facility that identifies imposed variables (variables appearing alone
on one side of an equation, e.g. `T_in = 25`), allows defining a sweep range,
and batch-solves the system:

| Item | Description |
|------|-------------|
| Sweep API | `POST /api/v1/sweep` — define input variable, range [min, max, steps], output variables; returns table of results |
| Frontend form | Swept-variable selector, range inputs, output-variable checkboxes |
| Plotly chart | Output variables vs. swept input; one line per output variable |
| CSV export | Download sweep results as CSV |
| 2D sweeps | Two swept variables → heatmap (future extension of Phase 3) |

### 14.4 Thermodynamic Diagram Enhancements

These are extensions to the already-implemented `ThermoDiagram.tsx`:

| Item | Description |
|------|-------------|
| Unit conversion | Display properties in user-friendly units (°C, kPa, kJ/kg, kJ/(kg·K)) instead of raw SI |
| Isolines | Isotherms, isobars, iso-quality lines on the saturation dome |
| Multiple fluids | Overlay domes for two fluids on the same chart for comparison |
| CoolProp.js (WASM) | Compute missing state properties in the browser without backend calls |
| Interactive hover | Show property values on hover over dome or state points |
| Diagram state persistence | Save/restore diagram configuration in the ZIP bundle |

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
- EES source code
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
