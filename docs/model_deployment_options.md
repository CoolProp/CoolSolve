# Embedding Predefined Models in an API or a UI - Options Evaluation

## Goal

Make predefined CoolSolve models callable from other applications, as a calculation service (API) or
compiled into a user interface, on Windows, macOS and Linux. A caller provides the inputs of a model,
optionally guess values and limits for its variables, and receives the solution together with warnings and
solver statistics.

**Hard requirement: the models are protected.** A user of the API or the UI must not be able to recover the
equations, constants, lookup tables or initial values of a predefined model. The models are either run where
the user has no access to them, encrypted, or compiled.

This page is an evaluation of the options; none of the new parts (worker mode, protected build, encrypted
bundles, code generation) is implemented yet. See [Section 9](#9-work-items) for the list.

---

## 1. What exists today

| Building block | Where | Use for a service |
|---|---|---|
| CLI with JSON output | `main.cpp`: `coolsolve model.eescode -o result.json` | Status, variables, per-block results, timing and diagnostics. Picks up `model.initials`, `coolsolve.conf` and the lookup CSVs from the model folder. Parses, analyses and warms up CoolProp on every call. |
| REST server | `coolsolve --gui [port] --no-browser`, routes in `src/server.cpp`, reference in [GUI & REST API](gui.md) | Works, but built for the browser GUI: session state, asynchronous solve (`POST /api/v1/solve`, then `GET /api/v1/solve/result`), several calls per result. Exposes the model source (`GET /api/v1/files/eescode`, debug files, bundle download), so it is not suitable for protected models. |
| C++ library | CMake target `coolsolve_lib`; `CoolSolveRunner` (`include/coolsolve/runner.h`) or `IR`, `StructuralAnalyzer`, `Solver` directly | Base of every option below. |
| Parse once, solve many | Parametric study fast path, `src/server.cpp` (`cachedAnalysis`) | The structural analysis is done once and reused; each point only changes the input values and warm-starts from the previous solution. This is the core of a "presolved model". |
| Parse from memory | `EESParser::parse(source)`, `include/coolsolve/parser.h` | Needed to load a decrypted model without writing it to disk. |
| Variable bounds | `lowerBound` / `upperBound` in the variable record, `include/coolsolve/ir.h` | The fields exist, but no solver reads them and the `.initials` format has no syntax for them. |
| Embedding files in the binary | `cmake/embed_assets.cmake`, `cmake/embed_library.cmake` | Already used for the GUI assets and the model library; reusable to embed model bundles. |

### Structural limitation of a presolved model

The block decomposition depends on **which** variables are specified, not on their values. A presolved model
therefore has a fixed list of input variables. Changing the values is free; specifying another variable (for
example the condensing pressure instead of the condensing temperature) requires a new structural analysis.
Each model declares its inputs and public outputs in a manifest.

### Limits

Two levels, the first one cheap:

1. **Post-solve check**: a solution outside the limits is reported as a warning or an error.
2. **Enforced limits**: Newton, Levenberg-Marquardt and trust-region keep trial points inside the box. More
   work, but also helps convergence.

---

## 2. Model protection

Protection depends on where the model is decrypted or executed, and who controls that machine.

| Approach | Casual user (opens the folder) | Technical user (hex editor, `strings`) | Reverse engineer (debugger, memory dump) |
|---|---|---|---|
| Plain model files | Not protected | Not protected | Not protected |
| Encrypted bundle, key built into the program | Protected | Protected | **Not protected**: the key and the decrypted model are in memory |
| Encrypted bundle, key from a licence file or dongle | Protected | Protected | Slowed down, but memory can be dumped after decryption |
| Compiled model (generated C++, symbols stripped) | Protected | Protected | Largely protected: requires decompiling optimised machine code |
| Model runs only on your server | Protected | Protected | Protected: the model never leaves your machine |

Encrypting JSON on a machine the user controls is a deterrent, not a guarantee. Only server-side execution or
compilation hold against a determined user. The one leak server-side execution cannot prevent is black-box
fitting of a surrogate model from many API calls; rate limits and licence terms cover that.

### 2.1 Leaks that encryption does not fix

These apply to every option, the server included:

- **All variables in the result.** A solve result today returns every variable, intermediate ones included.
  A protected model returns only the outputs listed in its manifest.
- **Diagnostics quote the model.** A `Diagnostic` carries a message and a source line
  (`include/coolsolve/diagnostic.h`), and many messages name equations, variables and CoolProp calls. For a
  protected model they are mapped to generic codes and public variable names.
- **Debug output and reports.** Debug mode (`-d`) writes the full analysis to disk; the LaTeX report prints
  every equation. Both are absent from a protected build.
- **GUI server routes.** `GET /api/v1/files/eescode`, `/api/v1/debug/file`, the bundle download and similar
  routes return the model. A protected service has its own minimal route set.
- **Files on disk.** A decrypted model is never written to disk. The parser reads from memory already, but
  `CoolSolveRunner` reads the initials, lookup tables and configuration from the model folder, so an
  in-memory loader is needed.

### 2.2 Protection mechanism (a): encrypted model bundle

- One file per model containing `.eescode`, `.initials`, lookup tables, `coolsolve.conf` and `manifest.json`.
- Encrypted with AES-256-GCM (authenticated: also detects tampering); optionally signed (for example
  Ed25519) so that only models issued by you load.
- Decrypted in memory only; the existing solver runs unchanged, with all its features (tearing, symbolic
  reduction, multi-start, Try Harder).
- Key, from weakest to strongest: built into the binary; per-customer licence file; online activation;
  hardware or software licence protector (commercial products such as Wibu CodeMeter or Thales Sentinel,
  which also encrypt the program code and resist debugging).
- Effort: medium.

### 2.3 Protection mechanism (b): compiled model

- A generator writes C++ for one model from its structural analysis: explicit equations become assignments,
  each implicit block becomes a residual function with exact derivatives from the existing automatic
  differentiation, called by the existing Newton, Levenberg-Marquardt and trust-region solvers.
- Block order, symbolic reduction and tearing are decided at generation time. No parser and no equation text
  ship. Lookup tables remain data and are encrypted as in (a). Constants remain in the binary, without names
  or structure.
- Costs: every model change is a rebuild and a release; features that work on the equation tree at solve
  time (symbolic reduction, tearing, Try Harder) are only available where generation can resolve them in
  advance; the generator must cover `FUNCTION` / `PROCEDURE`, `IF`, `DUPLICATE`, arrays and lookup tables,
  with `INTEGRAL` models last.
- Effort: high.

A cheaper intermediate step is to ship the analysed model in a binary format instead of source text. On its
own this is weak (expression trees are easily printed back as equations); combined with (a) it adds little.

---

## 3. Integration routes

| # | Route | How the caller uses it | Crash isolation | Protected models | Effort |
|---|---|---|---|---|---|
| A | CLI per call | Starts `coolsolve`, reads the JSON output | Yes | No (plain files) | None |
| B | Long-running worker process (JSON lines over stdin/stdout) | Any language that can start a process | Yes, a stuck worker can be killed | With (a) or (b) | Low |
| C | Local HTTP server | HTTP on localhost | Yes | With (a) or (b), GUI routes removed | Low to medium |
| D | **Remote HTTP service** (Python FastAPI plus CoolSolve workers) | HTTP over the network | Yes | **Yes, strongest**: plain models on the server | Medium |
| E | Shared library with a C interface (`.dll` / `.dylib` / `.so`) | Direct call from C#, Python, Excel, MATLAB, LabVIEW, Delphi, Rust, Swift | No | With (a) or (b) | Medium |
| F | Python extension module (pybind11 or nanobind, as a wheel) | `import coolsolve` | No | With (a) or (b), in the compiled part only | Medium |
| G | WebAssembly (Emscripten) | JavaScript in the browser, no server | Browser tab | **No**: downloaded code and any key are visible | Medium to high |
| H | Compiled model or FMU | Function call, or a simulation tool through FMI | No | **Yes** | High |

Notes:

- **E**: expose plain C functions (`cs_open_model`, `cs_solve(handle, request_json)`, `cs_free`), not C++
  classes. A C++ binary interface built with one compiler does not reliably work with another (for example
  MSVC and MinGW); a C interface does.
- **F**: release the Python global interpreter lock during a solve. Python code is readable, so nothing secret
  lives in it. `cibuildwheel` builds the Windows, macOS and Linux wheels with CoolProp linked in.
- **G**: CoolProp has a WebAssembly build. Excluded for protected models.
- **H**: an FMU (Functional Mock-up Unit, FMI standard) is a zip with one binary per platform. A binary-only
  FMU is the usual way to hand over protected models to Simulink, Dymola, OpenModelica or Python `fmpy`.

### 3.1 UI technology and route

| UI technology | Route |
|---|---|
| Web front end (React, Vue) with a server | D |
| Desktop app around a web front end (Tauri, Electron, native web view) | B or C as a helper process, or D |
| C++ desktop (Qt, wxWidgets, Dear ImGui) | Link `coolsolve_lib`, with (a) or (b) |
| .NET desktop (WPF, WinUI, Avalonia, MAUI) | E, or B for crash isolation, or D |
| Python desktop or notebook (PySide, Streamlit, Dash, NiceGUI, Jupyter) | F, or B, or D |
| Excel | E through Excel-DNA or VBA, or D through Power Query |
| MATLAB, Simulink, Modelica | E or H (FMU) |

Wrapping the existing CoolSolve GUI in a desktop window is evaluated separately in
[desktop_app_options.md](desktop_app_options.md).

---

## 4. Packaging the predefined models

1. **Model folders next to the binary** (unprotected, or on a server only):

   ```
   models/
     heat_pump_r290/
       model.eescode
       model.initials          presolved guess values
       model-compressor.csv    lookup table, LOOKUP('compressor', ...)
       coolsolve.conf
       manifest.json           inputs, public outputs, default limits, version
   ```

2. **Encrypted bundles**, one file per model, as files or embedded in the binary with the existing embedding
   mechanism (mechanism (a)).
3. **Compiled models** (mechanism (b)).

Parsing and analysing a model happens once at load time. A precomputed binary form is not worth it unless
measurements show that startup is too slow.

### One executable or several

One executable with the models as data (folders or encrypted bundles):

- CoolProp is most of the binary size; one binary per model would duplicate it.
- Models are updated without recompiling; one binary is built and tested.
- All models run on the same solver and CoolProp version, so results are comparable.
- Each `manifest.json` carries a model version; every response returns the model version and the CoolSolve
  version (`coolsolve --version`), which gives traceability without freezing a binary per model.

Separate executables only pay off for a single model handed to a third party as one file, or with compiled
models (mechanism (b)), where a model and its binary are released together.

---

## 5. Calculation service with Python FastAPI (route D)

```
client -> FastAPI (uvicorn, 1 process)
            |  pydantic request validation, one schema per model (from manifest.json)
            |  queue, timeout, output whitelist, public diagnostics
            +-- worker 1: coolsolve --worker models/   (all models preloaded)
            +-- worker 2: coolsolve --worker models/
            +-- worker N  (N = number of CPU cores)
```

Workers are separate processes rather than an in-process Python module: a crash in CoolProp or the solver
does not take the web server down, and a solve that hangs can be killed and the worker restarted. The solver
timeout (`src/solver.cpp`, thread-local deadline) is cooperative, so it is not enough on its own. A pipe
round trip is small compared with a solve.

### 5.1 Request and response

```json
{"inputs":  {"T_evap_C": -5, "T_cond_C": 45, "dT_sh_K": 5},
 "guesses": {"m_dot_ref_kg_s": 0.05},
 "limits":  {"p_evap_bar": [0.5, 20]},
 "options": {"timeout_s": 5, "tryHarder": false}}
```

```json
{"status": "Success",
 "outputs": {"COP": 3.9, "Q_dot_cond_kW": 12.1},
 "warnings": [], "limitViolations": [],
 "stats": {"iterations": 14, "blocks": 9, "maxResidual": 1e-11, "solve_ms": 12.3},
 "verification": "passed",
 "modelVersion": "1.2.0", "coolsolveVersion": "0.4"}
```

`outputs` holds only the public outputs of the manifest. `warnings` holds public diagnostics (Section 2.1).
The statistics come from the existing result JSON (`src/server.cpp`, block and timing fields) and the
solution check.

### 5.2 Concurrency

- The solver is designed for threads: the CoolProp state cache and the solve deadline are thread-local
  (`src/evaluator.cpp`, `src/solver.cpp`), and multi-start already runs candidates in parallel.
- A solve writes its working values into the variable records, so each concurrent solve needs its own copy;
  the structural analysis is shared read-only, as in the parametric study.
- A few timing counters in `src/variable_inference.cpp` are global and not thread-safe. They only affect
  statistics, and a presolved model skips variable inference.
- The thread safety of all of CoolProp has not been verified. The safe default is one single-threaded worker
  process per CPU core, with `multiStartNumCores = 1` so that workers do not oversubscribe the cores.
- **Deterministic results**: by default a solve starts from the presolved `.initials` plus the guesses of the
  request, so the same input always gives the same result. Warm-starting from the previous solution is
  faster in sweeps but makes results depend on request order; it is an explicit option only.
- FastAPI: an `async` endpoint waits on a free worker; 503 when the queue is full, 504 with partial
  statistics on timeout. One uvicorn process owns the pool; several uvicorn processes would multiply it.

### 5.3 Lookup tables

- **Read-only during a solve**: the solver receives a `const LookupTableStore*`
  (`include/coolsolve/solver.h`), and a cell write `LOOKUP('t', 1, 1) = x` is dropped by the parser with a
  warning. Each worker loads the tables once at startup and shares them across requests and threads.
- **Naming** as today: `model-compressor.csv` is `LOOKUP('compressor', ...)` (`src/lookup_table.cpp`,
  `loadLookupTableForModel`). In a bundle, the tables are part of the encrypted file.
- **Tables sent with a request** (for example a customer compressor map): the worker builds a table store for
  that request only; the preloaded tables are untouched. Table values are data, so the presolved analysis
  stays valid, unless a table value sets a `DUPLICATE` bound or otherwise changes the equation count; such a
  request is rejected or triggers a new analysis.
- Every response can return a hash of each table used, for traceability.

---

## 6. Platform notes

| | Windows | macOS | Linux |
|---|---|---|---|
| Build | MSVC (supported); NSIS installer exists (`coolsolve.nsi`) | Clang; universal binary (arm64 and x86_64) with `CMAKE_OSX_ARCHITECTURES` | GCC; build on the oldest supported distribution, or in a manylinux container for wheels |
| Library file | `.dll` plus import `.lib` | `.dylib`, set `install_name` / `rpath` | `.so`, set `rpath` |
| Distribution | Authenticode signing, or SmartScreen warns | Signing and notarisation required; every helper binary in the `.app` is signed | AppImage, `.deb` or wheel |
| Protected build | Release build without PDB files; commercial protectors (CodeMeter, Sentinel, VMProtect) are most mature here | Hardened runtime blocks attaching a debugger to a signed app | `strip` symbols; weakest on a user-controlled machine (root, `ptrace`); best as the server platform for D |
| Web view for a desktop UI | WebView2 (built into Windows 10 and 11) | WKWebView (built in) | WebKitGTK (may need installing) |

Licences: CoolSolve and CoolProp are MIT, so static linking and closed-source UIs are fine. Eigen (MPL-2.0)
and KINSOL (BSD) are compatible; list them in the third-party notices.

---

## 7. Recommendation

1. **If calculations may run on a server: route D.** FastAPI plus CoolSolve workers; models as plain folders
   readable only by the service account; output whitelist; public diagnostics; no GUI, debug or report routes.
   Strongest protection and least work. Desktop UIs, Excel and web front ends all call this API.
2. **If calculations must also run offline on customer machines**: routes B, C or E with encrypted bundles
   (mechanism (a)), per-customer licence keys, and a commercial licence protector for valuable models. Keeps
   all solver features and can ship quickly.
3. **If offline models need strong protection, or must go into Simulink, Modelica or other simulation
   tools**: compiled models (mechanism (b)), shipped as a DLL or a binary-only FMU, starting with the most
   valuable models that the generator supports (`INTEGRAL` models last).

Route G (WebAssembly) is excluded for protected models.

---

## 8. Comparison summary

| Criterion | D: server | B/C/E + encrypted bundle | H: compiled / FMU | G: WebAssembly |
|---|---|---|---|---|
| Protection | Strong | Deterrent, stronger with a licence protector | Strong | None |
| Offline use | No | Yes | Yes | Yes |
| All solver features | Yes | Yes | Partial | Yes |
| Model update | Replace a folder on the server | Issue a new bundle | Rebuild and release | Rebuild |
| Effort | Medium | Medium | High | Medium to high |

---

## 9. Work items

Shared by all recommended routes (needed first, for the server):

- [ ] Model manifest (`manifest.json`): inputs, public outputs, default limits, model version.
- [ ] In-memory model loader (source, initials, lookup tables, configuration), no files read from the model
      folder.
- [ ] Output whitelist and public diagnostics for protected models.
- [ ] Post-solve limit check, reported as `limitViolations`.
- [ ] `coolsolve --worker <models_dir>`: load and analyse all models at startup, one JSON request and one
      JSON response per line, per-request copy of the variable records.
- [ ] Protected build option without source, debug, LaTeX report or GUI routes.

Route D:

- [ ] FastAPI app (`service/`): worker pool, per-model pydantic schemas from the manifests, timeouts,
      worker restart, health endpoint.

Routes B, C, E (offline):

- [ ] Encrypted and signed model bundles, bundle tool, key handling.
- [ ] Shared library with a C interface.

Route H:

- [ ] Code generator from the structural analysis; FMU packaging.

Optional:

- [ ] Limits enforced inside Newton, Levenberg-Marquardt and trust-region.
- [ ] Explicit warm-start option for parameter sweeps.
