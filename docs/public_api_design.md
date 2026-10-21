# Public REST API: Review and Design Proposal

**Status:** proposal (2026-10-21). Phase 0 (security fixes) is in pull request
[CoolProp/CoolSolve#6](https://github.com/CoolProp/CoolSolve/pull/6), not merged yet;
the later phases are not started.
**Scope:** the HTTP API served by `coolsolve --gui` (`src/server.cpp`, client in
`gui/src/api/`), and what it takes to offer a stable API that other applications can
build on, as a process or as a container image.

This document reviews the current API (section 2), proposes an architecture with two
API surfaces on one engine (section 3), designs the public API (section 4), including
its use behind reverse proxies (section 4.10) and as a container image (section 4.11),
describes how the GUI keeps its freedom (section 5), lays out the phases to get there
(section 6) and lists the decisions still open (section 7).

---

## 1. Summary

The current `/api/v1` is a good **GUI backend** and a poor **public API**. It was built
for one client, the embedded React app, and it shows: the state lives in a server-side
session, results are "the last solve of this session", progress and results arrive on a
shared event channel, and the request and response shapes are defined only by the
TypeScript types of the GUI. Nothing is wrong with that for the GUI, but an external
application needs something different: send a model and some inputs, get the outputs
back, and rely on the shape not changing under it.

Recommendations, in order:

1. **Fix the security findings first** (section 2.3, F1), independently of the API work.
   One of them is a remote command execution through `POST /api/v1/latex/compile`.
   Done in pull request #6 (Phase 0), together with guards that the public API reuses:
   loopback by default, a cross-site guard, a host-name allow-list against DNS
   rebinding.
2. **Extract a service layer** from `src/server.cpp`: plain C++ functions and structs
   (`parseModel`, `solveModel`, `runStudy`, one JSON serializer) that know nothing about
   HTTP or sessions. Today the logic sits inside the request handlers, so a second API
   would copy it.
3. **Split the HTTP surface in two**, served by the same binary:
   - `/api/gui/...`: the **internal** GUI API. Session based, unversioned, free to
     change at every release because the GUI ships inside the same binary.
   - `/api/v1/...`: the **public** API. Stateless requests plus asynchronous jobs,
     defined by a committed OpenAPI file, changed only by adding things.
4. **Design the public API around one idea: a model is a function.** A request carries
   the model (or a library id), the input values that replace imposed values, and the
   outputs wanted; the response carries values with their units and a clear status.
5. **Add a headless mode** (`coolsolve --serve`) that serves only the public API, bound
   to `127.0.0.1` unless told otherwise, with a token when bound to a network address,
   and documented rules for running it behind a reverse proxy (section 4.10).
6. **Publish the headless mode as a container image** (section 4.11), so that other
   applications can run the API next to their own services without building CoolSolve
   or CoolProp.

---

## 2. Review of the current API

### 2.1 What it is

`startServer()` in `src/server.cpp` (3,555 lines) registers about 45 routes on a
cpp-httplib server. Grouped by purpose:

| Group | Routes | Nature |
|---|---|---|
| Lifecycle | `health`, `warmup`, `keepalive`, `initial-file` | GUI plumbing |
| Session model state | `files/eescode`, `files/initials`, `files/conf`, `files/sol` (GET/PUT), `new`, `back`, `model-name`, `update-guesses` | editor buffers held on the server |
| Files | `files/open`, `files/bundle`, `files/upload` | ZIP round trip, file association |
| Analysis | `parse`, `solve`, `solve/stream` (SSE), `solve/result`, `solve/cancel`, `variables`, `variables/inferred` | the engine |
| Studies | `parametric`, `parametric/result`, `integral/result`, `integral/csv` | the engine |
| Tables | `tables`, `tables/:name` (GET/PUT/DELETE) | session data |
| Fluids | `coolprop/fluids`, `coolprop/saturation` | stateless, mostly |
| Library | `library/menu`, `library/open`, `/library/*` pages | embedded model collection |
| Reports | `latex/report`, `latex/compile`, `debug/files`, `debug/file` | GUI features |
| Operations | `stats/log` | usage dashboard |

A typical solve from a client is: `PUT /files/eescode`, `PUT /files/initials`,
`POST /solve` (returns `{"status":"started"}`), open `GET /solve/stream` and wait for
the `done` or `error` event, then `GET /solve/result` or read the result from the event.
The session is chosen by the `coolsolve_session` cookie, or by the
`X-CoolSolve-Session` header (or `?sid=`) for tabs that opted into their own session.

### 2.2 What works well

These should be kept and carried into the new design:

- **A versioned prefix** (`/api/v1`) and JSON bodies throughout.
- **Asynchronous solving with progress and cancellation**: the solve runs in a worker
  thread, block progress is pushed as server-sent events (SSE), and cancellation goes
  through `SolverOptions::cancelToken`. That is the right model for solves that take
  seconds to minutes.
- **Session isolation per tab**, with idle sessions swept after 24 h and session ids
  validated (32 hex digits) before they name a temp folder.
- **A careful file policy**: `files/open` refuses every path except the one given on
  the command line, `debug/file` rejects separators and `..`, plot file names in the
  LaTeX build are sanitised.
- **Library models with stable ids** (`CSL-0012`): a natural, citable input for a
  public API.
- **Diagnostics with codes** (`P002`, `C001`, ...) and line numbers: these are already a
  machine-readable error vocabulary.
- **Python API tests** (`tests/test_*_api.py`) for tables, library, parametric studies,
  usage log and integral results.

### 2.3 Findings

Each finding says why it matters for an external integration. Findings F1 to F3 are
defects today, for the GUI too; F4 to F10 are reasons why the current API cannot simply
be declared public.

**F1. Security (fix now, independently of the API work).** *Addressed in pull request
#6* (Phase 0): only `pdflatex`, `xelatex` and `lualatex` are run, plot names are plain
`.png` names, the server listens on `127.0.0.1` unless started with `--host`, the CORS
header is gone and cross-site `POST`/`PUT`/`DELETE` are refused, every request must
name the server in its `Host` header (`localhost`, an IP address, or a name given with
`--allowed-host`; an adversarial review of the first fix showed that DNS rebinding
passed the cross-site guard), and the statistics show client networks instead of
addresses. The description below is the state before that pull request.

- **Command injection in `POST /api/v1/latex/compile`** (`src/server.cpp`, the
  `compiler` variable around line 3144). The `compiler` field of the JSON body, or the
  `latexCompiler` key of the session's `coolsolve.conf` (which any client can set with
  `PUT /files/conf`), is concatenated into a shell command and run with `popen`. A
  client that first solves any model with `enableLatexReport = true` can then run an
  arbitrary command as the server user. Combined with the next two points, this is
  reachable from the local network and from any web page open in the user's browser.
  Fix: accept only a fixed list of compiler names (`pdflatex`, `xelatex`, `lualatex`),
  never a value from the request body, and run it without a shell.
- **The server listens on all interfaces** (`bind_to_port("0.0.0.0", ...)`). On a
  desktop this exposes the GUI API, without authentication, to everyone on the same
  network. The default should be `127.0.0.1`, with an explicit `--host` option for the
  server deployment (which sits behind Apache anyway, see
  `docs/deployment_ubuntu_apache.md`).
- **`Access-Control-Allow-Origin: *`** on every response. It exists for the Vite dev
  server (port 5173), but it also lets any web site call the local API from the user's
  browser. It should be enabled only in dev mode, and only for the dev origin.
- **`GET /api/v1/stats/log` returns the top client IP addresses** to any caller. On the
  public deployment this is personal data; it should require an operator token or be
  limited to loopback.

**F2. Results are implicit session state.** A result is "the last solve of the
session", with no id. A client cannot run two solves in one session (the second gets
409), cannot tell which request a result belongs to, and loses it as soon as anything
else solves in that session. Integrations that call the model many times (an
optimiser, a spreadsheet, a plant model) need either a synchronous call or a job id.

**F3. Concurrency defects.**

- The solve thread writes `solContent`, `inferredVariables`, `modelFluids`,
  `latexReportContent` and `debugDir` of the session without a lock, while request
  handlers (for example `GET /files/sol`, `GET /variables/inferred`) read them. The GUI
  rarely triggers this because it waits for the `done` event; a script that polls will.
- `POST /solve` and `POST /parametric` test `session.solving` at the top and set it
  much later, so two concurrent requests can both start.
- The SSE queue is consumed by whoever reads it (`waitForEvents` moves the events out).
  Two listeners on one session each see part of the events. Solves and parametric
  studies share the same channel.
- There is no limit on concurrent work across sessions: each solve or study starts a
  detached thread. Ten clients mean ten solver threads.

**F4. Process-wide settings change per request.** `applyCoolPropConfig()` sets the
CoolProp superancillary flag for the whole process from the session's configuration. Two
clients with different settings change each other's results. A public API must either
make this a per-job setting (if CoolProp allows it per `AbstractState`) or not expose it.

**F5. No written contract.** The request and response shapes exist only in
`gui/src/api/types.ts` and in the handler code. Some consequences:

- Naming is mixed: `totalIterations`, `blockResults` (camelCase) next to `parse_ms`,
  `solve_ms` (snake_case) in `timing`.
- Errors have three shapes: an HTTP error with `{"error": "..."}`, an HTTP 200 with
  `"success": false`, and an SSE event of type `error` that may carry a full result. A
  failed solve (the model did not converge) and a failed request (bad JSON) are not
  clearly separated. Status codes for the same situation differ (`update-guesses`
  answers 400, `back` answers 409 for "nothing to do").
- `status` is documented as the solver enum (`Success`, `MaxIterations`, ...) but the
  server also writes `SolutionCheckFailed`, which is not in the enum.
- Units are free text in brackets, and the inferred units of `/parse` contain non-ASCII
  characters (degree sign, middle dot, superscript 3). The values returned by `/solve`
  carry no unit at all (units are only in `/parse` and `/variables/inferred`).
- Arrays are flattened into names (`"T[3]"`), so a client has to parse names to rebuild
  them.

**F6. GUI-only routes are mixed with engine routes.** `initial-file`, `keepalive`,
`back`, `model-name`, `update-guesses`, `debug/*`, `latex/*`, `files/bundle` and
`stats/log` make sense only for the GUI or the operator. Publishing them would freeze
GUI features into a contract.

**F7. Parametric studies edit the source text.** `POST /parametric` replaces
`var = value` in the source with a case-insensitive regex for every grid point. This is
fragile (the imposed equation must have exactly that form) and is the one feature most
integrations want in a robust form: "evaluate this model for these inputs". It should
become an override applied to the IR, shared by studies and by the public API.

**F8. Business logic lives in the handlers.** Building the result JSON, writing the
`.sol` text, running the solution check, extracting inferred variables, detecting
imposed variables (in `/parse`), building and reordering the parametric grid: all of it
is inside lambdas in `startServer()`. A public API, a Python binding or the CLI would
each have to copy it, and the copies would drift.

**F9. The server cannot run without the GUI.** The HTTP server is compiled only with
`COOLSOLVE_GUI` and always serves the SPA and opens a browser unless told not to. An
integration wants a headless service with only the API.

**F10. Contract tests are missing for the core path.** The Python tests cover tables,
library, parametric studies, usage log and integral results, but not
`/solve` + `/solve/stream` + `/solve/result`, which is the path every integrator uses.

**F11. The client address is taken from any `X-Forwarded-For`.** `getClientIP()`
(`src/server.cpp`) uses the first entry of `X-Forwarded-For` whenever the header is
present, whoever sent it. Behind Apache this is the real client; without a proxy, any
client can write any address into the usage log. Harmless for statistics, but a public
API that limits or audits per client must believe the header only when the request
comes from a proxy it trusts (section 4.10).

---

## 3. Target architecture

```
                    +--------------------------------------------+
  React GUI  ---->  |  /api/gui/...   internal, session based    |
                    |  (unversioned, changes with the GUI)       |
                    +---------------------+----------------------+
                                          |
  Other apps ---->  +---------------------+----------------------+
  (Python, Excel,   |  /api/v1/...    public, stateless + jobs   |
   TESPy, web)      |  (OpenAPI contract, additive changes only) |
                    +---------------------+----------------------+
                                          |
                    +---------------------v----------------------+
                    |  Service layer (new, no HTTP, no sessions) |
                    |  parseModel, solveModel, runStudy,         |
                    |  JobManager, api_json (one serializer)     |
                    +---------------------+----------------------+
                                          |
                    +---------------------v----------------------+
                    |  Engine (exists): parser, IR, inference,   |
                    |  structural analysis, solvers, runner      |
                    +--------------------------------------------+
```

### 3.1 Service layer

A new header `include/coolsolve/service.h` (and `src/service/*.cpp`) with plain structs
in and out. Sketch:

```cpp
// What a model is, independent of where it came from (editor, ZIP, library, request).
struct ModelInput {
    std::string source;                              // .eescode text
    std::map<std::string, double> guesses;           // initial values
    std::map<std::string, std::string> tablesCsv;    // lookup tables, name -> CSV text
    SolverOptions options;                           // already validated
};

// Values that replace imposed values ("T_amb = 25") without editing the source.
using Overrides = std::map<std::string, double, CaseInsensitiveLess>;

struct SolveOutcome {
    std::string status;            // "solved", "failed", "cancelled", "timeout", "invalid_model"
    SolveResult result;            // engine result (values, blocks, solver status)
    DiagnosticCollector diagnostics;
    CoolSolveRunner::PipelineTiming timing;
    std::vector<VariableInfo> variables;   // name, unit, unit source, inferred property and fluid
    std::optional<IntegralSolveResult> integral;
};

ParseOutcome parseModel(const std::string& source);
SolveOutcome solveModel(const ModelInput& model, const Overrides& overrides,
                        const std::atomic<bool>* cancel, ProgressCallback progress);
StudyOutcome runStudy(const ModelInput& model, const StudyRequest& study, ...);
```

Rules for this layer:

- No file system access for model data: `CoolSolveRunner` currently takes a file path,
  so the service layer first gets a constructor or factory that takes the source text
  (the server writes `model.eescode` to a temp folder only because of this).
- No session, no HTTP, no JSON in the function signatures. JSON is produced by one
  module (`src/service/api_json.cpp`) that both HTTP surfaces use, so a field has one
  name everywhere.
- Overrides are applied to the IR (replace the right-hand side of the imposed
  equation), with a clear error when the name is not an imposed variable (F7).
- The solution check, the `.sol` text and the inferred variables are computed here, once.

The GUI handlers become thin: read the session, call the service, store the outcome
in the session under one lock, answer. This also fixes F3: the worker thread builds a
complete outcome and swaps it into the session in one locked step.

### 3.2 Job manager

Long work (a solve, a study, a dynamic model) runs as a **job** in a bounded worker
pool, shared by both surfaces:

- a fixed number of workers (default: number of cores, configurable), and a bounded
  queue; when the queue is full the API answers `429 Too Many Requests` with
  `Retry-After`;
- each job has an id, a state (`queued`, `running`, `succeeded`, `failed`,
  `cancelled`), a deadline (default timeout, configurable per request up to a server
  maximum), and its own event buffer, so any number of SSE listeners each get every
  event (F3);
- finished jobs are kept for a configurable time (default 1 h), then dropped.

GUI sessions use the same job manager: a session simply remembers the id of its
current job.

---

## 4. Public API design (v1)

### 4.1 Principles

1. **A model is a function.** The core call takes a model and input values and returns
   output values. Everything else is support (fluids, library, parsing for validation).
2. **Requests are self-contained.** No cookies, no server-side editor state. The only
   server-side state is jobs, and they have ids.
3. **The contract is a file.** `api/openapi-v1.yaml` is committed, reviewed like code,
   served at `GET /api/v1/openapi.json`, and checked by contract tests (section 4.8).
4. **A solver failure is a result, not an HTTP error.** HTTP 4xx/5xx means the request
   could not be processed; a model that did not converge is a `200` (or a finished
   job) with `"status": "failed"` and diagnostics.
5. **Units travel with values**, as ASCII strings taken from the model.
6. **Only add, never change.** Within v1, fields, endpoints and enum values are only
   added. Clients must ignore unknown fields and treat unknown enum values as
   "other". Anything else is v2 (section 4.7).

### 4.2 Resources

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/api/v1/info` | Versions (CoolSolve, CoolProp and its revision, API version, library snapshot) and server limits (max workers, max timeout, max study points, max body size) |
| `GET` | `/api/v1/health` | Liveness and readiness (`coolpropReady`) |
| `GET` | `/api/v1/openapi.json` | The contract |
| `POST` | `/api/v1/models/check` | Parse and analyse without solving: equation and variable counts, square or not, imposed variables (the possible inputs), variables with units, diagnostics |
| `POST` | `/api/v1/jobs` | Start a job (`kind`: `solve` or `study`). Answers `202` with the job, or `200` with the finished job when `?wait=<s>` is given and the job ends in time |
| `GET` | `/api/v1/jobs/{id}` | Job state and, once finished, its result |
| `GET` | `/api/v1/jobs/{id}/events` | SSE progress stream for that job (every listener gets every event) |
| `DELETE` | `/api/v1/jobs/{id}` | Cancel a queued or running job |
| `GET` | `/api/v1/fluids` | Fluid registry: name, CoolProp name, type, `hasDome` |
| `GET` | `/api/v1/fluids/{name}/saturation?points=200` | Saturation curve and critical / triple point data |
| `GET` | `/api/v1/library/models` | Library catalogue (id, title, category, status) |
| `GET` | `/api/v1/library/models/{id}` | One library model as a `Model` object (section 4.3), ready to send back in a job |

`POST /api/v1/jobs?wait=60` covers the simple synchronous case ("solve this and give me
the answer") without a second endpoint; a client that wants progress or long studies
leaves out `wait` and follows the job.

### 4.3 The `Model` object

```json
{
  "source": "T_amb = 25 [C]\nQ_dot_cond = ...",
  "guesses": { "T_cond": 45.0, "p_evap": 3.2e5 },
  "tables":  { "cop_map": "T_cond,COP\n30,4.1\n40,3.4\n" },
  "options": { "tolerance": 1e-9, "maxIterations": 200, "deepSearch": false }
}
```

or, for a library model, a reference instead of the source:

```json
{ "library": { "id": "CSL-0012" } }
```

(with an optional `"file"` for another `.eescode` of the same library folder, as in
`library/open` today)

- `options` is a **whitelist** of `SolverOptions` keys documented in the contract
  (tolerances, iteration limits, pipeline choice, `deepSearch`, timeout). Keys that
  touch the process (`coolpropConfig`, `latexCompiler`, debug output paths) are not
  accepted (F1, F4).
- `guesses` is structured, not `.initials` text. The text format stays a GUI and CLI
  concern.

### 4.4 Solve job

Request:

```json
POST /api/v1/jobs?wait=60
{
  "kind": "solve",
  "model": { "library": { "id": "CSL-0012" } },
  "inputs":  { "T_amb": 30.0, "m_dot_sink": 0.5 },
  "outputs": ["COP", "Q_dot_cond", "T[1..5]"],
  "timeout_s": 120
}
```

- `inputs` replace imposed values (`T_amb = 25` becomes `T_amb = 30`) through the
  service-layer overrides. A name that is not an imposed variable is a `422` with the
  list of valid inputs (from `models/check`). Values are in the units the model
  declares for that variable.
- `outputs` is optional; without it all variables are returned.

Response (finished job):

```json
{
  "id": "j_7f3c...",
  "kind": "solve",
  "state": "succeeded",
  "createdAt": "2026-10-21T09:12:03Z",
  "finishedAt": "2026-10-21T09:12:04Z",
  "result": {
    "status": "solved",
    "solverStatus": "Success",
    "variables": {
      "COP":        { "value": 3.42,  "unit": "-",  "unitSource": "declared" },
      "Q_dot_cond": { "value": 12.8,  "unit": "kW", "unitSource": "declared" },
      "T":          { "index": [1, 2, 3, 4, 5], "value": [4.1, 61.0, 45.0, 45.0, 4.1],
                      "unit": "C", "unitSource": "inferred" }
    },
    "strings": { "Fluid$": "R290" },
    "diagnostics": [
      { "severity": "warning", "code": "C001", "message": "...", "line": 12 }
    ],
    "statistics": { "equations": 84, "blocks": 31, "largestBlock": 6, "iterations": 57 },
    "timing": { "parseMs": 1.2, "analysisMs": 3.4, "solveMs": 812.0, "totalMs": 820.1 },
    "engine": { "coolsolve": "0.5.0", "coolprop": "8.1.0dev", "coolpropRevision": "..." }
  }
}
```

- `state` is the job (did it run to the end); `result.status` is the engineering answer
  (`solved`, `failed`, `cancelled`, `timeout`, `invalid_model`); `solverStatus` keeps
  the detailed solver enum for diagnosis, and the solution check gets its own documented
  value instead of the ad-hoc `SolutionCheckFailed`.
- Arrays are grouped (`index` and `value` lists) instead of flattened names (F5).
- Names keep the spelling of the model; lookups are case-insensitive, as in the
  language.
- Non-finite numbers are returned as `null` with a diagnostic, never as invalid JSON.
- Dynamic (`INTEGRAL`) models add `result.integral` with the columnar table already
  produced by `integralResultToJSON`.
- Block-level detail (per-block iterations and residuals) is available with
  `?detail=blocks`, so the default answer stays small.

### 4.5 Study job

```json
{
  "kind": "study",
  "model": { ... },
  "sweep": [
    { "name": "T_amb", "values": [0, 5, 10, 15] },
    { "name": "m_dot_sink", "values": [0.3, 0.5] }
  ],
  "outputs": ["COP"],
  "warmStart": true
}
```

The result holds one row per grid point: the inputs, `status`, and the requested
outputs. This is today's parametric study, moved onto overrides (F7) and onto the job
manager. The server limits the number of points (`info.limits.maxStudyPoints`).

A later, additive extension is `"points": [{...}, {...}]` for an explicit list of input
sets (design of experiments, Monte Carlo), which the grid cannot express.

### 4.6 Errors

HTTP errors use one shape, the problem details format of RFC 9457:

```json
{
  "type": "/api/v1/problems/unknown-input",
  "title": "Unknown input variable",
  "status": 422,
  "detail": "'T_ambient' is not an imposed variable of this model",
  "validInputs": ["T_amb", "m_dot_sink", "eta_is"]
}
```

| Status | When |
|---|---|
| `400` | Body is not valid JSON or does not match the schema |
| `401` / `403` | Missing or invalid token (headless mode on a network address) |
| `404` | Unknown job, fluid or library model |
| `409` | Cancel of a job that has already finished |
| `413` | Body larger than `info.limits.maxBodyBytes` |
| `422` | Valid JSON that cannot be run: parse errors (with the diagnostics), unknown input, non-square system |
| `429` | Job queue full, with `Retry-After` |
| `503` | CoolProp not ready yet, or the server is shutting down |

The `type` values are part of the contract and listed in the OpenAPI file.

### 4.7 Versioning and stability policy

- The path carries the major version (`/api/v1`). `GET /api/v1/info` carries the
  minor version (`"apiVersion": "1.3"`), raised when something is added.
- Allowed within v1: new endpoints, new optional request fields, new response fields,
  new enum values, new diagnostic codes, new error `type` values.
- Not allowed within v1: removing or renaming a field, changing a type or a unit
  convention, making an optional field required, changing the meaning of a status.
- A v2 runs next to v1 for at least two releases. Deprecated endpoints answer with
  `Deprecation` and `Sunset` headers (RFC 9745, RFC 8594) and are listed in
  `docs/versions.md`.
- The API version is independent of the CoolSolve version: CoolSolve 0.7 can still
  serve API 1.2.
- Before the first stable release, the public API is marked **beta** in `info` and in
  the docs for one release cycle, so the contract can be corrected with real users.

### 4.8 Contract and tests

- `api/openapi-v1.yaml` is the source of truth. A change to the public API is a change
  to this file first, reviewed with the code.
- `tests/test_public_api.py` validates every response against the schema (with
  `openapi-core` or `jsonschema`), on a set of example models (`examples/*.eescode`)
  and library ids. It covers the solve path that is untested today (F10).
- A **golden-response** test stores the JSON answer for a few models and fails on any
  change that is not an addition. This is what catches accidental breaking changes.
- The contract can generate a typed Python client; an example notebook and a `curl`
  walkthrough go into the docs.

### 4.9 Security and deployment

- **Headless mode**: `coolsolve --serve [--host 127.0.0.1] [--port 8551]` serves only
  `/api/v1`, no SPA, no browser, no GUI routes.
- **Binding**: `127.0.0.1` by default in both modes; a network address only with
  `--host` (as introduced in Phase 0).
- **Host names**: the Phase 0 guard applies unchanged: a request must name the server
  in its `Host` header (`localhost`, an IP address, the `--host` name, or a name given
  with `--allowed-host`). This keeps DNS rebinding out of the public API too.
- **Authentication**: a bearer token (`Authorization: Bearer ...`, read from
  `--api-token-file`, compared in constant time). `--serve` on a network address
  refuses to start without a token unless the operator says explicitly that a proxy
  authenticates (`--auth=proxy`, section 4.10). The GUI mode keeps its Phase 0 warning
  instead, so that existing deployments do not stop.
- **Cross-site requests and CORS**: requests without `Origin` (scripts, other
  services: the usual integration) are not affected. Browser applications on another
  origin need CORS: `--cors-origin <origin>` (repeatable) answers their preflight and
  adds them to the origins the Phase 0 cross-site guard accepts; every other origin is
  still refused. A token in the `Authorization` header cannot be sent cross-site without
  a preflight, so an allowed origin cannot be abused by a third site.
- **Limits**: maximum body size, workers, queue length, job timeout, study points and
  job retention, all reported in `GET /api/v1/info`.
- **No public endpoint runs an external program or reads a path chosen by the client.**
  LaTeX compilation, debug files, bundles and the usage statistics stay GUI or operator
  features.

### 4.10 Behind a reverse proxy

Most deployments put a reverse proxy (Apache, nginx, Caddy, a cloud load balancer) in
front of CoolSolve for TLS, authentication and rate limits, as
`docs/deployment_ubuntu_apache.md` does for the GUI. The public API sets these rules:

| Concern | Rule |
|---|---|
| Listening address | CoolSolve listens on `127.0.0.1` (same machine) or on the container network (section 4.11); the proxy is the only thing exposed. |
| Host name | The proxy passes the browser's or client's host name: either in `Host` (Apache `ProxyPreserveHost On`, nginx `proxy_set_header Host $host;`), with that name given to `--allowed-host`, or in `X-Forwarded-Host` (Apache adds it). Otherwise requests get `403` from the Phase 0 host guard, or from the cross-site guard for browser `POST`s. |
| Client address (F11) | `X-Forwarded-For` is believed only when the connection comes from an address given with `--trusted-proxy <address or CIDR>` (repeatable; for example `127.0.0.1`, or the proxy network of a container setup). Then the right-most entry that is not itself a trusted proxy is the client. Without `--trusted-proxy` the header is ignored and the connection address is used. The usage log and any per-client limit use this address. |
| Scheme and path prefix | The API never builds absolute URLs: `Location: /api/v1/jobs/{id}` is a path, and the OpenAPI `servers` entry is relative. So `X-Forwarded-Proto` is not needed, and a sub-path deployment (`/coolsolve/api/v1/...`, deployment guide section 10) works when the proxy strips the prefix, as the GUI already does. |
| Server-sent events | Job events (`/api/v1/jobs/{id}/events`) must not be buffered: CoolSolve sends `X-Accel-Buffering: no` (nginx) as today, Apache needs `SetEnv proxy-sendchunked 1` (as in the deployment guide). The stream sends a comment line (`: keep-alive`) every 15 s, so that proxies with an idle timeout (`ProxyTimeout`, `proxy_read_timeout`, load balancers, often 60 s) do not cut a long solve. Clients reconnect with `Last-Event-ID` and get the events they missed from the job's buffer. |
| Request size | The proxy limit (`LimitRequestBody`, `client_max_body_size`) is at least `info.limits.maxBodyBytes`, so that a too-large model gets CoolSolve's `413` problem details, not a proxy error page. |
| Long synchronous calls | `POST /api/v1/jobs?wait=<s>` keeps the connection open for up to `wait` seconds; the proxy timeout must be longer. The server caps `wait` (`info.limits.maxWaitSeconds`, default 60 s) so that it stays under common proxy defaults; longer work is followed as a job. |
| Authentication | Either CoolSolve's bearer token, or the proxy (basic auth, OAuth proxy, mutual TLS) with `--auth=proxy`. Not both missing: `--serve` on a network address refuses to start otherwise (section 4.9). With proxy authentication, CoolSolve listens only on loopback or on a network that only the proxy can reach. |
| Rate limits | Per-client limits belong in the proxy (`mod_ratelimit`, `limit_req`); CoolSolve's own limits (workers, queue, `429` with `Retry-After`) protect the solver whatever the proxy does. |

The deployment guide gets a section for the public API with an Apache and an nginx
example covering these rows, and a troubleshooting table for the `403` and `413`
answers.

### 4.11 Container image

An integration often wants CoolSolve as a service next to its own: in Docker Compose,
Kubernetes or a cloud container service, without building CoolSolve and CoolProp.
The headless mode is published as an image for that purpose.

**Image.** `ghcr.io/coolprop/coolsolve-api`, built from a multi-stage `Dockerfile` at
the root of the repository:

```dockerfile
# --- build stage: compilers, CMake, CoolProp sources (not shipped) ---
FROM ubuntu:24.04 AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git python3 ca-certificates && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
ARG COOLSOLVE_LIBRARY_TAG=main
# No Node.js: the API image does not embed the GUI assets
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOOLSOLVE_BUILD_TESTS=OFF \
          -DCOOLSOLVE_LIBRARY_TAG=${COOLSOLVE_LIBRARY_TAG} \
 && cmake --build build -j"$(nproc)" --target coolsolve

# --- runtime stage: the binary and the C++ runtime only ---
FROM ubuntu:24.04 AS api
RUN useradd --system --uid 10001 --home-dir /var/lib/coolsolve coolsolve
COPY --from=build /src/build/coolsolve /usr/local/bin/coolsolve
USER 10001
EXPOSE 8551
ENV COOLSOLVE_HOST=0.0.0.0 COOLSOLVE_PORT=8551
ENTRYPOINT ["coolsolve", "--serve"]
```

A `.dockerignore` keeps `build/`, `.fetchcontent_cache/`, `gui/node_modules/` and
`gui/dist/` out of the build context: a local `gui/dist/` would otherwise be embedded
into the API image.

A second target (`FROM api AS gui`, with the SPA embedded and `--gui`) replaces the
single-stage example of `docs/gui.md` section 12.2, so that the GUI and the API images
come from one file.

**Configuration through environment variables.** Every `--serve` option has an
environment variable, read when the option is not given on the command line:
`COOLSOLVE_HOST`, `COOLSOLVE_PORT`, `COOLSOLVE_ALLOWED_HOSTS` (comma-separated),
`COOLSOLVE_TRUSTED_PROXIES`, `COOLSOLVE_CORS_ORIGINS`, `COOLSOLVE_API_TOKEN_FILE`,
`COOLSOLVE_AUTH` (`token` or `proxy`), `COOLSOLVE_WORKERS`, and the limits of section
4.9. The usage log already follows this pattern (`COOLSOLVE_GUI_LOG`).

**Rules that are specific to containers:**

- **Listening address.** Inside the container the server must listen on `0.0.0.0`:
  `127.0.0.1` would be the container's own loopback, unreachable from outside. Exposure
  is decided when the container is started: `docker run -p 127.0.0.1:8551:8551 ...`
  keeps it on the host's loopback (for a proxy on the host), a Compose network keeps it
  reachable only by the proxy service.
- **Host names.** Clients reach the container by a service name (`coolsolve-api`), a
  host name or an IP address; names must be given in `COOLSOLVE_ALLOWED_HOSTS`, or
  requests get the Phase 0 `403`.
- **Secrets.** The token is a file (`/run/secrets/coolsolve_api_token`, a Docker or
  Kubernetes secret) given with `COOLSOLVE_API_TOKEN_FILE`, never an environment
  variable value (those end up in `docker inspect` and in logs).
- **Least privilege.** A non-root user (uid 10001), a read-only root file system and a
  `tmpfs` for the job folders (`--read-only --tmpfs /tmp`). No shell, compiler or LaTeX
  in the image: the API runs no external program (section 4.9).
- **CPU limits.** `std::thread::hardware_concurrency()` (used for the worker count and
  by the parallel solver pipeline, `src/solver.cpp`) reports the host's cores, not the
  container's CPU quota. The worker count defaults to the cgroup quota when one is set,
  and `COOLSOLVE_WORKERS` overrides it.
- **Health.** `GET /api/v1/health` is the liveness probe and the readiness probe
  (`coolpropReady`). The image has no `curl`, so Docker's `HEALTHCHECK` calls the
  binary itself: `coolsolve --health-check` requests `/api/v1/health` on the local port
  and exits 0 or 1. Kubernetes uses an `httpGet` probe directly.
- **Logs.** Requests, warnings and the usage log go to stdout/stderr by default in the
  image (`COOLSOLVE_GUI_LOG=-`), where the container runtime collects them; a file in a
  mounted volume stays possible.
- **Shutdown.** `SIGTERM` (sent by `docker stop` and Kubernetes) stops accepting jobs,
  cancels the running ones through their cancel tokens, and exits within the grace
  period.

**Example: Compose with a TLS proxy.**

```yaml
services:
  proxy:
    image: caddy:2
    ports: ["443:443"]
    volumes: ["./Caddyfile:/etc/caddy/Caddyfile:ro"]   # reverse_proxy coolsolve-api:8551
  coolsolve-api:
    image: ghcr.io/coolprop/coolsolve-api:1
    read_only: true
    tmpfs: ["/tmp"]
    environment:
      COOLSOLVE_ALLOWED_HOSTS: api.example.com,coolsolve-api
      COOLSOLVE_TRUSTED_PROXIES: 172.16.0.0/12
      COOLSOLVE_API_TOKEN_FILE: /run/secrets/coolsolve_api_token
      COOLSOLVE_WORKERS: "4"
    secrets: [coolsolve_api_token]
    healthcheck:
      test: ["CMD", "coolsolve", "--health-check"]
secrets:
  coolsolve_api_token:
    file: ./coolsolve_api_token.txt
```

**Tags and provenance.** Each release publishes `ghcr.io/coolprop/coolsolve-api` with
the CoolSolve version (`0.6.0`, `0.6`), the API major version (`1`, which moves only
within the stability policy of section 4.7) and the commit. OCI labels record the
CoolSolve version, the API version, the CoolProp revision (pinned, contributing guide
section 3.13) and the model library tag, the same facts as `coolsolve --version`. Images
are built for `linux/amd64` and `linux/arm64`.

**CI.** A GitHub Actions workflow builds the image on each pull request (not pushed),
runs `tests/test_public_api.py` against the running container, and on a release tag
pushes it to the GitHub Container Registry. The build cache keeps the CoolProp
compilation (the slow part) between runs.

**Size.** The stripped binary is about 14 MB with CoolProp and without the GUI assets
(measured on Linux x86-64, Release build), and it links only the C and C++ runtime
(`libstdc++`, `libgcc_s`, `libm`, `libc`). The runtime stage adds the `ubuntu:24.04`
base (about 80 MB). A distroless C++ base, which carries exactly these libraries, would
bring the image under 40 MB, at the cost of the tools for debugging inside it (open
decision D9).

---

## 5. Keeping the GUI flexible

The GUI keeps its own API and loses nothing:

- **Separate prefix, no promise.** GUI routes move to `/api/gui/...`. They are
  unversioned and may change at any release, because the SPA is embedded in the same
  binary (`COOLSOLVE_EMBEDDED_ASSETS`): the client and the server always match. This is
  written at the top of the GUI section of `docs/gui.md` so nobody builds on it by
  mistake.
- **Sessions are a GUI concept.** Editor buffers, `.initials` text, Back snapshot,
  model name, keepalive, ZIP bundles, debug files and LaTeX reports stay in
  `/api/gui`. None of them appears in the public contract.
- **Shared engine, shared vocabulary.** Both surfaces call the service layer and
  serialize through `api_json.cpp`. The GUI response for a solve is the public result
  object plus GUI fields (`isImposed`, `imposedValue`, block traces, `latexReportAvailable`).
  The GUI can add fields freely; the public object changes only by the rules of
  section 4.7.
- **One-way dependency.** The GUI may call public endpoints where they fit (fluids,
  saturation, library catalogue), which keeps them exercised. The public API never
  calls or depends on a GUI handler.
- **A promotion path.** A new capability lands in `/api/gui` first. When it has been
  stable through a release and an external use exists, it is added to the OpenAPI file
  and exposed under `/api/v1`, usually by mapping the same service function. Until then
  it can change freely.

Mapping of today's routes:

| Today (`/api/v1/...`) | Public `/api/v1` | GUI `/api/gui` |
|---|---|---|
| `health`, `warmup` | `health`, `info` | `health`, `warmup` |
| `parse` | `models/check` | `parse` (with editor extras) |
| `solve`, `solve/stream`, `solve/result`, `solve/cancel` | `jobs` (`kind: solve`), `jobs/{id}`, `jobs/{id}/events`, `DELETE jobs/{id}` | same names, on the session's current job |
| `parametric`, `parametric/result` | `jobs` (`kind: study`) | same names |
| `integral/result`, `integral/csv` | `result.integral` of the job | same names |
| `variables`, `variables/inferred` | `result.variables` (with units and inferred fields) | same names |
| `coolprop/fluids`, `coolprop/saturation` | `fluids`, `fluids/{name}/saturation` | GUI may use the public ones |
| `library/menu`, `library/open` | `library/models`, `library/models/{id}` | `library/menu`, `library/open` |
| `files/*`, `new`, `back`, `model-name`, `update-guesses`, `keepalive`, `initial-file`, `tables/*` | not public | unchanged |
| `debug/*`, `latex/*`, `stats/log` | not public | unchanged (after the F1 fixes) |

---

## 6. Phases

Each phase is a separate branch and pull request, reviewed before the next starts.

### Phase 0: security fixes (pull request #6, open)

Independent of everything else; worth a patch release.

- `latex/compile`: fixed compiler allow-list (`pdflatex`, `xelatex`, `lualatex`),
  plain `.png` plot names, `-no-shell-escape` on TeX Live.
- Listen on `127.0.0.1` by default, `--host` to change it.
- No CORS header; cross-site `POST`/`PUT`/`DELETE` refused (`Origin` against `Host` or
  `X-Forwarded-Host`; the Vite dev proxy sends `X-Forwarded-Host`).
- `Host` allow-list against DNS rebinding, `--allowed-host`. Server deployments behind
  Apache with `ProxyPreserveHost On` must add `--allowed-host <domain>`.
- `stats/log`: client networks instead of addresses.
- Still to do after the merge: the Docker example of `docs/gui.md` section 12.2 needs
  `--host 0.0.0.0 --no-browser`, since the default address is now unreachable from
  outside a container (fixed on the pull request branch).

### Phase 1: service layer (pure refactor)

- `CoolSolveRunner` accepts source text as well as a path.
- Move the logic of `/parse`, `/solve`, `/parametric` and the result JSON into
  `src/service/`. GUI behaviour and JSON unchanged.
- Fix the races of F3: one session lock, outcome swapped in one step, atomic
  check-and-set for "solving".
- C++ unit tests on the service functions; a Python test of the GUI solve path.

### Phase 2: overrides and jobs

- Input overrides applied to the IR; the parametric study uses them instead of the
  regex substitution (F7).
- `JobManager` with a worker pool, per-job event buffers, timeouts and retention. GUI
  solve and study run through it.

### Phase 3: route split

- GUI routes move to `/api/gui` (one constant, `API_BASE` in
  `gui/src/api/client.ts`, plus the Python tests). The old `/api/v1` GUI routes stay
  as deprecated aliases for one release, in GUI mode only.

### Phase 4: public API beta

- Write `api/openapi-v1.yaml` first; review it on its own.
- Implement `/api/v1` on the service layer, `--serve` headless mode, token, limits,
  the reverse-proxy rules of section 4.10 (`--trusted-proxy`, `--auth=proxy`, SSE
  keep-alive, capped `wait`), environment variables for every option.
- Contract and golden tests; docs page `docs/public_api.md` with `curl` and Python
  examples; entry in the README documentation table.
- Mark as beta in `info` and in `docs/versions.md`.

### Phase 4b: container image

- `Dockerfile` (targets `api` and `gui`), `coolsolve --health-check`, cgroup-aware
  worker count, log to stdout, `SIGTERM` handling (section 4.11).
- GitHub Actions workflow: build on pull requests, test the container with
  `tests/test_public_api.py`, push multi-arch images to `ghcr.io` on release tags.
- Docs: a container section in the deployment guide (`docker run`, Compose with a TLS
  proxy, Kubernetes probes), `docs/gui.md` section 12.2 pointing to the `gui` target.
- Can start as soon as Phase 4 has a working `--serve`; the beta label covers the image
  too.

### Phase 5: stable

- After one release cycle with feedback, set `apiVersion` to `1.0`, drop the beta label
  and apply the policy of section 4.7 from then on.

### Later, optional

The service layer is also the natural base for an **in-process** integration: a small C
API (`coolsolve_c.h`) and a Python module built on it, for applications that want to
call the solver without running a server (for example a plant simulation calling a
component model thousands of times). It would expose the same `Model`, inputs and
result concepts as the REST API, so a user can move between the two.

---

## 7. Open decisions

| # | Decision | Proposal |
|---|---|---|
| D1 | Path of the public API | `/api/v1` after the GUI moves to `/api/gui`. If external scripts are known to use today's `/api/v1`, start the public API at `/api/v2` instead and keep `/api/v1` as the deprecated GUI alias for longer |
| D2 | Synchronous calls | `POST /api/v1/jobs?wait=<s>` rather than a separate synchronous endpoint |
| D3 | Units of inputs and outputs | The units declared in the model. An option `"unitSystem": "SI"` that converts every value to SI could be added later without breaking v1 |
| D4 | Arrays | Grouped (`index`, `value`), as in section 4.4 |
| D5 | Stored models on the server (`POST /api/v1/models`, then jobs by model id) | Not in v1: requests stay self-contained. Library ids cover the "known model" case. Can be added later |
| D6 | Authentication | Bearer token in headless mode on a network address; reverse proxy otherwise |
| D7 | Per-request CoolProp settings | Not exposed in v1, until they can be set per job without affecting other jobs (F4) |
| D8 | `--serve` on a network address without a token | Refuse to start unless `--auth=proxy` is given (section 4.9); the GUI mode keeps its warning |
| D9 | Base image of the runtime stage | `ubuntu:24.04` first (same as the build, tools for debugging); distroless C++ base once the image is stable |
| D10 | Where the image is published | GitHub Container Registry of the `CoolProp` organisation, next to the releases; Docker Hub only if users ask |
