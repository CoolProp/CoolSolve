<p align="center">
  <img src="docs/logo_coolsolve.png" alt="CoolSolve" width="200">
</p>

# CoolSolve

CoolSolve is a parser, structural analyzer, and equation solver using a language compatible with EES (Engineering Equation Solver), designed as an open-source alternative that integrates with CoolProp for thermodynamic property calculations.

![CoolSolve GUI](docs/coolsolve_gui.png)

## Try CoolSolve

### Windows Desktop Application (v0.3)

Download the Windows installer: **[CoolSolve v0.3 for Windows](https://github.com/CoolProp/CoolSolve/releases/download/v0.3/CoolSolve_v0.3_Installer.exe)**

See [Version History](docs/versions.md) for all releases, changelogs, and
the CoolProp version bundled in each release.

### Online Demo

Try CoolSolve in your browser: **[https://coolsolve.squoilin.eu/](https://coolsolve.squoilin.eu/)**

> **⚠ Disclaimer**: The online demo runs on a private server and corresponds to the latest development version, which may not always be stable. Under high load, the server may respond slowly or be temporarily unavailable.

## Features

- **Parser**: Parses CoolSolve source code (`.eescode` files) into an AST — variables, equations, procedural blocks (`FUNCTION`, `PROCEDURE`, `CALL`), control flow (`IF-THEN-ELSE`, `DUPLICATE`, `REPEAT-UNTIL`), directives, array subscripts, and named function arguments.

- **Structural Analysis**: Decomposes the equation system into solvable blocks using the Hopcroft–Karp matching algorithm and Tarjan's SCC algorithm.

- **Automatic Differentiation**: Forward-mode AD with exact analytical derivatives for all arithmetic, transcendental, and CoolProp functions — no finite differences.

- **Configurable Solver Pipeline**: Multiple algorithms with automatic fallback:
  - **Newton + Line Search** (default workhorse)
  - **Trust-Region Dogleg** with optional hybrd-style Broyden reuse
  - **Levenberg-Marquardt** with geodesic acceleration
  - **Multi-dimensional Bisection** (derivative-free, for singular Jacobians)
  - **Homotopy Continuation** (for distant starting points)
  - **KINSOL** (opt-in: Dennis-Schnabel line search, Picard, Anderson-accelerated fixed-point)
  - **Structural Tearing** (equation tearing for large blocks)
  - **Symbolic Block Reduction** (explicit extraction, CoolProp call inversion, equation substitution)
  - **Multi-Start** fallback (CoolProp-consistent starting points for thermo blocks, scale factors for algebraic blocks; parallel execution available; engageable always, never, or only on Deep Search)

- **CoolProp Integration**: Thermodynamic property calculations via the low-level `AbstractState` API (2–5× faster than `PropsSI`), with analytical derivatives and thread-local caching.

- **Lookup Tables**: 1D/2D interpolation, cell access, and aggregate functions from external CSV files. A GUI panel lets you create and edit tables in-browser.

- **Equation-Based Dynamic Solving (`INTEGRAL`)**: Solve initial-value DAE models in the EES integral form (`y = y0 + INTEGRAL(dydt, t, t0, tf)`). Supports coupled ODEs, algebraic variables, RK4/RK45/Euler integrators, Richardson extrapolation, and an **Integral** GUI tab. See [Language Reference §12](docs/language_reference.md#12-equation-based-integration-dynamicdae-solving) and [docs/integral_table.md](docs/integral_table.md).

  - **Solution Verification**: Post-solve check that re-evaluates every equation (LHS ≈ RHS) to catch bugs the solver's convergence check may miss.

  - **Try Harder**: One-click recovery for failed solves. After a failed run the GUI's *Solve* button morphs into *Try Harder*; clicking it re-runs the model with the full Deep Search pipeline, tearing and symbolic reduction forced on, and the configured multi-start policy. Editing the model, the initials, or the configuration restores the normal *Solve* button.

- **Output Formats**: JSON, LaTeX, `.sol` files, CSV trajectory files, and a comprehensive **Debug Mode** (`-d`) folder.

- **GUI & REST API**: Embedded single-page app (React/TypeScript) with a code editor, variable table, parametric studies, thermodynamic diagrams, and ZIP bundle round-trip. See [GUI & REST API](docs/gui.md).

## Documentation

| Document | Description |
|----------|-------------|
| [Language Reference](docs/language_reference.md) | CoolSolve language syntax, built-in functions, and configuration keys |
| [Dynamic Solving](docs/integral_table.md) | Equation-based `INTEGRAL` / `$IntegralTable` — algorithms, architecture, and limitations |
| [Debugging Models](docs/debugging_models.md) | Diagnosing and fixing solver failures |
| [Solver Roadmap](docs/solver_roadmap.md) | Performance roadmap, algorithm details, and implementation status |
| [Symbolic Block Reduction](docs/symbolic_redecomposition.md) | Symbolic block reduction algorithm |
| [GUI & REST API](docs/gui.md) | Web interface, REST API, and parametric studies |
| [Deployment Guide](docs/deployment_ubuntu_apache.md) | Deploying CoolSolve on Ubuntu with Apache |
| [Contributing Guide](docs/contributing.md) | Workflow and integration checklist for new features and bug fixes |
| [Version History](docs/versions.md) | Release notes, download links, and CoolProp version per release |

## Performance and Roadmap

For details on current performance optimizations and the future roadmap for solver integration, see the [Solver Roadmap](docs/solver_roadmap.md).

## File Formats

CoolSolve uses several file formats for input and verification:

- **.eescode**: The EES source code to be parsed and solved.
- **.initials**: Initial values for variables, used to seed the solver or evaluator. Format: `variable=value` (one per line).
- **coolsolve.conf**: Optional solver configuration. Place in the **same folder** as your .eescode file (not in subfolders). Format: `key = value` per line; lines starting with `#` are comments. Only the options you set override the defaults (see `include/coolsolve/solver.h` for `SolverOptions`). An example with all keys and comments is in `examples/coolsolve.conf`. In debug mode (`-d`), this file is copied into the debug folder. Key options:
  - **Pipeline options**:
    - `solverPipeline`: Comma-separated list of solvers to try (e.g. `Newton, LM, TrustRegion, BisectionND, Homotopy, Partitioned`). Available solvers: `Newton`, `TrustRegion`, `LM` (or `LevenbergMarquardt`), `BisectionND`, `Homotopy`, `Partitioned`, `Kinsol`.
    - `pipelineMode`: `sequential` (default) or `parallel` (first-to-converge wins)
    - `deepSearchPipeline` / `deepSearchPipelineMode`: pipeline used when the GUI **Try Harder** button is clicked after a failed solve (defaults to the full sequential chain). Deep Search also forces tearing and symbolic reduction on.
    - `enableTearing`: When `true`, use structural tearing for blocks of size ≥ `tearingMinBlockSize`.
    - `enableSymbolicReduction`: When `true`, pre-process blocks to reduce their size via explicit extraction, CoolProp call inversion, and equation substitution, with automatic re-decomposition of the reduced block.
    - `bisectionNDMaxBlockSize`: Maximum block size for BisectionND (default: `8`).
    - `bisectionNDIterFactor`: Multiplier for BisectionND iteration budget (default: `1.0`).
    - `lsNonMonotoneMemory`: Number of recent merit values kept for non-monotone acceptance (default: `10`). Set to `1` for classic monotone line search.
    - `multiStartMode`: When to retry a failed multi-variable block from alternative starting points — `always`, `deepsearch` (default; only when running a Deep Search), or `never`. The legacy `multiStartEnabled = true/false` is still accepted and maps to `always`/`never`.
    - `multiStartMaxRestarts`: Number of alternative starting points to try on a failed block (default: `4`).
    - `multiStartNumCores`: Threads for concurrent multi-start candidates (default: `4`; `N`>1 or `0`=auto runs candidates in parallel, first-to-converge wins).
   - **Integration options** (equation-based `INTEGRAL` models — all inert by default; see [Language Reference §12](docs/language_reference.md#12-equation-based-integration-dynamicdae-solving)):
     - `integralMethod`: Integrator — `RK4` (default), `RK45` (adaptive), `EulerExplicit`, or `EulerImplicit`.
     - `integralFixedStep`: Fixed step size (`0` ⇒ derive from `integralMaxSteps` for fixed methods, or adapt for RK45).
     - `integralMaxSteps`: Upper bound on the number of integration steps (default: `1000`).
     - `integralRelTol` / `integralAbsTol`: RK45 error control (defaults `1e-6` / `1e-9`).
     - `integralMinStep` / `integralMaxStep`: Step-size bounds (`0` = auto).
     - `integralRichardson`: Richardson extrapolation on fixed-step methods (default: `false`).
     - `integralOutputInterval`: Default `$IntegralTable` row interval when `:n` is omitted (`0` = every step).
   - **CoolProp Integration options**:
    - `coolpropBackend`: CoolProp backend string (default: `HEOS`). Options: `HEOS`, `INCOMP`, `TTSE&HEOS`, `BICUBIC&HEOS`.
    - `coolpropUseAbstractState`: Use the low-level `AbstractState` API instead of `PropsSI` (default: `true`). Provides 2–5× speedup by caching fluid objects and avoiding string parsing.
    - `coolpropEnableAnalyticalDerivatives`: Use `first_partial_deriv()` for exact gradients with a forward-FD consistency check (default: `true`). Falls back to finite differences near phase boundaries where analytical derivatives are inaccurate.
    - `coolpropCacheEnabled`: Enable thread-local `AbstractState` caching (default: `true`).
    - `coolpropEnableSuperancillaries`: Enable CoolProp superancillary equations (default: `true`).

## Building

### Prerequisites

- C++17 compatible compiler (GCC 7+, Clang 5+, or MSVC 2019+)
- CMake 3.14 or later
- Git (for fetching dependencies)
- Python 3 (required by CoolProp's CMake scripts at configure time)
- Node.js 18+ and npm (required to build the React GUI frontend)

### How Dependencies Are Handled

CoolSolve is a **standalone library** that automatically downloads all its dependencies using CMake's FetchContent mechanism:

| Dependency | Purpose |
|------------|---------|
| **CoolProp** | Thermodynamic property calculations |
| **cpp-peglib** | PEG parser generator |
| **nlohmann/json** | JSON serialization |
| **Eigen** | Linear algebra |
| **Catch2** | Unit testing |

Dependencies are cached in `.fetchcontent_cache/` at the project root. This cache **persists across build folder deletions**, so you won't need to re-download dependencies if you clean and rebuild. The first `cmake` run will take several minutes to fetch CoolProp and its submodules; subsequent runs take only a few seconds.

### Build Steps

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

This will build:
- `coolsolve` - The main executable
- `coolsolve_tests` - The test suite

### Building on Windows

#### Prerequisites
- **[Visual Studio 2022](https://visualstudio.microsoft.com/)** with the **"Desktop development with C++"** workload (includes MSVC, CMake, and Windows SDK)
- **[Python 3](https://www.python.org/downloads/)** (or Anaconda) — required by CoolProp's CMake scripts; must be on `PATH`
- **[Node.js](https://nodejs.org/)** (v18+) — required to build the embedded React GUI
- **[NSIS](https://nsis.sourceforge.net/Download)** — to generate the single-file installer (optional, installer only)

#### Building the GUI frontend
The React frontend must be built before CMake can embed it into the binary:
```bat
cd gui
npm ci
npm run build
cd ..
```
This produces `gui/dist/` which CMake will embed into the `.exe`.

#### CMake build (Developer PowerShell for VS 2022)
```bat
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64 -DCOOLSOLVE_BUILD_GUI=ON -DCOOLSOLVE_BUILD_TESTS=OFF
cmake --build . --config Release --target coolsolve
```

The compiled binary is at `build\Release\coolsolve.exe` and is fully self-contained (no external runtime or DLLs needed).

#### Generating the installer
After a successful build:
```bat
makensis coolsolve.nsi
```
This produces `CoolSolve_Installer.exe`, a single-file installer that:
- Installs to `Program Files\CoolSolve`
- Adds a Start Menu shortcut
- Associates `.eescode` files with CoolSolve (double-click to open)
- Includes an uninstaller

#### One-click build script
A convenience script `build_installer.bat` in the project root automates all steps above. Edit the paths at the top of the file if your Visual Studio, Anaconda, or NSIS installations are in non-default locations, then right-click → **Run as Administrator**.

### Build Type: Release vs Debug

**Important**: CoolSolve defaults to **Release** mode for optimal performance. CoolProp property calculations are computationally intensive, and Debug mode can be **10-50x slower** than Release. **Do not switch to Debug** for routine work, benchmarks, or robustness reports — only when you genuinely need a debugger on the C++ code. **When you are done debugging, always switch back to Release and rebuild:**

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j$(nproc)    # or: make -j$(nproc)
```

```bash
# Default: Release mode (recommended)
cmake ..

# Explicitly set Release mode
cmake -DCMAKE_BUILD_TYPE=Release ..

# Debug mode (C++ debugger only — revert to Release afterwards!)
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

If you experience slow solve times, verify you're building in Release mode:
```bash
# Check current build type (must show Release)
grep CMAKE_BUILD_TYPE CMakeCache.txt
```

### Updating CoolProp Version

By default, CoolSolve fetches the latest `master` branch of CoolProp. To use a specific version:

```bash
# Use a specific tag or commit
cmake -DCOOLSOLVE_COOLPROP_TAG=v6.6.0 ..

# Or a specific commit hash
cmake -DCOOLSOLVE_COOLPROP_TAG=abc123def ..
```

### Cleaning the Build

```bash
# Clean build artifacts only (keeps dependency cache)
rm -rf build && mkdir build && cd build && cmake .. && make -j$(nproc)

# Full clean including dependency cache (will re-download everything)
rm -rf build .fetchcontent_cache
```

## Usage

### Basic Usage

```bash
# Parse a .eescode file and output JSON analysis
./coolsolve input.eescode

# Output LaTeX equations
./coolsolve -f latex input.eescode

# Save output to a file
./coolsolve -o output.json input.eescode
```

### Debug Mode

Debug mode creates a folder containing all analysis information, useful for understanding the equation structure and debugging. For a complete guide on diagnosing and fixing solver failures (including the simplified-model workflow), see [Debugging Models](docs/debugging_models.md).

```bash
# Create debug output in <input>_coolsolve/ folder
./coolsolve -d input.eescode

# Specify custom debug output directory
./coolsolve -d my_debug_folder input.eescode
```

The debug folder contains:
| File | Description |
|------|-------------|
| `README.md` | Index of all generated files |
| `coolsolve.conf` | Copy of solver config from source folder (if present) |
| `report.md` | Model statistics and block summary |
| `variables.md` | Variable mapping table (Markdown) |
| `variables.csv` | Variable mapping (CSV for external tools) |
| `equations.md` | Equations grouped by solution block |
| `analysis.json` | Full JSON analysis data |
| `residuals.txt` | CoolSolve residuals report |
| `equations.tex` | LaTeX formatted equations |
| `incidence.md` | Variable-equation incidence matrix |
| `evaluator.md` | Evaluator structure and block evaluation tests |
| `symbolic_reduction.md` | Symbolic reduction debug report (when `enableSymbolicReduction` is on) |
| `integral.md` | Integral solver report (dynamic/`INTEGRAL` models) |
| `integral_table.csv` | Full trajectory CSV (dynamic/`INTEGRAL` models) |
| `solution_check.md` | Post-solve equation verification (LHS vs RHS for every equation) |
| `original.eescode` | Copy of the original input |

### Command Line Options

| Option | Description |
|--------|--------------|
| `-o, --output <file>` | Output file (default: stdout) |
| `-f, --format <fmt>` | Output format: `json`, `latex` (default: json) |
| `-d, --debug [dir]` | Create debug output folder |
| `-g, --guess` | Update .initials file with solution on success |
| `--no-sol` | Disable generation of .sol file |
| `--no-superancillary` | Disable CoolProp superancillary functions (faster VLE solving) |
| `-h, --help` | Show help message |

### Performance Options

#### Superancillary Functions

CoolProp's superancillary functions provide high-accuracy VLE (vapor-liquid equilibrium) calculations, especially near the critical point. However, they add computational overhead during solving.

```bash
# Default: superancillaries enabled 
./coolsolve model.eescode

# superancillaries disabled
./coolsolve --no-superancillary model.eescode

# Alternative: use environment variable
COOLPROP_ENABLE_SUPERANCILLARIES=false ./coolsolve model.eescode
```

## Running Tests

```bash
cd build

# Run all unit tests (parser, evaluator, Newton, tearing, config, etc.)
./coolsolve_tests

# Run comprehensive example file tests (solves all 40 .eescode examples,
# checks expected values, and runs solution verification on each)
./coolsolve_tests "[examples-comprehensive]"

# Run solver robustness tests (stress-tests every example with multiple
# solver configurations: Newton-only, LM-only, tearing, symbolic reduction, etc.)
# Requires a Release build — Debug timings are ~10x slower and must not be committed as reports.
./coolsolve_tests "[solver-robustness]"
```

Or using CTest:
```bash
cd build
ctest --output-on-failure
```

The comprehensive test runs all `.eescode` files in the `examples/` folder, validates solutions against known expected values (with 1% tolerance), and performs **solution verification** — independently re-evaluating every equation to confirm LHS ≈ RHS. A detailed report is written to `examples/test_examples.md`.

The robustness test runs every solvable example under multiple solver configurations and reports a summary table of successes, failures, and iteration counts per configuration. A detailed report is written to `examples/solver_robustness_report.md`.

## GUI (Web Interface)

The CoolSolve GUI is a React/TypeScript single-page app served by an embedded HTTP server inside the `coolsolve` binary. The GUI is **optional**; the CLI continues to work as described above.

### GUI Prerequisites

- **Node.js + npm**: Required only to build or develop the frontend (no Node.js is needed at runtime once the binary is built).

### Run the GUI in development mode (hot reload)

Use this when working on the GUI itself:

```bash
# Terminal 1: build and run the backend with the embedded HTTP server
mkdir -p build
cd build
cmake -DCOOLSOLVE_BUILD_GUI=ON ..
make -j$(nproc) coolsolve
./coolsolve --gui --no-browser

# Terminal 2: run the frontend dev server with hot reload
cd gui
npm install
npm run dev
# Then open http://localhost:5173 in your browser
```

In this mode, Vite serves the frontend on port `5173` and proxies API calls to the CoolSolve server (default port `8550`).

### Run the GUI from the compiled binary (no dev server)

To use the GUI as an integrated part of the `coolsolve` binary (no Node.js or Vite needed at runtime):

```bash
# 1) Build the frontend once
cd gui
npm install
npm run build        # produces gui/dist/

# 2) Build or rebuild the backend so it embeds gui/dist/ into the binary
cd ../build
cmake -DCOOLSOLVE_BUILD_GUI=ON ..
make -j$(nproc) coolsolve

# 3) Run the GUI
./coolsolve --gui     # opens the default browser, serving on http://localhost:8550
```

For more detailed information on the GUI architecture and advanced workflows (ZIP bundles, online deployment, thermodynamic diagrams, etc.), see `docs/gui.md`.

## Project Structure

```
CoolSolve/
├── CMakeLists.txt              # Build configuration (handles all dependencies)
├── README.md                   # This file
├── main.cpp                    # CLI entry point
├── .fetchcontent_cache/        # Cached dependencies (auto-generated, git-ignored)
├── include/coolsolve/
│   ├── ast.h                   # Abstract Syntax Tree definitions
│   ├── parser.h                # Parser interface
│   ├── ir.h                    # Intermediate Representation
│   ├── structural_analysis.h   # Analysis algorithms
│   ├── symbolic_reduction.h    # Symbolic block reduction analysis
│   ├── autodiff_node.h         # Forward-mode AD types and operations
│   ├── evaluator.h             # Block and system evaluators
│   ├── solution_checker.h      # Post-solve solution verification
│   └── solver.h                # Solver pipeline, all SolverOptions declarations
├── include/coolsolve/integral/
│   └── *.h                     # Equation-based dynamic solver (IntegralProblem, IntegralTable, Integrator, IntegralSolver)
├── src/
    ├── parser.cpp              # CoolSolve language parser
│   ├── ir.cpp                  # IR building and LaTeX generation
│   ├── structural_analysis.cpp # Matching and SCC algorithms
│   ├── autodiff_node.cpp       # AD function implementations
│   ├── evaluator.cpp           # Evaluator implementations
│   ├── integral/               # Equation-based dynamic solver (INTEGRAL / $IntegralTable)
│   │   ├── integrator_euler_explicit.cpp  # Forward Euler (fixed step)
│   │   ├── integrator_euler_implicit.cpp  # Backward Euler, A-stable (fixed step)
│   │   ├── integrator_rk4.cpp             # Classic RK4 (fixed step, default)
│   │   ├── integrator_rk45.cpp            # Dormand-Prince RK45 (adaptive)
│   │   ├── richardson.cpp                 # Richardson extrapolation wrapper
│   │   ├── integral_table.cpp             # Trajectory storage + interpolation + CSV
│   │   ├── integral_extraction.cpp        # IR → IntegralProblem classification
│   │   └── integral_solver.cpp            # Time-march loop, reuses algebraic Solver
│   ├── solver.cpp              # Pipeline orchestrator, Newton, tearing, config loading
│   ├── solver_bisection_nd.cpp # Multi-dimensional bisection solver
│   ├── solver_homotopy.cpp     # Homotopy continuation solver
    │   ├── solver_lm.cpp           # Levenberg-Marquardt solver
    │   ├── solver_kinsol.cpp       # KINSOL (SUNDIALS-style) solver: line search / Picard / Anderson
    │   ├── solver_newton.cpp       # Newton + line search solver
│   ├── solver_symbolic.cpp     # Symbolic block reduction preprocessing
│   ├── solver_trust_region.cpp # Trust-region dogleg solver
│   └── solution_checker.cpp    # Post-solve equation-by-equation verification
├── tests/
│   ├── test_parser.cpp         # Parser/IR unit tests (Catch2)
│   ├── test_evaluator.cpp      # AD/Evaluator unit tests (Catch2)
│   ├── test_newton.cpp         # Newton solver unit tests
│   ├── test_solver_pipeline.cpp # Pipeline, LM, config tests
│   ├── test_solver_integration.cpp # Full-system solver integration tests
│   ├── test_new_solvers.cpp    # BisectionND, Homotopy, advantage tests
│   ├── test_solver_robustness.cpp  # Robustness tests (stiff/near-singular blocks)
│   ├── test_symbolic_reduction.cpp # Symbolic block reduction tests
│   ├── test_tearing.cpp        # Structural tearing unit tests
│   ├── test_config.cpp         # coolsolve.conf loading tests
│   ├── test_fluids.cpp         # CoolProp fluid property tests
│   ├── test_examples.cpp       # Integration tests with example files
│   └── test_kinsol.cpp         # KINSOL solver unit tests (3 modes + config)
└── examples/                   # Example .eescode files for testing
```

## How It Works

CoolSolve follows a standard equation-solving pipeline:

1. **Parsing** — Source code (`.eescode`) is parsed into an AST and then transformed into an Intermediate Representation (IR) that builds the variable-equation incidence matrix.

2. **Structural Analysis** — The Hopcroft–Karp algorithm matches equations to output variables, Tarjan's algorithm finds strongly connected components (SCCs), and blocks are topologically sorted.

3. **Automatic Differentiation** — Forward-mode AD computes exact derivatives alongside values, building the Jacobian matrix analytically.

4. **Evaluation** — The evaluator computes residuals F(x) and Jacobian J for each block.

5. **Solving** — Each block is solved by the configurable pipeline (`Newton → TrustRegion → LM → BisectionND → Homotopy → Partitioned`), with optional tearing and symbolic reduction. For detailed algorithm descriptions, see the [Solver Roadmap](docs/solver_roadmap.md).

6. **Dynamic Solving** — When the model contains `INTEGRAL(...)` calls, CoolSolve time-marches the state variables and re-solves the algebraic subsystem at each step. See [Dynamic Solving](docs/integral_table.md).

For implementation details, see the [Language Reference](docs/language_reference.md) and [Contributing Guide](docs/contributing.md).

## Example

```ees
T_in = 25            // Temperature in Celsius
P = 101325           // Pressure in Pa
h = enthalpy('Water', T=T_in, P=P)
s = entropy('Water', T=T_in, P=P)
```

CoolSolve parses these 4 equations, identifies the variable-equation structure, and evaluates the thermodynamic properties using CoolProp. Temperatures are written in **°C** in the model — CoolSolve's unit system is fixed (it cannot be changed through `coolsolve.conf`) and the conversion `°C → K` required by CoolProp's SI interface is performed automatically behind the scenes.

For the full language syntax, built-in functions, and CoolProp integration details, see the [Language Reference](docs/language_reference.md).

## Future Work

Future features include a complete library of models embedded into coolsolve and improved plotting capabilities. 

Planned improvements for the solver include a stiff ODE integrator (BDF), pseudo-arclength continuation, and an improved plotting interface. See [docs/solver_roadmap.md](docs/solver_roadmap.md) for the full prioritized roadmap.

## Acknowledgements & License

CoolSolve is released under the [MIT License](LICENSE).

**Main contributor**: [Sylvain Quoilin](https://www.uliege.be/cms/c_9054334/en/directory?uid=U203754) — [ISES Research Group](https://www.ises.uliege.be/), Université de Liège.

This work is part of the broader development of [CoolProp](https://github.com/CoolProp/CoolProp), an open-source thermophysical property library.


