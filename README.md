# CoolSolve

CoolSolve is a parser, structural analyzer, and equation evaluator for the EES (Engineering Equation Solver) language, designed to be the foundation of an open-source equation solver that integrates with CoolProp for thermodynamic property calculations.

## Features

- **EES Parser**: Parses EES source code (.eescode files) into an Abstract Syntax Tree
  - Variables (scalars and arrays like `P[1]`)
  - Equations with operators (`+`, `-`, `*`, `/`, `^`)
  - **Procedural Statements**: Support for `=` and `:=` in procedural blocks
  - **Functions and Procedures**: `FUNCTION` and `PROCEDURE` blocks with local scoping
  - **Procedure Calls**: `CALL name(inputs : outputs)` syntax
  - **Conditional Logic**: `IF-THEN-ELSE` support within procedural blocks
  - Comments (`"..."`, `{...}`, `//...`)
  - Directives (`$ifnot`, `$endif`, etc.)
  - Units annotations (`"[Pa]"`, `"[kJ/kg]"`)
  - Function calls with named arguments (`enthalpy(R134a, T=300, P=1E5)`)

- **Structural Analysis**: Decomposes the equation system into solvable blocks
  - Hopcroft-Karp algorithm for variable-equation matching
  - Tarjan's algorithm for Strongly Connected Components (SCCs)
  - Block decomposition (identifies algebraic loops)

- **Automatic Differentiation**: Forward-mode AD for efficient Jacobian computation
  - Full support for arithmetic operators (`+`, `-`, `*`, `/`, `^`)
  - Transcendental functions (`sin`, `cos`, `exp`, `log`, `sqrt`, etc.)
  - Exact analytical derivatives (no finite differences)
  - Prepares for Newton-based numerical solvers

- **Equation Evaluator**: Evaluates residuals and Jacobians for each block
  - Block-level evaluation with external variable support
  - System-level orchestration for sequential block solving
  - **CoolProp integration** for thermodynamic property calculations
  - Automatic temperature conversion (Celsius ↔ Kelvin)

- **Output Formats**:
  - JSON (for automated testing and integration)
  - EES-compatible residuals format
  - LaTeX equations
  - Evaluator report (block summary and state)
  - **Solution/Initials files**: Support for loading initial values and comparing solutions

- **Debug Mode**: Creates a comprehensive output folder with all analysis information

## Performance and Roadmap

For details on current performance optimizations and the future roadmap for solver integration, see the [Performance Improvement and Solver Integration Plan](performance_plan.md).

## File Formats

CoolSolve uses several file formats for input and verification:

- **.eescode**: The EES source code to be parsed and solved.
- **.initials**: Initial values for variables, used to seed the solver or evaluator. Format: `variable=value` (one per line).

## Building

### Prerequisites

- C++17 compatible compiler (GCC 7+, Clang 5+)
- CMake 3.14 or later
- Git (for fetching dependencies)

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

### Build Type: Release vs Debug

**Important**: CoolSolve defaults to **Release** mode for optimal performance. CoolProp property calculations are computationally intensive, and Debug mode can be **10-50x slower** than Release.

```bash
# Default: Release mode (recommended)
cmake ..

# Explicitly set Release mode
cmake -DCMAKE_BUILD_TYPE=Release ..

# Debug mode (for development only - very slow!)
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

If you experience slow solve times, verify you're building in Release mode:
```bash
# Check current build type
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
# Parse an EES file and output JSON analysis
./coolsolve input.eescode

# Output in EES residuals format
./coolsolve -f residuals input.eescode

# Output LaTeX equations
./coolsolve -f latex input.eescode

# Save output to a file
./coolsolve -o output.json input.eescode
```

### Debug Mode

Debug mode creates a folder containing all analysis information, useful for understanding the equation structure and debugging:

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
| `report.md` | Model statistics and block summary |
| `variables.md` | Variable mapping table (Markdown) |
| `variables.csv` | Variable mapping (CSV for external tools) |
| `equations.md` | Equations grouped by solution block |
| `analysis.json` | Full JSON analysis data |
| `residuals.txt` | EES-compatible residuals format |
| `equations.tex` | LaTeX formatted equations |
| `incidence.md` | Variable-equation incidence matrix |
| `evaluator.md` | Evaluator structure and block evaluation tests |
| `original.eescode` | Copy of the original input |

### Compare with EES

To validate the structural analysis against EES output:

```bash
./coolsolve -c reference.residuals input.eescode
```

### Command Line Options

| Option | Description |
|--------|-------------|
| `-o, --output <file>` | Output file (default: stdout) |
| `-f, --format <fmt>` | Output format: `json`, `residuals`, `latex` (default: json) |
| `-c, --compare <file>` | Compare with EES .residuals file |
| `-d, --debug [dir]` | Create debug output folder |
| `-v, --verbose` | Verbose progress output |
| `--no-sol` | Disable generation of .sol file |
| `-g, --guess` | Update .initials file with solution on success |
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

# Run all unit tests
./coolsolve_tests

# Run comprehensive example file tests
./coolsolve_tests "[examples-comprehensive]"
```

Or using CTest:
```bash
cd build
ctest --output-on-failure
```

The comprehensive test runs all `.eescode` files in the `examples/` folder and validates solutions against known expected values (with 1% tolerance). A detailed report is written to `examples/test_examples.md`.

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
│   ├── autodiff_node.h         # Forward-mode AD types and operations
│   └── evaluator.h             # Block and system evaluators
├── src/
│   ├── parser.cpp              # EES parser implementation
│   ├── ir.cpp                  # IR building and LaTeX generation
│   ├── structural_analysis.cpp # Matching and SCC algorithms
│   ├── autodiff_node.cpp       # AD function implementations
│   └── evaluator.cpp           # Evaluator implementations
├── tests/
│   ├── test_parser.cpp         # Parser/IR unit tests (Catch2)
│   ├── test_evaluator.cpp      # AD/Evaluator unit tests (Catch2)
│   └── test_examples.cpp       # Integration tests with example files
└── examples/                   # Example .eescode files for testing
```

## How It Works

### 1. Parsing

The parser reads EES source code and builds an Abstract Syntax Tree (AST). It handles EES-specific syntax including:
- Case-insensitive keywords
- Multiple comment styles
- Thermodynamic function calls with named parameters
- Array subscript notation

### 2. IR Building

The AST is transformed into an Intermediate Representation (IR) that:
- Extracts all variables from equations
- Builds the incidence matrix (which variables appear in which equations)
- Identifies thermodynamic function calls (first argument is fluid name, not a variable)

### 3. Structural Analysis

The equation system is analyzed using graph algorithms:

1. **Matching**: The Hopcroft-Karp algorithm finds a maximum bipartite matching, assigning each equation to one "output" variable it will determine.

2. **SCC Detection**: Tarjan's algorithm finds Strongly Connected Components in the dependency graph. Each SCC becomes a "block" that must be solved simultaneously.

3. **Block Ordering**: Blocks are topologically sorted so they can be solved in sequence.

### 4. Automatic Differentiation

Forward-mode automatic differentiation computes exact derivatives alongside values:

1. **ADValue**: A dual number storing both `value` and `gradient` (partial derivatives w.r.t. all block variables).

2. **Propagation Rules**: Each operation propagates gradients using the chain rule:
   - Addition: `(x + y).grad = x.grad + y.grad`
   - Multiplication: `(x * y).grad = x.grad * y + y.grad * x`
   - Power: `(x^n).grad = n * x^(n-1) * x.grad`
   - Functions: `sin(x).grad = cos(x) * x.grad`

3. **Jacobian Construction**: For each equation residual F_i, the gradient gives row i of the Jacobian matrix.

### 5. Evaluation

The evaluator system provides numerical computation:

1. **ExpressionEvaluator**: Traverses the AST and computes ADValues for any expression.

2. **BlockEvaluator**: Evaluates a block of equations, returning residuals F(x) and Jacobian J.

3. **SystemEvaluator**: Orchestrates evaluation across all blocks with proper variable handling.

### 6. Output

The analysis results can be exported in various formats for:
- Integration with numerical solvers (JSON)
- Documentation (LaTeX)

## Example

Given this EES code:
```
T_in = 25            // Temperature in Celsius
P = 101325           // Pressure in Pa
h = enthalpy(Water, T=T_in, P=P)
s = entropy(Water, T=T_in, P=P)
```

CoolSolve will:
1. Parse 4 equations
2. Identify 4 variables: `T_in`, `P`, `h`, `s`
3. Create 4 blocks (all explicit, no algebraic loops)
4. Determine solution order: T_in → P → h, s
5. Evaluate thermodynamic properties using CoolProp (with automatic °C→K conversion)

## CoolProp Integration

CoolSolve integrates with CoolProp for thermodynamic property calculations:

- **Supported functions**: `enthalpy`, `entropy`, `density`, `pressure`, `temperature`, `quality`, `cp`, `cv`
- **Supported input pairs**: `T,P`, `T,x`, `P,h`, `P,s`, `H,P`, and others
- **Temperature units**: Automatically converts Celsius (EES convention) to Kelvin (CoolProp)
- **Derivatives**: Computed using finite differences from CoolProp's `PropsSI` function
- **Fluids**: All CoolProp fluids are supported (Water, R134a, Nitrogen, etc.)

### Default Units

CoolSolve assumes the following default units:
- Temperature: **°C** (automatically converted to K for CoolProp calls)
- Pressure: **Pa**
- Energy: **J**
- Mass: **kg**
- Force: **N**
- Length: **m**

Built-in constants like `g#` use these SI units (e.g., `g# = 9.80665 m/s^2`).

Example usage in EES code:
```
T_ev = -10           // Celsius
P_ev = pressure(R134a, T=T_ev, x=1)   // Saturation pressure
h_1 = enthalpy(R134a, P=P_ev, T=T_1)  // Enthalpy at state 1
```

### Functions and Procedures

CoolSolve supports user-defined functions and procedures for modular code:

```ees
PROCEDURE single_phase_HX(cf$, hf$, t_cf_su : Q_dot)
    "Procedural block with local scope"
    cp_cf = specheat(cf$, T=t_cf_su, P=101325)
    Q_dot := alpha * cp_cf * (t_hf_su - t_cf_su)
END

CALL single_phase_HX('Air_ha', 'Water', 20 : Q_total)
```

- **Local Scope**: Variables inside functions and procedures are local and do not interfere with the main program.
- **Assignment**: Use `:=` for assignment inside procedural blocks (though `=` is also supported for compatibility).
- **Control Flow**: `IF-THEN-ELSE` statements are supported within these blocks.
- **Automatic Differentiation**: CoolSolve automatically propagates derivatives through procedural calls, ensuring accurate Jacobians.

### Reference State Warning

**Important**: EES and CoolProp use different reference states for thermodynamic properties:

| Property | EES (IIR Reference) | CoolProp (ASHRAE Reference) |
|----------|---------------------|----------------------------|
| Enthalpy | h = 200 kJ/kg at 0°C sat. liquid | h = 0 at -40°C |
| Entropy  | s = 1 kJ/kg/K at 0°C sat. liquid | s varies |

This means:
- Enthalpy values from EES may differ from CoolProp by ~100-200 kJ/kg
- Entropy values from EES may differ from CoolProp by ~0.8-1 kJ/kg/K
- **Differences** (e.g., h2 - h1) are the same in both reference states
- When loading initial values from `.initials` files, some CoolProp function calls may fail if the enthalpy/entropy is outside the expected range for the given pressure

For new equation systems, use CoolProp-computed values for initial guesses to ensure consistency.

## Future Work

The next steps in the implementation plan include:
- **Numerical Solver**: Integration with KINSOL or Ceres for Newton iteration on algebraic loops
- **Analytical Derivatives**: Use CoolProp's `AbstractState::first_partial_deriv` for exact derivatives
- **Profiling and optimization**: Improve the solving time
