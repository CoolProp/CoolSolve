# CoolSolve Performance Improvement and Solver Integration Plan

This document outlines the plan to address the performance regression identified in commit `f4063bb4` and the roadmap for integrating alternative solvers and monitoring tools.

## 1. Monitoring and Profiling Tools

To ensure long-term performance and code quality, the following tools will be implemented:

### 1.1. Solving Time Tracker
- **Tool**: A utility to track and log the execution time of individual blocks and the total model solution time.
- **Integration**: Results will be included in the `report.md` generated in debug mode and optionally logged to a performance history file.
- **Purpose**: Detect performance regressions in future changes.

### 1.2. Code Complexity Monitor
- **Tool**: A script integrated into the test suite (e.g., using `cloc` or a custom counter) to track the total lines of code (LOC) for the `coolsolve` directory.
- **Purpose**: Identify significant code duplication or unnecessary complexity growth.

### 1.3. Detailed Profiling
- **Action**: Use system profilers (like `gprof`, `perf`, or Valgrind's `callgrind`) to identify remaining bottlenecks in the parser, structural analyzer, and solver.

### 1.4. Profiling Results (refrigeration2.eescode)

**Test Case**: R22 heat pump with 39 equations, all explicit blocks

**Timing Breakdown** (total: ~35 seconds):

| Component | Time | Notes |
|-----------|------|-------|
| Parse | 0.3 ms | Negligible |
| IR Build | 0.5 ms | Negligible |
| **Variable Inference** | **32,600 ms** | **95% of total time!** |
| Structural Analysis | 0.3 ms | Negligible |
| Solver | 1,550 ms | Includes 120 PropsSI calls |

**Root Cause Analysis**:

The first `CoolProp::PropsSI` call for any fluid triggers expensive one-time initialization:
- Pre-computing "super ancillary" equations for saturation properties
- Building equation of state infrastructure
- Loading and validating fluid parameters

For R22: **First call: 32,608 ms**, subsequent calls: **~13 ms each**

**Key Insight**: Once a fluid is initialized, CoolProp is fast (13ms/call). The variable inference phase makes the first PropsSI call, paying the initialization cost. The solver phase then benefits from the cached fluid.

**Recommendations**:
1. **Lazy fluid initialization**: Defer PropsSI calls in variable inference until actually needed
2. **AbstractState caching**: Use `CoolProp::AbstractState` objects instead of PropsSI
3. **Warm-up on demand**: Only initialize fluids that will be used in the solve

## 2. Performance Optimization: Thermodynamic Property Evaluations

The current overhead is primarily due to redundant and expensive calls to `CoolProp::PropsSI` and the use of finite differences for Jacobian computation.

### 2.1. AbstractState Caching
- **Issue**: `PropsSI` performs fluid lookup, string parsing, and EOS initialization on every call.
- **Fix**: Implement a cache of `shared_ptr<CoolProp::AbstractState>` in the `ExpressionEvaluator`.
- **Benefit**: Reusing `AbstractState` objects and calling `update()` is significantly faster for repeated evaluations at the same or nearby states.

### 2.2. Analytical Derivatives
- **Issue**: Current Jacobians use finite differences (3+ `PropsSI` calls per property).
- **Fix**: Use `AbstractState::first_partial_deriv` to obtain exact analytical derivatives.
- **Benefit**: Reduces the number of EOS evaluations to 1 per property call and improves solver robustness with exact gradients.

### 2.3. Solver Loop Optimization
- **Issue**: The Newton loop re-computes the Jacobian during line searches.
- **Fix**: Parameterize `BlockEvaluator::evaluate` to skip Jacobian computation when only the residual norm is required (e.g., during line search).


## 3. Alternative Solvers Integration

As described in the initial plan, `coolsolve` will support multiple solver backends to handle different types of equation systems:

- **Newton-Raphson (Internal)**: The current default for small to medium algebraic loops.
- **KINSOL (SUNDIALS)**: For large-scale non-linear systems requiring robust preconditioning.
- **Ceres Solver**: For systems that can be framed as optimization or least-squares problems.
- **Levenberg-Marquardt**: To improve convergence when the initial guess is far from the solution.

## 4. Implementation Status and Roadmap

| Task | Priority | Status |
|------|----------|--------|
| **Lazy Fluid Initialization** | **Critical** | **Next** |
| AbstractState Caching | High | Planned |
| Analytical Derivatives | High | Planned |
| Solving Time Tracker | Medium | Implemented |
| Code Complexity Monitor | Medium | Planned |
| KINSOL Integration | Low | Planned |
| Detailed Profiling | Low | **Complete** |

### Implemented: `--no-superancillary` Option

Disabling CoolProp's superancillary functions significantly speeds up VLE calculations during solving:

```bash
# 22x faster solving with this option
./coolsolve --no-superancillary model.eescode
```

| Mode | Solver Time |
|------|-------------|
| Default (with superancillary) | 1612 ms |
| `--no-superancillary` | 73 ms |

Trade-off: Superancillaries provide higher accuracy near the critical point and faster phase determination. For most engineering calculations, disabling them is acceptable.

### Next: Lazy Fluid Initialization

The most impactful remaining optimization is to defer the first CoolProp call:

1. **Skip PropsSI in `initializeVariables()`**: Use heuristic-based default values instead of CoolProp calls for initial guesses. Only call CoolProp during the actual solve.

2. **Pre-warm fluids lazily**: On first solver evaluation, warm up CoolProp for all fluids used in the model, then proceed with solving.

3. **Background initialization**: Start fluid initialization in a background thread while parsing/analyzing, so it completes before solving starts.
