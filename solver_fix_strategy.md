# Strategy to Fix orc_co2 Newton Solver Convergence

## Problem Analysis

### Current Situation
- **Model**: orc_co2.eescode with 159 equations, 132 blocks
- **Failing Block**: Block 77 (size 28 variables) - algebraic loop with thermodynamic equations
- **Expected W_dot_t**: 1.1744e+05, **Actual**: 9.6429e+04
- **Error**: 17.9% (unacceptable for engineering accuracy)

### Root Cause Diagnosis

From the solver trace analysis, the fundamental issue is:

1. **Highly Ill-Conditioned Jacobian**: Block 77 has condition number > 1e12
2. **Stalled Convergence**: With regularization, residual reduces from 4.9e5 to ~30 over 500 iterations
3. **Linear Convergence Rate**: Only ~0.1% reduction per iteration near the end
4. **The 99.9% Acceptance Trap**: The smart acceptance criterion accepts at 99.9% reduction, but the last 0.1% contains critical accuracy

### Why Current Solutions Fail

**Regularization (current approach)**:
- Prevents huge Newton steps (good)
- But converts quadratic convergence to linear convergence (bad)
- With damping, the solver takes 1000s of iterations to reach tight tolerance
- After 500 iterations, residual is still ~30, far from 1e-9 tolerance

**Problem**: The solver is stuck in a "stiff" region of the thermodynamic surface where:
- The linear model (Jacobian) doesn't match the nonlinear reality
- CoolProp property evaluations have discontinuities/near-singularities
- Polynomial maps for compressor/turbine create multiple local minima

---

## Proposed Solutions (Out of the Box Thinking)

### Solution 1: Variable Scaling / Nondimensionalization - IMPLEMENTED

**Concept**: The variables in Block 77 have vastly different scales:
- Temperatures: ~300 K
- Pressures: ~1e7 Pa  
- Enthalpies: ~3e5 J/kg
- Efficiencies: ~0.6 (dimensionless)

This creates a poorly scaled Jacobian. By scaling variables to O(1), we improve conditioning.

**Implementation**:
```cpp
// Scale variables before solving
T_scaled = T / 300.0;        // ~1.0
P_scaled = P / 1e7;          // ~1.0
h_scaled = h / 3e5;          // ~1.0
```

**Pros**: Simple, doesn't change algorithm, improves conditioning
**Cons**: Requires identifying characteristic scales for each variable type

---

### Solution 2: Block Variable Partitioning (DAE-style)

**Concept**: Block 77 has 28 variables but maybe not all need to be solved simultaneously.

**Analysis of Block 77**:
The block contains:
- Isentropic relations: h_is, P_is, eta_s (compressor)
- Heat exchanger: h[3], T[3], P[3], h[8], h[9]
- Solar heating: T[5], T[4], T[6]
- Turbine: eta_t, h[6], h[7], h_is_t, T[6], T[7], s[6], Dh_turb

**Idea**: Partition into 2-3 smaller blocks using tear variables or hierarchical solving.

**Implementation**: 
1. Solve compressor section first (eta_s, Dh)
2. Use results as inputs for heat exchanger + turbine
3. Iterate if needed

**Pros**: Reduces problem size, faster convergence
**Cons**: Requires structural analysis changes

---

### Solution 3: Homotopy / Continuation Method

**Concept**: Instead of solving the hard problem directly, start from an easier problem and gradually transform to the target.

**Homotopy formulation**:
```
H(x, λ) = λ * F(x) + (1-λ) * G(x) = 0
```
Where:
- F(x) = target equations (orc_co2)
- G(x) = easier equations (e.g., linearized version, or simpler model)
- λ goes from 0 to 1

**Specific for orc_co2**:
1. Start with linearized thermodynamic properties
2. Gradually introduce CoolProp nonlinearities
3. At each λ step, use solution as initial guess for next step

**Pros**: Very robust for highly nonlinear problems
**Cons**: Requires solving multiple problems, more computational effort

---

### Solution 4: Trust-Region Dogleg Method (replace line search) - IMPLEMENTED

**Concept**: Line search only works along the Newton direction. Trust regions adaptively choose between Newton and steepest descent based on model accuracy.

**Algorithm**:
1. Compute Newton step dx_n
2. Compute steepest descent step dx_sd = -gradient
3. If ||dx_n|| > trust_radius, use dogleg path
4. Adapt trust radius based on actual vs predicted reduction

**Why better**: When Newton direction points to infeasible region (CoolProp error), trust region uses smaller step or switches to gradient descent.

**Pros**: More robust than line search, better handling of nonlinearities
**Cons**: More complex implementation

---

### Solution 5: Nonlinear Preconditioning (Residual Weighting)

**Concept**: Different equations have different "importance" and scales. Weight residuals to balance their influence.

**For Block 77**:
- Thermodynamic property equations (h = Enthalpy(...)): weight ~1e-5
- Energy balance: weight ~1
- Isentropic relations: weight ~1

**Implementation**:
```cpp
// Weighted residual: W * F(x) = 0
W = diag(1/max(|F_i|, 1))
```

**Pros**: Balances equation contributions
**Cons**: Requires careful tuning

---

### Solution 6: Better Initial Guess via Simplified Model

**Concept**: The .initials file has T[3]=158.19 but the current solution has different values. Run a simplified version first.

**Simplified orc_co2**:
1. Remove polynomial maps (use constant efficiency)
2. Remove pressure drops
3. Solve this simplified model
4. Use solution as initial guess for full model

**Implementation**:
- Add preprocessing step in runner
- Create "orchestration" capability to solve sub-models

**Pros**: Uses problem structure, very effective
**Cons**: Requires model-specific knowledge

---

### Solution 7: Anderson Acceleration (Fixed-Point Iteration)

**Concept**: Instead of Newton, use fixed-point iteration with acceleration.

**Why for orc_co2**:
- Thermodynamic properties are often evaluated as: h = f(T,P)
- This is naturally a fixed-point form
- Anderson acceleration converges faster than Picard

**Algorithm**:
```
x_{k+1} = g(x_k)  // fixed point
with acceleration using m previous iterates
```

**Pros**: Doesn't need Jacobian, robust for property evaluations
**Cons**: Slower than Newton when Newton works

---

### Solution 8: Explicit Solve for Size-1 Blocks (Long-Term Plan)

**Concept**: For blocks with exactly one equation and one variable—structurally explicit assignments of the form `x = expr(known_vars)`—bypass Newton entirely and solve by direct evaluation.

**How it works**:
1. **Detection**: When solving a block of size 1, check that the equation has the form `lhs_var = rhs_expr` where the LHS is a single variable and the RHS depends only on variables already solved (external to this block).
2. **Evaluation**: Evaluate the RHS expression with current external variable values (no Jacobian, no iteration).
3. **Assignment**: Assign the result to the block variable.
4. **Validation**: Optionally verify that the residual (LHS − RHS) is acceptably small as a consistency check.

**Implementation sketch**: `tryExplicitSolve(blockIndex)` currently returns `false`. Extend it to detect the explicit pattern (variable on LHS, expression on RHS with no cycle) and, if applicable, evaluate the RHS and assign to the output variable.

**Benefits for robustness**:
- **No Jacobian exposure**: Explicit blocks never compute or invert a Jacobian, so singular/NaN Jacobians (e.g. from AD edge cases like `pow(x,y)` at x=0) cannot occur for these blocks.
- **No Newton failure modes**: No line search failures, no convergence stalls, no invalid Newton steps.
- **Deterministic**: A single evaluation always succeeds when the RHS is well-defined; no iteration and no numerical solver instability.
- **Reduced solver surface**: Fewer blocks use the nonlinear solver, reducing the chance of hitting solver pathologies.

**Benefits for computational efficiency**:
- **One evaluation vs many**: Newton typically needs several F and J evaluations per block; explicit solve needs one RHS evaluation.
- **No Jacobian computation**: Avoids building and storing the Jacobian and computing AD gradients for that block.
- **No linear algebra**: No QR or linear solve for these blocks.
- **Scalability**: In typical EES models, many blocks are size 1 (e.g. 40 of 41 in scroll_compressor). Solving them explicitly can cut solver cost substantially (often 50% or more of block evaluations).

**Typical impact**: Models like `scroll_compressor.eescode` have dozens of explicit blocks. Treating them as explicit would eliminate most Newton calls and Jacobian computations, speeding up solves and avoiding the singular Jacobian issue that currently affects block 19.

**Status**: Planned for later implementation. Does not replace the short-term AD fix (guard `pow` at x=0) but complements it by avoiding the nonlinear solver for many blocks.

---

## Recommended Strategy

### Phase 1: Variable Scaling (Quick Win)
Implement automatic variable scaling based on initial guesses:
- Scale temperatures by 300 K
- Scale pressures by 1e7 Pa  
- Scale enthalpies by 3e5 J/kg
- Keep dimensionless quantities (efficiencies) as-is

This should improve conditioning from 1e12 to ~1e4, allowing Newton to work without heavy regularization.

### Phase 2: Trust Region Method (Robustness)
Replace line search with trust region dogleg method:
- Better handling of infeasible regions
- Automatic switching between Newton and gradient descent
- More robust for CoolProp evaluation failures

### Phase 3: Hierarchical Solving (Structure Exploitation)
If Phase 1-2 insufficient, implement block partitioning:
- Identify tear variables in Block 77
- Solve in 2-3 stages
- Iterate between stages

### Phase 4: Explicit Solve for Size-1 Blocks (Long-Term)
Implement `tryExplicitSolve()` for structurally explicit blocks:
- Detect `x = expr(known_vars)` pattern in size-1 blocks
- Evaluate RHS directly and assign; skip Newton and Jacobian
- Improves both robustness and efficiency for models with many explicit equations

---

## Implementation Priority

1. **Variable Scaling** (easiest, most impact)
2. **Trust Region** (medium difficulty, robustness)
3. **Hierarchical Solving** (hardest, best accuracy)
4. **Explicit Solve for Size-1 Blocks** (medium difficulty, robustness + efficiency)

The key insight: **Don't just fix the solver algorithm - fix the problem formulation** through scaling and structure exploitation.


## Status of coolsolve_test

The best implementation leads to the following results:

| File | Parse | IR | Analysis | Solve | Value Check | Eqs | Blocks | Time (s) |
|------|-------|----|---------|-------|-------------|----:|-------:|--------:|
| condenser_3zones.eescode | OK | OK | OK | OK | OK | 99 | 50 | 0.03 |
| exchangers1.eescode | OK | OK | OK | OK | OK | 20 | 20 | 0.01 |
| exchangers2.eescode | OK | OK | OK | OK | OK | 29 | 26 | 0.00 |
| exchangers3.eescode | OK | OK | OK | OK | OK | 35 | 33 | 0.03 |
| humidair1.eescode | OK | OK | OK | OK | OK | 25 | 25 | 0.00 |
| humidair2.eescode | OK | OK | OK | OK | OK | 37 | 33 | 0.02 |
| orc_co2.eescode | OK | OK | OK | OK | OK | 139 | 112 | 0.05 |
| orc_complex.eescode | FAIL | FAIL | FAIL | FAIL | - | 0 | 0 | 0.00 |
| orc_extraction.eescode | OK | OK | OK | OK | OK | 133 | 113 | 0.08 |
| orc_r245fa.eescode | OK | OK | OK | OK | OK | 173 | 151 | 0.30 |
| orc_simple.eescode | OK | OK | OK | OK | OK | 172 | 150 | 0.10 |
| pressuredrop.eescode | OK | OK | OK | OK | OK | 26 | 26 | 0.01 |
| rankine1.eescode | OK | OK | OK | OK | OK | 30 | 30 | 0.04 |
| rankine2.eescode | OK | OK | OK | OK | OK | 45 | 42 | 0.09 |
| refrigeration1.eescode | OK | OK | OK | OK | OK | 39 | 39 | 0.05 |
| refrigeration2.eescode | OK | OK | OK | OK | OK | 39 | 39 | 0.04 |
| refrigeration3.eescode | OK | OK | OK | OK | OK | 38 | 38 | 0.01 |
| scroll_compressor.eescode | OK | OK | OK | OK | OK | 99 | 66 | 0.42 |

Any modfication of coolsolve should improve these results, or at leat not deteriorate them. This table will be updated as the code improves.