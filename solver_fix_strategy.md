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

---

## Implementation Priority

1. **Variable Scaling** (easiest, most impact)
2. **Trust Region** (medium difficulty, robustness)
3. **Hierarchical Solving** (hardest, best accuracy)

The key insight: **Don't just fix the solver algorithm - fix the problem formulation** through scaling and structure exploitation.


## Status of coolsolve_test

The best implementation leads to the following results:
| File | Parse | IR | Analysis | Solve | Value Check | Eqs | Blocks |
|------|-------|----|---------|-------|-------------|----:|-------:|
| condenser_3zones.eescode | OK | OK | OK | OK | OK | 99 | 50 |
| exchangers1.eescode | OK | OK | OK | OK | OK | 20 | 20 |
| exchangers2.eescode | OK | OK | OK | OK | OK | 29 | 26 |
| exchangers3.eescode | OK | OK | OK | OK | OK | 35 | 33 |
| humidair1.eescode | OK | OK | OK | OK | OK | 25 | 25 |
| humidair2.eescode | OK | OK | OK | OK | OK | 37 | 33 |
| orc_co2.eescode | OK | OK | OK | FAIL | - | 159 | 132 |
| orc_complex.eescode | FAIL | FAIL | FAIL | FAIL | - | 0 | 0 |
| orc_extraction.eescode | OK | OK | OK | OK | OK | 133 | 113 |
| orc_r245fa.eescode | OK | OK | OK | OK | OK | 173 | 151 |
| orc_simple.eescode | OK | OK | OK | OK | OK | 172 | 150 |
| pressuredrop.eescode | OK | OK | OK | OK | OK | 26 | 26 |
| rankine1.eescode | OK | OK | OK | OK | OK | 30 | 30 |
| rankine2.eescode | OK | OK | OK | OK | OK | 45 | 42 |
| refrigeration1.eescode | OK | OK | OK | OK | OK | 39 | 39 |
| refrigeration2.eescode | OK | OK | OK | OK | OK | 39 | 39 |
| refrigeration3.eescode | OK | OK | OK | OK | OK | 38 | 38 |
| scroll_compressor.eescode | OK | OK | OK | FAIL | - | 99 | 66 |

Any modfication of coolsolve should improve these results, or at leat not deteriorate them. This table will be updated as the code improves.