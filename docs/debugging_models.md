# Debugging Models in CoolSolve

This guide explains how to diagnose and fix solver failures in CoolSolve. It is based on the workflow used to fix the `condenser_3zones.eescode` model and applies to other models that fail at the solver level.

## Overview

When a model fails to solve, the typical causes are:

1. **Convergence problems** – Newton solver cannot find a solution (line search failure, max iterations, poor initial guess)
2. **Algebraic loops** – Large strongly connected components (SCCs) that are difficult to solve
3. **Invalid physics** – CoolProp errors from out-of-range temperatures, pressures, or enthalpies
4. **Procedure/function issues** – Bugs or unit mismatches in user-defined procedures

## Step 1: Run in Debug Mode

Always start by running the model in debug mode to gather diagnostic information:

```bash
./coolsolve -d input.eescode
```

This creates a folder `<input>_coolsolve/` (or `<input>_coolsolve/` next to the input file) containing:

| File | Use for diagnosis |
|------|-------------------|
| `report.md` | Block summary, solver status, largest block size |
| `solver_errors.md` | Error type (LineSearchFailed, MaxIterations, EvaluationError) |
| `solver_trace.md` | Per-block Newton iterations, residual norms |
| `evaluator.md` | Block evaluation status, max residuals at current state |
| `coolsolve_residuals.md` | Equation-by-equation LHS, RHS, residual values |
| `equations.md` | Which equations are in each block |

## Step 2: Interpret the Error

### LineSearchFailed

- **Meaning**: The Newton step direction is valid, but no step size satisfies the line search (e.g. step would push variables into invalid region, or Armijo condition not met).
- **Typical causes**: Poor initial guess, variables near physical limits, ill-conditioned Jacobian.
- **Check**: `solver_trace.md` – residual at failure, step norm. If residual is small (e.g. < 1), the solution may be very close.

### MaxIterations

- **Meaning**: Newton did not converge within the iteration limit (default 100).
- **Typical causes**: Large algebraic loop, poor initials, slow convergence.
- **Check**: `solver_trace.md` – does the residual decrease or stall? A stall near a small residual suggests a near-solution.

### EvaluationError (CoolProp)

- **Meaning**: CoolProp rejected an input (e.g. enthalpy outside valid range, pressure too low).
- **Typical causes**: Incorrect initial values, wrong units (e.g. P=1 interpreted as 1 Pa instead of 1 bar), cascaded errors from a previous block.
- **Check**: Error message for the exact invalid input. Inspect `coolsolve_residuals.md` for suspicious values (e.g. enthalpy = 1.0 instead of ~4e5 for air).

### SingularJacobian

- **Meaning**: Jacobian is rank-deficient; Newton cannot compute a unique step.
- **Typical causes**: Redundant equations, degenerate system, bad structural analysis.

## Step 3: Identify the Failing Block

From `report.md` and `solver_errors.md`:

- Note the **block index** and **size** (e.g. Block 38, size 62).
- Large blocks (tens of equations) are algebraic loops and are hardest to solve.

From `equations.md`:

- List the variables and equations in that block.
- Look for coupling: shared variables, procedure calls, thermodynamic functions.

## Step 4: Check Initial Values

- Variables without values in `.initials` default to 1.0, which is often wrong (e.g. enthalpies, mass fractions).
- Compare `evaluator.md` “Current State” with physically reasonable values.
- If you have an EES solution or reference, use those values in `.initials`.

## Step 5: Break the Loop (Simplified Model)

For large algebraic loops that fail to converge, a practical approach is to **fix one or more variables** to reduce the loop size. Use a simplified model to obtain good initials, then solve the full model.

### 5.1 Create a Simplified Model

1. Copy the original model to `model_simple.eescode`.
2. Pick a variable that strongly couples the loop (e.g. a pressure, temperature, or mass fraction) and **fix it** to a known value from `.initials` or EES.
3. Replace the equation that originally determined it with an explicit assignment.

Example (condenser_3zones):

```ees
" Original (in implicit block): P_r_su_cd is unknown "
" Simplified: fix condenser inlet pressure from initials "
P_r_su_cd=949215
```

4. Ensure the system stays **square** (same number of equations and unknowns). Fixing one variable replaces one implicit equation with one explicit one, so the count stays the same.

### 5.2 Solve the Simplified Model

```bash
cp model.initials model_simple.initials
./coolsolve -g model_simple.eescode
```

- If it solves, `model_simple.initials` is updated with the solution.

### 5.3 Transfer Initials to the Full Model

```bash
cp model_simple.initials model.initials
./coolsolve model.eescode
```

- The full model often converges with these improved initials.

### 5.4 Choosing Variables to Fix

- **Pressures** (e.g. condenser/evaporator inlet): often good candidates in HVAC/refrigeration.
- **Mass/energy fractions** (e.g. zone splits): can decouple zone balances.
- **Temperatures** at key points: useful when they are central to the loop.

Use `equations.md` and `report.md` to see which variables appear in the largest blocks.

## Step 6: Fix Model Issues

### Unit Mismatches (CoolProp)

- CoolProp expects **Pa** for pressure. If the model uses bar or kPa, convert (e.g. `P=101325` for 1 atm).
- Example: `specheat(cf$, T=t_bar_cf, P=1)` may mean 1 bar in EES but 1 Pa in CoolProp → use `P=101325`.

### Procedure Fixes

- Procedures are a common source of issues. Check:
  - Correct mapping of inputs/outputs in `CALL`.
  - Pressure and temperature units in thermodynamic calls.
  - Use of `:=` vs `=` for assignment inside procedures.

## Step 7: Solver Options (Advanced)

For difficult models, you can adjust `SolverOptions` in code (e.g. in `runner.cpp` or via a future CLI):

- `maxIterations` – increase for large blocks (e.g. 500).
- `tolerance` – loosen if near-convergence is acceptable.
- `lsRelaxedTolerance` – accept convergence when line search fails but residual is small.

## Example: condenser_3zones

The `condenser_3zones.eescode` model originally failed with `LineSearchFailed` in Block 38 (62 equations). The workflow was:

1. Run `./coolsolve -d condenser_3zones.eescode` → saw Block 38, size 62, LineSearchFailed/MaxIterations.
2. Create `condenser_3zones_simple.eescode` with `P_r_su_cd=949215` fixed.
3. Solve simple model with `-g` → it converged.
4. Copy `condenser_3zones_simple.initials` to `condenser_3zones.initials`.
5. Run full model → it converged with the new initials.

Additionally, `P=1` was corrected to `P=101325` in the `two_phase_CD` procedure for the `specheat` call (CoolProp pressure in Pa).

## Summary Checklist

- [ ] Run with `-d` to generate debug output
- [ ] Identify error type and failing block from `solver_errors.md`, `report.md`
- [ ] Inspect `solver_trace.md`, `evaluator.md`, `coolsolve_residuals.md`
- [ ] Check `.initials` and default values for unphysical guesses
- [ ] Create simplified model by fixing one or more coupling variables
- [ ] Solve simplified model with `-g` and copy initials to full model
- [ ] Fix unit mismatches (especially pressure) in procedures
- [ ] Re-run full model with improved initials
