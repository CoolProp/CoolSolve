import { useState, useMemo, useCallback } from 'react';
import { useModelStore } from '../stores/modelStore';
import { api } from '../api/client';
import { RotateCcw, Info } from 'lucide-react';
import Tooltip from './Tooltip';
import ConfigSelect from './ConfigSelect';

interface ConfigField {
  key: string;
  label: string;
  type: 'number' | 'boolean' | 'string';
  defaultVal: string;
  description: string;
}

interface ConfigGroup {
  title: string;
  fields: ConfigField[];
  /** Render this group with a custom component instead of the standard field list */
  custom?: 'pipeline';
}

// ---------------------------------------------------------------------------
// Pipeline presets — each preset sets both solverPipeline and pipelineMode.
// ---------------------------------------------------------------------------
const ALL_SOLVERS = 'Newton, TrustRegion, LevenbergMarquardt, BisectionND, Homotopy, Partitioned, Kinsol';

const PIPELINE_PRESETS: Array<{
  id: string;
  label: string;
  pipeline: string;
  mode: string;
  description: string;
}> = [
  {
    id: 'sequential',
    label: 'All Solvers (Sequential)',
    pipeline: ALL_SOLVERS,
    mode: 'sequential',
    description:
      'Tries all solvers in order: Newton → TrustRegion → LM → BisectionND → Homotopy → Partitioned → Kinsol. '
      + 'Each solver warm-starts from the best result found so far. '
      + 'BisectionND is automatically skipped for blocks larger than bisectionNDMaxBlockSize. '
      + 'Robust and covers the widest range of problem types.',
  },
  {
    id: 'parallel',
    label: 'All Solvers (Parallel)',
    pipeline: ALL_SOLVERS,
    mode: 'parallel',
    description:
      'Launches all solvers concurrently in separate threads; the first solver to converge wins. '
      + 'Useful when the best solver is unknown and you have spare CPU cores. '
      + 'Uses more memory than sequential mode.',
  },
  {
    id: 'newton-only',
    label: 'Newton (Default)',
    pipeline: 'Newton',
    mode: 'sequential',
    description:
      'Newton solver with backtracking line search. '
      + 'Fastest per block — ideal when the Jacobian is well-conditioned and '
      + 'initial guesses are reasonably close to the solution. No fallback.',
  },
  {
    id: 'trustregion-only',
    label: 'TrustRegion',
    pipeline: 'TrustRegion',
    mode: 'sequential',
    description:
      'Trust-Region Dogleg solver. Blends steepest-descent and Newton steps, '
      + 'keeping updates inside a safe radius. Helps when trial points would drive '
      + 'thermodynamic calls into physically invalid regions. No fallback.',
  },
  {
    id: 'lm-only',
    label: 'LevenbergMarquardt',
    pipeline: 'LevenbergMarquardt',
    mode: 'sequential',
    description:
      'Levenberg-Marquardt solver. Uses adaptive damping: large λ → gradient descent '
      + '(safe, slow); small λ → Gauss-Newton (fast, quadratic near solution). '
      + 'Effective when the initial guess is far from the solution or the Jacobian is near-singular. No fallback.',
  },
  {
    id: 'bisectionnd-only',
    label: 'BisectionND',
    pipeline: 'BisectionND',
    mode: 'sequential',
    description:
      'Multi-dimensional Bisection (derivative-free). Converges even when J(x₀)=0. '
      + 'Automatically skipped for blocks larger than bisectionNDMaxBlockSize (default: 8). '
      + 'Cost is O(2ⁿ) function evaluations per iteration — use only for small blocks. No fallback.',
  },
  {
    id: 'homotopy-only',
    label: 'Homotopy',
    pipeline: 'Homotopy',
    mode: 'sequential',
    description:
      'Homotopy continuation. Deforms the problem from an easy auxiliary system to the real '
      + 'system via a parameter t ∈ [0,1]. Can reach solutions from starting points that are '
      + 'too far away for Newton-based methods. No fallback.',
  },
  {
    id: 'partitioned-only',
    label: 'Partitioned',
    pipeline: 'Partitioned',
    mode: 'sequential',
    description:
      'Partitioned Block Updates. Updates each variable using its matched equation diagonal: '
      + 'xᵢ ← xᵢ − w·Fᵢ/(∂Fᵢ/∂xᵢ). Acts like a DAE-style tear without restructuring the block. '
      + 'Last-resort stabilizer for stiff or highly nonlinear loops. No fallback.',
  },
  {
    id: 'kinsol-only',
    label: 'KINSOL',
    pipeline: 'Kinsol',
    mode: 'sequential',
    description:
      'SUNDIALS-KINSOL-style solver. Three globalisation modes (set via kinsolGlobalStrategy): '
      + 'inexact Newton + Dennis-Schnabel line search (default), Picard/Richardson fixed-point, '
      + 'or Anderson-accelerated fixed point. The derivative-free modes can solve blocks where '
      + 'the Jacobian is zero or singular. No fallback.',
  },
  // Sentinel: shown only when the conf contains values that match no preset above.
  {
    id: 'custom',
    label: 'Custom…',
    pipeline: '',
    mode: 'sequential',
    description:
      'The current coolsolve.conf contains a custom pipeline or mode that does not match '
      + 'any preset. Edit the values below, or select a preset to replace them.',
  },
];

/** Infer which preset matches the current conf entries, or "custom" if none match. */
function detectPipelinePreset(confMap: Map<string, string>): string {
  // Nothing set → Newton (default)
  if (!confMap.has('solverPipeline') && !confMap.has('pipelineMode')) return 'newton-only';
  const pipeline = (confMap.get('solverPipeline') ?? '').trim().toLowerCase();
  const mode = (confMap.get('pipelineMode') ?? 'sequential').trim().toLowerCase();

  // 1) Check single-solver presets (exact pipeline match + mode)
  for (const p of PIPELINE_PRESETS) {
    if (p.id === 'custom' || p.id === 'sequential' || p.id === 'parallel') continue;
    if (p.pipeline.toLowerCase() === pipeline && p.mode.toLowerCase() === mode) return p.id;
  }
  // 2) Sequential / parallel: match on mode only (solver list is user-editable)
  if (mode === 'parallel') return 'parallel';
  // If mode is sequential and pipeline has more than one solver, keep "sequential"
  if (mode === 'sequential') {
    const solverCount = pipeline.split(',').filter((s) => s.trim()).length;
    if (solverCount !== 1) return 'sequential';
  }
  return 'custom';
}

// ---------------------------------------------------------------------------
// Config schema — all groups except pipeline use the standard field renderer.
// ---------------------------------------------------------------------------
const CONFIG_SCHEMA: ConfigGroup[] = [
  {
    title: 'Main Iteration',
    fields: [
      { key: 'maxIterations', label: 'Max iterations', type: 'number', defaultVal: '100',
        description: 'Maximum Newton (or bisection) iterations per block. Also used as the base for bisectionNDIterFactor.' },
      { key: 'tolerance', label: 'Tolerance', type: 'number', defaultVal: '1e-9',
        description: 'Convergence: stop when max residual ‖F‖∞ is below this value.' },
      { key: 'relativeTolerance', label: 'Relative tolerance', type: 'number', defaultVal: '1e-9',
        description: 'Stop when ‖F‖/‖F₀‖ falls below this (relative to the initial residual).' },
      { key: 'stepTolerance', label: 'Step tolerance', type: 'number', defaultVal: '1e-12',
        description: 'Minimum step size treated as converged (prevents infinite looping on stuck iterates).' },
      { key: 'verbose', label: 'Verbose', type: 'boolean', defaultVal: 'false',
        description: 'Print per-iteration residual info to stderr.' },
    ],
  },
  {
    title: 'Line Search',
    fields: [
      { key: 'lsAlpha', label: 'Armijo α', type: 'number', defaultVal: '1e-4',
        description: 'Armijo sufficient-decrease parameter. Smaller values accept more aggressive steps.' },
      { key: 'lsRho', label: 'Step reduction ρ', type: 'number', defaultVal: '0.5',
        description: 'Backtracking step reduction factor (0 < ρ < 1). Applied each time a step is rejected.' },
      { key: 'lsMaxIterations', label: 'Max LS iterations', type: 'number', defaultVal: '100',
        description: 'Max line search iterations per Newton step.' },
      { key: 'lsMinStep', label: 'Min step', type: 'number', defaultVal: '1e-10',
        description: 'Minimum step size in line search; smaller steps are rejected.' },
      { key: 'lsRelaxedTolerance', label: 'Relaxed tolerance', type: 'number', defaultVal: '1e-2',
        description: 'Accept as converged if line search fails but ‖F‖ is below this. Prevents discarding near-converged solutions.' },
      { key: 'lsNonMonotoneMemory', label: 'Non-monotone memory M', type: 'number', defaultVal: '10',
        description:
          'Non-monotone line search memory (Grippo et al. 1986). '
          + 'Instead of requiring the merit function to decrease at every step (monotone Armijo), '
          + 'the solver compares against the maximum merit over the last M iterations, '
          + 'bounded at 10× the current merit to prevent catastrophic step acceptance. '
          + 'This helps escape narrow curved valleys and saddle points that trap monotone methods. '
          + 'Applied to Newton (Armijo condition), TrustRegion and LM (step acceptance). '
          + 'M=1 = classic monotone behavior. M=10 = default non-monotone (recommended). '
          + 'Increase M if the solver oscillates near a valley; decrease toward 1 if it accepts bad steps too aggressively.' },
    ],
  },
  {
    title: 'Variable Scaling',
    fields: [
      { key: 'enableScaling', label: 'Enable scaling', type: 'boolean', defaultVal: 'true',
        description: 'Automatic variable scaling for better Jacobian conditioning. Recommended for models mixing variables of very different magnitudes (e.g. Pa and m³/kg).' },
      { key: 'broydenRecomputeInterval', label: 'Broyden interval', type: 'number', defaultVal: '0',
        description:
          'Broyden quasi-Newton enhancement for the Newton solver. '
          + 'When > 0, the Newton solver computes a full Jacobian every K iterations '
          + 'and uses Broyden rank-1 updates in between, saving expensive Jacobian evaluations. '
          + 'If the Broyden approximation fails line search, a full Jacobian is recomputed automatically. '
          + '0 = disabled (classic full-Jacobian Newton, default). '
          + '5 = recompute every 5 iterations (good starting point for large CoolProp models).' },
    ],
  },
  {
    title: 'Trust Region',
    fields: [
      { key: 'trInitialRadius', label: 'Initial radius', type: 'number', defaultVal: '10.0',
        description: 'Initial trust region radius. Larger → more aggressive initial steps.' },
      { key: 'trMaxRadius', label: 'Max radius', type: 'number', defaultVal: '1000.0',
        description: 'Maximum trust region radius. Caps step size to avoid driving variables into unphysical ranges.' },
      { key: 'trEta', label: 'Acceptance η', type: 'number', defaultVal: '0.05',
        description: 'Step acceptance threshold (ρ ≥ η). Lower → accept more steps; higher → require better agreement between model and quadratic.' },
      { key: 'trShrinkFactor', label: 'Shrink factor', type: 'number', defaultVal: '0.5',
        description: 'Multiply trust radius by this factor when a step is rejected.' },
      { key: 'trGrowFactor', label: 'Grow factor', type: 'number', defaultVal: '2.0',
        description: 'Multiply trust radius by this factor when a step is accepted with ρ ≥ 0.9.' },
      { key: 'trAdaptiveRadius', label: 'Adaptive radius', type: 'boolean', defaultVal: 'true',
        description:
          'When true, the initial trust radius is set from the Cauchy step norm on the first iteration '
          + 'instead of using the fixed trInitialRadius. This scales the trust region to the problem '
          + 'geometry automatically. Also enables smoother rho-based radius adaptation '
          + 'and gradient-based recovery after consecutive rejections. Recommended for most models.' },
      { key: 'trBroydenRecomputeInterval', label: 'Broyden interval', type: 'number', defaultVal: '0',
        description:
          'hybrd-style (Powell 1970) Broyden quasi-Newton enhancement for TrustRegion. '
          + 'When > 0, computes a full Jacobian every K iterations and maintains a Broyden rank-1 '
          + 'approximation (via an incrementally-updated QR factorization) in between, saving expensive '
          + 'Jacobian evaluations. A full Jacobian is also recomputed after consecutive rejected steps '
          + '(see Broyden restart rejects) or if the resulting step is non-finite. '
          + '0 = disabled (classic full-Jacobian TrustRegion, default). '
          + '5 = recompute every 5 iterations (good starting point for large CoolProp models).' },
      { key: 'trBroydenRestartRejects', label: 'Broyden restart rejects', type: 'number', defaultVal: '2',
        description:
          'Number of consecutive rejected steps that triggers a full Jacobian recompute (Powell restart '
          + 'criterion) when Broyden interval > 0. Lower values recover from a stale approximation sooner '
          + 'at the cost of more Jacobian evaluations.' },
    ],
  },
  {
    title: 'Levenberg-Marquardt',
    fields: [
      { key: 'lmInitialLambda', label: 'Initial λ', type: 'number', defaultVal: '1e-3',
        description: 'Initial damping parameter. Large λ → gradient descent (safe, slow). Small λ → Gauss-Newton (fast, quadratic near solution).' },
      { key: 'lmLambdaIncrease', label: 'λ increase', type: 'number', defaultVal: '10.0',
        description: 'Multiply λ by this when a step increases ‖F‖ (bad step).' },
      { key: 'lmLambdaDecrease', label: 'λ decrease', type: 'number', defaultVal: '0.1',
        description: 'Multiply λ by this when a step decreases ‖F‖ (good step).' },
      { key: 'lmMinLambda', label: 'Min λ', type: 'number', defaultVal: '1e-12',
        description: 'Floor for the damping parameter.' },
      { key: 'lmMaxLambda', label: 'Max λ', type: 'number', defaultVal: '1e8',
        description: 'Ceiling for the damping parameter. If λ exceeds this, the solver gives up.' },
      { key: 'lmNielsenUpdate', label: 'Nielsen update', type: 'boolean', defaultVal: 'true',
        description:
          'Nielsen\'s smooth lambda adaptation (Madsen et al. 2004). '
          + 'Uses λ = λ × max(1/3, 1 − (2ρ−1)³) on acceptance and λ = λ × ν with ν doubling '
          + 'on consecutive rejections. Provides smoother, faster-converging lambda transitions '
          + 'than the legacy step-wise increase/decrease. Recommended.' },
      { key: 'lmGeodesicAcceleration', label: 'Geodesic accel.', type: 'boolean', defaultVal: 'true',
        description:
          'Geodesic acceleration (Transtrum & Sethna 2012). Adds a second-order correction '
          + 'to the LM step by evaluating the directional second derivative of F along the velocity step. '
          + 'Costs 1 extra residual evaluation per iteration but can halve iterations on curved problems. '
          + 'Automatically disabled when the acceleration term dominates the velocity.' },
    ],
  },
  {
    title: 'BisectionND',
    fields: [
      { key: 'bisectionNDMaxBlockSize', label: 'Max block size', type: 'number', defaultVal: '8',
        description:
          'Maximum number of unknowns for which BisectionND is attempted. '
          + 'Blocks larger than this are automatically skipped (the next pipeline solver is tried). '
          + 'BisectionND requires 2ⁿ probe evaluations per phase — cost grows exponentially. '
          + 'Default: 8. Increase with caution (n=10 → 1024 probes, n=12 → 4096 probes).' },
      { key: 'bisectionNDIterFactor', label: 'Iteration factor', type: 'number', defaultVal: '1.0',
        description:
          'Multiplier for the BisectionND iteration budget: actual steps = round(maxIterations × bisectionNDIterFactor). '
          + 'BisectionND converges linearly (much slower than Newton), so it typically needs more steps. '
          + 'Example: maxIterations=100, bisectionNDIterFactor=5 → 500 bisection steps. '
          + 'Increase if BisectionND reports "max iterations reached" on a problem you know has a solution.' },
    ],
  },
  {
    title: 'Partitioned Solver',
    fields: [
      { key: 'partitionedMaxIterations', label: 'Max iterations', type: 'number', defaultVal: '300',
        description: 'Max iterations for the partitioned (per-variable diagonal update) solver.' },
      { key: 'partitionedRelaxation', label: 'Relaxation', type: 'number', defaultVal: '0.6',
        description: 'Relaxation factor w (0 < w ≤ 1) applied to each variable update. Smaller values improve stability on ill-conditioned loops.' },
      { key: 'partitionedMinDiagonal', label: 'Min diagonal', type: 'number', defaultVal: '1e-12',
        description: 'Minimum |∂Fᵢ/∂xᵢ| to update a variable. Prevents division by near-zero diagonals.' },
      { key: 'partitionedMinBlockSize', label: 'Min block size', type: 'number', defaultVal: '4',
        description: 'Only use the partitioned solver for blocks with at least this many unknowns.' },
    ],
  },
  {
    title: 'KINSOL (SUNDIALS-style)',
    fields: [
      { key: 'kinsolGlobalStrategy', label: 'Global strategy', type: 'string', defaultVal: 'linesearch',
        description:
          'Globalisation strategy for the KINSOL solver (add "Kinsol" to the pipeline). '
          + 'One of: linesearch (inexact Newton + Dennis-Schnabel line search; default), '
          + 'picard (fixed-point/Richardson iteration, no Jacobian), '
          + 'fp (Anderson-accelerated fixed point, derivative-free).' },
      { key: 'kinsolLineSearchAlpha', label: 'LS sufficient decrease α', type: 'number', defaultVal: '1e-4',
        description:
          'Armijo sufficient-decrease coefficient α for the KINSOL line search (linesearch mode). '
          + 'Range (0,1). Smaller → stricter decrease.' },
      { key: 'kinsolLineSearchMaxIters', label: 'LS max backtracks', type: 'number', defaultVal: '30',
        description: 'Maximum number of backtracking trials per Newton step in the KINSOL line search.' },
      { key: 'kinsolPicardOmega', label: 'Picard ω', type: 'number', defaultVal: '1.0',
        description:
          'Relaxation ω for Picard and Anderson modes: x ← x − ω·F(x). Must be > 0. '
          + 'Use ω < 1 to under-relax stiff (near-divergent) fixed-point iterations.' },
      { key: 'kinsolAndersonDepth', label: 'Anderson depth m', type: 'number', defaultVal: '5',
        description:
          'History depth m for Anderson acceleration (fp mode). '
          + 'm = 0 disables acceleration (degenerates to plain fixed point).' },
      { key: 'kinsolAndersonRelaxation', label: 'Anderson damping θ', type: 'number', defaultVal: '1.0',
        description:
          'Damping θ for Anderson acceleration: x_next = x + θ·(x_Anderson − x). '
          + 'θ = 1 → pure Anderson (default); θ < 1 → more conservative. Range (0,1].' },
    ],
  },
  {
    title: 'Tearing',
    fields: [
      { key: 'enableTearing', label: 'Enable tearing', type: 'boolean', defaultVal: 'true',
        description:
          'When true, blocks of size ≥ tearingMinBlockSize are pre-processed by structural tearing '
          + 'BEFORE the pipeline solvers. Tearing selects a minimal "tear set" (feedback vertex set) '
          + 'so that the remaining equations form an acyclic system solvable sequentially. '
          + 'Newton is then applied only to the small tear-variable residual. '
          + 'If tearing succeeds, pipeline solvers are skipped for that block. '
          + 'If tearing fails, the pipeline runs normally. '
          + 'Enable for large, densely-coupled loops (ORC, heat exchangers). '
          + 'Tearing applies regardless of which solvers are in solverPipeline.' },
      { key: 'tearingMaxIterations', label: 'Max iterations', type: 'number', defaultVal: '100',
        description: 'Max Newton iterations on the tear-variable residual system.' },
      { key: 'tearingMinBlockSize', label: 'Min block size', type: 'number', defaultVal: '3',
        description: 'Only apply tearing to blocks with at least this many unknowns.' },
      { key: 'tearingInnerIterations', label: 'Inner iterations', type: 'number', defaultVal: '5',
        description: 'Max solver iterations per equation in the acyclic (sequential) solve pass.' },
    ],
  },
  {
    title: 'Symbolic Reduction',
    fields: [
      { key: 'enableSymbolicReduction', label: 'Enable reduction', type: 'boolean', defaultVal: 'false',
        description:
          'Pre-process blocks of size ≥ 2 to reduce their size before the iterative solver. '
          + 'Techniques: (1) extract explicit equations where the output is computable from known values, '
          + '(2) invert CoolProp calls so that an unknown input becomes the output '
          + '(e.g. T = temperature(Water, H=h, P=P) instead of h = enthalpy(Water, T=T, P=P)), '
          + '(3) substitute variables that appear only in their own defining equation. '
          + 'Can dramatically reduce block sizes for CoolProp-heavy models (refrigeration, ORC, heat exchangers). '
          + 'Off by default; when disabled, zero overhead is added.' },
    ],
  },
  {
    title: 'Multi-Start',
    fields: [
      { key: 'multiStartEnabled', label: 'Enable multi-start', type: 'boolean', defaultVal: 'true',
        description:
          'When a multi-variable block fails the entire solver pipeline, retry it from a few '
          + 'alternative starting points. Thermo blocks re-evaluate every property via CoolProp at '
          + 'a coherent reference (T, P) state (pressure regimes); purely-algebraic blocks scale the '
          + 'default guess (wrong magnitude is the typical failure, e.g. C ~ 0.05 guessed as 1.0). '
          + 'Zero overhead when every block converges on the first try: the search only triggers after '
          + 'a block failure. Size-1 blocks are skipped (Newton1D already does its own multi-probe).' },
      { key: 'multiStartMaxRestarts', label: 'Max restarts', type: 'number', defaultVal: '4',
        description:
          'Number of alternative starting points to try on a failed block. Each candidate replays '
          + 'the full solver pipeline, so large values increase the worst-case cost of a failure. '
          + 'Range 1..6.' },
      { key: 'multiStartNumCores', label: 'Num cores', type: 'number', defaultVal: '1',
        description:
          'Number of threads used to evaluate multi-start candidates concurrently. '
          + '1 = sequential (default, no threading overhead); N>1 = run up to N candidates '
          + 'concurrently (first-to-converge wins); 0 = auto (min(hardware_concurrency, candidates)). '
          + 'Only consulted when multi-start engages (a block fails), so successful solves are unaffected.' },
    ],
  },
  {
    title: 'Solver Pipeline',
    fields: [],
    custom: 'pipeline',
  },
  {
    title: 'Integration',
    fields: [
      { key: 'integralMethod', label: 'Method', type: 'string', defaultVal: 'RK4',
        description:
          'Time-marching integrator for equation-based INTEGRAL / $IntegralTable models. '
          + 'One of: EulerExplicit (order 1, cheap, conditionally stable), '
          + 'EulerImplicit (order 1, A-stable, good for stiff systems), '
          + 'RK4 (classic 4th-order Runge-Kutta; fixed step, recommended default), '
          + 'RK45 (Dormand-Prince adaptive embedded pair; variable step, error-controlled). '
          + 'Ignored for non-integral models.' },
      { key: 'integralFixedStep', label: 'Fixed step', type: 'number', defaultVal: '0',
        description:
          'Fixed time step. 0 ⇒ derive the step from integralMaxSteps and the integration interval. '
          + 'Used by EulerExplicit, EulerImplicit, RK4 (and by RK45 as the initial step guess).' },
      { key: 'integralMaxSteps', label: 'Max steps', type: 'number', defaultVal: '1000',
        description: 'Upper bound on the number of integration steps. The march aborts if exceeded.' },
      { key: 'integralRelTol', label: 'Relative tol.', type: 'number', defaultVal: '1e-6',
        description:
          'RK45 local-error control (relative component). '
          + 'The step is accepted when the weighted error ≤ relTol·|y| + absTol. Ignored for fixed-step methods.' },
      { key: 'integralAbsTol', label: 'Absolute tol.', type: 'number', defaultVal: '1e-9',
        description: 'RK45 absolute-error floor. Prevents over-tightening the step when |y| ≈ 0.' },
      { key: 'integralMinStep', label: 'Min step', type: 'number', defaultVal: '0',
        description: 'Minimum step size for RK45 (0 = auto, derived from the interval). ' },
      { key: 'integralMaxStep', label: 'Max step', type: 'number', defaultVal: '0',
        description: 'Maximum step size for RK45 (0 = auto = the full interval).' },
      { key: 'integralRichardson', label: 'Richardson', type: 'boolean', defaultVal: 'false',
        description:
          'Richardson extrapolation (fixed-step methods only). Runs one h step + two h/2 steps and '
          + 'combines them as (2^p·I_{h/2} − I_h)/(2^p − 1) where p is the method order, '
          + 'reducing truncation error by 1–2 orders at ~3× cost.' },
      { key: 'integralOutputInterval', label: 'Output interval', type: 'number', defaultVal: '0',
        description:
          'Default $IntegralTable output interval when the directive omits the ":n". '
          + '0 = record every step. Does not change the integration — only which rows are tabulated.' },
    ],
  },
  {
    title: 'CoolProp Integration',
    fields: [
      { key: 'coolpropBackend', label: 'Backend', type: 'string', defaultVal: 'HEOS',
        description:
          'CoolProp AbstractState backend. HEOS = full Helmholtz EOS (default). '
          + 'TTSE&HEOS or BICUBIC&HEOS = tabular interpolation backends that build lookup '
          + 'tables on first use (slower startup, much faster evaluations for large models).' },
      { key: 'coolpropUseAbstractState', label: 'Use AbstractState', type: 'boolean', defaultVal: 'true',
        description:
          'Use the low-level AbstractState API instead of the high-level PropsSI function. '
          + 'Eliminates string parsing and fluid lookup overhead on every call. '
          + 'When false, the legacy PropsSI path is used with zero additional overhead.' },
      { key: 'coolpropEnableAnalyticalDerivatives', label: 'Analytical derivatives', type: 'boolean', defaultVal: 'true',
        description:
          'Compute CoolProp derivatives via first_partial_deriv() instead of 4 central '
          + 'finite-difference calls per property evaluation. Much faster and more accurate. '
          + 'Requires Use AbstractState = true. When false, finite differences are used.' },
      { key: 'coolpropCacheEnabled', label: 'Cache enabled', type: 'boolean', defaultVal: 'true',
        description:
          'Cache AbstractState instances across evaluations (thread-local). '
          + 'Avoids re-creating the state object for every CoolProp call.' },
      { key: 'coolpropEnableSuperancillaries', label: 'Superancillaries', type: 'boolean', defaultVal: 'true',
        description:
          'Enable CoolProp superancillary functions for faster VLE initialisation. '
          + 'Adds a small startup cost but accelerates two-phase lookups.' },
    ],
  },
  {
    title: 'LaTeX Report',
    fields: [
      { key: 'enableLatexReport', label: 'Enable LaTeX report', type: 'boolean', defaultVal: 'true',
        description:
          'Generate a LaTeX (.tex) report after every successful solve. '
          + 'The report contains model equations, variable solutions, and plots. '
          + 'The backend only produces the LaTeX source — no LaTeX installation is required on the solver side.' },
      { key: 'latexCompiler', label: 'LaTeX compiler', type: 'string', defaultVal: 'pdflatex',
        description:
          'LaTeX compiler command used for PDF compilation. '
          + 'Common alternatives: xelatex, lualatex.' },
    ],
  },
  {
    title: 'Safety',
    fields: [
      { key: 'timeoutSeconds', label: 'Timeout (s)', type: 'number', defaultVal: '0',
        description: 'Global timeout in seconds; 0 = no timeout. Aborts the entire solve if exceeded.' },
    ],
  },
];

/** Parse conf text into key→value map */
function parseConf(text: string): Map<string, string> {
  const map = new Map<string, string>();
  for (const line of text.split('\n')) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith('#')) continue;
    const eq = trimmed.indexOf('=');
    if (eq < 0) continue;
    map.set(trimmed.substring(0, eq).trim(), trimmed.substring(eq + 1).trim());
  }
  return map;
}

/** Serialize key→value map to conf text */
function serializeConf(map: Map<string, string>): string {
  const lines: string[] = ['# CoolSolve solver configuration'];
  for (const [key, val] of map) {
    lines.push(`${key} = ${val}`);
  }
  return lines.join('\n') + '\n';
}

// ---------------------------------------------------------------------------
// Pipeline group — custom dropdown component
// ---------------------------------------------------------------------------
interface PipelineGroupProps {
  confMap: Map<string, string>;
  onChange: (key: string, value: string) => void;
  onBatchChange: (changes: Record<string, string>) => void;
}

function PipelineGroup({ confMap, onChange, onBatchChange }: PipelineGroupProps) {
  const currentPresetId = detectPipelinePreset(confMap);
  const isCustom = currentPresetId === 'custom';
  const showPipelineList = currentPresetId === 'sequential' || currentPresetId === 'parallel';
  const currentPreset = PIPELINE_PRESETS.find((p) => p.id === currentPresetId) ?? PIPELINE_PRESETS[0];

  const handlePresetChange = (newId: string) => {
    if (newId === 'custom') return;
    const preset = PIPELINE_PRESETS.find((p) => p.id === newId);
    if (!preset) return;
    onBatchChange({ solverPipeline: preset.pipeline, pipelineMode: preset.mode });
  };

  return (
    <div className="config-group-body">
      {/* Preset dropdown */}
      <div className="config-field">
        <label title="Choose the solver pipeline preset">
          <span className="config-field-label">Pipeline</span>
          <span className="config-field-default">default: Sequential</span>
        </label>
        <ConfigSelect
          value={currentPresetId}
          onChange={handlePresetChange}
          options={[
            ...PIPELINE_PRESETS.filter((p) => p.id !== 'custom').map((p) => ({
              value: p.id,
              label: p.label,
            })),
            ...(isCustom ? [{ value: 'custom', label: 'Custom (current conf)' }] : []),
          ]}
        />
      </div>

      {/* Description */}
      <div style={{ fontSize: '0.78rem', color: 'var(--color-text-muted, #888)', padding: '2px 2px 6px 2px', lineHeight: '1.45' }}>
        {currentPreset.description}
      </div>

      {/* Solver list — only for multi-solver presets (sequential / parallel) or custom */}
      {(showPipelineList || isCustom) && (
        <div style={{ padding: '4px 0' }}>
          <div style={{ fontSize: '0.82rem', marginBottom: '4px' }}>
            <span className="config-field-label" title="Comma-separated list of solvers: Newton, TrustRegion, LM, BisectionND, Homotopy, Partitioned">
              Solver list
            </span>
          </div>
          <textarea
            className="config-input"
            rows={2}
            style={{ resize: 'vertical', width: '100%', boxSizing: 'border-box', fontFamily: 'inherit', fontSize: 'inherit', display: 'block' }}
            value={confMap.get('solverPipeline') ?? ''}
            placeholder="Newton, TrustRegion, LevenbergMarquardt, BisectionND, Homotopy, Partitioned"
            onChange={(e) => onChange('solverPipeline', e.target.value)}
          />
        </div>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------
export default function ConfigEditor() {
  const conf = useModelStore((s) => s.conf);
  const setConf = useModelStore((s) => s.setConf);
  const [collapsed, setCollapsed] = useState<Set<string>>(new Set());

  const confMap = useMemo(() => parseConf(conf), [conf]);

  const toggleGroup = useCallback((title: string) => {
    setCollapsed((prev) => {
      const next = new Set(prev);
      if (next.has(title)) next.delete(title);
      else next.add(title);
      return next;
    });
  }, []);

  const handleChange = useCallback(
    (key: string, value: string) => {
      const newMap = new Map(confMap);
      if (value === '' || value === undefined) {
        newMap.delete(key);
      } else {
        newMap.set(key, value);
      }
      const text = serializeConf(newMap);
      setConf(text);
      api.putConf(text).catch(() => {});
    },
    [confMap, setConf]
  );

  const handleBatchChange = useCallback(
    (changes: Record<string, string>) => {
      const newMap = new Map(confMap);
      for (const [key, value] of Object.entries(changes)) {
        if (value === '' || value === undefined) {
          newMap.delete(key);
        } else {
          newMap.set(key, value);
        }
      }
      const text = serializeConf(newMap);
      setConf(text);
      api.putConf(text).catch(() => {});
    },
    [confMap, setConf]
  );

  const handleReset = useCallback(() => {
    setConf('');
    api.putConf('').catch(() => {});
  }, [setConf]);

  return (
    <div className="config-editor">
      <div className="config-header">
        <span className="config-title">coolsolve.conf</span>
        <button className="toolbar-btn" onClick={handleReset} title="Reset all to defaults">
          <RotateCcw size={14} /> Reset
        </button>
      </div>
      {CONFIG_SCHEMA.map((group) => (
        <div key={group.title} className="config-group">
          <div className="config-group-header" onClick={() => toggleGroup(group.title)}>
            <span>{collapsed.has(group.title) ? '▸' : '▾'} {group.title}</span>
          </div>
          {!collapsed.has(group.title) && (
            group.custom === 'pipeline' ? (
              <PipelineGroup confMap={confMap} onChange={handleChange} onBatchChange={handleBatchChange} />
            ) : (
              <div className="config-group-body">
                {group.fields.map((field) => {
                  const current = confMap.get(field.key) ?? '';
                  const isSet = confMap.has(field.key);
                  return (
                    <div key={field.key} className="config-field">
                      <label>
                        <span className="config-field-label">
                          {field.label}
                          {field.description && (
                            <Tooltip content={field.description}>
                              <Info size={12} className="config-info-icon" />
                            </Tooltip>
                          )}
                        </span>
                        <span className="config-field-default">default: {field.defaultVal}</span>
                      </label>
                      {field.type === 'boolean' ? (
                        <ConfigSelect
                          value={isSet ? current : ''}
                          onChange={(v) => handleChange(field.key, v)}
                          options={[
                            { value: '', label: '— default —' },
                            { value: 'true', label: 'true' },
                            { value: 'false', label: 'false' },
                          ]}
                        />
                      ) : (
                        <input
                          type="text"
                          className="config-input"
                          value={isSet ? current : ''}
                          placeholder={field.defaultVal}
                          onChange={(e) => handleChange(field.key, e.target.value)}
                        />
                      )}
                    </div>
                  );
                })}
              </div>
            )
          )}
        </div>
      ))}
    </div>
  );
}
