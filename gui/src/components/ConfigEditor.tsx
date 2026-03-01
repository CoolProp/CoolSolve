import { useState, useMemo, useCallback } from 'react';
import { useModelStore } from '../stores/modelStore';
import { api } from '../api/client';
import { RotateCcw } from 'lucide-react';

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
const ALL_SOLVERS = 'Newton, TrustRegion, LevenbergMarquardt, BisectionND, Homotopy, Partitioned';

const PIPELINE_PRESETS: Array<{
  id: string;
  label: string;
  pipeline: string;
  mode: string;
  description: string;
}> = [
  {
    id: 'sequential',
    label: 'Sequential',
    pipeline: ALL_SOLVERS,
    mode: 'sequential',
    description:
      'Tries all solvers in order: Newton → TrustRegion → LM → BisectionND → Homotopy → Partitioned. '
      + 'Each solver warm-starts from the best result found so far. '
      + 'BisectionND is automatically skipped for blocks larger than bisectionNDMaxBlockSize. '
      + 'Recommended for most models — robust and covers the widest range of problem types.',
  },
  {
    id: 'parallel',
    label: 'Parallel',
    pipeline: ALL_SOLVERS,
    mode: 'parallel',
    description:
      'Launches all solvers concurrently in separate threads; the first solver to converge wins. '
      + 'Useful when the best solver is unknown and you have spare CPU cores. '
      + 'Uses more memory than sequential mode.',
  },
  {
    id: 'newton-only',
    label: 'Newton',
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
  // Nothing set → sequential (default)
  if (!confMap.has('solverPipeline') && !confMap.has('pipelineMode')) return 'sequential';
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
    ],
  },
  {
    title: 'Variable Scaling',
    fields: [
      { key: 'enableScaling', label: 'Enable scaling', type: 'boolean', defaultVal: 'true',
        description: 'Automatic variable scaling for better Jacobian conditioning. Recommended for models mixing variables of very different magnitudes (e.g. Pa and m³/kg).' },
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
    title: 'Solver Pipeline',
    fields: [],
    custom: 'pipeline',
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
        <select
          className="config-input"
          value={currentPresetId}
          onChange={(e) => handlePresetChange(e.target.value)}
        >
          {PIPELINE_PRESETS.filter((p) => p.id !== 'custom').map((p) => (
            <option key={p.id} value={p.id}>{p.label}</option>
          ))}
          {isCustom && <option value="custom">Custom (current conf)</option>}
        </select>
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
                      <label title={field.description}>
                        <span className="config-field-label">{field.label}</span>
                        <span className="config-field-default">default: {field.defaultVal}</span>
                      </label>
                      {field.type === 'boolean' ? (
                        <select
                          value={isSet ? current : ''}
                          onChange={(e) => handleChange(field.key, e.target.value)}
                          className="config-input"
                        >
                          <option value="">— default —</option>
                          <option value="true">true</option>
                          <option value="false">false</option>
                        </select>
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
