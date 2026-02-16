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
}

const CONFIG_SCHEMA: ConfigGroup[] = [
  {
    title: 'Main Iteration',
    fields: [
      { key: 'maxIterations', label: 'Max iterations', type: 'number', defaultVal: '100', description: 'Maximum Newton iterations per block' },
      { key: 'tolerance', label: 'Tolerance', type: 'number', defaultVal: '1e-9', description: 'Convergence: stop when max residual below this' },
      { key: 'relativeTolerance', label: 'Relative tolerance', type: 'number', defaultVal: '1e-9', description: 'Stop when ||F||/||F0|| below this' },
      { key: 'stepTolerance', label: 'Step tolerance', type: 'number', defaultVal: '1e-12', description: 'Minimum step size treated as converged' },
      { key: 'verbose', label: 'Verbose', type: 'boolean', defaultVal: 'false', description: 'Print iteration info' },
    ],
  },
  {
    title: 'Line Search',
    fields: [
      { key: 'lsAlpha', label: 'Armijo α', type: 'number', defaultVal: '1e-4', description: 'Armijo condition parameter' },
      { key: 'lsRho', label: 'Step reduction ρ', type: 'number', defaultVal: '0.5', description: 'Backtracking step reduction factor' },
      { key: 'lsMaxIterations', label: 'Max LS iterations', type: 'number', defaultVal: '100', description: 'Max line search iterations per Newton step' },
      { key: 'lsMinStep', label: 'Min step', type: 'number', defaultVal: '1e-10', description: 'Minimum step size in line search' },
      { key: 'lsRelaxedTolerance', label: 'Relaxed tolerance', type: 'number', defaultVal: '1e-2', description: 'Accept as converged if LS fails and ||F|| below this' },
    ],
  },
  {
    title: 'Variable Scaling',
    fields: [
      { key: 'enableScaling', label: 'Enable scaling', type: 'boolean', defaultVal: 'true', description: 'Automatic variable scaling for better conditioning' },
    ],
  },
  {
    title: 'Trust Region',
    fields: [
      { key: 'trInitialRadius', label: 'Initial radius', type: 'number', defaultVal: '10.0', description: 'Initial trust region radius' },
      { key: 'trMaxRadius', label: 'Max radius', type: 'number', defaultVal: '1000.0', description: 'Maximum trust region radius' },
      { key: 'trEta', label: 'Acceptance η', type: 'number', defaultVal: '0.05', description: 'Step acceptance threshold' },
      { key: 'trShrinkFactor', label: 'Shrink factor', type: 'number', defaultVal: '0.5', description: 'Shrink factor when step is rejected' },
      { key: 'trGrowFactor', label: 'Grow factor', type: 'number', defaultVal: '2.0', description: 'Grow factor when step is accepted' },
    ],
  },
  {
    title: 'Levenberg-Marquardt',
    fields: [
      { key: 'lmInitialLambda', label: 'Initial λ', type: 'number', defaultVal: '1e-3', description: 'Initial damping parameter' },
      { key: 'lmLambdaIncrease', label: 'λ increase', type: 'number', defaultVal: '10.0', description: 'Factor to increase λ on bad step' },
      { key: 'lmLambdaDecrease', label: 'λ decrease', type: 'number', defaultVal: '0.1', description: 'Factor to decrease λ on good step' },
      { key: 'lmMinLambda', label: 'Min λ', type: 'number', defaultVal: '1e-12', description: 'Minimum damping parameter' },
      { key: 'lmMaxLambda', label: 'Max λ', type: 'number', defaultVal: '1e8', description: 'Maximum damping parameter' },
    ],
  },
  {
    title: 'Partitioned Solver',
    fields: [
      { key: 'partitionedMaxIterations', label: 'Max iterations', type: 'number', defaultVal: '300', description: 'Max iterations for partitioned solver' },
      { key: 'partitionedRelaxation', label: 'Relaxation', type: 'number', defaultVal: '0.6', description: 'Relaxation factor 0 < w ≤ 1' },
      { key: 'partitionedMinDiagonal', label: 'Min diagonal', type: 'number', defaultVal: '1e-12', description: 'Minimum |dF_i/dx_i| to update variable' },
      { key: 'partitionedMinBlockSize', label: 'Min block size', type: 'number', defaultVal: '4', description: 'Only use for blocks of this size or larger' },
    ],
  },
  {
    title: 'Tearing',
    fields: [
      { key: 'enableTearing', label: 'Enable tearing', type: 'boolean', defaultVal: 'true', description: 'Decompose algebraic loops via tearing' },
      { key: 'tearingMaxIterations', label: 'Max iterations', type: 'number', defaultVal: '100', description: 'Max tearing iterations' },
      { key: 'tearingMinBlockSize', label: 'Min block size', type: 'number', defaultVal: '3', description: 'Min block size for tearing' },
      { key: 'tearingInnerIterations', label: 'Inner iterations', type: 'number', defaultVal: '5', description: 'Inner iterations per tearing step' },
    ],
  },
  {
    title: 'Solver Pipeline',
    fields: [
      { key: 'solverPipeline', label: 'Pipeline', type: 'string', defaultVal: 'Newton, TrustRegion, LevenbergMarquardt, Partitioned', description: 'Comma-separated list of solvers to try' },
      { key: 'pipelineMode', label: 'Mode', type: 'string', defaultVal: 'sequential', description: 'sequential or parallel' },
    ],
  },
  {
    title: 'Safety',
    fields: [
      { key: 'timeoutSeconds', label: 'Timeout (s)', type: 'number', defaultVal: '0', description: 'Timeout in seconds; 0 = no timeout' },
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
          )}
        </div>
      ))}
    </div>
  );
}
