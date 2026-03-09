/**
 * ParametricStudy.tsx — Parametric sensitivity analysis panel.
 *
 * Allows users to:
 *  1. Select 1 or 2 imposed variables to sweep
 *  2. Define sweep ranges (min, max, steps)
 *  3. Run the study (backend solves each grid point with warm-starting)
 *  4. View results in a table (swept variables highlighted)
 *  5. Plot 1D line charts or 2D contour/surface plots
 *  6. Save/load studies in the model bundle
 */
import { useState, useMemo, useCallback, useRef, useEffect } from 'react';
import Plot from './PlotlyChart';
import { useModelStore } from '../stores/modelStore';
import { api } from '../api/client';
import type {
  VariableInfo,
  SweepVariable,
  ParametricStudyResult,
  SavedParametricStudy,
} from '../api/types';

// ============================================================================
// Helpers
// ============================================================================

/** Generate a linear range of values */
function linspace(min: number, max: number, n: number): number[] {
  if (n <= 1) return [min];
  const step = (max - min) / (n - 1);
  return Array.from({ length: n }, (_, i) => min + i * step);
}

/** Generate a unique ID for saved studies */
function makeId(): string {
  return Date.now().toString(36) + Math.random().toString(36).substring(2, 6);
}

/** Format a number for display */
function fmt(val: number): string {
  if (Number.isInteger(val) && Math.abs(val) < 1e9) return val.toString();
  if (Math.abs(val) < 0.01 || Math.abs(val) >= 1e6) return val.toExponential(4);
  return val.toPrecision(5);
}

// ============================================================================
// Sub-components
// ============================================================================

interface SweepInputProps {
  label: string;
  varName: string;
  setVarName: (n: string) => void;
  min: string;
  setMin: (v: string) => void;
  max: string;
  setMax: (v: string) => void;
  steps: string;
  setSteps: (v: string) => void;
  imposedVars: VariableInfo[];
  disabled: boolean;
}

function SweepInput({
  label, varName, setVarName, min, setMin, max, setMax, steps, setSteps,
  imposedVars, disabled,
}: SweepInputProps) {
  // When user selects a variable, pre-fill min/max from its imposed value
  const handleVarChange = (name: string) => {
    setVarName(name);
    const info = imposedVars.find((v) => v.name === name);
    if (info?.imposedValue !== undefined) {
      const val = info.imposedValue;
      // Default range: ±30% around imposed value (or ±1 if val is 0)
      const delta = val !== 0 ? Math.abs(val) * 0.3 : 1;
      setMin((val - delta).toPrecision(4));
      setMax((val + delta).toPrecision(4));
    }
  };

  return (
    <div className="sweep-input-row">
      <label className="sweep-label">{label}</label>
      <select
        className="sweep-select"
        value={varName}
        onChange={(e) => handleVarChange(e.target.value)}
        disabled={disabled}
      >
        <option value="">— select variable —</option>
        {imposedVars.map((v) => (
          <option key={v.name} value={v.name}>
            {v.name} {v.imposedValue !== undefined ? `(= ${fmt(v.imposedValue)})` : ''}
          </option>
        ))}
      </select>
      <input
        className="sweep-num" type="text" placeholder="min"
        value={min} onChange={(e) => setMin(e.target.value)} disabled={disabled}
      />
      <input
        className="sweep-num" type="text" placeholder="max"
        value={max} onChange={(e) => setMax(e.target.value)} disabled={disabled}
      />
      <input
        className="sweep-num sweep-steps" type="text" placeholder="steps"
        value={steps} onChange={(e) => setSteps(e.target.value)} disabled={disabled}
      />
    </div>
  );
}

// ============================================================================
// Main Component
// ============================================================================

export default function ParametricStudy() {
  const parsedVariables = useModelStore((s) => s.parsedVariables);
  const parametricStudies = useModelStore((s) => s.parametricStudies);
  const addParametricStudy = useModelStore((s) => s.addParametricStudy);
  const removeParametricStudy = useModelStore((s) => s.removeParametricStudy);

  // Sweep configuration
  const [var1Name, setVar1Name] = useState('');
  const [var1Min, setVar1Min] = useState('');
  const [var1Max, setVar1Max] = useState('');
  const [var1Steps, setVar1Steps] = useState('10');
  const [var2Name, setVar2Name] = useState('');
  const [var2Min, setVar2Min] = useState('');
  const [var2Max, setVar2Max] = useState('');
  const [var2Steps, setVar2Steps] = useState('10');

  // Options
  const [timeout, setTimeout] = useState('1');
  const [updateGuesses, setUpdateGuesses] = useState(false);

  // Run state
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState('');
  const [error, setError] = useState('');

  // Active result (from run or loaded study)
  const [activeResult, setActiveResult] = useState<ParametricStudyResult | null>(null);
  const [activeStudyId, setActiveStudyId] = useState<string | null>(null);

  // Plot config
  const [plotYVar, setPlotYVar] = useState('');
  const [plotType, setPlotType] = useState<'line' | 'contour'>('line');

  const sseRef = useRef<EventSource | null>(null);

  // Imposed variables available for sweep
  const imposedVars = useMemo(
    () => parsedVariables.filter((v) => v.isImposed),
    [parsedVariables]
  );

  // All output variable names from the active result
  const outputVarNames = useMemo(() => {
    if (!activeResult || activeResult.results.length === 0) return [];
    const firstSuccess = activeResult.results.find((r) => r.success && r.variables);
    if (!firstSuccess?.variables) return [];
    return Object.keys(firstSuccess.variables).sort((a, b) =>
      a.localeCompare(b, undefined, { sensitivity: 'base' })
    );
  }, [activeResult]);

  // Determine number of sweep dimensions
  const sweepDim = useMemo(() => {
    if (!activeResult) return 0;
    return activeResult.sweepVariables.length;
  }, [activeResult]);

  // Reset plotYVar when result changes
  useEffect(() => {
    if (outputVarNames.length > 0 && !outputVarNames.includes(plotYVar)) {
      setPlotYVar(outputVarNames[0]);
    }
  }, [outputVarNames, plotYVar]);

  // Cleanup SSE on unmount
  useEffect(() => {
    return () => {
      if (sseRef.current) sseRef.current.close();
    };
  }, []);

  // ---- Run parametric study ----
  const handleRun = useCallback(async () => {
    const sweeps: SweepVariable[] = [];

    const n1 = parseInt(var1Steps, 10) || 10;
    const min1 = parseFloat(var1Min);
    const max1 = parseFloat(var1Max);
    if (!var1Name || isNaN(min1) || isNaN(max1)) {
      setError('Please select a variable and enter valid min/max for Sweep 1.');
      return;
    }
    sweeps.push({ name: var1Name, values: linspace(min1, max1, n1) });

    if (var2Name) {
      const n2 = parseInt(var2Steps, 10) || 10;
      const min2 = parseFloat(var2Min);
      const max2 = parseFloat(var2Max);
      if (isNaN(min2) || isNaN(max2)) {
        setError('Please enter valid min/max for Sweep 2.');
        return;
      }
      sweeps.push({ name: var2Name, values: linspace(min2, max2, n2) });
    }

    setError('');
    setRunning(true);
    setProgress('Starting parametric study...');
    setActiveResult(null);
    setActiveStudyId(null);

    try {
      // Start the parametric study
      const timeoutSec = Math.max(0, Math.round(parseFloat(timeout) || 0));
      await api.runParametricStudy(sweeps, {
        timeout: timeoutSec,
        updateGuesses,
      });

      // Subscribe to SSE for progress
      const es = new EventSource('/api/v1/solve/stream');
      sseRef.current = es;

      es.onmessage = (msg) => {
        try {
          const data = JSON.parse(msg.data);
          if (data.type === 'progress' || data.type === 'block') {
            setProgress(data.message || `Point ${data.block || '?'}...`);
          } else if (data.type === 'done') {
            es.close();
            sseRef.current = null;
            // Fetch result
            api.getParametricResult()
              .then((result) => {
                setActiveResult(result);
                setRunning(false);
                setProgress('');
                // Auto-save
                const study: SavedParametricStudy = {
                  id: makeId(),
                  name: sweeps.map((s) => s.name).join(' × '),
                  sweepVariables: sweeps,
                  result,
                  timestamp: Date.now(),
                };
                addParametricStudy(study);
                setActiveStudyId(study.id);
              })
              .catch((err) => {
                setError('Failed to fetch results: ' + err.message);
                setRunning(false);
              });
          } else if (data.type === 'error') {
            es.close();
            sseRef.current = null;
            // Partial success: some points may have solved — fetch result anyway
            setError(data.message || 'Parametric study failed');
            api.getParametricResult()
              .then((result) => {
                setActiveResult(result);
                setRunning(false);
                if (result.successCount > 0) {
                  const study: SavedParametricStudy = {
                    id: makeId(),
                    name: sweeps.map((s) => s.name).join(' × '),
                    sweepVariables: sweeps,
                    result,
                    timestamp: Date.now(),
                  };
                  addParametricStudy(study);
                  setActiveStudyId(study.id);
                }
              })
              .catch(() => {
                setRunning(false);
              });
          }
        } catch { /* ignore parse errors */ }
      };

      es.onerror = () => {
        es.close();
        sseRef.current = null;
        // If still running, try to fetch result (SSE may close prematurely)
        if (running) {
          api.getParametricResult()
            .then((result) => {
              setActiveResult(result);
              setRunning(false);
            })
            .catch(() => {
              setError('Lost connection during parametric study.');
              setRunning(false);
            });
        }
      };
    } catch (err: unknown) {
      setError(err instanceof Error ? err.message : 'Failed to start parametric study');
      setRunning(false);
    }
  }, [var1Name, var1Min, var1Max, var1Steps, var2Name, var2Min, var2Max, var2Steps, running, addParametricStudy, timeout, updateGuesses]);

  // ---- Stop parametric study ----
  const handleStop = useCallback(async () => {
    try {
      await api.cancelSolve();
    } catch { /* ignore */ }
  }, []);

  // ---- Load a saved study ----
  const loadStudy = useCallback((study: SavedParametricStudy) => {
    setActiveResult(study.result);
    setActiveStudyId(study.id);
    // Restore sweep config
    if (study.sweepVariables[0]) {
      const s = study.sweepVariables[0];
      setVar1Name(s.name);
      setVar1Min(s.values[0]?.toString() || '');
      setVar1Max(s.values[s.values.length - 1]?.toString() || '');
      setVar1Steps(s.values.length.toString());
    }
    if (study.sweepVariables[1]) {
      const s = study.sweepVariables[1];
      setVar2Name(s.name);
      setVar2Min(s.values[0]?.toString() || '');
      setVar2Max(s.values[s.values.length - 1]?.toString() || '');
      setVar2Steps(s.values.length.toString());
    } else {
      setVar2Name('');
    }
  }, []);

  // ---- Build Plot ----
  const plotData = useMemo((): { data: Plotly.Data[]; layout: Partial<Plotly.Layout> } => {
    if (!activeResult || !plotYVar) return { data: [], layout: {} };

    const sweepNames = activeResult.sweepVariables.map((s) => s.name);
    const dim = sweepNames.length;

    if (dim === 1) {
      // 1D line chart
      const sweep = activeResult.sweepVariables[0];
      const xVals: number[] = [];
      const yVals: number[] = [];
      const markers: string[] = [];

      for (const pt of activeResult.results) {
        if (!pt.success || !pt.variables) continue;
        const xVal = pt.overrides[sweep.name];
        const yVal = pt.variables[plotYVar];
        if (xVal !== undefined && yVal !== undefined) {
          xVals.push(xVal);
          yVals.push(yVal);
          markers.push(`${sweep.name}=${fmt(xVal)}, ${plotYVar}=${fmt(yVal)}`);
        }
      }

      return {
        data: [{
          x: xVals,
          y: yVals,
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          marker: { size: 6, color: '#3b82f6' },
          line: { color: '#3b82f6', width: 2 },
          text: markers,
          hoverinfo: 'text' as const,
          name: plotYVar,
        }],
        layout: {
          xaxis: { title: { text: sweep.name }, gridcolor: '#e5e7eb' },
          yaxis: { title: { text: plotYVar }, gridcolor: '#e5e7eb' },
          margin: { t: 30, r: 30, b: 50, l: 60 },
          paper_bgcolor: 'transparent',
          plot_bgcolor: 'transparent',
          font: { family: 'system-ui, sans-serif', size: 11, color: 'var(--text-primary)' },
          hovermode: 'closest' as const,
        },
      };
    }

    if (dim === 2) {
      // 2D contour or heatmap
      const s0 = activeResult.sweepVariables[0];
      const s1 = activeResult.sweepVariables[1];
      const xVals = s0.values;
      const yVals = s1.values;

      // Build z-matrix (yVals × xVals)
      const z: (number | null)[][] = [];
      for (let j = 0; j < yVals.length; j++) {
        const row: (number | null)[] = [];
        for (let i = 0; i < xVals.length; i++) {
          const pt = activeResult.results.find(
            (r) => r.overrides[s0.name] === xVals[i] && r.overrides[s1.name] === yVals[j]
          );
          row.push(pt?.success && pt.variables ? (pt.variables[plotYVar] ?? null) : null);
        }
        z.push(row);
      }

      if (plotType === 'contour') {
        return {
          data: [{
            x: xVals,
            y: yVals,
            z,
            type: 'contour' as const,
            colorscale: 'Viridis',
            colorbar: { title: plotYVar },
          } as Plotly.Data],
          layout: {
            xaxis: { title: { text: s0.name } },
            yaxis: { title: { text: s1.name } },
            margin: { t: 30, r: 30, b: 50, l: 60 },
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: { family: 'system-ui, sans-serif', size: 11 },
          },
        };
      }

      // Heatmap (default for 2D)
      return {
        data: [{
          x: xVals,
          y: yVals,
          z,
          type: 'heatmap' as const,
          colorscale: 'Viridis',
          colorbar: { title: plotYVar },
        } as Plotly.Data],
        layout: {
          xaxis: { title: { text: s0.name } },
          yaxis: { title: { text: s1.name } },
          margin: { t: 30, r: 30, b: 50, l: 60 },
          paper_bgcolor: 'transparent',
          plot_bgcolor: 'transparent',
          font: { family: 'system-ui, sans-serif', size: 11 },
        },
      };
    }

    return { data: [], layout: {} };
  }, [activeResult, plotYVar, plotType]);

  // ---- Render ----
  return (
    <div className="parametric-study">
      {/* ---- SETUP SECTION ---- */}
      <div className="parametric-setup">
        <div className="parametric-section-title">Parametric Sweep</div>

        {imposedVars.length === 0 ? (
          <div className="parametric-notice">
            No imposed variables detected. Write equations like <code>T_in = 300</code> in your
            model, then parse to enable parametric studies.
          </div>
        ) : (
          <>
            <SweepInput
              label="Sweep 1"
              varName={var1Name} setVarName={setVar1Name}
              min={var1Min} setMin={setVar1Min}
              max={var1Max} setMax={setVar1Max}
              steps={var1Steps} setSteps={setVar1Steps}
              imposedVars={imposedVars}
              disabled={running}
            />
            <SweepInput
              label="Sweep 2"
              varName={var2Name} setVarName={setVar2Name}
              min={var2Min} setMin={setVar2Min}
              max={var2Max} setMax={setVar2Max}
              steps={var2Steps} setSteps={setVar2Steps}
              imposedVars={imposedVars}
              disabled={running}
            />

            <div className="parametric-options-row">
              <label className="parametric-option" title="Maximum time per solve point (seconds). 0 = no limit.">
                Timeout (s):
                <input
                  className="sweep-num sweep-timeout"
                  type="text"
                  value={timeout}
                  onChange={(e) => setTimeout(e.target.value)}
                  disabled={running}
                />
              </label>
              <label className="parametric-option" title="When checked, each solved point updates the initial guesses for the next point. Points are reordered to spiral outward from the value closest to the original initials.">
                <input
                  type="checkbox"
                  checked={updateGuesses}
                  onChange={(e) => setUpdateGuesses(e.target.checked)}
                  disabled={running}
                />
                Update guesses
              </label>
            </div>

            <div className="parametric-actions">
              <button
                className="parametric-run-btn"
                disabled={running || !var1Name}
                onClick={handleRun}
              >
                {running ? (
                  <>
                    <span className="spinner" style={{ width: 14, height: 14 }} />
                    Running...
                  </>
                ) : (
                  'Run Study'
                )}
              </button>
              {running && (
                <button
                  className="parametric-stop-btn"
                  onClick={handleStop}
                  title="Stop the parametric study"
                >
                  Stop
                </button>
              )}
              {running && (
                <span className="parametric-progress">{progress}</span>
              )}
            </div>
          </>
        )}

        {error && <div className="parametric-error">{error}</div>}
      </div>

      {/* ---- SAVED STUDIES ---- */}
      {parametricStudies.length > 0 && (
        <div className="parametric-saved">
          <div className="parametric-section-title">Saved Studies</div>
          <div className="parametric-study-list">
            {parametricStudies.map((s) => (
              <div
                key={s.id}
                className={`parametric-study-item ${activeStudyId === s.id ? 'active' : ''}`}
                onClick={() => loadStudy(s)}
              >
                <span className="study-name">{s.name}</span>
                <span className="study-meta">
                  {s.result.successCount}/{s.result.totalPoints} pts
                </span>
                <button
                  className="study-delete"
                  onClick={(e) => { e.stopPropagation(); removeParametricStudy(s.id); }}
                  title="Delete study"
                >×</button>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* ---- RESULTS ---- */}
      {activeResult && (
        <div className="parametric-results">
          <div className="parametric-section-title">
            Results — {activeResult.successCount}/{activeResult.totalPoints} converged
          </div>

          {/* Plot controls */}
          <div className="parametric-plot-controls">
            <label>Y variable:</label>
            <select
              className="sweep-select"
              value={plotYVar}
              onChange={(e) => setPlotYVar(e.target.value)}
            >
              {outputVarNames.map((n) => (
                <option key={n} value={n}>{n}</option>
              ))}
            </select>
            {sweepDim === 2 && (
              <>
                <label>Plot type:</label>
                <select
                  className="sweep-select"
                  value={plotType}
                  onChange={(e) => setPlotType(e.target.value as 'line' | 'contour')}
                >
                  <option value="contour">Contour</option>
                  <option value="line">Heatmap</option>
                </select>
              </>
            )}
          </div>

          {/* Plot */}
          {plotData.data.length > 0 && (
            <div className="parametric-plot">
              <Plot
                data={plotData.data}
                layout={plotData.layout}
                config={{ responsive: true, displaylogo: false }}
                useResizeHandler
                style={{ width: '100%', height: 320 }}
              />
            </div>
          )}

          {/* Results table */}
          <div className="parametric-table-container">
            <table className="variable-table parametric-table">
              <thead>
                <tr>
                  <th>#</th>
                  {activeResult.sweepVariables.map((s) => (
                    <th key={s.name} className="sweep-col">{s.name}</th>
                  ))}
                  <th>Status</th>
                  {outputVarNames.slice(0, 20).map((n) => (
                    <th key={n}>{n}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {activeResult.results.map((pt, idx) => (
                  <tr key={idx} className={pt.success ? '' : 'failed-row'}>
                    <td>{idx + 1}</td>
                    {activeResult.sweepVariables.map((s) => (
                      <td key={s.name} className="sweep-col">{fmt(pt.overrides[s.name])}</td>
                    ))}
                    <td className={pt.success ? 'status-ok' : 'status-fail'}>
                      {pt.success ? '✓' : '✗'}
                    </td>
                    {outputVarNames.slice(0, 20).map((n) => (
                      <td key={n} className="value-cell">
                        {pt.success && pt.variables ? fmt(pt.variables[n]) : '—'}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
}
