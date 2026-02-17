import { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import Plot from './PlotlyChart';
import { api } from '../api/client';
import { useModelStore } from '../stores/modelStore';
import type {
  FluidInfo,
  SaturationResponse,
  InferredVariable,
} from '../api/types';

type DiagramType = 'T-s' | 'P-h' | 'h-s' | 'T-h';

interface DiagramConfig {
  xProp: 'T' | 'P' | 'H' | 'S' | 'D';
  yProp: 'T' | 'P' | 'H' | 'S' | 'D';
  xLabel: string;
  yLabel: string;
  yLog: boolean;
}

const DIAGRAM_CONFIGS: Record<DiagramType, DiagramConfig> = {
  'T-s': { xProp: 'S', yProp: 'T', xLabel: 'Entropy [J/(kg·K)]', yLabel: 'Temperature [°C]', yLog: false },
  'P-h': { xProp: 'H', yProp: 'P', xLabel: 'Enthalpy [J/kg]', yLabel: 'Pressure [Pa]', yLog: true },
  'h-s': { xProp: 'S', yProp: 'H', xLabel: 'Entropy [J/(kg·K)]', yLabel: 'Enthalpy [J/kg]', yLog: false },
  'T-h': { xProp: 'H', yProp: 'T', xLabel: 'Enthalpy [J/kg]', yLabel: 'Temperature [°C]', yLog: false },
};

const PROP_UNITS: Record<string, string> = {
  T: '[°C]',
  P: '[Pa]',
  H: '[J/kg]',
  S: '[J/(kg·K)]',
  D: '[kg/m³]',
  Q: '',
};

export default function ThermoDiagram() {
  const [fluids, setFluids] = useState<FluidInfo[]>([]);
  const [modelFluids, setModelFluids] = useState<string[]>([]);
  const [selectedFluid, setSelectedFluid] = useState('');
  const [diagramType, setDiagramType] = useState<DiagramType>('T-s');
  const [satData, setSatData] = useState<SaturationResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [computeTime, setComputeTime] = useState<number | null>(null);

  // State overlay
  const [showStates, setShowStates] = useState(true);
  const [inferredVars, setInferredVars] = useState<InferredVariable[]>([]);

  // Array overlay (from solved array variables)
  const [showArrayOverlay, setShowArrayOverlay] = useState(false);
  const [arrayXVar, setArrayXVar] = useState('');
  const [arrayYVar, setArrayYVar] = useState('');
  const [closeCycle, setCloseCycle] = useState(false);

  const lastResult = useModelStore((s) => s.lastResult);
  const plotRef = useRef<HTMLDivElement>(null);

  // Load fluids list
  useEffect(() => {
    api.getCoolPropFluids().then((res) => {
      setFluids(res.fluids);
      setModelFluids(res.modelFluids);
      // Default: first model fluid with dome, or first fluid with dome
      const withDome = res.fluids.filter((f) => f.hasDome);
      const modelWithDome = withDome.filter((f) => res.modelFluids.includes(f.name));
      if (modelWithDome.length > 0) {
        setSelectedFluid(modelWithDome[0].name);
      } else if (withDome.length > 0) {
        setSelectedFluid(withDome[0].name);
      }
    }).catch(() => {});
  }, []);

  // Reload fluids when model changes (after solve)
  useEffect(() => {
    if (lastResult) {
      api.getCoolPropFluids().then((res) => {
        setModelFluids(res.modelFluids);
      }).catch(() => {});
    }
  }, [lastResult]);

  // Load inferred variables when we have a solve result
  useEffect(() => {
    if (lastResult?.success) {
      api.getInferredVariables().then((res) => {
        setInferredVars(res.variables);
      }).catch(() => {});
    } else {
      setInferredVars([]);
    }
  }, [lastResult]);

  // Array columns from solve result (all array variables, not just inferred)
  const arrayColumns = useMemo(() => {
    if (!lastResult?.variables) return new Map<string, number[]>();
    const cols = new Map<string, { indices: number[]; values: number[] }>();
    for (const [name, val] of Object.entries(lastResult.variables)) {
      const match = name.match(/^(.+)\[(\d+)\]$/);
      if (match) {
        const base = match[1];
        const idx = parseInt(match[2]);
        if (!cols.has(base)) cols.set(base, { indices: [], values: [] });
        cols.get(base)!.indices.push(idx);
        cols.get(base)!.values.push(val);
      }
    }
    // Sort each column by index
    const result = new Map<string, number[]>();
    for (const [base, data] of cols) {
      const sorted = data.indices
        .map((idx, i) => ({ idx, val: data.values[i] }))
        .sort((a, b) => a.idx - b.idx);
      result.set(base, sorted.map((s) => s.val));
    }
    return result;
  }, [lastResult]);

  // Sorted fluids: model fluids first, then alphabetical
  const sortedFluids = useMemo(() => {
    const withDome = fluids.filter((f) => f.hasDome);
    const model = withDome.filter((f) => modelFluids.includes(f.name));
    const other = withDome.filter((f) => !modelFluids.includes(f.name));
    return [...model, ...other];
  }, [fluids, modelFluids]);

  const generate = useCallback(async () => {
    if (!selectedFluid) return;
    setLoading(true);
    setError('');
    setSatData(null);
    setComputeTime(null);
    try {
      const res = await api.getSaturationDome(selectedFluid);
      setSatData(res);
      setComputeTime(res.computeTime_ms);
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setLoading(false);
    }
  }, [selectedFluid]);

  // Build Plotly traces
  const traces = useMemo(() => {
    if (!satData) return [];
    const cfg = DIAGRAM_CONFIGS[diagramType];
    const t: Plotly.Data[] = [];

    // Saturation dome: liquid line + vapor line (connected)
    const liqX = satData.liquid[cfg.xProp];
    const liqY = satData.liquid[cfg.yProp];
    const vapX = [...satData.vapor[cfg.xProp]].reverse();
    const vapY = [...satData.vapor[cfg.yProp]].reverse();

    t.push({
      x: [...liqX, ...vapX],
      y: [...liqY, ...vapY],
      type: 'scatter',
      mode: 'lines',
      name: 'Saturation dome',
      line: { color: '#3b82f6', width: 2 },
      hovertemplate: `${cfg.xLabel}: %{x:.4g}<br>${cfg.yLabel}: %{y:.4g}<extra>Saturation</extra>`,
    });

    // Critical point
    t.push({
      x: [satData.critical[cfg.xProp]],
      y: [satData.critical[cfg.yProp]],
      type: 'scatter',
      mode: 'markers',
      name: 'Critical point',
      marker: { color: '#ef4444', size: 10, symbol: 'diamond' },
      hovertemplate: `Critical<br>${cfg.xLabel}: %{x:.4g}<br>${cfg.yLabel}: %{y:.4g}<extra></extra>`,
    });

    // State overlay: solved variables for this fluid
    if (showStates && inferredVars.length > 0) {
      const statesForFluid = inferredVars.filter(
        (v) => !v.isArray && v.inferredFluid === selectedFluid
      );

      // Group by state point (variables with same base name)
      const stateMap = new Map<string, Record<string, { value: number; varName: string }>>();
      for (const v of statesForFluid) {
        // Try to group by base name: e.g., T_1, P_1, H_1 -> state "1"
        const baseName = v.name.replace(/^[A-Za-z]+_?/, '');
        const key = baseName || v.name;
        if (!stateMap.has(key)) stateMap.set(key, {});
        stateMap.get(key)![v.inferredProperty] = { value: v.value, varName: v.name };
      }

      // Collect state points that have both x and y properties
      const stateX: number[] = [];
      const stateY: number[] = [];
      const stateLabels: string[] = [];
      const stateHovers: string[] = [];

      for (const [key, props] of stateMap) {
        if (cfg.xProp in props && cfg.yProp in props) {
          stateX.push(props[cfg.xProp].value);
          stateY.push(props[cfg.yProp].value);
          stateLabels.push(key);
          // Build tooltip showing all variable names and values for this state
          const lines: string[] = [];
          for (const [prop, info] of Object.entries(props)) {
            const unit = PROP_UNITS[prop] || '';
            lines.push(`${info.varName} = ${info.value.toPrecision(5)} ${unit}`);
          }
          stateHovers.push(lines.join('<br>'));
        }
      }

      if (stateX.length > 0) {
        t.push({
          x: stateX,
          y: stateY,
          type: 'scatter',
          mode: 'text+markers',
          name: 'State points',
          marker: { color: '#22c55e', size: 10, symbol: 'circle' },
          text: stateLabels,
          textposition: 'top center',
          textfont: { size: 11, color: '#22c55e' },
          hovertemplate: stateHovers.map(
            (h) => `${h}<extra></extra>`
          ),
        });
      }
    }

    // Array overlay (from solve result arrays)
    if (showArrayOverlay && arrayXVar && arrayYVar) {
      const xData = arrayColumns.get(arrayXVar);
      const yData = arrayColumns.get(arrayYVar);

      if (xData && yData) {
        const len = Math.min(xData.length, yData.length);
        const pathX = xData.slice(0, len);
        const pathY = yData.slice(0, len);

        if (closeCycle && pathX.length > 0) {
          pathX.push(pathX[0]);
          pathY.push(pathY[0]);
        }

        if (pathX.length > 0) {
          t.push({
            x: pathX,
            y: pathY,
            type: 'scatter',
            mode: 'lines+markers',
            name: 'Array data',
            line: { color: '#f59e0b', width: 2.5, dash: 'dot' },
            marker: { color: '#f59e0b', size: 7 },
            hovertemplate: `${arrayXVar}: %{x:.4g}<br>${arrayYVar}: %{y:.4g}<extra>Array</extra>`,
          });
        }
      }
    }

    return t;
  }, [satData, diagramType, showStates, inferredVars, selectedFluid, showArrayOverlay, arrayXVar, arrayYVar, closeCycle, arrayColumns]);

  const layout = useMemo((): Partial<Plotly.Layout> => {
    const cfg = DIAGRAM_CONFIGS[diagramType];
    return {
      title: { text: satData ? `${satData.fluid} — ${diagramType} diagram` : '' },
      xaxis: { title: { text: cfg.xLabel }, gridcolor: '#ddd', zerolinecolor: '#bbb' },
      yaxis: {
        title: { text: cfg.yLabel },
        type: cfg.yLog ? 'log' : 'linear',
        gridcolor: '#ddd',
        zerolinecolor: '#bbb',
      },
      paper_bgcolor: 'rgba(255,255,255,0)',
      plot_bgcolor: '#ffffff',
      font: { color: '#333', size: 12 },
      legend: { x: 0.02, y: 0.98, bgcolor: 'rgba(255,255,255,0.8)' },
      margin: { l: 70, r: 30, t: 50, b: 60 },
      autosize: true,
    };
  }, [diagramType, satData]);

  return (
    <div className="thermo-diagram" style={{ display: 'flex', flexDirection: 'column', height: '100%', gap: 8, padding: 8 }}>
      {/* Controls */}
      <div style={{ display: 'flex', gap: 8, alignItems: 'center', flexWrap: 'wrap' }}>
        <select
          value={selectedFluid}
          onChange={(e) => {
            setSelectedFluid(e.target.value);
            setSatData(null);
          }}
          style={{ padding: '4px 8px', minWidth: 140 }}
        >
          {modelFluids.length > 0 && (
            <optgroup label="Model fluids">
              {sortedFluids
                .filter((f) => modelFluids.includes(f.name))
                .map((f) => (
                  <option key={f.name} value={f.name}>
                    {f.name}
                  </option>
                ))}
            </optgroup>
          )}
          <optgroup label="All fluids">
            {sortedFluids
              .filter((f) => !modelFluids.includes(f.name))
              .map((f) => (
                <option key={f.name} value={f.name}>
                  {f.name}
                </option>
              ))}
          </optgroup>
        </select>

        <select
          value={diagramType}
          onChange={(e) => setDiagramType(e.target.value as DiagramType)}
          style={{ padding: '4px 8px' }}
        >
          <option value="T-s">T-s</option>
          <option value="P-h">P-h</option>
          <option value="h-s">h-s</option>
          <option value="T-h">T-h</option>
        </select>

        <button
          onClick={generate}
          disabled={loading || !selectedFluid}
          className="toolbar-btn"
          style={{ padding: '4px 12px' }}
        >
          {loading ? 'Computing...' : 'Generate'}
        </button>

        {computeTime !== null && (
          <span style={{ fontSize: 11, color: '#888' }}>
            {computeTime.toFixed(0)} ms
          </span>
        )}
      </div>

      {/* Overlays */}
      {satData && lastResult?.success && (
        <div style={{ display: 'flex', gap: 16, alignItems: 'center', flexWrap: 'wrap', fontSize: 13 }}>
          <label style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
            <input
              type="checkbox"
              checked={showStates}
              onChange={(e) => setShowStates(e.target.checked)}
            />
            State points
          </label>

          <label style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
            <input
              type="checkbox"
              checked={showArrayOverlay}
              onChange={(e) => setShowArrayOverlay(e.target.checked)}
            />
            Array overlay
          </label>

          {showArrayOverlay && arrayColumns.size > 0 && (
            <>
              <select
                value={arrayXVar}
                onChange={(e) => setArrayXVar(e.target.value)}
                style={{ padding: '2px 6px', fontSize: 12 }}
              >
                <option value="">X column...</option>
                {[...arrayColumns.keys()].sort().map((base) => (
                  <option key={base} value={base}>
                    {base} ({arrayColumns.get(base)!.length} pts)
                  </option>
                ))}
              </select>
              <select
                value={arrayYVar}
                onChange={(e) => setArrayYVar(e.target.value)}
                style={{ padding: '2px 6px', fontSize: 12 }}
              >
                <option value="">Y column...</option>
                {[...arrayColumns.keys()].sort().map((base) => (
                  <option key={base} value={base}>
                    {base} ({arrayColumns.get(base)!.length} pts)
                  </option>
                ))}
              </select>
              <label style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
                <input
                  type="checkbox"
                  checked={closeCycle}
                  onChange={(e) => setCloseCycle(e.target.checked)}
                />
                Close loop
              </label>
            </>
          )}
        </div>
      )}

      {/* Error */}
      {error && (
        <div style={{ color: '#ef4444', fontSize: 13, padding: '4px 8px' }}>
          Error: {error}
        </div>
      )}

      {/* Plot */}
      <div ref={plotRef} style={{ flex: 1, minHeight: 200 }}>
        {satData ? (
          <Plot
            data={traces}
            layout={layout}
            config={{
              responsive: true,
              displayModeBar: true,
              toImageButtonOptions: {
                format: 'png',
                filename: `${selectedFluid}_${diagramType}`,
                height: 800,
                width: 1200,
              },
              modeBarButtonsToAdd: ['toImage'],
            }}
            style={{ width: '100%', height: '100%' }}
            useResizeHandler
          />
        ) : !loading ? (
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              height: '100%',
              color: '#666',
              fontSize: 14,
            }}
          >
            Select a fluid and click Generate to create a diagram
          </div>
        ) : null}
      </div>
    </div>
  );
}
