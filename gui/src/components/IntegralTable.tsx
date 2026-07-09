/**
 * IntegralTable.tsx — Bottom-panel tab for the equation-based dynamic
 * (INTEGRAL / $IntegralTable) trajectory.
 *
 * Mirrors the Parametric Study panel's UX: a scrollable columnar table
 * (integration variable always first), a row count, an "Export CSV"
 * button, an optional line plot of one tabulated variable versus the
 * integration variable, and an empty-state message when no integral has
 * been solved. Read-only — EES Integral Tables are not user-editable.
 *
 * Data sources (in priority order):
 *   1. `integralTable` (columnar JSON from a live solve response)
 *   2. `integralCSV`   (raw CSV text, restored from a ZIP bundle when the
 *                       columnar JSON is unavailable)
 *
 * The CSV round-trip is the persistence contract: download bundle →
 * upload bundle → the same columns and rows reappear here (see
 * `tests/test_integral_api.py`).
 */
import { useMemo, useState, useCallback, useEffect } from 'react';
import Plot from './PlotlyChart';
import { useModelStore } from '../stores/modelStore';
import type { IntegralTableData } from '../api/types';

// ============================================================================
// Helpers
// ============================================================================

/** Format a number for display (matches ParametricStudy's fmt). */
function fmt(val: number): string {
  if (!Number.isFinite(val)) return '—';
  if (Number.isInteger(val) && Math.abs(val) < 1e9) return val.toString();
  if (Math.abs(val) < 0.01 || Math.abs(val) >= 1e6) return val.toExponential(4);
  return val.toPrecision(6);
}

/** Parse a CSV string (header row + numeric rows) into IntegralTableData. */
function parseCSVToTable(csv: string, csvName: string): IntegralTableData | null {
  if (!csv || !csv.trim()) return null;
  const lines = csv.split(/\r?\n/).filter((l) => l.trim().length > 0);
  if (lines.length < 2) return null;

  // Split a CSV line respecting simple quoting (no embedded commas in quotes).
  const splitLine = (line: string): string[] => {
    const out: string[] = [];
    let cur = '';
    let inQuote = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (ch === '"') {
        inQuote = !inQuote;
      } else if (ch === ',' && !inQuote) {
        out.push(cur);
        cur = '';
      } else {
        cur += ch;
      }
    }
    out.push(cur);
    return out.map((s) => s.trim());
  };

  const header = splitLine(lines[0]);
  if (header.length === 0) return null;

  const data: Record<string, number[]> = {};
  for (const col of header) data[col] = [];

  for (let i = 1; i < lines.length; i++) {
    const cells = splitLine(lines[i]);
    for (let c = 0; c < header.length; c++) {
      const raw = cells[c] ?? '';
      const num = Number(raw);
      data[header[c]].push(Number.isFinite(num) ? num : NaN);
    }
  }

  const numRows = data[header[0]].length;
  return {
    integrationVar: header[0],
    columns: header,
    data,
    numRows,
    csvName,
  };
}

/** Build a CSV string from an IntegralTableData. */
function tableToCSV(table: IntegralTableData): string {
  const lines: string[] = [table.columns.join(',')];
  for (let r = 0; r < table.numRows; r++) {
    lines.push(table.columns.map((c) => table.data[c]?.[r] ?? '').join(','));
  }
  return lines.join('\n') + '\n';
}

// ============================================================================
// Main component
// ============================================================================

export default function IntegralTable() {
  const integralTable = useModelStore((s) => s.integralTable);
  const integralCSV = useModelStore((s) => s.integralCSV);

  // Resolve the active table: prefer the live columnar JSON; fall back to
  // parsing the restored CSV (bundle round-trip).
  const table = useMemo<IntegralTableData | null>(() => {
    if (integralTable && integralTable.numRows > 0) return integralTable;
    if (integralCSV) {
      const name = integralTable?.csvName ?? 'model-integral.csv';
      return parseCSVToTable(integralCSV, name);
    }
    return null;
  }, [integralTable, integralCSV]);

  // Plot Y-variable selector — default to the first non-integration column.
  const [plotYVar, setPlotYVar] = useState('');

  const plotCandidates = useMemo(() => {
    if (!table) return [];
    return table.columns.filter((c) => c !== table.integrationVar);
  }, [table]);

  // Keep the plot selector valid as the table changes.
  useEffect(() => {
    if (plotCandidates.length > 0 && !plotCandidates.includes(plotYVar)) {
      setPlotYVar(plotCandidates[0]);
    }
  }, [plotCandidates, plotYVar]);

  const plotData = useMemo<{ data: Plotly.Data[]; layout: Partial<Plotly.Layout> }>(() => {
    if (!table || !plotYVar) return { data: [], layout: {} };
    const x = table.data[table.integrationVar] ?? [];
    const y = table.data[plotYVar] ?? [];
    const n = Math.min(x.length, y.length);
    return {
      data: [{
        x: x.slice(0, n),
        y: y.slice(0, n),
        type: 'scatter' as const,
        mode: 'lines+markers' as const,
        marker: { size: 5, color: '#3b82f6' },
        line: { color: '#3b82f6', width: 2 },
        name: plotYVar,
      }],
      layout: {
        xaxis: { title: { text: table.integrationVar }, gridcolor: '#e5e7eb' },
        yaxis: { title: { text: plotYVar }, gridcolor: '#e5e7eb' },
        margin: { t: 20, r: 20, b: 45, l: 60 },
        paper_bgcolor: 'transparent',
        plot_bgcolor: 'transparent',
        font: { family: 'system-ui, sans-serif', size: 11 },
      },
    };
  }, [table, plotYVar]);

  // ---- Export CSV ----
  const handleExport = useCallback(() => {
    if (!table) return;
    const csv = tableToCSV(table);
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = table.csvName || 'integral.csv';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, [table]);

  // ---- Empty state ----
  if (!table) {
    return (
      <div className="integral-table-panel">
        <div className="parametric-section-title">Integral Table</div>
        <div className="parametric-notice">
          No integral trajectory available. Write an equation-based integral model
          (using the <code>INTEGRAL(integrand, t, t0, tf)</code> function and a{' '}
          <code>$IntegralTable</code> directive) and solve it to populate this tab.
        </div>
      </div>
    );
  }

  const stepInfo = [
    table.totalSteps !== undefined && `${table.totalSteps} steps`,
    table.rejectedSteps !== undefined && table.rejectedSteps > 0
      ? `${table.rejectedSteps} rejected` : null,
    table.csvName,
  ].filter(Boolean).join(' · ');

  return (
    <div className="integral-table-panel">
      <div className="integral-header">
        <div className="parametric-section-title">
          Integral Table — {table.numRows} rows
          {stepInfo && <span className="integral-meta"> · {stepInfo}</span>}
        </div>
        <button
          className="parametric-run-btn"
          onClick={handleExport}
          title={`Export the trajectory as ${table.csvName || 'integral.csv'}`}
        >
          Export CSV
        </button>
      </div>

      {/* Plot */}
      {plotCandidates.length > 0 && (
        <>
          <div className="parametric-plot-controls">
            <label>Y variable:</label>
            <select
              className="sweep-select"
              value={plotYVar}
              onChange={(e) => setPlotYVar(e.target.value)}
            >
              {plotCandidates.map((n) => (
                <option key={n} value={n}>{n}</option>
              ))}
            </select>
          </div>
          {plotData.data.length > 0 && (
            <div className="parametric-plot">
              <Plot
                data={plotData.data}
                layout={plotData.layout}
                config={{ responsive: true, displaylogo: false }}
                useResizeHandler
                style={{ width: '100%', height: 240 }}
              />
            </div>
          )}
        </>
      )}

      {/* Columnar table (integration variable first) */}
      <div className="parametric-table-container">
        <table className="variable-table parametric-table">
          <thead>
            <tr>
              {table.columns.map((c) => (
                <th
                  key={c}
                  className={c === table.integrationVar ? 'sweep-col' : ''}
                >
                  {c}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {Array.from({ length: table.numRows }, (_, r) => (
              <tr key={r}>
                {table.columns.map((c) => (
                  <td
                    key={c}
                    className={c === table.integrationVar ? 'sweep-col' : 'value-cell'}
                  >
                    {fmt(table.data[c]?.[r] ?? NaN)}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
