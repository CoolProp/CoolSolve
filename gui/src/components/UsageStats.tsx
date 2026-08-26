/**
 * UsageStats.tsx — Solve usage statistics dashboard.
 *
 * Standalone page served at the hidden URL `/stats` (no link points to it).
 * It fetches aggregated statistics from GET /api/v1/stats/log, which the
 * backend computes from the GUI solve-attempt log file (see
 * include/coolsolve/usage_log.h). Charts are rendered with Plotly, which is
 * already bundled for the parametric study — no additional dependency.
 */
import { useCallback, useEffect, useMemo, useState } from 'react';
import Plot from './PlotlyChart';
import { api } from '../api/client';
import { useUIStore } from '../stores/uiStore';
import type { UsageStatsResponse } from '../api/types';

// ============================================================================
// Helpers
// ============================================================================

/** Format a duration in ms using human-friendly units. */
function fmtDuration(ms: number): string {
  if (!isFinite(ms) || ms <= 0) return '—';
  if (ms < 1) return `${ms.toFixed(2)} ms`;
  if (ms < 1000) return `${ms.toFixed(0)} ms`;
  if (ms < 60000) return `${(ms / 1000).toFixed(1)} s`;
  return `${(ms / 60000).toFixed(1)} min`;
}

/** Geometric center of a histogram bin (for plotting on a log axis). */
function binCenter(lo: number, hi: number): number {
  return Math.sqrt(lo * hi);
}

// ============================================================================
// Main Component
// ============================================================================

export default function UsageStats() {
  const theme = useUIStore((s) => s.theme);
  const [stats, setStats] = useState<UsageStatsResponse | null>(null);
  const [error, setError] = useState<string>('');
  const [loading, setLoading] = useState(true);

  // This page is rendered outside the main editor layout, so it must apply
  // the theme attribute itself.
  useEffect(() => {
    document.documentElement.setAttribute('data-theme', theme);
  }, [theme]);

  const refresh = useCallback(async () => {
    setLoading(true);
    setError('');
    try {
      setStats(await api.getUsageStats());
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    refresh();
  }, [refresh]);

  const dark = theme === 'dark';

  const successRate = stats && stats.totalAttempts > 0
    ? (100 * stats.successes) / stats.totalAttempts
    : null;

  // --- Chart data ----------------------------------------------------------

  const dailyData = useMemo((): Plotly.Data[] => {
    if (!stats?.daily.length) return [];
    const dates = stats.daily.map((d) => d.date);
    const successes = stats.daily.map((d) => d.successes);
    const other = stats.daily.map((d) => d.attempts - d.successes);
    return [
      { type: 'bar', name: 'Success', x: dates, y: successes, marker: { color: '#4ade80' } },
      {
        type: 'bar', name: 'Failed', x: dates, y: other,
        marker: { color: '#f87171' },
      },
    ];
  }, [stats]);

  const outcomeData = useMemo((): Plotly.Data[] => {
    if (!stats) return [];
    const o = stats.outcomes;
    return [
      {
        type: 'pie',
        labels: ['Success', 'Failed', 'Parse error'],
        values: [o.success ?? 0, o.failed ?? 0, o.parse_error ?? 0],
        marker: { colors: ['#4ade80', '#f87171', '#fbbf24'] },
        textinfo: 'value',
      },
    ];
  }, [stats]);

  const durationData = useMemo((): Plotly.Data[] => {
    if (!stats?.durationHistogram.counts.length) return [];
    const edges = stats.durationHistogram.edges;
    const counts = stats.durationHistogram.counts;
    const centers = counts.map((_, i) => binCenter(edges[i], edges[i + 1]));
    const labels = counts.map((_, i) =>
      `${fmtDuration(edges[i])} – ${fmtDuration(edges[i + 1])}`,
    );
    return [
      {
        type: 'bar',
        x: centers,
        y: counts,
        text: labels,
        hovertemplate: '%{text}<br>%{y} solves<extra></extra>',
        marker: { color: dark ? '#4dabf7' : '#0066cc' },
      },
    ];
  }, [stats, dark]);

  const rankingData = (entries: { name: string; count: number }[]): Plotly.Data[] => [
    {
      type: 'bar',
      orientation: 'h',
      // Reverse so the most frequent entry appears at the top
      y: entries.map((e) => e.name).reverse(),
      x: entries.map((e) => e.count).reverse(),
      marker: { color: dark ? '#4dabf7' : '#0066cc' },
    },
  ];

  const baseLayout = useMemo(
    (): Partial<Plotly.Layout> => ({
      margin: { t: 12, r: 16, b: 40, l: 48 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: {
        size: 11,
        color: dark ? '#d4d4d4' : '#1a1a1a',
      },
      bargap: 0.15,
      xaxis: { gridcolor: dark ? '#3c3c3c' : '#e0e0e0', automargin: true },
      yaxis: { gridcolor: dark ? '#3c3c3c' : '#e0e0e0', automargin: true },
    }),
    [dark],
  );

  const dailyLayout = useMemo(
    (): Partial<Plotly.Layout> => ({
      ...baseLayout,
      barmode: 'stack',
      legend: { orientation: 'h', y: 1.12 },
      xaxis: { ...baseLayout.xaxis as NonNullable<Partial<Plotly.Layout>['xaxis']>, type: 'date' },
      yaxis: {
        ...(baseLayout.yaxis as NonNullable<Partial<Plotly.Layout>['yaxis']>),
        title: { text: 'Attempts' },
        rangemode: 'nonnegative',
      },
    }),
    [baseLayout],
  );

  const durationLayout = useMemo(
    (): Partial<Plotly.Layout> => ({
      ...baseLayout,
      xaxis: {
        ...(baseLayout.xaxis as NonNullable<Partial<Plotly.Layout>['xaxis']>),
        type: 'log',
        title: { text: 'Duration (ms)' },
      },
      yaxis: {
        ...(baseLayout.yaxis as NonNullable<Partial<Plotly.Layout>['yaxis']>),
        title: { text: 'Solves' },
        rangemode: 'nonnegative',
      },
    }),
    [baseLayout],
  );

  const rankingLayout = useMemo(
    (): Partial<Plotly.Layout> => ({
      ...baseLayout,
      margin: { t: 8, r: 16, b: 32, l: 150 },
      xaxis: {
        ...(baseLayout.xaxis as NonNullable<Partial<Plotly.Layout>['xaxis']>),
        rangemode: 'nonnegative',
      },
    }),
    [baseLayout],
  );

  // --- Render ---------------------------------------------------------------

  return (
    <div className="app stats-page" data-theme={theme}>
      <div className="stats-header">
        <h1>Solve Usage Statistics</h1>
        <button className="btn" onClick={refresh} disabled={loading}>
          {loading ? 'Loading…' : 'Refresh'}
        </button>
      </div>

      {error && (
        <div className="stats-empty">
          <p>Failed to load statistics: {error}</p>
        </div>
      )}

      {!error && stats && !stats.valid && (
        <div className="stats-empty">
          <p>No usage log found yet.</p>
          <p className="stats-empty-hint">
            Statistics appear here once at least one model has been solved
            through the web interface.
          </p>
        </div>
      )}

      {!error && stats?.valid && (
        <>
          {/* Summary cards */}
          <div className="stats-cards">
            <div className="stat-card">
              <div className="stat-value">{stats.totalAttempts}</div>
              <div className="stat-label">Total attempts</div>
            </div>
            <div className="stat-card">
              <div className="stat-value">{successRate === null ? '—' : `${successRate.toFixed(1)}%`}</div>
              <div className="stat-label">Success rate</div>
            </div>
            <div className="stat-card">
              <div className="stat-value">{fmtDuration(stats.medianMs)}</div>
              <div className="stat-label">Median duration</div>
            </div>
            <div className="stat-card">
              <div className="stat-value">{fmtDuration(stats.p95Ms)}</div>
              <div className="stat-label">P95 duration</div>
            </div>
            <div className="stat-card">
              <div className="stat-value">{stats.tryHarderAttempts}</div>
              <div className="stat-label">Try Harder runs</div>
            </div>
            <div className="stat-card">
              <div className="stat-value">{stats.uniqueIps}</div>
              <div className="stat-label">Unique clients</div>
            </div>
          </div>

          {/* Charts */}
          <div className="stats-charts">
            <div className="stats-chart-card stats-chart-wide">
              <h2>Attempts per day</h2>
              <Plot data={dailyData} layout={dailyLayout}
                config={{ responsive: true, displaylogo: false }}
                useResizeHandler style={{ width: '100%', height: 300 }} />
            </div>

            <div className="stats-chart-card">
              <h2>Outcomes</h2>
              <Plot data={outcomeData} layout={{ ...baseLayout, showlegend: false }}
                config={{ responsive: true, displaylogo: false }}
                useResizeHandler style={{ width: '100%', height: 300 }} />
            </div>

            <div className="stats-chart-card stats-chart-wide">
              <h2>Solve duration distribution</h2>
              <Plot data={durationData} layout={durationLayout}
                config={{ responsive: true, displaylogo: false }}
                useResizeHandler style={{ width: '100%', height: 280 }} />
            </div>

            <div className="stats-chart-card">
              <h2>Most solved models</h2>
              {stats.topModels.length > 0 ? (
                <Plot data={rankingData(stats.topModels)} layout={rankingLayout}
                  config={{ responsive: true, displaylogo: false }}
                  useResizeHandler style={{ width: '100%', height: 280 }} />
              ) : (
                <p className="stats-chart-empty">No named models yet</p>
              )}
            </div>

            <div className="stats-chart-card">
              <h2>Clients by attempts</h2>
              {stats.topIps.length > 0 ? (
                <Plot data={rankingData(stats.topIps)} layout={rankingLayout}
                  config={{ responsive: true, displaylogo: false }}
                  useResizeHandler style={{ width: '100%', height: 280 }} />
              ) : (
                <p className="stats-chart-empty">No data</p>
              )}
            </div>
          </div>
        </>
      )}
    </div>
  );
}
