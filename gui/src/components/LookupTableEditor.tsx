import React, { useState, useCallback, useRef } from 'react';
import { api } from '../api/client';
import { useModelStore } from '../stores/modelStore';

// ============================================================================
// Helpers
// ============================================================================

function parseCSV(csv: string): string[][] {
  const rows: string[][] = [];
  const lines = csv.split(/\r?\n/);
  for (const line of lines) {
    if (line.trim() === '') continue;
    // Simple CSV parser: handles quoted fields
    const cells: string[] = [];
    let current = '';
    let inQuotes = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (ch === '"') {
        inQuotes = !inQuotes;
      } else if (ch === ',' && !inQuotes) {
        cells.push(current.trim());
        current = '';
      } else {
        current += ch;
      }
    }
    cells.push(current.trim());
    rows.push(cells);
  }
  return rows;
}

function serializeCSV(rows: string[][]): string {
  return rows
    .map((row) =>
      row
        .map((cell) => (cell.includes(',') || cell.includes('"') ? `"${cell.replace(/"/g, '""')}"` : cell))
        .join(',')
    )
    .join('\n');
}

// ============================================================================
// TableGrid: editable grid for one table
// ============================================================================

interface TableGridProps {
  tableName: string;
  onClose: () => void;
  onSaved: () => void;
}

function TableGrid({ tableName, onClose, onSaved }: TableGridProps) {
  const lookupTableCSVs = useModelStore((s) => s.lookupTableCSVs);
  const setLookupTableCSV = useModelStore((s) => s.setLookupTableCSV);

  const rawCSV = lookupTableCSVs[tableName] ?? '';
  const [rows, setRows] = useState<string[][]>(() => {
    const parsed = parseCSV(rawCSV);
    // Ensure at least 1 header row + 2 data rows and 2 columns
    if (parsed.length === 0) return [['Column1', 'Column2'], ['', ''], ['', '']];
    return parsed;
  });
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const setCell = (r: number, c: number, val: string) => {
    setRows((prev) => {
      const next = prev.map((row) => [...row]);
      next[r][c] = val;
      return next;
    });
  };

  const addRow = () =>
    setRows((prev) => [...prev, Array(prev[0]?.length ?? 1).fill('')]);

  const addColumn = () =>
    setRows((prev) => prev.map((row, i) => [...row, i === 0 ? `Column${row.length + 1}` : '']));

  const deleteRow = (r: number) => {
    if (r === 0 || rows.length <= 2) return; // keep header + 1 data row
    setRows((prev) => prev.filter((_, i) => i !== r));
  };

  const deleteColumn = (c: number) => {
    if (rows[0]?.length <= 1) return;
    setRows((prev) => prev.map((row) => row.filter((_, i) => i !== c)));
  };

  const handleSave = async () => {
    setSaving(true);
    setError(null);
    try {
      const csv = serializeCSV(rows);
      await api.putTable(tableName, csv);
      setLookupTableCSV(tableName, csv);
      onSaved();
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setSaving(false);
    }
  };

  const handleImportCSV = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      const text = ev.target?.result as string;
      const parsed = parseCSV(text);
      if (parsed.length > 0) setRows(parsed);
    };
    reader.readAsText(file);
    e.target.value = '';
  };

  const handleExportCSV = () => {
    const csv = serializeCSV(rows);
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${tableName}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const numCols = rows[0]?.length ?? 0;

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%', gap: 6 }}>
      {/* Toolbar */}
      <div style={{ display: 'flex', gap: 6, alignItems: 'center', flexShrink: 0, padding: '4px 0' }}>
        <span style={{ fontWeight: 600, fontSize: 13 }}>{tableName}</span>
        <button className="toolbar-btn" onClick={addRow} title="Add row">+ Row</button>
        <button className="toolbar-btn" onClick={addColumn} title="Add column">+ Col</button>
        <label className="toolbar-btn" style={{ cursor: 'pointer' }} title="Import CSV">
          Import CSV
          <input type="file" accept=".csv" style={{ display: 'none' }} onChange={handleImportCSV} />
        </label>
        <button className="toolbar-btn" onClick={handleExportCSV} title="Export CSV">Export CSV</button>
        <div style={{ flex: 1 }} />
        {error && <span style={{ color: '#f87171', fontSize: 12 }}>{error}</span>}
        <button className="toolbar-btn" onClick={handleSave} disabled={saving}>
          {saving ? 'Saving…' : 'Save'}
        </button>
        <button className="toolbar-btn" onClick={onClose}>✕</button>
      </div>

      {/* Grid */}
      <div style={{ flex: 1, overflow: 'auto' }}>
        <table style={{ borderCollapse: 'collapse', fontSize: 12, tableLayout: 'fixed', minWidth: '100%' }}>
          <thead>
            <tr>
              <th style={{ width: 28 }} />
              {Array.from({ length: numCols }).map((_, c) => (
                <th key={c} style={{ border: '1px solid #3a3a3a', padding: '2px 4px', background: '#1e1e2e' }}>
                  <div style={{ display: 'flex', alignItems: 'center', gap: 2 }}>
                    <input
                      value={rows[0]?.[c] ?? ''}
                      onChange={(e) => setCell(0, c, e.target.value)}
                      style={{ flex: 1, background: 'transparent', border: 'none', color: '#cdd6f4', fontWeight: 600, fontSize: 12, outline: 'none' }}
                    />
                    {numCols > 1 && (
                      <button
                        onClick={() => deleteColumn(c)}
                        style={{ background: 'none', border: 'none', cursor: 'pointer', color: '#585b70', fontSize: 11, padding: 0, lineHeight: 1 }}
                        title="Delete column"
                      >×</button>
                    )}
                  </div>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.slice(1).map((row, ri) => (
              <tr key={ri + 1}>
                <td style={{ textAlign: 'center', color: '#6c7086', fontSize: 11, border: '1px solid #2a2a3a' }}>
                  <button
                    onClick={() => deleteRow(ri + 1)}
                    style={{ background: 'none', border: 'none', cursor: 'pointer', color: '#585b70', fontSize: 12 }}
                    title="Delete row"
                  >×</button>
                </td>
                {Array.from({ length: numCols }).map((_, c) => (
                  <td key={c} style={{ border: '1px solid #2a2a3a', padding: 0 }}>
                    <input
                      value={row[c] ?? ''}
                      onChange={(e) => setCell(ri + 1, c, e.target.value)}
                      style={{ width: '100%', background: 'transparent', border: 'none', color: '#cdd6f4', fontSize: 12, padding: '2px 4px', outline: 'none', boxSizing: 'border-box' }}
                    />
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

// ============================================================================
// LookupTableEditor: list of tables + optional grid view
// ============================================================================

export default function LookupTableEditor() {
  const lookupTables = useModelStore((s) => s.lookupTables);
  const lookupTableCSVs = useModelStore((s) => s.lookupTableCSVs);
  const setLookupTables = useModelStore((s) => s.setLookupTables);
  const setLookupTableCSV = useModelStore((s) => s.setLookupTableCSV);
  const removeLookupTable = useModelStore((s) => s.removeLookupTable);

  const [selectedTable, setSelectedTable] = useState<string | null>(null);
  const [newTableName, setNewTableName] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const newTableRef = useRef<HTMLInputElement>(null);

  // Refresh table list from server
  const refreshTables = useCallback(async () => {
    try {
      const resp = await api.getTables();
      setLookupTables(resp.tables);
    } catch {
      // Non-fatal; list may be stale
    }
  }, [setLookupTables]);

  // Open a table (fetch CSV if not cached)
  const openTable = useCallback(async (name: string) => {
    if (!lookupTableCSVs[name]) {
      try {
        const csv = await api.getTableCSV(name);
        setLookupTableCSV(name, csv);
      } catch (e) {
        setError(e instanceof Error ? e.message : String(e));
        return;
      }
    }
    setSelectedTable(name);
  }, [lookupTableCSVs, setLookupTableCSV]);

  const handleCreateTable = async () => {
    const name = newTableName.trim();
    if (!name) return;
    const defaultCSV = 'Column1,Column2\n1,\n2,\n';
    setLoading(true);
    setError(null);
    try {
      await api.putTable(name, defaultCSV);
      setLookupTableCSV(name, defaultCSV);
      await refreshTables();
      setNewTableName('');
      setSelectedTable(name);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setLoading(false);
    }
  };

  const handleDeleteTable = async (name: string, e: React.MouseEvent) => {
    e.stopPropagation();
    if (!confirm(`Delete table "${name}"?`)) return;
    try {
      await api.deleteTable(name);
      removeLookupTable(name);
      if (selectedTable === name) setSelectedTable(null);
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    }
  };

  // If editing a table, show the grid
  if (selectedTable !== null) {
    return (
      <div style={{ padding: 8, height: '100%', boxSizing: 'border-box', display: 'flex', flexDirection: 'column' }}>
        <TableGrid
          tableName={selectedTable}
          onClose={() => setSelectedTable(null)}
          onSaved={() => { refreshTables(); }}
        />
      </div>
    );
  }

  return (
    <div style={{ padding: 8, display: 'flex', flexDirection: 'column', gap: 8, height: '100%', boxSizing: 'border-box' }}>
      <div style={{ fontSize: 12, color: '#a6adc8', flexShrink: 0 }}>
        Lookup tables are CSV files loaded from the model's directory using the
        naming convention <code style={{ color: '#cba6f7' }}>&lt;modelname&gt;-&lt;tablename&gt;.csv</code>.{' '}
        <code style={{ color: '#cba6f7' }}>LOOKUP('table', row, col)</code> — direct cell access;{' '}
        <code style={{ color: '#cba6f7' }}>INTERPOLATE('table', 'xcol', 'ycol', x)</code> — 1-D interpolation.
      </div>

      {/* Create new table */}
      <div style={{ display: 'flex', gap: 6, alignItems: 'center', flexShrink: 0 }}>
        <input
          ref={newTableRef}
          value={newTableName}
          onChange={(e) => setNewTableName(e.target.value)}
          onKeyDown={(e) => e.key === 'Enter' && handleCreateTable()}
          placeholder="New table name…"
          style={{ flex: 1, padding: '3px 6px', background: '#1e1e2e', border: '1px solid #3a3a4a', borderRadius: 3, color: '#cdd6f4', fontSize: 12 }}
        />
        <button className="toolbar-btn" onClick={handleCreateTable} disabled={loading || !newTableName.trim()}>
          Create
        </button>
      </div>

      {error && <div style={{ color: '#f87171', fontSize: 12 }}>{error}</div>}

      {/* Table list */}
      <div style={{ flex: 1, overflow: 'auto' }}>
        {lookupTables.length === 0 ? (
          <div style={{ color: '#585b70', fontSize: 12, paddingTop: 8 }}>
            No lookup tables. Create one above or load a model with .csv files.
          </div>
        ) : (
          <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: 12 }}>
            <thead>
              <tr style={{ color: '#a6adc8', textAlign: 'left', borderBottom: '1px solid #3a3a4a' }}>
                <th style={{ padding: '3px 6px' }}>Name</th>
                <th style={{ padding: '3px 6px' }}>Rows</th>
                <th style={{ padding: '3px 6px' }}>Columns</th>
                <th style={{ width: 40 }} />
              </tr>
            </thead>
            <tbody>
              {lookupTables.map((tbl) => (
                <tr
                  key={tbl.name}
                  onClick={() => openTable(tbl.name)}
                  style={{ cursor: 'pointer', borderBottom: '1px solid #2a2a3a' }}
                  className="lookup-table-row"
                >
                  <td style={{ padding: '4px 6px', color: '#cba6f7' }}>{tbl.name}</td>
                  <td style={{ padding: '4px 6px', color: '#cdd6f4' }}>{tbl.rows}</td>
                  <td style={{ padding: '4px 6px', color: '#a6adc8' }}>{tbl.columns.join(', ')}</td>
                  <td style={{ padding: '4px 6px' }}>
                    <button
                      onClick={(e) => handleDeleteTable(tbl.name, e)}
                      style={{ background: 'none', border: 'none', cursor: 'pointer', color: '#585b70', fontSize: 14 }}
                      title="Delete table"
                    >🗑</button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </div>
    </div>
  );
}
