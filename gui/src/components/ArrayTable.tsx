import { useState, useMemo } from 'react';
import { useModelStore } from '../stores/modelStore';

function formatValue(val: number): string {
  if (Number.isInteger(val) && Math.abs(val) < 1e9) return val.toString();
  if (Math.abs(val) < 0.01 || Math.abs(val) >= 1e6) return val.toExponential(6);
  return val.toPrecision(7);
}

interface ArrayEntry {
  baseName: string;
  index: number;
  value: number;
}

export default function ArrayTable() {
  const [filter, setFilter] = useState('');
  const lastResult = useModelStore((s) => s.lastResult);

  const { baseNames, indices, grid } = useMemo(() => {
    if (!lastResult) return { baseNames: [] as string[], indices: [] as number[], grid: new Map<string, Map<number, number>>() };

    const entries: ArrayEntry[] = [];
    for (const [name, val] of Object.entries(lastResult.variables)) {
      const match = name.match(/^(.+)\[(\d+)\]$/);
      if (match) {
        entries.push({ baseName: match[1], index: parseInt(match[2], 10), value: val });
      }
    }

    // Collect unique base names and indices
    const baseSet = new Set<string>();
    const indexSet = new Set<number>();
    const grid = new Map<string, Map<number, number>>();

    for (const e of entries) {
      baseSet.add(e.baseName);
      indexSet.add(e.index);
      if (!grid.has(e.baseName)) grid.set(e.baseName, new Map());
      grid.get(e.baseName)!.set(e.index, e.value);
    }

    const baseNames = Array.from(baseSet).sort((a, b) =>
      a.localeCompare(b, undefined, { sensitivity: 'base' })
    );
    const indices = Array.from(indexSet).sort((a, b) => a - b);

    return { baseNames, indices, grid };
  }, [lastResult]);

  const filteredBaseNames = useMemo(() => {
    if (!filter) return baseNames;
    const lower = filter.toLowerCase();
    return baseNames.filter((n) => n.toLowerCase().includes(lower));
  }, [baseNames, filter]);

  if (!lastResult || baseNames.length === 0) {
    return (
      <div style={{ padding: 16, color: 'var(--text-muted)', fontSize: 13 }}>
        {lastResult
          ? 'No array variables found in the solution.'
          : 'No solve result yet. Click Solve to compute variable values.'}
      </div>
    );
  }

  return (
    <>
      <div className="search-filter">
        <input
          type="text"
          placeholder="Filter array names..."
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
        />
      </div>
      <div className="variable-table-container">
        <table className="variable-table array-table">
          <thead>
            <tr>
              <th>Index</th>
              {filteredBaseNames.map((name) => (
                <th key={name} style={{ textAlign: 'right' }}>{name}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {indices.map((idx) => (
              <tr key={idx}>
                <td style={{ fontWeight: 600 }}>{idx}</td>
                {filteredBaseNames.map((name) => {
                  const val = grid.get(name)?.get(idx);
                  return (
                    <td key={name} className="value-cell solved">
                      {val !== undefined ? formatValue(val) : ''}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
        <div style={{ padding: '4px 8px', fontSize: 11, color: 'var(--text-muted)' }}>
          {filteredBaseNames.length} array{filteredBaseNames.length !== 1 ? 's' : ''} × {indices.length} indices
        </div>
      </div>
    </>
  );
}
