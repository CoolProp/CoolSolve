import { useState, useMemo, useCallback } from 'react';
import { useModelStore } from '../stores/modelStore';
import { api } from '../api/client';

function formatValue(val: number): string {
  if (Number.isInteger(val) && Math.abs(val) < 1e9) return val.toString();
  if (Math.abs(val) < 0.01 || Math.abs(val) >= 1e6) return val.toExponential(6);
  return val.toPrecision(7);
}

/** Parse initials text (name=value # comment) into a Map */
function parseInitials(text: string): Map<string, string> {
  const map = new Map<string, string>();
  for (const line of text.split('\n')) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith('#')) continue;
    const eqIdx = trimmed.indexOf('=');
    if (eqIdx < 0) continue;
    const name = trimmed.substring(0, eqIdx).trim();
    // Value is everything after = up to # comment
    let rest = trimmed.substring(eqIdx + 1);
    const hashIdx = rest.indexOf('#');
    if (hashIdx >= 0) rest = rest.substring(0, hashIdx);
    map.set(name, rest.trim());
  }
  return map;
}

/** Rebuild initials text from the map */
function serializeInitials(map: Map<string, string>): string {
  const lines: string[] = [];
  for (const [name, value] of map) {
    lines.push(`${name}=${value}`);
  }
  return lines.join('\n') + '\n';
}

export default function VariableTable() {
  const [filter, setFilter] = useState('');
  const [editingCell, setEditingCell] = useState<string | null>(null);
  const [editValue, setEditValue] = useState('');
  const lastResult = useModelStore((s) => s.lastResult);
  const initials = useModelStore((s) => s.initials);
  const setInitials = useModelStore((s) => s.setInitials);

  const initialsMap = useMemo(() => parseInitials(initials), [initials]);

  const variables = useMemo(() => {
    if (!lastResult) {
      // If no solve result, show initials only
      if (initialsMap.size === 0) return [];
      const vars: { name: string; value: string; initial: string; isArray: boolean }[] = [];
      for (const [name, val] of initialsMap) {
        vars.push({ name, value: '', initial: val, isArray: name.includes('[') });
      }
      vars.sort((a, b) => a.name.localeCompare(b.name, undefined, { sensitivity: 'base' }));
      return vars;
    }

    const vars: { name: string; value: string; initial: string; isArray: boolean; isString: boolean }[] = [];

    // Numeric variables
    for (const [name, val] of Object.entries(lastResult.variables)) {
      vars.push({
        name,
        value: formatValue(val),
        initial: initialsMap.get(name) ?? '',
        isArray: name.includes('['),
        isString: false,
      });
    }

    // String variables
    for (const [name, val] of Object.entries(lastResult.stringVariables)) {
      vars.push({
        name,
        value: `'${val}'`,
        initial: '',
        isArray: false,
        isString: true,
      });
    }

    // Sort: scalars first alphabetically, then arrays
    vars.sort((a, b) => {
      if (a.isArray !== b.isArray) return a.isArray ? 1 : -1;
      return a.name.localeCompare(b.name, undefined, { sensitivity: 'base' });
    });

    return vars;
  }, [lastResult, initialsMap]);

  const filtered = useMemo(() => {
    if (!filter) return variables;
    const lower = filter.toLowerCase();
    return variables.filter((v) => v.name.toLowerCase().includes(lower));
  }, [variables, filter]);

  const handleInitialClick = useCallback((name: string, current: string) => {
    setEditingCell(name);
    setEditValue(current);
  }, []);

  const commitEdit = useCallback((name: string) => {
    setEditingCell(null);
    const newMap = new Map(initialsMap);
    const trimmed = editValue.trim();
    if (trimmed === '') {
      newMap.delete(name);
    } else {
      newMap.set(name, trimmed);
    }
    const newText = serializeInitials(newMap);
    setInitials(newText);
    api.putInitials(newText).catch(() => {});
  }, [editValue, initialsMap, setInitials]);

  const handleKeyDown = useCallback((e: React.KeyboardEvent, name: string) => {
    if (e.key === 'Enter') {
      commitEdit(name);
    } else if (e.key === 'Escape') {
      setEditingCell(null);
    }
  }, [commitEdit]);

  return (
    <>
      <div className="search-filter">
        <input
          type="text"
          placeholder="Filter variables..."
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
        />
      </div>
      <div className="variable-table-container">
        {variables.length === 0 ? (
          <div style={{ padding: 16, color: 'var(--text-muted)', fontSize: 13 }}>
            No solve result yet. Click <strong>Solve</strong> to compute variable values.
          </div>
        ) : (
          <table className="variable-table">
            <thead>
              <tr>
                <th>Name</th>
                <th style={{ textAlign: 'right' }}>Initial</th>
                <th style={{ textAlign: 'right' }}>Value</th>
              </tr>
            </thead>
            <tbody>
              {filtered.map((v) => (
                <tr key={v.name}>
                  <td>{v.name}</td>
                  <td
                    className="value-cell initial-cell"
                    onClick={() => handleInitialClick(v.name, v.initial)}
                  >
                    {editingCell === v.name ? (
                      <input
                        className="initial-input"
                        type="text"
                        value={editValue}
                        onChange={(e) => setEditValue(e.target.value)}
                        onBlur={() => commitEdit(v.name)}
                        onKeyDown={(e) => handleKeyDown(e, v.name)}
                        autoFocus
                      />
                    ) : (
                      <span className="initial-value">{v.initial || '—'}</span>
                    )}
                  </td>
                  <td className="value-cell solved">{v.value}</td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
        {variables.length > 0 && (
          <div style={{ padding: '4px 8px', fontSize: 11, color: 'var(--text-muted)' }}>
            {filtered.length} of {variables.length} variables
          </div>
        )}
      </div>
    </>
  );
}
