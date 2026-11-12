/**
 * DiagramOverlays — the overlays of a diagram definition (docs/coolplot_integration.md §4.1): the
 * inferred states (the states of the solved model as the server groups its variables, drawn as numbered
 * markers) and the arrays (state i from X[i] and Y[i], joined in index order, closed on request, with
 * paths along the processes or straight). Both are cycle definitions; the server reads their variables
 * from the last solve.
 */

import type { CycleDefinition, CycleStateDefinition, DiagramDefinition } from '../../api/types';
import {
  STATE_KEYS,
  arraysCycle,
  isArraysCycle,
  isInferredCycle,
  newInferredCycle,
  overlayComplete,
  overlayOfCycle,
  type ArraysOverlay,
} from './definition';
import { symbolOf } from './options';

interface Props {
  def: DiagramDefinition;
  onChange: (def: DiagramDefinition) => void;
  /** A solve succeeded: its variables can be drawn */
  solved: boolean;
  /** The inferred states of the definition's fluid in the last solve; null while unknown */
  inferredStates: CycleStateDefinition[] | null;
  /** The arrays of the last solve and their indices */
  arrays: Map<string, number[]>;
  /** The property key of an array, when the variable inference found it */
  arrayKeys: Map<string, string>;
  /** The arrays overlay while it names no two arrays (it is then not in the definition) */
  draft: ArraysOverlay;
  setDraft: (draft: ArraysOverlay) => void;
}

/** The cycles with `cycle` in place of the one `isIt` finds (appended if none; removed if `cycle` is null). */
function replaceCycle(cycles: CycleDefinition[], isIt: (c: CycleDefinition) => boolean, cycle: CycleDefinition | null) {
  const index = cycles.findIndex(isIt);
  if (index < 0) return cycle ? [...cycles, cycle] : cycles;
  return cycle ? cycles.map((c, i) => (i === index ? cycle : c)) : cycles.filter((_, i) => i !== index);
}

export default function DiagramOverlays({
  def, onChange, solved, inferredStates, arrays, arrayKeys, draft, setDraft,
}: Props) {
  const inferred = def.cycles.find(isInferredCycle);
  const storedArrays = def.cycles.find(isArraysCycle);
  // A complete overlay is always in the definition: a complete draft without it belongs to a model
  // left since (New, Open, Back), and shows off, its arrays kept for the next time
  const overlay = storedArrays
    ? overlayOfCycle(storedArrays)
    : overlayComplete(draft) ? { ...draft, on: false } : draft;

  const setInferred = (on: boolean) =>
    onChange({ ...def, cycles: replaceCycle(def.cycles, isInferredCycle, on ? newInferredCycle(inferredStates ?? []) : null) });

  const setOverlay = (next: ArraysOverlay) => {
    setDraft(next);
    if (overlayComplete(next)) {
      onChange({ ...def, cycles: replaceCycle(def.cycles, isArraysCycle, arraysCycle(next, arrays, storedArrays)) });
    } else if (storedArrays) {
      onChange({ ...def, cycles: replaceCycle(def.cycles, isArraysCycle, null) });
    }
  };

  const chooseArray = (i: 0 | 1, name: string) => {
    const next: ArraysOverlay = { ...overlay, arrays: [...overlay.arrays] as [string, string], keys: [...overlay.keys] as [string, string] };
    next.arrays[i] = name;
    const key = arrayKeys.get(name);
    if (key && key !== next.keys[1 - i]) next.keys[i] = key;
    setOverlay(next);
  };

  const chooseKey = (i: 0 | 1, key: string) => {
    const keys = [...overlay.keys] as [string, string];
    keys[i] = key;
    setOverlay({ ...overlay, keys });
  };

  const arrayNames = [...arrays.keys()].sort((a, b) => a.localeCompare(b));
  for (const name of overlay.arrays) if (name && !arrays.has(name)) arrayNames.push(name);

  const stateCount = inferred?.states?.length ?? 0;

  return (
    <div className="diagram-row">
      <span className="diagram-row-title">Overlays</span>
      <label
        className="diagram-field"
        title="States of the solved model: the variables of the fluid that the same property calls share"
      >
        <input
          type="checkbox"
          checked={!!inferred}
          disabled={!solved && !inferred}
          onChange={(e) => setInferred(e.target.checked)}
        />
        inferred states{inferred && solved ? ` (${stateCount})` : ''}
      </label>

      <label className="diagram-field" title="A path through the elements of two arrays: state i from X[i] and Y[i]">
        <input
          type="checkbox"
          checked={overlay.on}
          disabled={!solved && !overlay.on}
          onChange={(e) => setOverlay({ ...overlay, on: e.target.checked })}
        />
        arrays
      </label>
      {overlay.on && ([0, 1] as const).map((i) => (
        <span key={i} className="diagram-field">
          <select
            className="diagram-select"
            value={overlay.keys[i]}
            title="Property held by the array"
            onChange={(e) => chooseKey(i, e.target.value)}
          >
            {STATE_KEYS.map((k) => <option key={k} value={k}>{symbolOf(k)}</option>)}
          </select>
          <select
            className="diagram-select"
            value={overlay.arrays[i]}
            title="Array"
            onChange={(e) => chooseArray(i, e.target.value)}
          >
            <option value="">array…</option>
            {arrayNames.map((name) => (
              <option key={name} value={name}>
                {name}{arrays.has(name) ? ` (${arrays.get(name)!.length})` : ''}
              </option>
            ))}
          </select>
        </span>
      ))}
      {overlay.on && (
        <select
          className="diagram-select"
          value={overlay.paths}
          title="Paths between the states: along the processes (auto) or straight"
          onChange={(e) => setOverlay({ ...overlay, paths: e.target.value as 'auto' | 'straight' })}
        >
          <option value="auto">auto paths</option>
          <option value="straight">straight</option>
        </select>
      )}
      {overlay.on && (
        <label className="diagram-field" title="Join the last state to the first">
          <input
            type="checkbox"
            checked={overlay.closed}
            onChange={(e) => setOverlay({ ...overlay, closed: e.target.checked })}
          />
          close loop
        </label>
      )}
      {overlay.on && overlay.keys[0] === overlay.keys[1] && (
        <span className="diagram-error">two different properties needed</span>
      )}
    </div>
  );
}
