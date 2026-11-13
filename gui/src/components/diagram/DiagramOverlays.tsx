/**
 * DiagramOverlays — the overlays of a diagram definition (docs/coolplot_integration.md §4.1): the
 * inferred states (the states of the solved model as the server groups its variables, drawn as numbered
 * markers) and any number of arrays overlays, one row each: the arrays of the chart's horizontal and
 * vertical quantities (state i from their elements i, joined in index order with paths along the
 * processes, closed on request). All are cycle definitions; the server reads their variables from the
 * last solve. The first arrays overlay switches the inferred states off, which would number the same states
 * a second time (the checkbox brings them back). The tab gives a new overlay its default arrays and maps the overlays to the axes of
 * another diagram type (normalizeDefinition).
 */

import { Plus, X } from 'lucide-react';
import type { CycleDefinition, CycleStateDefinition, DiagramDefinition, DiagramOptionsResponse, DiagramTypeOption } from '../../api/types';
import {
  arraysCycle,
  axisKeys,
  diagramOption,
  isArraysCycle,
  isInferredCycle,
  newArraysCycle,
  newInferredCycle,
  overlayOf,
  type ArraysOverlay,
  type SolvedArrays,
} from './definition';
import { axisSymbol } from './options';

interface Props {
  def: DiagramDefinition;
  options: DiagramOptionsResponse;
  onChange: (def: DiagramDefinition) => void;
  /** A solve succeeded: its variables can be drawn */
  solved: boolean;
  /** The inferred states of the definition's fluid in the last solve; null while unknown */
  inferredStates: CycleStateDefinition[] | null;
  /** The arrays of the last solve */
  arrays: SolvedArrays;
}

/** One arrays overlay: an array per axis, close loop, remove. */
function ArraysRow({ cycle, type, arrays, onChange, onRemove }: {
  cycle: CycleDefinition;
  type: DiagramTypeOption;
  arrays: SolvedArrays;
  onChange: (cycle: CycleDefinition) => void;
  onRemove: () => void;
}) {
  const keys = axisKeys(type);
  const overlay = overlayOf(cycle, keys);
  const set = (change: Partial<ArraysOverlay>) => onChange(arraysCycle({ ...overlay, ...change }, keys, cycle));
  return (
    <div className="diagram-row diagram-arrays">
      {(['x', 'y'] as const).map((axis, i) => {
        const chosen = overlay[axis];
        const names = [...arrays.sizes.keys()];
        if (chosen && !arrays.sizes.has(chosen)) names.push(chosen);
        return (
          <label key={axis} className="diagram-field" title={`Array of the ${axis === 'x' ? 'horizontal' : 'vertical'} axis quantity`}>
            {axisSymbol(type, keys[i])}:
            <select className="diagram-select" value={chosen} onChange={(e) => set(axis === 'x' ? { x: e.target.value } : { y: e.target.value })}>
              {!chosen && <option value="" disabled>array…</option>}
              {names.map((name) => (
                <option key={name} value={name}>
                  {name}{arrays.sizes.has(name) ? ` (${arrays.sizes.get(name)})` : ''}
                </option>
              ))}
            </select>
          </label>
        );
      })}
      <label className="diagram-field" title="Join the last state to the first">
        <input type="checkbox" checked={overlay.closed} onChange={(e) => set({ closed: e.target.checked })} />
        close loop
      </label>
      <button className="diagram-icon-btn" title="Remove these arrays" onClick={onRemove}><X size={12} /></button>
    </div>
  );
}

export default function DiagramOverlays({ def, options, onChange, solved, inferredStates, arrays }: Props) {
  const inferred = def.cycles.find(isInferredCycle);
  const type = diagramOption(options, def.spec.diagram);

  const setInferred = (on: boolean) => {
    const others = def.cycles.filter((c) => !isInferredCycle(c));
    onChange({ ...def, cycles: on ? [...others, newInferredCycle(inferredStates ?? [])] : others });
  };
  const setCycle = (index: number, cycle: CycleDefinition | null) => onChange({
    ...def,
    cycles: cycle ? def.cycles.map((c, i) => (i === index ? cycle : c)) : def.cycles.filter((_, i) => i !== index),
  });

  return (
    <div className="diagram-overlays">
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
          inferred states{inferred && solved ? ` (${inferred.states?.length ?? 0})` : ''}
        </label>
        <button
          className="diagram-link"
          disabled={arrays.sizes.size === 0}
          title={arrays.sizes.size > 0
            ? 'A path through the elements of two arrays, those of the quantities of the axes: state i from their elements i'
            : 'The last solve has no arrays'}
          onClick={() => {
            const first = !def.cycles.some(isArraysCycle);
            const kept = first ? def.cycles.filter((c) => !isInferredCycle(c)) : def.cycles;
            onChange({ ...def, cycles: [...kept, newArraysCycle()] });
          }}
        >
          <Plus size={11} /> arrays
        </button>
      </div>
      {type && def.cycles.map((cycle, i) => isArraysCycle(cycle) && (
        <ArraysRow
          key={i}
          cycle={cycle}
          type={type}
          arrays={arrays}
          onChange={(c) => setCycle(i, c)}
          onRemove={() => setCycle(i, null)}
        />
      ))}
    </div>
  );
}
