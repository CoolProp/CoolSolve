/**
 * DiagramTab — the Diagrams tab: a property diagram of the model drawn by CoolPlot on the server
 * (docs/coolplot_integration.md, phase 1).
 *
 * It edits the first diagram of the model's diagram definitions (`diagrams` in the model store, saved in
 * the session and the ZIP bundle as `diagrams.json`). Until the user changes something, it shows a
 * default chart of the model's first fluid that is not stored. The chart is redrawn by `DiagramView`
 * when the definition, the GUI theme or the solution changes.
 */

import { useCallback, useEffect, useMemo, useState } from 'react';
import { api } from '../../api/client';
import type {
  DiagramDefinition,
  DiagramOptionsResponse,
  DiagramViewSI,
  FluidInfo,
  InferredVariable,
  SolveResponse,
} from '../../api/types';
import { useModelStore } from '../../stores/modelStore';
import { useUIStore } from '../../stores/uiStore';
import { watchDiagramEdits } from '../../utils/diagrams';
import { toErrMsg } from '../../utils/errors';
import DiagramControls from './DiagramControls';
import DiagramOverlays from './DiagramOverlays';
import DiagramView from './DiagramView';
import {
  FALLBACK_FLUID,
  arrayIndices,
  arraysCycle,
  defaultDefinition,
  inferredStates,
  isArraysCycle,
  isInferredCycle,
  newInferredCycle,
  overlayOfCycle,
  propertyKey,
  sameFluid,
  sameJson,
  sceneInput,
  withView,
  type ArraysOverlay,
} from './definition';

export default function DiagramTab() {
  const lastResult = useModelStore((s) => s.lastResult);
  const diagrams = useModelStore((s) => s.diagrams);
  const setDiagramDefinition = useModelStore((s) => s.setDiagramDefinition);
  const theme = useUIStore((s) => s.theme);

  const [options, setOptions] = useState<DiagramOptionsResponse | null>(null);
  const [optionsError, setOptionsError] = useState('');
  const [optionsAttempt, setOptionsAttempt] = useState(0);
  const [fluids, setFluids] = useState<FluidInfo[]>([]);
  const [modelFluids, setModelFluids] = useState<string[]>([]);
  // The inferred variables and the solve they belong to
  const [inferred, setInferred] = useState<{ result: SolveResponse | null; vars: InferredVariable[] }>({ result: null, vars: [] });
  const [arraysDraft, setArraysDraft] = useState<ArraysOverlay>({ on: false, keys: ['h', 'p'], arrays: ['', ''], closed: false });

  // Edits of the definitions go to the session
  useEffect(() => watchDiagramEdits(), []);

  useEffect(() => {
    let cancelled = false;
    api.getDiagramOptions()
      .then((res) => {
        if (cancelled) return;
        setOptions(res);
        setOptionsError('');
      })
      .catch((err: unknown) => { if (!cancelled) setOptionsError(toErrMsg(err)); });
    return () => { cancelled = true; };
  }, [optionsAttempt]);

  // The fluids, and those of the model (known after a solve)
  useEffect(() => {
    let cancelled = false;
    api.getCoolPropFluids()
      .then((res) => {
        if (cancelled) return;
        setFluids(res.fluids);
        setModelFluids(res.modelFluids);
      })
      .catch(() => {});
    return () => { cancelled = true; };
  }, [lastResult]);

  // The inferred variables of a successful solve
  useEffect(() => {
    if (!lastResult?.success) return;
    let cancelled = false;
    api.getInferredVariables()
      .then((res) => { if (!cancelled) setInferred({ result: lastResult, vars: res.variables }); })
      .catch(() => { if (!cancelled) setInferred({ result: lastResult, vars: [] }); });
    return () => { cancelled = true; };
  }, [lastResult]);

  const solved = !!lastResult?.success;
  const currentVars = solved && inferred.result === lastResult ? inferred.vars : null;
  const indices = useMemo(() => arrayIndices(lastResult?.variables), [lastResult]);
  const arrayKeys = useMemo(() => {
    const keys = new Map<string, string>();
    for (const v of currentVars ?? []) {
      const match = v.name.match(/^(.+)\[\d+\]$/);
      const key = propertyKey(v.inferredProperty);
      if (match && key && !keys.has(match[1])) keys.set(match[1], key);
    }
    return keys;
  }, [currentVars]);

  // The fluids offered: pure fluids, and the model's fluids of the other types the server draws (its
  // CoolProp mixtures, which are not listed: they are made from their string); the server refuses the
  // others (ideal gases, humid air, ...) for now
  const drawableFluids = useMemo(() => {
    const inModel = (name: string) => modelFluids.some((m) => sameFluid(m, name));
    const types = options?.fluid_types ?? ['Real'];
    const listed = fluids.filter((f) => f.type === 'Real' || (types.includes(f.type) && inModel(f.name)));
    const unlisted: FluidInfo[] = modelFluids
      .filter((m) => !fluids.some((f) => sameFluid(f.name, m)))
      .map((name) => ({ name, coolpropName: name, type: 'Mixture', hasDome: false }));
    return [...unlisted, ...listed];
  }, [fluids, modelFluids, options]);
  const defaultFluid = useMemo(() => {
    const model = drawableFluids.find((f) => modelFluids.some((m) => sameFluid(m, f.name)));
    return (model ?? drawableFluids.find((f) => f.name === FALLBACK_FLUID) ?? drawableFluids[0])?.name ?? '';
  }, [drawableFluids, modelFluids]);

  const stored = diagrams.diagrams[0] ?? null;
  // The definition shown: the stored one, or a default chart (stored at its first change)
  const def = useMemo<DiagramDefinition | null>(() => {
    if (stored) return stored;
    if (!options || !defaultFluid) return null;
    const states = currentVars ? inferredStates(currentVars, defaultFluid) : [];
    return defaultDefinition(defaultFluid, options, [newInferredCycle(states)]);
  }, [stored, options, defaultFluid, currentVars]);

  // After a solve (or a change of fluid), the stored overlays follow the model's variables
  useEffect(() => {
    if (!stored || !solved || currentVars === null) return;
    const cycles = stored.cycles.map((c) => {
      if (isInferredCycle(c)) return { ...c, states: inferredStates(currentVars, stored.fluid) };
      if (isArraysCycle(c)) return arraysCycle(overlayOfCycle(c), indices, c);
      return c;
    });
    if (!sameJson(cycles, stored.cycles)) setDiagramDefinition({ ...stored, cycles });
  }, [stored, solved, currentVars, indices, setDiagramDefinition]);

  const input = useMemo(() => (def ? sceneInput(def, theme) : null), [def, theme]);

  const onViewChange = useCallback(
    (view: DiagramViewSI | null) => { if (def) setDiagramDefinition({ ...def, spec: withView(def.spec, view) }); },
    [def, setDiagramDefinition],
  );

  const liveStates = useMemo(
    () => (currentVars && def ? inferredStates(currentVars, def.fluid) : null),
    [currentVars, def],
  );

  return (
    <div className="diagram-tab">
      {optionsError && (
        <div className="diagram-error">
          Diagram options not available: {optionsError}{' '}
          <button className="diagram-retry" onClick={() => setOptionsAttempt((n) => n + 1)}>Retry</button>
        </div>
      )}
      {def && options && (
        <>
          <DiagramControls
            def={def}
            options={options}
            fluids={drawableFluids}
            modelFluids={modelFluids}
            onChange={setDiagramDefinition}
          />
          <DiagramOverlays
            def={def}
            onChange={setDiagramDefinition}
            solved={solved}
            inferredStates={liveStates}
            arrays={indices}
            arrayKeys={arrayKeys}
            draft={arraysDraft}
            setDraft={setArraysDraft}
          />
        </>
      )}
      <DiagramView input={input} revision={lastResult} onViewChange={onViewChange} />
    </div>
  );
}
