/**
 * DiagramTab — the Diagrams tab: a property diagram of the model drawn by CoolPlot on the server
 * (docs/coolplot_integration.md, phases 1 and 2).
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
  DiagramStatesResponse,
  DiagramViewSI,
  FluidInfo,
  InferredVariable,
  SolveResponse,
} from '../../api/types';
import { useModelStore } from '../../stores/modelStore';
import { useUIStore } from '../../stores/uiStore';
import { watchDiagramEdits } from '../../utils/diagrams';
import { toErrMsg } from '../../utils/errors';
import CycleEditor, { type VariableChoice } from './CycleEditor';
import DiagramControls from './DiagramControls';
import DiagramOverlays from './DiagramOverlays';
import DiagramView from './DiagramView';
import { loadDiagramOptions } from './options';
import {
  FALLBACK_FLUID,
  arrayIndices,
  arraysCycle,
  defaultDefinition,
  defaultFractions,
  fluidKind,
  isArraysCycle,
  isInferredCycle,
  modelStates,
  newInferredCycle,
  overlayOfCycle,
  propertyKey,
  sameFluid,
  sameJson,
  sceneInput,
  withView,
  type ArraysOverlay,
} from './definition';

/** What the server knows of the last solve: its inferred variables and its states, with the solve they belong to. */
interface SolveInfo {
  result: SolveResponse | null;
  vars: InferredVariable[];
  states: DiagramStatesResponse | null;
}

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
  const [info, setInfo] = useState<SolveInfo>({ result: null, vars: [], states: null });
  const [arraysDraft, setArraysDraft] = useState<ArraysOverlay>(
    { on: false, keys: ['h', 'p'], arrays: ['', ''], closed: false, paths: 'auto' });

  // Edits of the definitions go to the session
  useEffect(() => watchDiagramEdits(), []);

  useEffect(() => {
    let cancelled = false;
    loadDiagramOptions()
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

  // The inferred variables and the states of a successful solve
  useEffect(() => {
    if (!lastResult?.success) return;
    let cancelled = false;
    Promise.all([
      api.getInferredVariables().then((r) => r.variables).catch(() => [] as InferredVariable[]),
      api.getDiagramStates().catch(() => null),
    ]).then(([vars, states]) => {
      if (!cancelled) setInfo({ result: lastResult, vars, states });
    });
    return () => { cancelled = true; };
  }, [lastResult]);

  const solved = !!lastResult?.success;
  const current = solved && info.result === lastResult ? info : null;
  const indices = useMemo(() => arrayIndices(lastResult?.variables), [lastResult]);
  const arrayKeys = useMemo(() => {
    const keys = new Map<string, string>();
    for (const v of current?.vars ?? []) {
      const match = v.name.match(/^(.+)\[\d+\]$/);
      const key = propertyKey(v.inferredProperty);
      if (match && key && !keys.has(match[1])) keys.set(match[1], key);
    }
    return keys;
  }, [current]);

  // The fluids offered: those of the types the server draws (pure fluids, incompressible solutions, ideal
  // gases), and the model's CoolProp mixtures, which are not listed: they are made from their string; the
  // server refuses the others (humid air, NH3H2O) for now
  const drawableFluids = useMemo(() => {
    const types = options?.fluid_types ?? ['Real'];
    const listed = fluids.filter((f) => types.includes(f.type));
    const unlisted: FluidInfo[] = modelFluids
      .filter((m) => !fluids.some((f) => sameFluid(f.name, m)))
      .map((name) => ({ name, coolpropName: name, type: 'Mixture', hasDome: false }));
    return [...unlisted, ...listed];
  }, [fluids, modelFluids, options]);
  const defaultFluid = useMemo(() => {
    const model = drawableFluids.find((f) => modelFluids.some((m) => sameFluid(m, f.name)));
    return model ?? drawableFluids.find((f) => f.name === FALLBACK_FLUID) ?? drawableFluids[0];
  }, [drawableFluids, modelFluids]);
  // The concentration of the model's solutions, from its property calls
  const modelFractions = current?.states?.fractions;

  const stored = diagrams.diagrams[0] ?? null;
  // The definition shown: the stored one, or a default chart (stored at its first change); a solution drawn
  // at the concentration of the model (or the middle of its range) until one is chosen
  const def = useMemo<DiagramDefinition | null>(() => {
    if (stored) {
      if (stored.fractions || !options) return stored;
      const given = defaultFractions(drawableFluids.find((f) => sameFluid(f.name, stored.fluid)), modelFractions);
      return given ? { ...stored, fractions: given } : stored;
    }
    if (!options || !defaultFluid) return null;
    return defaultDefinition(defaultFluid.name, options, [newInferredCycle(modelStates(current?.states ?? null, defaultFluid.name))],
      fluidKind(options, defaultFluid.type), defaultFractions(defaultFluid, modelFractions));
  }, [stored, options, defaultFluid, drawableFluids, modelFractions, current]);

  // After a solve (or a change of fluid), the stored overlays follow the model's variables; the user's
  // cycles name their variables and are left as they are
  useEffect(() => {
    if (!stored || !current) return;
    const cycles = stored.cycles.map((c) => {
      if (isInferredCycle(c)) return { ...c, states: modelStates(current.states, stored.fluid) };
      if (isArraysCycle(c)) return arraysCycle(overlayOfCycle(c), indices, c);
      return c;
    });
    if (!sameJson(cycles, stored.cycles)) setDiagramDefinition({ ...stored, cycles });
  }, [stored, current, indices, setDiagramDefinition]);

  const input = useMemo(() => (def ? sceneInput(def, theme) : null), [def, theme]);

  const onViewChange = useCallback(
    (view: DiagramViewSI | null) => { if (def) setDiagramDefinition({ ...def, spec: withView(def.spec, view) }); },
    [def, setDiagramDefinition],
  );

  const liveStates = useMemo(
    () => (current && def ? modelStates(current.states, def.fluid) : null),
    [current, def],
  );

  // The variables the editor offers: those the inference ties to the chart's fluid first (with their
  // property), then the other solved variables
  const variables = useMemo<VariableChoice[]>(() => {
    if (!lastResult?.variables || !def) return [];
    const inferred = new Map((current?.vars ?? []).map((v) => [v.name, v]));
    const choices = Object.keys(lastResult.variables).map((name): VariableChoice => {
      const v = inferred.get(name);
      const ofFluid = !!v && sameFluid(v.inferredFluid, def.fluid);
      return { name, key: ofFluid ? propertyKey(v.inferredProperty) : undefined, ofFluid };
    });
    return choices.sort((a, b) => Number(b.ofFluid) - Number(a.ofFluid) || a.name.localeCompare(b.name));
  }, [lastResult, current, def]);

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
            modelFractions={modelFractions}
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
          <CycleEditor
            def={def}
            onChange={setDiagramDefinition}
            inferredStates={liveStates}
            variables={variables}
          />
        </>
      )}
      <DiagramView input={input} revision={lastResult} onViewChange={onViewChange} />
    </div>
  );
}
