/**
 * DiagramTab — the Diagrams tab: a property diagram of the model drawn by CoolPlot on the server
 * (docs/coolplot_integration.md).
 *
 * It edits the first diagram of the model's diagram definitions (`diagrams` in the model store, saved in
 * the session and the ZIP bundle as `diagrams.json`); every definition it shows or writes is in the
 * model's units, with its overlays on the chart's axes (normalizeDefinition). Until the user changes
 * something, it shows a default chart of the model's first fluid that is not stored. The chart is redrawn
 * by `DiagramView` when the definition, the GUI theme or the solution changes.
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
import DiagramControls from './DiagramControls';
import DiagramOverlays from './DiagramOverlays';
import DiagramView from './DiagramView';
import { loadDiagramOptions } from './options';
import {
  FALLBACK_FLUID,
  defaultDefinition,
  defaultFractions,
  fluidKind,
  isHumidAirKind,
  isInferredCycle,
  modelPressure,
  modelStates,
  newInferredCycle,
  normalizeDefinition,
  sameFluid,
  sameJson,
  sceneInput,
  solvedArrays,
  withView,
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
  // The arrays of the last successful solve, once its states and inferred variables are known
  const arrays = useMemo(
    () => solvedArrays(current ? lastResult?.variables : undefined, current?.states ?? null, current?.vars ?? []),
    [current, lastResult],
  );

  // The fluids of the model: those of the variable inference, and those whose states the server found in its
  // property calls (humid air is AirH2O there; the inference names it after its registered fluid, Air_ha)
  const allModelFluids = useMemo(() => {
    const names = [...modelFluids];
    for (const name of Object.keys(current?.states?.fluids ?? {}))
      if (!names.some((m) => sameFluid(m, name))) names.push(name);
    return names;
  }, [modelFluids, current]);

  // The fluids offered: those of the types the server draws (pure fluids, incompressible solutions, ideal
  // gases, humid air), and the model's CoolProp mixtures, which are not listed: they are made from their string;
  // the server refuses the others (NH3H2O) for now
  const drawableFluids = useMemo(() => {
    const types = options?.fluid_types ?? ['Real'];
    const listed = fluids.filter((f) => types.includes(f.type));
    const unlisted: FluidInfo[] = modelFluids
      .filter((m) => !fluids.some((f) => sameFluid(f.name, m)))
      .map((name) => ({ name, coolpropName: name, type: 'Mixture', hasDome: false }));
    return [...unlisted, ...listed];
  }, [fluids, modelFluids, options]);
  const defaultFluid = useMemo(() => {
    const model = drawableFluids.find((f) => allModelFluids.some((m) => sameFluid(m, f.name)));
    return model ?? drawableFluids.find((f) => f.name === FALLBACK_FLUID) ?? drawableFluids[0];
  }, [drawableFluids, allModelFluids]);
  // The concentration of the model's solutions, from its property calls
  const modelFractions = current?.states?.fractions;
  // The total pressure of the model's humid air (the default of its charts), from its states
  const humidAirFluid = useMemo(
    () => drawableFluids.find((f) => options && isHumidAirKind(fluidKind(options, f.type))),
    [drawableFluids, options],
  );
  const airPressure = useMemo(
    () => (humidAirFluid ? modelPressure(current?.states ?? null, humidAirFluid.name, lastResult?.variables) : undefined),
    [humidAirFluid, current, lastResult],
  );

  // Every definition the tab shows or writes: in the model's units, its overlays on the chart's axes
  const normalize = useCallback(
    (d: DiagramDefinition) => (options ? normalizeDefinition(d, options, arrays) : d),
    [options, arrays],
  );
  const edit = useCallback((d: DiagramDefinition) => setDiagramDefinition(normalize(d)), [normalize, setDiagramDefinition]);

  const stored = diagrams.diagrams[0] ?? null;
  // The definition shown: the stored one, or a default chart (stored at its first change); a solution drawn
  // at the concentration of the model (or the middle of its range) until one is chosen, humid air at the
  // total pressure of the model until one is chosen
  const def = useMemo<DiagramDefinition | null>(() => {
    const atModelPressure = (d: DiagramDefinition, type: string | undefined) =>
      d.p_Pa === undefined && airPressure !== undefined && options && isHumidAirKind(fluidKind(options, type))
        ? { ...d, p_Pa: airPressure } : d;
    if (stored) {
      if (!options) return stored;
      const info = drawableFluids.find((f) => sameFluid(f.name, stored.fluid));
      const given = stored.fractions ? undefined : defaultFractions(info, modelFractions);
      return normalize(atModelPressure(given ? { ...stored, fractions: given } : stored, info?.type));
    }
    if (!options || !defaultFluid) return null;
    const created = defaultDefinition(defaultFluid.name, options, [newInferredCycle(modelStates(current?.states ?? null, defaultFluid.name))],
      fluidKind(options, defaultFluid.type), defaultFractions(defaultFluid, modelFractions));
    return created && normalize(atModelPressure(created, defaultFluid.type));
  }, [stored, options, defaultFluid, drawableFluids, modelFractions, current, airPressure, normalize]);

  // After a solve (or a change of fluid), the stored inferred states follow the model's states; the arrays
  // overlays name their arrays and are drawn with the new solve as they are
  useEffect(() => {
    if (!stored || !current) return;
    const cycles = stored.cycles.map((c) => (isInferredCycle(c) ? { ...c, states: modelStates(current.states, stored.fluid) } : c));
    if (!sameJson(cycles, stored.cycles)) edit({ ...stored, cycles });
  }, [stored, current, edit]);

  const input = useMemo(() => (def ? sceneInput(def, theme) : null), [def, theme]);

  const onViewChange = useCallback(
    (view: DiagramViewSI | null) => { if (def) edit({ ...def, spec: withView(def.spec, view) }); },
    [def, edit],
  );

  const liveStates = useMemo(
    () => (current && def ? modelStates(current.states, def.fluid) : null),
    [current, def],
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
            modelFluids={allModelFluids}
            modelFractions={modelFractions}
            modelPressure={airPressure}
            onChange={edit}
          />
          <DiagramOverlays
            def={def}
            options={options}
            onChange={edit}
            solved={solved}
            inferredStates={liveStates}
            arrays={arrays}
          />
        </>
      )}
      <DiagramView input={input} revision={lastResult} onViewChange={onViewChange} />
    </div>
  );
}
