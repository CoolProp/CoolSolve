/**
 * definition.ts — the diagram definition the Diagrams tab edits (docs/coolplot_integration.md §3.4, §4.1).
 *
 * A definition is a fluid, a CoolPlot spec record (without "fluid" and "cycles") and cycle definitions
 * that name the model's variables. The GUI writes the fields below and keeps every other field of a
 * spec written elsewhere. Its cycles are of three kinds: the inferred states (the states of the solved
 * model as the server groups them, `GET /api/v1/diagram/states`, drawn as markers without lines), the
 * arrays (state i from X[i], Y[i], joined in index order) and the user's cycles of the state editor
 * (states, loops and paths chosen by the user).
 */

import type {
  CoolPlotSpec,
  CycleDefinition,
  CycleStateDefinition,
  DiagramDefinition,
  DiagramName,
  DiagramOptionsResponse,
  DiagramSceneInput,
  DiagramTypeOption,
  DiagramUnits,
  CyclePathEntry,
  CyclePathKind,
  DiagramModelState,
  DiagramStatesResponse,
  DiagramViewSI,
  IsolineSpec,
} from '../../api/types';

/** Name of the cycle of the inferred states (the former "state points" overlay). */
export const INFERRED_CYCLE = 'Inferred states';
/** Name of the cycle of the arrays (the former "array path" overlay). */
export const ARRAYS_CYCLE = 'Arrays';

export const DEFAULT_DIAGRAM: DiagramName = 'ph';
/** Decision D7: printed charts' units (deg C, bar, kJ/kg) */
export const DEFAULT_UNITS = 'engineering';
export const DEFAULT_REGION = 'DEF';
export const DEFAULT_LABEL_MODE = 'rails';
/** The fluid drawn when the model has none */
export const FALLBACK_FLUID = 'Water';

/** The fluid quantities a cycle state can name (CoolPlot keys). */
export const STATE_KEYS = ['T', 'p', 'h', 's', 'rho', 'u', 'v', 'Q'];

/** Property codes of CoolSolve's variable inference (src/variable_inference.cpp) -> CoolPlot keys. */
const PROPERTY_KEYS: Record<string, string> = { T: 'T', P: 'p', H: 'h', S: 's', D: 'rho', U: 'u', VOL: 'v', Q: 'Q' };

/** The CoolPlot key of an inferred property code, or undefined (viscosity, specific volume, ...). */
export const propertyKey = (code: string): string | undefined => PROPERTY_KEYS[code];

/** The name of a `units` value: a preset name, or the name of a unit record. */
export const unitsName = (units: DiagramUnits): string => (typeof units === 'string' ? units : units.name);

export const sameFluid = (a: string, b: string) => a.toLowerCase() === b.toLowerCase();

/** The option of a diagram type, or undefined. */
export function diagramOption(options: DiagramOptionsResponse, name: string): DiagramTypeOption | undefined {
  return options.diagram_types.find((d) => d.name === name);
}

/** The `units` value of a unit system offered by the server (a preset name, or the "coolsolve" record). */
export function unitsRecord(options: DiagramOptionsResponse, name: string): DiagramUnits | undefined {
  return options.unit_systems.find((u) => u.name === name)?.record;
}

/** Whether coolplot-core draws a saturation scale in this chart: p or T vertical, neither horizontal. */
export function hasSaturationScale(type: { x: string; y: string }): boolean {
  return (type.y === 'p' || type.y === 'T') && type.x !== 'p' && type.x !== 'T';
}

const familyLine = (key: string, num: number): IsolineSpec => ({ key, labels: true, num, rounding: false });

/** The lines of a new chart: every family but the internal energy (CoolPlot's `set_default_lines`). */
function defaultIsolines(type: DiagramTypeOption): IsolineSpec[] {
  return type.families.filter((f) => f.key !== 'u').map((f) => familyLine(f.key, f.default_num));
}

/** A new spec for a diagram type. */
export function defaultSpec(type: DiagramTypeOption, options: DiagramOptionsResponse): CoolPlotSpec {
  const regions = options.regions.map((r) => r.name);
  return {
    diagram: type.name,
    units: unitsRecord(options, options.default_unit_system || DEFAULT_UNITS)
      ?? unitsRecord(options, DEFAULT_UNITS) ?? options.unit_systems[0]?.record ?? DEFAULT_UNITS,
    limits: regions.includes(DEFAULT_REGION) || regions.length === 0 ? DEFAULT_REGION : regions[0],
    isolines: defaultIsolines(type),
    marked_points: ['critical'],
    settings: {
      // No "look": the GUI's theme (sceneInput). No "branding": CoolPlot's badge, the references of the
      // equation of state and the versions (definitions saved before keep the branding they have)
      labels: { mode: DEFAULT_LABEL_MODE, symbols: 'unicode' },
    },
  };
}

/** A new definition: the default chart of `fluid`, with the cycles given. */
export function defaultDefinition(
  fluid: string,
  options: DiagramOptionsResponse,
  cycles: CycleDefinition[],
): DiagramDefinition | null {
  const type = diagramOption(options, DEFAULT_DIAGRAM) ?? options.diagram_types[0];
  if (!type) return null;
  return { name: 'Diagram 1', fluid, spec: defaultSpec(type, options), cycles };
}

/**
 * The spec for another diagram type: the families of both keep their settings, the new ones are on
 * (but the internal energy), the zoom is dropped and the saturation scale kept where the chart has one.
 */
export function withDiagramType(spec: CoolPlotSpec, from: DiagramTypeOption | undefined, to: DiagramTypeOption): CoolPlotSpec {
  const before = new Map(spec.isolines.map((s) => [s.key, s]));
  const oldFamilies = new Set(from?.families.map((f) => f.key) ?? []);
  const isolines: IsolineSpec[] = [];
  for (const f of to.families) {
    const kept = before.get(f.key);
    if (kept) isolines.push(kept);
    else if (!oldFamilies.has(f.key) && f.key !== 'u') isolines.push(familyLine(f.key, f.default_num));
  }
  const { view_SI: _view, saturation_scale, ...rest } = spec;
  const next: CoolPlotSpec = { ...rest, diagram: to.name, isolines };
  if (saturation_scale && hasSaturationScale(to)) next.saturation_scale = saturation_scale;
  return next;
}

/** The spec with a zoom window, or without one (null). */
export function withView(spec: CoolPlotSpec, view: DiagramViewSI | null): CoolPlotSpec {
  const { view_SI: _view, ...rest } = spec;
  return view ? { ...rest, view_SI: view } : rest;
}

/** A family of lines on or off; `num` and `rounding` change one that is on. */
export function withFamily(
  spec: CoolPlotSpec,
  key: string,
  change: { on?: boolean; num?: number; rounding?: boolean },
  defaultNum: number,
): CoolPlotSpec {
  const existing = spec.isolines.find((s) => s.key === key);
  if (change.on === false) return { ...spec, isolines: spec.isolines.filter((s) => s.key !== key) };
  if (!existing) {
    if (change.on !== true) return spec;
    return { ...spec, isolines: [...spec.isolines, familyLine(key, defaultNum)] };
  }
  const updated: IsolineSpec = { ...existing };
  if (change.num !== undefined) updated.num = change.num;
  if (change.rounding !== undefined) updated.rounding = change.rounding;
  return { ...spec, isolines: spec.isolines.map((s) => (s.key === key ? updated : s)) };
}

/** The request of a definition: the theme follows the GUI when the spec names none; cycles without states are left out. */
export function sceneInput(def: DiagramDefinition, guiTheme: 'light' | 'dark'): DiagramSceneInput {
  const settings = def.spec.settings ?? {};
  const look = settings.look || (guiTheme === 'dark' ? 'dark' : 'default');
  return {
    fluid: def.fluid,
    spec: { ...def.spec, settings: { ...settings, look } },
    cycles: def.cycles.filter(isDrawable),
  };
}

function isDrawable(cycle: CycleDefinition): boolean {
  if (cycle.states && cycle.states.length > 0) return true;
  return !!cycle.arrays && Object.values(cycle.arrays).filter(Boolean).length >= 2;
}

// ---------------------------------------------------------------------------
// The inferred states: the server's grouping of the solved variables
// ---------------------------------------------------------------------------

/** The states of `fluid` in the server's grouping (`GET /api/v1/diagram/states`), as cycle states. */
export function modelStates(response: DiagramStatesResponse | null, fluid: string): CycleStateDefinition[] {
  const entry = Object.entries(response?.fluids ?? {}).find(([name]) => sameFluid(name, fluid));
  return (entry?.[1] ?? []).map((s: DiagramModelState) => ({
    number: s.number,
    ...(s.label ? { label: s.label } : {}),
    vars: { ...s.vars },
  }));
}

/** A new cycle of the inferred states: markers and numbers, no lines (`loops: []`). */
export function newInferredCycle(states: CycleStateDefinition[]): CycleDefinition {
  return { name: INFERRED_CYCLE, numbers: true, states, loops: [] };
}

export const isInferredCycle = (c: CycleDefinition) => c.name === INFERRED_CYCLE && !c.arrays;

// ---------------------------------------------------------------------------
// The arrays (former "array path" overlay)
// ---------------------------------------------------------------------------

/** The arrays of a solve (`X[i]` variables): name -> indices in ascending order. */
export function arrayIndices(variables: Record<string, number> | undefined): Map<string, number[]> {
  const out = new Map<string, number[]>();
  for (const name of Object.keys(variables ?? {})) {
    const match = name.match(/^(.+)\[(\d+)\]$/);
    if (!match) continue;
    const list = out.get(match[1]) ?? [];
    list.push(parseInt(match[2], 10));
    out.set(match[1], list);
  }
  for (const list of out.values()) list.sort((a, b) => a - b);
  return out;
}

/** The arrays overlay as the controls show it. */
export interface ArraysOverlay {
  on: boolean;
  /** CoolPlot property keys of the two arrays */
  keys: [string, string];
  /** Array names ('' until chosen) */
  arrays: [string, string];
  /** Join the last state to the first */
  closed: boolean;
  /** Paths between the states: along the processes ("auto") or straight */
  paths: 'auto' | 'straight';
}

export const isArraysCycle = (c: CycleDefinition) => c.name === ARRAYS_CYCLE && !!c.arrays;

/** The overlay of a stored arrays cycle. */
export function overlayOfCycle(cycle: CycleDefinition): ArraysOverlay {
  const entries = Object.entries(cycle.arrays ?? {});
  const chain = cycle.loops?.[0];
  return {
    on: true,
    keys: [entries[0]?.[0] ?? 'h', entries[1]?.[0] ?? 'p'],
    arrays: [entries[0]?.[1] ?? '', entries[1]?.[1] ?? ''],
    closed: chain ? chain.length > 2 && chain[0] === chain[chain.length - 1] : cycle.loops === undefined && cycle.closed !== false,
    paths: cycle.paths === 'straight' ? 'straight' : 'auto',
  };
}

/** Whether the overlay names two arrays of two different quantities. */
export const overlayComplete = (o: ArraysOverlay) =>
  o.on && !!o.arrays[0] && !!o.arrays[1] && o.keys[0] !== o.keys[1];

/**
 * The chain of the arrays' common indices, closed on request; undefined when the indices are not known
 * (no solve): `previous` (the stored chain) is then kept.
 */
export function arraysLoops(
  o: ArraysOverlay,
  indices: Map<string, number[]>,
  previous?: number[][],
): number[][] | undefined {
  const a = indices.get(o.arrays[0]);
  const b = indices.get(o.arrays[1]);
  let common: number[] | undefined;
  if (a && b) common = a.filter((i) => b.includes(i));
  else if (previous?.[0]) {
    const chain = previous[0];
    common = chain.length > 2 && chain[0] === chain[chain.length - 1] ? chain.slice(0, -1) : chain;
  }
  if (!common) return undefined;
  if (common.length < 2) return [];
  return [o.closed && common.length > 2 ? [...common, common[0]] : common];
}

/** The arrays cycle of a complete overlay; `previous` keeps the CoolPlot fields of the stored cycle. */
export function arraysCycle(
  o: ArraysOverlay,
  indices: Map<string, number[]>,
  previous?: CycleDefinition,
): CycleDefinition {
  const cycle: CycleDefinition = {
    ...(previous ?? { name: ARRAYS_CYCLE, numbers: true }),
    arrays: { [o.keys[0]]: o.arrays[0], [o.keys[1]]: o.arrays[1] },
    paths: o.paths,
  };
  // The indices known (a solve): the chain itself; else the server joins the states in index order
  const loops = arraysLoops(o, indices, previous?.loops);
  delete cycle.loops;
  delete cycle.closed;
  if (loops) cycle.loops = loops;
  else cycle.closed = o.closed;
  return cycle;
}

// ---------------------------------------------------------------------------
// The user's cycles (state editor)
// ---------------------------------------------------------------------------

/** The kinds of path between two states; "auto" lets CoolPlot choose from the two states. */
export const PATH_KINDS: CyclePathKind[] = ['auto', 'straight', 'isobaric', 'isenthalpic', 'isentropic', 'isothermal', 'polytropic'];

/** A cycle of the state editor: neither the inferred states nor the arrays overlay. */
export const isUserCycle = (c: CycleDefinition) => !isInferredCycle(c) && !isArraysCycle(c);

/** A name not used by another cycle: `base`, else `base 2`, `base 3`, ... */
export function freeCycleName(cycles: CycleDefinition[], base: string): string {
  const used = new Set(cycles.map((c) => c.name));
  if (!used.has(base)) return base;
  let n = 2;
  while (used.has(`${base} ${n}`)) n++;
  return `${base} ${n}`;
}

/** A user cycle through `states`, one closed loop in the order given, paths along the processes. */
export function newUserCycle(name: string, states: CycleStateDefinition[]): CycleDefinition {
  const numbers = states.map((s) => s.number);
  const loops = numbers.length >= 2 ? [[...numbers, ...(numbers.length > 2 ? [numbers[0]] : [])]] : [];
  return { name, numbers: true, states, loops, paths: 'auto' };
}

/** Loops as the editor writes them: "1-2-3-4-1; 2-5". */
export const formatLoops = (loops: number[][] | undefined) => (loops ?? []).map((chain) => chain.join('-')).join('; ');

/** Loops from their text; an error names what is wrong (a chain of one state, an unknown state). */
export function parseLoops(text: string, numbers: number[]): { loops: number[][] } | { error: string } {
  const loops: number[][] = [];
  for (const part of text.split(/[;\n]/).map((p) => p.trim()).filter(Boolean)) {
    const chain = part.split(/[-,\s>]+/).filter(Boolean).map((n) => Number(n));
    if (chain.some((n) => !Number.isInteger(n))) return { error: `"${part}" is not a chain of state numbers` };
    if (chain.length < 2) return { error: `"${part}": a chain joins at least two states` };
    const unknown = chain.filter((n) => !numbers.includes(n));
    if (unknown.length > 0) return { error: `no state ${unknown[0]}` };
    loops.push(chain);
  }
  return { loops };
}

/** The connections the loops draw, each once, in order: "a-b". */
export function connectionsOf(loops: number[][] | undefined): [number, number][] {
  const seen = new Set<string>();
  const out: [number, number][] = [];
  for (const chain of loops ?? []) {
    for (let i = 0; i + 1 < chain.length; i++) {
      const key = `${chain[i]}-${chain[i + 1]}`;
      if (!seen.has(key)) {
        seen.add(key);
        out.push([chain[i], chain[i + 1]]);
      }
    }
  }
  return out;
}

/** The paths of a cycle as the editor shows them: one kind for the cycle, and the connections that differ. */
export function pathsOf(cycle: CycleDefinition): { kind: CyclePathKind; overrides: Record<string, CyclePathKind> } {
  const paths = cycle.paths;
  if (!Array.isArray(paths)) return { kind: paths ?? 'auto', overrides: {} };
  const listed = new Map(paths.map(([[a, b], kind]) => [`${a}-${b}`, kind] as const));
  const connections = connectionsOf(cycle.loops);
  // The kind of the cycle: "auto" when a connection is not listed, else the most frequent one
  let kind: CyclePathKind = 'auto';
  if (connections.length > 0 && connections.every(([a, b]) => listed.has(`${a}-${b}`))) {
    const counts = new Map<CyclePathKind, number>();
    for (const k of listed.values()) counts.set(k, (counts.get(k) ?? 0) + 1);
    kind = [...counts].sort((x, y) => y[1] - x[1])[0][0];
  }
  const overrides: Record<string, CyclePathKind> = {};
  for (const [key, k] of listed) if (k !== kind) overrides[key] = k;
  return { kind, overrides };
}

/**
 * The `paths` of a cycle from a kind and the connections that differ: the kind alone without
 * overrides; else CoolPlot's record, listing every connection of the loops unless the kind is "auto"
 * (connections not listed are "auto").
 */
export function pathsValue(
  kind: CyclePathKind,
  overrides: Record<string, CyclePathKind>,
  loops: number[][] | undefined,
): CyclePathKind | CyclePathEntry[] {
  const connections = connectionsOf(loops);
  const differing = connections.filter(([a, b]) => overrides[`${a}-${b}`] && overrides[`${a}-${b}`] !== kind);
  if (differing.length === 0) return kind;
  const listed = kind === 'auto' ? differing : connections;
  return listed.map(([a, b]): CyclePathEntry => [[a, b], overrides[`${a}-${b}`] ?? kind]);
}

/** JSON equality of two plain values (definitions are JSON). */
export const sameJson = (a: unknown, b: unknown) => JSON.stringify(a) === JSON.stringify(b);
