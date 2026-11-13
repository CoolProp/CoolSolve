/**
 * definition.ts — the diagram definition the Diagrams tab edits (docs/coolplot_integration.md §3.4, §4.1).
 *
 * A definition is a fluid, a CoolPlot spec record (without "fluid" and "cycles") and cycle definitions
 * that name the model's variables. The GUI writes the fields below and keeps every other field of a
 * spec written elsewhere; the chart is always in the model's units (`model_units` of the options). Its
 * cycles are of two kinds: the inferred states (the states of the solved model as the server groups them,
 * `GET /api/v1/diagram/states`, drawn as markers without lines) and the arrays overlays (state i from the
 * elements X[i] and Y[i] of the arrays of the chart's two axis quantities, joined in index order).
 */

import type {
  CoolPlotSpec,
  CycleDefinition,
  CycleStateDefinition,
  DiagramDefinition,
  DiagramFamilyOption,
  DiagramFluidKind,
  DiagramName,
  DiagramOptionsResponse,
  DiagramSceneInput,
  DiagramTypeOption,
  DiagramModelState,
  DiagramStatesResponse,
  DiagramViewSI,
  FluidInfo,
  InferredVariable,
  IsolineSpec,
} from '../../api/types';

/** Name of the cycle of the inferred states (the former "state points" overlay). */
export const INFERRED_CYCLE = 'Inferred states';

export const DEFAULT_DIAGRAM: DiagramName = 'ph';
export const DEFAULT_REGION = 'DEF';
export const DEFAULT_LABEL_MODE = 'rails';
/** The fluid drawn when the model has none */
export const FALLBACK_FLUID = 'Water';

/** Property codes of CoolSolve's variable inference (src/variable_inference.cpp) -> CoolPlot keys. */
const PROPERTY_KEYS: Record<string, string> = { T: 'T', P: 'p', H: 'h', S: 's', D: 'rho', U: 'u', VOL: 'v', Q: 'Q' };

/** The CoolPlot key of an inferred property code, or undefined (viscosity, conductivity, ...). */
const propertyKey = (code: string): string | undefined => PROPERTY_KEYS[code];

export const sameFluid = (a: string, b: string) => a.toLowerCase() === b.toLowerCase();

/** The option of a diagram type, or undefined. */
export function diagramOption(options: DiagramOptionsResponse, name: string): DiagramTypeOption | undefined {
  return options.diagram_types.find((d) => d.name === name);
}

/** Whether coolplot-core draws a saturation scale in this chart: p or T vertical, neither horizontal. */
export function hasSaturationScale(type: { x: string; y: string }): boolean {
  return (type.y === 'p' || type.y === 'T') && type.x !== 'p' && type.x !== 'T';
}

/**
 * A family of lines as the GUI adds it. A fluid without saturation curve gets round values: the families spaced
 * evenly over the region of an incompressible liquid begin and end at its corners, where a line is a single point.
 */
const familyLine = (key: string, num: number, round = false): IsolineSpec => ({ key, labels: true, num, rounding: round });

/** The lines of a new chart: every family but the internal energy (CoolPlot's `set_default_lines`). */
function defaultIsolines(type: DiagramTypeOption): IsolineSpec[] {
  return type.families.filter((f) => f.key !== 'u').map((f) => familyLine(f.key, f.default_num));
}

// ---------------------------------------------------------------------------
// Kinds of fluids (phase 3): what a pure fluid, a mixture, an incompressible solution or an ideal gas offers
// ---------------------------------------------------------------------------

/** What a type of fluid (`FluidInfo.type`) offers; undefined when the server does not say (everything then). */
export function fluidKind(options: DiagramOptionsResponse, type: string | undefined): DiagramFluidKind | undefined {
  return type ? options.fluid_kinds?.[type] : undefined;
}

/** Whether a kind of fluid is humid air, drawn in the psychrometric and Mollier charts at a total pressure. */
export const isHumidAirKind = (kind: DiagramFluidKind | undefined) => kind?.chart_kind === 'humid_air';

/** PsychrometricChart's region: -10 to 50 °C, 0 to 30 g/kg (SI). */
const DEFAULT_AIR_REGION: Record<string, [number, number]> = { T: [263.15, 323.15], W: [0, 0.03] };

/** The bounds of a region record ({"bounds": {key: [lo, hi]}}), or undefined for a preset name. */
export function regionBounds(limits: CoolPlotSpec['limits'] | null | undefined): Record<string, [number, number]> | undefined {
  if (!limits || typeof limits === 'string') return undefined;
  const bounds = (limits as { bounds?: unknown }).bounds;
  return bounds && typeof bounds === 'object' ? bounds as Record<string, [number, number]> : undefined;
}

/**
 * The region of a chart of humid air: bounds on its temperature (the dry-bulb temperature T, or the temperature scale
 * Thx of the Mollier chart) and on the humidity ratio W, taken from `limits` (either temperature), else from `fallback`
 * (the kind's default region), else PsychrometricChart's.
 */
export function humidAirRegion(
  limits: CoolPlotSpec['limits'] | null | undefined,
  type: { x: string; y: string },
  fallback?: CoolPlotSpec['limits'] | null,
): Record<string, unknown> {
  const temperature = type.y === 'Thx' ? 'Thx' : 'T';
  const given = regionBounds(limits) ?? regionBounds(fallback) ?? DEFAULT_AIR_REGION;
  const T = given.T ?? given.Thx ?? DEFAULT_AIR_REGION.T;
  const W = given.W ?? DEFAULT_AIR_REGION.W;
  return { bounds: { [temperature]: [...T], W: [...W] } };
}

/** The diagram types a kind of fluid offers. */
export function diagramTypesOf(options: DiagramOptionsResponse, kind: DiagramFluidKind | undefined): DiagramTypeOption[] {
  return kind ? options.diagram_types.filter((d) => kind.diagram_types.includes(d.name)) : options.diagram_types;
}

/** The families of a diagram type a kind of fluid can draw. */
export function familiesOf(type: DiagramTypeOption | undefined, kind: DiagramFluidKind | undefined): DiagramFamilyOption[] {
  const all = type?.families ?? [];
  return kind ? all.filter((f) => kind.families.includes(f.key)) : all;
}

/**
 * The spec for a fluid of another kind: a diagram type the kind offers, its families (those set kept, the others
 * on but the internal energy), the kind's default region preset (`DEF`, `ACHP` for mixtures; a fluid without
 * saturation curve is drawn over its default box), and without saturation curve no critical point nor saturation
 * scale.
 */
export function fitToKind(spec: CoolPlotSpec, kind: DiagramFluidKind | undefined, options: DiagramOptionsResponse): CoolPlotSpec {
  if (!kind) return spec;
  const type = diagramOption(options, kind.diagram_types.includes(spec.diagram) ? spec.diagram : kind.diagram_types[0]);
  if (!type) return spec;
  const before = new Map(spec.isolines.map((s) => [s.key, s]));
  const isolines = familiesOf(type, kind)
    .filter((f) => before.has(f.key) || f.key !== 'u')
    .map((f) => before.get(f.key) ?? familyLine(f.key, f.default_num, !kind.saturation));
  const { view_SI: _view, saturation_scale, marked_points, ...rest } = spec;
  const limits = isHumidAirKind(kind) ? humidAirRegion(kind.default_region, type) : kind.default_region || DEFAULT_REGION;
  const next: CoolPlotSpec = { ...rest, diagram: type.name, isolines, limits };
  if (kind.saturation && saturation_scale && hasSaturationScale(type)) next.saturation_scale = saturation_scale;
  if (kind.saturation) next.marked_points = marked_points ?? ['critical'];
  return next;
}

/**
 * The concentration of a solution's charts: that of the model's calls, else the middle of the fluid's range;
 * undefined for a fluid that takes none.
 */
export function defaultFractions(fluid: FluidInfo | undefined, modelFractions?: Record<string, number[]>): number[] | undefined {
  if (!fluid?.fractionRange) return undefined;
  const fromModel = Object.entries(modelFractions ?? {}).find(([name]) => sameFluid(name, fluid.name))?.[1];
  if (fromModel?.length === 1) return fromModel;
  const [lo, hi] = fluid.fractionRange;
  return [Math.round(500 * (lo + hi)) / 1000];
}

/**
 * The light lines of the Diagrams tab: every family of lines thinner and partly transparent, so that the saturation
 * curve (1.5 pt) and the cycles (2 pt) of CoolPlot's themes stand out. CoolPlot's themes draw isolines at 0.75 pt
 * and lines of constant quality at 0.5 pt, fully opaque.
 */
export const LIGHT_LINE = { width_pt: 0.5, alpha: 0.55 };

/**
 * CoolPlot's `style_overrides` of the light lines: LIGHT_LINE for every family a chart of the options can draw
 * (the saturation curve keeps the theme's style). A CoolPlot spec field, so that an exported spec draws the same.
 */
export function lightLineOverrides(options: DiagramOptionsResponse): Record<string, typeof LIGHT_LINE> {
  const keys = new Set<string>();
  for (const type of options.diagram_types) for (const f of type.families) keys.add(f.key);
  return Object.fromEntries([...keys].sort().map((k) => [k, { ...LIGHT_LINE }]));
}

/** The line weight of a spec: "light" (the overrides of lightLineOverrides), "theme" (none), or "custom". */
export function lineWeight(spec: CoolPlotSpec): 'light' | 'theme' | 'custom' {
  const overrides = spec.style_overrides as Record<string, Record<string, unknown>> | undefined;
  if (!overrides || Object.keys(overrides).length === 0) return 'theme';
  // Field by field: the server keeps JSON objects with their keys sorted
  const light = Object.values(overrides).every(
    (o) => Object.keys(o).length === 2 && o.width_pt === LIGHT_LINE.width_pt && o.alpha === LIGHT_LINE.alpha,
  );
  return light ? 'light' : 'custom';
}

/** The spec with the light lines ("light") or the theme's ("theme"). */
export function withLineWeight(spec: CoolPlotSpec, weight: 'light' | 'theme', options: DiagramOptionsResponse): CoolPlotSpec {
  const { style_overrides: _overrides, ...rest } = spec;
  return weight === 'light' ? { ...rest, style_overrides: lightLineOverrides(options) } : rest;
}

/** A new spec for a diagram type. */
export function defaultSpec(type: DiagramTypeOption, options: DiagramOptionsResponse): CoolPlotSpec {
  const regions = options.regions.map((r) => r.name);
  return {
    diagram: type.name,
    units: options.model_units.record,
    limits: regions.includes(DEFAULT_REGION) || regions.length === 0 ? DEFAULT_REGION : regions[0],
    isolines: defaultIsolines(type),
    marked_points: ['critical'],
    style_overrides: lightLineOverrides(options),
    settings: {
      // No "look": the GUI's theme (sceneInput). No "branding": CoolPlot's badge, the references of the
      // equation of state and the versions (definitions saved before keep the branding they have)
      labels: { mode: DEFAULT_LABEL_MODE, symbols: 'unicode' },
    },
  };
}

/** A new definition: the default chart of `fluid` (of the kind given, at the concentration given), with the cycles given. */
export function defaultDefinition(
  fluid: string,
  options: DiagramOptionsResponse,
  cycles: CycleDefinition[],
  kind?: DiagramFluidKind,
  fractions?: number[],
): DiagramDefinition | null {
  const type = diagramOption(options, DEFAULT_DIAGRAM) ?? options.diagram_types[0];
  if (!type) return null;
  const spec = kind ? fitToKind(defaultSpec(type, options), kind, options) : defaultSpec(type, options);
  return { name: 'Diagram 1', fluid, ...(fractions ? { fractions } : {}), spec, cycles };
}

/**
 * The spec for another diagram type: the families of both keep their settings, the new ones are on
 * (but the internal energy, and those the kind of fluid cannot draw), the zoom is dropped and the
 * saturation scale kept where the chart has one.
 */
export function withDiagramType(
  spec: CoolPlotSpec,
  from: DiagramTypeOption | undefined,
  to: DiagramTypeOption,
  kind?: DiagramFluidKind,
): CoolPlotSpec {
  const before = new Map(spec.isolines.map((s) => [s.key, s]));
  const oldFamilies = new Set(from?.families.map((f) => f.key) ?? []);
  const isolines: IsolineSpec[] = [];
  for (const f of familiesOf(to, kind)) {
    const kept = before.get(f.key);
    if (kept) isolines.push(kept);
    else if (!oldFamilies.has(f.key) && f.key !== 'u') isolines.push(familyLine(f.key, f.default_num, kind ? !kind.saturation : false));
  }
  const { view_SI: _view, saturation_scale, ...rest } = spec;
  const next: CoolPlotSpec = { ...rest, diagram: to.name, isolines };
  if (saturation_scale && hasSaturationScale(to)) next.saturation_scale = saturation_scale;
  // The psychrometric chart is bounded in T, the Mollier chart in its temperature scale Thx: the same values
  if (to.kind === 'humid_air') next.limits = humidAirRegion(spec.limits, to, kind?.default_region);
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
  round = false,
): CoolPlotSpec {
  const existing = spec.isolines.find((s) => s.key === key);
  if (change.on === false) return { ...spec, isolines: spec.isolines.filter((s) => s.key !== key) };
  if (!existing) {
    if (change.on !== true) return spec;
    return { ...spec, isolines: [...spec.isolines, familyLine(key, defaultNum, round)] };
  }
  const updated: IsolineSpec = { ...existing };
  if (change.num !== undefined) updated.num = change.num;
  if (change.rounding !== undefined) updated.rounding = change.rounding;
  return { ...spec, isolines: spec.isolines.map((s) => (s.key === key ? updated : s)) };
}

/**
 * The request of a definition: the theme follows the GUI when the spec names none; cycles with nothing to draw (no
 * states, an overlay without its two arrays) are left out.
 */
export function sceneInput(def: DiagramDefinition, guiTheme: 'light' | 'dark'): DiagramSceneInput {
  const settings = def.spec.settings ?? {};
  const look = settings.look || (guiTheme === 'dark' ? 'dark' : 'default');
  return {
    fluid: def.fluid,
    ...(def.fractions ? { fractions: def.fractions } : {}),
    ...(def.p_Pa !== undefined ? { p_Pa: def.p_Pa } : {}),
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

/**
 * The total pressure (Pa) of the model's humid air: that of the first state of `fluid` in the server's grouping whose
 * pressure is a solved variable (CoolSolve's pressures are in Pa); undefined when none is.
 */
export function modelPressure(
  response: DiagramStatesResponse | null,
  fluid: string,
  variables: Record<string, number> | undefined,
): number | undefined {
  const entry = Object.entries(response?.fluids ?? {}).find(([name]) => sameFluid(name, fluid));
  for (const state of entry?.[1] ?? []) {
    const value = state.vars.p !== undefined ? variables?.[state.vars.p] : undefined;
    if (value !== undefined && Number.isFinite(value) && value > 0) return value;
  }
  return undefined;
}

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
// The arrays overlays: state i from the elements i of two arrays, those of the chart's axis quantities
// ---------------------------------------------------------------------------

/** The arrays of a solve: their number of elements, and the quantity (CoolPlot key) of those the model ties to one. */
export interface SolvedArrays {
  /** Array name -> number of elements, by name */
  sizes: Map<string, number>;
  /** Array name -> CoolPlot key of its elements */
  quantities: Map<string, string>;
}

/** The array of an element variable (`h[2]` -> `h`), or undefined. */
const arrayOf = (name: string) => name.match(/^(.+)\[\d+\]$/)?.[1];

/** The axis quantities of the fluid charts, which an array can be named after (`P[i]`, `h[i]`, `T[i]`, `s[i]`). */
const NAMED_QUANTITIES = ['T', 'p', 'h', 's', 'rho'];

/**
 * The arrays of a solve. The quantity of an array is that of the first of its elements found in the server's states
 * (`GET /api/v1/diagram/states`: CoolPlot keys, those of humid air included), else in the variable inference
 * (`GET /api/v1/variables/inferred`), else the quantity it is named after: the arrays a model fills for its
 * diagrams (`P[1] = P_1`, as the library models do) are in no property call.
 */
export function solvedArrays(
  variables: Record<string, number> | undefined,
  states: DiagramStatesResponse | null,
  inferred: InferredVariable[],
): SolvedArrays {
  const sizes = new Map<string, number>();
  for (const name of Object.keys(variables ?? {}).sort((a, b) => a.localeCompare(b))) {
    const array = arrayOf(name);
    if (array) sizes.set(array, (sizes.get(array) ?? 0) + 1);
  }
  const quantities = new Map<string, string>();
  const add = (name: string, key: string | undefined) => {
    const array = arrayOf(name);
    if (array && key && !quantities.has(array)) quantities.set(array, key);
  };
  for (const list of Object.values(states?.fluids ?? {}))
    for (const state of list) for (const [key, name] of Object.entries(state.vars)) add(name, key);
  for (const v of inferred) add(v.name, propertyKey(v.inferredProperty));
  for (const array of sizes.keys())
    if (!quantities.has(array)) {
      const named = NAMED_QUANTITIES.find((key) => key.toLowerCase() === array.toLowerCase());
      if (named) quantities.set(array, named);
    }
  return { sizes: new Map([...sizes].sort(([a], [b]) => a.localeCompare(b))), quantities };
}

/**
 * The quantities a state of an overlay gives on a chart, [horizontal, vertical]: the axis keys, but for the temperature
 * scale of the Mollier chart (Thx), which no model computes: its arrays are dry-bulb temperatures (T).
 */
export function axisKeys(type: DiagramTypeOption): [string, string] {
  const stateKey = (key: string) => (key === 'Thx' ? 'T' : key);
  return [stateKey(type.x), stateKey(type.y)];
}

/** An arrays overlay as its row shows it: the arrays of the horizontal and vertical quantities ('' until chosen). */
export interface ArraysOverlay {
  x: string;
  y: string;
  /** Join the last state to the first */
  closed: boolean;
}

export const isArraysCycle = (c: CycleDefinition) => !!c.arrays;

/** The overlay of an arrays cycle on a chart whose states give the quantities `keys` (axisKeys). */
export function overlayOf(cycle: CycleDefinition, keys: [string, string]): ArraysOverlay {
  return { x: cycle.arrays?.[keys[0]] ?? '', y: cycle.arrays?.[keys[1]] ?? '', closed: cycle.closed !== false };
}

/**
 * The cycle of an overlay: its arrays by quantity, one chain through the common elements in index order (closed on
 * request), paths along the processes. The CoolPlot fields of `previous` (styles) are kept; normalizeDefinition names it.
 */
export function arraysCycle(o: ArraysOverlay, keys: [string, string], previous?: CycleDefinition): CycleDefinition {
  const { loops: _loops, ...kept } = previous ?? { name: '', numbers: true };
  const arrays: Record<string, string> = {};
  if (o.x) arrays[keys[0]] = o.x;
  if (o.y) arrays[keys[1]] = o.y;
  return { ...kept, arrays, closed: o.closed, paths: 'auto' };
}

/** A new overlay, its arrays to be chosen (normalizeDefinition gives it the default ones). */
export const newArraysCycle = (): CycleDefinition => ({ name: '', numbers: true, arrays: {}, closed: false, paths: 'auto' });

/**
 * The default array of a quantity: the first array of that quantity not in `taken` (the arrays of the overlays
 * before), else the first one, else ''.
 */
function defaultArray(arrays: SolvedArrays, key: string, taken: string[]): string {
  const candidates = [...arrays.sizes.keys()].filter((name) => arrays.quantities.get(name) === key);
  return candidates.find((name) => !taken.includes(name)) ?? candidates[0] ?? '';
}

/** `base`, else `base 2`, `base 3`, ...: the first name not in `used`. */
function freeName(base: string, used: Set<string>): string {
  let name = base;
  for (let n = 2; used.has(name); n++) name = `${base} ${n}`;
  return name;
}

/**
 * The definition as the Diagrams tab shows and writes it: in the model's units; with the cycles the tab shows, the
 * inferred states and the arrays overlays (others, the cycles of the former state editor, are dropped); each overlay on
 * the axes of the chart, an axis quantity it has no array for taking the default one (an array of that quantity, so
 * that a change of diagram type maps the overlays to the new axes), and named after its arrays ("T, s"), uniquely.
 */
export function normalizeDefinition(def: DiagramDefinition, options: DiagramOptionsResponse, arrays: SolvedArrays): DiagramDefinition {
  const type = diagramOption(options, def.spec.diagram);
  const keys = type ? axisKeys(type) : null;
  const names = new Set([INFERRED_CYCLE]);
  const taken: [string[], string[]] = [[], []];
  const cycles: CycleDefinition[] = [];
  for (const cycle of def.cycles) {
    if (isInferredCycle(cycle)) cycles.push(cycle);
    if (!isArraysCycle(cycle)) continue;
    if (!keys) {  // a diagram type the server does not offer: the overlay is left as it is
      cycles.push(cycle);
      continue;
    }
    const o = overlayOf(cycle, keys);
    o.x = o.x || defaultArray(arrays, keys[0], taken[0]);
    o.y = o.y || defaultArray(arrays, keys[1], taken[1]);
    taken[0].push(o.x);
    taken[1].push(o.y);
    const name = freeName(`${o.y || '?'}, ${o.x || '?'}`, names);
    names.add(name);
    cycles.push({ ...arraysCycle(o, keys, cycle), name });
  }
  return { ...def, spec: { ...def.spec, units: options.model_units.record }, cycles };
}

/** JSON equality of two plain values (definitions are JSON). */
export const sameJson = (a: unknown, b: unknown) => JSON.stringify(a) === JSON.stringify(b);
