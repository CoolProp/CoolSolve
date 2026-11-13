// API types matching the CoolSolve REST API responses

export interface HealthResponse {
  status: string;
  coolpropReady: boolean;
}

export interface ParseError {
  line: number;
  column: number;
  message: string;
  context: string;
}

export interface VariableInfo {
  name: string;
  units: string;
  unitSource: 'code' | 'inferred' | '';  // Source of unit annotation
  isArray: boolean;
  isImposed: boolean;        // True if variable = constant in the code
  imposedValue?: number;     // The imposed constant value (if isImposed)
}

/** One `$INCLUDE` directive of the model, as resolved by the server (POST /api/v1/parse). */
export interface IncludeInfo {
  line: number;
  name: string;                     // as written: "hx_ees_function"
  file: string;                     // the file read: "hx_ees_function.eescode"
  origin: 'library' | 'project';    // a library model, or a file of the model folder
  libraryId: string;                // "CSL-0331" for a library model, else ''
  description: string;              // "library model hx_ees_function (CSL-0331)" / "my_props.eescode"
  routines: string[];               // "FUNCTION HX", "PROCEDURE p"
  nested: string[];                 // the files it includes in turn
}

export interface ParseResponse {
  success: boolean;
  equationCount: number;
  totalLines: number;
  variableCount?: number;
  isSquare?: boolean;
  errors: ParseError[];
  variables?: VariableInfo[];
  diagnostics?: Diagnostic[];
  includes?: IncludeInfo[];
}

export interface BlockResult {
  id: number;
  success: boolean;
  status: string;
  iterations: number;
  maxResidual: number;
  errorMessage: string;
}

export interface Diagnostic {
  severity: 'Error' | 'Warning' | 'Info';
  code: string;
  message: string;
  source: string;
  line?: number;
  column?: number;
}

export interface TimingInfo {
  coolprop_warmup_ms: number;
  parse_ms: number;
  ir_ms: number;
  infer_ms: number;
  analysis_ms: number;
  solve_ms: number;
  total_ms: number;
}

export interface SolveResponse {
  success: boolean;
  status: string;
  errorMessage: string;
  totalIterations: number;
  blocksEvaluated: number;
  totalTimeMs: number;
  timing: TimingInfo;
  variables: Record<string, number>;
  stringVariables: Record<string, string>;
  blockResults: BlockResult[];
  detailedError?: string;
  equationCount?: number;
  variableCount?: number;
  isSquare?: boolean;
  totalBlocks?: number;
  largestBlock?: number;
  latexReportAvailable?: boolean;
  diagnostics?: Diagnostic[];
  // Integral (INTEGRAL/$IntegralTable) trajectory — present only when the
  // model is a dynamic/DAE solve (see IntegralTableData).
  integralTable?: IntegralTableData;
  integralCsvName?: string;
}

export interface SolvedVariable {
  name: string;
  value: number | string;
  isArray?: boolean;
  isString?: boolean;
}

export interface VariablesResponse {
  variables: SolvedVariable[];
}

export interface FileOpenResponse {
  success: boolean;
  filePath: string;
  modelName: string;
  hasInitials: boolean;
  hasSol: boolean;
  hasConf: boolean;
}

/** A model of the library menu (verified models only). */
export interface LibraryMenuModel {
  id: string;     // permanent ID, e.g. "CSL-0012"
  name: string;
  title: string;
  level: number;  // 1 introductory ... 4 research
  kind: string;   // steady | dynamic | optimisation | functions
  functions: number;  // FUNCTIONs and PROCEDUREs: what `$INCLUDE <name>` imports
}

/** A node of the library taxonomy holding models (a category or a sub-category). */
export interface LibraryMenuGroup {
  key: string;
  title: string;
  models: LibraryMenuModel[];
}

export interface LibraryMenuCategory extends LibraryMenuGroup {
  count: number;  // models of the category, sub-categories included
  children: LibraryMenuGroup[];
}

/** Response of GET /api/v1/library/menu */
export interface LibraryMenuResponse {
  available: boolean;
  reason?: string;  // why the library is not available
  total: number;
  categories: LibraryMenuCategory[];
  snapshot: { commit?: string; models?: number; content_hash?: string };
}

/** Response of POST /api/v1/library/open */
export interface LibraryOpenResponse {
  success: boolean;
  libraryId: string;
  title: string;
  status: string;
  file: string;
  variants: string[];
  modelName: string;
  hasInitials: boolean;
  hasSol: boolean;
  hasConf: boolean;
  tables: number;
}

export interface UploadResponse {
  success: boolean;
  files: string[];
  modelName: string;
  mainFile?: string;  // the model, when the ZIP holds several .eescode files (the others can be included)
}

/** Response of POST /api/v1/includes/open: a new session holding the included file. */
export interface IncludeOpenResponse {
  sessionId: string;
  modelName: string;
  description: string;
}

// Lookup table types
export interface LookupTableInfo {
  name: string;
  columns: string[];
  rows: number;
}

export interface LookupTablesResponse {
  tables: LookupTableInfo[];
}

// SSE progress events from /api/v1/solve/stream
export interface SSEBlockEvent {
  type: 'block';
  block: number;
  totalBlocks: number;
  event: 'start' | 'done' | 'fail';
  iterations: number;
  residualNorm: number;
  message: string;
}

export interface SSEProgressEvent {
  type: 'progress' | 'start';
  message: string;
  /** Present on 'start' events when this run is a Deep Search ("Try Harder"). */
  deepSearch?: boolean;
}

export interface SSEDoneEvent {
  type: 'done' | 'error';
  message: string;
  result?: SolveResponse;
}

export type SSEEvent = SSEBlockEvent | SSEProgressEvent | SSEDoneEvent;

// Debug output
export interface DebugFile {
  name: string;
  size: number;
}

export interface DebugFilesResponse {
  files: DebugFile[];
}

// Inferred variables (thermodynamic property inference)
export interface InferredVariable {
  name: string;
  value: number;
  inferredProperty: string; // T, P, H, S, D, Q, V (viscosity), VOL (specific volume), L, C
  inferredFluid: string;
  units: string;
  isArray: boolean;
}

export interface InferredVariablesResponse {
  variables: InferredVariable[];
}

// CoolProp fluids
export interface FluidInfo {
  name: string;
  coolpropName: string;
  type: string;
  hasDome: boolean;
  /** An incompressible solution: the mass fractions CoolProp has it for (the concentration of its charts) */
  fractionRange?: [number, number];
}

export interface CoolPropFluidsResponse {
  fluids: FluidInfo[];
  modelFluids: string[];
}

// ============================================================================
// Property diagrams (CoolPlot): docs/coolplot_integration.md §3.1, §3.4, §4.1
// ============================================================================

/**
 * The diagram types of coolplot-core, vertical axis first ("ph" is log(p)-h): those of a fluid, and the charts of
 * humid air, the psychrometric chart "WT" (humidity ratio over dry-bulb temperature) and the Mollier h,x chart "ThxW"
 * (its temperature scale over the humidity ratio).
 */
export type DiagramName = 'ph' | 'Ts' | 'hs' | 'ps' | 'prho' | 'Trho' | 'pT' | 'Th' | 'WT' | 'ThxW';

/** A family of lines a diagram type can draw: its quantity key and default number of lines. */
export interface DiagramFamilyOption {
  key: string;
  label: string;
  symbol: string;
  default_num: number;
}

/** A diagram type offered by the server: its name, short title and axis quantity keys. */
export interface DiagramTypeOption {
  name: DiagramName;
  /** The short name shown: the axis symbols, vertical first ("p-h", "T-s", "p-ρ"); "psychrometric", "Mollier h-x" */
  title: string;
  kind: string;
  x: string;
  y: string;
  x_log: boolean;
  y_log: boolean;
  families: DiagramFamilyOption[];
}

/**
 * A custom unit system of CoolPlot: quantity key -> [label, scale, offset, difference label?],
 * a display value being SI * scale + offset.
 */
export interface UnitSystemRecord {
  name: string;
  units: Record<string, [string, number, number] | [string, number, number, string]>;
}

/** The `units` of a CoolPlot spec: a preset name ("SI", "KSI", "engineering", "IP") or a unit record. */
export type DiagramUnits = string | UnitSystemRecord;

/** The unit system of the model, in which every chart is drawn: `record` goes into the spec's `units`. */
export interface DiagramModelUnits {
  record: UnitSystemRecord;
  /** Its unit table: quantity key -> [label, scale, offset, difference label?] */
  units: Record<string, unknown>;
}

/** A region preset (T-p box of the lines). */
export interface DiagramRegionOption {
  name: string;
  limits: Record<string, unknown>;
}

/** What a kind of fluid ("Real", "Mixture", "Incompressible", "IdealGas", "HumidAir") offers (`fluid_kinds` of the options). */
export interface DiagramFluidKind {
  /** The diagram types offered */
  diagram_types: DiagramName[];
  /** The keys of the families its lines can have */
  families: string[];
  /** A saturation curve: dome, lines of constant quality, critical point, saturation scale, region presets */
  saturation: boolean;
  /** A concentration (mass fraction) is part of the fluid: an incompressible solution */
  fractions: boolean;
  /** The T-p box of a fluid without saturation curve where the spec gives no absolute limit (null: the fluid's own) */
  default_box: (number | null)[] | null;
  /** The region of a new chart: a preset (ACHP for mixtures, DEF for pure fluids); null without saturation curve;
   *  for humid air a region record {"bounds": {"T": [lo, hi], "W": [lo, hi]}} in SI */
  default_region?: string | Record<string, unknown> | null;
  /** CoolPlot's kind of the charts: "fluid", or "humid_air" (psychrometric and Mollier charts) */
  chart_kind?: 'fluid' | 'humid_air';
  /** Humid air: the total pressure of a new chart, Pa */
  total_pressure?: number;
  /** Humid air: the total pressures HAPropsSI covers, Pa */
  total_pressure_range?: [number, number];
  /** Humid air: the temperatures HAPropsSI covers, K */
  temperature_range?: [number, number];
}

/** Response of GET /api/v1/diagram/options */
export interface DiagramOptionsResponse {
  diagram_types: DiagramTypeOption[];
  model_units: DiagramModelUnits;
  regions: DiagramRegionOption[];
  themes: string[];
  label_modes: string[];
  marked_points: string[];
  /** Fluid types (GET /api/v1/coolprop/fluids) the server draws */
  fluid_types: string[];
  /** What each fluid type offers */
  fluid_kinds?: Record<string, DiagramFluidKind>;
  coolplot: { core_version: string; commit: string };
}

/** A family of lines of a CoolPlot spec (IsolineSet), as the GUI writes it. */
export interface IsolineSpec {
  key: string;
  labels: boolean;
  num: number;
  rounding: boolean;
  [field: string]: unknown;  // other IsolineSet fields of a spec written elsewhere are kept
}

/** Zoom window of a CoolPlot spec, in SI: [[x_min, x_max], [y_min, y_max]], null for the whole axis. */
export type DiagramViewSI = [[number, number] | null, [number, number] | null];

/** Settings of a CoolPlot spec, as far as the GUI writes them. */
export interface DiagramSettingsSpec {
  /** A theme name ("default", "dark", "print", "coolpack"); absent: the GUI's theme ("default" or "dark") */
  look?: string;
  labels?: { mode?: string; symbols?: 'unicode' | 'ascii'; [field: string]: unknown };
  branding?: { references?: 'none' | 'short' | 'full'; versions?: boolean; [field: string]: unknown };
  [field: string]: unknown;
}

/**
 * A CoolPlot chart spec record (DiagramSpec.to_dict) without "fluid" and "cycles", which the
 * server writes. Fields the GUI does not edit are kept as they are.
 */
export interface CoolPlotSpec {
  diagram: DiagramName;
  units: DiagramUnits;
  /** A region preset ("DEF", "ACHP", "ORC"), or the bounds of a spec written elsewhere */
  limits: string | Record<string, unknown>;
  isolines: IsolineSpec[];
  marked_points?: string[];
  saturation_scale?: { Q: number; step_SI?: number | null };
  view_SI?: DiagramViewSI;
  settings: DiagramSettingsSpec;
  [field: string]: unknown;
}

/** One state of a cycle definition: CoolPlot property keys ("p", "h", "T", ...) -> model variable names. */
export interface CycleStateDefinition {
  number: number;
  label?: string;
  vars: Record<string, string>;
}

/**
 * A cycle of a diagram definition (§4.1): states as variable names, resolved by the server with
 * the session's last solve. CoolPlot's CycleSpec fields ("numbers", styles) pass through.
 */
export interface CycleDefinition {
  name: string;
  numbers?: boolean;
  /** Explicit states... */
  states?: CycleStateDefinition[];
  /** ...or arrays: property key -> array name, state i from X[i] */
  arrays?: Record<string, string>;
  /** Chains of state numbers; [] draws the states without lines; absent: one loop through the states
   *  in ascending order, closed unless `closed` is false */
  loops?: number[][];
  closed?: boolean;
  /** One kind for every connection, or CoolPlot's record [[[a, b], kind], ...] (connections not
   *  listed: "auto"); absent: "auto" */
  paths?: CyclePathKind | CyclePathEntry[];
  /** Segments per drawn path */
  steps?: number;
  [field: string]: unknown;
}

/** How the path between two states is drawn ("auto": CoolPlot chooses from the two states). */
export type CyclePathKind = 'auto' | 'straight' | 'isobaric' | 'isenthalpic' | 'isentropic' | 'isothermal' | 'polytropic';

/** The path of one connection of a cycle: [[from, to], kind]. */
export type CyclePathEntry = [[number, number], CyclePathKind];

/** One diagram of a model (`diagrams.json`, §3.4). */
export interface DiagramDefinition {
  name: string;
  /** A CoolSolve fluid name (FluidRegistry) */
  fluid: string;
  /** The concentration of an incompressible solution: one mass fraction (CoolPlot's fluid record) */
  fractions?: number[];
  /** Humid air: the total pressure of the chart in Pa (CoolPlot's HumidAirSpec); absent: 101325 Pa */
  p_Pa?: number;
  spec: CoolPlotSpec;
  cycles: CycleDefinition[];
  [field: string]: unknown;
}

/** The diagrams of a model: `diagrams.json` in the ZIP bundle, GET/PUT /api/v1/diagram/definitions. */
export interface DiagramsDocument {
  format: 'coolsolve-diagrams';
  version: 1;
  diagrams: DiagramDefinition[];
}

/** What a diagram shows: the body of POST /api/v1/diagram/scene without the figure size and labels. */
export interface DiagramSceneInput {
  fluid: string;
  /** The mass fraction of an incompressible solution */
  fractions?: number[];
  /** Humid air: the total pressure of the chart, Pa */
  p_Pa?: number;
  spec: CoolPlotSpec;
  cycles: CycleDefinition[];
}

/** Body of POST /api/v1/diagram/scene */
export interface DiagramSceneRequest extends DiagramSceneInput {
  /** The whole figure, in pt */
  figure_size_pt: [number, number];
  /** The `labels` of the previous response (same fluid and diagram), so that labels keep their places */
  previous_labels?: Record<string, unknown>[];
}

/** Body of POST /api/v1/diagram/export: a chart request and the format of its file */
export interface DiagramExportRequest extends DiagramSceneRequest {
  format: 'svg' | 'pdf';
}

/** A display unit of an axis: display = SI * scale + offset. */
export interface DiagramAxisUnit {
  label: string;
  scale: number;
  offset: number;
}

/** The axes of a scene; limits in display units. */
export interface DiagramAxes {
  x_key: string;
  y_key: string;
  x_label: string;
  y_label: string;
  x_limits: [number, number];
  y_limits: [number, number];
  x_log: boolean;
  y_log: boolean;
  x_unit: DiagramAxisUnit;
  y_unit: DiagramAxisUnit;
  /** Size of the plot area, pt */
  plot_size_pt?: [number, number];
}

/** A state of a cycle as drawn: its variables and its full state in SI. */
export interface DiagramSceneState {
  cycle: number;
  number: number;
  label?: string;
  vars: Record<string, string>;
  values_SI: Record<string, number>;
}

/** A state of the solved model as the server groups it (the property calls its variables share). */
export interface DiagramModelState {
  number: number;
  label: string | null;
  /** CoolPlot property key -> variable name */
  vars: Record<string, string>;
}

/** Response of GET /api/v1/diagram/states: the states of the last solve per fluid (empty without a solve). */
export interface DiagramStatesResponse {
  fluids: Record<string, DiagramModelState[]>;
  /** The concentration (mass fraction) of the model's incompressible solutions, from its property calls */
  fractions?: Record<string, number[]>;
}

/** Body of POST /api/v1/diagram/state: the point of a chart, in SI. */
export interface DiagramPointRequest {
  fluid: string;
  fractions?: number[];
  p_Pa?: number;
  diagram: string;
  x_SI: number;
  y_SI: number;
}

/**
 * The full state at a point of a chart, SI (Q = -1 in a single phase); a property may be null. Humid air: T, W, RH,
 * Twb, Tdp, h, v (per kg of dry air), Thx and p.
 */
export type DiagramPointState = Partial<Record<
  'T' | 'p' | 'h' | 's' | 'rho' | 'u' | 'v' | 'Q' | 'W' | 'RH' | 'Twb' | 'Tdp' | 'Thx', number | null>>;

/** Response of POST /api/v1/diagram/state (null: no state there). */
export interface DiagramPointResponse {
  state: DiagramPointState | null;
}

/** Response of POST /api/v1/diagram/scene */
export interface DiagramSceneResponse {
  /** coolplot-core's SVG document */
  svg: string;
  axes: DiagramAxes;
  /** Isoline labels, sent back as `previous_labels` */
  labels: Record<string, unknown>[];
  states: DiagramSceneState[];
  /** The spec as composed (fluid and cycles filled in) */
  spec: Record<string, unknown>;
  diagnostics: { curves: Record<string, unknown>[]; warnings: string[] };
  timing_ms: { curves: number; compose: number; svg: number };
}

// ============================================================================
// Parametric Study Types
// ============================================================================

/** A sweep variable definition for parametric study */
export interface SweepVariable {
  name: string;
  values: number[];
}

/** One grid point result from a parametric study */
export interface ParametricGridPoint {
  index: number;
  success: boolean;
  overrides: Record<string, number>;
  variables?: Record<string, number>;
  errorMessage?: string;
}

/** Full parametric study response */
export interface ParametricStudyResult {
  success: boolean;
  totalPoints: number;
  successCount: number;
  failCount: number;
  sweepVariables: SweepVariable[];
  results: ParametricGridPoint[];
}

/** Saved parametric study (for persistence in bundle and UI state) */
export interface SavedParametricStudy {
  id: string;
  name: string;
  sweepVariables: SweepVariable[];
  result: ParametricStudyResult;
  timestamp: number;
}

// ============================================================================
// Integral Table Types (equation-based dynamic / DAE solving)
// ============================================================================

/**
 * Columnar trajectory table produced by an `INTEGRAL` / `$IntegralTable`
 * solve. Mirrors the JSON emitted by the backend
 * (`integralResultToJSON` in src/server.cpp):
 *
 *   {
 *     "integrationVar": "t",
 *     "columns": ["t", "y", "dydt"],
 *     "data": { "t": [...], "y": [...], "dydt": [...] },
 *     "numRows": 9,
 *     "csvName": "model-integral.csv",
 *     "totalSteps": 40,
 *     "rejectedSteps": 0
 *   }
 *
 * The `data` field is keyed by column name; every column has the same
 * length (`numRows`). The integration variable is always `columns[0]`.
 */
export interface IntegralTableData {
  integrationVar: string;
  columns: string[];
  data: Record<string, number[]>;
  numRows: number;
  csvName: string;
  totalSteps?: number;
  rejectedSteps?: number;
}

/** User-defined unit override */
export interface UserUnitOverride {
  variableName: string;
  units: string;
}

// ============================================================================
// LaTeX Report Types
// ============================================================================

/** Response from GET /api/v1/latex/report */
export interface LatexReportResponse {
  available: boolean;
  content: string;
}

/** A plot image to embed in the LaTeX report */
export interface LatexPlotImage {
  /** Filename to use in the .tex (e.g. "thermo_diagrams.png") */
  name: string;
  /** Base64-encoded PNG data (data URI prefix is stripped automatically) */
  data: string;
}

/** Request body for POST /api/v1/latex/compile */
export interface LatexCompileRequest {
  compiler?: string;
  plots?: LatexPlotImage[];
}

// ============================================================================
// Usage Log Statistics (hidden dashboard at /stats)
// ============================================================================

/** Per-day attempt counters (UTC) */
export interface UsageDailyCount {
  date: string;
  attempts: number;
  successes: number;
}

/** A "name → count" entry in a top-ranking (models, IPs) */
export interface UsageNameCount {
  name: string;
  count: number;
}

/** Aggregated solve-attempt statistics from GET /api/v1/stats/log */
export interface UsageStatsResponse {
  valid: boolean;
  totalAttempts: number;
  successes: number;
  failures: number;
  parseErrors: number;
  tryHarderAttempts: number;
  uniqueIps: number;
  malformedLines: number;
  meanMs: number;
  medianMs: number;
  p95Ms: number;
  daily: UsageDailyCount[];
  outcomes: Record<string, number>;
  durationHistogram: { edges: number[]; counts: number[] };
  topModels: UsageNameCount[];
  topIps: UsageNameCount[];
}