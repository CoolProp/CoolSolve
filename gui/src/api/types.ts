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

export interface ParseResponse {
  success: boolean;
  equationCount: number;
  totalLines: number;
  variableCount?: number;
  isSquare?: boolean;
  errors: ParseError[];
  variables?: VariableInfo[];
}

export interface BlockResult {
  id: number;
  success: boolean;
  status: string;
  iterations: number;
  maxResidual: number;
  errorMessage: string;
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

export interface ExampleFile {
  name: string;
  path: string;
}

export interface ExamplesResponse {
  examples: ExampleFile[];
}

export interface UploadResponse {
  success: boolean;
  files: string[];
  modelName: string;
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
  inferredProperty: string; // T, P, H, S, D, Q, V, L, C
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
}

export interface CoolPropFluidsResponse {
  fluids: FluidInfo[];
  modelFluids: string[];
}

// Saturation dome
export interface SaturationData {
  T: number[];
  P: number[];
  H: number[];
  S: number[];
  D: number[];
}

export interface CriticalPoint {
  T: number;
  P: number;
  H: number;
  S: number;
  D: number;
}

export interface SaturationResponse {
  fluid: string;
  critical: CriticalPoint;
  triplePoint: { T: number };
  liquid: SaturationData;
  vapor: SaturationData;
  nPoints: number;
  computeTime_ms: number;
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