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
  isArray: boolean;
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
