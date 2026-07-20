// REST API client for CoolSolve backend

import type {
  HealthResponse,
  ParseResponse,
  SolveResponse,
  VariablesResponse,
  FileOpenResponse,
  ExamplesResponse,
  UploadResponse,
  SSEEvent,
  InferredVariablesResponse,
  CoolPropFluidsResponse,
  SaturationResponse,
  ParametricStudyResult,
  SweepVariable,
  LatexReportResponse,
  LatexCompileRequest,
  LookupTablesResponse,
  IntegralTableData,
} from './types';

const API_BASE = '/api/v1';

async function request<T>(path: string, options?: RequestInit): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, {
    headers: { 'Content-Type': 'application/json' },
    ...options,
  });
  if (!res.ok) {
    const body = await res.json().catch(() => ({ error: res.statusText }));
    throw new Error(body.error || `HTTP ${res.status}`);
  }
  return res.json();
}

async function requestOptional<T>(path: string): Promise<T | null> {
  const res = await fetch(`${API_BASE}${path}`, {
    headers: { 'Content-Type': 'application/json' },
  });
  if (res.status === 204 || res.status === 404) return null;
  if (!res.ok) return null;
  return res.json();
}

export const api = {
  // Health
  health: () => request<HealthResponse>('/health'),

  // Warmup
  warmup: () => request<{ warmupMs: number }>('/warmup', { method: 'POST' }),

  // Initial file (set by CLI when launched with a file argument)
  getInitialFile: () => requestOptional<{ path: string }>('/initial-file'),

  // Files
  getEescode: () => request<{ content: string; filePath: string }>('/files/eescode'),
  putEescode: (content: string) =>
    request<{ success: boolean }>('/files/eescode', {
      method: 'PUT',
      body: JSON.stringify({ content }),
    }),

  getInitials: () => request<{ content: string }>('/files/initials'),
  putInitials: (content: string) =>
    request<{ success: boolean }>('/files/initials', {
      method: 'PUT',
      body: JSON.stringify({ content }),
    }),

  getSol: () => request<{ content: string }>('/files/sol'),
  getConf: () => request<{ content: string }>('/files/conf'),
  putConf: (content: string) =>
    request<{ success: boolean }>('/files/conf', {
      method: 'PUT',
      body: JSON.stringify({ content }),
    }),

  openFile: (path: string) =>
    request<FileOpenResponse>('/files/open', {
      method: 'POST',
      body: JSON.stringify({ path }),
    }),

  // New model (clear session)
  newModel: () =>
    request<{ success: boolean; hadContent: boolean }>('/new', { method: 'POST' }),

  // Back (restore previous model)
  goBack: () =>
    request<{ success: boolean; modelName: string }>('/back', { method: 'POST' }),

  // Model name
  getModelName: () =>
    request<{ modelName: string }>('/model-name'),
  setModelName: (modelName: string) =>
    request<{ success: boolean }>('/model-name', {
      method: 'PUT',
      body: JSON.stringify({ modelName }),
    }),

  // Parse
  parse: (source: string) =>
    request<ParseResponse>('/parse', {
      method: 'POST',
      body: JSON.stringify({ source }),
    }),

  // Solve (async — returns immediately, results come via SSE)
  // `deepSearch: true` triggers a "Try Harder" run: the deep-search pipeline
  // is used and tearing + symbolic reduction are forced on.
  solve: (options?: { eescode?: string; initials?: string; debug?: boolean; deepSearch?: boolean }) =>
    request<{ status: string }>('/solve', {
      method: 'POST',
      body: JSON.stringify(options || {}),
    }),

  // Subscribe to solve progress via SSE
  subscribeSolveProgress: (onEvent: (event: SSEEvent) => void): EventSource => {
    const es = new EventSource(`${API_BASE}/solve/stream`);
    es.onmessage = (msg) => {
      try {
        const data = JSON.parse(msg.data) as SSEEvent;
        onEvent(data);
      } catch {
        // Ignore unparseable events
      }
    };
    es.onerror = () => {
      es.close();
    };
    return es;
  },

  getSolveResult: () => request<SolveResponse>('/solve/result'),

  cancelSolve: () =>
    request<{ success: boolean; message: string }>('/solve/cancel', { method: 'POST' }),

  // Variables
  getVariables: () => request<VariablesResponse>('/variables'),

  // Update guesses
  updateGuesses: () =>
    request<{ success: boolean }>('/update-guesses', { method: 'POST' }),

  // Examples
  getExamples: () => request<ExamplesResponse>('/examples'),

  // Debug output
  getDebugFiles: () => request<{ files: { name: string; size: number }[] }>('/debug/files'),
  getDebugFile: (name: string) => request<{ name: string; content: string }>(`/debug/file?name=${encodeURIComponent(name)}`),

  // Upload ZIP file
  uploadFiles: async (formData: FormData): Promise<UploadResponse> => {
    const res = await fetch(`${API_BASE}/files/upload`, { method: 'POST', body: formData });
    if (!res.ok) {
      const body = await res.json().catch(() => ({ error: res.statusText }));
      throw new Error(body.error || `HTTP ${res.status}`);
    }
    return res.json();
  },

  // Inferred variables (thermodynamic property inference)
  getInferredVariables: () =>
    request<InferredVariablesResponse>('/variables/inferred'),

  // CoolProp fluids
  getCoolPropFluids: () =>
    request<CoolPropFluidsResponse>('/coolprop/fluids'),

  // Saturation dome
  getSaturationDome: (fluid: string, nPoints?: number) =>
    request<SaturationResponse>('/coolprop/saturation', {
      method: 'POST',
      body: JSON.stringify({ fluid, nPoints: nPoints || 200 }),
    }),

  // Parametric study (async — returns immediately, results come via SSE)
  runParametricStudy: (
    sweepVariables: SweepVariable[],
    options?: { timeout?: number; updateGuesses?: boolean },
  ) =>
    request<{ status: string; totalPoints: number }>('/parametric', {
      method: 'POST',
      body: JSON.stringify({
        sweepVariables,
        timeout: options?.timeout ?? 0,
        updateGuesses: options?.updateGuesses ?? false,
      }),
    }),

  // Get last parametric study result
  getParametricResult: () =>
    request<ParametricStudyResult>('/parametric/result'),

  // Integral (INTEGRAL/$IntegralTable) — columnar JSON of the last trajectory
  getIntegralResult: () =>
    requestOptional<{ integralTable?: IntegralTableData; integralCsvName?: string }>('/integral/result'),

  // Integral trajectory as raw CSV text (survives ZIP bundle round-trip;
  // the GUI tab parses this on bundle load when no live JSON is present)
  getIntegralCSV: async (): Promise<string | null> => {
    const res = await fetch(`${API_BASE}/integral/csv`);
    if (res.status === 204 || res.status === 404) return null;
    if (!res.ok) return null;
    return res.text();
  },

  // LaTeX report
  getLatexReport: () =>
    request<LatexReportResponse>('/latex/report'),

  // Compile LaTeX report to PDF (returns binary PDF blob)
  compileLatexReport: async (body: LatexCompileRequest): Promise<Blob> => {
    const res = await fetch(`${API_BASE}/latex/compile`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });
    if (!res.ok) {
      let errMsg = `HTTP ${res.status}`;
      try {
        const errBody = await res.json();
        if (errBody.error) errMsg = errBody.error;
        if (errBody.log) console.warn('LaTeX log:', errBody.log);
      } catch {
        errMsg += ` ${res.statusText}`;
      }
      throw new Error(errMsg);
    }
    return res.blob();
  },

  // Lookup Tables
  getTables: () => request<LookupTablesResponse>('/tables'),

  getTableCSV: async (name: string): Promise<string> => {
    const res = await fetch(`${API_BASE}/tables/${encodeURIComponent(name)}`);
    if (!res.ok) {
      const body = await res.json().catch(() => ({ error: res.statusText }));
      throw new Error(body.error || `HTTP ${res.status}`);
    }
    return res.text();
  },

  putTable: async (name: string, csvContent: string): Promise<{ success: boolean; name: string; created: boolean }> => {
    const res = await fetch(`${API_BASE}/tables/${encodeURIComponent(name)}`, {
      method: 'PUT',
      headers: { 'Content-Type': 'text/csv' },
      body: csvContent,
    });
    if (!res.ok) {
      const body = await res.json().catch(() => ({ error: res.statusText }));
      throw new Error(body.error || `HTTP ${res.status}`);
    }
    return res.json();
  },

  deleteTable: (name: string) =>
    request<{ success: boolean; name: string }>(
      `/tables/${encodeURIComponent(name)}`,
      { method: 'DELETE' }
    ),
};
