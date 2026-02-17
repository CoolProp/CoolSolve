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

export const api = {
  // Health
  health: () => request<HealthResponse>('/health'),

  // Warmup
  warmup: () => request<{ warmupMs: number }>('/warmup', { method: 'POST' }),

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
  solve: (options?: { eescode?: string; initials?: string; debug?: boolean }) =>
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
};
