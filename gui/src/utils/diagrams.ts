/**
 * diagrams.ts — the diagram definitions of the model (`diagrams.json`, docs/coolplot_integration.md §3.4).
 *
 * The document lives in the model store (`diagrams`); the server session keeps a copy, which goes into
 * the ZIP bundle. `showSessionModel()` reads it back after every operation that replaces the model; an
 * edit made in the Diagrams tab is written to the session with `PUT /api/v1/diagram/definitions`, debounced
 * (`watchDiagramEdits`), and sent at once before a Save or an operation that replaces the model
 * (`flushDiagramsPut`).
 */

import { api } from '../api/client';
import type { DiagramDefinition, DiagramsDocument } from '../api/types';
import { useModelStore } from '../stores/modelStore';
import { toErrMsg } from './errors';

/** A document without diagrams: what a model has until the Diagrams tab is changed. */
export function emptyDiagramsDocument(): DiagramsDocument {
  return { format: 'coolsolve-diagrams', version: 1, diagrams: [] };
}

const isObject = (v: unknown): v is Record<string, unknown> => typeof v === 'object' && v !== null && !Array.isArray(v);

/**
 * The document as read from the session, or an empty one when it is not a diagrams document. Diagrams
 * without a fluid or a spec are left out; every other field is kept as it is.
 */
export function normalizeDiagramsDocument(raw: unknown): DiagramsDocument {
  if (!isObject(raw) || raw.format !== 'coolsolve-diagrams' || !Array.isArray(raw.diagrams)) {
    return emptyDiagramsDocument();
  }
  const diagrams = raw.diagrams
    .filter((d): d is Record<string, unknown> => isObject(d) && typeof d.fluid === 'string' && isObject(d.spec))
    .map((d, i) => ({
      ...d,
      name: typeof d.name === 'string' && d.name ? d.name : `Diagram ${i + 1}`,
      cycles: Array.isArray(d.cycles) ? d.cycles : [],
    }) as DiagramDefinition);
  return { format: 'coolsolve-diagrams', version: 1, diagrams };
}

// ---------------------------------------------------------------------------
// Writing the document to the session
// ---------------------------------------------------------------------------

const PUT_DEBOUNCE_MS = 400;
let timer: ReturnType<typeof setTimeout> | undefined;
let lastError = '';

async function putNow(): Promise<void> {
  timer = undefined;
  const { diagrams, addConsoleLine } = useModelStore.getState();
  try {
    await api.putDiagramDefinitions(diagrams);
    lastError = '';
  } catch (err: unknown) {
    const msg = toErrMsg(err);
    if (msg !== lastError) addConsoleLine(`>>> Diagram definitions not saved in the session: ${msg}`);
    lastError = msg;
  }
}

/**
 * Write the document to the session, debounced, whenever the GUI edits it (`diagramsEdits` grows); a
 * load (`diagramsEdits` back to 0) cancels a pending write. Returns the unsubscribe function.
 */
export function watchDiagramEdits(): () => void {
  return useModelStore.subscribe((state, prev) => {
    if (state.diagramsEdits === prev.diagramsEdits) return;
    if (timer !== undefined) clearTimeout(timer);
    timer = state.diagramsEdits > 0 ? setTimeout(() => { void putNow(); }, PUT_DEBOUNCE_MS) : undefined;
  });
}

/** Send a pending write now: before a Save, and before the session's model is replaced (New, Open, Back). */
export async function flushDiagramsPut(): Promise<void> {
  if (timer === undefined) return;
  clearTimeout(timer);
  await putNow();
}
