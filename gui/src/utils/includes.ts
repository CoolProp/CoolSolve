/**
 * includes.ts — the files a model includes (`$INCLUDE`) and the include line of a library model.
 */

import { api } from '../api/client';
import { adoptTabSession, getTabSessionId } from '../api/session';
import { useModelStore } from '../stores/modelStore';
import { toErrMsg } from './errors';
import { showSessionModel } from './sessionModel';

/** The directive that imports the functions and procedures of a library model. */
export const includeDirective = (name: string) => `$INCLUDE ${name}`;

// The window that opens an included file is a plain link (`target="_blank"`), which no pop-up blocker
// stops; the click leaves this note in localStorage (shared by the windows of the site) to tell the new
// window which session includes the file.
const OPENER_KEY = 'coolsolve_include_opener';
const OPENER_MAX_AGE_MS = 60_000;

/** Address of a new window showing the file of `$INCLUDE name`. */
export const includeWindowUrl = (name: string) => `/?include=${encodeURIComponent(name)}`;

/** On the click of an include link: note this window's session for the window it opens. */
export function noteIncludeOpener(name: string): void {
  try {
    localStorage.setItem(OPENER_KEY, JSON.stringify({ name, from: getTabSessionId(), time: Date.now() }));
  } catch {
    // Storage blocked: the new window falls back to the shared cookie session
  }
}

/** The session that asked to open `name` (null: the cookie session), consumed once. */
function takeIncludeOpener(name: string): string | null {
  try {
    const note = JSON.parse(localStorage.getItem(OPENER_KEY) ?? 'null');
    if (note && note.name === name && Date.now() - note.time < OPENER_MAX_AGE_MS) {
      localStorage.removeItem(OPENER_KEY);
      return typeof note.from === 'string' ? note.from : null;
    }
  } catch {
    // No note: the cookie session
  }
  return null;
}

/**
 * In a window opened on `/?include=name`: ask the server, as the including session, for a new session
 * holding the included file, make it this window's session and show its model. Errors are reported in
 * the console.
 */
export async function openIncludeInThisWindow(name: string): Promise<void> {
  const { addConsoleLine } = useModelStore.getState();
  try {
    const res = await api.openInclude(name, takeIncludeOpener(name));
    if (!adoptTabSession(res.sessionId)) throw new Error('this browser keeps no per-window session');
    await showSessionModel();
    addConsoleLine(`>>> Opened ${res.description} (included by the model of the other window)`);
  } catch (err: unknown) {
    addConsoleLine(`>>> ERROR opening $INCLUDE ${name}: ${toErrMsg(err)}`);
  }
}

/** Copy `text` to the clipboard; false when the browser refuses. */
export async function copyToClipboard(text: string): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch {
    // No clipboard API (a page served over plain http from another host): the selection fallback
    const area = document.createElement('textarea');
    area.value = text;
    area.style.position = 'fixed';
    area.style.opacity = '0';
    document.body.appendChild(area);
    area.select();
    const ok = document.execCommand('copy');
    document.body.removeChild(area);
    return ok;
  }
}
