/**
 * Per-tab session id.
 *
 * By default the server identifies a browser by a cookie, so all the tabs of a browser share
 * one model. A tab that must be an independent instance (the "Open in CoolSolve" links of the
 * model library) calls `enableTabSession()`: it gets its own id, kept in `sessionStorage`
 * (per tab, survives a reload), and sends it with every API call. Tabs that never call it
 * keep the cookie session, so nothing changes for them.
 */

const STORAGE_KEY = 'coolsolve_tab_session';

function readStoredId(): string | null {
  try {
    return sessionStorage.getItem(STORAGE_KEY);
  } catch {
    return null; // storage blocked: fall back to the cookie session
  }
}

/** The id of this tab's own session, or null when the tab uses the shared cookie session. */
export function getTabSessionId(): string | null {
  return readStoredId();
}

/** Give this tab its own server-side session (idempotent). Returns the id. */
export function enableTabSession(): string | null {
  const existing = readStoredId();
  if (existing) return existing;
  // 32 lowercase hex digits: the format the server accepts
  const bytes = crypto.getRandomValues(new Uint8Array(16));
  const id = Array.from(bytes, (b) => b.toString(16).padStart(2, '0')).join('');
  try {
    sessionStorage.setItem(STORAGE_KEY, id);
  } catch {
    return null;
  }
  return id;
}

/** Headers identifying this tab's session (empty for the cookie session). */
export function sessionHeaders(): Record<string, string> {
  const id = getTabSessionId();
  return id ? { 'X-CoolSolve-Session': id } : {};
}

/** URL with the session id, for what cannot set headers (EventSource, downloads). */
export function withSession(url: string): string {
  const id = getTabSessionId();
  if (!id) return url;
  return url + (url.includes('?') ? '&' : '?') + 'sid=' + id;
}
