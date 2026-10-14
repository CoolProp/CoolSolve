/**
 * library.ts — open a model of the embedded CoolSolve Library into the current session.
 *
 * A library model is a whole folder (`.eescode`, `.initials`, `coolsolve.conf`, lookup tables,
 * `.sol` baseline): the server loads it all into the session, and the editor state is then
 * read back through the same calls as for any opened model.
 */

import { api } from '../api/client';
import { useModelStore } from '../stores/modelStore';
import { toErrMsg } from './errors';

/** Address of the library dashboard (all models, blocked ones included). */
export const LIBRARY_DASHBOARD_URL = '/library/';

/** Dashboard page of one model (summary, verification, README, gaps). */
export const libraryModelPageUrl = (id: string) => `/library/#${encodeURIComponent(id)}`;

/**
 * Open a library model (by permanent ID or name, optionally one of its variant files).
 * `allowBack` arms the Back button, which restores the model that was open before.
 * Returns true on success; errors are reported in the console panel.
 */
export async function openLibraryModel(id: string, file?: string, allowBack = true): Promise<boolean> {
  const { loadFile, addConsoleLine, setCanGoBack } = useModelStore.getState();
  try {
    const res = await api.openLibraryModel(id, file);
    const ees = await api.getEescode();
    const init = await api.getInitials();
    const solRes = await api.getSol();
    const confRes = await api.getConf();
    setCanGoBack(allowBack);
    loadFile('', ees.content, init.content, solRes.content, confRes.content, res.modelName, res.libraryId);
    addConsoleLine(`>>> Opened library model: ${res.title} (${res.libraryId})`);
    if (res.status !== 'verified') {
      addConsoleLine(`>>> Note: this model is "${res.status}" in the library; it may not solve as is.`);
    }
    return true;
  } catch (err: unknown) {
    addConsoleLine(`>>> ERROR opening library model ${id}: ${toErrMsg(err)}`);
    return false;
  }
}
