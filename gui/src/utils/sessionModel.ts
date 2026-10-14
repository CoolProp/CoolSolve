/**
 * sessionModel.ts — show in the editor the model held by this tab's server session.
 *
 * The server session holds the reference state of the model (files, lookup tables, last
 * result). After an operation that replaces it (opening a library model or the file given on
 * the command line, Back, ZIP upload) and on the reload of a tab that has its own session,
 * the editor state is read back from it by this one function.
 */

import { api } from '../api/client';
import { useModelStore } from '../stores/modelStore';

/**
 * Load the session's model into the store: its files, name and library ID, its lookup tables
 * (the Lookup Tables panel does not fetch them by itself) and, when the model was solved, its
 * result and integral trajectory.
 */
export async function showSessionModel(): Promise<void> {
  const [ees, init, sol, conf, tables] = await Promise.all([
    api.getEescode(),
    api.getInitials(),
    api.getSol(),
    api.getConf(),
    api.getTables(),
  ]);
  const store = useModelStore.getState();
  store.loadFile(ees.filePath, ees.content, init.content, sol.content, conf.content, ees.modelName, ees.libraryId);
  store.setLookupTables(tables.tables);

  const [result, integral, integralCSV] = await Promise.all([
    api.getSolveResult().catch(() => null),  // 404 when the model was not solved
    api.getIntegralResult().catch(() => null),
    api.getIntegralCSV().catch(() => null),
  ]);
  if (result) store.setSolveResult(result);
  store.setIntegralTable(integral?.integralTable ?? null);
  store.setIntegralCSV(integralCSV ?? '');
}
