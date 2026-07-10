/**
 * errors.ts — small helper for safe `catch`-clause message extraction.
 *
 * TypeScript's `useUnknownInCatchVariables` makes caught values `unknown`,
 * which is correct but verbose to narrow at every catch site. This helper
 * reproduces the historical `err.message` convenience without resorting to
 * `catch (err: any)`.
 */

/** Extract a human-readable message from a caught value of unknown type. */
export function toErrMsg(err: unknown): string {
  if (err instanceof Error) return err.message;
  if (typeof err === 'string') return err;
  if (err && typeof err === 'object' && 'message' in err) {
    const m = (err as { message: unknown }).message;
    if (typeof m === 'string') return m;
  }
  return String(err);
}
