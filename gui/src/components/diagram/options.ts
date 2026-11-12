/**
 * options.ts — what the server draws (`GET /api/v1/diagram/options`), loaded once and shared by the
 * Diagrams tab and every DiagramView, and the display units of a spec (for the hover).
 */

import { api } from '../../api/client';
import type { DiagramAxisUnit, DiagramOptionsResponse, DiagramUnits } from '../../api/types';

let cached: Promise<DiagramOptionsResponse> | null = null;

/** The diagram options; a failed load is tried again at the next call. */
export function loadDiagramOptions(): Promise<DiagramOptionsResponse> {
  if (!cached) {
    cached = api.getDiagramOptions().catch((err: unknown) => {
      cached = null;
      throw err;
    });
  }
  return cached;
}

/** Units written with their own characters, as coolplot-core writes them with symbols "unicode". */
const UNICODE_UNITS: Record<string, string> = {
  'deg C': '°C', 'deg F': '°F', R: '°R', 'Btu/lb/R': 'Btu/lb/°R',
  'kg/m3': 'kg/m³', 'm3/kg': 'm³/kg', 'lb/ft3': 'lb/ft³', 'ft3/lb': 'ft³/lb',
};

/** Display units by quantity key: display = SI * scale + offset. */
export type UnitTable = Record<string, DiagramAxisUnit>;

function tableOf(raw: Record<string, unknown> | undefined): UnitTable {
  const table: UnitTable = {};
  for (const [key, entry] of Object.entries(raw ?? {})) {
    if (!Array.isArray(entry) || entry.length < 3) continue;
    const [label, scale, offset] = entry as [string, number, number];
    if (typeof scale !== 'number' || typeof offset !== 'number') continue;
    table[key] = { label: UNICODE_UNITS[label] ?? label, scale, offset };
  }
  return table;
}

/** The display units of a spec's `units` (a preset name, or a unit record). */
export function unitTableOf(options: DiagramOptionsResponse | null, units: DiagramUnits): UnitTable {
  const name = typeof units === 'string' ? units : units.name;
  const offered = options?.unit_systems.find((u) => u.name === name);
  if (offered) return tableOf(offered.units);
  return typeof units === 'string' ? {} : tableOf(units.units);
}

/** An SI value in display units, or null when it has no unit in the table. */
export function toDisplay(value: number, unit: DiagramAxisUnit | undefined): number {
  return unit ? value * unit.scale + unit.offset : value;
}

/** A value with four significant digits, in exponent form when very large or small. */
export function formatValue(value: number): string {
  if (!Number.isFinite(value)) return '—';
  const size = Math.abs(value);
  if (size !== 0 && (size >= 1e6 || size < 1e-3)) return value.toExponential(3).replace(/\.?0+e/, 'e');
  return String(parseFloat(value.toPrecision(4)));
}

/** How a quantity key is written. */
export const SYMBOLS: Record<string, string> = { rho: 'ρ', Q: 'x' };
export const symbolOf = (key: string) => SYMBOLS[key] ?? key;
