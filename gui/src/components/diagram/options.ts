/**
 * options.ts — what the server draws (`GET /api/v1/diagram/options`), loaded once and shared by the
 * Diagrams tab and every DiagramView, and the display units of the charts: the model's (for the hover and
 * the range fields).
 */

import { api } from '../../api/client';
import type { DiagramAxisUnit, DiagramOptionsResponse } from '../../api/types';

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

/** The display units of the model (`model_units`), those of every chart the GUI draws. */
export function modelUnitTable(options: DiagramOptionsResponse): UnitTable {
  const table: UnitTable = {};
  for (const [key, entry] of Object.entries(options.model_units.units)) {
    if (!Array.isArray(entry) || entry.length < 3) continue;
    const [label, scale, offset] = entry as [string, number, number];
    if (typeof scale !== 'number' || typeof offset !== 'number') continue;
    table[key] = { label: UNICODE_UNITS[label] ?? label, scale, offset };
  }
  return table;
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

/** How a quantity key is written (as coolplot-core writes the symbols with "unicode"). */
export const SYMBOLS: Record<string, string> = { rho: 'ρ', Q: 'x', RH: 'φ', Twb: 'T_wb', Tdp: 'T_dp' };
export const symbolOf = (key: string) => SYMBOLS[key] ?? key;

/** How a chart writes the quantity of one of its axes: the Mollier chart writes its humidity ratio x and its temperature scale T. */
export function axisSymbol(type: { y: string } | undefined, key: string): string {
  if (type?.y === 'Thx' && key === 'W') return 'x';
  return symbolOf(key === 'Thx' ? 'T' : key);
}
