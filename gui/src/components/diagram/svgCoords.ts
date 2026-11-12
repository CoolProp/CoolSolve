/**
 * svgCoords.ts — from the pointer to data values on coolplot-core's SVG (docs/coolplot_integration.md §3.2).
 *
 * The root `<svg>` carries the axis limits in display units (`data-x-limits="lo hi"`) and the log flags
 * (`data-x-log="0|1"`); the plot area is `rect.cp-background`, in viewBox units. CoolPlot's
 * `render/web.py` (pointer client) maps the pointer the same way; here the screen transformation of the
 * SVG is used, which also holds when the drawing is scaled with its aspect kept.
 */

import type { DiagramAxes, DiagramAxisUnit, DiagramViewSI } from '../../api/types';

/** The plot area of a drawn scene and its axes. */
export interface PlotFrame {
  /** Plot area, viewBox units */
  x: number;
  y: number;
  width: number;
  height: number;
  /** Axis limits, display units */
  xLimits: [number, number];
  yLimits: [number, number];
  xLog: boolean;
  yLog: boolean;
}

export interface SvgPoint {
  x: number;
  y: number;
}

function limitsAttr(svg: Element, name: string, fallback: [number, number] | undefined): [number, number] | null {
  const parts = (svg.getAttribute(name) ?? '').trim().split(/\s+/).map(Number);
  if (parts.length === 2 && parts.every(Number.isFinite)) return [parts[0], parts[1]];
  return fallback ?? null;
}

function flagAttr(svg: Element, name: string, fallback: boolean): boolean {
  const value = svg.getAttribute(name);
  return value === null ? fallback : value === '1' || value === 'true';
}

/** The plot area and axes of a scene's SVG (the response's axes stand in for missing attributes), or null. */
export function readPlotFrame(svg: SVGSVGElement, axes?: DiagramAxes): PlotFrame | null {
  const bg = svg.querySelector('rect.cp-background');
  if (!bg) return null;
  const attr = (name: string) => parseFloat(bg.getAttribute(name) ?? '');
  const [x, y, width, height] = [attr('x'), attr('y'), attr('width'), attr('height')];
  if (![x, y, width, height].every(Number.isFinite) || width <= 0 || height <= 0) return null;
  const xLimits = limitsAttr(svg, 'data-x-limits', axes?.x_limits);
  const yLimits = limitsAttr(svg, 'data-y-limits', axes?.y_limits);
  if (!xLimits || !yLimits) return null;
  return {
    x, y, width, height, xLimits, yLimits,
    xLog: flagAttr(svg, 'data-x-log', axes?.x_log ?? false),
    yLog: flagAttr(svg, 'data-y-log', axes?.y_log ?? false),
  };
}

/** The viewBox point under a position of the screen (client coordinates), or null. */
export function clientToSvg(svg: SVGSVGElement, clientX: number, clientY: number): SvgPoint | null {
  const ctm = svg.getScreenCTM();
  if (!ctm) return null;
  const p = new DOMPoint(clientX, clientY).matrixTransform(ctm.inverse());
  return { x: p.x, y: p.y };
}

/** Whether a viewBox point lies in the plot area. */
export function inPlot(frame: PlotFrame, p: SvgPoint): boolean {
  return p.x >= frame.x && p.x <= frame.x + frame.width && p.y >= frame.y && p.y <= frame.y + frame.height;
}

/** The point moved into the plot area. */
export function clampToPlot(frame: PlotFrame, p: SvgPoint): SvgPoint {
  return {
    x: Math.min(Math.max(p.x, frame.x), frame.x + frame.width),
    y: Math.min(Math.max(p.y, frame.y), frame.y + frame.height),
  };
}

/** The value at a fraction of an axis (0 at its lower limit, 1 at its upper limit), linear or log. */
function valueAt(fraction: number, [lo, hi]: [number, number], log: boolean): number {
  if (log) {
    const a = Math.log10(lo);
    const b = Math.log10(hi);
    return 10 ** (a + fraction * (b - a));
  }
  return lo + fraction * (hi - lo);
}

/** Display values [x, y] at a viewBox point. */
export function svgToDisplay(frame: PlotFrame, p: SvgPoint): [number, number] {
  const fx = (p.x - frame.x) / frame.width;
  const fy = 1 - (p.y - frame.y) / frame.height;  // SVG y grows downwards
  return [valueAt(fx, frame.xLimits, frame.xLog), valueAt(fy, frame.yLimits, frame.yLog)];
}

/** SI value of a display value: display = SI * scale + offset. */
export function displayToSI(value: number, unit: DiagramAxisUnit | undefined): number {
  if (!unit || !Number.isFinite(unit.scale) || unit.scale === 0) return value;
  return (value - (unit.offset ?? 0)) / unit.scale;
}

function sortedRange(a: number, b: number): [number, number] | null {
  const lo = Math.min(a, b);
  const hi = Math.max(a, b);
  if (!Number.isFinite(lo) || !Number.isFinite(hi) || !(hi > lo)) return null;
  // A span too small to draw (a click, or a zoom far past the precision of the numbers)
  if (hi - lo <= 1e-9 * Math.max(Math.abs(lo), Math.abs(hi))) return null;
  return [lo, hi];
}

/** The zoom window (view_SI) of the rectangle between two viewBox points, clamped to the plot area; null if unusable. */
export function viewFromRectangle(frame: PlotFrame, axes: DiagramAxes, a: SvgPoint, b: SvgPoint): DiagramViewSI | null {
  const [x0, y0] = svgToDisplay(frame, clampToPlot(frame, a));
  const [x1, y1] = svgToDisplay(frame, clampToPlot(frame, b));
  const xs = sortedRange(displayToSI(x0, axes.x_unit), displayToSI(x1, axes.x_unit));
  const ys = sortedRange(displayToSI(y0, axes.y_unit), displayToSI(y1, axes.y_unit));
  if (!xs || !ys) return null;
  if ((frame.xLog && xs[0] <= 0) || (frame.yLog && ys[0] <= 0)) return null;
  return [xs, ys];
}
