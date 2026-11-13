/**
 * useDiagramHover — the hover of a drawn diagram (docs/coolplot_integration.md §5, phase 2).
 *
 * Near a state marker: the cycle, the state and its variables with their values. The marker positions
 * come from the scene's `states` (values in SI) through the axes' units and the plot area of the SVG
 * (log axes included), so nothing depends on how coolplot-core draws a marker. Elsewhere in the plot
 * area: the full state under the cursor, from `POST /api/v1/diagram/state`, debounced; a newer request
 * cancels the older (of humid air: T, W, φ, T_wb, h and v at the chart's total pressure; nothing above
 * saturation). Values are shown in the chart's display units (`GET /api/v1/diagram/options`).
 */

import { useEffect, useMemo, useRef, useState } from 'react';
import type { MouseEvent as ReactMouseEvent, RefObject } from 'react';
import { api } from '../../api/client';
import type { DiagramOptionsResponse, DiagramPointState, DiagramSceneInput, DiagramSceneResponse } from '../../api/types';
import { formatValue, loadDiagramOptions, symbolOf, toDisplay, unitTableOf, type UnitTable } from './options';
import { clientToSvg, displayToSI, inPlot, readPlotFrame, svgToDisplay, type PlotFrame } from './svgCoords';

/** Wait this long after the pointer stops before asking for the state under it */
const POINT_DEBOUNCE_MS = 80;
/** A marker is hovered within this distance (CSS px) */
const MARKER_RADIUS_PX = 9;
/** The quantities of the state under the cursor: of a fluid, of humid air (per kg of dry air) */
const POINT_KEYS = ['T', 'p', 'h', 's', 'rho'] as const;
const AIR_POINT_KEYS = ['T', 'W', 'RH', 'Twb', 'h', 'v'] as const;
/** The order of a state's variables in its tooltip */
const KEY_ORDER = ['T', 'W', 'RH', 'Twb', 'Tdp', 'p', 'h', 's', 'rho', 'u', 'v', 'Q'];

/** A tooltip, positioned in the drawing area (CSS px); `flipX`/`flipY` put it left of / above the pointer. */
export interface HoverTip {
  kind: 'state' | 'point';
  left: number;
  top: number;
  flipX: boolean;
  flipY: boolean;
  title?: string;
  lines: string[];
}

/** What is on screen: the scene, its SVG node and the request it answers. */
export interface DrawnScene {
  scene: DiagramSceneResponse;
  node: SVGSVGElement;
  input: DiagramSceneInput;
}

function withUnit(value: number | null | undefined, key: string, units: UnitTable): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return '—';
  const unit = units[key];
  const label = unit && unit.label !== '-' ? ` ${unit.label}` : '';
  return `${formatValue(toDisplay(value, unit))}${label}`;
}

/** The point of the viewBox of a value pair in SI, or null outside the plot area or a log axis. */
function svgPointOf(frame: PlotFrame, scene: DiagramSceneResponse, values: Record<string, number>) {
  const { axes } = scene;
  const x = toDisplay(values[axes.x_key], axes.x_unit);
  const y = toDisplay(values[axes.y_key], axes.y_unit);
  const fraction = (v: number, [lo, hi]: [number, number], log: boolean) =>
    log ? (v > 0 ? (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo)) : NaN) : (v - lo) / (hi - lo);
  const fx = fraction(x, frame.xLimits, frame.xLog);
  const fy = fraction(y, frame.yLimits, frame.yLog);
  if (!(fx >= 0 && fx <= 1 && fy >= 0 && fy <= 1)) return null;
  return { x: frame.x + fx * frame.width, y: frame.y + (1 - fy) * frame.height };
}

export function useDiagramHover(drawn: DrawnScene | null, canvasRef: RefObject<HTMLDivElement>) {
  const [tip, setTip] = useState<HoverTip | null>(null);
  const [options, setOptions] = useState<DiagramOptionsResponse | null>(null);
  const timerRef = useRef<ReturnType<typeof setTimeout> | undefined>(undefined);
  const requestRef = useRef<{ controller: AbortController; seq: number } | null>(null);
  const seqRef = useRef(0);
  const draggingRef = useRef(false);

  useEffect(() => {
    let cancelled = false;
    loadDiagramOptions().then((o) => { if (!cancelled) setOptions(o); }).catch(() => {});
    return () => { cancelled = true; };
  }, []);

  const units = useMemo(() => (drawn ? unitTableOf(options, drawn.input.spec.units) : {}), [options, drawn]);

  const cancel = () => {
    if (timerRef.current !== undefined) clearTimeout(timerRef.current);
    timerRef.current = undefined;
    requestRef.current?.controller.abort();
    requestRef.current = null;
    seqRef.current++;
  };

  // A new scene: positions and units change
  useEffect(() => {
    setTip(null);
    if (timerRef.current !== undefined) clearTimeout(timerRef.current);
    requestRef.current?.controller.abort();
  }, [drawn]);

  // Nothing on its way when the view goes
  useEffect(() => {
    const timer = timerRef;
    const request = requestRef;
    return () => {
      if (timer.current !== undefined) clearTimeout(timer.current);
      request.current?.controller.abort();
    };
  }, []);

  const stateTip = (base: Omit<HoverTip, 'kind' | 'lines' | 'title'>, index: number): HoverTip | null => {
    if (!drawn) return null;
    const state = drawn.scene.states[index];
    const cycle = drawn.input.cycles[state.cycle]?.name ?? `Cycle ${state.cycle + 1}`;
    const rank = (key: string) => (KEY_ORDER.indexOf(key) + KEY_ORDER.length + 1) % (KEY_ORDER.length + 1);
    const lines = Object.entries(state.vars)
      .sort(([a], [b]) => rank(a) - rank(b))
      .map(([key, name]) => `${name} = ${withUnit(state.values_SI[key], key, units)}`);
    return { ...base, kind: 'state', title: `${cycle}: state ${state.number}${state.label ? ` (${state.label})` : ''}`, lines };
  };

  const pointTip = (base: Omit<HoverTip, 'kind' | 'lines' | 'title'>, state: DiagramPointState): HoverTip => {
    const humidAir = state.W !== undefined;  // the state of humid air: its own quantities, no quality
    const lines = (humidAir ? AIR_POINT_KEYS : POINT_KEYS).map((key) => `${symbolOf(key)} = ${withUnit(state[key], key, units)}`);
    const Q = state.Q;
    if (Q !== null && Q !== undefined && Q >= 0 && Q <= 1) lines.push(`x = ${formatValue(Q)}`);
    return { ...base, kind: 'point', lines };
  };

  const onMouseMove = (e: ReactMouseEvent<HTMLDivElement>) => {
    const canvas = canvasRef.current;
    if (draggingRef.current || !drawn || !canvas) return;
    const svg = drawn.node;
    const frame = readPlotFrame(svg, drawn.scene.axes);
    const p = clientToSvg(svg, e.clientX, e.clientY);
    if (!frame || !p || !inPlot(frame, p)) {
      cancel();
      setTip(null);
      return;
    }
    const box = canvas.getBoundingClientRect();
    const left = e.clientX - box.left;
    const top = e.clientY - box.top;
    const base = { left, top, flipX: left > box.width * 0.6, flipY: top > box.height * 0.6 };

    // The nearest state marker within reach
    const ctm = svg.getScreenCTM();
    let nearest = -1;
    let best = MARKER_RADIUS_PX;
    if (ctm) {
      drawn.scene.states.forEach((state, i) => {
        const at = svgPointOf(frame, drawn.scene, state.values_SI);
        if (!at) return;
        const screen = new DOMPoint(at.x, at.y).matrixTransform(ctm);
        const d = Math.hypot(screen.x - e.clientX, screen.y - e.clientY);
        if (d <= best) {
          best = d;
          nearest = i;
        }
      });
    }
    if (nearest >= 0) {
      cancel();
      setTip(stateTip(base, nearest));
      return;
    }

    // Elsewhere: the state under the cursor, once the pointer rests
    setTip((t) => (t?.kind === 'state' ? null : t));
    const [xd, yd] = svgToDisplay(frame, p);
    const body = {
      fluid: drawn.input.fluid,
      ...(drawn.input.fractions ? { fractions: drawn.input.fractions } : {}),
      ...(drawn.input.p_Pa !== undefined ? { p_Pa: drawn.input.p_Pa } : {}),
      diagram: drawn.input.spec.diagram,
      x_SI: displayToSI(xd, drawn.scene.axes.x_unit),
      y_SI: displayToSI(yd, drawn.scene.axes.y_unit),
    };
    if (timerRef.current !== undefined) clearTimeout(timerRef.current);
    timerRef.current = setTimeout(() => {
      timerRef.current = undefined;
      requestRef.current?.controller.abort();
      const controller = new AbortController();
      const seq = ++seqRef.current;
      requestRef.current = { controller, seq };
      api.getDiagramPointState(body, controller.signal)
        .then((res) => {
          if (seq !== seqRef.current) return;
          setTip(res.state ? pointTip(base, res.state) : null);
        })
        .catch(() => {
          if (seq === seqRef.current) setTip(null);
        });
    }, POINT_DEBOUNCE_MS);
  };

  const onMouseLeave = () => {
    cancel();
    setTip(null);
  };

  /** A zoom drag starts (no hover until it ends) or ends. */
  const setDragging = (dragging: boolean) => {
    draggingRef.current = dragging;
    if (dragging) {
      cancel();
      setTip(null);
    }
  };

  return { tip, onMouseMove, onMouseLeave, setDragging };
}
