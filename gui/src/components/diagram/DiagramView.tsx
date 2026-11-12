/**
 * DiagramView — a property diagram of the session's model, drawn by the server (coolplot-core's SVG).
 *
 * Given what to draw (fluid, CoolPlot spec, cycle definitions), it asks `POST /api/v1/diagram/scene`
 * for a figure the size of its container and shows the SVG inline, as DOM nodes (so that a hover layer
 * can find the items later), scaled to the container with its aspect kept. Requests are debounced, a
 * newer one cancels the older, and the labels of the scene on screen are sent back when only the size or
 * the zoom changed. Drag a rectangle on the plot area to zoom (a new `view_SI`); double-click or the
 * button resets it. Independent of the Diagrams tab's controls, so that diagrams can be shown elsewhere.
 */

import { useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react';
import type { MouseEvent as ReactMouseEvent } from 'react';
import { AlertTriangle, ZoomOut } from 'lucide-react';
import { api } from '../../api/client';
import type { CoolPlotSpec, DiagramSceneInput, DiagramSceneResponse, DiagramViewSI } from '../../api/types';
import { toErrMsg } from '../../utils/errors';
import { clientToSvg, inPlot, readPlotFrame, viewFromRectangle } from './svgCoords';

/** CSS px to pt */
const PT_PER_PX = 0.75;
/** Wait this long after the last change (resize, edit, solve) before asking for a new scene */
const DEBOUNCE_MS = 150;
/** No request below this size (a hidden tab measures 0) */
const MIN_SIZE_PX = 120;
/** A smaller drag is a click, not a zoom */
const MIN_BAND_PX = 6;
const SVG_NS = 'http://www.w3.org/2000/svg';

export interface DiagramViewProps {
  /** What to draw (the spec may hold a zoom, `view_SI`); null draws nothing. */
  input: DiagramSceneInput | null;
  /** Changes when the session's solution changes, so that the same input is drawn again with it. */
  revision?: unknown;
  /** Controlled zoom: called with the new view (null: the whole chart), which the owner puts in the
   *  input's spec. Without it the zoom is kept here, for the fluid and diagram type it was made on. */
  onViewChange?: (view: DiagramViewSI | null) => void;
  /** Called with every scene shown. */
  onScene?: (scene: DiagramSceneResponse) => void;
  className?: string;
}

interface Drawn {
  scene: DiagramSceneResponse;
  node: SVGSVGElement;
}

interface Shown {
  contentKey: string;
  requestKey: string;
  revision: unknown;
  labels: Record<string, unknown>[];
}

interface Pending {
  requestKey: string;
  revision: unknown;
  controller: AbortController;
}

let instanceCount = 0;

/**
 * The SVG document as a node of this page. Its ids (the clip path of the plot area) get the suffix of
 * this view, so that two diagrams on one page do not clip each other.
 */
function svgNode(text: string, suffix: string): SVGSVGElement {
  const doc = new DOMParser().parseFromString(text, 'image/svg+xml');
  let root: Element | null = doc.documentElement;
  if (doc.getElementsByTagName('parsererror').length > 0 || !root || root.localName !== 'svg') {
    throw new Error('the server sent no SVG document');
  }
  if (root.namespaceURI !== SVG_NS) {
    // Without xmlns the XML parser makes plain elements: the HTML parser puts <svg> in its namespace
    const template = document.createElement('template');
    template.innerHTML = text;
    root = template.content.querySelector('svg');
    if (!root) throw new Error('the server sent no SVG document');
  }
  const ids = new Map<string, string>();
  root.querySelectorAll('[id]').forEach((el) => {
    const renamed = `${el.id}-${suffix}`;
    ids.set(el.id, renamed);
    el.id = renamed;
  });
  if (ids.size > 0) {
    const rename = (value: string) =>
      value.replace(/url\(#([^)]+)\)/g, (match, id: string) => (ids.has(id) ? `url(#${ids.get(id)})` : match));
    for (const el of [root, ...Array.from(root.querySelectorAll('*'))]) {
      for (const attr of Array.from(el.attributes)) {
        if (!attr.value.includes('#')) continue;
        let value = rename(attr.value);
        if ((attr.localName === 'href') && value.startsWith('#') && ids.has(value.slice(1))) {
          value = `#${ids.get(value.slice(1))}`;
        }
        if (value !== attr.value) el.setAttributeNS(attr.namespaceURI, attr.name, value);
      }
    }
  }
  return document.importNode(root, true) as unknown as SVGSVGElement;
}

const round1 = (v: number) => Math.round(v * 10) / 10;

function specWithView(spec: CoolPlotSpec, view: DiagramViewSI | null): CoolPlotSpec {
  const { view_SI: _view, ...rest } = spec;
  return view ? { ...rest, view_SI: view } : rest;
}

export default function DiagramView({ input, revision, onViewChange, onScene, className }: DiagramViewProps) {
  const [suffix] = useState(() => `cs-diagram-${++instanceCount}`);
  const canvasRef = useRef<HTMLDivElement>(null);
  const hostRef = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState<[number, number]>([0, 0]);
  const [drawn, setDrawn] = useState<Drawn | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [band, setBand] = useState<{ left: number; top: number; width: number; height: number } | null>(null);
  // The zoom when the owner does not control it, for one fluid and diagram type
  const [localView, setLocalView] = useState<{ chart: string; view: DiagramViewSI | null } | null>(null);

  const shownRef = useRef<Shown | null>(null);
  const pendingRef = useRef<Pending | null>(null);
  const seqRef = useRef(0);
  const onSceneRef = useRef(onScene);
  useEffect(() => {
    onSceneRef.current = onScene;
  });

  const chart = input ? `${input.fluid}\n${input.spec.diagram}` : '';
  const view: DiagramViewSI | null = onViewChange || !localView || localView.chart !== chart
    ? input?.spec.view_SI ?? null
    : localView.view;

  const spec = useMemo(() => (input ? specWithView(input.spec, view) : null), [input, view]);
  // Everything but the zoom and the size: the labels of the scene shown are reused while it stays the same
  const contentKey = useMemo(
    () => (input ? JSON.stringify({ fluid: input.fluid, spec: specWithView(input.spec, null), cycles: input.cycles }) : ''),
    [input],
  );

  // The size of the drawing area
  useEffect(() => {
    const el = canvasRef.current;
    if (!el) return;
    const observer = new ResizeObserver((entries) => {
      const box = entries[0].contentRect;
      const w = Math.floor(box.width);
      const h = Math.floor(box.height);
      setSize((s) => (s[0] === w && s[1] === h ? s : [w, h]));
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  // Ask for a scene when what is drawn, the zoom, the size or the solution changes
  useEffect(() => {
    if (!input || !spec || size[0] < MIN_SIZE_PX || size[1] < MIN_SIZE_PX) return;
    const figure: [number, number] = [round1(size[0] * PT_PER_PX), round1(size[1] * PT_PER_PX)];
    const requestKey = `${contentKey}|${JSON.stringify(spec.view_SI ?? null)}|${figure[0]}x${figure[1]}`;
    const matches = (k: { requestKey: string; revision: unknown } | null) =>
      !!k && k.requestKey === requestKey && Object.is(k.revision, revision);
    if (matches(pendingRef.current)) return;
    if (matches(shownRef.current)) {
      // Back to the scene on screen: drop a request on its way
      if (pendingRef.current) {
        pendingRef.current.controller.abort();
        pendingRef.current = null;
        seqRef.current++;
        setLoading(false);
      }
      return;
    }
    const timer = setTimeout(() => {
      pendingRef.current?.controller.abort();
      const controller = new AbortController();
      const seq = ++seqRef.current;
      pendingRef.current = { requestKey, revision, controller };
      const shown = shownRef.current;
      const sameContent = !!shown && shown.contentKey === contentKey && Object.is(shown.revision, revision);
      setLoading(true);
      api.getDiagramScene(
        {
          fluid: input.fluid,
          spec,
          cycles: input.cycles,
          figure_size_pt: figure,
          ...(sameContent && shown.labels.length > 0 ? { previous_labels: shown.labels } : {}),
        },
        controller.signal,
      )
        .then((scene) => {
          if (seq !== seqRef.current) return;
          const node = svgNode(scene.svg, suffix);
          shownRef.current = { contentKey, requestKey, revision, labels: scene.labels ?? [] };
          setDrawn({ scene, node });
          setError('');
          onSceneRef.current?.(scene);
        })
        .catch((err: unknown) => {
          if (seq !== seqRef.current || controller.signal.aborted) return;
          shownRef.current = null;
          setDrawn(null);
          setError(toErrMsg(err));
        })
        .finally(() => {
          if (seq !== seqRef.current) return;
          pendingRef.current = null;
          setLoading(false);
        });
    }, DEBOUNCE_MS);
    return () => clearTimeout(timer);
  }, [input, spec, contentKey, size, revision, suffix]);

  // A request still on its way when the view goes
  useEffect(() => {
    const pending = pendingRef;
    const seq = seqRef;
    return () => {
      pending.current?.controller.abort();
      seq.current++;
    };
  }, []);

  // The SVG as DOM nodes
  useLayoutEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    if (drawn) host.replaceChildren(drawn.node);
    else host.replaceChildren();
  }, [drawn]);

  const applyView = (next: DiagramViewSI | null) => {
    if (onViewChange) onViewChange(next);
    else setLocalView({ chart, view: next });
  };

  // Rubber-band zoom on the plot area
  const onMouseDown = (e: ReactMouseEvent<HTMLDivElement>) => {
    if (e.button !== 0 || !drawn || !canvasRef.current) return;
    const svg = drawn.node;
    const axes = drawn.scene.axes;
    const frame = readPlotFrame(svg, axes);
    const start = clientToSvg(svg, e.clientX, e.clientY);
    const plot = svg.querySelector('rect.cp-background');
    if (!frame || !start || !plot || !inPlot(frame, start)) return;
    e.preventDefault();
    const canvasBox = canvasRef.current.getBoundingClientRect();
    const plotBox = plot.getBoundingClientRect();
    const clampX = (x: number) => Math.min(Math.max(x, plotBox.left), plotBox.right);
    const clampY = (y: number) => Math.min(Math.max(y, plotBox.top), plotBox.bottom);
    const x0 = clampX(e.clientX);
    const y0 = clampY(e.clientY);
    const move = (ev: MouseEvent) => {
      const x1 = clampX(ev.clientX);
      const y1 = clampY(ev.clientY);
      setBand({
        left: Math.min(x0, x1) - canvasBox.left,
        top: Math.min(y0, y1) - canvasBox.top,
        width: Math.abs(x1 - x0),
        height: Math.abs(y1 - y0),
      });
    };
    const up = (ev: MouseEvent) => {
      window.removeEventListener('mousemove', move);
      window.removeEventListener('mouseup', up);
      setBand(null);
      if (Math.abs(clampX(ev.clientX) - x0) < MIN_BAND_PX || Math.abs(clampY(ev.clientY) - y0) < MIN_BAND_PX) return;
      const end = clientToSvg(svg, ev.clientX, ev.clientY);
      const next = end ? viewFromRectangle(frame, axes, start, end) : null;
      if (next) applyView(next);
    };
    window.addEventListener('mousemove', move);
    window.addEventListener('mouseup', up);
  };

  const scene = drawn?.scene;
  const warnings = scene?.diagnostics?.warnings ?? [];
  const timing = scene?.timing_ms;

  return (
    <div className={`diagram-view${className ? ` ${className}` : ''}`}>
      <div
        ref={canvasRef}
        className={`diagram-canvas${drawn ? ' zoomable' : ''}`}
        onMouseDown={onMouseDown}
        onDoubleClick={() => { if (view) applyView(null); }}
      >
        <div ref={hostRef} className="diagram-svg-host" />
        {band && <div className="diagram-band" style={band} />}
        {view && drawn && (
          <button
            className="diagram-reset"
            title="Reset the zoom (or double-click the chart)"
            onClick={() => applyView(null)}
            onMouseDown={(e) => e.stopPropagation()}
          >
            <ZoomOut size={14} />
          </button>
        )}
        {!drawn && !loading && !error && (
          <div className="diagram-empty">{input ? '' : 'Nothing to draw'}</div>
        )}
      </div>
      <div className="diagram-status">
        {loading && <span className="spinner diagram-spinner" />}
        {error ? (
          <span className="diagram-error" title={error}>Diagram not drawn: {error}</span>
        ) : scene ? (
          <>
            {warnings.length > 0 && (
              <span className="diagram-warnings" title={warnings.join('\n')}>
                <AlertTriangle size={11} /> {warnings.length} warning{warnings.length > 1 ? 's' : ''}
              </span>
            )}
            {timing && (
              <span title="Server time: lines, composition, SVG">
                {Math.round(timing.curves)} + {Math.round(timing.compose)} + {Math.round(timing.svg)} ms
              </span>
            )}
            <span className="diagram-hint">Drag to zoom{view ? ', double-click to reset' : ''}</span>
          </>
        ) : null}
      </div>
    </div>
  );
}
