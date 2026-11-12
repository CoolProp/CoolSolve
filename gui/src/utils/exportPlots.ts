/**
 * exportPlots.ts — Capture the visible charts as PNG data URIs for the LaTeX report.
 *
 * Returns an array of { name, data } suitable for the /api/v1/latex/compile endpoint.
 *
 *   - The property diagram (coolplot-core's SVG in the Diagrams tab, `.thermo-diagram`) is rasterised
 *     in the browser (Image + canvas, 2×, white background) → "thermo_diagrams.png".
 *   - Plotly charts (every <div class="js-plotly-plot">) go through Plotly.toImage():
 *     inside a ".parametric-plot" ancestor → "parametric_1.png", "parametric_2.png", …;
 *     any other chart → "plot_<n>.png".
 */

import Plotly from 'plotly.js-cartesian-dist-min';
import type { LatexPlotImage } from '../api/types';

/**
 * The `plotly.js-cartesian-dist-min` type definitions omit the imperative
 * `toImage` method that ships with the runtime bundle. Define the minimal
 * slice we call so the cast stays scoped instead of `any`.
 */
interface PlotlyImperative {
  toImage(gd: HTMLDivElement, opts: Record<string, unknown>): Promise<string>;
}
type PlotlyWithImperative = typeof Plotly & PlotlyImperative;

/** Minimal typing for a Plotly-rendered DOM element (has .data attached). */
interface PlotlyGraphDiv extends HTMLDivElement {
  data?: unknown[];
}

/** The diagram's SVG as drawn by DiagramView (not the icons of its buttons). */
const DIAGRAM_SVG_SELECTOR = '.thermo-diagram .diagram-svg-host > svg';

/** Rasterise an SVG element to a PNG data URI at `scale` times its own size, on white. */
async function rasteriseSvg(svg: SVGSVGElement, scale: number): Promise<string> {
  const box = svg.viewBox?.baseVal;
  const width = parseFloat(svg.getAttribute('width') ?? '') || box?.width || svg.getBoundingClientRect().width;
  const height = parseFloat(svg.getAttribute('height') ?? '') || box?.height || svg.getBoundingClientRect().height;
  if (!(width > 0 && height > 0)) throw new Error('the diagram has no size');

  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute('width', String(width));
  clone.setAttribute('height', String(height));
  if (!clone.getAttribute('xmlns')) clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  const text = new XMLSerializer().serializeToString(clone);
  const url = URL.createObjectURL(new Blob([text], { type: 'image/svg+xml;charset=utf-8' }));
  try {
    const image = new Image();
    await new Promise<void>((resolve, reject) => {
      image.onload = () => resolve();
      image.onerror = () => reject(new Error('the diagram SVG could not be read as an image'));
      image.src = url;
    });
    const canvas = document.createElement('canvas');
    canvas.width = Math.round(width * scale);
    canvas.height = Math.round(height * scale);
    const ctx = canvas.getContext('2d');
    if (!ctx) throw new Error('no 2D canvas');
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(image, 0, 0, canvas.width, canvas.height);
    return canvas.toDataURL('image/png');
  } finally {
    URL.revokeObjectURL(url);
  }
}

export async function exportAllPlots(): Promise<LatexPlotImage[]> {
  const plots: LatexPlotImage[] = [];

  // The property diagram, when it is shown (a hidden tab may hold a chart of an older solve)
  const diagram = document.querySelector<SVGSVGElement>(DIAGRAM_SVG_SELECTOR);
  if (diagram && diagram.getBoundingClientRect().width > 0) {
    try {
      plots.push({ name: 'thermo_diagrams.png', data: await rasteriseSvg(diagram, 2) });
    } catch {
      // Skip a diagram that fails to rasterise
    }
  }

  // Plotly attaches the class "js-plotly-plot" to every rendered chart div.
  const chartDivs = document.querySelectorAll<HTMLDivElement>('.js-plotly-plot');

  let parametricIdx = 0;
  let genericIdx = 0;

  for (const div of chartDivs) {
    // Skip charts that are not visible (e.g. tab is hidden, display:none parents)
    if (div.offsetParent === null && div.offsetWidth === 0) continue;
    // Also skip if no data traces
    const gd = div as unknown as PlotlyGraphDiv;
    if (!gd.data || gd.data.length === 0) continue;

    try {
      const dataUrl: string = await (Plotly as PlotlyWithImperative).toImage(gd, {
        format: 'png',
        width: 1200,
        height: 800,
        scale: 2, // 2× for crisp text in the PDF
      });

      // Determine a filename based on ancestry
      let name: string;
      if (div.closest('.parametric-plot')) {
        parametricIdx++;
        name = `parametric_${parametricIdx}.png`;
      } else {
        genericIdx++;
        name = `plot_${genericIdx}.png`;
      }

      plots.push({ name, data: dataUrl });
    } catch {
      // Skip charts that fail to export (e.g. blank canvas)
    }
  }

  return plots;
}
