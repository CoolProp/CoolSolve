/**
 * exportPlots.ts — Capture all visible Plotly charts as PNG data URIs.
 *
 * Uses Plotly.toImage() on every <div class="js-plotly-plot"> in the DOM.
 * Returns an array of { name, data } suitable for the /api/v1/latex/compile
 * endpoint.
 *
 * Naming convention:
 *   - The first chart inside a ".thermo-diagram" ancestor → "thermo_diagrams.png"
 *   - Charts inside a ".parametric-plot" ancestor → "parametric_1.png", "parametric_2.png", …
 *   - Any other chart → "plot_<n>.png"
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

export async function exportAllPlots(): Promise<LatexPlotImage[]> {
  const plots: LatexPlotImage[] = [];

  // Plotly attaches the class "js-plotly-plot" to every rendered chart div.
  const chartDivs = document.querySelectorAll<HTMLDivElement>('.js-plotly-plot');
  if (chartDivs.length === 0) return plots;

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
      if (div.closest('.thermo-diagram')) {
        name = 'thermo_diagrams.png';
      } else if (div.closest('.parametric-plot')) {
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
