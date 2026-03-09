// Use the cartesian plotly.js distribution (scatter, bar, heatmap, contour, etc.)
// This avoids the massive full plotly.js bundle while supporting 2D plots.
import Plotly from 'plotly.js-cartesian-dist-min';
import createPlotlyComponent from 'react-plotly.js/factory';

const Plot = createPlotlyComponent(Plotly);
export default Plot;
