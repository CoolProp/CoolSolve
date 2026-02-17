// Use the minimal plotly.js distribution (scatter charts only, ~800KB vs ~5MB)
// This avoids the massive full plotly.js bundle.
import Plotly from 'plotly.js-basic-dist-min';
import createPlotlyComponent from 'react-plotly.js/factory';

const Plot = createPlotlyComponent(Plotly);
export default Plot;
