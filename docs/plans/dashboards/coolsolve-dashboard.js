/*
 * CoolSolve dashboard runtime — PROTOTYPE
 *
 * A dependency-free subset of the v1 dashboard contract planned in
 * docs/plans/dashboards.md, written so that the examples of this folder run
 * standalone (double-click the .html file, no server, no network). The real
 * runtime is planned as gui/src/dashboard/ (TypeScript, built by Vite, served
 * by the GUI and inlined in standalone exports); this file fixes the authoring
 * contract that it must keep:
 *
 *   <cs-stage width height>   design-coordinate canvas, scaled as one piece
 *   <cs-value name|expr>      a value (or expression) with its unit
 *   <cs-input name>           an input: a `name = constant` line of the model
 *   <cs-button action>        solve | update-guesses | stop | reset-inputs
 *   <cs-status>               solve state of the data shown
 *   <cs-plot ref | x y>       a saved plot, or a plot of integral columns
 *   <cs-player>               animation controls of the integral table
 *   data-cs-*                 bindings on any element, SVG included
 *   CoolSolve.onFrame(fn)     custom drawing, called for every frame
 *
 * Data: a host page (the CoolSolve GUI) posts it with postMessage; standalone,
 * it comes from window.COOLSOLVE_DATA, from <script type="application/json"
 * id="coolsolve-data">, or from the results file named by
 * <cs-dashboard results="..."> (default results.js), loaded as a script so that
 * it also works from file://.
 */
(function () {
  'use strict';

  const CS = {};
  window.CoolSolve = CS;

  let data = null;                 // the data contract (schema coolsolve-dashboard-data/1)
  let hosted = false;              // true once a CoolSolve host has answered
  let lookup = new Map();          // lower-case name -> variable record
  let columns = new Map();         // lower-case name -> integral column
  let integralVar = null;          // lower-case name of the integration variable
  const updaters = [];             // functions called at every frame
  const frameHooks = [];           // CoolSolve.onFrame callbacks
  const problems = new Set();      // binding errors shown in the problem badge
  const player = { t: NaN, t0: NaN, t1: NaN, playing: false, speed: 1, loop: true, duration: 10 };

  // ---------------------------------------------------------------- units --
  // Display units: a scale from the model unit (CoolSolve SI: Pa, J, W, kg, m).
  const UNIT_SCALE = {
    'Pa>kPa': 1e-3, 'Pa>bar': 1e-5, 'Pa>MPa': 1e-6,
    'J/kg>kJ/kg': 1e-3, 'J/kg-K>kJ/kg-K': 1e-3, 'J>kJ': 1e-3, 'J>MJ': 1e-6,
    'W>kW': 1e-3, 'W>MW': 1e-6, 'kg/s>g/s': 1e3, 'kg>g': 1e3, 'kg>mg': 1e6,
    'm^3>cm3': 1e6, 'm^3>L': 1e3, 'm>mm': 1e3, 'm>cm': 1e2, '->%': 100,
  };

  function prettyUnit(u) {
    if (!u || u === '-') return '';
    return u.replace(/\^2/g, '²').replace(/\^3/g, '³').replace(/cm3/g, 'cm³');
  }

  /** Scale that converts a value from the model unit `from` to the display unit `to`. */
  function unitScale(from, to, name) {
    if (!to || to === from) return 1;
    const s = UNIT_SCALE[(from || '-') + '>' + to];
    if (s === undefined) {
      problems.add(`unit "${to}" of ${name}: no conversion from "${from || 'no unit'}" (use scale="…")`);
      return 1;
    }
    return s;
  }

  // ---------------------------------------------------------- expressions --
  // A small safe expression language for bindings (no eval): numbers, variable
  // names (case-insensitive, h[3] allowed), + - * / ^, comparisons, && || !,
  // parentheses and a few functions. Trigonometric functions take degrees, as
  // most CoolSolve models do.
  const FUNCS = {
    abs: Math.abs, sqrt: Math.sqrt, exp: Math.exp, ln: Math.log, log10: Math.log10,
    min: Math.min, max: Math.max, floor: Math.floor, round: Math.round,
    sin: (d) => Math.sin(d * Math.PI / 180), cos: (d) => Math.cos(d * Math.PI / 180),
    tan: (d) => Math.tan(d * Math.PI / 180),
    clamp: (x, a, b) => Math.min(Math.max(x, a), b),
  };

  function compile(src) {
    const toks = src.match(/\s*(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+|[A-Za-z_][A-Za-z0-9_]*(?:\[\d+\])?|&&|\|\||[<>!=]=|[-+*/^(),<>!])/g) || [];
    if (toks.join('').replace(/\s/g, '') !== src.replace(/\s/g, '')) throw new Error(`cannot read "${src}"`);
    let i = 0;
    const names = new Set();
    const peek = () => (toks[i] || '').trim();
    const next = () => (toks[i++] || '').trim();
    const expect = (t) => { if (next() !== t) throw new Error(`"${t}" expected in "${src}"`); };
    const bin = (sub, ops) => () => {
      let a = sub();
      while (ops[peek()]) { const f = ops[next()]; const l = a, r = sub(); a = (c) => f(l(c), r(c)); }
      return a;
    };
    let orExpr;
    const atom = () => {
      const t = next();
      if (t === '(') { const e = orExpr(); expect(')'); return e; }
      if (/^[\d.]/.test(t)) { const v = parseFloat(t); return () => v; }
      if (/^[A-Za-z_]/.test(t)) {
        if (peek() === '(') {
          const f = FUNCS[t.toLowerCase()];
          if (!f) throw new Error(`unknown function ${t}() in "${src}"`);
          next();
          const args = [];
          if (peek() !== ')') { args.push(orExpr()); while (peek() === ',') { next(); args.push(orExpr()); } }
          expect(')');
          return (c) => f(...args.map((a) => a(c)));
        }
        if (t.toLowerCase() === 'pi') return () => Math.PI;
        names.add(t);
        return (c) => c(t);
      }
      throw new Error(`unexpected "${t}" in "${src}"`);
    };
    const unary = () => {
      if (peek() === '-') { next(); const e = unary(); return (c) => -e(c); }
      if (peek() === '!') { next(); const e = unary(); return (c) => (e(c) ? 0 : 1); }
      const a = atom();
      if (peek() === '^') { next(); const b = unary(); return (c) => Math.pow(a(c), b(c)); }
      return a;
    };
    const mul = bin(unary, { '*': (a, b) => a * b, '/': (a, b) => a / b });
    const add = bin(mul, { '+': (a, b) => a + b, '-': (a, b) => a - b });
    const cmp = bin(add, {
      '<': (a, b) => +(a < b), '>': (a, b) => +(a > b), '<=': (a, b) => +(a <= b),
      '>=': (a, b) => +(a >= b), '==': (a, b) => +(a === b), '!=': (a, b) => +(a !== b),
    });
    const and = bin(cmp, { '&&': (a, b) => +(a && b) });
    orExpr = bin(and, { '||': (a, b) => +(a || b) });
    const fn = orExpr();
    if (i < toks.length) throw new Error(`unexpected "${peek()}" in "${src}"`);
    return { fn, names: [...names] };
  }

  /** Compiles a binding expression, recording unknown names as problems. */
  function binding(src, where) {
    try {
      const { fn, names } = compile(src);
      for (const n of names) if (!known(n)) problems.add(`${where}: unknown variable "${n}"`);
      return () => fn(value);
    } catch (e) {
      problems.add(`${where}: ${e.message}`);
      return () => NaN;
    }
  }

  // --------------------------------------------------------------- values --
  function known(name) {
    const k = name.toLowerCase();
    return lookup.has(k) || columns.has(k);
  }

  function unitsOf(name) {
    const k = name.toLowerCase();
    if (lookup.has(k)) return lookup.get(k).units || '';
    return (data && data.integral && data.integral.units && data.integral.units[columnName(k)]) || '';
  }

  function columnName(k) {
    return Object.keys(data.integral.columns).find((c) => c.toLowerCase() === k);
  }

  /** Linear interpolation of an integral column at the integration variable t. */
  function interpolate(col, t) {
    const x = columns.get(integralVar);
    const n = x.length;
    if (!(t > x[0])) return col[0];
    if (t >= x[n - 1]) return col[n - 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (x[mid] <= t) lo = mid; else hi = mid; }
    const w = (t - x[lo]) / (x[hi] - x[lo]);
    return col[lo] + w * (col[hi] - col[lo]);
  }

  /**
   * Value of a variable at the current frame: an integral column at the player
   * position, otherwise the solved value. NaN when unknown.
   */
  function value(name) {
    const k = name.toLowerCase();
    if (columns.has(k) && !isNaN(player.t)) return interpolate(columns.get(k), player.t);
    const v = lookup.get(k);
    return v ? v.value : NaN;
  }

  function format(v, decimals) {
    if (typeof v !== 'number' || !isFinite(v)) return '—';
    if (decimals !== null && decimals !== undefined && decimals !== '') return v.toFixed(+decimals);
    const a = Math.abs(v);
    if (a !== 0 && (a >= 1e6 || a < 1e-3)) return v.toExponential(3);
    return String(parseFloat(v.toPrecision(4)));
  }

  // ------------------------------------------------------------- colours --
  function parseColormap(spec) {
    return spec.split(',').map((s) => {
      const [v, c] = s.split(':');
      const h = c.trim().replace('#', '');
      return [parseFloat(v), [0, 2, 4].map((j) => parseInt(h.substr(j, 2), 16))];
    }).sort((a, b) => a[0] - b[0]);
  }

  function colormap(stops, v) {
    if (!isFinite(v)) return 'transparent';
    if (v <= stops[0][0]) return rgb(stops[0][1]);
    for (let j = 1; j < stops.length; j++) {
      if (v <= stops[j][0]) {
        const w = (v - stops[j - 1][0]) / (stops[j][0] - stops[j - 1][0]);
        return rgb(stops[j - 1][1].map((c, k) => c + w * (stops[j][1][k] - c)));
      }
    }
    return rgb(stops[stops.length - 1][1]);
  }

  const rgb = (c) => `rgb(${c.map(Math.round).join(',')})`;

  // ---------------------------------------------------------------- stage --
  // A stage is drawn in design coordinates (width × height) and scaled as one
  // piece, like a slide: a label placed at (x, y) stays on the same point of
  // the schematic at any window size.
  function setupStage(stage) {
    const W = +stage.getAttribute('width') || 1000;
    const H = +stage.getAttribute('height') || 600;
    const inner = document.createElement('div');
    inner.className = 'cs-stage-inner';
    inner.style.width = W + 'px';
    inner.style.height = H + 'px';
    while (stage.firstChild) inner.appendChild(stage.firstChild);
    stage.appendChild(inner);
    for (const el of inner.children) {
      if (!el.hasAttribute('x')) continue;
      el.classList.add('cs-placed');
      placeChild(el);
    }
    const fit = () => {
      const avail = stage.clientWidth;
      let s = avail / W;
      if (stage.getAttribute('fit') === 'contain') {
        const maxH = stage.parentElement.clientHeight || window.innerHeight;
        s = Math.min(s, maxH / H);
      }
      inner.style.transform = `scale(${s})`;
      stage.style.height = H * s + 'px';
      stage.dataset.scale = s;
    };
    new ResizeObserver(fit).observe(stage);
    fit();
  }

  const ANCHORS = {
    'top-left': [0, 0], top: [-50, 0], 'top-right': [-100, 0], left: [0, -50], center: [-50, -50],
    right: [-100, -50], 'bottom-left': [0, -100], bottom: [-50, -100], 'bottom-right': [-100, -100],
  };

  function placeChild(el) {
    const a = ANCHORS[el.getAttribute('anchor') || 'top-left'] || ANCHORS['top-left'];
    el.style.left = el.getAttribute('x') + 'px';
    el.style.top = el.getAttribute('y') + 'px';
    el.style.transform = `translate(${a[0]}%, ${a[1]}%)`;
  }

  // ------------------------------------------------------------- elements --
  function setupValue(el) {
    const expr = el.getAttribute('expr');
    const name = el.getAttribute('name') || expr;
    const where = `<cs-value ${expr ? 'expr' : 'name'}="${name}">`;
    const get = binding(name, where);
    const from = el.getAttribute('units') || (expr ? '' : unitsOf(name));
    const to = el.getAttribute('unit');
    const scale = el.hasAttribute('scale') ? +el.getAttribute('scale') : unitScale(from, to, name);
    const shown = prettyUnit(to || from);
    const decimals = el.getAttribute('decimals');
    const label = el.getAttribute('label');
    const warnAbove = el.getAttribute('warn-above'), warnBelow = el.getAttribute('warn-below');
    el.innerHTML = (label ? `<span class="cs-label">${label}</span>` : '') +
      '<span class="cs-num"></span>' + (shown && !el.hasAttribute('hide-unit') ? `<span class="cs-unit">${shown}</span>` : '');
    const num = el.querySelector('.cs-num');
    const rec = expr ? null : lookup.get(name.toLowerCase());
    el.title = expr ? `${expr}` : `${name}${rec && rec.input ? ' (input)' : ''}`;
    updaters.push(() => {
      const v = get() * scale;
      num.textContent = format(v, decimals);
      el.classList.toggle('cs-warn', (warnAbove !== null && v > +warnAbove) || (warnBelow !== null && v < +warnBelow));
    });
  }

  function setupInput(el) {
    const name = el.getAttribute('name');
    const rec = lookup.get((name || '').toLowerCase());
    const label = el.getAttribute('label') || name;
    const shown = prettyUnit(el.getAttribute('unit') || (rec && rec.units) || '');
    const v = rec ? rec.value : NaN;
    const min = el.getAttribute('min') ?? (isFinite(v) ? v - Math.abs(v) * 0.5 - 1 : 0);
    const max = el.getAttribute('max') ?? (isFinite(v) ? v + Math.abs(v) * 0.5 + 1 : 1);
    const step = el.getAttribute('step') || 'any';
    el.innerHTML = `<label><span class="cs-label">${label}</span>` +
      `<input type="range" min="${min}" max="${max}" step="${step}">` +
      `<input type="number" step="${step}"><span class="cs-unit">${shown}</span></label>`;
    const [range, box] = el.querySelectorAll('input');
    range.value = box.value = isFinite(v) ? v : '';
    let reason = '';
    if (!rec) reason = `${name} is not a variable of the model`;
    else if (!rec.input) reason = `${name} is not a directly defined input (only "${name} = value" lines can be dashboard inputs)`;
    if (reason) problems.add(`<cs-input name="${name}">: ${reason}`);
    const enable = () => {
      const why = reason || (hosted ? '' : 'Open this dashboard in CoolSolve to change inputs and solve');
      range.disabled = box.disabled = !!why;
      el.title = why;
    };
    const send = (x) => { range.value = box.value = x; post({ type: 'set-input', name, value: +x }); };
    range.addEventListener('input', () => send(range.value));
    box.addEventListener('change', () => send(box.value));
    enable();
    statusListeners.push(enable);
  }

  const ACTIONS = { solve: 'Solve', 'update-guesses': 'Update guesses', stop: 'Stop', 'reset-inputs': 'Reset inputs' };

  function setupButton(el) {
    const action = el.getAttribute('action');
    if (!ACTIONS[action]) problems.add(`<cs-button action="${action}">: unknown action`);
    const b = document.createElement('button');
    b.type = 'button';
    b.textContent = el.textContent.trim() || ACTIONS[action] || action;
    el.textContent = '';
    el.appendChild(b);
    b.addEventListener('click', () => post({ type: 'action', action }));
    const enable = () => {
      b.disabled = !hosted;
      b.title = hosted ? '' : 'Available when the dashboard is opened from the CoolSolve GUI';
    };
    enable();
    statusListeners.push(enable);
  }

  const statusListeners = [];

  function setupStatus(el) {
    const show = () => {
      const s = (data && data.status) || { state: 'unsolved' };
      const when = s.solvedAt ? new Date(s.solvedAt).toLocaleString() : '';
      el.innerHTML = `<span class="cs-pill cs-${s.state}">${s.state}</span> ` +
        `<span class="cs-muted">${hosted ? (s.message || '') : 'snapshot' + (when ? ' of ' + when : '')}</span>`;
    };
    show();
    statusListeners.push(show);
  }

  // [data-cs-*] bindings, on HTML and SVG elements alike
  function setupBindings(root) {
    for (const el of root.querySelectorAll('[data-cs-value]')) {
      const name = el.dataset.csValue;
      const get = binding(name, `data-cs-value="${name}"`);
      const from = unitsOf(name), to = el.dataset.csUnit;
      const scale = el.dataset.csScale ? +el.dataset.csScale : unitScale(from, to, name);
      const unit = 'csHideUnit' in el.dataset ? '' : prettyUnit(to || from);
      updaters.push(() => { el.textContent = format(get() * scale, el.dataset.csDecimals) + (unit ? ' ' + unit : ''); });
    }
    for (const prop of ['fill', 'stroke']) {
      for (const el of root.querySelectorAll(`[data-cs-${prop}]`)) {
        const src = el.getAttribute(`data-cs-${prop}`);
        const get = binding(src, `data-cs-${prop}="${src}"`);
        const stops = parseColormap(el.dataset.csColormap || '0:#3b82f6,100:#ef4444');
        updaters.push(() => { el.style[prop] = colormap(stops, get()); });
      }
    }
    for (const el of root.querySelectorAll('[data-cs-show]')) {
      const get = binding(el.dataset.csShow, `data-cs-show="${el.dataset.csShow}"`);
      updaters.push(() => { el.style.visibility = get() ? '' : 'hidden'; });
    }
    for (const el of root.querySelectorAll('[data-cs-rotate]')) {
      const get = binding(el.dataset.csRotate, `data-cs-rotate="${el.dataset.csRotate}"`);
      const origin = (el.dataset.csOrigin || '0 0').trim().split(/[\s,]+/).join(' ');
      updaters.push(() => { el.setAttribute('transform', `rotate(${get()} ${origin})`); });
    }
    // Moving dashes along a pipe; the period is 1 s when the value equals data-cs-flow-ref.
    for (const el of root.querySelectorAll('[data-cs-flow]')) {
      const get = binding(el.dataset.csFlow, `data-cs-flow="${el.dataset.csFlow}"`);
      const ref = +(el.dataset.csFlowRef || 1);
      el.classList.add('cs-flow');
      updaters.push(() => {
        const v = get();
        el.style.animationDuration = v > 0 ? Math.min(20, ref / v) + 's' : '0s';
        el.style.animationPlayState = v > 0 ? 'running' : 'paused';
      });
    }
  }

  // ---------------------------------------------------------------- plots --
  // Minimal SVG charts (the real runtime uses Plotly, as the GUI): line series,
  // log axes, labelled state points and a cursor that follows the player.
  const PW = 520, PH = 300, M = { l: 58, r: 14, t: 26, b: 42 };
  let plotCount = 0;

  function niceTicks(lo, hi, log) {
    if (log) {
      const out = [];
      for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
        for (const m of (Math.log10(hi / lo) < 2.5 ? [1, 2, 5] : [1])) {
          const v = m * Math.pow(10, e);
          if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
        }
      }
      return out;
    }
    const raw = (hi - lo) / 5, p = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw);
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
    return out;
  }

  function extent(arrs, log) {
    let lo = Infinity, hi = -Infinity;
    for (const a of arrs) for (const v of a) if (isFinite(v) && (!log || v > 0)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    if (log) return [lo / 1.15, hi * 1.15];
    const pad = (hi - lo) * 0.05 || 1;
    return [lo - pad, hi + pad];
  }

  function setupPlot(el) {
    const ns = 'http://www.w3.org/2000/svg';
    const svg = document.createElementNS(ns, 'svg');
    svg.setAttribute('viewBox', `0 0 ${PW} ${PH}`);
    svg.setAttribute('class', 'cs-plot-svg');
    el.appendChild(svg);
    const ref = el.getAttribute('ref');
    const spec = ref ? plotFromRef(ref, el) : plotFromAttributes(el);
    if (!spec) return;
    const logx = el.hasAttribute('log-x') || spec.logx, logy = el.hasAttribute('log-y') || spec.logy;
    const xs = [], ys = [];
    for (const s of spec.series) { xs.push(s.x); ys.push(s.y); }
    // Axis ranges in display units: x-range="lo, hi" / y-range="lo, hi", else the data
    const range = (attr, arrs, log) => {
      const r = (el.getAttribute(attr) || '').split(',').map(Number);
      return r.length === 2 && r.every(isFinite) ? r : extent(arrs, log);
    };
    const [x0, x1] = range('x-range', xs, logx), [y0, y1] = range('y-range', ys, logy);
    const clipId = 'cs-clip-' + (++plotCount);
    const fx = (v) => M.l + (PW - M.l - M.r) * (logx ? Math.log(v / x0) / Math.log(x1 / x0) : (v - x0) / (x1 - x0));
    const fy = (v) => PH - M.b - (PH - M.t - M.b) * (logy ? Math.log(v / y0) / Math.log(y1 / y0) : (v - y0) / (y1 - y0));
    let g = `<text class="cs-plot-title" x="${M.l}" y="16">${el.getAttribute('title') || spec.title || ''}</text>`;
    for (const band of spec.bands) {
      g += `<rect class="cs-plot-band" x="${fx(band[0])}" y="${M.t}" width="${fx(band[1]) - fx(band[0])}" height="${PH - M.t - M.b}"/>` +
        `<text class="cs-plot-bandlabel" x="${fx(band[0]) + 3}" y="${M.t + 11}">${band[2] || ''}</text>`;
    }
    for (const v of niceTicks(x0, x1, logx)) {
      g += `<line class="cs-grid" x1="${fx(v)}" x2="${fx(v)}" y1="${M.t}" y2="${PH - M.b}"/>` +
        `<text class="cs-tick" x="${fx(v)}" y="${PH - M.b + 14}" text-anchor="middle">${format(v)}</text>`;
    }
    for (const v of niceTicks(y0, y1, logy)) {
      g += `<line class="cs-grid" x1="${M.l}" x2="${PW - M.r}" y1="${fy(v)}" y2="${fy(v)}"/>` +
        `<text class="cs-tick" x="${M.l - 5}" y="${fy(v) + 4}" text-anchor="end">${format(v)}</text>`;
    }
    g += `<rect class="cs-frame" x="${M.l}" y="${M.t}" width="${PW - M.l - M.r}" height="${PH - M.t - M.b}"/>` +
      `<text class="cs-axis" x="${(M.l + PW - M.r) / 2}" y="${PH - 6}" text-anchor="middle">${spec.xlabel}</text>` +
      `<text class="cs-axis" transform="rotate(-90)" x="${-(M.t + PH - M.b) / 2}" y="14" text-anchor="middle">${spec.ylabel}</text>` +
      `<clipPath id="${clipId}"><rect x="${M.l}" y="${M.t}" width="${PW - M.l - M.r}" height="${PH - M.t - M.b}"/></clipPath><g clip-path="url(#${clipId})">`;
    for (const s of spec.series) {
      let d = '', pen = false;
      s.x.forEach((xv, j) => {
        const ok = isFinite(xv) && isFinite(s.y[j]) && (!logx || xv > 0) && (!logy || s.y[j] > 0);
        if (ok) { d += (pen ? 'L' : 'M') + fx(xv).toFixed(1) + ' ' + fy(s.y[j]).toFixed(1); }
        pen = ok;
      });
      if (s.close && d) d += 'Z';
      g += `<path class="cs-series ${s.cls || ''}" d="${d}"/>`;
      if (s.labels) {
        s.x.forEach((xv, j) => {
          g += `<circle class="cs-point" cx="${fx(xv)}" cy="${fy(s.y[j])}" r="3.5"/>` +
            `<text class="cs-pointlabel" x="${fx(xv) + 6}" y="${fy(s.y[j]) - 6}">${s.labels[j]}</text>`;
        });
      }
    }
    g += '</g><g class="cs-cursor"><line class="cs-cursor-line"/><circle r="5"/></g>';
    svg.innerHTML = g;
    if (!spec.cursor) { svg.querySelector('.cs-cursor').remove(); return; }
    const cur = svg.querySelector('.cs-cursor'), line = cur.querySelector('line'), dot = cur.querySelector('circle');
    const vertical = spec.cursor.x.toLowerCase() === integralVar;
    updaters.push(() => {
      const xv = value(spec.cursor.x) * spec.xscale, yv = value(spec.cursor.y) * spec.yscale;
      if (!isFinite(xv) || !isFinite(yv)) { cur.style.display = 'none'; return; }
      cur.style.display = '';
      dot.setAttribute('cx', fx(xv)); dot.setAttribute('cy', fy(yv));
      if (vertical) { line.setAttribute('x1', fx(xv)); line.setAttribute('x2', fx(xv)); line.setAttribute('y1', M.t); line.setAttribute('y2', PH - M.b); }
    });
  }

  function axisLabel(name, unit) {
    return name + (prettyUnit(unit) ? ` [${prettyUnit(unit)}]` : '');
  }

  /** A saved plot of the session (here: a thermodynamic diagram with the cycle on it). */
  function plotFromRef(ref, el) {
    const p = data.plots && data.plots[ref];
    if (!p) { problems.add(`<cs-plot ref="${ref}">: no saved plot of that name`); return null; }
    const xu = el.getAttribute('x-unit'), yu = el.getAttribute('y-unit');
    const xscale = unitScale(p.x.units, xu, p.x.name), yscale = unitScale(p.y.units, yu, p.y.name);
    const series = [{ x: p.dome.x.map((v) => v * xscale), y: p.dome.y.map((v) => v * yscale), cls: 'cs-dome' }];
    if (p.path) {
      const pts = arrayPoints(p.path.x, p.path.y);
      series.push({ x: pts.map((q) => q[1] * xscale), y: pts.map((q) => q[2] * yscale), labels: pts.map((q) => q[0]), close: p.path.close, cls: 'cs-cycle' });
    }
    return { title: p.title, series, bands: [], logx: !!p.x.log, logy: !!p.y.log, xlabel: axisLabel(p.x.name, xu || p.x.units), ylabel: axisLabel(p.y.name, yu || p.y.units) };
  }

  /** [index, x, y] of the solved array elements X[i], Y[i]. */
  function arrayPoints(X, Y) {
    const out = [];
    for (let i = 1; lookup.has(`${X}[${i}]`.toLowerCase()); i++) {
      out.push([i, value(`${X}[${i}]`), value(`${Y}[${i}]`)]);
    }
    if (!out.length) problems.add(`path ${X}[i], ${Y}[i]: no such arrays`);
    return out;
  }

  /** A plot of integral columns: x="V" y="p" (y may list several columns). */
  function plotFromAttributes(el) {
    if (!data.integral) { problems.add('<cs-plot>: the model has no integral table'); return null; }
    const xname = el.getAttribute('x') || data.integral.var;
    const ynames = (el.getAttribute('y') || '').split(',').map((s) => s.trim()).filter(Boolean);
    for (const n of [xname, ...ynames]) if (!columns.has(n.toLowerCase())) { problems.add(`<cs-plot>: "${n}" is not a column of the integral table`); return null; }
    const xu = el.getAttribute('x-unit'), yu = el.getAttribute('y-unit');
    const xscale = unitScale(unitsOf(xname), xu, xname), yscale = unitScale(unitsOf(ynames[0]), yu, ynames[0]);
    const series = ynames.map((n, j) => ({ x: columns.get(xname.toLowerCase()).map((v) => v * xscale), y: columns.get(n.toLowerCase()).map((v) => v * yscale), cls: 'cs-s' + j }));
    const bands = [];
    const band = el.getAttribute('band-x');
    if (band) {
      const [a, b] = band.split(';');
      bands.push([binding(a, 'band-x')() * xscale, binding(b, 'band-x')() * xscale, el.getAttribute('band-label') || '']);
    }
    return {
      series, bands, xscale, yscale,
      xlabel: axisLabel(xname, xu || unitsOf(xname)), ylabel: axisLabel(ynames.join(', '), yu || unitsOf(ynames[0])),
      cursor: el.hasAttribute('cursor') ? { x: xname, y: ynames[0] } : null,
    };
  }

  // --------------------------------------------------------------- player --
  // Plays the integral table: the full range of the integration variable lasts
  // `duration` seconds at speed ×1; values between rows are interpolated, so
  // the animation is smooth whatever the number of rows.
  function setupPlayer(el) {
    if (!data.integral) { el.innerHTML = '<span class="cs-muted">No integral table: nothing to animate.</span>'; return; }
    const x = columns.get(integralVar);
    player.t0 = x[0]; player.t1 = x[x.length - 1]; player.t = player.t0;
    player.duration = +(el.getAttribute('duration') || 10);
    player.loop = !el.hasAttribute('no-loop');
    const unit = prettyUnit(unitsOf(integralVar));
    el.innerHTML =
      '<button type="button" data-a="start" title="Start (Home)">⏮</button>' +
      '<button type="button" data-a="back" title="Previous row (←)">◀︎◀︎</button>' +
      '<button type="button" data-a="play" class="cs-play" title="Play / pause (space)">▶</button>' +
      '<button type="button" data-a="fwd" title="Next row (→)">▶︎▶︎</button>' +
      '<input type="range" min="0" max="1000" step="1" value="0" aria-label="position">' +
      `<span class="cs-player-time"></span>` +
      '<select title="Speed"><option>0.25</option><option>0.5</option><option selected>1</option><option>2</option><option>4</option></select>' +
      `<label title="Loop"><input type="checkbox" ${player.loop ? 'checked' : ''}> loop</label>`;
    const play = el.querySelector('.cs-play'), range = el.querySelector('input[type=range]');
    const time = el.querySelector('.cs-player-time'), speed = el.querySelector('select');
    const loop = el.querySelector('input[type=checkbox]');
    // One row forward or back: the next row after the position, or the last one before it
    const step = (dir) => {
      const eps = 1e-9 * (player.t1 - player.t0);
      let k = 0;
      while (k + 1 < x.length && x[k + 1] <= player.t + eps) k++;
      if (dir > 0) k = Math.min(x.length - 1, k + 1);
      else if (Math.abs(x[k] - player.t) <= eps) k = Math.max(0, k - 1);
      seek(x[k]);
    };
    let last = null;
    const toggle = () => { player.playing = !player.playing; if (player.playing && player.t >= player.t1) player.t = player.t0; last = null; };
    el.addEventListener('click', (e) => {
      const a = e.target.closest('button') && e.target.closest('button').dataset.a;
      if (a === 'play') toggle();
      if (a === 'start') seek(player.t0);
      if (a === 'back') { player.playing = false; step(-1); }
      if (a === 'fwd') { player.playing = false; step(1); }
    });
    range.addEventListener('input', () => { player.playing = false; seek(player.t0 + (player.t1 - player.t0) * range.value / 1000); });
    speed.addEventListener('change', () => { player.speed = +speed.value; });
    loop.addEventListener('change', () => { player.loop = loop.checked; });
    document.addEventListener('keydown', (e) => {
      if (/^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName)) return;
      if (e.key === ' ') { e.preventDefault(); toggle(); }
      if (e.key === 'ArrowRight') { player.playing = false; step(1); }
      if (e.key === 'ArrowLeft') { player.playing = false; step(-1); }
      if (e.key === 'Home') seek(player.t0);
    });
    updaters.push(() => {
      play.textContent = player.playing ? '❚❚' : '▶';
      range.value = Math.round(1000 * (player.t - player.t0) / (player.t1 - player.t0));
      time.textContent = `${data.integral.var} = ${format(player.t, 1)}${unit ? ' ' + unit : ''}`;
    });
    const tick = (now) => {
      if (player.playing) {
        if (last !== null) {
          player.t += (player.t1 - player.t0) * (now - last) / 1000 / player.duration * player.speed;
          if (player.t >= player.t1) {
            if (player.loop) player.t = player.t0 + (player.t - player.t1); else { player.t = player.t1; player.playing = false; }
          }
        }
        last = now;
        render();
      } else last = null;
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
    if (el.hasAttribute('autoplay')) player.playing = true;
  }

  function seek(t) {
    player.t = Math.min(player.t1, Math.max(player.t0, t));
    render();
  }

  // ------------------------------------------------------------ rendering --
  const ctx = { get: value, format, get t() { return player.t; }, get playing() { return player.playing; }, get data() { return data; } };

  function render() {
    for (const u of updaters) u();
    for (const f of frameHooks) {
      try { f(ctx); } catch (e) { problems.add(`onFrame: ${e.message}`); showProblems(); }
    }
  }

  /**
   * Registers a drawing function, called after the data arrive and at every
   * frame of the player, with ctx.get(name), ctx.t, ctx.playing, ctx.format.
   */
  CS.onFrame = (fn) => { frameHooks.push(fn); if (data) fn(ctx); };
  CS.value = value;
  CS.format = format;

  // ------------------------------------------------------- problem badge --
  // Unknown variables and other binding errors are listed in a badge, so that
  // a dashboard written for another version of the model shows what to fix.
  function showProblems() {
    let badge = document.querySelector('.cs-problems');
    if (!problems.size) { if (badge) badge.remove(); return; }
    if (!badge) {
      badge = document.createElement('details');
      badge.className = 'cs-problems';
      document.body.appendChild(badge);
    }
    badge.innerHTML = `<summary>${problems.size} binding problem${problems.size > 1 ? 's' : ''}</summary><ul>` +
      [...problems].map((p) => `<li>${p.replace(/</g, '&lt;')}</li>`).join('') + '</ul>';
  }

  // ------------------------------------------------------- layout helper --
  // Press L (or open the page with #layout): the stage shows the design
  // coordinates under the pointer, and placed elements can be dragged; the new
  // x/y attributes are shown and copied, ready to paste into the HTML file.
  function setupLayoutHelper() {
    let on = false;
    const tip = document.createElement('div');
    tip.className = 'cs-layout-tip';
    document.body.appendChild(tip);
    const toggle = () => { on = !on; document.body.classList.toggle('cs-layout', on); tip.style.display = on ? 'block' : 'none'; };
    document.addEventListener('keydown', (e) => {
      if (/^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName) || e.ctrlKey || e.metaKey) return;
      if (e.key === 'l' || e.key === 'L') toggle();
    });
    if (location.hash === '#layout') toggle();
    const stageXY = (inner, e) => {
      const r = inner.getBoundingClientRect(), s = r.width / inner.offsetWidth;
      return [Math.round((e.clientX - r.left) / s), Math.round((e.clientY - r.top) / s)];
    };
    let drag = null;
    document.addEventListener('pointermove', (e) => {
      if (!on) return;
      const inner = (drag && drag.inner) || (e.target.closest && e.target.closest('.cs-stage-inner'));
      if (!inner) { tip.textContent = 'Layout mode (L): point at a stage'; return; }
      const [x, y] = stageXY(inner, e);
      tip.textContent = drag ? `x="${x - drag.dx}" y="${y - drag.dy}"` : `x=${x}  y=${y}`;
      if (drag) { drag.el.setAttribute('x', x - drag.dx); drag.el.setAttribute('y', y - drag.dy); placeChild(drag.el); }
    });
    document.addEventListener('pointerdown', (e) => {
      if (!on) return;
      const el = e.target.closest && e.target.closest('.cs-placed');
      if (!el) return;
      e.preventDefault();
      const inner = el.closest('.cs-stage-inner');
      const [x, y] = stageXY(inner, e);
      drag = { el, inner, dx: x - +el.getAttribute('x'), dy: y - +el.getAttribute('y') };
    });
    document.addEventListener('pointerup', () => {
      if (!drag) return;
      const el = drag.el;
      drag = null;
      // The start tag as written by the author: without what the runtime added
      const copy = el.cloneNode(false);
      copy.classList.remove('cs-placed', 'cs-warn');
      if (!copy.classList.length) copy.removeAttribute('class');
      copy.removeAttribute('style');
      copy.removeAttribute('title');
      const tag = copy.outerHTML.replace(/>[\s\S]*$/, '>');
      tip.textContent = `${tag} — copied`;
      if (navigator.clipboard) navigator.clipboard.writeText(tag).catch(() => {});
      console.log('[CoolSolve layout]', tag);
    });
  }

  // ---------------------------------------------------------------- host --
  // Protocol with a CoolSolve host page (planned GUI side): the dashboard posts
  // {source: 'coolsolve-dashboard', type: 'ready' | 'action' | 'set-input'} to
  // its parent; the host answers with {source: 'coolsolve-host', type: 'data',
  // data, theme}. Without an answer the dashboard stays a standalone snapshot.
  function post(msg) {
    if (hosted) window.parent.postMessage(Object.assign({ source: 'coolsolve-dashboard', version: 1 }, msg), '*');
  }

  window.addEventListener('message', (e) => {
    const m = e.data;
    if (!m || m.source !== 'coolsolve-host' || e.source !== window.parent) return;
    if (m.theme) document.documentElement.dataset.theme = m.theme;
    if (m.type === 'data') { hosted = true; start(m.data); }
  });

  // ---------------------------------------------------------------- start --
  let built = false;

  function start(d) {
    data = d;
    lookup = new Map(Object.entries(data.variables || {}).map(([k, v]) => [k.toLowerCase(), v]));
    columns = new Map(data.integral ? Object.entries(data.integral.columns).map(([k, v]) => [k.toLowerCase(), v]) : []);
    integralVar = data.integral ? data.integral.var.toLowerCase() : null;
    if (!built) {
      built = true;
      for (const s of document.querySelectorAll('cs-stage')) setupStage(s);
      for (const el of document.querySelectorAll('cs-player')) setupPlayer(el);
      for (const el of document.querySelectorAll('cs-value')) setupValue(el);
      for (const el of document.querySelectorAll('cs-input')) setupInput(el);
      for (const el of document.querySelectorAll('cs-button')) setupButton(el);
      for (const el of document.querySelectorAll('cs-status')) setupStatus(el);
      for (const el of document.querySelectorAll('cs-plot')) setupPlot(el);
      setupBindings(document);
      setupLayoutHelper();
      for (const el of document.querySelectorAll('[data-cs-model-title]')) el.textContent = (data.model && data.model.title) || '';
    }
    for (const f of statusListeners) f();
    render();
    showProblems();
  }

  function loadStandalone() {
    if (window.COOLSOLVE_DATA) return start(window.COOLSOLVE_DATA);
    const inline = document.getElementById('coolsolve-data');
    if (inline) return start(JSON.parse(inline.textContent));
    const root = document.querySelector('cs-dashboard');
    const s = document.createElement('script');
    s.src = (root && root.getAttribute('results')) || 'results.js';
    s.onload = () => (window.COOLSOLVE_DATA ? start(window.COOLSOLVE_DATA) : null);
    s.onerror = () => { problems.add(`no results: ${s.src} not found (open the dashboard from CoolSolve, or export it as a standalone file)`); showProblems(); };
    document.head.appendChild(s);
  }

  function boot() {
    injectStyle();
    if (window.parent !== window) {
      window.parent.postMessage({ source: 'coolsolve-dashboard', version: 1, type: 'ready' }, '*');
      setTimeout(() => { if (!hosted) loadStandalone(); }, 400);
    } else loadStandalone();
  }

  // ---------------------------------------------------------------- style --
  // Defaults only: a dashboard restyles everything with its own CSS, through
  // the --cs-* variables or the class names.
  function injectStyle() {
    const css = `
:root{--cs-bg:#f8fafc;--cs-card:#ffffff;--cs-fg:#0f172a;--cs-muted:#64748b;--cs-line:#e2e8f0;--cs-accent:#2563eb;--cs-warn:#dc2626;--cs-ok:#16a34a;--cs-s0:#2563eb;--cs-s1:#dc2626;--cs-s2:#16a34a;--cs-cycle:#ea580c}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--cs-bg:#0b1220;--cs-card:#111a2e;--cs-fg:#e2e8f0;--cs-muted:#94a3b8;--cs-line:#24314d;--cs-accent:#60a5fa;--cs-warn:#f87171;--cs-ok:#4ade80;--cs-s0:#60a5fa;--cs-s1:#f87171;--cs-s2:#4ade80;--cs-cycle:#fb923c}}
:root[data-theme=dark]{--cs-bg:#0b1220;--cs-card:#111a2e;--cs-fg:#e2e8f0;--cs-muted:#94a3b8;--cs-line:#24314d;--cs-accent:#60a5fa;--cs-warn:#f87171;--cs-ok:#4ade80;--cs-s0:#60a5fa;--cs-s1:#f87171;--cs-s2:#4ade80;--cs-cycle:#fb923c}
cs-stage{display:block;position:relative;overflow:hidden;width:100%}
.cs-stage-inner{position:absolute;left:0;top:0;transform-origin:0 0}
.cs-stage-inner>svg{position:absolute;left:0;top:0}
.cs-placed{position:absolute;white-space:nowrap}
cs-value{display:inline-flex;align-items:baseline;gap:.3em;font-variant-numeric:tabular-nums}
cs-value .cs-label{color:var(--cs-muted)}
cs-value .cs-unit,cs-input .cs-unit{color:var(--cs-muted);font-size:.85em}
cs-value.cs-warn .cs-num{color:var(--cs-warn);font-weight:600}
cs-input{display:block}
cs-input label{display:grid;grid-template-columns:minmax(5em,auto) minmax(0,1fr) 5em 2.5em;gap:.5em;align-items:center}
cs-input input[type=range]{min-width:0;width:100%}
cs-input input[type=number]{width:100%;font:inherit;padding:.15em .3em;border:1px solid var(--cs-line);border-radius:4px;background:var(--cs-card);color:var(--cs-fg)}
cs-input input:disabled{opacity:.6;cursor:not-allowed}
cs-button button,cs-player button,cs-player select{font:inherit;padding:.35em .8em;border-radius:6px;border:1px solid var(--cs-line);background:var(--cs-card);color:var(--cs-fg);cursor:pointer}
cs-button button:not(:disabled){background:var(--cs-accent);border-color:var(--cs-accent);color:#fff}
cs-button button:disabled{opacity:.55;cursor:not-allowed}
.cs-pill{display:inline-block;padding:.05em .55em;border-radius:999px;font-size:.8em;font-weight:600;text-transform:uppercase;letter-spacing:.04em;background:var(--cs-line)}
.cs-pill.cs-solved{background:color-mix(in srgb,var(--cs-ok) 20%,transparent);color:var(--cs-ok)}
.cs-pill.cs-failed{background:color-mix(in srgb,var(--cs-warn) 20%,transparent);color:var(--cs-warn)}
.cs-muted{color:var(--cs-muted)}
cs-plot{display:block}
.cs-plot-svg{width:100%;height:auto;display:block;font:11px system-ui,sans-serif}
.cs-plot-title{font-weight:600;fill:var(--cs-fg);font-size:12px}
.cs-tick{fill:var(--cs-muted)} .cs-axis{fill:var(--cs-fg)}
.cs-grid{stroke:var(--cs-line);stroke-width:1} .cs-frame{fill:none;stroke:var(--cs-muted);stroke-width:1}
.cs-series{fill:none;stroke:var(--cs-s0);stroke-width:2;stroke-linejoin:round}
.cs-series.cs-s1{stroke:var(--cs-s1)} .cs-series.cs-s2{stroke:var(--cs-s2)}
.cs-series.cs-dome{stroke:var(--cs-muted);stroke-width:1.5} .cs-series.cs-cycle{stroke:var(--cs-cycle);stroke-width:2.5}
.cs-point{fill:var(--cs-cycle)} .cs-pointlabel{fill:var(--cs-fg);font-weight:600}
.cs-plot-band{fill:var(--cs-cycle);opacity:.12} .cs-plot-bandlabel{fill:var(--cs-cycle);font-size:10px}
.cs-cursor circle{fill:var(--cs-accent);stroke:var(--cs-card);stroke-width:2}
.cs-cursor-line{stroke:var(--cs-accent);stroke-width:1;stroke-dasharray:3 3}
cs-player{display:flex;align-items:center;gap:.5em;flex-wrap:wrap}
cs-player input[type=range]{flex:1;min-width:10em}
.cs-player-time{font-variant-numeric:tabular-nums;min-width:9em}
.cs-flow{stroke-dasharray:7 9;animation:cs-flow 1s linear infinite}
@keyframes cs-flow{to{stroke-dashoffset:-16}}
@media (prefers-reduced-motion:reduce){.cs-flow{animation:none}}
.cs-problems{position:fixed;right:12px;bottom:12px;max-width:min(32em,90vw);background:#fef2f2;color:#7f1d1d;border:1px solid #fca5a5;border-radius:8px;padding:.4em .8em;font:13px system-ui,sans-serif;z-index:10}
.cs-problems summary{cursor:pointer;font-weight:600}
.cs-layout .cs-placed{outline:1px dashed var(--cs-accent);cursor:move}
.cs-layout cs-stage{outline:2px solid var(--cs-accent)}
.cs-layout-tip{display:none;position:fixed;left:12px;bottom:12px;background:var(--cs-fg);color:var(--cs-bg);padding:.3em .7em;border-radius:6px;font:12px ui-monospace,monospace;z-index:10}
`;
    const style = document.createElement('style');
    style.textContent = css;
    document.head.prepend(style);
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot);
  else boot();
})();
