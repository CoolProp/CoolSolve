/**
 * DiagramControls — the chart settings of a diagram definition: fluid (and the concentration of a
 * solution), diagram type, units, region (or the temperature and pressure range of a fluid without
 * saturation curve; for humid air the total pressure and the temperature and humidity ratio range), theme,
 * label mode, line families (on/off, count, rounding), critical point and saturation scale. What can be
 * chosen comes from `GET /api/v1/diagram/options`, per kind of fluid.
 */

import { useEffect, useState } from 'react';
import type { DiagramDefinition, DiagramOptionsResponse, FluidInfo } from '../../api/types';
import {
  DEFAULT_REGION,
  defaultFractions,
  diagramOption,
  diagramTypesOf,
  familiesOf,
  fitToKind,
  fluidKind,
  hasSaturationScale,
  humidAirRegion,
  isHumidAirKind,
  lineWeight,
  regionBounds,
  sameFluid,
  sameJson,
  unitsName,
  unitsRecord,
  withDiagramType,
  withFamily,
  withLineWeight,
  withView,
} from './definition';
import { formatValue, unitTableOf } from './options';

interface Props {
  def: DiagramDefinition;
  options: DiagramOptionsResponse;
  /** Fluids the server draws */
  fluids: FluidInfo[];
  /** Fluids of the model, offered first */
  modelFluids: string[];
  /** The concentration of the model's solutions (mass fractions), the default of their charts */
  modelFractions?: Record<string, number[]>;
  /** The total pressure of the model's humid air (Pa), the default of its charts */
  modelPressure?: number;
  onChange: (def: DiagramDefinition) => void;
}

/**
 * A number typed in a field, committed when valid (inside [min, max] when given): as it is typed, or when
 * `lazy` once the field is left (or Enter pressed); the field then shows the value committed. An empty
 * field commits undefined when `optional`.
 */
function NumberField({ value, onCommit, min, max, step, optional, lazy, placeholder, title, className }: {
  value: number | undefined;
  onCommit: (v: number | undefined) => void;
  min?: number;
  max?: number;
  step?: number;
  optional?: boolean;
  lazy?: boolean;
  placeholder?: string;
  title?: string;
  className?: string;
}) {
  const shown = value === undefined ? '' : formatValue(value);
  const [text, setText] = useState(shown);
  useEffect(() => setText(shown), [shown]);
  const commit = (t: string) => {
    if (t.trim() === '') {
      if (optional && value !== undefined) onCommit(undefined);
      return;
    }
    const v = Number(t);
    if (Number.isFinite(v) && v !== value && (min === undefined || v >= min) && (max === undefined || v <= max)) onCommit(v);
  };
  return (
    <input
      type="number"
      className={className ?? 'diagram-num'}
      value={text}
      min={min}
      max={max}
      step={step}
      placeholder={placeholder}
      title={title}
      onChange={(e) => {
        setText(e.target.value);
        if (!lazy) commit(e.target.value);
      }}
      onBlur={() => {
        if (lazy) commit(text);
        setText(shown);
      }}
      onKeyDown={(e) => { if (e.key === 'Enter') e.currentTarget.blur(); }}
    />
  );
}

/** The bounds of a region record, in the order they are shown. */
const BOUNDS = ['T_min', 'T_max', 'p_min', 'p_max'] as const;

/** How a quantity key (or a family's symbol) is written in the controls. */
const SYMBOLS: Record<string, string> = { rho: 'ρ', Q: 'x', RH: 'φ', phi: 'φ' };
const symbol = (key: string) => SYMBOLS[key] ?? key;

const UNIT_LABELS: Record<string, string> = { coolsolve: 'model (°C, Pa, J/kg)' };
const THEME_AUTO = '';
const MAX_LINES = 60;

export default function DiagramControls({ def, options, fluids, modelFluids, modelFractions, modelPressure, onChange }: Props) {
  const spec = def.spec;
  const type = diagramOption(options, spec.diagram);
  const setSpec = (next: typeof spec) => onChange({ ...def, spec: next });

  // What the kind of the fluid offers (a pure fluid, a mixture, an incompressible solution, an ideal gas)
  const fluidInfo = fluids.find((f) => sameFluid(f.name, def.fluid));
  const kind = fluidKind(options, fluidInfo?.type);
  const saturation = kind?.saturation ?? true;
  const diagramTypes = diagramTypesOf(options, kind);
  const fractions = def.fractions ?? defaultFractions(fluidInfo, modelFractions);
  const fractionRange = fluidInfo?.fractionRange;
  const humidAir = isHumidAirKind(kind);

  const setFluid = (name: string) => {
    const info = fluids.find((f) => f.name === name);
    const next = fluidKind(options, info?.type);
    const unzoomed = withView(spec, null);
    const { fractions: _fractions, p_Pa: _pressure, ...rest } = def;
    const given = defaultFractions(info, modelFractions);
    // Humid air at the model's total pressure, else at the kind's default (the server's when absent)
    const pressure = isHumidAirKind(next) ? modelPressure : undefined;
    onChange({
      ...rest,
      fluid: name,
      ...(given ? { fractions: given } : {}),
      ...(pressure !== undefined ? { p_Pa: pressure } : {}),
      spec: sameJson(next, kind) ? unzoomed : fitToKind(unzoomed, next, options),
    });
  };

  // The temperature and pressure range of a fluid without saturation curve: absolute bounds in the chart's
  // units, the others those of the preset DEF (the fluid's default box)
  const units = unitTableOf(options, spec.units);
  const presetLimits = options.regions.find((r) => r.name === DEFAULT_REGION)?.limits ?? {};
  const limits: Record<string, unknown> = typeof spec.limits === 'string' ? presetLimits : spec.limits;
  const boundShown = (name: string): number | undefined => {
    const v = limits[name];
    const unit = units[name[0] === 'T' ? 'T' : 'p'];
    return typeof v === 'number' && unit ? v * unit.scale + unit.offset : undefined;
  };
  const setBound = (name: string, shown: number | undefined) => {
    const unit = units[name[0] === 'T' ? 'T' : 'p'];
    if (!unit) return;
    const next: Record<string, unknown> = {};
    for (const b of BOUNDS) next[b] = typeof limits[b] === 'number' ? limits[b] : presetLimits[b] ?? null;
    next[name] = shown === undefined ? presetLimits[name] ?? null : (shown - unit.offset) / unit.scale;
    const asPreset = BOUNDS.every((b) => sameJson(next[b], presetLimits[b] ?? null));
    setSpec({ ...withView(spec, null), limits: asPreset ? DEFAULT_REGION : next });
  };

  // Humid air: the total pressure and the region of the chart (bounds on T, or Thx on the Mollier chart, and W) in the
  // chart's units
  const temperatureKey = type?.y === 'Thx' ? 'Thx' : 'T';
  const airBounds = regionBounds(spec.limits) ?? regionBounds(humidAirRegion(kind?.default_region, type ?? { x: 'T', y: 'W' })) ?? {};
  const pressure = def.p_Pa ?? modelPressure ?? kind?.total_pressure ?? 101325;
  const pressureRange = kind?.total_pressure_range;
  const inUnits = (key: string, value: number | undefined) => {
    const unit = units[key];
    return value !== undefined && unit ? value * unit.scale + unit.offset : value;
  };
  const toSI = (key: string, value: number) => {
    const unit = units[key];
    return unit ? (value - unit.offset) / unit.scale : value;
  };
  const setAirBound = (key: string, end: 0 | 1, value: number | undefined) => {
    if (value === undefined || !type) return;
    const bounds = { ...airBounds };
    const pair: [number, number] = [...(bounds[key] ?? [0, 0])] as [number, number];
    pair[end] = toSI(key, value);
    if (!(pair[1] > pair[0])) return;  // the region needs two different values, low then high
    bounds[key] = pair;
    setSpec({ ...withView(spec, null), limits: humidAirRegion({ bounds }, type) });
  };

  const inModel = (name: string) => modelFluids.some((m) => sameFluid(m, name));
  const modelList = fluids.filter((f) => inModel(f.name));
  const otherList = fluids.filter((f) => !inModel(f.name));
  const fluidListed = fluids.some((f) => f.name === def.fluid);

  const unitNames = options.unit_systems.map((u) => u.name);
  const regionNames = options.regions.map((r) => r.name);
  const currentUnits = unitsName(spec.units);
  const regionName = typeof spec.limits === 'string' ? spec.limits : '';
  const look = spec.settings?.look ?? THEME_AUTO;
  const weight = lineWeight(spec);
  const labelMode = spec.settings?.labels?.mode ?? 'rails';
  const marksCritical = spec.marked_points?.includes('critical') ?? false;
  const scaleAllowed = type ? hasSaturationScale(type) : false;

  const setUnits = (name: string) => {
    const units = unitsRecord(options, name);
    if (units !== undefined) setSpec({ ...spec, units });
  };

  const setSettings = (change: Record<string, unknown>) =>
    setSpec({ ...spec, settings: { ...spec.settings, ...change } });

  const setLook = (value: string) => {
    const { look: _look, ...rest } = spec.settings ?? {};
    setSpec({ ...spec, settings: value === THEME_AUTO ? rest : { ...rest, look: value } });
  };

  const setCritical = (on: boolean) => {
    const others = (spec.marked_points ?? []).filter((p) => p !== 'critical');
    const marked = on ? ['critical', ...others] : others;
    const { marked_points: _marked, ...rest } = spec;
    setSpec(marked.length > 0 ? { ...rest, marked_points: marked } : rest);
  };

  const setScale = (on: boolean) => {
    const { saturation_scale: _scale, ...rest } = spec;
    setSpec(on ? { ...rest, saturation_scale: { Q: 1 } } : rest);
  };

  return (
    <div className="diagram-controls">
      <div className="diagram-row">
        <label className="diagram-field" title="Fluid of the chart">
          <span>Fluid</span>
          <select
            className="diagram-select"
            value={def.fluid}
            onChange={(e) => setFluid(e.target.value)}
          >
            {!fluidListed && <option value={def.fluid}>{def.fluid}</option>}
            {modelList.length > 0 && (
              <optgroup label="Model fluids">
                {modelList.map((f) => <option key={f.name} value={f.name}>{f.name}</option>)}
              </optgroup>
            )}
            <optgroup label="All fluids">
              {otherList.map((f) => <option key={f.name} value={f.name}>{f.name}</option>)}
            </optgroup>
          </select>
        </label>

        {kind?.fractions && fractions && (
          <label
            className="diagram-field diagram-fraction"
            title={`Concentration of the solution: mass fraction${fractionRange ? ` between ${fractionRange[0]} and ${fractionRange[1]}` : ''}`}
          >
            <span>Mass fraction</span>
            <NumberField
              value={fractions[0]}
              min={fractionRange?.[0]}
              max={fractionRange?.[1]}
              step={0.01}
              onCommit={(v) => { if (v !== undefined) onChange({ ...def, fractions: [v] }); }}
            />
          </label>
        )}

        <label className="diagram-field" title="Diagram type">
          <span>Diagram</span>
          <select
            className="diagram-select"
            value={spec.diagram}
            onChange={(e) => {
              const to = diagramOption(options, e.target.value);
              if (to) setSpec(withDiagramType(spec, type, to, kind));
            }}
          >
            {!diagramTypes.some((d) => d.name === spec.diagram) && <option value={spec.diagram}>{type?.title || spec.diagram}</option>}
            {diagramTypes.map((d) => (
              <option key={d.name} value={d.name}>{d.title || d.name}</option>
            ))}
          </select>
        </label>

        <label className="diagram-field" title="Units of the axes and labels">
          <span>Units</span>
          <select className="diagram-select" value={currentUnits} onChange={(e) => setUnits(e.target.value)}>
            {!unitNames.includes(currentUnits) && <option value={currentUnits}>{currentUnits}</option>}
            {unitNames.map((name) => (
              <option key={name} value={name}>{UNIT_LABELS[name] ?? name}</option>
            ))}
          </select>
        </label>

        {humidAir ? (
          <span className="diagram-field diagram-range" title="Total pressure of the humid air, and the region of the chart">
            <span className="diagram-field" title="Total pressure of the humid air (the model's when it has one)">
              <span>p [{units.p?.label ?? 'Pa'}]</span>
              <NumberField
                className="diagram-num diagram-pressure"
                value={inUnits('p', pressure)}
                lazy
                min={pressureRange ? inUnits('p', pressureRange[0]) : undefined}
                max={pressureRange ? inUnits('p', pressureRange[1]) : undefined}
                title="Total pressure of the humid air"
                onCommit={(v) => { if (v !== undefined) onChange({ ...def, p_Pa: toSI('p', v) }); }}
              />
            </span>
            {[temperatureKey, 'W'].map((key) => (
              <span key={key} className="diagram-field">
                <span>{key === 'W' ? (type?.y === 'Thx' ? 'x' : 'W') : 'T'} [{units[key]?.label ?? (key === 'W' ? 'kg/kg' : 'K')}]</span>
                {([0, 1] as const).map((end) => (
                  <NumberField
                    key={end}
                    className="diagram-num diagram-bound"
                    value={inUnits(key, airBounds[key]?.[end])}
                    lazy
                    placeholder={end === 0 ? 'min' : 'max'}
                    title={`${end === 0 ? 'Lowest' : 'Highest'} ${key === 'W' ? 'humidity ratio' : 'temperature'} of the chart`}
                    onCommit={(v) => setAirBound(key, end, v)}
                  />
                ))}
              </span>
            ))}
          </span>
        ) : saturation ? (
          <label className="diagram-field" title="Temperature and pressure range of the lines (CoolPlot region preset)">
            <span>Region</span>
            <select
              className="diagram-select"
              value={regionName}
              onChange={(e) => setSpec({ ...withView(spec, null), limits: e.target.value })}
            >
              {regionName === '' && <option value="" disabled>custom</option>}
              {regionName !== '' && !regionNames.includes(regionName) && <option value={regionName}>{regionName}</option>}
              {regionNames.map((r) => <option key={r} value={r}>{r}</option>)}
            </select>
          </label>
        ) : (
          <span
            className="diagram-field diagram-range"
            title="Temperature and pressure range of the lines; an empty field takes the default of the fluid"
          >
            {(['T', 'p'] as const).map((key) => (
              <span key={key} className="diagram-field">
                <span>{key} [{units[key]?.label ?? (key === 'T' ? 'K' : 'Pa')}]</span>
                {(['min', 'max'] as const).map((end) => (
                  <NumberField
                    key={end}
                    className="diagram-num diagram-bound"
                    value={boundShown(`${key}_${end}`)}
                    optional
                    lazy
                    placeholder={end === 'min' ? 'min' : 'max'}
                    title={`${end === 'min' ? 'Lowest' : 'Highest'} ${key === 'T' ? 'temperature' : 'pressure'} of the lines`}
                    onCommit={(v) => setBound(`${key}_${end}`, v)}
                  />
                ))}
              </span>
            ))}
          </span>
        )}

        <label className="diagram-field" title="Colours and fonts of the chart">
          <span>Theme</span>
          <select className="diagram-select" value={look} onChange={(e) => setLook(e.target.value)}>
            <option value={THEME_AUTO}>as the GUI</option>
            {look !== THEME_AUTO && !options.themes.includes(look) && <option value={look}>{look}</option>}
            {options.themes.map((t) => <option key={t} value={t}>{t}</option>)}
          </select>
        </label>

        <label className="diagram-field" title="Light: thin, partly transparent isolines, so that the saturation curve and the cycles stand out. As the theme: CoolPlot's line widths">
          <span>Line weight</span>
          <select
            className="diagram-select"
            value={weight}
            onChange={(e) => setSpec(withLineWeight(spec, e.target.value as 'light' | 'theme', options))}
          >
            <option value="light">light</option>
            <option value="theme">as the theme</option>
            {weight === 'custom' && <option value="custom" disabled>custom</option>}
          </select>
        </label>

        <label className="diagram-field" title="Where the lines are labelled">
          <span>Labels</span>
          <select
            className="diagram-select"
            value={labelMode}
            onChange={(e) => setSettings({ labels: { ...spec.settings?.labels, mode: e.target.value, symbols: spec.settings?.labels?.symbols ?? 'unicode' } })}
          >
            {!options.label_modes.includes(labelMode) && <option value={labelMode}>{labelMode}</option>}
            {options.label_modes.map((m) => <option key={m} value={m}>{m}</option>)}
          </select>
        </label>
      </div>

      <div className="diagram-row">
        <span className="diagram-row-title">Lines</span>
        {familiesOf(type, kind).map((family) => {
          const name = symbol(family.symbol || family.key);
          const set = spec.isolines.find((s) => s.key === family.key);
          return (
            <span key={family.key} className={`diagram-family${set ? ' on' : ''}`}>
              <label className="diagram-field" title={`Lines of constant ${(family.label || name).toLowerCase()}`}>
                <input
                  type="checkbox"
                  checked={!!set}
                  onChange={(e) => setSpec(withFamily(spec, family.key, { on: e.target.checked }, family.default_num, !saturation))}
                />
                {name}
              </label>
              {set && (
                <>
                  <input
                    type="number"
                    className="diagram-num"
                    min={1}
                    max={MAX_LINES}
                    value={set.num}
                    title="Number of lines"
                    onChange={(e) => {
                      const num = parseInt(e.target.value, 10);
                      if (num >= 1 && num <= MAX_LINES) setSpec(withFamily(spec, family.key, { num }, family.default_num));
                    }}
                  />
                  <label className="diagram-field" title="Round the values of the lines">
                    <input
                      type="checkbox"
                      checked={!!set.rounding}
                      onChange={(e) => setSpec(withFamily(spec, family.key, { rounding: e.target.checked }, family.default_num))}
                    />
                    round
                  </label>
                </>
              )}
            </span>
          );
        })}
        {saturation && (
          <>
            <label className="diagram-field" title="Mark the critical point">
              <input type="checkbox" checked={marksCritical} onChange={(e) => setCritical(e.target.checked)} />
              critical point
            </label>
            <label
              className="diagram-field"
              title={scaleAllowed
                ? 'Saturation temperatures (or pressures) next to the vertical axis'
                : 'Only with pressure or temperature on the vertical axis and neither on the horizontal one'}
            >
              <input
                type="checkbox"
                checked={!!spec.saturation_scale && scaleAllowed}
                disabled={!scaleAllowed}
                onChange={(e) => setScale(e.target.checked)}
              />
              saturation scale
            </label>
          </>
        )}
      </div>
    </div>
  );
}
