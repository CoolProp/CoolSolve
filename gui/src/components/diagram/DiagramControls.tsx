/**
 * DiagramControls — the chart settings of a diagram definition: fluid, diagram type, units, region,
 * theme, label mode, line families (on/off, count, rounding), critical point and saturation scale.
 * What can be chosen comes from `GET /api/v1/diagram/options`.
 */

import type { DiagramDefinition, DiagramOptionsResponse, FluidInfo } from '../../api/types';
import {
  diagramOption,
  hasSaturationScale,
  sameFluid,
  unitsName,
  unitsRecord,
  withDiagramType,
  withFamily,
  withView,
} from './definition';

interface Props {
  def: DiagramDefinition;
  options: DiagramOptionsResponse;
  /** Fluids the server draws */
  fluids: FluidInfo[];
  /** Fluids of the model, offered first */
  modelFluids: string[];
  onChange: (def: DiagramDefinition) => void;
}

/** How a quantity key is written in the controls. */
const SYMBOLS: Record<string, string> = { rho: 'ρ', Q: 'x' };
const symbol = (key: string) => SYMBOLS[key] ?? key;

const UNIT_LABELS: Record<string, string> = { coolsolve: 'model (°C, Pa, J/kg)' };
const THEME_AUTO = '';
const MAX_LINES = 60;

export default function DiagramControls({ def, options, fluids, modelFluids, onChange }: Props) {
  const spec = def.spec;
  const type = diagramOption(options, spec.diagram);
  const setSpec = (next: typeof spec) => onChange({ ...def, spec: next });

  const inModel = (name: string) => modelFluids.some((m) => sameFluid(m, name));
  const modelList = fluids.filter((f) => inModel(f.name));
  const otherList = fluids.filter((f) => !inModel(f.name));
  const fluidListed = fluids.some((f) => f.name === def.fluid);

  const unitNames = options.unit_systems.map((u) => u.name);
  const regionNames = options.regions.map((r) => r.name);
  const currentUnits = unitsName(spec.units);
  const regionName = typeof spec.limits === 'string' ? spec.limits : '';
  const look = spec.settings?.look ?? THEME_AUTO;
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
            onChange={(e) => onChange({ ...def, fluid: e.target.value, spec: withView(spec, null) })}
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

        <label className="diagram-field" title="Diagram type">
          <span>Diagram</span>
          <select
            className="diagram-select"
            value={spec.diagram}
            onChange={(e) => {
              const to = diagramOption(options, e.target.value);
              if (to) setSpec(withDiagramType(spec, type, to));
            }}
          >
            {!type && <option value={spec.diagram}>{spec.diagram}</option>}
            {options.diagram_types.map((d) => (
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

        <label className="diagram-field" title="Colours and fonts of the chart">
          <span>Theme</span>
          <select className="diagram-select" value={look} onChange={(e) => setLook(e.target.value)}>
            <option value={THEME_AUTO}>as the GUI</option>
            {look !== THEME_AUTO && !options.themes.includes(look) && <option value={look}>{look}</option>}
            {options.themes.map((t) => <option key={t} value={t}>{t}</option>)}
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
        {(type?.families ?? []).map((family) => {
          const name = family.key === 'rho' ? 'ρ' : family.symbol || symbol(family.key);
          const set = spec.isolines.find((s) => s.key === family.key);
          return (
            <span key={family.key} className={`diagram-family${set ? ' on' : ''}`}>
              <label className="diagram-field" title={`Lines of constant ${(family.label || name).toLowerCase()}`}>
                <input
                  type="checkbox"
                  checked={!!set}
                  onChange={(e) => setSpec(withFamily(spec, family.key, { on: e.target.checked }, family.default_num))}
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
      </div>
    </div>
  );
}
