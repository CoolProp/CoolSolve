/**
 * CycleEditor — the user's cycles of a diagram definition (docs/coolplot_integration.md §4.1, source 3,
 * `CS-FEAT-OVERLAY-STATES`): numbered states, each a set of the model's variables by CoolPlot property
 * key, chains of states ("1-2-3-4-1; 2-5") and the paths between them (one kind for the cycle, and the
 * connections that differ). "From inferred states" starts a cycle from the server's grouping of the
 * solved variables. Text fields are applied on Enter or when they lose the focus.
 */

import { useState } from 'react';
import type { KeyboardEvent } from 'react';
import { ChevronDown, ChevronRight, Plus, Trash2, X } from 'lucide-react';
import type { CycleDefinition, CyclePathKind, CycleStateDefinition, DiagramDefinition } from '../../api/types';
import {
  ARRAYS_CYCLE,
  INFERRED_CYCLE,
  PATH_KINDS,
  STATE_KEYS,
  connectionsOf,
  formatLoops,
  freeCycleName,
  isInferredCycle,
  isUserCycle,
  newUserCycle,
  parseLoops,
  pathsOf,
  pathsValue,
} from './definition';
import { symbolOf } from './options';

/** A solved variable the editor offers, with the property the inference found for the chart's fluid. */
export interface VariableChoice {
  name: string;
  key?: string;
  ofFluid: boolean;
}

interface Props {
  def: DiagramDefinition;
  onChange: (def: DiagramDefinition) => void;
  /** The server's states of the chart's fluid in the last solve; null while unknown */
  inferredStates: CycleStateDefinition[] | null;
  variables: VariableChoice[];
}

const listId = (key: string) => `cs-cycle-vars-${key}`;
/** Properties in the order of STATE_KEYS (T, p, h, s, ...), unknown ones last. */
const keyRank = (key: string) => (STATE_KEYS.includes(key) ? STATE_KEYS.indexOf(key) : STATE_KEYS.length);

/**
 * A text field applied on Enter or blur. `onCommit` returns an error to show (the text is kept), or
 * null. The parent gives it `key={value}`, so that it shows the value again after a change.
 */
function CommitInput({ value, onCommit, className, placeholder, list, title }: {
  value: string;
  onCommit: (text: string) => string | null;
  className?: string;
  placeholder?: string;
  list?: string;
  title?: string;
}) {
  const [text, setText] = useState(value);
  const [error, setError] = useState('');
  const commit = () => {
    if (text === value) { setError(''); return; }
    setError(onCommit(text.trim()) ?? '');
  };
  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter') (e.target as HTMLInputElement).blur();
    if (e.key === 'Escape') { setText(value); setError(''); }
  };
  return (
    <input
      className={`diagram-input${className ? ` ${className}` : ''}${error ? ' invalid' : ''}`}
      value={text}
      placeholder={placeholder}
      list={list}
      title={error || title}
      onChange={(e) => setText(e.target.value)}
      onBlur={commit}
      onKeyDown={onKeyDown}
    />
  );
}

/** Replace state `from` by `to` in the loops and in the keys of the path overrides. */
function renumber(cycle: CycleDefinition, from: number, to: number): CycleDefinition {
  const { kind, overrides } = pathsOf(cycle);
  const loops = (cycle.loops ?? []).map((chain) => chain.map((n) => (n === from ? to : n)));
  const moved: Record<string, CyclePathKind> = {};
  for (const [key, k] of Object.entries(overrides)) {
    const [a, b] = key.split('-').map(Number);
    moved[`${a === from ? to : a}-${b === from ? to : b}`] = k;
  }
  return { ...cycle, loops, paths: pathsValue(kind, moved, loops) };
}

function StateRow({ state, onChange, onRemove, numbers, variables }: {
  state: CycleStateDefinition;
  onChange: (state: CycleStateDefinition) => string | null;
  onRemove: () => void;
  numbers: number[];
  variables: VariableChoice[];
}) {
  const used = Object.keys(state.vars);
  const free = STATE_KEYS.filter((k) => !used.includes(k));
  const [slotKey, setSlotKey] = useState('');
  const newKey = free.includes(slotKey) ? slotKey : free[0];

  const setVar = (key: string, name: string): string | null => {
    const vars = { ...state.vars };
    if (name) vars[key] = name;
    else delete vars[key];
    return onChange({ ...state, vars });
  };
  const renameKey = (from: string, to: string) => {
    const vars: Record<string, string> = {};
    for (const [k, v] of Object.entries(state.vars)) vars[k === from ? to : k] = v;
    onChange({ ...state, vars });
  };
  const addVar = (name: string): string | null => {
    if (!name || !newKey) return null;
    // The property the inference found for the variable, when the slot's own is not chosen
    const inferred = variables.find((v) => v.name === name)?.key;
    const key = !slotKey && inferred && free.includes(inferred) ? inferred : newKey;
    setSlotKey('');
    return setVar(key, name);
  };

  return (
    <div className={`diagram-state${used.length < 2 ? ' incomplete' : ''}`}
         title={used.length < 2 ? 'A state needs two properties' : undefined}>
      <CommitInput
        key={`n${state.number}`}
        className="narrow"
        value={String(state.number)}
        title="State number"
        onCommit={(text) => {
          const n = Number(text);
          if (!Number.isInteger(n)) return 'A state number is an integer';
          if (n !== state.number && numbers.includes(n)) return `State ${n} exists`;
          return onChange({ ...state, number: n });
        }}
      />
      <CommitInput
        key={`l${state.label ?? ''}`}
        className="label"
        value={state.label ?? ''}
        placeholder="label"
        title="Written next to the state's number"
        onCommit={(text) => {
          const { label: _label, ...rest } = state;
          return onChange(text ? { ...rest, label: text } : rest);
        }}
      />
      <div className="diagram-pairs">
        {Object.entries(state.vars).sort(([a], [b]) => keyRank(a) - keyRank(b)).map(([key, name]) => (
          <span key={key} className="diagram-pair">
            <select className="diagram-select" value={key} title="Property" onChange={(e) => renameKey(key, e.target.value)}>
              {[key, ...free].map((k) => <option key={k} value={k}>{symbolOf(k)}</option>)}
            </select>
            <CommitInput key={`${key}=${name}`} value={name} list={listId(key)} title="Variable of the model"
                         onCommit={(text) => setVar(key, text)} />
            <button className="diagram-icon-btn" title="Remove this property" onClick={() => setVar(key, '')}>
              <X size={11} />
            </button>
          </span>
        ))}
        {newKey && (
          <span className="diagram-pair new">
            <select className="diagram-select" value={newKey} title="Property" onChange={(e) => setSlotKey(e.target.value)}>
              {free.map((k) => <option key={k} value={k}>{symbolOf(k)}</option>)}
            </select>
            <CommitInput key={`new-${used.join()}`} value="" list={listId(newKey)} placeholder="variable"
                         title="Add a variable of the model" onCommit={addVar} />
          </span>
        )}
      </div>
      <button className="diagram-icon-btn" title="Remove this state" onClick={onRemove}><Trash2 size={12} /></button>
    </div>
  );
}

function CycleBlock({ cycle, others, onChange, onRemove, variables }: {
  cycle: CycleDefinition;
  others: string[];
  onChange: (cycle: CycleDefinition) => void;
  onRemove: () => void;
  variables: VariableChoice[];
}) {
  const [perConnection, setPerConnection] = useState(false);
  const states = cycle.states ?? [];
  const numbers = states.map((s) => s.number);
  const { kind, overrides } = pathsOf(cycle);
  const connections = connectionsOf(cycle.loops);

  const setPaths = (k: CyclePathKind, o: Record<string, CyclePathKind>) =>
    onChange({ ...cycle, paths: pathsValue(k, o, cycle.loops) });

  const setState = (index: number, next: CycleStateDefinition): string | null => {
    const before = states[index];
    let updated: CycleDefinition = { ...cycle, states: states.map((s, i) => (i === index ? next : s)) };
    if (next.number !== before.number) updated = renumber(updated, before.number, next.number);
    onChange(updated);
    return null;
  };
  const removeState = (index: number) => {
    const n = states[index].number;
    const loops = (cycle.loops ?? []).map((chain) => chain.filter((m) => m !== n)).filter((chain) => chain.length >= 2);
    onChange({ ...cycle, states: states.filter((_, i) => i !== index), loops, paths: pathsValue(kind, overrides, loops) });
  };
  const addState = () => {
    const number = numbers.length ? Math.max(...numbers) + 1 : 1;
    onChange({ ...cycle, states: [...states, { number, vars: {} }] });
  };

  return (
    <div className="diagram-cycle">
      <div className="diagram-row">
        <CommitInput
          key={`name-${cycle.name}`}
          className="name"
          value={cycle.name}
          title="Name of the cycle (legend)"
          onCommit={(text) => {
            if (!text) return 'A cycle needs a name';
            if (others.includes(text)) return 'Another cycle has this name';
            if (text === INFERRED_CYCLE || text === ARRAYS_CYCLE) return 'This name belongs to an overlay';
            onChange({ ...cycle, name: text });
            return null;
          }}
        />
        <label className="diagram-field" title="How the states are joined; auto: isobar, isenthalp, isentrope, isotherm or polytropic, from the two states">
          paths
          <select className="diagram-select" value={kind} onChange={(e) => setPaths(e.target.value as CyclePathKind, overrides)}>
            {PATH_KINDS.map((k) => <option key={k} value={k}>{k}</option>)}
          </select>
        </label>
        <label className="diagram-field" title="Write the state numbers">
          <input type="checkbox" checked={cycle.numbers !== false}
                 onChange={(e) => onChange({ ...cycle, numbers: e.target.checked })} />
          numbers
        </label>
        <button className="diagram-icon-btn" title="Delete this cycle" onClick={onRemove}><Trash2 size={12} /></button>
      </div>
      <div className="diagram-row">
        <span>loops</span>
        <CommitInput
          key={`loops-${formatLoops(cycle.loops)}`}
          className="loops"
          value={formatLoops(cycle.loops)}
          placeholder="1-2-3-4-1; 2-5 (empty: markers only)"
          title="Chains of state numbers, separated by ';'"
          onCommit={(text) => {
            const parsed = parseLoops(text, numbers);
            if ('error' in parsed) return parsed.error;
            onChange({ ...cycle, loops: parsed.loops, paths: pathsValue(kind, overrides, parsed.loops) });
            return null;
          }}
        />
        {connections.length > 0 && (
          <button className="diagram-link" onClick={() => setPerConnection((v) => !v)}>
            {perConnection ? 'hide paths per connection' : `paths per connection${Object.keys(overrides).length ? ` (${Object.keys(overrides).length})` : ''}`}
          </button>
        )}
      </div>
      {perConnection && connections.length > 0 && (
        <div className="diagram-row">
          {connections.map(([a, b]) => {
            const key = `${a}-${b}`;
            return (
              <label key={key} className="diagram-field">
                {a}→{b}
                <select
                  className="diagram-select"
                  value={overrides[key] ?? ''}
                  onChange={(e) => {
                    const next = { ...overrides };
                    if (e.target.value) next[key] = e.target.value as CyclePathKind;
                    else delete next[key];
                    setPaths(kind, next);
                  }}
                >
                  <option value="">({kind})</option>
                  {PATH_KINDS.filter((k) => k !== kind).map((k) => <option key={k} value={k}>{k}</option>)}
                </select>
              </label>
            );
          })}
        </div>
      )}
      {states.map((state, i) => (
        <StateRow
          key={state.number}
          state={state}
          numbers={numbers}
          variables={variables}
          onChange={(next) => setState(i, next)}
          onRemove={() => removeState(i)}
        />
      ))}
      <button className="diagram-link" onClick={addState}><Plus size={11} /> state</button>
    </div>
  );
}

export default function CycleEditor({ def, onChange, inferredStates, variables }: Props) {
  const [open, setOpen] = useState(false);
  const cycles = def.cycles.map((cycle, index) => ({ cycle, index })).filter(({ cycle }) => isUserCycle(cycle));

  // A cycle of the user supersedes the inferred states, which would mark its states twice: adding the
  // first one switches them off (the checkbox brings them back)
  const add = (cycle: CycleDefinition) => {
    const kept = cycles.length === 0 ? def.cycles.filter((c) => !isInferredCycle(c)) : def.cycles;
    onChange({ ...def, cycles: [...kept, cycle] });
    setOpen(true);
  };
  const setCycle = (index: number, cycle: CycleDefinition) =>
    onChange({ ...def, cycles: def.cycles.map((c, i) => (i === index ? cycle : c)) });
  const removeCycle = (index: number) => onChange({ ...def, cycles: def.cycles.filter((_, i) => i !== index) });

  // Per property: the variables the inference found for it first, then the fluid's, then the others
  const lists = open ? STATE_KEYS.map((key) => {
    const rank = (v: VariableChoice) => (v.key === key ? 0 : v.ofFluid ? 1 : 2);
    return { key, names: [...variables].sort((a, b) => rank(a) - rank(b)).map((v) => v.name) };
  }) : [];

  return (
    <div className="diagram-cycles">
      <div className="diagram-row">
        <button className="diagram-link title" onClick={() => setOpen((v) => !v)}>
          {open ? <ChevronDown size={12} /> : <ChevronRight size={12} />}
          Cycles{cycles.length ? ` (${cycles.length})` : ''}
        </button>
        <button className="diagram-link" title="A new cycle; add its states below"
                onClick={() => add(newUserCycle(freeCycleName(def.cycles, 'Cycle'), []))}>
          <Plus size={11} /> new
        </button>
        <button
          className="diagram-link"
          disabled={!inferredStates || inferredStates.length === 0}
          title="A cycle through the inferred states of the solved model, in order, to edit"
          onClick={() => add(newUserCycle(freeCycleName(def.cycles, 'Cycle'), inferredStates ?? []))}
        >
          <Plus size={11} /> from inferred states
        </button>
      </div>
      {open && (
        <>
          {lists.map(({ key, names }) => (
            <datalist key={key} id={listId(key)}>
              {names.map((name) => <option key={name} value={name} />)}
            </datalist>
          ))}
          {cycles.map(({ cycle, index }) => (
            <CycleBlock
              key={index}
              cycle={cycle}
              others={def.cycles.filter((_, i) => i !== index).map((c) => c.name)}
              variables={variables}
              onChange={(next) => setCycle(index, next)}
              onRemove={() => removeCycle(index)}
            />
          ))}
        </>
      )}
    </div>
  );
}
