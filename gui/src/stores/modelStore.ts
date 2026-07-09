import { create } from 'zustand';
import type { SolveResponse, ParseError, SolvedVariable, VariableInfo, SavedParametricStudy, UserUnitOverride, LookupTableInfo, IntegralTableData } from '../api/types';

interface ModelState {
  // File state
  filePath: string;
  modelName: string;
  eescode: string;
  initials: string;
  sol: string;
  conf: string;
  dirty: boolean;
  canGoBack: boolean;

  // Parse state
  parseErrors: ParseError[];
  equationCount: number;
  variableCount: number;
  isSquare: boolean;
  parsedVariables: VariableInfo[];  // Variable info from last parse (units, imposed, etc.)

  // Solve state
  solving: boolean;
  lastResult: SolveResponse | null;
  solvedVariables: SolvedVariable[];

  // Parametric studies (saved results)
  parametricStudies: SavedParametricStudy[];

  // Integral (INTEGRAL/$IntegralTable) trajectory
  integralTable: IntegralTableData | null;  // columnar table (live solve)
  integralCSV: string;                       // raw CSV text (round-trip via ZIP)

  // Lookup tables
  lookupTables: LookupTableInfo[];
  lookupTableCSVs: Record<string, string>;  // name → raw CSV content

  // User unit overrides (GUI-specified units)
  userUnitOverrides: UserUnitOverride[];

  // Console log
  consoleLines: string[];

  // Actions
  setEescode: (code: string) => void;
  setInitials: (initials: string) => void;
  setConf: (conf: string) => void;
  setFilePath: (path: string) => void;
  setModelName: (name: string) => void;
  setParseResult: (errors: ParseError[], eqCount: number, varCount: number, isSquare: boolean) => void;
  setParsedVariables: (vars: VariableInfo[]) => void;
  setSolving: (solving: boolean) => void;
  setSolveResult: (result: SolveResponse) => void;
  setSolvedVariables: (vars: SolvedVariable[]) => void;
  setSol: (sol: string) => void;
  setCanGoBack: (canGoBack: boolean) => void;
  addConsoleLine: (line: string) => void;
  clearConsole: () => void;
  loadFile: (path: string, eescode: string, initials: string, sol: string, conf: string, modelName?: string) => void;
  clearModel: () => void;
  addParametricStudy: (study: SavedParametricStudy) => void;
  removeParametricStudy: (id: string) => void;
  setParametricStudies: (studies: SavedParametricStudy[]) => void;
  setIntegralTable: (table: IntegralTableData | null) => void;
  setIntegralCSV: (csv: string) => void;
  setUserUnitOverrides: (overrides: UserUnitOverride[]) => void;
  setUserUnit: (variableName: string, units: string) => void;
  setLookupTables: (tables: LookupTableInfo[]) => void;
  setLookupTableCSV: (name: string, csv: string) => void;
  removeLookupTable: (name: string) => void;
}

export const useModelStore = create<ModelState>((set) => ({
  filePath: '',
  modelName: '',
  eescode: '',
  initials: '',
  sol: '',
  conf: '',
  dirty: false,
  canGoBack: false,

  parseErrors: [],
  equationCount: 0,
  variableCount: 0,
  isSquare: true,
  parsedVariables: [],

  solving: false,
  lastResult: null,
  solvedVariables: [],

  parametricStudies: [],
  userUnitOverrides: [],
  lookupTables: [],
  lookupTableCSVs: {},
  integralTable: null,
  integralCSV: '',

  consoleLines: [],

  setEescode: (code) => set({ eescode: code, dirty: true }),
  setInitials: (initials) => set({ initials }),
  setConf: (conf) => set({ conf }),
  setFilePath: (path) => set({ filePath: path }),
  setModelName: (name) => set({ modelName: name }),
  setParseResult: (errors, eqCount, varCount, isSquare) =>
    set({ parseErrors: errors, equationCount: eqCount, variableCount: varCount, isSquare }),
  setParsedVariables: (vars) => set({ parsedVariables: vars }),
  setSolving: (solving) => set({ solving }),
  setSolveResult: (result) => set({ lastResult: result }),
  setSolvedVariables: (vars) => set({ solvedVariables: vars }),
  setSol: (sol) => set({ sol }),
  setCanGoBack: (canGoBack) => set({ canGoBack }),
  addConsoleLine: (line) =>
    set((state) => ({ consoleLines: [...state.consoleLines, line] })),
  clearConsole: () => set({ consoleLines: [] }),
  loadFile: (path, eescode, initials, sol, conf, modelName) =>
    set({ filePath: path, eescode, initials, sol, conf, dirty: false, lastResult: null, solvedVariables: [], consoleLines: [], modelName: modelName ?? '', parametricStudies: [], parsedVariables: [], lookupTables: [], lookupTableCSVs: {}, integralTable: null, integralCSV: '' }),
  clearModel: () =>
    set({ filePath: '', modelName: '', eescode: '', initials: '', sol: '', conf: '', dirty: false, lastResult: null, solvedVariables: [], consoleLines: [], parseErrors: [], equationCount: 0, variableCount: 0, isSquare: true, parametricStudies: [], parsedVariables: [], userUnitOverrides: [], lookupTables: [], lookupTableCSVs: {}, integralTable: null, integralCSV: '' }),
  addParametricStudy: (study) =>
    set((state) => ({ parametricStudies: [...state.parametricStudies, study] })),
  removeParametricStudy: (id) =>
    set((state) => ({ parametricStudies: state.parametricStudies.filter((s) => s.id !== id) })),
  setParametricStudies: (studies) => set({ parametricStudies: studies }),
  setIntegralTable: (table) => set({ integralTable: table }),
  setIntegralCSV: (csv) => set({ integralCSV: csv }),
  setUserUnitOverrides: (overrides) => set({ userUnitOverrides: overrides }),
  setUserUnit: (variableName, units) =>
    set((state) => {
      const existing = state.userUnitOverrides.filter((o) => o.variableName !== variableName);
      if (units) existing.push({ variableName, units });
      return { userUnitOverrides: existing };
    }),
  setLookupTables: (tables) => set({ lookupTables: tables }),
  setLookupTableCSV: (name, csv) =>
    set((state) => ({
      lookupTableCSVs: { ...state.lookupTableCSVs, [name]: csv },
      lookupTables: state.lookupTables.some((t) => t.name === name)
        ? state.lookupTables
        : [...state.lookupTables, { name, columns: [], rows: 0 }],
    })),
  removeLookupTable: (name) =>
    set((state) => {
      const { [name]: _removed, ...rest } = state.lookupTableCSVs;
      return { lookupTableCSVs: rest, lookupTables: state.lookupTables.filter((t) => t.name !== name) };
    }),
}));