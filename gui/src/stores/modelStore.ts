import { create } from 'zustand';
import type { SolveResponse, ParseError, SolvedVariable } from '../api/types';

interface ModelState {
  // File state
  filePath: string;
  eescode: string;
  initials: string;
  sol: string;
  conf: string;
  dirty: boolean;

  // Parse state
  parseErrors: ParseError[];
  equationCount: number;
  variableCount: number;
  isSquare: boolean;

  // Solve state
  solving: boolean;
  lastResult: SolveResponse | null;
  solvedVariables: SolvedVariable[];

  // Console log
  consoleLines: string[];

  // Actions
  setEescode: (code: string) => void;
  setInitials: (initials: string) => void;
  setConf: (conf: string) => void;
  setFilePath: (path: string) => void;
  setParseResult: (errors: ParseError[], eqCount: number, varCount: number, isSquare: boolean) => void;
  setSolving: (solving: boolean) => void;
  setSolveResult: (result: SolveResponse) => void;
  setSolvedVariables: (vars: SolvedVariable[]) => void;
  setSol: (sol: string) => void;
  addConsoleLine: (line: string) => void;
  clearConsole: () => void;
  loadFile: (path: string, eescode: string, initials: string, sol: string, conf: string) => void;
}

export const useModelStore = create<ModelState>((set) => ({
  filePath: '',
  eescode: '',
  initials: '',
  sol: '',
  conf: '',
  dirty: false,

  parseErrors: [],
  equationCount: 0,
  variableCount: 0,
  isSquare: true,

  solving: false,
  lastResult: null,
  solvedVariables: [],

  consoleLines: [],

  setEescode: (code) => set({ eescode: code, dirty: true }),
  setInitials: (initials) => set({ initials }),
  setConf: (conf) => set({ conf }),
  setFilePath: (path) => set({ filePath: path }),
  setParseResult: (errors, eqCount, varCount, isSquare) =>
    set({ parseErrors: errors, equationCount: eqCount, variableCount: varCount, isSquare }),
  setSolving: (solving) => set({ solving }),
  setSolveResult: (result) => set({ lastResult: result }),
  setSolvedVariables: (vars) => set({ solvedVariables: vars }),
  setSol: (sol) => set({ sol }),
  addConsoleLine: (line) =>
    set((state) => ({ consoleLines: [...state.consoleLines, line] })),
  clearConsole: () => set({ consoleLines: [] }),
  loadFile: (path, eescode, initials, sol, conf) =>
    set({ filePath: path, eescode, initials, sol, conf, dirty: false, lastResult: null, solvedVariables: [], consoleLines: [] }),
}));
