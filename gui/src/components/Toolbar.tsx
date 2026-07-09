import { useState, useEffect, useCallback, useRef } from 'react';
import {
  Play, FolderOpen, Save, RefreshCw, Bug, Sun, Moon,
  ChevronDown, ChevronUp, BookOpen, Square,
  Braces, Quote, FilePlus, Undo, Pencil, FileText, Library,
} from 'lucide-react';
import { useModelStore } from '../stores/modelStore';
import { useUIStore } from '../stores/uiStore';
import { api } from '../api/client';
import { editorInstance, toggleBraceComment, toggleQuoteComment } from './CodeEditor';
import { exportAllPlots } from '../utils/exportPlots';
import type { ExampleFile, SSEEvent, SolveResponse } from '../api/types';

export default function Toolbar() {
  const solving = useModelStore((s) => s.solving);
  const setSolving = useModelStore((s) => s.setSolving);
  const setSolveResult = useModelStore((s) => s.setSolveResult);
  const addConsoleLine = useModelStore((s) => s.addConsoleLine);
  const clearConsole = useModelStore((s) => s.clearConsole);
  const eescode = useModelStore((s) => s.eescode);
  const initials = useModelStore((s) => s.initials);
  const sol = useModelStore((s) => s.sol);
  const conf = useModelStore((s) => s.conf);
  const modelName = useModelStore((s) => s.modelName);
  const setModelName = useModelStore((s) => s.setModelName);
  const setLookupTables = useModelStore((s) => s.setLookupTables);
  const setIntegralTable = useModelStore((s) => s.setIntegralTable);
  const setIntegralCSV = useModelStore((s) => s.setIntegralCSV);
  const canGoBack = useModelStore((s) => s.canGoBack);
  const setSol = useModelStore((s) => s.setSol);
  const setInitials = useModelStore((s) => s.setInitials);
  const setCanGoBack = useModelStore((s) => s.setCanGoBack);
  const loadFile = useModelStore((s) => s.loadFile);
  const clearModel = useModelStore((s) => s.clearModel);
  const lastResult = useModelStore((s) => s.lastResult);

  const theme = useUIStore((s) => s.theme);
  const toggleTheme = useUIStore((s) => s.toggleTheme);
  const bottomPanelOpen = useUIStore((s) => s.bottomPanelOpen);
  const toggleBottomPanel = useUIStore((s) => s.toggleBottomPanel);
  const setBottomPanelOpen = useUIStore((s) => s.setBottomPanelOpen);

  const [examples, setExamples] = useState<ExampleFile[]>([]);
  const [showExamples, setShowExamples] = useState(false);
  const [editingName, setEditingName] = useState(false);
  const [nameInput, setNameInput] = useState('');
  const nameInputRef = useRef<HTMLInputElement>(null);
  const eventSourceRef = useRef<EventSource | null>(null);

  // Update page title when model name changes
  useEffect(() => {
    document.title = modelName ? `CoolSolve: ${modelName}` : 'CoolSolve';
  }, [modelName]);

  // Focus name input when editing starts
  useEffect(() => {
    if (editingName && nameInputRef.current) {
      nameInputRef.current.focus();
      nameInputRef.current.select();
    }
  }, [editingName]);

  // Load examples on mount
  useEffect(() => {
    api.getExamples().then((res) => setExamples(res.examples)).catch(() => {});
  }, []);

  // Cleanup SSE on unmount
  useEffect(() => {
    return () => {
      if (eventSourceRef.current) eventSourceRef.current.close();
    };
  }, []);

  // Solve handler — async with SSE progress streaming
  const handleSolve = useCallback(
    async (debug = false) => {
      if (!eescode.trim()) return;
      setSolving(true);
      clearConsole();
      setBottomPanelOpen(true);
      addConsoleLine('>>> Starting solve...');

      try {
        // Save current content to backend
        await api.putEescode(eescode);
        if (initials) await api.putInitials(initials);

        // Trigger async solve
        await api.solve({ debug });

        // Subscribe to SSE progress stream
        if (eventSourceRef.current) eventSourceRef.current.close();
        const es = api.subscribeSolveProgress((event: SSEEvent) => {
          switch (event.type) {
            case 'start':
              addConsoleLine(`>>> ${event.message}`);
              break;
            case 'progress':
              addConsoleLine(`    ${event.message}`);
              break;
            case 'block':
              if (event.event === 'start') {
                addConsoleLine(`    ${event.message}`);
              } else if (event.event === 'done') {
                addConsoleLine(`    ${event.message}`);
              } else if (event.event === 'fail') {
                addConsoleLine(`    ${event.message}`);
              }
              break;
            case 'done': {
              addConsoleLine(`>>> Solve SUCCESS`);
              const result = event.result as SolveResponse;
              if (result) {
                setSolveResult(result);
                addConsoleLine(`    Equations: ${result.equationCount}, Variables: ${result.variableCount}`);
                addConsoleLine(`    Blocks: ${result.totalBlocks}, Largest: ${result.largestBlock}`);
                addConsoleLine(`    Total iterations: ${result.totalIterations}`);
                addConsoleLine(`    Time: ${result.timing.total_ms.toFixed(0)} ms (solve: ${result.timing.solve_ms.toFixed(0)} ms)`);
                addConsoleLine(`    Variables solved: ${Object.keys(result.variables).length}`);
                // Display diagnostics (warnings/info)
                if (result.diagnostics && result.diagnostics.length > 0) {
                  addConsoleLine(`--- Diagnostics (${result.diagnostics.length}) ---`);
                  for (const d of result.diagnostics) {
                    const loc = d.line ? ` (line ${d.line})` : '';
                    addConsoleLine(`    [${d.severity}]${loc}: ${d.message}`);
                  }
                }
                // Fetch .sol content
                api.getSol().then((s) => setSol(s.content)).catch(() => {});
                // Refresh lookup table list (tables may have been created/loaded)
                api.getTables().then((r) => setLookupTables(r.tables)).catch(() => {});
                // Integral trajectory: prefer the columnar JSON embedded in the
                // solve response; fall back to fetching the CSV (e.g. if the
                // result was large and stripped, or for restored sessions).
                if (result.integralTable) {
                  setIntegralTable(result.integralTable);
                  setIntegralCSV('');
                } else {
                  setIntegralTable(null);
                  api.getIntegralCSV().then((csv) => setIntegralCSV(csv ?? '')).catch(() => {});
                }
                // Auto-switch to the Integral tab when a trajectory was produced,
                // so the user immediately sees the result.
                if (result.integralTable) {
                  setBottomPanelOpen(true);
                }
              }
              setSolving(false);
              es.close();
              break;
            }
            case 'error': {
              const result = event.result as SolveResponse | undefined;
              if (result) {
                setSolveResult(result);
                // Provide user-friendly failure messages
                if (result.isSquare === false) {
                  addConsoleLine(`>>> Solve FAILED: System is not square (${result.equationCount} equations, ${result.variableCount} variables)`);
                } else if (result.status === 'ParseFailed') {
                  addConsoleLine(`>>> Solve FAILED: Parse error`);
                } else {
                  addConsoleLine(`>>> Solve FAILED: ${result.errorMessage || result.status}`);
                }
                if (result.detailedError) addConsoleLine(`    ${result.detailedError}`);
                // Show failed blocks (but not for parse/structural failures)
                if (result.status !== 'ParseFailed' && result.status !== 'InvalidInput') {
                  for (const block of result.blockResults) {
                    if (!block.success) {
                      addConsoleLine(
                        `    Block ${block.id}: ${block.status} (res=${block.maxResidual.toExponential(2)}) ${block.errorMessage}`
                      );
                    }
                  }
                }
                // Display diagnostics
                if (result.diagnostics && result.diagnostics.length > 0) {
                  addConsoleLine(`--- Diagnostics (${result.diagnostics.length}) ---`);
                  for (const d of result.diagnostics) {
                    const loc = d.line ? ` (line ${d.line})` : '';
                    addConsoleLine(`    [${d.severity}]${loc}: ${d.message}`);
                  }
                }
              } else {
                addConsoleLine(`>>> ERROR: ${event.message}`);
              }
              setSolving(false);
              es.close();
              break;
            }
          }
        });
        eventSourceRef.current = es;

      } catch (err: any) {
        addConsoleLine(`>>> ERROR: ${err.message}`);
        setSolving(false);
      }
    },
    [eescode, initials, setSolving, clearConsole, addConsoleLine, setSolveResult, setSol, setBottomPanelOpen, setLookupTables, setIntegralTable, setIntegralCSV]
  );

  // ================================================================
  // New — clear everything
  // ================================================================
  const handleNew = useCallback(async () => {
    try {
      const res = await api.newModel();
      clearModel();
      setCanGoBack(res.hadContent);
      addConsoleLine('>>> New model');
    } catch (err: any) {
      addConsoleLine(`>>> ERROR: ${err.message}`);
    }
  }, [clearModel, addConsoleLine, setCanGoBack]);

  // ================================================================
  // Open — upload ZIP file
  // ================================================================
  const handleOpen = useCallback(() => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.zip';
    input.onchange = async (e) => {
      const files = (e.target as HTMLInputElement).files;
      if (!files || files.length === 0) return;
      const file = files[0];
      const formData = new FormData();
      formData.append('file', file);
      try {
        const res = await api.uploadFiles(formData);
        if (res.success) {
          const ees = await api.getEescode();
          const init = await api.getInitials();
          const solRes = await api.getSol();
          const confRes = await api.getConf();
          setCanGoBack(true);
          loadFile(file.name.replace(/\.zip$/, ''), ees.content, init.content, solRes.content, confRes.content, res.modelName);
          // Restore integral trajectory: the columnar JSON is not round-tripped
          // through the bundle, but the CSV is — fetch it so the Integral tab
          // repopulates on bundle load.
          api.getIntegralResult().then((t) => setIntegralTable(t?.integralTable ?? null)).catch(() => setIntegralTable(null));
          api.getIntegralCSV().then((csv) => setIntegralCSV(csv ?? '')).catch(() => setIntegralCSV(''));
          addConsoleLine(`>>> Opened: ${res.files.join(', ')}`);
        }
      } catch (err: any) {
        addConsoleLine(`>>> ERROR: ${err.message}`);
      }
    };
    input.click();
  }, [loadFile, addConsoleLine, setCanGoBack, setIntegralTable, setIntegralCSV]);

  // ================================================================
  // Save — download ZIP bundle
  // ================================================================
  const handleSave = useCallback(async () => {
    if (!eescode.trim()) {
      addConsoleLine('>>> Nothing to save');
      return;
    }
    try {
      // Sync current content to backend before downloading
      await api.putEescode(eescode);
      if (initials) await api.putInitials(initials);

      // If no model name, prompt for one
      let name = modelName;
      if (!name) {
        const input = prompt('Enter model name:', 'model');
        if (!input) return; // cancelled
        name = input.trim() || 'model';
        setModelName(name);
        await api.setModelName(name);
      }

      // Build file list for logging
      const fileNames: string[] = [name + '.eescode'];
      if (initials) fileNames.push(name + '.initials');
      if (sol) fileNames.push(name + '.sol');
      if (conf) fileNames.push('coolsolve.conf');
      // Check for debug files
      try {
        const debugRes = await api.getDebugFiles();
        if (debugRes.files.length > 0) fileNames.push('debug_output/');
      } catch { /* no debug */ }

      addConsoleLine(`>>> Saved: ${fileNames.join(', ')}`);
      window.location.href = '/api/v1/files/bundle';
    } catch (err: any) {
      addConsoleLine(`>>> ERROR: ${err.message}`);
    }
  }, [eescode, initials, sol, conf, modelName, setModelName, addConsoleLine]);

  // ================================================================
  // Back — restore previous model
  // ================================================================
  const handleBack = useCallback(async () => {
    try {
      const res = await api.goBack();
      if (res.success) {
        const ees = await api.getEescode();
        const init = await api.getInitials();
        const solRes = await api.getSol();
        const confRes = await api.getConf();
        loadFile(ees.filePath || '', ees.content, init.content, solRes.content, confRes.content, res.modelName);
        setCanGoBack(false);
        // Try to restore solve result
        try {
          const result = await api.getSolveResult();
          setSolveResult(result);
        } catch { /* no result */ }
        // Restore integral trajectory (CSV-backed; survives the model swap).
        api.getIntegralResult().then((t) => setIntegralTable(t?.integralTable ?? null)).catch(() => setIntegralTable(null));
        api.getIntegralCSV().then((csv) => setIntegralCSV(csv ?? '')).catch(() => setIntegralCSV(''));
        addConsoleLine('>>> Restored previous model');
      }
    } catch (err: any) {
      addConsoleLine(`>>> ERROR: ${err.message}`);
    }
  }, [loadFile, addConsoleLine, setCanGoBack, setSolveResult, setIntegralTable, setIntegralCSV]);

  // Open example by server path
  const openFileByPath = useCallback(async (path: string) => {
    try {
      const res = await api.openFile(path);
      if (res.success) {
        const ees = await api.getEescode();
        const init = await api.getInitials();
        const solRes = await api.getSol();
        const confRes = await api.getConf();
        setCanGoBack(true);
        loadFile(res.filePath, ees.content, init.content, solRes.content, confRes.content, res.modelName);
        addConsoleLine(`>>> Opened example: ${res.modelName || res.filePath.split('/').pop()}`);
      }
    } catch (err: any) {
      addConsoleLine(`>>> ERROR opening example: ${err.message}`);
    }
  }, [loadFile, addConsoleLine, setCanGoBack]);

  // Stop/cancel solve
  const handleStopSolve = useCallback(async () => {
    try {
      await api.cancelSolve();
      addConsoleLine('>>> Cancel requested...');
    } catch (err: any) {
      addConsoleLine(`>>> ${err.message}`);
    }
  }, [addConsoleLine]);

  // Update guesses
  const handleUpdateGuesses = useCallback(async () => {
    try {
      await api.updateGuesses();
      const init = await api.getInitials();
      setInitials(init.content);
      addConsoleLine('>>> Guesses updated from last solution');
    } catch (err: any) {
      addConsoleLine(`>>> ERROR: ${err.message}`);
    }
  }, [setInitials, addConsoleLine]);

  // Comment toggles
  const handleBraceComment = useCallback(() => {
    toggleBraceComment(editorInstance);
  }, []);

  const handleQuoteComment = useCallback(() => {
    toggleQuoteComment(editorInstance);
  }, []);

  // Helper: trigger a file download from a Blob
  const downloadBlob = useCallback((blob: Blob, filename: string) => {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.style.display = 'none';
    document.body.appendChild(a);
    a.click();
    // Small delay before cleanup so the download starts
    setTimeout(() => {
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    }, 100);
  }, []);

  // Compile LaTeX report to PDF (with embedded plot images)
  const handleLatexReport = useCallback(async () => {
    try {
      addConsoleLine('>>> Exporting plots...');
      let plots: { name: string; data: string }[] = [];
      try {
        plots = await exportAllPlots();
      } catch (plotErr: any) {
        addConsoleLine(`    Plot export error: ${plotErr.message} (continuing without plots)`);
      }
      if (plots.length > 0) {
        addConsoleLine(`    Captured ${plots.length} plot(s): ${plots.map((p) => p.name).join(', ')}`);
      } else {
        addConsoleLine('    No visible plots to embed');
      }

      addConsoleLine('>>> Compiling LaTeX report to PDF...');
      const pdfBlob = await api.compileLatexReport({ plots });
      addConsoleLine(`    Received ${pdfBlob.size} bytes (type: ${pdfBlob.type})`);

      const filename = (modelName || 'model') + '_report.pdf';
      downloadBlob(pdfBlob, filename);
      addConsoleLine(`>>> PDF report downloaded: ${filename}`);
    } catch (err: any) {
      // If compilation fails (e.g. pdflatex not installed), offer the .tex file
      addConsoleLine(`>>> PDF compilation failed: ${err.message || err}`);
      addConsoleLine('>>> Falling back to .tex download...');
      try {
        const res = await api.getLatexReport();
        if (res.available && res.content) {
          const blob = new Blob([res.content], { type: 'application/x-tex' });
          const filename = (modelName || 'model') + '_report.tex';
          downloadBlob(blob, filename);
          addConsoleLine(`>>> LaTeX source downloaded: ${filename}`);
        } else {
          addConsoleLine('>>> No LaTeX report available');
        }
      } catch {
        addConsoleLine('>>> Could not download .tex fallback either');
      }
    }
  }, [modelName, addConsoleLine, downloadBlob]);

  // Open example
  const handleOpenExample = useCallback(
    async (ex: ExampleFile) => {
      setShowExamples(false);
      await openFileByPath(ex.path);
    },
    [openFileByPath]
  );

  // Keyboard shortcuts
  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
        e.preventDefault();
        if (!solving) handleSolve(e.shiftKey);
      } else if ((e.ctrlKey || e.metaKey) && e.key === 's') {
        e.preventDefault();
        handleSave();
      } else if ((e.ctrlKey || e.metaKey) && e.key === 'o') {
        e.preventDefault();
        handleOpen();
      } else if ((e.ctrlKey || e.metaKey) && e.key === 'n') {
        e.preventDefault();
        handleNew();
      } else if ((e.ctrlKey || e.metaKey) && e.key === 'g') {
        e.preventDefault();
        handleUpdateGuesses();
      } else if (e.key === 'Escape' && solving) {
        e.preventDefault();
        handleStopSolve();
      }
    };
    window.addEventListener('keydown', handler);
    return () => window.removeEventListener('keydown', handler);
  }, [solving, handleSolve, handleSave, handleOpen, handleNew, handleUpdateGuesses, handleStopSolve]);

  return (
    <div className="toolbar">
      {/* File group */}
      <div className="toolbar-group">
        <button className="toolbar-btn" onClick={handleNew} title="New model (Ctrl+N)">
          <FilePlus size={16} /> New
        </button>
        <button className="toolbar-btn" onClick={handleOpen} title="Open ZIP file (Ctrl+O)">
          <FolderOpen size={16} /> Open
        </button>
        <button className="toolbar-btn" onClick={handleSave} disabled={!eescode.trim()} title="Save as ZIP (Ctrl+S)">
          <Save size={16} /> Save
        </button>
        <button className="toolbar-btn" onClick={handleBack} disabled={!canGoBack} title="Restore previous model">
          <Undo size={16} /> Back
        </button>
      </div>

      <div className="toolbar-separator" />

      {/* Examples */}
      <div className="toolbar-group example-selector">
        <button className="toolbar-btn" onClick={() => setShowExamples(!showExamples)}>
          <BookOpen size={16} /> Examples
        </button>
        {showExamples && examples.length > 0 && (
          <div className="example-dropdown">
            {examples.map((ex) => (
              <button
                key={ex.name}
                className="example-dropdown-item"
                onClick={() => handleOpenExample(ex)}
              >
                {ex.name}
              </button>
            ))}
          </div>
        )}
      </div>

      <div className="toolbar-separator" />

      {/* Solve group */}
      <div className="toolbar-group">
        <button
          className="toolbar-btn primary"
          onClick={() => handleSolve(false)}
          disabled={solving || !eescode.trim()}
          title="Solve (Ctrl+Enter)"
        >
          {solving ? <span className="spinner" /> : <Play size={16} />}
          {solving ? 'Solving...' : 'Solve'}
        </button>
        <button
          className="toolbar-btn"
          onClick={() => handleSolve(true)}
          disabled={solving || !eescode.trim()}
          title="Debug Solve (Ctrl+Shift+Enter)"
        >
          <Bug size={16} /> Debug
        </button>
        {solving && (
          <button
            className="toolbar-btn danger"
            onClick={handleStopSolve}
            title="Stop solve (Escape)"
          >
            <Square size={16} /> Stop
          </button>
        )}
      </div>

      <div className="toolbar-separator" />

      {/* Guesses */}
      <div className="toolbar-group">
        <button className="toolbar-btn" onClick={handleUpdateGuesses} title="Update guesses (Ctrl+G)">
          <RefreshCw size={16} /> Update Guesses
        </button>
      </div>

      <div className="toolbar-separator" />

      {/* Comment toggles */}
      <div className="toolbar-group">
        <button className="toolbar-btn" onClick={handleBraceComment} title="Toggle { } comment (Ctrl+/)">
          <Braces size={16} /> {'{}'}
        </button>
        <button className="toolbar-btn" onClick={handleQuoteComment} title='Toggle &quot; &quot; comment (Ctrl+Shift+/)'>
          <Quote size={16} /> &quot;&quot;
        </button>
      </div>

      <div className="toolbar-separator" />

      {/* LaTeX report */}
      <div className="toolbar-group">
        <button
          className="toolbar-btn"
          onClick={handleLatexReport}
          disabled={!lastResult?.success || !lastResult?.latexReportAvailable}
          title={
            !lastResult?.success
              ? 'Solve a model first to generate a LaTeX report'
              : !lastResult?.latexReportAvailable
                ? 'Enable enableLatexReport in coolsolve.conf to generate LaTeX reports'
                : 'Compile and download PDF report (with plots)'
          }
        >
          <FileText size={16} /> LaTeX
        </button>
      </div>

      <div className="toolbar-separator" />

      {/* Language Reference */}
      <div className="toolbar-group">
        <a
          className="toolbar-btn"
          href="/docs/#docs%2Flanguage_reference.md"
          target="_blank"
          rel="noopener noreferrer"
          title="Open Language Reference (syntax, functions, fluids)"
        >
          <Library size={16} /> Language Reference
        </a>
      </div>

      <div className="toolbar-spacer" />

      {/* Status — editable model name */}
      <div className="toolbar-status" title="Click to rename model">
        {editingName ? (
          <input
            ref={nameInputRef}
            className="model-name-input"
            value={nameInput}
            onChange={(e) => setNameInput(e.target.value)}
            onBlur={() => {
              const trimmed = nameInput.trim();
              if (trimmed && trimmed !== modelName) {
                setModelName(trimmed);
                api.setModelName(trimmed).catch(() => {});
              }
              setEditingName(false);
            }}
            onKeyDown={(e) => {
              if (e.key === 'Enter') {
                (e.target as HTMLInputElement).blur();
              } else if (e.key === 'Escape') {
                setEditingName(false);
              }
            }}
          />
        ) : (
          <span
            className="model-name-display"
            onClick={() => {
              setNameInput(modelName || '');
              setEditingName(true);
            }}
          >
            {modelName || 'Untitled'} <Pencil size={12} />
          </span>
        )}
      </div>

      <div className="toolbar-separator" />

      {/* Theme & panels */}
      <button className="toolbar-btn" onClick={toggleTheme} title="Toggle theme">
        {theme === 'dark' ? <Sun size={16} /> : <Moon size={16} />}
      </button>
      <button className="toolbar-btn" onClick={toggleBottomPanel} title="Toggle console">
        {bottomPanelOpen ? <ChevronDown size={16} /> : <ChevronUp size={16} />}
      </button>
    </div>
  );
}
