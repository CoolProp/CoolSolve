import { useState, useEffect, useCallback, useRef } from 'react';
import {
  Play, FolderOpen, Save, RefreshCw, Bug, Sun, Moon,
  ChevronDown, ChevronUp, BookOpen, Square, SaveAll,
  Braces, Quote,
} from 'lucide-react';
import { useModelStore } from '../stores/modelStore';
import { useUIStore } from '../stores/uiStore';
import { api } from '../api/client';
import { editorInstance, toggleBraceComment, toggleQuoteComment } from './CodeEditor';
import type { ExampleFile, SSEEvent, SolveResponse } from '../api/types';

export default function Toolbar() {
  const solving = useModelStore((s) => s.solving);
  const setSolving = useModelStore((s) => s.setSolving);
  const setSolveResult = useModelStore((s) => s.setSolveResult);
  const addConsoleLine = useModelStore((s) => s.addConsoleLine);
  const clearConsole = useModelStore((s) => s.clearConsole);
  const eescode = useModelStore((s) => s.eescode);
  const initials = useModelStore((s) => s.initials);
  const filePath = useModelStore((s) => s.filePath);
  const setSol = useModelStore((s) => s.setSol);
  const setInitials = useModelStore((s) => s.setInitials);
  const setEescode = useModelStore((s) => s.setEescode);
  const setConf = useModelStore((s) => s.setConf);
  const setFilePath = useModelStore((s) => s.setFilePath);
  const loadFile = useModelStore((s) => s.loadFile);

  const theme = useUIStore((s) => s.theme);
  const toggleTheme = useUIStore((s) => s.toggleTheme);
  const bottomPanelOpen = useUIStore((s) => s.bottomPanelOpen);
  const toggleBottomPanel = useUIStore((s) => s.toggleBottomPanel);
  const setBottomPanelOpen = useUIStore((s) => s.setBottomPanelOpen);

  const [examples, setExamples] = useState<ExampleFile[]>([]);
  const [showExamples, setShowExamples] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const eventSourceRef = useRef<EventSource | null>(null);

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
                // Fetch .sol content
                api.getSol().then((sol) => setSol(sol.content)).catch(() => {});
              }
              setSolving(false);
              es.close();
              break;
            }
            case 'error': {
              const result = event.result as SolveResponse | undefined;
              if (result) {
                setSolveResult(result);
                addConsoleLine(`>>> Solve FAILED: ${result.status}`);
                if (result.errorMessage) addConsoleLine(`    ${result.errorMessage}`);
                if (result.detailedError) addConsoleLine(`    ${result.detailedError}`);
                for (const block of result.blockResults) {
                  if (!block.success) {
                    addConsoleLine(
                      `    Block ${block.id}: ${block.status} (res=${block.maxResidual.toExponential(2)}) ${block.errorMessage}`
                    );
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
    [eescode, initials, setSolving, clearConsole, addConsoleLine, setSolveResult, setSol, setBottomPanelOpen]
  );

  // Open file via backend path (for examples and known paths)
  const openFileByPath = useCallback(async (path: string) => {
    try {
      const res = await api.openFile(path);
      if (res.success) {
        const ees = await api.getEescode();
        const init = await api.getInitials();
        const sol = await api.getSol();
        const conf = await api.getConf();
        loadFile(res.filePath, ees.content, init.content, sol.content, conf.content);
        addConsoleLine(`>>> Opened: ${res.filePath}`);
      }
    } catch (err: any) {
      addConsoleLine(`>>> ERROR opening file: ${err.message}`);
    }
  }, [loadFile, addConsoleLine]);

  // Open file via browser file picker
  const handleOpenFile = useCallback(async () => {
    // Try File System Access API first (Chrome/Edge)
    if ('showOpenFilePicker' in window) {
      try {
        const [handle] = await (window as any).showOpenFilePicker({
          types: [{
            description: 'EES Code Files',
            accept: { 'text/plain': ['.eescode'] },
          }],
          multiple: false,
        });
        const file = await handle.getFile();
        const content = await file.text();
        
        // Load into editor and send to backend
        setEescode(content);
        await api.putEescode(content);
        
        // Try to load companion files from same directory
        try {
          const dirHandle = await handle.getParent?.();
          if (dirHandle) {
            const stem = file.name.replace(/\.eescode$/, '');
            // Try .initials
            try {
              const initHandle = await dirHandle.getFileHandle(stem + '.initials');
              const initFile = await initHandle.getFile();
              const initContent = await initFile.text();
              setInitials(initContent);
              await api.putInitials(initContent);
            } catch { /* no initials file */ }
            // Try coolsolve.conf
            try {
              const confHandle = await dirHandle.getFileHandle('coolsolve.conf');
              const confFile = await confHandle.getFile();
              const confContent = await confFile.text();
              setConf(confContent);
              await api.putConf(confContent);
            } catch { /* no conf file */ }
          }
        } catch { /* getParent not supported */ }
        
        loadFile(file.name, content, useModelStore.getState().initials, '', useModelStore.getState().conf);
        addConsoleLine(`>>> Opened: ${file.name}`);
        return;
      } catch (e: any) {
        if (e.name === 'AbortError') return; // User cancelled
        // Fall through to input element
      }
    }
    
    // Fallback: use hidden file input
    fileInputRef.current?.click();
  }, [loadFile, addConsoleLine, setEescode, setInitials, setConf]);

  // Handle file input change (fallback for browsers without File System Access API)
  const handleFileInputChange = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = e.target.files;
    if (!files || files.length === 0) return;

    let eesContent = '';
    let initContent = '';
    let confContent = '';
    let fileName = '';

    for (const file of Array.from(files)) {
      const content = await file.text();
      if (file.name.endsWith('.eescode')) {
        eesContent = content;
        fileName = file.name;
      } else if (file.name.endsWith('.initials')) {
        initContent = content;
      } else if (file.name === 'coolsolve.conf') {
        confContent = content;
      }
    }

    if (eesContent) {
      await api.putEescode(eesContent);
      if (initContent) await api.putInitials(initContent);
      if (confContent) await api.putConf(confContent);
      loadFile(fileName, eesContent, initContent, '', confContent);
      addConsoleLine(`>>> Opened: ${fileName}`);
    }

    // Reset input so same file can be selected again
    e.target.value = '';
  }, [loadFile, addConsoleLine]);

  // Save file
  const handleSave = useCallback(async () => {
    try {
      await api.putEescode(eescode);
      await api.saveFile();
      addConsoleLine(`>>> Saved: ${filePath || '<no path>'}`);
    } catch (err: any) {
      addConsoleLine(`>>> ERROR saving: ${err.message}`);
    }
  }, [eescode, filePath, addConsoleLine]);

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

  // Save As — prompt for path or use File System Access API
  const handleSaveAs = useCallback(async () => {
    // Try File System Access API (Chrome/Edge)
    if ('showSaveFilePicker' in window) {
      try {
        const handle = await (window as any).showSaveFilePicker({
          suggestedName: filePath ? filePath.split('/').pop() : 'model.eescode',
          types: [{
            description: 'EES Code Files',
            accept: { 'text/plain': ['.eescode'] },
          }],
        });
        const writable = await handle.createWritable();
        await writable.write(eescode);
        await writable.close();
        addConsoleLine(`>>> Saved: ${handle.name}`);
        setFilePath(handle.name);
        return;
      } catch (e: any) {
        if (e.name === 'AbortError') return;
        // Fall through to prompt
      }
    }
    // Fallback: prompt for path
    const newPath = prompt('Save as (full path):', filePath || 'model.eescode');
    if (!newPath) return;
    try {
      const res = await api.saveFileAs(newPath);
      addConsoleLine(`>>> Saved as: ${res.filePath}`);
      setFilePath(res.filePath);
    } catch (err: any) {
      addConsoleLine(`>>> ERROR: ${err.message}`);
    }
  }, [eescode, filePath, addConsoleLine, setFilePath]);

  // Stop/cancel solve
  const handleStopSolve = useCallback(async () => {
    try {
      await api.cancelSolve();
      addConsoleLine('>>> Cancel requested...');
    } catch (err: any) {
      addConsoleLine(`>>> ${err.message}`);
    }
  }, [addConsoleLine]);

  // Comment toggles
  const handleBraceComment = useCallback(() => {
    toggleBraceComment(editorInstance);
  }, []);

  const handleQuoteComment = useCallback(() => {
    toggleQuoteComment(editorInstance);
  }, []);

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
        if (e.shiftKey) {
          handleSaveAs();
        } else {
          handleSave();
        }
      } else if ((e.ctrlKey || e.metaKey) && e.key === 'o') {
        e.preventDefault();
        handleOpenFile();
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
  }, [solving, handleSolve, handleSave, handleSaveAs, handleOpenFile, handleUpdateGuesses, handleStopSolve]);

  return (
    <div className="toolbar">
      {/* Hidden file input for fallback browser file picker */}
      <input
        ref={fileInputRef}
        type="file"
        accept=".eescode,.initials,.conf"
        multiple
        style={{ display: 'none' }}
        onChange={handleFileInputChange}
      />

      {/* File group */}
      <div className="toolbar-group">
        <button className="toolbar-btn" onClick={handleOpenFile} title="Open file (Ctrl+O)">
          <FolderOpen size={16} /> Open
        </button>
        <button className="toolbar-btn" onClick={handleSave} title="Save (Ctrl+S)">
          <Save size={16} /> Save
        </button>
        <button className="toolbar-btn" onClick={handleSaveAs} title="Save As (Ctrl+Shift+S)">
          <SaveAll size={16} /> Save As
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

      <div className="toolbar-spacer" />

      {/* Status */}
      <div className="toolbar-status">
        {filePath ? filePath.split('/').pop() : 'Untitled'}
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
