import { useEffect, useState } from 'react';
import Toolbar from './components/Toolbar';
import CodeEditor from './components/CodeEditor';
import VariableTable from './components/VariableTable';
import ArrayTable from './components/ArrayTable';
import ConfigEditor from './components/ConfigEditor';
import DebugViewer from './components/DebugViewer';
import ThermoDiagram from './components/ThermoDiagram';
import ParametricStudy from './components/ParametricStudy';
import Console from './components/Console';
import LookupTableEditor from './components/LookupTableEditor';
import IntegralTable from './components/IntegralTable';
import UsageStats from './components/UsageStats';
import SplitPane from './components/SplitPane';
import { useModelStore } from './stores/modelStore';
import { useUIStore } from './stores/uiStore';
import { api } from './api/client';
import { enableTabSession, getTabSessionId } from './api/session';
import { openLibraryModel } from './utils/library';
import { showSessionModel } from './utils/sessionModel';
import './App.css';

// Hidden utility page (solve usage statistics) — reachable at /stats only,
// deliberately not linked from anywhere in the UI.
const STATS_PATH = '/stats';

// An open tab tells the server that its session is in use at this interval: the server drops
// sessions idle for 24 h, which would otherwise lose the model of a tab left open overnight.
const KEEPALIVE_INTERVAL_MS = 15 * 60 * 1000;

export default function App() {
  if (window.location.pathname === STATS_PATH) {
    return <UsageStats />;
  }
  return <MainApp />;
}

function MainApp() {
  const theme = useUIStore((s) => s.theme);
  const rightTab = useUIStore((s) => s.rightTab);
  const setRightTab = useUIStore((s) => s.setRightTab);
  const bottomPanelOpen = useUIStore((s) => s.bottomPanelOpen);
  const toggleBottomPanel = useUIStore((s) => s.toggleBottomPanel);
  const bottomTab = useUIStore((s) => s.bottomTab);
  const setBottomTab = useUIStore((s) => s.setBottomTab);
  const equationCount = useModelStore((s) => s.equationCount);
  const variableCount = useModelStore((s) => s.variableCount);
  const isSquare = useModelStore((s) => s.isSquare);
  const parseErrors = useModelStore((s) => s.parseErrors);
  const lastResult = useModelStore((s) => s.lastResult);
  const solving = useModelStore((s) => s.solving);

  const [ready, setReady] = useState(false);

  // Initialize: warmup CoolProp, then open what was asked for: a library model
  // ("Open in CoolSolve" link: /?library=CSL-0012) or the file passed on the command line
  useEffect(() => {
    const params = new URLSearchParams(window.location.search);
    const libraryParam = params.get('library');
    // A library link opens in its own instance, so it never replaces the model of another tab
    if (libraryParam) enableTabSession();
    api.warmup()
      .catch(() => {})
      .finally(async () => {
        try {
          if (libraryParam) {
            await openLibraryModel(libraryParam, params.get('file') ?? undefined, false);
            // Consumed: a reload shows the tab's session (below) instead of reopening the model
            window.history.replaceState(null, '', window.location.pathname);
          } else if (getTabSessionId()) {
            // Reload of a tab that has its own session: the model it holds, edits included
            await showSessionModel();
          } else {
            const initial = await api.getInitialFile();
            if (initial?.path) {
              await api.openFile(initial.path);
              await showSessionModel();
            }
          }
        } catch {
          // Not critical — user can open files manually
        }
        setReady(true);
      });
  }, []);

  useEffect(() => {
    const timer = window.setInterval(() => { api.keepAlive().catch(() => {}); }, KEEPALIVE_INTERVAL_MS);
    return () => window.clearInterval(timer);
  }, []);

  // Apply theme
  useEffect(() => {
    document.documentElement.setAttribute('data-theme', theme);
  }, [theme]);

  if (!ready) {
    return (
      <div className="loading-splash" data-theme={theme}>
        <h1>CoolSolve</h1>
        <div className="spinner" style={{ width: 32, height: 32 }} />
        <p>Initializing thermodynamic engine...</p>
      </div>
    );
  }

  return (
    <div className="app">
      {/* Toolbar */}
      <Toolbar />

      {/* Main content: split vertically (top: editor+panel, bottom: console) */}
      <SplitPane
        direction="vertical"
        defaultSize={bottomPanelOpen ? 70 : 96}
        minSize={28}
        className="main-vertical-split"
      >
        {/* Top half: editor + right panel, split horizontally */}
        <SplitPane
          direction="horizontal"
          defaultSize={60}
          minSize={300}
        >
          {/* Left: Code editor */}
          <CodeEditor />

          {/* Right: Panel with tabs */}
          <div className="right-panel">
            <div className="tab-bar">
              <button
                className={`tab-btn ${rightTab === 'variables' ? 'active' : ''}`}
                onClick={() => setRightTab('variables')}
              >
                Variables
              </button>
              <button
                className={`tab-btn ${rightTab === 'arrays' ? 'active' : ''}`}
                onClick={() => setRightTab('arrays')}
              >
                Arrays
              </button>
              <button
                className={`tab-btn ${rightTab === 'config' ? 'active' : ''}`}
                onClick={() => setRightTab('config')}
              >
                Config
              </button>
              <button
                className={`tab-btn ${rightTab === 'debug' ? 'active' : ''}`}
                onClick={() => setRightTab('debug')}
              >
                Debug
              </button>
              <button
                className={`tab-btn ${rightTab === 'diagrams' ? 'active' : ''}`}
                onClick={() => setRightTab('diagrams')}
              >
                Diagrams
              </button>
            </div>
            <div className="tab-content">
              {rightTab === 'variables' && <VariableTable />}
              {rightTab === 'arrays' && <ArrayTable />}
              {rightTab === 'config' && <ConfigEditor />}
              {rightTab === 'debug' && <DebugViewer />}
              <div style={{ display: rightTab === 'diagrams' ? 'contents' : 'none' }}>
                <ThermoDiagram />
              </div>
            </div>
          </div>
        </SplitPane>

        {/* Bottom panel: Console + Parametric Study */}
        <div className="bottom-panel">
          <div className="tab-bar bottom-panel-header">
            <button
              className={`tab-btn ${bottomTab === 'console' ? 'active' : ''}`}
              onClick={() => { if (bottomTab === 'console' && bottomPanelOpen) toggleBottomPanel(); else setBottomTab('console'); }}
            >
              Console
            </button>
            <button
              className={`tab-btn ${bottomTab === 'sensitivity' ? 'active' : ''}`}
              onClick={() => { if (bottomTab === 'sensitivity' && bottomPanelOpen) toggleBottomPanel(); else setBottomTab('sensitivity'); }}
            >
              Parametric
            </button>
            <button
              className={`tab-btn ${bottomTab === 'lookuptables' ? 'active' : ''}`}
              onClick={() => { if (bottomTab === 'lookuptables' && bottomPanelOpen) toggleBottomPanel(); else setBottomTab('lookuptables'); }}
            >
              Lookup Tables
            </button>
            <button
              className={`tab-btn ${bottomTab === 'integral' ? 'active' : ''}`}
              onClick={() => { if (bottomTab === 'integral' && bottomPanelOpen) toggleBottomPanel(); else setBottomTab('integral'); }}
            >
              Integral
            </button>
          </div>
          {bottomPanelOpen && (
            <>
              {bottomTab === 'console' && <Console />}
              {bottomTab === 'sensitivity' && <ParametricStudy />}
              {bottomTab === 'lookuptables' && <LookupTableEditor />}
              {bottomTab === 'integral' && <IntegralTable />}
            </>
          )}
        </div>
      </SplitPane>

      {/* Status bar */}
      <div className="status-bar">
        <span>
          {equationCount > 0
            ? `${equationCount} equations, ${variableCount} variables`
            : 'No model loaded'}
        </span>
        {equationCount > 0 && (
          <span style={{ color: isSquare ? '#4ade80' : '#f87171' }}>
            {isSquare ? 'Square system' : 'Non-square system'}
          </span>
        )}
        {parseErrors.length > 0 && (
          <span
            style={{ color: '#f87171', cursor: 'help' }}
            title={parseErrors.map(e => `Line ${e.line}: ${e.message}`).join('\n')}
          >
            {parseErrors.length} error(s): Line {parseErrors[0].line}: {parseErrors[0].message}
          </span>
        )}
        {solving && <span style={{ color: '#facc15' }}>Solving...</span>}
        {lastResult && !solving && (
          <span style={{ color: lastResult.success ? '#4ade80' : '#f87171' }}>
            {lastResult.success ? 'Solved' : 'Failed'}
          </span>
        )}
      </div>
    </div>
  );
}
