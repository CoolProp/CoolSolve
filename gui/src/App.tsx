import { useEffect, useState } from 'react';
import Toolbar from './components/Toolbar';
import CodeEditor from './components/CodeEditor';
import VariableTable from './components/VariableTable';
import ArrayTable from './components/ArrayTable';
import ConfigEditor from './components/ConfigEditor';
import DebugViewer from './components/DebugViewer';
import Console from './components/Console';
import SplitPane from './components/SplitPane';
import { useModelStore } from './stores/modelStore';
import { useUIStore } from './stores/uiStore';
import { api } from './api/client';
import './App.css';

export default function App() {
  const theme = useUIStore((s) => s.theme);
  const rightTab = useUIStore((s) => s.rightTab);
  const setRightTab = useUIStore((s) => s.setRightTab);
  const bottomPanelOpen = useUIStore((s) => s.bottomPanelOpen);
  const toggleBottomPanel = useUIStore((s) => s.toggleBottomPanel);
  const equationCount = useModelStore((s) => s.equationCount);
  const variableCount = useModelStore((s) => s.variableCount);
  const isSquare = useModelStore((s) => s.isSquare);
  const parseErrors = useModelStore((s) => s.parseErrors);
  const lastResult = useModelStore((s) => s.lastResult);
  const solving = useModelStore((s) => s.solving);

  const [ready, setReady] = useState(false);

  // Initialize: warmup CoolProp
  useEffect(() => {
    api.warmup().then(() => setReady(true)).catch(() => setReady(true));
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
            </div>
            <div className="tab-content">
              {rightTab === 'variables' && <VariableTable />}
              {rightTab === 'arrays' && <ArrayTable />}
              {rightTab === 'config' && <ConfigEditor />}
              {rightTab === 'debug' && <DebugViewer />}
            </div>
          </div>
        </SplitPane>

        {/* Bottom panel: Console */}
        <div className="bottom-panel">
          <div className="tab-bar bottom-panel-header" onClick={toggleBottomPanel}>
            <button className="tab-btn active">Console</button>
          </div>
          {bottomPanelOpen && <Console />}
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
          <span style={{ color: '#f87171' }}>{parseErrors.length} error(s)</span>
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
