import { useState, useEffect, useCallback } from 'react';
import { api } from '../api/client';
import { FileText, RefreshCw } from 'lucide-react';

interface DebugFile {
  name: string;
  size: number;
}

export default function DebugViewer() {
  const [files, setFiles] = useState<DebugFile[]>([]);
  const [selected, setSelected] = useState<string | null>(null);
  const [content, setContent] = useState('');
  const [loading, setLoading] = useState(false);

  const refresh = useCallback(async () => {
    try {
      const res = await api.getDebugFiles();
      setFiles(res.files);
    } catch {
      setFiles([]);
    }
  }, []);

  useEffect(() => {
    refresh();
  }, [refresh]);

  const handleSelectFile = useCallback(async (name: string) => {
    setSelected(name);
    setLoading(true);
    try {
      const res = await api.getDebugFile(name);
      setContent(res.content);
    } catch (err: any) {
      setContent(`Error loading file: ${err.message}`);
    }
    setLoading(false);
  }, []);

  if (files.length === 0) {
    return (
      <div className="debug-viewer">
        <div className="debug-empty">
          <p>No debug output available.</p>
          <p>Click <strong>Debug</strong> (Ctrl+Shift+Enter) to run a solve with tracing enabled.</p>
          <button className="toolbar-btn" onClick={refresh} style={{ marginTop: 8 }}>
            <RefreshCw size={14} /> Refresh
          </button>
        </div>
      </div>
    );
  }

  // Priority order for file display
  const fileOrder = ['report.md', 'solver_errors.md', 'equations.md', 'variables.md', 'variables.csv', 'analysis.json', 'residuals.md', 'ees_residuals.txt', 'original.eescode', 'coolsolve.conf'];
  const sortedFiles = [...files].sort((a, b) => {
    const ai = fileOrder.indexOf(a.name);
    const bi = fileOrder.indexOf(b.name);
    if (ai >= 0 && bi >= 0) return ai - bi;
    if (ai >= 0) return -1;
    if (bi >= 0) return 1;
    return a.name.localeCompare(b.name);
  });

  return (
    <div className="debug-viewer">
      <div className="debug-file-list">
        <div className="debug-file-list-header">
          <span>Debug Files</span>
          <button className="toolbar-btn" onClick={refresh} title="Refresh" style={{ padding: '2px 4px' }}>
            <RefreshCw size={12} />
          </button>
        </div>
        {sortedFiles.map((f) => (
          <button
            key={f.name}
            className={`debug-file-item ${selected === f.name ? 'active' : ''}`}
            onClick={() => handleSelectFile(f.name)}
          >
            <FileText size={12} />
            <span>{f.name}</span>
            <span className="debug-file-size">{formatSize(f.size)}</span>
          </button>
        ))}
      </div>
      <div className="debug-content">
        {loading ? (
          <div className="debug-loading"><span className="spinner" /></div>
        ) : selected ? (
          <pre className="debug-pre">{content}</pre>
        ) : (
          <div className="debug-empty">Select a file to view its content.</div>
        )}
      </div>
    </div>
  );
}

function formatSize(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}
