import { useEffect, useRef } from 'react';
import { useModelStore } from '../stores/modelStore';

export default function Console() {
  const consoleLines = useModelStore((s) => s.consoleLines);
  const containerRef = useRef<HTMLDivElement>(null);

  // Auto-scroll to bottom
  useEffect(() => {
    if (containerRef.current) {
      containerRef.current.scrollTop = containerRef.current.scrollHeight;
    }
  }, [consoleLines]);

  const getLineClass = (line: string) => {
    if (line.includes('ERROR') || line.includes('FAIL') || line.includes('Error') || line.includes('[Error]'))
      return 'console-line error';
    if (line.includes('[Warning]'))
      return 'console-line warning';
    if (line.includes('SUCCESS') || line.includes('converged'))
      return 'console-line success';
    if (line.startsWith('>>>') || line.startsWith('---'))
      return 'console-line info';
    return 'console-line';
  };

  return (
    <div className="console" ref={containerRef}>
      {consoleLines.length === 0 ? (
        <div className="console-line" style={{ color: 'var(--text-muted)' }}>
          Ready. Open a file or paste EES code to begin.
        </div>
      ) : (
        consoleLines.map((line, i) => (
          <div key={i} className={getLineClass(line)}>
            {line}
          </div>
        ))
      )}
    </div>
  );
}
