import { BookOpen, ExternalLink, FileText } from 'lucide-react';
import { useModelStore } from '../stores/modelStore';
import { includeWindowUrl, noteIncludeOpener } from '../utils/includes';
import type { IncludeInfo } from '../api/types';

function tooltip(inc: IncludeInfo): string {
  const lines = [`$INCLUDE ${inc.name} (line ${inc.line}): ${inc.description}`];
  lines.push(inc.routines.length ? `Imports ${inc.routines.join(', ')}` : 'Defines no function or procedure');
  if (inc.nested.length) lines.push(`Includes in turn: ${inc.nested.join(', ')}`);
  lines.push('Click to open it in a new CoolSolve window');
  return lines.join('\n');
}

/**
 * The files the model includes, above the editor. Rendered only when the model has an `$INCLUDE`
 * (otherwise nothing changes); each file opens in a new CoolSolve window.
 */
export default function IncludeBar() {
  const includes = useModelStore((s) => s.includes);
  if (includes.length === 0) return null;
  return (
    <div className="include-bar">
      <span className="include-bar-label">Includes</span>
      {includes.map((inc) => (
        <a
          key={`${inc.line}:${inc.name}`}
          className="include-chip"
          title={tooltip(inc)}
          href={includeWindowUrl(inc.name)}
          target="_blank"
          rel="noopener"
          onClick={() => noteIncludeOpener(inc.name)}
          onAuxClick={() => noteIncludeOpener(inc.name)}
        >
          {inc.origin === 'library' ? <BookOpen size={12} /> : <FileText size={12} />}
          <span>{inc.origin === 'library' ? inc.name : inc.file}</span>
          {inc.origin === 'library' && <span className="include-chip-id">{inc.libraryId}</span>}
          <ExternalLink size={11} />
        </a>
      ))}
    </div>
  );
}
