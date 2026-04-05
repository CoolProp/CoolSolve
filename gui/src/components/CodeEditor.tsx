import { useCallback, useRef, useEffect } from 'react';
import Editor, { OnMount } from '@monaco-editor/react';
import { registerEESLanguage, EES_LANGUAGE_ID } from '../languages/ees';
import { useModelStore } from '../stores/modelStore';
import { useUIStore } from '../stores/uiStore';
import { api } from '../api/client';

let languageRegistered = false;

// Expose editor ref globally so Toolbar can invoke comment actions
export let editorInstance: any = null;

/**
 * Toggle EES invisible comment: wraps/unwraps selected lines in { ... }
 */
export function toggleBraceComment(editor: any) {
  if (!editor) return;
  const model = editor.getModel();
  const selection = editor.getSelection();
  if (!model || !selection) return;

  const startLine = selection.startLineNumber;
  const endLine = selection.endLineNumber;
  const firstLineText = model.getLineContent(startLine);
  const lastLineText = model.getLineContent(endLine);

  // Check if already wrapped: first line starts with { and last line ends with }
  const trimFirst = firstLineText.trimStart();
  const trimLast = lastLineText.trimEnd();
  const isWrapped = trimFirst.startsWith('{') && trimLast.endsWith('}');

  const edits: any[] = [];
  if (isWrapped) {
    // Unwrap: remove leading { from the first line and trailing } from the last line
    if (startLine === endLine) {
      // Single line: remove first { and last }
      const openIdx = firstLineText.indexOf('{');
      const closeIdx = lastLineText.lastIndexOf('}');
      if (openIdx >= 0 && closeIdx > openIdx) {
        const newText = firstLineText.substring(0, openIdx) + firstLineText.substring(openIdx + 1, closeIdx) + firstLineText.substring(closeIdx + 1);
        edits.push({ range: { startLineNumber: startLine, startColumn: 1, endLineNumber: startLine, endColumn: firstLineText.length + 1 }, text: newText });
      }
    } else {
      // Multi-line: remove { from first line, } from last line
      const openIdx = firstLineText.indexOf('{');
      if (openIdx >= 0) {
        const newFirst = firstLineText.substring(0, openIdx) + firstLineText.substring(openIdx + 1);
        edits.push({ range: { startLineNumber: startLine, startColumn: 1, endLineNumber: startLine, endColumn: firstLineText.length + 1 }, text: newFirst });
      }
      const closeIdx = lastLineText.lastIndexOf('}');
      if (closeIdx >= 0) {
        const newLast = lastLineText.substring(0, closeIdx) + lastLineText.substring(closeIdx + 1);
        edits.push({ range: { startLineNumber: endLine, startColumn: 1, endLineNumber: endLine, endColumn: lastLineText.length + 1 }, text: newLast });
      }
    }
  } else {
    // Wrap: add { at start of first line, } at end of last line
    edits.push({ range: { startLineNumber: startLine, startColumn: 1, endLineNumber: startLine, endColumn: 1 }, text: '{' });
    const lastLen = model.getLineContent(endLine).length;
    edits.push({ range: { startLineNumber: endLine, startColumn: lastLen + 1, endLineNumber: endLine, endColumn: lastLen + 1 }, text: '}' });
  }

  if (edits.length > 0) {
    editor.executeEdits('comment-toggle', edits);
  }
}

/**
 * Toggle EES visible comment: wraps/unwraps selected text in "..."
 */
export function toggleQuoteComment(editor: any) {
  if (!editor) return;
  const model = editor.getModel();
  const selection = editor.getSelection();
  if (!model || !selection) return;

  const selectedText = model.getValueInRange(selection);
  if (!selectedText) return;

  const trimmed = selectedText.trim();
  const isWrapped = trimmed.startsWith('"') && trimmed.endsWith('"') && trimmed.length >= 2;

  if (isWrapped) {
    // Unwrap: remove surrounding quotes
    const newText = selectedText.replace(/^(\s*)"/, '$1').replace(/"(\s*)$/, '$1');
    editor.executeEdits('comment-toggle', [{
      range: selection,
      text: newText,
    }]);
  } else {
    // Wrap: add quotes around selected text
    editor.executeEdits('comment-toggle', [{
      range: selection,
      text: '"' + selectedText + '"',
    }]);
  }
}

export default function CodeEditor() {
  const eescode = useModelStore((s) => s.eescode);
  const setEescode = useModelStore((s) => s.setEescode);
  const setParseResult = useModelStore((s) => s.setParseResult);
  const theme = useUIStore((s) => s.theme);
  const editorRef = useRef<any>(null);
  const parseTimeoutRef = useRef<ReturnType<typeof setTimeout>>();

  const handleMount: OnMount = (editor, monaco) => {
    editorRef.current = editor;
    editorInstance = editor;

    if (!languageRegistered) {
      registerEESLanguage(monaco);
      languageRegistered = true;
    }

    // Set model language
    const model = editor.getModel();
    if (model) {
      monaco.editor.setModelLanguage(model, EES_LANGUAGE_ID);
    }

    // Register comment toggle actions
    editor.addAction({
      id: 'ees.toggleBraceComment',
      label: 'Toggle { } Comment',
      keybindings: [monaco.KeyMod.CtrlCmd | monaco.KeyCode.Slash],
      run: () => toggleBraceComment(editor),
    });
    editor.addAction({
      id: 'ees.toggleQuoteComment',
      label: 'Toggle " " Comment',
      keybindings: [monaco.KeyMod.CtrlCmd | monaco.KeyMod.Shift | monaco.KeyCode.Slash],
      run: () => toggleQuoteComment(editor),
    });

    // Focus editor
    editor.focus();
  };

  // Debounced parse (does not call setEescode — call separately when needed)
  const doParse = useCallback(
    (code: string) => {
      if (parseTimeoutRef.current) clearTimeout(parseTimeoutRef.current);
      parseTimeoutRef.current = setTimeout(async () => {
        if (!code.trim()) {
          setParseResult([], 0, 0, true);
          useModelStore.getState().setParsedVariables([]);
          return;
        }
        try {
          const result = await api.parse(code);
          setParseResult(
            result.errors,
            result.equationCount,
            result.variableCount || 0,
            result.isSquare ?? true
          );

          // Store parsed variable info (units, imposed flags) for the VariableTable
          if (result.variables) {
            useModelStore.getState().setParsedVariables(result.variables);
          }

          // Set Monaco markers for errors and diagnostics
          if (editorRef.current) {
            const monaco = (window as any).monaco;
            if (monaco) {
              const model = editorRef.current.getModel();
              if (model) {
                const markers = result.errors.map((err) => ({
                  severity: monaco.MarkerSeverity.Error,
                  startLineNumber: err.line,
                  startColumn: err.column || 1,
                  endLineNumber: err.line,
                  endColumn: (err.column || 1) + 10,
                  message: err.message,
                }));
                // Add diagnostic markers (warnings/info from P002, P003, etc.)
                if (result.diagnostics) {
                  for (const d of result.diagnostics) {
                    const sev = d.severity === 'Error' ? monaco.MarkerSeverity.Error
                              : d.severity === 'Warning' ? monaco.MarkerSeverity.Warning
                              : monaco.MarkerSeverity.Info;
                    markers.push({
                      severity: sev,
                      startLineNumber: d.line || 1,
                      startColumn: d.column || 1,
                      endLineNumber: d.line || 1,
                      endColumn: (d.column || 1) + 10,
                      message: `[${d.code}] ${d.message}`,
                    });
                  }
                }
                monaco.editor.setModelMarkers(model, 'ees', markers);
              }
            }
          }
        } catch {
          // Parse endpoint unavailable, no big deal
        }
      }, 500);
    },
    [setParseResult]
  );

  // Monaco onChange handler — user typing
  const handleChange = useCallback(
    (value: string | undefined) => {
      const code = value || '';
      setEescode(code);
      doParse(code);
    },
    [setEescode, doParse]
  );

  // Re-parse when eescode changes externally (file load, examples, back button, etc.)
  // This also handles the initial mount.  The debounce in doParse prevents double-work
  // when the change originates from handleChange (same code → same debounced timer).
  useEffect(() => {
    if (eescode) {
      doParse(eescode);
    }
  }, [eescode, doParse]);

  return (
    <div className="editor-pane">
      <Editor
        height="100%"
        defaultLanguage={EES_LANGUAGE_ID}
        value={eescode}
        onChange={handleChange}
        onMount={handleMount}
        theme={theme === 'dark' ? 'vs-dark' : 'vs'}
        options={{
          fontSize: 13,
          minimap: { enabled: true },
          lineNumbers: 'on',
          wordWrap: 'on',
          automaticLayout: true,
          scrollBeyondLastLine: false,
          renderWhitespace: 'selection',
          bracketPairColorization: { enabled: true },
          tabSize: 4,
        }}
      />
    </div>
  );
}
