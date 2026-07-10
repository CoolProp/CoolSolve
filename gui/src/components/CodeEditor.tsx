import { useCallback, useRef, useEffect } from 'react';
import type * as Monaco from 'monaco-editor';
import Editor, { OnMount } from '@monaco-editor/react';
import { registerEESLanguage, EES_LANGUAGE_ID } from '../languages/ees';
import { useModelStore } from '../stores/modelStore';
import { useUIStore } from '../stores/uiStore';
import { api } from '../api/client';
import { setEditorInstance, toggleBraceComment, toggleQuoteComment } from './editorUtils';

let languageRegistered = false;

export default function CodeEditor() {
  const eescode = useModelStore((s) => s.eescode);
  const setEescode = useModelStore((s) => s.setEescode);
  const setParseResult = useModelStore((s) => s.setParseResult);
  const theme = useUIStore((s) => s.theme);
  const editorRef = useRef<Monaco.editor.IStandaloneCodeEditor | null>(null);
  const parseTimeoutRef = useRef<ReturnType<typeof setTimeout>>();

  const handleMount: OnMount = (editor, monaco) => {
    editorRef.current = editor;
    setEditorInstance(editor);

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
            const monaco = window.monaco;
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
