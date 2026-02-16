import { useCallback, useRef, useEffect } from 'react';
import Editor, { OnMount } from '@monaco-editor/react';
import { registerEESLanguage, EES_LANGUAGE_ID } from '../languages/ees';
import { useModelStore } from '../stores/modelStore';
import { useUIStore } from '../stores/uiStore';
import { api } from '../api/client';

let languageRegistered = false;

export default function CodeEditor() {
  const eescode = useModelStore((s) => s.eescode);
  const setEescode = useModelStore((s) => s.setEescode);
  const setParseResult = useModelStore((s) => s.setParseResult);
  const theme = useUIStore((s) => s.theme);
  const editorRef = useRef<any>(null);
  const parseTimeoutRef = useRef<ReturnType<typeof setTimeout>>();

  const handleMount: OnMount = (editor, monaco) => {
    editorRef.current = editor;

    if (!languageRegistered) {
      registerEESLanguage(monaco);
      languageRegistered = true;
    }

    // Set model language
    const model = editor.getModel();
    if (model) {
      monaco.editor.setModelLanguage(model, EES_LANGUAGE_ID);
    }

    // Focus editor
    editor.focus();
  };

  // Debounced parse on change
  const handleChange = useCallback(
    (value: string | undefined) => {
      const code = value || '';
      setEescode(code);

      // Debounced parse
      if (parseTimeoutRef.current) clearTimeout(parseTimeoutRef.current);
      parseTimeoutRef.current = setTimeout(async () => {
        if (!code.trim()) {
          setParseResult([], 0, 0, true);
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

          // Set Monaco markers for errors
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
                monaco.editor.setModelMarkers(model, 'ees', markers);
              }
            }
          }
        } catch {
          // Parse endpoint unavailable, no big deal
        }
      }, 500);
    },
    [setEescode, setParseResult]
  );

  // Initial parse when eescode changes externally (e.g., file load)
  useEffect(() => {
    if (eescode) {
      handleChange(eescode);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

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
