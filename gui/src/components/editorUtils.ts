/**
 * editorUtils.ts — Monaco editor helpers shared between the CodeEditor
 * component and the Toolbar.
 *
 * Kept in a sibling non-component module so that `CodeEditor.tsx` only
 * exports the React component (required by react-refresh / Fast Refresh
 * for clean hot-module replacement).
 */

import type * as Monaco from 'monaco-editor';

/**
 * `@monaco-editor/react`'s loader exposes the Monaco namespace on
 * `window`. Augment the DOM lib so we can read it without an `any` cast.
 */
declare global {
  interface Window {
    monaco?: typeof Monaco;
  }
}

/** Live reference to the currently mounted editor instance (or null). */
export let editorInstance: Monaco.editor.IStandaloneCodeEditor | null = null;

/** Bind the shared editor reference. Called from the editor's onMount. */
export function setEditorInstance(editor: Monaco.editor.IStandaloneCodeEditor | null): void {
  editorInstance = editor;
}

/**
 * Toggle EES invisible comment: wraps/unwraps selected lines in { ... }
 */
export function toggleBraceComment(editor: Monaco.editor.IStandaloneCodeEditor | null): void {
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

  const edits: Monaco.editor.IIdentifiedSingleEditOperation[] = [];
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
export function toggleQuoteComment(editor: Monaco.editor.IStandaloneCodeEditor | null): void {
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
