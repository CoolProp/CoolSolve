import type * as Monaco from 'monaco-editor';

export const EES_LANGUAGE_ID = 'ees';

export function registerEESLanguage(monaco: typeof Monaco) {
  // Register language
  monaco.languages.register({ id: EES_LANGUAGE_ID });

  // Token provider (Monarch)
  monaco.languages.setMonarchTokensProvider(EES_LANGUAGE_ID, {
    ignoreCase: true,
    keywords: [
      'FUNCTION', 'PROCEDURE', 'END', 'CALL',
      'IF', 'THEN', 'ELSE', 'ENDIF',
      'AND', 'OR', 'NOT',
      'REPEAT', 'UNTIL',
      'GOTO',
      'SUBPROGRAM', 'MODULE',
    ],
    typeKeywords: [
      'Water', 'R134a', 'R245fa', 'R22', 'R410A', 'R32',
      'CO2', 'Nitrogen', 'Oxygen', 'Air', 'Air_ha',
      'Ammonia', 'Propane', 'Butane', 'Isobutane',
      'R1234yf', 'R1234ze', 'R407C', 'R404A', 'R507A',
    ],
    builtinFunctions: [
      'enthalpy', 'entropy', 'density', 'pressure', 'temperature',
      'quality', 'cp', 'cv', 'specheat', 'volume',
      'viscosity', 'conductivity', 'prandtl', 'surfacetension',
      'soundspeed', 'internalenergy',
      'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'atan2',
      'exp', 'ln', 'log', 'log10', 'sqrt', 'abs',
      'min', 'max', 'round', 'trunc', 'ceil', 'floor',
      'DELTAh_s', 'T_sat', 'P_sat', 'h_sat', 's_sat',
      'convert',
    ],
    operators: [
      '=', ':=', '+', '-', '*', '/', '^',
      '<', '>', '<=', '>=', '<>', '==',
    ],
    symbols: /[=><!~?:&|+\-*/^%]+/,

    tokenizer: {
      root: [
        // Comments
        [/\/\/.*$/, 'comment'],
        [/\{/, 'comment', '@bracketComment'],
        [/"[^"]*"/, 'string'],

        // Directives
        [/\$\w+/, 'keyword.directive'],

        // Numbers
        [/\d+\.?\d*([eE][+-]?\d+)?/, 'number'],
        [/\.\d+([eE][+-]?\d+)?/, 'number'],

        // Identifiers
        [/[a-zA-Z_]\w*/, {
          cases: {
            '@keywords': 'keyword',
            '@typeKeywords': 'type.identifier',
            '@builtinFunctions': 'predefined',
            '@default': 'identifier',
          },
        }],

        // Array subscripts
        [/\[/, 'delimiter.bracket', '@array'],

        // Operators
        [/@symbols/, {
          cases: {
            '@operators': 'operator',
            '@default': '',
          },
        }],

        // Whitespace
        [/\s+/, 'white'],

        // Delimiters
        [/[(),:;]/, 'delimiter'],
      ],

      bracketComment: [
        [/[^}]+/, 'comment'],
        [/\}/, 'comment', '@pop'],
      ],

      array: [
        [/\d+/, 'number'],
        [/[a-zA-Z_]\w*/, 'identifier'],
        [/\]/, 'delimiter.bracket', '@pop'],
        [/[+\-*/]/, 'operator'],
      ],
    },
  });

  // Language configuration
  monaco.languages.setLanguageConfiguration(EES_LANGUAGE_ID, {
    comments: {
      lineComment: '//',
      blockComment: ['{', '}'],
    },
    brackets: [
      ['{', '}'],
      ['[', ']'],
      ['(', ')'],
    ],
    autoClosingPairs: [
      { open: '{', close: '}' },
      { open: '[', close: ']' },
      { open: '(', close: ')' },
      { open: '"', close: '"' },
    ],
    surroundingPairs: [
      { open: '{', close: '}' },
      { open: '[', close: ']' },
      { open: '(', close: ')' },
      { open: '"', close: '"' },
    ],
    folding: {
      markers: {
        start: /^\s*(FUNCTION|PROCEDURE|IF)\b/i,
        end: /^\s*(END|ENDIF)\b/i,
      },
    },
    indentationRules: {
      increaseIndentPattern: /^\s*(FUNCTION|PROCEDURE|IF|ELSE)\b/i,
      decreaseIndentPattern: /^\s*(END|ENDIF|ELSE)\b/i,
    },
  });
}
