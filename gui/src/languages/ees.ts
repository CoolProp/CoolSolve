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
      'nDay_', 'NumDay_', 'AveDay_', 'EqnTime_',
      'DELTAh_s', 'T_sat', 'P_sat', 'h_sat', 's_sat',
      'convert', 'unitsystem', 'error', 'warning', 'colebrook',
      // Lookup table functions
      'interpolate', 'interpolate1', 'interpolate2', 'interpolate2dm',
      'lookup', 'lookup$', 'lookupcol', 'lookupcol1', 'lookupcellempty',
      'tablevalue', 'nlookuprows', 'nlookupcolumns',
      'sumlookup', 'avglookup', 'maxlookup', 'minlookup', 'stddevlookup',
      // Equation-based integration (INTEGRAL / $IntegralTable)
      'integral', 'integralvalue',
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
        [/"/, 'string', '@quoteComment'],

        // Directives
        [/\$\w+/, 'keyword.directive'],

        // Numbers
        [/\d+\.?\d*([eE][+-]?\d+)?/, 'number'],
        [/\.\d+([eE][+-]?\d+)?/, 'number'],

        // Identifiers: letters (accented ones too), digits, _ and the symbols % | # \ ° of EES names,
        // with the $ that ends a string variable or a string function (LOOKUP$)
        [/[a-zA-Z_\u00C0-\u00FF][\w%|#\\\u00B0\u00B2\u00B3\u00B5\u00B9\u00C0-\u00FF]*\$?/, {
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

      quoteComment: [
        [/[^"]+/, 'string'],
        [/"/, 'string', '@pop'],
      ],

      array: [
        [/\d+/, 'number'],
        [/[a-zA-Z_\u00C0-\u00FF][\w%|#\\\u00B0\u00B2\u00B3\u00B5\u00B9\u00C0-\u00FF]*/, 'identifier'],
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
