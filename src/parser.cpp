#include "coolsolve/parser.h"
#include "coolsolve/conditional_directives.h"
#include "coolsolve/constants.h"
#include "coolsolve/diagnostic.h"
#include "coolsolve/ees_string.h"
#include "coolsolve/ir.h"
#include "coolsolve/units.h"
#include "coolsolve/user_hints.h"
#include <peglib.h>
#include <fstream>
#include <functional>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <set>
#include <unordered_set>

namespace coolsolve {

// ============================================================================
// EES Thermophysical Functions
// ============================================================================
// In EES, thermophysical functions have a fluid name as their first positional
// argument, which can be specified without quotes. This list is used to identify
// such functions so that unquoted first arguments are treated as string literals.

static const std::set<std::string> EES_THERMOPHYSICAL_FUNCTIONS = {
    // Standard thermodynamic properties
    "enthalpy", "entropy", "temperature", "pressure", "density",
    "quality", "volume", "cp", "cv", "viscosity", "conductivity",
    "t_sat", "p_sat", "t_crit", "p_crit", "v_crit", "t_triple",
    "molarmass", "soundspeed", "specheat", "intenergy",
    // Psychrometric / humid air functions
    "humrat", "wetbulb", "dewpoint", "relhum",
    // Additional EES functions
    "surfacetension", "prandtl", "thermaldiffusivity", "kinematicviscosity",
    "compressibilityfactor", "fugacity", "isentropicexponent", "volexpcoef",
    "isothermalcompress", "joulethomsoncoef", "dipole", "acentricfactor",
    "enthalpy_fusion", "enthalpy_vaporization", "normalboilingpt", "freezingpt",
    "gibbsfreeenergy", "helmholtzfreeenergy", "phase$",
    "higherheatingvalue", "lowerheatingvalue", "enthalpy_formation",
    "massfraction", "molarmass_soln", "saturation%",
    "p_melting", "p_sublimation", "t_melting"
};

// Check if a function name is an EES thermophysical function (case-insensitive)
static bool isEESThermophysicalFunction(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return EES_THERMOPHYSICAL_FUNCTIONS.find(lower) != EES_THERMOPHYSICAL_FUNCTIONS.end();
}

// ============================================================================
// Characters of variable names (EES help, "Variable names")
// ============================================================================
// A name starts with a letter or an underscore and holds letters, digits, `_` and the
// symbols `%`, `|` and `#` (`C%`, `K|star_r`, `Nu#_f`); `\` may join two parts of a name
// (`dU\dtau`), `$` ends a string variable. Letters are ASCII or the letters of Latin-1 in UTF-8
// (`réf`), and `°`, `µ`, `²`, `³`, `¹` may follow (`n°_90°elbow`), see ees_string.h.

/// True when a name can start with the character at s[i]: a letter or an underscore.
static bool nameStartsAt(const std::string& s, size_t i) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c < 0x80) return std::isalpha(c) || c == '_';
    bool letter = false;
    return eesNonAsciiNameCharLength(s, i, &letter) > 0 && letter;
}

/// Position of the first character of `s` that a name cannot hold; npos when `s` is a valid name.
static size_t firstInvalidNameChar(const std::string& s) {
    if (s.empty() || !nameStartsAt(s, 0)) return 0;
    for (size_t i = 0; i < s.size();) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c >= 0x80) {
            const size_t len = eesNonAsciiNameCharLength(s, i);
            if (len == 0) return i;
            i += len;
            continue;
        }
        if (c == '$') {
            if (i != s.size() - 1) return i;                  // only at the end (string variable)
        } else if (c == '\\') {
            if (i + 1 >= s.size() || !nameStartsAt(s, i + 1)) return i;   // joins two parts of a name
        } else if (!std::isalnum(c) && c != '_' && c != '%' && c != '|' && c != '#') {
            return i;
        }
        ++i;
    }
    return std::string::npos;
}

/// "the character 'x'" for the character of `s` at `i`, or the byte when it is not valid UTF-8.
static std::string describeCharAt(const std::string& s, size_t i) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c < 0x80) return std::string("the character '") + s[i] + "'";
    const size_t len = c >= 0xF0 && c <= 0xF4 ? 4 : c >= 0xE0 && c < 0xF0 ? 3 : c >= 0xC2 && c < 0xE0 ? 2 : 0;
    bool valid = len > 0 && i + len <= s.size();
    unsigned cp = len == 2 ? (c & 0x1Fu) : len == 3 ? (c & 0x0Fu) : (c & 0x07u);
    for (size_t k = 1; valid && k < len; ++k) {
        const unsigned char cc = static_cast<unsigned char>(s[i + k]);
        valid = (cc & 0xC0) == 0x80;
        cp = (cp << 6) | (cc & 0x3Fu);
    }
    char buf[128];
    if (!valid) {
        std::snprintf(buf, sizeof buf, "the byte 0x%02X, which is not valid UTF-8 (save the model as a UTF-8 file),", c);
        return buf;
    }
    std::snprintf(buf, sizeof buf, "' (U+%04X)", cp);
    return "the character '" + s.substr(i, len) + buf;
}

/// Why `token` (no blank in it) is neither a number nor a valid name; empty when the token does
/// not look like a name at all (the caller then keeps its generic message).
static std::string invalidNameMessage(const std::string& token) {
    if (token.empty() || token.find_first_of(" \t") != std::string::npos || token.find("..") != std::string::npos)
        return "";   // not a name: an expression, or a range i..j of an array
    const size_t bad = firstInvalidNameChar(token);
    if (bad == std::string::npos) return "";
    bool wordLike = false;   // holds a letter, a digit or a non-ASCII character
    for (unsigned char c : token) wordLike = wordLike || std::isalnum(c) || c >= 0x80;
    if (!wordLike) return "";
    const unsigned char first = static_cast<unsigned char>(token[0]);
    const std::string quoted = "'" + token + "'";
    if (bad == 0 && (std::isdigit(first) || first == '.')) return "";   // a number, or a range 1..5 of an array
    if (bad == 0 && first < 0x80 && !std::isalpha(first) && first != '_') {
        if (std::string("%|#\\$@&~?!").find(token[0]) == std::string::npos) return "";
        return "Invalid variable name " + quoted + ": a name must start with a letter or an underscore";
    }
    const unsigned char c = static_cast<unsigned char>(token[bad]);
    if (c == '$') return "Invalid variable name " + quoted + ": '$' can only end a name (a string variable)";
    if (c == '\\')
        return "Invalid variable name " + quoted + ": '\\' must be followed by a letter or an underscore";
    const std::string what = describeCharAt(token, bad);
    if (c >= 0x80 && what.find("not valid UTF-8") != std::string::npos) return "Invalid variable name: " + what;
    return "Invalid variable name " + quoted + ": " + what + " is not allowed in a name";
}

// ============================================================================
// EES Grammar (PEG format)
// ============================================================================

// Note: EES is case-insensitive for keywords and functions, but variable names
// preserve case. We handle case-insensitivity at the semantic level.

static const char* EES_GRAMMAR = R"(
    # Top-level program
    Program         <- Spacing Statement* EndOfFile
    
    # Statements
    Statement       <- BlockComment / LineComment / Directive / Equation / EmptyLine
    
    # Comments
    BlockComment    <- '"' < (!'"' .)* > '"' Spacing
    BraceComment    <- '{' < (![}] .)* > '}' Spacing
    LineComment     <- '//' < (![\r\n] .)* > Spacing
    
    # Directives ($if, $endif, etc.)
    Directive       <- '$' < [a-zA-Z_][a-zA-Z0-9_]* > DirectiveContent Spacing
    DirectiveContent <- < (![\r\n] .)* >
    
    # Equations: LHS = RHS with optional units annotation
    Equation        <- Expression '=' Expression Units? InlineComment? Spacing
    Units           <- '"[' < (!']"' .)* > ']"'
    InlineComment   <- '"' < (!'"' .)* > '"' / BraceComment
    
    # Expressions with operator precedence
    Expression      <- Additive
    Additive        <- Multiplicative (('+' / '-') Multiplicative)*
    Multiplicative  <- Power (('*' / '/') Power)*
    Power           <- Unary ('^' Unary)*
    Unary           <- ('-' / '+')? Primary
    
    # Primary expressions
    Primary         <- FunctionCall / ArrayAccess / Variable / Number / StringLiteral / '(' Expression ')'
    
    # Function calls: func(args) or func(Fluid$, T=val, P=val)
    FunctionCall    <- Identifier '(' ArgList? ')' Spacing
    ArgList         <- Arg (',' Arg)*
    Arg             <- NamedArg / Expression
    NamedArg        <- Identifier '=' Expression
    
    # Array access: var[index] or var[index1, index2]
    ArrayAccess     <- Identifier '[' Expression (',' Expression)* ']' Spacing
    
    # Variables
    Variable        <- Identifier Spacing
    
    # Identifiers (case-insensitive comparison for keywords, preserve for vars)
    Identifier      <- < [a-zA-Z_$] [a-zA-Z0-9_$#]* > Spacing
    
    # Numbers: integer, float, scientific notation
    Number          <- < '-'? [0-9]+ ('.' [0-9]*)? ([eE] [+-]? [0-9]+)? > Spacing
    
    # String literals (for fluid names, etc.)
    StringLiteral   <- "'" < (!"'" .)* > "'" Spacing
    
    # Whitespace and empty lines
    Spacing         <- ([ \t] / Comment)*
    Comment         <- BraceComment
    EmptyLine       <- [\r\n]+ Spacing
    EndOfFile       <- !.
    
    # Handle newlines
    ~_              <- [ \t\r\n]*
)";

// ============================================================================
// Parser Implementation
// ============================================================================

class EESParser::Impl {
public:
    Impl() {
        initializeGrammar();
    }
    
    ParseResult parse(const std::string& source, const std::string& filename) {
        presetConstants_.clear();
        ParseResult result = parsePass(source, filename, false);
        // `$IF INTEGRALTABLE` is true when an $IntegralTable directive is in the file, above or
        // below the test: EES compiles in two passes for it, so do we. A $CONSTANT is global: the
        // body of a routine written above its directive uses it too (second pass, constants known).
        if ((integralTableTested_ && integralTableSeen_ != compileState_.integralTable) || constantAfterRoutine_) {
            for (const auto& entry : userConstants_) presetConstants_.insert(entry.first);
            result = parsePass(source, filename, integralTableSeen_);
        }
        return result;
    }

    std::string getLastError() const { return lastError_; }

private:
    ParseResult parsePass(const std::string& source, const std::string& filename, bool integralTableAssumed) {
        conditionalStack_.clear();
        compileState_ = CompileTimeState{};
        compileState_.integralTable = integralTableAssumed;
        integralTableSeen_ = false;
        integralTableTested_ = false;
        procedureNames_.clear();
        userConstants_.clear();
        routineChecked_ = constantAfterRoutine_ = rangeSeen_ = duplicateSeen_ = indexedAggregateSeen_ = unitIndexSeen_ = false;
        lastIgnoredEndifLine_ = 0;
        ParseResult result;
        result.success = false;
        result.program.sourceFilename = filename;
        
        if (!grammarValid_) {
            result.errors.push_back({0, 0, "Grammar initialization failed: " + lastError_, ""});
            return result;
        }
        
        // Preprocess: normalize line endings and handle EES-specific quirks
        std::string preprocessed = preprocess(source);
        
        // Parse line by line for better error recovery
        std::istringstream stream(preprocessed);
        std::string line;
        int lineNumber = 0;
        
        bool unsupportedFound = false;
        comments_ = CommentState{};
        while (std::getline(stream, line)) {
            lineNumber++;

            // Skip comments (also those spanning several lines) and blank lines
            if (skipCommentLine(line, lineNumber, result)) continue;
            std::string trimmed = trim(line);

            // $IF / $IFNOT / $ELSE / $ENDIF; the lines of a removed branch are not looked at
            if (handleConditionalDirective(trimmed, lineNumber, line, result)) continue;
            if (!conditionalActive()) continue;

            // Check for unsupported constructs (procedures, modules, etc.)
            std::string upper = trimmed;
            std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
            
            // Match keywords at start of line with word boundary
            auto isKeyword = [&](const std::string& kw) {
                if (upper.find(kw) != 0) return false;
                if (upper.size() == kw.size()) return true;
                return !isIdentChar(upper[kw.size()]);
            };

            if (isKeyword("MODULE") || 
                isKeyword("SUBPROGRAM")) {
                result.errors.push_back({lineNumber, 0, "Construct '" + trimmed.substr(0, trimmed.find_first_of(" (")) + "' is not yet handled by coolsolve", line});
                unsupportedFound = true;
                continue;
            }

            // Handle Function Definition
            if (isKeyword("FUNCTION")) {
                std::string name;
                std::vector<std::string> params;
                if (tryParseFunctionHeader(line, name, params)) {
                    std::vector<StmtPtr> body;
                    int startLine = lineNumber;
                    procedureBodyDepth_++;
                    routineKind_ = "FUNCTION"; routineName_ = name;
                    const size_t errorsBefore = result.errors.size();
                    std::string term = parseBodyLines(stream, lineNumber, body, result);
                    procedureBodyDepth_--;
                    if (term != "END") {
                        result.errors.push_back({startLine, 0, "Function '" + name + "' missing END", line});
                    } else {
                        if (result.errors.size() == errorsBefore)
                            checkRoutineBody(params, {name}, body, startLine, preprocessed, result);
                        result.program.statements.push_back(makeFunctionDefinition(name, params, body, startLine));
                    }
                    continue;
                }
            }

            // Handle Procedure Definition
            if (isKeyword("PROCEDURE")) {
                std::string name;
                std::vector<std::string> inputs, outputs;
                if (tryParseProcedureHeader(line, name, inputs, outputs)) {
                    procedureNames_.insert(name);
                    std::vector<StmtPtr> body;
                    int startLine = lineNumber;
                    procedureBodyDepth_++;
                    routineKind_ = "PROCEDURE"; routineName_ = name;
                    const size_t errorsBefore = result.errors.size();
                    std::string term = parseBodyLines(stream, lineNumber, body, result);
                    procedureBodyDepth_--;
                    if (term != "END") {
                        result.errors.push_back({startLine, 0, "Procedure '" + name + "' missing END", line});
                    } else {
                        if (result.errors.size() == errorsBefore)
                            checkRoutineBody(inputs, outputs, body, startLine, preprocessed, result);
                        result.program.statements.push_back(makeProcedureDefinition(name, inputs, outputs, body, startLine));
                    }
                    continue;
                }
            }

            // IF ... THEN statements are not valid in the main program (EES)
            if (isKeyword("IF") && findKeywordAtDepth0(trim(removeInlineComments(line)), "THEN") != std::string::npos) {
                std::vector<StmtPtr> discarded;
                const int ifLine = lineNumber;
                parseIfStatement(trim(removeInlineComments(line)), stream, lineNumber, discarded, result, line);
                result.errors.push_back({ifLine, 0, kIfOutsideProcedureMessage, line});
                continue;
            }

            // GOTO and statement labels belong to the body of a FUNCTION or PROCEDURE (CS-GAP-GOTO)
            if (isKeyword("GOTO") || std::isdigit(static_cast<unsigned char>(trimmed.empty() ? ' ' : trimmed[0]))) {
                std::vector<StmtPtr> discarded;
                const std::string code = trim(removeInlineComments(line));
                int labelNumber = 0;
                std::string labelRest;
                if (parseGotoSegment(code, lineNumber, discarded, result, line)) continue;
                if (splitStatementLabel(code, labelNumber, labelRest)) {
                    result.errors.push_back({lineNumber, 0, kLabelOutsideRoutineMessage, line});
                    continue;
                }
            }

            // Handle Procedure Call
            if (isKeyword("CALL")) {
                std::string callError;
                if (auto stmt = tryParseProcedureCall(line, lineNumber, &callError)) {
                    result.program.statements.push_back(stmt);
                    noteComputedStrings(stmt);
                    continue;
                }
                if (!callError.empty()) {
                    result.errors.push_back({lineNumber, 0, callError, line});
                    continue;
                }
            }

            // Handle DUPLICATE loop: DUPLICATE i = start, end
            if (isKeyword("DUPLICATE")) {
                int startLine = lineNumber;
                // Parse header: DUPLICATE varName = startExpr, endExpr
                std::string header = trim(trimmed.substr(9));  // skip "DUPLICATE"
                // Find '=' to get iterator variable
                auto eqPos = header.find('=');
                if (eqPos == std::string::npos) {
                    result.errors.push_back({lineNumber, 0, "DUPLICATE missing '=' in header", line});
                    continue;
                }
                std::string iterVar = trim(header.substr(0, eqPos));
                std::string rangeStr = trim(header.substr(eqPos + 1));
                // Find comma separating start and end
                auto commaPos = rangeStr.find(',');
                if (commaPos == std::string::npos) {
                    result.errors.push_back({lineNumber, 0, "DUPLICATE missing ',' between start and end", line});
                    continue;
                }
                std::string startStr = trim(rangeStr.substr(0, commaPos));
                std::string endStr = trim(rangeStr.substr(commaPos + 1));
                auto startExpr = parseExpression(startStr, lineNumber);
                auto endExpr = parseExpression(endStr, lineNumber);
                if (!startExpr || !endExpr) {
                    result.errors.push_back({lineNumber, 0, "DUPLICATE could not parse range expressions", line});
                    continue;
                }
                // Parse body until END
                std::vector<StmtPtr> body;
                std::string term = parseBodyLines(stream, lineNumber, body, result);
                if (term != "END") {
                    result.errors.push_back({startLine, 0, "DUPLICATE missing END", line});
                } else {
                    result.program.statements.push_back(makeDuplicate(iterVar, startExpr, endExpr, body, startLine));
                    duplicateSeen_ = true;
                }
                continue;
            }

            // REPEAT-UNTIL loops are not valid in the main program (EES): parsed to skip them, then reported
            if (isKeyword("REPEAT")) {
                std::vector<StmtPtr> discarded;
                parseRepeatStatement(stream, lineNumber, discarded, result, line);
                continue;
            }
            
            // Handle comments
            if (auto stmt = tryParseComment(line, lineNumber)) {
                result.program.statements.push_back(stmt);
                result.commentCount++;
                continue;
            }
            
            // Handle directives
            if (auto stmt = tryParseDirective(line, lineNumber, &result.diagnostics)) {
                result.program.statements.push_back(stmt);
                result.directiveCount++;
                noteConstantDirective(stmt->as<Directive>(), lineNumber, line, result);
                if (upper.rfind("$INTEGRALTABLE", 0) == 0) integralTableSeen_ = true;
                continue;
            }
            
            // Handle equations — split on semicolons first (EES uses ; as
            // a statement separator, e.g. `a = 1 ; b = 2 ; c = 3`).
            {
                auto segments = splitOnSemicolons(line);
                if (segments.size() > 1) {
                    // Multiple statements on one line
                    bool allOk = true;
                    for (const auto& seg : segments) {
                        std::string segError;
                        if (auto stmt = tryParseEquationOrAssignment(seg, lineNumber, &segError)) {
                            result.program.statements.push_back(stmt);
                            result.equationCount++;
                            noteStringAssignment(stmt);
                        } else {
                            allOk = false;
                            std::string msg = segError.empty() ? ("Could not parse segment: " + seg) : segError;
                            result.errors.push_back({lineNumber, 0, msg, line});
                        }
                    }
                    if (allOk) continue;
                    // Even if some segments failed, we already pushed the
                    // successful ones; skip the fallback single-parse path.
                    continue;
                }
            }

            {
                std::string parseError;
                if (auto stmt = tryParseEquationOrAssignment(line, lineNumber, &parseError)) {
                    result.program.statements.push_back(stmt);
                    result.equationCount++;
                    noteStringAssignment(stmt);
                    continue;
                }
                if (!parseError.empty()) {
                    result.errors.push_back({lineNumber, 0, parseError, line});
                    continue;
                }
            }
            
            // If we couldn't parse the line, record an error but continue
            // We don't set unsupportedFound here, as it might just be a line we don't understand yet
            result.errors.push_back({lineNumber, 0, "Could not parse line", line});
        }
        
        for (auto it = conditionalStack_.rbegin(); it != conditionalStack_.rend(); ++it) {
            const auto& info = result.conditionalDirectives[it->infoIndex];
            result.errors.push_back({it->line, 0,
                std::string(it->ifnot ? "$IFNOT" : "$IF") + " without a matching $ENDIF", info.text});
        }
        if (comments_.inQuote) {
            result.diagnostics.push(DiagnosticSeverity::Warning, "P007",
                "The comment opened by a '\"' at line " + std::to_string(comments_.quoteOpenLine) +
                " is never closed: everything after it is ignored",
                "parser", comments_.quoteOpenLine, 0);
        }
        result.totalLines = lineNumber;
        checkRoutineCycles(result, preprocessed);
        if (result.errors.empty() && (rangeSeen_ || duplicateSeen_ || indexedAggregateSeen_ || unitIndexSeen_ ||
                                      !userConstants_.empty()))
            resolveCompileTimeValues(result, preprocessed);
        // Success if no errors AND we found something valid
        bool somethingParsed = result.equationCount > 0 || result.commentCount > 0 || result.directiveCount > 0 || !result.program.statements.empty();
        result.success = !unsupportedFound && result.errors.empty() && somethingParsed;
        
        // Post-parse validation: check function names against known builtins
        validateFunctionNames(result);
        validateThermoUserPatterns(result);
        
        return result;
    }
    
    peg::parser parser_;
    bool grammarValid_ = false;
    std::string lastError_;
    // Set when an argument of a call, or an index, cannot be parsed (or when a named
    // argument is given to a function that takes none); consumed (and cleared) by the
    // statement parsers to emit a specific diagnostic instead of the generic
    // "Could not parse line". The first message is kept: it is that of the innermost call.
    std::string pendingExpressionError_;
    // Why a name was refused (invalidNameMessage): added to the specific message above, or
    // alone, so that a character EES takes in a name and CoolSolve does not is named.
    std::string pendingNameHint_;

    void setExpressionError(const std::string& message) {
        if (pendingExpressionError_.empty()) pendingExpressionError_ = message;
    }
    void clearExpressionError() {
        pendingExpressionError_.clear();
        pendingNameHint_.clear();
    }
    /// The diagnostic to report for an expression that cannot be parsed; empty when nothing specific is known.
    std::string expressionErrorText() const {
        if (pendingNameHint_.empty()) return pendingExpressionError_;
        if (pendingExpressionError_.empty()) return pendingNameHint_;
        return pendingExpressionError_ + "; " + pendingNameHint_;
    }
    // > 0 while the body of a FUNCTION/PROCEDURE is being parsed (IF ... THEN
    // statements are only valid there).
    int procedureBodyDepth_ = 0;
    // The routine whose body is being parsed ("FUNCTION" or "PROCEDURE", and its name).
    std::string routineKind_, routineName_;
    // Text that follows `ENDIF;` on a body line, returned by nextBodyLine before the stream is read.
    std::string pendingBodyText_;
    int lastIgnoredEndifLine_ = 0;  // the line of the last ENDIF ignored in a body (one warning per line)
    // The PROCEDUREs defined so far (EES: a routine is defined before it is called); one of them
    // named like a built-in procedure replaces it (see isBuiltinColebrook).
    std::set<std::string, CaseInsensitiveLess> procedureNames_;

    // $CONSTANT and the other values known before the solve (see resolveCompileTimeValues)
    struct UserConstant {
        double value;
        int line;
    };
    std::map<std::string, UserConstant, CaseInsensitiveLess> userConstants_;  // the $CONSTANTs parsed so far
    std::set<std::string, CaseInsensitiveLess> presetConstants_;  // all those of the file, on a second pass
    bool routineChecked_ = false;        // the body of a FUNCTION/PROCEDURE has been checked
    bool constantAfterRoutine_ = false;  // a $CONSTANT defined below such a body: parse again
    bool rangeSeen_ = false;             // an array range X[a..b] was parsed (a header included)
    bool duplicateSeen_ = false;         // a DUPLICATE loop was parsed
    bool indexedAggregateSeen_ = false;  // SUM(expr, i=a, b) or PRODUCT(expr, i=a, b) was parsed
    bool unitIndexSeen_ = false;         // an element indexed by a bare unit name (a[kW]) was parsed

    // ------------------------------------------------------------------------
    // Compile-time directives $IF / $IFNOT / $ELSE / $ENDIF (CS-GAP-IF-DIRECTIVE)
    //
    // EES compiles the Equations window top to bottom and ignores the lines of the
    // branch that is not selected. The parse loops (main program and bodies) call
    // handleConditionalDirective() on every line that is not a comment and skip the
    // line while conditionalActive() is false: a removed line is never parsed (it may
    // hold anything, even text CoolSolve does not support), so it never reaches the IR,
    // the debug reports or the solver, whichever entry point built the program.
    // The string constants that the conditions test are recorded as the statements
    // that set them are parsed: EES needs the variable defined before the directive.
    // The conditions themselves are evaluated by conditional_directives.cpp.
    // ------------------------------------------------------------------------
    struct ConditionalFrame {
        int line;              // line of the $IF / $IFNOT
        size_t infoIndex;      // its entry in ParseResult::conditionalDirectives
        bool ifnot;
        bool parentActive;     // the enclosing branch is kept
        bool failed;           // the condition cannot be evaluated (reported): both branches are skipped
        bool condition;        // the THEN branch is selected
        bool inElse;           // $ELSE seen
    };
    std::vector<ConditionalFrame> conditionalStack_;
    CompileTimeState compileState_;
    bool integralTableSeen_ = false;     // a kept $IntegralTable directive was parsed
    bool integralTableTested_ = false;   // a condition asked for INTEGRALTABLE

    static bool frameActive(const ConditionalFrame& f) {
        return f.parentActive && !f.failed && (f.inElse ? !f.condition : f.condition);
    }
    bool conditionalActive() const { return conditionalStack_.empty() || frameActive(conditionalStack_.back()); }

    // Consume the line if it is $IF, $IFNOT, $ELSE or $ENDIF (any letter case). Text after
    // the condition other than a comment or a trailing `then` is an error; text after
    // $ELSE / $ENDIF is ignored (EES files annotate them: `$else fuel$='ch4'`). A bare
    // `$ELSEIF` is a `$ELSE`, as in the EES help.
    bool handleConditionalDirective(const std::string& trimmed, int lineNumber, const std::string& sourceLine,
                                    ParseResult& result) {
        if (trimmed.size() < 3 || trimmed[0] != '$') return false;
        size_t nameEnd = 1;
        while (nameEnd < trimmed.size() &&
               (std::isalnum(static_cast<unsigned char>(trimmed[nameEnd])) || trimmed[nameEnd] == '_'))
            nameEnd++;
        std::string name = trimmed.substr(1, nameEnd - 1);
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        if (name != "if" && name != "ifnot" && name != "else" && name != "elseif" && name != "endif") return false;
        result.directiveCount++;

        std::string upperName = name;
        std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
        const std::string rest = trim(removeInlineComments(trimmed.substr(nameEnd)));
        auto fail = [&](const std::string& message) {
            result.errors.push_back({lineNumber, 0, message, sourceLine});
        };

        if (name == "if" || name == "ifnot") {
            ConditionalFrame frame{lineNumber, result.conditionalDirectives.size(), name == "ifnot",
                                   conditionalActive(), false, false, false};
            ConditionalDirectiveInfo info;
            info.line = lineNumber;
            info.text = trim(removeInlineComments(trimmed));
            if (frame.parentActive) {
                const ConditionResult r = evaluateDirectiveCondition(rest, compileState_);
                if (r.usedIntegralTable) integralTableTested_ = true;
                if (r.ok) {
                    frame.condition = frame.ifnot ? !r.value : r.value;
                    info.evaluated = true;
                    info.thenKept = frame.condition;
                } else {
                    frame.failed = true;
                    fail("$" + upperName + ": " + r.error);
                }
            }
            result.conditionalDirectives.push_back(info);
            conditionalStack_.push_back(frame);
            return true;
        }

        if (conditionalStack_.empty()) {
            fail("$" + upperName + " without a matching $IF or $IFNOT");
            return true;
        }
        ConditionalFrame& frame = conditionalStack_.back();
        ConditionalDirectiveInfo& info = result.conditionalDirectives[frame.infoIndex];
        if (name == "endif") {
            info.endLine = lineNumber;
            conditionalStack_.pop_back();
            return true;
        }
        // $ELSE or $ELSEIF
        if (frame.inElse) {
            fail("Second $" + upperName + " in the " + (frame.ifnot ? "$IFNOT" : "$IF") + " of line " + std::to_string(frame.line));
            return true;
        }
        if (name == "elseif" && !rest.empty()) {
            frame.failed = true;
            fail("$ELSEIF with a condition is not supported: the EES help documents $IF, $IFNOT, $ELSE and $ENDIF "
                 "(write $ELSE followed by a nested $IF)");
            return true;
        }
        frame.inElse = true;
        info.elseLine = lineNumber;
        return true;
    }

    // What a main-program statement tells the $IF conditions about string variables: a
    // constant (`X$ = 'text'`, `X$ = Y$` with Y$ a constant) or a value unknown at compile time.
    void noteStringAssignment(const StmtPtr& stmt) {
        if (!stmt->is<Equation>()) return;
        const auto& eq = stmt->as<Equation>();
        if (!eq.lhs || !eq.rhs || !eq.lhs->is<Variable>()) return;
        const auto& var = eq.lhs->as<Variable>();
        if (var.name.empty() || var.name.back() != '$' || !var.indices.empty()) return;
        const std::string* constant = nullptr;
        if (eq.rhs->is<StringLiteral>()) {
            constant = &eq.rhs->as<StringLiteral>().value;
        } else if (eq.rhs->is<Variable>() && eq.rhs->as<Variable>().indices.empty()) {
            const auto it = compileState_.strings.find(eq.rhs->as<Variable>().name);
            if (it != compileState_.strings.end()) constant = &it->second;
        }
        if (constant) {
            const std::string value = *constant;   // may alias the entry being replaced
            compileState_.strings[var.name] = value;
            compileState_.computedStrings.erase(var.name);
        } else {
            compileState_.strings.erase(var.name);
            compileState_.computedStrings.insert(var.name);
        }
    }

    // The string outputs of a CALL are computed, not constants.
    void noteComputedStrings(const StmtPtr& stmt) {
        if (!stmt->is<ProcedureCall>()) return;
        for (const auto& out : stmt->as<ProcedureCall>().outputVars) {
            if (out.name.empty() || out.name.back() != '$' || !out.indices.empty()) continue;
            compileState_.strings.erase(out.name);
            compileState_.computedStrings.insert(out.name);
        }
    }

    // Known built-in math/utility functions (case-insensitive)
    static const std::unordered_set<std::string>& knownBuiltinFunctions() {
        static const std::unordered_set<std::string> fns = {
            "sin", "cos", "tan", "exp", "ln", "log", "log10", "sqrt", "abs",
            "asin", "arcsin", "acos", "arccos", "atan", "arctan", "atan2",
            "sinh", "cosh", "tanh", "arcsinh", "asinh", "arccosh", "acosh",
            "arctanh", "atanh",
            "ceil", "floor", "round", "trunc", "sign", "mod", "pow",
            "nday_", "numday_", "aveday_", "eqntime_",  // EES time functions (SETP.LIB)
            "min", "max", "sum", "sum2d", "average", "product", "stddev",
            "if", "convert", "converttemp", "pi",
            "tablevalue", "tablevalue#", "tablerun#", "lookup", "lookup$", "uppercase$", "lowercase$",
            "interpolate", "interpolate1", "interpolate2", "interpolate2dm",
            "nlookuprows", "nlookupcolumns",
            "sumlookup", "avglookup", "maxlookup", "minlookup", "stddevlookup",
            "lookupcol", "lookupcol1", "lookupcellempty",
            "unitsystem", "error", "warning",  // UNITSYSTEM('K') and the built-in procedures CALL error('msg'), CALL warning('msg')
            "colebrook",  // the built-in procedure CALL Colebrook(Re, RR : f)
            "ntumethod", "effectivenessntu",
            "lmtd",
            "integral", "integralvalue",  // equation-based dynamic solver
        };
        return fns;
    }
    
    // Built-in functions that take positional arguments only: for them `name=value` is
    // not a named argument (typically the comparison of IF(a=5, 1, 2), which EES does not
    // have). The indexed aggregations SUM(X[i], i=1, N) and the EES form of INTERPOLATE
    // carry `name=value` arguments: the former is expanded into the elements once the file is parsed
    // (resolveIndexedAggregate; indexedAggregateError refuses it in any other aggregate, which are
    // listed here only so that it is the one to answer), the latter is evaluated from the column names
    // of the call (ExpressionEvaluator::resolveEesInterpolation). The thermophysical functions,
    // user functions and unknown names are not listed here.
    static bool isPositionalOnlyBuiltin(const std::string& name) {
        static const std::unordered_set<std::string> takingNamedArguments = {
            "sum", "sum2d", "average", "product", "stddev",
            "interpolate", "interpolate1", "interpolate2", "interpolate2dm",
        };
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        return knownBuiltinFunctions().count(lower) > 0 && takingNamedArguments.count(lower) == 0;
    }

    static std::string namedArgumentRefusedMessage(const std::string& function, const std::string& arg,
                                                   const std::string& paramName, const std::string& paramValue) {
        std::string lower = function;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        std::string message = function + " does not take named arguments ('" + arg + "'); ";
        if (lower == "if")
            return message + "to test an equality use the 5-argument form IF(A, B, X, Y, Z) "
                   "(X if A<B, Y if A=B, Z if A>B), here IF(" + paramName + ", " + paramValue
                   + ", X, Y, Z); IF(cond, a, b) tests cond > 0";
        return message + "only the thermophysical property functions do, e.g. enthalpy(Water, T=20, P=1E5)";
    }

    /// Receives the calls found by forEachCall: the name, whether it is a CALL statement (false: a call in an
    /// expression, `f(x)`) and the line.
    using CallVisitor = std::function<void(const std::string& name, bool isCallStatement, int line)>;

    /// Every call of an expression, the calls in the arguments of a call and in the indices of a variable
    /// included, outermost first.
    static void forEachCall(const ExprPtr& expr, const CallVisitor& visit) {
        if (!expr) return;
        std::visit([&](const auto& node) {
            using T = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<T, FunctionCall>) {
                visit(node.name, false, expr->sourceLineNumber);
                for (const auto& arg : node.args) forEachCall(arg, visit);
                for (const auto& [k, v] : node.namedArgs) forEachCall(v, visit);
            } else if constexpr (std::is_same_v<T, UnaryOp>) {
                forEachCall(node.operand, visit);
            } else if constexpr (std::is_same_v<T, BinaryOp>) {
                forEachCall(node.left, visit);
                forEachCall(node.right, visit);
            } else if constexpr (std::is_same_v<T, Variable>) {
                for (const auto& idx : node.indices) forEachCall(idx, visit);
            }
        }, expr->node);
    }

    /// Every call of a statement and of the statements nested in it (IF branches, loop bodies, bodies of
    /// routines): the calls of its expressions and its CALL statements.
    static void forEachCall(const StmtPtr& stmt, const CallVisitor& visit) {
        if (!stmt) return;
        std::visit([&](const auto& node) {
            using T = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<T, Equation> || std::is_same_v<T, Assignment>) {
                forEachCall(node.lhs, visit);
                forEachCall(node.rhs, visit);
            } else if constexpr (std::is_same_v<T, IfThenElse>) {
                forEachCall(node.condition, visit);
                for (const auto& s : node.thenBranch) forEachCall(s, visit);
                for (const auto& s : node.elseBranch) forEachCall(s, visit);
            } else if constexpr (std::is_same_v<T, Duplicate>) {
                forEachCall(node.start, visit);
                forEachCall(node.end, visit);
                for (const auto& s : node.body) forEachCall(s, visit);
            } else if constexpr (std::is_same_v<T, RepeatUntil>) {
                forEachCall(node.condition, visit);
                for (const auto& s : node.body) forEachCall(s, visit);
            } else if constexpr (std::is_same_v<T, ProcedureCall>) {
                visit(node.name, true, stmt->sourceLineNumber);
                for (const auto& arg : node.inputArgs) forEachCall(arg, visit);
                for (const auto& out : node.outputExprs) forEachCall(out, visit);
                for (const auto& out : node.outputVars)
                    for (const auto& idx : out.indices) forEachCall(idx, visit);
            } else if constexpr (std::is_same_v<T, FunctionDefinition> || std::is_same_v<T, ProcedureDefinition>) {
                for (const auto& s : node.body) forEachCall(s, visit);
            }
        }, stmt->node);
    }

    // ------------------------------------------------------------------------
    // Call cycles between routines (CS-BUG-RECURSION-CRASH)
    //
    // EES does not allow recursion ("Recursive calls in which a function calls itself are, however, not
    // allowed", Internal Functions), and a body that calls itself, directly or through other routines,
    // overflowed the stack and killed the process. A cycle in the graph "routine R calls routine S" is
    // therefore a parse error, one per cycle found, with the line of the call that closes it. The nodes
    // are the FUNCTIONs and PROCEDUREs of the model; a call that names none of them (a built-in, a
    // property function, an unknown name) is no edge. A call goes where the evaluator sends it: a call in
    // an expression to the FUNCTION of that name or else to the PROCEDURE (the inline form `y = p(x)`),
    // a CALL statement to the PROCEDURE, so that a FUNCTION and a PROCEDURE of the same name are distinct
    // routines and a model PROCEDURE Colebrook (which replaces the built-in) is a node like any other.
    // ------------------------------------------------------------------------
    void checkRoutineCycles(ParseResult& result, const std::string& source) {
        struct Routine {
            std::string label;                           // "FUNCTION f"
            std::vector<std::pair<int, int>> calls;      // (callee, line of its first call here)
            int state = 0;                               // 0 not visited, 1 on the search path, 2 done
        };
        std::vector<Routine> routines;
        std::map<std::string, int, CaseInsensitiveLess> functions, procedures;   // name -> index in `routines`
        auto defined = [&](const StmtPtr& stmt) -> int {   // the routine a definition declares (added if new), -1 if none
            const bool isFunction = stmt->is<FunctionDefinition>();
            if (!isFunction && !stmt->is<ProcedureDefinition>()) return -1;
            const std::string& name = isFunction ? stmt->as<FunctionDefinition>().name : stmt->as<ProcedureDefinition>().name;
            auto& byName = isFunction ? functions : procedures;
            auto [it, isNew] = byName.emplace(name, static_cast<int>(routines.size()));
            if (isNew) routines.push_back({std::string(isFunction ? "FUNCTION " : "PROCEDURE ") + name, {}, 0});
            return it->second;
        };
        for (const auto& stmt : result.program.statements) defined(stmt);
        if (routines.empty()) return;

        for (const auto& stmt : result.program.statements) {
            const int caller = defined(stmt);
            if (caller < 0) continue;
            const auto& body = stmt->is<FunctionDefinition>() ? stmt->as<FunctionDefinition>().body
                                                              : stmt->as<ProcedureDefinition>().body;
            for (const auto& inner : body) {
                forEachCall(inner, [&](const std::string& name, bool isCallStatement, int line) {
                    int callee = -1;
                    auto function = functions.find(name);
                    if (!isCallStatement && function != functions.end()) {
                        callee = function->second;
                    } else {
                        std::string lower = name;
                        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                        auto procedure = procedures.find(name);
                        // the built-in psychprops always runs, whatever the model defines under that name
                        if (lower != "psychprops" && procedure != procedures.end()) callee = procedure->second;
                    }
                    if (callee < 0) return;
                    auto& calls = routines[caller].calls;
                    const bool known = std::any_of(calls.begin(), calls.end(), [&](const auto& c) { return c.first == callee; });
                    if (!known) calls.emplace_back(callee, line);
                });
            }
        }

        std::vector<int> path;   // the routines being searched, the last one calls the next
        std::function<void(int)> search = [&](int node) {
            routines[node].state = 1;
            path.push_back(node);
            for (const auto& [callee, line] : routines[node].calls) {
                if (routines[callee].state == 0) {
                    search(callee);
                } else if (routines[callee].state == 1) {   // back to a routine of the path: a cycle
                    auto first = std::find(path.begin(), path.end(), callee);
                    std::string message;
                    if (first + 1 == path.end()) {
                        message = "Recursive call: " + routines[callee].label + " calls itself";
                    } else {
                        message = "Recursive calls: ";
                        for (auto it = first; it != path.end(); ++it) message += routines[*it].label + " → ";
                        message += routines[callee].label;
                    }
                    message += "; recursion is not allowed in EES";
                    result.errors.push_back({line, 0, message, sourceLineText(source, line)});
                }
            }
            path.pop_back();
            routines[node].state = 2;
        };
        for (int i = 0; i < static_cast<int>(routines.size()); ++i)
            if (routines[i].state == 0) search(i);
    }

    // Post-parse validation: check that all function calls reference either a
    // built-in, a CoolProp thermo function, or a user-defined function/procedure.
    void validateFunctionNames(ParseResult& result) {
        // Collect user-defined function and procedure names
        std::set<std::string, CaseInsensitiveLess> userDefined;
        for (const auto& stmt : result.program.statements) {
            if (stmt->is<FunctionDefinition>()) {
                userDefined.insert(stmt->as<FunctionDefinition>().name);
            } else if (stmt->is<ProcedureDefinition>()) {
                userDefined.insert(stmt->as<ProcedureDefinition>().name);
            }
        }
        
        // Collect all function call names from the AST
        struct FuncCallInfo { std::string name; int line; };
        std::vector<FuncCallInfo> callsFound;
        for (const auto& stmt : result.program.statements) {
            forEachCall(stmt, [&](const std::string& name, bool, int line) { callsFound.push_back({name, line}); });
        }
        
        // Check each function call
        const auto& builtins = knownBuiltinFunctions();
        std::set<std::string, CaseInsensitiveLess> alreadyWarned;
        for (const auto& call : callsFound) {
            if (alreadyWarned.count(call.name)) continue;
            
            std::string lower = call.name;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            
            if (builtins.count(lower)) continue;
            if (isEESThermophysicalFunction(call.name)) continue;
            if (userDefined.count(call.name)) continue;
            
            const std::string typo = suggestThermoFunctionTypo(call.name);
            if (!typo.empty()) {
                result.diagnostics.push(DiagnosticSeverity::Warning, "P004",
                    "Unknown function '" + call.name + "'. " + typo,
                    "parser", call.line, 0);
                alreadyWarned.insert(call.name);
                continue;
            }
            
            result.diagnostics.push(DiagnosticSeverity::Warning, "P003",
                "Unknown function '" + call.name + "'",
                "parser", call.line, 0);
            alreadyWarned.insert(call.name);
        }
    }

    static bool tryGetNumericLiteral(const ExprPtr& expr, double& out) {
        if (!expr) return false;
        if (expr->is<NumberLiteral>()) {
            out = expr->as<NumberLiteral>().value;
            return true;
        }
        if (expr->is<UnaryOp>()) {
            const auto& u = expr->as<UnaryOp>();
            if (u.op == "-" && u.operand && u.operand->is<NumberLiteral>()) {
                out = -u.operand->as<NumberLiteral>().value;
                return true;
            }
        }
        return false;
    }

    static ThermoInputKind thermoKindFromArgName(const std::string& argName) {
        std::string n = argName;
        std::transform(n.begin(), n.end(), n.begin(), ::tolower);
        if (n == "p" || n == "pressure") return ThermoInputKind::Pressure;
        if (n == "t" || n == "temperature") return ThermoInputKind::Temperature;
        if (n == "h" || n == "enthalpy") return ThermoInputKind::Enthalpy;
        if (n == "s" || n == "entropy") return ThermoInputKind::Entropy;
        if (n == "d" || n == "rho" || n == "density") return ThermoInputKind::Density;
        if (n == "v" || n == "volume") return ThermoInputKind::SpecificVolume;
        if (n == "q" || n == "x" || n == "quality") return ThermoInputKind::Quality;
        if (n == "w") return ThermoInputKind::HumidityRatio;
        if (n == "r") return ThermoInputKind::RelativeHumidity;
        return ThermoInputKind::Temperature;
    }

    void validateThermoUserPatterns(ParseResult& result) {
        resetThermoInputHints();
        UnitSystem units;
        std::set<std::string, CaseInsensitiveLess> fluidWarningsDone;

        auto inspectThermoCall = [&](const FunctionCall& call, int line) {
            if (!isEESThermophysicalFunction(call.name)) return;

            // First positional argument: fluid name patterns
            if (!call.args.empty() && call.args[0]->is<Variable>()) {
                const auto& var = call.args[0]->as<Variable>();
                const bool isStringVar = !var.name.empty() && var.name.back() == '$';
                const std::string hint = unquotedFluidArgumentHint(
                    call.name, var.name, isStringVar);
                if (!hint.empty() && !fluidWarningsDone.count(var.name)) {
                    result.diagnostics.push(DiagnosticSeverity::Warning, "P005",
                                            hint, "parser", line, 0);
                    fluidWarningsDone.insert(var.name);
                }
                if (!isStringVar) {
                    const std::string fluidHint = suggestFluidNameHint(var.name);
                    if (!fluidHint.empty() && !fluidWarningsDone.count("fluid:" + var.name)) {
                        result.diagnostics.push(DiagnosticSeverity::Warning, "P005",
                            call.name + "(): fluid argument '" + var.name + "'. " + fluidHint,
                            "parser", line, 0);
                        fluidWarningsDone.insert("fluid:" + var.name);
                    }
                }
            } else if (!call.args.empty() && call.args[0]->is<StringLiteral>()) {
                const std::string fluid = call.args[0]->as<StringLiteral>().value;
                const std::string fluidHint = suggestFluidNameHint(fluid);
                if (!fluidHint.empty() && !fluidWarningsDone.count("str:" + fluid)) {
                    result.diagnostics.push(DiagnosticSeverity::Warning, "P005",
                        call.name + "(): fluid '" + fluid + "'. " + fluidHint,
                        "parser", line, 0);
                    fluidWarningsDone.insert("str:" + fluid);
                }
            }

            // Literal named arguments: static unit hints with line numbers
            for (const auto& [argName, argExpr] : call.namedArgs) {
                double lit = 0.0;
                if (!tryGetNumericLiteral(argExpr, lit)) continue;
                checkThermoInputValueHint(&result.diagnostics, call.name,
                                          units, thermoKindFromArgName(argName),
                                          argName, lit, line);
            }
        };

        auto walkExpr = [&](const ExprPtr& expr, auto& self) -> void {
            if (!expr) return;
            std::visit([&](const auto& node) {
                using T = std::decay_t<decltype(node)>;
                if constexpr (std::is_same_v<T, FunctionCall>) {
                    inspectThermoCall(node, expr->sourceLineNumber);
                    for (const auto& arg : node.args) self(arg, self);
                    for (const auto& [k, v] : node.namedArgs) self(v, self);
                } else if constexpr (std::is_same_v<T, UnaryOp>) {
                    self(node.operand, self);
                } else if constexpr (std::is_same_v<T, BinaryOp>) {
                    self(node.left, self); self(node.right, self);
                } else if constexpr (std::is_same_v<T, Variable>) {
                    for (const auto& idx : node.indices) self(idx, self);
                }
            }, expr->node);
        };

        auto walkStmt = [&](const StmtPtr& stmt, auto& self) -> void {
            if (!stmt) return;
            std::visit([&](const auto& node) {
                using T = std::decay_t<decltype(node)>;
                if constexpr (std::is_same_v<T, Equation>) {
                    walkExpr(node.lhs, walkExpr);
                    walkExpr(node.rhs, walkExpr);
                } else if constexpr (std::is_same_v<T, Assignment>) {
                    walkExpr(node.lhs, walkExpr);
                    walkExpr(node.rhs, walkExpr);
                } else if constexpr (std::is_same_v<T, IfThenElse>) {
                    walkExpr(node.condition, walkExpr);
                    for (const auto& s : node.thenBranch) self(s, self);
                    for (const auto& s : node.elseBranch) self(s, self);
                } else if constexpr (std::is_same_v<T, Duplicate>) {
                    walkExpr(node.start, walkExpr);
                    walkExpr(node.end, walkExpr);
                    for (const auto& s : node.body) self(s, self);
                } else if constexpr (std::is_same_v<T, RepeatUntil>) {
                    walkExpr(node.condition, walkExpr);
                    for (const auto& s : node.body) self(s, self);
                } else if constexpr (std::is_same_v<T, ProcedureCall>) {
                    for (const auto& arg : node.inputArgs) walkExpr(arg, walkExpr);
                    for (const auto& out : node.outputExprs) walkExpr(out, walkExpr);
                } else if constexpr (std::is_same_v<T, FunctionDefinition>) {
                    for (const auto& s : node.body) self(s, self);
                } else if constexpr (std::is_same_v<T, ProcedureDefinition>) {
                    for (const auto& s : node.body) self(s, self);
                }
            }, stmt->node);
        };

        for (const auto& stmt : result.program.statements) {
            walkStmt(stmt, walkStmt);
        }
    }
    
    void initializeGrammar() {
        // Use a simpler grammar approach - parse expressions directly
        grammarValid_ = true;  // We'll use manual parsing
    }

    // State of comments that span several lines. Shared by the main parse loop
    // and the FUNCTION/PROCEDURE body parser, so that both skip them alike.
    struct CommentState {
        bool inQuote = false;      // inside a multi-line "..." comment
        bool standalone = false;   // block opened by a line that is exactly `"`
        int quoteOpenLine = 0;     // line that opened the multi-line "..." comment
        int braceDepth = 0;        // nesting depth of { } comments
    };
    CommentState comments_;

    // Position of a `"` that opens a comment not closed on the same line (the comment
    // continues on the next lines), or npos. Mirrors removeInlineComments(): single-quoted
    // strings, `"[units]"` annotations and `//` comments are skipped.
    static size_t findUnterminatedQuote(const std::string& line) {
        bool inString = false;
        for (size_t i = 0; i < line.size(); ++i) {
            const char c = line[i];
            if (c == '\'') { inString = !inString; continue; }
            if (inString) continue;
            if (c == '/' && i + 1 < line.size() && line[i + 1] == '/') return std::string::npos;
            if (c != '"') continue;
            if (i + 1 < line.size() && line[i + 1] == '[') {
                const size_t units = line.find("]\"", i);
                if (units != std::string::npos) { i = units + 1; continue; }
            }
            const size_t close = line.find('"', i + 1);
            if (close == std::string::npos) return i;
            i = close;
        }
        return std::string::npos;
    }

    // Remove from `line` (in place) what belongs to multi-line comments and
    // update comments_. Returns true when nothing is left to parse on the line:
    // it is blank or entirely inside a comment. `lineNumber` is the line of `line`.
    bool skipCommentLine(std::string& line, int lineNumber, ParseResult& result) {
        // Multi-line quote comment " ... "
        if (comments_.inQuote) {
            if (trim(line) == "\"") {
                // Standalone `"` always ends a quote block (and starts one when not in block).
                comments_.inQuote = false;
                result.commentCount++;
                return true;
            }
            // Block started with a standalone `"` -> only a standalone `"` can close it.
            if (comments_.standalone) return true;
            // Classic style: block started with `"` on a line with no closing `"` -> first `"` on a line closes.
            const size_t quoteEnd = line.find('"');
            if (quoteEnd == std::string::npos) return true;
            comments_.inQuote = false;
            result.commentCount++;
            line = line.substr(quoteEnd + 1);
        }

        // Nested brace comments that can span multiple lines. Done before any
        // further processing so that code inside `{ ... }` is completely ignored.
        line = stripBraceComments(line, comments_.braceDepth, result);
        if (comments_.braceDepth > 0 && line.empty()) return true;  // still inside a block comment

        if (isEmptyOrWhitespace(line)) return true;
        const std::string trimmed = trim(line);

        // Start of multi-line quote: (1) line is exactly `"` (standalone delimiter), or
        // (2) a `"` with no closing `"` on the same line (classic EES style): the comment
        // runs up to the first `"` of a later line. Code before the `"` is kept, so that
        // `x = 1 "a comment <newline> continued here"` is the equation x = 1.
        if (trimmed == "\"") {
            comments_.inQuote = true;
            comments_.standalone = true;
            comments_.quoteOpenLine = lineNumber;
            return true;
        }
        const size_t open = findUnterminatedQuote(line);
        if (open != std::string::npos) {
            comments_.inQuote = true;
            comments_.standalone = false;
            comments_.quoteOpenLine = lineNumber;
            line.erase(open);
            return isEmptyOrWhitespace(line);
        }
        return false;
    }

    // Strip brace-delimited comments `{ ... }` from a line, tracking nested
    // comments across lines via `braceDepth`. Text inside comments is removed,
    // while code outside comments is preserved. Braces appearing inside
    // single-quoted string literals are ignored (treated as normal characters).
    std::string stripBraceComments(const std::string& line, int& braceDepth, ParseResult& result) {
        if (line.empty()) return line;
        
        std::string output;
        output.reserve(line.size());
        
        bool inString = false; // single-quoted EES string literal
        
        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];
            
            // Handle single-quoted strings (used for literals like 'R134a')
            if (c == '\'') {
                // Toggle string state, but only when not inside a brace comment
                if (braceDepth == 0) {
                    inString = !inString;
                    output += c;
                } else {
                    // Inside a comment, we do not need to keep string delimiters
                    // or contents; just skip them.
                }
                continue;
            }
            
            if (inString) {
                // Inside a string literal, copy characters verbatim
                if (braceDepth == 0) {
                    output += c;
                }
                continue;
            }
            
            // Outside of string literals, handle brace comments
            if (c == '{') {
                // Starting a new brace comment when not already in one:
                // count it once as a comment. Nested `{` while braceDepth>0
                // simply increase depth.
                if (braceDepth == 0) {
                    result.commentCount++;
                }
                braceDepth++;
                continue; // Do not copy '{'
            }
            
            if (c == '}') {
                if (braceDepth > 0) {
                    braceDepth--;
                    // We don't copy the closing brace either.
                    continue;
                }
                // Unmatched '}' outside a comment is treated as a normal char
                // to avoid being overly strict on slightly malformed input.
                output += c;
                continue;
            }
            
            // Only keep characters that are outside of any brace comment
            if (braceDepth == 0) {
                output += c;
            }
        }
        
        return output;
    }
    
    std::string preprocess(const std::string& source) {
        std::string result;
        result.reserve(source.size());
        
        for (size_t i = 0; i < source.size(); ++i) {
            // Convert CRLF to LF
            if (source[i] == '\r' && i + 1 < source.size() && source[i+1] == '\n') {
                result += '\n';
                ++i;
            } else if (source[i] == '\r') {
                result += '\n';
            } else {
                result += source[i];
            }
        }
        return result;
    }
    
    bool isEmptyOrWhitespace(const std::string& line) {
        return std::all_of(line.begin(), line.end(), 
            [](char c) { return std::isspace(static_cast<unsigned char>(c)); });
    }

    // Next line of a body: the text kept after an `ENDIF;` (see parseBodyLines), else the next line of
    // the stream. `resumed` is set for that text: its comments and directives are already dealt with.
    bool nextBodyLine(std::istream& stream, int& lineNumber, std::string& line, bool& resumed) {
        resumed = !pendingBodyText_.empty();
        if (resumed) {   // same line number as the ENDIF it follows
            line.swap(pendingBodyText_);
            pendingBodyText_.clear();
            return true;
        }
        if (!std::getline(stream, line)) return false;
        lineNumber++;
        return true;
    }

    static constexpr const char* kRepeatOutsideProcedureMessage =
        "REPEAT ... UNTIL loops are only allowed inside FUNCTION and PROCEDURE bodies "
        "(in the main program an iteration is solved by the solver as equations; to run a loop, "
        "put it in a FUNCTION or PROCEDURE)";

    /// Parses the REPEAT ... UNTIL(condition) loop whose REPEAT is on the current line (`line`) and appends it to
    /// `body`. Like IF ... THEN, the loop exists in the body of a FUNCTION or PROCEDURE only (EES): anywhere in
    /// the main program (a DUPLICATE body included) it is an error, with the line of the REPEAT. It used to be
    /// parsed and then silently dropped by the IR (CS-BUG-MAIN-REPEAT).
    void parseRepeatStatement(std::istream& stream, int& lineNumber, std::vector<StmtPtr>& body,
                              ParseResult& result, const std::string& line) {
        const int startLine = lineNumber;
        std::vector<StmtPtr> repBody;
        ExprPtr condition;
        std::string term = parseBodyLines(stream, lineNumber, repBody, result, &condition);
        if (term == "UNTIL" && condition) {
            body.push_back(makeRepeatUntil(repBody, condition, startLine));
        } else {
            result.errors.push_back({startLine, 0, "REPEAT missing UNTIL(condition)", line});
        }
        if (procedureBodyDepth_ == 0)
            result.errors.push_back({startLine, 0, kRepeatOutsideProcedureMessage, line});
    }

    // ========================================================================
    // Unified body parser for function/procedure/DUPLICATE/REPEAT/IF bodies.
    // Reads lines from 'stream', parsing comments, directives, IF-THEN-ELSE,
    // DUPLICATE, REPEAT-UNTIL, equations/assignments, and procedure calls.
    // Returns "END" if END keyword terminates, "UNTIL" if UNTIL(condition)
    // terminates (conditionOut set), "" on EOF, and — only when called for the
    // block of an IF statement (inIfBlock) — "ELSE" (the text of the ELSE line
    // is stored in *terminatorLine) or "ENDIF". `;` separates statements: what follows
    // `ENDIF;` on a line (another ENDIF, an ELSE, statements) is read as the next line, so
    // that `ENDIF;ENDIF;ENDIF` closes three IF blocks.
    // ========================================================================
    std::string parseBodyLines(std::istream& stream, int& lineNumber,
                               std::vector<StmtPtr>& body, ParseResult& result,
                               ExprPtr* conditionOut = nullptr,
                               bool inIfBlock = false,
                               std::string* terminatorLine = nullptr) {
        std::string bodyLine;
        bool resumed;
        while (nextBodyLine(stream, lineNumber, bodyLine, resumed)) {
            if (!resumed && skipCommentLine(bodyLine, lineNumber, result)) continue;
            std::string trimmedBody = trim(bodyLine);
            if (!resumed && handleConditionalDirective(trimmedBody, lineNumber, bodyLine, result)) continue;
            if (!conditionalActive()) continue;
            std::string upperBody = trimmedBody;
            std::transform(upperBody.begin(), upperBody.end(), upperBody.begin(), ::toupper);

            auto isBodyKeyword = [&](const std::string& kw) {
                if (upperBody.find(kw) != 0) return false;
                if (upperBody.size() == kw.size()) return true;
                return !isIdentChar(upperBody[kw.size()]);
            };

            // END, possibly followed by a comment ("END // end of the function")
            std::string codeBody = trim(removeInlineComments(bodyLine));
            std::transform(codeBody.begin(), codeBody.end(), codeBody.begin(), ::toupper);
            if (codeBody == "END") return "END";

            // Check for UNTIL(condition)
            if (isBodyKeyword("UNTIL")) {
                if (conditionOut) {
                    size_t parenOpen = trimmedBody.find('(');
                    size_t parenClose = trimmedBody.rfind(')');
                    if (parenOpen != std::string::npos && parenClose != std::string::npos && parenClose > parenOpen) {
                        std::string condStr = trim(trimmedBody.substr(parenOpen + 1, parenClose - parenOpen - 1));
                        *conditionOut = parseExpression(condStr, lineNumber);
                    }
                }
                return "UNTIL";
            }

            if (auto stmt = tryParseComment(bodyLine, lineNumber)) {
                body.push_back(stmt);
            } else if (auto stmt = tryParseDirective(bodyLine, lineNumber)) {
                body.push_back(stmt);
                noteConstantDirective(stmt->as<Directive>(), lineNumber, bodyLine, result);
                if (procedureBodyDepth_ > 0 && isCommonDirective(stmt->as<Directive>()))
                    result.errors.push_back({lineNumber, 0, commonDirectiveMessage(stmt->as<Directive>().content), bodyLine});
            } else if (isBodyKeyword("IF")) {
                const int ifLine = lineNumber;
                std::string term = parseIfStatement(trim(removeInlineComments(bodyLine)),
                                                    stream, lineNumber, body, result, bodyLine);
                if (procedureBodyDepth_ == 0)  // e.g. inside a DUPLICATE of the main program
                    result.errors.push_back({ifLine, 0, kIfOutsideProcedureMessage, bodyLine});
                if (!term.empty()) return term;  // block ended prematurely (END/UNTIL)
            } else if (isBodyKeyword("ELSE") || isBodyKeyword("ENDIF")) {
                const bool isElse = isBodyKeyword("ELSE");
                if (!isElse) {  // `ENDIF;text`: the text is read as the next line
                    const std::string code = trim(removeInlineComments(bodyLine));
                    const size_t semicolon = code.find(';');
                    const std::string stray = trim(code.substr(5, semicolon == std::string::npos ? semicolon
                                                                                               : semicolon - 5));
                    if (!stray.empty())
                        result.errors.push_back({lineNumber, 0, "Unexpected text '" + stray +
                                                 "' after ENDIF (separate statements with ';')", bodyLine});
                    const size_t next = semicolon == std::string::npos ? semicolon
                                                                        : code.find_first_not_of("; \t", semicolon);
                    if (next != std::string::npos) pendingBodyText_ = code.substr(next);
                }
                if (!inIfBlock) {
                    if (!isElse && procedureBodyDepth_ > 0) {
                        // EES ignores an ENDIF that closes nothing (the ladders of its libraries have some)
                        if (lastIgnoredEndifLine_ != lineNumber)
                            result.diagnostics.push(DiagnosticSeverity::Warning, "P008",
                                "ENDIF without a matching IF ... THEN: ignored (EES accepts it)",
                                "parser", lineNumber, 0);
                        lastIgnoredEndifLine_ = lineNumber;
                    } else {
                        result.errors.push_back({lineNumber, 0,
                            std::string(isElse ? "ELSE" : "ENDIF") + " without a matching IF ... THEN", bodyLine});
                    }
                    continue;
                }
                if (terminatorLine) *terminatorLine = trimmedBody;
                return isElse ? "ELSE" : "ENDIF";
            } else if (isBodyKeyword("DUPLICATE")) {
                // Parse nested DUPLICATE: DUPLICATE varName = startExpr, endExpr
                std::string header = trim(trimmedBody.substr(9));
                auto eqPos = header.find('=');
                if (eqPos == std::string::npos) {
                    result.errors.push_back({lineNumber, 0, "DUPLICATE missing '=' in header", bodyLine});
                    continue;
                }
                std::string iterVar = trim(header.substr(0, eqPos));
                std::string rangeStr = trim(header.substr(eqPos + 1));
                auto commaPos = rangeStr.find(',');
                if (commaPos == std::string::npos) {
                    result.errors.push_back({lineNumber, 0, "DUPLICATE missing ',' between start and end", bodyLine});
                    continue;
                }
                std::string startStr = trim(rangeStr.substr(0, commaPos));
                std::string endStr = trim(rangeStr.substr(commaPos + 1));
                auto startExpr = parseExpression(startStr, lineNumber);
                auto endExpr = parseExpression(endStr, lineNumber);
                if (!startExpr || !endExpr) {
                    result.errors.push_back({lineNumber, 0, "DUPLICATE could not parse range expressions", bodyLine});
                    continue;
                }
                int dupStartLine = lineNumber;
                std::vector<StmtPtr> dupBody;
                std::string term = parseBodyLines(stream, lineNumber, dupBody, result);
                if (term == "END") {
                    body.push_back(makeDuplicate(iterVar, startExpr, endExpr, dupBody, dupStartLine));
                    duplicateSeen_ = true;
                } else {
                    result.errors.push_back({dupStartLine, 0, "DUPLICATE missing END", bodyLine});
                }
            } else if (isBodyKeyword("REPEAT")) {
                parseRepeatStatement(stream, lineNumber, body, result, bodyLine);
            } else if (std::isdigit(static_cast<unsigned char>(trimmedBody.empty() ? ' ' : trimmedBody[0])) &&
                       parseLabelLine(bodyLine, lineNumber, body, result)) {
                // a statement label (CS-GAP-GOTO)
            } else {
                parseSimpleStatements(bodyLine, lineNumber, body, result, bodyLine);
            }
        }
        return "";  // EOF
    }

    // ------------------------------------------------------------------------
    // GOTO and numeric statement labels (CS-GAP-GOTO)
    //
    //   GOTO 100      the body goes on after the label `100:` of the same body
    //   100:          a label alone on its line, or followed by a statement (`90 : j = j - 1`)
    //
    // EES has them in the body of a FUNCTION or PROCEDURE only (the statements of the main program are
    // equations, solved in no particular order): elsewhere they are errors. GOTO is a statement like
    // another: a line, a segment after `;`, a branch of a single-line IF (`IF (c) THEN GOTO 10 ELSE GOTO 20`)
    // or of a block IF. The label must be in the body of the GOTO, in its block or in a block that contains
    // it: a jump into a block is refused, like a missing or a duplicated label (checkBodyJumps). The
    // evaluator executes the jumps (ProceduralEvaluator).
    // ------------------------------------------------------------------------
    static constexpr const char* kGotoOutsideRoutineMessage =
        "GOTO is only allowed inside FUNCTION and PROCEDURE bodies "
        "(the statements of the main program are equations, solved in no particular order)";
    static constexpr const char* kLabelOutsideRoutineMessage =
        "statement labels are only allowed inside FUNCTION and PROCEDURE bodies "
        "(the statements of the main program are equations, solved in no particular order)";

    /// The statement label at the start of `code` (comments removed, trimmed): `100:`, `10 :`, `90 : j = j - 1`.
    /// Sets its number and the text after the colon; false when `code` does not start with a label
    /// (`n := 3` is an assignment).
    bool splitStatementLabel(const std::string& code, int& number, std::string& rest) {
        size_t digits = 0;
        while (digits < code.size() && std::isdigit(static_cast<unsigned char>(code[digits]))) ++digits;
        if (digits == 0 || digits > 9) return false;
        const size_t colon = code.find_first_not_of(" \t", digits);
        if (colon == std::string::npos || code[colon] != ':' || (colon + 1 < code.size() && code[colon + 1] == '='))
            return false;
        number = std::stoi(code.substr(0, digits));
        rest = trim(code.substr(colon + 1));
        return true;
    }

    /// Parses a body line that starts with a statement label: appends the label to `body`, and the text after it
    /// is read as the next line. False (nothing done) when the line starts with no label.
    bool parseLabelLine(const std::string& bodyLine, int lineNumber, std::vector<StmtPtr>& body, ParseResult& result) {
        int number = 0;
        std::string rest;
        if (!splitStatementLabel(trim(removeInlineComments(bodyLine)), number, rest)) return false;
        if (procedureBodyDepth_ == 0)  // e.g. inside a DUPLICATE of the main program
            result.errors.push_back({lineNumber, 0, kLabelOutsideRoutineMessage, bodyLine});
        else
            body.push_back(makeLabel(number, lineNumber));
        pendingBodyText_ = rest;
        return true;
    }

    /// Parses `segment` as `GOTO n` and appends it to `body`; false when it is no GOTO (the keyword followed by
    /// `=` is a variable named goto). A GOTO that is refused (no label number, main program) is reported here.
    bool parseGotoSegment(const std::string& segment, int lineNumber, std::vector<StmtPtr>& body,
                          ParseResult& result, const std::string& errorLine) {
        const std::string text = trim(segment);
        if (text.empty() || (text[0] != 'g' && text[0] != 'G')) return false;
        const std::string code = trim(removeInlineComments(text));
        if (findKeywordAtDepth0(code, "GOTO") != 0) return false;
        const std::string target = trim(code.substr(4));
        if (!target.empty() && (target[0] == '=' || target[0] == ':')) return false;
        if (procedureBodyDepth_ == 0) {
            result.errors.push_back({lineNumber, 0, kGotoOutsideRoutineMessage, errorLine});
        } else if (target.empty() || target.size() > 9 ||
                   !std::all_of(target.begin(), target.end(), [](unsigned char c) { return std::isdigit(c) != 0; })) {
            result.errors.push_back({lineNumber, 0, "GOTO must be followed by a statement label number (GOTO 100), "
                                     "found '" + target + "'", errorLine});
        } else {
            body.push_back(makeGoto(std::stoi(target), lineNumber));
        }
        return true;
    }

    // Parse the equations/assignments and procedure calls of `text` (one line
    // of a body, or one branch of a single-line IF), split on semicolons.
    // `errorLine` is the source line reported with parse errors. A statement that
    // cannot be parsed is an error, as in the main program: it is never dropped
    // (CS-BUG-BODY-LINE-DROP); only a segment that holds nothing but comments (or
    // separators) is skipped.
    void parseSimpleStatements(const std::string& text, int lineNumber,
                               std::vector<StmtPtr>& body, ParseResult& result,
                               const std::string& errorLine) {
        const auto segments = splitOnSemicolons(text);
        for (const auto& seg : segments) {
            if (parseGotoSegment(seg, lineNumber, body, result, errorLine)) continue;
            std::string segErr;
            StmtPtr stmt = tryParseEquationOrAssignment(seg, lineNumber, &segErr);
            if (!stmt) stmt = tryParseProcedureCall(seg, lineNumber, &segErr);
            if (stmt) {
                body.push_back(stmt);
            } else if (!segErr.empty()) {
                result.errors.push_back({lineNumber, 0, segErr, errorLine});
            } else if (!trim(removeInlineComments(seg)).empty()) {
                result.errors.push_back({lineNumber, 0, segments.size() > 1
                                         ? "Could not parse segment: " + seg : "Could not parse line", errorLine});
            }
        }
    }

    // ------------------------------------------------------------------------
    // A unit annotation after a variable (CS-BUG-VAR-ANNOTATION)
    //
    // `y = a [kW]` has the form of the array element `a[kW]`, so the unit was read as its index: a variable
    // `kW` that nothing sets, a structural error of the main program or the read of an unassigned local in a
    // body, neither saying why. An element whose only index is a bare identifier that names a unit
    // (isUnitName) is therefore an error, with the line, unless that identifier is a valid index in its scope:
    // the index of an enclosing DUPLICATE or indexed SUM, or a variable of the scope (main program: one of its
    // equations uses it, see walkExpression; body: an input or a local that the body assigns, see
    // checkBodyReads). What EES does with a unit after a variable is not established. The unit after a number,
    // a parenthesised group, a call or an element is accepted (parsePrimary, CS-BUG-TRAILING-TEXT) and is no
    // element of that form.
    // ------------------------------------------------------------------------
    /// An element `variable[unit]` whose index is a bare unit name, at `line`.
    struct UnitSubscript {
        int line;
        std::string variable;
        std::string unit;
    };

    /// The bare unit name that is the only index of `var` (`a[kW]`); null for any other variable.
    static const std::string* unitIndexOf(const Variable& var) {
        if (var.indices.size() != 1 || !var.indices[0]->is<Variable>()) return nullptr;
        const Variable& index = var.indices[0]->as<Variable>();
        return index.indices.empty() && isUnitName(index.name) ? &index.name : nullptr;
    }

    static std::string unitAnnotationMessage(const UnitSubscript& element) {
        return "'" + element.variable + " [" + element.unit + "]': a unit annotation after a variable is not "
               "supported ('" + element.unit + "' would be read as an array index, and nothing sets '" + element.unit +
               "' here); put the unit in the comment of the equation (\"[" + element.unit + "]\") or after a number";
    }

    // ------------------------------------------------------------------------
    // Check of a FUNCTION / PROCEDURE body (CS-BUG-IMPLICIT-PROC, CS-BUG-LOCAL-UNASSIGNED,
    // CS-GAP-STRING-LOCAL-UNSET; maintainer decision of 2026-10-12)
    //
    // A body is a sequence of assignments executed in order, as in EES (`=` and `:=` alike), and its
    // variables other than the inputs are local and not initialised (EES help, Internal Functions).
    // Two errors, reported with their line like the parse errors so that the model is not solved
    // (the evaluator would skip the statement, or read the variable as 1):
    //  - an implicit statement: the left-hand side is not a variable or an array element. A
    //    lookup-table cell write `LOOKUP(...) = ...` is left alone (the open item CS-BUG-LOOKUP-WRITE);
    //  - a read of a variable that is no input, no built-in constant (pi, g#...) and that no earlier
    //    statement of the body assigns, in text order. Flow-insensitive: an assignment in either
    //    branch of an IF counts, an array counts as assigned once one of its elements is; inside a
    //    REPEAT loop a variable read before its assignment in the loop is an error. The variables
    //    read are those the equations of the main program count (IR::extractVariables);
    //  - an output (the FUNCTION's own name, a PROCEDURE output) that no statement assigns, reported
    //    on the routine's header line (flow-insensitive too: one branch or a CALL output suffices).
    // Each variable is reported once, at its first read. Called once per routine, at parse time,
    // and only when its body parsed without error: a statement that failed to parse would make the
    // variables it assigns look unassigned.
    // ------------------------------------------------------------------------
    struct BodyCheck {
        std::string routine;                                  // "FUNCTION 'f'" / "PROCEDURE 'p'"
        std::set<std::string, CaseInsensitiveLess> assigned;  // the inputs, then every name assigned so far
        std::set<std::string, CaseInsensitiveLess> reported;  // variables already reported as unassigned
        const std::string& source;                            // the preprocessed source (error context)
        ParseResult& result;
        std::set<std::string, CaseInsensitiveLess> constants; // the $CONSTANTs (global, above or below)
        std::vector<UnitSubscript> unitSubscripts;            // the elements a[kW] whose index is not assigned (yet)
    };

    void checkRoutineBody(const std::vector<std::string>& inputs, const std::vector<std::string>& outputs,
                          const std::vector<StmtPtr>& body, int headerLine, const std::string& source,
                          ParseResult& result) {
        BodyCheck check{routineKind_ + " '" + routineName_ + "'", {}, {}, source, result, presetConstants_, {}};
        for (const auto& entry : userConstants_) check.constants.insert(entry.first);
        routineChecked_ = true;
        for (const auto& input : inputs) check.assigned.insert(arrayBaseName(input));   // an array input `SP[1..5]`: SP
        checkBodyJumps(body, check);
        checkBodyStatements(body, check);
        // `a [kW]` is the element a[kW] read before `kW` is assigned, unless the body assigns `kW` further down
        for (const auto& element : check.unitSubscripts) {
            if (!check.assigned.count(element.unit))
                bodyError(check, element.line, unitAnnotationMessage(element));
            else if (check.reported.insert(element.unit).second)
                bodyError(check, element.line, unassignedReadMessage(element.unit, check));
        }
        for (const auto& output : outputs)
            if (!check.assigned.count(arrayBaseName(output)))
                bodyError(check, headerLine, check.routine + " never assigns its output '" + arrayBaseName(output) +
                          "': every output of a FUNCTION or PROCEDURE must be assigned by its body");
    }

    // ------------------------------------------------------------------------
    // The GOTOs and labels of a body (CS-GAP-GOTO), reported with their line like the other errors of a body:
    //  - a GOTO whose label the body does not have;
    //  - a label defined twice in the body (the GOTO could not tell which one it means);
    //  - a jump into a block: the label is in an IF branch or a loop body that does not contain the GOTO.
    //    Jumps within a block, backward jumps (loops) and jumps out of blocks are valid.
    // The reads and assignments of the body are checked in text order whatever the jumps, as for an IF.
    // ------------------------------------------------------------------------
    struct JumpScan {
        using Block = const std::vector<StmtPtr>*;
        /// The statement of which a block is a branch or the body (null for the body of the routine).
        struct Owner {
            const Statement* statement;
            bool isElse;
        };
        struct LabelAt {
            int line;
            Block block;   // the block that holds the label
            Owner owner;
        };
        struct GotoAt {
            int line;
            int label;
            std::vector<Block> enclosing;   // the blocks that hold the GOTO, outermost first
        };
        std::map<int, LabelAt> labels;
        std::vector<GotoAt> gotos;
    };

    static void scanJumps(const std::vector<StmtPtr>& block, JumpScan::Owner owner,
                          std::vector<JumpScan::Block>& enclosing, JumpScan& scan, BodyCheck& check) {
        enclosing.push_back(&block);
        for (const auto& stmt : block) {
            const int line = stmt->sourceLineNumber;
            if (stmt->is<Label>()) {
                const int number = stmt->as<Label>().number;
                const auto added = scan.labels.emplace(number, JumpScan::LabelAt{line, &block, owner});
                if (!added.second)
                    bodyError(check, line, "label " + std::to_string(number) + " is defined twice in " + check.routine +
                              " (lines " + std::to_string(added.first->second.line) + " and " + std::to_string(line) +
                              "): a GOTO could not tell which one it means");
            } else if (stmt->is<Goto>()) {
                scan.gotos.push_back({line, stmt->as<Goto>().label, enclosing});
            } else if (stmt->is<IfThenElse>()) {
                const auto& ite = stmt->as<IfThenElse>();
                scanJumps(ite.thenBranch, {stmt.get(), false}, enclosing, scan, check);
                scanJumps(ite.elseBranch, {stmt.get(), true}, enclosing, scan, check);
            } else if (stmt->is<Duplicate>()) {
                scanJumps(stmt->as<Duplicate>().body, {stmt.get(), false}, enclosing, scan, check);
            } else if (stmt->is<RepeatUntil>()) {
                scanJumps(stmt->as<RepeatUntil>().body, {stmt.get(), false}, enclosing, scan, check);
            }
        }
        enclosing.pop_back();
    }

    /// "the THEN branch of the IF at line 4", "the REPEAT loop at line 7": the block a label is in.
    static std::string describeJumpBlock(const JumpScan::Owner& owner) {
        const std::string at = " at line " + std::to_string(owner.statement->sourceLineNumber);
        if (owner.statement->is<IfThenElse>())
            return std::string("the ") + (owner.isElse ? "ELSE" : "THEN") + " branch of the IF" + at;
        return std::string("the ") + (owner.statement->is<Duplicate>() ? "DUPLICATE" : "REPEAT") + " loop" + at;
    }

    static void checkBodyJumps(const std::vector<StmtPtr>& body, BodyCheck& check) {
        JumpScan scan;
        std::vector<JumpScan::Block> enclosing;
        scanJumps(body, {nullptr, false}, enclosing, scan, check);
        for (const auto& jump : scan.gotos) {
            const std::string number = std::to_string(jump.label);
            const auto target = scan.labels.find(jump.label);
            if (target == scan.labels.end()) {
                bodyError(check, jump.line, "GOTO " + number + " has no target: " + check.routine +
                          " has no statement label '" + number + ":'");
            } else if (std::find(jump.enclosing.begin(), jump.enclosing.end(), target->second.block) == jump.enclosing.end()) {
                bodyError(check, jump.line, "GOTO " + number + " would jump into a block: the label '" + number + ":' (line " +
                          std::to_string(target->second.line) + ") is inside " + describeJumpBlock(target->second.owner) +
                          "; a GOTO can only jump to a label of its own block or of a block that contains it");
            }
        }
    }

    static void checkBodyStatements(const std::vector<StmtPtr>& statements, BodyCheck& check) {
        for (const auto& stmt : statements) {
            const int line = stmt->sourceLineNumber;
            if (stmt->is<Equation>() || stmt->is<Assignment>()) {
                const bool isEquation = stmt->is<Equation>();
                const ExprPtr& lhs = isEquation ? stmt->as<Equation>().lhs : stmt->as<Assignment>().lhs;
                const ExprPtr& rhs = isEquation ? stmt->as<Equation>().rhs : stmt->as<Assignment>().rhs;
                if (isLookupCellWrite(lhs)) continue;
                if (!lhs->is<Variable>()) {
                    bodyError(check, line, "implicit equation in " + check.routine + ": the statements of a FUNCTION "
                              "or PROCEDURE are assignments executed in order, so the left-hand side must be a variable. "
                              "Solve the equation in the main program, or iterate explicitly in the body (REPEAT ... UNTIL)");
                    markBodyAssigned(lhs, check);   // its variables are not reported again at each later read
                    continue;
                }
                checkBodyReads(rhs, line, check);
                for (const auto& index : lhs->as<Variable>().indices) checkBodyReads(index, line, check);
                check.assigned.insert(lhs->as<Variable>().name);
            } else if (stmt->is<IfThenElse>()) {
                const auto& ite = stmt->as<IfThenElse>();
                checkBodyReads(ite.condition, line, check);
                checkBodyStatements(ite.thenBranch, check);
                checkBodyStatements(ite.elseBranch, check);
            } else if (stmt->is<ProcedureCall>()) {
                const auto& call = stmt->as<ProcedureCall>();
                for (const auto& arg : call.inputArgs) checkBodyReads(arg, line, check);
                for (const auto& out : call.outputVars) {
                    for (const auto& index : out.indices) checkBodyReads(index, line, check);
                    check.assigned.insert(out.name);
                }
            } else if (stmt->is<Duplicate>()) {
                const auto& dup = stmt->as<Duplicate>();
                checkBodyReads(dup.start, line, check);
                checkBodyReads(dup.end, line, check);
                check.assigned.insert(dup.iteratorVar);
                checkBodyStatements(dup.body, check);
            } else if (stmt->is<RepeatUntil>()) {
                const auto& rep = stmt->as<RepeatUntil>();
                checkBodyStatements(rep.body, check);   // the condition is tested after the loop body
                checkBodyReads(rep.condition, rep.condition->sourceLineNumber > 0 ? rep.condition->sourceLineNumber
                                                                                   : line, check);
            }
        }
    }

    // `LOOKUP('table', row, col) = value` / `LOOKUP$(...) = 'text'`: a lookup-table cell write
    static bool isLookupCellWrite(const ExprPtr& lhs) {
        if (!lhs->is<FunctionCall>()) return false;
        std::string name = lhs->as<FunctionCall>().name;
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        return name == "lookup" || name == "lookup$";
    }

    // The name of an array without the indices that IR::extractVariables flattens into it (x[i] -> x).
    static std::string arrayBaseName(const std::string& flattened) {
        return flattened.substr(0, flattened.find('['));
    }

    static void checkBodyReads(const ExprPtr& expr, int line, BodyCheck& check) {
        std::set<std::string, CaseInsensitiveLess> read, annotated;
        IR::extractVariables(expr, read);
        const auto unset = [&](const std::string& name) {
            return !check.assigned.count(name) && !Constants::isConstant(name) && !check.constants.count(name);
        };
        for (const auto& flattened : read) {   // `a [kW]` reads the element a[kW]: kW is an index (or the unit)
            const size_t open = flattened.find('[');
            if (open == std::string::npos || flattened.back() != ']') continue;
            const std::string unit = flattened.substr(open + 1, flattened.size() - open - 2);
            if (!isUnitName(unit) || !unset(unit)) continue;
            check.unitSubscripts.push_back({line, flattened.substr(0, open), unit});   // reported once the body is read
            annotated.insert(unit);
        }
        for (const auto& flattened : read) {
            const std::string name = arrayBaseName(flattened);
            if (unset(name) && !annotated.count(name) && check.reported.insert(name).second)
                bodyError(check, line, unassignedReadMessage(name, check));
        }
    }

    static std::string unassignedReadMessage(const std::string& name, const BodyCheck& check) {
        return "variable '" + name + "' is read in " + check.routine + " before it is assigned: "
               "the variables of a FUNCTION or PROCEDURE other than its inputs are local and not initialised "
               "(pass main-program values as input arguments)";
    }

    static void markBodyAssigned(const ExprPtr& expr, BodyCheck& check) {
        std::set<std::string, CaseInsensitiveLess> names;
        IR::extractVariables(expr, names);
        for (const auto& flattened : names) check.assigned.insert(arrayBaseName(flattened));
    }

    static void bodyError(BodyCheck& check, int line, const std::string& message) {
        // the source line as context, like the other parse errors
        check.result.errors.push_back({line, 0, message, sourceLineText(check.source, line)});
    }

    /// Line `line` (1-based) of `source`; empty when there is none.
    static std::string sourceLineText(const std::string& source, int line) {
        size_t start = 0;
        for (int i = 1; i < line && start != std::string::npos; ++i) {
            start = source.find('\n', start);
            if (start != std::string::npos) ++start;
        }
        return start == std::string::npos ? std::string() : source.substr(start, source.find('\n', start) - start);
    }

    // ------------------------------------------------------------------------
    // $CONSTANT Name# = ValueSI [UnitsSI] ValueEng [UnitsEng] (EES help, $Constant)
    //
    // A value known before the solve and global: the main program and every FUNCTION / PROCEDURE see
    // it, above or below the directive. ValueSI is a number or an expression without spaces of numbers
    // and of the constants defined above; the units and the English value are ignored (CoolSolve works
    // in SI). Its reads are replaced by the value once the file is parsed (resolveCompileTimeValues).
    // A string constant (`Name#$ = 'text'`) is not supported.
    // ------------------------------------------------------------------------
    void noteConstantDirective(const Directive& directive, int lineNumber, const std::string& sourceLine,
                               ParseResult& result) {
        std::string kind = directive.name;
        std::transform(kind.begin(), kind.end(), kind.begin(), ::tolower);
        if (kind != "constant") return;
        auto fail = [&](const std::string& message) { result.errors.push_back({lineNumber, 0, message, sourceLine}); };
        const std::string text = trim(removeInlineComments(directive.content));
        const size_t eq = text.find('=');
        const std::string name = trim(text.substr(0, eq));
        if (eq == std::string::npos || name.empty()) return fail("$CONSTANT: expected $CONSTANT Name# = value [units]");
        if (name.size() > 2 && name.compare(name.size() - 2, 2, "#$") == 0)
            return fail("$CONSTANT " + name + ": string constants are not supported (set a string variable in the "
                        "main program and pass it to the routines as an argument)");
        if (name.back() != '#' || !isIdentifier(name))
            return fail("$CONSTANT: the name of a constant ends with # ('" + name + "')");
        const std::string rest = trim(text.substr(eq + 1));
        const size_t valueEnd = rest.find_first_of(" \t[");
        const std::string valueText = rest.substr(0, valueEnd);
        const std::string after = valueEnd == std::string::npos ? "" : trim(rest.substr(valueEnd));
        if (valueText.empty() || (!after.empty() && after[0] != '['))
            return fail("$CONSTANT " + name + ": the value is a number or an expression without spaces, followed "
                        "by its [units] if any (EES), not '" + rest + "'");
        const auto existing = userConstants_.find(name);
        if (existing != userConstants_.end())
            return fail("$CONSTANT " + name + " is already defined (line " + std::to_string(existing->second.line) + ")");
        clearExpressionError();
        std::string unknown;
        const auto value = evaluateConstantExpression(parseExpression(valueText, lineNumber),
            [this](const std::string& n) { return constantValue(n); }, &unknown);
        clearExpressionError();
        if (!value)
            return fail("$CONSTANT " + name + ": " + (unknown.empty() ? "'" + valueText + "' is not a constant expression"
                                                                      : "'" + unknown + "' is not a constant defined above"));
        userConstants_[name] = {*value, lineNumber};
        if (routineChecked_ && !presetConstants_.count(name)) constantAfterRoutine_ = true;
    }

    /// Value of a $CONSTANT defined so far or of a built-in constant (pi, g#...); nullopt otherwise.
    std::optional<double> constantValue(const std::string& name) const {
        const auto it = userConstants_.find(name);
        if (it != userConstants_.end()) return it->second.value;
        if (Constants::isConstant(name)) return Constants::getValue(name);
        return std::nullopt;
    }

    // ------------------------------------------------------------------------
    // Values known before the solve (CS-BUG-DUPLICATE-VAR-BOUND, CS-GAP-ARRAY-ARG)
    //
    // One walk of the program, once the whole file is parsed without error (a variable may be set
    // below the loop that uses it), and only when the file has a $CONSTANT, a DUPLICATE or an array range:
    //  - every read of a $CONSTANT becomes its value (a routine that declares the name in its header
    //    keeps it as its own variable); assigning a constant is an error;
    //  - the bounds of the DUPLICATE loops of the main program and of the array ranges X[a..b] are
    //    resolved by evaluateConstantExpression, a name standing for a constant or for a variable of the
    //    main program set to a constant expression by an equation `n = <expression>`, possibly through
    //    other such variables (n = 2*m, m = 3). A DUPLICATE bound that uses the index of an enclosing
    //    DUPLICATE stays an expression of it, evaluated by the IR for each copy; any other bound becomes
    //    its integer value. In a FUNCTION / PROCEDURE (header and body) the bounds of a range are numbers
    //    or constants: a body sees no variable of the main program and its own are known only when it
    //    runs (the DUPLICATE loops of a body are run by the evaluator and keep their bounds);
    //  - an indexed aggregation SUM(expr, i=a, b) / PRODUCT(expr, i=a, b) becomes the variadic SUM / PRODUCT of
    //    the copies of expr with i = a, ..., b (the index is local to it), its bounds being resolved as those
    //    of a DUPLICATE; one that uses the index of an enclosing DUPLICATE is expanded by the IR in each copy
    //    (resolveIndexedAggregate);
    //  - a range stands for its elements, in order (the last index varying fastest), in the header of a
    //    routine and in the arguments of a function or of a CALL, outputs included; anywhere else it is
    //    an error. A call that passes a range, or calls a routine whose header has one, must give as many
    //    values as the header declares (checked here, with the line, rather than at the solve).
    // Each error names the bound and is reported with its line, like a parse error.
    // ------------------------------------------------------------------------
    using NameSet = std::set<std::string, CaseInsensitiveLess>;

    struct RoutineSignature {
        std::string kind;   // "FUNCTION" / "PROCEDURE"
        std::string name;
        int line = 0;
        size_t inputs = 0, outputs = 0;   // after the ranges of the header are expanded
        bool hasRange = false;            // the header has a range
    };

    struct CompileTimeContext {
        ParseResult& result;
        const std::string& source;
        /// The variables of the main program set by an equation `name = <expression>` (the first one).
        std::map<std::string, std::pair<ExprPtr, int>, CaseInsensitiveLess> definitions;
        std::map<std::string, std::optional<double>, CaseInsensitiveLess> values;   // resolved so far
        NameSet resolving;                                                         // cycle guard
        std::map<std::string, RoutineSignature, CaseInsensitiveLess> routines;
        NameSet iterators;     // the indices of the enclosing DUPLICATE loops (main program)
        std::string routine;   // "FUNCTION 'f'" while a routine is walked, empty in the main program
        NameSet locals;        // the names its header declares (they hide a $CONSTANT)
        int line = 0;          // line of the statement being walked
        std::vector<UnitSubscript> unitSubscripts;   // the elements a[kW] of the main program (CS-BUG-VAR-ANNOTATION)
        NameSet used;          // the variables of the main program that the walk met (filled when unitIndexSeen_)

        void error(const std::string& message) {
            result.errors.push_back({line, 0, message, sourceLineText(source, line)});
        }
    };

    static constexpr const char* kBoundRule =
        ". A bound is a number, a $CONSTANT, or a variable set to a constant expression by an equation of the "
        "main program (n = 17), combined with + - * / ^";
    static constexpr size_t kMaxRangeElements = 100000;

    static bool isRangeIndex(const ExprPtr& index) {
        return index && index->is<BinaryOp>() && index->as<BinaryOp>().op == "..";
    }
    static bool hasRange(const Variable& var) {
        return std::any_of(var.indices.begin(), var.indices.end(), isRangeIndex);
    }

    static std::string formatNumber(double value) {
        std::ostringstream os;
        os << std::setprecision(15) << value;
        return os.str();
    }

    /// The expression as written in a message: `x[1..n]`, `y + 1`, `2*m`.
    static std::string formatExpr(const ExprPtr& expr) {
        if (!expr) return "";
        auto precedence = [](const std::string& op) {
            return op == "+" || op == "-" ? 1 : op == "*" || op == "/" ? 2 : op == "^" ? 3 : 0;
        };
        auto operand = [&](const ExprPtr& e, int parent, bool right) {
            const std::string text = formatExpr(e);
            if (!e->is<BinaryOp>() || e->as<BinaryOp>().op == "..") return text;
            const int own = precedence(e->as<BinaryOp>().op);
            return own < parent || (right && own == parent) ? "(" + text + ")" : text;
        };
        if (expr->is<NumberLiteral>()) return formatNumber(expr->as<NumberLiteral>().value);
        if (expr->is<StringLiteral>()) return "'" + expr->as<StringLiteral>().value + "'";
        if (expr->is<Variable>()) {
            const auto& var = expr->as<Variable>();
            std::string text = var.name;
            for (size_t i = 0; i < var.indices.size(); ++i)
                text += (i == 0 ? "[" : ",") + formatExpr(var.indices[i]) + (i + 1 == var.indices.size() ? "]" : "");
            return text;
        }
        if (expr->is<UnaryOp>()) return expr->as<UnaryOp>().op + operand(expr->as<UnaryOp>().operand, 3, false);
        if (expr->is<BinaryOp>()) {
            const auto& op = expr->as<BinaryOp>();
            if (op.op == "..") return formatExpr(op.left) + ".." + formatExpr(op.right);
            const int own = precedence(op.op);
            const std::string sep = own == 1 || own == 0 ? " " : "";
            return operand(op.left, own, false) + sep + op.op + sep + operand(op.right, own, true);
        }
        const auto& call = expr->as<FunctionCall>();
        std::string text = call.name + "(";
        for (size_t i = 0; i < call.args.size(); ++i) text += (i ? ", " : "") + formatExpr(call.args[i]);
        for (size_t i = 0; i < call.namedArgs.size(); ++i)
            text += (i || !call.args.empty() ? ", " : "") + call.namedArgs[i].first + "=" + formatExpr(call.namedArgs[i].second);
        return text + ")";
    }

    void resolveCompileTimeValues(ParseResult& result, const std::string& source) {
        auto& statements = result.program.statements;
        CompileTimeContext ctx{result, source, {}, {}, {}, {}, {}, {}, {}, 0};
        for (const auto& stmt : statements) {
            const ExprPtr* lhs = nullptr;
            const ExprPtr* rhs = nullptr;
            if (stmt->is<Equation>()) { lhs = &stmt->as<Equation>().lhs; rhs = &stmt->as<Equation>().rhs; }
            if (stmt->is<Assignment>()) { lhs = &stmt->as<Assignment>().lhs; rhs = &stmt->as<Assignment>().rhs; }
            if (lhs && (*lhs)->is<Variable>() && (*lhs)->as<Variable>().indices.empty())
                ctx.definitions.emplace((*lhs)->as<Variable>().name, std::make_pair(*rhs, stmt->sourceLineNumber));
        }
        for (auto& stmt : statements) {   // the headers first: the calls are checked against them
            ctx.line = stmt->sourceLineNumber;
            RoutineSignature signature;
            signature.line = stmt->sourceLineNumber;
            const size_t errorsBefore = result.errors.size();
            if (stmt->is<FunctionDefinition>()) {
                auto& def = stmt->as<FunctionDefinition>();
                signature.kind = "FUNCTION";
                signature.name = def.name;
                ctx.routine = "FUNCTION '" + def.name + "'";
                expandHeader(def.parameters, signature.hasRange, ctx);
                signature.inputs = def.parameters.size();
                signature.outputs = 1;
            } else if (stmt->is<ProcedureDefinition>()) {
                auto& def = stmt->as<ProcedureDefinition>();
                signature.kind = "PROCEDURE";
                signature.name = def.name;
                ctx.routine = "PROCEDURE '" + def.name + "'";
                expandHeader(def.inputs, signature.hasRange, ctx);
                expandHeader(def.outputs, signature.hasRange, ctx);
                signature.inputs = def.inputs.size();
                signature.outputs = def.outputs.size();
            } else {
                continue;
            }
            if (result.errors.size() == errorsBefore)   // a header in error is not checked against its calls
                ctx.routines.emplace(signature.name, signature);
        }
        ctx.routine.clear();
        walkStatements(statements, ctx);
        for (const auto& element : ctx.unitSubscripts)   // `a [kW]`, unless the main program uses a variable `kW`
            if (!ctx.used.count(element.unit)) {
                ctx.line = element.line;
                ctx.error(unitAnnotationMessage(element));
            }
    }

    // --- bounds ---------------------------------------------------------------

    /// Value of `name` in a bound or a constant expression: a constant (not hidden by a name of the
    /// routine's header) and, in the main program, a variable set to a constant expression.
    std::optional<double> knownValue(CompileTimeContext& ctx, const std::string& name) {
        if (ctx.locals.count(name)) return std::nullopt;
        if (!ctx.definitions.count(name) || userConstants_.count(name))
            if (const auto constant = constantValue(name)) return constant;
        if (!ctx.routine.empty()) return std::nullopt;
        const auto known = ctx.values.find(name);
        if (known != ctx.values.end()) return known->second;
        const auto def = ctx.definitions.find(name);
        if (def == ctx.definitions.end() || !ctx.resolving.insert(name).second) return std::nullopt;   // n = m, m = n
        const auto value = evaluateConstantExpression(def->second.first,
                                                      [&](const std::string& n) { return knownValue(ctx, n); });
        ctx.resolving.erase(name);
        ctx.values[name] = value;
        return value;
    }

    /// Why `name` (a name without value in a bound of the main program) is not known before the solve.
    static std::string explainUnknown(const CompileTimeContext& ctx, const std::string& name) {
        if (name.find('[') != std::string::npos)
            return "'" + name + "' is an array element, known only at the solve";
        const auto def = ctx.definitions.find(name);
        if (def == ctx.definitions.end())
            return "'" + name + "' is not set by an equation '" + name + " = <constant expression>' of the main program";
        return "'" + name + "' is set by '" + name + " = " + formatExpr(def->second.first) + "' (line " +
               std::to_string(def->second.second) + "), which is not a constant expression";
    }

    /// Resolve in place the `which` ("first"/"last") bound of `owner` ("DUPLICATE i", "SUM index i"), a loop
    /// of the main program or an indexed aggregation: its integer value, or (a bound using an enclosing
    /// index) the expression with the known names replaced by their values. False, reported, when it is
    /// neither. In a FUNCTION or PROCEDURE only constants are known: the other names are values of the run.
    bool resolveDuplicateBound(ExprPtr& bound, const char* which, const std::string& owner, CompileTimeContext& ctx) {
        std::string unknown;
        bool usesIndex = false;
        const auto value = evaluateConstantExpression(bound, [&](const std::string& n) -> std::optional<double> {
            if (ctx.iterators.count(n) && ctx.routine.empty()) {
                usesIndex = true;
                return 1.0;   // the form is checked here, the value in each copy (IR)
            }
            return knownValue(ctx, n);
        }, &unknown);
        const std::string what = owner + ": the " + which + " bound '" + formatExpr(bound) + "'";
        if (!value) {
            if (unknown.empty())
                ctx.error(what + " is not an arithmetic expression" + kBoundRule);
            else if (!ctx.routine.empty())
                ctx.error(what + " is known only when " + ctx.routine + " runs: in a FUNCTION or PROCEDURE the "
                          "bounds of an indexed SUM or PRODUCT are numbers or $CONSTANTs");
            else
                ctx.error(what + " is not known before the solve: " + explainUnknown(ctx, unknown) + kBoundRule);
            return false;
        }
        if (usesIndex) {
            replaceKnownNames(bound, ctx);
            return true;
        }
        const auto integer = asIntegerBound(*value);
        if (!integer) {
            ctx.error(what + " is " + formatNumber(*value) + ", not an integer");
            return false;
        }
        bound = makeNumber(*integer, bound->sourceLineNumber);
        return true;
    }

    void replaceKnownNames(ExprPtr& expr, CompileTimeContext& ctx) {
        if (expr->is<Variable>()) {
            const auto& var = expr->as<Variable>();
            if (!var.indices.empty() || ctx.iterators.count(var.name)) return;
            if (const auto value = knownValue(ctx, var.name)) expr->node = NumberLiteral{*value};
        } else if (expr->is<UnaryOp>()) {
            replaceKnownNames(expr->as<UnaryOp>().operand, ctx);
        } else if (expr->is<BinaryOp>()) {
            replaceKnownNames(expr->as<BinaryOp>().left, ctx);
            replaceKnownNames(expr->as<BinaryOp>().right, ctx);
        }
    }

    /// "array range 'x[1..n]'", with the routine when it is in one.
    static std::string rangeLabel(const std::string& text, const CompileTimeContext& ctx, bool inHeader) {
        std::string label = "array range '" + text + "'";
        if (!ctx.routine.empty()) label += (inHeader ? " in the header of " : " in ") + ctx.routine;
        return label;
    }

    /// The integer value of the `which` bound of the range `text`; nullopt (reported) when it is not known.
    std::optional<int> rangeBound(const ExprPtr& bound, const char* which, const std::string& text,
                                  CompileTimeContext& ctx, bool inHeader) {
        std::string unknown, index;
        const auto value = evaluateConstantExpression(bound, [&](const std::string& n) -> std::optional<double> {
            if (ctx.iterators.count(n) && ctx.routine.empty()) {
                if (index.empty()) index = n;
                return std::nullopt;
            }
            return knownValue(ctx, n);
        }, &unknown);
        const std::string label = rangeLabel(text, ctx, inHeader);
        if (!value) {
            if (!ctx.routine.empty() && !unknown.empty())
                ctx.error(label + ": the bounds of a range in a FUNCTION or PROCEDURE are numbers or $CONSTANTs ('" +
                          unknown + "' is known only when the routine runs)");
            else if (!index.empty())
                ctx.error(label + ": its bounds cannot depend on the index '" + index + "' of the enclosing DUPLICATE "
                          "(not supported: write the elements one by one)");
            else
                ctx.error(label + ": the " + which + " index " +
                          (unknown.empty() ? std::string("is not an arithmetic expression")
                                           : "is not known before the solve: " + explainUnknown(ctx, unknown)) + kBoundRule);
            return std::nullopt;
        }
        const auto integer = asIntegerBound(*value);
        if (!integer) ctx.error(label + ": the " + which + " index is " + formatNumber(*value) + ", not an integer");
        return integer;
    }

    // --- ranges and constants ---------------------------------------------------

    /// The elements of `var` (an index at least is a range a..b), the last index varying fastest; empty
    /// when a bound is not known (reported).
    std::vector<Variable> rangeElements(Variable& var, CompileTimeContext& ctx, bool inHeader) {
        const std::string text = formatExpr(makeArrayVariable(var.name, var.indices));
        std::vector<Variable> elements{Variable{var.name, {}}};
        for (auto& index : var.indices) {
            std::vector<ExprPtr> values;
            if (!isRangeIndex(index)) {
                walkExpression(index, ctx);
                values.push_back(index);
            } else {
                const auto& range = index->as<BinaryOp>();
                const auto first = rangeBound(range.left, "first", text, ctx, inHeader);
                const auto last = first ? rangeBound(range.right, "last", text, ctx, inHeader) : std::nullopt;
                if (!last) return {};
                if (*first > *last) {
                    ctx.error(rangeLabel(text, ctx, inHeader) + ": the first index (" + std::to_string(*first) +
                              ") is greater than the last (" + std::to_string(*last) + ")");
                    return {};
                }
                if (static_cast<size_t>(*last - *first + 1) * elements.size() > kMaxRangeElements) {
                    ctx.error(rangeLabel(text, ctx, inHeader) + ": more than " + std::to_string(kMaxRangeElements) +
                              " elements");
                    return {};
                }
                for (int k = *first; k <= *last; ++k) values.push_back(makeNumber(k, index->sourceLineNumber));
            }
            std::vector<Variable> next;
            for (const auto& element : elements)
                for (const auto& value : values) {
                    next.push_back(element);
                    next.back().indices.push_back(value);
                }
            elements = std::move(next);
        }
        return elements;
    }

    /// Replace the ranges of an argument list by their elements; the other arguments are walked. True when
    /// the list held a range.
    bool expandArguments(std::vector<ExprPtr>& args, CompileTimeContext& ctx) {
        bool hadRange = false;
        std::vector<ExprPtr> expanded;
        for (auto& arg : args) {
            if (!arg->is<Variable>() || !hasRange(arg->as<Variable>())) {
                walkExpression(arg, ctx);
                expanded.push_back(arg);
                continue;
            }
            hadRange = true;
            for (auto& element : rangeElements(arg->as<Variable>(), ctx, false))
                expanded.push_back(makeArrayVariable(element.name, std::move(element.indices), arg->sourceLineNumber));
        }
        if (hadRange) args = std::move(expanded);
        return hadRange;
    }

    /// Replace the ranges of a routine header (`SP[1..5]`) by the names of their elements as the body
    /// reads them (`SP[1]`, ..., `SP[5]`). The bounds are numbers or constants.
    void expandHeader(std::vector<std::string>& names, bool& hasRangeOut, CompileTimeContext& ctx) {
        std::vector<std::string> expanded;
        for (const auto& name : names) {
            if (name.find("..") == std::string::npos) {
                expanded.push_back(name);
                continue;
            }
            clearExpressionError();
            ExprPtr parsed = parseExpression(name, ctx.line);
            clearExpressionError();
            if (!parsed || !parsed->is<Variable>() || !hasRange(parsed->as<Variable>())) {
                ctx.error("Could not parse the argument '" + name + "' in the header of " + ctx.routine);
                continue;
            }
            hasRangeOut = true;
            for (const auto& element : rangeElements(parsed->as<Variable>(), ctx, true))
                expanded.push_back(element.flattenedName());
        }
        names = std::move(expanded);
    }

    /// A call that uses a range, or calls a routine whose header has one, gives as many values as declared.
    static void checkCallSize(const std::string& call, const std::string& routine, size_t given, bool callHasRange,
                              bool outputs, CompileTimeContext& ctx) {
        const auto it = ctx.routines.find(routine);
        if (it == ctx.routines.end() || (!callHasRange && !it->second.hasRange)) return;
        const auto& sig = it->second;
        const size_t declared = outputs ? sig.outputs : sig.inputs;
        if (given == declared) return;
        const std::string what = call.rfind("CALL ", 0) == 0 ? (outputs ? " outputs" : " inputs") : " arguments";
        ctx.error(call + " gives " + std::to_string(given) + what + " (an array range counts its elements), but " +
                  sig.kind + " " + sig.name + " (line " + std::to_string(sig.line) + ") " +
                  (outputs ? "returns " : "takes ") + std::to_string(declared));
    }

    /// A variable assigned by a statement (left-hand side, CALL output) must not be a $CONSTANT.
    void checkNotConstant(const Variable& var, CompileTimeContext& ctx) {
        if (!var.indices.empty() || ctx.locals.count(var.name)) return;
        const auto it = userConstants_.find(var.name);
        if (it != userConstants_.end())
            ctx.error("'" + var.name + "' is a constant ($CONSTANT, line " + std::to_string(it->second.line) +
                      "): its value cannot be changed");
    }

    /// SUM(expr, i=a, b) and PRODUCT(expr, i=a, b): replaced in place by the variadic call over the copies of
    /// `expr` with i = a, ..., b (expandIndexedAggregate, the substitution of the DUPLICATE copies), then
    /// walked. The bounds are resolved like those of a DUPLICATE; one that uses the index of an enclosing
    /// DUPLICATE leaves the call as it is, the IR expands it in each copy of that loop.
    void resolveIndexedAggregate(ExprPtr& expr, CompileTimeContext& ctx) {
        auto& call = expr->as<FunctionCall>();
        std::string owner = call.name;
        std::transform(owner.begin(), owner.end(), owner.begin(), ::toupper);
        owner += " index " + call.namedArgs[0].first;
        if (!resolveDuplicateBound(call.namedArgs[0].second, "first", owner, ctx) ||
            !resolveDuplicateBound(call.args[1], "last", owner, ctx))
            return;
        try {
            if (const auto expanded = expandIndexedAggregate(expr)) {
                *expr = *expanded;
                return walkExpression(expr, ctx);
            }
        } catch (const ModelError& e) {
            return ctx.error(e.what());
        }
        const std::string& index = call.namedArgs[0].first;
        const bool added = ctx.iterators.insert(index).second;
        walkExpression(call.args[0], ctx);
        if (added) ctx.iterators.erase(index);
    }

    /// The main program uses `name` as a variable (and it is no loop index): it is a real index when a unit
    /// is its name (CS-BUG-VAR-ANNOTATION). Only recorded when such an index was parsed.
    void noteUsed(const std::string& name, CompileTimeContext& ctx) const {
        if (unitIndexSeen_ && ctx.routine.empty() && !ctx.iterators.count(name)) ctx.used.insert(name);
    }

    void walkExpression(ExprPtr& expr, CompileTimeContext& ctx) {
        if (!expr) return;
        if (expr->is<Variable>()) {
            auto& var = expr->as<Variable>();
            if (var.indices.empty()) {   // a $CONSTANT becomes its value
                noteUsed(var.name, ctx);
                const auto it = ctx.locals.count(var.name) ? userConstants_.end() : userConstants_.find(var.name);
                if (it != userConstants_.end()) expr->node = NumberLiteral{it->second.value};
                return;
            }
            if (const std::string* unit = unitIndexOf(var)) {   // `a [kW]`, unless `kW` is a loop index
                if (ctx.routine.empty() && !ctx.iterators.count(*unit)) ctx.unitSubscripts.push_back({ctx.line, var.name, *unit});
                return;
            }
            if (hasRange(var)) {
                ctx.error("array range '" + formatExpr(expr) + "' is only allowed as an argument of a function or of "
                          "a CALL, or in the header of a FUNCTION or PROCEDURE (EES array range notation)");
                return;
            }
            for (auto& index : var.indices) walkExpression(index, ctx);
        } else if (expr->is<UnaryOp>()) {
            walkExpression(expr->as<UnaryOp>().operand, ctx);
        } else if (expr->is<BinaryOp>()) {
            walkExpression(expr->as<BinaryOp>().left, ctx);
            walkExpression(expr->as<BinaryOp>().right, ctx);
        } else if (expr->is<FunctionCall>()) {
            auto& call = expr->as<FunctionCall>();
            if (isIndexedAggregate(call)) return resolveIndexedAggregate(expr, ctx);
            const bool hadRange = expandArguments(call.args, ctx);
            for (auto& named : call.namedArgs) walkExpression(named.second, ctx);
            checkCallSize(call.name + "(...)", call.name, call.args.size(), hadRange, false, ctx);
        }
    }

    /// The left-hand side of an equation or an assignment.
    void walkAssigned(ExprPtr& lhs, CompileTimeContext& ctx) {
        if (!lhs->is<Variable>()) return walkExpression(lhs, ctx);   // an implicit equation of the main program
        checkNotConstant(lhs->as<Variable>(), ctx);
        noteUsed(lhs->as<Variable>().name, ctx);
        if (!lhs->as<Variable>().indices.empty()) walkExpression(lhs, ctx);
    }

    void walkCall(ProcedureCall& call, CompileTimeContext& ctx) {
        const bool inputRange = expandArguments(call.inputArgs, ctx);
        bool outputRange = false;
        std::vector<Variable> outputs;
        std::vector<ExprPtr> outputExprs;   // kept aligned with the outputs when the CALL has expression outputs
        for (size_t i = 0; i < call.outputVars.size(); ++i) {
            ExprPtr expression = call.outputExprs.empty() ? nullptr : call.outputExprs[i];
            Variable& var = call.outputVars[i];
            if (expression || !hasRange(var)) {
                if (expression) walkExpression(expression, ctx);
                checkNotConstant(var, ctx);
                noteUsed(var.name, ctx);
                for (auto& index : var.indices) walkExpression(index, ctx);
                outputs.push_back(var);
                outputExprs.push_back(expression);
                continue;
            }
            outputRange = true;
            for (auto& element : rangeElements(var, ctx, false)) {
                outputs.push_back(std::move(element));
                outputExprs.push_back(nullptr);
            }
        }
        if (outputRange) {
            call.outputVars = std::move(outputs);
            if (!call.outputExprs.empty()) call.outputExprs = std::move(outputExprs);
        }
        const auto sig = ctx.routines.find(call.name);
        if (sig == ctx.routines.end() || sig->second.kind != "PROCEDURE") return;
        checkCallSize("CALL " + call.name, call.name, call.inputArgs.size(), inputRange || outputRange, false, ctx);
        checkCallSize("CALL " + call.name, call.name, call.outputVars.size(), inputRange || outputRange, true, ctx);
    }

    void walkStatements(std::vector<StmtPtr>& statements, CompileTimeContext& ctx) {
        for (auto& stmt : statements) {
            ctx.line = stmt->sourceLineNumber;
            if (stmt->is<Equation>()) {
                walkAssigned(stmt->as<Equation>().lhs, ctx);
                walkExpression(stmt->as<Equation>().rhs, ctx);
            } else if (stmt->is<Assignment>()) {
                walkAssigned(stmt->as<Assignment>().lhs, ctx);
                walkExpression(stmt->as<Assignment>().rhs, ctx);
            } else if (stmt->is<ProcedureCall>()) {
                walkCall(stmt->as<ProcedureCall>(), ctx);
            } else if (stmt->is<Duplicate>()) {
                auto& dup = stmt->as<Duplicate>();
                if (ctx.routine.empty()) {   // a loop of a body is run by the evaluator, with its bounds
                    const std::string owner = "DUPLICATE " + dup.iteratorVar;
                    if (resolveDuplicateBound(dup.start, "first", owner, ctx))
                        resolveDuplicateBound(dup.end, "last", owner, ctx);
                } else {
                    walkExpression(dup.start, ctx);
                    walkExpression(dup.end, ctx);
                }
                const bool added = ctx.iterators.insert(dup.iteratorVar).second;
                walkStatements(dup.body, ctx);
                if (added) ctx.iterators.erase(dup.iteratorVar);
            } else if (stmt->is<IfThenElse>()) {
                auto& ite = stmt->as<IfThenElse>();
                walkExpression(ite.condition, ctx);
                walkStatements(ite.thenBranch, ctx);
                walkStatements(ite.elseBranch, ctx);
            } else if (stmt->is<RepeatUntil>()) {
                auto& rep = stmt->as<RepeatUntil>();
                walkStatements(rep.body, ctx);
                if (rep.condition && rep.condition->sourceLineNumber > 0) ctx.line = rep.condition->sourceLineNumber;
                walkExpression(rep.condition, ctx);
            } else if (stmt->is<FunctionDefinition>()) {
                auto& def = stmt->as<FunctionDefinition>();
                ctx.routine = "FUNCTION '" + def.name + "'";
                ctx.locals = {def.name};
                for (const auto& p : def.parameters) ctx.locals.insert(arrayBaseName(p));
                walkStatements(def.body, ctx);
                ctx.routine.clear();
                ctx.locals.clear();
            } else if (stmt->is<ProcedureDefinition>()) {
                auto& def = stmt->as<ProcedureDefinition>();
                ctx.routine = "PROCEDURE '" + def.name + "'";
                for (const auto& p : def.inputs) ctx.locals.insert(arrayBaseName(p));
                for (const auto& p : def.outputs) ctx.locals.insert(arrayBaseName(p));
                walkStatements(def.body, ctx);
                ctx.routine.clear();
                ctx.locals.clear();
            }
        }
    }

    // ------------------------------------------------------------------------
    // IF ... THEN ... ELSE ... ENDIF (EES: only inside FUNCTION/PROCEDURE bodies)
    //
    //   IF (cond) THEN stmts [ELSE stmts]        single line, no ENDIF
    //   IF (cond) THEN stmts ELSE                single line ending with ELSE: the ELSE part
    //       stmts                                is on the following lines, up to the ENDIF
    //   ENDIF
    //   IF (cond) THEN                           block form, closed by ENDIF
    //       stmts
    //   [ELSE [stmt]
    //       stmts]
    //   ENDIF
    //
    // ENDIFs may share a line (`ENDIF;ENDIF`: one per open IF) and be followed by statements after a `;`.
    // An ENDIF that closes no IF is a warning (P008) in a body, as EES accepts it; an error elsewhere.
    // The condition may use = <> < > <= >= (numbers or strings) combined with
    // AND / OR and parentheses; the parentheses around the whole condition are
    // optional. Branch statements may themselves be IF statements.
    // ------------------------------------------------------------------------
    static constexpr const char* kIfOutsideProcedureMessage =
        "IF ... THEN statements are only allowed inside FUNCTION and PROCEDURE bodies "
        "(in the main program use the IF(cond, a, b) or IF(A, B, X, Y, Z) function)";

    // A character that can be part of a name (see firstInvalidNameChar): a keyword is a keyword only
    // when none follows it. Bytes >= 0x80 are the UTF-8 letters and symbols of names.
    static bool isIdentChar(char c) {
        const unsigned char u = static_cast<unsigned char>(c);
        return std::isalnum(u) || u >= 0x80 || c == '_' || c == '$' || c == '#' || c == '%' || c == '|' || c == '\\';
    }

    // Position (>= from) of the whole-word, case-insensitive keyword `kw`
    // (given in upper case) in `s`, at parenthesis depth 0 and outside 'strings',
    // "comments" and {comments}; npos when absent.
    static size_t findKeywordAtDepth0(const std::string& s, const std::string& kw, size_t from = 0) {
        bool inSingle = false, inDouble = false;
        int braceDepth = 0, depth = 0;
        for (size_t i = 0; i < s.size(); ++i) {
            const char c = s[i];
            if (inSingle) { if (c == '\'') inSingle = false; continue; }
            if (inDouble) { if (c == '"') inDouble = false; continue; }
            if (braceDepth > 0) {
                if (c == '{') braceDepth++;
                else if (c == '}') braceDepth--;
                continue;
            }
            if (c == '\'') { inSingle = true; continue; }
            if (c == '"') { inDouble = true; continue; }
            if (c == '{') { braceDepth++; continue; }
            if (c == '(' || c == '[') { depth++; continue; }
            if (c == ')' || c == ']') { depth--; continue; }
            if (depth != 0 || i < from) continue;
            if (i > 0 && isIdentChar(s[i - 1])) continue;
            if (i + kw.size() > s.size()) continue;
            bool match = true;
            for (size_t k = 0; k < kw.size() && match; ++k)
                match = std::toupper(static_cast<unsigned char>(s[i + k])) == kw[k];
            if (match && (i + kw.size() == s.size() || !isIdentChar(s[i + kw.size()])))
                return i;
        }
        return std::string::npos;
    }

    // True when the '(' at s[0] is closed by the last character of `s`.
    static bool isFullyParenthesized(const std::string& s) {
        if (s.size() < 2 || s.front() != '(' || s.back() != ')') return false;
        bool inSingle = false;
        int depth = 0;
        for (size_t i = 0; i < s.size(); ++i) {
            if (s[i] == '\'') { inSingle = !inSingle; continue; }
            if (inSingle) continue;
            if (s[i] == '(') depth++;
            else if (s[i] == ')') {
                depth--;
                if (depth == 0 && i + 1 < s.size()) return false;
            }
        }
        return depth == 0;
    }

    // Parse the condition of an IF statement. Relations produce BinaryOp nodes
    // with op "=", "<>", "<", ">", "<=" or ">="; "and"/"or" combine them
    // (OR binds weaker than AND). Evaluates to +1 (true) / -1 (false).
    ExprPtr parseCondition(const std::string& text, int lineNum) {
        const std::string s = trim(text);
        if (s.empty()) return nullptr;

        for (const char* kw : {"OR", "AND"}) {
            const size_t pos = findKeywordAtDepth0(s, kw);
            if (pos == std::string::npos) continue;
            auto left = parseCondition(s.substr(0, pos), lineNum);
            auto right = parseCondition(s.substr(pos + std::string(kw).size()), lineNum);
            if (!left || !right) return nullptr;
            return makeBinaryOp(std::string(kw) == "OR" ? "or" : "and", left, right, lineNum);
        }

        if (isFullyParenthesized(s))
            return parseCondition(s.substr(1, s.size() - 2), lineNum);

        // Leftmost relational operator outside parentheses and strings.
        bool inSingle = false;
        int depth = 0;
        for (size_t i = 0; i < s.size(); ++i) {
            const char c = s[i];
            if (c == '\'') { inSingle = !inSingle; continue; }
            if (inSingle) continue;
            if (c == '(' || c == '[') { depth++; continue; }
            if (c == ')' || c == ']') { depth--; continue; }
            if (depth != 0 || (c != '<' && c != '>' && c != '=')) continue;
            std::string op(1, c);
            if (i + 1 < s.size() && ((c == '<' && (s[i + 1] == '=' || s[i + 1] == '>')) ||
                                     (c == '>' && s[i + 1] == '=')))
                op += s[i + 1];
            auto left = parseExpression(s.substr(0, i), lineNum);
            auto right = parseExpression(s.substr(i + op.size()), lineNum);
            if (!left || !right) return nullptr;
            return makeBinaryOp(op, left, right, lineNum);
        }
        return parseExpression(s, lineNum);  // no relation: numeric truth value
    }

    // Parse one statement of a single-line IF branch: a nested IF statement or
    // simple statements. Returns the terminator of a prematurely ended nested
    // block (see parseIfStatement), "" otherwise.
    std::string parseBranchText(const std::string& text, std::istream& stream, int& lineNumber,
                                std::vector<StmtPtr>& branch, ParseResult& result,
                                const std::string& srcLine) {
        const std::string t = trim(text);
        if (t.empty()) return "";
        if (findKeywordAtDepth0(t, "IF") == 0)
            return parseIfStatement(t, stream, lineNumber, branch, result, srcLine);
        parseSimpleStatements(t, lineNumber, branch, result, srcLine);
        return "";
    }

    // Parse an IF statement whose first line is `ifText` (comments removed,
    // trimmed, starting with IF) and append it to `body`. For the block form, and for
    // a single-line form ending with ELSE, the following lines are read from `stream`
    // up to the matching ENDIF.
    // Returns "" when the statement is complete, or the terminator ("END",
    // "UNTIL") that ended the enclosing block before the ENDIF, which the caller
    // must propagate.
    std::string parseIfStatement(const std::string& ifText, std::istream& stream, int& lineNumber,
                                 std::vector<StmtPtr>& body, ParseResult& result,
                                 const std::string& srcLine) {
        const int ifLine = lineNumber;
        const size_t thenPos = findKeywordAtDepth0(ifText, "THEN", 2);
        if (thenPos == std::string::npos) {
            result.errors.push_back({ifLine, 0, "IF statement without THEN", srcLine});
            return "";
        }
        const std::string condText = trim(ifText.substr(2, thenPos - 2));
        clearExpressionError();
        ExprPtr cond = parseCondition(condText, ifLine);
        if (!cond)
            result.errors.push_back({ifLine, 0, "Could not parse IF condition '" + condText + "'"
                                     + (expressionErrorText().empty() ? "" : ": " + expressionErrorText()), srcLine});

        std::vector<StmtPtr> thenBranch, elseBranch;
        std::string premature;
        const std::string rest = trim(ifText.substr(thenPos + 4));
        std::vector<StmtPtr>* current = &thenBranch;
        bool sawElse = false;
        bool readBlock = rest.empty();
        if (!rest.empty()) {
            // Single-line form: IF (cond) THEN stmts [ELSE stmts]
            const size_t elsePos = findKeywordAtDepth0(rest, "ELSE");
            premature = parseBranchText(elsePos == std::string::npos ? rest : rest.substr(0, elsePos),
                                        stream, lineNumber, thenBranch, result, srcLine);
            if (premature.empty() && elsePos != std::string::npos) {
                const std::string elseText = rest.substr(elsePos + 4);
                premature = parseBranchText(elseText, stream, lineNumber, elseBranch, result, srcLine);
                if (premature.empty() && trim(elseText).empty()) {
                    // The line ends with ELSE: the ELSE part is on the following lines, up to the ENDIF
                    current = &elseBranch;
                    sawElse = true;
                    readBlock = true;
                }
            }
        }
        if (readBlock) {
            // Block form: lines up to ENDIF, with an optional ELSE
            for (;;) {
                std::string elseLine;
                const std::string term = parseBodyLines(stream, lineNumber, *current, result,
                                                        nullptr, /*inIfBlock=*/true, &elseLine);
                if (term == "ENDIF") break;
                if (term == "ELSE") {
                    if (sawElse)
                        result.errors.push_back({lineNumber, 0, "IF statement with two ELSE parts", elseLine});
                    sawElse = true;
                    current = &elseBranch;
                    // Statements after ELSE on the same line, e.g. `ELSE x = 1` or `ELSE IF (c) THEN`
                    premature = parseBranchText(removeInlineComments(elseLine).substr(4), stream,
                                                lineNumber, *current, result, elseLine);
                    if (!premature.empty()) break;
                    continue;
                }
                result.errors.push_back({ifLine, 0, "IF ... THEN without a matching ENDIF", srcLine});
                premature = term;  // END / UNTIL end the enclosing block; EOF ("") ends the program
                break;
            }
        }
        if (cond) body.push_back(makeIfThenElse(cond, std::move(thenBranch), std::move(elseBranch), ifLine));
        return premature;
    }

    StmtPtr tryParseComment(const std::string& line, int lineNum) {
        std::string trimmed = trim(line);
        
        // Check for quote-delimited comment ("..." display comments)
        if (trimmed.size() >= 2 && trimmed.front() == '"' && trimmed.back() == '"') {
            // Find the position of the second quote (end of first comment)
            size_t firstQuoteEnd = trimmed.find('"', 1);
            if (firstQuoteEnd == trimmed.size() - 1) {
                // Only one pair of quotes - this is a pure comment
                std::string content = trimmed.substr(1, trimmed.size() - 2);
                return makeComment(content, true, /*isDisplay=*/true, lineNum);
            }
            // Multiple quote pairs - check if there's code between them
            // If there's an '=' outside of quotes, it's likely an equation with inline comments
            std::string stripped = removeInlineComments(trimmed);
            stripped = trim(stripped);
            if (!stripped.empty() && stripped.find('=') != std::string::npos) {
                // This is an equation with inline comments, not a pure comment
                return nullptr;
            }
            // No equation found, treat as comment
            std::string content = trimmed.substr(1, trimmed.size() - 2);
            return makeComment(content, true, /*isDisplay=*/true, lineNum);
        }
        
        // Check for brace-delimited comment (only if it starts with { and is a whole line)
        if (trimmed.size() >= 2 && trimmed.front() == '{' && trimmed.back() == '}') {
            // Make sure it's just a comment, not embedded in an equation
            if (trimmed.find('=') == std::string::npos) {
                std::string content = trimmed.substr(1, trimmed.size() - 2);
                return makeComment(content, true, /*isDisplay=*/false, lineNum);
            }
        }
        
        // Check for C-style comment
        if (trimmed.size() >= 2 && trimmed[0] == '/' && trimmed[1] == '/') {
            std::string content = trimmed.substr(2);
            return makeComment(content, false, /*isDisplay=*/false, lineNum);
        }
        
        return nullptr;
    }

    // Parse the body of a `$IntegralTable` directive into an IntegralTableSpec.
    // EES syntax:  $IntegralTable t:0.1  y  dydt  X[1..5]
    //   - first name: integration variable, optionally followed by ":interval" (a number, or
    //     the name of a variable that holds a constant, read at the initial point)
    //   - remaining names: tabulated variables; "X[lo..hi]" ranges are expanded
    //     to X[1], X[2], ... at parse time so the rest of the pipeline sees a
    //     flat column list (additive — only active inside this directive parser).
    //   - names are separated by spaces or by the EES list delimiters `,` and `;`
    //     (`t:0.5; y; dydt`, `t:0.5, y, dydt`); a delimiter inside brackets belongs to the
    //     array subscript (`X[1,2]`).
    static IntegralTableSpec parseIntegralTableContent(const std::string& content) {
        IntegralTableSpec spec;
        std::vector<std::string> tokens;
        {
            std::string cur;
            int depth = 0;
            for (char c : content) {
                if (c == '[') ++depth;
                else if (c == ']' && depth > 0) --depth;
                const bool space = std::isspace(static_cast<unsigned char>(c)) != 0;
                if (depth == 0 && (space || c == ',' || c == ';')) {
                    if (!cur.empty()) tokens.push_back(cur);
                    cur.clear();
                } else if (!(depth > 0 && space)) {   // blanks inside a subscript are dropped
                    cur += c;
                }
            }
            if (!cur.empty()) tokens.push_back(cur);
        }
        if (tokens.empty()) {
            spec.errorMessage = "no integration variable specified";
            return spec;
        }
        // First token: var[:interval].
        const std::string& first = tokens[0];
        size_t colon = first.find(':');
        if (colon != std::string::npos) {
            spec.integrationVar = first.substr(0, colon);
            const std::string interval = first.substr(colon + 1);
            size_t used = 0;
            double value = 0.0;
            bool isNumber = false;
            try {
                value = std::stod(interval, &used);
                isNumber = used == interval.size();
            } catch (...) {
                isNumber = false;
            }
            if (isNumber) {
                spec.outputInterval = value;
            } else if (!interval.empty() && nameStartsAt(interval, 0)) {
                spec.outputIntervalVar = interval;   // a variable set to a constant (EES)
            } else {
                spec.errorMessage = "invalid output interval '" + interval + "'";
                return spec;
            }
        } else {
            spec.integrationVar = first;
        }
        spec.columns.push_back(spec.integrationVar);

        // Remaining tokens: expand X[lo..hi] ranges, otherwise push as-is.
        for (size_t i = 1; i < tokens.size(); ++i) {
            const std::string& t = tokens[i];
            size_t lb = t.find('[');
            size_t rb = t.find(']');
            size_t dots = t.find("..");
            if (lb != std::string::npos && rb != std::string::npos && dots != std::string::npos
                && lb < dots && dots < rb) {
                std::string base = t.substr(0, lb);
                std::string loStr = t.substr(lb + 1, dots - (lb + 1));
                std::string hiStr = t.substr(dots + 2, rb - (dots + 2));
                try {
                    int lo = std::stoi(loStr);
                    int hi = std::stoi(hiStr);
                    if (lo > hi) std::swap(lo, hi);
                    for (int k = lo; k <= hi; ++k)
                        spec.columns.push_back(base + "[" + std::to_string(k) + "]");
                    continue;
                } catch (...) {
                    // fall through and record the raw token
                }
            }
            spec.columns.push_back(t);
        }
        spec.valid = true;
        return spec;
    }

    static bool isCommonDirective(const Directive& directive) {
        std::string name = directive.name;
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        return name == "common";
    }

    /// The error for a `$COMMON` line in a FUNCTION/PROCEDURE body (CS-BUG-COMMON-PROC): it is not
    /// supported, the shared variables are passed as input arguments (the body is unchanged).
    /// `content` is the text after the directive name: variables separated by commas or semicolons.
    std::string commonDirectiveMessage(const std::string& content) {
        std::string list = removeInlineComments(content);
        std::replace(list.begin(), list.end(), ';', ',');
        std::string names;
        size_t count = 0;
        for (const auto& item : splitArguments(list)) {
            const std::string name = trim(item);
            if (name.empty()) continue;
            names += (count++ ? ", " : "") + name;
        }
        const std::string what = count == 0 ? "the shared variables as input arguments"
                               : names + (count == 1 ? " as an input argument" : " as input arguments");
        return "$COMMON is not supported: pass " + what + " of " + routineKind_ + " " + routineName_ +
               (routineKind_ == "FUNCTION" ? " (and in each call)" : " (and in each CALL)");
    }

    StmtPtr tryParseDirective(const std::string& line, int lineNum, DiagnosticCollector* diag = nullptr) {
        std::string trimmed = trim(line);
        
        if (trimmed.empty() || trimmed[0] != '$') return nullptr;
        
        // Extract directive name
        size_t nameEnd = 1;
        while (nameEnd < trimmed.size() && 
               (std::isalnum(static_cast<unsigned char>(trimmed[nameEnd])) || trimmed[nameEnd] == '_')) {
            nameEnd++;
        }
        
        std::string name = trimmed.substr(1, nameEnd - 1);
        std::string content = trimmed.substr(nameEnd);
        content = trim(content);
        
        // Validate directive name against known EES directives
        if (diag) {
            std::string lower = name;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            static const std::unordered_set<std::string> knownDirectives = {
                "if", "ifnot", "else", "endif",
                "include", "export", "import",
                "tabstops", "unitsystem", "opentable",
                "arrays", "common", "integraltable",
                "integralautostep", "integralstop",  // recognised-but-not-interpreted
                "diagrams", "bookmark", "hiddenvariables",
                "checkunits", "warnings", "complex",
                "reference", "savelookup", "constant",
            };
            if (knownDirectives.find(lower) == knownDirectives.end()) {
                diag->push(DiagnosticSeverity::Warning, "P002",
                           "Unknown directive '$" + name + "'",
                           "parser", lineNum, 1);
            }
            // $IntegralAutoStep / $IntegralStop are recognised for compatibility
            // but not interpreted — integration parameters live in coolsolve.conf
            // (see docs/integral_table.md §5 Configuration).
            if (lower == "integralautostep" || lower == "integralstop") {
                diag->push(DiagnosticSeverity::Warning, "P005",
                           "'$" + name + "' is recognised but not interpreted; "
                           "configure integration via the 'integral*' keys in coolsolve.conf.",
                           "parser", lineNum, 1);
            }
        }

        StmtPtr stmt = makeDirective(name, content, lineNum);
        // Capture the structured payload for $IntegralTable (variable list +
        // optional output interval, with X[1..5] range expansion).
        if (diag) {
            std::string lower = name;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            if (lower == "integraltable") {
                auto& d = stmt->as<Directive>();
                d.integralTableSpec = parseIntegralTableContent(content);
                d.hasIntegralTableSpec = true;
                if (!d.integralTableSpec.valid) {
                    diag->push(DiagnosticSeverity::Warning, "P006",
                               "$IntegralTable: " + d.integralTableSpec.errorMessage,
                               "parser", lineNum, 1);
                }
            }
        }

        return stmt;
    }
    
    bool tryParseFunctionHeader(const std::string& line, std::string& name, std::vector<std::string>& params) {
        std::string cleaned = removeInlineComments(line);
        cleaned = trim(cleaned);
        std::string upper = cleaned;
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        
        if (upper.find("FUNCTION ") != 0) return false;
        
        std::string content = trim(cleaned.substr(9));
        size_t parenStart = content.find('(');
        if (parenStart == std::string::npos) return false;
        
        name = trim(content.substr(0, parenStart));
        
        size_t parenEnd = content.rfind(')');
        if (parenEnd == std::string::npos || parenEnd <= parenStart) return false;
        
        std::string argsStr = content.substr(parenStart + 1, parenEnd - parenStart - 1);
        std::vector<std::string> argStrs = splitArguments(argsStr);
        for (const auto& s : argStrs) {
            params.push_back(trim(s));
            if (s.find("..") != std::string::npos) rangeSeen_ = true;   // X[1..5]: CS-GAP-ARRAY-ARG
        }
        return true;
    }

    bool tryParseProcedureHeader(const std::string& line, std::string& name, std::vector<std::string>& inputs, std::vector<std::string>& outputs) {
        std::string cleaned = removeInlineComments(line);
        cleaned = trim(cleaned);
        std::string upper = cleaned;
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        
        if (upper.find("PROCEDURE ") != 0) return false;
        
        std::string content = trim(cleaned.substr(10));
        size_t parenStart = content.find('(');
        if (parenStart == std::string::npos) return false;
        
        name = trim(content.substr(0, parenStart));
        
        size_t parenEnd = content.rfind(')');
        if (parenEnd == std::string::npos || parenEnd <= parenStart) return false;
        
        std::string argsStr = content.substr(parenStart + 1, parenEnd - parenStart - 1);
        
        // Split by colon to separate inputs and outputs
        size_t colonPos = std::string::npos;
        int depth = 0;
        bool inString = false;
        for (size_t i = 0; i < argsStr.size(); ++i) {
            if (argsStr[i] == '\'') inString = !inString;
            if (inString) continue;
            if (argsStr[i] == '(' || argsStr[i] == '[') depth++;
            else if (argsStr[i] == ')' || argsStr[i] == ']') depth--;
            else if (argsStr[i] == ':' && depth == 0) {
                colonPos = i;
                break;
            }
        }
        
        std::string inputsStr, outputsStr;
        if (colonPos != std::string::npos) {
            inputsStr = argsStr.substr(0, colonPos);
            outputsStr = argsStr.substr(colonPos + 1);
        } else {
            inputsStr = argsStr;
        }
        
        std::vector<std::string> inputStrs = splitArguments(inputsStr);
        for (const auto& s : inputStrs) {
            inputs.push_back(trim(s));
        }
        
        std::vector<std::string> outputStrs = splitArguments(outputsStr);
        for (const auto& s : outputStrs) {
            outputs.push_back(trim(s));
        }
        if (argsStr.find("..") != std::string::npos) rangeSeen_ = true;   // X[1..5]: CS-GAP-ARRAY-ARG
        return true;
    }

    // The reason why a CALL of the built-in procedures ERROR or WARNING is refused; "" when the CALL is valid
    // or is that of another procedure. EES help (ERROR and WARNING procedures): they are used only within a
    // FUNCTION or a PROCEDURE, with no outputs, as CALL P(X), CALL P('message XXXF1', X) or CALL P('message');
    // the message is a string, so a numeric first argument is the value and cannot be followed by another one.
    static std::string messageProcedureError(const std::string& name, const std::vector<ExprPtr>& args,
                                             bool hasOutputs, bool inRoutine) {
        std::string upper = name;
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        if (upper != "ERROR" && upper != "WARNING") return "";
        const std::string call = "CALL " + upper;
        if (!inRoutine) return call + " can only be used inside a FUNCTION or PROCEDURE (as in EES)";
        if (hasOutputs) return call + " has no outputs: " + call + "('message XXXA1', value)";
        if (args.empty())
            return call + " needs an argument: " + call + "('message XXXA1', value), " + call + "('message') or " + call + "(value)";
        if (args.size() > 2)
            return call + " takes a message and a value at most: " + call + "('message XXXA1', value)";
        if (args.size() == 2 && !isStringExpression(args[0]))
            return "the first argument of " + call + " must be the message (a string) when a value follows: " +
                   call + "('message XXXA1', value)";
        return "";
    }

    // True for a CALL of the built-in procedure Colebrook(Re, RR : f) (Darcy friction factor, see
    // colebrookFrictionFactor), which the colon-less form CALL Colebrook(Re, RR, f) also calls. A PROCEDURE
    // Colebrook defined by the model before the call replaces it: the call is then that of any procedure.
    bool isBuiltinColebrook(const std::string& name) const {
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        return lower == "colebrook" && !procedureNames_.count(name);
    }

    // `errorMsg`, when given, receives the reason why a CALL statement was refused
    // (an argument that cannot be parsed, an output that is not a variable).
    StmtPtr tryParseProcedureCall(const std::string& line, int lineNum, std::string* errorMsg = nullptr) {
        clearExpressionError();
        std::string cleaned = removeInlineComments(line);
        cleaned = trim(cleaned);
        
        std::string upper = cleaned;
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        
        if (upper.find("CALL ") != 0) return nullptr;
        
        std::string callContent = trim(cleaned.substr(5));
        size_t parenStart = callContent.find('(');
        if (parenStart == std::string::npos) return nullptr;
        
        std::string name = trim(callContent.substr(0, parenStart));
        
        const size_t parenEnd = findClosingBracket(callContent, parenStart);
        if (parenEnd == std::string::npos || callContent[parenEnd] != ')') return nullptr;

        // Nothing may follow the closing parenthesis but a statement separator or a "[unit]"
        // comment, which the comment remover keeps (the text was ignored: CS-BUG-TRAILING-TEXT)
        std::string trailing = trim(callContent.substr(parenEnd + 1));
        for (bool more = true; more && !trailing.empty();) {
            more = false;
            if (trailing.back() == ';') {
                trailing = trim(trailing.substr(0, trailing.size() - 1));
                more = true;
            } else if (trailing.front() == '"') {
                const size_t close = trailing.find('"', 1);
                if (close != std::string::npos) {
                    trailing = trim(trailing.substr(close + 1));
                    more = true;
                }
            }
        }
        if (!trailing.empty()) {
            if (errorMsg) *errorMsg = unexpectedTextMessage(trailing, "the arguments of", name);
            return nullptr;
        }
        
        std::string argsStr = callContent.substr(parenStart + 1, parenEnd - parenStart - 1);
        
        // Split by colon to separate inputs and outputs
        size_t colonPos = std::string::npos;
        int depth = 0;
        bool inString = false;
        for (size_t i = 0; i < argsStr.size(); ++i) {
            if (argsStr[i] == '\'') inString = !inString;
            if (inString) continue;
            if (argsStr[i] == '(' || argsStr[i] == '[') depth++;
            else if (argsStr[i] == ')' || argsStr[i] == ']') depth--;
            else if (argsStr[i] == ':' && depth == 0) {
                colonPos = i;
                break;
            }
        }
        
        std::string inputsStr, outputsStr;
        if (colonPos != std::string::npos) {
            inputsStr = argsStr.substr(0, colonPos);
            outputsStr = argsStr.substr(colonPos + 1);
        } else {
            inputsStr = argsStr;
        }
        
        // EES writes the built-in Colebrook without the colon, CALL Colebrook(Re, RR, f): its last argument is the output
        const bool colebrook = isBuiltinColebrook(name);
        bool colonLessColebrook = false;
        if (colebrook && colonPos == std::string::npos) {
            const std::vector<std::string> args = splitArguments(inputsStr);
            if (args.size() == 3) {
                inputsStr = args[0] + "," + args[1];
                outputsStr = args[2];
                colonLessColebrook = true;
            }
        }

        // A faulty argument refuses the whole statement (it is never skipped)
        auto refuse = [&](const std::string& fallback) -> StmtPtr {
            if (errorMsg) *errorMsg = expressionErrorText().empty() ? fallback : expressionErrorText();
            return nullptr;
        };

        std::vector<ExprPtr> inputArgs;
        std::vector<std::string> inputStrs = splitArguments(inputsStr);
        for (size_t i = 0; i < inputStrs.size(); ++i) {
            auto expr = parseExpression(inputStrs[i], lineNum);
            if (!expr)
                return refuse("Could not parse argument " + std::to_string(i + 1) + " of '" + name
                              + "': '" + trim(inputStrs[i]) + "'");
            inputArgs.push_back(expr);
        }
        
        // The built-in procedures ERROR and WARNING have a grammar of their own (EES help)
        {
            const std::string why = messageProcedureError(name, inputArgs, colonPos != std::string::npos,
                                                          procedureBodyDepth_ > 0);
            if (!why.empty()) return refuse(why);
        }

        // An output is a variable (possibly an array element). In the main program EES also accepts an
        // expression (`CALL p(x : c/1000)`): the output of the routine is equated to it. Inside a
        // FUNCTION/PROCEDURE the statements run in order and nothing solves such an equation.
        std::vector<Variable> outputVars;
        std::vector<ExprPtr> outputExprs;   // filled only when an output is an expression
        std::vector<std::string> outputStrs = splitArguments(outputsStr);
        for (size_t i = 0; i < outputStrs.size(); ++i) {
            std::string ts = trim(outputStrs[i]);
            if (ts.empty()) continue;
            
            auto expr = parseExpression(ts, lineNum);
            if (!expr)
                return refuse("Could not parse output " + std::to_string(i + 1) + " of '" + name
                              + "': '" + ts + "'");
            if (expr->is<Variable>()) {
                outputVars.push_back(expr->as<Variable>());
                if (!outputExprs.empty()) outputExprs.push_back(nullptr);
                continue;
            }
            if (colonLessColebrook)   // read as an input, so an expression there is more likely a mistake
                return refuse("the third argument of CALL Colebrook(Re, RR, f) is its output f and must be a variable, "
                              "got '" + ts + "'");
            if (procedureBodyDepth_ > 0)
                return refuse("Output " + std::to_string(i + 1) + " of '" + name + "' is the expression '" + ts
                              + "': inside a FUNCTION or PROCEDURE the statements run in order and the expression "
                                "cannot be solved for; write the result into a variable and compute the expression from it");
            if (outputExprs.empty()) outputExprs.assign(outputVars.size(), nullptr);
            outputVars.push_back(Variable{});   // placeholder, replaced when the IR is built
            outputExprs.push_back(expr);
        }

        if (colebrook && (inputArgs.size() != 2 || outputVars.size() != 1))
            return refuse("CALL Colebrook takes 2 inputs (Re, RR) and 1 output (f): CALL Colebrook(Re, RR : f) or "
                          "CALL Colebrook(Re, RR, f); got " + std::to_string(inputArgs.size()) + " input(s) and " +
                          std::to_string(outputVars.size()) + " output(s)");
        
        return makeProcedureCall(name, std::move(inputArgs), std::move(outputVars), lineNum, std::move(outputExprs));
    }

    StmtPtr tryParseEquationOrAssignment(const std::string& line, int lineNum,
                                        std::string* errorMsg = nullptr) {
        // Clear any stale expression error from a previous call.
        clearExpressionError();

        // Remove inline comments first, extracting any "" comment
        std::string inlineComment;
        std::string cleaned = removeInlineComments(line, &inlineComment);
        cleaned = trim(cleaned);
        
        if (cleaned.empty()) return nullptr;
        
        // Strip trailing semicolon (EES uses semicolons as statement terminators/separators)
        if (!cleaned.empty() && cleaned.back() == ';') {
            cleaned.pop_back();
            cleaned = trim(cleaned);
        }
        
        if (cleaned.empty()) return nullptr;
        
        // Find the main '=' or ':=' sign (not inside function calls or arrays)
        std::string op;
        int opPos = findMainOperator(cleaned, op);

        // Reject '==' — EES uses only '=' for equations
        if (op == "==") {
            if (errorMsg) *errorMsg = "Invalid operator '=='. EES uses '=' for equations, not '=='";
            return nullptr;
        }

        if (opPos < 0) return nullptr;
        
        std::string lhsStr = trim(cleaned.substr(0, opPos));
        std::string rhsStr = trim(cleaned.substr(opPos + op.size()));
        
        // Extract units if present
        std::string units;
        size_t unitsStart = rhsStr.find("\"[");
        if (unitsStart != std::string::npos) {
            size_t unitsEnd = rhsStr.find("]\"", unitsStart);
            if (unitsEnd != std::string::npos) {
                units = rhsStr.substr(unitsStart + 2, unitsEnd - unitsStart - 2);
                rhsStr = trim(rhsStr.substr(0, unitsStart));
            }
        }
        
        if (lhsStr.empty() || rhsStr.empty()) return nullptr;

        // Validate expression syntax on both sides
        std::string syntaxErr;
        if (!validateExpressionSyntax(lhsStr, syntaxErr)) {
            if (errorMsg) *errorMsg = syntaxErr;
            return nullptr;
        }
        if (!validateExpressionSyntax(rhsStr, syntaxErr)) {
            if (errorMsg) *errorMsg = syntaxErr;
            return nullptr;
        }
        
        auto lhs = parseExpression(lhsStr, lineNum);
        auto rhs = parseExpression(rhsStr, lineNum);
        
        // If RHS parsing failed, it might be due to unquoted units (e.g., 10 [m])
        // Try stripping units in square brackets at the end. Not when a specific
        // reason was found (an argument or an index that cannot be parsed): `x[3 4]`
        // is a faulty index, not the variable x annotated with a unit.
        if (!rhs && pendingExpressionError_.empty()) {
            size_t lastBracket = rhsStr.rfind(']');
            if (lastBracket == rhsStr.size() - 1) {
                size_t openBracket = rhsStr.rfind('[');
                if (openBracket != std::string::npos) {
                    std::string noUnits = trim(rhsStr.substr(0, openBracket));
                    if (!noUnits.empty()) {
                        rhs = parseExpression(noUnits, lineNum);
                        if (rhs) {
                            // Found valid expression after stripping units
                            units = rhsStr.substr(openBracket + 1, lastBracket - openBracket - 1);
                        }
                    }
                }
            }
        }
        
        if (!lhs || !rhs) {
            // Propagate a specific message from parseFunctionCall if available.
            if (errorMsg && !expressionErrorText().empty())
                *errorMsg = expressionErrorText();
            return nullptr;
        }
        
        auto stmt = std::make_shared<Statement>();
        if (op == ":=") {
            stmt->node = Assignment{lhs, rhs};
        } else {
            stmt->node = Equation{lhs, rhs, units, inlineComment};
        }
        stmt->sourceLineNumber = lineNum;
        return stmt;
    }

    // Split a line on semicolons that are outside of strings, comments, and
    // parentheses.  EES uses `;` as a statement separator on a single line,
    // e.g.  `m = 86/12 ; n = 14/1 ; f = 0.07`.
    // Returns a vector of sub-strings (trimmed); empty segments are dropped.
    std::vector<std::string> splitOnSemicolons(const std::string& line) {
        std::vector<std::string> parts;
        std::string current;
        bool inSingleQuote = false;   // '...'  string literal
        bool inDoubleQuote = false;   // "..." comment / units
        int  parenDepth    = 0;
        int  bracketDepth  = 0;

        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];

            // Toggle single-quote string
            if (c == '\'' && !inDoubleQuote) {
                inSingleQuote = !inSingleQuote;
                current += c;
                continue;
            }
            if (inSingleQuote) { current += c; continue; }

            // Toggle double-quote comment/units
            if (c == '"' && !inSingleQuote) {
                inDoubleQuote = !inDoubleQuote;
                current += c;
                continue;
            }
            if (inDoubleQuote) { current += c; continue; }

            // Track parentheses and brackets
            if (c == '(' || c == '[') { parenDepth++; current += c; continue; }
            if (c == ')' || c == ']') { if (parenDepth > 0) parenDepth--; current += c; continue; }

            // Brace comment (skip entire block on this line)
            if (c == '{') { bracketDepth++; current += c; continue; }
            if (c == '}') { if (bracketDepth > 0) bracketDepth--; current += c; continue; }
            if (bracketDepth > 0) { current += c; continue; }

            // C-style line comment — rest of line is consumed
            if (c == '/' && i + 1 < line.size() && line[i + 1] == '/') {
                // Keep the rest as part of the current segment (so that
                // removeInlineComments can strip it later).
                current += line.substr(i);
                break;
            }

            // Semicolon at top level → split here
            if (c == ';' && parenDepth == 0) {
                std::string t = trim(current);
                if (!t.empty()) parts.push_back(t);
                current.clear();
                continue;
            }

            current += c;
        }
        std::string t = trim(current);
        if (!t.empty()) parts.push_back(t);
        return parts;
    }

    std::string removeInlineComments(const std::string& line, std::string* extractedQuoteComment = nullptr) {
        std::string result;
        bool inString = false;
        bool inQuoteComment = false;
        bool inBraceComment = false;
        int parenDepth = 0;
        
        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];
            
            if (inBraceComment) {
                if (c == '}') inBraceComment = false;
                continue;
            }
            
            if (c == '\'') {
                inString = !inString;
                result += c;
                continue;
            }
            
            if (inString) {
                result += c;
                continue;
            }
            
            if (c == '"') {
                // Check if this is a units annotation
                if (i + 1 < line.size() && line[i + 1] == '[') {
                    // Keep units annotation
                    size_t end = line.find("]\"", i);
                    if (end != std::string::npos) {
                        result += line.substr(i, end - i + 2);
                        i = end + 1;
                        continue;
                    }
                }
                // Skip quote comment to end of line or next quote
                size_t end = line.find('"', i + 1);
                if (end != std::string::npos) {
                    // Extract the first "" comment if requested
                    if (extractedQuoteComment && extractedQuoteComment->empty()) {
                        *extractedQuoteComment = line.substr(i + 1, end - i - 1);
                    }
                    i = end;
                }
                continue;
            }
            
            if (c == '{') {
                // Skip brace comment
                size_t end = line.find('}', i + 1);
                if (end != std::string::npos) {
                    i = end;
                }
                continue;
            }
            
            if (c == '/' && i + 1 < line.size() && line[i + 1] == '/') {
                // C-style comment - skip rest of line
                break;
            }
            
            result += c;
        }
        
        return result;
    }
    
    int findMainOperator(const std::string& expr, std::string& op) {
        int depth = 0;  // Parenthesis/bracket depth
        bool inString = false;
        
        for (size_t i = 0; i < expr.size(); ++i) {
            char c = expr[i];
            
            if (c == '\'') {
                inString = !inString;
                continue;
            }
            if (inString) continue;
            
            if (c == '(' || c == '[') depth++;
            else if (c == ')' || c == ']') depth--;
            else if (depth == 0) {
                if (c == ':' && i + 1 < expr.size() && expr[i+1] == '=') {
                    op = ":=";
                    return static_cast<int>(i);
                }
                if (c == '=' && (i == 0 || expr[i-1] != ':')) {
                    // Reject '==' — EES does not have a comparison equality operator
                    if (i + 1 < expr.size() && expr[i+1] == '=') {
                        op = "==";
                        return static_cast<int>(i);
                    }
                    op = "=";
                    return static_cast<int>(i);
                }
            }
        }
        return -1;
    }

    // Validate expression syntax before parsing.
    // Returns true if valid; sets errorMsg on failure.
    bool validateExpressionSyntax(const std::string& expr, std::string& errorMsg) {
        bool inString = false;
        int depth = 0;

        for (size_t i = 0; i < expr.size(); ++i) {
            char c = expr[i];
            if (c == '\'') { inString = !inString; continue; }
            if (inString) continue;
            if (c == '(' || c == '[') { depth++; continue; }
            if (c == ')' || c == ']') {
                depth--;
                if (depth < 0) {
                    errorMsg = std::string("Unmatched '") + c + "' — extra closing parenthesis or bracket";
                    return false;
                }
                continue;
            }

            // Detect consecutive arithmetic operators: ++, --, +-, -+, +++, etc.
            if ((c == '+' || c == '-') && depth == 0) {
                // Look ahead past optional whitespace for another +/-
                size_t j = i + 1;
                while (j < expr.size() && (expr[j] == ' ' || expr[j] == '\t')) j++;
                if (j < expr.size() && (expr[j] == '+' || expr[j] == '-')) {
                    // Allow if the first operator follows a binary operator (* / ^ =)
                    // which makes the +/- unary.  E.g. T * -2 is fine.
                    // Scan backwards past whitespace to find preceding token.
                    int k = static_cast<int>(i) - 1;
                    while (k >= 0 && (expr[k] == ' ' || expr[k] == '\t')) k--;
                    if (k >= 0) {
                        char prev = expr[k];
                        if (prev == '*' || prev == '/' || prev == '^' ||
                            prev == '=' || prev == '(' || prev == '[' || prev == ',') {
                            // First +/- is unary after a binary op — still reject
                            // double sign: use parentheses instead.
                        } else {
                            // Two consecutive additive operators like T + -2 or T+++2
                            errorMsg = "Consecutive operators '" + std::string(1, c)
                                     + std::string(1, expr[j])
                                     + "' are not allowed (use parentheses instead)";
                            return false;
                        }
                    }
                    // At the very start of the expression it could be unary too
                    // e.g. +-2 — still reject.
                    errorMsg = "Consecutive operators '" + std::string(1, c)
                             + std::string(1, expr[j])
                             + "' are not allowed (use parentheses instead)";
                    return false;
                }
            }

            // Detect consecutive multiplicative operators: **, //, */, /*
            if ((c == '*' || c == '/') && depth == 0) {
                size_t j = i + 1;
                while (j < expr.size() && (expr[j] == ' ' || expr[j] == '\t')) j++;
                if (j < expr.size() && (expr[j] == '*' || expr[j] == '/')) {
                    errorMsg = "Invalid operator '" + std::string(1, c) + std::string(1, expr[j]) + "'";
                    return false;
                }
            }
        }
        if (depth > 0) {
            errorMsg = "Unmatched '(' — missing closing parenthesis ')'";
            return false;
        }
        return true;
    }
    
    ExprPtr parseExpression(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        if (trimmed.empty()) return nullptr;
        
        // Parse as comparison expression (lowest precedence)
        return parseComparison(trimmed, lineNum);
    }
    
    ExprPtr parseComparison(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        
        // Look for comparison operators at depth 0: >=, <=, >, <
        // Scan left-to-right (left associative)
        int depth = 0;
        bool inString = false;
        for (size_t i = 0; i < trimmed.size(); ++i) {
            if (trimmed[i] == '\'') inString = !inString;
            if (inString) continue;
            if (trimmed[i] == '(' || trimmed[i] == '[') depth++;
            if (trimmed[i] == ')' || trimmed[i] == ']') depth--;
            if (depth != 0) continue;
            
            std::string op;
            if (i + 1 < trimmed.size() && trimmed[i] == '>' && trimmed[i+1] == '=') op = ">=";
            else if (i + 1 < trimmed.size() && trimmed[i] == '<' && trimmed[i+1] == '=') op = "<=";
            else if (trimmed[i] == '>' && (i == 0 || trimmed[i-1] != '<')) op = ">";
            else if (trimmed[i] == '<' && (i + 1 >= trimmed.size() || trimmed[i+1] != '=')) op = "<";
            
            if (!op.empty()) {
                std::string left = trim(trimmed.substr(0, i));
                std::string right = trim(trimmed.substr(i + op.size()));
                if (!left.empty() && !right.empty()) {
                    auto leftExpr = parseAdditive(left, lineNum);
                    auto rightExpr = parseAdditive(right, lineNum);
                    if (leftExpr && rightExpr) {
                        return makeBinaryOp(op, leftExpr, rightExpr, lineNum);
                    }
                }
            }
        }
        
        return parseAdditive(trimmed, lineNum);
    }
    
    ExprPtr parseAdditive(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        
        // Find + or - at depth 0 (right to left for left associativity)
        for (int i = static_cast<int>(trimmed.size()) - 1; i >= 0; --i) {
            if (trimmed[i] == '+' || trimmed[i] == '-') {
                int depth = 0;
                bool inString = false;
                
                // Check depth
                for (int j = 0; j < i; ++j) {
                    if (trimmed[j] == '\'') inString = !inString;
                    if (inString) continue;
                    if (trimmed[j] == '(' || trimmed[j] == '[') depth++;
                    if (trimmed[j] == ')' || trimmed[j] == ']') depth--;
                }
                
                if (depth == 0 && i > 0) {
                    // Skip if part of scientific notation (e.g., 1E+5 or 1.2E-5)
                    if ((trimmed[i] == '+' || trimmed[i] == '-') && i > 0) {
                        char prev = std::tolower(static_cast<unsigned char>(trimmed[i-1]));
                        if (prev == 'e' && i > 1) {
                            char prev2 = trimmed[i-2];
                            if (std::isdigit(static_cast<unsigned char>(prev2)) || prev2 == '.') {
                                continue; // This is the sign of an exponent, not a binary operator
                            }
                        }
                    }

                    // Check if this is a unary minus (preceded by operator or nothing)
                    char prev = trimmed[i-1];
                    if (prev == '(' || prev == '[' || prev == '+' || prev == '-' || 
                        prev == '*' || prev == '/' || prev == '^' || prev == '=' || prev == ',') {
                        continue;  // This is unary, not binary
                    }
                    
                    std::string left = trimmed.substr(0, i);
                    std::string right = trimmed.substr(i + 1);
                    std::string op(1, trimmed[i]);
                    
                    auto leftExpr = parseAdditive(left, lineNum);
                    auto rightExpr = parseMultiplicative(right, lineNum);
                    
                    if (leftExpr && rightExpr) {
                        return makeBinaryOp(op, leftExpr, rightExpr, lineNum);
                    }
                }
            }
        }
        
        return parseMultiplicative(trimmed, lineNum);
    }
    
    ExprPtr parseMultiplicative(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        
        // Find * or / at depth 0 (right to left)
        for (int i = static_cast<int>(trimmed.size()) - 1; i >= 0; --i) {
            if (trimmed[i] == '*' || trimmed[i] == '/') {
                int depth = 0;
                bool inString = false;
                
                for (int j = 0; j < i; ++j) {
                    if (trimmed[j] == '\'') inString = !inString;
                    if (inString) continue;
                    if (trimmed[j] == '(' || trimmed[j] == '[') depth++;
                    if (trimmed[j] == ')' || trimmed[j] == ']') depth--;
                }
                
                if (depth == 0) {
                    std::string left = trimmed.substr(0, i);
                    std::string right = trimmed.substr(i + 1);
                    std::string op(1, trimmed[i]);
                    
                    auto leftExpr = parseMultiplicative(left, lineNum);
                    auto rightExpr = parsePower(right, lineNum);
                    
                    if (leftExpr && rightExpr) {
                        return makeBinaryOp(op, leftExpr, rightExpr, lineNum);
                    }
                }
            }
        }
        
        return parsePower(trimmed, lineNum);
    }
    
    ExprPtr parsePower(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        
        // Find ^ at depth 0 (left to right for right associativity)
        for (size_t i = 0; i < trimmed.size(); ++i) {
            if (trimmed[i] == '^') {
                int depth = 0;
                bool inString = false;
                
                for (size_t j = 0; j < i; ++j) {
                    if (trimmed[j] == '\'') inString = !inString;
                    if (inString) continue;
                    if (trimmed[j] == '(' || trimmed[j] == '[') depth++;
                    if (trimmed[j] == ')' || trimmed[j] == ']') depth--;
                }
                
                if (depth == 0) {
                    std::string left = trimmed.substr(0, i);
                    std::string right = trimmed.substr(i + 1);
                    
                    auto leftExpr = parseUnary(left, lineNum);
                    auto rightExpr = parsePower(right, lineNum);  // Right associative
                    
                    if (leftExpr && rightExpr) {
                        return makeBinaryOp("^", leftExpr, rightExpr, lineNum);
                    }
                }
            }
        }
        
        return parseUnary(trimmed, lineNum);
    }
    
    ExprPtr parseUnary(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        
        if (trimmed.empty()) return nullptr;
        
        // Handle unary minus/plus
        if (trimmed[0] == '-' || trimmed[0] == '+') {
            std::string op(1, trimmed[0]);
            std::string rest = trimmed.substr(1);
            auto operand = parsePrimary(rest, lineNum);
            if (operand) {
                return makeUnaryOp(op, operand, lineNum);
            }
        }
        
        return parsePrimary(trimmed, lineNum);
    }
    
    ExprPtr parsePrimary(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        
        if (trimmed.empty()) return nullptr;
        
        // Parenthesized expression
        if (trimmed.front() == '(' && trimmed.back() == ')') {
            // Check if the parens match
            int depth = 0;
            bool matches = true;
            for (size_t i = 0; i < trimmed.size() - 1; ++i) {
                if (trimmed[i] == '(') depth++;
                else if (trimmed[i] == ')') depth--;
                if (depth == 0 && i > 0) {
                    matches = false;
                    break;
                }
            }
            if (matches) {
                return parseExpression(trimmed.substr(1, trimmed.size() - 2), lineNum);
            }
        }
        
        // String literal (the closing quote ends it)
        if (trimmed.front() == '\'') {
            size_t end = trimmed.find('\'', 1);
            if (end != std::string::npos) {
                std::string value = trimmed.substr(1, end - 1);
                if (end + 1 < trimmed.size()) {
                    setExpressionError(unexpectedTextMessage(trimmed.substr(end + 1), "the string", value));
                    return nullptr;
                }
                return makeString(value, lineNum);
            }
        }
        
        // Number
        if (isNumber(trimmed)) {
            try {
                double value = std::stod(trimmed);
                return makeNumber(value, lineNum);
            } catch (...) {
                return nullptr;
            }
        }
        
        // A unit annotation in brackets ends the term it follows: a number or a parenthesized
        // group (3600[s], 0.05 [kg/s], (25 + 273.15) [K]), a function call or an array element
        // (max(a, b) [kW], x[1] [K]). As for the annotation of a whole equation, the unit is
        // documentation only (values are not converted): the term is kept and the annotation
        // dropped. After a call or an element the annotation must hold a unit: a number
        // between brackets (f(x)[2], A[1][2]) is an index, not a unit, and is refused below.
        const size_t annotationPos = trimmed.rfind('[');
        if (annotationPos != std::string::npos && annotationPos > 0 && trimmed.back() == ']' &&
            trimmed.find(']', annotationPos) == trimmed.size() - 1) {
            const std::string annotated = trim(trimmed.substr(0, annotationPos));
            if (isNumber(annotated) || isFullyParenthesized(annotated))
                return parsePrimary(annotated, lineNum);
            const std::string unit = trim(trimmed.substr(annotationPos + 1, trimmed.size() - annotationPos - 2));
            if (!annotated.empty() && (annotated.back() == ')' || annotated.back() == ']') &&
                !unit.empty() && !isNumber(unit))
                return parseCallOrElement(annotated, lineNum);
        }

        // Function call or array element
        if (trimmed.find_first_of("([") != std::string::npos)
            return parseCallOrElement(trimmed, lineNum);

        // Simple variable
        if (isIdentifier(trimmed)) {
            return makeVariable(trimmed, lineNum);
        }
        
        // Say why a name is refused (a character EES does not take in a name, a Latin-1 file)
        const std::string why = invalidNameMessage(trimmed);
        if (!why.empty() && pendingNameHint_.empty()) pendingNameHint_ = why;
        return nullptr;
    }

    // `expr` starts with a function name followed by '(', or with a variable name followed by
    // '[' (whichever comes first): a function call or an array element. The text after the
    // closing bracket is refused here (see parseFunctionCall, parseArrayAccess).
    ExprPtr parseCallOrElement(const std::string& expr, int lineNum) {
        const size_t parenPos = expr.find('(');
        const size_t bracketPos = expr.find('[');
        if (parenPos != std::string::npos && (bracketPos == std::string::npos || parenPos < bracketPos))
            return parseFunctionCall(expr, lineNum);
        if (bracketPos != std::string::npos)
            return parseArrayAccess(expr, lineNum);
        return nullptr;
    }

    // Position of the bracket that closes the '(' or '[' at s[open], the brackets of both kinds
    // being counted and the 'strings' skipped; npos when it is never closed.
    static size_t findClosingBracket(const std::string& s, size_t open) {
        int depth = 0;
        bool inString = false;
        for (size_t i = open; i < s.size(); ++i) {
            const char c = s[i];
            if (c == '\'') { inString = !inString; continue; }
            if (inString) continue;
            if (c == '(' || c == '[') depth++;
            else if (c == ')' || c == ']') {
                if (--depth == 0) return i;
            }
        }
        return std::string::npos;
    }

    // The text found after the end of a call, an array element or a string, where nothing may
    // follow but a unit annotation (CS-BUG-TRAILING-TEXT): it was ignored.
    static std::string unexpectedTextMessage(const std::string& text, const std::string& what,
                                             const std::string& after) {
        std::string t = text;
        t.erase(0, t.find_first_not_of(" \t"));
        return "Unexpected text '" + t + "' after " + what + " '" + after + "'";
    }

    /// Whether `arg`, an argument of CONVERT or CONVERTTEMP as written, is the name of a unit given without quotes
    /// (CONVERT(lbm/hr, kg/s), CONVERTTEMP(C, K, T)), which EES allows. It is read from the text, before it becomes
    /// an expression (lbm/hr is no division of two variables, kJ/kg-K is no subtraction): a text made of the
    /// characters of a unit string only (names, digits, - * / ^ . and parentheses) with at least one letter. A string
    /// ('..', u$, a string function) and a number or an expression with other operators are not.
    static bool isUnquotedUnitText(const std::string& arg) {
        bool hasLetter = false;
        for (const char c : arg) {
            const unsigned char u = static_cast<unsigned char>(c);
            if (std::isalpha(u) || c == '_' || u >= 0x80) hasLetter = true;
            else if (!std::isdigit(u) && std::string(" \t-*/^.()").find(c) == std::string::npos) return false;
        }
        return hasLetter;
    }

    /// The aggregate built-ins take a named argument only in the EES indexed form `SUM(expr, i=first, last)`
    /// (and PRODUCT): an error message for any other use, empty when `args` (the argument texts of the call
    /// `call`) is fine. AVERAGE, STDDEV and SUM2D have no indexed form: the named argument would be taken
    /// for a variable of the model. The indexed form is noted: resolveCompileTimeValues expands it
    /// (resolveIndexedAggregate).
    std::string indexedAggregateError(const std::string& name, const std::vector<std::string>& args,
                                      const std::string& call) {
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        const bool indexable = lower == "sum" || lower == "product";
        if (!indexable && lower != "average" && lower != "stddev" && lower != "sum2d") return "";
        std::vector<bool> named;
        for (const auto& arg : args) {
            const size_t eq = arg.find('=');
            named.push_back(eq != std::string::npos && eq > 0 && isIdentifier(trim(arg.substr(0, eq))));
        }
        if (std::find(named.begin(), named.end(), true) == named.end()) return "";
        std::string upper = lower;
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        if (!indexable)
            return upper + ": the indexed form " + upper + "(expression, i=first, last) is implemented for SUM and "
                   "PRODUCT only (found '" + trim(call) + "'); for the other aggregates use an array range, as in " +
                   upper + "(x[1..n])";
        if (args.size() == 3 && !named[0] && named[1] && !named[2]) {
            indexedAggregateSeen_ = true;
            return "";
        }
        return upper + ": the indexed form is " + upper + "(expression, i=first, last), as in " + upper +
               "(x[i], i=1, 4); found '" + trim(call) + "'";
    }

    ExprPtr parseFunctionCall(const std::string& expr, int lineNum) {
        size_t parenStart = expr.find('(');
        if (parenStart == std::string::npos) return nullptr;
        
        std::string name = trim(expr.substr(0, parenStart));
        if (name.empty()) return nullptr;   // "(a + b) c": a group, not a call
        // "2*f(x)" is no call of a function "2*f": when f(x) is refused, the whole text must be refused with
        // the reason, not accepted as a call of an unknown function
        if (name.find_first_of("+-*/^<>=,[]") != std::string::npos) return nullptr;
        
        // The parenthesis that closes the call: nothing but a unit annotation (handled by
        // parsePrimary) may follow it
        const size_t parenEnd = findClosingBracket(expr, parenStart);
        if (parenEnd == std::string::npos || expr[parenEnd] != ')') return nullptr;
        if (parenEnd + 1 < expr.size()) {
            const std::string rest = trim(expr.substr(parenEnd + 1));
            if (!rest.empty()) {
                setExpressionError(unexpectedTextMessage(rest, "the call", trim(expr.substr(0, parenEnd + 1))));
                return nullptr;
            }
        }
        
        std::string argsStr = expr.substr(parenStart + 1, parenEnd - parenStart - 1);
        
        auto func = std::make_shared<Expression>();
        FunctionCall fc;
        fc.name = name;
        
        // Check if this is a thermophysical function (first arg is fluid name)
        bool isThermo = isEESThermophysicalFunction(name);
        bool isFirstPositionalArg = true;
        
        // Parse arguments
        std::vector<std::string> argStrings = splitArguments(argsStr);
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        // CONVERT(from, to) and CONVERTTEMP(from, to, value): `from` and `to` are unit names, the quotes of which are optional
        const bool unitArguments = lower == "convert" || lower == "converttemp";

        const std::string aggregateError = indexedAggregateError(name, argStrings, expr);
        if (!aggregateError.empty()) {
            setExpressionError(aggregateError);
            return nullptr;
        }
        
        for (size_t argIndex = 0; argIndex < argStrings.size(); ++argIndex) {
            std::string arg = trim(argStrings[argIndex]);
            
            // Check for named argument (e.g., T=value)
            size_t eqPos = arg.find('=');
            if (eqPos != std::string::npos && eqPos > 0) {
                // Check if the part before = is an identifier (not an expression)
                std::string paramName = trim(arg.substr(0, eqPos));
                std::string paramValue = trim(arg.substr(eqPos + 1));
                
                if (isIdentifier(paramName)) {
                    if (isPositionalOnlyBuiltin(name)) {
                        setExpressionError(namedArgumentRefusedMessage(name, arg, paramName, paramValue));
                        return nullptr;
                    }
                    auto valueExpr = parseExpression(paramValue, lineNum);
                    if (valueExpr) {
                        fc.namedArgs.push_back({paramName, valueExpr});
                        continue;
                    }
                    // Value expression is malformed — report a specific error and
                    // abort the whole function call so the equation fails to parse
                    // instead of silently continuing with a missing argument.
                    setExpressionError("Could not parse value for named argument '"
                        + paramName + "': '" + paramValue + "'");
                    return nullptr;
                }
            }
            
            // Regular argument
            if (unitArguments && argIndex < 2 && isUnquotedUnitText(arg)) {
                fc.args.push_back(makeString(arg, lineNum));   // the unit as written, checked below like a quoted one
                continue;
            }
            auto argExpr = parseExpression(arg, lineNum);
            if (!argExpr) {
                // Never skip an argument: the call would be evaluated with fewer
                // arguments and the model solved with a wrong value.
                setExpressionError("Could not parse argument " + std::to_string(argIndex + 1)
                    + " of '" + name + "': '" + arg + "'");
                return nullptr;
            }
            // For thermophysical functions, convert unquoted first argument to string literal
            // In EES, the first argument is always the fluid name (e.g., humrat(airH2O, ...) 
            // or enthalpy(R134a, ...)) and quotes are optional
            // Exception: string variables (ending with $) should remain as variables
            // because they need to be resolved at runtime
            if (isThermo && isFirstPositionalArg && argExpr->is<Variable>()) {
                const Variable& var = argExpr->as<Variable>();
                // Only convert if: simple variable (no array access) AND not a string variable (not ending with $)
                if (var.indices.empty() && !var.name.empty() && var.name.back() != '$') {
                    // Convert the variable to a string literal (fluid name)
                    argExpr = makeString(var.name, lineNum);
                }
            }
            fc.args.push_back(argExpr);
            isFirstPositionalArg = false;
        }
        
        // CONVERT('from', 'to') and CONVERTTEMP('from', 'to', value) with literal units: an unknown unit, units of
        // different dimensions (CONVERT) or a unit that is not a temperature scale (CONVERTTEMP) are an error of
        // this line, found before the solve (a unit string computed at run time is checked when evaluated)
        if (fc.args.size() == 2 || fc.args.size() == 3) {
            const bool fromLiteral = fc.args[0]->is<StringLiteral>();
            const bool toLiteral = fc.args[1]->is<StringLiteral>();
            if (fromLiteral || toLiteral) {
                try {
                    if (lower == "convert" && fc.args.size() == 2 && fromLiteral && toLiteral) {
                        unitConversionFactor(fc.args[0]->as<StringLiteral>().value, fc.args[1]->as<StringLiteral>().value);
                    } else if (lower == "converttemp" && fc.args.size() == 3) {
                        if (fromLiteral) checkTemperatureUnit(fc.args[0]->as<StringLiteral>().value);
                        if (toLiteral) checkTemperatureUnit(fc.args[1]->as<StringLiteral>().value);
                    }
                } catch (const std::runtime_error& e) {
                    clearExpressionError();   // what is pending comes from a rejected alternative reading of the text
                    setExpressionError(e.what());
                    return nullptr;
                }
            }
        }
        
        func->node = fc;
        func->sourceLineNumber = lineNum;
        return func;
    }
    
    ExprPtr parseArrayAccess(const std::string& expr, int lineNum) {
        size_t bracketStart = expr.find('[');
        if (bracketStart == std::string::npos) return nullptr;
        
        std::string name = trim(expr.substr(0, bracketStart));
        
        // Ensure the name is a valid identifier (to distinguish from units like 10[m])
        if (!isIdentifier(name)) return nullptr;

        // The bracket that closes the element: nothing but a unit annotation (handled by
        // parsePrimary) may follow it
        const size_t bracketEnd = findClosingBracket(expr, bracketStart);
        if (bracketEnd == std::string::npos || expr[bracketEnd] != ']') return nullptr;

        // A string array (U$[2]) is not supported, by decision: refused here, in the main program, in
        // FUNCTION/PROCEDURE bodies and in DUPLICATE alike, instead of being stored as the number 0
        // (main program) or failing at the read (bodies)
        if (name.back() == '$') {
            setExpressionError("String arrays are not supported: '" + trim(expr.substr(0, bracketEnd + 1)) +
                "' is an element of the string array " + name + ". Put the strings in a lookup table (one "
                "string per row, file <model>-<table>.csv), keep the row number in the array instead (a "
                "numeric array) and read the string with LOOKUP$('table', row, col), which can be used "
                "wherever a string is expected, for instance as the fluid of a property call: "
                "enthalpy(LOOKUP$('table', i, 1), T=T)");
            return nullptr;
        }

        if (bracketEnd + 1 < expr.size()) {
            const std::string rest = trim(expr.substr(bracketEnd + 1));
            if (!rest.empty()) {
                setExpressionError(unexpectedTextMessage(rest, "the array element", trim(expr.substr(0, bracketEnd + 1))));
                return nullptr;
            }
        }
        
        std::string indexStr = expr.substr(bracketStart + 1, bracketEnd - bracketStart - 1);
        
        std::vector<ExprPtr> indices;
        std::vector<std::string> indexStrings = splitArguments(indexStr);
        
        for (size_t i = 0; i < indexStrings.size(); ++i) {
            auto indexExpr = parseIndexOrRange(indexStrings[i], lineNum);
            if (!indexExpr) {
                setExpressionError("Could not parse index " + std::to_string(i + 1)
                    + " of '" + name + "': '" + trim(indexStrings[i]) + "'");
                return nullptr;
            }
            indices.push_back(indexExpr);
        }
        
        ExprPtr element = makeArrayVariable(name, indices, lineNum);
        if (unitIndexOf(element->as<Variable>())) unitIndexSeen_ = true;
        return element;
    }

    // An index, or a range `a..b` of indices (EES array range notation, CS-GAP-ARRAY-ARG): the binary
    // operator "..", which resolveCompileTimeValues() expands into the elements once the whole file is
    // parsed. nullptr when the text is neither.
    ExprPtr parseIndexOrRange(const std::string& text, int lineNum) {
        size_t dots = std::string::npos;
        int depth = 0;
        for (size_t i = 0; i + 1 < text.size() && dots == std::string::npos; ++i) {
            if (text[i] == '(' || text[i] == '[') depth++;
            else if (text[i] == ')' || text[i] == ']') depth--;
            else if (depth == 0 && text[i] == '.' && text[i + 1] == '.') dots = i;
        }
        if (dots == std::string::npos) return parseExpression(text, lineNum);
        const std::string first = trim(text.substr(0, dots));
        const std::string last = trim(text.substr(dots + 2));
        if (first.empty() || last.empty() || last[0] == '.' || first.back() == '.') return nullptr;   // 1...3
        auto firstExpr = parseExpression(first, lineNum);
        auto lastExpr = firstExpr ? parseExpression(last, lineNum) : nullptr;
        if (!lastExpr) return nullptr;
        rangeSeen_ = true;
        return makeBinaryOp("..", firstExpr, lastExpr, lineNum);
    }

    std::vector<std::string> splitArguments(const std::string& argsStr) {
        std::vector<std::string> result;
        std::string current;
        int depth = 0;
        bool inString = false;
        
        for (char c : argsStr) {
            if (c == '\'') inString = !inString;
            
            if (!inString) {
                if (c == '(' || c == '[') depth++;
                else if (c == ')' || c == ']') depth--;
                else if (c == ',' && depth == 0) {
                    result.push_back(current);
                    current.clear();
                    continue;
                }
            }
            current += c;
        }
        
        // No argument at all ("" or blanks: f(), f( )) gives an empty list; after a comma
        // the last piece is kept even if empty, so that f(a, ) is not read as f(a).
        if (!result.empty() || !trim(current).empty()) {
            result.push_back(current);
        }
        
        return result;
    }
    
    bool isNumber(const std::string& s) {
        if (s.empty()) return false;
        
        size_t i = 0;
        // Leading +/- is handled by parseUnary, not here
        
        bool hasDigits = false;
        while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
            hasDigits = true;
            i++;
        }
        
        if (i < s.size() && s[i] == '.') {
            i++;
            while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
                hasDigits = true;
                i++;
            }
        }
        
        if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
            i++;
            if (i < s.size() && (s[i] == '+' || s[i] == '-')) i++;
            bool hasExpDigits = false;
            while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
                hasExpDigits = true;
                i++;
            }
            if (!hasExpDigits) return false;
        }
        
        return hasDigits && i == s.size();
    }
    
    bool isIdentifier(const std::string& s) {
        return firstInvalidNameChar(s) == std::string::npos;
    }
    
    std::string trim(const std::string& s) {
        size_t start = 0;
        while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) start++;
        
        size_t end = s.size();
        while (end > start && std::isspace(static_cast<unsigned char>(s[end-1]))) end--;
        
        return s.substr(start, end - start);
    }
};

// ============================================================================
// EESParser Public Interface
// ============================================================================

EESParser::EESParser() : pImpl(std::make_unique<Impl>()) {}
EESParser::~EESParser() = default;

ParseResult EESParser::parse(const std::string& source, const std::string& filename) {
    return pImpl->parse(source, filename);
}

ParseResult EESParser::parseFile(const std::string& filepath) {
    auto content = readFile(filepath);
    if (!content) {
        ParseResult result;
        result.success = false;
        result.errors.push_back({0, 0, "Could not read file: " + filepath, ""});
        return result;
    }
    return parse(*content, filepath);
}

std::string EESParser::getLastError() const {
    return pImpl->getLastError();
}

// ============================================================================
// Utility Functions
// ============================================================================

std::optional<double> evaluateConstantExpression(
    const ExprPtr& expr, const std::function<std::optional<double>(const std::string&)>& lookup,
    std::string* unknown) {
    if (!expr) return std::nullopt;
    if (expr->is<NumberLiteral>()) return expr->as<NumberLiteral>().value;
    if (expr->is<Variable>()) {
        const auto& var = expr->as<Variable>();
        std::optional<double> value;
        if (var.indices.empty() && lookup) value = lookup(var.name);
        if (!value && unknown && unknown->empty()) *unknown = var.flattenedName();
        return value;
    }
    if (expr->is<UnaryOp>()) {
        const auto& op = expr->as<UnaryOp>();
        const auto operand = evaluateConstantExpression(op.operand, lookup, unknown);
        if (!operand || (op.op != "-" && op.op != "+")) return std::nullopt;
        return op.op == "-" ? -*operand : *operand;
    }
    if (expr->is<BinaryOp>()) {
        const auto& op = expr->as<BinaryOp>();
        const auto left = evaluateConstantExpression(op.left, lookup, unknown);
        if (!left) return std::nullopt;
        const auto right = evaluateConstantExpression(op.right, lookup, unknown);
        if (!right) return std::nullopt;
        double value;
        if (op.op == "+") value = *left + *right;
        else if (op.op == "-") value = *left - *right;
        else if (op.op == "*") value = *left * *right;
        else if (op.op == "/") value = *left / *right;
        else if (op.op == "^") value = std::pow(*left, *right);
        else return std::nullopt;   // a comparison, a range a..b
        if (!std::isfinite(value)) return std::nullopt;
        return value;
    }
    return std::nullopt;   // a function call, a string
}

std::optional<int> asIntegerBound(double value) {
    const double rounded = std::round(value);
    if (!std::isfinite(value) || std::fabs(rounded) > 1e9 ||
        std::fabs(value - rounded) > 1e-9 * std::max(1.0, std::fabs(value)))
        return std::nullopt;
    return static_cast<int>(rounded);
}

std::optional<std::string> readFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return std::nullopt;
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

void collectVariables(const ExprPtr& expr, std::vector<std::string>& vars) {
    if (!expr) return;
    
    std::visit([&vars](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        
        if constexpr (std::is_same_v<T, Variable>) {
            if (node.indices.empty()) {
                vars.push_back(node.name);
            } else {
                vars.push_back(node.flattenedName());
            }
            for (const auto& idx : node.indices) {
                collectVariables(idx, vars);
            }
        } else if constexpr (std::is_same_v<T, UnaryOp>) {
            collectVariables(node.operand, vars);
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            collectVariables(node.left, vars);
            collectVariables(node.right, vars);
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            for (const auto& arg : node.args) {
                collectVariables(arg, vars);
            }
            for (const auto& [name, arg] : node.namedArgs) {
                collectVariables(arg, vars);
            }
        }
    }, expr->node);
}

std::string astToString(const ExprPtr& expr) {
    if (!expr) return "<null>";
    
    return std::visit([](const auto& node) -> std::string {
        using T = std::decay_t<decltype(node)>;
        
        if constexpr (std::is_same_v<T, NumberLiteral>) {
            return std::to_string(node.value);
        } else if constexpr (std::is_same_v<T, StringLiteral>) {
            return "'" + node.value + "'";
        } else if constexpr (std::is_same_v<T, Variable>) {
            std::string result = node.name;
            if (!node.indices.empty()) {
                result += "[";
                for (size_t i = 0; i < node.indices.size(); ++i) {
                    if (i > 0) result += ",";
                    result += astToString(node.indices[i]);
                }
                result += "]";
            }
            return result;
        } else if constexpr (std::is_same_v<T, UnaryOp>) {
            return "(" + node.op + astToString(node.operand) + ")";
        } else if constexpr (std::is_same_v<T, BinaryOp>) {
            if (node.op == "..") return astToString(node.left) + ".." + astToString(node.right);   // an array range
            return "(" + astToString(node.left) + " " + node.op + " " + astToString(node.right) + ")";
        } else if constexpr (std::is_same_v<T, FunctionCall>) {
            std::string result = node.name + "(";
            bool first = true;
            for (const auto& arg : node.args) {
                if (!first) result += ", ";
                result += astToString(arg);
                first = false;
            }
            for (const auto& [name, arg] : node.namedArgs) {
                if (!first) result += ", ";
                result += name + "=" + astToString(arg);
                first = false;
            }
            result += ")";
            return result;
        }
        return "<unknown>";
    }, expr->node);
}

std::string astToString(const StmtPtr& stmt) {
    if (!stmt) return "<null>";
    
    return std::visit([](const auto& node) -> std::string {
        using T = std::decay_t<decltype(node)>;
        
        if constexpr (std::is_same_v<T, Equation>) {
            return astToString(node.lhs) + " = " + astToString(node.rhs);
        } else if constexpr (std::is_same_v<T, Assignment>) {
            return astToString(node.lhs) + " := " + astToString(node.rhs);
        } else if constexpr (std::is_same_v<T, ProcedureCall>) {
            return "CALL " + node.name + "(...)";
        } else if constexpr (std::is_same_v<T, FunctionDefinition>) {
            return "FUNCTION " + node.name + "(...)";
        } else if constexpr (std::is_same_v<T, ProcedureDefinition>) {
            return "PROCEDURE " + node.name + "(...)";
        } else if constexpr (std::is_same_v<T, Comment>) {
            return "// " + node.text;
        } else if constexpr (std::is_same_v<T, Directive>) {
            return "$" + node.name + " " + node.content;
        } else if constexpr (std::is_same_v<T, IfThenElse>) {
            return "IF " + astToString(node.condition) + " THEN ...";
        } else if constexpr (std::is_same_v<T, Duplicate>) {
            return "DUPLICATE " + node.iteratorVar + " = " + astToString(node.start) + ", " + astToString(node.end);
        } else if constexpr (std::is_same_v<T, RepeatUntil>) {
            return "REPEAT ... UNTIL " + astToString(node.condition);
        } else if constexpr (std::is_same_v<T, Label>) {
            return std::to_string(node.number) + ":";
        } else if constexpr (std::is_same_v<T, Goto>) {
            return "GOTO " + std::to_string(node.label);
        }
        return "<unknown>";
    }, stmt->node);
}

std::string astToString(const Program& program) {
    std::string result;
    for (const auto& stmt : program.statements) {
        result += astToString(stmt) + "\n";
    }
    return result;
}

// ============================================================================
// Variable::flattenedName implementation
// ============================================================================

std::string Variable::flattenedName() const {
    if (indices.empty()) return name;
    
    std::string result = name + "[";
    for (size_t i = 0; i < indices.size(); ++i) {
        if (i > 0) result += ",";
        if (indices[i]->is<NumberLiteral>()) {
            double val = indices[i]->as<NumberLiteral>().value;
            if (val == static_cast<int>(val)) {
                result += std::to_string(static_cast<int>(val));
            } else {
                result += std::to_string(val);
            }
        } else if (indices[i]->is<Variable>()) {
            result += indices[i]->as<Variable>().flattenedName();
        } else {
            result += "expr"; // Fallback for complex expressions in indices
        }
    }
    result += "]";
    return result;
}

}  // namespace coolsolve
