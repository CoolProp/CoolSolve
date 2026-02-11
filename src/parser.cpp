#include "coolsolve/parser.h"
#include <peglib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <set>

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
        bool inMultiLineQuote = false;
        bool inMultiLineQuoteStandalone = false;  // true if block started by a line that is exactly `"`
        int braceCommentDepth = 0; // Track nested { } comments across lines
        while (std::getline(stream, line)) {
            lineNumber++;
            
            // Handle multi-line quote comments " ... "
            if (inMultiLineQuote) {
                std::string innerTrimmed = trim(line);
                if (innerTrimmed == "\"") {
                    // Standalone `"` always ends a quote block (and starts one when not in block).
                    inMultiLineQuote = false;
                    result.commentCount++;
                    continue;
                }
                if (inMultiLineQuoteStandalone) {
                    // Block started with standalone `"` -> only a standalone `"` can close it.
                    continue;
                }
                // Classic style: block started with `"` on a line with no closing `"` -> first `"` on a line closes.
                size_t quoteEnd = line.find('"');
                if (quoteEnd != std::string::npos) {
                    inMultiLineQuote = false;
                    result.commentCount++;
                    line = line.substr(quoteEnd + 1);
                } else {
                    continue;
                }
            }

            // Strip nested brace comments that can span multiple lines.
            // This is done before any further processing so that code inside
            // `{ ... }` (including nested braces) is completely ignored.
            line = stripBraceComments(line, braceCommentDepth, result);
            if (braceCommentDepth > 0 && line.empty()) {
                // Still inside a block comment and nothing else on this line.
                continue;
            }

            // Skip empty lines
            if (isEmptyOrWhitespace(line)) continue;
            
            std::string trimmed = trim(line);
            if (trimmed.empty()) continue;

            // Start of multi-line quote: (1) line is exactly `"` (standalone delimiter), or
            // (2) line starts with `"` and has no closing `"` on the same line (classic EES style).
            if (trimmed == "\"") {
                inMultiLineQuote = true;
                inMultiLineQuoteStandalone = true;
                continue;
            }
            if (trimmed.front() == '"' && trimmed.find('"', 1) == std::string::npos) {
                inMultiLineQuote = true;
                inMultiLineQuoteStandalone = false;
                continue;
            }

            // Check for unsupported constructs (procedures, modules, etc.)
            std::string upper = trimmed;
            std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
            
            // Match keywords at start of line with word boundary
            auto isKeyword = [&](const std::string& kw) {
                if (upper.find(kw) != 0) return false;
                if (upper.size() == kw.size()) return true;
                char next = upper[kw.size()];
                return !std::isalnum(static_cast<unsigned char>(next)) && next != '_';
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
                    std::string bodyLine;
                    int startLine = lineNumber;
                    bool foundEnd = false;
                    while (std::getline(stream, bodyLine)) {
                        lineNumber++;
                        std::string trimmedBody = trim(bodyLine);
                        std::string upperBody = trimmedBody;
                        std::transform(upperBody.begin(), upperBody.end(), upperBody.begin(), ::toupper);
                        if (upperBody == "END") {
                            foundEnd = true;
                            break;
                        }
                        if (auto stmt = tryParseComment(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        } else if (auto stmt = tryParseDirective(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        } else if (auto stmt = tryParseEquationOrAssignment(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        } else if (auto stmt = tryParseProcedureCall(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        }
                    }
                    if (!foundEnd) {
                        result.errors.push_back({startLine, 0, "Function '" + name + "' missing END", line});
                    } else {
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
                    std::vector<StmtPtr> body;
                    std::string bodyLine;
                    int startLine = lineNumber;
                    bool foundEnd = false;
                    while (std::getline(stream, bodyLine)) {
                        lineNumber++;
                        std::string trimmedBody = trim(bodyLine);
                        std::string upperBody = trimmedBody;
                        std::transform(upperBody.begin(), upperBody.end(), upperBody.begin(), ::toupper);
                        if (upperBody == "END") {
                            foundEnd = true;
                            break;
                        }
                        if (auto stmt = tryParseComment(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        } else if (auto stmt = tryParseDirective(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        } else if (auto stmt = tryParseEquationOrAssignment(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        } else if (auto stmt = tryParseProcedureCall(bodyLine, lineNumber)) {
                            body.push_back(stmt);
                        }
                    }
                    if (!foundEnd) {
                        result.errors.push_back({startLine, 0, "Procedure '" + name + "' missing END", line});
                    } else {
                        result.program.statements.push_back(makeProcedureDefinition(name, inputs, outputs, body, startLine));
                    }
                    continue;
                }
            }

            // Handle Procedure Call
            if (isKeyword("CALL")) {
                if (auto stmt = tryParseProcedureCall(line, lineNumber)) {
                    result.program.statements.push_back(stmt);
                    continue;
                }
            }
            
            // Handle comments
            if (auto stmt = tryParseComment(line, lineNumber)) {
                result.program.statements.push_back(stmt);
                result.commentCount++;
                continue;
            }
            
            // Handle directives
            if (auto stmt = tryParseDirective(line, lineNumber)) {
                result.program.statements.push_back(stmt);
                result.directiveCount++;
                continue;
            }
            
            // Handle equations
            if (auto stmt = tryParseEquationOrAssignment(line, lineNumber)) {
                result.program.statements.push_back(stmt);
                result.equationCount++;
                continue;
            }
            
            // If we couldn't parse the line, record an error but continue
            // We don't set unsupportedFound here, as it might just be a line we don't understand yet
            result.errors.push_back({lineNumber, 0, "Could not parse line", line});
        }
        
        result.totalLines = lineNumber;
        // Success if no unsupported constructs AND we found something valid
        bool somethingParsed = result.equationCount > 0 || result.commentCount > 0 || result.directiveCount > 0 || !result.program.statements.empty();
        result.success = !unsupportedFound && somethingParsed;
        return result;
    }
    
    std::string getLastError() const { return lastError_; }

private:
    peg::parser parser_;
    bool grammarValid_ = false;
    std::string lastError_;
    
    void initializeGrammar() {
        // Use a simpler grammar approach - parse expressions directly
        grammarValid_ = true;  // We'll use manual parsing
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
    
    StmtPtr tryParseComment(const std::string& line, int lineNum) {
        std::string trimmed = trim(line);
        
        // Check for quote-delimited comment
        if (trimmed.size() >= 2 && trimmed.front() == '"' && trimmed.back() == '"') {
            // Find the position of the second quote (end of first comment)
            size_t firstQuoteEnd = trimmed.find('"', 1);
            if (firstQuoteEnd == trimmed.size() - 1) {
                // Only one pair of quotes - this is a pure comment
                std::string content = trimmed.substr(1, trimmed.size() - 2);
                return makeComment(content, true, lineNum);
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
            return makeComment(content, true, lineNum);
        }
        
        // Check for brace-delimited comment (only if it starts with { and is a whole line)
        if (trimmed.size() >= 2 && trimmed.front() == '{' && trimmed.back() == '}') {
            // Make sure it's just a comment, not embedded in an equation
            if (trimmed.find('=') == std::string::npos) {
                std::string content = trimmed.substr(1, trimmed.size() - 2);
                return makeComment(content, true, lineNum);
            }
        }
        
        // Check for C-style comment
        if (trimmed.size() >= 2 && trimmed[0] == '/' && trimmed[1] == '/') {
            std::string content = trimmed.substr(2);
            return makeComment(content, false, lineNum);
        }
        
        return nullptr;
    }
    
    StmtPtr tryParseDirective(const std::string& line, int lineNum) {
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
        
        return makeDirective(name, content, lineNum);
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
        return true;
    }

    StmtPtr tryParseProcedureCall(const std::string& line, int lineNum) {
        std::string cleaned = removeInlineComments(line);
        cleaned = trim(cleaned);
        
        std::string upper = cleaned;
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        
        if (upper.find("CALL ") != 0) return nullptr;
        
        std::string callContent = trim(cleaned.substr(5));
        size_t parenStart = callContent.find('(');
        if (parenStart == std::string::npos) return nullptr;
        
        std::string name = trim(callContent.substr(0, parenStart));
        
        size_t parenEnd = callContent.rfind(')');
        if (parenEnd == std::string::npos || parenEnd <= parenStart) return nullptr;
        
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
        
        std::vector<ExprPtr> inputArgs;
        std::vector<std::string> inputStrs = splitArguments(inputsStr);
        for (const auto& s : inputStrs) {
            auto expr = parseExpression(s, lineNum);
            if (expr) inputArgs.push_back(expr);
        }
        
        std::vector<Variable> outputVars;
        std::vector<std::string> outputStrs = splitArguments(outputsStr);
        for (const auto& s : outputStrs) {
            std::string ts = trim(s);
            if (ts.empty()) continue;
            
            // Outputs must be variables (potentially array access)
            auto expr = parseExpression(ts, lineNum);
            if (expr && expr->is<Variable>()) {
                outputVars.push_back(expr->as<Variable>());
            }
        }
        
        return makeProcedureCall(name, std::move(inputArgs), std::move(outputVars), lineNum);
    }

    StmtPtr tryParseEquationOrAssignment(const std::string& line, int lineNum) {
        // Remove inline comments first
        std::string cleaned = removeInlineComments(line);
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
        
        auto lhs = parseExpression(lhsStr, lineNum);
        auto rhs = parseExpression(rhsStr, lineNum);
        
        // If RHS parsing failed, it might be due to unquoted units (e.g., 10 [m])
        // Try stripping units in square brackets at the end
        if (!rhs) {
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
        
        if (!lhs || !rhs) return nullptr;
        
        auto stmt = std::make_shared<Statement>();
        if (op == ":=") {
            stmt->node = Assignment{lhs, rhs};
        } else {
            stmt->node = Equation{lhs, rhs, units, ""};
        }
        stmt->sourceLineNumber = lineNum;
        return stmt;
    }
    
    std::string removeInlineComments(const std::string& line) {
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
                    op = "=";
                    return static_cast<int>(i);
                }
            }
        }
        return -1;
    }
    
    ExprPtr parseExpression(const std::string& expr, int lineNum) {
        std::string trimmed = trim(expr);
        if (trimmed.empty()) return nullptr;
        
        // Parse as additive expression
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
        
        // String literal
        if (trimmed.front() == '\'') {
            size_t end = trimmed.find('\'', 1);
            if (end != std::string::npos) {
                std::string value = trimmed.substr(1, end - 1);
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
        
        // Function call or array access
        size_t parenPos = trimmed.find('(');
        size_t bracketPos = trimmed.find('[');
        
        if (parenPos != std::string::npos && 
            (bracketPos == std::string::npos || parenPos < bracketPos)) {
            return parseFunctionCall(trimmed, lineNum);
        }
        
        if (bracketPos != std::string::npos) {
            return parseArrayAccess(trimmed, lineNum);
        }
        
        // Simple variable
        if (isIdentifier(trimmed)) {
            return makeVariable(trimmed, lineNum);
        }
        
        return nullptr;
    }
    
    ExprPtr parseFunctionCall(const std::string& expr, int lineNum) {
        size_t parenStart = expr.find('(');
        if (parenStart == std::string::npos) return nullptr;
        
        std::string name = trim(expr.substr(0, parenStart));
        
        // Find matching close paren
        size_t parenEnd = expr.size() - 1;
        while (parenEnd > parenStart && expr[parenEnd] != ')') parenEnd--;
        
        if (expr[parenEnd] != ')') return nullptr;
        
        std::string argsStr = expr.substr(parenStart + 1, parenEnd - parenStart - 1);
        
        auto func = std::make_shared<Expression>();
        FunctionCall fc;
        fc.name = name;
        
        // Check if this is a thermophysical function (first arg is fluid name)
        bool isThermo = isEESThermophysicalFunction(name);
        bool isFirstPositionalArg = true;
        
        // Parse arguments
        std::vector<std::string> argStrings = splitArguments(argsStr);
        
        for (const auto& argStr : argStrings) {
            std::string arg = trim(argStr);
            
            // Check for named argument (e.g., T=value)
            size_t eqPos = arg.find('=');
            if (eqPos != std::string::npos && eqPos > 0) {
                // Check if the part before = is an identifier (not an expression)
                std::string paramName = trim(arg.substr(0, eqPos));
                std::string paramValue = trim(arg.substr(eqPos + 1));
                
                if (isIdentifier(paramName)) {
                    auto valueExpr = parseExpression(paramValue, lineNum);
                    if (valueExpr) {
                        fc.namedArgs.push_back({paramName, valueExpr});
                        continue;
                    }
                }
            }
            
            // Regular argument
            auto argExpr = parseExpression(arg, lineNum);
            if (argExpr) {
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
        
        size_t bracketEnd = expr.rfind(']');
        if (bracketEnd == std::string::npos || bracketEnd <= bracketStart) return nullptr;
        
        std::string indexStr = expr.substr(bracketStart + 1, bracketEnd - bracketStart - 1);
        
        std::vector<ExprPtr> indices;
        std::vector<std::string> indexStrings = splitArguments(indexStr);
        
        for (const auto& idx : indexStrings) {
            auto indexExpr = parseExpression(idx, lineNum);
            if (indexExpr) {
                indices.push_back(indexExpr);
            }
        }
        
        return makeArrayVariable(name, indices, lineNum);
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
        
        if (!current.empty()) {
            result.push_back(current);
        }
        
        return result;
    }
    
    bool isNumber(const std::string& s) {
        if (s.empty()) return false;
        
        size_t i = 0;
        if (s[i] == '-' || s[i] == '+') i++;
        
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
        if (s.empty()) return false;
        
        // First character must be a letter or underscore
        char first = s[0];
        if (!std::isalpha(static_cast<unsigned char>(first)) && first != '_') {
            return false;
        }
        
        // Check remaining characters
        for (size_t i = 1; i < s.size(); ++i) {
            char c = s[i];
            // $ is only allowed as the last character (for string variables)
            if (c == '$') {
                if (i != s.size() - 1) return false;
            } 
            // # is only allowed as the last character (for constants)
            else if (c == '#') {
                if (i != s.size() - 1) return false;
            }
            else if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') {
                return false;
            }
        }
        
        return true;
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
