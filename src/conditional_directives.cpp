#include "coolsolve/conditional_directives.h"

#include "coolsolve/ees_string.h"

#include <cctype>
#include <stdexcept>
#include <vector>

namespace coolsolve {

namespace {

struct ConditionError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

enum class TokenKind { Word, String, LParen, RParen, Operator };

struct Token {
    TokenKind kind;
    std::string text;   ///< For a string: the text between the quotes.
};

std::string upperCase(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

/// Words are runs of characters other than blanks, parentheses, quotes and the
/// characters of operators, so that `Tablerun#<13`, `X$='a'` and `MIN/MAX` split as expected.
bool isOperatorChar(char c) { return c == '=' || c == '<' || c == '>' || c == '!'; }

std::vector<Token> tokenize(const std::string& s) {
    std::vector<Token> tokens;
    size_t i = 0;
    while (i < s.size()) {
        const char c = s[i];
        if (std::isspace(static_cast<unsigned char>(c))) { ++i; continue; }
        if (c == '(') { tokens.push_back({TokenKind::LParen, "("}); ++i; continue; }
        if (c == ')') { tokens.push_back({TokenKind::RParen, ")"}); ++i; continue; }
        if (c == '\'') {
            const size_t close = s.find('\'', i + 1);
            if (close == std::string::npos)
                throw ConditionError("the string constant " + s.substr(i) + " has no closing quote");
            tokens.push_back({TokenKind::String, s.substr(i + 1, close - i - 1)});
            i = close + 1;
            continue;
        }
        const size_t start = i;
        if (isOperatorChar(c)) {
            while (i < s.size() && isOperatorChar(s[i])) ++i;
            tokens.push_back({TokenKind::Operator, s.substr(start, i - start)});
            continue;
        }
        while (i < s.size() && !std::isspace(static_cast<unsigned char>(s[i])) && s[i] != '(' && s[i] != ')' &&
               s[i] != '\'' && !isOperatorChar(s[i]))
            ++i;
        tokens.push_back({TokenKind::Word, s.substr(start, i - start)});
    }
    return tokens;
}

/// Recursive-descent evaluation of `or-expr := and-expr {OR and-expr}`,
/// `and-expr := term {AND term}`, `term := (or-expr) | test`.
class ConditionParser {
public:
    ConditionParser(std::vector<Token> tokens, const CompileTimeState& state)
        : tokens_(std::move(tokens)), state_(state) {}

    bool usedIntegralTable() const { return usedIntegralTable_; }

    bool parse() {
        if (tokens_.empty())
            throw ConditionError("the directive has no condition (e.g. $IF fluid$='Water' or $IFNOT PARAMETRICTABLE)");
        const bool value = parseOr();
        if (peekWord("THEN")) ++pos_;   // `$if x$='a' then` is accepted
        if (pos_ < tokens_.size())
            throw ConditionError("unexpected text '" + tokens_[pos_].text + "' after the condition");
        return value;
    }

private:
    std::vector<Token> tokens_;
    const CompileTimeState& state_;
    size_t pos_ = 0;
    bool usedIntegralTable_ = false;

    bool atEnd() const { return pos_ >= tokens_.size(); }
    bool peekWord(const char* upper) const {
        return !atEnd() && tokens_[pos_].kind == TokenKind::Word && upperCase(tokens_[pos_].text) == upper;
    }
    bool peekOperator() const { return !atEnd() && tokens_[pos_].kind == TokenKind::Operator; }

    bool parseOr() {
        bool value = parseAnd();
        while (peekWord("OR")) {
            ++pos_;
            const bool right = parseAnd();   // every test is checked: no short-circuit
            value = value || right;
        }
        return value;
    }

    bool parseAnd() {
        bool value = parseTerm();
        while (peekWord("AND")) {
            ++pos_;
            const bool right = parseTerm();
            value = value && right;
        }
        return value;
    }

    bool parseTerm() {
        if (atEnd()) throw ConditionError("the condition is incomplete: a test is missing");
        if (tokens_[pos_].kind == TokenKind::LParen) {
            ++pos_;
            const bool value = parseOr();
            if (atEnd() || tokens_[pos_].kind != TokenKind::RParen) throw ConditionError("missing ')' in the condition");
            ++pos_;
            return value;
        }
        if (tokens_[pos_].kind != TokenKind::Word)
            throw ConditionError("expected a test, found '" + tokens_[pos_].text + "'");
        return parseTest();
    }

    /// The text of a string constant: quoted or, as EES allows, a single unquoted word.
    std::string takeStringConstant(const std::string& what) {
        if (atEnd() || (tokens_[pos_].kind != TokenKind::String && tokens_[pos_].kind != TokenKind::Word))
            throw ConditionError("a string constant is expected after " + what);
        return tokens_[pos_++].text;
    }

    bool parseTest() {
        const std::string written = tokens_[pos_++].text;
        const std::string name = upperCase(written);

        if (name == "UNITSYSTEM") return parseUnitSystem();
        if (name == "VERSION")
            throw ConditionError("VERSION conditions are not supported: CoolSolve is not EES and has no EES version number");
        if (name == "TABLERUN#") return parseTableRun();
        if (!written.empty() && written.back() == '$') return parseStringTest(written);

        // Flags: false in a CoolSolve solve, but PROFESSIONAL and INTEGRALTABLE
        std::string flag = name;
        if (flag == "MIN/MAX") flag = "MINMAX";
        if (flag == "MIN/MAXTABLE") flag = "MINMAXTABLE";
        bool value = false;
        bool takesTableName = false;
        if (flag == "PARAMETRICTABLE" || flag == "MINMAXTABLE" || flag == "UNCERTAINTYTABLE") takesTableName = true;
        else if (flag == "PROFESSIONAL") value = true;
        else if (flag == "INTEGRALTABLE") { value = state_.integralTable; usedIntegralTable_ = true; }
        else if (flag != "MINMAX" && flag != "UNCERTAINTY" && flag != "DIAGRAMWINDOW" && flag != "MACRO")
            throw ConditionError("unknown condition '" + written + "': expected a string test (X$='text'), PARAMETRICTABLE, "
                "MINMAX, MINMAXTABLE, UNCERTAINTY, UNCERTAINTYTABLE, DIAGRAMWINDOW, INTEGRALTABLE, MACRO, PROFESSIONAL, "
                "TABLERUN#, UNITSYSTEM('..'), or a combination with AND / OR");

        if (peekOperator()) {
            if (!takesTableName || tokens_[pos_].text != "=")
                throw ConditionError("'" + written + "' cannot be followed by '" + tokens_[pos_].text + "'");
            ++pos_;
            if (atEnd() || tokens_[pos_].kind != TokenKind::String)
                throw ConditionError("the table name after '" + written + "=' must be a string constant in single quotes "
                    "(a string variable is not accepted)");
            ++pos_;
            value = false;   // CoolSolve has no named parametric / min-max / uncertainty table
        }
        return value;
    }

    /// `UNITSYSTEM('XX')` or `UNITSYSTEM('XX')=n`.
    bool parseUnitSystem() {
        if (atEnd() || tokens_[pos_].kind != TokenKind::LParen) throw ConditionError("UNITSYSTEM needs a setting: UNITSYSTEM('K')");
        ++pos_;
        const std::string setting = takeStringConstant("UNITSYSTEM(");
        if (atEnd() || tokens_[pos_].kind != TokenKind::RParen) throw ConditionError("missing ')' after UNITSYSTEM('" + setting + "'");
        ++pos_;
        bool active = false;
        try {
            active = unitSettingActive(setting, state_.units);
        } catch (const std::runtime_error& e) {
            throw ConditionError(e.what());
        }
        if (peekOperator()) {
            if (tokens_[pos_].text != "=") throw ConditionError("UNITSYSTEM('" + setting + "') can only be compared with '='");
            ++pos_;
            if (atEnd() || tokens_[pos_].kind != TokenKind::Word) throw ConditionError("a number is expected after UNITSYSTEM('" + setting + "')=");
            const std::string number = tokens_[pos_++].text;
            try {
                size_t used = 0;
                const double n = std::stod(number, &used);
                if (used != number.size()) throw std::invalid_argument(number);
                return (active ? 1.0 : 0.0) == n;
            } catch (const std::logic_error&) {
                throw ConditionError("a number is expected after UNITSYSTEM('" + setting + "')=, found '" + number + "'");
            }
        }
        return active;
    }

    /// `TABLERUN# = n`, `< n`, `> n`: true only for a run of a parametric table, never in a solve.
    bool parseTableRun() {
        if (!peekOperator() || (tokens_[pos_].text != "=" && tokens_[pos_].text != "<" && tokens_[pos_].text != ">"))
            throw ConditionError("TABLERUN# must be followed by =, < or > and an integer number");
        ++pos_;
        if (atEnd() || tokens_[pos_].kind != TokenKind::Word) throw ConditionError("an integer number is expected after TABLERUN#");
        const std::string number = tokens_[pos_++].text;
        const bool integer = !number.empty() && number.find_first_not_of("0123456789") == std::string::npos;
        if (!integer) throw ConditionError("the value compared with TABLERUN# must be an integer number, found '" + number + "'");
        return false;
    }

    /// `X$ = 'text'`: X$ must hold a string constant set before the directive.
    bool parseStringTest(const std::string& variable) {
        if (!peekOperator())
            throw ConditionError("a test on the string variable '" + variable + "' must be written " + variable + "='text'");
        const std::string op = tokens_[pos_].text;
        if (op != "=")
            throw ConditionError("'" + op + "' is not an EES condition on strings; only " + variable +
                                 "='text' is (use $IFNOT for the opposite)");
        ++pos_;
        const std::string constant = takeStringConstant(variable + "=");

        const auto it = state_.strings.find(variable);
        if (it != state_.strings.end()) return eesStringEqual(it->second, constant);
        if (state_.computedStrings.count(variable))
            throw ConditionError("the string variable '" + variable + "' is set by a CALL or an expression: $IF is "
                "evaluated when the equations are compiled, so its value must be a string constant ('" + variable + " = 'text'') "
                "set before the directive");
        throw ConditionError("the string variable '" + variable + "' has no value when the directive is compiled: it must be "
            "set to a string constant (" + variable + " = 'text') before the $IF (a value given in the EES Diagram window "
            "is not available in CoolSolve)");
    }
};

}  // namespace

ConditionResult evaluateDirectiveCondition(const std::string& condition, const CompileTimeState& state) {
    ConditionResult result;
    try {
        ConditionParser parser(tokenize(condition), state);
        result.value = parser.parse();
        result.usedIntegralTable = parser.usedIntegralTable();
        result.ok = true;
    } catch (const ConditionError& e) {
        result.error = e.what();
    }
    return result;
}

}  // namespace coolsolve
