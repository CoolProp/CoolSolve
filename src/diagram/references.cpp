// The badge lines of a chart (include/coolsolve/diagram/references.h).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/references.py
// (latex_to_text, parse_bibtex, cite, fluid_references, blend_name, composition_text, versions_line) and
// engine/engine.py (Engine.source_lines). Each function names the Python function it follows; the regular
// expressions of the Python code are written out as scans that match the same text.
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#include "coolsolve/diagram/references.h"

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"

#include <coolplot/capi.h>
#include <coolplot/core.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string_view>

namespace CoolProp {
// Backends/Helmholtz/MixtureParameters.h of the CoolProp sources, which CoolProp's Python module exposes; declared here
// rather than included, since that header pulls in the internals of the Helmholtz backend
std::string get_mixture_binary_pair_data(const std::string& CAS1, const std::string& CAS2, const std::string& param);
}  // namespace CoolProp

namespace coolsolve {
namespace diagram {

// CoolPropBibTeXLibrary.bib of the CoolProp sources, written into the build by cmake/embed_bibtex.cmake
extern const unsigned char COOLPROP_BIBTEX[];
extern const std::size_t COOLPROP_BIBTEX_SIZE;

const char* const HUMID_AIR_KEY = "Herrmann-2009";

namespace {

// ---------------------------------------------------------------------------
// Characters
// ---------------------------------------------------------------------------

bool isSpace(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

bool isLetter(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

/// Python's \w on ASCII text.
bool isWordChar(char c) {
    return isLetter(c) || (c >= '0' && c <= '9') || c == '_';
}

std::string lower(std::string text) {
    for (char& c : text)
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    return text;
}

std::string strip(const std::string& text) {
    size_t begin = 0, end = text.size();
    while (begin < end && isSpace(text[begin])) ++begin;
    while (end > begin && isSpace(text[end - 1])) --end;
    return text.substr(begin, end - begin);
}

bool endsWith(const std::string& text, const std::string& tail) {
    return text.size() >= tail.size() && text.compare(text.size() - tail.size(), tail.size(), tail) == 0;
}

std::string utf8(char32_t code) {
    std::string out;
    if (code < 0x80) {
        out += static_cast<char>(code);
    } else if (code < 0x800) {
        out += static_cast<char>(0xC0 | (code >> 6));
        out += static_cast<char>(0x80 | (code & 0x3F));
    } else {
        out += static_cast<char>(0xE0 | (code >> 12));
        out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
        out += static_cast<char>(0x80 | (code & 0x3F));
    }
    return out;
}

/// The first code point of UTF-8 text (0 for empty text).
char32_t firstCodePoint(const std::string& text) {
    if (text.empty()) return 0;
    const unsigned char c = static_cast<unsigned char>(text[0]);
    if (c < 0x80 || text.size() < 2) return c;
    if (c < 0xE0) return ((c & 0x1F) << 6) | (static_cast<unsigned char>(text[1]) & 0x3F);
    if (text.size() < 3) return c;
    return ((c & 0x0F) << 12) | ((static_cast<unsigned char>(text[1]) & 0x3F) << 6) | (static_cast<unsigned char>(text[2]) & 0x3F);
}

// ---------------------------------------------------------------------------
// LaTeX of a BibTeX field (latex_to_text)
// ---------------------------------------------------------------------------

/// The combining accent of each LaTeX accent command (references._ACCENTS).
char32_t accentMark(char command) {
    switch (command) {
        case '\'': return 0x0301;
        case '`': return 0x0300;
        case '"': return 0x0308;
        case '^': return 0x0302;
        case '~': return 0x0303;
        case 'v': return 0x030C;
        case 'c': return 0x0327;
        case 'r': return 0x030A;
        case '=': return 0x0304;
        case '.': return 0x0307;
        case 'u': return 0x0306;
        case 'H': return 0x030B;
        default: return 0;
    }
}

/// The letters of Unicode made of an ASCII letter and a combining accent: what unicodedata.normalize("NFC")
/// composes (generated with Python's unicodedata for every ASCII letter and the accents above).
struct CompositionRow {
    char32_t mark;
    const char* letters;
    std::vector<char32_t> codes;  ///< in the order of the letters
};

const std::map<std::pair<char32_t, char>, char32_t>& compositions() {
    static const std::map<std::pair<char32_t, char>, char32_t> table = [] {
        const std::vector<CompositionRow> rows = {
            {0x0301, "acegiklmnoprsuwyzACEGIKLMNOPRSUWYZ",
             {0xE1, 0x107, 0xE9, 0x1F5, 0xED, 0x1E31, 0x13A, 0x1E3F, 0x144, 0xF3, 0x1E55, 0x155, 0x15B, 0xFA, 0x1E83, 0xFD, 0x17A,
              0xC1, 0x106, 0xC9, 0x1F4, 0xCD, 0x1E30, 0x139, 0x1E3E, 0x143, 0xD3, 0x1E54, 0x154, 0x15A, 0xDA, 0x1E82, 0xDD, 0x179}},
            {0x0300, "aeinouwyAEINOUWY",
             {0xE0, 0xE8, 0xEC, 0x1F9, 0xF2, 0xF9, 0x1E81, 0x1EF3, 0xC0, 0xC8, 0xCC, 0x1F8, 0xD2, 0xD9, 0x1E80, 0x1EF2}},
            {0x0308, "aehiotuwxyAEHIOUWXY",
             {0xE4, 0xEB, 0x1E27, 0xEF, 0xF6, 0x1E97, 0xFC, 0x1E85, 0x1E8D, 0xFF, 0xC4, 0xCB, 0x1E26, 0xCF, 0xD6, 0xDC, 0x1E84,
              0x1E8C, 0x178}},
            {0x0302, "aceghijosuwyzACEGHIJOSUWYZ",
             {0xE2, 0x109, 0xEA, 0x11D, 0x125, 0xEE, 0x135, 0xF4, 0x15D, 0xFB, 0x175, 0x177, 0x1E91, 0xC2, 0x108, 0xCA, 0x11C,
              0x124, 0xCE, 0x134, 0xD4, 0x15C, 0xDB, 0x174, 0x176, 0x1E90}},
            {0x0303, "aeinouvyAEINOUVY",
             {0xE3, 0x1EBD, 0x129, 0xF1, 0xF5, 0x169, 0x1E7D, 0x1EF9, 0xC3, 0x1EBC, 0x128, 0xD1, 0xD5, 0x168, 0x1E7C, 0x1EF8}},
            {0x030C, "acdeghijklnorstuzACDEGHIKLNORSTUZ",
             {0x1CE, 0x10D, 0x10F, 0x11B, 0x1E7, 0x21F, 0x1D0, 0x1F0, 0x1E9, 0x13E, 0x148, 0x1D2, 0x159, 0x161, 0x165, 0x1D4, 0x17E,
              0x1CD, 0x10C, 0x10E, 0x11A, 0x1E6, 0x21E, 0x1CF, 0x1E8, 0x13D, 0x147, 0x1D1, 0x158, 0x160, 0x164, 0x1D3, 0x17D}},
            {0x0327, "cdeghklnrstCDEGHKLNRST",
             {0xE7, 0x1E11, 0x229, 0x123, 0x1E29, 0x137, 0x13C, 0x146, 0x157, 0x15F, 0x163, 0xC7, 0x1E10, 0x228, 0x122, 0x1E28,
              0x136, 0x13B, 0x145, 0x156, 0x15E, 0x162}},
            {0x030A, "auwyAU", {0xE5, 0x16F, 0x1E98, 0x1E99, 0xC5, 0x16E}},
            {0x0304, "aegiouyAEGIOUY", {0x101, 0x113, 0x1E21, 0x12B, 0x14D, 0x16B, 0x233, 0x100, 0x112, 0x1E20, 0x12A, 0x14C, 0x16A, 0x232}},
            {0x0307, "abcdefghmnoprstwxyzABCDEFGHIMNOPRSTWXYZ",
             {0x227, 0x1E03, 0x10B, 0x1E0B, 0x117, 0x1E1F, 0x121, 0x1E23, 0x1E41, 0x1E45, 0x22F, 0x1E57, 0x1E59, 0x1E61, 0x1E6B,
              0x1E87, 0x1E8B, 0x1E8F, 0x17C, 0x226, 0x1E02, 0x10A, 0x1E0A, 0x116, 0x1E1E, 0x120, 0x1E22, 0x130, 0x1E40, 0x1E44,
              0x22E, 0x1E56, 0x1E58, 0x1E60, 0x1E6A, 0x1E86, 0x1E8A, 0x1E8E, 0x17B}},
            {0x0306, "aegiouAEGIOU", {0x103, 0x115, 0x11F, 0x12D, 0x14F, 0x16D, 0x102, 0x114, 0x11E, 0x12C, 0x14E, 0x16C}},
            {0x030B, "ouOU", {0x151, 0x171, 0x150, 0x170}},
        };
        std::map<std::pair<char32_t, char>, char32_t> out;
        for (const CompositionRow& row : rows)
            for (size_t i = 0; row.letters[i]; ++i) out[{row.mark, row.letters[i]}] = row.codes.at(i);
        return out;
    }();
    return table;
}

/// unicodedata.normalize("NFC", letter + mark): the composed letter, else the letter and the combining accent.
std::string composeAccent(char letter, char32_t mark) {
    const auto& table = compositions();
    const auto found = table.find({mark, letter});
    if (found != table.end()) return utf8(found->second);
    return std::string(1, letter) + utf8(mark);
}

/// Letters written as commands (references._LETTERS), in the order of the alternatives of latex_to_text.
const std::vector<std::pair<std::string, char32_t>>& letterCommands() {
    static const std::vector<std::pair<std::string, char32_t>> commands = {
        {"ss", 0xDF}, {"AA", 0xC5}, {"aa", 0xE5}, {"AE", 0xC6}, {"ae", 0xE6}, {"o", 0xF8},
        {"O", 0xD8},  {"l", 0x142}, {"L", 0x141}, {"i", 'i'}};
    return commands;
}

/// Python's str.islower() of the first character of a text latexToText wrote: an ASCII lower-case letter, a letter
/// composed from one, or a lower-case letter command (ß, ø, å, æ, ł).
bool startsLower(const std::string& text) {
    static const std::set<char32_t> lowerLetters = [] {
        std::set<char32_t> out = {0xDF, 0xF8, 0xE5, 0xE6, 0x142};
        for (const auto& [key, code] : compositions())
            if (key.second >= 'a' && key.second <= 'z') out.insert(code);
        return out;
    }();
    const char32_t first = firstCodePoint(text);
    return (first >= 'a' && first <= 'z') || lowerLetters.count(first) > 0;
}

}  // namespace

std::string latexToText(const std::string& input) {
    // {\"u}, \"{u}, \"u: \\(['`"^~=.])\s*(?:\{([A-Za-z])\}|([A-Za-z]))
    std::string text;
    for (size_t i = 0; i < input.size();) {
        if (input[i] == '\\' && i + 1 < input.size() && std::string_view("'`\"^~=.").find(input[i + 1]) != std::string_view::npos) {
            size_t j = i + 2;
            while (j < input.size() && isSpace(input[j])) ++j;
            if (j + 2 < input.size() && input[j] == '{' && isLetter(input[j + 1]) && input[j + 2] == '}') {
                text += composeAccent(input[j + 1], accentMark(input[i + 1]));
                i = j + 3;
                continue;
            }
            if (j < input.size() && isLetter(input[j])) {
                text += composeAccent(input[j], accentMark(input[i + 1]));
                i = j + 1;
                continue;
            }
        }
        text += input[i++];
    }
    // {\v s}, \r{A}: \\([vcruH])\s*(?:\{([A-Za-z])\}|\s([A-Za-z])) (a letter needs a space before it)
    std::string out;
    for (size_t i = 0; i < text.size();) {
        if (text[i] == '\\' && i + 1 < text.size() && std::string_view("vcruH").find(text[i + 1]) != std::string_view::npos) {
            size_t j = i + 2;
            while (j < text.size() && isSpace(text[j])) ++j;
            if (j + 2 < text.size() && text[j] == '{' && isLetter(text[j + 1]) && text[j + 2] == '}') {
                out += composeAccent(text[j + 1], accentMark(text[i + 1]));
                i = j + 3;
                continue;
            }
            if (j > i + 2 && j < text.size() && isLetter(text[j])) {
                out += composeAccent(text[j], accentMark(text[i + 1]));
                i = j + 1;
                continue;
            }
        }
        out += text[i++];
    }
    // \ss, \o, \aa, ...: \\(ss|AA|aa|AE|ae|o|O|l|L|i)(?![A-Za-z])
    text.clear();
    for (size_t i = 0; i < out.size();) {
        bool matched = false;
        if (out[i] == '\\') {
            for (const auto& [command, code] : letterCommands()) {
                const size_t end = i + 1 + command.size();
                if (out.compare(i + 1, command.size(), command) == 0 && (end >= out.size() || !isLetter(out[end]))) {
                    text += utf8(code);
                    i = end;
                    matched = true;
                    break;
                }
            }
        }
        if (!matched) text += out[i++];
    }
    // Inline math such as $_3$ or $^{13}$: the characters, without the markup
    out.clear();
    for (size_t i = 0; i < text.size();) {
        const size_t close = text[i] == '$' ? text.find('$', i + 1) : std::string::npos;
        if (close == std::string::npos) {
            out += text[i++];
            continue;
        }
        for (size_t k = i + 1; k < close; ++k)
            if (std::string_view("_^{}\\").find(text[k]) == std::string_view::npos) out += text[k];
        i = close + 1;
    }
    // Escaped special characters (\& -> &), then any other command (\\[A-Za-z]+\s*)
    text.clear();
    for (size_t i = 0; i < out.size();) {
        if (out[i] == '\\' && i + 1 < out.size() && std::string_view("&%$#_").find(out[i + 1]) != std::string_view::npos) {
            text += out[i + 1];
            i += 2;
        } else {
            text += out[i++];
        }
    }
    out.clear();
    for (size_t i = 0; i < text.size();) {
        if (text[i] == '\\' && i + 1 < text.size() && isLetter(text[i + 1])) {
            size_t j = i + 1;
            while (j < text.size() && isLetter(text[j])) ++j;
            while (j < text.size() && isSpace(text[j])) ++j;
            i = j;
        } else {
            out += text[i++];
        }
    }
    // Braces removed, "--" as "-", white space collapsed
    text.clear();
    for (char c : out)
        if (c != '{' && c != '}') text += c;
    out.clear();
    for (size_t i = 0; i < text.size();) {
        if (text.compare(i, 2, "--") == 0) {
            out += '-';
            i += 2;
        } else {
            out += text[i++];
        }
    }
    text.clear();
    for (size_t i = 0; i < out.size();) {
        if (isSpace(out[i])) {
            while (i < out.size() && isSpace(out[i])) ++i;
            text += ' ';
        } else {
            text += out[i++];
        }
    }
    return strip(text);
}

namespace {

// ---------------------------------------------------------------------------
// The bibliography (parse_bibtex)
// ---------------------------------------------------------------------------

/// key -> {"type": entry type, field name (lower case): raw value}
using Entry = std::map<std::string, std::string>;

/// (value, position after it) of one field value starting at `pos` (references._bibtex_value).
std::pair<std::string, size_t> bibtexValue(std::string_view text, size_t pos) {
    if (pos < text.size() && text[pos] == '{') {
        int depth = 0;
        const size_t start = pos;
        for (; pos < text.size(); ++pos) {
            if (text[pos] == '{') {
                ++depth;
            } else if (text[pos] == '}') {
                if (--depth == 0) return {std::string(text.substr(start + 1, pos - start - 1)), pos + 1};
            }
        }
        return {std::string(text.substr(start + 1)), text.size()};
    }
    if (pos < text.size() && text[pos] == '"') {
        size_t end = text.find('"', pos + 1);
        if (end == std::string_view::npos) end = text.size();
        return {std::string(text.substr(pos + 1, end - pos - 1)), end + 1};
    }
    size_t end = pos;
    while (end < text.size() && text[end] != ',' && text[end] != '}' && !isSpace(text[end])) ++end;
    return {std::string(text.substr(pos, end - pos)), end};
}

/// Entries of a BibTeX file: fields in braces (nested braces allowed) or quotes, or bare (references.parse_bibtex).
std::map<std::string, Entry> parseBibtex(std::string_view text) {
    std::map<std::string, Entry> entries;
    const size_t n = text.size();
    auto skipSpace = [&](size_t i) {
        while (i < n && isSpace(text[i])) ++i;
        return i;
    };
    // The entries start at @(\w+)\s*\{\s*([^,\s]+)\s*, and the search goes on after each start, as re.finditer does
    for (size_t at = text.find('@'); at != std::string_view::npos; at = text.find('@', at)) {
        size_t i = at + 1;
        while (i < n && isWordChar(text[i])) ++i;
        const std::string kind = lower(std::string(text.substr(at + 1, i - at - 1)));
        i = skipSpace(i);
        if (kind.empty() || i >= n || text[i] != '{') {
            ++at;
            continue;
        }
        i = skipSpace(i + 1);
        const size_t keyStart = i;
        while (i < n && text[i] != ',' && !isSpace(text[i])) ++i;
        const std::string key(text.substr(keyStart, i - keyStart));
        i = skipSpace(i);
        if (key.empty() || i >= n || text[i] != ',') {
            ++at;
            continue;
        }
        at = i + 1;
        if (kind == "comment" || kind == "string" || kind == "preamble") continue;
        Entry fields = {{"type", kind}};
        for (size_t pos = at;;) {
            // \s*(\w+)\s*=\s*
            size_t p = skipSpace(pos);
            const size_t nameStart = p;
            while (p < n && isWordChar(text[p])) ++p;
            if (p == nameStart) break;
            const std::string name = lower(std::string(text.substr(nameStart, p - nameStart)));
            p = skipSpace(p);
            if (p >= n || text[p] != '=') break;
            auto [value, next] = bibtexValue(text, skipSpace(p + 1));
            fields[name] = value;
            p = skipSpace(next);
            if (p >= n || text[p] != ',') break;
            pos = p + 1;
        }
        entries[key] = fields;
    }
    return entries;
}

const std::map<std::string, Entry>& bibliography() {
    static const std::map<std::string, Entry> entries =
        parseBibtex(std::string_view(reinterpret_cast<const char*>(COOLPROP_BIBTEX), COOLPROP_BIBTEX_SIZE));
    return entries;
}

// ---------------------------------------------------------------------------
// Authors (_split_top, _last_name, _authors_short)
// ---------------------------------------------------------------------------

/// Where a separator found at `i` ends, or npos.
using Separator = size_t (*)(const std::string&, size_t);

size_t commaAt(const std::string& text, size_t i) {
    return text[i] == ',' ? i + 1 : std::string::npos;
}

size_t spacesAt(const std::string& text, size_t i) {
    if (!isSpace(text[i])) return std::string::npos;
    while (i < text.size() && isSpace(text[i])) ++i;
    return i;
}

/// \s+and\s+
size_t andAt(const std::string& text, size_t i) {
    size_t j = spacesAt(text, i);
    if (j == std::string::npos || text.compare(j, 3, "and") != 0 || j + 3 >= text.size()) return std::string::npos;
    return spacesAt(text, j + 3);
}

/// `text` split at `separator` outside braces, the parts stripped, empty ones dropped (references._split_top).
std::vector<std::string> splitTop(const std::string& text, Separator separator) {
    std::vector<std::string> parts;
    int depth = 0;
    size_t start = 0;
    for (size_t i = 0; i < text.size();) {
        if (text[i] == '{') {
            ++depth;
            ++i;
        } else if (text[i] == '}') {
            depth = std::max(depth - 1, 0);
            ++i;
        } else if (const size_t end = separator(text, i); end != std::string::npos) {
            if (depth == 0) {
                parts.push_back(text.substr(start, i - start));
                start = end;
            }
            i = end;
        } else {
            ++i;
        }
    }
    parts.push_back(text.substr(start));
    std::vector<std::string> out;
    for (const std::string& part : parts)
        if (!strip(part).empty()) out.push_back(strip(part));
    return out;
}

/// Name suffixes that are not a family name: "Rachford Jr." is cited as "Rachford" (references._SUFFIXES).
bool isSuffix(const std::string& word) {
    static const std::set<std::string> suffixes = {"jr", "jr.", "sr", "sr.", "ii", "iii", "iv"};
    return suffixes.count(lower(word)) > 0;
}

std::string join(const std::vector<std::string>& parts, const std::string& separator) {
    std::string out;
    for (size_t i = 0; i < parts.size(); ++i) out += (i ? separator : "") + parts[i];
    return out;
}

/// The family name of one BibTeX author, with its particles ("de Vries") (references._last_name).
std::string lastName(const std::string& author) {
    const std::vector<std::string> parts = splitTop(author, commaAt);
    std::vector<std::string> words;
    if (parts.size() > 1) {
        words = splitTop(parts[0], spacesAt);
    } else {
        words = splitTop(author, spacesAt);
        while (words.size() > 1 && isSuffix(latexToText(words.back()))) words.pop_back();
        // From the first lower-case word after the first one ("von", "de", "van der"), else the last word
        size_t start = words.size();
        for (size_t i = 1; i + 1 < words.size(); ++i)
            if (startsLower(latexToText(words[i]))) {
                start = i;
                break;
            }
        if (start < words.size()) words.erase(words.begin(), words.begin() + static_cast<std::ptrdiff_t>(start));
        else if (!words.empty()) words.erase(words.begin(), words.end() - 1);
    }
    // text.split(): the words, without empty ones
    std::vector<std::string> tokens = splitTop(latexToText(join(words, " ")), spacesAt);
    while (tokens.size() > 1 && isSuffix(tokens.back())) tokens.pop_back();
    return join(tokens, " ");
}

std::string authorsShort(const std::string& authors) {
    std::vector<std::string> names;
    for (const std::string& author : splitTop(authors, andAt)) names.push_back(lastName(author));
    if (names.empty()) return "";
    if (names.size() == 1) return names[0];
    if (names.size() == 2) return names[0] + " and " + names[1];
    return names[0] + " et al.";
}

}  // namespace

bool hasBibliography() {
    return !bibliography().empty();
}

std::string cite(const std::string& key, const std::string& style) {
    const auto found = bibliography().find(key);
    if (found == bibliography().end() || !found->second.count("author")) return key;
    const Entry& entry = found->second;
    auto field = [&entry](const std::string& name) {
        const auto value = entry.find(name);
        return value == entry.end() ? std::string() : value->second;
    };
    const std::string authors = authorsShort(field("author"));
    const std::string year = latexToText(field("year"));
    if (style == "short") return year.empty() ? authors : authors + " (" + year + ")";
    if (field("type") == "article" && entry.count("journal")) {
        std::string source = latexToText(field("journal"));
        if (!field("volume").empty()) source += " " + latexToText(field("volume"));
        if (!year.empty()) source += " (" + year + ")";
        if (!field("pages").empty()) source += " " + latexToText(field("pages"));
        return authors + ", " + source;
    }
    // Reports, theses, books: the title up to its first colon names them best
    const std::string title = latexToText(field("title"));
    std::string source = strip(title.substr(0, title.find(':')));
    if (source.empty()) source = latexToText(entry.count("institution") ? field("institution") : field("publisher"));
    return source.empty() ? authors + " (" + year + ")" : authors + ", " + source + " (" + year + ")";
}

namespace {

// ---------------------------------------------------------------------------
// Fluids (fluid_references, _fluid_lines, mixtures)
// ---------------------------------------------------------------------------

/// The short name of a CoolProp backend, as in fluid strings: "HelmholtzEOSBackend" -> "HEOS"; a tabular backend
/// ("BICUBIC&HEOS") is named by the equation of state behind it (references._short_backend).
std::string shortBackend(const std::string& backend) {
    const std::string name = backend.substr(backend.rfind('&') == std::string::npos ? 0 : backend.rfind('&') + 1);
    for (const auto& [part, shortName] : std::vector<std::pair<std::string, std::string>>{
             {"Helmholtz", "HEOS"}, {"REFPROP", "REFPROP"}, {"PengRobinson", "PR"}, {"SRK", "SRK"}, {"Incompressible", "INCOMP"}})
        if (name.find(part) != std::string::npos) return shortName;
    return name;
}

/// The citations of several keys separated by "; ", one that reads as an earlier one left out
/// (references._join_citations).
std::string joinCitations(const std::vector<std::string>& keys, const std::string& style) {
    std::vector<std::string> texts;
    for (const std::string& key : keys) {
        const std::string text = cite(key, style);
        if (std::find(texts.begin(), texts.end(), text) == texts.end()) texts.push_back(text);
    }
    return join(texts, "; ");
}

/// CoolProp's own name of a fluid, e.g. "Nitrogen" for the alias "N2"; the name itself if unknown.
std::string canonicalName(const std::string& name) {
    try {
        return CoolProp::get_fluid_param_string(name, "name");
    } catch (const std::exception&) {
        return name;
    }
}

/// Mass fractions from mole fractions, w_i = x_i M_i / sum(x_j M_j); none if CoolProp lacks a molar mass.
std::optional<std::vector<double>> massFractions(const std::vector<std::string>& components, const std::vector<double>& fractions) {
    std::vector<double> weights;
    try {
        for (size_t i = 0; i < std::min(components.size(), fractions.size()); ++i) {
            std::unique_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory("HEOS", components[i]));
            weights.push_back(fractions[i] * state->molar_mass());
        }
    } catch (const std::exception&) {
        return std::nullopt;
    }
    double total = 0.0;
    for (double w : weights) total += w;
    if (!(total > 0.0)) return std::nullopt;
    for (double& w : weights) w /= total;
    return weights;
}

/// A predefined blend: its name without ".mix" and the mass fraction of each component.
using Blend = std::pair<std::string, std::map<std::string, double>>;

/// The predefined blends of CoolProp by their set of components (references._predefined_blends); a blend whose
/// components CoolProp cannot load is left out.
const std::map<std::set<std::string>, std::vector<Blend>>& predefinedBlends() {
    static const std::map<std::set<std::string>, std::vector<Blend>> blends = [] {
        std::map<std::set<std::string>, std::vector<Blend>> out;
        const std::string list = CoolProp::get_global_param_string("predefined_mixtures");
        for (size_t start = 0; start <= list.size();) {
            const size_t comma = std::min(list.find(',', start), list.size());
            const std::string name = list.substr(start, comma - start);
            start = comma + 1;
            if (!endsWith(name, ".mix")) continue;
            try {
                std::unique_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory("HEOS", name));
                const std::vector<std::string> components = state->fluid_names();
                const auto mass = state->get_mass_fractions();
                std::map<std::string, double> byComponent;
                for (size_t i = 0; i < components.size() && i < mass.size(); ++i) byComponent[components[i]] = static_cast<double>(mass[i]);
                out[std::set<std::string>(components.begin(), components.end())].push_back({name.substr(0, name.size() - 4), byComponent});
            } catch (const std::exception&) {
            }
        }
        return out;
    }();
    return blends;
}

/// "R32 70 %, R125 30 % by mole": every component with its fraction in per cent (references.composition_text).
std::string compositionText(const std::vector<std::string>& components, const std::vector<double>& fractions, const std::string& basis) {
    std::vector<std::string> parts;
    for (size_t i = 0; i < std::min(components.size(), fractions.size()); ++i)
        parts.push_back(components[i] + " " + coolplot::format_number(100.0 * fractions[i], 3) + " %");
    static const std::map<std::string, std::string> words = {{"mole", "by mole"}, {"mass", "by mass"}, {"volu", "by volume"}};
    const auto word = words.find(basis);
    return join(parts, ", ") + " " + (word == words.end() ? "by " + basis : word->second);
}

/// "R407C: R32 23 %, R125 25 %, R134a 52 % by mass", or "Mixture: ..." without a blend of that composition; a named
/// refrigerant blend by mass, as its R number is defined (references._composition_line).
std::string compositionLine(const std::string& backend, const std::vector<std::string>& components, std::vector<double> fractions,
                            std::string basis, std::string blend) {
    std::optional<std::vector<double>> mass;
    if (backend == "HEOS") {
        if (basis == "mass") mass = fractions;
        else if (basis == "mole") mass = massFractions(components, fractions);
    }
    if (mass) {
        if (blend.empty()) blend = blendName(components, *mass);
        if (blend.size() >= 2 && blend[0] == 'R' && blend[1] >= '0' && blend[1] <= '9') {
            fractions = *mass;
            basis = "mass";
        }
    }
    return (blend.empty() ? "Mixture" : blend) + ": " + compositionText(components, fractions, basis);
}

/// The BibTeX key of the binary pair of two components, in either order, or "" (references._pair_key).
std::string pairKey(const std::string& a, const std::string& b) {
    std::string casA, casB;
    try {
        casA = CoolProp::get_fluid_param_string(a, "CAS");
        casB = CoolProp::get_fluid_param_string(b, "CAS");
    } catch (const std::exception&) {
        return "";
    }
    for (const auto& [first, second] : {std::make_pair(casA, casB), std::make_pair(casB, casA)}) {
        try {
            return CoolProp::get_mixture_binary_pair_data(first, second, "BibTeX");
        } catch (const std::exception&) {
            // CoolProp throws for an unknown pair, or for one given in the other order
        }
    }
    return "";
}

std::vector<std::string> fluidLines(const std::string& backend, const std::vector<std::string>& components,
                                    const std::vector<double>& fractions, const std::string& basis, const std::string& blend,
                                    const std::string& style) {
    if (backend != "HEOS") {
        std::string names = fractions.empty() ? join(components, ", ") : compositionText(components, fractions, basis);
        static const std::map<std::string, std::string> labels = {{"REFPROP", "REFPROP"},
                                                                  {"PR", "Peng-Robinson equation of state"},
                                                                  {"SRK", "Soave-Redlich-Kwong equation of state"},
                                                                  {"INCOMP", "CoolProp incompressible fluid"}};
        const auto label = labels.find(backend);
        return {names + ": " + (label == labels.end() ? backend + " backend of CoolProp" : label->second)};
    }
    std::vector<std::string> lines;
    if (components.size() > 1 && fractions.size() == components.size())
        lines.push_back(compositionLine(backend, components, fractions, basis, blend));
    for (const std::string& name : components) {
        std::string key;
        try {
            key = CoolProp::get_fluid_param_string(name, "BibTeX-EOS");
        } catch (const std::exception&) {
            continue;  // an unknown name: nothing to cite
        }
        // Some fluids name several sources, separated by commas
        std::vector<std::string> sources;
        for (size_t start = 0; start <= key.size();) {
            const size_t comma = std::min(key.find(',', start), key.size());
            const std::string source = strip(key.substr(start, comma - start));
            if (!source.empty()) sources.push_back(source);
            start = comma + 1;
        }
        lines.push_back(name + ": " + joinCitations(sources, style));
    }
    if (components.size() > 1) {
        std::vector<std::string> keys;
        for (size_t i = 0; i < components.size(); ++i)
            for (size_t j = i + 1; j < components.size(); ++j) {
                const std::string pair = pairKey(components[i], components[j]);
                if (!pair.empty() && std::find(keys.begin(), keys.end(), pair) == keys.end()) keys.push_back(pair);
            }
        if (!keys.empty()) lines.push_back("Mixture model: " + joinCitations(keys, style));
    }
    return lines;
}

/// Truth of a JSON value as Python's bool() reads the value of a record.
bool truthy(const json& value) {
    if (value.is_boolean()) return value.get<bool>();
    if (value.is_number()) return value.get<double>() != 0.0;
    if (value.is_string() || value.is_array() || value.is_object()) return !value.empty();
    return false;
}

}  // namespace

std::string blendName(const std::vector<std::string>& components, const std::vector<double>& massFractions) {
    std::vector<std::string> names;
    for (const std::string& component : components) names.push_back(canonicalName(component));
    const auto blends = predefinedBlends().find(std::set<std::string>(names.begin(), names.end()));
    if (blends == predefinedBlends().end()) return "";
    std::string best;
    double bestDistance = BLEND_TOLERANCE;
    for (const auto& [name, blend] : blends->second) {
        double distance = 0.0;
        for (size_t i = 0; i < std::min(names.size(), massFractions.size()); ++i)
            distance = std::max(distance, std::fabs(blend.at(names[i]) - massFractions[i]));
        if (distance <= bestDistance) {
            best = name;
            bestDistance = distance;
        }
    }
    return best;
}

std::vector<std::string> fluidReferences(const json& fluid, const std::string& style) {
    if (style != "short" && style != "full")
        throw DiagramError("Unknown reference style '" + style + "', expected 'short' or 'full'.");
    if (!fluid.is_object()) return {};
    if (fluid.value("model", std::string()) == "humid_air") return {"Humid air: " + cite(HUMID_AIR_KEY, style)};

    // The lines of a fluid are kept, as CoolPlot's lru_cache keeps them: a blend asks CoolProp for every predefined blend
    static std::mutex mutex;
    static std::map<std::string, std::vector<std::string>> kept;
    const std::string cacheKey = fluid.dump() + "|" + style;
    {
        std::lock_guard<std::mutex> lock(mutex);
        const auto found = kept.find(cacheKey);
        if (found != kept.end()) return found->second;
    }
    std::string backend = shortBackend(fluid.value("backend", std::string("HEOS")));
    std::vector<std::string> components;
    std::vector<double> fractions;
    try {
        components = fluid.value("components", std::vector<std::string>());
        fractions = fluid.value("fractions", std::vector<double>());
    } catch (const json::exception&) {
        throw DiagramError("The fluid record has no list of components or fractions: " + fluid.dump() + ".");
    }
    std::string basis = fluid.value("basis", std::string("mole"));
    std::string blend;
    if (backend == "HEOS" && components.size() == 1 && endsWith(lower(components[0]), ".mix")) {
        // A predefined blend such as "R407C.mix": described by its components (references._blend_composition)
        try {
            std::unique_ptr<CoolProp::AbstractState> state(CoolProp::AbstractState::factory("HEOS", components[0]));
            const std::vector<std::string> names = state->fluid_names();
            const auto moles = state->get_mole_fractions();
            blend = components[0].substr(0, components[0].size() - 4);
            components = names;
            fractions.assign(moles.begin(), moles.end());
            basis = "mole";
        } catch (const std::exception&) {
        }
    }
    std::vector<std::string> lines = fluidLines(backend, components, fractions, basis, blend, style);
    std::lock_guard<std::mutex> lock(mutex);
    if (kept.size() >= 256) kept.clear();
    kept[cacheKey] = lines;
    return lines;
}

std::string versionsLine() {
    std::string version = CoolProp::get_global_param_string("version");
    if (version.find("dev") != std::string::npos) {
        const std::string revision = CoolProp::get_global_param_string("gitrevision").substr(0, 7);
        if (!revision.empty()) version += " (" + revision + ")";
    }
    return "CoolProp " + version + ", CoolPlot " + coolplot_version();
}

std::vector<std::string> sourceLines(const json& spec) {
    const json settings = spec.is_object() ? spec.value("settings", json::object()) : json::object();
    json branding = settings.is_object() ? settings.value("branding", json()) : json();
    // A record before 1.2 always holds its other parts, never the branding: it keeps the chart without badge it was
    // drawn as (Settings.from_dict)
    if (branding.is_null() && settings.is_object() && !settings.value("calc", json()).is_null())
        branding = {{"references", "none"}, {"versions", false}};
    if (!branding.is_object()) branding = json::object();
    const json references = branding.value("references", json("short"));
    if (!references.is_string() || (references != "none" && references != "short" && references != "full"))
        throw DiagramError("Unknown references style '" + (references.is_string() ? references.get<std::string>() : references.dump()) +
                           "', expected one of ('none', 'short', 'full').");
    std::vector<std::string> lines;
    if (references != "none") lines = fluidReferences(spec.value("fluid", json()), references.get<std::string>());
    if (truthy(branding.value("versions", json(true)))) lines.push_back(versionsLine());
    return lines;
}

}  // namespace diagram
}  // namespace coolsolve
