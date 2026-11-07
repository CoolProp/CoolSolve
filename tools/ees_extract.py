#!/usr/bin/env python3
"""
ees_extract.py -- extract the content of a binary EES file (.ees) for CoolSolve.

EES (Engineering Equation Solver) stores a model in a binary file written by a
Delphi program.  This tool recovers, without needing EES itself:

  * the Equations-window text (plain text or RTF), decoded from Windows-1252
    and normalised to UTF-8 with LF line endings;
  * the unit-system settings (mass/molar basis, C/K, kPa/bar/Pa, deg/rad, J/kJ);
  * the Variable Information records: stored solution value, guess value,
    lower/upper bounds and units of every variable of the main program;
  * the embedded Lookup and Parametric tables (column names, units, values).

The format knowledge was reverse-engineered on ~480 files written by EES 6.x
to 10.x (see docs/ees_import.md, section "EES binary format").  Equations are
read reliably for all versions; variable records and tables are decoded for
EES 7.x-10.x (EES 6.x files: equations only).

Typical use (creates a CoolSolve-ready folder):

    python3 tools/ees_extract.py model.ees -o out_dir [--name my_model]
    python3 tools/ees_extract.py --info model.ees        # summary only

Output folder content:

    <name>.eescode                    equations (EES syntax, UTF-8)
    <name>.initials                   guess values (stored EES solution)
    <name>-<table>.csv                one file per Lookup table (CoolSolve naming)
    reference/ees_variables.csv       EES stored solution: name,value,guess,lower,upper,units
                                      (`#` header lines: the stored values can be stale, see below)
    reference/ees_parametric-<t>.csv  Parametric tables
    reference/ees_extract_report.md   provenance, unit system, warnings

The stored solution is NOT guaranteed to be the solution of the equations: EES does not always
rewrite its Variable Info records when it saves.  The tool warns about the cases it can detect (inf/NaN
values, values that contradict an explicit equation, variables without record); cross-check the
reference before relying on it (docs/ees_import.md, "Stored solution may be stale").

Companion script: tools/compare_solution.py compares a CoolSolve .sol file with
reference/ees_variables.csv.
"""

import argparse
import ast
import collections
import csv
import datetime
import math
import os
import re
import struct
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Low-level decoding helpers
# ---------------------------------------------------------------------------

TRIG_DEG_FACTOR = math.pi / 180.0
EES_INF_EXP = 0x73E3          # EES stores +/-"infinity" as +/-1E+3999 (exp 0x73E3)
EES_NO_VALUE = -9999.0        # EES marks cleared / not-yet-computed values with -9999


def read_extended(buf, off):
    """Decode an 80-bit little-endian x87 'Extended' float (Delphi)."""
    mant = int.from_bytes(buf[off:off + 8], "little")
    se = int.from_bytes(buf[off + 8:off + 10], "little")
    sign = -1.0 if se & 0x8000 else 1.0
    exp = se & 0x7FFF
    if exp == 0 and mant == 0:
        return 0.0
    if exp == 0x7FFF:
        return float("nan")
    e = exp - 16383
    if e > 1023:
        return sign * math.inf
    if e < -1074:
        return 0.0
    return sign * math.ldexp(mant / 2.0 ** 63, e)


def plausible_extended(buf, off):
    """True if the 10 bytes at `off` look like a value EES would store."""
    if off + 10 > len(buf):
        return False
    mant = int.from_bytes(buf[off:off + 8], "little")
    exp = int.from_bytes(buf[off + 8:off + 10], "little") & 0x7FFF
    if exp == 0 and mant == 0:
        return True
    if not mant >> 63:                       # explicit integer bit must be set
        return False
    return abs(exp - 16383) < 400 or exp == EES_INF_EXP


def pascal_string(buf, off, maxlen=255):
    """Return (string, length) of a Delphi short string, or (None, 0)."""
    if off >= len(buf):
        return None, 0
    n = buf[off]
    if n > maxlen or off + 1 + n > len(buf):
        return None, 0
    return buf[off + 1:off + 1 + n].decode("cp1252", errors="replace"), n


def fmt_num(x):
    """Compact, round-trippable text for a float (empty for None/nan)."""
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return f"{x:.10g}"


# ---------------------------------------------------------------------------
# RTF -> plain text (EES stores 'formatted equations' as RTF)
# ---------------------------------------------------------------------------

_RTF_SKIP_DEST = {"fonttbl", "colortbl", "stylesheet", "info", "pict", "object",
                  "header", "footer", "listtable", "listoverridetable",
                  "rsidtbl", "generator", "themedata", "datastore", "latentstyles"}


def rtf_to_text(rtf):
    """Minimal RTF to text converter (enough for EES equation windows)."""
    out = []
    stack = []
    skip = False
    i, n = 0, len(rtf)
    word_re = re.compile(r"\\([a-zA-Z]+)(-?\d+)? ?")
    while i < n:
        c = rtf[i]
        if c == "{":
            stack.append(skip)
            m = re.match(r"\{\\(\*\\)?([a-zA-Z]+)", rtf[i:i + 40])
            if m and (m.group(1) or m.group(2) in _RTF_SKIP_DEST):
                skip = True
            i += 1
        elif c == "}":
            skip = stack.pop() if stack else False
            i += 1
        elif c == "\\":
            nxt = rtf[i + 1:i + 2]
            if nxt in ("\\", "{", "}"):
                if not skip:
                    out.append(nxt)
                i += 2
            elif nxt == "'":
                if not skip:
                    out.append(bytes.fromhex(rtf[i + 2:i + 4]).decode("cp1252", errors="replace"))
                i += 4
            elif nxt == "~":
                if not skip:
                    out.append("\u00a0")
                i += 2
            else:
                m = word_re.match(rtf, i)
                if not m:
                    i += 1
                    continue
                word, arg = m.group(1), m.group(2)
                i = m.end()
                if skip:
                    continue
                if word in ("par", "line"):
                    out.append("\n")
                elif word == "tab":
                    out.append("\t")
                elif word == "u" and arg is not None:
                    out.append(chr(int(arg) % 65536))
                    if i < n and rtf[i] == "?":    # skip the ANSI replacement char
                        i += 1
        elif c in "\r\n":
            i += 1                               # raw newlines are not content in RTF
        else:
            if not skip:
                out.append(c)
            i += 1
    return "".join(out)


# ---------------------------------------------------------------------------
# File sections
# ---------------------------------------------------------------------------

def parse_header(b):
    """Return (version, text_offset, text_length)."""
    if len(b) < 24:
        raise ValueError("file too short to be an EES file")
    version = b[1:1 + b[0]].decode("latin-1", errors="replace")
    n = struct.unpack_from("<I", b, 16)[0]
    if not version.startswith(("X", "W")) or not 0 < n <= len(b) - 20:
        raise ValueError("unrecognised EES header (version %r, text length %d)" % (version, n))
    return version, 20, n


_CONTROL_RE = re.compile("[\x00-\x08\x0b-\x1f\x7f]")      # every control character except tab and newline


def strip_control_chars(text):
    """EES pads its text fields with NUL bytes (and the RTF converter keeps them): CoolSolve cannot read them."""
    return _CONTROL_RE.sub("", text)


def decode_equations(raw):
    """Return (text, is_rtf) from the raw equation bytes."""
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):                 # EES 10.8 can store Unicode text (UTF-16 with BOM)
        text = raw.decode("utf-16", errors="replace")
    else:
        text = raw.decode("cp1252", errors="replace")
    is_rtf = text.lstrip().startswith("{\\rtf")
    if is_rtf:
        text = rtf_to_text(text)
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return strip_control_chars(text), is_rtf


# EES-only tags embedded in the text that have no meaning for CoolSolve.
_TAG_RES = [
    (re.compile(r"\{\$ID\$[^}]*\}"), "license ID tag {$ID$...}"),
    (re.compile(r"\{\$[A-Z]{2}\$[^}]*\}"), "EES display tag {$XX$...}"),
    (re.compile(r"\{\$DS[.,]\}"), "decimal-separator tag {$DS.}"),
]


def _split_code(text):
    """Yield (is_code, segment): comments "..", {..}, //.. and strings '..' are not code."""
    i, n, start = 0, len(text), 0
    while i < n:
        c = text[i]
        end = None
        if c == '"':
            end = text.find('"', i + 1)
        elif c == "{":
            end = text.find("}", i + 1)
        elif c == "'":
            end = text.find("'", i + 1)
        elif text.startswith("//", i):
            end = text.find("\n", i)
            end = n - 1 if end < 0 else end - 1
        if end is None:
            i += 1
            continue
        end = n - 1 if end < 0 else end
        if i > start:
            yield True, text[start:i]
        yield False, text[i:end + 1]
        i = start = end + 1
    if start < n:
        yield True, text[start:]


# The comma of a range clause separates its bounds: `SUM(expr, i=1,4)`, `PRODUCT(expr, i=1,4)`, `DUPLICATE i=1,4`.
# With the decimal comma EES writes these bounds with `;`, so a comma between two digits there is never a number.
_RANGE_COMMA_RES = [re.compile(r"(,\s*[A-Za-z_]\w*\s*=\s*[^,;()=\n]*?\d),(?=\d)"),
                    re.compile(r"^(\s*DUPLICATE\s+[A-Za-z_]\w*\s*=\s*[^,;\n]*?\d),(?=\d)", flags=re.I | re.M)]
_RANGE_MARK = "\ue000"


def _protect_range_commas(code):
    """Replace the comma between the bounds of a range clause by a mark (restored by the caller)."""
    for rx in _RANGE_COMMA_RES:
        code = rx.sub(lambda m: m.group(1) + _RANGE_MARK, code)
    return code


def normalize_decimal_comma(text):
    """Convert text typed with the EES 'comma' decimal separator to the dot convention.

    With the European setting EES writes `0,5` for numbers and `;` between the
    arguments of functions, procedures, arrays and DUPLICATE bounds.  The mode
    is detected by a `;` inside parentheses/brackets (impossible with dots) or
    a number with a decimal comma after `=`.  Comments and strings are kept.
    Returns (text, converted).
    """
    code = _protect_range_commas("".join(seg for is_code, seg in _split_code(text) if is_code))
    depth, semicolon_in_args = 0, False
    for ch in code:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth = max(0, depth - 1)
        elif ch == ";" and depth > 0:
            semicolon_in_args = True
            break
    comma_numbers = re.search(r"^(?!\s*DUPLICATE)[^\n]*=\s*-?\d+,\d", code, flags=re.I | re.M)
    if not (semicolon_in_args or comma_numbers):
        return text, False
    out, depth = [], 0
    for is_code, seg in _split_code(text):
        if not is_code:
            out.append(seg)
            continue
        seg = re.sub(r"(?<=\d),(?=\d)", ".", _protect_range_commas(seg))
        chars = []
        for line in seg.split("\n"):
            dup = re.match(r"\s*DUPLICATE\b", line, flags=re.I) is not None
            for ch in line:
                if ch in "([":
                    depth += 1
                elif ch in ")]":
                    depth = max(0, depth - 1)
                chars.append("," if ch == ";" and (depth > 0 or dup) else ch)
            chars.append("\n")
        out.append("".join(chars)[:-1].replace(_RANGE_MARK, ","))
    return "".join(out), True


def strip_ees_tags(text):
    removed = []
    for rx, label in _TAG_RES:
        found = rx.findall(text)
        if found:
            removed.append((label, found))
            text = rx.sub("", text)
    return text.rstrip() + "\n", removed


_UNIT_T = {0: "C", 1: "K", 2: "F", 3: "R"}
_UNIT_P = {0: "kPa", 1: "bar", 2: "Pa", 3: "MPa"}


def parse_unit_settings(b, off):
    """Decode the unit-system bytes that follow the equation text.

    Layout (EES 7+): [system, basis, T, P, trig, energy] + trig factor (Extended).
    Older files have no energy byte (5 bytes).  Codes verified against
    $UnitSystem directives: basis 0=mass 1=molar; T 0=C 1=K; P 0=kPa 1=bar
    2=Pa; trig 0=deg 1=rad; energy 0=J 1=kJ.  P 3=MPa is assumed.
    """
    for nbytes in (6, 5):
        f = read_extended(b, off + nbytes)
        if abs(f - TRIG_DEG_FACTOR) < 1e-8 or f == 1.0:     # EES stores 0.017453293
            s = b[off:off + nbytes]
            energy = ("J" if s[5] == 0 else "kJ") if nbytes == 6 else "kJ"
            system = "SI" if s[0] == 0 else "ENG"
            return {
                "system": system,
                "basis": "MOLE" if s[1] == 1 else "MASS",
                "T": _UNIT_T.get(s[2], "?%d" % s[2]),
                "P": _UNIT_P.get(s[3], "?%d" % s[3]),
                "trig": "RAD" if s[4] == 1 else "DEG",
                "energy": energy,
                "raw": s.hex(" "),
                "confident": system == "SI" and s[2] in (0, 1) and s[3] in (0, 1, 2) and nbytes == 6,
            }
    return None


def unit_directive(u):
    if not u:
        return None
    return "$UnitSystem %s %s %s %s %s %s" % (u["system"], u["basis"], u["trig"],
                                              u["P"].upper(), u["T"], u["energy"].upper())


# Variable names as CoolSolve reads them (docs/language_reference.md, "Numbers and identifiers"; closed gaps
# CS-GAP-NAME-SYMBOL, CS-GAP-NAME-PIPE, CS-GAP-ACCENTED-IDENT, CS-BUG-HASH-IDENT): a name starts with a letter
# or `_`, holds letters, digits, `_` and the symbols `%`, `|`, `#`; a `\` may join two parts of a name
# (`M\M_ex_duct_n`); the letters are ASCII or the Latin-1 letters (`é`, `ß`, `µ`, but neither `×` nor `÷`) and
# the signs `°`, `²`, `³`, `¹` may follow them.  EES stores names in Windows-1252, whose bytes 0xA0-0xFF are the
# Latin-1 characters.  A trailing `$` marks a string variable: its record holds no number and is skipped.
_NAME_LETTERS = "A-Za-z_\u00c0-\u00d6\u00d8-\u00ff\u00b5"
_NAME_CHARS = _NAME_LETTERS + "0-9%|#\u00b0\u00b2\u00b3\u00b9"
_NAME_PART = "[%s][%s]*" % (_NAME_LETTERS, _NAME_CHARS)
_NAME_BODY = r"%s(?:\\%s)*" % (_NAME_PART, _NAME_PART)
_NAME_RE = re.compile(_NAME_BODY + r"(\[[0-9, ]+\])?")
# EES's built-in constants (`g#`, `R#`...): EES keeps a record for them, but it holds a placeholder, not their value.
_BUILTIN_CONSTANTS = {"g#", "r#", "k#", "sigma#", "po#", "to#", "c#", "h#", "e#", "me#", "mp#", "mn#", "na#",
                      "avogadro#"}
_WORD_RE = re.compile(r"(?<![%s.])%s" % (_NAME_CHARS, _NAME_BODY))


def _code_text(text):
    """The code of the equations: no comments (`"..."`, `{...}`, `//...`) and no 'strings'."""
    return "".join(seg if is_code else " " for is_code, seg in _split_code(text))


def _words(text, code_only=False):
    """Lower-case set of the names that occur in `text` (in its code only when `code_only`)."""
    return {w.lower() for w in _WORD_RE.findall(_code_text(text) if code_only else text)}


def _record_at(b, i, words, code_words):
    """Decode the Variable Information record whose name length byte is at `i`, or return None.

    The record is recognised by its structure, not by one format byte: a Windows-1252 name that CoolSolve
    accepts and that occurs in the equations, four plausible Extended floats (value, guess, lower, upper)
    with lower <= upper, a display-format code (@+74: 0 to 5, EES 7.x writes 1 or 3 and EES 10.x also 0 and
    5) and a printable units string.  Blocks of zero bytes that start with a word of the text (a fluid name,
    a table name) have all-zero bounds, which no variable has (the default bounds are -inf and +inf).
    A record whose display format is not 3 (the only one this tool accepted before) must also name a
    variable of the *code*: EES keeps the records of deleted variables, whose names may survive in a comment
    (a legend of symbols), and their cleared values must not be taken for a solution.
    """
    n = b[i]
    if not 1 <= n <= 40 or b[i + 74] > 5:
        return None
    try:
        name = b[i + 1:i + 1 + n].decode("cp1252")
    except UnicodeDecodeError:
        return None
    if (not _NAME_RE.fullmatch(name) or name.lower() in _BUILTIN_CONSTANTS
            or not all(plausible_extended(b, i + k) for k in (31, 41, 51, 61))):
        return None
    units, ul = pascal_string(b, i + 77, 40)
    lower, upper = read_extended(b, i + 51), read_extended(b, i + 61)
    if units is None or not units.isprintable() or lower > upper or (lower == 0.0 and upper == 0.0):
        return None
    base = name.split("[")[0].lower()
    if base not in words or (b[i + 74] != 3 and base not in code_words):
        return None
    value = read_extended(b, i + 31)
    return {
        "name": name,
        "value": None if value == EES_NO_VALUE else value,
        "guess": read_extended(b, i + 41),
        "lower": lower,
        "upper": upper,
        "units": units.strip(),
        "offset": i,
        "end": i + 77 + 1 + ul,
    }


def scan_records(b, start, text):
    """Every Variable Information record from `start` on, duplicates included (EES 7.x-10.x).

    Each record starts with the variable name as a Delphi short string; relative
    to the length byte: value @+31, guess @+41, lower @+51, upper @+61 (Extended),
    format bytes @+71..76 (@+74 is the display format: 3 for most records, 1 for
    fixed-format variables, 0 or 5 in EES 8-10 files), units short string @+77
    (see `_record_at`).  The records of one file follow each other on a regular
    grid (122 bytes in EES 7-9, 338 in EES 9.7+).
    """
    words, code_words = _words(text), _words(text, code_only=True)
    records = []
    i, end = start, len(b) - 120
    while i < end:
        rec = _record_at(b, i, words, code_words)
        if rec:
            records.append(rec)
            i = rec["end"]
        else:
            i += 1
    return records


def records_block_end(records):
    """End of the first block of records, where the variables of the main program are stored.

    The records follow each other on the grid of the file (122 or 338 bytes); the records of deleted variables
    that are not decoded leave holes of a few strides.  Tables never lie inside this block, but they can lie
    between it and later blocks of records (copies for the parametric tables, modules).
    """
    steps = collections.Counter(nxt["offset"] - cur["offset"] for cur, nxt in zip(records, records[1:]))
    stride = steps.most_common(1)[0][0] if steps else 122
    end = records[0]["end"]
    for rec in records[1:]:
        if rec["offset"] - end > 12 * stride:
            break
        end = rec["end"]
    return end


def unique_variables(records):
    """The first record of each name (EES is case-insensitive), in file order."""
    seen, found = set(), []
    for rec in records:
        if rec["name"].lower() not in seen:
            seen.add(rec["name"].lower())
            found.append(rec)
    return found


def extract_variables(b, start, text):
    """The stored value, guess, bounds and units of each variable (see `scan_records`)."""
    return unique_variables(scan_records(b, start, text))


_HEADER_RE = re.compile(r"(?P<name>[A-Za-z_][A-Za-z0-9_$#]*(\[[0-9, ]*\])?)[^\n\w]?\n")


def _header_at(b, p):
    """Column header = short string 'name\n<units or junk>' in a 31-byte field (one junk byte may precede the \n)."""
    s, n = pascal_string(b, p, 40)
    m = _HEADER_RE.match(s) if s and n >= 2 else None
    if not m or not plausible_extended(b, p + 31):
        return None
    name, units = m.group("name"), s[m.end():]
    units = units.strip().strip("[]").strip()
    if not units.replace(" ", "").replace("/", "").replace("-", "").replace("^", "").isalnum():
        units = ""                     # EES leaves junk after the newline when units are empty
    return name, units


def _string_cells(b, q, nrows, stop):
    """Parse nrows string cells (int32 n+1, short string) starting at q; return end or None."""
    for _ in range(nrows):
        if q + 5 > stop:
            return None
        k = struct.unpack_from("<I", b, q)[0]
        s, n = pascal_string(b, q + 4, 255)
        if s is None or k != n + 1:
            return None
        q += 4 + 1 + n
    return q


def _cell_strings(b, q, nrows):
    """The nrows strings that `_string_cells` has validated at q."""
    out = []
    for _ in range(nrows):
        s, n = pascal_string(b, q + 4, 255)
        out.append(strip_control_chars(s))
        q += 5 + n
    return out


# After its header (31 bytes) and its nrows Extended cells, a column holds `per_row` state bytes (0 or 1) per row
# and `fixed` format bytes: (1, 3) in the Lookup tables of EES 7.9, (1, 0) and (0, 3) in other files.  A string
# column is then followed by its nrows strings.
_TAILS = ((1, 3), (1, 0), (0, 3))
_MAX_ROWS = 20000


def _numeric_cells_ok(b, q, nrows, per_row):
    """True if the nrows cells at q are plausible numbers and the state bytes after them are 0 or 1."""
    return (all(plausible_extended(b, q + 10 * r) for r in range(nrows))
            and all(x <= 1 for x in b[q + 10 * nrows:q + (10 + per_row) * nrows]))


def _link(b, p, p1, nrows):
    """How the column at p reaches the next header p1: (nrows, per_row, fixed, is_string) or None.

    `nrows` is None for the first column of a table: the row count is then solved from the gap.
    """
    gap = p1 - p - 31
    for per_row, fixed in _TAILS:                                  # string column (stricter test)
        for n in ([nrows] if nrows else range(1, min(_MAX_ROWS, gap // 10) + 1)):
            if _string_cells(b, p + 31 + (10 + per_row) * n + fixed, n, p1) == p1:
                return n, per_row, fixed, True
    for per_row, fixed in _TAILS:
        n, rest = divmod(gap - fixed, 10 + per_row)
        if rest == 0 and 0 < n <= _MAX_ROWS and nrows in (None, n) and _numeric_cells_ok(b, p + 31, n, per_row):
            return n, per_row, fixed, False
    return None


def extract_tables(b, start, rejected=None):
    """Find embedded Lookup/Parametric table grids (EES 7.x-10.x).

    A grid is stored column-wise: header short string 'name\nunits' in a 31-byte
    field, nrows Extended values, their state bytes (see `_TAILS`) and, for a
    string column, nrows strings.  Consecutive headers are chained by solving
    for the common row count.  `start` is the end of the variable records: the
    bytes before it are records, never a grid.  A chain whose cells are not all
    plausible numbers or strings is not a table: it is not returned, and is
    appended to `rejected` (offset, column names) when that list is given.
    """
    heads = [p for p in range(start, len(b) - 41)
             if 2 <= b[p] <= 40 and b"\n" in b[p + 1:p + 1 + b[p]] and _header_at(b, p)]
    tables, used = [], set()
    for k, p0 in enumerate(heads):
        if p0 in used:
            continue
        chain, links, nrows, p = [p0], [], None, p0
        for p1 in heads[k + 1:]:
            link = _link(b, p, p1, nrows)
            if not link:
                break
            nrows = link[0]
            links.append(link)
            chain.append(p1)
            p = p1
        if nrows is None or len(chain) < 2:
            continue                                            # single columns are ambiguous
        # The last column has no next header to measure it: try the layout of the other columns first.
        for per_row, fixed in dict.fromkeys([links[-1][1:3]] + list(_TAILS)):
            end = chain[-1] + 31 + (10 + per_row) * nrows + fixed
            string_end = _string_cells(b, end, nrows, len(b))
            if string_end is not None:
                links.append((nrows, per_row, fixed, True))
                end = string_end
                break
            if _numeric_cells_ok(b, chain[-1] + 31, nrows, per_row):
                links.append((nrows, per_row, fixed, False))
                break
        else:
            if rejected is not None:
                rejected.append((p0, [_header_at(b, c)[0] for c in chain]))
            continue
        cols = []
        for c, (_, per_row, fixed, is_string) in zip(chain, links):
            name, units = _header_at(b, c)
            if is_string:
                values = _cell_strings(b, c + 31 + (10 + per_row) * nrows + fixed, nrows)
            else:
                values = [read_extended(b, c + 31 + 10 * r) for r in range(nrows)]
            cols.append({"name": name, "units": units, "values": values})
        used.update(chain)
        tables.append({"offset": p0, "end": end, "nrows": nrows, "columns": cols})
    return tables


_TABLE_NAME_RE = re.compile(rb"([\x01-\x28])((?:Lookup|Table|Parametric|Integral)[ _A-Za-z0-9#\-]{0,30})")


def name_tables(b, tables, start, text):
    """Attach a name and a kind (lookup/parametric/integral) to each table."""
    lookup_refs = {m.lower(): m for m in re.findall(
        r"(?:LOOKUP\$?|INTERPOLATE\w*|DIFFERENTIATE\w*|NLOOKUP\w+|\w+LOOKUP|LOOKUPCOL\w*)\s*\(\s*'([^']+)'",
        text, flags=re.I)}
    candidates = []
    for m in _TABLE_NAME_RE.finditer(b, start):
        n = m.group(1)[0]
        s = b[m.start(2):m.start(2) + n]
        if len(s) == n and re.fullmatch(rb"(Lookup|Table|Parametric|Integral)[ _A-Za-z0-9#\-]*", s):
            candidates.append((m.start(), s.decode("cp1252", errors="replace")))
    for ref in lookup_refs.values():                  # user-renamed lookup tables
        key = bytes([len(ref)]) + ref.encode("cp1252", errors="replace")
        for m in re.finditer(re.escape(key), b[start:]):
            candidates.append((start + m.start(), ref))
    candidates.sort()
    prev_end, taken = start, set()
    for t in tables:
        before = [c for c in candidates if c[0] < t["offset"]]
        name = before[-1][1] if before else "table%d" % (tables.index(t) + 1)
        if name.lower() in taken:
            # The settings of a table can repeat the name of the previous one (its key column): the right
            # name is then the first one still free between the two grids.
            free = [c for c in before if c[0] >= prev_end and c[1].lower() not in taken]
            name = free[0][1] if free else name
        t["name"] = name
        taken.add(name.lower())
        prev_end = t["end"]
        low = t["name"].lower()
        if low in lookup_refs or low.startswith("lookup"):
            t["kind"] = "lookup"
        elif low.startswith("integral"):
            t["kind"] = "integral"
        else:
            t["kind"] = "parametric"
    return lookup_refs


def safe_table_name(name):
    s = re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_").lower()
    return s or "table"


# ---------------------------------------------------------------------------
# Feature scan (what will need attention during the CoolSolve import)
# ---------------------------------------------------------------------------

_FEATURE_PATTERNS = [
    ("PROCEDURE/FUNCTION definitions", r"^\s*(PROCEDURE|FUNCTION)\b"),
    ("MODULE/SUBPROGRAM (not supported)", r"^\s*(MODULE|SUBPROGRAM)\b"),
    ("DUPLICATE loops", r"\bDUPLICATE\b"),
    ("INTEGRAL (dynamic model)", r"\bINTEGRAL\s*\("),
    ("$IntegralTable", r"\$IntegralTable"),
    ("Lookup-table functions", r"\b(LOOKUP|INTERPOLATE\w*|TABLEVALUE|NLOOKUPROWS)\s*\("),
    ("Parametric-table functions", r"\b(TABLERUN#|\w+PARAMETRIC|NPARAMETRICROWS)\b"),
    ("$UnitSystem directive", r"\$UnitSystem"),
    ("Other $ directives", r"^\s*\$(?!UnitSystem|IntegralTable|if|ifnot|endif)\w+"),
    ("$Include / $Load (external library)", r"\$(Include|Load)\b"),
    ("Unit conversion functions", r"\b(CONVERT|CONVERTTEMP)\s*\("),
    ("Humid-air properties", r"\b(AirH2O|PSYCHPROPS)\b"),
    ("Unsupported fluids (NH3H2O, LiBrH2O, ...)", r"\b(NH3H2O|LiBrH2O|NH3-H2O|LiBr-H2O)\b"),
    ("Units annotations [..]", r"\[[^\]]*(kPa|bar|MPa|kJ|K\b|C\b|J\b|Pa\b)[^\]]*\]"),
    ("ERROR/WARNING calls", r"\b(CALL\s+)?(ERROR|WARNING)\s*\("),
    ("String variables", r"\w\$\s*="),
    ("Complex numbers", r"\$Complex"),
]


def scan_features(text):
    hits = []
    for label, rx in _FEATURE_PATTERNS:
        n = len(re.findall(rx, text, flags=re.I | re.M))
        if n:
            hits.append((label, n))
    return hits


def called_functions(text):
    """Names called like functions, minus those defined in the file."""
    code = re.sub(r'"[^"]*"|\{[^}]*\}|//[^\n]*', " ", text)
    calls = {m.lower() for m in re.findall(r"\b([A-Za-z_][A-Za-z0-9_]*\$?)\s*\(", code)}
    calls |= {m.lower() for m in re.findall(r"\bCALL\s+([A-Za-z_][A-Za-z0-9_]*)", code, flags=re.I)}
    defined = {m.lower() for m in re.findall(r"^\s*(?:FUNCTION|PROCEDURE)\s+([A-Za-z_][A-Za-z0-9_]*\$?)",
                                             code, flags=re.I | re.M)}
    return sorted(calls - defined), sorted(defined)


# ---------------------------------------------------------------------------
# Checks of the decoded variable records
# ---------------------------------------------------------------------------

def _equations(text):
    """Code that defines the stored variables: no comments, no `$` directives, no FUNCTION/PROCEDURE body
    (their variables are local and have no record); the equations of a MODULE do define records."""
    lines, in_body = [], False
    for line in _code_text(text).split("\n"):
        word = line.strip().lower()
        if re.match(r"(function|procedure)\b", word):
            in_body = True
        if not in_body and not word.startswith("$"):
            lines.append(line)
        elif re.fullmatch(r"end\s*", word):
            in_body = False
    return "\n".join(lines)


_ASSIGNMENT_RE = re.compile(r"^[ \t]*(?P<lhs>%s(?:\[[^\]\n]*\])?)[ \t]*=(?P<rhs>[^\n]*)$" % _NAME_BODY, re.M)


def defined_variables(text):
    """Lower-case names of the variables the main program defines with an equation `name = ...`."""
    return {m.group("lhs").split("[")[0].lower() for m in _ASSIGNMENT_RE.finditer(_equations(text))}


def missing_records(text, variables):
    """Variables defined in the equations that have no Variable Information record."""
    have = {v["name"].split("[")[0].lower() for v in variables}
    return sorted(defined_variables(text) - have)


_OPERAND_RE = re.compile(r"(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)|(?P<name>%s(?:\[\d+\])?)|(?P<op>[-+*/^()])|\s+"
                         % _NAME_BODY)


def _evaluate(rhs, values):
    """Value of a right-hand side made of numbers, `+ - * / ^`, parentheses and variables of `values`.

    Returns (value, largest operand) or None for anything else (function call, unknown or non-finite
    variable, division by zero).  The expression is evaluated by walking its syntax tree, never with eval().
    """
    rhs = re.sub(r"(?<=[\d)\s])\s*\[[^\]\n]*\]\s*$", "", rhs)          # a trailing units annotation [kPa]
    pieces, pos, scale = [], 0, 0.0
    while pos < len(rhs):
        m = _OPERAND_RE.match(rhs, pos)
        if not m:
            return None
        if m.group("name"):
            x = values.get(m.group("name").lower())
            if x is None:
                return None
            scale = max(scale, abs(x))
            pieces.append(repr(x))
        elif m.group("num"):
            scale = max(scale, float(m.group("num")))
            pieces.append(m.group("num"))
        elif m.group("op"):
            pieces.append("**" if m.group("op") == "^" else m.group("op"))
        pos = m.end()
    try:
        tree = ast.parse("".join(pieces), mode="eval")
        return _walk(tree.body), scale
    except (SyntaxError, ValueError, ZeroDivisionError, OverflowError, TypeError):
        return None


_BINARY = {ast.Add: lambda a, b: a + b, ast.Sub: lambda a, b: a - b, ast.Mult: lambda a, b: a * b,
           ast.Div: lambda a, b: a / b, ast.Pow: lambda a, b: a ** b}


def _walk(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        x = _walk(node.operand)
        return x if isinstance(node.op, ast.UAdd) else -x
    if isinstance(node, ast.BinOp) and type(node.op) in _BINARY:
        x = _BINARY[type(node.op)](_walk(node.left), _walk(node.right))
        if isinstance(x, float) and math.isfinite(x):
            return x
    raise ValueError("not a finite arithmetic expression")


def check_stored_solution(variables, text):
    """Stored values that cannot be the solution of the equations: [(name, reason)].

    EES does not always rewrite the Variable Information records when it saves a file, so a record can hold the
    value of an older run (docs/ees_import.md, 'Stored solution may be stale').  Only two cases can be
    detected here: a value that is inf or NaN, and a variable defined by an explicit equation of numbers and
    other stored variables (`x = 5`, `r = p1/p2`) whose stored value differs from the value of that equation.
    A stored solution that passes this check is not thereby correct.
    """
    suspect = [(v["name"], "stored value is %s" % fmt_num(v["value"])) for v in variables
               if v["value"] is not None and not math.isfinite(v["value"])]
    flagged = {name for name, _ in suspect}
    values = {v["name"].lower(): v["value"] for v in variables
              if v["value"] is not None and math.isfinite(v["value"])}
    assignments = [(m.group("lhs").replace(" ", ""), m.group("rhs")) for m in _ASSIGNMENT_RE.finditer(_equations(text))]
    counts = collections.Counter(lhs.lower() for lhs, _ in assignments)
    for lhs, rhs in assignments:
        key = lhs.lower()
        if counts[key] > 1 or key not in values or lhs in flagged:     # defined twice ($if branches): skip
            continue
        result = _evaluate(rhs.strip(), values)
        if result and abs(values[key] - result[0]) > 1e-4 * max(abs(values[key]), abs(result[0]), 1e-3 * result[1]):
            suspect.append((lhs, "stored value %s, but the equation `%s = %s` gives %s" % (
                fmt_num(values[key]), lhs, rhs.strip(), fmt_num(result[0]))))
    return suspect


# ---------------------------------------------------------------------------
# Main extraction
# ---------------------------------------------------------------------------

def extract(path):
    return extract_bytes(Path(path).read_bytes(), str(path))


def extract_bytes(b, label="<bytes>"):
    """Same as extract() for the content of a file already in memory (e.g. a zip member)."""
    version, off, n = parse_header(b)
    text, is_rtf = decode_equations(b[off:off + n])
    text, removed_tags = strip_ees_tags(text)
    text, decimal_comma = normalize_decimal_comma(text)
    settings_off = off + n
    units = parse_unit_settings(b, settings_off)
    records = scan_records(b, settings_off, text)
    variables = unique_variables(records)
    rejected_tables = []
    tables = extract_tables(b, records_block_end(records) if records else settings_off, rejected_tables)
    lookup_refs = name_tables(b, tables, settings_off, text)
    calls, defined = called_functions(text)
    return {
        "path": label, "size": len(b), "version": version, "is_rtf": is_rtf,
        "decimal_comma": decimal_comma,
        "text": text, "removed_tags": removed_tags, "units": units,
        "variables": variables, "tables": tables, "rejected_tables": rejected_tables, "lookup_refs": lookup_refs,
        "missing_records": missing_records(text, variables) if variables else [],
        "suspect_values": check_stored_solution(variables, text),
        "features": scan_features(text), "calls": calls, "defined": defined,
    }


STORED_SOLUTION_NOTE = ("The stored values are the Variable Information records of the .ees file. EES does not always "
                        "rewrite them when it saves, so a value can be stale (an older run) even when this tool "
                        "reports nothing: cross-check a few results by hand before using them as the reference.")


def _free_file(path, taken):
    """`path`, or `path` with a numeric suffix when an earlier table of this run already uses that file name."""
    free, k = path, 2
    while free.name in taken:
        free = path.with_name("%s_%d%s" % (path.stem, k, path.suffix))
        k += 1
    taken.add(free.name)
    return free


def write_outputs(res, outdir, name):
    outdir = Path(outdir)
    ref = outdir / "reference"
    ref.mkdir(parents=True, exist_ok=True)
    warnings = []
    text = res["text"]

    # Unit system: make the EES settings explicit with a $UnitSystem directive.
    directive = unit_directive(res["units"])
    if directive and not re.search(r"\$UnitSystem", text, flags=re.I):
        text = directive + "\n\n" + text
    u = res["units"]
    if u is None:
        warnings.append("Unit-system settings not decoded: check units manually.")
    elif (u["T"], u["P"], u["energy"], u["basis"], u["trig"]) != ("C", "Pa", "J", "MASS", "DEG"):
        warnings.append("EES unit system is %s, CoolSolve uses SI MASS DEG PA C J: the model must be "
                        "converted (see docs/ees_import.md, 'Unit conversion')." % directive)
    if u and not u["confident"]:
        warnings.append("Unit-system decoding is not fully verified for this file (raw bytes %s)." % u["raw"])

    # Lookup tables -> <name>-<table>.csv (CoolSolve convention), references renamed.
    lookup_files, param_files, taken = [], [], set()
    for t in res["tables"]:
        cols = t["columns"]
        if t["kind"] == "lookup":
            safe = safe_table_name(t["name"])
            if safe != t["name"]:
                text = re.sub(r"'%s'" % re.escape(t["name"]), "'%s'" % safe, text, flags=re.I)
                warnings.append("Lookup table '%s' renamed '%s' (file name safety); equation "
                                "references updated." % (t["name"], safe))
            fname = outdir / ("%s-%s.csv" % (name, safe))
        else:
            fname = ref / ("ees_parametric-%s.csv" % safe_table_name(t["name"]))
        wanted, fname = fname, _free_file(fname, taken)
        if fname != wanted:
            warnings.append("Several tables are named '%s': the next one is written to %s instead of "
                            "overwriting the first -- check which one the equations use." % (t["name"], fname.name))
        (lookup_files if t["kind"] == "lookup" else param_files).append((t, fname.name))
        with open(fname, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow([c["name"] for c in cols])
            for r in range(t["nrows"]):
                w.writerow([fmt_num(c["values"][r]) for c in cols])
    for offset, columns in res["rejected_tables"]:
        warnings.append("A table-like structure at byte %d (columns %s) failed validation and was not "
                        "extracted: check the file in EES if it holds a table." % (offset, ", ".join(columns[:6])))
    for ref_name in res["lookup_refs"].values():
        if not any(t["name"].lower() == ref_name.lower() for t, _ in lookup_files):
            warnings.append("Lookup table '%s' is used but was not found in the file: it was probably "
                            "loaded from an external .lkt/.csv/.txt file -- look for it next to the .ees "
                            "file and convert it to %s-%s.csv." % (ref_name, name, safe_table_name(ref_name)))

    (outdir / ("%s.eescode" % name)).write_text(text, encoding="utf-8")

    # Initial guesses: the stored EES value (last solution) or, if zero, the guess.
    variables = res["variables"]
    if variables:
        with open(outdir / ("%s.initials" % name), "w", encoding="utf-8") as f:
            for v in variables:
                x = v["value"] if v["value"] not in (None, 0.0) and math.isfinite(v["value"]) else v["guess"]
                if math.isfinite(x) and not (v["value"] is None and x == 1.0):   # skip EES default guess
                    f.write("%s=%s\n" % (v["name"], fmt_num(x)))
        with open(ref / "ees_variables.csv", "w", newline="", encoding="utf-8") as f:
            f.write("# %s\n" % STORED_SOLUTION_NOTE)
            if res["suspect_values"]:
                f.write("# suspect: %s\n" % ", ".join(n for n, _ in res["suspect_values"]))
            w = csv.writer(f)
            w.writerow(["name", "value", "guess", "lower", "upper", "units"])
            for v in variables:
                w.writerow([v["name"], fmt_num(v["value"]), fmt_num(v["guess"]),
                            fmt_num(v["lower"]), fmt_num(v["upper"]), v["units"]])
    else:
        warnings.append("No variable records decoded (EES %s layout?): no initials / reference "
                        "solution available." % res["version"])
    if variables and res["missing_records"]:
        shown = ", ".join(res["missing_records"][:10]) + (", ..." if len(res["missing_records"]) > 10 else "")
        warnings.append("%d variable(s) defined in the equations have no Variable Information record (%s): "
                        "reference/ees_variables.csv does not cover them." % (len(res["missing_records"]), shown))
    if res["suspect_values"]:
        shown = "; ".join("%s: %s" % sv for sv in res["suspect_values"][:3]) + (
            "; ..." if len(res["suspect_values"]) > 3 else "")
        warnings.append("Stored solution may be STALE: %d value(s) are inconsistent with the equations (%s), either "
                        "stale records or the result of a run of a Parametric table. %s"
                        % (len(res["suspect_values"]), shown, STORED_SOLUTION_NOTE))
    if variables and sum(v["value"] is None for v in variables) > len(variables) / 2:
        warnings.append("Most stored values are -9999 (cleared): the file was last run from a Parametric "
                        "table or not solved -- use the parametric table(s) as reference instead.")
    elif variables and all(v["value"] == 0 or v["value"] == v["guess"] for v in variables):
        warnings.append("Stored values equal the guesses: the file was probably saved without "
                        "being solved -- the reference solution is NOT meaningful.")

    for label, n in res["features"]:
        if label.startswith(("MODULE", "Unsupported", "$Include", "Complex")):
            warnings.append("Feature needing attention: %s (%d occurrence(s))." % (label, n))

    report = render_report(res, name, directive, lookup_files, param_files, warnings)
    (ref / "ees_extract_report.md").write_text(report, encoding="utf-8")
    return warnings


def render_report(res, name, directive, lookup_files, param_files, warnings):
    u = res["units"]
    L = ["# EES extraction report: `%s`" % name, "",
         "| Item | Value |", "|---|---|",
         "| Source file | `%s` |" % os.path.basename(res["path"]),
         "| File size | %d bytes |" % res["size"],
         "| EES version | %s |" % res["version"],
         "| Equations format | %s%s |" % ("RTF (converted to text)" if res["is_rtf"] else "plain text",
                                         ", decimal comma converted to dot (`0,5` → `0.5`, `;` → `,` in "
                                         "argument lists)" if res["decimal_comma"] else ""),
         "| Equation-window lines | %d |" % res["text"].count("\n"),
         "| Unit system | `%s`%s |" % (directive, "" if (u and u["confident"]) else " (unverified)"),
         "| Variables decoded | %d |" % len(res["variables"]),
         "| Lookup tables | %s |" % (", ".join("%s (%d×%d) → `%s`" % (t["name"], t["nrows"], len(t["columns"]), f)
                                         for t, f in lookup_files) or "none"),
         "| Parametric tables | %s |" % (", ".join("%s (%d×%d) → `reference/%s`" % (t["name"], t["nrows"], len(t["columns"]), f)
                                             for t, f in param_files) or "none"),
         "| Extracted on | %s with `tools/ees_extract.py` |" % datetime.date.today().isoformat(),
         ""]
    if res["removed_tags"]:
        L += ["## Removed EES tags", ""]
        for label, found in res["removed_tags"]:
            L.append("- %s: `%s`" % (label, "`, `".join(x.replace("|", "\\|") for x in found)))
        L.append("")
    if res["features"]:
        L += ["## Features found in the equations", "", "| Feature | Count |", "|---|---|"]
        L += ["| %s | %d |" % (label, n) for label, n in res["features"]]
        L.append("")
    if res["defined"]:
        L += ["## Functions/procedures defined in the file", "", ", ".join("`%s`" % d for d in res["defined"]), ""]
    if res["calls"]:
        L += ["## Functions called (built-in or external)", "",
              "Check that each one exists in CoolSolve (docs/language_reference.md) or is provided by "
              "a library model:", "", ", ".join("`%s`" % c for c in res["calls"]), ""]
    L += ["## Stored solution", "", STORED_SOLUTION_NOTE, ""]
    if res["suspect_values"]:
        L += ["Values that cannot be the solution of the equations (also listed in the header of "
              "`ees_variables.csv`):", ""] + ["- `%s`: %s" % sv for sv in res["suspect_values"]] + [""]
    L += ["## Warnings", ""]
    L += ["- %s" % w for w in warnings] if warnings else ["- none"]
    L.append("")
    return "\n".join(L)


def print_info(res):
    u = res["units"]
    print("%s: EES %s, %d bytes, %s, %d lines" % (res["path"], res["version"], res["size"],
                                                  "RTF" if res["is_rtf"] else "text", res["text"].count("\n")))
    print("  units      : %s%s" % (unit_directive(u), "" if (u and u["confident"]) else " (unverified)"))
    if res["decimal_comma"]:
        print("  format     : decimal comma converted to dot")
    print("  variables  : %d decoded" % len(res["variables"]))
    if res["missing_records"]:
        print("  no record  : %d defined variable(s): %s" % (
            len(res["missing_records"]), ", ".join(res["missing_records"][:8])))
    if res["suspect_values"]:
        print("  STALE?     : %d stored value(s) inconsistent with the equations: %s" % (
            len(res["suspect_values"]), ", ".join(n for n, _ in res["suspect_values"][:8])))
    for t in res["tables"]:
        print("  table      : %-12s %-10s %d rows x %d cols: %s" % (
            t["name"], t["kind"], t["nrows"], len(t["columns"]), ", ".join(c["name"] for c in t["columns"])))
    for label, n in res["features"]:
        print("  feature    : %s (%d)" % (label, n))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ees_file", nargs="+", help="binary EES file(s) (.ees)")
    ap.add_argument("-o", "--outdir", help="output folder (default: <stem>_extracted next to the file)")
    ap.add_argument("--name", help="base name of the generated files (default: sanitised file stem)")
    ap.add_argument("--info", action="store_true", help="print a summary only, write nothing")
    args = ap.parse_args(argv)
    status = 0
    for p in args.ees_file:
        try:
            res = extract(p)
        except (OSError, ValueError) as e:
            print("ERROR %s: %s" % (p, e), file=sys.stderr)
            status = 1
            continue
        if args.info:
            print_info(res)
            continue
        name = args.name or safe_table_name(Path(p).stem)
        outdir = args.outdir or str(Path(p).with_name(Path(p).stem + "_extracted"))
        if len(args.ees_file) > 1 and args.outdir:
            outdir = os.path.join(args.outdir, name)
        warnings = write_outputs(res, outdir, name)
        print("%s -> %s/ (%d variables, %d tables, %d warnings)" % (
            p, outdir, len(res["variables"]), len(res["tables"]), len(warnings)))
        for w in warnings:
            print("  WARNING: " + w)
    return status


if __name__ == "__main__":
    sys.exit(main())
