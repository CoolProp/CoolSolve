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
    reference/ees_parametric-<t>.csv  Parametric tables
    reference/ees_extract_report.md   provenance, unit system, warnings

Companion script: tools/compare_solution.py compares a CoolSolve .sol file with
reference/ees_variables.csv.
"""

import argparse
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


def decode_equations(raw):
    """Return (text, is_rtf) from the raw equation bytes."""
    text = raw.decode("cp1252", errors="replace")
    is_rtf = text.lstrip().startswith("{\\rtf")
    if is_rtf:
        text = rtf_to_text(text)
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return text, is_rtf


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


def normalize_decimal_comma(text):
    """Convert text typed with the EES 'comma' decimal separator to the dot convention.

    With the European setting EES writes `0,5` for numbers and `;` between the
    arguments of functions, procedures, arrays and DUPLICATE bounds.  The mode
    is detected by a `;` inside parentheses/brackets (impossible with dots) or
    a number with a decimal comma after `=`.  Comments and strings are kept.
    Returns (text, converted).
    """
    code = "".join(seg for is_code, seg in _split_code(text) if is_code)
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
        seg = re.sub(r"(?<=\d),(?=\d)", ".", seg)
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
        out.append("".join(chars)[:-1])
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


_NAME_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*(\[[0-9, ]+\])?")


def extract_variables(b, start, text):
    """Decode the Variable Information records (EES 7.x-10.x).

    Each record starts with the variable name as a Delphi short string; relative
    to the length byte: value @+31, guess @+41, lower @+51, upper @+61 (Extended),
    format bytes @+71..76 (byte @+74 is always 3), units short string @+77.
    """
    words = {w.lower() for w in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", text)}
    found, seen = [], set()
    i, end = start, len(b) - 120
    while i < end:
        n = b[i]
        if 1 <= n <= 40 and b[i + 74] == 3:
            try:
                name = b[i + 1:i + 1 + n].decode("ascii")
            except UnicodeDecodeError:
                name = None
            if name and _NAME_RE.fullmatch(name) and all(plausible_extended(b, i + k) for k in (31, 41, 51, 61)):
                units, ul = pascal_string(b, i + 77, 40)
                lower, upper = read_extended(b, i + 51), read_extended(b, i + 61)
                base = name.split("[")[0].lower()
                if units is not None and units.isprintable() and lower <= upper and base in words:
                    if name.lower() not in seen:
                        seen.add(name.lower())
                        value = read_extended(b, i + 31)
                        found.append({
                            "name": name,
                            "value": None if value == EES_NO_VALUE else value,
                            "guess": read_extended(b, i + 41),
                            "lower": lower,
                            "upper": upper,
                            "units": units.strip(),
                        })
                    i += 77 + 1 + ul
                    continue
        i += 1
    return found


_HEADER_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_$#]*(\[[0-9, ]*\])?\n")


def _header_at(b, p):
    """Column header = short string 'name\n<units or junk>' in a 31-byte field."""
    s, n = pascal_string(b, p, 40)
    if not s or n < 2 or not _HEADER_RE.match(s) or not plausible_extended(b, p + 31):
        return None
    name, _, units = s.partition("\n")
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


def extract_tables(b, start):
    """Find embedded Lookup/Parametric table grids (EES 7.x-10.x).

    A grid is stored column-wise: header short string 'name\nunits' in a 31-byte
    field, nrows Extended values, then 6 trailing bytes (Lookup tables) or 3
    (Parametric tables; string columns are followed by their nrows strings).
    Consecutive headers are chained by solving for the common row count.
    """
    heads = [p for p in range(start, len(b) - 41)
             if 2 <= b[p] <= 40 and b"\n" in b[p + 1:p + 1 + b[p]] and _header_at(b, p)]
    tables, used = [], set()
    for k, p0 in enumerate(heads):
        if p0 in used:
            continue
        chain, nrows, p = [p0], None, p0
        for p1 in heads[k + 1:]:
            link = None
            if nrows is None:                                  # string column first (stricter test)
                for n in range(1, min(20000, (p1 - p - 34) // 10) + 1):
                    if _string_cells(b, p + 31 + 10 * n + 3, n, p1) == p1:
                        link = (n, "str")
                        break
            elif _string_cells(b, p + 31 + 10 * nrows + 3, nrows, p1) == p1:
                link = (nrows, "str")
            for trail in (6, 3):
                if link:
                    break
                gap = p1 - p - 31 - trail
                if nrows is None and gap > 0 and gap % 10 == 0:
                    link = (gap // 10, None)
                elif nrows is not None and gap == 10 * nrows:
                    link = (nrows, None)
            if not link or link[0] > 20000:
                break
            nrows = link[0]
            chain.append(p1)
            p = p1
        if nrows is None or len(chain) < 2:
            continue                                            # single columns are ambiguous
        cols = []
        for q in chain:
            name, units = _header_at(b, q)
            vals = [read_extended(b, q + 31 + 10 * r) if plausible_extended(b, q + 31 + 10 * r) else None
                    for r in range(nrows)]
            strs = None
            end_str = _string_cells(b, q + 31 + 10 * nrows + 3, nrows, len(b))
            if end_str is not None:
                strs, qq = [], q + 31 + 10 * nrows + 3
                for _ in range(nrows):
                    s, n = pascal_string(b, qq + 4, 255)
                    strs.append(s)
                    qq += 5 + n
            cols.append({"name": name, "units": units, "values": strs or vals})
        used.update(chain)
        tables.append({"offset": p0, "nrows": nrows, "columns": cols})
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
    for t in tables:
        before = [c for c in candidates if c[0] < t["offset"]]
        t["name"] = before[-1][1] if before else "table%d" % (tables.index(t) + 1)
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
    variables = extract_variables(b, settings_off, text)
    tables = extract_tables(b, settings_off)
    lookup_refs = name_tables(b, tables, settings_off, text)
    calls, defined = called_functions(text)
    return {
        "path": label, "size": len(b), "version": version, "is_rtf": is_rtf,
        "decimal_comma": decimal_comma,
        "text": text, "removed_tags": removed_tags, "units": units,
        "variables": variables, "tables": tables, "lookup_refs": lookup_refs,
        "features": scan_features(text), "calls": calls, "defined": defined,
    }


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
    lookup_files, param_files = [], []
    for t in res["tables"]:
        cols = t["columns"]
        if t["kind"] == "lookup":
            safe = safe_table_name(t["name"])
            if safe != t["name"]:
                text = re.sub(r"'%s'" % re.escape(t["name"]), "'%s'" % safe, text, flags=re.I)
                warnings.append("Lookup table '%s' renamed '%s' (file name safety); equation "
                                "references updated." % (t["name"], safe))
            fname = outdir / ("%s-%s.csv" % (name, safe))
            lookup_files.append((t, fname.name))
        else:
            fname = ref / ("ees_parametric-%s.csv" % safe_table_name(t["name"]))
            param_files.append((t, fname.name))
        with open(fname, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow([c["name"] for c in cols])
            for r in range(t["nrows"]):
                w.writerow([fmt_num(c["values"][r]) for c in cols])
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
            w = csv.writer(f)
            w.writerow(["name", "value", "guess", "lower", "upper", "units"])
            for v in variables:
                w.writerow([v["name"], fmt_num(v["value"]), fmt_num(v["guess"]),
                            fmt_num(v["lower"]), fmt_num(v["upper"]), v["units"]])
    else:
        warnings.append("No variable records decoded (EES %s layout?): no initials / reference "
                        "solution available." % res["version"])
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
