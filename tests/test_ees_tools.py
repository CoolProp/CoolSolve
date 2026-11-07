#!/usr/bin/env python3
"""
test_ees_tools.py -- unit tests of the EES import tools (tools/ees_extract.py, tools/compare_solution.py).

The tests build small synthetic EES files byte by byte (a few hundred bytes, no copy of a real EES file), so
they need neither EES files, nor a CoolSolve build, nor pytest.  Each `test_*` function covers one entry of the
"CS-BUG-EXTRACT-*" / "CS-BUG-COMPARE-*" family of docs/model_library_support.md; the negative tests check that
a malformed or truncated file gives a clear error or warning and never a crash or a silent loss.

Usage:
    python3 tests/test_ees_tools.py [-v]

Exit code 0 if all tests pass, 1 otherwise.
"""

import contextlib
import io
import math
import os
import struct
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import compare_solution as cmp  # noqa: E402
import ees_extract as ex  # noqa: E402

VERBOSE = "-v" in sys.argv or "--verbose" in sys.argv

# ── Synthetic EES file builder ───────────────────────────────────────────────


def ext(x):
    """80-bit x87 Extended (Delphi) encoding of a float, as EES stores it."""
    if x == 0:
        return bytes(10)
    if math.isinf(x):
        return (1 << 63).to_bytes(8, "little") + (0x73E3 | (0x8000 if x < 0 else 0)).to_bytes(2, "little")
    m, e = math.frexp(abs(x))
    return int(m * 2 ** 64).to_bytes(8, "little") + ((e - 1 + 16383) | (0x8000 if x < 0 else 0)).to_bytes(2, "little")


def pascal(s, size):
    """Delphi short string in a field of `size` bytes."""
    raw = s.encode("cp1252")
    return (bytes([len(raw)]) + raw).ljust(size, b"\0")


def record(name, value, guess=1.0, lower=-math.inf, upper=math.inf, units="", fmt=3, digits=3):
    """One 122-byte Variable Information record (EES 7)."""
    head = pascal(name, 31) + ext(value) + ext(guess) + ext(lower) + ext(upper)
    return (head + bytes([0, 0, 0, fmt, digits, 0]) + pascal(units, 45)).ljust(122, b"\0")


def column(name, values=None, strings=None, nrows=None, per_row=1, fixed=3):
    """One table column: 31-byte header, nrows cells, `per_row` state bytes per row, `fixed` format bytes and,
    for a string column, nrows strings (uint32 len+1, short string)."""
    n = nrows or len(values or strings)
    cells = [ext(v) for v in (values if values is not None else [0.0] * n)]
    out = pascal(name + "\n", 31) + b"".join(cells) + bytes([1]) * (per_row * n) + bytes([3, 3, 0][:fixed])
    for s in strings or []:
        raw = s.encode("cp1252")
        out += struct.pack("<I", len(raw) + 1) + bytes([len(raw)]) + raw
    return out


def build(text, records=(), tables=b"", before_records=b"", raw_text=None, version="X7.888"):
    """A whole file: header, equations, unit settings, [before_records], records, [tables], padding."""
    body = raw_text if raw_text is not None else text.replace("\n", "\r\n").encode("cp1252")
    head = (bytes([len(version)]) + version.encode()).ljust(16, b"\0") + struct.pack("<I", len(body))
    settings = bytes([0, 0, 0, 2, 0, 0]) + ext(math.pi / 180)           # SI, mass, C, Pa, deg, J
    return head + body + settings + bytes(30) + before_records + b"".join(records) + tables + bytes(200)


# A file with a unit system, three explicit equations and their records.
TEXT = "a = 5\nb = 2*a\nc = a*b\n"
RECORDS = [record("a", 5.0), record("b", 10.0), record("c", 50.0)]

passed = 0
failed = 0
errors = []


def extract(b):
    return ex.extract_bytes(b, "fixture.ees")


def written(res, name="m"):
    """Run write_outputs in a temp folder; return (warnings, {relative path: text})."""
    with tempfile.TemporaryDirectory() as d:
        warnings = ex.write_outputs(res, d, name)
        files = {str(p.relative_to(d)): p.read_text(encoding="utf-8")
                 for p in Path(d).rglob("*") if p.is_file()}
    return warnings, files


# ── CS-BUG-EXTRACT-FMT-BYTE ──────────────────────────────────────────────────

def test_records_with_any_valid_format_byte_are_kept():
    text = "a = 1\nb = 2\nc = 3\nd = 4\ne = 5\n"
    recs = [record("a", 1.0, fmt=3), record("b", 2.0, fmt=1), record("c", 3.0, fmt=0, digits=3),
            record("d", 4.0, fmt=5), record("e", 5.0, fmt=79)]
    names = [v["name"] for v in extract(build(text, recs))["variables"]]
    assert names == ["a", "b", "c", "d"], names          # 79 is no display format: the block is not a record


def test_junk_blocks_are_not_records():
    text = "a = 1\nb = 2\nc = 3\n"
    recs = [record("a", 1.0),
            record("b", 2.0, lower=0.0, upper=0.0, fmt=0),     # a word of the text followed by zero bytes
            record("c", 3.0, lower=5.0, upper=1.0)]            # lower > upper
    assert [v["name"] for v in extract(build(text, recs))["variables"]] == ["a"]


def test_records_without_name_in_the_code_are_ignored():
    text = '{ legend: ghost is not a variable }\na = 1\n'
    recs = [record("a", 1.0), record("ghost", 3.0, fmt=1)]
    assert [v["name"] for v in extract(build(text, recs))["variables"]] == ["a"]


def test_missing_record_is_reported():
    text = "a = 1\nb = a + 1\nfluid$ = 'R22'\nPROCEDURE p(x : y)\ny = x\nEND\n"
    res = extract(build(text, [record("a", 1.0)]))
    assert res["missing_records"] == ["b"], res["missing_records"]      # not fluid$, not the procedure's y
    warnings, _ = written(res)
    assert any("no Variable Information record (b)" in w for w in warnings), warnings
    assert extract(build("a = 1\n", [record("a", 1.0)]))["missing_records"] == []


def test_equations_of_a_module_define_records():
    text = "MODULE m(x : y)\ny = 2*x\nEND\nCALL m(1 : z)\nz = 3\n"
    assert extract(build(text, [record("y", 2.0)]))["missing_records"] == ["z"]


# ── CS-BUG-EXTRACT-NAME-CHARS ────────────────────────────────────────────────

def test_names_with_symbols_are_records():
    names = ["M\\M_ex", "C%", "K|star", "Nu#_f", "T_a_réf", "n°_90°elbow", "x_µ", "s[1]"]
    text = "\n".join("%s = 1" % n.replace("[1]", "[i]") for n in names) + "\n"
    recs = [record(n, float(i + 1)) for i, n in enumerate(names)]
    got = [v["name"] for v in extract(build(text, recs))["variables"]]
    assert got == names, got


def test_names_coolsolve_refuses_are_not_records():
    for bad in ["a b", "1abc", "a\\", "a\\1", "×x", "a×", "x$", "a.b", "%x", "a'b"]:
        assert not ex._NAME_RE.fullmatch(bad), bad
    for good in ["a", "_a", "a%", "a|b", "a#1", "é", "ÿ", "a\\b_c", "T[1,2]"]:
        assert ex._NAME_RE.fullmatch(good), good
    res = extract(build("x = 1\n", [record("x$", 1.0), record("x", 1.0)]))
    assert [v["name"] for v in res["variables"]] == ["x"]               # the record of a string variable


def test_builtin_constants_are_not_records():
    res = extract(build("g = g#\nx = 2\n", [record("g#", 1.0), record("g", 9.80665), record("x", 2.0)]))
    assert [v["name"] for v in res["variables"]] == ["g", "x"]


# ── CS-BUG-EXTRACT-DECIMAL-RANGE ─────────────────────────────────────────────

def test_range_comma_is_not_a_decimal_point():
    text = ("MM = sum(x[i]*MM[i],i=1,4)\ns = sum(x[i]*s[i], i=1,4)/MM\nduplicate i = 1,4\nx[i] = 0.25\nend\n"
            "r = 2.16\nC% = 86.9\n")
    out, converted = ex.normalize_decimal_comma(text)
    assert (out, converted) == (text, False)


def test_decimal_comma_files_still_convert():
    text = ("A = 0,5\nMM = sum(x[i]*y[i];i=1;4)\nduplicate i = 1;4\nx[i] = 0,25\nend\n"
            "h = enthalpy(water;T=20;P=1,5)\n\"comment 3,5 and i=1,4\"\n")
    out, converted = ex.normalize_decimal_comma(text)
    assert converted
    assert "A = 0.5" in out and "x[i] = 0.25" in out
    assert "sum(x[i]*y[i],i=1,4)" in out and "duplicate i = 1,4" in out
    assert "enthalpy(water,T=20,P=1.5)" in out
    assert '"comment 3,5 and i=1,4"' in out                              # comments are kept as they are


def test_range_is_kept_in_a_file_that_mixes_both_conventions():
    out, converted = ex.normalize_decimal_comma("A = 0,5\ns = sum(x[i],i=1,4)\n")
    assert converted and out == "A = 0.5\ns = sum(x[i],i=1,4)\n", out


# ── CS-BUG-EXTRACT-TABLE-FP ──────────────────────────────────────────────────

def test_no_phantom_table_in_the_records():
    # Two header-like columns inside the record area (a variable name and a string): not a table.
    fake = column("v", [1.0, 2.0, 3.0]) + column("Regime$", [4.0, 5.0, 6.0])
    b = build(TEXT, RECORDS, before_records=fake)
    assert extract(b)["tables"] == []
    start = 20 + len(TEXT.replace("\n", "\r\n"))
    assert len(ex.extract_tables(b, start)) == 1, "the scan of the whole file finds the fake grid"


def test_real_table_after_the_records_is_found():
    tab = column("x", [1.0, 2.0, 3.0]) + column("y", [4.0, 5.5, 7.0])
    (t,) = extract(build(TEXT, RECORDS, tables=tab))["tables"]
    assert (t["nrows"], [c["name"] for c in t["columns"]]) == (3, ["x", "y"])
    assert t["columns"][1]["values"] == [4.0, 5.5, 7.0]


def test_table_between_two_blocks_of_records_is_kept():
    tab = column("x", [1.0, 2.0, 3.0]) + column("y", [4.0, 5.5, 7.0])
    b = build(TEXT, RECORDS, tables=bytes(100) + tab + bytes(3000) + b"".join(RECORDS))
    assert [t["nrows"] for t in extract(b)["tables"]] == [3]
    assert len(extract(b)["variables"]) == 3                               # the second block repeats the names


def test_table_that_fails_validation_is_rejected_with_a_warning():
    bad = column("x", [1.0, 2.0, 3.0]) + column("y", [4.0, 5.0, 6.0])
    bad = bad[:-13] + bytes([0x01] * 10) + bad[-3:]            # third cell of the last column: no float
    res = extract(build(TEXT, RECORDS, tables=bad))
    assert res["tables"] == [] and res["rejected_tables"][0][1] == ["x", "y"]
    warnings, files = written(res)
    assert any("failed validation" in w for w in warnings), warnings
    assert not any("parametric" in f for f in files)


def test_tables_with_the_same_name_do_not_overwrite_each_other():
    name = bytes([7]) + b"Table 1"
    tabs = (name + column("x", [1.0, 2.0]) + column("y", [3.0, 4.0]) + bytes(40) + name
            + column("u", [5.0, 6.0]) + column("w", [7.0, 8.0]))
    res = extract(build(TEXT, RECORDS, tables=tabs))
    assert len(res["tables"]) == 2
    warnings, files = written(res)
    assert sorted(f for f in files if "parametric" in f) == ["reference/ees_parametric-table_1.csv",
                                                             "reference/ees_parametric-table_1_2.csv"]
    assert any("Several tables are named" in w for w in warnings), warnings


# ── CS-BUG-EXTRACT-LOOKUP-STRING ─────────────────────────────────────────────

def lookup_tables(per_row, fixed):
    kw = dict(per_row=per_row, fixed=fixed)
    first = (bytes([8]) + b"activite" + bytes(20) + column("activite", strings=["Sommeil", "Repos_assis", "Sport"], **kw)
             + column("alpha_E", [41.0, 60.0, 350.0], **kw))
    second = (bytes(30) + bytes([6]) + b"veture" + bytes(20)
              + column("veture", strings=["Tropicale", "Eté", "Intérieur (légère)", "Hiver"], **kw)
              + column("Column2", [0.35, 0.5, 0.7, 1.5], **kw))
    return first + second + bytes(30)


LOOKUP_TEXT = ("a = 1\nrow = LOOKUP$ROW('activite','activite',act$)\nalpha = lookup('activite',row,2)\n"
               "I = lookup('veture',2,2)\n")


def test_lookup_tables_with_a_string_key_column():
    for per_row, fixed in ((1, 3), (0, 3)):                     # EES 7.9 and older layouts
        res = extract(build(LOOKUP_TEXT, [record("a", 1.0)], tables=lookup_tables(per_row, fixed)))
        assert [(t["name"], t["kind"], t["nrows"]) for t in res["tables"]] == [
            ("activite", "lookup", 3), ("veture", "lookup", 4)], (per_row, res["tables"])
        assert res["tables"][0]["columns"][0]["values"] == ["Sommeil", "Repos_assis", "Sport"]
        assert res["tables"][1]["columns"][0]["values"][1:3] == ["Eté", "Intérieur (légère)"]
        assert res["tables"][1]["columns"][1]["values"] == [0.35, 0.5, 0.7, 1.5]
        warnings, files = written(res)
        assert files["m-activite.csv"] == "activite,alpha_E\nSommeil,41\nRepos_assis,60\nSport,350\n"
        assert files["m-veture.csv"].splitlines()[3] == "Intérieur (légère),0.7"
        assert not any("not found" in w for w in warnings), warnings


def test_truncated_string_column_is_not_decoded_as_numbers():
    tab = lookup_tables(1, 3)
    res = extract(build(LOOKUP_TEXT, [record("a", 1.0)], tables=tab[:tab.index(b"veture") + 60]))
    assert [t["name"] for t in res["tables"]] == ["activite"]            # the cut table is dropped, not guessed
    assert res["tables"][0]["columns"][0]["values"][0] == "Sommeil"


# ── CS-BUG-EXTRACT-NUL ───────────────────────────────────────────────────────

def test_control_characters_are_stripped():
    assert ex.strip_control_chars("a\0b\x0c\tc\nd\x7f") == "ab\tc\nd"
    text, is_rtf = ex.decode_equations(b"a = 1\r\n\0b = 2\r\n\0")
    assert (text, is_rtf) == ("a = 1\nb = 2\n", False)
    rtf = b"{\\rtf1\\ansi a = 1\\par b = 2\\par c\\'e9 = 3\\par}\r\n\0"
    text, is_rtf = ex.decode_equations(rtf)
    assert is_rtf and text == "a = 1\nb = 2\ncé = 3\n", repr(text)


def test_no_nul_in_any_output_file():
    raw = b"{\\rtf1\\ansi a = 5\\par b = 2*a\\par c = a*b\\par}\r\n\0"
    res = extract(build("", RECORDS, raw_text=raw))
    _, files = written(res)
    assert files["m.eescode"].endswith("c = a*b\n") and all("\0" not in t for t in files.values())


def test_unicode_equations_are_decoded():
    raw = "﻿a = 5\r\nb = 2*a\r\n".encode("utf-16-le")
    res = extract(build("", RECORDS[:2], raw_text=raw))
    assert res["text"] == "a = 5\nb = 2*a\n" and len(res["variables"]) == 2


# ── CS-BUG-EXTRACT-STALE ─────────────────────────────────────────────────────

def stale(text, recs):
    return dict(extract(build(text, recs))["suspect_values"])


def test_stale_values_are_detected():
    text = "x = 5 [kPa]\ny = 2\nz = x*y + 3^2\nw = z/x\nq = 7\n"
    got = stale(text, [record("x", 6.0), record("y", 2.0), record("z", 19.0), record("w", 19 / 6),
                       record("q", math.inf)])
    assert "stored value 6, but the equation `x = 5 [kPa]` gives 5" in got["x"], got
    assert got["q"] == "stored value is inf"
    assert "z" in got and "w" not in got                                 # z is wrong, w is consistent with it


def test_consistent_solution_gives_no_warning():
    text = "x = 5 [kPa]\nr = x/2\nh = enthalpy(water, T=r, P=x)\nk = 1/0\nm = 2^0.5*x\nn = x^-1\n"
    recs = [record("x", 5.0), record("r", 2.5), record("h", 123.0), record("k", 7.0), record("m", 2 ** 0.5 * 5),
            record("n", 0.2)]
    assert extract(build(text, recs))["suspect_values"] == []
    warnings, files = written(extract(build(text, recs)))
    assert not any("STALE" in w for w in warnings)
    assert files["reference/ees_variables.csv"].startswith("# The stored values are")    # the limitation is always stated
    assert "## Stored solution" in files["reference/ees_extract_report.md"]


def test_equations_not_to_be_trusted_are_skipped():
    text = "$if opt$='a'\nx = 1\n$else\nx = 2\n$endif\ny = 3\nz[i] = 4\nt = -9999\n"
    recs = [record("x", 9.0), record("y", 3.0), record("t", -9999.0)]
    assert extract(build(text, recs))["suspect_values"] == []             # x is defined twice, t was cleared


def test_stale_warning_csv_header_and_comparison():
    res = extract(build(TEXT, [record("a", 5.0), record("b", 10.0), record("c", 7.0)]))
    warnings, files = written(res)
    assert any("STALE" in w and "c: stored value 7" in w for w in warnings), warnings
    csv_text = files["reference/ees_variables.csv"]
    assert csv_text.splitlines()[0].startswith("# ") and "# suspect: c" in csv_text
    with tempfile.TemporaryDirectory() as d:
        ref, sol = os.path.join(d, "ees_variables.csv"), os.path.join(d, "m.sol")
        Path(ref).write_text(csv_text, encoding="utf-8")
        Path(sol).write_text("a = 5\nb = 10\nc = 50\n", encoding="utf-8")
        assert cmp.read_suspect(ref) == {"c"} and set(cmp.read_reference(ref)) == {"a", "b", "c"}
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            status = cmp.main([sol, ref])
        assert status == 1 and "**DIFF** (suspect stored value)" in out.getvalue(), out.getvalue()


# ── CS-BUG-COMPARE-UNIT-CASE ─────────────────────────────────────────────────

def test_ees_units_are_matched_ignoring_the_case():
    conv = cmp.to_coolsolve_units
    assert conv("R_air", 0.28705, "kj/(kg-C)")[0] == 287.05
    assert conv("h", 1.0, "kJ/kg")[0] == 1000.0 and conv("h", 1.0, "KJ/KG")[0] == 1000.0
    assert conv("p", 5.0, "kpa")[0] == 5000.0 and conv("p", 5.0, "BAR")[0] == 5e5 and conv("p", 2.0, "mpa")[0] == 2e6
    assert conv("W", 3.0, "KW")[0] == 3000.0
    assert abs(conv("T_1", 300.0, "k")[0] - 26.85) < 1e-9
    assert conv("DELTAT_1", 5.0, "k")[0] == 5.0                          # a difference of temperatures: unchanged


def test_milli_units_and_other_units_are_not_converted():
    conv = cmp.to_coolsolve_units
    for units in ["mW", "mJ", "mPa-s", "mbar", "kgf", "", "m", "eur/kW", "W/kj"]:
        assert conv("x", 7.0, units)[0] == 7.0, units
    assert conv("x", 7.0, "kW-h")[0] == 7000.0                          # unchanged behaviour for the numerator


def test_compare_with_lower_case_units_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        Path(d, "bug.sol").write_text('R_air = 2.870500000000e+02 "J/kg-K"\nh_a = 1.000000000000e+03 "J/kg"\n')
        Path(d, "ref.csv").write_text("name,value,guess,lower,upper,units\nR_air,0.28705,1,-inf,inf,kj/(kg-C)\n"
                                      "h_a,1,1,-inf,inf,kJ/kg\n")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            status = cmp.main([os.path.join(d, "bug.sol"), os.path.join(d, "ref.csv"), "--ees-units"])
        assert status == 0 and "2 common variables, 0 differ" in out.getvalue(), out.getvalue()


# ── Negative tests: malformed and truncated files ────────────────────────────

def test_garbage_and_short_files_give_a_clear_error():
    for blob in [b"", b"X7.888", os.urandom(300), b"\x06X7.888" + bytes(10) + struct.pack("<I", 10 ** 6) + bytes(50)]:
        try:
            extract(blob)
        except ValueError as e:
            assert "EES" in str(e), str(e)
        else:
            raise AssertionError("no error for %r" % blob[:20])


def test_truncated_file_never_crashes_nor_invents_data():
    full = build(LOOKUP_TEXT, [record("a", 1.0)], tables=lookup_tables(1, 3))
    for cut in range(0, len(full), 5):
        try:
            res = extract(full[:cut])
        except ValueError:
            continue                                                    # cut inside the header or the equations
        warnings, files = written(res)
        for t in res["tables"]:                                         # whatever is extracted is complete
            assert all(len(c["values"]) == t["nrows"] for c in t["columns"]), cut
        assert [v["name"] for v in res["variables"]] in ([], ["a"]), cut
        if not res["variables"]:
            assert any("No variable records decoded" in w for w in warnings), (cut, warnings)


def test_main_reports_errors_and_continues():
    with tempfile.TemporaryDirectory() as d:
        good, bad = os.path.join(d, "good.ees"), os.path.join(d, "bad.ees")
        Path(good).write_bytes(build(TEXT, RECORDS))
        Path(bad).write_bytes(b"not an ees file at all, not an ees file at all")
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = ex.main([bad, good, "-o", os.path.join(d, "o")])
        assert status == 1 and "ERROR" in err.getvalue() and "bad.ees" in err.getvalue()
        assert "3 variables" in out.getvalue() and os.path.exists(os.path.join(d, "o", "good", "good.eescode"))


def test_unbalanced_rtf_does_not_crash():
    text, is_rtf = ex.decode_equations(b"{\\rtf1 a = 1\\par {b = 2\\par")
    assert is_rtf and "a = 1" in text


def test_compare_ignores_unreadable_rows_and_reports_missing_variables():
    with tempfile.TemporaryDirectory() as d:
        Path(d, "r.csv").write_text("# note\nname,value,guess,lower,upper,units\na,x,1,,,\nb,2,1,,,\nc\n")
        assert set(cmp.read_reference(os.path.join(d, "r.csv"))) == {"b"}
        assert cmp.read_suspect(os.path.join(d, "r.csv")) == set()


# ── Runner ───────────────────────────────────────────────────────────────────

def main():
    global passed, failed
    for name, fn in sorted((n, f) for n, f in globals().items() if n.startswith("test_") and callable(f)):
        try:
            fn()
        except (Exception, SystemExit) as e:                            # AssertionError or a crash: both fail
            failed += 1
            errors.append("%s: %s: %s" % (name, type(e).__name__, e))
            print("FAIL %s" % name)
        else:
            passed += 1
            if VERBOSE:
                print("ok   %s" % name)
    for e in errors:
        print("  " + e)
    print("%d passed, %d failed" % (passed, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
