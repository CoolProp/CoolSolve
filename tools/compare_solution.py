#!/usr/bin/env python3
"""
compare_solution.py -- compare a CoolSolve solution with an EES reference.

Usage:
    python3 tools/compare_solution.py model.sol reference/ees_variables.csv [--rtol 1e-3]
    python3 tools/compare_solution.py model.sol ref.csv --row 2     # one row of a parametric table

The reference is either the `ees_variables.csv` written by tools/ees_extract.py
(columns name,value,...) or a parametric table CSV (one column per variable,
one row per run; select the run with --row, 1-based).  Variable names are
matched case-insensitively, as in EES.

--ees-units converts the EES reference values to CoolSolve units (Pa, J, degC)
using the units stored by EES for each variable (column `units` of
ees_variables.csv): kPa/bar/MPa -> Pa, kJ -> J, kW -> W, and K -> degC for
variables whose name starts with T (not DELTAT/DT: temperature differences are
unchanged).  The unit names are matched ignoring the case (EES stores declared
units in lower case).  Use it after converting a kPa/kJ/K model
(docs/ees_import.md §6).

The stored solution of an EES file can be stale (docs/ees_import.md, "Stored
solution may be stale"): a difference on a variable listed in the `# suspect:`
line of ees_variables.csv is marked, but a reference without that mark is not
thereby correct.

Exit status: 0 if every common variable agrees within tolerance, 1 otherwise.
Absolute enthalpies/entropies may legitimately differ between EES and CoolProp
(different reference states): check that the *derived* results agree.
"""

import argparse
import csv
import math
import re
import sys


def read_sol(path):
    sol = {}
    rx = re.compile(r'^\s*([^=\s]+)\s*=\s*([-+0-9.eEinfINFnaN]+)\s*(?:"([^"]*)")?')
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = rx.match(line)
            if m:
                try:
                    sol[m.group(1).lower()] = (m.group(1), float(m.group(2)), m.group(3) or "")
                except ValueError:
                    pass
    return sol


def read_reference(path, row=None):
    with open(path, encoding="utf-8", newline="") as f:
        rows = [r for r in csv.reader(f) if not (r and r[0].startswith("#"))]     # `#` lines are notes
    ref = {}
    if rows and [h.lower() for h in rows[0][:2]] == ["name", "value"]:
        units_col = rows[0].index("units") if "units" in rows[0] else None
        for r in rows[1:]:
            try:
                ref[r[0].lower()] = (r[0], float(r[1]), r[units_col] if units_col is not None else "")
            except (ValueError, IndexError):
                pass
    else:                                     # parametric table: header = variable names
        if row is None:
            raise SystemExit("reference looks like a parametric table: select a run with --row N")
        values = rows[row]
        for name, val in zip(rows[0], values):
            try:
                ref[name.lower()] = (name, float(val), "")
            except ValueError:
                pass                          # string column (e.g. fluid$) or empty cell
    return ref


def read_suspect(path):
    """Names listed in the `# suspect:` line of an ees_variables.csv: stored values that look stale."""
    with open(path, encoding="utf-8", newline="") as f:
        for line in f:
            if line.startswith("# suspect:"):
                return {n.strip().lower() for n in line[len("# suspect:"):].split(",") if n.strip()}
            if not line.startswith("#"):
                break
    return set()


_SCALE = [("MPa", 1e6), ("kPa", 1e3), ("bar", 1e5), ("kJ", 1e3), ("kW", 1e3), ("MJ", 1e6), ("MW", 1e6)]


def to_coolsolve_units(name, value, units):
    """Convert one EES reference value to CoolSolve units; return (value, note).

    EES stores the units a variable was declared with in lower case (`kj/(kg-C)`, `k`), so the unit names are
    matched ignoring the case.  The milli spellings `mW`, `mJ`, `mPa` are not mega units and stay unconverted;
    the all-lower-case `mw`, `mj`, `mpa` cannot be told apart from them and are read as MW, MJ, MPa.
    """
    u = (units or "").strip()
    if u.lower() == "k" and name[:1] in "Tt":    # absolute temperature (DELTAT_... starts with D)
        return value - 273.15, "K->C"
    numerator = u.split("/")[0]   # a prefixed unit in a denominator (eur/kW) is kept as written by the model
    for unit, factor in _SCALE:
        m = re.search(r"(?<![A-Za-z])%s(?![A-Za-z])" % unit, numerator, flags=re.IGNORECASE)
        if m and not re.match(r"m[A-Z]", m.group(0)):
            return value * factor, "%s x%g" % (unit, factor)
    return value, ""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sol", help="CoolSolve .sol file")
    ap.add_argument("reference", help="ees_variables.csv or parametric table CSV")
    ap.add_argument("--rtol", type=float, default=1e-3, help="relative tolerance (default 1e-3)")
    ap.add_argument("--atol", type=float, default=1e-9, help="absolute tolerance (default 1e-9)")
    ap.add_argument("--row", type=int, help="run number (1-based) for a parametric-table reference")
    ap.add_argument("--all", action="store_true", help="list matching variables too")
    ap.add_argument("--ees-units", action="store_true",
                    help="convert EES reference values (kPa, bar, kJ, kW, K) to CoolSolve units first")
    args = ap.parse_args(argv)

    sol, ref = read_sol(args.sol), read_reference(args.reference, args.row)
    suspect = read_suspect(args.reference)
    if args.ees_units:
        converted = []
        for k, (name, val, units) in list(ref.items()):
            v2, note = to_coolsolve_units(name, val, units)
            if note:
                ref[k] = (name, v2, units)
                converted.append("%s (%s)" % (name, note))
        if converted:
            print("Converted EES reference values: " + ", ".join(converted) + "\n")
    common = sorted(set(sol) & set(ref))
    bad = []
    lines = ["| Variable | EES | CoolSolve | rel. diff | |", "|---|---:|---:|---:|---|"]
    for k in common:
        name, a, units = ref[k]
        b = sol[k][1]
        if math.isinf(a) or math.isnan(a):
            continue
        rel = abs(b - a) / max(abs(a), abs(b), 1e-300)
        ok = abs(b - a) <= args.atol + args.rtol * max(abs(a), abs(b))
        if not ok:
            bad.append(k)
        if args.all or not ok:
            verdict = "ok" if ok else "**DIFF**" + (" (suspect stored value)" if k in suspect else "")
            lines.append("| %s | %.6g | %.6g | %.2e | %s |" % (name, a, b, rel, verdict))
    print("\n".join(lines))
    print("\n%d common variables, %d differ (rtol=%g); only in EES: %d; only in CoolSolve: %d" % (
        len(common), len(bad), args.rtol, len(set(ref) - set(sol)), len(set(sol) - set(ref))))
    only = sorted(set(ref) - set(sol))
    if only:
        print("Only in EES reference: " + ", ".join(ref[k][0] for k in only[:40]) + (" ..." if len(only) > 40 else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
