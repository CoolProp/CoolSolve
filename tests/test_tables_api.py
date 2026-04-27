#!/usr/bin/env python3
"""
test_tables_api.py — Integration tests for the lookup-table REST API.

Tests the four /api/v1/tables routes (GET list, GET CSV, PUT, DELETE) and
verifies that lookup tables are correctly used in a solve.

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_tables_api.py [base_url]

    # Or run with auto-start (from build/ directory):
    python3 test_tables_api.py --auto

Exit code 0 if all tests pass, 1 otherwise.
"""

import json
import os
import signal
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import http.cookiejar

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv

# ── Helpers ──────────────────────────────────────────────────────────────────

passed = 0
failed = 0
errors: list[str] = []

_cookie_jar = http.cookiejar.CookieJar()
_opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cookie_jar))


def api(method: str, path: str, body=None, content_type: str = "application/json"):
    """Send a request to the API.  Returns parsed JSON or raw text."""
    url = f"{BASE_URL}/api/v1{path}"
    if body is None:
        data = None
    elif isinstance(body, (bytes, bytearray)):
        data = body
    elif content_type == "text/csv":
        data = body.encode() if isinstance(body, str) else body
    else:
        data = json.dumps(body).encode()
    headers = {"Content-Type": content_type}
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        resp = _opener.open(req, timeout=15)
        raw = resp.read()
        ct = resp.headers.get("Content-Type", "")
        if "json" in ct:
            return json.loads(raw)
        return raw.decode()
    except urllib.error.HTTPError as e:
        body_text = e.read().decode()
        raise RuntimeError(f"HTTP {e.code} on {method} {path}: {body_text}") from e


def check(description: str, condition: bool, detail: str = ""):
    global passed, failed
    if condition:
        passed += 1
        if VERBOSE:
            print(f"  \u2713 {description}")
    else:
        failed += 1
        msg = f"  \u2717 {description}"
        if detail:
            msg += f" \u2014 {detail}"
        print(msg)
        errors.append(description)


def check_approx(description: str, actual, expected, tol: float = 0.01):
    if expected == 0:
        ok = abs(actual) < tol
    else:
        ok = abs(actual - expected) / abs(expected) < tol
    check(description, ok, f"got {actual}, expected {expected}")


# ── Test groups ───────────────────────────────────────────────────────────────

def test_tables_crud():
    """Create, read, update, and delete a lookup table via REST."""
    print("\n\u2500\u2500 Lookup Tables: CRUD \u2500\u2500")

    # Start fresh session
    api("POST", "/new")

    # Initially empty
    result = api("GET", "/tables")
    check("GET /tables returns list", isinstance(result, (list, dict)))
    tables_list = result if isinstance(result, list) else result.get("tables", [])
    check("No tables initially", len(tables_list) == 0, f"got {len(tables_list)}")

    # PUT a simple table
    csv = "T,h\n100,2675.6\n150,2745.9\n200,2826.8\n"
    api("PUT", "/tables/mysteam", csv, content_type="text/csv")

    # List should now have 1 table
    result = api("GET", "/tables")
    tables_list = result if isinstance(result, list) else result.get("tables", [])
    check("1 table after PUT", len(tables_list) == 1, f"got {len(tables_list)}")
    if tables_list:
        tbl = tables_list[0]
        check("table name is 'mysteam'", tbl.get("name", "").lower() == "mysteam")
        check("3 rows", tbl.get("rows") == 3, f"got {tbl.get('rows')}")
        check("2 columns", len(tbl.get("columns", [])) == 2,
              f"got {tbl.get('columns')}")

    # GET CSV back
    csv_back = api("GET", "/tables/mysteam")
    check("GET /tables/:name returns CSV text", isinstance(csv_back, str))
    check("CSV contains header", "T" in csv_back and "h" in csv_back)
    check("CSV contains data", "2675" in csv_back)

    # PUT a second table
    csv2 = "x,y\n0,0\n1,1\n"
    api("PUT", "/tables/linear", csv2, content_type="text/csv")
    result = api("GET", "/tables")
    tables_list = result if isinstance(result, list) else result.get("tables", [])
    check("2 tables after second PUT", len(tables_list) == 2,
          f"got {len(tables_list)}")

    # DELETE first table
    api("DELETE", "/tables/mysteam")
    result = api("GET", "/tables")
    tables_list = result if isinstance(result, list) else result.get("tables", [])
    check("1 table after DELETE", len(tables_list) == 1, f"got {len(tables_list)}")
    if tables_list:
        check("remaining table is 'linear'",
              tables_list[0].get("name", "").lower() == "linear")

    # GET deleted table should fail
    try:
        api("GET", "/tables/mysteam")
        check("GET deleted table returns 404", False, "expected 404")
    except RuntimeError as e:
        check("GET deleted table returns 404", "404" in str(e), str(e))

    # Invalid table name should be rejected.  Percent-encode the path so the
    # request actually reaches the server (urllib refuses raw spaces).
    try:
        bad_path = "/tables/" + urllib.parse.quote("bad name!", safe="")
        api("PUT", bad_path, "a,b\n1,2\n", content_type="text/csv")
        check("PUT with invalid name returns 400", False, "expected 400")
    except RuntimeError as e:
        check("PUT with invalid name returns 400", "400" in str(e), str(e))


def test_tables_cleared_on_new():
    """POST /new must clear lookup tables from the previous session.

    Regression: previously the New endpoint cleared the .eescode source and
    other model state but left session.lookupTableCSVs untouched, so tables
    from a previously loaded model resurfaced as soon as the GUI refreshed
    its table list (e.g. after creating a new table).
    """
    print("\n\u2500\u2500 Lookup Tables: cleared by POST /new \u2500\u2500")

    # Start from a known empty session
    api("POST", "/new")

    # Pre-populate with two tables (mimicking what loading lookup_demo would do)
    api("PUT", "/tables/data",    "T,h\n100,2675\n200,2826\n", content_type="text/csv")
    api("PUT", "/tables/watercp", "T,Cp\n100,4.216\n",          content_type="text/csv")

    result = api("GET", "/tables")
    tables_list = result if isinstance(result, list) else result.get("tables", [])
    check("2 tables loaded before New", len(tables_list) == 2,
          f"got {len(tables_list)}")

    # Click "New" — must drop both tables
    api("POST", "/new")

    result = api("GET", "/tables")
    tables_list = result if isinstance(result, list) else result.get("tables", [])
    check("0 tables after POST /new", len(tables_list) == 0,
          f"got {len(tables_list)}: " +
          ", ".join(t.get("name", "?") for t in tables_list))

    # Creating a new table on the fresh session must NOT show stale tables
    api("PUT", "/tables/fresh", "x,y\n0,0\n1,1\n", content_type="text/csv")
    result = api("GET", "/tables")
    tables_list = result if isinstance(result, list) else result.get("tables", [])
    check("only the new table is listed", len(tables_list) == 1,
          f"got {len(tables_list)}: " +
          ", ".join(t.get("name", "?") for t in tables_list))
    if tables_list:
        check("listed table is 'fresh'",
              tables_list[0].get("name", "").lower() == "fresh")


def test_tables_solve():
    """Verify that a table loaded via PUT is used correctly in a solve."""
    print("\n\u2500\u2500 Lookup Tables: solve with INTERPOLATE \u2500\u2500")

    # Start fresh session
    api("POST", "/new")

    # Upload a linear table: y = x
    csv = "x,y\n0,0\n10,10\n20,20\n"
    api("PUT", "/tables/linear", csv, content_type="text/csv")

    # Model: interpolate at x=5, expect y=5
    code = "x_in = 5\ny_out = INTERPOLATE('linear', 'x', 'y', x_in)\n"
    api("PUT", "/files/eescode", {"content": code})

    # Trigger solve and wait
    api("POST", "/solve", {})
    time.sleep(3)

    try:
        result = api("GET", "/solve/result")
        check("solve succeeded", result.get("status") == "Success",
              f"status={result.get('status')}, err={result.get('errorMessage','')}")
        if result.get("status") == "Success":
            y_val = result.get("variables", {}).get("y_out")
            if y_val is not None:
                check_approx("INTERPOLATE('linear','x','y',5) = 5", y_val, 5.0, tol=1e-4)
            else:
                check("y_out in variables", False, "not found")
    except Exception as e:
        check("solve result accessible", False, str(e))


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    global BASE_URL
    server_proc = None

    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if args:
        BASE_URL = args[0].rstrip("/")

    auto_start = "--auto" in sys.argv
    if auto_start:
        binary = "./coolsolve"
        if not os.path.exists(binary):
            binary = "../build/coolsolve"
        if not os.path.exists(binary):
            print(f"ERROR: CoolSolve binary not found at {binary}")
            sys.exit(1)
        print(f"Starting CoolSolve from {binary}...")
        server_proc = subprocess.Popen(
            [binary, "--gui", "18600", "--no-browser"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        BASE_URL = "http://localhost:18600"
        for _ in range(20):
            time.sleep(0.5)
            try:
                api("GET", "/health")
                break
            except Exception:
                pass
        else:
            print("ERROR: Server did not start in 10 s")
            server_proc.terminate()
            sys.exit(1)

    try:
        health = api("GET", "/health")
        print(f"Connected to {BASE_URL} (status={health.get('status')})")

        test_tables_crud()
        test_tables_cleared_on_new()
        test_tables_solve()

    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        if server_proc:
            server_proc.send_signal(signal.SIGTERM)
            server_proc.wait(timeout=5)

    total = passed + failed
    print(f"\n{'=' * 60}")
    print(f"Results: {passed}/{total} passed, {failed} failed")
    if errors:
        print("Failed:")
        for e in errors:
            print(f"  - {e}")
    print(f"{'=' * 60}")
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
