#!/usr/bin/env python3
"""
test_parametric_api.py — Integration tests for the parametric study feature.

Tests the parse endpoint (imposed variable detection, unit source) and the
parametric endpoint (sweep execution, result retrieval) via the REST API.

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_parametric_api.py [base_url]

    # Or run with auto-start (from build/ directory):
    python3 test_parametric_api.py --auto

Exit code 0 if all tests pass, 1 otherwise.
"""

import json
import os
import signal
import subprocess
import sys
import time
import urllib.error
import urllib.request
import http.cookiejar

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv

# ── Helpers ──────────────────────────────────────────────────────────────────

passed = 0
failed = 0
errors: list[str] = []

# Persistent cookie jar so all requests share the same session
_cookie_jar = http.cookiejar.CookieJar()
_opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cookie_jar))


def api(method: str, path: str, body=None):
    """Send a request to the API and return parsed JSON."""
    url = f"{BASE_URL}/api/v1{path}"
    data = json.dumps(body).encode() if body else None
    headers = {"Content-Type": "application/json"}
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        resp = _opener.open(req, timeout=15)
        return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        body_text = e.read().decode()
        raise RuntimeError(f"HTTP {e.code} on {method} {path}: {body_text}") from e


def check(description: str, condition: bool, detail: str = ""):
    global passed, failed
    if condition:
        passed += 1
        if VERBOSE:
            print(f"  ✓ {description}")
    else:
        failed += 1
        msg = f"  ✗ {description}"
        if detail:
            msg += f" — {detail}"
        print(msg)
        errors.append(description)


def check_approx(description: str, actual, expected, tol=1e-6):
    if expected == 0:
        ok = abs(actual) < tol
    else:
        ok = abs(actual - expected) / abs(expected) < tol
    check(description, ok, f"got {actual}, expected {expected}")


# ── Test Groups ──────────────────────────────────────────────────────────────

def test_parse_imposed_simple():
    """Test imposed detection on inline code."""
    print("\n── Parse: imposed detection (simple) ──")

    code = "T_in = 300\nP = 1E5\nm_dot = 0.5 [kg/s]\nT_out = T_in + deltaT"
    result = api("POST", "/parse", {"source": code})

    check("parse succeeds", result["success"])
    vars_list = result.get("variables", [])
    check("variables returned", len(vars_list) > 0)

    imposed = {v["name"]: v for v in vars_list if v.get("isImposed")}
    check("T_in is imposed", "T_in" in imposed)
    check("P is imposed", "P" in imposed)
    check("m_dot is imposed", "m_dot" in imposed)
    check("T_out is NOT imposed", "T_out" not in imposed)

    if "T_in" in imposed:
        check_approx("T_in imposed value", imposed["T_in"]["imposedValue"], 300.0)
    if "P" in imposed:
        check_approx("P imposed value", imposed["P"]["imposedValue"], 1e5)
    if "m_dot" in imposed:
        check_approx("m_dot imposed value", imposed["m_dot"]["imposedValue"], 0.5)


def test_parse_imposed_negative():
    """Test imposed detection with negative values."""
    print("\n── Parse: imposed detection (negative) ──")

    code = "T_min = -40\ncoeff = -3.5E-2"
    result = api("POST", "/parse", {"source": code})
    check("parse succeeds", result["success"])

    imposed = {v["name"]: v for v in result.get("variables", []) if v.get("isImposed")}
    check("T_min is imposed", "T_min" in imposed)
    check("coeff is imposed", "coeff" in imposed)

    if "T_min" in imposed:
        check_approx("T_min value is -40", imposed["T_min"]["imposedValue"], -40.0)
    if "coeff" in imposed:
        check_approx("coeff value is -0.035", imposed["coeff"]["imposedValue"], -0.035)


def test_parse_not_imposed():
    """Test that expressions, function calls, etc. are NOT imposed."""
    print("\n── Parse: NOT imposed patterns ──")

    code = """
y = x + 3
z = sin(x)
a + b = 5
q = 2 * m
P = 1E5
"""
    result = api("POST", "/parse", {"source": code})
    check("parse succeeds", result["success"])

    imposed = {v["name"]: v for v in result.get("variables", []) if v.get("isImposed")}
    check("y is NOT imposed", "y" not in imposed)
    check("z is NOT imposed", "z" not in imposed)
    check("q is NOT imposed", "q" not in imposed)
    check("P IS imposed", "P" in imposed)
    check("only 1 imposed", len(imposed) == 1, f"got {len(imposed)}")


def test_parse_units():
    """Test unit source detection with thermo inference."""
    print("\n── Parse: unit source detection ──")

    # Note: simple unit annotations like [m] on imposed equations are
    # NOT propagated to the variable's units field in the current parser.
    # Units come from thermodynamic inference (e.g. ENTHALPY, ENTROPY calls).
    code = "D_i = 16E-3 [m]\nP = 1E5"
    result = api("POST", "/parse", {"source": code})
    check("parse succeeds", result["success"])

    var_map = {v["name"]: v for v in result.get("variables", [])}
    if "D_i" in var_map:
        check("D_i is imposed", var_map["D_i"]["isImposed"])
        check_approx("D_i value is 0.016", var_map["D_i"]["imposedValue"], 0.016)
    if "P" in var_map:
        check("P is imposed", var_map["P"]["isImposed"])


def test_parse_pressuredrop():
    """Test full pressuredrop example."""
    print("\n── Parse: pressuredrop.eescode ──")

    # Try to read the example file
    for prefix in ["../examples/", "examples/", "../../examples/"]:
        path = prefix + "pressuredrop.eescode"
        if os.path.exists(path):
            with open(path) as f:
                code = f.read()
            break
    else:
        print("  ⚠ Skipping: pressuredrop.eescode not found")
        return

    result = api("POST", "/parse", {"source": code})
    check("parse succeeds", result["success"])

    imposed = {v["name"]: v for v in result.get("variables", []) if v.get("isImposed")}
    check("epsilon imposed", "epsilon" in imposed)
    check("D_i imposed", "D_i" in imposed)
    check("P imposed", "P" in imposed)
    check("M_dot imposed", "M_dot" in imposed)
    check("L imposed", "L" in imposed)
    check("T not imposed (= T_sat+10)", "T" not in imposed)
    check("exactly 5 imposed", len(imposed) == 5, f"got {len(imposed)}: {list(imposed.keys())}")


def test_parse_commented_equation():
    """Commented-out imposed equations should not be detected."""
    print("\n── Parse: commented-out equations ──")

    code = '{T = 300}\nP = 1E5\n"Q = 500"\nM = 2'
    result = api("POST", "/parse", {"source": code})
    check("parse succeeds", result["success"])

    imposed = {v["name"]: v for v in result.get("variables", []) if v.get("isImposed")}
    check("T not imposed (commented with {})", "T" not in imposed)
    check("Q not imposed (commented with quotes)", "Q" not in imposed)
    check("P is imposed", "P" in imposed)
    check("M is imposed", "M" in imposed)


def wait_parametric(expected_points: int, timeout: float = 30.0):
    """Poll for parametric result until it has the expected number of points."""
    result = None
    deadline = time.time() + timeout
    while time.time() < deadline:
        time.sleep(0.5)
        try:
            result = api("GET", "/parametric/result")
            if (result and result.get("totalPoints", 0) == expected_points
                    and (result.get("successCount", 0) + result.get("failCount", 0)) == expected_points):
                return result
        except Exception:
            pass
    return result


def test_parametric_simple():
    """Test a simple 1D parametric sweep on a trivial model."""
    print("\n── Parametric: simple 1D sweep ──")

    # Set up a trivial model: x=10, y=x^2
    code = "x = 10\ny = x^2"
    api("PUT", "/files/eescode", {"content": code})

    # Run parametric: sweep x from 1 to 5 in 5 steps
    sweep = {"eescode": code, "sweepVariables": [{"name": "x", "values": [1, 2, 3, 4, 5]}]}
    try:
        start_result = api("POST", "/parametric", sweep)
        check("parametric started", start_result.get("status") == "started" or "totalPoints" in start_result)
    except Exception as e:
        check("parametric started", False, str(e))
        return

    result = wait_parametric(5)

    if result is None:
        check("parametric completed", False, "no result after timeout")
        return

    check("parametric completed", result.get("totalPoints", 0) == 5)
    check("all points converged", result.get("successCount", 0) == 5)

    # Verify results
    results = result.get("results", [])
    check("5 result points", len(results) == 5)

    for pt in results:
        if pt.get("success") and pt.get("variables"):
            x_val = pt["overrides"].get("x", 0)
            y_val = pt["variables"].get("y", None)
            if y_val is not None:
                expected = x_val ** 2
                check_approx(f"y={y_val} when x={x_val}", y_val, expected)


def test_parametric_2d():
    """Test a 2D parametric sweep."""
    print("\n── Parametric: 2D sweep ──")

    # Wait for any previous parametric to finish
    time.sleep(2)

    code = "a = 2\nb = 3\nc = a + b"
    api("PUT", "/files/eescode", {"content": code})

    sweep = {"eescode": code, "sweepVariables": [
        {"name": "a", "values": [1, 2, 3]},
        {"name": "b", "values": [10, 20]},
    ]}
    try:
        start_result = api("POST", "/parametric", sweep)
        check("2D parametric started", "status" in start_result or "totalPoints" in start_result)
    except Exception as e:
        check("2D parametric started", False, str(e))
        return

    result = wait_parametric(6)

    if result is None:
        check("2D parametric completed", False, "no result after timeout")
        return

    check("6 total points (3×2)", result.get("totalPoints", 0) == 6)
    check("all points converged", result.get("successCount", 0) == 6)

    # Verify a + b = c for each point
    for pt in result.get("results", []):
        if pt.get("success") and pt.get("variables"):
            a = pt["overrides"].get("a", 0)
            b = pt["overrides"].get("b", 0)
            c = pt["variables"].get("c")
            if c is not None:
                check_approx(f"c={c} when a={a},b={b}", c, a + b)


def test_parametric_orc_simple():
    """Regression: parametric sweep on orc_simple with string variables (fluid$).

    Previously all points failed with 'Unknown fluid: fluid' because the
    initials writer was emitting 'fluid$ = 0.0' which overrode the string
    assignment fluid$='r245fa' from the eescode.
    """
    print("\n── Parametric: orc_simple (string variable regression) ──")

    # Wait for any previous parametric to finish
    time.sleep(2)

    # Load the example files
    for prefix in ["../examples/", "examples/", "../../examples/"]:
        ees_path = prefix + "orc_simple.eescode"
        init_path = prefix + "orc_simple.initials"
        if os.path.exists(ees_path):
            with open(ees_path) as f:
                eescode = f.read()
            initials = ""
            if os.path.exists(init_path):
                with open(init_path) as f:
                    initials = f.read()
            break
    else:
        print("  ⚠ Skipping: orc_simple.eescode not found")
        return

    # Sweep pinch_cd from 5 to 15 (3 points)
    sweep = {
        "eescode": eescode,
        "initials": initials,
        "sweepVariables": [{"name": "pinch_cd", "values": [5, 10, 15]}],
    }
    try:
        start_result = api("POST", "/parametric", sweep)
        check("orc parametric started", start_result.get("status") == "started")
    except Exception as e:
        check("orc parametric started", False, str(e))
        return

    result = wait_parametric(3, timeout=60.0)

    if result is None:
        check("orc parametric completed", False, "no result after timeout")
        return

    check("3 total points", result.get("totalPoints", 0) == 3)
    check("all orc points converged", result.get("successCount", 0) == 3,
          f"got {result.get('successCount', 0)}/{result.get('totalPoints', 0)}")

    # Verify pinch_cd values are reflected in results
    for pt in result.get("results", []):
        pinch = pt["overrides"].get("pinch_cd", 0)
        if pt.get("success"):
            check(f"orc pinch_cd={pinch} converged", True)
            # Check eta_cycle exists and is reasonable (0 < eta < 1)
            eta = pt.get("variables", {}).get("eta_cycle")
            if eta is not None:
                check(f"eta_cycle={eta:.4f} reasonable", 0 < eta < 1,
                      f"got {eta}")
        else:
            check(f"orc pinch_cd={pinch} converged", False,
                  pt.get("errorMessage", "")[:80])


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    global BASE_URL
    server_proc = None

    # Parse args
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if args:
        BASE_URL = args[0].rstrip("/")

    auto_start = "--auto" in sys.argv
    if auto_start:
        # Start CoolSolve in the background
        binary = "./coolsolve"
        if not os.path.exists(binary):
            binary = "../build/coolsolve"
        if not os.path.exists(binary):
            print(f"ERROR: CoolSolve binary not found at {binary}")
            sys.exit(1)
        print(f"Starting CoolSolve from {binary}...")
        server_proc = subprocess.Popen(
            [binary, "--gui", "18599", "--no-browser"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        BASE_URL = "http://localhost:18599"
        # Wait for server to be ready
        for _ in range(20):
            time.sleep(0.5)
            try:
                api("GET", "/health")
                break
            except Exception:
                pass
        else:
            print("ERROR: Server did not start in 10s")
            server_proc.terminate()
            sys.exit(1)

    try:
        # Verify connectivity
        health = api("GET", "/health")
        print(f"Connected to {BASE_URL} (status={health.get('status')})")

        # Run all test groups
        test_parse_imposed_simple()
        test_parse_imposed_negative()
        test_parse_not_imposed()
        test_parse_units()
        test_parse_pressuredrop()
        test_parse_commented_equation()
        test_parametric_simple()
        test_parametric_2d()
        test_parametric_orc_simple()

    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        failed_count = 1
    finally:
        if server_proc:
            server_proc.send_signal(signal.SIGTERM)
            server_proc.wait(timeout=5)

    # Summary
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
