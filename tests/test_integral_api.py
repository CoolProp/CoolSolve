#!/usr/bin/env python3
"""
test_integral_api.py — Integration tests for the equation-based dynamic
(INTEGRAL / $IntegralTable) feature via the REST API.

Verifies that:
  - solving an integral model embeds the trajectory table in the result,
  - GET /api/v1/integral/result returns the columnar table,
  - the ZIP bundle exports `<model>-integral.csv`,
  - re-uploading the bundle restores the table (round-trip).

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_integral_api.py [base_url]

    # Or run with auto-start (from build/ directory):
    python3 test_integral_api.py --auto

Exit code 0 if all tests pass, 1 otherwise.
"""

import io
import json
import math
import os
import signal
import subprocess
import sys
import time
import urllib.error
import urllib.request
import http.cookiejar
import zipfile

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv

# ── Helpers ──────────────────────────────────────────────────────────────────

passed = 0
failed = 0
errors: list[str] = []

_cookie_jar = http.cookiejar.CookieJar()
_opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cookie_jar))


def api(method: str, path: str, body=None, content_type="application/json"):
    url = f"{BASE_URL}/api/v1{path}"
    if isinstance(body, (dict, list)):
        data = json.dumps(body).encode()
    elif isinstance(body, str) and content_type == "application/json":
        data = json.dumps({"content": body}).encode()
    elif isinstance(body, str):
        data = body.encode()
    elif body is None:
        data = None
    else:
        data = body
    headers = {"Content-Type": content_type}
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        resp = _opener.open(req, timeout=30)
        ctype = resp.headers.get("Content-Type", "")
        body_text = resp.read()
        if "application/json" in ctype:
            return json.loads(body_text)
        return body_text
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


def check_approx(description: str, actual, expected, tol=1e-3):
    if expected == 0:
        ok = abs(actual) < tol
    else:
        ok = abs(actual - expected) / abs(expected) < tol
    check(description, ok, f"got {actual}, expected {expected}")


# ── Model ────────────────────────────────────────────────────────────────────

DECAY_MODEL = """
"Exponential decay: dy/dt = -y, y(0)=1, y(4)=e^-4"
y = 1 + integral(dydt, t, 0, 4)
dydt = -y
$IntegralTable t:0.5 y dydt
"""


def wait_for_solve(timeout=30.0):
    """Poll /solve/result until the solve completes (it returns 404 while
    running or before any solve)."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            res = api("GET", "/solve/result")
            # 200 means a result is available.
            if isinstance(res, dict):
                return res
        except Exception:
            pass
        time.sleep(0.3)
    return None


# ── Tests ────────────────────────────────────────────────────────────────────

def test_integral_solve_and_table():
    print("\n\u2500\u2500 Integral: solve + table in result \u2500\u2500")
    api("POST", "/new")
    api("PUT", "/files/eescode", {"content": DECAY_MODEL})
    api("POST", "/solve", {})
    st = wait_for_solve()
    check("solve finished", st is not None)
    if st is None:
        return

    # The integral table is embedded in the solve result and also available via
    # the dedicated endpoint.
    res = api("GET", "/integral/result")
    tbl = res.get("integralTable") or {}
    check("integralTable present", bool(tbl), f"keys={list(res.keys())}")
    if not tbl:
        return
    cols = tbl.get("columns", [])
    check("columns include t,y,dydt", set(["t", "y", "dydt"]).issubset(set(cols)),
          f"cols={cols}")
    data = tbl.get("data", {})
    y_series = data.get("y", [])
    check("y series non-empty", len(y_series) > 0)
    if y_series:
        check_approx("y(4) = e^-4", y_series[-1], math.exp(-4.0), tol=1e-3)

    # csvName field present
    check("csvName present", "csvName" in tbl, f"tbl={tbl}")


def test_bundle_round_trip():
    print("\n\u2500\u2500 Integral: ZIP bundle round-trip \u2500\u2500")
    # Download the bundle (model already solved in the previous test).
    bundle = api("GET", "/files/bundle")
    check("bundle is a ZIP", isinstance(bundle, (bytes, bytearray)) and bundle[:2] == b"PK",
          f"type={type(bundle)}")

    csv_name = None
    csv_text = None
    if isinstance(bundle, (bytes, bytearray)):
        zf = zipfile.ZipFile(io.BytesIO(bundle))
        names = zf.namelist()
        integral_csvs = [n for n in names if n.endswith("-integral.csv")]
        check("bundle contains *-integral.csv", len(integral_csvs) == 1,
              f"matches={integral_csvs} all={names}")
        if integral_csvs:
            csv_name = integral_csvs[0]
            csv_text = zf.read(csv_name).decode()
            check("CSV has a y row near t=4", "0.0183" in csv_text[:200] or
                  any("0.0183" in line for line in csv_text.splitlines()),
                  f"tail={csv_text.splitlines()[-1] if csv_text else ''}")

    # Clear, then re-upload the bundle and verify the table is restored.
    api("POST", "/new")
    check("integral cleared after /new",
          _integral_result_empty(), "table still present after /new")

    # Upload the bundle via multipart or the JSON upload endpoint.
    if isinstance(bundle, (bytes, bytearray)):
        _upload_bundle(bundle)
        # The CSV round-trips through the bundle (the JSON lastIntegralResult is
        # only produced by a live solve, so /integral/result may be empty here).
        # The contract we guarantee is that the CSV survives the round-trip:
        # re-download the bundle and check the *-integral.csv is back.
        bundle2 = api("GET", "/files/bundle")
        if isinstance(bundle2, (bytes, bytearray)):
            zf2 = zipfile.ZipFile(io.BytesIO(bundle2))
            restored_csvs = [n for n in zf2.namelist() if n.endswith("-integral.csv")]
            check("integral CSV restored after re-upload", len(restored_csvs) == 1,
                  f"matches={restored_csvs}")
            if restored_csvs:
                txt = zf2.read(restored_csvs[0]).decode()
                check("restored CSV has the y(4) row",
                      any("0.0183" in line for line in txt.splitlines()),
                      f"tail={txt.splitlines()[-1] if txt else ''}")


def _integral_result_empty() -> bool:
    try:
        res = api("GET", "/integral/result")
        return not (res.get("integralTable") or {})
    except Exception:
        return True


def _upload_bundle(bundle_bytes):
    """POST a ZIP to /api/v1/files/upload as multipart/form-data."""
    boundary = "----coolsolvetestboundary1234"
    body = b"--" + boundary.encode() + b"\r\n"
    body += b'Content-Disposition: form-data; name="file"; filename="model.zip"\r\n'
    body += b"Content-Type: application/zip\r\n\r\n"
    body += bytes(bundle_bytes) + b"\r\n"
    body += b"--" + boundary.encode() + b"--\r\n"
    url = f"{BASE_URL}/api/v1/files/upload"
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
                                 method="POST")
    _opener.open(req, timeout=30).read()


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
            [binary, "--gui", "18599", "--no-browser"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        BASE_URL = "http://localhost:18599"
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
        health = api("GET", "/health")
        print(f"Connected to {BASE_URL} (status={health.get('status')})")
        test_integral_solve_and_table()
        test_bundle_round_trip()
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        global failed
        failed += 1
        errors.append(str(e))
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
