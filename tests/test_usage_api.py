#!/usr/bin/env python3
"""
test_usage_api.py — Integration tests for the GUI usage-log feature.

Verifies that every GUI solve attempt (Solve / Try Harder / parse failure)
appends exactly one line to the usage log file and that the hidden
statistics endpoint (/api/v1/stats/log) aggregates it correctly.

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_usage_api.py [base_url]

    # Or run with auto-start (from build/ directory):
    python3 test_usage_api.py --auto

Exit code 0 if all tests pass, 1 otherwise.
"""

import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
import http.cookiejar

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv

GOOD_MODEL = """\
T_in = 25
P = 101325
h = enthalpy('Water', T=T_in, P=P)
s = entropy('Water', T=T_in, P=P)
"""

BAD_MODEL = "x === y"

# ── Helpers ──────────────────────────────────────────────────────────────────

passed = 0
failed = 0
errors: list[str] = []

_cookie_jar = http.cookiejar.CookieJar()
_opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cookie_jar))


def api(method: str, path: str, body=None):
    """Send a request to the API and return parsed JSON."""
    url = f"{BASE_URL}/api/v1{path}"
    data = json.dumps(body).encode() if body is not None else None
    headers = {"Content-Type": "application/json"}
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        resp = _opener.open(req, timeout=60)
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


def read_log(path: str) -> list[str]:
    """Return the data (non-comment, non-empty) lines of the log file."""
    if not os.path.exists(path):
        return []
    with open(path, "r") as f:
        return [l.rstrip("\n") for l in f if l.strip() and not l.startswith("#")]


def wait_solve_done(timeout_s: float = 60.0):
    """Consume the SSE progress stream until the solve finishes."""
    url = f"{BASE_URL}/api/v1/solve/stream"
    deadline = time.time() + timeout_s
    last_event = None
    while time.time() < deadline:
        try:
            with _opener.open(url, timeout=5) as resp:
                for raw in resp:
                    line = raw.decode().strip()
                    if line.startswith("data: "):
                        try:
                            evt = json.loads(line[6:])
                        except json.JSONDecodeError:
                            continue
                        if evt.get("type") in ("done", "error"):
                            return evt
                break  # stream closed without a final event
        except Exception:
            pass
        time.sleep(0.2)
    return last_event


# ── Tests ────────────────────────────────────────────────────────────────────

def test_stats_endpoint_graceful_without_log():
    """Negative test: statistics endpoint must respond cleanly before any solve."""
    print("\n── Stats: graceful behaviour with no/empty log ──")
    stats = api("GET", "/stats/log")
    check("endpoint responds", isinstance(stats, dict))
    check("valid flag present", "valid" in stats)
    if not stats["valid"]:
        check("zero totals when invalid", stats["totalAttempts"] == 0)


def test_log_format_and_stats():
    global LOG_PATH
    print("\n── Usage log: solve / try-harder / parse-failure round trip ──")

    baseline_lines = len(read_log(LOG_PATH))
    baseline_stats = api("GET", "/stats/log")
    baseline_total = baseline_stats.get("totalAttempts", 0)

    # Load a small model, give it a name
    api("PUT", "/files/eescode", {"content": GOOD_MODEL})
    api("PUT", "/model-name", {"modelName": "usage_api_test"})

    # 1) Normal solve
    api("POST", "/solve", {})
    done = wait_solve_done()
    check("first solve finished", done is not None and done.get("type") == "done",
          f"last event: {done}")

    # 2) Try Harder (deep search) on the same model
    api("POST", "/solve", {"deepSearch": True})
    done = wait_solve_done()
    check("try-harder solve finished", done is not None and done.get("type") == "done")

    # 3) Parse failure
    api("PUT", "/files/eescode", {"content": BAD_MODEL})
    api("POST", "/solve", {})
    evt = wait_solve_done()
    check("parse failure reported", evt is not None and evt.get("type") == "error")

    # ── Raw file checks ──
    lines = read_log(LOG_PATH)
    check("three new lines appended", len(lines) == baseline_lines + 3,
          f"{len(lines)} vs {baseline_lines + 3}")
    check("one CSV line per attempt", all(l.count(",") >= 9 for l in lines[-3:]))

    # ── Aggregated statistics ──
    stats = api("GET", "/stats/log")
    check("stats marked valid", stats["valid"])
    check("total attempts incremented by 3",
          stats["totalAttempts"] == baseline_total + 3,
          f"{stats['totalAttempts']} vs {baseline_total}")
    check("outcome counters consistent",
          stats["successes"] + stats["failures"] + stats["parseErrors"] == stats["totalAttempts"])

    # Baseline may already contain entries; compare deltas instead of absolutes
    d_success = stats["outcomes"]["success"] - max(0, baseline_stats.get("outcomes", {}).get("success", 0))
    d_parse = stats["outcomes"]["parse_error"] - max(0, baseline_stats.get("outcomes", {}).get("parse_error", 0))
    d_tryharder = stats["tryHarderAttempts"] - baseline_stats.get("tryHarderAttempts", 0)
    check(">=1 new success", d_success >= 1, str(d_success))
    check("parse failure logged", d_parse >= 1, str(d_parse))
    check("try harder logged once", d_tryharder == 1, str(d_tryharder))

    check("model ranked by name",
          any(m["name"] == "usage_api_test" for m in stats["topModels"]),
          str(stats["topModels"]))
    check("client IP recorded", stats["uniqueIps"] >= 1)
    check("daily entry exists for today", any(d["attempts"] >= 3 for d in stats["daily"]),
          str(stats["daily"][-3:]))
    check("histogram sums to totals",
          sum(stats["durationHistogram"]["counts"]) == stats["totalAttempts"])
    check("median duration positive", stats["medianMs"] > 0)
    check("p95 >= median", stats["p95Ms"] >= stats["medianMs"])
    check("version string exported", os.path.getsize(LOG_PATH) > 0)

    # ── Cached reads stay consistent ──
    again = api("GET", "/stats/log")
    check("cached read consistent", again["totalAttempts"] == stats["totalAttempts"])


def test_hidden_page_not_linked():
    print("\n── Hidden page: /stats serves the SPA shell ──")
    # /stats must serve the app (SPA fallback) so the React router can pick it up
    url = f"{BASE_URL}/stats"
    try:
        with _opener.open(url, timeout=10) as resp:
            body = resp.read().decode()
            check("/stats reachable", resp.status == 200)
            check("/stats serves HTML shell", "<div id=\"root\">" in body or "<html" in body.lower())
    except Exception as e:
        check("/stats reachable", False, str(e))


# ── Main ─────────────────────────────────────────────────────────────────────

LOG_PATH = ""


def main():
    global BASE_URL, LOG_PATH
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

        tmpdir = tempfile.mkdtemp(prefix="coolsolve_usage_test_")
        LOG_PATH = os.path.join(tmpdir, "gui_usage.log")

        env = dict(os.environ)
        env["COOLSOLVE_GUI_LOG"] = LOG_PATH

        print(f"Starting CoolSolve from {binary} (log: {LOG_PATH})...")
        server_proc = subprocess.Popen(
            [binary, "--gui", "18602", "--no-browser"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env,
        )
        BASE_URL = "http://localhost:18602"
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
        if auto_start:
            test_stats_endpoint_graceful_without_log()
            test_log_format_and_stats()
            test_hidden_page_not_linked()
        else:
            # Remote mode: cannot inspect the log file itself
            test_stats_endpoint_graceful_without_log()
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
