#!/usr/bin/env python3
"""
test_security_api.py - Integration tests for the security guards of the GUI server.

Covers:
  - POST /api/v1/latex/compile runs only pdflatex, xelatex or lualatex, whether the
    name comes from the request body or from latexCompiler in coolsolve.conf, and
    accepts only plain PNG names for the plot images;
  - no CORS header is sent, and a POST, PUT or DELETE whose Origin names another
    site is refused (403), while same-origin, proxied (X-Forwarded-Host) and
    non-browser (no Origin) requests pass;
  - a request whose Host header is not localhost, an IP address or a name given with
    --allowed-host is refused (403), whatever its method (DNS rebinding);
  - GET /api/v1/stats/log shows client networks, never full client addresses;
  - the server listens on 127.0.0.1 only by default (auto mode).

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_security_api.py [base_url]

    # Or run with auto-start (from the build/ directory):
    python3 test_security_api.py --auto

Exit code 0 if all tests pass, 1 otherwise.
"""

import http.cookiejar
import json
import os
import signal
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

BASE_URL = "http://localhost:8550"
ALLOWED_HOST = "coolsolve.test"   # passed with --allowed-host in auto mode
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv

MODEL = """\
T_in = 25
P = 101325
h = enthalpy('Water', T=T_in, P=P)
"""

# The commands an attacker would try: each creates this file if it ever runs
MARKER = os.path.join(tempfile.gettempdir(), "coolsolve_security_test_marker")

# -- Helpers -------------------------------------------------------------------

passed = 0
failed = 0
errors = []

_cookie_jar = http.cookiejar.CookieJar()
_opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cookie_jar))


def request(method, path, body=None, headers=None, raw_body=None):
    """Send a request; return (status, headers, body text). Never raises on HTTP errors."""
    url = f"{BASE_URL}{path}"
    data = raw_body
    all_headers = dict(headers or {})
    if body is not None:
        data = json.dumps(body).encode()
        all_headers.setdefault("Content-Type", "application/json")
    req = urllib.request.Request(url, data=data, headers=all_headers, method=method)
    try:
        with _opener.open(req, timeout=60) as resp:
            return resp.status, resp.headers, resp.read().decode(errors="replace")
    except urllib.error.HTTPError as e:
        return e.code, e.headers, e.read().decode(errors="replace")


def check(description, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        if VERBOSE:
            print(f"  ok   {description}")
    else:
        failed += 1
        msg = f"  FAIL {description}"
        if detail:
            msg += f" - {detail}"
        print(msg)
        errors.append(description)


def host_header():
    """The Host header urllib sends for BASE_URL (host[:port])."""
    return urllib.parse.urlparse(BASE_URL).netloc


def wait_solve_done(timeout_s=60.0):
    """Consume the SSE progress stream until the solve finishes."""
    url = f"{BASE_URL}/api/v1/solve/stream"
    deadline = time.time() + timeout_s
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
                break
        except Exception:
            pass
        time.sleep(0.2)
    return None


def remove_marker():
    if os.path.exists(MARKER):
        os.remove(MARKER)


# -- Tests ---------------------------------------------------------------------

def test_latex_compiler_allow_list():
    print("\n-- LaTeX compile: compiler allow-list and plot names --")
    remove_marker()
    status, _, _ = request("PUT", "/api/v1/files/conf", {"content": ""})
    check("conf cleared", status == 200, str(status))
    status, _, _ = request("PUT", "/api/v1/files/eescode", {"content": MODEL})
    check("model stored", status == 200, str(status))
    status, _, _ = request("POST", "/api/v1/solve", {})
    check("solve started", status == 200, str(status))
    evt = wait_solve_done()
    check("solve finished", evt is not None and evt.get("type") == "done", str(evt)[:200])
    status, _, text = request("GET", "/api/v1/latex/report")
    available = status == 200 and json.loads(text).get("available") is True
    check("LaTeX report available", available, text[:200])
    if not available:
        # Every compile would answer 400 "No LaTeX report available": the checks below
        # would pass or fail for the wrong reason
        return

    not_allowed = "is not allowed"

    # A command in the request body is refused before anything runs
    injected = f"touch {MARKER} #"
    status, _, text = request("POST", "/api/v1/latex/compile", {"compiler": injected})
    check("command in body refused (400)", status == 400 and not_allowed in text, f"{status} {text[:200]}")
    check("error names the allowed compilers", "pdflatex, xelatex or lualatex" in text, text[:200])
    check("command in body not run", not os.path.exists(MARKER))

    for name in ["pdflatex; touch " + MARKER, "/usr/bin/pdflatex", "PDFLATEX", "pdflatex "]:
        status, _, text = request("POST", "/api/v1/latex/compile", {"compiler": name})
        check(f"compiler {name!r} refused", status == 400 and not_allowed in text, f"{status} {text[:200]}")
    check("no refused compiler ran", not os.path.exists(MARKER))

    # The same command through coolsolve.conf
    status, _, _ = request("PUT", "/api/v1/files/conf", {"content": f"latexCompiler = {injected}\n"})
    check("conf with command stored", status == 200, str(status))
    status, _, text = request("POST", "/api/v1/latex/compile", {})
    check("command in conf refused (400)", status == 400 and not_allowed in text, f"{status} {text[:200]}")
    check("error names coolsolve.conf", "coolsolve.conf" in text, text[:200])
    check("command in conf not run", not os.path.exists(MARKER))

    # An allowed name passes the check: the answer is the PDF, or a compile error when no
    # TeX installation is present, never the 400 of the allow-list
    status, _, _ = request("PUT", "/api/v1/files/conf", {"content": "latexCompiler = xelatex \r\n"})
    check("conf with xelatex stored", status == 200, str(status))
    status, _, text = request("POST", "/api/v1/latex/compile", {})
    check("xelatex from conf (CRLF, trailing blank) accepted", status != 400, f"{status} {text[:200]}")
    request("PUT", "/api/v1/files/conf", {"content": ""})
    status, _, text = request("POST", "/api/v1/latex/compile", {"compiler": "pdflatex"})
    check("pdflatex from body accepted", status != 400, f"{status} {text[:200]}")

    # Plot names that would overwrite the report or leave the build folder
    png = "data:image/png;base64,iVBORw0KGgo="
    for name in ["report.tex", "../x.png", "a.png.tex", "x/y.png", ".png", "plot.PNG", "p q.png"]:
        status, _, text = request("POST", "/api/v1/latex/compile",
                                  {"plots": [{"name": name, "data": png}]})
        check(f"plot name {name!r} refused", status == 400 and "Invalid plot name" in text,
              f"{status} {text[:200]}")
    status, _, text = request("POST", "/api/v1/latex/compile",
                              {"plots": [{"name": "parametric_1.png", "data": png}]})
    check("plot name 'parametric_1.png' accepted", status != 400, f"{status} {text[:200]}")
    remove_marker()


def test_cross_site_guard():
    print("\n-- Cross-site guard and CORS --")
    status, headers, _ = request("GET", "/api/v1/health", headers={"Origin": "http://evil.example"})
    check("cross-site GET answered", status == 200, str(status))
    check("no Access-Control-Allow-Origin header", headers.get("Access-Control-Allow-Origin") is None,
          str(headers.get("Access-Control-Allow-Origin")))

    status, headers, _ = request("OPTIONS", "/api/v1/files/conf", headers={
        "Origin": "http://evil.example",
        "Access-Control-Request-Method": "PUT",
        "Access-Control-Request-Headers": "content-type",
    })
    check("cross-site preflight refused", status == 403, str(status))
    check("preflight grants no origin", headers.get("Access-Control-Allow-Origin") is None)

    # A "simple" cross-site POST (text/plain, no preflight in a browser)
    status, _, _ = request("POST", "/api/v1/new", raw_body=b"{}",
                           headers={"Origin": "http://evil.example", "Content-Type": "text/plain"})
    check("cross-site POST refused", status == 403, str(status))
    status, _, _ = request("PUT", "/api/v1/model-name", {"modelName": "x"},
                           headers={"Origin": "http://evil.example"})
    check("cross-site PUT refused", status == 403, str(status))
    status, _, _ = request("DELETE", "/api/v1/tables/x", headers={"Origin": "http://evil.example"})
    check("cross-site DELETE refused", status == 403, str(status))
    status, _, _ = request("POST", "/api/v1/new", {}, headers={"Origin": "null"})
    check("Origin 'null' refused", status == 403, str(status))
    status, _, _ = request("POST", "/api/v1/new", {}, headers={"Origin": "http://" + host_header() + ".evil"})
    check("look-alike origin refused", status == 403, str(status))

    status, _, _ = request("POST", "/api/v1/new", {}, headers={"Origin": "http://" + host_header()})
    check("same-origin POST accepted", status == 200, str(status))
    status, _, _ = request("POST", "/api/v1/new", {},
                           headers={"Origin": "http://" + host_header().upper()})
    check("same-origin POST accepted (host case ignored)", status == 200, str(status))
    status, _, _ = request("POST", "/api/v1/new", {}, headers={
        "Origin": "http://localhost:5173", "X-Forwarded-Host": "localhost:5173, proxy.local"})
    check("proxied POST accepted (Vite dev, X-Forwarded-Host)", status == 200, str(status))
    status, _, _ = request("POST", "/api/v1/new", {})
    check("POST without Origin accepted (scripts)", status == 200, str(status))


def test_host_allow_list():
    print("\n-- Host header allow-list (DNS rebinding) --")
    port = urllib.parse.urlparse(BASE_URL).port
    # A rebound page: Host and Origin both name the attacker's domain
    rebound = {"Host": f"attacker.example:{port}", "Origin": f"http://attacker.example:{port}"}
    status, _, text = request("POST", "/api/v1/new", {}, headers=rebound)
    check("rebound POST refused", status == 403 and "--allowed-host" in text, f"{status} {text[:200]}")
    status, _, _ = request("GET", "/api/v1/files/eescode", headers=rebound)
    check("rebound GET refused", status == 403, str(status))
    status, _, _ = request("GET", "/", headers={"Host": f"attacker.example:{port}"})
    check("rebound page load refused", status == 403, str(status))
    status, _, _ = request("GET", "/api/v1/health", headers={"Host": "localhost.attacker.example"})
    check("look-alike host refused", status == 403, str(status))
    for host in ["[attacker.example", "[attacker.example]", "999.999.999.999", "1.2.3.4.attacker.example"]:
        status, _, _ = request("GET", "/api/v1/health", headers={"Host": host})
        check(f"Host {host!r} refused", status == 403, str(status))

    for host in [f"localhost:{port}", f"LOCALHOST:{port}", f"127.0.0.1:{port}", f"[::1]:{port}",
                 f"192.168.1.5:{port}", "localhost."]:
        status, _, _ = request("GET", "/api/v1/health", headers={"Host": host})
        check(f"Host {host!r} accepted", status == 200, str(status))
    # Given with --allowed-host in auto mode
    status, _, _ = request("POST", "/api/v1/new", {}, headers={
        "Host": ALLOWED_HOST, "Origin": f"https://{ALLOWED_HOST}"})
    check(f"Host {ALLOWED_HOST!r} from --allowed-host accepted", status == 200, str(status))


def test_stats_hide_addresses():
    print("\n-- Usage statistics: networks, not addresses --")
    status, _, text = request("GET", "/api/v1/stats/log")
    check("stats answered", status == 200, str(status))
    stats = json.loads(text)
    check("no topIps field", "topIps" not in stats)
    check("topNetworks field", isinstance(stats.get("topNetworks"), list))
    if stats.get("valid") and stats.get("totalAttempts", 0) > 0:
        names = [n["name"] for n in stats["topNetworks"]]
        check("networks listed", len(names) >= 1, str(names))
        check("every entry is a network", all(n.endswith("/24") or n.endswith("/48") or n == "(other)"
                                              for n in names), str(names))
        check("no full loopback address in the answer", "127.0.0.1" not in text)


def local_non_loopback_ip():
    """An address of this machine that is not loopback, or None (no packet is sent)."""
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("10.255.255.255", 1))
        ip = s.getsockname()[0]
        s.close()
        return None if ip.startswith("127.") else ip
    except OSError:
        return None


def test_loopback_only(port):
    print("\n-- Listening address --")
    ip = local_non_loopback_ip()
    if ip is None:
        print("  (skipped: no non-loopback address on this machine)")
        return
    try:
        with socket.create_connection((ip, port), timeout=3):
            reachable = True
    except OSError:
        reachable = False
    check(f"not reachable on {ip}:{port}", not reachable)
    with socket.create_connection(("127.0.0.1", port), timeout=3):
        pass
    check("reachable on 127.0.0.1", True)


# -- Main ----------------------------------------------------------------------

def main():
    global BASE_URL, failed
    server_proc = None
    port = 18603

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
        tmpdir = tempfile.mkdtemp(prefix="coolsolve_security_test_")
        env = dict(os.environ)
        env["COOLSOLVE_GUI_LOG"] = os.path.join(tmpdir, "gui_usage.log")
        print(f"Starting CoolSolve from {binary}...")
        server_proc = subprocess.Popen([binary, "--gui", str(port), "--no-browser",
                                        "--allowed-host", ALLOWED_HOST],
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        BASE_URL = f"http://127.0.0.1:{port}"
        for _ in range(40):
            time.sleep(0.5)
            try:
                if request("GET", "/api/v1/health")[0] == 200:
                    break
            except Exception:
                pass
        else:
            print("ERROR: Server did not start in 20 s")
            server_proc.terminate()
            sys.exit(1)

    try:
        print(f"Connected to {BASE_URL}")
        test_latex_compiler_allow_list()
        test_cross_site_guard()
        test_stats_hide_addresses()
        if auto_start:
            test_host_allow_list()
            test_loopback_only(port)
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        failed += 1
        errors.append(str(e))
    finally:
        remove_marker()
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
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
