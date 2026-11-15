#!/usr/bin/env python3
"""
test_bundle_api.py — Opening a ZIP bundle (POST /api/v1/files/upload) through the REST API.

Verifies that:
  - the ZIP files of the usual zip tools open (tests/fixtures/zip/, tools/zip_fixtures.py: Info-ZIP
    Deflate and streamed, 7-Zip, bsdtar, ZIP64, the ZIP of a folder, macOS' Compress, Windows' code
    page 437 names) and their files are returned as they were (the GUI reads them back as JSON),
  - a ZIP that cannot be opened is refused with a message that says why (encrypted, BZIP2/LZMA, not
    a ZIP, damaged, a name outside the folder, no .eescode file, a file that is not UTF-8), and the
    model of the session is left as it was,
  - Save then Open gives the same files back,
  - an exception that a route does not catch reaches the GUI as a JSON message, not a bare
    "Internal Server Error" (a Latin-1 model given on the command line).

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_bundle_api.py [base_url]

    # Or run with auto-start (from build/ directory):
    python3 test_bundle_api.py --auto

Exit code 0 if all tests pass, 1 otherwise.
"""

import io
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
import zipfile

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv
FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "zip")

# ── Helpers ──────────────────────────────────────────────────────────────────

passed = 0
failed = 0
errors: list[str] = []

_cookie_jar = http.cookiejar.CookieJar()
_opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cookie_jar))


def request(method: str, path: str, data=None, headers=None):
    """(status, body bytes, Content-Type) of a request to the API; HTTP errors are returned, not raised."""
    req = urllib.request.Request(f"{BASE_URL}/api/v1{path}", data=data, headers=headers or {}, method=method)
    try:
        resp = _opener.open(req, timeout=30)
        return resp.status, resp.read(), resp.headers.get("Content-Type", "")
    except urllib.error.HTTPError as e:
        return e.code, e.read(), e.headers.get("Content-Type", "")


def api(method: str, path: str, body=None):
    """(status, JSON) of a JSON request; the JSON is {} when the body is not JSON."""
    data = json.dumps(body).encode() if body is not None else None
    status, content, _ = request(method, path, data, {"Content-Type": "application/json"})
    try:
        return status, json.loads(content)
    except ValueError:
        return status, {}


def upload(bundle: bytes, filename: str):
    """POST a ZIP to /api/v1/files/upload (multipart, as the Open button); return (status, JSON)."""
    boundary = "----coolsolvebundletest"
    body = (f'--{boundary}\r\nContent-Disposition: form-data; name="file"; filename="{filename}"\r\n'
            "Content-Type: application/zip\r\n\r\n").encode("utf-8") + bundle + f"\r\n--{boundary}--\r\n".encode()
    status, content, _ = request("POST", "/files/upload", body,
                                 {"Content-Type": f"multipart/form-data; boundary={boundary}"})
    try:
        return status, json.loads(content)
    except ValueError:
        return status, {"raw": content[:200]}


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


def fixture(name: str) -> bytes:
    with open(os.path.join(FIXTURES, name), "rb") as f:
        return f.read()


def make_zip(files: dict, method=zipfile.ZIP_DEFLATED) -> bytes:
    """A ZIP of {name: str or bytes} written by Python's zipfile."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", method) as zf:
        for name, content in files.items():
            zf.writestr(name, content)
    return buf.getvalue()


def session_file(kind: str) -> str:
    """The content of the session's .eescode / .initials / .sol (GET /api/v1/files/<kind>)."""
    status, res = api("GET", f"/files/{kind}")
    return res.get("content", f"<HTTP {status}: {res}>")


# The model of every fixture (tools/zip_fixtures.py)
FIXTURE_MODEL = zipfile.ZipFile(os.path.join(FIXTURES, "infozip_deflate.zip")).read("model.eescode").decode()
FIXTURE_INITIALS = zipfile.ZipFile(os.path.join(FIXTURES, "infozip_deflate.zip")).read("model.initials").decode()

# ── Tests ────────────────────────────────────────────────────────────────────


def test_zip_tools():
    print("\n── ZIP files of the usual zip tools ──")
    for name in ["infozip_deflate.zip", "infozip_stream.zip", "7zip_deflate.zip", "bsdtar.zip",
                 "python_zip64.zip", "infozip_folder.zip", "macos_folder.zip"]:
        status, res = upload(fixture(name), name)
        check(f"{name} opens", status == 200 and res.get("success"), str(res))
        check(f"{name}: the model is model.eescode", res.get("mainFile") == "model.eescode", str(res))
        check(f"{name}: .eescode read back as it was", session_file("eescode") == FIXTURE_MODEL)
        check(f"{name}: .initials read back as it was", session_file("initials") == FIXTURE_INITIALS)

    status, res = upload(fixture("windows_cp437.zip"), "windows_cp437.zip")
    check("Windows' code page 437 names become UTF-8",
          status == 200 and res.get("mainFile") == "Échangeur.eescode", str(res))
    check("windows_cp437.zip: .eescode read back as it was", session_file("eescode") == FIXTURE_MODEL)

    model = "x = 2\ny = x^2 \"in °C\"\n"
    status, res = upload(make_zip({"upper.eescode": model}), "UPPER.ZIP")
    check("a .ZIP file name in capitals opens", status == 200 and res.get("success"), str(res))
    check("UPPER.ZIP: .eescode read back as it was", session_file("eescode") == model)


def test_refused():
    print("\n── ZIP files that cannot be opened: the message says why ──")
    reference = "x = 3\n"
    status, _ = upload(make_zip({"reference.eescode": reference}), "reference.zip")
    check("the reference model opens", status == 200)

    truncated = fixture("infozip_deflate.zip")[:300]
    cases = [
        ("password-protected", fixture("infozip_encrypted.zip"), "secret.zip", "is password-protected"),
        ("BZIP2 (7-Zip)", fixture("7zip_bzip2.zip"), "bzip2.zip", "compressed with BZIP2"),
        ("LZMA", fixture("python_lzma.zip"), "lzma.zip", "compressed with LZMA"),
        ("not a ZIP archive", b"x = 2\n" * 20, "model.zip", "not a ZIP archive"),
        ("a 7z archive", b"7z\xbc\xaf\x27\x1c" + bytes(64), "model.zip", "is a 7z archive"),
        ("damaged (truncated)", truncated, "truncated.zip", "Cannot open truncated.zip"),
        ("a name outside the folder", make_zip({"../evil.eescode": "x=1"}), "evil.zip", "outside its folder"),
        ("no .eescode file", make_zip({"notes.txt": "hello"}), "notes.zip",
         "it contains no .eescode file (files: notes.txt)"),
        ("an empty archive", make_zip({}), "empty.zip", "(the archive is empty)"),
        ("a model saved by EES (Windows ANSI)", make_zip({"ees.eescode": "x = 2\nT = 20 \"\xb0C\"\n".encode("latin-1")}),
         "ees.zip", "'ees.eescode' is not UTF-8 text (byte 0xB0 on line 2)"),
        ("not UTF-8 lookup table", make_zip({"m.eescode": "x=1", "m-t.csv": "a;b\n\xe9;1\n".encode("latin-1")}),
         "m.zip", "'m-t.csv' is not UTF-8 text"),
        ("not a .zip file name", make_zip({"m.eescode": "x=1"}), "model.rar", "Please upload a .zip file"),
    ]
    for label, data, filename, expected in cases:
        status, res = upload(data, filename)
        message = res.get("error", "")
        check(f"{label} → 400 with a message", status == 400 and expected in message, f"{status} {res}")
        check(f"{label}: the session keeps its model", session_file("eescode") == reference)
    if VERBOSE:
        print("  e.g. " + upload(fixture("7zip_bzip2.zip"), "bzip2.zip")[1].get("error", ""))


def test_round_trip():
    print("\n── Save then Open ──")
    status, _ = upload(fixture("infozip_deflate.zip"), "model.zip")
    check("a Deflate ZIP opens", status == 200)
    status, bundle, ctype = request("GET", "/files/bundle")
    check("Save gives a ZIP", status == 200 and "zip" in ctype)
    names = zipfile.ZipFile(io.BytesIO(bundle)).namelist()
    check("the bundle has the model and its initials", "model.eescode" in names and "model.initials" in names, str(names))
    api("POST", "/new")
    status, res = upload(bundle, "model.zip")
    check("the bundle opens again", status == 200 and res.get("success"), str(res))
    check("the same .eescode", session_file("eescode") == FIXTURE_MODEL)
    check("the same .initials", session_file("initials") == FIXTURE_INITIALS)


def test_uncaught_exception(binary: str):
    """A Latin-1 model opened from the command line: the route that returns it cannot write it as JSON."""
    global BASE_URL
    print("\n── An uncaught exception reaches the GUI as a JSON message ──")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "latin1.eescode")
        with open(path, "wb") as f:
            f.write("x = 2\nT = 20 \"\xb0C\"\n".encode("latin-1"))
        main_url, BASE_URL = BASE_URL, "http://localhost:18631"
        proc = start_server(binary, 18631, [path])
        try:
            status, res = api("POST", "/files/open", {"path": path})
            check("the command-line file opens", status == 200, f"{status} {res}")
            status, content, ctype = request("GET", "/files/eescode")
            message = json.loads(content).get("error", "") if "application/json" in ctype else ""
            check("its content cannot be returned → 500 with a JSON message",
                  status == 500 and "Internal error (GET /api/v1/files/eescode)" in message and "UTF-8" in message,
                  f"{status} {ctype} {content[:200]!r}")
        finally:
            stop_server(proc)
            BASE_URL = main_url


# ── Runner ───────────────────────────────────────────────────────────────────

def start_server(binary: str, port: int, extra_args=()):
    """Start `coolsolve --gui <port>` (BASE_URL must already point to that port); wait until it answers."""
    print(f"Starting CoolSolve from {binary} on port {port}...")
    proc = subprocess.Popen([binary, "--gui", str(port), "--no-browser", *extra_args],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(40):
        time.sleep(0.25)
        if api("GET", "/health")[0] == 200:
            return proc
    print("ERROR: Server did not start in 10s")
    proc.terminate()
    sys.exit(1)


def stop_server(proc):
    proc.send_signal(signal.SIGTERM)
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()


def main():
    global BASE_URL, failed
    server_proc = None

    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if args:
        BASE_URL = args[0].rstrip("/")

    binary = "./coolsolve"
    if not os.path.exists(binary):
        binary = "../build/coolsolve"
    if "--auto" in sys.argv:
        if not os.path.exists(binary):
            print(f"ERROR: CoolSolve binary not found at {binary}")
            sys.exit(1)
        BASE_URL = "http://localhost:18630"
        server_proc = start_server(binary, 18630)

    try:
        status, health = api("GET", "/health")
        print(f"Connected to {BASE_URL} (status={health.get('status')})")
        test_zip_tools()
        test_refused()
        test_round_trip()
        if os.path.exists(binary):
            test_uncaught_exception(binary)
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        failed += 1
        errors.append(str(e))
    finally:
        if server_proc:
            stop_server(server_proc)

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
