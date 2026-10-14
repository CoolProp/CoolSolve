#!/usr/bin/env python3
"""
test_library_api.py — Integration tests for the embedded model library (GUI server).

Tests the Library menu (verified models only), opening a model folder into a session,
the per-tab session id, and the /library/ pages and files, via the REST API.

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_library_api.py [base_url]

    # Or run with auto-start (from the build/ directory):
    python3 test_library_api.py --auto

Needs a build with the library embedded (COOLSOLVE_EMBED_LIBRARY=ON, the default).
Exit code 0 if all tests pass, 1 otherwise.
"""

import http.cookiejar
import json
import os
import signal
import subprocess
import sys
import time
import urllib.error
import urllib.request

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv

passed = 0
failed = 0
errors: list[str] = []

# The default client: one cookie session, like a browser tab that did not opt into its own instance
_cookie_jar = http.cookiejar.CookieJar()
_opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cookie_jar))
# A client without cookies: only an explicit session id identifies it
_plain_opener = urllib.request.build_opener()


def raw(method: str, path: str, body=None, headers=None, opener=None):
    """Send a request; return (status, content_type, body bytes). HTTP errors are returned, not raised."""
    data = json.dumps(body).encode() if body is not None else None
    hdrs = {"Content-Type": "application/json"}
    hdrs.update(headers or {})
    req = urllib.request.Request(BASE_URL + path, data=data, headers=hdrs, method=method)
    try:
        resp = (opener or _opener).open(req, timeout=15)
        return resp.status, resp.headers.get("Content-Type", ""), resp.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers.get("Content-Type", ""), e.read()


def api(method: str, path: str, body=None, headers=None, opener=None):
    """Call /api/v1; return (status, parsed JSON)."""
    status, _, content = raw(method, "/api/v1" + path, body, headers, opener)
    return status, json.loads(content) if content else {}


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


def menu_models(menu):
    """All models of the menu, flattened."""
    out = []
    for cat in menu["categories"]:
        out += cat["models"]
        for sub in cat["children"]:
            out += sub["models"]
    return out


# ── Test groups ──────────────────────────────────────────────────────────────

def test_menu():
    print("\n── Library menu ──")
    status, menu = api("GET", "/library/menu")
    check("menu returns 200", status == 200)
    check("library is available", menu.get("available") is True, str(menu.get("reason")))
    models = menu_models(menu)
    check("menu lists models", len(models) > 0)
    check("total matches the listed models", menu["total"] == len(models))
    check("snapshot information present", "content_hash" in menu["snapshot"] and "commit" in menu["snapshot"])
    check("categories carry a count and a title",
          all("count" in c and c["title"] for c in menu["categories"]))

    # Only verified models: none of the blocked ones of the full index is in the menu
    _, _, content = raw("GET", "/library/library.json")
    index = json.loads(content)
    verified = {m["id"] for m in index["models"] if m["status"] == "verified"}
    blocked = {m["id"] for m in index["models"] if m["status"] != "verified"}
    listed = {m["id"] for m in models}
    check("menu = the verified models of library.json", listed == verified)
    check("no blocked/failing model in the menu", not (listed & blocked))
    return models, index


def test_open(models, index):
    print("\n── Open a library model ──")
    model = models[0]
    status, res = api("POST", "/library/open", {"id": model["id"]})
    check("open verified model succeeds", status == 200 and res.get("success"), str(res))
    check("response names the model", res.get("libraryId") == model["id"] and res.get("modelName") == model["name"])

    status, ees = api("GET", "/files/eescode")
    check("eescode loaded into the session", len(ees.get("content", "")) > 50)
    check("session remembers the library ID", ees.get("libraryId") == model["id"])
    check("model name set", ees.get("modelName") == model["name"])

    status, parsed = api("POST", "/parse", {"source": ees["content"]})
    check("the loaded model parses", status == 200 and parsed.get("success"), str(parsed)[:200])

    # A model is a folder: its companion files come with it (initials, solution baseline, tables)
    by_id = {m["id"]: m for m in index["models"]}
    for m in models:
        folder = "/library/" + by_id[m["id"]]["path"] + "/"
        st_i, _, initials = raw("GET", folder + m["name"] + ".initials")
        st_s, _, _ = raw("GET", folder + m["name"] + ".sol")
        if st_i == 200 and st_s == 200:
            _, res = api("POST", "/library/open", {"id": m["id"]})
            check("hasInitials / hasSol reported", res.get("hasInitials") and res.get("hasSol"), str(res))
            _, got = api("GET", "/files/initials")
            check("the .initials file is loaded", got["content"].strip() == initials.decode().strip())
            _, got = api("GET", "/files/sol")
            check("the .sol baseline is loaded", len(got["content"]) > 0)
            break
    else:
        check("some verified model has .initials and .sol", False)
    tabled = None
    for m in models:
        _, res = api("POST", "/library/open", {"id": m["id"]})
        if res.get("tables"):
            tabled = res
            break
    if tabled:
        _, listing = api("GET", "/tables")
        check("lookup tables of the model folder are loaded", len(listing.get("tables", [])) == tabled["tables"],
              str(listing)[:200])

    api("POST", "/library/open", {"id": model["id"]})  # the Back test below starts from this model
    # Opening a different model replaces everything, and Back restores the first
    other = models[1]
    status, res2 = api("POST", "/library/open", {"id": other["id"]})
    check("second model opens", status == 200 and res2.get("libraryId") == other["id"])
    status, back = api("POST", "/back")
    check("Back restores the previous model", status == 200 and back.get("modelName") == model["name"])
    _, ees = api("GET", "/files/eescode")
    check("Back restores its library ID", ees.get("libraryId") == model["id"])

    # By name, and a model that is not verified (blocked models open too)
    status, res = api("POST", "/library/open", {"id": other["name"]})
    check("open by name", status == 200 and res.get("libraryId") == other["id"])
    blocked = next((m for m in index["models"] if m["status"] == "blocked"), None)
    if blocked:
        status, res = api("POST", "/library/open", {"id": blocked["id"]})
        check("a blocked model opens", status == 200 and res.get("status") == "blocked", str(res))

    # A variant: the folder's other .eescode file
    with_variant = None
    for m in models:
        status, res = api("POST", "/library/open", {"id": m["id"]})
        if res.get("variants"):
            with_variant = (m, res["variants"][0])
            break
    if with_variant:
        status, res = api("POST", "/library/open", {"id": with_variant[0]["id"], "file": with_variant[1]})
        check("a variant opens", status == 200 and res.get("file") == with_variant[1], str(res))

    # Opening a plain file from disk clears the library ID
    api("POST", "/new")
    _, ees = api("GET", "/files/eescode")
    check("New model clears the library ID", ees.get("libraryId") == "")


def test_open_errors(models):
    print("\n── Open: errors ──")
    status, res = api("POST", "/library/open", {"id": "CSL-9999"})
    check("unknown model → 404", status == 404 and "Unknown library model" in res.get("error", ""), str(res))
    status, res = api("POST", "/library/open", {"id": models[0]["id"], "file": "../../../etc/passwd"})
    check("a path is not accepted as file → 404", status == 404, str(res))
    status, res = api("POST", "/library/open", {"id": models[0]["id"], "file": "README.md"})
    check("only .eescode files open → 404", status == 404, str(res))
    status, _, _ = raw("POST", "/api/v1/library/open", None)
    check("missing body → 400", status == 400)
    status, res = api("POST", "/library/open", {})
    check("missing id → 404", status == 404)


def test_tab_sessions(models):
    print("\n── Per-tab sessions ──")
    sid_a = "a" * 32
    sid_b = "b" * 32
    a, b = models[0], models[1]

    # The cookie session works on its own model
    api("POST", "/library/open", {"id": b["id"]})

    status, res = api("POST", "/library/open", {"id": a["id"]}, {"X-CoolSolve-Session": sid_a}, _plain_opener)
    check("a tab session opens a model", status == 200 and res.get("libraryId") == a["id"])

    _, ees_a = api("GET", "/files/eescode", None, {"X-CoolSolve-Session": sid_a}, _plain_opener)
    _, ees_cookie = api("GET", "/files/eescode")
    check("the tab session holds its own model", ees_a["libraryId"] == a["id"])
    check("the cookie session is untouched", ees_cookie["libraryId"] == b["id"])

    # ?sid= (EventSource and downloads cannot set headers) reaches the same session
    _, ees_q = api("GET", "/files/eescode?sid=" + sid_a, None, None, _plain_opener)
    check("?sid= selects the tab session", ees_q["libraryId"] == a["id"])

    # A second tab session is independent
    _, ees_b = api("GET", "/files/eescode", None, {"X-CoolSolve-Session": sid_b}, _plain_opener)
    check("another tab session starts empty", ees_b["libraryId"] == "" and ees_b["content"] == "")

    # The id names a temp folder: anything that is not 32 lowercase hex digits is ignored
    for bad in ["../../etc", "A" * 32, "a" * 31, "g" * 32, "a" * 33]:
        status, res = api("GET", "/files/eescode", None, {"X-CoolSolve-Session": bad}, _plain_opener)
        check(f"invalid session id {bad[:12]!r} is not used", status == 200 and res["libraryId"] == "")
    status, _ = api("GET", "/files/eescode?sid=../../etc", None, None, _plain_opener)
    check("invalid ?sid= does not break the request", status == 200)


def test_pages(index):
    print("\n── /library/ pages and files ──")
    status, ctype, content = raw("GET", "/library/")
    check("/library/ serves the dashboard", status == 200 and "text/html" in ctype and b"CoolSolve Library" in content)
    status, _, content = raw("GET", "/library")
    check("/library redirects to the dashboard", status == 200 and b"CoolSolve Library" in content)

    meta = index["models"][0]
    status, ctype, content = raw("GET", "/library/" + meta["path"] + "/README.md")
    check("a model README is served", status == 200 and len(content) > 100)
    status, ctype, content = raw("GET", "/library/" + meta["path"] + "/" + meta["main_file"])
    check("a model's .eescode is served", status == 200 and "text/plain" in ctype)
    status, _, content = raw("GET", "/docs/raw/library/" + meta["path"] + "/README.md")
    check("the docs viewer can fetch a model README", status == 200 and len(content) > 100)

    fig = next((m for m in index["models"] if m.get("figures")), None)
    if fig:
        status, ctype, _ = raw("GET", "/library/" + fig["path"] + "/figures/" + "x.png")
        check("a missing figure → 404", status == 404)

    check("not embedded: sources/ → 404", raw("GET", "/library/sources/thermo_models/inventory.csv")[0] == 404)
    check("not embedded: tools/ → 404", raw("GET", "/library/tools/build_index.py")[0] == 404)
    check("path traversal rejected", raw("GET", "/library/models/../../etc/passwd")[0] in (400, 404))
    check("the documentation links to the library",
          b"/library/" in raw("GET", "/docs/")[2])


# ── Runner ───────────────────────────────────────────────────────────────────

def main():
    global BASE_URL
    server_proc = None

    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if args:
        BASE_URL = args[0].rstrip("/")

    if "--auto" in sys.argv:
        binary = "./coolsolve"
        if not os.path.exists(binary):
            binary = "../build/coolsolve"
        if not os.path.exists(binary):
            print(f"ERROR: CoolSolve binary not found at {binary}")
            sys.exit(1)
        print(f"Starting CoolSolve from {binary}...")
        server_proc = subprocess.Popen([binary, "--gui", "18600", "--no-browser"],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        BASE_URL = "http://localhost:18600"
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
        status, health = api("GET", "/health")
        print(f"Connected to {BASE_URL} (status={health.get('status')})")
        models, index = test_menu()
        if not models:
            print("No library in this build: the remaining tests are skipped")
        else:
            test_open(models, index)
            test_open_errors(models)
            test_tab_sessions(models)
            test_pages(index)
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        check("no fatal error", False, str(e))
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
