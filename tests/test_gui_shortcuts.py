#!/usr/bin/env python3
"""
test_gui_shortcuts.py — End-to-end test of the GUI keyboard shortcuts in a real browser.

Opens the GUI in headless Chromium, driven over the Chrome DevTools Protocol, presses keys
as a user would (trusted key events, so the browser and Monaco see them exactly like real
key presses) and checks the requests that reach the CoolSolve server. Covers F2 (Solve,
or Try Harder after a failed solve, as in EES); see the comment on the F2 handler in
gui/src/components/Toolbar.tsx.

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_gui_shortcuts.py [base_url]

    # Or run with auto-start (from the build/ directory):
    python3 test_gui_shortcuts.py --auto

Needs a GUI build (gui/dist embedded), Chromium or Chrome on the PATH (or in $CHROME) and the
Python package `websockets`. Without the browser or the package the test is skipped (exit 0).
Exit code 0 if all tests pass, 1 otherwise.
"""

import asyncio
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import urllib.request

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv
DEVTOOLS_PORT = 18611
# The Solve button shows "Solving..." / "Searching..." while a solve runs
SOLVING = "[...document.querySelectorAll('button')].some(b => /Solving|Searching/.test(b.textContent))"

passed = 0
failed = 0
errors: list[str] = []


def check(name: str, cond: bool, detail: str = ""):
    global passed, failed
    if cond:
        passed += 1
        print(f"  PASS  {name}")
    else:
        failed += 1
        errors.append(name)
        print(f"  FAIL  {name}" + (f"  ({detail})" if detail else ""))


def find_browser():
    for name in (os.environ.get("CHROME"), "chromium", "chromium-browser", "google-chrome", "chrome"):
        if name and shutil.which(name):
            return shutil.which(name)
    return None


class Page:
    """A minimal CDP client on one page target: commands, and the solve requests it sends."""

    def __init__(self, ws):
        self.ws = ws
        self.next_id = 0
        self.pending: dict[int, asyncio.Future] = {}
        self.solve_requests: list[dict] = []   # bodies of POST /api/v1/solve
        self.reader = asyncio.ensure_future(self._read())

    async def _read(self):
        async for message in self.ws:
            msg = json.loads(message)
            if "id" in msg and msg["id"] in self.pending:
                self.pending.pop(msg["id"]).set_result(msg)
            elif msg.get("method") == "Network.requestWillBeSent":
                req = msg["params"]["request"]
                if req["method"] == "POST" and req["url"].split("?")[0].endswith("/api/v1/solve"):
                    self.solve_requests.append(json.loads(req.get("postData") or "{}"))
                    if VERBOSE:
                        print(f"        POST /solve {req.get('postData')}")

    async def send(self, method: str, **params):
        self.next_id += 1
        fut = asyncio.get_running_loop().create_future()
        self.pending[self.next_id] = fut
        await self.ws.send(json.dumps({"id": self.next_id, "method": method, "params": params}))
        msg = await asyncio.wait_for(fut, 30)
        if "error" in msg:
            raise RuntimeError(f"{method}: {msg['error']}")
        return msg.get("result", {})

    async def eval(self, expression: str):
        res = await self.send("Runtime.evaluate", expression=expression, returnByValue=True, awaitPromise=True)
        return res.get("result", {}).get("value")

    async def wait_for(self, expression: str, timeout: float = 15.0) -> bool:
        end = time.time() + timeout
        while time.time() < end:
            if await self.eval(expression):
                return True
            await asyncio.sleep(0.2)
        return False

    async def press(self, key: str, code: str, vk: int, modifiers: int = 0):
        """A trusted key press (modifiers: 1 Alt, 2 Ctrl, 4 Meta, 8 Shift)."""
        for kind in ("rawKeyDown", "keyUp"):
            await self.send("Input.dispatchKeyEvent", type=kind, key=key, code=code,
                            windowsVirtualKeyCode=vk, nativeVirtualKeyCode=vk, modifiers=modifiers)

    async def press_f2(self, modifiers: int = 0):
        await self.press("F2", "F2", 113, modifiers)

    async def set_model(self, code: str):
        """Replace the editor content (as typing would: the GUI store follows Monaco) and focus it."""
        await self.eval(f"window.monaco.editor.getModels()[0].setValue({json.dumps(code)})")
        await self.eval("document.querySelector('.monaco-editor textarea').focus()")
        await asyncio.sleep(0.3)

    async def solve_count_after(self, action, settle: float = 1.0) -> int:
        """Number of POST /solve requests sent within `settle` seconds of `action`."""
        before = len(self.solve_requests)
        await action
        await asyncio.sleep(settle)
        return len(self.solve_requests) - before

    def console_text(self):
        return "document.body.innerText"


async def run_tests(page: Page):
    await page.send("Network.enable")
    await page.send("Page.enable")
    await page.send("Page.navigate", url=BASE_URL + "/")
    check("the GUI loads with the Monaco editor",
          await page.wait_for("!!(window.monaco && document.querySelector('.monaco-editor textarea'))", 30))

    print("\n=== F2: negative cases ===")
    await page.set_model("")
    n = await page.solve_count_after(page.press_f2())
    check("F2 on an empty model sends no solve", n == 0, f"{n} request(s)")

    await page.set_model("x = 3\ny = x^2 + 1")
    for name, mod in (("Shift+F2", 8), ("Ctrl+F2", 2), ("Alt+F2", 1)):
        n = await page.solve_count_after(page.press_f2(mod))
        check(f"{name} sends no solve", n == 0, f"{n} request(s)")

    print("\n=== F2 in the editor: Solve ===")
    n = await page.solve_count_after(page.press_f2(), settle=0.5)
    check("F2 with the focus in the editor sends one POST /solve", n == 1, f"{n} request(s)")
    body = page.solve_requests[-1] if page.solve_requests else {}
    check("the solve is a normal one (no debug, no deep search)",
          not body.get("debug") and not body.get("deepSearch"), json.dumps(body))
    check("the GUI reports the solve success",
          await page.wait_for(f"{page.console_text()}.includes('Solve SUCCESS')", 30))
    # Read back in the page, so that the request carries the page's session (cookie or tab id)
    result = await page.eval("fetch('/api/v1/solve/result').then(r => r.json())") or {}
    y = result.get("variables", {}).get("y", {})
    y = y.get("value") if isinstance(y, dict) else y
    check("the server solved the model of the editor (y = 10)", y is not None and abs(float(y) - 10) < 1e-9,
          f"y={y}")
    check("Monaco did not take F2 (no rename box)", not await page.eval("!!document.querySelector('.rename-box')"))
    check("the model text is unchanged",
          await page.eval("window.monaco.editor.getModels()[0].getValue()") == "x = 3\ny = x^2 + 1")

    print("\n=== F2 outside the editor ===")
    await page.eval("document.activeElement && document.activeElement.blur()")
    n = await page.solve_count_after(page.press_f2(), settle=0.5)
    check("F2 with the focus outside the editor sends one POST /solve", n == 1, f"{n} request(s)")
    await page.wait_for(f"{page.console_text()}.includes('Solve SUCCESS')", 30)

    print("\n=== F2 after a failed solve: Try Harder ===")
    await page.set_model("x^2 = -1")
    await page.solve_count_after(page.press_f2(), settle=0.5)
    check("the impossible model fails and arms Try Harder",
          await page.wait_for("[...document.querySelectorAll('button')].some(b => b.textContent.includes('Try Harder'))",
                              60))
    n = await page.solve_count_after(page.press_f2(), settle=0.5)
    body = page.solve_requests[-1]
    check("F2 then sends one Try Harder solve (deepSearch: true)",
          n == 1 and body.get("deepSearch") is True, f"{n} request(s), {json.dumps(body)}")
    await page.press("Escape", "Escape", 27)   # stop the deep search
    await page.wait_for(f"!{SOLVING}", 60)

    print("\n=== F2 during a solve ===")
    # A long time march (10^6 steps), stopped with Escape once checked
    await page.set_model("y = 1 + INTEGRAL(dydt, t, 0, 1000, 0.001)\ndydt = -y")
    started = len(page.solve_requests)
    await page.press_f2()
    # The button shows "Solving..." before its POST /solve leaves (the model is PUT first)
    end = time.time() + 15
    while len(page.solve_requests) == started and time.time() < end:
        await asyncio.sleep(0.1)
    running = len(page.solve_requests) == started + 1 and await page.wait_for(SOLVING, 15)
    n = await page.solve_count_after(page.press_f2(), settle=0.5)
    check("F2 while a solve runs sends nothing", running and n == 0 and await page.eval(SOLVING),
          f"running={running}, {n} request(s)")
    await page.press("Escape", "Escape", 27)
    check("Escape stops it", await page.wait_for(f"!{SOLVING}", 30))


async def run_browser(browser: str):
    import websockets
    profile = tempfile.mkdtemp(prefix="coolsolve_shortcuts_")
    proc = subprocess.Popen([browser, "--headless=new", f"--remote-debugging-port={DEVTOOLS_PORT}",
                             f"--user-data-dir={profile}", "--no-first-run", "--no-default-browser-check",
                             "--disable-gpu", "about:blank"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        targets = None
        for _ in range(60):
            time.sleep(0.5)
            try:
                with urllib.request.urlopen(f"http://127.0.0.1:{DEVTOOLS_PORT}/json/list", timeout=2) as r:
                    targets = [t for t in json.loads(r.read()) if t.get("type") == "page"]
                if targets:
                    break
            except Exception:
                pass
        if not targets:
            check("the browser starts with its DevTools endpoint", False, browser)
            return
        async with websockets.connect(targets[0]["webSocketDebuggerUrl"], max_size=None) as ws:
            page = Page(ws)
            try:
                await run_tests(page)
            finally:
                page.reader.cancel()
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        shutil.rmtree(profile, ignore_errors=True)


# ── Runner ───────────────────────────────────────────────────────────────────

def start_server(binary: str, port: int):
    print(f"Starting CoolSolve from {binary} on port {port}...")
    proc = subprocess.Popen([binary, "--gui", str(port), "--no-browser"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(20):
        time.sleep(0.5)
        try:
            urllib.request.urlopen(BASE_URL + "/api/v1/health", timeout=2)
            return proc
        except Exception:
            pass
    proc.terminate()
    print("ERROR: Server did not start in 10s")
    sys.exit(1)


def main():
    global BASE_URL
    browser = find_browser()
    try:
        import websockets  # noqa: F401
    except ImportError:
        browser = None
    if not browser:
        print("SKIPPED: needs Chromium/Chrome on the PATH (or $CHROME) and the Python package websockets")
        sys.exit(0)

    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if args:
        BASE_URL = args[0].rstrip("/")

    server_proc = None
    if "--auto" in sys.argv:
        binary = "./coolsolve"
        if not os.path.exists(binary):
            binary = "../build/coolsolve"
        if not os.path.exists(binary):
            print(f"ERROR: CoolSolve binary not found at {binary}")
            sys.exit(1)
        BASE_URL = "http://localhost:18610"
        server_proc = start_server(binary, 18610)

    try:
        print(f"Browser: {browser}\nGUI: {BASE_URL}")
        asyncio.run(run_browser(browser))
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
