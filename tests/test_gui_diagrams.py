#!/usr/bin/env python3
"""
test_gui_diagrams.py — End-to-end test of the Diagrams tab in a real browser.

Opens the GUI in headless Chromium (Chrome DevTools Protocol, as test_gui_shortcuts.py does) and checks
that the tab draws coolplot-core's SVG served by POST /api/v1/diagram/scene: the default chart after
load, the states of a solved model as cycle markers (GET /api/v1/diagram/states), a cycle of the state
editor drawn with paths along the processes, the hover (a state marker's variables; elsewhere the state
under the cursor from POST /api/v1/diagram/state), a rubber-band zoom and its reset, the downloads of the
chart (SVG and PDF from POST /api/v1/diagram/export, the chart spec for CoolPlot) and the export of the
diagram for the LaTeX report (thermo_diagrams.png, the fallback of the server's PDF figure). See
docs/coolplot_integration.md, phases 1, 2 and 4.

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_gui_diagrams.py [base_url]

    # Or run with auto-start (from the build/ directory):
    python3 test_gui_diagrams.py --auto

    # Also save screenshots of the tab (PNG) to a folder:
    python3 test_gui_diagrams.py --auto --screenshots /tmp/diagrams

Needs a GUI build (gui/dist embedded), Chromium or Chrome on the PATH (or in $CHROME) and the
Python package `websockets`. Without the browser or the package the test is skipped (exit 0).
Exit code 0 if all tests pass, 1 otherwise.
"""

import asyncio
import base64
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
DEVTOOLS_PORT = 18621
SERVER_PORT = 18620
SCREENSHOTS = None
WIDTH, HEIGHT = 1400, 900

# The chart drawn by DiagramView (not the icons of its buttons)
SVG = "document.querySelector('.thermo-diagram .diagram-svg-host > svg')"

# A vapour-compression cycle: four states, each the variables its property calls share
CYCLE_MODEL = """P_1 = 3e5
T_1 = 10
h_1 = enthalpy(R134a, P=P_1, T=T_1)
s_1 = entropy(R134a, P=P_1, h=h_1)
P_2 = 12e5
s_2 = s_1
h_2 = enthalpy(R134a, P=P_2, s=s_2)
T_2 = temperature(R134a, P=P_2, h=h_2)
P_3 = P_2
x_3 = 0
h_3 = enthalpy(R134a, P=P_3, x=x_3)
T_3 = temperature(R134a, P=P_3, x=x_3)
P_4 = P_1
h_4 = h_3
T_4 = temperature(R134a, P=P_4, h=h_4)
"""

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
    """A minimal CDP client on one page target: commands, the scene requests, and the paused report request."""

    def __init__(self, ws):
        self.ws = ws
        self.next_id = 0
        self.pending: dict[int, asyncio.Future] = {}
        self.scene_requests: list[dict] = []   # bodies of POST /api/v1/diagram/scene
        self.point_requests: list[dict] = []   # bodies of POST /api/v1/diagram/state
        self.export_requests: list[dict] = []  # bodies of POST /api/v1/diagram/export
        self.report_plots: list | None = None   # the plots of POST /api/v1/latex/compile
        self.reader = asyncio.ensure_future(self._read())

    async def _read(self):
        async for message in self.ws:
            msg = json.loads(message)
            if "id" in msg and msg["id"] in self.pending:
                self.pending.pop(msg["id"]).set_result(msg)
            elif msg.get("method") == "Network.requestWillBeSent":
                req = msg["params"]["request"]
                if req["method"] == "POST" and req["url"].split("?")[0].endswith("/api/v1/diagram/scene"):
                    body = json.loads(req.get("postData") or "{}")
                    self.scene_requests.append(body)
                    if VERBOSE:
                        print(f"        POST /diagram/scene {body.get('fluid')} {body.get('figure_size_pt')}")
                elif req["method"] == "POST" and req["url"].split("?")[0].endswith("/api/v1/diagram/state"):
                    self.point_requests.append(json.loads(req.get("postData") or "{}"))
                elif req["method"] == "POST" and req["url"].split("?")[0].endswith("/api/v1/diagram/export"):
                    self.export_requests.append(json.loads(req.get("postData") or "{}"))
            elif msg.get("method") == "Fetch.requestPaused":
                # The report: read the plots, then fail the request so that nothing is compiled or downloaded
                params = msg["params"]
                url = params["request"]["url"]
                if url.split("?")[0].endswith("/api/v1/latex/compile"):
                    try:
                        self.report_plots = json.loads(params["request"].get("postData") or "{}").get("plots", [])
                    except ValueError:
                        self.report_plots = []
                asyncio.ensure_future(self.send("Fetch.failRequest", requestId=params["requestId"], errorReason="Aborted"))

    async def send(self, method: str, **params):
        self.next_id += 1
        fut = asyncio.get_running_loop().create_future()
        self.pending[self.next_id] = fut
        await self.ws.send(json.dumps({"id": self.next_id, "method": method, "params": params}))
        msg = await asyncio.wait_for(fut, 60)
        if "error" in msg:
            raise RuntimeError(f"{method}: {msg['error']}")
        return msg.get("result", {})

    async def eval(self, expression: str):
        res = await self.send("Runtime.evaluate", expression=expression, returnByValue=True, awaitPromise=True)
        return res.get("result", {}).get("value")

    async def wait_for(self, expression: str, timeout: float = 30.0) -> bool:
        end = time.time() + timeout
        while time.time() < end:
            if await self.eval(expression):
                return True
            await asyncio.sleep(0.2)
        return False

    async def idle(self, settle: float = 0.5):
        """Wait until no scene request is on its way (the spinner of DiagramView is gone)."""
        await asyncio.sleep(settle)
        await self.wait_for("!document.querySelector('.diagram-spinner')", 60)
        await asyncio.sleep(0.2)

    async def set_model(self, code: str):
        await self.eval(f"window.monaco.editor.getModels()[0].setValue({json.dumps(code)})")
        await asyncio.sleep(0.3)

    async def mouse(self, kind: str, x: float, y: float, buttons: int):
        await self.send("Input.dispatchMouseEvent", type=kind, x=x, y=y, button="left", buttons=buttons, clickCount=1)

    async def screenshot(self, name: str, selector: str = ".diagram-tab"):
        if not SCREENSHOTS:
            return
        box = await self.eval(f"(() => {{ const r = document.querySelector({json.dumps(selector)}).getBoundingClientRect();"
                              " return [r.left, r.top, r.width, r.height]; })()")
        shot = await self.send("Page.captureScreenshot", format="png",
                               clip={"x": box[0], "y": box[1], "width": box[2], "height": box[3], "scale": 2})
        path = os.path.join(SCREENSHOTS, name)
        with open(path, "wb") as f:
            f.write(base64.b64decode(shot["data"]))
        print(f"        screenshot: {path}")


async def run_tests(page: Page):
    await page.send("Network.enable")
    await page.send("Page.enable")
    await page.send("Emulation.setDeviceMetricsOverride", width=WIDTH, height=HEIGHT, deviceScaleFactor=1, mobile=False)
    await page.send("Page.navigate", url=BASE_URL + "/")
    check("the GUI loads with the Monaco editor",
          await page.wait_for("!!(window.monaco && document.querySelector('.monaco-editor textarea'))", 30))

    print("\n=== The default chart ===")
    await page.eval("[...document.querySelectorAll('.tab-btn')].find(b => b.textContent === 'Diagrams').click()")
    check("the Diagrams tab shows coolplot-core's SVG with the saturation curve (iso/Q/0)",
          await page.wait_for(f"!!{SVG}?.querySelector('[data-item=\"iso/Q/0\"]')", 60))
    check("the SVG is inline (an <svg> element, not an image)",
          await page.eval(f"{SVG} instanceof SVGSVGElement"))
    size = await page.eval(f"(() => {{ const c = document.querySelector('.diagram-canvas').getBoundingClientRect();"
                           f" return [c.width, c.height, +{SVG}.getAttribute('width'), +{SVG}.getAttribute('height')]; }})()")
    req = page.scene_requests[-1] if page.scene_requests else {}
    fig = req.get("figure_size_pt") or [0, 0]
    check("the figure is the size of its container (pt = px * 0.75)",
          size and abs(fig[0] - size[0] * 0.75) < 2 and abs(fig[1] - size[1] * 0.75) < 2
          and abs(size[2] - size[0]) < 3 and abs(size[3] - size[1]) < 3, f"canvas/svg {size}, figure {fig}")
    spec = req.get("spec", {})
    check("the spec is CoolPlot's: ph, engineering units, DEF, CoolPlot's default badge (no branding written)",
          spec.get("diagram") == "ph" and spec.get("units") == "engineering" and spec.get("limits") == "DEF"
          and "branding" not in spec.get("settings", {}), json.dumps(spec)[:300])
    check("the badge names the equation of state and the versions",
          await page.wait_for(f"(() => {{ const t = {SVG}?.textContent || ''; return t.includes('Wagner and Pru') "
                              "&& t.includes('CoolPlot 1.5.0'); })()", 10))
    check("no error shown", not await page.eval("!!document.querySelector('.diagram-status .diagram-error')"),
          await page.eval("document.querySelector('.diagram-status')?.textContent"))
    await page.idle()
    await page.screenshot("diagram_default.png")

    print("\n=== States of a solved model ===")
    await page.set_model(CYCLE_MODEL)
    await page.eval("[...document.querySelectorAll('button')].find(b => b.textContent.trim() === 'Solve').click()")
    check("the model solves", await page.wait_for("document.body.innerText.includes('Solve SUCCESS')", 60))
    check("the model's fluid is drawn (R134a)",
          await page.wait_for("document.querySelector('.diagram-controls select')?.value === 'R134a'", 30))
    check("the cycle's state markers appear (cp-role-state_points)",
          await page.wait_for(f"!!{SVG}?.querySelector('.cp-role-state_points')", 60))
    await page.idle()
    numbers = await page.eval(f"[...{SVG}.querySelectorAll('.cp-role-state_number')].map(g => g.textContent.trim())")
    check("the four states of the model are numbered on the chart", len(numbers or []) >= 4, f"{numbers}")
    check("markers alone: no line through the inferred states",
          not await page.eval(f"!!{SVG}.querySelector('.cp-role-cycle')"))
    cycles = page.scene_requests[-1].get("cycles", []) if page.scene_requests else []
    inferred = cycles[0] if cycles else {}
    check("the request names the server's states by their variables, loops []",
          len(cycles) == 1 and inferred.get("loops") == []
          and any(st["vars"].get("p") == "P_1" for st in inferred.get("states", [])), json.dumps(cycles)[:300])
    await page.screenshot("diagram_cycle_states.png")

    print("\n=== A cycle of the state editor ===")
    await page.eval("[...document.querySelectorAll('.diagram-cycles button')]"
                    ".find(b => b.textContent.includes('from inferred states')).click()")
    end = time.time() + 20
    while time.time() < end and not any(c.get("loops") for c in (page.scene_requests[-1].get("cycles", [])
                                                                 if page.scene_requests else [])):
        await asyncio.sleep(0.1)
    cycles = page.scene_requests[-1].get("cycles", []) if page.scene_requests else []
    user = cycles[0] if cycles else {}
    check("'from inferred states' makes an editable cycle through the states, closed, paths auto; the overlay goes",
          len(cycles) == 1 and user.get("name") == "Cycle" and user.get("paths") == "auto"
          and len(user.get("loops", [[]])[0]) == len(user.get("states", [])) + 1
          and user["loops"][0][0] == user["loops"][0][-1], json.dumps(cycles)[:400])
    check("the editor shows the cycle with its states",
          await page.wait_for("document.querySelectorAll('.diagram-cycle .diagram-state').length >= 4", 10))
    check("the cycle is drawn with lines (cp-role-cycle)",
          await page.wait_for(f"!!{SVG}?.querySelector('.cp-role-cycle path')", 60))
    await page.idle()
    points = await page.eval(f"[...{SVG}.querySelectorAll('.cp-role-cycle path')]"
                             ".map(p => (p.getAttribute('d').match(/[ML]/g) || []).length).reduce((a, b) => a + b, 0)")
    check("auto paths follow the processes (more points than straight segments)",
          (points or 0) > len(user.get("states", [])) + 4, f"{points} points")
    await page.screenshot("diagram_editor_cycle.png")

    print("\n=== Hover ===")
    marker = await page.eval(f"(() => {{ const g = {SVG}.querySelector('.cp-role-state_points');"
                             " const m = g && g.firstElementChild; if (!m) return null;"
                             " const r = m.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()")
    if marker:
        await page.mouse("mouseMoved", marker[0], marker[1], 0)
    check("hovering a state marker shows the cycle, the state and its variables with their values",
          await page.wait_for("(() => { const t = document.querySelector('.diagram-tooltip[data-kind=\"state\"]');"
                              " return !!t && t.textContent.includes('Cycle: state') && /[PhTs]_\\d = /.test(t.textContent); })()", 10),
          await page.eval("document.querySelector('.diagram-tooltip')?.textContent"))
    box = await page.eval(f"(() => {{ const r = {SVG}.querySelector('rect.cp-background').getBoundingClientRect();"
                          " return [r.left, r.top, r.width, r.height]; })()")
    count = len(page.point_requests)
    await page.mouse("mouseMoved", box[0] + 0.12 * box[2], box[1] + 0.55 * box[3], 0)
    check("hovering elsewhere asks for the state under the cursor and shows T, p, h, s, rho",
          await page.wait_for("(() => { const t = document.querySelector('.diagram-tooltip[data-kind=\"point\"]');"
                              " return !!t && ['T = ', 'p = ', 'h = ', 's = ', 'ρ = '].every(k => t.textContent.includes(k)); })()", 15)
          and len(page.point_requests) > count,
          await page.eval("document.querySelector('.diagram-tooltip')?.textContent"))
    body = page.point_requests[-1] if page.point_requests else {}
    check("the state request names the fluid, the diagram and the point in SI",
          body.get("fluid") == "R134a" and body.get("diagram") == "ph" and body.get("x_SI", 0) > 1000, json.dumps(body))
    await page.screenshot("diagram_hover.png")
    await page.mouse("mouseMoved", 5, 5, 0)
    check("leaving the plot hides the tooltip", await page.wait_for("!document.querySelector('.diagram-tooltip')", 5))

    print("\n=== Zoom and reset ===")
    before = await page.eval(f"{SVG}.getAttribute('data-x-limits')")
    box = await page.eval(f"(() => {{ const r = {SVG}.querySelector('rect.cp-background').getBoundingClientRect();"
                          " return [r.left, r.top, r.width, r.height]; })()")
    x0, y0 = box[0] + 0.25 * box[2], box[1] + 0.25 * box[3]
    x1, y1 = box[0] + 0.75 * box[2], box[1] + 0.75 * box[3]
    count = len(page.scene_requests)
    await page.mouse("mouseMoved", x0, y0, 0)
    await page.mouse("mousePressed", x0, y0, 1)
    for k in range(1, 6):
        await page.mouse("mouseMoved", x0 + (x1 - x0) * k / 5, y0 + (y1 - y0) * k / 5, 1)
    await page.mouse("mouseReleased", x1, y1, 0)
    end = time.time() + 15   # the request leaves after the debounce
    while time.time() < end and not (len(page.scene_requests) > count and page.scene_requests[-1]["spec"].get("view_SI")):
        await asyncio.sleep(0.1)
    check("a dragged rectangle asks for a zoomed scene (view_SI)",
          len(page.scene_requests) > count and page.scene_requests[-1]["spec"].get("view_SI") is not None
          and await page.wait_for("!!document.querySelector('.diagram-reset')", 15))
    await page.idle()
    after = await page.eval(f"{SVG}.getAttribute('data-x-limits')")
    view = page.scene_requests[-1]["spec"].get("view_SI") if page.scene_requests else None
    lo, hi = (float(v) for v in after.split())
    # h in kJ/kg on the chart, J/kg in view_SI
    check("the chart shows the zoom window", after != before and view is not None
          and abs(lo - view[0][0] / 1000) < 1e-6 * max(1, abs(lo)) and abs(hi - view[0][1] / 1000) < 1e-6 * max(1, abs(hi)),
          f"{before} -> {after}, view {view}")
    check("the labels of the scene on screen are sent back", bool(page.scene_requests[-1].get("previous_labels")))
    await page.eval("document.querySelector('.diagram-reset').click()")
    check("reset: the whole chart again",
          await page.wait_for(f"{SVG}?.getAttribute('data-x-limits') === {json.dumps(before)}", 30))
    await page.idle()

    print("\n=== Downloads: SVG, PDF, chart spec ===")
    downloads = tempfile.mkdtemp(prefix="coolsolve_diagram_downloads_")
    try:
        await page.send("Page.setDownloadBehavior", behavior="allow", downloadPath=downloads)
        check("the status line offers SVG, PDF and Spec",
              await page.eval("[...document.querySelectorAll('.diagram-download')].map(b => b.dataset.format).join() === 'svg,pdf,spec'"))
        for fmt in ("svg", "pdf", "spec"):
            await page.eval(f"document.querySelector('.diagram-download[data-format=\"{fmt}\"]').click()")
            await page.wait_for("[...document.querySelectorAll('.diagram-download')].every(b => !b.disabled)", 30)
        end = time.time() + 30
        names = []
        while time.time() < end:
            names = sorted(n for n in os.listdir(downloads) if not n.endswith(".crdownload"))
            if len(names) >= 3:
                break
            await asyncio.sleep(0.2)
        check("three files downloaded, named after the model and the diagram",
              names == ["diagram_ph.json", "diagram_ph.pdf", "diagram_ph.svg"], f"{names}")
        body = page.export_requests[-1] if page.export_requests else {}
        check("the export asks for the chart on screen: its fluid, cycles, figure size and labels",
              [r.get("format") for r in page.export_requests[-2:]] == ["svg", "pdf"] and body.get("fluid") == "R134a"
              and body.get("cycles") and body.get("figure_size_pt") == page.scene_requests[-1].get("figure_size_pt")
              and body.get("previous_labels") is not None, json.dumps(body)[:300])
        if len(names) >= 3:
            with open(os.path.join(downloads, "diagram_ph.svg"), "rb") as f:
                svg = f.read()
            with open(os.path.join(downloads, "diagram_ph.pdf"), "rb") as f:
                pdf = f.read()
            with open(os.path.join(downloads, "diagram_ph.json"), encoding="utf-8") as f:
                spec = json.load(f)
            check("the SVG is coolplot-core's, with the cycle", svg.startswith(b"<svg") and b'data-item="cycle/Cycle/loop/0"' in svg)
            check("the PDF is coolplot-core's, with its fonts", pdf.startswith(b"%PDF-") and b"/FontFile2" in pdf)
            check("the spec is the one composed: fluid and cycles written, as CoolPlot reads it",
                  spec.get("diagram") == "ph" and spec.get("fluid", {}).get("components") == ["R134a"]
                  and spec.get("cycles", [{}])[0].get("name") == "Cycle" and "states" in spec["cycles"][0], json.dumps(spec)[:300])
        check("no download error shown", not await page.eval("!!document.querySelector('.diagram-status .diagram-error')"),
              await page.eval("document.querySelector('.diagram-status')?.textContent"))
    finally:
        shutil.rmtree(downloads, ignore_errors=True)

    print("\n=== LaTeX report export ===")
    await page.send("Fetch.enable", patterns=[{"urlPattern": "*/api/v1/latex/*"}])
    await page.eval("document.querySelector('button[title*=\"PDF report\"]').click()")
    end = time.time() + 30
    while page.report_plots is None and time.time() < end:
        await asyncio.sleep(0.2)
    plots = page.report_plots or []
    diagram = next((p for p in plots if p.get("name") == "thermo_diagrams.png"), None)
    check("the report gets the diagram as thermo_diagrams.png (a PNG)",
          diagram is not None and diagram.get("data", "").startswith("data:image/png;base64,")
          and len(diagram["data"]) > 10000, f"{[p.get('name') for p in plots]}")
    await asyncio.sleep(1.0)
    await page.send("Fetch.disable")

    if SCREENSHOTS:
        print("\n=== Screenshots: dark theme, T-s chart ===")
        await page.eval("document.querySelector('button[title=\"Toggle theme\"]').click()")
        await page.wait_for(f"{SVG}?.querySelector('rect.cp-figure')?.getAttribute('fill') !== '#fcfcfb'", 30)
        await page.idle()
        await page.screenshot("diagram_dark.png")
        await page.eval("document.querySelector('button[title=\"Toggle theme\"]').click()")
        await page.eval("(() => { const s = document.querySelectorAll('.diagram-controls select')[1];"
                        " Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set.call(s, 'Ts');"
                        " s.dispatchEvent(new Event('change', {bubbles: true})); })()")
        await page.wait_for(f"{SVG}?.getAttribute('data-x-limits') !== {json.dumps(before)}", 30)
        await page.idle(1.0)
        await page.screenshot("diagram_ts.png")


async def run_browser(browser: str):
    import websockets
    profile = tempfile.mkdtemp(prefix="coolsolve_diagrams_")
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
    # In a temporary folder: the server writes its log and the session files where it runs
    workdir = tempfile.mkdtemp(prefix="coolsolve_diagrams_srv_")
    proc = subprocess.Popen([os.path.abspath(binary), "--gui", str(port), "--no-browser"], cwd=workdir,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(20):
        time.sleep(0.5)
        try:
            urllib.request.urlopen(BASE_URL + "/api/v1/health", timeout=2)
            return proc, workdir
        except Exception:
            pass
    proc.terminate()
    print("ERROR: Server did not start in 10s")
    sys.exit(1)


def main():
    global BASE_URL, SCREENSHOTS
    browser = find_browser()
    try:
        import websockets  # noqa: F401
    except ImportError:
        browser = None
    if not browser:
        print("SKIPPED: needs Chromium/Chrome on the PATH (or $CHROME) and the Python package websockets")
        sys.exit(0)

    argv = sys.argv[1:]
    if "--screenshots" in argv:
        i = argv.index("--screenshots")
        SCREENSHOTS = os.path.abspath(argv[i + 1])
        os.makedirs(SCREENSHOTS, exist_ok=True)
        del argv[i:i + 2]
    args = [a for a in argv if not a.startswith("-")]
    if args:
        BASE_URL = args[0].rstrip("/")

    server_proc, workdir = None, None
    if "--auto" in argv:
        binary = "./coolsolve"
        if not os.path.exists(binary):
            binary = "../build/coolsolve"
        if not os.path.exists(binary):
            print(f"ERROR: CoolSolve binary not found at {binary}")
            sys.exit(1)
        BASE_URL = f"http://localhost:{SERVER_PORT}"
        server_proc, workdir = start_server(binary, SERVER_PORT)

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
            shutil.rmtree(workdir, ignore_errors=True)

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
