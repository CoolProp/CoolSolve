#!/usr/bin/env python3
"""
test_gui_diagrams.py — End-to-end test of the Diagrams tab in a real browser.

Opens the GUI in headless Chromium (Chrome DevTools Protocol, as test_gui_shortcuts.py does) and checks
that the tab draws coolplot-core's SVG served by POST /api/v1/diagram/scene: the default chart after
load (in the model's units, short diagram names, no unit selector), the states of a solved model as
cycle markers (GET /api/v1/diagram/states), the hover (a state marker's variables; elsewhere the state
under the cursor from POST /api/v1/diagram/state), a rubber-band zoom and its reset, two arrays overlays
drawn with paths along the processes (their default arrays those of the axis quantities, mapped to the
axes of another diagram type), the downloads of the chart (SVG and PDF from POST /api/v1/diagram/export,
the chart spec for CoolPlot), the export of the diagram for the LaTeX report (thermo_diagrams.png, the
fallback of the server's PDF figure), and the other kinds of fluids: an ideal gas with its solved cycle,
an incompressible solution at a concentration, a blend of the model, humid air (a solved
air-conditioning process on the psychrometric chart at the model's total pressure, the hover, the
Mollier chart). See docs/coolplot_integration.md, phases 1 to 4, and docs/gui.md §6.8.

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

# Two vapour-compression cycles of R134a in arrays. P, h, T and s are filled from the scalars of the cycle above, as
# the library models do (the tab knows their quantities from their names); Pb, hb, Tb and sb are arguments and results
# of property calls (the variable inference gives their quantities)
ARRAYS_MODEL = CYCLE_MODEL + """s_3 = entropy(R134a, P=P_3, x=x_3)
s_4 = entropy(R134a, P=P_4, h=h_4)
""" + "".join(f"P[{i}] = P_{i}\nh[{i}] = h_{i}\nT[{i}] = T_{i}\ns[{i}] = s_{i}\n" for i in range(1, 5)) + "".join(f"""P{b}[1] = {p_ev}
T{b}[1] = 10
h{b}[1] = enthalpy(R134a, P=P{b}[1], T=T{b}[1])
s{b}[1] = entropy(R134a, P=P{b}[1], h=h{b}[1])
P{b}[2] = {p_cd}
s{b}[2] = s{b}[1]
h{b}[2] = enthalpy(R134a, P=P{b}[2], s=s{b}[2])
T{b}[2] = temperature(R134a, P=P{b}[2], h=h{b}[2])
P{b}[3] = P{b}[2]
h{b}[3] = enthalpy(R134a, P=P{b}[3], x=0)
T{b}[3] = temperature(R134a, P=P{b}[3], x=0)
s{b}[3] = entropy(R134a, P=P{b}[3], x=0)
P{b}[4] = P{b}[1]
h{b}[4] = h{b}[3]
T{b}[4] = temperature(R134a, P=P{b}[4], h=h{b}[4])
s{b}[4] = entropy(R134a, P=P{b}[4], h=h{b}[4])
""" for b, p_ev, p_cd in (("b", "2e5", "9e5"),))

# Phase 3: an air-standard Brayton cycle of the ideal gas N2 in arrays, and a blend of the model
BRAYTON_MODEL = """T[1] = 25
P[1] = 1e5
P[2] = 8*P[1]
s[1] = entropy(N2, T=T[1], P=P[1])
s[2] = s[1]
h[2] = enthalpy(N2, P=P[2], s=s[2])
T[2] = temperature(N2, P=P[2], s=s[2])
P[3] = P[2]
T[3] = 1000
s[3] = entropy(N2, T=T[3], P=P[3])
P[4] = P[1]
s[4] = s[3]
h[4] = enthalpy(N2, P=P[4], s=s[4])
T[4] = temperature(N2, P=P[4], s=s[4])
"""
BLEND_MODEL = "T_b = temperature('R32[0.7]&R125[0.3]', P=10e5, x=0)\n"

# Humid air at 1 bar: outdoor air cooled and dehumidified, then reheated at constant humidity ratio
AIR_CONDITIONING_MODEL = """P = 1e5
T[1] = 32
RH[1] = 0.5
w[1] = humrat(AirH2O, T=T[1], P=P, R=RH[1])
h[1] = enthalpy(AirH2O, T=T[1], P=P, w=w[1])
T[2] = 13
RH[2] = 0.95
w[2] = humrat(AirH2O, T=T[2], P=P, R=RH[2])
h[2] = enthalpy(AirH2O, T=T[2], P=P, R=RH[2])
w[3] = w[2]
T[3] = 20
h[3] = enthalpy(AirH2O, T=T[3], P=P, w=w[3])
RH[3] = relhum(AirH2O, T=T[3], P=P, w=w[3])
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
    units = spec.get("units") if isinstance(spec.get("units"), dict) else {}
    check("the spec is CoolPlot's: ph, the model's units (°C, Pa, J/kg), DEF, CoolPlot's default badge (no branding written)",
          spec.get("diagram") == "ph" and units.get("name") == "coolsolve" and units.get("units", {}).get("p") == ["Pa", 1.0, 0.0]
          and units["units"].get("h") == ["J/kg", 1.0, 0.0] and spec.get("limits") == "DEF"
          and "branding" not in spec.get("settings", {}), json.dumps(spec)[:300])
    check("no unit selector", not await page.eval(
        "[...document.querySelectorAll('.diagram-controls .diagram-field > span')].some(s => s.textContent === 'Units')"))
    titles = await page.eval("[...document.querySelectorAll('.diagram-controls select')[1].options].map(o => o.textContent)")
    check("the diagram types have short names", titles == ["p-h", "T-s", "h-s", "p-s", "p-ρ", "T-ρ", "p-T", "T-h"], str(titles))
    overrides = spec.get("style_overrides", {})
    check("light lines by default: every family 0.5 pt at 55 % (CoolPlot's style_overrides)",
          {"T", "s", "rho", "Q"} <= set(overrides)
          and all(o == {"width_pt": 0.5, "alpha": 0.55} for o in overrides.values()), json.dumps(overrides)[:300])
    check("the isolines are drawn light (stroke-opacity 0.55)",
          await page.wait_for(f"(() => {{ const p = {SVG}?.querySelector('g.cp-role-isoline path'); "
                              "return !!p && p.getAttribute('stroke-opacity') === '0.55'; })()", 10))
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

    print("\n=== Hover ===")
    marker = await page.eval(f"(() => {{ const g = {SVG}.querySelector('.cp-role-state_points');"
                             " const m = g && g.firstElementChild; if (!m) return null;"
                             " const r = m.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()")
    if marker:
        await page.mouse("mouseMoved", marker[0], marker[1], 0)
    check("hovering a state marker shows the cycle, the state and its variables with their values",
          await page.wait_for("(() => { const t = document.querySelector('.diagram-tooltip[data-kind=\"state\"]');"
                              " return !!t && t.textContent.includes('Inferred states: state') && /[PhTs]_\\d = /.test(t.textContent); })()", 10),
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
    # h in J/kg on the chart (the model's units) and in view_SI
    check("the chart shows the zoom window", after != before and view is not None
          and abs(lo - view[0][0]) < 1e-6 * max(1, abs(lo)) and abs(hi - view[0][1]) < 1e-6 * max(1, abs(hi)),
          f"{before} -> {after}, view {view}")
    check("the labels of the scene on screen are sent back", bool(page.scene_requests[-1].get("previous_labels")))
    await page.eval("document.querySelector('.diagram-reset').click()")
    check("reset: the whole chart again",
          await page.wait_for(f"{SVG}?.getAttribute('data-x-limits') === {json.dumps(before)}", 30))
    await page.idle()

    await arrays_overlays(page)

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
            check("the SVG is coolplot-core's, with the arrays", svg.startswith(b"<svg") and b'data-item="cycle/P, h/loop/0"' in svg)
            check("the PDF is coolplot-core's, with its fonts", pdf.startswith(b"%PDF-") and b"/FontFile2" in pdf)
            arrays = next((c for c in spec.get("cycles", []) if c.get("name") == "P, h"), {})
            check("the spec is the one composed: fluid and cycles written (their states), as CoolPlot reads it",
                  spec.get("diagram") == "ph" and spec.get("fluid", {}).get("components") == ["R134a"]
                  and len(arrays.get("states", [])) == 4, json.dumps(spec)[:300])
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

    await more_fluids(page)
    await humid_air(page)
    await theme_lines(page)

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


def select_js(index: int, value: str) -> str:
    """Choose `value` in the index-th select of the controls, as a user does (React sees the change)."""
    return ("(() => { const s = document.querySelectorAll('.diagram-controls select')[" + str(index) + "];"
            " Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set.call(s, " + json.dumps(value) + ");"
            " s.dispatchEvent(new Event('change', {bubbles: true})); })()")


def theme_js(value: str) -> str:
    """Choose `value` ('' for "as the GUI", or a CoolPlot theme) in the Theme select of the controls."""
    return ("(() => { const s = [...document.querySelectorAll('.diagram-controls label')]"
            ".find(l => l.querySelector('span')?.textContent === 'Theme')?.querySelector('select');"
            " Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set.call(s, " + json.dumps(value) + ");"
            " s.dispatchEvent(new Event('change', {bubbles: true})); })()")


async def theme_lines(page: Page):
    """The light lines belong to the tab's own look ("as the GUI"); CoolPlot's themes keep their own lines."""
    print("\n=== Themes and lines ===")
    check("no Line weight control: the light lines belong to the look of the tab",
          not await page.eval("[...document.querySelectorAll('.diagram-controls span')].some(e => e.textContent === 'Line weight')"))
    count = len(page.scene_requests)
    await page.eval(theme_js("default"))
    end = time.time() + 20
    while time.time() < end and not (len(page.scene_requests) > count
                                     and page.scene_requests[-1].get("spec", {}).get("settings", {}).get("look") == "default"):
        await asyncio.sleep(0.1)
    spec = page.scene_requests[-1].get("spec", {})
    check("CoolPlot's own theme draws its own, thicker lines: no style overrides",
          spec.get("settings", {}).get("look") == "default" and "style_overrides" not in spec, json.dumps(spec)[:300])
    check("the isolines of CoolPlot's theme are opaque",
          await page.wait_for(f"(() => {{ const p = {SVG}?.querySelector('g.cp-role-isoline path'); "
                              "return !!p && p.getAttribute('stroke-opacity') === '1'; })()", 20))
    count = len(page.scene_requests)
    await page.eval(theme_js(""))
    end = time.time() + 20
    while time.time() < end and not (len(page.scene_requests) > count
                                     and "style_overrides" in page.scene_requests[-1].get("spec", {})):
        await asyncio.sleep(0.1)
    spec = page.scene_requests[-1].get("spec", {})
    check("back to the look of the tab: the light lines again",
          all(o == {"width_pt": 0.5, "alpha": 0.55} for o in spec.get("style_overrides", {}).values())
          and spec.get("style_overrides"), json.dumps(spec)[:300])


async def solve_model(page: Page, code: str, fluid: str):
    """Solve `code` in the editor; wait until `fluid` is offered among the model's fluids."""
    await page.set_model(code)
    await page.eval("[...document.querySelectorAll('button')].find(b => b.textContent.trim() === 'Solve').click()")
    return await page.wait_for("[...document.querySelectorAll('.diagram-controls select optgroup[label=\"Model fluids\"] option')]"
                               f".some(o => o.value === {json.dumps(fluid)})", 60)


async def last_scene(page: Page, fluid: str, timeout: float = 30.0) -> dict:
    """The last scene request for `fluid`, once one has gone."""
    end = time.time() + timeout
    while time.time() < end:
        found = [r for r in page.scene_requests if r.get("fluid") == fluid]
        if found:
            return found[-1]
        await asyncio.sleep(0.1)
    return {}


def row_select_js(row: int, index: int, value: str) -> str:
    """Choose `value` in the index-th select of the row-th arrays overlay, as a user does."""
    return ("(() => { const s = document.querySelectorAll('.diagram-arrays')[" + str(row) + "].querySelectorAll('select')["
            + str(index) + "]; Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set.call(s, "
            + json.dumps(value) + "); s.dispatchEvent(new Event('change', {bubbles: true})); })()")


# The labels of the arrays overlays' fields (the axis symbols) and the arrays chosen, row by row
ROWS_JS = ("[...document.querySelectorAll('.diagram-arrays')].map(r => [...r.querySelectorAll('select')]"
           ".map(s => s.parentElement.textContent.split(':')[0] + '=' + s.value).join(' '))")
ADD_ARRAYS_JS = "[...document.querySelectorAll('.diagram-overlays button')].find(b => b.textContent.trim() === 'arrays').click()"
REMOVE_ARRAYS_JS = "document.querySelector('.diagram-arrays button[title=\"Remove these arrays\"]').click()"


async def arrays_cycles(page: Page, test, timeout: float = 20.0) -> list:
    """The arrays cycles of the last scene request, once `test` holds for them (or the timeout ends)."""
    end = time.time() + timeout
    while True:
        cycles = [c for c in (page.scene_requests[-1].get("cycles", []) if page.scene_requests else []) if "arrays" in c]
        if test(cycles) or time.time() > end:
            return cycles
        await asyncio.sleep(0.1)


async def arrays_overlays(page: Page):
    """Two arrays overlays: default arrays of the axis quantities, paths along the processes, another diagram type."""
    print("\n=== Arrays overlays: two cycles of R134a in arrays ===")
    await page.set_model(ARRAYS_MODEL)
    await page.eval("[...document.querySelectorAll('button')].find(b => b.textContent.trim() === 'Solve').click()")
    check("the arrays model solves: '+ arrays' is offered",
          await page.wait_for("[...document.querySelectorAll('.diagram-overlays button')].some(b => b.textContent.trim() === 'arrays' && !b.disabled)", 60))
    for _ in range(2):
        await page.eval(ADD_ARRAYS_JS)
        await asyncio.sleep(0.3)
    cycles = await arrays_cycles(page, lambda c: len(c) == 2)
    check("two overlays, each with arrays of the axis quantities by default (h and p: from their names, from the property"
          " calls), the second ones not taken by the first",
          [c.get("arrays") for c in cycles] == [{"h": "h", "p": "P"}, {"h": "hb", "p": "Pb"}], json.dumps(cycles)[:400])
    others = [c for c in page.scene_requests[-1].get("cycles", []) if "arrays" not in c]
    check("the first arrays overlay switched the inferred states off (they would number the same states twice)",
          others == [], json.dumps(others)[:300])
    check("named after their arrays, paths along the processes, an open chain, no loops written",
          [c.get("name") for c in cycles] == ["P, h", "Pb, hb"] and all(c.get("paths") == "auto" and c.get("closed") is False
                                                                     and "loops" not in c for c in cycles), json.dumps(cycles)[:400])
    rows = await page.eval(ROWS_JS)
    check("one row per overlay: the horizontal then the vertical axis symbol", rows == ["h=h p=P", "h=hb p=Pb"], str(rows))
    await page.eval("document.querySelectorAll('.diagram-arrays input[type=\"checkbox\"]').forEach(c => c.click())")
    cycles = await arrays_cycles(page, lambda c: len(c) == 2 and all(x.get("closed") for x in c))
    check("close loop", len(cycles) == 2 and all(c.get("closed") is True for c in cycles), json.dumps(cycles)[:300])
    check("both cycles are drawn", await page.wait_for(
        f"(() => {{ const items = [...{SVG}?.querySelectorAll('[data-item]') ?? []].map(e => e.dataset.item);"
        " return items.includes('cycle/P, h/loop/0') && items.includes('cycle/Pb, hb/loop/0'); })()", 60))
    await page.idle()
    points = await page.eval(f"[...{SVG}.querySelectorAll('.cp-role-cycle path')]"
                             ".map(p => (p.getAttribute('d').match(/[ML]/g) || []).length).reduce((a, b) => a + b, 0)")
    check("auto paths follow the processes (more points than two closed loops of straight segments)", (points or 0) > 2 * 5,
          f"{points} points")
    await page.screenshot("diagram_arrays_ph.png")

    print("\n=== Another diagram type: the overlays follow the axes ===")
    await page.eval(select_js(1, "Ts"))
    cycles = await arrays_cycles(page, lambda c: [x.get("name") for x in c] == ["T, s", "Tb, sb"])
    check("on T-s the overlays take the arrays of s and T", [c.get("arrays") for c in cycles] == [{"s": "s", "T": "T"}, {"s": "sb", "T": "Tb"}]
          and all(c.get("closed") is True for c in cycles), json.dumps(cycles)[:400])
    rows = await page.eval(ROWS_JS)
    check("the rows name s and T", rows == ["s=s T=T", "s=sb T=Tb"], str(rows))
    await page.idle(1.0)
    await page.screenshot("diagram_arrays_ts.png")
    await page.eval(select_js(1, "ph"))
    cycles = await arrays_cycles(page, lambda c: [x.get("name") for x in c] == ["P, h", "Pb, hb"])
    check("back on p-h: the arrays of h and p again", [c.get("arrays") for c in cycles] == [{"h": "h", "p": "P"}, {"h": "hb", "p": "Pb"}],
          json.dumps(cycles)[:400])
    # An array chosen by hand
    await page.eval(row_select_js(1, 1, "P"))
    cycles = await arrays_cycles(page, lambda c: len(c) == 2 and c[1].get("arrays", {}).get("p") == "P")
    check("an array chosen by hand; the names stay unique", [c.get("name") for c in cycles] == ["P, h", "P, hb"], json.dumps(cycles)[:300])
    await page.eval(row_select_js(1, 1, "Pb"))
    await arrays_cycles(page, lambda c: len(c) == 2 and c[1].get("arrays", {}).get("p") == "Pb")
    check("no error shown", not await page.eval("!!document.querySelector('.diagram-status .diagram-error')"),
          await page.eval("document.querySelector('.diagram-status')?.textContent"))
    await page.idle(1.0)


async def more_fluids(page: Page):
    """Phase 3: an ideal gas with its solved cycle, an incompressible solution at a concentration, a blend."""
    print("\n=== An ideal gas: a solved Brayton cycle of N2 ===")
    check("the Brayton model solves and N2 is a model fluid", await solve_model(page, BRAYTON_MODEL, "N2"))
    await page.eval(select_js(0, "N2"))
    req = await last_scene(page, "N2")
    check("a chart of the ideal gas is asked for", req.get("fluid") == "N2" and "fractions" not in req, json.dumps(req)[:200])
    names = await page.eval("[...document.querySelectorAll('.diagram-controls select')[1].options].map(o => o.value)")
    check("its diagram types: all but T-h", names == ["ph", "Ts", "hs", "ps", "prho", "Trho", "pT"], str(names))
    labels = await page.eval("document.querySelector('.diagram-controls').textContent")
    check("no critical point, saturation scale or region preset; a T and p range instead",
          "critical point" not in labels and "saturation scale" not in labels and "Region" not in labels
          and await page.eval("!!document.querySelector('.diagram-range')"), labels[:300])
    families = await page.eval("[...document.querySelectorAll('.diagram-family')].map(f => f.textContent.trim())")
    check("no lines of constant quality", not any(f.startswith("x") for f in families or []), str(families))
    check("no quality family in the spec", all(i["key"] != "Q" for i in req.get("spec", {}).get("isolines", [])),
          json.dumps(req.get("spec", {}).get("isolines")))
    # The cycle of the model's arrays on a T-s chart: one overlay (the other removed), mapped to the arrays of s and T
    await page.eval(REMOVE_ARRAYS_JS)
    await asyncio.sleep(0.3)
    check("an overlay removed", await page.wait_for("document.querySelectorAll('.diagram-arrays').length === 1", 10))
    await page.eval(select_js(1, "Ts"))
    cycles = await arrays_cycles(page, lambda c: [x.get("arrays") for x in c] == [{"s": "s", "T": "T"}])
    check("on the T-s chart of N2 the overlay takes the arrays of s and T", [c.get("arrays") for c in cycles] == [{"s": "s", "T": "T"}]
          and page.scene_requests[-1].get("fluid") == "N2", json.dumps(cycles)[:300])
    check("the Brayton cycle is drawn on the T-s chart of N2",
          await page.wait_for(f"!!{SVG}?.querySelector('[data-item=\"cycle/T, s/loop/0\"]')", 60))
    await page.idle(1.0)
    check("no error shown", not await page.eval("!!document.querySelector('.diagram-status .diagram-error')"),
          await page.eval("document.querySelector('.diagram-status')?.textContent"))
    # A T range set by hand: absolute limits in the spec
    count = len(page.scene_requests)
    await page.eval("(() => { const i = document.querySelectorAll('.diagram-range input')[1]; i.focus();"
                    " Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set.call(i, '1200');"
                    " i.dispatchEvent(new Event('input', {bubbles: true})); i.blur(); })()")
    end = time.time() + 15
    while time.time() < end and not (len(page.scene_requests) > count and isinstance(page.scene_requests[-1]["spec"].get("limits"), dict)):
        await asyncio.sleep(0.1)
    limits = page.scene_requests[-1]["spec"].get("limits") if page.scene_requests else None
    check("the highest temperature typed gives absolute limits (°C to K)", isinstance(limits, dict)
          and abs(limits.get("T_max", 0) - 1473.15) < 1e-6 and limits.get("p_min") == {"relative": 1.01}, json.dumps(limits))
    await page.idle(1.0)
    await page.screenshot("diagram_ideal_gas_ts.png")

    print("\n=== An incompressible solution at a concentration ===")
    await page.eval(select_js(0, "MEG"))
    req = await last_scene(page, "MEG")
    check("the solution is drawn at a mass fraction (the middle of its range: 0.3)", req.get("fractions") == [0.3],
          json.dumps(req)[:200])
    check("a concentration field shows it", await page.wait_for("document.querySelector('.diagram-fraction input')?.value === '0.3'", 10))
    names = await page.eval("[...document.querySelectorAll('.diagram-controls select')[1].options].map(o => o.value)")
    check("its diagram types: the pressure charts", names == ["ph", "ps", "prho", "pT"], str(names))
    check("its lines: T, s (no density, no quality), rounded", sorted(i["key"] for i in req.get("spec", {}).get("isolines", []))
          == ["T", "s"] and all(i.get("rounding") for i in req["spec"]["isolines"]), json.dumps(req.get("spec", {}).get("isolines")))
    check("the chart shows its isotherms", await page.wait_for(f"!!{SVG}?.querySelector('[data-item^=\"iso/T/\"]')", 60))
    count = len(page.scene_requests)
    await page.eval("(() => { const i = document.querySelector('.diagram-fraction input');"
                    " Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set.call(i, '0.4');"
                    " i.dispatchEvent(new Event('input', {bubbles: true})); })()")
    end = time.time() + 15
    while time.time() < end and not (len(page.scene_requests) > count and page.scene_requests[-1].get("fractions") == [0.4]):
        await asyncio.sleep(0.1)
    check("another concentration asks for another chart", page.scene_requests[-1].get("fractions") == [0.4])
    await page.idle(1.0)
    check("no error shown", not await page.eval("!!document.querySelector('.diagram-status .diagram-error')"),
          await page.eval("document.querySelector('.diagram-status')?.textContent"))
    await page.screenshot("diagram_incompressible.png")

    print("\n=== A blend of the model ===")
    check("the blend model solves", await solve_model(page, BLEND_MODEL, "R32[0.7]&R125[0.3]"))
    await page.eval(select_js(0, "R32[0.7]&R125[0.3]"))
    req = await last_scene(page, "R32[0.7]&R125[0.3]")
    # A mixture opens on the ACHP region (its kind's default_region): DEF reaches the triple point of the lightest
    # component, where the lines of a blend are slow to compute
    check("the blend is drawn with its dome, on the ACHP region", req.get("spec", {}).get("limits") == "ACHP"
          and any(i["key"] == "Q" for i in req.get("spec", {}).get("isolines", []))
          and await page.wait_for(f"!!{SVG}?.querySelector('[data-item=\"iso/Q/0.5\"]')", 60), json.dumps(req)[:300])
    await page.idle(1.0)
    await page.screenshot("diagram_blend.png")


async def humid_air(page: Page):
    """Humid air: a solved air-conditioning process on the psychrometric chart, the hover, the Mollier chart."""
    print("\n=== Humid air: a solved air-conditioning process ===")
    check("the air-conditioning model solves and AirH2O is a model fluid", await solve_model(page, AIR_CONDITIONING_MODEL, "AirH2O"))
    await page.eval(select_js(0, "AirH2O"))
    req = await last_scene(page, "AirH2O")
    spec = req.get("spec", {})
    check("a psychrometric chart of humid air at the model's total pressure (1 bar)", spec.get("diagram") == "WT"
          and req.get("p_Pa") == 1e5 and set(spec.get("limits", {}).get("bounds", {})) == {"T", "W"}, json.dumps(req)[:300])
    names = await page.eval("[...document.querySelectorAll('.diagram-controls select')[1].options].map(o => o.value)")
    check("its diagram types: the psychrometric and Mollier charts", names == ["WT", "ThxW"], str(names))
    labels = await page.eval("document.querySelector('.diagram-controls').textContent")
    check("no critical point, saturation scale or region preset; the total pressure in Pa (the model's units) instead",
          "critical point" not in labels and "saturation scale" not in labels and "Region" not in labels
          and await page.eval("document.querySelector('.diagram-pressure')?.value") == "100000", labels[:300])
    check("its lines: φ, Twb, h, v, round and labelled", sorted(i["key"] for i in spec.get("isolines", [])) == ["RH", "Twb", "h", "v"]
          and all(i.get("rounding") and i.get("labels") for i in spec["isolines"]), json.dumps(spec.get("isolines")))
    check("the chart shows the saturation line (iso/RH/1)", await page.wait_for(f"!!{SVG}?.querySelector('[data-item=\"iso/RH/1\"]')", 60))
    # The process through the elements of the arrays: a new overlay takes those of T and W (w, from the states
    # of humid air), an open chain
    while await page.eval("!!document.querySelector('.diagram-arrays')"):
        await page.eval(REMOVE_ARRAYS_JS)
        await asyncio.sleep(0.3)
    await page.eval(ADD_ARRAYS_JS)
    cycles = await arrays_cycles(page, lambda c: len(c) == 1 and len(c[0].get("arrays", {})) == 2)
    check("a new overlay of the psychrometric chart takes the arrays of T and of the humidity ratio, an open chain",
          [c.get("arrays") for c in cycles] == [{"T": "T", "W": "w"}] and cycles[0].get("closed") is False
          and page.scene_requests[-1].get("fluid") == "AirH2O", json.dumps(cycles)[:300])
    check("the process is drawn through its three states",
          await page.wait_for(f"!!{SVG}?.querySelector('.cp-role-cycle path') && "
                              f"{SVG}.querySelectorAll('.cp-role-state_number').length >= 3", 60))
    await page.idle(1.0)
    check("no error shown", not await page.eval("!!document.querySelector('.diagram-status .diagram-error')"),
          await page.eval("document.querySelector('.diagram-status')?.textContent"))
    await page.screenshot("diagram_psychrometric.png")
    # The hover: the state of humid air under the cursor
    box = await page.eval(f"(() => {{ const r = {SVG}.querySelector('rect.cp-background').getBoundingClientRect();"
                          " return [r.left, r.top, r.width, r.height]; })()")
    count = len(page.point_requests)
    await page.mouse("mouseMoved", box[0] + 0.6 * box[2], box[1] + 0.85 * box[3], 0)
    check("hovering the chart shows T, W, φ, T_wb, h and v",
          await page.wait_for("(() => { const t = document.querySelector('.diagram-tooltip[data-kind=\"point\"]');"
                              " return !!t && ['T = ', 'W = ', 'φ = ', 'T_wb = ', 'h = ', 'v = '].every(k => t.textContent.includes(k)); })()", 15)
          and len(page.point_requests) > count, await page.eval("document.querySelector('.diagram-tooltip')?.textContent"))
    body = page.point_requests[-1] if page.point_requests else {}
    check("the state request names humid air at its total pressure", body.get("fluid") == "AirH2O" and body.get("p_Pa") == 1e5
          and body.get("diagram") == "WT", json.dumps(body))
    await page.screenshot("diagram_psychrometric_hover.png")
    await page.mouse("mouseMoved", 5, 5, 0)
    # The Mollier h,x chart: the same region in its temperature scale, isotherms
    await page.eval(select_js(1, "ThxW"))
    check("the Mollier chart shows its isotherms", await page.wait_for(f"!!{SVG}?.querySelector('[data-item^=\"iso/T/\"]')", 60))
    await page.idle(1.0)
    req = await last_scene(page, "AirH2O")
    check("the Mollier chart is bounded in Thx and W", set(req.get("spec", {}).get("limits", {}).get("bounds", {})) == {"Thx", "W"}
          and any(i["key"] == "T" for i in req["spec"].get("isolines", [])), json.dumps(req.get("spec", {}))[:300])
    check("the process is drawn on the Mollier chart", await page.wait_for(f"!!{SVG}?.querySelector('.cp-role-cycle path')", 30))
    rows = await page.eval(ROWS_JS)
    check("the overlay keeps its arrays, named as the Mollier chart writes them (x, T)", rows == ["x=w T=T"], str(rows))
    await page.screenshot("diagram_mollier.png")


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
