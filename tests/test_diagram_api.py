#!/usr/bin/env python3
"""
test_diagram_api.py — Integration tests of the Diagrams endpoints of the GUI server (docs/coolplot_integration.md §3).

Tests GET /api/v1/diagram/options, POST /api/v1/diagram/scene (a default p-h chart of R134a, the badge naming the
equation of state and the versions, a chart of a solved model with its cycle from the model's variables and the paths
between its states, the errors), POST /api/v1/diagram/export (the chart as SVG and PDF), the diagrams of the LaTeX
report as vector PDF figures (when pdflatex is installed), POST
/api/v1/diagram/state (the state under the cursor), GET /api/v1/diagram/states (the states of the solved model),
GET/PUT /api/v1/diagram/definitions, the definitions
in the ZIP bundle (diagrams.json), their reset by New and their return by Back, and that the former route
POST /api/v1/coolprop/saturation is gone. Phase 3: incompressible solutions at a mass fraction, ideal gases, a tabular
backend set in coolsolve.conf, blends, humid air (psychrometric and Mollier charts at a total pressure, a solved
air-conditioning process, the state under the cursor; docs/coolplot_integration.md §5).

Usage:
    # Start CoolSolve with --gui, then:
    python3 test_diagram_api.py [base_url]

    # Or run with auto-start (from the build/ directory) on a free port:
    python3 test_diagram_api.py --auto

Exit code 0 if all tests pass, 1 otherwise.
"""

import io
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zipfile
import zlib
from typing import Optional

BASE_URL = "http://localhost:8550"
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv
EXAMPLE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "examples", "refrigeration1.eescode"))
SESSION = {"X-CoolSolve-Session": "e" * 32}  # one tab of its own

passed = 0
failed = 0
errors: list[str] = []
_opener = urllib.request.build_opener()


def raw(method: str, path: str, body=None, headers=None, data: Optional[bytes] = None):
    """Send a request; return (status, content_type, body bytes). HTTP errors are returned, not raised."""
    if data is None and body is not None:
        data = json.dumps(body).encode()
    hdrs = {"Content-Type": "application/json", **SESSION}
    hdrs.update(headers or {})
    req = urllib.request.Request(BASE_URL + path, data=data, headers=hdrs, method=method)
    try:
        resp = _opener.open(req, timeout=60)
        return resp.status, resp.headers.get("Content-Type", ""), resp.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers.get("Content-Type", ""), e.read()


def raw_headers(method: str, path: str, body=None, headers=None):
    """Send a JSON request; return (status, response headers, body bytes). HTTP errors are returned, not raised."""
    hdrs = {"Content-Type": "application/json", **SESSION}
    hdrs.update(headers or {})
    req = urllib.request.Request(BASE_URL + path, data=json.dumps(body).encode() if body is not None else None,
                                 headers=hdrs, method=method)
    try:
        resp = _opener.open(req, timeout=120)
        return resp.status, resp.headers, resp.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers, e.read()


def api(method: str, path: str, body=None, data: Optional[bytes] = None):
    """Call /api/v1; return (status, parsed JSON or None)."""
    status, _, content = raw(method, "/api/v1" + path, body, None, data)
    try:
        return status, json.loads(content) if content else None
    except ValueError:
        return status, None


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


def default_spec(diagram="ph", families=("T", "s", "rho", "Q")):
    return {"diagram": diagram, "units": "engineering", "limits": "DEF", "show_saturation": True,
            "isolines": [{"key": k, "num": 11 if k == "Q" else 15} for k in families],
            "settings": {"look": "default", "branding": {"references": "none", "versions": False}}}


# ── Test groups ──────────────────────────────────────────────────────────────

def test_options():
    print("\n── Options ──")
    status, options = api("GET", "/diagram/options")
    check("options return 200", status == 200, str(status))
    names = [t["name"] for t in options["diagram_types"]]
    check("the fluid diagram types, T-h included, and the charts of humid air",
          names == ["ph", "Ts", "hs", "ps", "prho", "Trho", "pT", "Th", "WT", "ThxW"], str(names))
    check("each type lists its families with a default count",
          all(f["default_num"] > 0 for t in options["diagram_types"] for f in t["families"]))
    units = {u["name"]: u for u in options["unit_systems"]}
    check("unit systems: the presets and coolsolve", list(units) == ["SI", "KSI", "engineering", "IP", "coolsolve"], str(list(units)))
    record = units["coolsolve"]["record"]
    check("the coolsolve record is a custom unit record (°C, Pa, J/kg)",
          record["name"] == "coolsolve" and record["units"]["T"] == ["deg C", 1.0, -273.15, "K"]
          and record["units"]["p"] == ["Pa", 1.0, 0.0] and record["units"]["h"] == ["J/kg", 1.0, 0.0], json.dumps(record))
    check("region presets", [r["name"] for r in options["regions"]] == ["DEF", "ACHP", "ORC", "NONE"])
    check("themes and label modes", options["themes"] == ["default", "dark", "print", "coolpack"]
          and options["label_modes"] == ["rails", "edge", "margin", "inline"])
    check("coolplot-core version and commit", options["coolplot"]["core_version"] == "1.5.0"
          and len(options["coolplot"]["commit"]) >= 7, str(options["coolplot"]))
    return options


def test_scene(options):
    print("\n── Scene: default p-h chart of R134a ──")
    body = {"fluid": "R134a", "spec": default_spec(), "cycles": [], "figure_size_pt": [600, 420]}
    status, scene = api("POST", "/diagram/scene", body)
    check("scene returns 200", status == 200, str(scene)[:300])
    if status != 200:
        return
    check("an SVG document", scene["svg"].startswith("<svg") and 'data-item="iso/Q/0"' in scene["svg"])
    axes = scene["axes"]
    check("axes: keys, limits, log flags, units", axes["x_key"] == "h" and axes["y_key"] == "p" and axes["y_log"] is True
          and len(axes["x_limits"]) == 2 and axes["x_unit"]["scale"] == 1e-3 and axes["y_unit"]["label"] == "bar", json.dumps(axes))
    check("the spec as composed has the fluid and the cycles",
          scene["spec"]["fluid"]["components"] == ["R134a"] and scene["spec"]["cycles"] == [])
    check("diagnostics of every curve", len(scene["diagnostics"]["curves"]) == 15 * 3 + 11
          and all("item_id" in c and "n_failed" in c for c in scene["diagnostics"]["curves"]))
    check("no warning", scene["diagnostics"]["warnings"] == [], str(scene["diagnostics"]["warnings"]))
    check("timings", set(scene["timing_ms"]) == {"curves", "compose", "svg"})

    # Labels, units of the model, a zoom window, the labels sent back
    spec = default_spec()
    spec["isolines"][0]["labels"] = True
    spec["units"] = next(u for u in options["unit_systems"] if u["name"] == "coolsolve")["record"]
    spec["view_SI"] = [[2.0e5, 4.5e5], [1.0e5, 3.0e6]]
    status, zoomed = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": spec})
    check("a zoomed chart in the model's units", status == 200 and zoomed["axes"]["x_unit"]["label"] == "J/kg"
          and abs(zoomed["axes"]["x_limits"][0] - 2.0e5) < 1.0, str(zoomed.get("axes") if zoomed else status))
    if status == 200:
        check("isoline labels returned", len(zoomed["labels"]) > 0 and all(t["role"] == "isoline_label" for t in zoomed["labels"]))
        status, again = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": spec, "previous_labels": zoomed["labels"]})
        check("the labels sent back keep their places", status == 200 and again["labels"] == zoomed["labels"])


def test_badge():
    print("\n── Badge: CoolPlot's defaults ──")
    spec = default_spec()
    del spec["settings"]["branding"]
    status, scene = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": spec, "figure_size_pt": [600, 420]})
    check("a chart with CoolPlot's default badge", status == 200, str(scene)[:300])
    if status != 200:
        return
    check("the badge names the equation of state (CoolPlot's short reference)",
          "R134a: Tillner-Roth and Baehr (1994)" in scene["svg"], scene["svg"][-600:])
    check("the badge names CoolProp and CoolPlot", "CoolProp 8." in scene["svg"] and "CoolPlot 1.5.0" in scene["svg"])
    check("no warning about the badge", scene["diagnostics"]["warnings"] == [], str(scene["diagnostics"]["warnings"]))
    spec["settings"]["branding"] = {"references": "full", "versions": False}
    status, scene = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": spec})
    check("full references", status == 200 and "Tillner-Roth and Baehr, J. Phys. Chem. Ref. Data 23 (1994) 657-729" in scene["svg"]
          and "CoolPlot 1.5.0" not in scene["svg"], str(status))
    spec["settings"]["branding"] = {"references": "long"}
    status, result = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": spec})
    check("an unknown references style: 400", status == 400 and "Unknown references style" in (result or {}).get("error", ""),
          f"{status} {result}")


def test_export():
    print("\n── Export: SVG, PDF ──")
    body = {"fluid": "R134a", "spec": default_spec(), "cycles": [], "figure_size_pt": [600, 420]}
    status, content_type, svg = raw("POST", "/api/v1/diagram/export", {**body, "format": "svg"})
    check("SVG export returns 200 with image/svg+xml", status == 200 and content_type.startswith("image/svg+xml"),
          f"{status} {content_type} {svg[:200]!r}")
    check("the SVG of the scene", svg.startswith(b"<svg") and b'data-item="iso/Q/0"' in svg and b'width="800"' in svg)
    status, headers, pdf = raw_headers("POST", "/api/v1/diagram/export", {**body, "format": "pdf"})
    check("PDF export returns 200 with application/pdf", status == 200 and headers.get("Content-Type", "").startswith("application/pdf"),
          f"{status} {pdf[:300]!r}")
    check("a PDF with an embedded TrueType font, the figure size as its page",
          pdf.startswith(b"%PDF-") and b"/FontFile2" in pdf and b"/MediaBox [0 0 600 420]" in pdf, repr(pdf[:100]))
    disposition = headers.get("Content-Disposition", "")
    check("a file name for the download", disposition.startswith("attachment;") and disposition.endswith('_ph.pdf"'), disposition)
    status, result = api("POST", "/diagram/export", {**body, "format": "png"})
    check("an unknown format: 400", status == 400 and "svg" in (result or {}).get("error", ""), f"{status} {result}")
    status, result = api("POST", "/diagram/export", {**body, "fluid": "Unobtainium", "format": "svg"})
    check("an unknown fluid: 400", status == 400 and "Unknown fluid" in (result or {}).get("error", ""), f"{status} {result}")
    # The session of the tab names the file: the model of this session, not of another one
    api("PUT", "/model-name", {"modelName": "my model"})
    _, headers_, _ = raw_headers("POST", "/api/v1/diagram/export", {**body, "format": "svg"})
    check("the file is named after the model of the tab's session", 'filename="my_model_ph.svg"' in headers_.get("Content-Disposition", ""),
          headers_.get("Content-Disposition", ""))
    _, headers_, _ = raw_headers("POST", "/api/v1/diagram/export", {**body, "format": "svg"}, {"X-CoolSolve-Session": "f" * 32})
    check("another session names it after its own model", 'filename="diagram_ph.svg"' in headers_.get("Content-Disposition", ""),
          headers_.get("Content-Disposition", ""))


def test_report_figures(definitions_with_cycle):
    """The LaTeX report of a solved model takes its diagrams as vector PDF (needs pdflatex on the machine)."""
    print("\n── Report: the diagrams as PDF figures ──")
    if shutil.which("pdflatex") is None:
        print("  (pdflatex not found: skipped)")
        return
    api("PUT", "/diagram/definitions", definitions_with_cycle)
    status, report = api("GET", "/latex/report")
    check("a LaTeX report of the solve", status == 200 and report.get("available") is True, str(report)[:200])
    if not report.get("available"):
        return
    check("the report includes thermo_diagrams.pdf under \\IfFileExists",
          "\\IfFileExists{thermo_diagrams.pdf}" in report["content"])
    status, headers, pdf = raw_headers("POST", "/api/v1/latex/compile", {"plots": []})
    check("the report compiles", status == 200 and pdf.startswith(b"%PDF-"), f"{status} {pdf[:300]!r}")
    # The figures are coolplot-core's PDFs, with their fonts: pdfTeX copies their font subsets ("COOLPL+...") into
    # its compressed object streams, once per figure
    fonts = 0
    for stream in re.findall(rb"stream\r?\n(.*?)\r?\nendstream", pdf, re.S):
        try:
            fonts += zlib.decompress(stream).count(b"/BaseFont /COOLPL+")
        except zlib.error:
            pass
    check("both diagrams are the server's vector figures", fonts >= 2 * 2, f"{fonts} font dictionaries of coolplot-core")


def test_errors():
    print("\n── Errors ──")
    cases = [
        ("unknown fluid", {"fluid": "Unobtainium", "spec": default_spec()}, "Unknown fluid"),
        ("no fluid", {"spec": default_spec()}, "fluid"),
        ("a p-h chart of humid air", {"fluid": "AirH2O", "spec": default_spec()}, "cannot show humid air"),
        ("unknown diagram", {"fluid": "R134a", "spec": {**default_spec(), "diagram": "xy"}}, "Unknown diagram type"),
        ("bad theme", {"fluid": "R134a", "spec": {**default_spec(), "settings": {"look": "neon"}}}, "Unknown theme"),
        ("bad figure size", {"fluid": "R134a", "spec": default_spec(), "figure_size_pt": [1, 2]}, "figure size"),
    ]
    for name, body, message in cases:
        status, result = api("POST", "/diagram/scene", body)
        check(f"{name}: 400 with an error", status == 400 and result is not None and message in result.get("error", ""),
              f"{status} {result}")
    status, result = api("POST", "/diagram/scene", data=b"{not json")
    check("invalid JSON: 400 with an error", status == 400 and "error" in (result or {}), f"{status} {result}")
    status, _, _ = raw("POST", "/api/v1/coolprop/saturation", {"fluid": "R134a"})
    check("the former saturation route is gone", status == 404, str(status))


def test_definitions():
    print("\n── Definitions ──")
    api("POST", "/new")
    status, doc = api("GET", "/diagram/definitions")
    check("an empty document by default", status == 200 and doc == {"format": "coolsolve-diagrams", "version": 1, "diagrams": []},
          str(doc))
    definitions = {"format": "coolsolve-diagrams", "version": 1, "diagrams": [
        {"name": "Cycle on log(p)-h", "fluid": "R134a", "spec": default_spec(),
         "cycles": [{"name": "Heat pump", "states": []}]}]}
    status, result = api("PUT", "/diagram/definitions", definitions)
    check("PUT stores the document", status == 200 and result.get("success") is True, str(result))
    status, doc = api("GET", "/diagram/definitions")
    check("GET returns it", doc == definitions, str(doc)[:200])
    for name, bad in (("another format", {"format": "other"}), ("a version that is no number", {"version": "1"}),
                      ("diagrams that are no list", {"diagrams": {}}), ("a list", [1, 2])):
        status, result = api("PUT", "/diagram/definitions", bad)
        check(f"PUT refuses {name} with a JSON error", status == 400 and "error" in (result or {}), f"{status} {result}")
    status, result = api("PUT", "/diagram/definitions", data=b"[1,")
    check("PUT refuses invalid JSON with a JSON error", status == 400 and "error" in (result or {}), f"{status} {result}")
    status, doc = api("GET", "/diagram/definitions")
    check("a refused PUT changes nothing", doc == definitions)
    return definitions


def test_bundle(definitions):
    print("\n── ZIP bundle, New and Back ──")
    api("PUT", "/files/eescode", {"content": "x = 1\n"})
    api("PUT", "/model-name", {"modelName": "diagrams"})
    api("PUT", "/diagram/definitions", definitions)
    status, _, bundle = raw("GET", "/api/v1/files/bundle")
    check("the bundle downloads", status == 200 and bundle[:2] == b"PK", str(status))
    if bundle[:2] != b"PK":
        return
    names = zipfile.ZipFile(io.BytesIO(bundle)).namelist()
    check("the bundle holds diagrams.json", "diagrams.json" in names, str(names))
    stored = json.loads(zipfile.ZipFile(io.BytesIO(bundle)).read("diagrams.json"))
    check("diagrams.json is the document", stored == definitions)

    api("POST", "/new")
    _, doc = api("GET", "/diagram/definitions")
    check("New clears the definitions", doc["diagrams"] == [], str(doc))
    status, _ = api("POST", "/back")
    _, doc = api("GET", "/diagram/definitions")
    check("Back brings them back", status == 200 and doc == definitions, str(doc)[:200])

    api("POST", "/new")
    boundary = "----coolsolvediagramtest"
    body = (f'--{boundary}\r\nContent-Disposition: form-data; name="file"; filename="diagrams.zip"\r\n'
            "Content-Type: application/zip\r\n\r\n").encode() + bundle + f"\r\n--{boundary}--\r\n".encode()
    status, _, content = raw("POST", "/api/v1/files/upload", None, {"Content-Type": f"multipart/form-data; boundary={boundary}"},
                             body)
    result = json.loads(content)
    check("the bundle opens", status == 200 and "diagrams.json" in result.get("files", []), str(result))
    _, doc = api("GET", "/diagram/definitions")
    check("opening the bundle restores the definitions", doc == definitions, str(doc)[:200])


def test_solved_cycle():
    print("\n── A solved model and its cycle ──")
    api("POST", "/new")
    with open(EXAMPLE, encoding="utf-8") as f:
        api("PUT", "/files/eescode", {"content": f.read()})
    status, _ = api("POST", "/solve", {})
    result = None
    for _ in range(120):
        time.sleep(0.25)
        status, result = api("GET", "/solve/result")
        if status == 200 and result.get("status") not in (None, "", "Running"):
            break
    check("the model solves", status == 200 and result.get("status") == "Success", str(result)[:200])
    if status != 200:
        return
    states = [{"number": i, "vars": {"p": f"P_{i}", "h": f"h_{i}"}} for i in range(1, 5)]
    states[0]["vars"]["T"] = "T_1"
    states[0]["label"] = "su,cp"
    cycles = [{"name": "Refrigeration", "states": states, "paths": "straight"},
              {"name": "Missing", "states": [{"number": 1, "vars": {"p": "P_99", "h": "h_1"}},
                                             {"number": 2, "vars": {"p": "P_2", "h": "h_2"}}]},
              {"name": "Markers", "states": states, "loops": []}]
    status, scene = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": default_spec(), "cycles": cycles})
    check("the scene of the solved model", status == 200, str(scene)[:300])
    if status != 200:
        return
    drawn = [s for s in scene["states"] if s["cycle"] == 0]
    check("the four states of the cycle", [s["number"] for s in drawn] == [1, 2, 3, 4], str(drawn)[:300])
    h1 = result["variables"]["h_1"]
    check("the state holds the solved enthalpy (SI)", abs(drawn[0]["values_SI"]["h"] - h1) <= 1e-9 * abs(h1),
          f"{drawn[0]['values_SI']['h']} vs {h1}")
    check("the solved temperature in K", abs(drawn[0]["values_SI"]["T"] - (result["variables"]["T_1"] + 273.15)) < 1e-9)
    check("full states", set(drawn[2]["values_SI"]) == {"T", "p", "h", "s", "rho", "v", "u", "Q"}, str(drawn[2]))
    check("the label and the variables of a state", drawn[0]["label"] == "su,cp" and drawn[0]["vars"] == {"p": "P_1", "h": "h_1",
                                                                                                         "T": "T_1"})
    check("the cycle is drawn", 'data-item="cycle/Refrigeration/loop/0"' in scene["svg"])
    check("markers alone for loops []", 'data-item="cycle/Markers/states"' in scene["svg"]
          and 'cycle/Markers/loop' not in scene["svg"])
    cycle = scene["spec"]["cycles"][0]
    check("the spec writes a CoolPlot cycle with straight paths",
          cycle["loops"] == [[1, 2, 3, 4, 1]] and all(kind == "straight" for _, kind in cycle["paths"])
          and cycle["labels"] == [[1, "su,cp"]] and "state" in cycle["states"][0][1], json.dumps(cycle)[:300])
    warnings = scene["diagnostics"]["warnings"]
    check("a cycle with a missing variable is left out with a warning",
          [c["name"] for c in scene["spec"]["cycles"]] == ["Refrigeration", "Markers"]
          and any("Missing" in w and "P_99" in w for w in warnings), str(warnings))

    # The paths between the states (phase 2): "auto" by default, CoolPlot's record, a kind for every connection
    auto = {"name": "Auto", "states": states}
    record = {"name": "Record", "states": states, "paths": [[[1, 2], "isentropic"]], "steps": 10}
    status, scene = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": default_spec(), "cycles": [auto, record]})
    check("cycles with paths", status == 200, str(scene)[:300])
    if status == 200:
        cycles = scene["spec"]["cycles"]
        check("auto paths: CoolPlot's default, nothing written", cycles[0]["paths"] == [] and cycles[0]["steps"] == 20,
              json.dumps(cycles[0]["paths"]))
        check("a record of paths and steps pass through", cycles[1]["paths"] == [[[1, 2], "isentropic"]] and cycles[1]["steps"] == 10)
        check("the paths are drawn", 'data-item="cycle/Auto/loop/0"' in scene["svg"])
        check("no warning about paths", not any("path" in w for w in scene["diagnostics"]["warnings"]),
              str(scene["diagnostics"]["warnings"]))
    status, result = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": default_spec(),
                                                    "cycles": [{"name": "Bad", "states": states, "paths": "zigzag"}]})
    check("an unknown path kind: 400", status == 400 and "zigzag" in result.get("error", ""), f"{status} {result}")

    # The states of the solved model from its property calls (§4.2)
    status, model_states = api("GET", "/diagram/states")
    check("states of the solved model", status == 200 and "R134a" in model_states.get("fluids", {}), str(model_states)[:300])
    if status == 200 and "R134a" in model_states.get("fluids", {}):
        found = model_states["fluids"]["R134a"]
        by_number = {s["number"]: s for s in found}
        check("unique state numbers", len(by_number) == len(found))
        # h_1 = enthalpy(R134a, T=T_1, P=P_1) and s_1 = entropy(R134a, T=T_1, P=P_1) share T_1 and P_1: two calls
        # without quality at the same (T, p) are one single-phase state (the rule as decided, plan §4.2)
        check("state 1 holds P_1, T_1, h_1 and s_1 (calls at the same T and p merge)",
              by_number.get(1, {}).get("vars") == {"p": "P_1", "T": "T_1", "h": "h_1", "s": "s_1"}
              and by_number[1]["label"] is None, str(by_number.get(1)))
        # h_2_s = enthalpy(R134a, P=P_2, s=s_1) shares only P_2 with state 2 and only s_1 with state 1: a state of
        # its own (the isentropic outlet of the compressor)
        check("the isentropic outlet is a state of its own",
              any(s["vars"] == {"p": "P_2", "s": "s_1", "h": "h_2_s"} for s in found), str(found)[:400])
        check("state 2: the outlet of the compressor (P_2, h_2, T_2)", by_number.get(2, {}).get("vars") == {"p": "P_2", "h": "h_2", "T": "T_2"},
              str(by_number.get(2)))
        check("every state has two independent properties",
              all(len(set(s["vars"]) - {"T", "p"}) >= 1 for s in found), str(found)[:300])
        # The states as GET /diagram/states writes them (label null) are a valid cycle definition
        status, drawn = api("POST", "/diagram/scene", {"fluid": "R134a", "spec": default_spec(),
                                                       "cycles": [{"name": "Model states", "states": found, "loops": []}]})
        check("the states of the model, labels null included, draw as a cycle",
              status == 200 and len([s for s in drawn.get("states", []) if s["cycle"] == 0]) == len(found), str(drawn)[:300])

    spec = default_spec()
    del spec["settings"]["branding"]
    test_report_figures({"format": "coolsolve-diagrams", "version": 1, "diagrams": [
        {"name": "Refrigeration", "fluid": "R134a", "spec": spec, "cycles": [{"name": "Cycle", "states": states}]},
        {"name": "T-s", "fluid": "R134a", "spec": {**default_spec("Ts", ("p", "h", "Q"))}, "cycles": []}]})


BRAYTON_MODEL = """T_1 = 25
P_1 = 1e5
P_2 = 10*P_1
h_1 = enthalpy(N2, T=T_1)
s_1 = entropy(N2, T=T_1, P=P_1)
s_2 = s_1
h_2 = enthalpy(N2, P=P_2, s=s_2)
T_2 = temperature(N2, P=P_2, s=s_2)
P_3 = P_2
T_3 = 1100
h_3 = enthalpy(N2, T=T_3)
s_3 = entropy(N2, T=T_3, P=P_3)
P_4 = P_1
s_4 = s_3
h_4 = enthalpy(N2, P=P_4, s=s_4)
T_4 = temperature(N2, P=P_4, s=s_4)
"""

GLYCOL_MODEL = """c = 30
T_1 = 5
T_2 = 15
h_1 = enthalpy(MEG, T=T_1, C=c)
h_2 = enthalpy(MEG, T=T_2, C=c)
q = h_2 - h_1
"""


def solve(source: str):
    """Solve `source` in the session; the result of GET /solve/result."""
    api("POST", "/new")
    api("PUT", "/files/eescode", {"content": source})
    api("POST", "/solve", {})
    result = None
    for _ in range(240):
        time.sleep(0.25)
        status, result = api("GET", "/solve/result")
        if status == 200 and result.get("status") not in (None, "", "Running"):
            break
    return result or {}


def test_more_fluids(options):
    """Phase 3: incompressible solutions, ideal gases, tabular backends, blends computed on several threads."""
    print("\n── More fluids (phase 3) ──")
    check("the fluid types drawn", options["fluid_types"] == ["Real", "Mixture", "Incompressible", "IdealGas", "HumidAir"],
          str(options["fluid_types"]))
    kinds = options.get("fluid_kinds", {})
    check("what each kind offers", kinds.get("Incompressible", {}).get("diagram_types") == ["ph", "ps", "prho", "pT"]
          and kinds["Incompressible"]["fractions"] is True and kinds["IdealGas"]["saturation"] is False
          and "Th" not in kinds["IdealGas"]["diagram_types"] and kinds["IdealGas"]["default_box"] == [200, 2000, 1e4, 1e7],
          json.dumps(kinds)[:300])
    check("the default region of a new chart: DEF for pure fluids, ACHP for mixtures",
          kinds["Real"].get("default_region") == "DEF" and kinds["Mixture"].get("default_region") == "ACHP",
          json.dumps(kinds)[:300])
    _, fluids = api("GET", "/coolprop/fluids")
    meg = next((f for f in fluids["fluids"] if f["name"] == "MEG"), {})
    check("a solution lists its range of mass fractions", meg.get("fractionRange") == [0, 0.6], str(meg))

    # An incompressible solution at a mass fraction
    spec = default_spec("ph", ("T", "s"))
    status, scene = api("POST", "/diagram/scene", {"fluid": "MEG", "fractions": [0.3], "basis": "mass", "spec": spec})
    check("an incompressible chart", status == 200 and 'data-item="iso/T/' in scene.get("svg", ""), str(scene)[:300])
    if status == 200:
        check("its absolute limits and CoolPlot's fluid record in the spec",
              scene["spec"]["fluid"] == {"backend": "INCOMP", "components": ["MEG"], "fractions": [0.3], "basis": "mass",
                                         "interactions": []}
              and scene["spec"]["limits"]["p_min"] == 1e5 and scene["spec"]["limits"]["p_max"] == 1e7,
              json.dumps(scene["spec"]["fluid"]) + json.dumps(scene["spec"]["limits"]))
    status, result = api("POST", "/diagram/state", {"fluid": "MEG", "fractions": [0.3], "diagram": "ph", "x_SI": 5e4, "y_SI": 2e5})
    point = (result or {}).get("state") or {}
    check("the state under the cursor of a solution", status == 200 and point.get("p") == 2e5 and point.get("T", 0) > 250
          and point.get("Q") is None, str(result))
    for name, body, message in (
            ("a solution without its concentration", {"fluid": "MEG", "spec": spec}, "give its concentration"),
            ("a concentration out of range", {"fluid": "MEG", "fractions": [0.9], "spec": spec}, "not 0.9"),
            ("fractions that are no numbers", {"fluid": "MEG", "fractions": ["a"], "spec": spec}, "list of numbers"),
            ("a basis other than mass", {"fluid": "MEG", "fractions": [0.3], "basis": "volu", "spec": spec}, "mass fraction"),
            ("fractions of a pure fluid", {"fluid": "R134a", "fractions": [0.3], "spec": spec}, "takes no fractions"),
            ("an unknown incompressible", {"fluid": "NoSuchBrine", "fractions": [0.3], "spec": spec}, "Unknown fluid"),
            ("an ideal gas on a T-h chart", {"fluid": "N2", "spec": default_spec("Th", ("p", "s"))}, "no T-h chart")):
        status, result = api("POST", "/diagram/scene", body)
        check(f"{name}: 400 with an error", status == 400 and message in (result or {}).get("error", ""), f"{status} {result}")
    status, result = api("POST", "/diagram/state", {"fluid": "MEG", "diagram": "ph", "x_SI": 5e4, "y_SI": 2e5})
    check("the state of a solution without its concentration: 400", status == 400, f"{status} {result}")

    # A glycol loop: the concentration of the model's calls, its states at the pressure of the solver
    result = solve(GLYCOL_MODEL)
    check("the glycol model solves", result.get("status") == "Success", str(result)[:200])
    status, states = api("GET", "/diagram/states")
    check("the concentration of the model's solution", status == 200 and states.get("fractions", {}).get("MEG") == [0.3],
          str(states)[:300])
    meg_states = states.get("fluids", {}).get("MEG", [])
    cycles = [{"name": "Loop", "states": [{"number": s["number"], "vars": s["vars"]} for s in meg_states], "loops": []}]
    status, scene = api("POST", "/diagram/scene", {"fluid": "MEG", "fractions": [0.3], "spec": spec, "cycles": cycles})
    drawn = scene.get("states", []) if status == 200 else []
    check("the states of the glycol loop, at 1 atm", len(drawn) == 2 and all(s["values_SI"]["p"] == 101325 for s in drawn),
          str(scene)[:300])

    # An ideal gas: a solved Brayton cycle of N2 on its T-s chart (EES/JANAF reference states)
    result = solve(BRAYTON_MODEL)
    check("the Brayton model solves", result.get("status") == "Success", str(result)[:200])
    if result.get("status") == "Success":
        states = [{"number": i, "vars": {"T": f"T_{i}", "p": f"P_{i}", "h": f"h_{i}", "s": f"s_{i}"}} for i in range(1, 5)]
        spec = default_spec("Ts", ("p", "rho"))
        status, scene = api("POST", "/diagram/scene", {"fluid": "N2", "spec": spec, "cycles": [{"name": "Brayton", "states": states}]})
        check("an ideal-gas T-s chart with its cycle", status == 200 and 'data-item="cycle/Brayton/loop/0"' in scene.get("svg", "")
              and scene["diagnostics"]["warnings"] == [], str(scene)[:300])
        if status == 200:
            h1 = result["variables"]["h_1"]
            check("N2 on the JANAF basis: h(25 °C) is its enthalpy of formation, 0",
                  abs(h1) < 1e-6 and scene["states"][0]["values_SI"]["h"] == h1, str(scene["states"][0]))
            check("the lines of CoolSolve's provider", all(c["method"] == "ideal_gas" for c in scene["diagnostics"]["curves"]))
        status, result = api("POST", "/diagram/state", {"fluid": "N2", "diagram": "Ts", "x_SI": 7000, "y_SI": 400})
        check("the state under the cursor of an ideal gas", status == 200 and abs(result["state"]["T"] - 400) < 1e-9
              and result["state"]["Q"] is None, str(result))

    # A tabular backend set in coolsolve.conf: the lines of constant quality from the equation of state
    api("POST", "/new")
    api("PUT", "/files/conf", {"content": "coolpropBackend = BICUBIC&HEOS\n"})
    status, scene = api("POST", "/diagram/scene", {"fluid": "Water", "spec": default_spec("Ts", ("p", "Q"))})
    check("a chart with the model's tabular backend", status == 200 and scene["spec"]["fluid"]["backend"] == "BICUBIC&HEOS"
          and 'data-item="iso/Q/0.5"' in scene["svg"], str(scene)[:300])
    api("PUT", "/files/conf", {"content": ""})

    # Blends: the lines of one chart on several threads
    status, scene = api("POST", "/diagram/scene", {"fluid": "R32[0.7]&R125[0.3]", "spec": default_spec()})
    check("a blend chart", status == 200 and 'data-item="iso/Q/0.5"' in scene.get("svg", ""), str(scene)[:200])
    if status == 200:
        print(f"  R32[0.7]&R125[0.3] p-h cold: {scene['timing_ms']['curves']:.0f} ms of lines")
    status, opened = api("POST", "/library/open", {"id": "CSL-0324"})
    if status == 200 and opened.get("success"):
        start = time.time()
        status, scene = api("POST", "/diagram/scene", {"fluid": "Propane[0.5]&IsoButane[0.5]", "spec": default_spec(),
                                                       "figure_size_pt": [600, 420]})
        check("the default p-h chart of CSL-0324 (propane/isobutane)", status == 200, str(scene)[:200])
        print(f"  CSL-0324 Propane[0.5]&IsoButane[0.5] p-h cold: {(time.time() - start) * 1000:.0f} ms over HTTP")


AIR_CONDITIONING_MODEL = """P = 101325
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


def humid_air_spec(diagram="WT"):
    families = ("T", "RH", "Twb", "h", "v") if diagram == "ThxW" else ("RH", "Twb", "h", "v")
    return {"diagram": diagram, "units": "engineering", "show_saturation": True,
            "limits": {"bounds": {"Thx" if diagram == "ThxW" else "T": [263.15, 323.15], "W": [0, 0.03]}},
            "isolines": [{"key": k, "num": 10, "rounding": True, "labels": True} for k in families],
            "settings": {"look": "default", "branding": {"references": "none", "versions": False}}}


def test_humid_air(options):
    """Humid air (CS-FEAT-PSYCHRO): the psychrometric and Mollier charts at a total pressure, a solved process."""
    print("\n── Humid air ──")
    kind = options.get("fluid_kinds", {}).get("HumidAir", {})
    check("what humid air offers", kind.get("diagram_types") == ["WT", "ThxW"] and kind.get("chart_kind") == "humid_air"
          and kind.get("total_pressure") == 101325 and kind.get("saturation") is False, json.dumps(kind))
    _, fluids = api("GET", "/coolprop/fluids")
    listed = [f for f in fluids["fluids"] if f["type"] == "HumidAir"]
    check("humid air is listed as AirH2O, once", [f["name"] for f in listed] == ["AirH2O"], str(listed))
    names = [f["name"] for f in fluids["fluids"]]
    check("every fluid is listed once", len(names) == len(set(names)), str([n for n in names if names.count(n) > 1]))
    for diagram, x_key in (("WT", "T"), ("ThxW", "W")):
        status, scene = api("POST", "/diagram/scene", {"fluid": "AirH2O", "spec": humid_air_spec(diagram)})
        check(f"a {diagram} chart of humid air", status == 200 and 'data-item="iso/RH/1"' in scene.get("svg", "")
              and scene["axes"]["x_key"] == x_key and scene["diagnostics"]["warnings"] == [], str(scene)[:300])
        if status == 200:
            check(f"{diagram}: CoolPlot's humid-air record in the spec", scene["spec"]["fluid"] == {"model": "humid_air", "p_Pa": 101325}
                  and "W" in scene["spec"]["limits"]["bounds"], json.dumps(scene["spec"]["fluid"]))
            print(f"  {diagram} cold: {scene['timing_ms']['curves']:.0f} ms of lines")
    status, scene = api("POST", "/diagram/scene", {"fluid": "AirH2O", "p_Pa": 9e4, "spec": humid_air_spec()})
    check("a chart at another total pressure", status == 200 and scene["spec"]["fluid"]["p_Pa"] == 9e4, str(scene)[:200])
    for name, body, message in (
            ("a total pressure out of range", {"fluid": "AirH2O", "p_Pa": 1e3, "spec": humid_air_spec()}, "between 10000.0"),
            ("a total pressure that is no number", {"fluid": "AirH2O", "p_Pa": "1 atm", "spec": humid_air_spec()}, "p_Pa"),
            ("a total pressure of R134a", {"fluid": "R134a", "p_Pa": 1e5, "spec": default_spec()}, "takes no total pressure"),
            ("a psychrometric chart of R134a", {"fluid": "R134a", "spec": humid_air_spec()}, "cannot show R134a"),
            ("a region preset on a psychrometric chart", {"fluid": "AirH2O", "spec": {**humid_air_spec(), "limits": "DEF"}},
             "needs a Region"),
            ("a region in the fog", {"fluid": "AirH2O", "spec": {**humid_air_spec(), "limits": {"bounds": {"T": [263.15, 273.15],
                                                                                                         "W": [0.02, 0.03]}}}},
             "above saturation")):
        status, result = api("POST", "/diagram/scene", body)
        check(f"{name}: 400 with an error", status == 400 and message in (result or {}).get("error", ""), f"{status} {result}")

    # The state under the cursor: T, W, RH, Twb, Tdp, h, v, Thx, p; nothing above saturation
    status, result = api("POST", "/diagram/state", {"fluid": "AirH2O", "diagram": "WT", "x_SI": 298.15, "y_SI": 0.01})
    point = (result or {}).get("state") or {}
    check("the state of humid air under the cursor", status == 200 and set(point) == {"T", "W", "RH", "Twb", "Tdp", "h", "v", "Thx", "p"}
          and point["T"] == 298.15 and point["W"] == 0.01 and 0 < point["RH"] < 1, str(result))
    status, result = api("POST", "/diagram/state", {"fluid": "AirH2O", "diagram": "WT", "x_SI": 298.15, "y_SI": 0.021})
    check("fog has no state", status == 200 and result == {"state": None}, str(result))
    status, result = api("POST", "/diagram/state", {"fluid": "AirH2O", "p_Pa": 9e4, "diagram": "WT", "x_SI": 298.15, "y_SI": 0.021})
    check("at 0.9 bar the same point is humid air", status == 200 and (result.get("state") or {}).get("p") == 9e4, str(result))
    status, result = api("POST", "/diagram/state", {"fluid": "AirH2O", "diagram": "ThxW", "x_SI": 0.01, "y_SI": 300.0})
    check("a point of the Mollier chart", status == 200 and abs(((result or {}).get("state") or {}).get("Thx", 0) - 300) < 1e-6, str(result))
    start = time.time()
    for i in range(20):
        api("POST", "/diagram/state", {"fluid": "AirH2O", "diagram": "WT", "x_SI": 290 + i, "y_SI": 0.005})
    print(f"  state of humid air under the cursor: {(time.time() - start) / 20 * 1000:.2f} ms per request (HTTP)")

    # A solved air-conditioning process: its states grouped from the humid-air calls, drawn with straight paths
    result = solve(AIR_CONDITIONING_MODEL)
    check("the air-conditioning model solves", result.get("status") == "Success", str(result)[:200])
    status, states = api("GET", "/diagram/states")
    air = states.get("fluids", {}).get("AirH2O", [])
    check("its three states, by the humid-air calls they share",
          [s["number"] for s in air] == [1, 2, 3] and air[0]["vars"] == {"T": "T[1]", "p": "P", "RH": "RH[1]", "W": "w[1]", "h": "h[1]"},
          str(states)[:400])
    cycle = {"name": "Air conditioning", "states": [{"number": s["number"], "vars": s["vars"]} for s in air], "loops": [[1, 2, 3]]}
    status, scene = api("POST", "/diagram/scene", {"fluid": "AirH2O", "spec": humid_air_spec(), "cycles": [cycle]})
    drawn = scene.get("states", []) if status == 200 else []
    check("the process on the psychrometric chart", len(drawn) == 3 and 'data-item="cycle/Air conditioning/loop/0"' in scene.get("svg", "")
          and scene["diagnostics"]["warnings"] == [], str(scene)[:300])
    if len(drawn) == 3:
        check("the states where the solver put them", drawn[0]["values_SI"]["W"] == result["variables"]["w[1]"]
              and drawn[2]["values_SI"]["W"] == drawn[1]["values_SI"]["W"], str(drawn[0]))
    status, scene = api("POST", "/diagram/scene", {"fluid": "AirH2O", "p_Pa": 1e5, "spec": humid_air_spec(), "cycles": [cycle]})
    check("at another total pressure the states are left out with a warning",
          status == 200 and scene["states"] == [] and any("101325 Pa" in w for w in scene["diagnostics"]["warnings"]), str(scene)[:300])


def test_state():
    print("\n── State under the cursor ──")
    status, result = api("POST", "/diagram/state", {"fluid": "R134a", "diagram": "ph", "x_SI": 3.0e5, "y_SI": 2.0e5})
    state = (result or {}).get("state") or {}
    check("a point inside the dome", status == 200 and 0.0 < state.get("Q", -1) < 1.0, str(result))
    check("the axis values come back", state.get("h") == 3.0e5 and state.get("p") == 2.0e5, str(state))
    check("a full state", set(state) == {"T", "p", "h", "s", "rho", "u", "v", "Q"}, str(state))
    status, result = api("POST", "/diagram/state", {"fluid": "R134a", "diagram": "ph", "x_SI": 4.5e5, "y_SI": 2.0e5})
    check("single phase: Q = -1", status == 200 and result["state"]["Q"] == -1.0, str(result))
    status, result = api("POST", "/diagram/state", {"fluid": "R134a", "diagram": "ph", "x_SI": 3.0e5, "y_SI": 1e12})
    check("outside the fluid's range: null", status == 200 and result == {"state": None}, str(result))
    for name, body, message in (
            ("unknown fluid", {"fluid": "Unobtainium", "diagram": "ph", "x_SI": 1, "y_SI": 1}, "Unknown fluid"),
            ("unknown diagram", {"fluid": "R134a", "diagram": "xy", "x_SI": 1, "y_SI": 1}, "Unknown diagram type"),
            ("no axis value", {"fluid": "R134a", "diagram": "ph", "x_SI": 3e5}, "y_SI"),
            ("an axis value that is no number", {"fluid": "R134a", "diagram": "ph", "x_SI": "a", "y_SI": 1}, "x_SI")):
        status, result = api("POST", "/diagram/state", body)
        check(f"{name}: 400 with an error", status == 400 and message in (result or {}).get("error", ""), f"{status} {result}")
    start = time.time()
    for i in range(50):
        api("POST", "/diagram/state", {"fluid": "R134a", "diagram": "ph", "x_SI": 2.0e5 + 1000 * i, "y_SI": 2.0e5})
    print(f"  state under the cursor: {(time.time() - start) / 50 * 1000:.2f} ms per request (HTTP)")


# ── Runner ───────────────────────────────────────────────────────────────────

def free_port() -> int:
    with socket.socket() as s:
        s.bind(("localhost", 0))
        return s.getsockname()[1]


def start_server(binary: str, port: int):
    """Start `coolsolve --gui <port>` (BASE_URL must already point to that port); wait until it answers."""
    print(f"Starting CoolSolve from {binary} on port {port}...")
    proc = subprocess.Popen([binary, "--gui", str(port), "--no-browser"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(40):
        time.sleep(0.25)
        try:
            api("GET", "/health")
            return proc
        except Exception:
            pass
    proc.terminate()
    print("ERROR: Server did not start in 10s")
    sys.exit(1)


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
        port = free_port()
        BASE_URL = f"http://localhost:{port}"
        server_proc = start_server(binary, port)
    try:
        options = test_options()
        test_scene(options)
        test_badge()
        test_export()
        test_errors()
        definitions = test_definitions()
        test_bundle(definitions)
        status, empty = api("GET", "/diagram/states")
        check("no states without a solve", status == 200 and empty == {"fluids": {}}, str(empty))
        test_state()
        test_solved_cycle()
        test_more_fluids(options)
        test_humid_air(options)
    except Exception as e:
        print(f"\nFATAL ERROR: {e!r}")
        check("no fatal error", False, repr(e))
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
