#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write the CoolPlot fixtures of the [diagram] tests (tests/fixtures/coolplot/).

Modelled on CoolPlot's cpp/tests/data/make_data.py (compose_case): CoolPlot (Python, at the pinned commit
COOLSOLVE_COOLPLOT_PINNED_TAG of CMakeLists.txt) draws a few chart specs, and for each the compose record
(the input of coolplot-core: spec, domain, curves) and the scene it composes are written. The C++ tests
(tests/test_diagram.cpp) read them without Python: the curve requests must be the same, the domain and the
curves within a tolerance (the CoolProp wheel is not the CoolProp commit CoolSolve is pinned to), and the
scene must have the same items, groups, roles, styles and texts (docs/coolplot_integration.md §6.1).

Run it in a virtual environment with CoolPlot at the pinned commit:

    python3 -m venv venv
    venv/bin/pip install -e /path/to/CoolPlot          # pulls coolprop>=8.0 from PyPI
    venv/bin/python -I tools/coolplot_fixtures.py [--out tests/fixtures/coolplot]

header.json records the CoolPlot version and commit, the CoolProp version of the wheel, whether CoolPlot used
CoolProp's C++ port of its numerics (CoolProp.CoolProp.Plot) or its Python numerics, and the date.

To keep the files small the charts use CoolPlot's "interactive" calculation settings (60 points per line, 100
per quality line), the curves of the compose record keep every third point and the last one ("stride",
"n_points") with 8 significant digits, and the lines of the scene are written without their points (the tests
compare the curves of the compose record instead). "input" is therefore the compose record with thinned curves.
"""
import argparse
import datetime
import json
import os
import subprocess
import sys

import CoolProp
from coolplot import PropertyDiagram, __version__ as coolplot_version
from coolplot.compose.builder import compose
from coolplot.engine import plan
from coolplot.scene import Scene
from coolplot.thermo import native

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.join(HERE, "..", "tests", "fixtures", "coolplot")

# The fixtures compare the lines and the cycles, so their charts have no badge: its lines name the CoolProp of the
# wheel, not the pinned one (the badge lines are tested against CoolPlot's text in tests/test_diagram_outputs.cpp)
NO_BADGE = {"references": "none", "versions": False}


STRIDE = 3


def rounded(value, digits=8):
    """A number with `digits` significant digits (None stays None)."""
    return None if value is None else float("{0:.{1}g}".format(value, digits))


def curve_record(request, curve):
    record = curve.to_dict()
    n = len(record["x"])
    kept = sorted(set(range(0, n, STRIDE)) | ({n - 1} if n else set()))
    record["x"] = [rounded(record["x"][i]) for i in kept]
    record["y"] = [rounded(record["y"][i]) for i in kept]
    record["stride"] = STRIDE
    record["n_points"] = n
    return {"request": {"item_id": request.item_id, "key": request.key, "value_SI": request.value_SI,
                        "points": request.points, "labels": request.labels, "anchor_SI": request.anchor_SI},
            "curve": record}


def light_scene(scene):
    """The scene record in drawing order, the points of its lines replaced by their number and the styles of the
    items written once each, in "styles", the items giving the index of theirs."""
    record = scene.to_dict()
    styles = []
    for item in record["items"]:
        if item.get("type") == "line":
            item["n_points"] = len(item.pop("x"))
            item.pop("y")
        if item["style"] not in styles:
            styles.append(item["style"])
        item["style"] = styles.index(item["style"])
    record["styles"] = styles
    return record


def compose_case(out, name, chart):
    """The compose record of a chart and the scene the Python composer makes of it (make_data.compose_case)."""
    chart.settings = chart.settings.with_(branding=NO_BADGE)
    spec = chart.spec
    domain = chart.engine.domain(spec)
    curves = chart.engine.curves(spec, domain, plan(spec, domain))
    cycles = [chart.engine.cycle(spec, c) for c in spec.cycles]
    sources = chart.engine.source_lines(spec)
    assert not sources, "the fixtures have no badge"
    items, axes = compose(spec, domain, curves, None, None, cycles, sources)
    scene = Scene(axes, items)
    data = {
        "spec": spec.to_dict(),
        "domain": {"tp_box": list(domain.tp_box), "ranges_SI": {k: list(v) for k, v in domain.ranges_SI}},
        "curves": [curve_record(r, c) for r, c in curves],
        "cycles": [{"states": {str(n): st.to_dict() for n, st in states.items()},
                    "chains": [[st.to_dict() for st in chain] for chain in chains]} for states, chains in cycles],
        "plot_size_pt": None,
        "previous_labels": None,
    }
    if spec.marked_points and domain.points:
        data["domain"]["points"] = {n: dict(values) for n, values in domain.points}
    if spec.saturation_scale is not None and domain.scale:
        data["domain"]["scale"] = [list(pair) for pair in domain.scale]
    if spec.settings.labels.dome and domain.dome:
        data["domain"]["dome"] = [[value, dict(bubble), dict(dew)] for value, bubble, dew in domain.dome]
    # The requests of every curve, also those left out for want of valid points, as plan() gives them
    requests = [{"item_id": r.item_id, "key": r.key, "value_SI": r.value_SI, "points": r.points, "labels": r.labels,
                 "anchor_SI": r.anchor_SI} for r in plan(spec, domain)]
    with open(os.path.join(out, name + ".json"), "w", encoding="ascii", newline="\n") as f:
        json.dump({"input": data, "requests": requests, "scene": light_scene(scene)}, f, sort_keys=True,
                  allow_nan=False, separators=(",", ":"))
        f.write("\n")
    print("{0}: {1} curves, {2} items".format(name, len(curves), len(items)))


def cases(out):
    # The default p-h chart of R134a (no golden file upstream)
    chart = PropertyDiagram("R134a", "ph", settings="interactive")
    chart.set_default_lines()
    compose_case(out, "r134a_ph_default", chart)

    # Water T-s with labels
    chart = PropertyDiagram("Water", "Ts", settings="interactive")
    chart.set_default_lines(labels=True)
    compose_case(out, "water_ts", chart)

    # R134a T-h (no golden file upstream), SI units
    chart = PropertyDiagram("R134a", "Th", units="SI", settings="interactive")
    chart.set_default_lines(labels=True)
    compose_case(out, "r134a_th", chart)

    # R290 p-h, ACHP region, engineering units, round values, isentropes through the saturated vapour every 10 K
    chart = PropertyDiagram("R290", "ph", tp_limits="ACHP", units="engineering", settings="interactive")
    chart.set_default_lines(rounding=True, labels=True)
    chart.set_lines("s", at_saturation={"T": 10}, labels=True)
    compose_case(out, "r290_ph_achp", chart)

    # An R134a refrigeration cycle on p-h, its states given by p and h (SI) and joined by "auto" paths: a polytropic
    # compression, isobars through the dome with the saturation states, an isenthalp (CoolPlot's cycle_paths)
    chart = PropertyDiagram("R134a", "ph", units="SI", settings="interactive")
    chart.set_lines("T", num=5, labels=True)
    chart.set_lines("Q", num=3)
    chart.add_cycle("Refrigeration", [{"p": 2.0e5, "h": 4.0e5}, {"p": 1.0e6, "h": 4.4e5}, {"p": 1.0e6, "h": 2.4e5},
                                      {"p": 2.0e5, "h": 2.4e5}], labels={1: "su,cp"})
    compose_case(out, "r134a_ph_cycle", chart)

    # CO2 p-h with the critical and triple points marked, a saturation scale and labels inside the dome
    chart = PropertyDiagram("CO2", "ph", units="KSI", settings="interactive")
    chart.set_default_lines(labels=True)
    chart.marked_points = ("critical", "triple")
    chart.saturation_scale = True
    chart.settings = chart.settings.with_(labels={"dome": True})
    compose_case(out, "co2_ph_marked", chart)


def coolplot_commit():
    """The git commit of the CoolPlot checkout that is installed (editable install), or "unknown"."""
    import coolplot
    try:
        return subprocess.run(["git", "-C", os.path.dirname(coolplot.__file__), "rev-parse", "HEAD"], check=True,
                              capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", default=DEFAULT_OUT, help="the folder of the fixtures")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cases(args.out)
    header = {"coolplot_version": coolplot_version, "coolplot_commit": coolplot_commit(),
              "coolprop_version": CoolProp.__version__, "coolprop_revision": getattr(CoolProp, "__gitrevision__", ""),
              "coolprop_port_used": native.enabled(), "python": sys.version.split()[0],
              "date": datetime.date.today().isoformat(), "generator": "tools/coolplot_fixtures.py"}
    with open(os.path.join(args.out, "header.json"), "w", encoding="ascii", newline="\n") as f:
        json.dump(header, f, indent=2, sort_keys=True)
        f.write("\n")
    print(json.dumps(header, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
