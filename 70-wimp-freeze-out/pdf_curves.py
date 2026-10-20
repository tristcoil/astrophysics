"""Read curves exactly from a matplotlib vector figure (the published limits of Ackermann et al. 2015).

The Fermi-LAT collaboration's six-year dwarf limits are not in a table, but the figures in the
paper's arXiv source are matplotlib PDFs: every curve is a path of exact coordinates. We convert a
page to SVG with poppler's ``pdftocairo``, take the path coordinates, and calibrate both log axes
on the major tick marks, each identified by the tick label next to it (``pdftotext -bbox`` gives
the label words and their boxes). Curves are recognised by their drawing style (colour, width,
dash). The figure's own thermal-relic line, drawn from Steigman et al. (2012)'s table, checks the
calibration (tests).
"""

from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np

_NUM = r"[-\d.]+"


def _svg_paths(pdf: Path) -> list[dict]:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "fig.svg"
        subprocess.run(["pdftocairo", "-svg", str(pdf), str(out)], check=True, capture_output=True)
        s = out.read_text()
    body = s[s.find("</defs>"):]
    paths = []
    for attrs in re.findall(r"<path ([^>]*)/>", body):
        d = re.search(r' d="([^"]*)"', attrs)
        if not d:
            continue
        pts = [(float(x), float(y)) for x, y in re.findall(rf"[ML] ({_NUM}) ({_NUM})", d.group(1))]
        tr = re.search(r'transform="matrix\(([^)]*)\)"', attrs)
        m = [float(v) for v in tr.group(1).split(",")] if tr else [1, 0, 0, 1, 0, 0]
        P = np.array([(m[0] * x + m[2] * y + m[4], m[1] * x + m[3] * y + m[5]) for x, y in pts])

        def attr(name):
            a = re.search(rf' {name}="([^"]*)"', attrs)
            return a.group(1) if a else None

        paths.append({"xy": P, "stroke": attr("stroke"), "width": attr("stroke-width"),
                      "dash": attr("stroke-dasharray"), "fill": attr("fill"),
                      "opacity": attr("stroke-opacity")})
    return paths


def _words(pdf: Path) -> list[dict]:
    r = subprocess.run(["pdftotext", "-bbox", str(pdf), "-"], check=True, capture_output=True, text=True)
    out = []
    for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>',
                         r.stdout):
        out.append({"x0": float(m[1]), "y0": float(m[2]), "x1": float(m[3]), "y1": float(m[4]),
                    "text": m[5].replace("−", "-")})
    return out


def _decade_labels(words):
    """'10' followed immediately by an integer exponent word -> (centre x, centre y, exponent)."""
    labs = []
    for i, w in enumerate(words):
        if w["text"] != "10":
            continue
        for v in words:
            if v is w or not re.fullmatch(r"-?\d+", v["text"]):
                continue
            # the exponent starts where '10' ends and is raised a little above it
            if abs(v["x0"] - w["x1"]) < 1.0 and v["y0"] < w["y0"] and abs(v["y0"] - w["y0"]) < 6.0:
                cx = 0.5 * (w["x0"] + v["x1"])
                cy = 0.5 * (w["y0"] + w["y1"])
                labs.append((cx, cy, int(v["text"])))
                break
    return labs


def read_figure(pdf: Path) -> dict:
    """Calibrated curves of a log-log matplotlib figure: {'curves': [...], 'calibration': {...}}."""
    paths = _svg_paths(pdf)
    words = _words(pdf)
    segs = [p for p in paths if len(p["xy"]) == 2]
    # the axes frame: the longest horizontal and vertical segments
    horiz = [s for s in segs if abs(s["xy"][0, 1] - s["xy"][1, 1]) < 1e-6]
    vert = [s for s in segs if abs(s["xy"][0, 0] - s["xy"][1, 0]) < 1e-6]
    hlen = max(abs(s["xy"][1, 0] - s["xy"][0, 0]) for s in horiz)
    vlen = max(abs(s["xy"][1, 1] - s["xy"][0, 1]) for s in vert)
    frame_h = [s for s in horiz if abs(abs(s["xy"][1, 0] - s["xy"][0, 0]) - hlen) < 1e-6]
    frame_v = [s for s in vert if abs(abs(s["xy"][1, 1] - s["xy"][0, 1]) - vlen) < 1e-6]
    y_top, y_bot = sorted({round(float(s["xy"][0, 1]), 4) for s in frame_h})[:: len(frame_h) and 1][0], \
        sorted({round(float(s["xy"][0, 1]), 4) for s in frame_h})[-1]
    x_left = min(float(s["xy"][0, 0]) for s in frame_v)
    x_right = max(float(s["xy"][0, 0]) for s in frame_v)
    # major ticks: the longest short segments perpendicular to the bottom and left axes
    tick_v = [s for s in vert if abs(s["xy"][0, 1] - y_bot) < 1e-3 or abs(s["xy"][1, 1] - y_bot) < 1e-3]
    tick_v = [s for s in tick_v if abs(s["xy"][1, 1] - s["xy"][0, 1]) < 10]
    lmax = max(abs(s["xy"][1, 1] - s["xy"][0, 1]) for s in tick_v)
    xmaj = sorted({round(float(s["xy"][0, 0]), 4) for s in tick_v if abs(abs(s["xy"][1, 1] - s["xy"][0, 1]) - lmax) < 1e-3})
    tick_h = [s for s in horiz if abs(s["xy"][0, 0] - x_left) < 1e-3 or abs(s["xy"][1, 0] - x_left) < 1e-3]
    tick_h = [s for s in tick_h if abs(s["xy"][1, 0] - s["xy"][0, 0]) < 10]
    lmax = max(abs(s["xy"][1, 0] - s["xy"][0, 0]) for s in tick_h)
    ymaj = sorted({round(float(s["xy"][0, 1]), 4) for s in tick_h if abs(abs(s["xy"][1, 0] - s["xy"][0, 0]) - lmax) < 1e-3})
    labs = _decade_labels(words)
    xl = [(cx, e) for cx, cy, e in labs if cy > y_bot]
    yl = [(cy, e) for cx, cy, e in labs if cx < x_left]
    xpairs = [(t, min(xl, key=lambda q: abs(q[0] - t))) for t in xmaj]
    ypairs = [(t, min(yl, key=lambda q: abs(q[0] - t))) for t in ymaj]
    xpairs = [(t, e) for t, (c, e) in xpairs if abs(c - t) < 12]
    ypairs = [(t, e) for t, (c, e) in ypairs if abs(c - t) < 12]
    ax = np.polyfit([t for t, _ in xpairs], [e for _, e in xpairs], 1)
    ay = np.polyfit([t for t, _ in ypairs], [e for _, e in ypairs], 1)
    resid = max(np.max(np.abs(np.polyval(ax, [t for t, _ in xpairs]) - [e for _, e in xpairs])),
                np.max(np.abs(np.polyval(ay, [t for t, _ in ypairs]) - [e for _, e in ypairs])))
    curves = []
    for p in paths:
        if len(p["xy"]) < 5 or p["fill"] not in (None, "none"):
            continue
        X, Y = p["xy"][:, 0], p["xy"][:, 1]
        inside = (X >= x_left - 1e-6) & (X <= x_right + 1e-6)
        curves.append({"stroke": p["stroke"], "width": p["width"], "dash": p["dash"], "opacity": p["opacity"],
                       "log10x": np.polyval(ax, X[inside]), "log10y": np.polyval(ay, Y[inside])})
    return {"curves": curves, "calibration": {"x_ticks": xpairs, "y_ticks": ypairs, "max_tick_residual": float(resid),
                                              "x_range": [float(np.polyval(ax, x_left)), float(np.polyval(ax, x_right))],
                                              "y_range": [float(np.polyval(ay, y_bot)), float(np.polyval(ay, y_top))]}}


def ackermann2015(pdf: Path) -> dict:
    """The six-year Pass 8 limit, its median expectation, the 4-year Pass 7 limit and the figure's
    thermal-relic line, as (mass GeV, <sigma v> cm^3/s)."""
    fig = read_figure(pdf)
    out = {"calibration": fig["calibration"]}
    for c in fig["curves"]:
        key = None
        if c["stroke"] == "rgb(0%, 0%, 0%)" and c["width"] == "1.5" and not c["dash"]:
            key = "limit_6yr_pass8"
        elif c["stroke"] == "rgb(0%, 0%, 100%)" and c["width"] == "1.5":
            key = "limit_4yr_pass7"
        elif c["stroke"] == "rgb(0%, 0%, 0%)" and c["dash"] and c["opacity"] in ("1", None):
            key = "median_expected"
        elif c["stroke"] == "rgb(0%, 0%, 0%)" and c["dash"] and c["opacity"] == "0.5":
            key = "thermal_sdb12"
        if key:
            out[key] = {"mass": (10 ** c["log10x"]).tolist(), "sigmav": (10 ** c["log10y"]).tolist()}
    return out
