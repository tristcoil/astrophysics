"""Read the published plots of the CMS four-lepton outreach analysis back into numbers.

Record 12360 ships its final plots as PDFs written by ROOT. They are vector
graphics: every histogram is a polyline and every data point a filled circle,
with coordinates printed to six significant digits. Mapping the page
coordinates back through the axis ticks recovers, bin by bin, the observed
event counts, the stacked Z gamma* / ZZ background and the Higgs signal the
analysis published -- an exact answer key for our replica of its code.
"""

from __future__ import annotations

import json
import re
import zlib
from pathlib import Path

import numpy as np

import config as cfg

#: Histogram ranges of the published plots (histograms.py: 36 bins each).
RANGES = {"Higgs_mass": (70.0, 180.0), "Z1_mass": (40.0, 160.0), "Z2_mass": (12.0, 160.0)}
NBINS = 36
STATES = {"FourMuons": "4mu", "FourElectrons": "4e", "TwoMuonsTwoElectrons": "2e2mu",
          "combined": "combined"}

_NUM = r"(-?\d+(?:\.\d+)?)"


def _content_stream(path: Path) -> str:
    raw = path.read_bytes()
    streams = []
    for m in re.finditer(rb"stream\r?\n", raw):
        end = raw.find(b"endstream", m.end())
        try:
            streams.append(zlib.decompress(raw[m.end():end]))
        except zlib.error:
            continue
    return max(streams, key=len).decode("latin-1")


def _frame(text: str) -> tuple[float, float, float, float]:
    """The plot frame, drawn as 'x y w h re S' right after the white pads."""
    m = re.search(_NUM + " " + _NUM + " " + _NUM + " " + _NUM + r" re S", text)
    return tuple(float(v) for v in m.groups())


def _y_axis(text: str, x0: float, y0: float) -> tuple[float, float]:
    """Linear map value -> page y from the labelled major ticks of the y axis."""
    ticks = []
    for m in re.finditer(_NUM + " " + _NUM + " m " + _NUM + " " + _NUM + " l S", text):
        xa, ya, xb, yb = (float(v) for v in m.groups())
        if abs(xb - x0) < 1e-3 and abs(ya - yb) < 1e-6 and 12.0 < xa - x0 < 15.0:
            ticks.append(ya)
    ticks = sorted(set(ticks))
    labels = []
    for m in re.finditer(r"1 0 0 1 " + _NUM + " " + _NUM + r" cm BT /F\d+ [\d.]+ Tf \(([^)]*)\) Tj",
                         text):
        x, y, s = float(m.group(1)), float(m.group(2)), m.group(3)
        # left of the frame and not below it (the x axis's first label is both)
        if x < x0 and y > y0 - 9.0 and re.fullmatch(r"-?\d+(\.\d+)?", s):
            labels.append((y, float(s)))
    # each label sits ~6 pt below the tick it names: pair them by order
    labels.sort()
    if len(labels) != len(ticks):
        raise ValueError(f"{len(labels)} y labels for {len(ticks)} major ticks")
    vals = np.array([v for _, v in labels])
    a, b = np.polyfit(vals, np.array(ticks), 1)
    resid = np.array(ticks) - (a * vals + b)
    if np.max(np.abs(resid)) > 0.01:
        raise ValueError("y-axis ticks are not linear in their labels")
    return b, a                       # y_page = b + a * value


def _step_values(points: list[tuple[float, float]], x0: float, width: float) -> np.ndarray:
    """Bin heights of a ROOT step polyline (page coordinates)."""
    segs = [(min(p[0], q[0]), max(p[0], q[0]), p[1])
            for p, q in zip(points[:-1], points[1:]) if abs(p[1] - q[1]) < 1e-9]
    out = np.full(NBINS, np.nan)
    for k in range(NBINS):
        xc = x0 + (k + 0.5) * width
        for xa, xb, y in segs:
            if xa - 1e-6 <= xc <= xb + 1e-6:
                out[k] = y
                break
    return out


def _path_after(text: str, marker: str, stop: str) -> list[tuple[float, float]]:
    i = text.find(marker)
    if i < 0:
        raise ValueError(f"no path after {marker!r}")
    j = text.find(" m ", i)
    start = text.rfind(" ", i, j - 1)
    start = text.rfind(" ", i, start - 1)
    k = text.find(stop, j)
    nums = [float(v) for v in re.findall(_NUM, text[start:k])]
    return list(zip(nums[0::2], nums[1::2]))


def read_plot(path: Path, variable: str) -> dict:
    text = _content_stream(path)
    x0, y0, w, h = _frame(text)
    lo, hi = RANGES[variable]
    width = w / NBINS
    yb, ya = _y_axis(text, x0, y0)

    def to_events(y):
        return (np.asarray(y) - yb) / ya

    # stacked signal + background: the red line ("1 0 0 RG"), stroked with S
    red = _path_after(text, "1 0 0 RG 1 0 0 rg", " S")
    total = to_events(_step_values(red, x0, width))
    # background: the light-blue fill (100, 192, 232)/255
    blue = _path_after(text, "0.392157 0.752941 0.909804 rg", " f*")
    background = to_events(_step_values(blue, x0, width))
    # data: filled circles drawn as two Bezier halves
    circ = re.compile(_NUM + " " + _NUM + " m " + r"(\S+) (\S+) " + _NUM + r" \S+ \5 \S+ c "
                      + r"\S+ \S+ \S+ \S+ \S+ \S+ c h f")
    data = np.zeros(NBINS)
    for m in circ.finditer(text):
        xa, yc, xb = float(m.group(1)), float(m.group(2)), float(m.group(5))
        xc = 0.5 * (xa + xb)
        k = int(np.floor((xc - x0) / width))
        # data points sit exactly on bin centres; the legend's sample marker does not
        if 0 <= k < NBINS and abs(xc - (x0 + (k + 0.5) * width)) < 0.01:
            data[k] = to_events(yc)
    if np.max(np.abs(data - np.round(data))) > 2e-3:
        raise ValueError(f"{path.name}: data points are not integers")
    edges = np.linspace(lo, hi, NBINS + 1)
    return {"edges": edges.tolist(), "data": np.round(data).astype(int).tolist(),
            "background": background.tolist(), "total": total.tolist(),
            "signal": (total - background).tolist()}


def read_all(plots_dir: Path = cfg.OUTREACH_PLOTS_DIR) -> dict:
    """{state: {variable: histogram dict}} for all twelve published plots."""
    out = {}
    for fname, state in STATES.items():
        for var in RANGES:
            p = plots_dir / f"{fname}_{var}.pdf"
            out.setdefault(state, {})[var] = read_plot(p, var)
    return out


def load(refresh: bool = False) -> dict:
    cache = cfg.CACHE_DIR / "outreach_published_histograms.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())
    out = read_all()
    cache.write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    h = load(refresh=True)["combined"]["Higgs_mass"]
    for lo, n, b, s in zip(h["edges"], h["data"], h["background"], h["signal"]):
        print(f"{lo:7.2f}  data {n:3d}   ZZ {b:7.4f}   H {s:7.4f}")
    print("total data", sum(h["data"]), " background", round(sum(h["background"]), 3),
          " signal", round(sum(h["signal"]), 3))
