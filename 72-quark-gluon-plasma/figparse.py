"""Read published vector figures back into numbers.

Andronic et al. 2018 made their figures with ROOT and put them in the arXiv
source as PDF. Converted to SVG (poppler's ``pdftocairo``), every marker,
error bar, model line and axis tick is a path with exact coordinates, so the
plotted values can be recovered to a small fraction of a marker's size:

* Fig. 1: the measured yields they fitted and their model's yields at the
  best fit (T = 156.5 MeV), particle by particle;
* Fig. 4: the chemical freeze-out points (T, mu_B) at all energies, theirs
  and other groups', drawn on the QCD phase diagram.

Axes are calibrated on their own tick marks (major ticks at known values), so
nothing is read by eye.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

import config as cfg
import fetch_data

_TOKEN = re.compile(r"[MLCZ]|-?[0-9.]+(?:e-?[0-9]+)?")


@dataclass
class Path_:
    stroke: str
    fill: str
    width: float
    subpaths: list            # list of (k, 2) arrays in PDF points, y upward

    def bbox_centres(self) -> np.ndarray:
        return np.array([[0.5 * (s[:, 0].min() + s[:, 0].max()),
                          0.5 * (s[:, 1].min() + s[:, 1].max())] for s in self.subpaths])


def _parse_d(d: str) -> list:
    toks = _TOKEN.findall(d)
    subs, cur, i, cmd = [], [], 0, None
    while i < len(toks):
        t = toks[i]
        if t in "MLCZ":
            cmd = t
            i += 1
            if t == "Z":
                if cur:
                    subs.append(cur)
                cur = []
                continue
            if t == "M" and cur:
                subs.append(cur)
                cur = []
            continue
        n = 6 if cmd == "C" else 2
        vals = [float(x) for x in toks[i:i + n]]
        cur.append(vals[-2:])           # end point (control points ignored)
        if cmd == "C":
            cur.append(vals[0:2])
            cur.append(vals[2:4])
        i += n
    if cur:
        subs.append(cur)
    return [np.array(s) for s in subs if len(s)]


@lru_cache(maxsize=None)
def read_svg_paths(pdf: Path) -> tuple:
    """All drawn paths of a one-page PDF, in PDF coordinates (y up)."""
    svg = pdf.with_suffix(".svg")
    cache = cfg.DERIVED_DIR / "figparse" / (pdf.parent.parent.name + "_" + svg.name)
    cache.parent.mkdir(parents=True, exist_ok=True)
    if not cache.exists():
        subprocess.run(["pdftocairo", "-svg", str(pdf), str(cache)], check=True)
    s = cache.read_text()
    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', s)
    H = float(m.group(2))
    body = s[s.index("</defs>"):]
    out = []
    for p in re.findall(r"<path ([^>]*)/>", body):
        a = dict(re.findall(r'([\w-]+)="([^"]*)"', p))
        subs = _parse_d(a.get("d", ""))
        tr = a.get("transform", "")
        if tr.startswith("matrix(1, 0, 0, -1, 0,"):
            pass                                   # already PDF coordinates
        elif tr == "":
            subs = [np.column_stack([s[:, 0], H - s[:, 1]]) for s in subs]
        else:
            mm = [float(x) for x in re.findall(r"-?[0-9.]+", tr)]
            a_, b_, c_, d_, e_, f_ = mm
            # SVG point = M (x, y); convert to PDF y-up afterwards
            subs = [np.column_stack([a_ * s[:, 0] + c_ * s[:, 1] + e_,
                                     H - (b_ * s[:, 0] + d_ * s[:, 1] + f_)]) for s in subs]
        out.append(Path_(a.get("stroke", "none"), a.get("fill", "none"),
                         float(a.get("stroke-width", "0") or 0), subs))
    return tuple(out)


def segments(paths, stroke: str, width: float | None = None) -> np.ndarray:
    """Straight two-point segments (x0, y0, x1, y1) of one stroke colour."""
    segs = []
    for p in paths:
        if p.stroke != stroke or (width is not None and abs(p.width - width) > 1e-6):
            continue
        for s in p.subpaths:
            for a, b in zip(s[:-1], s[1:]):
                segs.append([a[0], a[1], b[0], b[1]])
    return np.array(segs)


def filled_markers(paths, fill: str) -> np.ndarray:
    pts = [p.bbox_centres() for p in paths if p.fill == fill and p.stroke == "none"]
    return np.vstack(pts) if pts else np.zeros((0, 2))


def frames(paths, stroke: str = "rgb(0%, 0%, 0%)") -> list:
    """Axis-aligned rectangles (x0, x1, y0, y1) drawn in one colour, tallest first."""
    out = set()
    for p in paths:
        if p.stroke != stroke:
            continue
        for s in p.subpaths:
            if len(s) == 4:
                xs, ys = np.unique(np.round(s[:, 0], 3)), np.unique(np.round(s[:, 1], 3))
                if len(xs) == 2 and len(ys) == 2 and np.ptp(xs) > 50 and np.ptp(ys) > 50:
                    out.add((xs[0], xs[1], ys[0], ys[1]))
    return sorted(out, key=lambda r: -(r[3] - r[2]))


def _axis_ticks(segs: np.ndarray, x_axis: float, y_range: tuple, tol: float = 0.05) -> tuple:
    """Horizontal tick marks starting on a vertical axis at x_axis; returns
    (y positions, lengths)."""
    ys, ls = [], []
    for x0, y0, x1, y1 in segs:
        if abs(y0 - y1) < tol and (abs(x0 - x_axis) < tol or abs(x1 - x_axis) < tol) \
                and y_range[0] - 1 <= y0 <= y_range[1] + 1:
            length = abs(x1 - x0)
            if 0.5 < length < 30:
                ys.append(y0)
                ls.append(length)
    order = np.argsort(ys)
    return np.array(ys)[order], np.array(ls)[order]


# ---------------------------------------------------------------------------
# Andronic et al. 2018, Fig. 1: yields and model
# ---------------------------------------------------------------------------

FIG1_ORDER = ("pi+", "pi-", "K+", "K-", "K0S", "phi", "p", "pbar", "Lambda", "Lambdabar",
              "Xi-", "Xibar+", "Omega-", "Omegabar+", "d", "dbar", "He3", "He3bar",
              "H3L", "H3Lbar", "He4", "He4bar")
RED, BLUE, BLACK = "rgb(100%, 0%, 0%)", "rgb(0%, 0%, 100%)", "rgb(0%, 0%, 0%)"


def andronic_fig1() -> dict:
    """Data and model yields read from Andronic et al. 2018, Fig. 1 (upper panel).

    The y axis is logarithmic with major ticks at 10^-6 ... 10^3 (labels read
    with pdftotext); data are red filled circles with red error bars, the
    model is a blue horizontal bar per particle.
    """
    pdf = fetch_data.arxiv_dir("1710.09425") / "yields_fit2_276.pdf"
    paths = read_svg_paths(pdf)
    black = segments(paths, BLACK, 0.75)
    x_lo, x_hi, y_lo, y_hi = frames(paths)[0]       # upper panel: the tallest frame
    ty, tl = _axis_ticks(black, x_lo, (y_lo, y_hi))
    major = ty[tl > 0.75 * tl.max()]
    major = np.unique(np.round(major, 3))
    decades = _lowest_decade(pdf) + np.arange(len(major))
    slope, icpt = np.polyfit(major, decades, 1)
    resid = np.abs(np.polyval([slope, icpt], major) - decades).max()
    to_val = lambda y: 10.0 ** (slope * y + icpt)   # noqa: E731

    # Data markers are the large red circles (the tiny red dots sit on the
    # model values); the legend's marker is off the particles' regular x grid.
    big, small = [], []
    for p in paths:
        if p.fill != RED or p.stroke != "none":
            continue
        for s in p.subpaths:
            c = [0.5 * (s[:, 0].min() + s[:, 0].max()), 0.5 * (s[:, 1].min() + s[:, 1].max())]
            if not (y_lo < c[1] < y_hi):
                continue
            (big if np.ptp(s[:, 1]) > 5 else small).append(c)
    big = np.unique(np.round(np.array(big), 2), axis=0)
    xu = np.unique(np.round(big[:, 0], 1))
    step = np.median(np.diff(xu))
    on_grid = np.abs((big[:, 0] - xu.min()) / step - np.round((big[:, 0] - xu.min()) / step)) < 0.05
    xs = big[on_grid][np.argsort(big[on_grid][:, 0])]
    red = segments(paths, RED, 1.5)
    blue = segments(paths, BLUE, 1.5)
    out = {"axis_fit_residual_decades": float(resid), "entries": {}}
    if len(xs) != len(FIG1_ORDER):
        raise RuntimeError(f"found {len(xs)} data markers, expected {len(FIG1_ORDER)}")
    for key, (x, y) in zip(FIG1_ORDER, xs):
        # error bar: vertical red segments at this x
        v = red[(np.abs(red[:, 0] - x) < 0.3) & (np.abs(red[:, 2] - x) < 0.3)
                & (np.abs(red[:, 0] - red[:, 2]) < 0.05)]
        lo = min([y] + list(v[:, 1]) + list(v[:, 3])) if len(v) else y
        hi = max([y] + list(v[:, 1]) + list(v[:, 3])) if len(v) else y
        # model: horizontal blue segments whose x-range covers this x
        h = blue[(np.abs(blue[:, 1] - blue[:, 3]) < 0.05)
                 & (np.minimum(blue[:, 0], blue[:, 2]) <= x + 0.5)
                 & (np.maximum(blue[:, 0], blue[:, 2]) >= x - 0.5)]
        model_y = float(np.median(h[:, 1])) if len(h) else float("nan")
        dots = [c[1] for c in small if abs(c[0] - x) < 0.3]
        out["entries"][key] = {"data": to_val(y), "data_lo": to_val(lo), "data_hi": to_val(hi),
                               "model": to_val(model_y), "x": float(x),
                               "model_dot": to_val(dots[0]) if dots else float("nan")}
    return out


def andronic_fig4() -> dict:
    """Chemical freeze-out points from Andronic et al. 2018, Fig. 4 (T vs mu_B).

    Axes: mu_B logarithmic (major ticks at 1, 10, 100, 1000 MeV), T linear
    (major ticks every 20 MeV from 0). Markers are closed outlines in one
    colour per analysis; their error bars are two-point segments in the same
    colour. Returns one record per marker with its group (colour + number of
    outline vertices), centre, and error-bar extent.
    """
    pdf = fetch_data.arxiv_dir("1710.09425") / "t-mu_2017_lqcd3.pdf"
    paths = read_svg_paths(pdf)
    x0, x1, y0, y1 = frames(paths)[0]
    black = segments(paths, BLACK)
    # x ticks: vertical segments rising from the bottom axis
    xt, xl = [], []
    for a, b, c, d in black:
        if abs(a - c) < 0.01 and (abs(b - y0) < 0.01 or abs(d - y0) < 0.01) and x0 - 1 <= a <= x1 + 1:
            length = abs(d - b)
            if 0.5 < length < 30:
                xt.append(a)
                xl.append(length)
    xt, xl = np.array(xt), np.array(xl)
    xmaj = np.unique(np.round(xt[xl > 0.75 * xl.max()], 3))
    xdec = np.arange(len(xmaj)) + 0.0          # 1, 10, 100, 1000 MeV
    sx, ix = np.polyfit(xmaj, xdec, 1)
    ty, tl = _axis_ticks(black, x0, (y0, y1))
    ymaj = np.unique(np.round(ty[tl > 0.75 * tl.max()], 3))
    yval = 20.0 * np.arange(len(ymaj))
    sy, iy = np.polyfit(ymaj, yval, 1)
    to_mu = lambda x: 10.0 ** (sx * x + ix)     # noqa: E731
    to_T = lambda y: sy * y + iy                 # noqa: E731
    out = {"x_axis_residual": float(np.abs(np.polyval([sx, ix], xmaj) - xdec).max()),
           "y_axis_residual_MeV": float(np.abs(np.polyval([sy, iy], ymaj) - yval).max()),
           "points": []}
    bars: dict = {}
    for p in paths:
        if p.fill == "none" and 5 <= p.width <= 10:
            bars.setdefault(p.stroke, []).extend(s for s in p.subpaths if len(s) == 2)
    for p in paths:
        if p.fill != "none" or p.stroke in ("none",) or p.width < 5 or p.width > 10:
            continue
        segs = bars.get(p.stroke, [])
        for s in p.subpaths:
            if len(s) < 3:
                continue
            cx = 0.5 * (s[:, 0].min() + s[:, 0].max())
            cy = 0.5 * (s[:, 1].min() + s[:, 1].max())
            if not (x0 < cx < x1 and y0 < cy < y1):
                continue
            # error bars of this marker: vertical/horizontal segments through its centre
            hx, hy = 0.5 * np.ptp(s[:, 0]) + 0.3, 0.5 * np.ptp(s[:, 1]) + 0.3
            up = [max(q[:, 1]) for q in segs if np.ptp(q[:, 0]) < 0.01
                  and abs(q[0, 0] - cx) < 0.05 and cy <= min(q[:, 1]) <= cy + hy]
            dn = [min(q[:, 1]) for q in segs if np.ptp(q[:, 0]) < 0.01
                  and abs(q[0, 0] - cx) < 0.05 and cy - hy <= max(q[:, 1]) <= cy]
            lf = [min(q[:, 0]) for q in segs if np.ptp(q[:, 1]) < 0.01
                  and abs(q[0, 1] - cy) < 0.05 and cx - hx <= max(q[:, 0]) <= cx]
            rt = [max(q[:, 0]) for q in segs if np.ptp(q[:, 1]) < 0.01
                  and abs(q[0, 1] - cy) < 0.05 and cx <= min(q[:, 0]) <= cx + hx]
            out["points"].append({
                "group": f"{p.stroke}|{len(s)}", "muB": float(to_mu(cx)), "T": float(to_T(cy)),
                "T_hi": float(to_T(max(up))) if up else float(to_T(cy)),
                "T_lo": float(to_T(min(dn))) if dn else float(to_T(cy)),
                "muB_lo": float(to_mu(min(lf))) if lf else float(to_mu(cx)),
                "muB_hi": float(to_mu(max(rt))) if rt else float(to_mu(cx)),
                "on_left_axis": bool(lf and min(lf) <= x0 + 0.05)})
    return out


def _lowest_decade(pdf: Path) -> int:
    """Exponent of the lowest labelled decade of a log axis: ROOT writes the
    negative ones inline ("10−6"), so the most negative inline exponent."""
    txt = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True,
                         text=True, check=True).stdout
    exps = [-int(e) for e in re.findall(r"10[−-](\d+)", txt)]
    if not exps:
        raise RuntimeError("no negative decade labels found")
    return min(exps)
