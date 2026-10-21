"""Read Davis & Lineweaver's (2004) own figures back out of their vector PostScript.

The arXiv source of Davis & Lineweaver, "Expanding Confusion" (PASA 21, 97;
astro-ph/0310808 v2), contains the EPS files that Grace (xmgrace 5.1) wrote for
Figures 1-6. Grace writes every curve as a polyline in page coordinates with
four decimals, every axis tick as a short segment and every label as a string
at a known position. So the curves can be recovered exactly, without
digitising pixels: calibrate each axis from its labelled major ticks, then map
every polyline vertex to physical units.

Figure 1 (the spacetime diagrams, (Omega_M, Omega_Lambda) = (0.3, 0.7), H0 = 70):
    fig1a  time vs proper distance      fig1b  time vs comoving distance
    fig1c  conformal time vs comoving distance
Figure 5: the Perlmutter et al. (1999) supernovae they used for the SR test,
as filled circles with error bars.

The curves are identified by their drawing style, which is the same in all
three panels: our past light cone and the event horizon are thick solid lines,
the particle horizon thick dashed, the Hubble sphere thin solid, worldlines
dotted.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_NUM = re.compile(r"^-?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


@dataclass
class Path2D:
    points: np.ndarray            # (n, 2) page coordinates
    lw: float | None
    dash: str
    color: str | None
    filled: bool = False


@dataclass
class Text:
    x: float
    y: float
    s: str
    size: float


@dataclass
class Circle:
    x: float
    y: float
    r: float


@dataclass
class EPS:
    paths: list = field(default_factory=list)
    texts: list = field(default_factory=list)
    circles: list = field(default_factory=list)


def parse_grace_eps(path: Path) -> EPS:
    """Polylines (with line width / dash / colour), text strings and EARC circles."""
    raw = Path(path).read_text(errors="replace")
    body = raw.split("%%EndSetup", 1)[1]
    # Strings may contain spaces: protect them before splitting.
    strings = []

    def _keep(m):
        strings.append(m.group(1))
        return f" __STR{len(strings) - 1}__ "

    body = re.sub(r"\(((?:[^()\\]|\\.)*)\)", _keep, body)
    toks = body.split()
    out = EPS()
    lw, dash, color = None, "[]", None
    cur: list | None = None
    stack: list = []
    last_point = None
    text_size = 0.0
    i = 0
    while i < len(toks):
        t = toks[i]
        if _NUM.match(t):
            stack.append(float(t))
        elif t == "SLW":
            lw = stack[-1] if stack else lw
            stack.clear()
        elif t == "SD":
            k = i - 1
            while k >= 0 and not toks[k].startswith("["):
                k -= 1
            dash = " ".join(toks[k:i - 1]).replace("[ ", "[").strip()
            stack.clear()
        elif t == "SC":
            color = toks[i - 1]
            stack.clear()
        elif t == "n":
            cur = []
            stack.clear()
        elif t in ("m", "l"):
            if len(stack) >= 2:
                last_point = (stack[-2], stack[-1])
                if cur is not None:
                    cur.append(last_point)
            stack.clear()
        elif t == "EARC":
            if len(stack) >= 6:
                x, y, rx = stack[-6], stack[-5], stack[-4]
                out.circles.append(Circle(x, y, rx))
            cur = None
            stack.clear()
        elif t == "s":
            if cur:
                out.paths.append(Path2D(np.array(cur), lw, dash, color))
            cur = None
            stack.clear()
        elif t == "fill":
            if cur:
                out.paths.append(Path2D(np.array(cur), lw, dash, color, filled=True))
            cur = None
            stack.clear()
        elif t == "CC":
            # "[a 0 0 d 0 0] CC" sets the font size for the following show
            m = re.findall(r"-?\d*\.\d+|-?\d+", " ".join(toks[max(0, i - 7):i]))
            if m:
                text_size = abs(float(m[0])) or abs(float(m[1]))
            stack.clear()
        elif t.startswith("__STR") and i + 1 < len(toks) and toks[i + 1] == "show":
            s = strings[int(t[5:-2])].replace("\\(", "(").replace("\\)", ")")
            if last_point is not None:
                out.texts.append(Text(last_point[0], last_point[1], s, text_size))
            i += 1
            stack.clear()
        elif t in ("c", "GS", "GR"):
            pass
        else:
            stack.clear()
        i += 1
    return out


# --- Axis calibration -------------------------------------------------------------


def _frame(eps: EPS) -> tuple[float, float, float, float]:
    """The plotting frame (xmin, xmax, ymin, ymax) in page units.

    The bottom and left axes are where the tick marks start: many short solid
    segments share that base, while error bars and other short lines do not.
    The right axis is where the right-hand ticks end; the top is the highest
    point of the long curves, which Grace clips to the frame.
    """
    from collections import Counter

    base_y, base_x, end_x = Counter(), Counter(), Counter()
    for p in eps.paths:
        if len(p.points) != 2 or p.filled or p.dash != "[]":
            continue
        (x0, y0), (x1, y1) = p.points
        length = max(abs(x1 - x0), abs(y1 - y0))
        if not 0 < length < 0.05:
            continue
        if abs(x0 - x1) < 1e-6:
            base_y[round(min(y0, y1), 4)] += 1
        elif abs(y0 - y1) < 1e-6:
            base_x[round(min(x0, x1), 4)] += 1
            end_x[round(max(x0, x1), 4)] += 1
    ymin = base_y.most_common(1)[0][0]
    xmin = base_x.most_common(1)[0][0]
    right = [(n, x) for x, n in end_x.items() if x > xmin + 0.5]
    long_pts = np.vstack([p.points for p in eps.paths if len(p.points) > 50 and not p.filled])
    xmax = max(right)[1] if right else long_pts[:, 0].max()
    return xmin, xmax, ymin, float(long_pts[:, 1].max())


def _major_ticks(eps: EPS, axis: str, where: float, tol: float = 2e-4) -> np.ndarray:
    """Positions of the longest ticks on an axis line ('x' ticks sit on y = where)."""
    ticks = []
    for p in eps.paths:
        if len(p.points) != 2 or p.filled:
            continue
        (x0, y0), (x1, y1) = p.points
        if max(abs(x1 - x0), abs(y1 - y0)) > 0.05:
            continue  # a line across the panel (Grace's "now" line), not a tick
        if axis == "x" and abs(x0 - x1) < 1e-6 and abs(min(y0, y1) - where) < tol:
            ticks.append((x0, abs(y1 - y0)))
        if axis == "y" and abs(y0 - y1) < 1e-6 and abs(min(x0, x1) - where) < tol:
            ticks.append((y0, abs(x1 - x0)))
    if not ticks:
        raise ValueError(f"no {axis} ticks at {where}")
    ticks = np.array(ticks)
    longest = ticks[:, 1].max()
    return np.sort(ticks[np.abs(ticks[:, 1] - longest) < 1e-4, 0])


def _numeric_labels(eps: EPS, axis: str, near: float, band: float,
                    span: tuple[float, float]) -> list[tuple[float, float]]:
    """(page position along the axis, value) of numeric labels in a band beside the
    axis and within ``span`` (page coordinates along it, padded a little)."""
    out = []
    lo, hi = span[0] - 0.02, span[1] + 0.02
    for t in eps.texts:
        try:
            v = float(t.s.replace("−", "-"))
        except ValueError:
            continue
        if axis == "x" and near - band < t.y < near and lo < t.x < hi:
            out.append((t.x, v))
        if axis == "y" and near - band < t.x < near and lo < t.y < hi:
            out.append((t.y, v))
    return sorted(out)


@dataclass
class Calibration:
    x0: float
    xs: float  # value = (page - x0) / xs
    y0: float
    ys: float
    x_ticks: list
    y_ticks: list
    rms_x: float
    rms_y: float

    def to_world(self, pts: np.ndarray) -> np.ndarray:
        return np.column_stack([(pts[:, 0] - self.x0) / self.xs, (pts[:, 1] - self.y0) / self.ys])


def calibrate(eps: EPS) -> Calibration:
    """Linear page -> data maps from the bottom and left axes' labelled major ticks."""
    xmin, xmax, ymin, ymax = _frame(eps)
    xt = _major_ticks(eps, "x", ymin)
    yt = _major_ticks(eps, "y", xmin)
    xl = _numeric_labels(eps, "x", ymin, 0.06, (xt.min(), xt.max()))
    yl = _numeric_labels(eps, "y", xmin, 0.08, (yt.min(), yt.max()))
    if len(xl) != len(xt) or len(yl) != len(yt):
        raise ValueError(f"tick/label mismatch: x {len(xt)}/{len(xl)}, y {len(yt)}/{len(yl)}")
    xv = np.array([v for _, v in xl])
    yv = np.array([v for _, v in yl])
    bx = np.polyfit(xv, xt, 1)
    by = np.polyfit(yv, yt, 1)
    rx = xt - np.polyval(bx, xv)
    ry = yt - np.polyval(by, yv)
    return Calibration(x0=bx[1], xs=bx[0], y0=by[1], ys=by[0],
                       x_ticks=list(zip(xv.tolist(), xt.tolist())),
                       y_ticks=list(zip(yv.tolist(), yt.tolist())),
                       rms_x=float(np.sqrt(np.mean(rx**2))), rms_y=float(np.sqrt(np.mean(ry**2))))


def right_axis_ticks(eps: EPS) -> list[tuple[float, float]]:
    """(page y, value) of the labelled major ticks on the right axis (the scale factor)."""
    xmin, xmax, ymin, ymax = _frame(eps)
    yt = _major_ticks_right(eps, xmax)
    labels = []
    for t in eps.texts:
        try:
            v = float(t.s)
        except ValueError:
            continue
        if xmax < t.x < xmax + 0.08:
            labels.append((t.y, v))
    labels.sort()
    if len(labels) != len(yt):
        raise ValueError(f"right axis: {len(yt)} ticks, {len(labels)} labels")
    return [(float(y), v) for y, (_, v) in zip(yt, labels)]


def _major_ticks_right(eps: EPS, xmax: float) -> np.ndarray:
    ticks = []
    for p in eps.paths:
        if len(p.points) != 2 or p.filled:
            continue
        (x0, y0), (x1, y1) = p.points
        if abs(y0 - y1) < 1e-6 and abs(max(x0, x1) - xmax) < 2e-4 and 0 < abs(x1 - x0) < 0.05:
            ticks.append((y0, abs(x1 - x0)))
    ticks = np.array(ticks)
    longest = ticks[:, 1].max()
    return np.sort(ticks[np.abs(ticks[:, 1] - longest) < 1e-4, 0])


# --- Figure 1: which polyline is which ----------------------------------------------


def _style(p: Path2D) -> str:
    thick = p.lw is not None and p.lw > 0.002
    if p.dash.startswith("[0.0015"):
        return "dotted"
    if p.dash != "[]":
        return "thick_dashed" if thick else "thin_dashed"
    return "thick" if thick else "thin"


def figure1_curves(path: Path) -> dict:
    """The named curves of one Figure 1 panel, in data units, split into the halves
    left (x < 0) and right (x > 0) of our worldline.

    Returns {"calibration", "light_cone", "event_horizon", "particle_horizon",
    "hubble_sphere", "worldlines": [array, ...], "frame": (xmin, xmax, ymin, ymax)}
    with each curve a list of (n, 2) arrays [x, y] (one per half).
    """
    eps = parse_grace_eps(path)
    cal = calibrate(eps)
    xmin, xmax, ymin, ymax = _frame(eps)
    centre = cal.x0  # page x of distance 0
    groups: dict[str, list] = {"thick": [], "thick_dashed": [], "thin": [], "dotted": []}
    for p in eps.paths:
        if p.filled or len(p.points) < 20:
            continue
        st = _style(p)
        if st in groups:
            groups[st].append(p)

    def world(p):
        return cal.to_world(p.points)

    out = {"calibration": cal, "frame": tuple(cal.to_world(np.array([[xmin, ymin], [xmax, ymax]])).ravel()),
           "worldlines": [], "light_cone": [], "event_horizon": [], "particle_horizon": [],
           "hubble_sphere": []}
    # Thick solid: the light cone (closes on our worldline below the top: apex at
    # x = 0) and the event horizon (runs to the top of the frame).
    for p in groups["thick"]:
        w = world(p)
        top_reached = p.points[:, 1].max() > ymax - 1e-3
        out["event_horizon" if top_reached else "light_cone"].append(w)
    out["particle_horizon"] = [world(p) for p in groups["thick_dashed"]]
    # Thin solid long curves: the Hubble sphere, plus the edges of the shaded
    # regions (which duplicate the event horizon); keep the ones that are not
    # within a hair of an event-horizon curve.
    eh_pts = np.vstack([p.points for p in groups["thick"] if p.points[:, 1].max() > ymax - 1e-3]) \
        if out["event_horizon"] else np.zeros((0, 2))
    for p in groups["thin"]:
        if len(p.points) < 50:
            continue
        if np.ptp(cal.to_world(p.points)[:, 0]) > 30.0:
            continue  # the edge of a shaded region, not the Hubble sphere (|D_HS| < 17 Gly)
        if eh_pts.size:
            d = np.min(np.hypot(p.points[::25, None, 0] - eh_pts[None, ::5, 0],
                                p.points[::25, None, 1] - eh_pts[None, ::5, 1]), axis=1)
            if np.median(d) < 5e-4:
                continue
        xs = p.points[:, 0]
        if np.all(np.abs(xs - centre) < 2e-3):
            continue
        out["hubble_sphere"].append(world(p))
    for p in groups["dotted"]:
        w = world(p)
        if np.all(np.abs(p.points[:, 0] - centre) < 2e-3):
            continue  # our own worldline
        out["worldlines"].append(w)
    out["right_axis"] = right_axis_ticks(eps)
    out["right_axis_y_to_value"] = [((y - cal.y0) / cal.ys, v) for y, v in out["right_axis"]]
    return out


# --- Figure 5: the supernova points ----------------------------------------------------


def figure5_points(path: Path) -> dict:
    """Perlmutter et al. (1999) supernovae as plotted in D&L's Fig. 5: z, m_B, and the
    upper and lower error bars recovered from the vertical segments at each circle."""
    eps = parse_grace_eps(path)
    cal = calibrate(eps)
    circles = [c for c in eps.circles if c.r < 0.01]
    seen = set()
    pts = []
    for c in circles:  # each circle is drawn twice (fill + outline)
        key = (round(c.x, 4), round(c.y, 4))
        if key in seen:
            continue
        seen.add(key)
        up = dn = None
        for p in eps.paths:
            if len(p.points) != 2 or p.filled:
                continue
            (x0, y0), (x1, y1) = p.points
            if abs(x0 - c.x) < 2e-4 and abs(x1 - c.x) < 2e-4 and abs(y0 - c.y) < 2e-4:
                if y1 > y0:
                    up = y1 - y0
                else:
                    dn = y0 - y1
        pts.append((c.x, c.y, up if up is not None else 0.0, dn if dn is not None else 0.0))
    pts = np.array(sorted(pts))
    zm = cal.to_world(pts[:, :2])
    return {"z": zm[:, 0], "m": zm[:, 1], "err_up": pts[:, 2] / cal.ys, "err_dn": pts[:, 3] / cal.ys,
            "calibration": cal,
            "curves": [cal.to_world(p.points) for p in eps.paths if len(p.points) > 50 and not p.filled],
            "curve_styles": [(_style(p)) for p in eps.paths if len(p.points) > 50 and not p.filled]}
