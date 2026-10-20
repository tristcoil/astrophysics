"""Keeper animation: four years of LAGEOS and LARES 2, one frame per week.

    python animation.py      # writes figures/orbital_planes.mp4 (~1-2 min)

Left: the two real orbital planes (weekly mean node and inclination from the
ILRS orbits, to scale, Earth radius = 1) seen from a fixed direction in space.
The Earth's equatorial bulge swings them round in opposite directions, ~125 deg
a year each, LAGEOS eastward and LARES 2 westward. Right: what is left of the
sum of their nodes once the bulge (and its J2 errors) and the fitted tide-model
errors are taken out: frame dragging, building up week by week against general
relativity's 61.4 mas/yr. At the satellites, 1 mas of node is 6 cm.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.patches import Circle

import config as C
import frames as F
import framedrag as FD
import nodes
import plotstyle as ps

ps.use_style()
RE_KM = 6378.137


def _basis(view_lon_deg=-35.0, view_lat_deg=22.0):
    lon, lat = np.radians(view_lon_deg), np.radians(view_lat_deg)
    d = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])   # towards viewer
    e1 = np.cross([0.0, 0.0, 1.0], d)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    return d, e1, e2


def _orbit_points(a, inc, node, n=400):
    u = np.linspace(0, 2 * np.pi, n)
    en = np.array([np.cos(node), np.sin(node), 0.0])
    h = np.array([np.sin(inc) * np.sin(node), -np.sin(inc) * np.cos(node), np.cos(inc)])
    em = np.cross(h, en)
    return a * (np.outer(np.cos(u), en) + np.outer(np.sin(u), em)), a * en


def _draw_curve(ax, P, d, e1, e2, color, lw, zorder):
    """Draw a 3-D curve in projection; the part hidden behind the Earth faint and dotted."""
    x, y, depth = P @ e1, P @ e2, P @ d
    hidden = (depth < 0) & (np.hypot(x, y) < 1.0)
    for mask, style in ((~hidden, dict(lw=lw, alpha=1.0, ls="-")),
                        (hidden, dict(lw=lw * 0.7, alpha=0.35, ls=":"))):
        xs, ys = np.where(mask, x, np.nan), np.where(mask, y, np.nan)
        ax.plot(xs, ys, color=color, zorder=zorder if style["ls"] == "-" else 1, **style)


def make(out=None, fps: int = 14, dpi: int = 110) -> str:
    out = out or C.FIGURE_DIR / "orbital_planes.mp4"
    we = {s: nodes.weekly_elements(s, C.PRIMARY_AC, C.BASE_TAG) for s in ("lageos1", "lares2")}
    curves = FD.curves(C.PRIMARY_AC, C.BASE_TAG, "baseline")
    yr_c = F.mjd_to_year(curves["mjd"])
    # frames = the weeks both satellites have, matched by date (LAGEOS starts three weeks earlier)
    m1 = np.round(we["lageos1"]["mjd"], 1)
    m2 = np.round(we["lares2"]["mjd"], 1)
    common = np.intersect1d(m1, m2)
    idx = {"lageos1": np.searchsorted(m1, common), "lares2": np.searchsorted(m2, common)}
    n_frames = common.size
    d, e1, e2 = _basis()

    fig = plt.figure(figsize=(12.8, 6.4))
    ax = fig.add_axes([0.0, 0.02, 0.5, 0.9])
    ax2 = fig.add_axes([0.58, 0.16, 0.38, 0.68])
    lt = curves["lt_sum"]

    def frame(k):
        ax.clear()
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_xlim(-2.45, 2.45)
        ax.set_ylim(-2.3, 2.3)
        ax.add_patch(Circle((0, 0), 1.0, color=ps.BLUE_150, zorder=2))
        # equator and spin axis
        u = np.linspace(0, 2 * np.pi, 300)
        eq = np.stack([np.cos(u), np.sin(u), np.zeros_like(u)], 1)
        x, y, dep = eq @ e1, eq @ e2, eq @ d
        ax.plot(np.where(dep >= 0, x, np.nan), np.where(dep >= 0, y, np.nan), color=ps.BLUE_550,
                lw=0.8, zorder=3)
        pole = np.array([0, 0, 1.0])
        # the axis only outside the globe (the north pole faces us)
        ax.plot([1.0 * (pole @ e1), 1.4 * (pole @ e1)], [1.0 * (pole @ e2), 1.4 * (pole @ e2)],
                color=ps.INK, lw=1.2, zorder=5)
        ax.plot(pole @ e1, pole @ e2, "o", color=ps.INK, ms=2.5, zorder=5)
        ax.text(1.45 * (pole @ e1), 1.45 * (pole @ e2), "N", fontsize=9, ha="center", va="bottom")
        for sat in ("lageos1", "lares2"):
            w = we[sat]
            st = ps.SAT_STYLE[sat]
            kk = idx[sat][k]
            a = w["a"][kk] / 1e3 / RE_KM
            P, nodevec = _orbit_points(a, w["inc_cip"][kk], w["node_cip"][kk])
            _draw_curve(ax, P, d, e1, e2, st["color"], 2.0, 6)
            nx, ny, nd = nodevec @ e1, nodevec @ e2, nodevec @ d
            behind = nd < 0 and np.hypot(nx, ny) < 1.0
            ax.plot(nx, ny, st["marker"], color=st["color"], ms=7, mec=ps.SURFACE,
                    zorder=1 if behind else 7, alpha=0.35 if behind else 1.0)
        date = F.mjd_to_year(common[k])
        ax.text(-2.3, 2.12, f"{date:.2f}", fontsize=14, fontweight="bold")
        ax.text(-2.3, -1.78, "LAGEOS  (node swings east, +125°/yr)", fontsize=9, color=ps.BLUE)
        ax.text(-2.3, -1.98, "LARES 2 (node swings west, −125°/yr)", fontsize=9, color=ps.ORANGE)
        ax.text(-2.3, -2.2, "orbits to scale; markers: where each orbit crosses the equator northward",
                fontsize=7.5, color=ps.MUTED)
        # right panel: frame dragging accumulating
        ax2.clear()
        tt = curves["t_yr"]
        j = int(np.clip(np.searchsorted(curves["mjd"], common[k] + 3.5) - 1, 0, len(tt) - 1))
        ax2.plot(yr_c, lt * tt, color=ps.FRAME, lw=1.6, label=f"general relativity, {lt:.1f} mas/yr")
        ax2.plot(yr_c[:j + 1], curves["sum"]["clean"][:j + 1], "D", color=ps.INK, ms=2.6,
                 label="LAGEOS + LARES 2, measured")
        ax2.axhline(0, color=ps.MUTED, lw=0.8, ls="--")
        ax2.set_xlim(yr_c[0] - 0.05, yr_c[-1] + 0.05)
        ax2.set_ylim(-15, 275)
        ax2.set_ylabel("node residual of the sum (mas)")
        ax2.set_title("Frame dragging, building up", fontsize=11)
        ax2.legend(loc="upper left", fontsize=8)
        val = curves["sum"]["clean"][j]
        ax2.text(0.98, 0.04, f"{val:.0f} mas = {val * 4.848e-9 * 12_268e3:.1f} m at the satellites",
                 transform=ax2.transAxes, ha="right", fontsize=9, color=ps.INK)
        return []

    anim = FuncAnimation(fig, frame, frames=n_frames, blit=False)
    anim.save(str(out), writer=FFMpegWriter(fps=fps, bitrate=1800, codec="libx264",
                                            extra_args=["-pix_fmt", "yuv420p", "-crf", "28"]), dpi=dpi)
    plt.close(fig)
    return str(out)


if __name__ == "__main__":
    print(make())
