"""All-sky maps in galactic coordinates: a Hammer-Aitoff projection on plain
matplotlib axes, so arrows can be drawn with directions that stay true.

Convention (as in Gaia's papers): l = 0, b = 0 at the centre, north up,
longitude increasing to the LEFT, as the sky looks from inside the sphere.

Proper-motion arrows: a vector with east (increasing l) and north (increasing
b) components is drawn at its point by stepping a fraction of a degree along
the great circle it points to, projecting both ends, and scaling the projected
direction to the arrow's length -- so arrows keep their true direction on the
sky even near the edges of the map, where meridians converge.
"""

from __future__ import annotations

import math

import numpy as np

SQ2 = math.sqrt(2.0)
X_MAX, Y_MAX = 2 * SQ2, SQ2


def wrap180(l_deg):
    return (np.asarray(l_deg, float) + 180.0) % 360.0 - 180.0


def project(l_deg, b_deg):
    """Galactic (l, b) in degrees -> Hammer map (x, y)."""
    lam = np.radians(-wrap180(l_deg))
    phi = np.radians(np.asarray(b_deg, float))
    den = np.sqrt(1.0 + np.cos(phi) * np.cos(lam / 2))
    return 2 * SQ2 * np.cos(phi) * np.sin(lam / 2) / den, SQ2 * np.sin(phi) / den


def unproject(x, y):
    """Hammer map (x, y) -> (l, b) degrees; NaN outside the map."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    inside = (x / X_MAX) ** 2 + (y / Y_MAX) ** 2 <= 1.0
    z = np.sqrt(np.clip(1 - (x / 4) ** 2 - (y / 2) ** 2, 0, None))
    lam = 2 * np.arctan2(z * x, 2 * (2 * z**2 - 1))
    phi = np.arcsin(np.clip(z * y, -1, 1))
    l = (-np.degrees(lam)) % 360.0
    b = np.degrees(phi)
    return np.where(inside, l, np.nan), np.where(inside, b, np.nan)


def unit(l_deg, b_deg) -> np.ndarray:
    l, b = np.radians(l_deg), np.radians(b_deg)
    return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1)


def lonlat(u: np.ndarray):
    return (np.degrees(np.arctan2(u[..., 1], u[..., 0])) % 360.0,
            np.degrees(np.arcsin(np.clip(u[..., 2], -1, 1))))


def frame(ax, *, labels: bool = True, grid_color="#e1e0d9", edge_color="#c3c2b7",
          label_color="#898781", fontsize=7) -> None:
    """Outline, graticule (every 30 deg) and longitude labels."""
    ax.set_aspect("equal")
    ax.set_xlim(-X_MAX * 1.02, X_MAX * 1.02)
    ax.set_ylim(-Y_MAX * 1.04, Y_MAX * 1.04)
    ax.axis("off")
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(X_MAX * np.cos(t), Y_MAX * np.sin(t), color=edge_color, lw=0.8, zorder=1)
    bb = np.linspace(-90, 90, 181)
    for l in range(-150, 181, 30):
        lv = np.full_like(bb, l if l != 180 else 179.999)
        ax.plot(*project(lv, bb), color=grid_color, lw=0.5, zorder=0.5)
    ll = np.linspace(-179.999, 179.999, 361)
    for b in (-60, -30, 0, 30, 60):
        ax.plot(*project(ll, np.full_like(ll, b)), color=grid_color, lw=0.5, zorder=0.5)
    if labels:
        for l in (90, 180 - 1e-6, 270, 0):
            x, y = project(l, 0)
            txt = {0: "l = 0°", 90: "90°", 270: "270°"}.get(int(round(l)) % 360, "180°")
            if abs(l - 180) < 1:
                x = X_MAX * 0.985
                ax.text(x + 0.05, 0.02, "180°", fontsize=fontsize, color=label_color, ha="left", va="bottom")
                ax.text(-x - 0.05, 0.02, "180°", fontsize=fontsize, color=label_color, ha="right", va="bottom")
                continue
            ax.text(x, y - 0.05, txt, fontsize=fontsize, color=label_color, ha="center", va="top")
        for b in (-60, -30, 30, 60):
            x, y = project(0, b)
            ax.text(x + 0.03, y, f"{b:+d}°", fontsize=fontsize - 0.5, color=label_color, va="center")


def raster(ax, value_of_lb, *, n=900, **imshow_kw):
    """Fill the map with value_of_lb(l, b) (vectorised), e.g. a HEALPix lookup."""
    xs = np.linspace(-X_MAX, X_MAX, n)
    ys = np.linspace(-Y_MAX, Y_MAX, n // 2)
    X, Y = np.meshgrid(xs, ys)
    L, B = unproject(X, Y)
    ok = np.isfinite(L)
    val = np.full(L.shape, np.nan)
    val[ok] = value_of_lb(L[ok], B[ok])
    return ax.imshow(np.ma.masked_invalid(val), origin="lower",
                     extent=(-X_MAX, X_MAX, -Y_MAX, Y_MAX), interpolation="bilinear", **imshow_kw)


def arrows(ax, l_deg, b_deg, east, north, *, scale: float, step_deg: float = 0.3, **quiver_kw):
    """Draw tangent vectors (east = +l, north = +b) of any units; ``scale`` is
    the map length of one unit."""
    l_deg, b_deg = np.asarray(l_deg, float), np.asarray(b_deg, float)
    east, north = np.asarray(east, float), np.asarray(north, float)
    amp = np.hypot(east, north)
    u = unit(l_deg, b_deg)
    lr, br = np.radians(l_deg), np.radians(b_deg)
    e_l = np.stack([-np.sin(lr), np.cos(lr), np.zeros_like(lr)], -1)
    e_b = np.stack([-np.sin(br) * np.cos(lr), -np.sin(br) * np.sin(lr), np.cos(br)], -1)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = (east[:, None] * e_l + north[:, None] * e_b) / amp[:, None]
    eps = math.radians(step_deg)
    u2 = math.cos(eps) * u + math.sin(eps) * d
    x0, y0 = project(l_deg, b_deg)
    x1, y1 = project(*lonlat(u2))
    dx, dy = x1 - x0, y1 - y0
    norm = np.hypot(dx, dy)
    ok = np.isfinite(norm) & (norm > 0) & (norm < 0.5) & np.isfinite(amp)
    ux, uy = np.zeros_like(dx), np.zeros_like(dy)
    ux[ok], uy[ok] = dx[ok] / norm[ok], dy[ok] / norm[ok]
    kw = {"angles": "xy", "scale_units": "xy", "scale": 1.0, "width": 0.0028,
          "headwidth": 3.6, "headlength": 4.2, "headaxislength": 3.8, "pivot": "tail", "zorder": 5}
    kw.update(quiver_kw)
    return ax.quiver(x0[ok], y0[ok], (amp * scale * ux)[ok], (amp * scale * uy)[ok], **kw)


def ellipse_points(centre_lb, cov_lb_deg, nsig: float = 1.0, n: int = 200):
    """Points of a covariance ellipse around (l, b) given in degrees, with the
    longitude error as a coordinate error (scaled by cos b on the sky)."""
    l0, b0 = centre_lb
    vals, vecs = np.linalg.eigh(np.asarray(cov_lb_deg))
    t = np.linspace(0, 2 * np.pi, n)
    circle = np.stack([np.cos(t), np.sin(t)])
    pts = vecs @ (nsig * np.sqrt(np.maximum(vals, 0))[:, None] * circle)
    return l0 + pts[0], b0 + pts[1]
