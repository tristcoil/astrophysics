"""The keeper animation: TIC 290061484's three suns, as seen from Earth.

A labelled simulation: REBOUND integrates the four stars from Kostov et al.
(2024)'s solution (or our refit of the outer orbits, when available); the
three B stars of the triple are drawn to scale in the plane of the sky,
nearer stars over farther ones; the light curve is computed from the
simulated positions (limb-darkened spheres; the fourth star, C, adds
constant light) and our TESS measurements appear as time passes.

    python animate.py
"""

from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, PillowWriter
from matplotlib.patches import Circle

import config as C
import eclipse_model as E
import etv
import plotstyle as PS
import tic290_nbody as N

PS.use_style()

SYS = C.SYSTEMS["TIC 290061484"]
RS = C.R_SUN_AU
NAMES = ("Aa", "Ab", "B")
#: TESS-band limb darkening for ~22,000 K, log g 4.3 (Claret 2017)
LD_TESS = (0.073, 0.207)


def simulate(t_start, t_end, dt, params=None):
    sim = N.build(params)
    sim.integrator = "ias15"
    if t_start > sim.t:
        sim.integrate(t_start)
    ts = np.arange(t_start, t_end, dt)
    xyz = np.empty((len(ts), 4, 3))
    for k, t in enumerate(ts):
        sim.integrate(t)
        for j in range(4):
            p = sim.particles[j]
            xyz[k, j] = (p.x, p.y, p.z)
    return ts, xyz


def light_curve(xyz):
    """Normalised flux of the four stars (C constant) from sky positions."""
    R = np.array([SYS["radii"][n] for n in NAMES]) * RS
    L = np.array([SYS["flux_tess"][n] for n in NAMES])
    LC = SYS["flux_tess"]["C"]
    u1, u2 = LD_TESS
    f = np.empty(len(xyz))
    for k in range(len(xyz)):
        tot = LC
        for i in range(3):
            hidden = 0.0
            for j in range(3):
                if j == i or xyz[k, j, 2] <= xyz[k, i, 2]:
                    continue                      # star j is not in front of star i
                d = np.hypot(*(xyz[k, j, :2] - xyz[k, i, :2]))
                if d < R[i] + R[j]:
                    hidden += E.occulted_fraction(d / R[i], R[j] / R[i], u1, u2)
            tot += L[i] * (1 - min(1.0, hidden))
        f[k] = tot / (L.sum() + LC)
    return f


def make(t_start=2458712.0, days=25.6, dt=0.02, stem="tic290061484_three_suns", fps=30, params=None,
         label="Kostov et al. (2024) solution"):
    ts, xyz = simulate(t_start, t_start + days, dt, params)
    flux = light_curve(xyz)
    # view centred on the triple's own barycentre (star C swings the triple
    # around the system's centre by ~600 solar radii)
    m = np.array([SYS["masses"][n] for n in NAMES])
    bary = (xyz[:, :3, :] * m[None, :, None]).sum(axis=1) / m.sum()
    xyz = xyz - bary[:, None, :]
    # observed: our photometry of sector 15 (and 16), detrended
    lcs = [x for x in etv.light_curves() if x["sector"] in (15, 16)]
    eph = json.loads((C.DERIVED / "tic290061484" / "ephemeris.json").read_text())
    obs_t = np.concatenate([x["jd"] for x in lcs])
    obs_f = np.concatenate([etv.detrend(x, eph["P"], eph["T0"]) for x in lcs])
    sel = (obs_t >= ts[0]) & (obs_t <= ts[-1])
    obs_t, obs_f = obs_t[sel], obs_f[sel]
    # the model's out-of-eclipse level matches the data's normalisation
    flux = flux / np.median(flux)
    R = np.array([SYS["radii"][n] for n in NAMES]) * RS
    col = [PS.STAR_COLOR[n] for n in NAMES]
    scale = 1 / RS                         # plot in solar radii
    fig = plt.figure(figsize=(12.8, 6.0), dpi=100)
    ax = fig.add_axes([0.04, 0.10, 0.40, 0.80])
    bx = fig.add_axes([0.52, 0.16, 0.45, 0.68])
    lim = 85
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim * 0.95, lim * 0.95)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("TIC 290061484 as seen from Earth", loc="left")
    # full paths, faint
    pair_c = xyz[:, :2, :].mean(axis=1)
    ax.plot(xyz[:, 2, 0] * scale, xyz[:, 2, 1] * scale, color=PS.CONTEXT, lw=0.8)
    ax.plot(pair_c[:, 0] * scale, pair_c[:, 1] * scale, color=PS.CONTEXT, lw=0.8)
    ax.plot([-75, -55], [-75, -75], color=PS.INK, lw=1.2)
    ax.text(-65, -71, "20 R☉", ha="center", fontsize=8)
    circles = [Circle((0, 0), R[i] * scale, color=col[i], ec=PS.INK, lw=0.4) for i in range(3)]
    for c in circles:
        ax.add_patch(c)
    texts = [ax.text(0, 0, n, fontsize=9, ha="left", va="bottom", color=PS.INK) for n in NAMES]
    clock = ax.text(0.02, 0.97, "", transform=ax.transAxes, fontsize=9, va="top")
    bx.set_xlim(0, days)
    lo = min(flux.min(), obs_f.min()) - 0.02
    bx.set_ylim(lo, 1.03)
    bx.set_xlabel(f"days from BJD {t_start:.1f} (2019 Aug)")
    bx.set_ylabel("brightness")
    bx.set_title("Light curve: simulation (line) and TESS (dots)", loc="left")
    model_line, = bx.plot([], [], color=PS.BLUE, lw=1.4)
    obs_pts, = bx.plot([], [], ".", color=PS.INK, ms=3)
    fig.text(0.52, 0.035, f"Labelled simulation: REBOUND, {label}.\nDots: our photometry of TESS sector 15 "
                          "(30-minute exposures).", fontsize=7.5, color=PS.MUTED)

    def draw(k):
        order = np.argsort(xyz[k, :3, 2])      # farthest first
        for z, i in enumerate(order):
            circles[i].center = (xyz[k, i, 0] * scale, xyz[k, i, 1] * scale)
            circles[i].set_zorder(3 + z)
            texts[i].set_position((xyz[k, i, 0] * scale + R[i] * scale + 1.0, xyz[k, i, 1] * scale + 1.0))
        clock.set_text(f"day {ts[k] - t_start:5.2f}")
        model_line.set_data(ts[:k + 1] - t_start, flux[:k + 1])
        m = obs_t <= ts[k]
        obs_pts.set_data(obs_t[m] - t_start, obs_f[m])
        return circles + texts + [clock, model_line, obs_pts]

    out = []
    w = FFMpegWriter(fps=fps, bitrate=4000, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "28"])
    mp4 = C.FIGURES / f"{stem}.mp4"
    with w.saving(fig, str(mp4), dpi=100):
        for k in range(len(ts)):
            draw(k)
            w.grab_frame()
    out.append(mp4)
    # a lighter GIF for the README: every 4th frame
    gif = C.FIGURES / f"{stem}.gif"
    gw = PillowWriter(fps=15)
    fig.set_dpi(70)
    with gw.saving(fig, str(gif), dpi=70):
        for k in range(0, len(ts), 4):
            draw(k)
            gw.grab_frame()
    out.append(gif)
    plt.close(fig)
    return [str(o) for o in out]


if __name__ == "__main__":
    print(make())
