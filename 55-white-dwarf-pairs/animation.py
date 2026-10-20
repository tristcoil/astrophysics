"""Animation (labelled visualisation): ZTF J1539+5027 through two orbits.

Left: the two white dwarfs to scale, as seen from Earth at Burdge et al.'s inclination
(84.15 deg), sizes and separation from their Table 1; the cool star's side facing its
hot companion is drawn brighter (the irradiation that dominates the light curve).
Right, top: the ZTF g-band light curve, binned data and our fitted model, with a cursor.
Right, bottom: the two polarisations of the gravitational wave reaching Earth
(h+ ~ A (1 + cos^2 i), hx ~ 2 A cos i), with A from our measured chirp mass at the
spectroscopic distance; they oscillate twice per orbit. Two orbits play in 24 s,
~35 x faster than real time; the real orbit lasts 6.91 minutes.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Circle, Polygon

import config as C
import gw
import plotstyle as S
import timing as T

S.use_style()
R_EARTH_KM = 6371.0


def make(res, art, fps=30, seconds=24, out_name="j1539_orbit.mp4"):
    P = res["joint_ztf_ptf"]["P_s"]
    m1, m2 = C.BURDGE19["M1"][0], C.BURDGE19["M2"][0]
    a_km = C.BURDGE19["a_rsun"][0] * C.R_SUN / 1e3
    r1_km = C.BURDGE19["R1_rsun"][0] * C.R_SUN / 1e3
    r2_km = C.BURDGE19["R2_rsun"][0] * C.R_SUN / 1e3
    inc = np.radians(C.BURDGE19["incl_deg"][0])
    a1, a2 = a_km * m2 / (m1 + m2), a_km * m1 / (m1 + m2)  # orbits about the centre of mass
    # data: binned g-band light curve from the fit
    fit, shapes, th = art["ztf_fit"], art["shapes"], art["joint_theta"]
    d = fit["data"][fit["data"]["band"] == "zg"]
    off = shapes["zg"]["offsets"]
    y = d["flux"].to_numpy() - d["group"].map(off).to_numpy()
    ph = T.phase(d["bjd_tdb"].to_numpy() - T.T_REF, th) % 1.0
    bins = np.linspace(0, 1, 41)
    idx = np.digitize(ph, bins) - 1
    w = 1 / d["err"].to_numpy() ** 2
    mb = np.array([np.sum(y[idx == k] * w[idx == k]) / np.sum(w[idx == k]) for k in range(40)])
    eb = np.array([1 / np.sqrt(np.sum(w[idx == k])) for k in range(40)])
    cb = 0.5 * (bins[1:] + bins[:-1])
    grid = np.linspace(0, 2, 2001)
    basis = T.shape_basis(grid, 4, False, True)
    model = basis @ shapes["zg"]["coef"]
    # gravitational wave, h+ and hx at Earth
    mc = res["chirp_mass"]["Mc_gr_only"]
    A = gw.amplitude(2 / P, mc, C.BURDGE19["d_kpc"][0] * C.KPC)
    hp = A * (1 + np.cos(inc) ** 2) * np.cos(4 * np.pi * grid)
    hx = 2 * A * np.cos(inc) * np.sin(4 * np.pi * grid)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    ax0 = fig.add_axes([0.03, 0.12, 0.44, 0.74])
    ax1 = fig.add_axes([0.55, 0.52, 0.42, 0.34])
    ax2 = fig.add_axes([0.55, 0.12, 0.42, 0.28])
    lim = 1.25 * a_km
    ax0.set_xlim(-lim, lim)
    ax0.set_ylim(-0.6 * lim, 0.6 * lim)
    ax0.set_aspect("equal")
    ax0.grid(False)
    for s in ax0.spines.values():
        s.set_visible(False)
    ax0.set_xticks([])
    ax0.set_yticks([])
    # scale bar: the Earth's diameter
    ax0.plot([-lim * 0.9, -lim * 0.9 + 2 * R_EARTH_KM], [-0.52 * lim] * 2, color=S.INK, lw=2.5)
    ax0.text(-lim * 0.9 + R_EARTH_KM, -0.49 * lim, "Earth's diameter", ha="center", fontsize=9, color=S.INK_SECONDARY)
    ax0.set_title("as seen from Earth (to scale)", fontsize=11)
    ax1.errorbar(cb, mb, eb, fmt="o", color=S.BLUE, ms=4.5, mec=S.SURFACE, mew=0.6, elinewidth=0.8,
                 label="ZTF g, 2018-2026 (binned)")
    ax1.errorbar(cb + 1, mb, eb, fmt="o", color=S.BLUE, ms=4.5, mec=S.SURFACE, mew=0.6, elinewidth=0.8)
    ax1.plot(grid, model, color=S.INK, lw=1.4, label="model")
    ax1.set_xlim(0, 2)
    ax1.set_ylabel("brightness (microjansky)")
    ax1.legend(loc="lower left", fontsize=8)
    ax1.set_title("what ZTF sees", fontsize=11)
    ax2.plot(grid, hp / 1e-22, color=S.ORANGE, lw=1.4, label="h$_+$")
    ax2.plot(grid, hx / 1e-22, color=S.AQUA, lw=1.4, label="h$_\\times$")
    ax2.set_xlim(0, 2)
    ax2.set_xlabel("orbits")
    ax2.set_ylabel("strain at Earth (10$^{-22}$)")
    ax2.legend(loc="lower left", fontsize=8, ncol=2)
    ax2.set_title("the gravitational wave LISA will hear (2 cycles per orbit)", fontsize=11)
    cur1 = ax1.axvline(0, color=S.INK_SECONDARY, lw=1)
    cur2 = ax2.axvline(0, color=S.INK_SECONDARY, lw=1)
    hot = Circle((0, 0), r1_km, color="#dbe8ff", ec="#8fb3e8", lw=1, zorder=3)
    cool = Circle((0, 0), r2_km, color="#c96d3a", ec="#a0522d", lw=1, zorder=2)
    glow = Polygon(np.zeros((4, 2)), closed=True, color="#ffc58a", zorder=2.5, alpha=0.95, lw=0)
    for p in (cool, glow, hot):
        ax0.add_patch(p)
    ttl = fig.suptitle("", x=0.03, ha="left", fontsize=13)
    fig.text(0.03, 0.05, "Visualisation of a model: sizes, separation and inclination from Burdge et al. (2019); light "
                         "curve and orbital period from this project's ZTF fit; wave amplitude from our chirp mass at "
                         "2.34 kpc.", fontsize=8, color=S.MUTED)
    n_frames = fps * seconds
    orbits = 2.0
    speed = orbits * P / seconds

    def draw(k):
        phi = orbits * k / n_frames  # orbital phase, 0 = primary eclipse (cool star in front)
        ang = 2 * np.pi * phi
        # observer looks along -y_orbit; x on sky = sin(ang), depth = cos(ang)
        x1, z1 = -a1 * np.sin(ang), -a1 * np.cos(ang)
        x2, z2 = a2 * np.sin(ang), a2 * np.cos(ang)
        y1, y2 = z1 * np.cos(inc), z2 * np.cos(inc)
        hot.center = (x1, y1)
        cool.center = (x2, y2)
        # the cool star's hemisphere facing the hot one is heated; we see it like a lunar
        # phase: phase angle alpha between the hot star and us, seen from the cool star
        cos_alpha = -np.cos(ang) * np.sin(inc)
        dx, dy = x1 - x2, y1 - y2
        e1 = np.array([dx, dy]) / max(np.hypot(dx, dy), 1e-9)  # towards the hot star (sky plane)
        e2 = np.array([-e1[1], e1[0]])
        t = np.linspace(np.pi / 2, -np.pi / 2, 60)
        limb = r2_km * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
        term = r2_km * (-cos_alpha * np.cos(t[::-1])[:, None] * e1 + np.sin(t[::-1])[:, None] * e2)
        glow.set_xy(np.vstack([limb, term]) + np.array([x2, y2]))
        glow.set_visible(cos_alpha > -0.985)
        # depth order: the nearer star drawn on top (z > 0 is towards the observer)
        hot.set_zorder(4 if z1 > z2 else 2)
        cool.set_zorder(3 if z2 > z1 else 1)
        glow.set_zorder(cool.get_zorder() + 0.5)
        cur1.set_xdata([phi, phi])
        cur2.set_xdata([phi, phi])
        ttl.set_text(f"ZTF J1539+5027: two white dwarfs, one orbit every {P / 60:.2f} minutes   "
                     f"t = {phi * P / 60:4.1f} min   (shown {speed:.0f}x faster than real time)")
        return hot, cool, glow, cur1, cur2, ttl

    out = C.FIGURE_DIR / out_name
    writer = FFMpegWriter(fps=fps, codec="libx264", bitrate=-1,
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "26", "-preset", "slow"])
    with writer.saving(fig, str(out), dpi=100):
        for k in range(n_frames):
            draw(k)
            writer.grab_frame()
    plt.close(fig)
    return str(out)
