"""Every figure in the project, from ``analysis.run_all()`` output.

    01  S2 on the sky: 27 years of positions from two observatories
    02  The radial velocity curve and the 2018 pericentre swing
    03  Mass and distance: the "banana", and the published values
    04  Einstein's test: the gravitational redshift, and the precession these data cannot see
    05  The scale: S2's orbit, the Solar System, and the event horizon
    06  Why a black hole: how long could a dark cluster survive?
    07  A prediction for the Event Horizon Telescope: the shadow of Sgr A*
    08  Grading our code against both teams' own fits

Each figure is saved as PNG and SVG in ``figures/``. To redraw them without
re-running any fit: ``python -c "import analysis, figures; figures.make_all(*analysis.load_cached())"``.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patches  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402

import black_hole_case as bh  # noqa: E402
import orbit_physics as op  # noqa: E402
import plotstyle as ps  # noqa: E402
import s2data as sd  # noqa: E402
from config import FIGURE_DIR, PUBLISHED  # noqa: E402

ps.use_style()

PHYS_KEYS = ("gm6", "r0", "period", "t0", "e", "inc", "omega", "Omega")
#: viridis without its palest yellows, so the newest points stay visible on a light page
TIME_CMAP = ListedColormap(plt.get_cmap(ps.TIME_CMAP)(np.linspace(0.0, 0.86, 256)))
SUPERSCRIPT = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def best_params(art: dict, name: str) -> dict:
    fit = art["fits"][name]
    model = art["models"][name]
    p = dict(zip(fit["names"], fit["best_theta"]))
    for k, v in model.fixed.items():
        p.setdefault(k, v)
    return p


def sky_track(p: dict, t, *, roemer=True) -> tuple[np.ndarray, np.ndarray]:
    """Model sky offsets of S2 from the black hole (mas, x = East), 1-D arrays."""
    pred = op.predict({k: p[k] for k in (*PHYS_KEYS, "f", "f_sp")}, np.atleast_1d(t), None, roemer=roemer)
    return np.ravel(pred["x"]), np.ravel(pred["y"])


def sky_point(p: dict, t: float) -> tuple[float, float]:
    x, y = sky_track(p, [t])
    return float(x[0]), float(y[0])


def total_err(v) -> float:
    """Quadrature sum of a published (value, err1, err2, ...) tuple's errors."""
    return float(np.sqrt(np.sum(np.square(v[1:]))))


def rv_track(p: dict, t, *, f=None, roemer=True) -> np.ndarray:
    q = {k: p[k] for k in (*PHYS_KEYS, "f", "f_sp")}
    if f is not None:
        q["f"] = f
    return np.ravel(op.predict(q, None, np.atleast_1d(t), roemer=roemer)["rv"]) + p.get("vz0", 0.0)


def relative_positions(frame, p: dict):
    """Data positions with that team's fitted frame offset and drift removed."""
    k = frame.key
    dt = frame.data.t - frame.t_ref
    return (frame.data.x - p[f"x0_{k}"] - p[f"vx0_{k}"] * dt,
            frame.data.y - p[f"y0_{k}"] - p[f"vy0_{k}"] * dt)


def hpd_levels(h: np.ndarray, fractions=(0.68, 0.95)) -> list[float]:
    """Density thresholds enclosing the given probability fractions."""
    flat = np.sort(h.ravel())[::-1]
    cum = np.cumsum(flat)
    cum /= cum[-1]
    return [float(flat[np.searchsorted(cum, f)]) for f in fractions]


def contour2d(ax, x, y, *, weights=None, color, lw=1.4, ls="-", bins=80, smooth=1.6,
              fill=True, zorder=3, extent=None):
    """68% and 95% credible contours of a 2D sample (optionally weighted)."""
    h, xe, ye = np.histogram2d(x, y, bins=bins, range=extent, weights=weights)
    h = gaussian_filter(h, smooth)
    lv = hpd_levels(h)
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    if fill:
        ax.contourf(xc, yc, h.T, levels=[lv[0], h.max() * 1.01], colors=[color], alpha=0.16,
                    zorder=zorder - 1)
    ax.contour(xc, yc, h.T, levels=[lv[1], lv[0]], colors=[color], linewidths=[lw * 0.7, lw],
               linestyles=ls, zorder=zorder)


def save(fig, name: str) -> list[str]:
    return ps.save(fig, FIGURE_DIR / name)


def fmt_pm(s: dict, digits=2) -> str:
    lo, hi, med = s["lo68"], s["hi68"], s["median"]
    return f"{med:.{digits}f} ± {0.5 * (hi - lo):.{digits}f}"


def sci(x: float, digits: int = 1) -> str:
    """6.3e7 -> '6.3 × 10⁷' (readable powers of ten for labels)."""
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10.0**e
    if round(m, digits) >= 10:
        m, e = m / 10.0, e + 1
    return f"{m:.{digits}f} × 10{str(e).translate(SUPERSCRIPT)}"


def leader(ax, text, xy, xytext, *, color=ps.INK_SECONDARY, fontsize=8, **kw):
    return ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
                       arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.6, "shrinkA": 2, "shrinkB": 3},
                       **kw)


# --------------------------------------------------------------------------
# 01  The orbit on the sky
# --------------------------------------------------------------------------


def fig_sky(res: dict, art: dict) -> list[str]:
    p = best_params(art, "joint_gr")
    model = art["models"]["joint_gr"]
    d = res["fits"]["joint_gr"]["derived"]
    fig = plt.figure(figsize=(12.0, 7.6))
    ax = fig.add_axes([0.06, 0.09, 0.46, 0.83])
    az = fig.add_axes([0.60, 0.36, 0.37, 0.56])
    cax = fig.add_axes([0.62, 0.20, 0.33, 0.022])
    norm = plt.Normalize(1992, 2019)
    markers = {"keck": "o", "vlt": "s"}

    def draw_data(axis, sel_fn, size_keck=24, size_vlt=15):
        for fr in model.frames:
            xr, yr = relative_positions(fr, p)
            sel = sel_fn(fr)
            confused = (fr.key == "vlt") & (fr.data.t > 2002.0) & (fr.data.t < 2003.0)
            axis.errorbar(xr[sel], yr[sel], xerr=fr.data.x_err[sel], yerr=fr.data.y_err[sel], fmt="none",
                          ecolor=ps.AXIS, elinewidth=0.6, zorder=3)
            for mask, alpha in ((sel & ~confused, 1.0), (sel & confused, 0.35)):
                if np.any(mask):
                    axis.scatter(xr[mask], yr[mask], c=fr.data.t[mask], cmap=TIME_CMAP, norm=norm,
                                 marker=markers[fr.key], s=size_keck if fr.key == "keck" else size_vlt,
                                 edgecolors=ps.SURFACE, linewidths=0.6, zorder=4, alpha=alpha)

    # --- the whole orbit -------------------------------------------------
    t_orb = np.linspace(p["t0"] - p["period"], p["t0"], 6000)
    xo, yo = sky_track(p, t_orb)
    ax.plot(xo, yo, color=ps.INK, lw=0.9, zorder=2)
    draw_data(ax, lambda fr: np.ones(fr.data.n, bool))
    ax.scatter([0], [0], s=80, facecolor=ps.INK, edgecolor=ps.SURFACE, lw=1.5, zorder=6)
    xp, yp = sky_point(p, p["t0"])
    xa, ya = sky_point(p, p["t0"] - 0.5 * p["period"])
    xc, yc = 0.5 * (xp + xa), 0.5 * (yp + ya)  # projection keeps midpoints: the ellipse's centre
    ax.scatter([xc], [yc], marker="+", s=70, color=ps.INK_SECONDARY, lw=1.1, zorder=6)
    leader(ax, "Sgr A*: the mass sits\nat a focus of the ellipse", (0, 0), (-8, 44), fontsize=8.5, color=ps.INK)
    ax.text(xc - 4, yc + 4, "the ellipse's centre:\nnothing there", fontsize=8, color=ps.INK_SECONDARY,
            va="bottom")
    rp_au, ra_au = d["pericentre_au"]["median"], d["apocentre_au"]["median"]
    leader(ax, f"pericentre (2002.3 and 2018.4):\n{rp_au:.0f} AU from Sgr A*", (xp, yp), (42, -38))
    leader(ax, f"apocentre ({p['t0'] - 0.5 * p['period']:.1f}): {ra_au:,.0f} AU", (xa, ya), (-18, 196))
    k = np.searchsorted(t_orb, 2012.3)
    ax.annotate("", xy=(xo[k + 70], yo[k + 70]), xytext=(xo[k], yo[k]),
                arrowprops={"arrowstyle": "-|>", "color": ps.INK, "lw": 0.9}, zorder=5)
    ax.text(-74, 166, "moving\nclockwise", fontsize=8, color=ps.INK_SECONDARY)
    bar = 500.0 / p["r0"]  # mas
    x_bar = -45.0
    ax.plot([x_bar, x_bar - bar], [-52, -52], color=ps.INK, lw=1.5)
    ax.text(x_bar - bar / 2, -48, f"500 AU = {bar:.0f} mas", ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlim(78, -112)
    ax.set_ylim(-62, 206)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("offset in right ascension (mas; East is left, as on the sky)")
    ax.set_ylabel("offset in declination (mas, North up)")
    ax.set_title("S2 around Sagittarius A*, 1992–2018")
    n_ep = {fr.key: fr.data.n for fr in model.frames}
    handles = [Line2D([], [], marker="o", ls="", color=ps.MUTED, ms=5, label=f"Keck (UCLA): {n_ep['keck']} nights"),
               Line2D([], [], marker="s", ls="", color=ps.MUTED, ms=4.5,
                      label=f"NTT + VLT (MPE): {n_ep['vlt']} nights"),
               Line2D([], [], color=ps.INK, lw=0.9, label="best-fitting orbit, both teams")]
    ax.legend(handles=handles, loc="upper left", fontsize=8)

    # --- zoom: the 2018 pericentre ----------------------------------------
    t_in = np.linspace(2016.9, 2018.95, 3000)
    xi, yi = sky_track(p, t_in)
    az.plot(xi, yi, color=ps.INK, lw=0.9, zorder=2)
    draw_data(az, lambda fr: (fr.data.t > 2017.0) & (fr.data.t < 2019.0), size_keck=34)
    az.scatter([0], [0], s=90, facecolor=ps.INK, edgecolor=ps.SURFACE, lw=1.5, zorder=6)
    az.text(1.5, 1.5, "Sgr A*", fontsize=8.5, color=ps.INK, ha="right", va="bottom")
    for tt, lab, off in ((2017.35, "May 2017", (-8, 2)), (2017.63, "Aug 2017", (-9, -3)),
                         (2018.21, "Mar 2018", (-3, -8)), (p["t0"], "pericentre\n18 May 2018", (9, -9)),
                         (2018.67, "Sep 2018\n(model)", (8, 4))):
        xx, yy = sky_point(p, tt)
        leader(az, lab, (xx, yy), (xx + off[0], yy + off[1]), fontsize=7.5)
    az.set_xlim(40, -72)
    az.set_ylim(-22, 48)
    az.set_aspect("equal", adjustable="datalim")
    az.set_xlabel("offset in right ascension (mas)")
    az.set_title("zoom: the closest approach, 2017–2018", fontsize=10)
    swing = abs(np.degrees(np.arctan2(*sky_point(p, 2018.21)) - np.arctan2(*sky_point(p, 2018.67))))
    az.text(38, 46, f"Between March and September 2018 S2 swung {swing:.0f}°\n"
            f"around the black hole, passing {d['pericentre_light_hours']['median']:.0f} light-hours\n"
            f"from it at {d['pericentre_speed_kms']['median']:,.0f} km/s.", fontsize=8, color=ps.INK_SECONDARY,
            va="top")

    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=TIME_CMAP), cax=cax, orientation="horizontal")
    cb.set_label("year of observation", color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    ps.caption(fig, "Positions relative to the fitted position of the mass: each team's reference frame gets its own offset "
                    "and drift. Faded squares: the ten 2002 NTT/VLT nights, when S2 was blended with another star near "
                    "pericentre (their published errors are inflated ~10×). Data: Do et al. 2019 (Keck, via the Internet "
                    "Archive); Gillessen et al. 2017 (NTT/VLT, CDS J/ApJ/837/30).")
    return save(fig, "01_sky_orbit")


# --------------------------------------------------------------------------
# 02  Radial velocities
# --------------------------------------------------------------------------


def fig_rv(res: dict, art: dict) -> list[str]:
    p = best_params(art, "joint_gr")
    model = art["models"]["joint_gr"]
    rv = model.rv
    is_vlt = rv.instrument == "VLT"
    corr = np.where(rv.instrument == "NIRC2", p.get("nirc2_offset", 0.0), 0.0)
    v_plot = rv.v + corr  # apply the fitted NIRC2 correction for display
    t = np.linspace(1999.5, 2020.0, 8000)
    vm = rv_track(p, t)
    resid = v_plot - rv_track(p, rv.t)

    fig = plt.figure(figsize=(11.0, 7.2))
    ax = fig.add_axes([0.08, 0.36, 0.60, 0.56])
    axr = fig.add_axes([0.08, 0.10, 0.60, 0.20], sharex=ax)
    axz = fig.add_axes([0.74, 0.36, 0.24, 0.56])
    axzr = fig.add_axes([0.74, 0.10, 0.24, 0.20], sharex=axz)

    for a, ar in ((ax, axr), (axz, axzr)):
        a.plot(t, vm, color=ps.INK, lw=1.0, zorder=2)
        ar.axhline(0, color=ps.INK, lw=0.8)
        for sel, col, lab in ((~is_vlt, ps.C_KECK, "Keck, Subaru, Gemini (UCLA)"), (is_vlt, ps.C_VLT, "VLT (MPE)")):
            a.errorbar(rv.t[sel], v_plot[sel], yerr=rv.v_err[sel], fmt="o", ms=3.2, color=col,
                       ecolor=col, elinewidth=0.7, mec=ps.SURFACE, mew=0.5, zorder=3, label=lab)
            ar.errorbar(rv.t[sel], resid[sel], yerr=rv.v_err[sel], fmt="o", ms=3.0, color=col,
                        ecolor=col, elinewidth=0.7, mec=ps.SURFACE, mew=0.5, zorder=3)
    ax.set_xlim(1999.5, 2019.8)
    ax.set_ylim(-2500, 4700)
    axr.set_ylim(-160, 160)
    axz.set_xlim(2018.1, 2018.8)
    axz.set_ylim(-2500, 4700)
    axzr.set_ylim(-160, 160)
    ax.set_ylabel("line-of-sight velocity (km/s, + = moving away)")
    axr.set_ylabel("data − model\n(km/s)")
    axr.set_xlabel("year")
    axzr.set_xlabel("2018")
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(axz.get_xticklabels(), visible=False)
    axz.set_yticklabels([])
    axzr.set_yticklabels([])
    axzr.set_xticks([2018.17, 2018.33, 2018.5, 2018.67])
    axzr.set_xticklabels(["Mar", "May", "Jul", "Sep"])
    ax.set_title(f"S2's line-of-sight velocity: {rv.n} measurements, 2000–2019")
    axz.set_title("the 2018 pericentre swing", fontsize=9.5, loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.48, 1.0), fontsize=8)

    d = res["fits"]["joint_gr"]["derived"]
    in18 = (t > 2017.8) & (t < 2018.9)
    k = np.flatnonzero(in18)[np.argmax(vm[in18])]
    k2 = np.flatnonzero(in18)[np.argmin(vm[in18])]
    leader(ax, f"+{vm[k]:,.0f} km/s:\nrushing away from us", (t[k], vm[k]), (2011.3, 3650))
    leader(ax, f"{vm[k2]:,.0f} km/s: rushing towards us", (t[k2], vm[k2]), (2008.8, -2150))
    leader(ax, "2002: the previous\npericentre", (2002.37, 900), (2003.4, 2300))
    axz.text(2018.52, 3300, f"{vm[k] - vm[k2]:,.0f} km/s\nswing in\n{(t[k2] - t[k]) * 12:.1f} months",
             fontsize=8, color=ps.INK_SECONDARY, va="top")
    axz.text(2018.12, -1150, f"top speed {d['pericentre_speed_kms']['median']:,.0f} km/s\n"
             f"= {100 * d['pericentre_speed_frac_c']['median']:.1f}% of light speed", fontsize=7.5,
             color=ps.INK_SECONDARY, va="top")
    rms = float(np.sqrt(np.mean(resid**2)))
    chi = res["fits"]["joint_gr"]["chi2"]
    axr.text(1999.8, 135, f"rms {rms:.0f} km/s; χ² {chi['rv']:.0f} for {chi['n_rv']} velocities",
             fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    ps.caption(fig, "Model: the joint Keck + VLT orbit with general relativity (light-travel delay, relativistic redshift, "
                    "precession). Velocities are relative to the Local Standard of Rest; the six Keck/NIRC2 points "
                    f"(2002–2005) include the fitted instrumental correction (+{p.get('nirc2_offset', 0):.0f} km/s, as Do et "
                    "al. found). Data: Do et al. 2019, whose file includes the VLT velocities of Gillessen et al. 2017.")
    return save(fig, "02_radial_velocity")


# --------------------------------------------------------------------------
# 03  Mass and distance
# --------------------------------------------------------------------------


def fig_mass_distance(res: dict, art: dict) -> list[str]:
    fits = art["fits"]
    fig, ax = plt.subplots(figsize=(9.0, 6.8))
    extent = [[7.3, 9.0], [3.5, 5.0]]
    for name, col in (("vlt_gr", ps.C_VLT), ("keck_gr", ps.C_KECK), ("joint_gr", ps.C_JOINT)):
        f = fits[name]
        names = list(f["names"])
        contour2d(ax, f["chain"][:, names.index("r0")], f["chain"][:, names.index("gm6")], color=col,
                  extent=extent, bins=170, smooth=2.0, lw=1.6 if name == "joint_gr" else 1.2,
                  fill=name != "keck_gr")

    j = res["fits"]["joint_gr"]
    alpha = j["mass_distance_exponent"]
    r_ref, m_ref = j["parameters"]["r0"]["median"], j["parameters"]["gm6"]["median"]
    rr = np.linspace(*extent[0], 50)
    ax.plot(rr, m_ref * (rr / r_ref) ** alpha, color=ps.MUTED, lw=0.8, zorder=1)
    ax.text(8.93, m_ref * (8.93 / r_ref) ** alpha + 0.03, f"M ∝ R₀^{alpha:.2f}", color=ps.INK_SECONDARY,
            fontsize=8, ha="right", va="bottom", rotation=31, rotation_mode="anchor")
    ax.plot(rr, m_ref * (rr / r_ref) ** 3, color=ps.GRID, lw=0.9, zorder=1)
    ax.text(7.645, m_ref * (7.645 / r_ref) ** 3 + 0.015, "positions alone: M ∝ R₀³", color=ps.MUTED, fontsize=7,
            rotation=41, rotation_mode="anchor", ha="left", va="bottom")

    pubs = [("do2019_gr", "Do+ 2019 (Keck)", (7.45, 3.72)),
            ("gravity2018", "GRAVITY 2018", (8.36, 3.86)),
            ("gravity2019", "GRAVITY 2019", (8.46, 3.97)),
            ("gravity2020", "GRAVITY 2020", (8.56, 4.08)),
            ("gravity2022", "GRAVITY 2022", (8.64, 4.20)),
            ("g17_multistar", "Gillessen+ 2017\n(17 stars)", (8.52, 4.72))]
    for key, lab, xy in pubs:
        v = PUBLISHED[key]
        r, rerr = v["r0"][0], total_err(v["r0"])
        m, merr = v["gm6"][0], total_err(v["gm6"])
        ax.errorbar(r, m, xerr=rerr, yerr=merr, fmt="D", ms=4, color=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED,
                    elinewidth=0.8, capsize=0, mec=ps.SURFACE, mew=0.6, zorder=6)
        leader(ax, lab, (r, m), xy, fontsize=7.5)

    mg = res["headline"]["mass_at_gravity_distance"]["joint_gr"]
    ax.text(7.33, 4.93, "Slide along our banana to GRAVITY's distance, 8.277 kpc:\n"
            f"M = {mg['gm6']:.3f} ± {mg['gm6_err']:.3f} million suns (GRAVITY: "
            f"{mg['gravity2022_gm6']} ± {mg['gravity2022_gm6_err']})", fontsize=8, color=ps.INK, va="top")
    ax.set_xlim(*extent[0])
    ax.set_ylim(*extent[1])
    ax.set_xlabel("distance to the Galactic Centre, R₀ (kpc)")
    ax.set_ylabel("mass of Sgr A* (million solar masses)")
    ax.set_title("Weighing Sgr A*, and measuring the distance to it, with one star")

    def lab(name, text):
        q = res["fits"][name]["parameters"]
        return f"{text}: M = {fmt_pm(q['gm6'])}, R₀ = {fmt_pm(q['r0'])} kpc"

    handles = [Line2D([], [], color=ps.C_KECK, lw=1.2, label=lab("keck_gr", "Keck alone")),
               Line2D([], [], color=ps.C_VLT, lw=1.2, label=lab("vlt_gr", "VLT alone")),
               Line2D([], [], color=ps.C_JOINT, lw=1.6, label=lab("joint_gr", "Keck + VLT")),
               Line2D([], [], marker="D", ls="", color=ps.C_PUBLISHED, ms=4, label="published values (±1σ)")]
    leg = ax.legend(handles=handles, loc="lower right", fontsize=7.8,
                    title="our fits, general relativity (68% and 95% regions)", title_fontsize=8)
    leg._legend_box.align = "left"
    ps.caption(fig, "Why a banana: the positions measure the orbit's size as an angle, the velocities in km/s. Put the star "
                    "further away and its orbit must be bigger, so the black hole heavier; the velocities push back, "
                    "leaving M ∝ R₀^~1.9. Each team's reference frame gets its own offset and drift. GRAVITY's "
                    "interferometric data (2017 onwards) are not public; its distance sits 4% beyond the Keck-based value, "
                    "a known disagreement between the teams.")
    return save(fig, "03_mass_distance")


# --------------------------------------------------------------------------
# 04  The relativity tests
# --------------------------------------------------------------------------


def fig_redshift(res: dict, art: dict) -> list[str]:
    p = best_params(art, "joint_redshift")
    model = art["models"]["joint_redshift"]
    rv = model.rv
    corr = np.where(rv.instrument == "NIRC2", p.get("nirc2_offset", 0.0), 0.0)
    v_plot = rv.v + corr
    is_vlt = rv.instrument == "VLT"
    t = np.linspace(2017.0, 2019.6, 6000)
    newton = rv_track(p, t, f=0.0)
    einstein = rv_track(p, t, f=1.0) - newton
    best = rv_track(p, t, f=p["f"]) - newton
    resid = v_plot - rv_track(p, rv.t, f=0.0)

    fig = plt.figure(figsize=(11.5, 5.8))
    ax = fig.add_axes([0.07, 0.14, 0.50, 0.74])
    axf = fig.add_axes([0.65, 0.58, 0.33, 0.30])
    axp = fig.add_axes([0.65, 0.14, 0.33, 0.30])

    ax.axhline(0, color=ps.C_NEWTON, lw=1.2, ls="--", zorder=1)
    ax.plot(t, einstein, color=ps.C_GR, lw=2.0, zorder=2)
    ax.plot(t, best, color=ps.INK, lw=0.9, zorder=2)
    sel = (rv.t > 2017.0) & (rv.t < 2019.6)
    for s, col, lab in ((sel & ~is_vlt, ps.C_KECK, "Keck, Subaru, Gemini"), (sel & is_vlt, ps.C_VLT, "VLT")):
        if np.any(s):
            ax.errorbar(rv.t[s], resid[s], yerr=rv.v_err[s], fmt="o", ms=4, color=col, ecolor=col, elinewidth=0.8,
                        mec=ps.SURFACE, mew=0.6, zorder=3, label=lab)
    ax.set_xlim(2017.0, 2019.6)
    ax.set_ylim(-150, 300)
    ax.set_xlabel("year")
    ax.set_ylabel("measured velocity − Newtonian orbit (km/s)")
    ax.set_title("Clocks run slow near a black hole: the 2018 pericentre")
    k = np.argmax(einstein)
    ax.text(t[k] + 0.04, einstein[k] + 6, f"Einstein: +{einstein[k]:.0f} km/s at pericentre\n"
            "(time dilation + gravitational redshift)", color=ps.INK_SECONDARY, fontsize=8, va="bottom")
    ax.text(2017.03, -9, "Newton: no shift", color=ps.MUTED, fontsize=8, va="top")
    rt = res["relativity"]["joint_redshift"]
    kb = np.argmax(best)
    leader(ax, f"best fit: f = {rt['f_median']:.2f} ± {rt['f_sd']:.2f}", (t[kb] + 0.08, best[kb] * 0.62),
           (2018.72, 150), fontsize=8, color=ps.INK)
    ax.text(2017.03, -95, f"Newton ruled out at {rt['sigma_from_newton']:.1f}σ; Einstein "
            f"{abs(rt['sigma_from_einstein']):.1f}σ away", color=ps.INK, fontsize=8.5)
    ax.legend(loc="upper left", fontsize=8)

    # Posterior of f.
    bins = np.linspace(-0.2, 1.8, 70)
    fits = art["fits"]
    for name, col, lab in (("keck_redshift", ps.C_KECK, "Keck positions"),
                           ("joint_redshift", ps.C_JOINT, "Keck + VLT positions")):
        names = list(fits[name]["names"])
        h, e = np.histogram(fits[name]["chain"][:, names.index("f")], bins=bins, density=True)
        axf.stairs(h, e, color=col, lw=1.3, label=lab)
    do = sd.load_do2019_chain("redshift")
    h, e = np.histogram(do["f"], bins=bins, weights=do["weight"], density=True)
    axf.stairs(h, e, color=ps.C_PUBLISHED, lw=1.0, label="Do+ 2019's chain")
    axf.set_ylim(0, axf.get_ylim()[1] * 1.45)  # headroom for the legend
    axf.axvline(0, color=ps.C_NEWTON, ls="--", lw=1.0)
    axf.axvline(1, color=ps.C_GR, lw=1.4)
    g = PUBLISHED["gravity2018"]["f"]
    y_g = axf.get_ylim()[1] * 0.10
    axf.errorbar(g[0], y_g, xerr=total_err(g), fmt="D", ms=4, color=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED,
                 elinewidth=0.8, mec=ps.SURFACE)
    axf.text(g[0] + total_err(g) + 0.04, y_g, "GRAVITY 2018", fontsize=7, color=ps.INK_SECONDARY, va="center")
    axf.set_yticks([])
    axf.spines["left"].set_visible(False)
    axf.set_xlabel("redshift factor f (0 = Newton, 1 = Einstein)", fontsize=8.5)
    axf.set_title("how much redshift the data want", fontsize=9.5)
    axf.legend(fontsize=7, loc="upper left", handlelength=1.2)

    # The precession factor: a profile likelihood, all positions vs adaptive-optics era only.
    prof = res.get("precession_profile")
    g = PUBLISHED["gravity2020"]["f_sp"]
    if prof:
        for key, col, lw, lab in (("all_positions", ps.C_JOINT, 1.5, "all positions, 1992–2018"),
                                  ("ao_only", ps.C_KECK, 1.3, "adaptive-optics era only")):
            pr_ = prof[key]
            axp.plot(pr_["f_sp"], pr_["delta_chi2"], color=col, lw=lw, marker="o", ms=3, label=lab)
        top = max(prof["all_positions"]["delta_chi2"])
        axp.set_ylim(-1.5, top * 1.08)
        axp.set_ylabel("Δχ² from best", fontsize=8)
        axp.errorbar(g[0], 0.12 * top, xerr=g[1], fmt="D", ms=4, color=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED,
                     elinewidth=0.8, mec=ps.SURFACE)
        axp.text(g[0] + g[1] + 0.2, 0.12 * top, "GRAVITY 2020", fontsize=7, color=ps.INK_SECONDARY, va="center")
        axp.legend(fontsize=7, loc="upper right", handlelength=1.5)
        fbest = prof["all_positions"]["f_sp"][int(np.argmin(prof["all_positions"]["delta_chi2"]))]
        axp.set_title(f"the orbit's turning: old positions fake {fbest:.0f}× Einstein", fontsize=9.5)
    else:  # fall back to the MCMC posterior
        names = list(fits["joint_precession"]["names"])
        h, e = np.histogram(fits["joint_precession"]["chain"][:, names.index("f_sp")], bins=np.linspace(-1.5, 9, 70),
                            density=True)
        axp.stairs(h, e, color=ps.C_JOINT, lw=1.3)
        axp.set_yticks([])
        axp.set_title("the orbit's turning (posterior)", fontsize=9.5)
    axp.axvline(0, color=ps.C_NEWTON, ls="--", lw=1.0)
    axp.axvline(1, color=ps.C_GR, lw=1.4)
    axp.set_xlim(-2.4, 10.4)
    axp.set_xlabel("precession factor (0 = Newton, 1 = Einstein: 12′ per orbit)", fontsize=8.5)
    ps.caption(fig, "Left: velocities near the 2018 pericentre minus the same orbit computed without the relativistic "
                    "redshift. Top right: the redshift factor f, with Do et al.'s own chain for the same Keck data. Bottom "
                    "right: how well the orbit fits when its precession is forced to each value. With every position a "
                    "precession ~5× Einstein's fits best: a false signal from the 1990s speckle-imaging positions. With "
                    "adaptive-optics positions alone the data cannot tell. GRAVITY's data, which measured it, are not public.")
    return save(fig, "04_relativity")


# --------------------------------------------------------------------------
# 05  Scale
# --------------------------------------------------------------------------


def fig_scale(res: dict, art: dict) -> list[str]:
    p = best_params(art, "joint_gr")
    d = res["fits"]["joint_gr"]["derived"]
    a = d["semi_major_axis_au"]["median"]
    e = p["e"]
    rp = a * (1 - e)
    rs = d["schwarzschild_radius_au"]["median"]
    b_ax = a * np.sqrt(1 - e**2)
    th = np.linspace(0, 2 * np.pi, 1000)

    fig = plt.figure(figsize=(12.0, 10.6))
    ax_a = fig.add_axes([0.06, 0.575, 0.90, 0.35])
    ax_b = fig.add_axes([0.06, 0.07, 0.42, 0.42])
    ax_c = fig.add_axes([0.55, 0.07, 0.42, 0.42])

    def ellipse(ax, lw=1.1):
        ax.plot(a * np.cos(th) - a * e, b_ax * np.sin(th), color=ps.INK, lw=lw)  # focus at the origin

    def circle(ax, r, color, lw=0.9, ls="-"):
        ax.plot(r * np.cos(th), r * np.sin(th), color=color, lw=lw, ls=ls)

    # (a) the whole orbit, seen face-on.
    ellipse(ax_a)
    ax_a.scatter([0], [0], s=18, color=ps.INK, zorder=5)
    circle(ax_a, 169.0, ps.C_VLT, lw=0.9, ls=(0, (4, 2)))
    leader(ax_a, "Voyager 1's distance from the Sun\ntoday (169 AU), for scale", (-120, 119), (-520, 330))
    ax_a.text(-a * e, 0, f"S2's orbit seen face-on: {2 * a:,.0f} AU end to end, one lap every "
              f"{p['period']:.1f} years", ha="center", va="center", fontsize=9, color=ps.INK)
    leader(ax_a, f"apocentre: {a * (1 + e):,.0f} AU (2010)", (-a * (1 + e), 0), (-a * (1 + e) - 60, -440))
    leader(ax_a, f"pericentre: {rp:.0f} AU\n(2002, 2018)", (rp, 0), (260, -330))
    ax_a.set_xlim(-2.0 * a - 60, 460)
    ax_a.set_ylim(-1.12 * b_ax - 30, 1.12 * b_ax + 30)
    ax_a.set_aspect("equal", adjustable="box")
    ax_a.set_xlabel("AU from Sgr A*")
    ax_a.set_title("(a) the whole orbit, true shape", fontsize=10)

    # (b) zoom on the pericentre: the outer Solar System.
    ax = ax_b
    ellipse(ax, lw=1.3)
    ax.scatter([0], [0], s=16, color=ps.INK, zorder=5)
    circle(ax, 30.07, ps.C_KECK)
    circle(ax, 50.0, ps.GRID, lw=0.8)
    circle(ax, 169.0, ps.C_VLT, ls=(0, (4, 2)))
    ax.text(0, -38, "Neptune, 30 AU", fontsize=7.5, color=ps.INK_SECONDARY, ha="center", va="top")
    ax.text(0, 57, "Kuiper belt's edge, ~50 AU", fontsize=7.5, color=ps.INK_SECONDARY, ha="center")
    ax.text(-195, -190, "dashed: Voyager 1 today,\n169 AU from the Sun", fontsize=7.5, color=ps.INK_SECONDARY,
            va="bottom")
    ax.scatter([rp], [0], s=30, color=ps.C_GR, zorder=6, edgecolor=ps.SURFACE)
    leader(ax, f"S2 at pericentre:\n{rp:.0f} AU = {d['pericentre_light_hours']['median']:.0f} light-hours,\n"
               f"{d['pericentre_speed_kms']['median']:,.0f} km/s", (rp, 0), (124, 42), color=ps.INK, fontsize=7.5)
    ax.set_xlim(-200, 200)
    ax.set_ylim(-200, 200)
    ax.set_aspect("equal")
    ax.set_xlabel("AU")
    ax.set_title("(b) ×10: S2's closest approach against the Solar System", fontsize=10)

    # (c) zoom on the black hole.
    ax = ax_c
    ax.add_patch(patches.Circle((0, 0), rs, color=ps.INK, zorder=4))
    circle(ax, 1.5 * rs, ps.MUTED, lw=0.8, ls=(0, (3, 2)))
    circle(ax, 0.387, ps.C_KECK)
    circle(ax, 1.0, ps.C_KECK)
    ax.text(0, 0.41, "Mercury's orbit, 0.39 AU", fontsize=7.5, color=ps.INK_SECONDARY, ha="center")
    ax.text(0, 1.03, "Earth's orbit, 1 AU", fontsize=7.5, color=ps.INK_SECONDARY, ha="center")
    ax.add_patch(patches.Circle((0.70, -0.55), 0.00465, color=ps.YELLOW, zorder=4))
    leader(ax, "the Sun, to scale", (0.70, -0.55), (0.35, -0.85), fontsize=7.5)
    rs_km = rs * op.AU_KM
    leader(ax, f"event horizon: {2 * rs:.2f} AU across\n({rs_km / 1e6:.0f} million km radius, "
               f"{rs_km / 695_700:.0f}× the Sun's)", (rs * 0.7, -rs * 0.7), (-1.05, -0.62), color=ps.INK)
    leader(ax, "photon sphere: light\ncan orbit here", (-1.06 * rs, 1.06 * rs), (-0.95, 0.42), fontsize=7.5)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect("equal")
    ax.set_xlabel("AU")
    ax.set_title(f"(c) ×170 more: S2's closest approach is {d['pericentre_in_rs']['median']:,.0f} horizon radii",
                 fontsize=10)
    fig.suptitle(f"Four million suns inside {rp:.0f} AU: a mean density above "
                 f"{sci(d['density_bound_msun_pc3']['median'], 0)} suns per cubic parsec",
                 x=0.06, ha="left", y=0.975, fontsize=12)
    ps.caption(fig, "Every panel is to scale. Solar System distances are from the Sun, S2's from the black hole; the Solar "
                    "System is drawn only for scale (Sgr A* is 26,000 light-years away). Horizon: the Schwarzschild radius "
                    "2GM/c² for our mass (spin could shrink it by up to half). A cubic parsec around the Sun holds about "
                    "0.1 solar masses of stars.")
    return save(fig, "05_scale")


# --------------------------------------------------------------------------
# 06  The dark-cluster test
# --------------------------------------------------------------------------


def fig_dark_cluster(res: dict, art: dict) -> list[str]:
    case = res["black_hole_case"]
    mh = case["m_half_msun"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.4), gridspec_kw={"width_ratios": [1.25, 1]})

    maoz = bh.MAOZ_1997_GALAXY
    for m_half, rho, col, lw in ((maoz["m_half"], maoz["rho_half"], ps.MUTED, 1.0),
                                 (mh, case["rho_half_msun_pc3"], ps.INK, 1.7)):
        curves = bh.lifetime_curves(m_half, rho)
        for cls, c in curves.items():
            if cls != "primordial black holes":
                ax.plot(c["m"], c["lifetime"], color=col, lw=lw)
    now_curves = bh.lifetime_curves(mh, case["rho_half_msun_pc3"])
    c = now_curves["primordial black holes"]
    sel = c["m"] < 3.0
    ax.plot(c["m"][sel], c["lifetime"][sel], color=ps.C_PUBLISHED, lw=0.9, ls=(0, (4, 2)))
    ax.text(1.3e-5, 4e7, "violet dashes: hypothetical\nprimordial black holes", fontsize=7.5,
            color=ps.INK_SECONDARY)
    for cls, x, dy in (("brown dwarfs", 0.045, 0.22), ("white dwarfs", 0.3, 4.0), ("neutron stars", 1.45, 4.0),
                       ("stellar black holes", 9.0, 3.0)):
        cc = now_curves[cls]
        y = np.interp(np.log10(x), np.log10(cc["m"]), cc["lifetime"])
        ax.text(x, y * dy, cls, fontsize=7.5, color=ps.INK, ha="center")
    ax.text(3e-4, 3e1, "planets: < 10 years", fontsize=7.5, color=ps.INK, ha="center")
    ax.axhline(bh.GALAXY_AGE_YR, color=ps.C_GR, lw=1.2, ls="--")
    ax.text(1.3e-5, 1.5e10, "age of the Galaxy, 10 billion years", fontsize=8, color=ps.INK_SECONDARY)
    then, now = case["max_lifetime_1997"], case["max_lifetime_now"]
    leader(ax, f"best case in 1997: {sci(then['lifetime_yr'], 0)} years", (then["m_star"], then["lifetime_yr"]),
           (0.05, 1.2e9))
    leader(ax, f"best case now:\n{sci(now['lifetime_yr'], 0)} years\n({now['class'][:-1]}s)",
           (now["m_star"], now["lifetime_yr"]), (3.0, 1.5e3), color=ps.INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-5, 100)
    ax.set_ylim(1e1, 1e13)
    ax.set_xlabel("mass of each dark object (M☉)")
    ax.set_ylabel("lifetime of the cluster (years)")
    ax.set_title("Could 4 million suns of dark objects hide inside S2's orbit?", fontsize=10.5)
    ax.legend(handles=[Line2D([], [], color=ps.MUTED, lw=1.0,
                              label=f"1997: {sci(maoz['rho_half'], 0)} M☉/pc³ (Maoz 1998)"),
                       Line2D([], [], color=ps.INK, lw=1.7,
                              label=f"S2's pericentre: {sci(case['rho_half_msun_pc3'], 1)} M☉/pc³")],
              loc="upper right", fontsize=7.8, title="half-mass density of the cluster", title_fontsize=7.8)

    rho = np.logspace(6, 18, 90)
    tmax = np.array([bh.max_lifetime(mh, r)["lifetime_yr"] for r in rho])
    ax2.plot(rho, tmax, color=ps.INK, lw=1.5)
    ax2.axhline(bh.GALAXY_AGE_YR, color=ps.C_GR, lw=1.2, ls="--")
    need = case["density_needed_for_10_gyr"]
    ax2.axvline(need, color=ps.AXIS, lw=0.8)
    ax2.text(need * 1.6, 2e12, f"only below {sci(need, 0)} M☉/pc³\ncould a cluster last\n10 billion years",
             fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    for r, lab, col, off in ((maoz["rho_half"], "1997: the nearest\nstars' motions", ps.MUTED, (1 / 3e3, 1 / 60)),
                             (case["rho_half_msun_pc3"], "2019: S2's\npericentre", ps.INK, (1 / 2e3, 1 / 300))):
        tt = bh.max_lifetime(mh, r)["lifetime_yr"]
        ax2.scatter([r], [tt], s=40, color=col, zorder=5, edgecolor=ps.SURFACE)
        leader(ax2, lab, (r, tt), (r * off[0], tt * off[1]))
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(1e6, 1e18)
    ax2.set_ylim(1e3, 1e14)
    ax2.set_xlabel("half-mass density of the would-be cluster (M☉/pc³)")
    ax2.set_ylabel("longest possible lifetime (years)")
    ax2.set_title("The denser it must be, the faster it dies", fontsize=10.5)
    ps.caption(fig, "A cluster loses stars by evaporation (encounters kick members out: ~300 relaxation times) or destroys "
                    "itself in collisions, whichever is faster; every assumption favours its survival (Maoz 1998, whose "
                    "recipe and 1997 numbers we reproduce). The cluster holds half our mass inside S2's pericentre. Only "
                    f"exotic options survive: primordial black holes lighter than {sci(case['primordial_bh_mass_limit_msun'], 0)} "
                    "M☉ (about 20 Earth masses), or a star made of unknown particles (a boson star).")
    return save(fig, "06_dark_cluster")


# --------------------------------------------------------------------------
# 07  The shadow prediction
# --------------------------------------------------------------------------


def fig_shadow(res: dict, art: dict) -> list[str]:
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    rows = []
    for name, lab, col in (("keck_gr", "Keck alone (ours)", ps.C_KECK), ("vlt_gr", "VLT alone (ours)", ps.C_VLT),
                           ("joint_gr", "Keck + VLT (ours)", ps.C_JOINT)):
        s = res["fits"][name]["derived"]["shadow_diameter_uas"]
        rows.append((lab, s["median"], 0.5 * (s["hi68"] - s["lo68"]), col))
    for key, lab in (("do2019_gr", "Do+ 2019"), ("gravity2022", "GRAVITY 2022")):
        v = PUBLISHED[key]
        m, r = v["gm6"][0], v["r0"][0]
        dval = float(op.shadow_diameter_uas(m, r))
        rel = float(np.hypot(total_err(v["gm6"]) / m, total_err(v["r0"]) / r))
        rows.append((f"{lab} (from their M, R₀)", dval, dval * rel, ps.C_PUBLISHED))
    eht = PUBLISHED["eht2022"]
    ax.axvspan(eht["shadow_uas"][0] - eht["shadow_uas"][1], eht["shadow_uas"][0] + eht["shadow_uas"][1],
               color=ps.GRID, alpha=0.45, lw=0, zorder=0)
    ax.axvspan(eht["ring_uas"][0] - eht["ring_uas"][1], eht["ring_uas"][0] + eht["ring_uas"][1],
               color=ps.C_GR, alpha=0.16, lw=0, zorder=0)
    ax.axvline(eht["ring_uas"][0], color=ps.C_GR, lw=1.3, zorder=1)
    y = np.arange(len(rows))[::-1]
    for yy, (lab, v, err, col) in zip(y, rows):
        spin = 0.04 * v
        ax.plot([v - spin, v + spin], [yy, yy], color=col, lw=4, alpha=0.2, solid_capstyle="butt")
        ax.errorbar(v, yy, xerr=err, fmt="o", ms=5.5, color=col, ecolor=col, elinewidth=1.4, mec=ps.SURFACE)
        ax.text(36.4, yy, lab, va="center", fontsize=8, color=ps.INK_SECONDARY)
        ax.text(63.6, yy, f"{v:.1f} ± {err:.1f} µas", va="center", ha="right", fontsize=8, color=ps.INK)
    ax.text(eht["ring_uas"][0] + eht["ring_uas"][1] + 0.3, len(rows) - 0.4,
            f"EHT 2022, the ring:\n{eht['ring_uas'][0]} ± {eht['ring_uas'][1]} µas", fontsize=7.5,
            color=ps.INK_SECONDARY, va="top")
    ax.text(eht["shadow_uas"][0] - eht["shadow_uas"][1] + 0.3, -0.72,
            f"EHT 2022, the shadow itself: {eht['shadow_uas'][0]} ± {eht['shadow_uas'][1]} µas (grey)",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlim(36, 64)
    ax.set_ylim(-1.05, len(rows) - 0.25)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("angular diameter of Sgr A*'s shadow (micro-arcseconds)")
    ax.set_title("A prediction for the Event Horizon Telescope, from one star's orbit")
    ps.caption(fig, "Predicted diameter of the shadow of a non-spinning black hole, 2√27·GM/(c²R₀): it depends only on "
                    "M/R₀, which S2 pins down better than M or R₀ separately. Pale bars: ±4%, the most any spin or viewing "
                    "angle can change it. The EHT measured the bright ring around the shadow; converting ring to shadow "
                    "needs simulations of the glowing gas, hence the wider grey band.")
    return save(fig, "07_shadow_prediction")


# --------------------------------------------------------------------------
# 08  Validation against the teams' own fits
# --------------------------------------------------------------------------


def fig_validation(res: dict, art: dict) -> list[str]:
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.8), gridspec_kw={"width_ratios": [1.1, 1]})
    labels = {"gm6": "mass", "r0": "distance R₀", "f": "redshift factor f", "period": "period", "e": "eccentricity",
              "inc": "inclination", "Omega": "node Ω", "vz0": "velocity zero point", "nirc2_offset": "NIRC2 offset",
              "x0_keck": "frame x₀", "y0_keck": "frame y₀", "vx0_keck": "frame drift vx₀",
              "vy0_keck": "frame drift vy₀", "corr_mix_keck": "noise mixing p"}
    order = list(labels)
    #: The relativistic chains' period and omega are osculating elements at J2000 for a
    #: numerically integrated orbit; ours refer to the 2018 pericentre. Not comparable.
    skip = {"keck_gr": {"period", "omega"}, "keck_redshift": {"period", "omega"}, "keck_newton": set()}
    fits = (("keck_newton", "Newtonian model", "o", ps.C_KECK), ("keck_gr", "GR model", "s", ps.C_PUBLISHED),
            ("keck_redshift", "GR, f free", "D", ps.INK))
    for j, (name, lab, mk, col) in enumerate(fits):
        cmp = res["fits"][name]["vs_do2019_chain"]
        ys, xs = [], []
        for i, par in enumerate(order):
            if par in cmp and par not in skip[name] and cmp[par]["difference_in_sigma"] is not None:
                ys.append(len(order) - i + (1 - j) * 0.22)
                xs.append(cmp[par]["difference_in_sigma"])
        ax.scatter(xs, ys, marker=mk, s=24, color=col, edgecolor=ps.SURFACE, lw=0.5, label=lab, zorder=3)
    ax.axvspan(-1, 1, color=ps.GRID, alpha=0.5, lw=0, zorder=0)
    ax.axvline(0, color=ps.INK, lw=0.8)
    ax.set_yticks([len(order) - i for i in range(len(order))])
    ax.set_yticklabels([labels[k] for k in order], fontsize=8)
    ax.set_ylim(0.3, len(order) + 0.7)
    ax.set_xlim(-2.0, 2.0)
    ax.set_xlabel("(ours − their chain) / their posterior σ")
    ax.set_title("Keck data: our posteriors vs Do et al.'s own chains", fontsize=10.5)
    ax.legend(fontsize=7.5, loc="lower left", bbox_to_anchor=(0.0, -0.01))
    ax.grid(axis="y", visible=False)

    cmp = res["fits"]["vlt_kepler"]["vs_published"]
    items = [("gm6", "mass (10⁶ M☉)"), ("r0", "distance R₀ (kpc)"), ("alpha_mas", "frame α₀ at J2000 (mas)"),
             ("delta_mas", "frame δ₀ at J2000 (mas)"), ("valpha_uas_yr", "frame drift vα (µas/yr)"),
             ("vdelta_uas_yr", "frame drift vδ (µas/yr)"), ("vz_kms", "velocity zero point (km/s)")]
    for i, (k, lab) in enumerate(items):
        c = cmp[k]
        y = len(items) - i
        ax2.errorbar(0, y + 0.14, xerr=1, fmt="D", ms=4, color=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED, elinewidth=1.0)
        z = (c["ours"] - c["published"]) / c["published_err"]
        ze = c["ours_err"] / c["published_err"]
        ax2.errorbar(z, y - 0.14, xerr=ze, fmt="o", ms=4, color=ps.C_VLT, ecolor=ps.C_VLT, elinewidth=1.0,
                     mec=ps.SURFACE)
        ax2.text(2.2, y, f"{c['published']:g} ± {c['published_err']:g}\nvs {c['ours']:.3g} ± {c['ours_err']:.2g}",
                 fontsize=7.2, color=ps.INK_SECONDARY, va="center")
    ax2.set_yticks([len(items) - i for i in range(len(items))])
    ax2.set_yticklabels([lab for _, lab in items], fontsize=8)
    ax2.axvspan(-1, 1, color=ps.GRID, alpha=0.5, lw=0, zorder=0)
    ax2.axvline(0, color=ps.INK, lw=0.8)
    ax2.set_xlim(-2.0, 3.4)
    ax2.set_ylim(-0.4, len(items) + 0.6)
    ax2.set_xticks([-2, -1, 0, 1, 2])
    ax2.set_xlabel("(ours − published) / published σ")
    ax2.set_title("VLT data: our fit vs Gillessen et al. 2017, Table 1", fontsize=10.5)
    ax2.legend(handles=[Line2D([], [], marker="D", ls="", color=ps.C_PUBLISHED, ms=4, label="published (S2, VLT, no priors)"),
                        Line2D([], [], marker="o", ls="", color=ps.C_VLT, ms=4, label="ours: same data, same model")],
               fontsize=7.5, loc="lower left")
    ax2.grid(axis="y", visible=False)
    ps.caption(fig, "Before trusting anything new, we re-ran each team's own analysis on its own data. Grey band: within "
                    "1σ. Not compared: the period and pericentre angle of the relativistic fits, which Do et al. define as "
                    "osculating elements at J2000 while ours refer to the 2018 pericentre. Do et al. sampled with nested "
                    "sampling (MultiNest); we use an ensemble MCMC (emcee).")
    return save(fig, "08_validation")


# --------------------------------------------------------------------------


def make_all(res: dict, art: dict) -> list[str]:
    written = []
    for fn in (fig_sky, fig_rv, fig_mass_distance, fig_redshift, fig_scale, fig_dark_cluster, fig_shadow,
               fig_validation):
        written += fn(res, art)
    return written
