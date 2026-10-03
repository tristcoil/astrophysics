"""Static figures (PNG + SVG), built from the analysis results and the cached data.

Nothing here re-measures anything: every number printed on a figure comes from
``analysis.run_all()``, so the figures and the README cannot drift apart.

Design notes. Each kind of object keeps one colour throughout (see
``plotstyle``): spacecraft blue, interstellar objects orange, comets aqua,
planets neutral ink, the Sun drawn as its astronomical symbol. Individual
objects are always named by a direct label, so no figure relies on colour
alone. Top-down maps look down on the plane of Earth's orbit from the north,
with the x-axis pointing to the March-equinox direction; bodies far above or
below that plane appear closer to the Sun than they are, so their labels carry
the true distance.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import horizons as hz
import orbit_physics as op
import plotstyle as ps
from config import FIGURE_DIR, PLANETS, PUBLISHED, SNAPSHOT_DATE, body

ps.use_style()

PROVENANCE = "Data: JPL Horizons ephemeris service (Solar System Dynamics, JPL/Caltech)."

YEAR = 365.25
ORBIT_LW = 0.8


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _jd(date: str) -> float:
    from astropy.time import Time

    return float(Time(date, scale="tdb").jd)


def _year(jd):
    """Decimal year from a Julian date (TDB) -- close enough for axes."""
    return 2000.0 + (np.asarray(jd) - 2451545.0) / YEAR


SNAP_JD = _jd(SNAPSHOT_DATE)


def _sun(ax, size=7.0, zorder=8):
    """The Sun as the astronomical symbol: a ring with a central dot."""
    ax.plot(0, 0, marker="o", ms=size, mfc=ps.SURFACE, mec=ps.C_SUN, mew=1.1, zorder=zorder)
    ax.plot(0, 0, marker="o", ms=size * 0.3, color=ps.C_SUN, zorder=zorder + 1)


def _one_orbit(key: str, around_jd: float = SNAP_JD) -> np.ndarray:
    """Positions over one orbital period centred on a date, closed into a loop."""
    d = hz.load(hz.orbit_key(body(key)))
    period = body(key).period_days
    m = np.abs(d["jd"] - around_jd) <= period / 2 + 1e-9
    r = d["r"][m]
    return np.vstack([r, r[:1]])


def _label(ax, x, y, text, *, dist=7, color=ps.INK_SECONDARY, size=7.5, angle=None,
           weight="normal", **kw):
    """Label a point, nudged radially away from the Sun so labels fan outwards."""
    ang = np.arctan2(y, x) if angle is None else np.radians(angle)
    dx, dy = dist * np.cos(ang), dist * np.sin(ang)
    ha = "left" if dx > 1 else ("right" if dx < -1 else "center")
    va = "bottom" if dy > 1 else ("top" if dy < -1 else "center")
    return ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                       ha=ha, va=va, fontsize=size, color=color, weight=weight, **kw)


def _square(ax, half):
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")


def _map_axes(ax, half, title):
    _square(ax, half)
    ax.set_title(title)
    ax.set_xlabel("x (au)  →  March-equinox direction")
    ax.set_ylabel("y (au)")


# --------------------------------------------------------------------------
# 1. The solar system today
# --------------------------------------------------------------------------

def _callout(ax, xy, text_xy, text, *, ha="left", va="center", color=ps.INK, size=7.2):
    """A label placed in free space, tied to its point by a hairline leader."""
    return ax.annotate(text, xy, xytext=text_xy, textcoords="data", ha=ha, va=va,
                       fontsize=size, color=color,
                       arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.6,
                                   "shrinkA": 1, "shrinkB": 3})


def figure_today(res):
    snap = res["snapshot"]["bodies"]
    sc = res["spacecraft"]["craft"]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4), gridspec_kw={"wspace": 0.3})
    fig.subplots_adjust(top=0.9)
    ax1, ax2, ax3 = axes

    def xy(key):
        return snap[key]["r"][0], snap[key]["r"][1]

    # (a) Inner system.
    for key in ("mercury", "venus", "earth", "mars", "ceres"):
        orb = _one_orbit(key)
        ax1.plot(orb[:, 0], orb[:, 1], color=ps.C_ORBIT, lw=ORBIT_LW, zorder=1)
        filled = key != "ceres"
        ax1.plot(*xy(key), "o", ms=6, mfc=ps.C_PLANET if filled else ps.SURFACE,
                 mec=ps.C_PLANET, mew=1.1, zorder=5)
        _label(ax1, *xy(key), snap[key]["name"] + ("" if filled else " (dwarf planet)"),
               color=ps.INK)
    _sun(ax1)
    _map_axes(ax1, 3.2, "Inner solar system")
    ax1.text(0.02, 0.02, "planets move anticlockwise in this view", transform=ax1.transAxes,
             fontsize=7, color=ps.MUTED)

    # (b) Planets and visitors.
    for key in ("jupiter", "saturn", "uranus", "neptune", "pluto"):
        orb = _one_orbit(key)
        ax2.plot(orb[:, 0], orb[:, 1], color=ps.C_ORBIT, lw=ORBIT_LW, zorder=1)
        filled = key != "pluto"
        ax2.plot(*xy(key), "o", ms=5.5, mfc=ps.C_PLANET if filled else ps.SURFACE,
                 mec=ps.C_PLANET, mew=1.1, zorder=5)
    for key, (dx, dy, ha, va) in {"jupiter": (-2.5, -1.5, "right", "top"),
                                  "saturn": (2.5, -1.5, "left", "top"),
                                  "uranus": (2.5, 1.5, "left", "bottom"),
                                  "neptune": (2.5, -1.5, "left", "top"),
                                  "pluto": (2.5, -1.5, "left", "top")}.items():
        x, y = xy(key)
        ax2.text(x + dx, y + dy, snap[key]["name"], ha=ha, va=va, fontsize=7.5, color=ps.INK)

    halley = hz.load(hz.track_key(body("halley")))
    m = (halley["jd"] >= _jd("1986-02-09")) & (halley["jd"] <= _jd("2061-07-29"))
    ax2.plot(halley["r"][m, 0], halley["r"][m, 1], color=ps.C_COMET, lw=1.0, zorder=2,
             label="comet (Halley: whole orbit)")
    hb = hz.load(hz.track_key(body("hale_bopp")))
    m = hb["jd"] <= SNAP_JD
    ax2.plot(hb["r"][m, 0], hb["r"][m, 1], color=ps.C_COMET, lw=1.0, zorder=2)
    for key in ("oumuamua", "borisov", "atlas"):
        tr = hz.load(hz.track_key(body(key)))
        m = (tr["jd"] >= _jd("2010-01-01")) & (tr["jd"] <= SNAP_JD)
        ax2.plot(tr["r"][m, 0], tr["r"][m, 1], color=ps.C_INTERSTELLAR, lw=1.1, zorder=3,
                 label="interstellar object (path since 2010)" if key == "atlas" else None)
    for key, color in (("halley", ps.C_COMET), ("hale_bopp", ps.C_COMET),
                       ("oumuamua", ps.C_INTERSTELLAR), ("borisov", ps.C_INTERSTELLAR),
                       ("atlas", ps.C_INTERSTELLAR), ("new_horizons", ps.C_SPACECRAFT)):
        ax2.plot(*xy(key), "o", ms=5.5, color=color, mec=ps.SURFACE, mew=1.0, zorder=6)
    s_ = snap
    _callout(ax2, xy("halley"), (-30, 42),
             f"1P/Halley · {s_['halley']['dist_au']:.0f} au\npast aphelion, returning",
             ha="center", va="bottom")
    _callout(ax2, xy("atlas"), (-24, 14),
             f"3I/ATLAS · {s_['atlas']['dist_au']:.0f} au · {s_['atlas']['speed_kms']:.0f} km/s",
             ha="right")
    _callout(ax2, xy("oumuamua"), (68, 27),
             f"1I/ʻOumuamua · {s_['oumuamua']['dist_au']:.0f} au · "
             f"{s_['oumuamua']['speed_kms']:.0f} km/s", ha="right", va="bottom")
    _callout(ax2, xy("borisov"), (-6, -48),
             f"2I/Borisov · {s_['borisov']['dist_au']:.0f} au · "
             f"{s_['borisov']['speed_kms']:.0f} km/s", ha="right", va="top")
    _callout(ax2, xy("hale_bopp"), (-4, -24),
             f"Hale-Bopp · {s_['hale_bopp']['dist_au']:.0f} au\n"
             f"({abs(s_['hale_bopp']['ecl_lat_deg']):.0f}° below the plane)", ha="right",
             va="top")
    _callout(ax2, xy("new_horizons"), (26, -58),
             f"New Horizons · {s_['new_horizons']['dist_au']:.0f} au", ha="left", va="bottom")
    _sun(ax2, size=5)
    _map_axes(ax2, 72, "Planets, comets and interstellar visitors")
    ax2.legend(loc="lower left", fontsize=7)

    # (c) The escaping spacecraft.
    for key in ("neptune", "pluto"):
        orb = _one_orbit(key)
        ax3.plot(orb[:, 0], orb[:, 1], color=ps.C_ORBIT, lw=ORBIT_LW, zorder=1)
    helio = PUBLISHED["heliopause"]
    for b in ("voyager1", "voyager2", "pioneer10", "pioneer11", "new_horizons"):
        tr = hz.load(hz.track_key(body(b)))
        past = tr["jd"] <= SNAP_JD
        ax3.plot(tr["r"][past, 0], tr["r"][past, 1], color=ps.C_SPACECRAFT, lw=1.2, zorder=3,
                 label="since launch" if b == "voyager1" else None)
        fut = (tr["jd"] >= SNAP_JD) & (tr["jd"] <= _jd("2050-01-01"))
        ax3.plot(tr["r"][fut, 0], tr["r"][fut, 1], color=ps.C_SPACECRAFT, lw=1.0,
                 ls=(0, (2, 2)), alpha=0.8, zorder=3,
                 label="predicted to 2050" if b == "voyager1" else None)
        ax3.plot(*xy(b), "o", ms=6, color=ps.C_SPACECRAFT, mec=ps.SURFACE, mew=1.0, zorder=6)
        if b in helio:
            p = hz.load(hz.track_key(body(b)))
            i = int(np.argmin(np.abs(p["jd"] - _jd(helio[b]["date"]))))
            ax3.plot(p["r"][i, 0], p["r"][i, 1], marker="o", ms=8, mfc="none",
                     mec=ps.INK_SECONDARY, mew=1.0, ls="none", zorder=7,
                     label="left the heliosphere" if b == "voyager1" else None)

    def sc_text(b, where):
        c = sc[b]
        return (f"{c['name']}  {c['today_r_au']:.1f} au · {c['today_light_hours']:.1f} light-h"
                f"\n{abs(c['today_ecliptic_lat_deg']):.0f}° {where} the plane")

    _callout(ax3, xy("voyager1"), (-45, -178), sc_text("voyager1", "above"), ha="right")
    _callout(ax3, xy("voyager2"), (62, -100), sc_text("voyager2", "below"))
    _callout(ax3, xy("pioneer11"), (62, -150), sc_text("pioneer11", "above"))
    _callout(ax3, xy("new_horizons"), (62, -45), sc_text("new_horizons", "above"))
    _callout(ax3, xy("pioneer10"), (40, 150), sc_text("pioneer10", "above"))
    _callout(ax3, (-21.2, 21.2), (-70, 62), "orbits of Neptune\nand Pluto", ha="right",
             color=ps.MUTED, size=6.5)
    _sun(ax3, size=4)
    _map_axes(ax3, 205, "The five spacecraft leaving the solar system")
    ax3.legend(loc="center left", fontsize=7)

    fig.suptitle(f"The solar system on {SNAPSHOT_DATE} — three scales", x=0.125, ha="left",
                 y=0.995)
    ps.caption(fig, f"{PROVENANCE} Positions at 00:00 TDB, seen from above the plane of "
                    "Earth's orbit (the ecliptic). The Voyagers, Hale-Bopp and 2I/Borisov are far "
                    "above or below that plane, so they appear closer to the Sun than they are; "
                    "labels give true distances. Light-hours: one-way signal time from Earth.")
    return ps.save(fig, FIGURE_DIR / "01_solar_system_today")


# --------------------------------------------------------------------------
# 2. Kepler's third law
# --------------------------------------------------------------------------

def figure_kepler3(res):
    k = res["kepler3"]
    rows = k["bodies"]
    fig = plt.figure(figsize=(13.5, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1], hspace=0.55, wspace=0.28)
    ax = fig.add_subplot(gs[:, 0])
    ax_in = fig.add_subplot(gs[0, 1])
    ax_out = fig.add_subplot(gs[1, 1])

    # (a) The law.
    a_line = np.logspace(np.log10(0.3), np.log10(60), 200)
    p_line = op.kepler_period(a_line, k["gm_sun_measured_au3_d2"]) / YEAR
    ax.plot(a_line, p_line, color=ps.C_MODEL, lw=1.4, zorder=1,
            label="P = 2π √(a³ / GM☉), with our measured GM☉")
    for key, r in rows.items():
        dwarf = key in ("pluto", "ceres")
        ax.plot(r["a_au"], r["period_days"] / YEAR, "o", ms=7,
                mfc=ps.SURFACE if dwarf else ps.C_PLANET, mec=ps.C_PLANET, mew=1.2, zorder=4)
        _label(ax, r["a_au"], r["period_days"] / YEAR,
               f"{r['name']}  ({r['n_orbits']} orbit{'s' if r['n_orbits'] > 1 else ''})",
               angle={"pluto": -35, "neptune": 180}.get(key, 150), color=ps.INK, size=7.5,
               dist=8)
    h = k["halley"]
    ax.plot(h["a_au"], h["period_days"] / YEAR, "o", ms=7, color=ps.C_COMET, mec=ps.SURFACE,
            zorder=5)
    _label(ax, h["a_au"], h["period_days"] / YEAR, "Halley (1986–2061)",
           angle=-30, color=ps.INK, size=7.5, dist=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([0.3, 1, 3, 10, 30], ["0.3", "1", "3", "10", "30"])
    ax.set_yticks([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
    ax.minorticks_off()
    ax.set_xlabel("orbit size: semi-major axis a (au)")
    ax.set_ylabel("orbital period P (years)")
    ax.set_title("Kepler's third law, measured")
    ax.legend(loc="lower right")
    txt = (f"fitted slope of log P vs log a: {k['loglog_slope']:.4f}   (Kepler: 3/2)\n"
           f"the Sun's GM from the inner planets: {k['gm_sun_measured_si']:.7e} m³/s²\n"
           f"   = JPL's value {k['gm_sun_vs_de440_ppm']:+.1f} ± "
           f"{k['gm_sun_measured_err_ppm']:.1f} parts per million\n"
           f"the Sun's mass (÷ G): {k['sun_mass_kg']:.4e} kg")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=7.8, color=ps.INK,
            linespacing=1.5)

    # (b), (c) The residuals: how far each orbit's GM sits above the Sun's.
    def residual_panel(axr, keys, title, ylim):
        x = np.arange(len(keys))
        for i, key in enumerate(keys):
            r = rows[key]
            err = r["excess_err_ppm"]
            axr.errorbar(i, r["excess_ppm"], yerr=err if err else None, fmt="o", ms=6.5,
                         color=ps.C_MEASURED, mec=ps.SURFACE, mew=1.0, capsize=0, lw=1.2,
                         zorder=5, label="measured from the orbit" if i == 0 else None)
            axr.plot(i + 0.22, r["own_mass_ppm"], marker="o", ms=6, mfc=ps.SURFACE,
                     mec=ps.C_CONTEXT, mew=1.1, ls="none", zorder=3,
                     label="Sun + the planet's own mass" if i == 0 else None)
            axr.plot(i - 0.22, r["predicted_ppm"], marker="D", ms=5, color=ps.C_MODEL,
                     mec=ps.SURFACE, mew=0.8, ls="none", zorder=4,
                     label="+ every other planet, smeared into a ring" if i == 0 else None)
        axr.axhline(0, color=ps.AXIS, lw=0.8, zorder=1)
        axr.set_xticks(x, [rows[k_]["name"] for k_ in keys])
        axr.set_xlim(-0.6, len(keys) - 0.4)
        axr.set_ylim(*ylim)
        axr.set_ylabel("GM from orbit ÷ GM☉ − 1 (ppm)")
        axr.set_title(title)
        axr.grid(axis="x", visible=False)

    residual_panel(ax_in, ["mercury", "venus", "earth", "mars"],
                   "Inner planets: the Sun alone, to a few parts per million", (-9, 7))
    handles, labels = ax_in.get_legend_handles_labels()
    order = [labels.index(x) for x in ("measured from the orbit",
                                       "+ every other planet, smeared into a ring",
                                       "Sun + the planet's own mass")]
    ax_in.legend([handles[i] for i in order], [labels[i] for i in order], loc="lower left",
                 fontsize=7)
    residual_panel(ax_out, ["ceres", "jupiter", "saturn", "uranus", "neptune", "pluto"],
                   "Beyond Mars: the other planets show up", (-250, 1750))
    j = rows["jupiter"]
    ratio = k["jupiter_sun_over_jupiter_measured"]
    ratio_err = k["jupiter_mass_ratio_err"] / k["jupiter_mass_ratio_measured"] * ratio
    ps.annotate(ax_out, f"Jupiter's own mass, from its orbit:\nthe Sun is {ratio:.0f} ± "
                        f"{ratio_err:.0f} Jupiters\n(JPL: {k['jupiter_sun_over_jupiter_jpl']:.1f})",
                (1, j["excess_ppm"]), (1.35, 330), color=ps.INK_SECONDARY)
    ax_out.text(4.0, 1060, "Uranus, Neptune and Pluto orbit\nthe Sun plus Jupiter and Saturn",
                fontsize=7, color=ps.INK_SECONDARY, ha="center", va="center")
    ps.annotate(ax_out, "Jupiter, as a ring,\npulls Ceres outward", (0, rows["ceres"]["excess_ppm"]),
                (0.1, 560), color=ps.INK_SECONDARY, ha="center")

    ps.caption(fig, f"{PROVENANCE} Heliocentric positions of each planet system's barycentre, "
                    "1950–2050 (inner planets) and 1750–2250 (outer). Period: time to return to "
                    "the same direction, averaged over every complete orbit. Size: an ellipse with "
                    "the Sun at its focus fitted to every position (Halley: half the sum of its "
                    "perihelion and aphelion distances). Error bars: orbit-to-orbit "
                    "scatter ÷ √(number of orbits); none for Pluto (2 orbits). Planet masses in "
                    "the models are JPL's (DE441).")
    return ps.save(fig, FIGURE_DIR / "02_kepler_third_law")


# --------------------------------------------------------------------------
# 3. Kepler's second law: Halley's comet
# --------------------------------------------------------------------------

def figure_halley(res, art):
    k2 = res["kepler2"]
    a = art["kepler2"]
    jd, r, v = a["jd"], a["r"], a["v"]
    t_p1, t_p2 = a["t_peri"]
    t_a = a["t_aph"]

    # In-plane coordinates: x towards the 1986 perihelion, y 90 deg ahead.
    i_p = int(np.argmin(np.abs(jd - t_p1)))
    P = op.unit(op.eccentricity_vector(r[i_p], v[i_p]))
    Q = np.cross(op.unit(np.cross(r[i_p], v[i_p])), P)
    x, y = op.dot(r, P), op.dot(r, Q)
    orbit = (jd >= t_p1) & (jd <= t_p2)
    half = YEAR / 2

    def wedge_pts(t_mid):
        m = (jd >= t_mid - half) & (jd <= t_mid + half)
        return np.column_stack([np.r_[0, x[m], 0], np.r_[0, y[m], 0]])

    fig = plt.figure(figsize=(13.5, 8.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], width_ratios=[0.8, 1, 1],
                          hspace=0.45, wspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    ax_z = fig.add_subplot(gs[1, 0])
    ax_v = fig.add_subplot(gs[1, 1])
    ax_h = fig.add_subplot(gs[1, 2])

    # (a) The whole orbit, to scale, in its own plane.
    for name, dist in (("Jupiter", 5.20), ("Saturn", 9.54), ("Uranus", 19.19),
                       ("Neptune", 30.07)):
        th = np.linspace(0, 2 * np.pi, 400)
        ax.plot(dist * np.cos(th), dist * np.sin(th), color=ps.GRID, lw=0.8, zorder=0)
        if dist > 6:
            ax.text(-np.sqrt(dist**2 - 5.9**2), 5.95, f"{name}'s distance", fontsize=6.3,
                    color=ps.MUTED, ha="center", va="top")
    ax.plot(x[orbit], y[orbit], color=ps.C_COMET, lw=1.3, zorder=3)
    ax.add_patch(Polygon(wedge_pts(t_p1), closed=True, facecolor=ps.C_COMET, alpha=0.22,
                         edgecolor="none", zorder=2))
    ax.add_patch(Polygon(wedge_pts(t_a), closed=True, facecolor=ps.C_COMET, alpha=0.75,
                         edgecolor="none", zorder=2))
    yrs = np.arange(np.ceil(_year(t_p1)), _year(t_p2))
    tj = np.array([_jd(f"{int(yy)}-01-01") for yy in yrs])
    ax.plot(np.interp(tj, jd, x), np.interp(tj, jd, y), "o", ms=3.2, color=ps.C_COMET,
            mec=ps.SURFACE, mew=0.5, zorder=4)
    for yy, t in zip(yrs, tj):
        if int(yy) in (1990, 2000, 2010, 2040, 2050, 2060):
            xx, yv = np.interp(t, jd, x), np.interp(t, jd, y)
            ax.annotate(str(int(yy)), (xx, yv), xytext=(0, 7 if yv > 0 else -7),
                        textcoords="offset points", ha="center",
                        va="bottom" if yv > 0 else "top", fontsize=7, color=ps.INK_SECONDARY)
    xs, ys = np.interp(SNAP_JD, jd, x), np.interp(SNAP_JD, jd, y)
    ax.plot(xs, ys, "o", ms=7, mfc=ps.SURFACE, mec=ps.INK, mew=1.3, zorder=6)
    ax.annotate(f"today: {k2['today_r_au']:.1f} au,\n{k2['today_speed_kms']:.2f} km/s, "
                "falling back in", (xs, ys), xytext=(-33.6, -0.5), textcoords="data",
                fontsize=7.3, color=ps.INK, va="top",
                arrowprops={"arrowstyle": "-", "color": ps.INK, "lw": 0.6})
    _sun(ax, size=6)
    pa, aa = k2["year_at_perihelion"], k2["year_at_aphelion"]
    ps.annotate(ax, f"the year around perihelion (1986)\nsweeps {pa['angle_deg']:.0f}° — "
                    f"area {pa['area_au2']:.3f} au²", (-0.9, 1.6), (-12.5, 2.2), color=ps.INK)
    ps.annotate(ax, f"the year around aphelion (Dec 2023)\nsweeps {aa['angle_deg']:.2f}° — "
                    f"area {aa['area_au2']:.3f} au²", (-34.6, 0.12), (-27.5, 2.2), color=ps.INK)
    ax.text(-15.0, -0.35, "↑ that year is this sliver, drawn to scale: 35 au long, "
                          f"{aa['angle_deg'] / 180 * np.pi * k2['Q_au']:.2f} au wide",
            fontsize=7, color=ps.INK_SECONDARY, ha="center", va="top")
    ax.set_xlim(-37.2, 3.0)
    ax.set_ylim(-6.0, 6.1)
    ax.set_aspect("equal")
    ax.set_xlabel("au, in the plane of Halley's orbit (perihelion to the right)")
    ax.set_ylabel("au")
    ax.set_title("Kepler's second law: equal areas in equal times — the two shaded years "
                 f"differ by {abs(k2['equal_area_ratio'] - 1) * 100:.1f}%")

    # (b) Perihelion up close.
    ax_z.add_patch(Polygon(wedge_pts(t_p1), closed=True, facecolor=ps.C_COMET, alpha=0.22,
                           edgecolor="none"))
    mm = (jd >= t_p1 - 1.3 * YEAR) & (jd <= t_p1 + 1.3 * YEAR)
    ax_z.plot(x[mm], y[mm], color=ps.C_COMET, lw=1.3)
    months = t_p1 + np.arange(-6, 7) * YEAR / 12
    ax_z.plot(np.interp(months, jd, x), np.interp(months, jd, y), "o", ms=3.5,
              color=ps.C_COMET, mec=ps.SURFACE, mew=0.5)
    th = np.linspace(0, 2 * np.pi, 200)
    ax_z.plot(np.cos(th), np.sin(th), color=ps.GRID, lw=0.8, zorder=0)
    ax_z.text(0.74, -0.8, "1 au", fontsize=6.5, color=ps.MUTED, ha="left", va="top")
    ax_z.annotate(f"perihelion {k2['perihelion_1_utc'][:10]}\n{k2['q1_au']:.3f} au, "
                  f"{k2['speed_max_kms']:.1f} km/s", (k2["q1_au"], 0), xytext=(-8, 26),
                  textcoords="offset points", fontsize=7, color=ps.INK, ha="right",
                  arrowprops={"arrowstyle": "-", "color": ps.INK, "lw": 0.6})
    _sun(ax_z, size=5)
    ax_z.set_xlim(-3.0, 1.4)
    ax_z.set_ylim(-2.7, 3.3)
    ax_z.set_aspect("equal")
    ax_z.set_title("The perihelion year up close")
    ax_z.set_xlabel("au")
    ax_z.text(0.03, 0.03, "dots one month apart", transform=ax_z.transAxes, fontsize=7,
              color=ps.MUTED)

    # (c) Speed.
    yr = _year(jd[orbit])
    ax_v.plot(yr, a["speed_kms"][orbit], color=ps.C_COMET, lw=1.3)
    ax_v.set_yscale("log")
    ax_v.set_yticks([1, 3, 10, 30], ["1", "3", "10", "30"])
    ax_v.minorticks_off()
    ax_v.set_ylim(0.6, 90)
    ax_v.set_ylabel("speed around the Sun (km/s)")
    ax_v.set_title(f"Speed changes {k2['speed_ratio']:.0f}-fold")
    ax_v.annotate(f"{k2['speed_max_kms']:.1f} km/s", (yr[0], k2["speed_max_kms"]),
                  xytext=(6, -2), textcoords="offset points", fontsize=7, va="top",
                  color=ps.INK_SECONDARY)
    ax_v.annotate(f"{k2['speed_min_kms']:.2f} km/s at\naphelion, {k2['aphelion_utc'][:4]}",
                  (_year(t_a), k2["speed_min_kms"]), xytext=(-12, -2), textcoords="offset points",
                  ha="right", va="top", fontsize=7, color=ps.INK_SECONDARY)
    ax_v.set_xlim(1984, 2064)

    # (d) Areal velocity.
    dev = (a["areal"][orbit] / a["mean_areal"] - 1) * 100
    ax_h.plot(yr, dev, color=ps.C_COMET, lw=1.2)
    ax_h.axhline(0, color=ps.AXIS, lw=0.8)
    ax_h.set_ylim(-3.5, 3.5)
    ax_h.set_xlim(1984, 2064)
    ax_h.set_ylabel("area swept per day, % from average")
    ax_h.set_title(f"Area swept per day: steady to ±{k2['areal_velocity_max_dev_percent']:.1f}%")
    ax_h.text(0.5, 0.97,
              f"The ripple repeats every {k2['ripple_period_years']:.1f} years: Jupiter's year.\n"
              f"Jupiter swings the Sun itself around at {k2['sun_reflex_speed_ms']:.0f} m/s;\n"
              f"against a comet crawling at 0.9 km/s, that predicts ±{k2['expected_ripple_percent']:.1f}%.",
              transform=ax_h.transAxes, ha="center", va="top", fontsize=7, color=ps.INK_SECONDARY)
    for axx in (ax_v, ax_h):
        axx.axvline(_year(SNAP_JD), color=ps.INK_SECONDARY, lw=0.8, zorder=0)
        axx.set_xlabel("year")
    ax_v.text(_year(SNAP_JD) + 1, 50, "today", fontsize=7, color=ps.INK_SECONDARY)

    ps.caption(fig, f"{PROVENANCE} Comet orbit solution {k2['source'].get('target_source')}, "
                    f"from {k2['source'].get('n_obs'):,} observations {k2['source'].get('obs_span')}. "
                    "Heliocentric positions every day from perihelion in 1986 to perihelion in 2061; "
                    "dots in the top panel mark each 1 January. Rings show planets' distances "
                    "only — Halley's orbit is tilted 162° to theirs, so it runs the other way round.")
    return ps.save(fig, FIGURE_DIR / "03_halley_equal_areas")


# --------------------------------------------------------------------------
# 4. Interstellar visitors
# --------------------------------------------------------------------------

def figure_interstellar(res, art):
    iso = res["interstellar"]["objects"]
    ia = art["interstellar"]
    fig, (ax, axv) = plt.subplots(1, 2, figsize=(14, 6.2),
                                  gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.25})
    fig.subplots_adjust(top=0.87)

    for key in ("mercury", "venus", "earth", "mars", "jupiter"):
        orb = _one_orbit(key)
        ax.plot(orb[:, 0], orb[:, 1], color=ps.C_ORBIT, lw=ORBIT_LW, zorder=1)
    for name, dist in (("Earth", 1.0), ("Mars", 1.52), ("Jupiter", 5.2)):
        ax.text(0.0, -dist - 0.06, name, fontsize=6.5, color=ps.MUTED, ha="center", va="top")
    for key, o in iso.items():
        tr = ia[key]["passage"]
        tp = ia[key]["tp"]
        long = ia[key]["track"]
        m = np.abs(long["jd"] - tp) <= 4 * YEAR
        ax.plot(long["r"][m, 0], long["r"][m, 1], color=ps.C_INTERSTELLAR, lw=1.4, zorder=3)
        # Arrow showing direction of travel, just after perihelion.
        i0 = int(np.argmin(np.abs(tr["jd"] - (tp + 25))))
        p0, p1 = tr["r"][i0], tr["r"][i0 + 12]
        ax.annotate("", xy=(p1[0], p1[1]), xytext=(p0[0], p0[1]),
                    arrowprops={"arrowstyle": "-|>", "color": ps.C_INTERSTELLAR, "lw": 1.2,
                                "mutation_scale": 11}, zorder=4)
        ip = int(np.argmin(np.abs(tr["jd"] - tp)))
        ax.plot(tr["r"][ip, 0], tr["r"][ip, 1], "o", ms=5, color=ps.C_INTERSTELLAR,
                mec=ps.SURFACE, mew=0.8, zorder=5)
        spot, ha = {"oumuamua": ((2.2, 2.6), "left"), "borisov": ((-2.4, 2.9), "right"),
                    "atlas": ((-2.6, -1.9), "right")}[key]
        _callout(ax, (tr["r"][ip, 0], tr["r"][ip, 1]), spot,
                 f"{o['name']}\nperihelion {o['perihelion_utc'][:10]}\n{o['q_au']:.2f} au "
                 f"from the Sun, {o['speed_at_perihelion_kms']:.0f} km/s", ha=ha, size=7.2)
    _sun(ax, size=6)
    _map_axes(ax, 6.2, "Three paths through the inner solar system")
    ax.text(0.02, 0.02, "each path shown ±4 years from perihelion; arrows show direction",
            transform=ax.transAxes, fontsize=7, color=ps.MUTED)

    # Speed against distance: bound or not, at a glance.
    rr = np.logspace(np.log10(0.2), np.log10(120), 300)
    axv.plot(rr, op.to_kms(op.escape_speed(rr)), color=ps.INK, lw=1.1, ls=(0, (5, 3)), zorder=2)
    for key, o in iso.items():
        tr = ia[key]["track"]
        pas = ia[key]["passage"]
        jd = np.r_[tr["jd"][tr["jd"] < pas["jd"][0]], pas["jd"], tr["jd"][tr["jd"] > pas["jd"][-1]]]
        rv = np.vstack([tr["r"][tr["jd"] < pas["jd"][0]], pas["r"],
                        tr["r"][tr["jd"] > pas["jd"][-1]]])
        vv = np.vstack([tr["v"][tr["jd"] < pas["jd"][0]], pas["v"],
                        tr["v"][tr["jd"] > pas["jd"][-1]]])
        dist, spd = op.norm(rv), op.to_kms(op.norm(vv))
        inb = jd <= ia[key]["tp"]
        axv.plot(dist[inb], spd[inb], color=ps.C_INTERSTELLAR, lw=1.5, zorder=4)
        axv.plot(dist[~inb], spd[~inb], color=ps.C_INTERSTELLAR, lw=1.0, alpha=0.45, zorder=3)
        below = key == "oumuamua"
        axv.text(118, o["v_inf_in_kms"] * (0.96 if below else 1.04),
                 f"{o['name']}: {o['v_inf_in_kms']:.2f} km/s at infinity",
                 ha="right", va="top" if below else "bottom", fontsize=7.3, color=ps.INK)
    for key, lab in (("halley", "Halley"), ("hale_bopp", "Hale-Bopp")):
        tr = hz.load(hz.track_key(body(key)))
        axv.plot(op.norm(tr["r"]), op.to_kms(op.norm(tr["v"])), color=ps.C_COMET, lw=1.2,
                 zorder=3)
    tr = hz.load(hz.track_key(body("halley")))
    d_h, s_h = op.norm(tr["r"]), op.to_kms(op.norm(tr["v"]))
    i = int(np.argmin(np.abs(d_h - 20)))
    axv.annotate("Halley and Hale-Bopp: bound comets,\nalways below escape speed",
                 (d_h[i], s_h[i]), xytext=(-20, -34), textcoords="offset points", fontsize=7.3,
                 color=ps.INK_SECONDARY, ha="right",
                 arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "lw": 0.7})
    snap = res["snapshot"]["bodies"]
    for p in PLANETS:
        s = snap[p.key]
        axv.plot(s["dist_au"], s["speed_kms"], "o", ms=4.5, color=ps.C_PLANET, zorder=5)
    axv.annotate("planets", (snap["neptune"]["dist_au"], snap["neptune"]["speed_kms"]),
                 xytext=(0, -12), textcoords="offset points", ha="center", fontsize=7,
                 color=ps.INK_SECONDARY)
    axv.set_xscale("log")
    axv.set_yscale("log")
    axv.set_xlim(0.2, 120)
    axv.set_ylim(0.8, 110)
    axv.set_xticks([0.3, 1, 3, 10, 30, 100], ["0.3", "1", "3", "10", "30", "100"])
    axv.set_yticks([1, 3, 10, 30, 100], ["1", "3", "10", "30", "100"])
    axv.minorticks_off()
    # Label the escape line along its own slope, measured in display space.
    x0, x1 = 4.5, 9.0
    p0 = axv.transData.transform((x0, op.to_kms(op.escape_speed(x0))))
    p1 = axv.transData.transform((x1, op.to_kms(op.escape_speed(x1))))
    axv.text(x0, op.to_kms(op.escape_speed(x0)) * 1.07, "escape speed from the Sun",
             rotation=np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0])),
             rotation_mode="anchor", fontsize=7.5, color=ps.INK_SECONDARY, ha="left",
             va="bottom")
    axv.set_xlabel("distance from the Sun (au)")
    axv.set_ylabel("speed relative to the Sun (km/s)")
    axv.set_title("Above the escape line: never coming back")
    axv.text(0.02, 0.02, "bright: inbound · faded: outbound", transform=axv.transAxes,
             fontsize=7, color=ps.MUTED)

    fig.suptitle("The three interstellar objects", x=0.125, ha="left", y=0.99)
    rows = []
    for key, o in iso.items():
        c = o["comparison"]
        rows.append(f"{o['name']}: e = {o['e']:.3f}, v∞ = {o['v_inf_in_kms']:.2f} km/s "
                    f"(published {c['v_inf_kms']['published']:.2f}), radiant "
                    f"{c['radiant_separation_deg'] * 60:.1f}′ from the published one")
    ps.caption(fig, f"{PROVENANCE} Orbit solutions: " + ", ".join(
        f"{o['name']} {o['orbit_solution']}" for o in iso.values()) + ". "
        + "; ".join(rows) + ". v∞ measured from barycentric states more than 100 au out. "
        "ʻOumuamua's published radiant (2017) predates the discovery of its non-gravitational "
        "acceleration, which the current orbit includes.")
    return ps.save(fig, FIGURE_DIR / "04_interstellar_visitors")


# --------------------------------------------------------------------------
# 5. The sky: where visitors come from, where our probes go
# --------------------------------------------------------------------------

def _galactic_xy(ra_deg, dec_deg):
    """Mollweide coordinates in Galactic longitude/latitude, longitude increasing leftwards."""
    icrs = np.array([np.cos(np.radians(dec_deg)) * np.cos(np.radians(ra_deg)),
                     np.cos(np.radians(dec_deg)) * np.sin(np.radians(ra_deg)),
                     np.sin(np.radians(dec_deg))])
    l, b = op.spherical_deg(op.icrs_to_galactic(icrs.T))
    l = np.where(l > 180, l - 360, l)
    return -np.radians(l), np.radians(b)


def figure_sky(res):
    iso = res["interstellar"]["objects"]
    apex = res["interstellar"]["solar_apex"]
    sc = res["spacecraft"]["craft"]
    fig = plt.figure(figsize=(12.5, 6.8))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.grid(True, color=ps.GRID, lw=0.6)
    ax.set_xticks(np.radians([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]))
    ax.set_xticklabels([])
    ax.tick_params(axis="y", labelsize=7, labelcolor=ps.MUTED)
    for lon in (90, 0, 270):
        ax.text(-np.radians(lon if lon <= 180 else lon - 360), np.radians(-83),
                f"l = {lon}°", fontsize=6.5, color=ps.MUTED, ha="center")

    # The ecliptic, as a reference line.
    lam = np.linspace(0, 360, 721)
    ecl = np.column_stack([np.cos(np.radians(lam)), np.sin(np.radians(lam)), np.zeros_like(lam)])
    ra, dec = op.radec_from_ecliptic(ecl)
    ex, ey = _galactic_xy(ra, dec)
    order = np.argsort(ex)
    ax.plot(ex[order], ey[order], color=ps.AXIS, lw=1.0, zorder=1)
    ax.text(*_galactic_xy(*op.radec_from_ecliptic(np.array([np.cos(np.radians(200)),
                                                            np.sin(np.radians(200)), 0.0]))),
            "  ecliptic (plane of the planets)", fontsize=7, color=ps.MUTED, va="bottom")
    ax.text(0.0, np.radians(2), "Galactic centre", fontsize=7, color=ps.MUTED, ha="center",
            va="bottom")
    ax.plot(0, 0, "+", ms=8, color=ps.MUTED, mew=1.0)
    ax.text(np.radians(-178), np.radians(-4), "← plane of the Milky Way →", fontsize=7,
            color=ps.MUTED, va="top")

    ax_x, ax_y = _galactic_xy(apex["ra_deg"], apex["dec_deg"])
    ax.plot(ax_x, ax_y, marker="*", ms=13, color=ps.INK, mec=ps.SURFACE, mew=0.8, zorder=6)
    ax.annotate(f"solar apex ({apex['constellation']}): where the Sun is\nheading relative "
                f"to the nearby stars, {apex['speed_kms']:.0f} km/s", (ax_x, ax_y),
                xytext=(-14, 16), textcoords="offset points", fontsize=7.3, color=ps.INK,
                ha="right")

    offsets = {"oumuamua": (-10, -5), "borisov": (12, -10), "atlas": (-16, -26)}
    for key, o in iso.items():
        x, y = _galactic_xy(o["radiant_ra_deg"], o["radiant_dec_deg"])
        ax.plot(x, y, "o", ms=13 if key == "atlas" else 9, mfc="none" if key == "atlas"
                else ps.C_INTERSTELLAR, mec=ps.C_INTERSTELLAR, mew=2.0 if key == "atlas" else 1.0,
                zorder=7)
        dx, dy = offsets[key]
        extra = ("\n— and where New Horizons is heading, "
                 f"{res['spacecraft']['new_horizons_vs_atlas_radiant_deg']:.1f}° away"
                 if key == "atlas" else "")
        ax.annotate(f"{o['name']} came from here ({o['radiant_constellation']}),\n"
                    f"{o['speed_vs_lsr_kms']:.0f} km/s relative to the local stars{extra}",
                    (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=7.3,
                    color=ps.INK, ha="left" if dx > 0 else "right", va="top")
    offs = {"voyager1": (10, 6, "left"), "voyager2": (10, -12, "left"),
            "pioneer10": (-10, 10, "right"), "pioneer11": (10, -2, "left"),
            "new_horizons": None}
    for key, c in sc.items():
        x, y = _galactic_xy(c["heading_ra_deg"], c["heading_dec_deg"])
        ax.plot(x, y, "s", ms=6.5, color=ps.C_SPACECRAFT, mec=ps.SURFACE, mew=1.0, zorder=8)
        if offs[key] is None:
            continue
        dx, dy, ha = offs[key]
        ax.annotate(f"{c['name']} is heading here ({c['heading_constellation']})", (x, y),
                    xytext=(dx, dy), textcoords="offset points", fontsize=7.3, color=ps.INK,
                    ha=ha, va="center")
    ax.set_title("In and out: where the interstellar visitors came from, where our probes are going",
                 pad=16)
    sep = res["spacecraft"]["new_horizons_vs_atlas_radiant_deg"]
    ps.caption(fig, f"{PROVENANCE} The whole sky in Galactic coordinates: the Milky Way's "
                    "plane runs across the middle, the Galactic centre is at the centre, and "
                    "longitude increases to the left (the astronomers' convention). Orange: radiants "
                    "of the interstellar objects, from their barycentric orbits in 1900. Blue: "
                    "directions of the spacecraft's outgoing asymptotes. Solar apex from the Sun's "
                    "motion relative to the local standard of rest, (U, V, W) = (11.1, 12.24, 7.25) "
                    f"km/s (Schönrich et al. 2010). New Horizons' heading and 3I/ATLAS's radiant "
                    f"are {sep:.1f}° apart — a coincidence: both lie near the ecliptic in Sagittarius.")
    return ps.save(fig, FIGURE_DIR / "05_sky_in_and_out")


# --------------------------------------------------------------------------
# 6. The Voyagers' grand tour
# --------------------------------------------------------------------------

def figure_voyager(res, art):
    sa = art["spacecraft"]
    fly = res["spacecraft"]["flybys"]
    fig, (ax, axv) = plt.subplots(1, 2, figsize=(14.5, 6.6),
                                  gridspec_kw={"width_ratios": [1, 1.2], "wspace": 0.22})
    colors = {"voyager1": ps.BLUE, "voyager2": ps.ORANGE}
    names = {"voyager1": "Voyager 1", "voyager2": "Voyager 2"}
    half = 36

    # (a) The paths, launch to 1995, with the planets where the probes met them.
    for key in ("jupiter", "saturn", "uranus", "neptune"):
        orb = _one_orbit(key, around_jd=_jd("1983-01-01"))
        ax.plot(orb[:, 0], orb[:, 1], color=ps.C_ORBIT, lw=ORBIT_LW, zorder=1)
    for key in ("pioneer10", "pioneer11", "new_horizons"):
        d = sa[key]
        m = d["jd"] <= _jd("2016-01-01")
        ax.plot(d["r"][m, 0], d["r"][m, 1], color=ps.C_CONTEXT, lw=0.9, zorder=2)
        # Name each grey path where it leaves the frame.
        out = np.nonzero(np.max(np.abs(d["r"][m, :2]), axis=1) > half * 0.9)[0]
        if len(out):
            ex, ey = d["r"][m][out[0], :2]
            side = {"pioneer10": (-1.0, 0, "right"), "pioneer11": (-1.0, 0, "right"),
                    "new_horizons": (1.0, 0, "left")}[key]
            ax.text(ex + side[0], ey + side[1], body(key).name, fontsize=6.8, color=ps.MUTED,
                    ha=side[2], va="center")
    for key in ("voyager1", "voyager2"):
        d = sa[key]
        m = d["jd"] <= _jd("1995-01-01")
        ax.plot(d["r"][m, 0], d["r"][m, 1], color=colors[key], lw=1.6, zorder=4,
                label=names[key])
    spots = {}
    for f in fly:
        if f["craft"] not in colors:
            continue
        t = _jd(f["closest_approach_utc"][:10])
        d = sa[f["craft"]]
        p = np.array([np.interp(t, d["jd"], d["r"][:, k_]) for k_ in range(3)])
        ax.plot(p[0], p[1], "o", ms=5.5, mfc=ps.SURFACE, mec=ps.C_PLANET, mew=1.1, zorder=6)
        spots.setdefault(f["planet"], []).append((p, f["closest_approach_utc"][:4], f["craft"]))
    text_at = {"jupiter": (-2, 9, "center"), "saturn": (-17, -6, "right"),
               "uranus": (-12, -21, "right"), "neptune": (9, -26, "left")}
    for planet, hits in spots.items():
        yrs = " & ".join(sorted({h[1] for h in hits}))
        who = "both" if len(hits) > 1 else names[hits[0][2]].replace("oyager ", "")
        tx, ty, ha = text_at[planet]
        for p, _, _ in hits:
            ax.annotate("", (p[0], p[1]), xytext=(tx, ty), textcoords="data",
                        arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.6,
                                    "shrinkA": 0, "shrinkB": 3})
        ax.text(tx, ty, f"{planet.capitalize()} {yrs}\n({who})", fontsize=7, color=ps.INK,
                ha=ha, va="bottom" if ty > 0 else "top")
    _sun(ax, size=4.5)
    _map_axes(ax, half, "The grand tour, 1977–1995")
    ax.legend(loc="upper left")

    # (b) Speed against distance: each flyby a step up (or down).
    rr = np.logspace(0, np.log10(200), 300)
    axv.plot(rr, op.to_kms(op.escape_speed(rr)), color=ps.INK, lw=1.1, ls=(0, (5, 3)), zorder=2)
    for key in ("pioneer10", "pioneer11", "new_horizons"):
        d = sa[key]
        axv.plot(op.norm(d["r"]), op.to_kms(op.norm(d["v"])), color=ps.C_CONTEXT, lw=0.9,
                 zorder=3)
    for key in ("voyager1", "voyager2"):
        d = sa[key]
        axv.plot(op.norm(d["r"]), op.to_kms(op.norm(d["v"])), color=colors[key], lw=1.6,
                 zorder=4, label=names[key])
    # Kick labels sit in the plot's empty regions, tied to the middle of each jump.
    kick_text = {("voyager1", "jupiter"): (6.5, 30.0, "left"),
                 ("voyager1", "saturn"): (12.5, 26.5, "left"),
                 ("voyager2", "uranus"): (23.0, 24.0, "left"),
                 ("voyager2", "jupiter"): (2.5, 7.0, "center"),
                 ("voyager2", "saturn"): (10.5, 3.0, "center"),
                 ("voyager2", "neptune"): (62.0, 9.4, "left")}
    for f in fly:
        if f["craft"] not in colors:
            continue
        d = sa[f["craft"]]
        t = _jd(f["closest_approach_utc"][:10])
        i0 = int(np.searchsorted(d["jd"], t - 12))
        i1 = int(np.searchsorted(d["jd"], t + 12))
        x_ = 0.5 * (op.norm(d["r"][i0]) + op.norm(d["r"][i1]))
        y_ = 0.5 * (op.to_kms(op.norm(d["v"][i0])) + op.to_kms(op.norm(d["v"][i1])))
        tx, ty, ha = kick_text[(f["craft"], f["planet"])]
        dv = f["helio_speed_after_kms"] - f["helio_speed_before_kms"]
        axv.annotate(f"{names[f['craft']]} at {f['planet'].capitalize()}: {dv:+.1f} km/s",
                     (x_, y_), xytext=(tx, ty),
                     textcoords="data", fontsize=7.2, color=ps.INK, ha=ha, va="center",
                     arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.6,
                                 "shrinkA": 1, "shrinkB": 2})
    axv.set_xscale("log")
    axv.set_xlim(1, 200)
    axv.set_ylim(0, 45)
    axv.set_xticks([1, 2, 5, 10, 20, 50, 100, 200], ["1", "2", "5", "10", "20", "50", "100", "200"])
    axv.minorticks_off()
    x0, x1 = 24.0, 48.0
    p0 = axv.transData.transform((x0, op.to_kms(op.escape_speed(x0))))
    p1 = axv.transData.transform((x1, op.to_kms(op.escape_speed(x1))))
    axv.text(x0, op.to_kms(op.escape_speed(x0)) - 0.6, "escape speed from the Sun",
             rotation=np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0])),
             rotation_mode="anchor", fontsize=7.5, color=ps.INK, ha="left", va="top")
    axv.set_xlabel("distance from the Sun (au)")
    axv.set_ylabel("speed relative to the Sun (km/s)")
    axv.set_title("Each flyby is a kick; Jupiter pushed both past escape speed")
    handles, labels = axv.get_legend_handles_labels()
    from matplotlib.lines import Line2D
    handles.append(Line2D([], [], color=ps.C_CONTEXT, lw=0.9))
    labels.append("Pioneer 10 & 11, New Horizons")
    axv.legend(handles, labels, loc="upper right", bbox_to_anchor=(1.0, 0.8))
    v2j = next(f for f in fly if f["craft"] == "voyager2" and f["planet"] == "jupiter")
    axv.text(0.2, 0.985,
             f"How a flyby adds speed: in Jupiter's own frame Voyager 2 arrived at "
             f"{v2j['v_inf_in_kms']:.2f} km/s\nand left at {v2j['v_inf_out_kms']:.2f} km/s — it was "
             f"only turned, by {v2j['turn_observed_deg']:.1f}° (the hyperbola formula\nfor its "
             f"{v2j['closest_approach_km'] / 1e3:,.0f},000 km pass gives "
             f"{v2j['turn_predicted_deg']:.1f}°). Seen from the Sun, the same turn\nadds "
             f"Jupiter's own orbital motion: {v2j['helio_speed_after_kms'] - v2j['helio_speed_before_kms']:+.1f} km/s.",
             transform=axv.transAxes, fontsize=7.2, color=ps.INK_SECONDARY, va="top",
             linespacing=1.35)

    ps.caption(fig, f"{PROVENANCE} Heliocentric speeds every 5 days. Through their planetary "
                    "encounters (Voyager 1 before 1981, Voyager 2 before 1989-08-29) Horizons' "
                    "Voyager trajectories are patched-conic reconstructions of 'rough accuracy', "
                    "built to match the encounters; later arcs are refits of tracking data. Kick "
                    "sizes compare the speed 40 days before and after closest approach.")
    return ps.save(fig, FIGURE_DIR / "06_voyager_gravity_assists")


# --------------------------------------------------------------------------
# 7. Voyager 1's light-day
# --------------------------------------------------------------------------

def _dates(jd):
    """Julian dates (TDB) -> datetimes, for date axes (the 69 s TDB-UTC offset is irrelevant here)."""
    from astropy.time import Time

    return Time(np.atleast_1d(jd), format="jd", scale="tdb").utc.to_datetime()


def figure_light_day(res, art):
    import matplotlib.dates as mdates

    ld = res["light_day"]
    a = art["light_day"]
    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.plot(_dates(a["jd"]), a["rg"], color=ps.C_SPACECRAFT, lw=1.6, zorder=3,
            label="distance from Earth")
    tr = hz.load(hz.track_key(body("voyager1")))
    m = (tr["jd"] >= a["jd"][0]) & (tr["jd"] <= a["jd"][-1])
    ax.plot(_dates(tr["jd"][m]), op.norm(tr["r"][m]), color=ps.C_CONTEXT, lw=1.1, zorder=2,
            label="distance from the Sun")
    L = ld["light_day_au"]
    ax.axhline(L, color=ps.INK_SECONDARY, lw=0.9, zorder=1)
    x0 = _dates(a["jd"][0])[0]
    ax.text(x0, L - 0.1, f"  one light-day = {L:.2f} au = 25,902,068,371 km", fontsize=7.5,
            color=ps.INK_SECONDARY, va="top")
    tc = _dates(a["t_cross"])[0]
    ax.plot(tc, L, "o", ms=8, color=ps.C_SPACECRAFT, mec=ps.SURFACE, mew=1.2, zorder=5)
    diff_ms = abs(ld["difference_s"]) * 1000
    ps.annotate(ax, f"crosses it at {ld['crossing_utc'].replace('T', ' ')} UTC\n"
                    f"NASA's announced time: {ld['published_utc'].replace('T', ' ')} UTC\n"
                    f"difference: {diff_ms:.0f} milliseconds",
                (tc, L), (_dates(a["t_cross"] - 250)[0], L + 2.2), color=ps.INK)
    tm = _dates(a["t_min"])[0]
    ps.annotate(ax, f"closest it comes back: {a['r_min']:.2f} au (April 2027),\n"
                    f"still {ld['margin_light_minutes']:.1f} light-minutes beyond —\n"
                    "the crossing is for good",
                (tm, a["r_min"]), (_dates(a["t_min"] + 40)[0], a["r_min"] - 2.4), color=ps.INK)
    t_sun = _jd(ld["sun_light_day_utc"][:10])
    ps.annotate(ax, f"one light-day from the Sun:\n{ld['sun_light_day_utc'][:10]}",
                (_dates(t_sun)[0], L), (_dates(t_sun - 60)[0], L - 3.6),
                color=ps.INK_SECONDARY)
    ax.set_xlim(_dates(a["jd"][0])[0], _dates(a["jd"][-1])[0])
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.set_ylabel("distance (au)")
    ax.set_title("Voyager 1 becomes the first object one light-day from Earth")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93))
    ax.text(0.99, 0.03, f"Earth's own orbit swings the distance by {ld['annual_swing_au']:.2f} au "
            "a year:\nit shrinks every January to April, while Earth moves towards the spacecraft.",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.3, color=ps.INK_SECONDARY)
    ps.caption(fig, f"{PROVENANCE} Voyager 1 trajectory: 2022 refit of 1981–1992 tracking data "
                    "(R. Jacobson), extrapolated; there has been no two-way tracking since 1992. "
                    "Distance is geometric, Earth's centre to the spacecraft — the definition that "
                    "reproduces NASA's time. Counting light-travel time instead, a command sent from "
                    f"Earth takes 24 h to arrive from {ld['uplink_utc'].replace('T', ' ')} UTC, and a "
                    f"signal received takes 24 h from {ld['downlink_utc'].replace('T', ' ')} UTC.")
    return ps.save(fig, FIGURE_DIR / "07_voyager1_light_day")


# --------------------------------------------------------------------------
# 8. The sky tonight
# --------------------------------------------------------------------------

def figure_sky_tonight(res, art):
    sky = res["sky"]
    sa = art["sky"]
    site = sky["site"]
    off = sa["offset"]
    # Julian days start at noon, so shift by half a day to count from local midnight.
    t_local = (sa["jd_utc"] + off + 0.5 - np.floor(sa["jd_utc"][0] + off + 0.5)) * 24
    fig, ax = plt.subplots(figsize=(11, 5.6))

    # Twilight: contiguous bands between the moments the Sun crosses each depth.
    sun_el = sa["sun_el"]

    def crossing_times(level):
        g = sun_el - level
        idx = np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]
        return [t_local[i] + (t_local[i + 1] - t_local[i]) * g[i] / (g[i] - g[i + 1]) for i in idx]

    levels = [-0.833, -6, -12, -18]
    shades = [0.03, 0.055, 0.08, 0.11]  # civil, nautical, astronomical twilight, night
    for (lo, hi), shade in zip(zip(levels, levels[1:] + [None]), shades):
        a_, b_ = crossing_times(lo), (crossing_times(hi) if hi is not None else [])
        edges = sorted(a_ + b_)
        if hi is None:  # full night: between the evening and morning -18 deg crossings
            edges = sorted(a_)
        for x0, x1 in zip(edges[0::2], edges[1::2]):
            ax.axvspan(x0, x1, color=ps.INK, alpha=shade, lw=0, zorder=0)
            if hi is not None:
                pass
    order = ["saturn", "jupiter", "mars", "uranus", "neptune"]
    colors = dict(zip(order, ps.SERIES[:5]))
    d = sa["moon"]
    ax.plot(t_local, d["el"], color=ps.C_CONTEXT, lw=1.2, zorder=2,
            label=f"Moon — in {sky['bodies']['moon']['constellation']}, up all night")
    i = int(np.nanargmax(d["el"]))
    ax.annotate("Moon", (t_local[i], d["el"][i]), xytext=(0, 5), textcoords="offset points",
                ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    name_offset = {"neptune": (-40, 3, "right"), "saturn": (0, 6, "center"),
                   "uranus": (0, 6, "center"), "mars": (0, 6, "center"),
                   "jupiter": (0, 6, "center")}
    for key in order:
        b = sky["bodies"][key]
        d = sa[key]
        dark_el = np.where(sun_el < -12, d["el"], -99)
        i = int(np.argmax(dark_el))
        ax.plot(t_local, d["el"], color=colors[key], lw=1.7, zorder=3,
                label=f"{b['name']} — magnitude {b['magnitude']:.1f}, in {b['constellation']}; "
                      f"highest in the dark {b['max_alt_dark_deg']:.0f}° at {b['max_alt_time_local']}")
        ax.plot(t_local[i], d["el"][i], "o", ms=5, color=colors[key], mec=ps.SURFACE, zorder=4)
        dx, dy, ha = name_offset[key]
        ax.annotate(b["name"], (t_local[i], d["el"][i]), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, fontsize=7.5, color=ps.INK,
                    arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.5} if dx else None)
    hidden = [b["name"] for k_, b in sky["bodies"].items() if not b["visible"]]
    ax.axhline(0, color=ps.INK_SECONDARY, lw=0.9)
    ax.set_ylim(-5, 90)
    ax.set_xlim(17.5, 31.5)
    ticks = np.arange(18, 32, 1)
    ax.set_xticks(ticks, [f"{int(h % 24):02d}:00" for h in ticks])
    ax.set_xlabel(f"local time ({site['tz_label']}), evening of {SNAPSHOT_DATE} to the next morning")
    ax.set_ylabel("altitude above the horizon (°)")
    ax.set_title(f"What's up tonight from {site['name']} "
                 f"({site['lat']:.2f}° N, {site['lon']:.2f}° E)")
    tw = sky["twilight_local"]
    ax.text(0.01, 0.985, f"sunset {tw['sunset']['evening']} · fully dark (Sun 18° down) "
            f"{tw['astronomical']['evening']}–{tw['astronomical']['morning']} · sunrise "
            f"{tw['sunset']['morning']} · grey shading deepens from twilight to night"
            + (f" · below the horizon all night: {', '.join(hidden)}" if hidden else ""),
            transform=ax.transAxes, va="top", fontsize=7.3, color=ps.INK_SECONDARY)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.94), fontsize=7.2)
    chk = sky["astropy_check"]
    ps.caption(fig, f"{PROVENANCE} Apparent (refracted) altitudes every 10 minutes for the site "
                    f"above. Cross-check against astropy's independent ephemeris for Saturn at "
                    f"{chk['time_utc']} UTC: altitude differs by {abs(chk['alt_diff_arcsec']):.0f}″, "
                    f"azimuth by {abs(chk['az_diff_arcsec']):.0f}″. Magnitudes: lower is brighter; "
                    "the naked-eye limit under a dark sky is about 6.")
    return ps.save(fig, FIGURE_DIR / "08_sky_tonight")


def make_all(results, art) -> list[str]:
    written = []
    written += figure_today(results)
    written += figure_kepler3(results)
    written += figure_halley(results, art)
    written += figure_interstellar(results, art)
    written += figure_sky(results)
    written += figure_voyager(results, art)
    written += figure_light_day(results, art)
    written += figure_sky_tonight(results, art)
    return written
