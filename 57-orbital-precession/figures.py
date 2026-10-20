"""Every static figure, from ``data/results.json`` and ``data/derived/series.npz``.

    python figures.py            # all figures
    python figures.py ladder s2  # only some
"""

from __future__ import annotations

import json
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import plotstyle as PS
from analysis import SERIES_PATH
from config import FIGURE_DIR, PUBLISHED, RESULTS_PATH

PS.use_style()

PLANET_NAMES = {"venus": "Venus", "emb": "Earth-Moon", "mars": "Mars", "jupiter": "Jupiter",
                "saturn": "Saturn", "uranus": "Uranus", "neptune": "Neptune", "pluto": "Pluto",
                "mercury": "Mercury"}


def load():
    res = json.loads(RESULTS_PATH.read_text())
    with np.load(SERIES_PATH) as z:
        ser = {k: z[k] for k in z.files}
    return res, ser


def year(t_cy):
    return 2000.0 + 100.0 * np.asarray(t_cy)


# --------------------------------------------------------------------------
# 1. The budget
# --------------------------------------------------------------------------

def fig_budget(res, ser):
    c = res["mercury"]["contributions"]
    an = res["mercury"]["analytic"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [1.65, 1]})
    p = c["planets"]
    steps = [("Venus", p["venus"]["rate"], PS.C_PLANETS),
             ("Jupiter", p["jupiter"]["rate"], PS.C_PLANETS),
             ("Earth-Moon", p["emb"]["rate"], PS.C_PLANETS),
             ("Saturn", p["saturn"]["rate"], PS.C_PLANETS),
             ("Mars", p["mars"]["rate"], PS.C_PLANETS),
             ("Uranus +\nNeptune", p["uranus"]["rate"] + p["neptune"]["rate"] + p["pluto"]["rate"],
              PS.C_PLANETS),
             ("all at\nonce", c["nonadditivity"], PS.C_PLANETS),
             ("general\nrelativity", c["gr_in_system"]["gr_full"]["rate"], PS.C_GR),
             ("Sun's\nflattening", c["j2"]["rate"], PS.C_SUN),
             ("Sun's\nspin", c["lense_thirring"]["rate"], PS.C_SUN)]
    cum = 0.0
    for k, (label, v, col) in enumerate(steps):
        ax.bar(k, v, bottom=cum, width=0.72, color=col, linewidth=0)
        txt = f"{v:+.3f}" if abs(v) < 1 else f"{v:.2f}"
        ax.text(k, cum + max(v, 0) + 6, txt, ha="center", va="bottom", fontsize=7.5,
                color=PS.INK_SECONDARY)
        if k < len(steps) - 1:
            ax.plot([k + 0.36, k + 0.64], [cum + v, cum + v], color=PS.AXIS, lw=0.8)
        cum += v
    kt = len(steps)
    ax.bar(kt, cum, width=0.72, color=PS.INK, linewidth=0)
    ax.text(kt, cum + 6, f"{cum:.3f}", ha="center", va="bottom", fontsize=8, color=PS.INK,
            fontweight="bold")
    ax.set_xticks(range(kt + 1), [s[0] for s in steps] + ["our\ntotal"], fontsize=7.5)
    ax.set_ylabel("turn of Mercury's perihelion (arcsec per century)")
    ax.set_ylim(0, 640)
    ax.set_xlim(-0.6, kt + 0.6)
    de, pk = c["de440_total"], c["park2017_total"]
    ax.plot([kt - 0.5, kt + 0.5], [de, de], color=PS.C_PUBLISHED, lw=1.2)
    pr = res["mercury"].get("park_recipe", {})
    lines = [f"ours, all forces above: {cum:.4f}\"/cy",
             f"+ asteroid belt (estimate): {c['full_model_plus_asteroids']:.4f}",
             f"JPL's DE440, same fit: {de:.4f}"]
    if pr:
        lines.append(f"ours with Park et al.'s fit (AD 1000-3000): {pr['plus_asteroid_estimate']:.4f}")
    lines.append(f"MESSENGER (Park et al. 2017): {pk:.4f} +- 0.0015")
    ax.text(0.30, 0.60, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=7.5,
            color=PS.INK_SECONDARY, linespacing=1.5)
    newton = c["all_newton"]["rate"]
    PS.annotate(ax, f"Newton's planets: {newton:.2f}\"/cy", (6.3, newton), (3.9, 470))
    ax.set_title("Mercury's perihelion budget, measured force by force", loc="left")
    ax.text(0, 1.01, "Each bar = one force switched on in our own N-body runs (in-plane rate at J2000)",
            transform=ax.transAxes, fontsize=8, color=PS.INK_SECONDARY)

    # Right: each planet's share, defined four ways, against Park et al. (2017)'s Table 3.
    from config import PARK2017_TABLE3 as T3
    keys = ["venus", "jupiter", "emb", "saturn", "mars", "uranus", "neptune"]
    y = np.arange(len(keys))[::-1]
    ref = np.array([T3["planets"][k] for k in keys])
    inter = {k: sum(v for kk, v in T3["interactions"].items() if k in kk.split("+")) for k in keys}
    series = [
        ("our N-body: that planet alone (Sun + Mercury + it)", np.array([p[k]["rate"] for k in keys]),
         dict(marker="o", s=34, color=PS.INK)),
        ("our N-body: all planets minus that one", np.array([res["mercury"]["leave_one_out"][k] for k in keys]),
         dict(marker="o", s=40, facecolor="none", edgecolor=PS.INK, lw=1.2)),
        ("Gauss's averaging (1818): first order in the masses", np.array([an["gauss"][k] for k in keys]),
         dict(marker="D", s=30, color=PS.BLUE)),
        ("Laplace-Lagrange: also low order in eccentricity",
         np.array([an["laplace_lagrange"][k]["total"] for k in keys]),
         dict(marker="s", s=30, facecolor="none", edgecolor=PS.MUTED, lw=1.2)),
    ]
    ax2.axvline(0, color=PS.C_PUBLISHED, lw=1.0)
    for label, vals, kw in series:
        ax2.scatter(100 * (vals / ref - 1), y, zorder=3, label=label, **kw)
    ax2.scatter(100 * np.array([inter[k] for k in keys]) / ref, y, marker="|", s=120,
                color=PS.C_PUBLISHED, zorder=2, label="Park et al.: + indirect (via other orbits)")
    ax2.set_yticks(y, [f"{PLANET_NAMES[k]}  {T3['planets'][k]:.4f}\"" for k in keys], fontsize=7.5)
    ax2.tick_params(axis="y", length=0)
    ax2.set_xlim(-2.2, 5.0)
    ax2.set_xlabel("difference from Park et al. 2017's share (%)")
    ax2.set_title("Each planet's share, defined four ways", loc="left")
    ax2.text(0.0, 1.01, "violet line = Park et al. 2017, Table 3 (MESSENGER-era JPL ephemeris)",
             transform=ax2.transAxes, fontsize=7.5, color=PS.C_PUBLISHED)
    ax2.legend(loc="lower center", fontsize=7, bbox_to_anchor=(0.42, -0.42))
    ax2.grid(axis="y", visible=False)
    PS.caption(fig, "Simulation from DE440's J2000 state (REBOUND IAS15 + REBOUNDx), 1550-2650. "
               "Shares: Sun + Mercury + one planet; 'all at once' = all planets together minus the sum "
               "of the shares. GR from the Einstein-Infeld-Hoffmann equations; the Sun's J2 and spin "
               "from DE440's constants. Rates are the turn of the perihelion within Mercury's own "
               "orbital plane (the frame-independent definition used by MESSENGER).")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "01_budget")


# --------------------------------------------------------------------------
# 2. Newton vs Einstein vs JPL
# --------------------------------------------------------------------------

def fig_vs_de440(res, ser):
    yr = year(ser["t_cy"])
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.6), sharex=True,
                             gridspec_kw={"height_ratios": [1.2, 1, 1]})
    ax = axes[0]
    ax.plot(yr, ser["offset_all_newton"], color=PS.C_NEWTON, lw=1.4, label="Newton only")
    ax.plot(yr, ser["offset_all_gr_full"], color=PS.C_GR, lw=1.4, label="Newton + general relativity")
    ax.axhline(0, color=PS.AXIS, lw=0.8)
    ax.set_ylabel("perihelion direction\nminus DE440 (arcsec)")
    ax.set_title("Our Mercury against JPL's, 1550-2650", loc="left")
    ax.legend(loc="upper right")
    PS.annotate(ax, "falls behind by 43\" per century", (2500, ser["offset_all_newton"][
        np.argmin(np.abs(yr - 2500))]), (2190, -200))
    ax = axes[1]
    ax.plot(yr, ser["offset_all_gr_full"], color=PS.C_GR, lw=1.0, label="+ general relativity")
    ax.plot(yr, ser["offset_full_model"], color=PS.C_SUN, lw=1.0,
            label="+ the Sun's flattening (J2) and spin (Lense-Thirring)")
    ax.axhline(0, color=PS.AXIS, lw=0.8)
    ax.set_ylabel("zoom (arcsec)")
    ax.legend(loc="lower left")
    ax.set_ylim(-0.22, 0.17)
    ax = axes[2]
    for k, col, lab in (("all_newton", PS.C_NEWTON, "Newton only"),
                        ("all_gr_full", PS.C_GR, "+ GR"),
                        ("full_model", PS.C_SUN, "+ Sun's J2 and spin")):
        ax.semilogy(yr, np.maximum(ser[f"km_{k}"], 1e-3), color=col, lw=1.0, label=lab)
    ax.text(1600, 2.5e5, "Newton only", color=PS.INK_SECONDARY, fontsize=8)
    ax.text(1600, 1.5e2, "+ GR", color=PS.C_GR, fontsize=8)
    ax.text(1600, 1.2e0, "+ the Sun's J2 and spin", color=PS.C_SUN, fontsize=8)
    ax.set_ylabel("Mercury's distance from\nits DE440 position (km)")
    ax.set_ylim(1e-2, 1e6)
    ax.set_xlabel("year")
    ax.set_xlim(1550, 2650)
    PS.caption(fig, "Simulation started from DE440 at J2000.0 and run backwards to 1550 and forwards "
               "to 2650. What is left once GR, J2 and the Sun's spin are in (0.01\", tens of km) "
               "is mostly the asteroid belt, which DE440 includes and our runs do not.")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "02_newton_vs_einstein")


# --------------------------------------------------------------------------
# 3. Which frame?
# --------------------------------------------------------------------------

def fig_frames(res, ser):
    yr = year(ser["t_cy"])
    fr = res["mercury"]["frames"]["de440"]
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    specs = [("frame_psi", PS.INK, f"within Mercury's own orbital plane: {fr['psi']:.2f}\"/cy"),
             ("frame_varpi_ecl", PS.BLUE,
              f"longitude of perihelion, fixed J2000 ecliptic: {fr['varpi_ecliptic_j2000']:.2f}\"/cy"),
             ("frame_varpi_icrf", PS.ORANGE,
              f"longitude of perihelion, ICRF equator: {fr['varpi_icrf_equator']:.2f}\"/cy")]
    for k, col, lab in specs:
        ax.plot(yr, ser[k], color=col, lw=1.0, label=lab)
    ax.axhline(0, color=PS.AXIS, lw=0.8)
    ax.set_xlim(1550, 2650)
    ax.set_xlabel("year")
    ax.set_ylabel("angle minus 575\" per century (arcsec)")
    ax.set_title("One orbit, three precession rates: it depends on where you measure from", loc="left")
    ax.legend(loc="lower left", fontsize=7.5)
    ax.text(0.99, 0.97, f"Mercury's node slides along the ecliptic at {fr['Omega_ecliptic']:.0f}\"/cy;\n"
            f"the 'dogleg' angle Omega + omega then gains (1 - cos i) x dOmega/dt\n"
            f"= {fr['dogleg_predicted']:.2f}\"/cy relative to the in-plane turn (i = 7.0 deg).",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=PS.INK_SECONDARY)
    PS.caption(fig, "Mercury in DE440 (with all of GR), the same orbit measured three ways. The in-plane "
               "angle is frame-independent; the textbook numbers (43\"/cy, ~532\"/cy, 575.31\"/cy) "
               "refer to it.")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "03_which_frame")


# --------------------------------------------------------------------------
# 4. The perihelion on the march
# --------------------------------------------------------------------------

def fig_perihelia(res, ser):
    r0, v0 = ser["orbit_j2000"][0], ser["orbit_j2000"][1]
    h = np.cross(r0, v0)
    from orbit_physics import elements
    import ephemeris as eph
    el = elements(r0, v0, eph.gm("sun") + eph.gm("mercury"))
    P = el["ehat"][:, 0]
    Q = np.cross(h / np.linalg.norm(h), P)
    a, e = float(el["a"][0]), float(el["e"][0])
    km = eph.au_km() / 1e6  # million km per au
    fig, (ax0, ax) = plt.subplots(1, 2, figsize=(12, 5.4), gridspec_kw={"width_ratios": [1, 1.35]})
    E = np.linspace(0, 2 * np.pi, 400)
    x = a * (np.cos(E) - e) * km
    y = a * np.sqrt(1 - e**2) * np.sin(E) * km
    ax0.plot(x, y, color=PS.AXIS, lw=1.2)
    ax0.plot(0, 0, "o", color=PS.YELLOW, ms=11)
    ax0.plot(a * (1 - e) * km, 0, "o", color=PS.INK, ms=4)
    ax0.annotate("perihelion\n(46 million km)", (a * (1 - e) * km, 0), (20, 22), fontsize=8,
                 color=PS.INK_SECONDARY, arrowprops={"arrowstyle": "-", "color": PS.INK_SECONDARY,
                                                     "lw": 0.8})
    ax0.set_aspect("equal")
    ax0.set_xlabel("million km")
    ax0.set_ylabel("million km")
    ax0.set_title("Mercury's orbit (2000)", loc="left")

    # Right: how far the perihelion point travels along the orbit, year by year.
    tracks = {}
    for label in ("de440", "newton"):
        eh = ser[f"peri_e_{label}"]
        theta = np.unwrap(np.arctan2(eh @ Q, eh @ P))
        tracks[label] = (2000 + ser[f"peri_t_{label}"] / 365.25, theta * ser[f"peri_q_{label}"] * km)
    for label, col, lab in (("newton", PS.C_NEWTON, "our Newton-only simulation (planets only)"),
                            ("de440", PS.C_GR, "JPL's DE440 (with general relativity)")):
        yrs, s_ = tracks[label]
        ax.scatter(yrs, s_ * 1000, s=2, color=col, lw=0, label=lab)
        ax0.scatter(a * (1 - e) * km * np.cos(s_ / (a * (1 - e) * km)), s_, s=1, color=col, lw=0)
    ax.axhline(0, color=PS.AXIS, lw=0.8)
    y_de, s_de = tracks["de440"]
    y_nw, s_nw = tracks["newton"]
    gap = np.interp(2650, y_de, s_de) - np.interp(2650, y_nw, s_nw)
    PS.annotate(ax, f"{gap * 1000:,.0f} thousand km apart in 2650:\nEinstein's 43\" per century",
                (2640, np.interp(2640, y_de, s_de) * 1000 - 0.5 * gap * 1000), (2240, 200))
    ax.set_xlim(1550, 2650)
    ax.set_xlabel("year")
    ax.set_ylabel("perihelion position along the orbit,\nrelative to 2000 (thousand km)")
    ax.set_title("The perihelion point travels 1.4 million km in 1,100 years", loc="left")
    ax.legend(loc="upper left", markerscale=5)
    PS.caption(fig, "Every perihelion passage 1550-2650 (4,566 of them), found from the 2-day samples "
               "(parabolic minimum of the distance) and projected on Mercury's J2000 orbital plane. The "
               "perihelion swings through 1.8 degrees; without GR it would fall 0.13 degrees short by "
               "2650. The scatter is the planets' short-period tugs (~7\").")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "04_perihelion_march")


# --------------------------------------------------------------------------
# 5. The anomaly through history
# --------------------------------------------------------------------------

def fig_history(res, ser):
    from config import MERCURY_HISTORY as H
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10, 8.2), sharex=True,
                                  gridspec_kw={"height_ratios": [1.6, 1]})
    gr = 42.98
    ax.axhline(gr, color=PS.C_GR, lw=1.0)
    ax.text(1852, gr + 0.25, "general relativity: 42.98\"/cy", color=PS.C_GR, fontsize=8)
    offsets = {"Le Verrier": (8, -2), "Newcomb": (8, 6)}
    special = {("Clemence", 1943): (-8, 4, "right"), ("Clemence", 1947): (8, -12, "left")}
    for h in H:
        filled = h["quality"] == "primary"
        kw = dict(fmt="o", ms=6, lw=1.1, capsize=0, color=PS.INK,
                  mfc=PS.INK if filled else PS.SURFACE, mec=PS.INK)
        ax.errorbar(h["year"], h["value"], yerr=h["err"] or None, **kw)
        if h["err"]:
            ax2.semilogy(h["year"], h["err"], "o", ms=6, color=PS.INK,
                         mfc=PS.INK if filled else PS.SURFACE, mec=PS.INK)
    for h in H:
        who = h["who"].split(" (")[0]
        dx, dy = offsets.get(who, (7, 4))
        ha = "left"
        if (who, h["year"]) in special:
            dx, dy, ha = special[(who, h["year"])]
        if h["year"] in (1976, 1987, 1991, 1992, 2001, 2017, 1958, 1972):
            continue
        ax.annotate(f"{who} {h['year']}", (h["year"], h["value"]), (dx, dy), textcoords="offset points",
                    fontsize=7.5, color=PS.INK_SECONDARY, ha=ha)
    ax.annotate("radar era: Shapiro, Anderson et al.,\nPitjeva's ephemerides, MESSENGER (Park et al. 2017)",
                (1990, 43.2), (1905, 45.2), fontsize=7.5, color=PS.INK_SECONDARY,
                arrowprops={"arrowstyle": "-", "color": PS.INK_SECONDARY, "lw": 0.6})
    ax.set_ylim(36.5, 46.5)
    ax.set_ylabel("Mercury's unexplained perihelion advance\n(arcsec per century)")
    ax.set_title("The anomaly that Newton could not explain, measured for 160 years", loc="left")
    ax.text(1915, 37.2, "1915: Einstein\npredicts 43\"", fontsize=7.5, color=PS.C_GR, ha="center")
    ax.axvline(1915, color=PS.C_GR, lw=0.6)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=PS.INK, label="checked in the original paper"),
                       Line2D([], [], marker="o", ls="", mfc=PS.SURFACE, mec=PS.INK,
                              label="from a later compilation")], loc="lower right", fontsize=7.5)
    ax2.set_ylabel("uncertainty (arcsec per century)")
    ax2.set_xlabel("year of publication")
    ax2.set_ylim(3e-4, 5)
    ax2.set_xlim(1850, 2025)
    ax2.axvline(1915, color=PS.C_GR, lw=0.6)
    ax2.text(2016, 0.0009, "MESSENGER\n0.0009", fontsize=7.5, color=PS.INK_SECONDARY, ha="right", va="center")
    ax2.text(1897, 2.09, "Newcomb 2.1", fontsize=7.5, color=PS.INK_SECONDARY, va="center")
    PS.caption(fig, "Le Verrier (1859) and Newcomb (1882) gave no errors. Newcomb 1895 is his e x excess "
               "(8.48 +- 0.43) divided by e; Clemence 1943 is GR + his O - C (-0.07 +- 0.41 probable error). "
               "Radar results quoted as a ratio to GR are multiplied by 42.98. Park et al.'s value is the "
               "relativistic term of their fit. Sources in config.MERCURY_HISTORY.")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "05_history")


# --------------------------------------------------------------------------
# 6. GR across the planets
# --------------------------------------------------------------------------

def fig_planets(res, ser):
    gr = res["mercury"]["gr_by_planet"]
    pred = res["mercury"]["analytic"]["schwarzschild_by_planet"]
    import ephemeris as eph
    from config import J2000_JD
    from orbit_physics import elements
    keys = ["mercury", "venus", "emb", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto"]
    a = []
    for k in keys:
        r, v = eph.heliocentric_state(k, J2000_JD)
        a.append(float(elements(r, v, eph.gm("sun") + eph.gm(k))["a"][0]))
    a = np.array(a)
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    aa = np.logspace(np.log10(0.3), np.log10(45), 200)
    import orbit_physics as op
    c = eph.c_au_per_day()
    ax.loglog(aa, [op.schwarzschild_rate(eph.gm("sun"), x, 0.0, c) for x in aa], color=PS.C_GR,
              lw=1.0, label="GR, circular orbits: 6 pi GM/(c^2 a) per orbit")
    ax.scatter(a, [pred[k] for k in keys], s=60, facecolor="none", edgecolor=PS.INK, lw=1.0,
               zorder=3, label="GR formula for the actual orbits")
    ax.scatter(a, [max(gr[k]["rate"], 1e-5) for k in keys], s=16, color=PS.C_GR, zorder=4,
               label="our N-body: (with GR) minus (without)")
    for k, x in zip(keys, a):
        ax.annotate(f"{PLANET_NAMES[k]}  {pred[k]:.3g}", (x, pred[k]), (8, 4), textcoords="offset points",
                    fontsize=7.5, color=PS.INK_SECONDARY)
    ax.set_xlabel("semi-major axis (au)")
    ax.set_ylabel("relativistic perihelion advance (arcsec per century)")
    ax.set_title("General relativity turns every planet's orbit, Mercury's most", loc="left")
    ax.legend(loc="lower left", fontsize=7.5)
    PS.caption(fig, "N-body values from the 1550-2650 runs sampled every 20 days; nearly circular "
                    "orbits (Venus, Neptune) make the perihelion direction noisy, hence small scatter.")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "06_gr_every_planet")


# --------------------------------------------------------------------------
# 7. The ladder
# --------------------------------------------------------------------------

def fig_ladder(res, ser):
    import ladder
    rungs = res["ladder"]["rungs"]
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(10.5, 9.0), sharex=True,
                                  gridspec_kw={"height_ratios": [1.6, 1]})
    xx = np.logspace(-10.6, np.log10(0.166), 400)
    ax.loglog(xx, 6 * np.pi * xx, color=PS.INK, lw=1.0, label="GR, weak field: 6 pi x")
    ax.loglog(xx, ladder.exact_circular(xx), color=PS.MUTED, lw=1.0,
              label="GR, exact (circular orbits): 2 pi [(1 - 6x)^-1/2 - 1]")
    ax.axvline(1 / 6, color=PS.MUTED, lw=0.8)
    ax.text(1 / 6 * 0.8, 3e-9, "innermost stable orbit (x = 1/6)", rotation=90, ha="right",
            va="bottom", fontsize=7.5, color=PS.INK_SECONDARY)
    for p in res["ladder"]["planets_predicted"]:
        ax.scatter(p["x"], p["gr_per_orbit"], s=22, facecolor="none", edgecolor=PS.C_GROUP["Solar System"],
                   lw=0.9, zorder=3)
        ax.annotate(p["name"], (p["x"], p["gr_per_orbit"]), (-6, 6), textcoords="offset points",
                    fontsize=7, color=PS.MUTED, ha="right")
    groups_seen = set()
    label_pos = {"Mercury": (8, -12), "Icarus (asteroid)": (8, 6), "LAGEOS II (satellite)": (8, -3),
                 "Hulse-Taylor pulsar": (-8, 8), "Double pulsar": (8, -10), "S2 at Sgr A*": (-10, 8),
                 "DI Her": (8, -4)}
    for r in rungs:
        col = PS.C_GROUP[r["group"]]
        y = r["ratio"] * r["gr_per_orbit"]
        yerr = r["ratio_err"] * r["gr_per_orbit"]
        lab = r["group"] if r["group"] not in groups_seen else None
        groups_seen.add(r["group"])
        ax.errorbar(r["x"], y, yerr=yerr, fmt="o", ms=5, color=col, mec=PS.SURFACE, mew=0.6, lw=1,
                    zorder=4, label=lab)
        if r["name"] in label_pos or r["group"] != "Binary stars":
            off = label_pos.get(r["name"], (8, 0))
            ax.annotate(r["name"], (r["x"], y), off, textcoords="offset points", fontsize=7.5,
                        color=PS.INK_SECONDARY, ha="right" if off[0] < 0 else "left")
    bs = [r for r in rungs if r["group"] == "Binary stars"]
    ax.annotate("eclipsing binaries\n(TESS)", (np.median([r["x"] for r in bs]), 4e-6), (-60, 26),
                textcoords="offset points", fontsize=7.5, color=PS.INK_SECONDARY)
    ax.set_ylabel("relativistic turn of the orbit\nper revolution (radians)")
    ax.set_ylim(2e-9, 3)
    ax.set_title("The precession ladder: one formula from an Earth satellite to a black hole", loc="left")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.text(0.30, 0.03, "open circles: other planets, predicted (not yet measured)",
            transform=ax.transAxes, ha="left", fontsize=7.5, color=PS.MUTED)

    axr.axhline(1, color=PS.INK, lw=0.9)
    axr.axhline(0, color=PS.MUTED, lw=0.8)
    axr.text(2.5e-10, 0.04, "Newton (no relativistic precession)", fontsize=7.5, color=PS.MUTED,
             va="bottom")
    for r in rungs:
        col = PS.C_GROUP[r["group"]]
        axr.errorbar(r["x"], r["ratio"], yerr=r["ratio_err"], fmt="o", ms=5, color=col,
                     mec=PS.SURFACE, mew=0.6, lw=1, zorder=4)
    best = [r for r in rungs if r["ratio_err"] < 0.03]
    axr.text(0.99, 0.66, "most precise:\n" + "\n".join(
        f"{r['name']}: {r['ratio']:.5f} +- {r['ratio_err']:.5f}" for r in sorted(best, key=lambda r: r["ratio_err"])),
        transform=axr.transAxes, ha="right", va="top", fontsize=7.5, color=PS.INK_SECONDARY, linespacing=1.4)
    axr.set_xscale("log")
    axr.set_ylim(-0.25, 1.6)
    axr.set_ylabel("measured / GR prediction")
    axr.set_xlabel("depth in gravity, x = GM / (c^2 a (1 - e^2))")
    axr.set_xlim(2e-10, 0.3)
    PS.caption(fig, "Measured relativistic precession of each orbit (classical parts removed by the "
               "observers): Mercury from MESSENGER (Park et al. 2017, with Cassini's gamma); Icarus "
               "(Shapiro et al. 1971); LAGEOS II (Lucchesi & Peron 2010); eclipsing binaries from "
               "Baroch et al. 2021 and Claret et al. 2021 (DI Her with its tilted spins); binary pulsars "
               "with masses from other relativistic effects (this project); S2 from GRAVITY 2020.")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "07_precession_ladder")


# --------------------------------------------------------------------------
# 8. S2
# --------------------------------------------------------------------------

def fig_s2(res, ser):
    import s2 as S2
    fig = plt.figure(figsize=(13, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], width_ratios=[1.05, 1, 1])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[1, :2])
    ax4 = fig.add_subplot(gs[1, 2])
    yrs = ser["s2_years"]
    N, E = ser["s2_newton"], ser["s2_einstein"]
    one = (yrs >= 2010.35) & (yrs <= 2010.35 + 16.05)
    ax0.plot(E[one, 0], E[one, 1], color=PS.INK, lw=1.0)
    ax0.plot(0, 0, "+", color=PS.INK, ms=10)
    ax0.set_xlim(60, -120)
    ax0.set_aspect("equal")
    ax0.set_xlabel("RA offset (mas, east to the left)")
    ax0.set_ylabel("Dec offset (mas)")
    ax0.set_title("S2's orbit (GRAVITY 2020)", loc="left")
    ia = np.argmin(np.abs(yrs - res["s2"]["apocentre_year"]))
    ax0.add_patch(plt.Circle((E[ia, 0], E[ia, 1]), 6, fill=False, ec=PS.ORANGE, lw=1))
    ax0.annotate("2026 apocentre", (E[ia, 0], E[ia, 1]), (-30, -60), textcoords="offset points",
                 fontsize=7.5, color=PS.INK_SECONDARY,
                 arrowprops={"arrowstyle": "-", "color": PS.INK_SECONDARY, "lw": 0.6})
    ax0.annotate("2018 pericentre\n(120 au from the black hole)", (0, 0), (30, 50),
                 textcoords="offset points", fontsize=7.5, color=PS.INK_SECONDARY,
                 arrowprops={"arrowstyle": "-", "color": PS.INK_SECONDARY, "lw": 0.6})

    # Zoom on the 2026 apocentre: Newton vs GR, with error circles.
    win = (yrs > 2025.6) & (yrs < 2027.2)
    ax1.plot(N[win, 0], N[win, 1], color=PS.C_NEWTON, lw=1.4, label="Newton")
    ax1.plot(E[win, 0], E[win, 1], color=PS.C_GR, lw=1.4, label="general relativity")
    ax1.plot(N[ia, 0], N[ia, 1], "o", color=PS.C_NEWTON, ms=4)
    ax1.plot(E[ia, 0], E[ia, 1], "o", color=PS.C_GR, ms=4)
    cx, cy = E[ia, 0], E[ia, 1]
    ax1.add_patch(plt.Circle((cx + 0.9, cy - 1.0), 0.40, fill=False, ec=PS.INK_SECONDARY, lw=0.8))
    ax1.text(cx + 0.9, cy - 1.5, "NACO\nerror", ha="center", va="top", fontsize=7, color=PS.INK_SECONDARY)
    ax1.add_patch(plt.Circle((cx - 0.4, cy - 1.0), 0.065, fill=True, color=PS.INK_SECONDARY))
    ax1.text(cx - 0.4, cy - 1.2, "GRAVITY\nerror", ha="center", va="top", fontsize=7, color=PS.INK_SECONDARY)
    ax1.set_xlim(cx + 1.9, cx - 1.6)
    ax1.set_ylim(cy - 2.3, cy + 1.2)
    ax1.set_aspect("equal")
    ax1.set_xlabel("RA offset (mas)")
    ax1.set_title(f"Zoom: GR moves the apocentre {res['s2']['shift_at_apocentre_mas']:.2f} mas", loc="left")
    ax1.legend(loc="upper left", fontsize=7.5)

    # Rosette: S2 with GR boosted.
    ros = ser["s2_rosette_sky"]
    ax2.plot(ros[:, 0], ros[:, 1], color=PS.C_GR, lw=0.6)
    ax2.plot(0, 0, "+", color=PS.INK, ms=8)
    pad = 10
    ax2.set_xlim(ros[:, 0].max() + pad, ros[:, 0].min() - pad)
    ax2.set_ylim(ros[:, 1].min() - pad, ros[:, 1].max() + pad)
    ax2.set_aspect("equal")
    ax2.set_xlabel("RA offset (mas)")
    ax2.set_title(f"12 orbits with GR x{float(ser['s2_rosette_boost']):.0f} (simulation)", loc="left")

    d = ser["s2_diff_mas"]
    ax3.plot(yrs, d[:, 0], color=PS.BLUE, lw=1.2, label="RA (GR - Newton)")
    ax3.plot(yrs, d[:, 1], color=PS.ORANGE, lw=1.2, label="Dec (GR - Newton)")
    for name, (t0, t1, err) in S2.INSTRUMENT_ERRORS_MAS.items():
        t1 = min(t1, 2030)
        ax3.fill_between([t0, t1], -err, err, color=PS.GRID, alpha=0.9, lw=0, zorder=0)
        ax3.text(0.5 * (t0 + t1), (err if err < 1.5 else 1.45) + 0.03, name.split(",")[0], ha="center",
                 va="bottom", fontsize=7, color=PS.INK_SECONDARY)
    for tp in (2002.33, 2018.38, 2034.4):
        ax3.axvline(tp, color=PS.AXIS, lw=0.8)
    ax3.set_xlim(1992, 2036)
    ax3.set_ylim(-1.6, 1.6)
    ax3.set_xlabel("year")
    ax3.set_ylabel("precession signal on the sky (mas)")
    ax3.set_title("The signal against each instrument's per-point error (grey bands: +-1 sigma)", loc="left")
    ax3.legend(loc="lower left", fontsize=7.5)

    g = res["s2"]["gravity"]
    ys = list(np.arange(len(g))[::-1].astype(float))
    labels = []
    for yy, m in zip(ys, g):
        col = PS.YELLOW if m["stars"] == "S2" else PS.VIOLET
        ax4.errorbar(m["f_sp"], yy, xerr=m["err"], fmt="o", color=col, ms=5, lw=1)
        labels.append(f"{m['label']}, {'S2 alone' if m['stars'] == 'S2' else 'S2 + 3 stars'}")
    h1 = res["s2"].get("h1_public_data")
    if h1 and h1.get("f_sp") is not None:
        ys.append(-1.3)
        ax4.errorbar(h1["f_sp"], -1.3, xerr=h1["f_sp_err"], fmt="o", color=PS.MUTED, ms=5, lw=1)
        labels.append("public Keck + VLT data\n(project H1): spurious")
    ax4.axvline(1, color=PS.C_GR, lw=0.9)
    ax4.axvline(0, color=PS.C_NEWTON, lw=0.9)
    ax4.set_xlim(-0.6, 6.3)
    ax4.set_yticks(ys, labels, fontsize=7)
    ax4.tick_params(axis="y", length=0)
    ax4.set_xlabel("f_SP (0 = Newton, 1 = GR)")
    ax4.set_title("Measured precession factor", loc="left")
    ax4.grid(axis="y", visible=False)
    PS.caption(fig, "Our own Newtonian and first-post-Newtonian integrations of GRAVITY's 2020 orbit "
               "(osculating at the 2010.35 apocentre; distance 8.247 kpc). GRAVITY quotes 0.78 mas at "
               "apocentre. Per-point errors from GRAVITY Collaboration 2020 (speckle 3.8 mas, NACO 0.4 "
               "mas, GRAVITY 65 microarcsec).")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "08_s2")


# --------------------------------------------------------------------------
# 9. Eclipsing binaries and DI Her
# --------------------------------------------------------------------------

def fig_binaries(res, ser):
    import binaries as B
    b = res["binaries"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1.15, 1, 1.1]})
    ax = axes[0]
    names = ["V541 Cyg", "EY Cep", "V459 Cas", "AL Dor"]
    cols = [PS.BLUE, PS.ORANGE, PS.AQUA, PS.YELLOW]
    for name, col in zip(names, cols):
        tim = np.array(b["tess_timings"][name], dtype=object)
        n = np.array([t[0] for t in tim], float)
        d = np.array([t[1] for t in tim], float)
        e = np.array([t[2] for t in tim], float)
        fit = b["refits"][name]
        d0 = np.average(d, weights=1 / e**2)
        n0 = np.average(n, weights=1 / e**2)
        ax.errorbar(n, (d - d0) * 1440, yerr=e * 1440, fmt="o", ms=3.5, color=col, lw=0.8,
                    label=f"{name}: {fit['wdot_deg_per_cycle'] * 1e4:.2f} +- {fit['wdot_err'] * 1e4:.2f}")
        nn = np.linspace(n.min(), n.max(), 2)
        ax.plot(nn, fit["slope_days_per_cycle"] * (nn - n0) * 1440, color=col, lw=0.9)
    ax.set_xlabel("orbital cycle")
    ax.set_ylabel("T2 - T1, minus its mean (minutes)")
    ax.set_title("Secondary eclipses drift as the orbit turns", loc="left")
    ax.legend(fontsize=7, title="our rate (1e-4 deg per orbit)", title_fontsize=7, loc="lower left")

    ax = axes[1]
    sysd = b["systems"]
    xs = np.array([s["wdot_gr"] for s in sysd.values()])
    lim = [0, 0.00052]
    ax.plot(lim, lim, color=PS.INK, lw=0.9, label="GR")
    for name, s in sysd.items():
        used = s["used_in_test"]
        ax.errorbar(s["wdot_gr"], s["gr_measured_ours"], yerr=s["gr_measured_err"], fmt="o",
                    color=PS.ORANGE if used else PS.MUTED, ms=5, lw=1, mfc=None if used else "none")
        if used:
            off = {"KW Hya": (-36, 8), "EY Cep": (7, 2), "V541 Cyg": (7, -9)}.get(name, (6, -3))
            ax.annotate(name, (s["wdot_gr"], s["gr_measured_ours"]), off, textcoords="offset points",
                        fontsize=7, color=PS.INK_SECONDARY)
    A = b["combined_A"]
    ax.text(0.03, 0.97, f"slope A = (2 + 2 gamma - beta)/3\n= {A['A']:.3f} +- {A['err']:.3f}\n"
            f"(Baroch et al. 2021: {A['published']['A']:.3f} +- {A['published']['err']:.3f})",
            transform=ax.transAxes, va="top", fontsize=7.5, color=PS.INK_SECONDARY)
    ax.set_xlim(lim)
    ax.set_ylim(-0.00015, 0.00075)
    ax.set_xlabel("GR prediction (deg per orbit)")
    ax.set_ylabel("observed minus classical (deg per orbit)")
    ax.set_title("Binaries where GR does most of the turning", loc="left")

    ax = axes[2]
    d = PUBLISHED["di_her"]
    dh = res["di_her"]
    cpc = dh["cycles_per_century"]
    for k, m in enumerate(d["history"]):
        err = m.get("err")
        lo, hi = m.get("err_minus", err), m.get("err_plus", err)
        ax.errorbar(m["year"], m["obs"] * cpc, yerr=[[lo * cpc], [hi * cpc]], fmt="o", color=PS.INK,
                    ms=5, lw=1)
    ax.axhline(d["gr_claret10"] * cpc, color=PS.C_GR, lw=1.0)
    ax.text(1984, d["gr_claret10"] * cpc + 0.07, "general relativity alone", fontsize=7.5, color=PS.C_GR)
    ax.axhspan((d["theory_aligned_gm85"] - d["theory_aligned_gm85_err"]) * cpc,
               (d["theory_aligned_gm85"] + d["theory_aligned_gm85_err"]) * cpc, color=PS.GRID, lw=0)
    ax.text(1984, d["theory_aligned_gm85"] * cpc + 0.33, "GR + tides + spins, if spins were upright (1985)",
            fontsize=7.5, color=PS.INK_SECONDARY)
    lo = (d["theory_misaligned_claret21"] - d["theory_misaligned_claret21_err"]) * cpc
    hi = (d["theory_misaligned_claret21"] + d["theory_misaligned_claret21_err"]) * cpc
    ax.axhspan(lo, hi, color=PS.AQUA, alpha=0.25, lw=0)
    ax.text(1981, hi - 0.05, "... with the spins tipped\n~75-80 deg (Albrecht et al. 2009)",
            fontsize=7.5, color=PS.INK_SECONDARY, va="top")
    ax.set_xlim(1980, 2026)
    ax.set_ylim(0, 4.9)
    ax.set_xlabel("year of the measurement")
    ax.set_ylabel("apsidal motion (degrees per century)")
    ax.set_title("DI Herculis: the 'GR problem' that was spin", loc="left")
    PS.caption(fig, "Left: Baroch et al.'s tabulated eclipse-time differences (archival + TESS) refitted "
               "here with Sterne's formula. Middle: classical (tidal + rotational) terms from Baroch "
               "et al.'s stellar models; filled = the six systems of their GR test. Right: published "
               "DI Her rates (Guinan & Maloney 1985 to Liang et al. 2022) and predictions.")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "09_binaries")


# --------------------------------------------------------------------------
# 10. Numerical honesty
# --------------------------------------------------------------------------

def fig_numerics(res, ser):
    w = res["mercury"]["whfast"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    dts = np.array([v["dt_days"] for v in w["runs"].values()])
    err = np.array([abs(v["error"]) for v in w["runs"].values()])
    ax.loglog(dts, err, "o-", color=PS.BLUE, lw=1, ms=5, label="WHFast (fixed step), |rate - IAS15|")
    floor = abs(res["mercury"]["contributions"]["numerical_floor"])
    ax.axhline(max(floor, 1e-12), color=PS.ORANGE, lw=1, label=f"IAS15 floor: lone Mercury precesses {floor:.0e}\"/cy")
    ax.axhline(0.0015, color=PS.MUTED, lw=0.8)
    ax.text(0.55, 0.0019, "MESSENGER's quoted error, 0.0015\"/cy", fontsize=7.5, color=PS.INK_SECONDARY)
    ax.set_ylim(1e-10, 3.0)
    ax.annotate("step ~ 1/11 of Mercury's year:\na step-orbit resonance", (8, err[dts == 8][0]),
                (2.2, 0.12), fontsize=7.5, color=PS.INK_SECONDARY,
                arrowprops={"arrowstyle": "-", "color": PS.INK_SECONDARY, "lw": 0.6})
    ax.text(0.8, 4e-7, "error ~ step^2", fontsize=7.5, color=PS.INK_SECONDARY, rotation=17)
    ax.set_xlabel("time step (days); Mercury's year = 88 days")
    ax.set_ylabel("error in Mercury's precession (arcsec per century)")
    ax.set_title("How fine must the integration be?", loc="left")
    ax.legend(loc="lower right", fontsize=7.5, bbox_to_anchor=(1.0, 0.08))
    rev = res["mercury"]["reversibility"]
    ax.text(0.02, 0.30, f"Round trip 2000 -> 2650 -> 2000 (IAS15):\nMercury back within "
            f"{rev['position_error_km'] * 1000:.0f} m, perihelion within {rev['perihelion_direction_error_arcsec']:.0e}\"",
            transform=ax.transAxes, va="top", fontsize=7.5, color=PS.INK_SECONDARY)
    ms = res["mercury"].get("mass_scaling")
    if ms:
        f = np.array([m["factor"] for m in ms])
        r = np.array([m["rate_scaled_back"] for m in ms])
        g = res["mercury"]["analytic"]["gauss"]["jupiter"]
        ax2.semilogx(f, r, "o-", color=PS.BLUE, ms=6, lw=1, label="N-body, Jupiter's mass x f, rate / f")
        ax2.axhline(g, color=PS.ORANGE, lw=1, label=f"Gauss's first-order rate: {g:.3f}\"/cy")
        for fi, ri in zip(f, r):
            ax2.annotate(f"{ri:.3f}", (fi, ri), (0, 8), textcoords="offset points", ha="center",
                         fontsize=7.5, color=PS.INK_SECONDARY)
        ax2.set_xlabel("Jupiter's mass, as a fraction of the real one")
        ax2.set_ylabel("Jupiter's share, scaled back (arcsec per century)")
        ax2.set_title("The N-body sees Jupiter's second-order pull", loc="left")
        ax2.legend(loc="upper left", fontsize=7.5)
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "10_numerics")


# --------------------------------------------------------------------------
# 11. Six millennia
# --------------------------------------------------------------------------

def fig_millennia(res, ser):
    m = res["horizons"].get("millennia")
    if not m:
        return []
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.6), sharex=True)
    f, n = m["millennia_full"], m["millennia_newton"]
    ax = axes[0]
    ax.plot(n["years"], np.array(n["perihelion_offset_arcsec"]) / 60, "o-", ms=2.5, lw=0.8,
            color=PS.C_NEWTON, label="Newton only")
    ax.plot(f["years"], np.array(f["perihelion_offset_arcsec"]) / 60, "o-", ms=2.5, lw=0.8,
            color=PS.C_GR, label="with GR, J2 and the Sun's spin")
    ax.set_ylabel("perihelion minus\nDE441 (arcmin)")
    ax.legend(loc="upper right")
    ax.set_title("Six millennia: our integration against JPL's DE441, every 25 years", loc="left")
    ax = axes[1]
    ax.plot(f["years"], f["perihelion_offset_arcsec"], "o-", ms=2.5, lw=0.8, color=PS.C_GR)
    ax.axhline(0, color=PS.AXIS, lw=0.8)
    ax.set_ylabel("zoom (arcsec)")
    ast = res["mercury"]["asteroids"]["main_belt_rate"]
    ax.text(0.02, 0.08, f"drift {f['offset_rate_arcsec_per_cy']:+.4f}\"/cy: the asteroid belt, which DE441 "
            f"includes (estimate {-ast:+.4f})", transform=ax.transAxes, fontsize=7.5,
            color=PS.INK_SECONDARY)
    ax = axes[2]
    rv = m["rate_vs_epoch"]
    ax.plot(rv["year"], rv["rate"], color=PS.INK, lw=1.2)
    ax.set_ylabel("Mercury's precession\n(arcsec per century)")
    ax.set_xlabel("year")
    ax.text(0.02, 0.08, "in-plane rate in 300-year windows (our full model): the 575\"/cy is today's "
            "value, not a constant", transform=ax.transAxes, fontsize=7.5, color=PS.INK_SECONDARY)
    PS.caption(fig, "Simulation from DE440's J2000 state run to 1000 BC and AD 5000; DE441 positions "
               "from JPL Horizons at 241 epochs.")
    fig.tight_layout()
    return PS.save(fig, FIGURE_DIR / "11_six_millennia")


FIGURES = {"budget": fig_budget, "vs_de440": fig_vs_de440, "frames": fig_frames,
           "perihelia": fig_perihelia, "history": fig_history, "planets": fig_planets, "ladder": fig_ladder, "s2": fig_s2,
           "binaries": fig_binaries, "numerics": fig_numerics, "millennia": fig_millennia}


def main(names=None):
    res, ser = load()
    for name, fn in FIGURES.items():
        if names and name not in names:
            continue
        out = fn(res, ser)
        print(f"  {name}: {', '.join(str(p).split('/')[-1] for p in out) if out else 'skipped'}")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
