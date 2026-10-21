"""Every figure of the project, from cached results (see run.py)."""

from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch

import config as C
import literature as L
import plotstyle as ps

ps.use_style()
OUT = C.FIGURE_DIR


# ---------------------------------------------------------------------------
# Gravity Probe B
# ---------------------------------------------------------------------------

def gpb_result(pred: dict) -> list:
    """GP-B's four gyroscopes against general relativity's prediction (NS vs WE)."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.9), gridspec_kw={"width_ratios": [1.25, 1]})
    # --- left: the measurements in the drift-rate plane ------------------------------
    offsets = {"Gyro 1": (-4, 9, "right"), "Gyro 2": (7, 4, "left"), "Gyro 3": (8, 4, "left"),
               "Gyro 4": (-8, 4, "right"), "Joint": (9, -13, "left")}
    for name, (ns, ens, we, ewe) in L.GPB_RESULTS.items():
        joint = name == "Joint"
        ax.errorbar(we, ns, xerr=ewe, yerr=ens, fmt="D" if joint else "o",
                    color=ps.INK if joint else ps.MUTED, ms=7 if joint else 5,
                    mfc=ps.INK if joint else ps.SURFACE, capsize=0, lw=1.2 if joint else 0.8,
                    zorder=4 if joint else 3)
        dx, dy, ha = offsets[name]
        ax.annotate(name, (we, ns), xytext=(dx, dy), textcoords="offset points", fontsize=7.5,
                    ha=ha, color=ps.INK if joint else ps.INK_SECONDARY)
    jn = L.GPB_RESULTS["Joint"]
    ax.add_patch(Ellipse((jn[2], jn[0]), 2 * 2.4477 * jn[3], 2 * 2.4477 * jn[1], fill=False,
                         ec=ps.INK, lw=0.7, alpha=0.5))
    gr = L.GPB_GR_PREDICTION
    ax.plot(gr["we"], gr["ns"], **ps.ANCHOR_STYLE, zorder=6)
    ax.annotate("general relativity\n(GP-B team, from the\nGPS-tracked orbit)", (gr["we"], gr["ns"]),
                xytext=(-95, -6655), textcoords="data", fontsize=7.5, color=ps.INK,
                arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.7})
    ours_ns = pred["geodetic"]["ns_mas_yr"]
    ours_fd = pred["frame_dragging"]["we_mas_yr"]
    ours_we = pred["geodetic"]["we_mas_yr"] + ours_fd
    ax.plot(ours_we, ours_ns, "*", color=ps.BLUE, ms=13, mec=ps.SURFACE, mew=0.8, zorder=7)
    ax.annotate("our integration along\nthe published mean orbit", (ours_we, ours_ns),
                xytext=(-95, -6545), textcoords="data", fontsize=7.5, color=ps.BLUE,
                arrowprops={"arrowstyle": "-", "color": ps.BLUE, "lw": 0.7})
    ax.axvline(0, color=ps.AXIS, lw=0.8)
    ax.text(1.5, -6775, "no frame\ndragging", fontsize=7, color=ps.MUTED, ha="left", va="bottom")
    ax.set_xlim(-100, 22)
    ax.set_ylim(-6790, -6515)
    ax.set_xlabel("west-east drift (mas/yr); negative = eastward (frame dragging)")
    ax.set_ylabel("north-south drift (mas/yr); negative = southward,\nwithin the orbital plane (geodetic)")
    ax.set_title("Gravity Probe B: four gyroscopes, 2004-2005")
    ax.text(0.02, 0.03, "Error bars 1σ; ellipse 95% (joint). Everitt et al. 2011, PRL 106, 221101",
            transform=ax.transAxes, fontsize=7, color=ps.MUTED)
    # --- right: the prediction depends on the orbit's radius --------------------------
    scan = pred["radius_scan"]
    r = np.array([s["mean_r_km"] for s in scan])
    g = np.array([s["geodetic_ns"] for s in scan])
    gp = np.array([s["point_mass_ns"] for s in scan])
    o = np.argsort(r)
    ax2.plot(r[o], -g[o], color=ps.BLUE, lw=1.6, label="with the Earth's flattening (J2)")
    ax2.plot(r[o], -gp[o], color=ps.BLUE, lw=1.0, ls="--", label="point-mass Earth")
    ax2.axhline(6606.1, color=ps.INK, lw=1.0)
    ax2.text(7018.3, 6607.3, "GP-B team: 6,606.1 (from the GPS-tracked orbit)", fontsize=7.5, va="bottom")
    ax2.axhspan(6601.8 - 18.3, 6601.8 + 18.3, color=ps.WASH, zorder=0)
    ax2.text(7029.8, 6601.8 - 16.5, "measured ±1σ", fontsize=7, color=ps.MUTED, ha="right")
    for rk, lab in ((7018.0, "design mean\nsemi-major axis"), (7027.4, "'Quick facts'\nsemi-major axis")):
        ax2.axvline(rk, color=ps.AXIS, lw=0.8)
        ax2.text(rk + 0.15, 6633, lab, fontsize=7, color=ps.INK_SECONDARY, va="top")
    rp = pred["radius_for_published_km"]
    ax2.plot(rp, 6606.1, "o", color=ps.INK, ms=4)
    ax2.annotate(f"{rp:.1f} km", (rp, 6606.1), xytext=(4, -14), textcoords="offset points", fontsize=7.5)
    ax2.set_xlim(7017, 7031)
    ax2.set_ylim(6575, 6636)
    ax2.set_xlabel("mean orbital radius (km)")
    ax2.set_ylabel("predicted geodetic drift (mas/yr)")
    ax2.set_title("The prediction is set by the orbit's size")
    ax2.legend(loc="lower left", fontsize=7.5)
    ps.caption(fig, "Left: GP-B's measured drift of each gyroscope's spin axis relative to the guide star "
               "IM Pegasi, after removing the star's proper motion and the Sun's geodetic effect. "
               f"Our integration of de Sitter's and Schiff's formulas along the orbit gives "
               f"{-ours_ns:,.0f} (geodetic) and {-ours_fd:.1f} (frame dragging) mas/yr, plus "
               f"{-pred['geodetic']['we_mas_yr']:.1f} mas/yr of geodetic turning that leaks west-east as the "
               "orbital plane wanders ~100\" from the star. Right: the geodetic rate scales as "
               f"r^-5/2, {abs(pred['sensitivity_mas_yr_per_km']):.1f} mas/yr per km; the Earth's flattening "
               "lowers it by ~4.4 mas/yr. The published orbit elements are ambiguous at the km level.",
               y=-0.035)
    fig.tight_layout()
    return ps.save(fig, OUT / "gpb_result")


def gpb_diagram() -> list:
    """What the two effects do to a gyroscope in GP-B's polar orbit (not to scale).

    Left: viewed from local west at the guide star (north up, IM Peg to the right),
    GP-B's orbit runs clockwise and geodetic precession, prograde, tips the spin
    axis south. Right: viewed from the guide star (north up, east to the right),
    frame dragging pushes the spin axis east over the poles (and west, half as
    strongly, over the equator; eastward on average).
    """
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 5.6))
    cy, R, RE = -0.35, 1.45, 0.6
    for ax in (a1, a2):
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_xlim(-2.7, 3.0)
        ax.set_ylim(-2.75, 2.95)
        ax.add_patch(Circle((0, cy), RE, color=ps.BLUE_150, zorder=1))
        ax.plot([0, 0], [cy - RE - 0.22, cy + RE + 0.22], color=ps.INK, lw=1.1, zorder=2)
        ax.text(0.06, cy + RE + 0.24, "N", fontsize=8.5, color=ps.INK, va="bottom")
    # --- left: face-on to the orbital plane ------------------------------------------------
    a1.add_patch(Circle((0, cy), R, fill=False, ec=ps.AXIS, lw=1.2, zorder=1))
    for a0, a1_ in ((60, 30), (240, 210)):        # clockwise arrows on the orbit
        a1.annotate("", xy=(R * np.cos(np.radians(a1_)), cy + R * np.sin(np.radians(a1_))),
                    xytext=(R * np.cos(np.radians(a0)), cy + R * np.sin(np.radians(a0))),
                    arrowprops={"arrowstyle": "->", "color": ps.MUTED, "lw": 1.1,
                                "connectionstyle": "arc3,rad=-0.12"})
    a1.text(1.12, cy + 1.22, "GP-B's polar orbit,\n97.6 min", fontsize=8, color=ps.MUTED)
    tilt = np.radians(10)           # exaggerated about 5 million times
    for ang in np.radians([90, 180, 270, 0]):
        p0 = np.array([R * np.cos(ang), cy + R * np.sin(ang)])
        a1.plot(*p0, "o", color=ps.INK, ms=5, zorder=4)
        a1.annotate("", xy=p0 + [0.65, 0.0], xytext=p0,
                    arrowprops={"arrowstyle": "-|>", "color": ps.MUTED, "lw": 1.0, "ls": "--"})
        a1.annotate("", xy=p0 + 0.65 * np.array([np.cos(-tilt), np.sin(-tilt)]), xytext=p0,
                    arrowprops={"arrowstyle": "-|>", "color": ps.GEODETIC, "lw": 1.9})
    a1.annotate("", xy=(2.95, -2.45), xytext=(1.85, -2.45),
                arrowprops={"arrowstyle": "-|>", "color": ps.INK, "lw": 1.0})
    a1.text(1.8, -2.45, "to the guide star IM Pegasi", fontsize=8.5, ha="right", va="center")
    a1.text(-2.65, 2.85, "Geodetic precession: curved space", fontsize=10.5, fontweight="bold",
            va="top")
    a1.text(-2.65, 2.45, "Face-on to the orbit. Each gyroscope's spin axis (grey dashed: at the start)\n"
                         "turns within the orbital plane, in the sense of the orbit, so it\n"
                         "drifts south of the star: 6,606 mas/yr (blue; turning exaggerated)",
            fontsize=8, color=ps.INK_SECONDARY, va="top")
    # --- right: from the guide star, the orbit edge-on ------------------------------------------
    a2.plot([0, 0], [cy - R, cy + R], color=ps.AXIS, lw=1.2, zorder=0)
    a2.annotate("", xy=(0.3, cy + RE + 0.1), xytext=(-0.3, cy + RE + 0.1),
                arrowprops={"arrowstyle": "->", "color": ps.INK, "lw": 1.0,
                            "connectionstyle": "arc3,rad=0.45"})
    a2.text(0.38, cy + RE + 0.14, "Earth's rotation\n(near side moves east)", fontsize=7.5, va="center")
    for y, lab in ((cy + R, "over the North Pole"), (cy - R, "over the South Pole")):
        a2.plot(0, y, "o", ms=12, mfc=ps.SURFACE, mec=ps.INK, mew=1.0, zorder=3)
        a2.plot(0, y, "o", ms=3, color=ps.INK, zorder=4)
        a2.annotate("", xy=(0.75, y), xytext=(0.12, y),
                    arrowprops={"arrowstyle": "-|>", "color": ps.FRAME, "lw": 1.9})
        a2.text(-0.2, y, lab, fontsize=7.5, color=ps.INK_SECONDARY, ha="right", va="center")
    a2.text(0.82, cy + R, "east", fontsize=8, color=ps.FRAME, va="center")
    a2.text(0.82, cy - R, "east", fontsize=8, color=ps.FRAME, va="center")
    a2.text(-0.12, cy - R - 0.42, "orbit edge-on: it contains the line of sight to the star",
            fontsize=7.5, color=ps.MUTED, ha="center")
    a2.text(-2.65, 2.85, "Frame dragging: the spinning Earth", fontsize=10.5, fontweight="bold", va="top")
    a2.text(-2.65, 2.45, "Seen from the guide star. The spin axis (pointing at you, ⊙) is dragged\n"
                         "east, across the orbital plane, strongest over the poles (it is pushed\n"
                         "west, half as hard, over the equator): 39 mas/yr on average",
            fontsize=8, color=ps.INK_SECONDARY, va="top")
    ps.caption(fig, "Gravity Probe B's gyroscopes pointed at IM Pegasi, a star in the plane of the polar "
               "orbit. Curved space around the Earth turns them within the plane (north-south on the sky); "
               "the Earth's rotation drags them across it (west-east). The effects are at right angles, "
               "so one experiment separates them. Not to scale.", y=0.02)
    return ps.save(fig, OUT / "gpb_diagram")


# ---------------------------------------------------------------------------
# LAGEOS and LARES 2: the bulge-driven drift and frame dragging
# ---------------------------------------------------------------------------

def node_drift(we: dict, rates: dict) -> list:
    """The two orbital planes swing round in opposite directions; their sum barely moves.

    we: {sat: weekly_elements dict}; rates: analytic/fitted drift rates for the labels.
    """
    import frames as F
    fig = plt.figure(figsize=(11.5, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.25)
    ax = fig.add_subplot(gs[0, :])
    for sat in ("lageos1", "lares2"):
        w = we[sat]
        yr = F.mjd_to_year(w["mjd"])
        node = np.degrees(w["node_cip"]) % 360
        st = ps.SAT_STYLE[sat]
        ax.plot(yr, node, st["marker"], color=st["color"], ms=2.6, label=st["label"] + f": {rates[sat]:+.2f}°/yr")
    ax.set_ylim(0, 360)
    ax.set_yticks([0, 90, 180, 270, 360])
    ax.set_ylabel("orbital node (degrees,\nalong the true equator)")
    ax.set_title("Two orbital planes, turned in opposite directions by the Earth's equatorial bulge")
    ax.legend(loc="upper center", ncol=2, fontsize=8, bbox_to_anchor=(0.5, 1.0), frameon=True,
              facecolor=ps.SURFACE, edgecolor=ps.GRID)
    ax.set_xlim(2022.5, 2026.85)
    # --- bottom left: the sum ------------------------------------------------------------
    ax2 = fig.add_subplot(gs[1, 0])
    m1, m2 = we["lageos1"]["mjd"], we["lares2"]["mjd"]
    common = np.intersect1d(np.round(m1, 1), np.round(m2, 1))
    i1 = np.searchsorted(np.round(m1, 1), common)
    i2 = np.searchsorted(np.round(m2, 1), common)
    n1 = np.unwrap(we["lageos1"]["node_cip"])[i1]
    n2 = np.unwrap(we["lares2"]["node_cip"])[i2]
    ssum = np.degrees(n1 + n2)
    yr = F.mjd_to_year(m1[i1])
    ax2.plot(yr, ssum - ssum[0], "D", color=ps.INK, ms=2.4)
    slope = np.polyfit(yr, ssum, 1)[0]
    ax2.set_ylabel("sum of the two nodes (deg)")
    ax2.set_title(f"Their sum drifts only {slope:+.3f}°/yr", fontsize=10)
    ax2.text(0.03, 0.06, f"{abs(slope) / abs(rates['lageos1']) * 100:.2f}% of each drift: the inclinations sum to "
                         f"{rates['inc_sum']:.3f}°,\nnot 180°, and the orbits' sizes differ by {rates['da_km']:.2f} km",
             transform=ax2.transAxes, fontsize=7.5, color=ps.INK_SECONDARY)
    ax2.set_xlim(2022.5, 2026.85)
    # --- bottom right: scale ladder ------------------------------------------------------
    ax3 = fig.add_subplot(gs[1, 1])
    labels = ["bulge drift,\neach satellite", "residual of\nthe sum", "frame dragging\nin the sum",
              "de Sitter\nin the sum"]
    vals = [abs(rates["lageos1"]) * 3.6e6, abs(slope) * 3.6e6, rates["lt_sum"], rates["ds_sum"]]
    cols = [ps.CONTEXT, ps.INK, ps.FRAME, ps.GEODETIC]
    ax3.barh(range(len(vals))[::-1], vals, color=cols, height=0.6)
    ax3.set_xscale("log")
    ax3.set_yticks(range(len(vals))[::-1])
    ax3.set_yticklabels(labels, fontsize=8)
    for k, v in enumerate(vals):
        ax3.text(v * 1.4, len(vals) - 1 - k, f"{v:,.0f} mas/yr" if v > 100 else f"{v:.1f} mas/yr",
                 va="center", fontsize=8)
    ax3.set_xlim(10, 3e10)
    ax3.set_xlabel("node drift (mas/yr, log scale)")
    ax3.set_title("Seven orders of magnitude", fontsize=10)
    ax3.grid(axis="y", visible=False)
    ps.caption(fig, "Weekly mean nodes of LAGEOS and LARES 2 from the ILRS analysis centre ASI's orbits, "
               "measured from the Earth's true equator (CIP/CIO frame). The Earth's bulge (J2) turns each "
               "orbital plane by ~125° a year, in opposite directions because the inclinations are "
               "supplementary (109.84° and 70.16°). Frame dragging turns both the same way, by 31 mas/yr.",
               y=0.0)
    return ps.save(fig, OUT / "node_drift")


AC_NAMES = {"asi": "ASI (Italy)", "dgfi": "DGFI-TUM (Germany)", "jcet": "JCET/GSFC (USA)",
            "gfz": "GFZ (Germany)", "esa": "ESA/ESOC", "bkg": "BKG (Germany)", "nsgf": "NSGF (UK)",
            "cnes": "CNES/GRGS (France)"}


def frame_dragging(c: dict, summary: dict | None = None, ac_label: str = "ASI") -> list:
    """The measurement: summed node residuals of LAGEOS + LARES 2 against GR's frame dragging."""
    import frames as F
    yr = F.mjd_to_year(c["mjd"])
    fig = plt.figure(figsize=(11.5, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1], hspace=0.38, wspace=0.24)
    ax = fig.add_subplot(gs[0, :])
    t = c["t_yr"]
    s = c["sum"]
    lt = c["lt_sum"]
    ax.fill_between(yr, (c["mu"] - c["mu_err"]) * lt * t, (c["mu"] + c["mu_err"]) * lt * t,
                    color=ps.BLUE_150, alpha=0.6, lw=0, zorder=1)
    ax.plot(yr, s["raw"], "o", color=ps.CONTEXT, ms=2.6, zorder=2,
            label="as measured (summed weekly jumps)")
    ax.plot(yr, s["clean"], "D", color=ps.INK, ms=2.8, zorder=4,
            label="after removing the fitted tide-model errors")
    ax.plot(yr, lt * t, color=ps.FRAME, lw=2.0, zorder=3,
            label=f"general relativity: frame dragging, {lt:.1f} mas/yr")
    ax.plot(yr, c["mu"] * lt * t, color=ps.BLUE, lw=1.2, zorder=3,
            label=f"our fit ({ac_label}'s orbits): μ = {c['mu']:.3f} ± {c['mu_err']:.3f} (stat.)")
    ax.axhline(0, color=ps.MUTED, lw=1.0, ls="--")
    ax.text(yr[-1], 6, "no frame dragging", ha="right", va="bottom", fontsize=8, color=ps.MUTED)
    ax.set_ylabel("node residual, LAGEOS + LARES 2 (mas)")
    ax.set_title("Frame dragging, measured from four years of public laser-ranging orbits")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_xlim(yr[0] - 0.05, yr[-1] + 0.05)
    # --- every analysis centre ------------------------------------------------------------
    ax2 = fig.add_subplot(gs[1, 0])
    if summary is not None:
        acs = list(summary["ac_mus"])
        lo_x, hi_x = 0.9, 1.1
        for k, ac in enumerate(acs):
            mu, e = summary["ac_mus"][ac], summary["ac_errs"][ac]
            col = ps.BLUE if ac == "asi" else ps.INK_SECONDARY
            if e > 0.1:            # imprecise centre: say so instead of clipping it silently
                ax2.annotate(f"{mu:.2f} ± {e:.2f} (weekly frame rotations)", (lo_x + 0.003, k),
                             fontsize=7, color=ps.MUTED, va="center")
                continue
            ax2.errorbar(mu, k, xerr=e, fmt="D" if ac == "asi" else "o", color=col, ms=5,
                         mfc=col if ac == "asi" else ps.SURFACE, capsize=0, lw=1.0)
        ax2.set_yticks(range(len(acs)))
        ax2.set_yticklabels([AC_NAMES.get(a, a) for a in acs], fontsize=7.5)
        ax2.invert_yaxis()
        ax2.axvline(1.0, color=ps.FRAME, lw=1.6)
        ax2.axvspan(summary["mu"] - summary["total_err"], summary["mu"] + summary["total_err"],
                    color=ps.BLUE_150, alpha=0.5, lw=0)
        ax2.set_xlim(lo_x, hi_x)
        ax2.set_xlabel("μ = measured / general relativity")
        ax2.set_title("Eight analysis centres, eight sets of orbits", fontsize=10)
        ax2.grid(axis="y", visible=False)
    # --- the tidal terms we removed --------------------------------------------------------
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(yr, s["nuisance"], color=ps.INK, lw=1.4, label="sum")
    for sat in ("lageos1", "lares2"):
        st = ps.SAT_STYLE[sat]
        ax3.plot(yr, c[sat]["nuisance"] - c[sat]["nuisance"][0], color=st["color"], lw=0.9,
                 alpha=0.85, label=st["label"].split(" (")[0])
    ax3.axhline(0, color=ps.MUTED, lw=0.8)
    ax3.set_ylabel("fitted periodic + J2 terms (mas)")
    ax3.set_title("What the fit removed: tide-model errors", fontsize=10)
    lo, hi = ax3.get_ylim()
    ax3.set_ylim(lo, hi + 0.35 * (hi - lo))
    ax3.legend(loc="upper left", fontsize=7.5, ncol=3)
    ax3.set_xlim(yr[0] - 0.05, yr[-1] + 0.05)
    ps.caption(fig, f"Every week each ILRS analysis centre fits a 7-day orbit to each satellite's laser "
               "ranges. We refitted every week with our own orbit model (frame dragging switched off, the "
               "radiation-pressure scale held fixed) and measured how far the next week's orbital plane had "
               "turned beyond our prediction. Summed over the two satellites, the Earth's bulge and errors "
               "in J2 cancel; what accumulates is frame dragging. The fitted tide-model terms (K1, K2, S2, "
               "P1, radiation pressure, annual) are removed in black. Bands: ±1σ (top: statistical; "
               "bottom left: total).", y=0.0)
    return ps.save(fig, OUT / "frame_dragging")


def systematics(summary: dict, rows: list) -> list:
    """Forest plot: mu for every analysis choice, centre and model variation."""
    fig, ax = plt.subplots(figsize=(8.6, 0.32 * len(rows) + 1.6))
    y = 0
    labels = []
    headers = []
    for r in rows:
        if r.get("header"):
            labels.append(r["label"])
            headers.append(y)
            y += 1
            continue
        col = ps.BLUE if r.get("primary") else (ps.MUTED if r.get("test") else ps.INK_SECONDARY)
        ax.errorbar(r["mu"], y, xerr=r["err"], fmt="D" if r.get("primary") else "o", color=col, ms=5,
                    mfc=col if r.get("primary") else ps.SURFACE, capsize=0, lw=1.0)
        if r.get("note"):
            ax.text(r["mu"] + r["err"] + 0.012, y, r["note"], fontsize=6.8, va="center", color=ps.MUTED)
        labels.append(r["label"])
        y += 1
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.5)
    for k, tick in enumerate(ax.get_yticklabels()):
        if k in headers:
            tick.set_fontweight("bold")
            tick.set_fontsize(8.5)
            tick.set_color(ps.INK)
    ax.invert_yaxis()
    ax.axvline(1.0, color=ps.FRAME, lw=1.6, label="general relativity")
    ax.axvspan(summary["mu"] - summary["total_err"], summary["mu"] + summary["total_err"],
               color=ps.BLUE_150, alpha=0.45, lw=0, label="our result ±1σ (total)")
    ax.set_xlim(0.55, 1.45)
    ax.set_xlabel("μ = measured frame dragging / general relativity")
    ax.set_title("How the result depends on the analysis choices")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=7.5)
    fig.tight_layout()
    return ps.save(fig, OUT / "systematics")


def j2_figure(j: dict, nodal: dict) -> list:
    """The Earth's flattening as LAGEOS + LARES 2 measure it, against CSR's laser-ranging C20.

    j: framedrag.j2_comparison(); nodal: {'year','solid','ocean_iers','ocean_flip','cosN'} on a
    long grid, plus 'corrected_offset' (pair - CSR after the in-phase ocean tide).
    """
    ref = -4.84169e-4
    sc = 1e10
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
    ax.errorbar(j["csr_year"], (j["csr"] - ref) * sc, yerr=j["csr_err"] * sc, fmt="o", ms=3.2,
                color=ps.MUTED, mfc=ps.SURFACE, lw=0.7, capsize=0, label="CSR monthly C20 from SLR (5 satellites)")
    ax.plot(j["year"], (j["pair_smooth"] - ref) * sc, color=ps.BLUE, lw=1.6,
            label="LAGEOS - LARES 2 (this project, IERS tide file as published)")
    ax.plot(j["year"], (j["pair_smooth"] - nodal["fix_on_pair"] - ref) * sc, color=ps.ORANGE, lw=1.4,
            label="... with the 18.6-yr ocean tide put in phase")
    ax.plot(j["year"], (j["goco"] - ref) * sc, color=ps.INK, lw=1.0, ls="--", label="GOCO06s model (trend + annual, extrapolated)")
    ax.set_ylabel("C̄20 + 4.84169×10⁻⁴   (×10⁻¹⁰)")
    ax.set_title("The Earth's flattening (J2), seen by the same two satellites")
    ax.set_ylim(-12.2, -3.4)
    ax.legend(loc="lower center", fontsize=7.2, ncol=2, columnspacing=1.0)
    ax.text(0.99, 0.985, f"offset from CSR: {j['offset'] * sc:+.2f} ± {j['offset_err'] * sc:.2f}; "
                         f"after the phase fix: {nodal['corrected_offset'] * sc:+.2f} ± "
                         f"{nodal['corrected_err'] * sc:.2f}; correlation {j['corr']:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.3, color=ps.INK_SECONDARY,
            bbox={"facecolor": ps.SURFACE, "edgecolor": "none", "pad": 1.5})
    # --- the 18.6-year tide ---------------------------------------------------------------------
    ax2.plot(nodal["year"], nodal["solid"] * sc, color=ps.INK, lw=1.4, label="solid Earth (IERS steps 1-2)")
    ax2.plot(nodal["year"], nodal["ocean_iers"] * sc * 3, color=ps.BLUE, lw=1.3,
             label="ocean, IERS FES2004 file (×3)")
    ax2.plot(nodal["year"], nodal["ocean_flip"] * sc * 3, color=ps.ORANGE, lw=1.3, ls="--",
             label="ocean, in phase with the solid tide (×3)")
    ax2.axvspan(j["year"][0], j["year"][-1], color=ps.WASH, zorder=0)
    ax2.text(j["year"][0] + 0.1, -4.6, "our data", fontsize=7.5, color=ps.MUTED)
    ax2.axhline(0, color=ps.AXIS, lw=0.8)
    ax2.set_ylabel("18.6-year term in C̄20 (×10⁻¹⁰)")
    ax2.set_title("The lunar nodal tide", fontsize=10)
    ax2.legend(loc="upper right", fontsize=7)
    ax2.set_ylim(-5, 5.5)
    ps.caption(fig, "Left: frame dragging moves the two orbital planes alike, but an error in J2 moves them "
               "oppositely, so the difference of the weekly node jumps measures J2 (smoothed over five weeks). "
               "Right: when the Moon's orbit tilts furthest from the equator (major lunar standstill, 2025) "
               "the tidal bulge weakens; the IERS file's equilibrium ocean tide moves the other way. Putting it in "
               f"phase removes the {j['offset'] / j['offset_err']:.0f}σ offset from CSR's series.", y=-0.03)
    fig.tight_layout()
    return ps.save(fig, OUT / "j2_from_the_pair")


def ladder(rows: list, ours: dict | None = None) -> list:
    """Measured / predicted for spin and orbit precession, from the Moon to pulsars."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 5.0), gridspec_kw={"width_ratios": [1.2, 1]})
    # label positions (axes data coordinates) for each measurement, chosen by hand
    geo_pos = {"Moon's orbit (LLR)": (5.3e-6, 0.84), "GP-B geodetic": (1.8e-3, 0.84),
               "Stairs et al. 2004": (0.03, 0.62), "Fonseca et al. 2014": (0.025, 1.42),
               "Desvignes et al. 2019": (0.9, 0.70), "Breton et al. 2008": (9.5, 0.80),
               "Lower et al. 2024": (9.5, 1.30)}
    geo = [r for r in rows if r["kind"] == "geodetic"]
    for k, r in enumerate(sorted(geo, key=lambda r: (r["rate"], r["year"]))):
        x = r["rate"] * (1.18 if r["ref"] in ("Fonseca et al. 2014", "Lower et al. 2024") else 1.0)
        filled = "GP-B" in r["label"] or "LLR" in r["label"]
        a1.errorbar(x, r["ratio"], yerr=[[r["err_lo"]], [r["err_hi"]]], fmt="o", ms=5,
                    color=ps.GEODETIC, mfc=ps.GEODETIC if filled else ps.SURFACE, capsize=0, lw=1.0)
        key = r["label"] if r["label"] in geo_pos else r["ref"]
        tx, ty = geo_pos.get(key, (x, r["ratio"] + 0.1))
        a1.annotate(f"{r['label']}\n{r['ref']}", (x, r["ratio"]), xytext=(tx, ty), textcoords="data",
                    fontsize=6.8, ha="center", va="center", color=ps.INK_SECONDARY,
                    arrowprops={"arrowstyle": "-", "color": ps.AXIS, "lw": 0.6, "shrinkA": 6, "shrinkB": 3})
    a1.axhline(1.0, color=ps.INK, lw=1.0)
    a1.set_xscale("log")
    a1.set_xlim(2e-6, 40)
    a1.set_ylim(0.5, 1.9)
    a1.set_xlabel("general relativity's precession rate (degrees per year)")
    a1.set_ylabel("measured / general relativity")
    a1.set_title("Geodetic precession: the Moon, a gyroscope, pulsars")
    for x, lab in ((19.2 / 3.6e6, "19 mas/yr"), (6606 / 3.6e6, "6.6\"/yr"), (5.07, "5°/yr")):
        a1.text(x, 0.52, lab, fontsize=7, color=ps.MUTED, ha="center")
    # --- frame dragging ------------------------------------------------------------------------
    fd_pos = {2004: (2005.0, 0.80), 2011: (2011.0, 0.68), 2016: (2014.6, 1.13), 2019: (2018.2, 1.20),
              2020: (2020.6, 0.83), 2026: (2024.0, 1.12)}
    fd = [r for r in rows if r["kind"] == "frame"]
    for r in sorted(fd, key=lambda r: r["year"]):
        is_ours = r.get("ours", False)
        col = ps.BLUE if is_ours else ps.FRAME
        x = r["year"] + (0.45 if is_ours else 0)
        a2.errorbar(x, r["ratio"], yerr=[[r["err_lo"]], [r["err_hi"]]], fmt="D" if is_ours else "o",
                    ms=6 if is_ours else 5, color=col, capsize=0, lw=1.3 if is_ours else 1.0,
                    mfc=col if is_ours or "GP-B" in r["label"] else ps.SURFACE)
        tx, ty = (2026.5, 0.80) if is_ours else fd_pos[r["year"]]
        a2.annotate(r["label"].replace(" (", "\n(") if not is_ours else "this project:\nLAGEOS + LARES 2",
                    (x, r["ratio"]), xytext=(tx, ty), textcoords="data", fontsize=6.6, ha="center",
                    va="center", color=ps.BLUE if is_ours else ps.INK_SECONDARY,
                    arrowprops={"arrowstyle": "-", "color": ps.AXIS, "lw": 0.6, "shrinkA": 8, "shrinkB": 3})
    a2.axhline(1.0, color=ps.INK, lw=1.0)
    a2.set_xlim(2002.0, 2029.5)
    a2.set_ylim(0.6, 1.25)
    a2.set_xlabel("year of publication")
    a2.set_ylabel("measured / general relativity")
    a2.set_title("Frame dragging by the spinning Earth")
    ps.caption(fig, "Left: spin (geodetic) precession, measured against GR's prediction from each system's "
               "masses and orbit (Barker-O'Connell; our own computation for the pulsars): the Moon's orbit "
               "round the Sun (lunar laser ranging), Gravity Probe B, and pulsars whose beams change as they "
               "precess. Right: frame dragging, from GP-B's gyroscopes and the orbital planes of laser-ranged "
               "satellites. Error bars are total 1σ (statistical and systematic in quadrature).", y=-0.03)
    fig.tight_layout()
    return ps.save(fig, OUT / "precession_ladder")

