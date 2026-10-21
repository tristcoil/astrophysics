"""Every static figure (PNG + SVG), plus the gallery animation and the interactive page."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import config as C
import literature as L
import orbits
from plotstyle import (AQUA, AXIS, BLUE, BLUE_150, BLUE_300, BLUE_550, CONTEXT, GRID, INK, INK_SECONDARY, MUTED,
                       ORANGE, SURFACE, WASH, YELLOW, annotate, caption, save, use_style)

use_style()
FIG = C.FIGURE_DIR

METHOD_STYLE = {
    "M": {"color": CONTEXT, "marker": ".", "ms": 3.5, "label": "visual micrometer (1862-1980)", "zorder": 2},
    "P": {"color": BLUE, "marker": "o", "ms": 3.2, "label": "photographic plates (1930s-1986)", "zorder": 3},
    "C": {"color": AQUA, "marker": "D", "ms": 4.0, "label": "CCD and mid-infrared", "zorder": 4},
    "H": {"color": ORANGE, "marker": "s", "ms": 4.2, "label": "Hubble (1997-2016)", "zorder": 5},
}


def _east_north(rho, pa_deg):
    t = np.radians(pa_deg)
    return rho * np.sin(t), rho * np.cos(t)


def _star(ax, x=0.0, y=0.0, size=16, label=None):
    ax.plot([x], [y], marker="*", ms=size, color=INK, mec=SURFACE, mew=1.0, zorder=8, linestyle="none",
            label=label)


# ---------------------------------------------------------------------------
# 1. Sirius: the orbit and 154 years of data
# ---------------------------------------------------------------------------


def fig_sirius_orbit(res):
    fit = res["sirius"]
    df = fit["df"]
    el = np.median(fit["hier"]["elements"], axis=0)
    use = df["ours_used"]
    flagged = fit["screen"]["flagged"]
    fig = plt.figure(figsize=(12.5, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.0], hspace=0.12, wspace=0.18)
    ax = fig.add_subplot(gs[:, 0])
    for m in ("M", "P", "C", "H"):
        st = METHOD_STYLE[m]
        d = df[use & ~flagged & (df["method"] == m) & df["has_pa"] & df["has_sep"]]
        e, n = _east_north(d["sep"], d["pa"])
        ax.plot(e, n, linestyle="none", marker=st["marker"], ms=st["ms"], color=st["color"],
                mec=SURFACE if m != "M" else st["color"], mew=0.4, label=f"{st['label']}  ({len(d)})",
                zorder=st["zorder"], alpha=0.85 if m == "M" else 1.0, rasterized=(m == "M"))
    d = df[flagged & df["has_pa"] & df["has_sep"]]
    e, n = _east_north(d["sep"], d["pa"])
    ax.plot(e, n, linestyle="none", marker="x", ms=3, color=MUTED, mew=0.6, alpha=0.6,
            label=f"flagged as outliers ({int(flagged.sum())})", zorder=1)
    t = el[1] + el[0] * np.linspace(0, 1, 1441)
    xn, ye = orbits.relative_xy(t, *el)
    ax.plot(ye, xn, color=INK, lw=1.4, zorder=6, label="our orbit")
    _star(ax, label="Sirius A")
    # one cycle of years along the orbit (positions repeat every P), labels pushed outward
    cx, cy = np.mean(ye), np.mean(xn)
    for yr, lab in ((el[1], f"{el[1]:.1f}\nperiastron"), (2000.0, "2000"), (2005.0, "2005"), (2010.0, "2010"),
                    (el[1] + el[0] / 2, f"{el[1] + el[0] / 2:.1f}\napastron"),
                    (2030.0, "2030"), (2037.0, "2037")):
        xn1, ye1 = orbits.relative_xy(yr, *el)
        ax.plot([ye1], [xn1], marker="o", ms=4.5, color=SURFACE, mec=INK, mew=1.0, zorder=7)
        dx, dy = ye1 - cx, xn1 - cy
        nrm = np.hypot(dx, dy)
        # data x runs right-to-left on screen: flip the x offset
        ax.annotate(lab, (ye1, xn1), xytext=(-18 * dx / nrm, 16 * dy / nrm), textcoords="offset points",
                    fontsize=7.5, color=INK_SECONDARY, ha="center", va="center", linespacing=0.95)
    ax.text(0.98, 0.02, "positions repeat every %.2f years:\nB was discovered in 1862 near\nits 2012 position"
            % el[0], transform=ax.transAxes, fontsize=7.5, color=MUTED, ha="right", va="bottom")
    ax.set_xlim(12.5, -9.5)   # east to the left, as on the sky
    ax.set_ylim(-9.0, 13.0)
    ax.set_aspect("equal")
    ax.set_xlabel("Δα cos δ  (arcsec, east ←)")
    ax.set_ylabel("Δδ  (arcsec, north ↑)")
    ax.set_title("Sirius B around Sirius A, 1862-2016")
    ax.legend(loc="lower left", fontsize=7.5, handletextpad=0.4, borderaxespad=0.2)
    # right: separation and position angle against time
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 1], sharex=ax1)
    tt = np.linspace(1855, 2035, 3000)
    rho, pa = orbits.rho_theta(tt, *el)
    for m in ("M", "P", "C", "H"):
        st = METHOD_STYLE[m]
        d = df[use & ~flagged & (df["method"] == m)]
        ax1.plot(d["year"], d["sep"], linestyle="none", marker=st["marker"], ms=st["ms"] * 0.8,
                 color=st["color"], mec=SURFACE if m != "M" else st["color"], mew=0.3, zorder=st["zorder"])
        ax2.plot(d["year"], d["pa"], linestyle="none",
                 marker=st["marker"], ms=st["ms"] * 0.8, color=st["color"],
                 mec=SURFACE if m != "M" else st["color"], mew=0.3, zorder=st["zorder"])
    ax1.plot(tt, rho, color=INK, lw=1.0, zorder=1)
    jump = np.where(np.abs(np.diff(pa)) > 180)[0]
    pa_plot = pa.copy()
    pa_plot[jump + 1] = np.nan
    ax2.plot(tt, pa_plot, color=INK, lw=0.9, zorder=1)
    ax1.set_ylabel("separation (arcsec)")
    ax2.set_ylabel("position angle (deg)")
    ax2.set_xlabel("year")
    ax1.set_title("Three orbits in the data: P = %.2f yr" % el[0])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.set_ylim(0, 360)
    ax2.set_yticks([0, 90, 180, 270, 360])
    for ax_ in (ax1, ax2):
        for tp in (1894.6, 1944.7, 1994.6):
            ax_.axvline(tp, color=GRID, lw=0.8, zorder=0)
    ax1.text(1894.6, 12.3, "periastron", fontsize=7, color=MUTED, ha="center")
    ax1.set_ylim(0, 13)
    caption(fig, "Data: Bond et al. 2017 (ApJ 840, 70; VizieR J/ApJ/840/70), all published measures corrected to "
                 "J2000. Orbit: our hierarchical fit (each observer's personal equation marginalised), MCMC median.",
            y=0.0)
    return save(fig, FIG / "sirius_orbit")


# ---------------------------------------------------------------------------
# 2. Sirius: residuals, personal equations, fits compared
# ---------------------------------------------------------------------------


def fig_sirius_systematics(res):
    fit = res["sirius"]
    hl = fit["hier"]["like"]
    p_med = np.median(fit["hier"]["chain"], axis=0)
    el = hl.unpack(p_med)[0]
    df = fit["df"]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1.0, 1.1, 1.0], "wspace": 0.32})
    # (a) HST residuals in mas
    ax = axs[0]
    h = df[df["method"].eq("H") & df["ours_used"]]
    rho, pa = orbits.rho_theta(h["year"].to_numpy(), *el)
    dsep = (h["sep"] - rho) * 1e3
    dtan = h["sep"] * np.radians(orbits.wrap_angle_deg(h["pa"] - pa)) * 1e3
    s_h = np.exp(p_med[7])
    ax.errorbar(h["year"], dsep, yerr=h["e_sep"] * 1e3 * s_h, fmt="s", ms=4.5, color=ORANGE, mec=SURFACE,
                ecolor=ORANGE, elinewidth=0.9, capsize=0, label="along the line A→B (separation)")
    ax.errorbar(h["year"] + 0.15, dtan, yerr=h["sep"] * np.radians(h["e_pa"]) * 1e3 * s_h, fmt="o", ms=4.5,
                color=BLUE, mec=SURFACE, ecolor=BLUE, elinewidth=0.9, capsize=0, label="across it (position angle)")
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xlabel("year")
    ax.set_ylabel("Hubble measure - orbit (mas)")
    ax.set_title("Hubble residuals: a few mas")
    ax.set_ylim(-24, 22)
    ax.legend(loc="lower center", fontsize=7.5, ncol=2, columnspacing=1.0, handletextpad=0.3)
    ax.text(0.02, 0.97, f"errors scaled by {s_h:.2f} (fitted)", transform=ax.transAxes, fontsize=7.5,
            color=MUTED, va="top")
    # (b) personal equations of the most prolific observers
    ax = axs[1]
    off = fit["hier"]["offsets"]["pa"]
    off = off[off["group"].str.startswith("M:")].copy()
    off["name"] = off["group"].str[2:]
    years = df[df["method"].eq("M")].assign(name=df["obs"].map(lambda s: __import__("sirius").person(s))) \
        .groupby("name")["year"].median()
    off["year"] = off["name"].map(years)
    top = off.sort_values("n", ascending=False).head(18).sort_values("year")
    mean_pa = p_med[15]
    ypos = np.arange(len(top))
    ax.errorbar(top["offset"] + mean_pa, ypos, xerr=top["err"], fmt="o", ms=4.5, color=BLUE, mec=SURFACE,
                ecolor=BLUE, elinewidth=0.9, capsize=0)
    ax.axvline(mean_pa, color=BLUE_300, lw=1.0, ls="-", zorder=0)
    ax.axvline(0, color=AXIS, lw=0.8, zorder=0)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{n} ({int(c)}, {y:.0f})" for n, c, y in zip(top["name"], top["n"], top["year"])],
                       fontsize=7)
    ax.set_xlabel("personal equation in position angle (deg)")
    ax.set_title("Visual observers' systematic errors")
    ax.text(mean_pa, len(top) - 0.3, f" population mean {mean_pa:+.2f}°", fontsize=7.5, color=BLUE_550,
            va="bottom")
    ax.set_ylim(-0.8, len(top) + 0.6)
    ax.grid(axis="y", visible=False)
    # (c) period from different data choices
    ax = axs[2]
    rows = [("Bond et al. 2017 (Table 4)", L.BOND17_ELEMENTS["P"][0] * 365.242198781 / 365.25,
             L.BOND17_ELEMENTS["P"][1], MUTED, "^"),
            ("their weights, our code", fit["replication"]["theta"][0], fit["replication"]["err"][0], INK_SECONDARY, "v"),
            ("all data, observers modelled (ours)", fit["hier"]["summary"]["P"]["median"],
             fit["hier"]["summary"]["P"]["std"], BLUE, "o")]
    for lab, v in fit["alternatives"].items():
        rows.append((lab + " only" if lab != "HST only" else "Hubble only", v["summary"]["P"]["median"],
                     v["summary"]["P"]["std"], ORANGE if "HST only" in lab else AQUA, "s" if "HST" in lab else "D"))
    for k, (lab, v, e, col, mk) in enumerate(rows):
        ax.errorbar(v, k, xerr=e, fmt=mk, ms=6, color=col, mec=SURFACE, ecolor=col, elinewidth=1.0, capsize=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.set_xlim(49.95, 50.30)
    ax.set_xlabel("orbital period (Julian years)")
    ax.set_title("How long is a Sirius year?")
    ax.grid(axis="y", visible=False)
    ax.text(0.98, 0.02, "Hubble only: 50.06 ± 0.60 yr\n(19 years of data)", transform=ax.transAxes, fontsize=7,
            color=MUTED, ha="right", va="bottom")
    caption(fig, "Left: our orbit against Bond et al.'s Hubble astrometry. Middle: posterior mean offset of each "
                 "visual observer with 20+ measures (name, number of measures, median year). Right: the period "
                 "depends on how the visual data are modelled; Bond et al.'s P is converted from Besselian to "
                 "Julian years.", y=-0.03)
    return save(fig, FIG / "sirius_systematics")


# ---------------------------------------------------------------------------
# 3. Sirius: masses
# ---------------------------------------------------------------------------


def fig_sirius_masses(res):
    ms = res["masses"]
    s = ms["samples"]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.5), gridspec_kw={"wspace": 0.34})
    ax = axs[0]
    ax.plot(s["MA"][::4], s["MB"][::4], ".", ms=1.5, color=BLUE_300, alpha=0.35, zorder=1, rasterized=True)
    ax.errorbar(np.median(s["MA"]), np.median(s["MB"]), xerr=np.std(s["MA"]), yerr=np.std(s["MB"]), fmt="o",
                ms=7, color=BLUE, mec=SURFACE, ecolor=BLUE, elinewidth=1.2, capsize=0, zorder=4, label="this project")
    for lab, mm, col, mk in (("Bond et al. 2017", L.BOND17_MASSES, ORANGE, "s"),
                             ("Gatewood & Gatewood 1978", L.GG78_MASSES, AQUA, "D"),
                             ("van den Bos 1960", L.VDB60_MASSES, MUTED, "^")):
        ax.errorbar(mm["MA"][0], mm["MB"][0], xerr=mm["MA"][1], yerr=mm["MB"][1], fmt=mk, ms=6.5, color=col,
                    mec=SURFACE, ecolor=col, elinewidth=1.0, capsize=0, label=lab, zorder=3)
    ax.set_xlabel("mass of Sirius A (M☉)")
    ax.set_ylabel("mass of Sirius B (M☉)")
    ax.set_title("Weighing Sirius A and B")
    ax.legend(loc="lower right", fontsize=7.5)
    # (b) error budget
    ax = axs[1]
    rows = ms["budget"]["rows"]
    labs = list(rows.keys())
    y = np.arange(len(labs))
    ax.barh(y + 0.18, [rows[k][0] for k in labs], height=0.34, color=BLUE, label="Sirius A")
    ax.barh(y - 0.18, [rows[k][1] for k in labs], height=0.34, color=ORANGE, label="Sirius B")
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("uncertainty it causes (M☉)")
    ax.set_title("The parallax is everything")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.grid(axis="y", visible=False)
    # (c) M_B for each parallax
    ax = axs[2]
    gaia = res["gaia"]
    P_, a_ = np.median(res["sirius"]["hier"]["elements"][:, 0]), np.median(res["sirius"]["hier"]["elements"][:, 3])
    fB = L.BOND17_AA["value"] / a_
    entries = [("ground-based plates (GG78)", L.BOND17_PLX["GG78"]["value"] * 1e3,
                L.BOND17_PLX["GG78"]["err"] * 1e3, AQUA, "D"),
               ("Hipparcos (van Leeuwen 2007)", L.BOND17_PLX["HIP2"]["value"] * 1e3, L.BOND17_PLX["HIP2"]["err"] * 1e3,
                ORANGE, "s"),
               ("weighted mean (adopted)", L.BOND17_PLX["adopted"]["value"] * 1e3,
                L.BOND17_PLX["adopted"]["err"] * 1e3, BLUE, "o"),
               ("Gaia DR3 (Sirius B), as published", gaia["obs"]["plx"], gaia["obs"]["e_plx"], MUTED, "v"),
               ("Gaia DR3, curvature-corrected (ours)", *gaia["plx_corrected"], INK, "*")]
    for k, (lab, p, e, col, mk) in enumerate(entries):
        M = orbits.total_mass(a_, p / 1e3, P_) * fB
        dM = M * 3 * e / p
        ax.errorbar(M, k, xerr=dM, fmt=mk, ms=13 if mk == "*" else 7, color=col, mec=SURFACE, ecolor=col,
                    elinewidth=1.2, capsize=0)
        ax.text(M + dM + 0.004, k, f"ϖ = {p:.1f} ± {e:.1f}", fontsize=7, color=INK_SECONDARY, va="center")
    ax.set_yticks(range(len(entries)))
    ax.set_yticklabels([e[0] for e in entries], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("mass of Sirius B (M☉)")
    ax.set_title("Which parallax?")
    ax.set_xlim(0.97, 1.11)
    ax.set_ylim(len(entries) - 0.4, -0.6)
    ax.grid(axis="y", visible=False)
    caption(fig, "Masses from our orbit, Sirius A's orbit size a_A = 2.4761\" ± 0.0045\" (Bond et al. 2017, from "
                 "Gatewood & Gatewood's plates) and a parallax. Mass scales as parallax⁻³: the 1.4 mas parallax error "
                 "is 1.1% in mass.", y=-0.03)
    return save(fig, FIG / "sirius_masses")


# ---------------------------------------------------------------------------
# 4. Sirius B in Gaia DR3
# ---------------------------------------------------------------------------


def fig_sirius_gaia(res):
    g = res["gaia"]
    el = np.median(res["sirius"]["hier"]["elements"], axis=0)
    fA = 1 - L.BOND17_AA["value"] / el[3]
    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.6), gridspec_kw={"wspace": 0.3})
    # (a) curvature of B's path over the DR3 window
    ax = axs[0]
    t = np.linspace(2014.57, 2017.39, 400)
    xn, ye = orbits.relative_xy(t, *el)
    e_, n_ = fA * ye * 1e3, fA * xn * 1e3
    for comp, col, lab, mk in ((e_, BLUE, "east", "o"), (n_, ORANGE, "north", "s")):
        p = np.polyfit(t - 2016, comp, 1)
        r = comp - np.polyval(p, t - 2016)
        ax.plot(t, r, color=col, lw=1.6, label=f"{lab}: departs by up to {np.abs(r).max():.0f} mas")
    tr = g["transits"]
    vp = __import__("gaia_sirius").visibility_periods(tr)
    ax.plot(tr["year"], np.full(len(tr), -21.0), "|", ms=9, color=INK_SECONDARY, mew=1.0)
    ax.text(2017.38, 20.5, f"ticks: {len(tr)} transits in {vp.max() + 1} visibility periods (scanning law);\n"
                           f"Gaia used {int(g['record']['matched_transits'])} transits in "
                           f"{int(g['record']['visibility_periods_used'])}",
            fontsize=7, color=MUTED, ha="right", va="top")
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_ylim(-23, 22)
    ax.set_xlabel("year")
    ax.set_ylabel("Sirius B minus uniform motion (mas)")
    ax.set_title("Sirius B curves during Gaia DR3")
    ax.legend(loc="upper left", fontsize=7.5, bbox_to_anchor=(0.0, 0.88))
    # (b) simulated parallax bias
    ax = axs[1]
    mc = g["mc"]
    bins = np.linspace(-15, 12, 55)
    ax.hist(mc["plx_bias"], bins=bins, color=BLUE_150, label="all 20,000 simulated DR3 solutions")
    ax.hist(mc["plx_bias"], bins=bins, weights=g["weights"] * len(mc) * 0.25, color=BLUE,
            label="those that also reproduce Gaia's\nproper motion (scaled ×0.25)")
    ax.axvline(g["d_obs"]["plx"], color=ORANGE, lw=2.0)
    ax.text(g["d_obs"]["plx"] - 0.3, ax.get_ylim()[1] * 0.92, f"observed\n{g['d_obs']['plx']:+.1f} mas", ha="right",
            fontsize=7.5, color=ORANGE, va="top")
    ax.set_xlabel("parallax bias of a 5-parameter fit (mas)")
    ax.set_ylabel("simulations")
    ax.set_title("A straight-line model shifts the parallax")
    ax.legend(loc="upper right", fontsize=7, bbox_to_anchor=(1.0, 0.86))
    b, sb = g["bias_given_pm"]
    ax.text(0.98, 0.55, f"given the proper motion:\nbias = {b:+.1f} ± {sb:.1f} mas", transform=ax.transAxes,
            ha="right", fontsize=7.5, color=BLUE_550)
    # (c) parallaxes compared
    ax = axs[2]
    entries = [("Gatewood & Gatewood 1978 (plates)", L.BOND17_PLX["GG78"]["value"] * 1e3, L.BOND17_PLX["GG78"]["err"] * 1e3,
                AQUA, "D"),
               ("Hipparcos (1989-93)", L.BOND17_PLX["HIP2"]["value"] * 1e3, L.BOND17_PLX["HIP2"]["err"] * 1e3, ORANGE, "s"),
               ("Gaia DR2 (Sirius B)", L.JOYCE18["plx_dr2"][0] * 1e3, L.JOYCE18["plx_dr2"][1] * 1e3, MUTED, "^"),
               ("Gaia DR3 (Sirius B), as published", g["obs"]["plx"], g["obs"]["e_plx"], MUTED, "v"),
               ("Gaia DR3 corrected (this project)", *g["plx_corrected"], INK, "*")]
    for k, (lab, p, e, col, mk) in enumerate(entries):
        ax.errorbar(p, k, xerr=e, fmt=mk, ms=13 if mk == "*" else 7, color=col, mec=SURFACE, ecolor=col,
                    elinewidth=1.2, capsize=0)
        ax.text(float(np.clip(p, 377.2, 379.0)), k - 0.28, f"{lab}: {p:.1f} ± {e:.1f}", ha="center", va="bottom",
                fontsize=7.5, color=INK_SECONDARY)
    ax.axvspan(L.BOND17_PLX["adopted"]["value"] * 1e3 - L.BOND17_PLX["adopted"]["err"] * 1e3,
               L.BOND17_PLX["adopted"]["value"] * 1e3 + L.BOND17_PLX["adopted"]["err"] * 1e3, color=WASH, zorder=0)
    ax.set_yticks([])
    ax.set_ylim(len(entries) - 0.5, -0.9)
    ax.set_xlim(372.5, 383.5)
    ax.set_xlabel("parallax of Sirius (mas)")
    ax.set_title("Corrected, Gaia agrees with Hipparcos")
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    ax.text(0.02, 0.02, "band: adopted 378.9 ± 1.4 mas", transform=ax.transAxes, fontsize=7, color=MUTED)
    caption(fig, "Sirius B's path from our orbit and Sirius A's orbit size; transits of its position from Gaia's DR3 "
                 "commanded scanning law minus DPAC's data gaps. Each simulation picks 10 of the 19 visibility "
                 "periods at random and fits Gaia's 5-parameter model to along-scan positions.", y=-0.03)
    return save(fig, FIG / "sirius_gaia")


# ---------------------------------------------------------------------------
# 5. White dwarf mass-radius
# ---------------------------------------------------------------------------


def fig_mass_radius(res):
    wd = res["wd"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.25})
    for k, ax in enumerate(axs):
        ch, fe = wd["chandra"], wd["chandra_fe"]
        ax.plot(ch["M"], ch["R"] * 100, color=INK, lw=1.3, ls="--", label="Chandrasekhar: cold, μₑ = 2 (C, O)")
        ax.plot(fe["M"], fe["R"] * 100, color=MUTED, lw=1.0, ls=":", label="cold iron (μₑ = 2.15)")
        mr = wd["montreal"]
        temps = sorted(mr)
        for T, col, lw in zip(temps, (BLUE_150, BLUE, BLUE_300), (1.0, 1.8, 1.0)):
            lab = f"C/O model at {T:,.0f} K" if abs(T - L.BOND17_WD["Teff"][0]) > 1 else \
                f"C/O model at Sirius B's {T:,.0f} K"
            ax.plot(mr[T]["M"], mr[T]["R"] * 100, color=col, lw=lw, label=lab)
        p = wd["parsons"]
        ax.errorbar(p["M"], p["R"] * 100, xerr=p["e_M"], yerr=p["e_R"] * 100, fmt="o", ms=3.5, color=CONTEXT,
                    ecolor=CONTEXT, elinewidth=0.7, capsize=0, label="eclipsing white dwarfs (Parsons et al. 2017)",
                    zorder=2)
        v = wd["visual_wds"]
        v = v[v["name"] != "Sirius B"]
        ax.errorbar(v["M"], v["R"] * 100, xerr=v["e_M"], yerr=v["e_R"] * 100, fmt="D", ms=5.5, color=AQUA,
                    mec=SURFACE, ecolor=AQUA, elinewidth=1.0, capsize=0, label="Procyon B, 40 Eri B, Stein 2051 B",
                    zorder=3)
        MB = wd["MB"]
        for lab, mk, col in (("V (HST/STIS, Barstow et al. 2005)", "o", ORANGE), ("Gaia G (DR3)", "s", YELLOW)):
            R = wd["radii"][lab]["R"]
            ax.errorbar(MB["median"], R["median"] * 100, xerr=MB["std"], yerr=R["std"] * 100, fmt=mk,
                        ms=7, color=col, mec=SURFACE, ecolor=col, elinewidth=1.2, capsize=0, zorder=5,
                        label=f"Sirius B: our mass, radius from {lab.split(' (')[0] if 'Gaia' in lab else 'HST V'}")
        bw = L.BOND17_WD
        ax.errorbar(L.BOND17_MASSES["MB"][0], bw["R"][0] * 100, xerr=L.BOND17_MASSES["MB"][1], yerr=bw["R"][1] * 100,
                    fmt="^", ms=6, color=MUTED, mec=SURFACE, ecolor=MUTED, elinewidth=1.0, capsize=0, zorder=4,
                    label="Sirius B, Bond et al. 2017")
        if k == 0:
            ax.set_xlim(0.0, 1.5)
            ax.set_ylim(0.0, 3.2)
            ax.set_title("White dwarfs: heavier is smaller")
            ax.legend(loc="upper right", fontsize=7)
            for name, x, y in (("Sirius B", 1.02, 0.81),):
                annotate(ax, name, (x, y), (x + 0.12, y + 0.6))
            ax.text(1.47, 0.62, "Chandrasekhar\nlimit 1.45 M☉", fontsize=7.5, color=INK_SECONDARY, ha="right")
        else:
            ax.set_xlim(0.93, 1.10)
            ax.set_ylim(0.70, 0.92)
            ax.set_title("Sirius B up close")
        ax.set_xlabel("mass (M☉)")
        ax.set_ylabel("radius (R☉ / 100;  1 ≈ 7,000 km)")
    vgr = wd["vgr"]
    axs[1].text(0.03, 0.04, f"gravitational redshift GM/Rc from our M and R:\n{vgr['median']:.1f} ± {vgr['std']:.1f} km/s "
                            f"(measured: {L.JOYCE18['vgr'][0]} ± {L.JOYCE18['vgr'][1]}, Joyce et al. 2018)",
                transform=axs[1].transAxes, fontsize=7.5, color=INK_SECONDARY)
    caption(fig, "Chandrasekhar curves integrated here (exact degenerate-electron equation of state). C/O models: "
                 "Montreal thick-H sequences (Bédard et al. 2020). Sirius B's radius from its magnitude, our parallax "
                 "and the Montreal surface brightness at 25,369 K (Bond et al. 2017's temperature).", y=-0.02)
    return save(fig, FIG / "whitedwarf_mass_radius")


# ---------------------------------------------------------------------------
# 6-8. The mass-luminosity relation
# ---------------------------------------------------------------------------


def fig_mass_luminosity(res):
    ml = res["mlr"]
    s = ml["sample"]
    like, best = ml["like"], ml["best"]
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.3), gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.25})
    ax = axs[0]
    ek = ml["eker_stars"]
    sel = ek["selected"]
    ax.plot(ek.loc[sel, "MV"], ek.loc[sel, "M"], ".", ms=3.5, color=CONTEXT, zorder=1,
            label=f"eclipsing-binary stars ({int(sel.sum())}, Eker et al. 2018)")
    MV, band = ml["MV"], ml["band"]
    ax.fill_between(MV, 10 ** band[0], 10 ** band[4], color=BLUE_150, lw=0, zorder=2, label="95% band")
    ax.fill_between(MV, 10 ** band[1], 10 ** band[3], color=BLUE_300, lw=0, zorder=3, label="68% band")
    ax.plot(MV, 10 ** band[2], color=BLUE_550, lw=2.0, zorder=4,
            label=f"from {ml['n_fit']} visual-binary orbits (this project)")
    ec = ml["eker_curve"]
    ok = np.isfinite(ec["MV"])
    ax.plot(ec["MV"][ok], ec["M"][ok], color=ORANGE, lw=1.6, zorder=5, label="Eker et al. 2018 relation (→ M_V)")
    mm = ml["mamajek"].dropna(subset=["Mv", "Msun"])
    mm = mm[mm["Mv"].between(-2, 16)]
    ax.plot(mm["Mv"], mm["Msun"], color=AQUA, lw=1.4, ls="--", zorder=5, label="Pecaut & Mamajek dwarf sequence")
    ax.set_yscale("log")
    ax.set_xlim(-1.5, 15.5)
    ax.set_ylim(0.09, 8)
    ax.invert_xaxis()
    ax.set_xlabel("absolute visual magnitude M_V  (brighter →)")
    ax.set_ylabel("mass of one star (M☉)")
    ax.set_title("The mass-luminosity relation, from orbits alone")
    ax.legend(loc="lower right", fontsize=7.5)
    for mv, lab in ((4.81, "the Sun"),):
        m_sun = 10 ** np.interp(mv, MV, band[2])
        annotate(ax, f"a star as bright as {lab}\n(M_V = {mv}): {m_sun:.2f} M☉", (mv, m_sun), (10.5, 2.6))
    ax.set_yticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_yticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    # residuals of the totals
    ax = axs[1]
    d = s[s["in_fit"]]
    gaia = d["plx_source"].str.startswith("Gaia")
    out = d["p_outlier"] > 0.5
    for m, col, mk, lab in ((~gaia & ~out, ORANGE, "s", "Hipparcos parallax"), (gaia & ~out, BLUE, "o", "Gaia DR3 parallax")):
        ax.errorbar(d.loc[m, "MV1"], d.loc[m, "resid"], yerr=d.loc[m, "resid_err"], fmt=mk, ms=4, color=col,
                    mec=SURFACE, ecolor=col, elinewidth=0.6, capsize=0, alpha=0.9, label=f"{lab} ({int(m.sum())})")
    ax.plot(d.loc[out, "MV1"], d.loc[out, "resid"], "x", ms=6, color=MUTED, mew=1.2,
            label=f"outliers ({int(out.sum())})")
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_ylim(-0.6, 0.6)
    ax.invert_xaxis()
    ax.set_xlabel("M_V of the brighter star")
    ax.set_ylabel("log(dynamical mass / relation's prediction)")
    ax.set_title(f"Totals scatter by {ml['sigma_int'] * 100 * np.log(10):.0f}% beyond their errors")
    ax.legend(loc="upper right", fontsize=7.5)
    caption(fig, "Each orbit gives M1 + M2 = (a/ϖ)³/P². The relation log M = f(M_V) (piecewise linear, monotonic) is "
                 "fitted so that 10^f(M_V1) + 10^f(M_V2) matches every total, with parallax errors moving masses and "
                 "magnitudes together, an intrinsic scatter and an outlier mixture. Eclipsing stars: L from R and Teff, "
                 "BC_V from Mamajek's table.", y=-0.03)
    return save(fig, FIG / "mass_luminosity")


def fig_total_mass(res):
    ml = res["mlr"]
    s = ml["sample"]
    like, best = ml["like"], ml["best"]
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    ctx = s[(s["plx_source"] != "") & s["v_ok"] & ~s["in_fit"]]
    ax.plot(ctx["Mtot"], ctx["MV_comb"], ".", ms=4, color=CONTEXT, label=f"not in the fit ({len(ctx)}: evolved, "
                                                                         f"triples, >100 pc, poor mass)")
    d = s[s["in_fit"]]
    sc = ax.scatter(d["Mtot"], d["MV_comb"], c=np.clip(d["dm"], 0, 4), cmap=matplotlib.colors.LinearSegmentedColormap
                    .from_list("b", [BLUE_150, BLUE_550]), s=22, edgecolor=SURFACE, linewidth=0.4, zorder=3,
                    label=f"in the fit ({len(d)})")
    cb = fig.colorbar(sc, ax=ax, pad=0.01)
    cb.set_label("magnitude difference V₂ - V₁")
    mv = np.linspace(-1.0, 15.5, 300)
    m1 = 10 ** like.f(best, mv)
    ax.plot(2 * m1, mv - 2.5 * np.log10(2), color=INK, lw=1.4, label="identical twins (relation)")
    ax.plot(m1 + 10 ** like.f(best, mv + 5), -2.5 * np.log10(10 ** (-0.4 * mv) + 10 ** (-0.4 * (mv + 5))),
            color=INK, lw=1.2, ls="--", label="one star dominates (Δm = 5)")
    ax.set_xscale("log")
    ax.set_xlim(0.1, 12)
    ax.set_ylim(15.5, -2)
    ax.set_xlabel("total dynamical mass M₁ + M₂ (M☉)")
    ax.set_ylabel("total absolute magnitude M_V of the pair")
    ax.set_title("Heavier pairs shine brighter: every grade 1-2 orbit with a parallax")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10"])
    ax.legend(loc="upper left", fontsize=7.5)
    caption(fig, "ORB6 (October 2026) grade 1-2 orbits; parallaxes from Gaia DR3 where clean (RUWE < 1.4), else "
                 "Hipparcos. Twins sit on the solid line, pairs with a faint companion on the dashed one.", y=-0.01)
    return save(fig, FIG / "total_mass_luminosity")


def fig_gaia_vs_hip(res):
    ml = res["mlr"]
    b = ml["gaia_vs_hip"]
    s = ml["sample"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"wspace": 0.25})
    ax = axs[0]
    ax.scatter(b["ruwe_p"], b["dplx_sig"], s=14, color=BLUE, edgecolor=SURFACE, linewidth=0.3)
    ax.axvline(C.GAIA_RUWE_MAX, color=ORANGE, lw=1.2)
    ax.axhspan(-3, 3, color=WASH, zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=3)
    ax.set_yticks([-30, -10, -3, 0, 3, 10, 30])
    ax.set_yticklabels(["-30", "-10", "-3", "0", "3", "10", "30"])
    ax.set_xlabel("Gaia DR3 RUWE of the primary (or of the blended pair)")
    ax.set_ylabel("(ϖ Gaia - ϖ Hipparcos) / combined error")
    ax.set_title("Wobbling binaries: Gaia parallaxes fail more often")
    good = b["ruwe_p"] < C.GAIA_RUWE_MAX
    for m, x in ((good, 0.3), (~good, 0.98)):
        frac = np.mean(np.abs(b.loc[m, "dplx_sig"]) < 3)
        ax.text(x, 0.04, f"{frac:.0%} agree within 3σ\n(n = {int(m.sum())})", transform=ax.transAxes,
                ha="right", fontsize=7.5, color=INK_SECONDARY)
    ax.text(C.GAIA_RUWE_MAX * 1.05, ax.get_ylim()[1] * 0.6, "RUWE = 1.4", color=ORANGE, fontsize=7.5)
    ax = axs[1]
    g = s.copy()
    g["sep_bin"] = pd.cut(g["a_arcsec"], [0, 0.05, 0.1, 0.2, 0.5, 1, 30])
    src = g.groupby("sep_bin", observed=True)["plx_source"].agg(
        gaia=lambda x: np.mean(x.str.startswith("Gaia")), hip=lambda x: np.mean(x == "Hipparcos"),
        none=lambda x: np.mean(x == ""), n="size")
    x = np.arange(len(src))
    ax.bar(x, src["gaia"], color=BLUE, width=0.7, label="clean Gaia DR3 parallax")
    ax.bar(x, src["hip"], bottom=src["gaia"], color=ORANGE, width=0.7, label="Hipparcos instead")
    ax.bar(x, src["none"], bottom=src["gaia"] + src["hip"], color=CONTEXT, width=0.7, label="no usable parallax")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{iv.left:g}-{iv.right:g}\"\n(n={n})" for iv, n in zip(src.index, src["n"])], fontsize=7.5)
    ax.set_xlabel("angular semi-major axis of the orbit")
    ax.set_ylabel("fraction of grade 1-2 orbits")
    ax.set_title("Where each parallax comes from")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.grid(axis="x", visible=False)
    caption(fig, "Left: ORB6 grade 1-2 systems with both a Gaia DR3 (5/6-parameter) and a Hipparcos parallax. "
                 "Right: parallax source chosen per orbit (Gaia needs RUWE < 1.4 and 5% precision).", y=-0.02)
    return save(fig, FIG / "gaia_vs_hipparcos")


# ---------------------------------------------------------------------------
# 9. The gallery of famous orbits
# ---------------------------------------------------------------------------

PLANETS = (("Jupiter", 5.20), ("Saturn", 9.54), ("Uranus", 19.19), ("Neptune", 30.07))


def gallery_xy(entry, t):
    """Barycentric positions (au, face-on orbital plane) of both stars at times t."""
    P, T, e, a, i, O, w = entry["elements"]
    a_au = entry["a_au"]
    X, Y = orbits.unit_orbit(t, P, T, e)
    xr, yr = a_au * X, a_au * Y          # B relative to A, periastron along +x
    fB = entry["fB"]
    return (-fB * xr, -fB * yr), ((1 - fB) * xr, (1 - fB) * yr)


def gallery_extent(gal):
    """Largest distance of either star from its centre of mass, over all six orbits (au)."""
    m = 0.0
    for entry in gal:
        P = entry["elements"][0]
        t = entry["elements"][1] + P * np.linspace(0, 1, 721)
        (xa, ya), (xb, yb) = gallery_xy(entry, t)
        m = max(m, np.hypot(xa, ya).max(), np.hypot(xb, yb).max())
    return m


def fig_gallery(res):
    gal = res["gallery"]
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.6), gridspec_kw={"hspace": 0.22, "wspace": 0.08})
    lim = 1.05 * gallery_extent(gal)
    for ax, entry in zip(axs.flat, gal):
        P = entry["elements"][0]
        t = entry["elements"][1] + P * np.linspace(0, 1, 721)
        (xa, ya), (xb, yb) = gallery_xy(entry, t)
        for name, r in PLANETS:
            ax.add_patch(plt.Circle((0, 0), r, fill=False, color=GRID, lw=0.8, zorder=0))
        ax.plot(xa, ya, color=BLUE, lw=1.6, label="A")
        ax.plot(xb, yb, color=ORANGE, lw=1.6, label="B")
        ax.plot([0], [0], "+", color=MUTED, ms=7)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim * 0.8, lim * 0.8)
        ax.set_aspect("equal")
        ax.set_title(entry["name"])
        ax.text(0.03, 0.04, f"P = {P:.1f} yr   a = {entry['a_au']:.1f} au\n"
                            f"M = {entry['MA']:.2f} + {entry['MB']:.2f} M☉",
                transform=ax.transAxes, fontsize=8, color=INK_SECONDARY)
        ax.tick_params(labelsize=7)
        ax.grid(False)
    axs[0, 0].text(30.07 * 0.71, 30.07 * 0.71, "Neptune", fontsize=7, color=MUTED)
    axs[0, 0].text(5.2, 5.6, "Jupiter", fontsize=6.5, color=MUTED)
    for ax in axs[1]:
        ax.set_xlabel("au")
    fig.suptitle("Six famous pairs to one scale: primary (blue) and companion (orange) about their centre of mass",
                 fontsize=11, fontweight="bold", x=0.5, y=0.96)
    caption(fig, "Elements from ORB6 (October 2026; Sirius: this project). Distances: Sirius and Procyon from Bond et "
                 "al.; α Cen from Akeson et al. 2021; others Gaia DR3 or Hipparcos. Mass ratios: measured for Sirius, "
                 "Procyon, α Cen; from our mass-luminosity relation for ξ UMa, Castor (Castor A and B are each "
                 "themselves binaries) and 61 Cyg.", y=0.04)
    return save(fig, FIG / "orbit_gallery")


def make_movie(res, path=FIG / "orbit_gallery.mp4", years_per_second=12.0, fps=30, duration=40.0):
    """Animation: all six pairs moving together, to one scale, with a shared clock."""
    from matplotlib.animation import FFMpegWriter
    gal = res["gallery"]
    fig, axs = plt.subplots(2, 3, figsize=(12.8, 7.2), dpi=100)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.04, wspace=0.08, hspace=0.18)
    lim = 1.03 * gallery_extent(gal)
    arts = []
    t0 = 2026.0
    for ax, entry in zip(axs.flat, gal):
        P = entry["elements"][0]
        t = entry["elements"][1] + P * np.linspace(0, 1, 721)
        (xa, ya), (xb, yb) = gallery_xy(entry, t)
        for name, r in PLANETS:
            ax.add_patch(plt.Circle((0, 0), r, fill=False, color=GRID, lw=0.7))
        ax.plot(xa, ya, color=BLUE_150, lw=1.0)
        ax.plot(xb, yb, color="#f6c3ad", lw=1.0)
        ax.plot([0], [0], "+", color=MUTED, ms=6)
        sa = 7 + 3 * np.log10(max(entry["MA"], 0.1) * 10)
        sb = 7 + 3 * np.log10(max(entry["MB"], 0.1) * 10)
        pa, = ax.plot([], [], "o", color=BLUE, ms=sa, mec=SURFACE)
        pb, = ax.plot([], [], "o", color=ORANGE, ms=sb, mec=SURFACE)
        ta, = ax.plot([], [], color=BLUE, lw=1.8, alpha=0.7)
        tb, = ax.plot([], [], color=ORANGE, lw=1.8, alpha=0.7)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim * 0.7, lim * 0.7)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"{entry['name']}   P = {P:.0f} yr", fontsize=10)
        arts.append((entry, pa, pb, ta, tb))
    clock = fig.text(0.5, 0.93, "", ha="center", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.5, 0.965, "Six famous binaries to one scale (au); grey rings: orbits of Jupiter, Saturn, Uranus, Neptune",
             ha="center", fontsize=9, color=INK_SECONDARY)
    nframes = int(duration * fps)
    writer = FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "28"])
    with writer.saving(fig, str(path), dpi=100):
        for k in range(nframes):
            t = t0 + k * years_per_second / fps
            clock.set_text(f"year {t:.0f}")
            for entry, pa, pb, ta, tb in arts:
                (xa, ya), (xb, yb) = gallery_xy(entry, np.array([t]))
                pa.set_data(xa, ya)
                pb.set_data(xb, yb)
                tt = np.linspace(t - min(8.0, entry["elements"][0] * 0.15), t, 40)
                (xa2, ya2), (xb2, yb2) = gallery_xy(entry, tt)
                ta.set_data(xa2, ya2)
                tb.set_data(xb2, yb2)
            writer.grab_frame()
    plt.close(fig)
    return path


def make_all(res, movie=True):
    paths = []
    for fn in (fig_sirius_orbit, fig_sirius_systematics, fig_sirius_masses, fig_sirius_gaia, fig_mass_radius,
               fig_mass_luminosity, fig_total_mass, fig_gaia_vs_hip, fig_gallery):
        paths += fn(res)
    import interactive
    paths.append(interactive.write_page(res))
    if movie:
        paths.append(make_movie(res))
    return paths
