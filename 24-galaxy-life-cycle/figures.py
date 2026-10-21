"""Every figure, drawn from what analysis.analyse() returns.

Nothing here measures anything; it only draws. Each function writes one
figure (PNG + SVG) into figures/ and returns the paths.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm

import galaxies as gx
import literature as lit
import plotstyle as ps
import sfh
from config import FIGURE_DIR, IMF_TO_SALPETER_SFR

ps.use_style()

DENSITY_CMAP = LinearSegmentedColormap.from_list("blue_ramp", ps.BLUE_RAMP)


def _out(name: str):
    return FIGURE_DIR / name


def _time_axis(ax, cosmo=sfh.COSMO, ages=(0.6, 1, 2, 3, 5, 8)):
    """Top axis: age of the universe at the redshifts on the bottom log(1+z) axis."""
    top = ax.twiny()  # same y-axis: a second labelling of x, not a second scale
    top.set_xlim(ax.get_xlim())
    ticks = [np.log10(1 + sfh.find_redshift_at_time(a, cosmo)) for a in ages] + [0.0]
    t0 = float(sfh.cosmic_time_gyr(0.0, cosmo))
    top.set_xticks(ticks)
    top.set_xticklabels([f"{a:g}" for a in ages] + [f"{t0:.1f}\n(today)"])
    top.set_xlabel("age of the universe (billion years)")
    top.grid(False)
    for s in ("top",):
        top.spines[s].set_visible(True)
        top.spines[s].set_color(ps.AXIS)
    return top


def _z_ticks(ax, zs=(0, 0.5, 1, 2, 3, 4, 6, 8)):
    ax.set_xticks([np.log10(1 + z) for z in zs])
    ax.set_xticklabels([f"{z:g}" for z in zs])


# ---------------------------------------------------------------------------
# Cosmic star-formation history
# ---------------------------------------------------------------------------


def fig_cosmic_sfh(res: dict) -> list:
    """The keeper: when the universe made its stars."""
    h = res["sfh"]
    pts = h["points"]
    pub = h["published"]
    ours = h["ours"]
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    x = np.log10(1 + pts["z"])
    xlo = np.log10(1 + pts["z_lo"])
    xhi = np.log10(1 + pts["z_hi"])
    for band, color, label, marker in (("UV", ps.BLUE, "ultraviolet surveys (dust-corrected)", "o"),
                                       ("IR", ps.RED, "infrared surveys", "s")):
        m = pts["band"] == band
        ax.errorbar(x[m], pts["logpsi"][m], xerr=[x[m] - xlo[m], xhi[m] - x[m]],
                    yerr=[pts["err_minus"][m], pts["err_plus"][m]], fmt=marker, ms=4.5,
                    color=color, ecolor=color, elinewidth=0.8, capsize=0, alpha=0.85,
                    label=label, zorder=3, mec=ps.SURFACE, mew=0.6)
    zz = np.linspace(0, 9, 600)
    xz = np.log10(1 + zz)
    # bootstrap band of our fit
    bp = ours["boot_params"]
    curves = np.array([np.log10(sfh.psi_md(zz, *b)) for b in bp])
    lo, hi = np.percentile(curves, [16, 84], axis=0)
    ax.fill_between(xz, lo, hi, color=ps.ORANGE, alpha=0.18, lw=0, zorder=1)
    ax.plot(xz, np.log10(sfh.psi_md(zz, ours["A"], ours["B"], ours["C"], ours["D"])),
            color=ps.ORANGE, lw=2, zorder=2, label="our fit to their Table 1 (band: 68%)")
    ax.plot(xz, np.log10(sfh.psi_md(zz, pub["A"], pub["B"], pub["C"], pub["D"])),
            color=ps.VIOLET, lw=1.6, ls="--", zorder=2,
            label="Madau & Dickinson's published curve")
    # our SDSS anchor: likelihood mean (star); the bar spans the three estimates
    cm, cmed, cg = res["census_mean"], res["census"], res["census_gswlc"]
    vals = [cmed["rho_sfr_salpeter"], cg["rho_sfr_salpeter"], cm["rho_sfr_salpeter"]]
    x01 = np.log10(1.1)
    ax.plot([x01, x01], [np.log10(min(vals)), np.log10(max(vals))], color=ps.INK, lw=1.2,
            zorder=4)
    ax.plot([x01], [np.log10(cg["rho_sfr_salpeter"])], "D", ms=5, color=ps.INK,
            mec=ps.SURFACE, mew=0.6, zorder=5, label="SDSS with GSWLC-2 SFRs")
    ax.plot([x01], [np.log10(cmed["rho_sfr_salpeter"])], "_", ms=10, mew=1.5, color=ps.INK,
            zorder=5)
    ax.plot([x01], [np.log10(cm["rho_sfr_salpeter"])], "*", ms=13, color=ps.INK,
            mec=ps.SURFACE, mew=0.8, zorder=6, label="this project: SDSS, z = 0.1 (MPA-JHU)")
    zp = ours["summary"]["z_peak"]
    ps.annotate(ax, f"peak: z = {zp:.2f}\n{ours['summary']['t_peak_gyr']:.1f} billion years "
                "after the Big Bang", (np.log10(1 + zp), np.log10(ours["summary"]["psi_peak"])),
                (np.log10(1 + zp) + 0.12, -0.55), fontsize=8)
    ax.set_xlim(0, np.log10(10))
    ax.set_ylim(-2.45, -0.35)
    _z_ticks(ax)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("log star-formation rate density  (M$_\\odot$ yr$^{-1}$ Mpc$^{-3}$, Salpeter IMF)")
    _time_axis(ax)
    ax.set_title("When the universe made its stars", pad=30)
    ax.legend(loc="lower center", ncol=2, fontsize=7.5)
    ps.caption(fig, "Each point is one survey's count of newborn stars per year in a typical cubic "
               "megaparsec (Madau & Dickinson 2014, Table 1, parsed from the paper's LaTeX "
               "source; 28 ultraviolet and 23 infrared measurements). Horizontal bars span "
               "each survey's redshift bin. The universe formed stars ~9x faster at the peak "
               "than today. The SDSS star is our own 1/Vmax sum over "
               f"{cm['n']:,} galaxies (means of the MPA-JHU SFR likelihoods, Kroupa IMF "
               "converted to Salpeter by /0.67); the bar runs down to the sum of their medians, "
               "the diamond uses GSWLC-2 SFRs. Statistical errors (jackknife) are 4-6%; the "
               "choice of SFR estimate matters more.", y=-0.03)
    return ps.save(fig, _out("06_cosmic_star_formation_history"))


def fig_mass_buildup(res: dict) -> list:
    """Integrating the history: how today's stars built up, vs direct counts."""
    h = res["sfh"]
    hp_, ho = h["history_pub"], h["history_ours"]
    smd = h["smd"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.6), gridspec_kw={"wspace": 0.28})
    # (a) stellar mass density vs redshift
    zs = h["smd_z"]
    x = np.log10(1 + zs)
    yerr = np.vstack([np.nan_to_num(smd["err_minus"].to_numpy(), nan=0.0),
                      np.nan_to_num(smd["err_plus"].to_numpy(), nan=0.0)])
    a1.errorbar(x, smd["log_value"], yerr=yerr, fmt="o", ms=3.5, color=ps.INK_SECONDARY,
                ecolor=ps.MUTED, elinewidth=0.7, capsize=0, label="direct counts of stellar mass "
                "(Table 2, 71 points)", zorder=3)
    for hist, color, ls, label in ((ho, ps.ORANGE, "-", "integral of our fit"),
                                   (hp_, ps.VIOLET, "--", "integral of their curve")):
        sel = hist["z"] <= 9
        a1.plot(np.log10(1 + hist["z"][sel]), np.log10(hist["rho_star"][sel]), color=color,
                ls=ls, lw=1.8, label=label)
    c = res["census"]
    a1.plot([np.log10(1 + c["z_eff_mass"])], [np.log10(c["rho_star_salpeter"])], "*", ms=13,
            color=ps.INK, mec=ps.SURFACE, mew=0.8, zorder=5, label="this project: SDSS")
    a1.set_xlim(0, np.log10(10))
    a1.set_ylim(5.5, 9.0)
    _z_ticks(a1)
    a1.set_xlabel("redshift z")
    a1.set_ylabel("log stellar mass density  (M$_\\odot$ Mpc$^{-3}$, Salpeter)")
    a1.set_title("Stars ever formed, minus what they return", pad=8)
    a1.legend(loc="lower left", fontsize=7.5)

    # (b) cumulative fraction formed vs cosmic time
    t = hp_["t_gyr"]
    f = hp_["formed_before"] / hp_["total_formed"]
    a2.plot(t, f, color=ps.VIOLET, ls="--", lw=1.6, label="their curve")
    to = ho["t_gyr"]
    fo = ho["formed_before"] / ho["total_formed"]
    a2.plot(to, fo, color=ps.ORANGE, lw=2, label="our fit")
    s = h["published"]["summary"]
    so = h["ours"]["summary"]
    marks = [(sfh.cosmic_time_gyr(6.0), so["frac_before_reionization_end"],
              f"end of reionization (z = 6): {100 * so['frac_before_reionization_end']:.1f}%"),
             (so["t_peak_gyr"], so["frac_before_peak"],
              f"peak of star formation: {100 * so['frac_before_peak']:.0f}% formed"),
             (sfh.cosmic_time_gyr(so["z_half_mass"]), 0.5,
              f"half of today's stars by z = {so['z_half_mass']:.2f}"),
             (sfh.cosmic_time_gyr(so["z_last_quarter"]), 0.75,
              f"last quarter after z = {so['z_last_quarter']:.2f}")]
    offsets = [(2.2, 0.16), (1.5, -0.02), (1.6, -0.08), (-6.5, 0.12)]
    for (tm, fm, text), (dx, dy) in zip(marks, offsets):
        a2.plot([tm], [fm], "o", ms=5, color=ps.INK, zorder=4)
        ps.annotate(a2, text, (tm, fm), (tm + dx, fm + dy), fontsize=7.5)
    a2.set_xlim(0, 14)
    a2.set_ylim(0, 1.02)
    a2.set_xlabel("age of the universe (billion years)")
    a2.set_ylabel("fraction of today's stars already formed")
    a2.set_title("The build-up, in cosmic time", pad=8)
    a2.legend(loc="lower right", fontsize=7.5)
    ps.caption(fig, f"Left: integrating the star-formation history with a return fraction R = 0.27 "
               f"(stars give 27% of their mass back as winds and supernovae) predicts today's "
               f"stellar mass density: {s['rho_star_z0'] / 1e8:.2f} x 10^8 M_sun/Mpc^3 from "
               f"their curve (MD14 text: 5.8), {so['rho_star_z0'] / 1e8:.2f} from our fit; the SDSS "
               f"census gives {c['rho_star_salpeter'] / 1e8:.2f}. Right: the same integral as a "
               "fraction, with the milestones Madau & Dickinson quote.", y=-0.03)
    return ps.save(fig, _out("07_stellar_mass_buildup"))


# ---------------------------------------------------------------------------
# The SFR-mass plane
# ---------------------------------------------------------------------------


def _plane_hist(logm, logsfr, w, xe, ye):
    h, _, _ = np.histogram2d(logm, logsfr, bins=[xe, ye], weights=w)
    return h


def fig_sfr_mass_plane(res: dict) -> list:
    """The keeper: the main sequence, the green valley and the quenched cloud."""
    g, w, ms = res["_g"], res["_w_local"], res["ms"]
    ok = ms["ok"]
    fig = plt.figure(figsize=(11.5, 6.4))
    gs_ = fig.add_gridspec(3, 2, width_ratios=[2.15, 1], wspace=0.22, hspace=0.12)
    ax = fig.add_subplot(gs_[:, 0])
    xe = np.arange(8.5, 12.01, 0.05)
    ye = np.arange(-3.6, 1.81, 0.05)
    h = _plane_hist(g["logm"][ok], g["logsfr"][ok], w[ok], xe, ye)
    n_cell, _, _ = np.histogram2d(g["logm"][ok], g["logsfr"][ok], bins=[xe, ye])
    # cells holding a handful of galaxies are dominated by one large 1/Vmax
    # weight (a nearby dwarf); show only cells with at least 5 galaxies
    h = np.where(n_cell >= 5, h, np.nan)
    h = h / np.nanmax(h)
    im = ax.pcolormesh(xe, ye, h.T, cmap=DENSITY_CMAP, norm=LogNorm(1e-3, 1), rasterized=True)
    ridge = ms["ridge"]
    good = np.isfinite(ridge["peak"])
    xx = np.linspace(8.6, 11.9, 50)
    xb = np.linspace(9.2, 11.9, 50)  # the band only where there are galaxies to divide
    v = ms["valley"]["position"]
    ax.fill_between(xb, gx.ms_sfr(xb, ridge) + v - 0.25, gx.ms_sfr(xb, ridge) + v + 0.25,
                    color=ps.C_GREEN_VALLEY, alpha=0.22, lw=0, zorder=2)
    ax.plot(xx, gx.ms_sfr(xx, ridge), color=ps.C_BLUE_GALAXY, lw=2.0, zorder=3)
    ax.plot(ridge["logm"][good], ridge["peak"][good], "o", ms=5, color=ps.C_BLUE_GALAXY,
            mec=ps.SURFACE, mew=0.8, zorder=4)
    rp = lit.RENZINI_PENG_2015
    ax.plot(xx, rp["slope"][0] * xx + rp["intercept"][0], color=ps.VIOLET, lw=1.3, ls="--",
            zorder=3)
    ax.text(11.85, gx.ms_sfr(11.85, ridge) + 0.12, "star-forming\nmain sequence",
            color=ps.C_BLUE_GALAXY, ha="right", va="bottom", fontsize=9, weight="bold")
    ax.text(11.85, gx.ms_sfr(11.85, ridge) + v + 0.30, "green valley", color="#8a5d00",
            ha="right", va="bottom", fontsize=9, weight="bold")
    ax.text(11.85, -2.9, "quenched", color=ps.C_RED_GALAXY, ha="right", fontsize=9,
            weight="bold")
    ps.annotate(ax, f"our ridge line: slope {ridge['slope']:.2f} "
                f"({ms['ridge_lowmass']['slope']:.2f} below 10$^{{10.5}}$)",
                (9.2, gx.ms_sfr(9.2, ridge)),
                (8.65, 1.05), color=ps.C_BLUE_GALAXY)
    ps.annotate(ax, f"Renzini & Peng 2015: slope {rp['slope'][0]:.2f}",
                (11.2, rp["slope"][0] * 11.2 + rp["intercept"][0]), (8.65, 1.5),
                color=ps.VIOLET)
    ax.set_xlim(8.5, 12.0)
    ax.set_ylim(-3.6, 1.8)
    ax.set_xlabel("log stellar mass  (M$_\\odot$)")
    ax.set_ylabel("log star-formation rate  (M$_\\odot$ per year)")
    ax.set_title("Two kinds of galaxy: forming stars, or done", pad=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.015)
    cb.set_label("relative number of galaxies per Mpc$^3$ (1/V$_{max}$)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    # right column: Delta_MS distributions in three mass slices
    dms = ms["dms"]
    for k, (a, b) in enumerate(((9.5, 10.0), (10.25, 10.75), (11.0, 11.5))):
        sub = fig.add_subplot(gs_[k, 1])
        m = ok & (g["logm"] >= a) & (g["logm"] < b)
        edges = np.arange(-3.5, 1.51, 0.05)
        hh, _ = np.histogram(dms[m], bins=edges, weights=w[m])
        hh = hh / hh.sum() / 0.05
        c = 0.5 * (edges[1:] + edges[:-1])
        for lo_, hi_, col in ((-9, v - 0.25, ps.C_RED_GALAXY),
                              (v - 0.25, v + 0.25, ps.C_GREEN_VALLEY),
                              (v + 0.25, 9, ps.C_BLUE_GALAXY)):
            zone = (c >= lo_ - 0.026) & (c <= hi_ + 0.026)
            sub.fill_between(c[zone], hh[zone], step="mid", color=col, lw=0, alpha=0.85)
        sub.step(c, hh, where="mid", color=ps.INK_SECONDARY, lw=0.6)
        fq = np.sum(w[m] * ms["quenched"][m]) / np.sum(w[m])
        sub.text(0.98, 0.92, f"log M = {a:g}–{b:g}\n{100 * fq:.0f}% quenched",
                 transform=sub.transAxes, ha="right", va="top", fontsize=7.5,
                 color=ps.INK_SECONDARY)
        sub.set_xlim(-3.5, 1.5)
        sub.set_yticks([])
        if k < 2:
            sub.set_xticklabels([])
        else:
            sub.set_xlabel("distance from the main sequence  (dex in SFR)")
        if k == 0:
            sub.set_title("Below the sequence, by mass", pad=6)
    ps.caption(fig, f"{ms['n']:,} SDSS galaxies at 0.02 < z < 0.085 with MPA-JHU stellar masses "
               "and star-formation rates, each weighted by 1/Vmax (the volume it could have been "
               "seen in) and for fibre collisions, so the map shows galaxies per unit volume, "
               "not per unit of survey. Dots: the mode of the star-forming peak in each 0.2-dex "
               "mass column (AGN excluded), as defined by Renzini & Peng (2015); line: our fit. "
               f"The amber band (0.5 dex wide) is centred on the valley floor, {-v:.2f} dex below "
               "the sequence.", y=-0.02)
    return ps.save(fig, _out("01_sfr_mass_plane"))


def fig_3d_surface(res: dict) -> list:
    """Renzini & Peng's 3D SFR-mass-number plot, and number x SFR."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    g, w, ms = res["_g"], res["_w_local"], res["ms"]
    ok = ms["ok"]
    xe = np.arange(9.0, 11.81, 0.2)
    ye = np.arange(-3.2, 1.61, 0.2)
    h = _plane_hist(g["logm"][ok], g["logsfr"][ok], w[ok] * 10 ** g["logm"][ok], xe, ye)
    hs = _plane_hist(g["logm"][ok], g["logsfr"][ok],
                     w[ok] * 10 ** np.clip(g["logsfr"][ok], -5, 4), xe, ye)
    X, Y = np.meshgrid(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), indexing="ij")
    fig = plt.figure(figsize=(12, 5.4))
    for k, (H, title, zl) in enumerate(((h, "Number x mass: where the stars are",
                                         "stellar mass per Mpc$^3$ (relative)"),
                                        (hs, "Number x SFR: where stars are being born",
                                         "star formation per Mpc$^3$ (relative)"))):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        Z = H / H.max()
        ax.plot_surface(X, Y, Z, cmap=DENSITY_CMAP, linewidth=0.2, edgecolor=ps.BLUE_RAMP[6],
                        antialiased=True, rstride=1, cstride=1, alpha=0.95)
        ax.view_init(elev=32, azim=-128)
        ax.set_xlabel("log M$_*$", labelpad=4)
        ax.set_ylabel("log SFR", labelpad=4)
        ax.set_zlabel(zl, labelpad=6, fontsize=7.5)
        ax.set_title(title, pad=2, loc="center")
        ax.set_facecolor(ps.SURFACE)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor(ps.SURFACE)
            axis.pane.set_edgecolor(ps.GRID)
            axis._axinfo["grid"]["color"] = ps.GRID
        ax.tick_params(labelsize=7)
    ps.caption(fig, "The galaxies of Fig. 01 in 0.2 x 0.2 dex cells, as in Renzini & Peng (2015, "
               "their Figs. 2 and 3), each weighted by 1/Vmax and by its stellar mass (left) or "
               "its SFR (right). Left: most of the stellar mass in today's universe sits in the "
               "quenched peak, in galaxies of ~10^10.5-10^11 solar masses. Right: almost all new "
               "stars form on the star-forming ridge, a few tenths of a dex lower in mass. Fig. 01 "
               "shows the plain number of galaxies.", y=0.02)
    return ps.save(fig, _out("02_sfr_mass_3d"))


# ---------------------------------------------------------------------------
# Quenching: mass and environment
# ---------------------------------------------------------------------------


def _spread_labels(ax, items, min_gap):
    """Line-end labels: (x, y, text, color) nudged apart vertically by min_gap."""
    items = sorted(items, key=lambda t: t[1])
    ys = [t[1] for t in items]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + min_gap)
    for (x, _, text, col), y in zip(items, ys):
        ax.text(x, y, text, fontsize=7, color=col, va="center")


def fig_quenching(res: dict) -> list:
    ms, msg, msm = res["ms"], res["ms_gswlc"], res["ms_gswlc_m2"]
    sep = res["env"]["quench"]
    fig, axs = plt.subplots(2, 2, figsize=(11.0, 8.4), gridspec_kw={"hspace": 0.38,
                                                                    "wspace": 0.25})
    # (a) quenched fraction vs mass, three methods
    a = axs[0, 0]
    for m_, color, ls, label in ((ms, ps.C_BLUE_GALAXY, "-", "MPA-JHU"),
                                 (msg, ps.ORANGE, "--", "GSWLC-2 (X)"),
                                 (msm, ps.MAGENTA, ":", "GSWLC-2 (M, deeper UV)")):
        q = m_["fq_vs_mass"]
        a.errorbar(q["logm"], q["fq"], yerr=q["err"], color=color, ls=ls, lw=1.6, marker="o",
                   ms=3.5, capsize=0, label=label)
    a.axhline(0.5, color=ps.MUTED, lw=0.8)
    a.set_xlabel("log stellar mass  (M$_\\odot$)")
    a.set_ylabel("fraction quenched (below the valley)")
    a.set_ylim(0, 1)
    a.set_title("(a) Massive galaxies have stopped", pad=6)
    a.legend(fontsize=7.5, loc="upper left")

    # (b) quenched fraction vs density, per mass bin
    b = axs[0, 1]
    dc = np.array(sep["dens_centres"])
    mc = np.array(sep["mass_centres"])
    cmap = LinearSegmentedColormap.from_list("mass", ps.BLUE_RAMP[2:])
    labels = []
    for i, m in enumerate(mc):
        col = cmap(i / max(len(mc) - 1, 1))
        good = np.isfinite(sep["f"][i])
        b.errorbar(dc[good], sep["f"][i][good], yerr=sep["err"][i][good], color=col, marker="o",
                   ms=3, lw=1.3, capsize=0)
        if good.any():
            j = np.flatnonzero(good)[-1]
            labels.append((dc[j] + 0.06, sep["f"][i][j], f"{m:.2f}", ps.INK_SECONDARY))
    _spread_labels(b, labels, 0.03)
    b.set_xlabel("log (1 + δ)   local overdensity (5th neighbour)")
    b.set_ylabel("fraction quenched")
    b.set_ylim(0, 1)
    b.set_title("(b) ...and so have crowded ones, at every mass", pad=6)
    b.text(0.02, 0.96, "labels: log M$_*$", transform=b.transAxes, fontsize=7.5,
           color=ps.INK_SECONDARY, va="top")

    # (c) eps_rho vs mass for each density bin (Peng Fig 5 bottom)
    c = axs[1, 0]
    dcmap = LinearSegmentedColormap.from_list("dens", ["#f3c9b5", ps.RED, "#7a1f1f"])
    clabels = []
    for j in range(1, len(dc)):
        col = dcmap((j - 1) / max(len(dc) - 2, 1))
        e = sep["eps_rho"][:, j]
        good = np.isfinite(e)
        if good.sum() < 2:
            continue
        c.plot(mc[good], e[good], color=col, marker="o", ms=3, lw=1.3)
        clabels.append((mc[good][-1] + 0.05, e[good][-1], f"{dc[j]:.2f}", ps.INK_SECONDARY))
        fit = sep["fit"]
        model = 1 - np.exp(-(10 ** (dc[j] - fit["log_p1"])) ** fit["p2"])
        c.axhline(model, color=col, lw=0.8, alpha=0.6)
    _spread_labels(c, clabels, 0.025)
    c.set_xlabel("log stellar mass  (M$_\\odot$)")
    c.set_ylabel("ε$_\\rho$  environment-quenching efficiency")
    c.set_title("(c) Environment's effect does not depend on mass", pad=6)
    c.text(0.02, 0.96, "labels: log (1+δ); thin lines: separable fit", transform=c.transAxes,
           fontsize=7.5, color=ps.INK_SECONDARY, va="top")

    # (d) eps_m vs density for each mass bin
    d = axs[1, 1]
    dlabels = []
    for i in range(1, len(mc)):
        col = cmap(i / max(len(mc) - 1, 1))
        e = sep["eps_m"][i, :]
        good = np.isfinite(e)
        if good.sum() < 2:
            continue
        d.plot(dc[good], e[good], color=col, marker="o", ms=3, lw=1.3)
        dlabels.append((dc[good][-1] + 0.06, e[good][-1], f"{mc[i]:.2f}", ps.INK_SECONDARY))
        fit = sep["fit"]
        model = 1 - np.exp(-(10 ** (mc[i] - fit["log_p3"])) ** fit["p4"])
        d.axhline(model, color=col, lw=0.8, alpha=0.6)
    _spread_labels(d, dlabels, 0.03)
    d.set_xlabel("log (1 + δ)")
    d.set_ylabel("ε$_m$  mass-quenching efficiency")
    d.set_title("(d) Mass's effect does not depend on environment", pad=6)
    d.text(0.02, 0.96, "labels: log M$_*$; thin lines: separable fit", transform=d.transAxes,
           fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    f = sep["fit"]
    pg = lit.PENG_2010
    ps.caption(fig, "Peng et al. (2010)'s test, redone with sSFR-defined quenching. "
               "ε_ρ = the fraction of galaxies that would be star-forming in the emptiest "
               "environments but are quenched in denser ones (same mass); ε_m the same for mass. "
               "Flat curves = the two effects are separable: f = ε_m + ε_ρ − ε_m ε_ρ. Our fit of "
               f"Peng's model: log p1 = {f['log_p1']:.2f}, p2 = {f['p2']:.2f}, log p3 = "
               f"{f['log_p3']:.2f}, p4 = {f['p4']:.2f} (Peng et al., colour-defined: "
               f"{pg['log_p1'][0]}, {pg['p2'][0]}, {pg['log_p3'][0]}, {pg['p4'][0]}).", y=0.0)
    return ps.save(fig, _out("03_quenching_mass_environment"))


# ---------------------------------------------------------------------------
# AGN
# ---------------------------------------------------------------------------


def fig_agn(res: dict) -> list:
    g, w, ms = res["_g"], res["_w_local"], res["ms"]
    ok = ms["ok"]
    fig, axs = plt.subplots(1, 3, figsize=(14.0, 4.6), gridspec_kw={"wspace": 0.42})
    # (a) BPT diagram
    a = axs[0]
    full = ok & np.isfinite(g["bpt_x"])
    xe = np.arange(-1.6, 0.61, 0.02)
    ye = np.arange(-1.2, 1.31, 0.02)
    h, _, _ = np.histogram2d(g["bpt_x"][full], g["bpt_y"][full], bins=[xe, ye])
    a.pcolormesh(xe, ye, h.T, cmap=DENSITY_CMAP, norm=LogNorm(1, h.max()), rasterized=True)
    xx = np.linspace(-1.6, 0.0, 200)
    a.plot(xx, gx.kauffmann03(xx), color=ps.INK, lw=1.2)
    xx2 = np.linspace(-1.6, 0.4, 200)
    a.plot(xx2, gx.kewley01(xx2), color=ps.INK, lw=1.2, ls="--")
    xx3 = np.linspace(-0.18, 0.6, 50)
    a.plot(xx3, gx.schawinski07(xx3), color=ps.INK_SECONDARY, lw=0.9, ls=":")
    a.text(-1.1, -0.7, "star-forming", fontsize=8, color=ps.C_BLUE_GALAXY, weight="bold")
    a.text(-0.12, -0.9, "composite", fontsize=8, color=ps.INK_SECONDARY, weight="bold")
    a.text(-0.6, 1.05, "Seyfert", fontsize=8, color=ps.C_AGN, weight="bold")
    a.text(0.2, -0.35, "LINER", fontsize=8, color=ps.C_AGN, weight="bold")
    a.set_xlim(-1.6, 0.6)
    a.set_ylim(-1.2, 1.3)
    a.set_xlabel("log [N II] 6584 / Hα")
    a.set_ylabel("log [O III] 5007 / Hβ")
    a.set_title("(a) What lights the gas", pad=6)

    # (b) AGN fraction across the plane
    b = axs[1]
    xe = np.arange(9.0, 11.81, 0.2)
    ye = np.arange(-3.0, 1.61, 0.2)
    agn = gx.is_agn(g["bpt"])
    hw = _plane_hist(g["logm"][ok], g["logsfr"][ok], w[ok], xe, ye)
    ha = _plane_hist(g["logm"][ok & agn], g["logsfr"][ok & agn], w[ok & agn], xe, ye)
    hn, _, _ = np.histogram2d(g["logm"][ok], g["logsfr"][ok], bins=[xe, ye])
    frac = np.where(hn >= 30, ha / np.where(hw > 0, hw, 1), np.nan)
    acmap = LinearSegmentedColormap.from_list("agn", ["#f1effa", "#b3a9e3", ps.VIOLET, "#231a5c"])
    im = b.pcolormesh(xe, ye, frac.T, cmap=acmap, vmin=0, vmax=0.6, rasterized=True)
    ridge = ms["ridge"]
    xx = np.linspace(9.0, 11.8, 20)
    b.plot(xx, gx.ms_sfr(xx, ridge), color=ps.C_BLUE_GALAXY, lw=1.6)
    v = ms["valley"]["position"]
    b.plot(xx, gx.ms_sfr(xx, ridge) + v, color="#8a5d00", lw=1.2, ls="--")
    b.text(9.05, gx.ms_sfr(9.05, ridge) + 0.15, "main sequence", fontsize=7.5,
           color=ps.C_BLUE_GALAXY)
    b.text(9.05, gx.ms_sfr(9.05, ridge) + v + 0.1, "valley floor", fontsize=7.5, color="#8a5d00")
    cb = fig.colorbar(im, ax=b, fraction=0.045, pad=0.02)
    cb.set_label("fraction with AGN-lit gas (composite + Seyfert + LINER)", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    b.set_xlim(9.0, 11.8)
    b.set_ylim(-3.0, 1.6)
    b.set_xlabel("log stellar mass  (M$_\\odot$)")
    b.set_ylabel("log SFR  (M$_\\odot$/yr)")
    b.set_title("(b) Where AGN live on the plane", pad=6)

    # (c) AGN fraction vs distance from MS, MPA vs GSWLC SFRs
    c = axs[2]
    mcol = {(10.0, 10.5): ps.BLUE_RAMP[4], (10.5, 11.0): ps.BLUE_RAMP[6],
            (11.0, 11.5): ps.BLUE_RAMP[8]}
    from matplotlib.lines import Line2D
    for row in res["agn"]["by_dms"]:
        key = tuple(row["logm"])
        if key in mcol:
            c.plot(row["dms"], row["f_agn"], color=mcol[key], lw=1.6, marker="o", ms=3)
    for row in res["agn_gswlc"]["by_dms"]:
        key = tuple(row["logm"])
        if key in mcol:
            c.plot(row["dms"], row["f_agn"], color=mcol[key], lw=1.2, ls="--")
    handles = [Line2D([], [], color=col, lw=1.8, label=f"log M {k[0]:g}–{k[1]:g}")
               for k, col in mcol.items()]
    handles += [Line2D([], [], color=ps.INK_SECONDARY, lw=1.6, marker="o", ms=3,
                       label="MPA-JHU SFRs"),
                Line2D([], [], color=ps.INK_SECONDARY, lw=1.2, ls="--", label="GSWLC-2 SFRs")]
    c.axvspan(v - 0.25, v + 0.25, color=ps.C_GREEN_VALLEY, alpha=0.22, lw=0)
    c.text(v, 0.02, "green\nvalley", ha="center", fontsize=7.5, color="#8a5d00")
    c.set_xlabel("distance from the main sequence (dex)")
    c.set_ylabel("fraction with AGN-lit gas")
    c.set_ylim(0, None)
    c.set_title("(c) AGN peak just below the sequence", pad=6)
    c.legend(handles=handles, fontsize=6.5, loc="upper left")
    ps.caption(fig, "Line ratios from MPA-JHU fluxes, errors scaled up by the team's repeat-"
               "observation factors; all four lines at S/N >= 3 for panel (a). Boundaries: "
               "Kauffmann et al. 2003 (solid), Kewley et al. 2001 (dashed), Schawinski et al. "
               "2007 (dotted). Galaxies whose lines are too weak to classify count as non-AGN, "
               "so the fractions are lower limits. Panel (c) repeats the measurement with GSWLC's "
               "UV-to-IR SFRs (dashed), which do not use the emission lines or the 4000 Å "
               "break that MPA-JHU relies on for AGN hosts.", y=-0.02)
    return ps.save(fig, _out("04_agn_on_the_plane"))


# ---------------------------------------------------------------------------
# Two methods for one universe
# ---------------------------------------------------------------------------


def fig_methods(res: dict) -> list:
    cmp_ = res["gswlc_cmp"]
    arr = cmp_["arrays"]
    fig, axs = plt.subplots(1, 4, figsize=(14.5, 3.9), gridspec_kw={"wspace": 0.3})
    groups = (("star_forming", arr["bpt"] == gx.BPT_SF, "star-forming lines"),
              ("agn", gx.is_agn(arr["bpt"]), "AGN-lit lines"),
              ("weak_or_none", (arr["bpt"] == gx.BPT_WEAK) | (arr["bpt"] == gx.BPT_NONE),
               "weak or no lines"))
    e = np.arange(-13.5, -8.49, 0.05)
    for ax, (key, m, title) in zip(axs[:3], groups):
        h, _, _ = np.histogram2d(arr["ssfr_gswlc"][m], arr["ssfr_mpa"][m], bins=[e, e])
        ax.pcolormesh(e, e, h.T, cmap=DENSITY_CMAP, norm=LogNorm(1, max(h.max(), 2)),
                      rasterized=True)
        ax.plot([-13.5, -8.5], [-13.5, -8.5], color=ps.INK, lw=0.8)
        st = cmp_[f"dssfr_{key}"]
        ax.text(0.03, 0.97, f"{st['n']:,} galaxies\nmedian offset {st['median']:+.2f} dex\n"
                f"scatter {st['robust_sigma']:.2f} dex", transform=ax.transAxes, va="top",
                fontsize=7.5, color=ps.INK_SECONDARY)
        ax.set_xlim(-13.5, -8.5)
        ax.set_ylim(-13.5, -8.5)
        ax.set_xlabel("log sSFR, GSWLC-2 (UV-to-IR SED fit)")
        ax.set_title(title, pad=6)
    axs[0].set_ylabel("log sSFR, MPA-JHU (lines / D4000 + optical)")
    d = axs[3]
    dv = cmp_["dlogm_vs_ssfr"]
    d.plot(dv["ssfr"], dv["median"], color=ps.C_BLUE_GALAXY, marker="o", ms=4)
    lo, hi = lit.SALIM_2016["mass_offset_range_dex"]
    d.axhspan(-hi, -lo, color=ps.VIOLET, alpha=0.12, lw=0)
    d.text(-12.9, -0.5 * (lo + hi), "Salim et al. 2016:\nGSWLC 0.03–0.13 dex heavier",
           fontsize=7, color=ps.VIOLET, va="center")
    d.axhline(0, color=ps.MUTED, lw=0.8)
    d.set_xlabel("log sSFR (GSWLC)")
    d.set_ylabel("log M (MPA-JHU) − log M (GSWLC)")
    d.set_title("stellar masses", pad=6)
    ps.caption(fig, "The same galaxies measured two ways. As in Salim et al. (2016), MPA-JHU values "
                "are moved from the Kroupa to the Chabrier IMF (−0.025 dex) and both sSFRs use "
                "the GSWLC mass, so only the SFRs differ. Where the lines trace young stars the "
                "methods agree; for AGN hosts and line-less galaxies MPA-JHU has no direct "
                "tracer and the two disagree by up to ~1 dex: the green valley's population "
                "depends on which you believe.", y=-0.03)
    return ps.save(fig, _out("05_two_methods"))


# ---------------------------------------------------------------------------
# Who forms the stars today
# ---------------------------------------------------------------------------


def fig_census(res: dict) -> list:
    c = res["census_mean"]
    cl = res["census_local"]
    clm = res["census_mean_local"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.0, 4.3), gridspec_kw={"wspace": 0.28})
    bm = c["by_mass"]
    x = np.array(bm["logm"])
    sfr = np.array(bm["rho_sfr"]) / 0.2
    mass = np.array(bm["rho_star"]) / 0.2
    a.bar(x, sfr / sfr.sum() / 0.2 * 0.2 * 5, width=0.18, color=ps.C_BLUE_GALAXY,
          label="where stars are forming")
    a.bar(x, -(mass / mass.sum()) * 5, width=0.18, color=ps.C_RED_GALAXY,
          label="where stars already are")
    a.axhline(0, color=ps.INK_SECONDARY, lw=0.8)
    a.axvline(c["logm_half_sfr"], color=ps.C_BLUE_GALAXY, lw=0.9, ls="--")
    a.axvline(c["logm_half_mass"], color=ps.C_RED_GALAXY, lw=0.9, ls="--")
    a.text(8.05, 0.92 * a.get_ylim()[1],
           f"half of new stars\nform in galaxies\nbelow 10$^{{{c['logm_half_sfr']:.1f}}}$ M$_\\odot$",
           ha="left", fontsize=7.5, color=ps.C_BLUE_GALAXY, va="top")
    a.text(c["logm_half_mass"] + 0.55, 0.62 * a.get_ylim()[0],
           f"half of all stars\nlive above\n10$^{{{c['logm_half_mass']:.1f}}}$ M$_\\odot$",
           ha="left", fontsize=7.5, color=ps.C_RED_GALAXY, va="center")
    a.set_xlim(8.0, 12.2)
    yt = a.get_yticks()
    a.set_yticks(yt)
    a.set_yticklabels([f"{abs(t) / 5 * 100:.0f}%" for t in yt])
    a.set_xlabel("log stellar mass of the host galaxy (M$_\\odot$)")
    a.set_ylabel("share per 0.2 dex of mass")
    a.set_title("Today's star formation vs today's stars", pad=6)
    a.legend(fontsize=7.5, loc="lower left")

    # (b) shares by zone
    z = {k: {"sfr_share": clm["zones"][k]["sfr_share"], "mass_share": cl["zones"][k]["mass_share"]}
         for k in cl["zones"]}
    names = ["star_forming", "green_valley", "quenched"]
    labels = ["main sequence\n& above", "green valley", "quenched"]
    colors = [ps.C_BLUE_GALAXY, ps.C_GREEN_VALLEY, ps.C_RED_GALAXY]
    yy = np.arange(3)
    for i, nm in enumerate(names):
        b.barh(i + 0.18, 100 * z[nm]["sfr_share"], height=0.34, color=colors[i])
        b.barh(i - 0.18, 100 * z[nm]["mass_share"], height=0.34, color=colors[i], alpha=0.45)
        b.text(100 * z[nm]["sfr_share"] + 1, i + 0.18, f"{100 * z[nm]['sfr_share']:.0f}% of "
               "star formation", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
        b.text(100 * z[nm]["mass_share"] + 1, i - 0.18, f"{100 * z[nm]['mass_share']:.0f}% of "
               "stellar mass", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    b.set_yticks(yy)
    b.set_yticklabels(labels)
    b.invert_yaxis()
    b.set_xlim(0, 115)
    b.set_xlabel("share (%)")
    b.set_title("The quenched hold the stars; the sequence makes new ones", pad=6)
    b.grid(axis="y", visible=False)
    ps.caption(fig, f"1/Vmax census of SDSS galaxies (MPA-JHU, 0.02 < z < 0.2 left; "
               f"0.02 < z < 0.085 right, where every zone is well sampled). SFRs are the means of "
               f"each galaxy's likelihood. Total: {c['rho_sfr']:.4f} M_sun/yr/Mpc^3 of star "
               f"formation and {c['rho_star'] / 1e8:.2f} x 10^8 M_sun/Mpc^3 of stars (Kroupa IMF, "
               "H0 = 70).", y=-0.03)
    return ps.save(fig, _out("08_who_forms_stars_today"))


# ---------------------------------------------------------------------------
# Gallery: what the life cycle looks like
# ---------------------------------------------------------------------------

GALLERY_MASSES = (9.6, 10.1, 10.6, 11.1, 11.5)
GALLERY_ROWS = (("above the sequence", 0.45), ("on the main sequence", 0.0),
                ("green valley", None), ("quenched", -2.1))
GALLERY_FRAME_KPC = 50.0
GALLERY_PIXELS = 160


def gallery_selection(res: dict) -> list:
    """One galaxy per (mass, row) cell, chosen by rule, not by eye.

    In each cell (+/-0.12 dex in mass, +/-0.15 dex in distance from the main
    sequence), among galaxies at 0.02 < z < 0.05 (near enough to be resolved),
    take the one whose rest-frame g - r colour is closest to the cell's median:
    a typical member rather than an extreme. Cells with no such galaxy are
    searched again out to z = 0.085, the edge of the sample.
    """
    g, ms = res["_g"], res["ms"]
    dms = ms["dms"]
    v = ms["valley"]["position"]
    base = ms["ok"] & np.isfinite(g["gr0"])
    picks = []
    for row, (name, off) in enumerate(GALLERY_ROWS):
        target = v if off is None else off
        for col, m in enumerate(GALLERY_MASSES):
            cell = base & (np.abs(g["logm"] - m) < 0.12) & (np.abs(dms - target) < 0.15)
            near = cell & (g["z"] < 0.05)
            if not near.any():
                near = cell
            if not near.any():
                picks.append(None)
                continue
            idx = np.flatnonzero(near)
            i = idx[np.argmin(np.abs(g["gr0"][idx] - np.median(g["gr0"][idx])))]
            picks.append({"row": row, "col": col, "i": int(i), "objID": int(g["objID"][i]),
                          "ra": float(g["ra"][i]), "dec": float(g["dec"][i]),
                          "z": float(g["z"][i]), "logm": float(g["logm"][i]),
                          "logsfr": float(g["logsfr"][i]), "dms": float(dms[i]),
                          "bpt": gx.BPT_NAMES[int(g["bpt"][i])]})
    return picks


def fetch_cutout(p: dict):
    """SDSS DR17 colour cutout (g, r, i composite by SkyServer), cached as JPEG."""
    import io

    import requests

    import cosmo_physics as cp
    from config import CUTOUT_DIR, SKYSERVER_CUTOUT_URL

    path = CUTOUT_DIR / f"{p['objID']}_{int(GALLERY_FRAME_KPC)}kpc.jpg"
    if not path.exists():
        da_kpc = float(cp.angular_diameter_distance(p["z"])) * 1000.0
        size_arcsec = GALLERY_FRAME_KPC / da_kpc * 206264.806
        scale = size_arcsec / GALLERY_PIXELS
        for attempt in range(4):
            try:
                r = requests.get(SKYSERVER_CUTOUT_URL, timeout=(15, 120),
                                 params={"ra": p["ra"], "dec": p["dec"], "scale": f"{scale:.4f}",
                                         "width": GALLERY_PIXELS, "height": GALLERY_PIXELS})
                r.raise_for_status()
                path.write_bytes(r.content)
                break
            except requests.RequestException:
                import time
                time.sleep(5 * (attempt + 1))
    return plt.imread(io.BytesIO(path.read_bytes()), format="jpg")


def fig_gallery(res: dict) -> list:
    picks = gallery_selection(res)
    nr, nc = len(GALLERY_ROWS), len(GALLERY_MASSES)
    fig, axs = plt.subplots(nr, nc, figsize=(1.9 * nc + 1.4, 1.95 * nr + 0.9),
                            gridspec_kw={"wspace": 0.04, "hspace": 0.16})
    row_colors = [ps.C_BLUE_GALAXY, ps.C_BLUE_GALAXY, "#8a5d00", ps.C_RED_GALAXY]
    for p, ax in zip(picks, axs.ravel()):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
        if p is None:
            ax.set_facecolor(ps.SURFACE)
            ax.text(0.5, 0.5, "none in cell", ha="center", va="center", fontsize=7,
                    color=ps.MUTED, transform=ax.transAxes)
            continue
        ax.imshow(fetch_cutout(p))
        ax.text(0.03, 0.03, f"z={p['z']:.3f}  SFR={10 ** p['logsfr']:.2g}",
                transform=ax.transAxes, fontsize=6, color="white", va="bottom")
    for r, (name, _) in enumerate(GALLERY_ROWS):
        axs[r, 0].set_ylabel(name, fontsize=9, color=row_colors[r], weight="bold",
                             rotation=90, labelpad=6)
    for c, m in enumerate(GALLERY_MASSES):
        axs[0, c].set_title(f"10$^{{{m:g}}}$ M$_\\odot$", fontsize=9, loc="center", pad=4)
    fig.suptitle("The life cycle in pictures: SDSS images across the SFR–mass plane",
                 x=0.12, ha="left", y=0.995)
    ps.caption(fig, f"Every frame shows the same {GALLERY_FRAME_KPC:.0f} kpc (~160,000 light-years) "
               "at the galaxy's distance, so apparent sizes are true sizes. Columns: stellar mass. "
               "Rows: distance from the main sequence (+0.45 dex; 0; the valley floor; −2.1 dex). "
               "In each cell the galaxy at 0.02 < z < 0.05 (or further, if none) whose colour is "
               "closest to the cell's median is shown: picked by rule, not by eye. SFR in solar "
               "masses per year. Images: SDSS DR17 SkyServer g, r, i composites.",
               y=0.03)
    return ps.save(fig, _out("09_gallery"), formats=("png",))


# ---------------------------------------------------------------------------
# Where the quenched galaxies live: a slice of the cosmic web
# ---------------------------------------------------------------------------


def fig_web_slice(res: dict) -> list:
    """D2's wedge, coloured by star-formation status."""
    import cosmo_physics as cp

    g, ms = res["_g"], res["ms"]
    ok = ms["ok"] & (np.abs(g["dec"]) < 2.0) & (g["ra"] > 128) & (g["ra"] < 252)
    ra_mid = 190.0
    ang = np.radians(g["ra"] - ra_mid)
    d = g["d_c"]
    x, y = d * np.sin(ang), d * np.cos(ang)
    fig, ax = plt.subplots(figsize=(11.5, 6.6))
    groups = ((ms["star_forming"], ps.C_BLUE_GALAXY, "star-forming"),
              (ms["green_valley"], ps.C_GREEN_VALLEY, "green valley"),
              (ms["quenched"], ps.C_RED_GALAXY, "quenched"))
    # one scatter call in random order, so no colour is drawn on top of another
    idx = np.flatnonzero(ok & (ms["star_forming"] | ms["green_valley"] | ms["quenched"]))
    rng = np.random.default_rng(1)
    rng.shuffle(idx)
    cols = np.where(ms["quenched"][idx], ps.C_RED_GALAXY,
                    np.where(ms["green_valley"][idx], ps.C_GREEN_VALLEY, ps.C_BLUE_GALAXY))
    size = 0.6 + 2.2 * np.clip(g["logm"][idx] - 9.0, 0, 3) ** 1.5
    ax.scatter(x[idx], y[idx], s=size, c=list(cols), lw=0, alpha=0.75, rasterized=True)
    for flag, color, label in groups:
        ax.scatter([], [], s=20, color=color, label=f"{label} ({int((ok & flag).sum()):,})")
    dmax = float(cp.comoving_distance(0.085))
    dmin = float(cp.comoving_distance(0.02))
    for dd in (dmin, dmax):
        tt = np.radians(np.linspace(128 - ra_mid, 252 - ra_mid, 200))
        ax.plot(dd * np.sin(tt), dd * np.cos(tt), color=ps.AXIS, lw=0.8)
    for r_ in (128, 160, 190, 220, 252):
        t_ = np.radians(r_ - ra_mid)
        ax.plot([dmin * np.sin(t_), dmax * np.sin(t_)], [dmin * np.cos(t_), dmax * np.cos(t_)],
                color=ps.GRID, lw=0.6)
        ax.text(1.03 * dmax * np.sin(t_), 1.03 * dmax * np.cos(t_), f"RA {r_}°", fontsize=7,
                color=ps.MUTED, ha="center", va="bottom")
    for zz in (0.02, 0.05, 0.085):
        dd = float(cp.comoving_distance(zz))
        t_ = np.radians(252 - ra_mid)
        ax.text(dd * np.sin(t_) + 6, dd * np.cos(t_), f"z = {zz}", fontsize=7, color=ps.MUTED,
                va="center")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Quenched galaxies crowd into clusters: a 4° slice of the nearby universe",
                 loc="left")
    ax.legend(loc="lower left", fontsize=8, markerscale=1)
    ps.caption(fig, "Every galaxy of the 0.02 < z < 0.085 sample within 2° of the celestial equator "
               "(RA 128°–252°), placed at the distance its redshift implies (Milky Way at the "
               f"bottom; the outer arc is {dmax:.0f} Mpc away). Dot size grows with stellar mass; "
               "galaxies are drawn in random order. "
               "Clusters appear as radial spikes ('fingers of God', D2) and are dominated by red, "
               "quenched galaxies; star-forming galaxies trace the filaments and the void edges.",
               y=0.06)
    return ps.save(fig, _out("10_quenched_galaxies_in_the_web"))
