"""Every figure of the project, from the analysis results.

Conventions (plotstyle.py): the prediction of an expanding universe is blue
and solid, tired light is orange and dashed, measurements are ink, told apart
by marker shape. Redshift, a magnitude, uses one blue hue stepped light to
dark.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plotstyle as ps
import tcmb
import tolman
from config import FIGURE_DIR, T0_FIRAS, TOLMAN_CLUSTERS

ps.use_style()

T0 = T0_FIRAS[0]


def _binned(x, y, edges, stat=np.median):
    idx = np.digitize(x, edges) - 1
    xc, yc, ye, n = [], [], [], []
    for k in range(len(edges) - 1):
        m = idx == k
        if m.sum() < 5:
            continue
        xc.append(np.median(x[m]))
        yc.append(stat(y[m]))
        ye.append(1.2533 * np.std(y[m]) / np.sqrt(m.sum()))
        n.append(m.sum())
    return np.array(xc), np.array(yc), np.array(ye), np.array(n)


def _prediction_legend(ax, loc="upper left", extra=()):
    handles = [plt.Line2D([], [], color=ps.C_EXPANSION, ls=ps.LS_EXPANSION, lw=1.8,
                          label="expanding universe"),
               plt.Line2D([], [], color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8,
                          label="tired light (static universe)")]
    ax.legend(handles=list(extra) + handles, loc=loc)


# ---------------------------------------------------------------------------
# Test 1: supernova time dilation
# ---------------------------------------------------------------------------

def fig_slow_motion(td, path):
    """Light curves at the same rest-frame colour, in observed and in de-dilated time."""
    p = td["template"]["fit_points"].copy()
    lam = 10 ** p["loglam"]
    p["norm"] = p["flux"] / p["peak"]
    groups = [("r", (0.27, 0.45)), ("i", (0.56, 0.76)), ("z", (0.84, 1.06))]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for (band, (z0, z1)), col in zip(groups, ps.Z_RAMP):
        g = p[(p["band"] == band) & (p["z"] >= z0) & (p["z"] < z1) & (lam > 440) & (lam < 510)]
        nsn = g["cid"].nunique()
        zmed = g.groupby("cid")["z"].first().median()
        for ax, x in ((axes[0], g["t"].to_numpy()),
                      (axes[1], (g["t"] / (1 + g["z"])).to_numpy())):
            edges = np.arange(-30, 90, 3.0) if ax is axes[0] else np.arange(-30, 90, 2.0)
            xc, yc, ye, _ = _binned(x, g["norm"].to_numpy(), edges)
            ax.fill_between(xc, yc - ye, yc + ye, color=col, alpha=0.18, lw=0)
            ax.plot(xc, yc, color=col, lw=2,
                    label=f"{band} band, z = {z0:.2f}-{z1:.2f} ({nsn} supernovae)")
    tmpl = td["template"]["template"]
    j = np.argmin(np.abs(tmpl["lam"] - 470))
    S = tmpl["S"][:, j]
    axes[1].plot(tmpl["tau"], S / S.max(), color=ps.INK, lw=1.0, ls=":",
                 label="fitted template at 470 nm")
    axes[0].set_xlim(-25, 85)
    axes[1].set_xlim(-20, 60)
    axes[0].set_xlabel("days since peak, as observed on Earth")
    axes[1].set_xlabel("days since peak divided by (1 + z)")
    axes[0].set_ylabel("brightness / peak brightness (median)")
    axes[0].set_title("As observed: the farther, the slower", loc="left")
    axes[1].set_title("Time divided by (1 + z): one light curve", loc="left")
    axes[0].legend(loc="upper right", fontsize=7.5)
    axes[1].legend(loc="upper right", fontsize=7.5)
    for ax in axes:
        ax.axhline(0, color=ps.AXIS, lw=0.8)
    fig.suptitle("Type Ia supernovae seen in the same rest-frame light (440-510 nm), "
                 "sorted by distance", x=0.06, ha="left")
    ps.caption(fig, "Dark Energy Survey 5-year sample. Each band is chosen so that it samples blue "
                    "rest-frame light (~470 nm) at its redshift: r at z ~ 0.36, i at z ~ 0.66, "
                    "z at z ~ 0.95. Fluxes are divided by each light curve's own observed peak; lines "
                    "are medians in time bins with their uncertainties.")
    fig.tight_layout()
    return ps.save(fig, path)


def fig_widths(td, path):
    """White et al.'s method: light-curve width against 1 + z."""
    w = td["white"]
    ours = w["averaged_table"]
    theirs = w["white_average_table"]
    null = w["null_table"]
    fit = w["average"]["width"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    ax.scatter(1 + ours["z"], ours["width"], s=4, color=ps.C_DATA_SOFT, alpha=0.35, lw=0,
               rasterized=True, label=f"each supernova ({len(ours)}), ours")
    edges = np.quantile(ours["z"], np.linspace(0, 1, 15))
    xc, yc, ye, _ = _binned(ours["z"].to_numpy(), ours["width"].to_numpy(), edges)
    ax.errorbar(1 + xc, yc, ye, fmt="o", color=ps.C_DATA, ms=6, lw=1.2,
                label="ours, binned (median)")
    xc2, yc2, ye2, _ = _binned(theirs["z"].to_numpy(), theirs["width"].to_numpy(), edges)
    ax.errorbar(1 + xc2 + 0.006, yc2, ye2, fmt="s", color=ps.C_THIRD, ms=5, lw=1.0, mfc="none",
                label="White et al. 2024's own widths, binned")
    x = np.linspace(1.0, 2.25, 50)
    ax.plot(x, x, color=ps.C_EXPANSION, lw=1.8)
    ax.plot(x, np.ones_like(x), color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    ax.set_xlabel("1 + z")
    ax.set_ylabel("light-curve width relative to the reference")
    ax.set_ylim(0.5, 3.0)
    ax.set_xlim(1.0, 2.25)
    ax.set_title(r"Widths grow as $(1+z)^b$" + f" with b = {fit['b']:.3f} ± {fit['bootstrap_err']:.3f}",
                 loc="left")
    handles, labels = ax.get_legend_handles_labels()
    _prediction_legend(ax, extra=handles)

    ax = axes[1]
    xc, yc, ye, _ = _binned(null["z"].to_numpy(), null["width"].to_numpy(),
                            np.quantile(null["z"], np.linspace(0, 1, 12)))
    ax.scatter(1 + null["z"], null["width"], s=4, color=ps.C_DATA_SOFT, alpha=0.3, lw=0,
               rasterized=True)
    ax.errorbar(1 + xc, yc, ye, fmt="o", color=ps.C_DATA, ms=5, lw=1.2)
    ax.axhline(1.0, color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    nt = w["null_test"]
    ax.set_title("Null test: reference curves not de-dilated", loc="left")
    ax.set_xlabel("1 + z")
    ax.set_ylabel("width relative to an undilated reference")
    ax.set_ylim(0.3, 2.2)
    lin = nt["linear"]
    ps.annotate(ax, f"still rises: slope {lin['slope']:.2f} ± {lin['slope_err']:.2f}\n"
                    "per unit of (1 + z); tired light: flat", (1.75, 1.2), (1.12, 1.85))
    ps.caption(fig, "Replication of White et al. (2024) on the DES-SN5YR light curves: each "
                    "supernova's light curve is compared with a reference built from other supernovae "
                    "seen at the same rest-frame wavelength, and stretched in time until it matches. "
                    "Widths are averaged over the g, r, i, z bands. The error on b is from bootstrapping "
                    "supernovae.")
    fig.tight_layout()
    return ps.save(fig, path)


def fig_b_measurements(td, path):
    """Our global fit's chi-square against b, and the scatter-minimisation method."""
    t = td["template"]
    w = td["white"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    sc = t["scan"]
    d = sc["chi2"] - t["chi2_best"]
    ax.plot(sc["b"], d, "o-", color=ps.C_DATA, ms=5, lw=1.4)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylim(-1, 2e5)
    ax.axvline(1.0, color=ps.C_EXPANSION, lw=1.8)
    ax.axvline(0.0, color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    ax.set_xlabel(r"time-dilation power b in $(1+z)^b$")
    ax.set_ylabel("Δχ² against the best fit")
    ax.set_title("One template for all supernovae: χ² of each b", loc="left")
    ps.annotate(ax, f"tired light: Δχ² = {d[0]:,.0f}", (0.0, d[0]), (0.12, 6e4))
    ps.annotate(ax, f"best: b = {t['b']:.3f} ± {t['b_err_bootstrap']:.3f}\n(bootstrap)",
                (t["b"], 0.5), (1.12, 4))
    _prediction_legend(ax, loc="upper right")

    ax = axes[1]
    disp = w["dispersion"]
    from timedilation import B_DISPERSION
    ours = disp["riz_mean_curve"]
    theirs = disp["white_riz_mean_curve"]
    ax.plot(B_DISPERSION, ours / ours.min(), "o-", color=ps.C_DATA, ms=4, lw=1.4,
            label=f"ours (minimum at b = {disp['riz_b_min']:.2f})")
    ax.plot(B_DISPERSION, theirs / theirs.min(), "s--", color=ps.C_THIRD, ms=4, lw=1.0, mfc="none",
            label=f"White et al.'s curves (b = {disp['white_riz_b_min']:.2f})")
    ax.axvline(1.0, color=ps.C_EXPANSION, lw=1.8)
    ax.axvline(0.0, color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    ax.set_xlabel(r"b used to scale every reference point's time by $1/(1+z)^b$")
    ax.set_ylabel("scatter of the reference light curves (relative)")
    ax.set_title("No light-curve model at all: scatter is smallest at b ≈ 1", loc="left")
    ax.legend(loc="upper center")
    ps.caption(fig, "Left: our own measurement. A smooth light-curve surface in (rest-frame time, "
                    "rest-frame wavelength) is fitted to all 1,542 supernovae at once with the time "
                    "scaled by (1 + z)^b, for each b; the observer-frame data window is the same for "
                    "every b. Right: White et al.'s first method (median flux scatter in time bins of each "
                    "reference curve; r, i, z bands averaged).")
    fig.tight_layout()
    return ps.save(fig, path)


def fig_dilation_bins(td, path):
    """One dilation factor per redshift bin, no power law assumed."""
    b = td["template"]["bins"]
    ref = b["ref"]
    x = b["one_plus_z"] / b["one_plus_z"][ref]
    y, e = b["D_ratio"], b["D_ratio_err"]
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    xx = np.linspace(0.65, 1.5, 10)
    ax.plot(xx, xx, color=ps.C_EXPANSION, lw=1.8)
    ax.plot(xx, np.ones_like(xx), color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    keep = np.arange(len(x)) != ref
    ax.errorbar(x[keep], y[keep], e[keep], fmt="o", color=ps.C_DATA, ms=6, lw=1.2,
                label="measured (each redshift bin)")
    ax.plot(x[ref], y[ref], "D", color=ps.C_DATA, mfc=ps.SURFACE, ms=7,
            label="reference bin (z = 0.45-0.55)")
    for xi, yi, n in zip(x, y, b["n_sne"]):
        ax.annotate(f"{n}", (xi, yi), xytext=(6, -12), textcoords="offset points",
                    fontsize=7, color=ps.MUTED)
    ax.set_xlabel("(1 + z) of the bin / (1 + z) of the reference bin")
    ax.set_ylabel("measured time stretch / time stretch of the reference bin")
    ax.set_title("How much slower each redshift bin runs", loc="left")
    handles, _ = ax.get_legend_handles_labels()
    _prediction_legend(ax, extra=handles)
    ps.caption(fig, "Global template fit with one free time-scale factor per redshift bin instead "
                    "of a power law. Only ratios are measurable (the template's own time axis absorbs "
                    f"an overall factor). Error bars are the fit's, enlarged {b['error_scale']:.1f}x as "
                    "the bootstrap enlarged the error of b, to include supernova-to-supernova "
                    "differences. Grey numbers: supernovae per bin.")
    fig.tight_layout()
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Test 2: the CMB thermometer
# ---------------------------------------------------------------------------

MARKERS = [
    ("Planck clusters (Hurier+14)", lambda d: d["object"].str.startswith("Planck clusters"), "o"),
    ("SPT clusters (Saro+14)", lambda d: d["object"].str.startswith("SPT"), "s"),
    ("molecules, radio (Kotani+25, +26)", lambda d: d["method"].str.startswith("radio"), "D"),
    ("CO in quasar absorbers (Klimenko+20)", lambda d: d["method"].str.startswith("CO"), "^"),
    ("C I / C II in quasar absorbers", lambda d: d["method"].str.startswith(("CI", "CII")), "v"),
    ("water at z = 6.34 (Riechers+22)", lambda d: d["method"].str.startswith("H2O"), "*"),
]


def _tcmb_points(ax, base, yfun):
    for label, sel, m in MARKERS:
        d = base[sel(base) & (base["kind"] != "upper")]
        if d.empty:
            continue
        x = 1 + d["z"].to_numpy()
        y, lo, hi = yfun(d)
        ax.errorbar(x, y, [lo, hi], fmt=m, color=ps.C_DATA, mfc=ps.SURFACE if m in "Dv" else ps.C_DATA,
                    ms=7 if m == "*" else 5, lw=1.0, capsize=0, label=label)
    lim = base[base["kind"] == "upper"]
    x = 1 + lim["z"].to_numpy()
    y = yfun(lim.assign(err_plus=0, err_minus=0, err_sys=0))[0]
    ax.plot(x, y, marker=r"$\downarrow$", ls="none", color=ps.MUTED, ms=11,
            label="upper limits (not fitted)")


def fig_tcmb(tc, path):
    base = tcmb.select(tc["data"], "baseline")
    base = base[base["z"] > 0]
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 8.4), sharex=True,
                             gridspec_kw={"height_ratios": [1.5, 1]})
    x = np.geomspace(1, 8, 200)
    for ax, scale in zip(axes, ("T", "T/(1+z)")):
        div = (lambda d: 1.0) if scale == "T" else (lambda d: 1 + d["z"].to_numpy())

        def yfun(d, div=div):
            e_lo = np.hypot(d["err_minus"], d["err_sys"]).to_numpy()
            e_hi = np.hypot(d["err_plus"], d["err_sys"]).to_numpy()
            k = div(d)
            return d["T"].to_numpy() / k, e_lo / k, e_hi / k

        if scale == "T":
            ax.plot(x, T0 * x, color=ps.C_EXPANSION, lw=1.8)
            ax.plot(x, T0 * np.ones_like(x), color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
        else:
            ax.plot(x, T0 * np.ones_like(x), color=ps.C_EXPANSION, lw=1.8)
            ax.plot(x, T0 / x, color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
        _tcmb_points(ax, base, yfun)
        ax.plot([1], [T0], "*", color=ps.C_DATA, ms=11)
        ax.set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_yticks([3, 4, 6, 10, 20, 30])
    axes[0].set_yticklabels(["3", "4", "6", "10", "20", "30"])
    axes[0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[0].set_ylim(2.3, 36)
    for ax in axes:
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[0].set_ylabel("CMB temperature at that redshift (K)")
    axes[1].set_ylabel("temperature / (1 + z)  (K)")
    axes[1].set_ylim(0.5, 4.5)
    axes[1].set_xlabel("1 + z")
    axes[1].set_xticks([1, 1.5, 2, 3, 4, 5, 7])
    axes[1].set_xticklabels(["1", "1.5", "2", "3", "4", "5", "7"])
    b = tc["baseline"]
    axes[0].set_title(r"The CMB was hotter in the past: $T = T_0\,(1+z)^{1-\beta}$, "
                      f"β = {b['beta']:.3f} ± {max(b['plus'], b['minus']):.3f}", loc="left")
    handles, _ = axes[0].get_legend_handles_labels()
    _prediction_legend(axes[0], extra=handles)
    ps.annotate(axes[0], "COBE/FIRAS, today", (1.0, T0), (1.08, 2.45))
    m = tc["baseline_models"]
    ps.annotate(axes[1], f"tired light: χ² = {m['chi2_tired']:,.0f} for {m['n']} measurements\n"
                         f"expansion: χ² = {m['chi2_expansion']:.0f}", (3.0, T0 / 3.0), (3.2, 1.2))
    ps.caption(fig, "Published measurements (literature/tcmb_measurements.csv): Sunyaev-Zel'dovich "
                    "spectra of galaxy clusters, excitation of molecules and atoms in gas clouds in front "
                    "of quasars, and water absorbing the CMB at z = 6.34. Absorber values are corrected for "
                    "collisional excitation where available (Klimenko et al. 2020). Error bars include quoted "
                    "systematics.")
    fig.tight_layout()
    return ps.save(fig, path)


def fig_sz_null(tc, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    nu = np.linspace(20, 600, 800)
    z = 0.5
    x_exp = tcmb.H_PLANCK * nu * 1e9 / (tcmb.K_BOLTZ * T0)
    x_tl = x_exp * (1 + z)
    g_exp = tcmb.sz_spectrum(x_exp)
    g_tl = tcmb.sz_spectrum(x_tl)
    ax.plot(nu, g_exp / np.abs(g_exp).max(), color=ps.C_EXPANSION, lw=1.8)
    ax.plot(nu, g_tl / np.abs(g_tl).max(), color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    for f in (100, 143, 217, 353, 545):
        ax.axvline(f, color=ps.GRID, lw=3, zorder=0)
        ax.annotate(f"{f}", (f, -0.9), ha="center", fontsize=7, color=ps.MUTED)
    ax.annotate("Planck bands (GHz)", (390, -0.8), fontsize=7, color=ps.MUTED)
    ax.set_ylim(-0.98, 1.1)
    ax.set_xlabel("observed frequency (GHz)")
    ax.set_ylabel("change of the CMB seen through the cluster (normalised)")
    ax.set_title("A cluster at z = 0.5: shadow below the null, glow above", loc="left")
    _prediction_legend(ax, loc="upper left")
    ps.annotate(ax, "null 217 GHz", (217.5, 0), (240, 0.45), color=ps.C_EXPANSION)
    ps.annotate(ax, "null 145 GHz", (145.0, 0), (60, 0.45), color=ps.C_TIRED)

    ax = axes[1]
    d = tc["sz_null"]
    for label, sel, m in MARKERS[:2]:
        s = d[sel(d)]
        ax.errorbar(1 + s["z"], s["nu_null"], [s["nu_err_minus"], s["nu_err_plus"]], fmt=m,
                    color=ps.C_DATA, ms=5, lw=1.0, mfc=ps.C_DATA if m == "o" else ps.SURFACE,
                    label=label)
    xx = np.linspace(1, 2.05, 50)
    nu0 = tcmb.null_frequency_ghz(T0, 0)
    ax.plot(xx, nu0 * np.ones_like(xx), color=ps.C_EXPANSION, lw=1.8)
    ax.plot(xx, nu0 / xx, color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    ax.set_xlabel("1 + z of the clusters")
    ax.set_ylabel("observed frequency of the SZ null (GHz)")
    ax.set_title("Measured: the null stays at 217 GHz", loc="left")
    handles, _ = ax.get_legend_handles_labels()
    _prediction_legend(ax, loc="lower left", extra=handles)
    ps.caption(fig, "Left: non-relativistic thermal SZ spectrum, x⁴eˣ/(eˣ−1)² (x coth(x/2) − 4) with "
                    "x = hν/kT at the cluster. If the CMB there is (1 + z) hotter, the redshift of the "
                    "light cancels it and the null is observed at 217 GHz for every cluster. Right: each "
                    "published cluster-stack temperature re-expressed as the observed null frequency "
                    "(relativistic corrections, which shift it by a few GHz for the hottest clusters, are "
                    "included in the published analyses).")
    fig.tight_layout()
    return ps.save(fig, path)


def fig_co_thermometer(tc, path):
    co = tc["co_excitation"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    keep = co["n_h2"] <= 500
    for (z, ex), col in zip(sorted(co["excess"].items()), ps.Z_RAMP):
        ax.plot(co["n_h2"][keep], ex[keep], color=col, lw=2,
                label=f"z = {z:.2f}: CMB at {T0 * (1 + z):.2f} K")
    ax.axvspan(8, 200, color=ps.GRID, alpha=0.6, lw=0, zorder=0)
    ax.annotate("hydrogen densities of the CO\nabsorbers, log n_H = 0.9-2.3\n(Klimenko et al. 2020)",
                (9, 2.5), fontsize=8, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.005, 10)
    ax.set_xlabel("density of the colliding H₂ molecules (per cm³)")
    ax.set_ylabel("CO excitation temperature − CMB temperature (K)")
    ax.set_title("Why molecular thermometers read slightly warm", loc="left")
    ax.legend(loc="lower right")
    ps.caption(fig, "Our statistical-equilibrium calculation for carbon monoxide's J = 1-0 transition "
                    f"in an optically thin cloud at {co['t_kin']:.0f} K, with radiative rates and collision "
                    "rates with H₂ from the Leiden LAMDA database. Without collisions the levels follow "
                    "the CMB exactly; collisions add a fraction of a kelvin to a few kelvin at the absorbers' "
                    "densities, which is why published excitation temperatures are corrected downwards "
                    "(Klimenko et al. 2020 also include collisions with H atoms and He).")
    fig.tight_layout()
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Test 3: the Tolman test
# ---------------------------------------------------------------------------

def fig_tolman_diagram(to, path, eta=2.0):
    local, distant = to["local"], to["distant"]
    nl = tolman.nonlinearity_correction()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    for ax, model, col, ls in ((axes[0], "lcdm", ps.C_EXPANSION, ps.LS_EXPANSION),
                               (axes[1], "tired", ps.C_TIRED, ps.LS_TIRED)):
        fits, _ = tolman.fit_local(local, model)
        f = fits[eta]
        l = local[local["eta"] == eta]
        lr = l["logR"] - l["f"] + tolman.log_a(l["z"], model)
        sb0 = l["SB"] - 2.5 * tolman.MODELS[model]["n_true"] * np.log10(1 + l["z"])
        ax.scatter(lr, sb0, s=10, color=ps.C_DATA_SOFT, lw=0, alpha=0.7,
                   label="nearby ellipticals (Postman & Lauer), corrected to z = 0")
        x = np.linspace(3.2, 5.4, 200)
        line, _ = tolman.local_line(f, x, "R", nl)
        ax.plot(x, line, color=ps.INK, lw=1.4, label="z = 0 relation (Sandage & Lubin 2001)")
        z = TOLMAN_CLUSTERS["Cl 1604+4321"]["z"]
        n_true = tolman.MODELS[model]["n_true"]
        ax.plot(x, line + 2.5 * n_true * np.log10(1 + z), color=col, ls=ls, lw=1.8,
                label=f"prediction at z = {z:.2f} without evolution, " + rf"$(1+z)^{n_true}$")
        d = distant[(distant["eta"] == eta)].copy()
        d["logR"] = d["logr_arcsec"] + tolman.log_a(d["z"], model)
        r = d[d["band"] == "R"]
        ax.plot(r["logR"], r["SB"], "o", color=ps.C_DATA, ms=6,
                label="Cl 1604+4321, z = 0.92 (HST, R band)")
        ax.set_xlabel("log₁₀ radius (parsecs) at η = 2.0, " +
                      ("flat ΛCDM distances" if model == "lcdm" else "static tired-light distances"))
        ax.legend(loc="lower right", fontsize=7.5)
    axes[0].set_ylabel("mean surface brightness ⟨SB⟩ (mag per arcsec², R band)")
    res = to["results"]
    axes[0].set_title(f"Expanding universe: n = {res['lcdm']['n_R']:.2f} ± {res['lcdm']['n_R_err']:.2f}",
                      loc="left")
    axes[1].set_title(f"Tired light: n = {res['tired']['n_R']:.2f} ± {res['tired']['n_R_err']:.2f}",
                      loc="left")
    axes[0].annotate("brighter than predicted without evolution:\n"
                     "their stars were younger then", (3.35, 25.3), fontsize=8.5,
                     color=ps.INK_SECONDARY)
    axes[1].annotate("fainter than predicted without evolution:\n"
                     "tired light needs galaxies that were\n"
                     "dimmer in the past, which aging stars cannot do", (3.25, 25.3), fontsize=8.5,
                     color=ps.INK_SECONDARY)
    ps.caption(fig, "Fainter is up (magnitudes). Surface brightness needs no distance; the physical "
                    "radius does, so each panel uses its own world model's distances. Galaxies below their "
                    "model's no-evolution line are brighter than it predicts, as expected if their stars were "
                    "younger then; galaxies above it would have to have been fainter in the past, which aging "
                    "stars cannot do. Data: Sandage & Lubin 2001 (Paper I), Lubin & Sandage 2001 (Paper III).")
    fig.tight_layout()
    return ps.save(fig, path)


def fig_tolman_exponents(to, path):
    res = to["results"]
    order = ["q0=0", "q0=0.5", "q0=1", "lcdm", "tired"]
    labels = ["empty\n(q₀ = 0)", "matter only\n(q₀ = ½)", "closed\n(q₀ = 1)", "flat ΛCDM", "tired light"]
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    x = np.arange(len(order))
    # Expected n after luminosity evolution 2.5 log (1+z)^p with the p ranges
    # of Bruzual & Charlot models quoted in Paper IV (R: 0.85-2.36, I: 0.76-2.07).
    # R band: p = 0.85-2.36; I band: p = 0.76-2.07 (Paper IV, Bruzual & Charlot models).
    for off, (plo, phi) in ((-0.09, (0.85, 2.36)), (0.09, (0.76, 2.07))):
        for xi in range(4):
            ax.fill_between([xi + off - 0.08, xi + off + 0.08], 4 - phi, 4 - plo,
                            color=ps.BLUE_150, alpha=0.8, lw=0,
                            label="expected when expanding: 4 minus stellar evolution"
                            if (xi == 0 and off < 0) else None)
    ax.plot([-0.4, 3.4], [4, 4], color=ps.C_EXPANSION, lw=1.8,
            label="expanding, no evolution (n = 4)")
    ax.fill_between([3.6, 4.4], -0.5, 1.0, color=ps.C_TIRED, alpha=0.15, lw=0,
                    label="allowed for tired light: 1 minus evolution")
    ax.plot([3.6, 4.4], [1, 1], color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8,
            label="tired light, no evolution (n = 1)")
    for k, (m, mark) in enumerate(zip(("R", "I"), ("o", "s"))):
        v = [res[o][f"n_{m}"] for o in order]
        e = [res[o][f"n_{m}_err"] for o in order]
        ax.errorbar(x + (k - 0.5) * 0.18, v, e, fmt=mark, color=ps.C_DATA, ms=6, lw=1.2,
                    mfc=ps.C_DATA if mark == "o" else ps.SURFACE, label=f"measured, {m} band")
    pub = to["published"]
    ax.errorbar([1 - 0.3, 1 + 0.3], [pub["n_R"][0], pub["n_I"][0]], [pub["n_R"][1], pub["n_I"][1]],
                fmt="x", color=ps.C_THIRD, ms=7, mew=1.6, lw=1.2,
                label="Lubin & Sandage 2001 (q₀ = ½; R left, I right)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-0.5, 4.6)
    ax.set_ylabel("exponent n of the dimming, 2.5 log₁₀(1 + z)ⁿ")
    ax.set_title("Surface-brightness dimming of galaxies at z = 0.76-0.92, in five world models",
                 loc="left")
    ax.legend(loc="lower left", fontsize=7.5)
    ps.caption(fig, "Each world model gives its own physical radii and hence its own measured n "
                    "(points: eta = 1.7 and 2.0 combined). Expanding models predict n = 4 before "
                    "evolution; tired light predicts n = 1. Galaxies were brighter in the past (younger "
                    "stars), which can only lower n, so tired light needs n <= 1.")
    fig.tight_layout()
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def fig_summary(td, tc, to, path):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    # 1. Time dilation (ratios to the z ~ 0.5 bin: no model assumed for the anchor).
    ax = axes[0]
    b = td["template"]["bins"]
    ref = b["ref"]
    x = b["one_plus_z"] / b["one_plus_z"][ref]
    xx = np.linspace(0.7, 1.4, 10)
    ax.plot(xx, xx, color=ps.C_EXPANSION, lw=1.8)
    ax.plot(xx, np.ones_like(xx), color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    ax.errorbar(x, b["D_ratio"], b["D_ratio_err"], fmt="o", color=ps.C_DATA, ms=5, lw=1.2)
    ax.set_xlabel("(1 + z) relative to supernovae at z ≈ 0.5")
    ax.set_ylabel("how much slower they evolve (relative)")
    ax.set_title("1. Supernovae run in slow motion", loc="left")
    t = td["template"]
    ax.annotate(f"b = {t['b']:.3f} ± {t['b_err_bootstrap']:.3f}\n(expansion: 1, tired light: 0)",
                (0.72, 1.22), fontsize=8, color=ps.INK_SECONDARY)
    _prediction_legend(ax, loc="lower right")
    # 2. CMB temperature.
    ax = axes[1]
    base = tcmb.select(tc["data"], "baseline")
    base = base[(base["z"] > 0)]
    xx = np.geomspace(1, 8, 100)
    ax.plot(xx, T0 * xx, color=ps.C_EXPANSION, lw=1.8)
    ax.plot(xx, T0 * np.ones_like(xx), color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)

    def yfun(d):
        return (d["T"].to_numpy(), np.hypot(d["err_minus"], d["err_sys"]).to_numpy(),
                np.hypot(d["err_plus"], d["err_sys"]).to_numpy())
    for label, sel, m in MARKERS:
        d = base[sel(base) & (base["kind"] != "upper")]
        y, lo, hi = yfun(d)
        ax.errorbar(1 + d["z"], y, [lo, hi], fmt=m, color=ps.C_DATA, ms=6 if m == "*" else 4,
                    lw=0.9, mfc=ps.SURFACE if m in "Dv" else ps.C_DATA)
    ax.plot([1], [T0], "*", color=ps.C_DATA, ms=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 3, 5, 7])
    ax.set_xticklabels(["1", "2", "3", "5", "7"])
    ax.set_yticks([3, 5, 10, 20, 30])
    ax.set_yticklabels(["3", "5", "10", "20", "30"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("1 + z")
    ax.set_ylabel("CMB temperature (K)")
    ax.set_title("2. The CMB was hotter in the past", loc="left")
    bt = tc["baseline"]
    ax.annotate(f"β = {bt['beta']:.3f} ± {max(bt['plus'], bt['minus']):.3f}\n"
                f"(expansion: 0, tired light: 1)", (1.03, 18), fontsize=8, color=ps.INK_SECONDARY)
    # 3. Tolman.
    ax = axes[2]
    res = to["results"]
    cats = [("flat ΛCDM", "lcdm"), ("tired light", "tired")]
    for k, (lab, m) in enumerate(cats):
        for j, (band, mark) in enumerate((("R", "o"), ("I", "s"))):
            ax.errorbar(k + (j - 0.5) * 0.2, res[m][f"n_{band}"], res[m][f"n_{band}_err"], fmt=mark,
                        color=ps.C_DATA, ms=6, lw=1.2, mfc=ps.C_DATA if mark == "o" else ps.SURFACE,
                        label=f"measured, {band} band" if k == 0 else None)
    ax.fill_between([-0.35, 0.35], 4 - 2.36, 4 - 0.76, color=ps.BLUE_150, alpha=0.7, lw=0)
    ax.plot([-0.35, 0.35], [4, 4], color=ps.C_EXPANSION, lw=1.8)
    ax.fill_between([0.65, 1.35], -0.3, 1.0, color=ps.C_TIRED, alpha=0.15, lw=0)
    ax.plot([0.65, 1.35], [1, 1], color=ps.C_TIRED, ls=ps.LS_TIRED, lw=1.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([c[0] for c in cats])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.3, 4.5)
    ax.set_ylabel("surface-brightness dimming exponent n")
    ax.set_title("3. Distant galaxies are too dim", loc="left")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.annotate("shaded: allowed once galaxies'\nstars are allowed to age", (0.55, 2.7),
                fontsize=7.5, color=ps.INK_SECONDARY)
    fig.suptitle("Is the universe really expanding? Three tests that tired light fails",
                 x=0.04, ha="left")
    ps.caption(fig, "Left: DES 5-year Type Ia supernovae, time stretch per redshift bin (our global "
                    "template fit). Middle: published CMB temperatures (SZ clusters, quasar absorbers, "
                    "water at z = 6.34). Right: Lubin & Sandage's HST galaxies at z = 0.76-0.92 against "
                    "nearby ellipticals, radii computed in each world model.")
    fig.tight_layout()
    return ps.save(fig, path)


def make_all(td, tc, to):
    written = []
    jobs = [
        ("00_summary", lambda p: fig_summary(td, tc, to, p)),
        ("01_slow_motion", lambda p: fig_slow_motion(td, p)),
        ("02_light_curve_widths", lambda p: fig_widths(td, p)),
        ("03_measuring_b", lambda p: fig_b_measurements(td, p)),
        ("04_dilation_by_redshift", lambda p: fig_dilation_bins(td, p)),
        ("05_cmb_temperature", lambda p: fig_tcmb(tc, p)),
        ("06_sz_null", lambda p: fig_sz_null(tc, p)),
        ("07_molecular_thermometer", lambda p: fig_co_thermometer(tc, p)),
        ("08_tolman_diagram", lambda p: fig_tolman_diagram(to, p)),
        ("09_tolman_exponents", lambda p: fig_tolman_exponents(to, p)),
    ]
    for name, job in jobs:
        written += job(FIGURE_DIR / name)
        print(f"  {name}")
    return written
