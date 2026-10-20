"""Every figure of the project, drawn from the cached analysis results."""

from __future__ import annotations

import numpy as np
import healpy as hp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import catalog as K  # noqa: E402
import config as C  # noqa: E402
import cosmology as CO  # noqa: E402
import footprint as F  # noqa: E402
import literature as L  # noqa: E402
import model as MD  # noqa: E402
import plotstyle as PS  # noqa: E402
import selection as S  # noqa: E402

PS.use_style()
OUT = C.FIGURE_DIR


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _levels(p, fracs=(0.68, 0.95)):
    """Density thresholds enclosing the given posterior fractions."""
    flat = np.sort(p.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    return [flat[min(np.searchsorted(cum, f), flat.size - 1)] for f in sorted(fracs, reverse=True)]


def _grid_post(grid):
    return np.exp(grid["lnl"] - grid["lnl"].max())


def _chain_density(x, y, xr, yr, bins=70, smooth=1.4):
    from scipy.ndimage import gaussian_filter

    H, xe, ye = np.histogram2d(x, y, bins=bins, range=[xr, yr])
    H = gaussian_filter(H, smooth)
    return 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H


def _gauss_ellipse(ax, x, y, cov, color, label=None, lw=1.4, ls="-"):
    vals, vecs = np.linalg.eigh(cov)
    t = np.linspace(0, 2 * np.pi, 200)
    for k, chi in enumerate((2.30, 6.18)):  # 68% and 95% for 2 parameters
        r = np.sqrt(chi * vals)
        pts = vecs @ np.vstack([r[0] * np.cos(t), r[1] * np.sin(t)])
        ax.plot(x + pts[0], y + pts[1], color=color, lw=lw if k == 0 else lw * 0.7, ls=ls,
                alpha=1 if k == 0 else 0.7, label=label if k == 0 else None)


def _bin_centres(res):
    cr_edges = res["counts"]["cr_edges"]
    z_edges = res["counts"]["z_edges"]
    lcr = np.exp(cr_edges)
    return np.sqrt(lcr[1:] * lcr[:-1]), 0.5 * (z_edges[1:] + z_edges[:-1]), lcr, z_edges


def _galactic_xy(ra, dec):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    g = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg).galactic
    lon = np.where(g.l.deg > 180, g.l.deg - 360, g.l.deg)
    return -np.radians(lon), np.radians(g.b.deg)  # longitude increases to the left


# ---------------------------------------------------------------------------
# 01  sky map
# ---------------------------------------------------------------------------


def draw_sky(ax, fp, df, primary, psz2, small=False):
    # exposure map in Galactic coordinates, sampled on a lon/lat grid
    nl, nb = (720, 360) if not small else (480, 240)
    lon = np.linspace(-np.pi, np.pi, nl)
    lat = np.linspace(-np.pi / 2, np.pi / 2, nb)
    LON, LAT = np.meshgrid(lon, lat)
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    eq = SkyCoord(l=(-np.degrees(LON)) % 360 * u.deg, b=np.degrees(LAT) * u.deg, frame="galactic").icrs
    pix = hp.ang2pix(F.NSIDE, eq.ra.deg, eq.dec.deg, lonlat=True, nest=True)
    expo = np.array(fp["exposure"], float)[pix]
    expo[~fp["mask_de"][pix]] = np.nan
    cosmo = fp["mask_cosmo"][pix]
    im = ax.pcolormesh(LON, LAT, np.log10(expo), cmap=PS.SEQ_GREY, vmin=1.8, vmax=3.3, shading="auto",
                       rasterized=True, zorder=0)
    cs = ax.contour(LON, LAT, cosmo.astype(float), levels=[0.5], colors=[PS.INK_SECONDARY], linewidths=0.5, zorder=1)
    cs.set_rasterized(True)  # ~300k vertices: far too heavy as SVG paths
    in_cos = set(df.detuid)
    other = ~primary.detuid.isin(in_cos)
    x, y = _galactic_xy(primary.ra[other], primary.dec[other])
    ax.scatter(x, y, s=0.6 if small else 1.0, color=PS.CONTEXT, lw=0, zorder=2, rasterized=True)
    x, y = _galactic_xy(df.ra, df.dec)
    sc = ax.scatter(x, y, s=1.2 if small else 2.2, c=df.z, cmap=PS.SEQ_BLUE, vmin=0.1, vmax=0.8, lw=0, zorder=3,
                    rasterized=True)
    if psz2 is not None:
        x, y = _galactic_xy(psz2.ra, psz2.dec)
        ax.scatter(x, y, s=4 if small else 7, facecolors="none", edgecolors=PS.ORANGE, lw=0.5, zorder=4,
                   rasterized=True)
    ax.grid(True, color=PS.GRID, lw=0.4)
    ax.set_xticklabels([])
    ax.tick_params(labelsize=6.5, colors=PS.MUTED)
    return im, sc


def _gx(l):
    """Plot x (radians) of Galactic longitude l in degrees (l increases to the left)."""
    return -np.radians(l - 360 if l > 180 else l)


def fig_sky(res):
    fp = F.build()
    df = K.load_cosmology_sample(fp)
    primary = K.load_primary()
    psz2 = K.load_psz2()
    fig = plt.figure(figsize=(10.5, 6.6))
    ax = fig.add_axes([0.02, 0.27, 0.96, 0.66], projection="mollweide")
    im, sc = draw_sky(ax, fp, df, primary, psz2)
    ax.set_title("eRASS1: the galaxy clusters of the western Galactic hemisphere (Galactic coordinates)", pad=12)
    for txt, (l, b), (lt, bt) in (
        ("Large Magellanic Cloud (masked),\nnext to the south ecliptic pole,\nwhere the exposure is deepest", (280.5, -32.9), (238, -6)),
        ("Galactic plane: no optical\ncoverage, strong absorption", (330.0, 0.0), (318, 9)),
        ("Virgo cluster\n(masked)", (283.8, 74.5), (225, 52)),
    ):
        PS.annotate(ax, txt, (_gx(l), np.radians(b)), (_gx(lt), np.radians(bt)), fontsize=7, ha="center")
    ax.text(0.86, 0.985, "eROSITA-DE half:\nreleased", transform=ax.transAxes, ha="center", va="top", fontsize=8,
            color=PS.INK_SECONDARY)
    ax.text(0.14, 0.985, "Russian half: not released\n(only Planck clusters shown)", transform=ax.transAxes, ha="center",
            va="top", fontsize=8, color=PS.MUTED)
    cax1 = fig.add_axes([0.08, 0.13, 0.36, 0.022])
    cb = fig.colorbar(sc, cax=cax1, orientation="horizontal")
    cb.set_label("photometric redshift (cosmology sample, 5,259 clusters)", color=PS.INK_SECONDARY, fontsize=8)
    cb.outline.set_visible(False)
    cax2 = fig.add_axes([0.56, 0.13, 0.36, 0.022])
    cb2 = fig.colorbar(im, cax=cax2, orientation="horizontal")
    cb2.set_label("eRASS1 exposure (s)", color=PS.INK_SECONDARY, fontsize=8)
    cb2.set_ticks([2, 2.5, 3])
    cb2.set_ticklabels(["100", "320", "1000"])
    cb2.outline.set_visible(False)
    handles = [Line2D([], [], ls="", marker="o", ms=4, color=PS.BLUE_550, label="cosmology sample"),
               Line2D([], [], ls="", marker="o", ms=4, color=PS.CONTEXT, label="other eRASS1 clusters (12,247 in all)"),
               Line2D([], [], ls="", marker="o", ms=5, mfc="none", mec=PS.ORANGE, label="Planck PSZ2 (whole sky)"),
               Line2D([], [], color=PS.INK_SECONDARY, lw=0.8, label="cosmology footprint (our reconstruction)")]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.215), ncol=4, fontsize=7.5)
    fig.text(0.5, 0.03, "Background: eRASS1 0.2-2.3 keV exposure (HiPS map, eROSITA-DE DR1). Footprint: eROSITA-DE x Legacy Surveys "
             "DR10-South x Bulbul et al.'s masks, read from their Fig. 2.", ha="center", fontsize=7.5, color=PS.MUTED)
    return PS.save(fig, OUT / "01_sky")


# ---------------------------------------------------------------------------
# 02  survey depth and selection function
# ---------------------------------------------------------------------------


def fig_survey(res):
    sv = res["survey"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    ax = axes[0]
    e, cum = sv["cum_total"]
    ax.plot(e, cum, color=PS.CONTEXT, lw=1.6, label="all eROSITA-DE sky (ours)")
    e, cum = sv["cum_clean"]
    ax.plot(e, cum, color=PS.BLUE, lw=1.8, label="cleaned footprint (ours)")
    for curve, col, lab in ((sv["fig2b_total"], PS.INK_SECONDARY, "Bulbul+ 2024 Fig. 2 (read off)"),
                            (sv["fig2b_clean"], PS.INK, None)):
        k = curve[::12]
        ax.plot(k[:, 0], k[:, 1], ls="", marker="o", ms=3.5, mfc="none", mec=col, mew=0.8, label=lab)
    ax.set(xscale="log", yscale="log", xlim=(55, 5000), ylim=(0.3, 3e4), xlabel="exposure T (s)",
           ylabel="sky area deeper than T (deg$^2$)")
    ax.set_title("How deep is the survey?")
    ax.legend(loc="lower left", fontsize=7.5)
    PS.annotate(ax, "south ecliptic pole", (2500, 35), (900, 600), fontsize=7.5)

    ax = axes[1]
    model = S.SelectionModel()
    cr = np.geomspace(0.01, 20, 300)
    cols = {50.0: PS.BLUE_300, 150.0: PS.BLUE, 500.0: PS.BLUE_750}
    for (z, T), (lx, P) in sorted(res["fig4_curves"].items()):
        if z != 0.3:
            continue
        ax.plot(10**lx, P, color=PS.CONTEXT, lw=4.5, solid_capstyle="round", zorder=1)
        ax.plot(cr, model(cr, z, T), color=cols[T], lw=1.5, zorder=2)
        ax.text(model.cr50(z, T)[0] * 1.15, 0.47, f"{T:.0f} s", color=PS.INK_SECONDARY, fontsize=7.5)
    for (z, T), pts in res["clerc_points"].items():
        if T == 2000.0 and z == 0.2:
            P = sorted(pts)
            ax.plot([pts[p] for p in P], P, ls="", marker="D", ms=4.5, mfc="none", mec=PS.ORANGE, mew=1)
            ax.plot(cr, model(cr, 0.2, 2000.0, nh=3e20), color=PS.ORANGE, lw=1.2)
            ax.text(0.03, 0.82, "2000 s, z = 0.2", color=PS.INK_SECONDARY, fontsize=7.5)
    ax.set(xscale="log", xlabel="unabsorbed count rate (counts/s, 0.2-2.3 keV)", ylabel="probability of detection",
           ylim=(-0.02, 1.05))
    ax.set_title("Which clusters get detected? (z = 0.3)")
    handles = [Line2D([], [], color=PS.CONTEXT, lw=4.5, label="published model (Ghirardini+ Fig. 4)"),
               Line2D([], [], color=PS.BLUE, lw=1.5, label="our fit"),
               Line2D([], [], ls="", marker="D", ms=4.5, mfc="none", mec=PS.ORANGE, label="Clerc+ 2024 (raster)")]
    ax.legend(handles=handles, loc="upper left", fontsize=7)

    ax = axes[2]
    zz = np.linspace(0.1, 0.8, 71)
    for T, col in ((80.0, PS.BLUE_300), (170.0, PS.BLUE), (800.0, PS.BLUE_750)):
        c50 = np.array([model.cr50(z, T)[0] for z in zz])
        ax.plot(zz, c50, color=col, lw=1.6)
        ax.text(0.81, c50[-1], f"{T:.0f} s", va="center", fontsize=7.5, color=PS.INK_SECONDARY)
    ax.axvspan(0.5, 0.8, color=PS.WASH, zorder=0)
    ax.text(0.65, 1.4, "extrapolated\nbeyond the\npublished curves", ha="center", fontsize=7, color=PS.MUTED)
    ax.set(yscale="log", xlim=(0.1, 0.86), ylim=(0.02, 3), xlabel="redshift",
           ylabel="count rate detected half of the time")
    ax.set_title("Detection threshold vs redshift")
    PS.annotate(ax, "5th, 50th, 95th\npercentile exposure", (0.3, model.cr50(0.3, 170.0)[0]), (0.16, 0.05),
                fontsize=7)
    PS.caption(fig, "Left: cumulative area vs exposure of our reconstructed footprints against the curves of Bulbul et al. 2024. "
               "Middle/right: our fit to the published selection-function curves (extracted from the vector PDF) and its "
               "check at 2000 s. All for the reference absorption N_H = 3.5e20 cm^-2.", y=-0.04)
    fig.tight_layout()
    return PS.save(fig, OUT / "02_survey_selection")


# ---------------------------------------------------------------------------
# 03  data vs model
# ---------------------------------------------------------------------------


def draw_counts_z(ax, res, show_pub=True, curves=None):
    crc, zc, _, z_edges = _bin_centres(res)
    d = res["counts"]["data"].sum(0)
    if curves:
        for lab, N, col, lw in curves:
            ax.stairs(N.sum(0), z_edges, color=col, lw=lw, label=lab, baseline=None)
    else:
        if "pred_mcmc" in res["counts"]:
            ax.stairs(res["counts"]["pred_mcmc"].sum(0), z_edges, color=PS.BLUE, lw=1.9, label="best fit, relation free",
                      baseline=None)
        ax.stairs(res["counts"]["pred_best"].sum(0), z_edges, color=PS.BLUE_300, lw=1.3,
                  label="best fit, published relation fixed", baseline=None)
        if show_pub:
            ax.stairs(res["counts"]["pred_published"].sum(0), z_edges, color=PS.ORANGE, lw=1.2, ls="--",
                      label="at Ghirardini+ 2024 best fit", baseline=None)
    ax.errorbar(zc, d, yerr=np.sqrt(d), ls="", marker="o", ms=4, color=PS.INK, mec=PS.SURFACE, mew=0.6,
                elinewidth=0.9, label="eRASS1 (weighted by purity)", zorder=5)
    ax.set(xlabel="photometric redshift", ylabel="clusters per 0.05 in z", xlim=(0.1, 0.8))


def draw_counts_cr(ax, res, show_pub=True, curves=None, cumulative=False):
    crc, zc, lcr, _ = _bin_centres(res)
    d = res["counts"]["data"].sum(1)
    edges = lcr.copy()
    if cumulative:
        dd = np.cumsum(d[::-1])[::-1]
        for lab, N, col, lw in curves:
            ax.plot(edges[:-1], np.cumsum(N.sum(1)[::-1])[::-1], color=col, lw=lw, label=lab, drawstyle="steps-post")
        ax.errorbar(edges[:-1], dd, yerr=np.sqrt(dd), ls="", marker="o", ms=4, color=PS.INK, mec=PS.SURFACE,
                    mew=0.6, elinewidth=0.9, label="eRASS1", zorder=5)
        ax.set(xscale="log", yscale="log", xlim=(0.02, 20), ylim=(0.8, 8000),
               xlabel="count rate (counts/s)", ylabel="clusters brighter than this")
        return
    if curves:
        for lab, N, col, lw in curves:
            ax.stairs(N.sum(1), edges, color=col, lw=lw, label=lab, baseline=None)
    else:
        if "pred_mcmc" in res["counts"]:
            ax.stairs(res["counts"]["pred_mcmc"].sum(1), edges, color=PS.BLUE, lw=1.9, label="relation free",
                      baseline=None)
        ax.stairs(res["counts"]["pred_best"].sum(1), edges, color=PS.BLUE_300, lw=1.3, label="relation fixed",
                  baseline=None)
        if show_pub:
            ax.stairs(res["counts"]["pred_published"].sum(1), edges, color=PS.ORANGE, lw=1.2, ls="--",
                      label="at Ghirardini+ best fit", baseline=None)
    ax.errorbar(crc, d, yerr=np.sqrt(np.maximum(d, 1e-9)), ls="", marker="o", ms=4, color=PS.INK, mec=PS.SURFACE,
                mew=0.6, elinewidth=0.9, label="eRASS1", zorder=5)
    ax.set(xscale="log", xlim=(0.004, 40), xlabel="measured count rate (counts/s, 0.2 dex bins)",
           ylabel="clusters per bin")


def fig_counts(res):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.1), gridspec_kw={"width_ratios": [1, 1, 1.1]})
    draw_counts_z(axes[0], res)
    axes[0].set_title("Clusters by redshift")
    axes[0].legend(fontsize=7.5, loc="upper right")
    draw_counts_cr(axes[1], res)
    axes[1].set_title("Clusters by X-ray brightness")
    ax = axes[2]
    key = "pred_mcmc" if "pred_mcmc" in res["counts"] else "pred_best"
    d, m = res["counts"]["data"], res["counts"][key]
    pr = np.where(m > 1, (d - m) / np.sqrt(np.maximum(m, 1e-9)), np.nan)
    _, _, lcr, z_edges = _bin_centres(res)
    edges = lcr.copy()
    edges[0], edges[-1] = 0.0025, 60
    im = ax.pcolormesh(z_edges, edges, pr, cmap=PS.DIVERGING, vmin=-4, vmax=4, shading="flat")
    ax.set(yscale="log", ylim=(0.01, 30), xlabel="photometric redshift", ylabel="count rate (counts/s)")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("(observed - predicted) / sqrt(predicted)", color=PS.INK_SECONDARY, fontsize=7.5)
    cb.outline.set_visible(False)
    g = res["counts"]["gof_mcmc" if key == "pred_mcmc" else "gof"]
    ax.set_title("Where the fit misses" + (" (relation free)" if key == "pred_mcmc" else ""))
    ax.text(0.02, 0.97, f"chi^2 = {g['chi2_2d']:.0f} for {g['nbins_2d']} bins", transform=ax.transAxes,
            fontsize=7.5, va="top", color=PS.INK_SECONDARY)
    PS.caption(fig, "Counts weighted by each cluster's probability of being real (1 - PCONT); Poisson errors. Models: our forward "
               "model at the MCMC medians (relation free), at the best fit with the published relation fixed, and at "
               "Ghirardini et al.'s published best fit. Residual map: relation-free model, bins with > 1 predicted cluster.",
               y=-0.04)
    fig.tight_layout()
    return PS.save(fig, OUT / "03_counts_vs_model")


# ---------------------------------------------------------------------------
# 04  how clumpy?
# ---------------------------------------------------------------------------


def fig_clumpy(res):
    cols = {0.76: PS.BLUE_150, 0.81: PS.BLUE_300, 0.86: PS.BLUE_550}
    curves = [(f"sigma_8 = {s8:.2f}", N, cols[s8], 1.8) for s8, N in sorted(res["clumpy"].items())]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    draw_counts_z(axes[0], res, curves=curves)
    axes[0].set_title("The same survey in three universes (Omega_m = 0.3)")
    axes[0].legend(fontsize=7.5, loc="upper right")
    draw_counts_cr(axes[1], res, curves=curves, cumulative=True)
    axes[1].set_title("Clusters brighter than a given count rate")
    for s8, N in sorted(res["clumpy"].items()):
        n = N.sum()
        axes[1].text(0.03, 0.06 + 0.07 * list(cols).index(s8), f"sigma_8 = {s8:.2f}: {n:,.0f} clusters predicted",
                     transform=axes[1].transAxes, fontsize=7.5, color=PS.INK_SECONDARY)
    axes[1].text(0.03, 0.27, f"observed: {res['counts']['data'].sum():,.0f}", transform=axes[1].transAxes,
                 fontsize=7.5, color=PS.INK)
    axes[1].legend(fontsize=7.5, loc="upper right")
    PS.caption(fig, "Predictions with the published count rate-mass relation. A 13% change in sigma_8 (0.76 -> 0.86) "
               "changes the number of detectable clusters by more than half: the most massive halos are exponentially rare.",
               y=-0.04)
    fig.tight_layout()
    return PS.save(fig, OUT / "04_how_clumpy")


# ---------------------------------------------------------------------------
# 05  Omega_m - sigma_8
# ---------------------------------------------------------------------------


def draw_om_s8(ax, res, compact=False):
    om_r, s8_r = (0.18, 0.42), (0.70, 1.06)
    # late-universe lensing bands (S8 = const)
    om = np.linspace(*om_r, 200)
    for lab, S8, lo, hi, kind in L.S8_MEASUREMENTS:
        if lab.startswith(("KiDS", "DES Y6 3x2pt")):
            ax.fill_between(om, (S8 - lo) * np.sqrt(0.3 / om), (S8 + hi) * np.sqrt(0.3 / om), color=PS.WASH, lw=0,
                            zorder=0)
    # CMB: D3 samples (ellipse from their covariance) and Planck 2018
    d3 = res["d3"]["samples"]
    cov = np.cov(np.vstack([d3["Om"], d3["sigma8"]]))
    _gauss_ellipse(ax, np.mean(d3["Om"]), np.mean(d3["sigma8"]), cov, PS.AQUA, "CMB: project D3 -> CAMB")
    p = L.PLANCK2018
    ax.errorbar(p["Om"][0], p["sigma8"][0], xerr=p["Om"][1], yerr=p["sigma8"][1], ls="", marker="D", ms=5,
                color=PS.AQUA, mec=PS.SURFACE, mew=0.6, elinewidth=1, label="Planck 2018 CMB", zorder=6)
    # ours: fixed-relation grid (line) and MCMC (filled)
    g = res["grid"]
    P = _grid_post(g)
    ax.contour(g["om"], g["s8"], P.T, levels=_levels(P), colors=[PS.BLUE_550], linewidths=[0.8, 1.2], zorder=4)
    if "mcmc" in res:
        ch = res["mcmc"]["chain"]
        x, y, H = _chain_density(ch[:, 0], ch[:, 1], om_r, s8_r)
        lev = _levels(H)
        ax.contourf(x, y, H.T, levels=[lev[0], lev[1], H.max() * 1.01], colors=[PS.BLUE_150, PS.BLUE_300],
                    alpha=0.85, zorder=3)
    gh = L.GHIRARDINI2024
    ax.errorbar(gh["Om"][0], gh["sigma8"][0], xerr=[[gh["Om"][1]], [gh["Om"][2]]], yerr=gh["sigma8"][1], ls="",
                marker="s", ms=5, color=PS.ORANGE, mec=PS.SURFACE, mew=0.6, elinewidth=1, zorder=7,
                label="eRASS1 published (Ghirardini+ 2024)")
    ax.set(xlim=om_r, ylim=s8_r, xlabel="Omega_m (matter density)", ylabel="sigma_8 (clumpiness today)")
    handles = [Line2D([], [], color=PS.BLUE_300, lw=6, label="ours: relation free (MCMC, 68/95%)"),
               Line2D([], [], color=PS.BLUE_550, lw=1.2, label="ours: published relation fixed"),
               Line2D([], [], ls="", marker="s", ms=5, color=PS.ORANGE, label="eRASS1 published (Ghirardini+ 2024)"),
               Line2D([], [], color=PS.AQUA, lw=1.4, label="CMB: project D3, grown to today"),
               Line2D([], [], ls="", marker="D", ms=5, color=PS.AQUA, label="Planck 2018 CMB"),
               Line2D([], [], color=PS.WASH, lw=6, label="cosmic shear S_8 (KiDS-Legacy, DES Y6)")]
    ax.legend(handles=handles, loc="upper right", fontsize=7 if not compact else 6.3)


def fig_om_s8(res):
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    draw_om_s8(ax, res)
    ax.set_title("How much matter, how clumpy?")
    PS.caption(fig, "Counts constrain the combination S_8 = sigma_8 (Omega_m/0.3)^0.5 best: the banana runs along it. "
               "Shear bands: S_8 +-1 sigma. D3 ellipse: 68/95% from 300 CAMB runs drawn from D3's fit.", y=-0.02)
    return PS.save(fig, OUT / "05_omegam_sigma8")


# ---------------------------------------------------------------------------
# 06  S8 comparison
# ---------------------------------------------------------------------------


def draw_s8(ax, res, compact=False):
    rows = []
    if "mcmc" in res:
        m = res["mcmc"]["summary"]["S8"]
        rows.append(("This project: eRASS1 counts, relation free", m[0], m[0] - m[1], m[2] - m[0], "ours"))
    f = res["grid"]["S8"]
    rows.append(("This project: published relation fixed", f[0], f[0] - f[1], f[2] - f[0], "ours"))
    gh = L.GHIRARDINI2024["S8"]
    rows.append(("eRASS1 clusters (Ghirardini+ 2024)", gh[0], gh[1], gh[2], "published"))
    for lab, v, lo, hi, kind in L.S8_MEASUREMENTS:
        if kind == "clusters":
            rows.append((lab, v, lo, hi, "clusters"))
    d3 = res["d3"]["summary"]["S8"]
    rows.append(("CMB: project D3, grown to today", d3[0], d3[1], d3[1], "cmb"))
    p = L.PLANCK2018["S8"]
    rows.append(("Planck 2018 CMB", p[0], p[1], p[2], "cmb"))
    for lab, v, lo, hi, kind in L.S8_MEASUREMENTS:
        if kind == "cmb-lensing":
            rows.append((lab, v, lo, hi, "cmb"))
    for lab, v, lo, hi, kind in L.S8_MEASUREMENTS:
        if kind == "lensing":
            rows.append((lab, v, lo, hi, "lensing"))
    style = {"ours": (PS.BLUE, "o"), "published": (PS.ORANGE, "s"), "clusters": (PS.CONTEXT, "s"),
             "cmb": (PS.AQUA, "D"), "lensing": (PS.INK_SECONDARY, "^")}
    short = {"This project: eRASS1 counts, relation free": "ours: relation free",
             "This project: published relation fixed": "ours: relation fixed",
             "eRASS1 clusters (Ghirardini+ 2024)": "eRASS1 (Ghirardini+24)",
             "SPT clusters + DES/HST lensing (Bocquet+ 2024)": "SPT clusters (Bocquet+24)",
             "CMB: project D3, grown to today": "CMB: project D3",
             "Planck 2018 CMB": "Planck 2018 CMB",
             "ACT DR6 CMB lensing + BAO (Madhavacheril+ 2024)": "ACT DR6 lensing+BAO",
             "KiDS-Legacy cosmic shear (Wright+ 2025)": "KiDS-Legacy shear",
             "DES Y6 cosmic shear (2026)": "DES Y6 shear",
             "DES Y6 3x2pt (2026)": "DES Y6 3x2pt"}
    ax.axvspan(p[0] - p[1], p[0] + p[2], color=PS.WASH, zorder=0)
    for i, (lab, v, lo, hi, kind) in enumerate(rows):
        col, mk = style[kind]
        y = len(rows) - 1 - i
        ax.errorbar(v, y, xerr=[[lo], [hi]], ls="", marker=mk, ms=6, color=col, mec=PS.SURFACE, mew=0.6,
                    elinewidth=1.4, capsize=0)
        if not compact:
            ax.text(0.952, y, f"{v:.3f}", va="center", ha="right", fontsize=7, color=PS.INK_SECONDARY)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([(short[r[0]] if compact else r[0]) for r in rows[::-1]], fontsize=6.6 if compact else 7.5,
                       color=PS.INK)
    ax.tick_params(axis="y", length=0)
    ax.set(xlim=(0.74, 0.96) if not compact else (0.74, 0.95), ylim=(-0.7, len(rows) - 0.3),
           xlabel="S_8 = sigma_8 (Omega_m / 0.3)^0.5" if not compact else "S_8")
    ax.grid(axis="y", visible=False)
    ax.text(p[0], len(rows) - 0.55, "Planck", ha="center", fontsize=6.8, color=PS.MUTED)


def fig_s8(res):
    fig, ax = plt.subplots(figsize=(8.2, 4.9))
    draw_s8(ax, res)
    ax.set_title("Is the late universe less clumpy than the CMB predicts?  (S_8 tension)")
    PS.caption(fig, "Clusters (blue/orange), CMB and CMB lensing (aqua), galaxy shapes (grey triangles). Shaded: Planck 2018 +-1 sigma. "
               "Asymmetric errors as published.", y=-0.02)
    return PS.save(fig, OUT / "06_s8_comparison")


# ---------------------------------------------------------------------------
# 07  mass scale: the hydrostatic bias
# ---------------------------------------------------------------------------


def fig_mass_scale(res):
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), gridspec_kw={"width_ratios": [1.05, 0.95, 0.85]})
    ax = axes[0]
    ms = sorted(res["mass_scale"])
    s8 = np.array([res["mass_scale"][k]["S8"] for k in ms])
    ax.fill_between(ms, s8[:, 1], s8[:, 2], color=PS.BLUE_150, lw=0)
    ax.plot(ms, s8[:, 0], color=PS.BLUE, lw=1.8, marker="o", ms=4, label="eRASS1 counts (published relation)")
    p = L.PLANCK2018["S8"]
    ax.axhspan(p[0] - p[1], p[0] + p[2], color=PS.WASH, zorder=0)
    ax.text(1.14, p[0] + 0.004, "Planck CMB", ha="right", fontsize=7.5, color=PS.MUTED)
    for lab, v, err in L.HYDROSTATIC_BIAS:
        short = lab.split(" (")[0]
        ax.axvline(v, color=PS.CONTEXT, lw=0.8, zorder=1)
        ax.text(v + 0.006, 0.945, short.replace("Weighing the Giants lensing", "WtG lensing")
                .replace("CMB + Planck SZ counts need", "Planck SZ counts need")
                .replace("Planck's default assumption", "Planck default"),
                rotation=90, fontsize=6.6, color=PS.MUTED, va="top")
    ax.set(xlabel="assumed M_true / M_lensing  (1 = lensing; < 1 = hydrostatic-like)",
           ylabel="S_8 from the same counts", ylim=(0.66, 0.95), xlim=(0.55, 1.17))
    ax.set_title("A 10% mass error moves S_8 a lot")
    ax.legend(loc="lower right", fontsize=7.5)

    ax = axes[1]
    mt = res["hydrostatic"]["matches"]
    sel = mt.sz_snr >= 6
    ax.scatter(mt.x_m500[~sel], mt.sz_msz[~sel], s=9, color=PS.CONTEXT, lw=0, label="Planck S/N < 6")
    ax.scatter(mt.x_m500[sel], mt.sz_msz[sel], s=12, color=PS.ORANGE, edgecolor=PS.SURFACE, lw=0.4,
               label="Planck S/N >= 6")
    x = np.geomspace(1e14, 3e15, 10)
    ax.plot(x, x, color=PS.INK_SECONDARY, lw=0.9)
    ax.text(2.2e15, 2.2e15 * 1.07, "equal", ha="right", fontsize=7, color=PS.INK_SECONDARY)
    t = next(r for r in res["hydrostatic"]["table"] if r["estimator"] == "catalogue M500 column" and r["snr_min"] == 6)
    ax.plot(x, t["one_minus_b"] * x, color=PS.ORANGE, lw=1.2)
    ax.text(2.4e15, 2.4e15 * t["one_minus_b"] * 0.86, f"x {t['one_minus_b']:.2f}", ha="right", fontsize=7.5,
            color=PS.INK_SECONDARY)
    ax.set(xscale="log", yscale="log", xlim=(1.2e14, 2.5e15), ylim=(1.2e14, 2.5e15),
           xlabel="eRASS1 mass, lensing-calibrated (catalogue M500, Msun)", ylabel="Planck SZ mass (hydrostatic-based, Msun)")
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(FixedLocator([2e14, 5e14, 1e15, 2e15]))
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 1e14:.0f}e14" if v < 1e15 else f"{v / 1e15:.0f}e15"))
        axis.set_minor_formatter(NullFormatter())
    ax.set_title("Same clusters, two mass scales")
    ax.legend(loc="upper left", fontsize=7.5)

    ax = axes[2]
    rows = [(r["estimator"], r["one_minus_b"], r["err"]) for r in res["hydrostatic"]["table"] if r["snr_min"] == 6.0]
    rows += [(lab + "  [literature]", v, e) for lab, v, e in L.HYDROSTATIC_BIAS if e > 0]
    for i, (lab, v, e) in enumerate(rows):
        y = len(rows) - 1 - i
        lit = "literature" in lab
        ax.errorbar(v, y, xerr=e, ls="", marker="D" if lit else "o", ms=5.5, color=PS.AQUA if lit else PS.BLUE,
                    mec=PS.SURFACE, mew=0.6, elinewidth=1.2)
    names = {"catalogue M500 column": "eRASS1 catalogue M500",
             "catalogue PDF medians": "eRASS1 catalogue PDF medians",
             "our X-ray-only posteriors": "our X-ray-only posteriors",
             "forward: count rate predicted from Planck mass": "count rate predicted from Planck mass"}
    labels = [names.get(r[0], r[0].replace("  [literature]", "")) for r in rows]
    ax.set(xlim=(0.25, 1.4), ylim=(-0.7, len(rows) - 0.3), xlabel="(1 - b) = Planck mass / reference mass")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels[::-1], fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.axhline(len(rows) - 4 - 0.5, color=PS.GRID, lw=0.8)
    ax.text(1.38, len(rows) - 1, "this project\n(S/N >= 6)", ha="right", va="center", fontsize=6.8, color=PS.MUTED)
    ax.text(1.38, 1, "published", ha="right", va="center", fontsize=6.8, color=PS.MUTED)
    ax.grid(axis="y", visible=False)
    ax.set_title("Per-cluster comparison depends on the estimator")
    PS.caption(fig, "Left: the counts refitted with every mass rescaled; vertical lines: published hydrostatic biases. "
               "Middle/right: 374 eRASS1 clusters matched to Planck PSZ2 by Bulbul et al.; ratios for S/N >= 6. With a factor-2.7 "
               "scatter between count rate and mass, a single cluster's mass depends on the statistical treatment (see README).",
               y=-0.04)
    fig.tight_layout()
    return PS.save(fig, OUT / "07_mass_scale_hydrostatic")


# ---------------------------------------------------------------------------
# 08  masses and redshifts
# ---------------------------------------------------------------------------


def detection_mass(z, T, cosmo, sr, sel, nh=C.NH_REF, frac=0.5):
    lnm = np.log(np.geomspace(1e13, 5e15, 300))
    lncr = np.log(np.geomspace(1e-4, 300, 400))
    mu = sr.mean_lncr(cosmo, lnm, np.array([z]))[0]
    pdet = sel(np.exp(lncr), z, T, nh)
    g = np.exp(-0.5 * ((lncr[None, :] - mu[:, None]) / sr.sigma_X) ** 2)
    g /= g.sum(1, keepdims=True)
    p = g @ pdet
    return float(np.exp(np.interp(frac, p, lnm)))


def fig_masses(res):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    z, m = res["masses"]["z"], res["masses"]["cat_m500"]
    ax.scatter(z, m, s=2.5, color=PS.CONTEXT, lw=0, rasterized=True, label="cosmology sample (catalogue M500)")
    cosmo, sr, sel = CO.fiducial(), MD.Scaling(), S.SelectionModel()
    zz = np.linspace(0.1, 0.8, 36)
    for T, col, lab in ((80.0, PS.BLUE_300, "80 s"), (170.0, PS.BLUE, "170 s (median)"), (800.0, PS.BLUE_750, "800 s")):
        mm = [detection_mass(zi, T, cosmo, sr, sel) for zi in zz]
        ax.plot(zz, mm, color=col, lw=1.6, label=f"50% detected, exposure {lab}")
    ax.set(yscale="log", xlim=(0.1, 0.8), ylim=(1.5e13, 3e15), xlabel="photometric redshift",
           ylabel="M500 (Msun)")
    ax.set_title("What eROSITA can see")
    ax.legend(fontsize=7.2, loc="lower right")
    ax = axes[1]
    lnm = res["masses"]["lnm"]
    ax.plot(np.exp(lnm), res["masses"]["sum_pdf"], color=PS.BLUE, lw=1.6, label="sum of our posteriors")
    h, e = np.histogram(np.log(m[np.isfinite(m) & (m > 0)]), bins=40, range=(lnm[0], lnm[-1]))
    ax.stairs(h / np.diff(e), np.exp(e), color=PS.ORANGE, lw=1.2, label="catalogue M500 values")
    ax.set(xscale="log", xlim=(5e12, 5e15), xlabel="M500 (Msun)", ylabel="clusters per unit ln M")
    ax.set_title("The sample's mass distribution")
    vc = res["masses"]["vs_catalogue"]
    ax.text(0.03, 0.95, f"our posterior medians / catalogue M500:\nmedian {vc['median_ratio']:.2f}, "
                        f"68% range {vc['p16']:.2f}-{vc['p84']:.2f}", transform=ax.transAxes, va="top", fontsize=7.5,
            color=PS.INK_SECONDARY)
    ax.legend(fontsize=7.5, loc="upper right")
    PS.caption(fig, "Left: lines are the mass at which half of all clusters are detected, from our selection model and the published "
               "count rate-mass relation (N_H = 3.5e20). Right: the posterior uses the mass function as prior; the catalogue's values "
               "include an additional correction by Bulbul et al. (2024).", y=-0.04)
    fig.tight_layout()
    return PS.save(fig, OUT / "08_masses")


# ---------------------------------------------------------------------------
# 09  dark energy's brake
# ---------------------------------------------------------------------------


def fig_brake(res):
    b = res["brake"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    ax = axes[0]
    a = b["a"]
    z = 1 / a - 1
    ratio = b["Da_lcdm"] / b["Da_noL"]
    ax.plot(1 + z, ratio, color=PS.BLUE, lw=1.9)
    ax.axhline(1.0, color=PS.CONTEXT, lw=1.4)
    ax.text(28, 1.004, "a universe without dark energy", fontsize=7.5, color=PS.MUTED, va="bottom")
    om, ol = 0.3158, 1 - 0.3158
    z_eq = (ol / om) ** (1 / 3) - 1
    z_acc = (2 * ol / om) ** (1 / 3) - 1
    for zz, txt in ((z_acc, f"expansion starts to\naccelerate (z = {z_acc:.2f})"),
                    (z_eq, f"dark energy outweighs\nmatter (z = {z_eq:.2f})")):
        ax.axvline(1 + zz, color=PS.GRID, lw=1)
        ax.text(1 + zz, 0.80 if zz > 0.5 else 0.835, txt, fontsize=7, color=PS.INK_SECONDARY, ha="right")
    r0 = b["ratio_since_cmb"][0]
    ax.plot(1.0, r0, marker="o", ms=5, color=PS.BLUE, mec=PS.SURFACE)
    ax.text(1.0, r0 - 0.008, f"today: {r0:.2f}", fontsize=8, color=PS.INK, ha="right", va="top")
    ax.set(xscale="log", xlim=(30, 0.93), ylim=(0.76, 1.02), xlabel="1 + redshift (time runs to the right)",
           ylabel="growth since z = 1100, relative to no dark energy")
    ax.set_title("Dark energy slows the growth of structure")
    ax = axes[1]
    M = b["M"]
    for k, (zv, ls) in enumerate(zip(b["zz"], ("-", "--", ":"))):
        ax.plot(M, b["n_gt_noL"][k], color=PS.CONTEXT, lw=1.8, ls=ls)
        ax.plot(M, b["n_gt_lcdm"][k], color=PS.BLUE, lw=1.6, ls=ls, label=f"z = {zv:.1f}")
    ax.set(xscale="log", yscale="log", xlim=(1e14, 4e15), ylim=(1e-2, 2e5), xlabel="halo mass M200m (Msun)",
           ylabel="halos heavier than M, per comoving Gpc^3")
    ax.set_title("Same seeds, fewer giants")
    handles = [Line2D([], [], color=PS.BLUE, lw=1.6, label="with Lambda"), Line2D([], [], color=PS.CONTEXT, lw=1.8,
                                                                                 label="without")]
    handles += [Line2D([], [], color=PS.INK_SECONDARY, lw=1.2, ls=ls, label=f"z = {zv:.1f}")
                for zv, ls in zip(b["zz"], ("-", "--", ":"))]
    ax.legend(handles=handles, fontsize=7.5, loc="lower left")
    ax.text(0.97, 0.95, f"sigma_8 today: {b['sigma8_lcdm']:.2f} with Lambda,\n{b['sigma8_noL_phys']:.2f} without\n"
                        f"x{b['boost_gt_1e+15'][0]:.1f} more halos above 1e15 Msun",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=PS.INK_SECONDARY)
    PS.caption(fig, "Two flat universes with Planck 2018's matter and baryon densities and primordial fluctuations, identical at "
               "recombination; one has a cosmological constant. Linear growth from the growth equation (radiation included); "
               "abundances from the Tinker et al. (2008) mass function.", y=-0.04)
    fig.tight_layout()
    return PS.save(fig, OUT / "09_dark_energy_brake")


# ---------------------------------------------------------------------------
# 10  robustness
# ---------------------------------------------------------------------------


def fig_robustness(res):
    base = res["grid"]
    rows = [("baseline (published relation fixed)", base["S8"], base["sigma8"])]
    rows += [(v["label"], v["S8"], v["sigma8"]) for v in res["variants"]]
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.32 * len(rows) + 1.4), sharey=True)
    for ax, k, name in ((axes[0], 1, "S_8"), (axes[1], 2, "sigma_8")):
        b = rows[0][k]
        ax.axvspan(b[1], b[2], color=PS.WASH, zorder=0)
        for i, row in enumerate(rows):
            v = row[k]
            y = len(rows) - 1 - i
            ax.errorbar(v[0], y, xerr=[[v[0] - v[1]], [v[2] - v[0]]], ls="", marker="o", ms=4.5,
                        color=PS.BLUE if i else PS.INK, mec=PS.SURFACE, mew=0.5, elinewidth=1.1)
        ax.set_xlabel(name)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(range(len(rows)))
    axes[0].set_yticklabels([r[0] for r in rows[::-1]], fontsize=7.5)
    axes[0].set_title("What if one ingredient changes?")
    PS.caption(fig, "Each row repeats the consistency scan with one change. Shaded: baseline 68% interval.", y=-0.01)
    fig.tight_layout()
    return PS.save(fig, OUT / "10_robustness")


# ---------------------------------------------------------------------------
# 00  summary
# ---------------------------------------------------------------------------


def fig_summary(res):
    fp = F.build()
    df = K.load_cosmology_sample(fp)
    primary = K.load_primary()
    fig = plt.figure(figsize=(14, 10.2))
    ax = fig.add_axes([0.14, 0.50, 0.72, 0.46], projection="mollweide")
    draw_sky(ax, fp, df, primary, None, small=True)
    ax.set_title("5,259 galaxy clusters found by eROSITA in half the sky (colour: redshift, darker = farther)", pad=8,
                 loc="center")
    ax2 = fig.add_axes([0.05, 0.06, 0.25, 0.36])
    cols = {0.76: PS.BLUE_150, 0.81: PS.BLUE_300, 0.86: PS.BLUE_550}
    curves = [(f"sigma_8 = {s8:.2f}", N, cols[s8], 1.8) for s8, N in sorted(res["clumpy"].items())]
    draw_counts_z(ax2, res, curves=curves)
    ax2.set_title("Clumpier universe, more clusters")
    ax2.legend(fontsize=7, loc="upper right")
    ax3 = fig.add_axes([0.37, 0.06, 0.25, 0.36])
    draw_om_s8(ax3, res, compact=True)
    ax3.set_title("Matter density vs clumpiness")
    ax4 = fig.add_axes([0.80, 0.06, 0.18, 0.36])
    draw_s8(ax4, res, compact=True)
    ax4.set_title(SUMMARY_S8_TITLE)
    return PS.save(fig, OUT / "00_summary")


#: Set after looking at the result (see README): what the S_8 panel shows.
SUMMARY_S8_TITLE = "S_8: no deficit of clumpiness"


def make_all(res) -> list:
    paths = []
    for f in (fig_summary, fig_sky, fig_survey, fig_counts, fig_clumpy, fig_om_s8, fig_s8, fig_mass_scale,
              fig_masses, fig_brake, fig_robustness):
        paths += f(res)
    return paths
