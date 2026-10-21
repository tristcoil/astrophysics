"""Every figure, drawn from the cached analysis results (data/derived/results.pkl).

    python figures.py            # all figures
    python figures.py 3 7        # only figures 3 and 7
"""

from __future__ import annotations

import pickle
import sys

import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D

import literature as lit
import physics
from config import FIGURES, LUMINOUS_M_CUT, RESULTS_PICKLE, SEED_CLASSES
from plotstyle import (AXIS, BLUE_RAMP, C_MODEL, C_OURS, C_PUBLISHED, GRID, INK, INK_SECONDARY,
                       LINE_COLOURS, LINE_LABELS, MUTED, ORANGE, ORDINAL_BLUE, PROGRAMME_COLOURS, SURFACE,
                       annotate, caption, save, use_style)

use_style()
SEQ = LinearSegmentedColormap.from_list("seq", BLUE_RAMP[1:])
SEQ.set_bad((0, 0, 0, 0))
GREY_SEQ = LinearSegmentedColormap.from_list("grey", ["#eeede8", "#b4b2aa", "#6f6d66"])


def load_results() -> dict:
    with open(RESULTS_PICKLE, "rb") as fh:
        return pickle.load(fh)


def density(ax, x, y, xbins, ybins, cmap=SEQ, vmin=1, rasterized=True, **kw):
    h, xe, ye = np.histogram2d(x, y, [xbins, ybins])
    h = np.ma.masked_less(h.T, vmin)
    return ax.pcolormesh(xe, ye, h, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=max(h.max(), vmin + 1)),
                         rasterized=rasterized, **kw)


def age_top_axis(ax, zs=(0, 0.5, 1, 2, 3, 5, 10, 20)):
    """A redshift scale along the top of an axis whose x is cosmic age in Gyr."""
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    ages = physics.age_myr(np.array(zs, float)) / 1e3
    keep = (ages >= ax.get_xlim()[0]) & (ages <= ax.get_xlim()[1])
    top.set_xticks(ages[keep])
    top.set_xticklabels([f"{z:g}" for z in np.array(zs)[keep]])
    top.set_xlabel("redshift z", color=INK_SECONDARY)
    top.grid(False)
    for s in ("top",):
        top.spines[s].set_visible(True)
        top.spines[s].set_color(AXIS)
    top.spines["right"].set_visible(False)
    return top


def panel_label(ax, text, x=-0.02, y=1.02):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK,
            ha="right", va="bottom")


# ---------------------------------------------------------------------------
# 1. Census: sky and redshifts
# ---------------------------------------------------------------------------

def fig_census(res):
    c = res["census"]
    fig = plt.figure(figsize=(13.0, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.2)
    ax = fig.add_subplot(gs[0], projection="mollweide")
    nside = 64
    ok = c["z"] > 0
    pix = hp.ang2pix(nside, c["ra"][ok].astype(float), c["dec"][ok].astype(float), lonlat=True)
    counts = np.bincount(pix, minlength=hp.nside2npix(nside)).astype(float)
    dens = counts / hp.nside2pixarea(nside, degrees=True)
    lon_e = np.linspace(-np.pi, np.pi, 721)
    lat_e = np.linspace(-np.pi / 2, np.pi / 2, 361)
    lon_c, lat_c = 0.5 * (lon_e[1:] + lon_e[:-1]), 0.5 * (lat_e[1:] + lat_e[:-1])
    LON, LAT = np.meshgrid(lon_c, lat_c)
    ra = (180.0 - np.degrees(LON)) % 360.0
    grid = dens[hp.ang2pix(nside, ra, np.degrees(LAT), lonlat=True)]
    grid = np.ma.masked_less(grid, 0.5)
    m = ax.pcolormesh(lon_e, lat_e, grid, cmap=SEQ, norm=LogNorm(vmin=1, vmax=200), rasterized=True)
    # the Milky Way's plane, where dust hides everything behind it
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    gl = SkyCoord(l=np.linspace(0, 360, 721) * u.deg, b=np.zeros(721) * u.deg, frame="galactic").icrs
    glon = np.radians(180.0 - gl.ra.deg)
    glon = (glon + np.pi) % (2 * np.pi) - np.pi
    order = np.argsort(glon)
    ax.plot(glon[order], np.radians(gl.dec.deg[order]), color=MUTED, lw=0.8)
    ax.text(np.radians(180 - 290), np.radians(12), "Milky Way\nplane", color=MUTED, fontsize=7.5,
            ha="center")
    ticks = np.radians(np.arange(-150, 151, 30))
    ax.set_xticks(ticks)
    ax.set_xticklabels([])
    for d in (-120, -60, 0, 60, 120):
        ax.text(np.radians(d), np.radians(-36), f"RA {(180 - d) % 360:.0f}°", fontsize=7,
                color=MUTED, ha="center", va="center", zorder=10)
    ax.tick_params(axis="y", labelsize=7, colors=MUTED)
    ax.grid(True, color=GRID, lw=0.5)
    cb = fig.colorbar(m, ax=ax, orientation="horizontal", fraction=0.05, pad=0.06, aspect=40)
    cb.set_label("quasars per square degree", color=INK_SECONDARY)
    cb.outline.set_visible(False)
    ax.set_title(f"{c['n']:,} quasars of SDSS DR16Q on the sky", loc="left")

    ax2 = fig.add_subplot(gs[1])
    edges = c["z_edges"]
    base = np.zeros(len(edges) - 1)
    order = ["eBOSS CORE", "BOSS CORE", "SDSS-I/II", "Other BOSS/eBOSS"]
    for prog in order:
        h = c["z_hist"][prog]
        ax2.fill_between(edges[:-1], base, base + h, step="post", color=PROGRAMME_COLOURS[prog],
                         lw=0, label=prog)
        ax2.stairs(base + h, edges, color=SURFACE, lw=0.6)
        base = base + h
    ax2.set_xlim(0, 5.5)
    ax2.set_ylim(0, base.max() * 1.12)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel(f"quasars per Δz = {edges[1] - edges[0]:.2f}")
    leg = ax2.legend(loc="upper right", title="chosen by", title_fontsize=8)
    leg._legend_box.align = "left"
    rows = {r["programme"]: r for r in c["by_programme"]}
    annotate(ax2, f"eBOSS CORE aimed at\n0.9 < z < 2.2 ({rows['eBOSS CORE']['n']:,})",
             (1.55, base[int(1.55 / 0.05)] * 0.62), (0.05, base.max() * 0.93), fontsize=7.5)
    kb = int(2.45 / 0.05)
    blue_top = c["z_hist"]["eBOSS CORE"][kb]
    orange_mid = blue_top + 0.5 * c["z_hist"]["BOSS CORE"][kb]
    annotate(ax2, f"BOSS CORE aimed at\nz > 2.15 ({rows['BOSS CORE']['n']:,})",
             (2.45, orange_mid), (2.95, base.max() * 0.55), fontsize=7.5)
    hz = c["n_z_gt"]
    ax2.text(5.45, base.max() * 0.18, f"z > 4: {hz[4.0]:,}\nz > 5: {hz[5.0]:,}\nz > 6: {hz[6.0]:,}",
             ha="right", va="bottom", fontsize=7.5, color=INK_SECONDARY)
    ax2.set_title("A census shaped by who was looking for what", loc="left")
    caption(fig, "Left: sky density of DR16Q quasars (HEALPix, 0.84 deg² pixels; RA 180° at the "
            "centre, increasing to the left). The densest regions are eBOSS's footprint; the "
            "sparser cap is SDSS-I/II and BOSS only. Right: redshift distribution, coloured by the "
            "programme whose target bits selected each quasar (DR16Q, Lyke et al. 2020; redshifts "
            "from Wu & Shen 2022). The bumps are targeting choices, not features of the universe: "
            "eBOSS CORE chased 0.9 < z < 2.2 for clustering, BOSS CORE chased z > 2.15 for the "
            "Lyman-alpha forest. Counting quasars per volume therefore needs one uniformly selected "
            "subsample (figure 5).", y=0.02)
    return save(fig, FIGURES / "01_quasar_census")


# ---------------------------------------------------------------------------
# 2. Luminosity vs redshift
# ---------------------------------------------------------------------------

def fig_luminosity(res):
    from analysis import kcorr_r06
    lz = res["lum_z"]
    ph = res["photometry"]
    k06 = kcorr_r06()
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ok = np.isfinite(lz["mi"]) & (lz["mi"] < -18) & (lz["mi"] > -32)
    m = density(ax, lz["z"][ok], lz["mi"][ok], np.linspace(0, 7.2, 241), np.linspace(-31, -18, 261))
    z = np.linspace(0.05, 7.6, 600)
    dm = physics.distance_modulus(z, "d3")
    kz = k06(np.minimum(z, 5.49))

    def mlim(mag, colour=0.0):
        return mag - colour - dm - kz

    zlo = z < 2.9
    ax.plot(z[zlo], mlim(19.1)[zlo], color=INK, lw=1.4)
    ax.plot(z[~zlo & (z < 5.5)], mlim(20.2)[~zlo & (z < 5.5)], color=INK, lw=1.4)
    annotate(ax, "SDSS-I/II limit: i = 19.1,\nthen 20.2 above z = 2.9", (2.4, mlim(19.1)[np.argmin(abs(z - 2.4))]),
             (2.75, -30.3), fontsize=8, color=INK)
    # BOSS: g <= 22.0 or r <= 21.85 (eBOSS: g or r < 22), so the fainter of the two
    # limits (in i, via our composite's colours) applies at each redshift
    gi = np.interp(z, ph["zgrid"], ph["g_minus_i"], right=np.nan)
    ri = np.interp(z, ph["zgrid"], ph["r_minus_i"], right=np.nan)
    gl = np.maximum(mlim(22.0, gi), mlim(21.85, ri))
    ax.plot(z, gl, color=INK_SECONDARY, lw=1.2)
    annotate(ax, "BOSS / eBOSS limit: g ≈ 22 or r ≈ 21.85", (0.6, gl[np.argmin(abs(z - 0.6))]),
             (0.15, -18.6), fontsize=8)
    ax.plot(z[zlo], mlim(15.0)[zlo], color=MUTED, lw=0.9)
    ax.text(0.33, mlim(15.0)[np.argmin(abs(z - 0.33))] - 0.25, "bright limit i = 15", color=MUTED,
            fontsize=7.5, rotation=-22)
    # the luminous cut used for the quasar era, over its complete range
    era = res["qlf"]["era"][LUMINOUS_M_CUT]
    zc = era[era["complete"]]
    ax.plot([zc["z_lo"].min(), zc["z_hi"].max()], [LUMINOUS_M_CUT] * 2, color=C_MODEL, lw=1.6)
    annotate(ax, f"M_i < {LUMINOUS_M_CUT}: complete over\n{zc['z_lo'].min():.1f} < z < "
             f"{zc['z_hi'].max():.1f}, used in figure 6", (1.6, LUMINOUS_M_CUT), (0.95, -30.25),
             color=C_MODEL, fontsize=7.5)
    # z > 5.7 quasars from the literature
    hz = res["growth"]["hz"]
    ax.scatter(hz["z"], hz["M1450"] - 1.486, s=14, color=C_PUBLISHED, edgecolor=SURFACE, lw=0.5,
               zorder=5)
    ax.text(6.55, -24.0, f"{len(hz)} quasars at z > 5.7\n(Farina et al. 2022\ncompilation)",
            fontsize=7.5, color=C_PUBLISHED, ha="center", va="top")
    ax.set_xlim(0, 7.8)
    ax.set_ylim(-18, -31)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("absolute magnitude M_i(z=2)   (brighter ↑)")
    sec = ax.secondary_yaxis("right", functions=(lambda mi: physics.luminosity_2500(mi) + np.log10(2.99792458e18 / 2500),
                                                 lambda l: -2.5 * (l - np.log10(2.99792458e18 / 2500)
                                                                   - np.log10(4 * np.pi * (3.0857e19) ** 2)) - 48.60 - 2.5 * np.log10(3.0)))
    sec.set_ylabel("log νL_ν at 2500 Å  (erg/s)", color=INK_SECONDARY)
    cb = fig.colorbar(m, ax=ax, pad=0.09, fraction=0.035)
    cb.set_label("quasars per cell", color=INK_SECONDARY)
    cb.outline.set_visible(False)
    ax.set_title("Only the brightest quasars are seen far away", loc="left")
    caption(fig, "Each DR16Q quasar's luminosity (M_i(z=2), the i-band absolute magnitude K-corrected "
            "to z = 2 as in Richards et al. 2006, in D3's cosmology) against its redshift. Every "
            "survey sees down to a fixed apparent brightness, so the faintest quasar it can find "
            "grows more luminous with distance: the lower edge of the cloud is the survey limit, not "
            "a property of quasars. The BOSS/eBOSS curve converts g = 22 or r = 21.85 (whichever "
            "reaches fainter) into M_i with our composite spectrum's g - i and r - i colours. Right axis: the same luminosity as νL_ν at 2500 Å "
            "(Richards et al. 2006, eq. 4).", y=-0.01)
    return save(fig, FIGURES / "02_luminosity_redshift")


# ---------------------------------------------------------------------------
# 3. The composite spectrum
# ---------------------------------------------------------------------------

LINE_IDS = [
    (1033, "Lyβ/O VI"), (1216, "Lyα"), (1240, "N V"), (1400, "Si IV"), (1549, "C IV"),
    (1640, "He II"), (1909, "C III]"), (2326, "C II]"), (2798, "Mg II"), (3426, "[Ne V]"),
    (3727, "[O II]"), (4102, "Hδ"), (4340, "Hγ"), (4861, "Hβ"), (5007, "[O III]"),
    (5876, "He I"), (6563, "Hα"),
]


def fig_composite(res):
    c = res["composite"]
    w = c["wave"]
    fig = plt.figure(figsize=(13.0, 8.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[3.2, 0.9, 0.9], hspace=0.08, left=0.06, right=0.64)
    gsr = fig.add_gridspec(3, 1, hspace=0.42, left=0.70, right=0.99)
    ax = fig.add_subplot(gs[0, 0])
    good = c["n_per_pixel"] >= 10
    f = np.where(good, c["median"], np.nan)
    ax.plot(c["vb_wave"], c["vb_flux"], color=C_PUBLISHED, lw=0.8, label=f"Vanden Berk et al. 2001 ({lit.VB01['n_spectra']:,} spectra)")
    ax.plot(w, f, color=INK, lw=0.7, label=f"this project ({c['n_spectra']:,} spectra)")
    pl = c["powerlaw_median"]
    ax.plot(w, pl, color=C_MODEL, lw=1.0, label=f"power law, α_ν = {c['alpha_median']:.2f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(650, 9600)
    ax.set_ylim(0.08, 30)
    ax.set_ylabel("f_λ  (normalised at 3020-3100 Å)")
    ax.tick_params(labelbottom=False)
    for lam, name in LINE_IDS:
        if name in ("N V", "He I"):
            continue  # blended with their neighbours at this scale
        k = np.argmin(abs(w - lam))
        top = np.nanmax(f[max(k - 6, 0):k + 7])
        if not np.isfinite(top):
            continue
        ax.text(lam, top * 1.18, name, rotation=90, fontsize=7, ha="center", va="bottom",
                color=INK_SECONDARY)
    annotate(ax, "Lyman-α forest: hydrogen clouds\nbetween us and each quasar", (960, 1.2),
             (1350, 0.17), fontsize=7.5)
    annotate(ax, "host-galaxy starlight\nabove the power law", (7400, 0.43), (5200, 0.16),
             fontsize=7.5)
    ax.legend(loc="upper right")
    ax.set_title("The average quasar, from 4,446 SDSS spectra", loc="left")
    xt = [700, 1000, 1500, 2000, 3000, 5000, 7000, 9000]

    axr = fig.add_subplot(gs[1, 0], sharex=ax)
    ratio = c["ratio_to_vb"]
    axr.axhspan(0.95, 1.05, color=GRID, lw=0)
    axr.plot(w, np.where(good, ratio, np.nan), color=INK, lw=0.6)
    axr.axhline(1, color=AXIS, lw=0.8)
    axr.set_ylim(0.75, 1.45)
    axr.set_ylabel("ours ÷ VB01")
    axr.tick_params(labelbottom=False)
    axn = fig.add_subplot(gs[2, 0], sharex=ax)
    axn.fill_between(w, 0, c["n_per_pixel"], color=BLUE_RAMP[3], lw=0)
    axn.set_ylabel("spectra\nper pixel")
    axn.set_xlabel("rest-frame wavelength (Å)")
    axn.set_xticks(xt)
    axn.set_xticklabels([str(x) for x in xt])
    axn.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    # the luminosity sequence
    bins = c["luminosity_bins"]
    zooms = [((1420, 1720), "C IV region", (1440, 1460)), ((2650, 2950), "Mg II", (2650, 2670)),
             ((4700, 5100), "Hβ and [O III]", (4740, 4760))]
    for j, ((lo, hi), title, nw) in enumerate(zooms):
        a = fig.add_subplot(gsr[j])
        for b, col in zip(bins, ORDINAL_BLUE[::-1]):
            ok = (b["n_pix"] >= 20) & (w >= lo) & (w <= hi)
            if ok.sum() < 0.8 * (hi - lo):
                continue
            fn = b["flux"] / np.nanmedian(b["flux"][(w >= nw[0]) & (w <= nw[1])])
            lab = f"M_i {b['m_lo']:.0f} to {b['m_hi']:.0f}" if b["m_lo"] > -29 else f"M_i < {b['m_hi']:.0f}"
            if b["m_hi"] > -18.5:
                lab = f"M_i > {b['m_lo']:.0f}"
            a.plot(w[ok], fn[ok], color=col, lw=1.0, label=lab)
        a.set_xlim(lo, hi)
        a.set_title(title, loc="left", fontsize=9)
        a.tick_params(labelsize=7)
        if j == 0:
            a.text(0.02, 0.96, "brighter quasars:\nweaker C IV\n(the Baldwin effect)", transform=a.transAxes,
                   fontsize=7, va="top", color=INK_SECONDARY)
        if j == 1:
            handles = []
            for b, col in zip(bins, ORDINAL_BLUE[::-1]):
                if b["m_lo"] < -29:
                    lab = f"M_i < {b['m_hi']:.0f}"
                elif b["m_hi"] > -18.5:
                    lab = f"M_i > {b['m_lo']:.0f}"
                else:
                    lab = f"{b['m_lo']:.0f} < M_i < {b['m_hi']:.0f}"
                handles.append(Line2D([], [], color=col, lw=1.4, label=f"{lab}  ({b['n']:,})"))
            a.legend(handles=handles, fontsize=6.5, loc="upper right", title="luminosity (spectra)",
                     title_fontsize=7)
        if j == 2:
            a.set_xlabel("rest wavelength (Å)")
    caption(fig, f"Median composite of {c['n_spectra']:,} BOSS/eBOSS spectra chosen evenly in log(1+z) "
            f"from z = 0.06 to 4.6 (left), stacked as Vanden Berk et al. (2001) did: dereddened, "
            f"shifted to the rest frame with Wu & Shen's systemic redshifts, rebinned to 1 Å, scaled "
            f"in their overlap and median-combined. S/N per pixel ≈ {c['snr_2800']:.0f} at Mg II. "
            f"From 1,250 to 3,500 Å ours and VB01's agree to within "
            f"{min(abs(1 - c['vb_comparison'][k]['median_ratio']) for k in ('1250-2000', '2000-3500')) * 100:.0f}-"
            f"{max(abs(1 - c['vb_comparison'][k]['median_ratio']) for k in ('1250-2000', '2000-3500')) * 100:.0f}%. "
            f"Redward, ours rises 10-20% higher with stronger narrow lines: BOSS reaches fainter quasars, whose "
            f"host galaxies add starlight. Blueward of Lyα ours is 10-40% higher: the forest is intervening "
            f"gas whose depth depends on the redshift mix stacked there (and BOSS fibres placed for 4000 Å "
            f"make blue spectra 3-6% bluer). Right: the same stack split by luminosity (normalised beside each "
            f"line). Our most luminous quasars (M_i < -25) give α_ν = {c['alpha_gm_luminous']:.2f} "
            f"(geometric mean), VB01's {lit.VB01['alpha_nu_gm']}.", y=0.03)
    return save(fig, FIGURES / "03_composite_spectrum")


# ---------------------------------------------------------------------------
# 4. Colours, the stellar locus and the selection function
# ---------------------------------------------------------------------------

def fig_colours(res):
    ph = res["photometry"]
    zg = ph["zgrid"]
    tr = ph["track"]
    fig = plt.figure(figsize=(13.0, 8.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.0, 1.0], height_ratios=[1.0, 1.0],
                          hspace=0.32, wspace=0.28)
    ax = fig.add_subplot(gs[:, 0])
    st = ph["star_colours"]
    density(ax, st[:, 1], st[:, 0], np.linspace(-0.6, 2.2, 141), np.linspace(-0.6, 4.5, 256),
            cmap=GREY_SEQ)
    ok = np.isfinite(tr["u-g"]) & np.isfinite(tr["g-r"]) & (tr["u-g"] < 4.5)
    ax.plot(tr["g-r"][ok], tr["u-g"][ok], color=C_MODEL, lw=1.6)
    obs = ph["observed"]
    om = np.isfinite(obs["u-g"]["median"]) & np.isfinite(obs["g-r"]["median"])
    ax.scatter(obs["g-r"]["median"][om], obs["u-g"]["median"][om], s=9, color=C_OURS, zorder=4,
               edgecolor=SURFACE, lw=0.4)
    for zz, off in ((0.5, (-34, 6)), (2.0, (10, -12)), (2.5, (8, -4)), (2.7, (8, -2)), (3.0, (8, -2)),
                    (3.5, (8, -2))):
        k = np.argmin(abs(zg - zz))
        if np.isfinite(tr["u-g"][k]) and tr["u-g"][k] < 4.4:
            ax.annotate(f"z={zz:g}", (tr["g-r"][k], tr["u-g"][k]), xytext=off,
                        textcoords="offset points", fontsize=7.5, color=INK)
            ax.plot(tr["g-r"][k], tr["u-g"][k], "o", ms=3.5, color=C_MODEL)
    ax.text(1.25, 2.7, "ordinary stars\n(61,045 in a 100 deg²\npatch of sky)", fontsize=7.5,
            color=INK_SECONDARY)
    ax.set_xlim(-0.6, 2.0)
    ax.set_ylim(-0.6, 4.4)
    ax.set_xlabel("g - r")
    ax.set_ylabel("u - g")
    ax.legend(handles=[Line2D([], [], color=C_MODEL, lw=1.6, label="our composite moved through redshift"),
                       Line2D([], [], color=C_OURS, marker="o", ls="", ms=4,
                              label="DR16Q quasars: median colour per Δz = 0.05")],
              loc="upper left")
    ax.set_title("Quasars hiding among the stars", loc="left")

    # colours vs redshift, synthetic vs observed
    sub = gs[0, 1].subgridspec(2, 2, hspace=0.1, wspace=0.32)
    for k, col in enumerate(("u-g", "g-r", "r-i", "i-z")):
        a = fig.add_subplot(sub[k // 2, k % 2])
        o = obs[col]
        a.fill_between(o["z"], o["p16"], o["p84"], color=BLUE_RAMP[1], lw=0)
        a.plot(o["z"], o["median"], color=C_OURS, lw=1.0)
        a.plot(zg, tr[col], color=C_MODEL, lw=1.0)
        a.set_xlim(0, 5)
        a.set_ylim(-0.6, 3.0 if col in ("u-g", "g-r") else 1.0)
        a.text(0.04, 0.9, col, transform=a.transAxes, fontsize=8, fontweight="bold")
        a.tick_params(labelsize=7)
        if k >= 2:
            a.set_xlabel("z", fontsize=8)
        else:
            a.tick_params(labelbottom=False)
    fig.text(0.445, 0.915, "Colours vs redshift", fontsize=11, color=INK, fontweight="bold")
    fig.text(0.445, 0.897, "observed (blue, 16-84% band), synthetic (orange)", fontsize=7.5,
             color=INK_SECONDARY)

    # selection function and distance from the stellar locus
    sub2 = gs[0, 2].subgridspec(2, 1, hspace=0.1)
    a1 = fig.add_subplot(sub2[0])
    a1.plot(zg, ph["selfunc_marginal"], color=C_PUBLISHED, lw=1.4)
    a1.set_ylim(0, 1.05)
    a1.set_xlim(0, 5)
    a1.set_ylabel("fraction\nselected", fontsize=8)
    a1.tick_params(labelbottom=False, labelsize=7)
    a1.set_title("Selection vs nearness to stars", loc="left")
    a2 = fig.add_subplot(sub2[1], sharex=a1)
    a2.plot(zg, ph["locus_distance_ugri"], color=C_MODEL, lw=1.4)
    a2.plot(zg, ph["locus_distance_griz"], color=C_OURS, lw=1.4)
    a2.text(0.3, 0.42, "in ugri colours\n(low-z branch)", color=C_MODEL, fontsize=7)
    a2.text(3.6, 0.05, "in griz\n(high-z\nbranch)", color=C_OURS, fontsize=7)
    a2.set_ylim(0, 0.75)
    a2.set_ylabel("colour distance\nto stars (mag)", fontsize=8)
    a2.set_xlabel("z", fontsize=8)
    a2.tick_params(labelsize=7)
    for zz in (2.7, 3.5):
        for a in (a1, a2):
            a.axvline(zz, color=AXIS, lw=0.8)
    a1.text(2.62, 0.1, "2.7", fontsize=7, color=INK_SECONDARY, ha="right")
    a1.text(3.55, 0.3, "3.5", fontsize=7, color=INK_SECONDARY)
    a1.text(0.1, 0.12, "Richards et al. 2006,\n15 < i < 19.1", fontsize=7, color=C_PUBLISHED)

    # K-corrections
    a3 = fig.add_subplot(gs[1, 1:])
    a3.plot(ph["kz"], ph["k_r06"], color=C_PUBLISHED, lw=1.6, label="Richards et al. 2006 (their DR1 composite)")
    a3.plot(ph["kz"], ph["k_ours"], color=C_OURS, lw=1.2, label="ours (luminous-quasar template)")
    from photometry import kcorr_continuum
    a3.plot(ph["kz"], kcorr_continuum(ph["kz"]), color=MUTED, lw=0.9, label="continuum only (α_ν = -0.5)")
    a3.set_xlim(0, 5.5)
    a3.set_xlabel("redshift z")
    a3.set_ylabel("K(z), i band, pegged at z = 2 (mag)")
    a3.legend(loc="lower left")
    kc = ph["k_comparison"]
    a3.set_title("The K-correction: from observed i to M_i(z=2)", loc="left")
    for zz, name in ((0.27, "Hα"), (0.55, "Hβ, [O III]"), (1.67, "Mg II"), (3.82, "C IV")):
        a3.annotate(name, (zz, np.interp(zz, ph["kz"], ph["k_r06"])), xytext=(0, -16),
                    textcoords="offset points", fontsize=7, ha="center", color=INK_SECONDARY,
                    arrowprops=dict(arrowstyle="-", color=INK_SECONDARY, lw=0.6))
    ag = ph["agreement"]
    caption(fig, "Left: the u - g / g - r plane. Grey: stars from a SkyServer patch at high Galactic "
            "latitude. Orange: our composite spectrum moved to each redshift through the SDSS filters, "
            "with the mean intergalactic absorption of Madau (1995) applied (a model). Blue: the "
            "median colours of DR16Q quasars brighter than i = 19.5. Near z ≈ 2.7 the quasar track "
            "runs straight through the stars, so a colour cut that rejects stars also rejects quasars. "
            "Top right: the completeness of the SDSS colour selection (Richards et al. 2006, Table 1, "
            "averaged over 15 < i < 19.1) against how close our track comes to the stars (distance to the "
            "20th-nearest star). Its deep dip at z ≈ 2.7 is where the track, in the ugri colours that select "
            "z < 3 quasars, enters the stellar locus; in griz it runs along the stars over 2.8 < z < 3.3, "
            "just below the selection's second dip at 3.5. Synthetic and observed "
            f"colours agree to {ag['g-r 0.3-2.2']:.2f}-{ag['u-g 0.3-2.2']:.2f} mag (median) at z < 2.2. "
            "Bottom: our emission-line K-correction vs Richards et al.'s: median differences "
            f"{kc['0.5-2.2']['median']:+.2f} mag (0.5 < z < 2.2) and {kc['2.2-4.0']['median']:+.2f} "
            "(2.2 < z < 4). The dips are emission lines passing through the i band.", y=0.03)
    return save(fig, FIGURES / "04_colours_and_selection")


# ---------------------------------------------------------------------------
# 5. The luminosity function
# ---------------------------------------------------------------------------

def fig_qlf(res):
    q = res["qlf"]
    rt = q["r06_table"]          # Richards et al.'s bins with our rebuilt values
    d7 = q["dr7"]
    zbins = sorted(rt["zmid"].unique())
    fig = plt.figure(figsize=(13.0, 9.0))
    gs = fig.add_gridspec(3, 4, hspace=0.06, wspace=0.05, top=0.9)
    axes = []
    for k in range(11):
        share = axes[0] if axes else None
        axes.append(fig.add_subplot(gs[k // 4, k % 4], sharex=share, sharey=share))
    axes.append(fig.add_subplot(gs[2, 3]))
    for k in range(11):
        if k % 4:
            axes[k].tick_params(labelleft=False)
        if k < 8:
            axes[k].tick_params(labelbottom=False)
    for k, zm in enumerate(zbins):
        a = axes[k]
        t = rt[rt["zmid"] == zm]
        filled = t["Fill"] == 1
        sig = (t["errPhi"] * 1e-9) / (10 ** t["logPhi"]) / np.log(10)
        a.errorbar(t["iMAG"][filled], t["logPhi"][filled], yerr=sig[filled], fmt="o", ms=5,
                   mfc="none", mec=C_PUBLISHED, ecolor=C_PUBLISHED, elinewidth=0.8, lw=0, mew=1.0)
        a.errorbar(t["iMAG"][~filled], t["logPhi"][~filled], yerr=sig[~filled], fmt="o", ms=4,
                   mfc="none", mec=C_PUBLISHED, ecolor=C_PUBLISHED, elinewidth=0.6, lw=0, mew=0.6,
                   alpha=0.5)
        ours = t[np.isfinite(t["phi"])]
        a.plot(ours["iMAG"], np.log10(ours["phi"]), "+", color=C_OURS, ms=6, mew=1.0)
        g = d7[np.isclose(d7["z_mid"], zm)]
        g = g[g["fill"] > 0.02]
        full = g["fill"] > 0.999
        ly = np.log10(g["phi"])
        ey = g["phi_err"] / g["phi"] / np.log(10)
        a.errorbar(g["m_mid"][full], ly[full], yerr=ey[full], fmt="o", ms=3.2, color=INK,
                   elinewidth=0.7, lw=0)
        a.errorbar(g["m_mid"][~full], ly[~full], yerr=ey[~full], fmt="o", ms=3.2, mfc=SURFACE,
                   mec=INK, ecolor=INK, elinewidth=0.6, lw=0, mew=0.8)
        zlo, zhi = float(g["z_lo"].iloc[0]) if len(g) else zm, float(g["z_hi"].iloc[0]) if len(g) else zm
        a.text(0.95, 0.92, f"{zlo:.2f} < z < {zhi:.2f}", transform=a.transAxes, ha="right",
               fontsize=8, color=INK)
        a.set_xlim(-22.2, -30.3)
        a.set_ylim(-9.6, -5.0)
        a.tick_params(labelsize=7)
    a = axes[11]
    rep = rt["dlog"].dropna()
    skn = q["sk12_table_noimg"]["dlog"].dropna()
    bins = np.linspace(-0.3, 0.3, 61)
    a.hist(rep, bins, histtype="step", color=C_PUBLISHED, lw=1.4)
    a.hist(skn, bins, histtype="step", color=INK, lw=1.4)
    a.set_xlim(-0.3, 0.3)
    a.set_ylim(0, a.get_ylim()[1] * 1.9)
    a.set_xlabel("ours − published (dex)", fontsize=8)
    a.set_ylabel("bins", fontsize=8)
    a.tick_params(labelsize=7, labelleft=True, labelbottom=True)
    a.yaxis.set_label_position("right")
    a.yaxis.tick_right()
    r06s, sks = q["r06_reproduction"], q["sk12_comparison_without_image_factor"]
    a.text(0.03, 0.97, f"DR3 rebuilt vs Richards+06:\n{r06s['n_matched']} bins, median "
           f"{r06s['median_dlog']:+.3f}\n\nDR7 vs Shen & Kelly 12:\n{sks['n_matched']} bins, median "
           f"{sks['median_dlog']:+.4f}", transform=a.transAxes, va="top", fontsize=7, color=INK_SECONDARY)
    for k in (0, 4, 8):
        axes[k].set_ylabel("log Φ  (Mpc⁻³ mag⁻¹)", fontsize=8)
    for k in (8, 9, 10):
        axes[k].set_xlabel("M_i(z=2)", fontsize=8)
        axes[k].tick_params(labelbottom=True)
    handles = [Line2D([], [], marker="o", ls="", mfc="none", mec=C_PUBLISHED, color=C_PUBLISHED,
                      label="Richards et al. 2006 (DR3, published)"),
               Line2D([], [], marker="+", ls="", color=C_OURS, ms=7,
                      label="our rebuild from their 15,343 quasars"),
               Line2D([], [], marker="o", ls="", color=INK, ms=4,
                      label=f"ours: DR7 uniform sample ({q['n_dr7_uniform']:,} quasars)"),
               Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=INK, ms=4,
                      label="bin only partly inside the flux limits")]
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 0.955), fontsize=8)
    fig.suptitle("How many quasars of each luminosity, per cubic megaparsec, at each epoch",
                 x=0.06, y=0.995, ha="left")
    caption(fig, "The binned quasar luminosity function (number per Mpc³ per magnitude) by the Page & "
            "Carrera (2000) 1/V_a method, with each quasar weighted by 1/(selection probability) from "
            "Richards et al. (2006, Table 1; floor 1/3, global 95% image completeness), "
            "h = 0.7, Ω_m = 0.3 as in the papers compared with. Violet: Richards et al.'s DR3 result; "
            "blue crosses: our rebuild of it from their own published sample; black: the four-times-"
            "larger DR7 uniform sample (6,248 deg²) through the same code. Last panel: the rebuild "
            "matches Richards et al. bin by bin, and the DR7 result matches Shen & Kelly (2012), who "
            "did not apply the 5% image-completeness factor, to within 0.01 dex.", y=0.02)
    return save(fig, FIGURES / "05_luminosity_function")


# ---------------------------------------------------------------------------
# 6. The quasar era
# ---------------------------------------------------------------------------

def fig_era(res):
    q = res["qlf"]
    era, pk = q["era"], q["era_peak"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.6), gridspec_kw=dict(wspace=0.25))
    styles = {-28.5: (ORDINAL_BLUE[4], "M_i < -28.5"), -27.6: (C_OURS, "M_i < -27.6")}
    for cut, (col, lab) in styles.items():
        d = era[cut]
        age = physics.age_myr(d["z_mean"]) / 1e3
        y = d["density"] * 1e9
        e = d["err"] * 1e9
        c = d["complete"].to_numpy(bool)
        ax.errorbar(age[c], y[c], yerr=e[c], fmt="o", color=col, ms=5, elinewidth=0.8, lw=0, label=lab)
        ax.errorbar(age[~c], y[~c], yerr=e[~c], fmt="o", mfc=SURFACE, mec=col, ecolor=col, ms=5,
                    elinewidth=0.6, lw=0, mew=0.9)
        last = np.where(c)[0][0]
        ax.text(age[last] + 0.12, y[last] * 0.92, lab, color=col, fontsize=8, va="top")
    r = q["r06_luminous_density_d3"]
    r = r[r["complete"]]
    rz = physics.age_myr(r["z_mean"]) / 1e3
    # Richards et al.'s bins, moved into D3's cosmology (analysis.r06_to_d3)
    ax.errorbar(rz, r["density"] * 1e9, yerr=r["err"] * 1e9, fmt="s", ms=4.5, mfc="none",
                mec=C_PUBLISHED, ecolor=C_PUBLISHED, elinewidth=0.7, lw=0, mew=1.0,
                label="Richards et al. 2006 (DR3), M_i < -27.6")
    # the parabola fitted through the complete bins
    b = pk["bins"]
    x = np.log10(1 + b["z_mean"].to_numpy())
    yv = np.log10(b["density"].to_numpy())
    sig = b["err"].to_numpy() / b["density"].to_numpy() / np.log(10)
    cpar = np.polyfit(x, yv, 2, w=1 / sig)
    zz = np.linspace(b["z_mean"].min(), b["z_mean"].max(), 200)
    ax.plot(physics.age_myr(zz) / 1e3, 10 ** np.polyval(cpar, np.log10(1 + zz)) * 1e9, color=C_MODEL,
            lw=1.4)
    ap = pk["age_at_peak_gyr"]
    ax.axvline(ap, color=AXIS, lw=0.8)
    v = pk["variants"]
    vz = np.array([val for k, val in v.items() if k.startswith(("floor", "parabola", "z <", "z >", "no floor"))])
    annotate(ax, f"peak: z = {pk['z_peak']:.2f} ± {pk['z_peak_p84'] - pk['z_peak']:.2f} "
             f"(choices: {vz.min():.2f}-{vz.max():.2f})\n{ap:.2f} billion years after the Big Bang",
             (ap, pk["density_peak_per_gpc3"]), (ap + 1.2, 140), fontsize=8, color=INK)
    ax.set_yscale("log")
    ax.set_xlim(0.0, 6.5)
    ax.set_ylim(0.08, 300)
    ax.set_xlabel("age of the universe (billion years)")
    ax.set_ylabel("luminous quasars per Gpc³ (comoving)")
    ax.legend(loc="lower left", fontsize=7.5)
    age_top_axis(ax, zs=(1, 1.5, 2, 3, 4, 5, 6, 8))
    ax.set_title("The quasar era", loc="left", pad=28)

    # indexed: quasars vs star formation
    t = np.linspace(0.4, 13.8, 400)
    zt = physics.z_at_age(t * 1e3)
    sfr = lit.md14_sfrd(zt)
    ax2.plot(t, sfr / sfr.max(), color=MUTED, lw=1.6)
    imax = np.argmax(sfr)
    ax2.text(t[imax] + 0.3, 1.08, f"star formation (Madau & Dickinson 2014 fit)\npeaks at z = {zt[imax]:.1f}, "
             f"{t[imax]:.1f} Gyr", color=INK_SECONDARY, fontsize=8)
    for cut, (col, lab) in styles.items():
        d = era[cut]
        c = d["complete"].to_numpy(bool)
        y = (d["density"] / d["density"][c].max()).to_numpy()
        ax2.plot(physics.age_myr(d["z_mean"][c]) / 1e3, y[c], "o-", color=col, ms=4, lw=1.0, label=lab)
    ax2.set_yscale("log")
    ax2.set_ylim(0.05, 2)
    ax2.set_xlim(0, 13.8)
    ax2.set_xlabel("age of the universe (billion years)")
    ax2.set_ylabel("density relative to its own peak")
    ax2.legend(loc="lower right", fontsize=7.5, title="luminous quasars (complete bins)", title_fontsize=7.5)
    age_top_axis(ax2, zs=(0, 0.5, 1, 2, 3, 5))
    ax2.set_title("Quasars peaked before star formation", loc="left", pad=28)
    caption(fig, "Left: the comoving number density of luminous quasars (M_i(z=2) brighter than "
            "-27.6 and -28.5) from the SDSS DR7 uniform sample, summed over the luminosity function in "
            "D3's cosmology. Filled points: every quasar that luminous would have been inside the "
            "survey's flux limits; open: partly outside, so lower limits. The orange curve is a parabola "
            "in log n vs log(1+z) fitted to the complete points; its maximum, with a bootstrap over "
            "quasars, gives the peak; the range in brackets is how far it moves with other reasonable "
            "choices (selection-function floor, fitting form, redshift range). Squares: the same "
            "quantity from Richards et al.'s published DR3 bins, moved into D3's cosmology (a 0.05-0.09 "
            "mag shift in luminosity and a 5-10% change in volume). Right: the shapes, each scaled to its "
            "own maximum. The most luminous quasars peak slightly earlier (z ≈ 2.8 for M_i < -28.2 "
            "against 2.6), and all of them about a billion years before cosmic star formation "
            "(z ≈ 1.9).", y=0.02)
    return save(fig, FIGURES / "06_quasar_era")


# ---------------------------------------------------------------------------
# 7. Black hole masses
# ---------------------------------------------------------------------------

def fig_masses(res):
    m = res["masses"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.8), gridspec_kw=dict(width_ratios=[1.4, 1.0], wspace=0.22))
    ok = np.isfinite(m["m"]) & np.isfinite(m["z"])
    density(ax, m["z"][ok], m["m"][ok], np.linspace(0, 5.2, 209), np.linspace(6.0, 11.0, 201))
    pz = m["per_z"]
    zc = 0.5 * (pz["z_lo"] + pz["z_hi"])
    for line, (lo, hi) in (("Hbeta", (0, 0.75)), ("MgII", (0.6, 2.05)), ("CIV", (1.95, 5.1))):
        s = (zc > lo) & (zc < hi)
        ax.plot(zc[s], pz["m_p50"][s], color=LINE_COLOURS[line], lw=2.0)
    ax.text(0.12, 8.5, "median, by Hβ", color=LINE_COLOURS["Hbeta"], fontsize=8)
    ax.text(1.15, 8.95, "Mg II", color=LINE_COLOURS["MgII"], fontsize=8)
    ax.text(3.05, 9.15, "C IV", color=LINE_COLOURS["CIV"], fontsize=8)
    lim = m["limits"]
    zl = lim["z"]
    ax.plot(zl, lim["m_min_eboss"], color=INK, lw=1.1)
    ax.plot(zl, lim["m_min_eboss"] + 1.0, color=INK_SECONDARY, lw=0.7)
    kl = np.where(np.isfinite(lim["m_min_eboss"]))[0][-1]
    ax.text(zl[kl] + 0.08, lim["m_min_eboss"][kl], "at L_Edd", fontsize=7.5, color=INK, va="center")
    ax.text(zl[kl] + 0.08, lim["m_min_eboss"][kl] + 1.0, "at 0.1 L_Edd", fontsize=7.5,
            color=INK_SECONDARY, va="center")
    annotate(ax, "lightest black hole BOSS/eBOSS could see,\nfeeding at the stated Eddington ratio",
             (1.2, np.interp(1.2, zl, lim["m_min_eboss"])), (1.6, 6.35), fontsize=7.5, color=INK)
    hz = res["growth"]["hz"]
    ax.scatter(hz["z"], hz["log_mbh"], s=12, color=C_PUBLISHED, edgecolor=SURFACE, lw=0.4, zorder=5)
    ax.set_xlim(0, 7.8)
    ax.set_ylim(6.0, 11.0)
    ax.text(6.0, 10.4, "z > 5.7 quasars\n(Mg II, same recipe)", color=C_PUBLISHED, fontsize=7.5, ha="center")
    step = m["z2_step"]
    s2 = step[np.isclose(step["z_mid"], 2.05)].iloc[0]
    s1 = step[np.isclose(step["z_mid"], 1.95)].iloc[0]
    annotate(ax, f"C IV takes over at z = 2: the median\ndrops {s2['fiducial'] - s1['fiducial']:+.2f} dex "
             "(an artefact, figure 9)", (2.05, s2["fiducial"]), (2.6, 7.0), fontsize=7.5)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("log black hole mass  (M☉)")
    ax.set_title(f"Black hole masses of {m['n']:,} quasars", loc="left")

    density(ax2, m["lbol"][ok], m["m"][ok], np.linspace(44.0, 48.4, 177), np.linspace(6.0, 11.0, 201))
    led = np.log10(physics.eddington_per_msun())
    xx = np.linspace(44, 48.4, 10)
    ax2.set_xlim(44.0, 48.4)
    ax2.set_ylim(6.0, 11.0)
    for lam, lab in ((1.0, "L = L_Edd"), (0.1, "0.1 L_Edd"), (0.01, "0.01 L_Edd")):
        ax2.plot(xx, xx - led - np.log10(lam), color=INK if lam == 1 else INK_SECONDARY, lw=1.0 if lam == 1 else 0.7)
        xl = 47.75 if lam != 0.01 else 45.6
        yl = xl - led - np.log10(lam)
        ang = ax2.transData.transform_angles(np.array([45.0]), np.array([[xl, yl]]))[0]
        ax2.text(xl, yl + 0.1, lab, rotation=ang, rotation_mode="anchor", fontsize=7.5,
                 color=INK if lam == 1 else INK_SECONDARY, ha="center", va="bottom")
    ax2.set_xlim(44.0, 48.4)
    ax2.set_ylim(6.0, 11.0)
    ax2.set_xlabel("log bolometric luminosity  (erg/s)")
    ax2.set_ylabel("log black hole mass  (M☉)")
    ax2.set_title("Feeding rates: few quasars beat the Eddington limit", loc="left")
    ax2.text(0.97, 0.05, f"{100 * m['frac_super_edd_all']:.1f}% above L_Edd\nmedian L/L_Edd = "
             f"{10 ** m['median_edd_all']:.2f}", transform=ax2.transAxes, ha="right", fontsize=8,
             color=INK_SECONDARY)
    caption(fig, "Single-epoch virial masses from one spectrum each: line width (gas speed) and continuum "
            "luminosity (via the radius-luminosity relation of reverberation mapping, project I2). Recipes: "
            "Vestergaard & Peterson 2006 for Hβ and C IV, Shen et al. 2011 for Mg II, recomputed from Wu & "
            "Shen's measurements (identical to their catalogue). Left: the lower edge of the cloud is set by "
            "the survey: black lines show the mass of a quasar at the flux limit radiating at the Eddington "
            "limit (thin: at 10% of it). Most quasars sit well above the line, so most are found well below "
            "Eddington. Each mass is uncertain by ~0.4 dex. Right: the mass-luminosity plane with lines of "
            "constant Eddington ratio (L_Edd = 1.26 × 10³⁸ erg/s per solar mass).", y=0.02)
    return save(fig, FIGURES / "07_black_hole_masses")


# ---------------------------------------------------------------------------
# 8. Eddington ratios
# ---------------------------------------------------------------------------

def fig_eddington(res):
    m = res["masses"]
    kol = m["kollmeier"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.0), gridspec_kw=dict(wspace=0.22))
    sel = (m["z"] > 0.3) & (m["z"] < 4) & (m["lbol"] > 45) & (m["lbol"] < 47) & np.isfinite(m["edd"])
    e = m["edd"][sel]
    bins = np.linspace(-3.0, 1.0, 81)
    h, _ = np.histogram(e, bins)
    ax.stairs(h / h.sum() / np.diff(bins), bins, color=C_OURS, lw=1.5, fill=False)
    x = np.linspace(-3, 1, 400)
    mu, sg = np.log10(lit.KOLLMEIER06["peak"]), lit.KOLLMEIER06["sigma_dex"]
    ax.plot(x, np.exp(-0.5 * ((x - mu) / sg) ** 2) / (sg * np.sqrt(2 * np.pi)), color=C_PUBLISHED, lw=1.4)
    ax.axvline(0, color=INK, lw=0.8)
    ax.text(0.04, 1.3, "Eddington\nlimit", fontsize=7.5, color=INK, va="top")
    ax.text(-2.9, 1.15, f"DR16Q, 0.3 < z < 4,\n10⁴⁵ < L_bol < 10⁴⁷ erg/s\n({kol['n']:,} quasars):\n"
            f"median {10 ** kol['median']:.2f}, spread {kol['robust_std']:.2f} dex", color=C_OURS, fontsize=8, va="top")
    ax.text(0.06, 0.95, "Kollmeier et al. 2006\n(407 AGN): peak 1/4,\nspread 0.3 dex", color=C_PUBLISHED,
            fontsize=8, va="top", ha="left")
    ax.set_xlim(-3, 1)
    ax.set_xlabel("log Eddington ratio  L_bol / L_Edd")
    ax.set_ylabel("fraction per dex")
    ax.set_title("How hard quasars feed", loc="left")

    pz = m["per_z"]
    zc = 0.5 * (pz["z_lo"] + pz["z_hi"])
    for line, (lo, hi) in (("Hbeta", (0, 0.7)), ("MgII", (0.7, 2.0)), ("CIV", (2.0, 5.1))):
        s = (zc > lo) & (zc < hi)
        ax2.fill_between(zc[s], (pz["edd_p50"] - pz["edd_std"])[s], (pz["edd_p50"] + pz["edd_std"])[s],
                         color=LINE_COLOURS[line], alpha=0.18, lw=0)
        ax2.plot(zc[s], pz["edd_p50"][s], "o-", color=LINE_COLOURS[line], ms=4, lw=1.4,
                 label=LINE_LABELS[line])
    ax2.axhline(0, color=INK, lw=0.8)
    ax2.set_xlim(0, 5)
    ax2.set_ylim(-2.2, 0.6)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("log L_bol / L_Edd  (median ± spread)")
    ax2.legend(loc="lower right", title="mass from", title_fontsize=8)
    s1 = pz[(zc > 1.8) & (zc < 2.0)].iloc[0]
    s2 = pz[(zc > 2.0) & (zc < 2.2)].iloc[0]
    annotate(ax2, f"jump of {s2['edd_p50'] - s1['edd_p50']:+.2f} dex where C IV takes over:\n"
             "the C IV masses are too low here (figure 9)", (2.12, s2["edd_p50"]), (2.4, 0.35), fontsize=7.5)
    ax2.set_title("Eddington ratio across cosmic time", loc="left")
    caption(fig, "Left: Eddington ratios of DR16Q quasars in the luminosity and redshift range of "
            "Kollmeier et al. (2006). The distribution is a lognormal pile-up below the Eddington limit, "
            f"with {100 * kol['frac_above_1']:.1f}% formally above it (a lognormal of the observed width "
            f"would put {100 * kol['gauss_frac_above_1']:.1f}% there; without the ~0.3 dex mass errors, "
            f"{100 * kol['intrinsic_frac_above_1']:.1f}%). Its centre sits a factor ~2 below "
            "Kollmeier et al.'s (different recipes and bolometric corrections; a known systematic of this "
            f"size), and its width ({kol['robust_std']:.2f} dex) includes ~0.3 dex of mass error. Right: "
            "the median rises with redshift partly because flux-limited surveys see only the most "
            "luminous, fastest-feeding quasars far away; the step at z = 2 is the change of emission line.",
            y=0.02)
    return save(fig, FIGURES / "08_eddington_ratios")


# ---------------------------------------------------------------------------
# 9. Three rulers
# ---------------------------------------------------------------------------

def fig_rulers(res):
    m = res["masses"]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2), gridspec_kw=dict(hspace=0.32, wspace=0.24))
    ax = axes[0, 0]
    hb = m["hb_mg"]
    density(ax, hb["b"], hb["a"], np.linspace(6.5, 10.5, 161), np.linspace(6.5, 10.5, 161))
    ax.plot([6.5, 10.5], [6.5, 10.5], color=INK, lw=0.8)
    ax.set_xlim(6.5, 10.5)
    ax.set_ylim(6.5, 10.5)
    ax.set_xlabel("log M from Hβ  (Vestergaard & Peterson 2006)")
    ax.set_ylabel("log M from Mg II  (Shen et al. 2011)")
    ax.text(0.03, 0.97, f"{hb['n']:,} quasars at 0.4 < z < 0.9\nmedian offset {hb['median']:+.2f} dex\n"
            f"scatter {hb['robust_std']:.2f} dex", transform=ax.transAxes, va="top", fontsize=8,
            color=INK_SECONDARY)
    ax.set_title("Two lines, one black hole", loc="left")

    ax = axes[0, 1]
    rec = m["recipes"]
    mh = rec["mg_hb"]
    ax.plot(mh["l_mid"], mh["s11_minus_hb"], "o-", color=C_OURS, label="Mg II (Shen 2011, L^0.62) − Hβ")
    ax.plot(mh["l_mid"], mh["vo09_minus_hb"], "s-", color=C_PUBLISHED,
            label="Mg II (Vestergaard & Osmer 2009, L^0.5) − Hβ")
    c4 = rec["c4_mg"]
    ax.plot(c4["l_mid"], c4["civ_minus_s11"], "^-", color=C_MODEL, label="C IV − Mg II (vs L_1350)")
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_ylim(-0.52, 0.12)
    ax.set_xlabel("log λL_λ at 3000 Å (Mg II − Hβ) or 1350 Å (C IV − Mg II)   (erg/s)")
    ax.set_ylabel("median mass offset (dex)")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("The recipes agree where they were calibrated", loc="left")
    ax.text(0.03, 0.97, "luminous quasars (right):\nlines agree.\nfainter ones (left): up to\n0.3 dex apart",
            transform=ax.transAxes, va="top", fontsize=8, color=INK_SECONDARY)

    ax = axes[1, 0]
    bs = m["blueshift"]
    b = bs["bins"]
    ax.fill_between(b["bs_mid"], b["raw_q16"], b["raw_q84"], color=ORANGE, alpha=0.15, lw=0)
    ax.plot(b["bs_mid"], b["raw"], "o-", color=C_MODEL, ms=4, label="C IV as published recipe (VP06)")
    ax.plot(b["bs_mid"], b["c17"], "s-", color=C_PUBLISHED, ms=4,
            label="corrected with Coatman et al. 2017")
    ax.plot(b["bs_mid"], b["ours"], "o-", color=C_OURS, ms=4, label="our recalibration (cross-validated)")
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xlim(-1000, 4000)
    ax.set_ylim(-1.1, 1.0)
    ax.set_xlabel("C IV blueshift relative to Mg II  (km/s; outflowing gas)")
    ax.set_ylabel("log M(C IV) − log M(Mg II)  (median per bin)")
    ax.legend(loc="lower right", fontsize=7.5)
    cal = bs["calibration"]
    ax.text(0.03, 0.97, f"{bs['n']:,} quasars with both lines\nslope: {bs['slope_raw']:+.2f} → "
            f"{bs['slope_ours']:+.2f} dex per 1000 km/s\nscatter: {bs['std_raw']:.2f} → {bs['std_ours']:.2f} dex\n"
            f"FWHM(C IV)/FWHM(Mg II) = {cal['slope']:.2f}·Δv/1000 + {cal['intercept']:.2f}",
            transform=ax.transAxes, va="top", fontsize=8, color=INK_SECONDARY)
    ax.set_title("Winds inflate C IV widths", loc="left")

    ax = axes[1, 1]
    st = m["z2_step"]
    st = st[st["n"] > 100]
    ax.plot(st["z_mid"], st["fiducial"], "-", color=MUTED, lw=4, alpha=0.6, label="catalogue's fiducial mass")
    ax.plot(st["z_mid"], st["mgii"], "o-", color=LINE_COLOURS["MgII"], ms=4, label="Mg II")
    ax.plot(st["z_mid"], st["civ"], "^-", color=LINE_COLOURS["CIV"], ms=4, label="C IV, published recipe")
    ax.plot(st["z_mid"], st["civ_recal"], "^--", color=LINE_COLOURS["CIV"], mfc=SURFACE, ms=5,
            label="C IV, our recalibration")
    ax.axvline(2.0, color=AXIS, lw=0.8)
    ax.text(2.015, 9.0, "switch\nMg II → C IV", fontsize=7.5, color=INK_SECONDARY, va="top")
    ax.set_xlabel("redshift z  (same quasars for every curve)")
    ax.set_ylabel("median log M  (M☉)")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("An artificial step in cosmic history", loc="left")
    caption(fig, "Masses from different emission lines in the same spectrum. Top left: Mg II vs Hβ at "
            "0.4 < z < 0.9 (BOSS spectra). Top right: the offsets depend on luminosity, because the "
            "recipes scale differently with it and were matched on luminous SDSS-I/II quasars; the "
            "Vestergaard & Osmer Mg II recipe, which scales like Hβ's, stays within 0.12 dex. Bottom "
            "left: C IV is broadened by outflowing gas, and its mass relative to Mg II climbs with the "
            "blueshift that traces that wind (at fixed luminosity too). Coatman et al.'s correction, "
            "calibrated on ten-times-more-luminous quasars against Hα, over-corrects these; our own "
            "relation, fitted on one random half and tested on the other, flattens the trend. Bottom "
            "right: for quasars with both lines (1.6 < z < 2.4), Mg II and the recalibrated C IV run "
            "smoothly through z = 2, while the published C IV recipe sits ~0.25 dex lower.", y=0.02)
    return save(fig, FIGURES / "09_three_rulers")


# ---------------------------------------------------------------------------
# 10. The growth problem
# ---------------------------------------------------------------------------

def fig_growth(res):
    g = res["growth"]
    hz, seeds = g["hz"], g["seeds"]
    fig = plt.figure(figsize=(13.0, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 0.62], wspace=0.28)
    ax = fig.add_subplot(gs[0])
    band_cols = ["#eeede8", "#e1e0d9", "#d3d1c8"]
    for (name, (lo, hi)), col in zip(SEED_CLASSES.items(), band_cols):
        ax.axhspan(np.log10(lo), np.log10(hi), color=col, lw=0, zorder=0)
        ax.text(0.103, np.log10(np.sqrt(lo * hi)), name + " seeds", fontsize=7.5, color=INK_SECONDARY,
                va="center")
    for name, (t, lm) in g["tracks"].items():
        if name == "J0313-1806":
            continue
        ax.plot(t / 1e3, lm, color="#b4b2aa", lw=0.5, zorder=2)
    t, lm = g["tracks"]["J0313-1806"]
    ax.plot(t / 1e3, lm, color=C_OURS, lw=2.0, zorder=4)
    ax.scatter(hz["age_myr"] / 1e3, hz["log_mbh"], s=12, color=INK, zorder=5, edgecolor=SURFACE, lw=0.4)
    j = hz[hz["name"] == "J0313-1806"].iloc[0]
    ax.scatter([j["age_myr"] / 1e3], [j["log_mbh"]], s=40, color=C_OURS, zorder=6, edgecolor=SURFACE, lw=0.8)
    annotate(ax, f"J0313-1806, z = {j['z']:.2f}\n{10 ** (j['log_mbh'] - 9):.1f} billion M☉, "
             f"{j['age_myr']:.0f} Myr after the Big Bang", (j["age_myr"] / 1e3, j["log_mbh"]),
             (0.25, 10.25), fontsize=8, color=C_OURS)
    jw = g["jwst"]
    for _, r in jw.iterrows():
        tj, lj = growth_track_simple(r["log_mbh"], r["z"])
        ax.plot(tj / 1e3, lj, color=C_PUBLISHED, lw=0.8, alpha=0.6, zorder=3)
        ax.scatter([r["age_myr"] / 1e3], [r["log_mbh"]], marker="D", s=26, color=C_PUBLISHED, zorder=6,
                   edgecolor=SURFACE, lw=0.6)
        ax.annotate(r["name"], (r["age_myr"] / 1e3, r["log_mbh"]), xytext=(5, -10), textcoords="offset points",
                    fontsize=7, color=C_PUBLISHED)
    ax.set_xscale("log")
    ax.set_xlim(0.098, 1.1)
    ax.set_ylim(1.0, 10.8)
    ticks = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{x:g}" for x in ticks])
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlabel("age of the universe (billion years)")
    ax.set_ylabel("log black hole mass  (M☉)")
    top = ax.twiny()
    top.set_xscale("log")
    top.set_xlim(ax.get_xlim())
    zt = [30, 20, 15, 10, 8, 7, 6]
    top.set_xticks(physics.age_myr(np.array(zt, float)) / 1e3)
    top.set_xticklabels([str(z) for z in zt])
    top.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    top.set_xlabel("redshift z", color=INK_SECONDARY)
    top.grid(False)
    ax.set_title("Run the clock backwards", loc="left", pad=30)

    ax2 = fig.add_subplot(gs[1])
    for (name, (lo, hi)), col in zip(SEED_CLASSES.items(), band_cols):
        ax2.axhspan(np.log10(lo), np.log10(hi), color=col, lw=0, zorder=0)
    ax2.scatter(seeds["z"], seeds["log_seed_20"], s=14, facecolor="none", edgecolor=INK_SECONDARY, lw=0.7,
                label="growth started at z = 20", zorder=3)
    ax2.scatter(seeds["z"], seeds["log_seed_30"], s=14, color=C_OURS, edgecolor=SURFACE, lw=0.4,
                label="growth started at z = 30", zorder=4)
    for _, r in g["jwst"].iterrows():
        ax2.scatter([r["z"]], [r["log_seed_edd"]], marker="D", s=26, color=C_PUBLISHED, edgecolor=SURFACE,
                    lw=0.6, zorder=5)
    ax2.axhline(2.0, color=INK, lw=0.8)
    ax2.text(10.4, 2.08, "100 M☉", fontsize=7.5, color=INK, ha="right")
    ax2.set_xlim(5.6, 10.9)
    ax2.set_ylim(0, 6.6)
    ax2.set_xlabel("redshift of the quasar")
    ax2.set_ylabel("log minimum seed mass  (M☉)")
    ax2.legend(loc="upper left", fontsize=7.5)
    f = g["frac_seed_above_100"]
    ax2.text(0.97, 0.03, f"needing a seed > 100 M☉:\n{100 * f[30]:.0f}% (from z = 30)\n{100 * f[20]:.0f}% (from z = 20)",
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=8, color=INK_SECONDARY)
    ax2.set_title("The seed each quasar needs", loc="left", pad=30)

    ax3 = fig.add_subplot(gs[2])
    eff = g["efficiency"]
    ax3.plot([e["eps"] for e in eff], [e["log_seed"] for e in eff], "o-", color=C_OURS, ms=4)
    for e in eff:
        ax3.annotate(f"{e['efold_myr']:.0f} Myr", (e["eps"], e["log_seed"]), xytext=(6, -3),
                     textcoords="offset points", fontsize=7, color=INK_SECONDARY)
    ax3.set_xlabel("radiative efficiency ε")
    ax3.set_ylabel("log seed for J0313-1806, from z = 30")
    ax3.set_xlim(0, 0.5)
    ax3.set_title("Spin makes it worse", loc="left", pad=30)
    ax3.text(0.03, 0.97, "labels: mass\ne-folding time", transform=ax3.transAxes, va="top", fontsize=7,
             color=INK_SECONDARY)
    caption(fig, f"Left: every quasar beyond z = 5.7 with a Mg II width (Farina et al. 2022 compilation, "
            f"masses recomputed with one recipe), traced back in time assuming it grew without pause at "
            f"the Eddington limit with 10% radiative efficiency, an e-folding every {g['efold_myr']:.0f} Myr. "
            f"Shaded: the masses expected for black holes left by the first stars, by runaway collisions in "
            f"dense star clusters, and by gas clouds collapsing directly. Diamonds: four black holes found "
            f"by JWST (published masses; methods differ, see README). Middle: the seed each quasar needs; "
            f"J0313-1806 needs at least ~{round(10 ** seeds.loc[seeds['name'] == 'J0313-1806', 'log_seed_30'].iloc[0], -3):,.0f} "
            f"M☉ even if it started at z = 30 (its mass here, 3.2 billion M☉, uses our recipe with Yang et al. "
            f"2021's line width; Wang et al. 2021 give 1.6 billion, which halves the seed). Right: a spinning "
            f"black hole converts more of its inflow into light (higher ε) and so grows more slowly; "
            f"ε = 0.057-0.42 spans no spin to maximal spin.",
            y=0.02)
    return save(fig, FIGURES / "10_growth_problem")


def growth_track_simple(log_m, z, z_start=30.0):
    import highz
    return highz.growth_track(log_m, z, z_start)


FIGURE_FUNCS = {1: fig_census, 2: fig_luminosity, 3: fig_composite, 4: fig_colours, 5: fig_qlf,
                6: fig_era, 7: fig_masses, 8: fig_eddington, 9: fig_rulers, 10: fig_growth}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    FIGURES.mkdir(parents=True, exist_ok=True)
    res = load_results()
    which = [int(a) for a in argv] if argv else sorted(FIGURE_FUNCS)
    for k in which:
        out = FIGURE_FUNCS[k](res)
        print(f"figure {k}: {', '.join(str(p) for p in out)}")


if __name__ == "__main__":
    main()
