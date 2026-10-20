"""Every figure of the project, from the results and arrays of ``analysis.run_all``.

Maps use supergalactic coordinates: SGZ points to the supergalactic pole and
the SGX-SGY plane is the sheet in which the nearby galaxies lie (de Vaucouleurs
1953). Distances are in h^-1 Mpc (1 h^-1 Mpc = 100 km/s of expansion).
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm, Normalize, TwoSlopeNorm

import catalogs
import config
import cosmology as cg
import density as dn
import plotstyle as ps

ps.use_style()
FIG = config.FIGURE_DIR

#: Landmarks drawn on supergalactic-plane maps (the ones within |SGZ| < 30).
LABEL_OFFSETS = {"Virgo": (6, 6), "Fornax": (6, -12), "Centaurus": (2, -13), "Hydra": (5, 4), "Norma": (-40, -14),
                 "Perseus": (6, 6), "Coma": (6, 6), "Ophiuchus": (6, 6), "Hercules": (6, 6), "Shapley": (6, 6)}


# --------------------------------------------------------------------------
# Helpers: sampling the 3D grids on the supergalactic plane
# --------------------------------------------------------------------------

def sg_plane_grid(half: float = 160.0, step: float = 1.0, sgz: float = 0.0):
    x = np.arange(-half, half + step / 2, step)
    X, Y = np.meshgrid(x, x, indexing="xy")
    pts_sg = np.stack([X.ravel(), Y.ravel(), np.full(X.size, sgz)], axis=1)
    pts_gal = catalogs.sg_to_gal(pts_sg)
    return x, X, Y, pts_gal


def sample_scalar_sg(field, pts_gal, shape):
    inside = np.all(np.abs(pts_gal) < config.RECON_HALF - 2, axis=1)
    out = np.full(len(pts_gal), np.nan)
    out[inside] = dn.sample_field(field, pts_gal[inside])
    return out.reshape(shape)


def sample_labels_sg(labels, pts_gal, shape):
    c, cell = dn.grid_axes(labels.shape[0])
    idx = np.clip(np.round((pts_gal + config.RECON_HALF) / cell - 0.5).astype(int), 0, labels.shape[0] - 1)
    return labels[idx[:, 0], idx[:, 1], idx[:, 2]].reshape(shape)


def sample_velocity_sg(v, pts_gal, shape):
    inside = np.all(np.abs(pts_gal) < config.RECON_HALF - 2, axis=1)
    out = np.full((len(pts_gal), 3), np.nan)
    out[inside] = catalogs.gal_to_sg(dn.sample_field(v, pts_gal[inside]))
    return out[:, 0].reshape(shape), out[:, 1].reshape(shape), out[:, 2].reshape(shape)


def zoa_wedge_sg(ax, rmax, color=ps.MUTED):
    """Mark where the Milky Way hides the sky (|b| < 5 deg) on a supergalactic-plane map: the
    Galactic plane cuts the supergalactic plane along a line, and the hidden strip is a thin
    double wedge around it, drawn as its two edges."""
    ang = np.linspace(0, 2 * np.pi, 7201)
    pts = np.stack([np.cos(ang), np.sin(ang), np.zeros_like(ang)], axis=1)
    _, b = catalogs.lonlat(catalogs.sg_to_gal(pts))
    edge = np.flatnonzero(np.diff(np.sign(np.abs(b) - config.ZOA_B)) != 0)
    for i in edge:
        a = ang[i]
        ax.plot([0, rmax * np.cos(a)], [0, rmax * np.sin(a)], color=color, lw=0.7, ls=(0, (1.5, 2.5)), zorder=5)
    return ang[edge]


def basin_indicator_sg(labels, label, pts_gal, shape):
    """Membership of one basin, interpolated trilinearly so its outline is smooth."""
    ind = (labels == label).astype(np.float32)
    return sample_scalar_sg(ind, pts_gal, shape)


def draw_landmarks(ax, lm, *, zmax=30.0, rmax=160.0, color=ps.INK, fontsize=7.5):
    for _, row in lm.iterrows():
        if abs(row["SGZ"]) > zmax or row["r"] > rmax:
            continue
        ax.plot(row["SGX"], row["SGY"], marker="o", ms=3.5, mfc="none", mec=color, mew=0.9, zorder=6)
        dx, dy = LABEL_OFFSETS.get(row["name"], (5, 5))
        ax.annotate(row["name"], (row["SGX"], row["SGY"]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=fontsize, color=color, zorder=7)


def milky_way(ax, color=ps.INK):
    ax.plot(0, 0, marker="*", ms=9, color=color, mec=ps.SURFACE, mew=0.6, zorder=8)
    ax.annotate("Milky Way", (0, 0), xytext=(5, -11), textcoords="offset points", fontsize=7.5, color=color, zorder=8)


def sg_axes(ax, half):
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    ax.set_xlabel("SGX  (h$^{-1}$ Mpc)")
    ax.set_ylabel("SGY  (h$^{-1}$ Mpc)")
    ax.grid(False)


# --------------------------------------------------------------------------
# 00. Laniakea
# --------------------------------------------------------------------------

def figure_laniakea(res, art) -> list:
    rec = art["rec"]
    beta = art["vc"]["beta_grid"]
    delta = dn.field_for(rec, beta, smooth=5.0)
    _, cell = dn.grid_axes()
    bkey = f"{config.BASIN_SMOOTH_MAIN:g}"
    bas = art["basins"][bkey]
    rep = bas["report"]
    extra = np.sqrt(config.BASIN_SMOOTH_MAIN**2 - config.RECON_SMOOTH**2)
    dsm = dn.gaussian_smooth(dn.field_for(rec, beta), extra, cell)
    v = dn.velocity_field(dsm, beta, cell)
    half = 158.0
    x, X, Y, pts = sg_plane_grid(half, 1.0)
    dens = sample_scalar_sg(delta, pts, X.shape)
    xs, Xs, Ys, pts_s = sg_plane_grid(half, 2.5)
    U, V, _ = sample_velocity_sg(v, pts_s, Xs.shape)
    R = np.hypot(X, Y)
    fig = plt.figure(figsize=(13.2, 8.8))
    ax = fig.add_axes([0.05, 0.08, 0.62, 0.84])
    tax = fig.add_axes([0.70, 0.08, 0.29, 0.84])
    tax.axis("off")
    one = np.clip(1 + dens, 0.15, 12)
    one[R > config.RECON_RMAX] = np.nan
    ax.pcolormesh(x, x, one, cmap=ps.density_cmap(), norm=LogNorm(0.15, 12), shading="nearest", zorder=0,
                  rasterized=True)
    inside = R < config.RECON_RMAX - 2
    ours = rep["ours_label"]
    # Outlines of the other large basins, then ours on top.
    for b in rep["largest"][:10]:
        if b["is_ours"]:
            continue
        ind = basin_indicator_sg(bas["labels"], b["label"], pts, X.shape)
        ax.contour(x, x, np.where(inside, ind, np.nan), levels=[0.5], colors=[ps.INK_SECONDARY], linewidths=0.9,
                   alpha=0.8, zorder=3)
    ind = basin_indicator_sg(bas["labels"], ours, pts, X.shape)
    ind = np.where(inside, ind, np.nan)
    ax.contourf(x, x, ind, levels=[0.5, 1.5], colors=[ps.C_BASIN], alpha=0.13, zorder=3)
    ax.contour(x, x, ind, levels=[0.5], colors=[ps.C_BASIN], linewidths=2.4, zorder=4)
    Rs = np.hypot(Xs, Ys)
    U = np.where(Rs < config.RECON_RMAX - 4, U, np.nan)
    V = np.where(Rs < config.RECON_RMAX - 4, V, np.nan)
    ax.streamplot(xs, xs, U, V, density=2.1, color=ps.C_STREAM, linewidth=0.55, arrowsize=0.7, zorder=5)
    zoa_wedge_sg(ax, config.RECON_RMAX)
    ax.add_patch(plt.Circle((0, 0), config.RECON_RMAX, fill=False, color=ps.MUTED, lw=0.8, ls=(0, (2, 3)), zorder=5))
    att = rep["ours_attractor_sg"]
    ax.plot(att[0], att[1], marker="X", ms=11, color=ps.C_BASIN, mec=ps.INK, mew=0.8, zorder=9)
    ax.annotate("attractor", (att[0], att[1]), xytext=(-14, 12), textcoords="offset points", fontsize=8.5,
                color=ps.INK, fontweight="bold", zorder=9, ha="right")
    draw_landmarks(ax, art["landmarks"])
    milky_way(ax)
    sg_axes(ax, half)
    ax.set_title("Where the Milky Way is falling")
    # The side panel: what the map shows and how it compares with Laniakea.
    b = res["basins"]
    rows = [(k, b[k]) for k in ("6", "10", "15")]
    lines = [("Our home basin of attraction", True),
             ("The orange region drains to the same attractor as the", False),
             ("Milky Way: every point in it, followed downhill along", False),
             ("the gravitational flow, ends at the orange cross, near", False),
             ("the Centaurus cluster and the Great Attractor.", False), ("", False),
             ("Mass: the 2MRS survey's 43,285 galaxies, weighted", False),
             ("for what a flux-limited survey misses (grey, darker =", False),
             ("denser). Flow: linear theory, with the strength of", False),
             (f"gravity beta* = {res['velocity_comparison']['beta']:.2f} fitted to Cosmicflows-4. Lines", False),
             ("follow the flow of the mass inside 160 h-1 Mpc", False),
             ("(dashed circle); grey outlines are other basins.", False), ("", False),
             ("Size depends on how much the field is smoothed", True)]
    y = 0.98
    for text, bold in lines:
        tax.text(0.0, y, text, fontsize=9 if bold else 8.4, fontweight="bold" if bold else "normal",
                 color=ps.INK if bold else ps.INK_SECONDARY, va="top", transform=tax.transAxes)
        y -= 0.034 if text else 0.02
    y -= 0.01
    hdr = ["smoothing", "diameter", "mass"]
    xs_ = [0.0, 0.36, 0.72]
    for xx, hh in zip(xs_, hdr):
        tax.text(xx, y, hh, fontsize=8.2, color=ps.MUTED, va="top", transform=tax.transAxes)
    y -= 0.034
    for k, rep_k in rows:
        d_kms = rep_k["ours_equiv_diameter"] * 100
        mark = " (map)" if k == bkey else ""
        vals = [f"{k} h-1 Mpc{mark}", f"{d_kms:,.0f} km/s", f"{rep_k['ours_mass_msun'] / 1e17:.1f} x 1e17 Msun"]
        for xx, vv in zip(xs_, vals):
            tax.text(xx, y, vv, fontsize=8.4, color=ps.INK, va="top", transform=tax.transAxes)
        y -= 0.032
    for xx, vv in zip(xs_, ["Tully+2014", "12,000 km/s", "1 x 1e17 Msun"]):
        tax.text(xx, y, vv, fontsize=8.4, color=ps.C_PUBLISHED, va="top", transform=tax.transAxes)
    y -= 0.05
    more = [("Diameter: of a sphere of equal volume; 100 km/s =", False),
            ("1 h-1 Mpc, i.e. 1.33 Mpc for H0 = 75. Mass: mean", False),
            ("matter density times volume (h = 0.75).", False), ("", False),
            ("At 6 h-1 Mpc our basin matches Laniakea's size and", False),
            ("mass; at 15 h-1 Mpc it swallows Coma's. Adding the", False),
            ("pull of mass beyond 2MRS tips us into the Shapley", False),
            ("basin (figure 09), as Valade et al. (2024) found", False),
            ("plausible with Cosmicflows-4.", False)]
    for text, bold in more:
        tax.text(0.0, y, text, fontsize=8.4, color=ps.INK_SECONDARY, va="top", transform=tax.transAxes)
        y -= 0.034 if text else 0.02
    ps.caption(fig, "2MRS (Huchra et al. 2012) reconstructed in real space; dotted lines bound the strip the Milky Way's "
                    "disc hides (|b| < 5 deg), filled by cloning the strips beside it. Supergalactic coordinates, SGZ = 0.")
    return ps.save(fig, FIG / "00_laniakea")


# --------------------------------------------------------------------------
# 01. The Hubble diagram
# --------------------------------------------------------------------------

def figure_hubble(res, art) -> list:
    h = art["hubble_points"]
    h = h[(h["V3k"] > 0) & np.isfinite(h["logH"])]
    hb = res["hubble"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), gridspec_kw={"width_ratios": [1.15, 1.15, 0.9]})
    ax = axes[0]
    sn = h["method"] == "SN Ia"
    ax.scatter(h.loc[~sn, "d_mpc"], h.loc[~sn, "fcz"], s=1.2, color=ps.C_CONTEXT, rasterized=True, lw=0,
               label=f"{(~sn).sum():,} groups (TF, FP, SBF, SN II, calibrators)")
    ax.scatter(h.loc[sn, "d_mpc"], h.loc[sn, "fcz"], s=5, color=ps.BLUE, rasterized=True, lw=0,
               label=f"{sn.sum():,} groups with a type Ia supernova")
    d = np.geomspace(0.5, 800, 200)
    H0c = hb["corrected"]["baseline"]["H0"]
    ax.plot(d, H0c * d, color=ps.C_OURS, lw=1.6, label=f"H0 = {H0c:.1f} (CF4, our fit)")
    ax.plot(d, config.PLANCK18.H0 * d, color=ps.C_PUBLISHED, lw=1.2, ls=(0, (4, 2)),
            label=f"H0 = {config.PLANCK18.H0:.1f} (Planck CMB)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 800)
    ax.set_ylim(30, 4e4)
    ax.set_xlabel("distance (Mpc, Cepheid / TRGB / maser scale)")
    ax.set_ylabel("expansion velocity f(z) cz$_{CMB}$  (km/s)")
    ax.set_title("38,053 groups with distances that ignore redshift")
    ax.legend(loc="upper left", markerscale=3)
    # (b) log H against velocity
    ax = axes[1]
    for version, color, lab in (("published", ps.C_PUBLISHED_TABLE, "CF4 table as published"),
                                ("corrected", ps.C_CORRECTED, "SDSS sign corrected")):
        b = hb[version]["binned"]
        keep = np.array(b["n"]) >= 30
        v = np.array(b["v_mid"])[keep]
        m = np.array(b["mean_logH"])[keep]
        e = np.array(b["sem_logH"])[keep]
        ax.errorbar(v / 1000, 10**m, yerr=10**m * np.log(10) * e, fmt="o", ms=3.5, color=color, lw=1, capsize=0,
                    label=lab)
    ax.axhline(H0c, color=ps.C_OURS, lw=1.2, label=f"our H0: {H0c:.2f}")
    ax.axhline(74.6, color=ps.C_PUBLISHED, lw=1.0, ls=(0, (4, 2)), label="Tully et al. 2023: 74.6")
    ax.axvline(config.HUBBLE_VMIN / 1000, color=ps.MUTED, lw=0.8)
    ax.text(config.HUBBLE_VMIN / 1000 + 0.3, 63.2, "fit beyond 4,000 km/s", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_ylim(64, 82)
    ax.set_xlim(0, 31)
    ax.set_xlabel("CMB-frame velocity (1,000 km/s)")
    ax.set_ylabel("mean Hubble parameter in 1,000 km/s bins (km/s/Mpc)")
    ax.set_title("Nearby, our neighbours' motions dominate")
    ax.legend(loc="upper right")
    # (c) the Hubble constant in context
    ax = axes[2]
    rows = [
        ("CF4, Tully et al. 2023", 74.6, 0.8, 3.0, ps.C_PUBLISHED),
        ("CF4 table as published", hb["published"]["baseline"]["H0"], 0.76, None, ps.C_PUBLISHED_TABLE),
        ("CF4, SDSS sign corrected (ours)", H0c, hb["corrected"]["zero_point"]["H0_zp_err"], None, ps.C_OURS),
        ("D1: Pantheon+ with SH0ES (this atlas)", 73.55, 1.02, None, ps.INK_SECONDARY),
        ("SH0ES, Riess et al. 2022", 73.04, 1.04, None, ps.INK_SECONDARY),
        ("D3: our Planck CMB fit", 66.66, 0.90, None, ps.INK_SECONDARY),
        ("Planck 2018", 67.66, 0.42, None, ps.INK_SECONDARY),
    ]
    for i, (lab, val, err, sys_, col) in enumerate(rows):
        y = len(rows) - 1 - i
        if sys_:
            ax.plot([val - sys_, val + sys_], [y, y], color=col, lw=1, alpha=0.4)
        ax.errorbar(val, y, xerr=err, fmt="o", color=col, ms=5, lw=1.6, capsize=0)
        ax.text(61.2, y + 0.18, lab, fontsize=7.5, color=ps.INK, va="bottom")
    ax.set_yticks([])
    ax.set_xlim(61, 79)
    ax.set_ylim(-0.6, len(rows) - 0.2)
    ax.set_xlabel("H0 (km/s/Mpc)")
    ax.set_title("The distance ladder against the CMB")
    ax.spines["left"].set_visible(False)
    ps.caption(fig, "Cosmicflows-4 groups (Tully et al. 2023). H0 is the mean of log10(f cz / d) beyond 4,000 km/s; "
                    "its error is dominated by the Cepheid/TRGB zero point (+-0.022 mag -> +-0.76). The faint bar on "
                    "the CF4 point is its quoted ~3 km/s/Mpc systematic. D1 and D3 are earlier projects of this atlas.")
    return ps.save(fig, FIG / "01_hubble_diagram")


# --------------------------------------------------------------------------
# 02. The peculiar-velocity estimator and CF4's published column
# --------------------------------------------------------------------------

def figure_estimators(res, art) -> list:
    g = art["groups_published"]
    ok = g["ok"] & (g["V3k"] > config.CF4_RAMP_KMS)
    est = res["estimators"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    ax = axes[0]
    ax.hexbin(g.loc[ok, "v_direct"], g.loc[ok, "Vpwf"], gridsize=110, extent=(-4000, 4000, -4000, 4000),
              bins="log", cmap=ps.density_cmap(), mincnt=1, rasterized=True)
    xx = np.array([-4000, 4000])
    ax.plot(xx, xx, color=ps.C_OURS, lw=1.4, label="natural log (Watkins & Feldman): slope 1")
    ax.plot(xx, xx / np.log(10), color=ps.C_PUBLISHED, lw=1.4, ls=(0, (4, 2)), label="base-10 log: slope 1/2.303")
    ax.set_xlim(-4000, 4000)
    ax.set_ylim(-4000, 4000)
    ax.set_xlabel("direct estimate cz - H0 d (km/s)   [needs no logarithm]")
    ax.set_ylabel("CF4 Table 4 'Vpwf' column (km/s)")
    r = est["ratio_to_direct_small_v"]
    ax.set_title("CF4's published velocities are 2.3 times too small")
    ax.legend(loc="upper left")
    ax.text(0.98, 0.04, f"groups beyond 3,000 km/s: median ratio to the direct estimate\n"
                        f"CF4 column {r['cf4_column']:.3f}; natural-log estimator {r['natural_log']:.3f}",
            transform=ax.transAxes, ha="right", fontsize=7.5, color=ps.INK_SECONDARY)
    ax = axes[1]
    good = ok & (g["e_DMzp"] <= config.HUBBLE_GOOD_EDM)
    bins = np.arange(-3000, 3001, 150)
    wm = est["well_measured_beyond_ramp"]
    ax.hist(g.loc[good, "Vpec"], bins=bins, histtype="step", color=ps.C_PUBLISHED_TABLE, lw=1.6,
            label=f"CF4 'Vpec' column: rms {wm['cf4_Vpec_column']['std']:.0f} km/s")
    ax.hist(g.loc[good, "v_wf"], bins=bins, histtype="step", color=ps.C_OURS, lw=1.6,
            label=f"natural-log estimator: rms {wm['wf_natural_log']['std']:.0f} km/s")
    sig = wm["expected_error_median"]
    xx = np.linspace(-3000, 3000, 400)
    n = good.sum()
    ax.plot(xx, n * 150 / (np.sqrt(2 * np.pi) * sig) * np.exp(-0.5 * (xx / sig) ** 2), color=ps.INK_SECONDARY,
            lw=1, ls=(0, (2, 2)), label=f"distance errors alone: {sig:.0f} km/s")
    ax.set_xlabel("peculiar velocity (km/s)")
    ax.set_ylabel("groups per 150 km/s")
    ax.set_title(f"{n} groups with distances better than 5%")
    ax.legend(loc="upper right")
    ps.caption(fig, "CF4 Table 4 reproduced to 0.5 km/s with Eq. 10 evaluated in base 10; Watkins & Feldman (2015) derive "
                    "the estimator with the natural logarithm. The published scatter (~300 km/s) is smaller than the "
                    "distance errors alone would produce, which no unbiased estimator can do.")
    return ps.save(fig, FIG / "02_velocity_estimators")


# --------------------------------------------------------------------------
# 03. The sign of the SDSS distances
# --------------------------------------------------------------------------

def figure_sdss_sign(res, art) -> list:
    import cf4fix

    ev = res["sdss_sign"]["evidence"]
    gal = cf4fix._matched()
    s = gal[gal["in_sdss"]]
    dmz = cg.dm_from_redshift(s["zcmb"], s["zhelio"], 75.0)
    h = cf4fix.howlett()
    A = np.column_stack([h["s"], h["i"], np.ones(len(h))])
    coef, *_ = np.linalg.lstsq(A, h["r"], rcond=None)
    fpres = h["r"] - A @ coef
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0))
    ax = axes[0]
    ax.hexbin(s["logdist_corr"], s["DMfp_published"] - dmz if "DMfp_published" in s else s["DMfp"] - dmz,
              gridsize=90, extent=(-0.4, 0.4, -2.2, 2.2), bins="log", cmap=ps.density_cmap(), mincnt=1,
              rasterized=True)
    xx = np.array([-0.4, 0.4])
    c = ev["cf4_vs_eta"]
    ax.plot(xx, c["intercept"] + c["slope"] * xx, color=ps.C_PUBLISHED_TABLE, lw=1.6,
            label=f"CF4 table: DM - DM$_z$ = {c['slope']:+.2f} $\\eta$ (r = {c['r']:.4f})")
    ax.plot(xx, -5 * xx, color=ps.C_OURS, lw=1.6, label="what $\\eta$ = log$_{10}$(d$_z$/d$_H$) requires: -5 $\\eta$")
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-2.2, 2.2)
    ax.set_xlabel("Howlett et al. (2022) log-distance ratio $\\eta$")
    ax.set_ylabel("CF4 distance modulus - redshift modulus (mag)")
    ax.set_title(f"{ev['n_matched']:,} SDSS galaxies in CF4")
    ax.legend(loc="lower left", fontsize=7.5)
    ax = axes[1]
    ax.hexbin(fpres, h["logdist_corr"], gridsize=90, extent=(-0.4, 0.4, -0.4, 0.4), bins="log",
              cmap=ps.density_cmap(), mincnt=1, rasterized=True)
    e2 = ev["eta_vs_fp_residual"]
    ax.plot(xx, e2["slope"] * xx, color=ps.C_OURS, lw=1.6, label=f"slope {e2['slope']:.2f}, r = {e2['r']:.2f}")
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.4, 0.4)
    ax.set_xlabel("size residual from the Fundamental Plane, log$_{10}$ (R$_{obs}$ / R$_{FP}$)")
    ax.set_ylabel("Howlett et al. $\\eta$")
    ax.set_title("$\\eta$ > 0: looks too big at its redshift distance")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.text(0.98, 0.04, "sizes are computed at the group's redshift distance\n"
                        f"(reproduced to {ev['size_at_zcmb_group']['mad']:.4f} dex), so a galaxy that looks too\n"
                        "big is closer than its redshift says: d$_H$ < d$_z$, DM < DM$_z$",
            transform=ax.transAxes, ha="right", fontsize=7.2, color=ps.INK_SECONDARY)
    ax = axes[2]
    comp = res["nusser_tully_comparison"]
    for version, color, mk in (("published", ps.C_PUBLISHED_TABLE, "s"), ("corrected", ps.C_OURS, "o")):
        rows = comp[version]["rows"]
        x = np.concatenate([r["nt"] for r in rows])
        y = np.concatenate([r["ours"] for r in rows])
        lab = ("CF4 table as published" if version == "published" else "SDSS sign corrected (ours)")
        ax.scatter(x, y, s=22, marker=mk, color=color, label=f"{lab}: rms difference {comp[version]['rms_diff_sigma']:.2f} sigma",
                   zorder=3, edgecolors=ps.SURFACE, linewidths=0.6)
    ax.plot([-4.5, 3.5], [-4.5, 3.5], color=ps.MUTED, lw=0.8)
    ax.set_xlim(-4.5, 3.5)
    ax.set_ylim(-4.5, 3.5)
    ax.set_xlabel("Nusser & Tully (2026): standardised deviation")
    ax.set_ylabel("ours: standardised deviation")
    ax.set_title("Reproducing a 2026 bulk-flow analysis")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.text(0.03, 0.97, "monopole and three dipole components\nin six radial bins, 40-300 h$^{-1}$ Mpc",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    ps.caption(fig, "Left: CF4's SDSS moduli follow +5 eta. Middle: geometry fixes the sign of eta (Howlett et al. 2022 "
                    "catalogue, Zenodo 6824749). Right: with the sign corrected, our fits match Nusser & Tully's, who used a "
                    "later CF4 file; with the table as published they do not.")
    return ps.save(fig, FIG / "03_sdss_sign")


# --------------------------------------------------------------------------
# 04. Sky maps of the flows
# --------------------------------------------------------------------------

def _moll_xy(lon_deg, lat_deg):
    """Mollweide axes coordinates with longitude increasing to the left, as on the sky."""
    lon = (np.asarray(lon_deg) + 180.0) % 360.0 - 180.0
    return -np.radians(lon), np.radians(np.asarray(lat_deg))


def _moll_axes(ax):
    ax.grid(True, color=ps.GRID, lw=0.5)
    ticks = np.radians([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150])
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{int(-np.degrees(t)) % 360}" for t in ticks], fontsize=6.5, color=ps.MUTED)
    ax.set_yticks(np.radians([-60, -30, 0, 30, 60]))
    ax.set_yticklabels(["-60", "-30", "0", "30", "60"], fontsize=6.5, color=ps.MUTED)


def figure_sky_maps(res, art) -> list:
    import healpy as hp

    maps = art["sky_maps"]
    lm = art["landmarks"]
    lon = np.linspace(-180, 180, 721)
    lat = np.linspace(-90, 90, 361)
    LON, LAT = np.meshgrid(lon, lat)
    vec = catalogs.unit_vector(LON.ravel(), LAT.ravel())
    pix = hp.vec2pix(config.SKY_NSIDE, vec[:, 0], vec[:, 1], vec[:, 2])
    gl = np.linspace(0, 360, 721)
    gp_sg = catalogs.gal_to_sg(catalogs.unit_vector(gl, np.zeros_like(gl)))
    gp_lon, gp_lat = catalogs.lonlat(gp_sg)
    lg = catalogs.local_group_motion()
    lg_sgl, lg_sgb = catalogs.lonlat(catalogs.gal_to_sg(lg["vector_gal"]))
    fig = plt.figure(figsize=(15.5, 10.2))
    cmap = ps.velocity_cmap()
    cmap.set_bad(ps.GRID)
    titles = ["0-30 h$^{-1}$ Mpc: our neighbourhood", "30-60: Great Attractor, Perseus-Pisces",
              "60-100: Coma, Ophiuchus", "100-160: out to Shapley", "3-30 h$^{-1}$ Mpc, Local Sheet frame"]
    limits = [1000, 1000, 1000, 1000, 600]
    norm_last = None
    for i, (m, title, lim) in enumerate(zip(maps, titles, limits)):
        ax = fig.add_subplot(3, 2, i + 1, projection="mollweide")
        val = m["mean"].copy()
        val[(m["err"] > 400) | ~np.isfinite(val)] = np.nan
        img = val[pix].reshape(LON.shape)
        X, Y = _moll_xy(LON, LAT)
        # pcolormesh needs monotonic x: flip the left-increasing longitude grid.
        norm = TwoSlopeNorm(0, -lim, lim)
        ax.pcolormesh(X[:, ::-1], Y[:, ::-1], img[:, ::-1], cmap=cmap, norm=norm, shading="auto", rasterized=True)
        o = np.argsort(_moll_xy(gp_lon, gp_lat)[0])
        gx, gy = _moll_xy(gp_lon, gp_lat)
        ax.plot(gx[o], gy[o], color=ps.INK_SECONDARY, lw=0.8, ls=(0, (2, 2)))
        lo, hi = m["shell"]
        for _, row in lm.iterrows():
            if lo - 5 <= row["r"] <= hi + 5:
                x0, y0 = _moll_xy(row["SGL"], row["SGB"])
                ax.plot(x0, y0, "o", ms=4, mfc="none", mec=ps.INK, mew=0.9)
                ax.annotate(row["name"], (x0, y0), xytext=(4, 3), textcoords="offset points", fontsize=7, color=ps.INK)
        x0, y0 = _moll_xy(lg_sgl, lg_sgb)
        ax.plot(x0, y0, marker="*", ms=10, color=ps.INK, mec=ps.SURFACE, mew=0.7)
        if i == 0:
            ax.annotate("Local Group's\nmotion (CMB)", (x0, y0), xytext=(6, -16), textcoords="offset points",
                        fontsize=7, color=ps.INK)
        _moll_axes(ax)
        ax.set_title(f"{title}   ({m['n']:,} groups)", fontsize=9.5, pad=8)
    tax = fig.add_subplot(3, 2, 6)
    tax.axis("off")
    for k, (lim, lab) in enumerate(((1000, "CMB frame (first four maps)"), (600, "Local Sheet frame (bottom left)"))):
        cax = tax.inset_axes([0.05, 0.92 - 0.27 * k, 0.6, 0.05])
        sm = plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -lim, lim), cmap=cmap)
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(f"mean peculiar velocity (km/s), {lab}", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    txt = ("Each map is one shell of distance around us, the whole sky in supergalactic\n"
           "coordinates; the dotted line is the Milky Way's disc. Colour: the average\n"
           f"peculiar velocity of Cosmicflows-4 groups within {config.SKY_SMOOTH_DEG:g} degrees; red recedes\n"
           "faster than the expansion, blue lags. Grey: too few groups for a mean\n"
           "better than 400 km/s.\n\n"
           "CMB frame: the nearest shell is a dipole, because we and our neighbours\n"
           "stream together toward the star (the Local Group's 620 km/s): galaxies\n"
           "ahead of us look red, those behind blue. Deeper shells keep a dipole of\n"
           "the same sign: the whole volume shares the flow toward Centaurus and\n"
           "Shapley (the bulk flow of figure 08).\n\n"
           "Local Sheet frame: our shared motion removed. A large dipole remains,\n"
           "groups toward supergalactic south approaching and a patch of the north\n"
           "receding; Tully et al. (2008) traced this to our sheet being pushed\n"
           "away from the Local Void (259 km/s) and pulled toward Virgo (185 km/s).")
    tax.text(0.02, 0.42, txt, fontsize=8.4, va="top", color=ps.INK_SECONDARY, linespacing=1.3, transform=tax.transAxes)
    ps.caption(fig, "Cosmicflows-4 groups (SDSS sign corrected) placed at their CMB-frame redshift distance; peculiar "
                    "velocities from the Watkins & Feldman estimator with H0 = 75 (CF4's scale); error-weighted means with "
                    "250 km/s of intrinsic scatter; 59 catastrophic distance errors excluded.")
    return ps.save(fig, FIG / "04_sky_maps")


# --------------------------------------------------------------------------
# 05. Observed against predicted flows in the supergalactic plane
# --------------------------------------------------------------------------

def _smooth_plane(xg, px, py, vals, errs, sigma=6.0, sig_star=250.0):
    from scipy.spatial import cKDTree

    tree = cKDTree(np.column_stack([px, py]))
    XX, YY = np.meshgrid(xg, xg, indexing="xy")
    grid = np.column_stack([XX.ravel(), YY.ravel()])
    mean = np.full(len(grid), np.nan)
    err = np.full(len(grid), np.nan)
    w0 = 1.0 / (errs**2 + sig_star**2)
    for k, nb in enumerate(tree.query_ball_point(grid, 3 * sigma)):
        if len(nb) < 3:
            continue
        d2 = np.sum((np.column_stack([px[nb], py[nb]]) - grid[k]) ** 2, axis=1)
        kw = np.exp(-0.5 * d2 / sigma**2)
        w = kw * w0[nb]
        mean[k] = np.sum(w * vals[nb]) / np.sum(w)
        err[k] = np.sqrt(np.sum(w**2 / w0[nb])) / np.sum(w)
    return mean.reshape(XX.shape), err.reshape(XX.shape)


def figure_sg_flows(res, art) -> list:
    import gravity

    g = art["groups"]
    rec = art["rec"]
    vc = art["vc"]
    beta = vc["beta_grid"]
    _, cell = dn.grid_axes()
    v1 = dn.velocity_field(dn.field_for(rec, beta), 1.0, cell)
    s = g[g["map_ok"] & (g["SGZ_z"].abs() < config.SG_SLAB_HALF) & (g["r_z"] < config.RECON_RMAX)].copy()
    s["u_pred"] = gravity.predicted_at_redshift_position(s, v1, beta, np.asarray(vc["vext"]))
    xg = np.arange(-150, 151, 2.0)
    obs, eobs = _smooth_plane(xg, s["SGX_z"].to_numpy(), s["SGY_z"].to_numpy(), s["v_wf"].to_numpy(),
                              s["e_v"].to_numpy())
    pred, _ = _smooth_plane(xg, s["SGX_z"].to_numpy(), s["SGY_z"].to_numpy(), s["u_pred"].to_numpy(),
                            s["e_v"].to_numpy())
    bad = ~np.isfinite(eobs) | (eobs > 300)
    obs[bad] = np.nan
    pred[bad] = np.nan
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.4))
    cmap = ps.velocity_cmap()
    cmap.set_bad(ps.SURFACE)
    norm = TwoSlopeNorm(0, -900, 900)
    for ax, img, title in ((axes[0], obs, f"Measured: {len(s):,} Cosmicflows-4 groups within 15 h$^{{-1}}$ Mpc of the plane"),
                           (axes[1], pred, "Predicted from 2MRS gravity (+ the flow from beyond)")):
        ax.pcolormesh(xg, xg, img, cmap=cmap, norm=norm, shading="nearest", rasterized=True)
        zoa_wedge_sg(ax, 150)
        draw_landmarks(ax, art["landmarks"], rmax=150)
        milky_way(ax)
        sg_axes(ax, 150)
        ax.set_title(title, fontsize=10)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.025, pad=0.02)
    cb.set_label("radial peculiar velocity, CMB frame (km/s): red = receding, blue = approaching")
    ok = np.isfinite(obs) & np.isfinite(pred)
    r = np.corrcoef(obs[ok], pred[ok])[0, 1]
    axes[1].text(0.02, 0.02, f"pixel-by-pixel correlation with the measured map: r = {r:.2f}\n"
                             f"beta* = {res['velocity_comparison']['beta']:.2f}, V_ext = {res['velocity_comparison']['Vext']:.0f} km/s",
                 transform=axes[1].transAxes, fontsize=7.8, color=ps.INK_SECONDARY)
    ps.caption(fig, "Both maps use the same groups, positions (CMB-frame redshift distance) and Gaussian smoothing "
                    "(6 h-1 Mpc, error-weighted); blank where the measured mean is uncertain by more than 300 km/s. "
                    "Infall shows as red in front of a mass concentration and blue behind it.")
    return ps.save(fig, FIG / "05_measured_vs_predicted")


# --------------------------------------------------------------------------
# 06. The pull of the galaxies: the 2MRS flux dipole
# --------------------------------------------------------------------------

def figure_dipole(res, art) -> list:
    fd = res["flux_dipole"]
    pub = fd["published"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9), gridspec_kw={"width_ratios": [1, 1, 1.2]})
    ax = axes[0]
    for fr, color, lab in (("lg", ps.C_OURS, "distances from Local Group-frame redshifts"),
                           ("cmb", ps.BLUE, "... from CMB-frame redshifts")):
        ax.plot(fd[fr]["R"], fd[fr]["misalignment"], color=color, lw=1.8, label=lab)
    ax.errorbar([50, 130], [pub["misalign_50"][0], pub["misalign_130"][0]], yerr=[pub["misalign_50"][1], pub["misalign_130"][1]],
                fmt="s", color=ps.C_PUBLISHED, ms=5, lw=1.2, label="Erdogdu et al. 2006 (Ks < 11.25)")
    ax.set_xlabel("galaxies included out to R (h$^{-1}$ Mpc)")
    ax.set_ylabel("angle between the pull and our CMB motion (deg)")
    ax.set_ylim(0, 60)
    ax.set_title("The pull points where we move")
    ax.legend(loc="upper right")
    ax = axes[1]
    for fr, color in (("lg", ps.C_OURS), ("cmb", ps.BLUE)):
        R = np.array(fd[fr]["R"])
        amp = np.array(fd[fr]["amp"])
        ax.plot(R, fd[fr]["lg_cmb_speed"] / amp, color=color, lw=1.8)
    ax.axhspan(pub["beta_L"][0] - pub["beta_L"][1], pub["beta_L"][0] + pub["beta_L"][1], color=ps.C_PUBLISHED, alpha=0.12)
    ax.axhline(pub["beta_L"][0], color=ps.C_PUBLISHED, lw=1, label="Erdogdu et al. 2006: 0.40 +- 0.09")
    ax.set_ylim(0, 1.5)
    ax.set_xlabel("galaxies included out to R (h$^{-1}$ Mpc)")
    ax.set_ylabel("beta$_L$ = f / b$_L$ needed to give 620 km/s")
    ax.set_title("How strong must gravity be?")
    ax.legend(loc="upper right")
    ax.text(0.97, 0.05, "the curve flattens once the pull stops growing\n(the dipole converges)", transform=ax.transAxes,
            ha="right", fontsize=7.5, color=ps.INK_SECONDARY)
    ax = axes[2]
    lg = catalogs.local_group_motion()
    vc = res["velocity_comparison"]
    R = np.array(fd["lg"]["R"])
    sc = ax.scatter(fd["lg"]["l"], fd["lg"]["b"], c=R, cmap=ps.density_cmap(), norm=Normalize(0, 200), s=18,
                    edgecolors=ps.INK_SECONDARY, linewidths=0.3, zorder=3)
    ax.plot(fd["lg"]["l"], fd["lg"]["b"], color=ps.MUTED, lw=0.6, zorder=2)
    for lab, (l_, b_), mk, col in (("Local Group's motion (CMB)", (lg["l"], lg["b"]), "*", ps.INK),
                                   ("Sun's motion (CMB)", (config.SUN_CMB.l, config.SUN_CMB.b), "P", ps.INK_SECONDARY),
                                   ("predicted: 2MRS gravity + external flow", tuple(vc["v_lg_pred_lb"]), "X", ps.C_OURS),
                                   ("2MRS gravity alone", tuple(vc["v_lg_internal_lb"]), "D", ps.C_OURS)):
        ax.plot(l_, b_, marker=mk, ms=11 if mk == "*" else 8, color=col, mec=ps.SURFACE, mew=0.7, ls="none", label=lab, zorder=5)
    for Rlab, off in ((20, (4, 4)), (60, (6, 2)), (130, (6, 2)), (200, (-22, -10))):
        k = int(np.argmin(np.abs(R - Rlab)))
        ax.annotate(f"{Rlab}", (fd["lg"]["l"][k], fd["lg"]["b"][k]), xytext=off, textcoords="offset points",
                    fontsize=7, color=ps.INK_SECONDARY)
    ax.set_xlim(320, 220)
    ax.set_ylim(5, 75)
    ax.set_xlabel("Galactic longitude l (deg)")
    ax.set_ylabel("Galactic latitude b (deg)")
    ax.set_title("Direction of the pull as the survey deepens")
    ax.legend(loc="lower left", fontsize=7.2)
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("galaxies out to R (h$^{-1}$ Mpc)", fontsize=8)
    ps.caption(fig, f"Flux-weighted dipole of {fd['lg']['n_galaxies']:,} 2MRS galaxies (Ks <= 11.75) plus "
                    f"{fd['lg']['n_clones']:,} clones filling the Milky Way's shadow; each galaxy's flux is weighted up for "
                    "fainter galaxies missed at its distance. Because light and gravity both fall as 1/r^2, no distances are "
                    "needed for the direction.")
    return ps.save(fig, FIG / "06_dipole")


# --------------------------------------------------------------------------
# 07. Gravity against the measured motions
# --------------------------------------------------------------------------

def figure_velocity_comparison(res, art) -> list:
    import gravity

    vc = res["velocity_comparison"]
    pub = res["published"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    ax = axes[0]
    prof = vc["profile"]
    b = np.array([p["beta"] for p in prof])
    lnl = np.array([p["lnL"] for p in prof])
    ax.plot(b, lnl.max() - lnl, "o-", color=ps.C_OURS, ms=3.5, lw=1.4, label="Cosmicflows-4 vs 2MRS (ours)")
    ax.axvspan(pub["carrick2015"]["beta"][0] - pub["carrick2015"]["beta"][1],
               pub["carrick2015"]["beta"][0] + pub["carrick2015"]["beta"][1], color=ps.C_PUBLISHED, alpha=0.15,
               label="Carrick et al. 2015 (2M++): 0.431 +- 0.021")
    ax.set_ylim(0, 60)
    ax.set_xlim(0.2, 0.7)
    ax.set_xlabel("beta* = f / b* (growth rate / bias of L* galaxies)")
    ax.set_ylabel("- delta ln(likelihood)")
    ax.set_title(f"beta* = {vc['beta']:.3f} +- {vc['beta_err']:.3f} (statistical)")
    ax.legend(loc="upper center", fontsize=7.5)
    ax = axes[1]
    rec = art["rec"]
    beta = art["vc"]["beta_grid"]
    _, cell = dn.grid_axes()
    v1 = dn.velocity_field(dn.field_for(rec, beta), 1.0, cell)
    import cf4fix

    sdss_ids = set(cf4fix.corrected_groups().query("sdss_corrected")["1PGC"])
    edges = np.arange(-700, 701, 100)
    for version, color, mk, lab in (("corrected", ps.C_OURS, "o", "SDSS sign corrected"),
                                    ("published", ps.C_PUBLISHED_TABLE, "s", "CF4 table as published")):
        smp = gravity.cf4_sample(version)
        up = gravity.predicted_at_redshift_position(smp, v1, beta, np.asarray(art["vc"]["vext"]))
        for subset, alpha, sfx in ((~smp["1PGC"].isin(sdss_ids), 1.0, "groups without SDSS"),
                                   (smp["1PGC"].isin(sdss_ids), 1.0, "SDSS groups")):
            if version == "published" and sfx == "groups without SDSS":
                continue
            m = subset.to_numpy()
            w = 1 / (smp["e_v"].to_numpy() ** 2 + vc["sigma_v"] ** 2)
            idx = np.digitize(up[m], edges)
            xs_, ys_, es_ = [], [], []
            for k in range(1, len(edges)):
                sel = idx == k
                if sel.sum() < 20:
                    continue
                ww = w[m][sel]
                xs_.append(0.5 * (edges[k - 1] + edges[k]))
                ys_.append(np.sum(ww * smp["v_wf"].to_numpy()[m][sel]) / ww.sum())
                es_.append(1 / np.sqrt(ww.sum()))
            filled = sfx == "groups without SDSS"
            ax.errorbar(xs_, ys_, yerr=es_, fmt=mk, color=color, ms=4.5, lw=1, capsize=0,
                        mfc=color if filled else ps.SURFACE, label=f"{sfx}" + ("" if filled else f" ({lab})"))
    ax.plot([-700, 700], [-700, 700], color=ps.MUTED, lw=0.8)
    ax.set_xlim(-700, 700)
    ax.set_ylim(-900, 900)
    ax.set_xlabel("predicted radial velocity (km/s)")
    ax.set_ylabel("measured (error-weighted mean, km/s)")
    ax.set_title("Galaxies move as the 2MRS mass predicts")
    ax.legend(loc="upper left", fontsize=7.2)
    ax = axes[2]
    dec = res["lg_decomposition"]
    lg = catalogs.local_group_motion()
    nhat = lg["vector_gal"] / np.linalg.norm(lg["vector_gal"])
    labels, along, perp = [], [], []
    names = {"[0, 10]": "0-10 h-1 Mpc", "[10, 30]": "10-30 (Virgo side)", "[30, 60]": "30-60 (Great Attr.)",
             "[60, 100]": "60-100 (Coma)", "[100, 160]": "100-160 (Shapley)", "external": "beyond 2MRS"}
    for d in dec:
        v = np.asarray(d["v"])
        labels.append(names.get(str(d["shell"]), str(d["shell"])))
        along.append(float(v @ nhat))
        perp.append(float(np.linalg.norm(v - (v @ nhat) * nhat)))
    tot = np.sum([np.asarray(d["v"]) for d in dec], axis=0)
    labels.append("total predicted")
    along.append(float(tot @ nhat))
    perp.append(float(np.linalg.norm(tot - (tot @ nhat) * nhat)))
    y = np.arange(len(labels))[::-1]
    ax.barh(y, along, color=[ps.C_OURS] * (len(labels) - 1) + [ps.INK], height=0.6)
    ax.axvline(lg["speed"], color=ps.C_PUBLISHED, lw=1.2)
    ax.text(lg["speed"] + 8, y[-1] + 0.45, f"CMB: {lg['speed']:.0f} km/s", fontsize=7.5, color=ps.C_PUBLISHED, va="center")
    for yy, a, p_, lab in zip(y, along, perp, labels):
        ax.text(-10, yy, lab, ha="right", va="center", fontsize=7.6, color=ps.INK)
        ax.text(max(a, 0) + 8, yy, f"{a:.0f} (+{p_:.0f} sideways)", va="center", fontsize=7, color=ps.INK_SECONDARY)
    ax.set_yticks([])
    ax.set_xlim(-360, 860)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("contribution along the Local Group's CMB motion (km/s)")
    ax.set_title("Where our 620 km/s comes from")
    ps.caption(fig, f"Forward likelihood of {vc['n_groups']:,} CF4 groups (redshift distance < {config.VCOMP_RMAX:g} h-1 Mpc), "
                    "integrating each over its true distance with CF4's radial distribution and the 2MRS density as prior. "
                    f"Fitted V_ext = {vc['Vext']:.0f} km/s toward (l, b) = ({vc['Vext_l']:.0f}, {vc['Vext_b']:.0f}), "
                    f"sigma_v = {vc['sigma_v']:.0f} km/s. Predicted Local Group motion {vc['v_lg_pred_speed']:.0f} km/s, "
                    f"{vc['v_lg_misalignment']:.1f} deg from the CMB direction.")
    return ps.save(fig, FIG / "07_gravity_vs_flows")


# --------------------------------------------------------------------------
# 08. Bulk flows against Lambda-CDM
# --------------------------------------------------------------------------

def figure_bulk_flow(res, art) -> list:
    bf = res["bulk_flow"]
    pub = res["published"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.1), gridspec_kw={"width_ratios": [1.2, 1.1, 0.75]})
    ax = axes[0]
    sp = bf["corrected"]["spheres"]
    R = np.array([d["R"] for d in sp])
    band = np.array([d["lcdm_band_V"] for d in sp])
    ax.fill_between(R, band[:, 0], band[:, 2], color=ps.C_LCDM, alpha=0.35, lw=0, label="Lambda-CDM, this window: 68%")
    ax.fill_between(R, band[:, 2], band[:, 3], color=ps.C_LCDM, alpha=0.15, lw=0, label="... 68-97.7%")
    ideal = bf["ideal"]
    ax.plot(ideal["R"], ideal["rms_tophat"], color=ps.INK_SECONDARY, lw=1, ls=(0, (3, 2)),
            label="Lambda-CDM rms, ideal full sphere")
    for version, color, lab in (("published", ps.C_PUBLISHED_TABLE, "CF4 table as published"),
                                ("corrected", ps.C_OURS, "SDSS sign corrected (ours)")):
        s_ = bf[version]["spheres"]
        ax.errorbar([d["R"] for d in s_], [d["Vamp"] for d in s_], yerr=[d["Vamp_err_meas"] for d in s_], color=color,
                    fmt="o-", ms=3, lw=1.2, capsize=0, label=lab)
    for key, lab, off in (("watkins2023", "Watkins+2023", (-66, 6)), ("whitford2023_mle", "Whitford+2023 (MLE)", (5, 12)),
                          ("whitford2023_mve", "Whitford+2023 (MVE)", (6, 10)), ("courtois2025", "Courtois+2025", (-70, -14))):
        p_ = pub[key]
        ax.errorbar(p_["R"], p_["V"], yerr=p_["err"], fmt="s", color=ps.C_PUBLISHED, ms=4.5, lw=1, capsize=0)
        ax.annotate(lab, (p_["R"], p_["V"]), xytext=off, textcoords="offset points", fontsize=7, color=ps.C_PUBLISHED)
    ax.set_xlim(40, 305)
    ax.set_ylim(0, 800)
    ax.set_xlabel("all groups out to R (h$^{-1}$ Mpc; fitted beyond 4,000 km/s)")
    ax.set_ylabel("bulk flow |V| (km/s)")
    ax.set_title("Spheres of growing radius")
    ax.legend(loc="upper right", fontsize=7)
    ax.text(0.02, 0.03, f"inverse-variance weights favour nearby groups: even at R = 300 the\n"
                        f"weighted mean distance is {sp[-1]['reff']:.0f} h-1 Mpc. Published values use other windows.",
            transform=ax.transAxes, fontsize=7, color=ps.INK_SECONDARY)
    ax = axes[1]
    for version, color, mk, off in (("published", ps.C_PUBLISHED_TABLE, "s", -2.5), ("corrected", ps.C_OURS, "o", 2.5)):
        bins = bf[version]["bins"]
        x = np.array([d["reff"] for d in bins.values()]) + off
        y = [d["Vamp"] for d in bins.values()]
        e = [d["Vamp_err_meas"] for d in bins.values()]
        if version == "corrected":
            bands = np.array([d["lcdm_band_V"] for d in bins.values()])
            for xx, bb in zip(x, bands):
                ax.plot([xx - off, xx - off], [bb[0], bb[2]], color=ps.C_LCDM, lw=7, alpha=0.5, solid_capstyle="butt")
                ax.plot([xx - off, xx - off], [bb[2], bb[3]], color=ps.C_LCDM, lw=7, alpha=0.2, solid_capstyle="butt")
            for xx, yy, d in zip(x, y, bins.values()):
                ax.annotate(f"{d['pV_sigma']:.1f} s", (xx, yy), xytext=(6, -3), textcoords="offset points",
                            fontsize=7, color=ps.INK_SECONDARY)
        ax.errorbar(x, y, yerr=e, fmt=mk, color=color, ms=5, lw=1.2, capsize=0,
                    label="SDSS sign corrected (ours)" if version == "corrected" else "CF4 table as published")
    nt = pub["nusser_tully2026"]["bin_120_160"]
    ax.errorbar(138, nt["Vdip"], yerr=nt["err"], fmt="D", color=ps.C_PUBLISHED, ms=5, lw=1.2, capsize=0,
                label="Nusser & Tully 2026 (120-160)")
    ax.set_xlabel("radial bin, weighted mean distance (h$^{-1}$ Mpc)")
    ax.set_ylabel("dipole |V| in the bin (km/s)")
    ax.set_ylim(0, 1000)
    ax.set_title("Separate shells (Nusser & Tully's bins)")
    ax.legend(loc="upper left", fontsize=7)
    ax.text(0.98, 0.03, "grey bars: Lambda-CDM 68% / 97.7% for each bin's own window;\nlabels: deviation of the dipole (two-sided sigma)",
            transform=ax.transAxes, ha="right", fontsize=7, color=ps.INK_SECONDARY)
    ax = axes[2]
    bins = bf["corrected"]["bins"]
    mat = np.array([d["sigma_dev"] for d in bins.values()])
    im = ax.imshow(mat, cmap=ps.velocity_cmap(), norm=TwoSlopeNorm(0, -4, 4), aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:+.1f}", ha="center", va="center", fontsize=7.5,
                    color=ps.SURFACE if abs(mat[i, j]) > 2.5 else ps.INK)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["monopole", "SGX", "SGY", "SGZ"], fontsize=8)
    ax.set_yticks(range(len(bins)))
    ax.set_yticklabels(list(bins.keys()), fontsize=8)
    ax.grid(False)
    ax.set_title("Each component vs Lambda-CDM (sigma)")
    gl = bf["corrected"]["global"]["40-300"]
    ax.set_xlabel(f"all six bins together: dipole {gl['pV_sigma']:.1f} sigma, joint {gl['pmV_sigma']:.1f} sigma", fontsize=8)
    ps.caption(fig, "Joint monopole + dipole fits to 51,538 CF4 galaxies at their redshift distance (Nusser & Tully 2026 "
                    "method); the Lambda-CDM covariance (Planck 2018, linear theory) is evaluated at the actual positions and "
                    "weights of every galaxy pair. Negative SGX points toward the Great Attractor and Shapley.")
    return ps.save(fig, FIG / "08_bulk_flows")


# --------------------------------------------------------------------------
# 09. Basins at different smoothing
# --------------------------------------------------------------------------

def figure_basins(res, art) -> list:
    rec = art["rec"]
    beta = art["vc"]["beta_grid"]
    delta = dn.field_for(rec, beta, smooth=5.0)
    half = 158.0
    x, X, Y, pts = sg_plane_grid(half, 1.0)
    dens = sample_scalar_sg(delta, pts, X.shape)
    R = np.hypot(X, Y)
    one = np.clip(1 + dens, 0.15, 12)
    one[R > config.RECON_RMAX] = np.nan
    keys = [("6", "6 h$^{-1}$ Mpc smoothing"), ("10", "10 h$^{-1}$ Mpc"), ("15", "15 h$^{-1}$ Mpc"),
            ("15_ext", "15 h$^{-1}$ Mpc + the pull from beyond 2MRS")]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 12.0))
    fig.subplots_adjust(top=0.94, bottom=0.05, hspace=0.16, wspace=0.12)
    for ax, (k, title) in zip(axes.ravel(), keys):
        bas = art["basins"][k]
        rep = bas["report"]
        ax.pcolormesh(x, x, one, cmap=ps.density_cmap(), norm=LogNorm(0.15, 12), shading="nearest", alpha=0.6,
                      rasterized=True)
        inside = R < config.RECON_RMAX - 2
        for b in rep["largest"][:10]:
            if b["is_ours"]:
                continue
            ind = basin_indicator_sg(bas["labels"], b["label"], pts, X.shape)
            ax.contour(x, x, np.where(inside, ind, np.nan), levels=[0.5], colors=[ps.INK_SECONDARY], linewidths=0.8)
        ind = np.where(inside, basin_indicator_sg(bas["labels"], rep["ours_label"], pts, X.shape), np.nan)
        ax.contourf(x, x, ind, levels=[0.5, 1.5], colors=[ps.C_BASIN], alpha=0.25)
        ax.contour(x, x, ind, levels=[0.5], colors=[ps.C_BASIN], linewidths=2)
        att = rep["ours_attractor_sg"]
        ax.plot(att[0], att[1], marker="X", ms=10, color=ps.C_BASIN, mec=ps.INK, mew=0.8, zorder=9)
        ax.add_patch(plt.Circle((0, 0), config.RECON_RMAX, fill=False, color=ps.MUTED, lw=0.8, ls=(0, (2, 3))))
        draw_landmarks(ax, art["landmarks"], fontsize=7)
        milky_way(ax)
        sg_axes(ax, half)
        members = [m for m, inn in rep["members"].items() if inn]
        ax.set_title(title, fontsize=10.5)
        ax.text(0.01, 0.01, f"our basin: {rep['ours_equiv_diameter'] * 100:,.0f} km/s across, "
                            f"{rep['ours_mass_msun'] / 1e17:.1f} x 10$^{{17}}$ M$_\\odot$\n"
                            f"attractor {rep['ours_attractor_r']:.0f} h$^{{-1}}$ Mpc toward (l, b) = "
                            f"({rep['ours_attractor_lb'][0]:.0f}, {rep['ours_attractor_lb'][1]:.0f})\n"
                            f"contains: {', '.join(members)}",
                transform=ax.transAxes, fontsize=7.5, color=ps.INK, va="bottom",
                bbox={"facecolor": ps.SURFACE, "edgecolor": "none", "alpha": 0.85, "pad": 2})
    fig.suptitle("Which basin are we in? It depends on the scale you ask about", x=0.06, ha="left", y=0.985)
    ps.caption(fig, "Basins are catchment areas of the linear-theory velocity potential of the 2MRS mass inside "
                    f"{config.RECON_RMAX:g} h-1 Mpc (watershed on a 3.1 h-1 Mpc grid; supergalactic plane shown). "
                    "Bottom right adds the uniform flow fitted from beyond 2MRS, which tilts the landscape toward Shapley.")
    return ps.save(fig, FIG / "09_basins")


ALL = [figure_laniakea, figure_hubble, figure_estimators, figure_sdss_sign, figure_sky_maps, figure_sg_flows,
       figure_dipole, figure_velocity_comparison, figure_bulk_flow, figure_basins]
