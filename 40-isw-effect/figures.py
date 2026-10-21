"""Every figure of the project, from data/results.json and the cached derived arrays.

    python figures.py

All figures are PNG + SVG in figures/. Nothing here recomputes an analysis
step except cheap map projections.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter

import plotstyle as P
from config import (DERIVED_DIR, FIGURE_DIR, LITERATURE, RESULTS_PATH, SCAN_N, SCAN_RADII_DEG,
                    STACK_IMAGE_DEG, STACK_IMAGE_PIX)

P.use_style()
CMAP = P.diverging_cmap()
FIG = FIGURE_DIR


def load():
    return json.loads(RESULTS_PATH.read_text())


def _cbar(fig, im, ax, label):
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(label, color=P.INK_SECONDARY, fontsize=8)
    cb.ax.tick_params(labelsize=7, colors=P.MUTED)
    cb.outline.set_visible(False)
    return cb


def _img_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def ortho_patch(m, centre_lb, radius_deg, npix=400, mask=None):
    """Orthographic view of a HEALPix map around a centre (masked pixels NaN)."""
    import healpy as hp

    nside = hp.get_nside(m)
    s = np.sin(np.radians(radius_deg))
    x = np.linspace(-s, s, npix)
    xx, yy = np.meshgrid(x, x)
    rr = xx**2 + yy**2
    inside = rr <= s**2
    zz = np.sqrt(np.clip(1 - rr, 0, 1))
    lo, la = np.radians(centre_lb[0]), np.radians(centre_lb[1])
    er = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    el = np.array([-np.sin(lo), np.cos(lo), 0.0])
    eb = np.cross(er, el)
    v = xx[..., None] * (-el) + yy[..., None] * eb + zz[..., None] * er  # east to the left
    pix = hp.vec2pix(nside, v[..., 0], v[..., 1], v[..., 2])
    out = m[pix].astype(float)
    if mask is not None:
        out[~mask[pix]] = np.nan
    out[~inside] = np.nan
    return out


# --------------------------------------------------------------------------
# 1. The idea
# --------------------------------------------------------------------------


def fig_idea(r):
    import coldspot as CS
    import cosmology as C

    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.8))
    z = np.linspace(0, 3, 301)
    ax = axs[0]
    g = C.LCDM.g(z) / C.LCDM.g(3.0)
    ax.plot(z, g, color=P.C_MODEL, label="Our universe (D3's ΛCDM)")
    ax.plot(z, np.ones_like(z), color=P.C_EDS, ls="--", label="Flat, matter only")
    ax.set_xlim(3, 0)
    ax.set_xlabel("Redshift z (time runs to the right)")
    ax.set_ylabel("Depth of a potential well\n(relative to z = 3)")
    ax.set_title("Dark energy makes wells shallower")
    ax.legend(loc="lower left")
    zt = 0.0
    P.annotate(ax, f"today: {g[0]:.0%} of its depth at z = 3", (zt, g[0]), (1.6, 0.86))
    ax.set_ylim(0.7, 1.05)

    ax = axs[1]
    th = np.linspace(0, 40, 401)
    d0, r0, zc, q = CS.SUPERVOID_MODELS["Finelli+16"]
    prof = CS.supervoid_profile(d0, r0, zc, th, q)
    ax.plot(th, prof, color=P.C_MODEL, label="ΛCDM: photons leave colder")
    ax.plot(th, np.zeros_like(th), color=P.C_EDS, ls="--", label="Matter only: no effect")
    ax.set_xlabel("Angle from the supervoid's centre (deg)")
    ax.set_ylabel("ISW temperature change (µK)")
    ax.set_title(f"A supervoid's ISW shadow (r₀ = {r0:.0f} h⁻¹ Mpc, z = {zc})")
    ax.legend(loc="lower right")
    P.caption(fig, "Left: g = D(z)(1+z), the potential's growth, from our own growth calculation "
                   "with D3's parameters. Right: linear ISW profile of the compensated Eridanus "
                   "supervoid model of Finelli et al. (2016), integrated along each line of sight.")
    return P.save(fig, FIG / "01_the_isw_idea")


# --------------------------------------------------------------------------
# 2. Granett's stack, replicated, and what LCDM predicts there
# --------------------------------------------------------------------------


def fig_granett_stack(r):
    im = np.load(DERIVED_DIR / "stack_images.npz")
    st = r["stacking"]["stacks"]
    fig, axs = plt.subplots(2, 3, figsize=(11, 7.4))
    half = STACK_IMAGE_DEG / 2
    pix = STACK_IMAGE_DEG / STACK_IMAGE_PIX
    sig = 1.4 / 2.355 / pix  # Granett et al. displayed with a 1.4 deg FWHM smoothing
    lim = 15
    for row, sky in enumerate(("smica", "pred")):
        v = im[f"granett_voids__{sky}"]
        c = im[f"granett_clusters__{sky}"]
        panels = [(v, "50 supervoids"), (c, "50 superclusters"), ((c - v) / 2, "Combined: (clusters − voids)/2")]
        for col, (img, title) in enumerate(panels):
            ax = axs[row, col]
            sm = gaussian_filter(np.nan_to_num(img), sig)
            h = ax.imshow(sm, cmap=CMAP, vmin=-lim, vmax=lim, origin="lower",
                          extent=[half, -half, -half, half])
            for rad, ls in ((4.0, "-"), (4.0 * np.sqrt(2), "--")):
                ax.add_patch(Circle((0, 0), rad, fill=False, color=P.INK, lw=0.8, ls=ls))
            _img_axes(ax)
            key = ["granett_voids", "granett_clusters", "granett_combined"][col]
            if sky == "smica":
                e = st[key]["smica"]
                ax.set_title(f"{title}\nPlanck: {e['value']:+.1f} ± {e['sigma']:.1f} µK", fontsize=9.5)
            else:
                p = st[key]["predicted_uK"]["fiducial"]
                ax.set_title(f"ΛCDM prediction from BOSS: {p:+.2f} µK", fontsize=9.5)
    fig.subplots_adjust(left=0.08, right=0.88, top=0.86, bottom=0.08, hspace=0.32, wspace=0.12)
    cax = fig.add_axes([0.905, 0.25, 0.014, 0.5])
    cb = fig.colorbar(h, cax=cax)
    cb.set_label("Temperature (µK), multipoles below 10 removed", color=P.INK_SECONDARY, fontsize=8)
    cb.ax.tick_params(labelsize=7, colors=P.MUTED)
    cb.outline.set_visible(False)
    axs[0, 0].text(-0.08, 0.5, "Planck SMICA (real sky)", transform=axs[0, 0].transAxes,
                   rotation=90, va="center", ha="right", fontsize=10, fontweight="bold")
    axs[1, 0].text(-0.08, 0.5, "Predicted ISW sky (BOSS galaxies)", transform=axs[1, 0].transAxes,
                   rotation=90, va="center", ha="right", fontsize=10, fontweight="bold")
    fig.suptitle("Granett et al. (2008)'s 50 supervoids and 50 superclusters, stacked on the CMB",
                 x=0.02, ha="left")
    P.caption(fig, "24° × 24° patches, smoothed by 1.4° for display (as in Granett et al.). Circles: "
                   "the 4° compensated filter (solid) and its equal-area ring (dashed). Values: mean "
                   "filter output ± scatter in 1,000 simulated skies. Same colour scale in both rows.")
    return P.save(fig, FIG / "02_granett_stack")


# --------------------------------------------------------------------------
# 3. Filter size, and the price of choosing it
# --------------------------------------------------------------------------


def fig_scan(r):
    sc = r["stacking"]["scans"]["granett_combined"]
    le = r["stacking"]["look_elsewhere"]
    x = np.array(sc["x"])
    d, sd, pr = (np.array(sc[k]) for k in ("data", "sigma", "pred"))
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axs[0]
    ax.fill_between(x, -2 * sd, 2 * sd, color=P.GRID, lw=0, label="Chance (±2σ, simulated skies)")
    ax.fill_between(x, -sd, sd, color=P.AXIS, lw=0, alpha=0.6, label="Chance (±1σ)")
    ax.plot(x, d, color=P.C_OURS, marker="o", ms=4, label="Planck SMICA (this project)")
    ax.plot(x, sc["data_wmap"], color=P.C_WMAP, ls="--", label="WMAP 9-yr ILC (this project)")
    t1 = LITERATURE["granett2008"]["table1"]
    gx = sorted(k[1] for k in t1 if k[0] == 50)
    ax.plot(gx, [t1[(50, g)] for g in gx], ls="none", marker="D", ms=5, color=P.C_PUBLISHED,
            label="Granett+08 (WMAP5), their Table 1")
    ax.plot(x, pr, color=P.C_MODEL, label="ΛCDM prediction (BOSS)")
    ax.axhline(0, color=P.C_EDS, ls=":", lw=1.2)
    ax.text(9.9, 0.4, "matter only: 0", color=P.C_EDS, ha="right", fontsize=8)
    ax.axvline(4.0, color=P.INK_SECONDARY, lw=0.6)
    ax.set_xlabel("Compensated-filter radius (deg)")
    ax.set_ylabel("Combined signal (µK)")
    ax.set_title("The signal peaks at the radius Granett et al. chose")
    ax.legend(loc="lower left", fontsize=7.5)

    ax = axs[1]
    ms = np.load(DERIVED_DIR / "lookelsewhere_max_snr_sims.npy")
    ax.hist(ms, bins=40, color=P.AXIS, edgecolor=P.SURFACE)
    ax.axvline(le["best_snr"], color=P.C_OURS, lw=1.6)
    ax.text(le["best_snr"] + 0.05, ax.get_ylim()[1] * 0.85,
            f"Planck, best of the grid\nS/N {le['best_snr']:.2f}\n"
            f"chance does this {le['p_global']:.1%}\nof the time ({le['global_sigma']:.1f}σ)",
            color=P.C_OURS, fontsize=8, va="top")
    ax.set_xlabel(f"Largest S/N found on the grid ({len(SCAN_RADII_DEG)} radii × {len(SCAN_N)} sample sizes)")
    ax.set_ylabel("Simulated skies")
    ax.set_title("Searching many choices finds 'signals' in noise")
    P.caption(fig, f"Right: every one of 1,000 pure-noise skies was searched over the same grid of "
                   f"filter radii (1–10°) and numbers of structures (10–50) as the data. At the single "
                   f"published choice (4°, 50+50) the Planck S/N is {le['fixed_snr']:.2f}.")
    return P.save(fig, FIG / "03_filter_scan_look_elsewhere")


# --------------------------------------------------------------------------
# 4. All stacks
# --------------------------------------------------------------------------

STACK_LABELS = {
    "granett_voids": "Granett+08: 50 supervoids (4°)",
    "granett_clusters": "Granett+08: 50 superclusters (4°)",
    "nh14_lrgdim_voids": "NH14 lrgdim: 70 voids",
    "nh14_lrgdim_clusters": "NH14 lrgdim: 196 superclusters",
    "nh14_lrgbright_voids": "NH14 lrgbright: 13 voids",
    "nh14_lrgbright_clusters": "NH14 lrgbright: 39 superclusters",
    "mao_voids": "Mao+17 BOSS: 1,228 voids",
}


def fig_stacks(r):
    st = r["stacking"]["stacks"]
    hot = LITERATURE["hotchkiss2015"]
    pub = {"granett_voids": LITERATURE["granett2008"]["voids_uK"],
           "granett_clusters": LITERATURE["granett2008"]["clusters_uK"],
           "nh14_lrgdim_voids": hot["lrgdim_voids_uK"][:2],
           "nh14_lrgdim_clusters": hot["lrgdim_clusters_uK"][:2],
           "nh14_lrgbright_voids": hot["lrgbright_voids_uK"][:2],
           "nh14_lrgbright_clusters": hot["lrgbright_clusters_uK"][:2]}
    keys = list(STACK_LABELS)
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    y = np.arange(len(keys))[::-1]
    for yi, k in zip(y, keys):
        e = st[k]["smica"]
        col = P.C_VOID if "void" in k else P.C_CLUSTER
        ax.errorbar(e["value"], yi + 0.12, xerr=e["sigma"], fmt="o", color=col, ms=6, capsize=0,
                    lw=1.6)
        ax.plot(st[k]["predicted_uK"]["fiducial"], yi + 0.12, marker="|", ms=14, mew=2,
                color=P.C_MODEL)
        if k in pub:
            ax.errorbar(pub[k][0], yi - 0.18, xerr=pub[k][1], fmt="D", ms=4, color=P.C_PUBLISHED,
                        lw=1, capsize=0)
    ax.axvline(0, color=P.C_EDS, ls=":", lw=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels([STACK_LABELS[k] for k in keys])
    ax.set_xlabel("Mean compensated-filter temperature (µK)")
    ax.set_title("Every catalogue, its ΛCDM prediction, and the published value")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=P.C_VOID, marker="o", ls="none", label="Planck: voids"),
               Line2D([], [], color=P.C_CLUSTER, marker="o", ls="none", label="Planck: superclusters"),
               Line2D([], [], color=P.C_MODEL, marker="|", ms=12, mew=2, ls="none",
                      label="ΛCDM prediction (BOSS)"),
               Line2D([], [], color=P.C_PUBLISHED, marker="D", ls="none",
                      label="Published: Granett+08 (WMAP5),\nHotchkiss+15 (Planck 2013)"),
               Line2D([], [], color=P.C_EDS, ls=":", label="Matter only (0)")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7.5)
    ax.grid(axis="y", visible=False)
    P.caption(fig, "Granett+08 at 4°; NH14 (Nadathur & Hotchkiss 2014) and Mao+17 with the filter "
                   "scaled to each structure, R = 0.6 × its angular radius. Error bars: scatter in "
                   "1,000 simulated skies. Voids are expected cold (negative), superclusters hot.")
    return P.save(fig, FIG / "04_all_stacks")


# --------------------------------------------------------------------------
# 0. Summary: the ISW amplitude from every probe
# --------------------------------------------------------------------------


def fig_summary(r):
    mf = r["stacking"]["matched_filter"]
    cc = r["crosscorr"]
    comb = r["combined"]
    rows = [
        ("Granett+08 structures (matched filter)", mf["granett"]["fiducial"], P.C_OURS),
        ("NH14 DR7 LRG voids + superclusters", mf["nh14"]["fiducial"], P.C_OURS),
        ("Mao+17 BOSS voids", mf["mao"]["fiducial"], P.C_OURS),
        ("BOSS predicted-ISW template", cc["template"]["fiducial"], P.C_OURS),
        ("BOSS LOWZ × CMB (0.15 < z < 0.43)", cc["tracers"]["boss_0.15_0.43"], P.C_OURS),
        ("BOSS CMASS × CMB (0.43 < z < 0.75)", cc["tracers"]["boss_0.43_0.75"], P.C_OURS),
        ("NVSS radio sources × CMB", cc["tracers"]["nvss"], P.C_OURS),
        ("Combined: NVSS + BOSS cross-spectra", comb["all_crosscorr"], P.INK),
        ("Combined: NVSS + BOSS template", comb["nvss_template"], P.INK),
        ("Combined: + void matched filters", comb["nvss_template_voids"], P.INK),
    ]
    lit = [("Planck 2015 XXI: NVSS", LITERATURE["planck2015_xxi"]["NVSS_A"]),
           ("Planck 2015 XXI: all tracers", LITERATURE["planck2015_xxi"]["all_A"]),
           ("Nadathur & Crittenden 2016: BOSS voids", LITERATURE["nc2016"]["A_ISW"]),
           ("Kovács+ BOSS/DES supervoids", LITERATURE["kovacs2022"]["supervoid_A_ISW"])]
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    n = len(rows) + len(lit)
    y = np.arange(n)[::-1].astype(float)
    labels = []
    for i, (lab, e, col) in enumerate(rows):
        a, s = e["A"], e["sigma"]
        ax.errorbar(a, y[i], xerr=s, fmt="o", color=col, ms=6 if col == P.INK else 5, lw=1.6,
                    capsize=0)
        ax.text(1.01, y[i], f"{a:+.2f} ± {s:.2f}", va="center", ha="left", fontsize=7.5,
                color=P.INK_SECONDARY, transform=ax.get_yaxis_transform())
        labels.append(lab)
    for j, (lab, (a, s)) in enumerate(lit):
        i = len(rows) + j
        ax.errorbar(a, y[i], xerr=s, fmt="D", color=P.C_PUBLISHED, ms=4, lw=1, capsize=0)
        ax.text(1.01, y[i], f"{a:+.2f} ± {s:.2f}", va="center", ha="left", fontsize=7.5,
                color=P.C_PUBLISHED, transform=ax.get_yaxis_transform())
        labels.append(lab)
    ax.axhline(y[len(rows) - 1] - 0.5, color=P.GRID, lw=0.8)
    ax.axvline(0, color=P.C_EDS, ls="--", lw=1.2)
    ax.axvline(1, color=P.C_MODEL, lw=1.2)
    ax.text(0.03, y[0] + 0.7, "matter only", color=P.C_EDS, fontsize=8)
    ax.text(1.03, y[0] + 0.7, "ΛCDM", color=P.C_MODEL, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(-3.2, 7.2)
    ax.set_ylim(y[-1] - 0.7, y[0] + 1.3)
    ax.set_xlabel("ISW amplitude A (1 = the ΛCDM prediction, 0 = no dark energy)")
    n = cc["tracers"]["nvss"]
    ax.set_title(f"The ISW amplitude from every probe (NVSS alone: {n['A']:.2f} ± {n['sigma']:.2f}, "
                 f"{n['snr']:.1f}σ from no dark energy)")
    ax.grid(axis="y", visible=False)
    P.caption(fig, "Our measurements (circles) with errors from 1,000 simulated CMB skies; "
                   "combinations use the probes' joint covariance from the same skies. Diamonds: "
                   "published values for comparison (different catalogues and methods).")
    return P.save(fig, FIG / "00_summary")


# --------------------------------------------------------------------------
# 5. The BOSS reconstruction
# --------------------------------------------------------------------------


def fig_reconstruction(r):
    import healpy as hp

    import boss
    import catalogs
    import cmb
    from config import GRID_CELL, GRID_SHAPE

    grid = boss.build_grid()
    bias = boss.measure_bias()
    path = DERIVED_DIR / "wf_slab.npy"
    x0 = boss.grid_axes()
    zslab = 900.0  # Mpc/h above the equatorial plane: through the northern cap
    k = int((zslab - x0[2]) / GRID_CELL)
    if path.exists():
        slab = np.load(path)
    else:
        d0 = boss.wiener_filter(grid, bias)
        slab = d0[:, :, k - 1:k + 2].mean(axis=2)
        np.save(path, slab)
    raw = boss.overdensity(grid)[:, :, k - 1:k + 2].mean(axis=2)
    pred = boss.predicted_isw("fiducial")
    d = cmb.prepare_maps()
    import superstructures as SS
    sky = SS.highpass(d["smica"], d["mask"])
    pred_hp = SS.highpass(hp.ud_grade(pred, 512).astype(np.float32), d["mask"])

    fig = plt.figure(figsize=(13.5, 9.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], wspace=0.42, hspace=0.32)
    ext = [x0[0], x0[0] + GRID_SHAPE[0] * GRID_CELL, x0[1], x0[1] + GRID_SHAPE[1] * GRID_CELL]
    cover = grid["cover"].reshape(GRID_SHAPE)[:, :, k] > 0
    for j, (img, title, lim) in enumerate([(raw, "Galaxies: raw overdensity δg (noisy)", 1.5),
                                          (slab, "Wiener-filtered matter δ₀ (today)", 0.5)]):
        ax = fig.add_subplot(gs[0, j])
        show = np.where(cover, img, np.nan) if j == 0 else img
        h = ax.imshow(show.T, origin="lower", extent=ext, cmap=CMAP, vmin=-lim, vmax=lim)
        ax.set_xlim(-2100, 400)
        ax.set_ylim(-1900, 1900)
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("x (h⁻¹ Mpc)")
        ax.set_ylabel("y (h⁻¹ Mpc)")
        ax.grid(False)
        _cbar(fig, h, ax, "overdensity (red = dense)")
        if j == 1:
            m = catalogs.mao()
            xyz = boss.radec_to_xyz(m.ra, m.dec, boss.cosmology.LCDM.chi(m.z))
            sel = (np.abs(xyz[:, 2] - zslab) < 0.5 * m.r_eff) & (m.r_eff < 120)
            for xx, yy, rr in zip(xyz[sel, 0], xyz[sel, 1], m.r_eff[sel]):
                ax.add_patch(Circle((xx, yy), rr, fill=False, color=P.INK, lw=0.6))
    ax = fig.add_subplot(gs[0, 2])
    bres = r["boss"]["bias"]
    for key, col, lab in (("lo", P.BLUE, "LOWZ"), ("hi", P.ORANGE, "CMASS")):
        b = boss.measure_bias()[key]
        kk, pk, pl = b["k"], b["pk"], b["plin"]
        sel = kk < 0.25
        f = b["f"]
        amp = b["b"] ** 2 + 2 * b["b"] * f / 3 + f**2 / 5
        ax.loglog(kk[sel], pk[sel], color=col, marker="o", ms=3, ls="none", label=f"{lab} data")
        ax.loglog(kk[sel], amp * pl[sel], color=col, lw=1.2,
                  label=f"{lab} model: linear ΛCDM, b = {bres[key]['b']:.2f}")
    ax.axvspan(0.02, 0.08, color=P.GRID, lw=0)
    ax.text(0.04, 2.2e3, "fitted", ha="center", fontsize=7.5, color=P.INK_SECONDARY)
    ax.set_ylim(2e3, 5e5)
    ax.set_xlabel("k (h Mpc⁻¹)")
    ax.set_ylabel("P(k) (h⁻³ Mpc³)")
    ax.set_title("Galaxy bias from the power spectrum", fontsize=9.5)
    ax.legend(fontsize=7, loc="upper right")

    # the BOSS footprint in Galactic pixels, for the spectra
    import crosscorr as X
    bmask = hp.ud_grade(X.tracers()["boss_0.15_0.43"]["mask"].astype(float), 512) > 0.5
    w = (bmask & d["mask"]).astype(float)
    centre = (200.0, 75.0)  # Galactic: middle of the BOSS northern cap
    for j, (m, title) in enumerate([(pred_hp, "Predicted ISW sky from BOSS (ΛCDM)"),
                                    (sky, "Planck SMICA (real sky)")]):
        ax = fig.add_subplot(gs[1, j])
        img = ortho_patch(m, centre, 55, mask=d["mask"])
        lim = 8 if j == 0 else 150
        h = ax.imshow(img, origin="lower", cmap=CMAP, vmin=-lim, vmax=lim)
        _img_axes(ax)
        ax.set_title(title, fontsize=9.5)
        _cbar(fig, h, ax, "µK (ℓ ≥ 10)")
    ax = fig.add_subplot(gs[1, 2])
    ell = np.arange(61)
    for m, col, lab in ((sky, P.C_OURS, "Planck sky"), (pred_hp, P.C_MODEL, "predicted ISW")):
        cl = hp.anafast(m * w, lmax=60) / np.mean(w**2)
        ax.plot(ell[10:], (ell * (ell + 1) * cl / (2 * np.pi))[10:], color=col, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("Multipole ℓ")
    ax.set_ylabel("ℓ(ℓ+1)C_ℓ / 2π (µK²)")
    ax.set_title("Inside BOSS, the ISW is ~1–2% of the CMB's power", fontsize=9.5)
    ax.legend(fontsize=7.5)
    fig.suptitle("From 1.2 million BOSS galaxies to a predicted ISW sky", x=0.02, ha="left")
    P.caption(fig, f"Top: a {3 * GRID_CELL:.0f} h⁻¹ Mpc thick slab {zslab:.0f} h⁻¹ Mpc north of the "
                   "equatorial plane (circles: Mao+17 voids centred within half a radius of it). The "
                   "excess power at k < 0.02 h/Mpc (window smoothing and large-scale systematics) "
                   "only reaches ℓ ≲ 10, which every statistic here excludes. Bottom: orthographic "
                   "views of the northern cap; the colour scales differ by ×20.")
    return P.save(fig, FIG / "05_boss_reconstruction")


# --------------------------------------------------------------------------
# 6. Cross-correlations
# --------------------------------------------------------------------------


def fig_crosscorr(r):
    cc = r["crosscorr"]["tracers"]
    names = [("nvss", "NVSS radio sources (z̄ ≈ 1.2)"), ("boss_0.15_0.43", "BOSS LOWZ (0.15–0.43)"),
             ("boss_0.43_0.75", "BOSS CMASS (0.43–0.75)")]
    fig, axs = plt.subplots(2, 3, figsize=(12, 7.0), gridspec_kw={"height_ratios": [1.4, 1],
                                                                  "hspace": 0.5})
    for j, (k, title) in enumerate(names):
        e = cc[k]
        ell = np.array(e["ell"])
        scale = ell * (ell + 1) / (2 * np.pi)
        ax = axs[0, j]
        ax.fill_between(ell, -2 * scale * np.array(e["cl_gT_sigma"]),
                        2 * scale * np.array(e["cl_gT_sigma"]), color=P.GRID, lw=0)
        ax.errorbar(ell, scale * np.array(e["cl_gT_data"]), yerr=scale * np.array(e["cl_gT_sigma"]),
                    fmt="o", color=P.C_OURS, ms=5, capsize=0, label="Planck × galaxies")
        ax.plot(ell, scale * np.array(e["cl_gT_model"]), color=P.C_MODEL, label="ΛCDM (CAMB)")
        ax.axhline(0, color=P.C_EDS, ls="--", lw=1.1, label="matter only")
        ax.set_title(f"{title}\nA = {e['A']:.2f} ± {e['sigma']:.2f}", fontsize=9.5)
        ax.set_xlabel("Multipole ℓ")
        if j == 0:
            ax.set_ylabel("ℓ(ℓ+1) C_ℓ^gT / 2π (µK)")
            ax.legend(fontsize=7.5, loc="upper right")
        ax = axs[1, j]
        sc = ell
        ax.errorbar(ell, sc * np.array(e["cl_gg_data"]), yerr=sc * np.array(e["cl_gg_err"]), fmt="o",
                    color=P.C_RAW, ms=4, capsize=0, label="galaxy auto-spectrum")
        ax.plot(ell, sc * np.array(e["cl_gg_model"]), color=P.C_MODEL, lw=1.2, label="fitted bias")
        ax.set_xlabel("Multipole ℓ")
        if j == 0:
            ax.set_ylabel("ℓ C_ℓ^gg")
        ax.set_title(e["bias"], fontsize=8)
        ax.legend(fontsize=7)
    fig.suptitle("Are there more galaxies where the CMB is warmer?", x=0.02, ha="left")
    P.caption(fig, "Top: cross-spectra of each galaxy map with Planck SMICA (grey band: ±2σ of "
                   "chance, from 1,000 simulated skies against the real galaxy maps). The ΛCDM curve "
                   "includes the exact effect of both masks. Bottom: each tracer's own clustering, "
                   "which fixes its bias.")
    return P.save(fig, FIG / "06_cross_correlations")


# --------------------------------------------------------------------------
# 7. The Cold Spot
# --------------------------------------------------------------------------


def fig_coldspot(r):
    import coldspot as CS

    cs = r["coldspot"]
    an = CS.Analyzer()
    filt = an.filtered(an.map128)
    centre = tuple(cs["centre"])
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), gridspec_kw={"width_ratios": [1, 1.25, 1],
                                                                 "wspace": 0.42})
    ax = axs[0]
    import cmb
    d = cmb.prepare_maps()
    img = ortho_patch(d["smica"], centre, 30, mask=d["mask"])
    # 400 px over 60 deg: sigma 5.7 px = 2 deg FWHM, enough to show a 5-deg spot
    sm = gaussian_filter(np.nan_to_num(img), 5.7)
    sm[np.isnan(img)] = np.nan
    h = ax.imshow(sm, origin="lower", cmap=CMAP, vmin=-150, vmax=150)
    ax.add_patch(Circle((200, 200), 200 * np.sin(np.radians(5)) / np.sin(np.radians(30)),
                        fill=False, color=P.INK, lw=0.8, ls="--"))
    _img_axes(ax)
    ax.set_title(f"The Cold Spot, (l, b) = ({centre[0]:.0f}°, {centre[1]:.0f}°)\n"
                 "60° across, 2° smoothing; circle 5°", fontsize=9.5)
    _cbar(fig, h, ax, "µK")

    ax = axs[1]
    pr = cs["profile"]
    edges = np.array(pr["edges"])
    mid = 0.5 * (edges[1:] + edges[:-1])
    ax.fill_between(mid, -np.array(pr["null_sigma"]), np.array(pr["null_sigma"]), color=P.GRID, lw=0,
                    label="±1σ around a random point")
    ax.plot(mid, pr["data"], color=P.C_OURS, marker="o", ms=3.5, label="Planck: Cold Spot")
    first = True
    for name, m in cs["supervoid_models"].items():
        th, p = np.array(m["theta"]), np.array(m["profile"])
        ax.plot(th, p, color=P.C_MODEL, lw=1.4 if first else 0.9, ls="-" if first else ":",
                label=f"Supervoid ISW, ΛCDM ({name})" if first else name)
        if first:
            a5 = LITERATURE["kovacs2022"]["supervoid_A_ISW"][0]
            ax.plot(th, a5 * p, color=P.C_MODEL, ls="--", lw=1.0, label=f"same × {a5} (A_ISW of other supervoids)")
        first = False
    ax.axhline(0, color=P.C_EDS, ls="--", lw=1)
    ax.set_xlim(0, 25)
    ax.set_xlabel("Angle from the Cold Spot centre (deg)")
    ax.set_ylabel("Mean temperature in rings (µK)")
    ax.set_title("Could the Eridanus supervoid explain it?", fontsize=9.5)
    ax.legend(fontsize=6.8, loc="lower right")

    ax = axs[2]
    sims = CS.run_sims()
    k5 = ("smhw", 5.0)
    kg = ("gauss", CS.GAUSS_FWHM_DEG)
    f5 = cs["filters"]["smhw_5"]
    fg = cs["filters"][f"gauss_{CS.GAUSS_FWHM_DEG:g}"]
    bins = np.linspace(-6, -2, 41)
    ax.hist(sims[k5]["min_sigma"], bins=bins, color=P.AXIS, edgecolor=P.SURFACE,
            label="1,000 simulated skies: SMHW 5°")
    ax.axvline(f5["min_sigma"], color=P.C_OURS, lw=1.6)
    ax.text(f5["min_sigma"] - 0.05, ax.get_ylim()[1] * 0.9,
            f"Cold Spot\n{f5['min_sigma']:.2f}σ\np = {f5['p_min']:.3f}", ha="right", va="top",
            color=P.C_OURS, fontsize=8)
    ax.set_xlabel("Coldest spot on the sky (σ of the filtered map)")
    ax.set_ylabel("Simulated skies")
    ax.set_title(f"Unusual? Wavelet: p = {f5['p_min']:.3f}; 10° Gaussian: p = {fg['p_min']:.2f}",
                 fontsize=9.5)
    ax.legend(fontsize=7, loc="upper left")
    P.caption(fig, "Middle: ring-averaged temperature around the Cold Spot vs the linear ISW profile "
                   "of published Eridanus supervoid models computed with D3's cosmology. Right: how "
                   "often chance produces a spot as cold anywhere on the analysed sky, with the "
                   "spherical Mexican-hat wavelet; with a Gaussian filter the Cold Spot is unremarkable.")
    return P.save(fig, FIG / "07_cold_spot")


def main():
    r = load()
    for f in (fig_summary, fig_idea, fig_granett_stack, fig_scan, fig_stacks, fig_reconstruction,
              fig_crosscorr, fig_coldspot):
        out = f(r)
        print("  wrote", ", ".join(str(o).split("/")[-1] for o in out))


if __name__ == "__main__":
    main()
