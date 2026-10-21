"""Every figure of project H4, drawn from the cached results.

01 quasar proper motions: the haystack (raw) and the error bars (normalised)
02 Gaia's error bars on quasars versus magnitude, colour and quality
03 THE KEEPER: the sky of averaged quasar drifts converging on the Galactic Centre
04 the drift towards the Galactic Centre versus angular distance from it
05 where the acceleration points, against the Galactic Centre
06 how big it is, against the Galaxy's rotation (H1, H2 and others)
07 the glide versus the number of VSH terms fitted alongside it (l_max)
08 the same acceleration from every kind of quasar
09 the large-scale patterns in Gaia's quasar proper motions (VSH power)
10 how the signal climbs out of the noise as quasars are added
"""

from __future__ import annotations

import math

import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import config as C
import plotstyle as ps
import skymap

ps.use_style()
FIG = C.FIGURE_DIR
A = np.array(C.A_G)
K21 = C.REFERENCE["klioner2021"]


def _gal_lb_of_icrs(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    u = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    v = A @ u
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(v[2]))


SGRA_LB = _gal_lb_of_icrs(C.SGRA_RA_DEG, C.SGRA_DEC_DEG)


# --------------------------------------------------------------------------
# 01 the haystack
# --------------------------------------------------------------------------
def fig01_distributions(df, errs: dict) -> None:
    five = df["five_param"].to_numpy()
    pmra = df.loc[five, "pmra"].to_numpy() * 1000
    pmdec = df.loc[five, "pmdec"].to_numpy() * 1000
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.05, 1]})
    ax = axes[0]
    edges = np.linspace(-6000, 6000, 241)
    for vals, lab, col, mk in ((pmra, "east–west (µ_α*)", ps.BLUE, "o"),
                               (pmdec, "north–south (µ_δ)", ps.ORANGE, "s")):
        h, _ = np.histogram(vals, edges)
        c = 0.5 * (edges[1:] + edges[:-1])
        ax.step(c, h, where="mid", color=col, lw=1.4, label=lab)
    ax.set_yscale("log")
    ax.set_ylim(1, None)
    ax.axvspan(-5, 5, color=ps.RED, alpha=0.9, lw=0, zorder=4)
    ps.annotate(ax, "the signal we are after:\n±5 µas/yr — the thin red line", (5, 2e3), (1900, 5e3), color=ps.INK)
    med = float(np.median(np.hypot(df.loc[five, "pmra_error"], df.loc[five, "pmdec_error"]) / math.sqrt(2))) * 1000
    ax.set_xlabel("proper motion of each quasar (µas per year)")
    ax.set_ylabel("quasars per 50 µas/yr bin")
    ax.set_title("1.2 million quasars, none of which should move")
    ax.legend(loc="upper left")
    ax.text(0.98, 0.97, f"median error\n{med:.0f} µas/yr per component", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color=ps.INK_SECONDARY)

    ax = axes[1]
    for key, lab, col, mk in (("parallax", "parallax ÷ σ", ps.AQUA, "D"),
                              ("pmra", "µ_α* / σ", ps.BLUE, "o"), ("pmdec", "µ_δ / σ", ps.ORANGE, "s")):
        e = errs[key]
        ed = np.asarray(e["hist_edges"])
        cnt = np.asarray(e["hist_counts"])
        c = 0.5 * (ed[1:] + ed[:-1])
        ax.step(c, cnt, where="mid", color=col, lw=1.3, label=f"{lab}: width {e['fit_sigma']:.3f}")
        w = ed[1] - ed[0]
        model = cnt.sum() * w * np.exp(-0.5 * ((c - e["fit_mean"]) / e["fit_sigma"]) ** 2) / (
            e["fit_sigma"] * math.sqrt(2 * math.pi))
        ax.plot(c, model, color=ps.INK, lw=0.7, ls=(0, (2, 2)))
    ax.plot([], [], color=ps.INK, lw=0.7, ls=(0, (2, 2)), label="best-fit Gaussians")
    ax.set_yscale("log")
    ax.set_ylim(1, None)
    ax.set_xlabel("measurement ÷ Gaia's stated error")
    ax.set_ylabel("quasars per 0.05 bin")
    ax.set_title("Divided by its error bar, each is noise")
    ax.legend(loc="lower center", fontsize=7.5)
    ps.caption(fig, "Gaia DR3 (= EDR3 astrometry), 1,215,942 five-parameter Gaia-CRF3 sources. Widths of the "
                    f"normalised distributions: Klioner et al. 2021 found {K21['normalised_widths']['parallax']}, "
                    f"{K21['normalised_widths']['pmra']}, {K21['normalised_widths']['pmdec']} "
                    "(parallax, µ_α*, µ_δ). Gaia's errors are ~5% too small for most quasars, more for the brightest "
                    "(Fig. 02). The parallaxes centre on −0.05σ: Gaia's −17 µas zero point.",
               y=-0.03)
    fig.tight_layout()
    ps.save(fig, FIG / "01_quasar_proper_motions")


# --------------------------------------------------------------------------
# 02 the error model
# --------------------------------------------------------------------------
def fig02_error_model(errs: dict) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.6))
    panels = (("k_vs_g", "G magnitude", "brightness"), ("k_vs_nu", "ν_eff (µm⁻¹, bluer →)", "colour"),
              ("k_vs_ruwe", "RUWE (astrometric goodness of fit)", "quality of the fit"),
              ("k_vs_ecl_lat", "|ecliptic latitude| (°)", "Gaia's scanning pattern"),
              ("k_vs_aen", "astrometric excess noise (mas)", "excess noise"),
              ("k_vs_bprp_excess", "BP/RP flux excess (1.2 = point source; more = extended)", "extended light"))
    for ax, (key, xlabel, title) in zip(axes.ravel(), panels):
        rows = errs[key]
        x = np.array([0.5 * (r["lo"] + r["hi"]) for r in rows])
        if key in ("k_vs_ruwe", "k_vs_aen", "k_vs_bprp_excess"):
            x = np.arange(len(rows))
            ax.set_xticks(x, [f"{r['lo']:g}–{r['hi']:g}" for r in rows], rotation=45, fontsize=7)
        k = np.array([r["k"] for r in rows])
        ke = np.array([r["k_err"] for r in rows])
        n = np.array([r["n"] for r in rows])
        ax.errorbar(x, k, ke, fmt="o", ms=4, color=ps.BLUE, ecolor=ps.BLUE, elinewidth=1, capsize=0)
        ax.axhline(1.0, color=ps.MUTED, lw=0.8)
        ax.axhline(errs["k_overall"], color=ps.ORANGE, lw=1.0, ls=(0, (4, 2)))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("true scatter ÷ Gaia's error")
        for xi, ni, ki in zip(x, n, k):
            if ni < 3000:
                ax.text(xi, ki + 0.02, f"{ni}", fontsize=6, color=ps.MUTED, ha="center")
    axes[0, 0].text(0.98, 0.95, f"all quasars: k = {errs['k_overall']:.3f}", transform=axes[0, 0].transAxes,
                    color=ps.ORANGE, fontsize=8, va="top", ha="right")
    fig.suptitle("Gaia's proper-motion errors for quasars: how honest are they?", x=0.01, ha="left")
    ps.caption(fig, "k = sqrt(median X² / 2 ln 2), X² = µᵀC⁻¹µ of each five-parameter quasar against zero "
                    "motion (the 5 µas/yr signal is negligible here). k = 1 means Gaia's covariance is "
                    "exactly right; the orange line is the whole sample. Small counts are printed above "
                    "sparse bins.", y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "02_error_model")


# --------------------------------------------------------------------------
# 03 the keeper
# --------------------------------------------------------------------------
def _weight_map(df, nside=32):
    five = df["five_param"].to_numpy()
    l = df.loc[five, "l"].to_numpy()
    b = df.loc[five, "b"].to_numpy()
    w = 1 / (df.loc[five, "pmra_error"].to_numpy() ** 2) + 1 / (df.loc[five, "pmdec_error"].to_numpy() ** 2)
    pix = hp.ang2pix(nside, l, b, lonlat=True)
    m = np.bincount(pix, weights=w, minlength=hp.nside2npix(nside))
    return m / hp.nside2pixarea(nside, degrees=True), nside


def _flow_lines(ax, color, lw=0.6, alpha=0.5):
    """Streamlines of a glide towards the Galactic Centre: great circles from
    the anticentre to the centre, every 20 degrees of position angle."""
    t = np.radians(np.linspace(0.5, 179.5, 240))
    for pa in range(0, 360, 20):
        a = np.radians(pa)
        # start at the anticentre, head off along position angle pa
        u = np.stack([-np.cos(t), np.sin(t) * np.cos(a), np.sin(t) * np.sin(a)], -1)
        l, b = skymap.lonlat(u)
        x, y = skymap.project(l, b)
        jump = np.abs(np.diff(x)) > 1.0
        x, y = x.astype(float), y.astype(float)
        x[1:][jump] = np.nan
        ax.plot(x, y, color=color, lw=lw, alpha=alpha, zorder=1.5)


def fig03_keeper(df, skymaps: dict, base: dict, key: str = "clean_nside4_r30") -> None:
    sm = skymaps[key]
    wmap, nside = _weight_map(df)
    fig = plt.figure(figsize=(12.5, 8.6))
    ax = fig.add_axes([0.02, 0.12, 0.96, 0.77])
    vmax = np.percentile(wmap[wmap > 0], 99)

    def lookup(l, b):
        v = wmap[hp.ang2pix(nside, l, b, lonlat=True)]
        return np.where(v > 0, v / vmax, np.nan)

    skymap.raster(ax, lookup, cmap=ps.SEQ_BLUE, vmin=-0.25, vmax=1.9, zorder=0.2)
    skymap.frame(ax)
    _flow_lines(ax, ps.ORANGE, lw=0.7, alpha=0.35)
    l, b = np.array(sm["l"]), np.array(sm["b"])
    east, north = np.array(sm["east"]), np.array(sm["north"])
    me, mn = np.array(sm["model_east"]), np.array(sm["model_north"])
    n = np.array(sm["n"], float)
    full = n >= 0.5 * np.nanmedian(n)
    scale = 0.05  # map units per uas/yr
    skymap.arrows(ax, l, b, me, mn, scale=scale, color=ps.MODEL_ARROW, width=0.0048, alpha=0.95, zorder=4)
    skymap.arrows(ax, l[~full], b[~full], east[~full], north[~full], scale=scale, color=ps.CONTEXT,
                  width=0.0018, zorder=5)
    skymap.arrows(ax, l[full], b[full], east[full], north[full], scale=scale, color=ps.DATA_ARROW,
                  width=0.0021, zorder=6)
    gx, gy = skymap.project(*SGRA_LB)
    ax.plot(gx, gy, marker="*", ms=19, color=ps.YELLOW, mec=ps.INK, mew=0.9, zorder=8)
    ax.text(gx, gy - 0.13, "Galactic Centre", ha="center", va="top", fontsize=9, color=ps.INK, zorder=8,
            fontweight="bold", bbox={"boxstyle": "round,pad=0.15", "fc": ps.SURFACE, "ec": "none", "alpha": 0.8})
    # legend: one row beneath the sky (the y-range is extended to make room)
    ax.set_ylim(-1.66, skymap.Y_MAX * 1.04)
    for x0, (col, lw, txt) in zip((-2.85, -0.92, 1.0), (
            (ps.DATA_ARROW, 1.1, "5 µas/yr: Gaia DR3 quasars\naveraged within 30°"),
            (ps.CONTEXT, 1.1, "same, where fewer than half\nthe usual quasars are seen"),
            (ps.MODEL_ARROW, 2.4, "the fitted acceleration's\npattern, averaged the same way"))):
        y = -1.57
        ax.annotate("", xy=(x0 + 5 * scale, y), xytext=(x0, y),
                    arrowprops={"arrowstyle": "-|>", "color": col, "lw": lw}, annotation_clip=False)
        ax.text(x0 + 5 * scale + 0.06, y, txt, va="center", fontsize=8.2)
    cov = np.array(sm["cov"], float)
    typical = float(np.nanmedian(np.sqrt(0.5 * (cov[:, 0] + cov[:, 2]))))
    rms = skymaps["nonglide_rms"]
    fig.text(0.02, 0.105,
             "Shown after subtracting Gaia's own large-scale patterns, fitted at the same time as the acceleration: a spin of "
             f"the frame and VSH degrees 2–3\n(rms {rms['all_nonglide']:.1f} µas/yr over the quasars, more than the "
             f"acceleration's {rms['glide']:.1f}; the raw averages are in Fig. 03b). Each arrow averages "
             f"~{int(np.nanmedian(sm['n'])):,} quasars (error ±{typical:.1f} µas/yr per\ncomponent); neighbouring "
             "arrows share quasars. Shading: Gaia's statistical weight per square degree; the Galactic plane and bulge "
             "hide most quasars.",
             fontsize=8, color=ps.INK_SECONDARY, va="top")
    amp = base["bootstrap"]["glide_icrs"]["amplitude"]
    g = base["bootstrap"]["glide_galactic"]
    fig.suptitle("Every quasar in the sky is drifting towards the centre of the Milky Way", x=0.02, ha="left",
                 y=0.975, fontsize=15)
    fig.text(0.02, 0.925,
             f"The Sun falls around the Galaxy; its changing velocity tilts the light of 1.2 million quasars by "
             f"{amp['value']:.2f} ± {amp['err']:.2f} µas per year towards (l, b) = "
             f"({skymap.wrap180(g['lon']):.1f}°, {g['lat']:.1f}°), {base['bootstrap']['glide_angle_from_sgra_deg']:.1f}° from Sgr A*.\n"
             "Faint orange curves: the paths a pure acceleration towards the centre would draw.",
             fontsize=9.5, color=ps.INK_SECONDARY, va="top")
    ps.save(fig, FIG / "03_quasars_drift_to_galactic_centre")


def fig03b_raw_vs_clean(skymaps: dict, base: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.6))
    rms = skymaps["nonglide_rms"]
    for ax, key, title in ((axes[0], "raw_nside4_r30", "What Gaia measured (frame spin and all)"),
                           (axes[1], "clean_nside4_r30", "After subtracting the fitted spin and VSH degrees 2–3")):
        sm = skymaps[key]
        skymap.frame(ax, fontsize=6)
        _flow_lines(ax, ps.ORANGE, lw=0.5, alpha=0.3)
        scale = 0.045
        skymap.arrows(ax, sm["l"], sm["b"], sm["model_east"], sm["model_north"], scale=scale,
                      color=ps.MODEL_ARROW, width=0.0045, alpha=0.9)
        skymap.arrows(ax, sm["l"], sm["b"], sm["east"], sm["north"], scale=scale, color=ps.DATA_ARROW,
                      width=0.0022)
        gx, gy = skymap.project(*SGRA_LB)
        ax.plot(gx, gy, "*", ms=14, color=ps.YELLOW, mec=ps.INK, mew=0.8, zorder=8)
        ax.set_title(title, fontsize=10)
    fig.suptitle("Gaia's own large-scale patterns are bigger than the signal; the fit separates them by shape",
                 x=0.01, ha="left", y=1.02)
    ps.caption(fig, f"Arrows: quasars averaged within 30° (black) and the fitted acceleration pattern averaged the "
                    f"same way (orange). RMS over the quasars of the fitted acceleration (glide): "
                    f"{rms['glide']:.1f} µas/yr; of the frame spin: {rms['rotation']:.1f}; of the degree 2–3 patterns: "
                    f"{rms['l2_l3']:.1f}. Near the Galactic Centre, where few quasars are seen through the bulge, "
                    "the raw averages swirl at ~10 µas/yr.", y=0.0)
    fig.tight_layout()
    ps.save(fig, FIG / "03b_raw_versus_cleaned")


# --------------------------------------------------------------------------
# 04 the drift profile
# --------------------------------------------------------------------------
def fig04_profile(profiles: dict, base: dict) -> None:
    clean, raw = profiles["clean"], profiles["raw_minus_rotation"]
    gpar = clean["glide_towards_gc"]

    def arr(p, key):
        return np.array([r[key] for r in p["rings"]])

    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    t = np.linspace(0, 180, 300)
    ax = axes[0]
    ax.plot(t, gpar * np.sin(np.radians(t)), color=ps.ORANGE, lw=2.2,
            label=f"acceleration towards the centre: {gpar:.2f} µas/yr × sin θ")
    th = arr(clean, "theta_mean")
    ax.errorbar(arr(raw, "theta_mean") + 1.2, arr(raw, "along"), arr(raw, "along_err"), fmt="s", ms=4,
                color=ps.MUTED, mfc=ps.SURFACE, elinewidth=1, capsize=0,
                label="raw (only the frame's spin removed)")
    ax.errorbar(th, arr(clean, "along"), arr(clean, "along_err"), fmt="o", color=ps.BLUE, ms=5.5,
                elinewidth=1.3, capsize=0, label="after removing the fitted spin and degree 2–3 patterns")
    ax.axhline(0, color=ps.MUTED, lw=0.8)
    ax.set_ylabel("drift towards the\nGalactic Centre (µas/yr)")
    ax.set_title("The pull: quasars 90° from the centre drift fastest, those near it hardly at all")
    ax.legend(loc="upper right", fontsize=7.8)
    ax.set_ylim(-6, 14)
    ax = axes[1]
    ax.errorbar(th, arr(clean, "across"), arr(clean, "across_err"), fmt="o", color=ps.BLUE, ms=4.5,
                elinewidth=1.1, capsize=0)
    ax.axhline(0, color=ps.ORANGE, lw=2.2)
    ax.set_ylabel("sideways drift\n(µas/yr)")
    ax.set_xlabel("angular distance θ from the Galactic Centre (degrees)")
    ax.set_xticks(range(0, 181, 30))
    ax.set_ylim(-15, 15)
    al, ae = arr(clean, "along"), arr(clean, "along_err")
    chi_al = float(np.sum(((al - gpar * np.sin(np.radians(th))) / ae) ** 2))
    chi_ac = float(np.sum((arr(clean, "across") / arr(clean, "across_err")) ** 2))
    ax.text(0.01, 0.94, f"expected 0; χ² = {chi_ac:.1f} for {len(th)} rings", transform=ax.transAxes,
            fontsize=8, va="top", color=ps.INK_SECONDARY)
    axes[0].text(0.01, 0.04, f"blue points vs the curve: χ² = {chi_al:.1f} for {len(th)} rings",
                 transform=axes[0].transAxes, fontsize=8, color=ps.INK_SECONDARY)
    ps.caption(fig, "Five-parameter Gaia-CRF3 quasars minus the baseline fit's outliers, in 10° rings around "
                    "Sgr A*: weighted least-squares mean of the proper-motion component along (top) and across "
                    "(bottom) the great circle through the centre. No quasars are seen within 10° of it. The curve "
                    "is the glide's component towards Sgr A* from the full fit.", y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "04_drift_profile")


# --------------------------------------------------------------------------
# 05 the direction
# --------------------------------------------------------------------------
def fig05_direction(base: dict, boot_gal: np.ndarray, preds: dict, lmax_rows: list | None = None) -> None:
    gal = base["bootstrap"]["glide_galactic"]
    lb = np.array([skymap.lonlat(v / np.linalg.norm(v)) for v in boot_gal])
    l_b = skymap.wrap180(lb[:, 0])
    fig, ax = plt.subplots(figsize=(7.8, 6.8))
    from scipy.ndimage import gaussian_filter
    H, xe, ye = np.histogram2d(l_b, lb[:, 1], bins=100, range=[[-25, 25], [-25, 25]])
    H = gaussian_filter(H, 2.5)
    flat = np.sort(H.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    levels = [flat[np.searchsorted(cum, q)] for q in (0.95, 0.68)]
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    ax.contourf(xc, yc, H.T, levels=[levels[0], levels[1], H.max() * 1.01], colors=[ps.BLUE_RAMP[0], ps.BLUE_RAMP[2]],
                zorder=0)
    ax.contour(xc, yc, H.T, levels=levels, colors=[ps.BLUE], linewidths=[0.8, 1.2], zorder=1)
    kl, kb = K21["l_b"]
    ax.errorbar(skymap.wrap180(kl), kb, xerr=K21["l_b_err"][0], yerr=K21["l_b_err"][1], fmt="o", mfc="none",
                mec=ps.INK, ecolor=ps.INK, ms=13, mew=1.2, capsize=0, elinewidth=0.9,
                label=f"Klioner et al. 2021: ({kl - 360:.1f}°, {kb}°)", zorder=6)
    ax.plot(skymap.wrap180(gal["lon"]), gal["lat"], "o", color=ps.BLUE, ms=7, mec=ps.SURFACE, mew=1.0,
            label=f"this project: ({skymap.wrap180(gal['lon']):.1f}°, {gal['lat']:.1f}°)", zorder=7)
    for r in lmax_rows or []:
        if r["label"] == "1":
            ax.plot(skymap.wrap180(r["l"]), r["b"], "s", color=ps.MUTED, ms=6, zorder=6,
                    label=f"this project, dipole only (l_max = 1): ({skymap.wrap180(r['l']):.1f}°, {r['b']:.1f}°)")
    dl, db = _gal_lb_of_icrs(*C.REFERENCE["darling2025"]["ra_dec"])
    ax.plot(skymap.wrap180(dl), db, "s", color=ps.SURFACE, mec=ps.MUTED, mew=1.2, ms=6, zorder=6,
            label=f"Darling 2025, Quaia, dipole only: ({skymap.wrap180(dl):.1f}°, {db:.1f}°)")
    ch = C.REFERENCE["charlot2020"]
    cl, cb = _gal_lb_of_icrs(*ch["ra_dec"])
    ax.plot(skymap.wrap180(cl), cb, "D", color=ps.AQUA, ms=7, mec=ps.INK, mew=0.6,
            label=f"VLBI, ICRF3 (Charlot et al. 2020): ({skymap.wrap180(cl):.1f}°, {cb:.1f}°)", zorder=6)
    ax.plot(skymap.wrap180(SGRA_LB[0]), SGRA_LB[1], "*", color=ps.YELLOW, ms=18, mec=ps.INK, mew=0.9,
            label="Sgr A*, the Galactic Centre", zorder=7)
    lo, hi = K21["expected_with_lmc_range_kms_myr"]
    for vec, mk in ((K21["expected_kms_myr"], "P"), (lo, "^"), (hi, "v")):
        v = np.array(vec, float)
        pl, pb = skymap.lonlat(v / np.linalg.norm(v))
        ax.plot(skymap.wrap180(pl), pb, mk, color=ps.ORANGE, ms=7, mec=ps.INK, mew=0.5, zorder=6)
    ax.plot([], [], "P", color=ps.ORANGE, ms=7, mec=ps.INK, mew=0.5,
            label="expected: Galaxy model (+), with the LMC's pull (▲▼ range)")
    ax.set_xlim(25, -25)
    ax.set_ylim(-25, 25)
    ax.set_xlabel("galactic longitude l (°)")
    ax.set_ylabel("galactic latitude b (°)")
    ax.set_title("Where the acceleration points")
    ax.legend(loc="upper left", fontsize=7.5)
    ps.caption(fig, f"Shading: where 68% (dark) and 95% (light) of {len(boot_gal)} bootstrap resamplings of the "
                    "quasars point (each refitted with clipping). Klioner's error bars are their published 1σ. Expected "
                    "directions from Klioner et al. 2021 Sect. 3 (centripetal + bar/spiral + disc; LMC 1–2.5 × 10¹¹ M☉).",
               y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "05_direction")


# --------------------------------------------------------------------------
# 06 the amplitude against the Galaxy
# --------------------------------------------------------------------------
def fig06_amplitude(base: dict, preds: dict, implied: dict) -> None:
    b = base["bootstrap"]
    rows = []
    rows.append(("this project: towards the Galactic Centre (g_X)", b["glide_galactic"]["vec"][0],
                 b["glide_galactic"]["err"][0], "ours"))
    amp = b["glide_icrs"]["amplitude"]
    rows.append(("this project: total |g| (bias-corrected)", amp["corrected"], amp["err"], "ours"))
    rows.append(("Gaia EDR3, Klioner et al. 2021: |g|", *K21["amplitude_uas_yr"], "published"))
    rows.append(("Gaia EDR3, Klioner et al. 2021: g_X", K21["g_gal"][0], K21["g_gal_err"][0], "published"))
    for key, lab in (("tsigkas_kouvelis2026", "same data, Tsigkas Kouvelis et al. 2026: |g|"),
                     ("nguyen_murray2026", "same data, Nguyen & Murray 2026: |g|"),
                     ("darling2025", "Quaia quasars, dipole only, Darling 2025: |g|")):
        rows.append((lab, *C.REFERENCE[key]["amplitude_uas_yr"], "published"))
    ch = C.REFERENCE["charlot2020"]
    rows.append(("VLBI, ICRF3 data (Charlot et al. 2020)", *ch["amplitude_uas_yr"], "vlbi"))
    tk = C.REFERENCE["titov_krasna2018"]
    rows.append(("VLBI, Titov & Krásná 2018", *tk["amplitude_uas_yr"], "vlbi"))
    order = ["h1_sgra", "h1_h2", "h2_reid", "gravity2019", "gravity2022", "eilers2019"]
    for key in order:
        p = preds["constants"][key]
        rows.append((f"V_c²/R₀: {p['label']}", p["g_uas_yr"]["median"], p["g_uas_yr"]["sd"], "predicted"))
    mc = preds["mcmillan17"]
    rows.append(("McMillan 2017 Milky Way model (R₀ = 8.21 kpc)", mc["z_25pc"]["a_R_uas_yr"], 0.0, "predicted"))
    style = {"ours": (ps.OURS, "o", "full"), "published": (ps.PUBLISHED, "o", "none"),
             "vlbi": (ps.VLBI, "D", "full"), "predicted": (ps.PREDICTED, "s", "full")}
    fig, ax = plt.subplots(figsize=(10, 6.2))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, v, e, kind) in zip(y, rows):
        col, mk, fill = style[kind]
        ax.errorbar(v, yi, xerr=e if e else None, fmt=mk, color=col, ms=7, mfc=col if fill == "full" else ps.SURFACE,
                    mec=col, elinewidth=1.4, capsize=0)
        ax.text(3.02, yi, lab, va="center", fontsize=8.2, color=ps.INK)
    gx, gxe = b["glide_galactic"]["vec"][0], b["glide_galactic"]["err"][0]
    ax.axvspan(gx - gxe, gx + gxe, color=ps.BLUE, alpha=0.12, lw=0)
    ax.set_yticks([])
    ax.set_xlim(3.0, 7.0)
    ax.set_xlabel("drift of the quasars = the Sun's acceleration ÷ c (µas per year)")
    secax = ax.secondary_xaxis("top", functions=(lambda u: u * C.UAS_YR_TO_M_S2 * 1e10,
                                                  lambda a: a / (C.UAS_YR_TO_M_S2 * 1e10)))
    secax.set_xlabel("acceleration (10⁻¹⁰ m/s²)")
    secax.tick_params(colors=ps.MUTED)
    ax.grid(axis="y", visible=False)
    vc = implied["vc_kms_with_h1_r0"]
    ax.set_title("How hard the Galaxy pulls on the Sun", pad=28)
    ps.caption(fig, "Predictions are V_c²/R₀ with V_c = Ω☉R₀ − V☉ (Sgr A*'s proper motion, Reid & Brunthaler 2020; "
                    "V☉ = 12.24 km/s, Schönrich et al. 2010), or the constants named; they leave out the LMC's "
                    f"0.1–0.3 µas/yr pull. Our g_X with H1's R₀ gives V_c = {vc['median']:.0f} ± {vc['sd']:.0f} km/s.",
               y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "06_amplitude_vs_galaxy")


# --------------------------------------------------------------------------
# 07 l_max scan
# --------------------------------------------------------------------------
def fig07_lmax(lmax_rows: list) -> None:
    labels = [r["label"] for r in lmax_rows]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=False)
    for ax, frame, comps, ref, referr in ((axes[0], "icrs", ("x", "y", "z"), K21["g_icrs"], K21["g_icrs_err"]),
                                          (axes[1], "galactic", ("X", "Y", "Z"), K21["g_gal"], K21["g_gal_err"])):
        for i, (comp, col, mk) in enumerate(zip(comps, (ps.BLUE, ps.ORANGE, ps.AQUA), ("o", "s", "D"))):
            v = np.array([r[f"glide_{frame}"][i] for r in lmax_rows])
            e = np.array([r[f"glide_{frame}_err"][i] for r in lmax_rows])
            ax.errorbar(x + (i - 1) * 0.18, v, e, fmt=mk, color=col, ms=5, capsize=0, label=f"g_{comp}")
            ax.axhspan(ref[i] - referr[i], ref[i] + referr[i], color=col, alpha=0.12, lw=0)
        ax.set_xticks(x, labels)
        ax.set_xlabel("highest VSH degree fitted with the glide (l_max)")
        ax.set_ylabel("glide component (µas/yr)")
        ax.axhline(0, color=ps.MUTED, lw=0.8)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncols=3, fontsize=8)
        ax.set_title("equatorial (ICRS)" if frame == "icrs" else "galactic: X → centre, Y → rotation, Z → north")
    fig.suptitle("The glide is stable once the large-scale systematics are fitted alongside it", x=0.01, ha="left")
    ps.caption(fig, "Bands: Klioner et al. 2021 Table 2 (l_max = 3, bootstrap 1σ). 'alone' = the three glide "
                    "components with no other terms. Error bars: fit errors × unit weight error.", y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "07_lmax_scan")


# --------------------------------------------------------------------------
# 08 subsets
# --------------------------------------------------------------------------
def _pretty(label: str) -> str:
    rep = {"G -inf-18.0": "G ≤ 18", "G 18.0-19.0": "G 18–19", "G 19.0-20.0": "G 19–20", "G 20.0-inf": "G > 20",
           "kappa=inf": "no clipping (κ = ∞)", "errors x k(G, nu_eff)": "errors × k(G, ν_eff) (Fig. 02)",
           "conditioned, zero point -17 uas": "conditioned, zero point −17 µas",
           "conditioned, Lindegren+21 zero point": "conditioned, Lindegren+21 zero point"}
    if label in rep:
        return rep[label]
    if label.startswith("kappa="):
        return "clip at κ = " + label.split("=")[1]
    if label.startswith("nu_eff "):
        lo, hi = label.split(" ")[1].split("-")
        return f"ν_eff {float(lo):.2f}–{float(hi):.2f}"
    if label.startswith("z_QSOC "):
        lo, hi = label.split(" ")[1].split("-")
        return f"Gaia redshift {lo}–{hi}"
    return label


def fig08_subsets(subsets: list, choices: list, dr16q: list, contam: dict, errfit: dict, base: dict,
                  cond: list) -> None:
    groups: dict = {"analysis choices": []}
    for r in choices:
        groups["analysis choices"].append(r)
    groups["analysis choices"].append(errfit)
    for r in cond:
        if "glide_galactic" in r:
            groups["analysis choices"].append(r)
    names = {"magnitude": "brightness", "colour": "colour", "redshift (Gaia)": "redshift (Gaia's own, flag 0)",
             "catalogue": "which catalogue found it"}
    for r in subsets:
        groups.setdefault(names.get(r["group"], r["group"]), []).append(r)
    # layout: one row per fit, one blank row (holding the header) per group
    rows, headers, y = [], [], 0.0
    for gi, (g, rs) in enumerate(groups.items()):
        headers.append((g, y))
        y -= 1.0
        for r in rs:
            rows.append((gi, r, y))
            y -= 1.0
        y -= 0.4
    fig, axes = plt.subplots(1, 3, figsize=(13, 0.25 * (-y) + 1.8), sharey=True)
    b = base["bootstrap"]["glide_galactic"]
    cols = [ps.INK_SECONDARY, ps.ORANGE, ps.AQUA, ps.BLUE, ps.VIOLET]
    for i, (ax, comp) in enumerate(zip(axes, ("X (towards the centre)", "Y (direction of rotation)", "Z (north)"))):
        for gi, r, yy in rows:
            ax.errorbar(r["glide_galactic"][i], yy, xerr=r["glide_galactic_err"][i], fmt="o",
                        color=cols[gi % len(cols)], ms=4, elinewidth=1.1, capsize=0)
        ax.axvspan(b["vec"][i] - b["err"][i], b["vec"][i] + b["err"][i], color=ps.BLUE, alpha=0.12, lw=0)
        ax.axvline(b["vec"][i], color=ps.BLUE, lw=0.8)
        ax.axvline(0, color=ps.MUTED, lw=0.8)
        ax.set_title(f"g_{comp}", fontsize=10)
        ax.set_xlabel("µas/yr")
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks([yy for _, _, yy in rows],
                       [f"{_pretty(r['label'])}  [{r.get('n_sample', r['n_used']):,}]" for _, r, _ in rows], fontsize=7)
    for g, yy in headers:
        axes[0].text(-0.02, yy, g, transform=axes[0].get_yaxis_transform(), fontsize=8, fontweight="bold",
                     color=ps.INK, ha="right", va="center")
    axes[0].set_ylim(y + 0.6, 0.6)
    fig.suptitle("Slicing the quasars and changing the analysis: mostly the same answer (redshift: Fig. 11)",
                 x=0.01, ha="left")
    ps.caption(fig, "Each row is a separate VSH fit (l_max = 3, κ = 3; l_max = 1 for the SDSS-only catalogue and the "
                    "small remainder), errors = fit errors × unit weight error. Band: the baseline (bootstrap 1σ). "
                    "Brackets: quasars used or in the subset. Outliers: the naive zero-parallax conditioning and the "
                    "six-parameter sources (README §6).", y=-0.005)
    fig.tight_layout()
    ps.save(fig, FIG / "08_subsets")


# --------------------------------------------------------------------------
# 09 VSH power
# --------------------------------------------------------------------------
def fig09_power(lmax_rows: list) -> None:
    top = max((r for r in lmax_rows if r["label"] != "alone"), key=lambda r: r["lmax"])
    l = np.arange(1, top["lmax"] + 1)
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    for kind, col, mk, lab in (("S", ps.BLUE, "o", "spheroidal (flows: glide at l = 1)"),
                               ("T", ps.ORANGE, "s", "toroidal (swirls: rotation at l = 1)")):
        p = np.array(top["power"][kind])
        nz = np.array(top["noise_power"][kind])
        ax.plot(l + (0.08 if kind == "T" else -0.08), p, mk, color=col, ms=7, label=lab)
        ax.plot(l, nz, color=col, lw=1, ls=(0, (3, 2)))
    ax.plot([], [], color=ps.INK_SECONDARY, lw=1, ls=(0, (3, 2)), label="expected from noise alone")
    ax.set_yscale("log")
    ax.set_yticks([0.5, 1, 2, 3, 4, 6], ["0.5", "1", "2", "3", "4", "6"])
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_ylim(0.5, 8)
    ax.set_xticks(l)
    ax.annotate("the acceleration (glide)", (1 - 0.08, top["power"]["S"][0]), (1.4, 7.0), fontsize=8,
                color=ps.INK_SECONDARY, arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.7})
    ax.annotate("spin of Gaia's frame", (1 + 0.08, top["power"]["T"][0]), (1.5, 1.6), fontsize=8,
                color=ps.INK_SECONDARY, arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.7})
    ax.set_xlabel("VSH degree l (pattern size ~180°/l)")
    ax.set_ylabel("RMS proper motion over the sky (µas/yr)")
    ax.set_title(f"Patterns in Gaia's quasar proper motions (fit with l_max = {top['lmax']})")
    ax.legend(loc="upper right")
    ps.caption(fig, "Power per degree = sum of squared orthonormal VSH coefficients (Mignard & Klioner 2012), shown "
                    "as the RMS proper motion it contributes over the whole sphere. Noise level: the same from the "
                    "fit's covariance. Points far above their dashed line are real patterns.", y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "09_vsh_power")


# --------------------------------------------------------------------------
# 10 emergence
# --------------------------------------------------------------------------
def fig10_emergence(emerge: list, base: dict) -> None:
    n = np.array([r["n_sample"] for r in emerge])
    gx = np.array([r["glide_galactic"][0] for r in emerge])
    gxe = np.array([r["glide_galactic_err"][0] for r in emerge])
    ang = np.array([r["angle_from_sgra_deg"] for r in emerge])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax = axes[0]
    rng = np.random.default_rng(0)
    jitter = np.exp(rng.uniform(-0.08, 0.08, n.size))
    ax.errorbar(n * jitter, gx, gxe, fmt="o", color=ps.BLUE, ms=4, elinewidth=1, capsize=0, alpha=0.9)
    b = base["bootstrap"]["glide_galactic"]
    ax.axhline(b["vec"][0], color=ps.BLUE, lw=0.8)
    ax.axhline(0, color=ps.MUTED, lw=0.8)
    ns = np.logspace(3, np.log10(1.25e6), 100)
    sig = b["err"][0] * np.sqrt(base["n_sample"] / ns)
    ax.fill_between(ns, -sig, sig, color=ps.MUTED, alpha=0.15, lw=0, label="±1σ of pure noise")
    for k, ls, ha, f in ((3, (0, (3, 2)), "right", 0.92), (5, "-", "left", 1.08)):
        nk = base["n_sample"] * (k * b["err"][0] / b["vec"][0]) ** 2
        ax.axvline(nk, color=ps.ORANGE, lw=1, ls=ls)
        ax.text(nk * f, -33, f"{k}σ detection\n~{nk / 1000:.0f}k quasars", fontsize=7.5, color=ps.ORANGE,
                va="bottom", ha=ha)
    ax.set_xscale("log")
    ax.set_ylim(-40, 40)
    ax.set_xlabel("number of quasars used (random draws)")
    ax.set_ylabel("g_X, drift towards the centre (µas/yr)")
    ax.set_title("The signal climbs out of the noise as 1/√N")
    ax.legend(loc="upper right")
    ax = axes[1]
    ax.plot(n * jitter, ang, "o", color=ps.BLUE, ms=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1, 180)
    ax.axhline(90, color=ps.MUTED, lw=0.8)
    ax.text(1.2e3, 100, "median for a random direction: 90°", fontsize=7.5, color=ps.MUTED)
    ax.set_xlabel("number of quasars used")
    ax.set_ylabel("angle between fitted acceleration and Sgr A* (°)")
    ax.set_title("…and swings round to point at the Galactic Centre")
    ps.caption(fig, "Five random draws at each size from the five-parameter sample, each fitted like the baseline "
                    "(l_max = 3, κ = 3). Detection thresholds scale the full-sample error by √N.", y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "10_signal_emergence")


def make_all(df, results, skymaps, emerge, boot_gal, errs) -> None:
    base = results["baseline"]
    fig01_distributions(df, errs)
    fig02_error_model(errs)
    fig03_keeper(df, skymaps, base)
    fig03b_raw_vs_clean(skymaps, base)
    fig04_profile(results["drift_profile"], base)
    fig05_direction(base, boot_gal, results["predictions"], results["lmax_scan"])
    fig06_amplitude(base, results["predictions"], results["implications"])
    fig07_lmax(results["lmax_scan"])
    fig08_subsets(results["subsets"], results["choices"], results["dr16q"], results["contamination"],
                  results["error_model_refit"], base, results["conditioning"])
    fig09_power(results["lmax_scan"])
    fig10_emergence(emerge, base)
    fig11_redshift(results["redshift_quaia"], results["redshift_null"], results["dr16q"], base)
    fig12_frame_and_purity(results["frame_rotation"], results["lookalikes"], base)


# --------------------------------------------------------------------------
# 11 does the drift depend on distance?
# --------------------------------------------------------------------------
def fig11_redshift(redshift: dict, null: dict, dr16q: list, base: dict) -> None:
    b = base["bootstrap"]["glide_galactic"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.9), gridspec_kw={"width_ratios": [1.35, 1.05, 0.9]})
    # A: the drift in redshift bins
    ax = axes[0]
    ax.axhspan(b["vec"][0] - b["err"][0], b["vec"][0] + b["err"][0], color=ps.BLUE, alpha=0.12, lw=0)
    ax.axhline(b["vec"][0], color=ps.BLUE, lw=0.8)
    series = (("raw", ps.BLUE, "o", 0.97, "Quaia redshifts, all sky"),
              ("conditioned_l21", ps.ORANGE, "s", 1.03, "same, proper motions conditioned on zero parallax"))
    for variant, col, mk, off, lab in series:
        rows = redshift[variant]["bins"]
        z = np.array([r["z_median"] for r in rows])
        v = np.array([r["glide_galactic"][0] for r in rows])
        e = np.array([r["glide_galactic_err"][0] for r in rows])
        cons = redshift[variant]["consistency_gX"]
        ax.errorbar(z * off, v, e, fmt=mk, color=col, ms=6, capsize=0, elinewidth=1.3,
                    label=f"{lab} (χ² = {cons['chi2']:.0f} for {cons['dof']} dof)")
    rows = [r for r in dr16q if r["label"].startswith("z_SDSS") and "removed" in r["group"]]
    ax.errorbar([r["z_median"] for r in rows], [r["glide_galactic"][0] for r in rows],
                [r["glide_galactic_err"][0] for r in rows], fmt="D", color=ps.AQUA, ms=5, capsize=0,
                elinewidth=1.1, mec=ps.INK, mew=0.4, label="SDSS spectra (SDSS sky only; all-sky patterns removed)")
    ax.axhline(0, color=ps.MUTED, lw=0.8)
    ax.axvline(1.25, color=ps.MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.set_xlabel("redshift (median of each bin)")
    ax.set_ylabel("g_X, drift towards the centre (µas/yr)")
    ax.set_title("The drift should not depend on distance — but here it does", fontsize=10)
    ax.legend(loc="upper left", fontsize=7)
    ax.set_ylim(-6, 13.5)
    # B: the split at z = 1.25 under every variant
    ax = axes[1]
    variants = [("l_max = 3 (baseline)", null["real"]), ("l_max = 6", null["real_lmax6"]),
                ("conditioned on zero parallax", null["real_conditioned_l21"])]
    mag = redshift["z_by_magnitude"]
    for glo, ghi in (("0", "19.5"), ("19.5", "20.0"), ("20.0", "20.5")):
        lo = next(r for r in mag if r["label"] == f"z 0.0-1.25, G {glo}-{ghi}")
        hi = next(r for r in mag if r["label"] == f"z 1.25-5.0, G {glo}-{ghi}")
        variants.append((f"G {glo if glo != '0' else '<'}{'–' + ghi if glo != '0' else ghi}",
                         {"low": lo["glide_galactic"][0], "low_err": lo["glide_galactic_err"][0],
                          "high": hi["glide_galactic"][0], "high_err": hi["glide_galactic_err"][0]}))
    y = np.arange(len(variants))[::-1]
    for yi, (lab, d) in zip(y, variants):
        ax.errorbar(d["low"], yi + 0.12, xerr=d["low_err"], fmt="o", color=ps.BLUE, ms=5, capsize=0)
        ax.errorbar(d["high"], yi - 0.12, xerr=d["high_err"], fmt="s", color=ps.ORANGE, ms=5, capsize=0)
    ax.plot([], [], "o", color=ps.BLUE, label="z < 1.25")
    ax.plot([], [], "s", color=ps.ORANGE, label="z ≥ 1.25")
    ax.set_yticks(y, [v[0] for v in variants], fontsize=7.5)
    ax.axvline(0, color=ps.MUTED, lw=0.8)
    ax.axvspan(b["vec"][0] - b["err"][0], b["vec"][0] + b["err"][0], color=ps.BLUE, alpha=0.12, lw=0)
    ax.set_xlabel("g_X (µas/yr)")
    ax.set_title("Nearer vs farther quasars, every way", fontsize=10)
    ax.legend(loc="lower right", fontsize=7.5)
    ax.grid(axis="y", visible=False)
    # C: the null test
    ax = axes[2]
    d = np.array([s["diff"] for s in null["shuffled"]])
    ax.hist(d, bins=np.linspace(-3, 3, 25), color=ps.CONTEXT, edgecolor=ps.SURFACE)
    real = null["real"]["diff"]
    ax.axvline(real, color=ps.RED, lw=2)
    ax.text(real, ax.get_ylim()[1] * 0.92, f" real: {real:+.2f} ± {null['real']['diff_err']:.2f}",
            color=ps.INK, fontsize=8, ha="right" if real > 0 else "left")
    ax.set_xlabel("g_X(z ≥ 1.25) − g_X(z < 1.25)  (µas/yr)")
    ax.set_ylabel("shuffles")
    ax.set_title(f"Redshifts shuffled within 7° sky pixels ({len(d)}×)", fontsize=10)
    ax.text(0.02, 0.97, f"shuffled spread: ±{null['shuffled_diff_sd']:.2f}\nreal = {null['real_over_shuffled_sd']:.1f}× that",
            transform=ax.transAxes, va="top", fontsize=8, color=ps.INK_SECONDARY)
    ps.caption(fig, "Each point is a separate VSH fit of the five-parameter Gaia-CRF3 quasars with Quaia redshifts "
                    "(l_max = 3, κ = 3 unless noted; errors = fit errors × unit weight error). Shuffling keeps each "
                    "redshift bin's sky coverage but destroys its redshifts. Tsigkas Kouvelis et al. (2026) found the "
                    "same split at 2–3σ. A real acceleration cannot depend on distance: this is a systematic effect "
                    "that tracks redshift.", y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "11_redshift")


# --------------------------------------------------------------------------
# 12 Gaia's frame, and why the sample must be pure
# --------------------------------------------------------------------------
def fig12_frame_and_purity(frame: dict, look: dict, base: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.9), gridspec_kw={"width_ratios": [1.25, 0.8, 1.0]})
    # A: the frame's spin, by group and magnitude
    ax = axes[0]
    rows = [r for r in frame["fits"] if r["lmax"] == 1 and r["label"] != "frame-rotator sources (5p + 6p)"]
    names = [r["label"] for r in rows]
    y = np.arange(len(rows))[::-1]
    for i, (comp, col, mk) in enumerate(zip("xyz", (ps.BLUE, ps.ORANGE, ps.AQUA), ("o", "s", "D"))):
        v = np.array([r["rotation_icrs"][i] for r in rows])
        e = np.array([r["rotation_icrs_err"][i] for r in rows])
        ax.errorbar(v, y + (1 - i) * 0.2, xerr=e, fmt=mk, color=col, ms=5, capsize=0, label=f"ω_{comp}")
    ax.axvline(0, color=ps.MUTED, lw=0.8)
    ax.set_yticks(y, [f"{n}\n[{r['n_sample']:,}]" for n, r in zip(names, rows)], fontsize=7.5)
    ax.set_xlabel("spin of the quasars relative to Gaia's frame (µas/yr, ICRS axes)")
    ax.set_title("Gaia's frame was made not to spin — for its 428,034 quasars", fontsize=10)
    ax.legend(loc="lower right", ncols=3, fontsize=7.5)
    ax.grid(axis="y", visible=False)
    # B: look-alikes versus quasars
    ax = axes[1]
    lk = look["fit_lmax3"]
    q_amp = base["bootstrap"]["glide_icrs"]["amplitude"]["value"]
    vals = [q_amp, lk["amplitude"]["value"]]
    errs = [base["bootstrap"]["glide_icrs"]["amplitude"]["err"], lk["amplitude"]["err"]]
    ax.bar([0, 1], vals, color=[ps.BLUE, ps.ORANGE], width=0.6, yerr=errs, error_kw={"ecolor": ps.INK, "lw": 1})
    ax.set_yscale("log")
    ax.set_ylim(1, max(vals) * 4)
    ax.set_xticks([0, 1], [f"Gaia-CRF3 quasars\n(l, b) = ({skymap.wrap180(base['bootstrap']['glide_galactic']['lon']):.0f}°, "
                           f"{base['bootstrap']['glide_galactic']['lat']:.0f}°)",
                           f"look-alikes (mostly stars)\n(l, b) = ({lk['l']:.0f}°, {lk['b']:.0f}°)"], fontsize=7.5)
    for x, v in zip([0, 1], vals):
        ax.text(x, v * 1.35, f"{v:.1f}", ha="center", fontsize=9, color=ps.INK)
    ax.set_ylabel("fitted glide |g| (µas/yr)")
    ax.set_title(f"Same cuts, not quasars: {look['n']:,} sources", fontsize=10)
    ax.grid(axis="x", visible=False)
    # C: mixing look-alikes into the quasars
    ax = axes[2]
    mix = look["mixing"]
    f = np.array([r["fraction"] for r in mix]) * 100
    b = base["bootstrap"]["glide_galactic"]
    for i, (comp, col, mk) in enumerate(zip(("X", "Y", "Z"), (ps.BLUE, ps.ORANGE, ps.AQUA), ("o", "s", "D"))):
        v = np.array([r["glide_galactic"][i] for r in mix])
        e = np.array([r["glide_galactic_err"][i] for r in mix])
        ax.errorbar(f * (1 + (i - 1) * 0.06), v, e, fmt=mk, color=col, ms=5, capsize=0, label=f"g_{comp}",
                    ls="-", lw=1)
        ax.axhspan(b["vec"][i] - b["err"][i], b["vec"][i] + b["err"][i], color=col, alpha=0.10, lw=0)
    ax.set_xscale("log")
    ax.axhline(0, color=ps.MUTED, lw=0.8)
    ax.set_xlabel("look-alikes mixed into the quasars (%)")
    ax.set_ylabel("fitted glide component (µas/yr)")
    ax.set_title("A 1% stellar contamination would show", fontsize=10)
    ax.legend(loc="lower left", fontsize=7.5)
    ps.caption(fig, "Left: l = 1 rotation (glide + rotation fit, κ = 3); frame-rotator sources from "
                    "gaiadr3.frame_rotator_source (used_for_reference_frame_spin). Middle and right: a random 0.5% of "
                    "Gaia DR3 five-parameter sources with G > 16, parallax and proper motion within 5σ of zero, not in "
                    "the Gaia-CRF3; bands are the baseline result.", y=-0.01)
    fig.tight_layout()
    ps.save(fig, FIG / "12_frame_and_purity")
