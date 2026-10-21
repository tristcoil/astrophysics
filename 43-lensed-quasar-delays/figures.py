"""Every figure of the project, from the cached results.

    python figures.py            # all figures
    python figures.py flicker    # only those whose name contains "flicker"

Inputs: ``results/delays/*.json`` (``measure_delays.py``), ``results/h0.json``
and ``results/h0_posteriors.npz`` (``analysis_h0.py``), plus the light curves
and lens-model products in ``data/``.
"""

from __future__ import annotations

import json
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import config as C
import delays as D
import lensmodel as LM
import lightcurves as L
import measure_delays as M
import plotstyle as P

P.use_style()
DL = C.RESULTS / "delays"


def _j(name):
    return json.loads((DL / f"{name}.json").read_text())


def final_model(key: str):
    """The joint-spline model at the BIC-best setting and the combined delays."""
    comb = _j(f"{key}_combined")
    ks = _j(f"{key}_knots")
    ml = (M.ML_SINGLE if key in M.SINGLE_SEASON else M.ML_MULTI)[0]
    run = _j(f"{key}_spline_{ks['best']}_{ml}")
    scale = {tuple(k.split("|")): v for k, v in run["err_scale"].items()}
    lc = L.read_lens(key)
    m = D.JointSpline(lc, comb["ref"], D.SplineSettings(knot_int=float(ks["best"]), knot_ml=float(ml)), err_scale=scale)
    fit = m.fit(comb["taus"])
    return lc, m, fit, comb


def _corrected(lc, m, fit, im):
    """Image ``im`` with its telescope offsets and microlensing removed."""
    return lc.mag[im] - m.offsets(fit.coef, im) - m.microlensing(fit.coef, im, lc.t)


def _tau(comb, im):
    return 0.0 if im == comb["ref"] else comb["taus"][im]


def _spread(ys: dict, min_sep: float) -> dict:
    """Nudge label positions apart so no two are closer than ``min_sep`` (keeps their order)."""
    keys = sorted(ys, key=ys.get)
    pos = [ys[k] for k in keys]
    for _ in range(50):
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < min_sep:
                mid = 0.5 * (pos[i] + pos[i - 1])
                pos[i - 1], pos[i] = mid - min_sep / 2, mid + min_sep / 2
                moved = True
        if not moved:
            break
    return dict(zip(keys, pos))


# ---------------------------------------------------------------------------
# 1-2. The flicker, and the flicker lined up
# ---------------------------------------------------------------------------

ZOOM = {"RXJ1131": (2006.0, 2008.6), "HE0435": (2010.6, 2014.6)}


def fig_flicker(key: str, number: int):
    lc, m, fit, comb = final_model(key)
    name = C.LENSES[key]["name"]
    yr = L.decimal_year(lc.t)
    fig = plt.figure(figsize=(11, 9.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1, 1], hspace=0.42)

    # (a) the raw curves, offset for clarity
    ax = fig.add_subplot(gs[0])
    base = np.median(lc.mag[lc.images[0]])
    offs = {}
    for i, im in enumerate(lc.images):
        offs[im] = base - np.median(lc.mag[im]) + 0.9 * i
        y = lc.mag[im] - m.offsets(fit.coef, im) + offs[im]
        ax.plot(yr, y, "o", ms=2.2, color=P.image_color(im), mec="none", alpha=0.8)
        ax.annotate(im, (yr.max() + 0.15, np.median(y[-30:])), color=P.INK, fontsize=9, fontweight="bold", va="center")
    ax.invert_yaxis()
    ax.set_ylabel("R magnitude (offset per image)")
    ax.set_title(f"{name}: {len(lc.images)} images of one quasar, {int(yr.min())}–{int(yr.max())}")
    ax.axvspan(*ZOOM[key], color=P.EXCLUDED, zorder=0)
    ax.text(ZOOM[key][0] + 0.05, ax.get_ylim()[1] + 0.02, "zoom below", fontsize=7.5, color=P.MUTED, va="bottom")
    ax.set_xlim(yr.min() - 0.3, yr.max() + 0.6)

    # (b) zoom, raw
    z0, z1 = ZOOM[key]
    sel = (yr >= z0) & (yr <= z1)
    ax = fig.add_subplot(gs[1])
    seg = np.concatenate([[0], np.cumsum(np.diff(yr[sel]) > 30 / 365.25)])  # break lines at season gaps
    ends = {}
    for im in lc.images:
        y = lc.mag[im] - m.offsets(fit.coef, im)
        y = (y - np.median(y[sel]))[sel]
        for sgi in np.unique(seg):
            k = seg == sgi
            ax.plot(yr[sel][k], y[k], "o-", ms=3, lw=0.6, color=P.image_color(im), mec="none")
        ends[im] = float(np.median(y[-12:]))
    lim = ax.get_ylim()
    for im, yy in _spread(ends, 0.075 * abs(lim[1] - lim[0])).items():
        ax.annotate(im, (z1 + 0.02, yy), color=P.INK, fontsize=9, fontweight="bold", va="center", annotation_clip=False)
    ax.invert_yaxis()
    ax.set_xlim(z0, z1)
    ax.set_ylabel("Δmag")
    ax.set_title("As observed: each image shows the same flicker, at a different time")

    # (c) shifted by the delays, microlensing removed
    ax = fig.add_subplot(gs[2])
    tt = np.linspace(lc.t.min() - 200, lc.t.max(), 6000)
    for im in lc.images:
        tau = _tau(comb, im)
        y = _corrected(lc, m, fit, im)
        ys = L.decimal_year(lc.t - tau)
        s2 = (ys >= z0) & (ys <= z1)
        ax.plot(ys[s2], y[s2], "o", ms=3.2, color=P.image_color(im), mec="none",
                label=f"{im} shifted by {-tau:+.1f} d" if im != comb["ref"] else f"{im} (reference)")
    s_curve = m.intrinsic(fit.coef, tt)
    yt = L.decimal_year(tt)
    # only where some shifted image constrains it (in gaps the spline is just the prior)
    shifted = np.sort(np.concatenate([lc.t - _tau(comb, im) for im in lc.images]))
    j = np.clip(np.searchsorted(shifted, tt), 1, len(shifted) - 1)
    near = np.minimum(np.abs(tt - shifted[j - 1]), np.abs(tt - shifted[j])) < 0.3 * m.s.knot_int
    s_curve = np.where(near, s_curve, np.nan)
    k = (yt >= z0) & (yt <= z1)
    ax.plot(yt[k], s_curve[k], color=P.INK, lw=1.0, label="fitted intrinsic quasar curve")
    ax.invert_yaxis()
    ax.set_xlim(z0, z1)
    ax.set_xlabel("year (arrival time of image " + comb["ref"] + ")")
    ax.set_ylabel("R magnitude")
    ax.set_title("Shifted by our delays, microlensing removed: one light curve")
    ax.legend(ncol=2, loc="lower right")
    P.caption(fig, f"Data: {lc.source}. Delays from our joint-spline fit (this project). "
                   "Microlensing = slow spline per image relative to image " + comb["ref"] + "; telescope zero points removed.")
    P.save(fig, C.FIG / f"{number:02d}_flicker_{key}")


# ---------------------------------------------------------------------------
# 3. Microlensing: stars in the lens galaxy
# ---------------------------------------------------------------------------


def fig_microlensing():
    fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=False)
    for ax, key in zip(axs, ["RXJ1131", "HE0435"]):
        lc, m, fit, comb = final_model(key)
        yr = L.decimal_year(lc.t)
        tt = np.linspace(lc.t.min(), lc.t.max(), 3000)
        for j, im in enumerate([i for i in lc.images if i != comb["ref"]]):
            tau = _tau(comb, im)
            # data minus the shifted intrinsic curve and offsets = this image's microlensing
            pts = lc.mag[im] - m.offsets(fit.coef, im) - m.intrinsic(fit.coef, lc.t - tau)
            mlc = m.microlensing(fit.coef, im, tt)
            zero = np.median(mlc) - 0.7 * j  # remove the constant flux ratio; stack the curves
            ax.plot(yr, pts - zero, "o", ms=1.6, color=P.image_color(im), alpha=0.35, mec="none")
            ax.plot(L.decimal_year(tt), mlc - zero, color=P.image_color(im), lw=2)
            ax.annotate(f"{im} vs {comb['ref']}", (L.decimal_year(tt[-1]) + 0.1, mlc[-1] - zero), color=P.INK,
                        fontsize=8.5, fontweight="bold", va="center")
        ax.invert_yaxis()
        ax.set_ylabel("brightness change (mag; curves offset)")
        ax.set_xlim(yr.min() - 0.2, yr.max() + 0.9)
        ax.set_title(f"{C.LENSES[key]['name']}: how each image brightened or faded relative to image {comb['ref']}")
    axs[1].set_xlabel("year")
    P.caption(fig, "Points: each image minus the shared intrinsic quasar curve (shifted by its delay), with the constant "
                   "brightness ratio between images (set by the lens's magnification) removed; curves offset by 0.7 mag. "
                   "Lines: the fitted microlensing splines. Changes of up to a magnitude over years come from stars in the lens galaxy "
                   "drifting across the line of sight; they differ between images, so they are not the quasar.")
    P.save(fig, C.FIG / "03_microlensing")


# ---------------------------------------------------------------------------
# 9. The mass-sheet degeneracy
# ---------------------------------------------------------------------------


def he0435_delay():
    """Our A->D delay of HE 0435 (days), falling back to the published one."""
    p = DL / "HE0435_combined.json"
    if p.exists():
        c = json.loads(p.read_text())
        return c["taus"]["D"], c["err"]["D"], "ours"
    return 13.8, 0.9, "Bonvin et al. 2017"


def slope_family(gammas=np.round(np.arange(1.70, 2.31, 0.05), 3)):
    """EPL + shear fits to HE 0435's image positions for fixed slopes; cached."""
    p = C.RESULTS / "he0435_slope_family.json"
    if p.exists():
        return json.loads(p.read_text())
    img = LM.he0435_positions()
    rows = []
    for g in gammas:
        mdl, chi2 = LM.fit_positions(img, float(g))
        d = LM.fermat_differences(mdl, img)
        rows.append(dict(gamma=float(g), chi2=chi2, b=mdl.b, q=mdl.q, phi_q=mdl.phi_q, g1=mdl.g1, g2=mdl.g2,
                         dphi={k: float(v) for k, v in d.items()}))
    C.RESULTS.mkdir(exist_ok=True)
    p.write_text(json.dumps(rows, indent=1))
    return rows


def fig_mass_sheet():
    import h0 as H

    img = LM.he0435_positions()
    rows = slope_family()
    iso = next(r for r in rows if abs(r["gamma"] - 2.0) < 1e-6)
    mdl = LM.EPLShear(iso["b"], iso["q"], iso["phi_q"], 2.0, 0.0, 0.0, iso["g1"], iso["g2"])
    bx, by = LM.source_of(mdl, img)
    zd, zs = C.LENSES["HE0435"]["z_d"], C.LENSES["HE0435"]["z_s"]
    _, ddt70 = H.distances(zd, zs, 0.3, 70.0)
    ddt70 = float(ddt70[0])
    day_per_arcsec2 = ddt70 * C.MPC_KM / C.C_KM_S / C.DAY_S * C.ARCSEC_RAD**2  # days per arcsec^2 at H0 = 70
    dt_AD, e_AD, src = he0435_delay()

    fig = plt.figure(figsize=(13.5, 8.6))
    gs = fig.add_gridspec(2, 3, wspace=0.32, hspace=0.38)

    # (a) arrival-time surface
    ax = fig.add_subplot(gs[0, 0])
    g = np.linspace(-2.0, 2.0, 401)
    X, Y = np.meshgrid(g, g)
    tau = mdl.fermat(X, Y, bx, by)
    tA = mdl.fermat(*img["A"], bx, by)
    T = (tau - tA) * day_per_arcsec2
    levels = np.arange(-4, 31, 2.0)
    cs = ax.contour(X, Y, T, levels=levels, colors=P.BLUE_300, linewidths=0.6)
    imdl = {k: (mdl.fermat(*xy, bx, by) - tA) * day_per_arcsec2 for k, xy in img.items()}
    for k, (x, y) in img.items():
        kind = "minimum" if k in "AC" else "saddle"
        ax.plot(x, y, "o", ms=8, color=P.image_color(k), mec=P.SURFACE, mew=1.5, zorder=5)
        ax.annotate(f"{k}  {imdl[k]:+.1f} d\n({kind})", (x, y), xytext=(6, 6), textcoords="offset points", fontsize=7.5, color=P.INK)
    ax.plot(0, 0, "+", color=P.INK, ms=9)
    ax.plot(bx, by, "*", color=P.INK_SECONDARY, ms=9)
    ax.annotate("source", (bx, by), xytext=(-30, -14), textcoords="offset points", fontsize=7.5, color=P.INK_SECONDARY)
    ax.set_aspect("equal")
    ax.set_xlabel("x (arcsec, west)")
    ax.set_ylabel("y (arcsec, north)")
    ax.set_title("(a) HE 0435: the arrival-time surface")
    ax.text(0.02, 0.02, "isothermal model fitted to the 4 HST positions;\ncontours every 2 days (H0 = 70)",
            transform=ax.transAxes, fontsize=7, color=P.MUTED)

    # cut through A and D (the first and last image)
    A, Dd = np.array(img["A"]), np.array(img["D"])
    u = np.linspace(-0.25, 1.25, 600)
    path = A[None, :] + u[:, None] * (Dd - A)[None, :]
    lams = [0.85, 1.0, 1.15]
    cols = [P.BLUE_300, P.INK, P.ORANGE]

    # (b) convergence profiles under the MST
    ax = fig.add_subplot(gs[0, 1])
    r = np.linspace(0.05, 3.0, 300)
    for lam, c in zip(lams, cols):
        mm = mdl.mst(lam)
        k_major = mm.convergence(r * np.cos(mdl.phi_q), r * np.sin(mdl.phi_q))
        ax.plot(r, k_major, color=c, lw=1.8, label=f"\u03bb = {lam:.2f}")
    for k in "ABCD":
        rk = np.hypot(*img[k])
        ax.axvline(rk, color=P.image_color(k), lw=0.8, alpha=0.6)
    ax.set_yscale("log")
    ax.set_xlabel("radius along the major axis (arcsec)")
    ax.set_ylabel("surface density \u03ba (critical units)")
    ax.set_title("(b) Three lenses: \u03ba \u2192 \u03bb\u03ba + (1 \u2212 \u03bb)")
    ax.legend(loc="upper right")
    ax.text(0.03, 0.04, "vertical lines: radii of images A\u2013D", transform=ax.transAxes, fontsize=7, color=P.MUTED)

    # (c) the same images, different delays
    ax = fig.add_subplot(gs[0, 2])
    for lam, c in zip(lams, cols):
        mm = mdl.mst(lam)
        bxl, byl = lam * bx, lam * by
        tcut = mm.fermat(path[:, 0], path[:, 1], bxl, byl)
        t0 = mm.fermat(*img["A"], bxl, byl)
        tD = mm.fermat(*img["D"], bxl, byl)
        ax.plot(u, (tcut - t0) * day_per_arcsec2, color=c, lw=1.8,
                label=f"\u03bb = {lam:.2f}: D arrives {(tD - t0) * day_per_arcsec2:.1f} d after A")
        # images stay put: check the lens equation at the observed positions
        sx, sy = mm.source(np.array([img[k][0] for k in "ABCD"]), np.array([img[k][1] for k in "ABCD"]))
        assert np.ptp(sx) < 0.05 and np.ptp(sy) < 0.05
    ax.axvline(0, color=P.image_color("A"), lw=0.8)
    ax.axvline(1, color=P.image_color("D"), lw=0.8)
    ax.text(0.0, 29.5, " A", color=P.INK, fontsize=8, va="top")
    ax.text(1.0, 29.5, " D", color=P.INK, fontsize=8, va="top")
    ax.set_xlabel("position along the line from A to D")
    ax.set_ylabel("arrival time relative to A (days, H0 = 70)")
    ax.set_title("(c) Same images, delays \u00d7 \u03bb")
    ax.legend(loc="lower right", fontsize=7)
    ax.set_ylim(-2, 30)

    # (d) H0 against the density slope (all fit the image positions)
    ax = fig.add_subplot(gs[1, 0])
    gam = np.array([r_["gamma"] for r_ in rows])
    ddt = np.array([LM.ddt_from(dt_AD, r_["dphi"]["D"]) for r_ in rows])
    h0v = 70.0 * ddt70 / ddt
    h0e = h0v * e_AD / dt_AD
    ax.fill_between(gam, h0v - h0e, h0v + h0e, color=P.BLUE_150, lw=0)
    ax.plot(gam, h0v, color=P.BLUE, lw=2, label=f"simple models, A\u2192D = {dt_AD:.1f} d ({src})")
    w = C.LENSES["HE0435"]["wong20_h0"]
    ax.errorbar([1.93], [w[0]], yerr=[[w[1]], [w[2]]], xerr=[[0.02], [0.02]], fmt="s", color=P.INK, ms=6, capsize=0,
                label="H0LiCOW model (\u03b3\u2032 from the arcs)")
    ax.set_xlabel("density slope \u03b3\u2032  (\u03c1 \u221d r^\u2212\u03b3\u2032; 2 = isothermal)")
    ax.set_ylabel("H0 from HE 0435 (km/s/Mpc)")
    ax.set_title("(d) Every slope fits the four images")
    ax.legend(loc="upper left", fontsize=7)
    chi = [r_["chi2"] for r_ in rows]
    ax.text(0.98, 0.04, f"position \u03c7\u00b2 = {min(chi):.1f}\u2013{max(chi):.1f} for 1 d.o.f.\n(\u03c3 = 3 mas per coordinate)",
            transform=ax.transAxes, ha="right", fontsize=7, color=P.MUTED)

    # (e) TDCOSMO 2025: H0 against lambda_int
    ax = fig.add_subplot(gs[1, 1])
    with h5py.File(C.DATA / "tdcosmo2025" / "chains_export" / "ULambdaCDM1.h5") as f:
        smp = f["samples"][:]
        names = [x.decode() for x in f["parameters"][:]]
    h = smp[:, names.index("h0")]
    lam = smp[:, names.index("lambda_mst")]
    Hh, xe, ye = np.histogram2d(lam, h, bins=60, range=[[0.8, 1.2], [55, 95]])
    ax.imshow(Hh.T, origin="lower", extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect="auto", cmap="Blues")
    med = np.median(h / lam)
    xx = np.linspace(0.8, 1.2, 10)
    ax.plot(xx, med * xx, color=P.ORANGE, lw=1.5, label=f"H0 = \u03bb \u00d7 {med:.1f}")
    ax.set_xlabel("\u03bb_int (internal mass sheet, population mean)")
    ax.set_ylabel("H0 (km/s/Mpc)")
    ax.set_title("(e) TDCOSMO 2025: H0 tracks \u03bb")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.text(0.98, 0.04, f"corr(H0, \u03bb) = {np.corrcoef(h, lam)[0, 1]:.2f}\n8 lenses, no external prior",
            transform=ax.transAxes, ha="right", fontsize=7, color=P.MUTED)

    # (f) what breaks it
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    ax.set_title("(f) What breaks the degeneracy")
    txt = ("Images, flux ratios and delay ratios are all\n"
           "unchanged by the mass-sheet transform; only the\n"
           "absolute delays (and hence H0) scale with \u03bb.\n\n"
           "\u2022 H0LiCOW assumed simple profiles (power law,\n"
           "   or stars + NFW halo), which fixes \u03bb = 1:\n"
           "   73.3 \u00b1 1.8 (2.4%).\n"
           "\u2022 TDCOSMO 2025 lets \u03bb float and measures it\n"
           "   with the stars' velocity dispersion (JWST, Keck,\n"
           "   VLT): \u03bb ~ 0.97 \u00b1 0.05, H0 = 71.6 (+3.9/\u22123.3).\n"
           "\u2022 Mass along the line of sight acts as an\n"
           "   external sheet, \u03bb_ext = 1 \u2212 \u03ba_ext, estimated\n"
           "   from galaxy counts around each lens.")
    ax.text(0.0, 0.98, txt, va="top", fontsize=8.2, color=P.INK, transform=ax.transAxes, linespacing=1.35)
    P.caption(fig, "HE 0435\u22121223 image positions: Courbin et al. 2011 (HST/NICMOS). Models: our elliptical power law + "
                   "external shear (Tessore & Metcalf 2015 deflections). TDCOSMO 2025 chain: ULambdaCDM1 (github.com/TDCOSMO/TDCOSMO2025_public).")
    P.save(fig, C.FIG / "09_mass_sheet_degeneracy")


# ---------------------------------------------------------------------------
# 10. The gallery: all 23 COSMOGRAIL XIX lenses
# ---------------------------------------------------------------------------

#: Millon et al. 2020a, Table 3 (COSMOGRAIL convention: Delta t_XY = t_X - t_Y,
#: negative when X leads): the best-measured delay of each system.
MILLON20_TABLE3 = {
    "HE0047_Euler": ("A", "B", -10.4, 3.5, 3.5), "UM673_Euler": ("A", "B", -97.7, 15.5, 16.1),
    "J0158_Euler": ("A", "B", -22.7, 3.6, 3.6), "HE0230_Euler": ("A+B", "C", 15.7, 3.6, 4.2),
    "J0246_Euler": ("A", "B", 0.8, 5.2, 5.0), "HE0435_all": ("A", "D", -13.8, 0.8, 0.8),
    "HS0818_Euler": ("A", "B", -153.8, 14.6, 13.2), "J0832_Euler": ("A", "B", -125.3, 23.4, 12.8),
    "J0924_all": ("A+D", "B", 2.4, 3.8, 3.8), "RXJ1131_all": ("A", "D", -92.5, 1.8, 1.9),
    "J1226_Euler": ("A", "B", 33.7, 2.7, 2.7), "J1335_Euler": ("A", "B", -56.0, 6.1, 5.7),
    "Q1355_Euler": ("A", "B", -81.5, 12.0, 10.8), "J1455_Euler": ("A", "B", -47.2, 7.8, 7.5),
    "J1515_Euler": ("A", "B", -210.2, 5.7, 5.5), "J1620_Euler": ("A", "B", -171.5, 8.7, 8.7),
    "WFI2026_Euler": ("A", "B", 18.7, 4.9, 4.7), "HE2149_Euler": ("A", "B", -39.0, 16.7, 14.9),
}


def fig_gallery():
    lcs = L.read_gallery()
    ncol = 4
    nrow = int(np.ceil(len(lcs) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(13.5, 2.05 * nrow))
    for ax, lc in zip(axs.ravel(), lcs):
        yr = L.decimal_year(lc.t)
        for i, im in enumerate(lc.images):
            y = lc.mag[im] - np.median(lc.mag[im]) + 0.45 * i
            ax.plot(yr, y, "o", ms=1.1, color=P.image_color(im.replace("+", "")[0]), mec="none")
            ax.text(2019.3, np.median(y[-25:]), im, fontsize=6.5, color=P.INK, va="center", clip_on=False)
        ax.invert_yaxis()
        ax.set_xlim(2004, 2019)
        ax.tick_params(labelsize=6.5)
        ax.grid(False)
        ax.set_title(lc.name.replace("2017WFI", "WFI").replace("-", "\u2212"), fontsize=8.5, loc="left", pad=3)
        t3 = MILLON20_TABLE3.get(lc.key)
        if t3:
            a, b, v, em, ep = t3
            first, last = (a, b) if v < 0 else (b, a)
            txt = f"{last} trails {first} by {abs(v):.1f} d"
        else:
            txt = "no delay measured"
        ax.text(0.01, 0.03, txt, transform=ax.transAxes, fontsize=6.8, color=P.INK_SECONDARY)
    for ax in axs.ravel()[len(lcs):]:
        ax.axis("off")
    fig.suptitle("Fifteen years of flickering: the 23 lensed quasars of COSMOGRAIL XIX", x=0.01, ha="left", y=1.0)
    fig.tight_layout()
    P.caption(fig, "R-band light curves from the 1.2 m Euler telescope (plus SMARTS, Mercator, Maidanak and Liverpool data "
                   "for some systems), Millon et al. 2020a, VizieR J/A+A/640/A105. Images offset by 0.45 mag; delays from "
                   "the paper's Table 3. Q 2237+0305 (the Einstein Cross) has delays under a day.")
    P.save(fig, C.FIG / "10_gallery")


# ---------------------------------------------------------------------------
# Delay bookkeeping shared by figures 4-6
# ---------------------------------------------------------------------------

#: Rows of the delay comparison: (lens, first image, second image).
DELAY_ROWS = [("RXJ1131", "A", "B"), ("RXJ1131", "A", "C"), ("RXJ1131", "A", "D"),
              ("HE0435", "A", "B"), ("HE0435", "A", "C"), ("HE0435", "A", "D"),
              ("PG1115", "A", "B"), ("PG1115", "A", "C"),
              ("WFI2033", "B", "A"), ("WFI2033", "B", "C"),
              ("J1206", "A", "B"),
              ("DES0408", "A", "B"), ("DES0408", "A", "D")]


def published_delay(key, a, b):
    """(value, err_minus, err_plus, source) of t_b - t_a from the literature."""
    Lk = C.LENSES[key]
    for table, src in ((Lk.get("cosmograil19"), "Millon et al. 2020a"), (Lk.get("model_delays"), Lk["model_delay_ref"])):
        if not table:
            continue
        if (a, b) in table:
            v, em, ep = table[(a, b)]
            return v, em, ep, src
        if (b, a) in table:
            v, em, ep = table[(b, a)]
            return -v, ep, em, src
    if key == "WFI2033" and (a, b) == ("B", "A"):  # A1 and A2 averaged (merged in most data)
        return 36.75, 2.1, 2.65, "Bonvin et al. 2019 (A1, A2 averaged)"
    if key == "PG1115" and (a, b) == ("A", "B"):
        v, em, ep = Lk["model_delays"][("A", "B")]
        return v, em, ep, Lk["model_delay_ref"]
    raise KeyError((key, a, b))


def pair_from(res_taus, ref, a, b):
    ta = 0.0 if a == ref else res_taus[a]
    tb = 0.0 if b == ref else res_taus[b]
    return tb - ta


def pair_err(comb, a, b):
    ref, others = comb["ref"], comb["others"]
    cov = np.array(comb["cov"])
    v = np.zeros(len(others))
    if b != ref:
        v[others.index(b)] += 1
    if a != ref:
        v[others.index(a)] -= 1
    return float(np.sqrt(v @ cov @ v))


def method_values(key, a, b):
    """Our delay for (a, b) from each method: dict name -> (value, err_minus, err_plus) or value."""
    comb = _j(f"{key}_combined")
    ref = comb["ref"]
    out = {"joint spline (ours)": (pair_from(comb["taus"], ref, a, b),) + (pair_err(comb, a, b),) * 2}
    cc = _j(f"{key}_ccf")["delays"]
    if a == ref or b == ref:
        im, sgn = (b, 1) if a == ref else (a, -1)
        d = cc[im]
        out["ICCF (ours)"] = (sgn * d["med"], d["minus"] if sgn > 0 else d["plus"], d["plus"] if sgn > 0 else d["minus"])
    for f in sorted(DL.glob(f"{key}_pycs_*.json")):
        r = json.loads(f.read_text())
        lab = "PyCS3 spline" if "spl" in r["method"] else "PyCS3 regdiff"
        out.setdefault(lab, [])
        out[lab].append(pair_from(r["taus"], ref, a, b))
    return out


# ---------------------------------------------------------------------------
# 4. Our delays against the published ones
# ---------------------------------------------------------------------------


XLIM = 16.0


def _offscale(ax, dx, y, color):
    """An arrow at the edge of the axis for a point beyond it, labelled with its value."""
    x = np.sign(dx) * (XLIM - 0.6)
    ax.annotate("", (x + np.sign(dx) * 0.5, y), (x - np.sign(dx) * 0.4, y),
                arrowprops=dict(arrowstyle="->", color=color, lw=1.2))
    ax.text(x - np.sign(dx) * 0.6, y, f"{dx:+.0f}", color=P.INK_SECONDARY, fontsize=7, va="center",
            ha="right" if dx > 0 else "left")


def fig_delays():
    fig, ax = plt.subplots(figsize=(10.5, 8.2))
    rows = DELAY_ROWS[::-1]
    ylabels = []
    for i, (key, a, b) in enumerate(rows):
        v, em, ep, src = published_delay(key, a, b)
        ax.fill_betweenx([i - 0.38, i + 0.38], -em, ep, color=P.EXCLUDED, lw=0)
        mv = method_values(key, a, b)
        y0 = mv["joint spline (ours)"]
        ax.errorbar(y0[0] - v, i + 0.18, xerr=[[y0[1]], [y0[2]]], fmt="o", color=P.C_OURS, ms=6, capsize=0, lw=1.6, zorder=4)
        if "ICCF (ours)" in mv:
            c_ = mv["ICCF (ours)"]
            if abs(c_[0] - v) < XLIM:
                ax.errorbar(c_[0] - v, i - 0.02, xerr=[[c_[1]], [c_[2]]], fmt="D", color=P.C_CCF, ms=4.5, capsize=0, lw=1.0, zorder=3)
            else:
                _offscale(ax, c_[0] - v, i - 0.02, P.C_CCF)
        for lab, mk in (("PyCS3 spline", "^"), ("PyCS3 regdiff", "v")):
            for val in mv.get(lab, []):
                if abs(val - v) < XLIM:
                    ax.plot(val - v, i - 0.2, mk, color=P.C_PYCS, ms=5, mec=P.SURFACE, mew=0.6, zorder=3)
                else:
                    _offscale(ax, val - v, i - 0.2, P.C_PYCS)
        name = C.LENSES[key]["name"]
        ylabels.append(f"{name}  {a}\u2192{b}:  {v:+.1f} d")
    ax.axvline(0, color=P.MUTED, lw=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("our delay \u2212 published delay (days)")
    ax.set_xlim(-XLIM, XLIM)
    ax.grid(axis="y", visible=False)
    ax.set_title("Thirteen delays, three methods, against the published values")
    handles = [Line2D([], [], marker="o", color=P.C_OURS, lw=1.6, label="joint spline + microlensing (ours, total error)"),
               Line2D([], [], marker="D", color=P.C_CCF, lw=1.0, ms=4.5, label="high-passed ICCF (ours, FR/RSS error)"),
               Line2D([], [], marker="^", color=P.C_PYCS, lw=0, label="PyCS3 free-knot spline (two settings)"),
               Line2D([], [], marker="v", color=P.C_PYCS, lw=0, label="PyCS3 regression difference"),
               plt.Rectangle((0, 0), 1, 1, color=P.EXCLUDED, label="published 1\u03c3 range")]
    ax.legend(handles=handles, loc="upper left", fontsize=7.5)
    P.caption(fig, "Published: Millon et al. 2020a (RX J1131, HE 0435, same data), Bonvin et al. 2018 (PG 1115), 2019 "
                   "(WFI 2033), Birrer et al. 2019 (SDSS J1206), Courbin et al. 2018 (DES J0408). Delay t(second) \u2212 t(first): "
                   "positive = the second image arrives later.")
    P.save(fig, C.FIG / "04_delays_compared")


# ---------------------------------------------------------------------------
# 5. Finding the delay: wide cross-correlation, then the evidence profile
# ---------------------------------------------------------------------------

FIND = [("RXJ1131", "A", "D"), ("HE0435", "A", "D"), ("WFI2033", "B", "C"), ("J1206", "A", "B")]


def evidence_profile(key, im, half=12.0, step=0.5):
    """-2 ln evidence against one delay (others re-optimized), at the BIC-best setting; cached."""
    p = DL / f"{key}_profile_{im}.json"
    if p.exists():
        return json.loads(p.read_text())
    lc, m, fit, comb = final_model(key)
    centre = comb["taus"][im]
    grid = np.arange(centre - half, centre + half + step / 2, step)
    vals, cur = [], dict(comb["taus"])
    rest = [o for o in m.others if o != im]
    for g in grid:
        if rest:
            sub = D.JointSplineProfile(m, im, float(g))
            cur.update(sub.optimize({o: cur[o] for o in rest}, width=3.0))
        vals.append(m.chi2({**cur, im: float(g)}))
    out = dict(grid=grid.tolist(), obj=(np.array(vals) - min(vals)).tolist())
    p.write_text(json.dumps(out))
    return out


def fig_find():
    fig, axs = plt.subplots(2, len(FIND), figsize=(13.5, 6.4))
    for j, (key, a, b) in enumerate(FIND):
        st = _j(f"{key}_start")
        k = f"{min(a, b)}{max(a, b)}"
        pr = st["pairs"][k]
        sgn = 1 if a < b else -1
        lags = sgn * np.array(pr["lags"])
        ax = axs[0, j]
        ax.plot(lags, pr["r"], color=P.C_CCF, lw=1.4)
        v = published_delay(key, a, b)[0]
        ax.axvline(v, color=P.MUTED, lw=1, ls="-")
        ax.axvline(0, color=P.GRID, lw=0.8)
        ax.set_title(f"{C.LENSES[key]['name']}  {a}\u2192{b}", fontsize=9.5)
        ax.set_xlabel("trial delay (days)")
        if j == 0:
            ax.set_ylabel("cross-correlation r")
        comb = _j(f"{key}_combined")
        im = b if a == comb["ref"] else a
        prof = evidence_profile(key, im)
        g = np.array(prof["grid"])
        if a != comb["ref"]:
            g = -g
        ax = axs[1, j]
        ax.plot(g, prof["obj"], color=P.C_OURS, lw=1.8)
        best = pair_from(comb["taus"], comb["ref"], a, b)
        e = pair_err(comb, a, b)
        ax.axvspan(best - e, best + e, color=P.BLUE_150, lw=0, zorder=0)
        ax.axvline(v, color=P.MUTED, lw=1)
        ax.set_ylim(0, max(60, 1.1 * min(max(prof["obj"]), 400)))
        ax.set_xlabel("delay (days)")
        if j == 0:
            ax.set_ylabel("\u0394(\u22122 ln evidence)")
        ax.text(0.5, 0.97, f"ours {best:.1f} \u00b1 {e:.1f} d\npublished {v:.1f} d", transform=ax.transAxes, fontsize=7.5,
                va="top", ha="center", color=P.INK)
    axs[0, 0].text(0.02, 0.97, "blind search,\n\u00b1200 days", transform=axs[0, 0].transAxes, fontsize=7.5, va="top", color=P.MUTED)
    fig.suptitle("Finding a delay: a blind cross-correlation (top), then the joint spline fit (bottom)", x=0.01, ha="left")
    fig.tight_layout()
    P.caption(fig, "Top: correlation of the high-passed light curves at every trial delay (\u00b1200 d, no prior knowledge); "
                   "the peaks are broad because the filter that removes microlensing keeps only month-scale wiggles. Bottom: "
                   "the joint spline's \u22122 ln evidence at its preferred setting (other delays re-fitted). Its formal width "
                   "(\u0394 = 1 within 0.4\u20131 d) is too optimistic because the residuals are correlated (fast microlensing); "
                   "the errors therefore come from simulated light curves. Blue band = our final \u00b11\u03c3 (mean of six "
                   "settings, with their spread); grey line = published.")
    P.save(fig, C.FIG / "05_finding_the_delay")


# ---------------------------------------------------------------------------
# 6. The error budget
# ---------------------------------------------------------------------------


def fig_errors():
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw=dict(width_ratios=[1, 1.5]))
    ax = axs[0]
    ks = _j("RXJ1131_knots")
    run = _j(f"RXJ1131_spline_{ks['best']}_{M.ML_MULTI[0]}")
    dev = np.array(run["mock_dev"]["D"])
    ax.hist(dev, bins=25, color=P.BLUE_300, edgecolor=P.SURFACE)
    ax.axvline(0, color=P.INK, lw=1)
    ax.axvline(dev.mean(), color=P.ORANGE, lw=1.5)
    ax.set_xlabel("measured \u2212 true delay on mock curves (days)")
    ax.set_ylabel("mock light curves")
    ax.set_title("RX J1131 A\u2192D: 100 simulated data sets")
    ax.text(0.97, 0.95, f"bias {dev.mean():+.2f} d\nscatter {dev.std(ddof=1):.2f} d", transform=ax.transAxes, ha="right",
            va="top", fontsize=8.5)
    ax = axs[1]
    rows = DELAY_ROWS[::-1]
    for i, (key, a, b) in enumerate(rows):
        comb = _j(f"{key}_combined")
        mean = pair_from(comb["taus"], comb["ref"], a, b)
        for sset in comb["settings"]:
            val = pair_from(sset["taus"], comb["ref"], a, b)
            ml_long = sset["knot_ml"] == max(st["knot_ml"] for st in comb["settings"])
            ax.plot(val - mean, i, "o" if not ml_long else "s", ms=4.5, color=P.BLUE_300 if not ml_long else P.BLUE_700,
                    mec="none", alpha=0.9)
        e_stat = pair_err(dict(comb, cov=(np.array(comb["cov"]) - np.diag([comb["err_sys"][o] ** 2 for o in comb["others"]])).tolist()), a, b)
        e_tot = pair_err(comb, a, b)
        ax.plot([-e_tot, e_tot], [i, i], color=P.INK, lw=0.8, alpha=0.5)
        ax.plot([-e_stat, e_stat], [i, i], color=P.INK, lw=2.4, alpha=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{C.LENSES[k]['name']} {a}\u2192{b}" for k, a, b in rows], fontsize=8)
    ax.axvline(0, color=P.MUTED, lw=1)
    ax.set_xlabel("estimate at each setting \u2212 our final delay (days)")
    ax.set_title("Six estimator settings per lens")
    ax.grid(axis="y", visible=False)
    handles = [Line2D([], [], marker="o", lw=0, color=P.BLUE_300, label="short microlensing knots"),
               Line2D([], [], marker="s", lw=0, color=P.BLUE_700, label="long microlensing knots"),
               Line2D([], [], color=P.INK, lw=2.4, alpha=0.6, label="\u00b1 random error (mocks)"),
               Line2D([], [], color=P.INK, lw=0.8, alpha=0.5, label="\u00b1 total (with setting scatter)")]
    ax.legend(handles=handles, fontsize=7.5, loc="upper right")
    fig.tight_layout()
    P.caption(fig, "Left: the estimator run blind on mock curves (our fitted model resampled at the real dates, with the real "
                   "residuals rolled in time) whose true delays are known. Right: the three intrinsic knot steps around the "
                   "data's preferred one (lowest BIC), each with two microlensing knot steps.")
    P.save(fig, C.FIG / "06_error_budget")


# ---------------------------------------------------------------------------
# 7-8. H0
# ---------------------------------------------------------------------------


def fig_h0_lenses():
    post = np.load(C.RESULTS / "h0_posteriors.npz")
    res = json.loads((C.RESULTS / "h0.json").read_text())
    h = post["h0"]
    keys = C.H0LICOW_SIX
    fig, axs = plt.subplots(2, 4, figsize=(13.5, 6.0), sharex=True)
    for ax, key in zip(axs.ravel(), keys + ["ALL"]):
        pp = post[f"pub_{key}"] if key != "ALL" else post["pub_ALL"]
        po = post[f"ours_{key}"] if key != "ALL" else post["ours_ALL"]
        ax.fill_between(h, pp, color=P.C_REST, lw=0, label="published delays")
        if key != "B1608":
            ax.plot(h, po, color=P.C_OURS, lw=1.8, label="our delays")
        a = res["published"][key]
        o = res["ours"][key]
        title = "All six combined" if key == "ALL" else C.LENSES[key]["name"]
        ax.set_title(title, fontsize=9.5)
        txt = f"published {a[0]:.1f} \u2212{a[1]:.1f}/+{a[2]:.1f}"
        if key != "B1608":
            txt += f"\nours {o[0]:.1f} \u2212{o[1]:.1f}/+{o[2]:.1f}"
        else:
            txt += "\n(radio delays, not re-measured)"
        ax.text(0.97, 0.95, txt, transform=ax.transAxes, fontsize=7.4, va="top", ha="right", color=P.INK)
        ax.set_yticks([])
        ax.set_xlim(50, 115)
    ax = axs.ravel()[-1]
    ax.axis("off")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=P.C_REST, label="published delays"),
                       Line2D([], [], color=P.C_OURS, lw=1.8, label="our delays")], loc="lower left", fontsize=8)
    rows = [f"{C.LENSES[k]['name']}: delay scale r = {res['swap'][k]['r']:.3f}" for k in C.MEASURED if k in res["swap"]]
    ax.text(0.0, 0.95, "Our delays / the model's delays\n(D_dt and H0 scale as r and 1/r):\n\n" + "\n".join(rows),
            va="top", fontsize=8, color=P.INK, transform=ax.transAxes, linespacing=1.4)
    for ax in axs[1]:
        ax.set_xlabel("H0 (km/s/Mpc)")
    fig.suptitle("H0 from each lens: the published lens models with the published delays, and with ours", x=0.01, ha="left")
    fig.tight_layout()
    P.caption(fig, "Flat \u039bCDM, uniform priors H0 \u2208 [0, 150], \u03a9m \u2208 [0.05, 0.5] (as Wong et al. 2020). Distance "
                   "posteriors: H0LiCOW public release (Suyu+2010/2014, Jee+2019, Wong+2017, Chen+2019, Birrer+2019, Rusu+2020).")
    P.save(fig, C.FIG / "07_h0_per_lens")


def fig_h0_routes():
    res = json.loads((C.RESULTS / "h0.json").read_text())
    td = res["tdcosmo2025"]
    rows = [
        ("Lensed quasars: H0LiCOW models + our delays", res["ours"]["ALL"], P.C_LENS, "o", "this project"),
        ("H0LiCOW 2020, reproduced (power law / composite)", res["published"]["ALL"], P.C_LENS, "s", "this project"),
        ("TDCOSMO 2025 + our delays (estimate)", td["LambdaCDM1a_with_our_delays_estimate"], P.C_LENS, "o", "estimate"),
        ("TDCOSMO 2025 (\u03bb free, kinematics) + Pantheon+ \u03a9m", td["LambdaCDM1a"]["h0"], P.C_LENS, "D", "published chain"),
        ("TDCOSMO 2025 + SLACS + SL2S + Pantheon+", td["LambdaCDM1d"]["h0"], P.C_LENS, "D", "published chain"),
        ("D1: Cepheids + Type Ia supernovae", C.H0_OTHER["D1 Cepheids + SNe Ia (this repo)"], P.C_LADDER, "o", "this repo"),
        ("D3: Planck CMB (early universe)", C.H0_OTHER["D3 Planck CMB (this repo)"], P.C_CMB, "o", "this repo"),
        ("C4: GW170817 standard siren", C.H0_OTHER["C4 GW170817 standard siren (Abbott+2017)"], P.C_SIREN, "^", "Abbott+2017"),
        ("I9: dark + bright sirens (GWTC-5.0)", C.H0_OTHER["I9 GWTC-5.0 dark + bright sirens (LVK 2026)"], P.C_SIREN, "^", "LVK 2026"),
    ]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    d1 = C.H0_OTHER["D1 Cepheids + SNe Ia (this repo)"]
    d3 = C.H0_OTHER["D3 Planck CMB (this repo)"]
    ax.axvspan(d1[0] - d1[1], d1[0] + d1[2], color=P.ORANGE, alpha=0.12, lw=0)
    ax.axvspan(d3[0] - d3[1], d3[0] + d3[2], color=P.VIOLET, alpha=0.12, lw=0)
    n = len(rows)
    for i, (lab, (v, em, ep), col, mk, src) in enumerate(rows):
        y = n - 1 - i
        ax.errorbar(v, y, xerr=[[em], [ep]], fmt=mk, color=col, ms=7 if i == 0 else 5.5, lw=2 if i == 0 else 1.3, capsize=0)
        nd = 1 if all(abs(round(x, 1) - x) < 1e-6 for x in (v, em, ep)) or i < 5 else 2
        ax.text(101.5, y, f"{v:.{nd}f} \u2212{em:.{nd}f}/+{ep:.{nd}f}", va="center", fontsize=8, color=P.INK)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=8.5)
    ax.set_xlim(55, 100)
    ax.set_xlabel("H0 (km/s/Mpc)")
    ax.grid(axis="y", visible=False)
    ax.set_title("Four ways to measure the expansion rate today")
    P.caption(fig, "Blue: time-delay lenses (no distance ladder). Orange band: D1's Cepheid-calibrated supernovae (this "
                   "repository). Violet band: D3's Planck fit in \u039bCDM (this repository). Sirens: C4 and I9 are not built "
                   "yet, so their published values are shown. The TDCOSMO 2025 + our delays row rescales its H0 by our "
                   "delay ratios (an estimate, not a re-run of their hierarchical model).")
    P.save(fig, C.FIG / "08_h0_routes")


FIGURES = {
    "01_flicker_RXJ1131": lambda: fig_flicker("RXJ1131", 1),
    "02_flicker_HE0435": lambda: fig_flicker("HE0435", 2),
    "03_microlensing": fig_microlensing,
    "04_delays_compared": fig_delays,
    "05_finding_the_delay": fig_find,
    "06_error_budget": fig_errors,
    "07_h0_per_lens": fig_h0_lenses,
    "08_h0_routes": fig_h0_routes,
    "09_mass_sheet": fig_mass_sheet,
    "10_gallery": fig_gallery,
}


def main(argv):
    C.FIG.mkdir(exist_ok=True)
    pick = argv[1:] or [""]
    for name, fn in FIGURES.items():
        if any(p in name for p in pick):
            print("figure", name, flush=True)
            fn()


if __name__ == "__main__":
    main(sys.argv)
