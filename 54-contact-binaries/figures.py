"""Every figure of the project, drawn from data/results.json and data/derived/.

    python figures.py            # all figures
    python figures.py 07 12      # only those numbers
"""

from __future__ import annotations

import json
import os
import sys
import warnings

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
warnings.filterwarnings("ignore")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm

import config as C
import literature as L
import plotstyle as ps

ps.use_style()

#: Sequential ramp (one hue, light to dark) for densities.
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#eef4fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
#: Grey-to-warm-white ramp for rendered stellar surfaces on a dark panel.
STAR = LinearSegmentedColormap.from_list("star", ["#0b0b0b", "#3d2a1c", "#a8754a", "#f2c99a", "#fff6ea"])


def _out(stem: str) -> str:
    return str(C.FIGURE_DIR / stem)


def _res() -> dict:
    return json.loads(C.RESULTS_PATH.read_text())


def _panel(ax, letter, title=None):
    ax.set_title(f"{letter}  {title}" if title else letter, loc="left")


# ---------------------------------------------------------------------------
# 01  Roche lobes
# ---------------------------------------------------------------------------


def _stars_only(X, Y, F, level):
    """F with everything outside the stars' own connected region (Omega >= level)
    pushed below ``level``, so a contour at ``level`` draws only the stellar
    surface and not its twin far out in the orbital plane, where the
    centrifugal term makes Omega large again."""
    from scipy import ndimage
    lab, _ = ndimage.label(F >= level)
    keep = set()
    for cx in (0.0, 1.0):
        j = np.abs(Y[:, 0]).argmin()
        i = np.abs(X[0, :] - cx).argmin()
        if lab[j, i] > 0:
            keep.add(lab[j, i])
    inside = np.isin(lab, list(keep))
    G = F.copy()
    G[(F >= level) & ~inside] = level - 1.0
    return G, inside


def fig01_roche_lobes():
    import roche as R
    res = _res()["roche"]
    fig = plt.figure(figsize=(11.5, 7.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1.25, 1.0], height_ratios=[1, 1], wspace=0.28, hspace=0.32)
    for row, q in enumerate((1.0, 0.1)):
        ax = fig.add_subplot(gs[row, 0:2])
        lp = R.lagrange_points(q)
        X, Y, F = R.orbital_plane_contours(q, None, extent=(-1.55, 2.45, -1.25, 1.25), n=900)
        om1, om2, om3 = lp["omega_L1"], lp["omega_L2"], lp["omega_L3"]
        # background equipotentials (recessive): a few inside each star, a few outside
        for lev in np.concatenate([np.linspace(om1, om1 + 2.5, 5)[1:], [om3 - 0.15, om3 - 0.4]]):
            G, _ = _stars_only(X, Y, F, lev)
            ax.contour(X, Y, G, levels=[lev], colors=ps.AXIS, linewidths=0.6)
        # a contact envelope at fill-out 0.5, shaded
        omf = R.fillout_potential(q, 0.5)
        _, inside = _stars_only(X, Y, F, omf)
        ax.contourf(X, Y, inside.astype(float), levels=[0.5, 1.5], colors=["#cde2fb"])
        G1, _ = _stars_only(X, Y, F, om1)
        ax.contour(X, Y, G1, levels=[om1], colors=[ps.BLUE], linewidths=1.6)
        # a hair inside the outer critical surface: for q = 1 the L2 and L3 saddles
        # have the same potential and the region would leak out through L3
        G2, _ = _stars_only(X, Y, F, om2 + 2e-3)
        ax.contour(X, Y, G2, levels=[om2 + 2e-3], colors=[ps.ORANGE], linewidths=1.4)
        # Lagrange points
        for name in ("L1", "L2", "L3"):
            ax.plot(lp[name], 0, "o", ms=4.5, color=ps.INK, zorder=5)
            ax.annotate(name, (lp[name], 0), xytext=(4, 6), textcoords="offset points", fontsize=8, color=ps.INK)
        for sgn, name in ((1, "L4"), (-1, "L5")):
            ax.plot(0.5, sgn * np.sqrt(3) / 2, "o", ms=4.5, color=ps.INK_SECONDARY, zorder=5)
            ax.annotate(name, (0.5, sgn * np.sqrt(3) / 2), xytext=(5, -3), textcoords="offset points", fontsize=8,
                        color=ps.INK_SECONDARY)
        ax.plot([0, 1], [0, 0], "+", ms=9, color=ps.INK, mew=1.0)
        cm = q / (1 + q)
        ax.plot(cm, 0, "x", ms=6, color=ps.INK_SECONDARY, mew=1.0)
        ax.set_aspect("equal")
        ax.set_xlim(-1.55, 2.45)
        ax.set_ylim(-1.15, 1.15)
        ax.set_xlabel("x / separation")
        ax.set_ylabel("y / separation")
        ax.grid(False)
        _panel(ax, "ab"[row], f"mass ratio q = {q:g}  (orbital plane)")
        if row == 0:
            ps.annotate(ax, "Roche lobes: touch at L1", (0.15, 0.36), (-1.45, 0.92), color=ps.BLUE)
            ps.annotate(ax, "outer critical surface (through L2)", (1.62, 0.2), (1.25, 0.98), color=ps.ORANGE)
            ps.annotate(ax, "a contact binary's common envelope (fill-out 0.5)", (0.75, -0.33), (-1.45, -1.0),
                        color=ps.INK_SECONDARY)
        else:
            ps.annotate(ax, "gas lost here leaves the binary", (lp["L2"], 0.0), (1.05, -0.95), color=ps.INK_SECONDARY)
            ax.text(-1.45, 0.98, "+ star centres   x centre of mass", fontsize=7.5, color=ps.MUTED)
    # (c) lobe radius vs mass ratio: Eggleton's formula and our integration
    ax = fig.add_subplot(gs[:, 2])
    qq = np.logspace(-2, 2, 400)
    egg = np.array([R.eggleton_radius(x) for x in qq])
    ax.plot(qq, egg, color=ps.MUTED, lw=1.6, label="Eggleton (1983) formula")
    lobes = res["lobes"]
    q_star = [d["q"] for d in lobes] + [1 / d["q"] for d in lobes]
    r_num = [d["r2"] for d in lobes] + [d["r1"] for d in lobes]
    ax.plot(q_star, r_num, "o", ms=6, color=ps.BLUE, mec=ps.SURFACE, mew=1.2, label="our numerical lobes", zorder=4)
    ax.set_xscale("log")
    ax.set_xlabel("q = M(star) / M(companion)")
    ax.set_ylabel("Roche-lobe radius / separation  (equal volume)")
    _panel(ax, "c", "how big can a star get?")
    worst = max(abs(d["r1"] / d["r1_eggleton"] - 1) for d in lobes) * 100
    worst = max(worst, max(abs(d["r2"] / d["r2_eggleton"] - 1) for d in lobes) * 100)
    ax.text(0.04, 0.95, f"numerical volumes within {worst:.1f}%\nof the 1983 fitting formula", transform=ax.transAxes,
            va="top", fontsize=8, color=ps.INK_SECONDARY)
    ax.legend(loc="lower right")
    ps.caption(fig, "Kopal's Roche potential for a circular, synchronised orbit (separation 1, heavier star at the "
                    "origin). A star that outgrows its lobe spills gas through L1; stars that both fill theirs share "
                    "an envelope. Fill-out f = (Ω_L1 − Ω)/(Ω_L1 − Ω_L2).", y=0.01)
    return ps.save(fig, _out("01_roche_lobes"))


# ---------------------------------------------------------------------------
# 02  What a contact binary looks like, and its light curve
# ---------------------------------------------------------------------------


def fig02_contact_light_curve():
    import roche as R
    d = np.load(C.DERIVED_DIR / "roche_kic_model.npz")
    m = L.MOLNAR2017
    cb = R.ContactBinary(m["q"], m["fillout"], t1=m["t1"], t2=m["t2"], lam_nm=640.0, ld=(0.40, 0.26))
    phases = [0.0, 0.125, 0.25, 0.375, 0.5]
    fig = plt.figure(figsize=(11.5, 6.4))
    gs = fig.add_gridspec(2, 5, height_ratios=[1.0, 1.35], hspace=0.35, wspace=0.05)
    vmax = None
    for k, ph in enumerate(phases):
        ax = fig.add_subplot(gs[0, k])
        img, hit, _ = cb.render(m["incl"], ph, npix=260)
        if vmax is None:
            vmax = img.max()
        ax.imshow(np.where(hit, img, np.nan), origin="lower", cmap=STAR, vmin=0, vmax=vmax,
                  extent=(-cb.rb, cb.rb, -cb.rb, cb.rb))
        ax.set_facecolor("#0b0b0b")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlabel(f"phase {ph:g}", fontsize=9, color=ps.INK_SECONDARY, labelpad=3)
        if k == 0:
            ax.set_title("a  KIC 9832227 as our ray tracer sees it, from phase 0 to 0.5", loc="left", fontsize=11)
    ax = fig.add_subplot(gs[1, :])
    ph = d["phase"]
    for shift in (0.0, 1.0):
        ax.plot(ph + shift, -2.5 * np.log10(d["kepler"]), color=ps.BLUE, lw=2.0,
                label="Kepler, 2009-2013 mean light curve" if shift == 0 else None)
        ax.plot(ph + shift, -2.5 * np.log10(d["model"]), color=ps.ORANGE, lw=1.4, ls="-",
                label="Roche model with Molnar et al.'s parameters" if shift == 0 else None)
    ax.set_ylim(0.215, -0.045)
    ax.set_xlim(0, 2)
    ax.set_xlabel("orbital phase (0 = heavier star eclipsed; the pattern repeats every orbit)")
    ax.set_ylabel("brightness below maximum (mag)")
    _panel(ax, "b", "it never stops varying: two stars, one peanut-shaped envelope")
    r = _res()["roche"]["kic_model"]
    ax.legend(loc="upper center", ncol=2)
    ps.caption(fig, f"Model: q = {m['q']}, inclination {m['incl']}°, fill-out {m['fillout']}, T = {m['t1']:.0f} / "
                    f"{m['t2']:.0f} K (Molnar et al. 2017); amplitude {r['amp_model_mag']:.3f} mag against Kepler's "
                    f"{r['amp_kepler_mag']:.3f} mag, rms difference {100 * r['rms_diff']:.1f}% (the model has no "
                    "starspots, which make Kepler's two maxima unequal). Images: brightness of the Roche surface with "
                    "gravity darkening (T ∝ g^0.08) and limb darkening. Kepler's mean curve is our 14-harmonic template "
                    "of 59,451 long-cadence points.", y=0.0)
    return ps.save(fig, _out("02_contact_light_curve"))



# ---------------------------------------------------------------------------
# 03  The period distribution and its edge
# ---------------------------------------------------------------------------


def fig03_period_distribution():
    import census as CE
    res = _res()["census"]
    h = np.load(C.DERIVED_DIR / "census_hist.npz")
    e = h["edges"]
    mid = 0.5 * (e[1:] + e[:-1])
    fig = plt.figure(figsize=(11.5, 7.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1.0], height_ratios=[1, 1], hspace=0.38, wspace=0.25)
    ax = fig.add_subplot(gs[:, 0])
    dlog = e[1] - e[0]

    def coarse(n, k=2):
        return n.reshape(-1, k).sum(axis=1)

    ec = e[::2]
    midc = 0.5 * (ec[1:] + ec[:-1])
    dlc = ec[1] - ec[0]

    def norm(n):
        return n / (n.sum() * dlc)

    ax.step(10 ** midc, norm(coarse(h["ztf"])), where="mid", color=ps.BLUE, lw=1.6,
            label=f"ZTF contact (EW) binaries, all {res['n_ew']:,}")
    ax.step(10 ** midc, norm(coarse(h["near"])), where="mid", color=ps.ORANGE, lw=1.3,
            label=f"those within 1 kpc (Gaia parallaxes), {res['edge_within_1kpc']['n']:,}")
    ax.step(10 ** midc, norm(coarse(h["gaia"])), where="mid", color=ps.AQUA, lw=1.2,
            label=f"Gaia DR3 contact-like eclipsing binaries, {res['gaia']['n_contact_like']:,}")
    ax.set_xscale("log")
    ax.set_xlim(0.15, 1.5)
    ax.set_xticks([0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5])
    ax.set_xticklabels(["0.2", "0.25", "0.3", "0.4", "0.5", "0.7", "1", "1.5"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.axvspan(L.RUCINSKI2007["cutoff_lo"], L.RUCINSKI2007["cutoff_hi"], color=ps.GRID, alpha=0.9, lw=0)
    ax.axvline(0.2, color=ps.AQUA, lw=0.8, alpha=0.6)
    ps.annotate(ax, "Gaia's period search\nstops at 0.2 d", (0.2, 0.55), (0.153, 1.15), color=ps.INK_SECONDARY)
    ed = res["edge"]
    ps.annotate(ax, f"the edge: {ed['p_edge']:.3f} d\n(Rucinski 2007: 0.215-0.22)", (ed["p_edge"], 1.75),
                (0.153, 2.6), color=ps.INK)
    ax.set_xlabel("orbital period (days)")
    ax.set_ylabel("fraction of binaries per unit log P")
    _panel(ax, "a", "no contact binary orbits faster than ~5.3 hours")
    ax.legend(loc="upper right", fontsize=8)
    # (b) the edge at every apparent magnitude
    ax = fig.add_subplot(gs[0, 1])
    sl = res["mag_slices"]
    xm = [0.5 * (s["lo"] + s["hi"]) for s in sl]
    pe = [s["p_edge"] for s in sl]
    pe_err = [np.log(10) * s["p_edge"] * s["logp_edge_err"] for s in sl]
    ax.errorbar(xm, pe, yerr=pe_err, xerr=[0.5 * (s["hi"] - s["lo"]) for s in sl], fmt="o", ms=6, color=ps.BLUE,
                mec=ps.SURFACE, mew=1.2, elinewidth=1.0, capsize=0)
    ax.axhspan(L.RUCINSKI2007["cutoff_lo"], L.RUCINSKI2007["cutoff_hi"], color=ps.GRID, lw=0)
    ax.set_ylim(0.205, 0.245)
    ax.set_xlabel("ZTF r magnitude")
    ax.set_ylabel("edge period (d)")
    _panel(ax, "b", "the same edge at every depth")
    ax.text(0.03, 0.06, "a survey limit would move with depth;\nthe edge stays put", transform=ax.transAxes,
            fontsize=8, color=ps.INK_SECONDARY)
    # (c) what lies below the edge: luminosity and colour
    ax = fig.add_subplot(gs[1, 1])
    s = pd.read_parquet(C.DERIVED_DIR / "ztf_gaia_sample.parquet",
                        columns=["type", "blended", "ruwe", "plx_snr", "period", "g", "bp_rp", "av_layer", "dm"])
    s = s[(s["type"] == "EW") & ~s["blended"] & (s["ruwe"] < 1.4) & (s["plx_snr"] > 10)]
    s = s.assign(M_G=s["g"] - 0.789 * s["av_layer"] - s["dm"], c0=s["bp_rp"] - 0.413 * s["av_layer"])
    norm_ = s[s["period"] >= 0.22]
    sub = s[s["period"] < 0.21]
    ax.hist2d(norm_["c0"], norm_["M_G"], bins=[np.linspace(0, 3.2, 120), np.linspace(-1, 11, 120)], cmap=SEQ_BLUE,
              norm=LogNorm(), rasterized=True)
    ax.scatter(sub["c0"], sub["M_G"], s=9, color=ps.ORANGE, edgecolors=ps.SURFACE, linewidths=0.4, zorder=3,
               label=f"P < 0.21 d ({len(sub)})")
    ax.set_xlim(0, 3.2)
    ax.set_ylim(11, -1)
    ax.set_xlabel("colour (BP - RP)$_0$")
    ax.set_ylabel("absolute magnitude $M_G$")
    _panel(ax, "c", "below the edge: red dwarf pairs")
    ax.legend(loc="lower left", fontsize=8)
    ps.caption(fig, "Edge = midpoint of a logistic step fitted to the log P histogram. ZTF: Chen, Wang & Deng (2020), "
                    "amplitude > 0.08 mag. Gaia: DR3 eclipsing-binary candidates (Mowlavi et al. 2023) with an "
                    "ellipsoidal or broad-eclipse model and equal-ish minima. Within 1 kpc the distribution peaks at "
                    f"{res['peak_period_within_1kpc']:.2f} d, close to the edge; the magnitude-limited sample peaks at "
                    f"{res['peak_period_logp_hist']:.2f} d because bright, long-period systems are seen farther away. "
                    "Panel c: blue = all contact binaries with parallax/error > 10, extinction from a dust layer.", y=0.0)
    return ps.save(fig, _out("03_period_distribution"))


# ---------------------------------------------------------------------------
# 04  Period-luminosity relation
# ---------------------------------------------------------------------------


def fig04_period_luminosity():
    res = _res()["pl"]
    w = pd.read_parquet(C.DERIVED_DIR / "pl_sample.parquet")
    ch = pd.read_parquet(C.DERIVED_DIR / "chen2018_dr3.parquet")
    fig = plt.figure(figsize=(11.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 0.9], wspace=0.32)
    ax = fig.add_subplot(gs[0])
    lp = np.linspace(-0.6, -0.22, 10)
    ax.hist2d(w["logp"], w["M_W1"], bins=[np.linspace(-0.6, -0.22, 90), np.linspace(-0.4, 4.2, 110)],
              cmap=SEQ_BLUE, norm=LogNorm(), rasterized=True)
    good = np.isfinite(ch["dm_dr3"]) & (ch["RUWE"] < 1.4)
    m_dr3 = ch["w1_app_chen"] - 0.05 * ch["Av"] - ch["dm_dr3"]
    ax.scatter(ch["logp"][good], m_dr3[good], s=14, color=ps.ORANGE, edgecolors=ps.SURFACE, linewidths=0.6, zorder=3,
               label="Chen et al.'s 153 calibrators, Gaia DR3 distances")
    c18 = L.CHEN2018_PL["W1_mean"]
    ax.plot(lp, c18[0] * lp + c18[2], color=ps.VIOLET, lw=1.6, label="Chen et al. 2018 (TGAS)")
    rc = res["chen2018"]["W1_mean_dr3_ruwe"]
    ax.plot(lp, rc["slope"] * lp + rc["b0"], color=ps.ORANGE, lw=1.6, label="same stars, Gaia DR3 (this work)")
    zf = res["ztf_precise"]["W1_ols"]
    ax.plot(lp, zf["slope"] * lp + zf["b0"], color=ps.BLUE, lw=1.2, ls="--",
            label=f"{zf['n']:,} ZTF contact binaries (this work)")
    ax.set_ylim(4.2, -0.4)
    ax.set_xlim(-0.6, -0.22)
    ax.set_xlabel("log period (days)")
    ax.set_ylabel("absolute W1 magnitude (3.4 µm)")
    _panel(ax, "a", "the period tells the luminosity")
    ax.legend(loc="upper left", fontsize=7.5)
    # (b) the slope against distance: the survey's selection at work
    ax = fig.add_subplot(gs[1])
    b = res["distance_bins"]
    xm = [0.5 * (x["d0"] + x["d1"]) / 1000 for x in b]
    ax.errorbar(xm, [x["slope"] for x in b], yerr=[x["slope_err"] for x in b],
                xerr=[0.5 * (x["d1"] - x["d0"]) / 1000 for x in b], fmt="o", ms=6, color=ps.BLUE, mec=ps.SURFACE,
                mew=1.2, elinewidth=1.0, capsize=0)
    ax.axhspan(rc["slope"] - rc["slope_err"], rc["slope"] + rc["slope_err"], color="#fbe1d6", lw=0)
    ax.axhline(rc["slope"], color=ps.ORANGE, lw=1.2)
    ax.set_xlabel("distance (kpc)")
    ax.set_ylabel("fitted W1 slope")
    _panel(ax, "b", "the survey's fingerprint")
    ax.text(0.30, 0.97, "near: ZTF saturates on the bright,\nlong-period stars. Far: only the\nbright ones have "
                        "good parallaxes.", transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    ax.text(4.0, rc["slope"] + rc["slope_err"] + 0.06, "Gaia DR3 recalibration\nof Chen's stars", ha="right",
            va="bottom", fontsize=7.5, color=ps.ORANGE)
    # (c) scatter: TGAS vs DR3 for the same stars
    ax = fig.add_subplot(gs[2])
    tg = res["chen2018"]["W1_mean_tgas_same_stars"]
    labels = ["TGAS\n(2018)", "Gaia DR3", "ZTF sample"]
    vals = [tg["scatter"], rc["scatter"], zf["scatter"]]
    cols = [ps.VIOLET, ps.ORANGE, ps.BLUE]
    ax.bar(range(3), vals, width=0.55, color=cols)
    for k, v in enumerate(vals):
        ax.text(k, v + 0.004, f"{v:.3f}\n({100 * 0.4605 * v:.1f}%)", ha="center", va="bottom", fontsize=8,
                color=ps.INK)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 0.27)
    ax.set_ylabel("scatter about the relation (mag)")
    _panel(ax, "c", "distance precision")
    ax.grid(axis="x", visible=False)
    ps.caption(fig, "W1 from AllWISE (time-averaged), extinction A_W1 = 0.039 A_V from a dust-layer model. Gaia DR3 "
                    "parallaxes with the Lindegren et al. (2021) zero point; ZTF sample: parallax/error > 10, RUWE < "
                    f"1.4, no bright neighbour, Chen's period range. Fits by ordinary least squares (Chen's method; "
                    f"their Table 2 is reproduced exactly from their Table 1). Percentages in c: the distance error "
                    "of one star.", y=0.0)
    return ps.save(fig, _out("04_period_luminosity"))


# ---------------------------------------------------------------------------
# 05  A distance indicator: test and map
# ---------------------------------------------------------------------------


def fig05_distance_map():
    res = _res()["pl"]
    mp = pd.read_parquet(C.DERIVED_DIR / "pl_map.parquet")
    fig = plt.figure(figsize=(11.5, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.25)
    # (a) PL distance against Gaia distance
    ax = fig.add_subplot(gs[0])
    chk = mp[mp["plx_snr"] > 10]
    ax.hist2d(chk["dist_pc"] / 1000, chk["d_pl_pc"] / 1000, bins=[np.linspace(0, 3.5, 120)] * 2, cmap=SEQ_BLUE,
              norm=LogNorm(), rasterized=True)
    ax.plot([0, 3.5], [0, 3.5], color=ps.INK_SECONDARY, lw=0.9)
    wv = res["wesenheit_vs_parallax"]
    ax.text(0.04, 0.95, f"{wv['n']:,} stars with parallax/error > 10\nmedian ratio {wv['median_ratio']:.3f}, "
                        f"scatter {wv['scatter_percent']:.0f}%", transform=ax.transAxes, va="top", fontsize=8,
            color=ps.INK)
    ax.set_xlabel("distance from Gaia parallax (kpc)")
    ax.set_ylabel("distance from period + Gaia colours (kpc)")
    ax.set_aspect("equal")
    _panel(ax, "a", "the period predicts the distance")
    # (b) the map: face-on Milky Way around the Sun
    ax = fig.add_subplot(gs[1])
    lr, br = np.radians(mp["l"]), np.radians(mp["b"])
    d = mp["d_pl_pc"] / 1000
    x = d * np.cos(br) * np.cos(lr)
    y = d * np.cos(br) * np.sin(lr)
    ax.hist2d(x, y, bins=[np.linspace(-9, 9, 260), np.linspace(-9, 9, 260)], cmap=SEQ_BLUE, norm=LogNorm(),
              rasterized=True)
    ax.plot(0, 0, "*", ms=11, color=ps.ORANGE, mec=ps.SURFACE, mew=0.8)
    ax.annotate("Sun", (0, 0), xytext=(6, 6), textcoords="offset points", fontsize=8, color=ps.INK)
    ax.annotate("", xy=(8.6, -7.6), xytext=(5.6, -7.6),
                arrowprops={"arrowstyle": "->", "color": ps.INK_SECONDARY, "lw": 0.9})
    ax.text(5.6, -7.25, "Galactic centre (8.2 kpc)", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_aspect("equal")
    ax.set_xlim(-9, 9)
    ax.set_ylim(-9, 9)
    ax.set_xlabel("x (kpc, towards l = 0°)")
    ax.set_ylabel("y (kpc, towards l = 90°)")
    _panel(ax, "b", f"{wv['n_mapped']:,} contact binaries placed in the Galaxy")
    ax.text(0.02, 0.98, "ZTF sees the northern sky (dec > -30°): the\nempty wedge is the southern Milky Way. Dust\n"
                        "hides the disc beyond a few kpc in the plane.", transform=ax.transAxes, fontsize=7.5,
            color=ps.INK_SECONDARY, va="top")
    cal = res["wesenheit_calibration"]
    ps.caption(fig, f"Distances from the Wesenheit index W_G = G - 1.90 (BP - RP), which cancels extinction to first "
                    f"order, calibrated on Chen et al.'s stars with Gaia DR3 parallaxes: M(W_G) = "
                    f"{cal['slope']:.2f} log P {cal['b0']:+.2f} (scatter {cal['scatter']:.2f} mag, {cal['n']} stars). "
                    f"Median distance {wv['median_d_pl'] / 1000:.1f} kpc; 90% within {wv['p90_d_pl'] / 1000:.1f} kpc.",
               y=0.0)
    return ps.save(fig, _out("05_distance_map"))



# ---------------------------------------------------------------------------
# 06-07  KIC 9832227: the merger that was not
# ---------------------------------------------------------------------------

#: One colour per data source, in the palette's fixed order.
KIC_SOURCES = [("Kepler", ps.BLUE, "s"), ("TESS", ps.ORANGE, "D"), ("AAVSO", ps.AQUA, "o"),
               ("ASAS-SN", ps.YELLOW, "^"), ("SuperWASP", ps.MAGENTA, "v"), ("SuperWASP (VeSPA)", ps.MAGENTA, "v"),
               ("NSVS", ps.GREEN, "o")]


def _kic_models():
    import kic9832227 as K
    res = _res()["kic"]
    df = pd.read_csv(C.DERIVED_DIR / "kic9832227_oc.csv")
    p = res["fits_all"]["quad+ltte"]["params"]

    def ltte(t):
        y = (np.asarray(t) - K.T_REF) / C.YEAR_D
        return (p[0] + p[1] * y + p[2] * y * y + p[3] * np.sin(2 * np.pi * (y - p[5]) / p[4])) * 1440.0

    # Molnar et al.'s published law, its zero point fitted to the data they had (with the NSVS time as labelled)
    sel = df[(df["year"] <= 2016.75) & df["source"].isin([s for s in K.MOLNAR_SOURCES if s != "NSVS"] +
                                                          ["NSVS (time as labelled)"])]
    sel = sel[~((sel["source"] == "SuperWASP") & (sel["year"] < 2006))]
    base = K.molnar_model(sel["t"].to_numpy(), 0.0)
    w = 1 / (sel["err_min"].to_numpy() ** 2 + 1.0)
    off = np.sum(w * (sel["oc_min"].to_numpy() - base)) / np.sum(w)

    def molnar(t):
        return K.molnar_model(np.asarray(t), 0.0) + off

    return df, ltte, molnar, K.kovacs_model, res


def fig06_kic_oc():
    import kic9832227 as K
    df, ltte, molnar, kovacs, res = _kic_models()
    fig = plt.figure(figsize=(11.5, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    t_grid = K.jd_of(np.linspace(1998.5, 2026.5, 2000))
    yr_grid = K.year_of(t_grid)
    ax.plot(yr_grid, ltte(t_grid), color=ps.INK, lw=1.3, zorder=2,
            label=f"steady period decrease + {res['fits_all']['quad+ltte']['period_yr']:.0f}-yr wave (fit, this work)")
    m_t = t_grid[yr_grid < 2022.15]
    mm = molnar(m_t)
    era = K.year_of(m_t) <= 2016.75
    ax.plot(K.year_of(m_t)[era], mm[era], color=ps.VIOLET, lw=1.6, zorder=2,
            label="Molnar et al. 2017: exponential decay, merger 2022.2")
    ax.plot(K.year_of(m_t)[~era], mm[~era], color=ps.VIOLET, lw=1.6, ls=(0, (4, 2)), zorder=2)
    kv = yr_grid >= 2018.7
    ax.plot(yr_grid[~kv], kovacs(t_grid[~kv]), color=ps.INK_SECONDARY, lw=1.0, zorder=2)
    ax.plot(yr_grid[kv], kovacs(t_grid[kv]), color=ps.INK_SECONDARY, lw=1.0, ls=(0, (1.5, 1.5)), zorder=2,
            label="Kovács et al. 2019: parabola + 13.5-yr wave")
    ax.axvspan(2022.2 - 0.7, 2022.2 + 0.7, color="#ece9f7", lw=0, zorder=0)
    ax.text(2022.2, 60, "predicted\nmerger", ha="center", va="top", fontsize=8, color=ps.VIOLET)
    lit = df[df["kind"] == "literature"]
    ax.errorbar(lit["year"], lit["oc_min"], yerr=lit["err_min"], fmt="o", ms=4.5, mfc=ps.SURFACE, mec=ps.MUTED,
                ecolor=ps.MUTED, elinewidth=0.8, label="other surveys (Kovács et al. 2019 compilation)", zorder=3)
    shown = set()
    for src_name, col, mk in KIC_SOURCES:
        d = df[df["source"] == src_name]
        lab = src_name.replace(" (VeSPA)", "")
        ax.errorbar(d["year"], d["oc_min"], yerr=d["err_min"], fmt=mk, ms=5.5 if mk != "D" else 5, color=col,
                    mec=ps.SURFACE, mew=0.8, elinewidth=0.9, capsize=0, zorder=4,
                    label=None if lab in shown else f"{lab} (our timing)")
        shown.add(lab)
    bad = df[df["source"] == "NSVS (time as labelled)"].iloc[0]
    ax.plot(bad["year"], bad["oc_min"], "X", ms=10, color=ps.RED, mec=ps.SURFACE, mew=0.8, zorder=5,
            label="NSVS with the archive's time label taken literally")
    good = df[df["source"] == "NSVS"].iloc[0]
    ax.text(bad["year"] + 0.45, 0.5 * (bad["oc_min"] + good["oc_min"]), "half a day minus\none orbit = 60.5 min",
            fontsize=8, color=ps.RED, va="center")
    ax.annotate("", xy=(bad["year"], bad["oc_min"] + 3), xytext=(good["year"], good["oc_min"] - 3),
                arrowprops={"arrowstyle": "->", "color": ps.RED, "lw": 1.0})
    ax.set_xlim(1998.5, 2026.5)
    ax.set_ylim(-245, 65)
    ax.set_ylabel("O - C: eclipse time minus a constant-period prediction (minutes)")
    ax.tick_params(labelbottom=False)
    _panel(ax, "a", "KIC 9832227: 26 years of eclipse timing, and the merger that did not happen")
    ax.legend(loc="lower left", fontsize=7.5, ncol=2)
    # residuals from our model
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    ours = df[df["source"] != "NSVS (time as labelled)"]
    ax2.axhline(0, color=ps.INK, lw=1.0)
    ax2.errorbar(lit["year"], lit["oc_min"] - ltte(lit["t"]), yerr=lit["err_min"], fmt="o", ms=4, mfc=ps.SURFACE,
                 mec=ps.MUTED, ecolor=ps.MUTED, elinewidth=0.8)
    for src_name, col, mk in KIC_SOURCES:
        d = ours[ours["source"] == src_name]
        ax2.errorbar(d["year"], d["oc_min"] - ltte(d["t"]), yerr=d["err_min"], fmt=mk, ms=4.5, color=col,
                     mec=ps.SURFACE, mew=0.6, elinewidth=0.8, capsize=0)
    late = ours[ours["year"] > 2017]
    mres = late["oc_min"] - molnar(np.minimum(late["t"], K.jd_of(2022.1)))
    ax2.set_ylim(-11, 11)
    ax2.set_ylabel("residual (min)")
    ax2.set_xlabel("year")
    jit = res["fits_all"]["quad+ltte"]["jitter_min"]
    ax2.text(0.01, 0.92, f"scatter beyond the error bars: {jit:.1f} min (starspots). Molnar's curve is off this "
                         f"panel by 2018 and diverges at 2022.2; Kovács et al.'s misses 2024 by "
                         f"{abs(res['kovacs_prediction']['miss_2024_min']):.0f} min.", transform=ax2.transAxes,
             fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    ps.caption(fig, "O - C against Socia et al.'s reference ephemeris (BJD 2454953.949183 + 0.4579515 E), so values "
                    "compare one-to-one with Kovács et al. (2019). Our timings: one Kepler-based template, shifted, "
                    "scaled and offset per season, survey camera or observer; all times converted to BJD(TDB). The "
                    "wave's period is barely constrained (1.5 cycles in 26 years).", y=0.02)
    return ps.save(fig, _out("06_kic9832227_oc"))


def fig07_kic_mistake():
    import kic9832227 as K
    import timing as T
    import pickle
    res = _res()["kic"]
    tpl = pickle.loads((C.DERIVED_DIR / "kic9832227_template.pkl").read_bytes())
    fig = plt.figure(figsize=(11.5, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.25, 0.95], wspace=0.32)
    # (a) the NSVS light curve, folded both ways, with the template at each fitted shift
    ax = fig.add_subplot(gs[0])
    ph_t = np.linspace(-0.25, 1.6, 500)
    df_oc = pd.read_csv(C.DERIVED_DIR / "kic9832227_oc.csv")
    for conv, src_name, col, lab, dy in (("correct", "NSVS", ps.GREEN, "time = MJD + 2400000.5", 0.0),
                                         ("as_labelled", "NSVS (time as labelled)", ps.RED,
                                          "time = 'JD-2400000' + 2400000", 0.45)):
        ns = K.load_nsvs(conv)
        ph = np.mod((ns["t"] - K.T_REF) / K.P_REF, 1.0)
        fl = ns["flux"] / np.median(ns["flux"])
        mag = -2.5 * np.log10(fl) + dy
        ax.plot(ph, mag, ".", ms=2.5, color=col, alpha=0.55)
        ax.plot(ph + 1, mag, ".", ms=2.5, color=col, alpha=0.2)
        oc = df_oc.loc[df_oc["source"] == src_name, "oc_min"].iloc[0]
        shift = oc / 1440.0 / K.P_REF
        model = -2.5 * np.log10(tpl(ph_t - shift) / np.median(tpl(ph_t))) + dy
        ax.plot(ph_t, model, color=ps.INK, lw=1.0)
        for k in (0, 1):
            ax.plot([k + shift] * 2, [dy + 0.17, dy + 0.27], color=col, lw=1.6)
        ax.text(-0.22, dy - 0.27, lab, fontsize=7.5, color=col)
        ax.text(1.0 + shift + 0.03, dy + 0.27, f"{oc:+.0f} min", fontsize=7.5, color=col, va="center")
    ax.axvline(0, color=ps.AXIS, lw=0.8)
    ax.axvline(1, color=ps.AXIS, lw=0.8)
    ax.set_xlim(-0.25, 1.6)
    ax.set_ylim(0.8, -0.38)
    ax.set_xlabel("phase predicted by the Kepler-era ephemeris")
    ax.set_ylabel("NSVS brightness, 1999-2000 (mag, offset)")
    _panel(ax, "a", "same photometry, two clocks")
    # (b) the scatter each trial merger date needs
    ax = fig.add_subplot(gs[1])
    for key, col, lab in (("molnar_like_labelled", ps.RED, "with the mislabelled 1999 time"),
                          ("molnar_like_correct", ps.GREEN, "with the correct 1999 time")):
        p = res["exponential_profiles"][key]
        g = np.array(p["grid"])
        jg = np.array(p["jitter_grid"])
        k = (g < 2060) & np.isfinite(jg)
        ax.plot(g[k], jg[k], color=col, lw=1.6, label=f"merger law, {lab}")
        if key == "molnar_like_correct":
            ax.axhline(p["ltte_jitter_min"], color=col, lw=0.9, ls=(0, (1.5, 1.5)),
                       label="smoothly changing period, correct time")
    ax.axvspan(2022.2 - 0.7, 2022.2 + 0.7, color="#ece9f7", lw=0)
    ax.text(2023.3, 5.35, "Molnar et al.\n2022.2 ± 0.7", ha="left", va="top", fontsize=7.5, color=ps.VIOLET)
    pl = res["exponential_profiles"]["molnar_like_labelled"]
    pc = res["exponential_profiles"]["molnar_like_correct"]
    ax.set_ylim(0, 7.2)
    ax.set_xlim(2016.5, 2060)
    ax.set_xlabel("merger date t0 of the exponential law (year)")
    ax.set_ylabel("scatter the fit needs (minutes)")
    _panel(ax, "b", "a merger date made by one number")
    ax.text(0.97, 0.40, f"best: t0 = {pl['t0_best']:.1f} (+{pl['t0_hi'] - pl['t0_best']:.1f}/"
                        f"-{pl['t0_best'] - pl['t0_lo']:.1f}) with the bad time;\nwith the right one no date "
                        f"fits better than {pc['jitter_min']:.1f} min,\nwhile a smoothly changing period needs "
                        f"{pc['ltte_jitter_min']:.1f} min.", transform=ax.transAxes, ha="right", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ax.legend(loc="upper right", fontsize=7.5)
    # (c) the bookkeeping
    ax = fig.add_subplot(gs[2])
    ax.axis("off")
    _panel(ax, "c", "how half a day disappears")
    lines = [
        ("NSVS frame table, column 'JD':", ps.INK, True),
        ("labelled 'Julian date (JD-2400000)'", ps.INK_SECONDARY, False),
        ("holds MJD = JD - 2400000.5", ps.INK_SECONDARY, False),
        ("(1999-04-06T08:53:47 UT -> 51274.37069)", ps.MUTED, False),
        ("", ps.INK, False),
        ("Molnar et al. 2017 use", ps.INK, True),
        ("'MJD = JD - 2400000.0' (Kepler EB catalogue)", ps.INK_SECONDARY, False),
        ("", ps.INK, False),
        ("half a day = 1.09 orbits of 0.458 d", ps.INK, True),
        ("counted as 1 orbit, the rest is", ps.INK_SECONDARY, False),
        (f"0.5 - 0.45795 d = {res['nsvs']['expected_shift_min']:.1f} min early", ps.RED, False),
        (f"measured shift: {res['nsvs']['shift_min']:.2f} min", ps.INK_SECONDARY, False),
        ("", ps.INK, False),
        ("Socia et al. (2018) found the error;", ps.INK_SECONDARY, False),
        ("our NSVS time agrees with theirs", ps.INK_SECONDARY, False),
        (f"(O - C {res['nsvs']['oc_min']:+.1f} ± {res['nsvs']['err_min']:.1f} min; Socia +42, Kovács +38)",
         ps.INK_SECONDARY, False),
    ]
    y = 0.92
    for text, col, bold in lines:
        ax.text(0.0, y, text, transform=ax.transAxes, fontsize=8, color=col, fontweight="bold" if bold else "normal")
        y -= 0.058
    ps.caption(fig, "a: NSVS unfiltered photometry of both survey IDs (5597755, 5620022; Wozniak et al. 2004, VizieR "
                    "II/287), folded on the reference ephemeris; black: the Kepler template at each version's fitted shift. "
                    "b: the exponential period law of Tylenda et al. (2011) refitted for each trial merger date to "
                    "the data sets Molnar et al. had (2014-2016 ASAS-SN and AAVSO seasons standing in for their own "
                    "photometry); 'scatter' is the jitter added to the error bars in the maximum-likelihood fit.",
               y=0.0)
    return ps.save(fig, _out("07_kic9832227_mistake"))



# ---------------------------------------------------------------------------
# 08-10  V1309 Sco: a merger watched
# ---------------------------------------------------------------------------


def fig08_v1309_lightcurve():
    import v1309sco as V
    res = _res()["v1309"]
    d = V.load_ogle()
    I = d["I"]
    yr = V.jd_to_year(I["t"])
    fig = plt.figure(figsize=(11.5, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.22)
    ax = fig.add_subplot(gs[0])
    ax.plot(yr, I["mag"], ".", ms=3, color=ps.BLUE)
    on = res["onset"]
    t_disc = V.jd_to_year(L.TYLENDA2011["discovery_jd"])
    ax.axvline(on["t_on10_year"], color=ps.INK_SECONDARY, lw=0.8)
    ax.axvline(t_disc, color=ps.INK_SECONDARY, lw=0.8)
    ax.text(on["t_on10_year"] - 0.07, 9.0, f"the rise\nbegins\n{on['t_on10_year']:.2f}\n(late March)",
            ha="right", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.text(t_disc + 0.07, 7.6, "discovered as\nNova Sco 2008\n(Sep 2)", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.axhline(L.TYLENDA2011["peak_i"], xmin=0, xmax=1, color=ps.AXIS, lw=0.8)
    ax.text(2001.6, L.TYLENDA2011["peak_i"] + 0.35, "peak I ≈ 6.8 (AAVSO); OGLE's camera saturated above I ≈ 11",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.text(2009.6, 13.4, "decline:\nthe merged\nremnant", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.invert_yaxis()
    ax.set_ylim(17.4, 6.2)
    ax.set_xlim(2001.4, 2010.9)
    ax.set_xlabel("year")
    ax.set_ylabel("OGLE I magnitude")
    _panel(ax, "a", "V1309 Sco: ten magnitudes in half a year")
    # (b) the progenitor, zoomed
    ax = fig.add_subplot(gs[1])
    pre = I["t"] < 2454600
    ax.plot(yr[pre], I["mag"][pre], ".", ms=2.5, color=ps.BLUE)
    ax.invert_yaxis()
    ax.set_xlim(2001.5, 2008.45)
    ax.set_ylim(17.0, 15.35)
    ax.set_xlabel("year")
    ax.set_ylabel("OGLE I magnitude")
    _panel(ax, "b", "the progenitor: a 1.4-day binary")
    ax.text(0.03, 0.05, "the 0.3-0.5 mag spread each season is the\nbinary's light curve, sampled once a night",
            transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY)
    ps.annotate(ax, "brightest: April 2007", (2007.3, 15.5), (2004.4, 15.5), color=ps.INK_SECONDARY)
    ps.annotate(ax, "a year of fading,\nthen the rise", (2008.25, 16.55), (2005.7, 16.72), color=ps.INK_SECONDARY)
    ps.caption(fig, f"OGLE-III (2001-2009) and OGLE-IV (2010) photometry published with Tylenda et al. (2011): "
                    f"{res['n_I']:,} I-band points, {res['n_pre_discovery']:,} before the discovery. Start of the "
                    f"rise: when a fitted exponential brightening (e-folding {res['onset']['tau_d']:.0f} days) "
                    f"reaches 10% of the progenitor's flux.", y=0.0)
    return ps.save(fig, _out("08_v1309_lightcurve"))


def fig09_v1309_period():
    import v1309sco as V
    res = _res()["v1309"]
    per = res["periods"]
    t = np.array([p["t"] for p in per])
    p = np.array([p["period"] for p in per])
    e = np.array([max(x["err_boot"], x["err_curv"]) for x in per])
    ty = res["vs_tylenda"]
    fig = plt.figure(figsize=(11.5, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.24)
    ax = fig.add_subplot(gs[0])
    tg = np.linspace(2452200, 2455200, 1500)
    lp = res["law_from_periods"]
    lph = res["law_from_phases"]
    ax.plot(V.jd_to_year(tg), V.exp_period(tg, L.TYLENDA2011["p0"], L.TYLENDA2011["c"], L.TYLENDA2011["t0"]),
            color=ps.VIOLET, lw=1.4, label=f"Tylenda et al. 2011: t0 = {V.jd_to_year(L.TYLENDA2011['t0']):.2f}")
    ax.plot(V.jd_to_year(tg), V.exp_period(tg, lp["p0"], lp["c"], lp["t0"]), color=ps.BLUE, lw=1.4,
            label=f"our periods: t0 = {lp['t0_year']:.2f} ± {lp['t0_year_err']:.2f}")
    ax.plot(V.jd_to_year(tg), V.exp_period(tg, lph["p0"], lph["c"], lph["t0"]), color=ps.AQUA, lw=1.4,
            label=f"our eclipse timing: t0 = {lph['t0_year']:.2f} ± {lph['t0_year_err']:.2f}")
    ax.errorbar(V.jd_to_year(np.array(ty["t"])) + 0.03, ty["periods"], yerr=ty["err"], fmt="o", ms=4.5,
                mfc=ps.SURFACE, mec=ps.VIOLET, ecolor=ps.VIOLET, elinewidth=0.8, capsize=0,
                label="Tylenda et al.'s 17 periods (read from their figure)")
    ax.errorbar(V.jd_to_year(t), p, yerr=e, fmt="o", ms=5, color=ps.BLUE, mec=ps.SURFACE, mew=0.8,
                elinewidth=1.0, capsize=0, label="our 17 periods (bootstrap errors)")
    ax.axvline(res["onset"]["t_on10_year"], color=ps.INK_SECONDARY, lw=0.8)
    ax.text(res["onset"]["t_on10_year"] + 0.04, 1.4475, "rise\nbegins", fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    ax.set_xlim(2002.0, 2009.0)
    ax.set_ylim(1.4145, 1.4485)
    ax.set_xlabel("year")
    ax.set_ylabel("orbital period (days)")
    _panel(ax, "a", f"the orbit shrank by {res['period_drop_percent']:.1f}% in six years")
    ax.legend(loc="lower left", fontsize=7.5)
    # (b) O - C against a constant period, with the exponential law
    ax = fig.add_subplot(gs[1])
    oc = res["oc"]
    toc = np.array([o["t"] for o in oc])
    dph = np.array([o["dphase"] for o in oc])
    sdp = np.array([o["dphase_err"] for o in oc])
    a2 = np.array([o["a2"] for o in oc])
    p0, c, t0, z = lph["p0"], lph["c"], lph["t0"], lph["zero"]
    law0 = (lp["p0"], lp["c"], lp["t0"])
    t_ref = 2453500.0
    n_obs = V.cycles(toc, *law0, t_ref) + dph - z
    p_lin = 1.4329
    oc_lin = (toc - t_ref) - p_lin * np.round(n_obs) - p_lin * (n_obs - np.round(n_obs))
    # observed eclipse times: t_obs where the law's cycle count equals the observed count
    tg2 = np.linspace(2452300, 2454550, 3000)
    n_law = V.cycles(tg2, p0, c, t0, t_ref)
    t_ecl = np.interp(n_obs, n_law, tg2)
    e_lin = np.round((t_ecl - t_ref) / p_lin)
    ocd = t_ecl - (t_ref + p_lin * e_lin)
    tt = np.linspace(2452300, 2454550, 600)
    nn = V.cycles(tt, p0, c, t0, t_ref)
    ax.plot(V.jd_to_year(tt), (tt - (t_ref + p_lin * nn)) * 1.0, color=ps.AQUA, lw=1.4, label="exponential law")
    good = a2 > 0.03
    ax.errorbar(V.jd_to_year(toc[good]), (toc[good] - (t_ref + p_lin * n_obs[good])), yerr=sdp[good] * p_lin,
                fmt="o", ms=5, color=ps.BLUE, mec=ps.SURFACE, mew=0.8, elinewidth=1.0, capsize=0,
                label="eclipse timing, each subsample")
    if (~good).any():
        ax.errorbar(V.jd_to_year(toc[~good]), (toc[~good] - (t_ref + p_lin * n_obs[~good])),
                    yerr=sdp[~good] * p_lin, fmt="o", ms=5, mfc=ps.SURFACE, mec=ps.BLUE, ecolor=ps.BLUE,
                    elinewidth=1.0, capsize=0, label="eclipses too shallow to time (not fitted)")
    ax.set_xlabel("year")
    ax.set_ylabel(f"O - C against a constant {p_lin} d period (days)")
    _panel(ax, "b", "eclipses arriving ever earlier")
    ax.legend(loc="lower left", fontsize=7.5)
    ps.caption(fig, f"Periods: our multi-harmonic periodogram (Fourier series + linear trend) per subsample, split as in "
                    f"Tylenda et al.'s Fig. 2 (median difference {abs(ty['median_diff']) * 86400:.2f} s). Law: P = P0 exp"
                    f"(c / (t - t0)), formally reaching zero period at t0. Eclipse timing: phase of the second Fourier "
                    f"harmonic in each subsample; formal errors are scaled by the fit's chi^2 "
                    f"({lph['chi2']:.0f} for {lph['dof']} dof) because the light curve itself was changing.", y=0.0)
    return ps.save(fig, _out("09_v1309_period"))


def fig10_v1309_shape():
    res = _res()["v1309"]
    cur = list(np.load(C.DERIVED_DIR / "v1309_curves.npz", allow_pickle=True)["curves"])
    fig = plt.figure(figsize=(11.5, 6.4))
    n = len(cur)
    gs = fig.add_gridspec(2, n, height_ratios=[1.25, 1.0], hspace=0.42, wspace=0.08)
    for k, c in enumerate(cur):
        ax = fig.add_subplot(gs[0, k])
        ph = np.array(c["phase"])
        dm = np.array(c["dmag"])
        ax.plot(ph, dm, ".", ms=2.2, color=ps.MUTED)
        ax.plot(ph + 1, dm, ".", ms=2.2, color=ps.MUTED, alpha=0.4)
        fp, fm = np.array(c["fit_phase"]), np.array(c["fit_dmag"])
        ax.plot(np.concatenate([fp, fp + 1]), np.concatenate([fm, fm]), color=ps.BLUE, lw=1.6)
        ax.text(0.04, 0.04, f"{c['n']} pts", transform=ax.transAxes, fontsize=7, color=ps.MUTED)
        ax.set_ylim(0.45, -0.45)
        ax.set_xlim(0, 2)
        ax.set_xticks([0, 1, 2])
        if k:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("I minus season median (mag)")
        ax.set_title(f"{int(c['year'])}" if c["year"] < 2008 else "2008\n(Jan-Mar)", loc="center", fontsize=9,
                     fontweight="bold", color=ps.INK)
        ax.set_xlabel("phase", fontsize=8)
    fig.text(0.125, 0.965, "a  from two equal humps to one: the progenitor's light curve, season by season", fontsize=11,
             fontweight="bold", color=ps.INK)
    ax = fig.add_subplot(gs[1, :])
    sh = res["shape"]
    import v1309sco as V
    ts = V.jd_to_year(np.array([s["t"] for s in sh]))
    ax.plot(ts, [s["a2"] for s in sh], "o-", ms=5, color=ps.BLUE, mec=ps.SURFACE, mew=0.8, lw=1.2,
            label="2nd harmonic: two eclipses per orbit")
    ax.plot(ts, [s["a1"] for s in sh], "s-", ms=5, color=ps.ORANGE, mec=ps.SURFACE, mew=0.8, lw=1.2,
            label="1st harmonic: one hump per orbit")
    ax.set_xlabel("year")
    ax.set_ylabel("semi-amplitude (mag)")
    ax.set_ylim(0, 0.22)
    _panel(ax, "b", "the eclipses fade while one side of the binary brightens")
    ax.legend(loc="upper left", fontsize=8)
    ps.caption(fig, "Each season folded on our exponential ephemeris (phase 0 = deeper minimum), a linear trend in "
                    "time removed; blue: four-harmonic Fourier fit. Amplitudes from four-harmonic fits to each of the "
                    "17 subsamples. A contact binary "
                    "shows two equal maxima; from 2006 the maximum after the deeper minimum grew, then the eclipses "
                    "all but vanished in 2007-2008, as matter flowing out of the system (through L2) began to "
                    "shroud it (Tylenda et al. 2011; Pejcha et al. 2017).", y=0.0)
    return ps.save(fig, _out("10_v1309_shape"))



# ---------------------------------------------------------------------------
# 11-14  T CrB
# ---------------------------------------------------------------------------


def _tcrb():
    import tcrb as TC
    return TC, TC.load_aavso()


def _year(jd):
    return 2000.0 + (np.asarray(jd) - C.JD_J2000) / C.YEAR_D


def fig11_tcrb_history():
    TC, df = _tcrb()
    res = _res()["tcrb"]
    fig = plt.figure(figsize=(11.5, 7.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.32)
    ax = fig.add_subplot(gs[0])
    vis = TC.binned(TC.clean_daily(df, "Vis."), 30.0, min_n=3)
    ax.plot(vis["year"], vis["mag"], ".", ms=2.5, color=ps.BLUE, label="visual estimates, 30-day medians")
    for yr in (1866, 1946):
        er = TC.eruption_curve(df, yr, before=5, after=60)["Vis."]
        er = er[~er["rejected"]]
        ax.plot(_year(er["jd"]), er["mag"], ".", ms=3, color=ps.ORANGE)
        ax.text(_year(L.SCHAEFER2023["peak_jd"][yr]) + 1.2, 3.2, f"{yr}\neruption", fontsize=8, color=ps.ORANGE)
    ax.set_ylim(11.3, 1.6)
    ax.set_xlim(1862, 2028)
    ax.set_xlabel("year")
    ax.set_ylabel("visual magnitude")
    _panel(ax, "a", f"T CrB: 160 years from the AAVSO International Database ({res['n_obs']:,} observations)")
    ax.text(1868, 10.9, "AID's 1866 record starts ten days after the peak (V = 2.0 by Birmingham, May 12); "
                        "it holds almost nothing from 1876 to 1934", fontsize=7.5, color=ps.INK_SECONDARY)
    # (b) the last two decades in B and V
    ax = fig.add_subplot(gs[1])
    for band, col, lab in (("B", ps.BLUE, "B (blue: the accretion disc)"), ("V", ps.ORANGE, "V (mostly the red giant)")):
        m = TC.seasonal_levels(df[df["jd"] > 2451500], band, 15.0)
        ax.plot(m["year"], m["mag"], ".", ms=3.5, color=col, label=lab)
    lv = res["B_levels"]
    ax.axhline(lv["quiescence_2005_14"], color=ps.BLUE, lw=0.8, ls=(0, (1.5, 1.5)))
    ax.text(2003.2, lv["quiescence_2005_14"] + 0.08, "B in quiescence", fontsize=7.5, color=ps.BLUE, va="top")
    for x0, x1, lab in ((2015.3, 2023.2, "super-active (high) state"), (2023.3, 2024.3, "fading"),
                        (2024.4, 2025.9, "recovery"), (2026.0, 2026.8, "decline")):
        ax.plot([x0, x1], [9.12, 9.12], color=ps.INK_SECONDARY, lw=1.0)
        ax.text(0.5 * (x0 + x1), 9.07, lab, fontsize=7.5, color=ps.INK_SECONDARY, ha="center", va="bottom")
    last = res["last30"]
    ax.annotate(f"today: V = {last['V']['median']:.1f}, B = {last['B']['median']:.1f}\n(no eruption)",
                (_year(res["last_jd"]), last["B"]["median"]), xytext=(2019.3, 12.35), fontsize=8, color=ps.INK,
                arrowprops={"arrowstyle": "-", "color": ps.INK, "lw": 0.8})
    ax.set_ylim(12.45, 8.85)
    ax.set_xlim(2003, 2027.2)
    ax.set_xlabel("year")
    ax.set_ylabel("magnitude (15-day medians)")
    _panel(ax, "b", "the decade before the next eruption")
    ax.legend(loc="lower left", fontsize=8)
    ps.caption(fig, "Visual estimates and Johnson B/V (CCD and photoelectric) from the AAVSO International Database "
                    "through the public VSX interface; observations flagged discrepant by AAVSO, and estimates more "
                    "than 1 mag from the same day's consensus of at least three observers, are excluded. The scatter "
                    "within a year is the red giant's ellipsoidal variation (Figure 14) and flickering.", y=0.0)
    return ps.save(fig, _out("11_tcrb_history"))


def fig12_tcrb_eruptions():
    TC, df = _tcrb()
    res = _res()["tcrb"]["eruptions"]
    fig = plt.figure(figsize=(11.5, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.22)
    for k, (lo, hi, title) in enumerate(((-3, 30, "the eruption: down 3 magnitudes in 5 days"),
                                         (20, 320, "and a second, fainter brightening months later"))):
        ax = fig.add_subplot(gs[k])
        for yr, col, mk in ((1946, ps.BLUE, "o"), (1866, ps.ORANGE, "s")):
            ec = TC.eruption_curve(df, yr, before=10, after=400)["Vis."]
            ec = ec[~ec["rejected"] & (ec["days"] > lo) & (ec["days"] < hi)]
            if k == 1:
                b = TC.binned(ec.assign(rejected=False), 3.0, min_n=2)
                b["days"] = b["jd"] - L.SCHAEFER2023["peak_jd"][yr]
                ax.plot(b["days"], b["mag"], mk, ms=4, color=col, mec=ps.SURFACE, mew=0.5,
                        label=f"{yr} ({res[str(yr)]['n_vis']:,} visual estimates)")
            else:
                ax.plot(ec["days"], ec["mag"], mk, ms=3.5, color=col, mec=ps.SURFACE, mew=0.4, alpha=0.9,
                        label=f"{yr}")
        if k == 0:
            ax.plot(0, L.SCHAEFER2023["peak_v"], "*", ms=12, color=ps.VIOLET, mec=ps.SURFACE)
            ax.text(0.6, 2.0, "peak V = 2.0 (Schaefer 2023;\nbefore AID's records)", fontsize=7.5, color=ps.VIOLET,
                    va="center")
            r = res["1946"]
            for drop, t_ in ((2.0, r["t2"]), (3.0, r["t3"])):
                ax.plot([t_, t_], [2.0 + drop - 0.5, 2.0 + drop + 0.5], color=ps.INK, lw=1.0)
                ax.text(t_ + 0.4, 2.0 + drop - 0.35, f"t{int(drop)} = {t_:.1f} d", fontsize=7.5, color=ps.INK)
            ax.set_ylim(11, 1.5)
            ax.set_xlabel("days after the peak")
        else:
            for yr, col in ((1946, ps.BLUE), (1866, ps.ORANGE)):
                r = res[str(yr)]
                ax.annotate(f"{yr}: day {r['secondary_peak_day']:.0f}", (r["secondary_peak_day"], r["secondary_peak_mag"]),
                            xytext=(r["secondary_peak_day"] + (30 if yr == 1946 else -105),
                                    r["secondary_peak_mag"] - (0.25 if yr == 1946 else 0.12)),
                            fontsize=8, color=col, arrowprops={"arrowstyle": "-", "color": col, "lw": 0.8})
            ax.set_ylim(10.6, 7.3)
            ax.set_xlabel("days after the peak (3-day medians)")
        ax.set_ylabel("visual magnitude")
        _panel(ax, "ab"[k], title)
        ax.legend(loc="lower right" if k == 0 else "lower left", fontsize=8)
    ps.caption(fig, f"Days counted from Schaefer's (2023) times of peak light, JD 2402734.5 (1866) and 2431860.0 (1946). "
                    f"t2 and t3: the 1946 visual light curve crosses 4.0 and 5.0 mag (Schaefer: 3.0 and 5.0 days). "
                    f"The two eruptions eighty years apart are indistinguishable; the secondary maximum (Schaefer: "
                    f"~160 days, unique among novae) repeats too. One misdated 1946 estimate is excluded "
                    f"(J. Meek, 5.1 mag on Feb 5, four days before the rise; see tcrb.SUSPECT_RECORDS).", y=0.0)
    return ps.save(fig, _out("12_tcrb_eruptions"))


def fig13_tcrb_before():
    TC, df = _tcrb()
    res = _res()["tcrb"]
    fig = plt.figure(figsize=(11.5, 6.6))
    gs = fig.add_gridspec(2, 1, hspace=0.42)
    # (a) the 1930s and 40s: almost all of AID's record is Leslie Peltier's
    ax = fig.add_subplot(gs[0])
    old = TC.clean_daily(df[(df["jd"] > 2426700) & (df["jd"] < 2431859.5)], "Vis.")
    old = old[~old["rejected"]]
    pel = old[old["obscode"] == "P"]
    oth = old[old["obscode"] != "P"]
    rng = np.random.default_rng(1)
    jit = rng.uniform(-0.03, 0.03, len(pel))
    ax.plot(pel["year"], pel["mag"] + jit, "o", ms=3.2, color=ps.ORANGE, mec=ps.SURFACE, mew=0.3,
            label=f"Leslie Peltier ({len(pel)} estimates; small vertical jitter added)")
    ax.plot(oth["year"], oth["mag"], "s", ms=3.0, color=ps.BLUE, mec=ps.SURFACE, mew=0.3,
            label=f"all other observers ({len(oth)})")
    ax.axvline(_year(L.SCHAEFER2023["peak_jd"][1946]), color=ps.INK, lw=1.0)
    ax.text(_year(L.SCHAEFER2023["peak_jd"][1946]) - 0.05, 8.75, "eruption\n1946 Feb 9", fontsize=7.5, color=ps.INK,
            ha="right", va="top")
    ax.annotate("1945: 10.2, the fading Peltier noticed", (1945.35, 10.22), xytext=(1938.4, 10.55), fontsize=8,
                color=ps.ORANGE, arrowprops={"arrowstyle": "-", "color": ps.ORANGE, "lw": 0.8})
    ax.set_xlim(1932, 1946.6)
    ax.set_ylim(10.7, 8.6)
    ax.set_xlabel("year")
    ax.set_ylabel("visual magnitude")
    _panel(ax, "a", "before 1946: one observer's record")
    ax.legend(loc="upper left", fontsize=7.5)
    # (b) the last fourteen years
    ax = fig.add_subplot(gs[1])
    new = TC.binned(TC.clean_daily(df[df["jd"] > 2455900], "Vis."), 30.0, min_n=5)
    ax.plot(new["year"], new["mag"], "o", ms=3.5, color=ps.BLUE, mec=ps.SURFACE, mew=0.4,
            label="visual estimates, 30-day medians (100-200 observers a year)")
    vl = res["visual_levels"]
    ax.axhline(vl["quiescence_2005_14"], color=ps.MUTED, lw=0.9, ls=(0, (1.5, 1.5)))
    ax.text(2012.1, vl["quiescence_2005_14"] + 0.03, "quiescence 2005-14", fontsize=7.5, color=ps.MUTED, va="top")
    ax.axvline(_year(res["last_jd"]), color=ps.INK, lw=1.0)
    ax.text(_year(res["last_jd"]) - 0.05, 9.25, "today", fontsize=7.5, color=ps.INK, ha="right", va="top")
    ax.set_xlim(2012, 2026.95)
    ax.set_ylim(10.7, 9.2)
    ax.set_xlabel("year")
    ax.set_ylabel("visual magnitude")
    _panel(ax, "b", "eighty years later: high state, fading, recovery, a new decline, no eruption yet")
    ax.legend(loc="lower left", fontsize=7.5)
    hz = res["hazard"]
    ps.caption(fig, "a: AID's visual record from 1932 to 1945 is almost entirely L. Peltier's, pinned near 9.8 (his "
                    "comparison star) until he saw the star fade in 1945 (Peltier 1945); the rise to the high state "
                    "in 1936 shows only in blue light on Harvard plates (Schaefer 2023), not in these estimates, so "
                    "levels in a and b cannot be compared. b: the 2023-24 fading recovered (Merc et al. 2025). The "
                    f"recurrence record alone (three intervals, mean {hz['mean_interval']:.2f} ± {hz['sd_interval']:.2f} yr; "
                    f"Pei et al. 2026) gives {100 * hz['p_rest_of_2026']:.0f}% for the rest of 2026 and "
                    f"{100 * hz['p_next_12_months']:.0f}% within twelve months: an empirical indicator, not physics.",
               y=0.0)
    return ps.save(fig, _out("13_tcrb_before"))


def fig14_tcrb_orbit():
    TC, df = _tcrb()
    res = _res()["tcrb"]
    p_orb, t_ref = L.SCHAEFER2023["p_orb_fekel"], L.SCHAEFER2023["t_vmax_fekel"]
    fig = plt.figure(figsize=(11.5, 4.8))
    gs = fig.add_gridspec(1, 2, wspace=0.22)
    for k, (key, lo, hi, title) in enumerate((("quiescence", 2447800.0, 2456700.0, "quiescence, 1990-2014"),
                                              ("high_state", 2457400.0, 2460000.0, "high state, 2016-2023"))):
        ax = fig.add_subplot(gs[k])
        ph, resid, day = TC.ellipsoidal_fold(df, "V", lo, hi)
        coef, model = TC.fourier_fit(ph, resid, 2)
        for s in (0.0, 1.0):
            ax.plot(ph + s, resid, ".", ms=2.5, color=ps.MUTED, alpha=0.6)
        fph = np.linspace(0, 2, 400)
        fm = coef[0] + sum(coef[2 * j - 1] * np.cos(2 * np.pi * j * fph) + coef[2 * j] * np.sin(2 * np.pi * j * fph)
                           for j in (1, 2))
        ax.plot(fph, fm, color=ps.BLUE, lw=1.8)
        for yr, col in ((1866, ps.ORANGE), (1946, ps.VIOLET)):
            phe = np.mod((L.SCHAEFER2023["peak_jd"][yr] - t_ref) / p_orb, 1.0)
            for s in (0.0, 1.0):
                ax.axvline(phe + s, color=col, lw=1.0)
            ax.text(phe + 0.02, -0.33 if yr == 1866 else -0.27, f"{yr}", fontsize=7.5, color=col)
        ax.set_ylim(0.38, -0.38)
        ax.set_xlim(0, 2)
        ax.set_xlabel("orbital phase (227.57-day orbit; 0 = red giant receding fastest)")
        ax.set_ylabel("V minus a 400-day running median (mag)")
        r = res[f"ellipsoidal_{key}"]
        _panel(ax, "ab"[k], f"{title}: {r['full_amp_2nd']:.2f} mag peak to peak")
    ps.caption(fig, "The red giant fills its Roche lobe, so the white dwarf's tide stretches it into a teardrop: brightest "
                    "twice per orbit, seen side-on, faintest at the conjunctions. At phase 0.25 the giant is behind the "
                    "white dwarf and shows us its inner face, heated by the white dwarf and disc: that minimum is "
                    f"{res['ellipsoidal_quiescence']['min_diff']:.2f} ± {res['ellipsoidal_quiescence']['min_diff_err']:.2f} mag "
                    "shallower than the one at phase 0.75 (Schaefer 2023 found no significant difference in visual "
                    "estimates). Fekel et al.'s (2000) ephemeris, P = 227.5687 d; daily V medians; blue: two-harmonic fit. "
                    "The lines mark the orbital phases of the two eruptions.", y=0.0)
    return ps.save(fig, _out("14_tcrb_orbit"))


# ---------------------------------------------------------------------------
# 15  Predictions against outcomes
# ---------------------------------------------------------------------------


def fig15_predictions():
    import v1309sco as V
    res = _res()
    v = res["v1309"]
    k = res["kic"]
    t = res["tcrb"]
    today = _year(t["last_jd"])
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    rows = [
        ("V1309 Sco: period law (fitted 2011, after the fact)", v["law_from_periods"]["t0_year"],
         v["law_from_periods"]["t0_year_err"], v["onset"]["t_on10_year"], "merged; the rise began 2008.2", ps.VIOLET),
        ("KIC 9832227: Molnar et al. 2015 (AAS)", L.MOLNAR2017["me15"]["t0_year"], L.MOLNAR2017["me15"]["t0_year_err"],
         None, "still orbiting in 2025", ps.VIOLET),
        ("KIC 9832227: Molnar et al. 2017", L.MOLNAR2017["t0_year"], L.MOLNAR2017["t0_year_err"], None,
         "still orbiting in 2025", ps.VIOLET),
        ("KIC 9832227: refit, with the mislabelled 1999 time", k["exponential_profiles"]["molnar_like_labelled"]["t0_best"],
         0.5 * (k["exponential_profiles"]["molnar_like_labelled"]["t0_hi"] -
                k["exponential_profiles"]["molnar_like_labelled"]["t0_lo"]), None, "", ps.RED),
        ("T CrB: Schaefer 2023 (template of 1946)", L.SCHAEFER2023["prediction"][0], L.SCHAEFER2023["prediction"][1],
         None, "window closes 2026.8: not erupted by Oct 9", ps.VIOLET),
        ("T CrB: 1945-dip analogy, if the 2023 fading was the dip", 2024.2, 0.15, None, "recovered instead", ps.VIOLET),
    ]
    for i, (lab, c, err, outcome, note, col) in enumerate(rows):
        y = len(rows) - i
        ax.plot([c - err, c + err], [y, y], color=col, lw=6, solid_capstyle="butt", alpha=0.35)
        ax.plot(c, y, "o", ms=6, color=col, mec=ps.SURFACE, mew=1.0)
        if outcome is not None:
            ax.plot(outcome, y, "*", ms=12, color=ps.ORANGE, mec=ps.SURFACE, mew=0.6)
        if note:
            ax.text(min(c + err, 2030) + 0.4, y + 0.28, note, fontsize=7.5, color=ps.INK_SECONDARY, va="bottom")
    ax.axvline(today, color=ps.INK, lw=1.0)
    ax.text(today + 0.15, 0.55, f"today\n{today:.2f}", fontsize=7.5, color=ps.INK)
    ax.set_xlim(2005.5, 2031)
    ax.set_ylim(0.3, len(rows) + 0.8)
    ax.set_yticks([len(rows) - i for i in range(len(rows))])
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.5, color=ps.INK)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xticks(np.arange(2006, 2031, 2))
    ax.set_xlabel("predicted date (dot) and stated uncertainty (bar); star: what happened")
    _panel(ax, "", "Predicting stellar catastrophes from light curves: the record so far")
    ax.grid(axis="y", visible=False)
    ps.caption(fig, "V1309 Sco's law reaches zero period at t0 by construction: the real merger began about two years "
                    "earlier, so a period law gives at best the order of the time left. The 1945-dip analogy: Schaefer "
                    "(2023) gave about one year from the start of the dip to the eruption; the 2023 fading began in "
                    "March 2023.", y=0.0)
    return ps.save(fig, _out("15_predictions"))


FIGURES = {"01": fig01_roche_lobes, "02": fig02_contact_light_curve, "03": fig03_period_distribution,
           "04": fig04_period_luminosity, "05": fig05_distance_map, "06": fig06_kic_oc, "07": fig07_kic_mistake,
           "08": fig08_v1309_lightcurve, "09": fig09_v1309_period, "10": fig10_v1309_shape,
           "11": fig11_tcrb_history, "12": fig12_tcrb_eruptions, "13": fig13_tcrb_before, "14": fig14_tcrb_orbit,
           "15": fig15_predictions}


def main(which=None):
    for key, fn in FIGURES.items():
        if which and key not in which:
            continue
        print(f"figure {key}: {fn.__name__}", flush=True)
        fn()


if __name__ == "__main__":
    main(sys.argv[1:] or None)
