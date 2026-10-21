"""Every figure of the project, drawn from cached data and results.

    python figures.py            # all figures that have their inputs cached
    python figures.py overview   # one figure by name

Colour roles are fixed across figures: data in ink, Newtonian predictions in
blue, MOND predictions in orange, a third series (another method or sample)
in aqua; context recedes in grey. Each figure is saved as PNG and SVG.
"""

from __future__ import annotations

import json
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

import plotstyle as PS
from config import DERIVED_DIR, FIG_DIR

PS.use_style()
C_NEWTON = PS.BLUE
C_MOND = PS.ORANGE
C_THIRD = PS.AQUA
C_DATA = PS.INK
RESULTS = DERIVED_DIR / "results"


def _res(name):
    p = RESULTS / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


# --------------------------------------------------------------------------
# 1. The catalogue
# --------------------------------------------------------------------------

def fig_overview():
    t = Table.read(DERIVED_DIR / "elbadry_all.fits")
    sh = Table.read(DERIVED_DIR / "elbadry_shift_all.fits")
    R = np.asarray(t["R_chance_align"])
    good = R < 0.1
    l = np.asarray(t["l1"])
    b = np.asarray(t["b1"])
    fig = plt.figure(figsize=(11.5, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.38, wspace=0.25)
    ax = fig.add_subplot(gs[0, :], projection="mollweide")
    lon = np.radians(np.where(l > 180, l - 360, l))[good]
    lat = np.radians(b)[good]
    H, xe, ye = np.histogram2d(-lon, lat, bins=[np.linspace(-np.pi, np.pi, 241), np.linspace(-np.pi / 2, np.pi / 2, 121)])
    area = np.outer(np.diff(xe), np.sin(ye[1:]) - np.sin(ye[:-1])) * (180 / np.pi) ** 2
    dens = H / area
    X, Y = np.meshgrid(xe, ye, indexing="ij")
    ld = np.log10(np.where(dens > 0, dens, np.nan))
    im = ax.pcolormesh(X, Y, ld, cmap=PS.CMAP_BLUE, rasterized=True,
                       vmin=np.nanpercentile(ld, 1), vmax=np.nanpercentile(ld, 99.7))
    ax.set_xticks(np.radians([-150, -90, -30, 30, 90, 150]))
    ax.set_xticklabels(["150°", "90°", "30°", "330°", "270°", "210°"], color=PS.MUTED, fontsize=7)
    ax.set_yticks(np.radians([-60, -30, 0, 30, 60]))
    ax.set_yticklabels(["−60°", "−30°", "0°", "30°", "60°"], color=PS.MUTED, fontsize=7)
    ax.grid(color=PS.GRID, linewidth=0.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02, aspect=25)
    cb.set_label("log10 pairs per deg²", color=PS.INK_SECONDARY)
    cb.outline.set_visible(False)
    ax.set_title(f"A million binaries: {good.sum():,} pairs with R < 0.1 (El-Badry, Rix & Heintz 2021), "
                 "Galactic coordinates", loc="center")

    ax = fig.add_subplot(gs[1, 0])
    d = 1000.0 / np.asarray(t["parallax1"])
    bins = np.linspace(0, 1100, 111)
    ax.hist(d[good], bins=bins, histtype="step", color=C_DATA, linewidth=1.4, label="all, R < 0.1")
    for dm, lab, c in ((130, "130 pc: Cookson et al. (2026)", C_MOND), (200, "200 pc: Chae (2023)", C_THIRD),
                       (250, "250 pc: Banik et al. (2024)", C_NEWTON)):
        ax.axvline(dm, color=c, linewidth=1.2, linestyle="--", label=lab)
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_xlabel("distance of the brighter star, 1/parallax (pc)")
    ax.set_ylabel("pairs per 10 pc")
    ax.set_xlim(0, 1100)
    ax.set_title("Distances: the tests use only the nearest few per cent", loc="left")

    ax = fig.add_subplot(gs[1, 1])
    s_all = np.log10(np.asarray(t["sep_AU"]))
    s_sh = np.log10(np.asarray(sh["sep_AU"]))
    sb = np.linspace(1.0, 5.4, 89)
    w = 1.0 / np.diff(sb)[0]
    for x, lab, c, lw in ((s_all, "all candidate pairs", C_DATA, 1.4),
                          (s_sh, "chance alignments (stars shifted 0.5°)", C_MOND, 1.4),
                          (s_all[good], "R < 0.1 (>90% bound)", C_NEWTON, 1.4)):
        h = np.histogram(x, bins=sb)[0] * w
        ax.step(sb[:-1], np.where(h > 0, h, np.nan), where="post", color=c, linewidth=lw, label=lab)
    ax.set_yscale("log")
    ax.axvspan(np.log10(2000), np.log10(30000), color=PS.GRID, alpha=0.5, linewidth=0, zorder=0)
    ax.text(np.log10(7700), 3e1, "range of the\ngravity tests", ha="center", color=PS.INK_SECONDARY, fontsize=7.5)
    ax.set_xlabel("log10 projected separation s (AU)")
    ax.set_ylabel("dN / dlog s")
    ax.set_ylim(1e1, 3e7)
    ax.legend(loc="upper left", fontsize=7.5, ncol=1)
    ax.set_title("Chance alignments take over beyond ~30 kAU", loc="left")
    PS.caption(fig, "Data: Gaia EDR3 pairs within 1 kpc from El-Badry et al. (2021), Zenodo 4435257; R is their "
                    "chance-alignment ratio. The shifted catalogue pairs each star with stars displaced 0.5° on "
                    "the sky, so all its pairs are chance alignments with the real catalogue's statistics.")
    return PS.save(fig, FIG_DIR / "fig01_catalogue")


# --------------------------------------------------------------------------
# 2. What MOND predicts for a binary in the Galaxy's field
# --------------------------------------------------------------------------

def fig_mond_force():
    import qumond as Q
    gf = Q.galactic_field()
    f = Q.cached_field(gf["gNe_a0"])
    r = np.exp(f.lnr)
    fig = plt.figure(figsize=(12.0, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.05], wspace=0.42)
    ax = fig.add_subplot(gs[0])
    ax.axhline(1.0, color=C_NEWTON, linewidth=1.8, label="Newton")
    ax.plot(r, f.eta, color=C_MOND, linewidth=2.0, label="MOND (QUMOND) in the Galactic field")
    # isolated MOND for contrast (no external field): eta = nu(1/r^2) for a point mass
    ax.plot(r, Q.nu_simple(1.0 / r**2), color=C_MOND, linewidth=1.2, linestyle=":",
            label="MOND without the Galaxy (isolated star)")
    ax.axhline(gf["eta_inf"], color=PS.MUTED, linewidth=0.8, linestyle="--")
    ax.text(0.012, gf["eta_inf"] + 0.03, f"far-field limit ν_e(1 + K_e/3) = {gf['eta_inf']:.3f}", color=PS.INK_SECONDARY,
            fontsize=7.5)
    ax.axvline(1.0, color=PS.MUTED, linewidth=0.8)
    ax.text(1.25, 1.95, "MOND radius\nr_M = √(GM/a₀)\n= 7.0 kAU × √(M/M☉)", color=PS.INK_SECONDARY, fontsize=7.5, va="top")
    ax.set_xscale("log")
    ax.set_xlim(1e-2, 1e2)
    ax.set_ylim(0.9, 2.6)
    ax.set_xlabel("separation r / r_M")
    ax.set_ylabel("η = mean inward gravity / Newtonian")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("Gravity boost felt by a wide binary", loc="left")

    ax = fig.add_subplot(gs[1])
    # field map in the (x, z) plane: the Galaxy pulls towards +z
    th = np.arccos(np.clip(f.mu, -1, 1))
    xs = np.linspace(-3, 3, 301)
    X, Z = np.meshgrid(xs, xs)
    RR = np.hypot(X, Z)
    TT = np.arccos(np.clip(Z / np.maximum(RR, 1e-9), -1, 1))
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator((f.lnr, th), f.ratio_r, bounds_error=False, fill_value=None)
    val = interp(np.stack([np.log(np.maximum(RR, 1e-3)), TT], -1))
    im = ax.pcolormesh(X, Z, val, cmap=PS.CMAP_BLUE, vmin=1.0, vmax=1.6, rasterized=True, shading="auto")
    cs = ax.contour(X, Z, val, levels=[1.1, 1.3, 1.5], colors=[PS.INK_SECONDARY], linewidths=0.6)
    ax.clabel(cs, fmt="%.1f", fontsize=6.5)
    ax.plot(0, 0, marker="*", color=PS.INK, markersize=9)
    ax.annotate("", xy=(0, 2.85), xytext=(0, 2.0), arrowprops=dict(arrowstyle="-|>", color=PS.INK, lw=1.2))
    ax.text(-0.15, 2.42, "to the Galactic\nCentre", fontsize=7, color=PS.INK, ha="right", va="center")
    ax.set_aspect("equal")
    ax.set_xlabel("x / r_M")
    ax.set_ylabel("z / r_M")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.9, pad=0.16, aspect=30)
    cb.set_label("radial gravity / Newtonian", color=PS.INK_SECONDARY)
    cb.outline.set_visible(False)
    ax.set_title("…and it is lopsided", loc="left")

    ax = fig.add_subplot(gs[2])
    # radial acceleration relation: galaxies (isolated, D4) vs a binary in the Galactic field
    gN = np.logspace(-12.5, -7.5, 200)
    a0 = 1.2e-10
    ax.plot(gN, gN, color=C_NEWTON, linewidth=1.8, label="Newton: g = g_N")
    ax.plot(gN, Q.nu_simple(gN / a0) * gN, color=C_MOND, linewidth=1.2, linestyle=":",
            label="MOND, isolated (galaxy outskirts)")
    # binary: g_N = GM/r^2 -> r/r_M = sqrt(a0/g_N)
    rr = np.sqrt(a0 / gN)
    ax.plot(gN, np.interp(np.log(rr), f.lnr, f.eta) * gN, color=C_MOND, linewidth=2.0,
            label="MOND, binary in the Galactic field")
    ax.axvline(gf["ge_a0"] * a0, color=PS.MUTED, linewidth=0.8, linestyle="--")
    ax.text(gf["ge_a0"] * a0 * 1.15, 1.5e-12, f"Galaxy's pull\n{gf['ge_a0']:.2f} a₀", fontsize=7, color=PS.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Newtonian acceleration g_N (m s⁻²)")
    ax.set_ylabel("predicted acceleration g (m s⁻²)")
    ax.legend(loc="upper left", fontsize=7)
    ax.set_title("Why only ~20% faster", loc="left")
    PS.caption(fig, f"Our QUMOND solution (simple interpolating function, a₀ = 1.2×10⁻¹⁰ m s⁻²) for a star in the "
                    f"Galaxy's field g_e = {gf['ge_a0']:.3f} a₀ (v_c = 232.8 km/s at R₀ = 8.2 kpc; Newtonian-equivalent "
                    f"{gf['gNe_a0']:.3f} a₀). η → {gf['eta_inf']:.4f}, so orbital speeds rise by at most "
                    f"√η = {np.sqrt(gf['eta_inf']):.3f}. Middle: the ratio at each point of the x–z plane.")
    return PS.save(fig, FIG_DIR / "fig02_mond_force")


FIGURES = {"overview": fig_overview, "mond_force": fig_mond_force}

# --------------------------------------------------------------------------
# 3. The forward-model test on our rebuild of Banik et al.'s sample
# --------------------------------------------------------------------------
BANIK_ALPHA = (-0.021, 0.065, 0.045)      # Banik et al. (2024) table 4: mode, +, -


def _binom_ll(k, N, p):
    from scipy.special import gammaln
    p = np.clip(p, 1e-300, 1 - 1e-15)
    return gammaln(N + 1) - gammaln(k + 1) - gammaln(N - k + 1) + k * np.log(p) + (N - k) * np.log1p(-p)


def fig_forward(name: str = "forward_banik2024", out: str = "fig04_forward_test",
                title: str = "our rebuild of Banik et al. (2024)'s sample"):
    r = _res(name)
    if r is None:
        return None
    px = r["pixels"]
    k = np.array(px["counts"])
    pN, pM = np.array(px["p_newton"]), np.array(px["p_mond"])
    re, ve = np.array(px["r_edges"]), np.array(px["v_edges"])
    N = k.sum()
    groups = [(2, 3), (3, 5), (5, 12), (12, 30)]
    fig = plt.figure(figsize=(13.5, 7.8))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.95], hspace=0.45, wspace=0.32)
    vc = 0.5 * (ve[1:] + ve[:-1])
    dv = np.diff(ve)
    for j, (lo, hi) in enumerate(groups):
        ax = fig.add_subplot(gs[0, j])
        rows = (re[:-1] >= lo - 1e-9) & (re[1:] <= hi + 1e-9)
        kk = k[rows].sum(0)
        n_g = kk.sum()
        dens = kk / dv / n_g
        err = np.sqrt(np.maximum(kk, 1)) / dv / n_g
        for p, c, lab in ((pN, C_NEWTON, "best Newtonian model"), (pM, C_MOND, "best MOND model")):
            m = p[rows].sum(0) * N
            ax.step(ve, np.r_[m / dv / n_g, (m / dv / n_g)[-1]], where="post", color=c, linewidth=1.6, label=lab)
        ax.errorbar(vc, dens, yerr=err, fmt="o", color=C_DATA, markersize=3.2, elinewidth=0.8, capsize=0,
                    label="observed", zorder=5)
        ax.set_xlim(0, 3.0)
        ax.set_xlabel("v~ = v_sky / Newtonian v_c")
        if j == 0:
            ax.set_ylabel("fraction of pairs per unit v~")
            handles, labels = ax.get_legend_handles_labels()
        ax.set_title(f"r_sky {lo}–{hi} kAU ({int(n_g):,} pairs)", loc="left", fontsize=9.5)
    # per-separation likelihood difference
    gs2 = gs[1, :].subgridspec(1, 3, wspace=0.34, width_ratios=[1.15, 1.0, 1.1])
    ax = fig.add_subplot(gs2[0])
    dl = (_binom_ll(k, N, pN) - _binom_ll(k, N, pM)).sum(1)
    rc = np.sqrt(re[1:] * re[:-1])
    width = np.diff(np.log10(re)) * 0.85
    cols = np.where(dl >= 0, C_NEWTON, C_MOND)
    ax.bar(np.log10(rc), dl, width=width, color=cols, edgecolor=PS.SURFACE, linewidth=1.0)
    ax.axhline(0, color=PS.AXIS, linewidth=0.8)
    ax.set_xticks(np.log10([2, 3, 5, 10, 20, 30]))
    ax.set_xticklabels(["2", "3", "5", "10", "20", "30"])
    ax.set_xlabel("projected separation r_sky (kAU)")
    ax.set_ylabel("Δ ln L  (Newton − MOND)")
    tot = r["delta_lnP_newton_minus_mond"]
    ax.set_title(f"Eccentricity law free: Newton by Δ ln L = {tot:.0f} (≈ {r['sigma_newton_over_mond']:.0f}σ)",
                 loc="left", fontsize=9)
    ax.text(0.99, 0.95, "blue: prefers Newton\norange: prefers MOND", transform=ax.transAxes, ha="right",
            va="top", fontsize=7.2, color=PS.INK_SECONDARY)
    # Newton's advantage against a fixed eccentricity law, without and with the v-r angles
    ax = fig.add_subplot(gs2[1])
    gp = _res("gamma_profile")
    va = _res("vr_angle_test")
    if gp:
        g = np.array(gp["gamma"])
        lnN, lnM = np.array(gp["newton"]), np.array(gp["mond"])
        ax.axvspan(1.32 - 0.09, 1.32 + 0.09, color=PS.GRID, zorder=0)
        ax.text(1.32, 0.55, "Hwang et al.", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=7, color=PS.INK_SECONDARY, rotation=90)
        ax.plot(g, lnN - lnN.max(), "o-", color=C_NEWTON, markersize=4, linewidth=1.6, label="Newton, (r_sky, v~)")
        ax.plot(g, lnM - lnN.max(), "s-", color=C_MOND, markersize=4, linewidth=1.6, label="MOND, (r_sky, v~)")
        if va:
            cN = lnN + np.array(va["newton"])
            cM = lnM + np.array(va["mond"])
            ax.plot(g, cN - cN.max(), "o:", color=C_NEWTON, markersize=4, markerfacecolor=PS.SURFACE, linewidth=1.4,
                    label="Newton, + v–r angles")
            ax.plot(g, cM - cN.max(), "s:", color=C_MOND, markersize=4, markerfacecolor=PS.SURFACE, linewidth=1.4,
                    label="MOND, + v–r angles")
        ax.set_xlabel("eccentricity index γ, held fixed  (P(e) ∝ e^γ)")
        ax.set_ylabel("ln L − best Newtonian ln L")
        ax.set_xlim(0.85, 4.15)
        ax.legend(loc="lower right", fontsize=6.8)
        ax.set_title("Fixed eccentricity law", loc="left", fontsize=9)
    # alpha posterior
    ax = fig.add_subplot(gs2[2])
    ch = np.array(r["mcmc"]["chain_thin"])
    ia = r["mcmc"]["names"].index("alpha")
    ax.hist(ch[:, ia], bins=60, density=True, histtype="stepfilled", color=PS.GRID, edgecolor=C_DATA,
            linewidth=1.2, label="ours, eccentricity law free")
    a = r["posterior"]["alpha"]
    ax.axvline(0, color=C_NEWTON, linewidth=1.6, label="Newton (α = 0)")
    ax.axvline(1, color=C_MOND, linewidth=1.6, label="MOND (α = 1)")
    yt = ax.get_ylim()[1]
    vh = _res("variant_ecc_hwang")
    if vh:
        q = vh["posterior"]["alpha"]
        ax.errorbar([q["median"]], [0.72 * yt], xerr=[[q["median"] - q["lo"]], [q["hi"] - q["median"]]], fmt="o",
                    color=C_DATA, markersize=4.5, capsize=3, label="ours, with Hwang et al.'s γ prior")
    ax.errorbar([BANIK_ALPHA[0]], [0.5 * yt], xerr=[[BANIK_ALPHA[2]], [BANIK_ALPHA[1]]], fmt="D", color=C_THIRD,
                markersize=4.5, markerfacecolor=PS.SURFACE, capsize=3, label="Banik et al. 2024, γ free")
    ax.errorbar([0.024], [0.36 * yt], xerr=[[0.054], [0.042]], fmt="D", color=C_THIRD, markersize=4.5, capsize=3,
                label="Banik et al. 2024, γ prior")
    ax.set_xlim(-0.35, 1.15)
    ax.set_xlabel("gravity dial α (0 = Newton, 1 = MOND)")
    ax.set_ylabel("posterior density")
    ax.legend(loc="upper right", fontsize=6.8, frameon=True, facecolor=PS.SURFACE, edgecolor=PS.SURFACE,
              framealpha=1.0)
    ax.set_title(f"Gravity free: α = {a['median']:.3f} (+{a['hi'] - a['median']:.3f} −{a['median'] - a['lo']:.3f})",
                 loc="left", fontsize=9)
    fig.suptitle(f"The wide binary test, forward-modelled: {title}", x=0.01, ha="left", y=0.995)
    fig.legend(handles, labels, loc="upper right", ncol=3, fontsize=8, bbox_to_anchor=(0.99, 0.985))
    PS.caption(fig, f"{int(N):,} pairs, 2–30 kAU, v~ < 5, in Banik et al.'s 540 pixels. Models: our orbit libraries "
                    "(Newton; QUMOND with the Galactic external field) plus undetected close companions and "
                    "line-of-sight pairs, fitted by maximum likelihood with the gravity law fixed (top, bottom left), "
                    "with the eccentricity index also fixed (bottom middle; dotted: adding the likelihood of the "
                    "observed v–r angles, treated as independent), and by MCMC with gravity free (bottom right). "
                    "Error bars: Poisson (top), 68% intervals (bottom right).")
    return PS.save(fig, FIG_DIR / out)


FIGURES["forward"] = fig_forward


# --------------------------------------------------------------------------
# 5. Chae's acceleration-plane test, and what it reports for mock universes
# --------------------------------------------------------------------------
INJ_ROWS = [
    ("data|original", "Gaia data, Chae's method as published", "data"),
    ("data|corrected", "Gaia data, Newtonian expectation given the data's noise + perspective", "data"),
    ("newton_full|original", "Newtonian mock: companions + noise + perspective", "newton"),
    ("newton_noise_only|original", "Newtonian mock: companions + noise", "newton"),
    ("newton_persp_only|original", "Newtonian mock: companions + perspective", "newton"),
    ("newton_cb_only|original", "Newtonian mock: companions only", "newton"),
    ("newton_no_cb|original", "Newtonian mock: noise + perspective, no companions", "newton"),
    ("newton_full|corrected", "Newtonian mock, expectation given its noise + perspective", "newton"),
    ("mond_full|original", "MOND mock: companions + noise + perspective", "mond"),
    ("mond_cb_only|original", "MOND mock: companions only", "mond"),
    ("mond_full|corrected", "MOND mock, expectation given its noise + perspective", "mond"),
]


def fig_chae():
    r = _res("injection_chae")
    if r is None:
        return None
    import chae as CH
    d = CH.chae_sample(CH.load_chae_base())
    rng = np.random.default_rng(3)
    f = r["runs"]["data|original"]["f_multi"]
    m = CH.companion_masses(d["MA"], d["MB"], f, rng)
    mt = m["A1"] + m["A2"] + m["B1"] + m["B2"]
    gN, g = CH.deproject(d["s"], CH.observed_vp(d), mt, CH._draw_e(d["e"], d["e0"], d["e1"], rng), rng)
    fig = plt.figure(figsize=(13.5, 9.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.35], wspace=0.25)
    ax = fig.add_subplot(gs[0])
    H, xe, ye = np.histogram2d(gN, g, bins=[np.linspace(-12.5, -6.5, 121), np.linspace(-12.5, -6.5, 121)])
    ax.pcolormesh(xe, ye, np.log10(np.where(H.T > 0, H.T, np.nan)), cmap=PS.CMAP_CONTEXT, rasterized=True)
    ax.plot([-12.5, -6.5], [-12.5, -6.5], color=PS.MUTED, linewidth=0.8, linestyle="--")
    ax.text(-7.2, -6.9, "g = g_N", color=PS.INK_SECONDARY, fontsize=7.5, rotation=45, ha="center")
    # bin medians: data and Newtonian expectation (from the 200-set ensemble)
    run = r["runs"]["data|original"]
    for k, (lo, hi) in enumerate(zip(CH.X0_BINS[:-1], CH.X0_BINS[1:])):
        x0 = 0.5 * (lo + hi)
        for med, c, lab in ((run["obs_medians"][k], C_DATA, "data median"),
                            (run["newton_medians"][k], C_NEWTON, "Newtonian expectation")):
            # a point at orthogonal offset med from the line, at x0
            xx, yy = x0 - med / np.sqrt(2), x0 + med / np.sqrt(2)
            ax.plot(xx, yy, "o", color=c, markersize=7, markeredgecolor=PS.SURFACE, markeredgewidth=1.5,
                    label=lab if k == 0 else None, zorder=5)
    for b in CH.X0_BINS:
        ax.plot([b + 1.5, b - 1.5], [b - 1.5, b + 1.5], color=PS.AXIS, linewidth=0.6)
    ax.set_xlim(-12.5, -6.5)
    ax.set_ylim(-12.5, -6.5)
    ax.set_aspect("equal")
    ax.set_xlabel("log10 Newtonian acceleration g_N = GM/r² (m s⁻²)")
    ax.set_ylabel("log10 kinematic acceleration g = v²/r (m s⁻²)")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("Chae's acceleration plane (one Monte Carlo set)", loc="left")

    ax = fig.add_subplot(gs[1])
    runs = {f"inj:{k}": v for k, v in r["runs"].items()}
    for name, pre in (("injection_controls", "ctl"), ("injection_mechanism", "mech"),
                      ("injection_mechanism_opik", "mech2"), ("injection_chae_v2", "v2")):
        extra = _res(name)
        if extra:
            extra = extra.get("runs", extra)
            runs.update({f"{pre}:{k}": v for k, v in extra.items() if isinstance(v, dict) and "gamma_g" in v})
    groups = [
        ("The real sky (Chae's 26,615 pairs)", [
            ("inj:data|original", "as published", "data"),
            ("inj:data|corrected", "noise + perspective added to the Newtonian expectation", "data")]),
        ("Newtonian mock universes (the same pairs)", [
            ("inj:newton_full|original", "companions + noise + perspective", "newton"),
            ("inj:newton_noise_only|original", "companions + noise", "newton"),
            ("inj:newton_persp_only|original", "companions + perspective", "newton"),
            ("inj:newton_cb_only|original", "companions only", "newton"),
            ("inj:newton_no_cb|original", "noise + perspective, no companions", "newton"),
            ("inj:newton_full|corrected", "all three, noise + perspective in the expectation", "newton")]),
        ("MOND mock universes", [
            ("inj:mond_full|original", "companions + noise + perspective", "mond"),
            ("inj:mond_cb_only|original", "companions only", "mond"),
            ("inj:mond_full|corrected", "all three, noise + perspective in the expectation", "mond")]),
        ("Mocks from a forward fit to Chae's own sample (0.5–30 kAU)", [
            ("v2:newton|original", "Newton", "newton"),
            ("v2:newton|corrected", "Newton, noise + perspective in the expectation", "newton"),
            ("v2:mond_same_population|original", "MOND, same population", "mond"),
            ("v2:mond_refit|original", "MOND, its own best-fitting population", "mond")]),
        ("Pure Newtonian binaries: no companions, noise or perspective", [
            ("ctl:c1_pure_hwang_e", "eccentricities estimated per pair, as Hwang et al. (Chae's default)", "newton"),
            ("mech:true_eccentricity", "  ... each pair's true eccentricity instead", "newton"),
            ("mech:statistical_eccentricity", "  ... eccentricities from the population law", "newton"),
            ("mech2:opik_true_eccentricity", "  ... true eccentricity, scale-free separation law", "newton"),
            ("mech2:opik_hwang_individual", "  ... per-pair estimates, scale-free separation law", "newton"),
            ("ctl:c2_pure_superthermal_e", "more eccentric orbits (γ = 2.1), per-pair estimates", "newton"),
            ("ctl:c3_companions_hwang_e", "with hidden companions, per-pair estimates", "newton")]),
    ]
    cols = {"data": C_DATA, "newton": C_NEWTON, "mond": C_MOND}
    lines = []
    for head, rows in groups:
        rows = [row for row in rows if row[0] in runs]
        if rows:
            lines.append(("head", head))
            lines += [("row", row) for row in rows]
    n = len(lines)
    for j, (kind_, item) in enumerate(lines):
        y = n - 1 - j
        if kind_ == "head":
            ax.text(2.42, y - 0.05, item, ha="left", va="center", fontsize=7.6, color=PS.INK, fontweight="bold")
            continue
        k, lab, kind = item
        run = runs[k]
        g0 = run["gamma_g"][0]
        lo = 10 ** (np.sqrt(2) * (run["delta"][0] - run["delta_sd"][0]))
        hi = 10 ** (np.sqrt(2) * (run["delta"][0] + run["delta_sd"][0]))
        ax.errorbar([g0], [y], xerr=[[g0 - lo], [hi - g0]], fmt="o", color=cols[kind], markersize=5.5,
                    capsize=2, elinewidth=1.2, zorder=4)
        ax.text(2.48, y, lab, ha="left", va="center", fontsize=7.2, color=PS.INK)
        ax.text(hi + 0.045, y, f"{g0:.2f}", ha="left", va="center", fontsize=7.0, color=PS.INK_SECONDARY)
    ax.axvline(1.0, color=C_NEWTON, linewidth=1.0)
    ax.axvline(Q_ETA, color=C_MOND, linewidth=1.0, linestyle="--")
    ax.text(1.0, n - 0.2, "Newton", color=PS.INK_SECONDARY, fontsize=7.5, ha="center")
    ax.text(Q_ETA, n - 0.2, "MOND limit", color=PS.INK_SECONDARY, fontsize=7.5, ha="center")
    ax.set_xlim(0.8, 4.3)
    ax.set_ylim(-0.7, n + 0.3)
    ax.set_yticks([])
    ax.set_xticks([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(0.8, 2.3)
    ax.grid(False)
    for gx in (0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2):
        ax.axvline(gx, color=PS.GRID, linewidth=0.6, zorder=0)
    ax.set_xlabel("reported gravity boost g_obs / g_pred, lowest-acceleration bin (x₀ ≈ −10.3)", x=0.22)
    ax.set_title("What the method reports for real data and for mock universes", loc="left")
    PS.caption(fig, f"Chae's (2023, corrected 2024) sample of {r['n_pairs']:,} pairs, reproduced exactly from his "
                    "public files, and his method reimplemented from his code. Mocks: the same pairs with "
                    "relative motions drawn from our forward model (Newtonian or MOND; hidden companions as fitted "
                    "to our Banik-style sample, or in the fourth group the whole population as fitted to Chae's own "
                    "sample), Gaia-like noise, each pair's perspective effect and the catalogue's selection; every "
                    "mock is analysed exactly like the data. Bars: spread over 200 Monte Carlo sets.")
    return PS.save(fig, FIG_DIR / "fig06_chae_method")


def _q_eta():
    try:
        import qumond as Q
        return Q.galactic_field()["eta_inf"]
    except Exception:  # noqa: BLE001
        return 1.42


Q_ETA = _q_eta()
FIGURES["chae"] = fig_chae


# --------------------------------------------------------------------------
# 6. Which choices move the answer?
# --------------------------------------------------------------------------
CHAE_2026 = (1.60, 0.17, 0.14)      # Chae et al. (2026) 3D result, as quoted by Saad & Ting (2026)


def _row(ax, y, med, lo, hi, color, label, val_fmt="{:.2f}", marker="o"):
    ax.errorbar([med], [y], xerr=[[med - lo], [hi - med]], fmt=marker, color=color, markersize=5.5, capsize=2,
                elinewidth=1.2)
    ax.text(-0.02, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.4,
            color=PS.INK)


BANIK_TABLE4 = {      # Banik et al. (2024) table 4: alpha_grav mode, +, -
    "nominal": (-0.021, 0.065, 0.045), "ecc_hwang": (0.024, 0.042, 0.054), "fcb_03": (0.149, 0.035, 0.048),
    "mass_1_2": (-0.038, 0.076, 0.052), "mstar_gt_05": (-0.177, 0.064, 0.081)}


def _gg_interval(run):
    """Chae's statistic and its Monte Carlo spread as (median, lo, hi)."""
    d0, sd = run["delta"][0], run["delta_sd"][0]
    return run["gamma_g"][0], 10 ** (np.sqrt(2) * (d0 - sd)), 10 ** (np.sqrt(2) * (d0 + sd))


def fig_choices():
    fb = _res("forward_banik2024")
    if fb is None:
        return None
    # ---- A: forward model, ours (ink) against Banik et al.'s (aqua) for the same choice ----
    rows_a = [("nominal: eccentricity law free", fb["posterior"]["alpha"], BANIK_TABLE4["nominal"])]
    labels_a = [("ecc_hwang", "eccentricity prior γ = 1.32 ± 0.09 (Hwang et al.)"),
                ("fcb_03", "companion fraction fixed at 0.30"), ("mass_1_2", "total mass 1–2 M☉"),
                ("mstar_gt_05", "both stars > 0.5 M☉"), ("noise", "measurement noise modelled"),
                ("absolute_cb", "companion orbits independent of the wide orbit"),
                ("vt_lt_25", "only v~ < 2.5 fitted"), ("d_lt_150", "distance < 150 pc"),
                ("ruwe_lt_125", "RUWE < 1.25 for both stars"), ("both_rv", "both stars have Gaia RVs"),
                ("no_perspective", "no perspective correction")]
    for name, lab in labels_a:
        v = _res(f"variant_{name}")
        if v and "posterior" in v:
            rows_a.append((lab, v["posterior"]["alpha"], BANIK_TABLE4.get(name)))
        if name == "ecc_hwang":
            ag = _res("alpha_at_fixed_gamma")
            if ag and 2.0 in ag["gamma"]:
                a_ = ag["alpha"][ag["gamma"].index(2.0)]
                rows_a.append(("γ held at 2, where the v–r angles peak (best fit only)", dict(median=a_, lo=a_, hi=a_),
                               None))
    for tag, lab in (("forward_chae2023_2to30_relative", "Chae's sample, 2–30 kAU, its noise + selection modelled"),
                     ("forward_chae2023_2to30_absolute", "  ... companion orbits independent of the wide orbit"),
                     ("forward_chae2023_05to30_absolute", "  ... and fitted down to 0.5 kAU (best fit only)")):
        v = _res(tag)
        if v and "posterior" in v:
            rows_a.append((lab, v["posterior"]["alpha"], None))
        elif v and "free" in v.get("best", {}):
            a_ = v["best"]["free"]["theta"]["alpha"]
            rows_a.append((lab, dict(median=a_, lo=a_, hi=a_), None))
    # ---- B: Chae's statistic on the sky, a Newtonian mock and a MOND mock, one choice per row ----
    cv = _res("chae_variants") or {}
    inj = (_res("injection_chae") or {}).get("runs", {})
    sub = _res("injection_subsets") or {}
    mech2 = _res("injection_mechanism_opik") or {}
    v2 = (_res("injection_chae_v2") or {}).get("runs", {})
    rows_b = [("as published", inj.get("data|original", cv.get("published")), inj.get("newton_full|original"),
               inj.get("mond_full|original")),
              ("noise + perspective in the Newtonian expectation",
               inj.get("data|corrected", cv.get("noise_and_perspective_in_expectation")),
               inj.get("newton_full|corrected"), inj.get("mond_full|corrected")),
              ("eccentricities from the population law, not per pair", cv.get("statistical_eccentricities"),
               mech2.get("newton_full_statistical"), None),
              ("within 80 pc (smaller proper-motion errors)", cv.get("within_80pc"), sub.get("within_80pc|newton"),
               sub.get("within_80pc|mond")),
              ("precise half: σ(v~) < 0.05", cv.get("precise_only"), sub.get("precise_only|newton"),
               sub.get("precise_only|mond"))]
    if v2:
        rows_b.append(("mocks fitted to Chae's own sample (0.5–30 kAU)", inj.get("data|original", cv.get("published")),
                       v2.get("newton|original"), v2.get("mond_same_population|original")))
    rows_b = [rw for rw in rows_b if rw[1] is not None]
    # ---- C: 3D ----
    rows_c = []
    td = _res("three_d")
    if td:
        names = {"opik_drop0": "all 36 pairs, log-flat P(a)", "bpl_drop0": "all 36, broken power-law P(a)",
                 "opik_drop2": "  without the 2 fastest", "opik_drop4": "  without the 4 fastest",
                 "opik_drop6": "  without the 6 pairs too fast to be bound", "bpl_drop6": "  ... broken power-law P(a)"}
        for k, lab in names.items():
            if k in td["runs"]:
                q = td["runs"][k]
                rows_c.append((lab, dict(median=q["median"], lo=q["q16"], hi=q["q84"]), PS.INK))
        st = td.get("saad_ting", {})
        for k, lab in (("baseline", "Saad & Ting 2026, free semi-major axis"),
                       ("geometric_deprojection", "Saad & Ting 2026, geometric deprojection")):
            if k in st:
                rows_c.append((lab, dict(median=st[k]["median"], lo=st[k]["q16"], hi=st[k]["q84"]), C_THIRD))
        rows_c.append(("Chae et al. 2026 (published)", dict(median=CHAE_2026[0], lo=CHAE_2026[0] - CHAE_2026[2],
                                                    hi=CHAE_2026[0] + CHAE_2026[1]), C_THIRD))
    heights = [len(rows_a) + 1.4, 1.35 * len(rows_b) + 1.4] + ([len(rows_c) + 1.4] if rows_c else [])
    H = 0.31 * sum(heights) + 2.2
    fig = plt.figure(figsize=(11.5, H))
    gs = fig.add_gridspec(len(heights), 1, height_ratios=heights, hspace=0.5, left=0.42, right=0.97,
                          top=1 - 0.95 / H, bottom=0.75 / H)

    def frame(ax, n, xl, x_newton, x_mond, title):
        ax.axvline(x_newton, color=C_NEWTON, linewidth=1.2, zorder=1)
        ax.axvline(x_mond, color=C_MOND, linewidth=1.2, linestyle="--", zorder=1)
        ax.set_xlim(*xl)
        ax.set_ylim(-0.8, n - 0.2)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_title(title, loc="left", fontsize=9.5, x=-0.72)

    def label(ax, y, text):
        ax.text(-0.01, y, text, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.6,
                color=PS.INK)

    # A
    ax = fig.add_subplot(gs[0])
    n = len(rows_a)
    for j, (lab, q, bk) in enumerate(rows_a):
        y = n - 1 - j
        c = C_THIRD if lab.startswith(("Chae", "  ...")) else PS.INK
        med = q["median"]
        ax.errorbar([med], [y + (0.14 if bk else 0)], xerr=[[med - q["lo"]], [q["hi"] - med]], fmt="o", color=c,
                    markersize=5.0, capsize=2, elinewidth=1.2, zorder=4)
        ax.text(q["hi"] + 0.03, y + (0.14 if bk else 0), f"{med:.2f}", ha="left", va="center", fontsize=6.8,
                color=PS.INK_SECONDARY)
        if bk:
            ax.errorbar([bk[0]], [y - 0.2], xerr=[[bk[2]], [bk[1]]], fmt="D", color=C_THIRD, markersize=4.0,
                        markerfacecolor=PS.SURFACE, capsize=2, elinewidth=1.0, zorder=4)
        label(ax, y, lab)
    ax.errorbar([], [], xerr=[], fmt="o", color=PS.INK, markersize=5, label="this work (median, 68%)")
    ax.errorbar([], [], xerr=[], fmt="D", color=C_THIRD, markersize=4, markerfacecolor=PS.SURFACE,
                label="Banik et al. 2024, same choice (mode, 68%)")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=7.2, frameon=False)
    frame(ax, n, (-0.35, 1.15), 0.0, 1.0, "A. Forward model (Banik et al.'s approach): gravity dial α, 0 = Newton, "
                                         "1 = MOND")
    # B
    ax = fig.add_subplot(gs[1])
    n = len(rows_b)
    for j, (lab, dat, mN, mM) in enumerate(rows_b):
        y = n - 1 - j
        for run, c, mk, dy in ((dat, C_DATA, "o", 0.3), (mN, C_NEWTON, "s", 0.0), (mM, C_MOND, "^", -0.3)):
            if run is None:
                continue
            g0, lo, hi = _gg_interval(run)
            ax.errorbar([g0], [y + dy], xerr=[[g0 - lo], [hi - g0]], fmt=mk, color=c, markersize=4.6, capsize=2,
                        elinewidth=1.1, zorder=4)
            ax.text(hi + 0.03, y + dy, f"{g0:.2f}", ha="left", va="center", fontsize=6.6, color=PS.INK_SECONDARY)
        label(ax, y, lab)
    for c, mk, lab in ((C_DATA, "o", "the real sky"), (C_NEWTON, "s", "Newtonian mock"), (C_MOND, "^", "MOND mock")):
        ax.errorbar([], [], xerr=[], fmt=mk, color=c, markersize=4.6, label=lab)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=7.2, frameon=False)
    frame(ax, n, (0.85, 2.45), 1.0, Q_ETA, "B. Deprojection (Chae's approach): g_obs / g_pred, lowest-acceleration bin")
    # C
    if rows_c:
        ax = fig.add_subplot(gs[2])
        n = len(rows_c)
        for j, (lab, q, c) in enumerate(rows_c):
            y = n - 1 - j
            med = q["median"]
            ax.errorbar([med], [y], xerr=[[max(med - q["lo"], 0)], [max(q["hi"] - med, 0)]], fmt="o", color=c,
                        markersize=5.0, capsize=2, elinewidth=1.2, zorder=4)
            ax.text(min(q["hi"], 2.5) + 0.02, y, f"{med:.2f}", ha="left", va="center", fontsize=6.8,
                    color=PS.INK_SECONDARY)
            label(ax, y, lab)
        frame(ax, n, (0.6, 2.6), 1.0, Q_ETA, "C. 3D velocities, 36 pairs: gravity boost γ = G_eff / G")
    fig.suptitle("Which choices move the answer?", x=0.01, ha="left", y=1 - 0.18 / H)
    PS.caption(fig, "Blue lines: Newton; orange dashed: MOND (α = 1, or its asymptotic boost 1.42). A: one choice "
                    "changed per row; aqua circles: fits to Chae's sample. B: Chae's statistic with one choice "
                    "changed per row, for the real sky and for Newtonian and MOND mock skies analysed identically "
                    "(bars: spread of 200 Monte Carlo sets). C: the 3D test; aqua: published values.")
    return PS.save(fig, FIG_DIR / "fig08_choices")


FIGURES["choices"] = fig_choices


# --------------------------------------------------------------------------
# 2b. Kepler's envelope
# --------------------------------------------------------------------------

def fig_envelope():
    import chae as CH
    import kinematics as KIN
    import qumond as Q
    d = CH.chae_sample(CH.load_chae_base())
    vp = CH.observed_vp(d)
    s_ = d["s"]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw=dict(wspace=0.25))
    ax = axes[0]
    H, xe, ye = np.histogram2d(np.log10(s_), np.log10(vp), bins=[np.linspace(-0.75, 1.5, 120), np.linspace(-2.2, 1.0, 120)])
    ax.pcolormesh(10 ** xe, 10 ** ye, np.log10(np.where(H.T > 0, H.T, np.nan)), cmap=PS.CMAP_CONTEXT, rasterized=True)
    r = np.logspace(-0.8, 1.55, 200)
    gf = Q.galactic_field()
    f = Q.cached_field(gf["gNe_a0"])
    for M, ls in ((0.6, ":"), (1.2, "-"), (2.4, "--")):
        vesc = np.sqrt(2 * KIN.GM_KMS2_AU * M / (r * 1000))
        ax.plot(r, vesc, color=C_NEWTON, linewidth=1.5, linestyle=ls, label=f"Newtonian limit, M = {M} M☉")
    M = 1.2
    rm = 7.029 * np.sqrt(M)
    eta = np.interp(np.log(r / rm), f.lnr, f.eta)
    ax.plot(r, np.sqrt(2 * eta * KIN.GM_KMS2_AU * M / (r * 1000)), color=C_MOND, linewidth=1.8,
            label="MOND limit with the Galactic field, M = 1.2 M☉")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.18, 32)
    ax.set_ylim(10 ** -2.2, 10)
    ax.set_xlabel("projected separation s (kAU)")
    ax.set_ylabel("projected relative velocity v_sky (km/s)")
    ax.legend(loc="lower left", fontsize=7.2)
    ax.set_title(f"Chae (2023)'s {len(s_):,} pairs: speeds fall as s^-1/2", loc="left")
    ax = axes[1]
    vc = np.sqrt(KIN.GM_KMS2_AU * (d["M1"] + d["M2"]) / (s_ * 1000))
    edges = np.logspace(np.log10(0.2), np.log10(30), 13)
    xm, med, q25, q75 = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (s_ > lo) & (s_ < hi)
        xm.append(np.sqrt(lo * hi))
        med.append(np.median(vp[k] / vc[k]))
        q25.append(np.percentile(vp[k] / vc[k], 25))
        q75.append(np.percentile(vp[k] / vc[k], 75))
    ax.fill_between(xm, q25, q75, color=PS.GRID, linewidth=0, label="interquartile range")
    ax.plot(xm, med, "o-", color=C_DATA, markersize=4, label="median v~ = v_sky / sqrt(GM/s), all pairs")
    # the same for the precisely measured half (propagated sigma(v~) < 0.05)
    sig = (KIN.K_PM * np.hypot(np.hypot(d["mu1ra_err"], d["mu2ra_err"]), np.hypot(d["mu1dec_err"], d["mu2dec_err"]))
           / np.sqrt(2) * d["dM"] / 1000 / vc)
    xp, mp_ = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (s_ > lo) & (s_ < hi) & (sig < 0.05)
        if k.sum() >= 40:
            xp.append(np.sqrt(lo * hi))
            mp_.append(np.median(vp[k] / vc[k]))
    ax.plot(xp, mp_, "s-", color=C_THIRD, markersize=4, label="the same, only pairs with sigma(v~) < 0.05")
    ax.axhline(np.median(med[:3]), color=C_NEWTON, linewidth=1.2, label="Newton, if nothing else changes with s")
    ratio = np.sqrt(np.interp(np.log(np.array(xm) * 1.2 / rm), f.lnr, f.eta))
    ax.plot(xm, np.median(med[:3]) * ratio, color=C_MOND, linewidth=1.5, linestyle="--",
            label="MOND: up to ~19% faster (rough, 1.2 M☉)")
    ax.set_xscale("log")
    ax.set_xlim(0.18, 32)
    ax.set_ylim(0.2, 1.3)
    ax.set_xlabel("projected separation s (kAU)")
    ax.set_ylabel("v~")
    ax.legend(loc="upper left", fontsize=7.2)
    ax.set_title("Scaled by the Newtonian speed: the naive test", loc="left")
    PS.caption(fig, "v_sky from the proper-motion difference at the pair's mean distance (no perspective "
                    "correction, as in Chae 2023). A bound Newtonian orbit cannot exceed sqrt(2GM/s) on the sky; "
                    "MOND raises the limit by up to sqrt(1.42). Right: the MOND curve is only indicative "
                    "(it ignores projection and the spread of masses). The rise of the black curve comes mostly "
                    "from poorly measured pairs (section 5); the forward model and the mocks do this properly.")
    return PS.save(fig, FIG_DIR / "fig03_envelope")


# --------------------------------------------------------------------------
# 5b. Median v~ against r_sky / r_M for every sample
# --------------------------------------------------------------------------

def fig_profiles():
    mp = _res("median_profiles")
    if mp is None:
        return None
    panels = [("banik2024", "2.0", "Banik et al. rebuild, v~ < 2"),
              ("cookson2026", "2.5", "Cookson et al. rebuild, v~ < 2.5"),
              ("chae2023", "2.0", "Chae's sample, v~ < 2"),
              ("chae2023_precise", "2.0", "Chae's precise half, v~ < 2")]
    pt = _res("profile_tests") or {}
    from matplotlib.ticker import FuncFormatter, NullFormatter
    plain = FuncFormatter(lambda v, _: f"{v:g}")

    def logx(ax):
        ax.set_xscale("log")
        lo, hi = ax.get_xlim()
        ax.set_xticks([t for t in (0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2) if lo <= t <= hi])
        ax.xaxis.set_major_formatter(plain)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlim(lo, hi)
    fig, axs = plt.subplots(2, 4, figsize=(14.0, 8.6), sharey=False, gridspec_kw=dict(wspace=0.3, hspace=0.5))
    axes = axs[0]
    for ax, (k, cut, title) in zip(axes, panels):
        if k not in mp:
            ax.set_visible(False)
            continue
        c = mp[k]["cuts"][cut]
        for lab, col, ls in (("newton", C_NEWTON, "-"), ("mond", C_MOND, "--")):
            m = c[lab]
            ax.plot(m["x"], m["median"], color=col, linewidth=1.8, linestyle=ls,
                    label="Newton (fitted population)" if lab == "newton" else "same population, MOND")
        dd = c["data"]
        ax.errorbar(dd["x"], dd["median"], yerr=dd["err"], fmt="o", color=C_DATA, markersize=5, capsize=0,
                    elinewidth=1.0, label="data", zorder=5)
        logx(ax)
        ax.set_xlabel("r_sky / r_M")
        q = pt.get("chi2", {}).get(f"{k}|vt<{cut}")
        sub_ = (f"\nχ² ({q['n_bins']} bins): flat {q['flat']['chi2']:.1f}, MOND step {q['mond_step']['chi2']:.1f}"
                if q else "")
        ax.set_title(title + sub_, loc="left", fontsize=8.6)
        if ax is axes[0]:
            ax.set_ylabel("median v~")
            ax.legend(loc="upper left", fontsize=7.2)
    # second row: Chae's sample split by distance and by velocity precision
    import chae as CH
    import kinematics as KIN
    d = CH.chae_sample(CH.load_chae_base())
    vp = CH.observed_vp(d)
    mt = d["M1"] + d["M2"]
    vc = np.sqrt(KIN.GM_KMS2_AU * mt / (d["s"] * 1000))
    vt = vp / vc
    sig = (KIN.K_PM * np.hypot(np.hypot(d["mu1ra_err"], d["mu2ra_err"]), np.hypot(d["mu1dec_err"], d["mu2dec_err"]))
           / np.sqrt(2) * d["dM"] / 1000 / vc)
    x = d["s"] / (7.029 * np.sqrt(mt))
    ref = mp["chae2023"]["cuts"]["2.0"]
    rx = np.array(ref["newton"]["x"])
    rr = np.array(ref["mond"]["median"]) / np.array(ref["newton"]["median"])
    splits = [("Chae's sample by distance", [("d < 80 pc", d["dM"] < 80, C_NEWTON), ("80-140 pc", (d["dM"] > 80) & (d["dM"] < 140), C_THIRD),
                                             ("d > 140 pc", d["dM"] > 140, C_MOND)]),
              ("Chae's sample by velocity error", [("sigma(v~) < 0.05", sig < 0.05, C_NEWTON), ("0.05-0.1", (sig > 0.05) & (sig < 0.1), C_THIRD),
                                             ("> 0.1", sig > 0.1, C_MOND)]),
              ("Chae's precise pairs by distance", [("sigma < 0.05, d < 100 pc", (sig < 0.05) & (d["dM"] < 100), C_NEWTON),
                                              ("sigma < 0.05, d > 100 pc", (sig < 0.05) & (d["dM"] > 100), C_MOND)])]
    edges = np.quantile(x[vt < 2], np.linspace(0, 1, 9))
    for ax, (title, groups) in zip(axs[1][:3], splits):
        for lab, m, col in groups:
            k = m & (vt < 2)
            xm, med, err = [], [], []
            for lo, hi in zip(edges[:-1], edges[1:]):
                kk = k & (x >= lo) & (x < hi)
                if kk.sum() < 25:
                    continue
                xm.append(np.exp(np.mean(np.log(x[kk]))))
                med.append(np.median(vt[kk]))
                err.append(1.2533 * np.std(vt[kk]) / np.sqrt(kk.sum()))
            ax.errorbar(xm, med, yerr=err, fmt="o-", color=col, markersize=4, linewidth=1.2, capsize=0, label=lab)
        logx(ax)
        ax.set_xlabel("r_sky / r_M")
        ax.set_ylabel("median v~ (v~ < 2)")
        ax.set_title(title, loc="left", fontsize=8.6)
        ax.legend(loc="upper left", fontsize=7.2)
    ax = axs[1][3]
    ax.plot(rx, rr, color=C_MOND, linewidth=1.8, linestyle="--", label="MOND / Newton, our forward model")
    ax.axhline(1.0, color=C_NEWTON, linewidth=1.8, label="Newton")
    logx(ax)
    ax.set_xlabel("r_sky / r_M")
    ax.set_ylabel("predicted ratio of medians")
    ax.set_title("The step MOND predicts", loc="left", fontsize=9)
    ax.legend(loc="upper left", fontsize=7.2)
    PS.caption(fig, "Top: equal-count bins in projected separation over MOND radius; curves are our forward model at "
                    "the Newtonian best fit to the Banik-style sample (with hidden companions and, for Chae's sample, "
                    "its noise), with only gravity switched between Newton and MOND; below 2 kAU (r_sky/r_M < ~0.25) "
                    "that population is extrapolated. Bottom: Chae's 26,615 pairs split by distance and by each pair's "
                    "propagated velocity error; gravity cannot depend on either. Error bars: bootstrap (top), "
                    "1.25 sigma/sqrt(N) (bottom).")
    return PS.save(fig, FIG_DIR / "fig05_median_profiles")


# --------------------------------------------------------------------------
# 7. The 3D test
# --------------------------------------------------------------------------

def fig_threed():
    td = _res("three_d")
    if td is None:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), gridspec_kw=dict(wspace=0.25, width_ratios=[1, 1.1]))
    ax = axes[0]
    x = np.array(td["r_perp_kau"])
    y = np.array(td["v_over_vesc"])
    fast = y > 1
    ax.scatter(x[~fast], y[~fast], s=28, color=C_DATA, edgecolor=PS.SURFACE, linewidth=1.2, label="bound possible", zorder=4)
    ax.scatter(x[fast], y[fast], s=34, color=C_MOND, edgecolor=PS.SURFACE, linewidth=1.2,
               label=f"too fast to be bound in Newton ({fast.sum()})", zorder=5)
    ax.axhline(1.0, color=C_NEWTON, linewidth=1.2)
    ax.axhline(np.sqrt(Q_ETA), color=C_MOND, linewidth=1.2, linestyle="--")
    ax.text(0.99, 1.0 - 0.012, "Newtonian escape speed", color=PS.INK_SECONDARY, fontsize=7.5, va="top", ha="right",
            transform=ax.get_yaxis_transform())
    ax.text(0.5, np.sqrt(Q_ETA) + 0.012, "MOND escape speed (√η × Newtonian, far field)", color=PS.INK_SECONDARY,
            fontsize=7.5, va="bottom", ha="center", transform=ax.get_yaxis_transform())
    ax.set_xscale("log")
    ax.set_xticks([3, 5, 10, 20, 30])
    ax.set_xticklabels(["3", "5", "10", "20", "30"])
    ax.minorticks_off()
    ax.set_ylim(0.1, 1.36)
    ax.set_xlabel("projected separation r_⊥ (kAU)")
    ax.set_ylabel("3D relative speed / Newtonian escape speed at r_⊥")
    ax.legend(loc="upper left", fontsize=7.4)
    ax.set_title("36 pairs with precise RVs (Chae et al. 2026)", loc="left")
    ax = axes[1]
    for k, lab, c, ls in (("opik_drop0", "all 36, log-flat P(a)", C_DATA, "-"),
                          ("bpl_drop0", "all 36, broken power-law P(a)", C_DATA, "--"),
                          ("opik_drop6", "without the 6 too-fast pairs", C_THIRD, "-")):
        if k in td["runs"]:
            r = td["runs"][k]
            g = np.array(r["gammas"])
            ax.plot(g, np.array(r["posterior_log10"]), color=c, linewidth=1.8, linestyle=ls,
                    label=f"ours: {lab} ({r['median']:.2f})")
    try:
        import threed as T
        st = np.load(T.RAW_DIR / "chae2026_3d" / "posterior_samples_saad_ting.npz")
        for k, lab, ls in (("baseline", "Saad & Ting: free a, robust errors", "-"),
                           ("geometric_deprojection", "Saad & Ting: geometric", ":")):
            if k in st.files:
                h, e = np.histogram(np.log10(st[k]), bins=60, density=True)
                ax.step(10 ** e[:-1], h, where="post", color=PS.MUTED, linewidth=1.4, linestyle=ls,
                        label=f"{lab} ({np.median(st[k]):.2f})")
    except Exception:  # noqa: BLE001
        pass
    ax.axvline(1.0, color=C_NEWTON, linewidth=1.2)
    ax.axvline(Q_ETA, color=C_MOND, linewidth=1.2, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlim(0.5, 3.5)
    ax.set_xticks([0.6, 0.8, 1.0, 1.5, 2.0, 3.0])
    ax.set_xticklabels(["0.6", "0.8", "1", "1.5", "2", "3"])
    ax.minorticks_off()
    ax.set_ylim(0, 13.8)
    ax.set_xlabel("gravity boost γ = G_eff / G")
    ax.set_ylabel("posterior density per dex")
    ax.legend(loc="upper left", fontsize=7.2, frameon=True, facecolor=PS.SURFACE, edgecolor=PS.SURFACE,
              framealpha=1.0)
    ax.set_title("The answer hinges on a handful of pairs", loc="left")
    PS.caption(fig, "Our semi-analytic likelihood: each pair's phase-space density for an isotropic population of "
                    "Kepler orbits, integrated over the unknown line-of-sight separation and the velocity errors. "
                    "Pairs above the blue line cannot be bound Newtonian orbits at any line-of-sight separation.")
    return PS.save(fig, FIG_DIR / "fig07_three_d")


FIGURES["envelope"] = fig_envelope
FIGURES["profiles"] = fig_profiles
FIGURES["threed"] = fig_threed


if __name__ == "__main__":
    names = sys.argv[1:] or list(FIGURES)
    for n in names:
        print(n, FIGURES[n]())
