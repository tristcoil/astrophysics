"""Every figure of the project, drawn from the analysis results.

Each ``fig_*`` function takes the results dict from ``analysis.run()`` and
returns the paths it wrote (PNG + SVG into ``figures/``).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

import config as C
import lightcurves as LC
import ogle
import plotstyle as ps

ps.use_style()

#: Sequential blue ramp (palette steps 100 -> 700) and the blue <-> red diverging map.
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
DIVERGING = LinearSegmentedColormap.from_list(
    "near_far", ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f3a3a3", "#e34948", "#a32828"])
MODE_STYLE = {"F": {"color": ps.BLUE, "marker": "o", "label": "fundamental mode (F)"},
              "1O": {"color": ps.ORANGE, "marker": "s", "label": "first overtone (1O)"}}
GAL_STYLE = {"lmc": {"color": ps.BLUE, "marker": "o", "label": "LMC"},
             "smc": {"color": ps.ORANGE, "marker": "s", "label": "SMC"}}
BAND_STYLE = {"I": {"color": ps.BLUE, "marker": "o", "label": "I band (0.8 µm)"},
              "V": {"color": ps.ORANGE, "marker": "s", "label": "V band (0.55 µm)"}}


def _out(name):
    return str(C.FIGURE_DIR / name)


# ---------------------------------------------------------------------------
# 1. Light curves and the Hertzsprung progression
# ---------------------------------------------------------------------------

def fig_light_curves(res) -> list:
    ex = res["lc_examples"]
    hp = res["hertzsprung"]
    fig = plt.figure(figsize=(12.5, 7.6))
    gs = fig.add_gridspec(len(ex), 3, width_ratios=[1.25, 1.0, 1.0], hspace=0.12, wspace=0.32)
    for k, e in enumerate(ex):
        ax = fig.add_subplot(gs[k, 0])
        for band in ("V", "I"):
            d = e[band]
            if d is None:
                continue
            st = BAND_STYLE[band]
            dm = d["mag"] - d["mean"]
            mm = d["model"] - d["mean"]
            for shift in (0, 1):
                ax.scatter(d["phase"] + shift, dm, s=4, color=st["color"], alpha=0.35, linewidths=0,
                           marker=st["marker"])
                ax.plot(d["grid"] + shift, mm, color=st["color"], lw=1.4)
        ax.invert_yaxis()
        ax.set_xlim(0, 2)
        ax.text(0.01, 0.05, f"P = {e['P']:.2f} d  ({e['mode']})", transform=ax.transAxes, fontsize=8,
                color=ps.INK, va="bottom")
        ax.text(0.99, 0.05, e["id"].replace("OGLE-", ""), transform=ax.transAxes, fontsize=7,
                color=ps.MUTED, ha="right", va="bottom")
        ax.tick_params(labelsize=7)
        if k < len(ex) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("pulsation phase (two cycles; 0 = maximum light)")
        if k == 0:
            ax.set_title("a  Folded light curves", pad=18)
            handles = [plt.Line2D([], [], color=BAND_STYLE[b]["color"], marker=BAND_STYLE[b]["marker"], ms=4,
                                  lw=1.4, label=BAND_STYLE[b]["label"]) for b in ("V", "I")]
            ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=7.5,
                      handletextpad=0.3, borderaxespad=0.2)
    fig.text(0.075, 0.5, "brightness relative to the mean (mag; brighter is up)", rotation=90, va="center",
             ha="center", color=ps.INK_SECONDARY, fontsize=9)

    # Hertzsprung progression: median normalised I light curve per period bin, stacked
    ax = fig.add_subplot(gs[:, 1])
    grid = hp["grid"]
    step = 0.62
    for j, rung in enumerate(hp["ladder"]):
        y0 = j * step
        for shift in (0, 1):
            ax.plot(grid + shift, y0 + rung["curve"] * 0.55, color=ps.BLUE, lw=1.3)
        ax.text(-0.04, y0 + 0.27, f"{rung['P']:g} d", ha="right", va="center", fontsize=7.5, color=ps.INK)
        ax.text(2.04, y0 + 0.27, f"n={rung['n']}", ha="left", va="center", fontsize=6.5, color=ps.MUTED)
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.15, len(hp["ladder"]) * step + 0.05)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("phase (two cycles; 0 = maximum light)")
    ax.set_title("b  The bump moves with period", pad=18)
    ax.text(0.0, 1.005, "median I-band shape, LMC fundamental mode; up = brighter", transform=ax.transAxes,
            fontsize=7, color=ps.MUTED, va="bottom")
    # point at the two unambiguous bumps: falling branch at 8 d, rising branch at 16 d
    rungs = {r["P"]: (j, r) for j, r in enumerate(hp["ladder"])}
    for p0, text, dx, dy in ((8.0, "bump on the\nfalling branch", 0.38, -0.55), (16.0, "bump on the\nrising branch", -0.62, 0.42)):
        if p0 in rungs and rungs[p0][1]["bump_kind"] == "bump":
            j, r = rungs[p0]
            ph = r["bump_phase"] + (1 if r["bump_phase"] < 0.5 else 0)
            yb = j * step + np.interp(ph % 1.0, grid, r["curve"]) * 0.55
            ps.annotate(ax, text, (ph, yb), (ph + dx, yb + dy), color=ps.INK, fontsize=7)

    # Fourier phase phi21 vs period: the progression and the two modes
    ax = fig.add_subplot(gs[:, 2])
    for mode in ("F", "1O"):
        d = hp["phi21"][mode]
        st = MODE_STYLE[mode]
        ax.scatter(d["logP"], d["phi21"], s=3, color=st["color"], marker=st["marker"], alpha=0.45,
                   linewidths=0, label=st["label"])
    ax.set_xlabel("period (days)")
    ax.set_ylabel("Fourier phase φ₂₁ of the I-band light curve (rad)")
    xt = [0.3, 0.5, 1, 2, 5, 10, 20, 50]
    ax.set_xticks(np.log10(xt))
    ax.set_xticklabels([str(v) for v in xt])
    ax.set_ylim(2.5, 6.6)
    ax.set_xlim(np.log10(0.22), np.log10(60))
    ax.legend(loc="upper left", markerscale=3)
    ax.set_title("c  Light-curve shape sorts the modes", pad=18)
    ps.annotate(ax, "Hertzsprung progression:\nthe bump crosses\nmaximum light", (np.log10(9.5), 5.9),
                (np.log10(13), 6.3), color=ps.INK, fontsize=7.5)
    fig.suptitle("A Cepheid's heartbeat: ten years of OGLE photometry, folded on the pulsation period", y=1.0, x=0.08,
                 ha="left")
    ps.caption(fig, "Data: OGLE-IV Collection of Variable Stars (Soszyński et al. 2015-2023), LMC classical Cepheids; "
               "our own Fourier fits to the raw I- and V-band photometry (2010-2020). Panel b: each curve is the "
               "median normalised shape of all stars within ±5% of the period shown (n stars); the secondary bump "
               "rides the falling branch below ~9 days, merges with maximum light near 10 days and rides the rising "
               "branch above.", y=0.02)
    return ps.save(fig, _out("01_light_curves"))


# ---------------------------------------------------------------------------
# 3. The Leavitt law in V, I and W_I (LMC), and the Wesenheit's physics
# ---------------------------------------------------------------------------

def _logp_axis(ax, ticks=(1, 2, 5, 10, 20, 50, 100)):
    ax.set_xticks(np.log10(ticks))
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel("pulsation period (days)")


def fig_leavitt_laws(res) -> list:
    lmc = res["tabs"]["lmc"]
    laws = res["laws"]["lmc"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.6), sharex=True)
    bands = (("V", "V magnitude (0.55 µm)"), ("I", "I magnitude (0.8 µm)"),
             ("W_I", "Wesenheit W_I = I − 1.55 (V − I)"))
    for ax, (band, label) in zip(axes, bands):
        for mode in ("1O", "F"):
            st = MODE_STYLE[mode]
            f = laws[(mode, band)]
            m = f["mask"]
            sel = (lmc["mode"] == mode) & np.isfinite(lmc[band])
            out = sel & ~m
            ax.scatter(lmc["logP"][m], lmc[band][m], s=3, color=st["color"], marker=st["marker"], alpha=0.35,
                       linewidths=0, label=st["label"] if band == "V" else None, rasterized=True)
            ax.scatter(lmc["logP"][out], lmc[band][out], s=4, color=ps.CONTEXT, marker="x", linewidths=0.5,
                       label="clipped (> 3σ)" if (band == "V" and mode == "F") else None, rasterized=True)
            lp = np.linspace(np.nanmin(lmc["logP"][m]), np.nanmax(lmc["logP"][m]), 50)
            ax.plot(lp, f["a"] + f["b"] * (lp - 1), color=ps.INK, lw=1.0)
        fF, f1 = laws[("F", band)], laws[("1O", band)]
        ax.text(0.97, 0.03, f"scatter about the law\nF:  {fF['sigma']:.3f} mag\n1O: {f1['sigma']:.3f} mag\n"
                            f"slope (F): {fF['b']:.3f} ± {fF['eb']:.3f}",
                transform=ax.transAxes, fontsize=8, color=ps.INK, va="bottom", ha="right")
        ax.invert_yaxis()
        ax.set_ylabel(label)
        _logp_axis(ax, (0.3, 1, 3, 10, 30, 100))
        ax.set_xlim(np.log10(0.2), np.log10(150))
    axes[0].set_title("a  V band: blurred by dust")
    axes[1].set_title("b  I band: less dust")
    axes[2].set_title("c  Wesenheit: dust-free")
    axes[0].legend(loc="upper right", markerscale=3, fontsize=7.5)
    fig.suptitle(f"The Leavitt law in the Large Magellanic Cloud: {np.sum(lmc['mode'] == 'F'):,} fundamental-mode and "
                 f"{np.sum(lmc['mode'] == '1O'):,} first-overtone Cepheids", x=0.06, ha="left", y=1.01)
    ps.caption(fig, "Data: OGLE-IV (Soszyński et al. 2015-2023); mean magnitudes are our own intensity means from the "
               "raw light curves. Lines: least-squares fits with iterative 3σ clipping, m = a + b (log P − 1). The "
               "brighter a Cepheid, the longer its period: magnitudes run downward because smaller magnitude means brighter.",
               y=-0.03)
    return ps.save(fig, _out("03_leavitt_law_lmc"))


def fig_wesenheit(res) -> list:
    """Colour-magnitude diagram: the strip, lines of constant period and the dust vector."""
    lmc = res["tabs"]["lmc"]
    w = res["wesenheit"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.8), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    s = (lmc["mode"] == "F") & np.isfinite(lmc["VI"]) & np.isfinite(lmc["I"])
    sc = ax.scatter(lmc["VI"][s], lmc["I"][s], c=lmc["logP"][s], s=4, cmap=SEQ_BLUE, vmin=0.0, vmax=1.6,
                    linewidths=0, rasterized=True)
    # lines of constant period from the period-luminosity-colour fit: I = a + b (logP - 1) + c (V - I)
    plc = w["plc"]
    vi = np.linspace(0.45, 1.45, 10)
    for P in (2, 5, 10, 20, 40):
        I = plc["a"] + plc["b"] * (np.log10(P) - 1) + plc["c"] * vi
        ax.plot(vi, I, color=ps.INK_SECONDARY, lw=0.8)
        ax.text(vi[-1] + 0.01, I[-1], f"{P} d", fontsize=7, color=ps.INK_SECONDARY, va="center")
    # dust vector: E(V-I) = 0.3 moves a star 0.3 redder and 1.55 x 0.3 fainter
    x0, y0 = 0.5, 12.0
    ax.annotate("", xy=(x0 + 0.3, y0 + 1.55 * 0.3), xytext=(x0, y0),
                arrowprops={"arrowstyle": "-|>", "color": ps.ORANGE, "lw": 1.6})
    ax.text(x0 + 0.02, y0 - 0.12, "dust: E(V−I) = 0.3\n(A_I = 1.55 E)", color=ps.INK, fontsize=7.5, va="bottom")
    ax.invert_yaxis()
    ax.set_xlim(0.35, 1.6)
    ax.set_ylim(17.2, 11.3)
    ax.set_xlabel("colour V − I (mag; redder to the right)")
    ax.set_ylabel("I magnitude")
    cb = fig.colorbar(sc, ax=ax, pad=0.01)
    cb.set_label("period (days)")
    cb.set_ticks(np.log10([1, 2, 5, 10, 20, 40]))
    cb.set_ticklabels(["1", "2", "5", "10", "20", "40"])
    ax.set_title("The instability strip: at one period, bluer = brighter")
    ax.text(0.02, 0.02, f"constant-period lines: I = … + {plc['c']:.2f} (V − I)\ndust vector: slope 1.55",
            transform=ax.transAxes, fontsize=8, color=ps.INK)
    # right: residual from the PL law in I vs the colour residual
    ax = axes[1]
    fI = res["laws"]["lmc"][("F", "I")]
    pc = w["period_colour"]
    dI = lmc["I"] - (fI["a"] + fI["b"] * (lmc["logP"] - 1))
    dc = lmc["VI"] - (pc["a"] + pc["b"] * (lmc["logP"] - 1))
    s = (lmc["mode"] == "F") & np.isfinite(dI) & np.isfinite(dc) & (np.abs(dI) < 0.8) & (np.abs(dc) < 0.5)
    ax.scatter(dc[s], dI[s], s=3, color=ps.BLUE, alpha=0.35, linewidths=0, rasterized=True)
    x = np.linspace(-0.35, 0.45, 10)
    ax.plot(x, plc["c"] * x, color=ps.INK, lw=1.2, label=f"best fit: slope {plc['c']:.2f} ± {plc['err'][2]:.2f}")
    ax.plot(x, 1.55 * x, color=ps.ORANGE, lw=1.2, label="dust: slope 1.55 (the Wesenheit coefficient)")
    ax.invert_yaxis()
    ax.set_xlabel("colour offset from the period-colour relation (mag)")
    ax.set_ylabel("brightness offset from the I-band Leavitt law (mag)")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("Why W_I works: both blurs point the same way")
    ax.text(0.98, 0.03, f"σ(I) = {fI['sigma']:.3f} → σ(W_I) = {res['laws']['lmc'][('F', 'W_I')]['sigma']:.3f} mag",
            transform=ax.transAxes, ha="right", fontsize=8, color=ps.INK)
    fig.suptitle("Why the Wesenheit index tightens the law (LMC fundamental-mode Cepheids)", x=0.06, ha="left", y=1.0)
    ps.caption(fig, "A Cepheid that is redder than average at its period is either cooler (intrinsically fainter) or behind "
               "more dust (dimmed): both move it down and to the right along nearly the same line. The Wesenheit "
               "index measures brightness along that line, so both blurs cancel. Data: OGLE-IV; our intensity means.",
               y=-0.02)
    return ps.save(fig, _out("04_why_wesenheit"))


# ---------------------------------------------------------------------------
# 5. SMC vs LMC; 6. the shapes of the Clouds
# ---------------------------------------------------------------------------

def fig_smc_lmc(res) -> list:
    import literature as LIT
    tabs, sl, geo = res["tabs"], res["smc_lmc"], res["geometry"]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.4), gridspec_kw={"width_ratios": [1.2, 1, 0.9], "wspace": 0.45})
    ax = axes[0]
    for gal in ("lmc", "smc"):
        t = tabs[gal]
        st = GAL_STYLE[gal]
        s = (t["mode"] == "F") & (t["logP"] > 0.4) & np.isfinite(t["W_I"])
        ax.scatter(t["logP"][s], t["W_I"][s] - geo[gal]["geo_corr"][s], s=3, color=st["color"], marker=st["marker"],
                   alpha=0.3, linewidths=0, label=f"{st['label']} ({s.sum():,} F-mode)", rasterized=True)
    fl = res["laws"]["lmc"][("F", "W_I")]
    lp = np.linspace(0.4, 1.9, 20)
    d = sl["F"]["smc_central_1deg"]["dmu"]
    ax.plot(lp, fl["a"] + fl["b"] * (lp - 1), color=ps.INK, lw=1.0)
    ax.plot(lp, fl["a"] + d + fl["b"] * (lp - 1), color=ps.INK, lw=1.0)
    x = 1.55
    y1 = fl["a"] + fl["b"] * (x - 1)
    ax.annotate("", xy=(x, y1 + d), xytext=(x, y1), arrowprops={"arrowstyle": "<->", "color": ps.INK, "lw": 0.9})
    ax.text(x + 0.03, y1 + d + 0.35, f"Δ = {d:.3f} mag:\nthe SMC is {10 ** (0.2 * d):.2f}×\nas far away", fontsize=8,
            va="top")
    ax.invert_yaxis()
    _logp_axis(ax, (3, 10, 30, 100))
    ax.set_ylabel("Wesenheit W_I (mag), moved to each Cloud's centre")
    ax.legend(loc="upper left", markerscale=3, fontsize=7.5)
    ax.set_title("a  One law, two distances")
    # residual histograms: the SMC is deeper
    ax = axes[1]
    bins = np.linspace(-0.6, 0.6, 61)
    for gal, key in (("lmc", "sigma_lmc"), ("smc", "sigma_smc")):
        t = tabs[gal]
        r = geo[gal]["dmu"]
        s = (t["mode"] == "F") & np.isfinite(r)
        st = GAL_STYLE[gal]
        ax.hist(r[s], bins=bins, histtype="step", color=st["color"], lw=1.6, density=True,
                label=f"{st['label']}: σ = {np.nanstd(r[s]):.3f} mag")
    ax.set_xlabel("offset from the law (mag; − = nearer)")
    ax.set_ylabel("fraction of stars per mag")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title("b  The SMC is deep")
    ax.text(0.03, 0.80, f"extra spread ⇒ depth\nσ ≈ {sl['depth_sigma_kpc']:.1f} kpc along\nour line of sight",
            transform=ax.transAxes, fontsize=8, va="top")
    # gamma from this test, against the literature
    ax = axes[2]
    rows = [("this work: SMC−LMC, W_I", sl["F"]["smc_central_1deg"]["gamma"], sl["F"]["smc_central_1deg"]["gamma_err"],
             "ours"),
            ("this work: SMC−LMC, Gaia W_G", sl["gaia_WG_R1.90"]["gamma"], sl["gaia_WG_R1.90"]["gamma_err"], "ours"),
            ("this work: MW parallaxes, W_G", res["gaia"]["baseline"]["gamma"]["median"]
             - res["gaia"]["injection"]["gamma_bias"], res["gaia"]["baseline"]["gamma"]["std"], "ours")]
    for key in ("wielgorski17_WI", "breuval21_WVI", "breuval22_WVI", "riess22_WH", "breuval24_WH"):
        v, e, lab = LIT.GAMMA_OTHER[key]
        rows.append((lab.split(",")[0] + ", " + lab.split(",")[1].strip().split(" ")[0], v, e, "lit"))
    for key in ("breuval22", "ripepi22", "trentin24", "cmetall9"):
        g = LIT.GAIA_PWZ[key]
        rows.append((g["ref"].split("(")[0].strip() + ", W_G", g["gamma"], g["e_gamma"], "lit"))
    for j, (lab, v, e, kind) in enumerate(rows):
        col = ps.BLUE if kind == "ours" else ps.INK_SECONDARY
        ax.errorbar(v, -j, xerr=e, fmt="o" if kind == "ours" else "s", color=col, ms=5, capsize=2, lw=1.2)
    ax.axvline(0, color=ps.AXIS, lw=0.8)
    ax.set_xlim(-0.8, 0.15)
    ax.set_yticks(-np.arange(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    for tl, r in zip(ax.get_yticklabels(), rows):
        tl.set_color(ps.INK if r[3] == "ours" else ps.INK_SECONDARY)
    ax.set_xlabel("γ (mag/dex); < 0: metal-rich are brighter")
    ax.set_title("c  Do metal-poor Cepheids shine less?")
    ax.grid(axis="y", visible=False)
    fig.suptitle("The Small Magellanic Cloud against the Large: relative distance, depth and the metallicity term",
                 x=0.06, ha="left", y=1.01)
    ps.caption(fig, f"Fundamental-mode Cepheids longer than 2.5 days (the SMC's law bends below that); SMC within 1° of its "
               f"centre for the distance. Eclipsing binaries put the SMC {LIT.SMC_DEB['dmu_smc_lmc']:.3f} ± "
               f"{LIT.SMC_DEB['dmu_err']:.3f} mag beyond the LMC (Graczyk et al. 2020); with [Fe/H] = −0.75 (SMC) and "
               f"−0.41 (LMC), any extra offset of the Cepheids measures γ. Panel c: grey = literature values in various "
               f"Wesenheit indices (they need not agree exactly).", y=-0.03)
    return ps.save(fig, _out("05_smc_vs_lmc"))


def fig_maps(res) -> list:
    import geometry as GEO
    tabs, geo = res["tabs"], res["geometry"]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), gridspec_kw={"wspace": 0.38})
    for ax, gal, centre, lim, title in ((axes[0], "lmc", (80.05, -69.30), 6.2, "a  LMC: a tilted disk"),
                                        (axes[1], "smc", (12.54, -73.11), 3.3, "b  SMC: stretched along our line of sight")):
        t = tabs[gal]
        g = geo[gal]
        x, y, _, _, _ = GEO.sky_offsets(t["ra"], t["dec"], *centre)
        X, Y = np.degrees(x), np.degrees(y)
        r = g["dmu"]
        ok = np.isfinite(r)
        # median residual on a grid (stars are noisy individually)
        nb = 26
        edges = np.linspace(-lim, lim, nb + 1)
        H = np.full((nb, nb), np.nan)
        for i in range(nb):
            for j in range(nb):
                m = ok & (X >= edges[i]) & (X < edges[i + 1]) & (Y >= edges[j]) & (Y < edges[j + 1])
                if m.sum() >= 4:
                    H[j, i] = np.median(r[m])
        vmax = 0.12 if gal == "lmc" else 0.3
        im = ax.imshow(H, origin="lower", extent=[-lim, lim, -lim, lim], cmap=DIVERGING,
                       norm=TwoSlopeNorm(0, -vmax, vmax), interpolation="nearest")
        ax.scatter(X[~ok], Y[~ok], s=1, color=ps.CONTEXT, linewidths=0)
        # fitted plane: contours of distance offset
        gx, gy = np.meshgrid(np.linspace(-lim, lim, 200), np.linspace(-lim, lim, 200))
        rho = np.radians(np.hypot(gx, gy))
        pa = np.arctan2(gx, gy)
        dm = 5 * np.log10(GEO.disk_distance_ratio(rho, pa, g["disk"]["incl"], g["disk"]["theta"]))
        step = 0.05 if gal == "lmc" else 0.1
        levels = [v for v in np.arange(-0.3, 0.31, step) if abs(v) > 1e-9]
        inside = np.hypot(gx, gy) < 0.92 * lim
        cs = ax.contour(gx, gy, np.where(inside, dm, np.nan), levels=levels, colors=ps.INK, linewidths=0.6)
        ax.clabel(cs, fmt=lambda v: f"{v:+.2f}" if gal == "lmc" else f"{v:+.1f}", fontsize=6.5)
        # line of nodes
        th = np.radians(g["disk"]["theta"])
        ax.plot([-lim * np.sin(th), lim * np.sin(th)], [-lim * np.cos(th), lim * np.cos(th)], color=ps.INK_SECONDARY,
                lw=0.8)
        ax.set_xlim(lim, -lim)  # east to the left, as on the sky
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("← east    offset from the centre (deg)    west →")
        ax.set_ylabel("offset (deg); north up")
        ax.set_title(title)
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label("median offset from the law (mag; blue = nearer)")
        d = g["disk"]
        ax.text(0.02, 0.02, f"plane fit: i = {d['incl']:.1f}°, nodes at PA {d['theta']:.0f}°\nnear side toward PA "
                            f"{d['near_side_pa']:.0f}° (north-east)", transform=ax.transAxes, fontsize=7.5,
                color=ps.INK, bbox={"facecolor": ps.SURFACE, "edgecolor": "none", "alpha": 0.8, "pad": 2})
    fig.suptitle("Where the Clouds lie nearer and farther: Cepheid offsets from the Leavitt law (W_I, F and 1O)",
                 x=0.06, ha="left", y=1.0)
    ps.caption(fig, "Each pixel: the median offset of the Cepheids in it (0.1 mag = 4.7% in distance). Contours: our "
                    "fitted inclined plane (mag relative to the centre); grey line: its line of nodes. LMC literature: i = 24.2° ± 0.7°, PA 151.4° ± 1.7° "
                    "(Jacyszyn-Dobrzeniecka et al. 2016); i = 25.05°, PA 150.8° (Inno et al. 2016); north-east nearer. "
                    "The SMC is not a disk: the 'plane' only captures its north-east-near, south-west-far elongation.",
               y=-0.02)
    return ps.save(fig, _out("06_magellanic_geometry"))


# ---------------------------------------------------------------------------
# 7. The Gaia zero point; 8. the LMC distance and H0
# ---------------------------------------------------------------------------

def fig_gaia_calibration(res) -> list:
    g = res["gaia"]
    pts = g["points"]
    b = g["baseline"]
    mu = b["mu"]["median"]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2), gridspec_kw={"hspace": 0.38, "wspace": 0.28})
    plx, e = pts["plx"], pts["e_plx"] * pts["k"]
    snr = plx / e
    hi = pts["kept"] & (snr > 10)
    M = pts["W"] + 5 * np.log10(np.where(plx > 0, plx, np.nan) / 100.0)
    eM = 2.5 / np.log(10) / snr
    # a: absolute Wesenheit from parallaxes; the LMC law shifted by its fitted distance
    ax = axes[0, 0]
    for mode in ("F", "1O"):
        st = MODE_STYLE[mode]
        s = hi & (pts["mode"] == mode)
        ax.errorbar(pts["logP"][s], M[s], yerr=eM[s], fmt=st["marker"], ms=3, color=st["color"], alpha=0.55,
                    elinewidth=0.6, capsize=0, label=f"{st['label']}: {s.sum()} with parallax S/N > 10")
    shp = g["lmc_shapes"]
    R = b["R"]["median"] if "R" in b else C.R_WG
    for mode, sh in shp.items():
        pass
    lp = np.linspace(0.0, 1.8, 30)
    # model curves at the median calibrator metallicity
    feh_med = float(np.nanmedian(pts["feh"][hi]))
    for mode in ("F", "1O"):
        s = pts["kept"] & (pts["mode"] == mode)
        # the model absolute magnitude depends on metallicity; fit a line to the model values for drawing
        coef = np.polyfit(pts["logP"][s] - 1, pts["M_model"][s] - b["gamma"]["median"] * (pts["feh"][s] - feh_med), 1)
        x = lp if mode == "F" else lp[lp < 0.8]
        ax.plot(x, np.polyval(coef, x - 1), color=ps.INK, lw=1.1)
    ax.invert_yaxis()
    _logp_axis(ax, (1, 2, 5, 10, 20, 50))
    ax.set_ylabel(f"absolute Wesenheit M_W (mag), W = G − {R:.2f} (BP − RP)")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("a  Milky Way Cepheids weighed by parallax")
    ax.text(0.97, 0.04, f"lines: the LMC's own law, shifted by\nthe fitted μ_LMC = {mu:.3f} (raw fit)\n"
                        f"at [Fe/H] = {feh_med:+.2f}",
            transform=ax.transAxes, ha="right", fontsize=8)
    # b: residuals vs metallicity
    resid = M - pts["M_model"] + b["gamma"]["median"] * (pts["feh"] - (-0.409))
    ax = axes[0, 1]
    s = hi
    ax.scatter(pts["feh"][s], resid[s], s=9, color=ps.BLUE, alpha=0.6, linewidths=0)
    x = np.linspace(-1.1, 0.5, 10)
    gm = b["gamma"]["median"]
    gcor = gm - g["injection"]["gamma_bias"]
    ax.plot(x, gm * (x + 0.409) - 0.0, color=ps.INK, lw=1.1,
            label=f"fit: γ = {gm:.2f} ± {b['gamma']['std']:.2f} mag/dex ({gcor:.2f} after bias correction)")
    ax.axvline(-0.409, color=ps.ORANGE, lw=1.0)
    ax.text(-0.39, -0.72, "LMC Cepheids:\n[Fe/H] = −0.41", color=ps.INK, fontsize=7.5, va="top")
    ax.invert_yaxis()
    ax.set_ylim(0.8, -0.8)
    ax.set_xlabel("spectroscopic [Fe/H] of the Milky Way Cepheid (dex)")
    ax.set_ylabel("M_W − LMC law at LMC metallicity (mag)")
    ax.legend(loc="lower left", fontsize=7.5)
    ax.set_title("b  Metal-rich Cepheids are brighter")
    # c: residuals vs reddening (the dust coefficient)
    ax = axes[1, 0]
    r2 = M - pts["M_model"]
    ax.scatter(pts["E"][s], r2[s], s=9, color=ps.BLUE, alpha=0.6, linewidths=0, label=f"fitted R = {R:.3f}")
    # what the residuals would look like with R = 1.90 instead
    r190 = r2 - (C.R_WG - R) * (pts["E"])  # approximate: W changes by -(dR) * (BP-RP), mostly via E
    ax.scatter(pts["E"][s], r190[s] + 0.0, s=5, color=ps.ORANGE, alpha=0.45, linewidths=0, marker="s",
               label="same stars with R = 1.90")
    ax.axhline(0, color=ps.INK, lw=0.8)
    ax.invert_yaxis()
    ax.set_ylim(0.8, -0.8)
    ax.set_xlabel("reddening E(BP − RP) of the Milky Way Cepheid (mag)")
    ax.set_ylabel("residual (mag)")
    ax.legend(loc="lower left", fontsize=7.5)
    ax.set_title("c  Dust: Gaia's broad bands need R ≈ 1.74, not 1.90")
    # d: posterior of mu vs parallax offset
    ax = axes[1, 1]
    ch = g["chain"]
    inj = g["injection"]
    ax.hist2d((ch["delta"] - inj["delta_bias"]) * 1e3, ch["mu"] - inj["mu_bias"], bins=60, cmap=SEQ_BLUE, cmin=1)
    ax.axhline(18.477, color=ps.ORANGE, lw=1.2)
    ax.axhspan(18.477 - 0.026, 18.477 + 0.026, color=ps.ORANGE, alpha=0.12, lw=0)
    ax.text(ax.get_xlim()[0] + 0.3, 18.477 + 0.03, "eclipsing binaries (Pietrzyński+19)", color=ps.INK, fontsize=7.5)
    import literature as LIT
    for key, lab in (("riess21", "Riess+21"), ("cra23", "Cruz Reyes & Anderson 23"), ("cmetall9", "C-MetaLL IX")):
        v = LIT.PARALLAX_OFFSET[key]["delta_uas"]
        ax.axvline(v, color=ps.INK_SECONDARY, lw=0.7)
        ax.text(v, ax.get_ylim()[1], lab, rotation=90, fontsize=6.5, va="top", ha="right", color=ps.INK_SECONDARY)
    ax.set_xlabel("residual parallax offset δ (µas; true = Lindegren-corrected − δ)")
    ax.set_ylabel("LMC distance modulus μ")
    ax.set_title("d  Distance and parallax offset trade off")
    ax.text(0.98, 0.03, "posterior, shifted by the bias\nfound in the injection test", transform=ax.transAxes,
            ha="right", fontsize=7, color=ps.MUTED)
    ax.grid(False)
    fig.suptitle(f"The zero point: {b['n']} Milky Way Cepheids with Gaia DR3 parallaxes and spectroscopic metallicities",
                 x=0.06, ha="left", y=0.98)
    ps.caption(fig, "The fit is done in parallax space (no parallax cuts, no Lutz-Kelker bias); panels a-c show the stars "
               "with S/N > 10 only, for clarity. Parallaxes: Gaia DR3 with the Lindegren et al. (2021) zero point and a "
               f"fitted residual offset δ = {1e3 * (b['delta']['median'] - g['injection']['delta_bias']):.1f} ± "
               f"{1e3 * b['delta']['std']:.1f} µas (bias-corrected); catalogue errors "
               f"inflated by k = {b['k']['median']:.2f} (fitted). Metallicities: C-MetaLL IX and Trentin et al. 2024 (high-resolution "
               "spectra only). Photometry: Gaia light-curve (intensity) means.", y=0.03)
    return ps.save(fig, _out("07_gaia_zero_point"))


def fig_lmc_distance(res) -> list:
    import analysis
    import literature as LIT
    s = analysis.summarize(res)
    d = s["lmc_distance"]
    fig, ax = plt.subplots(figsize=(11, 7.4))
    rows = [("THIS WORK: Milky Way Cepheid parallaxes (Gaia DR3)", d["mu"], d["total"], "ours")]
    var = s["gaia"]["variants"]
    for name in ("fundamental mode only", "first overtones only", "C-MetaLL IX metallicities only",
                 "OGLE W_I route (literature V, I; R fitted)", "dust coefficient R fixed at 1.90",
                 "catalogue parallax errors (k = 1)", "parallax offset fixed: Riess+21 (14 uas)",
                 "gamma fixed: Breuval+22 (-0.384)", "gamma fixed: C-MetaLL IX (-0.527)",
                 "no metallicity term", "no residual parallax offset (L21 only)"):
        rows.append(("   variant: " + name, var[name]["mu"], np.nan, "variant"))
    rows += [("Eclipsing binaries (Pietrzyński et al. 2019)", LIT.LMC_DEB["mu"], np.hypot(LIT.LMC_DEB["stat"],
                                                                                         LIT.LMC_DEB["sys"]), "lit"),
             ("Cepheids, Gaia W_G (Ripepi et al. 2022)", LIT.GAIA_PWZ["ripepi22"]["mu_lmc"],
              LIT.GAIA_PWZ["ripepi22"]["e_mu_lmc"], "lit"),
             ("Cepheids, Gaia W_G (Trentin et al. 2024)", LIT.GAIA_PWZ["trentin24"]["mu_lmc"],
              LIT.GAIA_PWZ["trentin24"]["e_mu_lmc"], "lit"),
             ("Cepheids, Gaia W_G (C-MetaLL IX, 2026)", LIT.GAIA_PWZ["cmetall9"]["mu_lmc"],
              LIT.GAIA_PWZ["cmetall9"]["e_mu_lmc"], "lit"),
             ("What Planck's H0 would need", s["h0"]["mu_needed_for_planck"], np.nan, "planck")]
    for j, (lab, v, e, kind) in enumerate(rows):
        y = -j
        if kind == "ours":
            ax.errorbar(v, y, xerr=e, fmt="o", ms=7, color=ps.BLUE, capsize=3, lw=1.6)
            ax.errorbar(v, y, xerr=d["stat"], fmt="none", color=ps.BLUE, capsize=0, lw=3.5)
            ax.plot(d["mu_raw"], y, marker="o", ms=7, mfc="none", mec=ps.BLUE, mew=1.2)
            ps.annotate(ax, "raw fit, before the\nbias correction", (d["mu_raw"], y), (d["mu_raw"] + 0.06, y - 1.3),
                        fontsize=7)
        elif kind == "variant":
            ax.plot(v, y, marker="o", ms=4, color=ps.CONTEXT, mec=ps.INK_SECONDARY, mew=0.6)
        elif kind == "lit":
            ax.errorbar(v, y, xerr=e, fmt="s", ms=5, color=ps.INK_SECONDARY, capsize=2, lw=1.2)
        else:
            ax.plot(v, y, marker="D", ms=6, color=ps.ORANGE)
    ax.axvspan(LIT.LMC_DEB["mu"] - 0.026, LIT.LMC_DEB["mu"] + 0.026, color=ps.ORANGE, alpha=0.10, lw=0)
    ax.set_yticks(-np.arange(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    for tl, r in zip(ax.get_yticklabels(), rows):
        tl.set_color(ps.INK if r[3] in ("ours", "lit", "planck") else ps.MUTED)
    ax.set_xlabel("LMC distance modulus μ (mag)")
    ax.grid(axis="y", visible=False)
    sec = ax.secondary_xaxis("top", functions=(lambda m: 10 ** (0.2 * m + 1) / 1e3, lambda dk: 5 * np.log10(dk * 1e3) - 5))
    sec.set_xlabel("distance (kpc)")
    ax.set_xlim(18.2, 18.8)
    h = s["h0"]
    ax.text(0.01, 0.015, f"H0 with the LMC as the only anchor (SH0ES 2022 rescaled):\n"
                        f"{h['H0_lmc_anchor_ours']:.1f} ± {h['H0_err_from_mu']:.1f} km/s/Mpc with our distance;\n"
                        f"Planck's 67.4 needs the LMC at {h['d_needed_kpc']:.1f} kpc ({h['sigma_from_needed']:.1f}σ from ours)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, color=ps.INK,
            bbox={"facecolor": ps.SURFACE, "edgecolor": ps.GRID, "pad": 4})
    ax.set_title(f"How far is the Large Magellanic Cloud?  This work: {d['d_kpc']:.1f} ± {d['d_err_kpc']:.1f} kpc "
                 f"(μ = {d['mu']:.3f} ± {d['stat']:.3f} stat ± {d['sys']:.3f} sys)", pad=34)
    ps.caption(fig, "Thick blue bar: statistical error; thin: statistical and systematic combined (the RMS shift of seven "
               f"alternative assumptions, plus the bias correction's dependence on the assumed truth). Every fit is corrected for "
               f"its own bias measured on simulated data built from the real calibrators (regression dilution: the raw fit "
               f"returns μ {d['bias']:+.3f} too far). Grey dots: single-change variants (the last two drop physics the data "
               "require). Orange band: the eclipsing-binary distance, a method independent of ours. Rescaling H0 assumes the "
               "rest of the SH0ES ladder is unchanged.", y=0.0)
    return ps.save(fig, _out("08_lmc_distance"))


def fig_petersen(res) -> list:
    pt = res["petersen"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    styles = {"mw": {"color": ps.AQUA, "marker": "D", "label": "Milky Way ([Fe/H] ≈ 0)"},
              "lmc": GAL_STYLE["lmc"] | {"label": "LMC ([Fe/H] ≈ −0.41)"},
              "smc": GAL_STYLE["smc"] | {"label": "SMC ([Fe/H] ≈ −0.75)"}}
    for ax, kind, xlab, ylab, title in (
            (axes[0], "F1O", "fundamental period P_F (days)", "period ratio P_1O / P_F",
             "a  Double-mode F/1O Cepheids"),
            (axes[1], "1O2O", "first-overtone period P_1O (days)", "period ratio P_2O / P_1O",
             "b  Double-mode 1O/2O Cepheids")):
        for key in ("mw", "lmc", "smc"):
            d = pt[key][kind]
            st = styles[key]
            ax.scatter(np.log10(d["P"]), d["ratio"], s=14, color=st["color"], marker=st["marker"], alpha=0.7,
                       linewidths=0.4, edgecolors=ps.SURFACE, label=f"{st['label']}: median {d['median']:.4f} (n={d['n']})")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ticks = (0.3, 0.5, 1, 2, 3, 5) if kind == "F1O" else (0.2, 0.3, 0.5, 1, 1.5)
        ax.set_xticks(np.log10(ticks))
        ax.set_xticklabels([str(t) for t in ticks])
        ax.set_title(title)
        ax.legend(loc="lower left" if kind == "F1O" else "upper right", fontsize=7.5)
    axes[0].set_ylim(0.69, 0.75)
    axes[1].set_ylim(0.795, 0.815)
    fig.suptitle("Petersen diagram: the ratio of two pulsation periods sorts Cepheids by metal content",
                 x=0.06, ha="left", y=1.02)
    ps.caption(fig, "Periods from OGLE (Magellanic Clouds; Soszyński et al. 2015-2023) and OGLE's list of Galactic "
               "Cepheids (Pietrukowicz et al. 2021, 2026 update). The ratio of two pulsation periods depends on how sound "
               "travels through the star, which depends on the opacity of its envelope and so on its metal content: at a "
               "given period the F/1O ratio rises by ~0.01 from the Milky Way to the LMC and again to the SMC.",
               y=-0.04)
    return ps.save(fig, _out("09_petersen_diagram"))


# ---------------------------------------------------------------------------
# 2. Leavitt's 25 stars
# ---------------------------------------------------------------------------

def fig_leavitt(res) -> list:
    lv = res["leavitt"]
    t = lv["table"]
    smc = res["tabs"]["smc"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.4), gridspec_kw={"width_ratios": [1, 1.15, 1], "wspace": 0.3})
    # a: her Figure 2, redrawn from her Table I
    ax = axes[0]
    for key, st, lab in (("max", {"color": ps.BLUE, "marker": "o"}, "at maximum light"),
                         ("min", {"color": ps.ORANGE, "marker": "s"}, "at minimum light")):
        f = lv[f"fit_{key}"]
        ax.scatter(t["logP"], t[key], s=26, color=st["color"], marker=st["marker"], edgecolors=ps.SURFACE,
                   linewidths=0.6, zorder=3, label=f"{lab}: slope {f['b']:.2f} ± {f['eb']:.2f}")
        x = np.array([0.0, 2.2])
        ax.plot(x, f["a"] + f["b"] * (x - 1), color=st["color"], lw=1.1)
    ax.invert_yaxis()
    ax.set_xlim(-0.05, 2.25)
    ax.set_xlabel("log₁₀ (period in days)  — as Leavitt plotted it")
    ax.set_ylabel("photographic magnitude (Harvard provisional scale)")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("a  Leavitt 1912: 25 stars")
    ax.text(0.97, 0.03, "\"A straight line can readily be drawn\namong each of the two series of points\"\n"
                        "— H. S. Leavitt, Harvard Circular 173",
            transform=ax.transAxes, ha="right", fontsize=7.5, style="italic", color=ps.INK_SECONDARY)
    # b: OGLE today
    ax = axes[1]
    for mode in ("1O", "F"):
        st = MODE_STYLE[mode]
        s = (smc["mode"] == mode) & np.isfinite(smc["V"])
        ax.scatter(smc["logP"][s], smc["V"][s], s=2, color=st["color"], marker=st["marker"], alpha=0.3, linewidths=0,
                   label=f"{st['label']} ({s.sum():,})", rasterized=True)
    m = np.isfinite(t["V"])
    ax.scatter(t["logP"][m], t["V"][m], s=60, marker="*", color=ps.INK, edgecolors=ps.SURFACE, linewidths=0.6,
               zorder=4, label="Leavitt's 25 (OGLE V today; 2 saturated)")
    ax.invert_yaxis()
    ax.set_xlim(-0.7, 2.4)
    ax.set_ylim(19.6, 11.5)
    ax.set_xlabel("log₁₀ (period in days)")
    ax.set_ylabel("mean V magnitude (OGLE-IV)")
    ax.legend(loc="upper left", fontsize=7.5, markerscale=1.0)
    ax.set_title(f"b  OGLE today: {len(smc['id']):,} SMC Cepheids")
    # c: periods then and now
    ax = axes[2]
    frac = t["dP_over_P"]
    ali = {a["H"]: a["alias"] for a in lv["aliases"] if a["alias"]}
    for j in range(len(t["H"])):
        h = int(t["H"][j])
        col = ps.ORANGE if h in ali else ps.BLUE
        mk = "s" if h in ali else "o"
        ax.scatter(t["logP"][j], frac[j], s=30, color=col, marker=mk, edgecolors=ps.SURFACE, linewidths=0.6, zorder=3)
        if h in ali or abs(frac[j]) > 5e-3:
            lab = f"HV {h}" + (f"\n({ali[h]})" if h in ali else "")
            ax.annotate(lab, (t["logP"][j], frac[j]), xytext=(6, -3), textcoords="offset points", fontsize=6.5,
                        color=ps.INK)
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.set_xlabel("log₁₀ (period in days)")
    ax.set_ylabel("(P_OGLE − P_1912) / P_1912")
    ax.set_title("c  Her periods, 110 years on")
    ax.set_ylim(-0.03, 0.1)
    ax.scatter([], [], s=30, color=ps.BLUE, marker="o", label="same period, or a small change")
    ax.scatter([], [], s=30, color=ps.ORANGE, marker="s", label="an alias of her observing calendar")
    ax.legend(loc="lower right", fontsize=7)
    fig.suptitle("Henrietta Leavitt's discovery, 1912, and the same stars in OGLE's survey", x=0.06, ha="left", y=1.01)
    ps.caption(fig, "Panel a: Table I of Harvard College Observatory Circular 173 (Leavitt & Pickering 1912), transcribed by "
               f"us from the ADS scan. Her slope is shallower than the same stars' in V today ({lv['fit_v']['b']:.2f}): her "
               "photographic scale compresses faint magnitudes (plate non-linearity, and crowding that brightens the faint "
               "short-period stars; Breuval et al. 2025 find the same −2.03). "
               "Panel b: all 25 stars are fundamental-mode Cepheids in OGLE (identified through their Harvard numbers); "
               "the two longest-period ones are saturated in OGLE's V images. Panel c: three periods differ by one cycle "
               "per year or per lunar month, the classic aliases of plates taken only at night, in season and away from the "
               "full Moon.", y=-0.04)
    return ps.save(fig, _out("02_leavitt_1912"))


def make_all(res) -> list:
    import interactive
    paths = []
    for f in (fig_light_curves, fig_leavitt, fig_leavitt_laws, fig_wesenheit, fig_smc_lmc, fig_maps,
              fig_gaia_calibration, fig_lmc_distance, fig_petersen):
        paths += f(res)
    paths.append(interactive.write(res))
    return paths
