"""Every figure of L4, from the results and the cached data (PNG + SVG)."""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import plotstyle as ps  # noqa: E402
from config import ELEMENTS, FIGURES, REFERENCE  # noqa: E402

ps.use_style()
C_PROMPT = ps.ORANGE       # core-collapse supernovae / massive stars
C_DELAYED = ps.BLUE        # Type Ia supernovae / AGB stars
C_GSE = ps.MAGENTA         # Gaia-Sausage-Enceladus
C_DISC = ps.BLUE
GSE_MAIN = "GALAH orbits (Feuillet+21)"


def _out(name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    return FIGURES / name


def _density(ax, x, y, xr, yr, bins=(160, 140), cmap=ps.CMAP_CONTEXT, log=True):
    s = np.isfinite(x) & np.isfinite(y)
    H, xe, ye = np.histogram2d(x[s], y[s], bins=bins, range=[xr, yr])
    H = np.where(H > 0, np.log10(H) if log else H, np.nan)
    ax.pcolormesh(xe, ye, H.T, cmap=cmap, rasterized=True, shading="flat")


def _spread(ys, gap):
    """Nudge label positions apart so none are closer than ``gap``."""
    ys = np.asarray(ys, float)
    order = np.argsort(ys)
    out = ys.copy()
    for i in range(1, len(order)):
        a, b = order[i - 1], order[i]
        out[b] = max(out[b], out[a] + gap)
    return out


def _two_seq_lines(fit, x):
    hi = fit["plateau"] + fit["slope"] * np.clip(x - fit["knee"], 0, None)
    c0, c1, c2 = fit["low_c"]
    return hi, c0 + c1 * x + c2 * x ** 2


# ---------------------------------------------------------------------------


def _column_density(ax, x, y, xr, yr, bins=(110, 110), cmap=ps.CMAP_CONTEXT):
    """Density normalised within each [Fe/H] column, so both disc sequences
    show at every metallicity however many stars each column holds."""
    s = np.isfinite(x) & np.isfinite(y)
    H, xe, ye = np.histogram2d(x[s], y[s], bins=bins, range=[xr, yr])
    col = H.max(axis=1, keepdims=True)
    Hn = np.where((H > 0) & (col >= 30), H / np.maximum(col, 1), np.nan)
    ax.pcolormesh(xe, ye, Hn.T, cmap=cmap, rasterized=True, shading="flat", vmin=0, vmax=1)


def fig_alpha_knee(res, art, ctx):
    d = ctx.d
    kr = res["knee"]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.4), gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.62))
    ax = axes[0]
    m = ctx.giants & ctx.disc
    _column_density(ax, d["fe_h"][m], d["mg_fe"][m], (-2.2, 0.55), (-0.2, 0.6))
    x = np.linspace(-1.6, 0.45, 300)
    fit = art["knee"]["disc_giants"]["fit"]
    hi, lo = _two_seq_lines(fit, x)
    w = 1 / (1 + np.exp((x - fit["w_x0"]) / fit["w_width"]))
    ax.plot(x[w > 0.3], hi[w > 0.3], color=C_DISC, lw=2.4, label="thick (high-α) disc, fitted (GALAH giants)")
    ax.plot(x[w < 0.9], lo[w < 0.9], color=C_DISC, lw=1.4, ls="--", label="thin (low-α) disc, fitted")
    if isinstance(kr.get("apogee_disc_knee"), dict):
        a = kr["apogee_disc_knee"]
        xa = np.linspace(-1.6, 0.0, 200)
        ax.plot(xa, a["plateau"] + a["slope"] * np.clip(xa - a["knee"], 0, None), color=ps.AQUA, lw=1.6,
                ls="-.", label=f"thick disc in APOGEE giants (knee {a['knee']:+.2f})")
    g = ctx.gse[GSE_MAIN] & ctx.giants
    ax.scatter(d["fe_h"][g], d["mg_fe"][g], s=6, color=C_GSE, alpha=0.6, lw=0,
               label=f"Gaia-Sausage-Enceladus giants ({g.sum():,})", rasterized=True)
    gf = art["knee"][f"gse[{GSE_MAIN}][giants]"]
    xg = np.linspace(-2.2, np.percentile(d["fe_h"][g], 98), 200)
    ax.plot(xg, gf["plateau"] + gf["slope"] * np.clip(xg - gf["knee"], 0, None), color=C_GSE, lw=2.4,
            label="GSE, fitted")
    for kx, ky, col, lab, dx in ((fit["knee"], fit["plateau"], C_DISC, "disc knee", 0.25),
                                 (gf["knee"], gf["plateau"], C_GSE, "GSE knee", -0.25)):
        ax.plot([kx], [ky], marker="o", ms=9, mfc="white", mec=col, mew=2, zorder=5)
        ps.annotate(ax, f"{lab}\n[Fe/H] = {kx:+.2f}", (kx, ky), (kx + dx, ky + 0.17), color=col, ha="center")
    ax.set_xlim(-2.2, 0.55)
    ax.set_ylim(-0.2, 0.6)
    ax.set_xlabel("[Fe/H]  (iron, log relative to the Sun)")
    ax.set_ylabel("[Mg/Fe]  (magnesium per iron)")
    ax.set_title("The α-knee: where Type Ia supernovae start adding iron")
    ax.legend(loc="lower left", fontsize=7.3)

    ax = axes[1]
    rows = [
        ("Disc · GALAH giants · 2-sequence fit", kr["disc_giants"], C_DISC),
        ("Disc · GALAH dwarfs · 2-sequence fit", kr["disc_dwarfs"], C_DISC),
        ("Disc · GALAH giants · ridge-line fit", kr["disc_giants"]["ridge_method"], C_DISC),
        ("Disc · GALAH dwarfs · ridge-line fit", kr["disc_dwarfs"]["ridge_method"], C_DISC),
    ]
    if isinstance(kr.get("apogee_disc_knee"), dict):
        rows.append(("Disc · APOGEE giants · 2-sequence fit", kr["apogee_disc_knee"], ps.AQUA))
    for name in ctx.gse:
        rows.append((f"GSE · {name}", kr[f"gse[{name}][giants]"], C_GSE))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, f, col) in zip(y, rows):
        k = f["knee"]
        ax.errorbar(k, yi, xerr=[[max(k - f["knee_lo"], 0)], [max(f["knee_hi"] - k, 0)]], fmt="o", color=col,
                    ms=6, capsize=3)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=7.5)
    lit = REFERENCE["gse_knee"]
    ax.axvspan(*lit["literature_range"], color=C_GSE, alpha=0.08, lw=0)
    ax.axvline(lit["mackereth2019"], color=C_GSE, lw=1, ls=":")
    ax.text(lit["mackereth2019"] + 0.02, len(rows) - 0.55, "GSE, Mackereth+19 (APOGEE)", fontsize=6.5,
            color=C_GSE, va="top")
    ax.text(-2.0, len(rows) - 0.55, "other published\nGSE knees", fontsize=6.5, color=C_GSE, ha="center",
            va="top")
    ax.set_xlim(-2.6, -0.3)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("[Fe/H] of the knee")
    ax.set_title("Every knee we measured")
    ax.grid(axis="y", visible=False)
    ps.caption(fig, "Left: GALAH DR4 giants on disc-like orbits, density normalised in each [Fe/H] column. "
                    "Right: 16–84% ranges are statistical only (bootstrap; profile likelihood for the disc fits) — "
                    "the spread between methods and surveys is the real uncertainty of the disc's knee. "
                    "Whatever the method, Gaia-Sausage-Enceladus turns down at lower iron: it formed stars more slowly.")
    return ps.save(fig, _out("01_alpha_knee"))


def fig_fingerprint(res):
    fp = res["ncapture"]["fingerprint"]
    els = [e for e in ELEMENTS if e in fp and fp[e]["offset"] is not None and fp[e]["n"] >= 30]
    fig, ax = plt.subplots(figsize=(12, 4.2))
    x = np.arange(len(els))
    off = np.array([fp[e]["offset"] for e in els])
    err = np.array([fp[e]["err"] for e in els])
    cols = [C_GSE if o > 0 else ps.VIOLET for o in off]
    ax.bar(x, off, color=cols, width=0.7, alpha=0.85)
    ax.errorbar(x, off, yerr=err, fmt="none", color=ps.INK, lw=1, capsize=2)
    ax.axhline(0, color=ps.AXIS, lw=1)
    ax.set_xticks(x, els)
    ax.set_ylabel("GSE − disc at the same [Fe/H] (dex)")
    ax.set_title("The chemical fingerprint of another galaxy: Gaia-Sausage-Enceladus versus the Milky Way's disc")
    for xi, e in zip(x, els):
        ax.text(xi, -0.62, f"{fp[e]['n']}", ha="center", fontsize=6, color=ps.MUTED, rotation=90)
    ax.set_ylim(-0.65, max(0.35, off.max() + 0.08))
    low = ", ".join(e for e in els if fp[e]["offset"] < -0.08 and fp[e]["offset"] + 3 * fp[e]["err"] < 0)
    high = ", ".join(e for e in els if fp[e]["offset"] > 0.08 and fp[e]["offset"] - 3 * fp[e]["err"] > 0)
    same = ", ".join(e for e in els if abs(fp[e]["offset"]) < 0.02)
    ps.caption(fig, "GALAH DR4 giants, −1.5 < [Fe/H] < −0.7, medians compared bin by bin (bootstrap errors); "
                    f"numbers = GSE stars measured; Eu from Kane et al. (2026)'s VAC. Clearly low in GSE: {low}. "
                    f"Clearly high: {high}. Indistinguishable: {same}.")
    return ps.save(fig, _out("02_gse_fingerprint"))


def fig_two_process(res):
    pr = res["process"]
    g = pr["giants_ours"]
    dg = pr["dwarfs_griffith_recipe"]
    dw = pr["dwarfs_ours"]
    gr = REFERENCE["griffith2022_fcc"]
    order = [e for e in ELEMENTS if e not in ("Mg",)]
    keys = [("eu" if e == "Eu" else e.lower()) for e in order]
    fig, ax = plt.subplots(figsize=(13, 4.8))
    x = np.arange(len(order))
    ax.axhspan(1.0, 1.4, color=ps.GRID, alpha=0.6, lw=0)
    ax.text(len(order) - 0.6, 1.33, "> 100%: unphysical (no Type Ia share)", fontsize=6.5, color=ps.MUTED,
            ha="right", va="top")

    def series(src, dx, **kw):
        xs, ys, es = [], [], []
        for xi, k in zip(x, keys):
            r = src.get(k)
            if isinstance(r, dict) and r.get("f_cc") is not None:
                xs.append(xi + dx), ys.append(r["f_cc"]), es.append(r.get("f_cc_err") or 0)
        ax.errorbar(xs, ys, yerr=es, fmt="o", capsize=2, **kw)

    series(g, -0.18, color=ps.BLUE, ms=6, label="this work: GALAH DR4 giants (measured plateau)")
    series(dg, 0.0, color=ps.AQUA, ms=5, label="this work: GALAH DR4 dwarfs, Griffith et al.'s recipe")
    series(dw, 0.18, color=ps.MUTED, ms=4, mfc="white", label="DR4 dwarfs with their own (0.11 dex higher) plateau")
    gx = [xi for xi, e in zip(x, order) if e in gr]
    ax.scatter(gx, [gr[e] for e in order if e in gr], marker="D", s=30, color=ps.ORANGE, zorder=5,
               label="Griffith et al. 2022 (GALAH DR3 dwarfs)")
    # Europium: the VAC, DR3's own, and DR4's flawed values.
    xe = order.index("Eu")
    for k, mk, lab in (("eu_dr3", "s", "Eu: DR3 measurements"), ("eu_dr4", "x", "Eu: DR4's own (unreliable)")):
        r = g.get(k)
        if r and r.get("f_cc") is not None:
            ax.scatter([xe + 0.36], [r["f_cc"]], marker=mk, s=30, color=ps.INK_SECONDARY, zorder=6, label=lab)
    ax.set_xticks(x, order)
    ax.set_ylim(-0.05, 1.4)
    ax.set_xlim(-0.7, len(order) - 0.3)
    ax.set_ylabel("prompt share (core-collapse timescale)")
    ax.set_title("How much of each element arrived promptly? The two-process decomposition of GALAH's elements")
    ax.legend(loc="upper left", fontsize=7, ncol=2)
    ps.caption(fig, "Prompt share f_cc on the low-α disc sequence at [Mg/H] = 0 (≈ the Sun), Weinberg et al. "
                    "(2019)'s two-process model. Error bars are bootstrap (statistical) only and mostly smaller than "
                    "the dots; the spread between giants and dwarfs is the real uncertainty. Missing points fail the "
                    "detection test (a weak line detected in < 70% of stars biases the medians).")
    return ps.save(fig, _out("03_prompt_fractions"))


def fig_periodic(res):
    import periodic as P

    pr = res["process"]["giants_ours"]
    dg = res["process"]["dwarfs_griffith_recipe"]
    vals, uncertain = {"Mg": 1.0}, set()
    for e in ELEMENTS:
        k = "eu" if e == "Eu" else e.lower()
        r = pr.get(k)
        if isinstance(r, dict) and r.get("f_cc") is not None:
            vals[e] = r["f_cc"]
            rd = dg.get(k)
            if isinstance(rd, dict) and rd.get("f_cc") is not None and abs(rd["f_cc"] - r["f_cc"]) > 0.15:
                uncertain.add(e)
    s = res["ncapture"]["s_fraction"]
    compare = dict(REFERENCE["griffith2022_fcc"])
    compare.pop("O", None)
    compare.update({"Eu": 1 - s["Eu"]})
    for e in ("Ba",):
        compare[e] = 1 - s[e]
    fig, ax = plt.subplots(figsize=(14, 8.2))
    P.draw(ax, vals, compare=compare, uncertain=uncertain, bigbang={"H": 1.0, "He": 0.91},
           prompt_color=C_PROMPT, delayed_color=C_DELAYED, bigbang_color=ps.VIOLET)
    fig.suptitle("Where your atoms were made — measured in 160,000 Milky Way giants", x=0.06, ha="left", y=0.97)
    from matplotlib.patches import Rectangle

    # The explanation in the empty top row, the key in the gap between Be/Mg and B/Al.
    ax.text(1.25, -0.25, "Colour = share of each element that arrived promptly (bottom) or late (top) in "
                         "solar-metallicity disc stars,\nmeasured with the two-process model in GALAH DR4. "
                         "Hatched: dwarfs and giants disagree by > 15 points. Grey: not measured here.",
            fontsize=8, color=ps.INK_SECONDARY, va="top")
    key = [(C_PROMPT, "prompt: core-collapse supernovae of massive stars (< ~40 Myr)"),
           (C_DELAYED, "delayed: Type Ia supernovae, AGB stars (0.1–10 Gyr)"),
           (ps.VIOLET, "the Big Bang (H; 91% of He by mass)")]
    for i, (col, lab) in enumerate(key):
        y = -1.25 - 0.45 * i
        ax.add_patch(Rectangle((2.2, y - 0.17), 0.42, 0.34, color=col))
        ax.text(2.75, y, lab, fontsize=7.8, va="center")
    y = -1.25 - 0.45 * 3
    ax.plot([2.22, 2.6], [y, y], color=ps.INK, lw=1.6)
    ax.text(2.75, y, "published: Griffith+22 prompt share; r-process share of Ba, Eu (Bisterzo+14)",
            fontsize=7.8, va="center")
    ps.caption(fig, "Mg defines the prompt scale (100%). Europium arrived as promptly as magnesium — a puzzle if "
                    "neutron-star mergers, which need time to spiral in, made it. 'Delayed' for Ba, Y, Zr, La is "
                    "the slow neutron-capture process in AGB stars.")
    return ps.save(fig, _out("04_periodic_table_of_origins"))


def fig_ncapture(res, ctx, jina):
    d, eu = ctx.d, ctx.eu
    nr = res["ncapture"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    m = ctx.giants & ctx.disc
    specs = [("ba_fe", d["ba_fe"], "[Ba/Fe]  (barium: mostly s-process)", (-1.0, 1.2), "ba"),
             ("eu_fe", eu["eu_best"], "[Eu/Fe]  (europium: mostly r-process)", (-0.6, 1.4), "eu"),
             ("ba_eu", d["ba_fe"] - eu["eu_best"], "[Ba/Eu]", (-1.4, 1.0), None)]
    for ax, (q, y, lab, yr, jel) in zip(axes, specs):
        if jel is not None:
            det = ~jina[f"{jel}_upper"] & np.isfinite(jina[f"{jel}_fe"])
            ax.scatter(jina["feh"][det], jina[f"{jel}_fe"][det], s=4, color=ps.MUTED, alpha=0.5, lw=0,
                       label="JINAbase metal-poor stars", rasterized=True)
            ul = jina[f"{jel}_upper"]
            ax.scatter(jina["feh"][ul], jina[f"{jel}_fe"][ul], s=8, marker="v", color=ps.MUTED, alpha=0.35,
                       lw=0, label="  … upper limits", rasterized=True)
        _density(ax, d["fe_h"][m], y[m], (-2.0, 0.5), yr, bins=(100, 90), cmap=ps.CMAP_BLUE)
        t = nr["trends"][q]
        ax.plot([r["x"] for r in t["disc"]], [r["median"] for r in t["disc"]], color=ps.INK, lw=1.8,
                label="disc giants (median)")
        gs = t["gse"]
        ax.errorbar([r["x"] for r in gs], [r["median"] for r in gs],
                    yerr=[[r["median"] - r["p16"] for r in gs], [r["p84"] - r["median"] for r in gs]],
                    fmt="o", color=C_GSE, ms=6, capsize=2, label="GSE giants (median, 16–84%)")
        ax.axhline(0, color=ps.AXIS, lw=0.8)
        ax.set_xlim(-4.0, 0.5)
        ax.set_ylim(*yr)
        ax.set_xlabel("[Fe/H]")
        ax.set_ylabel(lab)
    both = (~jina["ba_upper"] & ~jina["eu_upper"] & np.isfinite(jina["ba_fe"]) & np.isfinite(jina["eu_fe"]))
    axes[2].scatter(jina["feh"][both], (jina["ba_fe"] - jina["eu_fe"])[both], s=5, color=ps.MUTED, alpha=0.6,
                    lw=0, rasterized=True)
    axes[2].axhline(nr["ba_eu_pure_r"], color=C_PROMPT, lw=1.5, ls="--")
    axes[2].text(-3.9, nr["ba_eu_pure_r"] - 0.05, f"pure r-process (solar pattern): {nr['ba_eu_pure_r']:+.2f}",
                 color=C_PROMPT, fontsize=7.5, va="top")
    axes[2].text(-3.9, 0.04, "solar mix", color=ps.MUTED, fontsize=7.5)
    axes[0].set_title("Barium", loc="left")
    axes[1].set_title("Europium", loc="left")
    axes[2].set_title("Barium per europium: r-process first, s-process later", loc="left")
    axes[0].legend(loc="upper left", fontsize=6.5)
    off = nr["gse_minus_disc"]
    axes[1].text(-3.9, -0.45, f"GSE − disc at the same [Fe/H]:\n[Eu/Fe] {off['eu_fe']['offset']:+.2f} ± "
                              f"{off['eu_fe']['err']:.2f},  [Eu/Mg] {off['eu_mg']['offset']:+.2f} ± "
                              f"{off['eu_mg']['err']:.2f}", fontsize=7.5, color=C_GSE, va="bottom")
    ps.caption(fig, "Blue density: GALAH DR4 disc giants (Eu from Kane et al. 2026's VAC; DR3 values where it has none). "
                    "Grey: JINAbase (Abohalima & Frebel 2018; 2021 update), Milky Way field stars, CEMP-s stars removed.")
    return ps.save(fig, _out("05_neutron_capture"))


def fig_scatter(res):
    nr = res["ncapture"]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    cols = {"mg_fe": ps.ORANGE, "ba_fe": ps.AQUA, "eu_fe": ps.MAGENTA}
    names = {"mg_fe": "Mg (core-collapse, frequent)", "ba_fe": "Ba", "eu_fe": "Eu (r-process, rare events?)"}
    for q in ("mg_fe", "ba_fe", "eu_fe"):
        j = nr["scatter_jinabase"][q]
        ax.plot([r["x"] for r in j], [r["robust_sigma"] for r in j], "o--", color=cols[q], mfc="white",
                ms=6, lw=1.2)
        gl = nr["scatter_galah"][q]
        ax.plot([r["x"] for r in gl], [r["intrinsic"] if q != "eu_fe" else r["robust_sigma"] for r in gl],
                "o-", color=cols[q], ms=5, lw=1.8, label=names[q])
    ax.plot([], [], "o--", color=ps.MUTED, mfc="white", label="JINAbase (literature compilation)")
    ax.plot([], [], "o-", color=ps.MUTED, label="GALAH DR4 giants")
    ax.set_xlabel("[Fe/H]")
    ax.set_ylabel("star-to-star spread (dex, robust σ)")
    ax.set_title("Europium is patchy in the oldest stars, magnesium is not")
    ax.set_ylim(0, 0.75)
    ax.legend(fontsize=7.5)
    ps.caption(fig, "Robust σ = 1.48 × median absolute deviation, detections only (JINAbase also lists upper limits, "
                    "which are not used: true Eu and Ba spreads at the lowest [Fe/H] are larger). GALAH Mg and Ba: "
                    "measurement errors removed. JINAbase mixes studies, and r-process-rich stars were targeted on purpose.")
    return ps.save(fig, _out("06_scatter"))


def fig_eu_check(res, ctx):
    """DR4's own europium against DR3 for the same stars, and the trends of
    the three europium sources."""
    import ncapture as nc

    d, eu = ctx.d, ctx.eu
    nr = res["ncapture"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    s = np.isfinite(d["eu_fe"]) & np.isfinite(eu["eu_dr3"])
    ax.hexbin(eu["eu_dr3"][s], d["eu_fe"][s], gridsize=70, extent=(-0.5, 1.2, -0.5, 1.2), bins="log",
              cmap=ps.CMAP_BLUE, mincnt=1, rasterized=True)
    ax.plot([-0.5, 1.2], [-0.5, 1.2], color=ps.ORANGE, lw=1)
    st = nr["eu_dr4_vs_dr3"]
    ax.text(-0.45, 1.15, f"{st['n']:,} stars measured by both\nr = {st['r']:.2f}, spread {st['sigma']:.2f} dex",
            va="top", fontsize=8)
    ax.set_xlabel("GALAH DR3 [Eu/Fe]")
    ax.set_ylabel("GALAH DR4's own [Eu/Fe] (same stars)")
    ax.set_title("DR4's pipeline lost europium", loc="left")
    ax = axes[1]
    edges = np.round(np.arange(-1.2, 0.51, 0.1), 2)
    g = ctx.giants
    for y, lab, col, ls in ((eu["eu_dr3"], "DR3 measurements (52k training stars)", ps.ORANGE, "-"),
                            (eu["eu_golden"], "Kane et al. CNN, golden sample (new stars)", ps.BLUE, "-"),
                            (np.where(np.isfinite(eu["eu_dr3"]) | np.isfinite(eu["eu_golden"]), np.nan, d["eu_fe"]),
                             "DR4's own, all other giants", ps.MUTED, "--")):
        rows = nc.binned(d["fe_h"][g], y[g], edges, min_n=30)
        xs = [r["x"] for r in rows]
        ax.plot(xs, [r["median"] for r in rows], ls, color=col, lw=2, label=lab)
        ax.fill_between(xs, [r["p16"] for r in rows], [r["p84"] for r in rows], color=col, alpha=0.12, lw=0)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.set_xlabel("[Fe/H]")
    ax.set_ylabel("[Eu/Fe] (giants: median, 16–84%)")
    ax.set_title("The re-derived europium follows DR3's trend", loc="left")
    ax.legend(fontsize=7.5)
    ps.caption(fig, "DR4's neural-network pipeline lost precision for europium (Buder et al. 2025); Kane et al. "
                    "(2026) re-derived it from the 6645 Å line with a network trained on DR3's measurements, for stars "
                    "DR3 never measured. This project uses the CNN golden sample, then DR3's values.")
    return ps.save(fig, _out("07_europium_check"))


def fig_kilonova(res, art):
    import kilonova as K

    kn, best, vil = art["kilonova"]["kn"], art["kilonova"]["best"], art["kilonova"]["villar"]
    bands = [("W2", "UV (Swift W2)"), ("g", "g"), ("r", "r"), ("i", "i"), ("z", "z"), ("J", "J"),
             ("H", "H"), ("Ks", "Ks")]
    offsets = np.arange(len(bands)) * 1.5 - 3
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax = axes[0]
    cmap = matplotlib.colormaps["plasma"]
    tt = np.logspace(np.log10(0.3), np.log10(30), 120)
    ends, labels, colors = [], [], []
    for i, ((b, lab), off) in enumerate(zip(bands, offsets)):
        col = cmap(0.05 + 0.85 * i / (len(bands) - 1))
        sel = (kn.band == b) & ~kn.upper
        ul = (kn.band == b) & kn.upper
        if not (sel | ul).any():
            continue
        ax.errorbar(kn.t_obs[sel], kn.mag[sel] + off, yerr=kn.err[sel], fmt="o", ms=3, color=col, lw=0.6)
        ax.scatter(kn.t_obs[ul], kn.mag[ul] + off, marker="v", s=12, color=col, alpha=0.5)
        idx = np.full(tt.size, kn.bands.ids.index(kn.svo[sel | ul][0]))
        mod = kn.model_mag(best, t_rest=tt / (1 + K.KN_REDSHIFT), idx=idx) + off
        ax.plot(tt, mod, color=col, lw=1.6)
        ax.plot(tt, kn.model_mag(vil, t_rest=tt / (1 + K.KN_REDSHIFT), idx=idx) + off, color=col, lw=1,
                ls=":")
        ends.append(mod[-1]), labels.append(f"{lab} {off:+.1f}"), colors.append(col)
    for yv, lab, col in zip(_spread(ends, 0.45), labels, colors):
        ax.text(31, yv, lab, fontsize=7, color=col, va="center")
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlim(0.3, 60)
    ax.set_xlabel("days after the merger")
    ax.set_ylabel("AB magnitude (+ offset)")
    ax.set_title("AT2017gfo: 646 measurements, 65 filters")
    ax.plot([], [], color=ps.INK, lw=1.6, label="our best fit")
    ax.plot([], [], color=ps.INK, lw=1, ls=":", label="Villar et al. 2017's parameters, our model")
    ax.legend(loc="lower left", fontsize=7.5)

    ax = axes[1]
    curves = K.component_curves(best)
    tg = K.T_GRID
    col = {"blue": ps.BLUE, "purple": ps.VIOLET, "red": ps.RED}
    total = sum(curves[c]["lum"] for c in K.COMPONENTS)
    p = dict(zip(K.PARAMS, best))
    heat = sum(K.heating_rate(tg, p[f"m_{c}"]) * K.thermalisation(tg, p[f"m_{c}"], p[f"v_{c}"])
               for c in K.COMPONENTS)
    raw = sum(K.heating_rate(tg, p[f"m_{c}"]) for c in K.COMPONENTS)
    ax.plot(tg, raw, color=ps.MUTED, lw=1, ls="--", label="radioactive power of the r-process nuclei")
    ax.plot(tg, heat, color=ps.MUTED, lw=1.4, label="… of which heats the gas (thermalised)")
    for c in K.COMPONENTS:
        ax.plot(tg, curves[c]["lum"], color=col[c], lw=1.4,
                label=f"{c}: {p[f'm_{c}']:.3f} M☉ at {p[f'v_{c}'] / 2.998e5:.2f}c (κ = {K.KN_KAPPA[c]:g})")
    ax.plot(tg, total, color=ps.INK, lw=2.2, label="total light")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 40)
    ax.set_ylim(1e39, 3e42)
    ax.set_xlabel("days after the merger")
    ax.set_ylabel("luminosity (erg/s)")
    ax.set_title("What powers it: decaying r-process nuclei")
    ax.legend(fontsize=6.8, loc="lower left")
    ps.caption(fig, "Data: Villar et al. (2017)'s unified set (rows they used). Model: their three-component "
                    "radioactive kilonova (Metzger 2017), our implementation, every point integrated over its filter. "
                    "Best fit by nested sampling, D = 39.5 Mpc.")
    return ps.save(fig, _out("08_kilonova_lightcurves"))


def fig_kilonova_masses(res, art):
    s = art["kilonova"]["samples"]
    v = REFERENCE["villar2017"]
    kr = res["kilonova"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axes[0]
    ax.scatter(s[:, 3], s[:, 6], s=3, color=ps.BLUE, alpha=0.25, lw=0, rasterized=True)
    ax.errorbar(v["purple"]["m"], v["red"]["m"], fmt="D", color=ps.ORANGE, ms=8, label="Villar et al. 2017")
    ax.set_xlabel("purple ejecta (κ = 3) [M☉]")
    ax.set_ylabel("red, lanthanide-rich ejecta (κ = 10) [M☉]")
    ax.set_title("Two ways to split the slow ejecta", loc="left")
    ax.legend(fontsize=7.5)
    for ax, (k, lab, ref) in zip(axes[1:], (("total", "all ejecta [M☉]", v["blue"]["m"] + v["purple"]["m"] + v["red"]["m"]),
                                            ("purple_red", "lanthanide-bearing ejecta, κ ≥ 3 [M☉]",
                                             v["purple"]["m"] + v["red"]["m"]))):
        x = s[:, 0] + s[:, 3] + s[:, 6] if k == "total" else s[:, 3] + s[:, 6]
        ax.hist(x, bins=60, color=ps.BLUE, alpha=0.8)
        ax.axvline(ref, color=ps.ORANGE, lw=2, label=f"Villar et al.: {ref:.3f}")
        q = np.percentile(x, [16, 50, 84])
        ax.set_title(f"{lab.split('[')[0].strip()}: {q[1]:.3f} (+{q[2] - q[1]:.3f} −{q[1] - q[0]:.3f}) M☉",
                     loc="left", fontsize=9)
        ax.set_xlabel(lab)
        ax.legend(fontsize=7.5)
    modes = kr["modes"]
    txt = ", ".join(f"{n}: {100 * m['weight']:.0f}%" for n, m in modes.items() if n != "all")
    ps.caption(fig, f"Posterior from nested sampling (dynesty). Posterior weight of the two modes — {txt}. "
                    "The split between purple and red is uncertain; their sum is not.")
    return ps.save(fig, _out("09_kilonova_masses"))


def fig_budget(res):
    b = res["budget"]
    import budget as B

    rates = np.logspace(0, 3.5, 200)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axes[0]
    from config import J3_LOG_RHO_STAR, MW_AGE_GYR, MW_STELLAR_MASS

    per_mass = 1e-9 / 10 ** J3_LOG_RHO_STAR       # mergers per yr per Msun of stars, for R0 = 1 Gpc^-3 yr^-1

    cols = {"all": ps.BLUE, "purple+red": ps.VIOLET, "red only": ps.RED}
    for reading, (lo, mid, hi) in b["per_event_eu"].items():
        supply = rates * per_mass * MW_STELLAR_MASS[0] * b["history_factor"] * MW_AGE_GYR * 1e9 / b["galaxy"]["m_eu"]
        ax.fill_between(rates, supply * lo, supply * hi, color=cols[reading], alpha=0.15, lw=0)
        ax.plot(rates, supply * mid, color=cols[reading], lw=2,
                label=f"Eu-making ejecta = {reading}: {mid:.1e} M☉ of Eu per merger")
    ax.axhline(1, color=ps.INK, lw=1)
    ax.text(1.1, 1.08, "enough for all the Galaxy's europium", fontsize=7.5)
    lo, hi = REFERENCE["gwtc4_bns_rate"]
    ax.axvspan(lo, hi, color=ps.AQUA, alpha=0.15, lw=0)
    ax.text(np.sqrt(lo * hi), 0.012, "GWTC-4.0\nmerger rate\n(90%)", ha="center", fontsize=7.5, color=ps.GREEN)
    ax.axvline(REFERENCE["bns_rate_2017"], color=ps.MUTED, ls=":", lw=1)
    ax.text(REFERENCE["bns_rate_2017"] * 1.07, 0.012, "2017 rate\n(after\nGW170817)", fontsize=7, color=ps.MUTED)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.01, 100)
    ax.set_xlabel("neutron-star merger rate today (Gpc⁻³ yr⁻¹)")
    ax.set_ylabel("europium supplied over 12 Gyr ÷ europium in the Galaxy")
    ax.set_title("Could GW170817-like mergers have made all our europium?")
    ax.legend(fontsize=7, loc="upper left")
    ax = axes[1]
    hc = b["history_curve"]
    ax.plot(hc["t_gyr"], hc["rate_rel"], color=ps.BLUE, lw=2)
    ax.axhline(b["history_factor"], color=ps.ORANGE, ls="--", lw=1)
    ax.text(0.5, b["history_factor"] * 1.05, f"mean over the last 12 Gyr: {b['history_factor']:.2f}× today",
            color=ps.ORANGE, fontsize=7.5)
    ax.set_xlabel("cosmic time (Gyr)")
    ax.set_ylabel("merger rate ÷ today's")
    ax.set_title("Mergers were commoner in the past", loc="left")
    del B
    g = b["galaxy"]
    ps.caption(fig, f"Galaxy: (5 × 10¹⁰ M☉ stars + 8 × 10⁹ M☉ gas) × solar Eu mass fraction {g['x_eu_sun']:.2e} × "
                    f"mean disc-giant Eu level {b['mean_eu_h_linear_disc_giants']:.2f} × r-share 0.94 = {g['m_eu']:.0f} ± "
                    f"{g['m_eu_err']:.0f} M☉ of Eu. Per merger: ejecta × Eu share of solar r-process matter "
                    f"(Kr–U, {100 * b['shares']['eu_share_kr_up']:.2f}%). Rate scaled by stellar mass (J3's local density); "
                    "history: Madau & Dickinson (2014) ⊗ 1/t delays from 20 Myr.")
    return ps.save(fig, _out("10_europium_budget"))


def fig_overview(res, art, ctx):
    """One-page summary: knee, periodic table, kilonova, budget."""
    import periodic as P

    d = ctx.d
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1.05], hspace=0.32, wspace=0.28)
    ax = fig.add_subplot(gs[0, 0])
    m = ctx.giants & ctx.disc
    _density(ax, d["fe_h"][m], d["mg_fe"][m], (-2.2, 0.5), (-0.2, 0.6), bins=(110, 90))
    g = ctx.gse[GSE_MAIN] & ctx.giants
    ax.scatter(d["fe_h"][g], d["mg_fe"][g], s=3, color=C_GSE, alpha=0.5, lw=0, rasterized=True)
    kr = res["knee"]
    for f, col, xmax in ((art["knee"]["disc_giants"]["fit"], C_DISC, 0.0),
                         (art["knee"][f"gse[{GSE_MAIN}][giants]"], C_GSE, -0.6)):
        x = np.linspace(-2.2, xmax, 200)
        ax.plot(x, f["plateau"] + f["slope"] * np.clip(x - f["knee"], 0, None), color=col, lw=2)
        ax.plot(f["knee"], f["plateau"], "o", mfc="white", mec=col, mew=2, ms=8)
    ax.set_title(f"1. The α-knee: disc {kr['disc_giants']['knee']:+.2f}, GSE "
                 f"{kr[f'gse[{GSE_MAIN}][giants]']['knee']:+.2f}", loc="left", fontsize=10)
    ax.set_xlabel("[Fe/H]")
    ax.set_ylabel("[Mg/Fe]")
    ax = fig.add_subplot(gs[0, 1:])
    pr = res["process"]["giants_ours"]
    vals = {"Mg": 1.0}
    for e in ELEMENTS:
        r = pr.get("eu" if e == "Eu" else e.lower())
        if isinstance(r, dict) and r.get("f_cc") is not None:
            vals[e] = r["f_cc"]
    P.draw(ax, vals, bigbang={"H": 1.0, "He": 0.91}, prompt_color=C_PROMPT, delayed_color=C_DELAYED,
           bigbang_color=ps.VIOLET)
    ax.set_title("2. Prompt (orange) vs delayed (blue) share of each element, from GALAH DR4", loc="left",
                 fontsize=10)
    ax = fig.add_subplot(gs[1, 0:2])
    kn, best = art["kilonova"]["kn"], art["kilonova"]["best"]
    import kilonova as K

    tt = np.logspace(np.log10(0.3), np.log10(30), 100)
    cmap = matplotlib.colormaps["plasma"]
    ends, labs, cols = [], [], []
    for i, (b, off) in enumerate((("g", -1.5), ("r", 0.0), ("i", 1.5), ("z", 3.0), ("J", 4.5), ("Ks", 6.0))):
        col = cmap(0.05 + 0.85 * i / 5)
        sel = (kn.band == b) & ~kn.upper
        ax.plot(kn.t_obs[sel], kn.mag[sel] + off, "o", ms=2.5, color=col)
        idx = np.full(tt.size, kn.bands.ids.index(kn.svo[sel][0]))
        mod = kn.model_mag(best, t_rest=tt / (1 + K.KN_REDSHIFT), idx=idx) + off
        ax.plot(tt, mod, color=col, lw=1.5)
        ends.append(mod[-1]), labs.append(b), cols.append(col)
    for yv, b, col in zip(_spread(ends, 0.5), labs, cols):
        ax.text(32, yv, b, color=col, fontsize=8, va="center")
    ax.set_xscale("log")
    ax.invert_yaxis()
    kq = res["kilonova"]["modes"]["all"]["quantiles"]
    ax.set_title(f"3. GW170817's kilonova: {kq['m_total'][1]:.3f} M☉ of r-process ejecta, "
                 f"{kq['m_purple_red'][1]:.3f} M☉ of it lanthanide-bearing", loc="left", fontsize=10)
    ax.set_xlabel("days after the merger")
    ax.set_ylabel("AB mag + offset")
    ax = fig.add_subplot(gs[1, 2])
    b = res["budget"]
    labels = ["red only", "purple+red", "all"]
    mid = [b["budgets"][k]["supply_over_demand"][1] for k in labels]
    lo = [b["budgets"][k]["supply_over_demand"][0] for k in labels]
    hi = [b["budgets"][k]["supply_over_demand"][2] for k in labels]
    y = np.arange(3)
    ax.errorbar(mid, y, xerr=[np.array(mid) - lo, np.array(hi) - mid], fmt="o", color=ps.BLUE, capsize=3, ms=7)
    ax.axvline(1, color=ps.INK, lw=1)
    ax.set_xscale("log")
    ax.set_yticks(y, ["only red ejecta\nmake Eu", "purple + red", "all ejecta"])
    ax.set_xlabel("Eu made by GW170817-like mergers ÷ Eu in the Galaxy")
    ax.set_title("4. Enough europium? (bars: GWTC-4.0 rate range)", loc="left", fontsize=10)
    ax.grid(axis="y", visible=False)
    fig.suptitle("L4 — Where your atoms were made: 530,000 GALAH stars and one kilonova", x=0.06, ha="left")
    return ps.save(fig, _out("00_overview"))


def make_all(res, art):
    ctx = art["ctx"]
    written = []
    written += fig_alpha_knee(res, art, ctx)
    written += fig_fingerprint(res)
    written += fig_two_process(res)
    written += fig_periodic(res)
    written += fig_ncapture(res, ctx, art["jinabase"])
    written += fig_scatter(res)
    written += fig_eu_check(res, ctx)
    written += fig_kilonova(res, art)
    written += fig_kilonova_masses(res, art)
    written += fig_budget(res)
    written += fig_overview(res, art, ctx)
    return written
