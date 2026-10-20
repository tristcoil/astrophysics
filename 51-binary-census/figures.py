"""Every figure of the project, drawn from the analysis results.

Each ``fig_*`` function takes the results dict from ``analysis.run()`` and
returns the paths it wrote (PNG + SVG into ``figures/``).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm

import config as C
import literature as L
import plotstyle as ps
import wide

ps.use_style()

#: Sequential blue ramp (palette steps 100 -> 700), for densities.
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#e8f1fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
REF_STYLE = {"color": ps.INK, "marker": "s", "mfc": "none", "mec": ps.INK, "ms": 5, "lw": 0}


def _out(name):
    return str(C.FIGURE_DIR / name)


def _logp_axis(ax, which="x"):
    """Secondary labels in years/AU-friendly units are added per figure."""
    getattr(ax, f"set_{which}label")("orbital period, log10(P / day)")


# ---------------------------------------------------------------------------
# Gaia's eyesight for close pairs
# ---------------------------------------------------------------------------

def fig_eyesight(res) -> list:
    import samples
    pairs = samples.wide_pairs()
    cur = res["contrast"]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.28})
    ax = axes[0]
    lt = np.log10(np.asarray(pairs["theta_arcsec"]))
    dg = np.asarray(pairs["dg"])
    h = ax.hist2d(lt, dg, bins=[np.arange(-0.7, 1.31, 0.025), np.arange(0, 13.01, 0.125)],
                  cmap=SEQ_BLUE, norm=LogNorm(vmin=1, vmax=400), rasterized=True)
    for key, style, label in (("all", dict(color=ps.INK, lw=1.8), "50% detection: all pairs"),
                              ("colour", dict(color=ps.ORANGE, lw=1.8, ls="--"),
                               "50% detection: both stars with BP-RP colours")):
        c = cur[key]
        ax.plot(c["log_theta0"], c["dg_mid"], **style, label=label)
    ax.plot([0.0, np.log10(2.0)], [4.0, 7.0], "o", color=ps.INK, mfc=ps.SURFACE, ms=7, mew=1.4,
            label="El-Badry et al. 2021: contrast limit")
    ax.set_xlabel("angular separation, log10(θ / arcsec)")
    ax.set_ylabel("brightness difference ΔG (mag)")
    ax.set_xlim(-0.7, 1.3)
    ax.set_ylim(13, 0)
    ax.set_title("a  How close can Gaia split a pair?")
    ax.legend(loc="lower right", fontsize=7.5)
    cb = fig.colorbar(h[3], ax=ax, pad=0.01, fraction=0.05)
    cb.set_label("pairs per cell (R < 0.1, within 200 pc)", color=ps.INK_SECONDARY, fontsize=8)
    cb.outline.set_visible(False)

    ax = axes[1]
    show = [(0.0, 0.5), (3.0, 4.0), (8.0, 9.5)]
    colours = [ps.BLUE, ps.ORANGE, ps.AQUA]
    markers = ["o", "s", "D"]
    for (lo, hi), col, mk in zip(show, colours, markers):
        b = next(b for b in cur["all"]["bins"] if b["dg_lo"] == lo and b["dg_hi"] == hi)
        x, n, model = b["curve"]
        ax.plot(x, n, mk, color=col, ms=3.5, alpha=0.8, mec="none")
        ax.plot(x, model, "-", color=col, lw=1.6, label=f"ΔG {lo:g}–{hi:g}: θ₀ = {10 ** b['log_theta0']:.2f}″")
    ax.set_yscale("log")
    ax.set_ylim(0.7, 3e3)
    ax.set_xlabel("angular separation, log10(θ / arcsec)")
    ax.set_ylabel("pairs per 0.05 dex")
    ax.set_title("b  The shortfall at small separation is the incompleteness")
    ax.legend(loc="lower right", fontsize=7.5, title="data (points) and fit (lines)", title_fontsize=7.5)
    ps.caption(fig, "Pairs from El-Badry, Rix & Heintz (2021) within 200 pc with chance-alignment ratio R < 0.1. "
               "In each ΔG range the number of pairs per log θ is fitted as a power law times "
               "S(θ) = 1/(1 + (θ/θ₀)^−β) (El-Badry et al. 2019, Eq. D5).", y=-0.04)
    return ps.save(fig, _out("08_gaia_eyesight"))


# ---------------------------------------------------------------------------
# Twins among wide binaries
# ---------------------------------------------------------------------------

def fig_twins(res) -> list:
    tw = res["twins"]
    fig = plt.figure(figsize=(13.0, 7.4))
    gs = fig.add_gridspec(2, 5, height_ratios=[1.15, 1.0], hspace=0.42, wspace=0.12)
    ax = fig.add_subplot(gs[0, :3])
    styles = [(ps.BLUE, "o", "-"), (ps.ORANGE, "s", "-"), (ps.AQUA, "D", "-")]
    for ((slo, shi), (h, e)), (col, mk, ls) in zip(tw["hists"].items(), styles):
        c = 0.5 * (e[1:] + e[:-1])
        dens = h / h.sum() / np.diff(e)
        ax.step(e[:-1], dens, where="post", color=col, lw=1.6)
        ax.plot(c, dens, mk, color=col, ms=4, label=f"{slo:,}–{shi:,} au ({h.sum():,} pairs)")
    ax.axvspan(0.95, 1.0, color=ps.WASH, zorder=0)
    ax.text(0.975, ax.get_ylim()[1] * 0.97 if ax.get_ylim()[1] else 1, "twins\nq > 0.95", ha="center", va="top",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlabel("mass ratio q = m₂ / m₁ (from both stars' absolute magnitudes)")
    ax.set_ylabel("fraction of pairs per unit q")
    ax.set_xlim(0.3, 1.0)
    ax.set_title("a  Sun-like primaries (0.8–1.2 M☉): the spike of twins fades with separation")
    ax.set_ylim(0, 3.3)
    ax.legend(loc="upper center", fontsize=7.5, ncol=3, bbox_to_anchor=(0.5, 1.0))

    ax = fig.add_subplot(gs[0, 3:])
    ax.axis("off")
    ax.text(0.0, 0.98,
            "What the panels below show\n\n"
            "Each mass and separation bin is fitted with El-Badry et al.'s (2019)\n"
            "model: a broken power law in q plus a step of extra 'twins' above\n"
            "q ≈ 0.95. F_twin is that excess as a fraction of all pairs with\n"
            "q > 0.3. Every pair's chance of being resolved by Gaia is folded\n"
            "into the likelihood, and only pairs wider than 3.5″ are used, where\n"
            "both stars have colours whatever their brightness difference.\n\n"
            "Filled circles: this work (median, 68% interval).\n"
            "Open squares: El-Badry et al. 2019 (Gaia DR2), their Table G1\n"
            "(median, 68% interval taken as half their 2σ bars).",
            va="top", fontsize=8, color=ps.INK_SECONDARY, linespacing=1.35)

    for i, (mlo, mhi) in enumerate(wide.TWIN_MASS_BINS):
        ax = fig.add_subplot(gs[1, i])
        xs, ys, lo, hi = [], [], [], []
        for j, (slo, shi) in enumerate(wide.TWIN_SEP_BINS):
            f = tw["fits"].get((i, j), {})
            if "f_twin" not in f:
                continue
            xs.append(np.sqrt(slo * shi))
            ys.append(f["f_twin"])
            lo.append(f["f_twin"] - f["f_twin_16"])
            hi.append(f["f_twin_84"] - f["f_twin"])
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o", color=ps.BLUE, ms=5, lw=1.2, capsize=0, label="this work")
        eb = L.ELBADRY19_FTWIN[:, i]
        xe = [np.sqrt(a * b) * 1.12 for a, b in wide.TWIN_SEP_BINS]
        ax.errorbar(xe, eb[:, 0], yerr=[eb[:, 2] / 2, eb[:, 1] / 2], fmt="s", mfc="none", mec=ps.INK, color=ps.INK,
                    ms=5, lw=1.0, capsize=0, label="El-Badry+2019")
        ax.axhline(0, color=ps.AXIS, lw=0.8)
        ax.set_xscale("log")
        ax.set_ylim(-0.04, 0.24)
        ax.set_title(f"{mlo}–{mhi} M☉", fontsize=9, loc="center")
        ax.set_xlabel("separation (au)")
        if i == 0:
            ax.set_ylabel("F_twin (excess twins)")
            ax.legend(loc="upper right", fontsize=7)
        else:
            ax.set_yticklabels([])
    fig.text(0.5, 0.465, "b  The twin excess by primary mass and separation", ha="center", fontsize=11,
             fontweight="bold", color=ps.INK)
    ps.caption(fig, "Main-sequence pairs within 200 pc from El-Badry, Rix & Heintz (2021), R < 0.1; masses from "
               "Mamajek's dwarf sequence. Panel a shows raw counts: at small separations faint companions are "
               "harder to resolve, which tilts the raw histogram towards q = 1; the fits in panel b correct for it.",
               y=0.0)
    return ps.save(fig, _out("06_twins"))


# ---------------------------------------------------------------------------
# Gaia's binaries on the HR diagram
# ---------------------------------------------------------------------------

HRD_PANELS = [
    ("sb2", "Double-lined (SB2)", ps.BLUE),
    ("sb1", "Single-lined (SB1)", ps.BLUE),
    ("astrometric", "Astrometric orbit", ps.ORANGE),
    ("acceleration", "Acceleration / RV trend", ps.YELLOW),
    ("ruwe", "RUWE > 1.4, no solution", ps.VIOLET),
]


def fig_hr_diagram(res) -> list:
    h = res["hrd"]
    fig = plt.figure(figsize=(14.0, 8.6))
    gs = fig.add_gridspec(2, 5, height_ratios=[1.35, 1.0], hspace=0.34, wspace=0.08)
    x, y, lab = h["bp_rp"], h["abs_g"], h["label"]
    for k, (ch, title, col) in enumerate(HRD_PANELS):
        ax = fig.add_subplot(gs[0, k])
        ax.hist2d(x, y, bins=[np.arange(0.3, 3.6, 0.02), np.arange(-1, 15.5, 0.05)], cmap="Greys",
                  norm=LogNorm(vmin=1, vmax=300), alpha=0.45, rasterized=True)
        m = lab == ch
        ax.scatter(x[m], y[m], s=2.5, color=col, alpha=0.75, linewidths=0, rasterized=True)
        ax.plot(h["ridge_c"], h["ridge_g"], color=ps.INK, lw=0.8)
        ax.plot(h["ridge_c"], h["ridge_g"] - 0.753, color=ps.INK, lw=0.8, ls="--")
        ax.set_xlim(0.3, 3.5)
        ax.set_ylim(15, -0.5)
        ax.set_title(f"{title}\n{int(m.sum()):,} stars", fontsize=9, loc="left")
        ax.set_xlabel("BP − RP")
        if k == 0:
            ax.set_ylabel("absolute magnitude M_G")
            ps.annotate(ax, "single-star ridge", (2.6, float(np.interp(2.6, h["ridge_c"], h["ridge_g"]))), (2.0, 13.8))
            ps.annotate(ax, "two equal stars\n(+0.75 mag)", (1.6, float(np.interp(1.6, h["ridge_c"], h["ridge_g"])) - 0.753),
                        (0.45, 6.2))
        else:
            ax.set_yticklabels([])

    # offsets above the ridge, K and early-M dwarfs (where the single-star sequence is narrowest)
    ax = fig.add_subplot(gs[1, :3])
    sel = h["ms_band"] & (x >= 1.0) & (x < 1.8)
    edges = np.arange(-1.6, 0.81, 0.04)
    rows = [("none", "no Gaia sign of a companion", ps.CONTEXT, "-"),
            ("sb1", "single-lined", ps.BLUE, "-"), ("astrometric", "astrometric orbit", ps.ORANGE, "-"),
            ("sb2", "double-lined", ps.INK, "-")]
    for ch, name, col, ls_ in rows:
        d = h["delta"][sel & (lab == ch)]
        d = d[np.isfinite(d)]
        hh, _ = np.histogram(d, edges)
        dens = hh / max(len(d), 1) / np.diff(edges)
        ax.step(edges[:-1], dens, where="post", color=col, lw=1.7, ls=ls_, label=f"{name} ({len(d):,} stars)")
    ax.axvline(-0.753, color=ps.MUTED, lw=0.8, ls="--")
    ax.text(-0.76, ax.get_ylim()[1] * 0.55, "two equal\nstars ", fontsize=7.5, color=ps.MUTED, va="top", ha="right")
    ax.set_xlim(-1.6, 0.8)
    ax.set_xlabel("height above the single-star ridge, ΔM_G (mag; negative = brighter)")
    ax.set_ylabel("fraction of stars per mag")
    ax.set_title("b  K and early-M dwarfs: height above the ridge")
    ax.legend(loc="upper left", fontsize=7.5)

    ax = fig.add_subplot(gs[1, 3:])
    st = h["stats"]["KM"]["twin_zone_composition"]
    names = {"none": "no sign", "ruwe": "RUWE only", "acceleration": "acceleration", "astrometric": "astrometric",
             "sb1": "single-lined", "sb2": "double-lined"}
    order = ["ruwe", "none", "sb1", "sb2", "acceleration", "astrometric"]
    vals = [st[k] for k in order]
    ax.barh(range(len(order)), vals, color=ps.BLUE, height=0.6)
    for i, v in enumerate(vals):
        ax.text(v + 0.01, i, f"{100 * v:.0f}%", va="center", fontsize=8, color=ps.INK_SECONDARY)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([names[k] for k in order])
    ax.invert_yaxis()
    ax.set_xlim(0, max(vals) * 1.25)
    ax.set_xlabel("share of stars in the twin zone (0.6–0.9 mag above the ridge)")
    ax.set_title(f"c  Who is on the twin sequence? ({h['stats']['KM']['twin_zone_n']:,} stars)", loc="left")
    ps.caption(fig, "Gaia DR3, every source within 100 pc with parallax/error > 10 (no RUWE cut). Grey: all stars. "
               "Solid line: the most common M_G at each colour among stars with no sign of a companion; dashed: "
               "0.753 mag brighter, where two identical stars appear as one. 'Single-lined' includes DR3's combined "
               "astrometric + spectroscopic orbits. Each star appears in one panel only, in this order of precedence: "
               "double-lined, single-lined, astrometric, acceleration, RUWE.", y=0.02)
    return ps.save(fig, _out("03_hr_diagram"))


# ---------------------------------------------------------------------------
# Detection windows and the period distribution
# ---------------------------------------------------------------------------

def _au_axis(ax, m_tot=1.5):
    """Top axis: semi-major axis (au) for a 1.5 Msun pair, from Kepler's law."""
    import binaries as B
    ticks_au = [0.03, 0.1, 1, 10, 100, 1000, 10_000, 100_000]
    ticks = np.log10(B.period_days_from_a(np.array(ticks_au), m_tot))
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    keep = (ticks > ax.get_xlim()[0]) & (ticks < ax.get_xlim()[1])
    top.set_xticks(ticks[keep])
    top.set_xticklabels([f"{t:g}" for t, k in zip(ticks_au, keep) if k])
    top.set_xlabel("orbit size for a Sun-like pair (au)", color=ps.INK_SECONDARY)
    top.grid(False)
    top.spines["top"].set_visible(True)
    top.spines["top"].set_color(ps.AXIS)
    return top


def fig_detection_windows(res) -> list:
    sc = res["solar"]
    comp, ideal = sc["completeness"], sc["completeness_ideal"]
    x = comp["centres"]
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 8.2), sharex=True, gridspec_kw={"hspace": 0.12})
    ax = axes[0]
    rows = [("spectroscopic", ps.BLUE, "o", "spectroscopic orbit"), ("astrometric", ps.ORANGE, "s", "astrometric orbit"),
            ("wide", ps.AQUA, "D", "resolved wide pair")]
    for key, col, mk, lab in rows:
        ax.plot(x, comp[key], "-", color=col, lw=2.0, marker=mk, ms=4, label=lab)
    ax.fill_between(x, comp["spectroscopic"], ideal["spectroscopic"], color=ps.BLUE, alpha=0.12, lw=0,
                    label="spectroscopic, if DR3 had kept every orbit passing its documented cuts")
    ax.plot(x, comp["period_known"], color=ps.INK, lw=1.0, ls="--", label="any channel that measures a period")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("fraction of companions detected")
    ax.set_title("a  Channels that measure the period", pad=34)
    ax.legend(loc="center right", fontsize=7.5)
    ax.set_xlim(-0.5, 9.0)
    _au_axis(ax)
    for xx, yy, txt in ((0.5, 1.0, "close binaries\n(spectroscopy)"), (2.5, 1.0, "Gaia's 34 months\n(astrometry)"),
                        (4.1, 1.0, "the gap"), (8.0, 0.9, "wide pairs\n(resolved by Gaia)")):
        ax.text(xx, yy, txt, ha="center", va="top", fontsize=7.5, color=ps.INK_SECONDARY)

    ax = axes[1]
    ax.plot(x, comp["acceleration"], "-", color=ps.BLUE, lw=2.0, marker="^", ms=4,
            label="astrometric acceleration (curved proper motion)")
    ax.plot(x, comp["ruwe_only"], "-", color=ps.ORANGE, lw=2.0, marker="v", ms=4,
            label="RUWE > 1.4 but no published solution")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("fraction of companions detected")
    ax.set_xlabel("orbital period, log10(P / day)")
    ax.set_title("b  Signs of a companion without a period")
    ax.legend(loc="upper right", fontsize=7.5)
    ps.caption(fig, f"Simulated companions (300 per 0.25 dex of period) around the {sc['n_primaries']:,} F6–K3 dwarfs "
               "within 100 pc, with Moe & Di Stefano's (2017) mass ratios and eccentricities and random orientations. "
               "Astrometry: gaiamock (El-Badry et al. 2024) with DR3's cuts; spectroscopy: DR3's documented SB1 "
               "selection, scaled by the acceptance rate measured on SB9 binaries; wide pairs: the contrast curve "
               "of figure 8.", y=0.0)
    return ps.save(fig, _out("04_detection_windows"))


def fig_period_distribution(res) -> list:
    sc = res["solar"]
    r10 = res["raghavan"]
    comp = sc["completeness"]
    x = comp["centres"]
    fig = plt.figure(figsize=(13.0, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1.0], wspace=0.25)
    ax = fig.add_subplot(gs[0, 0])
    n, width = sc["n_primaries"], np.diff(comp["edges"])
    raw = (sc["h_orbit"] + sc["h_wide"]) / n / width
    ax.step(comp["edges"][:-1], raw, where="post", color=ps.CONTEXT, lw=1.2, label="seen by Gaia DR3 (no correction)")
    f, fe, fi = sc["f_corr"], sc["f_err"], sc["f_ideal"]
    lo_band = np.fmin(f, fi)
    hi_band = np.fmax(f, fi)
    ok = np.isfinite(f)
    # shade the efficiency range only within contiguous runs of measured bins
    runs = np.split(np.where(ok)[0], np.where(np.diff(np.where(ok)[0]) > 1)[0] + 1)
    for k, run in enumerate(runs):
        ax.fill_between(x[run], lo_band[run], hi_band[run], color=ps.BLUE, alpha=0.18, lw=0,
                        label="range: DR3 spectroscopic acceptance as measured on SB9 ↔ as documented"
                        if k == 0 else None)
    cpk = comp["period_known"]
    solid = ok & (cpk >= 0.5)
    hollow = ok & (cpk < 0.5)
    ax.errorbar(x[solid], f[solid], yerr=fe[solid], fmt="o", color=ps.BLUE, ms=4.5, lw=1.0, capsize=0,
                label="Gaia DR3, corrected for completeness (this work)")
    ax.errorbar(x[hollow], f[hollow], yerr=fe[hollow], fmt="o", mfc=ps.SURFACE, color=ps.BLUE, ms=4.5, lw=1.0,
                capsize=0, label="… where Gaia catches fewer than half (less certain)")
    e = r10["hist_edges"]
    ax.step(e[:-1], r10["per_primary_per_dex"], where="post", color=ps.INK, lw=1.6,
            label=f"Raghavan et al. 2010: their {r10['n_companions']} companions, 25 pc")
    xx = np.linspace(-1, 10, 400)
    gauss = lambda m, s: np.exp(-0.5 * ((xx - m) / s) ** 2) / (s * np.sqrt(2 * np.pi))
    ax.plot(xx, 259 / 454 * gauss(L.RAGHAVAN_LOGP_MEAN, L.RAGHAVAN_LOGP_SIGMA), color=ps.INK, lw=1.2, ls="--",
            label="Raghavan's log-normal (log P = 5.03, σ = 2.28)")
    mds = np.where((x >= 0.2) & (x < 8.0), sc["models"]["mds"], np.nan)
    ax.plot(x, mds, color=ps.MUTED, lw=1.4, ls="-.", label="Moe & Di Stefano 2017 (q > 0.1; fitted to log P = 8)")
    ax.plot(xx, L.TOKOVININ_ALL["eps"] * gauss(L.TOKOVININ_ALL["x0"], L.TOKOVININ_ALL["sigma"]), color=ps.MUTED,
            lw=1.2, ls=":", label="Tokovinin 2014 (67 pc)")
    gap = ~ok & (x > 3) & (x < 5)
    if gap.any():
        ax.axvspan(x[gap].min() - 0.125, x[gap].max() + 0.125, color=ps.WASH, zorder=0)
        ax.text(x[gap].mean() - 0.12, 0.006, "Gaia sees these\nonly as accelerations\nand RUWE excess\n(panel b)",
                ha="center", va="bottom", fontsize=7.2, color=ps.INK_SECONDARY)
    ax.set_xlim(-0.5, 9.0)
    ax.set_ylim(0, 0.135)
    ax.set_xlabel("orbital period, log10(P / day)")
    ax.set_ylabel("companions per star per decade of period")
    ax.set_title("a  How far apart are the companions of Sun-like stars?", pad=34)
    ax.legend(loc="upper left", fontsize=7.0, ncol=2, columnspacing=1.0)
    _au_axis(ax)

    ax = fig.add_subplot(gs[0, 1])
    chans = [("orbit", "orbits (all types)"), ("acceleration", "accelerations"), ("ruwe_only", "RUWE > 1.4 only"),
             ("wide", "resolved pairs")]
    obs = sc["observed_counts"]
    styles = [("raghavan", ps.BLUE, "o", "Raghavan 2010"), ("mds", ps.ORANGE, "s", "Moe & Di Stefano 2017"),
              ("tokovinin", ps.AQUA, "D", "Tokovinin 2014")]
    for i, (ch, lab) in enumerate(chans):
        y0 = len(chans) - 1 - i
        for j, (name, col, mk, mlab) in enumerate(styles):
            a = sc["predicted"][name][ch]
            b = sc["predicted_ideal"][name][ch]
            yy = y0 + 0.18 * (j - 1)
            ax.plot([min(a, b), max(a, b)], [yy, yy], color=col, lw=2.5, alpha=0.5)
            ax.plot(a, yy, mk, color=col, ms=6, label=mlab if i == 0 else None)
        ax.plot(obs[ch], y0, "|", color=ps.INK, ms=22, mew=2.4)
        ax.text(obs[ch], y0 + 0.36, f"{obs[ch]:,}", ha="center", fontsize=7.5, color=ps.INK)
    ax.set_yticks(range(len(chans)))
    ax.set_yticklabels([c[1] for c in chans][::-1])
    ax.set_xlabel("number of the 32,487 stars' companions")
    ax.set_xlim(0, 8200)
    ax.set_ylim(-0.6, len(chans) - 0.3)
    ax.set_title("b  Gaia's raw counts vs. published distributions")
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], color=ps.INK, marker="|", ls="none", ms=10, mew=2.2))
    labels.append("observed in Gaia DR3")
    ax.legend(handles, labels, loc="lower right", fontsize=7.2)
    ps.caption(fig, "Sun-like stars (F6–K3 dwarfs) within 100 pc. Orbits carry Gaia's period; resolved pairs a period "
               "from their separation, as Raghavan et al. estimated theirs (a = 1.35 × separation). Panel b pushes "
               "each published period distribution through the same simulated Gaia selection; bars span the two "
               "assumptions about DR3's spectroscopic acceptance.", y=-0.02)
    return ps.save(fig, _out("05_period_distribution"))


# ---------------------------------------------------------------------------
# Period-eccentricity diagrams and the circularization clock
# ---------------------------------------------------------------------------

def _mm_curve(p_prime, p=None):
    import binaries as B
    if p is None:
        p = np.logspace(-0.5, 3.5, 400)
    return p, B.meibom_mathieu_e(p, p_prime)


def fig_period_eccentricity(res) -> list:
    import binaries as B
    s9 = res["sb9"]
    g = res["gaia_circ"]
    t = s9["table"]
    fig = plt.figure(figsize=(13.5, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.32, wspace=0.18)
    pp = np.logspace(-0.6, 3.6, 400)
    env = B.e_max_envelope(pp)

    ax = fig.add_subplot(gs[0, 0])
    sel = s9["fgk_dwarf"]
    P, E = np.asarray(t["period"])[sel], np.asarray(t["ecc"])[sel]
    fixed = np.asarray(t["ecc_fixed"])[sel]
    ax.scatter(np.log10(P[~fixed]), E[~fixed], s=9, color=ps.BLUE, alpha=0.7, linewidths=0, label="fitted eccentricity")
    ax.scatter(np.log10(P[fixed]), E[fixed], s=12, facecolors="none", edgecolors=ps.BLUE, linewidths=0.8,
               label="eccentricity fixed at 0 by the authors")
    ax.plot(np.log10(pp), env, color=ps.MUTED, lw=1.0, ls="--", label="e_max = 1 − (P/2 d)^−2/3 (Moe & Di Stefano 2017)")
    p_c, e_c = _mm_curve(s9["mm"]["p_prime"], pp)
    ax.plot(np.log10(p_c), e_c, color=ps.INK, lw=1.6, label=f"Meibom & Mathieu fit: P_circ = {s9['p_circ']:.1f} d")
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(-0.03, 1.0)
    ax.set_xlabel("orbital period, log10(P / day)")
    ax.set_ylabel("eccentricity")
    ax.set_title(f"a  SB9: {sel.sum():,} binaries with F, G, K dwarf primaries")
    ax.legend(loc="upper left", fontsize=7.3)

    ax = fig.add_subplot(gs[0, 1])
    lp = np.log10(g["P"])
    e_plot = np.where(g["circ"], 0.0, g["e"])
    h = ax.hist2d(lp, e_plot, bins=[np.arange(-0.6, 3.61, 0.04), np.arange(-0.02, 1.0, 0.012)], cmap=SEQ_BLUE,
                  norm=LogNorm(vmin=1, vmax=300), rasterized=True)
    ax.plot(np.log10(pp), env, color=ps.MUTED, lw=1.0, ls="--")
    # the fit stores P_circ = e(P) = 0.01 crossing; the curve needs P' = P_circ - (P_circ - P') at P' = 0
    p_prime = g["all"]["p_circ"] - B.meibom_mathieu_pcirc(0.0)
    p_c, e_c = _mm_curve(p_prime, pp)
    ax.plot(np.log10(p_c), e_c, color=ps.INK, lw=1.6,
            label=f"Meibom & Mathieu fit: P_circ = {g['all']['p_circ']:.1f} d")
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(-0.03, 1.0)
    ax.set_xlabel("orbital period, log10(P / day)")
    ax.set_title(f"b  Gaia DR3: {g['n']:,} single-lined orbits, dwarfs < 500 pc")
    ax.legend(loc="upper left", fontsize=7.3)
    cb = fig.colorbar(h[3], ax=ax, pad=0.01, fraction=0.04)
    cb.set_label("orbits per cell", color=ps.INK_SECONDARY, fontsize=8)
    cb.outline.set_visible(False)
    ax.text(3.55, 0.012, "orbits consistent with circular (e < 2.45σₑ) drawn at e = 0", ha="right", va="bottom",
            fontsize=7, color=ps.INK_SECONDARY)

    ax = fig.add_subplot(gs[1, :])
    for cf, col, mk, lab in ((s9["circ_frac"], ps.BLUE, "o", "SB9 F, G, K dwarfs"),
                             (g["circ_frac"], ps.ORANGE, "s", "Gaia DR3 dwarfs within 500 pc")):
        e = cf["edges"]
        c = 0.5 * (e[1:] + e[:-1])
        ok = cf["n"] >= 8
        err = np.sqrt(np.clip(cf["frac"] * (1 - cf["frac"]), 1e-4, None) / np.maximum(cf["n"], 1))
        ax.errorbar(c[ok], cf["frac"][ok], yerr=err[ok], fmt=mk + "-", color=col, ms=4.5, lw=1.4, capsize=0,
                    label=f"{lab}: half circular at {cf['p50']:.1f} d")
    ax.axvline(np.log10(L.RAGHAVAN_PCIRC_DAYS), color=ps.MUTED, lw=0.9, ls=":")
    ax.text(np.log10(L.RAGHAVAN_PCIRC_DAYS) + 0.02, 0.55, "Raghavan et al.\n2010: ~12 d", fontsize=7.5,
            color=ps.INK_SECONDARY, va="top")
    ax.set_xlim(-0.4, 2.0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("orbital period, log10(P / day)")
    ax.set_ylabel("fraction of orbits consistent with circular")
    ax.set_title("c  Tides have made the closest orbits circular")
    ax.legend(loc="upper right", fontsize=7.5)
    ps.caption(fig, "SB9 (Pourbaix et al. 2004, VizieR B/sb9): one orbit per system, dwarfs by Gaia absolute magnitude "
               "and spectral class. Gaia DR3 single-lined and astrometric-spectroscopic orbits with significance > 20 "
               "(DR3's double-lined eccentricities are unreliable below ~3 days). "
               "An eccentricity counts as consistent with circular if e < 2.45 σₑ (Lucy & Sweeney 1971).", y=0.02)
    return ps.save(fig, _out("01_period_eccentricity"))


def fig_circularization_clock(res) -> list:
    g = res["gaia_circ"]
    import analysis as A
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.4), gridspec_kw={"wspace": 0.22})
    ax = axes[0]
    for k, (name, age, cp, hi, lo) in enumerate(A.MM05_TABLE3):
        ax.errorbar(age, cp, yerr=[[lo], [hi]], fmt="s", mfc="none", mec=ps.INK, color=ps.INK, ms=6, lw=1.0, capsize=0,
                    label="Meibom & Mathieu 2005 (clusters, field, halo)" if k == 0 else None)
        ax.annotate(name, (age, cp), xytext=(6, -3), textcoords="offset points", fontsize=7, color=ps.INK_SECONDARY)
    gc = g["clusters"]["< 1 Gyr, 5300-6300 K"]
    if np.isfinite(gc["p_circ"]):
        ax.errorbar(gc["median_age_gyr"], gc["p_circ"], yerr=[[gc["p_circ"] - gc["lo"]], [gc["hi"] - gc["p_circ"]]],
                    fmt="o", color=ps.AQUA, ms=7, lw=1.4, capsize=0,
                    label=f"this work: open clusters < 1 Gyr ({gc['n']} binaries)")
    pops = g["populations"]
    for key, (age_lo, age_hi), mk, col, lab in (("thin disc (Teff-matched)", (1.0, 8.0), "o", ps.BLUE, "thin disc"),
                                                ("thick disc", (8.0, 12.0), "o", ps.ORANGE, "thick disc")):
        r = pops[key]
        age = np.sqrt(age_lo * age_hi)
        ax.errorbar(age, r["p_circ"], xerr=[[age - age_lo], [age_hi - age]],
                    yerr=[[r["p_circ"] - r["lo"]], [r["hi"] - r["p_circ"]]], fmt=mk, color=col, ms=7, lw=1.4, capsize=0,
                    label=f"this work: {lab} ({r['n']:,} binaries, 5300–6300 K)")
    ax.set_xscale("log")
    ax.set_xlabel("age (Gyr)")
    ax.set_ylabel("circularization period (days)")
    ax.set_ylim(0, 20)
    ax.set_title("a  Is there a circularization clock?")
    ax.legend(loc="upper left", fontsize=7.2)

    ax = axes[1]
    tb, te = g["teff_bins"], g["teff_bins_envelope"]
    tc = np.array([(b["teff_lo"] + b["teff_hi"]) / 2 for b in tb])
    def _plot(rows, col, mk, lab, dx):
        y = np.array([r["p_circ"] if "p_circ" in r else r["p_cut"] for r in rows])
        lo = np.array([r["lo"] for r in rows])
        hi = np.array([r["hi"] for r in rows])
        ok = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi)
        ax.errorbar(tc[ok] + dx, y[ok], yerr=[np.clip(y[ok] - lo[ok], 0, None), np.clip(hi[ok] - y[ok], 0, None)],
                    fmt=mk, color=col, ms=5, lw=1.1, capsize=0, label=lab)
    _plot(tb, ps.BLUE, "o", "this work: Meibom–Mathieu circularization period", -12)
    _plot(te, ps.ORANGE, "s", "this work: envelope base P_cut (Bashi et al.'s method, significance > 40)", 12)
    ax.errorbar(tc, A.BASHI_PCUT, yerr=A.BASHI_PCUT_ERR, fmt="D", mfc="none", mec=ps.INK, color=ps.INK, ms=5, lw=1.0,
                capsize=0, label="Bashi, Mazeh & Faigler 2023 (P_cut)")
    ax.set_xlabel("effective temperature of the primary (K)")
    ax.set_ylabel("period (days)")
    ax.set_ylim(0, 16)
    ax.set_title("b  Cooler stars, with deeper convective envelopes, circularize wider orbits")
    ax.legend(loc="upper right", fontsize=7.2)
    ps.caption(fig, "Gaia DR3 single-lined orbits of dwarfs within 500 pc (significance > 20), all with primaries of "
               "5300–6300 K in panel a. Populations by space velocity relative to the local standard of rest (thin "
               "disc < 50 km/s, thick disc 70–180 km/s); the thin disc is resampled to the thick disc's temperatures. "
               "Clusters: Hunt & Reffert (2023) members. The two estimators in panel b measure different things "
               "(where a fitted curve reaches e = 0.01, vs. where eccentric orbits begin).", y=-0.03)
    return ps.save(fig, _out("02_circularization_clock"))


# ---------------------------------------------------------------------------
# Multiplicity against primary mass
# ---------------------------------------------------------------------------

def _plain_mass_ticks(ax, ticks):
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    lo, hi = ax.get_xlim()
    t = [x for x in ticks if lo <= x <= hi]
    ax.xaxis.set_major_locator(FixedLocator(t))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{x:g}" for x in t]))
    ax.xaxis.set_minor_locator(NullLocator())


def fig_multiplicity_vs_mass(res) -> list:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.2), gridspec_kw={"wspace": 0.28})
    ax = axes[0]
    for k, ((lo, hi), mf, mfe, cf, cfe, src) in enumerate(L.OFFNER23_MF):
        m = np.sqrt(lo * hi)
        ax.errorbar(m, mf, xerr=[[m - lo], [hi - m]], yerr=mfe, fmt="o", color=ps.BLUE, ms=6, lw=1.2, capsize=0,
                    label="multiplicity fraction: share of stars with ≥ 1 companion" if k == 0 else None)
        ax.errorbar(m, cf, yerr=cfe, fmt="s", color=ps.ORANGE, ms=6, lw=1.2, capsize=0,
                    label="companion frequency: companions per star" if k == 0 else None)
    ax.axhline(1.0, color=ps.AXIS, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("primary mass (M☉)")
    ax.set_ylabel("fraction or companions per star")
    ax.set_title("a  Literature: companions rise with mass")
    ax.legend(loc="upper left", fontsize=7.2)
    for m, lab in ((0.2, "M"), (0.7, "K"), (1.0, "G"), (1.3, "F"), (2.0, "A"), (5, "B"), (25, "O")):
        ax.text(m, 0.02, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8,
                color=ps.MUTED)
    ax.set_ylim(0, 2.5)
    _plain_mass_ticks(ax, [0.1, 0.3, 1, 3, 10, 30])

    ax = axes[1]
    rows = res["close"]
    m = np.array([r["m_med"] for r in rows])
    a = np.array([r["close_sb9"] for r in rows])
    b = np.array([r["close_ideal"] for r in rows])
    ea = np.array([r["close_sb9_err"] for r in rows])
    for i in range(len(rows)):
        ax.plot([m[i], m[i]], [min(a[i], b[i]), max(a[i], b[i])], color=ps.BLUE, lw=6, alpha=0.25, solid_capstyle="butt")
    ax.errorbar(m, a, yerr=ea, fmt="o", color=ps.BLUE, ms=6, lw=1.2, capsize=0,
                label="this work: Gaia DR3 orbits, corrected (bar: acceptance range)")
    mm = np.array([r["mds_close"] for r in rows])
    ok = m >= 0.8
    ax.plot(m[ok], mm[ok], "s--", color=ps.INK, mfc="none", ms=6, lw=1.0,
            label="Moe & Di Stefano 2017 fit (valid ≥ 0.8 M☉)")
    for r in rows:
        ax.annotate(f"{r['name']}\n{r['n_primaries']:,} stars", (r["m_med"], max(r["close_sb9"], r["close_ideal"])),
                    xytext=(0, 8), textcoords="offset points", ha="center", fontsize=7, color=ps.INK_SECONDARY)
    rag = res["raghavan"]["per_primary_per_dex"]
    ax.errorbar(1.0, rag[1:4].sum(), yerr=np.sqrt(rag[1:4].sum() * 454) / 454, fmt="D", mfc="none", mec=ps.INK,
                color=ps.INK, ms=6, lw=1.0, capsize=0, label="Raghavan et al. 2010 (their companions, 25 pc)")
    ax.set_xscale("log")
    ax.set_xlabel("primary mass (M☉)")
    ax.set_ylabel("companions per star with 1 d < P < 1000 d (q > 0.1)")
    ax.set_title("b  Close companions (1–1,000 days)")
    ax.set_ylim(0, None)
    _plain_mass_ticks(ax, [0.3, 0.5, 0.7, 1, 1.5, 2])
    ax.legend(loc="upper left", fontsize=7.0)

    ax = axes[2]
    w = res["wbf"]
    for k, r in enumerate(w):
        mid = r["m_med"]
        lower_limit = r["m_lo"] >= 1.5
        style = dict(fmt="o", color=ps.BLUE, ms=6, lw=1.2, capsize=0)
        if lower_limit:
            style.update(mfc=ps.SURFACE)
        ax.errorbar(mid, r["wbf"], yerr=r["wbf_err"], label=("this work (El-Badry+2021 pairs, corrected)" if k == 0 else
                                                            ("this work: lower limit (see text)" if lower_limit and
                                                             r is next(x for x in w if x["m_lo"] >= 1.5) else None)),
                    **style)
        if lower_limit:
            ax.annotate("", xy=(mid, r["wbf"] + 0.05), xytext=(mid, r["wbf"]),
                        arrowprops={"arrowstyle": "->", "color": ps.BLUE, "lw": 1.0})
    for k, ((lo, hi), v, e, src) in enumerate(L.OFFNER23_WBF):
        mm_ = np.sqrt(lo * hi)
        ax.errorbar(mm_ * 1.04, v, xerr=[[mm_ - lo], [hi - mm_]], yerr=e, fmt="s", mfc="none", mec=ps.INK, color=ps.INK,
                    ms=6, lw=1.0, capsize=0, label="Offner et al. 2023 compilation" if k == 0 else None)
    ax.set_xscale("log")
    ax.set_xlabel("primary mass (M☉)")
    ax.set_ylabel("share of stars with a companion at 100–10,000 au")
    ax.set_title("c  Wide companions (100–10,000 au)")
    ax.set_ylim(0, 0.7)
    _plain_mass_ticks(ax, [0.2, 0.5, 1, 2, 5, 10])
    ax.legend(loc="upper left", fontsize=7.0)
    ps.caption(fig, "Panel a: Offner et al. (2023), Table 1 (Winters+2019, Raghavan+2010, Tokovinin 2014, Moe & Kratter "
               "2021, Moe & Di Stefano 2017, Sana et al.). Panels b and c: this work. Above ~1.5 M☉ the wide pairs are "
               "lower limits: El-Badry et al.'s catalogue drops pairs in crowded regions of phase space, i.e. in the "
               "young associations where most nearby A and B stars live.", y=-0.04)
    return ps.save(fig, _out("07_multiplicity_vs_mass"))


def make_all(res) -> list:
    import interactive
    paths = []
    for fn in (fig_period_eccentricity, fig_circularization_clock, fig_hr_diagram, fig_detection_windows,
               fig_period_distribution, fig_twins, fig_multiplicity_vs_mass, fig_eyesight):
        paths += fn(res)
    paths.append(interactive.build(res))
    return paths
