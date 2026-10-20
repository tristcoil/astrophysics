"""Every static figure, drawn from the analysis results (PNG + SVG in figures/)."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator

import config as C
import literature as LIT
import physics as P
import plotstyle as S

S.use_style()

LOG_PIVOT = np.log10(C.SIGMA_PIVOT)
SIG_TICKS = [20, 30, 50, 70, 100, 150, 200, 300, 400]


def _sigma_axis(ax, lo=25, hi=480):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(FixedLocator([t for t in SIG_TICKS if lo <= t <= hi]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("stellar velocity dispersion σ (km/s)")


def _mass_axis(ax, lo=1e5, hi=4e10, label="black hole mass (M☉)"):
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"$10^{{{int(round(np.log10(v)))}}}$"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel(label)


def _fit_line(ax, fit, lo=25, hi=480, *, band=True, color=S.INK_SECONDARY, lw=2.0, zorder=2, **kw):
    s = np.geomspace(lo, hi, 200)
    y = fit.a + fit.b * (np.log10(s) - LOG_PIVOT)
    if band:
        ax.fill_between(s, 10 ** (y - fit.eps), 10 ** (y + fit.eps), color=S.WASH, lw=0, zorder=0.5)
    return ax.plot(s, 10**y, color=color, lw=lw, solid_capstyle="round", zorder=zorder, **kw)[0]


def _anchor(ax, sigma, mass, label, xytext, *, ha="left"):
    ax.plot([sigma], [mass], ls="none", zorder=6, **S.ANCHOR_STYLE)
    S.annotate(ax, label, (sigma, mass), xytext, ha=ha, va="center", fontsize=8, color=S.INK)


def _points(ax, df, style, *, xcol="sigma", ycol="M", s=26, errors=True, zorder=3, label=None):
    if len(df) == 0:
        return None
    x, y = df[xcol].to_numpy(float), df[ycol].to_numpy(float)
    if errors and "elogM" in df and "elogsig" in df:
        ey = np.nan_to_num(df["elogM"].to_numpy(float))
        ex = np.nan_to_num(df["elogsig"].to_numpy(float))
        ax.errorbar(x, y, yerr=[y - y / 10**ey, y * 10**ey - y], xerr=[x - x / 10**ex, x * 10**ex - x],
                    fmt="none", ecolor=style["color"], elinewidth=0.6, alpha=0.35, zorder=zorder - 1)
    return ax.scatter(x, y, s=s, c=style["color"], marker=style["marker"], edgecolors=S.SURFACE,
                      linewidths=0.8, zorder=zorder, label=label or style["label"])


def _upper_limits(ax, df, *, color=S.MUTED, s=22, label="upper limits"):
    """Carets at the limit itself (the top of the allowed range)."""
    if len(df) == 0:
        return None
    y = df["M_lim"] if "M_lim" in df else df["M"]
    return ax.scatter(df["sigma"], y, s=s, marker="v", c=color, edgecolors="none", alpha=0.75,
                      zorder=2.5, label=label)


def _with_linear(df):
    d = df.copy()
    d["sigma"] = 10 ** d["logsig"]
    d["M"] = 10 ** d["logM"]
    if "logM_limit" in d:
        d["M_lim"] = np.where(d["upper"], 10 ** d["logM_limit"], d["M"])
    return d


# ---------------------------------------------------------------------------
# 1. The keeper: the M-sigma relation
# ---------------------------------------------------------------------------

def fig_m_sigma(res, out):
    cat, fits = res["cat"], res["fits"]
    v = _with_linear(cat["vdb16"])
    vs = v[v["selected"] == 1]
    fit = fits["vdb16"]
    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    _fit_line(ax, fit, 16, 480)
    handles = []
    for meth in ("stars", "gas", "maser", "reverb"):
        d = vs[(vs["method"] == meth) & ~vs["upper"]]
        handles.append(_points(ax, d, S.METHOD_STYLE[meth], zorder=3 if meth != "reverb" else 2.6))
    handles.append(_upper_limits(ax, vs[vs["upper"]], label="upper limits (stars or gas)"))
    mw, m87 = C.MILKY_WAY, C.M87
    ax.errorbar([mw["sigma"]], [mw["M"]], xerr=[[mw["sigma_err"]], [mw["sigma_err"]]], fmt="none",
                ecolor=S.INK, elinewidth=0.9, zorder=5)
    ax.errorbar([m87["sigma"]], [m87["M"]], xerr=[[m87["sigma_err_lo"]], [m87["sigma_err_hi"]]],
                yerr=[[m87["M_err"]], [m87["M_err"]]], fmt="none", ecolor=S.INK, elinewidth=0.9, zorder=5)
    _anchor(ax, mw["sigma"], mw["M"], "Milky Way: Sgr A*,\nweighed by star S2's orbit (H1)", (250, 2.5e5),
            ha="center")
    _anchor(ax, m87["sigma"], m87["M"], "M87*: weighed by the size\nof its shadow (EHT)", (150, 1.2e10), ha="right")
    for name, xy in (("NGC 4889", (210, 2.4e10)), ("NGC 4258", (55, 1.3e8)), ("NGC 221", (34, 4e6)),
                     ("NGC 598", (22, 3e4))):
        r = vs[vs["name"] == name]
        if len(r):
            lab = {"NGC 221": "M32", "NGC 4258": "NGC 4258 (maser disk)", "NGC 4889": "NGC 4889",
                   "NGC 598": "M33: no black hole\nabove 1,400 M☉"}[name]
            S.annotate(ax, lab, (r["sigma"].iloc[0], r["M_lim"].iloc[0]), xy, ha="center", fontsize=7.5)
    handles.append(Line2D([], [], color=S.INK_SECONDARY, lw=2,
                          label=f"our fit: slope {fit.b:.2f} ± {fit.eb:.2f}, scatter {fit.eps:.2f} dex"))
    handles.append(Line2D([], [], ls="none", **{k: v for k, v in S.ANCHOR_STYLE.items()},
                          label="weighed without a galaxy model"))
    ax.legend(handles=handles, loc="upper left", fontsize=7.5, handletextpad=0.4, borderaxespad=0.2)
    _sigma_axis(ax, 15, 480)
    _mass_axis(ax, 5e2, 6e10)
    ax.set_title("Black holes are heavier where their galaxy's stars move faster")
    S.caption(fig, f"{fit.n} galaxies (van den Bosch 2016 compilation; {fit.n - fit.n_upper} masses, {fit.n_upper} upper "
                   f"limits). Fit: log M = {fit.a:.2f} + {fit.b:.2f} log(σ/200 km/s), intrinsic scatter {fit.eps:.2f} dex "
                   "(shaded), Bayesian with censored limits. σ: stars' velocity dispersion within the half-light radius.")
    return S.save(fig, out / "01_m_sigma")


# ---------------------------------------------------------------------------
# 2. Four ways to weigh a black hole
# ---------------------------------------------------------------------------

def fig_methods(res, out):
    cat, fits = res["cat"], res["fits"]
    v = _with_linear(cat["vdb16"])
    vs = v[(v["selected"] == 1)]
    fit = fits["vdb16"]
    stats_ = res["methods"]
    notes = {
        "stars": "orbits of stars near the centre",
        "gas": "a rotating disk of gas",
        "maser": "radio-bright water clouds in a thin disk:\nmostly disk galaxies with pseudobulges",
        "reverb": "light echoes in active galaxies; the factor f\nis tuned so they sit on M-σ (not a test)",
    }
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.6), sharex=True, sharey=True)
    for ax, meth in zip(axes.flat, ("stars", "gas", "maser", "reverb")):
        _fit_line(ax, fit, 16, 480)
        others = vs[(vs["method"] != meth) & ~vs["upper"]]
        ax.scatter(others["sigma"], others["M"], s=12, c=S.GRID, edgecolors="none", zorder=1.5)
        st = dict(S.METHOD_STYLE[meth])
        if meth == "reverb":
            st["color"] = S.INK_SECONDARY
        _points(ax, vs[(vs["method"] == meth) & ~vs["upper"]], st, s=30)
        _upper_limits(ax, vs[(vs["method"] == meth) & vs["upper"]])
        m = stats_[meth]
        ax.set_title(S.METHOD_STYLE[meth]["label"].split(" (")[0].capitalize(), fontsize=10)
        ax.text(0.03, 0.97, f"{m['n']} masses; {notes[meth]}", transform=ax.transAxes, va="top", fontsize=7.5,
                color=S.INK_SECONDARY)
        ax.text(0.97, 0.04, f"offset {m['mean']:+.2f} ± {m['sem']:.2f} dex\nrms {m['rms']:.2f} dex",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=S.INK)
        _sigma_axis(ax, 15, 480)
        _mass_axis(ax, 5e2, 6e10)
        ax.label_outer()
    fig.suptitle("Four ways to weigh a black hole, one relation", x=0.06, ha="left", y=0.995)
    S.caption(fig, "Offsets are mean residuals (log M) from the fit to all 230 galaxies (line; shaded = intrinsic "
                   "scatter). That fit includes upper limits, so detections alone sit slightly above it on average.")
    fig.tight_layout()
    return S.save(fig, out / "02_methods")


# ---------------------------------------------------------------------------
# 3. Classical bulges vs pseudobulges
# ---------------------------------------------------------------------------

def fig_bulges(res, out):
    cat, fits = res["cat"], res["fits"]
    kh = _with_linear(cat["kh13"])
    fit = fits["kh13"]
    b = res["bulges"]
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(10.6, 5.2), gridspec_kw={"width_ratios": [1.55, 1]})
    _fit_line(ax, fit)
    fitted = kh[kh["kh13_fit"]]
    for kind in ("elliptical", "classical"):
        _points(ax, fitted[fitted["bulge"] == kind], S.BULGE_STYLE[kind], s=30)
    pb = kh[(kh["bulge"] == "pseudo")]
    _points(ax, pb[~pb["upper"]], S.BULGE_STYLE["pseudo"], s=30)
    _upper_limits(ax, pb[pb["upper"]], color=S.AQUA, label="pseudobulge upper limits")
    rest = kh[kh["bulge"].isin(["elliptical", "classical"]) & ~kh["kh13_fit"] & kh["sigma"].notna()]
    ax.scatter(rest["sigma"], rest["M"], s=22, marker="x", c=S.MUTED, linewidths=1.0, zorder=2.5,
               label="set aside by Kormendy & Ho")
    for name, lab, xy in (("NGC 1277", "NGC 1277 'monster'", (190, 4e10)),
                          ("NGC 2960", "NGC 2960 (merger)", (330, 2.5e7))):
        r = kh[kh["name"] == name].iloc[0]
        S.annotate(ax, lab, (r["sigma"], r["M"]), xy, ha="center", fontsize=7.5)
    mw = C.MILKY_WAY
    ax.errorbar([mw["sigma"]], [mw["M"]], xerr=[[mw["sigma_err"]], [mw["sigma_err"]]], fmt="none",
                ecolor=S.INK, elinewidth=0.9, zorder=5)
    _anchor(ax, mw["sigma"], mw["M"], "Milky Way (a pseudobulge)", (230, 6e5), ha="left")
    ax.legend(loc="upper left", fontsize=7.5, handletextpad=0.3)
    _sigma_axis(ax, 50, 480)
    _mass_axis(ax, 2e5, 6e10)
    ax.set_title("Ellipticals and classical bulges follow it; pseudobulges don't")

    # Right: residuals from the classical relation, one strip per kind.
    kh["res"] = kh["logM"] - fit.predict(kh["logsig"] - LOG_PIVOT)
    groups = [("ellipticals", fitted[fitted["bulge"] == "elliptical"], S.BULGE_STYLE["elliptical"]),
              ("classical\nbulges", fitted[fitted["bulge"] == "classical"], S.BULGE_STYLE["classical"]),
              ("pseudo-\nbulges", pb[~pb["upper"]], S.BULGE_STYLE["pseudo"]),
              ("mergers in\nprogress", kh[kh["merger"] & kh["sigma"].notna() & ~kh["upper"]],
               {"color": S.MUTED, "marker": "x"})]
    rng = np.random.default_rng(3)
    axr.axhspan(-fit.eps, fit.eps, color=S.WASH, lw=0, zorder=0)
    axr.axhline(0, color=S.INK_SECONDARY, lw=1.2, zorder=1)
    for i, (lab, d, st) in enumerate(groups):
        r = (d["logM"] - fit.predict(d["logsig"] - LOG_PIVOT)).to_numpy()
        jit = rng.uniform(-0.18, 0.18, len(r))
        kw = {"edgecolors": S.SURFACE, "linewidths": 0.7} if st["marker"] != "x" else {"linewidths": 1.0}
        axr.scatter(i + jit, r, s=24, c=st["color"], marker=st["marker"], zorder=3, **kw)
        axr.errorbar([i + 0.33], [r.mean()], yerr=[r.std(ddof=1) / np.sqrt(len(r))], fmt="o", ms=5,
                     color=S.INK, elinewidth=1.2, capsize=0, zorder=4)
        axr.text(i + 0.42, r.mean(), f"{r.mean():+.2f}", ha="left", va="center", fontsize=7.5, color=S.INK)
    axr.set_xticks(range(len(groups)))
    axr.set_xticklabels([g[0] for g in groups], fontsize=8)
    axr.set_xlim(-0.5, len(groups) - 0.15)
    axr.set_ylabel("log M − classical relation (dex)")
    axr.set_ylim(-2.0, 1.4)
    axr.grid(axis="x", visible=False)
    p = b["pseudo"]
    axr.set_title("Offsets from the classical relation", fontsize=10)
    axr.text(0.02, 0.02, f"pseudobulges: black holes {1 / p['median_mass_ratio_to_classical']:.0f}× lighter (median);\n"
                         f"{p['frac_below'] * 100:.0f}% below the line; no M-σ trend of their own "
                         f"(r = {p['pearson_r']:.2f}, p = {p['pearson_p']:.2f})",
             transform=axr.transAxes, fontsize=7.3, color=S.INK_SECONDARY, va="bottom")
    S.caption(fig, "Kormendy & Ho (2013) Tables 2-3, with their classifications. Line: our fit to their 45 ellipticals "
                   f"and classical bulges, log M = {fit.a:.2f} + {fit.b:.2f} log(σ/200), scatter {fit.eps:.2f} dex "
                   "(shaded). Black dots: group means ± standard error. Set aside: mergers, ionised-gas masses "
                   "without line widths, 'monsters' and the two largest black holes, as in the paper.")
    fig.tight_layout()
    return S.save(fig, out / "03_bulges")


# ---------------------------------------------------------------------------
# 4. Three published slopes, and where the difference comes from
# ---------------------------------------------------------------------------

def fig_slopes(res, out):
    cat, fits = res["cat"], res["fits"]
    sb = res["slope_budget"]
    fig = plt.figure(figsize=(11.0, 9.4))
    top = fig.add_gridspec(1, 2, left=0.07, right=0.98, top=0.95, bottom=0.57, width_ratios=[1.35, 1.0], wspace=0.22)
    bot = fig.add_gridspec(1, 1, left=0.31, right=0.98, top=0.45, bottom=0.08)
    ax, axj = fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1])
    axb = fig.add_subplot(bot[0, 0])

    v = _with_linear(cat["vdb16"])
    vs = v[(v["selected"] == 1) & ~v["upper"]]
    ax.scatter(vs["sigma"], vs["M"], s=12, c=S.CONTEXT, edgecolors="none", zorder=1,
               label=f"all {len(vs)} masses (van den Bosch)")
    lines = [("kh13_eq3", "Kormendy & Ho 2013", S.BLUE), ("vdb16", "van den Bosch 2016", S.ORANGE),
             ("mm13", "McConnell & Ma 2013", S.AQUA)]
    for key, lab, col in lines:
        f = fits[key]
        _fit_line(ax, f, 40, 420, band=False, color=col, label=f"{lab}: slope {f.b:.2f} ± {f.eb:.2f}")
    ax.legend(loc="upper left", fontsize=7.5)
    _sigma_axis(ax, 40, 420)
    _mass_axis(ax, 1e5, 5e10)
    ax.set_title("Same physics, three slopes", fontsize=10.5)
    ax.text(0.98, 0.03, "each line fitted by us, with that\npaper's own sample and method",
            transform=ax.transAxes, ha="right", fontsize=7.5, color=S.INK_SECONDARY)

    jack = sb["jackknife"]
    jb = np.array([j["b"] for j in jack])
    full = sb["steps"][-1]["b"]
    rng = np.random.default_rng(1)
    axj.axvline(full, color=S.INK_SECONDARY, lw=1.2)
    axj.text(full, 1.12, f"all 45: {full:.2f}", ha="center", fontsize=7.5, color=S.INK_SECONDARY,
             bbox={"facecolor": S.SURFACE, "edgecolor": "none", "pad": 1.5})
    axj.scatter(jb, rng.uniform(-0.3, 0.3, len(jb)), s=20, c=S.BLUE, edgecolors=S.SURFACE, linewidths=0.6, zorder=3)
    for j, ytxt in ((jack[0], -0.85), (jack[-1], 0.75)):
        lab = "M32" if j["name"] == "NGC 221" else j["name"]
        S.annotate(axj, f"without {lab}\n(σ = {j['sigma']:.0f} km/s)", (j["b"], 0.0), (j["b"], ytxt),
                   ha="center", fontsize=7.3)
    axj.set_ylim(-1.25, 1.25)
    axj.set_yticks([])
    axj.spines["left"].set_visible(False)
    axj.set_xlabel("Kormendy & Ho's slope, leaving out one galaxy")
    axj.set_title("A few galaxies at the ends carry the slope", fontsize=10.5)

    steps = sb["steps"]
    yy = np.arange(len(steps))[::-1]
    axb.axvline(4.0, color=S.AXIS, lw=1.0, zorder=0)
    axb.text(4.0, yy[0] + 0.55, "slope 4", ha="center", fontsize=7.3, color=S.MUTED)
    for i, st in enumerate(steps):
        axb.errorbar([st["b"]], [yy[i]], xerr=[st["eb"]], fmt="o", ms=6, color=S.INK, elinewidth=1.0,
                     mec=S.SURFACE, mew=0.8, zorder=3)
        if i:
            d = st["b"] - steps[i - 1]["b"]
            axb.annotate("", xy=(st["b"], yy[i] + 0.15), xytext=(steps[i - 1]["b"], yy[i - 1] - 0.15),
                         arrowprops={"arrowstyle": "-|>", "color": S.MUTED, "lw": 0.9, "shrinkA": 3, "shrinkB": 3})
            axb.text(st["b"] + st["eb"] + 0.06, yy[i], f"{d:+.2f}", fontsize=8, color=S.INK, va="center")
    axb.set_yticks(yy)
    axb.set_yticklabels([f"{st['label']}  (N = {st['n']})" for st in steps], fontsize=8)
    axb.set_xlabel("M-σ slope (least squares with intrinsic scatter, as both papers fit)")
    axb.set_xlim(3.7, 6.2)
    axb.set_ylim(-0.6, len(steps) - 0.25)
    axb.grid(axis="y", visible=False)
    axb.set_title(f"From McConnell & Ma's {steps[0]['b']:.2f} to Kormendy & Ho's {steps[-1]['b']:.2f}, "
                  "one choice at a time", fontsize=10.5)
    S.caption(fig, "Top left: grey = all van den Bosch masses; lines fitted to each paper's own sample. Top right: "
                   "Kormendy & Ho's slope recomputed 45 times, each time without one galaxy. Bottom: start from "
                   "McConnell & Ma's 72 galaxies and adopt Kormendy & Ho's choices in turn; the steps do not commute "
                   "(see README for each change applied alone).", y=0.0)
    return S.save(fig, out / "04_three_slopes")


# ---------------------------------------------------------------------------
# 5. Black hole mass vs bulge mass
# ---------------------------------------------------------------------------

def fig_bulge_mass(res, out):
    cat, fits = res["cat"], res["fits"]
    kh = cat["kh13"].copy()
    kh["Mb"] = 10 ** kh["logMbulge"]
    kh["M"] = 10 ** kh["logM"]
    bm = res["bulge_mass"]
    fit = fits["mbulge"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 5.0), gridspec_kw={"width_ratios": [1.25, 1]})
    xb = np.geomspace(3e8, 3e12, 100)
    for pct in (0.01, 0.1, 1.0, 10.0):
        ax.plot(xb, xb * pct / 100, color=S.GRID, lw=1.0, zorder=0.6)
        xt = 2.2e12 if pct < 5 else 2.5e9
        ax.text(xt, xt * pct / 100 * 1.15, f"{pct:g}%", color=S.MUTED, fontsize=7.5, rotation=31, ha="right")
    y = fit.a + fit.b * (np.log10(xb) - 11)
    ax.fill_between(xb, 10 ** (y - fit.eps), 10 ** (y + fit.eps), color=S.WASH, lw=0, zorder=0.5)
    ax.plot(xb, 10**y, color=S.INK_SECONDARY, lw=2, label=f"our fit: slope {fit.b:.2f} ± {fit.eb:.2f}", zorder=2)
    fitted = kh[kh["kh13_fit"]]
    for kind in ("elliptical", "classical"):
        d = fitted[fitted["bulge"] == kind]
        st = S.BULGE_STYLE[kind]
        ax.scatter(d["Mb"], d["M"], s=30, c=st["color"], marker=st["marker"], edgecolors=S.SURFACE,
                   linewidths=0.8, zorder=3, label=st["label"])
    pb = kh[(kh["bulge"] == "pseudo") & ~kh["upper"]]
    st = S.BULGE_STYLE["pseudo"]
    ax.scatter(pb["Mb"], pb["M"], s=30, c=st["color"], marker=st["marker"], edgecolors=S.SURFACE,
               linewidths=0.8, zorder=3, label=st["label"])
    for name, lab, xy in (("NGC 1277", "NGC 1277", (9e9, 3e10)), ("NGC 4486B", "NGC 4486B", (1.3e9, 2.5e9))):
        r = kh[kh["name"] == name].iloc[0]
        ax.scatter([r["Mb"]], [r["M"]], s=26, marker="x", c=S.MUTED, linewidths=1.0, zorder=3)
        S.annotate(ax, f"{lab}: {bm['monsters_percent'][name]:.0f}%", (r["Mb"], r["M"]), xy, ha="center",
                   fontsize=7.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e9, 3e12)
    ax.set_ylim(3e5, 6e10)
    ax.set_xlabel("bulge stellar mass (M☉)")
    ax.set_ylabel("black hole mass (M☉)")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("A black hole is about half a percent of its bulge")

    ratio = 100 * fitted["M"] / fitted["Mb"]
    pr = 100 * pb["M"] / pb["Mb"]
    ax2.axhline(bm["ratio_at_1e11_percent_ours"], color=S.INK_SECONDARY, lw=1.2, zorder=1)
    ax2.text(2.8e12, bm["ratio_at_1e11_percent_ours"] / 1.25, f"fit at 10¹¹ M☉: {bm['ratio_at_1e11_percent_ours']:.2f}%",
             ha="right", va="top", fontsize=7.5, color=S.INK_SECONDARY)
    for kind in ("elliptical", "classical"):
        d = fitted[fitted["bulge"] == kind]
        st = S.BULGE_STYLE[kind]
        ax2.scatter(d["Mb"], 100 * d["M"] / d["Mb"], s=30, c=st["color"], marker=st["marker"],
                    edgecolors=S.SURFACE, linewidths=0.8, zorder=3)
    ax2.scatter(pb["Mb"], pr, s=30, c=S.AQUA, marker="D", edgecolors=S.SURFACE, linewidths=0.8, zorder=3)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(1e9, 3e12)
    ax2.set_ylim(0.003, 5)
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
    ax2.set_xlabel("bulge stellar mass (M☉)")
    ax2.set_ylabel("black hole mass / bulge mass")
    ax2.set_title("The ratio, galaxy by galaxy", fontsize=10)
    ax2.text(0.03, 0.97, f"ellipticals + classical bulges: median {np.median(ratio):.2f}%\n"
                         f"(68% between {np.percentile(ratio, 16):.2f}% and {np.percentile(ratio, 84):.1f}%)\n"
                         f"pseudobulges: median {np.median(pr):.3f}%",
             transform=ax2.transAxes, va="top", fontsize=7.5, color=S.INK_SECONDARY)
    S.caption(fig, f"Kormendy & Ho (2013) bulge masses (K-band light × mass-to-light ratio). Our fit to their 45 "
                   f"ellipticals and classical bulges gives {bm['ratio_at_1e11_percent_ours']:.2f}% at 10¹¹ M☉ "
                   f"(theirs: {LIT.KH13['ratio_percent']['value']}%); slope {fit.b:.2f} ± {fit.eb:.2f}, scatter "
                   f"{fit.eps:.2f} dex. Grey lines: constant mass ratios. Crosses: the two 'monsters' they set aside.")
    fig.tight_layout()
    return S.save(fig, out / "05_bulge_mass")


# ---------------------------------------------------------------------------
# 6. How small is the sphere of influence?
# ---------------------------------------------------------------------------

def fig_scales(res, out):
    cat = res["cat"]
    sc = cat["vdb16_scales"]
    det = sc[~sc["upper"]].copy()
    det["M"] = 10 ** det["logM"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.8, 5.0))
    m = np.geomspace(1e5, 3e10, 50)
    ax.plot(m, P.schwarzschild_radius_pc(m), color=S.INK, lw=2, label="event horizon, 2GM/c²")
    ax.scatter(det["M"], det["r_infl_pc"], s=16, c=S.BLUE, edgecolors="none", alpha=0.85, zorder=3,
               label="sphere of influence, GM/σ²")
    ax.scatter(det["M"], det["re_pc"], s=16, c=S.ORANGE, marker="s", edgecolors="none", alpha=0.85, zorder=3,
               label="galaxy half-light radius")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e5, 3e10)
    ax.set_ylim(1e-8, 1e6)
    ax.set_xlabel("black hole mass (M☉)")
    ax.set_ylabel("radius (parsecs)")
    for r_pc, lab in ((4.848e-6, "1 AU"), (1.0, "1 pc"), (1e3, "1 kpc")):
        ax.axhline(r_pc, color=S.GRID, lw=1.0, zorder=0.4)
        ax.text(1.15e5, r_pc * 1.5, lab, fontsize=7, color=S.MUTED)
    s = res["scales"]
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("Nine orders of magnitude from horizon to galaxy", fontsize=10.5)
    ax.text(0.03, 0.97, f"median sphere of influence {s['r_infl_pc_median']:.0f} pc\n"
                        f"= {100 * s['r_infl_over_re_median']:.1f}% of the galaxy's radius,\n"
                        f"one ten-millionth of its volume",
            transform=ax.transAxes, va="top", fontsize=7.5, color=S.INK_SECONDARY)

    for meth in ("stars", "gas", "maser", "reverb"):
        d = det[det["method"] == meth]
        st = S.METHOD_STYLE[meth]
        ax2.scatter(d["dist"], d["theta_infl"], s=20, c=st["color"], marker=st["marker"], edgecolors=S.SURFACE,
                    linewidths=0.6, zorder=3, label=st["label"].split(" (")[0])
    th_eff = res["resolution_sim"]["theta_effective"]
    for th, lab in ((C.RESOLUTION_ARCSEC, "Hubble's resolution, 0.1″"),
                    (th_eff, f"effective limit of the sample, {th_eff:.3f}″")):
        ax2.axhline(th, color=S.INK_SECONDARY, lw=1.0, zorder=1)
        ax2.text(5e-3, th * 1.15, lab, ha="left", fontsize=7.3, color=S.INK_SECONDARY)
    mw = res["scales"]["milky_way"]
    m87 = res["scales"]["m87"]
    ax2.plot([0.00793], [mw["theta_arcsec"]], ls="none", zorder=6, **S.ANCHOR_STYLE)
    S.annotate(ax2, "Milky Way", (0.00793, mw["theta_arcsec"]), (0.03, 8), fontsize=7.5)
    ax2.plot([C.M87["D_mpc"]], [m87["theta_arcsec"]], ls="none", zorder=6, **S.ANCHOR_STYLE)
    S.annotate(ax2, "M87", (C.M87["D_mpc"], m87["theta_arcsec"]), (2.0, 20), fontsize=7.5)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(4e-3, 1.3e3)
    ax2.set_ylim(3e-4, 200)
    ax2.set_xlabel("distance (Mpc)")
    ax2.set_ylabel("angular size of the sphere of influence (arcsec)")
    ax2.legend(loc="upper right", fontsize=7.3)
    ax2.set_title("To weigh a black hole by its pull, resolve that sphere", fontsize=10.5)
    S.caption(fig, "van den Bosch (2016) galaxies with measured masses. Masers and light echoes (reverberation) "
                   "do not need to resolve the sphere of influence: masers are mapped by radio interferometry at "
                   "milliarcseconds, and echoes measure size by light-travel time.")
    fig.tight_layout()
    return S.save(fig, out / "06_scales")


# ---------------------------------------------------------------------------
# 7. Selection: are black hole hosts ordinary galaxies?
# ---------------------------------------------------------------------------

def fig_selection(res, out):
    cat = res["cat"]
    t = cat["atlas3d_test"]
    hb = res["host_bias"]
    sim = res["resolution_sim"]
    fig = plt.figure(figsize=(11.6, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 0.95, 1.05], wspace=0.42)
    ax, axh, axs = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    others, hosts = t[~t["host"]], t[t["host"]]
    ax.scatter(10 ** others["logLk"], 10 ** others["logsig"], s=14, c=S.CONTEXT, edgecolors="none",
               label=f"other ATLAS³ᴰ galaxies ({len(others)})")
    ax.scatter(10 ** hosts["logLk"], 10 ** hosts["logsig"], s=28, c=S.BLUE, edgecolors=S.SURFACE, linewidths=0.7,
               zorder=3, label=f"with a black hole mass ({len(hosts)})")
    xl = np.linspace(others["logLk"].min(), others["logLk"].max(), 100)
    ax.plot(10**xl, 10 ** np.polyval(hb["fj_coef"], xl - 11.0), color=S.INK_SECONDARY, lw=2,
            label="typical σ at that luminosity")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("K-band luminosity (L☉), a proxy for stellar mass")
    ax.set_ylabel("σ within the half-light radius (km/s)")
    ax.yaxis.set_major_locator(FixedLocator([50, 70, 100, 150, 200, 300]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("All 260 early-type galaxies within 42 Mpc", fontsize=10.5)

    rng = np.random.default_rng(7)
    for yrow, d, col, lab in ((1, others, S.CONTEXT, "others"), (0, hosts, S.BLUE, "hosts")):
        r = d["dlogsig"].to_numpy()
        axh.scatter(r, yrow + rng.uniform(-0.22, 0.22, len(r)), s=14 if yrow else 22, c=col,
                    edgecolors=S.SURFACE if not yrow else "none", linewidths=0.6, zorder=3)
        axh.plot([r.mean()] * 2, [yrow - 0.3, yrow + 0.3], color=S.INK, lw=2.2, solid_capstyle="round", zorder=4)
        txt = "mean 0 (by construction)" if yrow else f"mean {r.mean():+.3f}"
        axh.text(r.mean() + 0.012, yrow + 0.33, txt, fontsize=7.3, color=S.INK, va="bottom")
    h = hb["hosts"]
    axh.axvline(0, color=S.INK_SECONDARY, lw=1.0, zorder=1)
    axh.set_yticks([1, 0])
    axh.set_yticklabels([f"other\ngalaxies\n({len(others)})", f"black hole\nhosts\n({len(hosts)})"], fontsize=8)
    axh.set_ylim(-0.45, 1.75)
    axh.set_xlim(-0.35, 0.35)
    axh.grid(axis="y", visible=False)
    axh.text(0.5, 0.985, f"hosts: σ higher by {h['mean_dlogsig']:+.3f} ± {h['sem']:.3f} dex\n"
                         f"({100 * (10 ** h['mean_dlogsig'] - 1):.0f}%); a random set of {h['n']}\n"
                         f"does this well {100 * h['perm_p']:.2f}% of the time",
             transform=axh.transAxes, ha="center", va="top", fontsize=7.4, color=S.INK,
             bbox={"facecolor": S.SURFACE, "edgecolor": "none", "pad": 1.5})
    axh.set_xlabel("log σ − typical σ at that luminosity (dex)")
    axh.set_title("Hosts are dynamically hotter", fontsize=10.5)

    rows = [("truth", None), ("detections only", "detections_only"), ("with upper limits", "censored")]
    truths = {k: sim[f"{k}__effective"]["truth"] for k in ("our_fit", "kh13_like")}
    labels = {k: f"true slope {t['b']:.2f},\nscatter {t['eps']:.2f} dex" for k, t in truths.items()}
    for j, (tkey, tlab, col) in enumerate((("our_fit", labels["our_fit"], S.ORANGE),
                                           ("kh13_like", labels["kh13_like"], S.BLUE))):
        for k, th in enumerate(("effective", "strict_0p1")):
            r = sim[f"{tkey}__{th}"]
            y0 = j * 4 + k * 1.5
            vals = [r["truth"]["b"], r["detections_only"]["b"], r["censored"]["b"]]
            errs = [0, r["detections_only"]["b_std"], r["censored"]["b_std"]]
            axs.errorbar(vals[1:], [y0 + 0.35, y0 - 0.35], xerr=errs[1:], fmt="none", ecolor=col, elinewidth=1.0)
            axs.scatter([vals[1]], [y0 + 0.35], s=34, marker="v", c=col, edgecolors=S.SURFACE, zorder=3)
            axs.scatter([vals[2]], [y0 - 0.35], s=34, marker="o", c=col, edgecolors=S.SURFACE, zorder=3)
            axs.plot([vals[0]] * 2, [y0 - 0.7, y0 + 0.7], color=S.INK, lw=1.6)
            lim = f"{r['theta']:.3f}″" if th == "effective" else "0.1″"
            axs.text(3.15, y0, f"{lim}: {100 * r['frac_resolved']:.0f}% resolved", fontsize=7, va="center",
                     color=S.INK_SECONDARY)
        axs.text(3.15, j * 4 + 2.55, tlab, fontsize=7.5, color=S.INK, va="center")
    axs.set_xlim(3.1, 6.2)
    axs.set_ylim(-1.1, 8.3)
    axs.set_yticks([])
    axs.spines["left"].set_visible(False)
    axs.set_xlabel("fitted slope")
    axs.legend(handles=[Line2D([], [], color=S.INK, lw=1.6, label="true slope"),
                        Line2D([], [], ls="none", marker="v", color=S.MUTED, label="fit, resolved masses only"),
                        Line2D([], [], ls="none", marker="o", color=S.MUTED, label="fit, with upper limits")],
               loc="upper right", fontsize=7.2, bbox_to_anchor=(1.0, 1.0))
    axs.set_title("SIMULATION: requiring a resolved sphere", fontsize=10.5)
    S.caption(fig, "Left, middle: ATLAS³ᴰ (Cappellari et al. 2011, 2013): σ_e from integral-field spectra, 2MASS "
                   "luminosity, the same for every galaxy. Right: a labelled simulation on the real σ and distances of "
                   f"the {sim['n_galaxies']} star- and gas-dynamical galaxies; masses drawn from an assumed relation "
                   f"count as measured only if the sphere of influence spans the limit ({sim['realisations']} "
                   "realisations; bars = spread).", y=-0.04)
    return S.save(fig, out / "07_selection")


# ---------------------------------------------------------------------------
# 8. What it implies: feedback
# ---------------------------------------------------------------------------

def fig_feedback(res, out):
    cat, fits = res["cat"], res["fits"]
    fb = res["feedback"]
    kh = _with_linear(cat["kh13"])
    fitted = kh[kh["kh13_fit"]]
    fit = fits["kh13"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.9), gridspec_kw={"width_ratios": [1.3, 1]})
    s = np.geomspace(60, 450, 100)
    lo = P.king_msigma(s, f_gas=0.1)
    hi = P.king_msigma(s, f_gas=0.2, kappa=P.kappa_es(1.0))
    ax.fill_between(s, lo, hi, color=S.BLUE_150, alpha=0.55, lw=0, zorder=0.6)
    ax.plot(s, P.king_msigma(s), color=S.BLUE_550, lw=2, zorder=2,
            label="King (2003): f κ σ⁴ / (π G²): slope 4, no free parameter")
    _fit_line(ax, fit, 60, 450, band=False, label=f"fit to the data: slope {fit.b:.2f} ± {fit.eb:.2f}")
    for kind in ("elliptical", "classical"):
        _points(ax, fitted[fitted["bulge"] == kind], S.BULGE_STYLE[kind], s=26, errors=False)
    ax.legend(loc="upper left", fontsize=7.5)
    _sigma_axis(ax, 60, 450)
    _mass_axis(ax, 1e6, 5e10)
    ax.set_title("A black hole that blows its galaxy's gas away stops growing", fontsize=10.5)
    ax.text(0.98, 0.04, f"at σ = 200 km/s: predicted {fb['king_M_at_200'] / 1e8:.1f}×10⁸ M☉,\n"
                        f"measured {10 ** fits['kh13'].a / 1e8:.1f}×10⁸ M☉ (ellipticals + classical bulges)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=S.INK)

    er = cat["energy_ratio"]
    bins = np.geomspace(150, 8000, 15)
    ax2.hist(er, bins=bins, color=S.BLUE, edgecolor=S.SURFACE, linewidth=1.0)
    ax2.axvline(1.0, color=S.INK, lw=1)
    ax2.axvline(fb["energy_ratio_median"], color=S.INK_SECONDARY, lw=1.2)
    ax2.set_xscale("log")
    ax2.set_xlim(0.5, 2e4)
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}×"))
    ax2.set_xlabel("energy released growing the black hole / bulge binding energy")
    ax2.set_ylabel("galaxies")
    ax2.text(1.15, ax2.get_ylim()[1] * 0.92, "1×: just\nenough", fontsize=7.5, color=S.INK_SECONDARY, va="top")
    ax2.text(fb["energy_ratio_median"] * 1.1, ax2.get_ylim()[1] * 0.92,
             f"median {fb['energy_ratio_median']:.0f}×", fontsize=7.5, color=S.INK_SECONDARY, va="top")
    ax2.set_title("More than enough energy to do it", fontsize=10.5)
    S.caption(fig, f"Left: Kormendy & Ho's 45 ellipticals and classical bulges. King's line uses the cosmic baryon "
                   f"fraction f = {C.F_GAS_COSMIC:.3f} (Planck 2018) and electron scattering κ = {fb['kappa']:.3f} cm²/g; "
                   "the band spans f = 0.1-0.2 and pure-hydrogen κ. Right: 0.1 M c² (accretion's radiative output) "
                   "over M_bulge σ², galaxy by galaxy.")
    fig.tight_layout()
    return S.save(fig, out / "08_feedback")


# ---------------------------------------------------------------------------
# 9. Active galaxies and the virial factor
# ---------------------------------------------------------------------------

def fig_agn(res, out):
    cat, fits = res["cat"], res["fits"]
    m = cat["agn_matched"].copy()
    agn = res["agn"]
    m["sigma"] = 10 ** m["logsig"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.9), gridspec_kw={"width_ratios": [1.2, 1]})
    q = _with_linear(cat["vdb16"])
    q = q[(q["selected"] == 1) & (q["method"] != "reverb") & ~q["upper"]]
    ax.scatter(q["sigma"], q["M"], s=12, c=S.CONTEXT, edgecolors="none", label="quiet galaxies (stars, gas, masers)")
    qfit = fits["vdb16_no_reverb"]
    _fit_line(ax, qfit, 40, 420, band=False, label="their fit")
    cal = agn["calibrations"][0]
    vp = 10 ** m["logVP"]
    ax.scatter(m["sigma"], vp, s=30, marker="^", facecolors="none", edgecolors=S.INK_SECONDARY, linewidths=1.0,
               zorder=3, label="AGN: echo size × speed² / G (the 'virial product')")
    ax.scatter(m["sigma"], vp * cal["f"], s=30, marker="^", c=S.ORANGE, edgecolors=S.SURFACE, linewidths=0.6,
               zorder=4, label=f"× f = {cal['f']:.1f}: what puts them on the line")
    for _, r in m.iterrows():
        ax.plot([r["sigma"]] * 2, [10 ** r["logVP"], 10 ** r["logVP"] * cal["f"]], color=S.GRID, lw=0.8, zorder=1)
    ax.legend(loc="upper left", fontsize=7.3)
    _sigma_axis(ax, 40, 420)
    _mass_axis(ax, 1e4, 3e10)
    ax.set_title("Light echoes need one number from M-σ", fontsize=10.5)

    rows = []
    short = ["ours, vs all quiet galaxies", "ours, vs ellipticals + classical", "ours, vs pseudobulges"]
    for lab, c in zip(short, agn["calibrations"]):
        rows.append((lab, c["f"], c["f_lo"], c["f_hi"], S.ORANGE))
    rows.append(("Ho & Kim 2014: classical", LIT.HO_KIM_2014["classical"][0],
                 LIT.HO_KIM_2014["classical"][0] - LIT.HO_KIM_2014["classical"][1],
                 LIT.HO_KIM_2014["classical"][0] + LIT.HO_KIM_2014["classical"][1], S.INK_SECONDARY))
    rows.append(("Ho & Kim 2014: pseudobulges", LIT.HO_KIM_2014["pseudo"][0],
                 LIT.HO_KIM_2014["pseudo"][0] - LIT.HO_KIM_2014["pseudo"][1],
                 LIT.HO_KIM_2014["pseudo"][0] + LIT.HO_KIM_2014["pseudo"][1], S.INK_SECONDARY))
    for k, (f, ef) in LIT.F_FACTOR.items():
        rows.append((k.replace("+", " et al. 20"), f, f - ef, f + ef, S.INK_SECONDARY))
    yy = np.arange(len(rows))[::-1]
    for y, (lab, f, lo_, hi_, col) in zip(yy, rows):
        ax2.errorbar([f], [y], xerr=[[f - lo_], [hi_ - f]], fmt="o", color=col, ms=5, elinewidth=1.0,
                     mec=S.SURFACE)
    ax2.set_yticks(yy)
    ax2.set_yticklabels([r[0] for r in rows], fontsize=7.6)
    ax2.set_xscale("log")
    ax2.set_xlim(1, 15)
    ax2.xaxis.set_major_locator(FixedLocator([1, 2, 3, 5, 10]))
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax2.xaxis.set_minor_locator(NullLocator())
    ax2.set_xlabel("virial factor ⟨f⟩")
    ax2.grid(axis="y", visible=False)
    ax2.axhline(len(rows) - 3.5, color=S.GRID, lw=1)
    ax2.set_title("⟨f⟩ depends on which galaxies you match", fontsize=10.5)
    S.caption(fig, f"{agn['n_agn']} reverberation-mapped AGN with σ in van den Bosch (2016) and virial products from "
                   "the AGN Black Hole Mass Database (Bentz & Katz 2015, exported with f = 1). Our ⟨f⟩: offset of the "
                   "AGN from each quiet-galaxy relation (fitted without any reverberation masses), ± its uncertainty.")
    fig.tight_layout()
    return S.save(fig, out / "09_agn_virial_factor")


# ---------------------------------------------------------------------------
# 10. Same galaxy, different numbers
# ---------------------------------------------------------------------------

def fig_systematics(res, out):
    cat = res["cat"]
    m = cat["compilations_matched"]
    sy = res["systematics"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.8), gridspec_kw={"width_ratios": [1, 1.1]})
    lim = (50, 450)
    ax.plot(lim, lim, color=S.AXIS, lw=1, zorder=0.5)
    for kind in ("elliptical", "classical", "pseudo"):
        d = m[m["bulge"] == kind]
        st = S.BULGE_STYLE[kind]
        ax.scatter(10 ** d["logsig_kh"], 10 ** d["logsig_mm"], s=26, c=st["color"], marker=st["marker"],
                   edgecolors=S.SURFACE, linewidths=0.7, zorder=3, label=st["label"])
    for name, xy in (("Circinus", (110, 300)), ("NGC 1300", (60, 260))):
        r = m[m["name"] == name].iloc[0]
        S.annotate(ax, name, (10 ** r["logsig_kh"], 10 ** r["logsig_mm"]), xy, fontsize=7.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    for a in (ax.xaxis, ax.yaxis):
        a.set_major_locator(FixedLocator([50, 70, 100, 150, 200, 300, 400]))
        a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        a.set_minor_locator(NullLocator())
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("σ in Kormendy & Ho 2013 (km/s)")
    ax.set_ylabel("σ in McConnell & Ma 2013 (km/s)")
    ax.legend(loc="lower right", fontsize=7.5)
    ss = sy["sigma_mm_vs_kh"]
    ax.set_title("One galaxy, two velocity dispersions", fontsize=10.5)
    ax.text(0.03, 0.97, f"rms difference {ss['rms']:.3f} dex; for pseudobulges\nσ itself is ill-defined",
            transform=ax.transAxes, va="top", fontsize=7.5, color=S.INK_SECONDARY)

    hist = LIT.M87_HISTORY
    yy = np.arange(len(hist))[::-1]
    for y, h in zip(yy, hist):
        scale = C.M87["D_mpc"] / h["D"]  # dynamical and shadow masses both scale with distance
        ax2.errorbar([h["M"] * scale / 1e9], [y], xerr=[[h["lo"] * scale / 1e9], [h["hi"] * scale / 1e9]],
                     fmt="o", ms=6, color=S.INK if "EHT" in h["ref"] else S.BLUE, elinewidth=1.2, mec=S.SURFACE)
    kh_pred = res["anchors"]["m87"]["kh13"]["predicted_M"] / 1e9
    ax2.axvspan(kh_pred / 10 ** res["fits"]["kh13"].eps, kh_pred * 10 ** res["fits"]["kh13"].eps,
                color=S.WASH, lw=0, zorder=0)
    ax2.axvline(kh_pred, color=S.INK_SECONDARY, lw=1.2)
    ax2.text(kh_pred * 1.04, len(hist) - 0.45, "M-σ prediction\n(classical relation)", fontsize=7.3,
             color=S.INK_SECONDARY, va="top")
    ax2.set_xlim(0, 9)
    ax2.set_ylim(-0.7, len(hist) - 0.3)
    ax2.set_yticks(yy)
    ax2.set_yticklabels([f"{h['ref']}\n{h['method']}" for h in hist], fontsize=7.6)
    ax2.grid(axis="y", visible=False)
    ax2.set_xlabel("mass of M87's black hole at 16.8 Mpc (billion M☉)")
    ax2.set_title("Four weighings of M87*", fontsize=10.5)
    S.caption(fig, "Left: the 72 galaxies in both reviews. Right: published masses rescaled to the EHT's distance "
                   "(16.8 Mpc); dynamical masses and the shadow-based mass both scale in proportion to distance. "
                   "Shaded: the classical relation's intrinsic scatter at M87's σ = 324 km/s.")
    fig.tight_layout()
    return S.save(fig, out / "10_systematics")


def make_all(res) -> list:
    out = C.FIGURE_DIR
    paths = []
    for f in (fig_m_sigma, fig_methods, fig_bulges, fig_slopes, fig_bulge_mass, fig_scales, fig_selection,
              fig_feedback, fig_agn, fig_systematics):
        paths += f(res, out)
    return paths
