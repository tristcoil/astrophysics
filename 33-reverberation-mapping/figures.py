"""Every figure, drawn from the cached analysis results.

    python run.py --figures-only
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import config as C
import literature as LIT
import plotstyle as S

S.use_style()
FIG = C.FIGURE_DIR
Z1 = 1.0 + C.Z_NGC5548


def _std(y):
    y = np.asarray(y, float)
    return (y - y.mean()) / y.std()


def _label(ax, x, y, text, color=S.INK, dx=4, dy=0, ha="left", va="center", **kw):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", color=color, fontsize=8, ha=ha,
                va=va, **kw)


def _pm(v, nd=2):
    return f"{v[0]:.{nd}f} (−{v[1]:.{nd}f}/+{v[2]:.{nd}f})"


# ---------------------------------------------------------------------------
# 1. The keeper: the echo of NGC 5548 in 2014
# ---------------------------------------------------------------------------


def fig_echo(r: dict) -> list:
    lcs = r["storm"]["lcs"]
    c, h = lcs["F5100"], lcs["Hb"]
    res = r["storm"]["raw"][("Hb", "F5100", "full")]
    lag = r["storm"]["lags"][("Hb", "F5100", "full")]["cent"]
    lag_obs = lag[0] * Z1
    fig = plt.figure(figsize=(11.5, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.5, 1], hspace=0.42, wspace=0.18)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1])
    import iccf

    yc, yh = _std(iccf.detrend(c.t, c.y)), _std(iccf.detrend(h.t, h.y))
    for ax in (ax1, ax2):
        ax.axvspan(C.STORM_SPLIT_THJD, c.t.max() + 2, color=S.EXCLUDED, zorder=0, lw=0)
    ax1.plot(c.t, yc, "-", color=S.C_CONT, lw=1.3, label="5100 Å continuum (the accretion disk's light)")
    ax1.plot(h.t, yh, "-", color=S.C_LINE, lw=1.3, label="broad Hβ line (gas around the black hole)")
    ax1.set_title("1. As observed: the line follows the continuum, a few days late")
    ax1.legend(loc="upper left", ncol=2, fontsize=7.5)
    ax2.plot(c.t, yc, "-", color=S.C_CONT, lw=1.3)
    ax2.plot(h.t - lag_obs, yh, "-", color=S.C_LINE, lw=1.3)
    ax2.set_title(f"2. The line moved back by the lag ({lag_obs:.1f} d): echo and source line up")
    ax2.set_xlabel("Date in 2014 (HJD − 2,450,000)")
    for ax in (ax1, ax2):
        ax.set_ylabel("flux (standardized)")
        ax.set_ylim(-2.9, 3.2)
    ax1.text(C.STORM_SPLIT_THJD + 2, -2.65, "second half: the lines decouple\n(the 'BLR holiday', figure 2)",
             fontsize=7, color=S.INK_SECONDARY, va="bottom")
    # ICCF
    lags_rest = res.lags / Z1
    ax3.plot(lags_rest, res.r, color=S.C_MODEL, lw=1.4)
    ax3.axvline(0, color=S.AXIS, lw=0.8)
    ax3.axhline(0.8 * res.observed.r_max, color=S.MUTED, lw=0.6)
    ax3.text(38, 0.8 * res.observed.r_max - 0.03, "0.8 r_max", fontsize=6.5, color=S.MUTED, va="top", ha="right")
    ax3.set_xlim(-20, 39)
    ax3.set_ylim(-0.85, 1.05)
    ax3.set_xlabel("trial lag τ (days, rest frame)")
    ax3.set_ylabel("correlation r")
    ax3.set_title("3. Slide and correlate (ICCF)")
    ax3.annotate(f"peak r = {res.observed.r_max:.2f}", (res.observed.tau_peak / Z1, res.observed.r_max),
                 xytext=(28, 4), textcoords="offset points", fontsize=7.5, color=S.INK_SECONDARY,
                 arrowprops=dict(arrowstyle="-", color=S.INK_SECONDARY, lw=0.7))
    # CCCD
    bins = np.arange(2.0, 6.75, 0.125)
    ax4.hist(res.cccd / Z1, bins=bins, color=S.C_LINE, alpha=0.85, edgecolor=S.SURFACE, linewidth=0.6)
    pub = LIT.PEI17_LAGS[("Hb", "F5100", "cent")]
    ymax = ax4.get_ylim()[1]
    ax4.errorbar([pub[0]], [ymax * 1.08], xerr=[[pub[1]], [pub[2]]], fmt="o", color=S.C_LIT, ms=4, capsize=2)
    ax4.errorbar([lag[0]], [ymax * 1.18], xerr=[[lag[1]], [lag[2]]], fmt="o", color=S.INK, ms=4, capsize=2)
    _label(ax4, pub[0] + pub[2], ymax * 1.08, f"Pei et al. 2017: {pub[0]:.2f} ± {pub[1]:.2f}", S.INK_SECONDARY)
    _label(ax4, lag[0] + lag[2], ymax * 1.18, f"ours: {_pm(lag)} d", S.INK)
    ax4.set_ylim(0, ymax * 1.3)
    ax4.set_xlabel("centroid lag (days, rest frame)")
    ax4.set_ylabel(f"of {res.n_mc} resamplings")
    ax4.set_title("4. Lag and its error (FR/RSS)")
    R_au = r["storm_mass"]["R_au"]
    fig.suptitle(f"NGC 5548 in 2014: the broad-line gas echoes the disk {lag[0]:.1f} days later, "
                 f"so it sits {lag[0]:.1f} light-days (~{R_au:,.0f} AU) from the black hole", x=0.06, ha="left",
                 y=0.99)
    S.caption(fig, "Data: AGN STORM ground-based spectra (Pei et al. 2017), MDM Observatory, 133 nights, "
                   "4 Jan – 6 Jul 2014. Both light curves linearly detrended and scaled to zero mean, unit variance "
                   "(H-beta varies by ±4%, the continuum by ±7%). Lag: interpolated cross-correlation centroid, "
                   "1000 flux-randomization / random-subset resamplings.")
    return S.save(fig, FIG / "01_the_echo")


# ---------------------------------------------------------------------------
# 2. The BLR holiday
# ---------------------------------------------------------------------------


def fig_holiday(r: dict) -> list:
    h = r["holiday"]
    order = [("CIV", "C IV 1549 Å"), ("SiIV", "Si IV + O IV] 1400 Å"), ("HeII", "He II + O III] 1640 Å"),
             ("Lya", "Lyα 1216 Å"), ("Hb", "Hβ 4861 Å (ground)")]
    fig = plt.figure(figsize=(11.5, 9.0))
    gs = fig.add_gridspec(len(order), 2, width_ratios=[3.2, 1], hspace=0.35, wspace=0.12)
    seg = {k: (v[0] - 50000, v[1] - 50000) for k, v in LIT.GOAD16_SEGMENTS.items()}
    axes = []
    for i, (key, name) in enumerate(order):
        ax = fig.add_subplot(gs[i, 0], sharex=axes[0] if axes else None)
        axes.append(ax)
        d = h[key]
        t, fl, fp = d["t_shift"], d["f_line"], d["f_pred"]
        use = d["usable"]
        m = np.mean(fp[use])
        if key == "Hb":
            b = C.HB_SEGMENTS_THJD
            lo, hi = b[1], b[2]
        else:
            lo, hi = seg["red"]
        ax.axvspan(lo, hi, color=S.HOLIDAY_WASH, lw=0, zorder=0)
        o = np.argsort(t)
        ax.plot(t[o][use[o]], 100 * fp[o][use[o]] / m, "-", color=S.C_MODEL, lw=1.2, zorder=2)
        ax.plot(t[use], 100 * fl[use] / m, "o", color=S.C_LINE, ms=2.8, mec=S.SURFACE, mew=0.3, zorder=3)
        if (~use).any():
            ax.plot(t[~use], 100 * fl[~use] / m, "o", mfc=S.SURFACE, mec=S.C_LINE, ms=2.8, mew=0.6, zorder=3)
        ax.set_ylabel("% of mean", fontsize=7.5)
        pub = LIT.GOAD16_FLOST.get(key, LIT.PEI17_HB_FLOST)
        ax.set_title(name, loc="left", fontsize=9, pad=4)
        ax.set_title(f"lost {100 * d['flost_red']:.0f}% in the shaded weeks (published {100 * pub:.0f}%)",
                     loc="right", fontsize=8, fontweight="normal", color=S.INK_SECONDARY, pad=4)
        if i < len(order) - 1:
            plt.setp(ax.get_xticklabels(), visible=False)
    axes[-1].set_xlabel("Date of the continuum driving the line (HJD − 2,450,000; line moved back by its lag)")
    handles = [Line2D([], [], color=S.C_MODEL, lw=1.2, label="predicted: F_line ∝ F_cont^η, η fitted on normal epochs"),
               Line2D([], [], marker="o", ls="none", color=S.C_LINE, ms=4, label="observed broad-line flux"),
               plt.Rectangle((0, 0), 1, 1, color=S.HOLIDAY_WASH, label="the anomaly (Goad et al. 2016; Pei et al. 2017)")]
    fig.legend(handles=handles, loc="upper left", fontsize=7.5, ncol=3, bbox_to_anchor=(0.06, 0.965), frameon=False)
    # Summary: fraction lost, ours vs published
    axs = fig.add_subplot(gs[1:4, 1])
    ys = np.arange(len(order))[::-1]
    for y, (key, name) in zip(ys, order):
        d = h[key]
        pub = LIT.GOAD16_FLOST.get(key, LIT.PEI17_HB_FLOST)
        axs.errorbar(100 * d["flost_red"], y, xerr=100 * 2 * d["flost_red_err"], fmt="o", color=S.C_LINE, ms=5,
                     capsize=2, zorder=3)
        axs.plot(100 * pub, y, "|", color=S.INK, ms=12, mew=1.5, zorder=4)
    axs.set_yticks(ys)
    axs.set_yticklabels(["C IV", "Si IV", "He II", "Lyα", "Hβ"], fontsize=8)
    axs.set_xlim(0, 30)
    axs.set_xlabel("line flux missing during the anomaly (%)")
    axs.set_title("Missing flux", fontsize=10)
    axs.legend(handles=[Line2D([], [], marker="o", ls="none", color=S.C_LINE, label="ours (±2σ)"),
                        Line2D([], [], marker="|", ls="none", color=S.INK, ms=10, mew=1.5, label="published")],
               loc="lower right", fontsize=7)
    axs.grid(axis="y", visible=False)
    fig.suptitle("The 'BLR holiday' of 2014: for about five weeks the broad lines stopped echoing the continuum",
                 x=0.06, ha="left", y=1.0)
    S.caption(fig, "HST/COS far-UV lines (De Rosa et al. 2015) against the 1158 Å continuum (Fausnaugh et al. 2016); "
                   "H-beta from MDM (Pei et al. 2017; open circles: last 10 nights, calibration unreliable). "
                   "Narrow-line flux removed; each line moved back by its lag over the first 75 days. "
                   "High-ionization lines (C IV, Si IV, He II) lost the most: the gas saw fewer of the most energetic photons.")
    return S.save(fig, FIG / "02_blr_holiday")


# ---------------------------------------------------------------------------
# 3. Twenty years of NGC 5548: a breathing broad-line region
# ---------------------------------------------------------------------------


def _jd_to_year(jd):
    """JD - 2,400,000 (or MJD) to decimal year."""
    return 1858.877 + (np.asarray(jd, float) - 0.5) / 365.25


def fig_watch(r: dict) -> list:
    w = r["watch"]
    lcs, seasons = w["lcs"], w["seasons"]
    br = r["breathing"]
    fig = plt.figure(figsize=(11.5, 9.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.55], hspace=0.38, wspace=0.22)
    a1 = fig.add_subplot(gs[0, :])
    a2 = fig.add_subplot(gs[1, :], sharex=a1)
    for k, s in enumerate(seasons):
        y0, y1 = _jd_to_year(s["jd_start"]), _jd_to_year(s["jd_end"])
        for ax in (a1, a2):
            ax.axvspan(y0, y1, color=S.EXCLUDED, lw=0, zorder=0)
        a2.text(0.5 * (y0 + y1), 1.02, f"{s['tau'][0]:.0f}", transform=a2.get_xaxis_transform(), ha="center",
                fontsize=7, color=S.INK_SECONDARY)
    c, h = lcs["F5100"], lcs["Hb"]
    a1.plot(_jd_to_year(c.t), c.y, ".", color=S.C_CONT, ms=2.5)
    a2.plot(_jd_to_year(h.t), h.y, ".", color=S.C_LINE, ms=2.5)
    a1.set_ylabel("AGN 5100 Å\n(1e-15 erg/s/cm²/Å)")
    a2.set_ylabel("broad Hβ\n(1e-13 erg/s/cm²)")
    a1.set_title("NGC 5548 for 20 years: continuum (host starlight removed) and the broad Hβ echo", loc="left")
    a2.set_title("lag of each monitoring season (rest-frame days, ours) ↓", loc="right", fontsize=8,
                 fontweight="normal", color=S.INK_SECONDARY, pad=14)
    plt.setp(a1.get_xticklabels(), visible=False)
    a2.set_xlim(1988.5, 2008.8)
    a2.set_xticks([1990, 1995, 2000, 2005])
    a2.set_xticklabels(["1990", "1995", "2000", "2005"])
    # Our lags vs published
    a3 = fig.add_subplot(gs[2, 0])
    pub = np.array([s["pub_tau"] for s in seasons])
    ours = np.array([s["tau"] for s in seasons])
    lamp = np.array([s["jd_start"] > 54500 for s in seasons])
    a3.plot([0, 32], [0, 32], color=S.AXIS, lw=0.8, zorder=0)
    a3.errorbar(pub[~lamp, 0], ours[~lamp, 0], xerr=pub[~lamp, 1:].T, yerr=ours[~lamp, 1:].T, fmt="o",
                color=S.C_N5548, ms=4.5, mec=S.SURFACE, mew=0.6, elinewidth=0.8, capsize=0)
    a3.errorbar(pub[lamp, 0], ours[lamp, 0], xerr=pub[lamp, 1:].T, yerr=ours[lamp, 1:].T, fmt="o", mfc=S.SURFACE,
                color=S.C_N5548, ms=4.5, mew=1.0, elinewidth=0.8, capsize=0)
    S.annotate(a3, "2008 (LAMP): the published lag\nused B/V photometry as driver", (pub[lamp, 0][0], ours[lamp, 0][0]),
               (11.0, 0.5), color=S.INK_SECONDARY)
    a3.set_xlim(0, 32)
    a3.set_ylim(-2, 32)
    a3.set_xlabel("published H-beta lag (days; Peterson et al. 2002, 2004; Bentz et al. 2007; Denney et al. 2010)",
                  fontsize=7.5)
    a3.set_ylabel("our lag from the same light curves (days)")
    a3.set_title("Our cross-correlation reproduces 16 published seasons", loc="left", fontsize=10)
    # Breathing: lag vs luminosity
    a4 = fig.add_subplot(gs[2, 1])
    d = br["sample"]
    fit = br["bentz07_no_y12"]
    xs = np.linspace(42.4, 43.6, 50)
    a4.plot(xs, 10 ** (fit["a"][0] + fit["b"][0] * (xs - 44.0)), color=S.C_MODEL, lw=1.3, zorder=1)
    loc = r.get("local_rl")
    if loc is not None:
        g = loc["all"]
        a4.plot(xs, 10 ** (g["a"][0] + g["b"][0] * (xs - 44.0)), color=S.C_LIT, lw=1.0, zorder=1)
    late = d.jd > 53500
    a4.errorbar(d.x[~late] + 44.0, d.tau[~late], yerr=[d.lo[~late], d.hi[~late]], xerr=d.sx[~late], fmt="o",
                color=S.C_N5548, ms=4.5, mec=S.SURFACE, mew=0.6, elinewidth=0.8, capsize=0, zorder=3)
    a4.errorbar(d.x[late] + 44.0, d.tau[late].clip(lower=0.3), yerr=[np.minimum(d.lo[late], d.tau[late] - 0.25),
                                                                     d.hi[late]],
                fmt="o", mfc=S.SURFACE, color=S.C_N5548, ms=4.5, mew=1.0, elinewidth=0.8, capsize=0, zorder=3)
    # 2012 (De Rosa et al. 2018) and 2014 (ours)
    lum, sm = r["storm_lum"], r["storm_mass"]
    a4.errorbar([lum["logL"]], [sm["tau"]], yerr=[[sm["tau_err"]], [sm["tau_err"]]], xerr=[lum["logL_err"]],
                fmt="D", color=S.C_LINE, ms=6, mec=S.SURFACE, mew=0.6, elinewidth=0.9, capsize=0, zorder=4)
    _label(a4, lum["logL"], sm["tau"], "2014 (AGN STORM)", S.INK, dx=8, dy=-2)
    a4.errorbar([43.07], [2.83], yerr=[[0.96], [0.88]], xerr=[0.06], fmt="s", mfc=S.SURFACE, color=S.C_LINE, ms=5,
                mew=1.0, elinewidth=0.8, capsize=0, zorder=4)
    _label(a4, 43.07, 2.83, "2012", S.INK_SECONDARY, dx=7, dy=-3)
    pred = 10 ** (fit["a"][0] + fit["b"][0] * (lum["logL"] - 44.0))
    a4.annotate("", (lum["logL"], sm["tau"] * 1.08), (lum["logL"], pred), arrowprops=dict(arrowstyle="->",
                color=S.C_LINE, lw=1.0))
    _label(a4, lum["logL"], np.sqrt(pred * sm["tau"]), f"{pred / sm['tau']:.1f}× shorter\nthan its own trend",
           S.C_LINE, dx=-6, ha="right")
    a4.set_yscale("log")
    a4.set_ylim(0.8, 60)
    a4.set_yticks([1, 2, 5, 10, 20, 50])
    a4.set_yticklabels(["1", "2", "5", "10", "20", "50"])
    a4.set_xlim(42.45, 43.62)
    a4.set_xlabel("log AGN luminosity λL(5100 Å) (erg/s)")
    a4.set_ylabel("H-beta lag = BLR radius (light-days)")
    a4.set_title(f"Brighter season, bigger BLR: R ∝ L^{fit['b'][0]:.2f}±{0.5 * (fit['b'][1] + fit['b'][2]):.2f}",
                 loc="left", fontsize=10)
    hdl = [Line2D([], [], marker="o", ls="none", color=S.C_N5548, label="AGN Watch 1989–2001, 2005 (fitted)"),
           Line2D([], [], marker="o", ls="none", mfc=S.SURFACE, color=S.C_N5548, label="2007, 2008"),
           Line2D([], [], color=S.C_MODEL, label="NGC 5548's own trend (Year 12 excluded)"),
           Line2D([], [], color=S.C_LIT, lw=1.0, label="all local AGN (figure 5), slope 0.53")]
    a4.legend(handles=hdl, loc="upper left", fontsize=7)
    fig.suptitle("One black hole, 25 years: NGC 5548's broad-line region grows and shrinks with the light that "
                 "powers it", x=0.06, ha="left", y=0.995)
    S.caption(fig, "Light curves: Peterson et al. 2013 (AGN Watch and later campaigns, recalibrated for [O III] "
                   "variability). Lags: our ICCF centroids, no detrending (as published), 1000 FR/RSS resamplings, "
                   "rest frame. Luminosity: mean host-subtracted flux of each season at 75 Mpc. 2012: De Rosa et "
                   "al. 2018 (from the AGN Black Hole Mass Database).")
    return S.save(fig, FIG / "03_twenty_years")


# ---------------------------------------------------------------------------
# 4. The virial test: gravity sets the line widths
# ---------------------------------------------------------------------------

LINE_GROUPS = {
    "Hb": ("Hβ", S.C_LOCAL),
    "HeII": ("He II 4686, Hα, Hγ", S.C_N5548), "Ha": ("He II 4686, Hα, Hγ", S.C_N5548),
    "Hg": ("He II 4686, Hα, Hγ", S.C_N5548),
    "CIV": ("UV lines", S.C_LINE), "SiIV": ("UV lines", S.C_LINE), "HeII1640": ("UV lines", S.C_LINE),
    "CIII]": ("UV lines", S.C_LINE), "Lya": ("UV lines", S.C_LINE),
}


def fig_virial(r: dict) -> list:
    v = r["virial"]
    d = v["sample"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.9), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.25))
    taus = np.logspace(np.log10(1.0), np.log10(40.0), 50)
    vp = 10 ** v["logvp_weighted"]
    sig = np.sqrt(vp * C.MSUN_G * C.G_CGS / (C.C_CGS * taus * C.DAY_S)) / 1e5
    a1.plot(taus, sig, color=S.C_MODEL, lw=1.3, zorder=1)
    i15 = np.argmin(np.abs(taus - 2.4))
    _label(a1, taus[i15], sig[i15], "pure gravity: σ ∝ τ^−1/2", S.C_MODEL, dx=8, dy=4, ha="left", va="bottom")
    f = v["fit_all"]
    a1.plot(taus, 10 ** (f["a"][0] + f["b"][0] * np.log10(taus)), color=S.C_LIT, lw=1.0, zorder=1)
    _label(a1, 30, 10 ** (f["a"][0] + f["b"][0] * np.log10(30)), "best fit", S.C_LIT, dx=0, dy=8, ha="left",
           va="bottom")
    seen = set()
    for _, row in d.iterrows():
        name, col = LINE_GROUPS.get(row.line, (row.line, S.C_REST))
        mk = "D" if (row.line == "Hb" and row.jd_start > 56600) else "o"
        lab = None if name in seen else name
        seen.add(name)
        a1.errorbar(row.tau_cent, row.sigma_line, xerr=[[row.tau_cent_lo], [row.tau_cent_hi]],
                    yerr=row.sigma_line_err, fmt=mk, color=col, ms=5 if mk == "o" else 6, mec=S.SURFACE, mew=0.6,
                    elinewidth=0.7, capsize=0, label=lab, zorder=3)
    hb14 = d[(d.line == "Hb") & (d.jd_start > 56600)].iloc[0]
    _label(a1, hb14.tau_cent, hb14.sigma_line, "Hβ in 2014\n(our lag)", S.INK, dx=8, dy=6)
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlim(1.0, 40)
    a1.set_ylim(1300, 6500)
    a1.set_xticks([1, 2, 5, 10, 20, 40])
    a1.set_xticklabels(["1", "2", "5", "10", "20", "40"])
    a1.set_yticks([1500, 2000, 3000, 4000, 5000, 6000])
    a1.set_yticklabels(["1500", "2000", "3000", "4000", "5000", "6000"])
    a1.minorticks_off()
    a1.set_xlabel("lag τ = distance from the black hole (light-days)")
    a1.set_ylabel("line width σ_line, rms spectrum (km/s)")
    a1.set_title(f"Closer gas moves faster (fitted slope {f['b'][0]:.2f} ± {0.5 * (f['b'][1] + f['b'][2]):.2f})",
                 loc="left", fontsize=10)
    a1.legend(loc="lower left", fontsize=7.5)
    # Implied masses
    fdyn = r["f_dyn"]["mean"][0] if "f_dyn" in r else C.LOG_F_SHEN24
    m = np.log10(d.vp) + fdyn
    me = 0.4343 * d.vp_err / d.vp
    yrs = _jd_to_year(d.jd_start)
    for (_, row), mm, ee, yy in zip(d.iterrows(), m, me, yrs):
        name, col = LINE_GROUPS.get(row.line, (row.line, S.C_REST))
        mk = "D" if (row.line == "Hb" and row.jd_start > 56600) else "o"
        a2.errorbar(yy, mm, yerr=ee, fmt=mk, color=col, ms=5 if mk == "o" else 6, mec=S.SURFACE, mew=0.6,
                    elinewidth=0.7, capsize=0, zorder=3)
    a2.axhline(v["logvp_weighted"] + fdyn, color=S.C_MODEL, lw=1.2, zorder=1)
    S.annotate(a2, f"weighted mean {10 ** (v['logvp_weighted'] + fdyn) / 1e7:.1f} × 10⁷ M☉",
               (2003.5, v["logvp_weighted"] + fdyn), (1997.5, 8.48), color=S.C_MODEL)
    a2.set_xlabel("year of the campaign")
    a2.set_ylabel("log black hole mass (M☉), f = 4.2")
    a2.set_xlim(1987, 2016)
    a2.set_ylim(6.9, 8.7)
    a2.set_title(f"…so every line and year gives one mass ({v['logvp_scatter']:.2f} dex scatter)",
                 loc="left", fontsize=10)
    fig.suptitle("The virial test: NGC 5548's broad lines orbit one point mass", x=0.06, ha="left", y=1.01)
    S.caption(fig, f"{len(d)} lag + line-width pairs for NGC 5548 from the AGN Black Hole Mass Database (Bentz & Katz "
                   "2015): UV lines from IUE (1989) and HST (1993), optical lines 1988–2014; the 2014 H-beta lag is "
                   f"ours. Over a factor of 10 in lag the fitted widths change by {10 ** -f['b'][0]:.1f}× "
                   f"(gravity alone: 3.2×; slope {f['b'][0]:.2f} vs −1/2). Mass = f c τ σ² / G with log f = "
                   "0.63 (dynamical models, reproduced from Shen et al. 2024).")
    return S.save(fig, FIG / "04_virial_test")


# ---------------------------------------------------------------------------
# 5. The radius-luminosity relation
# ---------------------------------------------------------------------------


def fig_rl(r: dict) -> list:
    loc = r["local_rl"]
    sd = r["sdssrm_rl"]["HB"]
    ds, dd = loc["sample"], sd["sample"]
    fl = loc["all"]
    fs = sd["fit"]
    fig = plt.figure(figsize=(11.5, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.4, 1], width_ratios=[3.2, 1], hspace=0.08, wspace=0.05)
    a1 = fig.add_subplot(gs[0, 0])
    a2 = fig.add_subplot(gs[1, 0], sharex=a1)
    a3 = fig.add_subplot(gs[1, 1], sharey=a2)
    xs = np.linspace(41.6, 45.6, 60)
    b13 = LIT.BENTZ13_RL
    a1.plot(xs, 10 ** (fl["a"][0] + fl["b"][0] * (xs - 44)), color=S.C_MODEL, lw=1.5, zorder=2,
            label=f"our fit, local AGN: slope {fl['b'][0]:.2f} ± {0.5 * (fl['b'][1] + fl['b'][2]):.2f}, "
                  f"scatter {fl['sig'][0]:.2f} dex")
    a1.plot(xs, 10 ** (b13["K"][0] + b13["alpha"][0] * (xs - 44)), color=S.C_LIT, lw=1.0, zorder=1,
            label=f"Bentz et al. 2013: slope {b13['alpha'][0]:.3f}, scatter {b13['scatter_dex']:.2f} dex")
    xs2 = np.linspace(dd.x.min(), dd.x.max(), 20) + 44
    a1.plot(xs2, 10 ** (fs["a"][0] + fs["b"][0] * (xs2 - 44)), color=S.C_SDSS, lw=1.3, zorder=2,
            label=f"our fit, SDSS-RM: slope {fs['b'][0]:.2f} ± {0.5 * (fs['b'][1] + fs['b'][2]):.2f}, "
                  f"scatter {fs['sig'][0]:.2f} dex")
    n5548 = ds.object == "NGC5548"
    for m, col, lab in ((~n5548, S.C_LOCAL, f"local AGN ({loc['n_obj'] - 1} galaxies, HST host correction)"),
                        (n5548, S.C_N5548, "NGC 5548, 18 campaigns")):
        a1.errorbar(ds.logL5100[m], ds.tau_cent[m], xerr=ds.sx[m],
                    yerr=[ds.tau_cent_lo[m].clip(upper=ds.tau_cent[m] * 0.95), ds.tau_cent_hi[m]], fmt="o", color=col,
                    ms=4.2, mec=S.SURFACE, mew=0.5, elinewidth=0.6, capsize=0, zorder=4, label=lab)
    a1.errorbar(dd.logL, dd.tau, yerr=[dd.tau_lo.clip(upper=dd.tau * 0.95), dd.tau_hi], fmt="s", color=S.C_SDSS,
                ms=3.6, mec=S.SURFACE, mew=0.5, elinewidth=0.6, capsize=0, zorder=3, alpha=0.9,
                label=f"SDSS-RM quasars, z = 0.1–1 ({sd['n']}; Shen et al. 2024 lags)")
    a1.set_yscale("log")
    a1.set_ylim(0.8, 700)
    a1.set_yticks([1, 3, 10, 30, 100, 300])
    a1.set_yticklabels(["1", "3", "10", "30", "100", "300"])
    a1.set_ylabel("H-beta lag = BLR radius (light-days, rest frame)")
    a1.set_xlim(41.6, 45.7)
    a1.legend(loc="upper left", fontsize=7.5)
    plt.setp(a1.get_xticklabels(), visible=False)
    a1.set_title("Bigger AGN, bigger broad-line region: R ∝ L^0.5 holds over four decades of luminosity",
                 loc="left")
    # Residuals from the local fit
    res_l = ds.y - (fl["a"][0] + fl["b"][0] * ds.x)
    res_s = dd.y - (fl["a"][0] + fl["b"][0] * dd.x)
    a2.axhline(0, color=S.C_MODEL, lw=1.2)
    a2.plot(ds.logL5100[~n5548], res_l[~n5548], "o", color=S.C_LOCAL, ms=3.5, mec=S.SURFACE, mew=0.4)
    a2.plot(ds.logL5100[n5548], res_l[n5548], "o", color=S.C_N5548, ms=3.5, mec=S.SURFACE, mew=0.4)
    a2.plot(dd.logL, res_s, "s", color=S.C_SDSS, ms=3.2, mec=S.SURFACE, mew=0.4, alpha=0.9)
    a2.set_ylim(-1.3, 1.0)
    a2.set_ylabel("log (lag / fit)")
    a2.set_xlabel("log AGN luminosity λL(5100 Å) (erg/s), host galaxy removed")
    bins = np.arange(-1.3, 1.01, 0.1)
    a3.hist(res_l, bins=bins, orientation="horizontal", histtype="step", color=S.C_LOCAL, lw=1.6)
    a3.hist(res_s, bins=bins, orientation="horizontal", histtype="step", color=S.C_SDSS, lw=1.6)
    a3.axhline(np.median(res_s), color=S.C_SDSS, lw=1.0)
    a3.axhline(np.median(res_l), color=S.C_LOCAL, lw=1.0)
    _label(a3, a3.get_xlim()[1] * 0.97, -1.05, f"median offset\nSDSS-RM {np.median(res_s):+.2f} dex\n"
           f"local {np.median(res_l):+.2f} dex", S.INK_SECONDARY, dx=0, dy=0, ha="right", va="center")
    plt.setp(a3.get_yticklabels(), visible=False)
    a3.set_xlabel("number")
    fig.suptitle("The radius–luminosity relation: the key to weighing distant quasars from one spectrum",
                 x=0.06, ha="left", y=0.995)
    S.caption(fig, f"Local AGN: every H-beta campaign with an HST-based host-galaxy correction in the AGN Black Hole "
                   f"Mass Database ({loc['n']} lags; Bentz & Katz 2015). SDSS-RM: Shen et al. 2024 PyROA lags "
                   "detected at >2σ, luminosity = median over 2014–2020 minus host starlight. Fits: LINMIX_ERR "
                   "(Kelly 2007), which reproduces Bentz et al.'s own fit exactly on their table. Residuals are "
                   "relative to our local fit.")
    return S.save(fig, FIG / "05_radius_luminosity")


# ---------------------------------------------------------------------------
# 6. SDSS-RM: measuring lags in survey data
# ---------------------------------------------------------------------------


def fig_sdssrm(r: dict) -> list:
    sc = r["showcase"]
    ours = r["sdssrm_ours"]
    import rl as RL

    t = RL.sdssrm_table().set_index("RMID")
    fig = plt.figure(figsize=(11.5, 9.0))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.35], hspace=0.55, wspace=0.3)
    a1 = fig.add_subplot(gs[0, :])
    a2 = fig.add_subplot(gs[1, :2])
    a3 = fig.add_subplot(gs[1, 2])
    a4 = fig.add_subplot(gs[2, 0])
    a5 = fig.add_subplot(gs[2, 1])
    a6 = fig.add_subplot(gs[2, 2])
    c, h = sc["cont"], sc["line_lc"]
    lname = {"HB": "Hβ", "MG2": "Mg II", "C4": "C IV", "HA": "Hα"}[sc["line"]]
    w = sc["weighted"]
    lag = w["lag"]
    yr = lambda mjd: 1858.877 + np.asarray(mjd) / 365.25  # noqa: E731
    a1.plot(yr(c.t), _std(c.y), ".", color=S.C_CONT, ms=3, label="continuum, 11 years of imaging (SDSS, Pan-STARRS, CFHT, Bok)")
    a1.plot(yr(h.t), _std(h.y), "o", color=S.C_LINE, ms=3.2, mec=S.SURFACE, mew=0.4,
            label=f"broad {lname}, 7 years of SDSS spectra")
    a1.plot(yr(h.t - lag), _std(h.y), "o", mfc=S.SURFACE, mec=S.C_LINE, ms=3.2, mew=0.7,
            label=f"{lname} moved back by our lag ({lag:.0f} d observed)")
    a1.set_ylabel("flux (standardized)")
    a1.legend(loc="lower left", fontsize=7, ncol=3, bbox_to_anchor=(0, 1.0), borderaxespad=0.3)
    a1.set_title(f"A distant quasar: SDSS-RM {sc['rmid']:03d} at z = {sc['z']:.3f} (light from "
                 f"{_lookback(sc['z']):.1f} billion years ago)", loc="left", pad=24)
    # Lag distribution and alias weights
    res = sc["res"]
    lo_b, hi_b = w["bounds"]
    bins = np.arange(res.lags[0], res.lags[-1] + 20, 20.0)
    h_all, edges = np.histogram(res.cccd, bins=bins)
    inside = res.cccd[(res.cccd >= lo_b) & (res.cccd <= hi_b)]
    h_in, _ = np.histogram(inside, bins=bins)
    norm = h_all.max()
    centers = 0.5 * (edges[1:] + edges[:-1])
    a2.bar(centers, h_all / norm, width=np.diff(edges), color=S.C_REST, edgecolor=S.SURFACE, lw=0.4,
           label="lag of each resampling (all)")
    a2.bar(centers, h_in / norm, width=np.diff(edges), color=S.C_LINE, edgecolor=S.SURFACE, lw=0.4,
           label="primary peak: kept")
    a2.plot(res.lags, w["weights"], color=S.C_MODEL, lw=1.2, label="alias weight [N(τ)/N(0)]² ⊗ ACF")
    a2.set_ylim(0, 1.15)
    a2.set_xlabel("observed-frame lag (days)")
    a2.set_ylabel("relative (peak = 1)")
    pub = sc["pub_pyroa"]
    a2.set_title(f"Two candidate echoes: the main peak is kept ({lag:.0f} d; Shen et al. {pub[0]:.0f} d),\n"
                 "the bump near half a year, where shifted spectra land in seasonal gaps, is rejected",
                 loc="left", fontsize=9.5)
    a2.legend(loc="upper left", fontsize=7)
    # N(tau)/N(0) zoom
    a3.plot(res.lags, w["ntau"] / w["n0"], color=S.C_MODEL, lw=1.2)
    for k in (-2, -1, 1, 2):
        a3.axvline(k * 365.25, color=S.AXIS, lw=0.7, zorder=0)
    a3.set_xlim(res.lags[0], res.lags[-1])
    a3.set_ylim(0, 1.18)
    a3.set_xlabel("observed-frame lag (days)")
    a3.set_ylabel("N(τ)/N(0)")
    a3.set_title("Epochs that overlap", loc="left", fontsize=10)
    a3.text(0.5, 0.97, "dips at ±½ year: shifted\nspectra fall in the gaps", transform=a3.transAxes, ha="center",
            va="top", fontsize=7, color=S.INK_SECONDARY)
    # Ours vs published ICCF lags
    for ln, col, mk, lab in (("HB", S.C_LOCAL, "o", "Hβ"), ("MG2", S.C_LINE, "s", "Mg II"), ("C4", S.C_N5548, "^", "C IV")):
        o = ours[(ours.line == ln) & ours.ok].set_index("RMID")
        pub_det = t.loc[o.index, f"{ln}_lag"] > 0
        o = o[pub_det.values]
        p = t.loc[o.index, f"{ln}_iccf"]
        good = (p > 0) & (o.lag > 0)
        a4.plot(p[good], o.lag[good], mk, color=col, ms=3.4, mec=S.SURFACE, mew=0.4, label=lab, alpha=0.9)
    a4.plot([5, 3000], [5, 3000], color=S.AXIS, lw=0.8, zorder=0)
    a4.set_xscale("log")
    a4.set_yscale("log")
    a4.set_xlim(5, 3000)
    a4.set_ylim(5, 3000)
    a4.set_xlabel("published ICCF lag (days)")
    a4.set_ylabel("our ICCF lag (days)")
    a4.set_title("Same light curves, our code", loc="left", fontsize=10)
    a4.legend(loc="upper left", fontsize=7)
    comp0 = r["sdssrm_compare"]
    within = [100 * comp0[k]["frac_within_1sig_iccf"] for k in ("HB", "MG2", "C4")]
    a4.text(0.97, 0.05, f"{min(within):.0f}–{max(within):.0f}% agree\nwithin 1σ", transform=a4.transAxes, ha="right",
            fontsize=7.5, color=S.INK_SECONDARY)
    # Asymmetry test (all attempted H-beta)
    hb = ours[(ours.line == "HB") & ours.ok & np.isfinite(ours.lag)]
    a5.axvline(0, color=S.AXIS, lw=0.8)
    a5.axhline(0.4, color=S.AXIS, lw=0.8)
    hi = hb.rmax_near > 0.4
    a5.plot(hb.lag[hi], hb.rmax_near[hi], "o", color=S.C_LOCAL, ms=3.2, mec=S.SURFACE, mew=0.4)
    a5.plot(hb.lag[~hi], hb.rmax_near[~hi], "o", color=S.C_REST, ms=3.0, mec=S.SURFACE, mew=0.4)
    fpr = r["sdssrm_fpr"]
    a5.text(0.97, 0.97, f"r > 0.4: {fpr['pos_hi']} positive, {fpr['neg_hi']} negative", transform=a5.transAxes,
            ha="right", va="top", fontsize=7.5, color=S.INK)
    a5.text(0.97, 0.06, f"r < 0.4: {fpr['pos_lo']} positive, {fpr['neg_lo']} negative", transform=a5.transAxes,
            ha="right", fontsize=7.5, color=S.INK_SECONDARY)
    a5.set_xscale("symlog", linthresh=100)
    a5.set_xlim(-3000, 3000)
    a5.set_xticks([-1000, -100, 0, 100, 1000])
    a5.set_xticklabels(["−1000", "−100", "0", "100", "1000"])
    a5.set_ylim(0, 1)
    a5.set_xlabel("our H-beta lag (days, observed)")
    a5.set_ylabel("r_max near the lag")
    a5.set_title(f"All {fpr['n']} attempted Hβ lags", loc="left", fontsize=10)
    # Detected counts
    comp = r["sdssrm_compare"]
    names = [("HB", "Hβ"), ("MG2", "Mg II"), ("C4", "C IV"), ("HA", "Hα")]
    yy = np.arange(len(names))[::-1]
    a6.barh(yy + 0.18, [comp[k]["n_pub_det"] for k, _ in names], height=0.34, color=S.C_LIT, label="Shen et al. 2024")
    a6.barh(yy - 0.18, [comp[k]["n_both_det"] for k, _ in names], height=0.34, color=S.C_LINE,
            label="also detected by our ICCF")
    for y, (k, _) in zip(yy, names):
        _label(a6, comp[k]["n_pub_det"], y + 0.18, str(comp[k]["n_pub_det"]), S.INK_SECONDARY, dx=3)
        _label(a6, comp[k]["n_both_det"], y - 0.18, str(comp[k]["n_both_det"]), S.INK_SECONDARY, dx=3)
    a6.set_yticks(yy)
    a6.set_yticklabels([n for _, n in names])
    a6.set_xlabel("detected lags")
    a6.set_title("Detections", loc="left", fontsize=10)
    a6.legend(loc="lower right", fontsize=6.5)
    a6.grid(axis="y", visible=False)
    a6.set_xlim(0, 150)
    fig.suptitle("Echoes in a survey: 849 quasars, 11 years of photometry, 7 years of spectra (SDSS-RM)", x=0.06,
                 ha="left", y=1.0)
    S.caption(fig, "Data: SDSS Reverberation Mapping final release (Shen et al. 2024). Our lags: ICCF centroids, 500 "
                   "FR/RSS resamplings, alias weighting of Grier et al. (2017, 2019) as in pyPETAL. 'Detected': "
                   "positive at ≥2σ, r > 0.4 within ±1σ of the lag, less than half the resamplings rejected. "
                   "Published lags are Shen et al.'s fiducial PyROA lags; the comparison panel uses their ICCF lags.")
    return S.save(fig, FIG / "06_survey_echoes")


def _lookback(z: float) -> float:
    from astropy.cosmology import FlatLambdaCDM

    return float(FlatLambdaCDM(H0=70, Om0=0.3).lookback_time(z).value)


# ---------------------------------------------------------------------------
# 7. Black hole masses: calibrating f, and NGC 5548 on M-sigma
# ---------------------------------------------------------------------------


def fig_masses(r: dict) -> list:
    fd, fm, on = r["f_dyn"], r["f_msigma"], r["ngc5548_msigma"]
    fig = plt.figure(figsize=(11.5, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.42, wspace=0.28)
    a1 = fig.add_subplot(gs[0, 0])
    a2 = fig.add_subplot(gs[0, 1])
    a3 = fig.add_subplot(gs[1, :])
    # (a) individual f from dynamical models
    t = fd["table"].sort_values("logf").reset_index(drop=True)
    y = np.arange(len(t))
    mu, sc = fd["mean"][0], fd["scatter"][0]
    a1.axvspan(mu - sc, mu + sc, color=S.EXCLUDED, lw=0, zorder=0)
    a1.axvline(mu, color=S.C_MODEL, lw=1.3)
    n55 = t.object.str.contains("5548")
    a1.errorbar(t.logf[~n55], y[~n55], xerr=t.logf_err[~n55], fmt="o", color=S.C_LOCAL, ms=3.8, mec=S.SURFACE,
                mew=0.5, elinewidth=0.7, capsize=0)
    a1.errorbar(t.logf[n55], y[n55], xerr=t.logf_err[n55], fmt="o", color=S.C_N5548, ms=4.5, mec=S.SURFACE,
                mew=0.5, elinewidth=0.7, capsize=0)
    for yy in y[n55]:
        _label(a1, t.logf[yy] + t.logf_err[yy], yy, "NGC 5548", S.INK, dx=3)
    a1.set_yticks([])
    a1.set_xlabel("log f (= log of dynamical mass / virial product)")
    a1.set_title(f"30 AGN weighed by modelling their echoes:\n⟨log f⟩ = {mu:.2f} ± {0.5 * (fd['mean'][1] + fd['mean'][2]):.2f}, "
                 f"scatter {sc:.2f} dex (Shen et al.: 0.62, 0.31)", loc="left", fontsize=9.5)
    a1.grid(axis="y", visible=False)
    # (b) M-sigma
    q, ag = fm["quiescent"], fm["agn"]
    logf = fm["logf"][0]
    s_grid = np.linspace(1.3, 2.65, 50)
    a2.plot(10 ** s_grid, fm["alpha"][0] + fm["beta"][0] * (s_grid - np.log10(200)), color=S.C_MODEL, lw=1.3,
            zorder=1, label=f"joint fit: slope {fm['beta'][0]:.2f}, log f = {logf:.2f}")
    a2.plot(10 ** q.logsig, q.logM, "o", color=S.C_REST, ms=3.2, mec=S.SURFACE, mew=0.4,
            label=f"{len(q)} quiescent galaxies (van den Bosch 2016)", zorder=2)
    a2.errorbar(10 ** ag.logsig, ag.logVP + logf, yerr=ag.logVP_err, fmt="o", color=S.C_LOCAL, ms=4, mec=S.SURFACE,
                mew=0.5, elinewidth=0.6, capsize=0, zorder=3, label=f"{len(ag)} reverberation-mapped AGN (Woo et al. 2015)")
    m55 = on["ours_msigma"]["logM"]
    a2.errorbar([C.SIGMA_STAR_NGC5548], [m55], xerr=[C.SIGMA_STAR_NGC5548_ERR], yerr=[on["ours_msigma"]["logM_err"]],
                fmt="D", color=S.C_N5548, ms=6.5, mec=S.INK, mew=0.6, elinewidth=1.0, capsize=0, zorder=5)
    _label(a2, C.SIGMA_STAR_NGC5548, m55, "NGC 5548\n(our 2014 lag)", S.INK, dx=8, dy=-12, va="top")
    a2.set_xscale("log")
    a2.set_xlim(30, 420)
    a2.set_xticks([30, 50, 100, 200, 400])
    a2.set_xticklabels(["30", "50", "100", "200", "400"])
    a2.minorticks_off()
    a2.set_ylim(4.5, 10.8)
    a2.set_xlabel("stellar velocity dispersion σ* (km/s)")
    a2.set_ylabel("log black hole mass (M☉)")
    a2.set_title("…or by matching AGN to galaxies on M–σ\n(Woo et al. 2015: log f = 0.65 with Kormendy & Ho's 84)",
                 loc="left", fontsize=9.5)
    a2.legend(loc="upper left", fontsize=6.8)
    # (c) NGC 5548's mass, every route
    rows = [
        ("Reverberation, 2014 lag (ours), f from dynamical models (ours)", on["ours_dyn"]["logM"], on["ours_dyn"]["logM_err"], S.C_N5548),
        ("Reverberation, 2014 lag (ours), f from M–σ (ours)", on["ours_msigma"]["logM"], on["ours_msigma"]["logM_err"], S.C_N5548),
        ("Reverberation, 2014 lag (ours), f = 4.47 (Woo et al. 2015)", on["woo15"]["logM"], on["woo15"]["logM_err"], S.C_N5548),
        ("Reverberation, all 29 lines and seasons (virial test), f = 4.2", r["virial"]["logvp_weighted"] + fd["mean"][0],
         r["virial"]["logvp_scatter"] / np.sqrt(len(r["virial"]["sample"])), S.C_LOCAL),
        ("Pei et al. 2017 (2014 data, f = 4.47)", np.log10(LIT.PEI17_MBH[0]), 0.4343 * LIT.PEI17_MBH[1] / LIT.PEI17_MBH[0], S.C_LIT),
        ("Dynamical model of the echo, 2008 data (Pancoast et al. 2014)", np.log10(LIT.PANCOAST14_MBH[0]),
         0.5 * (np.log10(LIT.PANCOAST14_MBH[0] + LIT.PANCOAST14_MBH[2]) - np.log10(LIT.PANCOAST14_MBH[0] - LIT.PANCOAST14_MBH[1])), S.C_LIT),
        ("Dynamical model of the echo, 2014 data (Williams et al. 2020)", 7.54, 0.29, S.C_LIT),
        ("M–σ prediction for σ* = 195 km/s (our fit; its scatter)", on["msigma_prediction"], fm["eps"][0], S.C_REST),
    ]
    yy = np.arange(len(rows))[::-1]
    for y0, (lab, v, e, col) in zip(yy, rows):
        a3.errorbar(v, y0, xerr=e, fmt="o", color=col, ms=6, mec=S.SURFACE, mew=0.6, elinewidth=1.2, capsize=0)
        a3.text(1.01, y0, f"{10 ** v / 1e7:.1f} × 10⁷ M☉", transform=a3.get_yaxis_transform(), fontsize=8,
                color=S.INK_SECONDARY, va="center", ha="left")
    a3.set_yticks(yy)
    a3.set_yticklabels([r[0] for r in rows], fontsize=8, color=S.INK)
    a3.tick_params(axis="y", length=0)
    a3.set_xlim(6.9, 9.0)
    a3.set_ylim(-0.7, len(rows) - 0.3)
    a3.grid(axis="y", visible=False)
    a3.set_xlabel("log black hole mass of NGC 5548 (M☉)")
    a3.set_title("NGC 5548's black hole from its echoes: 3–9 × 10⁷ M☉, the spread mostly from the factor f "
                 "(M–σ alone would say 2 × 10⁸)", loc="left", fontsize=10)
    fig.suptitle("From light-days to solar masses: calibrating the BLR geometry factor f", x=0.06, ha="left", y=1.0)
    S.caption(fig, "Top left: Shen et al. 2024's compilation of BLR dynamical models (Pancoast et al. 2014; Grier et al. "
                   "2017; Williams et al. 2018, 2020; Villafana et al. 2022; ...), mean and intrinsic scatter refitted. "
                   "Top right: joint fit of quiescent galaxies and AGN virial products with one intrinsic scatter "
                   "(Woo et al. 2013, 2015 method); NGC 5548's σ* = 195 ± 13 km/s (Woo et al. 2010).")
    return S.save(fig, FIG / "07_black_hole_masses")


# ---------------------------------------------------------------------------
# 8. An AGN mapped in light-days: disk and broad-line region
# ---------------------------------------------------------------------------


def fig_structure(r: dict) -> list:
    import disk as DK

    dk = r["disk"]
    lags, fits = dk["lags"], dk["fits"]
    h = r["holiday"]
    st = r["storm"]["lags"]
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    lam = np.array([v["lam"] for v in lags.values()], float)
    tau = np.array([v["tau"] for v in lags.values()])
    pub = np.array([v["published"] for v in lags.values()])
    xs = np.linspace(1100, 9600, 200)
    f43 = fits["no_u_43"]["alpha"][0]
    ax.plot(xs, f43 * ((xs / DK.LAMBDA0) ** (4 / 3) - 1), color=S.C_CONT, lw=1.4, zorder=2)
    _label(ax, 6900, f43 * ((6900 / DK.LAMBDA0) ** (4 / 3) - 1), f"disk, fitted: τ ∝ λ^4/3, α = {f43:.2f} d",
           S.C_CONT, dx=-6, dy=10, ha="right", va="bottom")
    for md, lab, dy in ((0.1, "thin disk at 0.1 × Eddington", -8), (1.0, "thin disk at the Eddington limit", -8)):
        a = dk["alpha_thin"][md]
        ax.plot(xs, a * ((xs / DK.LAMBDA0) ** (4 / 3) - 1), color=S.C_LIT, lw=1.0, zorder=1)
        _label(ax, 7300, a * ((7300 / DK.LAMBDA0) ** (4 / 3) - 1), f"{lab} (α = {a:.2f} d)", S.C_LIT, dx=6,
               dy=dy, ha="left", va="top")
    ax.errorbar(lam, tau[:, 0], yerr=tau[:, 1:].T, fmt="o", color=S.C_CONT, ms=5, mec=S.SURFACE, mew=0.6,
                elinewidth=0.9, capsize=0, zorder=4, label="continuum lag vs 1367 Å (ours)")
    ax.plot(lam + 60, pub[:, 0], "_", color=S.INK, ms=9, mew=1.4, zorder=5, label="Fausnaugh et al. 2016")
    # Broad lines at their rest wavelengths: lags vs the 1158 A continuum
    lines = [("HeII", 1640, "He II", h["HeII"]["tau"]), ("CIV", 1549, "C IV", h["CIV"]["tau"]),
             ("SiIV", 1400, "Si IV", h["SiIV"]["tau"]), ("Lya", 1216, "Lyα", h["Lya"]["tau"])]
    z1 = 1 + C.Z_NGC5548
    for key, wl, name, t in lines:
        ax.errorbar(wl, t[0] / z1, yerr=[[t[1] / z1], [t[2] / z1]], fmt="D", color=S.C_LINE, ms=6, mec=S.SURFACE,
                    mew=0.6, elinewidth=0.9, capsize=0, zorder=4)
        _label(ax, wl, t[0] / z1, name, S.INK, dx=7, dy=0)
    for key, wl, name in (("Hb", 4861, "Hβ"), ("HeII", 4686, "He II 4686")):
        t = st[(key, "F1158", "full")]["cent"]
        ax.errorbar(wl, t[0], yerr=[[t[1]], [t[2]]], fmt="D", color=S.C_LINE, ms=6, mec=S.SURFACE, mew=0.6,
                    elinewidth=0.9, capsize=0, zorder=4)
        _label(ax, wl, t[0], name, S.INK, dx=7, dy=0)
    ax.plot([], [], "D", color=S.C_LINE, ms=6, label="broad-line lag vs 1158 Å (ours)")
    ax.axhline(0, color=S.AXIS, lw=0.8)
    hs, ls = ax.get_legend_handles_labels()
    order = [ls.index(k) for k in ("continuum lag vs 1367 Å (ours)", "Fausnaugh et al. 2016", "broad-line lag vs 1158 Å (ours)")]
    ax.set_xlim(1000, 9700)
    ax.set_ylim(-0.8, 8.0)
    ax.set_xlabel("rest wavelength (Å)")
    ax.set_ylabel("delay after the far-UV continuum (light-days)")
    ax.legend([hs[i] for i in order], [ls[i] for i in order], loc="upper left", fontsize=8)
    ratio = f43 / dk["alpha_thin"][0.1]
    ax.set_title(f"The disk's optical light comes from as far out as the inner broad-line gas: "
                 f"a disk {ratio:.1f}× larger than thin-disk theory", loc="left", fontsize=10.5)
    fig.suptitle("NGC 5548 measured in light-days (2014): accretion disk and broad-line region", x=0.06,
                 ha="left", y=1.0)
    S.caption(fig, "Continuum: HST/COS (1158-1746 Å) and ground-based u B g V r R i I z (Fausnaugh et al. 2016) "
                   "against HST 1367 Å, second-order detrending, 500 FR/RSS resamplings; fit τ = α[(λ/1367 Å)^4/3 − 1] "
                   "without u (Balmer-continuum contamination). Thin disk: Fausnaugh et al.'s eq. 8 with our mass "
                   f"({dk['mbh'] / 1e7:.1f} × 10⁷ M☉), radiative efficiency 0.1. Lines: UV lags over the first 75 days "
                   "of HST (holiday analysis), H-beta and He II 4686 over the full campaign.")
    return S.save(fig, FIG / "08_agn_in_light_days")


def make_all(r: dict) -> list:
    out = []
    for fn in (fig_echo, fig_holiday, fig_watch, fig_virial, fig_rl, fig_sdssrm, fig_masses, fig_structure):
        out += fn(r)
    return out
