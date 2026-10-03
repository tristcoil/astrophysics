"""Every figure, drawn from analysis results and artefacts (never recomputed).

Velocity figures show each point with its instrument zero point (and, for
ELODIE, its seasonal term) removed, so every spectrograph shares one axis.
Error bars are the quoted error and the fitted jitter in quadrature.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mpl_ticker
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

import plotstyle as ps
import rv_physics as phys
from config import FAP_THRESHOLD, FIGURE_DIR

ps.use_style()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _years(t):
    from astropy.time import Time
    return Time(np.asarray(t, float), format="jd").decimalyear


def corrected(fit):
    """Velocities with the fitted linear terms (zero points etc.) removed, and
    display errors (quoted error and white jitter in quadrature)."""
    m = fit.model
    y = m.data.y - m.offsets_and_trend(fit.theta)
    err = np.sqrt(m.data.err ** 2 + m.jitter(fit.theta)[0, m.data.inst] ** 2)
    return y, err


def scatter_by_family(ax, data, x, y, err=None, *, mask=None, ms=3.2, alpha=0.85,
                      elinewidth=0.6, zorder=3, label_counts=True):
    """Plot points coloured and shaped by spectrograph family. Returns legend handles."""
    fams = {}
    for k, name in enumerate(data.inst_names):
        fams.setdefault(ps.family(name), []).append(k)
    handles = []
    for fam, ks in fams.items():
        sel = np.isin(data.inst, ks)
        if mask is not None:
            sel &= mask
        if not sel.any():
            continue
        color, marker = ps.FAMILY_STYLE[fam]
        if err is not None:
            ax.errorbar(x[sel], y[sel], yerr=err[sel], fmt="none", ecolor=color,
                        elinewidth=elinewidth, alpha=0.35, zorder=zorder - 1)
        ax.plot(x[sel], y[sel], ls="none", marker=marker, ms=ms, mfc=color, mec=ps.SURFACE,
                mew=0.4, alpha=alpha, zorder=zorder)
        label = ps.FAMILY_LABEL[fam] + (f"  ({sel.sum()})" if label_counts else "")
        handles.append(Line2D([], [], ls="none", marker=marker, ms=5, mfc=color, mec=ps.SURFACE,
                              label=label))
    return handles


def densest_window(t, width):
    """Start of the ``width``-day window holding the most distinct nights."""
    nights = np.unique(np.floor(np.sort(t) - 0.3))
    best, start = -1, nights[0]
    for n0 in nights:
        k = np.sum((nights >= n0) & (nights < n0 + width))
        if k > best:
            best, start = k, n0
    return start


def fap_level(fap, t, err, fmax, n_offsets, n_eff):
    """Periodogram power at which the false-alarm probability equals ``fap``."""
    from scipy.optimize import brentq
    f = lambda z: phys.fap_baluev(z, t, err, fmax, n_offsets, n_eff) - fap
    return brentq(f, 1e-6, 0.999)


def period_axis(ax, periods=(20, 10, 4, 2, 1.3, 1, 0.8, 0.6, 0.5), label="period (days)"):
    """A second x axis on top giving the same frequency scale in days."""
    top = ax.secondary_xaxis("top", functions=(lambda f: f, lambda f: f))
    top.set_xticks([1 / p for p in periods])
    top.set_xticklabels([f"{p:g}" for p in periods])
    top.set_xlabel(label, color=ps.INK_SECONDARY)
    top.tick_params(colors=ps.MUTED)
    return top


# ---------------------------------------------------------------------------
# 51 Pegasi b
# ---------------------------------------------------------------------------

def fig_51peg_velocities(res, art, outdir=FIGURE_DIR):
    data = art["data"]
    fit = art["fit"]
    y, err = corrected(fit)
    yr = _years(data.t)

    fig = plt.figure(figsize=(10.5, 7.4))
    gs = GridSpec(2, 3, height_ratios=[1.25, 1.0], hspace=0.45, wspace=0.28, figure=fig)
    ax = fig.add_subplot(gs[0, :])
    handles = scatter_by_family(ax, data, yr, y, None, ms=3.0, alpha=0.7)
    ax.set_xlim(1994, 2019)
    ax.set_ylim(-95, 95)
    ax.set_ylabel("velocity of 51 Peg (m/s)")
    ax.set_title("51 Pegasi: 788 velocities from four spectrographs over 24 years", pad=26)
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4,
              handletextpad=0.3, columnspacing=1.4, borderaxespad=0.2)
    ax.text(0.995, 0.03, "at this scale the 4.2-day wobble is just a band ±56 m/s wide",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=ps.INK_SECONDARY)

    # Two zooms on the densest 12-day stretch of the early eras, and one on
    # HARPS's densest single night, where the orbit can be watched moving.
    from astropy.time import Time
    zooms = [("ELODIE", "The discovery instrument, 1995"),
             ("Lick", "Lick, the first autumn of confirmation")]
    for j, (fam, title) in enumerate(zooms):
        axz = fig.add_subplot(gs[1, j])
        ks = [k for k, n in enumerate(data.inst_names) if ps.family(n) == fam]
        sel_f = np.isin(data.inst, ks)
        tf = data.t[sel_f]
        if fam == "Lick":
            tf = tf[tf < 2450400]  # autumn 1995
        t0 = densest_window(tf, 12.0)
        m = sel_f & (data.t >= t0) & (data.t < t0 + 12.0)
        tg = np.linspace(t0 - 0.5, t0 + 12.5, 2000)
        axz.plot(tg - t0, fit.model.keplerian(fit.theta, tg)[0], color=ps.C_MODEL, lw=1.2, zorder=2)
        scatter_by_family(axz, data, data.t - t0, y, err, mask=m, ms=4.0, label_counts=False)
        axz.set_xlim(-0.5, 12.5)
        axz.set_ylim(-80, 80)
        axz.set_title(title, fontsize=9.5)
        axz.set_xlabel(f"days after {Time(t0, format='jd').iso[:10]}")
        if j == 0:
            axz.set_ylabel("velocity (m/s)")
        n_nights = np.unique(np.floor(data.t[m] - 0.3)).size
        axz.text(0.02, 0.03, f"{m.sum()} velocities, {n_nights} nights", transform=axz.transAxes,
                 fontsize=7.5, color=ps.INK_SECONDARY)
    axz = fig.add_subplot(gs[1, 2])
    ks = [k for k, n in enumerate(data.inst_names) if ps.family(n) == "HARPS"]
    th = data.t[np.isin(data.inst, ks)]
    nights = np.floor(th - 0.3)
    uniq, counts = np.unique(nights, return_counts=True)
    spans = np.array([np.ptp(th[nights == u]) for u in uniq])
    best = uniq[np.argmax(spans * (counts >= 10))]
    m = np.isin(data.inst, ks) & (np.floor(data.t - 0.3) == best)
    t0 = data.t[m].min()
    tg = np.linspace(t0 - 0.02, data.t[m].max() + 0.02, 400)
    axz.plot((tg - t0) * 24, fit.model.keplerian(fit.theta, tg)[0], color=ps.C_MODEL, lw=1.2, zorder=2)
    scatter_by_family(axz, data, (data.t - t0) * 24, y, err, mask=m, ms=4.0, label_counts=False)
    axz.set_title(f"HARPS, one night in {Time(t0, format='jd').iso[:4]}: watching it move", fontsize=9.5)
    axz.set_xlabel(f"hours after {Time(t0, format='jd').iso[:16]} UT")
    axz.text(0.03, 0.97, f"{m.sum()} velocities, ±{np.median(err[m]):.1f} m/s each;\n"
             "the night's shared offset is not removed", transform=axz.transAxes, fontsize=7.5,
             color=ps.INK_SECONDARY, va="top")
    ps.caption(fig, "Each instrument's zero point is fitted and removed (and ELODIE's seasonal "
               "term; see README). Black curve: our one-planet fit. Error bars: quoted error and "
               "fitted jitter in quadrature. Data: ELODIE archive (OHP); Fischer et al. 2014 (Lick); "
               "Butler et al. 2017 (Keck); HARPS RVBank (Trifonov et al. 2020).", y=0.0)
    return ps.save(fig, outdir / "01_51peg_velocities")


def fig_51peg_periodograms(res, art, outdir=FIGURE_DIR):
    disc = art["discovery"]
    steps = art["search"]
    freqs = art["freqs"]
    data = art["data"]
    fit = art["fit"]
    f_planet = 1 / res["search"][0]["peak_period"]
    f_sid = 1 / phys.SIDEREAL_DAY_D

    fig, axes = plt.subplots(4, 1, figsize=(10, 9.6), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1.25, 0.8, 0.9], "hspace": 0.32})
    xmax = 2.0

    # (a) the discovery data
    ax = axes[0]
    ax.plot(disc["freqs"], disc["power"], color=ps.C_SIGNAL, lw=0.8)
    ds = res["discovery_season"]
    ax.set_ylim(0, 1)
    ax.set_ylabel("power")
    ax.set_title(f"(a) Before the announcement: the {ds['n']} archived ELODIE spectra, "
                 f"{ds['first'][:7]} to {ds['last'][:7]} (the discovery paper used 35 velocities)", fontsize=10)
    ps.annotate(ax, f"{ds['peak']['period']:.3f} d  (FAP {ds['peak']['fap']:.1%})",
                (1 / ds["peak"]["period"], ds["peak"]["power"]), (0.36, 0.93), relpos=(0, 0.5))
    ax.text(0.99, 0.95, "21 points over one year: every peak has near-equal twins,\n"
            "one year and one day away", transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color=ps.INK_SECONDARY)

    # (b) everything
    ax = axes[1]
    s0 = steps[0]
    ax.plot(s0.freqs, s0.power, color=ps.C_SIGNAL, lw=0.7)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("power")
    ax.set_title(f"(b) All {data.n} velocities, 1994–2018: one planet, and its echoes", fontsize=10)
    pk = res["search"][0]
    ps.annotate(ax, f"51 Peg b: {pk['peak_period']:.4f} d", (f_planet, pk["peak_power"]),
                (0.33, 0.97), relpos=(0, 0.5))
    for fa, lab, xt in [(f_sid - f_planet, "alias 1.3048 d\n= 1 sidereal day − 51 Peg b", 0.52),
                        (f_sid + f_planet, "alias 0.8070 d\n= 1 sidereal day + 51 Peg b", 1.16)]:
        i = np.argmin(np.abs(s0.freqs - fa))
        j = i - 40 + int(np.argmax(s0.power[i - 40:i + 40]))
        ps.annotate(ax, lab, (s0.freqs[j], s0.power[j]), (xt, 0.62), relpos=(0.5, 0), ha="center")
    # inset: the yearly aliases around the planet
    ins = ax.inset_axes([0.745, 0.47, 0.235, 0.42])
    sel = np.abs(s0.freqs - f_planet) < 0.0045
    ins.plot(s0.freqs[sel], s0.power[sel], color=ps.C_SIGNAL, lw=0.8)
    for sgn in (-1, 1):
        ins.axvline(f_planet + sgn / phys.YEAR_D, color=ps.MUTED, lw=0.6)
    ins.set_xticks([f_planet - 1 / phys.YEAR_D, f_planet, f_planet + 1 / phys.YEAR_D])
    ins.set_xticklabels(["−1/yr", "4.2308 d", "+1/yr"], fontsize=6.5)
    ins.tick_params(axis="y", labelsize=6.5)
    ins.set_title("zoom: yearly aliases", fontsize=7, loc="left", pad=2, fontweight="normal")

    # (c) the window
    ax = axes[2]
    win = art["window"]
    ax.plot(freqs, win, color=ps.MUTED, lw=0.7)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("window")
    ax.set_title("(c) The sampling alone (spectral window): peaks at 0 and at 1 cycle per sidereal day",
                 fontsize=10)
    ax.text(0.99, 0.9, "every true signal is copied at f ± (each window peak)\n"
            "b's aliases sit at 1.0027 − 0.2364 and 1.0027 + 0.2364 cycles/day",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=ps.INK_SECONDARY)

    # (d) residuals
    ax = axes[3]
    s1 = steps[1]
    ax.plot(s1.freqs, s1.power, color=ps.C_SIGNAL, lw=0.7)
    err_eff = fit.model.effective_errors(fit.theta)
    lvl = fap_level(FAP_THRESHOLD, data.t, err_eff, freqs.max(), data.n_inst, fit.model.n_effective())
    ax.axhline(lvl, color=ps.INK_SECONDARY, lw=0.8)
    ax.text(1.99, lvl + 0.004, f"false-alarm probability {FAP_THRESHOLD:.1%}", ha="right", va="bottom",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_ylim(0, max(lvl * 1.35, s1.peak_power * 1.3))
    ax.set_ylabel("power")
    ax.set_title(f"(d) After removing 51 Peg b: nothing left that clears the bar", fontsize=10)
    nx = res["next_signal"]
    ps.annotate(ax, f"highest: {nx['period']:.2f} d (FAP {nx['fap']:.0%})",
                (1 / nx["period"], nx["power"]), (0.25, lvl * 1.2), relpos=(0, 0.5))
    ax.set_xlim(0, xmax)
    ax.set_xlabel("frequency (cycles per day)")
    period_axis(axes[0])
    ps.caption(fig, "Generalised Lomb–Scargle periodograms (Zechmeister & Kürster 2009), fitting "
               "each instrument's zero point at every frequency. Power = fraction of the variance a "
               "sinusoid removes. False-alarm probabilities: Baluev (2008), counting each night of "
               "multi-exposure instruments once. Periods of 0.5 d to twice the baseline searched.", y=0.02)
    return ps.save(fig, outdir / "02_51peg_periodogram")


def _summary_panel(axt, rows, *, y0=0.98, dy=0.05, fs=8):
    """Two-column text block: bold headers (value None) and key/value rows."""
    yy = y0
    for key, val in rows:
        if val is None:
            yy -= dy * 0.35 if key else 0
            axt.text(0.0, yy, key, fontsize=fs + 1, fontweight="bold", color=ps.INK,
                     transform=axt.transAxes, va="top")
        else:
            axt.text(0.0, yy, key, fontsize=fs, color=ps.INK_SECONDARY, transform=axt.transAxes, va="top")
            axt.text(1.0, yy, val, fontsize=fs, color=ps.INK, ha="right", transform=axt.transAxes, va="top")
        yy -= dy
    return yy


def fig_51peg_orbit(res, art, outdir=FIGURE_DIR):
    data = art["data"]
    fit = art["fit"]
    y, err = corrected(fit)
    pl = fit.planets()[0]
    P, tc = pl["P"], pl["tc"]
    phase = np.mod((data.t - tc) / P + 0.25, 1.0) - 0.25  # conjunction at 0
    model_at = fit.model.keplerian(fit.theta)[0]
    resid = y - model_at

    fig = plt.figure(figsize=(10.8, 6.4))
    gs = GridSpec(2, 2, width_ratios=[2.25, 1], height_ratios=[3, 1.1], hspace=0.08, wspace=0.10,
                  figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)
    ph_grid = np.linspace(-0.25, 0.75, 600)
    curve = fit.model.keplerian(fit.theta, tc + ph_grid * P)[0]
    ax.axvline(0, color=ps.AXIS, lw=0.8, zorder=1)
    ax.plot(ph_grid, curve, color=ps.C_MODEL, lw=1.6, zorder=4)
    handles = scatter_by_family(ax, data, phase, y, err, ms=3.0, alpha=0.75)
    ax.set_ylabel("velocity of 51 Peg (m/s)")
    ax.set_title("51 Pegasi b: the star's 4.23-day orbit, folded", pad=8)
    ax.legend(handles=handles, loc="lower right", fontsize=7.5, handletextpad=0.3)
    ax.set_ylim(-100, 100)
    ax.tick_params(labelbottom=False)
    ax.axhline(0, color=ps.AXIS, lw=0.6, zorder=1)
    # Positive velocity = moving away from us (redshift).
    ax.text(0.25, 95, "↑ above zero: star moving away from us (redshift)", ha="center", va="top",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.text(0.25, -72, "↓ below zero: star moving towards us (blueshift)", ha="center", va="top",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.text(0.012, -92, "planet passes in front of the star", ha="left", va="bottom", fontsize=7.5,
            color=ps.INK_SECONDARY)

    scatter_by_family(axr, data, phase, resid, err, ms=2.6, alpha=0.7)
    axr.axhline(0, color=ps.C_MODEL, lw=1.0)
    axr.set_ylim(-45, 45)
    axr.set_xlim(-0.25, 0.75)
    axr.set_ylabel("data − fit (m/s)")
    axr.set_xlabel("orbital phase (0 = conjunction: planet between us and the star)")

    axt = fig.add_subplot(gs[:, 1])
    axt.axis("off")
    fc = res["fits"]["circular"]["planet"]
    fe = res["fits"]["eccentric"]["planet"]
    ls = res["fits"]["lucy_sweeney"]
    mq = res["discovery_season"]["mayor_queloz_1995"]
    brogi = next(iter(res["true_mass"].values()))
    harps_rms = [r["rms"] for r in res["instruments"] if r["instrument"].startswith("HARPS")]
    rows = [
        ("Our fit, circular orbit", None),
        ("period", f"{fc['P']['median']:.7f} ± {fc['P']['plus']:.7f} d"),
        ("semi-amplitude K", f"{fc['K']['median']:.1f} ± {fc['K']['plus']:.1f} m/s"),
        ("eccentricity", f"< {fe['e']['upper_95']:.3f}  (95%)"),
        ("orbit radius", f"{fc['a_au']['median']:.4f} AU"),
        ("m sin i", f"{fc['msini_mjup']['median']:.3f} ± {fc['msini_mjup']['plus']:.3f} M_Jup"),
        ("", f"= {fc['msini_mearth']['median']:.0f} ± {fc['msini_mearth']['plus']:.0f} Earth masses"),
        ("scatter about fit", f"{res['fits']['circular']['rms']:.1f} m/s (HARPS {min(harps_rms):.1f})"),
        ("Mayor & Queloz 1995", None),
        ("period", f"{mq['P']:.4f} ± {mq['P_err']:.4f} d"),
        ("semi-amplitude K", f"{mq['K']:.0f} ± {mq['K_err']:.0f} m/s"),
        ("m sin i", f"{mq['msini_mjup']:.2f} M_Jup"),
        ("velocities", f"{mq['N']}, scatter 13 m/s"),
        ("True mass, no sin i", None),
        ("K★ / K_planet × M★", f"{brogi['mass_mjup']['median']:.3f} ± {brogi['mass_mjup']['plus']:.3f} M_Jup"),
        ("inclination", f"> {brogi['incl_deg_95_lower']:.0f}°  (95%)"),
    ]
    yy = _summary_panel(axt, rows, dy=0.049)
    axt.text(0.0, yy - 0.02, f"The fitted eccentricity, {ls['e_hat']:.3f} ± {ls['sigma_e']:.3f}, is what a\n"
             f"circular orbit shows by chance {ls['p_circular']:.0%} of the time\n(Lucy & Sweeney 1971): "
             "the orbit is consistent\nwith circular, as tides predict for a 4-day orbit.",
             fontsize=7.5, color=ps.INK_SECONDARY, transform=axt.transAxes, va="top")
    ps.caption(fig, "m sin i uses M★ = 1.07 ± 0.05 M☉ (Cont et al. 2026). True mass: our K★ over the "
               "planet's own orbital velocity K_p = 134.1 ± 1.8 km/s, measured from its CO and water "
               "lines (Brogi et al. 2013); K★/K_p = m/M★ needs no inclination.", y=0.02)
    return ps.save(fig, outdir / "03_51peg_orbit")


def make_51peg(res, art, outdir=FIGURE_DIR):
    out = []
    out += fig_51peg_velocities(res, art, outdir)
    out += fig_51peg_periodograms(res, art, outdir)
    out += fig_51peg_orbit(res, art, outdir)
    return out


# ---------------------------------------------------------------------------
# 55 Cancri
# ---------------------------------------------------------------------------

PLANET_ORDER_55CNC = ["b", "c", "d", "e", "f"]


def _phase_fold_panel(ax, fit, data, j, *, y, err, n_bins=12, label=""):
    """Velocities with every other signal removed, folded on signal j."""
    P, K, e, w, tc = (a[0] for a in fit.model.planet_arrays(fit.theta))
    others = [k for k in range(P.size) if k != j]
    y_j = y - fit.model.keplerian(fit.theta, planets=others)[0]
    phase = np.mod((data.t - tc[j]) / P[j] + 0.25, 1.0) - 0.25
    ax.plot(phase, y_j, ls="none", marker="o", ms=1.8, color=ps.CONTEXT, zorder=1)
    edges = np.linspace(-0.25, 0.75, n_bins + 1)
    idx = np.digitize(phase, edges) - 1
    w_ = 1 / err ** 2
    mean = np.array([np.average(y_j[idx == b], weights=w_[idx == b]) if np.any(idx == b) else np.nan
                     for b in range(n_bins)])
    sem = np.array([np.std(y_j[idx == b]) / np.sqrt(max((idx == b).sum() - 1, 1)) if np.any(idx == b) else np.nan
                    for b in range(n_bins)])
    ctr = 0.5 * (edges[1:] + edges[:-1])
    grid = np.linspace(-0.25, 0.75, 400)
    curve = fit.model.keplerian(fit.theta, tc[j] + grid * P[j], planets=[j])[0]
    ax.plot(grid, curve, color=ps.C_MODEL, lw=1.4, zorder=3)
    ax.errorbar(ctr, mean, yerr=sem, fmt="o", ms=4, color=ps.BLUE, mec=ps.SURFACE, mew=0.6,
                elinewidth=0.9, zorder=4)
    lim = max(3.2 * K[j], 6.0)
    ax.set_ylim(-lim, lim)
    ax.set_xlim(-0.25, 0.75)
    ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="top")


def fig_55cnc_planets(res, art, outdir=FIGURE_DIR):
    data = art["data"]
    steps = art["search"]
    fit = art["fit"]
    y, err = corrected(fit)
    fig = plt.figure(figsize=(11.5, 12.8))
    gs = GridSpec(6, 2, width_ratios=[1.55, 1], hspace=0.62, wspace=0.22, figure=fig)
    names = [s["identified_as"] for s in res["search"]]
    blurbs = {"b": "b: 14.65 d, a giant planet", "d": "d: ~15 years, a Jupiter-like orbit",
              "c": "c: 44.4 d", "e": "e: 0.74 d, a hot super-Earth", "f": "f: 262 d"}
    for k, s in enumerate(steps):
        ax = fig.add_subplot(gs[k, 0])
        per = 1 / s.freqs
        ax.plot(per, s.power, color=ps.C_SIGNAL, lw=0.6)
        ax.set_xscale("log")
        ax.set_xlim(0.5, 2e4)
        ax.set_ylim(0, max(0.35, s.peak_power * 1.35))
        ax.set_ylabel("power")
        nm = names[k]
        if k < 5:
            title = f"step {k + 1}: {k} signal{'s' if k != 1 else ''} removed → strongest peak {s.peak_period:.4g} d"
            lab = blurbs.get(nm, nm)
        else:
            title = "step 6: after five planets"
            lab = f"{s.peak_period:.0f} d left over (FAP {s.fap:.0e}): see bottom right"
        ax.set_title(title, fontsize=9.5)
        if k < 5:
            ps.annotate(ax, lab, (s.peak_period, s.peak_power),
                        (s.peak_period * (0.05 if s.peak_period > 50 else 3.0), s.peak_power * 1.08),
                        relpos=(1, 0.5) if s.peak_period > 50 else (0, 0.5),
                        ha="right" if s.peak_period > 50 else "left")
        else:
            ps.annotate(ax, lab, (s.peak_period, s.peak_power), (2.0, s.peak_power * 1.12),
                        relpos=(1, 0.5), ha="left")
        if nm == "e":
            pa = res["e_step"]["best_period_if_search_stops_at_1d"]
            j = int(np.argmin(np.abs(per - pa)))
            j = j - 20 + int(np.argmax(s.power[j - 20:j + 21]))
            ps.annotate(ax, f"its 1-day alias, {pa:.3f} d:\nthe answer if you only search P > 1 d",
                        (per[j], s.power[j]), (25, s.peak_power * 0.62), relpos=(0, 0.5))
        if k == 5:
            ax.set_xlabel("period (days)")
    # right column: folded planets and the long-period signal
    labels = {0: "b", 1: "c", 3: "e", 4: "f"}
    pl = {p["label"].split()[-1]: p for p in res["fiducial"]["planets"][:5]}
    rows = [(0, "b"), (1, "c"), (3, "e"), (4, "f")]
    for r_i, (j, nm) in enumerate(rows):
        ax = fig.add_subplot(gs[r_i, 1])
        p = pl[nm]
        _phase_fold_panel(ax, fit, data, j, y=y, err=err,
                          label=f"55 Cnc {nm}: K = {p['K']['median']:.1f} ± {p['K']['plus']:.1f} m/s,"
                                f" m sin i = {_mass_str(p)}")
        ax.set_ylabel("m/s")
        if r_i == 3:
            ax.set_xlabel("orbital phase")
    # d and the magnetic cycle, against time
    ax = fig.add_subplot(gs[4:, 1])
    P, K, e, w, tc = (a[0] for a in fit.model.planet_arrays(fit.theta))
    y_long = y - fit.model.keplerian(fit.theta, planets=[0, 1, 3, 4])[0]
    yr = _years(data.t)
    handles = scatter_by_family(ax, data, yr, y_long, None, ms=2.2, alpha=0.5)
    tg = np.linspace(data.t.min() - 200, data.t.max() + 200, 1500)
    ax.plot(_years(tg), fit.model.keplerian(fit.theta, tg, planets=[2, 5])[0], color=ps.C_MODEL, lw=1.6,
            label="planet d + magnetic cycle")
    ax.plot(_years(tg), fit.model.keplerian(fit.theta, tg, planets=[5])[0], color=ps.INK_SECONDARY, lw=1.0,
            label="cycle term alone")
    dd = pl["d"]
    ax.set_title(f"55 Cnc d ({dd['P']['median'] / 365.25:.0f}-yr orbit) and the star's magnetic cycle",
                 fontsize=9.5, loc="left", pad=14)
    ax.set_xlabel("year")
    ax.set_ylabel("m/s (b, c, e, f removed)")
    ax.set_ylim(-110, 110)
    ax.legend(handles=handles + [Line2D([], [], color=ps.C_MODEL, lw=1.6, label="planet d + cycle"),
                                 Line2D([], [], color=ps.INK_SECONDARY, lw=1.0, label="cycle term alone")],
              fontsize=7, loc="lower left", ncol=2)
    fig.suptitle("55 Cancri: five planets, found one at a time from 25 years of Lick and Keck velocities",
                 x=0.06, ha="left", y=0.935)
    ps.caption(fig, "Left: each periodogram is of the residuals after fitting (jointly) every signal found so far; "
               "the strongest peak is added and the loop repeats. Right: velocities with all other signals removed, "
               "grey = individual, blue = phase-binned means. m sin i with M★ = 0.905 ± 0.015 M☉ "
               "(von Braun et al. 2011); with Moutou et al. 2026's 1.015 M☉ every mass is 7.9% higher.", y=0.08)
    return ps.save(fig, outdir / "04_55cnc_planets")


def _mass_str(p) -> str:
    me = p["msini_mearth"]["median"]
    if me > 100:
        return f"{p['msini_mjup']['median']:.2f} M_Jup"
    return f"{me:.1f} ± {p['msini_mearth']['plus']:.1f} M⊕"


def fig_55cnc_e(res, art, outdir=FIGURE_DIR):
    from fetch_data import load_zeng
    tess = art["tess"]
    at = art["alias_test"]
    fig = plt.figure(figsize=(11.5, 11.0))
    gs = GridSpec(3, 2, hspace=0.5, wspace=0.25, figure=fig)

    # (a) the alias pair in the periodogram of the data without planet e
    ax = fig.add_subplot(gs[0, 0])
    freqs = art["freqs"]
    pw = art["no_e_power"]
    sel = (freqs > 0.2) & (freqs < 1.6)
    ax.plot(freqs[sel], pw[sel], color=ps.C_SIGNAL, lw=0.6)
    for P, lab, xt, yt in [(0.736547, "0.7365 d:\nthe real planet", 1.13, 0.75),
                           (res["e_step"]["best_period_if_search_stops_at_1d"], "2.817 d:\nits alias", 0.52, 0.9)]:
        j = np.argmin(np.abs(freqs - 1 / P))
        j = j - 20 + int(np.argmax(pw[j - 20:j + 20]))
        ps.annotate(ax, lab, (freqs[j], pw[j]), (xt, pw[j] * yt), relpos=(0.5, 0), ha="center")
    ax.set_xlabel("frequency (cycles per day)")
    ax.set_ylabel("power")
    ax.set_ylim(0, pw[sel].max() * 1.35)
    ax.set_title("(a) Two peaks one sidereal day apart: which is 55 Cnc e?", fontsize=9.5)
    hist = [r for r in res["alias_history"]]
    txt = "published as the 2.8-day planet:\n" + "\n".join(
        f"  {r['ref']}: {r['pl_orbper']:.3f} d" for r in sorted(hist, key=lambda r: r["year"] or 0))
    txt += "\nflagged as its alias:\n  Butler et al. 2017: 2.817 d"
    ax.text(0.45, 0.97, txt, transform=ax.transAxes, ha="center", va="top", fontsize=7, color=ps.INK_SECONDARY)

    # (b) the alias test
    ax = fig.add_subplot(gs[0, 1])
    keys = list(at["hypotheses"].keys())
    cols = [ps.BLUE, ps.ORANGE]
    bins = np.linspace(0.3, 3.5, 60)
    for key, c in zip(keys, cols):
        r = np.clip(at["hypotheses"][key]["sim_ratios"], bins[0], bins[-1])
        ax.hist(r, bins=bins, color=c, alpha=0.55, label=f"simulated: true period {float(key):.4f} d")
    ax.axvline(at["observed_ratio"], color=ps.INK, lw=1.4)
    ax.text(at["observed_ratio"] - 0.04, ax.get_ylim()[1] * 0.92, f"observed\n{at['observed_ratio']:.2f}",
            fontsize=7.5, va="top", ha="right")
    ax.set_xlabel("peak height at 0.7365 d ÷ peak height at 2.817 d")
    ax.set_ylabel("simulations")
    ax.legend(fontsize=7, loc="center right")
    h1 = res["alias_test"]["hypotheses"][keys[1]]
    ax.set_title("(b) Simulated data with each period as the truth (Dawson & Fabrycky 2010)", fontsize=9.5)
    h0 = res["alias_test"]["hypotheses"][keys[0]]
    n1 = len(at["hypotheses"][keys[1]]["sim_ratios"])
    far = (f"{h1['fraction_as_extreme']:.1%}" if h1["fraction_as_extreme"] > 0 else f"none of {n1}")
    ax.text(0.99, 0.97, f"a ratio this large (or larger):\n{h0['fraction_as_extreme']:.1%} of simulations if 0.7365 d is real\n"
            f"{far} if 2.817 d is real", transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, color=ps.INK_SECONDARY)

    # (c), (d) TESS folds
    for k, (key, title) in enumerate([("true", "(c) TESS, folded on 0.7365 d: a transit"),
                                      ("alias", "(d) TESS, folded on 2.817 d: nothing")]):
        ax = fig.add_subplot(gs[1, k])
        fo = tess["folds"][key]
        hours = fo["phase"] * 24
        m = np.abs(hours) < 3.5 if key == "true" else np.abs(hours) < 12
        ax.errorbar(hours[m], (fo["flux"][m] - 1) * 1e6, yerr=fo["err"][m] * 1e6, fmt="o", ms=3,
                    color=ps.BLUE, mec=ps.SURFACE, mew=0.5, elinewidth=0.7)
        if key == "true":
            ax.plot(art["tess"]["model_phase"] * 24, (art["tess"]["model_flux"] - 1) * 1e6, color=ps.C_MODEL, lw=1.4)
        ax.set_ylim(-520, 220)
        ax.axhline(0, color=ps.AXIS, lw=0.6)
        ax.set_xlabel("hours from mid-transit" if key == "true" else "hours from predicted transit")
        ax.set_ylabel("brightness change (ppm)")
        ax.set_title(title, fontsize=9.5)
    tr_res = res["tess"]

    # (e) HARPS: two nights the orbit fit never saw
    ax = fig.add_subplot(gs[2, 0])
    h = art["harps"]
    hd = h["data"]
    nights = np.floor(hd.t - 0.3)
    others = h["all"] - h["e"]
    for n_i, u in enumerate(np.unique(nights)):
        m = nights == u
        v = hd.y[m] - others[m]
        v = v - np.average(v - h["e"][m], weights=1 / hd.err[m] ** 2)
        x = (hd.t[m] - hd.t[m].min()) * 24 + n_i * 11
        ax.errorbar(x, v, yerr=hd.err[m], fmt="o", ms=3, color=ps.BLUE, mec=ps.SURFACE, mew=0.4, elinewidth=0.6)
        tg = np.linspace(hd.t[m].min() - 0.02, hd.t[m].max() + 0.02, 200)
        from astropy.time import Time
        ax.text(n_i * 11, 9.5, Time(u + 0.5, format="jd").iso[:10], fontsize=7.5, color=ps.INK_SECONDARY)
    fitm = art["fit"]
    for n_i, u in enumerate(np.unique(nights)):
        m = nights == u
        tg = np.linspace(hd.t[m].min() - 0.02, hd.t[m].max() + 0.02, 200)
        ax.plot((tg - hd.t[m].min()) * 24 + n_i * 11, fitm.model.keplerian(fitm.theta, tg, planets=[3])[0],
                color=ps.C_MODEL, lw=1.3)
    hc = res["harps_check"]
    ax.set_ylim(-11, 11)
    ax.set_xlabel("hours (two nights, placed side by side)")
    ax.set_ylabel("velocity from e alone (m/s)")
    ax.set_title("(e) HARPS, 2012: planet e predicted from Keck and Lick alone", fontsize=9.5)
    ax.text(0.99, 0.04, f"scatter {hc['rms_with_e']:.2f} m/s with e, {hc['rms_without_e']:.2f} without",
            transform=ax.transAxes, ha="right", fontsize=7.5, color=ps.INK_SECONDARY)

    # (f) mass-radius
    ax = fig.add_subplot(gs[2, 1])
    styles = {"iron": ps.INK_SECONDARY, "earthlike": ps.GREEN, "rock": ps.ORANGE}
    from config import ZENG_MODELS
    for key, (fname, lab) in ZENG_MODELS.items():
        mm, rr = load_zeng(key)
        s = (mm > 0.5) & (mm < 30)
        ax.plot(mm[s], rr[s], color=styles[key], lw=1.2, label=lab)
    mt = tr_res["true_mass_earth"]; rp = tr_res["radius_earth"]
    ax.errorbar([mt["median"]], [rp["median"]], xerr=[[mt["minus"]], [mt["plus"]]],
                yerr=[[rp["minus"]], [rp["plus"]]], fmt="o", ms=7, color=ps.BLUE, mec=ps.SURFACE, mew=1.2,
                elinewidth=1.2, zorder=5)
    ax.plot([1], [1], marker="o", ms=6, color=ps.INK, ls="none")
    ax.text(1.08, 0.97, "Earth", fontsize=7.5, va="top")
    ps.annotate(ax, f"55 Cnc e\n{mt['median']:.1f} ± {mt['plus']:.1f} M⊕, {rp['median']:.2f} ± {rp['plus']:.2f} R⊕\n"
                f"density {tr_res['density_gcc']['median']:.1f} ± {tr_res['density_gcc']['plus']:.1f} g/cm³",
                (mt["median"], rp["median"]), (0.8, 2.3), relpos=(1, 0.5), ha="left", va="center")
    ax.legend(fontsize=7, loc="lower right", title="Zeng et al. (2019)", title_fontsize=7)
    ax.set_xscale("log")
    ax.set_xlim(0.7, 30)
    ax.set_ylim(0.8, 2.6)
    ax.set_xlabel("mass (Earth masses)")
    ax.set_ylabel("radius (Earth radii)")
    ax.set_title("(f) Mass from our velocities, radius from TESS", fontsize=9.5)
    fig.suptitle("55 Cancri e: the alias, the transit, and a lava world's density", x=0.06, ha="left", y=0.95)
    ps.caption(fig, f"TESS sectors {', '.join(map(str, tr_res['sectors']))} (SPOC 2-min), flattened with a "
               f"12-h running median, binned in 10 min. Transit fit: quadratic limb darkening from Claret (2017), "
               f"R★ = 0.943 ± 0.010 R☉ (von Braun et al. 2011). Our RV orbit (data 1989–2014) predicts the TESS "
               f"transits to {tr_res['rv_predicted_transit_offset_minutes']:+.0f} ± "
               f"{tr_res['rv_predicted_transit_sigma_minutes']:.0f} min. Composition curves: Zeng et al. (2019).",
               y=0.03)
    return ps.save(fig, outdir / "05_55cnc_e")


def make_55cnc(res, art, outdir=FIGURE_DIR):
    return fig_55cnc_planets(res, art, outdir) + fig_55cnc_e(res, art, outdir)


# ---------------------------------------------------------------------------
# Activity: Proxima Centauri and Gliese 581
# ---------------------------------------------------------------------------

def _binned_fold(ax, t, y, err, P, t0, *, n_bins=10, curve=None, color=ps.BLUE, label=None, shift=0.0):
    """Phase-fold, draw single points in grey and weighted bin means with
    their errors (sqrt(1/sum w), inflated by sqrt(reduced chi^2) when the bin
    scatters more than its error bars say)."""
    phase = np.mod((t - t0) / P + 0.25, 1.0) - 0.25
    ax.plot(phase, y, ls="none", marker="o", ms=1.8, color=ps.CONTEXT, zorder=1)
    edges = np.linspace(-0.25, 0.75, n_bins + 1)
    idx = np.digitize(phase, edges) - 1
    w = 1 / err ** 2
    ctr = 0.5 * (edges[1:] + edges[:-1])
    mean = np.full(n_bins, np.nan)
    se = np.full(n_bins, np.nan)
    for b in range(n_bins):
        m = idx == b
        if m.sum() < 2:
            continue
        mean[b] = np.average(y[m], weights=w[m])
        chi2r = np.sum(w[m] * (y[m] - mean[b]) ** 2) / (m.sum() - 1)
        se[b] = np.sqrt(1 / np.sum(w[m])) * np.sqrt(max(chi2r, 1.0))
    ax.errorbar(ctr + shift, mean, yerr=se, fmt="o", ms=4.5, color=color, mec=ps.SURFACE, mew=0.6, elinewidth=1.0,
                zorder=4, label=label)
    if curve is not None:
        ax.plot(curve[0], curve[1], color=ps.C_MODEL, lw=1.4, zorder=3)
    ax.set_xlim(-0.25, 0.75)


def fig_proxima(res, art, outdir=FIGURE_DIR):
    data = art["data"]
    fit = art["fit"]
    steps = art["search"]
    fig = plt.figure(figsize=(11.5, 10.8))
    gs = GridSpec(4, 4, height_ratios=[1.25, 1, 1, 1.35], hspace=0.75, wspace=0.32, figure=fig)
    P_b = res["fiducial"]["planet"]["P"]["median"]
    rot = (83.0, 92.0)

    ax = fig.add_subplot(gs[0, :])
    s0 = steps[0]
    per = 1 / s0.freqs
    ax.axvspan(*rot, color=ps.CONTEXT, lw=0)
    ax.plot(per, s0.power, color=ps.C_SIGNAL, lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.9, 6000)
    ax.set_ylim(0, 0.26)
    ax.set_ylabel("power")
    ax.set_xlabel("period (days)")
    ax.set_title(f"(a) Proxima Centauri's velocities ({data.n} HARPS spectra, 2004–2019): the strongest peak",
                 fontsize=10)
    ps.annotate(ax, f"Proxima b: {s0.peak_period:.3f} d  (false-alarm probability {s0.fap:.0e})",
                (s0.peak_period, s0.peak_power), (1.3, 0.235), relpos=(0, 0.5))
    ax.text(np.sqrt(rot[0] * rot[1]), 0.245, "star's rotation\n(83–92 d)", ha="center", va="top", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ax.text(0.99, 0.95, "near 1 day: daily aliases of long-period signals", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5, color=ps.INK_SECONDARY)

    from analysis import ACTIVITY_LABELS
    inds = [k for k in ("Halpha", "NaD1", "RHK", "CRX", "dLW", "FWHM", "BIS") if k in art["indicators"]]
    for i, key in enumerate(inds):
        ax = fig.add_subplot(gs[1 + i // 4, i % 4])
        ia = art["indicators"][key]
        per = 1 / ia["freqs"]
        ax.axvspan(*rot, color=ps.CONTEXT, lw=0)
        ax.axvline(P_b, color=ps.INK_SECONDARY, lw=0.8)
        ax.plot(per, ia["power"], color=ps.C_SIGNAL, lw=0.5)
        ax.set_xscale("log")
        ax.set_xlim(2, 6000)
        top = max(0.12, ia["power"].max() * 1.15)
        ax.set_ylim(0, top)
        at = list(res["indicators"][key]["power_at_planets"].values())[0]
        ax.set_title(ACTIVITY_LABELS[key].split(" (")[0], fontsize=8.5)
        ax.text(0.03, 0.93, f"at 11.19 d: {at:.3f}", transform=ax.transAxes, fontsize=7, va="top",
                color=ps.INK_SECONDARY)
        if i % 4 == 0:
            ax.set_ylabel("power")
        if i // 4 == 1:
            ax.set_xlabel("period (days)")
    axl = fig.add_subplot(gs[2, 3])
    axl.axis("off")
    axl.text(0.0, 0.95, "(b) The same periodogram for seven\nactivity indicators, from the same\n"
             "spectra. Grey band: rotation.\nVertical line: the planet's period.\n\n"
             "None shows the 11.19-day signal, so\nit is not the star: a planet moves\n"
             "the whole spectrum, spots and\nplage change the lines' shapes.",
             fontsize=8, va="top", color=ps.INK_SECONDARY)

    ax = fig.add_subplot(gs[3, :2])
    th = fit.theta
    y_b = fit.model.residuals(th) + fit.model.keplerian(th, planets=[0])[0]
    err = fit.model.effective_errors(th)
    P, K, e, w, tc = (a[0] for a in fit.model.planet_arrays(th))
    grid = np.linspace(-0.25, 0.75, 300)
    curve = (grid, fit.model.keplerian(th, tc[0] + grid * P[0], planets=[0])[0])
    _binned_fold(ax, data.t, y_b, err, P[0], tc[0], n_bins=10, curve=curve)
    pb = res["fiducial"]["planet"]
    ax.set_ylim(-5.5, 5.5)
    ax.set_xlabel("orbital phase")
    ax.set_ylabel("velocity (m/s)")
    ax.set_title(f"(c) Folded on {pb['P']['median']:.3f} d", fontsize=9.5)
    ax.text(0.02, 0.96, f"K = {pb['K']['median']:.2f} ± {pb['K']['plus']:.2f} m/s\n"
            f"m sin i = {pb['msini_mearth']['median']:.2f} ± {pb['msini_mearth']['plus']:.2f} Earth masses",
            transform=ax.transAxes, va="top", fontsize=8, color=ps.INK)
    ax = fig.add_subplot(gs[3, 2:])
    sp = art["stacked_b"]
    tN = data.t[sp["N"] - 1]
    x = _years(tN)
    ax.fill_between(x, sp["band_lo"], sp["band_hi"], color=ps.CONTEXT, lw=0, label="90% range for a steady planet")
    ax.plot(x, sp["dchi2"], color=ps.C_SIGNAL, lw=1.4, label="observed")
    ax.set_xlabel("year of the last spectrum included")
    ax.set_ylabel("χ² improvement at 11.19 d")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("(d) The signal builds up as data accumulate", fontsize=9.5)
    ax.text(0.99, 0.04, "the 2016 'Pale Red Dot' campaign:\nnightly spectra for 10 weeks",
            transform=ax.transAxes, ha="right", fontsize=7, color=ps.INK_SECONDARY)
    ps.caption(fig, "HARPS RVBank velocities (Trifonov et al. 2020) after removing flares (H-alpha) and "
               "spectra with distorted lines. Panel c: rotation signal (a 91.7-d sinusoid) removed; grey = "
               "single spectra, blue = phase-binned means. m sin i with M★ = 0.1221 ± 0.0022 M☉ (Faria et al. "
               "2022). Band in d: the same signal injected into white noise; correlated activity noise makes "
               "real data wander more.", y=0.04)
    return ps.save(fig, outdir / "06_proxima")


def fig_gj581(res, art, outdir=FIGURE_DIR):
    data = art["data"]
    freqs = art["freqs"]
    three = art["three"]
    fig = plt.figure(figsize=(11.5, 11.2))
    gs = GridSpec(3, 2, hspace=0.55, wspace=0.34, figure=fig)
    P_d = res["signal_d"]["P"]
    rot = res["halpha_rotation"]

    ax = fig.add_subplot(gs[0, 0])
    pw, _, _ = __import__("rv_fit").residual_periodogram(three, freqs)
    per = 1 / freqs
    ax.axvspan(126, 139, color=ps.CONTEXT, lw=0)
    ax.plot(per, pw, color=ps.C_SIGNAL, lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(1, 3000)
    ax.set_ylim(0, 0.3)
    j = np.argmin(np.abs(per - P_d)); j = j - 20 + int(np.argmax(pw[j - 20:j + 20]))
    ps.annotate(ax, f"'GJ 581 d': {per[j]:.1f} d", (per[j], pw[j]), (6, 0.27), relpos=(1, 0.5), ha="right")
    ya = 1 / (1 / P_d - 1 / 365.25)
    ja = np.argmin(np.abs(per - ya))
    ps.annotate(ax, f"its yearly alias {ya:.0f} d\n(first published: 83.4 d)", (per[ja], pw[ja]), (300, 0.2),
                relpos=(0, 0.5))
    ax.text(np.sqrt(126 * 139), 0.29, "rotation", ha="center", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlabel("period (days)")
    ax.set_ylabel("power")
    ax.set_title("(a) Velocities with planets b, c and e removed", fontsize=9.5)

    ax = fig.add_subplot(gs[0, 1])
    ha = art["halpha"]
    per = 1 / ha["freqs"]
    ax.plot(per, ha["power"], color=ps.BLUE, lw=0.8, label="H-alpha index")
    ax.plot(per, ha["power_after"], color=ps.ORANGE, lw=0.8, label="after removing the rotation")
    ax.set_xscale("log")
    ax.set_xlim(2, 3000)
    ax.set_ylim(0, 0.36)
    ax.axvline(P_d, color=ps.INK_SECONDARY, lw=0.8)
    ps.annotate(ax, f"rotation {rot['period']:.0f} d", (rot["period"], rot["power"]), (600, 0.33), relpos=(0, 0.5))
    hp = rot["half_rotation_peak"]
    ps.annotate(ax, f"half-rotation harmonics\n{hp['period']:.0f} and 66.7 d", (hp["period"], hp["power"]),
                (4, 0.26), relpos=(0, 0.5), ha="left")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlabel("period (days)")
    ax.set_ylabel("power")
    ax.set_title("(b) The star's chromosphere (H-alpha): rotation and its harmonic", fontsize=9.5)

    ax = fig.add_subplot(gs[1, 0])
    r3 = three.model.residuals(three.theta)
    hav = data.activity["Halpha"]
    ax.plot(hav, r3, ls="none", marker="o", ms=3, color=ps.BLUE, mec=ps.SURFACE, mew=0.4, alpha=0.8)
    b = np.polyfit(hav, r3, 1)
    xx = np.linspace(hav.min(), hav.max(), 10)
    ax.plot(xx, np.polyval(b, xx), color=ps.C_MODEL, lw=1.3)
    corr = res["rv_halpha_correlation"]
    seasons = ", ".join(f"{k}: {v:+.2f}" for k, v in corr.items() if k != "all")
    ax.set_xlabel("H-alpha index")
    ax.set_ylabel("velocity residual (m/s)")
    ax.set_title(f"(c) More activity, lower velocity: r = {corr['all']:+.2f}", fontsize=9.5)
    ax.text(0.02, 0.03, "r by season:\n" + seasons.replace(", ", "\n"), transform=ax.transAxes, ha="left",
            va="bottom", fontsize=7, color=ps.INK_SECONDARY)

    ax = fig.add_subplot(gs[1, 1])
    display = {"b, c, e": "planets b, c, e", "b, c, e + H-alpha": "+ H-alpha term",
               "b, c, e + H-alpha by season": "+ H-alpha term\nper season", "b, c, e, d": "+ a 4th planet 'd'"}
    labels = [k for k in display if k in res["models"]]
    vals = [res["models"][k]["power_near_66.8d"] for k in labels]
    ypos = np.arange(len(labels))[::-1]
    ax.barh(ypos, vals, height=0.45, color=ps.BLUE)
    for yp, v, k in zip(ypos, vals, labels):
        fap = res["models"][k]["fap_near_66.8d"]
        ftxt = f"{fap:.0e}" if fap < 0.01 else ("≈ 1" if fap > 0.99 else f"{fap:.0%}")
        ax.text(v + 0.004, yp, f"{v:.2f}  (FAP {ftxt})", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_yticks(ypos)
    ax.set_yticklabels([display[k] for k in labels], fontsize=8)
    ax.set_xlim(0, 0.32)
    ax.set_xlabel("periodogram power left at 66.8 d")
    ax.grid(axis="y", visible=False)
    ax.set_title("(d) A linear H-alpha correction removes only a sixth of it", fontsize=9.5)

    ax = fig.add_subplot(gs[2, 0])
    err3 = three.model.effective_errors(three.theta)
    yr = _years(data.t)
    early = yr < 2008.9
    pP, tc0 = P_d, data.t.min()
    ph_fit = np.linspace(-0.25, 0.75, 300)
    ep = res["epochs"]["epochs"]
    _binned_fold(ax, data.t[early], r3[early], err3[early], pP, tc0, n_bins=8, color=ps.BLUE,
                 label=f"{int(yr[early].min())}–{int(yr[early].max())}: {ep[0]['amplitude']:.1f} ± "
                       f"{ep[0]['amplitude_err']:.1f} m/s", shift=-0.01)
    _binned_fold(ax, data.t[~early], r3[~early], err3[~early], pP, tc0, n_bins=8, color=ps.ORANGE,
                 label=f"{int(yr[~early].min())}–{int(yr[~early].max())}: {ep[1]['amplitude']:.1f} ± "
                       f"{ep[1]['amplitude_err']:.1f} m/s", shift=0.01)
    ax.set_ylim(-7, 7)
    ax.legend(fontsize=7, loc="upper right")
    ax.set_xlabel(f"phase at {P_d:.2f} d")
    ax.set_ylabel("velocity residual (m/s)")
    ax.set_title(f"(e) Same amplitude and phase in both halves (p = {res['epochs']['p_consistent']:.2f})",
                 fontsize=9.5)

    ax = fig.add_subplot(gs[2, 1])
    sp = art["stacked_d"]
    x = _years(data.t[sp["N"] - 1])
    ax.fill_between(x, sp["band_lo"], sp["band_hi"], color=ps.CONTEXT, lw=0, label="90% range for a steady planet")
    ax.plot(x, sp["dchi2"], color=ps.C_SIGNAL, lw=1.4, label="observed")
    ax.set_xlabel("year of the last spectrum included")
    ax.set_ylabel(f"χ² improvement at {P_d:.1f} d")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("(f) ...and it builds up like a planet's would", fontsize=9.5)
    fig.suptitle("Gliese 581 'd': a habitable-zone planet that turned out to be the star",
                 x=0.06, ha="left", y=0.95)
    ps.caption(fig, "HARPS RVBank velocities 2004–2012. Circumstantial case for activity (a–d): the period is half "
               f"the star's {rot['period']:.0f}-d rotation (literature 132.5 ± 6.3 d), H-alpha itself varies at the "
               "harmonic, and velocities anticorrelate with H-alpha. HARPS alone cannot close the case (e–f): the signal "
               "stayed coherent over 2004–2012. Robertson et al. (2014) first attributed it to activity; adding new CARMENES "
               "velocities, von Stauffenberg et al. (2024) find only three planets.", y=0.04)
    return ps.save(fig, outdir / "07_gj581_activity")


# ---------------------------------------------------------------------------
# All our planets together
# ---------------------------------------------------------------------------

#: A light step (250) of the blue ramp: velocity-discovered planets, kept
#: recessive behind our own (orange) but distinguishable from the grey rest.
RV_POP = "#86b6ef"

#: Solar System planets: orbital period (days) and mass (Earth masses).
SOLAR_SYSTEM = {"Mercury": (87.97, 0.0553), "Venus": (224.70, 0.815), "Earth": (365.26, 1.0),
                "Mars": (686.98, 0.107), "Jupiter": (4332.6, 317.8), "Saturn": (10759, 95.2),
                "Uranus": (30687, 14.5), "Neptune": (60190, 17.1)}


def our_planets(results) -> list[dict]:
    """Every planet we measured: name, period, minimum (or true) mass, errors."""
    out = []
    pk = results["51peg"]["fits"]["circular"]["planet"]
    out.append(("51 Peg b", pk))
    for p in results["55cnc"]["fiducial"]["planets"][:5]:
        out.append((p["label"], p))
    out.append(("Proxima Cen b", results["proxima"]["fiducial"]["planet"]))
    for p in results["gj581"]["fiducial"]["planets"][:3]:
        out.append((p["label"], p))
    rows = []
    for name, p in out:
        m = p["msini_mearth"]
        rows.append({"name": name, "P": p["P"]["median"], "m": m["median"], "m_lo": m["minus"], "m_hi": m["plus"],
                     "K": p["K"]["median"], "K_err": 0.5 * (p["K"]["plus"] + p["K"]["minus"])})
    return rows


def fig_mass_period(results, outdir=FIGURE_DIR):
    import pandas as pd
    from fetch_data import NASA_POPULATION_PATH
    pop = pd.read_csv(NASA_POPULATION_PATH)
    good = (pop["pl_bmasse"].notna() & pop["pl_orbper"].notna() & (pop["pl_bmasselim"].fillna(0) == 0)
            & (pop["pl_bmassprov"] != "M-R relationship")
            & ~pop["pl_bmasse_reflink"].fillna("").str.contains("CALCULATED_VALUE")
            & (pop["pl_controv_flag"].fillna(0) == 0))
    pop = pop[good]
    e_pl = next(p for p in results["55cnc"]["fiducial"]["planets"] if p["label"].endswith(" e"))
    true_over_min = results["55cnc"]["tess"]["true_mass_earth"]["median"] / e_pl["msini_mearth"]["median"]
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    rv = pop["discoverymethod"] == "Radial Velocity"
    ax.plot(pop.loc[~rv, "pl_orbper"], pop.loc[~rv, "pl_bmasse"], ls="none", marker="o", ms=2.4,
            color=ps.CONTEXT, label=f"other methods ({(~rv).sum():,})")
    ax.plot(pop.loc[rv, "pl_orbper"], pop.loc[rv, "pl_bmasse"], ls="none", marker="o", ms=2.4,
            color=RV_POP, label=f"found by radial velocity ({rv.sum():,})")
    # what a given K means around a Sun-like star
    P = np.logspace(-0.7, 5.3, 200)
    for K in (0.1, 1.0, 10.0, 100.0):
        m = phys.msini_mearth(P, K, 0.0, 1.0)
        ax.plot(P, m, color=ps.MUTED, lw=0.7)
        i = np.argmin(np.abs(P - 3e4))
        ax.text(P[i], m[i] * 1.12, f"K = {K:g} m/s", fontsize=7, color=ps.INK_SECONDARY, rotation=10,
                ha="center", va="bottom")
    for name, (Pd, m) in SOLAR_SYSTEM.items():
        ax.plot(Pd, m, marker="o", ms=5, color=ps.INK, ls="none")
        ax.text(Pd * 1.12, m, name, fontsize=7, va="center", color=ps.INK)
    rows = our_planets(results)
    for r in rows:
        ax.errorbar(r["P"], r["m"], yerr=[[r["m_lo"]], [r["m_hi"]]], fmt="o", ms=7, color=ps.C_OURS,
                    mec=ps.SURFACE, mew=1.2, elinewidth=1.2, zorder=6)
    offsets = {"51 Peg b": (1.35, 1.0), "55 Cnc b": (1.3, 1.3), "55 Cnc c": (1.25, 0.8), "55 Cnc d": (1.25, 1.1),
               "55 Cnc e": (0.8, 1.35), "55 Cnc f": (1.25, 0.85), "Proxima Cen b": (1.25, 0.72),
               "GJ 581 b": (1.2, 1.25), "GJ 581 c": (1.25, 0.75), "GJ 581 e": (0.8, 0.6)}
    for r in rows:
        dx, dy = offsets.get(r["name"], (1.2, 1.0))
        ax.text(r["P"] * dx, r["m"] * dy, r["name"], fontsize=7.5, color=ps.INK, fontweight="bold",
                ha="left" if dx > 1 else "right", va="center", zorder=7)
    gd = results["gj581"]["signal_d"]
    md = phys.msini_mearth(gd["P"], gd["K"], 0.0, 0.295)
    ax.plot(gd["P"], md, marker="o", ms=7, mfc="none", mec=ps.C_OURS, mew=1.4, ls="none", zorder=6)
    ax.text(gd["P"] * 1.2, md * 1.05, "GJ 581 'd': stellar activity", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.2, 2.5e5)
    ax.set_ylim(0.03, 2e4)
    ax.set_xlabel("orbital period (days)")
    ax.set_ylabel("planet mass (Earth masses); minimum mass m sin i for velocity planets")
    ax.set_title(f"Our ten planets among the {len(pop):,} known planets with a measured mass", pad=10)
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=7, color=ps.C_OURS, label="measured here"),
                       Line2D([], [], ls="none", marker="o", ms=4, color=RV_POP,
                              label=f"found by radial velocity ({rv.sum():,})"),
                       Line2D([], [], ls="none", marker="o", ms=4, color=ps.CONTEXT,
                              label=f"found by other methods ({(~rv).sum():,})"),
                       Line2D([], [], ls="none", marker="o", ms=5, color=ps.INK, label="Solar System")],
              loc="upper left", fontsize=8)
    ps.caption(fig, "Background: NASA Exoplanet Archive composite table, masses measured (not estimated from radii), "
                    "controversial planets excluded. Grey lines: minimum mass giving a stellar wobble K around a "
                    "Sun-like star; our masses use each host's own mass. HARPS reaches ~1 m/s and ESPRESSO a few tenths, "
                    "so an Earth twin (0.09 m/s) is still out of reach. Our points are m sin i (for the "
                    f"transiting 55\u00a0Cnc\u00a0e the true mass is only {true_over_min - 1:.0%} higher).")
    return ps.save(fig, outdir / "08_mass_period")


def fig_comparison(results, outdir=FIGURE_DIR):
    """Our semi-amplitudes and minimum masses against every published value."""
    import analysis
    rows = our_planets(results)
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.8), sharey=True, gridspec_kw={"wspace": 0.08})
    ypos = np.arange(len(rows))[::-1]
    for ax, (key, our_key, title) in zip(axes, [("pl_rvamp", "K", "semi-amplitude K: published ÷ ours"),
                                                ("pl_msinie", "m", "minimum mass: published ÷ ours")]):
        for y, r in zip(ypos, rows):
            pub = analysis.literature_rows(r["name"])
            ours = r[our_key]
            err = (r["K_err"] / r["K"]) if our_key == "K" else 0.5 * (r["m_lo"] + r["m_hi"]) / r["m"]
            ax.fill_betweenx([y - 0.28, y + 0.28], 1 - err, 1 + err, color=ps.C_OURS, alpha=0.25, lw=0)
            ax.plot([1, 1], [y - 0.28, y + 0.28], color=ps.C_OURS, lw=2)
            for p in pub:
                v = p.get(key)
                if v is None or not np.isfinite(v):
                    continue
                alias = r["name"] == "55 Cnc e" and (p.get("pl_orbper") or 0) > 1.5
                ratio = v / ours
                ax.plot(np.clip(ratio, 0.45, 2.2), y + (0.12 if p["default"] else 0), marker="o", ms=5,
                        mfc="none" if alias else ps.C_PUBLISHED, mec=ps.C_PUBLISHED, mew=1.0, ls="none", alpha=0.85)
        ax.set_xscale("log")
        ax.set_xlim(0.45, 2.2)
        ax.set_xticks([0.5, 0.67, 0.8, 1, 1.25, 1.5, 2])
        ax.set_xticklabels(["0.5", "0.67", "0.8", "1", "1.25", "1.5", "2"])
        ax.xaxis.set_minor_locator(mpl_ticker.NullLocator())
        ax.set_title(title, fontsize=10)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels([f"{r['name']}   K = {r['K']:.2f} m/s" for r in rows], fontsize=8)
    axes[0].text(0.47, ypos[[r["name"] for r in rows].index("55 Cnc e")] - 0.32,
                 "open circles: the 2.8-d alias era (McArthur 2004, Wright 2009)", fontsize=7, ha="left",
                 color=ps.INK_SECONDARY, va="top")
    axes[0].legend(handles=[Line2D([], [], color=ps.C_OURS, lw=2, label="ours (band: 1σ)"),
                            Line2D([], [], ls="none", marker="o", ms=5, color=ps.C_PUBLISHED,
                                   label="each published solution")],
                   loc="lower left", fontsize=7.5)
    fig.suptitle("Our results against every published solution in the NASA Exoplanet Archive", x=0.08, ha="left")
    ps.caption(fig, "Published minimum masses also differ because each paper adopts its own stellar mass (55 Cnc: "
               "0.85–1.02 M☉ across papers); K is the direct measurement. 55 Cnc d depends on how its star's "
               "magnetic cycle is modelled (see README). Points beyond the axis are drawn at its edge.", y=0.0)
    return ps.save(fig, outdir / "00_comparison")
