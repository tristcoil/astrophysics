"""Every figure, drawn from the cached analysis (no network, no refitting).

    python figures.py        # after run.py has produced data/results.json
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

import background as B
import catalogs as K
import config as C
import deltanu as D
import lightcurves as L
import literature as LIT
import period_spacing as P
import plotstyle as ps
import spectrum as S
import stars

ps.use_style()
FIG = C.FIGURE_DIR

#: Single-hue sequential ramp (the palette's blue, light to dark), starting
#: at the chart surface so that "no power" recedes into the background.
BLUE_RAMP = LinearSegmentedColormap.from_list(
    "blue_ramp", [ps.SURFACE, "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

#: Neutral ramp for background densities (other people's samples), so that
#: our own stars keep the categorical colours.
GREY_RAMP = LinearSegmentedColormap.from_list("grey_ramp", [ps.SURFACE, "#d6d5cf", "#a9a8a1", "#6f6e69"])
#: Single-hue ramp for markers coloured by a magnitude (mass): starts at the
#: palette's step 250 so the lightest marker still stands off the surface.
MASS_RAMP = LinearSegmentedColormap.from_list("mass_ramp", ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

C_RGB = ps.BLUE  # hydrogen-shell burning (red giant branch)
C_HEB = ps.ORANGE  # helium-core burning (red clump / secondary clump)
C_THIRD = ps.AQUA
STATE_LABEL = {"RGB": "Red giant branch (shell H burning)", "RC": "Red clump (core He burning)"}

#: Showcase stars for the single-star figures
HERO = 9145955  # bright (Kp 9.8) red giant branch star, Vrard et al. 2016's Fig. 3 star
RGB_STAR = 6779699  # Bedding et al. 2011's red giant branch star
RC_STAR = 12645107  # a red-clump star from the random sample with clean ridges (DPi1 = Vrard's)


def _results() -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(C.DATA_DIR / "results_stars.csv", index_col="kic")
    res = json.loads(C.RESULTS_PATH.read_text())
    return df, res


def _plain_log(ax, axis="both", ticks=None):
    """Log axis with plain-number tick labels (4, 6, 10, 20 ...), not 4x10^0."""
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator, ScalarFormatter

    for a in ([ax.xaxis, ax.yaxis] if axis == "both" else [getattr(ax, f"{axis}axis")]):
        if ticks is not None:
            a.set_major_locator(FixedLocator(ticks))
        fmt = ScalarFormatter()
        fmt.set_scientific(False)
        a.set_major_formatter(fmt)
        a.set_minor_formatter(NullFormatter())
        a.set_minor_locator(NullLocator())


def _state_color(state: str) -> str:
    return C_HEB if state == "RC" else C_RGB


def _snr(arr):
    return arr["psd"] / arr["bg"]


def _echelle(ax, nu, snr, dnu, lo, hi, *, x_offset=0.0, smooth_frac=0.004, vmax_pct=99.0,
             nx=180, cmap=BLUE_RAMP):
    """Draw an echelle image; x in uHz (nu mod Dnu, shifted by x_offset Dnu)."""
    df = nu[1] - nu[0]
    sm = S.smooth_boxcar(snr, max(1, int(round(smooth_frac * dnu / df))))
    xe, rows, img = D.echelle_image(nu, sm, dnu, lo, hi, x_offset=x_offset, nx=nx)
    vmax = np.percentile(img[img > 0], vmax_pct) if (img > 0).any() else 1
    ax.imshow(img, aspect="auto", origin="lower", cmap=cmap, vmin=0, vmax=vmax,
              extent=[xe[0], xe[-1], rows[0], rows[-1] + dnu], interpolation="nearest")
    ax.set_xlim(xe[0], xe[-1])
    ax.set_ylim(rows[0], rows[-1] + dnu)
    ax.grid(False)
    return xe, rows


# ---------------------------------------------------------------------------
# 1. From starlight to a power spectrum
# ---------------------------------------------------------------------------


def fig_light_to_spectrum(df, res):
    kic = HERO
    js, arr = stars.load_star(kic)
    ts = L.kepler_long_cadence(kic)
    numax, dnu = js["numax"], js["dnu"]
    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.25], hspace=0.55, wspace=0.22)

    ax = fig.add_subplot(gs[0, :])
    t = ts.time
    ax.plot(t, ts.flux, ",", color=ps.MUTED, alpha=0.6, rasterized=True)
    # one-day means, to show there is no visible trend left
    nb = np.floor(t).astype(int)
    u, inv = np.unique(nb, return_inverse=True)
    m = np.bincount(inv, ts.flux) / np.bincount(inv)
    ax.plot(u + 0.5, m, color=ps.BLUE, lw=0.8)
    ax.set_xlabel("Time (days since 2009 January 1, barycentric)")
    ax.set_ylabel("Brightness (ppm)")
    ax.set_ylim(np.percentile(ts.flux, 0.2) * 1.4, np.percentile(ts.flux, 99.8) * 1.4)
    ax.set_title(f"a. Four years of Kepler photometry of KIC {kic}: {len(t):,} points, "
                 f"each 29.4 minutes long")

    ax = fig.add_subplot(gs[1, 0])
    # the oscillations alone: a Gaussian band-pass around the hump (as in sonify.py)
    idx = np.round((t - t[0]) * 86400.0 / ts.cadence_s).astype(int)
    grid = np.zeros(idx.max() + 1)
    grid[idx] = ts.flux
    f = np.fft.rfftfreq(len(grid), ts.cadence_s) * 1e6
    sig = js["envelope_fwhm"] / 2.3548
    bp = np.fft.irfft(np.fft.rfft(grid) * np.exp(-0.5 * ((f - numax) / sig) ** 2), n=len(grid))[idx]
    t0 = 700.0
    m = (t > t0) & (t < t0 + 2.5)
    ax.plot((t[m] - t0) * 24, ts.flux[m], "o", ms=2.2, color=ps.MUTED, label="Measured")
    ax.plot((t[m] - t0) * 24, bp[m], "-", lw=1.3, color=ps.BLUE, label="Band-passed around νmax")
    ax.set_xlabel("Hours")
    ax.set_ylabel("Brightness (ppm)")
    period_h = 1e6 / numax / 3600
    ax.set_title(f"b. Sixty hours: ringing every ~{period_h:.1f} h")
    ax.legend(loc="lower right", fontsize=7, ncol=2)

    ax = fig.add_subplot(gs[1, 1])
    nu_b, p_b = arr["bin_nu"], arr["bin_psd"]
    th = arr["theta_post"]
    eta2 = S.apodization(nu_b, C.NYQUIST_LC_UHZ)
    tr = L.filter_transfer(nu_b, C.TREND_SIGMA_DAYS_LC)
    comp = B.components(th, nu_b, eta2, tr)
    ax.loglog(nu_b, p_b, color=ps.MUTED, lw=0.6, label="Power spectrum (log-binned)")
    ax.loglog(nu_b, comp["gran1"] + comp["gran2"] + comp["slow"], color=ps.INK_SECONDARY, lw=1.0,
              label="Granulation + slow variability")
    ax.loglog(nu_b, comp["white"], color=ps.INK_SECONDARY, lw=0.8, ls=(0, (1, 1.5)),
              label="Photon noise")
    ax.loglog(nu_b, comp["envelope"], color=ps.ORANGE, lw=1.2, label="Oscillation hump")
    ax.loglog(nu_b, sum(comp.values()), color=ps.BLUE, lw=1.4, label="Total model")
    lo = np.percentile(p_b, 0.5)
    ax.set_ylim(lo * 5e-3, p_b.max() * 3)
    ax.set_xlim(nu_b[0], C.NYQUIST_LC_UHZ)
    ps.annotate(ax, f"νmax = {numax:.1f} µHz", (numax, comp["envelope"].max() * 1.3),
                (numax * 0.25, p_b.max() * 1.2))
    ax.set_xlabel("Frequency (µHz)")
    ax.set_ylabel("Power density (ppm²/µHz)")
    ax.set_title("c. The power spectrum: background plus a hump")
    ax.legend(loc="lower left", fontsize=6.5, ncol=2)

    ax = fig.add_subplot(gs[2, :])
    nu, snr = arr["nu"], _snr(arr)
    m = (nu > numax - 3.3 * dnu) & (nu < numax + 3.3 * dnu)
    ax.plot(nu[m], S.smooth_boxcar(snr[m], 3), color=ps.BLUE, lw=0.6)
    ax.set_xlim(numax - 3.3 * dnu, numax + 3.3 * dnu)
    top = np.percentile(S.smooth_boxcar(snr[m], 3), 99.95)
    ax.set_ylim(0, top * 1.25)
    eps = js["epsilon"]
    up = D.universal_pattern(dnu)
    n_lo = int(np.floor((numax - 3.3 * dnu) / dnu - eps))
    for n in range(n_lo, n_lo + 8):
        f0 = (n + eps) * dnu
        if numax - 3.3 * dnu < f0 < numax + 3.3 * dnu:
            ax.text(f0, top * 1.12, "ℓ=0", ha="center", fontsize=7, color=ps.INK_SECONDARY)
        f2 = (n + eps - up["d02"]) * dnu
        if numax - 3.3 * dnu < f2 < numax + 3.3 * dnu:
            ax.text(f2, top * 1.03, "2", ha="center", fontsize=7, color=ps.MUTED)
        f1 = (n + 0.5 + eps - up["d01"]) * dnu
        if numax - 3.3 * dnu < f1 < numax + 3.3 * dnu:
            ax.text(f1, top * 1.12, "ℓ=1", ha="center", fontsize=7, color=ps.INK_SECONDARY)
    f0 = (round(numax / dnu - eps) + eps) * dnu
    ax.annotate("", xy=(f0 + dnu, top * 0.92), xytext=(f0, top * 0.92),
                arrowprops={"arrowstyle": "<->", "color": ps.ORANGE, "lw": 1.0})
    ax.text(f0 + dnu / 2, top * 0.95, f"Δν = {dnu:.2f} µHz", ha="center", fontsize=8,
            color=ps.INK)
    ax.set_xlabel("Frequency (µHz)")
    ax.set_ylabel("Signal / background")
    ax.set_title("d. Inside the hump: a comb of modes, the same pattern repeating every Δν")
    ps.caption(fig, f"KIC {kic}, Kepler long cadence Q0-Q17 (PDCSAP, MAST). Spectrum calibrated so that "
                     "its integral equals the variance of the light curve (Parseval). ℓ is the angular "
                     "degree of the mode: 0 = the whole star breathing, 1 = sloshing, 2 = quadrupole.")
    return ps.save(fig, FIG / "01_light_to_spectrum")


# ---------------------------------------------------------------------------
# 2. The ladder: big stars ring at low notes
# ---------------------------------------------------------------------------


def _ladder_stars(df, n=9):
    """One star per numax bin, the one with the strongest envelope detection."""
    out = []
    edges = np.asarray(C.NUMAX_BIN_EDGES)
    for i in range(len(edges) - 1):
        m = (df["numax"] >= edges[i]) & (df["numax"] < edges[i + 1])
        if m.any():
            out.append(df[m]["delta_lnL"].idxmax())
    return sorted(out, key=lambda k: -df.loc[k, "numax"])


def fig_ladder(df, res):
    kics = _ladder_stars(df)
    fig, ax = plt.subplots(figsize=(10, 8.5))
    for i, kic in enumerate(kics):
        js, arr = stars.load_star(int(kic))
        nu_b, p_b = arr["bin_nu"], arr["bin_psd"]
        th = arr["theta_post"]
        eta2 = S.apodization(nu_b, C.NYQUIST_LC_UHZ)
        tr = L.filter_transfer(nu_b, C.TREND_SIGMA_DAYS_LC)
        bg = B.background_only(th, nu_b, eta2, tr)
        # smooth by ~3% in frequency on a uniform log grid: keeps the hump, drops the chi^2 noise
        lg = np.arange(np.log(1.0), np.log(C.NYQUIST_LC_UHZ), 0.004)
        yy = np.interp(lg, np.log(nu_b), p_b / bg)
        from scipy.ndimage import gaussian_filter1d

        y = np.log10(np.maximum(gaussian_filter1d(yy, 8, mode="nearest"), 0.5))
        off = (len(kics) - 1 - i) * 1.0
        col = _state_color(df.loc[kic, "state_apokasc"])
        ax.plot(np.exp(lg), y + off, color=col, lw=1.1)
        ax.axhline(off, color=ps.GRID, lw=0.5, zorder=0)
        r = df.loc[kic, "R_self"]
        ax.text(1.05, off + 0.12, f"KIC {kic}   νmax {js['numax']:.0f} µHz   R = {r:.1f} R☉",
                fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_xlim(1, C.NYQUIST_LC_UHZ)
    ax.set_yticks([])
    ax.set_xlabel("Frequency (µHz)")
    ax.set_ylabel("Power / background (log scale, smoothed; stars offset)")
    ax.set_title(f"The bigger the star, the lower the note: oscillation humps of {len(kics)} red giants")
    from matplotlib.lines import Line2D

    ax.legend(handles=[Line2D([], [], color=C_RGB, label=STATE_LABEL["RGB"]),
                       Line2D([], [], color=C_HEB, label=STATE_LABEL["RC"])],
              loc="upper right", fontsize=7.5)
    ps.caption(fig, "One star from each νmax bin of the sample (the strongest detection). Each spectrum is "
                     "divided by its fitted granulation background, so the hump of oscillation power stands "
                     "out. Radii from our scaling relations, calibrated on our own measurement of the Sun. "
                     "Evolutionary states from APOKASC-3.")
    return ps.save(fig, FIG / "02_ladder")


# ---------------------------------------------------------------------------
# 3. Echelle diagrams
# ---------------------------------------------------------------------------


def _echelle_panel(ax, nu, snr, dnu, numax, sig, eps, title, *, n_sig=2.8, labels=True,
                   lund=None, smooth_frac=0.004, legend=True):
    up = D.universal_pattern(dnu)
    # put the l=0 ridge at 0.3 Dnu from the left edge so no ridge is split
    x_off = (eps - 0.3) % 1.0
    lo, hi = numax - n_sig * sig, numax + n_sig * sig
    xe, rows = _echelle(ax, nu, snr, dnu, lo, hi, x_offset=x_off, smooth_frac=smooth_frac)
    ax.set_xlabel(f"Frequency modulo Δν = {dnu:.3f} µHz (µHz)")
    ax.set_ylabel("Frequency (µHz)")
    ax.set_title(title)
    if labels:
        y_top = rows[-1] + dnu
        for name, x, ha in (("ℓ=0", 0.3, "left"), ("ℓ=2", 0.3 - up["d02"], "right"),
                            ("ℓ=1", 0.8 - up["d01"], "center")):
            xx = ((x % 1.0) + x_off) * dnu
            if xx > xe[-1]:
                xx -= dnu
            ax.text(xx, y_top - 0.015 * (y_top - rows[0]), name, ha=ha, va="top", fontsize=8,
                    color=ps.INK, bbox={"boxstyle": "round,pad=0.15", "fc": ps.SURFACE, "ec": "none",
                                        "alpha": 0.85})
    if lund is not None:
        markers = {0: "o", 1: "^", 2: "s", 3: "D"}
        for l, g in lund.groupby("l"):
            f = g["freq"].to_numpy()
            x = (f / dnu - x_off) % 1.0 * dnu + x_off * dnu
            ax.plot(x, f, markers[int(l)], ms=4, mfc="none", mec=ps.ORANGE, mew=0.9,
                    label=f"ℓ={int(l)}, Lund et al. 2017")
        if legend:
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), fontsize=7, ncol=2)


def fig_echelle(df, res):
    fig, axs = plt.subplots(1, 3, figsize=(13, 6.6))
    for ax, kic, kind, sm in ((axs[0], HERO, "a. Red giant branch", 0.004),
                              (axs[1], RC_STAR, "b. Red clump", 0.01)):
        js, arr = stars.load_star(kic)
        _echelle_panel(ax, arr["nu"], _snr(arr), js["dnu"], js["numax"], js["envelope_fwhm"] / 2.3548,
                       js["epsilon"], f"{kind}: KIC {kic}", smooth_frac=sm)
    ms = json.loads((C.CACHE_DIR / "sixteen_cyg.json").read_text())["16 Cyg A"]
    z = np.load(C.CACHE_DIR / "sixteen_cyg.npz")
    _, t6 = K.lund2017()
    lund = t6[t6["kic"] == "12069424"]
    _echelle_panel(axs[2], z["16_cyg_a_nu"], z["16_cyg_a_psd"] / z["16_cyg_a_bg"], ms["dnu"],
                   ms["numax"], ms["envelope_fwhm"] / 2.3548, ms["epsilon"],
                   "c. Main sequence: 16 Cyg A", labels=False, lund=lund, n_sig=1.6, smooth_frac=0.01)
    fig.suptitle("Échelle diagrams: cut the spectrum into Δν-wide strips and stack them", y=1.0)
    ps.caption(fig, "Each row is one Δν-wide slice of the spectrum (signal / background, lightly smoothed); "
                     "darker = more power. Radial (ℓ=0) and quadrupole (ℓ=2) modes form straight ridges. "
                     "In red giants the ℓ=1 ridge breaks into many 'mixed' modes, sound waves coupled to "
                     "gravity waves in the core (Fig. 8). The Sun-like 16 Cyg A has pure sound waves: "
                     "orange symbols are Lund et al. (2017)'s individually fitted frequencies.", y=-0.07)
    return ps.save(fig, FIG / "03_echelle")


# ---------------------------------------------------------------------------
# 4. Global parameters against the catalogues
# ---------------------------------------------------------------------------


def fig_vs_catalogues(df, res):
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1.5, 1.1], hspace=0.55, wspace=0.35)
    g = res["global"]
    letters = iter("abcdef")
    for row, q, name in ((0, "numax", "νmax"), (1, "dnu", "Δν")):
        for col, cat, label in ((0, "ak", "APOKASC-3"), (1, "yu", "Yu et al. 2018")):
            ax = fig.add_subplot(gs[row, col])
            r = 100 * (df[q] / df[f"{q}_{cat}"] - 1)
            err = 100 * np.hypot(df[f"{q}_err"] / df[q], df[f"{q}_err_{cat}"] / df[f"{q}_{cat}"])
            for st, mk in (("RGB", "o"), ("RC", "s")):
                m = df["state_apokasc"] == st
                c = _state_color(st)
                ax.errorbar(df["numax"][m], r[m], yerr=err[m], fmt=mk, ms=4, color=c,
                            mfc=c if st == "RGB" else ps.SURFACE, mew=1, elinewidth=0.6, alpha=0.9,
                            label=STATE_LABEL[st].split(" (")[0])
            ax.axhline(0, color=ps.INK_SECONDARY, lw=0.8)
            st = g[f"{q}_vs_{cat}"]
            lim = 6 if q == "numax" else 3
            ax.set_ylim(-lim, lim)
            ax.set_xscale("log")
            ax.set_xlabel("νmax (µHz)")
            ax.set_ylabel(f"{name}: ours / {label.split(' et al')[0]} − 1 (%)")
            ax.set_title(f"{next(letters)}. {name} vs {label}\nmedian {100 * st['median']:+.2f}%, "
                         f"robust scatter {100 * st['scatter']:.2f}%", fontsize=9.5)
            if row == 0 and col == 0:
                ax.legend(loc="lower left", fontsize=7)
        # where we sit among APOKASC-3's ten pipelines
        ax = fig.add_subplot(gs[row, 2])
        items = [(p, v[q]) for p, v in g["pipelines"].items() if q in v and v[q].get("n", 0) >= 10]
        items.sort(key=lambda t: t[1]["scatter"], reverse=True)
        for i, (p, stt) in enumerate(items):
            c = ps.INK if p == "this work" else ps.MUTED
            ax.errorbar(100 * stt["median"], i, xerr=100 * stt["scatter"], fmt="o", color=c, ms=5,
                        elinewidth=1.2 if p == "this work" else 0.9)
        ax.set_yticks(range(len(items)))
        ax.set_yticklabels([p for p, _ in items], fontsize=7.5)
        for t in ax.get_yticklabels():
            if t.get_text() == "this work":
                t.set_fontweight("bold")
                t.set_color(ps.INK)
        ax.axvline(0, color=ps.INK_SECONDARY, lw=0.8)
        ax.set_xlim((-4, 5) if q == "numax" else (-2.2, 2.2))
        ax.set_xlabel(f"{name}: pipeline / APOKASC-3 − 1 (%)")
        ax.set_title(f"{next(letters)}. Among APOKASC-3's pipelines\n(sorted by scatter)", fontsize=9.5)
        ax.grid(axis="y", visible=False)
    ps.caption(fig, "Our blind measurements of the 104 red giants compared star by star with the two catalogues "
                     "(bars: our and their uncertainties combined). Right: median offset (dot) and robust scatter "
                     "(bar) of each of APOKASC-3's ten pipelines relative to their consensus, on the same stars, "
                     "with ours in black; FLI and HON are coarse νmax-only methods and run off the scale. "
                     "Blue filled: red giant branch, orange open: red clump (APOKASC-3 states).")
    return ps.save(fig, FIG / "04_vs_catalogues")


# ---------------------------------------------------------------------------
# 5. Dnu-numax
# ---------------------------------------------------------------------------


def fig_dnu_numax(df, res):
    y = K.yu2018()
    sun = res["sun"]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.hexbin(y["numax"], y["dnu"], xscale="log", yscale="log", gridsize=90, cmap=GREY_RAMP,
              bins="log", mincnt=1, linewidths=0)
    # lines of constant mass from the scaling relations, at a typical giant's 4800 K:
    # Dnu = Dnu_sun (numax/numax_sun)^(3/4) M^(-1/4) (Teff/Teff_sun)^(3/8)
    xx = np.logspace(np.log10(4), np.log10(260), 60)
    tt = (4800.0 / C.TEFF_SUN_NOMINAL) ** 0.375
    for mass in (0.8, 1.2, 2.0, 3.0):
        yy = sun["dnu_sun"] * (xx / sun["numax_sun"]) ** 0.75 * mass ** -0.25 * tt
        ax.plot(xx, yy, color=ps.MUTED, lw=0.7, zorder=1)
        ax.text(xx[-1] * 1.03, yy[-1], f"{mass:g} M☉", fontsize=7, color=ps.INK_SECONDARY, va="center")
    for st, mk in (("RGB", "o"), ("RC", "s")):
        m = df["state_apokasc"] == st
        ax.plot(df["numax"][m], df["dnu"][m], mk, ms=5, color=_state_color(st), mec=ps.SURFACE,
                mew=0.7, label=f"This work: {STATE_LABEL[st].split(' (')[0].lower()}", zorder=3)
    f = res["dnu_numax"]["this work"]
    fy = res["dnu_numax"]["Yu et al. 2018 catalogue (our fit)"]
    ax.plot(xx, f["alpha"] * xx ** f["beta"], color=ps.INK, lw=1.1, zorder=2,
            label=f"Fit to our 104: Δν = {f['alpha']:.3f} νmax^{f['beta']:.3f} (±{f['beta_err']:.3f})")
    st = LIT.DNU_NUMAX["Stello et al. 2009"]
    ax.plot(xx, st["alpha"] * xx ** st["beta"], color=ps.INK_SECONDARY, lw=1.0, ls=(0, (4, 2)), zorder=2,
            label=f"Stello et al. 2009 (dwarfs to giants): {st['alpha']} νmax^{st['beta']}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(4, 330)
    ax.set_xlabel("νmax (µHz)")
    ax.set_ylabel("Δν (µHz)")
    ax.set_title(f"Δν against νmax: our 104 giants on Yu et al.'s {len(y):,} (grey)")
    ax.legend(loc="upper left", fontsize=7.5)
    ps.caption(fig, f"A fit to Yu et al.'s catalogue gives Δν = {fy['alpha']:.3f} νmax^{fy['beta']:.3f}, "
                     "like ours; Stello et al.'s steeper exponent comes from a fit that also spans main-sequence "
                     "stars. Thin grey lines: where stars of 0.8, 1.2, 2 and 3 M☉ at 4800 K fall according to the "
                     "scaling relations (Δν ∝ νmax^¾ M^−¼ Teff^⅜): the width of the band is mostly a spread in mass.")
    return ps.save(fig, FIG / "05_dnu_numax")


# ---------------------------------------------------------------------------
# 6. Masses and radii
# ---------------------------------------------------------------------------


def fig_masses_radii(df, res):
    fig, axs = plt.subplots(2, 2, figsize=(11, 10))
    fig.subplots_adjust(hspace=0.42, wspace=0.28)
    mr = res["masses_radii"]
    for ax, ours, theirs, label, lim, title in (
            (axs[0, 0], "M_ak", "mass_ak", "APOKASC-3 mass (M☉)", (0.6, 3.2),
             "a. Mass, on APOKASC-3's scale"),
            (axs[0, 1], "R_ak", "radius_ak", "APOKASC-3 radius (R☉)", (3, 45),
             "b. Radius, on APOKASC-3's scale"),
            (axs[1, 0], "M_yu", "mass_nocorr_yu", "Yu et al. mass, uncorrected (M☉)", (0.6, 3.2),
             "c. Mass, on Yu et al.'s scale")):
        for st, mk in (("RGB", "o"), ("RC", "s")):
            m = df["state_apokasc"] == st
            ax.errorbar(df[theirs][m], df[ours][m], yerr=df[f"{ours}_err"][m], fmt=mk, ms=4,
                        color=_state_color(st), mec=ps.SURFACE, mew=0.5, elinewidth=0.5, alpha=0.9,
                        label=STATE_LABEL[st].split(" (")[0])
        ax.plot(lim, lim, color=ps.INK_SECONDARY, lw=0.8)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        if "R_" in ours:
            ax.set_xscale("log")
            ax.set_yscale("log")
            _plain_log(ax, "both", [3, 4, 6, 10, 20, 30, 40])
        ax.set_xlabel(label)
        ax.set_ylabel("This work (" + ("M☉" if ours.startswith("M") else "R☉") + ")")
        key = "ak_scale_vs_apokasc" if "ak" in ours else "yu_scale_vs_yu"
        q = "mass" if ours.startswith("M") else "radius"
        s = mr[key][q]
        ax.set_title(f"{title}\nmedian {100 * s['median']:+.1f}%, scatter {100 * s['scatter']:.1f}%",
                     fontsize=9)
        ax.legend(loc="upper left", fontsize=7)
    ax = axs[1, 1]
    for st, mk in (("RGB", "o"), ("RC", "s")):
        m = df["state_apokasc"] == st
        ax.errorbar(df["M_self"][m], df["R_self"][m], xerr=df["M_self_err"][m], yerr=df["R_self_err"][m],
                    fmt=mk, ms=4, color=_state_color(st), mec=ps.SURFACE, mew=0.5, elinewidth=0.5,
                    label=STATE_LABEL[st])
    ax.set_yscale("log")
    _plain_log(ax, "y", [4, 6, 10, 20, 30])
    ax.set_xlabel("Mass (M☉)")
    ax.set_ylabel("Radius (R☉)")
    ax.set_title("d. 104 red giants weighed and measured by their sound\n(self-calibrated on our Sun)",
                 fontsize=9)
    ax.legend(loc="upper right", fontsize=7)
    ps.caption(fig, "a, b: our νmax and Δν with APOKASC-3's solar values (3076, 135.1416 µHz), their per-star "
                     "corrections f_νmax and f_Δν and APOGEE temperatures. c: with Yu et al.'s solar values "
                     "(3090, 135.1 µHz) and temperatures, no corrections. d: our own Sun as the reference, no "
                     "corrections. Error bars propagate νmax, Δν and Teff.")
    return ps.save(fig, FIG / "06_masses_radii")


# ---------------------------------------------------------------------------
# 7. Gaia radii
# ---------------------------------------------------------------------------


def fig_gaia(df, res):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.8), sharey=True)
    gr = res["gaia_radii"]
    titles = {"self": "a. Self-calibrated (our Sun), no corrections",
              "yu": "b. Yu et al.'s solar values, no corrections",
              "ak": "c. APOKASC-3 scale (with corrections)"}
    for ax, s in zip(axs, ("self", "yu", "ak")):
        r = df[f"R_{s}"] / df["radius_gaia_ak"]
        for st, mk in (("RGB", "o"), ("RC", "s")):
            m = df["state_apokasc"] == st
            ax.plot(df["numax"][m], r[m], mk, ms=4.5, color=_state_color(st), mec=ps.SURFACE, mew=0.5,
                    label=STATE_LABEL[st].split(" (")[0] +
                    f": median {gr[s][st]['median']:.3f}")
        ax.axhline(1, color=ps.INK_SECONDARY, lw=0.8)
        for st in ("RGB", "RC"):
            ax.axhline(gr[s][st]["median"], color=_state_color(st), lw=1.0, alpha=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("νmax (µHz)")
        ax.set_title(titles[s], fontsize=9)
        ax.legend(loc="lower right", fontsize=7)
        ax.set_ylim(0.8, 1.2)
    axs[0].set_ylabel("Seismic radius / Gaia radius")
    a, b = gr["self"], gr["yu"]
    ps.caption(fig, "Gaia radii (APOKASC-3) come from Gaia DR3 parallaxes, 2MASS K magnitudes and APOGEE "
                     "temperatures: no seismology. Coloured lines: medians. With our own measurement of the Sun as "
                     f"the reference, uncorrected radii are {100 * (a['RGB']['median'] - 1):+.1f}% (red giant "
                     f"branch) and {100 * (a['RC']['median'] - 1):+.1f}% (clump) off; with Yu et al.'s solar values "
                     f"{100 * (b['RGB']['median'] - 1):+.1f}% and {100 * (b['RC']['median'] - 1):+.1f}%. The split "
                     "between the two groups is the mean-density correction f_Δν from stellar models; APOKASC-3's "
                     "corrections, calibrated on Gaia radii, remove both offsets.")
    return ps.save(fig, FIG / "07_gaia_radii")


# ---------------------------------------------------------------------------
# 8. Period spacings
# ---------------------------------------------------------------------------


def _stretched_points(arr, js, dpi1, q, thresh=8.0):
    """One point per dipole-mode peak: (tau mod Delta Pi_1, frequency, height)."""
    from scipy.signal import find_peaks

    nu = arr["nu"]
    snr = _snr(arr)
    dnu, numax, eps = js["dnu"], js["numax"], js["epsilon"]
    up = D.universal_pattern(dnu)
    win = np.abs(nu - numax) < 3.5 * dnu
    nu_w, s_w = nu[win], S.smooth_boxcar(snr[win], 3)
    keep = P.dipole_mask(nu_w, dnu, eps, numax)
    tau, _ = P._stretch(nu_w, dnu, numax, eps, dpi1, q, up)
    df_ = nu_w[1] - nu_w[0]
    pk, props = find_peaks(np.where(keep, s_w, 0.0), height=thresh, distance=max(1, int(0.015 / df_)))
    return (tau[pk] % dpi1), nu_w[pk], props["peak_heights"]


def _stretched_fourier(arr, js, dpi1, q):
    """Fourier power of the stretched dipole spectrum at trial period spacings."""
    nu = arr["nu"]
    snr = _snr(arr)
    dnu, numax, eps = js["dnu"], js["numax"], js["epsilon"]
    up = D.universal_pattern(dnu)
    win = np.abs(nu - numax) < 3.5 * dnu
    nu_w = nu[win]
    keep = P.dipole_mask(nu_w, dnu, eps, numax)
    s_w = np.where(keep, np.clip(snr[win] - 1.0, -1.0, 50.0), 0.0)
    tau, dtau = P._stretch(nu_w, dnu, numax, eps, dpi1, q, up)
    trial = np.arange(30.0, 420.0, 0.25)
    pw = P.tau_spectrum(tau, s_w, dtau, trial)
    return trial, pw / np.mean(pw)


def fig_period_spacing(df, res):
    v = K.vrard2016()
    fig = plt.figure(figsize=(13.5, 7.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1, 1], hspace=0.5, wspace=0.3)
    ax = fig.add_subplot(gs[:, 0])
    ax.scatter(v["dnu"], v["dpi1"], s=3, color=ps.MUTED, alpha=0.35, lw=0,
               label=f"Vrard et al. 2016 ({len(v):,} giants)")
    det = df[df["dpi1_detected"]]
    for hyp, col, lab in (("RGB", C_RGB, "This work: shell H burning"),
                          ("helium burning", C_HEB, "This work: core He burning")):
        m = det["dpi1_hypothesis"] == hyp
        ax.errorbar(det["dnu"][m], det["dpi1"][m], yerr=det["dpi1_err"][m], fmt="o", ms=5, color=col,
                    mec=ps.SURFACE, mew=0.6, elinewidth=0.7, label=lab + f" ({int(m.sum())})")
    ax.set_xlabel("Large separation Δν (µHz)")
    ax.set_ylabel("Gravity-mode period spacing ΔΠ₁ (s)")
    ax.set_xlim(2.5, 18)
    ax.set_ylim(30, 420)
    ax.set_title("a. Period spacing against large separation")
    ax.legend(loc="center right", fontsize=7)
    for col_i, kic, kind in ((1, RGB_STAR, "Shell H burning"), (2, RC_STAR, "Core He burning")):
        js, arr = stars.load_star(kic)
        row = df.loc[kic]
        color = C_RGB if col_i == 1 else C_HEB
        # top: two orders of the spectrum around numax
        axs = fig.add_subplot(gs[0, col_i])
        nu, snr = arr["nu"], _snr(arr)
        dnu, eps = js["dnu"], js["epsilon"]
        n0 = np.floor(js["numax"] / dnu - eps)
        lo, hi = (n0 - 0.2 + eps) * dnu, (n0 + 1.85 + eps) * dnu
        m = (nu > lo) & (nu < hi)
        axs.plot(nu[m], S.smooth_boxcar(snr[m], 2), color=color, lw=0.6)
        up = D.universal_pattern(dnu)
        top = np.percentile(S.smooth_boxcar(snr[m], 2), 99.9) * 1.15
        for k in (0, 1):
            for name, x, ha in (("ℓ=0", 0.0, "left"), ("2", -up["d02"], "right"), ("ℓ=1", 0.5 - up["d01"], "center")):
                f = (n0 + k + eps + x) * dnu
                if lo < f < hi:
                    axs.text(f, top, name, ha=ha, fontsize=7, color=ps.INK_SECONDARY)
        axs.set_ylim(0, top * 1.12)
        axs.set_xlim(lo, hi)
        axs.set_xlabel("Frequency (µHz)")
        axs.set_ylabel("Signal / background")
        axs.set_title(f"{'bc'[col_i - 1]}. {kind}: KIC {kic}", fontsize=9.5)
        # bottom: Fourier spectrum of the stretched dipole spectrum
        axf = fig.add_subplot(gs[1, col_i])
        trial, pw = _stretched_fourier(arr, js, row["dpi1"], row["dpi1_q"])
        axf.plot(trial, pw, color=color, lw=0.9)
        axf.axvline(row["dpi1"], color=ps.INK_SECONDARY, lw=0.7)
        vv = row["dpi1_vrard"]
        label = f"ΔΠ₁ = {row['dpi1']:.1f} s" + (f"\n(Vrard et al.: {vv:.1f} s)" if np.isfinite(vv) else "")
        ps.annotate(axf, label, (row["dpi1"], pw[np.argmin(np.abs(trial - row['dpi1']))]),
                    (row["dpi1"] + 50 if col_i == 1 else 70, pw.max() * 0.8))
        axf.set_xlabel("Trial period spacing (s)")
        axf.set_ylabel("Fourier power (mean = 1)")
        axf.set_title(f"{'de'[col_i - 1]}. Its ℓ=1 modes, stretched\nand Fourier transformed", fontsize=9.5)
    fig.suptitle("Same surface, different cores: the period spacing of gravity modes tells hydrogen-shell "
                 "from helium-core burning", y=1.0, fontsize=11.5)
    ps.caption(fig, "a: ΔΠ₁ from our own stretched-period analysis of the dipole mixed modes (Mosser et al. 2015; "
                     "Vrard et al. 2016) on Vrard et al.'s measurements. b, c: two radial orders of each star. The "
                     "ℓ=1 modes are mixed with gravity waves in the core: dense in the shell-burning giant, sparser "
                     "around the helium-burning core. d, e: the dipole spectrum in 'stretched period', Fourier "
                     "transformed: a sharp peak at the core's gravity-mode period spacing. The broad bumps at long "
                     "periods in d come from the rise and fall of mode heights from one radial order to the next.")
    return ps.save(fig, FIG / "08_period_spacing")


# ---------------------------------------------------------------------------
# 9. HR diagrams
# ---------------------------------------------------------------------------


def fig_hr(df, res):
    import analysis as A

    bprp, gabs = A.b1_background()
    fig, axs = plt.subplots(1, 2, figsize=(12, 6.2))
    ax = axs[0]
    ax.hexbin(bprp, gabs, gridsize=220, cmap=GREY_RAMP, bins="log", mincnt=1, linewidths=0,
              extent=(-0.5, 4.5, -4, 17))
    for st, mk in (("RGB", "o"), ("RC", "s")):
        m = df["state_apokasc"] == st
        ax.plot(df["bp_rp0"][m], df["g_abs0"][m], mk, ms=4.5, color=_state_color(st), mec=ps.INK,
                mew=0.4, label=f"Kepler seismic giants: {STATE_LABEL[st].split(' (')[0].lower()}")
    ax.plot(0.82, 4.67, "*", ms=12, color=ps.YELLOW, mec=ps.INK, mew=0.6, label="The Sun")
    ax.set_xlim(-0.5, 4.0)
    ax.set_ylim(16, -3.5)
    ax.set_xlabel("BP − RP (dereddened for the giants)")
    ax.set_ylabel("Absolute G magnitude")
    ax.set_title("a. On B1's Gaia HR diagram (stars within 100 pc)")
    ax.legend(loc="lower left", fontsize=7)
    ax = axs[1]
    teff = df["teff_ak"]
    lum = df["L_self"]
    for st, mk in (("RGB", "o"), ("RC", "s")):
        m = (df["state_apokasc"] == st).to_numpy()
        sc = ax.scatter(teff[m], lum[m], c=df["M_self"][m], s=34, marker=mk, cmap=MASS_RAMP, vmin=0.8,
                        vmax=2.6, edgecolor=ps.INK_SECONDARY, linewidth=0.4,
                        label=STATE_LABEL[st].split(" (")[0])
    ax.legend(loc="upper left", fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(5300, 3900)
    ax.set_xticks([5200, 5000, 4800, 4600, 4400, 4200, 4000])
    ax.set_xticklabels(["5200", "5000", "4800", "4600", "4400", "4200", "4000"])
    ax.minorticks_off()
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Seismic mass (M☉)")
    ax.set_xlabel("Effective temperature (K, APOGEE)")
    ax.set_ylabel("Luminosity (L☉, from seismic radius)")
    ax.set_title("b. The seismic HR diagram, coloured by mass")
    ps.caption(fig, "a: B1's diagram (Gaia DR3, parallax > 10 mas, parallax/error > 10, RUWE < 1.4) with the "
                     "seismic giants at 0.4-4 kpc: parallaxes corrected for Gaia's zero point (Lindegren et al. "
                     "2021), reddening from Schlafly & Finkbeiner (2011) with Casagrande & VandenBerg (2018) "
                     "coefficients; the total dust column over-corrects the nearest few, which land too blue. "
                     "b: L = R² (Teff/5772 K)⁴ with the self-calibrated seismic radius.")
    return ps.save(fig, FIG / "09_hr_diagram")


# ---------------------------------------------------------------------------
# 10. The Sun
# ---------------------------------------------------------------------------


def fig_sun(df, res):
    sun = res["sun"]
    z = np.load(C.CACHE_DIR / "sun.npz")
    fig = plt.figure(figsize=(12.5, 8.4))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.32)
    ax = fig.add_subplot(gs[0, 0:2])
    nu_b, p_b = z["y2008_bin_nu"], z["y2008_bin_psd"]
    th = z["y2008_theta_post"]
    eta2 = S.apodization(nu_b, C.NYQUIST_VIRGO_UHZ)
    tr = L.filter_transfer(nu_b, C.TREND_SIGMA_DAYS_SC)
    comp = B.components(th, nu_b, eta2, tr)
    ax.loglog(nu_b, p_b, color=ps.MUTED, lw=0.6, label="VIRGO green, 2008")
    ax.loglog(nu_b, comp["gran1"] + comp["gran2"] + comp["slow"] + comp["white"], color=ps.INK_SECONDARY,
              lw=1.0, label="Background")
    ax.loglog(nu_b, comp["envelope"], color=ps.ORANGE, lw=1.2, label="Oscillation hump")
    ax.loglog(nu_b, sum(comp.values()), color=ps.BLUE, lw=1.3, label="Total model")
    ax.set_ylim(np.percentile(p_b, 0.5) * 0.5, p_b.max() * 2)
    ax.set_xlabel("Frequency (µHz)")
    ax.set_ylabel("Power density (ppm²/µHz)")
    ax.set_title("a. The Sun as a star, through the same pipeline (2008)")
    ax.legend(loc="lower left", fontsize=7)
    ax = fig.add_subplot(gs[0, 2])
    nu, snr = z["y2008_nu"], z["y2008_psd"] / z["y2008_bg"]
    y08 = sun["years"]["2008"]
    _echelle(ax, nu, snr, y08["dnu"], 2300, 4000, x_offset=0.0, smooth_frac=0.004)
    # label the ridges from the folded spectrum: the close pair is l=2 (left) and l=0 (right)
    from scipy.signal import find_peaks

    x, ce = D.collapsed_echelle(nu, S.smooth_boxcar(snr, 3), y08["dnu"], 2600, 3700, nbins=135)
    pk, pr = find_peaks(ce, distance=4)
    top3 = np.sort(pk[np.argsort(pr.get("peak_heights", ce[pk]) if "peak_heights" in pr else ce[pk])[-3:]])
    xs = x[top3] * y08["dnu"]
    gaps = np.diff(xs)
    j = int(np.argmin(gaps))  # the closest pair
    names = {j: "ℓ=2", j + 1: "ℓ=0"}
    for i, xx in enumerate(xs):
        ax.text(xx, 3990, names.get(i, "ℓ=1"), ha="center", va="top", fontsize=8, color=ps.INK,
                bbox={"boxstyle": "round,pad=0.15", "fc": ps.SURFACE, "ec": "none", "alpha": 0.85})
    ax.set_xlabel(f"Frequency mod Δν = {y08['dnu']:.2f} µHz (µHz)")
    ax.set_ylabel("Frequency (µHz)")
    ax.set_title("b. The Sun's échelle diagram", fontsize=10)
    years = sorted(int(y) for y in sun["years"])
    good = set(sun["years_used"])
    yr = np.array(years)
    nm = np.array([sun["years"][str(y)]["numax"] for y in years])
    nme = np.array([sun["years"][str(y)]["numax_err"] for y in years])
    dn = np.array([sun["years"][str(y)]["dnu"] for y in years])
    dne = np.array([sun["years"][str(y)]["dnu_err"] for y in years])
    used = np.array([y in good for y in years])
    ssn = pd.read_csv(C.DATA_DIR / "silso" / "SN_y_tot_V2.0.csv", sep=";", header=None)
    ssn = ssn[(ssn[0] > 1995) & (ssn[0] < 2024)]
    for k, (vals, errs, lab, ref) in enumerate((
            (nm, nme, "νmax☉ (µHz)", sun["numax_sun"]), (dn, dne, "Δν☉ (µHz)", sun["dnu_sun"]),
            (None, None, "Sunspot number", None))):
        ax = fig.add_subplot(gs[1, k])
        if vals is None:
            ax.bar(ssn[0] - 0.5, ssn[1], width=0.8, color=ps.MUTED)
            ax.set_title("e. The activity cycle (WDC-SILSO)", fontsize=9)
        else:
            ax.errorbar(yr[used], vals[used], yerr=errs[used], fmt="o", ms=4, color=ps.BLUE, elinewidth=0.8)
            ax.errorbar(yr[~used], vals[~used], yerr=errs[~used], fmt="o", ms=4, color=ps.BLUE, mfc=ps.SURFACE,
                        elinewidth=0.8)
            ax.axhline(ref, color=ps.INK_SECONDARY, lw=0.8)
            ax.set_title(("c. " if k == 0 else "d. ") + lab.split(" (")[0] + f" each year (mean {ref:.1f})"
                         if k == 0 else "d. " + lab.split(" (")[0] + f" each year (mean {ref:.2f})", fontsize=9.5)
        ax.set_xlabel("Year")
        ax.set_ylabel(lab)
        ax.set_xlim(1995, 2023.5)
    cyc = res["sun_cycle"]["all years"]
    ps.caption(fig, "SOHO/VIRGO sun-photometer data (green channel, 1-minute cadence), one calendar year at a "
                     "time, analysed like the Kepler giants (the narrow spike near 5.5 mHz is instrumental). Open "
                     "circle: 1998, 67% coverage, not used. νmax☉ rises with solar activity, "
                     f"{cyc['shift_0_to_170']:.0f} ± {cyc['shift_err']:.0f} µHz from no sunspots to 170 "
                     "(Howe et al. 2020 found ~25 µHz with BiSON); the low values after 2015 may be partly "
                     "instrumental. These numbers are the reference point of every scaling relation.")
    return ps.save(fig, FIG / "10_sun")


# ---------------------------------------------------------------------------
# 11. 16 Cyg A and B
# ---------------------------------------------------------------------------


def fig_16cyg(df, res):
    ms = json.loads((C.CACHE_DIR / "sixteen_cyg.json").read_text())
    z = np.load(C.CACHE_DIR / "sixteen_cyg.npz")
    _, t6 = K.lund2017()
    fig = plt.figure(figsize=(13, 6.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.1], wspace=0.34, hspace=0.6)
    for i, (name, kic) in enumerate((("16 Cyg A", "12069424"), ("16 Cyg B", "12069449"))):
        ax = fig.add_subplot(gs[:, i])
        key = name.replace(" ", "_").lower()
        r = ms[name]
        _echelle_panel(ax, z[f"{key}_nu"], z[f"{key}_psd"] / z[f"{key}_bg"], r["dnu"], r["numax"],
                       r["envelope_fwhm"] / 2.3548, r["epsilon"], f"{'ab'[i]}. {name}", labels=False,
                       lund=t6[t6["kic"] == kic], n_sig=1.6, smooth_frac=0.01, legend=(i == 0))
    sc = res["sixteen_cyg"]
    for k, (q, unit, lit_keys) in enumerate((("radius", "R☉", ("radius_interf", "radius_model")),
                                             ("mass", "M☉", (None, "mass_model")))):
        ax = fig.add_subplot(gs[k, 2])
        for j, name in enumerate(("16 Cyg A", "16 Cyg B")):
            s_ = sc[name]
            lit = s_["literature"]
            y = 1 - j
            v, e = s_["scaling"]["self"][q], s_["scaling"]["self"][f"{q}_err"]
            ax.errorbar(v, y + 0.15, xerr=e, fmt="o", color=ps.BLUE, ms=6,
                        label="This work: scaling relations, our Sun" if j == 0 else None)
            if lit_keys[0]:
                ax.errorbar(lit[lit_keys[0]], y, xerr=lit[f"{lit_keys[0]}_err"], fmt="s", color=ps.ORANGE,
                            ms=6, label="Interferometry (White et al. 2013)" if j == 0 else None)
            ax.errorbar(lit[lit_keys[1]], y - 0.15, xerr=lit[f"{lit_keys[1]}_err"], fmt="D", color=ps.AQUA,
                        ms=5.5, label="Models (Metcalfe et al. 2015)" if j == 0 else None)
        ax.set_yticks([1, 0])
        ax.set_yticklabels(["16 Cyg A", "16 Cyg B"])
        ax.set_ylim(-0.5, 1.5)
        ax.grid(axis="y", visible=False)
        ax.set_xlabel(f"{q.capitalize()} ({unit})")
        ax.set_title(f"{'cd'[k]}. {q.capitalize()}", fontsize=10)
        if k == 0:
            ax.legend(loc="center", fontsize=7, ncol=1)
    ps.caption(fig, "One year (Quarters 7-10) of 1-minute Kepler photometry from custom pixel masks. Orange "
                     "symbols: Lund et al. (2017)'s peak-bagged mode frequencies, which sit on our ridges. Radii "
                     "and masses from νmax and Δν alone (scaling relations calibrated on our own Sun, Teff from "
                     "Lund et al. 2017), against CHARA interferometric radii and detailed asteroseismic models.",
               y=-0.07)
    return ps.save(fig, FIG / "11_16cyg")


ALL = [fig_light_to_spectrum, fig_ladder, fig_echelle, fig_vs_catalogues, fig_dnu_numax, fig_masses_radii,
       fig_gaia, fig_period_spacing, fig_hr, fig_sun, fig_16cyg]


def make_all(*, quiet: bool = False):
    df, res = _results()
    written = []
    for f in ALL:
        out = f(df, res)
        written += out
        if not quiet:
            print(f"  {out[0]}")
    return written


if __name__ == "__main__":
    make_all()
