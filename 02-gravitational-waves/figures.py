"""Figure generation.

Every figure is drawn from the compact analysis artefacts, never by re-running
the analysis, so what is plotted is exactly what was measured.

Design notes. A colour means one thing everywhere: blue is Hanford, orange is
Livingston, black is the general-relativity waveform, violet a published value.
The palette is the colourblind-safe set used across these projects (validated
with a script, not by eye). Time runs in seconds after 2015-09-14 09:50:45 UTC,
the discovery paper's convention, so any panel can be set beside the paper's.
Each figure carries a caption with its data provenance and the caveat a reader
needs in order not to over-read it.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib.ticker import NullFormatter

import gw_physics as gp
import plotstyle as ps
from config import DISPLAY_BAND_HZ, FIGURE_DIR

ps.use_style()

PROVENANCE = (
    "Data: LIGO Hanford and Livingston strain, GWOSC (gwosc.org), GWTC-1 release."
)
T_LABEL = "Time (s) after 2015-09-14 09:50:45 UTC"
MSUN = r"$M_\odot$"
MC = r"$\mathcal{M}$"
WINDOW = (0.25, 0.46)  # the stretch of the discovery paper's Figure 1


def _clip(t, *arrays, lo=WINDOW[0], hi=WINDOW[1]):
    keep = (t >= lo) & (t <= hi)
    return (t[keep],) + tuple(a[keep] for a in arrays)


def _log_smooth(f, amp, per_decade=40):
    """RMS of ``amp`` in logarithmic frequency bins: the envelope, without the ripples."""
    edges = np.logspace(np.log10(f[f > 0].min()), np.log10(f.max()), int(per_decade * np.log10(f.max() / f[f > 0].min())))
    idx = np.digitize(f, edges)
    fc, ac = [], []
    for i in np.unique(idx):
        m = idx == i
        if m.sum():
            fc.append(np.exp(np.mean(np.log(f[m]))))
            ac.append(np.sqrt(np.mean(amp[m] ** 2)))
    return np.array(fc), np.array(ac)


def _label(ax, x, y, text, *, color=ps.INK_SECONDARY, ha="left", va="center", **kw):
    """Direct label in a text token (never the series colour)."""
    return ax.text(x, y, text, color=color, fontsize=8, ha=ha, va=va, **kw)


# --------------------------------------------------------------------------
# Reusable panels (each is used by a single-topic figure and by the summary)
# --------------------------------------------------------------------------

def panel_overlay(ax, res, art):
    """The discovery paper's Figure 1 view: L1, and H1 shifted and flipped onto it."""
    lag = art["arrival"]["lag_mf_ms"] / 1e3
    t = art["det"]["L1"]["prl"]["t"]
    l1 = art["det"]["L1"]["prl"]["data"]
    h1 = art["det"]["H1"]["prl"]["data"]
    h1_aligned = -np.interp(t + lag, t, h1)
    tt, l1c, h1c = _clip(t, l1, h1_aligned)
    ax.plot(tt, h1c * 1e21, color=ps.C_H1, lw=1.6,
            label=f"Hanford, shifted {lag * 1e3:.1f} ms earlier and inverted")
    ax.plot(tt, l1c * 1e21, color=ps.C_L1, lw=1.6, label="Livingston")
    ax.axhline(0, color=ps.AXIS, lw=0.8, zorder=0)
    ax.set_xlim(*WINDOW)
    ax.set_ylim(-1.45, 1.45)
    ax.set_xlabel(T_LABEL)
    ax.set_ylabel(r"Strain ($\times10^{-21}$)")
    ax.legend(loc="upper left", ncol=2)


def panel_qplane(ax, res, art, name, *, vmax=25.0, track=True):
    """Q-transform of one detector, with the fitted Newtonian chirp traced over it."""
    q = art["det"][name]["qplane"]
    f = q["frequencies"]
    mesh = ax.pcolormesh(q["times"], f, q["energy"], shading="nearest",
                         cmap=ps.sequential_cmap(), vmin=0.0, vmax=vmax, rasterized=True)
    ax.set_yscale("log")
    ax.set_ylim(20, 500)
    ax.set_xlim(0.20, 0.50)
    ax.set_yticks([20, 30, 50, 100, 200, 300, 500])
    ax.set_yticklabels(["20", "30", "50", "100", "200", "300", "500"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(False)
    if track:
        ch = art["chirp"]["data"]
        tc = ch["t_coalescence"] + (art["arrival"]["lag_mf_ms"] / 1e3 if name == "H1" else 0.0)
        tt = np.linspace(0.20, tc - 1e-3, 400)
        ff = gp.frequency_at(tc - tt, ch["chirp_mass"])
        keep = (ff > 20) & (ff < 350)
        ax.plot(tt[keep], ff[keep], color=ps.INK, lw=1.0, alpha=0.9)
    return mesh


def panel_chirp_fit(ax, res, art, *, compact=False):
    """f^(-8/3) against time: the chirp law as a straight line (the LIGO authors' Fig. 3)."""
    ch = art["chirp"]["data"]
    tm, fz, kept = ch["t_mid"], ch["freq"], ch["kept"]
    v = fz ** (-8.0 / 3.0) * 1e5
    others = np.setdiff1d(np.arange(len(fz)), kept)
    near = others[(tm[others] > WINDOW[0] + 0.05) & (tm[others] < ch["t_peak"] + 0.01)
                  & (v[others] < 12)]
    ax.plot(tm[near], v[near], "o", color=ps.C_CONTEXT, ms=4, alpha=0.6, zorder=2,
            label="half-cycles below the 2.5σ gate (not fitted)")
    ax.plot(tm[kept], v[kept], "o", color=ps.INK, ms=6.5, mec=ps.SURFACE, mew=1.5, zorder=4,
            label="half-cycles used in the fit")

    tc = ch["t_coalescence"]
    t0 = tm[kept].min() - 0.01
    tt = np.linspace(t0, tc, 50)
    ax.plot(tt, (ch["slope"] * tt + ch["intercept"]) * 1e5, color=ps.C_ACCENT, lw=2.0, zorder=3)
    xa = t0 + 0.35 * (tc - t0)
    ps.annotate(ax, f"straight-line fit:\n{MC} = {ch['chirp_mass']:.1f} {MSUN} (Newtonian)",
                xy=(xa, (ch["slope"] * xa + ch["intercept"]) * 1e5), xytext=(t0 - 0.02, 1.2),
                color=ps.INK, ha="left", va="center")
    for mc_ref in (30.0, 40.0):
        s = gp.chirp_slope(mc_ref)
        ax.plot(tt, s * (tt - tc) * 1e5, color=ps.C_CONTEXT, lw=0.9, zorder=1)
        _label(ax, t0, s * (t0 - tc) * 1e5, f"{MC} = {mc_ref:.0f} {MSUN}  ", ha="right",
               color=ps.MUTED)
    ax.axhline(0, color=ps.AXIS, lw=0.8, zorder=0)
    ax.set_xlim(t0 - (0.035 if not compact else 0.03), tc + 0.01)
    ax.set_ylim(-0.3, max(v[kept].max() * 1.25, 9.0))
    ax.set_xlabel(T_LABEL)
    ax.set_ylabel(r"$f^{-8/3}$  ($\times10^{-5}$ Hz$^{-8/3}$)")
    if not compact:
        ax.legend(loc="upper right")


def panel_background(ax, res, art):
    """Search results against the time-slide background, as in the discovery paper's Fig. 4."""
    sig = art["significance"]
    scale = sig["analysed_seconds"] / sig["background_livetime_s"]
    grid = np.linspace(4.0, 26.0, 900)

    def expected(stats):
        stats = np.sort(stats)
        return (len(stats) - np.searchsorted(stats, grid, side="left")) * scale

    for key, colour, text in (
        ("all_triggers", ps.C_CONTEXT, "background (all triggers)"),
        ("event_removed", ps.C_ACCENT, "background (event removed)"),
    ):
        y = expected(sig["background"][key])
        y = np.where(y > 0, y, np.nan)
        ax.step(grid, y, where="post", color=colour, lw=2.0, label=text)

    zl = np.sort(sig["zero_lag"]["statistic"])
    counts = len(zl) - np.arange(len(zl))
    ax.plot(zl, counts, "s", color=ps.INK, ms=5, mec=ps.SURFACE, mew=1.0, zorder=5,
            label="the real hour")
    ev = sig["event_statistic"]
    ax.annotate(f"GW150914, " + r"$\hat\rho_c$" + f" = {ev:.1f}", xy=(ev, 1), xytext=(ev - 1.5, 4.0),
                fontsize=8, color=ps.INK, ha="right", va="center",
                arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "lw": 0.8, "shrinkB": 4})
    ax.set_yscale("log")
    ax.set_xlim(4.0, 26.0)
    ax.set_ylim(1.0 / sig["slides"] / 3, 400)
    ax.set_xlabel(r"Ranking statistic $\hat\rho_c$ (network SNR, $\chi^2$-reweighted)")
    ax.set_ylabel("Events at least this loud,\nper 67 min of data")
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.43))


# --------------------------------------------------------------------------
# 1. The noise
# --------------------------------------------------------------------------

def figure_noise(res, art):
    """Amplitude spectral density of both detectors, with the signal's spectrum on top.

    The signal is drawn as 2 sqrt(f) |h(f)|, the convention in which its height
    above a noise curve shows where its SNR comes from: the squared ratio,
    integrated over log-frequency, is SNR^2. That integral is recomputed here
    as a consistency check against the matched filter.
    """
    fig, ax = plt.subplots(figsize=(10, 5.6))
    for name in ("H1", "L1"):
        n = art["det"][name]["noise"]
        f, asd = n["freqs"], np.sqrt(n["psd"])
        sel = (f >= 10) & (f <= 2000)
        ax.plot(f[sel], asd[sel], color=ps.DETECTOR_COLOUR[name], lw=1.0,
                label=ps.DETECTOR_NAME[name])

    s = art["det"]["H1"]["signal_spectrum"]
    f, amp = s["f"], 2.0 * np.sqrt(s["f"]) * s["abs_h"]
    sel = (f >= 20) & (f <= 400)
    fs_, as_ = _log_smooth(f[sel], amp[sel])
    ax.plot(fs_, as_, color=ps.C_MODEL, lw=2.0,
            label=r"GW150914 in Hanford (best-fit waveform), $2\sqrt{f}\,|\tilde h(f)|$")

    # Consistency: SNR^2 = integral (2 sqrt f |h|)^2 / S d(ln f) with the same PSD.
    n = art["det"]["H1"]["noise"]
    band = (f >= 20) & (f <= 2000)
    psd_i = np.interp(f[band], n["freqs"], n["psd"])
    df = f[1] - f[0]
    snr_check = np.sqrt(np.sum(4 * s["abs_h"][band] ** 2 / psd_i) * df)
    snr_mf = res["detectors"]["H1"]["matched_filter"]["snr"]

    ax.axvspan(*DISPLAY_BAND_HZ, color=ps.GRID, alpha=0.45, zorder=0, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, 2000)
    ax.set_ylim(3e-24, 6e-19)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(r"Strain noise (1/$\sqrt{\mathrm{Hz}}$)")
    ax.set_title("The noise: what each detector hears when nothing is happening")
    ax.legend(loc="upper right")

    _label(ax, 36, 4.2e-19, "35–350 Hz: the band\nthe discovery paper plots", va="top")
    ps.annotate(ax, "seismic and control noise:\nthe ground shaking the mirrors", xy=(12.5, 6e-21),
                xytext=(10.4, 6.5e-24), ha="left", va="bottom")
    ps.annotate(ax, "60 Hz mains\nand harmonics", xy=(60, 8e-22), xytext=(75, 2.5e-21))
    for x, y in ((36.7, 2e-21), (331.9, 1.5e-21)):
        ps.annotate(ax, "calibration lines,\ninjected on purpose" if x > 100 else "",
                    xy=(x, y), xytext=(190, 1.3e-20), ha="center")
    ps.annotate(ax, "violin resonances of the\nfibres the mirrors hang from", xy=(505, 5e-22),
                xytext=(420, 6e-21), ha="left")
    ps.annotate(ax, "quantum shot noise:\nthe laser light's graininess", xy=(1300, 2.3e-23),
                xytext=(1700, 6.5e-24), ha="right", va="bottom")

    ps.caption(fig, f"{PROVENANCE} Noise: median-averaged Welch spectrum over the hour around the "
               f"event (4 s segments). Where the black curve stands above a noise curve, the signal "
               f"out-shouts the noise at that frequency; the squared ratio integrated over log-frequency "
               f"is SNR² — recomputed from this plot it gives {snr_check:.1f}, against {snr_mf:.1f} from the "
               f"matched filter. Below ~20 Hz the noise climbs steeply; in the raw time series that "
               f"low-frequency noise dominates everything, which is why the raw strain shows nothing.")
    return ps.save(fig, FIGURE_DIR / "01_noise_spectrum"), {"snr_from_spectrum": float(snr_check)}


# --------------------------------------------------------------------------
# 2. From raw strain to a chirp
# --------------------------------------------------------------------------

def figure_reveal(res, art):
    """Raw strain, whitened strain, and the zoom where the chirp is plainly there."""
    fig, axes = plt.subplots(3, 2, figsize=(11, 8.6),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.2})
    for col, name in enumerate(("H1", "L1")):
        d = art["det"][name]
        c = ps.DETECTOR_COLOUR[name]

        ax = axes[0, col]
        ax.plot(d["raw"]["t"], d["raw"]["h"] * 1e19, color=c, lw=0.7, rasterized=True)
        ax.set_title(f"{ps.DETECTOR_NAME[name]}\nRaw strain, 32 s")
        ax.tick_params(labelbottom=False)
        ax.set_ylabel(r"Strain ($\times10^{-19}$)")
        ax.set_xlim(d["raw"]["t"][0], d["raw"]["t"][-1])

        ax = axes[1, col]
        ax.plot(d["white16"]["t"], d["white16"]["x"], color=c, lw=0.3, rasterized=True)
        ax.set_title("Whitened, 35–350 Hz, same 32 s")
        ax.set_ylabel("Noise σ")
        ax.set_xlim(d["white16"]["t"][0], d["white16"]["t"][-1])
        ax.set_ylim(-10, 10)
        cen = d["census"]
        louder = cen["louder_moments_in_hour"]
        ps.annotate(ax, f"GW150914: {cen['event_peak_sigma']:.1f}σ.\nLouder moments in the\n"
                        f"whole hour: {louder}", xy=(0.42, cen["event_peak_sigma"] * 0.95),
                    xytext=(4.0, 9.6), ha="left", va="top")

        ax = axes[2, col]
        tt, xx = _clip(d["white3"]["t"], d["white3"]["x"])
        ax.plot(tt, xx, color=c, lw=1.4)
        ax.axhline(0, color=ps.AXIS, lw=0.8, zorder=0)
        ax.set_title("Zoom on the last 0.2 s")
        ax.set_ylabel("Noise σ")
        ax.set_xlim(*WINDOW)
        ax.set_ylim(-10, 10)
        ax.set_xlabel(T_LABEL)
    for ax in axes[1]:
        ax.set_xlabel(T_LABEL)

    h1 = art["det"]["H1"]
    ratio = np.std(h1["raw"]["h"]) / res["detectors"]["H1"]["peak_strain_filtered"]
    ps.caption(fig, f"{PROVENANCE} Top: the calibrated strain as recorded; its scatter, dominated by noise "
               f"below 20 Hz, is ~{ratio:.0f} times the signal's peak. Middle: each frequency divided by its own noise level "
               f"(whitening), then restricted to 35–350 Hz; the y axis is in units of the noise's standard "
               f"deviation, measured over the rest of the hour. Bottom: the chirp — amplitude and frequency "
               f"rising together — in both detectors, 3000 km apart, before any template is used.")
    return ps.save(fig, FIGURE_DIR / "02_raw_to_chirp")


# --------------------------------------------------------------------------
# 3. Two detectors, and the check against the paper
# --------------------------------------------------------------------------

def figure_two_detectors(res, art):
    """PRL Figure 1 reproduced from raw data; our curves vs theirs; the arrival lag."""
    from fetch_data import load_prl_figure1

    paper = load_prl_figure1()
    fig = plt.figure(figsize=(11, 8.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1.0], hspace=0.55, wspace=0.32)

    ax0 = fig.add_subplot(gs[0, :])
    panel_overlay(ax0, res, art)
    ax0.set_title("One signal in two detectors 3,000 km apart — reproduced from the raw strain")

    for i, name in enumerate(("H1", "L1")):
        ax = fig.add_subplot(gs[1, i])
        d = art["det"][name]
        tt, ours = _clip(d["prl"]["t"], d["prl"]["data"])
        pp = paper[f"observed-{name[0]}"]
        tp, hp = _clip(pp["t"], pp["strain"])
        ax.plot(tt, ours * 1e21, color=ps.DETECTOR_COLOUR[name], lw=2.0, label="ours, from raw strain")
        ax.plot(tp, hp * 1e21, color=ps.C_PAPER, lw=0.9, label="discovery paper, Fig. 1")
        m = res["detectors"][name]["paper_figure1"]["observed"]
        ax.set_title(f"{name}: our processing vs the paper's")
        ax.text(0.03, 0.05, f"correlation {m['correlation']:.3f}\n"
                            f"amplitude ratio {m['rms_ratio_ours_over_paper']:.3f}",
                transform=ax.transAxes, fontsize=8, color=ps.INK_SECONDARY)
        ax.set_xlim(0.33, 0.44)
        ax.set_ylim(-1.45, 1.45)
        ax.set_xlabel("Time (s) after 09:50:45 UTC")
        ax.set_ylabel(r"Strain ($\times10^{-21}$)")
        ax.axhline(0, color=ps.AXIS, lw=0.8, zorder=0)
        if i == 0:
            ax.legend(loc="upper left")

    ax = fig.add_subplot(gs[1, 2])
    arr, pub = art["arrival"], art["published"]["predicted_lag_ms"]
    rows = [
        ("matched filter (ours)", arr["lag_mf_ms"], arr["lag_mf_err_ms"], arr["lag_mf_err_ms"], ps.INK),
        ("cross-correlation (ours)", arr["lag_xcorr_ms"], arr["lag_xcorr_err_ms"],
         arr["lag_xcorr_err_ms"], ps.INK),
        ("discovery paper", 6.9, 0.4, 0.5, ps.C_PAPER),
        ("GWTC-1 sky positions", pub["median"], pub["median"] - pub["p5"], pub["p95"] - pub["median"],
         ps.C_PAPER),
    ]
    for j, (label, v, lo, hi, colour) in enumerate(rows):
        y = len(rows) - 1 - j
        ax.errorbar(v, y, xerr=[[lo], [hi]], fmt="o", color=colour, ms=6, mec=ps.SURFACE, mew=1.5,
                    elinewidth=1.6, capsize=0)
        _label(ax, 4.15, y + 0.28, label, va="bottom")
    limit = res["geometry"]["light_travel_ms"]
    ax.axvspan(limit, 12, color=ps.GRID, alpha=0.8, lw=0)
    ax.axvline(limit, color=ps.INK_SECONDARY, lw=1.0)
    _label(ax, limit - 0.1, -0.75, f"light-travel\nlimit {limit:.1f} ms", ha="right", va="bottom")
    ax.set_xlim(4, 11)
    ax.set_ylim(-0.9, len(rows) - 0.2)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Hanford arrival after Livingston (ms)")
    ax.set_title("Arrival-time difference")

    ps.caption(fig, f"{PROVENANCE} Top: 35–350 Hz band-pass plus notches at the instrumental lines, the "
               f"paper's recipe; Hanford is advanced by our measured lag and flipped in sign, because the "
               f"two L-shaped detectors point their arms nearly at right angles to each other (their response "
               f"tensors overlap at {res['geometry']['tensor_overlap']:.2f}; exactly perpendicular would be "
               f"−1). Bottom: our filtered curves "
               f"overlaid on the collaboration's published ones. Error bars on our lags are the 68% spread "
               f"over {res['arrival']['injections']} injections into real noise; published intervals are 90%.")
    return ps.save(fig, FIGURE_DIR / "03_two_detectors")


# --------------------------------------------------------------------------
# 4. Time-frequency
# --------------------------------------------------------------------------

def figure_spectrogram(res, art):
    """Q-transforms of both detectors, the fitted Newtonian chirp traced over them."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    for ax, name in zip(axes, ("H1", "L1")):
        mesh = panel_qplane(ax, res, art, name)
        q = art["det"][name]["qplane"]
        ax.set_title(f"{ps.DETECTOR_NAME[name]}   (Q = {q['q']:.1f}, "
                     f"peak energy {q['peak_energy']:.0f})")
        ax.set_xlabel(T_LABEL)
    axes[0].set_ylabel("Frequency (Hz)")
    ch = art["chirp"]["data"]
    _label(axes[1], 0.205, 300, f"black line: Newtonian chirp with {MC} = {ch['chirp_mass']:.1f} {MSUN},\n"
                                 "fitted to zero crossings (Fig. 6), not to this plot",
           va="top", color=ps.INK)
    cbar = fig.colorbar(mesh, ax=axes, fraction=0.025, pad=0.015)
    cbar.set_label("Normalised energy")
    cbar.outline.set_visible(False)
    ps.caption(fig, f"{PROVENANCE} Q-transform of whitened strain: each row is a band-pass whose width grows "
               f"with frequency, so every tile holds the same number of wave cycles. Normalised energy is "
               f"about 1 for noise; a tile of energy Z has chance probability e^-Z. The colour scale "
               f"saturates at 25. The track climbs from ~35 Hz to ~250 Hz in 0.2 s — the chirp.")
    return ps.save(fig, FIGURE_DIR / "04_spectrogram")


# --------------------------------------------------------------------------
# 5. Matched filter
# --------------------------------------------------------------------------

def figure_matched_filter(res, art):
    """An hour of SNR in each detector, then the data against the best-fit waveform."""
    fig = plt.figure(figsize=(11, 9.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[0.75, 0.75, 1.0, 1.0], hspace=0.75)
    ev_rel = art["event_gps"] - art["ref"]

    for i, name in enumerate(("H1", "L1")):
        ax = fig.add_subplot(gs[i])
        d = art["det"][name]
        t_min = (d["snr_thin"]["t"] - ev_rel) / 60.0
        edge = d["mf"]["edge_s"] / 60.0
        keep = (t_min > t_min[0] + edge) & (t_min < t_min[-1] - edge)
        ax.plot(t_min[keep], d["snr_thin"]["snr"][keep], color=ps.DETECTOR_COLOUR[name], lw=0.6)
        ax.set_ylim(0, 22)
        ax.set_xlim(t_min[keep][0], t_min[keep][-1])
        ax.set_ylabel("SNR")
        ax.set_title(f"{ps.DETECTOR_NAME[name]}: matched-filter SNR across the hour")
        ps.annotate(ax, f"GW150914, SNR {d['mf']['snr']:.1f}\n"
                        f"χ²/dof {d['mf']['chisq'] / d['mf']['chisq_dof']:.2f} (signal-like)",
                    xy=(0.0, d["mf"]["snr"]), xytext=(3.0, 17.5), ha="left", va="top")
        tr = d["triggers"]
        away = np.abs(tr["t"] - ev_rel) > 8
        j = int(np.argmax(np.where(away, tr["snr"], 0)))
        if tr["snr"][j] > 10:
            tj = (tr["t"][j] - ev_rel) / 60.0
            ps.annotate(ax, f"glitch, SNR {tr['snr'][j]:.1f} — but χ²/dof {tr['chisq'][j] / d['mf']['chisq_dof']:.0f},\n"
                            f"so its reweighted SNR ρ̂ is just {tr['rho_hat'][j]:.1f}",
                        xy=(tj, tr["snr"][j]), xytext=(tj - 20, 19.5), ha="left", va="top")
        if i == 1:
            ax.set_xlabel("Minutes from GW150914")

    for i, name in enumerate(("H1", "L1")):
        ax = fig.add_subplot(gs[2 + i])
        d = art["det"][name]
        tt, xx, ff = _clip(d["white3"]["t"], d["white3"]["x"], d["white3"]["fit"])
        ax.plot(tt, xx, color=ps.DETECTOR_COLOUR[name], lw=1.2, label="whitened data")
        ax.plot(tt, ff, color=ps.C_MODEL, lw=1.6, label="best-fit general-relativity waveform")
        core = (tt > 0.30) & (tt < 0.44)
        resid = float(np.std(xx[core] - ff[core]))
        ax.axhline(0, color=ps.AXIS, lw=0.8, zorder=0)
        ax.set_xlim(*WINDOW)
        ax.set_ylim(-10, 10)
        ax.set_ylabel("Noise σ")
        ax.set_title(f"{name}: data and waveform (residual scatter {resid:.2f}σ; pure noise gives 1)")
        if i == 0:
            ax.legend(loc="upper left", ncol=2)
        else:
            ax.set_xlabel(T_LABEL)

    tm = art["template_chirp_mass"]
    ps.caption(fig, f"{PROVENANCE} Template: GWOSC's SEOBNRv2 waveform for GW150914 (detector-frame "
               f"masses 41.7 and 29.2 {MSUN}, {MC} = {tm:.1f} {MSUN}). Only its amplitude, phase and "
               f"arrival time are fitted to each detector. SNR is maximised over phase; χ²/dof tests "
               f"whether the power is spread across frequency the way the waveform predicts — near 1 for a "
               f"real signal, large for a glitch. Network SNR {res['network_snr']:.1f} "
               f"(the published searches: 23.6–25.2).")
    return ps.save(fig, FIGURE_DIR / "05_matched_filter")


# --------------------------------------------------------------------------
# 6. Chirp mass
# --------------------------------------------------------------------------

def figure_chirp_mass(res, art):
    """Zero crossings -> f^(-8/3) line -> chirp mass, calibrated by injection."""
    ch = art["chirp"]
    d = ch["data"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.6), gridspec_kw={"hspace": 0.52, "wspace": 0.28})

    # (a) the combined stream and its gated half-cycles
    ax = axes[0, 0]
    tt, yy = _clip(ch["time"], ch["combined"], lo=0.28, hi=0.45)
    ax.plot(tt, yy * 1e21, color=ps.INK, lw=1.3, label="combined strain")
    g = res["chirp"]["gate_sigma"] * ch["noise_sigma"] * 1e21
    for k, s in enumerate((-g, g)):
        ax.axhline(s, color=ps.MUTED, lw=0.9, ls=(0, (4, 3)),
                   label=f"gate: ±{res['chirp']['gate_sigma']:.1f} noise σ" if k else None)
    zc = d["crossings"]
    used = np.unique(np.concatenate([zc[d["kept"]], zc[d["kept"] + 1]]))
    ax.plot(used, np.zeros_like(used), "o", color=ps.C_ACCENT, ms=6, mec=ps.SURFACE, mew=1.2, zorder=4,
            label="zero crossings used")
    ax.axhline(0, color=ps.AXIS, lw=0.8, zorder=0)
    ax.set_ylim(-3.5, 2.4)
    ax.legend(loc="lower right", ncol=3)
    ax.set_xlim(0.28, 0.45)
    ax.set_title("1  Livingston plus shifted, inverted Hanford")
    ax.set_ylabel(r"Strain ($\times10^{-21}$)")
    ax.set_xlabel(T_LABEL)

    # (b) the straight line
    panel_chirp_fit(axes[0, 1], res, art)
    axes[0, 1].set_title(r"2  Each half-cycle's frequency, as $f^{-8/3}$")

    # (c) injections
    ax = axes[1, 0]
    rec = ch["injections"]["recovered"]
    rec = rec[np.isfinite(rec)]
    bins = np.arange(10, 56, 1.5)
    ax.hist(rec, bins=bins, color=ps.C_CONTEXT, alpha=0.55, rwidth=0.85, lw=0)
    top = ax.get_ylim()[1]
    ax.axvline(ch["template_chirp_mass"], color=ps.C_PAPER, lw=2.0)
    _label(ax, ch["template_chirp_mass"] - 0.8, top * 0.97,
           f"true {MC} of the\ninjected waveform:\n{ch['template_chirp_mass']:.1f}", ha="right", va="top")
    ax.axvline(d["chirp_mass"], color=ps.INK, lw=2.0)
    ps.annotate(ax, f"GW150914 itself: {d['chirp_mass']:.1f}\n— a typical draw", xy=(d["chirp_mass"], top * 0.8),
                xytext=(42.5, top * 0.8), ha="left", va="center", color=ps.INK)
    ax.set_xlim(10, 55)
    ax.set_xlabel(f"Chirp mass the method returns ({MSUN})")
    ax.set_ylabel("Injections")
    n_inj = len(ch["injections"]["recovered"])
    ax.set_title(f"3  The same method on {n_inj} injections")

    # (d) the comparison
    ax = axes[1, 1]
    inj = ch["injections"]
    pub = art["published"]["chirp_mass_detector_frame"]
    raw = d["chirp_mass"]
    rows = [
        ("Newtonian fit to zero crossings (ours)", raw, raw * (1 - inj["p16_p84"][0] / inj["median"]),
         raw * (inj["p16_p84"][1] / inj["median"] - 1), ps.INK),
        ("same method, LIGO authors (Ann. Phys. 2017)", 37.0, None, None, ps.C_PAPER),
        ("ours, bias-calibrated by injections", ch["calibrated"], ch["calibrated"] - ch["calibrated_68"][0],
         ch["calibrated_68"][1] - ch["calibrated"], ps.INK),
        ("GWTC-1 posterior, detector frame (90%)", pub["median"], pub["median"] - pub["p5"],
         pub["p95"] - pub["median"], ps.C_PAPER),
    ]
    for j, (label, v, lo, hi, colour) in enumerate(rows):
        y = len(rows) - 1 - j
        if lo is None:
            ax.plot(v, y, "o", color=colour, ms=6, mec=ps.SURFACE, mew=1.5)
        else:
            ax.errorbar(v, y, xerr=[[lo], [hi]], fmt="o", color=colour, ms=6, mec=ps.SURFACE, mew=1.5,
                        elinewidth=1.6, capsize=0)
        _label(ax, 21.3, y + 0.25, label, va="bottom")
    ax.set_xlim(21, 45)
    ax.set_ylim(-0.6, len(rows) - 0.1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(f"Detector-frame chirp mass ({MSUN})")
    ax.set_title("4  Calibrated vs the full analysis")

    ps.caption(fig, f"{PROVENANCE} Method after Abbott et al. 2017 (Ann. Phys. 529, 1600209): frequency = "
               f"1/(2 × gap between zero crossings); the leading-order chirp law makes f^(-8/3) fall "
               f"linearly to zero at merger with slope ∝ {MC}^(5/3). Half-cycles are kept, walking back from "
               f"the peak, while they stand above the gate (chosen on injections, not on the event). The "
               f"Newtonian law ignores the final orbits' extra speed-up, so it over-estimates {MC}; "
               f"injecting a known waveform into {n_inj} stretches of real noise ({n_inj - len(rec)} fits "
               f"failed) measures that bias "
               f"(×{res['chirp']['bias_factor']:.2f}) and the scatter (68% bars). Detector-frame masses are "
               f"(1+z) times the source's.")
    return ps.save(fig, FIGURE_DIR / "06_chirp_mass")


# --------------------------------------------------------------------------
# 7. Only black holes
# --------------------------------------------------------------------------

def figure_black_holes(res, art):
    """Contact frequency against total mass for every kind of star, and the binary to scale."""
    comp = res["compactness"]
    m_min, f_obs = comp["min_total_mass"], comp["f_gw_used"]
    fig = plt.figure(figsize=(11, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.18)
    ax = fig.add_subplot(gs[0])

    ax.axhspan(10, 2000, color=ps.GRID, alpha=0.5, lw=0, zorder=0)
    _label(ax, 0.12, 12, "LIGO's frequency band", va="bottom", color=ps.MUTED)

    m = np.logspace(np.log10(0.2), np.log10(200), 300)
    f_ms = gp.contact_frequency(m, 2 * gp.radius_main_sequence(m / 2))
    ax.plot(m, f_ms, color=ps.C_CONTEXT, lw=2.0)
    _label(ax, 0.25, f_ms[5] * 3, "two Sun-like stars", va="bottom")

    m = np.linspace(0.3, 2 * gp.M_CHANDRASEKHAR * 0.999, 300)
    f_wd = gp.contact_frequency(m, 2 * gp.radius_white_dwarf(m / 2))
    ax.plot(m, f_wd, color=ps.C_CONTEXT, lw=2.0)
    _label(ax, 0.3, f_wd[0] * 2.2, "two white dwarfs", va="bottom")

    m = np.linspace(2.4, 2 * gp.M_NEUTRON_STAR_MAX, 50)
    f_ns = gp.contact_frequency(m, 2 * gp.R_NEUTRON_STAR)
    ax.plot(m, f_ns, color=ps.C_CONTEXT, lw=2.0)
    _label(ax, 2.2, f_ns[0] * 0.55, "two neutron stars", va="top", ha="right")

    m = np.logspace(0, 3, 300)
    f_bh = gp.contact_frequency(m, gp.schwarzschild_radius(m))
    ax.plot(m, f_bh, color=ps.C_ACCENT, lw=2.4, zorder=4)
    ps.annotate(ax, "two black holes\n(event horizons touching)", xy=(12, gp.contact_frequency(12, gp.schwarzschild_radius(12))),
                xytext=(3.0, 9e4), ha="left", va="center", color=ps.INK)

    ax.add_patch(Rectangle((m_min, f_obs), 3000 - m_min, 3e5 - f_obs, facecolor=ps.C_PAPER,
                           alpha=0.10, edgecolor=ps.C_PAPER, lw=1.2, zorder=1))
    _label(ax, m_min * 1.1, 3.2e4,
           f"GW150914 must reach here:\nmass ≥ {np.floor(m_min):.0f} {MSUN} (chirp mass)\n"
           f"frequency ≥ {f_obs:.0f} Hz (observed)", va="top", color=ps.INK)

    bhns = comp["families"]["black hole + neutron star"]
    ax.plot(bhns["total_mass"], bhns["isco_frequency"], "X", color=ps.INK, ms=8, mec=ps.SURFACE, mew=1.2)
    _label(ax, bhns["total_mass"] / 1.3, bhns["isco_frequency"] / 2.2,
           f"neutron star + black hole with this\nchirp mass: {bhns['total_mass']:.0f} {MSUN} total,"
           f" merges\nat {bhns['isco_frequency']:.0f} Hz", ha="right", va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 3000)
    ax.set_ylim(1e-6, 3e5)
    ax.set_xlabel(f"Total mass of the binary ({MSUN})")
    ax.set_ylabel("Highest gravitational-wave frequency\nbefore the two bodies touch (Hz)")
    ax.set_title("Only one kind of object can make GW150914's signal")

    # To scale: the binary at the moment of peak emission, minimum-mass case.
    ax2 = fig.add_subplot(gs[1])
    sep = comp["separation_km"]
    r_bh = gp.schwarzschild_radius(m_min / 2) / 1e3
    ax2.add_patch(Circle((0, 0), sep / 2, fill=False, edgecolor=ps.C_CONTEXT, lw=1.0, ls=(0, (3, 3))))
    for x in (-sep / 2, sep / 2):
        ax2.add_patch(Circle((x, 0), r_bh, facecolor=ps.INK, edgecolor=ps.SURFACE, lw=2))
    ax2.add_patch(Circle((0, -sep / 2 - 60), gp.R_NEUTRON_STAR / 1e3, facecolor=ps.C_ACCENT,
                         edgecolor=ps.SURFACE, lw=1))
    _label(ax2, 20, -sep / 2 - 60, "a neutron star, to scale", va="center")
    ax2.annotate("", xy=(-sep / 2, r_bh + 22), xytext=(sep / 2, r_bh + 22),
                 arrowprops={"arrowstyle": "<->", "color": ps.INK_SECONDARY, "lw": 0.9})
    _label(ax2, 0, r_bh + 34, f"{sep:.0f} km", ha="center", va="bottom", color=ps.INK)
    _label(ax2, sep / 2, -r_bh - 18, f"horizon\n{2 * r_bh:.0f} km across", ha="center", va="top")
    ax2.plot([-200, -100], [-265, -265], color=ps.INK, lw=2)
    _label(ax2, -150, -275, "100 km", ha="center", va="top")
    ax2.set_xlim(-280, 280)
    ax2.set_ylim(-300, 190)
    ax2.set_aspect("equal")
    ax2.axis("off")
    ax2.set_title(f"To scale, at peak emission ({comp['orbital_frequency']:.0f} orbits/s, "
                  f"{comp['orbital_speed_c']:.2f} c)", loc="left")

    ps.caption(fig, f"{PROVENANCE} Newtonian estimates, deliberately conservative: {MC} ≥ "
               f"{comp['chirp_mass_used']:.1f} {MSUN} (lower 90% bound of our calibrated value) forces a total "
               f"mass ≥ 2^(6/5)·{MC} = {m_min:.1f} {MSUN}; the frequency is the half-cycle just before the "
               f"amplitude peak. Contact frequencies use simple mass–radius laws: main-sequence stars, "
               f"white dwarfs (Nauenberg 1972, up to the Chandrasekhar mass), neutron stars (12 km, "
               f"≤ {gp.M_NEUTRON_STAR_MAX} {MSUN}). The sketch shows the minimum-mass case; heavier holes "
               f"would be even more compact for their size.")
    return ps.save(fig, FIGURE_DIR / "07_only_black_holes")


# --------------------------------------------------------------------------
# 8. Significance
# --------------------------------------------------------------------------

def figure_significance(res, art):
    """Is it real? The event against 187 days of time-slid noise."""
    s = res["significance"]
    fig, ax = plt.subplots(figsize=(10, 5.2))
    panel_background(ax, res, art)
    ax.set_title("Is it real? GW150914 against half a year of noise-only coincidences")
    ax.text(0.205, 0.955,
            f"{s['slides']:,} time slides = {s['background_livetime_days']:.0f} days of background.\n"
            f"Loudest noise coincidence: {s['loudest_background_event_removed']:.1f}; or "
            f"{s['loudest_background']:.1f} if the event's own\nHanford trigger may pair with "
            f"Livingston noise. GW150914: {s['event_statistic']:.1f}.\n\n"
            f"False-alarm rate < {s['far_per_year_bound']:.1f} per year — a bound set only by how\n"
            f"much data we used (p < {s['false_alarm_probability_bound']:.1e}, ≥ {s['sigma_bound']:.1f}σ).\n"
            f"With 16 days of data and 10⁷ slides, the discovery\npaper reached 1 in 203,000 years: 5.1σ.",
            transform=ax.transAxes, fontsize=8, color=ps.INK_SECONDARY, va="top")
    zl = np.asarray(art["significance"]["zero_lag"]["statistic"])
    others = zl[zl < s["event_statistic"]]
    bg = np.asarray(art["significance"]["background"]["event_removed"])
    expect = np.sum(bg >= 5.5) * s["analysed_seconds"] / s["background_livetime_s"]
    ps.caption(fig, f"{PROVENANCE} The real hour's other coincidences are what noise predicts "
                    f"({np.sum(others >= 5.5)} at ≥ 5.5 against {expect:.1f} expected). "
                    f"Single-template search over {s['analysed_seconds'] / 60:.0f} minutes: "
                    f"triggers are SNR peaks above {art['det']['H1']['triggers']['snr'].min():.0f} "
                    f"({s['triggers']['H1']} in H1, {s['triggers']['L1']} in L1), ranked by χ²-reweighted "
                    f"SNR; coincidence window 15 ms (10 ms light travel + 5 ms timing slack) — all as in "
                    f"the discovery paper, except one template instead of ~250,000, which makes this "
                    f"background optimistic. Shifting Livingston by whole seconds destroys any real "
                    f"coincidence; what pairs up is noise with noise.")
    return ps.save(fig, FIGURE_DIR / "08_significance")


# --------------------------------------------------------------------------
# 0. Summary
# --------------------------------------------------------------------------

def figure_summary(res, art):
    """One page: the signal, the chirp, the chirp mass, the significance."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4), gridspec_kw={"hspace": 0.45, "wspace": 0.25})
    panel_overlay(axes[0, 0], res, art)
    axes[0, 0].set_title("Same signal, two detectors")
    mesh = panel_qplane(axes[0, 1], res, art, "H1")
    axes[0, 1].set_title("The chirp: frequency climbing (Hanford)")
    axes[0, 1].set_xlabel(T_LABEL)
    axes[0, 1].set_ylabel("Frequency (Hz)")
    cbar = fig.colorbar(mesh, ax=axes[0, 1], fraction=0.04, pad=0.02)
    cbar.set_label("Normalised energy")
    cbar.outline.set_visible(False)
    panel_chirp_fit(axes[1, 0], res, art, compact=True)
    ch = res["chirp"]
    axes[1, 0].set_title(f"Weighing it: {MC} = {ch['calibrated']:.0f} {MSUN} after calibration")
    panel_background(axes[1, 1], res, art)
    axes[1, 1].set_title("Louder than half a year of noise")
    fig.suptitle("GW150914, reprocessed from the raw LIGO data", x=0.07, ha="left", y=0.955)
    s = res["significance"]
    ps.caption(fig, f"{PROVENANCE} Arrival lag {res['arrival']['lag_mf_ms']:.1f} ± "
                    f"{res['arrival']['lag_mf_err_ms']:.1f} ms (published 6.9); network SNR "
                    f"{res['network_snr']:.1f} (published 23.6–25.2); detector-frame chirp mass "
                    f"{ch['calibrated']:.1f} {MSUN} (68%: {ch['calibrated_68'][0]:.1f}–"
                    f"{ch['calibrated_68'][1]:.1f}; published "
                    f"{res['published']['chirp_mass_detector_frame']['median']:.1f}); louder than every one "
                    f"of {s['slides']:,} time slides. Only black holes are compact enough to have made it.")
    return ps.save(fig, FIGURE_DIR / "00_summary")


def make_all(res, art) -> list:
    written = []
    paths, extra = figure_noise(res, art)
    written += paths
    for fn in (figure_reveal, figure_two_detectors, figure_spectrogram, figure_matched_filter,
               figure_chirp_mass, figure_black_holes, figure_significance, figure_summary):
        written += fn(res, art)
    return written, extra
