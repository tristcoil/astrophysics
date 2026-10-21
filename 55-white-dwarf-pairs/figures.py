"""Every figure, from the cached analysis (``run.py --figures-only`` redraws them)."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

import config as C
import fates
import gw
import plotstyle as S
import population as pop
import timing as T

S.use_style()
FIG = C.FIGURE_DIR
DAY = 86400.0


def _jyear(bjd):
    return 2000.0 + (np.asarray(bjd) - 2451545.0) / 365.25


def _model_curve(b_name, shapes, phi, exptime_s=None, f_orb_per_day=208.2986):
    """Shape model (no group offset) on a phase grid; smeared over an exposure if given."""
    sh = shapes[b_name]
    n_harm = {"zg": 4, "zr": 4, "zi": 3}.get(b_name, 4)
    if exptime_s:
        half = exptime_s / 2 / DAY * f_orb_per_day
        nodes = phi[:, None] + half * T.GL_X[None, :]
        basis = T.shape_basis(nodes, n_harm, False, True)
        basis = np.tensordot(basis, T.GL_W / 2, axes=([1], [0]))
    else:
        basis = T.shape_basis(phi, n_harm, False, True)
    return basis @ sh["coef"]


# ---------------------------------------------------------------------------
# 1. The light curve
# ---------------------------------------------------------------------------


def fig_light_curve(res, art):
    fit = art["ztf_fit"]
    th = art["joint_theta"]
    shapes = art["shapes"]
    d = fit["data"]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), sharey=True)
    grid = np.linspace(-0.6, 0.6, 1201)
    for ax, band in zip(axes, ("zg", "zr", "zi")):
        s = d[d["band"] == band]
        if s.empty:
            continue
        st = S.BAND_STYLE[band]
        off = shapes[band]["offsets"]
        y = s["flux"].to_numpy() - s["group"].map(off).to_numpy()
        phi = T.phase(s["bjd_tdb"].to_numpy() - T.T_REF, th)
        ph = (phi + 0.5) % 1.0 - 0.5
        ax.scatter(ph, y, s=4, color=st["color"], alpha=0.18, lw=0, zorder=1)
        bins = np.linspace(-0.5, 0.5, 41)
        idx = np.digitize(ph, bins) - 1
        w = 1 / s["err"].to_numpy() ** 2
        mb = np.array([np.sum(y[idx == k] * w[idx == k]) / np.sum(w[idx == k]) if np.any(idx == k) else np.nan
                       for k in range(40)])
        eb = np.array([1 / np.sqrt(np.sum(w[idx == k])) if np.any(idx == k) else np.nan for k in range(40)])
        cb = 0.5 * (bins[1:] + bins[:-1])
        ax.errorbar(cb, mb, eb, fmt=st["marker"], ms=4.5, color=st["color"], mec=S.SURFACE, mew=0.6,
                    elinewidth=0.8, zorder=3, label=f"{st['label']}: {len(s):,} exposures (binned)")
        sharp = _model_curve(band, shapes, grid)
        smeared = _model_curve(band, shapes, grid, exptime_s=30.0)
        ax.plot(grid, sharp, color=S.INK_SECONDARY, lw=0.8, zorder=2, label="model, instantaneous")
        ax.plot(grid, smeared, color=S.INK, lw=1.6, zorder=4, label="model, averaged over 30 s")
        ax.set_xlim(-0.5, 0.5)
        ax.set_xlabel("orbital phase (0 = hot white dwarf eclipsed)")
        ax.set_title(st["label"], fontsize=10)
        ax.legend(loc="lower right", fontsize=7)
    axes[0].set_ylabel("flux relative to the orbit average (microjansky)")
    S.annotate(axes[0], "primary eclipse:\nthe cool star hides\nthe hot one", xy=(0.0, -24), xytext=(-0.46, -27))
    S.annotate(axes[0], "irradiated face of\nthe cool star turned\ntowards us", xy=(0.38, 6.5), xytext=(0.12, 14))
    axes[0].plot([0.30, 0.30 + 30 / 414.79], [-30, -30], color=S.INK, lw=2.2, solid_capstyle="butt")
    axes[0].text(0.30 + 15 / 414.79, -28.4, "one 30-s\nexposure", ha="center", va="bottom", fontsize=7,
                 color=S.INK_SECONDARY)
    fig.suptitle("ZTF J1539+5027: a 6.9-minute orbit caught in 30-second snapshots", x=0.01, y=1.03, ha="left")
    S.caption(fig, "Every ZTF exposure 2018-2026 that covers the star, measured by forced PSF photometry on IPAC "
                   "difference images and folded on our own ephemeris (no published timing used). Points: single "
                   "exposures; markers: weighted means in 40 phase bins. The eclipse is ~54 s long, so 30-s "
                   "exposures smear it; the model is integrated over each exposure. Geometry (radii, inclination) "
                   "fixed to Burdge et al. (2019).", y=-0.04)
    return S.save(fig, FIG / "01_light_curve")


# ---------------------------------------------------------------------------
# 2. Finding the period, and watching it shrink
# ---------------------------------------------------------------------------


def fig_period(res, art):
    freq, power = art["periodogram"]
    seasons = art["seasons"]
    th, cov = art["joint_theta"], art["joint_cov"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.9), gridspec_kw={"width_ratios": [1.1, 1]})
    period_min = 1440.0 / freq
    a1.plot(period_min, power, color=S.BLUE, lw=0.5)
    a1.set_xscale("log")
    a1.set_xlim(3, 72)
    a1.set_xticks([3, 5, 7, 10, 20, 30, 60])
    a1.set_xticklabels(["3", "5", "7", "10", "20", "30", "60"])
    a1.xaxis.set_minor_formatter(plt.NullFormatter())
    a1.set_xlabel("trial orbital period (minutes)")
    a1.set_ylabel("Lomb-Scargle power")
    disc = res["discovery"]
    S.annotate(a1, f"P = {disc['P_best_s']:.3f} s = {disc['P_best_s'] / 60:.2f} min\n"
                   f"false-alarm probability {disc['fap']:.0e}",
               xy=(disc["P_best_s"] / 60, disc["power"]), xytext=(12, 0.9 * disc["power"]))
    a1.set_title("Hidden in 2018's survey data: one clear period", fontsize=10)
    # local period per season
    yr = _jyear(seasons["t_mid"])
    P = seasons["P_local_s"].to_numpy()
    Perr = seasons["P_err_s"].to_numpy()
    P0 = 414.79
    a2.errorbar(yr, (P - P0) * 1e3, Perr * 1e3, fmt="o", color=S.BLUE, ms=5, mec=S.SURFACE, mew=0.6,
                elinewidth=0.9, label="period measured within each season", zorder=3)
    tg = np.linspace(seasons["t_mid"].min() - 200, seasons["t_mid"].max() + 200, 200)
    Pm = DAY / (th[1] + th[2] * (tg - T.T_REF))
    a2.plot(_jyear(tg), (Pm - P0) * 1e3, color=S.INK, lw=1.4, label="our ephemeris (P shrinking at constant rate)",
            zorder=2)
    a2.axhline(0, color=S.AXIS, lw=0.6, zorder=0)
    a2.set_xlabel("year")
    a2.set_ylabel(f"orbital period - {P0} s (milliseconds)")
    e = res["joint_ztf_ptf"]
    a2.set_title(f"...and it shrinks: dP/dt = {e['Pdot'] * 1e11:.4f} x 10$^{{-11}}$ s/s", fontsize=10)
    S.annotate(a2, f"{abs(e['Pdot']) * C.YEAR * 1e6:.0f} microseconds\nshorter every year",
               xy=(_jyear(seasons["t_mid"].iloc[-1]), (P[-1] - P0) * 1e3), xytext=(2018.6, (P[-1] - P0) * 1e3 - 0.6))
    a2.legend(loc="upper right", fontsize=7.5)
    fig.suptitle("Finding the orbit, then timing it", x=0.01, y=1.03, ha="left")
    S.caption(fig, "Left: Lomb-Scargle periodogram of the 2018 ZTF g and r photometry (the high-cadence survey "
                   "in which the binary was discovered). Right: a constant period fitted to each observing season "
                   "on its own; the line is the phase-coherent fit to all seasons, which is ~100x more precise.",
              y=-0.04)
    return S.save(fig, FIG / "02_period_search")


# ---------------------------------------------------------------------------
# 3. The orbital-decay (O-C) diagram
# ---------------------------------------------------------------------------


def _oc_reference(th, t_ref):
    """Linear ephemeris tangent to ours at t_ref: returns f_ref (1/day)."""
    return th[1] + th[2] * (t_ref - T.T_REF)


def _model_deviation(th, t, t_ref):
    """Eclipse-time deviation (s) of the quadratic ephemeris from its tangent line at t_ref."""
    x, x0 = t - T.T_REF, t_ref - T.T_REF
    dphi = T.phase(x, th) - (T.phase(x0, th) + _oc_reference(th, t_ref) * (t - t_ref))
    return -dphi / _oc_reference(th, t_ref) * DAY


def _published_deviation(T0, P_d, pdot, th, t, t_ref):
    """Their eclipse-time prediction minus ours, plus our model deviation (s)."""
    f0 = 1 / P_d
    f1 = -pdot / P_d**2
    phi_theirs = f0 * (t - T0) + 0.5 * f1 * (t - T0) ** 2
    phi_ours = T.phase(t - T.T_REF, th) - T.phase(T0 - T.T_REF, th)
    # both count eclipses from their T0 (an eclipse of ours too, to within a second)
    t_ecl_ours, _ = T.predicted_eclipse_time(th, T0)
    phi_ours = phi_ours + (T0 - t_ecl_ours) * (th[1] + th[2] * (T0 - T.T_REF))
    diff_s = -(phi_theirs - phi_ours) / (th[1] + th[2] * (t - T.T_REF)) * DAY
    return diff_s


def fig_oc(res, art):
    th = art["joint_theta"]
    oc = art["oc"]
    t_ref = float(oc["t_mid"].min())  # the first PTF season
    tg = np.linspace(t_ref - 100, oc["t_mid"].max() + 250, 600)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.5, 6.6), sharex=True, gridspec_kw={"height_ratios": [2.6, 1]})
    dev = _model_deviation(th, tg, t_ref)
    a1.plot(_jyear(tg), dev, color=S.INK, lw=1.5, zorder=2, label="our fit: period shrinking at a constant rate")
    for band in ("ptf", "zg", "zr", "zi"):
        s = oc[oc["band"] == band]
        if s.empty:
            continue
        st = S.BAND_STYLE[band]
        y = _model_deviation(th, s["t_mid"].to_numpy(), t_ref) + s["dt_s"].to_numpy()
        a1.errorbar(_jyear(s["t_mid"].to_numpy()), y, s["dt_err_s"].to_numpy(), fmt=st["marker"], color=st["color"], ms=5.5,
                    mec=S.SURFACE, mew=0.6, elinewidth=0.9, zorder=3, label=f"{st['label']} (one point per season)")
        a2.errorbar(_jyear(s["t_mid"].to_numpy()), s["dt_s"].to_numpy(), s["dt_err_s"].to_numpy(), fmt=st["marker"], color=st["color"], ms=5,
                    mec=S.SURFACE, mew=0.6, elinewidth=0.9, zorder=3)
    a1.axhline(0, color=S.AXIS, lw=0.6)
    end = dev[-1]
    n_orbits = -end / (DAY / _oc_reference(th, t_ref))
    S.annotate(a1, f"by {_jyear(tg[-1]):.0f} the eclipses come {-end:,.0f} s\n= {n_orbits:.1f} orbits earlier "
                   "than a clock\nticking at the 2009 period", xy=(_jyear(tg[-1]), end), xytext=(2015.4, 0.93 * end))
    a1.set_ylabel("eclipse arrival time minus\nconstant-period prediction (s)")
    a1.legend(loc="lower left", fontsize=7.5, ncol=1)
    e = res["joint_ztf_ptf"]
    a1.set_title(f"Seventeen years of eclipses: dP/dt = ({e['Pdot'] * 1e11:.4f} $\\pm$ {e['Pdot_err'] * 1e11:.4f})"
                 " x 10$^{-11}$ s/s", fontsize=10)
    # published ephemerides relative to ours
    for name, T0, P_d, pdot, ls, lab in (
            ("burdge", C.BURDGE19["T0_bjd_tdb"][0], C.BURDGE19["P_s"][0] / DAY, C.BURDGE19["Pdot"][0], "-",
             "Burdge et al. 2019 ephemeris"),
            ("shaw", 59341.97656626 + 2400000.5, 0.0048008034992, -2.36812e-11, "--",
             "Burdge's 2024 update (Shaw et al. 2024)")):
        a2.plot(_jyear(tg), _published_deviation(T0, P_d, pdot, th, tg, t_ref), color=S.C_PUBLISHED, lw=1.1,
                ls=ls, label=lab, zorder=1)
    a2.axhline(0, color=S.INK, lw=1.0)
    a2.set_ylabel("minus our fit (s)")
    a2.set_xlabel("year")
    a2.set_xticks(np.arange(2009, 2027, 2))
    a2.legend(loc="lower left", fontsize=7.5)
    lim = np.nanmax(np.abs(oc["dt_s"]) + oc["dt_err_s"]) * 1.25
    a2.set_ylim(-lim, lim)
    fig.suptitle("Gravitational waves are draining the orbit of ZTF J1539+5027", x=0.01, y=1.03, ha="left")
    S.caption(fig, "Each point is the eclipse timing of one observing season and band, measured by fitting the light-"
                   "curve template to that season alone. PTF: our photometry on public Palomar Transient Factory "
                   "images (2009-2012); ZTF: forced photometry on difference images (2018-2026). The parabola is the "
                   "signature of a constant dP/dt, as in the Hulse-Taylor pulsar (C9).", y=-0.03)
    return S.save(fig, FIG / "03_orbital_decay")


# ---------------------------------------------------------------------------
# 4. Uniqueness of the solution
# ---------------------------------------------------------------------------


def fig_landscape(res, art):
    pd_grid, chi = art["landscape_ztf"]
    pred = art["ptf_prediction"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.6), gridspec_kw={"width_ratios": [1.6, 1]})
    a1.plot(pd_grid * 1e11, chi - chi.min(), color=S.BLUE, lw=1.4)
    a1.set_yscale("symlog", linthresh=10)
    a1.set_xlabel("trial dP/dt (10$^{-11}$ s/s, steps of 0.1), period and phase re-fitted at each")
    a1.set_ylabel("$\\Delta\\chi^2$ (ZTF only)")
    best = res["ztf_only"]["Pdot"]
    S.annotate(a1, f"one minimum: the fit gives {best * 1e11:.3f}", xy=(best * 1e11, 0),
               xytext=(best * 1e11 + 1.2, 3))
    a1.axvline(0, color=S.AXIS, lw=0.6)
    a1.text(0.05, 40, "constant period", fontsize=7.5, color=S.INK_SECONDARY, rotation=90, va="bottom")
    a1.set_title("Only one rate of decay fits the ZTF seasons", fontsize=10)
    yrs = pred["year"].to_numpy()
    for k in (-1, 1):
        a2.axhline(k, color=S.AXIS, lw=0.8)
        a2.text(yrs.min() - 0.5, k + 0.05, f"{k:+d} orbit", fontsize=7.5, color=S.INK_SECONDARY)
    a2.axhline(0, color=S.INK, lw=0.8)
    a2.errorbar(yrs, pred["dphi"], np.hypot(pred["dphi_err"], 0), fmt="D", color=S.INK, ms=6, mec=S.SURFACE,
                elinewidth=1.0, zorder=3, label="PTF season (measured)")
    a2.fill_between(np.r_[yrs.min() - 0.6, yrs, yrs.max() + 0.6],
                    -np.r_[pred["prediction_err"].iloc[0], pred["prediction_err"], pred["prediction_err"].iloc[-1]],
                    np.r_[pred["prediction_err"].iloc[0], pred["prediction_err"], pred["prediction_err"].iloc[-1]],
                    color=S.BLUE_150, alpha=0.8, lw=0, label="ZTF-only prediction, 1$\\sigma$")
    a2.set_xlim(yrs.min() - 0.7, yrs.max() + 0.7)
    a2.set_ylim(-1.4, 1.4)
    a2.set_xticks(yrs)
    a2.set_xlabel("PTF season")
    a2.set_ylabel("eclipse phase minus ZTF-only\nprediction (orbits)")
    a2.legend(loc="lower right", fontsize=7)
    a2.set_title("...and ZTF alone predicts\nPTF's eclipses a decade earlier", fontsize=10)
    fig.suptitle("Is the solution unique?", x=0.01, y=1.03, ha="left")
    one = pred["pdot_change_for_one_cycle"].min()
    S.caption(fig, "Left: chi-squared of the ZTF fit when dP/dt is fixed to each trial value and everything else is "
                   "re-fitted. Right: the eclipse phase of each PTF season (2009-2012) against the extrapolated ZTF-only "
                   f"ephemeris. Miscounting one orbit across the 6-year gap would need dP/dt to be off by "
                   f"{one * 1e11:.2f} x 10^-11, many times ZTF's own error.", y=-0.05)
    return S.save(fig, FIG / "04_uniqueness")


# ---------------------------------------------------------------------------
# 5. What the decay says about the masses
# ---------------------------------------------------------------------------


def fig_masses(res, art):
    cm = res["chirp_mass"]
    m1 = np.linspace(0.3, 1.0, 400)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))

    def m2_for_mc(mc):
        out = []
        for a in m1:
            g = np.linspace(0.05, 1.0, 4000)
            v = gw.chirp_mass(a, g) - mc
            i = np.where(np.diff(np.sign(v)))[0]
            out.append(np.interp(0, v[i[0]:i[0] + 2], g[i[0]:i[0] + 2]) if i.size else np.nan)
        return np.array(out)

    lo, mid, hi = (m2_for_mc(cm["Mc_gr_only"] + k * cm["Mc_gr_only_err"] * 3) for k in (-1, 0, 1))
    ax.fill_between(m1, lo, hi, color=S.BLUE_150, alpha=0.7, lw=0, label="chirp mass from our dP/dt if all of it is GR\n(its 3$\\sigma$ band is thinner than the line)")
    ax.plot(m1, mid, color=S.BLUE, lw=1.6)
    ax.plot(m1, m2_for_mc(cm["Mc_if_6p7pct_tides"]), color=S.BLUE, lw=1.1, ls="--",
            label="...if tides add 6.7% to the decay (Burdge et al.'s estimate)")
    b1, b2 = C.BURDGE19["M1"], C.BURDGE19["M2"]
    ax.errorbar([b1[0]], [b2[0]], xerr=[[b1[2]], [b1[1]]], yerr=[[b2[2]], [b2[1]]], fmt="*", ms=13,
                color=S.C_J1539, mec=S.SURFACE, mew=0.6, elinewidth=1.0, zorder=4,
                label="Burdge et al. 2019 masses (spectra + light-curve model + decay)")
    ax.set_xlabel("mass of the hot white dwarf M$_1$ (M$_\\odot$)")
    ax.set_ylabel("mass of the cool white dwarf M$_2$ (M$_\\odot$)")
    ax.set_xlim(0.3, 1.0)
    ax.set_ylim(0.1, 0.45)
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title(f"Chirp mass {cm['Mc_gr_only']:.4f} $\\pm$ {cm['Mc_gr_only_err']:.4f} M$_\\odot$ if the decay is "
                 f"pure GR;\nthe orbit decays {cm['ratio_obs_over_gr']:.2f} $\\pm$ {cm['ratio_err']:.2f} times as fast as GR "
                 "predicts from Burdge et al.'s masses", fontsize=10)
    S.caption(fig, "Gravitational-wave emission alone fixes the chirp mass (M1 M2)^(3/5) / (M1+M2)^(1/5) from P and "
                   "dP/dt. Tidal locking of the stars would take some orbital energy too and speed the decay, so the "
                   "true chirp mass would sit slightly lower (dashed).", y=-0.03)
    return S.save(fig, FIG / "05_chirp_mass")


# ---------------------------------------------------------------------------
# 6. LISA's view of the known binaries
# ---------------------------------------------------------------------------


def fig_lisa(res, art):
    b = art["binaries"]
    T4 = 4 * C.YEAR
    f = np.logspace(-4, -1, 600)
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    popd = art["pop_default"]
    # resolved population (sample)
    rf, rp = popd["resolved_f"], popd["resolved_pw"]
    rng = np.random.default_rng(3)
    k = rng.choice(rf.size, size=min(6000, rf.size), replace=False)
    ax.scatter(rf[k], np.sqrt(4 * rf[k] * T4 * rp[k]), s=2, color=S.CONTEXT, alpha=0.5, lw=0, zorder=1,
               label=f"model Galaxy: {rf.size:,} pairs LISA resolves (sample)")
    ax.plot(f, gw.characteristic_strain_noise(f, 4.0, confusion=False), color=S.C_INSTRUMENT, lw=1.6, zorder=3,
            label="LISA instrument noise (Robson et al. 2019)")
    s_ours = np.interp(f, popd["f"], popd["s_conf"], left=np.nan, right=0.0)
    ax.plot(f, np.sqrt(f * (gw.sn_instrument(f) + s_ours)), color=S.C_FOREGROUND_OURS, lw=1.3, zorder=3,
            label="+ unresolved Galactic pairs (our model)")
    ax.plot(f, gw.characteristic_strain_noise(f, 4.0, confusion=True), color=S.C_FOREGROUND_ROBSON, lw=1.1, ls="--",
            zorder=3, label="+ foreground in Robson et al.'s fit")
    for typ, st in S.TYPE_STYLE.items():
        for status, filled in (("verification", True), ("detectable", False)):
            s = b[(b["type"] == typ) & (b["status"] == status)]
            if s.empty:
                continue
            fac = np.where(s["incl_kind"] == "assumed", 0.8,
                           gw.inclination_factor(np.radians(s["incl"].to_numpy(float))))
            fs = s["f_Hz"].to_numpy(float)
            h = np.sqrt(4 * fac * fs * T4) * s["A"].to_numpy(float)
            lo = np.sqrt(4 * fac * fs * T4) * s["A_p16"].to_numpy(float)
            hi = np.sqrt(4 * fac * fs * T4) * s["A_p84"].to_numpy(float)
            ax.errorbar(fs, h, [h - lo, hi - h], fmt=st["marker"], ms=7, color=st["color"],
                        mfc=st["color"] if filled else S.SURFACE, mec=st["color"] if not filled else S.SURFACE,
                        mew=1.2 if not filled else 0.6, elinewidth=0.8, zorder=4,
                        label=f"{st['label']}, {'verification' if filled else 'detectable'}")
    j = b[b["name"] == "ZTFJ1539"].iloc[0]
    hj = float(np.sqrt(4 * gw.inclination_factor(np.radians(84.15)) * j["f_Hz"] * T4) * j["A"])
    ax.plot(float(j["f_Hz"]), hj, "*", ms=17, color=S.C_J1539, mec=S.SURFACE, mew=0.8, zorder=6)
    S.annotate(ax, f"ZTF J1539+5027\nSNR {j['snr_4yr']:.0f} in 4 yr", xy=(float(j["f_Hz"]), hj),
               xytext=(9e-3, 4e-19))
    # a few labels with leader lines (label position as a multiple of the point's f and h)
    for name, label, mx, my in (("HM Cnc", "HM Cnc", 1.5, 2.2), ("V407 Vul", "V407 Vul", 1.25, 3.2),
                                ("SDSSJ0651", "SDSS J0651", 1.9, 0.30),
                                ("4U1820-30", "4U 1820-30 (X-ray binary)", 1.6, 0.6),
                                ("HD265435", "HD 265435", 1.0, 0.45)):
        r = b[b["name"] == name]
        if r.empty:
            continue
        r = r.iloc[0]
        fac = 0.8 if r["incl_kind"] == "assumed" else float(gw.inclination_factor(np.radians(r["incl"])))
        hh = float(np.sqrt(4 * fac * r["f_Hz"] * T4) * r["A"])
        S.annotate(ax, label, xy=(float(r["f_Hz"]), hh), xytext=(float(r["f_Hz"]) * mx, hh * my))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-4, 1e-1)
    ax.set_ylim(1e-21, 2e-18)
    ax.set_xlabel("gravitational-wave frequency (Hz) = 2 / orbital period")
    ax.set_ylabel("characteristic strain after 4 years")
    ax.legend(loc="upper right", fontsize=6.6, ncol=2)
    ax.set_title("LISA's guaranteed sources: the known binaries against the sensitivity curve", fontsize=10)
    S.caption(fig, "Each binary is plotted at sqrt(4 F f T) A, the noise at sqrt(f S_n), so the height ratio is the "
                   "sky-averaged signal-to-noise ratio (F: inclination factor). Binaries from Kupfer et al. (2024) with "
                   "their masses, inclinations and Gaia-based distances; error bars from those uncertainties. Grey: "
                   "pairs our observationally driven Galaxy model predicts LISA will resolve; our foreground (violet) ends where "
                   "most frequency bins hold no unresolved pair.", y=-0.02)
    return S.save(fig, FIG / "06_lisa_sensitivity")


def fig_snr_growth(res, art):
    b = art["binaries"]
    ver = b[b["status"] == "verification"].copy()
    t_days = np.logspace(np.log10(3), np.log10(4 * 365.25), 80)
    ls = res["lisa_summary"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    for ax, with_fg, title in (
            (axes[0], False, f"Against LISA's instrument noise alone: "
                             f"{ls['n_verification_with_snr3mo_gt7_instrument_only']} of 18 by 3 months"),
            (axes[1], True, f"With the Galactic foreground (Robson et al.): "
                            f"{ls['n_verification_with_snr3mo_gt7']} of 18 by 3 months")):
        for r in ver.to_dict("records"):
            incl = None if r["incl_kind"] == "assumed" else np.radians(r["incl"])
            if with_fg:
                snr = [gw.snr(r["f_Hz"], r["A"], t * DAY, incl=incl, t_obs_yr=max(t / 365.25, 0.5)) for t in t_days]
            else:
                snr = [gw.snr(r["f_Hz"], r["A"], t * DAY, incl=incl, sn=gw.sn_instrument(r["f_Hz"])) for t in t_days]
            is_j = r["name"] == "ZTFJ1539"
            st = S.TYPE_STYLE[r["type"]]
            ax.plot(t_days, snr, color=S.C_J1539 if is_j else st["color"], lw=2.0 if is_j else 0.9,
                    alpha=1.0 if is_j else 0.75, zorder=3 if is_j else 2)
            if is_j or r["name"] in ("SDSSJ0935", "SDSSJ1351"):
                ax.text(t_days[-1] * 1.05, snr[-1], r["name"].replace("SDSSJ", "SDSS J").replace("ZTFJ", "ZTF J"),
                        fontsize=7, va="center", color=S.INK if is_j else S.INK_SECONDARY)
        ax.axhline(7, color=S.INK_SECONDARY, lw=0.8)
        ax.text(3.3, 7.6, "SNR 7", fontsize=7.5, color=S.INK_SECONDARY)
        ax.axvline(91.3, color=S.AXIS, lw=0.8)
        ax.text(96, 0.32, "3 months", fontsize=7.5, color=S.INK_SECONDARY)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(3, 4 * 365.25 * 3.2)
        ax.set_ylim(0.25, 400)
        ax.set_xlabel("LISA observing time (days)")
        ax.set_title(title, fontsize=9.5)
    axes[0].set_ylabel("sky-averaged signal-to-noise ratio")
    handles = [Line2D([], [], color=S.TYPE_STYLE[t]["color"], lw=1, label=S.TYPE_STYLE[t]["label"])
               for t in ("AMCVn", "DWD")] + [Line2D([], [], color=S.C_J1539, lw=2, label="ZTF J1539+5027")]
    axes[0].legend(handles=handles, loc="upper left", fontsize=7.5)
    fig.suptitle("How fast LISA picks up Kupfer et al.'s 18 verification binaries", x=0.01, y=1.03, ha="left")
    S.caption(fig, "SNR grows as the square root of observing time. Below ~3 mHz the unresolved Galactic pairs add "
                   "noise that fades as LISA resolves more of them (Robson et al.'s fits for 0.5-4 years; the 0.5-year "
                   "fit stands in for shorter times, so the right panel is pessimistic early on). Kupfer et al.'s "
                   "Bayesian test, which calls all 18 detectable within 3 months, sits between the two panels.", y=-0.03)
    return S.save(fig, FIG / "07_snr_growth")


# ---------------------------------------------------------------------------
# 8. Fates
# ---------------------------------------------------------------------------


def fig_fates(res, art):
    fat = art["fates"]
    b = art["binaries"]
    m1g = np.linspace(0.12, 1.4, 90)
    lo, hi = fates.boundary_curves(m1g)
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    ax.fill_between(m1g, 0, lo, color=S.BLUE_150, alpha=0.45, lw=0)
    ax.fill_between(m1g, hi, m1g, color=S.WASH, alpha=1.0, lw=0)
    ax.plot(m1g, lo, color=S.BLUE_450, lw=1.0)
    ax.plot(m1g, hi, color=S.INK_SECONDARY, lw=1.0)
    ax.plot(m1g, m1g, color=S.AXIS, lw=0.8)
    ax.text(1.05, 0.08, "stable transfer:\nbecomes an AM CVn", fontsize=8, color=S.BLUE_550)
    ax.text(0.47, 0.42, "unstable:\nthe stars merge", fontsize=8, color=S.INK_SECONDARY, rotation=38)
    ax.text(0.98, 0.36, "depends on how tides\ncouple spin and orbit", fontsize=8, color=S.INK_SECONDARY, rotation=26)
    mc_line = np.linspace(0.7, 1.4, 50)
    ax.plot(mc_line, 1.4 - mc_line, color=S.INK, lw=0.8, ls=":")
    ax.text(0.92, 0.50, "total mass 1.4 M$_\\odot$", fontsize=7.5, color=S.INK, rotation=-38)
    det = fat[fat["type"] == "DWD"]
    sc = ax.scatter(det["m_acc"], det["m_don"], c=np.log10(det["t_contact_yr"]), cmap="Blues_r", vmin=5, vmax=10.5,
                    s=60, marker="s", edgecolor=S.INK_SECONDARY, linewidth=0.6, zorder=4)
    cb = fig.colorbar(sc, ax=ax, pad=0.01)
    cb.set_label("log$_{10}$ years until the pair touches")
    sdb = fat[fat["type"] == "sdB"]
    ax.scatter(sdb["m_acc"], sdb["m_don"], marker="^", s=60, facecolor=S.SURFACE, edgecolor=S.AQUA, lw=1.2, zorder=4,
               label="hot subdwarf + white dwarf (Kupfer et al.)")
    am = b[b["type"] == "AMCVn"]
    ax.scatter(np.maximum(am["m1"], am["m2"]), np.minimum(am["m1"], am["m2"]), marker="o", s=26,
               facecolor=S.SURFACE, edgecolor=S.BLUE, lw=1.0, zorder=3, label="AM CVn: already transferring mass")
    j = fat[fat["name"] == "ZTFJ1539"].iloc[0]
    ax.plot(j["m_acc"], j["m_don"], "*", ms=17, color=S.C_J1539, mec=S.SURFACE, zorder=6)
    S.annotate(ax, f"ZTF J1539+5027: touches in ~{j['t_contact_yr'] / 1e3:,.0f} kyr\n(cold-core radius; its warm "
                   "envelope\nalready fills 96% of its Roche lobe)", xy=(float(j["m_acc"]), float(j["m_don"])),
               xytext=(0.15, 0.012))
    ax.scatter([], [], marker="s", s=50, color=S.BLUE_450, edgecolor=S.INK_SECONDARY, label="detached pair (Kupfer et al.)")
    ax.set_xlim(0.12, 1.4)
    ax.set_ylim(0.0, 0.62)
    ax.set_xlabel("more massive white dwarf (accretor), M$_\\odot$")
    ax.set_ylabel("less massive white dwarf (donor), M$_\\odot$")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("When the pairs touch: merge, or settle into an AM CVn?", fontsize=10)
    S.caption(fig, "Stability boundaries from Marsh, Nelemans & Steeghs (2004) for zero-temperature white dwarfs "
                   "(Eggleton's mass-radius relation): mass transfer is stable below the blue line whatever the tides "
                   "do, unstable above the grey one. Pairs above the dotted line weigh more than the Chandrasekhar mass.",
              y=-0.03)
    return S.save(fig, FIG / "08_fates")


# ---------------------------------------------------------------------------
# 9. The Galactic foreground
# ---------------------------------------------------------------------------


def fig_galaxy(res, art):
    popd = art["pop_default"]
    kor = art["pop_korol22_check"]
    sample = popd["sample"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.9), gridspec_kw={"width_ratios": [1, 1.25]})
    # face-on map of the sample: positions are not kept, so redraw positions for the sample size
    rng = np.random.default_rng(7)
    x, y, z = pop.sample_positions(200_000, rng, "mcmillan17")
    h = a1.hist2d(x, y, bins=180, range=[[-16, 16], [-16, 16]], cmap="Blues", norm=plt.matplotlib.colors.LogNorm(),
                  cmin=1)
    a1.plot(pop.MCMILLAN17["R0_kpc"], 0, "*", ms=12, color=S.ORANGE, mec=S.SURFACE)
    a1.text(pop.MCMILLAN17["R0_kpc"] + 0.6, 0.6, "Sun", color=S.INK, fontsize=8)
    a1.set_aspect("equal")
    a1.set_xlabel("x (kpc)")
    a1.set_ylabel("y (kpc)")
    a1.grid(False)
    p = res["population"]["default"]
    a1.set_title(f"{p['n_in_band'] / 1e6:.0f} million white-dwarf pairs\nin the LISA band (model)", fontsize=10)
    f = popd["f"]
    a2.plot(f, np.sqrt(popd["s_inst"]), color=S.C_INSTRUMENT, lw=1.5, label="LISA instrument noise")
    # all sources (before any are resolved): mean power in log-spaced bins, for display
    edges = np.logspace(-4, -2, 161)
    k = np.digitize(f, edges) - 1
    fm = np.array([np.sqrt(f[k == i].min() * f[k == i].max()) if np.any(k == i) else np.nan for i in range(160)])
    sm = np.array([popd["s_all_mean"][k == i].mean() if np.any(k == i) else np.nan for i in range(160)])
    a2.plot(fm, np.sqrt(sm), color=S.C_FOREGROUND_OURS, lw=1.0, alpha=0.45,
            label="all pairs, before any are resolved (mean)")
    a2.plot(f, np.sqrt(np.where(popd["s_conf"] > 0, popd["s_conf"], np.nan)), color=S.C_FOREGROUND_OURS, lw=1.6,
            label=f"unresolved foreground, our model (4 yr)")
    a2.plot(kor["f"], np.sqrt(np.where(kor["s_conf"] > 0, kor["s_conf"], np.nan)), color=S.C_FOREGROUND_OURS, lw=1.0,
            ls="--", label="same, Korol et al. 2022 set-up (check)")
    a2.plot(f, np.sqrt(gw.sn_confusion(f, 4.0)), color=S.C_FOREGROUND_ROBSON, lw=1.2, label="Robson et al. 2019 fit "
            "(population synthesis)")
    a2.set_xscale("log")
    a2.set_yscale("log")
    a2.set_xlim(1e-4, 1e-2)
    a2.set_ylim(1e-21, 1e-16)
    a2.set_xlabel("gravitational-wave frequency (Hz)")
    a2.set_ylabel("amplitude spectral density (1/$\\sqrt{\\rm Hz}$)")
    a2.legend(loc="upper right", fontsize=7.2)
    a2.set_title(f"A hum of unresolved pairs below {p['foreground_edge_Hz'] * 1e3:.1f} mHz; "
                 f"{p['n_resolved_snr7_4yr']:,} resolved", fontsize=10)
    fig.suptitle("The Galaxy's white-dwarf pairs as LISA will hear them (observationally driven model)", x=0.01,
                 ha="left")
    S.caption(fig, "Model: the 40 pc white-dwarf census (O'Brien et al. 2024) scaled to the McMillan (2017) Galaxy; "
                   "9.5% of white dwarfs in pairs closer than 4 au with separations from Maoz, Hallakoun & Badenes "
                   "(2018); gravitational-wave inspiral over 10 Gyr of constant star formation. Sources are resolved "
                   "iteratively at SNR > 7 (4 years); the foreground is the running median of the rest.", y=-0.03)
    return S.save(fig, FIG / "09_galactic_foreground")


# ---------------------------------------------------------------------------
# 10. The PTF seasons
# ---------------------------------------------------------------------------


def fig_ptf(res, art):
    joint = art["joint"]
    th = art["joint_theta"]
    shapes = art["shapes"]
    pt = art.get("ptf_data")
    if pt is None:
        return []
    years = sorted(pt["year"].unique())
    fig, axes = plt.subplots(1, len(years), figsize=(3.6 * len(years), 3.4), sharey=True)
    grid = np.linspace(-0.5, 1.5, 801)
    for ax, yr in zip(np.atleast_1d(axes), years):
        s = pt[pt["year"] == yr]
        b = T.make_band(s, "ptf", n_harm=4)
        beta, model, _ = T.solve_band(b, th, shapes["zr"]["coef"])
        off, scale = beta[0], beta[-1]
        phi = T.phase(s["bjd_tdb"].to_numpy() - T.T_REF, th) % 1.0
        y = s["flux"].to_numpy() - off
        for shift in (0, 1):
            ax.errorbar(phi + shift, y, s["err"], fmt="D", ms=4, color=S.INK, mec=S.SURFACE, mew=0.5,
                        elinewidth=0.7, alpha=0.85 if shift == 0 else 0.35)
        curve = scale * _model_curve("zr", shapes, grid, exptime_s=60.0)
        ax.plot(grid, curve, color=S.ORANGE, lw=1.4)
        ax.set_xlim(-0.05, 1.55)
        ax.set_title(f"PTF {yr}: {len(s)} images", fontsize=10)
        ax.set_xlabel("orbital phase")
    np.atleast_1d(axes)[0].set_ylabel("flux minus season mean (microjansky)")
    fig.suptitle("A decade earlier: the same star in Palomar Transient Factory images", x=0.01, y=1.03, ha="left")
    S.caption(fig, "Aperture photometry on public PTF images (60-s R-band exposures), calibrated on Gaia stars in each "
                   "image, folded on our final ephemeris. Orange: the ZTF r-band light curve averaged over 60 s, "
                   "scaled per season. Repeated over two cycles for clarity.", y=-0.04)
    return S.save(fig, FIG / "10_ptf_seasons")


def make_all(res, art, animation_too: bool = False) -> list:
    """Every figure and the interactive page; the movie only if missing or asked for."""
    written = []
    for fn in (fig_light_curve, fig_period, fig_oc, fig_landscape, fig_masses, fig_lisa, fig_snr_growth, fig_fates,
               fig_galaxy, fig_ptf):
        try:
            written += fn(res, art)
        except Exception as exc:  # report and continue: one broken figure should not hide the rest
            import traceback

            print(f"figure {fn.__name__} failed: {exc}")
            traceback.print_exc()
    try:
        import interactive

        written.append(interactive.build(res, art))
    except Exception as exc:
        import traceback

        print(f"interactive page failed: {exc}")
        traceback.print_exc()
    import shutil

    movie = FIG / "j1539_orbit.mp4"
    if (animation_too or not movie.exists()) and shutil.which("ffmpeg"):
        import animation

        written.append(animation.make(res, art))
    return written
