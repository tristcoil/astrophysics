"""Every figure, from the cached analysis (``run.py`` calls ``make_all``)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import config as cfg
import cosmology as cosmo
import plotstyle as ps

ps.use_style()
FIG = cfg.FIGURE_DIR
METHOD_LABEL = {"D": "4000 Å break", "F": "full-spectrum fits", "L": "Lick indices"}
BLUE_FILL_95, BLUE_FILL_68 = ps.BLUE_150, ps.BLUE_300


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _cc_points(ax, cc, y, yerr, *, size=4.5, zorder=5, label=True):
    """Chronometer points: ink, one marker shape per dating method."""
    handles = []
    for m, mk in ps.METHOD_MARKERS.items():
        sel = cc["method"] == m
        if not sel.any():
            continue
        face = ps.C_DATA if m == "D" else ps.SURFACE
        ax.errorbar(cc["z"][sel], y[sel], yerr=yerr[sel], fmt=mk, ms=size, mfc=face,
                    mec=ps.C_DATA, mew=0.9, ecolor=ps.C_DATA, elinewidth=0.8, capsize=0,
                    zorder=zorder)
        handles.append(Line2D([], [], marker=mk, ls="none", mfc=face, mec=ps.C_DATA,
                              ms=size, label=f"{METHOD_LABEL[m]} ({sel.sum()})"))
    return handles


def _band(ax, x, bands, color_line=ps.C_OURS, label=None, zorder=2):
    """bands = percentiles (2.3, 15.9, 50, 84.1, 97.7) along axis 0."""
    ax.fill_between(x, bands[0], bands[4], color=BLUE_FILL_95, lw=0, alpha=0.55, zorder=zorder)
    ax.fill_between(x, bands[1], bands[3], color=BLUE_FILL_68, lw=0, alpha=0.75, zorder=zorder)
    ax.plot(x, bands[2], color=color_line, lw=1.6, zorder=zorder + 1, label=label)


def _model_lines(ax, zz, d1, d3, *, kind="H", textbook=True):
    """D1 and D3 LCDM, plus matter-only and empty universes with D3's H0.
    Returns (line, linestyle) pairs: matplotlib reports custom dash patterns
    back as '--', so the pattern is passed on explicitly."""
    h3, om3 = d3["H0"][0], d3["Om"][0]
    h1, om1 = d1["H0"][0], d1["Om_calibrated"][0]
    curves = [(cosmo.hubble(zz, h3, om3), ps.C_CMB, "-", f"D3: Planck CMB, ΛCDM (H0 = {h3:.1f})"),
              (cosmo.hubble(zz, h1, om1), ps.C_SN, "-",
               f"D1: supernovae + Cepheids, ΛCDM (H0 = {h1:.1f})")]
    if textbook:
        curves += [(h3 * (1 + zz) ** 1.5, ps.C_OTHER_MODEL, (0, (5, 2)), "matter only (same H0)"),
                   (h3 * (1 + zz), ps.C_OTHER_MODEL, (0, (1, 1.5)), "empty, coasting (same H0)")]
    out = []
    for H, col, ls_, lab in curves:
        y = H if kind == "H" else H / (1 + zz)
        (ln,) = ax.plot(zz, y, color=col, ls=ls_, lw=1.3, zorder=3, label=lab)
        out.append((ln, ls_))
    return out


def _lookback_axis(ax, d3):
    """Top axis: lookback time in D3's LCDM (a second scale for the same x)."""
    zz = np.linspace(0, 6, 3001)
    lb = cosmo.lookback_time(zz, d3["H0"][0], d3["Om"][0])
    sec = ax.secondary_xaxis("top", functions=(lambda z: np.interp(z, zz, lb),
                                               lambda t: np.interp(t, lb, zz)))
    sec.set_xlabel("how long ago (billion years, in D3's ΛCDM)", color=ps.INK_SECONDARY)
    sec.tick_params(colors=ps.MUTED, labelsize=8)
    return sec


def _forest(ax, rows, *, ref_bands=(), xlabel="", xlim=None):
    """Rows of (label, centre, minus, plus, colour, kind); kind 'ours'/'published'."""
    y = np.arange(len(rows))[::-1]
    for (lab, c, lo, hi, col, kind), yy in zip(rows, y):
        if c is None:
            continue
        mfc = col if kind == "ours" else ps.SURFACE
        ax.errorbar([c], [yy], xerr=[[lo], [hi]], fmt="o", ms=5, mfc=mfc, mec=col, mew=1.1,
                    ecolor=col, elinewidth=1.4, capsize=0, zorder=4)
    for lo, hi, col, lab in ref_bands:
        ax.axvspan(lo, hi, color=col, alpha=0.16, lw=0, zorder=1)
        ax.text(0.5 * (lo + hi), len(rows) - 0.35, lab, color=col, ha="center", va="bottom",
                fontsize=7.5)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8, color=ps.INK)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.set_ylim(-0.7, len(rows) - 0.1 + 0.6)
    ax.set_xlabel(xlabel)
    if xlim:
        ax.set_xlim(*xlim)


def _save(fig, name):
    return ps.save(fig, FIG / name)


# ---------------------------------------------------------------------------
# 01: H(z)
# ---------------------------------------------------------------------------


def fig_hubble_rate(r, art):
    cc, cov = art["cc"], art["covs"]["full"]
    d1, d3 = r["inputs"]["D1"], r["inputs"]["D3"]
    zs = art["ZS"]
    g = art["gp"]["SE kernel, full covariance"]
    sd = np.sqrt(np.diag(cov))
    desi = art["desi"]
    zz = np.linspace(0, cfg.Z_MAX_PLOT, 400)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw={"wspace": 0.22})

    _band(ax, zs, g["bands"]["H"])
    lines = _model_lines(ax, zz, d1, d3)
    hd = _cc_points(ax, cc, cc["H"], sd)
    ax.errorbar(desi["z"], desi["H"], yerr=desi["err"], fmt="^", ms=5.5, mfc=ps.C_DESI,
                mec=ps.SURFACE, mew=0.8, ecolor=ps.C_DESI, elinewidth=1.2, zorder=6)
    hdesi = Line2D([], [], marker="^", ls="none", mfc=ps.C_DESI, mec=ps.SURFACE, ms=6,
                   label="DESI DR2 BAO, scaled by D3's sound horizon")
    ax.set_xlim(0, cfg.Z_MAX_PLOT)
    ax.set_ylim(40, 330)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("expansion rate H(z)  (km/s per Mpc)")
    ax.set_title("The expansion rate, measured epoch by epoch")

    # right: relative to the CMB's LCDM
    ref = lambda z_: cosmo.hubble(z_, d3["H0"][0], d3["Om"][0])
    _band(ax2, zs, 100 * (g["bands"]["H"] / ref(zs) - 1))
    for ln, ls_ in lines:
        x, y = ln.get_data()
        ax2.plot(x, 100 * (y / ref(x) - 1), color=ln.get_color(), ls=ls_, lw=1.3)
    for m, mk in ps.METHOD_MARKERS.items():
        sel = cc["method"] == m
        face = ps.C_DATA if m == "D" else ps.SURFACE
        ax2.errorbar(cc["z"][sel], 100 * (cc["H"][sel] / ref(cc["z"][sel]) - 1),
                     yerr=100 * sd[sel] / ref(cc["z"][sel]), fmt=mk, ms=4.5, mfc=face,
                     mec=ps.C_DATA, mew=0.9, ecolor=ps.C_DATA, elinewidth=0.8, zorder=5)
    ax2.errorbar(desi["z"], 100 * (desi["H"] / ref(desi["z"]) - 1),
                 yerr=100 * desi["err"] / ref(desi["z"]), fmt="^", ms=5.5, mfc=ps.C_DESI,
                 mec=ps.SURFACE, mew=0.8, ecolor=ps.C_DESI, elinewidth=1.2, zorder=6)
    ax2.axhline(0, color=ps.AXIS, lw=0.8)
    ax2.set_xlim(0, cfg.Z_MAX_PLOT)
    ax2.set_ylim(-45, 60)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("difference from the CMB's ΛCDM (%)")
    ax2.set_title("The same, relative to the CMB's prediction")
    z_m = 0.45
    ps.annotate(ax2, "matter only: far faster\nin the past", (z_m, 100 * ((1 + z_m) ** 1.5 * d3["H0"][0]
                / ref(z_m) - 1)), (0.75, 48))
    ps.annotate(ax2, "DESI: ~1–2% precision", (0.934, 100 * (desi["H"][2] / ref(0.934) - 1)),
                (1.3, -30), color=ps.C_DESI)
    band_h = Patch(color=BLUE_FILL_68, label="chronometers alone, smooth reconstruction (68%, 95%)")
    fig.legend(handles=hd + [hdesi, band_h] + [ln for ln, _ in lines], loc="upper center",
               bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=7.8)
    ps.caption(fig, "Chronometer error bars: the diagonal of the full covariance (statistical, "
               "metallicity and star-formation-history errors, plus a stellar-population-model "
               "systematic of 2–5% that is fully correlated between points). The blue band is a "
               "Gaussian-process reconstruction from the 32 chronometers alone. DESI's D_H/r_d "
               "converted with r_d = 146.9 Mpc from our own fit to Planck (D3). Data: Moresco "
               "et al. 2022 compilation; DESI DR2 (2025).", y=-0.17)
    return _save(fig, "01_hubble_rate")


# ---------------------------------------------------------------------------
# 02: the speedometer
# ---------------------------------------------------------------------------


def fig_speedometer(r, art):
    cc, cov = art["cc"], art["covs"]["full"]
    d1, d3 = r["inputs"]["D1"], r["inputs"]["D3"]
    zs = art["ZS"]
    g = art["gp"]["SE kernel, full covariance"]
    sd = np.sqrt(np.diag(cov)) / (1 + cc["z"])
    desi = art["desi"]
    zz = np.linspace(0, cfg.Z_MAX_PLOT, 400)
    sp = r["speedometer"]
    fig, ax = plt.subplots(figsize=(13.5, 6.2))
    fig.subplots_adjust(right=0.70)

    zt1, szt1 = d1["z_transition"]
    ax.axvspan(zt1 - szt1, zt1 + szt1, color=ps.C_SN, alpha=0.12, lw=0, zorder=0)
    _band(ax, zs, g["bands"]["adot"])
    lines = _model_lines(ax, zz, d1, d3, kind="adot")
    hd = _cc_points(ax, cc, cc["H"] / (1 + cc["z"]), sd)
    ax.errorbar(desi["z"], desi["H"] / (1 + desi["z"]), yerr=desi["err"] / (1 + desi["z"]),
                fmt="^", ms=6, mfc=ps.C_DESI, mec=ps.SURFACE, mew=0.8, ecolor=ps.C_DESI,
                elinewidth=1.2, zorder=6)
    hdesi = Line2D([], [], marker="^", ls="none", mfc=ps.C_DESI, mec=ps.SURFACE, ms=6,
                   label="DESI DR2 BAO (with D3's sound horizon)")
    zt = r["gp"]["SE kernel, full covariance"]["z_transition"]
    ax.errorbar([zt["median"]], [38.5], xerr=[[zt["minus"]], [zt["plus"]]], fmt="o", ms=6,
                color=ps.C_OURS, elinewidth=2, capsize=0, zorder=7, clip_on=False)
    ax.text(zt["median"] + zt["plus"] + 0.03, 38.5,
            f"chronometers alone: the dip is at z = {zt['median']:.2f} ± "
            f"{0.5 * (zt['plus'] + zt['minus']):.2f}", va="center", fontsize=8, color=ps.C_OURS)
    box = {"boxstyle": "round,pad=0.25", "fc": ps.SURFACE, "ec": "none", "alpha": 0.9}
    ax.text(0.02, 97.5, "recently: speeding up\n(dark energy wins)", ha="left", va="top",
            fontsize=8.5, color=ps.INK_SECONDARY, bbox=box, zorder=8)
    ax.text(2.48, 97.5, "earlier: slowing down\n(gravity wins)", ha="right", va="top",
            fontsize=8.5, color=ps.INK_SECONDARY, bbox=box, zorder=8)
    ax.set_xlim(0, cfg.Z_MAX_PLOT)
    ax.set_ylim(35, 100)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("speed at which two galaxies 1 Mpc apart today\nwere separating then:  "
                  "da/dt = H(z)/(1+z)  (km/s)")
    ax.set_title("The universe's speedometer: the expansion slowed, then sped up", pad=34)
    _lookback_axis(ax, d3)
    band_h = Patch(color=BLUE_FILL_68, label="chronometers alone (68%, 95%)")
    d1_band = Patch(color=ps.C_SN, alpha=0.25,
                    label=f"D1: acceleration began at z = {zt1:.2f} ± {szt1:.2f}")
    leg1 = ax.legend(handles=hd + [hdesi], loc="upper left", bbox_to_anchor=(1.02, 1.0),
                     fontsize=8, title="measurements", title_fontsize=8.5, alignment="left")
    ax.add_artist(leg1)
    ax.legend(handles=[band_h] + [ln for ln, _ in lines] + [d1_band], loc="upper left",
              bbox_to_anchor=(1.02, 0.62), fontsize=8, title="models", title_fontsize=8.5,
              alignment="left")
    ps.caption(fig, f"Reading the needle: in D3's ΛCDM two galaxies 1 Mpc apart today separate at "
               f"{sp['adot_today']:.0f} km/s now, were separating at only "
               f"{sp['adot_at_transition']:.0f} km/s at the slowest point "
               f"(z = {sp['z_transition']:.2f}, {sp['lookback_to_transition_gyr']:.1f} billion "
               f"years ago), and at {sp['adot_at_z2']:.0f} km/s at z = 2. A dip means the "
               "expansion decelerated and then accelerated. A matter-only universe only ever "
               "slows down, so its separation speed was always higher in the past (rising to the "
               "right); an empty one coasts (flat). Chronometer data: Moresco et al. 2022 "
               "compilation; error bars from the diagonal of the full covariance.")
    return _save(fig, "02_speedometer")


# ---------------------------------------------------------------------------
# 03: the transition and the deceleration parameter
# ---------------------------------------------------------------------------


def fig_transition(r, art):
    d1, d3 = r["inputs"]["D1"], r["inputs"]["D3"]
    fl = r["flat_lcdm"]["full"]["z_transition"]
    ol = r["open_lcdm"]["z_transition"]
    gpz = r["gp"]["SE kernel, full covariance"]
    pw = r["piecewise"]
    m16 = r["moresco2016_replication"]
    pub = cfg.PUBLISHED

    def prof(p):
        b = p["best"]
        return b, b - p["68"][0], p["68"][1] - b

    rows = [
        ("Flat ΛCDM fit", fl["median"], fl["minus"], fl["plus"], ps.C_OURS, "ours"),
        ("Curvature free (ΛCDM)", ol["median"], ol["minus"], ol["plus"], ps.C_OURS, "ours"),
        (f"Smooth reconstruction (GP; {100 * gpz['frac_transition_in_data']:.0f}% of curves)",
         gpz["z_transition"]["median"], gpz["z_transition"]["minus"],
         gpz["z_transition"]["plus"], ps.C_OURS, "ours"),
        ("Two straight lines", *prof(pw["none"]["profile"]), ps.C_OURS, "ours"),
        (f"Two lines + D1's H0 = {d1['H0'][0]:.1f}",
         *prof(pw["D1 H0 (Cepheids + supernovae)"]["profile"]), ps.C_SN, "ours"),
        (f"Two lines + D3's H0 = {d3['H0'][0]:.1f}",
         *prof(pw["D3 H0 (Planck CMB in LCDM)"]["profile"]), ps.C_CMB, "ours"),
        ("Moresco+2016, two lines + H0 = 73 ± 2.4 (published)", pub["m16_zt_piecewise"][0],
         pub["m16_zt_piecewise"][2], pub["m16_zt_piecewise"][1], ps.INK_SECONDARY, "published"),
        ("  redone here, same data and prior", *prof(m16["piecewise_profile"]), ps.INK_SECONDARY, "ours"),
        ("Moresco+2016, ΛCDM + H0 prior (published)", pub["m16_zt_olcdm"][0],
         pub["m16_zt_olcdm"][2], pub["m16_zt_olcdm"][1], ps.INK_SECONDARY, "published"),
        ("  redone here, same data and prior", *prof(m16["olcdm_profile"]), ps.INK_SECONDARY, "ours"),
    ]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.4),
                                  gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.42})
    zt1, s1 = d1["z_transition"]
    zt3 = float(cosmo.z_transition_flat(d3["Om"][0]))
    _forest(ax, rows, ref_bands=[(zt1 - s1, zt1 + s1, ps.C_SN, "D1 supernovae")],
            xlabel="redshift at which the expansion began to speed up, z_t", xlim=(0, 1.3))
    ax.axvline(zt3, color=ps.C_CMB, lw=1.1)
    ax.text(zt3 + 0.01, -0.6, f"D3 ΛCDM: {zt3:.2f}", color=ps.C_CMB, fontsize=7.5, va="bottom")
    ax.set_title("When did the acceleration begin?")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", color=ps.INK_SECONDARY,
                              label="this project (68%)"),
                       Line2D([], [], marker="o", ls="none", mfc=ps.SURFACE,
                              mec=ps.INK_SECONDARY, label="published")],
              loc="lower right", fontsize=7.5)

    zs = art["ZS"]
    g = art["gp"]["SE kernel, full covariance"]
    ax2.axhspan(-1.2, 0, color=ps.EXCLUDED, lw=0, zorder=0)
    ax2.text(2.35, -1.05, "accelerating", ha="right", fontsize=8, color=ps.MUTED)
    ax2.text(2.35, 0.9, "decelerating", ha="right", fontsize=8, color=ps.MUTED)
    _band(ax2, zs, g["bands"]["q"])
    zz = np.linspace(0, cfg.Z_MAX_PLOT, 300)
    ax2.plot(zz, cosmo.deceleration(zz, d3["Om"][0]), color=ps.C_CMB, lw=1.3,
             label="D3 ΛCDM (D1's is almost identical)")
    ax2.plot(zz, cosmo.deceleration(zz, d1["Om_calibrated"][0]), color=ps.C_SN, lw=1.0,
             ls=(0, (4, 2)), label="D1 ΛCDM")
    ax2.axhline(0.5, color=ps.C_OTHER_MODEL, ls=(0, (5, 2)), lw=1.1, label="matter only (q = ½)")
    ax2.axhline(0, color=ps.AXIS, lw=0.8)
    ax2.set_xlim(0, cfg.Z_MAX_PLOT)
    ax2.set_ylim(-1.2, 1.2)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("deceleration parameter q(z)")
    ax2.set_title("Slowing or speeding up, from the chronometers alone")
    h = [Patch(color=BLUE_FILL_68, label="chronometers alone (68%, 95%)")] + \
        ax2.get_legend_handles_labels()[0]
    ax2.legend(handles=h, loc="upper left", fontsize=7.5)
    q0 = r["gp"]["SE kernel, full covariance"]["q0"]
    pacc = r["gp"]["SE kernel, full covariance"]["p_accelerating_today"]
    ps.caption(fig, f"Left: estimates of the transition redshift (median and 68% interval; profile-χ² "
               f"intervals for the two-line fits). Right: q = (1+z) H'/H − 1 from the Gaussian-process "
               f"curves; today q0 = {q0['median']:.2f} +{q0['plus']:.2f} −{q0['minus']:.2f}, "
               f"and {100 * pacc:.0f}% of the curves are accelerating today. The derivative of noisy "
               "data is poorly determined, so the band widens quickly where data are sparse.")
    return _save(fig, "03_transition")


# ---------------------------------------------------------------------------
# 04: H0 and the cost of ignoring correlations
# ---------------------------------------------------------------------------


def _hpd_levels(p, probs=(0.954, 0.683)):
    flat = np.sort(p.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    return [flat[np.searchsorted(cum, q)] for q in probs]


def fig_hubble_constant(r, art):
    d1, d3 = r["inputs"]["D1"], r["inputs"]["D3"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0),
                                  gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.55})
    for name, style in (("full", "fill"), ("diagonal only", "line")):
        h, om, chi2 = art["planes"][name]
        p = np.exp(-0.5 * (chi2 - chi2.min()))
        levels = _hpd_levels(p)
        if style == "fill":
            ax.contourf(h, om, p, levels=[levels[0], levels[1], p.max() * 1.01],
                        colors=[BLUE_FILL_68, ps.BLUE], alpha=0.9)
        else:
            ax.contour(h, om, p, levels=levels, colors=[ps.INK_SECONDARY], linewidths=1.0,
                       linestyles=["-", "-"])
    lrr = cfg.PUBLISHED
    ax.errorbar([lrr["lrr_H0"][0]], [lrr["lrr_Om"][0]], xerr=[[5.4], [5.4]],
                yerr=[[lrr["lrr_Om"][2]], [lrr["lrr_Om"][1]]], fmt="o", mfc=ps.SURFACE,
                mec=ps.INK, ecolor=ps.INK, elinewidth=0.8, ms=5, zorder=6)
    ax.errorbar([d1["H0"][0]], [d1["Om_calibrated"][0]], xerr=[[d1["H0"][1]], [d1["H0"][1]]],
                yerr=[[d1["Om_calibrated"][1]], [d1["Om_calibrated"][1]]], fmt="s", color=ps.C_SN,
                ms=5, zorder=6)
    ax.errorbar([d3["H0"][0]], [d3["Om"][0]], xerr=[[d3["H0"][1]], [d3["H0"][1]]],
                yerr=[[d3["Om"][1]], [d3["Om"][1]]], fmt="D", color=ps.C_CMB, ms=5, zorder=6)
    ax.set_xlim(48, 90)
    ax.set_ylim(0.08, 0.75)
    ax.set_xlabel("Hubble constant H0 (km/s/Mpc)")
    ax.set_ylabel("matter density Ωm")
    ax.set_title("Flat ΛCDM from the chronometers")
    ax.legend(handles=[Patch(color=ps.BLUE, label="full covariance (68%, 95%)"),
                       Line2D([], [], color=ps.INK_SECONDARY, label="diagonal errors only"),
                       Line2D([], [], marker="o", ls="none", mfc=ps.SURFACE, mec=ps.INK,
                              label="published (Moresco+2022)"),
                       Line2D([], [], marker="s", ls="none", color=ps.C_SN, label="D1 supernovae"),
                       Line2D([], [], marker="D", ls="none", color=ps.C_CMB, label="D3 Planck CMB")],
              loc="upper right", fontsize=7.5)

    fl = r["flat_lcdm"]
    gpr = r["gp"]
    def row(lab, s, col=ps.C_OURS, kind="ours"):
        return (lab, s["median"], s["minus"], s["plus"], col, kind)
    rows = [
        row("Flat ΛCDM, full covariance", fl["full"]["H0"]),
        row("  + stellar-library systematic", fl["full + stellar library"]["H0"]),
        row("  systematics on Moresco's 15 only", fl["systematics on Moresco's 15 only"]["H0"]),
        row("  diagonal errors only", fl["diagonal only"]["H0"], ps.INK_SECONDARY),
        row("  published recipe (scaled by the data)", fl["full, published recipe"]["H0"],
            ps.INK_SECONDARY),
        row("Smooth reconstruction (GP), extrapolated", gpr["SE kernel, full covariance"]["H0"]),
        row("  hyperparameters marginalised", gpr["SE kernel, hyperparameters marginalised"]["H0"]),
        row("  Matérn kernel", gpr["Matern-9/2 kernel, full covariance"]["H0"]),
        row("  constant mean marginalised", gpr["SE kernel, constant mean marginalised"]["H0"]),
        row("  straight-line mean marginalised",
            gpr["SE kernel, straight-line mean marginalised"]["H0"]),
        ("Moresco+2022 (published)", 66.5, 5.4, 5.4, ps.INK_SECONDARY, "published"),
        ("Yu, Ratra & Wang 2018, GP (published)", 67.0, 4.0, 4.0, ps.INK_SECONDARY, "published"),
    ]
    h1, s1 = d1["H0"]
    h3, s3 = d3["H0"]
    _forest(ax2, rows, ref_bands=[(h3 - s3, h3 + s3, ps.C_CMB, "D3: CMB"),
                                  (h1 - s1, h1 + s1, ps.C_SN, "D1: Cepheids")],
            xlabel="H0 (km/s/Mpc)", xlim=(54, 84))
    ax2.set_title("H0 from the chronometers: too uncertain to settle the tension", loc="right")
    ps.caption(fig, "Left: ignoring the correlated stellar-population systematic (grey) shrinks the "
               f"H0 error from ±{0.5 * (fl['full']['H0']['plus'] + fl['full']['H0']['minus']):.1f} to "
               f"±{0.5 * (fl['diagonal only']['H0']['plus'] + fl['diagonal only']['H0']['minus']):.1f} "
               "km/s/Mpc: precision the data do not have. Right: medians and 68% intervals; the "
               "chronometers sit between the CMB's 66.7 and the Cepheids' 73.6, consistent with both. "
               "The reconstruction rows depend on the GP's prior mean, and fed an exact ΛCDM curve the "
               "zero-mean GP returns an H0 ~3 km/s/Mpc low (shrinkage); the flat-ΛCDM fit is unbiased.")
    return _save(fig, "04_hubble_constant")


# ---------------------------------------------------------------------------
# 05: where the error bars come from
# ---------------------------------------------------------------------------


def fig_error_budget(r, art):
    cc, cov, syst = art["cc"], art["covs"]["full"], art["syst"]
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.9), gridspec_kw={"wspace": 0.3})
    im = ax.imshow(corr, cmap="Blues", vmin=0, vmax=1, origin="lower")
    ticks = np.arange(0, cc["z"].size, 4)
    ax.set_xticks(ticks, [f"{cc['z'][i]:.2f}" for i in ticks], rotation=90, fontsize=7)
    ax.set_yticks(ticks, [f"{cc['z'][i]:.2f}" for i in ticks], fontsize=7)
    ax.grid(False)
    ax.set_xlabel("measurement (sorted by redshift)")
    ax.set_ylabel("measurement (sorted by redshift)")
    ax.set_title("Correlation between the 32 measurements")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("correlation coefficient", color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)

    zz = np.linspace(0.05, 2.0, 300)
    tot = np.sqrt(sum(np.interp(zz, syst["z"], syst[t]) ** 2 for t in cfg.SYSTEMATIC_TERMS))
    ax2.plot(zz, tot, color=ps.C_OURS, lw=1.8,
             label="shared: stellar-population model + IMF (used)")
    lib = np.sqrt(tot**2 + np.interp(zz, syst["z"], syst["stlib"]) ** 2)
    ax2.plot(zz, lib, color=ps.C_OURS, lw=1.2, ls=(0, (4, 2)),
             label="… + stellar library (left out; a variant)")
    for m, mk in ps.METHOD_MARKERS.items():
        sel = cc["method"] == m
        face = ps.C_DATA if m == "D" else ps.SURFACE
        ax2.plot(cc["z"][sel], 100 * cc["err"][sel] / cc["H"][sel], mk, ms=5, mfc=face,
                 mec=ps.C_DATA, mew=0.9, ls="none",
                 label=f"each point's own error: {METHOD_LABEL[m]}")
    ax2.set_xlim(0, 2.05)
    ax2.set_ylim(0, 75)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("error on H(z) (%)")
    ax2.set_title("Each point's own error vs the error they share")
    ax2.legend(loc="upper left", fontsize=7.5)
    ps.caption(fig, "The shared term (Moresco et al. 2020) is the change in H(z) if a different "
               "stellar-population-synthesis model were right: it would shift every measurement the "
               "same way, so it is fully correlated. It is smaller than most points' own errors, "
               "but it does not average down when points are combined. Mean off-diagonal correlation: "
               f"{r['data']['mean_offdiag_correlation']:.2f}.")
    return _save(fig, "05_error_budget")


# ---------------------------------------------------------------------------
# 06: inside a chronometer
# ---------------------------------------------------------------------------


def fig_inside_a_chronometer(r, art):
    cl = art["clock"]
    m12 = cl["m12"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.9), gridspec_kw={"wspace": 0.36,
                                                                     "width_ratios": [1.2, 1, 1]})
    ax, ax2, ax3 = axes
    marks = {"low mass": ("o", "log M = 11.0–11.25"), "high mass": ("s", "log M = 11.25–11.5")}
    for mass, (mk, lab) in marks.items():
        rows = np.array(m12["medians"][mass])
        ax.errorbar(rows[:, 0], rows[:, 1], yerr=rows[:, 2], fmt=mk + "-", color=ps.C_OURS,
                    ms=5.5, lw=1.0, elinewidth=1, label=f"Moresco+2012 selection, {lab}")
    for (mass, rows), mk in zip(cl["simple"].items(), ("o", "s")):
        rows = np.array(rows)
        ax.errorbar(rows[:, 0], rows[:, 1], yerr=rows[:, 2], fmt=mk + "--", color=ps.INK_SECONDARY,
                    mfc=ps.SURFACE, ms=5.5, lw=1.0, elinewidth=1,
                    label=f"J3's SDSS galaxies, simple cuts, log M = {mass.replace('-', '–')}")
    # what pure ageing would do between z = 0.16 and 0.22 (D3 LCDM, A = 0.066 / Gyr)
    d3 = r["inputs"]["D3"]
    zz = np.linspace(0.155, 0.225, 50)
    t = cosmo.lookback_time(zz, d3["H0"][0], d3["Om"][0])
    A = float(np.polyval(np.polyfit([0.4, 1, 2.5], [0.02893, 0.060, 0.193], 2), 1.1))
    for d0 in (1.961, 1.992):
        ax.plot(zz, d0 - A * (t - t[0]) + 0.0 * zz, color=ps.C_CMB, lw=1.0, ls=(0, (1, 1.5)))
    ax.text(0.226, 1.992 - A * (t[-1] - t[0]), "pure ageing\n(D3 ΛCDM, A = 0.066/Gyr)",
            color=ps.C_CMB, fontsize=7, va="center")
    ax.set_xlim(0.15, 0.255)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("4000 Å break, D4000ₙ (median)")
    ax.set_title("The clock: the 4000 Å break grows with age")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.14), fontsize=7.2)

    for model, col, ls_ in (("BC03", ps.C_OURS, "-"), ("MaStro", ps.C_SN, "-")):
        rows = np.array(m12["A_rows"][model])
        Z = np.linspace(0.4, 2.5, 200)
        rel = {k: np.polyfit(rows[:, 0], rows[:, 1] + s * rows[:, 2], 2)
               for k, s in (("mid", 0), ("lo", -1), ("hi", 1))}
        ax2.fill_between(Z, np.polyval(rel["lo"], Z), np.polyval(rel["hi"], Z), color=col,
                         alpha=0.2, lw=0)
        ax2.plot(Z, np.polyval(rel["mid"], Z), color=col, lw=1.5, ls=ls_,
                 label=f"{model} models" if model == "BC03" else "Maraston & Strömbäck 2011 models")
        ax2.plot(rows[:, 0], rows[:, 1], "o", color=col, ms=4.5)
    ax2.axvspan(1.066, 1.21, color=ps.MUTED, alpha=0.18, lw=0)
    ax2.text(1.14, 0.205, "measured\nmetallicity of\nthe SDSS bins", ha="center", fontsize=7,
             color=ps.INK_SECONDARY)
    ax2.set_xlabel("stellar metallicity Z / Z☉")
    ax2.set_ylabel("A: growth of D4000ₙ per Gyr")
    ax2.set_title("The calibration: how fast it ticks")
    ax2.set_ylim(0, 0.24)
    ax2.legend(loc="lower right", fontsize=7.5)

    pub = cfg.PUBLISHED["m12_sdss"]
    zpub = np.array(list(pub))
    ax3.errorbar(zpub - 0.002, [v[0] for v in pub.values()], yerr=[v[1] for v in pub.values()],
                 fmt="o", mfc=ps.SURFACE, mec=ps.INK, ecolor=ps.INK, ms=6,
                 label="Moresco+2012, published (BC03)")
    comb = r["clock"]["moresco2012_reread"]["BC03"]
    ax3.errorbar([c["z"] + 0.002 for c in comb], [c["H"] for c in comb],
                 yerr=[c["stat"] for c in comb], fmt="o", color=ps.C_OURS, ms=6,
                 label="re-read here from their medians")
    sh = r["clock"]["simple_cuts_H"].get("BC03", {})
    for (mass, pts), mk in zip(sh.items(), ("o", "s")):
        ax3.errorbar([p["z"] for p in pts], [p["H"] for p in pts], yerr=[p["stat"] for p in pts],
                     fmt=mk, mfc=ps.SURFACE, mec=ps.INK_SECONDARY, ecolor=ps.INK_SECONDARY,
                     ms=5.5, label=f"simple cuts, log M = {mass.replace('-', '–')}")
    zz = np.linspace(0.15, 0.23, 20)
    ax3.plot(zz, cosmo.hubble(zz, d3["H0"][0], d3["Om"][0]), color=ps.C_CMB, lw=1.2,
             label="D3 ΛCDM")
    ax3.set_xlim(0.165, 0.215)
    ax3.set_ylim(0, 185)
    ax3.set_xlabel("redshift z")
    ax3.set_ylabel("H(z) (km/s/Mpc)")
    ax3.set_title("The reading: H from the slope")
    ax3.legend(loc="lower left", fontsize=7)
    ps.caption(fig, "H(z) = −A/(1+z) · dz/dD4000ₙ. Left: median D4000ₙ in redshift bins. Blue: "
               "Moresco et al.'s (2012) published medians for their SDSS galaxies (SED-template "
               "selection with 2MASS photometry). Grey: the same SDSS spectra (J3's MPA-JHU table) "
               "with simple passive cuts (no Hα emission, u−r > 2.6, low specific SFR). The "
               "expected change is only ~0.03 between bins two apart, and the simple selection "
               "distorts it, giving H anywhere from ~20 to ~150 km/s/Mpc: picking the right "
               "galaxies is most of the method. Right: error bars are statistical only.", y=-0.26)
    return _save(fig, "06_inside_a_chronometer")


# ---------------------------------------------------------------------------
# 07: the future
# ---------------------------------------------------------------------------

FUTURE_STYLE = {
    "LCDM (D3's CMB fit)": (ps.C_CMB, "-", "ΛCDM (D3): exponential growth forever"),
    "Evolving dark energy (DESI DR2+CMB+Pantheon+, extrapolated)":
        (ps.C_DESI, "-", "DESI DR2 w0–wa, extrapolated (speculative)"),
    "Phantom, w = -1.24 (D1's 95% limit)": (ps.C_SN, "-", "w = −1.24 (D1's 95% limit)"),
    "Phantom, w = -1.5 (Caldwell et al. 2003)": (ps.C_PHANTOM, "-", "w = −1.5 (Caldwell et al. 2003)"),
    "Matter only (no dark energy)": (ps.C_OTHER_MODEL, (0, (5, 2)), "matter only (no dark energy)"),
}


def fig_future(r, art):
    fut, fart = r["future"], art["future"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 5.6),
                                  gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.55})
    ax.axvspan(0, 100, color=ps.EXCLUDED, lw=0, zorder=0)
    ax.text(50, 6e-3, "the future: a projection, not a measurement", ha="center",
            color=ps.MUTED, fontsize=8)
    for name, (col, ls_, lab) in FUTURE_STYLE.items():
        if name not in fart:
            continue
        t, a = fart[name]["t"], fart[name]["a"]
        ax.plot(t, a, color=col, ls=ls_, lw=1.7, label=lab, zorder=3)
        e = fut[name]
        if e.get("big_rip_from_now_gyr"):
            tr = e["big_rip_from_now_gyr"]
            ax.axvline(tr, color=col, lw=0.8, ls=(0, (1, 1.5)))
            ax.text(tr + 1.0, 2e2 if "1.5" in name else 2e3, f"Big Rip\n+{tr:.0f} Gyr", color=col,
                    fontsize=7.5)
        if e.get("acceleration_end"):
            te = e["acceleration_end"]["t_from_now_gyr"]
            ae = e["acceleration_end"]["a"]
            ax.plot([te], [ae], "o", color=col, ms=5, zorder=5)
            ps.annotate(ax, f"acceleration ends\n+{te:.0f} Gyr", (te, ae), (te + 9, ae / 9),
                        color=col)
    ax.axvline(0, color=ps.INK_SECONDARY, lw=0.9)
    ax.text(0.8, 1.5e-3, "today", ha="left", fontsize=8, color=ps.INK_SECONDARY)
    ax.set_yscale("log")
    ax.set_xlim(-15, 100)
    ax.set_ylim(1e-3, 1e4)
    ax.set_xlabel("time from today (billion years)")
    ax.set_ylabel("size of the universe relative to today, a(t)")
    ax.set_title("How the expansion could continue")
    ax.legend(loc="upper left", fontsize=7.5)

    tl = r["future"]["timeline"]
    labels = [k for k in tl["-1.5"]]
    y = np.arange(len(labels))[::-1]
    for w, col, off in (("-1.5", ps.C_PHANTOM, 0.17), (str(cfg.D1_W_95_LOW), ps.C_SN, -0.17)):
        vals = np.array([tl[w][k] for k in labels])
        ax2.barh(y + off, vals, height=0.32, color=col, left=0,
                 label=f"w = {w.replace('-', '−')}")
    def human(sec):
        for unit, s in (("Gyr", 3.15576e16), ("Myr", 3.15576e13), ("yr", 3.15576e7),
                        ("days", 86400), ("min", 60), ("s", 1)):
            if sec >= s:
                return f"{sec / s:.0f} {unit}" if sec / s >= 10 else f"{sec / s:.1f} {unit}"
        return f"{sec:.0e} s"
    for k, yy in zip(labels, y):
        ax2.text(tl["-1.5"][k] * 1.5, yy + 0.17, human(tl["-1.5"][k]), va="center", fontsize=7,
                 color=ps.INK_SECONDARY)
    ax2.set_xscale("log")
    ax2.set_xlim(1e-18, 1e19)
    ax2.set_yticks(y, [l.replace(" (electron's orbital period)", "") for l in labels], fontsize=8)
    ax2.tick_params(axis="y", length=0)
    ax2.grid(axis="y", visible=False)
    ax2.set_xlabel("torn apart this long before the Big Rip (seconds)")
    ax2.set_title("If w < −1: the order of destruction")
    ax2.legend(loc="lower right", fontsize=7.5)
    c03 = r["future"]["Caldwell 2003 check"]
    ps.caption(fig, "Left: scale factor for five universes, all starting from today's expansion. "
               "ΛCDM and the phantom models use D3's H0 and Ωm; the w0–wa curve uses DESI DR2's "
               "own best fit (with CMB and Pantheon+) and continues its w(a) = w0 + wa(1 − a) "
               "far beyond the range it was fitted on, which no data support. Right: Caldwell, "
               "Kamionkowski & Weinberg's (2003) rule, a system with orbital period P comes apart "
               "P·√(2|1+3w|)/(6π|1+w|) before the end; their check case gives "
               f"{c03['exact_gyr']:.1f} Gyr to the Rip by direct integration (their approximation: "
               f"{c03['approx_formula_gyr']:.1f}; quoted 22). For atoms the same rule with the electron's orbital "
               "period gives ~4e-17 s; Caldwell et al. quote ~1e-19 s.")
    return _save(fig, "07_future")


# ---------------------------------------------------------------------------
# 08: chronometers published since the compilation
# ---------------------------------------------------------------------------


def fig_since_2022(r, art):
    cc, cov = art["cc"], art["covs"]["full"]
    zs = art["ZS"]
    g = art["gp"]["SE kernel, full covariance"]
    d3 = r["inputs"]["D3"]
    sd = np.sqrt(np.diag(cov))
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    _band(ax, zs, g["bands"]["H"])
    zz = np.linspace(0, 1.45, 200)
    ax.plot(zz, cosmo.hubble(zz, d3["H0"][0], d3["Om"][0]), color=ps.C_CMB, lw=1.2,
            label="D3 ΛCDM")
    for m, mk in ps.METHOD_MARKERS.items():
        sel = cc["method"] == m
        ax.errorbar(cc["z"][sel], cc["H"][sel], yerr=sd[sel], fmt=mk, ms=3.5, mfc=ps.SURFACE,
                    mec=ps.MUTED, ecolor=ps.MUTED, elinewidth=0.7, zorder=3)
    new = sorted(r["since_2022"], key=lambda p: p["z"])
    key = []
    for k, p in enumerate(new, start=1):
        ax.errorbar([p["z"]], [p["H"]], yerr=[[p["err_minus"]], [p["err_plus"]]], fmt="D", ms=6,
                    color=ps.INK, elinewidth=1.3, zorder=6)
        ax.text(p["z"] + 0.008, p["H"] + 3, str(k), fontsize=8, color=ps.INK, weight="bold",
                zorder=7)
        key.append(f"{k}  z = {p['z']:.3g}: {p['label']}  ({p['pull_vs_gp']:+.1f}σ)")
    half = (len(key) + 1) // 2
    for col, chunk in ((0.0, key[:half]), (0.52, key[half:])):
        ax.text(col, -0.13, "\n".join(chunk), transform=ax.transAxes, ha="left", va="top",
                fontsize=7.8, color=ps.INK_SECONDARY, linespacing=1.5)
    ax.set_xlim(0, 1.45)
    ax.set_ylim(30, 230)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("H(z) (km/s/Mpc)")
    ax.set_title("Chronometers published since the 2022 compilation, against our reconstruction")
    ax.legend(handles=[Patch(color=BLUE_FILL_68, label="reconstruction from the 32 (68%, 95%)"),
                       Line2D([], [], marker="o", ls="none", mfc=ps.SURFACE, mec=ps.MUTED,
                              label="the 32 compiled measurements"),
                       Line2D([], [], marker="D", ls="none", color=ps.INK,
                              label="new measurements (2022–2026)"),
                       Line2D([], [], color=ps.C_CMB, label="D3 ΛCDM")], loc="upper left",
              fontsize=7.5)
    pulls = [p["pull_vs_gp"] for p in new]
    ps.caption(fig, f"Nine measurements from 2022–2026 papers (abstract values; statistical and "
               "systematic errors added in quadrature). Their distances from the reconstruction "
               f"range from {min(pulls):+.1f}σ to {max(pulls):+.1f}σ. They are compared, not "
               "fitted: the three DESI DR1 analyses share galaxies, Jiao et al. 2023 re-uses the "
               "LEGA-C galaxies behind the compilation's z = 0.75 point, and their shared "
               "model systematics are not published in a common form.", y=-0.2)
    return _save(fig, "08_since_2022")


# ---------------------------------------------------------------------------
# 09: validation against the published analyses
# ---------------------------------------------------------------------------


def fig_validation(r, art):
    v = r["validation"]
    pub = cfg.PUBLISHED
    m16 = r["moresco2016_replication"]
    rows = []

    def add(label, ours, pub_val, pub_err, nd=None):
        if nd is None:
            nd = 3 if abs(pub_val) < 2 else 2
        rows.append((label, ours, pub_val, pub_err, nd))

    add("Primer, 15 points, H0 (full)", v["primer_full"]["H0"]["median"], pub["primer_full_H0"][0],
        0.5 * (pub["primer_full_H0"][1] + pub["primer_full_H0"][2]))
    add("Primer, 15 points, H0 (diagonal)", v["primer_diag"]["H0"]["median"],
        pub["primer_diag_H0"][0], 0.5 * (pub["primer_diag_H0"][1] + pub["primer_diag_H0"][2]))
    add("Primer, 15 points, Ωm (full)", v["primer_full"]["Om"]["median"], pub["primer_full_Om"][0],
        0.5 * (pub["primer_full_Om"][1] + pub["primer_full_Om"][2]))
    add("Review, 32 points, H0", v["lrr_32"]["H0"]["median"], pub["lrr_H0"][0], pub["lrr_H0"][1])
    add("Review, 32 points, Ωm", v["lrr_32"]["Om"]["median"], pub["lrr_Om"][0],
        0.5 * (pub["lrr_Om"][1] + pub["lrr_Om"][2]))
    add("Moresco+2016, z_t, ΛCDM + H0 prior", m16["olcdm_profile"]["best"], pub["m16_zt_olcdm"][0],
        0.5 * (pub["m16_zt_olcdm"][1] + pub["m16_zt_olcdm"][2]))
    add("Moresco+2016, z_t, two lines + H0 prior", m16["piecewise_profile"]["best"],
        pub["m16_zt_piecewise"][0], pub["m16_zt_piecewise"][1])
    for c, (zp, (hp, sp)) in zip(r["clock"]["moresco2012_reread"]["BC03"], pub["m12_sdss"].items()):
        add(f"Moresco+2012, H({zp:.3f}), BC03", c["H"], hp, sp)
    for c, (zp, hp) in zip(r["clock"]["moresco2012_reread"]["MaStro"],
                           r["clock"]["moresco2012_published_M11"].items()):
        add(f"Moresco+2012, H({float(zp):.3f}), M11", c["H"], hp,
            dict(zip(["0.1791", "0.1993"], [4.78, 5.85]))[zp])
    c03 = r["future"]["Caldwell 2003 check"]
    add("Caldwell+2003, time to Big Rip (approx.)", c03["approx_formula_gyr"], 22.0, 0.5)

    fig, ax = plt.subplots(figsize=(10, 5.4))
    y = np.arange(len(rows))[::-1]
    ax.axvspan(-1, 1, color=ps.EXCLUDED, lw=0, zorder=0)
    for (lab, ours, p, e, nd), yy in zip(rows, y):
        pull = (ours - p) / e
        ax.plot([0, pull], [yy, yy], color=ps.AXIS, lw=1, zorder=1)
        ax.plot([pull], [yy], "o", color=ps.C_OURS, ms=6, zorder=3)
        ax.text(1.65, yy, f"ours {ours:.{nd}f}   published {p:.{nd}f} ± {e:.{nd}f}",
                va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.axvline(0, color=ps.INK_SECONDARY, lw=0.9)
    ax.set_yticks(y, [r_[0] for r_ in rows], fontsize=8, color=ps.INK)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(-1.5, 3.9)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel("(ours − published) / published 1σ error")
    ax.set_title("Checks against the published analyses")
    ps.caption(fig, "Each row repeats a published calculation with the same data and priors. "
               "Shaded: within one published standard deviation. Primer = Moresco's GitLab "
               "notebook (emcee); review = Moresco et al. 2022, Living Reviews in Relativity. "
               "M11 errors are the GitLab table's totals.")
    return _save(fig, "09_validation")


# ---------------------------------------------------------------------------
# 10: the scaling bias
# ---------------------------------------------------------------------------


def fig_scaling_bias(r, art):
    sb, est = r["scaling_bias"], art["scaling_bias"]
    fl = r["flat_lcdm"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={"wspace": 0.25})
    lab = "standard terms"
    truth = sb[lab]["truth_H0"]
    bins = np.linspace(45, 90, 46)
    styles = {"measured values (published recipe)": (ps.INK_SECONDARY, "step"),
              "smooth reconstruction (ours)": (ps.C_OURS, "step"),
              "true curve (ideal)": (ps.C_CMB, "step")}
    for key, (col, ht) in styles.items():
        v = est[lab][key]
        ax.hist(v, bins=bins, histtype=ht, color=col, lw=1.6, label=f"{key}: {v.mean():.1f}")
    ax.axvline(truth, color=ps.INK, lw=1.0)
    ax.set_xlim(45, 90)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
    ax.text(truth + 0.6, ax.get_ylim()[1] * 0.6, f"truth\n{truth:.1f}", fontsize=8)
    ax.set_xlabel("H0 fitted to a mock data set (km/s/Mpc)")
    ax.set_ylabel(f"mock data sets (of {sb[lab]['n_mock']})")
    ax.set_title("Mock chronometers: which covariance recovers the truth?")
    ax.legend(loc="upper left", fontsize=7.5, title="shared systematic scaled by the …  (mean H0)",
              title_fontsize=7.5, alignment="left")

    h = art["H0_AXIS"]
    dens = art["h0_densities"]
    curves = [("full", ps.C_OURS, "-", "scaled by the smooth curve"),
              ("full + stellar library", ps.C_OURS, (0, (4, 2)), "… + stellar library"),
              ("full, published recipe", ps.INK_SECONDARY, "-", "scaled by the data (published recipe)"),
              ("full + stellar library, published recipe", ps.INK_SECONDARY, (0, (4, 2)),
               "… + stellar library")]
    for key, col, ls_, lab2 in curves:
        s = fl[key]["H0"]
        ax2.plot(h, dens[key] / dens[key].max(), color=col, ls=ls_, lw=1.6,
                 label=f"{lab2}: {s['median']:.1f} ± {0.5 * (s['plus'] + s['minus']):.1f}")
    ax2.set_xlim(40, 95)
    ax2.set_ylim(0, 1.4)
    ax2.set_xlabel("H0 (km/s/Mpc), flat ΛCDM, the 32 real measurements")
    ax2.set_ylabel("posterior (peak = 1)")
    ax2.set_title("The real data: adding an error should not move the answer")
    ax2.legend(loc="upper left", fontsize=7.5)
    s = sb[lab]
    ps.caption(fig, f"Left: {s['n_mock']} mock data sets drawn from a flat ΛCDM (D3's parameters) with "
               "the full chronometer covariance, each fitted with the shared systematic scaled three "
               "ways. Scaling by the measured values gives points that scattered high a larger "
               "error, so fits lean on the low ones: mean bias "
               f"{s['measured values (published recipe)']['bias']:+.1f} km/s/Mpc, against "
               f"{s['smooth reconstruction (ours)']['bias']:+.1f} for our smooth scaling and "
               f"{s['true curve (ideal)']['bias']:+.1f} for the ideal (the small residual comes from "
               "quoting posterior medians). This is Peelle's pertinent puzzle of nuclear-data fitting "
               "(d'Agostini 1994). Right: on the real data, the stellar-library term moves the "
               "published recipe's H0 down by 4 km/s/Mpc; with smooth scaling it only widens it. (The "
               "step at 50 is the edge of the flat prior used by Moresco's primer.)")
    return _save(fig, "10_scaling_bias")


def make_all(results, art) -> list:
    written = []
    for fn in (fig_hubble_rate, fig_speedometer, fig_transition, fig_hubble_constant,
               fig_error_budget, fig_inside_a_chronometer, fig_future, fig_since_2022,
               fig_validation, fig_scaling_bias):
        written += fn(results, art)
    return written
