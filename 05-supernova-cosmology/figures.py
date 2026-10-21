"""Every figure, from the results of ``analysis.run``.

    01_hubble_diagram     Hubble's law nearby, and the full diagram out to z = 2.26
    02_acceleration       residuals from an empty, coasting universe: the discovery
    03_omega_plane        the (Omega_m, Omega_Lambda) plane: the 1998 figure, redone
    04_hubble_constant    the distance ladder's last rung, and the Hubble tension
    05_dark_energy_w      is dark energy a cosmological constant?
    06_expansion_history  what the numbers mean: the size of the universe over time
    07_validation         our likelihood vs the Pantheon+ team's own chains

Each is written as PNG and SVG into ``figures/``.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import NullFormatter
from scipy.ndimage import gaussian_filter

import cosmology as cosmo
import likelihood as lk
from config import FIGURE_DIR, PUBLISHED_PANTHEON, SURVEY_GROUPS, Z_MIN_COSMOLOGY
from plotstyle import (
    AQUA, AXIS, BLUE_550, C_BINNED, C_DATA, C_OTHER_MODEL, C_OURS, C_PLANCK, C_PUBLISHED,
    C_SURVEYS, CONTOUR_FILLS, EXCLUDED, GRID, INK, INK_SECONDARY, MUTED, ORANGE, SURFACE,
    annotate, caption, save, use_style,
)

use_style()

OM = r"$\Omega_\mathrm{m}$"
OL = r"$\Omega_\Lambda$"
H0 = r"$H_0$"
C_KM_S = 299_792.458

DATA_CREDIT = ("Data: Pantheon+ (Scolnic et al. 2022; Brout et al. 2022) with SH0ES Cepheid "
               "distances (Riess et al. 2022), github.com/PantheonPlusSH0ES/DataRelease.")
ERRORBAR_NOTE = ("Error bars are the square root of the covariance diagonal, for display; "
                 "every fit uses the full 1701 × 1701 covariance.")

SURVEY_SHORT = ["nearby searches", "SDSS, PS1, SNLS, DES", "Hubble Space\nTelescope"]


def _num(x: float, spec: str) -> str:
    """Format a number with a true minus sign, matching the axis tick labels."""
    return format(x, spec).replace("-", "−")


def _pm(s: dict, spec: str = ".3f") -> str:
    return f"{_num(s['mean'], spec)} ± {_num(s['std'], spec)}"


def _screen_angle(ax, x, y, i: int) -> float:
    """Angle on screen (degrees) of the curve (x, y) at index i."""
    p = ax.transData.transform(np.column_stack([x[i - 1:i + 2:2], y[i - 1:i + 2:2]]))
    return float(np.degrees(np.arctan2(p[1, 1] - p[0, 1], p[1, 0] - p[0, 0])))


def _label_curve(ax, x, y, x_at, text, *, dy=5, color=INK_SECONDARY, size=8, **kw):
    """Write ``text`` along a curve, offset ``dy`` points above it (negative: below)."""
    i = int(np.clip(np.searchsorted(x, x_at), 1, len(x) - 2))
    ax.annotate(text, xy=(x[i], y[i]), xytext=(0, dy), textcoords="offset points",
                rotation=_screen_angle(ax, x, y, i), rotation_mode="anchor",
                ha=kw.pop("ha", "center"), va="bottom" if dy >= 0 else "top",
                fontsize=size, color=color, **kw)


def _hpd_fill(ax, x, y, p, *, zorder=2):
    """Filled 68/95/99.7% highest-density regions of a gridded posterior."""
    levels = lk.credible_levels(p)  # [68, 95, 99.7] thresholds, descending
    ax.contourf(x, y, p.T, levels=[levels[2], levels[1], levels[0], p.max() * 1.01],
                colors=CONTOUR_FILLS, zorder=zorder)
    ax.contour(x, y, p.T, levels=[levels[0]], colors=[BLUE_550], linewidths=0.8,
               zorder=zorder + 0.1)


def _hpd_lines(ax, x, y, p, color, *, masses=(0.682689, 0.9545), lw=1.6, zorder=3):
    levels = sorted(lk.credible_levels(p, masses))
    return ax.contour(x, y, p.T, levels=levels, colors=[color], linewidths=lw, zorder=zorder)


def _chain_density_2d(xs, ys, x_edges, y_edges, smooth=1.5):
    h, _, _ = np.histogram2d(xs, ys, bins=[x_edges, y_edges])
    h = gaussian_filter(h, smooth)
    return h / h.sum()


def _model_minus_empty(z, om, ol, w=-1.0):
    """5 log10(d_L(model) / d_L(empty universe)) on a smooth redshift grid."""
    dl = cosmo.luminosity_distance_shape(z, z, om, ol, w)[0]
    return 5.0 * np.log10(dl / cosmo.dl_milne(z))


def _no_y_axis(ax):
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


# ---------------------------------------------------------------------------
# 01 -- The Hubble diagram
# ---------------------------------------------------------------------------


def figure_hubble_diagram(results: dict, art: dict) -> list:
    hb, ds = art["hubble"], results["dataset"]
    h0, om = hb["H0"], hb["Om"]
    z, mu, err, grp, is_cal, used = (hb["z"], hb["mu"], hb["err"], hb["group"],
                                     hb["is_cal"], hb["used"])
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.5, 5.8),
                                     gridspec_kw={"width_ratios": [1, 1.45], "wspace": 0.22})

    # (a) Nearby: velocity against distance -- Hubble's law.
    zmax_a = 0.06
    dist = 10 ** ((mu - 25) / 5)  # Mpc, from the supernova's brightness
    near = (z < zmax_a) & ~is_cal
    for k in range(len(C_SURVEYS)):
        m = near & (grp == k)
        if m.any():
            ax_a.scatter(dist[m], C_KM_S * z[m], s=9, color=C_SURVEYS[k], alpha=0.75, lw=0,
                         zorder=3)
    cal_d = 10 ** ((art["calibrated"]["cal_mu"] - 25) / 5)
    cal_v = C_KM_S * art["calibrated"]["cal_z"]
    ax_a.scatter(cal_d, cal_v, s=26, marker="D", facecolor=SURFACE, edgecolor=INK, lw=0.9,
                 zorder=4)
    zz = np.linspace(1e-4, zmax_a, 300)
    dl = cosmo.luminosity_distance_shape(zz, zz, om, 1 - om)[0] * cosmo.hubble_distance_mpc(h0)
    ax_a.plot(dl, C_KM_S * zz, color=C_OURS, lw=2, zorder=5)
    ax_a.axhspan(0, C_KM_S * Z_MIN_COSMOLOGY, color=EXCLUDED, zorder=0, lw=0)
    ax_a.text(282, C_KM_S * Z_MIN_COSMOLOGY / 2, "z < 0.01: not used — a galaxy's own\n"
              "motion is too big a share of its velocity", fontsize=7.5, color=MUTED,
              va="center", ha="right")
    ax_a.set_xlim(0, 290)
    ax_a.set_ylim(0, C_KM_S * zmax_a)
    ax_a.set_xlabel("Distance, from the supernova's brightness (Mpc)")
    ax_a.set_ylabel("Recession velocity  cz  (km/s)")
    ax_a.set_title("Nearby: velocity grows in step with distance", pad=36)
    annotate(ax_a, f"slope = {H0} = {results['calibrated']['H0']['mean']:.1f} km/s/Mpc",
             xy=(185, C_KM_S * 0.0435), xytext=(25, 15200), color=INK)
    far_cal = int(np.argmax(cal_d))
    annotate(ax_a, "calibrators: distances from their\nCepheid stars, not their supernovae",
             xy=(cal_d[far_cal] + 2, cal_v[far_cal] - 150), xytext=(140, 4300),
             color=INK_SECONDARY)

    # (b) The whole diagram, log redshift.
    for k in range(len(C_SURVEYS)):
        m = used & (grp == k) & ~is_cal
        ax_b.errorbar(z[m], mu[m], yerr=err[m], fmt="none", ecolor=C_SURVEYS[k], alpha=0.25,
                      lw=0.6, zorder=2)
        ax_b.scatter(z[m], mu[m], s=10, color=C_SURVEYS[k], lw=0, alpha=0.85, zorder=3)
    excl = ~used & ~is_cal
    ax_b.scatter(z[excl], mu[excl], s=12, facecolor="none", edgecolor=MUTED, lw=0.7, zorder=3)
    ax_b.scatter(z[is_cal], mu[is_cal], s=22, marker="D", facecolor=SURFACE, edgecolor=INK,
                 lw=0.8, zorder=4)
    zz = np.geomspace(1e-3, 2.5, 400)
    ax_b.plot(zz, cosmo.distance_modulus_curve(zz, om, 1 - om, H0=h0), color=C_OURS, lw=2,
              zorder=5)
    ax_b.set_xscale("log")
    ax_b.set_xlim(8e-4, 2.6)
    ax_b.set_ylim(27, 47)
    ax_b.set_xticks([0.001, 0.01, 0.1, 1])
    ax_b.set_xticklabels(["0.001", "0.01", "0.1", "1"])
    ax_b.set_xlabel("Redshift  z   (the universe has grown 1 + z times since the light set out)")
    ax_b.set_ylabel("Distance modulus  μ = m − M  (mag)      farther →")
    ax_b.set_title(f"All {ds['n_light_curves']:,} light curves, z = 0.001 to {ds['z_max']:.2f}")

    # Lookback time along the top: when the light we see set out.
    lb_z = np.geomspace(1e-4, 3.0, 400)
    lb_t = cosmo.lookback_time_gyr(lb_z, om, 1 - om, H0=h0)
    sec = ax_b.secondary_xaxis("top", functions=(
        lambda zv: np.interp(zv, lb_z, lb_t), lambda tv: np.interp(tv, lb_t, lb_z)))
    sec.set_xticks([0.1, 1, 3, 5, 7, 9])
    sec.set_xticklabels(["0.1", "1", "3", "5", "7", "9"])
    sec.set_xlabel("…light emitted this many billion years ago", color=MUTED, fontsize=8)
    sec.tick_params(colors=MUTED, labelsize=7.5)

    for label, xy, ha in zip(SURVEY_SHORT, [(0.028, 34.0), (0.4, 39.4), (0.75, 44.9)],
                             ["left", "left", "right"]):
        ax_b.text(*xy, label, fontsize=8, color=INK_SECONDARY, ha=ha, va="center")

    names = list(SURVEY_GROUPS)
    handles = [Line2D([], [], ls="", marker="o", ms=5, color=C_SURVEYS[k],
                      label=f"{names[k]}  ({ds['survey_groups'][names[k]]})")
               for k in range(len(C_SURVEYS))]
    handles += [
        Line2D([], [], ls="", marker="D", ms=5, markerfacecolor=SURFACE, markeredgecolor=INK,
               label=f"in a Cepheid-calibrated galaxy ({ds['n_calibrator_light_curves']})"),
        Line2D([], [], ls="", marker="o", ms=5, markerfacecolor="none", markeredgecolor=MUTED,
               label=f"z < 0.01, not used ({ds['n_excluded_low_z']})"),
        Line2D([], [], color=C_OURS, lw=2,
               label=f"best fit: flat ΛCDM, {OM} = {om:.3f}, {H0} = {h0:.1f}"),
    ]
    ax_b.legend(handles=handles, loc="lower right", fontsize=7.5)
    caption(fig, f"{DATA_CREDIT} Distances use the supernova absolute magnitude "
                 f"M = {_num(hb['M'], '.3f')} calibrated on the Cepheid hosts. "
                 f"{ds['n_cosmology']:,} light curves at z > 0.01 enter the cosmology fits. "
                 f"{ERRORBAR_NOTE}")
    return save(fig, FIGURE_DIR / "01_hubble_diagram")


# ---------------------------------------------------------------------------
# 02 -- The discovery: residuals from an empty universe
# ---------------------------------------------------------------------------


def figure_acceleration(results: dict, art: dict) -> list:
    hb, bn = art["hubble"], art["binned"]
    comp = results["model_comparison"]
    zero = art["model_offsets"]["best_flat"]  # one shared zero point for every curve
    fig, ax = plt.subplots(figsize=(11.5, 6.6))

    used = hb["used"]
    ax.scatter(hb["z"][used], hb["resid_milne"][used], s=5, color=C_DATA, alpha=0.35, lw=0,
               zorder=1)
    ax.axhline(0, color=INK_SECONDARY, lw=1.0, zorder=2)

    ax.set_xscale("log")
    ax.set_xlim(Z_MIN_COSMOLOGY * 0.9, 2.5)
    ax.set_ylim(-0.8, 0.66)

    zz = np.geomspace(Z_MIN_COSMOLOGY, 2.4, 500)
    curves = {k: _model_minus_empty(zz, comp[k]["Om"], comp[k]["OL"]) + zero
              for k in ("eds", "open", "de_sitter", "best_flat")}
    for key in ("eds", "open", "de_sitter"):
        ax.plot(zz, curves[key], color=C_OTHER_MODEL, lw=1.3, zorder=3)
    ax.plot(zz, curves["best_flat"], color=C_OURS, lw=2.2, zorder=4)

    ax.errorbar(bn["z"], bn["value"], yerr=bn["err"], fmt="o", color=C_BINNED, ms=6.5,
                mec=SURFACE, mew=1.2, elinewidth=1.2, capsize=0, zorder=6)

    ax.set_xticks([0.01, 0.03, 0.1, 0.3, 1, 2])
    ax.set_xticklabels(["0.01", "0.03", "0.1", "0.3", "1", "2"])
    ax.set_xlabel("Redshift  z")
    ax.set_ylabel("Δμ = observed − empty universe  (mag)\n"
                  "← brighter, nearer          fainter, farther →")
    ax.set_title("Distant supernovae are too faint for a universe that is slowing down")

    dchi = {k: comp[k]["delta_chi2"] for k in comp}
    _label_curve(ax, zz, curves["de_sitter"], 0.5, f"dark energy only  (Δχ² = +{dchi['de_sitter']:.0f})")
    _label_curve(ax, zz, curves["open"], 1.0, f"open, matter only, {OM} = 0.3  (Δχ² = +{dchi['open']:.0f})",
                 dy=-4, ha="right")
    _label_curve(ax, zz, curves["eds"], 0.42, f"Einstein–de Sitter: matter only, {OM} = 1  "
                                              f"(Δχ² = +{dchi['eds']:.0f})", dy=-4)
    ax.text(0.0112, -0.022, "empty universe: coasting, neither speeding up nor slowing",
            fontsize=8, color=INK_SECONDARY, va="top")
    annotate(ax, "beyond z ≈ 1 they turn\nbrighter again: back\nthen the expansion\nwas still slowing",
             xy=(bn["z"][-1], bn["value"][-1] + bn["err"][-1] + 0.01), xytext=(1.22, 0.33),
             color=INK_SECONDARY)
    ex = results["lcdm_plane"]["exclusion"]
    ax.text(0.0112, 0.62,
            f"The best universe that is not accelerating (an empty one) fits worse by "
            f"Δχ² = {ex['no_acceleration']['delta_chi2']:.0f}:\nthe expansion is speeding up, "
            f"at {ex['no_acceleration']['sigma']:.1f}σ.",
            fontsize=8.5, color=INK, va="top", linespacing=1.4)

    best = comp["best_flat"]
    handles = [
        Line2D([], [], ls="", marker="o", ms=3, color=C_DATA, alpha=0.6,
               label=f"individual supernovae ({results['dataset']['n_cosmology']:,})"),
        Line2D([], [], ls="", marker="o", ms=6.5, color=C_BINNED, mec=SURFACE,
               label="binned (full covariance)"),
        Line2D([], [], color=C_OURS, lw=2.2,
               label=f"our best fit: flat ΛCDM, {OM} = {best['Om']:.3f}, {OL} = {best['OL']:.3f}"),
        Line2D([], [], color=C_OTHER_MODEL, lw=1.3, label="other universes (labelled on the curves)"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=7.5)
    caption(fig, f"{DATA_CREDIT} Zero is an empty universe, expanding at a constant rate. All "
                 "curves share one brightness zero point, pinned by the nearby supernovae where "
                 "every model predicts the same distances; they differ only in how distance grows "
                 "with redshift. Δχ² is relative to our best fit, with each model given its own "
                 "best zero point, so no model is handicapped by calibration. Binned points are "
                 "generalised least-squares averages with the full covariance, for display only.")
    return save(fig, FIGURE_DIR / "02_acceleration")


# ---------------------------------------------------------------------------
# 03 -- The (Omega_m, Omega_Lambda) plane
# ---------------------------------------------------------------------------


def figure_omega_plane(results: dict, art: dict) -> list:
    pl = art["plane"]
    om, ol, p = pl["Om"], pl["OL"], pl["post_physical"]
    res = results["lcdm_plane"]
    fig, ax = plt.subplots(figsize=(9.4, 7.6))
    x0, x1, y0, y1 = 0.0, 1.2, -0.35, 1.6
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

    # No Big Bang region (upper left).
    xb = np.linspace(1e-4, x1, 600)
    edge = cosmo.no_big_bang_boundary(xb)
    ax.fill_between(xb, edge, y1 + 1, color=EXCLUDED, zorder=0, lw=0)
    ax.plot(xb, edge, color=AXIS, lw=0.8, zorder=1)
    annotate(ax, "no Big Bang: traced back, such a universe\nshrinks to a minimum size and "
                 "re-expands", xy=(0.05, 1.42), xytext=(0.26, 1.47), color=MUTED)

    _hpd_fill(ax, om, ol, p)

    xs = np.linspace(x0, x1, 200)
    ax.plot(xs, 1 - xs, color=INK_SECONDARY, lw=1.0, zorder=3)  # flat
    ax.plot(xs, xs / 2, color=INK_SECONDARY, lw=1.0, zorder=3)  # q0 = 0
    ax.axhline(0, color=AXIS, lw=1.0, zorder=1)

    _label_curve(ax, xs, 1 - xs, 0.83, "flat", size=8.5)
    ax.text(0.63, 0.24, "open", fontsize=8, color=MUTED, style="italic")
    ax.text(1.02, 0.03, "closed", fontsize=8, color=MUTED, style="italic")
    _label_curve(ax, xs, xs / 2, 1.02, "accelerating ↑", size=8.5)
    _label_curve(ax, xs, xs / 2, 1.02, "decelerating ↓", size=8.5, dy=-5)
    ax.text(0.5, -0.025, f"no dark energy ({OL} = 0)", fontsize=7.5, color=MUTED, ha="center",
            va="top")

    for (x, y, label, dx, dy, ha) in [(0, 0, "empty", 0.025, -0.06, "left"),
                                      (1, 0, "Einstein–de Sitter", -0.03, -0.075, "right"),
                                      (0, 1, "dark energy only", 0.025, 0.015, "left")]:
        ax.plot(x, y, "o", ms=5, color=INK, mec=SURFACE, mew=1.2, zorder=6, clip_on=False)
        ax.text(x + dx, y + dy, label, fontsize=8, color=INK_SECONDARY, ha=ha)
    b = res["best"]
    ax.plot(b["Om"], b["OL"], "+", ms=11, mew=1.8, color=INK, zorder=7)
    pk = results["planck18"]["Om"][0]
    ax.plot(pk, 1 - pk, "o", ms=6.5, color=C_PLANCK, mec=SURFACE, mew=1.2, zorder=7)
    annotate(ax, f"our best fit ({b['Om']:.3f}, {b['OL']:.3f})", xy=(b["Om"] + 0.012, b["OL"]),
             xytext=(0.5, 0.64), color=INK)
    annotate(ax, "Planck 2018 (CMB, assuming flat)", xy=(pk + 0.012, 1 - pk + 0.004),
             xytext=(0.52, 0.8), color=INK)

    ex = res["exclusion"]
    ax.text(0.62, 1.33,
            f"{OM} = {_pm(res['Om'])}\n{OL} = {_pm(res['OL'])}\n"
            f"dark energy detected at {ex['no_lambda']['sigma']:.1f}σ   (1998: 3–4σ)\n"
            f"expansion accelerating at {ex['no_acceleration']['sigma']:.1f}σ",
            fontsize=8.5, color=INK, va="top", linespacing=1.55)
    handles = [Patch(color=c, label=l) for c, l in
               zip(CONTOUR_FILLS[::-1], ["68.3% credible region", "95.4%", "99.7%"])]
    ax.legend(handles=handles, loc="lower left", fontsize=7.5)
    ax.set_xlabel(f"{OM}   matter density, as a fraction of the critical density")
    ax.set_ylabel(f"{OL}   dark-energy density")
    ax.set_title("What the universe contains, from supernovae alone")
    caption(fig, f"{DATA_CREDIT} Curvature left free; the only prior is {OM} ≥ 0. "
                 f"{results['dataset']['n_cosmology']:,} supernovae at z > 0.01, full statistical "
                 "+ systematic covariance. On the flat line the two densities add up to 1. "
                 f"Above the accelerating line {OL} > {OM}/2: dark energy's push beats matter's "
                 "pull and the deceleration parameter q₀ is negative.")
    return save(fig, FIGURE_DIR / "03_omega_plane")


# ---------------------------------------------------------------------------
# 04 -- H0: the distance ladder and the tension
# ---------------------------------------------------------------------------


def figure_hubble_constant(results: dict, art: dict) -> list:
    cb, res = art["calibrated"], results["calibrated"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.5, 5.4),
                                     gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.55})

    # (a) The standard candle: absolute magnitude of each calibrated supernova.
    d = 10 ** ((cb["cal_mu"] - 25) / 5)
    M, s = res["M"]["mean"], res["calibrator_M_scatter"]
    ax_a.axhspan(M - s, M + s, color=GRID, alpha=0.55, lw=0, zorder=0)
    ax_a.errorbar(d, cb["cal_M"], yerr=cb["cal_err"], fmt="none", ecolor=MUTED,
                  elinewidth=0.9, alpha=0.8, zorder=2)
    ax_a.scatter(d, cb["cal_M"], s=24, color=C_BINNED, edgecolor=SURFACE, lw=0.8, zorder=3)
    ax_a.axhline(M, color=C_OURS, lw=2, zorder=2.5)
    ax_a.set_xscale("log")
    ax_a.set_xlim(5, 60)
    ax_a.set_xticks([5, 10, 20, 40])
    ax_a.set_xticklabels(["5", "10", "20", "40"])
    ax_a.xaxis.set_minor_formatter(NullFormatter())
    ax_a.set_ylim(M + 0.75, M - 0.75)  # magnitudes: brighter is up
    ax_a.set_xlabel("Distance to the host galaxy, from its Cepheid stars (Mpc)")
    ax_a.set_ylabel(r"Supernova absolute magnitude  $M = m - \mu_\mathrm{Cepheid}$"
                    "\n(standardised)          brighter →")
    ax_a.set_title("Calibrating the candle")
    ax_a.text(5.3, M - 0.66,
              f"M = {_num(M, '.3f')} ± {res['M']['std']:.3f}   (blue line)\n"
              f"one supernova scatters by ±{s:.2f} mag (grey band),\n"
              f"i.e. ±{100 * (10 ** (s / 5) - 1):.0f}% in distance",
              fontsize=8.5, color=INK, va="top", linespacing=1.4)
    ax_a.text(58, M + 0.68,
              f"{results['dataset']['n_calibrator_light_curves']} light curves of "
              f"{results['dataset']['n_calibrator_supernovae']} supernovae",
              fontsize=7.5, color=MUTED, ha="right", va="bottom")

    # (b) H0 values, forest plot.
    ours = res["H0"]
    kin = results["cross_check_sh0es_hf"]
    pub = PUBLISHED_PANTHEON["shoes_flat_lcdm_H0"]
    rows = [
        ("This project: full fit", ours["mean"], ours["std"], C_OURS, "o"),
        (f"This project: only the {kin['n_hubble_flow']}\nSH0ES supernovae, no cosmology",
         kin["H0"], kin["H0_err"], C_OURS, "s"),
        ("Pantheon+ & SH0ES\n(Brout et al. 2022)", pub[0], pub[1], C_PUBLISHED, "o"),
        ("SH0ES (Riess et al. 2022)", res["sh0es_r22_H0"][0], res["sh0es_r22_H0"][1],
         C_PUBLISHED, "o"),
        ("Planck 2018: early universe\n+ ΛCDM, predicted", res["planck_H0"][0],
         res["planck_H0"][1], C_PLANCK, "o"),
    ]
    pv, pe = res["planck_H0"]
    ax_b.axvspan(pv - pe, pv + pe, color=C_PLANCK, alpha=0.12, lw=0, zorder=0)
    ax_b.axvspan(ours["mean"] - ours["std"], ours["mean"] + ours["std"], color=C_OURS,
                 alpha=0.12, lw=0, zorder=0)
    ys = np.arange(len(rows))[::-1].astype(float)
    for y, (label, val, err, color, marker) in zip(ys, rows):
        ax_b.errorbar(val, y, xerr=err, fmt=marker, color=color, ms=7, mec=SURFACE, mew=1.2,
                      elinewidth=2, capsize=0, zorder=3)
        ax_b.text(val, y + 0.2, f"{val:.2f} ± {err:.2f}", fontsize=7.5, color=INK_SECONDARY,
                  ha="center", va="bottom")
    ax_b.set_yticks(ys)
    ax_b.set_yticklabels([r[0] for r in rows])
    ax_b.tick_params(axis="y", length=0, labelsize=8.5, labelcolor=INK)
    ax_b.spines["left"].set_visible(False)
    ax_b.grid(axis="y", visible=False)
    ax_b.set_ylim(-1.0, len(rows) - 0.35)
    ax_b.set_xlim(64.5, 77)
    ax_b.set_xlabel(f"{H0}  (km/s/Mpc)")
    ax_b.set_title("The Hubble tension")
    gap = ours["mean"] - pv
    ax_b.annotate("", xy=(ours["mean"], -0.68), xytext=(pv, -0.68),
                  arrowprops={"arrowstyle": "<->", "color": INK_SECONDARY, "lw": 1})
    ax_b.text((ours["mean"] + pv) / 2, -0.62,
              f"{gap:.1f} km/s/Mpc apart = {res['tension_with_planck_sigma']:.1f}σ",
              fontsize=8.5, color=INK, ha="center", va="bottom")

    caption(fig, f"{DATA_CREDIT} Left: each point is one light curve's standardised peak "
                 "brightness minus the Cepheid distance modulus of its host (from SH0ES). Right: our "
                 "fit calibrates M on those hosts and reads H0 off the supernovae in the smooth "
                 "expansion (z > 0.01). The Cepheid distances -- and the geometric anchors beneath "
                 "them -- are taken from SH0ES, not re-derived here. Planck's value is inferred from "
                 "the young universe assuming ΛCDM: a prediction for today, not a local measurement.")
    return save(fig, FIGURE_DIR / "04_hubble_constant")


# ---------------------------------------------------------------------------
# 05 -- Dark energy's equation of state
# ---------------------------------------------------------------------------


def figure_dark_energy_w(results: dict, art: dict) -> list:
    wc, res = art["wcdm"], results["flat_wcdm"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 5.4),
                                     gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.22})
    _hpd_fill(ax_a, wc["Om"], wc["w"], wc["post_unboxed"])
    ax_a.axhline(-1, color=INK_SECONDARY, lw=1.0, zorder=3)
    ax_a.text(0.69, -0.985, "cosmological constant, w = −1", fontsize=8, color=INK_SECONDARY,
              ha="right", va="bottom")
    b = res["best"]
    ax_a.plot(b["Om"], b["w"], "+", ms=11, mew=1.8, color=INK, zorder=7)
    pk = results["planck18"]["Om"][0]
    ax_a.axvline(pk, color=C_PLANCK, lw=1.2, zorder=2)
    ax_a.text(pk + 0.007, -1.82, f"Planck's {OM}", fontsize=7.5, color=INK_SECONDARY,
              rotation=90, va="bottom")
    ax_a.set_xlim(0, 0.7)
    ax_a.set_ylim(-1.85, -0.35)
    ax_a.set_xlabel(OM)
    ax_a.set_ylabel("w = pressure / (energy density × c²)")
    ax_a.set_title("Flat universe, dark energy with a constant w")
    ax_a.text(0.69, -0.42, "w > −1: dark energy thins out as space grows\n"
                           "w < −1: it grows denser (\"phantom\")",
              fontsize=7.5, color=MUTED, va="top", ha="right")
    handles = [Patch(color=c, label=l) for c, l in
               zip(CONTOUR_FILLS[::-1], ["68.3%", "95.4%", "99.7%"])]
    ax_a.legend(handles=handles, loc="lower right", fontsize=7.5)

    # (b) The marginal posterior for w.
    w = wc["w"]
    pw = wc["post"].sum(0)
    pw = pw / np.trapezoid(pw, w)
    ax_b.fill_between(w, pw, color=C_OURS, alpha=0.12, lw=0)
    ax_b.plot(w, pw, color=C_OURS, lw=2)
    ax_b.axvline(-1, color=INK_SECONDARY, lw=1.0)
    pub = PUBLISHED_PANTHEON["flat_wcdm_w"]
    ytop = pw.max() * 1.2
    ax_b.errorbar(pub[0], ytop, xerr=pub[1], fmt="o", color=C_PUBLISHED, ms=6, mec=SURFACE,
                  elinewidth=2, capsize=0)
    ax_b.text(pub[0], ytop * 1.04, f"Brout et al. 2022: {_num(pub[0], '.2f')} ± {pub[1]:.2f}",
              fontsize=7.5, color=INK_SECONDARY, ha="center", va="bottom")
    ax_b.text(-1.97, pw.max() * 0.72,
              f"ours: w = {_num(res['w']['mean'], '.2f')} ± {res['w']['std']:.2f}\n\n"
              f"a cosmological constant\n(w = −1) is only Δχ² = "
              f"{results['flat_wcdm']['delta_chi2_w_minus_1']:.1f}\nfrom the best fit:\n"
              "fully consistent", fontsize=8.5, color=INK, va="top", linespacing=1.4)
    ax_b.set_xlim(-2.0, -0.3)
    ax_b.set_ylim(0, ytop * 1.15)
    _no_y_axis(ax_b)
    ax_b.set_xlabel("w")
    ax_b.set_title("How likely each value of w is")
    caption(fig, f"{DATA_CREDIT} Supernovae alone. Left: credible regions with only {OM} ≥ 0 "
                 "assumed. Right: marginalised within the published analysis's prior "
                 f"(0.1 < {OM} < 0.9, −2 < w < 0.5). Along the {OM}–w degeneracy the answer depends "
                 "on where a prior cuts it, which is why supernova-only w is always quoted with its "
                 "prior; adding the CMB and baryon acoustic oscillations breaks the degeneracy.")
    return save(fig, FIGURE_DIR / "05_dark_energy_w")


# ---------------------------------------------------------------------------
# 06 -- The expansion history implied by the fit
# ---------------------------------------------------------------------------


def figure_expansion_history(results: dict, art: dict) -> list:
    res = results["calibrated"]
    h0, om = res["best"]["H0"], res["best"]["Om"]
    fig, ax = plt.subplots(figsize=(11.5, 6.4))
    t_end = 8.0
    ax.set_xlim(-15.2, t_end)
    ax.set_ylim(-0.2, 1.85)
    a = np.concatenate([np.geomspace(1e-5, 0.02, 60), np.linspace(0.02, 2.2, 1200)])
    models = [("ours", om, 1 - om, C_OURS, 2.2, "our best fit: speeds up"),
              ("eds", 1.0, 0.0, INK_SECONDARY, 1.4, "matter only: keeps slowing"),
              ("empty", 0.0, 0.0, MUTED, 1.4, "empty: constant rate")]
    ends = []
    for key, m, l, color, lw, end_label in models:
        t = cosmo.cosmic_time_of_a(a, m, l, H0=h0)
        t0 = cosmo.age_gyr(m, l, H0=h0)
        tr = t - t0
        past = a <= 1.0
        z_ord = 4 if key == "ours" else 3
        ax.plot(tr[past], a[past], color=color, lw=lw, zorder=z_ord)
        ax.plot(tr[~past], a[~past], color=color, lw=lw, ls=(0, (4, 3)), zorder=z_ord)
        ax.plot(-t0, 0, "o", ms=5.5, color=color, mec=SURFACE, mew=1.2, zorder=5)
        ends.append([float(np.interp(t_end, tr, a)), end_label, color])
        label = {"ours": f"Big Bang\n{t0:.1f} Gyr ago", "eds": f"{t0:.1f} Gyr ago",
                 "empty": f"{t0:.1f} Gyr\nago"}[key]
        ax.text(-t0 + (0.12 if key != "empty" else -0.12), -0.035, label,
                fontsize=7.5 if key != "ours" else 8, va="top",
                ha="left" if key != "empty" else "right",
                color=INK if key == "ours" else INK_SECONDARY)

    # Right-margin labels at the line ends, nudged apart with leader lines.
    ends.sort(key=lambda e: e[0], reverse=True)
    placed = []
    for y, label, color in ends:
        y_lab = y if not placed else min(y, placed[-1] - 0.085)
        placed.append(y_lab)
        ax.annotate(label, xy=(t_end, y), xytext=(t_end + 0.35, y_lab), va="center",
                    fontsize=8, color=INK if color == C_OURS else INK_SECONDARY,
                    annotation_clip=False,
                    arrowprops={"arrowstyle": "-", "color": AXIS, "lw": 0.7, "shrinkA": 0,
                                "shrinkB": 1})

    # When the supernova light set out: a rug along the bottom.
    z_sn = art["hubble"]["z"][art["hubble"]["used"]]
    lb = cosmo.lookback_time_gyr(z_sn, om, 1 - om, H0=h0)
    ax.vlines(-lb, -0.195, -0.135, color=INK, alpha=0.12, lw=0.6, zorder=2)
    ax.text(-lb.max() - 0.2, -0.165, f"when the light from our {z_sn.size:,} supernovae set out →",
            fontsize=7.5, color=INK_SECONDARY, ha="right", va="center")

    zt = res["transition"]["z"]
    tt = -res["transition"]["lookback_gyr"]
    ax.plot(tt, 1 / (1 + zt), "o", ms=7, color=C_OURS, mec=SURFACE, mew=1.4, zorder=6)
    annotate(ax, f"acceleration begins: z = {zt:.2f},\n{-tt:.1f} billion years ago",
             xy=(tt, 1 / (1 + zt)), xytext=(-11.8, 0.95), color=INK)
    ax.plot(0, 1, "o", ms=7, color=INK, mec=SURFACE, mew=1.4, zorder=6)
    ax.text(0.25, 0.97, "today", fontsize=8.5, color=INK, va="top")
    ax.axvline(-4.57, color=GRID, lw=1.0, zorder=0)
    ax.text(-4.47, 1.8, "Sun and Earth\nform", fontsize=7.5, color=MUTED, va="top")
    ax.text(0.25, 1.8, "future (projection) →", fontsize=7.5, color=MUTED, va="top")

    ax.legend(handles=[
        Line2D([], [], color=C_OURS, lw=2.2, label=f"our best fit: flat ΛCDM, {OM} = {om:.2f}"),
        Line2D([], [], color=INK_SECONDARY, lw=1.4, label="matter only (Einstein–de Sitter)"),
        Line2D([], [], color=MUTED, lw=1.4, label="empty (coasting)"),
    ], loc="upper left", fontsize=7.5, title=f"all expanding at {H0} = {h0:.1f} km/s/Mpc today",
        title_fontsize=7.5)
    ax.set_xlabel("Time from today (billions of years)")
    ax.set_ylabel("Size of the universe relative to today  (scale factor a)")
    ax.set_title("The expansion history our fit implies")
    caption(fig, f"{DATA_CREDIT} Curves are models, not data: flat ΛCDM with our best fit and two "
                 "textbook universes with the same expansion rate today. They agree now and differ "
                 "in their pasts, hence their ages. The faster measured H0 gives an age of "
                 f"{res['age_gyr']['mean']:.1f} ± {res['age_gyr']['std']:.1f} Gyr -- younger than "
                 "Planck's 13.8 Gyr: the Hubble tension seen from another side. Radiation is "
                 "neglected; it matters only in the first ~100,000 years.")
    return save(fig, FIGURE_DIR / "06_expansion_history")


# ---------------------------------------------------------------------------
# 07 -- Validation against the published chains
# ---------------------------------------------------------------------------


def figure_validation(results: dict, art: dict) -> list:
    ch, val = art["chains"], results["validation"]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.2),
                             gridspec_kw={"hspace": 0.45, "wspace": 0.28})

    def _one_d(ax, x, dens, samples, label, key, spec):
        dens = dens / np.trapezoid(dens, x)
        ax.hist(samples, bins=60, range=(x.min(), x.max()), density=True, histtype="step",
                color=C_PUBLISHED, lw=1.5)
        ax.plot(x, dens, color=C_OURS, lw=2)
        ax.set_ylim(0, dens.max() * 1.55)
        v = val[key]
        ax.set_title(label)
        _no_y_axis(ax)
        ax.text(0.02, 0.97,
                f"ours   {_pm(v['ours'], spec)}\n"
                f"chain  {_pm(v['chain'], spec)}\n"
                f"shift {_num(v['diff_in_sigma'], '+.2f')}σ, width ×{v['width_ratio']:.2f}",
                transform=ax.transAxes, fontsize=7.5, color=INK, va="top",
                family="DejaVu Sans Mono")

    f = art["flat"]
    m = (f["Om"] > 0.25) & (f["Om"] < 0.42)
    _one_d(axes[0, 0], f["Om"][m], f["post"][m], ch["flat_lcdm"]["Om"],
           f"Flat ΛCDM: {OM}", "flat_lcdm:Om", ".3f")
    axes[0, 0].set_xlabel(OM)

    ax = axes[0, 1]
    pl = art["plane"]
    _hpd_lines(ax, pl["Om"], pl["OL"], pl["post"], C_OURS, lw=2)
    xe, ye = np.linspace(0.0, 1.2, 121), np.linspace(-0.4, 1.6, 201)
    dens = _chain_density_2d(ch["lcdm"]["Om"], ch["lcdm"]["OL"], xe, ye)
    _hpd_lines(ax, 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), dens, C_PUBLISHED, lw=1.5)
    ax.set_xlim(0.05, 0.6)
    ax.set_ylim(0.25, 1.0)
    ax.set_xlabel(OM)
    ax.set_ylabel(OL)
    ax.set_title("Curvature free: 68% and 95% regions")

    ax = axes[0, 2]
    pv = art["plane_variants"]
    for key, color, lw in (("statsys", C_OURS, 2), ("stat", ORANGE, 1.4), ("diag", AQUA, 1.4)):
        _hpd_lines(ax, pv["Om"], pv["OL"], pv[key], color, lw=lw, masses=(0.682689,))
    cv = results["covariance_variants"]["flat_Om"]
    ax.set_xlim(0.05, 0.6)
    ax.set_ylim(0.25, 1.0)
    ax.set_xlabel(OM)
    ax.set_ylabel(OL)
    ax.set_title("What the systematics do (68% regions)")
    ax.legend(title=f"flat ΛCDM {OM} with each:", title_fontsize=7, handles=[
        Line2D([], [], color=C_OURS, lw=2, label=f"statistical + systematic   {_pm(cv['statsys'])}"),
        Line2D([], [], color=ORANGE, lw=1.4, label=f"statistical only   {_pm(cv['stat'])}"),
        Line2D([], [], color=AQUA, lw=1.4, label=f"diagonal only, no correlations   {_pm(cv['diag'])}"),
    ], loc="lower right", fontsize=7)

    wc = art["wcdm"]
    _one_d(axes[1, 0], wc["w"], wc["post"].sum(0), ch["flat_wcdm"]["w"], "Flat wCDM: w",
           "flat_wcdm:w", ".3f")
    axes[1, 0].set_xlim(-1.7, -0.4)
    axes[1, 0].set_xlabel("w")

    cb = art["calibrated"]
    g, dh = cb["h0_grid"], cb["h0_dens"]
    m = (g > 70) & (g < 77.5)
    _one_d(axes[1, 1], g[m], dh[m], ch["shoes_flat_lcdm"]["H0"], f"With Cepheids: {H0}",
           "shoes_flat_lcdm:H0", ".2f")
    axes[1, 1].set_xlabel(f"{H0} (km/s/Mpc)")
    fit = cb["fit"]
    mg = np.linspace(-19.36, -19.13, 400)
    sig = np.sqrt(fit.cov_beta[0, 0])
    dm = (cb["om_wts"][:, None] * np.exp(-0.5 * ((mg[None, :] - fit.M[:, None]) / sig) ** 2)).sum(0)
    _one_d(axes[1, 2], mg, dm, ch["shoes_flat_lcdm"]["M"], "With Cepheids: supernova M",
           "shoes_flat_lcdm:M", ".3f")
    axes[1, 2].set_xlabel("M (mag)")

    fig.legend(handles=[Line2D([], [], color=C_OURS, lw=2, label="this project (grid likelihood)"),
                        Line2D([], [], color=C_PUBLISHED, lw=1.5, label="Pantheon+ MCMC chain")],
               loc="upper right", ncol=2, fontsize=8.5, bbox_to_anchor=(0.905, 0.985))
    fig.suptitle("Checking the likelihood against the Pantheon+ team's own analysis", x=0.125,
                 y=0.985, ha="left")
    caption(fig, "Our posteriors come from exact grid evaluation of a likelihood written for this "
                 "project; the chains are the published CosmoSIS/emcee runs from the data release "
                 "(first 30% discarded as burn-in), with the same flat priors applied to our grids. "
                 "Agreement to a few hundredths of a sigma means the two likelihoods are the same. "
                 "Top right: the same fit using only the covariance's statistical part, or only its "
                 "diagonal -- the second is the mistake the data release warns against.")
    return save(fig, FIGURE_DIR / "07_validation")


def make_all(results: dict, art: dict) -> list:
    written = []
    for fn in (figure_hubble_diagram, figure_acceleration, figure_omega_plane,
               figure_hubble_constant, figure_dark_energy_w, figure_expansion_history,
               figure_validation):
        written += fn(results, art)
        print(f"  {fn.__name__}")
    return written
