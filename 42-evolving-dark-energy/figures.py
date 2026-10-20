"""Every figure in the project, drawn from the cached results and chains.

    python run.py --figures-only

Figures (PNG + SVG in figures/):

01_bao_distances       the sound-horizon ruler from z = 0.3 to the CMB, and DESI vs LCDM
02_lcdm_omega_m        LCDM: DESI's low Omega_m against the CMB and the supernovae
03_w0wa_by_sample      the w0-wa plane, one panel per supernova sample
04_significance        how strongly LCDM is disfavoured, by data combination
05_dark_energy_history w(z) and the dark energy density over time
06_supernova_residuals where the supernova part of the signal comes from
07_validation          our likelihood and posteriors against DESI's published ones
08_cmb_prior           the compressed CMB prior, three ways (and the neutrino offset)
09_own_bao             our own BAO measurement in DESI DR1 luminous red galaxies
"""

from __future__ import annotations

import json
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from scipy.ndimage import gaussian_filter

import analysis as A
import background as bg
import fetch_data as fd
from config import (C_KM_S, DESI_DCHI2, DESI_PARAMS, DOVEKIE_PUBLISHED, FIGURE_DIR,
                    SN_OM_PUBLISHED, UNCALIBRATED_H0, UNCALIBRATED_OMBH2)
from plotstyle import (AXIS, C_BAO, C_CMB, C_COMB, C_DESI, C_OURS, C_SN, EXCLUDED, GRID, INK,
                       INK_SECONDARY, MUTED, ORANGE, RAMPS, SURFACE, VIOLET, BLUE, annotate,
                       caption, lighten, save, use_style)

use_style()
SN = A.SN_SAMPLES
SNL = A.SN_LABEL


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def hpd_levels(H, masses=(0.954, 0.683)):
    flat = np.sort(H.ravel())[::-1]
    cum = np.cumsum(flat)
    cum /= cum[-1]
    return [flat[min(np.searchsorted(cum, m), flat.size - 1)] for m in masses]


def contour2d(ax, x, y, color, *, w=None, xr, yr, bins=140, smooth=2.2, fill=True,
              alpha=0.85, lw=1.1, ls="-", zorder=2):
    H, xe, ye = np.histogram2d(x, y, bins=bins, range=[xr, yr], weights=w)
    H = gaussian_filter(H, smooth)
    l95, l68 = hpd_levels(H)
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    light, dark = RAMPS[color]
    if fill:
        ax.contourf(xc, yc, H.T, levels=[l95, l68, H.max() * 1.01], colors=[light, dark],
                    alpha=alpha, zorder=zorder)
    ax.contour(xc, yc, H.T, levels=[l95, l68], colors=[color, color],
               linewidths=[0.8 * lw, lw], linestyles=ls, zorder=zorder + 0.1)


def chain(model, combo):
    try:
        return A.load_chain(model, combo)
    except FileNotFoundError:
        return None


def col(ch, name):
    return A.chain_column(ch, name)


def model_curves(rec: dict, calibrated: bool, z: np.ndarray):
    """D_M/r_d, D_H/r_d, D_V/r_d of a best-fit universe on a redshift grid."""
    w0, wa = rec.get("w0", -1.0), rec.get("wa", 0.0)
    if calibrated:
        c = bg.Cosmo.from_Om(rec["Om"], rec["H0"] / 100, rec["ombh2"], w0, wa)
        rd = float(A.early_universe().rd(c)[0])
    else:
        c = bg.Cosmo.from_Om(rec["Om"], UNCALIBRATED_H0 / 100, UNCALIBRATED_OMBH2, w0, wa)
        rd = rec["hrd"] / (UNCALIBRATED_H0 / 100)
    grid = bg.DistanceGrid(z, z_max=float(z.max()))
    dm = grid.comoving(c)[0]
    dh = C_KM_S / bg.hubble(z, c)[0]
    return dm / rd, dh / rd, np.cbrt(z * dm**2 * dh) / rd


def bao_points(data: dict) -> dict:
    """Per redshift: D_M/r_d, D_H/r_d (with errors) and D_V/r_d (propagated)."""
    out = {}
    for zz in np.unique(data["z"]):
        idx = np.flatnonzero(data["z"] == zz)
        kinds = list(data["kind"][idx])
        rec = {}
        for k, i in zip(kinds, idx):
            rec[k] = (data["value"][i], np.sqrt(data["cov"][i, i]))
        if "DM_over_rs" in rec and "DH_over_rs" in rec:
            i, j = idx[kinds.index("DM_over_rs")], idx[kinds.index("DH_over_rs")]
            dm, dh = data["value"][i], data["value"][j]
            dv = np.cbrt(zz * dm**2 * dh)
            J = np.array([2 * dv / (3 * dm), dv / (3 * dh)])
            C = data["cov"][np.ix_([i, j], [i, j])]
            rec["DV_over_rs"] = (dv, float(np.sqrt(J @ C @ J)))
        out[float(zz)] = rec
    return out


def pm(x, s, d=3):
    return f"{x:.{d}f} ± {s:.{d}f}"


# ---------------------------------------------------------------------------
# 01  The ruler across cosmic time
# ---------------------------------------------------------------------------


def fig_bao_distances(r: dict) -> None:
    maps = r["maps"]
    ref = maps["cmb-desi"]["lcdm"]  # LCDM from the CMB prior alone
    data = fd.load_bao("dr2")
    pts = bao_points(data)
    early = A.early_universe()

    fig = plt.figure(figsize=(12.5, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.28)
    ax0 = fig.add_subplot(gs[0])

    # (a) one ruler from z = 0.1 to z = 1090
    z = np.geomspace(0.05, 1500, 600)
    dm_ref, _, _ = model_curves(ref, True, z)
    ax0.plot(z, dm_ref, color=C_CMB, lw=1.4, label="LCDM fitted to the CMB alone")
    zs = [k for k, v in pts.items() if "DM_over_rs" in v]
    ax0.errorbar(zs, [pts[k]["DM_over_rs"][0] for k in zs], [pts[k]["DM_over_rs"][1] for k in zs],
                 fmt="o", color=C_BAO, ms=4.5, label="DESI DR2 galaxies, quasars, Lyα forest")
    c_ref = bg.Cosmo.from_Om(ref["Om"], ref["H0"] / 100, ref["ombh2"])
    zstar = float(early("zstar", c_ref.ombh2, c_ref.omch2)[0])
    prior = A.cmb_prior("desi")
    # the CMB measures D_M(z*)/r_*; shown in units of r_d with the model's r_*/r_d
    rs_over_rd = float(early("rstar", c_ref.ombh2, c_ref.omch2)[0] / early.rd(c_ref)[0])
    y_cmb = rs_over_rd / prior.mean[0]
    ax0.errorbar([zstar], [y_cmb], [y_cmb * prior.sigma[0] / prior.mean[0]], fmt="D",
                 color=C_CMB, ms=5.5, label="Planck: acoustic angle θ* of the CMB")
    ax0.set_xscale("log"), ax0.set_yscale("log")
    ax0.set_xlim(0.05, 3000), ax0.set_ylim(1, 200)
    ax0.set_xlabel("redshift z"), ax0.set_ylabel("transverse distance  D_M / r_d")
    ax0.set_title("One ruler, from 3 to 13.8 billion years ago", fontsize=10)
    annotate(ax0, "380,000 years after\nthe Big Bang", (zstar, y_cmb), (40, 120))
    annotate(ax0, "the sound horizon r_d ≈ 147 Mpc\nseen 3–11 billion years ago",
             (0.706, pts[0.706]["DM_over_rs"][0]), (2.0, 3.0))
    ax0.legend(loc="lower right", fontsize=7.3)

    # (b), (c): relative to the CMB's LCDM prediction
    zz = np.linspace(0.02, 2.6, 300)
    base = model_curves(ref, True, zz)
    curves = [
        ("LCDM, DESI alone (Ωm = %.3f)" % maps["desi"]["lcdm"]["Om"], maps["desi"]["lcdm"],
         False, C_BAO, "-"),
        ("LCDM, DESI + CMB", maps["desi_cmb-desi"]["lcdm"], True, INK_SECONDARY, "--"),
        ("w0wa, DESI + CMB + Pantheon+", maps["desi_cmb-desi_pantheonplus"]["w0wa"], True,
         C_COMB, "-"),
        ("w0wa, DESI + CMB + DESY5", maps["desi_cmb-desi_desy5"]["w0wa"], True, C_COMB, ":"),
    ]
    for j, (kind, ylab) in enumerate([("DM_over_rs", "D_M / r_d  ÷  CMB-LCDM"),
                                      ("DH_over_rs", "D_H / r_d  ÷  CMB-LCDM")]):
        ax = fig.add_subplot(gs[j + 1])
        ax.axhline(1, color=C_CMB, lw=1.2)
        for lab, rec, cal, color, ls in curves:
            m = model_curves(rec, cal, zz)
            ax.plot(zz, m[j] / base[j], color=color, ls=ls, lw=1.3, label=lab)
        for k, v in pts.items():
            if kind in v:
                b = model_curves(ref, True, np.array([k]))[j][0]
                ax.errorbar(k, v[kind][0] / b, v[kind][1] / b, fmt="o", color=C_BAO, ms=4.5)
        if j == 0:  # BGS measures only the isotropic D_V
            b = model_curves(ref, True, np.array([0.295]))[2][0]
            v = pts[0.295]["DV_over_rs"]
            ax.errorbar(0.295, v[0] / b, v[1] / b, fmt="s", mfc=SURFACE, color=C_BAO, ms=4.5)
            annotate(ax, "BGS: D_V only", (0.295, v[0] / b), (0.35, 0.945))
        ax.set_xlim(0, 2.6), ax.set_ylim(0.93, 1.05)
        ax.set_xlabel("redshift z"), ax.set_ylabel(ylab)
        ax.set_title(["Across the line of sight", "Along the line of sight"][j], fontsize=10)
        if j == 1:
            ax.legend(loc="lower left", fontsize=6.8)
    caption(fig, "DESI DR2 BAO (DESI Collaboration 2025) divided by the prediction of flat LCDM "
                 "fitted to the compressed Planck CMB prior (violet line). DESI's distances sit "
                 "~2% short of it at z ≈ 0.5–1: within LCDM that asks for a lower matter "
                 "density than the CMB's. The CMB point is D_M(z*)/r_* converted to r_d units "
                 f"with the model's r_*/r_d = {rs_over_rd:.4f}.", y=-0.04)
    save(fig, FIGURE_DIR / "01_bao_distances")


# ---------------------------------------------------------------------------
# 02  LCDM: the matter density from each probe
# ---------------------------------------------------------------------------


def fig_lcdm_omega_m(r: dict) -> None:
    fig, (ax, axf) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                                  gridspec_kw={"width_ratios": [1, 1.05], "wspace": 0.3})
    xr, yr = (0.26, 0.355), (94.5, 105.5)
    sets = [("lcdm", A.Combo(), C_BAO, "DESI BAO"),
            ("lcdm", A.Combo(bao=False, cmb="desi"), C_CMB, "CMB (compressed Planck)"),
            ("lcdm", A.Combo(cmb="desi"), C_COMB, "DESI + CMB")]
    for model, combo, color, lab in sets:
        ch = chain(model, combo)
        if ch is None:
            continue
        contour2d(ax, col(ch, "Om"), col(ch, "hrd"), color, xr=xr, yr=yr,
                  alpha=0.75 if color != C_COMB else 0.9, zorder=3 if color == C_COMB else 2)
    ax.set_xlim(*xr), ax.set_ylim(*yr)
    ax.set_xlabel("matter density Ωm"), ax.set_ylabel("h · r_d  (Mpc)")
    ax.set_title("Flat LCDM: BAO and CMB disagree mildly", fontsize=10)
    ax.legend(handles=[Patch(facecolor=RAMPS[c][1], edgecolor=c, label=l) for _, _, c, l in sets],
              loc="upper right", fontsize=7.5)
    dsig = r["maps"]["desi_cmb-desi"]["lcdm"]
    ax.text(0.03, 0.04, "68% and 95% regions", transform=ax.transAxes, color=MUTED, fontsize=7)

    # forest of Omega_m
    rows = []
    post = r["posteriors"]

    def ours(key):
        p = post.get(key, {}).get("Om")
        return (p["mean"], p["std"]) if p else None

    rows.append(("DESI BAO", ours("lcdm__desi"), DESI_PARAMS[("lcdm", "DESI")]["Om"][:2]))
    rows.append(("CMB prior", ours("lcdm__cmb-desi"), (0.3169, 0.0065)))
    rows.append(("DESI + CMB", ours("lcdm__desi_cmb-desi"), (0.3027, 0.0036)))
    for s in SN:
        rows.append((SNL[s], ours(f"lcdm__{s.lower()}"), SN_OM_PUBLISHED[s][:2]))
    d1 = r.get("d1_flat_om")
    rows.append(("Pantheon+ in D1", tuple(d1) if d1 else None, None))
    y = np.arange(len(rows))[::-1]
    desi_om = rows[0][1]
    if desi_om:
        axf.axvspan(desi_om[0] - desi_om[1], desi_om[0] + desi_om[1], color=lighten(C_BAO, 0.8),
                    zorder=0)
    for yi, (lab, o, p) in zip(y, rows):
        if o:
            axf.errorbar(o[0], yi + 0.12, xerr=o[1], fmt="o", color=C_OURS, ms=5, capsize=0)
        if p:
            axf.errorbar(p[0], yi - 0.12, xerr=p[1], fmt="D", mfc=SURFACE, color=C_DESI, ms=4.5,
                         capsize=0)
    axf.set_yticks(y, [r_[0] for r_ in rows])
    axf.set_xlabel("matter density Ωm (flat LCDM)")
    axf.set_title("Ωm from each probe alone", fontsize=10)
    axf.legend(handles=[Line2D([], [], marker="o", color=C_OURS, ls="", label="this project"),
                        Line2D([], [], marker="D", mfc=SURFACE, color=C_DESI, ls="",
                               label="published (CMB rows: full Planck+ACT likelihood)")],
               loc="upper center", bbox_to_anchor=(0.45, -0.12), ncol=2, fontsize=7)
    axf.grid(axis="y", visible=False)
    caption(fig, "Left: DESI's BAO alone (blue) and the compressed CMB prior alone (violet), "
                 "each interpreted in flat LCDM. Right: Ωm from each probe; supernovae measure "
                 "only the shape of the expansion, so their Ωm needs no calibration. The blue "
                 "band is DESI's 68% range. Supernova samples overlap heavily (Pantheon+, Union3 "
                 "and DES share low-z supernovae), so their rows are not independent.", y=-0.1)
    save(fig, FIGURE_DIR / "02_lcdm_omega_m")


# ---------------------------------------------------------------------------
# 03  The w0-wa plane, one panel per supernova sample
# ---------------------------------------------------------------------------

W0R, WAR = (-1.45, 0.05), (-3.0, 1.4)


def _w0wa_frame(ax, show_labels=True):
    w0 = np.linspace(*W0R, 10)
    ax.fill_between(w0, -w0, WAR[1], color=EXCLUDED, zorder=0, lw=0)
    ax.axvline(-1, color=AXIS, lw=0.8, zorder=1), ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
    ax.plot(-1, 0, marker="*", ms=11, color=INK, mec=SURFACE, mew=0.8, zorder=6)
    if show_labels:
        ax.annotate("cosmological\nconstant", (-1, 0), (-1.42, 0.55), fontsize=7, color=INK,
                    arrowprops={"arrowstyle": "-", "color": INK, "lw": 0.6})
        ax.text(-0.35, 1.05, "excluded:\nw0 + wa > 0", fontsize=6.5, color=MUTED, ha="center")
    ax.set_xlim(*W0R), ax.set_ylim(*WAR)


def fig_w0wa_by_sample(r: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.2), sharex=True, sharey=True,
                             gridspec_kw={"hspace": 0.18, "wspace": 0.08})
    base = chain("w0wa", A.Combo(cmb="desi"))
    maps = r["maps"]
    published = {"PantheonPlus": "DESI+CMB(full)+Pantheon+", "Union3": "DESI+CMB(full)+Union3",
                 "DESY5": "DESI+CMB(full)+DESY5"}
    for ax, s in zip(axes.ravel(), SN):
        _w0wa_frame(ax, s == SN[0])
        sn_only = chain("w0wa", A.Combo(bao=False, sn=s))
        comb = chain("w0wa", A.Combo(sn=s, cmb="desi"))
        if sn_only is not None:
            contour2d(ax, col(sn_only, "w0"), col(sn_only, "wa"), C_SN, xr=W0R, yr=WAR, alpha=0.55)
        if base is not None:
            contour2d(ax, col(base, "w0"), col(base, "wa"), C_BAO, xr=W0R, yr=WAR, alpha=0.55)
        if comb is not None:
            contour2d(ax, col(comb, "w0"), col(comb, "wa"), C_COMB, xr=W0R, yr=WAR, alpha=0.9,
                      zorder=4)
        rec = maps[f"desi_cmb-desi_{s.lower()}"]
        title = f"{SNL[s]}:  ΛCDM disfavoured at {rec['sigma']:.1f}σ"
        ax.set_title(title, fontsize=10)
        pub = (f"DESI with the full CMB likelihood: {DESI_DCHI2[published[s]][1]}σ" if s in published
               else f"Popovic et al. 2025 (CMB: Planck+ACT+SPT): {DOVEKIE_PUBLISHED['sigma']}σ")
        ax.text(0.03, 0.03, f"Δχ² = {rec['dchi2']:.1f} for 2 extra parameters\n{pub}",
                transform=ax.transAxes, ha="left", fontsize=6.8, color=INK_SECONDARY, zorder=7,
                bbox={"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.85, "pad": 2})
    for ax in axes[1]:
        ax.set_xlabel("w0  (equation of state today)")
    for ax in axes[:, 0]:
        ax.set_ylabel("wa  (how fast it changes)")
    handles = [Patch(facecolor=RAMPS[C_BAO][1], edgecolor=C_BAO, label="DESI BAO + CMB prior"),
               Patch(facecolor=RAMPS[C_SN][1], edgecolor=C_SN, label="supernovae alone"),
               Patch(facecolor=RAMPS[C_COMB][1], edgecolor=C_COMB, label="all three combined")]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.955),
               fontsize=8.5)
    fig.suptitle("Is dark energy changing? The answer depends on the supernovae", y=0.995)
    caption(fig, "Dark energy's equation of state w(a) = w0 + wa(1 − a); a cosmological constant "
                 "is the star (w0 = −1, wa = 0). 68% and 95% posterior regions from our own "
                 "likelihood and emcee chains: DESI DR2 BAO with the compressed Planck PR4 prior "
                 "(as DESI ran it) and each supernova sample. Significances are from the best-fit "
                 "Δχ² (DESI's convention). DESI's published 2.8σ / 3.8σ / 4.2σ use the full CMB "
                 "likelihood with lensing, which adds evidence the compressed prior leaves out.",
            y=0.045)
    save(fig, FIGURE_DIR / "03_w0wa_by_sample")


# ---------------------------------------------------------------------------
# 04  Significance by data combination
# ---------------------------------------------------------------------------


def _combo_label(c) -> str:
    parts = (["DESI"] if c.bao else []) + (["CMB"] if c.cmb else []) + ([SNL[c.sn]] if c.sn else [])
    return " + ".join(parts)


def fig_significance(r: dict) -> None:
    maps = r["maps"]
    comp = {row["label"]: row for row in r["desi_comparison"]}
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4),
                                  gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.45})
    rows = [A.Combo()] + [A.Combo(sn=s) for s in SN] + [A.Combo(cmb="desi")] + \
           [A.Combo(sn=s, cmb="desi") for s in SN]
    y = np.arange(len(rows))[::-1].astype(float)
    y[: 5] += 0.6  # a gap between the two groups
    full = {"Pantheon+": 2.8, "Union3": 3.8, "DESY5": 4.2}
    d = 0.22
    for yi, c in zip(y, rows):
        rec = maps[c.key]
        ax.plot(rec["sigma"], yi + d, "o", color=C_OURS, ms=6.5, zorder=4)
        ax.text(rec["sigma"] + 0.1, yi + d, f"{rec['sigma']:.2f}", fontsize=7, color=INK,
                va="center")
        row = comp.get(c.label)
        if row and "desi_files_sigma" in row:
            ax.plot(row["desi_files_sigma"], yi, "D", mfc=SURFACE, color=C_DESI, ms=5.5, mew=1.3,
                    zorder=3)
        if c.cmb and c.sn:
            v = full.get(SNL[c.sn]) if c.sn != "DES-Dovekie" else DOVEKIE_PUBLISHED["sigma"]
            if v:
                ax.plot(v, yi - d, "s", mfc=SURFACE, color=MUTED, ms=5.5, mew=1.3, zorder=3)
                ax.text(v + 0.1, yi - d, f"{v}", fontsize=7, color=MUTED, va="center")
    for s in (2, 3, 4, 5):
        ax.axvline(s, color=GRID, lw=0.8, zorder=0)
    ax.set_yticks(y, [_combo_label(c) for c in rows])
    ax.set_xlim(0, 5.2), ax.set_xlabel("preference for w0wa over LCDM  (σ, 2 extra parameters)")
    ax.grid(axis="y", visible=False)
    ax.set_title("How strongly is a cosmological constant disfavoured?", fontsize=10)
    ax.legend(handles=[
        Line2D([], [], marker="o", color=C_OURS, ls="", label="this project"),
        Line2D([], [], marker="D", mfc=SURFACE, color=C_DESI, ls="", mew=1.3,
               label="DESI's own best fits, same data (released files)"),
        Line2D([], [], marker="s", mfc=SURFACE, color=MUTED, ls="", mew=1.3,
               label="published with the full CMB likelihood")],
        loc="upper center", bbox_to_anchor=(0.4, -0.13), ncol=2, fontsize=7)

    # (b) the CMB prior itself
    variants = [("desi", "PR4 prior as DESI ran it", C_OURS, "o", d),
                ("consistent", "same, neutrino-free mean", ORANGE, "s", 0.0),
                ("d3", "our own prior from D3 (TT only)", VIOLET, "^", -d)]
    rows2 = [None] + list(SN)
    y2 = np.arange(len(rows2))[::-1].astype(float)
    for yi, s in zip(y2, rows2):
        for v, lab, color, mk, off in variants:
            key = A.Combo(sn=s, cmb=v).key
            if key in maps and "sigma" in maps[key]:
                ax2.plot(maps[key]["sigma"], yi + off, mk, color=color, ms=6,
                         mfc=color if mk == "o" else SURFACE, mew=1.3)
                ax2.text(maps[key]["sigma"] + 0.08, yi + off, f"{maps[key]['sigma']:.2f}",
                         fontsize=6.5, color=INK_SECONDARY, va="center")
    ax2.set_yticks(y2, ["DESI + CMB"] + [f"DESI + CMB + {SNL[s]}" for s in SN])
    for s in (2, 3, 4):
        ax2.axvline(s, color=GRID, lw=0.8, zorder=0)
    ax2.set_xlim(1.5, 4.3), ax2.set_xlabel("σ")
    ax2.grid(axis="y", visible=False)
    ax2.set_title("…and which CMB prior is used", fontsize=10)
    ax2.legend(handles=[Line2D([], [], marker=mk, color=c, mfc=c if mk == "o" else SURFACE,
                               ls="", mew=1.3, label=l) for _, l, c, mk, _ in variants],
               loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=1, fontsize=7)
    caption(fig, "Significance = Δχ² between the best-fit w0waCDM and LCDM, converted to "
                 "Gaussian σ for two extra parameters (DESI's convention). Left: our "
                 "likelihood against DESI's released best-fit files for the same data "
                 "combinations, and DESI's headline numbers that use the full Planck+ACT CMB "
                 "likelihood (Popovic et al. 2025 for DES-Dovekie, with Planck+ACT+SPT). Right: "
                 "the same fits with three versions of the compressed CMB prior.", y=-0.2)
    save(fig, FIGURE_DIR / "04_significance")


# ---------------------------------------------------------------------------
# 05  Dark energy over time
# ---------------------------------------------------------------------------


def fig_dark_energy_history(r: dict) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.2), sharex=True, sharey="row",
                             gridspec_kw={"hspace": 0.12, "wspace": 0.08})
    z = np.linspace(0, 2.5, 126)
    for j, s in enumerate(SN):
        combo = A.Combo(sn=s, cmb="desi")
        if chain("w0wa", combo) is None:
            continue
        b = A.w_band("w0wa", combo, z=z)
        for i, (key, ref) in enumerate((("w", -1.0), ("rho", 1.0))):
            ax = axes[i, j]
            q = b[key]
            ax.fill_between(z, q[0], q[4], color=RAMPS[INK][0], lw=0, label="95%")
            ax.fill_between(z, q[1], q[3], color=RAMPS[INK][1], lw=0, label="68%")
            ax.plot(z, q[2], color=INK, lw=1.2, label="median")
            ax.axhline(ref, color=INK_SECONDARY, lw=1.0, ls="--")
            if key == "w":
                ax.axhspan(-3, -1, color=lighten(ORANGE, 0.9), zorder=0, lw=0)
                ax.set_ylim(-2.0, -0.2)
            else:
                ax.set_ylim(0.3, 1.6)
        axes[0, j].set_title(f"DESI + CMB + {SNL[s]}", fontsize=9.5)
        axes[1, j].set_xlabel("redshift z")
    axes[0, 0].set_ylabel("equation of state w(z)")
    axes[1, 0].set_ylabel("dark energy density ρ(z) / ρ(today)")
    axes[0, 0].text(1.25, -1.9, "'phantom': w < −1", fontsize=7, color=ORANGE, ha="center")
    axes[0, 0].text(2.45, -0.97, "cosmological constant", fontsize=6.5, color=INK_SECONDARY,
                    ha="right", va="bottom")
    axes[1, 0].legend(loc="upper left", fontsize=7)
    fig.suptitle("What the w0–wa fits say dark energy did", y=0.97)
    caption(fig, "Posterior bands of w(z) = w0 + wa z/(1+z) and of the dark energy density it "
                 "implies, ρ ∝ (1+z)^{3(1+w0+wa)} exp(−3 wa z/(1+z)). Read right to left (forward in "
                 "time): in every combination the density grows while w < −1, peaks at z ≈ 0.4–0.5 "
                 "(4–5 billion years ago) and has been fading since. The 'phantom' early phase may be "
                 "an artefact of the two-parameter form: w0–wa is a convenient parametrization, not a "
                 "physical model.", y=0.02)
    save(fig, FIGURE_DIR / "05_dark_energy_history")


# ---------------------------------------------------------------------------
# 06  Supernova residuals
# ---------------------------------------------------------------------------


def binned_gls(sn, resid, edges):
    which = np.digitize(sn.z_hd, edges) - 1
    used = [k for k in range(len(edges) - 1) if np.sum(which == k) > 0]
    B = np.zeros((sn.n, len(used)))
    for j, k in enumerate(used):
        B[which == k, j] = 1.0
    Aw = sn.whiten @ B
    cov = np.linalg.inv(Aw.T @ Aw)
    val = cov @ Aw.T @ (sn.whiten @ resid)
    zb = np.array([np.median(sn.z_hd[which == k]) for k in used])
    return zb, val, cov


def fig_supernova_residuals(r: dict) -> None:
    maps = r["maps"]
    ref = maps["desi_cmb-desi"]["lcdm"]
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.9), sharey=True,
                             gridspec_kw={"wspace": 0.07})
    edges = np.array([0.01, 0.02, 0.03, 0.045, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                      0.85, 1.0, 1.4, 2.4])
    zc = np.geomspace(0.012, 2.3, 300)
    for ax, s in zip(axes, SN):
        sn = A.sn_likelihood(s)
        P = A.problem("lcdm", A.Combo(cmb="desi"))
        c_ref = bg.Cosmo.from_Om(ref["Om"], ref["H0"] / 100, ref["ombh2"])
        grid = bg.DistanceGrid(np.unique(np.concatenate([sn.z_hd, zc])), z_max=2.4)
        dc_all = grid.comoving(c_ref)[0]
        zt = grid.z_targets
        dm_sn = np.interp(sn.z_hd, zt, dc_all)
        resid = sn.residuals(dm_sn)
        if s == "Union3":  # already binned (22 spline nodes of the UNITY1.5 model)
            zb, val = sn.z_hd, resid
            cov = np.linalg.inv(sn.whiten.T @ sn.whiten)
        else:
            zb, val, cov = binned_gls(sn, resid, edges)
        # errors after removing the fitted common offset, which cancels in any comparison
        icov = np.linalg.inv(cov)
        err = np.sqrt(np.clip(np.diag(cov) - 1.0 / icov.sum(), 0, None))
        ax.axhline(0, color=INK_SECONDARY, lw=1.0, ls="--")
        ax.errorbar(zb, val, err, fmt="o", color=C_SN, ms=4, label=f"{SNL[s]}, binned", zorder=3)
        # the w0wa best fit for DESI + CMB + this sample, as a residual from the reference
        rec = maps[f"desi_cmb-desi_{s.lower()}"]["w0wa"]
        c_w = bg.Cosmo.from_Om(rec["Om"], rec["H0"] / 100, rec["ombh2"], rec["w0"], rec["wa"])
        dc_w = grid.comoving(c_w)[0]
        _, off_ref = sn.chi2(dm_sn[None, :])
        _, off_w = sn.chi2(np.interp(sn.z_hd, zt, dc_w)[None, :])
        mu = lambda d: 5 * np.log10((1 + zc) * np.interp(zc, zt, d))
        ax.plot(zc, mu(dc_w) - mu(dc_all) + (off_w[0] - off_ref[0]), color=C_COMB, lw=1.4,
                label="best w0wa (DESI + CMB + this sample)")
        ax.set_xscale("log"), ax.set_xlim(0.01, 2.5), ax.set_ylim(-0.15, 0.15)
        ax.set_title(f"{SNL[s]} ({'22 spline nodes' if s == 'Union3' else f'{sn.n} supernovae'})",
                     fontsize=9.5)
        ax.set_xlabel("redshift z")
        if s == SN[0]:
            ax.set_ylabel("Hubble residual  Δμ  (mag)")
            ax.legend(loc="lower left", fontsize=6.8)
    fig.suptitle("Where the supernova evidence comes from: nearby vs distant supernovae", y=1.02)
    caption(fig, "Distance moduli minus the flat-LCDM best fit to DESI + CMB (dashed line), each "
                 "sample with its own best offset (the supernova magnitude is not calibrated "
                 "here). Bins are generalised-least-squares averages with the full covariance; error "
                 "bars exclude the common offset. Where nearby (z < 0.1) supernovae sit above the "
                 "distant ones (fainter, relative to LCDM) the fit is pulled toward w0 > −1, wa < 0: "
                 "strongly in DESY5, less after its Dovekie recalibration, least in Pantheon+. Bins "
                 "are correlated; fits use the individual supernovae.",
            y=-0.04)
    save(fig, FIGURE_DIR / "06_supernova_residuals")


# ---------------------------------------------------------------------------
# 07  Validation against DESI
# ---------------------------------------------------------------------------


def fig_validation(r: dict) -> None:
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(14, 5.2),
                                     gridspec_kw={"width_ratios": [1, 0.9, 1.1], "wspace": 0.55})
    rows = r["validation"]["rows"]
    lab = [f"{x['model'].replace('lcdm', 'LCDM').replace('w0wa', 'w0wa')}: {x['label']}"
           for x in rows]
    y = np.arange(len(rows))[::-1]
    for yi, x in zip(y, rows):
        a1.plot(x["chi2"] - x["desi_chi2"], yi, "o", color=C_OURS, ms=5)
        if "chi2_cmb_if_bcnu" in x:
            alt = x["chi2"] - x["chi2_cmb"] + x["chi2_cmb_if_bcnu"] - x["desi_chi2"]
            a1.plot(alt, yi, "x", color=ORANGE, ms=5, mew=1.3)
    a1.axvline(0, color=INK_SECONDARY, lw=0.8)
    a1.set_yticks(y, lab, fontsize=6.6)
    a1.set_xlim(-1.45, 0.15)
    a1.set_xlabel("our χ² − DESI's χ², at DESI's best fits")
    a1.set_title("Same parameters, same χ²", fontsize=10)
    a1.grid(axis="y", visible=False)
    a1.legend(handles=[Line2D([], [], marker="o", color=C_OURS, ls="", label="our likelihood"),
                       Line2D([], [], marker="x", color=ORANGE, ls="", mew=1.3,
                              label="if the CMB prior compared ω_b+ω_c+ω_ν")],
              loc="lower left", fontsize=6.8)

    comp = sorted([x for x in r["desi_comparison"] if "desi_files_dchi2" in x],
                  key=lambda x: x["dchi2"])
    for i, x in enumerate(comp):
        a2.plot(x["desi_files_dchi2"], x["dchi2"], "o", color=C_OURS, ms=6)
        a2.annotate(re.sub(r"\+(?=[A-Z(])", " + ", x["label"]),
                    (x["desi_files_dchi2"], x["dchi2"]),
                    (8, 6) if i % 2 else (8, -9), textcoords="offset points", fontsize=6.3,
                    color=INK_SECONDARY)
    lim = (0, 22)
    a2.plot(lim, lim, color=AXIS, lw=0.8, zorder=0)
    a2.set_xlim(*lim), a2.set_ylim(*lim)
    a2.set_xlabel("DESI's Δχ² (their best-fit files)"), a2.set_ylabel("our Δχ² (our own fits)")
    a2.set_title("Same fits, same significance", fontsize=10)

    # posterior means: ours vs DESI chains
    pairs = [("lcdm__desi", "base/desi-bao-all"), ("w0wa__desi", "base_w_wa/desi-bao-all"),
             ("w0wa__desi_pantheonplus", "base_w_wa/desi-bao-all_pantheonplus"),
             ("w0wa__desi_union3", "base_w_wa/desi-bao-all_union3"),
             ("w0wa__desi_desy5", "base_w_wa/desi-bao-all_desy5sn"),
             ("w0wa__desi_cmb-desi", "base_w_wa/desi-bao-all_CMB-compressed-theta-ombh2-ombch2"),
             ("w0wa__desi_cmb-desi_desy5",
              "base_w_wa/desi-bao-all_CMB-compressed-theta-ombh2-ombch2_desy5sn"),
             ("lcdm__pantheonplus", "base/pantheonplus"), ("lcdm__union3", "base/union3"),
             ("lcdm__desy5", "base/desy5sn")]
    names = {"lcdm__desi": "LCDM DESI", "w0wa__desi": "w0wa DESI"}
    post, dp = r["posteriors"], r["desi_posteriors"]
    yy = 0
    ylabels = []
    marks = {"Om": ("o", "Ωm"), "w0": ("s", "w0"), "wa": ("^", "wa"), "hrd": ("v", "h r_d"),
             "H0": ("P", "H0")}
    for ours_key, desi_key in pairs[::-1]:
        if ours_key not in post or desi_key not in dp:
            continue
        for p, (mk, _) in marks.items():
            o, d = post[ours_key].get(p), dp[desi_key].get(p)
            if o and d:
                a3.plot((o["mean"] - d["mean"]) / d["std"], yy, mk, color=INK, ms=5,
                        mfc=SURFACE if mk != "o" else INK)
        ylabels.append(names.get(ours_key, ours_key.split("__")[0].replace("lcdm", "LCDM") + " " +
                       post[ours_key]["label"]))
        yy += 1
    a3.set_yticks(range(yy), ylabels, fontsize=7)
    a3.axvline(0, color=INK_SECONDARY, lw=0.8)
    a3.axvspan(-0.1, 0.1, color=lighten(C_OURS, 0.85), zorder=0)
    a3.set_xlim(-0.5, 0.5)
    a3.set_xlabel("(our mean − DESI chain mean) / DESI σ")
    a3.set_title("Same posteriors", fontsize=10)
    a3.grid(axis="y", visible=False)
    a3.legend(handles=[Line2D([], [], marker=mk, color=INK, mfc=SURFACE if mk != "o" else INK,
                              ls="", label=l) for mk, l in marks.values()],
              loc="lower right", fontsize=6.8, ncol=2)
    caption(fig, "Left: our χ² at each of DESI's published best-fit points, minus DESI's. "
                 "Without the CMB they agree to 10⁻⁴; with it to a few hundredths (CAMB 2.0.4 "
                 "vs DESI's CAMB 1.5.4 put θ* ~10⁻⁵ apart). The crosses show the CMB prior "
                 "evaluated with the neutrinos included in ω_bc: clearly not what DESI did. "
                 "Middle: our own best-fit Δχ² against DESI's files. Right: our posterior means "
                 "against DESI's released chains (30% burn-in), in units of DESI's σ.", y=-0.03)
    save(fig, FIGURE_DIR / "07_validation")


# ---------------------------------------------------------------------------
# 08  The compressed CMB prior
# ---------------------------------------------------------------------------


def _ellipse(ax, mean, cov, color, ls="-", lw=1.4, nsig=1.52, fill=False):
    vals, vecs = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
    w, h = 2 * nsig * np.sqrt(vals[1]), 2 * nsig * np.sqrt(vals[0])
    e = Ellipse(mean, w, h, angle=ang, fill=fill, edgecolor=color, lw=lw, ls=ls,
                facecolor=lighten(color, 0.85) if fill else "none")
    ax.add_patch(e)


def fig_cmb_prior(r: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), gridspec_kw={"wspace": 0.3})
    ch = A.planck_early_chains()
    rng = np.random.default_rng(2)
    idx = rng.choice(ch["weight"].size, 2500, replace=False)
    th = ch["thetastar"][idx]  # stored as 100 theta_*
    obc = (ch["ombh2"] + ch["omch2"])[idx]
    ob = ch["ombh2"][idx]
    variants = [("desi", "as DESI ran it", BLUE, "-"), ("consistent", "neutrino-free mean", ORANGE, "--"),
                ("d3", "our D3 fit (TT only)", VIOLET, "-")]
    for ax, (xv, xlab) in zip(axes, [(0, "100 θ*"), (1, "ω_b")]):
        xs = th if xv == 0 else ob
        ax.scatter(xs, obc, s=2, color=AXIS, lw=0, label="Planck PR4 early-universe chain",
                   zorder=1)
        for v, lab, color, ls in variants:
            try:
                p = A.cmb_prior(v)
            except FileNotFoundError:
                continue
            sel = [xv, 2]
            mean = p.mean[sel] * np.array([100 if xv == 0 else 1, 1])
            scale = np.diag([100 if xv == 0 else 1, 1])
            cov = scale @ p.cov[np.ix_(sel, sel)] @ scale
            _ellipse(ax, mean, cov, color, ls=ls)
            ax.plot(*mean, "o", color=color, ms=3.5)
        ax.set_xlabel(xlab), ax.set_ylabel("ω_b + ω_c  (baryons + cold dark matter)")
    d = A.cmb_prior("desi")
    c = A.cmb_prior("consistent")
    annotate(axes[0], f"DESI's mean is {d.mean[2] - c.mean[2]:.5f} higher:\n"
                      "the density of a 0.06 eV neutrino", (100 * d.mean[0], d.mean[2]),
             (100 * d.mean[0] + 0.00032, d.mean[2] - 0.0037))
    axes[0].xaxis.set_major_locator(plt.MaxNLocator(5))
    axes[0].legend(handles=[Line2D([], [], color=cc, ls=ls, label=l) for _, l, cc, ls in variants]
                   + [Line2D([], [], marker="o", color=AXIS, ls="", ms=3,
                             label="Lemos & Lewis chain samples")], fontsize=7, loc="lower left")
    axes[0].set_title("The compressed CMB prior", fontsize=10)
    axes[1].set_title("…and its baryon direction", fontsize=10)
    caption(fig, "Ellipses: 68% regions of the three Gaussian CMB priors on (θ*, ω_b, ω_b + ω_c). "
                 "Grey points: the Planck PR4 'early-universe' chain of Lemos & Lewis (2023) "
                 "that DESI's prior summarises. DESI's third mean equals the chain average of "
                 "ω_m h², which includes the massive neutrino; its likelihood compares it with "
                 "ω_b + ω_c, which does not. Our D3 prior is wider (temperature only, τ fixed).",
            y=-0.03)
    save(fig, FIGURE_DIR / "08_cmb_prior")


# ---------------------------------------------------------------------------
# 09  Our own BAO measurement
# ---------------------------------------------------------------------------


def fig_own_bao(r: dict) -> None:
    import bao_measurement as BM

    path = BM.OUT_DIR / "bao_results.json"
    if not path.exists():
        print("  (skipping 09_own_bao: run `python bao_measurement.py` first)")
        return
    res = json.loads(path.read_text())
    fig = plt.figure(figsize=(14, 8.6))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.3)
    top = [fig.add_subplot(gs[0, i]) for i in range(3)]
    bot = [fig.add_subplot(gs[1, i]) for i in range(3)]
    have_post = all("post" in rr for rr in res.values())
    for ax, (key, rr) in zip(top[:2], res.items()):
        for stage, color, mfc, lab in (("pre", MUTED, SURFACE, "before reconstruction"),
                                       ("post", INK, INK, "after reconstruction")):
            f = BM.OUT_DIR / f"xi_{stage}_z{key}.npz"
            if stage not in rr or not f.exists():
                continue
            with np.load(f) as z:
                s, xi0, cov0, sm, model = z["s"], z["xi0"], z["cov0"], z["s_model"], z["model"]
            m = (s > 20) & (s < 190)
            err = np.sqrt(np.diag(cov0))
            off = -0.8 if stage == "pre" else 0.8
            ax.errorbar(s[m] + off, s[m] ** 2 * xi0[m], s[m] ** 2 * err[m], fmt="o", color=color,
                        mfc=mfc, ms=3.3, lw=1.0, label=lab)
            ax.plot(sm, sm**2 * model, color=C_BAO if stage == "post" else MUTED,
                    lw=1.6 if stage == "post" else 1.1, ls="-" if stage == "post" else "--")
        ax.axvspan(*BM.FIT_RANGE, color=lighten(C_BAO, 0.94), zorder=0)
        ax.set_xlabel("pair separation s  (h⁻¹ Mpc)"), ax.set_ylabel("s² ξ₀(s)  (h⁻² Mpc²)")
        ax.set_title(f"{rr['zmin']} < z < {rr['zmax']}:  {rr['n_data']:,} galaxies", fontsize=10)
        lines = [f"BAO peak detected at {rr['pre']['detection_sigma']:.1f}σ before"]
        if "post" in rr:
            lines.append(f"and {rr['post']['detection_sigma']:.1f}σ after reconstruction")
        ax.text(0.03, 0.04, "\n".join(lines), transform=ax.transAxes, fontsize=7.3, color=INK)
        if ax is top[0]:
            ax.legend(loc="upper right", fontsize=7)
    ax = top[2]
    for (key, rr), mk in zip(res.items(), ("o", "s")):
        stage = "post" if "post" in rr else "pre"
        with np.load(BM.OUT_DIR / f"xi_{stage}_z{key}.npz") as z:
            s8, xi2, cov, sa, m2 = z["s8"], z["xi2_8"], z["cov02_8"], z["s_aniso"], z["model2_aniso"]
        n = len(s8)
        err = np.sqrt(np.diag(cov))[n:]
        m = (s8 > 20) & (s8 < 190)
        ax.errorbar(s8[m], s8[m] ** 2 * xi2[m], s8[m] ** 2 * err[m], fmt=mk, color=INK,
                    mfc=SURFACE if mk == "s" else INK, ms=3.5, lw=1.0,
                    label=f"{rr['zmin']}–{rr['zmax']}")
        ax.plot(sa, sa**2 * m2, color=C_BAO, lw=1.4, ls="-" if mk == "o" else "--")
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.axvspan(*BM.ANISO_RANGE, color=lighten(C_BAO, 0.94), zorder=0)
    ax.set_xlabel("pair separation s  (h⁻¹ Mpc)"), ax.set_ylabel("s² ξ₂(s)  (h⁻² Mpc²)")
    ax.set_title("Quadrupole" + (" (after reconstruction)" if have_post else ""), fontsize=10)
    ax.legend(loc="lower left", fontsize=7, title="redshift slice", title_fontsize=7)

    dr1, dr2 = bao_points(fd.load_bao("dr1")), bao_points(fd.load_bao("dr2"))
    ref = r["maps"]["cmb-desi"]["lcdm"]
    specs = [("DV_over_rs", "dv_rd", 2, "D_V / r_d  (the average)"),
             ("DM_over_rs", "dm_rd", 0, "D_M / r_d  (across the line of sight)"),
             ("DH_over_rs", "dh_rd", 1, "D_H / r_d  (along the line of sight)")]
    for ax, (kind, ours, j, title) in zip(bot, specs):
        for i, (key, rr) in enumerate(res.items()):
            zz = rr["z_eff"]
            b = model_curves(ref, True, np.array([zz]))[j][0]
            for stage, mfc, dz in (("pre", SURFACE, -0.012), ("post", C_OURS, 0.0)):
                if stage not in rr:
                    continue
                q = rr[stage]
                key2 = "dv_rd_2d" if ours == "dv_rd" else ours  # D_V from the 2D fit, as DESI
                yerr = q[f"{key2}_err"] / b
                ax.errorbar(zz + dz, q[key2] / b, yerr, fmt="o", color=C_OURS, mfc=mfc, ms=5.5,
                            mew=1.3, label=(f"ours, {stage}-reconstruction") if i == 0 else None)
        for lab, pts, mk, dz in (("DESI DR1", dr1, "D", 0.014), ("DESI DR2", dr2, "s", 0.026)):
            first = True
            for zz in (0.51, 0.706):
                if zz in pts and kind in pts[zz]:
                    v, e = pts[zz][kind]
                    b = model_curves(ref, True, np.array([zz]))[j][0]
                    ax.errorbar(zz + dz, v / b, e / b, fmt=mk, color=C_DESI, mfc=SURFACE, ms=5,
                                mew=1.2, label=lab if first else None)
                    first = False
        ax.axhline(1, color=C_CMB, lw=1.0)
        ax.set_xlim(0.44, 0.78), ax.set_ylim(0.86, 1.10) if j == 1 else ax.set_ylim(0.92, 1.06)
        ax.set_xlabel("effective redshift"), ax.set_ylabel("÷ CMB-LCDM prediction")
        ax.set_title(title, fontsize=10)
        if j == 2:
            ax.legend(loc="lower left", fontsize=6.6)
    caption(fig, "Top: the correlation-function monopole of 1.28 million DESI DR1 luminous red "
                 "galaxies before (grey) and after (black) our reconstruction, with best-fit "
                 "templates (fit range shaded), and the quadrupole with the two-direction fit. "
                 "Bottom: the ruler in three forms, divided by the CMB-LCDM prediction, against "
                 "DESI's published (reconstructed) values; points shifted slightly in z for clarity. "
                 "Errors: delete-one jackknife over 100 sky regions (DESI's use simulations).",
            y=0.02)
    save(fig, FIGURE_DIR / "09_own_bao")


def make_all(results: dict | None = None) -> None:
    r = results or json.loads(A.RESULTS_PATH.read_text())
    for f in (fig_bao_distances, fig_lcdm_omega_m, fig_w0wa_by_sample, fig_significance,
              fig_dark_energy_history, fig_supernova_residuals, fig_validation, fig_cmb_prior,
              fig_own_bao):
        print(f"  {f.__name__}")
        f(r)


if __name__ == "__main__":
    make_all()
