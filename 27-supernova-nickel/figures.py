"""Every figure, drawn from the cached analysis results.

    python run.py --figures-only
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import config as C
import gammarays as G
import literature as LIT
import nickel as N
import plotstyle as S

S.use_style()
LABEL = {sn: C.SUPERNOVAE[sn]["label"] for sn in C.SUPERNOVAE}
FIG = C.FIGURE_DIR


def _log_axis(ax, lo, hi):
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)


def _end_label(ax, x, y, text, color, dx=4, dy=0, ha="left", va="center", **kw):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", color=S.INK,
                fontsize=8, ha=ha, va=va, **kw)


def _points(ax, sn, t, lum, err=None, open_mask=None, ms=3.2, zorder=3, alpha=1.0):
    col = S.SN_COLORS[sn]
    mk = S.SN_MARKERS[sn]
    t, lum = np.asarray(t), np.asarray(lum)
    filled = np.ones(len(t), bool) if open_mask is None else ~np.asarray(open_mask)
    ax.plot(t[filled], lum[filled], mk, color=col, ms=ms, mec=S.SURFACE, mew=0.4,
            zorder=zorder, alpha=alpha, ls="none")
    if open_mask is not None and np.any(~filled):
        ax.plot(t[~filled], lum[~filled], mk, mfc=S.SURFACE, mec=col, ms=ms, mew=0.8,
                zorder=zorder, alpha=alpha, ls="none")


# ---------------------------------------------------------------------------
# 1. The keeper: four supernovae against the radioactive power
# ---------------------------------------------------------------------------


def _slope_angle(ax, x0, y0, x1, y1):
    """Angle (deg) of the segment (x0, y0)-(x1, y1) as drawn on screen."""
    p0 = ax.transData.transform((x0, y0))
    p1 = ax.transData.transform((x1, y1))
    return float(np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0])))


def fig_radioactive_glow(r: dict):
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    _log_axis(ax, 2e38, 5e43)
    ax.set_xlim(0, 800)

    # Type Ia: data + Arnett's model with gamma-ray leakage.
    for sn in ("SN2011fe", "SN2014J"):
        tab = r[sn]["bolo"]["table"]
        af = r[sn]["fit"]["arnett_fit"]
        _points(ax, sn, tab["phase"], tab["L"], open_mask=tab["nir_source"] != "measured", ms=3.0)
        tt = np.linspace(0.5, tab["phase"].max(), 800)
        ax.plot(tt, N.arnett_light_curve(tt, af["m_ni"], af["tau_m"], af["t0"]),
                color=S.SN_COLORS[sn], lw=1.5, zorder=4)

    # Type II: published UVOIR of SN 1987A, our UVOIR of SN 2023ixf.
    sn87 = r["SN1987A"]["bolo"]
    late = sn87["published_late"][sn87["published_late"]["phase"] > sn87["published"]["phase"].max()]
    _points(ax, "SN1987A", pd.concat([sn87["published"]["phase"], late["phase"]]),
            pd.concat([sn87["published"]["L"], late["L"]]), ms=3.4)
    ixf = r["SN2023ixf"]["bolo"]["table"]
    _points(ax, "SN2023ixf", ixf["phase"], ixf["L"], open_mask=ixf["nir_source"] != "measured", ms=3.4)
    f87 = r["SN1987A"]["fit"]["tail_with_escape"]
    fix = r["SN2023ixf"]["fit"]["tail_measured_nir"]
    for sn, f, t_lo, t_hi in (("SN1987A", f87, 120.0, 800.0), ("SN2023ixf", fix, 95.0, 365.0)):
        tt = np.linspace(t_lo, t_hi, 400)
        ax.plot(tt, f["m_ni"] * N.deposition(tt, f["t0"]), color=S.SN_COLORS[sn], lw=1.5, zorder=4)

    # Radioactive power if every gamma ray were trapped, for the two families.
    m_ia = r["SN2011fe"]["fit"]["arnett_fit"]["m_ni"]
    m_ii = f87["m_ni"]
    t = np.linspace(0.5, 800, 1600)
    ax.plot(t, m_ia * N.deposition(t, None), color=S.C_MODEL, lw=0.9, ls=(0, (4, 3)), zorder=2)
    t2 = t[t >= 50]
    ax.plot(t2, m_ii * N.deposition(t2, None), color=S.C_MODEL, lw=0.9, ls=(0, (4, 3)), zorder=2)
    fig.canvas.draw()  # the on-screen slope of the lines needs the final layout
    for m, xa, txt in ((m_ia, 430.0, f"{m_ia:.2f} solar masses of nickel-56, every gamma ray trapped"),
                       (m_ii, 520.0, f"{m_ii:.3f} solar masses, every gamma ray trapped")):
        xb = xa + 60
        ya, yb = m * N.deposition(xa, None), m * N.deposition(xb, None)
        ax.text(xa, ya * 1.15, txt, rotation=_slope_angle(ax, xa, ya, xb, yb), rotation_mode="anchor",
                fontsize=7.8, color=S.INK_SECONDARY, ha="left", va="bottom", transform_rotates_text=False)

    ax.set_xlabel("Days since explosion")
    ax.set_ylabel("Luminosity (erg/s)")
    ax.set_title("What makes a supernova shine: the radioactive decay of nickel-56", pad=14)

    # Direct labels at the ends of each supernova's data.
    fe = r["SN2011fe"]["bolo"]["table"]
    j14 = r["SN2014J"]["bolo"]["table"]
    _end_label(ax, fe["phase"].iloc[-1] - 30, fe["L"].iloc[-1], "SN 2011fe\nType Ia", S.INK, dx=0, dy=-8,
               ha="right", va="top")
    _end_label(ax, j14["phase"].iloc[-1], j14["L"].iloc[-1], "SN 2014J\nType Ia", S.INK, dx=6)
    _end_label(ax, ixf["phase"].iloc[-1], ixf["L"].iloc[-1], "SN 2023ixf\nType II", S.INK, dx=6, dy=-4, va="top")
    _end_label(ax, 700, f87["m_ni"] * N.deposition(700.0, f87["t0"]), "SN 1987A\nType II", S.INK,
               dx=-6, dy=-10, ha="right", va="top")

    # What the shapes mean.
    pk = r["SN2011fe"]["fit"]["peak"]
    S.annotate(ax, "Arnett's rule: at peak, the light\ngiven out equals the radioactive\npower going in",
               (pk["t_peak"], pk["L_peak"]), (95, 2.2e43))
    S.annotate(ax, "SN 2023ixf's plateau: energy left by the\nshock (hydrogen recombining), not nickel",
               (72, 1.45e42), (150, 5.0e42))
    S.annotate(ax, "SN 1987A's slow rise: a compact blue\nstar whose shock heat faded fast", (32, 3.3e41), (75, 1.5e40))

    # The cobalt slope.
    x0, x1 = 600, 700
    y0 = 2.0e41
    ax.plot([x0, x1], [y0, y0 * np.exp(-(x1 - x0) / N.TAU_CO)], color=S.INK, lw=1.2)
    ax.text(x0 + 4, y0 * 1.25, "cobalt-56 decay:\n0.98 magnitudes per 100 days", fontsize=8, color=S.INK, va="bottom")

    handles = [
        Line2D([], [], color=S.C_MODEL, ls=(0, (4, 3)), lw=0.9, label="radioactive power, every gamma ray trapped"),
        Line2D([], [], color=S.INK_SECONDARY, lw=1.5, label="fitted: diffusion + gamma rays escaping"),
        Line2D([], [], marker="o", ls="none", color=S.INK_SECONDARY, ms=4, label="near-infrared measured"),
        Line2D([], [], marker="o", ls="none", mfc=S.SURFACE, mec=S.INK_SECONDARY, ms=4,
               label="near-infrared share estimated"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7.5)
    S.caption(fig, "Points: SN 2011fe, SN 2014J and SN 2023ixf from multi-band photometry (Open Supernova Catalog; "
              "Hsu et al. 2025; Singh et al. 2024), integrated over wavelength in this project; SN 1987A: published "
              "UVOIR (Suntzeff et al. 1991, 1992; Bouchet et al. 1991). Type Ia days count from first light.")
    return S.save(fig, FIG / "01_radioactive_glow")


# ---------------------------------------------------------------------------
# 2. The nickel clock
# ---------------------------------------------------------------------------


def fig_nickel_clock(r: dict):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    t = np.linspace(0, 400, 2001)
    n = N.nuclei(t, 1.0)
    tot = N.N_PER_MSUN
    for key, col, lab in (("Ni", S.BLUE, "nickel-56"), ("Co", S.ORANGE, "cobalt-56"), ("Fe", S.AQUA, "iron-56")):
        a1.plot(t, n[key] / tot, color=col, lw=1.8)
    a1.set_xlim(0, 400)
    a1.set_ylim(0, 1.02)
    a1.set_xlabel("Days since explosion")
    a1.set_ylabel("Fraction of the original nickel atoms")
    a1.set_title("The decay chain")
    _end_label(a1, 12, 0.22, "nickel-56\n(half-life 6.1 days)", S.INK, dx=6, dy=0)
    t_co = np.log(N.LAMBDA_NI / N.LAMBDA_CO) / (N.LAMBDA_NI - N.LAMBDA_CO)
    co_max = N.nuclei(t_co, 1.0)["Co"] / tot
    _end_label(a1, t_co, co_max, f"cobalt-56 peaks on day {t_co:.0f}\n(half-life 77.2 days)", S.INK, dx=8, dy=6, va="bottom")
    _end_label(a1, 300, 0.92, "iron-56 (stable)", S.INK, dx=-6, dy=8, ha="right", va="bottom")
    t2 = np.linspace(0.5, 600, 2400)
    p = N.power(t2)
    a2.plot(t2, p["total"], color=S.INK, lw=1.8)
    a2.plot(t2, p["ni_gamma"], color=S.BLUE, lw=1.2)
    a2.plot(t2, p["co_gamma"], color=S.ORANGE, lw=1.2)
    a2.plot(t2, p["co_positron"], color=S.MAGENTA, lw=1.2)
    for t0, lab, dy in ((35.0, "T0 = 35 d (Type Ia)", -10), (200.0, "T0 = 200 d (Type II)", 6)):
        a2.plot(t2, N.deposition(t2, t0), color=S.INK_SECONDARY, lw=1.0, ls=(0, (4, 3)))
    _log_axis(a2, 1e39, 2e44)
    a2.set_xlim(0, 600)
    a2.set_xlabel("Days since explosion")
    a2.set_ylabel("Heating per solar mass of nickel (erg/s)")
    a2.set_title("The power it releases")
    _end_label(a2, 8, 6e43, "total", S.INK, dx=8, dy=4)
    _end_label(a2, 98, 1.2e39, "nickel's gamma rays", S.INK, dx=6, dy=4, va="bottom")
    _end_label(a2, 300, 1.3e42, "cobalt's gamma rays", S.INK, dx=0, dy=10, va="bottom")
    _end_label(a2, 440, 1.0e40, "cobalt's positrons", S.INK, dx=0, dy=-10, va="top")
    _end_label(a2, 205, N.deposition(205.0, 35.0), "T0 = 35 d (Type Ia)", S.INK_SECONDARY, dx=0, dy=5,
               ha="left", va="bottom")
    _end_label(a2, 360, N.deposition(360.0, 200.0), "T0 = 200 d (Type II)", S.INK_SECONDARY, dx=0, dy=-4,
               ha="left", va="top")
    S.caption(fig, "Dashed: heat actually deposited when gamma rays start escaping after T0. "
              "Nuclear data: half-lives from NNDC/ENSDF; energies per decay from Nadyozhin (1994): 1.75 MeV in "
              "gamma rays per nickel decay, 3.61 MeV in gamma rays plus 0.12 MeV in positrons per cobalt decay. "
              "One solar mass of nickel-56 releases 1.9e50 erg in all: a fifth of a typical supernova's kinetic energy "
              "(1e51 erg), released slowly enough to be seen.")
    return S.save(fig, FIG / "02_nickel_clock")


# ---------------------------------------------------------------------------
# 3. Type Ia in detail
# ---------------------------------------------------------------------------


def fig_type_ia(r: dict):
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 10.5), sharex="col",
                             gridspec_kw={"height_ratios": [2.4, 1.3, 1.1], "hspace": 0.12})
    for col, sn in enumerate(("SN2011fe", "SN2014J")):
        tab = r[sn]["bolo"]["table"]
        f = r[sn]["fit"]
        af = f["arnett_fit"]
        c = S.SN_COLORS[sn]
        t_max = 360 if sn == "SN2011fe" else 130
        ax = axes[0, col]
        meas = (tab["nir_source"] == "measured").values
        ax.errorbar(tab["phase"], tab["L"], yerr=tab["L_err"], fmt="none", ecolor=c, alpha=0.35, lw=0.8)
        _points(ax, sn, tab["phase"], tab["L"], open_mask=~meas, ms=3.2)
        tt = np.linspace(0.3, t_max, 1200)
        ax.plot(tt, af["m_ni"] * N.deposition(tt, None), color=S.C_MODEL, lw=0.9, ls=(0, (4, 3)))
        ax.plot(tt, af["m_ni"] * N.deposition(tt, af["t0"]), color=S.C_MODEL, lw=0.9)
        ax.plot(tt, N.arnett_light_curve(tt, af["m_ni"], af["tau_m"], af["t0"]), color=c, lw=1.6)
        pk = f["peak"]
        ax.plot([pk["t_peak"]], [pk["L_peak"]], marker="*", ms=11, color=S.INK, mec=S.SURFACE, zorder=6)
        _end_label(ax, pk["t_peak"], pk["L_peak"],
                   f"peak {pk['L_peak'] / 1e43:.2f}e43 erg/s on day {pk['t_peak']:.1f}\n"
                   f"Arnett's rule: {f['arnett_rule']['m_ni']:.2f} solar masses", S.INK, dx=10, dy=6, va="bottom")
        _log_axis(ax, 1e40 if sn == "SN2011fe" else 1e41, 2.5e44)
        ax.set_xlim(0, t_max)
        ax.set_ylabel("Luminosity (erg/s)")
        ax.set_title(f"{LABEL[sn]}: Arnett's model fitted to days 0-{af['t_max_fit']:.0f}")
        txt = (f"M(Ni) = {af['m_ni']:.2f} solar masses\nT0 = {af['t0']:.1f} +/- {af['t0_err']:.1f} d "
               f"(gamma-ray escape)\ndiffusion time = {af['tau_m']:.1f} d")
        ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8, color=S.INK)
        if col == 0:
            handles = [Line2D([], [], color=c, lw=1.6, label="Arnett's model with gamma-ray escape"),
                       Line2D([], [], color=S.C_MODEL, lw=0.9, label="radioactive power deposited (same T0)"),
                       Line2D([], [], color=S.C_MODEL, lw=0.9, ls=(0, (4, 3)), label="radioactive power, all trapped"),
                       Line2D([], [], marker="o", ls="none", mfc=S.SURFACE, mec=S.INK_SECONDARY, ms=4,
                              label="near-infrared share estimated")]
            ax.legend(handles=handles, loc="lower left", fontsize=7.2)

        # Implied nickel mass, epoch by epoch.
        ax = axes[1, col]
        mi = f["m_inst"]
        ok = mi["phase"] <= t_max
        ax.axvspan(0, 40, color=S.EXCLUDED, zorder=0)
        ax.text(20, 0.6, "light still diffusing out", ha="center", va="center", fontsize=7.2, color=S.MUTED,
                rotation=90)
        ax.fill_between(mi["phase"][ok], (mi["m"] - mi["err"])[ok], (mi["m"] + mi["err"])[ok], color=c,
                        alpha=0.18, lw=0)
        ax.plot(mi["phase"][ok], mi["m"][ok], color=c, lw=1.5)
        ax.axhline(af["m_ni"], color=S.INK_SECONDARY, lw=0.9)
        ax.set_ylim(0, 1.2)
        ax.set_ylabel("Nickel implied by\nthe light (solar masses)")
        if sn == "SN2011fe":
            _end_label(ax, 340, 0.36, "after day ~150 light goes missing:\ninto the far infrared, "
                       "or positrons escape\n(Dimitriadis et al. 2017)", S.INK, dx=0, dy=-6, ha="right", va="top")
        else:
            _end_label(ax, 126, 0.98, "rises: M82's bright background and\ndust effects in broad filters",
                       S.INK, dx=-4, dy=8, ha="right", va="bottom")
        _end_label(ax, t_max, af["m_ni"], "Arnett fit", S.INK_SECONDARY, dx=-4, dy=-4, ha="right", va="top")

        # Where the infrared share came from.
        ax = axes[2, col]
        tmpl = r[sn]["bolo"]["nir"]
        last_meas = tab.loc[meas, "phase"].max()
        ax.plot(tab.loc[meas, "phase"], tab.loc[meas, "f_nir_used"], "o", color=c, ms=3, mec=S.SURFACE, mew=0.4)
        est = tab[(~meas) & (tab["phase"] > last_meas)]
        e_err = np.interp(est["phase"], tmpl["phase"], tmpl["err"])
        ax.fill_between(est["phase"], est["f_nir_used"] - e_err, est["f_nir_used"] + e_err, color=c, alpha=0.15, lw=0)
        ax.plot(est["phase"], est["f_nir_used"], color=c, lw=1.2)
        if sn == "SN2011fe":
            a = tmpl["anchor"]
            ax.errorbar([a["phase"]], [a["f"]], yerr=[[0.03], [0.03]], fmt="s", color=S.INK, ms=4)
            _end_label(ax, a["phase"], a["f"], "LBT J band, day 253\n(Shappee et al. 2017)", S.INK, dx=8, dy=12)
            _end_label(ax, 160, 0.27, "nothing observed in the infrared:\nstraight line, wide band", S.INK_SECONDARY,
                       dx=0, dy=8, ha="center", va="bottom")
        else:
            _end_label(ax, 90, 0.3, "SN 2011fe's share after day 48", S.INK_SECONDARY, dx=0, dy=14, ha="center", va="bottom")
        ax.set_ylim(0, 0.6)
        ax.set_ylabel("Share of light\nin the near-infrared")
        ax.set_xlabel("Days since first light")
    S.caption(fig, "Top: bolometric light curves (filled points: near-infrared measured; open: near-infrared share "
              "adopted, bottom row). Middle: the nickel mass each epoch's luminosity implies for the fitted T0 with "
              "positrons trapped; a flat line means the light is all accounted for. SN 2014J is dereddened with "
              "E(B-V) = 1.37, R_V = 1.4 (Amanullah et al. 2014).")
    return S.save(fig, FIG / "03_type_ia")


# ---------------------------------------------------------------------------
# 4. SN 2014J's gamma rays
# ---------------------------------------------------------------------------


def fig_gamma(r: dict):
    g = r["gamma"]
    d_cm = 10 ** (0.2 * C.SUPERNOVAE["SN2014J"]["mu"] + 1.0) * C.PC_CM
    fig = plt.figure(figsize=(11.5, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.38, wspace=0.25)
    t = np.linspace(5, 180, 200)
    dp = g["diehl"]
    m_u = g["uniform"]["combined"][0]
    m_e = g["exponential"]["combined"][0]
    ej_u = G.Ejecta("uniform", m_ni_msun=m_u)
    ej_e = G.Ejecta("exponential", m_ni_msun=m_e)
    for k, line in enumerate((847, 1238)):
        ax = fig.add_subplot(gs[0, k])
        ax.plot(t, G.line_flux(t, m_u, line, d_cm, None) * 1e4, color=S.C_MODEL, lw=0.9, ls=(0, (4, 3)))
        ax.plot(t, G.line_flux(t, m_u, line, d_cm, ej_u) * 1e4, color=S.VIOLET, lw=1.6)
        ax.plot(t, G.line_flux(t, m_e, line, d_cm, ej_e) * 1e4, color=S.VIOLET, lw=1.0, alpha=0.5)
        xerr = [dp["t_mid"] - dp["t_lo"], dp["t_hi"] - dp["t_mid"]]
        ax.errorbar(dp["t_mid"], dp[line] * 1e4, xerr=xerr, yerr=dp[f"{line}_err"] * 1e4, fmt="o",
                    color=S.INK_SECONDARY, ms=4, lw=0.9, capsize=0)
        lo, hi = G.CHURAZOV2014["t_range"]
        f, e = G.CHURAZOV2014[line]
        ax.fill_between([lo, hi], [f - e] * 2, [f + e] * 2, color=S.INK, alpha=0.10, lw=0)
        ax.plot([lo, hi], [f, f], color=S.INK, lw=2.0)
        ax.set_xlim(0, 180)
        ax.set_ylim(0, 8.0 if line == 847 else 6.0)
        ax.set_xlabel("Days since explosion")
        ax.set_ylabel("Line flux (1e-4 photons per cm² per s)")
        ax.set_title(f"Cobalt-56 gamma rays at {C.CO56_LINES[line]['energy_kev']:.0f} keV")
        _end_label(ax, 150, G.line_flux(150.0, m_u, line, d_cm, None)[0] * 1e4, "if the ejecta\nwere transparent",
                   S.INK_SECONDARY, dx=0, dy=14, ha="center", va="bottom")
        if k == 0:
            handles = [Line2D([], [], color=S.INK, lw=2.0, label="INTEGRAL, days 50-100 (Churazov et al. 2014)"),
                       Line2D([], [], color=S.INK_SECONDARY, marker="o", lw=0.9, label="INTEGRAL epochs (Diehl et al. 2015)"),
                       Line2D([], [], color=S.VIOLET, lw=1.6, label=f"model: {m_u:.2f} Msun, uniform ejecta"),
                       Line2D([], [], color=S.VIOLET, lw=1.0, alpha=0.5, label=f"model: {m_e:.2f} Msun, W7-like ejecta")]
            ax.legend(handles=handles, loc="upper right", fontsize=7)

    # Forest plot of every nickel estimate for SN 2014J.
    ax = fig.add_subplot(gs[1, :])
    s = r["summary"]
    v = s["SN2014J"]
    rows = [
        ("gamma rays, transparent ejecta (lower limit)", s["SN2014J_gamma"]["transparent"]["m_ni"], S.MUTED, "o"),
        ("gamma rays, uniform ejecta", s["SN2014J_gamma"]["uniform"]["m_ni"], S.VIOLET, "o"),
        ("gamma rays, W7-like ejecta", s["SN2014J_gamma"]["exponential"]["m_ni"], S.VIOLET, "o"),
        ("light: Arnett's model, days 0-128", v["m_ni_arnett_fit"], S.VIOLET, "s"),
        ("light: Arnett's model, days 0-60", v["m_ni_peak_only_t0_from_2011fe"], S.VIOLET, "s"),
        ("light: Arnett's rule, dust R_V = 1.4", v["m_ni_arnett_rule"], S.VIOLET, "s"),
        ("light: Arnett's rule, dust R_V = 1.6", {"value": v["m_ni_arnett_rule_alt_dust"], "err": None}, S.VIOLET, "s"),
    ]
    for lit in LIT.NICKEL["SN2014J"]:
        rows.append((lit["ref"].split(" (")[0] + " (gamma rays)", {"value": lit["m_ni"], "err": lit["err"]}, S.INK, "D"))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, val, col, mk) in zip(y, rows):
        if val["err"]:
            ax.errorbar(val["value"], yi, xerr=val["err"], fmt=mk, color=col, ms=5.5, lw=1.2, capsize=0)
        else:
            ax.plot(val["value"], yi, mk, color=col, ms=5.5)
        if "lower limit" in lab:
            ax.annotate("", xy=(val["value"] + 0.25, yi), xytext=(val["value"], yi),
                        arrowprops={"arrowstyle": "->", "color": S.MUTED, "lw": 1.0})
    ax.set_yticks(y)
    ax.set_yticklabels([row[0] for row in rows], fontsize=8)
    ax.set_xlim(0, 1.2)
    ax.set_xlabel("Nickel-56 mass (solar masses)")
    ax.set_title("SN 2014J, weighed by its gamma rays and by its light (squares: this project's optical; diamonds: published)")
    ax.grid(axis="y", visible=False)
    S.caption(fig, "Gamma-ray fluxes: Churazov et al. 2014 (Extended Data Table 1) and Diehl et al. 2015 (Table 1; "
              "errors for epochs other than the third scaled by epoch length, as only that one is quoted). Models: this "
              "project's escape probability for line photons through homologous ejecta (1.38 Msun, 1.2e51 erg). "
              "Distance 3.53 +/- 0.26 Mpc: +/-15% on every mass, included in the errors.")
    return S.save(fig, FIG / "04_gamma_rays")


# ---------------------------------------------------------------------------
# 5. Type II tails
# ---------------------------------------------------------------------------


def fig_type_ii(r: dict):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6), sharex="col",
                             gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.1})
    # SN 2023ixf
    tab = r["SN2023ixf"]["bolo"]["table"]
    fit = r["SN2023ixf"]["fit"]["tail_measured_nir"]
    ax = axes[0, 0]
    meas = tab["nir_source"] == "measured"
    ok = tab["L_err"] < 0.3 * tab["L"]  # the first days' hot-UV extrapolation is very uncertain
    ax.errorbar(tab["phase"][ok], tab["L"][ok], yerr=tab["L_err"][ok], fmt="none", ecolor=S.AQUA, alpha=0.3, lw=0.8)
    _points(ax, "SN2023ixf", tab["phase"], tab["L"], open_mask=~meas, ms=3.4)
    tt = np.linspace(90, 365, 600)
    ax.plot(tt, fit["m_ni"] * N.deposition(tt, None), color=S.C_MODEL, lw=0.9, ls=(0, (4, 3)))
    ax.plot(tt, fit["m_ni"] * N.deposition(tt, fit["t0"]), color=S.AQUA, lw=1.5)
    ax.axvspan(*fit["window"], color=S.EXCLUDED, zorder=0)
    _log_axis(ax, 2e40, 3e43)
    ax.set_xlim(0, 365)
    ax.set_ylabel("Luminosity (erg/s)")
    ax.set_title("SN 2023ixf: plateau, drop, radioactive tail")
    S.annotate(ax, "plateau: hydrogen recombining\nas the ejecta cool", (45, 2.0e42), (120, 6e42))
    S.annotate(ax, "drop onto the\nradioactive tail", (86, 9e41), (150, 1.6e42))
    fy = r["SN2023ixf"]["fit"]["tail_full_year"]
    ax.text(0.03, 0.03, f"tail, days {fit['window'][0]:.0f}-{fit['window'][1]:.0f} (shaded):\n"
            f"M(Ni) = {fit['m_ni']:.3f} solar masses\nT0 = {fit['t0']:.0f} +/- {fit['t0_err']:.0f} d\n"
            f"(to day {fy['window'][1]:.0f}, open points included:\n{fy['m_ni']:.3f} solar masses, T0 = {fy['t0']:.0f} d)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=S.INK)

    # SN 1987A
    sn87 = r["SN1987A"]["bolo"]
    ax = axes[0, 1]
    pub = sn87["published"]
    late = sn87["published_late"]
    late = late[late["phase"] > pub["phase"].max()]
    ours = sn87["ours"]
    ax.plot(ours["phase"], ours["L_uvoi"], color=S.ORANGE, lw=1.0, alpha=0.6)
    _points(ax, "SN1987A", pd.concat([pub["phase"], late["phase"]]), pd.concat([pub["L"], late["L"]]), ms=3.6)
    f87 = r["SN1987A"]["fit"]["tail_with_escape"]
    ft = r["SN1987A"]["fit"]["tail_full_trapping"]
    tt = np.linspace(60, 1400, 1000)
    ax.plot(tt, ft["m_ni"] * N.deposition(tt, None), color=S.C_MODEL, lw=0.9, ls=(0, (4, 3)))
    ax.plot(tt, f87["m_ni"] * N.deposition(tt, f87["t0"]), color=S.ORANGE, lw=1.5)
    ax.axvspan(*f87["window"], color=S.EXCLUDED, zorder=0)
    _log_axis(ax, 2e37, 3e42)
    ax.set_xlim(0, 1400)
    ax.set_title("SN 1987A: the textbook tail")
    ax.text(0.03, 0.03, f"days {ft['window'][0]:.0f}-{ft['window'][1]:.0f}, all trapped: {ft['m_ni']:.4f} Msun\n"
            f"days {f87['window'][0]:.0f}-{f87['window'][1]:.0f} (shaded): {f87['m_ni']:.4f} Msun,\n"
            f"T0 = {f87['t0']:.0f} d", transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=S.INK)
    _end_label(ax, 200, 3.4e41, "UVOIR, Suntzeff et al. 1991-92,\nBouchet et al. 1991", S.INK, dx=10, dy=6, va="bottom")
    _end_label(ax, 560, 3.0e39, "U to I only\n(this project)", S.INK, dx=-8, dy=-6, ha="right", va="top")
    _end_label(ax, 1050, 1.6e38, "57Co, 44Ti and dust\ntake over", S.INK, dx=4, dy=8)

    # Bottom: the light outside the optical.
    ax = axes[1, 0]
    nir = r["SN2023ixf"]["bolo"]["nir"]["measured"]
    ax.plot(nir["phase"], nir["f_nir"], "o", color=S.AQUA, ms=3, mec=S.SURFACE, mew=0.4)
    ax.set_ylim(0, 0.75)
    ax.set_ylabel("Share of light\nbeyond 0.9 micron")
    ax.set_xlabel("Days since explosion")
    _end_label(ax, 100, 0.52, "half the tail's light is infrared", S.INK, dx=0, dy=10, ha="center", va="bottom")
    ax = axes[1, 1]
    ok = ours["f_outside_ubvri"].notna()
    ax.plot(ours["phase"][ok], ours["f_outside_ubvri"][ok], color=S.ORANGE, lw=1.4)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Share of light\noutside U-I")
    ax.set_xlabel("Days since explosion")
    _end_label(ax, 780, 0.75, "dust forms, light\nmoves to the infrared", S.INK, dx=-10, dy=0, ha="right")
    S.caption(fig, "SN 2023ixf: Hsu et al. 2025 (UBVgriz), Singh et al. 2024 (Swift UV, BVRI, JHKs), integrated here; "
              "open points use the near-infrared share of days 120-145. SN 1987A: published UVOIR rescaled to the "
              "LMC's 49.59 kpc; the thin line is this project's own U-to-I integral of the SAAO/CTIO photometry. "
              "Dashed: every gamma ray trapped.")
    return S.save(fig, FIG / "05_type_ii")


# ---------------------------------------------------------------------------
# 6. Arnett's rule for hundreds of ZTF Type Ia
# ---------------------------------------------------------------------------


def fig_ztf(r: dict):
    from matplotlib.ticker import FixedLocator, NullLocator

    z = r["ztf"]["good"]
    low = z["low_dust"] & ~z["outlier"]
    out = z["outlier"]
    rest = ~z["low_dust"]
    s = r["summary"]["ztf"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), gridspec_kw={"wspace": 0.3})
    markers = {"normal": "o", "91T-like": "^", "91bg-like": "v"}

    ax = axes[0]
    for sel, kw, zo in ((rest, {"color": S.C_REST}, 2), (low, {"color": S.BLUE}, 3),
                        (out, {"mfc": S.SURFACE, "mec": S.BLUE, "mew": 1.0}, 4)):
        for typ, mk in markers.items():
            d = z[sel & (z["type"] == typ)]
            style = {"mec": S.SURFACE, "mew": 0.5, **kw}
            ax.plot(d["dm15_g"], d["M_g"], mk, ms=4.4, zorder=zo, ls="none", **style)
    ph = r["ztf"]["phillips"]
    xx = np.linspace(0.45, 1.8, 10)
    ax.plot(xx, ph["intercept"] + ph["slope"] * xx, color=S.INK_SECONDARY, lw=1.0)
    ax.invert_yaxis()
    ax.set_xlabel("Fading in 15 days after peak, Δm15(g) (mag)")
    ax.set_ylabel("Peak absolute magnitude, ZTF g")
    ax.set_title("Brighter supernovae fade more slowly")
    handles = [Line2D([], [], marker="o", ls="none", color=S.BLUE, ms=5,
                      label=f"blue at peak, little dust ({int(low.sum())})"),
               Line2D([], [], marker="o", ls="none", color=S.C_REST, ms=5, label=f"redder: host dust ({int(rest.sum())})"),
               Line2D([], [], marker="o", ls="none", mfc=S.SURFACE, mec=S.BLUE, ms=5,
                      label=f"> 3 sigma off the relation ({int(out.sum())})"),
               Line2D([], [], marker="^", ls="none", color=S.INK_SECONDARY, ms=5, label="91T-like (overluminous)")]
    ax.legend(handles=handles, loc="lower left", fontsize=7)

    ax = axes[1]
    b = z[low]
    for typ, mk in markers.items():
        d = b[b["type"] == typ]
        ax.errorbar(d["dm15_g"], d["m_ni"], yerr=d["m_ni_err"], fmt=mk, color=S.BLUE, ms=4.4, mec=S.SURFACE,
                    mew=0.5, elinewidth=0.6, alpha=0.9)
    ax.set_yscale("log")
    ax.set_ylim(0.15, 1.6)
    ax.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5, 0.7, 1.0, 1.5]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels(["0.2", "0.3", "0.5", "0.7", "1", "1.5"])
    ax.set_xlabel("Δm15(g) (mag)")
    ax.set_ylabel("Nickel-56 from Arnett's rule (solar masses)")
    ax.set_title("More nickel, slower fading")
    ax.text(0.97, 0.97, f"correlation r = {s['corr_dm15_logmni']:.2f}", transform=ax.transAxes, ha="right",
            va="top", fontsize=8, color=S.INK)

    ax = axes[2]
    bins = np.linspace(0.0, 1.4, 29)
    ax.hist(b["m_ni"], bins=bins, color=S.BLUE, alpha=0.85, rwidth=0.82)
    ymax = ax.get_ylim()[1]
    rs = r["summary"]
    for sn, frac, side in (("SN2011fe", 0.95, "right"), ("SN2014J", 0.83, "left")):
        v = rs[sn]["m_ni_arnett_rule"]["value"]
        ax.axvline(v, color=S.SN_COLORS[sn], lw=1.6)
        _end_label(ax, v, ymax * frac, LABEL[sn], S.INK, dx=-4 if side == "right" else 4, ha=side)
    from matplotlib.ticker import MaxNLocator

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Nickel-56 mass (solar masses)")
    ax.set_ylabel("Number of supernovae")
    ax.set_title(f"{int(low.sum())} low-dust Type Ia: median {s['m_ni_median']:.2f} Msun")
    S.caption(fig, f"ZTF Bright Transient Survey (z = 0.02-0.04), alert photometry from ALeRCE: {len(z)} of "
              f"{len(r['ztf']['all'])} candidates pass the sampling cuts. Peak luminosity = SN 2011fe's bolometric peak "
              "x the ratio of g+r fluxes at g-band peak (K-corrected with SN 2011fe's spectrum); rise time from a "
              "fireball fit to the first detections. H0 = 73 km/s/Mpc. Arnett's rule shown for SN 2011fe and 2014J.")
    return S.save(fig, FIG / "06_ztf_arnett")


# ---------------------------------------------------------------------------
# 7. The census
# ---------------------------------------------------------------------------


def fig_census(r: dict):
    s = r["summary"]
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    groups = [
        ("SN2011fe", [("Arnett's model fit", s["SN2011fe"]["m_ni_arnett_fit"]),
                      ("Arnett's rule", s["SN2011fe"]["m_ni_arnett_rule"])]),
        ("SN2014J", [("Arnett's model fit", s["SN2014J"]["m_ni_arnett_fit"]),
                     ("Arnett's rule", s["SN2014J"]["m_ni_arnett_rule"]),
                     ("gamma rays, uniform ejecta", s["SN2014J_gamma"]["uniform"]["m_ni"])]),
        ("SN2023ixf", [("radioactive tail", s["SN2023ixf"]["tail_measured_nir"]["m_ni"])]),
        ("SN1987A", [("radioactive tail", s["SN1987A"]["tail_with_escape"]["m_ni"])]),
    ]
    rows = []
    for sn, ours in groups:
        for lab, v in ours:
            rows.append((sn, f"{LABEL[sn]}: {lab}", v["value"], v["err"], "ours"))
        for lit in LIT.NICKEL[sn][:2]:
            rows.append((sn, f"{LABEL[sn]}: {lit['ref'].split(' (')[0]}", lit["m_ni"], lit["err"], "published"))
        rows.append(None)
    rows = rows[:-1]
    n = len(rows)
    ytick, ylab = [], []
    for i, row in enumerate(rows):
        if row is None:
            continue
        sn, lab, v, e, kind = row
        yi = n - i
        col = S.SN_COLORS[sn] if kind == "ours" else S.C_LIT
        mk = S.SN_MARKERS[sn] if kind == "ours" else "D"
        if e:
            ax.errorbar(v, yi, xerr=e, fmt=mk, color=col, ms=5.5, lw=1.2, capsize=0)
        else:
            ax.plot(v, yi, mk, color=col, ms=5.5)
        ytick.append(yi)
        ylab.append(lab)
    z = r["ztf"]["good"]
    zb = z[z["low_dust"] & ~z["outlier"]]
    q = np.percentile(zb["m_ni"], [16, 50, 84])
    yz = n + 1.5
    ax.plot([q[0], q[2]], [yz, yz], color=S.BLUE, lw=6, alpha=0.35, solid_capstyle="butt")
    ax.plot(q[1], yz, "|", color=S.BLUE, ms=14, mew=2.2)
    ytick.append(yz)
    ylab.append(f"{len(zb)} ZTF Type Ia (median, 68% range)")
    ax.set_yticks(ytick)
    ax.set_yticklabels(ylab, fontsize=7.8)
    ax.set_xscale("log")
    ax.set_xlim(0.02, 2)
    ax.set_xticks([0.03, 0.05, 0.1, 0.2, 0.5, 1.0])
    ax.set_xticklabels(["0.03", "0.05", "0.1", "0.2", "0.5", "1"])
    ax.set_xlabel("Nickel-56 made in the explosion (solar masses)")
    ax.set_title("Exploding white dwarfs make ten times more nickel than collapsing giants")
    ax.grid(axis="y", visible=False)
    handles = [Line2D([], [], marker="o", ls="none", color=S.INK_SECONDARY, ms=5, label="this project"),
               Line2D([], [], marker="D", ls="none", color=S.C_LIT, ms=5, label="published")]
    ax.legend(handles=handles, loc="lower right", fontsize=7.5)
    S.caption(fig, "Errors include distance; published values assume slightly different distances (see README). "
                   "Type Ia: thermonuclear explosions of white dwarfs. Type II: core collapse of massive stars.")
    return S.save(fig, FIG / "07_nickel_census")


def make_all(r: dict) -> list:
    out = []
    for fn in (fig_radioactive_glow, fig_nickel_clock, fig_type_ia, fig_gamma, fig_type_ii, fig_ztf, fig_census):
        out += fn(r)
    return out
