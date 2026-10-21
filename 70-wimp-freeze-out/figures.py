"""All figures, from data/results.json and data/derived/figdata.npz.

    python figures.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FixedLocator, LogLocator, NullFormatter  # noqa: E402

import config as C  # noqa: E402
import plotstyle as ps  # noqa: E402

F = None   # figdata
R = None   # results


def _load():
    global F, R
    F = dict(np.load(C.DERIVED_DIR / "figdata.npz"))
    R = json.loads(C.RESULTS_PATH.read_text())


def _fmt_sv(v):
    return f"{v / 1e-26:.2f}×10⁻²⁶"


def _sci(v, digits=1):
    e = int(math.floor(math.log10(abs(v))))
    m = v / 10 ** e
    sup = str(e).replace("-", "⁻").translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m:.{digits}f}×10{sup}"


def _top_temperature_axis(ax, m, ticks_T):
    """A twin x axis giving T = m/x for a log x axis in x = m/T."""
    top = ax.twiny()
    top.set_xscale("log")
    top.set_xlim(ax.get_xlim())
    top.xaxis.set_major_locator(FixedLocator([m / T for T in ticks_T]))
    top.xaxis.set_minor_locator(FixedLocator([]))
    top.set_xticklabels([(f"{T:g} GeV" if T >= 1 else f"{T * 1e3:g} MeV") for T in ticks_T])
    top.tick_params(axis="x", colors=ps.MUTED, labelsize=7.5)
    top.grid(False)
    for s in ("top", "right"):
        top.spines[s].set_visible(False)
    return top


# --- 1. freeze-out ------------------------------------------------------------------------------------


def fig_freeze_out(out: Path):
    ex = R["examples"]
    x = F["yx_x"]
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.loglog(x, F["yx_eq"], color=ps.C_MODEL, lw=1.3, ls="--")
    ps.annotate(ax, "equilibrium: Y_eq ∝ x^3/2 e^−x\n(annihilation keeps up)", (12.0, float(np.interp(12, x, F["yx_eq"]))),
                (1.6, 3e-8), color=ps.C_MODEL)
    greys = [ps.lighten(ps.INK, a) for a in (0.62, 0.5, 0.38, 0.0, 0.38, 0.5)]
    labels = []
    for k, row in enumerate(ex["rows"]):
        thermal = k == 3
        ax.loglog(x, F[f"yx_{k}"], color=ps.INK if thermal else greys[k], lw=2.2 if thermal else 1.1, zorder=3 if thermal else 2)
        lab = (f"⟨σv⟩ = {_fmt_sv(row['sigv'])} cm³/s → Ω h² = {row['omega_h2']:.3f}" if thermal else
               f"{_sci(row['sigv'], 0)} cm³/s → Ω h² = {row['omega_h2']:.3g}")
        labels.append([math.log10(F[f"yx_{k}"][-1]), lab, thermal])
    # spread the end labels at least 0.62 decades apart (one text line on this axis), keeping order
    labels.sort(key=lambda r: r[0])
    pos = [r[0] for r in labels]
    for _ in range(200):
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < 0.62:
                mid = 0.5 * (pos[i] + pos[i - 1])
                pos[i - 1], pos[i] = mid - 0.31, mid + 0.31
                moved = True
        if not moved:
            break
    for (y0, lab, thermal), y1 in zip(labels, pos):
        ax.annotate(lab, xy=(1e3, 10 ** y0), xytext=(1.25e3, 10 ** y1), fontsize=7.2, va="center",
                    color=ps.INK if thermal else ps.INK_SECONDARY, fontweight="bold" if thermal else "normal",
                    arrowprops={"arrowstyle": "-", "color": ps.AXIS, "lw": 0.6, "shrinkA": 0, "shrinkB": 1},
                    annotation_clip=False)
    ax.axhline(ex["Y_for_omega_0120"], color=ps.ORANGE, lw=0.9, ls=":")
    ax.text(1.15, ex["Y_for_omega_0120"] / 3.0, "the Y that makes Planck's\nΩ_c h² = 0.120 at 100 GeV",
            fontsize=7.5, color=ps.ORANGE, va="top")
    xf = ex["rows"][3]["x_f"]
    ax.axvline(xf, color=ps.AXIS, lw=0.8)
    ax.text(xf * 1.08, 3e-3, f"freeze-out\nx_f = {xf:.1f}\nT ≈ {100 / xf:.1f} GeV,\nt ≈ {_sci(ex['t_at_x_f_s'], 1)} s",
            fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    ax.set_xlim(1, 1e3)
    ax.set_ylim(1e-20, 1e-1)
    ax.set_xlabel("x = m / T   (time runs to the right)")
    ax.set_ylabel("Y = n / s   (WIMPs per unit entropy)")
    ax.set_title("Freeze-out of a 100 GeV WIMP: weaker annihilation leaves more behind", pad=26)
    _top_temperature_axis(ax, 100.0, [100, 10, 1, 0.1])
    ps.caption(fig, "Our Boltzmann solver (implicit Radau, Saikawa & Shirai 2020 g*(T), Majorana fermion, s-wave) — "
               "the plot Kolb & Turner made famous. Each grey curve is a 10× different annihilation cross section; "
               "the relic density scales almost exactly as 1/⟨σv⟩.", y=-0.01)
    return ps.save(fig, out / "01_freeze_out")


# --- 2. the thermal cross section ------------------------------------------------------------------


def fig_thermal(out: Path):
    m = F["th_m"]
    h = R["thermal"]["headline"]
    fig = plt.figure(figsize=(7.6, 8.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[3.2, 1.1, 1.4], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    ax.fill_between(m, F["th_eos_ss2020_lo"] / 1e-26, F["th_eos_ss2020_hi"] / 1e-26, color=ps.lighten(ps.INK, 0.85),
                    lw=0, label="equation-of-state uncertainty (Saikawa & Shirai)")
    ax.fill_between(m, F["th_om_hi"] / 1e-26, F["th_om_lo"] / 1e-26, color=ps.lighten(ps.C_MODEL, 0.6), lw=0,
                    label="Planck Ω_c h² = 0.1200 ± 0.0012")
    ax.semilogx(m, F["th_main"] / 1e-26, color=ps.INK, lw=2.0, label="this project (Ω_c h² = 0.120)")
    ax.semilogx(F["ss2020_m"], F["ss2020_sv"] / 1e-26, color=ps.ORANGE, lw=1.0, ls=(0, (4, 2)),
                label="Saikawa & Shirai 2020 (their table)")
    ax.semilogx(F["sdb12_m"], F["sdb12_sv"] / 1e-26, color=ps.AQUA, lw=1.2, marker="o", ms=2.5,
                label="Steigman, Dasgupta & Beacom 2012 (Ω h² = 0.11)")
    ax.axhline(3.0, color=ps.C_MODEL, lw=0.9, ls=":")
    ax.text(0.012, 3.05, "the 'canonical' 3×10⁻²⁶ cm³/s", fontsize=7.5, color=ps.C_MODEL, va="bottom")
    ax.set_ylabel("⟨σv⟩ needed  (10⁻²⁶ cm³/s)")
    ax.set_ylim(1.6, 5.4)
    ax.set_xlim(m[0], m[-1])
    ax.legend(loc="upper right", fontsize=7.3)
    ps.annotate(ax, f"{h['sigv_100GeV'] / 1e-26:.2f} at 100 GeV", (100, h["sigv_100GeV"] / 1e-26), (250, 1.85))
    ps.annotate(ax, "quarks and gluons bind into hadrons\n(QCD crossover at T ≈ 156 MeV)", (3.5, 3.3), (0.03, 2.05))
    ps.annotate(ax, "freeze-out during e⁺e⁻ annihilation\nand neutrino decoupling", (0.015, 3.94), (0.012, 3.35))
    ax.tick_params(labelbottom=False)
    ax.set_title("The annihilation cross section that makes all of the dark matter", pad=22)

    ax2 = fig.add_subplot(gs[1], sharex=ax)
    ax2.axhspan(0.998, 1.002, color=ps.WASH, lw=0)
    ax2.semilogx(F["val_ss_m"], F["val_ss_ratio"], color=ps.ORANGE, marker="o", ms=2.5, lw=0.8,
                 label="ours ÷ Saikawa & Shirai 2020")
    ax2.semilogx(F["sdb12_m"], F["sdb12_ratio"], color=ps.AQUA, marker="s", ms=2.5, lw=0.8,
                 label="ours in SDB12's set-up ÷ SDB12")
    ax2.set_ylim(0.88, 1.04)
    ax2.set_ylabel("ratio")
    ax2.legend(loc="lower right", fontsize=7, ncol=1)
    ax2.tick_params(labelbottom=False)
    v = R["thermal"]["validation_ss2020"]["majorana_s"]
    ax2.text(2e4, 1.025, f"median {v['median_ratio']:.4f}", fontsize=7, color=ps.ORANGE, ha="center")

    ax3 = fig.add_subplot(gs[2], sharex=ax)
    xf = F["th_xf"]
    Tf = xf[:, 0] / xf[:, 1]
    T = F["eos_T"]
    g = np.interp(np.log(Tf), np.log(T), F["eos_ss2020_grho"])
    ax3.semilogx(xf[:, 0], g, color=ps.INK, lw=1.4)
    ax3.set_ylabel("g*ρ at freeze-out")
    ax3.set_xlabel("WIMP mass m (GeV)")
    ax4 = ax3.twinx()
    ax4.semilogx(xf[:, 0], xf[:, 1], color=ps.MUTED, lw=1.0, ls="--")
    ax4.set_ylabel("x_f = m / T_f", color=ps.MUTED)
    ax4.tick_params(axis="y", colors=ps.MUTED)
    ax4.grid(False)
    ax4.spines["right"].set_visible(True)
    ax4.spines["right"].set_color(ps.AXIS)
    ax3.text(2e-2, 75, "freeze-out happens at T ≈ m/20, so the mass axis is also\n"
             "a thermometer: the curve above maps the particles in the plasma", fontsize=7.2, color=ps.INK_SECONDARY)
    ps.caption(fig, "Our solver with Saikawa & Shirai's (2020) g*(T) reproduces their published relic cross sections "
               f"to {100 * (v['median_ratio'] - 1):.2f}% at every mass (a constant offset); in SDB12's own set-up it "
               "reproduces their curve to <1% away from the QCD crossover, where their older equation of state differs.",
               y=0.02)
    return ps.save(fig, out / "02_thermal_cross_section")


# --- 3. what changes the number ------------------------------------------------------------------------


def fig_conventions(out: Path):
    m = F["th_m"]
    main = F["th_main"]
    fig, axs = plt.subplots(2, 2, figsize=(9.4, 7.0))
    ax = axs[0, 0]
    for key, col, lab in (("real_scalar", ps.MUTED, "real scalar (g = 1)"), ("vector", ps.AQUA, "real vector (g = 3)"),
                          ("dirac", ps.ORANGE, "Dirac fermion (g = 4), ⟨σv⟩ per n_DM²")):
        ax.semilogx(m, F[f"th_type_{key}"] / main, color=col, lw=1.4, label=lab)
    ax.axhline(1, color=ps.INK, lw=1.2, label="Majorana fermion or complex scalar (g = 2)")
    lo = float(min(np.nanmin(F[f"th_type_{k}"] / main) for k in ("real_scalar", "vector", "dirac")))
    ax.set_ylim(min(0.95, lo - 0.02), 1.075)
    ax.set_title("(a) Who is annihilating: internal states")
    ax.set_ylabel("⟨σv⟩ ÷ Majorana value")
    ax.legend(fontsize=7, loc="upper right")
    ax.text(0.012, min(0.95, lo - 0.02) + 0.004, "a Dirac WIMP's physical χχ̄ cross section is twice the orange value",
            fontsize=7, color=ps.ORANGE, va="bottom")

    ax = axs[0, 1]
    variants = [("lms2015", ps.AQUA, "Laine & Meyer 2015"), ("p2", ps.INK_SECONDARY, "P2's own g*(T) (68-)"),
                ("p2_ideal", ps.C_MODEL, "ideal gas, sharp QCD switch (P2)"), ("ideal", ps.C_MODEL, "ideal gas (fallback)")]
    drawn = set()
    for key, col, lab in variants:
        if f"th_eos_{key}" not in F or (key == "ideal" and "th_eos_p2_ideal" in F):
            continue
        ls = "--" if key.startswith("constant") and "10.75" in key else ("-." if key.startswith("constant") else "-")
        ax.semilogx(m, F[f"th_eos_{key}"] / main, color=col, lw=1.3, ls=ls, label=lab)
        drawn.add(key)
    ax.fill_between(m, F["th_eos_ss2020_lo"] / main, F["th_eos_ss2020_hi"] / main, color=ps.lighten(ps.INK, 0.85), lw=0,
                    label="Saikawa & Shirai's uncertainty")
    ax.axhline(1, color=ps.INK, lw=1.0)
    ax.set_ylim(0.7, 1.12)
    ax.set_title("(b) The plasma it freezes out of: g*(T)")
    ax.set_ylabel("⟨σv⟩ ÷ value with Saikawa & Shirai")
    ax.legend(fontsize=6.8, loc="lower right", ncol=1)
    er = R["thermal"]["eos_relative"]
    if "constant_106.75" in er and "constant_10.75" in er:
        r106 = F["th_eos_constant_106.75"] / main
        r10 = F["th_eos_constant_10.75"] / main
        ax.text(0.012, 1.11, f"ignoring the staircase does far worse: g* = 106.75 throughout gives\n"
                f"×{r106.min():.2f}–{r106.max():.2f}, g* = 10.75 throughout ×{r10.min():.2f}–{r10.max():.2f} (off this scale)",
                fontsize=6.9, color=ps.INK_SECONDARY, va="top")

    ax = axs[1, 0]
    ax.loglog(m, F["th_pwave"], color=ps.INK, lw=1.6, label="b, Majorana: ⟨σv⟩ = b ⟨v²⟩ = 6b/x")
    ax.loglog(m, F["th_main"], color=ps.MUTED, lw=1.0, ls="--", label="s-wave a, for comparison")
    xf = F["th_xf"]
    b100 = R["thermal"]["headline"]["pwave_b_100GeV"]
    ax.set_title("(c) Velocity-suppressed (p-wave) annihilation")
    ax.set_ylabel("cm³/s")
    ax.set_xlabel("WIMP mass (GeV)")
    ax.legend(fontsize=7.2, loc="lower right", bbox_to_anchor=(1.0, 0.12))
    ax.text(0.012, 1.15e-25, f"b ≈ {_sci(b100, 1)} cm³/s at 100 GeV. In today's halo\n(v ~ 10⁻³ c) it would annihilate "
            "~10⁵–10⁶× more slowly", fontsize=7, color=ps.INK_SECONDARY, va="top")

    ax = axs[1, 1]
    dec = R["thermal"]["decomposition_100GeV"]
    vals = [d["sigv_100GeV"] / 1e-26 for d in dec]
    labels = ["SDB12 set-up (Ω h² = 0.11, g*s = g*ρ,\nstop at x = 1000, Y_eq ∝ x^3/2 e^−x)",
              "+ Planck Ω h² = 0.120", "+ integrate to T = 1 eV", "+ separate g*s, g*ρ (Laine et al.)",
              "+ exact Y_eq (Bessel K₂)", "+ Saikawa & Shirai g*(T)"]
    y = np.arange(len(vals))[::-1]
    cols = [ps.AQUA] + [ps.lighten(ps.INK, 0.45)] * (len(vals) - 2) + [ps.INK]
    ax.plot(vals, y, color=ps.AXIS, lw=1.0, zorder=1)
    ax.scatter(vals, y, c=cols, s=46, zorder=2)
    for yy, v in zip(y, vals):
        ax.text(v + 0.012, yy, f"{v:.3f}", va="center", fontsize=7.4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(1.9, 2.5)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("⟨σv⟩ at 100 GeV (10⁻²⁶ cm³/s)")
    ax.set_title("(d) From SDB12's 2.2 to today's 2.1, one change at a time")
    for a in axs[0]:
        a.set_xlabel("WIMP mass (GeV)")
    fig.tight_layout()
    ps.caption(fig, "All curves are ours. The number that matters most is the one we measure: Ω_c h². The particle's spin "
               "and the plasma's equation of state move the answer by a few per cent; only a badly wrong g*(T) moves it more.",
               y=-0.005)
    return ps.save(fig, out / "03_what_changes_the_number")


# --- 4. analytic vs numerical ----------------------------------------------------------------------------


def fig_analytic(out: Path):
    a = R["analytic"]
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ax = axs[0]
    ax.semilogx(F["an_m"], F["an_kt"] / F["an_num"], color=ps.C_MODEL, lw=1.5, marker="o", ms=3,
                label="Kolb & Turner (1990): x_f from ln(0.038 g M_Pl m⟨σv⟩/√g*)")
    ax.semilogx(F["an_m"], F["an_im"] / F["an_num"], color=ps.ORANGE, lw=1.5, marker="s", ms=3,
                label="x_f + exact integral of λ(x) after it (Gondolo & Gelmini 1991)")
    ax.axhline(1, color=ps.INK, lw=1.0)
    lo = float(min(np.min(F["an_kt"] / F["an_num"]), np.min(F["an_im"] / F["an_num"])))
    ax.set_ylim(lo - 0.03, 1.04)
    ax.set_xlabel("WIMP mass (GeV)")
    ax.set_ylabel("thermal ⟨σv⟩: analytic ÷ numerical")
    ax.set_title("(a) The analytic estimates")
    ax.legend(fontsize=7, loc="lower left")
    ax = axs[1]
    rows = F["an_rows"]
    for sv, mk, col in ((2e-27, "v", ps.MUTED), (2.2e-26, "o", ps.INK), (2e-25, "^", ps.MUTED)):
        sel = np.isclose(rows[:, 1], sv, rtol=1e-6, atol=0.0)
        ax.semilogx(rows[sel, 0], rows[sel, 5], color=col, marker=mk, lw=1.2, label=f"numerical, ⟨σv⟩ = {_sci(sv, 1)}")
        ax.semilogx(rows[sel, 0], rows[sel, 6], color=ps.C_MODEL, marker=mk, lw=0.8, ls="--", ms=3)
    ax.plot([], [], color=ps.C_MODEL, ls="--", label="Kolb & Turner x_f")
    ax.set_xlabel("WIMP mass (GeV)")
    ax.set_ylabel("x_f = m / T at freeze-out")
    ax.set_title("(b) Freeze-out at x ≈ 15–30")
    ax.legend(fontsize=6.8, loc="upper left")
    fig.tight_layout()
    kt = a["thermal_ratio_kt"]
    ps.caption(fig, f"Numerical x_f: Y = 2 Y_eq; Kolb & Turner's: Δ(Δ+2) = 1. The logarithm makes x_f insensitive to such choices. "
               f"The textbook formula is good to a few per cent except where g* changes during freeze-out (the QCD era: "
               f"×{kt[0]:.2f} at worst); integrating λ(x) after x_f fixes that.", y=-0.02)
    return ps.save(fig, out / "04_analytic_vs_numerical")


# --- 5. the neutrons --------------------------------------------------------------------------------------


def fig_neutrons(out: Path):
    n = R["neutrons"]
    T = F["n_T"] * 1e3      # MeV
    fig, axs = plt.subplots(1, 3, figsize=(12.0, 4.0))
    ax = axs[0]
    rT = np.interp(np.log(F["n_rate_t"]), np.log(F["n_t"]), np.log(F["n_T"]))
    ax.loglog(np.exp(rT) * 1e3, F["n_rate"], color=ps.INK, lw=1.6, label="weak rate λ(n→p) + λ(p→n)")
    ax.loglog(np.exp(rT) * 1e3, F["n_H"], color=ps.C_MODEL, lw=1.4, ls="--", label="expansion rate H")
    ax.axhline(1 / C.TAU_N, color=ps.MUTED, lw=0.8, ls=":")
    ax.text(5, 1.3 / C.TAU_N, "free decay 1/τ_n", fontsize=7, color=ps.MUTED)
    Tf = n["runs"]["escudero"]["T_freeze_MeV"]
    ax.axvline(Tf, color=ps.AXIS, lw=0.8)
    ax.text(Tf * 0.92, 3e2, f"freeze-out\nT = {Tf:.2f} MeV\n(t = {n['runs']['escudero']['t_freeze_s']:.1f} s)",
            fontsize=7.2, ha="right", color=ps.INK_SECONDARY)
    ax.set_xlim(10, 0.03)
    ax.set_ylim(1e-4, 1e4)
    ax.set_xlabel("photon temperature (MeV)")
    ax.set_ylabel("rate (1/s)")
    ax.set_title("(a) Rate ∝ T⁵ loses to H ∝ T²")
    ax.legend(fontsize=7, loc="lower left")

    ax = axs[1]
    ax.loglog(T, F["n_eq_np"], color=ps.C_MODEL, lw=1.3, ls="--", label="equilibrium e^(−Q/T)")
    ax.loglog(T, F["n_free_np"], color=ps.INK, lw=1.8, label="free n/p: freeze-out, then decay")
    tot_n = F["n_Y_n"] + F["n_Y_d"] + 2 * F["n_Y_t"] + F["n_Y_h"] + 2 * F["n_Y_a"]
    ax.loglog(T, tot_n / (1 - tot_n), color=ps.ORANGE, lw=1.2, ls="-.", label="all neutrons (free + in nuclei)")
    ax.axvline(n["runs"]["escudero"]["T_nuc_keV"] / 1e3, color=ps.AXIS, lw=0.8)
    ax.text(n["runs"]["escudero"]["T_nuc_keV"] / 1e3 * 0.9, 0.5, "helium forms\n"
            f"T ≈ {n['runs']['escudero']['T_nuc_keV']:.0f} keV\nt ≈ {n['runs']['escudero']['t_nuc_s']:.0f} s",
            fontsize=7.2, ha="right", color=ps.INK_SECONDARY)
    T1 = float(np.exp(np.interp(math.log(1.0), np.log(F["n_t"]), np.log(T))))
    ps.annotate(ax, f"n/p = {n['runs']['escudero']['np_1s']:.2f} at t = 1 s\n(still adjusting)",
                (T1, n["runs"]["escudero"]["np_1s"]), (3.5, 0.045), color=ps.INK)
    ax.set_xlim(10, 0.03)
    ax.set_ylim(1e-3, 1.2)
    ax.set_xlabel("photon temperature (MeV)")
    ax.set_ylabel("neutron-to-proton ratio")
    ax.set_title("(b) n/p freezes near 1/5, decays to 1/7")
    ax.legend(fontsize=7, loc="lower left")

    ax = axs[2]
    for key, col, lab, mult in (("a", ps.INK, "⁴He (mass fraction Y = 4 Y_α)", 4), ("d", ps.ORANGE, "deuterium (D/H)", None),
                                ("h", ps.AQUA, "helium-3 (³He/H)", None), ("n", ps.MUTED, "free neutrons", 1)):
        y = F[f"n_Y_{key}"] * mult if mult else F[f"n_Y_{key}"] / F["n_Y_p"]
        ax.loglog(T, np.maximum(y, 1e-30), color=col, lw=1.6 if key == "a" else 1.2, label=lab)
    ax.set_xlim(0.3, 0.008)
    ax.set_ylim(1e-9, 1)
    ax.set_xlabel("photon temperature (MeV)")
    ax.set_ylabel("abundance")
    ax.set_title("(c) Where the neutrons end up")
    ax.legend(fontsize=7, loc="lower right")
    pr = n.get("prymordial", {})
    lines = [f"our Y_p = {n['runs']['escudero']['Y_p']:.4f}  (Born weak rates)"]
    if "born_small_instant_a" in pr:
        lines.append(f"PRyMordial, same physics: {pr['born_small_instant_a']['Y_p']:.4f}")
    if "full" in pr:
        lines.append(f"PRyMordial, full: {pr['full']['Y_p']:.4f}")
    lines.append(f"observed: {C.LIT['yp_aver2021']['value']:.4f} ± {C.LIT['yp_aver2021']['err']:.4f} (Aver+2021)")
    ax.text(0.06, 2.5e-2, "\n".join(lines), fontsize=6.9, color=ps.INK_SECONDARY, va="top")
    fig.tight_layout()
    ps.caption(fig, "Our weak rates (Born approximation, normalised to τ_n = 878.4 s) and an 8-reaction network (PRIMAT rates) "
               "on Escudero's (2020) background. The same freeze-out as the WIMP's, with a linear instead of a quadratic rate.",
               y=-0.02)
    return ps.save(fig, out / "05_neutron_freeze_out")


# --- 6. three freeze-outs -----------------------------------------------------------------------------------


def fig_three(out: Path):
    fig, axs = plt.subplots(1, 3, figsize=(12.2, 3.9))
    # WIMPs: Y vs time (S&S clock)
    ax = axs[0]
    x = F["yx_x"]
    Tw = 100.0 / x
    tw = np.exp(np.interp(np.log(Tw), np.log(F["time_T"][::-1]), np.log(F["time_t"][::-1])))
    ax.loglog(tw, F["yx_eq"], color=ps.C_MODEL, lw=1.3, ls="--", label="equilibrium")
    ax.loglog(tw, F["yx_3"], color=ps.INK, lw=1.8, label="actual (thermal ⟨σv⟩)")
    ex = R["examples"]
    T_geh = 100.0 / ex["rows"][3]["x_gamma_eq_H"]
    t_geh = float(np.exp(np.interp(math.log(T_geh), np.log(F["time_T"][::-1]), np.log(F["time_t"][::-1]))))
    ax.axvline(t_geh, color=ps.AXIS, lw=0.8)
    ax.set_ylim(1e-18, 1e-2)
    ax.set_xlim(tw[0], tw[-1])
    ax.set_xlabel("time since the Big Bang (s)")
    ax.set_ylabel("WIMPs per unit entropy")
    ax.set_title(f"WIMPs (100 GeV): t ≈ {t_geh * 1e9:.0f} ns")
    ax.legend(fontsize=7, loc="lower left")
    # neutrons
    ax = axs[1]
    n = R["neutrons"]
    ax.loglog(F["n_t"], F["n_eq_np"], color=ps.C_MODEL, lw=1.3, ls="--", label="equilibrium e^(−Q/T)")
    ax.loglog(F["n_t"], F["n_free_np"], color=ps.INK, lw=1.8, label="actual (free neutrons)")
    ax.axvline(n["runs"]["escudero"]["t_freeze_s"], color=ps.AXIS, lw=0.8)
    ax.set_ylim(1e-3, 1.2)
    ax.set_xlim(0.01, 1e4)
    ax.set_xlabel("time since the Big Bang (s)")
    ax.set_ylabel("neutrons per proton")
    ax.set_title(f"Neutrons: t ≈ {n['runs']['escudero']['t_freeze_s']:.1f} s")
    ax.legend(fontsize=7, loc="lower left")
    # electrons
    ax = axs[2]
    yr = 3.15576e7
    ax.loglog(F["rec_t"] / yr, F["rec_saha"], color=ps.C_MODEL, lw=1.3, ls="--", label="equilibrium (Saha)")
    ax.loglog(F["rec_t"] / yr, F["rec_xe"], color=ps.INK, lw=1.8, label="actual (Peebles atom)")
    zc, xc = F["rec_zc"], F["rec_xc"]
    tc = np.exp(np.interp(-zc, -F["rec_z"], np.log(F["rec_t"])))
    sel = (zc <= F["rec_z"][0]) & (zc >= F["rec_z"][-1])
    ax.loglog(tc[sel] / yr, xc[sel], color=ps.ORANGE, lw=1.0, ls=":", label="CAMB (RECFAST)")
    rec = R["recombination"]
    tz = float(np.exp(np.interp(-rec["z_freeze"], -F["rec_z"], np.log(F["rec_t"]))))
    ax.axvline(tz / yr, color=ps.AXIS, lw=0.8)
    ax.set_ylim(1e-6, 1.5)
    ax.set_xlim(2e5, 1e8)
    ax.set_xlabel("time since the Big Bang (years)")
    ax.set_ylabel("free electrons per hydrogen")
    ax.set_title(f"Electrons: z ≈ {rec['z_freeze']:.0f}, t ≈ {tz / yr / 1e6:.1f} Myr")
    ax.legend(fontsize=7, loc="lower left")
    fig.suptitle("Three freeze-outs: a reaction falls behind the expansion, and what is left stays", x=0.01,
                 ha="left", fontsize=11)
    fig.tight_layout()
    ps.caption(fig, "Vertical lines: where the reaction rate per particle equals H. In each case equilibrium (dashed) would "
               "leave almost nothing; the frozen remnant is the dark matter, the neutrons that became helium, and the "
               f"x_e ≈ {_sci(rec['x_e_200'], 1)} electrons per hydrogen that kept the gas coupled to the CMB.", y=-0.02)
    return ps.save(fig, out / "06_three_freeze_outs")


# --- 7. gamma rays from dwarfs ---------------------------------------------------------------------------


def fig_dwarfs(out: Path):
    L = R["dwarfs"]["limits"]
    fig, axs = plt.subplots(1, 2, figsize=(11.2, 4.6), sharey=True)
    for ax, ch, title in ((axs[0], "bb", "annihilation into b quarks"), (axs[1], "tau", "annihilation into τ leptons")):
        th_m, th = F["dw_th_m"], F["dw_th"]
        cross = L[f"{ch}_pppc_benchmark"]["crossing_GeV"]
        if cross:
            ax.axvspan(th_m[0], cross[0], color=ps.lighten(ps.C_MODEL, 0.88), lw=0)
            ax.text(th_m[0] * 1.08, 1.12e-27, f"thermal WIMPs excluded\nbelow {cross[0]:.0f} GeV",
                    fontsize=7.6, ha="left", va="bottom", color=ps.C_MODEL)
        ax.loglog(th_m, th, color=ps.C_MODEL, lw=2.4, label="thermal relic ⟨σv⟩ (this project)")
        ax.loglog(F[f"dw_{ch}_pppc_benchmark_m"], F[f"dw_{ch}_pppc_benchmark_ul"], color=ps.INK, lw=1.8,
                  label="our limit: 42 dwarfs, 14 yr (McDaniel et al. 2024 data)")
        ax.loglog(F[f"dw_{ch}_pppc_measured_m"], F[f"dw_{ch}_pppc_measured_ul"], color=ps.INK, lw=0.9, ls="--",
                  label="our limit: 30 dwarfs with measured J")
        m = F[f"dw_{ch}_mcdsum_benchmark_m"]
        ok = m >= C.GAMMA_MASS_MIN
        ax.loglog(m[ok], F[f"dw_{ch}_mcdsum_benchmark_ul"][ok], color=ps.ORANGE, lw=1.1, ls=(0, (4, 2)),
                  label="their published TS profiles, summed")
        ax.loglog(F[f"dw_{ch}_cosmixs_benchmark_m"], F[f"dw_{ch}_cosmixs_benchmark_ul"], color=ps.MUTED, lw=0.9, ls=":",
                  label="ours with CosmiXs spectra")
        ax.loglog(F[f"dw_{ch}_ack_m"], F[f"dw_{ch}_ack_ul"], color=ps.AQUA, lw=1.3, marker="o", ms=2.5,
                  label="Ackermann et al. 2015 (6 yr, 15 dwarfs)")
        g = R["dwarfs"]["gce_calore2015"][ch]
        ax.errorbar([g["mass"]], [g["sigv"]], xerr=[[g["mass_err"][0]], [g["mass_err"][1]]],
                    yerr=[[g["sigv_err"][0]], [g["sigv_err"][1]]], fmt="D", color=ps.INK_SECONDARY, ms=4, lw=0.9,
                    capsize=0, label="Galactic Centre excess fit (Calore et al. 2015)")
        if ch == "bb":
            d = R["dwarfs"]["gce_daylan2016"]
            ax.add_patch(plt.Rectangle((d["bb_mass_range"][0], d["sigv_range"][0]), d["bb_mass_range"][1] - d["bb_mass_range"][0],
                                       d["sigv_range"][1] - d["sigv_range"][0], fill=False, ec=ps.INK_SECONDARY, lw=0.7, ls="--"))
            ax.text(36, 3.3e-26, "Daylan+16", fontsize=6.5, color=ps.INK_SECONDARY)
        ax.set_xlim(C.GAMMA_MASS_MIN, 1e4)
        ax.set_ylim(1e-27, 1e-23)
        ax.set_xlabel("WIMP mass (GeV)")
        ax.set_title(title)
    axs[0].set_ylabel("⟨σv⟩ today, 95% upper limit (cm³/s)")
    axs[1].legend(fontsize=6.9, loc="upper left")
    fig.suptitle("Dwarf galaxies see no annihilation light: light thermal WIMPs that annihilate to quarks or leptons are ruled out",
                 x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout()
    ps.caption(fig, "Our limit combines McDaniel et al.'s per-dwarf, per-energy-bin likelihoods (fermipy SEDs) with PPPC 4 DM ID "
               "spectra integrated over each bin and log-normal J-factor priors profiled out. Assumes 100% annihilation into "
               "one channel and a Majorana WIMP; the Galactic Centre point depends on the Milky Way's halo profile.", y=-0.01)
    return ps.save(fig, out / "07_dwarf_galaxy_limits")


# --- 8. direct detection --------------------------------------------------------------------------------


def fig_direct(out: Path):
    D = R["direct"]
    fig, ax = plt.subplots(figsize=(7.8, 5.4))
    fm, fs = F["dd_floor_m"], F["dd_floor"]
    ax.fill_between(fm, 1e-52, fs, color=ps.lighten(ps.AQUA, 0.85), lw=0)
    ax.loglog(fm, fs, color=ps.AQUA, lw=1.0, ls="--")
    ax.text(2.5e3, 1.6e-49, "neutrino fog (xenon)", fontsize=7.4, color=ps.AQUA, rotation=14)
    for i in range(2):
        reg = F[f"dd_dama{i}"]
        ax.fill(reg[:, 0], reg[:, 1], color=ps.lighten(ps.MUTED, 0.6), lw=0)
    ax.text(19, 1.6e-40, "DAMA", fontsize=7, color=ps.MUTED)
    for k, lab_m in (("LUX 2017", 3e3), ("PandaX-II 2017", None), ("XENON1T 2018", 5e3), ("CRESST", None),
                     ("CDMSlite", None), ("DarkSide-50", None)):
        ax.loglog(F[f"dd_hist_{k}_m"], F[f"dd_hist_{k}"], color=ps.MUTED, lw=0.8)
    ax.text(900, 3.2e-44, "LUX 2017, PandaX-II 2017,\nXENON1T 2018", fontsize=6.8, color=ps.MUTED)
    ax.loglog(F["dd_xn_m"], F["dd_xn"], color=ps.ORANGE, lw=1.4, label="XENONnT 2025 (3.1 t·yr)")
    ax.loglog(F["dd_lz_m"], F["dd_lz"], color=ps.INK, lw=1.8, label="LZ 2024 (4.2 t·yr)")
    mx = np.array([F["dd_lz_m"][-1], 3e4])
    ax.loglog(mx, F["dd_lz"][-1] * mx / mx[0], color=ps.INK, lw=1.0, ls="--")
    ax.text(1.15e4, F["dd_lz"][-1] * 0.55, "∝ m", fontsize=7, color=ps.INK_SECONDARY)
    hm, hs, hl = F["hp_m"], F["hp_sigma"], F["hp_lambda"]
    ax.loglog(hm, hs, color=ps.C_MODEL, lw=2.2, label="Higgs-portal singlet that freezes out to Ω h² = 0.120")
    for w in D["higgs_portal"]["allowed_by_lz2024"]:
        sel = (hm >= w[0]) & (hm <= w[1])
        ax.loglog(hm[sel], hs[sel], color=ps.C_MODEL, lw=5.0, alpha=0.45, solid_capstyle="butt")
    ax.set_xlim(1, 3e4)
    ax.set_ylim(1e-50, 1e-39)
    ax.set_xlabel("WIMP mass (GeV)")
    ax.set_ylabel("spin-independent WIMP–nucleon cross section (cm²)")
    ax.set_title("Underground: no WIMP has hit a xenon nucleus", pad=8)
    ax.legend(fontsize=7.2, loc="upper right")
    win = D["higgs_portal"]["allowed_by_lz2024"]
    txt = ("the singlet survives LZ only at " + ", ".join(f"{a:.1f}–{b:.1f} GeV" for a, b in win)
           + "\n(the Higgs resonance, m ≈ m_h/2)") if win else "the singlet is excluded everywhere"
    ps.annotate(ax, txt, (61, 1.5e-49), (130, 2.2e-50), color=ps.C_MODEL)
    ps.caption(fig, "Limits: 90% CL, standard halo (ρ = 0.3 GeV/cm³). Blue: our freeze-out coupling λ_hs(m_S) "
               "(Cline et al. 2013 cross sections, LHC Higgs WG widths, exact thermal average) turned into a scattering rate "
               "with f_N = 0.30 ± 0.03 (±20% in σ). Thick: allowed by LZ 2024 and by the invisible Higgs width.", y=-0.01)
    return ps.save(fig, out / "08_direct_detection")


# --- 0. where the thermal WIMP stands (keeper) ----------------------------------------------------------


def fig_status(out: Path):
    L = R["dwarfs"]["limits"]
    D = R["direct"]
    fig = plt.figure(figsize=(9.4, 6.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.6], hspace=0.32)
    ax = fig.add_subplot(gs[0])
    m = F["th_m"]
    ax.fill_between(m, F["th_eos_ss2020_lo"], F["th_eos_ss2020_hi"], color=ps.lighten(ps.C_MODEL, 0.75), lw=0)
    ax.loglog(m, F["th_main"], color=ps.C_MODEL, lw=2.6, label="thermal relic: ⟨σv⟩ for Ω_c h² = 0.120")
    ax.loglog(F["dw_bb_pppc_benchmark_m"], F["dw_bb_pppc_benchmark_ul"], color=ps.INK, lw=1.6,
              label="Fermi-LAT dwarfs (our limit): → b b̄")
    ax.loglog(F["dw_tau_pppc_benchmark_m"], F["dw_tau_pppc_benchmark_ul"], color=ps.ORANGE, lw=1.6,
              label="… → τ⁺τ⁻")
    ax.axhline(3e-26, color=ps.MUTED, lw=0.8, ls=":")
    ax.text(1.2e4, 3.25e-26, "textbook 3×10⁻²⁶", fontsize=7, color=ps.MUTED, va="bottom")
    ax.set_xlim(0.01, 3e5)
    ax.set_ylim(1e-27, 1e-24)
    ax.set_ylabel("⟨σv⟩ (cm³/s)")
    ax.legend(fontsize=7.4, loc="lower left")
    ax.set_title("Where the thermal WIMP stands (October 2026)")
    ax.axvspan(1.4e5, 3e5, color=ps.lighten(ps.MUTED, 0.8), lw=0)
    ax.text(1.5e5, 1.5e-27, "unitarity\nlimit", fontsize=7, color=ps.INK_SECONDARY)

    ax2 = fig.add_subplot(gs[1], sharex=ax)
    rows = []
    cb = L["bb_pppc_benchmark"]["crossing_GeV"]
    ct = L["tau_pppc_benchmark"]["crossing_GeV"]
    rows.append(("γ rays, dwarfs: if WIMPs → b b̄", [(C.GAMMA_MASS_MIN, cb[0])] if cb else [], ps.INK))
    rows.append(("γ rays, dwarfs: if WIMPs → τ⁺τ⁻", [(C.GAMMA_MASS_MIN, ct[0])] if ct else [], ps.ORANGE))
    win = D["higgs_portal"]["allowed_by_lz2024"]
    hm = F["hp_m"]
    excl = []
    lo = float(hm[0])
    for a, b in win:
        if a > lo:
            excl.append((lo, a))
        lo = b
    if lo < hm[-1]:
        excl.append((lo, float(hm[-1])))
    rows.append(("LZ 2024: the Higgs-portal singlet", excl, ps.C_MODEL))
    rows.append(("unitarity (Smirnov & Beacom 2019)", [(1.4e5, 3e5)], ps.MUTED))
    survive = (f"; the Higgs resonance {win[0][0]:.1f}–{win[0][1]:.1f} GeV survives" if win else "")
    notes = {0: "above ~140 TeV (s-wave, self-conjugate)",
             1: f"excluded from {hm[0]:.0f} GeV to {hm[-1] / 1e3:.0f} TeV" + survive,
             2: f"excluded below {ct[0]:.0f} GeV" if ct else "", 3: f"excluded below {cb[0]:.0f} GeV" if cb else ""}
    for i, (lab, spans, col) in enumerate(rows[::-1]):
        for a, b in spans:
            ax2.barh(i, b - a, left=a, height=0.5, color=ps.lighten(col, 0.25), lw=0)
        ax2.text(0.012, i, lab, va="center", fontsize=7.6)
        if spans and notes.get(i):
            if i == 0:
                ax2.text(spans[0][0] * 0.9, i + 0.32, notes[i], va="bottom", ha="right", fontsize=6.9, color=ps.INK_SECONDARY)
            else:
                ax2.text(spans[0][0], i + 0.32, notes[i], va="bottom", ha="left", fontsize=6.9, color=ps.INK_SECONDARY)
    ax2.set_yticks([])
    ax2.set_ylim(-0.6, len(rows) - 0.4)
    ax2.set_xscale("log")
    ax2.set_xlabel("WIMP mass (GeV)")
    ax2.grid(axis="y", visible=False)
    ax2.set_title("Masses ruled out (bars), for the stated assumption", fontsize=9.5)
    ps.caption(fig, "Each bar assumes something: 100% annihilation into one channel (dwarfs), or the minimal Higgs-portal model "
                    "(LZ). The gaps are not a prediction: other final states, p-wave annihilation, co-annihilation or a dark "
                    "sector (and non-WIMP candidates such as axions) escape these searches. Unitarity: s-wave, self-conjugate, "
                    "with bound states.", y=0.0)
    return ps.save(fig, out / "00_wimp_status")


def make_all(results=None) -> list:
    _load()
    ps.use_style()
    out = C.FIG_DIR
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for f in (fig_status, fig_freeze_out, fig_thermal, fig_conventions, fig_analytic, fig_neutrons, fig_three,
              fig_dwarfs, fig_direct):
        written += f(out)
    return written


if __name__ == "__main__":
    for p in make_all():
        print("  wrote", Path(p).relative_to(C.HERE))
