"""The project's figures (PNG + SVG in figures/).

    python figures.py            # all figures from the cached analysis
    python figures.py 1 4        # only figures 1 and 4
"""

from __future__ import annotations

import json
import math
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import baryons
import config as cfg
import literature as lit
import plotstyle as ps
from analysis import load_pickle


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(v: float, digits: int = 1) -> str:
    """3.1e-19 -> '3.1 × 10⁻¹⁹'."""
    if v == 0:
        return "0"
    e = int(math.floor(math.log10(abs(v))))
    m = v / 10**e
    if round(m, digits) >= 10:
        m, e = m / 10, e + 1
    return f"{m:.{digits}f} × 10{str(e).translate(_SUP)}"


def _res() -> dict:
    return json.loads(cfg.RESULTS_PATH.read_text())


def _secondary_ombh2(ax, axis: str = "x"):
    """A top axis in Omega_b h^2 for an eta axis (the same quantity, converted)."""
    k = baryons.eta10_coefficient()
    if axis == "x":
        sec = ax.secondary_xaxis("top", functions=(lambda e: e / k, lambda o: o * k))
        sec.set_xlabel("Ω_b h²  (same quantity, converted)", color=ps.INK_SECONDARY)
    sec.tick_params(colors=ps.MUTED, labelsize=7.5)
    return sec


# ---------------------------------------------------------------------------
# 1. eta two ways: deuterium (3 minutes) and the CMB (380,000 years)
# ---------------------------------------------------------------------------


def fig_eta_two_ways():
    res = _res()["bbn"]
    curves = load_pickle("bbn_curves")
    fig = plt.figure(figsize=(14.5, 4.9))
    axes = [fig.add_axes([0.045, 0.14, 0.24, 0.66]), fig.add_axes([0.345, 0.14, 0.25, 0.66]),
            fig.add_axes([0.795, 0.14, 0.195, 0.66])]
    title_pad = 30

    # (a) the deuterium plateau
    ax = axes[0]
    dh = res["dh_measured"]["pdg2026_12"]
    for r in lit.DH_SYSTEMS:
        name, z, met, (v, lo, hi), ref, cooke = r
        ax.errorbar(met, v, yerr=[[lo], [hi]], fmt="o" if cooke else "s", ms=5,
                    color=ps.C_BBN if cooke else ps.lighten(ps.C_BBN, 0.0),
                    mfc=ps.C_BBN if cooke else ps.SURFACE, mec=ps.C_BBN, ecolor=ps.C_BBN,
                    elinewidth=0.9, capsize=0, zorder=3)
    ax.axhspan(dh["mean"] - dh["err_scaled"], dh["mean"] + dh["err_scaled"],
               color=ps.RAMPS[ps.C_BBN][0], zorder=1, lw=0)
    ax.axhline(dh["mean"], color=ps.C_BBN, lw=1.0, zorder=2)
    ax.set_xlabel("metallicity of the gas cloud  [O/H]  (log, relative to the Sun)")
    ax.set_ylabel("deuterium / hydrogen  (× 10⁻⁵)")
    ax.set_ylim(1.5, 3.4)
    ax.set_title("a  Primordial deuterium in 12 gas clouds", pad=title_pad)
    ax.text(0.03, 0.96, f"weighted mean {dh['mean']:.3f} ± {dh['err_scaled']:.3f}\n"
            f"(PDG 2025/26 selection: {dh['mean'] * 10:.2f} ± {dh['err_scaled'] * 10:.2f} ×10⁻⁶)",
            transform=ax.transAxes, va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.text(0.03, 0.06, "● Cooke+18 precision sample   □ later systems\n"
            "no trend with metallicity: this is the Big Bang's deuterium",
            transform=ax.transAxes, fontsize=7.0, color=ps.MUTED)

    # (b) D/H against eta from PRyMordial, both rate compilations
    ax = axes[1]
    eta_plot = np.geomspace(5.4e-10, 6.8e-10, 200)
    styles = {"primat": ("PRIMAT rates (theory d+d)", "-"),
              "nacreii": ("NACRE-II rates (data d+d)", "--")}
    for key, (lab, ls) in styles.items():
        c = curves[key]
        y = np.exp(np.interp(np.log(eta_plot), np.log(c["eta"]), np.log(c["dh"])))
        ax.fill_between(eta_plot * 1e10, y * (1 - c["frac_theory"]), y * (1 + c["frac_theory"]),
                        color=ps.RAMPS[ps.C_MODEL][0], alpha=0.55 if key == "primat" else 0.3, lw=0)
        ax.plot(eta_plot * 1e10, y, ls, color=ps.C_MODEL, lw=1.6,
                label=f"{lab}: ±{100 * c['frac_theory']:.1f}% theory")
    ax.axhspan(dh["mean"] - dh["err_scaled"], dh["mean"] + dh["err_scaled"],
               color=ps.RAMPS[ps.C_BBN][0], lw=0, zorder=0)
    ax.axhline(dh["mean"], color=ps.C_BBN, lw=1.0)
    cmb = [c for c in res["cmb"] if c["label"].startswith("Planck 2018 TT,TE")][0]
    ax.axvspan((cmb["eta"] - cmb["eta_err"]) * 1e10, (cmb["eta"] + cmb["eta_err"]) * 1e10,
               color=ps.RAMPS[ps.C_CMB][0], lw=0, zorder=0)
    ax.text(cmb["eta"] * 1e10 + 0.045, 2.13, "CMB (Planck)", color=ps.C_CMB, fontsize=7.5,
            va="bottom")
    ax.text(5.45, dh["mean"] + 0.035, "measured D/H", color=ps.C_BBN, fontsize=7.5)
    ax.set_xlim(5.4, 6.8)
    ax.set_ylim(2.1, 3.0)
    ax.set_xlabel("baryons per photon  η  (× 10⁻¹⁰)")
    ax.set_ylabel("predicted deuterium / hydrogen  (× 10⁻⁵)")
    ax.set_title("b  More baryons, less deuterium left over", pad=title_pad)
    ax.legend(loc="upper right", fontsize=7)
    sec = ax.secondary_xaxis("top", functions=(lambda e: e / baryons.eta10_coefficient(),
                                                lambda o: o * baryons.eta10_coefficient()))
    sec.set_xlabel("Ω_b h²", color=ps.INK_SECONDARY, fontsize=8)
    sec.tick_params(colors=ps.MUTED, labelsize=7)

    # (c) the two clocks on one axis
    ax = axes[2]
    rows = []
    for key, lab in (("primat", "BBN, our run: PRIMAT rates"), ("nacreii", "BBN, our run: NACRE-II rates")):
        r = res["rates"][key]["pdg2026_12"]
        rows.append((lab, r["eta_median"] * 1e10, r["eta_err"] * 1e10, "bbn", True))
    for r in res["literature"]:
        rows.append(("BBN, " + r["label"], baryons.eta_from_ombh2(r["ombh2"]) * 1e10,
                     baryons.eta_from_ombh2(r["err"]) * 1e10, "bbn", False))
    for c in res["cmb"]:
        mine = c["label"].startswith("D3")
        rows.append(("CMB, " + c["label"], c["eta"] * 1e10, c["eta_err"] * 1e10, "cmb", mine))
    ylab = []
    for i, (lab, v, e, kind, ours) in enumerate(rows):
        col = ps.C_BBN if kind == "bbn" else ps.C_CMB
        ax.errorbar(v, -i, xerr=e, fmt="*" if ours else "o", ms=8 if ours else 4.5, color=col,
                    mec=ps.INK if ours else col, mew=0.6 if ours else 0, elinewidth=1.1, capsize=0)
        ylab.append(lab)
    ax.set_yticks(-np.arange(len(rows)))
    ax.set_yticklabels(ylab, fontsize=7)
    for t, (_, _, _, kind, _) in zip(ax.get_yticklabels(), rows):
        t.set_color(ps.C_BBN if kind == "bbn" else ps.C_CMB)
    ax.set_xlabel("baryons per photon  η  (× 10⁻¹⁰)")
    ax.set_title("c  Two clocks, one number   (★ = computed in this atlas)", pad=title_pad,
                 x=-0.75, ha="left")
    ax.grid(axis="y", visible=False)
    ax.set_xlim(5.75, 6.4)
    sec = ax.secondary_xaxis("top", functions=(lambda e: e / baryons.eta10_coefficient(),
                                                lambda o: o * baryons.eta10_coefficient()))
    sec.set_xlabel("Ω_b h²", color=ps.INK_SECONDARY, fontsize=8)
    sec.tick_params(colors=ps.MUTED, labelsize=7)
    ps.caption(fig, "Deuterium made 3 minutes after the Big Bang (orange) and the CMB's sound waves "
                    "380,000 years later (violet) measure the same ratio: about 6.1 baryons per 10 billion "
                    "photons. The spread between the two BBN rate compilations is the main systematic.",
               y=0.035)
    return ps.save(fig, cfg.FIGURE_DIR / "01_eta_two_ways")


# ---------------------------------------------------------------------------
# 2. The asymmetry in numbers
# ---------------------------------------------------------------------------


def fig_asymmetry():
    r = _res()["baryons"]
    cur = load_pickle("baryons_curves")
    pl = r["planck"]
    fig = plt.figure(figsize=(14.5, 4.6))
    axes = [fig.add_axes([0.06, 0.15, 0.25, 0.66]), fig.add_axes([0.385, 0.15, 0.26, 0.66]),
            fig.add_axes([0.72, 0.15, 0.27, 0.66])]

    # (a) today, per cubic metre
    ax = axes[0]
    n_g = pl["n_gamma_per_cm3"] * 1e6
    n_nu = 3 * (3 / 11) * n_g                           # three species, particles + antiparticles
    items = [("CMB photons", n_g, ps.C_CMB), ("relic neutrinos", n_nu, ps.C_OTHER),
             ("protons + neutrons\n(baryons)", pl["n_b_per_m3"], ps.C_MODEL),
             ("antibaryons", None, None)]
    for i, (lab, v, col) in enumerate(items):
        if v is None:
            ax.text(1.4e-2, -i, "none found (but half of the relic\nneutrinos are antineutrinos)",
                    va="center", fontsize=7.5, color=ps.INK_SECONDARY)
            continue
        ax.barh(-i, v, left=1e-2, height=0.55, color=col, lw=0)
        txt = f"{v / 1e6:.0f} million" if v > 1e6 else f"{v:.2f}"
        ax.text(v * 1.6, -i, txt, va="center", fontsize=8.5, color=ps.INK)
    ax.set_xscale("log")
    ax.set_xlim(1e-2, 1e12)
    ax.set_yticks(-np.arange(len(items)))
    ax.set_ylim(-3.6, 0.6)
    ax.set_yticklabels([it[0] for it in items], fontsize=8)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("number per cubic metre of today's universe")
    ax.set_title("a  One baryon per 1.6 billion photons", pad=10)
    ax.text(0.98, 0.36, f"η = n_b / n_γ = {pl['eta'] * 1e10:.3f} × 10⁻¹⁰\n(Planck 2018)\n"
            f"{pl['photons_per_baryon'] / 1e9:.2f} billion photons\nper baryon",
            transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY, ha="right", va="center")

    # (b) the excess per quark before the quarks were bound into protons
    ax = axes[1]
    T = cur["T_quark"]
    ax.plot(T, cur["quark_excess"] * 1e9, color=ps.C_MODEL, lw=1.8)
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_ylim(0, 2.4)
    ax.set_xlabel("temperature of the universe  (GeV; time runs to the right →)")
    ax.set_ylabel("extra quarks per quark + antiquark  (× 10⁻⁹)")
    ax.set_title("b  Before 20 µs: quarks and antiquarks almost equal", pad=10)
    pairs = r["quark_excess"]["pairs_per_extra_quark_1GeV"]
    ax.text(3e3, 2.15, f"≈ 1 extra quark for every {pairs / 1e6:.0f} million\nquark–antiquark pairs",
            fontsize=8, color=ps.INK, va="top")
    for t, lab in ((159.5, "electroweak\ncrossover"), (0.156, "quarks bind\ninto hadrons")):
        ax.axvline(t, color=ps.AXIS, lw=0.8, ls=":")
        ax.text(t * 0.9, 0.12, lab, fontsize=7, color=ps.MUTED, ha="left")

    # (c) the annihilation catastrophe
    ax = axes[2]
    sym, real = cur["symmetric"], cur["real"]
    Tm = lambda d: d["T"] * 1e3
    ax.plot(Tm(real), real["Y_eq"], ":", color=ps.C_OTHER, lw=1.2, label="in equilibrium (no expansion)")
    ax.plot(Tm(sym), sym["Y_N"], color=ps.C_ALT, lw=2.0, label="symmetric universe (no excess)")
    ax.plot(Tm(real), real["Y_N"], color=ps.C_MODEL, lw=2.0, label="our universe: nucleons")
    ax.plot(Tm(real), np.maximum(real["Y_Nbar"], 1e-40), "--", color=ps.C_MODEL, lw=1.6,
            label="our universe: antinucleons")
    c = r["catastrophe"]
    ax.axhline(real["Y_B"], color=ps.INK, lw=0.8, ls="-", alpha=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(300, 2)
    ax.set_ylim(1e-22, 1e-1)
    ax.set_xlabel("temperature  (MeV; time runs to the right →)")
    ax.set_ylabel("nucleons per unit entropy  n / s")
    ax.set_title("c  Without the excess: the annihilation catastrophe", pad=10)
    ax.annotate(f"freeze-out at {c['T_freezeout_MeV']:.0f} MeV:\n"
                f"η = {sci(c['eta_symmetric'])},\n{c['observed_over_symmetric'] / 1e9:.1f} billion times\n"
                "fewer baryons than we see", xy=(8, sym["Y_N"][-1]), xytext=(9, 2e-17),
                fontsize=7.5, color=ps.INK_SECONDARY,
                arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "lw": 0.7})
    ax.text(2.3, real["Y_B"] * 2.5, "the excess, Y_B = 8.7 × 10⁻¹¹", fontsize=7.5, color=ps.INK,
            ha="right")
    ax.legend(loc="lower left", fontsize=7)
    ps.caption(fig, "Today's universe holds 1.6 billion photons for every proton or neutron (a). "
                    "Before 20 microseconds there was one extra quark for every ~270 million quark–antiquark pairs (b; "
                    "the wiggles are heavy quarks dropping out as the universe cools below their mass). "
                    "As the universe cooled, every antinucleon met a nucleon; only the excess survived (c). "
                    "With no excess, freeze-out at 21 MeV would leave a universe almost empty of matter.",
               y=0.03)
    return ps.save(fig, cfg.FIGURE_DIR / "02_asymmetry_in_numbers")


# ---------------------------------------------------------------------------
# 3. Where matter meets antimatter: the flow into a domain boundary
# ---------------------------------------------------------------------------


def fig_boundary():
    import domains as dm

    d = load_pickle("domain_curves")
    res = _res()["domains"]
    runs, ref = d["runs"], d["cdg_ref"]
    fig, axs = plt.subplots(1, 4, figsize=(15.5, 4.3))
    fig.subplots_adjust(wspace=0.36, left=0.05, right=0.99, bottom=0.17, top=0.86)
    styles = {"cdg": ("-", "CDG's heating shortcut (conservative)"),
              "full": ("--", "full electron heating"),
              "none": (":", "no heating")}

    # (a) the annihilation rate
    ax = axs[0]
    ax.plot(ref["J"]["y"], ref["J"]["J_cm2_s"], color=ps.INK, lw=2.4,
            alpha=0.35, label="CDG 1998, their universe (read from their figure)")
    ax.plot(d["cdg_run"]["y"], d["cdg_run"]["J"], color=ps.INK, lw=1.0,
            label="our solution, their universe")
    for k, (ls, lab) in styles.items():
        ax.plot(runs[k]["y"], runs[k]["J"], ls, color=ps.C_MODEL, lw=1.6, label=f"ours, Planck universe: {lab}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1100, 20)
    ax.set_xlabel("1 + z  (time runs to the right →)")
    ax.set_ylabel("annihilations per cm² of boundary per second")
    ax.set_title("a  Annihilation rate at a boundary")
    ax.legend(fontsize=6.3, loc="lower left")

    # (b) the temperature of the gas
    ax = axs[1]
    t = ref["T"]
    ax.plot(t["electron_heating"]["y"], t["electron_heating"]["T_eV"], color=ps.INK, lw=2.4, alpha=0.35,
            label="CDG 1998 (their figure)")
    ax.plot(d["cdg_run"]["y"], d["cdg_run"]["T"], color=ps.INK, lw=1.0, label="our solution, their universe")
    for k, (ls, lab) in styles.items():
        ax.plot(runs[k]["y"], runs[k]["T"], ls, color=ps.C_MODEL, lw=1.6)
    y = np.geomspace(20, 1100, 50)
    u = dm.planck2018()
    ax.plot(y, u.t_gamma_ev(y), color=ps.C_OTHER, lw=1.0, ls="-.", label="CMB temperature")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1100, 20)
    ax.set_xlabel("1 + z")
    ax.set_ylabel("gas temperature near the boundary (eV)")
    ax.set_title("b  Annihilation electrons heat the gas")
    ax.legend(fontsize=6.5, loc="upper left")

    # (c) the depletion zone
    ax = axs[2]
    prof = runs["cdg"]["profiles"]
    keys = sorted(prof, reverse=True)
    cols = ps.SEQ_BLUE[2:2 + len(keys)]
    for k, col in zip(keys, cols):
        l, n, T, v = prof[k]
        ax.plot(l / cfg.MPC * 1e3, n, color=col, lw=1.6, label=f"1 + z = {k}")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(1e-5, 1e2)
    ax.set_xlabel("distance from the boundary (kpc, proper)")
    ax.set_ylabel("gas density / density far away")
    ax.set_title("c  The gas thins out near the boundary")
    ax.legend(fontsize=6.8, loc="lower right")

    # (d) the length scales
    ax = axs[3]
    ys = np.geomspace(20, 1100, 60)
    rng = [dm.ElectronHeating(u, yy).range_m / cfg.MPC * 1e3 for yy in ys]
    ax.plot(ys, rng, color=ps.C_ALT, lw=1.6, label="range of annihilation electrons")
    ax.plot(runs["cdg"]["y"], runs["cdg"]["D"] / cfg.MPC * 1e3, color=ps.C_MODEL, lw=1.6,
            label="depletion zone (density halved)")
    ax.plot(ys, cfg.C_LIGHT / u.hubble(ys) / cfg.MPC * 1e3, color=ps.C_OTHER, lw=1.2,
            label="Hubble radius c/H")
    d0 = res["modern"]["headline_d0_min_gpc"]
    ax.plot(ys, d0 * 1e6 / ys, color=ps.INK, lw=1.0, ls="--",
            label=f"domain scale for d₀ = {d0:.0f} Gpc")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1100, 20)
    ax.set_xlabel("1 + z")
    ax.set_ylabel("proper length (kpc)")
    ax.set_title("d  Scales: why a 1D flow is enough")
    ax.legend(fontsize=6.5, loc="lower right")
    ps.caption(fig, "Our numerical solution of Cohen, De Rújula & Glashow's fluid equations reproduces their "
                    "annihilation rate within 2–31% (z ≥ 29) and their temperature within 1–17% (grey: theirs, ink: ours). "
                    "For today's universe (blue: 3× more baryons, slower expansion) the boundary annihilates 5–10× faster.",
               y=0.02)
    return ps.save(fig, cfg.FIGURE_DIR / "03_domain_boundary")


# ---------------------------------------------------------------------------
# 4. The MeV sky: no room for antimatter domains
# ---------------------------------------------------------------------------


def fig_mev_sky():
    import fetch_data as fd

    d = load_pickle("domain_curves")
    res = _res()["domains"]
    m = res["modern"]
    Eg, F = d["Eg"], d["spectra_1gpc"]["cdg"]
    fig = plt.figure(figsize=(15.5, 5.0))
    ax = fig.add_axes([0.05, 0.14, 0.42, 0.72])
    # data (direct labels)
    bat = lit.swift_bat()
    ax.errorbar(bat["E_mid"], bat["E2I"], yerr=bat["E2I_err"], fmt="o", ms=2.5, color=ps.C_OTHER,
                elinewidth=0.6)
    ax.text(0.03, 0.12, "Swift-BAT\n(X-ray background)", fontsize=7, color=ps.C_OTHER, ha="center")
    for w, mk in (("weidenspointner", "s"), ("kappadath", "D")):
        c = lit.comptel(w)
        e = c["E_mid"]
        filled = w == "weidenspointner"
        ax.errorbar(e, e**2 * c["I"], yerr=e**2 * c["err"], xerr=[e - c["E_lo"], c["E_hi"] - e],
                    fmt=mk, ms=4.5 if filled else 4, color=ps.INK,
                    mfc=ps.INK if filled else ps.SURFACE, elinewidth=0.8, zorder=5)
    ax.text(25, 2.5e-4, "COMPTEL\n■ Weidenspointner 1999\n◇ Kappadath 1998", fontsize=7, color=ps.INK,
            ha="center", va="center")
    eg = lit.egret()
    ok = ~eg["is_ul"]
    ax.errorbar(eg["E_mid"][ok], eg["E_mid"][ok] ** 2 * eg["I"][ok],
                yerr=eg["E_mid"][ok] ** 2 * eg["err"][ok], fmt="^", ms=3.5, color=ps.C_OTHER,
                elinewidth=0.7)
    ax.text(70, 3.4e-3, "EGRET", fontsize=7, color=ps.C_OTHER, ha="center")
    fl = {k: fd.load_fermi_igrb(k) for k in "ABC"}
    e = fl["A"]["E_mid"]
    lo = np.min([fl[k]["igrb"] - fl[k]["igrb_err_lo"] for k in "ABC"], axis=0)
    hi = np.max([fl[k]["igrb"] + fl[k]["igrb_err_up"] for k in "ABC"], axis=0)
    ax.fill_between(e, e**2 * np.maximum(lo, 1e-30), e**2 * hi, color=ps.RAMPS[ps.INK][0], lw=0,
                    step="mid")
    ax.text(600, 1.8e-4, "Fermi-LAT\nisotropic", fontsize=7, color=ps.MUTED, ha="center")
    # predictions, labelled at their peaks
    d0s = [(1.0, ps.SEQ_BLUE[4]), (10.0, ps.SEQ_BLUE[5]), (m["headline_d0_min_gpc"], ps.SEQ_BLUE[7])]
    k = int(np.argmax(Eg**2 * F))
    for i, (d0, col) in enumerate(d0s):
        ax.plot(Eg, Eg**2 * F / d0, color=col, lw=2.0)
        if i < 2:
            ax.text(Eg[k], Eg[k] ** 2 * F[k] / d0 * 1.35, f"d₀ = {d0:.0f} Gpc", fontsize=7.5,
                    color=col, ha="center")
        else:
            ax.text(3.0, 9e-5, f"d₀ = {d0:.0f} Gpc:\nthe bound", fontsize=7.5, color=col, ha="center")
    full = d["spectra_1gpc"]["full"]
    ax.plot(Eg, Eg**2 * full / m["headline_d0_min_gpc"], color=ps.SEQ_BLUE[7], lw=1.0, ls="--",
            label=f"d₀ = {m['headline_d0_min_gpc']:.0f} Gpc with full electron heating")
    g = d["cdg_ref"]["spec"]["computed"]
    ax.plot(g["E_MeV"], g["E_MeV"] ** 2 * g["flux"] * 20 / 1000, color=ps.INK, lw=1.0, ls=":",
            label="Cohen et al. 1998's prediction for 1 Gpc (their universe)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.012, 2000)
    ax.set_ylim(1e-5, 30)
    ax.set_xlabel("photon energy today (MeV)")
    ax.set_ylabel("E² × intensity  (MeV cm⁻² s⁻¹ sr⁻¹)")
    ax.set_title("a  Predicted annihilation glow against the measured sky")
    ax.legend(fontsize=6.8, loc="upper right")
    ax.text(0.02, 0.97, "blue: annihilations between z = 1100 and 19, Planck universe,\n"
            "conservative heating, 3.93 photons per annihilation", transform=ax.transAxes,
            fontsize=7, color=ps.MUTED, va="top")

    # (b) how robust is the bound?
    ax = fig.add_axes([0.62, 0.14, 0.16, 0.72])
    L = m["limits"]
    rows = [("conservative heating, COMPTEL (W99)", L["cdg/weidenspointner"]["d0_min_gpc"], True),
            ("conservative heating, COMPTEL (K98)", L["cdg/kappadath"]["d0_min_gpc"], False),
            ("full electron heating, COMPTEL (W99)", L["full/weidenspointner"]["d0_min_gpc"], False),
            ("no heating at all, COMPTEL (W99)", L["none/weidenspointner"]["d0_min_gpc"], False),
            ("annihilation lasts to z = 9", m["y_window"]["y_min_10"], False),
            ("annihilation stops at z = 29", m["y_window"]["y_min_30"], False),
            ("annihilation stops at z = 49", m["y_window"]["y_min_50"], False),
            ("Cohen et al.'s 1998 curve, final COMPTEL", res["cdg_1998_vs_final_comptel_d0_min_gpc"], False)]
    for i, (lab, v, head) in enumerate(rows):
        ax.barh(-i, v, color=ps.C_MODEL if head else ps.RAMPS[ps.C_MODEL][1], height=0.6, lw=0)
        ax.text(v * 1.15, -i, f"{v:.0f}" if v >= 10 else f"{v:.1f}", va="center", fontsize=7.5)
    ax.axvline(m["comoving_distance_y1100_gpc"], color=ps.INK, lw=0.8, ls="--")
    ax.text(m["comoving_distance_y1100_gpc"] * 1.08, 0.75,
            "distance to the CMB (13.9 Gpc)", fontsize=6.8, color=ps.INK_SECONDARY)
    ax.set_ylim(-len(rows) + 0.4, 1.1)
    ax.set_xscale("log")
    ax.set_xlim(1, 2e4)
    ax.set_yticks(-np.arange(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("smallest allowed domain scale d₀ (Gpc)")
    ax.set_title("b  The bound under other assumptions", x=-1.6, ha="left")

    # (c) a single boundary
    ax = fig.add_axes([0.83, 0.14, 0.155, 0.72])
    planes = m["single_plane_excess"]
    r = np.array([float(kk) for kk in planes])
    v = np.array([planes[kk] for kk in planes])
    o = np.argsort(r)
    ax.plot(r[o], v[o], "o-", color=ps.C_MODEL, lw=1.8, ms=3.5)
    ax.axhline(1.0, color=ps.INK, lw=0.8)
    ax.axhspan(1.0, 1e3, color=ps.EXCLUDED, lw=0, zorder=0)
    ax.text(0.4, 1.5, "excluded", fontsize=7.5, color=ps.INK_SECONDARY)
    for dist, lab in ((m["comoving_distance_y20_gpc"], "z = 19"),
                      (m["comoving_distance_y1100_gpc"], "CMB")):
        ax.axvline(dist, color=ps.AXIS, lw=0.8, ls=":")
        ax.text(dist - 0.2, 1.2e-3, lab, fontsize=7, color=ps.MUTED, ha="right")
    ax.set_yscale("log")
    ax.set_ylim(3e-4, 30)
    ax.set_xlim(0, 14.2)
    ax.set_xlabel("distance of ONE boundary (Gpc)")
    ax.set_ylabel("its MeV glow / COMPTEL background")
    ax.set_title("c  Even one boundary shows")
    ps.caption(fig, f"To hide their annihilation glow under COMPTEL's MeV background, domains of matter and antimatter "
                    f"would need d₀ > {m['headline_d0_min_gpc']:.0f} Gpc (1/d₀ = boundary area per volume): far larger "
                    "than the 13.9 Gpc to the CMB, under every assumption tried except two extremes: no heating of the gas "
                    "by annihilation electrons, or annihilation ending already at z ≈ 49 (b). "
                    "A single boundary within ~11.8 Gpc of us would outshine the measured sky (c).",
               y=0.02)
    return ps.save(fig, cfg.FIGURE_DIR / "04_mev_sky")


# ---------------------------------------------------------------------------
# 5. Antiprotons in cosmic rays: made in collisions, not left over
# ---------------------------------------------------------------------------


def fig_antiprotons():
    import cosmic_rays as cr
    import fetch_data as fd

    r = _res()["cosmic_rays"]
    cur = load_pickle("cr_curves")
    fig = plt.figure(figsize=(15.5, 4.8))
    axs = [fig.add_axes([0.05, 0.15, 0.25, 0.7]), fig.add_axes([0.37, 0.15, 0.33, 0.7]),
           fig.add_axes([0.77, 0.15, 0.21, 0.7])]

    # (a) boron / carbon: the grammage clock
    ax = axs[0]
    bc = cr.load_ratio("b_over_c", cr.AMS_PERIOD)
    ax.errorbar(bc.e_mean, bc.value, yerr=bc.err, fmt="o", ms=3, color=ps.INK, elinewidth=0.7,
                label="AMS-02 2011–2018")
    ax.plot(cur["bc_model_R"], cur["bc_model"], color=ps.C_MODEL, lw=1.8, label="leaky box, fitted")
    ax.axvline(4.0, color=ps.AXIS, lw=0.8, ls=":")
    ax.text(4.3, 0.36, "fit range →", fontsize=7, color=ps.MUTED)
    f = r["bc_fit"]
    ax.text(0.97, 0.97, f"grammage X = {f['X0_10GV']:.2f} g/cm² × (R / 10 GV)^−{f['delta']:.3f}\n"
            f"(hardening above 312 GV)\nχ² = {f['chi2']:.0f} for {f['ndof']} points (correlated\n"
            "systematics counted as independent)", transform=ax.transAxes, ha="right", va="top",
            fontsize=7, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("rigidity R (GV)")
    ax.set_ylabel("boron / carbon")
    ax.set_title("a  Boron measures the gas cosmic rays cross")
    ax.legend(fontsize=7, loc="lower left")

    # (b) antiprotons / protons
    ax = axs[1]
    ekn = fd.load_crdb("pbar_over_p_ekn")
    for sub, mk, lab in (("PAMELA (2006/07-2009/12)", "v", "PAMELA 2006–2009"),
                         ("BESS-PolarII (2007/12-2008/01)", "<", "BESS-Polar II 2007–2008")):
        dd = ekn[(ekn.subexp == sub) & ~ekn.is_ul]
        R = np.sqrt(dd.e_mean * (dd.e_mean + 2 * cr.MP))
        ax.errorbar(R, dd.value, yerr=dd.err, fmt=mk, ms=3.5, color=ps.C_OTHER, elinewidth=0.6,
                    label=lab, mfc=ps.SURFACE)
    for key, mk, lab in (("2011-2015/II", "o", "AMS-02 2011–2015 (PRL 2016)"),
                         ("2011-2022/II", "s", "AMS-02 2011–2022 (PRL 2025)")):
        v = cur["preds"][key]
        ax.errorbar(v["data_R"], v["data"], yerr=v["data_err"], fmt=mk, ms=3.2, color=ps.INK,
                    elinewidth=0.7, mfc=ps.INK if mk == "o" else ps.SURFACE, label=lab, zorder=5)
    v2, v1 = cur["preds"]["2011-2015/II"], cur["preds"]["2011-2015/I"]
    band = math.hypot(0.12, 0.10)
    ax.fill_between(v2["R"], v2["model"] * (1 - band), v2["model"] * (1 + band),
                    color=ps.RAMPS[ps.C_MODEL][0], lw=0, alpha=0.8)
    ax.plot(v2["R"], v2["model"], color=ps.C_MODEL, lw=1.8,
            label="secondaries, same grammage (cross-section fit II)")
    ax.plot(v1["R"], v1["model"], "--", color=ps.C_MODEL, lw=1.4, label="… cross-section fit I")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.9, 600)
    ax.set_ylim(2e-6, 6e-4)
    ax.set_xlabel("rigidity R (GV)")
    ax.set_ylabel("antiprotons / protons")
    ax.set_title("b  Antiprotons: just what collisions make")
    s15 = r["pbar_over_p"]["2011-2015/II"]
    ax.text(0.97, 0.05, f"above 20 GV: model / data = {s15['mean_model_over_data_above_20GV']:.2f}\n"
            f"band: ±{100 * band:.0f}% (production and boron cross-sections)\n"
            "below ~10 GV this simple model has no tertiaries or energy losses",
            transform=ax.transAxes, ha="right", fontsize=7, color=ps.INK_SECONDARY)
    ax.legend(fontsize=6.6, loc="upper left")

    # (c) antihelium
    ax = axs[2]
    ah = fd.load_crdb("hebar_over_he")
    show = ["BESS-COMBINED (1993-2007)", "BESS-PolarII (2007/12-2008/01)", "PAMELA (2006/07-2009/12)",
            "AMS01 (1998/06)"]
    names = {"BESS-COMBINED (1993-2007)": "BESS 1993–2007", "BESS-PolarII (2007/12-2008/01)": "BESS-Polar II",
             "PAMELA (2006/07-2009/12)": "PAMELA", "AMS01 (1998/06)": "AMS-01 (1998)"}
    for i, sub in enumerate(show):
        row = ah[ah.subexp == sub].iloc[0]
        col = ps.C_DATA if i == 0 else ps.C_OTHER
        ax.hlines(row.value, row.e_low, row.e_high, color=col, lw=2.0)
        mid = math.sqrt(row.e_low * row.e_high)
        ax.annotate("", xy=(mid, row.value / 3), xytext=(mid, row.value),
                    arrowprops={"arrowstyle": "-|>", "color": col, "lw": 1.0})
        ax.text(row.e_high * 1.15, row.value, f"{names[sub]}: < {sci(row.value)}", fontsize=7,
                color=col, va="center")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 3e4)
    ax.set_ylim(1e-8, 1e-5)
    ax.set_xlabel("rigidity R (GV)")
    ax.set_ylabel("antihelium / helium  (upper limit)")
    ax.set_title("c  No antihelium")
    ax.text(0.03, 0.04, "an antistar or antigalaxy would shed\nantihelium; none has been confirmed",
            transform=ax.transAxes, fontsize=7, color=ps.MUTED)
    ps.caption(fig, "One grammage fitted to boron (made by carbon and oxygen breaking up on interstellar gas) predicts "
                    "the antiprotons made by cosmic-ray protons and helium hitting the same gas: above 20 GV it matches "
                    "AMS-02 to 2–4%, leaving no room for a primary antimatter source. Antihelium, which only antistars "
                    "could supply in quantity, has never been seen.", y=0.02)
    return ps.save(fig, cfg.FIGURE_DIR / "05_antiprotons")


# ---------------------------------------------------------------------------
# 6. Sakharov's three conditions, and why the Standard Model fails them
# ---------------------------------------------------------------------------


def fig_sakharov():
    import sakharov as sk

    r = _res()["sakharov"]
    cur = load_pickle("sakharov_curves")
    lat = lit.LATTICE_EW
    fig = plt.figure(figsize=(15.5, 4.7))
    axs = [fig.add_axes([0.05, 0.17, 0.24, 0.68]), fig.add_axes([0.505, 0.17, 0.17, 0.68]),
           fig.add_axes([0.75, 0.17, 0.24, 0.68])]

    # (1) baryon number violation: sphalerons against the expansion
    ax = axs[0]
    T = np.linspace(125, 175, 300)
    b = lat["sphaleron_broken"]
    broken = np.exp(b["a"] * T - b["b"])                    # Gamma / T^4
    sym = lat["sphaleron_symmetric"][0]
    rate = np.where(T < lat["T_c"][0], np.minimum(broken, sym), sym)
    ax.plot(T, rate * T, color=ps.C_MODEL, lw=2.0, label="sphaleron rate Γ/T³ (lattice)")
    H = sk.hubble_rate_gev(T, lat["freezeout_gstar"])
    ax.plot(T, lat["freezeout_alpha"] * H, color=ps.INK, lw=1.2, ls="--",
            label="expansion rate (× 0.1015, the freeze-out criterion)")
    ax.axvline(r["T_star_ours"], color=ps.AXIS, lw=0.8, ls=":")
    ax.text(r["T_star_ours"] + 0.6, 1e-10, f"sphalerons stop\nat {r['T_star_ours']:.1f} GeV\n"
            f"(paper: {lat['T_star'][0]} ± {lat['T_star'][1]})", fontsize=7.2, color=ps.INK_SECONDARY,
            ha="right")
    ax.axvspan(lat["T_c"][0] - lat["T_c"][1], lat["T_c"][0] + lat["T_c"][1], color=ps.EXCLUDED, lw=0)
    ax.text(lat["T_c"][0] + 2, 3e-4, "Higgs field\nswitches on\n(159.5 GeV)", fontsize=7, color=ps.MUTED)
    ax.set_yscale("log")
    ax.set_ylim(1e-20, 1)
    ax.set_xlim(125, 175)
    ax.invert_xaxis()
    ax.set_xlabel("temperature (GeV; time runs to the right →)")
    ax.set_ylabel("rate (GeV)")
    ax.set_title("1  Baryon number can change ✓")
    ax.legend(fontsize=6.8, loc="lower right")

    # (2) CP violation: present, but tiny at the electroweak temperature
    ax = axs[1]
    jm = r["jarlskog_measure"]
    rows = [("observed excess Y_B", r["Y_B_needed"], ps.INK),
            ("Jarlskog invariant J (PDG 2026)", r["J_pdg"], ps.C_OTHER),
            ("SM CP measure, T = 100 GeV", jm["pdg2026"]["T100"], ps.C_MODEL),
            (f"… T = {r['T_star_ours']:.1f} GeV (freeze-out)", list(jm["pdg2026"].values())[1], ps.C_MODEL),
            ("… T = 159.5 GeV (crossover)", jm["pdg2026"]["T159.5"], ps.C_MODEL),
            ("… same, running masses", jm["running_mz"]["T159.5"], ps.C_MODEL)]
    for i, (lab, v, col) in enumerate(rows):
        ax.barh(-i, math.log10(v) + 26, left=-26, color=col, height=0.6, lw=0)
        ax.text(math.log10(v) + 0.4, -i, sci(v), va="center", fontsize=7.5, color=ps.INK)
    ax.set_yticks(-np.arange(len(rows)))
    ax.set_yticklabels([x[0] for x in rows], fontsize=7)
    ax.set_xlim(-26, 2)
    ax.set_xlabel("log₁₀ of the dimensionless measure")
    ax.grid(axis="y", visible=False)
    ax.set_title("2  C and CP violation: ✗ 9–13 orders too weak", x=-0.75, ha="left")
    ax.text(-25.5, -len(rows) + 0.35, "measure = J × Π(m_i² − m_j²) / T¹²: CP violation needs all six quark\n"
            "masses to differ, and at ~100 GeV all are tiny compared with T", fontsize=6.8,
            color=ps.MUTED, va="top")
    ax.set_ylim(-len(rows) - 0.6, 0.6)

    # (3) departure from equilibrium: the strength of the transition
    ax = axs[2]
    ax.axhspan(*lit.WASHOUT["range"], color=ps.EXCLUDED, lw=0)
    ax.axhline(1.0, color=ps.INK, lw=0.8)
    ax.text(21, 1.03, "needed: v/T ≳ 1 (sphalerons must switch off)", fontsize=7, color=ps.INK_SECONDARY)
    ax.plot(cur["mh"], cur["vt_simple"], color=ps.C_MODEL, lw=1.0, ls="--", label="one loop: 2E/λ")
    ax.plot(cur["mh"], cur["vt"], color=ps.C_MODEL, lw=1.8, label="one loop with thermal masses (λ_T)")
    lv = np.array(lat["v_over_T"])
    ax.plot(lv[:, 0], lv[:, 1], "*", ms=10, color=ps.INK, label="lattice, full SM (Kajantie et al. 1996)")
    ax.axvspan(lat["m_h_endpoint"][0], 135, color=ps.RAMPS[ps.C_ALT][0], lw=0, alpha=0.6)
    ax.text(lat["m_h_endpoint"][0] + 2, 1.72, "no phase transition at all:\na smooth crossover\n(lattice: m_H > 72 GeV)",
            fontsize=7, color=ps.INK_SECONDARY, va="top")
    ax.axvline(125.2, color=ps.C_ALT, lw=2.0)
    ax.text(124, 0.45, "measured\nHiggs\n125.2 GeV", fontsize=7.5, color=ps.INK, ha="right")
    ax.set_xlim(20, 135)
    ax.set_ylim(0, 2.2)
    ax.set_xlabel("Higgs boson mass (GeV)")
    ax.set_ylabel("v(T_c) / T_c  (strength of the transition)")
    ax.set_title("3  Out of equilibrium: ✗ no transition")
    ax.legend(fontsize=6.6, loc="upper right")
    ps.caption(fig, "Sakharov (1967): a matter excess needs (1) baryon-number violation, (2) C and CP violation and "
                    "(3) departure from equilibrium. The Standard Model has the first (sphalerons, active above 132 GeV), "
                    "but its CP violation is 9–13 orders of magnitude too weak and, with a 125 GeV Higgs, the universe "
                    "never left equilibrium. Something beyond the Standard Model made the matter we are made of.",
               y=0.02)
    return ps.save(fig, cfg.FIGURE_DIR / "06_sakharov")


# ---------------------------------------------------------------------------
# 7. How much antimatter is there? Limits from planets to the universe
# ---------------------------------------------------------------------------


def fig_ladder():
    res = _res()
    cl = res["clusters"]
    dom = res["domains"]["modern"]
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    fig.subplots_adjust(left=0.33, right=0.97, bottom=0.16, top=0.9)
    rows = [
        # (scale label, limit, note, ours)
        ("Moon, planets (~100 AU)", None, "spacecraft landed, the solar wind would make antiplanets glow", False),
        ("Milky Way gas (10 kpc)", 1e-15, "Steigman 1976: Galactic gamma rays", False),
        ("antistars in the Milky Way", 2.5e-6, "Dupourqué et al. 2021: no antistar among Fermi sources", False),
        ("cosmic-ray antihelium (kpc)", 6.9e-8, "BESS 1993–2007: He̅/He below this", False),
        ("Virgo cluster (2 Mpc)", 5e-9, "Steigman 2008, EGRET", False),
        ("Perseus cluster (2 Mpc)", 8e-9, "Steigman 2008, EGRET", False),
        ("Bullet Cluster, colliding gas (1 Mpc)", 3.1e-6, "Steigman 2008, EGRET", False),
        ("Bullet Cluster, colliding gas (1 Mpc)", cl["f_limit_steigman_FX"],
         "this project: Fermi-LAT limit, thermal ⟨σv⟩, our photon yield", True),
    ]
    for i, (lab, f, note, ours) in enumerate(rows):
        y = -i
        if f is None:
            ax.text(1e-16, y, note, va="center", fontsize=7.5, color=ps.INK_SECONDARY)
            continue
        col = ps.C_MODEL if ours else ps.INK
        ax.hlines(y, f, 1.0, color=ps.EXCLUDED, lw=8, zorder=1)
        ax.plot(f, y, "*" if ours else "|", ms=13 if ours else 14, mew=2, color=col, mec=ps.INK if ours else col,
                zorder=3)
        txt = "≲ 10⁻¹⁵" if f == 1e-15 else f"< {sci(f)}"
        if f < 1e-13:
            ax.text(f * 3, y + 0.28, txt, va="center", ha="left", fontsize=8, color=col)
        else:
            ax.text(f / 2.5, y, txt, va="center", ha="right", fontsize=8, color=col)
        ax.text(1.4, y, note, va="center", fontsize=7, color=ps.MUTED)
    i = len(rows)
    ax.text(1e-16, -i, f"observable universe (14 Gpc to the CMB): no boundary between matter and antimatter;\n"
            f"domains would need d₀ > {dom['headline_d0_min_gpc']:.0f} Gpc (this project, MeV gamma rays)",
            va="center", fontsize=7.6, color=ps.C_MODEL)
    labels = [r[0] for r in rows] + ["the whole visible universe"]
    ax.set_yticks(-np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(1e-17, 1e6)
    ax.set_ylim(-len(labels) + 0.4, 0.6)
    ax.set_xticks([1e-16, 1e-12, 1e-8, 1e-4, 1])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("largest fraction of antimatter allowed (grey: excluded)")
    ax.set_title("Wherever we look, there is no antimatter", x=-0.4, ha="left")
    ps.caption(fig, "Annihilation of antimatter with ordinary matter would make gamma rays (or antinuclei) that are "
                    "not seen, on every scale from the Solar System to the observable universe. Blue: computed in "
                    "this project.", y=0.02)
    return ps.save(fig, cfg.FIGURE_DIR / "07_antimatter_ladder")


FIGURES = {1: fig_eta_two_ways, 2: fig_asymmetry, 3: fig_boundary, 4: fig_mev_sky, 5: fig_antiprotons,
           6: fig_sakharov, 7: fig_ladder}


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ps.use_style()
    which = [int(a) for a in argv] if argv else sorted(FIGURES)
    for k in which:
        out = FIGURES[k]()
        print("   ", out[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
