"""Every figure of the project (PNG + SVG in figures/), from the cached analysis."""

from __future__ import annotations

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import config as cfg  # noqa: E402
import plotstyle as ps  # noqa: E402

D = cfg.DERIVED_DIR


def load(name: str) -> dict:
    return json.loads((D / f"{name}.json").read_text())


def out(name: str) -> str:
    return str(cfg.FIGURE_DIR / name)


def lab(key: str) -> str:
    import datasets
    return datasets.ALICE_LABELS[key][0]


# ---------------------------------------------------------------------------
# 1. Yields: data against the fitted hadron gas
# ---------------------------------------------------------------------------


def fig_yields() -> None:
    m = load("main_fit")
    ys = [y for y in m["yields"] if y["key"] != "Lambdabar"]
    p = m["params"]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9.2, 6.2), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.25], "hspace": 0.06})
    x = np.arange(len(ys))
    for i, y in enumerate(ys):
        ax.plot([i - 0.32, i + 0.32], [y["model"]] * 2, color=ps.C_MODEL, lw=2.2,
                solid_capstyle="round", zorder=2)
        fitted = y["fitted"]
        ax.errorbar(i, y["value"], yerr=y["err"], fmt="o" if fitted else "o", ms=5.5,
                    mfc=ps.INK if fitted else ps.SURFACE, mec=ps.INK, ecolor=ps.INK,
                    elinewidth=1, capsize=0, zorder=3)
        r, e = y["ratio"], y["err"] / y["model"]
        ax2.errorbar(i, r, yerr=e, fmt="o", ms=4.5, mfc=ps.INK if fitted else ps.SURFACE,
                     mec=ps.INK, ecolor=ps.INK, elinewidth=1, capsize=0, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 3e3)
    ax.set_ylabel("dN/dy  (particles per unit rapidity)")
    ax.set_title("Twenty-one hadron and nucleus yields from one temperature")
    T = p["T"]["value"]
    txt = (f"ALICE, Pb–Pb √s$_{{NN}}$ = 2.76 TeV, 0–10% central\n"
           f"Fit: T = {T:.1f} ± {p['T']['error']:.1f} MeV,  "
           f"μ$_B$ = {p['muB']['value']:.1f} ± {p['muB']['error']:.1f} MeV,\n"
           f"V = {p['V']['value']:.0f} ± {p['V']['error']:.0f} fm³,  "
           f"χ² = {m['chi2']:.1f} for {m['ndf']} degrees of freedom")
    ax.text(0.98, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8,
            color=ps.INK_SECONDARY)
    handles = [Line2D([], [], color=ps.C_MODEL, lw=2.2, label="hadron resonance gas (fit)"),
               Line2D([], [], marker="o", ls="", mfc=ps.INK, mec=ps.INK, label="measured, fitted"),
               Line2D([], [], marker="o", ls="", mfc=ps.SURFACE, mec=ps.INK,
                      label="measured, not fitted (K*: rescattered)")]
    ax.legend(handles=handles, loc="lower left")
    ax2.axhline(1, color=ps.C_MODEL, lw=1.2)
    ax2.set_ylim(0.3, 2.3)
    ax2.set_ylabel("data / model")
    ax2.set_xticks(x)
    ax2.set_xticklabels([lab(y["key"]) for y in ys], fontsize=8.5)
    ax2.tick_params(axis="x", colors=ps.INK_SECONDARY)
    ip = [y["key"] for y in ys].index("p")
    ax2.text(ip + 0.5, 0.45, "proton anomaly", fontsize=7, color=ps.INK_SECONDARY, ha="center")
    ps.caption(fig, "Data: HEPData and the ALICE papers (see README). Model: every PDG hadron "
               "(PDG 2025 list), decayed with ALICE's definition of primary particles; "
               "errors are statistical and systematic in quadrature.", y=0.02)
    ps.save(fig, out("01_yields"))


# ---------------------------------------------------------------------------
# 2. Yield per spin state against mass
# ---------------------------------------------------------------------------


def fig_mass() -> None:
    m = load("main_fit")
    T = m["params"]["T"]["value"]
    keys = ("pi+", "K+", "phi", "p", "Lambda", "Xi-", "Omega-", "d", "He3", "H3L", "He4")
    ys = [y for y in m["yields"] if y["key"] in keys]
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    for y in ys:
        g = y["g"] if y["key"] != "K0S" else 1.0
        ax.errorbar(y["mass"], y["value"] / g, yerr=y["err"] / g, fmt="o", ms=5.5,
                    color=ps.INK, elinewidth=1, zorder=3)
        ax.plot(y["mass"], y["model"] / g, "_", ms=16, mew=2.2, color=ps.C_MODEL, zorder=2)
        ax.plot(y["mass"], y["direct"] / g, "_", ms=16, mew=2.2,
                color=ps.lighten(ps.C_MODEL, 0.55), zorder=2)
        off = (-14, -12) if y["key"] == "phi" else (6, 4)
        ax.annotate(lab(y["key"]), (y["mass"], y["value"] / g), off, textcoords="offset points",
                    fontsize=8, color=ps.INK_SECONDARY)
    # a pure thermal (Boltzmann) line through the nuclei: n/g ~ m^2 T K2(m/T)
    from scipy.special import kn
    mm = np.linspace(0.9, 4.0, 200)
    nuc = [y for y in ys if y["key"] in ("d", "He3", "He4")]
    Vn = np.mean([y["direct"] / y["g"] / (y["mass"] ** 2 * kn(2, y["mass"] / (T / 1e3)))
                  for y in nuc])
    ax.plot(mm, Vn * mm ** 2 * kn(2, mm / (T / 1e3)), color=ps.C_MODEL, lw=1, ls="--", zorder=1)
    ax.text(2.2, Vn * 2.2 ** 2 * kn(2, 2.2 / (T / 1e3)) * 4, "∝ m² T K₂(m/T): a factor\n"
            "≈ 300 per extra nucleon", fontsize=7.5, color=ps.C_MODEL)
    ax.set_yscale("log")
    ax.set_xlabel("mass (GeV)")
    ax.set_ylabel("dN/dy per spin state, (2J + 1)⁻¹ dN/dy")
    ax.set_title(f"Heavier particles are rarer, at the rate a {T:.0f} MeV gas sets")
    handles = [Line2D([], [], marker="o", ls="", color=ps.INK, label="measured (particles)"),
               Line2D([], [], marker="_", ms=14, mew=2.2, ls="", color=ps.C_MODEL,
                      label="model, after resonance decays"),
               Line2D([], [], marker="_", ms=14, mew=2.2, ls="", color=ps.lighten(ps.C_MODEL, .55),
                      label="model, born directly (primordial)")]
    ax.legend(handles=handles, loc="lower left")
    ps.caption(fig, "The gap between the two model marks is the share of each species made in "
               "decays of heavier resonances (2/3 of pions, 3/4 of Λ); nuclei have none.")
    ps.save(fig, out("02_mass_dependence"))


# ---------------------------------------------------------------------------
# 3. The temperature: profile and every variant
# ---------------------------------------------------------------------------


def fig_temperature() -> None:
    prof, var, m = load("profile"), load("variants"), load("main_fit")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.8),
                                 gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.5})
    T, c2 = np.array(prof["T"]), np.array(prof["chi2"])
    a1.plot(T, c2 - c2.min(), color=ps.C_MODEL)
    a1.axhline(1, color=ps.MUTED, lw=0.8, ls="--")
    a1.text(T[1], 1.3, "Δχ² = 1 (68%)", fontsize=7.5, color=ps.INK_SECONDARY)
    hq = cfg.HOTQCD_2019
    a1.axvspan(hq["Tpc"] - hq["Tpc_err"], hq["Tpc"] + hq["Tpc_err"], color=ps.WASH[ps.VIOLET])
    a1.text(hq["Tpc"] + 0.4, 18, "lattice\nT$_{pc}$", fontsize=7.5, color=ps.C_LATTICE)
    a1.set_ylim(0, 25)
    a1.set_xlim(T.min(), T.max())
    a1.set_xlabel("chemical freeze-out temperature T (MeV)")
    a1.set_ylabel("χ² − χ²$_{min}$ (μ$_B$, V refitted)")
    a1.set_title("How sharply the yields fix T")
    rows = [(k, v) for k, v in var.items() if not k.startswith("_")]
    short = {"default": "default fit", "no_charm": "no charm feed-down",
             "charm5": "charm feed-down × 0.5", "charm15": "charm feed-down × 1.5",
             "zero_width": "zero-width resonances", "pdg2014": "PDG 2014 hadron list",
             "boltzmann": "Boltzmann statistics", "gammaS": "γ$_S$ free",
             "lambdabar": "Λ̄ added (as Andronic)", "andronic_data": "Andronic's data compilation",
             "no_protons": "without p, p̄", "hadrons_only": "hadrons only",
             "nuclei_only": "nuclei only", "d_he3_he4": "d, ³He, ⁴He (μ$_B$ = 0)",
             "p_d_he3": "p, d, ³He (μ$_B$ = 0)"}
    labels = [short.get(k, k) for k, _ in rows]
    pub = [("Andronic et al. 2018 (Nature)", cfg.ANDRONIC_2018["T"], cfg.ANDRONIC_2018["T_err"]),
           ("ALICE 2016: p, d, ³He", cfg.ALICE_NUCLEI_2016["T_p_d_3He"], cfg.ALICE_NUCLEI_2016["T_err"]),
           ("ALICE 2018: d, ³He, ⁴He", cfg.ALICE_HE4_2018["T_nuclei"], cfg.ALICE_HE4_2018["T_err"]),
           ("Andronic: nuclei only", cfg.ANDRONIC_2018["T_nuclei"], cfg.ANDRONIC_2018["T_nuclei_err"])]
    n = len(rows) + len(pub)
    yv = np.arange(n)[::-1]
    a2.axvspan(hq["Tpc"] - hq["Tpc_err"], hq["Tpc"] + hq["Tpc_err"], color=ps.WASH[ps.VIOLET])
    T0 = var["default"]["T"]
    a2.axvline(T0, color=ps.MUTED, lw=0.7, ls=":")
    for i, (k, v) in enumerate(rows):
        a2.errorbar(v["T"], yv[i], xerr=v["T_err"], fmt="*" if k == "default" else "o",
                    ms=9 if k == "default" else 4.5, color=ps.INK, elinewidth=1)
    for j, (l, t, e) in enumerate(pub):
        a2.errorbar(t, yv[len(rows) + j], xerr=e, fmt="s", ms=4.5, color=ps.C_ANDRONIC,
                    elinewidth=1)
    a2.set_yticks(yv)
    a2.set_yticklabels(labels + [p[0] for p in pub], fontsize=7.4)
    a2.tick_params(axis="y", colors=ps.INK_SECONDARY, length=0)
    a2.set_xlabel("T (MeV)")
    a2.set_title("Every way of doing the fit, and the published values")
    a2.text(hq["Tpc"], yv[0] + 1.0, "lattice crossover 156.5 ± 1.5", fontsize=7.2,
            color=ps.C_LATTICE, ha="center")
    a2.set_ylim(-1, n + 0.6)
    ps.caption(fig, "Ink: our fits (star = default). Orange: published thermal fits. Violet band: "
               "HotQCD 2019's chiral crossover temperature.", y=-0.03)
    ps.save(fig, out("03_temperature"))


# ---------------------------------------------------------------------------
# 4. The QCD phase diagram
# ---------------------------------------------------------------------------


def fig_phase() -> None:
    ph = load("phase_diagram")
    fig, (ax, axl) = plt.subplots(1, 2, figsize=(11.5, 5.4),
                                  gridspec_kw={"width_ratios": [1.45, 1], "wspace": 0.28})
    lat = ph["lattice"]
    hq, bz = lat["hotqcd"], lat["borsanyi"]
    cut = lambda d: {k: np.array(v)[np.array(d["mu_B"]) <= 300] for k, v in d.items()}  # noqa: E731
    hq, bz = cut(hq), cut(bz)
    ax.fill_between(hq["mu_B"], hq["lo"], hq["hi"], color=ps.WASH[ps.VIOLET], lw=0, zorder=1)
    ax.plot(hq["mu_B"], hq["T"], color=ps.C_LATTICE, lw=1.4, zorder=2)
    ax.plot(bz["mu_B"], bz["T"], color=ps.C_LATTICE, lw=1.2, ls="--", zorder=2)
    ax.text(310, 160, "lattice QCD crossover\n(HotQCD 2019 band; dashed:\nBorsányi et al. 2020)",
            fontsize=7.5, color=ps.C_LATTICE)
    # published points
    for g, mk, col, lab in (("other_a", "^", ps.C_OTHER, None), ("other_b", "D", ps.C_OTHER, None),
                            ("other_c", "v", ps.C_OTHER, None),
                            ("andronic", "s", ps.C_ANDRONIC, "Andronic et al. 2018")):
        pts = ph["published"].get(g, [])
        if not pts:
            continue
        ax.errorbar([p["muB"] for p in pts], [p["T"] for p in pts],
                    yerr=[p["T_err"] for p in pts], fmt=mk, ms=4.5, color=col, mfc=ps.SURFACE,
                    elinewidth=0.8, zorder=3, label=lab)
    sh = ph["star_hepdata"]
    ax.errorbar([v["muB"] for v in sh.values()], [v["T"] for v in sh.values()],
                xerr=[v["muB_err"] for v in sh.values()], yerr=[v["T_err"] for v in sh.values()],
                fmt="h", ms=5, color=ps.C_STAR, mfc=ps.SURFACE, elinewidth=0.8, zorder=3,
                label="STAR 2017 (published)")
    ours = ph["ours"]
    ax.errorbar([o["muB"] for o in ours], [o["T"] for o in ours],
                xerr=[o["muB_err"] for o in ours], yerr=[o["T_err"] for o in ours],
                fmt="*", ms=11, color=ps.INK, elinewidth=1, zorder=5, label="this project's fits")
    ax.text(420, 178, "★ left to right: LHC 2.76 TeV; RHIC 39, 27,\n19.6, 11.5, 7.7 GeV "
            "(STAR data, our fits)", fontsize=7, color=ps.INK_SECONDARY)
    # freeze-out at constant energy per hadron
    en = np.array(ph["e_over_n_line"])
    ax.plot(en[:, 0], en[:, 1], color=ps.MUTED, lw=1, ls=":", zorder=1)
    ax.text(en[-1, 0] - 30, en[-1, 1] - 14, f"⟨E⟩/⟨N⟩ = {ph['e_over_n_GeV']:.2f} GeV",
            fontsize=7, color=ps.INK_SECONDARY)
    # the universe
    up = ph["universe_path"]
    ax.plot(up["muB_MeV"], up["T_MeV"], color=ps.C_UNIVERSE, lw=2.2, zorder=4)
    ax.annotate("the early universe\n(μ$_B$ ≈ 2 eV at 156 MeV;\nantimatter gone below 40 MeV)",
                (3, 118), (60, 98), fontsize=7.5, color=ps.INK_SECONDARY,
                arrowprops={"arrowstyle": "-", "color": ps.C_UNIVERSE, "lw": 0.9})
    ax.plot(ph["nuclear_matter_muB"], 0, "v", ms=8, color=ps.INK, clip_on=False, zorder=5)
    ax.text(ph["nuclear_matter_muB"] - 10, 6, "atomic nuclei\n(923 MeV)", fontsize=7,
            ha="right", color=ps.INK_SECONDARY)
    ax.text(30, 186, "quark–gluon plasma", fontsize=10, color=ps.INK_SECONDARY, style="italic")
    ax.text(330, 40, "hadron gas", fontsize=10, color=ps.INK_SECONDARY, style="italic")
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 200)
    ax.set_xlabel("baryon chemical potential μ$_B$ (MeV)")
    ax.set_ylabel("temperature T (MeV)")
    ax.set_title("The QCD phase diagram: where hadrons froze out")
    ax.legend(loc="lower left", fontsize=7.2, handlelength=1.5)
    # right: log mu_B/T, universe vs experiments
    axl.plot(up["muB_over_T"], up["T_MeV"], color=ps.C_UNIVERSE, lw=2.2, label="early universe")
    axl.errorbar([max(o["muB"], 1e-9) / o["T"] for o in ours], [o["T"] for o in ours],
                 fmt="*", ms=10, color=ps.INK, label="heavy-ion collisions (ours)")
    lhc = ours[0]
    axl.errorbar(lhc["muB_err"] / lhc["T"], lhc["T"], fmt="none")
    axl.annotate("LHC: μ$_B$/T = %.3f ± %.3f\n(consistent with zero)" % (
        lhc["muB"] / lhc["T"], lhc["muB_err"] / lhc["T"]), (lhc["muB"] / lhc["T"], lhc["T"]),
        (1e-6, 172), fontsize=7.2, color=ps.INK_SECONDARY,
        arrowprops={"arrowstyle": "-", "color": ps.MUTED, "lw": 0.6})
    axl.set_xscale("log")
    axl.set_xlim(1e-9, 1e3)
    axl.set_ylim(0, 260)
    axl.axhspan(cfg.HOTQCD_2019["Tpc"] - 1.5, cfg.HOTQCD_2019["Tpc"] + 1.5,
                color=ps.WASH[ps.VIOLET])
    axl.set_xlabel("μ$_B$ / T  (logarithmic)")
    axl.set_ylabel("T (MeV)")
    axl.set_title("Matter excess: universe vs collider")
    axl.text(2e-8, 205, "μ$_B$/T ≈ 10⁻⁸: one extra\nquark per ~10⁸ quark–\nantiquark pairs",
             fontsize=7.2, color=ps.INK_SECONDARY)
    axl.legend(loc="lower right", fontsize=7.2)
    ps.caption(fig, "Published points read from Andronic et al. 2018 Fig. 4 (orange: theirs; grey: "
               "other analyses they compiled) and STAR's HEPData record. Universe: Planck's baryon "
               "density, lepton asymmetry zero.", y=-0.02)
    ps.save(fig, out("04_phase_diagram"))


# ---------------------------------------------------------------------------
# 5. Equation of state
# ---------------------------------------------------------------------------


def fig_eos() -> None:
    L = load("lattice")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"wspace": 0.25})
    hq, wb = L["hotqcd2014"], L["wb2014"]
    h = L["hrg"]
    hT = [r["T"] for r in h]
    for q, lab in (("e", "ε/T⁴"), ("p", "3p/T⁴")):
        f = 3.0 if q == "p" else 1.0
        T = np.array(hq["T"])
        v, e = np.array(hq[q]) * f, np.array(hq[q + "_err"]) * f
        a1.fill_between(T, v - e, v + e, color=ps.WASH[ps.VIOLET], lw=0)
        a1.plot(T, v, color=ps.C_LATTICE, lw=1.3)
        Tw = np.array(wb["T"])
        a1.plot(Tw, np.array(wb[q]) * f, color=ps.C_LATTICE, lw=1, ls="--")
        a1.plot(hT, np.array([r[q] for r in h]) * f, color=ps.C_MODEL, lw=1.8)
        a1.text(395, v[-1] + 0.25, lab, fontsize=8, color=ps.INK_SECONDARY, ha="right")
    sb = L["sb_p"]
    a1.axhline(3 * sb, color=ps.MUTED, lw=0.9, ls=":")
    a1.text(105, 3 * sb + 0.3, "free quarks and gluons: ε = 3p = 15.6 T⁴", fontsize=7.5,
            color=ps.INK_SECONDARY)
    a1.axvspan(155, 158, color=ps.WASH[ps.ORANGE])
    a1.set_xlim(100, 400)
    a1.set_ylim(0, 17)
    a1.set_xlabel("T (MeV)")
    a1.set_ylabel("energy density and 3 × pressure over T⁴")
    a1.set_title("Hadron gas below T$_{pc}$, deconfined quarks above")
    handles = [Line2D([], [], color=ps.C_LATTICE, lw=1.3, label="lattice QCD, HotQCD 2014 (band)"),
               Line2D([], [], color=ps.C_LATTICE, lw=1, ls="--", label="lattice QCD, Borsányi et al. 2014"),
               Line2D([], [], color=ps.C_MODEL, lw=1.8, label="our hadron resonance gas")]
    a1.legend(handles=handles, loc="lower right", fontsize=7.2)
    s = L["hotqcd2012_chi"]
    T = np.array(s["T"])
    a2.fill_between(T, np.array(s["BB"]) - s["BB_err"], np.array(s["BB"]) + s["BB_err"],
                    color=ps.WASH[ps.VIOLET], lw=0)
    a2.plot(T, s["BB"], color=ps.C_LATTICE, lw=1.3, label="lattice QCD (HotQCD 2012)")
    a2.plot(hT, [r["chiBB"] for r in h], color=ps.C_MODEL, lw=1.8, label="hadron resonance gas")
    a2.axhline(1 / 3, color=ps.MUTED, lw=0.9, ls=":")
    a2.text(102, 1 / 3 + 0.008, "free quarks: 1/3", fontsize=7.5, color=ps.INK_SECONDARY)
    a2.axvspan(155, 158, color=ps.WASH[ps.ORANGE])
    a2.text(156.5, 0.36, "T$_{CF}$", fontsize=7.5, color=ps.C_ANDRONIC, ha="center")
    a2.set_xlim(100, 250)
    a2.set_ylim(0, 0.4)
    a2.set_xlabel("T (MeV)")
    a2.set_ylabel("baryon-number susceptibility χ₂$^B$")
    a2.set_title("Baryon number: carried by 3-quark hadrons, then by free quarks")
    a2.legend(loc="lower right", fontsize=7.2)
    ps.caption(fig, "Orange band: our chemical freeze-out temperature. The hadron gas follows the "
               "lattice up to the crossover and then keeps rising: hadrons melt.", y=-0.02)
    ps.save(fig, out("05_equation_of_state"))


# ---------------------------------------------------------------------------
# 6. The universe through the crossover
# ---------------------------------------------------------------------------


def fig_universe() -> None:
    import universe as U
    un, m = load("universe"), load("main_fit")
    cal = U.p2_calendar()
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    T = np.logspace(np.log10(5), np.log10(2000), 300)
    t = U.age_at(T, cal)
    ax.plot(t * 1e6, T, color=ps.INK, lw=1.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time since the Big Bang (microseconds)")
    ax.set_ylabel("temperature (MeV)")
    ax.set_title("The universe cooled through our freeze-out temperature 20 µs after the Big Bang")
    ev = [("direct-photon slope ≈ 300 MeV\n(the collision's hottest light)", un["ages"]["T_eff photons (ours)"]),
          ("quarks bind into hadrons:\nT$_{CF}$ = %.1f MeV" % m["params"]["T"]["value"],
           un["ages"]["T_CF (our fit)"])]
    for lab, a in ev:
        age = a["age_s_P2"] or a["age_s_own"]
        ax.plot(age * 1e6, a["T_MeV"], "o", ms=6, color=ps.C_ANDRONIC, zorder=3)
        ax.annotate(lab, (age * 1e6, a["T_MeV"]), (12, 6), textcoords="offset points",
                    fontsize=7.5, color=ps.INK_SECONDARY)
    Ta, ta = un["T_antibaryons_gone_MeV"], un["age_antibaryons_gone_s"]
    ax.plot(ta * 1e6, Ta, "o", ms=6, color=ps.C_UNIVERSE, zorder=3)
    ax.annotate("antibaryons annihilated\n(μ$_B$ = T): %.0f MeV" % Ta, (ta * 1e6, Ta), (12, 6),
                textcoords="offset points", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.text(0.02, 0.04, "Calendar from project P2 (full Standard-Model equation of state).",
            transform=ax.transAxes, fontsize=7, color=ps.MUTED)
    ps.save(fig, out("06_universe"))


# ---------------------------------------------------------------------------
# 7. Centrality dependence
# ---------------------------------------------------------------------------


def fig_centrality() -> None:
    c = load("centrality")["classes"]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), gridspec_kw={"wspace": 0.32})
    npart = np.array([r["npart"] for r in c])
    for k, mk, col, lab in (("eq", "o", ps.INK, "equilibrium (γ$_S$ = 1)"),
                            ("gammaS", "s", ps.C_MODEL, "γ$_S$ free")):
        axs[0].errorbar(npart, [r[k]["T"] for r in c], yerr=[r[k]["T_err"] for r in c], fmt=mk,
                        ms=5, color=col, mfc=col if k == "eq" else ps.SURFACE, label=lab)
    axs[0].axhspan(155, 158, color=ps.WASH[ps.VIOLET])
    axs[0].set_ylabel("T (MeV)")
    axs[0].set_title("T does not depend on centrality")
    axs[0].legend(fontsize=7.2, loc="lower right")
    axs[1].errorbar(npart, [r["gammaS"]["gammaS"] for r in c],
                    yerr=[r["gammaS"]["gammaS_err"] for r in c], fmt="s", ms=5, color=ps.C_MODEL,
                    mfc=ps.SURFACE)
    axs[1].axhline(1, color=ps.MUTED, lw=0.8, ls=":")
    axs[1].set_ylabel("strangeness saturation γ$_S$")
    axs[1].set_title("Strangeness: at or slightly above equilibrium")
    for key, mk, col, lab in (("p", "o", ps.INK, "p"), ("Omega-", "D", ps.C_ANDRONIC, "Ω⁻"),
                              ("K-", "^", ps.C_MODEL, "K⁻")):
        axs[2].errorbar(npart, [r["eq"]["ratios"][key] for r in c],
                        yerr=[r["eq"]["errs"][key] for r in c], fmt=mk, ms=5, color=col, label=lab)
    axs[2].axhline(1, color=ps.MUTED, lw=0.8, ls=":")
    axs[2].set_ylabel("data / model (equilibrium fit)")
    axs[2].set_title("The proton deficit is everywhere")
    axs[2].legend(fontsize=7.2)
    for a in axs:
        a.set_xlabel("participating nucleons ⟨N$_{part}$⟩")
    ps.caption(fig, "ALICE Pb–Pb 2.76 TeV in five centrality classes (0–10% at right, 60–80% at "
               "left); each class fitted separately.", y=-0.03)
    ps.save(fig, out("07_centrality"))


# ---------------------------------------------------------------------------
# 8. Direct photons
# ---------------------------------------------------------------------------


def fig_photons() -> None:
    import photons as P
    ph = load("photons")
    sp = P.spectrum()
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ok = np.isfinite(sp["y"])
    ax.errorbar(sp["pt"][ok], sp["y"][ok], yerr=np.hypot(sp["stat"][ok], sp["sys"][ok]),
                fmt="o", ms=4.5, color=ps.INK, elinewidth=0.9, label="ALICE direct photons, 0–20%")
    x = np.linspace(0.9, 2.1, 50)
    ax.plot(x, ph["A"] * np.exp(-x / (ph["T_eff_MeV"] / 1e3)), color=ps.C_MODEL, lw=2,
            label=f"exp(−p$_T$/T$_{{eff}}$), T$_{{eff}}$ = {ph['T_eff_MeV']:.0f} ± "
                  f"{ph['stat_MeV']:.0f} (stat) ± {ph['syst_MeV']:.0f} (syst) MeV")
    xx = np.linspace(0.9, 2.6, 50)
    Tcf = load("main_fit")["params"]["T"]["value"]
    ax.plot(xx, ph["A"] * np.exp(-0.9 / (ph["T_eff_MeV"] / 1e3)) * np.exp(-(xx - 0.9) / (Tcf / 1e3)),
            color=ps.MUTED, lw=1, ls="--", label=f"slope of a {Tcf:.0f} MeV source")
    ax.set_yscale("log")
    ax.set_xlim(0.8, 6)
    ax.set_ylim(1e-4, 8)
    ax.set_xlabel("photon transverse momentum p$_T$ (GeV/c)")
    ax.set_ylabel("1/(2π p$_T$) d²N/(dp$_T$ dy)  (GeV/c)⁻²")
    ax.set_title("Light from the plasma: twice as hot as freeze-out")
    ax.legend(fontsize=7.2, loc="upper right")
    ps.caption(fig, "HEPData ins1394677 Table 1. T_eff averages over the whole evolution and is "
               "blue-shifted by the expansion; hydrodynamic models need 300–600 MeV at the start.")
    ps.save(fig, out("08_photons"))


# ---------------------------------------------------------------------------
# 9. Validation: code against code
# ---------------------------------------------------------------------------


def fig_validation() -> None:
    an, st = load("andronic_comparison"), load("star")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.6, 1],
                                                                       "wspace": 0.25})
    keys = list(an["fig1"])
    x = np.arange(len(keys))
    for lab, mk, col in (("default", "o", ps.INK), ("no charm", "s", ps.C_MODEL),
                         ("zero width", "^", ps.C_OTHER)):
        r = an["variants"][lab]["ratio"]
        a1.plot(x, [r[k] for k in keys], mk, ms=5, color=col,
                mfc=col if lab == "default" else ps.SURFACE, label=lab)
    a1.axhline(1, color=ps.MUTED, lw=0.8)
    a1.axhspan(0.99, 1.01, color=ps.WASH[ps.NEUTRAL])
    import datasets as DS
    a1.set_xticks(x)
    a1.set_xticklabels([DS.ALICE_LABELS[k][0] for k in keys], fontsize=7.5, rotation=0)
    a1.text(len(keys) - 0.5, 0.82, "default = no charm except for π, K", fontsize=7,
            color=ps.INK_SECONDARY, ha="right")
    a1.set_ylabel("our model / Andronic et al.'s model")
    a1.set_title("Same T and μ$_B$, two codes: equal where no decays feed in")
    a1.legend(fontsize=7.2, loc="lower left")
    a1.set_ylim(0.8, 1.12)
    pub = st["published"]
    for e, mk in (("7.7", "o"), ("39.0", "s")):
        r = st["replica"][e]
        p = pub[e]
        a2.errorbar(p["muB"], p["T"], xerr=p["muB_err"], yerr=p["T_err"], fmt=mk, ms=6,
                    color=ps.C_STAR, mfc=ps.SURFACE, label=f"STAR {float(e):g} GeV (THERMUS)")
        a2.errorbar(r["muB"], r["T"], xerr=r["muB_err"], yerr=r["T_err"], fmt=mk, ms=5,
                    color=ps.INK, label=f"our code, STAR's inputs and choices")
    a2.set_xlabel("μ$_B$ (MeV)")
    a2.set_ylabel("T (MeV)")
    a2.set_title("STAR's own fits reproduced")
    h, l_ = a2.get_legend_handles_labels()
    a2.legend(h[:3], l_[:3], fontsize=7, loc="lower left")
    ps.caption(fig, "Left: models evaluated at T = 156.5 MeV, μ_B = 0.7 MeV, volume matched on the "
               "nuclei (read from Andronic et al.'s Fig. 1). Right: STAR 2017 GCE fits (μ_Q = 0, μ_S "
               "and γ_S free, inclusive protons).", y=-0.03)
    ps.save(fig, out("09_validation"))


ALL = {"yields": fig_yields, "mass": fig_mass, "temperature": fig_temperature,
       "phase": fig_phase, "eos": fig_eos, "universe": fig_universe,
       "centrality": fig_centrality, "photons": fig_photons, "validation": fig_validation}


def main(argv=None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    ps.use_style()
    for name, fn in ALL.items():
        if argv and name not in argv:
            continue
        fn()
        print(f"  figure {name}")


if __name__ == "__main__":
    main()
