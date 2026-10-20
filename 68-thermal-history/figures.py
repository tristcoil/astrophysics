"""Every figure, from data/results.json and data/derived/curves.npz.

    python figures.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, FuncFormatter, LogLocator, NullFormatter  # noqa: E402

import config as C  # noqa: E402
import plotstyle as S  # noqa: E402

YEAR = C.YEAR_S


def _load():
    res = json.loads(C.RESULTS_PATH.read_text())
    d = dict(np.load(C.DERIVED_DIR / "curves.npz", allow_pickle=False))
    return res, d


def _energy_formatter(x, _pos=None):
    return S.human_energy(x)


def _T_of_t(d, t):
    """Photon temperature (GeV) at age t (s), from our calendar."""
    lt, lT = np.log(d["cal_t"][::-1]), np.log(d["cal_T"][::-1])
    return np.exp(np.interp(np.log(np.asarray(t, dtype=float)), lt, lT))


def _time_axis(ax, d, *, ticks_s=None, label="Age of the universe"):
    """A top axis labelled in time, for an axis in temperature: ticks placed at T(t)."""
    top = ax.twiny()
    top.set_xscale(ax.get_xscale())
    top.set_xlim(ax.get_xlim())
    if ticks_s is None:
        ticks_s = [1e-12, 1e-9, 1e-6, 1e-3, 1, 60, 3600]
    lo, hi = sorted(ax.get_xlim())
    pos = [(_T_of_t(d, t), t) for t in ticks_s]
    pos = [(T, t) for T, t in pos if lo <= T <= hi]
    top.xaxis.set_major_locator(FixedLocator([T for T, _ in pos]))
    top.xaxis.set_minor_locator(FixedLocator([]))
    top.xaxis.set_major_formatter(FixedFormatter([_tick_time(t) for _, t in pos]))
    top.set_xlabel(label, color=S.INK_SECONDARY, fontsize=8.5)
    top.tick_params(colors=S.MUTED, labelsize=7.5)
    top.grid(False)
    for side in ("left", "right", "bottom"):
        top.spines[side].set_visible(False)
    top.spines["top"].set_visible(True)
    top.spines["top"].set_color(S.AXIS)
    return top


def _tick_time(t):
    names = {1e-12: "1 ps", 1e-9: "1 ns", 1e-6: "1 µs", 1e-3: "1 ms", 1: "1 s", 60: "1 min",
             3600: "1 hour", 86400: "1 day", YEAR: "1 year", 1e3 * YEAR: "1,000 yr",
             1e6 * YEAR: "1 Myr", 1e9 * YEAR: "1 Gyr", 1e-15: "1 fs", 1e-14: "10 fs"}
    for k, v in names.items():
        if abs(np.log(t / k)) < 1e-6:
            return v
    return S.human_time(t)


# --- Figure 1: the staircase ----------------------------------------------------------------------


def fig_staircase(res, d) -> list:
    T = d["T"]
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(10.5, 7.8), sharex=True,
                                  gridspec_kw={"height_ratios": [3.2, 1.25], "hspace": 0.07})
    ssT = d["ss_T"]
    m = (ssT >= T[0]) & (ssT <= T[-1])
    ax.fill_between(ssT[m], (d["ss_g_rho"] - d["ss_g_rho_err"])[m], (d["ss_g_rho"] + d["ss_g_rho_err"])[m],
                    color=S.lighten(S.C_REF1, 0.78), lw=0, zorder=1)
    ax.plot(ssT[m], d["ss_g_rho"][m], color=S.C_REF1, lw=1.0, zorder=2)
    ax.plot(T, d["g_rho_ideal"], color=S.C_IDEAL, lw=1.2, zorder=3)
    ax.plot(T, d["g_s"], color=S.INK_SECONDARY, lw=0.9, zorder=4)
    ax.plot(T, d["g_rho"], color=S.C_OURS, lw=1.9, zorder=5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(T[-1], T[0])
    ax.set_ylim(2.6, 140)
    ax.set_ylabel("g*   (effective number of relativistic species, log scale)")
    ax.yaxis.set_major_locator(FixedLocator([3, 5, 10, 20, 50, 100]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=20))
    ax.xaxis.set_minor_formatter(NullFormatter())
    _time_axis(ax, d, ticks_s=[1e-14, 1e-12, 1e-9, 1e-6, 1e-3, 1, 60, 3600])
    ax.set_title("The staircase: how many kinds of particle the hot universe held", pad=30)

    # step labels: where each heavy species' share of the entropy has halved
    h = res["half_temperatures"]
    gT = lambda TT: float(np.interp(np.log(TT), np.log(T), d["g_rho_ideal"]))
    # (key, label, factor on g at the label: >1 above the curve, <1 below it)
    steps = [("top", "t", 0.90), ("higgs", "H", 0.83), ("Z", "Z", 0.765), ("W", "W", 0.70),
             ("bottom", "b", 1.16), ("tau", "τ", 1.16), ("charm", "c", 1.18), ("pions", "π", 1.45),
             ("muon", "μ", 1.32), ("electron", "e±", 1.6)]
    for key, lab, fac in steps:
        Th = h[key]["T_half_GeV"]
        g0 = gT(Th)
        above = fac > 1
        ax.annotate(lab, xy=(Th, g0 * (1.02 if above else 0.985)), xytext=(Th, g0 * fac), ha="center",
                    va="bottom" if above else "top", fontsize=8.5, color=S.C_IDEAL, fontweight="bold",
                    arrowprops={"arrowstyle": "-", "color": S.lighten(S.C_IDEAL, 0.4), "lw": 0.6,
                                "shrinkA": 0, "shrinkB": 0})
    T_pc = C.LIT["T_pc_hotqcd2019"]["value"]
    ax.axvspan(T_pc + 0.0015, T_pc - 0.0015, color=S.WASH, zorder=0)
    S.annotate(ax, "quarks and gluons bind into hadrons\n(lattice QCD crossover, 156.5 MeV)",
               (T_pc, 26.5), (2.4, 21.5), ha="center", va="center", fontsize=7.5)
    S.annotate(ax, "the ideal-gas picture swaps free quarks\nand gluons for hadrons at once",
               (0.162, 55), (0.012, 62), color=S.C_IDEAL, ha="center", va="center", fontsize=7.5)
    Td = res["neutrinos"]["decoupling_T_MeV"]
    S.annotate(ax, f"neutrinos decouple\n({Td['nu_mu']:.1f} and {Td['nu_e']:.1f} MeV)", (2.0e-3, 10.75),
               (2.2e-3, 6.0), ha="center", va="center", fontsize=7.5)
    # plateau values
    ck = res["checks"]
    tf = ck["today_full"]
    ax.text(9e3, 113, "ideal gases with PDG masses: 106.75,\nevery Standard Model particle",
            fontsize=7.6, color=S.C_IDEAL, ha="left", va="bottom")
    ax.text(9e3, 62, f"this work: ideal gases + hadron gas +\nlattice QCD + perturbative QCD\n"
            f"(g*ρ = {ck['g_at_10TeV']['lattice']:.1f} at 10 TeV)", fontsize=7.6, color=S.C_OURS,
            ha="left", va="center")
    ax.text(22, 46, "Saikawa & Shirai 2018\n(band: their uncertainty)", color=S.C_REF1, fontsize=7.6,
            ha="center", va="center")
    S.annotate(ax, "", (12, 80), (20, 52), color=S.C_REF1)
    ax.text(1.5e-3, 11.9, "10.75 = γ + e± + 3 ν", fontsize=7.5, color=S.INK_SECONDARY,
            ha="center", va="bottom")
    ax.text(1.1e-5, 4.55, f"today: g*s = {tf['g_s']:.3f}", fontsize=7.6, color=S.INK_SECONDARY,
            ha="right", va="bottom")
    ax.text(1.1e-5, 3.2, f"today: g*ρ = {tf['g_rho']:.3f}", fontsize=7.6, color=S.C_OURS, ha="right", va="top")

    # residuals against Saikawa & Shirai
    S_on = np.interp(np.log(T), np.log(ssT), d["ss_g_rho"])
    S_err = np.interp(np.log(T), np.log(ssT), d["ss_g_rho_err"])
    axr.fill_between(T, -100 * S_err / S_on, 100 * S_err / S_on, color=S.lighten(S.C_REF1, 0.78), lw=0)
    axr.axhline(0, color=S.C_REF1, lw=0.9)
    axr.fill_between(T, 100 * (d["band_lo"] / S_on - 1), 100 * (d["band_hi"] / S_on - 1),
                     color=S.lighten(S.C_OURS, 0.75), lw=0)
    axr.plot(T, 100 * (d["g_rho_ideal"] / S_on - 1), color=S.C_IDEAL, lw=1.1)
    axr.plot(T, 100 * (d["g_rho"] / S_on - 1), color=S.C_OURS, lw=1.6)
    axr.set_ylim(-14, 22)
    axr.set_ylabel("difference from\nSaikawa & Shirai (%)")
    axr.set_xlabel("Temperature of the universe (kT)   ←  hotter, earlier   ·   cooler, later  →")
    axr.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    mx = res["comparisons"]["ideal_max_ratio_to_SS"]
    axr.annotate(f"ideal gases: +{100 * (mx['ratio'] - 1):.0f}% at {1000 * mx['T_GeV']:.0f} MeV",
                 xy=(mx["T_GeV"] * 1.04, 21.0), xytext=(9.0, 15.0), color=S.C_IDEAL, fontsize=7.5,
                 arrowprops={"arrowstyle": "->", "color": S.C_IDEAL, "lw": 0.7})
    frac = res["comparisons"]["fraction_of_grid_within_SS_error"]
    axr.text(9e3, -12.0, f"this work: inside their band at {100 * frac:.0f}% of temperatures; "
             "grey: our own systematic range", fontsize=7.3, color=S.C_OURS, ha="left")
    S.caption(fig, "g*ρ = energy density / (π²T⁴/30), g*s = entropy density / (2π²T³/45), with T the photon "
              "temperature. Letters mark where each species' share of the entropy has halved (T ≈ m/2.5–2.7). "
              "Top axis: the age of the universe at that temperature, from this project's calendar.", y=0.035)
    return S.save(fig, C.FIG_DIR / "01_staircase")


# --- Figure 2: what the universe was made of -------------------------------------------------------

GROUPS = [  # bottom to top; colours in the validated slot order (adjacent pairs)
    ("photons", ("photon",), S.SLOTS[0]),
    ("neutrinos", ("nu_e", "nu_mu", "nu_tau"), S.SLOTS[1]),
    ("electrons & positrons", ("electron", "qed"), S.SLOTS[2]),
    ("muons & taus", ("muon", "tau"), S.SLOTS[3]),
    ("W, Z & Higgs bosons", ("W", "Z", "higgs"), S.SLOTS[4]),
    ("quarks & gluons, then hadrons", ("strong",), S.SLOTS[5]),
]


def _group_fractions(d):
    names = [str(x) for x in d["parts_names"]]
    parts = d["parts_g_rho"]
    tot = parts.sum(axis=0)
    out = []
    for label, keys, col in GROUPS:
        g = sum(parts[names.index(k)] for k in keys if k in names)
        out.append((label, g / tot, col))
    return out


def fig_composition(res, d) -> list:
    T = d["T"]
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    fr = _group_fractions(d)
    base = np.zeros_like(T)
    for label, f, col in fr:
        ax.fill_between(T, 100 * base, 100 * (base + f), color=S.lighten(col, 0.45), lw=0)
        ax.plot(T, 100 * (base + f), color=S.SURFACE, lw=1.2)  # 2-px surface gap between bands
        base = base + f
    ax.set_xscale("log")
    ax.set_xlim(T[-1], T[0])
    ax.set_ylim(0, 100)
    ax.set_ylabel("share of the energy density (%)")
    ax.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    ax.set_xlabel("Temperature of the universe (kT)   ←  hotter, earlier   ·   cooler, later  →")
    ax.grid(False)
    _time_axis(ax, d, ticks_s=[1e-14, 1e-12, 1e-9, 1e-6, 1e-3, 1, 60, 3600])
    ax.set_title("What the universe was made of: shares of the energy, 10 TeV to 10 keV", pad=30)
    # direct labels inside the bands, at a temperature where each band is wide
    spots = {"quarks & gluons, then hadrons": (3e2, None), "W, Z & Higgs bosons": (1.2e2, None),
             "muons & taus": (5.0e-2, None), "electrons & positrons": (1.5e-3, None),
             "neutrinos": (1.0e-4, None), "photons": (4e-5, None)}
    cum = np.zeros_like(T)
    for label, f, col in fr:
        Tl = spots[label][0]
        i = int(np.argmin(np.abs(np.log(T / Tl))))
        mid = 100 * (cum[i] + f[i] / 2)
        if f[i] > 0.035:
            ax.text(Tl, mid, label, ha="center", va="center", fontsize=8, color=S.INK)
        cum = cum + f
    # the strong sector's two faces
    sfrac = fr[-1][1]
    i = int(np.argmin(np.abs(np.log(T / 0.105))))
    ax.text(0.105, 100 * (1 - sfrac[i] / 2), "hadrons", ha="center", va="center", fontsize=7.5)
    T_pc = C.LIT["T_pc_hotqcd2019"]["value"]
    ax.axvline(T_pc, color=S.INK_SECONDARY, lw=0.8)
    ax.text(T_pc * 1.06, 2, "quark–hadron\ncrossover", fontsize=7.2, color=S.INK, ha="right", va="bottom")
    ax.axvline(159.5, color=S.INK_SECONDARY, lw=0.8)
    ax.text(159.5 * 1.06, 2, "electroweak\ncrossover", fontsize=7.2, color=S.INK, ha="right", va="bottom")
    hi = fr[-1][1][-1]
    nu = fr[1][1]
    i10 = int(np.argmin(np.abs(np.log(T / 1e-2))))
    S.caption(fig, f"From this project's equation of state (lattice build). At 10 TeV quarks and gluons carry "
              f"{100 * hi:.0f}% of the energy; between 1 GeV and 120 MeV the strongly interacting matter "
              f"is lattice QCD (neither free quarks nor separate hadrons). Neutrinos carry {100 * nu[i10]:.0f}% at "
              f"10 MeV and {100 * nu[0]:.0f}% after e± annihilation, whose energy goes to the photons.", y=-0.02)
    return S.save(fig, C.FIG_DIR / "02_composition")


# --- Figure 3: the quark-hadron crossover ------------------------------------------------------------


def fig_qcd(res, d) -> list:
    import particles
    import qcd
    import thermo

    hot, wb = qcd.hotqcd2014(), qcd.wb2014()
    Tf = np.linspace(0.100, 0.420, 321)
    Tg = np.logspace(np.log10(3e-3), np.log10(0.6), 1800)
    strong = qcd.strong_sector(Tg)
    hs = particles.hadrons()
    pions = [h for h in hs if h.pdgid in (111, 211, -211)]
    h_all, h_pi = qcd.hrg(Tf, hs), qcd.hrg(Tf, pions)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.9), gridspec_kw={"width_ratios": [1.55, 1], "wspace": 0.26})
    T_pc = C.LIT["T_pc_hotqcd2019"]["value"]
    ax.axvspan(1000 * (T_pc - 0.0015), 1000 * (T_pc + 0.0015), color=S.WASH, zorder=0)
    ax.axvline(1000 * T_pc, color=S.AXIS, lw=0.8, zorder=0)
    sb = 2 * np.pi ** 2 / 45 * 47.5
    ax.axhline(sb, color=S.C_IDEAL, lw=1.0)
    ax.text(418, sb - 0.3, "free u, d, s quarks + gluons (ideal gas): 20.8", color=S.C_IDEAL, fontsize=7.5,
            ha="right", va="top")
    ax.plot(1000 * Tf, h_all["s"], color=S.C_IDEAL, lw=1.4)
    ax.plot(1000 * Tf, h_pi["s"], color=S.lighten(S.C_IDEAL, 0.35), lw=1.0)
    m = (Tg > 0.1) & (Tg < 0.42)
    ax.plot(1000 * Tg[m], strong["s"][m], color=S.C_OURS, lw=1.8, zorder=5)
    ax.errorbar(1000 * hot["T"], hot["s"], yerr=[hot["s_lo"], hot["s_hi"]], fmt="o", ms=3.2, color=S.C_REF1,
                ecolor=S.lighten(S.C_REF1, 0.35), elinewidth=0.8, capsize=0, zorder=4,
                markeredgecolor=S.SURFACE, markeredgewidth=0.6)
    ax.errorbar(1000 * wb["T"], wb["s"], yerr=wb["s_err"], fmt="s", ms=3.0, color=S.C_REF2,
                ecolor=S.lighten(S.C_REF2, 0.35), elinewidth=0.8, capsize=0, zorder=3,
                markeredgecolor=S.SURFACE, markeredgewidth=0.6)
    ax.set_xlim(100, 420)
    ax.set_ylim(0, 22)
    ax.set_xlabel("Temperature (MeV)")
    ax.set_ylabel("entropy density of QCD matter, s / T³")
    ax.set_title("Quarks and gluons become hadrons: lattice QCD vs ideal gases")
    q = res["qcd"]
    T_lab = 0.2135
    s_lab = float(qcd.hrg(np.array([T_lab]), hs)["s"][0])
    S.annotate(ax, "hadron gas: every PDG hadron\n(this work; it runs away above ~170 MeV)",
               (1000 * T_lab, s_lab), (232, 18.9), color=S.C_IDEAL, ha="left", va="center", fontsize=7.6)
    ax.text(300, 0.55, "pions only (textbook)", color=S.C_IDEAL, fontsize=7.5, ha="center", va="center")
    S.annotate(ax, "HotQCD 2014 (lattice, 2+1 flavours)", (395, 17.0), (418, 13.3), color=S.C_REF1,
               ha="right", va="center", fontsize=7.6)
    S.annotate(ax, "Borsányi et al. 2014 (lattice, 2+1)", (405, 16.0), (418, 12.1), color=S.C_REF2,
               ha="right", va="center", fontsize=7.6)
    S.annotate(ax, "this work's strong sector:\nhadron gas → Borsányi et al. 2016 fit (2+1+1)", (300, 15.5),
               (418, 9.6), color=S.C_OURS, ha="right", va="center", fontsize=7.6)
    S.annotate(ax, f"chiral crossover\n{1000 * T_pc:.1f} ± 1.5 MeV\n(HotQCD 2019)", (1000 * T_pc, 6.0), (122, 8.0),
               ha="left", fontsize=7.4)
    st = q["steepest_entropy_rise_MeV"]
    ax.text(418, 3.4, f"s/T³ rises fastest at {st['B16_fit_our_integration']:.0f} MeV (2016 fit),\n"
            f"{st['HotQCD2014_table']:.0f} and {st['WB2014_table']:.0f} MeV (2014 tables)", fontsize=7.2,
            color=S.INK_SECONDARY, ha="right", va="center")

    # how many resonances the hadron gas needs (g*rho of QCD matter at 150 MeV)
    T150 = 0.150
    hs_sorted = sorted(hs, key=lambda h: h.mass)
    contrib = np.array([h.g * thermo.functions(np.array([h.mass / T150]), h.stat)["rho"][0] * 30 / np.pi ** 2
                        for h in hs_sorted])
    masses = np.array([h.mass for h in hs_sorted])
    states = np.cumsum([h.g for h in hs_sorted])
    cum = np.cumsum(contrib)
    ax2.step(np.concatenate([[0.1], masses]), np.concatenate([[0], cum]), where="post", color=S.C_IDEAL, lw=1.4)
    lat = q["lattice_g_rho_strong_150MeV"]
    ax2.axhline(lat, color=S.C_OURS, lw=1.6)
    ax2.text(2.95, lat + 0.25, f"lattice QCD: {lat:.1f}", color=S.C_OURS, fontsize=7.6, ha="right", va="bottom")
    ax2.set_xlim(0.1, 3.0)
    ax2.set_ylim(0, 15.5)
    ax2.set_xlabel("heaviest hadron included (GeV)")
    ax2.set_ylabel("g*ρ of QCD matter at 150 MeV")
    ax2.set_title("How many hadrons it takes")
    for mcut in (0.2, 1.0, 2.0):
        i = np.searchsorted(masses, mcut) - 1
        ax2.annotate(f"{states[i]:,.0f} states", xy=(mcut, cum[i]), xytext=(mcut + 0.12, cum[i] - 1.7),
                     fontsize=7.3, color=S.C_IDEAL, arrowprops={"arrowstyle": "-", "color": S.C_IDEAL, "lw": 0.6})
    ax2.text(2.95, 1.0, f"free quarks + gluons: {q['free_quark_gluon_g_rho_150MeV']:.0f} (off scale)",
             color=S.C_IDEAL, fontsize=7.3, ha="right")
    ax2.text(0.15, 15.3, "PDG hadrons added lightest first,\neach charge and spin state counted",
             fontsize=7.2, color=S.INK_SECONDARY, va="top")
    S.caption(fig, "Left: lattice tables of HotQCD (Bazavov et al. 2014, Table I) and Borsányi et al. (2014), with "
              "errors; the hadron gas computed here from the 2026 PDG table (454 hadrons, 1,680 states, f0(500) and "
              "K0*(700) left out). Right: the same hadron gas at 150 MeV, adding states in order of mass.", y=-0.03)
    return S.save(fig, C.FIG_DIR / "03_qcd_crossover")


# --- Figure 4: the cosmic calendar ---------------------------------------------------------------------


def fig_calendar(res, d) -> list:
    cal_T, cal_t = d["cal_T"], d["cal_t"]
    fig, ax = plt.subplots(figsize=(12.5, 8.6))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-15, 2e18)
    ax.set_ylim(1e-13, 2e4)
    # eras (from the computed events)
    ev = {e["key"]: e for e in res["calendar"]["events"]}
    le = res["calendar"]["late_epochs"]
    t_mL = le["matter_lambda_equality"]["t_Gyr"] * 1e9 * YEAR
    eras = [(1e-15, ev["ew_crossover"]["t_s"], "Higgs field\noff"),
            (ev["ew_crossover"]["t_s"], ev["qcd_crossover"]["t_s"], "quark–gluon\nplasma"),
            (ev["qcd_crossover"]["t_s"], ev["muon_annihilation"]["t_s"], "hadrons"),
            (ev["muon_annihilation"]["t_s"], ev["nuclei_form"]["t_s"], "leptons and\nneutrinos"),
            (ev["nuclei_form"]["t_s"], 20 * 60, "nuclei\nform"),
            (20 * 60, ev["equality"]["t_s"], "radiation era"),
            (ev["equality"]["t_s"], t_mL, "matter era"),
            (t_mL, 2e18, "dark\nenergy")]
    for k, (a, b, name) in enumerate(eras):
        if k % 2 == 0:
            ax.axvspan(a, b, color=S.WASH, zorder=0, lw=0)
        ax.text(np.sqrt(a * b), 1.0e4, name, ha="center", va="top", fontsize=7.6, color=S.INK_SECONDARY)
    ax.plot(cal_t, cal_T, color=S.C_OURS, lw=2.0, zorder=4)
    # events above the line, temperatures we know below it; label anchors in (log10 t, log10 T)
    lay = {
        "ew_crossover": ("Electroweak crossover: the Higgs field switches on", -9.3, 3.15, "left"),
        "top_annihilation": ("top quarks annihilate", -8.4, 2.45, "left"),
        "bottom_annihilation": ("bottom quarks gone", -5.6, 1.55, "left"),
        "charm_annihilation": ("charm quarks gone", -4.6, 0.85, "left"),
        "qcd_crossover": ("Quarks and gluons bind into protons, neutrons and pions", -3.4, 0.15, "left"),
        "pion_annihilation": ("pions gone", -2.1, -0.55, "left"),
        "muon_annihilation": ("muons gone", -1.3, -1.25, "left"),
        "nu_mu_decoupling": ("neutrinos stop interacting", 0.35, -1.95, "left"),
        "np_freeze_out": ("the neutron-to-proton ratio freezes (about 1/6)", 1.55, -2.65, "left"),
        "electron_annihilation": ("electrons and positrons annihilate", 2.75, -3.35, "left"),
        "nuclei_form": ("deuterium survives: nuclei form (done in ~3 min)", 3.75, -4.05, "left"),
        "equality": ("matter outweighs radiation", 11.3, -6.7, "left"),
        "last_scattering": ("atoms form, CMB released", 13.0, -7.65, "left"),
        "today": ("today, 2.7255 K", 16.95, -12.62, "right"),
    }
    for e in res["calendar"]["events"]:
        if e["key"] not in lay:
            continue
        lab, lx, ly, ha = lay[e["key"]]
        x, y = e["t_s"], e["T_GeV"]
        ax.plot([x], [y], "o", ms=5, color=S.C_OURS, markeredgecolor=S.SURFACE, markeredgewidth=1.0, zorder=5)
        txt = f"{lab}\n{S.human_time(x)} · {S.human_energy(y)}"
        ax.annotate(txt, xy=(x, y), xytext=(10 ** lx, 10 ** ly), fontsize=7.4, ha=ha, va="center", color=S.INK,
                    arrowprops={"arrowstyle": "-", "color": S.AXIS, "lw": 0.6, "shrinkA": 2, "shrinkB": 3})
    labs = {l["key"]: l for l in res["calendar"]["labs"]}
    lab_lay = {
        "lhc_proton": ("an LHC proton carries 6.8 TeV:\nthe average particle did at", -14.85, 0.75, "left"),
        "alice_qgp": ("as hot as the quark–gluon plasma at the LHC\n(ALICE photons, 297 ± 43 MeV)", -6.6, -1.55, "right"),
        "phenix_qgp": ("…and at RHIC (PHENIX photons, 221 ± 27 MeV)", -6.0, -2.45, "right"),
        "sun_core": ("as hot as the centre of the Sun\n(15.6 million K)", 4.5, -6.85, "right"),
        "sun_surface": ("as hot as the surface of the Sun (5,772 K)", 11.2, -9.95, "right"),
        "water_boils": ("water would boil everywhere (373 K)", 13.0, -10.85, "right"),
        "water_freezes": ("…and freeze (273 K)", 13.4, -11.75, "right"),
    }
    for key, (lab, lx, ly, ha) in lab_lay.items():
        l = labs[key]
        x, y = l["t_s"], l["T_GeV"]
        ax.plot([x], [y], "D", ms=4.6, color=S.C_REF1, markeredgecolor=S.SURFACE, markeredgewidth=0.8, zorder=6)
        txt = f"{lab}\n{S.human_time(x)}" if key != "lhc_proton" else f"{lab} {S.human_time(x)}"
        ax.annotate(txt, xy=(x, y), xytext=(10 ** lx, 10 ** ly), fontsize=7.2, ha=ha, va="center",
                    color=S.C_REF1, arrowprops={"arrowstyle": "-", "color": S.lighten(S.C_REF1, 0.4), "lw": 0.6,
                                                "shrinkA": 2, "shrinkB": 3})
    ax.set_xlabel("Age of the universe")
    ax.set_ylabel("Temperature (kT)")
    xt = [1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1, 60, 3600, 86400, YEAR, 1e3 * YEAR, 1e6 * YEAR, 1e9 * YEAR]
    ax.xaxis.set_major_locator(FixedLocator(xt))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.xaxis.set_major_formatter(FixedFormatter([_tick_time(x) for x in xt]))
    yt = [1e4, 1e2, 1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
    ax.yaxis.set_major_locator(FixedLocator(yt))
    ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: S.human_energy(v).replace("0.001 eV", "1 meV")
                                               if v >= 1e-9 else f"{v * 1e12:g} meV"))
    axk = ax.twinx()
    axk.set_yscale("log")
    axk.set_ylim(np.array(ax.get_ylim()) / C.K_B_GEV_PER_K)
    axk.set_ylabel("Temperature (kelvin)")
    axk.grid(False)
    axk.spines["right"].set_visible(True)
    axk.spines["right"].set_color(S.AXIS)
    axk.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
    axk.yaxis.set_minor_formatter(NullFormatter())
    axk.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"10$^{{{int(round(np.log10(v)))}}}$ K"
                                                if v > 0 else ""))
    ax.set_title("A calendar of the hot universe: temperature against age, from 10 TeV to today", pad=10)
    S.caption(fig, "The line is this project's integration of the Friedmann equation with its own g*(T) (Planck 2018 "
              "matter and dark energy at late times); it reproduces Escudero's (2020) clock to 0.01% and gives an "
              "age of 13.787 Gyr today. Dots: events computed here (half-points for annihilations; Γ = H for "
              "neutrinos) or taken from lattice QCD, the PDG and Planck. Diamonds: temperatures known on Earth "
              "and in the Sun. Ages count from the formal t = 0 of a hot, radiation-filled universe.", y=-0.01)
    return S.save(fig, C.FIG_DIR / "04_cosmic_calendar")


# --- Figure 5: neutrino decoupling ------------------------------------------------------------------------


def fig_neutrinos(res, d) -> list:
    import neutrinos as NU

    fig = plt.figure(figsize=(11.5, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.42, wspace=0.24)
    ax = fig.add_subplot(gs[0, 0])
    axe = fig.add_subplot(gs[1, 0], sharex=ax)
    axg = fig.add_subplot(gs[0, 1])
    axn = fig.add_subplot(gs[1, 1])
    Tg, Te, Tm = d["nu_T"], d["nu_Te"], d["nu_Tm"]
    Ti, Tni = d["nu_inst_T"], d["nu_inst_Tn"]
    inst = np.interp(np.log(Tg[::-1]), np.log(Ti[::-1]), (Tni / Ti)[::-1])[::-1]
    ax.plot(Tg, inst, color=S.C_IDEAL, lw=1.2)
    ax.plot(Tg, Te / Tg, color=S.C_OURS, lw=1.8)
    ax.plot(Tg, Tm / Tg, color=S.INK_SECONDARY, lw=1.0)
    esc_T, esc_r = d["esc_T"], d["esc_ratio"]
    k = slice(None, None, 40)
    ax.plot(esc_T[k], 1 / esc_r[k], "o", ms=3.4, color=S.C_REF1, markeredgecolor=S.SURFACE, markeredgewidth=0.6)
    ax.set_xscale("log")
    ax.set_xlim(1e-2, 5e-6)
    ax.set_ylim(0.69, 1.01)
    ax.set_ylabel("neutrino / photon temperature")
    ax.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    ax.set_title("e± annihilation heats the photons, not the neutrinos")
    ax.axhline((4 / 11) ** (1 / 3), color=S.AXIS, lw=0.7, zorder=0)
    ax.text(9e-3, (4 / 11) ** (1 / 3) - 0.004, "(4/11)^(1/3) = 0.7138", fontsize=7.4, color=S.INK_SECONDARY, va="top")
    S.annotate(ax, "textbook: neutrinos gone\nbefore annihilation began", (2.4e-5, 0.7145), (3.5e-5, 0.775),
               color=S.C_IDEAL, ha="center", va="center", fontsize=7.5)
    S.annotate(ax, r"this work: $\nu_e$ (heavy) and $\nu_{\mu,\tau}$ (thin)", (5e-5, 0.7175), (5.5e-5, 0.85),
               color=S.C_OURS, ha="center", va="center", fontsize=7.5)
    ax.text(7e-3, 0.93, "Escudero 2020 (dots,\none common ν temperature)", color=S.C_REF1, fontsize=7.4)
    # the extra heating, as a percentage over the textbook ratio
    for col, lw, arr in ((S.C_OURS, 1.8, Te), (S.INK_SECONDARY, 1.0, Tm)):
        axe.plot(Tg, 100 * (arr / Tg / inst - 1), color=col, lw=lw)
    cn = res["neutrinos"]["cases"]["fdme_qed3_separate"]
    axe.set_ylabel("neutrino temperature above\nthe textbook value (%)")
    axe.set_xlabel("photon temperature (kT)")
    axe.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    fe = 100 * (1.40102 / cn["T_gamma_over_T_nu_e"] - 1)
    fm = 100 * (1.40102 / cn["T_gamma_over_T_nu_mu"] - 1)
    axe.text(4e-5, fe + 0.03, rf"$\nu_e$: +{fe:.2f}%", fontsize=7.5, color=S.C_OURS, va="bottom")
    axe.text(4e-5, fm - 0.04, rf"$\nu_{{\mu,\tau}}$: +{fm:.2f}%", fontsize=7.5, color=S.INK_SECONDARY, va="top")
    axe.set_ylim(-0.05, 0.66)
    axe.text(9e-3, 0.64, "late annihilations still heat the neutrinos a little, and the\n"
             "photon–electron plasma is not quite an ideal gas (QED)", fontsize=7.3, color=S.INK_SECONDARY, va="top")

    # rates
    Tr = np.logspace(-4, -1.5, 300)
    rr = NU.relaxation_rates(Tr)
    H = NU.hubble_equilibrium(Tr)
    axg.plot(Tr, rr["nu_e"] / H, color=S.C_OURS, lw=1.8)
    axg.plot(Tr, rr["nu_mu"] / H, color=S.INK_SECONDARY, lw=1.0)
    axg.axhline(1, color=S.AXIS, lw=0.8)
    axg.set_xscale("log")
    axg.set_yscale("log")
    axg.set_xlim(3e-2, 1e-4)
    axg.set_ylim(1e-3, 1e4)
    axg.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    axg.set_ylabel("energy-exchange rate / expansion rate")
    axg.set_title("When neutrinos stop exchanging energy")
    td = res["neutrinos"]["decoupling_T_MeV"]
    for key, col, lab, xy_text in (("nu_e", S.C_OURS, r"$\nu_e$", (4.2e-4, 8.0)),
                                   ("nu_mu", S.INK_SECONDARY, r"$\nu_\mu$, $\nu_\tau$", (9e-3, 0.02))):
        x = td[key] / 1000
        axg.plot([x], [1], "o", ms=5, color=col, markeredgecolor=S.SURFACE)
        axg.annotate(f"{lab}: {td[key]:.2f} MeV", xy=(x, 1), xytext=xy_text, fontsize=7.6, color=col, ha="center",
                     arrowprops={"arrowstyle": "-", "color": col, "lw": 0.6})
    axg.text(2.5e-2, 1.4, "Γ = H", fontsize=7.4, color=S.INK_SECONDARY)
    axg.text(2.8e-2, 1.5e-3, r"rate $\propto G_F^2 T^5$, expansion $\propto T^2$:" + "\nthe ratio falls as T³",
             fontsize=7.3, color=S.INK_SECONDARY, va="bottom")

    # N_eff ladder
    cs = res["neutrinos"]["cases"]
    rows = [("textbook: instantaneous, ideal gas", 3.0, S.C_IDEAL),
            ("+ QED plasma corrections", cs["instantaneous_10MeV_qed"]["N_eff"], S.C_OURS),
            ("+ ongoing collisions (massless e)", cs["fd_qed3_separate"]["N_eff"], S.C_OURS),
            ("+ electron mass in the rates", cs["fdme_qed3_separate"]["N_eff"], S.C_OURS),
            ("Bennett et al. 2021, full calculation: 3.0440 ± 0.0002", 3.0440, S.C_REF1)]
    for i, (lab, v, col) in enumerate(rows):
        y = len(rows) - 1 - i
        axn.plot([3.0, v], [y, y], color=S.lighten(col, 0.6), lw=2.2, solid_capstyle="butt")
        axn.plot([v], [y], "o", ms=6, color=col, markeredgecolor=S.SURFACE)
        axn.text(3.0005, y + 0.28, lab, fontsize=7.4, color=S.INK, va="bottom")
        axn.text(v + 0.0007, y, f"{v:.4f}", fontsize=7.4, color=col, va="center")
    axn.set_xlim(2.998, 3.058)
    axn.set_ylim(-0.6, len(rows) - 0.1)
    axn.set_yticks([])
    axn.spines["left"].set_visible(False)
    axn.set_xlabel(r"$N_{\rm eff}$ (effective number of neutrino species)")
    axn.set_title(r"Why $N_{\rm eff}$ is 3.044, not 3")
    S.caption(fig, "This project solves Escudero's (2020) temperature equations with the QED corrections of Bennett "
              "et al. (2020) and the electron-mass rate coefficients of Escudero, Jackson, Laine & Sandner (2025). "
              f"It reproduces their published runs to every printed digit: N_eff = {cs['fdme_qed3_separate']['N_eff']:.4f}, "
              f"g*ρ = {cs['fdme_qed3_separate']['g_rho_today']:.4f} and g*s = {cs['fdme_qed3_separate']['g_s_today']:.4f} "
              "today. Rates: our linearised energy-exchange rate per degree of temperature difference.", y=0.015)
    return S.save(fig, C.FIG_DIR / "05_neutrino_decoupling")


# --- Figure 6: a softer universe and the horizon mass -----------------------------------------------------


def fig_horizon(res, d) -> list:
    T = d["T"]
    fig, (ax, axm) = plt.subplots(2, 1, figsize=(10.5, 8.2), sharex=True,
                                  gridspec_kw={"height_ratios": [1, 1.25], "hspace": 0.1})
    ssT = d["ss_w_T"]
    ax.fill_between(ssT, d["ss_w"] - d["ss_w_err"], d["ss_w"] + d["ss_w_err"], color=S.lighten(S.C_REF1, 0.78), lw=0)
    ax.plot(ssT, d["ss_w"], color=S.C_REF1, lw=1.0)
    ax.plot(d["lms_T"], d["lms_w"], color=S.C_REF2, lw=1.0)
    ax.plot(T, d["cs2"], color=S.INK_SECONDARY, lw=0.9)
    ax.plot(T, d["w"], color=S.C_OURS, lw=1.8)
    ax.axhline(1 / 3, color=S.AXIS, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(T[-1], T[0])
    ax.set_ylim(0.19, 0.345)
    ax.set_ylabel("pressure / energy density")
    ax.set_title("A softer universe: the pressure dips when particles annihilate or quarks bind")
    sof = res["eos_softening"]
    q = sof["qcd"]
    ax.annotate(f"QCD crossover: w = {q['w_min_ours']:.3f} at {1000 * q['T_w_min_ours_GeV']:.0f} MeV\n"
                f"(Saikawa & Shirai {q['w_min_SS']:.3f} at {1000 * q['T_w_min_SS_GeV']:.0f} MeV)",
                xy=(q["T_w_min_ours_GeV"], q["w_min_ours"]), xytext=(4e-3, 0.205), fontsize=7.5,
                arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.6}, color=S.INK)
    ax.text(5e3, 1 / 3 + 0.002, "w = 1/3: pure radiation", fontsize=7.4, color=S.INK_SECONDARY, va="bottom")
    ax.text(9e3, 0.226, "this work: w = p/ρ (heavy line)\nand sound speed² = dp/dρ (thin line)", fontsize=7.5,
            color=S.C_OURS, ha="left", va="center")
    ax.text(2.5, 0.272, "Saikawa & Shirai 2018 (w, band)", fontsize=7.5, color=S.C_REF1, ha="center")
    ax.text(0.03, 0.25, "Laine, Meyer &\nSchröder 2015", fontsize=7.5, color=S.C_REF2, ha="center")
    e = sof["electron"]
    ax.annotate(f"e± annihilation: w = {e['w_min_ours']:.3f}", xy=(e["T_w_min_ours_GeV"], e["w_min_ours"]),
                xytext=(2e-4, 0.24), fontsize=7.5, arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.6})
    w = sof["ew"]
    ax.annotate(f"W, Z, H, t annihilate: w = {w['w_min_ours']:.3f}", xy=(w["T_w_min_ours_GeV"], w["w_min_ours"]),
                xytext=(2e3, 0.296), fontsize=7.5, arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.6})

    cT, MH = d["cal_T"], d["cal_MH"]
    m = (cT >= T[0]) & (cT <= T[-1])
    axm.plot(cT[m], MH[m], color=S.C_OURS, lw=1.9)
    axm.set_yscale("log")
    axm.set_ylim(1e-11, 1e11)
    axm.set_ylabel("mass inside the Hubble sphere (solar masses)")
    axm.set_xlabel("Temperature of the universe (kT)   ←  hotter, earlier   ·   cooler, later  →")
    axm.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    refs = [(3.003e-6, "Earth"), (9.546e-4, "Jupiter"), (1.0, "Sun"), (4.3e6, "Sagittarius A*")]
    for mass, lab in refs:
        axm.axhline(mass, color=S.AXIS, lw=0.8, zorder=1)
        axm.text(9e3, mass * 1.35, lab, fontsize=7.3, color=S.INK_SECONDARY, va="bottom")
    pta = res["calendar"]["pta_band_T_GeV"]
    axm.axvspan(pta[1], pta[0], color=S.WASH, zorder=0.5)
    axm.text(np.sqrt(pta[0] * pta[1]), 3e9, "today's pulsar-timing band\n(2–100 nHz waves entered\nthe horizon here)",
             fontsize=7.3, color=S.INK_SECONDARY, ha="center", va="center")
    Mw = q["M_hubble_at_w_min_Msun"]
    axm.plot([q["T_w_min_ours_GeV"]], [Mw], "o", ms=5.5, color=S.C_OURS, markeredgecolor=S.SURFACE, zorder=5)
    axm.annotate(f"QCD softening: {Mw:.1f} M☉ inside the horizon\n"
                 "(primordial black holes forming then would peak\nat ≈0.7 M☉: Byrnes et al. 2018)",
                 xy=(q["T_w_min_ours_GeV"], Mw), xytext=(30, 2e-6), fontsize=7.4, color=S.INK,
                 arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.6})
    for key, lab, dx in (("ew_crossover", "electroweak crossover", 1.0), ("nu_e_decoupling", "neutrinos decouple", 1.0),
                         ("nuclei_form", "nuclei form", 1.0)):
        e2 = next(x for x in res["calendar"]["events"] if x["key"] == key)
        axm.plot([e2["T_GeV"]], [e2["M_hubble_Msun"]], "o", ms=4.5, color=S.C_OURS, markeredgecolor=S.SURFACE)
        axm.annotate(f"{lab}: {S.sci(e2['M_hubble_Msun'])} M☉", xy=(e2["T_GeV"], e2["M_hubble_Msun"]),
                     xytext=(10, -14), textcoords="offset points", fontsize=7.4, color=S.INK, va="top",
                     arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.6})
    _time_axis(ax, d, ticks_s=[1e-14, 1e-12, 1e-9, 1e-6, 1e-3, 1, 60, 3600])
    S.caption(fig, "Top: w = p/ρ and the sound speed squared from this project's equation of state, against Saikawa "
              "& Shirai's published w(T) and Laine, Meyer & Schröder's tables. Bottom: the mass-energy inside the "
              "Hubble sphere, c³/(2GH), from our calendar; it agrees with Byrnes et al.'s (2018) formula, "
              "1.5 × 10⁵ M☉ × (g/10.75)^(−1/2) × (T/1 MeV)^(−2), to 0.1%.", y=0.03)
    return S.save(fig, C.FIG_DIR / "06_softening_and_horizon_mass")


# --- Figure 7: checks against the literature ---------------------------------------------------------------


def fig_checks(res, d) -> list:
    T = d["T"]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10.5, 7.6), gridspec_kw={"hspace": 0.32, "height_ratios": [1.25, 1]})
    S_on = np.interp(np.log(T), np.log(d["ss_T"]), d["ss_g_rho"])
    S_err = np.interp(np.log(T), np.log(d["ss_T"]), d["ss_g_rho_err"])
    ax.fill_between(T, -100 * S_err / S_on, 100 * S_err / S_on, color=S.lighten(S.C_REF1, 0.78), lw=0)
    ax.axhline(0, color=S.C_REF1, lw=0.9)
    ax.fill_between(T, 100 * (d["band_lo"] / S_on - 1), 100 * (d["band_hi"] / S_on - 1),
                    color=S.lighten(S.C_OURS, 0.75), lw=0)
    L_on = np.interp(np.log(T), np.log(d["lms_T"]), d["lms_g_rho"])
    ax.plot(T, 100 * (L_on / S_on - 1), color=S.C_REF2, lw=1.0)
    Sb = np.interp(np.log(d["b16_T"]), np.log(d["ss_T"]), d["ss_g_rho"])
    ax.plot(d["b16_T"], 100 * (d["b16_g_rho"] / Sb - 1), "^", ms=5.5, color=S.C_OTHER,
            markeredgecolor=S.SURFACE, markeredgewidth=0.6)
    ax.plot(T, 100 * (d["g_rho"] / S_on - 1), color=S.C_OURS, lw=1.7)
    ax.set_xscale("log")
    ax.set_xlim(T[-1], T[0])
    ax.set_ylim(-10, 10)
    ax.set_ylabel("g*ρ: difference from\nSaikawa & Shirai 2018 (%)")
    ax.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    ax.set_title("Checks: this project's g* against three published Standard Model curves")
    ax.text(9e3, -3.6, "this work (grey band: our variants: switch temperatures, perturbative scale, hadron list)",
            fontsize=7.5, color=S.C_OURS, va="top")
    ax.text(9e3, -5.2, "Saikawa & Shirai 2018 (zero line; band = their error)", fontsize=7.5, color=S.C_REF1, va="top")
    ax.text(9e3, -6.8, "Laine, Meyer & Schröder 2015 (line; off scale near 150 MeV, older lattice input)",
            fontsize=7.5, color=S.C_REF2, va="top")
    ax.text(9e3, -8.4, "Borsányi et al. 2016 spline knots (triangles; they quote ~1% accuracy)", fontsize=7.5,
            color=S.INK_SECONDARY, va="top")
    # zoom: neutrino decoupling and e+- annihilation
    m = T <= 1.2e-2
    for key, col, lw in (("g_rho", S.C_OURS, 1.7), ("g_s", S.INK_SECONDARY, 1.0)):
        So = np.interp(np.log(T[m]), np.log(d["ss_T"]), d["ss_" + key])
        ax2.plot(T[m], 100 * (d[key][m] / So - 1), color=col, lw=lw)
        Se = np.interp(np.log(d["esc_T"]), np.log(d["ss_T"]), d["ss_" + key])
        k = slice(None, None, 30)
        ax2.plot(d["esc_T"][k], 100 * (d["esc_" + key][k] / Se[k] - 1), "s" if key == "g_rho" else "D",
                 ms=3.6, color=S.C_OTHER, markeredgecolor=S.SURFACE, markeredgewidth=0.5)
        Lo = np.interp(np.log(T[m]), np.log(d["lms_T"]), d["lms_" + key])
        ax2.plot(T[m], 100 * (Lo / So - 1), color=S.C_REF2, lw=0.9 if key == "g_rho" else 0.6)
    ax2.axhline(0, color=S.C_REF1, lw=0.9)
    ax2.set_xscale("log")
    ax2.set_xlim(1.2e-2, T[0])
    ax2.set_ylim(-0.75, 0.75)
    ax2.set_ylabel("difference from\nSaikawa & Shirai 2018 (%)")
    ax2.set_xlabel("Temperature of the universe (kT)")
    ax2.xaxis.set_major_formatter(FuncFormatter(_energy_formatter))
    ax2.set_title("Zoom: neutrino decoupling and e± annihilation", fontsize=9.5)
    ax2.text(1.1e-2, -0.22, "this work: g*ρ (heavy), g*s (thin)", fontsize=7.5, color=S.C_OURS, va="top")
    ax2.text(1.1e-2, -0.38, "Escudero 2020 (squares g*ρ, diamonds g*s; one common ν temperature)", fontsize=7.5,
             color=S.INK_SECONDARY, va="top")
    ax2.text(1.1e-2, -0.54, "Laine, Meyer & Schröder 2015 (sudden decoupling at 2 MeV)", fontsize=7.5,
             color=S.C_REF2, va="top")
    S.caption(fig, "All curves are divided by Saikawa & Shirai's. From 1 GeV down to 120 MeV we and they use the same "
              "lattice input (Borsányi et al. 2016) and agree to within 1%; above 1 GeV the spread is the poorly "
              "convergent perturbative QCD. Below 3 MeV our decoupling run follows Escudero's to 10⁻⁴; Saikawa & "
              "Shirai's curve runs up to 0.5% lower there, then reaches the same values today.", y=0.025)
    return S.save(fig, C.FIG_DIR / "07_checks")


def make_all(res=None) -> list:
    S.use_style()
    C.FIG_DIR.mkdir(parents=True, exist_ok=True)
    res2, d = _load()
    res = res or res2
    out = []
    for f in (fig_staircase, fig_composition, fig_qcd, fig_calendar, fig_neutrinos, fig_horizon, fig_checks):
        out += f(res, d)
    return out


if __name__ == "__main__":
    for p in make_all():
        print("wrote", p)
