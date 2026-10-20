"""Every figure, from the cached results (``data/results.json`` and
``data/derived/``). Run ``python analysis.py`` first.

    python figures.py            # all figures
    python figures.py 01 05      # only these
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

import config as cfg  # noqa: E402
import plotstyle as ps  # noqa: E402
import relic  # noqa: E402

FIG = cfg.FIGURE_DIR


def _results() -> dict:
    return json.loads(cfg.RESULTS_PATH.read_text())


# ---------------------------------------------------------------------------
# 01 The relic neutrinos today
# ---------------------------------------------------------------------------


def fig_relic(r: dict) -> list:
    rl = r["relic"]
    T_g, T_n = cfg.T_CMB_K, rl["T_nu_K"]
    kTg, kTn = relic.kT_ev(T_g), relic.kT_ev(T_n)
    hc3 = relic.HBARC_EV_CM**3
    p = np.linspace(1e-6, 3.0e-3, 2000)  # eV
    dn_g = (1 / np.pi**2) * p**2 / np.expm1(p / kTg) / hc3 * 1e-3  # per cm^3 per meV
    dn_n = (3 / np.pi**2) * p**2 / (np.exp(p / kTn) + 1) / hc3 * 1e-3

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.2), gridspec_kw={"wspace": 0.30})
    ax = axes[0]
    ax.fill_between(p * 1e3, dn_g, color=ps.lighten(ps.C_GAMMA, 0.8), lw=0)
    ax.fill_between(p * 1e3, dn_n, color=ps.lighten(ps.C_NU, 0.72), lw=0, alpha=0.85)
    ax.plot(p * 1e3, dn_g, color=ps.C_GAMMA, label=f"CMB photons, {T_g:.4f} K")
    ax.plot(p * 1e3, dn_n, color=ps.C_NU, label=f"relic neutrinos, {T_n:.3f} K")
    ax.text(0.47, 95, f"{rl['n_nu_total_cm3']:.0f} neutrinos\nper cm³\n(ν + ν̄, 3 flavours)",
            color=ps.INK, fontsize=7.6, ha="center")
    ax.text(1.12, 185, f"{rl['n_gamma_cm3']:.0f} photons\nper cm³", color=ps.INK, fontsize=7.6)
    ax.set_xlabel("momentum $pc$ today (meV)")
    ax.set_ylabel("number per cm³ per meV of momentum")
    ax.set_xlim(0, 2.2)
    ax.set_ylim(0, 640)
    ax.legend(loc="upper right")
    ax.set_title("a  Two relic gases share the sky")

    ax = axes[1]
    states = [(k, v) for k, v in rl["states"].items() if v["m_eV"] > 0]
    shades = [ps.SEQ_BLUE[2], ps.SEQ_BLUE[4], ps.SEQ_BLUE[6]]
    heights = [0.06, 0.20, 0.06]  # staggered so neighbouring labels do not touch
    for (label, s_), col, dy in zip(states, shades, heights):
        v, pdf = relic.speed_distribution(s_["m_eV"], T_n, n=6000)
        m = v > 30
        logpdf = pdf[m] * v[m] * np.log(10)
        ax.plot(v[m], logpdf, color=col)
        med = s_["quantiles_kms"]["0.50"]
        ax.plot([med], [np.interp(med, v[m], logpdf)], "o", color=col, ms=4,
                mec=ps.SURFACE, mew=0.8)
        vpk = v[m][np.argmax(logpdf)]
        mev = s_["m_eV"] * 1e3
        ax.text(vpk, logpdf.max() + dy, f"{mev:.1f} meV" if mev < 10 else f"{mev:.0f} meV",
                ha="center", va="bottom", fontsize=7.4, color=ps.INK)
    for x, txt in ((cfg.V_ESC_MW_KMS, "escape speed from\nthe Milky Way here"),
                   (relic.C_KM_S, "light")):
        ax.axvline(x, color=ps.MUTED, lw=0.8, ls=(0, (3, 3)))
        ax.text(x * 0.93, 0.97, txt, rotation=90, ha="right", va="top", fontsize=7,
                color=ps.INK_SECONDARY, transform=ax.get_xaxis_transform())
    ax.set_xscale("log")
    ax.set_xlim(150, 4e5)
    ax.set_ylim(0, 2.0)
    ax.set_xlabel("speed today (km/s)")
    ax.set_ylabel("fraction per decade of speed")
    ax.set_title("b  Massive relics are slow now")

    ax = axes[2]
    zp1 = np.geomspace(1, 2e4, 400)
    pmean = relic.mean_momentum_ev(T_n) * zp1
    ax.plot(zp1, pmean, color=ps.INK)
    ps.annotate(ax, "mean momentum $\\langle p\\rangle c$\nof a relic neutrino",
                xy=(4.0, relic.mean_momentum_ev(T_n) * 4.0), xytext=(25, 8e-4), color=ps.INK,
                fontsize=7.4)
    for (label, s_), col in zip(states, shades):
        ax.axhline(s_["m_eV"], color=col, lw=1.2)
        z_nr = s_["z_nr"]
        ax.plot([1 + z_nr], [s_["m_eV"]], "o", color=col, ms=5, mec=ps.SURFACE, mew=0.8)
        ax.text(1.6e4, s_["m_eV"] * 0.9, f"{s_['m_eV'] * 1e3:.0f} meV: since $z$ ≈ {z_nr:.0f}",
                fontsize=7.2, color=ps.INK_SECONDARY, ha="right", va="top")
    ax.axvline(1091, color=ps.MUTED, lw=0.8, ls=(0, (3, 3)))
    ax.text(1091 * 0.88, 2.6e-4, "CMB\nreleased", fontsize=7, color=ps.INK_SECONDARY, ha="right")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 2e4)
    ax.set_ylim(2e-4, 0.5)
    ax.set_xlabel("1 + redshift")
    ax.set_ylabel("energy (eV)")
    ax.set_title("c  When each became non-relativistic")
    ps.caption(fig, "Calculation from the entropy argument, $T_\\nu$ = (4/11)$^{1/3}$ × 2.7255 K. "
               "Masses: the two heavier states if the lightest is massless (normal ordering, "
               "NuFIT 6.0), and 0.1 eV for comparison. Dots in b mark median speeds; in c, the "
               "moment ⟨p⟩c fell below mc², i.e. the state became non-relativistic.", y=-0.04)
    return ps.save(fig, FIG / "01_relic_neutrinos")


# ---------------------------------------------------------------------------
# 02 Decoupling: why 3.044
# ---------------------------------------------------------------------------


def fig_decoupling(r: dict) -> list:
    runs = r["decoupling"]["runs"]
    off = runs[cfg.DECOUPLING_RUNS[0][0]]["track"]
    full = runs[cfg.DECOUPLING_HEADLINE]["track"]
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.3),
                             gridspec_kw={"wspace": 0.36, "width_ratios": [1, 1, 1.15]})

    def comoving(tr, key, t_ref=10.0):
        """a·T normalised to 1 at T_γ = t_ref (both runs share that reference)."""
        a = np.exp(np.array(tr["lna"]))
        aT = a * np.array(tr[key])
        T = np.array(tr["T"])
        ref = np.interp(np.log(t_ref), np.log(T[::-1]), (a * T)[::-1])
        return aT / ref

    ax = axes[0]
    T_off = np.array(off["T"])
    T = np.array(full["T"])
    m_off = T_off <= 10.0
    ax.plot(T_off[m_off], comoving(off, "T")[m_off], color=ps.C_GAMMA, lw=1.1, ls=(0, (4, 3)),
            label="photons, instantaneous decoupling")
    ax.plot(T, comoving(full, "T"), color=ps.C_GAMMA, label="photons, our calculation")
    ax.plot(T_off[m_off], comoving(off, "Tnue")[m_off], color=ps.C_NU, lw=1.1, ls=(0, (4, 3)),
            label="neutrinos, instantaneous")
    ax.plot(T, comoving(full, "Tnue"), color=ps.C_NU, label="$\\nu_e$, our calculation")
    ax.axhline((11 / 4) ** (1 / 3), color=ps.MUTED, lw=0.7, ls=":")
    ax.text(9, (11 / 4) ** (1 / 3) + 0.008, "(11/4)$^{1/3}$ = 1.401", fontsize=7,
            color=ps.INK_SECONDARY, ha="left", va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(10, 0.009)
    ax.set_ylim(0.97, 1.45)
    ax.set_xlabel("photon temperature (MeV)  →  later")
    ax.set_ylabel("comoving temperature $a\\,T$ (= 1 at 10 MeV)")
    ax.legend(loc="center left", fontsize=7.2)
    ax.set_title("a  e± annihilation heats the photons")

    ax = axes[1]
    for key, col, lab in (("Tnue", ps.C_NU, "$\\nu_e$ (also W exchange with electrons)"),
                          ("Tnumu", ps.lighten(ps.C_NU, 0.45), "$\\nu_\\mu$, $\\nu_\\tau$ (Z only)")):
        ax.plot(T, (comoving(full, key) - 1) * 1e3, color=col, label=lab)
    ge = np.array(full["gamma_e_over_H"])
    gm = np.array(full["gamma_mu_over_H"])
    for g, col in ((ge, ps.C_NU), (gm, ps.lighten(ps.C_NU, 0.45))):
        i = int(np.argmax(g < 1.0))
        if i > 0:
            ax.axvline(T[i], color=col, lw=0.8, ls=":")
    hd = runs[cfg.DECOUPLING_HEADLINE]
    ax.text(0.97, 0.10, f"dotted: interaction rate = expansion rate\nat $T$ = {hd['T_dec_mu_MeV']:.1f} MeV "
            f"($\\nu_\\mu$) and {hd['T_dec_e_MeV']:.1f} MeV ($\\nu_e$)", transform=ax.transAxes,
            fontsize=7, color=ps.INK_SECONDARY, ha="right")
    ax.set_xscale("log")
    ax.set_xlim(10, 0.009)
    ax.set_xlabel("photon temperature (MeV)  →  later")
    ax.set_ylabel("neutrino heating, ($a\\,T_\\nu$ − 1) × 1000")
    ax.legend(loc="upper left", fontsize=7.2)
    ax.set_title("b  …and the neutrinos catch a little")

    ax = axes[2]
    rows = []
    for label, key in cfg.DECOUPLING_WATERFALL:
        rows.append((f"{label}", runs[key]["neff"], runs[key].get("escudero_2020"), "ours"))
    lit = list(cfg.NEFF_SM_LITERATURE.items())
    y_ours = np.arange(len(rows))[::-1] + len(lit) + 1
    ax.plot([v for _, v, _, _ in rows], y_ours, color=ps.lighten(ps.C_NU, 0.5), lw=1, zorder=1)
    for yy, (lab, v, ref, _) in zip(y_ours, rows):
        ax.plot([v], [yy], "o", color=ps.C_NU, ms=6, mec=ps.SURFACE, mew=0.8, zorder=3)
        if ref is not None:
            ax.plot([ref], [yy], "o", mfc="none", mec=ps.INK, ms=9, mew=0.8, zorder=2)
        ax.text(v + 0.0015, yy + 0.18, f"{v:.4f}", fontsize=7, color=ps.INK, va="bottom")
    y_lit = np.arange(len(lit))[::-1]
    for yy, (lab, (v, e)) in zip(y_lit, lit):
        ax.errorbar([v], [yy], xerr=[[e], [e]], fmt="D", color=ps.INK, ms=4.5, capsize=0)
    ax.axhline(len(lit) - 0.0 + 0.4, color=ps.AXIS, lw=0.8)
    ax.set_yticks(list(y_ours) + list(y_lit))
    ax.set_yticklabels([row[0] for row in rows] + [k for k, _ in lit], fontsize=7)
    ax.axvline(3.0, color=ps.MUTED, lw=0.7, ls=":")
    ax.set_xlim(2.997, 3.058)
    ax.set_xlabel("$N_\\mathrm{eff}$")
    ax.grid(axis="y", visible=False)
    ax.set_title("c  Our $N_\\mathrm{eff}$, one ingredient at a time")
    ax.plot([], [], "o", color=ps.C_NU, label="ours")
    ax.plot([], [], "o", mfc="none", mec=ps.INK, label="Escudero 2020, Table 1")
    ax.plot([], [], "D", color=ps.INK, label="full calculations (literature)")
    ax.legend(loc="lower left", fontsize=6.8)
    ps.caption(fig, "Calculation: thermal temperatures for photons + e±, $\\nu_e$ and $\\nu_{\\mu,\\tau}$ "
               "evolved with the weak energy-transfer rates and finite-temperature QED (Escudero "
               "2020; electron-mass rate coefficients from Escudero et al. 2026).", y=-0.04)
    return ps.save(fig, FIG / "02_decoupling")


# ---------------------------------------------------------------------------
# 03 N_eff in the CMB
# ---------------------------------------------------------------------------


def fig_neff_spectra(r: dict) -> list:
    import fetch_data as fd

    s = dict(np.load(cfg.DERIVED_DIR / "neff_spectra.npz"))
    ell = s["ell"]
    neffs = s["neffs"]
    binned = fd.load_planck_binned_tt()
    cols = [ps.C_FEWER, ps.C_SM, ps.C_MORE]
    fig = plt.figure(figsize=(12.8, 7.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.38, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    m = ell >= 2
    for i, n in enumerate(neffs):
        ax.plot(ell[m], s[f"zeq_tt_{i}"][m], color=cols[i], lw=1.4 if i != 1 else 1.1,
                label=f"$N_\\mathrm{{eff}}$ = {n:g}" + (" (Standard Model)" if i == 1 else ""))
    ax.errorbar(binned[:, 0], binned[:, 1], yerr=binned[:, 2], fmt="o", ms=2.6, color=ps.INK,
                ecolor=ps.MUTED, elinewidth=0.8, capsize=0, label="Planck 2018 (binned)")
    ax.set_xlim(2, 2500)
    ax.set_ylim(0, 6500)
    ax.set_xlabel("multipole ℓ")
    ax.set_ylabel("$D_\\ell^{TT}$ (µK²)")
    ax.legend(loc="upper right")
    ax.set_title("More neutrino species, fainter small-scale ripples (acoustic angle and "
                 "matter–radiation equality held fixed)")
    titles = {"theta": "a  Only the acoustic angle held",
              "zeq": "b  + equality held: the damping tail",
              "thetad": "c  + damping angle held: the phase shift"}
    frac_err = binned[:, 2] / binned[:, 1]
    for j, tag in enumerate(("theta", "zeq", "thetad")):
        ax = fig.add_subplot(gs[1, j])
        key = "utt" if tag == "thetad" else "tt"
        ref = s[f"{tag}_{key}_1"]
        mm = ell >= 30
        ax.fill_between(binned[:, 0], -frac_err * 100, frac_err * 100, color=ps.WASH, lw=0,
                        label="Planck error per band of 30 ℓ")
        for i in (0, 2):
            ax.plot(ell[mm], (s[f"{tag}_{key}_{i}"][mm] / ref[mm] - 1) * 100, color=cols[i],
                    lw=1.1, label=f"$N_\\mathrm{{eff}}$ = {neffs[i]:g}")
        ax.axhline(0, color=ps.AXIS, lw=0.8)
        ax.set_xlim(30, 2500)
        lim = {"theta": 40, "zeq": 16, "thetad": 9}[tag]
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("multipole ℓ")
        if j == 0:
            ax.set_ylabel("change in $D_\\ell^{TT}$ (%)")
        ax.set_title(titles[tag], fontsize=9.5)
        if j == 0:
            ax.legend(loc="lower left", fontsize=7)
    ps.caption(fig, "Calculation with CAMB around Planck 2018's best fit. a: H0 adjusts to keep "
               "θ*; b: ω_c also rises with N_eff so z_eq stays put (ω_b, helium fixed); c: the "
               "helium fraction is also tuned to restore the damping angle (unlensed spectra).",
               y=0.03)
    return ps.save(fig, FIG / "03_neff_in_the_cmb")


# ---------------------------------------------------------------------------
# 07 Matter power suppression
# ---------------------------------------------------------------------------


def fig_matter_power(r: dict) -> list:
    m = dict(np.load(cfg.DERIVED_DIR / "matter_power.npz"))
    kh = m["kh"]
    masses = m["masses"]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.4), gridspec_kw={"wspace": 0.26})
    shades = ps.SEQ_BLUE[1:][: len(masses) - 1]
    ax = axes[0]
    ax.loglog(kh, m["fixed_h_lin_0"], color=ps.INK, label="massless neutrinos")
    for i, (mass, col) in enumerate(zip(masses[1:], shades), start=1):
        ax.loglog(kh, m[f"fixed_h_lin_{i}"], color=col, lw=1.2, label=f"$\\Sigma m_\\nu$ = {mass:g} eV")
    ax.set_xlabel("wavenumber k (h/Mpc)")
    ax.set_ylabel("$P(k)$ at $z$ = 0  ((Mpc/h)³)")
    ax.set_xlim(1e-4, 10)
    ax.legend(loc="lower left", fontsize=7.4)
    ax.set_title("a  The linear matter power spectrum")
    ax = axes[1]
    for i, (mass, col) in enumerate(zip(masses[1:], shades), start=1):
        fnu = float(m[f"fixed_h_fnu_{i}"])
        ax.semilogx(kh, (m[f"fixed_h_lin_{i}"] / m["fixed_h_lin_0"] - 1) * 100, color=col,
                    label=f"$\\Sigma m_\\nu$ = {mass:g} eV ($f_\\nu$ = {fnu * 100:.2f}%)")
        ax.semilogx(kh, (m[f"fixed_h_nl_{i}"] / m["fixed_h_nl_0"] - 1) * 100, color=col,
                    lw=1.0, ls=(0, (3, 2)))
        ax.plot([8.0], [-8 * fnu * 100], marker="<", color=col, ms=6, mec=ps.SURFACE, mew=0.6)
        knr = relic.k_nonrelativistic(mass / 3, cfg.T_CMB_K * (4 / 11) ** (1 / 3),
                                      float(cfg.OMEGA_M_REF))
        ax.plot([knr], [0.6], marker="v", color=col, ms=5, mec=ps.SURFACE, mew=0.6)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.set_xlim(1e-4, 10)
    ax.set_xlabel("wavenumber k (h/Mpc)")
    ax.set_ylabel("change in $P(k)$ (%)")
    ax.set_title("b  Neutrinos erase small-scale structure")
    ax.legend(loc="lower left", fontsize=7.0)
    ax.text(0.98, 0.98, "solid: linear   dashed: non-linear (HMcode 2020)\n"
            "◀ the −8 $f_\\nu$ rule   ▼ free-streaming scale $k_\\mathrm{nr}$", transform=ax.transAxes,
            ha="right", va="top", fontsize=7, color=ps.INK_SECONDARY)
    ps.caption(fig, "Calculation with CAMB: three degenerate neutrinos at fixed $\\Omega_m$, "
               "$\\omega_b$, $h$ and $A_s$ around Planck 2018's best fit; $f_\\nu = \\Omega_\\nu/\\Omega_m$.",
               y=-0.03)
    return ps.save(fig, FIG / "07_matter_power")


# ---------------------------------------------------------------------------
# 04 Counting neutrino species
# ---------------------------------------------------------------------------

#: Our rungs: (key, label, colour, dash). Three hues by data (see plotstyle).
NEFF_STYLE = [
    ("tt", "TT (D3)", ps.AQUA, "-"),
    ("plik_tt", "TT + low-ℓ", ps.ORANGE, "-"),
    ("plik", "TT,TE,EE + low-ℓ", ps.ORANGE, (0, (4, 2))),
    ("planck18", "+ lens + BAO 2018", ps.BLUE, (0, (4, 2))),
    ("cmb_desi", "+ ACT lens + DESI", ps.BLUE, "-"),
]


def _smooth_profile(x, y, n=400):
    """A smooth curve through a profile's points (shape-preserving cubic)."""
    from scipy.interpolate import PchipInterpolator

    x, y = np.asarray(x, float), np.asarray(y, float)
    o = np.argsort(x)
    xx = np.linspace(x[o][0], x[o][-1], n)
    return xx, PchipInterpolator(x[o], y[o])(xx)


def fig_neff(r: dict) -> list:
    import fetch_data as fd

    fits = r["fits"]["neff"]
    neff_sm = r["decoupling"]["runs"][cfg.DECOUPLING_HEADLINE]["neff"]
    fig = plt.figure(figsize=(14.6, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.1, 1], wspace=0.62)
    ax = fig.add_subplot(gs[0])
    for key, lab, col, ls in NEFF_STYLE:
        if key not in fits:
            continue
        xx, yy = _smooth_profile(fits[key]["profile"]["values"], fits[key]["profile"]["dchi2"])
        ax.plot(xx, yy, color=col, ls=ls, label=lab)
    for lev in (1, 4):
        ax.axhline(lev, color=ps.AXIS, lw=0.8, ls=":")
    ax.axvline(neff_sm, color=ps.INK, lw=1, label=f"Standard Model, {neff_sm:.3f}")
    ax.set_xlim(1.6, 5.6)
    ax.set_ylim(0, 10)
    ax.set_xlabel("$N_\\mathrm{eff}$")
    ax.set_ylabel("Δχ² (profile likelihood)")
    ax.legend(loc="lower right", fontsize=6.6)
    ax.set_title("a  How many neutrino species?")

    ax = fig.add_subplot(gs[1])
    rows = []
    forest = [("tt", "TT, D3's likelihood", ps.AQUA), ("plik_tt_tau", "Planck TT, no low-ℓ", ps.AQUA),
              ("plik_tt", "Planck TT + low-ℓ", ps.ORANGE), ("plik", "TT,TE,EE + low-ℓ", ps.ORANGE),
              ("planck18", "+ lensing + BAO (2018 data)", ps.BLUE),
              ("cmb_desi", "+ ACT lensing + DESI DR2", ps.BLUE)]
    shown = set()
    for key, lab, col in forest:
        if key not in fits:
            continue
        iv = fits[key]["interval68"]
        rows.append((f"ours · {lab}", iv["best"], iv["minus"], iv["plus"], col, "o"))
        cmp_ = fits[key].get("compare")
        if cmp_ and cmp_["value"] and cmp_["label"] not in shown and key not in ("tt", "plik_tt_tau"):
            m, e = cmp_["value"]
            rows.append((cmp_["label"], m, e, e, ps.INK, "s"))
            shown.add(cmp_["label"])
    for lab in ("ACT DR6 P-ACT-LB", "SPT-3G CMB-SPA"):
        m, e = cfg.NEFF_PUBLISHED[lab]
        rows.append((lab, m, e, e, ps.INK, "s"))
    m, e = cfg.NEFF_BBN_YEH2022
    rows.append(("BBN alone (Yeh et al. 2022)", m, e, e, ps.MUTED, "D"))
    y = np.arange(len(rows))[::-1]
    for yy, (lab, m, lo, hi, col, mk) in zip(y, rows):
        ax.errorbar([m], [yy], xerr=[[lo], [hi]], fmt=mk, color=col, ms=5, capsize=0,
                    elinewidth=1.4, mec=ps.SURFACE if mk == "o" else col, mew=0.6)
    ax.axvline(neff_sm, color=ps.INK, lw=1)
    ax.axvline(3.0, color=ps.MUTED, lw=0.7, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([row[0] for row in rows], fontsize=7)
    ax.set_xlim(2.0, 4.0)
    ax.set_xlabel("$N_\\mathrm{eff}$ (68% intervals)")
    ax.grid(axis="y", visible=False)
    ax.set_title("b  Ours against the published")

    ax = fig.add_subplot(gs[2])
    try:
        ch = fd.load_planck_chain("nnu_TTTEEE_lens_BAO")
        h0 = ch["H0*"] if "H0*" in ch else ch["H0"]
        _contour(ax, ch["nnu"], h0, ch["weight"], ps.INK, "Planck 2018 chain, + lensing + BAO")
    except FileNotFoundError:
        pass
    try:
        ch = fd.load_act_chain()
        _contour(ax, ch["nnu"], ch["H0"], ch["weight"], ps.MUTED, "ACT DR6 chain (P-ACT-LB)", ls="--")
    except FileNotFoundError:
        pass
    for key, lab, col, ls in NEFF_STYLE:
        if key not in fits:
            continue
        h = fits[key]["H0_at"]
        xs = np.array([float(k) for k in h])
        ys = np.array([h[k] for k in h])
        prof = dict(zip(np.round(fits[key]["profile"]["values"], 6), fits[key]["profile"]["dchi2"]))
        keep = np.array([prof.get(round(v, 6), 99) <= 9.0 for v in xs])
        o = np.argsort(xs)
        xs, ys, keep = xs[o], ys[o], keep[o]
        if keep.sum() >= 2:
            ax.plot(xs[keep], ys[keep], color=col, ls=ls, lw=1.6, marker="o", ms=2.5)
    ax.axhspan(cfg.H0_SH0ES[0] - cfg.H0_SH0ES[1], cfg.H0_SH0ES[0] + cfg.H0_SH0ES[1],
               color=ps.lighten(ps.ORANGE, 0.8), lw=0)
    ax.text(2.05, cfg.H0_SH0ES[0] + 1.25, "local H0 (SH0ES)", fontsize=7, color=ps.INK_SECONDARY)
    ax.plot([], [], color=ps.INK_SECONDARY, marker="o", ms=2.5, lw=1.4,
            label="ours (colours as in a)")
    ax.set_xlim(2.0, 4.6)
    ax.set_ylim(58, 78)
    ax.set_xlabel("$N_\\mathrm{eff}$")
    ax.set_ylabel("$H_0$ (km/s/Mpc)")
    ax.legend(loc="lower right", fontsize=6.6)
    ax.set_title("c  Extra radiation raises $H_0$")
    ps.caption(fig, "Our fits: CAMB with Planck 2018 data (D3's binned TT, or plik_lite with Planck's "
               "low-ℓ likelihoods), lensing and BAO as labelled. In c each line follows the best-fit H0 "
               "as N_eff is stepped, where the fit allows it (Δχ² ≤ 9); contours: 68% and 95%.", y=-0.03)
    return ps.save(fig, FIG / "04_counting_neutrinos")


def _contour(ax, x, y, w, color, label, ls="-", levels=(0.68, 0.95), bins=60):
    """Weighted 2D histogram contours enclosing the given probability mass."""
    from scipy.ndimage import gaussian_filter

    h, xe, ye = np.histogram2d(x, y, bins=bins, weights=w)
    h = gaussian_filter(h, 1.2)
    flat = np.sort(h.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    lev = sorted(float(flat[np.searchsorted(cum, q)]) for q in levels)
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    ax.contour(xc, yc, h.T, levels=lev, colors=[color], linestyles=[ls], linewidths=[0.9, 1.3])
    ax.plot([], [], color=color, ls=ls, label=label)


# ---------------------------------------------------------------------------
# 05 The neutrino mass plot
# ---------------------------------------------------------------------------


def fig_mass_plot(r: dict) -> list:
    import masses as ms

    nf = r["masses"]["nufit"]["TByes"]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), gridspec_kw={"wspace": 0.26})
    ml = np.geomspace(1e-4, 2.0, 600)
    limits = [("Planck 2018 + BAO", cfg.MNU_PUBLISHED["Planck TT,TE,EE+lowE+lensing+BAO"]),
              ("DESI DR1 + CMB", cfg.MNU_PUBLISHED["DESI DR1 + CMB"]),
              ("DESI DR2 + CMB", cfg.MNU_PUBLISHED["DESI DR2 + CMB (CamSpec)"])]
    ours = r.get("fits", {}).get("mnu", {}).get("cmb_desi")
    for order, col in (("NO", ps.C_NORMAL), ("IO", ps.C_INVERTED)):
        rng = {k: (nf[order][k]["best"], nf[order][k]["lo3"], nf[order][k]["hi3"])
               for k in ("dm21", "dm3l", "s12sq", "s13sq")}
        lo, hi = ms.mass_band(ml, rng, order, "sum")
        best = ms.sum_masses(ml, rng["dm21"][0], rng["dm3l"][0], order)
        ax = axes[0]
        ax.fill_between(ml, lo, hi, color=ps.lighten(col, 0.6), lw=0)
        ax.plot(ml, best, color=col, lw=1.4,
                label=f"{'normal' if order == 'NO' else 'inverted'} ordering (NuFIT 6.0, 3σ)")
        lo, hi = ms.mass_band(ml, rng, order, "mbeta")
        tot = ms.sum_masses(ml, rng["dm21"][0], rng["dm3l"][0], order)
        mb = ms.m_beta(ml, rng["dm21"][0], rng["dm3l"][0], order, rng["s12sq"][0], rng["s13sq"][0])
        ax = axes[1]
        ax.fill_between(tot, lo, hi, color=ps.lighten(col, 0.6), lw=0)
        ax.plot(tot, mb, color=col, lw=1.4)
    ax = axes[0]
    for i, (lab, lim) in enumerate(limits):
        ax.axhline(lim, color=ps.INK, lw=0.9, ls=[(0, (1, 1.5)), (0, (4, 2)), "-"][i],
                   label=f"{lab} < {lim:g} eV")
    if ours and ours.get("limit95_flat_ge_0"):
        lim = ours["limit95_flat_ge_0"]
        ax.axhspan(lim, 3.0, color=ps.WASH, lw=0, zorder=0)
        ax.axhline(lim, color=ps.BLUE, lw=1.6, label=f"ours, CMB + DESI DR2 < {lim:.3f} eV")
        ax.text(0.04, 0.0715, "↑ excluded at 95% by our fit", fontsize=7, color=ps.INK_SECONDARY,
                ha="left", va="bottom")
    if ours and ours.get("limit95_flat_ge_NO"):
        # Each ordering on its own (prior starting at its minimum sum): the limit that
        # bounds the lightest mass. Reading it off the Σm_ν ≥ 0 line would overstate it.
        for order, col in (("NO", ps.C_NORMAL), ("IO", ps.C_INVERTED)):
            lim = ours[f"limit95_flat_ge_{order}"]
            m_l = ms.lightest_from_sum(lim, nf[order]["dm21"]["best"], nf[order]["dm3l"]["best"],
                                       order)
            ax.plot([m_l], [lim], "o", ms=6.5, mfc=ps.SURFACE, mec=col, mew=1.7, zorder=5)
        ax.plot([], [], "o", ms=6.5, mfc=ps.SURFACE, mec=ps.INK, mew=1.4,
                label="ours, each ordering on its own")
    kat = r["masses"]["NO"]["katrin_lightest_limit_eV"]
    ax.axvspan(kat, 2.0, color=ps.lighten(ps.MUTED, 0.75), lw=0, zorder=0)
    ax.text(kat * 1.08, 0.03, "excluded by\nKATRIN\n($m_\\beta$ < 0.45 eV)", fontsize=7,
            color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-4, 2.0)
    ax.set_ylim(0.02, 3.0)
    ax.set_xlabel("mass of the lightest neutrino (eV)")
    ax.set_ylabel("sum of the three masses $\\Sigma m_\\nu$ (eV)")
    ax.legend(loc="upper left", fontsize=7.0)
    ax.set_title("a  Cosmology presses on the oscillation minimum")
    ax = axes[1]
    ax.axhline(cfg.KATRIN_MBETA_LIMIT, color=ps.MUTED, lw=1.2)
    ax.text(0.14, cfg.KATRIN_MBETA_LIMIT * 1.06, "KATRIN 2025: $m_\\beta$ < 0.45 eV (90%)", fontsize=7,
            color=ps.INK_SECONDARY)
    ax.axhline(cfg.KATRIN_FINAL_SENSITIVITY, color=ps.MUTED, lw=0.8, ls=(0, (4, 2)))
    ax.text(0.14, cfg.KATRIN_FINAL_SENSITIVITY * 0.86, "KATRIN's final goal, ~0.3 eV", fontsize=7,
            color=ps.INK_SECONDARY, va="top")
    styles = [(0, (1, 1.5)), (0, (4, 2)), "-"]
    for i, (lab, lim) in enumerate(limits):
        ax.axvline(lim, color=ps.INK, lw=0.9, ls=styles[i], label=f"{lab} < {lim:g} eV")
    if ours and ours.get("limit95_flat_ge_0"):
        ax.axvline(ours["limit95_flat_ge_0"], color=ps.BLUE, lw=1.6,
                   label=f"ours, CMB + DESI DR2 < {ours['limit95_flat_ge_0']:.3f} eV")
    ax.legend(loc="lower right", fontsize=6.8, title="cosmology, 95%", title_fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.05, 3.0)
    ax.set_ylim(0.005, 2.0)
    ax.set_xlabel("sum of the masses $\\Sigma m_\\nu$ (eV): what cosmology weighs")
    ax.set_ylabel("electron-neutrino mass $m_\\beta$ (eV): what KATRIN weighs")
    ax.set_title("b  Two scales, two experiments")
    each = ""
    if ours and ours.get("limit95_flat_ge_NO"):
        each = (f"\nOpen circles: our limits for each ordering on its own (prior starting at its "
                f"minimum sum), {ours['limit95_flat_ge_NO']:.3f} and "
                f"{ours['limit95_flat_ge_IO']:.3f} eV; these, not the lines, bound the lightest mass.")
    ps.caption(fig, "Bands: NuFIT 6.0 3σ ranges of the mass splittings and mixing angles (with SK "
               "atmospheric data), from NuFIT's Δχ² tables. Cosmological limits: 95%, flat prior "
               "Σm_ν ≥ 0, three degenerate masses; DESI DR1's 0.072 eV was later corrected to 0.082."
               + each, y=-0.03)
    return ps.save(fig, FIG / "05_neutrino_mass_plot")


# ---------------------------------------------------------------------------
# 06 Weighing neutrinos: our profile likelihoods
# ---------------------------------------------------------------------------

MNU_STYLE = [
    ("tt", "TT, D3's likelihood (+ τ prior)", ps.AQUA, "-"),
    ("tt_bao", "TT (D3) + DESI BAO", ps.AQUA, (0, (4, 2))),
    ("plik", "Planck TT,TE,EE + low-ℓ", ps.ORANGE, "-"),
    ("planck18", "Planck + lensing + BAO (Planck 2018's data)", ps.ORANGE, (0, (4, 2))),
    ("cmb_desi", "Planck + ACT lensing + DESI DR2 BAO", ps.BLUE, "-"),
    ("cmb_desi_sn_w0wa", "… + DES-SN5YR, w0wa dark energy", ps.BLUE, (0, (4, 2))),
]


#: Rungs whose profile is close enough to a parabola for the negative-mass extension to mean
#: something (the CMB + BAO combinations; temperature-only profiles are far from Gaussian).
EXTRAPOLATE = ("tt_bao", "planck18", "cmb_desi", "cmb_desi_sn_w0wa")


def fig_mass_bound(r: dict) -> list:
    import fetch_data as fd
    from fitting import parabola_fit

    fits = r["fits"]["mnu"]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.8), gridspec_kw={"wspace": 0.24})
    ax = axes[0]
    for key, lab, col, ls in MNU_STYLE:
        if key not in fits:
            continue
        x = np.array(fits[key]["profile"]["values"])
        y = np.array(fits[key]["profile"]["dchi2"])
        o = np.argsort(x)
        ax.plot(x[o], y[o], color=col, ls=ls, marker="o", ms=2.5, label=lab)
        par = fits[key]["parabola"]
        if key in EXTRAPOLATE and par.get("mu") is not None and par.get("sigma") is not None \
                and np.isfinite(par["mu"]) and par["mu"] < 0:
            # the parabola through the profile, continued to negative masses; its minimum
            # lies below zero when the data prefer less clustering than any mass allows
            xx = np.linspace(max(par["mu"] - 1.5 * par["sigma"], -0.3), 0, 80)
            ax.plot(xx, ((xx - par["mu"]) / par["sigma"]) ** 2 + par["offset"], color=col,
                    lw=1.0, ls=":")
            ax.plot([par["mu"]], [par["offset"]], marker="v", color=col, ms=5)
    ax.axvspan(-0.3, 0.0, color=ps.WASH, lw=0, zorder=0)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.text(-0.29, 9.6, "negative masses are unphysical:\ndotted = the parabola through each\n"
            "profile, continued (▼ its minimum)", fontsize=6.6, color=ps.INK_SECONDARY, va="top")
    ax.axvline(cfg.MIN_SUM_NO, color=ps.C_NORMAL, lw=0.8, ls=":")
    ax.axvline(cfg.MIN_SUM_IO, color=ps.C_INVERTED, lw=0.8, ls=":")
    ax.text(cfg.MIN_SUM_NO - 0.006, -1.45, "NO\nmin", fontsize=6.4, color=ps.INK_SECONDARY,
            ha="right", va="bottom")
    ax.text(cfg.MIN_SUM_IO + 0.006, -1.45, "IO\nmin", fontsize=6.4, color=ps.INK_SECONDARY,
            ha="left", va="bottom")
    ax.set_xlim(-0.3, 1.0)
    ax.set_ylim(-1.6, 10)
    ax.set_xlabel("$\\Sigma m_\\nu$ (eV)")
    ax.set_ylabel("Δχ² (profile likelihood)")
    ax.legend(loc="lower right", fontsize=6.6)
    ax.set_title("a  How each data set weighs the neutrinos")
    ax = axes[1]
    grid = np.linspace(0, 0.4, 801)
    for key, lab, col, ls in MNU_STYLE[4:]:
        if key not in fits:
            continue
        x = np.array(fits[key]["profile"]["values"])
        y = np.array(fits[key]["profile"]["dchi2"])
        o = np.argsort(x)
        from scipy.interpolate import PchipInterpolator

        like = np.exp(-0.5 * PchipInterpolator(x[o], y[o], extrapolate=False)(grid))
        like = np.nan_to_num(like)
        like /= np.trapezoid(like, grid)
        ax.plot(grid, like, color=col, ls=ls, label=f"ours: {lab}")
    for run, col, ls, lab in ((fd.DESI_CHAIN_RUNS[1], ps.INK, "-", "DESI chain: DESI + CMB (plik)"),
                              (fd.DESI_CHAIN_RUNS[2], ps.INK, (0, (4, 2)),
                               "DESI chain: + DESY5, w0wa")):
        try:
            ch = fd.load_desi_chain(run)
        except FileNotFoundError:
            continue
        h, e = np.histogram(ch["mnu"], bins=np.linspace(0, 0.4, 41), weights=ch["weight"],
                            density=True)
        ax.step(0.5 * (e[1:] + e[:-1]), h, where="mid", color=col, ls=ls, lw=0.9, label=lab)
    ax.axvline(cfg.MIN_SUM_NO, color=ps.C_NORMAL, lw=0.8, ls=":")
    ax.axvline(cfg.MIN_SUM_IO, color=ps.C_INVERTED, lw=0.8, ls=":")
    ax.set_xlim(0, 0.4)
    ax.set_ylim(0, None)
    ax.set_xlabel("$\\Sigma m_\\nu$ (eV)")
    ax.set_ylabel("posterior density, flat prior $\\Sigma m_\\nu$ ≥ 0 (1/eV)")
    ax.legend(loc="upper right", fontsize=6.8)
    ax.set_title("b  Our posteriors against DESI's own chains")
    ps.caption(fig, "Profile likelihoods: at each Σm_ν (three degenerate masses) every other "
               "parameter is re-fitted with CAMB. Posteriors in b: exp(−Δχ²/2) with a flat prior; "
               "histograms: DESI DR2's released MCMC chains (Planck's full plik and low-ℓ likelihoods).",
               y=-0.03)
    return ps.save(fig, FIG / "06_weighing_neutrinos")


# ---------------------------------------------------------------------------
# 00 Scorecard
# ---------------------------------------------------------------------------

SCORE_NEFF = [("plik_tt", "Planck TT + low-ℓ"), ("plik", "Planck TT,TE,EE + low-ℓ"),
              ("planck18", "+ lensing + BAO (Planck 2018)"),
              ("cmb_desi", "+ ACT lensing + DESI DR2 BAO")]
SCORE_MNU = [("tt", "TT, D3's likelihood"), ("tt_bao", "TT (D3) + DESI BAO"),
             ("plik_tt", "Planck TT + low-ℓ"), ("plik", "Planck TT,TE,EE + low-ℓ"),
             ("plik_lens18", "+ Planck lensing"), ("planck18", "+ lensing + BAO (Planck 2018)"),
             ("cmb_desi", "+ ACT lensing + DESI DR2 BAO"),
             ("cmb_desi_sn_w0wa", "+ DES-SN5YR, w0wa dark energy")]
#: Rows whose comparison is only the nearest published result (no one else uses D3's likelihood).
SCORE_NEAREST_ONLY = ("tt", "tt_bao")


def fig_scorecard(r: dict) -> list:
    fits = r["fits"]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.6), gridspec_kw={"wspace": 0.55})
    ax = axes[0]
    rows = [(k, lab) for k, lab in SCORE_NEFF if k in fits["neff"]]
    y = np.arange(len(rows))[::-1]
    for yy, (k, lab) in zip(y, rows):
        s_ = fits["neff"][k]
        iv = s_["interval68"]
        ax.errorbar([iv["best"]], [yy + 0.13], xerr=[[iv["minus"]], [iv["plus"]]], fmt="o",
                    color=ps.BLUE, ms=6, mec=ps.SURFACE, mew=0.8, capsize=0, elinewidth=1.6)
        if s_.get("compare") and s_["compare"]["value"]:
            m, e = s_["compare"]["value"]
            ax.errorbar([m], [yy - 0.13], xerr=[[e], [e]], fmt="s", color=ps.INK, ms=5,
                        capsize=0, elinewidth=1.2)
            ax.text(m + e + 0.03, yy - 0.13, s_["compare"]["label"], fontsize=6.4,
                    color=ps.INK_SECONDARY, va="center", ha="left")
    neff = r["decoupling"]["runs"][cfg.DECOUPLING_HEADLINE]["neff"]
    ax.axvline(neff, color=ps.INK, lw=1)
    ax.text(neff + 0.02, len(rows) - 0.45, f"our decoupling\ncalculation: {neff:.4f}", fontsize=7,
            color=ps.INK_SECONDARY, va="top")
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in rows], fontsize=7.6)
    ax.set_xlim(2.3, 4.6)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("$N_\\mathrm{eff}$ (68%)")
    ax.plot([], [], "o", color=ps.BLUE, label="ours")
    ax.plot([], [], "s", color=ps.INK, label="published (same data)")
    ax.legend(loc="upper right", fontsize=7)
    ax.set_title("a  How many neutrino species?")

    ax = axes[1]
    rows = [(k, lab) for k, lab in SCORE_MNU if k in fits["mnu"]]
    y = np.arange(len(rows))[::-1]
    for yy, (k, lab) in zip(y, rows):
        s_ = fits["mnu"][k]
        lim = s_.get("limit95_flat_ge_0")
        if lim:
            ax.plot([lim], [yy + 0.13], "o", color=ps.BLUE, ms=6, mec=ps.SURFACE, mew=0.8)
            ax.plot([0.02, lim], [yy + 0.13, yy + 0.13], color=ps.BLUE, lw=1.6)
            ax.text(lim * 1.08, yy + 0.2, f"{lim:.3g}", fontsize=6.8, color=ps.INK, va="bottom")
        if s_.get("compare") and s_["compare"]["value"]:
            v = s_["compare"]["value"]
            # No published analysis uses exactly D3's likelihood: show the nearest one, muted.
            near = k in SCORE_NEAREST_ONLY
            col = ps.MUTED if near else ps.INK
            ax.plot([v], [yy - 0.13], "s", color=col, mfc=ps.SURFACE if near else col, ms=5)
            ax.plot([0.02, v], [yy - 0.13, yy - 0.13], color=col, lw=1.0,
                    ls=(0, (3, 2)) if near else "-")
            ax.text(v * 1.08, yy - 0.2, ("nearest published: " if near else "")
                    + s_["compare"]["label"], fontsize=6.2, color=ps.INK_SECONDARY, va="top",
                    bbox=dict(facecolor=ps.SURFACE, edgecolor="none", pad=0.6), zorder=4)
    ax.axvline(cfg.MIN_SUM_NO, color=ps.C_NORMAL, lw=0.8, ls=":")
    ax.axvline(cfg.MIN_SUM_IO, color=ps.C_INVERTED, lw=0.8, ls=":")
    top = len(rows) - 0.42
    ax.text(cfg.MIN_SUM_NO * 0.97, top, "normal-ordering\nminimum", fontsize=6.4,
            color=ps.INK_SECONDARY, ha="right", va="bottom")
    ax.text(cfg.MIN_SUM_IO * 1.03, top, "inverted-ordering\nminimum", fontsize=6.4,
            color=ps.INK_SECONDARY, ha="left", va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(0.02, 3.0)
    ax.set_ylim(-0.6, len(rows) + 0.35)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in rows], fontsize=7.6)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("95% upper limit on $\\Sigma m_\\nu$ (eV), flat prior $\\Sigma m_\\nu$ ≥ 0")
    ax.set_title("b  How heavy can they be?")
    ps.caption(fig, "Blue: our fits (CAMB + the data on each row). Ink: the published result for the "
               "same data (Planck 2018 VI and its parameter tables; DESI DR2 Results II; Elbers et al. "
               "2025). Open grey squares: the nearest\npublished result, for D3's likelihood, which no "
               "published analysis uses. Our Planck rungs use plik_lite and Planck's low-ℓ likelihoods.",
               y=-0.03)
    return ps.save(fig, FIG / "00_scorecard")


FIGURES = {
    "00": fig_scorecard,
    "01": fig_relic,
    "02": fig_decoupling,
    "03": fig_neff_spectra,
    "04": fig_neff,
    "05": fig_mass_plot,
    "06": fig_mass_bound,
    "07": fig_matter_power,
}


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ps.use_style()
    r = _results()
    todo = argv or list(FIGURES)
    for key in todo:
        written = FIGURES[key](r)
        print(f"  {key}: " + ", ".join(Path(w).name for w in written))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
