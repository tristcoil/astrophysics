"""Every figure of the project, from the cached analysis results.

    python figures.py            # needs run.py (or its caches) to have been run once
"""

from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from scipy import stats

import config as cfg
import dimuon
import electroweak as ew
import fourlepton as fl
import higgs_search
import higgsfit as hf
import plotstyle as ps
import stability as st

ps.use_style()
FIG = cfg.FIGURE_DIR
C_DATA = ps.INK
C_ZZ = ps.BLUE
C_ZX = ps.AQUA
C_H = ps.ORANGE
C_OTHER = ps.VIOLET
WASH = {"stable": "#dbe8f9", "metastable": "#ecebe6", "unstable": "#fbe1d6",
        "nonpert": "#d9d8d2"}


def _save(fig, name):
    return ps.save(fig, str(FIG / name))


def _step(ax, edges, values, **kw):
    ax.stairs(values, edges, **kw)


def _poisson_err(n):
    """Central 68% Garwood interval for Poisson counts (asymmetric bars)."""
    n = np.asarray(n, dtype=float)
    lo = np.where(n > 0, n - stats.chi2.ppf(0.15865, 2 * n) / 2, 0.0)
    hi = stats.chi2.ppf(0.84135, 2 * (n + 1)) / 2 - n
    return lo, hi


# ---------------------------------------------------------------------------
# 1. The dimuon spectrum
# ---------------------------------------------------------------------------

LABELS = [("η", 0.5479, 1.0), ("ρ,ω", 0.78, 1.0), ("φ", 1.0195, 1.0), ("J/ψ", 3.0969, 1.0),
          ("ψ(2S)", 3.686, 1.18), ("ϒ(1S,2S,3S)", 9.9, 1.0), ("Z", 91.19, 1.0)]


def fig_dimuon_spectrum(d, R):
    sp = dimuon.spectrum(d, nbins=1800)
    e, y = sp["edges"], sp["per_gev"]
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    ax.stairs(y, e, color=C_DATA, linewidth=0.8, baseline=None)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.25, 300)
    ax.set_ylim(20, 3e8)
    ax.set_xlabel("Invariant mass of the muon pair, m(μμ)  [GeV]")
    ax.set_ylabel("Pairs per GeV")
    mids = np.sqrt(e[:-1] * e[1:])
    for lab, m, shift in LABELS:
        j = np.argmin(np.abs(mids - m))
        top = y[max(j - 3, 0):j + 4].max()
        ax.annotate(lab, (m, top * 1.25), xytext=(m * shift, top * 6), ha="center", fontsize=9,
                    color=ps.INK, arrowprops={"arrowstyle": "-", "color": ps.MUTED,
                                               "linewidth": 0.6})
    j = np.argmin(np.abs(mids - 28))
    ax.annotate("trigger threshold\n(not a particle)", (28, y[j]), xytext=(17, y[j] * 0.04),
                fontsize=7.5, color=ps.MUTED, ha="center",
                arrowprops={"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6})
    n = int(d["n_opposite_charge"])
    ax.set_title(f"The particle zoo in one plot: {n / 1e6:.1f} million muon pairs, "
                 "CMS 2012 open data")
    ps.caption(fig, "Events with exactly two muons of opposite charge (record 12341; 61.5 million "
               "events, Run 2012B+C, 11.6 fb⁻¹ at √s = 8 TeV). Each bump is a particle that decays "
               "to μ⁺μ⁻; the smooth continuum is mostly virtual photons (Drell–Yan) and heavy-quark "
               "decays.", y=-0.03)
    return _save(fig, "01_dimuon_spectrum")


# ---------------------------------------------------------------------------
# 2. Masses from the bumps, and the muon momentum scale
# ---------------------------------------------------------------------------


def _fit_panel(ax, fr, title, *, log=False, xlab="m(μμ)  [GeV]", scale=1e6,
               ylab="Million pairs per GeV"):
    e = np.array(fr["edges"])
    n = np.array(fr["counts"])
    mu = np.array(fr["expected"])
    w = np.diff(e)
    mid = 0.5 * (e[:-1] + e[1:])
    ax.errorbar(mid, n / w / scale, yerr=np.sqrt(n) / w / scale, fmt="o", ms=1.6, color=C_DATA,
                elinewidth=0.5, label="data")
    ax.plot(mid, mu / w / scale, color=C_H, linewidth=1.3, label="fit")
    ax.set_title(title, fontsize=9.5)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if log:
        ax.set_yscale("log")


def fig_resonance_masses(d, R):
    fits, cal = R["dimuon_fits"], R["higgs"]["calibration"]
    fig = plt.figure(figsize=(11, 7.2))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.3)
    ax = fig.add_subplot(gs[0, 0])
    _fit_panel(ax, fits["phi"], f"φ meson: {fits['phi']['peaks']['phi']['mass'] * 1000:.1f} MeV "
               f"(PDG {fits['phi']['peaks']['phi']['pdg_mass'] * 1000:.1f})")
    ax = fig.add_subplot(gs[0, 1])
    jp = fits["jpsi"]["peaks"]["jpsi"]
    _fit_panel(ax, fits["jpsi"], f"J/ψ: {jp['mass'] * 1000:.1f} MeV (PDG {jp['pdg_mass'] * 1000:.1f})")
    # Upsilons: barrel vs endcap resolution
    ax = fig.add_subplot(gs[0, 2])
    cats = dimuon.categories(d)
    edges = np.linspace(8.6, 11.2, 261)
    for name, c, lab in (("barrel", C_ZZ, "both muons central, |η| < 0.9"),
                         ("endcap", C_OTHER, "both forward, |η| > 1.2")):
        h = np.histogram(d["mass"][cats[name]], edges)[0].astype(float)
        ax.stairs(h / h.max(), edges, color=c, linewidth=1.1, label=lab, baseline=None)
    ax.set_ylim(0, 1.32)
    ax.set_title("ϒ(1S), ϒ(2S), ϒ(3S): resolution decides", fontsize=9.5)
    ax.set_xlabel("m(μμ)  [GeV]")
    ax.set_ylabel("Pairs (peak = 1)")
    ax.legend(loc="upper right", fontsize=7)
    # Z
    ax = fig.add_subplot(gs[1, 0])
    zb = R["dimuon_z_barrel"]
    _fit_panel(ax, zb, f"Z boson, central muons: {zb['mass']:.2f} GeV\n(PDG {zb['pdg_mass']:.2f})",
               scale=1e3, ylab="Thousand pairs per GeV")
    # the scale plot
    ax = fig.add_subplot(gs[1, 1:])
    rows = []
    for key, fr in fits.items():
        for pk, p in fr["peaks"].items():
            if pk in ("eta", "omega", "z"):     # the Z: data vs simulation below instead
                continue
            rows.append((p["pdg_mass"], p["ratio_minus_1"], p["ratio_err"], pk))
    offsets = {"phi": (0, 9), "jpsi": (-14, 9), "psi2s": (6, 9), "upsilon1s": (-22, 9),
               "upsilon2s": (-8, -16), "upsilon3s": (6, 9)}
    names = {"phi": "φ", "jpsi": "J/ψ", "psi2s": "ψ(2S)", "upsilon1s": "ϒ(1S)",
             "upsilon2s": "ϒ(2S)", "upsilon3s": "ϒ(3S)"}
    for m0, r, e_, k in rows:
        ax.errorbar(m0, r * 1e3, yerr=e_ * 1e3, fmt="o", color=C_DATA, ms=4)
        ax.annotate(names.get(k, k), (m0, r * 1e3), xytext=offsets.get(k, (4, 6)),
                    textcoords="offset points", fontsize=7.5, color=ps.INK_SECONDARY,
                    ha="center")
    for reg, mk in (("barrel", "s"), ("endcap", "D")):
        z = cal["mu"][reg]["data"]
        zs = cal["mu"][reg]["sim"]
        ax.errorbar(91.1876 * 1.04 if reg == "endcap" else 91.1876 * 0.96,
                    (z["mass"] / zs["mass"] - 1) * 1e3, yerr=z["mass_err"] / zs["mass"] * 1e3,
                    fmt=mk, color=C_ZZ, ms=4.5, mfc="white")
    for reg, mk in (("barrel", "s"), ("endcap", "D")):
        z = cal["e"][reg]["data"]
        zs = cal["e"][reg]["sim"]
        ax.errorbar(91.1876 * (1.08 if reg == "endcap" else 0.92), (z["mass"] / zs["mass"] - 1) * 1e3,
                    yerr=np.hypot(z["mass_err"], zs["mass_err"]) / zs["mass"] * 1e3, fmt=mk,
                    color=C_H, ms=4.5, mfc="white")
    ax.axhline(0, color=ps.AXIS, linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.7, 160)
    ax.set_xlabel("Particle mass  [GeV]")
    ax.set_ylabel("Measured / reference − 1  [per mille]")
    ax.set_title("The momentum scale of the outreach files is off by a few per mille",
                 fontsize=9.5)
    handles = [Line2D([], [], marker="o", ls="", color=C_DATA, label="μμ peak vs PDG 2026"),
               Line2D([], [], marker="s", ls="", color=C_ZZ, mfc="white",
                      label="Z→μμ, data vs simulation (barrel ■, endcap ◆)"),
               Line2D([], [], marker="s", ls="", color=C_H, mfc="white",
                      label="Z→ee, data vs simulation (barrel ■, endcap ◆)")]
    ax.legend(handles=handles, loc="upper left", fontsize=7)
    ax.set_ylim(-4, 13)
    fig.suptitle("Weighing particles with muon pairs", x=0.01, ha="left", y=0.995)
    ps.caption(fig, "Peaks fitted with Gaussian/Crystal Ball shapes (Voigt for φ, Breit–Wigner ⊗ "
               "Crystal Ball for the Z) over smooth backgrounds. Data muons come out ~0.1% light at "
               "J/ψ and ϒ and ~0.2–0.3% light at the Z, electrons ~0.4% (barrel) to 1% (endcap) "
               "heavy: the files carry no momentum calibration.", y=0.02)
    return _save(fig, "02_resonance_masses")


# ---------------------------------------------------------------------------
# 3. Forward-backward asymmetry
# ---------------------------------------------------------------------------


def fig_afb(d, R):
    edges = np.geomspace(40, 300, 41)
    mid = np.sqrt(edges[:-1] * edges[1:])
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.axvspan(70, 100, color=ps.EXCLUDED, linewidth=0)
    ax.text(71.5, 0.195, "Z decays that radiated a photon\nmigrate down into this band",
            ha="left", fontsize=7.2, color=ps.MUTED)
    for ymin, c, lab in ((0.0, C_ZZ, "all pairs"), (1.0, C_OTHER, "pairs with |y| > 1")):
        r = dimuon.afb(d, edges, ymin=ymin)
        f = dimuon.afb_with_dilution(r)
        ax.errorbar(mid, r["afb"], yerr=r["err"], fmt="o", ms=3.5, color=c, label=lab)
        hi = (f["curve_m"] >= 100) & (f["curve_m"] < 300)
        lo = (f["curve_m"] > 40) & (f["curve_m"] <= 100)
        ax.plot(f["curve_m"][hi], f["curve"][hi], color=c, linewidth=1.1, ls="--")
        ax.plot(f["curve_m"][lo], f["curve"][lo], color=c, linewidth=0.8, ls=":", alpha=0.6)
        R.setdefault("afb_fits", {})[lab] = {k: f[k] for k in ("D", "D_err", "chi2", "ndof",
                                                              "zero_crossing_prediction")}
    ax.axhline(0, color=ps.AXIS, linewidth=0.8)
    ax.axvline(91.19, color=ps.MUTED, linewidth=0.7, ls=":")
    ax.text(92.5, -0.06, "Z mass", fontsize=7.5, color=ps.MUTED)
    ax.set_xscale("log")
    ax.set_xlim(40, 300)
    ticks = [40, 50, 60, 70, 80, 90, 100, 150, 200, 300]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(-0.09, 0.23)
    ax.set_xlabel("m(μμ)  [GeV]")
    ax.set_ylabel("Forward–backward asymmetry  A_FB")
    ax.set_title("The photon and the Z interfere: the asymmetry changes sign at the Z")
    ax.legend(loc="upper left")
    ps.caption(fig, "A_FB = (N_F − N_B)/(N_F + N_B) with the Collins–Soper angle, the quark direction "
               "taken as the pair's boost. Dashed: leading-order quark-level prediction (u:d = 2:1, "
               "sin²θ = 0.2316) times one dilution factor fitted above 110 GeV, dotted where photon "
               "radiation (not modelled here) dominates; its zero is at 89.2 GeV, the data's near "
               "83 GeV.", y=-0.02)
    return _save(fig, "03_forward_backward_asymmetry")


# ---------------------------------------------------------------------------
# 4. The replica and the duplicates
# ---------------------------------------------------------------------------


def _stack_panel(ax, edges, bkg, sig, data, *, title, extra=None):
    ax.stairs(bkg, edges, fill=True, color=ps.RAMPS[ps.BLUE][0], label="Z/γ*, ZZ (simulation)")
    ax.stairs(bkg, edges, color=C_ZZ, linewidth=1.0)
    ax.stairs(bkg + sig, edges, baseline=bkg, color=C_H, linewidth=1.6,
              label="+ Higgs, m_H = 125 GeV (simulation)")
    mid = 0.5 * (edges[:-1] + edges[1:])
    lo, hi = _poisson_err(data)
    show = data > 0
    ax.errorbar(mid[show], data[show], yerr=[lo[show], hi[show]], fmt="o", ms=4, color=C_DATA,
                label="data", zorder=5)
    ax.set_xlabel("Four-lepton mass, m(4ℓ)  [GeV]")
    ax.set_ylabel("Events / 3.06 GeV")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_title(title, fontsize=10)
    if extra:
        extra(ax, mid)


def fig_replica(R):
    pub = __import__("outreach_plots").load()["combined"]["Higgs_mass"]
    ours = higgs_search.replica_histograms("")
    dedup = higgs_search.replica_histograms("", dedupe=True)
    edges = np.array(pub["edges"])
    comb = lambda h, k: sum(h[s]["Higgs_mass"][k] for s in higgs_search.STATES)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharey=True)

    def published_marks(ax, mid):
        n = np.array(pub["data"])
        s = n > 0
        ax.plot(mid[s], n[s], "o", ms=9, mfc="none", mec=C_OTHER, mew=1.0,
                label="published plot (read from its PDF)")

    _stack_panel(axes[0], edges, comb(ours, "background"), comb(ours, "signal"),
                 comb(ours, "data").astype(float),
                 title="Our replica = the published plot",
                 extra=published_marks)
    axes[0].legend(loc="upper right", fontsize=7)
    d = R["higgs"]["duplicates"]
    removed = comb(ours, "data") - comb(dedup, "data")

    def dup_marks(ax, mid):
        s = removed > 0
        ax.plot(mid[s], comb(dedup, "data")[s] + removed[s], "x", ms=7, color=C_H, mew=1.6,
                label=f"counted twice: {int(removed.sum())} events")

    _stack_panel(axes[1], edges, comb(dedup, "background"), comb(dedup, "signal"),
                 comb(dedup, "data").astype(float),
                 title="Each collision counted once", extra=dup_marks)
    axes[1].legend(loc="upper right", fontsize=7)
    b = R["higgs"]["binned"]
    for ax, key in zip(axes, ("as_published", "duplicates_removed")):
        ax.text(0.04, 0.95, f"local significance {b[key]['z_local']:.1f}σ\n"
                f"signal strength μ̂ = {b[key]['mu_hat']:.2f}", transform=ax.transAxes,
                va="top", fontsize=8, color=ps.INK_SECONDARY)
    ps.caption(fig, f"Left: all 9 published plots × 36 bins reproduced exactly (simulation to 10⁻⁵ "
               f"events). Right: {d['duplicated_events']} of the {d['unique']} unique 2e2μ "
               f"collisions are in both the DoubleMuParked and DoubleElectron streams and were "
               f"added twice; {d['duplicated_in_70_180']} of them fall in this plot, "
               f"{d['duplicated_in_peak_bins']} in the two bins at 122–128 GeV.", y=-0.02)
    return _save(fig, "04_replica_and_duplicates")


# ---------------------------------------------------------------------------
# 5. Our four-lepton spectrum
# ---------------------------------------------------------------------------


def fig_four_lepton(R, mc, ssd):
    data = R["higgs"]["data_events"]
    m = np.array([x["m4l"] for x in data])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), gridspec_kw={"width_ratios": [1.35, 1]})
    for ax, (lo, hi, w) in zip(axes, ((70, 180, 2.5), (180, 600, 20))):
        edges = np.arange(lo, hi + 1e-9, w)
        zz = np.histogram(mc["ZZ"]["m4l"], edges, weights=mc["ZZ"]["weight"])[0]
        hh = np.histogram(mc["H"]["m4l"], edges, weights=mc["H"]["weight"])[0]
        # Z+X: the same-sign control region's shape and size (Moyal fit), all channels
        zx_par = R["higgs"]["zx_shape"]
        cdf = stats.moyal.cdf(edges, zx_par["loc"], zx_par["scale"])
        zx = R["higgs"]["zx"]["total_all_masses"] * np.diff(cdf)
        base = np.zeros(len(edges) - 1)
        for vals, col, lab in ((zz, ps.RAMPS[ps.BLUE][0], "ZZ (simulation)"),
                               (zx, "#bfe8d8", "Z + misidentified leptons (same-sign data)")):
            ax.stairs(base + vals, edges, baseline=base, fill=True, color=col, label=lab)
            base = base + vals
        ax.stairs(base, edges, color=C_ZZ, linewidth=0.8)
        ax.stairs(base + hh, edges, baseline=base, color=C_H, linewidth=1.6,
                  label="+ Higgs boson (simulation, μ = 1)")
        n = np.histogram(m, edges)[0].astype(float)
        mid = 0.5 * (edges[:-1] + edges[1:])
        lo_, hi_ = _poisson_err(n)
        s = n > 0
        ax.errorbar(mid[s], n[s], yerr=[lo_[s], hi_[s]], fmt="o", ms=4, color=C_DATA,
                    label="data", zorder=5)
        ax.set_xlim(lo, hi)
        ax.set_xlabel("m(4ℓ)  [GeV]")
        ax.set_ylabel(f"Events / {w:g} GeV")
    axes[0].set_title("Our selection: the Z → 4ℓ peak and the Higgs boson", fontsize=10)
    axes[0].legend(loc="upper right", fontsize=7)
    y = R["higgs"]["yields_onshell_ZZ"]
    axes[1].set_title(f"Pairs of real Z bosons: data {y['data']['all']:.0f}, "
                      f"simulation {y['ZZ']['all']:.0f}", fontsize=10)
    ps.caption(fig, "Each collision counted once, lepton momenta calibrated on the Z peaks, "
                    "simulation normalised as in the outreach analysis (×1.386 for ZZ). The right "
                    "panel checks that normalisation where no Higgs signal is expected.", y=-0.02)
    return _save(fig, "05_four_lepton_mass")


# ---------------------------------------------------------------------------
# 6. Mass, signal strength, significance
# ---------------------------------------------------------------------------


def fig_higgs_fit(R):
    v = R["higgs"]["variants"]["ours_calibrated_with_zx"]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    ax = axes[0]
    pm = v["profile_mass"]
    ax.plot(pm["m"], pm["nll"], color=C_H)
    ax.axhline(1, color=ps.MUTED, linewidth=0.7, ls="--")
    ax.text(129.8, 1.15, "68%", fontsize=7.5, color=ps.MUTED, ha="right")
    mh_pdg, s_pdg = cfg.M_HIGGS
    ax.axvspan(mh_pdg - s_pdg, mh_pdg + s_pdg, color=ps.RAMPS[ps.BLUE][0], linewidth=0)
    ax.text(mh_pdg + 0.25, 5.2, "PDG 2026", fontsize=7.5, color=C_ZZ)
    ax.set_ylim(0, 6)
    ax.set_xlim(120, 130)
    ax.set_xlabel("Higgs mass m_H  [GeV]")
    ax.set_ylabel("−2 Δ ln L")
    syst = R["higgs"]["mass_systematic"]["total"]
    ax.set_title(f"Mass: {v['mh']:.2f} GeV\n+{v['mh_hi'] - v['mh']:.2f} −{v['mh'] - v['mh_lo']:.2f} "
                 f"(stat) ± {syst:.2f} (syst)", fontsize=9.5)
    ax = axes[1]
    pu = v["profile_mu"]
    ax.plot(pu["mu"], pu["nll"], color=C_H)
    for lvl, lab in ((1, "68%"), (3.84, "95%")):
        ax.axhline(lvl, color=ps.MUTED, linewidth=0.7, ls="--")
        ax.text(2.95, lvl + 0.15, lab, fontsize=7.5, color=ps.MUTED, ha="right")
    ax.axvline(1, color=C_ZZ, linewidth=0.8, ls=":")
    ax.text(1.03, 7.2, "Standard Model", fontsize=7.5, color=C_ZZ)
    ax.set_ylim(0, 8)
    ax.set_xlabel("Signal strength μ = σ / σ_SM")
    ax.set_ylabel("−2 Δ ln L")
    ax.set_title(f"Signal strength: μ = {v['mu']:.2f}\n+{v['mu_hi'] - v['mu']:.2f} "
                 f"−{v['mu'] - v['mu_lo']:.2f} (stat + syst)", fontsize=9.5)
    ax = axes[2]
    q = np.array(R["higgs"]["scan_q0"])
    masses = np.array(R["higgs"]["toys"]["scan_masses"])
    p = 0.5 * stats.chi2.sf(q, 1)
    p = np.where(q > 0, p, 0.5)
    ax.plot(masses, p, color=C_H, label="observed")
    t = R["higgs"]["toys"]
    ax.plot([mh_pdg], [stats.norm.sf(t["expected_z_median"])], "o", color=C_ZZ, ms=5,
            label="expected for a Standard Model Higgs")
    for k in range(1, 5):
        ax.axhline(stats.norm.sf(k), color=ps.MUTED, linewidth=0.6, ls="--")
        ax.text(159.5, stats.norm.sf(k) * 1.25, f"{k}σ", fontsize=7.5, color=ps.MUTED, ha="right")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    ax.set_xlim(110, 160)
    ax.set_xlabel("Hypothesised m_H  [GeV]")
    ax.set_ylabel("Local p-value (background only)")
    lee = t["lee"]
    ax.set_title(f"Local {v['z_local_known_mass']:.1f}σ at 125.13 GeV\n"
                 f"global {lee['z_global']:.1f}σ for a search over 110–160 GeV", fontsize=9.5)
    ax.legend(loc="lower left", fontsize=7)
    fig.suptitle("The Higgs boson in 11.6 fb⁻¹ of 2012 CMS data, four-lepton channel",
                 x=0.01, ha="left", y=1.03)
    ps.caption(fig, "Unbinned likelihood over 105–160 GeV in three channels: simulated signal "
               "shapes (resolution widened to match data), ZZ from simulation (±10%), Z+X from "
               "same-sign events (±50%). Global significance: Gross–Vitells upcrossings from "
               f"{t['n_background']:,} pseudo-experiments.", y=-0.04)
    return _save(fig, "06_higgs_fit")


def fig_evidence(R):
    b, v, t = R["higgs"]["binned"], R["higgs"]["variants"], R["higgs"]["toys"]
    rows = [("Published outreach plot", b["as_published"]["z_local"]),
            ("… each collision counted once", b["duplicates_removed"]["z_local"]),
            ("Our selection, calibrated, ZZ only",
             v["ours_calibrated_no_zx"]["z_local_known_mass"]),
            ("… + misidentified-lepton background", v["ours_calibrated_with_zx"]
             ["z_local_known_mass"])]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.9), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    ys = np.arange(len(rows))[::-1]
    for y, (lab, z) in zip(ys, rows):
        ax.barh(y, z, height=0.55, color=C_H if y == ys[-1] else ps.RAMPS[ps.ORANGE][0])
        ax.text(z + 0.06, y, f"{z:.1f}σ", va="center", fontsize=8.5)
    ax.axvline(t["expected_z_median"], color=C_ZZ, ls="--", linewidth=1)
    ax.text(t["expected_z_median"] + 0.05, ys[0] + 0.45, f"expected {t['expected_z_median']:.1f}σ",
            color=C_ZZ, fontsize=7.5)
    ax.set_yticks(ys, [r[0] for r in rows])
    ax.set_xlim(0, 5)
    ax.set_xlabel("Local significance of the excess at 125 GeV")
    ax.grid(axis="y", visible=False)
    ax.set_title("How strong is the evidence? It depends on the bookkeeping", fontsize=10)
    ax = axes[1]
    h = np.array(t["q0_bg_hist"], dtype=float)
    e = np.linspace(0, 16, 65)
    dens = h / h.sum() / np.diff(e)
    ax.stairs(dens, e, color=C_ZZ, fill=True, alpha=0.35, label="background-only pseudo-experiments")
    x = np.linspace(0.05, 16, 400)
    ax.plot(x, 0.5 * stats.chi2.pdf(x, 1), color=ps.INK_SECONDARY, ls="--", linewidth=1,
            label="asymptotic ½χ²(1)")
    q_obs = v["ours_calibrated_with_zx"]["q0_at_known_mass"]
    ax.axvline(q_obs, color=C_H, linewidth=1.5)
    ax.text(q_obs + 0.2, 0.3, f"data\np = {t['p_local_known_mass_toys']:.4f}", fontsize=7.5,
            color=C_H)
    ax.set_yscale("log")
    ax.set_xlim(0, 16)
    ax.set_xlabel("Test statistic q₀ at m_H = 125.13 GeV")
    ax.set_ylabel("Probability density")
    ax.set_title("Calibrating the significance with pseudo-experiments", fontsize=10)
    ax.legend(loc="upper right", fontsize=7)
    ps.caption(fig, "Left: binned fit to the published 36-bin histogram (top two bars), unbinned "
               "fit of our selection (bottom two). Right: the fraction of the background-only "
               "spike at q₀ = 0 is not drawn.", y=-0.03)
    return _save(fig, "07_evidence_steps")


# ---------------------------------------------------------------------------
# 7. The Higgs potential and the couplings
# ---------------------------------------------------------------------------


def fig_potential(R):
    pot = R["electroweak"]["potential"]
    v, lam = pot["v"], pot["lambda_tree"]
    fig = plt.figure(figsize=(11.5, 4.8))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    r = np.linspace(0, 1.36 * v, 80)
    th = np.linspace(0, 2 * np.pi, 120)
    Rr, T = np.meshgrid(r, th)
    V = lam * (Rr**2 - v**2) ** 2 / 4 - lam * v**4 / 4      # V(0) = 0
    X, Y = Rr * np.cos(T), Rr * np.sin(T)
    ax.plot_surface(X, Y, V / 1e8, cmap="Blues_r", linewidth=0, antialiased=True, alpha=0.95,
                    rstride=2, cstride=2)
    ax.plot(v * np.cos(th), v * np.sin(th), np.full_like(th, -lam * v**4 / 4 / 1e8),
            color=C_H, linewidth=2)
    ax.set_xlabel("Re H  [GeV]", fontsize=8)
    ax.set_ylabel("Im H  [GeV]", fontsize=8)
    ax.set_zticklabels([])
    ax.set_title("The Mexican hat: the Higgs potential", fontsize=10)
    ax.view_init(elev=33, azim=-60)
    ax2 = fig.add_subplot(1, 2, 2)
    h = np.linspace(-1.5 * v, 1.5 * v, 600)
    Vh = lam * (h**2 - v**2) ** 2 / 4 - lam * v**4 / 4
    ax2.plot(h, Vh / 1e8, color=C_ZZ)
    ax2.plot([v, -v], [-lam * v**4 / 4 / 1e8] * 2, "o", color=C_H, ms=6)
    ax2.annotate(f"our vacuum: |H| = v/√2,  v = {v:.2f} GeV", (v, -lam * v**4 / 4 / 1e8),
                 xytext=(v * 0.12, -1.38), fontsize=8, color=ps.INK_SECONDARY,
                 arrowprops={"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6})
    ax2.annotate("H = 0: symmetric, all particles massless\n(the early universe above 160 GeV)",
                 (0, 0), xytext=(-1.45 * v, 0.62), fontsize=8, color=ps.INK_SECONDARY,
                 arrowprops={"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6})
    ax2.set_xlabel("Higgs field value v + h  [GeV]")
    ax2.set_ylabel("V  [10⁸ GeV⁴]")
    ax2.set_ylim(-1.45, 1.0)
    ax2.set_title(f"λ = M_h²/2v² = {lam:.4f},  μ = M_h/√2 = {pot['mu']:.1f} GeV", fontsize=10)
    ps.caption(fig, f"V = −μ²|H|² + λ|H|⁴ (tree level), v from the Fermi constant, λ from the Higgs "
               f"mass. The vacuum sits ({pot['depth_fourth_root_GeV']:.0f} GeV)⁴ below the top of "
               f"the hat: {_sci(10 ** pot['log10_depth_over_rho_lambda'])} times the dark-energy density "
               f"({pot['rho_lambda_fourth_root_meV']:.2f} meV)⁴ — the cosmological-constant problem.",
               y=-0.02)
    return _save(fig, "08_higgs_potential")


def fig_couplings(R):
    at = R["electroweak"]["atlas"]
    sm = R["electroweak"]["yukawas"]
    v = R["electroweak"]["potential"]["v"]
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    m = np.geomspace(3e-4, 300, 50)
    ax.plot(m, m / v, color=ps.MUTED, linewidth=1, ls="--", label="coupling ∝ mass (Standard Model)")
    measured = {r["particle"] for r in at}
    for r in sm:
        if r["particle"] in measured:
            continue
        ax.plot(r["mass_GeV"], r["reduced_coupling"], "o", mfc="white", mec=ps.MUTED, ms=6)
        ax.annotate(r["particle"], (r["mass_GeV"], r["reduced_coupling"]), xytext=(5, -10),
                    textcoords="offset points", fontsize=8, color=ps.MUTED)
    for r in at:
        col = C_ZZ if r["particle"] in ("W", "Z") else C_H
        ax.errorbar(r["mass_GeV"], r["reduced"], yerr=[[-r["reduced_down"]], [r["reduced_up"]]],
                    fmt="o", color=col, ms=6)
        off = {"W": (-16, -4), "Z": (2, 9)}.get(r["particle"], (-12, 6))
        ax.annotate(r["particle"], (r["mass_GeV"], r["reduced"]), xytext=off,
                    textcoords="offset points", fontsize=9, color=ps.INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Particle mass  [GeV]")
    ax.set_ylabel("Coupling to the Higgs boson: κ_F m_F / v  or  √κ_V m_V / v")
    handles = [Line2D([], [], color=ps.MUTED, ls="--", label="Standard Model: coupling ∝ mass"),
               Line2D([], [], marker="o", ls="", color=C_H, label="fermions, measured (ATLAS 2022)"),
               Line2D([], [], marker="o", ls="", color=C_ZZ,
                      label="W and Z bosons, measured (ATLAS 2022)"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=ps.MUTED,
                      label="not yet measured (predicted)")]
    ax.legend(handles=handles, loc="upper left", fontsize=7.5)
    ax.set_title("The heavier the particle, the harder it pulls on the Higgs field")
    ps.caption(fig, "Measured coupling modifiers from ATLAS Run 2 (Nature 607, 52; HEPData "
               "ins2104706, Fig. 5). Hollow: not yet measured, drawn where the Standard Model puts "
               "them (PDG 2026 masses); the electron's, 2 × 10⁻⁶, is the weakest.", y=-0.02)
    return _save(fig, "09_couplings_vs_mass")


# ---------------------------------------------------------------------------
# 8. Running couplings, lambda, vacuum stability
# ---------------------------------------------------------------------------


def fig_running(R, curves):
    S = R["stability"]
    head = "smdr" if S["headline"].endswith("smdr") else "buttazzo"
    mu = curves[f"{head}_mu"]
    X = curves[f"{head}_X"]
    lmu = np.log10(mu)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    ax = axes[0]
    series = [(np.sqrt(X[:, 0]), "g₁ hypercharge", C_ZZ, "-", 0.0),
              (np.sqrt(X[:, 1]), "g₂ weak", C_ZZ, "--", 0.022),
              (np.sqrt(X[:, 2]), "g₃ strong", C_ZZ, ":", -0.022),
              (np.sqrt(X[:, 3]), "y_t top quark", C_H, "-", 0.0),
              (X[:, 6], "λ Higgs", C_OTHER, "-", 0.0)]
    for y, lab, c, ls, dy in series:
        ax.plot(lmu, y, color=c, ls=ls, linewidth=1.5)
        ax.text(lmu[-1] + 0.2, y[-1] + dy, lab, fontsize=7.5, va="center", color=ps.INK_SECONDARY)
    ax.axhline(0, color=ps.AXIS, linewidth=0.8)
    ax.axvspan(np.log10(st.M_PL), 20.5, color=ps.EXCLUDED, linewidth=0)
    ax.set_xlim(2, 23.6)
    ax.set_ylim(-0.1, 1.25)
    ax.set_xticks(range(2, 21, 2), [f"10{_sup(k)}" for k in range(2, 21, 2)])
    ax.set_xlabel("Energy scale μ  [GeV]")
    ax.set_ylabel("Coupling (MS-bar)")
    ax.set_title("All couplings, top mass to Planck mass", fontsize=10)
    ax.text(np.log10(st.M_PL) + 0.1, 1.15, "beyond\nPlanck", fontsize=7, color=ps.MUTED)
    ax = axes[1]
    for n, alpha in ((3, 0.18), (1, 0.4)):
        lo = curves[f"var_Mt{n:+d}_lam"]
        hi = curves[f"var_Mt{-n:+d}_lam"]
        mm = np.log10(curves[f"var_Mt{n:+d}_mu"])
        mh_ = np.log10(curves[f"var_Mt{-n:+d}_mu"])
        grid = np.linspace(2.3, 20, 800)
        ax.fill_between(grid, np.interp(grid, mm, lo), np.interp(grid, mh_, hi), color=C_OTHER,
                        alpha=alpha, linewidth=0, label=f"M_t ± {n}σ")
    ax.plot(lmu, X[:, 6], color=C_OTHER, linewidth=1.8, label="λ, 3-loop running")
    big = mu > 1e6
    leff = rge_lambda_eff(X[big])
    ax.plot(lmu[big], leff, color=C_OTHER, ls="--", linewidth=1.1, label="λ_eff (2-loop potential)")
    for n, ls in ((1, ":"), (2, "-.")):
        ax.plot(np.log10(curves[f"loops{n}_mu"]), curves[f"loops{n}_lam"], color=ps.MUTED,
                ls=ls, linewidth=1.0, label=f"{n}-loop running")
    if "smdr_Q" in curves:
        ax.plot(np.log10(curves["smdr_Q"])[::2], curves["smdr_lam"][::2], "o", ms=4.5, mfc="none",
                mec=C_H, mew=1.0, label="SMDR's own running (Martin & Robertson)")
    ax.axhline(0, color=ps.INK, linewidth=0.8)
    ax.axvspan(np.log10(st.M_PL), 20.5, color=ps.EXCLUDED, linewidth=0)
    r = S["results"][S["headline"]]
    ax.annotate(f"λ < 0 above {_sci(r['Lambda_lambda'])} GeV", (np.log10(r["Lambda_lambda"]), 0),
                xytext=(5.8, -0.035), fontsize=8, color=ps.INK_SECONDARY,
                arrowprops={"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6})
    ax.set_xlim(2, 20.5)
    ax.set_ylim(-0.05, 0.14)
    ax.set_xticks(range(2, 21, 2), [f"10{_sup(k)}" for k in range(2, 21, 2)])
    ax.set_xlabel("Energy scale μ  [GeV]")
    ax.set_ylabel("Higgs self-coupling λ")
    ax.set_title("The Higgs self-coupling turns negative", fontsize=10)
    ax.legend(loc="upper right", fontsize=7)
    ps.caption(fig, "PDG 2026 masses (M_h = 125.13, M_t = 172.60 GeV, α_s = 0.1180) matched to "
               "MS-bar couplings at M_t as in SMDR; three-loop RGEs (Buttazzo et al. 2013, "
               "appendix B); the small b and τ Yukawa couplings are included but not drawn. "
               "One-loop running puts the zero two orders of magnitude too low.", y=-0.02)
    return _save(fig, "10_running_couplings")


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _sup(n) -> str:
    return str(n).translate(_SUP)


def _sci(x: float, d: int = 1) -> str:
    e = int(np.floor(np.log10(x)))
    return f"{x / 10**e:.{d}f} × 10{_sup(e)}"


def rge_lambda_eff(Xrows):
    import rge
    return rge.lambda_eff(np.asarray(Xrows).T, loops=2)


def _phase_image(ax, g, tag, *, labels=True):
    mt, mh = g[f"{tag}_mt"], g[f"{tag}_mh"]
    leff = g[f"{tag}_lambda_eff_min"]
    p0 = g[f"{tag}_log10_p0"]
    lmpl = g[f"{tag}_lambda_MPl"]
    # smooth fields for each border: lambda_eff = 0 (stability), P = 1 (instability),
    # lambda(M_Pl) = 1 (the self-coupling grows too strong before the Planck mass)
    lmax = g.get(f"{tag}_lambda_max")
    big = np.where(np.isfinite(lmax), lmax, 50.0) if lmax is not None else \
        np.where(np.isfinite(lmpl), lmpl, 50.0)
    big = np.clip(big, -50.0, 50.0)        # lambda can grow huge (finite) near a Landau pole
    leff_c = np.where(np.isfinite(leff), leff, 1.0)
    ax.set_facecolor(WASH["metastable"])
    ax.contourf(mt, mh, leff_c, levels=[0, 1e9], colors=[WASH["stable"]])
    ax.contourf(mt, mh, np.where(np.isfinite(p0), p0, -1e9), levels=[0, 1e12],
                colors=[WASH["unstable"]])
    ax.contourf(mt, mh, big, levels=[1.0, 1e9], colors=[WASH["nonpert"]])
    ax.contour(mt, mh, np.where(big > 1.0, np.nan, leff_c), levels=[0], colors=[C_ZZ],
               linewidths=1.2)
    ax.contour(mt, mh, np.where(big > 1.0, np.nan, p0), levels=[0], colors=[C_H], linewidths=1.2)


def fig_phase_diagram(R, g):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
    ax = axes[0]
    _phase_image(ax, g, "wide")
    ax.set_xlabel("Top quark pole mass M_t  [GeV]")
    ax.set_ylabel("Higgs boson pole mass M_h  [GeV]")
    ax.set_title("The Standard Model's vacuum: phase diagram", fontsize=10)
    for txt, x, y, rot in (("stable", 100, 90, 0), ("metastable", 153, 82, 62),
                           ("unstable", 200, 70, 0), ("λ non-perturbative", 130, 188, 0)):
        ax.text(x, y, txt, fontsize=9, color=ps.INK_SECONDARY, ha="center", va="center",
                rotation=rot)
    zoom = (166, 180, 118, 132)
    ax.add_patch(plt.Rectangle((zoom[0], zoom[2]), zoom[1] - zoom[0], zoom[3] - zoom[2],
                               fill=False, edgecolor=ps.INK, linewidth=0.8))
    ax.plot(cfg.M_TOP[0], cfg.M_HIGGS[0], "o", color=ps.INK, ms=3)
    ax = axes[1]
    _phase_image(ax, g, "zoom")
    mt, mh = g["zoom_mt"], g["zoom_mh"]
    lI = g["zoom_log10_Lambda_I"]
    cs = ax.contour(mt, mh, lI, levels=[9, 10, 11, 12, 14, 16, 18], colors=[ps.MUTED],
                    linewidths=0.6, linestyles=":")
    ax.clabel(cs, fmt=lambda x: f"10{_sup(int(round(x)))} GeV", fontsize=6.5)
    t0, h0 = cfg.M_TOP[0], cfg.M_HIGGS[0]
    st_, sh = cfg.M_TOP[1], cfg.M_HIGGS[1]
    for k in (1, 2, 3):
        ax.add_patch(Ellipse((t0, h0), 2 * k * st_, 2 * k * sh, fill=False, edgecolor=ps.INK,
                             linewidth=0.9 if k == 1 else 0.6))
    ax.plot(t0, h0, "o", color=ps.INK, ms=4, label="PDG 2026 (1, 2, 3σ)")
    ax.plot(172.57, 125.20, "s", mfc="white", mec=ps.INK, ms=4, label="PDG 2024")
    S = R["stability"]
    b = S["results"][S["headline"]]["stability_bound_mt"]
    ax.annotate(f"absolutely stable only\nif M_t < {b:.2f} GeV", (b, h0), xytext=(166.6, 129.5),
                fontsize=8, color=C_ZZ, arrowprops={"arrowstyle": "-", "color": C_ZZ,
                                                    "linewidth": 0.6})
    for txt, x, y in (("stable", 168.0, 121.0), ("metastable", 176.0, 120.0)):
        ax.text(x, y, txt, fontsize=9, color=ps.INK_SECONDARY, ha="center")
    ax.set_xlim(mt.min(), mt.max())
    ax.set_ylim(mh.min(), mh.max())
    ax.set_xlabel("Top quark pole mass M_t  [GeV]")
    ax.set_ylabel("Higgs boson pole mass M_h  [GeV]")
    gap = S["stability_gap"]
    ax.set_title(f"We live {gap['gap_GeV']:.2f} GeV of top mass from stability "
                 f"({gap['n_sigma_exp_th']:.1f}σ)", fontsize=10)
    ax.legend(loc="lower right", fontsize=7)
    ps.caption(fig, "Stable: the Higgs potential (two-loop, λ_eff) stays positive up to the Planck "
               "mass. Unstable: our vacuum would probably have decayed inside our past light cone "
               "(leading-order bounce, Planck 2018 cosmology). Grey: λ exceeds 1 before the "
               "Planck mass. Dotted: the scale where the potential turns negative. "
               "α_s = 0.1180.", y=-0.02)
    return _save(fig, "11_vacuum_phase_diagram")


# ---------------------------------------------------------------------------
# 9. The electroweak calendar
# ---------------------------------------------------------------------------


def fig_calendar(R):
    lat = ew.lattice_condensate()
    T = lat["T"]
    t = ew.p2_time_of(T)["t_s"]
    cal = R["electroweak"]["calendar"]
    fig, ax = plt.subplots(figsize=(10, 5.0))
    phi = lat["phi2_over_T"]
    ax.fill_between(t * 1e12, phi - lat["err"], phi + lat["err"], color=C_H, alpha=0.3,
                    linewidth=0)
    ax.plot(t * 1e12, phi, color=C_H, linewidth=2,
            label="Higgs condensate ⟨φ†φ⟩/T (lattice, D'Onofrio & Rummukainen 2016)")
    ax.axhline(0, color=ps.AXIS, linewidth=0.8)
    c = cal["crossover"]
    ax.axvspan(cal["crossover_hi"]["t_s"] * 1e12, cal["crossover_lo"]["t_s"] * 1e12,
               color=ps.RAMPS[ps.ORANGE][0], linewidth=0)
    yc = float(np.interp(c["t_s"], t[::-1], phi[::-1]))
    ax.annotate(f"crossover at {c['T_GeV']:.1f} ± 1.5 GeV,\n{c['t_s'] * 1e12:.1f} ps after the Big Bang",
                (c["t_s"] * 1e12, yc), xytext=(10.4, 0.12), fontsize=8,
                ha="left", color=ps.INK_SECONDARY,
                arrowprops={"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6})
    s = cal["sphaleron_freeze_out"]
    ax.axvline(s["t_s"] * 1e12, color=C_OTHER, linewidth=1, ls="--")
    ax.text(s["t_s"] * 1e12 * 1.02, 0.42, f"sphalerons freeze:\nbaryon number fixed\n"
            f"({s['T_GeV']:.1f} GeV, {s['t_s'] * 1e12:.1f} ps)", fontsize=7.5, color=C_OTHER)
    ax.set_xlim(8.0, 14.6)
    ax.set_ylim(-0.03, 0.66)
    ax.set_xlabel("Time after the Big Bang  [picoseconds]")
    ax.set_ylabel("Higgs condensate ⟨φ†φ⟩ / T")
    ax.set_title("The Higgs field switches on: a smooth crossover, not a phase transition")
    ax.legend(loc="upper left", fontsize=7.5)
    to_T = lambda x: np.interp(np.asarray(x, float), (t * 1e12)[::-1], T[::-1])
    to_t = lambda y: np.interp(np.asarray(y, float), T, t * 1e12)
    top = ax.secondary_xaxis("top", functions=(to_T, to_t))
    Tt = [170, 165, 160, 155, 150, 145, 140]
    top.set_xticks(Tt, [str(v) for v in Tt])
    top.set_xlabel("Temperature  [GeV]")
    later = []
    for key, lab in (("top_half", "top"), ("higgs_half", "Higgs"), ("Z_half", "Z"),
                     ("W_half", "W")):
        if key in cal:
            later.append(f"{lab} {cal[key]['t_s'] * 1e12:.0f} ps")
    ps.caption(fig, "Continuum lattice results for the Higgs condensate in the 3D effective theory "
               "(the authors' data file, band = their error); in the broken phase v(T)² ≈ "
               "2⟨φ†φ⟩. Time from P2's equation of state. Later, half of each heavy species has "
               "annihilated by: " + ", ".join(later) + ".", y=-0.02)
    return _save(fig, "12_electroweak_calendar")


# ---------------------------------------------------------------------------


def make_all(R: dict) -> list:
    out = []
    d = dimuon.load()
    out += fig_dimuon_spectrum(d, R)
    out += fig_resonance_masses(d, R)
    out += fig_afb(d, R)
    out += fig_replica(R)
    mc = fl.mc_candidates()
    out += fig_four_lepton(R, mc, None)
    out += fig_higgs_fit(R)
    out += fig_evidence(R)
    out += fig_potential(R)
    out += fig_couplings(R)
    curves = dict(np.load(st.CURVES_CACHE))
    out += fig_running(R, curves)
    out += fig_phase_diagram(R, st.grids())
    out += fig_calendar(R)
    return out


if __name__ == "__main__":
    import sys
    R = json.loads((cfg.CACHE_DIR / "full_results.json").read_text())
    if len(sys.argv) > 1:          # draw only the named figure functions
        d = dimuon.load() if any("dimuon" in a or "resonance" in a or "afb" in a
                                 for a in sys.argv[1:]) else None
        for name in sys.argv[1:]:
            f = globals()[name]
            args = {"fig_dimuon_spectrum": (d, R), "fig_resonance_masses": (d, R),
                    "fig_afb": (d, R), "fig_four_lepton": (R, fl.mc_candidates(), None),
                    "fig_running": (R, dict(np.load(st.CURVES_CACHE))),
                    "fig_phase_diagram": (R, st.grids())}.get(name, (R,))
            print(f(*args))
        sys.exit(0)
    for f in make_all(R):
        print(f)
