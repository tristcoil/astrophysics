"""Every figure of the project, from the analysis results and artefacts."""

from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from astropy.time import Time

import gaia_orbits
import literature as lit
import masses
import plotstyle as ps
from config import FIGURE_DIR, J1820_HOMAN_TRANSITION, NS_MAX_CAUSAL, NS_MAX_KALOGERA_BAYM

ps.use_style()


def _label(row_or_adopted) -> str:
    name = row_or_adopted["name"]
    common = row_or_adopted.get("common_name") or ""
    short = common.split(" (")[0]
    return f"{name} ({short})" if short else name


def _wrap(text: str, width: int) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width)


def _log_mass_axis(ax, lo=1.0, hi=40.0, axis="x"):
    ticks = [t for t in (1, 1.5, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 200, 300) if lo <= t <= hi]
    if axis == "x":
        ax.set_xscale("log")
        ax.set_xlim(lo, hi)
        ax.xaxis.set_major_locator(mticker.FixedLocator(ticks))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_locator(mticker.NullLocator())
    else:
        ax.set_yscale("log")
        ax.set_ylim(lo, hi)
        ax.yaxis.set_major_locator(mticker.FixedLocator(ticks))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_locator(mticker.NullLocator())


def _gap_band(ax, axis="x", lo=2.5, hi=5.0, label=True):
    if axis == "x":
        ax.axvspan(lo, hi, color=ps.GAP_FILL, zorder=0, lw=0)
    else:
        ax.axhspan(lo, hi, color=ps.GAP_FILL, zorder=0, lw=0)


def _violin(ax, y, samples, color, *, width=0.8, lo=None, hi=None, alpha=0.85, logx=True):
    """Horizontal half-violin (density of log mass) centred on row y."""
    s = np.asarray(samples)
    s = s[np.isfinite(s) & (s > 0)]
    x = np.log10(s) if logx else s
    a, b = np.quantile(x, [0.002, 0.995])
    if lo is not None:
        a = max(a, np.log10(lo) if logx else lo)
    if hi is not None:
        b = min(b, np.log10(hi) if logx else hi)
    grid = np.linspace(a, b, 300)
    hist, edges = np.histogram(x, bins=120, range=(a, b), density=True)
    centres = 0.5 * (edges[1:] + edges[:-1])
    dens = np.interp(grid, centres, np.convolve(hist, np.ones(5) / 5, mode="same"))
    dens = dens / dens.max() * width / 2
    gx = 10**grid if logx else grid
    ax.fill_between(gx, y - dens, y + dens, color=color, alpha=alpha, lw=0, zorder=3)
    q16, q50, q84 = np.quantile(s, [0.16, 0.5, 0.84])
    ax.plot([q16, q84], [y, y], color=ps.SURFACE, lw=1.2, zorder=4, solid_capstyle="butt")
    ax.plot(q50, y, "o", ms=3.2, color=ps.SURFACE, zorder=5)
    return q50


# ---------------------------------------------------------------------------
# 1. Mass functions
# ---------------------------------------------------------------------------

def fig_mass_functions(res, art):
    w, wart = res["weighing"], art["weighing"]
    rows = sorted(w["systems"], key=lambda r: r["fm_ours"]["median"])
    adopted = {a["name"]: a for a in wart["adopted"]}
    pers = w["persistent"]
    n_t, n_p = len(rows), len(pers)
    fig, ax = plt.subplots(figsize=(8.6, 0.27 * (n_t + n_p) + 2.0))
    y_t = np.arange(n_t)[::-1] + n_p + 1.2
    y_p = np.arange(n_p)[::-1]
    ax.axvspan(NS_MAX_CAUSAL, 20, color=ps.GAP_FILL, lw=0, zorder=0)
    import matplotlib.transforms as mtransforms
    top = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for xv, lab, ha in ((2.08, "2.08", "center"), (NS_MAX_KALOGERA_BAYM, "2.9 ", "right"),
                        (NS_MAX_CAUSAL, " 3.2", "left")):
        ax.axvline(xv, color=ps.MUTED if xv == 2.08 else ps.AXIS, lw=0.8, zorder=1)
        ax.text(xv, 1.005, lab, transform=top, fontsize=7, color=ps.MUTED, ha=ha, va="bottom")
    for y, r in zip(y_t, rows):
        f = r["fm_ours"]
        ax.errorbar(f["median"], y, xerr=[[f["minus"]], [f["plus"]]], fmt="o", ms=4.5,
                    color=ps.C_OURS, ecolor=ps.C_OURS, elinewidth=1.2, capsize=0, zorder=4)
        pub = adopted[r["name"]]["fm_published"]
        if pub is not None:
            ax.errorbar(pub["value"], y - 0.28, xerr=[[pub.get("minus", 0)], [pub.get("plus", 0)]],
                        fmt="D", ms=3.6, mfc="none", mec=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED,
                        elinewidth=0.9, capsize=0, zorder=3)
    for y, p in zip(y_p, pers):
        ax.plot(p["fm_from_P_K"], y, "s", ms=4.5, color=ps.INK_SECONDARY, zorder=4)
    labels = [_label(r) for r in rows] + [f"{p['name']} ({p['host']})" for p in pers]
    ax.set_yticks(list(y_t) + list(y_p))
    ax.set_yticklabels(labels, fontsize=7.6)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 12)
    ax.xaxis.set_major_locator(mticker.FixedLocator([0.1, 0.2, 0.5, 1, 2, 3, 5, 10]))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_ylim(-0.8, y_t.max() + 1.3)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("mass function f(M) = P K³ / 2πG   [M☉]  — a strict lower limit on the unseen mass")
    ax.text(NS_MAX_CAUSAL * 1.16, y_t.max() + 0.2, "f(M) > 3.2 M☉:\nno neutron star\ncan be this heavy\n"
            "(Rhoades & Ruffini\n1974): a black hole\nwhatever the\ninclination",
            fontsize=7.2, color=ps.INK_SECONDARY, va="top")
    ax.text(0.105, n_p + 0.35, "Persistent black holes with massive companions: f(M) is small\n"
            "because the visible star is the heavier one", fontsize=7.2, color=ps.INK_SECONDARY, va="center")
    j = next((k for k, r in enumerate(rows) if r["name"] == "GRO J1655-40"), None)
    if j is not None:
        ps.annotate(ax, "BlackCAT pairs a 2008 K₂ (226 km/s) with f(M)\nfrom a 1999 K₂ (215.5 km/s): 4σ apart",
                    (adopted["GRO J1655-40"]["fm_published"]["value"] / 1.04, y_t[j] - 0.28),
                    (0.42, y_t[j] - 0.28), fontsize=7, relpos=(1, 0.5), ha="left", va="center")
    h1 = ax.errorbar([], [], xerr=[], fmt="o", color=ps.C_OURS, label="ours, from BlackCAT's P and K₂")
    h2 = ax.errorbar([], [], xerr=[], fmt="D", mfc="none", mec=ps.C_PUBLISHED, color=ps.C_PUBLISHED,
                     label="f(M) as listed in BlackCAT")
    h3, = ax.plot([], [], "s", color=ps.INK_SECONDARY, label="from P and K (literature)")
    ax.legend(handles=[h1, h2, h3], loc="lower right", fontsize=7.4)
    ax.set_title(f"The mass functions of {n_t} black-hole transients: {len(w['fm_above_3_2'])} exceed "
                 "any possible neutron star", pad=16)
    ps.caption(fig, "Data: BlackCAT online (Corral-Santana et al. 2016 and updates, Sept 2025), each P and K₂ "
                    "as listed with its paper. Persistent systems: Miller-Jones et al. 2021; Orosz et al. 2009, 2014, 2007.\n"
                    "Vertical lines: the heaviest neutron star weighed (PSR J0740+6620, 2.08 M☉); the maximum neutron-star "
                    "mass for an equation of state matched to nuclear data (2.9 M☉, Kalogera & Baym 1996) and from "
                    "causality alone (3.2 M☉, Rhoades & Ruffini 1974).", y=0.0)
    return ps.save(fig, FIGURE_DIR / "01_mass_functions")


# ---------------------------------------------------------------------------
# 2. From a mass function to a mass
# ---------------------------------------------------------------------------

EXAMPLES = [("GS 2023+338", "V404 Cyg", ps.BLUE), ("3A 0620-003", "A0620-00", ps.VIOLET),
            ("MAXI J1820+070", "MAXI J1820+070", ps.GREEN), ("GRO J0422+32", "GRO J0422+32", ps.AQUA)]


def fig_inclination(res, art, rng):
    draws = art["weighing"]["draws"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.05, 1]})
    incl = np.linspace(8, 90, 400)
    for ax in (ax1,):
        ax.axhspan(2.08, NS_MAX_CAUSAL, color=ps.GAP_FILL, lw=0, zorder=0)
        ax.text(9.5, 2.55, "neutron stars stop here (2.1–3.2 M☉)", fontsize=7, color=ps.MUTED, va="center")
    for name, label, col in EXAMPLES:
        d = draws[name]
        fm, q = np.median(d["fm"]), np.median(d["q"])
        m = masses.m1_from_fm(fm, q, incl)
        ax1.plot(incl, m, color=col, lw=1.0, alpha=0.55)
        lo, hi = np.quantile(d["incl"], [0.16, 0.84])
        sel = (incl >= lo) & (incl <= hi)
        ax1.plot(incl[sel], m[sel], color=col, lw=3.2, solid_capstyle="round")
        ax1.plot(90, masses.m1_from_fm(fm, q, 90.0), "o", ms=5, color=col, mec=ps.SURFACE, mew=1.5, zorder=5)
        offset = {"V404 Cyg": 1.22, "MAXI J1820+070": 0.82}.get(label, 1.0)
        ax1.text(91.5, masses.m1_from_fm(fm, q, 90.0) * offset, f"{label}\nf = {fm:.2f}", fontsize=7.2,
                 color=ps.INK_SECONDARY, va="center")
        bins = np.geomspace(1, 80, 90)
        h, e = np.histogram(d["m1"], bins=bins, density=True)
        c = np.sqrt(e[1:] * e[:-1])
        h = h * c  # density per log interval
        ax2.plot(c, h / h.max(), color=col, lw=1.8, label=f"{label}: {np.median(d['m1']):.1f} M☉")
        ax2.fill_between(c, 0, h / h.max(), color=col, alpha=0.10, lw=0)
    # V404 Cyg if nothing were known about the inclination (only: no eclipses)
    d = draws["GS 2023+338"]
    q = d["q"]
    imax = masses.max_inclination_no_eclipse(q)
    i_iso = np.degrees(np.arccos(rng.uniform(np.cos(np.radians(imax)), 1.0)))
    m_iso = masses.m1_from_fm(d["fm"], q, i_iso)
    bins = np.geomspace(1, 80, 90)
    h, e = np.histogram(m_iso, bins=bins, density=True)
    c = np.sqrt(e[1:] * e[:-1])
    h = h * c
    ax2.plot(c, h / h.max() * 0.45, color=ps.MUTED, lw=1.2,
             label=f"V404 Cyg if i were unknown: {np.median(m_iso):.0f} M☉ median")
    _log_mass_axis(ax1, 1, 100, axis="y")
    ax1.set_xlim(8, 104)
    ax1.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax1.set_xlabel("orbital inclination i  [deg]   (90° = edge-on)")
    ax1.set_ylabel("black hole mass  M₁ = f (1 + q)² / sin³ i   [M☉]")
    ax1.set_title("(a) One orbit, many possible masses")
    ax1.text(0.03, 0.04, "thick: the inclination range measured from the\ncompanion's ellipsoidal light curve "
             "(68%)\ndot: the minimum mass, edge-on", transform=ax1.transAxes, fontsize=7.2,
             color=ps.INK_SECONDARY)
    _log_mass_axis(ax2, 2, 60)
    _gap_band(ax2)
    ax2.set_ylim(0, 1.08)
    ax2.set_yticks([])
    ax2.set_xlabel("black hole mass  [M☉]")
    ax2.set_title("(b) The masses, marginalised over inclination")
    ax2.legend(loc="upper right", fontsize=7.2)
    ax2.text(3.5, 1.02, "2.5–5 M☉", fontsize=7, color=ps.MUTED, ha="center")
    ps.caption(fig, "Monte Carlo over the measured P, K₂, mass ratio q and inclination (BlackCAT online). "
                    "Isotropic inclinations (cos i uniform) where only a range is known. "
                    "Grey: V404 Cyg with only the no-eclipse limit (Eggleton 1983 Roche lobe).", y=-0.03)
    return ps.save(fig, FIGURE_DIR / "02_inclination")


# ---------------------------------------------------------------------------
# 3. The mass gap, system by system
# ---------------------------------------------------------------------------

def fig_masses(res, art):
    w = res["weighing"]
    draws = art["weighing"]["draws"]
    adopted = {a["name"]: a for a in art["weighing"]["adopted"]}
    g = art["graveyard"]
    ns = g["ns"]
    items = []
    for r in w["systems"]:
        if r["in_population"]:
            items.append((_label(r), draws[r["name"]]["m1"], ps.C_BH_LIGHT, r["name"]))
    for p in w["persistent"]:
        items.append((f"{p['name']} ({p['host']})", draws[p["name"]]["m1"], ps.C_BH_LIGHT, p["name"]))
    for gb in w["gaia"]:
        items.append((f"{gb['name']} (dormant)", draws[gb["name"]]["m1"], ps.C_BH_GAIA, gb["name"]))
    items.sort(key=lambda t: np.median(t[1]))
    n = len(items)
    fig = plt.figure(figsize=(9.2, 0.25 * n + 4.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.25 * n / 1.6], hspace=0.08)
    ax0 = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharex=ax0)
    # (a) neutron stars and black holes, one point each
    nsx = lit.NS_XRB
    ax0.errorbar(ns["mass"], np.full(len(ns), 2.0) + np.linspace(-0.25, 0.25, len(ns)),
                 xerr=[ns["minus"], ns["plus"]], fmt="o", ms=3, color=ps.C_NS_LIGHT, ecolor=ps.C_NS_LIGHT,
                 elinewidth=0.6, alpha=0.9, label=f"neutron stars, pulsar timing ({len(ns)})")
    ax0.errorbar([v[1] for v in nsx], np.full(len(nsx), 1.0) + np.linspace(-0.2, 0.2, len(nsx)),
                 xerr=[[v[3] for v in nsx], [v[2] for v in nsx]], fmt="s", ms=3.2, mfc="none",
                 mec=ps.C_NS_LIGHT, ecolor=ps.C_NS_LIGHT, elinewidth=0.6,
                 label=f"neutron stars weighed in X-ray binaries ({len(nsx)})")
    med = np.array([np.median(t[1]) for t in items])
    cols = [t[2] for t in items]
    for k, (m, c) in enumerate(zip(med, cols)):
        ax0.plot(m, 0.0 + (k % 5 - 2) * 0.06, "o" if c == ps.C_BH_LIGHT else "^", ms=4.5, color=c,
                 mec=ps.SURFACE, mew=0.8)
    ax0.plot([], [], "o", color=ps.C_BH_LIGHT, label="black holes in X-ray binaries (median)")
    ax0.plot([], [], "^", color=ps.C_BH_GAIA, label="dormant black holes, Gaia (median)")
    ax0.set_yticks([0, 1, 2])
    ax0.set_yticklabels(["black holes", "NS in X-ray\nbinaries", "NS, pulsar\ntiming"], fontsize=7.6)
    ax0.set_ylim(-0.6, 2.6)
    ax0.grid(axis="y", visible=False)
    _gap_band(ax0)
    ax0.legend(loc="upper right", fontsize=7, ncol=2)
    ax0.set_title("(a) Every compact object weighed by light: neutron stars stop near 2 M☉, "
                  "black holes cluster at 7–10 M☉", fontsize=10.5)
    plt.setp(ax0.get_xticklabels(), visible=False)
    # (b) every black hole's mass distribution
    _gap_band(ax)
    for k, (lab, s, c, key) in enumerate(items):
        _violin(ax, k, s, c, width=0.82, lo=1.0, hi=80)
        a = adopted.get(key)
        if a and a["m1_published"] is not None:
            p = a["m1_published"]
            if p["kind"] == "range":
                ax.plot([p["low"], p["high"]], [k + 0.36, k + 0.36], color=ps.C_PUBLISHED, lw=1.0)
            elif "value" in p and p["flag"] == "":
                ax.errorbar(p["value"], k + 0.36, xerr=[[p.get("minus", 0)], [p.get("plus", 0)]],
                            fmt="D", ms=3, mfc="none", mec=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED, elinewidth=0.9)
            elif "value" in p:
                ax.plot(p["value"], k + 0.36, marker=">" if p["flag"] in (">", "≥") else "<", ms=4,
                        color=ps.C_PUBLISHED)
    ax.set_yticks(range(n))
    ax.set_yticklabels([t[0] for t in items], fontsize=7.4)
    ax.set_ylim(-0.7, n - 0.3)
    ax.grid(axis="y", visible=False)
    _log_mass_axis(ax, 1, 40)
    ax.set_xlabel("mass  [M☉]")
    ax.set_title("(b) Each black hole's mass, marginalised over inclination (ours), with BlackCAT's value",
                 fontsize=10)
    h1, = ax.plot([], [], "s", color=ps.C_BH_LIGHT, ms=6, label="ours: density, 68% bar, median")
    h2 = ax.errorbar([], [], xerr=[], fmt="D", mfc="none", mec=ps.C_PUBLISHED, color=ps.C_PUBLISHED,
                     label="BlackCAT (value, range or limit)")
    ax.legend(handles=[h1, h2], loc="lower right", fontsize=7.2)
    ps.caption(fig, "Neutron stars: P. Freire's table of precise pulsar-timing masses (Oct 2026); Özel & Freire 2016 "
                    "(X-ray binaries). Black holes: our Monte Carlo from BlackCAT's P, K₂, q, i; persistent and Gaia "
                    "black holes as published (Miller-Jones et al. 2021; Orosz et al. 2009, 2014; Ramachandran et al. "
                    "2022; El-Badry et al. 2023a,b; Panuzzo et al. 2024).", y=0.04)
    return ps.save(fig, FIGURE_DIR / "03_masses")


# ---------------------------------------------------------------------------
# 4. Is there a gap? Population fits
# ---------------------------------------------------------------------------

SAMPLE_LABELS = {
    "farr2011": "Farr et al. 2011 inputs\n(15 transients)",
    "kreidberg2012": "Kreidberg et al. 2012 inputs\n(16, revised inclinations)",
    "blackcat_lmxb": "BlackCAT 2025, ours\n(transients)",
    "catnoir_lmxb": "CatNoir 2026 inputs\n(same systems)",
    "blackcat_combined": "BlackCAT 2025 + persistent\n+ Gaia (combined)",
}
MODEL_LABELS = {"power_law": "power law", "exponential": "decaying exponential", "gaussian": "Gaussian"}


def fig_population(res, art):
    fits = res["mass_gap"]["fits"]
    m1 = art["mass_gap"]["m1pct"]
    order = ["farr2011", "kreidberg2012", "blackcat_lmxb", "catnoir_lmxb", "blackcat_combined"]
    models = ["power_law", "exponential", "gaussian"]
    fig = plt.figure(figsize=(12, 7.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.42, wspace=0.08)
    farr = lit.CHECKS["farr2011_m1pct"]
    rocha = lit.CHECKS["rocha2026_m1pct_90"]
    for j, model in enumerate(models):
        ax = fig.add_subplot(gs[0, j])
        ax.axvspan(NS_MAX_KALOGERA_BAYM, NS_MAX_CAUSAL, color=ps.GAP_FILL, lw=0, zorder=0)
        for k, label in enumerate(order[::-1]):
            s = m1[(label, model)]
            s = s[np.isfinite(s)]
            neg = float(np.mean(s <= 0))
            if neg > 0.95:  # a Gaussian forced to span 6-33 M_sun puts its 1% point below zero
                ax.text(0.2, k, f"M₁% < 0 in {100 * neg:.0f}% of the posterior: one Gaussian\n"
                        "cannot fit black holes from 6 to 33 M☉", fontsize=6.8, color=ps.INK_SECONDARY,
                        va="center")
            else:
                _violin(ax, k, np.clip(s[s > 0], 0.05, 12), ps.C_OURS, width=0.75, lo=0.05, hi=12,
                        logx=False)
                q10 = np.quantile(s, 0.1)
                ax.plot([max(q10, 0.0)] * 2, [k - 0.3, k + 0.3], color=ps.INK, lw=1.2)
                if neg > 0.01:
                    ax.text(0.15, k - 0.38, f"{100 * neg:.0f}% below 0", fontsize=6.4, color=ps.MUTED)
            lnb = fits[label][model]["ln_bayes_vs_best"]
            ax.text(11.8, k + 0.27, f"ln B = {lnb:+.1f}" if lnb < -0.05 else "preferred",
                    fontsize=6.6, color=ps.MUTED, ha="right", va="center")
            if label == "farr2011":
                for qv in farr[model]:
                    ax.plot(qv, k - 0.36, "|", ms=7, mew=1.6, color=ps.C_PUBLISHED)
            if label == "blackcat_combined" and model == "power_law":
                ax.plot(rocha["combined_power_law"], k - 0.36, "|", ms=7, mew=1.6, color=ps.C_PUBLISHED)
            if label == "catnoir_lmxb" and model == "gaussian":
                ax.plot(rocha["lmxb_gaussian"], k - 0.36, "|", ms=7, mew=1.6, color=ps.C_PUBLISHED)
        ax.set_xlim(0, 12)
        ax.set_ylim(-0.6, len(order) - 0.4)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([SAMPLE_LABELS[l] for l in order[::-1]] if j == 0 else [], fontsize=7.4)
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("lightest black hole, M₁%  [M☉]")
        ax.set_title(f"({'abc'[j]}) {MODEL_LABELS[model]}", fontsize=10)
        if j == 0:
            ax.text(3.05, len(order) - 0.45, "max NS\n2.9–3.2", fontsize=6.6, color=ps.MUTED, ha="center",
                    va="top")
    h1, = plt.gca().plot([], [], "s", color=ps.C_OURS, ms=6, label="ours: posterior of M₁% (black bar: 90% lower bound)")
    h2, = plt.gca().plot([], [], "|", color=ps.C_PUBLISHED, ms=8, mew=1.6,
                         label="published: Farr et al. 2011 10/50/90%; Rocha et al. 2026 90% bound")
    fig.legend(handles=[h1, h2], loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=2, fontsize=7.6)
    fig.suptitle("Is there a mass gap? The 1% quantile of the black hole mass distribution, five samples × three models",
                 x=0.07, ha="left", y=1.0)
    # (d) posterior predictive mass distributions for our transients
    ax = fig.add_subplot(gs[1, :2])
    grid = np.geomspace(1, 40, 600)
    styles = {"power_law": "-", "exponential": "--", "gaussian": ":"}
    _gap_band(ax)
    for model in models:
        chain = art["mass_gap"]["chains"][("blackcat_lmxb", model)]
        mean, lo, hi = masses.posterior_predictive(chain, model, grid)
        ax.plot(grid, mean, ls=styles[model], color=ps.INK, lw=1.4,
                label=f"{MODEL_LABELS[model]}")
    smp = art["mass_gap"]["samples"]["blackcat_lmxb"]
    meds = [np.median(v) for v in smp.values()]
    ax.plot(meds, np.full(len(meds), -0.012), "|", ms=10, mew=1.4, color=ps.C_BH_LIGHT,
            label="the transients' median masses")
    _log_mass_axis(ax, 1, 40)
    ax.set_ylim(-0.03, None)
    ax.set_ylabel("p(M)  [per M☉]")
    ax.set_xlabel("black hole mass  [M☉]")
    ax.set_title("(d) The inferred mass distribution of black holes in transients (BlackCAT 2025), "
                 "averaged over each model's posterior", fontsize=10)
    ax.legend(loc="upper right", fontsize=7.4)
    # (e) table of numbers
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    short = {"farr2011": "Farr 2011 inputs", "kreidberg2012": "Kreidberg 2012 inputs",
             "blackcat_lmxb": "BlackCAT 2025 (ours)", "catnoir_lmxb": "CatNoir 2026 inputs",
             "blackcat_combined": "BlackCAT + persistent + Gaia"}
    lines = ["M₁% > X at 90% confidence", "", f"{'':28s} power  expon  Gauss"]
    for label in order:
        vals = "  ".join(f"{fits[label][m]['m1pct_q10_q50_q90'][0]:5.1f}" for m in models)
        lines.append(f"{short[label]:28s} {vals}")
    lines += ["", "published:", f"{'Farr et al. 2011':28s}   4.3    5.3    2.9",
              f"{'Rocha et al. 2026, LMXBs':28s}     –      –    2.56",
              f"{'Rocha et al. 2026, combined':28s}  4.06     –      –"]
    ax.text(0.0, 1.0, "\n".join(lines), fontsize=7.2, family="DejaVu Sans Mono", va="top",
            color=ps.INK_SECONDARY, transform=ax.transAxes)
    ps.caption(fig, "Hierarchical model of Farr et al. (2011): each system's Monte Carlo mass samples stand in for its "
                    "measurement; uniform priors; nested sampling (dynesty) with the evidence for each model (ln B relative "
                    "to the best). No selection effects modelled. M₁%: the 1% quantile of the fitted mass distribution.",
               y=0.02)
    return ps.save(fig, FIGURE_DIR / "04_mass_gap_population")


# ---------------------------------------------------------------------------
# 5. Inclination, the dominant systematic
# ---------------------------------------------------------------------------

def fig_inclination_bias(res, art):
    inc = res["inclination"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), gridspec_kw={"width_ratios": [1.1, 0.8, 1.25],
                                                                      "wspace": 0.42})
    ax = axes[0]
    a = inc["a0620"]
    i = np.linspace(25, 90, 300)
    d = art["weighing"]["draws"]["3A 0620-003"]
    band = np.array([np.quantile(masses.m1_from_fm(d["fm"], d["q"], x), [0.16, 0.84]) for x in i])
    ax.fill_between(i, band[:, 0], band[:, 1], color=ps.CONTEXT, lw=0)
    ax.plot(i, masses.m1_from_fm(a["fm"], a["q"], i), color=ps.INK, lw=1.4)
    _gap_band(ax, axis="y", lo=2.5, hi=5.0)
    for k, (who, v) in enumerate(sorted(a["history"].items(), key=lambda kv: kv[1]["i"])):
        ax.plot(v["i"], v["m_quoted"], "D", ms=5, mfc="none", mec=ps.C_PUBLISHED, mew=1.3, zorder=5)
        ax.plot(v["i"], v["m_ours"], "o", ms=4.5, color=ps.C_OURS, zorder=6)
        label = f"{who}\ni = {v['i']:.0f}°: {v['m_quoted']:g} M☉"
        if v["i"] > 60:
            ps.annotate(ax, label, (v["i"], v["m_quoted"]), (91, v["m_quoted"] * 1.4), fontsize=6.9,
                        relpos=(0.4, 0), ha="right", va="bottom")
        else:
            ax.text(v["i"] + 1.5, v["m_quoted"] * 1.12, label, fontsize=6.9, color=ps.INK_SECONDARY)
    _log_mass_axis(ax, 2, 30, axis="y")
    ax.set_xlim(25, 92)
    ax.set_xlabel("inclination adopted  [deg]")
    ax.set_ylabel("mass of A0620-00's black hole  [M☉]")
    ax.set_title("(a) A0620-00: one orbit,\nthree published masses", fontsize=10)
    h1, = ax.plot([], [], "D", mfc="none", mec=ps.C_PUBLISHED, label="as published")
    h2, = ax.plot([], [], "o", color=ps.C_OURS, label="ours: today's f(M) and q at that i")
    ax.legend(handles=[h1, h2], loc="lower left", fontsize=7)
    ax = axes[1]
    dd = np.linspace(0, 0.6, 200)
    ax.plot(dd * 100, masses.veiling_mass_bias(dd), color=ps.INK, lw=1.6)
    for dv in (0.2, 0.4):
        ax.plot(dv * 100, masses.veiling_mass_bias(dv), "o", ms=4.5, color=ps.INK)
        ax.text(dv * 100 + 2, masses.veiling_mass_bias(dv), f"{masses.veiling_mass_bias(dv):.2f}×",
                fontsize=7.4, color=ps.INK_SECONDARY, va="top")
    ax.set_xlabel("share of the optical light from the disc,\nignored in the fit  [%]")
    ax.set_ylabel("mass overestimated by a factor")
    ax.set_title("(b) Hidden disc light\ninflates masses", fontsize=10)
    ax.text(2, 2.95, "ellipsoidal amplitude ∝ sin² i, diluted\nby (1 − d); mass ∝ 1/sin³ i\n→ factor (1 − d)^−3/2",
            fontsize=7, color=ps.INK_SECONDARY, va="top")
    ax = axes[2]
    rows = [r for r in inc["distance_to_gap"] if r["disc_fraction_needed"] is not None]
    rows.sort(key=lambda r: r["disc_fraction_needed"])
    robust = [r for r in inc["distance_to_gap"] if r["disc_fraction_needed"] is None]
    y = np.arange(len(rows))
    vals = [100 * r["disc_fraction_needed"] for r in rows]
    ax.barh(y, vals, height=0.62, color=ps.C_BH_LIGHT)
    for k, (r, v) in enumerate(zip(rows, vals)):
        txt = (f"already below 5 M☉ ({r['m1_median']:.1f} M☉)" if v <= 0 else
               f"{r['m1_median']:.1f} M☉ at i = {r['i_median']:.0f}° → 5 M☉ at {r['i_for_5msun']:.0f}°")
        ax.text(max(v, 0) + 1, k, txt, fontsize=6.4, color=ps.INK_SECONDARY, va="center")
    ax.text(0.98, 0.03, f"+ {len(robust)} more that no inclination brings\nbelow 5 M☉: f(M)(1 + q)² > 5 M☉",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8, color=ps.INK_SECONDARY)
    ax.set_yticks(y)
    ax.set_yticklabels([r["common_name"].split(" (")[0] or r["name"] for r in rows], fontsize=7.2)
    ax.set_xlim(0, 100)
    ax.set_xlabel("ignored disc light that would move the black hole\ndown to 5 M☉  [% of the optical light]")
    ax.set_title("(c) How far each black hole\nsits from the gap", fontsize=10)
    ax.grid(axis="y", visible=False)
    ps.caption(fig, "(a) Inclinations and masses from Kreidberg et al. 2012, Sect. 1 (Shahbaz et al. 1994; "
                    "Froning & Robinson 2001; Cantrell et al. 2010). (b)–(c): first-order scaling, not a light-curve model; "
                    "median f(M), q and i of each system (BlackCAT online).", y=-0.1)
    return ps.save(fig, FIGURE_DIR / "05_inclination_bias")


# ---------------------------------------------------------------------------
# 6. The stellar graveyard
# ---------------------------------------------------------------------------

def fig_graveyard(res, art):
    g = art["graveyard"]
    w = res["weighing"]
    draws = art["weighing"]["draws"]
    ns, gw = g["ns"], g["gw"].copy()
    fig, ax = plt.subplots(figsize=(13, 6.2))
    x = 0
    ticks = []
    # neutron stars weighed by light
    nsm = np.sort(np.concatenate([ns["mass"].values, [v[1] for v in lit.NS_XRB]]))
    xs = np.arange(len(nsm)) * 0.5
    ax.plot(xs, nsm, "o", ms=3.6, color=ps.C_NS_LIGHT, mec=ps.SURFACE, mew=0.6, zorder=4)
    ticks.append((xs.mean(), f"neutron stars\nweighed by light\n({len(nsm)})"))
    x = xs.max() + 6
    # black holes weighed by light
    bh = [(np.median(v["m1"]), k) for k, v in draws.items() if "m1" in v and not k.startswith("Gaia")]
    bh.sort()
    xb = x + np.arange(len(bh)) * 0.9
    ax.plot(xb, [m for m, _ in bh], "o", ms=5, color=ps.C_BH_LIGHT, mec=ps.SURFACE, mew=0.8, zorder=4)
    ticks.append((xb.mean(), f"black holes,\nX-ray binaries\n({len(bh)})"))
    for nm in ("Cyg X-1", "1H 1659-487"):
        i_ = [k for _, k in bh].index(nm)
        ax.text(xb[i_] + 0.8, bh[i_][0], "GX 339-4" if nm.startswith("1H") else nm, fontsize=7,
                color=ps.INK_SECONDARY, va="center")
    x = xb.max() + 9
    gb = sorted((np.median(draws[k]["m1"]), k) for k in draws if k.startswith("Gaia"))
    xg = x + np.arange(len(gb)) * 1.2
    ax.plot(xg, [m for m, _ in gb], "^", ms=6.5, color=ps.C_BH_GAIA, mec=ps.SURFACE, mew=0.8, zorder=4)
    for xx, (m, k) in zip(xg, gb):
        ax.text(xx + 0.7, m * (0.9 if k == "Gaia BH2" else 1.08 if k == "Gaia BH1" else 1.0),
                k.replace("Gaia ", ""), fontsize=7, color=ps.INK_SECONDARY, va="center")
    ticks.append((xg.mean(), "Gaia,\ndormant\n(3)"))
    x = xg.max() + 7
    # gravitational-wave events, sorted by remnant mass
    gw["sort_key"] = gw["mf"].fillna(gw["m1"] + gw["m2"])  # a few events lack a published remnant mass
    gw = gw.sort_values("sort_key").reset_index(drop=True)
    xw = x + np.arange(len(gw)) * 0.55
    for k, r in gw.iterrows():
        for m in (r["m1"], r["m2"]):
            ax.plot([xw[k], xw[k]], [m, r["mf"]], color=ps.CONTEXT, lw=0.5, zorder=1)
    cols1 = np.where(gw["m1"] < 2.5, ps.C_NS_GW, ps.C_BH_GW)
    cols2 = np.where(gw["m2"] < 2.5, ps.C_NS_GW, ps.C_BH_GW)
    ax.scatter(xw, gw["m1"], s=7, c=cols1, lw=0, zorder=3)
    ax.scatter(xw, gw["m2"], s=7, c=cols2, lw=0, zorder=3)
    ax.scatter(xw, gw["mf"], s=9, facecolors="none", edgecolors=ps.C_BH_GW, lw=0.6, zorder=3)
    ticks.append((xw.mean(), f"gravitational-wave mergers ({len(gw)}): two objects merge into one"))
    for name, dx, dy in (("GW170817", -2, 0.65), ("GW230529", 2, 1.25), ("GW231123", -14, 1.0),
                         ("GW190814", 0, 0.55)):
        hit = gw[gw["name"].str.startswith(name)]
        if len(hit):
            k = hit.index[0]
            m = hit["m2"].iloc[0] if name == "GW190814" else hit["mf"].iloc[0] if name == "GW231123" else hit["m1"].iloc[0]
            ps.annotate(ax, name, (xw[k], m), (xw[k] + dx, m * dy), fontsize=7)
    ax.axhspan(2.5, 5.0, color=ps.GAP_FILL, lw=0, zorder=0)
    ax.axhspan(50, 130, color=ps.GAP_FILL, lw=0, zorder=0)
    ax.text(xs.min(), 3.45, "lower mass gap 2.5–5 M☉: sparsely populated", fontsize=7.4,
            color=ps.INK_SECONDARY, va="center")
    ax.text(xs.min(), 82, "predicted pair-instability gap\n(~50–130 M☉): no star should leave\na black hole here",
            fontsize=7.4, color=ps.INK_SECONDARY, va="center")
    _log_mass_axis(ax, 0.9, 300, axis="y")
    ax.set_ylabel("mass  [M☉]")
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], fontsize=7.6)
    ax.set_xlim(-3, xw.max() + 3)
    ax.grid(axis="x", visible=False)
    hs = [ax.plot([], [], "o", color=ps.C_NS_LIGHT, label="neutron stars, weighed by light")[0],
          ax.plot([], [], "o", color=ps.C_BH_LIGHT, label="black holes, X-ray binaries")[0],
          ax.plot([], [], "^", color=ps.C_BH_GAIA, label="black holes, Gaia astrometry")[0],
          ax.plot([], [], "o", ms=3, color=ps.C_NS_GW, label="neutron stars, gravitational waves")[0],
          ax.plot([], [], "o", ms=3, color=ps.C_BH_GW, label="black holes, gravitational waves")[0],
          ax.plot([], [], "o", ms=3.5, mfc="none", mec=ps.C_BH_GW, label="merger remnant")[0]]
    ax.legend(handles=hs, loc="upper left", fontsize=7.2, ncol=2)
    ax.set_title("The stellar graveyard: everything weighed so far")
    ps.caption(fig, f"GWOSC cumulative event list (GWTC-1 to GWTC-5.0; {len(gw)} candidates with p_astro ≥ 0.5 and "
                    "published masses; medians, source frame). Neutron stars: Freire (pulsar timing), Özel & Freire 2016. "
                    "Black holes from light: this project. After the LIGO-Virgo-KAGRA 'stellar graveyard' plots.", y=0.0)
    return ps.save(fig, FIGURE_DIR / "06_stellar_graveyard")


# ---------------------------------------------------------------------------
# 7. Why the gravitational-wave black holes look heavier
# ---------------------------------------------------------------------------

def fig_light_vs_gravity(res, art):
    g = art["graveyard"]
    bbh = g["bbh"]
    gst = res["graveyard"]["gw_stats"]
    comp = np.concatenate([bbh["m1"], bbh["m2"]])
    wv = np.concatenate([bbh["w_vmax"], bbh["w_vmax"]])
    wc = np.concatenate([bbh["w_chirp"], bbh["w_chirp"]])
    xrb = g["xrb_mix"]
    bins = np.geomspace(3, 200, 40)
    c = np.sqrt(bins[1:] * bins[:-1])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=False)

    def dens(x, w=None):
        h, _ = np.histogram(x, bins=bins, weights=w)
        return h / h.sum()

    for ax, title, gwd, gwlab in (
            (axes[0], "(a) As detected", dens(comp), f"gravitational waves, as detected ({len(bbh)} binaries)"),
            (axes[1], "(b) Weighted by the volume each could be heard in", dens(comp, wv),
             "gravitational waves, weighted by 1/V_max")):
        ax.step(c, dens(xrb), where="mid", color=ps.C_BH_LIGHT, lw=2.0,
                label=f"X-ray binaries (each system once; {len(res['graveyard']['xrb_stats']['names'])})")
        ax.fill_between(c, 0, dens(xrb), step="mid", color=ps.C_BH_LIGHT, alpha=0.12, lw=0)
        ax.step(c, gwd, where="mid", color=ps.C_BH_GW, lw=2.0, label=gwlab)
        if ax is axes[1]:
            ax.step(c, dens(comp, wc), where="mid", color=ps.C_BH_GW, lw=1.0, ls=":",
                    label="weighted by M_chirp^−5/2 instead")
        for gb in res["weighing"]["gaia"]:
            ax.plot(gb["m_bh"]["median"], 0.004, "^", ms=6, color=ps.C_BH_GAIA, zorder=5)
        ax.plot([], [], "^", color=ps.C_BH_GAIA, label="Gaia BH1, BH2, BH3 (dormant)")
        _log_mass_axis(ax, 3, 200)
        ax.set_xlabel("black hole mass  [M☉]   (both members of each merging pair)")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(0, None)
        ax.legend(loc="upper right", fontsize=7.2)
    axes[0].set_ylabel("fraction per bin")
    o, v = gst["observed"]["all"], gst["w_vmax"]["all"]
    xq = res["graveyard"]["xrb_stats"]["quantiles"]
    neff = gst["w_vmax"]["n_eff"]
    axes[0].text(0.97, 0.80, f"medians: X-ray binaries {xq[1]:.1f} M☉,\ngravitational waves {o[1]:.0f} M☉",
                 transform=axes[0].transAxes, fontsize=7.6, va="top", ha="right", color=ps.INK_SECONDARY)
    axes[1].text(0.97, 0.62, f"weighted median {v[1]:.1f} M☉; 90% below {v[2]:.0f} M☉\n"
                 f"(as detected: 90% below {o[2]:.0f} M☉;\nX-ray binaries: 90% below {xq[2]:.0f} M☉)\n"
                 f"the weights leave ~{neff:.0f} effective events,\nhence the jagged histogram",
                 transform=axes[1].transAxes, fontsize=7.6, va="top", ha="right", color=ps.INK_SECONDARY)
    fig.suptitle("Black holes heard merging look heavier mostly because heavy mergers are louder", x=0.07,
                 ha="left", y=1.02)
    ps.caption(fig, "V_max: comoving volume (Planck 2018) out to d_L × SNR/8 for each event; every observing run treated "
                    "alike. Approximate: it ignores how redshift changes the signal. X-ray binaries: our mass "
                    "distributions, every system weighted equally.", y=-0.03)
    return ps.save(fig, FIGURE_DIR / "07_light_vs_gravity")


# ---------------------------------------------------------------------------
# 8-9. Outbursts
# ---------------------------------------------------------------------------

def _dates(mjd):
    return Time(np.asarray(mjd, float), format="mjd").to_datetime()


def fig_outburst(res, art):
    df = art["outbursts"]["j1820"]
    o = res["outbursts"]["j1820"]
    tr = o["transitions"]
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.2), sharex=True,
                             gridspec_kw={"height_ratios": [1.5, 1, 1], "hspace": 0.12})
    ax = axes[0]
    t = _dates(df["mjd"])
    def broken(mjd, y):
        """Insert gaps where consecutive points are more than 3 days apart."""
        mjd, y = np.asarray(mjd, float), np.asarray(y, float)
        cut = np.where(np.diff(mjd) > 3)[0] + 1
        mjd = np.insert(mjd, cut, mjd[cut - 1] + 1)
        y = np.insert(y, cut, np.nan)
        return _dates(mjd), y

    for band in ("2-4", "4-10", "10-20", "bat"):
        y, e = df[f"crab_{band}"], df[f"err_{band}"]
        ok = np.isfinite(y) & (y > 2 * e)
        tt, yy = broken(df["mjd"][ok], y[ok])
        ax.plot(tt, yy, "o-", ms=2.2, lw=0.8, color=ps.BAND_COLOURS[band], label=ps.BAND_LABELS[band])
    ax.set_yscale("log")
    ax.set_ylim(0.01, 6)
    ax.yaxis.set_major_locator(mticker.FixedLocator([0.01, 0.03, 0.1, 0.3, 1, 3]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_ylabel("flux  [Crab]")
    ax.legend(loc="upper right", fontsize=7.4, ncol=2)
    ax.set_title("MAXI J1820+070, 2018: a black hole's disc switches from hard to soft and back")
    ps.annotate(ax, f"peak of the hard state: {o['peak_bat']['crab']:.1f} Crab in hard X-rays",
                (_dates([o["peak_bat"]["mjd"]])[0], o["peak_bat"]["crab"]),
                (_dates([o["peak_bat"]["mjd"] + 25])[0], 4.6), fontsize=7.4, relpos=(0, 0.5))
    ax = axes[1]
    hr = df["hr"]
    ax.plot(t, hr, "o", ms=2.4, color=ps.INK_SECONDARY)
    ax.axhline(tr["hr_hard"], color=ps.AXIS, lw=0.8)
    ax.axhline(tr["hr_soft"], color=ps.AXIS, lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 3)
    ax.yaxis.set_major_locator(mticker.FixedLocator([0.1, 0.2, 0.5, 1, 2]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_ylabel("hardness (4–10)/(2–4) keV\n[Crab = 1]")
    h0, h1 = _dates(J1820_HOMAN_TRANSITION)
    for a_ in axes:
        a_.axvspan(h0, h1, color=ps.GAP_FILL, lw=0, zorder=0)
        a_.axvline(_dates([tr["hard_to_soft_mjd"]])[0], color=ps.C_OURS, lw=1.1)
        a_.axvline(_dates([tr["soft_to_hard_mjd"]])[0], color=ps.C_OURS, lw=1.1)
    ax.text(_dates([tr["hard_to_soft_mjd"] + 2])[0], 2.2, f"hard → soft\n{Time(tr['hard_to_soft_mjd'], format='mjd').iso[:10]}",
            fontsize=7.2, color=ps.INK_SECONDARY, va="top")
    ax.text(_dates([tr["soft_to_hard_mjd"] - 2])[0], 1.3, f"soft → hard\n{Time(tr['soft_to_hard_mjd'], format='mjd').iso[:10]}",
            fontsize=7.2, color=ps.INK_SECONDARY, va="top", ha="right")
    ax.text(_dates([58200])[0], 0.13, "hard state: a hot corona dominates", fontsize=7.2, color=ps.INK_SECONDARY)
    ax.text(_dates([58318])[0], 0.13, "soft state: the thermal disc dominates", fontsize=7.2, color=ps.INK_SECONDARY)
    ax = axes[2]
    lam = df["eflux_2-50"] * (o["peak_2_50_kev"]["l_ledd"]["median"] / o["peak_2_50_kev"]["eflux_cgs"])
    okl = np.isfinite(lam) & (lam > 0)
    tt, yy = broken(df["mjd"][okl], lam[okl])
    ax.plot(tt, yy, "o-", ms=2.2, lw=0.8, color=ps.INK)
    ax.set_yscale("log")
    ax.set_ylim(3e-4, 0.3)
    ax.yaxis.set_major_locator(mticker.FixedLocator([0.001, 0.003, 0.01, 0.03, 0.1, 0.3]))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_ylabel("luminosity 2–50 keV\n[% of Eddington]")
    p = o["peak_2_50_kev"]["l_ledd"]
    s = o["soft_to_hard"]["l_ledd"]
    ax.text(_dates([58212])[0], p["median"] * 1.15, f"peak {100 * p['median']:.0f}% of Eddington",
            fontsize=7.4, color=ps.INK_SECONDARY)
    ps.annotate(ax, f"switches back at {100 * s['median']:.1f}%", (_dates([tr["soft_to_hard_mjd"]])[0], s["median"]),
                (_dates([tr["soft_to_hard_mjd"] + 12])[0], s["median"] * 4), fontsize=7.4, relpos=(0, 0.5))
    ax.set_xlabel("2018")
    import matplotlib.dates as mdates
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ps.caption(fig, "MAXI/GSC (RIKEN, JAXA) and Swift/BAT transient monitor (Krimm et al. 2013) daily fluxes, each divided "
                    "by the Crab measured the same year. Orange lines: our transitions from the hardness; grey band: the "
                    "rapid transition of Homan et al. 2020 (July 3–9). Eddington ratio for 2.96 kpc (VLBI parallax, Atri "
                    "et al. 2020) and 8.48 M☉ (Torres et al. 2020), 2–50 keV only.", y=0.03)
    return ps.save(fig, FIGURE_DIR / "08_outburst_lightcurve")


def fig_hid(res, art):
    df = art["outbursts"]["j1820"]
    gx = art["outbursts"]["gx339"]
    outs = [o for o in res["outbursts"]["gx339"] if o["kind"] == "full"]
    fig = plt.figure(figsize=(12.5, 6.2))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.7, 1, 1, 1], hspace=0.35, wspace=0.14)
    ax = fig.add_subplot(gs[:, 0])

    def log_ticks(a, xt, yt):
        for axis, t in ((a.xaxis, xt), (a.yaxis, yt)):
            axis.set_major_locator(mticker.FixedLocator(t))
            axis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
            axis.set_minor_formatter(mticker.NullFormatter())

    ok = np.isfinite(df["hr"]) & np.isfinite(df["crab_2-20"]) & (df["crab_2-20"] > 0)
    d = df[ok]
    tt = d["mjd"].values
    sc = ax.scatter(d["hr"], d["crab_2-20"], c=tt - tt.min(), cmap="Blues", vmin=-60, vmax=tt.max() - tt.min(),
                    s=14, lw=0.3, edgecolors=ps.INK_SECONDARY, zorder=3)
    ax.plot(d["hr"], d["crab_2-20"], color=ps.CONTEXT, lw=0.6, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.12, 2.5)
    ax.set_ylim(0.03, 4)
    log_ticks(ax, [0.2, 0.5, 1, 2], [0.05, 0.1, 0.2, 0.5, 1, 2])
    ax.set_xlabel("hardness (4–10)/(2–4) keV  [Crab = 1]")
    ax.set_ylabel("MAXI 2–20 keV flux  [Crab]")
    ax.set_title("(a) MAXI J1820+070, 2018", fontsize=10)
    for txt, xy in (("rise in the\nhard state", (1.5, 0.35)), ("hard → soft\nin days", (0.55, 2.6)),
                    ("soft state:\nslow decay", (0.21, 0.27)), ("back to hard\nand fading", (1.2, 0.06))):
        ax.text(*xy, txt, fontsize=7.4, color=ps.INK_SECONDARY, ha="center")
    cax = ax.inset_axes([0.06, 0.07, 0.36, 0.028])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("days since the outburst began (light → dark)", fontsize=7.2)
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_visible(False)
    for k, o in enumerate(outs[:6]):
        axk = fig.add_subplot(gs[k // 3, 1 + k % 3])
        seg = gx[(gx["mjd"] >= o["window_mjd"][0]) & (gx["mjd"] <= o["window_mjd"][1])]
        seg = seg[np.isfinite(seg["hr"]) & (seg["crab_2-20"] > 0)]
        axk.plot(seg["hr"], seg["crab_2-20"], color=ps.CONTEXT, lw=0.6, zorder=2)
        axk.scatter(seg["hr"], seg["crab_2-20"], c=seg["mjd"] - seg["mjd"].min(), cmap="Blues", vmin=-60,
                    vmax=seg["mjd"].max() - seg["mjd"].min(), s=7, lw=0.2, edgecolors=ps.INK_SECONDARY, zorder=3)
        tr = o.get("transitions")
        if tr:
            for key, mk in (("hard_to_soft", "o"), ("soft_to_hard", "s")):
                f = o.get(f"{key}_crab_2_20")
                if f is not None and np.isfinite(f):
                    h = seg.loc[(np.abs(seg["mjd"] - tr[f"{key}_mjd"]) <= 1), "hr"].mean()
                    axk.plot(h, f, mk, ms=6, mfc="none", mec=ps.C_OURS, mew=1.4, zorder=5)
        axk.set_xscale("log")
        axk.set_yscale("log")
        axk.set_xlim(0.12, 2.8)
        axk.set_ylim(0.008, 1.2)
        log_ticks(axk, [0.2, 0.5, 1, 2], [0.01, 0.1, 1])
        y0, y1 = Time(o["start_mjd"], format="mjd").datetime.year, Time(o["end_mjd"], format="mjd").datetime.year
        axk.set_title(f"GX 339-4, {y0}" + (f"–{str(y1)[2:]}" if y1 != y0 else ""), fontsize=8.6)
        if k % 3:
            axk.set_yticklabels([])
        if k < 3:
            axk.set_xticklabels([])
    axk = fig.axes[-1]
    h1, = axk.plot([], [], "o", mfc="none", mec=ps.C_OURS, mew=1.4, label="hard → soft")
    h2, = axk.plot([], [], "s", mfc="none", mec=ps.C_OURS, mew=1.4, label="soft → hard")
    fig.legend(handles=[h1, h2], loc="lower right", bbox_to_anchor=(0.9, 0.0), ncol=2, fontsize=7.4)
    fig.suptitle("Hardness–intensity diagrams: every outburst runs the same anticlockwise loop", x=0.07, ha="left")
    n_full = len(outs)
    n_all = len(res["outbursts"]["gx339"])
    ps.caption(fig, f"MAXI/GSC daily fluxes in Crab units. GX 339-4: the {n_full} outbursts of {n_all} since 2009 that "
                    "reached the soft state (the rest stayed hard). Open markers: our transition days. Each loop "
                    "switches back to hard at a flux 4–18 times lower than the one it left at: hysteresis.", y=0.01)
    return ps.save(fig, FIGURE_DIR / "09_hardness_intensity")


# ---------------------------------------------------------------------------
# 10. Gaia's dormant black holes
# ---------------------------------------------------------------------------

def fig_gaia(res, art):
    gres = res["gaia"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), gridspec_kw={"width_ratios": [1.15, 1, 1], "wspace": 0.25})
    ax = axes[0]
    a = art["gaia"]["Gaia BH1"]
    row, dr = a["row"], a["draws"]
    for k in range(120):
        xa, yd = gaia_orbits.orbit_on_sky(dr["A"][k], dr["B"][k], dr["F"][k], dr["G"][k], dr["e"][k])
        ax.plot(xa, yd, color=ps.CONTEXT, lw=0.5, alpha=0.6)
    xa, yd = gaia_orbits.orbit_on_sky(row["a_thiele_innes"], row["b_thiele_innes"], row["f_thiele_innes"],
                                      row["g_thiele_innes"], row["eccentricity"])
    ax.plot(xa, yd, color=ps.C_OURS, lw=2)
    m_star, m_bh = 0.93, gres["Gaia BH1"]["m2"]["median"]
    ax.plot(-xa * m_star / m_bh, -yd * m_star / m_bh, color=ps.INK, lw=1.4)
    ax.plot(0, 0, "+", color=ps.INK, ms=8)
    ax.set_aspect("equal")
    ax.invert_xaxis()
    ax.set_xlabel("Δα cos δ  [mas]  (east to the left)")
    ax.set_ylabel("Δδ  [mas]")
    ax.set_title("(a) Gaia BH1: the star's orbit on the sky", fontsize=10)
    ax.set_ylim(-3.6, 2.6)
    h1, = ax.plot([], [], color=ps.C_OURS, lw=2, label="the Sun-like star (DR3 solution)")
    h2, = ax.plot([], [], color=ps.CONTEXT, lw=1, label="120 draws from its covariance")
    h3, = ax.plot([], [], color=ps.INK, lw=1.4, label=f"the black hole, {m_bh / m_star:.0f}× closer to the centre")
    ax.legend(handles=[h1, h2, h3], loc="lower center", fontsize=6.8, ncol=1)
    for ax, name in ((axes[1], "Gaia BH1"), (axes[2], "Gaia BH2")):
        r = gres[name]
        m2 = art["gaia"][name]["draws"]["m2"]
        bins = np.linspace(4, 25, 120)
        ax.hist(m2, bins=bins, density=True, histtype="stepfilled", color=ps.C_OURS, alpha=0.22, lw=0)
        ax.hist(m2, bins=bins, density=True, histtype="step", color=ps.C_OURS, lw=1.6,
                label=f"ours, Gaia DR3 alone: {r['m2']['median']:.1f} +{r['m2']['plus']:.1f} −{r['m2']['minus']:.1f} M☉")
        chk = r["check"]
        peak = np.histogram(m2, bins=bins, density=True)[0].max()
        ax.set_ylim(0, peak * 1.95)
        top = peak / 0.62 * 1.0
        if "m2" in chk:
            ax.errorbar(chk["m2"][0], peak * 1.1, xerr=chk["m2"][1], fmt="D", ms=4.5, mfc="none",
                        mec=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED, capsize=0,
                        label=f"El-Badry et al., astrometry only: {chk['m2'][0]} ± {chk['m2'][1]}")
        if "m2_for_mstar_1" in chk:
            ax.errorbar(chk["m2_for_mstar_1"][0], peak * 1.1, xerr=chk["m2_for_mstar_1"][1], fmt="D", ms=4.5,
                        mfc="none", mec=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED, capsize=0,
                        label=f"El-Badry et al., DR3, M★ = 1: {chk['m2_for_mstar_1'][0]} ± {chk['m2_for_mstar_1'][1]}")
        ax.errorbar(chk["joint_m2"][0], peak * 1.24, xerr=chk["joint_m2"][1], fmt="D", ms=4.5,
                    color=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED, capsize=0,
                    label=f"with follow-up velocities: {chk['joint_m2'][0]} ± {chk['joint_m2'][1]}")
        ax.set_xlim(4, 25)
        ax.set_yticks([])
        ax.set_xlabel("mass of the dark companion  [M☉]")
        f = r["f_ast"]
        ax.set_title(f"({'b' if name.endswith('1') else 'c'}) {name}: f_ast = {f['median']:.1f} M☉ from DR3",
                     fontsize=10)
        ax.legend(loc="upper right", fontsize=6.8)
    ps.caption(fig, "Gaia DR3 nss_two_body_orbit (Gaia Collaboration 2023): Thiele-Innes constants with their full "
                    "covariance → a₀, i → astrometric mass function (a₀/ϖ)³/P² = M₂³/(M₁+M₂)². Star masses 0.93 ± 0.05 "
                    "and 1.07 ± 0.19 M☉ (El-Badry et al. 2023a,b). (a) Positions relative to the centre of mass; "
                    "1 mas = 0.48 AU at Gaia BH1's 480 pc.", y=-0.03)
    return ps.save(fig, FIGURE_DIR / "10_gaia_black_holes")
