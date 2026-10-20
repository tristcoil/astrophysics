"""Every figure of the project, from the cached analysis (run.py --figures-only)."""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

import agefit as A  # noqa: E402
import config as C  # noqa: E402
import cosmology as CO  # noqa: E402
import plotstyle as S  # noqa: E402
from fiducial import turnoff_and_sgb_point  # noqa: E402
from isochrones import load_grid  # noqa: E402

S.use_style()
OUT = C.FIGURE_DIR
MODELS = list(C.MODEL_SETS)
TITLE = {"M92": "M92 (NGC 6341)", "NGC 6397": "NGC 6397", "M15": "M15 (NGC 7078)", "47 Tuc": "47 Tucanae (NGC 104)"}


def _draw(inp: dict) -> A.Draw:
    return A.Draw(inp["feh"], inp["alpha"], inp["dm"], inp["ebv"])


def _iso(model: str, prep: dict, inp: dict, ages, dc: float = 0.0):
    """Apparent isochrones (F606W, colour) of one model set at the cluster's inputs."""
    f = A.ClusterFitter(load_grid(model), prep["ridge"], feh0=inp["feh"])
    f.m_lo, f.m_hi = -10, 40  # keep the whole isochrone for drawing
    m, c, _ = f.models(_draw(inp), np.atleast_1d(ages))
    return m, c + dc


# ---------------------------------------------------------------------------
# 1. The four clusters
# ---------------------------------------------------------------------------


def fig_cmd_gallery(results, art):
    fig, axes = plt.subplots(1, 4, figsize=(15, 6.4), sharey=False)
    for ax, name in zip(axes, C.CORE_CLUSTERS):
        prep, inp = art["prepared"][name], art["inputs"][name]
        res = results["clusters"][name]
        ax.scatter(prep["c"], prep["m606"], s=0.4, c=S.C_STARS, alpha=0.35, lw=0, rasterized=True)
        mist = res["models"]["MIST"]
        m, c = _iso("MIST", prep, inp, [mist["nominal_age"]], mist["nominal_dc"])
        bright = np.percentile(prep["m606"], 0.05) - 0.3  # stop where the unsaturated stars stop
        keep = m[0] > bright
        ax.plot(np.where(keep, c[0], np.nan), m[0], color=S.C_OURS, lw=1.4)
        r = prep["ridge"]
        ax.set_xlim(r.c_to - 0.35, r.c_to + 0.85)
        ax.set_ylim(r.m_to + 4.6, max(r.m_to - 6.8, bright - 0.6))
        ax.set_xlabel("colour  F606W − F814W")
        ax.set_title(TITLE[name], fontsize=10, pad=22)
        ax.text(0.0, 1.015, f"{inp['dist_kpc']:.2f} kpc · [Fe/H] = {inp['feh']:+.2f} · "
                f"{res['photometry']['n_selected']:,} stars · MIST {mist['nominal_age']:.1f} Gyr",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, color=S.INK_SECONDARY)
    axes[0].set_ylabel("brightness  F606W (mag, fainter ↓)")
    a = axes[0]
    r = art["prepared"]["M92"]["ridge"]
    S.annotate(a, "main sequence\n(burning hydrogen)", (r.c_to + 0.16, r.m_to + 2.0), (r.c_to + 0.42, r.m_to + 3.2))
    S.annotate(a, "turnoff", (r.c_to, r.m_to), (r.c_to - 0.30, r.m_to + 1.2))
    S.annotate(a, "subgiant branch", (r.c_to + 0.10, r.m_to - 0.45), (r.c_to + 0.28, r.m_to + 0.6))
    S.annotate(a, "red giant branch", (r.c_to + 0.24, r.m_to - 3.0), (r.c_to + 0.40, r.m_to - 2.2))
    S.annotate(a, "horizontal branch\n(burning helium)", (r.c_to - 0.27, r.m_to - 3.45), (r.c_to - 0.33, r.m_to - 4.25),
               va="center")
    S.annotate(a, "blue stragglers", (r.c_to - 0.12, r.m_to - 1.25), (r.c_to - 0.33, r.m_to - 1.95))
    fig.suptitle("Four globular clusters seen by Hubble, with the best-fitting MIST isochrone", x=0.01, ha="left")
    S.caption(fig, "HUGS photometry (Nardiello et al. 2018), method 1: well-measured proper-motion members, corrected "
              "for differential reddening. Blue: MIST v2.5 isochrone ([α/Fe] = +0.4) at the fitted age, placed at "
              "the Baumgardt & Vasiliev (2021) distance and Harris (2010) reddening, with the fitted colour shift. "
              "The brightest giants saturate the detectors and are not shown.", y=-0.03)
    return S.save(fig, OUT / "01_cmd_gallery")


# ---------------------------------------------------------------------------
# 2. How the turnoff tells the time
# ---------------------------------------------------------------------------


def fig_turnoff_clock(results, art):
    name = "M92"
    prep, inp = art["prepared"][name], art["inputs"][name]
    res = results["clusters"][name]
    r = prep["ridge"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw={"width_ratios": [1.1, 1]})
    a1.scatter(prep["c"], prep["m606"], s=0.5, c=S.C_STARS, alpha=0.3, lw=0, rasterized=True)
    ages = [10, 11.5, 13, 14.5, 16]
    m, c = _iso("MIST", prep, inp, ages, res["models"]["MIST"]["nominal_dc"])
    for k, (age, col) in enumerate(zip(ages, S.AGE_RAMP)):
        a1.plot(c[k], m[k], color=col, lw=1.4, label=f"{age:g} Gyr")
    a1.legend(loc="upper left", fontsize=8, title="MIST isochrone age", title_fontsize=8)
    a1.set_xlim(r.c_to - 0.12, r.c_to + 0.28)
    a1.set_ylim(r.m_to + 1.6, r.m_to - 2.2)
    a1.set_xlabel("colour  F606W − F814W")
    a1.set_ylabel("F606W (mag)")
    a1.set_title("Older clusters have fainter turnoffs", fontsize=10)

    # right: the subgiant point's absolute magnitude against age, per model set
    c_to, m_to, m_sgb = turnoff_and_sgb_point(r.c, r.m)
    a606 = res["photometry"]["a606_per_ebv"]
    M_obs = m_sgb - inp["dm"] - a606 * inp["ebv"]
    M_err = np.hypot(inp["dm_err"], a606 * inp["ebv_err"])
    grid_ages = np.arange(8.0, 16.01, 0.25)
    for model in MODELS:
        g = load_grid(model)
        f = A.ClusterFitter(g, r, feh0=inp["feh"])
        f.m_lo, f.m_hi = -10, 40  # absolute magnitudes: keep everything
        ages_ok = grid_ages[grid_ages <= f.age_max + 1e-6]
        mm, cc, _ = f.models(A.Draw(inp["feh"], inp["alpha"], 0.0, 0.0), ages_ok)
        msgb = np.array([turnoff_and_sgb_point(cc[k], mm[k])[2] for k in range(len(ages_ok))])
        a2.plot(ages_ok, msgb, color=S.MODEL_COLOURS[model], lw=1.8)
        a2.text(ages_ok[-1] + 0.1, msgb[-1], model, color=S.INK_SECONDARY, fontsize=8, va="center")
    a2.axhspan(M_obs - M_err, M_obs + M_err, color=S.BLUE_150, alpha=0.6, lw=0)
    a2.axhline(M_obs, color=S.INK, lw=1.0)
    a2.text(8.15, M_obs + M_err + 0.012, f"M92's subgiant point: M = {M_obs:.2f} ± {M_err:.2f}\n"
            "(band: distance and reddening errors)", fontsize=8, color=S.INK, va="top")
    a2.invert_yaxis()
    a2.set_xlim(8, 17.2)
    a2.set_xlabel("age (Gyr)")
    a2.set_ylabel("absolute F606W of the subgiant point")
    a2.set_title("…by about 0.07 mag per billion years", fontsize=10)
    fig.suptitle("The turnoff clock, illustrated on M92", x=0.01, ha="left")
    S.caption(fig, "Left: MIST isochrones of 10-16 Gyr at M92's distance, reddening and composition. Stars leave the "
              "main sequence when their core hydrogen runs out; lower-mass stars take longer, so the turnoff and "
              "subgiant branch fade as the cluster ages. Right: the classic age indicator, the brightness of the "
              "subgiant point 0.05 mag redder than the turnoff (Chaboyer et al. 1996), for each isochrone set. "
              "A 0.1 mag (5%) error in distance shifts the age by 1.2-1.4 Gyr.", y=-0.04)
    return S.save(fig, OUT / "02_turnoff_clock")


# ---------------------------------------------------------------------------
# 3. Three model sets on each ridge line
# ---------------------------------------------------------------------------


def fig_three_models(results, art):
    fig, axes = plt.subplots(1, 4, figsize=(15, 5.6))
    for ax, name in zip(axes, C.CORE_CLUSTERS):
        prep, inp = art["prepared"][name], art["inputs"][name]
        res = results["clusters"][name]
        r = prep["ridge"]
        ax.scatter(prep["c"], prep["m606"], s=0.4, c=S.MUTED, alpha=0.25, lw=0, rasterized=True)
        for model in MODELS:
            mres = res["models"][model]
            m, c = _iso(model, prep, inp, [mres["nominal_age"]], mres["nominal_dc"])
            ax.plot(c[0], m[0], color=S.MODEL_COLOURS[model], lw=1.3,
                    label=f"{model} {mres['nominal_age']:.1f} Gyr (Δc {mres['nominal_dc']:+.3f})")
        used = art["nominal"][(name, "MIST")]
        ax.plot(r.c, r.m, "o", ms=2.2, color=S.INK, mfc=S.INK, mew=0, zorder=5)
        ax.set_xlim(r.c_to - 0.07, r.c_to + 0.32)
        ax.set_ylim(r.m_to + 2.2, r.m_to - 2.0)
        ax.set_xlabel("F606W − F814W")
        ax.set_title(TITLE[name], fontsize=10)
        ax.legend(loc="upper left", fontsize=7, handlelength=1.2)
        del used
    axes[0].set_ylabel("F606W (mag)")
    fig.suptitle("Three stellar-evolution codes, fitted to each cluster's ridge line (black dots)", x=0.01, ha="left")
    S.caption(fig, "Each set's best age at the central inputs (distance, reddening, [Fe/H], [α/Fe] = +0.4). Models are "
              "shifted in colour by Δc (fitted); their brightness is never adjusted. MIST v2.5 (Dotter et al. 2026), "
              "BaSTI-IAC (Pietrinferni et al. 2021), DSED (Dotter et al. 2008).", y=-0.03)
    return S.save(fig, OUT / "03_three_models")


# ---------------------------------------------------------------------------
# 4. Age distributions
# ---------------------------------------------------------------------------


def fig_age_distributions(results, art):
    fig, axes = plt.subplots(4, 1, figsize=(10, 9.6), sharex=True)
    x = np.linspace(8, 17, 600)
    t_u = results["cosmology"]["universe"].get("D3", {}).get("age", 13.8)
    for ax, name in zip(axes, C.CORE_CLUSTERS):
        res = results["clusters"][name]
        top = 0
        for model in MODELS:
            ages = art["mc"][(name, model)]["ages"]
            dens = gaussian_kde(ages, bw_method=0.25)(x)
            ax.plot(x, dens, color=S.MODEL_COLOURS[model], lw=1.6, label=model)
            top = max(top, dens.max())
            edge = res["models"][model]["mc_edge_fraction"]
            if edge > 0.05:
                k = int(np.argmax(dens * (x > 14.5)))
                ax.annotate(f"{model} grid ends at 15 Gyr:\n{100 * edge:.0f}% of draws pile up here", (x[k], dens[k]),
                            (x[k] + 0.35, dens[k] * 0.9), fontsize=7, color=S.INK_SECONDARY,
                            arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.8})
        comb = res["combined"]
        y0 = -0.18 * top
        ax.errorbar(comb["age"], y0, xerr=comb["total_err"], fmt="o", color=S.C_OURS, ms=6, capsize=0, lw=1.6)
        ax.text(comb["age"] + comb["total_err"] + 0.08, y0, f"this work {comb['age']:.1f} ± {comb['total_err']:.1f}",
                fontsize=8, va="center", color=S.INK)
        lit = results["literature"][name]
        k = 0
        for src, v in lit.items():
            if src == "this work":
                continue
            k += 1
            y = y0 - 0.16 * top * k
            ax.errorbar(v["age"], y, xerr=v["err"], fmt=S.LIT_MARKERS.get(src, "o"), color=S.C_PUBLISHED,
                        ms=5, capsize=0, lw=1.0, mfc="white")
            ax.text(v["age"] + v["err"] + 0.08, y, f"{src} {v['age']:.1f}", fontsize=7, va="center",
                    color=S.INK_SECONDARY)
        ax.axvline(t_u, color=S.INK, lw=0.9)
        ax.set_ylim(y0 - 0.16 * top * (k + 1), top * 1.25)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_title(TITLE[name], fontsize=10, loc="left")
    axes[0].text(t_u + 0.06, axes[0].get_ylim()[1] * 0.97, "age of the universe (D3, Planck)", fontsize=7.5,
                 va="top", color=S.INK)
    axes[0].legend(loc="upper left", ncol=3, fontsize=8)
    axes[-1].set_xlabel("age (Gyr)")
    axes[-1].set_xlim(8.5, 17)
    fig.subplots_adjust(top=0.93, hspace=0.35)
    fig.suptitle("Monte Carlo age distributions, three isochrone sets", x=0.01, ha="left")
    S.caption(fig, "Curves: ages from 1,000 draws of distance, reddening, [Fe/H], [α/Fe] and ridge-line noise per set. "
              "Blue point: mean of the three sets; error = input spread ⊕ model spread. Grey: published ages "
              "(VandenBerg et al. 2013 quote relative errors; their absolute uncertainty is ~±1.5 Gyr).", y=-0.02)
    return S.save(fig, OUT / "04_age_distributions")


# ---------------------------------------------------------------------------
# 5. Error budget
# ---------------------------------------------------------------------------


def fig_error_budget(results, art):
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.6), sharex=True)
    terms = ["distance", "reddening", "[Fe/H]", "[alpha/Fe]", "ridge statistics"]
    labels = ["distance", "reddening", "[Fe/H]", "[α/Fe]", "ridge noise", "isochrone set", "total"]
    for ax, name in zip(axes, C.CORE_CLUSTERS):
        res = results["clusters"][name]
        vals = [float(np.mean([res["models"][m]["budget"][t] for m in MODELS])) for t in terms]
        vals += [res["combined"]["models_err"], res["combined"]["total_err"]]
        y = np.arange(len(vals))[::-1]
        cols = [S.BLUE] * 6 + [S.INK]
        ax.barh(y, vals, color=cols, height=0.62)
        for yi, v in zip(y, vals):
            ax.text(v + 0.02, yi, f"{v:.2f}", va="center", fontsize=7.5, color=S.INK_SECONDARY)
        ax.set_yticks(y, labels if ax is axes[0] else [""] * len(y))
        ax.set_title(TITLE[name], fontsize=10)
        ax.grid(axis="y", visible=False)
    for ax in axes:
        ax.set_xlabel("age uncertainty (Gyr, 1σ)")
    fig.suptitle("What limits each age", x=0.01, ha="left", y=1.04)
    S.caption(fig, "Each input varied alone (others at their central values), averaged over the three isochrone sets. "
              "Distance: Baumgardt & Vasiliev (2021) plus a 0.025 mag systematic floor; reddening ±10% (≥0.01 mag); "
              "[Fe/H] ±0.1; [α/Fe] ±0.1. Isochrone set: spread of the three model ages. Total: all inputs varied "
              "together, combined in quadrature with the isochrone-set spread.", y=-0.05)
    return S.save(fig, OUT / "05_error_budget")


# ---------------------------------------------------------------------------
# 6. Age against metallicity
# ---------------------------------------------------------------------------


def fig_age_metallicity(results, art):
    import literature as L

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    vdb, feh = L.vandenberg2013(), L.vandenberg2013_feh()
    ids = [k for k in vdb if k in feh]
    ax.scatter([feh[k] for k in ids], [vdb[k][0] for k in ids], marker="s", s=16, color=S.MUTED, alpha=0.6,
               lw=0, label=f"VandenBerg et al. 2013 ({len(ids)} NGC clusters)")
    for name in C.CORE_CLUSTERS:
        res = results["clusters"][name]
        comb = res["combined"]
        inp = res["inputs"]
        ax.errorbar(inp["feh"], comb["age"], xerr=inp["feh_err"], yerr=comb["total_err"], fmt="o", ms=7,
                    color=S.C_OURS, capsize=0, lw=1.5, zorder=5)
        ax.text(inp["feh"] + 0.04, comb["age"] + 0.12, name, fontsize=8.5, color=S.INK)
    t_u = results["cosmology"]["universe"].get("D3", {}).get("age", 13.8)
    ax.axhline(t_u, color=S.INK, lw=0.9)
    ax.text(-0.25, t_u + 0.08, f"age of the universe, {t_u:.2f} Gyr (D3)", fontsize=8, ha="right")
    ax.set_xlabel("[Fe/H]  (iron relative to the Sun, log scale)")
    ax.set_ylabel("age (Gyr)")
    ax.set_xlim(-2.6, -0.2)
    ax.set_ylim(8.5, 16)
    ax.legend(loc="lower left")
    fig.suptitle("Metal-rich 47 Tuc is younger than the metal-poor clusters", x=0.01, ha="left")
    S.caption(fig, "Blue: this work (mean of three isochrone sets, total error). Grey squares: VandenBerg et al. (2013) "
              "relative ages for context. The metal-poor clusters formed first, from gas that fewer generations "
              "of stars had yet enriched with iron.", y=-0.03)
    return S.save(fig, OUT / "06_age_metallicity")


# ---------------------------------------------------------------------------
# 7. Chromosome maps
# ---------------------------------------------------------------------------


def fig_chromosome_maps(results, art):
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6))
    for ax, name in zip(axes, C.CORE_CLUSTERS):
        a = art["chromosome_maps"][name]
        r = results["chromosome_maps"][name]
        g1 = a["is_1g"]
        ax.scatter(a["d_x"][g1], a["d_c"][g1], s=9, marker="o", color=S.C_1G, lw=0, alpha=0.8)
        ax.scatter(a["d_x"][~g1], a["d_c"][~g1], s=11, marker="^", color=S.C_2G, lw=0, alpha=0.8)
        lit = r["milone2017"].get("f1G", [np.nan, np.nan])
        ax.text(0.03, 0.03, f"1G fraction {r['f1G']:.2f} ± {r['f1G_err']:.2f}\nMilone+17: {lit[0]:.2f} ± {lit[1]:.2f}\n"
                f"{r['n_rgb']} giants", transform=ax.transAxes, va="bottom", fontsize=7.5, color=S.INK_SECONDARY)
        ax.set_title(TITLE[name], fontsize=10)
        ax.set_xlabel("Δ F275W,F814W")
        if name == "47 Tuc":
            x1, y1 = np.median(a["d_x"][g1]), np.median(a["d_c"][g1])
            x2, y2 = np.median(a["d_x"][~g1]), np.median(a["d_c"][~g1])
            S.annotate(ax, "first generation", (x1, y1), (x1 - 0.55, y1 + 0.02))
            S.annotate(ax, "second generation\n(N-rich, O-poor)", (x2, y2), (x2 - 0.62, y2 + 0.20))
    axes[0].set_ylabel("Δ C F275W,F336W,F438W")
    fig.suptitle("Chromosome maps: every cluster holds two generations of stars", x=0.01, ha="left")
    S.caption(fig, "Red giants 1-3.5 mag above the turnoff, placed between the blue and red envelopes of the RGB in two "
              "ultraviolet colours (Milone et al. 2017 recipe). Circles: first generation, triangles: second "
              "generation, assigned by a Gaussian-mixture fit across the 1G sequence (after an 18° rotation). "
              "NGC 6397's field holds few giants.", y=-0.04)
    return S.save(fig, OUT / "07_chromosome_maps")


# ---------------------------------------------------------------------------
# 8. Timeline and 9. the age crisis
# ---------------------------------------------------------------------------


def fig_timeline(results, art):
    co = results["cosmology"]
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    rows = []
    for key, col in (("D3", S.C_CMB), ("D8", S.C_CC), ("D1", S.C_SN)):
        u = co["universe"].get(key)
        if u:
            rows.append((u["label"], u["age"], u["err"], col, "s"))
    for name in C.CORE_CLUSTERS:
        comb = results["clusters"][name]["combined"]
        rows.append((f"{name} (this work)", comb["age"], comb["total_err"], S.C_OURS, "o"))
    y = np.arange(len(rows))[::-1]
    for yi, (label, age, err, col, mk) in zip(y, rows):
        ax.errorbar(age, yi, xerr=err, fmt=mk, color=col, ms=7, capsize=0, lw=1.8)
        ax.text(age + err + 0.1, yi, f"{age:.2f} ± {err:.2f}", va="center", fontsize=8, color=S.INK_SECONDARY)
    ax.set_yticks(y, [r[0].replace("LCDM", "ΛCDM") for r in rows])
    n_u = sum(1 for r in rows if "this work" not in r[0])
    ax.axhline(len(rows) - n_u - 0.5, color=S.AXIS, lw=0.8)
    ax.text(16.4, len(rows) - 0.75, "the universe", ha="right", fontsize=8, color=S.INK_SECONDARY)
    ax.text(16.4, len(rows) - n_u - 0.75, "its oldest stars", ha="right", fontsize=8, color=S.INK_SECONDARY)
    t_need = co["t_universe_min"]["median"]
    ax.axvspan(t_need - co["t_universe_min"]["std"], t_need + co["t_universe_min"]["std"], color=S.BLUE_150,
               alpha=0.45, lw=0)
    ax.axvline(t_need, color=S.C_OURS, lw=1.0)
    ax.text(t_need, len(rows) - 0.35, f"universe must be ≥ {t_need:.1f} Gyr\n(metal-poor clusters + "
            f"{co['delay_gyr']['median']:.2f} Gyr to form)", ha="center", fontsize=7.5, color=S.INK)
    ax.set_xlim(9.5, 16.5)
    ax.set_ylim(-0.8, len(rows) + 0.6)
    ax.set_xlabel("age (Gyr)")
    ax.grid(axis="y", visible=False)
    fig.suptitle("The oldest stars against the age of the universe", x=0.01, ha="left")
    S.caption(fig, "Universe ages from earlier atlas projects (all assume flat ΛCDM): D3 Planck CMB fit, D8 cosmic "
              "chronometers, D1 supernovae with the Cepheid H0 = 73.5. Clusters: this work. A universe younger than "
              "its oldest stars would be a contradiction.", y=-0.04)
    return S.save(fig, OUT / "08_timeline")


def fig_age_crisis(results, art):
    co = results["cosmology"]
    fig, ax = plt.subplots(figsize=(9, 5.6))
    h0 = np.linspace(45, 90, 300)
    om = co["universe"].get("D3", {}).get("Om", 0.315)
    ax.plot(h0, CO.age_flat_lcdm(h0, om), color=S.INK, lw=1.6)
    ax.text(59.0, float(CO.age_flat_lcdm(56, om)), f"flat ΛCDM, Ωm = {om:.2f}\n(with dark energy)", fontsize=8.5,
            ha="left", va="center")
    ax.plot(h0, CO.age_flat_lcdm(h0, 1.0), color=S.INK_SECONDARY, lw=1.2)
    ax.text(56, float(CO.age_flat_lcdm(56, 1.0)) - 0.55, "no dark energy (Ωm = 1)", fontsize=8.5,
            color=S.INK_SECONDARY, ha="left")
    t = co["t_universe_min"]
    ax.axhspan(t["median"] - t["std"], t["median"] + t["std"], color=S.BLUE_150, alpha=0.5, lw=0)
    ax.axhline(t["median"], color=S.C_OURS, lw=1.2)
    ax.text(89.5, t["median"] + t["std"] + 0.1, "oldest clusters + time to form them (this work):\n"
            "the universe must be at least this old", ha="right", va="bottom", fontsize=8, color=S.INK)
    offsets = {"D3": (1.3, 0.25, "left"), "D8": (-5.6, 0.0, "right"), "D1": (-1.2, -0.75, "right")}
    for key, col in (("D3", S.C_CMB), ("D8", S.C_CC), ("D1", S.C_SN)):
        u = co["universe"].get(key)
        if u:
            ax.errorbar(u["H0"], u["age"], xerr=u["H0_err"], yerr=u["err"], fmt="s", color=col, ms=6, capsize=0, lw=1.4)
            dx, dy, ha = offsets[key]
            ax.text(u["H0"] + dx, u["age"] + dy, u["label"].replace("LCDM", "ΛCDM"), fontsize=8,
                    color=S.INK_SECONDARY, ha=ha)
    lim = co["h0_limits"].get("D3's Omega_m")
    if lim:
        ax.axvline(lim["H0_max_95"], color=S.INK, lw=0.8)
        ax.text(lim["H0_max_95"] + 0.4, 9.6, f"H0 < {lim['H0_max_95']:.1f}\n(95%, D3's Ωm)", fontsize=8)
    ax.add_patch(plt.Rectangle((70, 15.4), 10, 2.0, fill=False, ec=S.MUTED, lw=1.0))
    ax.text(75, 17.5, "early 1990s: H0 ~ 70-80, clusters 14-16 Gyr", ha="center", fontsize=7.5, color=S.MUTED)
    ax.set_xlim(45, 90)
    ax.set_ylim(7, 18)
    ax.set_xlabel("Hubble constant H0 (km/s/Mpc)")
    ax.set_ylabel("age of the universe (Gyr)")
    fig.suptitle("The age crisis, restaged", x=0.01, ha="left")
    S.caption(fig, "Curves: the age a flat universe would have for each H0. Without dark energy, H0 ≈ 70 gives under "
              "10 Gyr, far younger than the clusters: the 1990s crisis. With dark energy the CMB universe (D3) is "
              "older than the clusters; the Cepheid-calibrated H0 (D1) would make it younger.", y=-0.03)
    return S.save(fig, OUT / "09_age_crisis")


# ---------------------------------------------------------------------------
# 10. Differential reddening, 11. robustness
# ---------------------------------------------------------------------------


def fig_differential_reddening(results, art):
    cmap = LinearSegmentedColormap.from_list("div", S.DIVERGING)
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.0))
    for ax, name in zip(axes, C.CORE_CLUSTERS):
        prep = art["prepared"][name]
        dE = prep["dE"]
        lim = 0.03
        sc = ax.scatter(prep["x"] * 0.04, prep["y"] * 0.04, c=np.clip(dE, -lim, lim), cmap=cmap, vmin=-lim, vmax=lim,
                        s=1.0, lw=0, rasterized=True)
        ax.set_aspect("equal")
        ax.set_title(f"{TITLE[name]}\nrms ΔE(B−V) = {np.std(dE):.3f}", fontsize=9)
        ax.set_xlabel("x (arcsec)")
        ax.grid(False)
    axes[0].set_ylabel("y (arcsec)")
    cb = fig.colorbar(sc, ax=axes, shrink=0.85, pad=0.01)
    cb.set_label("ΔE(B−V) (mag)")
    fig.suptitle("Patchy dust: differential reddening across each field", x=0.01, ha="left")
    S.caption(fig, "Local reddening relative to the cluster mean, from the median offset of the 50 nearest main-sequence "
              "stars along the reddening vector (Milone et al. 2012). 0.04 arcsec per ACS/WFC pixel.", y=-0.05)
    return S.save(fig, OUT / "10_differential_reddening")


def fig_robustness(results, art):
    fig, axes = plt.subplots(1, 4, figsize=(15, 5.4), sharey=True, sharex=True)
    labels = None
    for ax, name in zip(axes, C.CORE_CLUSTERS):
        res = results["clusters"][name]
        rows = []
        for model in MODELS:
            m = res["models"][model]
            base_mc, base_fit = m["mc"]["median"], m["nominal_age"]
            v = m["variants"]
            entries = []
            if "geometric_distance" in v:
                entries.append(("geometric distance only (MC)", v["geometric_distance"]["median"] - base_mc))
            entries.append(("subgiant point only (MC)", v["sgb_point"]["median"] - base_mc))
            entries.append(("no colour shift (MC)", v["no_colour_shift"]["median"] - base_mc))
            for k, age in m["region_tests"].items():
                if not k.startswith("baseline"):
                    entries.append((k, age - base_fit))
            rows.append((model, entries))
        if labels is None:
            labels = [e[0] for e in rows[0][1]] + ["helium Y = 0.275 (BaSTI, MC)"]
        y = np.arange(len(labels))[::-1]
        ax.axvline(0, color=S.AXIS, lw=1.0)
        for model, entries in rows:
            d = dict(entries)
            ax.plot([d.get(lbl, np.nan) for lbl in labels], y, "o", color=S.MODEL_COLOURS[model], ms=5, label=model)
        ax.plot([res["helium_test"]["shift"]], [y[-1]], "o", color=S.MODEL_COLOURS["BaSTI"], ms=5)
        ax.set_title(TITLE[name], fontsize=10)
        ax.set_xlabel("age change (Gyr)")
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(np.arange(len(labels))[::-1], labels)
    axes[-1].legend(loc="lower right", fontsize=7.5)
    fig.suptitle("How much do the ages depend on the analysis choices?", x=0.01, ha="left", y=1.02)
    S.caption(fig, "Change from each isochrone set's own baseline: Monte Carlo medians (MC; 300 draws) against the "
              "1,000-draw baseline, single fits at the central inputs against the central fit. 'No colour shift' "
              "forces the model colours, which drags ages young: the models are 0.01-0.04 mag redder than the data. "
              "The geometric distances (Gaia parallax and kinematics only) are shorter for the metal-poor clusters.",
              y=-0.03)
    return S.save(fig, OUT / "11_robustness")


# ---------------------------------------------------------------------------
# 12-13. The survey extension (only if survey.py has been run)
# ---------------------------------------------------------------------------

ORIGIN_COLOURS = {"in situ": S.BLUE, "accreted": S.ORANGE, "unassociated": S.AQUA}


def _survey():
    import json

    path = C.DATA_DIR / "survey.json"
    return json.loads(path.read_text()) if path.exists() else None


def fig_survey_amr(results, art):
    sv = _survey()
    if not sv:
        return []
    cl = sv["clusters"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.8), gridspec_kw={"width_ratios": [1.5, 1]})
    core_ids = {C.CLUSTERS[n]["hugs"]: n for n in C.CORE_CLUSTERS}
    for org, col in ORIGIN_COLOURS.items():
        pts = [e for e in cl.values() if "age" in e and e["origin"] == org and not e["flags"]]
        if pts:
            a1.errorbar([e["inputs"]["feh"] for e in pts], [e["age"] for e in pts], yerr=[e["age_err"] for e in pts],
                        fmt="o", ms=6, color=col, lw=1.0, capsize=0, label=f"{org} ({len(pts)})")
    flagged = [e for e in cl.values() if "age" in e and e["flags"]]
    a1.scatter([e["inputs"]["feh"] for e in flagged], [min(max(e["age"], 7.7), 16.3) for e in flagged], s=22,
               facecolors="none", edgecolors=S.MUTED, lw=0.9, label=f"flagged ({len(flagged)}), any origin")
    for hid, name in core_ids.items():
        e = cl.get(hid)
        if e and "age" in e:
            a1.annotate(name, (e["inputs"]["feh"], e["age"]), (e["inputs"]["feh"] - 0.02, e["age"] - 0.75),
                        fontsize=7.5, color=S.INK, ha="right",
                        arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.6})
    t_u = results["cosmology"]["universe"].get("D3", {}).get("age", 13.8)
    a1.axhline(t_u, color=S.INK, lw=0.9)
    a1.text(-0.42, t_u + 0.1, "age of the universe (D3)", fontsize=8, ha="right")
    a1.set_xlabel("[Fe/H]")
    a1.set_ylabel("age (Gyr), mean of three isochrone sets")
    a1.set_ylim(7.5, 16.5)
    a1.legend(loc="lower left", fontsize=7.5, ncol=2)
    a1.set_xlim(-2.5, -0.35)
    a1.set_title("Born here or captured? Age against metallicity", fontsize=10)

    for ref, mk, col, label in (("vandenberg2013", "s", S.INK, "VandenBerg et al. 2013"),
                                ("valcin2025", "^", S.MUTED, "Valcin et al. 2025")):
        pts = [(e[ref][0], e["age"]) for e in cl.values() if "age" in e and e.get(ref) and not e["flags"]]
        if pts:
            xr, yr = zip(*pts)
            a2.plot(xr, yr, mk, ms=5, color=col, mfc="none" if col == S.MUTED else col, label=f"{label} ({len(pts)})")
    a2.plot([8, 16], [8, 16], color=S.AXIS, lw=1.0)
    summ = sv["summary"]
    note = []
    for ref, label in (("vs_vandenberg2013", "VandenBerg"), ("vs_valcin2025", "Valcin")):
        if ref in summ:
            note.append(f"{label}: median +{summ[ref]['median_diff']:.2f} Gyr, scatter {summ[ref]['rms_diff']:.2f}")
    a2.text(0.97, 0.03, "this work − published\n" + "\n".join(note), transform=a2.transAxes, ha="right",
            va="bottom", fontsize=7.5, color=S.INK_SECONDARY)
    a2.set_xlim(8.5, 15.5)
    a2.set_ylim(8.5, 15.5)
    a2.set_aspect("equal")
    a2.set_xlabel("published age (Gyr)")
    a2.set_ylabel("this work (Gyr)")
    a2.legend(loc="upper left", fontsize=7.5)
    a2.set_title("Against published ages (unflagged clusters)", fontsize=10)
    fig.suptitle("All 56 globular clusters of the HUGS survey", x=0.01, ha="left")
    S.caption(fig, "Origins from Massari, Koppelman & Helmi (2019): formed in the Milky Way's disc or bulge (in situ), "
              "brought in by Gaia-Enceladus, Sagittarius, Sequoia or the Helmi streams (accreted), or on orbits not tied to "
              "a known progenitor (unassociated). Hollow: flagged (high reddening, split subgiant branch, or poor fit).",
              y=-0.03)
    return S.save(fig, OUT / "12_survey_age_metallicity")


def fig_survey_gallery(results, art):
    sv = _survey()
    if not sv:
        return []
    import catalog as K_
    import photometry as P

    ids = [i for i, e in sorted(sv["clusters"].items(), key=lambda kv: kv[1]["inputs"]["feh"] if "inputs" in kv[1] else 0)
           if "age" in e]
    ncol = 8
    nrow = int(np.ceil(len(ids) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(16, 2.6 * nrow))
    fig.subplots_adjust(hspace=0.55, wspace=0.12, top=0.95, bottom=0.03, left=0.02, right=0.99)
    grid = load_grid("MIST")
    for ax, hid in zip(axes.flat, ids):
        e = sv["clusters"][hid]
        inp = K_.cluster_inputs(hid)
        prep = P.prepare_hugs(hid, inp)
        r = prep["ridge"]
        ax.scatter(prep["c"], prep["m606"], s=0.2, c=S.C_STARS, alpha=0.25, lw=0, rasterized=True)
        f = A.ClusterFitter(grid, r, feh0=inp["feh"])
        f.m_lo, f.m_hi = -10, 40
        nm = e["nominal"]["MIST"]
        m, c, _ = f.models(A.Draw(inp["feh"], inp["alpha"], inp["dm"], inp["ebv"]), [nm["age"]])
        ax.plot(c[0] + nm["dc"], m[0], color=S.C_OURS, lw=1.0)
        ax.set_xlim(r.c_to - 0.15, r.c_to + 0.45)
        ax.set_ylim(r.m_to + 2.5, r.m_to - 3.0)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        flag = "  ⚑" if e["flags"] else ""
        ax.set_title(f"{K_.hugs_to_ngc(hid)} [Fe/H] {inp['feh']:+.2f}\n{e['age']:.1f} Gyr{flag}", fontsize=7,
                     color=S.MUTED if e["flags"] else S.INK)
    for ax in list(axes.flat)[len(ids):]:
        ax.axis("off")
    fig.suptitle("Every HUGS cluster with its best MIST isochrone (sorted by [Fe/H])", x=0.01, ha="left", y=0.985)
    S.caption(fig, "Turnoff region of each cluster's cleaned CMD (F606W vs F606W − F814W) with the MIST fit at the "
              "central inputs. ⚑: flagged in survey.json (high reddening, Type II split subgiant branch, age at the "
              "grid edge, large colour shift or poor fit).", y=0.012)
    return S.save(fig, OUT / "13_survey_gallery", formats=("png",))


def make_all(results: dict, art: dict) -> list:
    written = []
    for fn in (fig_cmd_gallery, fig_turnoff_clock, fig_three_models, fig_age_distributions, fig_error_budget,
               fig_age_metallicity, fig_chromosome_maps, fig_timeline, fig_age_crisis, fig_differential_reddening,
               fig_robustness, fig_survey_amr, fig_survey_gallery):
        written += fn(results, art)
        print(f"  {fn.__name__}")
    import interactive

    written.append(interactive.build(results, art))
    print("  interactive explorer")
    return written
