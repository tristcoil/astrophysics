"""Static figures (PNG + SVG in figures/).

    python figures.py            # every figure whose inputs are cached

Each function takes the results it needs, so run.py can call them after the
analyses; figures never refit anything.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.time import Time

import config as C
import kepler as K
import plotstyle as S

S.use_style()
warnings.filterwarnings("ignore", message=".*dubious year.*")   # UTC after 2026: leap seconds unknown, harmless
PL = C.PLANETS                      # b, c, d, e (outside in) for legends
IO = C.INNER_TO_OUTER


def _year(mjd):
    return Time(np.asarray(mjd, float), format="mjd").decimalyear


def _mjd(year):
    return Time(np.asarray(year, float), format="decimalyear").mjd


def _ring(ax, x, y, size=26, color=S.INK, zorder=5, **kw):
    """Data marker: filled ink dot with a 2px surface ring."""
    ax.scatter(x, y, s=size, color=color, edgecolors=S.SURFACE, linewidths=1.2, zorder=zorder, **kw)


# ---------------------------------------------------------------------------
# beta Pic b
# ---------------------------------------------------------------------------

def betapic_figure(bp: dict, n_draws: int = 120):
    import betapic as B
    s = bp["_samples"]
    data = B.load_data()
    pan_med = np.median(s["pan"])                  # the orbit's line of nodes on the sky
    c, sn = np.cos(pan_med), np.sin(pan_med)

    def rotate(ra, dec):
        # along = position along the node line (PA = Omega), perp = across it
        return ra * sn + dec * c, ra * c - dec * sn

    pa = np.radians(data.pa.to_numpy())
    d_ra, d_dec = data.sep.to_numpy() * np.sin(pa), data.sep.to_numpy() * np.cos(pa)
    d_al, d_pe = rotate(d_ra, d_dec)
    d_year = _year(data.epoch.to_numpy())

    years = np.linspace(2003, 2046, 1500)
    t = _mjd(years)
    rng = np.random.default_rng(1)
    pick = rng.choice(len(s["a"]), n_draws, replace=False)
    fig = plt.figure(figsize=(11, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], width_ratios=[1.6, 1], hspace=0.42, wspace=0.22)
    ax = fig.add_subplot(gs[0, :])
    for k in pick:
        x, y, z = K.orbit_xyz(t, s["a"][k], s["e"][k], s["inc"][k], s["aop"][k], s["pan"][k], s["tau"][k], s["m"][k])
        al, pe = rotate(x * s["plx"][k], y * s["plx"][k])
        ax.plot(years, al, color=S.ORANGE, lw=0.6, alpha=0.08)
    ax.errorbar(d_year, d_al, yerr=data.sep_err.to_numpy(), fmt="none", ecolor=S.INK_SECONDARY, lw=0.8, zorder=4)
    _ring(ax, d_year, d_al, size=22)
    ax.axhline(0, color=S.AXIS, lw=0.8)
    ymin, ymax = ax.get_ylim()
    for key, label, front in (("conj2017_mjd", "passes in front\n2017 Sep 12", True),
                              ("behind_year", "behind the star\n{:.1f}", False),
                              ("next_front_year", "in front again\n{:.1f}", True)):
        v = bp[key][0]
        yr = _year(v) if key.endswith("mjd") else v
        ax.axvline(yr, color=S.VIOLET if front else S.MUTED, lw=1.0, zorder=1)
        ax.text(yr + 0.3, ymax * 0.97, label.format(yr), va="top", fontsize=7.5, color=S.INK_SECONDARY)
    ax.set_xlim(2003, 2046)
    ax.set_ylabel("position along the orbit's line of nodes (mas)")
    ax.set_xlabel("year")
    ax.set_title("β Pictoris b: an orbit seen edge-on")
    ax.text(0.005, -0.2, f"{n_draws} orbits from the posterior (orange); measured positions 2003–2018 (dots, Nielsen et al. 2020). "
            f"P = {bp['P_yr'][0]:.1f} (+{bp['P_yr'][1]:.1f}/−{bp['P_yr'][2]:.1f}) yr, i = {bp['i_deg'][0]:.2f}°.",
            transform=ax.transAxes, fontsize=7.5, color=S.MUTED)

    # zoom on the 2017 passage: across vs along, the star and the Hill sphere
    ax2 = fig.add_subplot(gs[1, 0])
    tz = _mjd(np.linspace(2016.6, 2018.9, 600))
    for k in pick[:60]:
        x, y, z = K.orbit_xyz(tz, s["a"][k], s["e"][k], s["inc"][k], s["aop"][k], s["pan"][k], s["tau"][k], s["m"][k])
        al, pe = rotate(x * s["plx"][k], y * s["plx"][k])
        ax2.plot(al, pe, color=S.ORANGE, lw=0.6, alpha=0.15)
    rh = bp["hill_radius_mas"][0]
    tm = bp["conj2017_mjd"][0]
    x, y, z = K.orbit_xyz(tm, np.median(s["a"]), np.median(s["e"]), np.median(s["inc"]), np.median(s["aop"]),
                          np.median(s["pan"]), np.median(s["tau"]), np.median(s["m"]))
    al0, pe0 = rotate(x * np.median(s["plx"]), y * np.median(s["plx"]))
    ax2.add_patch(Circle((al0, pe0), rh, facecolor=S.ORANGE, alpha=0.08, edgecolor=S.ORANGE, lw=0.8))
    ax2.add_patch(Circle((0, 0), B.STAR_DIAMETER_MAS / 2 * 6, color=S.INK, zorder=6))
    ax2.annotate("star (disc ×6)", (0, 0), (-95, 35), fontsize=7.5, color=S.INK_SECONDARY,
                 arrowprops=dict(arrowstyle="-", color=S.INK_SECONDARY, lw=0.8))
    ax2.annotate(f"Hill sphere at closest approach\n(r = {rh:.0f} mas = {bp['hill_radius_au'][0]:.2f} au)",
                 (al0 + rh * 0.62, pe0 + rh * 0.78), (al0 + rh * 1.05, pe0 + rh * 1.15), fontsize=7.5,
                 color=S.INK_SECONDARY, arrowprops=dict(arrowstyle="-", color=S.INK_SECONDARY, lw=0.8))
    ax2.text(-125, -9.8 - 6, "the planet's track: it passes 9.8 mas from the star", fontsize=7,
             color=S.INK_SECONDARY, va="top")
    ax2.set_aspect("equal")
    ax2.set_xlim(-130, 130)
    ax2.set_ylim(-80, 80)
    ax2.set_xlabel("along the node line (mas)")
    ax2.set_ylabel("across (mas)")
    ax2.set_title("The 2017 passage, seen from Earth", fontsize=10)

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis("off")
    lines = [
        ("Closest approach", f"{bp['conj2017_iso']}  (± {bp['conj2017_mjd'][1]:.0f} d)", "Nielsen+20: 2017 Sep 13 ± 2.8 d"),
        ("Projected miss distance", f"{bp['conj2017_sep_mas'][0]:.1f} ± {bp['conj2017_sep_mas'][1]:.1f} mas = {bp['conj2017_sep_au'][0]:.2f} au",
         "Wang+16: 9.9 (+0.9/−0.8) mas"),
        ("", f"= {bp['conj2017_sep_stellar_radii'][0]:.0f} stellar radii: the planet itself never transits", ""),
        ("Hill sphere over the star", f"{bp['hill_ingress_iso']} → {bp['hill_egress_iso']}", "Nielsen+20: Apr 11 → Feb 16"),
        ("Passage in front (RV)", f"probability {bp['p_branch_front_2017']:.3f}", "CRIRES: approaching in 2013"),
        ("Next passages", f"behind {bp['behind_year'][0]:.1f}, in front {bp['next_front_year'][0]:.1f} (± {bp['next_front_year'][1]:.1f})", ""),
    ]
    y = 0.97
    for a, b, c_ in lines:
        if a:
            ax3.text(0.0, y, a, fontsize=8, color=S.INK_SECONDARY, va="top", transform=ax3.transAxes)
            y -= 0.075
        else:
            y -= 0.015
        ax3.text(0.0, y, b, fontsize=8.5, color=S.INK, va="top", transform=ax3.transAxes)
        y -= 0.07
        if c_:
            ax3.text(0.0, y, c_, fontsize=7.5, color=S.MUTED, va="top", transform=ax3.transAxes)
            y -= 0.085
        else:
            y -= 0.045
    S.caption(fig, "Fit: our nested sampling of orbitize!'s example data (Nielsen et al. 2020, Table 1); the CRIRES radial "
              "velocity (Snellen et al. 2014) decides which passage is in front. Hill radius for 12.8 MJup (Nielsen et al. 2020).", y=0.0)
    return S.save(fig, C.FIGURES / "betapic_edge_on")


# ---------------------------------------------------------------------------
# Reproduction of Wang et al. (2018)
# ---------------------------------------------------------------------------

def reproduction_figure(res: dict):
    """Forest plot: our posterior medians and 68% intervals vs Wang et al.'s Table 2."""
    models = list(C.MODELS)
    rows = []
    for p in PL:
        for par, lab in (("a", "a (au)"), ("e", "e"), ("i", "i (°)"), ("Omega", "Ω (°)")):
            rows.append((p, par, lab))
    fig, axes = plt.subplots(len(rows) // 4, 4, figsize=(11, 9.5), sharex=False)
    for r, (p, par, lab) in enumerate(rows):
        ax = axes[r // 4, r % 4]
        name = f"{par}_{p}"
        for k, kind in enumerate(models):
            o = res[kind]["posterior"][name]
            w = res[kind]["wang"].get(name)
            y = len(models) - k
            ax.errorbar(o["median"], y + 0.15, xerr=[[o["minus"]], [o["plus"]]], fmt="o", ms=4.5,
                        color=S.PLANET_COLOR[p], mec=S.SURFACE, mew=1, ecolor=S.PLANET_COLOR[p], elinewidth=1.6)
            if w:
                ax.errorbar(w[0], y - 0.15, xerr=[[w[2]], [w[1]]], fmt="s", ms=4, color=S.C_PUBLISHED,
                            mec=S.SURFACE, mew=1, ecolor=S.C_PUBLISHED, elinewidth=1.2)
        ax.set_yticks(range(1, len(models) + 1))
        ax.set_yticklabels([C.MODEL_LABEL[m].replace(", low e", "\n+ low e") for m in models[::-1]] if r % 4 == 0 else [],
                           fontsize=7)
        ax.set_ylim(0.4, len(models) + 0.6)
        ax.set_title(f"{p}: {lab}", fontsize=9)
        ax.grid(axis="y", visible=False)
    from matplotlib.lines import Line2D
    fig.legend([Line2D([], [], marker="o", color=S.INK, ls="", ms=5),
                Line2D([], [], marker="s", color=S.C_PUBLISHED, ls="", ms=5)],
               ["this work (nested sampling; planet colour)", "Wang et al. 2018, Table 2"],
               loc="upper left", ncol=2, bbox_to_anchor=(0.015, 0.968), fontsize=8.5)
    fig.suptitle("Reproducing Wang et al. (2018): four orbit models, same data, same priors", x=0.02, ha="left",
                 y=0.995)
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    S.caption(fig, "Medians and 16-84% intervals. Ω folded into 0-180° as in Wang et al. Data: Keck/NIRC2 2004-2014 "
              "(Konopacky et al. 2016) + GPI 2013-2016 (Wang et al. 2018). Wang's a_e for the two 1:2:4:8 models is "
              "printed as ±0.06 au, apparently a decimal slip (±0.6); drawn as printed.", y=0.0)
    return S.save(fig, C.FIGURES / "wang2018_reproduction")


def coplanarity_figure(fit_unc, fit_cop):
    """Wang et al.'s Fig. 2: Omega vs i of each planet, unconstrained fit, with the coplanar plane."""
    t = fit_unc.table()
    tc = fit_cop.table()
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.4), sharey=True)
    for ax, p in zip(axes, IO):
        O = K.wrap180(t[f"Omega_{p}"].to_numpy())
        i = t[f"i_{p}"].to_numpy()
        H, xe, ye = np.histogram2d(O, i, bins=[np.linspace(0, 180, 46), np.linspace(0, 70, 36)])
        ax.pcolormesh(xe, ye, H.T, cmap=_seq_cmap(S.PLANET_COLOR[p]), shading="auto")
        Oc, ic = K.wrap180(tc["Omega_e"].to_numpy()), tc["i_e"].to_numpy()
        ax.scatter(Oc[::20], ic[::20], s=2, color=S.INK, alpha=0.15, lw=0)
        ax.set_title(f"planet {p}", fontsize=9)
        ax.set_xlabel("Ω (°, folded)")
        ax.set_xlim(0, 180)
        ax.set_ylim(0, 70)
    axes[0].set_ylabel("inclination i (°)")
    fig.suptitle("Is the system flat? Each planet's orbital plane fitted on its own", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    S.caption(fig, "Coloured: posterior density of each planet's plane in the Unconstrained fit. Grey dots: the single "
              "plane of the Coplanar fit. Where the colour overlaps the grey, one flat system fits that planet too.", y=0.02)
    return S.save(fig, C.FIGURES / "coplanarity")


def _seq_cmap(color):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("seq", [S.SURFACE, color])



# ---------------------------------------------------------------------------
# Stability
# ---------------------------------------------------------------------------

#: model colours where models (not planets) are the series: slots 5, 6, 8 and ink
MODEL_COLOR = {"unconstrained": S.MAGENTA, "coplanar": S.GREEN, "near1248": S.RED, "near1248_coplanar": S.INK}


def stability_figure(summ: dict, wang: dict, stable_masses: pd.DataFrame, hot: pd.DataFrame | None = None,
                     extra: dict | None = None, all_masses: pd.DataFrame | None = None):
    """summ: model -> analysis.stability_summary (W18 posteriors); extra: fit key ->
    summary for the all-data fits (no published counterpart)."""
    import stability as St
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3), gridspec_kw=dict(width_ratios=[1.1, 1.2, 1.05]))
    # A: stable per million, ours vs Wang
    ax = axes[0]
    models = [m for m in C.MODELS if m in summ]
    rows = [(C.MODEL_LABEL[m].replace(", low e", "\n+ low e"), summ[m], wang[m]["stable_per_million"])
            for m in models]
    for key, sm in (extra or {}).items():
        kind = key.split("/")[1]
        rows.append(("all data: " + C.MODEL_LABEL[kind].replace(", low e", "\n+ low e"), sm, None))
    n = len(rows)
    for k, (lab, sm, w) in enumerate(rows):
        s = sm["stable"]
        y = n - k
        if s["k"] > 0:
            ax.errorbar(s["per_million"], y + 0.12, xerr=[[s["per_million"] - s["lo"]], [s["hi"] - s["per_million"]]],
                        fmt="o", color=S.INK, ms=5, mec=S.SURFACE, mew=1, elinewidth=1.4)
        else:
            ax.errorbar(s["upper95"], y + 0.12, xerr=[[s["upper95"] * 0.5], [0]], xuplims=True, fmt="o",
                        color=S.INK, ms=5, mec=S.SURFACE, mew=1, elinewidth=1.4)
        if w is not None:
            ax.scatter(max(w, 0.15), y - 0.12, marker="s", s=28, color=S.C_PUBLISHED, edgecolors=S.SURFACE, lw=1)
        ax.text(1.4e5, y, f"{s['k']} / {s['n']:,}", va="center", fontsize=7.5, color=S.INK_SECONDARY)
    if extra:
        ax.axhline(n - len(models) + 0.5, color=S.AXIS, lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(0.1, 3e6)
    ax.set_xticks([1, 10, 100, 1e3, 1e4, 1e5])
    ax.set_yticks(range(1, n + 1))
    ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=7)
    ax.set_ylim(0.4, n + 0.6)
    ax.set_xlabel("orbits stable for 40 Myr, per million draws")
    ax.set_title("Which orbits survive?", fontsize=10)
    ax.grid(axis="y", visible=False)
    from matplotlib.lines import Line2D
    ax.legend([Line2D([], [], marker="o", color=S.INK, ls=""), Line2D([], [], marker="s", color=S.C_PUBLISHED, ls="")],
              ["this work (Wang's set-up)", "Wang et al. 2018 (their 0 drawn at 0.15)"], loc="lower left", fontsize=7)
    # B: survival-time distribution (solid: Wang's data; dashed: all data)
    ax = axes[1]
    curves = [(f"W18/{m}", m, "-", C.MODEL_LABEL[m].replace(" coplanar, low e", " coplanar")) for m in models]
    curves += [(k, k.split("/")[1], "--", "all data: " + C.MODEL_LABEL[k.split("/")[1]].replace(" coplanar, low e",
                                                                                                    " coplanar"))
               for k in (extra or {})]
    for key, m, ls, lab in curves:
        df = pd.read_csv(St.results_path(key))
        n_all = len(df)
        t = np.sort(df.t_end_yr[df.outcome != "stable"].to_numpy())
        frac = 1 - np.arange(1, len(t) + 1) / n_all            # still intact after each break-up
        t, frac = np.append(t, 4e7), np.append(frac, frac[-1])  # survivors: flat to 40 Myr
        ax.plot(t, np.maximum(frac, 1e-9), color=MODEL_COLOR[m], lw=1.8, ls=ls)
        if ls == "-":                                     # dashed (all-data) curves: legend below
            j = np.searchsorted(t, 3e5)
            ax.text(3.2e5, frac[min(j, len(frac) - 1)] * 1.15, lab, fontsize=7, color=S.INK_SECONDARY)
    if extra:
        from matplotlib.lines import Line2D
        ax.legend([Line2D([], [], color=S.INK, lw=1.8), Line2D([], [], color=S.INK, lw=1.8, ls="--")],
                  ["2004–2016 data (Wang et al.'s)", "1998–2025 data, same model colours"],
                  loc="lower left", fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e2, 4.5e7)
    ax.set_ylim(3e-5, 1.1)
    ax.axvline(4e7, color=S.AXIS, lw=1)
    ax.set_xlabel("years survived (integrating backwards from 2013)")
    ax.set_ylabel("fraction of draws still intact")
    ax.set_title("Most configurations break up quickly", fontsize=10)
    # C: stable masses
    ax = axes[2]
    ax.fill_between([4, 11], [3, 3], [4, 11], color=S.CONTEXT, alpha=0.35, lw=0)
    ax.text(10.8, 3.3, "Wang's mass prior", ha="right", fontsize=7, color=S.MUTED)
    if all_masses is not None and len(all_masses):
        ax.scatter(all_masses.m_cde, all_masses.m_b, s=14, facecolors="none", edgecolors=S.MUTED, lw=0.7,
                   label=f"stable, all-data orbits ({len(all_masses)})", zorder=2)
    ax.scatter(stable_masses.m_cde, stable_masses.m_b, s=12, color=S.INK, edgecolors=S.SURFACE, lw=0.5,
               label=f"stable, Wang's data ({len(stable_masses)})", zorder=3)
    if hot is not None and len(hot):
        ax.scatter(hot.m_cde, hot.m_b, s=14, color=S.ORANGE, edgecolors=S.SURFACE, lw=0.5,
                   label=f"stable, hot-start masses ({len(hot)})", zorder=4)
    ax.set_xlim(4, 11)
    ax.set_ylim(3, 9)
    ax.set_xlabel("mass of c, d and e (MJup, each)")
    ax.set_ylabel("mass of b (MJup)")
    ax.set_title("Survivors are light", fontsize=10)
    ax.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    S.caption(fig, "N-body simulations (REBOUND WHFast): posterior orbits run 40 Myr back in time with Wang et al.'s "
              "set-up and instability criteria. Arrows: 95% upper limits where nothing survived.", y=-0.02)
    return S.save(fig, C.FIGURES / "stability")


def sky_orbits_figure(sets: list, data: pd.DataFrame, stem="sky_orbits", n=150, title=None):
    """Small multiples of orbits on the sky: sets = [(label, model, samples)]."""
    import predict as Pr
    fig, axes = plt.subplots(1, len(sets), figsize=(4.2 * len(sets), 4.6), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(C.RANDOM_SEED)
    for ax, (label, model, samples) in zip(axes, sets):
        s = samples[rng.choice(len(samples), min(n, len(samples)), replace=False)]
        phys = model.physical(s)
        for p in C.PLANETS:
            a, e, inc, aop, pan, tau = phys["elements"][p]
            P = np.median(K.period_days(a, phys["mstar"]))
            t = _mjd(2000.0) + np.linspace(0, 1.02 * P, 500)
            sky = Pr.predicted_positions(model, s, t)
            ra, dec = sky[p]
            for k in range(ra.shape[0]):
                ax.plot(ra[k], dec[k], color=S.PLANET_COLOR[p], lw=0.5, alpha=0.06)
        _ring(ax, data.ra, data.dec, size=6)
        ax.scatter([0], [0], marker="*", s=90, color=S.INK, zorder=6)
        ax.set_aspect("equal")
        ax.set_xlim(2100, -2100)
        ax.set_ylim(-2100, 2100)
        ax.set_title(label, fontsize=9.5)
        ax.set_xlabel("ΔRA (mas)")
        for p in C.PLANETS:          # letter just outside the middle of the planet's measured arc
            dp = data[data.planet == p]
            x0, y0 = np.median(dp.ra), np.median(dp.dec)
            r0 = np.hypot(x0, y0)
            ax.text(x0 * (1 + 120 / r0), y0 * (1 + 120 / r0), p, fontsize=10, fontweight="bold", color=S.INK,
                    ha="center", va="center")
    axes[0].set_ylabel("ΔDec (mas)")
    if title:
        fig.suptitle(title, x=0.02, ha="left")
    fig.tight_layout()
    return S.save(fig, C.FIGURES / stem)


# ---------------------------------------------------------------------------
# Prediction test, model evidence, GPI calibration
# ---------------------------------------------------------------------------

SOURCE_MARKER = {"IRDIS-Z22": "o", "IFS-Z22": "D", "LUCI-Z22": "^", "GRAVITY-L19": "*",
                 "NIRCam-B25": "s", "AMI-N26": "v", "NIRC2-N26": "P"}
SOURCE_LABEL = {"IRDIS-Z22": "SPHERE/IRDIS", "IFS-Z22": "SPHERE/IFS", "LUCI-Z22": "LBT/LUCI",
                "GRAVITY-L19": "VLTI/GRAVITY", "NIRCam-B25": "JWST/NIRCam", "AMI-N26": "JWST/AMI",
                "NIRC2-N26": "Keck/NIRC2 2025"}


def prediction_figure(tab: pd.DataFrame, fits_shown=("W18/near1248_coplanar", "W18R/near1248_coplanar")):
    """Miss distance of each later measurement vs time, per planet, for two fits."""
    fig, axes = plt.subplots(len(fits_shown), 4, figsize=(12, 3.0 * len(fits_shown)), sharex=True, sharey=True)
    labels = {"W18/near1248_coplanar": "fitted to 2004-2016, GPI as published (Wang et al.)",
              "W18R/near1248_coplanar": "same data, GPI with De Rosa et al.'s 2020 calibration",
              "W18/near1248_coplanar (stable)": "only the orbits stable for 40 Myr"}
    for r, key in enumerate(fits_shown):
        t = tab.loc[key] if key in tab.index.get_level_values(0) else None
        for c, p in enumerate(IO):
            ax = axes[r, c]
            d = t[t.planet == p]
            for src, g in d.groupby("source"):
                grav = src == "GRAVITY-L19"
                ax.errorbar(g.year, g.res_mas, yerr=np.hypot(g.sigma_obs_mas, g.sigma_pred_mas), fmt=SOURCE_MARKER[src],
                            ms=12 if grav else 5, color=S.PLANET_COLOR[p], mec=S.INK if grav else S.SURFACE,
                            mew=1.0 if grav else 0.8, ecolor=S.AXIS, elinewidth=0.8, label=SOURCE_LABEL[src],
                            zorder=6 if grav else 3)
            ax.axhline(0, color=S.AXIS, lw=0.8)
            if r == 0:
                ax.set_title(f"planet {p}", fontsize=9.5)
            if c == 0:
                ax.set_ylabel("miss distance (mas)")
            if r == len(fits_shown) - 1:
                ax.set_xlabel("year")
            ax.set_ylim(0, 72)
            ax.set_xlim(2016.5, 2026.3)
    handles, lab = [], []
    from matplotlib.lines import Line2D
    for src in SOURCE_MARKER:
        handles.append(Line2D([], [], marker=SOURCE_MARKER[src], color=S.INK, ls="", ms=6))
        lab.append(SOURCE_LABEL[src])
    fig.legend(handles, lab, loc="upper center", ncol=7, fontsize=7.5, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Predicting the planets: orbits fitted up to Sep 2016 against every later measurement",
                 x=0.02, ha="left", y=1.08)
    fig.tight_layout(h_pad=3.0)
    for r, key in enumerate(fits_shown):        # row labels above each row, clear of the panel titles
        bb = axes[r, 0].get_position()
        fig.text(bb.x0, bb.y1 + (0.045 if r == 0 else 0.012), labels.get(key, key), fontsize=9,
                 fontweight="bold", color=S.INK_SECONDARY, va="bottom")
    S.caption(fig, "Miss = distance between the measured position and the median prediction; error bars combine the "
              "measurement (with each instrument's fitted extra scatter) and the spread of the predictions.", y=-0.03)
    return S.save(fig, C.FIGURES / "prediction_test")


def evidence_figure(lnz: dict):
    """lnz: data set -> model -> (lnZ, err). Bayes factors relative to Unconstrained."""
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    sets = list(lnz)
    for k, ds in enumerate(sets):
        base = lnz[ds]["unconstrained"][0]
        for j, m in enumerate(C.MODELS):
            v, e = lnz[ds][m]
            y = len(C.MODELS) - j + (0.15 if k == 0 else -0.15)
            ax.errorbar(v - base, y, xerr=e, fmt="o" if k == 0 else "s", color=S.MUTED if k == 0 else S.INK,
                        ms=6, mec=S.SURFACE, mew=1, label=ds if j == 0 else None)
            ax.text(v - base + 0.6, y, f"{v - base:+.1f}", va="center", fontsize=7.5, color=S.INK_SECONDARY)
    ax.axvline(0, color=S.AXIS, lw=0.8)
    ax.set_yticks(range(1, len(C.MODELS) + 1))
    ax.set_yticklabels([C.MODEL_LABEL[m].replace(", low e", "\n+ low e") for m in C.MODELS[::-1]], fontsize=8)
    ax.set_xlabel("ln Bayes factor against Unconstrained (nested sampling)")
    ax.set_title("Which orbital architecture do the data prefer?", fontsize=10)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="upper right", fontsize=7.5)
    fig.tight_layout()
    S.caption(fig, "A difference of 5 in ln Z is odds of ~150:1. These Bayes factors use Wang et al.'s priors: they "
              "measure how much prior volume each model wastes on orbits the data exclude.", y=-0.03)
    return S.save(fig, C.FIGURES / "model_evidence")


def _offset_draws(fit, n=600, seed=C.RANDOM_SEED):
    """Posterior draws of each instrument's (scale, rotation) - conditional Gaussians."""
    rng = np.random.default_rng(seed)
    m = fit.model
    out = {s: [] for s in m.data.offset_sources}
    for th in fit.equal[rng.choice(len(fit.equal), n, replace=False)]:
        _, offs = m.loglike_phys(m.physical(th), return_offsets=True)
        for s, (sol, cov) in offs.items():
            out[s].append(rng.multivariate_normal(sol, cov))
    return {s: np.array(v) for s, v in out.items()}


def calibration_figure(fit_rev, fit_orig, derosa_change=0.36, derosa_err=0.12):
    """What the planets say about each instrument's north angle and plate scale."""
    o_rev, o_orig = _offset_draws(fit_rev), _offset_draws(fit_orig)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.9), gridspec_kw=dict(width_ratios=[1.1, 1.3]))
    ax = axes[0]
    corr_orig = -np.degrees(o_orig["GPI-W18"][:, 1])      # correction to add to the published PAs
    corr_rev = -np.degrees(o_rev["GPI-D20"][:, 1])
    bins = np.linspace(-0.35, 0.8, 47)
    ax.hist(corr_orig, bins, color=S.INK, alpha=0.85, label="GPI as published (2018)")
    ax.hist(corr_rev, bins, color=S.MUTED, alpha=0.6, label="GPI recalibrated (De Rosa+2020)")
    ax.axvspan(derosa_change - derosa_err, derosa_change + derosa_err, color=S.C_PUBLISHED, alpha=0.15, lw=0)
    ax.axvline(derosa_change, color=S.C_PUBLISHED, lw=1.5)
    ax.text(derosa_change + 0.02, ax.get_ylim()[1] * 0.92, "De Rosa et al.'s\ncorrection\n(binary stars)",
            fontsize=7.5, color=S.INK_SECONDARY, va="top")
    ax.axvline(0, color=S.AXIS, lw=0.8)
    ax.set_xlabel("north-angle correction the planets ask for (deg)")
    ax.set_ylabel("posterior draws")
    ax.set_title("HR 8799 itself re-derives GPI's calibration error", fontsize=10)
    ax.legend(loc="upper left", fontsize=7.5)
    ax.grid(axis="y", visible=False)
    # B: every instrument
    ax = axes[1]
    names = [s for s in ("NIRC2-K16", "GPI-D20", "IRDIS-Z22", "IFS-Z22", "LUCI-Z22", "NIRCam-B25", "NIRC2-N26")
             if s in o_rev]
    lab = {"NIRC2-K16": "Keck/NIRC2 2004-14", "GPI-D20": "Gemini/GPI (revised)", "IRDIS-Z22": "SPHERE/IRDIS",
           "IFS-Z22": "SPHERE/IFS", "LUCI-Z22": "LBT/LUCI", "NIRCam-B25": "JWST/NIRCam", "NIRC2-N26": "Keck/NIRC2 2025"}
    for k, s in enumerate(names):
        rot = np.degrees(o_rev[s][:, 1])
        y = len(names) - k
        lo, med, hi = np.percentile(rot, [16, 50, 84])
        ax.errorbar(med, y, xerr=[[med - lo], [hi - med]], fmt="o", color=S.INK, ms=5, mec=S.SURFACE, mew=1)
        prior = C.CALIBRATION[s][0]
        ax.plot([-prior, prior], [y - 0.25, y - 0.25], color=S.CONTEXT, lw=3, solid_capstyle="butt")
    ax.set_yticks(range(1, len(names) + 1))
    ax.set_yticklabels([lab[s] for s in names[::-1]], fontsize=8)
    ax.axvline(0, color=S.AXIS, lw=0.8)
    ax.set_xlabel("fitted north-angle offset (deg; grey bar: prior from the instrument's own calibration)")
    ax.set_xlim(-0.6, 0.6)
    ax.set_title("All instruments, fitted together with the orbits", fontsize=10)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    S.caption(fig, "Coplanar fits to all data 1998-2025; the offsets enter the model linearly and are integrated out "
              "analytically, then drawn back for this figure. GRAVITY (interferometric) and HST carry no offsets.", y=-0.03)
    return S.save(fig, C.FIGURES / "instrument_calibration")


def hgca_figure(h: dict, Df, mtot):
    """The star's measured acceleration and what the four planets predict."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[1, 1.2]))
    ax = axes[0]
    d = np.array(h["measured"]["delta_mu_masyr"])
    cov = np.array(h["measured"]["cov"])
    w, v = np.linalg.eigh(cov)
    from matplotlib.patches import Ellipse
    for nsig, a in ((1, 0.25), (2, 0.12)):
        ax.add_patch(Ellipse(d, 2 * nsig * np.sqrt(w[1]), 2 * nsig * np.sqrt(w[0]),
                             angle=np.degrees(np.arctan2(v[1, 1], v[0, 1])), color=S.INK, alpha=a, lw=0))
    ax.scatter(*d, s=46, color=S.INK, zorder=6, edgecolors=S.SURFACE, linewidths=1.2)
    ax.text(d[0] - 0.035, d[1] - 0.205, "measured\n(1σ, 2σ)", fontsize=7.5, color=S.INK_SECONDARY, ha="center")
    # the four planets' pulls, tip to tail, at hot-start masses (median orbit)
    m = np.array([8.7, 8.7, 8.7, 6.0]) * C.MJUP_MSUN
    contrib = np.median(Df * m[None, :, None] / mtot[:, None, None], axis=0)
    start = np.zeros(2)
    for k, p in enumerate(IO):
        end = start + contrib[k]
        ax.annotate("", end, start, arrowprops=dict(arrowstyle="-|>", color=S.PLANET_COLOR[p], lw=2,
                                                    shrinkA=0, shrinkB=0))
        mid = 0.5 * (start + end)
        ax.text(mid[0] + (0.018 if p != "e" else -0.035), mid[1] + 0.01, p, fontsize=9, fontweight="bold",
                color=S.INK)
        start = end
    for label, key, mk, dxy in (("hot-start masses\n(b 6.0, c d e 8.7 MJup)", "hotstart_42Myr", "D", (0.03, 0.02)),
                                ("Wang+18 combined\n(b 5.8, c d e 7.2)", "wang_combined", "s", (0.03, 0.03))):
        pm = h["final_orbits"][key]["predicted_masyr"]
        ax.scatter(*pm, marker=mk, s=40, color=S.INK_SECONDARY, edgecolors=S.SURFACE, linewidths=1, zorder=6)
        ax.annotate(label, pm, (pm[0] + dxy[0] + 0.1, pm[1] + dxy[1] + (0.0 if key == "wang_combined" else -0.07)),
                    fontsize=7.5, color=S.INK_SECONDARY, va="center",
                    arrowprops=dict(arrowstyle="-", color=S.INK_SECONDARY, lw=0.8))
    ax.axhline(0, color=S.AXIS, lw=0.8)
    ax.axvline(0, color=S.AXIS, lw=0.8)
    ax.set_xlim(-0.5, 0.15)
    ax.set_ylim(-0.6, 0.15)
    ax.set_aspect("equal")
    ax.set_xlabel("Δμ, east (mas/yr)")
    ax.set_ylabel("Δμ, north (mas/yr)")
    ax.set_title("The star is accelerating", fontsize=10)
    # B: mass of e
    ax = axes[1]
    for key, label, col in (("brandt_reproduction", "Wang's stable orbits (Brandt+21 set-up)", S.MUTED),
                            ("final_orbits", "orbits fitted to all data 1998-2025", S.INK)):
        me = h[key]["mass_e"]
        x = np.linspace(0, 20, 400)
        sig = 0.5 * (me["plus"] + me["minus"])
        ax.plot(x, np.exp(-0.5 * ((x - me["median"]) / sig) ** 2), color=col, lw=2,
                label=f"{label}: {me['median']:.1f} ± {sig:.1f}")
    b = h["brandt_reproduction"]["brandt_mass_e"]
    ax.errorbar(b[0], 1.07, xerr=[[b[2]], [b[1]]], fmt="s", color=S.C_PUBLISHED, ms=5, label="Brandt et al. 2021: 9.6 +1.9/−1.8")
    ax.set_ylim(0, 1.45)
    ax.axvspan(8.7 - 1.0, 8.7 + 1.0, color=S.CONTEXT, alpha=0.4, lw=0)
    ax.text(8.7, 0.05, "hot-start\n42 Myr", ha="center", fontsize=7, color=S.MUTED)
    ax.axvline(13, color=S.AXIS, lw=1)
    ax.text(13.2, 0.45, "deuterium-burning\nlimit (~13 MJup)", fontsize=7, color=S.MUTED)
    ax.set_xlabel("mass of HR 8799 e (MJup), with m_c = m_d = m_e = 1.25 m_b")
    ax.set_yticks([])
    ax.legend(loc="upper right", fontsize=7, bbox_to_anchor=(1.0, 1.0))
    ax.set_xlim(0, 22)
    ax.set_title("Weighing planet e by the star's wobble", fontsize=10)
    ax.grid(visible=False)
    fig.tight_layout()
    S.caption(fig, "HGCA EDR3 (Brandt 2021, CDS J/ApJS/254/42). Left: the measured change in the star's proper motion "
              "and the four planets' contributions at hot-start masses. Right: the mass of e if the others scale with it.",
              y=-0.03)
    return S.save(fig, C.FIGURES / "hgca_wobble")


def keeper_figure(fit, data: pd.DataFrame, today_mjd: float, n=120):
    """The orbits fitted to 1998-2025, every measurement, and where the planets are today."""
    import predict as Pr
    model, samples = fit.model, fit.equal
    t = fit.table()
    rng = np.random.default_rng(C.RANDOM_SEED)
    s = samples[rng.choice(len(samples), n, replace=False)]
    fig = plt.figure(figsize=(11.5, 8.2))
    ax = fig.add_axes([0.04, 0.06, 0.62, 0.86])
    # faint posterior orbits
    ext = []
    for p in C.PLANETS:
        P = float(np.median(t[f"P_{p}"]))
        tt = _mjd(1990.0) + np.linspace(0, 1.01 * P * 365.25, 600)
        ra, dec = Pr.predicted_positions(model, s, tt)[p]
        for k in range(ra.shape[0]):
            ax.plot(ra[k], dec[k], color=S.PLANET_COLOR[p], lw=0.5, alpha=0.05)
        ext.append([ra.min(), ra.max(), dec.min(), dec.max()])
    ext = np.array(ext)
    xlo, xhi, ylo, yhi = ext[:, 0].min(), ext[:, 1].max(), ext[:, 2].min(), ext[:, 3].max()
    pad = 0.04 * max(xhi - xlo, yhi - ylo)
    xlo, xhi, ylo, yhi = xlo - pad, xhi + pad, ylo - pad, yhi + pad
    # yearly ticks 1990-2050 on the median track, labelled every decade
    yrs = np.arange(1990, 2051, 1.0)
    trk = Pr.predicted_positions(model, s, _mjd(yrs))
    for p in C.PLANETS:
        ra, dec = np.median(trk[p][0], 0), np.median(trk[p][1], 0)
        ax.plot(ra, dec, color=S.PLANET_COLOR[p], lw=2.2, solid_capstyle="round")
        ax.scatter(ra, dec, s=6, color=S.PLANET_COLOR[p], zorder=3)
        for y in (2000, 2050):
            k = int(np.where(yrs == y)[0][0])
            r = np.hypot(ra[k], dec[k])
            ax.text(ra[k] * (1 + 95 / r), dec[k] * (1 + 95 / r), str(y), fontsize=7, color=S.INK_SECONDARY,
                    ha="center", va="center")
    _ring(ax, data.ra, data.dec, size=9)
    now = Pr.predicted_positions(model, s, np.array([today_mjd]))
    for p in C.PLANETS:
        x, y = np.median(now[p][0][:, 0]), np.median(now[p][1][:, 0])
        ax.scatter(x, y, s=120, color=S.PLANET_COLOR[p], edgecolors=S.SURFACE, linewidths=2, zorder=7)
        r = np.hypot(x, y)
        ax.text(x * (1 - 120 / r), y * (1 - 120 / r), p, fontsize=13, fontweight="bold", color=S.INK,
                ha="center", va="center", zorder=8)
    ax.scatter([0], [0], marker="*", s=260, color=S.INK, zorder=8, edgecolors=S.SURFACE, linewidths=1)
    plx = float(np.median(t.plx))
    x0, y0 = xlo + 0.04 * (xhi - xlo), ylo + 0.04 * (yhi - ylo)
    ax.plot([x0, x0 + 20 * plx], [y0, y0], color=S.INK_SECONDARY, lw=2.2)
    ax.text(x0 + 10 * plx, y0 + 0.015 * (yhi - ylo), "20 au", ha="center", va="bottom", fontsize=8,
            color=S.INK_SECONDARY)
    ax.set_aspect("equal")
    ax.set_xlim(xhi, xlo)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel("ΔRA (mas, east to the left)")
    ax.set_ylabel("ΔDec (mas)")
    ax.set_title("HR 8799: four giant planets on their orbits", fontsize=13)
    # side table
    tx = fig.add_axes([0.69, 0.06, 0.30, 0.86])
    tx.axis("off")
    y = 0.98
    tx.text(0, y, "Fitted to all data 1998-2025", fontsize=10, fontweight="bold", va="top")
    y -= 0.045
    tx.text(0, y, "(near 1:2:4:8 coplanar model; medians ± 1σ)", fontsize=8, color=S.INK_SECONDARY, va="top")
    y -= 0.07

    def pm(name, fmt="{:.1f}"):
        lo, med, hi = np.percentile(t[name], [16, 50, 84])
        return f"{fmt.format(med)} ± {fmt.format(0.5 * (hi - lo))}"

    for p in C.PLANETS:
        tx.scatter([0.01], [y - 0.012], s=60, color=S.PLANET_COLOR[p], clip_on=False)
        tx.text(0.06, y, f"planet {p}", fontsize=10, fontweight="bold", va="top")
        y -= 0.04
        tx.text(0.06, y, f"a = {pm(f'a_{p}')} au   P = {pm(f'P_{p}', '{:.0f}')} yr", fontsize=8.5, va="top")
        y -= 0.035
        tx.text(0.06, y, f"e = {pm(f'e_{p}', '{:.2f}')}", fontsize=8.5, va="top")
        y -= 0.06
    tx.text(0, y, f"shared plane: i = {pm('i_e')}°, Ω = {pm('Omega_e')}°", fontsize=8.5, va="top")
    y -= 0.04
    tx.text(0, y, f"period ratios {pm('P_b/P_c', '{:.2f}')} : {pm('P_c/P_d', '{:.2f}')} : {pm('P_d/P_e', '{:.2f}')}",
            fontsize=8.5, va="top")
    y -= 0.04
    tx.text(0, y, f"star: {pm('mstar', '{:.2f}')} Msun (from the orbits)", fontsize=8.5, va="top")
    y -= 0.07
    tx.text(0, y, f"Ink dots: {len(data)} measured positions (HST 1998,\nKeck, GPI, SPHERE, LBT, GRAVITY, JWST).\n"
            "Coloured dots: positions on 2026 Oct 8.\nTicks: one per year 1990-2050 (labels 2000, 2050).",
            fontsize=8, color=S.INK_SECONDARY, va="top")
    tx.set_xlim(0, 1)
    tx.set_ylim(0, 1)
    return S.save(fig, C.FIGURES / "hr8799_orbits")


# ---------------------------------------------------------------------------
# Resonances
# ---------------------------------------------------------------------------

ANGLE_LABEL = {"b:c,b": "θ b:c (ϖ_b)", "b:c,c": "θ b:c (ϖ_c)", "c:d,c": "θ c:d (ϖ_c)", "c:d,d": "θ c:d (ϖ_d)",
               "d:e,d": "θ d:e (ϖ_d)", "d:e,e": "θ d:e (ϖ_e)", "c:d:e": "Laplace c:d:e", "b:c:d:e": "four-body b:c:d:e"}


def resonance_figure(df: pd.DataFrame, recs: list):
    """Top: the resonant angles of one stable system over 40 Myr. Below: for every
    stable system, the fraction of 100-kyr windows in which each angle librates."""
    import stability as St
    st = df[df.outcome == "stable"].reset_index(drop=True)
    locked = st[st["lib_d:e,e"] >= 0.9]
    ex = (locked if len(locked) else st).sort_values("m_cde").iloc[-1]
    rec = next(r for r in recs if int(r["draw"]) == int(ex["draw"]))
    fig = plt.figure(figsize=(12, 9.2))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1], hspace=0.62, wspace=0.12)
    show = ["d:e,e", "c:d,d", "c:d:e", "b:c,b"]
    for c, name in enumerate(show):
        ax = fig.add_subplot(gs[0, c])
        th = np.degrees(rec["_series"][name].astype(float))
        th = np.where(th > 180, th - 360, th)
        t = np.arange(1, len(th) + 1) * 10 * C.SNAPSHOT_YEARS / 1e6
        ax.scatter(t, th, s=0.6, color=S.INK, alpha=0.35, lw=0, rasterized=True)
        ax.set_ylim(-180, 180)
        ax.set_yticks([-180, -90, 0, 90, 180])
        if c:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel("angle (°)")
        ax.set_xlim(0, 40)
        ax.set_xlabel("Myr into the past")
        f = ex[f"lib_{name}"]
        ax.set_title(f"{ANGLE_LABEL[name]}\nlibrates {100 * f:.0f}% of the time", fontsize=8.5)
    fig.text(0.125, 0.945, f"One stable system (M_cde = {ex.m_cde:.1f}, M_b = {ex.m_b:.1f} MJup): a librating angle stays "
             "in a band, a circulating one fills every value", fontsize=9.5, fontweight="bold", color=S.INK_SECONDARY)
    for k, name in enumerate(St.ANGLES):
        ax = fig.add_subplot(gs[1 + k // 4, k % 4])
        ax.scatter(st.m_cde, 100 * st[f"lib_{name}"], s=10, color=S.INK, edgecolors=S.SURFACE, lw=0.4, alpha=0.8)
        ax.scatter([ex.m_cde], [100 * ex[f"lib_{name}"]], s=60, facecolors="none", edgecolors=S.INK, lw=1.2)
        ax.set_ylim(-5, 105)
        ax.set_xlim(st.m_cde.min() - 0.2, st.m_cde.max() + 0.2)
        ax.set_title(ANGLE_LABEL[name], fontsize=8.5)
        if k % 4:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel("% of time librating")
        if k // 4 == 1:
            ax.set_xlabel("mass of c, d, e (MJup)")
    fig.text(0.125, 0.627, f"All {len(st)} stable systems: how much of the 40 Myr each angle spends librating "
             "(ring: the system above)", fontsize=9.5, fontweight="bold", color=S.INK_SECONDARY)
    fig.suptitle("Are the planets locked in resonance?", x=0.125, ha="left", y=0.985, fontsize=13)
    S.caption(fig, "θ for the pair i:j (j inside) = λ_j − 2λ_i + ϖ_k; Laplace = λ_e − 3λ_d + 2λ_c; four-body = "
              "λ_e − 2λ_d − λ_c + 2λ_b. Librating: within a 100-kyr window the angle never visits some 30° stretch of the circle. "
              "Angles stored every 2,000 yr (the streaks are stroboscopic). Orbits: 2004-2016 data, near 1:2:4:8 coplanar fit, "
              "Wang et al.'s set-up.", y=0.02)
    return S.save(fig, C.FIGURES / "resonances")


def future_figure(fut: pd.DataFrame, twins: pd.DataFrame, years: float = 100e6):
    """Each stable system run forwards from today (up to 100 Myr), against its mass;
    filled if its nudged twin also survived the 40 Myr into the past."""
    d = fut.merge(twins[["draw", "outcome"]].rename(columns={"outcome": "twin"}), on="draw", how="left")
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    t = np.where(d.outcome == "stable", years, d.t_end_yr) / 1e6
    same = d.twin == "stable"
    ax.scatter(d.m_cde[same], t[same], s=26, color=S.INK, edgecolors=S.SURFACE, lw=0.6, zorder=3,
               label=f"twin also stable 40 Myr ({int(same.sum())})")
    ax.scatter(d.m_cde[~same], t[~same], s=26, facecolors=S.SURFACE, edgecolors=S.INK, lw=1.0, zorder=3,
               label=f"twin broke up ({int((~same).sum())})")
    ax.axhline(40, color=S.AXIS, lw=1, ls=":")
    ax.set_yscale("log")
    ax.set_ylim(max(0.5 * t.min(), 0.01), years / 1e6 * 1.4)
    ax.set_xlabel("mass of c, d and e (MJup, each)")
    ax.set_ylabel("survives into the future (Myr)")
    n_surv = int((d.outcome == "stable").sum())
    ax.set_title(f"After today: {n_surv} of {len(d)} stable systems last another {years / 1e6:.0f} Myr", fontsize=10)
    ax.legend(loc="lower left", fontsize=7.5)
    fig.tight_layout()
    S.caption(fig, "Every system stable for 40 Myr in the past, integrated forwards from 2013 (top row: still intact "
              "at the end, 100 Myr; dotted: 40 Myr, about the system's age). Its twin is the same system with planet b "
              "moved by 10⁻¹² au (0.15 m), run 40 Myr into the past again.", y=-0.03)
    return S.save(fig, C.FIGURES / "future_stability")
