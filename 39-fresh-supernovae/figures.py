"""Every figure in D6, drawn from data/results.json and the cached tables.

    python run.py --figures-only
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import config  # noqa: E402
import cosmology as cosmo  # noqa: E402
from plotstyle import (  # noqa: E402
    AXIS, BAND_COLOR, BAND_LABEL, BAND_MARKER, BLUE, BLUE_150, BLUE_300, C_LIVE, C_MODEL,
    C_PANTHEON, C_REFERENCE, C_ZTF, EXCLUDED, GRID, INK, INK_SECONDARY, MUTED, ORANGE, SURFACE,
    VIOLET, annotate, caption, save, use_style,
)

use_style()
FIG = config.FIGURE_DIR


def _tables():
    t = pd.read_parquet(config.CACHE_DIR / "hubble_table.parquet")
    p = pd.read_parquet(config.CACHE_DIR / "pantheon_table.parquet")
    live_path = config.CACHE_DIR / "live_table.parquet"
    lt = pd.read_parquet(live_path) if live_path.exists() else pd.DataFrame()
    return t, p, lt


def _signed(v, spec=".3f"):
    """A signed number with a typographic minus sign."""
    return f"{v:+{spec}}".replace("-", "\u2212")


def _pm(v, e, spec=".3f"):
    return f"{v:{spec}} ± {e:{spec}}"


def _d1_curve(R, z):
    d1 = R["d1_universe"]
    return cosmo.distance_modulus(z, z, d1["Om"], d1["H0"])


def _volume_band(ax, top=True):
    ax.axvspan(0.0, config.Z_VOLUME_LIMIT, color=EXCLUDED, zorder=0, lw=0)
    if top:
        ax.text(0.0585, 0.97, "volume-limited\n(no selection bias)", transform=ax.get_xaxis_transform(),
                ha="right", va="top", fontsize=7, color=MUTED)


# ---------------------------------------------------------------------------
# 1. The Hubble diagram
# ---------------------------------------------------------------------------

def figure_hubble_diagram(R: dict) -> list:
    t, p, lt = _tables()
    land = R["landing"]["full"]
    off = land["offset_hubble_flow"]
    d1s = t["d1_test"]
    pm = (p["zHD"] > 0.01) & (p["IS_CALIBRATOR"] == 0)

    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.32)
    ax = fig.add_subplot(gs[0])
    zz = np.geomspace(0.003, 2.4, 400)
    ax.scatter(p.loc[pm, "zHD"], p.loc[pm, "mu"], s=7, color=C_PANTHEON, lw=0, alpha=0.8,
               label=f"D1's supernovae (Pantheon+, {int(pm.sum()):,})", zorder=2)
    ax.scatter(t.loc[d1s, "redshift"], t.loc[d1s, "mu_d1_full"], s=6, color=C_ZTF, lw=0, alpha=0.55,
               label=f"ZTF DR2, new to D1 ({int(d1s.sum()):,})", zorder=3)
    if not lt.empty and lt["d1_test"].any():
        lg = lt[lt["d1_test"]]
        ax.scatter(lg["redshift"], lg["mu_d1_full"], s=60, marker="*", color=C_LIVE,
                   edgecolor=SURFACE, lw=0.6, zorder=5,
                   label=f"classified in the last weeks ({len(lg)})")
    ax.plot(zz, _d1_curve(R, zz), color=C_MODEL, lw=1.4, zorder=4,
            label=f"D1's universe (Ωm = {R['d1_universe']['Om']:.3f}, H0 = {R['d1_universe']['H0']:.1f})")
    ax.set_xscale("log")
    ax.set_xlim(0.004, 2.5)
    ax.set_ylim(31, 47)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("distance modulus μ  [mag]")
    ax.set_title(f"D1's Hubble diagram, with {int(d1s.sum()):,} supernovae it never saw")
    ax.legend(loc="lower right", markerscale=1.6)
    ax.text(0.01, 0.97, "every new supernova standardised exactly as D1 would:\n"
            "D1's light-curve model, α = 0.148, β = 3.09, M = −19.244",
            transform=ax.transAxes, va="top", fontsize=8, color=INK_SECONDARY)

    ax2 = fig.add_subplot(gs[1])
    _volume_band(ax2)
    r = t.loc[d1s, "resid_d1_full"]
    ax2.scatter(t.loc[d1s, "redshift"], r, s=4, color=C_ZTF, lw=0, alpha=0.25, zorder=2)
    for key, color, label, dx in (("pantheon", INK_SECONDARY, "D1's supernovae, binned", 0.985),
                                   ("ztf_d1", C_ZTF, "ZTF DR2, binned", 1.015)):
        b = pd.DataFrame(R["binned"][key])
        ax2.errorbar(b["x"] * dx, b["mean"], yerr=b["err"], fmt="o-", ms=4, lw=1.4, color=color,
                     mec=SURFACE, mew=0.8, capsize=0, zorder=4, label=label)
    if not lt.empty and lt["d1_test"].any():
        lg = lt[lt["d1_test"]]
        ax2.scatter(lg["redshift"], lg["resid_d1_full"], s=60, marker="*", color=C_LIVE,
                    edgecolor=SURFACE, lw=0.6, zorder=5, label="classified in the last weeks")
    ax2.axhline(0, color=C_MODEL, lw=1.2, zorder=3)
    ax2.axhline(off["mean"], color=C_ZTF, lw=0.9, ls=(0, (4, 3)), zorder=3)
    annotate(ax2, f"ZTF's mean offset {_signed(off['mean'])} ± {off['err']:.3f} mag\n"
             "(DR2 fluxes are blinded: not cosmology)", (0.2, off["mean"]), (0.12, 0.62),
             color=INK_SECONDARY)
    ax2.set_xscale("log")
    ax2.set_xlim(0.008, 0.31)
    ax2.set_ylim(-0.95, 0.95)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("μ − D1's prediction  [mag]")
    ax2.set_title("Residuals from D1's universe")
    ax2.legend(loc="lower left", ncol=3, fontsize=7.5)
    caption(fig, "ZTF SN Ia DR2 (Rigault et al. 2025) light curves fitted here with D1's own SALT2 model; "
            "D1 = Pantheon+ (Brout et al. 2022) fitted in 05-supernova-cosmology. Binned points: weighted means. "
            "Grey band: z < 0.06, where ZTF classified essentially every normal SN Ia.")
    return save(fig, FIG / "01_hubble_diagram")


# ---------------------------------------------------------------------------
# 2. Standardisable, not standard
# ---------------------------------------------------------------------------

def figure_standardisable(R: dict) -> list:
    t, _, _ = _tables()
    zf = R["ztf_standardisation"]["with_mass_step"]
    a, b = zf["params"]["alpha"], zf["params"]["beta"]
    d1a, d1b = R["d1_recipe"]["alpha"], R["d1_recipe"]["beta"]
    vl = t["cosmo"] & t["redshift"].between(0.015, config.Z_VOLUME_LIMIT)
    raw = t.loc[vl, "resid_raw"]
    raw = raw - raw.median()
    x1, c = t.loc[vl, "x1"], t.loc[vl, "c"]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), gridspec_kw={"width_ratios": [1, 1, 1.05],
                                                                      "wspace": 0.38})
    # (a) brighter-slower, after removing colour
    ax = axes[0]
    y = raw - b * c
    y = y - y.median()
    ax.scatter(x1, y, s=5, color=C_ZTF, alpha=0.35, lw=0)
    edges = np.linspace(-3, 3, 13)
    bx = pd.Series(y.values).groupby(pd.cut(x1.values, edges), observed=True)
    mids = [iv.mid for iv in bx.median().index]
    ax.plot(mids, bx.median().values, "o", color=INK, ms=4, zorder=4, label="binned median")
    xx = np.linspace(-3, 3, 10)
    ax.plot(xx, -a * xx, color=C_ZTF, lw=1.6, label=f"fitted on ZTF: α = {a:.3f}")
    ax.plot(xx, -d1a * xx, color=C_REFERENCE, lw=1.2, ls=(0, (4, 3)), label=f"D1: α = {d1a:.3f}")
    ax.set_ylim(1.2, -1.2)
    ax.set_xlabel("stretch x1  (← faster fading · slower →)")
    ax.set_ylabel("peak magnitude − distance  [mag]\n(colour removed; up = brighter)")
    ax.set_title("Brighter supernovae fade more slowly")
    ax.legend(loc="lower left", fontsize=7.5)

    # (b) brighter-bluer, after removing stretch
    ax = axes[1]
    y = raw + a * x1
    y = y - y.median()
    ax.scatter(c, y, s=5, color=C_ZTF, alpha=0.35, lw=0)
    edges = np.linspace(-0.2, 0.6, 17)
    bc = pd.Series(y.values).groupby(pd.cut(c.values, edges), observed=True)
    ax.plot([iv.mid for iv in bc.median().index], bc.median().values, "o", color=INK, ms=4, zorder=4)
    cc = np.linspace(-0.25, 0.65, 10)
    c0 = float(np.median(c))
    ax.plot(cc, b * (cc - c0), color=C_ZTF, lw=1.6, label=f"fitted on ZTF: β = {b:.2f}")
    ax.plot(cc, d1b * (cc - c0), color=C_REFERENCE, lw=1.2, ls=(0, (4, 3)), label=f"D1: β = {d1b:.2f}")
    ax.set_ylim(1.6, -1.2)
    ax.set_xlim(-0.25, 0.65)
    ax.set_xlabel("colour c  (← bluer · redder →)")
    ax.set_ylabel("peak magnitude − distance  [mag]\n(stretch removed; up = brighter)")
    ax.set_title("Redder supernovae look fainter")
    ax.legend(loc="lower left", fontsize=7.5)

    # (c) the scatter collapses
    ax = axes[2]
    sc = R["scatter"]["dr2_clean_0.03-0.1"]
    m = t["cosmo"] & t["redshift"].between(*config.Z_DR2_CLEAN)
    before = t.loc[m, "resid_raw"] - t.loc[m, "resid_raw"].median()
    after = t.loc[m, "resid_ztf"] - t.loc[m, "resid_ztf"].median()
    bins = np.linspace(-1.5, 1.5, 61)
    ax.hist(before, bins=bins, color=MUTED, alpha=0.55, label=f"peak magnitude only: nMAD {sc['raw_salt3']['nmad']:.2f} mag")
    ax.hist(after, bins=bins, color=C_ZTF, alpha=0.85,
            label=f"standardised: nMAD {sc['ztf_fit_salt3']['nmad']:.3f} mag")
    ax.set_xlabel("Hubble residual  [mag]")
    ax.set_ylabel("supernovae per 0.05 mag")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
    ax.set_title("Two corrections: 30% → 7% in distance")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_xlim(-1.5, 1.5)
    caption(fig, f"ZTF DR2 normal SNe Ia new to D1, our SALT3 fits. (a, b): 0.015 < z < 0.06 "
            f"({int(vl.sum())} supernovae, no selection bias); the binned stretch points bend near x1 = −0.5, "
            f"the 'broken α' found by Ginolin et al. 2025a. (c): 0.03 < z < 0.1 "
            f"({sc['ztf_fit_salt3']['n']}); nMAD = robust standard deviation. A residual of 0.15 mag is 7% in distance.")
    return save(fig, FIG / "02_standardisable")


# ---------------------------------------------------------------------------
# 3. Light curves and their SALT3 fits
# ---------------------------------------------------------------------------

def _model_curve(row, band, times, kind="salt3"):
    import saltfit

    model = saltfit.get_model(kind)
    model.set(z=row["redshift"], t0=row["t0"], x0=row["x0"], x1=row["x1"], c=row["c"],
              mwebv=row["mwebv"])
    return model.bandflux(band, times, zp=30, zpsys="ab")


def _pick_examples(t: pd.DataFrame) -> list[tuple[str, str]]:
    ok = t["basic_cuts"] & (t["lc_n_phases"] >= 30) & t["z_group"].eq("galaxy")
    picks = []

    def take(mask, label, key, asc):
        cand = t[mask & ~t.index.isin([p for p, _ in picks])].sort_values(key, ascending=asc)
        if len(cand):
            picks.append((cand.index[0], label))

    norm = ok & t["normal"]
    take(norm & t["x1"].between(-0.3, 0.3) & t["c"].between(-0.08, 0.0), "a typical one", "x1_err", True)
    take(norm & (t["x1"] > 1.2), "slow (high stretch)", "x1_err", True)
    take(norm & (t["x1"] < -2.0), "fast (low stretch)", "x1_err", True)
    take(norm & (t["c"] > 0.35), "red (dusty)", "c_err", True)
    take(t["basic_cuts"] & t["normal"] & (t["redshift"] > 0.13) & (t["lc_n_phases"] >= 15),
         "far away (z > 0.13)", "x1_err", True)
    pec = t["ok"] & t["sub_type"].isin(["91bg"]) & t["good_sampling"].fillna(False).astype(bool)
    take(pec, "peculiar: 91bg-like", "lc_n_phases", False)
    return picks


def figure_light_curves(R: dict) -> list:
    from lightcurves import read_dr2_lightcurve

    t, _, _ = _tables()
    picks = _pick_examples(t)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    for ax, (name, label) in zip(axes.flat, picks):
        row = t.loc[name]
        lc, _ = read_dr2_lightcurve(name)
        phase = (lc["mjd"] - row["t0"]) / (1 + row["redshift"])
        lc = lc[lc["good"] & (phase > -25) & (phase < 70)]
        ph = (lc["mjd"] - row["t0"]) / (1 + row["redshift"])
        peak = 0.0
        # The model is defined for -20..+50 rest-frame days; sncosmo returns zero outside.
        tt = row["t0"] + np.linspace(-19.9, 49.9, 240) * (1 + row["redshift"])
        for band in ("ztfg", "ztfr", "ztfi"):
            m = lc["band"] == band
            if not m.any():
                continue
            mod = _model_curve(row, band, tt)
            peak = max(peak, float(np.max(mod)))
        for band in ("ztfg", "ztfr", "ztfi"):
            m = lc["band"] == band
            if not m.any():
                continue
            ax.errorbar(ph[m], lc.loc[m, "flux"] / peak, yerr=lc.loc[m, "fluxerr"] / peak,
                        fmt=BAND_MARKER[band], ms=3.2, color=BAND_COLOR[band], ecolor=BAND_COLOR[band],
                        elinewidth=0.6, alpha=0.75, mec=SURFACE, mew=0.4, zorder=3)
            mod = _model_curve(row, band, tt) / peak
            xs = (tt - row["t0"]) / (1 + row["redshift"])
            ax.plot(xs, mod, color=BAND_COLOR[band], lw=1.5, zorder=4)
            i = int(np.argmax(mod))
            ax.text(xs[i] + 3, mod[i] + 0.03, BAND_LABEL[band], color=INK, fontsize=8, va="bottom",
                    bbox={"boxstyle": "round,pad=0.15", "fc": SURFACE, "ec": "none", "alpha": 0.85})
        ax.axvspan(-10, 40, color=EXCLUDED, zorder=0, lw=0)
        ax.axhline(0, color=AXIS, lw=0.8)
        ax.set_title(label, fontsize=10)
        ax.text(0.98, 0.95, f"{row['iau_name'] if isinstance(row['iau_name'], str) else name}\n"
                f"z = {row['redshift']:.3f}\nx1 = {row['x1']:+.2f}, c = {row['c']:+.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK_SECONDARY)
        ax.set_xlim(-25, 70)
        ax.set_ylim(-0.15, 1.35)
    for ax in axes[-1]:
        ax.set_xlabel("days from peak, in the supernova's own time")
    for ax in axes[:, 0]:
        ax.set_ylabel("flux / peak flux")
    handles = [Line2D([], [], ls="-", marker=BAND_MARKER[b], color=BAND_COLOR[b], ms=5,
                      label=f"ZTF {BAND_LABEL[b]}") for b in ("ztfg", "ztfr", "ztfi")]
    fig.legend(handles=handles, loc="upper right", ncol=3, bbox_to_anchor=(0.98, 1.0))
    fig.suptitle("ZTF light curves and their SALT3 fits", x=0.06, ha="left", y=1.0)
    caption(fig, "Points: ZTF DR2 forced photometry (good-quality epochs). Lines: our SALT3 fit. Shaded: the "
            "−10 to +40 day window the fit uses. Time is divided by (1 + z), so all curves run on the "
            "supernova's own clock. Fluxes on DR2's blinded zero point, normalised to the fitted peak.")
    return save(fig, FIG / "03_light_curves")


# ---------------------------------------------------------------------------
# 4. Where the scatter comes from
# ---------------------------------------------------------------------------

def figure_scatter_budget(R: dict) -> list:
    sb = R["scatter_budget"]
    g = sb["galaxy_redshifts"]
    zf = sb["sn_feature_redshifts"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    zz = np.geomspace(0.012, 0.07, 200)
    sfit = sb["median_sigma_fit"]
    base = np.sqrt(g["sigma_sn"] ** 2 + sfit ** 2)
    vel = cosmo.sigma_mu_from_velocity(zz, g["sigma_v"])
    zerr = cosmo.sigma_mu_from_velocity(zz, zf["sigma_z"] * cosmo.C_KM_S)
    ax.plot(zz, np.sqrt(base ** 2 + vel ** 2), color=C_ZTF, lw=1.6,
            label=f"model: σ_SN = {g['sigma_sn']:.3f} mag ⊕ {g['sigma_v']:.0f} km/s")
    ax.plot(zz, np.sqrt(base ** 2 + vel ** 2 + zerr ** 2), color=ORANGE, lw=1.6,
            label=f"… ⊕ redshift error σz = {zf['sigma_z']:.4f}")
    ax.axhline(base, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax.text(0.0122, base - 0.012, "supernovae alone (+ fit errors)", fontsize=7.5, color=MUTED, va="top")
    for key, color, label in (("binned_galaxy", C_ZTF, "host-galaxy redshifts"),
                              ("binned_sn_features", ORANGE, "redshifts from SN spectra")):
        b = pd.DataFrame(sb[key])
        err = b["nmad"] / np.sqrt(2 * (b["n"] - 1))
        ax.errorbar(b["x"], b["nmad"], yerr=err, fmt="o", ms=5, color=color, mec=SURFACE, mew=0.8,
                    label=f"{label} (nMAD per bin)")
    ax.set_xscale("log")
    ax.set_xlim(0.012, 0.07)
    ax.set_ylim(0, 0.6)
    from matplotlib.ticker import FixedLocator, NullFormatter

    ax.xaxis.set_major_locator(FixedLocator([0.015, 0.02, 0.03, 0.04, 0.05, 0.06]))
    ax.set_xticklabels(["0.015", "0.02", "0.03", "0.04", "0.05", "0.06"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("redshift z")
    ax.set_ylabel("Hubble-residual scatter  [mag]")
    ax.set_title("Redshift errors show up; galaxy motions barely do")
    ax.legend(loc="upper right", fontsize=7.5)

    ax = axes[1]
    grid = np.array(g["profile_grid"])
    prof = np.array(g["profile_dchi2"])
    ax.plot(grid, prof, color=C_ZTF, lw=1.6, label=f"with an outlier component ({100 * g['outlier_fraction']:.1f}%)")
    gg = g["gaussian"]
    ax.plot(gg["profile_grid"], gg["profile_dchi2"], color=MUTED, lw=1.2, ls=(0, (4, 3)),
            label=f"pure Gaussian: {gg['sigma_v']:.0f} km/s (outliers win)")
    ax.axhline(1, color=MUTED, lw=0.8, ls=(0, (4, 3)))
    lo, hi = g["sigma_v_68"]
    ax.axvspan(lo, hi, color=BLUE_150, lw=0, zorder=0)
    ax.axvline(250, color=C_REFERENCE, lw=1.1, ls=(0, (1, 2)), label="250 km/s, Pantheon+'s assumption")
    ax.set_xlim(0, 900)
    ax.set_ylim(0, 7)
    ax.set_xlabel("random galaxy velocity σ_v  [km/s]")
    ax.set_ylabel("Δχ²  (profile likelihood)")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title(f"Random galaxy velocities: below {hi:.0f} km/s (68%)")
    caption(fig, f"Normal ZTF DR2 SNe Ia at z < 0.06 standardised on ZTF (SALT3). {g['n']} with host-galaxy "
            f"redshifts fix σ_SN and σ_v (few lie below z = 0.02, where velocities would dominate); {zf['n']} with "
            f"redshifts from the supernova's own spectrum then need an extra σz = {zf['sigma_z']:.4f} "
            f"({zf['sigma_v_equiv']:.0f} km/s; ZTF team: ~0.003). Fits allow a small fraction of outliers "
            f"({100 * g['outlier_fraction']:.1f}%); without that, a handful of them fake a large σ_v.")
    return save(fig, FIG / "04_scatter_budget")


# ---------------------------------------------------------------------------
# 5. Malmquist bias
# ---------------------------------------------------------------------------

def figure_selection(R: dict) -> list:
    s = R["selection"]
    mod = pd.DataFrame(s["model"]["grid"])
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    specs = (("data_resid_ztf", "bias", "mean Hubble residual  [mag]", "The far ones are the bright ones"),
             ("data_mean_c", "mean_c", "mean colour c", "… and the blue ones"),
             ("data_mean_x1", "mean_x1", "mean stretch x1", "… and the slow ones"))
    for ax, (dk, mk, ylab, title) in zip(axes, specs):
        _volume_band(ax, top=False)
        if dk == "data_resid_ztf":
            ax.text(0.003, 0.03, "volume-limited:\nno selection bias", transform=ax.get_xaxis_transform(),
                    fontsize=7, color=MUTED, va="bottom")
        b = pd.DataFrame(s[dk])
        ref = 0.0
        if dk == "data_resid_ztf":
            vl = b[b["hi"] <= config.Z_VOLUME_LIMIT + 1e-9]
            ref = float(np.average(vl["mean"], weights=1 / vl["err"] ** 2))
        ax.errorbar(b["x"], b["mean"] - ref, yerr=b["err"], fmt="o", ms=5, color=C_ZTF, mec=SURFACE,
                    mew=0.8, label="ZTF DR2 (binned)", zorder=4)
        ax.plot(mod["z"], mod[mk], color=INK, lw=1.4, label="selection model", zorder=3)
        ax.set_xlim(0, 0.2)
        ax.set_xlabel("redshift z")
        ax.set_ylabel(ylab)
        ax.set_title(title)
    axes[0].legend(loc="lower left", fontsize=7.5)
    axes[0].axhline(0, color=AXIS, lw=0.8)
    caption(fig, "Selection model (a simulation from real inputs): supernovae with the stretch and colour "
            "mix of ZTF's volume-limited sample and its intrinsic scatter, kept with the probability that ZTF "
            "classified one of that peak brightness, 1/(1+exp(4.5 (m − 18.8))) (Rigault et al. 2025). "
            "Pantheon+ removes this bias with simulations; ZTF DR2 has none, hence z < 0.06 for tests.")
    return save(fig, FIG / "05_selection_bias")


# ---------------------------------------------------------------------------
# 6. Peculiar types
# ---------------------------------------------------------------------------

SUBTYPE_LABELS = {
    "norm": "normal", "91t": "91T-like", "99aa": "99aa-like", "norm/99aa": "normal/99aa",
    "99aa/91t": "99aa/91T", "unspecified": "Ia (no sub-type)", "norm/04gs": "normal/04gs",
    "91bg": "91bg-like", "91bg/86G": "86G-like", "iax": "Iax (02cx-like)", "03fg": "03fg-like",
    "csm": "Ia-CSM", "02es": "02es-like", "00cx": "00cx-like", "18byg": "18byg-like",
}


def figure_peculiar(R: dict) -> list:
    t, _, _ = _tables()
    usable = t["ok"] & t["good_sampling"].fillna(False).astype(bool) & (t["c_err"] < 0.15) \
        & (t["x1_err"] < 1.5) & t["fresh"]
    t = t[usable].copy()
    t["subtype"] = t["sub_type"].fillna("unspecified")
    zero = float(np.median(t.loc[t["cosmo"], "resid_ztf"]))
    t["r"] = t["resid_ztf"] - zero
    order = [k for k in SUBTYPE_LABELS if (t["subtype"] == k).sum() >= 3]
    order = sorted(order, key=lambda k: (t["peculiar"][t["subtype"] == k].mean() > 0.5,
                                         -(t["subtype"] == k).sum()))
    fig, ax = plt.subplots(figsize=(11, 5.6))
    rng = np.random.default_rng(1)
    for i, k in enumerate(order):
        g = t[t["subtype"] == k]
        pec = g["peculiar"].mean() > 0.5
        color = ORANGE if pec else C_ZTF
        y = np.clip(g["r"], -2.4, 3.4)
        ax.scatter(i + rng.uniform(-0.28, 0.28, len(g)), y, s=6 if len(g) > 200 else 14, color=color,
                   alpha=0.35 if len(g) > 200 else 0.75, lw=0, zorder=2)
        med = float(np.median(g["r"]))
        ax.plot([i - 0.36, i + 0.36], [med, med], color=INK, lw=2, zorder=3)
        ax.text(i, -2.75, f"n={len(g)}", ha="center", fontsize=7, color=MUTED)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.axvline(sum(t["peculiar"][t["subtype"] == k].mean() <= 0.5 for k in order) - 0.5,
               color=AXIS, lw=0.8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([SUBTYPE_LABELS[k] for k in order], rotation=30, ha="right")
    ax.set_ylim(3.6, -3.0)
    ax.set_ylabel("standardised Hubble residual  [mag]\n(up = brighter than a normal SN Ia)")
    ax.set_title("Why cosmology drops peculiar Type Ia supernovae")
    handles = [Line2D([], [], ls="", marker="o", color=C_ZTF, label="kept for cosmology"),
               Line2D([], [], ls="", marker="o", color=ORANGE, label="peculiar: dropped"),
               Line2D([], [], color=INK, lw=2, label="median")]
    ax.legend(handles=handles, loc="lower right", fontsize=8)
    caption(fig, "ZTF DR2 spectroscopic sub-types (Dimitriadis et al. 2025; Burgaz et al. 2025), standardised "
              "with the stretch and colour relations fitted on normal ZTF supernovae. Points beyond the "
              "axes are drawn at the edge. Quality cuts on the fit only (no x1 or c range cut).")
    return save(fig, FIG / "06_peculiar_types")


# ---------------------------------------------------------------------------
# 7. Validation
# ---------------------------------------------------------------------------

def figure_validation(R: dict) -> list:
    t, _, _ = _tables()
    both = t["ok"] & np.isfinite(t["s3zt_x1"]) & t["basic_cuts"]
    v = R["validation"]
    conv = pd.read_csv(config.CACHE_DIR / "snana_convention_check.csv")
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.6), gridspec_kw={"wspace": 0.38})
    for ax, k, lim, label in ((axes[0], "x1", (-3.3, 3.3), "stretch x1"),
                              (axes[1], "c", (-0.25, 0.85), "colour c")):
        ax.scatter(t.loc[both, "s3zt_" + k], t.loc[both, k], s=4, color=C_ZTF, alpha=0.35, lw=0)
        ax.plot(lim, lim, color=INK, lw=0.9)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel(f"ZTF team's SALT3 {label}")
        ax.set_ylabel(f"our SALT3 {label}")
        ax.set_title(f"{label}\nmedian Δ {_signed(v[k]['median_diff'])}, nMAD {v[k]['nmad_diff']:.3f}",
                     fontsize=9.5)
    ax = axes[2]
    d = t.loc[both, "mB"] - t.loc[both, "s3zt_mB"]
    ax.hist(d.clip(-0.1, 0.1), bins=np.linspace(-0.1, 0.1, 41), color=C_ZTF)
    ax.set_xlabel("our mB − ZTF team's mB  [mag]")
    ax.set_ylabel("supernovae")
    ax.set_title(f"peak magnitude\nnMAD {v['mB']['nmad_diff']:.3f} mag", fontsize=9.5)
    ax = axes[3]
    cc = R["convention_check"]
    ax.hist(conv["d_fixed"], bins=np.linspace(0.2, 0.32, 25), color=VIOLET, alpha=0.85,
            label="D1's model, x1, c, t0 fixed")
    ax.axvline(cc["mag_offset_in_SALT2_INFO"], color=INK, lw=1.4)
    ax.text(cc["mag_offset_in_SALT2_INFO"] + 0.002, ax.get_ylim()[1] * 0.95, "MAG_OFFSET\n= 0.27",
            fontsize=7.5, va="top", color=INK_SECONDARY)
    ax.set_xlabel("sncosmo − SNANA peak magnitude  [mag]")
    ax.set_ylabel("Pantheon+ Foundation supernovae")
    ax.set_title(f"the hidden convention\n{_pm(cc['fixed_x1_c_t0']['mean'], cc['fixed_x1_c_t0']['err'])} mag",
                 fontsize=9.5)
    caption(fig, f"Left three: our SALT3 fits against the ZTF team's own SALT3 fits for the same {v['n']:,} "
            "supernovae (same light curves, same model). Right: Pantheon+'s own Foundation light curves "
            "refitted here with D1's model; SNANA adds the model's MAG_OFFSET to every magnitude and "
            "sncosmo does not, so D1's M only applies after subtracting 0.27 mag.")
    return save(fig, FIG / "07_validation")


# ---------------------------------------------------------------------------
# 8. The newest supernovae
# ---------------------------------------------------------------------------

def figure_live(R: dict) -> list:
    lv = R.get("live", {})
    if not lv.get("available"):
        return []
    t, _, lt = _tables()
    lcs = pd.read_parquet(config.CACHE_DIR / "live_lightcurves.parquet")
    good = lt[lt["d1_test"]].sort_values("days_since_peak", ascending=False)
    show = good.head(4).index.tolist()
    fig = plt.figure(figsize=(13, 7.4))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.15], hspace=0.45, wspace=0.28)
    import saltfit

    for i, name in enumerate(show):
        ax = fig.add_subplot(gs[0, i])
        row = lt.loc[name]
        lc = lcs[lcs["tns_name"] == name]
        lc = lc[lc["good"]]
        ph = (lc["mjd"] - row["t0"]) / (1 + row["redshift"])
        model = saltfit.get_model("salt3")
        model.set(z=row["redshift"], t0=row["t0"], x0=row["x0"], x1=row["x1"], c=row["c"], mwebv=row["mwebv"])
        tt = row["t0"] + np.linspace(-20, max(30, ph.max() + 5), 200) * (1 + row["redshift"])
        peak = max(float(np.max(model.bandflux(b, tt, zp=30, zpsys="ab")))
                   for b in lc["band"].unique())
        for band in ("ztfg", "ztfr", "ztfi"):
            m = (lc["band"] == band) & (ph > -30)
            if not m.any():
                continue
            ax.errorbar(ph[m], lc.loc[m, "flux"] / peak, yerr=lc.loc[m, "fluxerr"] / peak,
                        fmt=BAND_MARKER[band], ms=3.2, color=BAND_COLOR[band], elinewidth=0.6,
                        mec=SURFACE, mew=0.4, alpha=0.85)
            ax.plot((tt - row["t0"]) / (1 + row["redshift"]),
                    model.bandflux(band, tt, zp=30, zpsys="ab") / peak, color=BAND_COLOR[band], lw=1.4)
        ax.axhline(0, color=AXIS, lw=0.8)
        ax.set_title(f"{name}  (z = {row['redshift']:.3f})", fontsize=9)
        ax.set_xlim(-30, max(35, ph.max() + 5))
        ax.set_ylim(-0.2, 1.4)
        ax.set_xlabel("days from peak (rest frame)")
        if i == 0:
            ax.set_ylabel("flux / peak")
        disc = str(row["discovery_date"])[:10]
        ax.text(0.97, 0.95, f"discovered {disc}\nx1 = {row['x1']:+.2f}, c = {row['c']:+.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color=INK_SECONDARY)

    ax = fig.add_subplot(gs[1, :])
    _volume_band(ax)
    d1s = t["d1_test"]
    ax.scatter(t.loc[d1s, "redshift"], t.loc[d1s, "resid_d1_full"], s=4, color=C_ZTF, lw=0, alpha=0.18,
               label="ZTF DR2 (2018–2020)")
    b = pd.DataFrame(R["binned"]["ztf_d1"])
    ax.errorbar(b["x"], b["mean"], yerr=b["err"], fmt="o-", ms=4, color=C_ZTF, mec=SURFACE, mew=0.8,
                lw=1.2, label="ZTF DR2, binned")
    ax.errorbar(good["redshift"], good["resid_d1_full"], yerr=good["sigma_meas_d1"], fmt="none",
                ecolor=C_LIVE, elinewidth=0.8, alpha=0.7)
    o = lv.get("offset")
    lab = f"classified since {lv['meta']['discovered_from']} ({len(good)})"
    if o:
        lab += f": mean {_signed(o['mean'])} ± {o['err']:.3f}"
    ax.scatter(good["redshift"], good["resid_d1_full"], s=90, marker="*", color=C_LIVE, edgecolor=SURFACE,
               lw=0.6, zorder=5, label=lab)
    for name, r in good.iterrows():
        if abs(r["resid_d1_full"]) > 0.45 or name == show[0]:
            ax.annotate(name.replace("SN ", ""), (r["redshift"], r["resid_d1_full"]), xytext=(6, -3),
                        textcoords="offset points", fontsize=7.5, color=INK_SECONDARY)
    ax.axhline(0, color=C_MODEL, lw=1.2)
    if o:
        ax.axhline(o["mean"], color=C_LIVE, lw=0.9, ls=(0, (4, 3)))
    ax.set_xlim(0, 0.16)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("μ − D1's prediction  [mag]")
    ax.set_title("Supernovae a few weeks old, on D1's Hubble diagram")
    ax.legend(loc="upper right", fontsize=8, ncol=1)
    caption(fig, f"Snapshot {lv['meta']['snapshot']}: SNe Ia classified on the Transient Name Server and "
            "discovered since " + lv['meta']['discovered_from'] + "; ZTF forced photometry from alert packets "
            "(ALeRCE broker), on ZTF's real (unblinded) calibration. Same fits, cuts and D1 recipe as the DR2 "
            "points. Redshifts mostly from the classification spectra (±0.005); a wrong one moves a point a lot "
            "(2026znw's normal light curve sits 0.95 mag off).")
    return save(fig, FIG / "08_live")


def figure_flows(R: dict) -> list:
    from matplotlib.colors import LinearSegmentedColormap
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    fl = R.get("flows")
    if not fl:
        return []
    t, _, _ = _tables()
    m = t["d1_test"] & t["z_group"].eq("galaxy") & t["redshift"].between(0.01, config.Z_VOLUME_LIMIT)
    d = t[m]
    g = SkyCoord(d["ra"].to_numpy() * u.deg, d["dec"].to_numpy() * u.deg).galactic
    k = cosmo.sigma_mu_from_velocity(d["z_cmb"].to_numpy(), 1.0)
    v = -(d["resid_d1_full"].to_numpy() - fl["dipole"]["offset"]) / k
    cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])

    def xy(l, b):
        lw = (np.asarray(l) + 180) % 360 - 180
        return -np.radians(lw), np.radians(b)

    fig = plt.figure(figsize=(13.5, 5.0))
    ax = fig.add_axes([0.02, 0.1, 0.58, 0.8], projection="mollweide")
    ax.set_facecolor(SURFACE)
    # The sky ZTF cannot see (declination below -31 deg).
    ra, dec = np.meshgrid(np.linspace(0, 360, 361), np.linspace(-90, -31, 60))
    hid = SkyCoord(ra.ravel() * u.deg, dec.ravel() * u.deg).galactic
    hx, hy = xy(hid.l.deg, hid.b.deg)
    ax.scatter(hx, hy, s=3, color=EXCLUDED, lw=0, zorder=0.5)
    ax.xaxis.set_zorder(6)
    order = np.argsort(np.abs(v))
    x, y = xy(g.l.deg[order], g.b.deg[order])
    sc = ax.scatter(x, y, c=np.clip(v[order], -1500, 1500), cmap=cmap, vmin=-1500, vmax=1500, s=14,
                    lw=0.3, edgecolor=GRID, zorder=3)
    marks = [((fl["dipole"]["l"], fl["dipole"]["b"]), "*", INK, 14, "this fit"),
             ((304, 6), "P", C_REFERENCE, 9, "Carrick+15 bulk flow"),
             ((264.021, 48.253), "X", INK_SECONDARY, 8, "CMB dipole (our motion)")]
    for (l, b), mk, col, ms, lab in marks:
        mx, my = xy(l, b)
        ax.plot(mx, my, mk, color=col, ms=ms, mec=SURFACE, mew=1.0, zorder=5, label=lab)
    ax.set_xticks(np.radians([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]))
    ax.set_xticklabels(["150°", "120°", "90°", "60°", "30°", "0°", "330°", "300°", "270°", "240°", "210°"],
                       fontsize=7, color=MUTED)
    ax.tick_params(axis="y", labelsize=7, colors=MUTED)
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_title("ZTF supernovae at z < 0.06 on the sky (Galactic coordinates)", pad=30)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3, fontsize=7.5)
    cax = fig.add_axes([0.08, 0.06, 0.46, 0.025])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("implied radial velocity of the host galaxy  [km/s]   (blue: approaching · red: receding)",
                 fontsize=7.5, color=INK_SECONDARY)
    cb.ax.tick_params(labelsize=7, colors=MUTED)
    cb.outline.set_visible(False)

    ax2 = fig.add_axes([0.68, 0.16, 0.30, 0.68])
    apex = np.array(fl["dipole"]["vector"]) / fl["dipole"]["dipole_kms"]
    from analysis import _unit_vectors

    cos_t = _unit_vectors(d["ra"], d["dec"]) @ apex
    edges = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06])
    import analysis as A

    for sel, col, lab, dx in ((cos_t > 0.3, ORANGE, "towards the flow apex", 0.997),
                              (cos_t < -0.3, C_ZTF, "away from it", 1.003)):
        b = pd.DataFrame(A.binned(d.loc[sel, "z_cmb"], d.loc[sel, "resid_d1_full"], edges,
                                  np.sqrt(d.loc[sel, "sigma_fit_d1"] ** 2 + 0.15 ** 2)))
        ax2.errorbar(b["x"] * dx, b["mean"] - fl["dipole"]["offset"], yerr=b["err"], fmt="o-", ms=5,
                     color=col, mec=SURFACE, mew=0.8, label=f"{lab}")
    zz = np.linspace(0.012, 0.06, 100)
    amp = fl["dipole"]["dipole_kms"]
    for sgn, col in ((1, ORANGE), (-1, C_ZTF)):
        ax2.plot(zz, -sgn * 0.65 * amp * cosmo.sigma_mu_from_velocity(zz, 1.0), color=col, lw=1, ls=(0, (4, 3)))
    ax2.axhline(0, color=AXIS, lw=0.8)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("residual from D1  [mag]  (offset removed)")
    ax2.set_ylim(0.45, -0.45)
    ax2.set_title("The nearest are faint in every direction", fontsize=10)
    ax2.legend(loc="lower right", fontsize=7.5)
    dp = fl["dipole"]
    mono = fl["monopole"]
    caption(fig, f"{fl['n']} normal ZTF SNe Ia with host-galaxy redshifts, 0.01 < z < 0.06 ({fl['n_clipped']} "
            f"outliers removed). Fitted dipole: {dp['dipole_kms']:.0f} ± {dp['dipole_err']:.0f} km/s towards (l, b) = "
            f"({dp['l']:.0f}° ± {dp['l_err']:.0f}°, {dp['b']:.0f}° ± {dp['b_err']:.0f}°), Δχ² = "
            f"{dp['delta_chi2_vs_offset']:.1f} for 3 parameters; an isotropic 1/z term (all galaxies moving in, or a "
            f"brightness-dependent photometry bias) gives Δχ² = {mono['delta_chi2_vs_offset']:.1f} for 1. Carrick et al. "
            "(2015): 159 ± 23 km/s towards (304°, 6°) from beyond the 2M++ volume. Right: a flow would make the "
            "nearest supernovae towards the apex brighter (dashed), but below z = 0.03 they are as faint as those "
            "on the far side; the apex itself is below ZTF's horizon (grey).")
    return save(fig, FIG / "09_bulk_flow")


def make_all(R: dict) -> list:
    written = []
    for fn in (figure_hubble_diagram, figure_standardisable, figure_light_curves,
               figure_scatter_budget, figure_selection, figure_peculiar, figure_validation,
               figure_live, figure_flows):
        written += fn(R)
    return written
