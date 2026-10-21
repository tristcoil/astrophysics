"""Figure generation.

Every figure is drawn from the analysis results, never by re-running the
analysis, so what is plotted is exactly what was measured.

One colour rule runs through all of them: **ink** is what was measured (ring
masses, velocity dispersions, our own Einstein radii); **blue** is the stars
(for a Chabrier stellar mass function where one is needed) and **aqua** the
same stars with the heavier Salpeter mass function; **orange** is dark matter
and our own fits and models; **violet** is a published result we compare
against; **magenta** is the lensed background galaxy. The palette is the
atlas's colourblind-safe set, and its lighter hues always carry text labels.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import plotstyle as ps
from config import FIGURE_DIR, FIBRE_RADIUS_ARCSEC

ps.use_style()

C_SOURCE = ps.MAGENTA  # the lensed background galaxy
DASHED = (0, (5, 2.5))
DOTTED = (0, (1, 1.8))

PROV_SLACS = ("Data: SLACS (Bolton et al. 2008, ApJ 682, 964; Auger et al. 2009, ApJ 705, 1099) — "
              "HST/ACS imaging and SDSS spectroscopy.")


# --------------------------------------------------------------------------
# Building blocks
# --------------------------------------------------------------------------

def _title(ax, title: str, subtitle: str | None = None) -> None:
    lines = subtitle.count("\n") + 1 if subtitle else 0
    ax.set_title(title, pad=8 + 10.8 * lines)
    if subtitle:
        ax.annotate(subtitle, xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom", fontsize=8.2,
                    color=ps.INK_SECONDARY, linespacing=1.3)


def _label(ax, text, xy, xytext, *, ha="left", va="center", leader=True, fontsize=8,
           textcoords="data", color=None):
    ax.annotate(text, xy=xy, xytext=xytext, textcoords=textcoords, ha=ha, va=va,
                fontsize=fontsize, color=color or ps.INK_SECONDARY,
                arrowprops=({"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6,
                             "shrinkA": 1, "shrinkB": 3} if leader else None))


def _points(ax, x, y, yerr=None, xerr=None, *, color=ps.C_OBSERVED, ms=4.6, label=None,
            zorder=5, alpha=1.0, mfc=None):
    ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt="o", color=color, ecolor=ps.C_OBSERVED_ERR,
                elinewidth=0.8, capsize=0, markersize=ms, markeredgecolor=ps.SURFACE,
                markeredgewidth=0.8, zorder=zorder, label=label, alpha=alpha,
                markerfacecolor=mfc or color)


def _sci(value: float, digits: int = 1) -> str:
    exp = int(np.floor(np.log10(abs(value))))
    return f"{value / 10**exp:.{digits}f} × 10{str(exp).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))}"


def _pm(v, e, fmt=".3f"):
    return f"{v:{fmt}} ± {e:{fmt}}"


# --------------------------------------------------------------------------
# 3. Two redshifts in one spectrum
# --------------------------------------------------------------------------

def figure_two_redshifts(name, z_lens, z_source, search):
    import spectra

    sp = spectra.read_spectrum(name)
    r, w = spectra.clean_residual(sp)
    good = sp.ivar > 0
    fig = plt.figure(figsize=(11.6, 7.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.25, 0.85, 0.75], hspace=0.62, wspace=0.28)
    top = fig.add_subplot(gs[0, :])
    sm = lambda a: np.convolve(a, np.ones(5) / 5, mode="same")  # display smoothing only
    top.plot(sp.wave[good], sm(sp.flux)[good], color=ps.C_OBSERVED, lw=0.6,
             label="SDSS spectrum (5-pixel smoothed)")
    top.plot(sp.wave, sp.model, color=ps.C_STARS, lw=1.2,
             label=f"pipeline's model of the lens galaxy's starlight (z = {z_lens:.3f})")
    top.set_xlim(3800, 9200)
    ymax = np.percentile(sm(sp.flux)[good], 99.7) * 1.25
    top.set_ylim(0, ymax)
    absorption = {"Ca H+K": 0.5 * (spectra.ABSORPTION_LINES["Ca K"] + spectra.ABSORPTION_LINES["Ca H"]),
                  "G band": spectra.ABSORPTION_LINES["G band"],
                  "Mg b": spectra.ABSORPTION_LINES["Mg b"], "Na D": spectra.ABSORPTION_LINES["Na D"]}
    for lab, lam in absorption.items():
        x = lam * (1 + z_lens)
        if 3800 < x < 9200:
            top.annotate(lab, xy=(x, np.interp(x, sp.wave, sp.model)), xytext=(x, ymax * 0.10),
                         ha="center", fontsize=7, color=ps.INK_SECONDARY,
                         arrowprops={"arrowstyle": "-", "color": ps.C_STARS, "lw": 0.5})
    for lab, lams in spectra.EMISSION_LINES.items():
        x = lams[0] * (1 + z_source)
        if 3800 < x < 9200:
            top.axvline(x, ymin=0.80, ymax=0.93, color=C_SOURCE, lw=1.6)
    # Labels: neighbouring lines share one, staggered so none collide.
    em = {k: v[0] * (1 + z_source) for k, v in spectra.EMISSION_LINES.items()}
    groups = [("[O II]", [em["[O II]"]], "center"), ("Hβ", [em["H-beta"]], "right"),
              ("[O III]", [em["[O III] 4959"], em["[O III] 5007"]], "left"),
              ("Hα", [em["H-alpha"]], "center")]
    for lab, xs, ha in groups:
        xs = [x for x in xs if 3800 < x < 9200]
        if xs:
            top.annotate(lab, xy=(np.mean(xs), ymax * 0.95), ha=ha, va="bottom", fontsize=7.5,
                         color=ps.INK_SECONDARY, annotation_clip=False)
    top.set_ylabel("Flux (10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹)")
    top.set_xlabel("Observed wavelength (Å)")
    top.legend(loc="upper left", fontsize=7.6, ncol=2, bbox_to_anchor=(0.0, 0.80))
    _title(top, f"SDSS {name}: two galaxies, two redshifts, one 3″ fibre",
           f"The elliptical's starlight at z = {z_lens:.4f}; marked above it, emission lines of a "
           f"star-forming galaxy far behind at z = {z_source:.4f}, whose light the elliptical "
           "bends into a ring.")

    zoom = [("[O II]", 3728.5, 30), ("H-beta", 4862.7, 30), ("[O III] 4959", 4960.3, 30),
            ("[O III] 5007", 5008.2, 30)]
    for k, (lab, rest, half) in enumerate(zoom):
        ax = fig.add_subplot(gs[1, k])
        c = rest * (1 + z_source)
        sel = (sp.wave > c - half * (1 + z_source) / 1.4) & (sp.wave < c + half * (1 + z_source) / 1.4)
        err = 1 / np.sqrt(np.where(w > 0, w, np.nan))
        ax.fill_between(sp.wave[sel], -err[sel], err[sel], color=ps.GRID, lw=0, step="mid")
        ax.step(sp.wave[sel], r[sel], where="mid", color=ps.C_OBSERVED, lw=0.9)
        ax.axhline(0, color=ps.AXIS, lw=0.8)
        ax.axvline(c, color=C_SOURCE, lw=1.2, ls=DOTTED)
        sn = search.lines_at_best.get(lab, np.nan)
        ax.set_title(f"{lab}  ({sn:.0f}σ)", fontsize=9, pad=4)
        ax.set_xlabel("Å", fontsize=8)
        if k == 0:
            ax.set_ylabel("Data − galaxy model")
        ax.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[2, :])
    ax.plot(search.z_grid, np.sqrt(search.score), color=ps.C_OBSERVED, lw=0.8)
    ax.axvline(z_source, color=ps.C_PUBLISHED, lw=1.0, ls=DASHED)
    _label(ax, f"catalogued source redshift {z_source:.4f}", (z_source, np.sqrt(search.score_best) * 0.9),
           (z_source + 0.08, np.sqrt(search.score_best) * 0.85), color=ps.C_PUBLISHED)
    ax.set_xlim(search.z_grid[0], search.z_grid[-1])
    ax.set_ylim(0, None)
    ax.set_xlabel("Trial redshift of a background galaxy")
    ax.set_ylabel("Line evidence (√Σ S/N²)")
    _title(ax, "The blind search",
           f"every trial redshift scored by how strongly the residuals show [O II], Hβ, [O III] and Hα "
           f"there; the highest peak is at z = {search.z_best:.4f}")
    ps.caption(fig, "Data: SDSS DR17 spectrum (plate-MJD-fibre from Bolton et al. 2008). The galaxy "
                    "model is the SDSS pipeline's own fit at the lens redshift; residuals are "
                    "high-pass filtered (101-pixel running median) before the search; pixels near "
                    "bright night-sky lines and near the lens galaxy's own possible emission lines get "
                    "zero weight. [O III] 4959/5007 searched as a pair (ratio 1:2.98). Vacuum "
                    "wavelengths.", y=0.02)
    return ps.save(fig, FIGURE_DIR / "03_two_redshifts")


# --------------------------------------------------------------------------
# 4. How the lenses were found
# --------------------------------------------------------------------------

def figure_discovery(disc):
    df = disc["table"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.8), gridspec_kw={"wspace": 0.28})
    grade_style = {"A": (ps.C_OBSERVED, "definite lens (grade A)"),
                   "B": (ps.MUTED, "probable lens (B)"),
                   "X": (ps.C_CONTEXT, "not a lens (X)")}
    for grade in ("X", "B", "A"):
        m = df["grade"] == grade
        col, lab = grade_style[grade]
        a1.scatter(df.loc[m, "z_source"], df.loc[m, "z_found"], s=22, color=col,
                   edgecolor=ps.SURFACE, linewidth=0.6, zorder=4 if grade == "A" else 3,
                   label=f"{lab}: {int(df.loc[m, 'recovered'].sum())} of {int(m.sum())} found")
    lim = [0.05, 1.4]
    a1.plot(lim, lim, color=ps.AXIS, lw=0.9, zorder=1)
    a1.set_xlim(lim)
    a1.set_ylim(lim)
    for row in disc["missed"]:
        _label(a1, row["name"], (row["z_source"], row["z_found"]),
               (row["z_source"] + 0.05, row["z_found"] - 0.07), fontsize=7)
    a1.set_xlabel("Source redshift in the SLACS catalogue")
    a1.set_ylabel("Strongest background redshift we find")
    a1.legend(loc="upper left", fontsize=7.5)
    _title(a1, f"Second redshifts found blind: {disc['n_recovered']} of {disc['n']}",
           "every Hubble-observed SLACS candidate; on the line = same redshift to < 0.002")

    rate = disc["confirmation_vs_theta"]
    x = np.arange(len(rate))
    fa = np.array([r_["frac_a"] for r_ in rate])
    fab = np.array([r_["frac_ab"] for r_ in rate])
    n = np.array([r_["n"] for r_ in rate])
    a2.bar(x, fab * 100, width=0.62, color=ps.C_CONTEXT, zorder=2, label="grade A or B")
    a2.bar(x, fa * 100, width=0.62, color=ps.C_OBSERVED, zorder=3, label="grade A")
    for xi, f, nn in zip(x, fab, n):
        a2.annotate(f"{nn} systems", xy=(xi, f * 100), xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    labs = [f"{r_['lo']:.1f}–{r_['hi']:.1f}″" if r_["hi"] < 3 else f"> {r_['lo']:.1f}″"
            for r_ in rate]
    labs[0] = f"< {rate[0]['hi']:.1f}″"
    a2.set_xticks(x, labs)
    a2.set_ylim(0, 105)
    a2.set_ylabel("Share confirmed by Hubble (%)")
    a2.set_xlabel("Einstein radius predicted from the spectrum alone")
    a2.grid(axis="x", visible=False)
    a2.legend(loc="upper left", fontsize=7.6)
    _title(a2, "Bigger predicted rings are confirmed more often",
           f"prediction: isothermal sphere with the SDSS velocity dispersion;\n"
           f"the fibre radius is {FIBRE_RADIUS_ARCSEC}″")
    ps.caption(fig, "Data: SDSS DR17 spectra of all 131 SLACS candidates imaged with HST; grades, "
                    "redshifts and dispersions from Bolton et al. (2008). Candidates without a "
                    "reliable dispersion (low S/N, or two galaxies in the fibre) are left out of "
                    "the right panel.")
    return ps.save(fig, FIGURE_DIR / "04_how_found")


# --------------------------------------------------------------------------
# 6. Isothermal: lensing speed against the stars' speed
# --------------------------------------------------------------------------

def figure_isothermal(t, iso):
    s = t[t["name"].isin(iso["names"])]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.9), gridspec_kw={"wspace": 0.27})
    ey = np.maximum(s["e_sigma_star"], 0.05 * s["sigma_star"]) * s["sigma_e8"] / s["sigma_star"]
    lim = [160, 420]
    a1.plot(lim, lim, color=ps.AXIS, lw=1.0, zorder=1)
    _label(a1, "equal", (395, 395), (395, 395), leader=False, ha="center", va="bottom")
    f, ef = iso["f_e8"][0], iso["f_e8"][1]
    a1.plot(lim, np.array(lim) * f, color=ps.C_DARK, lw=1.4, zorder=2)
    _label(a1, f"our fit: σ★ = {f:.3f} σ_SIE", (330, 330 * f), (175, 395))
    _points(a1, s["sigma_sie"], s["sigma_e8"], yerr=ey)
    a1.set_xlim(lim)
    a1.set_ylim(lim)
    a1.set_aspect("equal")
    a1.set_xlabel("σ_SIE: isothermal velocity dispersion from the ring (km/s)")
    a1.set_ylabel("σ★: stars' velocity dispersion (km/s)")
    _title(a1, "The ring knows how fast the stars move",
           f"{iso['n']} early-type lenses; σ★ corrected to a radius of R_e/8\n"
           "(Jørgensen et al. 1995), errors floored at 5%")

    ratio = s["sigma_e8"] / s["sigma_sie"]
    x = s["theta_e"] / s["r_eff"]
    a2.axhspan(f - iso["f_e8"][2], f + iso["f_e8"][2], color=ps.C_DARK, alpha=0.10, lw=0)
    a2.axhline(f, color=ps.C_DARK, lw=1.4)
    a2.axhline(1.0, color=ps.AXIS, lw=1.0)
    _points(a2, x, ratio, yerr=ey / s["sigma_sie"])
    a2.set_xlabel("Einstein radius / effective radius of the light")
    a2.set_ylabel("σ★ / σ_SIE")
    a2.set_ylim(0.6, 1.45)
    _label(a2, f"mean {f:.3f} ± {ef:.3f}; band: ±1 rms ({iso['f_e8'][2]:.3f})", (x.max(), f),
           (x.max() * 0.62, 1.38), leader=False)
    _title(a2, "…at every radius the rings probe",
           f"no trend with how far out the ring lies (Spearman ρ = "
           f"{iso['spearman_ratio_vs_rein_over_reff']:+.2f})".replace("-", "−"))
    ps.caption(fig, f"{PROV_SLACS} σ_SIE = c·√(θ_E D_s / 4π D_ls). Published (Bolton et al. 2008b): "
                    "σ★ = (1.019 ± 0.008) σ_SIE for the same 53 lenses; ours "
                    f"{f:.3f} ± {ef:.3f}. Uncorrected SDSS fibre values give "
                    f"{iso['f_sdss'][0]:.3f} ± {iso['f_sdss'][1]:.3f} (published 0.948 ± 0.008).")
    return ps.save(fig, FIGURE_DIR / "06_isothermal")


# --------------------------------------------------------------------------
# 7. Density slopes from lensing + dynamics
# --------------------------------------------------------------------------

def figure_density_slopes(t, slopes, auger, variants):
    tab = slopes["table"].merge(t[["name", "theta_e", "r_eff"]], on="name")
    pop = slopes["population"]
    fig = plt.figure(figsize=(12.4, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 0.9], wspace=0.33)
    a1, a2, a3 = (fig.add_subplot(gs[0, i]) for i in range(3))

    x = tab["theta_e"] / tab["r_eff"]
    a1.axhspan(pop["mu"] - pop["scatter"], pop["mu"] + pop["scatter"], color=ps.C_DARK,
               alpha=0.10, lw=0)
    a1.axhline(pop["mu"], color=ps.C_DARK, lw=1.4)
    a1.axhline(2.0, color=ps.AXIS, lw=1.0, ls=DASHED)
    _points(a1, x, tab["gamma"], yerr=[tab["gamma_lo"], tab["gamma_hi"]])
    a1.set_ylim(1.1, 2.9)
    a1.set_xlabel("Einstein radius / effective radius")
    a1.set_ylabel("Density slope γ′   (ρ ∝ r^−γ′)")
    _label(a1, f"population mean {pop['mu']:.3f} ± {pop['mu_err']:.3f}\n"
               f"intrinsic scatter {pop['scatter']:.2f} (band)", (x.max(), pop["mu"]),
           (x.max() * 1.02, 1.28), leader=False, ha="right")
    _label(a1, "isothermal: γ′ = 2", (0.25, 2.0), (0.19, 1.32), ha="left")
    _title(a1, f"How steeply does the total mass fall? ({slopes['n']} lenses)",
           "lensing fixes the mass inside the ring; the stars' speed through the SDSS fibre\n"
           "fixes the slope (spherical Jeans equation, isotropic orbits)")

    common = tab.merge(auger, on="name", suffixes=("", "_pub"))
    lim = [1.4, 2.8]
    a2.plot(lim, lim, color=ps.AXIS, lw=1.0)
    _points(a2, common["gamma_pub"], common["gamma"], ms=4.2)
    d = common["gamma"] - common["gamma_pub"]
    a2.set_xlim(lim)
    a2.set_ylim(lim)
    a2.set_aspect("equal")
    a2.set_xlabel("γ′ published (Auger et al. 2010)")
    a2.set_ylabel("γ′ ours")
    _title(a2, "Lens by lens, against SLACS X",
           f"{len(common)} lenses in both: median difference {d.median():+.3f},\nrms {np.sqrt(np.mean(d**2)):.3f}")

    # Population mean: ours (baseline and variants) against the published ones.
    rows = [("ours, isotropic, 1.5″ seeing", pop["mu"], pop["mu_lo"], pop["mu_hi"], ps.C_DARK)]
    for lab, v in variants:
        p = v["population"]
        rows.append((lab, p["mu"], p["mu_lo"], p["mu_hi"], ps.MUTED))
    rows += [("Auger et al. 2010 (73 lenses)", 2.078, 0.027, 0.027, ps.C_PUBLISHED),
             ("Koopmans et al. 2009 (58)", 2.085, 0.018, 0.025, ps.C_PUBLISHED)]
    for i, (lab, mu, lo, hi, col) in enumerate(rows):
        a3.errorbar(mu, i, xerr=[[lo], [hi]], fmt="o", color=col, ecolor=col, elinewidth=1.2,
                    capsize=0, markersize=5, markeredgecolor=ps.SURFACE)
    a3.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=7.6, color=ps.INK_SECONDARY)
    a3.tick_params(axis="y", length=0)
    a3.yaxis.tick_right()
    a3.invert_yaxis()
    a3.axvline(2.0, color=ps.AXIS, lw=1.0, ls=DASHED)
    a3.set_xlabel("Population mean γ′")
    a3.grid(axis="y", visible=False)
    a3.spines["left"].set_visible(False)
    _title(a3, "The mean, and how it moves",
           "grey: our variants (orbits, seeing)")
    ps.caption(fig, f"{PROV_SLACS} Stars traced by a Hernquist profile with the measured effective "
                    "radius; dispersion averaged over the 3″ fibre after 1.5″ seeing (Koopmans et al. "
                    "2006, 2009). Population mean and intrinsic scatter from a hierarchical fit to "
                    "every lens's full likelihood.")
    return ps.save(fig, FIGURE_DIR / "07_density_slopes")


# --------------------------------------------------------------------------
# 8. How much of the mass inside the ring is dark?
# --------------------------------------------------------------------------

def figure_dark_matter(t, dm):
    s = dm["table"]
    from matplotlib.ticker import FixedLocator, NullFormatter

    fig = plt.figure(figsize=(13.4, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.0], wspace=0.36)
    a1, a2, a3 = (fig.add_subplot(gs[0, i]) for i in range(3))

    lim = [2e10, 7e11]
    a1.fill_between(lim, lim, [lim[1] * 10] * 2, color=ps.C_CONTEXT, alpha=0.5, lw=0, zorder=0)
    a1.plot(lim, lim, color=ps.AXIS, lw=1.0, zorder=1)
    _label(a1, "more stars than mass:\nimpossible", (4e10, 3e11), (4e10, 3e11), leader=False,
           ha="center")
    a1.scatter(s["m_ein"], s["mstar_ein_salp"], s=20, color=ps.C_STARS_HEAVY, edgecolor=ps.SURFACE,
               lw=0.6, zorder=3, label="Salpeter IMF")
    a1.scatter(s["m_ein"], s["mstar_ein_chab"], s=20, color=ps.C_STARS, edgecolor=ps.SURFACE,
               lw=0.6, zorder=4, label="Chabrier IMF")
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlim(lim)
    a1.set_ylim(lim)
    a1.set_aspect("equal")
    ticks = [3e10, 1e11, 3e11]
    for axis in (a1.xaxis, a1.yaxis):
        axis.set_major_locator(FixedLocator(ticks))
        axis.set_major_formatter(lambda v, _: {3e10: "3×10¹⁰", 1e11: "10¹¹", 3e11: "3×10¹¹"}.get(v, ""))
        axis.set_minor_formatter(NullFormatter())
    a1.set_xlabel("Total mass inside the ring, from lensing (M☉)")
    a1.set_ylabel("Mass in stars inside the ring (M☉)")
    a1.legend(loc="lower right", fontsize=7.6)
    _title(a1, f"Stars against the lensing mass ({dm['n']} lenses)",
           "stellar mass from Hubble colours (Auger et al. 2009),\nscaled to the light inside the ring")

    bins = np.linspace(-0.4, 1.0, 29)
    for imf, col, lab in (("chab", ps.C_STARS, "Chabrier"), ("salp", ps.C_STARS_HEAVY, "Salpeter")):
        f = 1 - s[f"fstar_ein_{imf}_ours"]
        a2.hist(f, bins=bins, color=col, alpha=0.22, histtype="stepfilled", lw=0)
        a2.hist(f, bins=bins, color=col, histtype="step", lw=1.8,
                label=f"{lab}: mean {f.mean():.2f}")
    a2.axvline(0, color=ps.INK_SECONDARY, lw=1.0)
    a2.set_xlabel("Dark share of the mass inside the ring")
    a2.set_ylabel("Lenses")
    a2.legend(loc="upper left", fontsize=7.6)
    _title(a2, "A quarter to three fifths is dark",
           f"depending on the stars' mass function;\n{dm['salp']['n_fstar_above_1']} of {dm['n']} "
           "lenses go negative with Salpeter")

    fdm = 1 - s["fstar_ein_chab_ours"]
    _points(a3, s["sigma_star"], fdm, ms=4.2, color=ps.C_DARK)
    a3.axhline(0, color=ps.AXIS, lw=1.0)
    a3.set_ylim(-0.05, 0.9)
    a3.set_xlabel("Stars' velocity dispersion (km/s)")
    a3.set_ylabel("Dark share inside the ring (Chabrier)")
    _title(a3, "Heavier galaxies are darker",
           f"Spearman ρ = {dm['spearman_fdm_vs_sigma']:+.2f} against σ★ (p = {dm['p_fdm_vs_sigma']:.2g}),\n"
           f"{dm['spearman_fdm_vs_mstar']:+.2f} against stellar mass (p = {dm['p_fdm_vs_mstar']:.2g})")
    ps.caption(fig, f"{PROV_SLACS} Stellar masses: Bayesian population-synthesis fits to HST "
                    "photometry (Bruzual & Charlot 2003) for two initial mass functions; Salpeter "
                    f"masses are {dm['salp_over_chab']:.2f}× Chabrier. Light inside the ring from the "
                    "de Vaucouleurs fit to each galaxy. H0 = 70 km/s/Mpc, Ωm = 0.3, as in SLACS.")
    return ps.save(fig, FIGURE_DIR / "08_dark_matter")


# --------------------------------------------------------------------------
# 9. Stars + dark halo: the conspiracy, and how heavy the stars are
# --------------------------------------------------------------------------

def figure_conspiracy(t, dc, example_row, prof, treu, slopes):
    tab = dc["table"]
    pop = dc["population"]
    fig = plt.figure(figsize=(13.4, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.25, 0.85], wspace=0.36)
    a1, a2, a3 = (fig.add_subplot(gs[0, i]) for i in range(3))

    x = prof["r_over_re"]
    norm = (prof["rho"] * x**2)[np.argmin(np.abs(x - prof["theta_e_over_re"]))]
    a1.plot(x, prof["rho_star"] * x**2 / norm, color=ps.C_STARS, lw=2.0)
    a1.plot(x, prof["rho_dm"] * x**2 / norm, color=ps.C_DARK, lw=2.0)
    a1.plot(x, prof["rho"] * x**2 / norm, color=ps.C_OBSERVED, lw=2.2)
    a1.axvspan(0.2, 1.3, color=ps.C_CONTEXT, alpha=0.35, lw=0)
    a1.axvline(prof["theta_e_over_re"], color=ps.MUTED, lw=0.8, ls=DOTTED)
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_ylim(0.03, 5)
    a1.set_xlim(0.03, 10)
    i_lab = np.argmin(np.abs(x - 3.0))
    _label(a1, "stars", (x[i_lab], prof["rho_star"][i_lab] * x[i_lab] ** 2 / norm),
           (4.5, 0.06), color=ps.INK_SECONDARY)
    _label(a1, "dark halo (NFW)", (x[i_lab], prof["rho_dm"][i_lab] * x[i_lab] ** 2 / norm),
           (2.2, 2.5))
    _label(a1, f"total: mean slope {prof['mean_slope']:.2f}", (0.5, 1.0), (0.04, 2.6))
    _label(a1, "ring", (prof["theta_e_over_re"], 0.04), (prof["theta_e_over_re"] * 1.1, 0.045),
           leader=False)
    a1.set_xlabel("Radius / effective radius")
    a1.set_ylabel("Density × r²   (flat = isothermal)")
    _title(a1, "Stars + dark halo ≈ isothermal",
           f"{example_row.name}, best stars + NFW model; shaded:\n0.2–1.3 R_e, where lensing and the fibre measure")

    e_lo, e_hi = tab["log_alpha_lo"], tab["log_alpha_hi"]
    a2.axhline(0.0, color=ps.C_STARS_HEAVY, lw=1.4)
    a2.axhline(-dc["log_salp_over_chab"], color=ps.C_STARS, lw=1.4)
    _label(a2, "Salpeter IMF", (152, 0.0), (152, 0.02), leader=False, ha="left", va="bottom")
    _label(a2, "Chabrier IMF", (152, -dc["log_salp_over_chab"]), (152, -dc["log_salp_over_chab"] - 0.02),
           leader=False, ha="left", va="top")
    common = tab.merge(treu, on="name", suffixes=("", "_pub"))
    a2.scatter(common["sigma_star"], common["log_alpha_pub"], s=24, facecolor="none",
               edgecolor=ps.C_PUBLISHED, lw=0.9, zorder=3, label="Treu et al. 2010")
    _points(a2, tab["sigma_star"], tab["log_alpha"], yerr=[e_lo, e_hi], ms=4.4, label="ours",
            color=ps.C_OBSERVED)
    xs = np.linspace(160, 410, 50)
    a2.plot(xs, dc["trend_at_250"] + dc["trend_slope"] * np.log10(xs / 250), color=ps.C_DARK,
            lw=1.4, zorder=4, label=f"our trend: slope {dc['trend_slope']:.2f} ± {dc['trend_slope_err']:.2f}")
    a2.set_ylim(-1.0, 0.6)
    a2.set_xlim(150, 415)
    a2.set_xlabel("Stars' velocity dispersion (km/s)")
    a2.set_ylabel("log α = log (M★ lensing+dynamics / M★ Salpeter)")
    a2.legend(loc="lower right", fontsize=7.5)
    _title(a2, f"How heavy are the stars? ({dc['n']} lenses)",
           f"population mean log α = {pop['mu']:+.3f} ± {pop['mu_err']:.3f}: Salpeter-like,\n"
           f"{(pop['mu'] + dc['log_salp_over_chab']) / pop['mu_err']:.0f}σ above Chabrier")

    d = common["log_alpha"] - common["log_alpha_pub"]
    lim = [-0.8, 0.4]
    a3.plot(lim, lim, color=ps.AXIS, lw=1.0)
    _points(a3, common["log_alpha_pub"], common["log_alpha"], ms=4.0)
    a3.set_xlim(lim)
    a3.set_ylim(lim)
    a3.set_aspect("equal")
    a3.set_xlabel("log α, Treu et al. 2010")
    a3.set_ylabel("log α, ours")
    _title(a3, "Lens by lens",
           f"{len(common)} in both: median difference\n{d.median():+.3f}, rms {np.sqrt(np.mean(d**2)):.3f}")
    ps.caption(fig, f"{PROV_SLACS} Model (Treu et al. 2010): Hernquist stars plus an NFW halo with a "
                    "30-kpc break radius; lensing fixes the sum inside the Einstein radius, the "
                    "fibre dispersion the split (isotropic orbits). α compares the stellar mass this "
                    "requires with a Salpeter-IMF population fit to the same galaxy.")
    return ps.save(fig, FIGURE_DIR / "09_conspiracy")


# --------------------------------------------------------------------------
# 1. The gallery
# --------------------------------------------------------------------------

def _ring_image(npz, *, stretch=None):
    """Lens-subtracted image in units of the noise, for display; outside the
    fit region the full-box light model (stage 1) is used."""
    data, var = npz["data"], npz["var"]
    light = np.where(npz["fit_mask"], npz["light"], npz["light_full"]) if "light_full" in npz else npz["light"]
    sig = np.sqrt(np.median(var[np.isfinite(var)]))
    resid = (data - light) / sig
    n = data.shape[0]
    yy, xx = np.mgrid[:n, :n] - n // 2
    resid[np.hypot(xx, yy) * 0.04 < 0.25] = 0.0  # the masked core: no light model there
    return np.arcsinh(resid / (stretch or 3.0))


def _crit_curve(ax, geo, color, ls="-", lw=1.0, alpha=0.95):
    import lensing as L
    cx, cy = L.sie_critical_curve(geo["b"], geo["q"], geo["phi"], geo["x0"], geo["y0"])
    ax.plot(cx, cy, color=color, lw=lw, ls=ls, alpha=alpha)


def _image_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def _scalebar(ax, half, color="white"):
    """1-arcsecond bar in the lower-left corner."""
    x0, y0 = -half * 0.88, -half * 0.86
    ax.plot([x0, x0 + 1.0], [y0, y0], color=color, lw=1.6, solid_capstyle="butt")
    ax.text(x0 + 0.5, y0 + 0.12 * half / 2.4, "1″", color=color, ha="center", va="bottom", fontsize=7)


def figure_gallery(t, names, results):
    import hstimages as H

    half_arcsec = 2.4
    half_pix = int(round(half_arcsec / 0.04))
    ncol_lens = 3
    nrow = int(np.ceil(len(names) / ncol_lens))
    fig, axes = plt.subplots(nrow, 2 * ncol_lens, figsize=(13.6, 2.42 * nrow + 0.9),
                             gridspec_kw={"wspace": 0.04, "hspace": 0.16, "top": 0.90,
                                          "bottom": 0.03})
    axes = np.atleast_2d(axes)
    row_of = t.set_index("name")
    for k, name in enumerate(names):
        r_, c_ = divmod(k, ncol_lens)
        a_img, a_ring = axes[r_, 2 * c_], axes[r_, 2 * c_ + 1]
        npz, res = H.load_model(name)
        info = row_of.loc[name]
        centre = tuple(int(v) for v in npz["centre_pix"])
        rgb, blue = H.composite(name, centre, half_pix)
        ext = [-half_arcsec, half_arcsec, -half_arcsec, half_arcsec]
        if rgb is not None:
            a_img.imshow(rgb, extent=ext, origin="lower", interpolation="lanczos")
        else:
            img = H.grey(name, centre, half_pix)
            peak = np.percentile(img, 99.9)
            a_img.imshow(np.arcsinh(img / (0.02 * peak)), extent=ext, origin="lower",
                         cmap="gray", vmin=0, interpolation="lanczos")
        ring = _ring_image(npz)
        h0 = ring.shape[0] // 2
        sub = ring[h0 - half_pix:h0 + half_pix, h0 - half_pix:h0 + half_pix]
        vmax = np.nanpercentile(sub[np.hypot(*np.meshgrid(np.arange(sub.shape[0]) - half_pix,
                                                           np.arange(sub.shape[0]) - half_pix)) * 0.04 < 2.2], 99.6)
        a_ring.imshow(sub, extent=ext, origin="lower", cmap=ps.CMAP_RING, vmin=0, vmax=vmax,
                      interpolation="lanczos")
        geo = H.model_geometry(npz)
        _crit_curve(a_ring, geo, "white", ls=(0, (3, 2)), lw=0.8, alpha=0.8)
        for a in (a_img, a_ring):
            _image_axes(a)
            # North up, east left, in every panel
        _scalebar(a_img, half_arcsec)
        import lensing as L
        m_e = float(L.einstein_mass(res["b"], info["d_l"], info["d_s"], info["d_ls"]))  # from our radius
        a_img.set_title(f"{name}", fontsize=8.6, pad=3, loc="left")
        a_ring.set_title(f"z = {info['z_lens']:.3f} → {info['z_source']:.3f}", fontsize=7.6, pad=3,
                         loc="right", fontweight="normal", color=ps.INK_SECONDARY)
        a_ring.text(0.03, 0.03, f"θ_E = {res['b']:.2f}″\nM = {m_e / 1e11:.1f} × 10¹¹ M☉",
                    transform=a_ring.transAxes, color="white", fontsize=6.8, va="bottom",
                    linespacing=1.2)
    for k in range(len(names), nrow * ncol_lens):
        r_, c_ = divmod(k, ncol_lens)
        axes[r_, 2 * c_].axis("off")
        axes[r_, 2 * c_ + 1].axis("off")
    fig.suptitle("Einstein rings: twelve SLACS lenses, and what is left when the lens galaxy is removed",
                 x=0.125, ha="left", y=0.985, fontsize=12)
    fig.text(0.125, 0.962, "Left of each pair: Hubble colour image (F814W red, a blue filter blue; "
             "yellow-orange = old stars of the lens, blue = young stars of the galaxy behind it). "
             "Right: F814W minus our model of the lens galaxy's light;\ndashed: our lens model's "
             "critical curve. Each panel 4.8″ across, north up, east left. θ_E: our Einstein radius; "
             "M: the mass inside it, from lensing.", fontsize=8.2, color=ps.INK_SECONDARY, va="top")
    ps.caption(fig, "Data: HST ACS/WFC (and WFC3/UVIS for some blue images) via MAST Hubble Advanced "
                    "Products mosaics; programs 10174, 10494, 10587, 10798, 10886 and others. Lens "
                    "models and lens-light subtraction: this project. Redshifts: Bolton et al. 2008; "
                    "Auger et al. 2009.", y=0.022)
    return ps.save(fig, FIGURE_DIR / "01_gallery")


# --------------------------------------------------------------------------
# 2. One lens taken apart
# --------------------------------------------------------------------------

def figure_anatomy(t, name):
    import hstimages as H
    import lensing as L

    npz, res = H.load_model(name)
    info = t.set_index("name").loc[name]
    geo = H.model_geometry(npz)
    data, light, model, var = npz["data"], npz["light"], npz["model"], npz["var"]
    sig = np.sqrt(np.median(var[np.isfinite(var)]))
    n = data.shape[0]
    half = n // 2 * 0.04
    ext = [-half, half, -half, half]
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 8.9), gridspec_kw={"wspace": 0.08, "hspace": 0.26,
                                                                       "top": 0.86})
    peak = np.percentile(data, 99.9)
    st = lambda im: np.arcsinh(im / (0.02 * peak))
    vmax = st(peak)
    fm = npz["fit_mask"]
    ring = np.where(fm, data - light, np.nan) / sig
    yy, xx = np.mgrid[:n, :n] - n // 2
    core = np.hypot(xx, yy) * 0.04 < 0.25  # masked in the fit: the light model is free there
    rr = np.hypot(xx, yy)
    in_fit = (rr <= rr[fm].max()) & ~core  # the fitted circle, holes (cosmic rays) filled
    light_disp = np.where(core, np.nan, npz["light_full"] if "light_full" in npz else light)
    rmax = np.nanpercentile(ring, 99.5)

    panels = [
        (axes[0, 0], st(data), dict(cmap="gray", vmin=0, vmax=vmax), "1. What Hubble sees",
         "F814W (near-infrared), lens galaxy plus arcs"),
        (axes[0, 1], st(light_disp), dict(cmap="gray", vmin=0, vmax=vmax), "2. Our model of the lens's light",
         "14 concentric elliptical Gaussians (centre masked)"),
        (axes[0, 2], np.arcsinh(ring / 3), dict(cmap=ps.CMAP_RING, vmin=0, vmax=np.arcsinh(rmax / 3)),
         "3. Image minus lens light", "the lensed galaxy; dashed: critical curve"),
        (axes[1, 0], np.arcsinh(np.where(in_fit, model - light, np.nan) / sig / 3),
         dict(cmap=ps.CMAP_RING, vmin=0, vmax=np.arcsinh(rmax / 3)), "4. Our lens model's prediction",
         f"SIE lens + one Sérsic source, θ_E = {res['b']:.3f}″"),
        (axes[1, 1], np.where(fm, (data - model) / np.sqrt(var), np.nan),
         dict(cmap=ps.CMAP_RESIDUAL, vmin=-5, vmax=5), "5. What the model misses",
         f"(data − model)/noise, χ²/pixel = {res['chi2_red']:.2f}"),
    ]
    import copy
    for ax, img, kw, title, sub in panels:
        cmap = copy.copy(plt.get_cmap(kw.pop("cmap")))
        cmap.set_bad(ps.GRID if cmap.name == ps.CMAP_RESIDUAL else "#141414")
        ax.imshow(img, extent=ext, origin="lower", interpolation="nearest", cmap=cmap, **kw)
        _image_axes(ax)
        ax.set_title(title, fontsize=10, loc="left", pad=16)
        ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=7.8, color=ps.INK_SECONDARY, va="bottom")
    for ax in (axes[0, 2], axes[1, 0]):
        _crit_curve(ax, geo, "white", ls=(0, (3, 2)), lw=0.9)
    _scalebar(axes[0, 0], half)

    ax = axes[1, 2]
    img, sext = H.unlens(npz, half_width=0.7, n=70)
    ax.imshow(img, extent=sext, origin="lower", cmap=ps.CMAP_RING, vmin=0,
              vmax=np.nanpercentile(img, 99), interpolation="nearest")
    (cx, cy), (ux, uy) = L.sie_caustics(geo["b"], geo["q"], geo["phi"], 0, 0)
    ax.plot(cx, cy, color=ps.C_DARK, lw=1.2)
    ax.plot(ux, uy, color=ps.C_DARK, lw=1.0, ls=(0, (3, 2)))
    ax.set_xlim(sext[:2])
    ax.set_ylim(sext[2:])
    _image_axes(ax)
    ax.set_facecolor("#0b0b0b")
    mu = H.magnification(npz)
    r_kpc = res["source_r"] * float(L.kpc_per_arcsec(L.distances(info["z_lens"], info["z_source"],
                                                                     70.0, 0.3)["d_s"]))
    ax.set_title("6. The galaxy behind, un-lensed", fontsize=10, loc="left", pad=16)
    ax.text(0, 1.012, f"traced back through the lens; magnified {mu:.0f}×, half-light radius "
            f"{r_kpc:.1f} kpc", transform=ax.transAxes, fontsize=7.8, color=ps.INK_SECONDARY,
            va="bottom")
    ax.text(0.97, 0.04, "orange: caustics", transform=ax.transAxes, fontsize=7.2, color="white",
            ha="right")
    x0 = sext[0] * 0.88
    ax.plot([x0, x0 + 0.2], [sext[2] * 0.86] * 2, color="white", lw=1.6)
    ax.text(x0 + 0.1, sext[2] * 0.80, "0.2″", color="white", ha="center", fontsize=7)

    m_ours = float(L.einstein_mass(res["b"], info["d_l"], info["d_s"], info["d_ls"]))
    s_ours = float(L.sigma_sie(res["b"], info["d_s"], info["d_ls"]))
    fig.suptitle(f"{name}: one Einstein ring, taken apart", x=0.125, ha="left", y=0.99, fontsize=13)
    fig.text(0.125, 0.958,
             f"Lens galaxy at z = {info['z_lens']:.4f}, background galaxy at z = {info['z_source']:.4f}. "
             f"Our Einstein radius {res['b']:.3f}″ ({res['b'] * info['kpc_per_arcsec']:.2f} kpc); "
             f"published {info['theta_e']:.2f}″ (Bolton et al. 2008).\nMass inside it: "
             f"{m_ours / 1e11:.2f} × 10¹¹ M☉. Isothermal speed σ_SIE = {s_ours:.0f} km/s; "
             f"the stars' measured speed {info['sigma_star']:.0f} ± {info['e_sigma_star']:.0f} km/s.",
             fontsize=8.6, color=ps.INK_SECONDARY, va="top", linespacing=1.4)
    ps.caption(fig, "Data: HST ACS/WFC F814W, MAST Hubble Advanced Products mosaic (0.04″ pixels). "
                    "North up, east left; 6.4″ across. Point-spread function approximated by a "
                    "Moffat profile (FWHM 0.10″). Fit: this project.", y=0.03)
    return ps.save(fig, FIGURE_DIR / "02_anatomy")


# --------------------------------------------------------------------------
# 5. Our Einstein radii against the published ones
# --------------------------------------------------------------------------

def figure_einstein_radii(ours, ours_pass1, shear):
    df = ours["table"]
    ok = (df["failed"] == "") & df["b_pub"].notna()
    d = df[ok].copy()
    d["frac"] = d["b_ours"] / d["b_pub"] - 1
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 5.0), gridspec_kw={"wspace": 0.27,
                                                                         "width_ratios": [1, 1.1]})
    lim = [0.35, 2.45]  # wide enough that no point, however far off, is hidden
    a1.plot(lim, lim, color=ps.AXIS, lw=1.0)
    a1.fill_between(lim, np.array(lim) * 0.98, np.array(lim) * 1.02, color=ps.C_PUBLISHED,
                    alpha=0.10, lw=0)
    v = d["b_pub_from"] == "SLACS V"
    fl = d["flag"] != ""
    _points(a1, d.loc[v & d["ring"] & ~fl, "b_pub"], d.loc[v & d["ring"] & ~fl, "b_ours"], ms=5.0,
            label="full or partial ring (SLACS 'ring subset')")
    a1.scatter(d.loc[v & ~d["ring"] & ~fl, "b_pub"], d.loc[v & ~d["ring"] & ~fl, "b_ours"], s=26,
               facecolor=ps.SURFACE, edgecolor=ps.C_OBSERVED, lw=1.0, zorder=6,
               label="arcs or image pairs")
    a1.scatter(d.loc[~v & ~fl, "b_pub"], d.loc[~v & ~fl, "b_ours"], s=30, marker="D", facecolor="none",
               edgecolor=ps.C_PUBLISHED, lw=1.0, zorder=6, label="published in kpc only (SLACS IX)")
    if fl.any():
        a1.scatter(d.loc[fl, "b_pub"], d.loc[fl, "b_ours"], s=30,
                   marker="x", color=ps.MUTED, lw=1.1, zorder=5,
                   label="fit flagged by its own diagnostics")
    for row in d[(d["frac"].abs() > 0.10) & (d["flag"] == "")].itertuples():
        dy = 0.10 if row.b_ours > row.b_pub else -0.14
        if row.b_ours + dy > lim[1] - 0.08:  # keep labels near the top inside the axes
            dy = -0.16
        _label(a1, row.name, (row.b_pub, row.b_ours), (row.b_pub + 0.12, row.b_ours + dy), fontsize=6.8)
    a1.set_xlim(lim)
    a1.set_ylim(lim)
    a1.set_aspect("equal")
    a1.set_xlabel("Published Einstein radius (″)")
    a1.set_ylabel("Our Einstein radius (″)")
    a1.legend(loc="upper left", fontsize=7.4, bbox_to_anchor=(0.0, 1.0))
    _title(a1, f"{ours['n_compared']} Einstein radii re-measured from the images",
           "violet band: ±2%, the uncertainty SLACS adopted")

    bins = np.arange(-0.16, 0.161, 0.01)
    a2.hist(d["frac"].clip(-0.155, 0.155) * 100, bins=bins * 100, color=ps.C_OBSERVED,
            edgecolor=ps.SURFACE, lw=0.6)
    a2.axvspan(-2, 2, color=ps.C_PUBLISHED, alpha=0.10, lw=0)
    a2.axvline(0, color=ps.INK_SECONDARY, lw=1.0)
    a2.set_xlabel("Ours − published (%)   (outermost bins hold everything beyond ±15%)")
    a2.set_ylabel("Lenses")
    lines = [f"all {ours['n_compared']}: median {100 * ours['median_frac_diff']:+.1f}%, robust scatter "
             f"{100 * ours['robust_scatter']:.1f}%; {ours['n_within_2pct']} within 2%, "
             f"{ours['n_within_5pct']} within 5%, {ours['n_off_by_10pct']} off by >10%"]
    if ours.get("n_unflagged"):
        lines.append(f"{ours['n_unflagged']} unflagged: median {100 * ours['unflagged_median_frac_diff']:+.1f}%, "
                     f"scatter {100 * ours['unflagged_robust_scatter']:.1f}%, "
                     f"{ours['unflagged_n_off_by_10pct']} off by >10%")
    if ours_pass1:
        lines.append(f"first blind pass ({ours_pass1['n_compared']}): {ours_pass1['n_within_5pct']} within 5%, "
                     f"{ours_pass1['n_off_by_10pct']} off by >10%")
    if shear:
        lines.append(f"with external shear ({shear['n_compared']} ring-subset lenses): median "
                     f"{100 * shear['median_frac_diff']:+.1f}%, scatter {100 * shear['robust_scatter']:.1f}%")
    _title(a2, "How far off are we?", "\n".join(lines))
    ps.caption(fig, "Published: SIE models of Bolton et al. (2008) (circles) and Auger et al. (2009) "
                    "(diamonds, converted from kpc). Ours: SIE lens + Sérsic source + multi-Gaussian "
                    "lens light, fitted to the F814W image with no knowledge of the published model; "
                    "started from the radius where the arcs, traced back to the source plane, "
                    "converge.")
    return ps.save(fig, FIGURE_DIR / "05_einstein_radii")


# --------------------------------------------------------------------------
# 0. Our numbers against the literature
# --------------------------------------------------------------------------

def figure_comparison(rows):
    labels = [r[0] for r in rows]
    diffs = np.array([100.0 * (r[1] / r[2] - 1.0) for r in rows])
    unc = [100.0 * r[3] if r[3] else None for r in rows]
    fig, ax = plt.subplots(figsize=(9.8, 0.44 * len(rows) + 2.3))
    ypos = np.arange(len(rows))
    for y, u in zip(ypos, unc):
        if u:
            ax.add_patch(plt.Rectangle((-u, y - 0.32), 2 * u, 0.64, facecolor=ps.C_PUBLISHED,
                                       alpha=0.12, edgecolor="none", zorder=1))
    ax.barh(ypos, diffs, height=0.36, color=ps.C_DARK, zorder=3)
    ax.axvline(0.0, color=ps.INK_SECONDARY, linewidth=1.0, zorder=4)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8.5, color=ps.INK_SECONDARY)
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    lim = max(16.0, np.nanmax(np.abs(diffs)) * 1.35)
    ax.set_xlim(-lim, lim)
    for y, v, r in zip(ypos, diffs, rows):
        off = 5 if v >= 0 else -5
        txt = f"{v:+.1f}%  ({r[1]:.4g} vs {r[2]:.4g})".replace("-", "−")
        ax.annotate(txt, xy=(v, y), xytext=(off, 0), textcoords="offset points", va="center",
                    ha="left" if v >= 0 else "right", fontsize=7.6, color=ps.INK_SECONDARY)
    ax.set_xlabel("Difference from the published value (%)")
    _title(ax, "Our numbers against the literature",
           "violet band: the published uncertainty, where one is quoted")
    ps.caption(fig, "Published: Bolton et al. 2008a,b (σ ratios, Einstein radii, ±2% adopted error); "
                    "Auger et al. 2009 (masses and stellar shares, their per-lens table averaged over "
                    "the same lenses); Auger et al. 2010 (γ′, dark share inside R_e/2); Treu et al. "
                    "2010 (α). For α and γ′ the comparison is of population means.")
    return ps.save(fig, FIGURE_DIR / "00_comparison")
