"""Every figure, from results.json plus the artefacts analysis.py returns."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.table import Table

import binary
import plotstyle as ps
from config import FIGURE_DIR, GAIA_DIR, R0_KPC, V_TAN_MIN
from literature import CENSUS, REFERENCE

ps.use_style()

CLASS_STYLE = {
    "D6": dict(color=ps.BLUE, marker="o", label="D⁶ stars (surviving donors)"),
    "LP40-365": dict(color=ps.ORANGE, marker="s", label="LP 40-365 stars (Type Iax remnants?)"),
    "He-star donor": dict(color=ps.AQUA, marker="^", label="US 708 (helium-star donor)"),
    "other": dict(color=ps.AQUA, marker="D", label="J1240+6710 (oxygen white dwarf)"),
}
CLS = {r.name: r.cls for r in CENSUS}


def _out(name):
    return FIGURE_DIR / name


def _background():
    t = Table.read(GAIA_DIR / "background_100pc.fits")
    ok = np.isfinite(t["bp_rp"]) & (t["ruwe"] < 1.4)
    g = np.asarray(t["phot_g_mean_mag"])[ok] + 5 * np.log10(np.asarray(t["parallax"])[ok] / 100)
    return np.asarray(t["bp_rp"])[ok], g


def _hr_background(ax, extent=(-0.9, 3.6, 17.5, -1.0)):
    c, m = _background()
    ax.hexbin(c, m, gridsize=(160, 120), extent=(extent[0], extent[1], extent[3], extent[2]),
              bins="log", cmap="ramp_context", mincnt=1, linewidths=0, rasterized=True, zorder=0)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("colour BP − RP (mag)")
    ax.set_ylabel("absolute magnitude M_G (mag)")


# --------------------------------------------------------------------------
# 1. The search: candidates in the HR diagram, and the mirror test
# --------------------------------------------------------------------------

def fig_search(results, art):
    import search as S
    s = art["search"]
    pos = s[np.asarray(s["plx_snr"]) > 3]
    clean = S.quality_masks(pos)["all of these"]
    rob = art["robust"]
    st = results["search"]["shift_test"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw=dict(width_ratios=[1.1, 1]))
    ax = axes[0]
    _hr_background(ax)
    ax.scatter(np.asarray(pos["bp_rp"], float)[clean], np.asarray(pos["abs_M_G"], float)[clean], s=3,
               color=ps.MUTED, alpha=0.45, lw=0, rasterized=True,
               label=f"> 600 km/s at face value ({clean.sum():,}, clean astrometry)")
    unb = np.asarray(rob["v_gc_min"]) > np.asarray(rob["v_esc"])
    ax.scatter(np.asarray(rob["bp_rp"], float), np.asarray(rob["abs_M_G_corr"], float), s=22, color=ps.BLUE,
               edgecolor=ps.SURFACE, lw=0.5, zorder=4,
               label=f"robust: parallax would have to be > 3σ off ({len(rob)})")
    ax.scatter(np.asarray(rob["bp_rp"], float)[unb], np.asarray(rob["abs_M_G_corr"], float)[unb], s=110,
               facecolor="none", edgecolor=ps.ORANGE, lw=1.6, zorder=5,
               label=f"unbound whatever its radial velocity ({unb.sum()})")
    for r in rob[unb]:
        nm = str(r["known_as"]).split(" (")[0] or str(r["source_id"])
        ax.annotate(nm, (r["bp_rp"], r["abs_M_G_corr"]), xytext=(8, -3), textcoords="offset points",
                    fontsize=7, color=ps.INK)
    ax.set_title("Faster than 600 km/s across the sky?")
    leg = ax.legend(loc="upper right", fontsize=7, markerscale=1.0, frameon=True)
    leg.get_frame().set(facecolor=ps.SURFACE, edgecolor=ps.GRID, alpha=0.92)

    ax = axes[1]
    styles = {"catalogue errors": (ps.MUTED, "--", "Gaia's quoted errors"),
              "catalogue errors, clean astrometry": (ps.INK_SECONDARY, "-", "quoted errors, clean astrometry"),
              "zero point + external errors, clean": (ps.BLUE, "-", "zero point + realistic errors, clean")}
    for key, (col, ls, lbl) in styles.items():
        c = st["counts"][key]
        k = np.array([float(x) for x in c])
        n = np.array(list(c.values()), float)
        ax.plot(k, np.maximum(n, 0.6), color=col, ls=ls, marker="o", ms=3.5, lw=1.6, label=lbl)
    ax.axvline(3, color=ps.GRID, lw=1)
    ax.set_yscale("log")
    ax.set_ylim(0.6, 3e4)
    from matplotlib.ticker import FixedLocator
    ax.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels(["1", "10", "100", "1,000", "10,000"])
    ax.set_xlabel("how far (σ) the parallax would have to be off for the star to be slower than 600 km/s")
    ax.set_ylabel("candidates needing more than that")
    ax.set_title("Most candidates are the noise tail of ordinary stars")
    neg = results["search"]["n_negative"]
    ax.text(0.97, 0.62, f"{st['frac_within_1sigma']:.0%} would drop below 600 km/s\n"
            f"with a parallax just 1σ larger.\n\nOnly {neg} sources show the mirror image:\n"
            f"a parallax more than 3σ negative.", transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.legend(loc="upper right", fontsize=7)
    ps.caption(fig, "Gaia DR3 (ARI Heidelberg mirror): every source with pm > 20 mas/yr, G < 20.5, parallax > 3σ and "
               "4.74 μ/ϖ > 600 km/s. Clean astrometry: RUWE < 1.4, clean image parameters, Rybizki et al. (2022) "
               "fidelity > 0.5, ≥ 10 visibility periods. Realistic errors: Lindegren et al. (2021) zero point, "
               "errors scaled as Maíz Apellániz et al. (2021). Unbound: the slowest Galactocentric speed allowed by "
               "any radial velocity exceeds the escape speed of McMillan (2017)'s potential.")
    return ps.save(fig, _out("fig1_search"))


# --------------------------------------------------------------------------
# 2. Parallax errors: El-Badry et al.'s candidates, and Shen et al.'s DR2 ones
# --------------------------------------------------------------------------

def fig_parallax_errors(results):
    eb = results["search"]["elbadry"]["rows"]
    shen = results["search"]["shen"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw=dict(wspace=0.28))
    ax = axes[0]
    colors = {"runaway": ps.BLUE, "false positive": ps.ORANGE, "unclassified": ps.MUTED}
    marks = {"runaway": "o", "false positive": "s", "unclassified": "x"}
    for cls in ("false positive", "unclassified", "runaway"):
        rr = [r for r in eb if r["cls"] == cls and r["plx"] + r["plx_err"] > 0]
        x = [r["v_tan_lower"] for r in rr]
        y = [r["sigma_to_300"] for r in rr]
        lbl = {"runaway": "spectrum: thermonuclear runaway",
               "false positive": "spectrum: ordinary white dwarf / star",
               "unclassified": "bright neighbour / not observed"}[cls]
        ax.scatter(x, y, s=36, marker=marks[cls], color=colors[cls], label=lbl,
                   edgecolor=ps.SURFACE, lw=0.6, zorder=3)
        if cls == "runaway":
            offs = {"J1332-3541": (4, -10), "J1109+0001": (-6, -12), "J0927-6335": (-10, 7),
                    "J0546+0836": (5, 2)}
            for r in rr:
                nm = r["verdict"].split("(")[-1].rstrip(")")
                ax.annotate(nm, (r["v_tan_lower"], r["sigma_to_300"]), xytext=offs.get(nm, (5, 2)),
                            textcoords="offset points", fontsize=7, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("4.74 μ / (ϖ + σ): speed across the sky with the parallax 1σ high (km/s)")
    ax.set_ylabel("how far off (σ) the parallax must be for a 300 km/s halo star")
    e = results["search"]["elbadry"]
    same = "the same" if e["identical_to_table1"] else "not the published"
    ax.set_title(f"El-Badry et al.'s query rerun: {same} {e['n_returned']} sources")
    from matplotlib.ticker import FixedLocator, NullLocator
    ax.xaxis.set_major_locator(FixedLocator([400, 600, 1000, 2000, 3000]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels(["400", "600", "1,000", "2,000", "3,000"])
    ax.yaxis.set_major_locator(FixedLocator([1.5, 2, 3, 5, 10, 20, 30]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels(["1.5", "2", "3", "5", "10", "20", "30"])
    ax.legend(loc="upper left", fontsize=7)

    ax = axes[1]
    for r in shen:
        real = r["name"].startswith("D6")
        col = ps.BLUE if real else ps.ORANGE
        ax.annotate("", xy=(r["plx_dr3"], r["vtan_dr3"]), xytext=(r["plx_dr2"], r["vtan_dr2"]),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.4, shrinkA=3, shrinkB=3))
        ax.scatter([r["plx_dr2"]], [r["vtan_dr2"]], s=18, color=col, alpha=0.45, zorder=3)
        ax.scatter([r["plx_dr3"]], [r["vtan_dr3"]], s=30, color=col, zorder=4,
                   edgecolor=ps.SURFACE, lw=0.6)
        sh = r["shift_sigma"] if abs(r["shift_sigma"]) >= 0.05 else 0.0
        note = f"{r['name']}: {sh:+.1f}σ, RUWE {r['ruwe_dr3']:.2f}"
        off = {"D6-2": (6, 6), "O4": (6, -11), "D6-3": (6, 3), "D6-1": (6, -6)}.get(r["name"], (6, -3))
        ax.annotate(note, (r["plx_dr3"], r["vtan_dr3"]), xytext=off, textcoords="offset points",
                    fontsize=7, color=ps.INK_SECONDARY)
    ax.axhline(V_TAN_MIN, color=ps.MUTED, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("parallax ϖ (mas)")
    ax.set_ylabel("face-value speed across the sky (km/s)")
    ax.set_title("Shen et al.'s DR2 candidates, remeasured by DR3")
    ax.plot([], [], color=ps.BLUE, label="D⁶ stars (spectra confirm)")
    ax.plot([], [], color=ps.ORANGE, label="O1–O4 (spectra: ordinary stars)")
    ax.legend(loc="upper right", fontsize=7)
    ps.caption(fig, "Left: El-Badry et al.'s ADQL query rerun on Gaia DR3 (two sources with negative ϖ + σ not shown); "
               "verdicts from their Table 1 spectra. Right: arrows run from Gaia DR2 (22 months) to DR3 (34 months); "
               "shift in units of the DR2 error. O3's DR2 motion came from a blend with a star 1.8 mag brighter 0.55″ away.")
    return ps.save(fig, _out("fig2_parallax_errors"))


# --------------------------------------------------------------------------
# 3. HR diagram of the known runaways
# --------------------------------------------------------------------------

#: Label positions (data coordinates) in the zoomed HR panel, where the
#: LP 40-365 stars and the cool D6 stars crowd together.
HR_LABELS = {
    "D6-1": (0.52, 5.45), "J0905+2510": (0.52, 5.85), "D6-3": (0.52, 6.25),
    "LP 40-365": (0.52, 6.65), "J1311-1846": (0.52, 7.05), "D6-2": (0.52, 7.45),
    "J1603-6613": (-0.10, 5.55), "US 708": (-0.62, 3.75), "J1022-3414": (-0.22, 3.85),
}


def _hr_points(ax, kin, labels=True, label_pos=None):
    for name, k in kin.items():
        st = CLASS_STYLE[k["cls"]]
        y, lo, hi = k["M_G0"]["median"], k["M_G0"]["minus"], k["M_G0"]["plus"]
        ax.errorbar(k["bp_rp0"], y, yerr=[[lo], [hi]], fmt=st["marker"], color=st["color"],
                    ms=7, mec=ps.SURFACE, mew=0.8, elinewidth=1.0, capsize=0, zorder=5)
        if not labels:
            continue
        pos = (label_pos or {}).get(name)
        if pos is None:
            ax.annotate(name, (k["bp_rp0"], y), xytext=(7, 2), textcoords="offset points", fontsize=7,
                        color=ps.INK_SECONDARY, zorder=6)
        else:
            ax.annotate(name, (k["bp_rp0"], y), xytext=pos, textcoords="data", fontsize=7,
                        color=ps.INK_SECONDARY, va="center", ha="left" if pos[0] > k["bp_rp0"] else "right",
                        arrowprops=dict(arrowstyle="-", color=ps.MUTED, lw=0.6, shrinkA=1, shrinkB=4), zorder=6)


def fig_hr(results):
    kin = results["kinematics"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.2), gridspec_kw=dict(width_ratios=[1, 1.15]))
    ax = axes[0]
    _hr_background(ax, extent=(-0.9, 3.0, 15.5, 1.0))
    _hr_points(ax, kin, labels=False)
    ax.set_xlim(-0.9, 2.6)
    ax.text(0.55, 13.9, "white dwarfs", color=ps.INK_SECONDARY, fontsize=8, rotation=-22)
    ax.text(1.45, 6.35, "main sequence", color=ps.INK_SECONDARY, fontsize=8, rotation=-35)
    zx, zy = (-0.78, 0.62), (9.6, 2.9)
    ax.add_patch(plt.Rectangle((zx[0], zy[1]), zx[1] - zx[0], zy[0] - zy[1], fill=False, ec=ps.INK_SECONDARY,
                               lw=0.8, ls="--"))
    for cls, st in CLASS_STYLE.items():
        ax.errorbar([], [], fmt=st["marker"], color=st["color"], label=st["label"], ms=6)
    ax.legend(loc="upper right", fontsize=7)
    ax.set_title("Between the white dwarfs and the main sequence")
    ax = axes[1]
    _hr_background(ax, extent=(-0.9, 3.0, 15.5, 1.0))
    _hr_points(ax, kin, labels=True, label_pos=HR_LABELS)
    ax.set_xlim(zx[0], zx[1] + 0.33)
    ax.set_ylim(*zy)
    ax.set_title("Zoom: the 17 runaways")
    ps.caption(fig, "Absolute magnitudes and dereddened colours: medians and 68% ranges from our distance "
               "posteriors (prior flat in absolute magnitude; US 708 uses its spectroscopic distance). "
               "SFD dust (Schlafly & Finkbeiner 2011). Grey: Gaia DR3 stars within 100 pc.")
    return ps.save(fig, _out("fig3_hr_diagram"))


# --------------------------------------------------------------------------
# 4. The velocity ladder
# --------------------------------------------------------------------------

def fig_velocities(results):
    kin = results["kinematics"]
    names = sorted(kin, key=lambda n: kin[n]["v_ej"]["median"])
    fig, ax = plt.subplots(figsize=(9.5, 6.6))
    y = np.arange(len(names))
    vesc = np.median([kin[n]["v_esc"]["median"] for n in names])
    ax.axvspan(0, vesc, color=ps.GRID, alpha=0.6, lw=0, zorder=0)
    ax.text(vesc / 2, -0.75, f"bound to the Galaxy\n(escape speed ≈ {vesc:.0f} km/s)",
            ha="center", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    for i, n in enumerate(names):
        k = kin[n]
        st = CLASS_STYLE[k["cls"]]
        v = k["v_ej"]
        ax.plot([v["p2"], v["p98"]], [i, i], color=st["color"], lw=1.2, alpha=0.45, solid_capstyle="round")
        ax.plot([v["median"] - v["minus"], v["median"] + v["plus"]], [i, i], color=st["color"], lw=4,
                solid_capstyle="round")
        ax.plot(v["median"], i, st["marker"], color=st["color"], ms=7, mec=ps.SURFACE, mew=1.0)
        ax.plot(abs(k["rv"]), i, "|", color=ps.INK, ms=9, mew=1.4)
        eb = k.get("eb23") or {}
        if "vej" in eb:
            ax.plot(eb["vej"][0], i, marker=7, color=ps.MAGENTA, ms=7, zorder=6)
    ax.set_yticks(y, [f"{n}" for n in names])
    ax.set_ylim(-1.6, len(names) - 0.5)
    ax.set_xlabel("launch speed relative to the rotating disc (km/s)")
    ax.set_xlim(0, 3100)
    ax.set_title("The fastest stars in the Galaxy")
    ax.plot([], [], "|", color=ps.INK, ms=9, mew=1.4, label="|radial velocity| (measured directly)")
    ax.plot([], [], marker=7, color=ps.MAGENTA, ls="none", label="El-Badry et al. (2023) median")
    for cls in ("D6", "LP40-365", "He-star donor", "other"):
        st = CLASS_STYLE[cls]
        ax.plot([], [], st["marker"], color=st["color"], label=st["label"])
    ax.legend(loc="lower right", fontsize=7)
    ps.caption(fig, "Thick bar: 68% range; thin: 95%. Launch speed = today's Galactocentric velocity minus the "
                    "disc's circular velocity where the straight-line path crosses the midplane (McMillan 2017 "
                    "potential); Gaia DR3 parallaxes with the Lindegren et al. (2021) zero point. The radial velocity "
                    "is heliocentric, so it bounds the speed relative to the Sun, not the launch speed (see US 708).")
    return ps.save(fig, _out("fig4_velocities"))


# --------------------------------------------------------------------------
# 5. Trajectories through the Galaxy
# --------------------------------------------------------------------------

def fig_trajectories(results, orbit_art):
    """Small multiples: each D6 star seen edge-on, in Galactocentric radius R
    and height z, 3 Myr back and forward, with where it left the disc."""
    kin = results["kinematics"]
    d6 = sorted([n for n in orbit_art if CLS[n] == "D6"], key=lambda n: kin[n]["d"]["median"])
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.6), sharex=True, sharey=True)
    for ax, n in zip(axes.ravel(), d6):
        a = orbit_art[n]
        ax.axhline(0, color=ps.AXIS, lw=0.8)
        for k in range(min(60, a["past"].shape[1])):
            ax.plot(np.hypot(a["past"][0, k], a["past"][1, k]), a["past"][2, k], color=ps.BLUE, lw=0.5,
                    alpha=0.3)
            ax.plot(np.hypot(a["future"][0, k], a["future"][1, k]), a["future"][2, k], color=ps.ORANGE,
                    lw=0.5, alpha=0.3)
        xc = a["xc"]
        if xc.shape[1]:
            ax.scatter(np.hypot(xc[0], xc[1]), xc[2], s=4, color=ps.INK, alpha=0.35, lw=0, zorder=4)
        R0 = np.hypot(a["past"][0, :, 0], a["past"][1, :, 0])
        ax.scatter(np.median(R0), np.median(a["past"][2, :, 0]), s=30, color=ps.INK, zorder=6,
                   edgecolor=ps.SURFACE, lw=0.8)
        ax.plot(R0_KPC, 0.0208, marker="*", color=ps.YELLOW, ms=11, mec=ps.INK, mew=0.5, zorder=7)
        t = kin[n]["t_disk"]
        tt = (f"left the disc {t['median']:.2f} Myr ago" if t["median"] > 0
              else f"heading for the disc (in {abs(t['median']):.2f} Myr)")
        v = kin[n]["v_ej"]
        ax.set_title(n, fontsize=10)
        ax.text(0.03, 0.04, f"{tt}\nlaunched at {v['median']:.0f} km/s", transform=ax.transAxes,
                fontsize=7, color=ps.INK_SECONDARY, va="bottom")
        ax.set_xlim(0, 22)
        ax.set_ylim(-7, 7)
    fig.supxlabel("distance from the Galactic Centre R (kpc)", fontsize=9, color=ps.INK_SECONDARY, y=0.04)
    for ax in axes[:, 0]:
        ax.set_ylabel("height above the disc z (kpc)")
    axes[0, 0].plot([], [], color=ps.BLUE, label="past 3 Myr")
    axes[0, 0].plot([], [], color=ps.ORANGE, label="next 3 Myr")
    axes[0, 0].scatter([], [], s=8, color=ps.INK, label="launch point (disc crossing)")
    axes[0, 0].plot([], [], "*", color=ps.YELLOW, mec=ps.INK, ms=9, ls="none", label="Sun")
    axes[0, 0].legend(loc="upper left", fontsize=6.5)
    fig.suptitle("The D⁶ stars leaving the Milky Way, seen edge-on", x=0.06, ha="left")
    ps.caption(fig, "Each thin line is one draw from the star's posterior (distance, proper motion, radial velocity) "
                    "integrated in McMillan (2017)'s potential with galpy; ordered by distance. Black dots: where draws "
                    "crossed the midplane in the last 15 Myr. D6-2 and J0927-6335 are moving toward the midplane, so their "
                    "paths never crossed it: they were launched below the plane (D6-2 from G70.0-21.5, ~0.3 kpc down), "
                    "as old thick-disc binaries can be. Interactive: figures/runaway_explorer.html.")
    return ps.save(fig, _out("fig5_trajectories"))


# --------------------------------------------------------------------------
# 6. D6-2 and the remnant G70.0-21.5
# --------------------------------------------------------------------------

def halpha_cutout(lc, bc, win, n=700):
    """Finkbeiner (2003) H-alpha (Rayleigh) on an l-b grid around (lc, bc)."""
    import healpy as hp
    from config import AUX_DIR
    path = AUX_DIR / "lambda_halpha_fwhm06_0512.fits"
    if not path.exists():
        return None
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = hp.read_map(path, nest=True)
    ls = np.linspace(lc - win, lc + win, n)
    bs = np.linspace(bc - 0.8 * win, bc + 0.8 * win, int(n * 0.8))
    L, B = np.meshgrid(ls, bs)
    pix = hp.ang2pix(512, np.radians(90 - B), np.radians(L % 360), nest=True)
    img = m[pix]
    img[img < -1e20] = np.nan
    return img, (ls[0], ls[-1], bs[0], bs[-1])


def fig_remnant(results, rem_art):
    res = results["remnants"]
    tr = rem_art["tracks"]["D6-2"]
    green = rem_art["green"]
    best = res["D6-2"]["best"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax = axes[0]
    j0 = int(np.flatnonzero(np.asarray(green["name"]) == best["name"])[0])
    lc, bc = float(green["l"][j0]) - 1.0, float(green["b"][j0]) - 0.5
    win = 9.5
    halpha = halpha_cutout(lc, bc, win)
    if halpha is not None:
        img, ext = halpha
        from matplotlib.colors import LogNorm
        lo, hi = np.nanpercentile(img, 2), np.nanpercentile(img, 99.5)
        ax.imshow(np.clip(img, lo, hi), extent=ext, origin="lower", cmap="ramp_context", aspect="auto",
                  norm=LogNorm(vmin=max(lo, 0.05), vmax=hi), zorder=0, interpolation="bilinear")
    for k in range(min(150, tr["l"].shape[0])):
        ax.plot(tr["l"][k], tr["b"][k], color=ps.BLUE, lw=0.5, alpha=0.12, zorder=2)
    ax.plot(tr["lm"], tr["bm"], color=ps.BLUE, lw=1.6, zorder=3, label="D6-2's explosion site, traced back")
    for t_ky in (0, 50, 100, 150, 200, 250, 300):
        j = int(np.argmin(np.abs(tr["t_yr"] - t_ky * 1e3)))
        ax.plot(tr["lm"][j], tr["bm"][j], "o", color=ps.BLUE, ms=3.5, zorder=4)
        if t_ky % 100 == 0:
            ax.annotate("today" if t_ky == 0 else f"{t_ky} kyr ago", (tr["lm"][j], tr["bm"][j]),
                        xytext=(4, -11), textcoords="offset points", fontsize=7, color=ps.INK,
                        zorder=5)
    near = [r for r in green if abs(((r["l"] - lc + 180) % 360) - 180) < win and abs(r["b"] - bc) < win]
    for r in near:
        # Green gives sizes, not orientations: drawn axis-aligned, with the
        # longitude extent stretched by 1/cos(b) so it spans the true angle
        w, h = r["maj_arcmin"] / 60 / np.cos(np.radians(r["b"])), r["min_arcmin"] / 60
        hit = r["name"] == best["name"]
        ax.add_patch(Ellipse((r["l"], r["b"]), w, h, fill=False, lw=1.6 if hit else 0.9,
                             ec=ps.ORANGE if hit else ps.MUTED, ls="-" if hit else "--", zorder=3))
        ax.annotate(r["name"], (r["l"] + w / 2, r["b"] + h / 2), fontsize=7.5,
                    color=ps.ORANGE if hit else ps.MUTED, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", fc=ps.SURFACE, ec="none", alpha=0.85))
    ax.set_xlim(lc + win, lc - win)
    ax.set_ylim(bc - win * 0.8, bc + win * 0.8)
    ax.set_xlabel("Galactic longitude l (deg)")
    ax.set_ylabel("Galactic latitude b (deg)")
    ax.set_title(f"D6-2 left {best['name']} {best['t_yr'] / 1e3:.0f},000 years ago")
    ax.set_aspect("equal")
    ax.legend(loc="lower left", fontsize=7)

    ax = axes[1]
    le = res.get("look_elsewhere", {})
    per = le.get("p_per_star", {})
    names = sorted(per, key=lambda n: -per[n])
    ax.barh(range(len(names)), [100 * per[n] for n in names],
            color=[CLASS_STYLE[CLS[n]]["color"] for n in names], height=0.7)
    ax.set_yticks(range(len(names)), names)
    ax.invert_yaxis()
    ax.set_xlabel("chance of a match as close by coincidence (%)")
    ax.set_title("How surprising is the match?")
    pd = res["D6-2"]["chance"]["p_chance"]
    n_sh = le.get("n_shuffles", 1)
    txt = (f"D6-2: {100 * pd:.1f} ± {100 * np.sqrt(pd * (1 - pd) / n_sh):.1f}%  (its own track)\n"
           f"any D⁶ star: {100 * le.get('p_any_d6', np.nan):.1f}%\n"
           f"any of the 17 runaways: {100 * le.get('p_any_all', np.nan):.1f}%")
    ax.text(0.97, 0.05, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color=ps.INK)
    ps.caption(fig, f"Left: 150 posterior draws (they overlap: D6-2's distance is precise); the curve marks where an explosion at each past time "
                    f"would appear today, carried along by the disc's rotation. Remnants: Green (2025). "
                    + ("Background: Finkbeiner (2003) H-alpha map. " if halpha is not None else "")
                    + f"Right: remnants' longitudes shifted at random ({res['D6-2']['chance']['n_shuffles']:,} times), "
                    f"latitudes and sizes kept; a match counts if the track passes within {le.get('ratio', np.nan):.2f} "
                    "of a remnant's radius from its centre, as D6-2's does, within 300,000 years.")
    return ps.save(fig, _out("fig6_remnant"))


# --------------------------------------------------------------------------
# 7. Orbital speed vs donor mass
# --------------------------------------------------------------------------

def fig_donors(results):
    kin = results["kinematics"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 6.0), sharey=True,
                                  gridspec_kw=dict(width_ratios=[2.3, 1], wspace=0.05))
    m = np.linspace(0.2, 1.2, 200)
    for macc, col in ((0.85, ps.MUTED), (1.05, ps.BLUE), (1.2, ps.VIOLET)):
        mm = m[m <= macc]
        ax.plot(mm, binary.ejection_speed(mm, macc), color=col, lw=2,
                label=f"accretor {macc:.2f} M☉ (Bédard et al. radii at 10,000 K)")
        mc = np.linspace(0.05, macc, 200)
        ax.plot(mc, binary.ejection_speed(mc, macc, "eggleton"), color=col, lw=1, ls="--")
        ax.annotate(f"{macc:.2f}", (macc, binary.ejection_speed(macc, macc)), xytext=(4, -2),
                    textcoords="offset points", fontsize=7, color=col)
    ax.plot([], [], color=ps.MUTED, ls="--", lw=1, label="same, zero-temperature radii (fastest possible)")
    ax.set_xlabel("donor white dwarf mass (M☉)")
    ax.set_ylabel("launch speed of the surviving donor (km/s)")
    ax.set_xlim(0.0, 1.25)
    ax.set_ylim(400, 3000)
    ax.set_title("How heavy were the exploding binaries?")
    ax.legend(loc="upper left", fontsize=7)
    ax.text(0.98, 0.03, "a donor heavier than its accretor is impossible:\neach curve ends where the two masses are equal",
            transform=ax.transAxes, ha="right", fontsize=7, color=ps.MUTED)
    # the measured launch speeds, one column per star
    d6 = sorted([n for n in kin if kin[n]["cls"] == "D6"], key=lambda n: kin[n]["v_ej"]["median"])
    names = d6 + ["US 708"]
    for i, n in enumerate(names):
        v = kin[n]["v_ej"]
        col = CLASS_STYLE[kin[n]["cls"]]["color"]
        mk = CLASS_STYLE[kin[n]["cls"]]["marker"]
        ax2.plot([i, i], [v["p2"], v["p98"]], color=col, lw=1.2, alpha=0.45)
        ax2.plot([i, i], [v["median"] - v["minus"], v["median"] + v["plus"]], color=col, lw=4,
                 solid_capstyle="round")
        ax2.plot(i, v["median"], mk, color=col, ms=7, mec=ps.SURFACE)
    vmax105 = binary.max_speed(1.05)
    for a in (ax, ax2):
        a.axhline(vmax105, color=ps.BLUE, lw=0.8, ls=":")
    ax.text(0.22, vmax105 + 30, f"two 1.05 M☉ white dwarfs at most: {vmax105:,.0f} km/s", fontsize=7,
            color=ps.INK_SECONDARY)
    ax2.set_xticks(range(len(names)), names, rotation=60, ha="right", fontsize=7.5)
    ax2.set_xlim(-0.7, len(names) - 0.3)
    ax2.set_title("measured (68% / 95%)", fontsize=9)
    ax2.grid(axis="x", visible=False)
    ps.caption(fig, "Donor fills its Roche lobe (Eggleton 1983) when the accretor explodes; launch speed = orbital speed "
                    "plus the debris kick of El-Badry et al. (2023, eq. 7) in quadrature. Accretors above ~1.05 M☉ are "
                    "unlikely to have carbon-oxygen cores. US 708 (a helium-star donor) is shown for scale only.")
    return ps.save(fig, _out("fig7_donor_masses"))


# --------------------------------------------------------------------------
# 8. Survivors vs the Type Ia rate
# --------------------------------------------------------------------------

def fig_population(results):
    pop = results["population"]
    taus = np.array(pop["taus_myr"])
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    styles = {"D6 (all)": (ps.BLUE, "-", 2.4, "all D⁶ stars found"),
              "D6 bright, fast": (ps.BLUE, "--", 1.4, "bright, fast D⁶ stars"),
              "D6 faint": (ps.BLUE, ":", 1.6, "faint D⁶ stars (D6-2, J1332-3541)"),
              "LP40-365 (per Type Iax)": (ps.ORANGE, "-", 2.0, "LP 40-365 stars, as % of Type Iax")}
    for g, (col, ls, lw, lbl) in styles.items():
        c = pop["groups"].get(g)
        if not c:
            continue
        f = np.array(c["fraction"])
        if g in ("D6 (all)", "LP40-365 (per Type Iax)"):
            ax.fill_between(taus, 100 * np.array(c["fraction_lo"]), 100 * np.array(c["fraction_hi"]),
                            color=col, alpha=0.15, lw=0)
        ax.plot(taus, 100 * f, color=col, ls=ls, lw=lw, label=f"{lbl} ({c['n_found']})")
        if c.get("oldest_age_myr"):
            t0, f0 = c["oldest_age_myr"], c["fraction_at_oldest"]
            ax.plot(t0, 100 * f0, "o", color=col, ms=6, mec=ps.SURFACE, zorder=5)
            ax.annotate(f"≤ {100 * f0:.1f}%", (t0, 100 * f0), xytext=(6, 4), textcoords="offset points",
                        fontsize=7.5, color=ps.INK)
    d62 = pop["per_star"].get("D6-2")
    if d62 and d62.get("fraction_of_parent_at_age"):
        t0, f0 = d62["age_myr"], d62["fraction_of_parent_at_age"]
        ax.plot(t0, 100 * f0, "D", color=ps.BLUE, ms=6, mec=ps.SURFACE, zorder=5)
        ax.annotate(f"stars like D6-2, if they fade as\nfast as its remnant's age allows: {100 * f0:.0f}%",
                    (t0, 100 * f0), xytext=(0.0115, 9.0), textcoords="data", fontsize=7.5, color=ps.INK,
                    arrowprops=dict(arrowstyle="-", color=ps.MUTED, lw=0.6))
    ax.axhline(100, color=ps.INK_SECONDARY, lw=1)
    ax.text(taus[-1], 118, "every supernova of the parent type", fontsize=7.5, color=ps.INK_SECONDARY, ha="right")
    ax.axhline(100 * REFERENCE["shen25_fraction"], color=ps.MAGENTA, lw=1, ls="--")
    ax.text(taus[0], 100 * REFERENCE["shen25_fraction"] * 0.72, "Shen (2025): ~2% of Type Ia",
            fontsize=7.5, color=ps.INK_SECONDARY, ha="left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.05, 3e4)
    ax.set_xlabel("how long a survivor stays this bright, τ (Myr)")
    ax.set_ylabel("survivors needed, as % of the parent supernova rate")
    ax.set_title("What fraction of thermonuclear supernovae leave a runaway?")
    ax.legend(loc="upper right", fontsize=7)
    ps.caption(fig, "Each star found by the search stands for look-alikes (its own launch speed, absolute magnitude "
                    "and colour) launched from the disc at the Milky Way's Type Ia rate (0.54 per century; Li et al. 2011) "
                    "and observed with real Gaia DR3 parallax errors through the search's cuts. Type Iax rate: 31% of "
                    "Type Ia (Foley et al. 2013). Bands: Poisson 68%. Dots: each group's oldest flight time, a lower "
                    "limit on τ, so the fraction there is an upper bound.")
    return ps.save(fig, _out("fig8_population"))


# --------------------------------------------------------------------------
# 9. The runaways on the sky
# --------------------------------------------------------------------------

def _to_aitoff(l_deg, b_deg):
    """Galactic (l, b) in degrees -> matplotlib aitoff radians, l increasing
    to the left with the Galactic Centre in the middle."""
    lw = ((np.asarray(l_deg) + 180) % 360) - 180
    return -np.radians(lw), np.radians(b_deg)


def fig_sky(results, orbit_art, rem_art):
    import astropy.units as u
    from astropy.coordinates import Galactic, SkyCoord
    from kinematics import galcen_frame
    green = rem_art["green"]
    fig = plt.figure(figsize=(12, 6.4))
    ax = fig.add_subplot(111, projection="aitoff")
    x, y = _to_aitoff(green["l"], green["b"])
    ax.scatter(x, y, s=np.clip(6 + 40 * np.asarray(green["radius_deg"]), 6, 160), facecolor="none",
               edgecolor=ps.MUTED, lw=0.6, label="supernova remnants (Green 2025)", zorder=2)
    for name, a in orbit_art.items():
        st = CLASS_STYLE[CLS[name]]
        med = np.median(a["past"], axis=1)[:, :51]          # today -> 1 Myr ago
        g = SkyCoord(x=med[0] * u.kpc, y=med[1] * u.kpc, z=med[2] * u.kpc,
                     frame=galcen_frame()).transform_to(Galactic())
        lx, by = _to_aitoff(g.l.deg, g.b.deg)
        # break the line where it wraps around the map edge
        jump = np.flatnonzero(np.abs(np.diff(lx)) > np.pi) + 1
        for seg_x, seg_y in zip(np.split(lx, jump), np.split(by, jump)):
            ax.plot(seg_x, seg_y, color=st["color"], lw=1.2, alpha=0.8, zorder=3)
        ax.plot(lx[0], by[0], st["marker"], color=st["color"], ms=7, mec=ps.SURFACE, mew=0.8, zorder=4)
        ax.annotate(name, (lx[0], by[0]), xytext=(5, 4), textcoords="offset points", fontsize=6.5,
                    color=ps.INK_SECONDARY, zorder=5)
    # x = -l (wrapped to +-180), so the tick at x = -150 deg is l = 150 deg
    ticks = np.array([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150])
    ax.set_xticks(np.radians(ticks), [f"{(-t) % 360:.0f}°" for t in ticks], fontsize=7, color=ps.MUTED)
    ax.set_yticks(np.radians([-60, -30, 0, 30, 60]), ["−60°", "−30°", "0°", "30°", "60°"], fontsize=7,
                  color=ps.MUTED)
    ax.grid(color=ps.GRID, lw=0.6)
    for cls in ("D6", "LP40-365", "He-star donor", "other"):
        st = CLASS_STYLE[cls]
        ax.plot([], [], st["marker"], color=st["color"], ls="-", lw=1.2, label=st["label"])
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=7)
    ax.set_title("Where they are, and where they were 1 Myr ago (lines)", pad=26)
    ps.caption(fig, "Galactic coordinates, longitude increasing to the left. Lines: the median past path over 1 Myr, "
                    "as it would have appeared from the Sun's present position. Remnant circles scale with size.")
    return ps.save(fig, _out("fig9_sky"))
