"""Every figure in the project, from the outputs of ``analysis.run_all``.

    01_sky          Where FRBs are on the sky, next to pulsars
    02_dispersion   What dispersion does to a burst, and a test of the 1/f^2 law
    03_dm_excess    More dispersion than the whole Milky Way can supply
    04_repeaters    Repeating sources emit differently shaped bursts
    05_macquart     DM vs redshift: weighing the diffuse baryons

Colour use follows one rule per job. Categorical identity uses the palette's
first three slots only (blue, orange, aqua), the set that stays distinguishable
for colour-blind readers in scatter plots; anything beyond that is grey context.
Intensity in dynamic spectra uses a single-hue blue ramp that fades into the
background at zero.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import frb_physics as fp
import plotstyle as ps
import waterfalls as wf
from config import (
    CHIME_DEC_LIMIT_DEG,
    DM_HALO,
    DM_HALO_MAX,
    F_DIFFUSE_FIDUCIAL,
    FIGURE_DIR,
    MACQUART_2020_LOG10_F,
    WRONG_DM_BURST,
    WRONG_DM_FRACTION,
)

ps.use_style()

C_ONEOFF = ps.BLUE
C_REPEATER = ps.ORANGE
C_THIRD = ps.AQUA
C_CONTEXT = ps.MUTED
C_BLIND = "#eeede8"
C_FLAGGED = "#e9e8e3"

SURVEY_STYLE = {
    "ASKAP (CRAFT)": (ps.BLUE, "o"),
    "DSA-110": (ps.ORANGE, "s"),
    "CHIME Outriggers": (ps.AQUA, "D"),
    "MeerKAT": (ps.MUTED, "^"),
    "Other": (ps.MUTED, "v"),
}

#: Sequential ramp for intensity: the palette's blue steps, starting at the
#: chart surface so that zero signal recedes into the background.
CMAP = LinearSegmentedColormap.from_list(
    "frb_blue",
    [ps.SURFACE, "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
)
CMAP.set_bad(C_FLAGGED)


def _out(name):
    return FIGURE_DIR / name


def _fmt_frb(tns: str) -> str:
    """'FRB20211005A' -> 'FRB 20211005A' (the conventional printed form)."""
    return tns.replace("FRB", "FRB ", 1) if tns.startswith("FRB") and " " not in tns else tns


# ---------------------------------------------------------------------------
# 01 -- the sky
# ---------------------------------------------------------------------------


def _mollweide_x(l_deg):
    """Galactic longitude -> Mollweide x (radians), l increasing to the left."""
    lw = (np.asarray(l_deg, float) + 180.0) % 360.0 - 180.0
    return -np.deg2rad(lw)


def figure_sky(results: dict, art: dict) -> list:
    pop, psr = art["pop"], art["psr"]
    frb = pop[pop["gl"].notna()]
    one, rep = frb[~frb["is_repeater"]], frb[frb["is_repeater"]]

    fig = plt.figure(figsize=(10, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.1, 1], hspace=0.32, top=0.93, bottom=0.1)
    ax = fig.add_subplot(gs[0], projection="mollweide")

    # The part of the sky CHIME cannot see: Dec < -11 deg.
    xs = np.linspace(-180, 180, 721)
    ys = np.linspace(-89.9, 89.9, 360)
    X, Y = np.meshgrid(xs, ys)
    eq = SkyCoord(l=((-X) % 360) * u.deg, b=Y * u.deg, frame="galactic").icrs
    blind = (eq.dec.deg < CHIME_DEC_LIMIT_DEG).astype(float)
    shade = ax.contourf(np.deg2rad(X), np.deg2rad(Y), blind, levels=[0.5, 1.5], colors=[C_BLIND], zorder=0)
    shade.set_rasterized(True)  # a jagged 721 x 360 boundary is megabytes as vector paths

    ax.scatter(_mollweide_x(psr["gl"]), np.deg2rad(psr["gb"]), s=1.2, color=C_CONTEXT, alpha=0.55,
               linewidths=0, zorder=2, rasterized=True)
    ax.scatter(_mollweide_x(one["gl"]), np.deg2rad(one["gb"]), s=3.0, color=C_ONEOFF, alpha=0.75,
               linewidths=0, zorder=3, rasterized=True)
    ax.scatter(_mollweide_x(rep["gl"]), np.deg2rad(rep["gb"]), s=26, color=C_REPEATER,
               edgecolors=ps.SURFACE, linewidths=1.0, zorder=4)

    # Longitude labels sit just below the plane, on a surface-coloured backing so
    # the band of pulsars behind them does not swallow them.
    ticks = np.arange(-150, 151, 30)
    ax.set_xticks(np.deg2rad(ticks))
    ax.set_xticklabels([])
    for lon in (60, 120, 240, 300):
        ax.text(_mollweide_x(lon), np.deg2rad(-7), f"l = {lon}°", ha="center", va="top", fontsize=7,
                color=ps.MUTED, bbox={"facecolor": ps.SURFACE, "edgecolor": "none", "alpha": 0.8, "pad": 0.6},
                zorder=6)
    ax.set_yticks(np.deg2rad(np.arange(-60, 61, 30)))
    ax.set_yticklabels([f"{t:+.0f}°" if t else "0°" for t in np.arange(-60, 61, 30)], fontsize=7, color=ps.MUTED)
    ax.grid(True, color=ps.GRID, linewidth=0.6)

    # Landmarks.
    ncp_x, ncp_y = _mollweide_x(122.93), np.deg2rad(27.13)
    ax.plot([ncp_x], [ncp_y], marker="+", color=ps.INK, markersize=9, mew=1.4, zorder=5)
    ps.annotate(ax, "north celestial pole:\nCHIME watches this\nsky the longest", (ncp_x, ncp_y),
                (0.0, 0.93), color=ps.INK_SECONDARY, textcoords="axes fraction")
    ax.text(np.deg2rad(40), np.deg2rad(-52), f"not visible to CHIME\n(declination below {CHIME_DEC_LIMIT_DEG:.0f}°)",
            ha="center", va="center", fontsize=8, color=ps.INK_SECONDARY)
    ax.text(0, np.deg2rad(4), "Galactic centre, l = 0°", ha="center", va="bottom", fontsize=7.5,
            color=ps.INK_SECONDARY, bbox={"facecolor": ps.SURFACE, "edgecolor": "none", "alpha": 0.8, "pad": 0.6},
            zorder=6)
    ax.set_title("Fast radio bursts come from every direction; pulsars trace the Milky Way", pad=18)

    handles = [
        Line2D([], [], ls="", marker="o", ms=4, color=C_ONEOFF, label=f"one-off FRB sources ({len(one):,})"),
        Line2D([], [], ls="", marker="o", ms=6.5, color=C_REPEATER, mec=ps.SURFACE,
               label=f"repeating FRB sources ({len(rep)})"),
        Line2D([], [], ls="", marker="o", ms=2.5, color=C_CONTEXT, label=f"pulsars, all telescopes ({len(psr):,})"),
        Patch(facecolor=C_BLIND, edgecolor="none", label="outside CHIME's view"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.17), ncol=4, fontsize=8,
              handletextpad=0.4, columnspacing=1.4)

    # Panel B: density relative to an isotropic sky seen through the same exposure.
    axb = fig.add_subplot(gs[1])
    sky = art["sky"]
    for key, colour, label in (("pulsar_chime_sky", C_CONTEXT, "pulsars in CHIME's part of the sky"),
                               ("frb", C_ONEOFF, "FRB sources")):
        h = sky[key]
        obs, exp, edges = h["hist_observed"], h["hist_isotropic"], h["edges"]
        centre = 0.5 * (edges[1:] + edges[:-1])
        ok = obs > 0
        ratio = obs[ok] / exp[ok]
        err = np.sqrt(obs[ok]) / exp[ok]
        axb.errorbar(centre[ok], ratio, yerr=err, color=colour, marker="o", ms=4.5, lw=1.6, capsize=0,
                     mec=ps.SURFACE, mew=1.0, label=label, zorder=3)
    axb.axhline(1.0, color=ps.INK_SECONDARY, lw=0.9, zorder=1)
    axb.set_yscale("log")
    axb.set_ylim(0.03, 12)
    axb.set_xlim(0, 90)
    axb.set_xticks(np.arange(0, 91, 10))
    axb.set_yticks([0.1, 0.3, 1, 3, 10])
    axb.set_yticklabels(["0.1", "0.3", "1", "3", "10"])
    axb.set_xlabel("distance from the Galactic plane, |b| (degrees)")
    axb.set_ylabel("sources per bin ÷ isotropic expectation")
    axb.text(62, 0.9, "isotropic sky, same exposure", ha="center", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    s = results["sky"]
    axb.set_title("How concentrated toward the Galactic plane?", fontsize=10)
    axb.legend(loc="lower left", fontsize=8)
    ps.annotate(axb, f"FRBs: {s['frb']['frac_abs_b_lt_10_observed']:.0%} within 10° of the plane vs "
                     f"{s['frb']['frac_abs_b_lt_10_isotropic']:.0%} if isotropic.\n"
                     "A mild deficit: the Milky Way's own radio glow and\nturbulent gas make bursts harder to detect there",
                (7.5, 0.84), (38, 3.2))
    ps.annotate(axb, f"pulsars: {s['pulsar_chime_sky']['frac_abs_b_lt_10_observed']:.0%} within 10° vs "
                     f"{s['pulsar_chime_sky']['frac_abs_b_lt_10_isotropic']:.0%}",
                (2.5, 5.0), (14, 6.5))

    ps.caption(fig, (
        "Galactic coordinates (Mollweide), longitude increasing to the left. FRBs: CHIME/FRB Catalog 2, one point per source "
        "passing CHIME's population cuts. Pulsars: ATNF Pulsar Catalogue, found by many telescopes in both hemispheres. "
        "Lower panel: isotropic expectation from 'RA scrambling' -- each source keeps its declination and is spread uniformly in "
        "right ascension, which reproduces the exposure of a transit telescope without modelling it. Bursts pile up near the "
        "celestial pole because CHIME sees circumpolar sky twice a day, not because FRBs cluster there."
    ))
    return ps.save(fig, _out("01_sky"))


# ---------------------------------------------------------------------------
# 02 -- dispersion
# ---------------------------------------------------------------------------


def _waterfall(ax, ds: wf.DynamicSpectrum, data=None, *, tlim=(-40, 40), vmax=None, label_y=True):
    d = ds.data if data is None else data
    t_ms = ds.t * 1e3
    sel = (t_ms >= tlim[0] - 1) & (t_ms <= tlim[1] + 1)
    if vmax is None:
        vmax = np.nanpercentile(ds.data[:, np.abs(ds.t) < 3 * ds.dt], 97)
    mesh = ax.pcolormesh(t_ms[sel], ds.freqs, d[:, sel], cmap=CMAP, vmin=0, vmax=vmax, shading="nearest",
                         rasterized=True)
    ax.set_xlim(*tlim)
    ax.set_ylim(400, 800)
    ax.grid(False)
    ax.set_xlabel("time (ms)")
    if label_y:
        ax.set_ylabel("frequency (MHz)")
    return mesh, vmax


def _profile(ax, ds: wf.DynamicSpectrum, *, tlim, color=C_ONEOFF):
    t_ms = ds.t * 1e3
    ax.plot(t_ms, ds.profile(), color=color, lw=1.2)
    ax.set_xlim(*tlim)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.tick_params(labelbottom=False)
    ax.grid(False)


def figure_dispersion(results: dict, art: dict) -> list:
    disp = art["dispersion"]
    names = list(disp.keys())

    fig = plt.figure(figsize=(11, 12.2))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.15, 1.0, 0.95], hspace=0.45, wspace=0.28, top=0.93, bottom=0.09)

    # A: the sweep on the true time axis.
    ax = fig.add_subplot(gs[0, :])
    colours = [C_ONEOFF, C_REPEATER]
    ff = np.linspace(400.2, 800.2, 400)
    for tns, colour in zip(names, colours):
        r = disp[tns]
        t_top = r["t0_alpha_fixed"] + r["k_dm"] * r["dm_alpha_fixed"] * 800.2 ** -2
        ax.plot(r["t0_alpha_fixed"] + r["k_dm"] * r["dm_alpha_fixed"] * ff ** -2 - t_top, ff,
                color=colour, lw=1.0, alpha=0.6, zorder=2)
        ax.scatter(r["t_arrival_s"] - t_top, r["freq_mhz"], s=22, color=colour, edgecolors=ps.SURFACE,
                   linewidths=1.0, zorder=3, label=f"{_fmt_frb(tns)}  (DM {r['dm_chime']:.0f} pc cm⁻³)")
        ax.text(r["sweep_s"] + 0.15, 402, f"arrives {r['sweep_s']:.1f} s\nafter 800 MHz", fontsize=8,
                color=ps.INK_SECONDARY, va="bottom", ha="left")
    ax.set_xlim(-0.3, max(disp[n]["sweep_s"] for n in names) + 1.6)
    ax.set_ylim(395, 805)
    ax.set_xlabel("seconds after the burst reached the top of the band (800 MHz)")
    ax.set_ylabel("frequency (MHz)")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("A: one instant of emission, spread over seconds -- each dot is when the burst arrived in a 6 MHz slice")

    # B-D: the same burst dedispersed with the right DM and with DM off by 1%.
    ds = wf.load(WRONG_DM_BURST)
    panels = [
        (ds.dm, f"B: dedispersed, DM = {ds.dm:.2f}"),
        (ds.dm * (1 - WRONG_DM_FRACTION), f"C: DM {WRONG_DM_FRACTION:.0%} too low"),
        (ds.dm * (1 + WRONG_DM_FRACTION), f"D: DM {WRONG_DM_FRACTION:.0%} too high"),
    ]
    vmax = None
    mesh = None
    for j, (dm_new, title) in enumerate(panels):
        axw = fig.add_subplot(gs[1, j])
        data = ds.data if dm_new == ds.dm else wf.redisperse(ds, dm_new)
        mesh, vmax = _waterfall(axw, ds, data, tlim=(-60, 60), vmax=vmax, label_y=(j == 0))
        axw.set_title(title, fontsize=9.5)
        if j:
            resid = fp.dispersion_delay(dm_new - ds.dm, 400.2, 800.2, k=ds.k_dm) * 1e3
            axw.text(0.03, 0.04, f"residual sweep {abs(resid):.0f} ms", transform=axw.transAxes,
                     fontsize=7.5, color=ps.INK_SECONDARY)
    cb = fig.colorbar(mesh, ax=fig.axes[-3:], fraction=0.025, pad=0.015)
    cb.set_label("signal-to-noise per pixel", fontsize=8)
    cb.outline.set_visible(False)

    # E: arrival-time residuals from the cold-plasma law.
    axe = fig.add_subplot(gs[2, :2])
    for k, (tns, colour) in enumerate(zip(names, colours)):
        r = disp[tns]
        axe.errorbar(r["freq_mhz"], r["residuals_s"] * 1e3, yerr=r["t_err_s"] * 1e3, ls="", marker="o", ms=4,
                     color=colour, mec=ps.SURFACE, mew=0.8, elinewidth=0.8, capsize=0,
                     label=(f"{_fmt_frb(tns)}: free index α = {r['alpha']:.5f} ± {r['alpha_err']:.5f}; "
                            f"rms {r['residual_rms_ms']:.2f} ms"))
    axe.axhline(0, color=ps.INK_SECONDARY, lw=0.9)
    axe.set_xlim(400, 800)
    lim = max(0.45, 1.15 * max(np.max(np.abs(disp[n]["residuals_s"])) * 1e3 for n in names))
    axe.set_ylim(-lim, lim)
    axe.set_xlabel("frequency (MHz)")
    axe.set_ylabel("arrival time − ν⁻² law (ms)")
    axe.set_title("E: what is left after removing a ν⁻² delay of seconds", fontsize=9.5)
    axe.legend(loc="upper right", fontsize=7.5)

    # F: scattering, the other fingerprint of the plasma.
    axf = fig.add_subplot(gs[2, 2])
    sc = wf.load("FRB20230905D")
    _waterfall(axf, sc, tlim=(-20, 150), label_y=False)
    axf.set_title("F: scattering, a second fingerprint", fontsize=9.5)
    tau400 = 84.3
    axf.text(0.97, 0.96, f"tail ∝ ν⁻⁴:\n{tau400:.0f} ms at 400 MHz,\n{tau400 * (400 / 800) ** 4:.0f} ms at 800 MHz",
             transform=axf.transAxes, ha="right", va="top", fontsize=7.5, color=ps.INK_SECONDARY)

    fig.suptitle("Dispersion: the delay that counts the electrons between us and the burst", x=0.06, ha="left",
                 y=0.975)
    ps.caption(fig, (
        "CHIME/FRB Catalog 2 dynamic spectra, one-off bursts FRB 20211005A and FRB 20200711F. A: arrival time in up to 64 sub-bands, "
        "measured by fitting each sub-band's pulse and converting to the recorded timestamps (lines: t ∝ DM/ν² with the fitted "
        "DM). B-D: 16384 channels averaged to 1.56 MHz; grey rows are channels removed for radio interference. Re-dedispersing "
        "with a DM off by 1% tilts the burst by tens of ms. E: residuals from the fixed ν⁻² fit; the free-index values quote "
        "errors scaled by sqrt(chi²/dof). F: FRB 20230905D, the highest-S/N strongly scattered one-off in the catalogue; the public "
        "window clips the lowest-frequency tail."
    ))
    return ps.save(fig, _out("02_dispersion"))


# ---------------------------------------------------------------------------
# 03 -- more DM than the Milky Way holds
# ---------------------------------------------------------------------------


def figure_dm_excess(results: dict, art: dict) -> list:
    frb, psr = art["dm"]["frb"], art["dm"]["psr"]
    r = results["dm"]
    fig, (ax, axh) = plt.subplots(1, 2, figsize=(11.5, 5.6),
                                  gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.28, "top": 0.86})

    gal = psr[~psr["magellanic"]]
    mc = psr[psr["magellanic"]]
    ax.scatter(gal["dm_mw_ne2001"], gal["dm"], s=4, color=C_CONTEXT, alpha=0.6, linewidths=0, rasterized=True,
               label=f"pulsars in the Milky Way ({len(gal):,})", zorder=2)
    ax.scatter(mc["dm_mw_ne2001"], mc["dm"], s=22, color=C_THIRD, edgecolors=ps.SURFACE, linewidths=0.8,
               label=f"pulsars in the Magellanic Clouds ({len(mc)})", zorder=4)
    ax.scatter(frb["dm_mw_ne2001"], frb["dm_fitb"], s=5, color=C_ONEOFF, alpha=0.55, linewidths=0, rasterized=True,
               label=f"FRB sources, CHIME ({len(frb):,})", zorder=3)
    x = np.geomspace(10, 3000, 200)
    ax.plot(x, x, color=ps.INK_SECONDARY, lw=1.0, zorder=5)
    ax.plot(x, x + DM_HALO_MAX, color=ps.INK_SECONDARY, lw=0.8, alpha=0.6, zorder=5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(12, 2500)
    ax.set_ylim(1.5, 6000)
    ax.set_xlabel("NE2001: DM through the entire Milky Way disk in that direction (pc cm⁻³)")
    ax.set_ylabel("observed DM (pc cm⁻³)")
    handles, _ = ax.get_legend_handles_labels()
    handles += [
        Line2D([], [], color=ps.INK_SECONDARY, lw=1.0, label="observed = all the disk holds"),
        Line2D([], [], color=ps.INK_SECONDARY, lw=0.8, alpha=0.6, label=f"... plus {DM_HALO_MAX:.0f} for the halo (generous)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7.5, markerscale=1.6)
    ax.set_title("A: observed DM against everything the Milky Way can supply")
    ps.annotate(ax, f"median FRB: {r['frb']['median_ratio']:.0f}× what the whole disk holds", (60, 700), (16, 3300))

    # B: the extragalactic excess, and what redshift it implies.
    dm_eg = frb["dm_eg"].to_numpy()
    bins = np.arange(-100, 3001, 50)
    axh.hist(dm_eg, bins=bins, color=C_ONEOFF, alpha=0.85, rwidth=0.88, zorder=3)
    axh.axvline(0, color=ps.INK_SECONDARY, lw=0.9)
    axh.set_xlim(-150, 2600)
    axh.set_xlabel(f"DM beyond the Milky Way (disk from NE2001, halo {DM_HALO:.0f}) (pc cm⁻³)")
    axh.set_ylabel("FRB sources per 50 pc cm⁻³")
    axh.set_title("B: the excess, and the distance it implies")
    top = axh.secondary_xaxis(
        "top",
        functions=(lambda d: np.nan_to_num(fp.z_from_mean_dm(d), nan=6.0),
                   lambda zz: fp.mean_dm_cosmic(np.clip(zz, 0, 5.9))),
    )
    top.set_xticks([0.25, 0.5, 1.0, 1.5, 2.0])
    top.set_xticklabels(["0.25", "0.5", "1", "1.5", "2"])
    top.set_xlabel("redshift if all of the excess were intergalactic (an upper limit)", fontsize=8, color=ps.INK_SECONDARY)
    top.tick_params(colors=ps.MUTED, labelsize=7.5)
    med = r["frb"]["median_dm_eg"]
    ps.annotate(axh, f"median {med:.0f} pc cm⁻³\n→ z ≲ {r['frb']['median_z_max']:.2f}, light that left\n"
                     f"its source ~{_lookback_gyr(r['frb']['median_z_max']):.0f} billion years ago",
                (med, axh.get_ylim()[1] * 0.55), (med + 450, axh.get_ylim()[1] * 0.7))

    ps.caption(fig, (
        "FRBs: CHIME/FRB Catalog 2, first burst of each source passing the population cuts; their Milky Way column is CHIME's "
        "own NE2001 value for the sightline. Pulsars: ATNF catalogue, NE2001 through the whole disk toward each (HEALPix map "
        f"from the FRBs/FRB repository). {r['pulsar']['milky_way_all']['frac_above_ne2001']:.0%} of Milky Way pulsars sit above the line, a measure of NE2001's errors on "
        "individual sightlines; Magellanic Cloud pulsars sit above it because their own galaxy adds DM, just as an FRB's host "
        "does. Top axis of B: mean intergalactic DM for Planck 2018 with 84% of baryons diffuse; a host galaxy's contribution "
        "would place each burst nearer."
    ))
    return ps.save(fig, _out("03_dm_excess"))


def _lookback_gyr(z: float) -> float:
    from astropy.cosmology import Planck18

    return float(Planck18.lookback_time(z).to(u.Gyr).value)


# ---------------------------------------------------------------------------
# 04 -- repeaters vs one-offs
# ---------------------------------------------------------------------------


def _ecdf(ax, values, *, color, label, lw=2.0, alpha=1.0):
    v = np.sort(np.asarray(values, float)[np.isfinite(values)])
    y = np.arange(1, v.size + 1) / v.size
    ax.step(v, y, where="post", color=color, lw=lw, alpha=alpha, label=label)


def figure_repeaters(results: dict, art: dict) -> list:
    samples = art["repeaters"]
    r = results["repeaters"]
    gallery = [
        ("FRB20211005A", "one-off", C_ONEOFF),
        ("FRB20200711F", "one-off", C_ONEOFF),
        ("FRB20181222A", "repeater", C_REPEATER),
        ("FRB20211114A", "repeater", C_REPEATER),
    ]
    fig = plt.figure(figsize=(12, 9.6))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.05, 1], hspace=0.36, top=0.9, bottom=0.1)
    top = outer[0].subgridspec(2, 4, height_ratios=[0.28, 1], hspace=0.05, wspace=0.12)
    tlim = (-30, 30)
    for j, (tns, kind, colour) in enumerate(gallery):
        ds = wf.load(tns)
        axp = fig.add_subplot(top[0, j])
        axw = fig.add_subplot(top[1, j], sharex=axp)
        _profile(axp, ds, tlim=tlim, color=colour)
        _waterfall(axw, ds, tlim=tlim, label_y=(j == 0))
        if j:
            axw.tick_params(labelleft=False)
        source = _fmt_frb(ds.repeater_name) if ds.repeater_name else _fmt_frb(tns)
        head = "one-off" if kind == "one-off" else "repeater"
        sub = f"burst {_fmt_frb(tns)}" if ds.repeater_name else "never seen to repeat"
        axp.set_title(f"{head}: {source}\n{sub}", fontsize=8.5, loc="left")

    bottom = outer[1].subgridspec(1, 2, wspace=0.25)
    one, rep, rep_all = samples["one_off"], samples["repeater_first"], samples["repeater_all"]
    q = r["quantities"]

    ax1 = fig.add_subplot(bottom[0])
    _ecdf(ax1, one["width_ms"], color=C_ONEOFF, label=f"one-off sources ({len(one):,})")
    _ecdf(ax1, rep["width_ms"], color=C_REPEATER, label=f"repeaters, first burst ({len(rep)})")
    _ecdf(ax1, rep_all["width_ms"], color=C_REPEATER, lw=1.0, alpha=0.45, label=f"repeaters, all bursts ({len(rep_all)})")
    ax1.set_xscale("log")
    ax1.set_xlim(0.03, 30)
    ax1.set_ylim(0, 1.02)
    ax1.set_xlabel("intrinsic burst width, Gaussian σ (ms)")
    ax1.set_ylabel("fraction of bursts narrower than this")
    w = q["width_ms"]
    ax1.set_title(f"Repeater bursts last longer: median {w['repeater_first']['median']:.2f} vs "
                  f"{w['one_off']['median']:.2f} ms", fontsize=10)
    ax1.legend(loc="upper left", fontsize=8)
    ax1.text(0.98, 0.04, f"KS test p = {w['test_one_off_vs_repeater_first']['ks_pvalue']:.0e}",
             transform=ax1.transAxes, ha="right", fontsize=7.5, color=ps.INK_SECONDARY)

    ax2 = fig.add_subplot(bottom[1])
    _ecdf(ax2, one["bandwidth_mhz"], color=C_ONEOFF, label="one-off sources")
    _ecdf(ax2, rep["bandwidth_mhz"], color=C_REPEATER, label="repeaters, first burst")
    _ecdf(ax2, rep_all["bandwidth_mhz"], color=C_REPEATER, lw=1.0, alpha=0.45, label="repeaters, all bursts")
    ax2.set_xlim(0, 410)
    ax2.set_ylim(0, 1.02)
    ax2.set_xlabel("emission bandwidth, full width at tenth maximum (MHz)")
    ax2.set_ylabel("fraction of bursts narrower than this")
    b = q["bandwidth_mhz"]
    bf = r["band_filling_fraction"]
    ax2.set_title(f"...and cover less of the band: median {b['repeater_first']['median']:.0f} vs "
                  f"{b['one_off']['median']:.0f} MHz", fontsize=10)
    ax2.legend(loc="upper left", fontsize=8)
    ps.annotate(ax2, f"{bf['one_off']:.0%} of one-offs and {bf['repeater_first']:.0%} of repeaters\n"
                     "fill CHIME's whole 400 MHz band:\ntheir true bandwidth is larger",
                (399, 0.8), (205, 0.14))
    ax2.text(0.98, 0.04, f"KS test p = {b['test_one_off_vs_repeater_first']['ks_pvalue']:.0e}",
             transform=ax2.transAxes, ha="right", fontsize=7.5, color=ps.INK_SECONDARY)

    fig.suptitle("Repeating sources emit different bursts: longer, narrower in frequency, drifting downward",
                 x=0.06, ha="left", y=0.965)
    ps.caption(fig, (
        "Top: CHIME dynamic spectra (1.56 MHz x 0.98 ms pixels, shaded by signal-to-noise; grey rows: interference removed) with "
        "the frequency-averaged pulse above. The repeater bursts show several sub-bursts whose frequency drifts down with time. "
        "Bottom: CHIME/FRB Catalog 2 fitburst results, excluding flagged events and unstable fits. Each repeater contributes "
        "one burst (its first usable one) so that FRB 20220912A's 373 bursts cannot dominate; the faint lines use every repeat "
        "burst. 'One-off' means not yet seen to repeat -- some will."
    ))
    return ps.save(fig, _out("04_repeaters"))


# ---------------------------------------------------------------------------
# 05 -- the Macquart relation
# ---------------------------------------------------------------------------


#: Display names for the Macquart fit variants in panel C.
FIT_LABELS = {
    "baseline": "this project (baseline)",
    "macquart2020_setup": "F ≤ 0.32, no MW scatter (Macquart+20 set-up)",
    "no_mw_scatter": "no Milky Way scatter term",
    "halo_30": "Milky Way halo 30, not 50",
    "halo_80": "Milky Way halo 80, not 50",
    "no_outlier": "without FRB 20190520B",
    "z_ge_0.1": "only z ≥ 0.1",
}


def _hpd_levels(p2d: np.ndarray, fracs=(0.68, 0.95)) -> list:
    """Density thresholds enclosing the given posterior mass (for contours)."""
    flat = np.sort(p2d.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    return [flat[np.searchsorted(cum, f)] for f in fracs]


def figure_macquart(results: dict, art: dict) -> list:
    m = art["macquart"]
    sample, fits, band = m["sample"], m["fits"], m["band"]
    base = fits["baseline"]
    mr = results["macquart"]
    bj = mr["fits"]["baseline"]

    fig = plt.figure(figsize=(12, 7.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], height_ratios=[1, 1], hspace=0.42, wspace=0.26,
                          top=0.9, bottom=0.12)
    ax = fig.add_subplot(gs[:, 0])

    z = band["z"]
    ax.fill_between(z, band["q05"], band["q95"], color=ps.BLUE, alpha=0.10, lw=0, zorder=1,
                    label="model: 90% of bursts")
    ax.fill_between(z, band["q16"], band["q84"], color=ps.BLUE, alpha=0.18, lw=0, zorder=1,
                    label="model: 68% of bursts")
    ax.plot(z, band["q50"], color=ps.BLUE, lw=2.0, zorder=2, label="model: median burst")
    zz = np.linspace(0.0, 1.2, 200)
    ax.plot(zz, fp.mean_dm_cosmic(zz, f_d=F_DIFFUSE_FIDUCIAL), color=ps.INK_SECONDARY, lw=1.0, zorder=2,
            label=f"mean intergalactic DM if {F_DIFFUSE_FIDUCIAL:.0%} of baryons are diffuse")
    for survey, (colour, marker) in SURVEY_STYLE.items():
        s = sample[sample["survey"] == survey]
        if s.empty:
            continue
        ax.errorbar(s["z"], s["dm_eg"], yerr=base.sigma_mw, ls="", marker=marker, ms=5.5, color=colour,
                    mec=ps.SURFACE, mew=0.8, elinewidth=0.7, ecolor=colour, alpha=0.95, capsize=0, zorder=4,
                    label=f"{survey} ({len(s)})")
    for name, dx, dy in (("FRB20190520B", 0.05, 60), ("FRB20220610A", -0.02, 280)):
        row = sample[sample["name"] == name]
        if not row.empty:
            ps.annotate(ax, _fmt_frb(name) + ("\nhost adds ~900 pc cm⁻³" if name == "FRB20190520B" else "\nfarthest, z = 1.02"),
                        (float(row["z"].iloc[0]), float(row["dm_eg"].iloc[0])),
                        (float(row["z"].iloc[0]) + dx, float(row["dm_eg"].iloc[0]) + dy))
    ax.set_xlim(0, 1.12)
    ax.set_ylim(-60, 2100)
    ax.set_xlabel("host-galaxy redshift z")
    ax.set_ylabel(f"extragalactic DM = observed − NE2001 − {DM_HALO:.0f} (pc cm⁻³)")
    ax.set_title(f"A: the Macquart relation, {bj['n_frb']} localized FRBs")
    ax.legend(loc="upper left", fontsize=7.5, ncol=1)

    # B: the f_d - F degeneracy.
    axb = fig.add_subplot(gs[0, 1])
    post = np.exp(base.log_post).sum(axis=(2, 3))  # (f_d, log10 F)
    levels = sorted(_hpd_levels(post))
    F = 10 ** base.axes["log10_F"]
    axb.contourf(F, base.axes["f_d"], post, levels=[levels[0], levels[1], post.max() * 1.01],
                 colors=[ps.BLUE, ps.BLUE], alpha=0.25)
    axb.contourf(F, base.axes["f_d"], post, levels=[levels[1], post.max() * 1.01], colors=[ps.BLUE], alpha=0.35)
    axb.contour(F, base.axes["f_d"], post, levels=levels, colors=[ps.BLUE], linewidths=1.0)
    lo, hi = 10 ** np.array(MACQUART_2020_LOG10_F[:2])
    axb.axvspan(lo, hi, color=ps.GRID, alpha=0.6, lw=0, zorder=0)
    axb.text(np.sqrt(lo * hi), 1.52, "range allowed by\nMacquart+2020", ha="center", va="top", fontsize=7,
             color=ps.INK_SECONDARY)
    axb.axhline(F_DIFFUSE_FIDUCIAL, color=ps.INK_SECONDARY, lw=0.9)
    axb.axhline(1.0, color=ps.INK_SECONDARY, lw=0.6, alpha=0.6)
    axb.text(2.9, F_DIFFUSE_FIDUCIAL - 0.02, "census: 84% diffuse", ha="right", va="top", fontsize=7, color=ps.INK_SECONDARY)
    axb.text(2.9, 1.02, "all baryons diffuse", ha="right", va="bottom", fontsize=7, color=ps.INK_SECONDARY)
    axb.set_xscale("log")
    axb.set_xlim(F.min(), F.max())
    axb.set_ylim(0.45, 1.6)
    axb.set_xlabel("F: how clumpy the intergalactic gas is")
    axb.set_ylabel("f_d: fraction of baryons diffuse")
    axb.set_title("B: clumpier web ↔ more gas needed (68%, 95%)", fontsize=9.5)

    # C: forest plot of Omega_b h70.
    axc = fig.add_subplot(gs[1, 1])
    ob = mr["omega_b_h70_planck"]
    rows = [(FIT_LABELS.get(k, v["label"]), v["omega_b_h70_if_fd_fiducial"]) for k, v in mr["fits"].items()]
    lit = [
        ("Macquart+2020, first localized FRBs (95%)", (0.051, 0.051 - 0.025, 0.051 + 0.021), True),
        ("Connor+2025, DSA-110 + literature", (0.051, 0.045, 0.057), False),
    ]
    y = 0
    labels = []
    for label, q in rows:
        colour = ps.BLUE if y == 0 else ps.INK_SECONDARY
        axc.plot([q["q025"], q["q975"]], [y, y], color=colour, lw=1.0)
        axc.plot([q["q160"], q["q840"]], [y, y], color=colour, lw=3.0 if y == 0 else 2.0, solid_capstyle="round")
        axc.plot([q["q500"]], [y], marker="o", ms=6 if y == 0 else 4.5, color=colour, mec=ps.SURFACE, mew=1.0)
        labels.append(label)
        y -= 1
    y -= 0.4
    for label, (mid, lo_, hi_), is95 in lit:
        axc.plot([lo_, hi_], [y, y], color=C_CONTEXT, lw=1.0 if is95 else 2.0)
        axc.plot([mid], [y], marker="s", ms=4.5, color=C_CONTEXT, mec=ps.SURFACE, mew=1.0)
        labels.append(label)
        y -= 1
    axc.axvline(ob, color=ps.INK, lw=1.0)
    axc.text(ob, 0.9, " Planck 2018\n (CMB)", fontsize=7, color=ps.INK_SECONDARY, va="bottom", ha="center")
    ticks = [0, -1, -2, -3, -4, -5, -6][: len(rows)] + [-len(rows) - 0.4 - i for i in range(len(lit))]
    axc.set_yticks(ticks)
    axc.set_yticklabels(labels, fontsize=7)
    axc.set_ylim(ticks[-1] - 0.7, 1.9)
    axc.set_xlim(0.02, 0.095)
    axc.grid(axis="y", visible=False)
    axc.set_xlabel("Ω_b h₇₀  (assuming 84% of baryons diffuse)")
    axc.set_title("C: the baryon density, 68% (thick) and 95% (thin)", fontsize=9.5)

    fd = bj["f_d"]
    fig.suptitle(
        f"Fast radio bursts weigh the missing baryons: diffuse fraction f_d = {fd['q500']:.2f} "
        f"(+{fd['q840'] - fd['q500']:.2f} / −{fd['q500'] - fd['q160']:.2f})",
        x=0.06, ha="left", y=0.975)
    ps.caption(fig, (
        "Localized FRBs with secure hosts (PATH P(O|x) ≥ 0.9) and z ≥ 0.01 from the FRBs/FRB community compilation; error bars "
        f"show the ±{base.sigma_mw:.0f} pc cm⁻³ Milky Way uncertainty the fit allows for. Model: Macquart et al. (2020) cosmic "
        "DM distribution plus a lognormal host term, Planck 2018 cosmology, all parameters marginalised on a grid; the band is "
        "posterior-predictive. C: our f_d converted to Ω_b h₇₀ assuming f_d = 0.844; the grey rows below the baseline change one "
        "assumption each. Selection effects of the surveys are not modelled."
    ))
    return ps.save(fig, _out("05_macquart"))


def make_all(results: dict, art: dict) -> list:
    written = []
    written += figure_sky(results, art)
    if "dispersion" in art:
        written += figure_dispersion(results, art)
    written += figure_dm_excess(results, art)
    if "dispersion" in art:
        written += figure_repeaters(results, art)
    written += figure_macquart(results, art)
    return written
