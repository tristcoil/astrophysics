"""Every figure of the project, from the analysis results and artefacts.

Each function takes what analysis.py produced and writes PNG + SVG files into
figures/. No figure recomputes science: the numbers quoted on the figures are
the ones in results.json.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

import galaxy as gx
import photometry as ph
import plotstyle as ps
from config import FIGURE_DIR, GD1_PHI1_RANGE, GD1_PHI2_RANGE, M92, REFERENCE

ps.use_style()
DIV = ps.diverging_cmap()

PHI1_LABEL = "φ₁ — along the stream (deg)"
PHI2_LABEL = "φ₂ — across (deg)"


def _out(name):
    return FIGURE_DIR / name


def _sky_axes(ax, ylabel=True):
    ax.set_xlim(*GD1_PHI1_RANGE)
    ax.set_ylim(*GD1_PHI2_RANGE)
    ax.set_aspect("equal")
    ax.grid(False)
    if ylabel:
        ax.set_ylabel(PHI2_LABEL)


def _hist_sky(phi1, phi2, b=0.25):
    e1 = np.arange(GD1_PHI1_RANGE[0], GD1_PHI1_RANGE[1] + b / 2, b)
    e2 = np.arange(GD1_PHI2_RANGE[0], GD1_PHI2_RANGE[1] + b / 2, b)
    H, _, _ = np.histogram2d(phi1, phi2, bins=[e1, e2])
    return e1, e2, H


def _band_polygon(track, phi1_range, off_lo, off_hi):
    x = np.linspace(*phi1_range, 60)
    y0 = track(x)
    return np.column_stack([np.concatenate([x, x[::-1]]),
                            np.concatenate([y0 + off_lo, (y0 + off_hi)[::-1]])])


# ==========================================================================
# GD-1
# ==========================================================================


def gd1_emerges(res, art):
    """Before and after: the stream appears out of a smooth sea of stars."""
    import gd1 as g1

    st, steps = art["stars"], art["steps"]
    e1, e2, H_all = g1.load_countmap()
    fig, axes = plt.subplots(3, 1, figsize=(12, 7.6), sharex=True)
    n_all = int(H_all.sum())

    ax = axes[0]
    ax.pcolormesh(e1, e2, np.log10(np.maximum(H_all, 1)).T, cmap=ps.CMAP_CONTEXT,
                  rasterized=True)
    ax.set_title(f"1. Every star Gaia sees here, brighter than G = 20.5: "
                 f"{n_all / 1e6:.1f} million", fontsize=10)
    _sky_axes(ax)
    ax.text(-79, 4.4, "denser towards the Galactic plane", fontsize=8, color=ps.SURFACE)

    ax = axes[1]
    pm = steps["proper_motion"]
    _, _, H = _hist_sky(st["phi1"][pm], st["phi2"][pm])
    ax.pcolormesh(e1, e2, H.T, cmap=ps.CMAP_CONTEXT, vmax=np.percentile(H, 99.7),
                  rasterized=True)
    ax.set_title(f"2. …moving across the sky like the stream, and far away "
                 f"(proper motion + parallax): {pm.sum():,}", fontsize=10)
    _sky_axes(ax)

    ax = axes[2]
    m = art["members"]
    ax.scatter(st["phi1"][m], st["phi2"][m], s=1.2, c=ps.C_OURS, lw=0, rasterized=True)
    ext = art["extent"]
    ax.axvspan(ext[1], GD1_PHI1_RANGE[1], color=ps.C_MASK, alpha=0.6, lw=0)
    ax.text(0.5 * (ext[1] + GD1_PHI1_RANGE[1]), -1, "no detection beyond here", rotation=90,
            ha="center", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_title(f"3. …and on the colour–magnitude sequence of an ancient cluster "
                 f"at the stream's distance: {m.sum():,}", fontsize=10)
    _sky_axes(ax)
    ax.set_xlabel(PHI1_LABEL)
    ref = REFERENCE["pwb2018"]
    lab = dict(fontsize=8, color=ps.INK_SECONDARY)
    ps.annotate(ax, "the spur", (-33, 1.6), (-40, 4.6), annotation_clip=False)
    ps.annotate(ax, "a gap", (-20.5, 0.3), (-24, -4.5), annotation_clip=False)
    ps.annotate(ax, "densest part: the lost\ncluster's last remains?", (-14, -0.4), (-10, -6.0),
                annotation_clip=False)
    ps.annotate(ax, "the stream fades into\nthe crowded Galactic plane", (-72, -0.4),
                (-79, -6.2), annotation_clip=False)
    fig.suptitle(f"GD-1: a star cluster torn into a stream {ext[1] - ext[0]:.0f}° long, picked "
                 f"out of {n_all / 1e6:.1f} million stars by motion and colour", x=0.06, ha="left",
                 y=0.995)
    fig.tight_layout(h_pad=0.6)
    ps.caption(fig, "Gaia DR3, sky coordinates rotated so the stream runs horizontally "
               "(Koposov et al. 2010 frame). Panel 1: server-side counts, 0.25° bins. "
               "Panels 2-3: stars in the downloaded sample (parallax < 1 mas, BP-RP < 1.6); "
               "panel 3 plots every selected star individually. The track, motion and distance "
               "were all found in the data; see gd1_02 and gd1_03.")
    return ps.save(fig, _out("gd1_01_stream_emerges"))


def gd1_method(res, art):
    tpl, st, track = art["tpl"], art["stars"], art["track"]
    fig, axes = plt.subplots(1, 4, figsize=(17.5, 5.2),
                             gridspec_kw={"width_ratios": [1, 1.15, 1.3, 1]})

    # (a) the template
    ax = axes[0]
    ax.scatter(tpl.members_c0, tpl.members_absmag, s=1, c=ps.C_CONTEXT, lw=0, rasterized=True)
    mm = np.linspace(*tpl.absmag_range, 300)
    ax.plot(tpl.colour_at(mm), mm, color=ps.C_OURS, lw=1.6, label="median sequence")
    ax.fill_betweenx(mm, tpl.colour_at(mm) - 2 * tpl.width_at(mm),
                     tpl.colour_at(mm) + 2 * tpl.width_at(mm), color=ps.C_OURS, alpha=0.15,
                     lw=0, label="accepted (±2σ)")
    ax.plot(tpl.hb_colour, tpl.hb_absmag, color=ps.C_OURS, lw=1.6)
    ax.set_xlim(-0.4, 1.5)
    ax.set_ylim(6.2, -3.4)
    ax.set_xlabel("colour (BP−RP)₀")
    ax.set_ylabel("absolute magnitude M_G")
    ax.set_title("a. The filter: globular cluster M92", fontsize=10)
    ax.text(0.03, 0.30, f"{tpl.n_members:,} M92 stars\n[Fe/H] = {M92.feh}, 8.50 kpc:\n"
            "same instrument,\nsame noise as GD-1", transform=ax.transAxes, fontsize=7.5,
            color=ps.INK_SECONDARY, va="center")
    ax.legend(loc="upper right", fontsize=7)

    # (b) the proper-motion knot in one segment
    ax = axes[1]
    key = min(track.maps, key=lambda k: abs(0.5 * (k[0] + k[1]) + 45))
    x, y, snr = track.maps[key]
    lim = np.nanpercentile(np.abs(snr), 99.9)
    im = ax.pcolormesh(x, y, snr.T, cmap=DIV, vmin=-lim, vmax=lim, rasterized=True)
    ax.grid(False)
    ax.set_xlabel("μφ₁, along the stream (mas/yr)")
    ax.set_ylabel("μφ₂, across (mas/yr)")
    ax.set_title(f"b. Motion search at φ₁ {key[0]:.0f}…{key[1]:.0f}°", fontsize=10)
    k = [kk for kk in track.knots if kk["pass"] == 1 and abs(kk["phi1"] - 0.5 * sum(key)) < 1e-6][0]
    ps.annotate(ax, f"GD-1: {k['snr']:.0f}σ", (k["pm1"], k["pm2"]), (k["pm1"] + 5, k["pm2"] - 5))
    cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cb.set_label("excess on the track (σ)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # (c) the proper-motion track
    ax = axes[2]
    p = np.linspace(*GD1_PHI1_RANGE, 200)
    t1, t2 = track.at(p)
    det = [kk for kk in track.knots if kk["pass"] == 2]
    for comp, t, lab in ((1, t1, "μφ₁ (along)"), (2, t2, "μφ₂ (across)")):
        ok = np.array([kk["detected"] for kk in det])
        ph1 = np.array([kk["phi1"] for kk in det])
        v = np.array([kk[f"pm{comp}"] for kk in det])
        ax.plot(p, t, color=ps.C_OURS, lw=1.4)
        ax.scatter(ph1[ok], v[ok], s=22, color=ps.C_OURS, zorder=3,
                   label="detected (> 6σ)" if comp == 1 else None)
        ax.scatter(ph1[~ok], v[~ok], s=22, facecolor="none", edgecolor=ps.C_OURS, zorder=3,
                   label="below 6σ" if comp == 1 else None)
        ax.text(p[5], t[5] + 0.7, lab, fontsize=8, color=ps.INK_SECONDARY)
    ref = REFERENCE["pwb2018"]
    ra, dec = ref["radec_at_m14"]
    _, _, rp1, rp2, _ = gx.gd1_coordinates(np.array([ra]), np.array([dec]),
                                           np.array([ref["pm_radec_at_m14"][0]]),
                                           np.array([ref["pm_radec_at_m14"][1]]))
    rphi1, _ = gx.gd1_coordinates(np.array([ra]), np.array([dec]))
    ax.scatter([rphi1[0]] * 2, [rp1[0], rp2[0]], marker="D", s=30, color=ps.C_REFERENCE,
               zorder=4, label="Price-Whelan & Bonaca 2018 (φ₁ −18…−10)")
    ax.set_xlabel(PHI1_LABEL)
    ax.set_ylabel("proper motion (mas/yr)")
    ax.set_title("c. The stream's motion, segment by segment", fontsize=10)
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.42), fontsize=7)
    ax.set_xlim(*GD1_PHI1_RANGE)

    # (d) GD-1's own colour-magnitude diagram, in absolute magnitude
    ax = axes[3]
    steps = art["steps"]
    d2 = st["phi2"] - art["phi2_track"](st["phi1"])
    base = steps["proper_motion"] & (st["phi1"] > -60) & (st["phi1"] < -5)
    on = base & (np.abs(d2) < 0.75)
    off = base & (np.abs(d2) > 2.5) & (np.abs(d2) < 5.5)
    absmag = st["g0"] - art["dm_of_phi1"](st["phi1"])
    ce = np.arange(-0.4, 1.5, 0.04)
    me = np.arange(-3.4, 6.3, 0.2)
    Hon, _, _ = np.histogram2d(st["c0"][on], absmag[on], bins=[ce, me])
    Hoff, _, _ = np.histogram2d(st["c0"][off], absmag[off], bins=[ce, me])
    D = Hon - (1.5 / 6.0) * Hoff
    ax.pcolormesh(ce, me, np.maximum(D, 0).T, cmap=ps.CMAP_BLUE, vmax=np.percentile(D, 99.5),
                  rasterized=True)
    ax.plot(tpl.colour_at(mm), mm, color=ps.C_REFERENCE, lw=1.0, label="M92 (panel a)")
    ax.plot(tpl.hb_colour, tpl.hb_absmag, color=ps.C_REFERENCE, lw=1.0)
    ax.set_xlim(-0.4, 1.5)
    ax.set_ylim(6.2, -3.4)
    ax.grid(False)
    ax.set_xlabel("colour (BP−RP)₀")
    ax.set_ylabel("absolute magnitude M_G")
    ax.set_title("d. GD-1's own stars, field subtracted", fontsize=10)
    ax.legend(loc="upper right", fontsize=7)
    ax.text(0.03, 0.97, "stars on the stream minus\nthe area-scaled control bands,\n"
            "shifted by the measured distance", transform=ax.transAxes, fontsize=7.5,
            color=ps.INK_SECONDARY, va="top")
    fig.suptitle("How GD-1 was picked out: an old cluster's colour–magnitude diagram as a "
                 "filter, and the stream's motion found in the data", x=0.02, ha="left")
    fig.tight_layout()
    ps.caption(fig, "b: stars within 0.75° of the stream's midline minus the scaled control "
               "bands 2-5.5° away, counted in 0.75 mas/yr apertures, divided by the Poisson "
               "noise. c: pass-2 searches in 10° segments overlapping by half; the curve is a "
               "cubic fit weighted by significance. Proper motions are as observed (the Sun's "
               "own motion is not removed).")
    return ps.save(fig, _out("gd1_02_how_it_was_found"))


def gd1_distance(res, art):
    ref = REFERENCE["koposov2010"]
    dist = res["distance"]
    model = art["dm_of_phi1"]
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.2), sharex=True,
                             gridspec_kw={"height_ratios": [1.3, 1]})
    ax = axes[0]
    lo, hi = res["distance_model"]["phi1_range"]
    p = np.linspace(lo, hi, 200)
    ax.plot(p, ph.distance_from_modulus(model(p)), color=ps.C_OURS, lw=1.6,
            label="this work: M92 template, smooth fit")
    ph1 = np.array([d["phi1"] for d in dist])
    dm = np.array([d["dm"] for d in dist])
    dme = np.array([d["dm_err"] for d in dist])
    used = (ph1 >= lo - 1e-6) & (ph1 <= hi + 1e-6) & (dme < 0.3)
    dk = ph.distance_from_modulus(dm)
    dk_lo = dk - ph.distance_from_modulus(dm - dme)
    dk_hi = ph.distance_from_modulus(dm + dme) - dk
    ax.errorbar(ph1[used], dk[used], yerr=[dk_lo[used], dk_hi[used]], fmt="o", color=ps.C_OURS,
                ms=5, capsize=0, label="this work: 10° segments (bootstrap 1σ)")
    ax.errorbar(ph1[~used], dk[~used], yerr=[dk_lo[~used], dk_hi[~used]], fmt="o", mfc="none",
                color=ps.C_OURS, ms=5, capsize=0, label="stream too faint to fit")
    ax.errorbar(ref["distance_phi1"], ref["distance_kpc"], yerr=ref["distance_err_kpc"], fmt="s",
                color=ps.C_REFERENCE, ms=5, capsize=0, label="Koposov et al. 2010")
    ax.set_ylabel("distance from the Sun (kpc)")
    ax.set_ylim(5.5, 14)
    ax.set_title("a. Distance along the stream: the far end is ~25% further away", fontsize=10)
    ax.legend(loc="upper left", fontsize=7.5)

    ax = axes[1]
    ext = art["extent"]
    p_all = np.linspace(ext[0], ext[1], 300)
    rows = res["parallax_check"]
    x = np.array([r["phi1"] for r in rows])
    y = np.array([r["parallax"] for r in rows])
    e = np.array([r["parallax_err"] for r in rows])
    ax.fill_between(p, 1 / ph.distance_from_modulus(model(p)) - 0.01,
                    1 / ph.distance_from_modulus(model(p)) + 0.01, color=ps.C_OURS, alpha=0.12,
                    lw=0, label="±0.01 mas: Gaia's small-scale parallax systematics")
    ax.plot(p, 1 / ph.distance_from_modulus(model(p)), color=ps.C_OURS, lw=1.6,
            label="predicted from the photometric distance (panel a)")
    for side in (p_all[p_all < lo], p_all[p_all > hi]):
        if len(side):
            ax.plot(side, 1 / ph.distance_from_modulus(model(side)), color=ps.C_OURS, lw=0.8,
                    alpha=0.6)
    ax.errorbar(x, y, yerr=e, fmt="D", color=ps.C_REFERENCE, ms=5, capsize=0,
                label="mean Gaia parallax of members, field-subtracted")
    ax.set_ylabel("parallax (mas)")
    ax.set_xlabel(PHI1_LABEL)
    ax.set_xlim(*GD1_PHI1_RANGE)
    pc = res["parallax_summary"]
    ax.set_title(f"b. An independent check: parallax, no stellar models — "
                 f"agreement {pc['mean_offset_mas']:+.3f} ± {pc['mean_offset_err_mas']:.3f} mas "
                 "overall", fontsize=10)
    ax.legend(loc="upper left", fontsize=7.5)
    fig.tight_layout()
    ps.caption(fig, "a: distance moduli from sliding M92's sequence (8.50 kpc; Baumgardt & "
               "Vasiliev 2021) over GD-1's field-subtracted stars. Koposov et al. used theoretical "
               "isochrones for [Fe/H] = -1.4 and reduced their distances by 10% by hand; they "
               "noted an M92-based fiducial would be more consistent. b: Lindegren et al. (2021) "
               "zero-point applied; each point averages 50-230 stars weighted by 1/σ².")
    return ps.save(fig, _out("gd1_03_distance"))


def gd1_density(res, art):
    dens = res["density"]
    x = np.array(dens["phi1"])
    y = np.array(dens["per_degree"])
    e = np.array(dens["err"])
    ext = art["extent"]
    inside = (x >= ext[0]) & (x <= ext[1])
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.8), sharex=True,
                             gridspec_kw={"height_ratios": [1.5, 1]})
    ax = axes[0]
    ref, bon = REFERENCE["pwb2018"], REFERENCE["bonaca2019"]
    for (a, b), lab in ((ref["spur_phi1"], "spur"), (ref["blob_phi1"], "blob")):
        ax.axvspan(a, b, color=ps.C_REFERENCE, alpha=0.10, lw=0)
        ax.text(0.5 * (a + b), 37.5, lab, ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    for g in bon["gaps_phi1"]:
        ax.plot([g, g], [34.2, 36.2], color=ps.C_REFERENCE, lw=2.2, solid_capstyle="butt")
    ax.text(bon["gaps_phi1"][0] - 1.0, 35.2, "published gaps G-40, G-20\n(Bonaca et al. 2019)",
            ha="right", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.fill_between(x[inside], (y - e)[inside], (y + e)[inside], color=ps.C_OURS, alpha=0.2,
                    lw=0)
    ax.plot(x[inside], y[inside], color=ps.C_OURS, lw=1.6, label="GD-1 members, field subtracted")
    ax.plot(x[inside], np.array(dens["background"])[inside], color=ps.C_CONTEXT, lw=1.0,
            label="field stars that pass every cut (subtracted)")
    for gp in res["gaps"]:
        if gp["significant"]:
            ps.annotate(ax, f"gap {gp['phi1']:.0f}°\n{100 * gp['depth_fraction']:.0f}% deep, "
                        f"{gp['sigma']:.1f}σ", (gp["phi1"], gp["density"]),
                        (gp["phi1"] - 3.5, max(gp["density"] - 7.5, -2.6)), ha="center")
        else:
            ax.text(gp["phi1"], gp["density"] - 2.6, f"{gp['sigma']:.1f}σ", ha="center",
                    va="top", fontsize=7, color=ps.INK_SECONDARY)
    top = max(res["peaks"], key=lambda q: q["density"]) if res["peaks"] else None
    if top:
        ps.annotate(ax, f"peak {top['phi1']:.1f}° (PWB2018's\nprogenitor candidate: −13.5°)",
                    (top["phi1"], top["density"]), (top["phi1"] + 3, top["density"] + 0.5))
    ax.set_ylabel("stars per degree")
    ax.set_ylim(-3, 40)
    ax.set_title("a. How many stars along the stream: gaps and clumps", fontsize=10)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.8), fontsize=7.5)

    ax = axes[1]
    tr = res["phi2_track"]
    ax.errorbar(tr["centres"], tr["median"], yerr=tr["err"], fmt="o", color=ps.C_OURS, ms=4,
                capsize=0, label="this work: median φ₂ in 4° windows")
    p = np.linspace(tr["centres"][0], tr["centres"][-1], 300)
    ax.plot(p, art["phi2_track"](p), color=ps.C_OURS, lw=1.4)
    k = REFERENCE["koposov2010"]
    ax.errorbar(k["track_phi1"], k["track_phi2"], yerr=k["track_phi2_err"], fmt="s",
                color=ps.C_REFERENCE, ms=4, capsize=0, label="Koposov et al. 2010")
    ax.set_ylabel("φ₂ of the midline (deg)")
    ax.set_xlabel(PHI1_LABEL)
    ax.set_ylim(-1.9, 1.0)
    ax.set_title("b. The stream's midline bends by ~1.4° over its length", fontsize=10)
    ax.legend(loc="lower center", fontsize=7.5, ncol=2)
    ax.set_xlim(*GD1_PHI1_RANGE)
    fig.tight_layout()
    ps.caption(fig, "a: members within 0.75° of the midline, in 3° windows stepped by 1°, minus "
               "members 2.5-5.5° from it scaled by area (±1σ Poisson band). Gaps are labelled by how far "
               "they lie below the lower of the highest points within 10° on either side (small "
               "labels: 2-3σ, marginal). "
               "Ripples at φ₁ < -60° partly trace Gaia's scanning pattern (PWB2018). Shaded: "
               "PWB2018's spur and blob fields.")
    return ps.save(fig, _out("gd1_04_density_and_track"))


def gd1_spur(res, art):
    st, m, track = art["stars"], art["members"], art["phi2_track"]
    fig, ax = plt.subplots(figsize=(11, 4.4))
    lo, hi = -62, -2
    s = m & (st["phi1"] > lo) & (st["phi1"] < hi)
    ax.scatter(st["phi1"][s], st["phi2"][s] - track(st["phi1"][s]), s=6, color=ps.C_OURS, lw=0,
               rasterized=True, label="GD-1 members")
    ref = REFERENCE["pwb2018"]
    boxes = [("spur", ref["spur_phi1"], ref["spur_phi2"], res["spur"]),
             ("mirror", ref["spur_phi1"], (-ref["spur_phi2"][1], -ref["spur_phi2"][0]), None),
             ("blob", ref["blob_phi1"], (-ref["spur_phi2"][1], -ref["spur_phi2"][0]), res["blob"]),
             ("mirror", ref["blob_phi1"], ref["spur_phi2"], None),
             ("control", ref["control_phi1"], ref["spur_phi2"], None),
             ("control", ref["control_phi1"], (-ref["spur_phi2"][1], -ref["spur_phi2"][0]), None)]
    for name, (a, b), (c, d), r in boxes:
        col = ps.C_REFERENCE if name in ("spur", "blob") else ps.C_CONTEXT
        ax.add_patch(Rectangle((a, c), b - a, d - c, fill=False, edgecolor=col, lw=1.2))
    sp_, bl, sc = res["spur"], res["blob"], res["spur_control"]
    ax.text(-33, 2.05, f"spur: {sp_['n_feature']} stars vs {sp_['n_mirror']} in the mirror box",
            ha="center", fontsize=8, color=ps.INK)
    ax.text(-15, -2.35, f"blob: {bl['n_feature']} vs {bl['n_mirror']}", ha="center", fontsize=8,
            color=ps.INK)
    ax.text(-48, 2.05, f"control: {sc['n_feature']} vs {sc['n_mirror']}", ha="center",
            fontsize=8, color=ps.INK_SECONDARY)
    ax.axhline(0, color=ps.C_CONTEXT, lw=0.8)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-3, 3)
    ax.set_xlabel(PHI1_LABEL)
    ax.set_ylabel("φ₂ relative to the midline (deg)")
    p_sp = sp_["p_value"]
    ax.set_title(f"The spur: stars pulled 1–1.7° off the stream (chance probability "
                 f"{p_sp:.0e}) — a hint of a dark-matter clump's fly-by?", fontsize=10)
    fig.tight_layout()
    ps.caption(fig, "Boxes follow Price-Whelan & Bonaca (2018): spur at φ₁ -36…-30°, blob "
               "at -18…-12°, control at -51…-45°, each compared with the same box mirrored "
               "across the stream (same area, same field density). Vertical scale stretched "
               "10x relative to the horizontal.")
    return ps.save(fig, _out("gd1_05_spur"))


def gd1_confirmation(res, art):
    st = art["stars"]
    hu, fate = art["huang"], art["huang_fate"]
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.0), sharex=True)
    ax = axes[0]
    m = art["members"]
    ax.scatter(st["phi1"][m], st["phi2"][m], s=1, color=ps.C_CONTEXT, lw=0, rasterized=True,
               label="our members")
    rec = fate == "member"
    ax.scatter(art["huang_phi1"][rec], art["huang_phi2"][rec], s=26, marker="D",
               color=ps.C_REFERENCE, label=f"spectroscopic members we recover ({rec.sum()})",
               zorder=3)
    miss = (fate != "member") & (fate != "outside our field")
    ax.scatter(art["huang_phi1"][miss], art["huang_phi2"][miss], s=26, marker="D",
               facecolor="none", edgecolor=ps.C_REFERENCE,
               label=f"missed ({miss.sum()}: {', '.join(sorted(set(fate[miss])))})", zorder=3)
    _sky_axes(ax)
    ax.set_ylim(-4, 4)
    ax.set_aspect("auto")
    ax.set_title(f"a. Huang et al. (2019) confirmed {len(fate)} GD-1 stars with spectra; "
                 f"{rec.sum()} of the {len(fate) - np.sum(fate == 'outside our field')} in our "
                 "field are in our list", fontsize=10)
    ax.legend(loc="lower left", fontsize=7.5, ncol=3)

    ax = axes[1]
    k = REFERENCE["koposov2010"]
    ax.errorbar(k["rv_phi1"], k["rv_kms"], yerr=k["rv_err_kms"], fmt="s", color=ps.C_REFERENCE,
                ms=4, capsize=0, label="Koposov et al. 2010 (Calar Alto)")
    ax.errorbar(art["huang_phi1"], hu["rv"], yerr=hu["rv_err"], fmt="D", color=ps.C_THIRD, ms=4,
                capsize=0, label="Huang et al. 2019 (SDSS/LAMOST)")
    r = art["has_rv"]
    rc = res["rv_check"]
    pp = np.linspace(-80, 5, 200)
    ax.plot(pp, np.polyval(art["rv_coef"], pp), color=ps.C_CONTEXT, lw=1.0,
            label="smooth fit to the published velocities")
    ax.errorbar(st["phi1"][r], st["radial_velocity"][r], yerr=st["radial_velocity_error"][r],
                fmt="o", color=ps.C_OURS, ms=6, capsize=0,
                label=f"Gaia DR3 itself: our {rc['n']} members bright enough for a velocity "
                      f"({rc['n_consistent']} on the curve)")
    ax.set_ylabel("radial velocity (km/s)")
    ax.set_xlabel(PHI1_LABEL)
    ax.set_title("b. The sixth dimension: GD-1 swings from +250 to −200 km/s along its length",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_xlim(*GD1_PHI1_RANGE)
    fig.tight_layout()
    ps.caption(fig, "Radial velocities were not used anywhere in the selection, so their smooth "
               "trend along φ₁ is an independent confirmation that the selected stars share one "
               "orbit. Heliocentric velocities: the gradient combines the stream's own orbit and "
               "the Sun's motion, seen along lines of sight that turn through 90° of sky.")
    return ps.save(fig, _out("gd1_06_confirmation"))


def gd1_all(res, art):
    written = []
    for f in (gd1_emerges, gd1_method, gd1_distance, gd1_density, gd1_spur, gd1_confirmation):
        written += f(res, art)
    return written


# ==========================================================================
# Gaia-Sausage-Enceladus
# ==========================================================================

from matplotlib.colors import LogNorm  # noqa: E402


def _density(ax, x, y, bins, cmap, weights=None, vmin=None, rasterized=True):
    H, xe, ye = np.histogram2d(x, y, bins=bins, weights=weights)
    H = np.where(H > 0, H, np.nan)
    im = ax.pcolormesh(xe, ye, H.T, cmap=cmap, norm=LogNorm(vmin=vmin or np.nanmin(H),
                                                             vmax=np.nanmax(H)),
                       rasterized=rasterized)
    ax.grid(False)
    return im


def gse_toomre(res, art):
    f = art["fast"]
    from config import HALO_MAX_DISTANCE_KPC, RANDOM_FRACTION, V_PECULIAR_KMS

    near = f["good"] & (f["distance"] < HALO_MAX_DISTANCE_KPC)
    x = f["V"] + V_PECULIAR_KMS[1]
    y = np.hypot(f["U"] + V_PECULIAR_KMS[0], f["W"] + V_PECULIAR_KMS[2])
    # Slow stars come only from the random 2% (each stands for 50); the fast
    # file has every star above 150 km/s. Weighted, the two join seamlessly.
    slow = near & art["rnd"] & (f["speed_lsr"] <= 150)
    fast = near & (f["speed_lsr"] > 150)
    xs = np.concatenate([x[slow], x[fast]])
    ys = np.concatenate([y[slow], y[fast]])
    ws = np.concatenate([np.full(slow.sum(), 1 / RANDOM_FRACTION), np.ones(fast.sum())])
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    im = _density(ax, xs, ys, [np.arange(-520, 200, 5), np.arange(0, 520, 5)], ps.CMAP_BLUE,
                  weights=ws, vmin=1)
    th = np.linspace(0, np.pi, 200)
    for r, lab, ang in ((100, "100", 60), (210, "210 km/s: the halo cut\n(Helmi et al. 2018)", 118)):
        ax.plot(r * np.cos(th), r * np.sin(th), color=ps.INK_SECONDARY, lw=0.9)
        ax.text(r * np.cos(np.radians(ang)) + 4, r * np.sin(np.radians(ang)) + 6, lab,
                fontsize=7.5, color=ps.INK_SECONDARY, ha="left" if ang < 90 else "right")
    ax.axvline(-233.4, color=ps.C_REFERENCE, lw=1.0)
    ax.text(-236, 470, "no rotation at all\n(V = −V_LSR)", ha="right", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ax.text(0, 25, "thin disc", ha="center", fontsize=8, color=ps.SURFACE)
    ax.text(-120, 55, "thick disc", ha="center", fontsize=8, color=ps.INK)
    ax.text(-330, 300, "halo", ha="center", fontsize=9, color=ps.INK, weight="bold")
    ax.set_xlim(-520, 200)
    ax.set_ylim(0, 520)
    ax.set_aspect("equal")
    ax.set_xlabel("V: speed along the Galaxy's rotation, relative to the local disc (km/s)")
    ax.set_ylabel("√(U² + W²): speed in the other two directions (km/s)")
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("stars per 5×5 km/s cell", fontsize=8)
    pct = 100 * res["halo_fraction_random"]
    ax.set_title(f"The Toomre diagram: {pct:.1f}% of stars near the Sun move like the halo",
                 fontsize=10)
    fig.tight_layout()
    ps.caption(fig, f"Gaia DR3 stars within {HALO_MAX_DISTANCE_KPC} kpc with radial velocities and "
               "parallaxes better than 20%: every star faster than 150 km/s (the fast sample) plus "
               "a random 2% of the rest, each of those counted 50 times. Velocities relative to the "
               "Local Standard of Rest (the Sun's peculiar motion from Schönrich et al. 2010 removed).")
    return ps.save(fig, _out("gse_01_toomre"))


def gse_sausage(res, art):
    f, mp = art["fast"], art["mp"]
    from config import HALO_MAX_DISTANCE_KPC, V_SUN_KMS

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.6), sharex=True, sharey=True)
    bins = [np.arange(-450, 451, 7.5), np.arange(-350, 401, 7.5)]
    lim = dict(xlim=(-450, 450), ylim=(-350, 400))

    ax = axes[0, 0]
    s = art["rnd_near"]
    _density(ax, f["v_r"][s], f["v_phi"][s], bins, ps.CMAP_CONTEXT, vmin=1)
    ax.set(**lim)
    ax.set_title(f"a. All stars (random 2%, {s.sum():,}):\nthe disc hides everything", fontsize=10)
    ax.set_ylabel("v_φ, rotation (km/s)")

    ax = axes[0, 1]
    s = art["halo"]
    _density(ax, f["v_r"][s], f["v_phi"][s], bins, ps.CMAP_BLUE, vmin=1)
    th = np.linspace(0, 2 * np.pi, 300)
    v_lsr = V_SUN_KMS[1] - 12.24
    ax.plot(210 * np.cos(th), v_lsr + 210 * np.sin(th), color=ps.INK_SECONDARY, lw=0.9)
    ax.text(0, v_lsr, "removed by\nthe velocity cut", ha="center", va="center", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ax.set(**lim)
    kb = res["halo_kinematic_ellipsoid"]
    ax.set_title(f"b. Halo by velocity (|v − v_LSR| > 210 km/s,\n{s.sum():,} stars): "
                 f"β = {kb['beta']:.2f}", fontsize=10)

    ax = axes[0, 2]
    ax.axis("off")
    ax.set_title("How to read these panels", fontsize=10)
    ax.text(0.0, 0.98,
            "x: speed towards (+) or away from (−)\nthe Galactic Centre.\n"
            "y: speed along the Galaxy's rotation.\n\n"
            "Disc stars circle at ~230 km/s (top).\n"
            "The halo barely rotates but plunges in\nand out at up to 400 km/s: a long\n"
            "horizontal 'sausage' -- orbits inherited\nfrom one dwarf galaxy that fell in on a\n"
            "near-radial path ~10 billion years ago.\n\n"
            "Bottom row: no velocity cut at all, only\niron abundance from Gaia's own spectra.",
            transform=ax.transAxes, fontsize=8.5, va="top", color=ps.INK_SECONDARY,
            linespacing=1.35)

    for ax, r, letter in zip(axes[1], res["beta_display_bins"], "cde"):
        lo, hi = r["feh_lo"], r["feh_hi"]
        s = art["chem_mp"] & (mp["mh"] >= lo) & (mp["mh"] < hi)
        # Few enough stars to show individually: a density map would break
        # into isolated pixels.
        ax.scatter(mp["v_r"][s], mp["v_phi"][s], s=2.5 if s.sum() < 5000 else 1, lw=0,
                   color=ps.C_OURS, alpha=0.7 if s.sum() < 5000 else 0.35, rasterized=True)
        ax.set(**lim)
        ax.set_title(f"{letter}. {lo} < [M/H] < {hi}, no velocity cut\n{s.sum():,} giants; halo "
                     f"β = {r['beta']:.2f} ± {r['beta_err']:.2f}", fontsize=10)
        ax.set_xlabel("v_r, radial (km/s)")
    axes[1, 0].set_ylabel("v_φ, rotation (km/s)")
    fig.tight_layout()
    ps.caption(fig, "Galactocentric spherical velocities of Gaia DR3 stars with radial velocities. "
               "a-b: within 2.5 kpc, parallax better than 20%. c-e: giants (4000-5500 K, log g < 3.5) "
               "within 3 kpc with medium-quality GSP-Spec metallicities (Recio-Blanco et al. 2023, "
               "calibrated) and parallaxes better than 10%. β from a disc + halo mixture model "
               "(halo: no net radial or polar motion), bootstrap errors.")
    return ps.save(fig, _out("gse_02_sausage"))


def gse_beta(res, art):
    rows = [r for r in res["beta_by_metallicity"] if "beta" in r]
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.8), sharex=True,
                             gridspec_kw={"height_ratios": [1.5, 1]})
    ax = axes[0]
    b = REFERENCE["belokurov2018"]
    ax.fill_between([b["feh_transition"], -0.7], 0.85, 0.95, color=ps.C_REFERENCE, alpha=0.18,
                    lw=0, label="Belokurov et al. 2018 (SDSS dwarfs, 1-10 kpc)")
    ax.fill_between([-3.0, b["feh_transition"]], *b["beta_metal_poor_halo"], color=ps.C_REFERENCE,
                    alpha=0.18, lw=0)
    x = np.array([0.5 * (r["feh_lo"] + r["feh_hi"]) for r in rows])
    xerr = np.array([0.5 * (r["feh_hi"] - r["feh_lo"]) for r in rows])
    y = np.array([r["beta"] for r in rows])
    e = np.array([r["beta_err"] for r in rows])
    ok = np.array([r["weight"] > 0.05 for r in rows])
    ax.errorbar(x[ok], y[ok], xerr=xerr[ok], yerr=e[ok], fmt="o", color=ps.C_OURS, ms=6, capsize=0,
                label="this work: Gaia DR3 giants, no velocity cut")
    ax.errorbar(x[~ok], y[~ok], xerr=xerr[~ok], yerr=e[~ok], fmt="o", mfc="none",
                color=ps.C_OURS, ms=6, capsize=0, label="halo < 5% of the stars: unreliable")
    kb = res["halo_kinematic_ellipsoid"]
    ax.axhline(kb["beta"], color=ps.C_CONTEXT, lw=1.0)
    ax.text(-2.95, kb["beta"] + 0.02, f"halo picked by velocity instead: β = {kb['beta']:.2f}",
            fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_ylabel("anisotropy β")
    ax.set_ylim(-0.1, 1.05)
    ax.set_title("a. Halo orbits turn radial above [Fe/H] ≈ −1.7: the sausage's metallicity",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=7.5)
    ax.text(-2.95, 0.02, "β = 1: purely radial orbits;  0: isotropic;  < 0: circular (disc)",
            fontsize=7.5, color=ps.INK_SECONDARY)

    ax = axes[1]
    comp = np.array([r["cut_completeness"] for r in rows])
    pur = np.array([r["cut_purity"] for r in rows])
    ax.plot(x[ok], 100 * comp[ok], "o-", color=ps.C_OURS,
            label="completeness: halo stars the cut keeps")
    ax.plot(x[ok], 100 * pur[ok], "s-", color=ps.C_CONTEXT,
            label="purity: kept stars that are halo")
    ax.set_ylabel("per cent")
    ax.set_ylim(0, 105)
    ax.set_xlabel("iron abundance [M/H] (0 = the Sun; −1 = one tenth of the Sun's)")
    ax.set_title("b. What the usual velocity cut (> 210 km/s) would have done in each bin",
                 fontsize=10)
    ax.legend(loc="lower left", fontsize=7.5)
    ax.set_xlim(-3.05, 0.55)
    fig.tight_layout()
    ps.caption(fig, "a: β = 1 - (σ_θ² + σ_φ²)/(2σ_r²) of the halo component of a disc + halo "
               "Gaussian mixture, fitted in bins of calibrated GSP-Spec [M/H] (bootstrap 1σ). "
               "Metal-rich bins use the random 2% sample. Belokurov et al. give β ~ 0.9 for "
               "[Fe/H] > -1.7 and 0.2-0.4 at the metal-poor end (bands).")
    return ps.save(fig, _out("gse_03_anisotropy"))


def gse_integrals(res, art):
    f = art["fast"]
    s = art["halo5"]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))
    ref = REFERENCE

    ax = axes[0]
    im = _density(ax, f["L_z"][s], f["E_h17"][s] / 1e5, [np.linspace(-4500, 4500, 150),
                                                          np.linspace(-2.2, -0.6, 120)],
                  ps.CMAP_BLUE, vmin=1)
    lz = ref["helmi2018"]["gse_lz"]
    emin = ref["helmi2018"]["gse_energy_min"] / 1e5
    ax.add_patch(Rectangle((lz[0], emin), lz[1] - lz[0], -0.6 - emin, fill=False,
                           edgecolor=ps.C_REFERENCE, lw=1.4))
    n = res["substructure_counts_5kpc"]["GSE (Helmi+18)"]
    ax.text(lz[0] + 60, -0.66, f"Gaia-Enceladus\n(Helmi et al. 2018 box)\n{n:,} stars",
            fontsize=7.5, va="top", color=ps.INK)
    ax.set_xlabel("L_z, angular momentum about the Galaxy's axis (kpc km/s)")
    ax.set_ylabel("orbital energy E (10⁵ km²/s²)")
    ax.set_title("a. Energy vs angular momentum\n(Helmi et al. 2017 potential)", fontsize=10)
    ax.text(3900, -2.12, "prograde\n(with the disc) →", ha="right", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ax.text(-4400, -2.12, "← retrograde", fontsize=7.5, color=ps.INK_SECONDARY)

    ax = axes[1]
    _density(ax, f["L_z"][s], f["L_perp"][s], [np.linspace(-4500, 4500, 150),
                                              np.linspace(0, 4500, 120)], ps.CMAP_BLUE, vmin=1)
    k = ref["koppelman2019"]
    ax.add_patch(Rectangle((k["hs_lz"][0], k["hs_lperp"][0]), k["hs_lz"][1] - k["hs_lz"][0],
                           k["hs_lperp"][1] - k["hs_lperp"][0], fill=False,
                           edgecolor=ps.C_REFERENCE, lw=1.4))
    ax.text(k["hs_lz"][1] + 80, k["hs_lperp"][1], f"Helmi streams\n(Koppelman+19, box A)",
            fontsize=7.5, va="top")
    ax.set_xlabel("L_z (kpc km/s)")
    ax.set_ylabel("L⊥, angular momentum in the plane (kpc km/s)")
    ax.set_title("b. The Helmi streams: a clump at high L⊥", fontsize=10)

    ax = axes[2]
    _density(ax, f["circ_m17"][s], f["E_m17"][s] / 1e5, [np.linspace(-1.05, 1.05, 120),
                                                         np.linspace(-2.1, -0.8, 120)],
             ps.CMAP_BLUE, vmin=1)
    kb = ref["koppelman2019b"]
    counts = res["substructure_counts_5kpc"]
    for key, name in (("gse", "GSE (Koppelman+19b)"), ("sequoia", "Sequoia"),
                      ("thamnos1", "Thamnos 1"), ("thamnos2", "Thamnos 2")):
        (e0, e1), (c0, c1) = kb[key]["energy"], kb[key]["circ"]
        ax.add_patch(Rectangle((c0, e0), c1 - c0, e1 - e0, fill=False, edgecolor=ps.C_REFERENCE,
                               lw=1.2))
        lab = name.split(" (")[0]
        ax.text(c0, e1 + 0.015, f"{lab}: {counts[name]:,}", fontsize=7, va="bottom")
    ax.set_xlabel("circularity (−1 retrograde circular, 0 radial, +1 disc)")
    ax.set_ylabel("orbital energy E (10⁵ km²/s²)")
    ax.set_title("c. Circularity vs energy (McMillan 2017 potential):\n"
                 "the published debris boxes", fontsize=10)
    fig.tight_layout()
    ps.caption(fig, f"{s.sum():,} halo stars (|v - v_LSR| > 210 km/s) within 5 kpc, as in Helmi "
               "et al.'s Extended Data Fig. 1. Energies depend on the assumed potential, so each "
               "published box is drawn in the potential it was defined in. Boxes are selections, "
               "not detections: the counts show how many stars each one catches.")
    return ps.save(fig, _out("gse_04_integrals_of_motion"))


def gse_chemistry(res, art):
    f = art["fast"]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    t = art["thick"]
    _density(ax, f["mh"][t], f["alpha"][t], [np.arange(-2.0, 0.41, 0.04),
                                            np.arange(-0.3, 0.71, 0.02)], ps.CMAP_CONTEXT, vmin=1)
    g = art["gse_sel"]["Koppelman+19b"]
    ax.scatter(f["mh"][g], f["alpha"][g], s=7, color=ps.C_OURS, lw=0, alpha=0.8,
               label=f"Gaia-Enceladus giants ({g.sum()})")
    rows = res["alpha_by_metallicity"]
    xm = [0.5 * (r["mh_lo"] + r["mh_hi"]) for r in rows]
    ax.errorbar(xm, [r["alpha_thick"] for r in rows], yerr=[r["alpha_thick_err"] for r in rows],
                fmt="s-", color=ps.INK_SECONDARY, ms=5, capsize=0, label="thick disc, median")
    ax.errorbar(xm, [r["alpha_gse"] for r in rows], yerr=[r["alpha_gse_err"] for r in rows],
                fmt="o-", color=ps.C_OURS, ms=6, capsize=0, label="Gaia-Enceladus, median")
    ax.set_xlim(-2.0, 0.4)
    ax.set_ylim(-0.3, 0.7)
    ax.set_xlabel("[M/H]")
    ax.set_ylabel("[α/Fe]")
    comb = res["alpha_difference_combined"]
    ax.set_title(f"a. Accreted stars are α-poorer: {comb['value']:+.3f} ± {comb['err']:.3f} dex "
                 "at fixed [M/H]", fontsize=10)
    ax.legend(loc="lower left", fontsize=7.5)
    ax.text(0.35, 0.62, "thick-disc giants\n(grey density)", ha="right", fontsize=7.5,
            color=ps.INK_SECONDARY)

    ax = axes[1]
    bins = np.arange(-2.6, 0.31, 0.1)
    for key, col, lab in (("Koppelman+19b", ps.C_OURS, "near-radial orbits (Koppelman+19b)"),
                          ("Helmi+18", ps.C_CONTEXT, "Helmi+18 box, |v − v_LSR| > 210")):
        m = art["gse_sel"][key]
        ax.hist(f["mh"][m], bins=bins, density=True, histtype="step", lw=1.6, color=col,
                label=f"{lab}: median {res['gse_chemistry'][key]['mh_median']:.2f}")
    fe = REFERENCE["feuillet2021"]
    ax.axvspan(fe["gse_feh_mean"] - fe["gse_feh_dispersion"], fe["gse_feh_mean"] + fe["gse_feh_dispersion"],
               color=ps.C_REFERENCE, alpha=0.12, lw=0)
    ax.axvline(fe["gse_feh_mean"], color=ps.C_REFERENCE, lw=1.0)
    ax.text(fe["gse_feh_mean"] + 0.04, 0.02, "Feuillet+21:\n−1.17 ± 0.34", transform=ax.get_xaxis_transform(),
            fontsize=7.5, color=ps.INK_SECONDARY, va="bottom")
    ax.axvline(REFERENCE["helmi2018"]["gse_feh_peak"], color=ps.C_REFERENCE, lw=1.0)
    ax.text(REFERENCE["helmi2018"]["gse_feh_peak"] - 0.04, 0.02, "Helmi+18\npeak", ha="right",
            transform=ax.get_xaxis_transform(), fontsize=7.5, color=ps.INK_SECONDARY, va="bottom")
    ax.set_xlabel("[M/H]")
    ax.set_ylabel("share of stars (normalised)")
    ax.set_title("b. How metal-rich was the galaxy we ate?", fontsize=10)
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
    fig.tight_layout()
    ps.caption(fig, "Giants (4000-5500 K, log g < 3.5) within 3 kpc with medium-quality GSP-Spec "
               "parameters and calibrated [M/H], [α/Fe] (Recio-Blanco et al. 2023). Thick disc: "
               "prograde giants 80-180 km/s from the LSR. In a small galaxy stars form slowly, so "
               "iron from type Ia supernovae arrives at lower metallicity and [α/Fe] drops earlier "
               "(Nissen & Schuster 2010). GSP-Spec separates the two sequences less than APOGEE "
               "does (Helmi et al. 2018: ~0.1-0.15 dex).")
    return ps.save(fig, _out("gse_05_chemistry"))


def gse_hrd(res, art):
    f = art["fast"]
    h = art["hrd"]
    g = h & art["masks"]["GSE (Helmi+18)"]
    o = h & ~art["masks"]["GSE (Helmi+18)"]
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    ax.scatter(f["bp_rp"][o], f["abs_g"][o], s=2, color=ps.C_CONTEXT, lw=0, rasterized=True,
               label=f"other halo stars ({o.sum():,})")
    ax.scatter(f["bp_rp"][g], f["abs_g"][g], s=2, color=ps.C_OURS, lw=0, rasterized=True,
               label=f"in the Gaia-Enceladus box ({g.sum():,})")
    ax.set_xlim(0.3, 1.6)
    ax.set_ylim(7.0, -3.5)
    ax.set_xlabel("colour BP−RP")
    ax.set_ylabel("absolute magnitude M_G")
    ax.set_title("Two halos in one diagram: the blue sequence\nis the swallowed galaxy's",
                 fontsize=10)
    ax.legend(loc="lower left", fontsize=7.5, markerscale=4)
    fig.tight_layout()
    ps.caption(fig, f"{h.sum():,} halo stars (|v - v_LSR| > 210 km/s, within 2.5 kpc) with parallaxes "
               "better than 10% and E(B-V) < 0.015 (SFD), as in Helmi et al. (2018) Fig. 2c. "
               "The bluer sequence is more metal-poor at the same age; the redder is the heated "
               "old disc (Gaia Collaboration, Babusiaux et al. 2018).")
    return ps.save(fig, _out("gse_06_hr_diagram"))


def gse_all(res, art):
    written = []
    for f in (gse_toomre, gse_sausage, gse_beta, gse_integrals, gse_chemistry, gse_hrd):
        written += f(res, art)
    return written


# ==========================================================================
# The phase spiral
# ==========================================================================


def spiral_snail(res, art):
    ez, ev, maps = art["ez"], art["ev"], art["maps"]
    zc = 0.5 * (ez[1:] + ez[:-1])
    vc = 0.5 * (ev[1:] + ev[:-1])
    centre = (np.abs(zc)[:, None] < 0.6) & (np.abs(vc)[None, :] < 30)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6), sharey=True)
    panels = (("contrast", "a. Star counts", "excess over a fully mixed disc", 1.0, "{:+.0%}"),
              ("dR", "b. Radial velocity", "v_R minus the mixed-disc average (km/s)", 1.0, None),
              ("dphi", "c. Rotation speed", "v_φ minus the mixed-disc average (km/s)", 1.0, None))
    for ax, (key, title, lab, _, fmt) in zip(axes, panels):
        m = maps[key]
        lim = float(np.nanpercentile(np.abs(m[centre]), 97))
        im = ax.pcolormesh(ez, ev, m.T, cmap=DIV, vmin=-lim, vmax=lim, rasterized=True)
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(lab, fontsize=8)
        if fmt:
            from matplotlib.ticker import PercentFormatter

            cb.ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        cb.ax.tick_params(labelsize=7)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("z, height above the Galactic plane (kpc)")
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-50, 50)
        ax.grid(False)
    axes[0].set_ylabel("v_z, vertical velocity (km/s)")
    ps.annotate(axes[0], "outer rim: fewer stars far\nfrom the plane (sample limit)", (0.85, 0),
                (0.3, -44), bbox=dict(boxstyle="round,pad=0.2", fc=ps.SURFACE, ec="none", alpha=0.85))
    fig.suptitle(f"The phase spiral in {res['n_ring_good'] / 1e6:.1f} million stars: the disc near "
                 "the Sun is still ringing from a kick", x=0.02, ha="left")
    fig.tight_layout()
    ps.caption(fig, "b also shows the tilt of the velocity ellipsoid (a four-leaf pattern that does not "
               "wind). Gaia DR3 stars with radial velocities, parallaxes better than 20% and RUWE < 1.4, "
               f"within 0.1 kpc of the Sun's Galactocentric radius (Antoja et al. 2018's selection; "
               f"they had {REFERENCE['antoja2018']['n_ring_stars_dr2']:,} DR2 stars). Each star is "
               "compared with a perfectly phase-mixed disc: stars of the same vertical amplitude "
               "spread evenly around their loop, with the same average velocities (vertical orbits in "
               "the McMillan 2017 potential). Bins of 0.025 kpc x 1.25 km/s with at least 20 stars.")
    return ps.save(fig, _out("spiral_01_snail_shell"))


def spiral_dating(res, art):
    import spiral as sp

    fig = plt.figure(figsize=(16.0, 5.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1, 1.05])

    # (a) the unwound spiral, in the rotation-speed residual
    ax = fig.add_subplot(gs[0])
    u = art["use"]
    from scipy.stats import binned_statistic_2d

    om = art["omega"][u] * 1000
    th = np.degrees(art["theta"][u])
    eo = np.quantile(om, np.linspace(0, 1, 41))
    et = np.linspace(0, 360, 49)
    m = binned_statistic_2d(th, om, art["dphi"][u], "mean", bins=[et, eo]).statistic
    lim = float(np.nanpercentile(np.abs(m), 97))
    ax.pcolormesh(et, eo, m.T, cmap=DIV, vmin=-lim, vmax=lim, rasterized=True)
    ax.grid(False)
    ax.set_xlabel("vertical phase θ_z (deg; 0 = crossing the plane upwards)")
    ax.set_ylabel("vertical frequency Ω_z (rad/Gyr)")
    ax.set_title("a. The spiral unwound (rotation-speed residual):\ntilted stripes, "
                 "tilt = time since the kick", fontsize=10)
    lo, hi = sp.DATING_RANGES_KPC[0]
    ax.text(4, eo[-1] - 0.3, f"orbits reaching {lo}-{hi} kpc from the plane", fontsize=7.5,
            color=ps.INK_SECONDARY, va="top",
            bbox=dict(boxstyle="round,pad=0.2", fc=ps.SURFACE, ec="none", alpha=0.85))

    # (b) the alignment curves
    ax = fig.add_subplot(gs[1])
    a = REFERENCE["antoja2018"]
    ax.axvspan(*a["perturbation_range_myr"], color=ps.C_REFERENCE, alpha=0.12, lw=0,
               label="Antoja et al. 2018: 300-900 Myr")
    t_grid = sp.T_GRID_MYR
    base = res["age_summary"]["baseline"]
    for key, col, lab in (("counts", ps.C_OURS, "star counts"),
                          ("v_phi", ps.C_THIRD, "rotation speed"),
                          ("v_R", ps.YELLOW, "radial velocity")):
        y = np.array(res["baseline_curves"][key]["amplitude"])
        ax.plot(t_grid, y / y[t_grid >= 150].max(), color=col, lw=1.5,
                label=f"{lab}: {base[key]['t_myr']:.0f} Myr")
    ax.set_xlim(0, 1500)
    ax.set_ylim(0, 1.3)
    ax.set_xlabel("trial age of the disturbance (Myr)")
    ax.set_ylabel("how well the spiral lines up (normalised)")
    ax.set_title("b. Dating the kick: three tracers", fontsize=10)
    ax.legend(loc="upper right", fontsize=7.5)
    ax.text(10, 0.05, "t < 150 Myr excluded: any pattern\nthat does not wind adds up there",
            fontsize=7, color=ps.INK_SECONDARY)

    # (c) the systematic spread
    ax = fig.add_subplot(gs[2])
    ax.axvspan(*a["perturbation_range_myr"], color=ps.C_REFERENCE, alpha=0.12, lw=0)
    rows = res["dating_table"]
    ranges = [f"{r0}-{r1} kpc" for r0, r1 in sp.DATING_RANGES_KPC]
    ymap = {r: i for i, r in enumerate(ranges)}
    for key, col, mk in (("counts", ps.C_OURS, "o"), ("v_phi", ps.C_THIRD, "s"),
                         ("v_R", ps.YELLOW, "D")):
        xs = [r["t_myr"] for r in rows if r["tracer"] == key]
        ys = [ymap[f"{r['amp_lo']}-{r['amp_hi']} kpc"] + {"counts": -0.15, "v_phi": 0, "v_R": 0.15}[key]
              for r in rows if r["tracer"] == key]
        ax.scatter(xs, ys, color=col, marker=mk, s=30, zorder=3,
                   label={"counts": "star counts", "v_phi": "rotation speed",
                          "v_R": "radial velocity"}[key])
    alt = res["age_other_potential"]["Helmi et al. (2017)"]
    ax.scatter(list(alt.values()), [len(ranges) + 0.3] * 3, facecolor="none",
               edgecolor=ps.INK_SECONDARY, s=30)
    sd = res["surface_density_11"]
    ax.set_yticks(range(len(ranges) + 1))
    ax.set_yticklabels(ranges + [f"baseline, Helmi+17 disc\n(Σ₁.₁ = {sd['Helmi et al. (2017)']:.0f} M☉/pc²)"],
                       fontsize=7.5)
    s_ = res["age_summary"]
    ax.axvline(s_["median_myr"], color=ps.INK_SECONDARY, lw=1.0)
    ax.set_xlim(0, 1000)
    ax.legend(loc="lower right", fontsize=7, title="orbits reaching (from the plane):",
              title_fontsize=7)
    ax.set_xlabel("age of the disturbance (Myr)")
    ax.set_title(f"c. {s_['min_myr']:.0f}-{s_['max_myr']:.0f} Myr across tracers and orbit "
                 f"ranges\n(median {s_['median_myr']:.0f} Myr)", fontsize=10)
    ax.invert_yaxis()
    fig.tight_layout()
    ps.caption(fig, "Vertical orbits in the McMillan (2017) potential, whose local disc "
               f"(Σ within 1.1 kpc = {sd['McMillan (2017)']:.0f} M☉/pc²) matches the measured "
               f"{SIGMA_TEXT}; the Helmi et al. (2017) disc is too heavy, so its stars oscillate too "
               "fast and the spiral looks younger (hollow circles). Bootstrap errors on each age are a "
               "few Myr: the spread between tracers and orbit ranges is the real uncertainty.")
    return ps.save(fig, _out("spiral_02_dating"))


SIGMA_TEXT = "68 ± 4 (Bovy & Rix 2013)"


def spiral_all(res, art):
    written = []
    for f in (spiral_snail, spiral_dating):
        written += f(res, art)
    return written


# ==========================================================================
# Overview
# ==========================================================================


def overview_arrays(part, res, art):
    """The few small arrays the overview figure needs from each part."""
    if part == "gd1":
        st, m = art["stars"], art["members"]
        return dict(phi1=st["phi1"][m].copy(), phi2=st["phi2"][m].copy())
    if part == "gse":
        mp, r = art["mp"], res["beta_display_bins"][1]
        sel = art["chem_mp"] & (mp["mh"] >= r["feh_lo"]) & (mp["mh"] < r["feh_hi"])
        return dict(v_r=mp["v_r"][sel].copy(), v_phi=mp["v_phi"][sel].copy())
    if part == "spiral":
        return dict(contrast=art["maps"]["contrast"], ez=art["ez"], ev=art["ev"])
    raise ValueError(part)


def overview(results, keep):
    g, h, s = keep["gd1"], keep["gse"], keep["spiral"]
    rg, rh, rs = results["gd1"], results["gse"], results["spiral"]
    fig = plt.figure(figsize=(16.5, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 1, 1])

    ax = fig.add_subplot(gs[0])
    ax.scatter(g["phi1"], g["phi2"], s=1.2, color=ps.C_OURS, lw=0, rasterized=True)
    ax.set_xlim(*GD1_PHI1_RANGE)
    ax.set_ylim(-5, 5)
    ax.set_xlabel(PHI1_LABEL)
    ax.set_ylabel(PHI2_LABEL)
    ax.grid(False)
    ax.set_title(f"1. GD-1: a dissolved star cluster, {rg['extent_phi1'][1] - rg['extent_phi1'][0]:.0f}° "
                 f"of sky\n{len(g['phi1']):,} stars found by motion and colour", fontsize=10)

    ax = fig.add_subplot(gs[1])
    r = rh["beta_display_bins"][1]
    ax.scatter(h["v_r"], h["v_phi"], s=2.5, lw=0, color=ps.C_OURS, alpha=0.7, rasterized=True)
    ax.set_xlim(-450, 450)
    ax.set_ylim(-350, 400)
    ax.set_xlabel("v_r, towards/away from the Galactic Centre (km/s)")
    ax.set_ylabel("v_φ, rotation (km/s)")
    ax.text(0, 290, "the disc", ha="center", fontsize=8, color=ps.INK_SECONDARY)
    ax.text(-400, -60, "the 'sausage': the galaxy we ate", fontsize=8, color=ps.INK_SECONDARY)
    ax.set_title(f"2. Gaia-Sausage-Enceladus: {len(h['v_r']):,} giants with\n"
                 f"{r['feh_lo']} < [M/H] < {r['feh_hi']}, halo β = {r['beta']:.2f}", fontsize=10)

    ax = fig.add_subplot(gs[2])
    mm = s["contrast"]
    zc = 0.5 * (s["ez"][1:] + s["ez"][:-1])
    vc = 0.5 * (s["ev"][1:] + s["ev"][:-1])
    centre = (np.abs(zc)[:, None] < 0.6) & (np.abs(vc)[None, :] < 30)
    lim = float(np.nanpercentile(np.abs(mm[centre]), 97))
    ax.pcolormesh(s["ez"], s["ev"], mm.T, cmap=DIV, vmin=-lim, vmax=lim, rasterized=True)
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-40, 40)
    ax.grid(False)
    ax.set_xlabel("z, height above the plane (kpc)")
    ax.set_ylabel("v_z, vertical velocity (km/s)")
    a = rs["age_summary"]
    ax.set_title(f"3. The phase spiral: the disc still ringing\nfrom a kick "
                 f"{a['min_myr']:.0f}-{a['max_myr']:.0f} Myr ago", fontsize=10)
    fig.suptitle("Galactic archaeology with Gaia DR3: what the Milky Way has eaten, and what "
                 "still shakes it", x=0.02, ha="left")
    fig.tight_layout()
    ps.caption(fig, "All panels: Gaia DR3 only. Details in gd1_*, gse_* and spiral_* figures.")
    return ps.save(fig, _out("00_overview"))


def make_all(results, artefacts):
    written = []
    if "gd1" in results:
        written += gd1_all(results["gd1"], artefacts["gd1"])
    if "gse" in results:
        written += gse_all(results["gse"], artefacts["gse"])
    if "spiral" in results:
        written += spiral_all(results["spiral"], artefacts["spiral"])
    return written
