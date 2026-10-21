"""Static figures, drawn from the analysis artefacts (never by re-running it).

Stars are drawn in their true colours -- the colour of a blackbody at the
temperature implied by each star's Gaia colour -- which is why the figures use
a dark surface: a hot star's pale blue-white would vanish on white paper. All
other colours are the documented dark-mode palette (see plotstyle.py), used in
fixed order, with direct labels wherever a colour carries identity.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

import plotstyle as ps
import stellar_physics as sp
from analysis import CLASS_LABELS, LITERATURE_SUN_HEIGHT_PC
from config import FIGURE_DIR, MAX_DISTANCE_PC, NAKED_EYE_V, NAMED_RADIUS_PC

ps.use_style()
COLOURS = sp.StarColours()

PROVENANCE = ("Data: ESA Gaia DR3, cleaned with the Rybizki et al. (2022) fidelity classifier; "
              "stars too bright for Gaia from Hipparcos (van Leeuwen 2007).")


def _arr(cat, name):
    col = np.asarray(cat[name])
    return col.astype(str) if col.dtype.kind in "SU" else col


def star_colours(cat) -> np.ndarray:
    """RGB for every star; grey where there is no temperature to go on."""
    t = np.asarray(cat["display_teff"], dtype=float)
    rgb = COLOURS(np.where(np.isfinite(t), t, 5800.0))
    grey = np.array([0.54, 0.53, 0.51])
    rgb[~np.isfinite(t)] = grey
    return rgb


def star_sizes(m_g, scale: float = 0.35) -> np.ndarray:
    """Marker area growing with the square root of luminosity.

    Area in the square root of luminosity (not luminosity itself), because a
    giant is ~10^6 times more luminous than a late M dwarf, which no marker
    scale could show honestly. Faint dwarfs become specks, giants become beacons.
    """
    m_g = np.asarray(m_g, dtype=float)
    s = scale * 10 ** (-0.2 * (np.where(np.isfinite(m_g), m_g, 14.0) - 8.0))
    return np.clip(s, 0.04, 45.0)


def _ink_on(fill: str) -> str:
    """White or near-black text for a label sitting inside a coloured fill."""
    from matplotlib.colors import to_rgb

    def lum(c):
        c = np.array(to_rgb(c))
        c = np.where(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
        return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]

    fl = lum(fill)
    white = (1.05) / (fl + 0.05)
    dark = (fl + 0.05) / (lum(ps.SURFACE) + 0.05)
    return ps.INK if white >= dark else ps.SURFACE


def _log_ticks(ax, values, axis="x"):
    ticks = np.log10(values)
    labels = [f"{v:g}" for v in values]
    if axis == "x":
        ax.set_xticks(ticks, labels)
    else:
        ax.set_yticks(ticks, labels)


# --------------------------------------------------------------------------
# 1. The map: from above and edge-on
# --------------------------------------------------------------------------

SLICE_HALF_THICKNESS_PC = 20.0


def _label(ax, text, xy, xytext, **kw):
    """A callout whose text sits on a patch of surface, so star fields cannot swallow it."""
    return ax.annotate(text, xy=xy, xytext=xytext, color=kw.pop("color", ps.INK_SECONDARY), fontsize=8,
                       ha=kw.pop("ha", "left"), va=kw.pop("va", "center"),
                       bbox={"boxstyle": "round,pad=0.3", "fc": ps.SURFACE, "ec": "none", "alpha": 0.9},
                       arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "lw": 0.8, "shrinkA": 0, "shrinkB": 3},
                       zorder=8, **kw)


def figure_projections(cat, res):
    """Two 40 pc-thick slices through the Sun: one in the Galactic plane, one edge-on.

    A projection of the whole 100 pc ball piles 300,000 stars on top of each
    other and hides everything; a slice keeps the crowding honest, and the
    Hyades -- which sits 18 pc below the plane of the Sun -- falls inside it.
    """
    x, y, z = (np.asarray(cat[k], dtype=float) for k in ("x", "y", "z"))
    size = star_sizes(cat["m_g"], scale=0.45)
    rgb = star_colours(cat)
    peaks = {p["match"]["name"]: p for p in res["structure"]["peaks"] if p.get("match")}
    t = SLICE_HALF_THICKNESS_PC

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 7.2), gridspec_kw={"wspace": 0.14})
    panels = [
        (np.abs(z) < t, x, y, "Y (pc)   \u2192  direction of Galactic rotation",
         "In the plane of the Galaxy", f"Looking down on a slice |Z| < {t:g} pc"),
        (np.abs(y) < t, x, z, "Z (pc)   \u2192  toward the North Galactic Pole",
         "Edge-on", f"A vertical slice |Y| < {t:g} pc, seen from the side"),
    ]
    lim = 112
    for k, (ax, (sl, u, v, vlab, title, sub)) in enumerate(zip(axes, panels)):
        order = np.flatnonzero(sl)[np.argsort(size[sl])]
        ax.scatter(u[order], v[order], s=size[order], c=rgb[order], alpha=0.85, linewidths=0, rasterized=True)
        for r in (25, 50, 75, 100):
            ax.add_patch(Circle((0, 0), r, fill=False, color=ps.INK_SECONDARY if r == 100 else ps.MUTED,
                                lw=0.7, zorder=6, alpha=0.8))
            ax.text(r * 0.707 + 1.5, -r * 0.707 - 1.5, f"{r} pc", color=ps.INK_SECONDARY, fontsize=7,
                    ha="left", va="top", zorder=8,
                    bbox={"boxstyle": "round,pad=0.15", "fc": ps.SURFACE, "ec": "none", "alpha": 0.85})
        ax.plot(0, 0, marker="o", ms=9, mfc="none", mec=ps.INK, mew=1.3, zorder=7)
        _label(ax, "the Sun", (1.5, 1.5), (12, 14), va="bottom", color=ps.INK)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_xlabel("X (pc)   \u2192  toward the Galactic centre")
        ax.set_ylabel(vlab)
        ax.set_title(title)
        ps.subtitle(ax, sub)
        found = {"Hyades": "the densest clump found", "Coma Berenices": "the second-densest clump found"}
        for name, p in peaks.items():
            cx, cy, cz = p["centre_xyz_pc"]
            cv = cy if k == 0 else cz
            if abs(cz if k == 0 else cy) > t:
                continue  # this cluster is outside the slice
            ax.add_patch(Circle((cx, cv), 10, fill=False, color=ps.INK, lw=1.0, zorder=7))
            text = f"{name}, {p['distance_pc']:.1f} pc\n{found[name]}"
            if name == "Hyades":
                _label(ax, text, (cx - 7, cv - 7), (-104, -98), va="bottom")
            else:
                _label(ax, text, (cx + 8, cv + 6), (40, 104), va="top")
    # For reference only -- the mid-plane is a literature value, not measured here.
    ax = axes[1]
    ax.plot([-100, 100], [-LITERATURE_SUN_HEIGHT_PC] * 2, color=ps.MUTED, lw=0.8, zorder=0.6)
    ax.text(lim + 2, -LITERATURE_SUN_HEIGHT_PC, "Galactic\nmid-plane\n(literature)", color=ps.MUTED,
            fontsize=7, va="center", ha="left", clip_on=False)

    ps.temperature_key(fig, [0.60, 0.955, 0.30, 0.02], COLOURS)
    n_gaia, n_hip = res["sample"]["n_gaia"], res["sample"]["n_hipparcos"]
    n_slices = int((np.abs(z) < t).sum()), int((np.abs(y) < t).sum())
    fig.suptitle(f"Our corner of the Galaxy: the stars within 100 parsecs (326 light years)",
                 x=0.125, y=1.03, ha="left")
    ps.caption(fig, f"{PROVENANCE} Of {n_gaia + n_hip:,} stars in the catalogue, the slices hold {n_slices[0]:,} "
               f"(left) and {n_slices[1]:,} (right). Marker area grows with the square root of luminosity; colour is "
               "the blackbody colour of each star's temperature. Coma Berenices lies almost straight overhead, "
               "so only the edge-on slice contains it. The mid-plane line is a literature value (Bennett & Bovy "
               "2019), shown for orientation.", y=0.02)
    return ps.save(fig, FIGURE_DIR / "01_the_neighbourhood")


# --------------------------------------------------------------------------
# 2. What is out there, versus what you can see
# --------------------------------------------------------------------------

def figure_census(res):
    """Stars within 20 pc by type: all of them, then only the naked-eye ones."""
    radius = str(int(NAMED_RADIUS_PC))
    counts = res["census"]["by_radius"][radius]["counts"]
    eye = res["census"]["naked_eye"][radius]["counts"]
    classes = [c for c in ("A", "F", "G", "K", "M", "L/T", "WD", "Giant", "No colour")]
    labels = [CLASS_LABELS[c] for c in classes]
    ypos = np.arange(len(classes))[::-1]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, data, title in ((axes[0], counts, f"Every star within {radius} pc"),
                            (axes[1], eye, f"Only those visible to the naked eye (V < {NAKED_EYE_V:g})")):
        vals = np.array([data[c] for c in classes])
        colours = [ps.C_EMPHASIS if c == "M" else ps.C_CONTEXT for c in classes]
        ax.barh(ypos, vals, height=0.62, color=colours, edgecolor=ps.SURFACE, linewidth=2)
        top = vals.max()
        for yy, v, c in zip(ypos, vals, classes):
            ax.text(v + top * 0.012, yy, f"{v:,}", va="center", ha="left", fontsize=8,
                    color=ps.INK if c == "M" else ps.INK_SECONDARY,
                    fontweight="bold" if c == "M" else "normal")
        ax.set_xlim(0, top * 1.14)
        ax.grid(axis="y", visible=False)
        ax.set_title(title)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(ypos, labels)
    total = res["census"]["by_radius"][radius]["total"]
    n_eye = res["census"]["naked_eye"][radius]["total"]
    ps.subtitle(axes[0], f"{total:,} objects")
    ps.subtitle(axes[1], f"{n_eye:,} of them: {100 * n_eye / total:.0f}%")
    axes[0].set_xlabel("number of stars")
    axes[1].set_xlabel("number of stars")
    share = res["census"]["m_dwarf_share_of_stars_20pc"]
    fig.suptitle(f"Red dwarfs are {100 * share:.0f}% of the stars around us, and not one of them "
                 "can be seen without a telescope", x=0.125, y=1.02, ha="left")
    ps.caption(fig, f"{PROVENANCE} Classes from each star's position in the HR diagram (see figure 05). "
               "The L/T count is only the brown dwarfs Gaia can detect, a fraction of the true number "
               "(figures 03 and 06). V magnitudes are measured by Hipparcos where it observed the star, "
               "otherwise estimated from Gaia photometry.")
    return ps.save(fig, FIGURE_DIR / "02_census")


# --------------------------------------------------------------------------
# 3. Selection effects: what Gaia can and cannot see
# --------------------------------------------------------------------------

def figure_selection(cat, res, seq):
    """Absolute magnitude against distance, with Gaia's bright and faint walls."""
    origin = _arr(cat, "origin")
    d = np.asarray(cat["dist_pc"], dtype=float)
    mg = np.asarray(cat["m_g"], dtype=float)
    gaia = (origin == "Gaia DR3") & np.isfinite(mg)
    hip = (origin == "Hipparcos") & np.isfinite(mg)
    lim = res["gaia_limits"]

    fig, ax = plt.subplots(figsize=(10, 7.4))
    ax.hist2d(np.log10(d[gaia]), mg[gaia], bins=[np.linspace(0.0, 2.0, 140), np.linspace(-6, 26, 200)],
              cmap=ps.DENSITY_CMAP, norm=LogNorm(), rasterized=True, cmin=1)
    ax.scatter(np.log10(d[hip]), mg[hip], s=9, color=ps.C_HIPPARCOS, linewidths=0, zorder=4,
               label=f"added from Hipparcos ({hip.sum():,})")

    pts = res["validation"]["census_10pc"]["points"]
    status = np.array(pts["status"])
    kind = np.array(pts["kind"])
    pd, pm = np.array(pts["dist_pc"], dtype=float), np.array(pts["m_g"], dtype=float)
    missed = (status == "missed") & np.isfinite(pm)
    bd = missed & (kind == "brown dwarf")
    other = missed & (kind != "brown dwarf")
    ax.scatter(np.log10(pd[bd]), pm[bd], s=22, marker="v", color=ps.AQUA, edgecolor=ps.SURFACE,
               linewidths=0.6, zorder=5, label=f"brown dwarfs within 10 pc that Gaia missed ({bd.sum()} with a magnitude to plot)")
    ax.scatter(np.log10(pd[other]), pm[other], s=22, marker="o", color=ps.AQUA, edgecolor=ps.SURFACE,
               linewidths=0.6, zorder=5, label=f"stars within 10 pc that Gaia missed ({other.sum()} plotted; all in binaries)")

    dd = np.geomspace(1.0, 100.0, 200)
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(26, -6)
    walls = ((lim["g_faint"], f"faint limit: apparent G = {lim['g_faint']:.1f}", 30.0, 0.5, "top"),
             (lim["g_bright"], f"bright limit: the brightest Gaia star here, G = {lim['g_bright']:.1f}", 1.6, -0.5, "bottom"))
    for g, label, d_at, dy, va in walls:
        m = g - 5 * np.log10(dd / 10)
        ax.plot(np.log10(dd), m, color=ps.INK_SECONDARY, lw=1.0, zorder=3)
        # A line of constant apparent magnitude has slope -5 in (log d, M) space;
        # turn that into the angle it makes on screen, y axis inverted and all.
        x0, y0 = np.log10(d_at), g - 5 * np.log10(d_at / 10)
        angle = ax.transData.transform_angles(np.array([np.degrees(np.arctan2(-5.0, 1.0))]),
                                              np.array([[x0, y0]]))[0]
        ax.text(x0, y0 + dy, label, rotation=angle, rotation_mode="anchor", color=ps.INK_SECONDARY,
                fontsize=7.5, va=va, ha="left", zorder=6)
    ax.fill_between(np.log10(dd), lim["g_faint"] - 5 * np.log10(dd / 10), 27, color=ps.MUTED, alpha=0.10, lw=0, zorder=0)
    ax.fill_between(np.log10(dd), -7, lim["g_bright"] - 5 * np.log10(dd / 10), color=ps.MUTED, alpha=0.10, lw=0, zorder=0)
    ax.text(np.log10(40), 24.3, "too faint for Gaia", color=ps.INK_SECONDARY, fontsize=9, ha="center")
    ax.text(np.log10(2.2), -4.6, "too bright for Gaia", color=ps.INK_SECONDARY, fontsize=9, ha="center")

    # Reference: where typical dwarfs sit, from the Mamajek sequence.
    for spt in ("G2V", "K5V", "M0V", "M5V", "M9V", "L5V"):
        i = int(np.flatnonzero(seq.spt == spt)[0])
        ax.axhline(seq.m_g[i], xmax=0.985, color=ps.GRID, lw=0.7, zorder=0.3)
        ax.text(2.005, seq.m_g[i], spt, color=ps.MUTED, fontsize=7, va="center", ha="left")

    _log_ticks(ax, [1, 2, 5, 10, 20, 50, 100])
    ax.set_xlabel("distance (pc)")
    ax.set_ylabel("absolute magnitude  M$_G$   (brighter ↑)")
    ax.grid(False)
    ax.set_title("Two walls: what Gaia cannot see")
    ps.subtitle(ax, "Every star lies between the diagonals. Below the lower one, an object is too faint to "
                    "register; above the upper one, too bright.")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), handletextpad=0.4)
    ps.caption(fig, f"{PROVENANCE} Blue shading: density of Gaia stars (log scale). Orange points between the walls "
               "are Hipparcos stars inside Gaia's range that Gaia still could not measure: almost all binaries. "
               "Missed objects are from the "
               "10 pc census of Reyle et al. (2021, 2023 update), plotted with their catalogued G or an estimate "
               "from other photometry; the coldest brown dwarfs have no optical magnitude at all and cannot be "
               "plotted. Right-hand labels mark the absolute magnitudes of typical dwarfs (Mamajek table).")
    return ps.save(fig, FIGURE_DIR / "03_selection_effects")


# --------------------------------------------------------------------------
# 4. Density with distance
# --------------------------------------------------------------------------

def figure_density(res):
    """Number density of each class in 10 pc shells."""
    dens = res["density_by_class"]
    edges = np.array(dens["edges_pc"])
    centres = 0.5 * (edges[1:] + edges[:-1])
    series = [("M", ps.SERIES[0]), ("K", ps.SERIES[1]), ("G", ps.SERIES[2]),
              ("WD", ps.SERIES[3]), ("L/T", ps.SERIES[4])]
    fig, ax = plt.subplots(figsize=(10, 6.2))
    for cls, colour in series:
        y = np.array(dens["classes"][cls]["density_per_pc3"])
        e = np.array(dens["classes"][cls]["err_per_pc3"])
        ok = y > 0
        ax.errorbar(centres[ok], y[ok], yerr=e[ok], color=colour, lw=2, marker="o", ms=5,
                    mec=ps.SURFACE, mew=1.5, capsize=0, elinewidth=1, label=CLASS_LABELS[cls])
        ax.text(centres[ok][-1] + 2, y[ok][-1], CLASS_LABELS[cls], color=ps.INK_SECONDARY,
                fontsize=8, va="center")
    lt = np.array(dens["classes"]["L/T"]["density_per_pc3"])
    last = int(np.flatnonzero(lt > 0)[-1])
    lt_fall = lt[0] / lt[last]
    lt_end = edges[last + 1]
    ax.set_yscale("log")
    ax.set_xlim(0, 113)
    ax.set_xticks(np.arange(0, 101, 10))
    ax.set_xlabel("distance from the Sun (pc), 10 pc shells")
    ax.set_ylabel("stars per cubic parsec")
    ax.set_title("How the census changes with distance")
    ps.subtitle(ax, f"Flat means Gaia sees every one. L/T dwarfs fall {lt_fall:,.0f}-fold by {lt_end:.0f} pc, and none "
                    "are seen beyond: Gaia running out of sensitivity, not brown dwarfs running out.")
    ax.legend(loc="lower left", ncol=5)
    ps.caption(fig, f"{PROVENANCE} Poisson error bars. The two detected clusters are masked. The gentle decline "
               "shared by every class has two ordinary causes: at larger distances more binary companions blur "
               "into single Gaia sources, and a larger sphere reaches further out of the thin Galactic disc.")
    return ps.save(fig, FIGURE_DIR / "04_density_with_distance")


# --------------------------------------------------------------------------
# 5. How each star is classified
# --------------------------------------------------------------------------

def figure_hr(cat, res, seq):
    """The HR diagram of the Gaia sample, with every classification boundary drawn."""
    import config as C

    origin = _arr(cat, "origin")
    gaia = origin == "Gaia DR3"
    bp_rp = np.asarray(cat["bp_rp"], dtype=float)[gaia]
    mg = np.asarray(cat["m_g"], dtype=float)[gaia]
    cls = _arr(cat, "cls")[gaia]
    bad = np.asarray(cat["colour_contaminated"], dtype=bool)[gaia]
    ok = np.isfinite(bp_rp) & np.isfinite(mg) & ~bad

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.hist2d(bp_rp[ok], mg[ok], bins=[np.linspace(-0.6, 5.2, 300), np.linspace(-3, 22, 320)],
              cmap=ps.DENSITY_CMAP, norm=LogNorm(), rasterized=True, cmin=1)
    ax.scatter(bp_rp[bad], mg[bad], s=0.6, color=ps.MUTED, alpha=0.5, linewidths=0, rasterized=True,
               label=f"colour contaminated by a neighbour ({bad.sum():,}): reclassified by brightness")

    c = np.linspace(-0.12, 5.1, 300)
    ms = seq.main_sequence_mg(c)
    ax.plot(c, ms, color=ps.INK_SECONDARY, lw=1.0, label="mean main sequence (Mamajek)")
    giant = ms - C.GIANT_OFFSET_MAG
    ax.plot(c, np.where(giant < C.GIANT_MAX_MG, giant, np.nan), color=ps.INK, lw=1.0, ls=(0, (4, 2)))
    wd_c = np.linspace(-0.6, C.WD_MAX_BP_RP, 50)
    ax.plot(wd_c, C.WD_LINE_SLOPE * wd_c + C.WD_LINE_INTERCEPT, color=ps.INK, lw=1.0, ls=(0, (4, 2)))
    ax.plot([C.WD_MAX_BP_RP] * 2, [C.WD_LINE_SLOPE * C.WD_MAX_BP_RP + C.WD_LINE_INTERCEPT, 22],
            color=ps.INK, lw=1.0, ls=(0, (4, 2)))
    ax.axhline(C.ULTRACOOL_MG, xmin=(C.WD_MAX_BP_RP + 0.6) / 5.8, color=ps.INK, lw=1.0, ls=(0, (4, 2)))

    n = {k: int((cls == k).sum()) for k in ("WD", "Giant", "L/T")}
    n_ms = int(np.isin(cls, ["B", "A", "F", "G", "K", "M"]).sum())
    ax.text(0.05, 15.2, f"White dwarfs\n{n['WD']:,}", color=ps.INK, fontsize=9, ha="center")
    ax.text(1.95, 0.9, f"Giants & subgiants\n{n['Giant']:,}", color=ps.INK, fontsize=9)
    ax.text(3.35, 19.6, f"L/T dwarfs\n{n['L/T']:,}", color=ps.INK, fontsize=9, ha="center")
    ax.text(3.45, 9.0, f"Main sequence\n{n_ms:,}", color=ps.INK, fontsize=9, ha="left")
    ax.plot(0.82, 4.67, marker="o", ms=7, mfc="none", mec=ps.YELLOW, mew=1.4, zorder=5)
    ps.annotate(ax, "the Sun", (0.82, 4.67), (-0.35, 6.4), color=ps.INK_SECONDARY)

    # Spectral-class boundaries along the top, converted from temperature to colour.
    bounds = seq.class_boundaries()
    x_of_t = lambda t: float(np.interp(t, seq.teff[::-1][np.isfinite(seq.bp_rp[::-1])],
                                       seq.bp_rp[::-1][np.isfinite(seq.bp_rp[::-1])]))
    edges = [x_of_t(bounds[k][1]) for k in ("F", "G", "K", "M")]
    for e in edges:
        ax.axvline(e, ymin=0.965, color=ps.MUTED, lw=0.8)
    mids = [(-0.35, "A"), ((edges[0] + edges[1]) / 2, "F"), ((edges[1] + edges[2]) / 2, "G"),
            ((edges[2] + edges[3]) / 2, "K"), (3.2, "M")]
    for xm, lab in mids:
        ax.text(xm, -2.5, lab, color=ps.INK_SECONDARY, fontsize=9, ha="center", va="center", fontweight="bold")

    ax.set_xlim(-0.6, 5.2)
    ax.set_ylim(22, -3)
    ax.set_xlabel("Gaia colour  BP − RP   (redder, cooler →)")
    ax.set_ylabel("absolute magnitude  M$_G$   (brighter ↑)")
    ax.grid(False)
    ax.set_title("How every star is classified")
    ps.subtitle(ax, "Dashed lines are the class boundaries; each sits in a gap between populations.")
    ax.legend(loc="lower left", markerscale=12, handletextpad=0.3)
    ps.caption(fig, f"{PROVENANCE} Gaia stars only (log density). Main-sequence dwarfs take the spectral letter of "
               "their temperature, itself read from colour via Mamajek's dwarf sequence. Giants must lie 1.5 mag "
               "above the main sequence and be brighter than M_G = 4.5; white dwarfs must fall below the lower "
               "diagonal and be bluer than BP-RP = 2; anything fainter than M_G = 16.5 is an L/T dwarf.")
    return ps.save(fig, FIGURE_DIR / "05_classification")


# --------------------------------------------------------------------------
# 6. Gaia against the 10 parsec census
# --------------------------------------------------------------------------

def figure_ten_parsecs(res):
    """Which of the known objects within 10 pc our catalogue contains, and why not."""
    v = res["validation"]["census_10pc"]["by_kind"]
    kinds = [("star", "Stars"), ("white dwarf", "White dwarfs"), ("brown dwarf", "Brown dwarfs")]
    segments = [
        ("gaia", "measured by Gaia", ps.SERIES[0]),
        ("hipparcos", "from Hipparcos (too bright for Gaia)", ps.SERIES[1]),
        ("binary: companion not resolved by Gaia", "missed: companion not resolved", ps.SERIES[2]),
        ("binary: Gaia parallax spoiled by orbital motion", "missed: orbit spoiled the parallax", ps.SERIES[3]),
        ("too faint for Gaia", "missed: too faint", ps.SERIES[4]),
    ]
    fig, ax = plt.subplots(figsize=(11, 3.9))
    ypos = np.arange(len(kinds))[::-1]
    for yy, (k, label) in zip(ypos, kinds):
        left = 0
        total = v[k]["total"]
        for key, seg_label, colour in segments:
            n = v[k][key] if key in ("gaia", "hipparcos") else v[k]["missed_reasons"][key]
            if n == 0:
                continue
            width = 100 * n / total
            ax.barh(yy, width, left=left, height=0.56, color=colour, edgecolor=ps.SURFACE, linewidth=2)
            if width > 5:
                ax.text(left + width / 2, yy, f"{n}", ha="center", va="center", fontsize=8,
                        color=_ink_on(colour), fontweight="bold")
            left += width
        other = v[k]["missed_reasons"].get("other", 0) + v[k]["missed_reasons"].get("single: no usable Gaia parallax", 0)
        found = v[k]["gaia"] + v[k]["hipparcos"]
        ax.text(101.5, yy, f"{found} of {total} found" + (f"  ({other} other)" if other else ""),
                va="center", fontsize=8, color=ps.INK_SECONDARY)
    ax.set_yticks(ypos, [label for _, label in kinds])
    ax.set_xlim(0, 125)
    ax.set_xticks(np.arange(0, 101, 20), [f"{t}%" for t in range(0, 101, 20)])
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in segments]
    ax.legend(handles, [s for _, s, _ in segments], loc="upper center", bbox_to_anchor=(0.42, -0.14),
              ncol=3, handlelength=1.2)
    total_all = sum(v[k]["total"] for k, _ in kinds)
    found_all = sum(v[k]["gaia"] + v[k]["hipparcos"] for k, _ in kinds)
    ax.set_title(f"Our catalogue against every object known within 10 pc: {found_all} of {total_all}")
    ps.subtitle(ax, "Stars are missed only when they have a close companion. Brown dwarfs are mostly too faint.")
    ps.caption(fig, "Reference: the 10 parsec catalogue of Reyle et al. (2021), 2023 update, which combines Gaia "
               "with ground-based parallaxes, infrared surveys and adaptive-optics imaging. Exoplanets excluded. "
               "An object counts as found if its Gaia DR3 source (or Hipparcos number) is in our catalogue.", y=-0.2)
    return ps.save(fig, FIGURE_DIR / "06_ten_parsec_census")


# --------------------------------------------------------------------------
# 7. Are we above or below the disc?
# --------------------------------------------------------------------------

def figure_north_south(res):
    ns = res["north_south"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"wspace": 0.26})

    edges = np.array(ns["in_plane_ring_edges_pc"])
    centres = 0.5 * (edges[1:] + edges[:-1])
    dens = np.array(ns["in_plane_density_per_pc3"])
    ax1.plot(centres, dens / dens[0], color=ps.SERIES[0], marker="o", ms=5, mec=ps.SURFACE, mew=1.5)
    ax1.axhline(1.0, color=ps.AXIS, lw=0.8)
    ax1.set_xlabel("distance within the Galactic plane (pc)")
    ax1.set_ylabel("density, relative to the 10-20 pc ring")
    ax1.set_title("Why a direct fit fails")
    ps.subtitle(ax1, f"Even within the plane, density falls {ns['in_plane_density_drop_percent']:.0f}% with distance.")
    ax1.set_ylim(0.9, 1.03)

    z = np.array(ns["bin_centres_pc"])
    ratio, err = np.array(ns["ratio_south_north"]), np.array(ns["ratio_err"])
    preds = ns["prediction"]["ratio"]
    lo, mid, hi = (np.array(preds[k]) for k in sorted(preds, key=float))
    ax2.fill_between(z, np.minimum(lo, hi), np.maximum(lo, hi), color=ps.MUTED, alpha=0.25, lw=0)
    ax2.plot(z, mid, color=ps.MUTED, lw=1.2)
    ps.annotate(ax2, "prediction: Sun 20.8 pc above an\nisothermal disc, scale height 250-350 pc",
                (40, float(np.interp(40, z, mid))), (14, 0.968), color=ps.INK_SECONDARY)
    ax2.axhline(1.0, color=ps.AXIS, lw=0.8)
    ax2.text(2, 0.994, "no asymmetry", color=ps.MUTED, fontsize=7.5, va="top")
    ax2.errorbar(z, ratio, yerr=err, color=ps.SERIES[0], marker="o", ms=5, mec=ps.SURFACE, mew=1.5,
                 capsize=0, elinewidth=1, lw=2)
    ax2.set_xlabel("height |z| above or below the Sun (pc), 10 pc slabs")
    ax2.set_ylabel("stars below us / stars above us")
    ax2.set_title("More stars below us than above")
    ps.subtitle(ax2, f"{100 * (ns['pooled_ratio_10_90pc'] - 1):.1f} ± {100 * ns['pooled_ratio_err']:.1f}% "
                     f"excess for |z| = 10-90 pc ({ns['pooled_significance_sigma']:.1f}σ)")
    ax2.set_ylim(0.95, 1.18)
    ps.caption(fig, f"{PROVENANCE} Dwarfs brighter than M_G = 13 (complete far beyond 100 pc); the two clusters "
               "are masked and their volume removed. A slab above and a slab below have identical volume and "
               "distance distribution, so the distance-dependent loss on the left cancels in the ratio on the "
               "right. The excess has the sign and rough size expected for a Sun above the mid-plane, but not "
               "the predicted shape, so this 100 pc bubble cannot measure the Sun's height itself.")
    return ps.save(fig, FIGURE_DIR / "07_north_south")


# --------------------------------------------------------------------------
# 8. Famous stars: a distance ladder
# --------------------------------------------------------------------------

def figure_famous(res):
    rows = [r for r in res["famous"]["inside"] + res["famous"]["beyond"] if r["dist_pc"]]
    rows.sort(key=lambda r: r["dist_pc"])
    fig, ax = plt.subplots(figsize=(10, 9.6))
    ax.axvspan(np.log10(MAX_DISTANCE_PC), np.log10(600), color=ps.MUTED, alpha=0.10, lw=0)
    ax.axvline(np.log10(MAX_DISTANCE_PC), color=ps.INK_SECONDARY, lw=0.9)
    ax.text(np.log10(MAX_DISTANCE_PC) + 0.02, len(rows) - 0.3, "edge of the map, 100 pc", color=ps.INK_SECONDARY,
            fontsize=8, va="top")
    markers = {"Gaia DR3": "o", "Hipparcos": "D"}
    for i, r in enumerate(rows):
        yy = len(rows) - 1 - i
        teff = r["teff"] if r.get("teff") else None
        colour = COLOURS(np.atleast_1d(teff))[0] if teff else ps.MUTED
        src = r["source"] if r["source"] in markers else "SIMBAD"
        mk = markers.get(src, "o")
        filled = src != "SIMBAD"
        ax.plot(np.log10(r["dist_pc"]), yy, marker=mk, ms=8 if mk == "o" else 6.5,
                mfc=colour if filled else "none", mec=colour if filled else ps.INK_SECONDARY, mew=1.2, zorder=4)
        ax.text(np.log10(r["dist_pc"]) + 0.035, yy, f"{r['label']}  {r['dist_pc']:.1f} pc", va="center",
                fontsize=8, color=ps.INK)
    ax.set_yticks([])
    ax.set_ylim(-1, len(rows))
    _log_ticks(ax, [1, 2, 5, 10, 20, 50, 100, 200, 500])
    ax.set_xlim(0.0, np.log10(600))
    ax.set_xlabel("distance (pc)   —   1 pc = 3.26 light years")
    ax.grid(axis="y", visible=False)
    for spine in ("left",):
        ax.spines[spine].set_visible(False)
    handles = [Line2D([], [], marker="o", ls="none", mfc=ps.INK_SECONDARY, mec=ps.INK_SECONDARY, label="Gaia DR3"),
               Line2D([], [], marker="D", ls="none", mfc=ps.INK_SECONDARY, mec=ps.INK_SECONDARY, ms=6,
                      label="Hipparcos (too bright for Gaia)"),
               Line2D([], [], marker="o", ls="none", mfc="none", mec=ps.INK_SECONDARY,
                      label="other parallax, via SIMBAD (not in our catalogue)")]
    ax.legend(handles=handles, loc="lower left")
    ax.set_title("Where the famous stars are")
    ps.subtitle(ax, "The night sky's brightest stars are mostly rare, highly luminous ones -- giants and hot blue "
                    "stars -- seen from far away, so many lie beyond the map.")
    ps.caption(fig, "Marker colour: each star's true colour where our catalogue has its temperature. Positions "
               "beyond 100 pc are SIMBAD's adopted parallaxes, mostly Hipparcos; for supergiants such as "
               "Betelgeuse and Deneb these are uncertain by tens of percent.", y=0.03)
    return ps.save(fig, FIGURE_DIR / "08_famous_stars")


# --------------------------------------------------------------------------
# 9. Finding the clusters
# --------------------------------------------------------------------------

def figure_clusters(cat, res):
    st = res["structure"]
    h = st["histogram"]
    n = np.array(h["n_neighbours"])
    real, null = np.array(h["real"], dtype=float), np.array(h["null_mean_of_3"], dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"wspace": 0.22, "width_ratios": [1.1, 1]})

    ax1.step(n, np.where(null > 0, null, np.nan), where="mid", color=ps.MUTED, lw=1.4,
             label="same stars, random directions (null)")
    ax1.step(n, np.where(real > 0, real, np.nan), where="mid", color=ps.SERIES[0], lw=1.6, label="real sky")
    ax1.axvline(max(st["null_max_neighbours"]), color=ps.MUTED, lw=0.8)
    ax1.text(max(st["null_max_neighbours"]) + 1, 3e4, "the best clump\nchance can make", color=ps.INK_SECONDARY,
             fontsize=7.5, va="top")
    for p in st["peaks"]:
        name = p["match"]["name"] if p.get("match") else "unidentified"
        ps.annotate(ax1, f"{name}\n{p['peak_neighbours']} neighbours", (p["peak_neighbours"], 1.2),
                    (p["peak_neighbours"] - 18, 40), color=ps.INK)
    ax1.set_yscale("log")
    ax1.set_xlabel(f"neighbours within {st['search_radius_pc']:g} pc of a star")
    ax1.set_ylabel("number of stars")
    ax1.set_title("Is the clump real?")
    ps.subtitle(ax1, "Randomising every star's direction, at its true distance, destroys real clustering.")
    ax1.legend(loc="upper right")

    hy = st["peaks"][0]
    cx, cy, cz = hy["centre_xyz_pc"]
    x, y, z = (np.asarray(cat[k], dtype=float) for k in ("x", "y", "z"))
    slab = np.abs(z - cz) < 12
    size = star_sizes(np.asarray(cat["m_g"], dtype=float)[slab], scale=1.4)
    rgb = star_colours(cat)[slab]
    order = np.argsort(size)
    ax2.scatter(x[slab][order], y[slab][order], s=size[order], c=rgb[order], alpha=0.85, linewidths=0, rasterized=True)
    ax2.add_patch(Circle((cx, cy), 10, fill=False, color=ps.INK_SECONDARY, lw=0.9))
    ax2.set_xlim(cx - 32, cx + 32)
    ax2.set_ylim(cy - 32, cy + 32)
    ax2.set_aspect("equal")
    ax2.grid(False)
    ax2.set_xlabel("X  (pc)")
    ax2.set_ylabel("Y  (pc)")
    ax2.set_title("The Hyades, 47 pc away")
    m = hy["match"]
    ps.subtitle(ax2, f"A 24 pc-thick slice through the cluster. Found at {hy['distance_pc']:.1f} pc; "
                     f"literature {m['literature_distance_pc']:.1f} pc.")
    ps.caption(fig, f"{PROVENANCE} Left: for every star, the number of other stars within {st['search_radius_pc']:g} pc. "
               "The null keeps each star's distance (and so all of Gaia's distance-dependent selection) and throws "
               "away its direction. Circle: the 10 pc radius used to count members and mask the cluster.")
    return ps.save(fig, FIGURE_DIR / "09_finding_the_clusters")


def make_all(res, art) -> list:
    cat, seq = art["catalogue"], art["sequence"]
    written = []
    for fn in (lambda: figure_projections(cat, res), lambda: figure_census(res),
               lambda: figure_selection(cat, res, seq), lambda: figure_density(res),
               lambda: figure_hr(cat, res, seq), lambda: figure_ten_parsecs(res),
               lambda: figure_north_south(res), lambda: figure_famous(res),
               lambda: figure_clusters(cat, res)):
        written += fn()
    return written
