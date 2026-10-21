"""Every static figure, from the measurements in ``run.build``.

00 sample       the fair census: flux limit, volume limit, mass completeness
01 tuning fork  the poster (``poster.py``), plus a gallery of all candidates
02 bias         distant galaxies look smoother, and how well it is corrected
03 colours      red spirals and blue ellipticals
04 bars         how common bars are, and in which galaxies
05 arms         how many spiral arms, and in which galaxies
06 environment  morphology-density versus colour-density, at fixed mass
07 experts      volunteers against a professional classifier

Colours follow ``plotstyle.py``: categorical slots in fixed order; galaxy
colours drawn in the colour they are named for; one-hue ramps for ordered
quantities (stellar mass bins), validated for colour-vision deficiency.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import fetch_data as fd
import literature as lit
import morphology as mo
import physics as ph
import poster
from analysis import faint_limit, is_red
from config import COMPLETENESS_R, FIGURE_DIR, GZ2_R_LIMIT
from plotstyle import (
    AQUA, AXIS, BLUE, BLUE_RAMP, C_BLUE_GALAXY, C_OURS, C_PUBLISHED, C_RED_GALAXY, GRID, INK,
    INK_SECONDARY, MUTED, ORANGE, SURFACE, VIOLET, annotate, caption, save, use_style,
)

use_style()

#: The three readings of the votes: raw, Willett et al. (2013), Hart et al. (2016).
READER_COLOURS = {"hart16": BLUE, "willett13": ORANGE, "raw": AQUA}
READER_LABELS = {"hart16": "Hart et al. 2016 (used here)", "willett13": "Willett et al. 2013",
                 "raw": "raw votes"}
#: Ordinal ramp for the stellar-mass bins (validated: monotone, one hue).
MASS_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
DENSITY_CMAP = LinearSegmentedColormap.from_list("blues", BLUE_RAMP)
GREY_CMAP = LinearSegmentedColormap.from_list("greys", ["#fcfcfb", "#c3c2b7", "#52514e", "#0b0b0b"])


def _out(name: str):
    return FIGURE_DIR / name


def _img(path):
    from PIL import Image

    return np.asarray(Image.open(path).convert("RGB"))


def _tile(ax, objid, title=None, sub=None, crop: float = 0.2, fontsize: float = 7.5) -> None:
    """One galaxy image (the GZ2 subject image, central part) in an axes."""
    path = fd.gz2_image_path(objid)
    if path.exists():
        img = _img(path)
        n = img.shape[0]
        k = int(round(n * crop))
        ax.imshow(img[k:n - k, k:n - k], interpolation="lanczos")
    else:  # not downloaded (offline run): say so rather than fail
        ax.set_facecolor("#1a1a19")
        ax.text(0.5, 0.5, "image not\ndownloaded", ha="center", va="center", color=MUTED,
                fontsize=7, transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.grid(False)
    if title:
        ax.set_title(title, fontsize=fontsize + 0.5, fontweight="bold", loc="center", pad=3,
                     color=INK)
    if sub:
        ax.text(0.5, -0.04, sub, transform=ax.transAxes, ha="center", va="top",
                fontsize=fontsize, color=INK_SECONDARY, linespacing=1.25)


def needed_images(state: dict) -> list[int]:
    """Every galaxy any figure shows, so they can be fetched in one go."""
    return sorted({*[e["dr7objid"] for v in poster.select(state["c"], state["d"]).values()
                     for e in v],
                   *bias_examples(state), *colour_examples(state)["red_spirals"],
                   *colour_examples(state)["blue_ellipticals"],
                   *[o for o in arm_examples(state) if o is not None]})


def make_all(state: dict, *, offline: bool = False, quiet: bool = False) -> list:
    c = state["c"]
    ids = needed_images(state)
    if not offline:
        fd.fetch_gz2_images(ids, ra=c["ra"], dec=c["dec"], r90=c["petroR90_r"], quiet=quiet)
    paths = []
    for fn in (fig_sample, fig_poster, fig_bias, fig_colours, fig_bars, fig_arms,
               fig_environment, fig_experts):
        if not quiet:
            print(f"  figure: {fn.__name__}", flush=True)
        paths += fn(state)
    return paths


# ---------------------------------------------------------------------------
# 00 The sample
# ---------------------------------------------------------------------------

def fig_sample(state: dict) -> list:
    d, vl, kz, mi = state["d"], state["vl"], state["kz"], state["mass_info"]
    res = state["results"]
    g = d["good"] & (d["z"] < 0.2)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={"wspace": 0.3})
    fig.subplots_adjust(right=0.8)

    h = a.hexbin(d["z"][g], d["M"][g], gridsize=(120, 70), bins="log", mincnt=1,
                 cmap=GREY_CMAP, extent=(0, 0.2, -24.5, -16), rasterized=True, linewidths=0)
    zs = np.linspace(0.005, 0.17, 400)
    a.plot(zs, faint_limit(zs, kz), color=ORANGE, lw=1.6)
    annotate(a, f"flux limit r = {COMPLETENESS_R}\nfor the reddest galaxies",
             (0.15, float(faint_limit(0.15, kz))), (0.12, -17.2), color=INK_SECONDARY)
    z0, z1 = vl["z_range"]
    a.add_patch(plt.Rectangle((z0, -24.5), z1 - z0, vl["M_lim"] + 24.5, fill=False, ec=BLUE,
                              lw=1.6))
    a.text(z0 + 0.002, -24.3, f"volume-limited sample\n{vl['n']:,d} galaxies, "
                              f"$M_r < {vl['M_lim']:.2f}$", color=BLUE, fontsize=8, va="top")
    a.set_xlim(0, 0.2)
    a.set_ylim(-16, -24.5)
    a.set_xlabel("redshift z")
    a.set_ylabel("absolute magnitude $M_r$ (brighter up)")
    a.set_title("A. Far away, only luminous galaxies make the cut")
    cb = fig.colorbar(h, ax=a, pad=0.01, fraction=0.04)
    cb.set_label("galaxies per cell", color=INK_SECONDARY)
    cb.outline.set_visible(False)

    gm = d["good"] & np.isfinite(d["mass"]) & (d["z"] < 0.12)
    b.hexbin(d["z"][gm], d["mass"][gm], gridsize=(100, 60), bins="log", mincnt=1, cmap=GREY_CMAP,
             extent=(0, 0.12, 9.0, 12.0), rasterized=True, linewidths=0)
    from matplotlib.transforms import blended_transform_factory

    right = blended_transform_factory(b.transAxes, b.transData)
    for k, bb in enumerate([x for x in mi if x["usable"]]):
        lo_z, hi_z = bb["z_range"]
        b.add_patch(plt.Rectangle((lo_z, bb["mass_lo"]), hi_z - lo_z, bb["mass_hi"] - bb["mass_lo"],
                                  fc=MASS_RAMP[k], ec=MASS_RAMP[k], lw=1.6, alpha=0.30))
        b.add_patch(plt.Rectangle((lo_z, bb["mass_lo"]), hi_z - lo_z, bb["mass_hi"] - bb["mass_lo"],
                                  fill=False, ec=MASS_RAMP[k], lw=1.6))
        b.text(1.03, 0.5 * (bb["mass_lo"] + bb["mass_hi"]),
               f"log M* {bb['mass_lo']:.1f}-{bb['mass_hi']:.1f}: {bb['n']:,d} galaxies,\n"
               f"z {lo_z:.2f}-{hi_z:.3f}", transform=right, fontsize=7.5, va="center",
               color=INK_SECONDARY)
    b.set_xlim(0, 0.12)
    b.set_ylim(9.0, 12.0)
    b.set_xlabel("redshift z")
    b.set_ylabel("stellar mass  log$_{10}$($M_*$ / $M_\\odot$)")
    b.set_title("B. Each mass bin, out to where it is complete")
    caption(fig, f"Galaxy Zoo 2 holds every SDSS galaxy brighter than r = {GZ2_R_LIMIT} (DR7). "
                 "Fractions are only fair in samples complete for the galaxies they include. "
                 "A: the main volume-limited sample. B: for 'at fixed mass' results each "
                 "stellar-mass bin keeps only the redshifts out to which even its reddest (faintest "
                 "for their mass) 98% of galaxies are bright enough. Stellar masses: MPA-JHU (SDSS DR17).")
    return save(fig, _out("00_sample"))


# ---------------------------------------------------------------------------
# 01 Poster
# ---------------------------------------------------------------------------

def fig_poster(state: dict) -> list:
    c, d = state["c"], state["d"]
    sel = poster.select(c, d)
    paths = poster.draw(sel, fd.gz2_image_path, _out("01_tuning_fork"),
                        n_volunteers=lit.WILLETT13["n_volunteers"], n_galaxies=len(c["dr7objid"]))
    paths += poster.draw_gallery(sel, fd.gz2_image_path, _out("01b_tuning_fork_candidates"))
    state["results"]["poster"] = {k: [{kk: vv for kk, vv in e.items() if kk != "index"}
                                      for e in v] for k, v in sel.items()}
    return paths


# ---------------------------------------------------------------------------
# 02 Classification bias
# ---------------------------------------------------------------------------

def bias_examples(state: dict) -> list[int]:
    """One typical featured galaxy per redshift step, all of one luminosity and size.

    In each step, among galaxies of the panel-B luminosity and size range that
    Hart et al.'s debiasing calls featured (p > 0.9), the one whose raw vote
    is closest to the median raw vote of that group: a typical case, not the
    most dramatic one.
    """
    c, d = state["c"], state["d"]
    ok = (d["good"] & (d["M"] >= -21.8) & (d["M"] < -21.4) & (d["R50_kpc"] >= 3) & (d["R50_kpc"] < 6)
          & (np.nan_to_num(c["p_features"]) > 0.9))
    out = []
    for zc in (0.03, 0.045, 0.06, 0.075, 0.09):
        m = ok & (np.abs(d["z"] - zc) < 0.004)
        idx = np.flatnonzero(m)
        raw = c["wf_features"][idx]
        out.append(int(c["dr7objid"][idx[np.argmin(np.abs(raw - np.median(raw)))]]))
    return out


def fig_bias(state: dict) -> list:
    c, d = state["c"], state["d"]
    rb = state["results"]["redshift_bias"]
    fig = plt.figure(figsize=(12, 8.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.25, 1.0], hspace=0.72, wspace=0.32)

    # Row 1: what the volunteers saw.
    ex = bias_examples(state)
    sub = fig.add_gridspec(1, 5, left=0.08, right=0.92, top=0.955, bottom=0.735, wspace=0.12)
    lookup = {int(o): i for i, o in enumerate(c["dr7objid"])}
    for k, oid in enumerate(ex):
        ax = fig.add_subplot(sub[0, k])
        i = lookup[oid]
        _tile(ax, oid, f"z = {d['z'][i]:.3f}",
              f"{100 * c['wf_features'][i]:.0f}% saw features\n(debiased: {100 * c['p_features'][i]:.0f}%)")
    fig.text(0.08, 0.985, "A. The same kind of galaxy, farther and farther away "
                          "(M$_r$ -21.8 to -21.4, half-light radius 3-6 kpc)",
             fontsize=11, fontweight="bold", va="top")

    # Row 2, left: the mean vote for galaxies of one luminosity and size.
    ax = fig.add_subplot(gs[1, :2])
    one = rb["one_kind"]
    z = np.array([r["z"] for r in one["bins"]])
    for name in ("raw", "willett13", "hart16"):
        y = np.array([r[name] for r in one["bins"]])
        e = np.array([r[name + "_err"] for r in one["bins"]])
        ax.errorbar(z, y, e, color=READER_COLOURS[name], marker="o", ms=5, lw=1.8, capsize=0,
                    mec=SURFACE, mew=1.2)
        dy = {"hart16": 0.035, "willett13": -0.035, "raw": 0.0}[name]
        ax.text(z[-1] + 0.003, y[-1] + dy, READER_LABELS[name], color=INK_SECONDARY, fontsize=8,
                va="center")
    ax.set_xlim(0.015, 0.13)
    ax.set_ylim(0, 0.85)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("mean vote fraction for\n'features or disk'")
    ax.set_title("B. Raw votes lose disks with distance;\nthe corrections restore them", fontsize=10)

    # Row 2, right + row 3: the flatness test.
    fl = rb["flatness"]
    ze = np.asarray(fl["z_edges"])
    zc = 0.5 * (ze[:-1] + ze[1:])
    names = ["smooth (p > 0.8)", "features (p > 0.8)", "spiral (product > 0.5)", "barred, of disks"]
    titles = ["Clean ellipticals", "Clean disks", "Spirals", "Barred, of disks"]
    slots = [gs[1, 2], gs[1, 3], gs[2, 2], gs[2, 3]]
    for k, (cn, tt, sl) in enumerate(zip(names, titles, slots)):
        ax = fig.add_subplot(sl)
        for r in ("raw", "willett13", "hart16"):
            f = np.array([x["f"] for x in fl["classes"][r][cn]["bins"]])
            ax.plot(zc, f, color=READER_COLOURS[r], marker="o", ms=4, lw=1.6, mec=SURFACE, mew=1)
        ax.set_title(tt, fontsize=9.5)
        ax.set_ylim(0, 0.65)
        ax.set_xlim(0.03, 0.085)
        ax.tick_params(labelsize=7.5)
        if k >= 2:
            ax.set_xlabel("redshift z")
        if k in (0, 2):
            ax.set_ylabel("fraction of the\nvolume-limited sample", fontsize=8)
    fig.text(0.54, 0.672, "C. The flatness test: in a volume-limited sample "
                          "every fraction should be flat", fontsize=10, fontweight="bold")
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color=READER_COLOURS[r], marker="o", ms=4, lw=1.6, mec=SURFACE)
               for r in ("raw", "willett13", "hart16")]
    fig.legend(handles, [READER_LABELS[r] for r in ("raw", "willett13", "hart16")],
               loc="upper left", bbox_to_anchor=(0.54, 0.665), ncol=3, fontsize=8, frameon=False,
               handlelength=1.6, columnspacing=1.4)

    # Row 3, left: what the debiasing does.
    ax = fig.add_subplot(gs[2, :2])
    ax.axis("off")
    h16 = fl["classes"]["hart16"]
    raw = fl["classes"]["raw"]
    w13 = fl["classes"]["willett13"]
    lines = ["Change from z = 0.03-0.04 to 0.07-0.085, in percentage points:", ""]
    for cn, tt in zip(names, titles):
        lines.append(f"{tt:<20} raw {100 * raw[cn]['drift']:+5.1f}   "
                     f"Willett {100 * w13[cn]['drift']:+5.1f}   Hart {100 * h16[cn]['drift']:+5.1f}")
    ax.text(0.0, 0.95, "\n".join(lines), family="DejaVu Sans Mono", fontsize=8.2, va="top",
            color=INK_SECONDARY, linespacing=1.6, transform=ax.transAxes)
    caption(fig, "Volunteers classified 424 x 424 pixel SDSS images scaled to each galaxy's size, so a "
                 "distant galaxy is drawn with fewer, noisier pixels, and its arms and bars fade. The "
                 "galaxies in A and B all have the same luminosity and size, so their true mix of shapes "
                 "should not depend on distance. Hart et al. (2016) correct each answer by matching the "
                 "distribution of votes to that of similar galaxies at z = 0.03; Willett et al. (2013) "
                 "rescaled ratios of answers. Neither is perfect: Hart's leaves drifts of 4-9 points "
                 "(C), a systematic uncertainty on every fraction in this project.", y=0.035)
    return save(fig, _out("02_redshift_bias"))


# ---------------------------------------------------------------------------
# 03 Colours by shape
# ---------------------------------------------------------------------------

def colour_examples(state: dict) -> dict:
    """Four red spirals and four blue ellipticals, by rule.

    Red spirals: the Masters et al. (2010)-style sample (face-on, disk-
    dominated, p_spiral path >= 0.8, redder than their red-sequence line),
    the highest spiral scores. Blue ellipticals: the Schawinski et al.
    (2009)-style sample (smooth, 3 sigma bluer than the red sequence in u - r),
    the highest p_smooth.
    """
    c = state["c"]
    wc = state["results"]["wrong_colours"]
    m10 = wc["masters10_their_cut"]
    rs = np.flatnonzero(m10["_sel"] & m10["_red"])
    pp = mo.path_probabilities(c)["spiral"]
    rs = rs[np.argsort(-pp[rs])][:4]
    be = np.flatnonzero(wc["schawinski09"]["_blue"])
    be = be[np.argsort(-np.nan_to_num(c["p_smooth"][be]))][:4]
    return {"red_spirals": [int(c["dr7objid"][i]) for i in rs],
            "blue_ellipticals": [int(c["dr7objid"][i]) for i in be]}


def fig_colours(state: dict) -> list:
    c, d, cls, vl = state["c"], state["d"], state["cls"], state["vl"]
    div = state["divider"]
    wc = state["results"]["wrong_colours"]
    fig = plt.figure(figsize=(12.5, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.42, wspace=0.42)
    g = d["good"] & (d["z"] < 0.085) & (d["z"] > 0.01)
    mm = np.linspace(-24, -18, 50)
    line = div["slope"] * mm + div["intercept"]
    for k, (key, title, frac_key, word) in enumerate([
            ("smooth", "A. Clean ellipticals ('smooth', p > 0.8)", "smooth_blue", "blue"),
            ("features", "B. Clean disks ('features', p > 0.8)", "features_red", "red")]):
        ax = fig.add_subplot(gs[0, k])
        m = g & cls[key]
        ax.hexbin(d["M"][m], d["gr0"][m], gridsize=(70, 55), bins="log", mincnt=1, cmap=DENSITY_CMAP,
                  extent=(-24, -18, 0.1, 1.1), rasterized=True, linewidths=0)
        ax.plot(mm, line, color=INK, lw=1.2)
        ax.text(-18.15, div["slope"] * -18.15 + div["intercept"] - 0.035, "D2's red/blue divider",
                fontsize=7.5, color=INK, va="top")
        f = wc["vl"][frac_key]
        ax.text(0.03, 0.97, f"{100 * f['f']:.0f}% {word} among the bright\n"
                            f"(M$_r$ < {state['vl']['M_lim']:.1f}, n = {f['n']:,d})",
                transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color=INK)
        ax.set_xlim(-18, -24)
        ax.set_ylim(0.1, 1.1)
        ax.set_xlabel("$M_r$ (brighter right)")
        ax.set_ylabel("rest-frame g - r")
        ax.set_title(title, fontsize=10)

    # C: how many spirals are red? It depends on the question.
    ax = fig.add_subplot(gs[0, 2])
    v = wc["vl"]
    m10 = wc["masters10_their_cut"]
    items = [("edge-on disks", v["edgeon_red"], ""),
             ("all clean disks", v["features_red"], ""),
             ("spirals", v["spiral_red"], ""),
             ("face-on spirals", v["faceon_spiral_red"], ""),
             ("face-on, disky\n(Masters' cuts)", m10["red"], "   published: 6%")]
    y = np.arange(len(items))[::-1]
    for yy, (lab, f, note) in zip(y, items):
        ax.barh(yy, 100 * f["f"], height=0.55, color=C_RED_GALAXY, edgecolor="none")
        ax.text(100 * f["f"] + 1.5, yy, f"{100 * f['f']:.0f}%{note}", va="center", fontsize=8.5,
                color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items], fontsize=8)
    ax.set_xlim(0, 110)
    ax.set_xlabel("% that are red")
    ax.set_title("C. How many spirals are red?\nIt depends what you ask", fontsize=10)
    ax.grid(axis="y", visible=False)

    # D: pictures.
    ex = colour_examples(state)
    lookup = {int(o): i for i, o in enumerate(c["dr7objid"])}
    sub = fig.add_gridspec(2, 4, left=0.06, right=0.62, top=0.43, bottom=0.07, hspace=0.42,
                           wspace=0.08)
    for k, oid in enumerate(ex["red_spirals"]):
        ax = fig.add_subplot(sub[0, k])
        i = lookup[oid]
        _tile(ax, oid, None, f"g-r {d['gr0'][i]:.2f}, z {d['z'][i]:.3f}")
    for k, oid in enumerate(ex["blue_ellipticals"]):
        ax = fig.add_subplot(sub[1, k])
        i = lookup[oid]
        _tile(ax, oid, None, f"u-r {d['ur0'][i]:.2f}, z {d['z'][i]:.3f}")
    fig.text(0.06, 0.445, "D. Red spirals (top) and blue ellipticals (bottom)", fontsize=10,
             fontweight="bold")

    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    s09 = wc["schawinski09"]
    m10t = wc["masters10_their_z"]
    txt = (f"Red spirals, Masters et al. (2010) cuts:\n"
           f"  {100 * m10['red']['f']:.1f} ± {100 * m10['red']['err']:.1f}% of {m10['n']:,d}"
           f"   (published 6%)\n"
           f"  their redshift range: {100 * m10t['red']['f']:.1f}%\n"
           f"  star formation rate per unit mass:\n"
           f"  red {10 ** m10['ssfr_median_red'] * 1e11:.1f}, blue {10 ** m10['ssfr_median_blue'] * 1e11:.1f}"
           f" × 10⁻¹¹ per yr\n\n"
           f"Blue ellipticals, Schawinski et al. (2009):\n"
           f"  {100 * s09['blue']['f']:.1f} ± {100 * s09['blue']['err']:.1f}% of {s09['n_early']:,d}"
           f"   (published 5.7 ± 0.4%)")
    ax.text(0.0, 1.0, txt, va="top", fontsize=8.6, color=INK, linespacing=1.45,
            transform=ax.transAxes)
    caption(fig, "Colour tracks star formation: young stars are blue, so a galaxy that stopped forming "
                 "stars turns red. Most ellipticals are red and most disks blue, but not all. How many "
                 "disks count as red depends on the definition (C): dust reddens disks seen edge-on, and "
                 "big bulges are red whatever the disk does; Masters et al. kept only face-on, disk-"
                 "dominated spirals to find the truly passive ones. Their red spirals form stars at a "
                 "quarter of the rate of blue ones (MPA-JHU, spectra). Images: the GZ2 subject images, "
                 "chosen by rule (highest vote fractions in each sample).", y=0.0)
    return save(fig, _out("03_colours_by_shape"))


# ---------------------------------------------------------------------------
# 04 Bars
# ---------------------------------------------------------------------------

def fig_bars(state: dict) -> list:
    b = state["results"]["bars"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), gridspec_kw={"wspace": 0.42,
                                                                     "width_ratios": [1.1, 1, 1.15]})
    ax = axes[0]
    for key, col, lab in (("masters11", AQUA, "Masters et al. (2011)'s rule, raw votes"),
                          ("vl", BLUE, "debiased votes, bright disks")):
        rows = [r for r in b[key]["vs_colour"] if r["n"] >= 100]
        x = np.array([r["gr"] for r in rows])
        y = np.array([r["f"] for r in rows])
        e = np.array([r["err"] for r in rows])
        ax.errorbar(x, 100 * y, 100 * e, color=col, marker="o", ms=4.5, lw=1.6, mec=SURFACE, mew=1,
                    label=lab)
    ax.axhline(100 * lit.MASTERS11["bar_fraction"][0], color=C_PUBLISHED, lw=1.0,
               label="Masters et al. (2011): 29.4% overall")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_xlim(0.2, 1.0)
    ax.set_ylim(0, 75)
    ax.set_xlabel("rest-frame g - r")
    ax.set_ylabel("% of disks with a bar")
    ax.set_title("A. Bars are rare in the bluest disks", fontsize=10)

    ax = axes[1]
    rows = b["vs_mass"]
    x = np.array([r["mass"] for r in rows])
    for key, col, lab in (("red", C_RED_GALAXY, "red disks"), ("blue", C_BLUE_GALAXY, "blue disks")):
        y = np.array([r[key]["f"] for r in rows])
        e = np.array([r[key]["err"] for r in rows])
        ax.errorbar(x, 100 * y, 100 * e, color=col, marker="o", ms=5, lw=1.8, mec=SURFACE, mew=1)
        ax.text(x[0] - 0.05, 100 * y[0], lab, color=INK_SECONDARY, fontsize=8, va="center",
                ha="right")
    ax.set_xlim(9.6, 11.5)
    ax.set_ylim(0, 75)
    ax.set_xlabel("stellar mass log$_{10}$($M_*/M_\\odot$)")
    ax.set_ylabel("% of disks with a bar")
    ax.set_title("B. At fixed mass: red disks are barred\nmore often, except the most massive",
                 fontsize=10)

    ax = axes[2]
    m11 = {r["min_votes"]: r for r in b["masters11_n20"]}
    v = b["vl"]
    items = [("Masters' cuts, raw majority,\n≥10 bar votes", m11[10], False),
             ("same, ≥20 bar votes", m11[20], False),
             ("bright disks, raw majority", v["barred_raw"], False),
             ("bright disks, debiased\nmajority (used here)", v["barred"], True),
             ("bright disks, debiased\np_bar ≥ 0.3 (inclusive)", v["barred_p03"], False)]
    y = np.arange(len(items))[::-1]
    for yy, (lab, f, ours) in zip(y, items):
        ax.barh(yy, 100 * f["f"], height=0.55, color=BLUE if ours else MUTED, edgecolor="none")
        ax.text(100 * f["f"] + 1.2, yy, f"{100 * f['f']:.0f}%", va="center", fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items], fontsize=7.6)
    ax.set_xlim(0, 85)
    ax.set_xlabel("% of disks with a bar")
    ax.set_title("C. How common are bars?\nIt depends how the votes are read", fontsize=10)
    ax.grid(axis="y", visible=False)
    caption(fig, "A bar is a straight band of stars across a disk galaxy's centre; it funnels gas "
                 "inwards. Disks: debiased p_features x p_not-edge-on > 0.5, inclination below ~60 deg "
                 "(b/a > 0.5). A: bar fraction against colour, with Masters et al. (2011)'s rule on the "
                 "final GZ2 votes (their redshifts, luminosities and inclinations) and with debiased votes "
                 "in the bright volume-limited sample. B: at fixed stellar mass (mass-complete samples). "
                 "C: the answer moves by a factor of two with the reading; the volunteers' bar icon showed "
                 "a strong bar, and weak bars are mostly missed (Fig. 07).", y=-0.03)
    return save(fig, _out("04_bars"))


# ---------------------------------------------------------------------------
# 05 Spiral arms
# ---------------------------------------------------------------------------

def arm_examples(state: dict) -> list:
    """The face-on spiral with the most decisive vote for each arm number."""
    c, d = state["c"], state["d"]
    arm, pmax = mo.arm_number(c)
    base = (d["good"] & state["cls"]["spiral_clean"] & (d["z"] < 0.05) & (d["R90_arcsec"] > 15)
            & (c["expAB_r"] > 0.7))
    out = []
    for k in range(5):
        idx = np.flatnonzero(base & (arm == k))
        out.append(int(c["dr7objid"][idx[np.argmax(pmax[idx])]]) if idx.size else None)
    return out


def fig_arms(state: dict) -> list:
    a = state["results"]["arms"]
    c, d = state["c"], state["d"]
    fig = plt.figure(figsize=(12.5, 7.6))
    gs = fig.add_gridspec(1, 3, left=0.07, right=0.97, top=0.53, bottom=0.12, wspace=0.32)

    ex = arm_examples(state)
    sub = fig.add_gridspec(1, 5, left=0.08, right=0.92, top=0.96, bottom=0.66, wspace=0.12)
    lookup = {int(o): i for i, o in enumerate(c["dr7objid"])}
    arm, pmax = mo.arm_number(c)
    for k, oid in enumerate(ex):
        if oid is None:
            continue
        ax = fig.add_subplot(sub[0, k])
        i = lookup[oid]
        _tile(ax, oid, f"{mo.ARM_LABELS[k]} arm" + ("" if k == 0 else "s"),
              f"{100 * pmax[i]:.0f}% of arm votes, z {d['z'][i]:.3f}")

    ax = fig.add_subplot(gs[0, 0])
    labs = mo.ARM_LABELS
    x = np.arange(5)
    for off, key, pub, col, lab in ((-0.17, "arm_fractions_lum", lit.HART16["arm_fractions_lum"], C_OURS,
                                     "this project"),
                                    (0.17, None, lit.HART16["arm_fractions_lum"], C_PUBLISHED,
                                     "Hart et al. 2016")):
        vals = [a["arm_fractions_lum"][l] for l in labs] if key else [pub[l] for l in labs]
        ax.bar(x + off, 100 * np.array(vals), width=0.32, color=col, edgecolor="none", label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_xlabel("number of spiral arms")
    ax.set_ylabel("% of spirals")
    ax.legend(loc="upper right")
    ax.set_title(f"A. Two arms is the norm\n({a['n_spirals_lum']:,d} spirals; Hart: 17,957)",
                 fontsize=10)
    ax.grid(axis="x", visible=False)

    ax = fig.add_subplot(gs[0, 1])
    rows = a["vs_mass"]
    xm = np.array([r["mass"] for r in rows])
    for lab, col, key in (("2 arms", BLUE, "2"), ("5+ arms", ORANGE, "5+")):
        y = np.array([r["f" + key] for r in rows])
        e = np.array([r["e" + key] for r in rows])
        ax.errorbar(xm, 100 * y, 100 * e, color=col, marker="o", ms=3.5, lw=1.4, mec=SURFACE, mew=0.8)
        ax.text(xm[0] - 0.02, 100 * y[0], lab, color=INK_SECONDARY, fontsize=8, ha="right", va="center")
    for val, err, mass in ((lit.HART16["m2_fraction_highest_mass"], None, 11.2),
                           (lit.HART16["m5_fraction_highest_mass"], None, 11.2)):
        ax.errorbar([mass], [100 * val[0]], [100 * val[1]], color=C_PUBLISHED, marker="D", ms=5,
                    lw=1.2, mec=SURFACE)
    ax.text(11.22, 100 * lit.HART16["m2_fraction_highest_mass"][0] - 6, "Hart et al.,\nhighest mass",
            color=C_PUBLISHED, fontsize=7.5, ha="center", va="top")
    ax.set_xlim(10.45, 11.45)
    ax.set_ylim(0, 85)
    ax.set_xlabel("stellar mass log$_{10}$($M_*/M_\\odot$)")
    ax.set_ylabel("% of spirals")
    ax.set_title("B. The most massive spirals have\nmore arms", fontsize=10)

    ax = fig.add_subplot(gs[0, 2])
    rows = a["vs_colour"]
    xg = np.array([r["gi"] for r in rows])
    y = np.array([r["f2"] for r in rows])
    e = np.array([r["e2"] for r in rows])
    ax.errorbar(xg, 100 * y, 100 * e, color=BLUE, marker="o", ms=3.5, lw=1.4, mec=SURFACE, mew=0.8)
    ax.errorbar([0.67, 1.17], [32, 84], [2, 2], xerr=[0.07, 0.05], color=C_PUBLISHED, marker="D",
                ms=5, lw=1.2, ls="none", mec=SURFACE)
    ax.text(0.69, 27, "Hart et al.\n(bluest, reddest bins)", color=C_PUBLISHED, fontsize=7.5, va="top")
    mg = a["mass_matched_gi"]
    ax.text(0.98, 0.04, f"many-armed spirals are {abs(mg['difference']):.2f} mag\nbluer in g-i than "
                        f"two-armed ones\nof the same mass", transform=ax.transAxes, ha="right",
            fontsize=7.5, color=INK_SECONDARY)
    ax.set_xlim(0.55, 1.35)
    ax.set_ylim(0, 100)
    ax.set_xlabel("rest-frame g - i")
    ax.set_ylabel("% of spirals with 2 arms")
    ax.set_title("C. Blue spirals have more arms", fontsize=10)
    caption(fig, "Spirals of Hart et al. (2016)'s samples (0.03 < z < 0.085, M_r < -21 for H0 = 70; "
                 "B and C also log M* > 10.6), arm number = the answer with the highest debiased vote "
                 "fraction. Two-armed 'grand design' spirals dominate; many-armed ones are bluer at the "
                 "same mass, so their arms go with recent star formation. Our colours are model-magnitude "
                 "colours, ~0.1 mag redder than Hart's, which shifts panel C sideways.", y=0.005)
    return save(fig, _out("05_spiral_arms"))


# ---------------------------------------------------------------------------
# 06 Environment
# ---------------------------------------------------------------------------

def fig_environment(state: dict) -> list:
    env = state["results"]["environment"]
    bins = env["bins"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6), gridspec_kw={"hspace": 0.48, "wspace": 0.3})
    for ax, key, ylab, title in ((axes[0, 0], "f_red", "fraction red",
                                  "A. Colour depends strongly on environment..."),
                                 (axes[0, 1], "w_smooth", "early-type (smooth) likelihood",
                                  "B. ...shape much less")):
        for k, b in enumerate(bins):
            x = [cell["log_sigma_median"] for cell in b["cells"]]
            y = [cell[key] for cell in b["cells"]]
            ax.plot(x, y, color=MASS_RAMP[k], marker="o", ms=5, lw=1.8, mec=SURFACE, mew=1,
                    label=f"{b['mass_range'][0]:.1f}-{b['mass_range'][1]:.1f}")
        ax.set_xlabel("local density log$_{10}$ Σ$_5$ (h$^2$ Mpc$^{-2}$)")
        ax.set_ylabel(ylab)
        ax.set_xlim(-1.3, 0.9)
        # The same vertical span (0.6) in A and B, so the slopes compare directly.
        ax.set_ylim((0.1, 0.7) if key == "w_smooth" else (0.4, 1.0))
        ax.set_title(title, fontsize=10)
    axes[0, 0].legend(title="stellar mass, log M*", loc="lower right", fontsize=7.5,
                      title_fontsize=7.5)

    ax = axes[1, 0]
    y = np.arange(len(bins))[::-1]
    for yy, b in zip(y, bins):
        dd, de = b["decomposition"], b["decomposition_err"]
        a1, a2 = 100 * dd["red_from_morphology_mix"], 100 * dd["red_at_fixed_morphology"]
        ax.barh(yy, a2, height=0.55, color=C_RED_GALAXY, edgecolor="none")
        ax.barh(yy, a1, left=a2 + (0.6 if a1 > 0 else 0), height=0.55, color=MUTED, edgecolor="none")
        ax.text(a1 + a2 + 1.5, yy, f"{100 * dd['delta_red']:+.0f} pts", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"log M* {b['mass_range'][0]:.1f}-{b['mass_range'][1]:.1f}" for b in bins],
                       fontsize=8)
    ax.set_xlabel("rise in the red fraction, sparsest to densest fifth (points)")
    ax.set_xlim(0, 32)
    ax.grid(axis="y", visible=False)
    ax.set_title("C. Why dense places are redder", fontsize=10)
    ax.text(0.98, 0.06, "■ same shapes, redder\n", color=C_RED_GALAXY, transform=ax.transAxes,
            ha="right", fontsize=8)
    ax.text(0.98, 0.02, "■ more ellipticals", color=MUTED, transform=ax.transAxes, ha="right",
            fontsize=8)

    ax = axes[1, 1]
    for k, b in enumerate(bins):
        x = [cell["log_sigma_median"] for cell in b["cells"]]
        yv = [cell["red_given_disk"] for cell in b["cells"]]
        ax.plot(x, yv, color=MASS_RAMP[k], marker="o", ms=5, lw=1.8, mec=SURFACE, mew=1)
    ax.set_xlabel("local density log$_{10}$ Σ$_5$ (h$^2$ Mpc$^{-2}$)")
    ax.set_ylabel("fraction of disk galaxies that are red")
    ax.set_xlim(-1.3, 1.05)
    ax.set_ylim(0.3, 0.9)
    ax.set_title("D. Disks redden in crowds, small ones most\n(lines: mass bins as in A)",
                 fontsize=10)
    caption(fig, "GZ2 galaxies in mass bins complete over their redshift ranges (Fig. 00B); density as "
                 "in project D2 (5th-nearest neighbour within ±1000 km/s among D2's 134,975 volume-"
                 "limited tracers), galaxies whose density circle crosses a survey edge dropped. Shape is "
                 "counted fractionally: each galaxy's p_smooth / (p_smooth + p_features). C splits the rise "
                 "in the red fraction exactly into 'more ellipticals, which are redder' (grey) and 'each "
                 "shape redder' (red): 75-90% is the second, as Bamford et al. (2009) found. A and B "
                 "share one vertical span, so their slopes compare directly. Red is D2's divider, which "
                 "counts green-valley galaxies as red.", y=0.0)
    return save(fig, _out("06_environment"))


# ---------------------------------------------------------------------------
# 07 Volunteers versus experts
# ---------------------------------------------------------------------------

def fig_experts(state: dict) -> list:
    c = state["c"]
    ex = state["results"]["experts"]
    m = ex["_match"]
    gi, tt, bar, obl = m["gz2_index"], m["tt"], m["bar"], m["oblique"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), gridspec_kw={"wspace": 0.34,
                                                                   "width_ratios": [1.2, 1, 1]})
    ax = axes[0]
    groups = [("E", tt == -5), ("S0", (tt >= -3) & (tt <= 0)), ("Sa-Sb", (tt >= 1) & (tt <= 3)),
              ("Sbc-Sd", (tt >= 4) & (tt <= 7)), ("Sdm-Im", (tt >= 8) & (tt <= 10))]
    pf = np.asarray(c["wf_features"][gi], float)
    parts = ax.violinplot([pf[g] for _, g in groups], positions=np.arange(5), widths=0.8,
                          showmedians=True, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(BLUE)
        body.set_edgecolor("none")
        body.set_alpha(0.55)
    parts["cmedians"].set_color(INK)
    ax.set_xticks(np.arange(5))
    ax.set_xticklabels([f"{n}\n({g.sum():,d})" for n, g in groups], fontsize=8)
    ax.set_ylabel("fraction of volunteers who saw\n'features or a disk' (raw)")
    ax.set_xlabel("expert's Hubble type (Nair & Abraham 2010)")
    ax.set_title("A. Volunteers sort galaxies\nlike the expert", fontsize=10)
    ax.grid(axis="x", visible=False)

    ax = axes[1]
    pb = np.nan_to_num(np.asarray(c["w13_bar"][gi], float))
    edges = np.linspace(0, 1, 21)
    for lab, sel, col in (("no bar", obl & (bar == 0), MUTED), ("weak", obl & ((bar & 8) > 0), AQUA),
                          ("intermediate", obl & ((bar & 4) > 0), ORANGE),
                          ("strong", obl & ((bar & 2) > 0), BLUE)):
        h, _ = np.histogram(pb[sel], edges)
        cum = 1 - np.concatenate([[0], np.cumsum(h)]) / max(h.sum(), 1)
        ax.plot(edges, 100 * cum, color=col, lw=1.8)
        k = 9
        ax.text(edges[k] + 0.02, 100 * cum[k] + 2, f"{lab} ({sel.sum():,d})", fontsize=7.5,
                color=INK_SECONDARY)
    ax.axvline(0.5, color=AXIS, lw=0.8)
    ax.set_xlabel("GZ2 bar vote fraction (Willett et al. debiased)")
    ax.set_ylabel("% of the expert's class above it")
    ax.set_title("B. Strong bars are seen,\nweak ones mostly missed", fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 102)

    ax = axes[2]
    ax.axis("off")
    raw, w13 = ex["raw"], ex["willett13"]
    bb = w13["bars"]
    L = lit.WILLETT13
    rows = [("matched galaxies", f"{ex['n_matched']:,d}", f"{L['na10_overlap']:,d}"),
            ("oblique disks", f"{ex['n_oblique']:,d}", f"{L['na10_oblique_overlap']:,d}"),
            ("E: median p_smooth", f"{raw['E_median_psmooth']:.2f}±{raw['E_psmooth_sd']:.2f}",
             "0.86±0.07"),
            ("E: highest p_features", f"{raw['E_max_pfeatures']:.3f}", "0.741"),
            ("Sa+: median p_features", f"{raw['late_median_pfeatures']:.2f}", "0.88"),
            ("p_smooth>0.8: late types", f"{100 * raw['smooth08_late']:.1f}%", "4%"),
            ("  ...late types or S0", f"{100 * raw['smooth08_late_or_s0']:.0f}%", "20%"),
            ("barred, p_bar < 0.5", f"{100 * bb['barred_below_05']:.0f}%", "39%"),
            ("unbarred, p_bar > 0.5", f"{100 * bb['unbarred_above_05']:.1f}%", "6%"),
            ("strong bars, p > 0.8", f"{100 * bb['strong_above_08']:.0f}%", "78%"),
            ("intermediate, p > 0.8", f"{100 * bb['intermediate_above_08']:.0f}%", "40%"),
            ("weak, p > 0.8 / ≥ 0.5", f"{100 * bb['weak_above_08']:.0f} / {100 * bb['weak_above_05']:.0f}%",
             "9 / 32%")]
    ax.text(0.0, 1.0, f"{'':<26}{'here':>12}{'Willett+13':>12}", family="DejaVu Sans Mono",
            fontsize=8, fontweight="bold", va="top", transform=ax.transAxes)
    for k, (a_, b_, c_) in enumerate(rows):
        ax.text(0.0, 0.92 - k * 0.072, f"{a_:<26}{b_:>12}{c_:>12}", family="DejaVu Sans Mono",
                fontsize=8, va="top", transform=ax.transAxes, color=INK)
    ax.set_title("C. Willett et al.'s checks,\nredone", fontsize=10)
    caption(fig, "Nair & Abraham (2010): one expert classified 14,034 SDSS galaxies by eye, in g-band "
                 "images. Matched to GZ2 within 2\". A: raw vote fractions by expert type; S0s split "
                 "between the volunteers' 'smooth' and 'features'. B: among oblique disks (Willett et al. "
                 "2013's Table 3 cuts), the share of each expert bar class above a given volunteer bar "
                 "fraction. C: the published comparison reproduced from the public tables (Willett et al. "
                 "used raw votes for types, their debiased votes for bars); one number is not reproduced "
                 "(late types or S0 in the clean smooth sample; see README).", y=-0.03)
    return save(fig, _out("07_volunteers_vs_experts"))
