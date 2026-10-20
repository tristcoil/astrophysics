"""All figures, drawn from the analysis caches (no fitting happens here).

    python figures.py          # every figure into figures/ (PNG + SVG)
    python figures.py 04 06    # just some of them
"""

from __future__ import annotations

import json
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter

import config as cfg
import fetch_data as fd
import friedmann as fr
import likelihoods as lk
import plotstyle as ps

ps.use_style()
FIG = cfg.FIGURE_DIR


def _load(name: str):
    with np.load(cfg.DERIVED_DIR / f"{name}.npz", allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    return arrays, json.loads((cfg.DERIVED_DIR / f"{name}.json").read_text())


def _results() -> dict:
    return json.loads(cfg.RESULTS_PATH.read_text())


def _fate_key(f: dict) -> str:
    if "a_min" in f or "static" in f.get("past", "") or "no beginning" in f.get("past", ""):
        return "no big bang" if "a_min" in f else "accelerating"
    if f["future"] in ("Big Crunch", "Big Rip"):
        return "ends"
    return "accelerating" if "accelerating" in f["future"] else "slowing"


def _hpd(p):
    t68, t95, t997 = lk.credible_levels(p)
    return t68, t95, t997


def filled(ax, X, Y, p, color, *, zorder=3, lines=True, alpha=1.0, three=False):
    """68% (dark) and 95% (light) credible regions of a grid posterior."""
    t68, t95, t997 = _hpd(p)
    light, dark = ps.RAMPS[color]
    if three:
        ax.contourf(X, Y, p, levels=[t997, t95], colors=[ps.lighten(light, 0.5)], zorder=zorder,
                    alpha=alpha)
    ax.contourf(X, Y, p, levels=[t95, t68], colors=[light], zorder=zorder, alpha=alpha)
    ax.contourf(X, Y, p, levels=[t68, p.max() * 1.01], colors=[dark], zorder=zorder, alpha=alpha)
    if lines:
        ax.contour(X, Y, p, levels=[t95, t68], colors=[color], linewidths=0.8, zorder=zorder + 0.1)


def outline(ax, X, Y, p, color, *, ls="--", lw=1.0, zorder=6, masses=(0.683, 0.954)):
    levels = sorted(lk.credible_levels(p, masses))
    ax.contour(X, Y, p, levels=levels, colors=[color], linestyles=[ls], linewidths=lw,
               zorder=zorder)


def chain_density(x, y, w, xr, yr, bins=120, smooth=1.6):
    H, xe, ye = np.histogram2d(x, y, bins=bins, range=[xr, yr], weights=w)
    H = gaussian_filter(H, smooth)
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    return X, Y, H / H.sum()


def _plane_lines(ax, om_max, Or=0.0, *, recollapse=True, nbb=True, labels=True):
    om = np.linspace(0, om_max, 400)
    ax.plot(om, 1 - om - Or, color=ps.INK_SECONDARY, lw=0.9, zorder=4)
    ax.plot(om, om / 2 + Or, color=ps.INK_SECONDARY, lw=0.9, zorder=4)
    if recollapse:
        ax.plot(om, fr.recollapse_boundary(om, Or), color=ps.INK_SECONDARY, lw=0.9, zorder=4)
    if nbb:
        ax.plot(om, fr.no_big_bang_boundary(om, Or), color=ps.INK_SECONDARY, lw=0.9, zorder=4)


# ---------------------------------------------------------------------------
# 01  Every universe, all expanding at today's rate
# ---------------------------------------------------------------------------


def fig_universes():
    arrays, g = _load("gallery")
    show = ["ours", "eds", "milne", "desitter", "closed", "negative", "bounce", "phantom"]
    items = {u["key"]: u for u in g["universes"]}
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    curves = {}
    for key in show:
        u = items[key]
        col = ps.FATE[_fate_key(u)]
        lw = 2.8 if key == "ours" else 1.5
        pieces = [(arrays[f"{key}_t{j}"], arrays[f"{key}_a{j}"]) for j in range(u["n_pieces"])]
        for t, a in pieces:
            keep = (t >= -36) & (t <= 46)
            ax.plot(t[keep], a[keep], color=col, lw=lw, zorder=5 if key == "ours" else 3,
                    solid_capstyle="round")
        curves[key] = pieces
    ax.axvline(0, color=ps.AXIS, lw=0.8, zorder=1)
    ax.plot(0, 1, "o", color=ps.INK, ms=6, zorder=7, mec=ps.SURFACE, mew=1.2)
    ps.annotate(ax, "today: every universe here is the same size\nand expanding at the same rate, "
                f"H0 = {g['H0']:.1f} km/s/Mpc", (0, 1), (2.5, 0.42), fontsize=8)

    def at(key, t):  # a(t) on the expanding branch through today
        tt, aa = curves[key][0] if key != "bounce" else curves[key][1]
        return float(np.interp(t, tt, aa))

    def lab(text, xy, **kw):
        ax.text(*xy, text, fontsize=8.2, color=ps.INK_SECONDARY, zorder=8, **kw)

    t_ours = float(np.interp(4.1, curves["ours"][0][1], curves["ours"][0][0]))
    ax.text(t_ours + 0.9, 4.1, "our universe", fontsize=9.5, color=ps.INK, fontweight="bold",
            ha="left", va="center", zorder=8)
    rip = items["phantom"]["t_rip_gyr"]
    ps.annotate(ax, f"phantom dark energy (w = −1.5):\nBig Rip in {rip:.0f} Gyr",
                (12.0, at("phantom", 12.0)), (-14.0, 3.85), fontsize=8.2)
    lab("bounce: shrank, stopped, re-expanded;\nno Big Bang", (-34.0, 1.12), va="bottom")
    lab("de Sitter (Λ only): no beginning", (-34.0, 0.29), va="bottom")
    lab("empty (Milne)", (45.0, at("milne", 45.0) + 0.06), ha="right", va="bottom")
    lab("Einstein–de Sitter", (45.0, at("eds", 45.0) - 0.08), ha="right", va="top")
    lab("negative Λ: Big Crunch", (31.0, at("negative", 31.0) - 0.07), ha="left", va="top")
    lab("closed: Big Crunch", (16.0, 1.42), ha="center", va="top")
    # Big Bang markers on the time axis
    for key in ("ours", "eds", "milne", "closed", "negative", "phantom"):
        u = items[key]
        if isinstance(u.get("t_bang_gyr"), (int, float)):
            ax.plot(u["t_bang_gyr"], 0, marker="|", color=ps.FATE[_fate_key(u)], ms=10, mew=1.6,
                    zorder=6, clip_on=False)
    ax.set_xlim(-35, 45)
    ax.set_ylim(0, 4.5)
    ax.set_xlabel("time from today (billions of years)")
    ax.set_ylabel("size of the universe, relative to today (scale factor a)")
    ax.set_title("Build a universe: the same expansion today, very different pasts and futures")
    handles = [Line2D([], [], color=ps.FATE[k], lw=2, label=lbl) for k, lbl in
               (("accelerating", "expands forever, speeding up"),
                ("slowing", "expands forever, slowing down"),
                ("ends", "ends: Big Crunch or Big Rip"),
                ("no big bang", "had no Big Bang"))]
    ax.legend(handles=handles, loc="upper left", fontsize=8)
    ps.caption(fig, "Solutions of the Friedmann equation, all scaled to the expansion rate we measure "
               "today (our combined fit of supernovae, BAO and the CMB). Ticks on the time axis mark "
               "each universe's Big Bang: its age is how far back the curve hits zero. Ages scale as "
               "1/H0. Labelled simulations of textbook universes; ours is the measured one.")
    return ps.save(fig, FIG / "01_universes")


# ---------------------------------------------------------------------------
# 02  The gallery, one universe per panel
# ---------------------------------------------------------------------------


def _fmt_t(x):
    return f"{x:.1f}" if abs(x) < 100 else f"{x:.0f}"


def fig_gallery():
    arrays, g = _load("gallery")
    fig, axes = plt.subplots(3, 4, figsize=(13.5, 9.2), sharex=True, sharey=True)
    for ax, u in zip(axes.ravel(), g["universes"]):
        key = u["key"]
        col = ps.FATE[_fate_key(u)]
        for j in range(u["n_pieces"]):
            ax.plot(arrays[f"{key}_t{j}"], arrays[f"{key}_a{j}"], color=col, lw=1.8)
        ax.axvline(0, color=ps.AXIS, lw=0.7)
        ax.plot(0, 1, "o", color=ps.INK, ms=3.5, zorder=5)
        title = u["name"]
        params = (f"Ωm = {u['Om']:.2f}, ΩΛ = {u['OL']:.2f}" + (f", w = {u['w']:.1f}"
                                                              if u["w"] != -1 else "")
                  + (", Ωr = 1" if u["Or"] == 1.0 else ""))
        ax.set_title(title, fontsize=9.5)
        lines = [params]
        past = u["past"]
        if past == "Big Bang":
            lines.append(f"Big Bang {_fmt_t(u['age_gyr'])} Gyr ago")
        elif "a_min" in u:
            lines.append(f"no Big Bang: bounced {abs(u['t_bounce_gyr']):.1f} Gyr ago")
            lines.append(f"at {u['a_min']:.2f}× today's size")
        else:
            lines.append("no beginning: infinitely old")
        fut = u["future"]
        if fut == "Big Crunch":
            lines.append(f"Big Crunch in {u['t_crunch_gyr']:.1f} Gyr")
        elif fut == "Big Rip":
            lines.append(f"Big Rip in {u['t_rip_gyr']:.1f} Gyr")
        else:
            lines.append(fut.replace("expands forever, ", "forever, "))
        ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes, va="top", ha="left",
                fontsize=7.6, color=ps.INK_SECONDARY, linespacing=1.4, zorder=1,
                bbox={"fc": ps.SURFACE, "ec": "none", "alpha": 0.85, "pad": 1.0})
        ax.set_xlim(-40, 60)
        ax.set_ylim(0, 4)
        ax.set_xticks([-40, -20, 0, 20, 40, 60])
    for ax in axes[-1]:
        ax.set_xlabel("time from today (Gyr)")
    for ax in axes[:, 0]:
        ax.set_ylabel("scale factor a")
    fig.suptitle("A gallery of Friedmann universes, all expanding at today's rate", x=0.02,
                 ha="left", y=0.995)
    ps.caption(fig, "Each panel solves the Friedmann equation for one set of contents. Colour: blue "
               "expands forever and speeds up, aqua expands forever but slows, orange ends (Big "
               "Crunch or Big Rip), grey had no Big Bang. The dot is today. Loitering and bounce "
               "examples use Ωm = 0.1, where Petrosian, Salpeter & Szekeres (1967) put a Lemaître "
               "universe to explain a pile-up of quasar redshifts near z = 2.", y=0.0)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    return ps.save(fig, FIG / "02_gallery")


# ---------------------------------------------------------------------------
# 03  The map of universes, and why there was a Big Bang
# ---------------------------------------------------------------------------


def fig_map():
    A, s = _load("map")
    res = _results()
    Or = s["Or"]
    fig = plt.figure(figsize=(13.2, 6.6))
    ax = fig.add_axes([0.05, 0.1, 0.52, 0.8])
    ax2 = fig.add_axes([0.66, 0.1, 0.32, 0.8])
    Om, OL, code = A["Om"], A["OL"], A["code"]
    names = ["accelerating", "slowing", "ends", "no big bang"]
    cols = [ps.FATE_WASH[n] for n in names]
    from matplotlib.colors import ListedColormap
    ax.pcolormesh(Om, OL, code.T, cmap=ListedColormap(cols), vmin=-0.5, vmax=3.5, shading="auto",
                  rasterized=True, zorder=0)
    om = A["line_Om"]
    ax.plot(om, 1 - om - Or, color=ps.INK_SECONDARY, lw=0.9)
    ax.plot(om, om / 2 + Or, color=ps.INK_SECONDARY, lw=0.9)
    ax.plot(om, A["recollapse"], color=ps.INK, lw=1.2)
    ax.plot(om, A["nbb"], color=ps.INK, lw=1.2)
    # sparkline glyphs: a(t) from 30 Gyr ago to 40 Gyr ahead, a from 0 to 3; each glyph's
    # "today" dot sits exactly on the universe it shows
    w_box, h_box = 0.40, 0.36
    for i, om0 in enumerate(A["lat_Om"]):
        for j, ol0 in enumerate(A["lat_OL"]):
            k = 0
            pieces = []
            while f"spark_{i}_{j}_t{k}" in A:
                pieces.append((A[f"spark_{i}_{j}_t{k}"], A[f"spark_{i}_{j}_a{k}"]))
                k += 1
            if not pieces:
                continue
            x0, y0 = om0 - 30 / 70 * w_box, ol0 - h_box / 3
            c = int(fr.classify_grid([om0], [ol0], Or)[0])
            col = ps.FATE[names[c]]
            ax.plot([x0, x0 + w_box], [y0, y0], color=ps.AXIS, lw=0.5, zorder=4)
            for t, a in pieces:
                xx = x0 + (t + 30) / 70 * w_box
                yy = y0 + np.where(a <= 3.0, a, np.nan) / 3 * h_box
                ax.plot(xx, yy, color=col, lw=1.1, zorder=5, solid_capstyle="round")
            ax.plot(om0, ol0, "o", ms=2.4, color=ps.INK, zorder=6)
    # region labels, in gaps between glyphs
    for txt, xy in (("expands forever,\nspeeding up today", (0.5, 1.5)),
                    ("expands forever,\nslowing today", (2.5, 0.5)),
                    ("Big Crunch", (2.0, -0.5)),
                    ("no Big Bang:\na bounce", (0.5, 2.5))):
        ax.text(*xy, txt, fontsize=8.5, color=ps.INK, ha="center", va="center", zorder=7,
                bbox={"boxstyle": "round,pad=0.25", "fc": ps.SURFACE, "ec": "none", "alpha": 0.9})
    # line labels, rotated to follow the lines on screen
    def along(x, y, slope, text):
        p1 = ax.transData.transform((x, y))
        p2 = ax.transData.transform((x + 0.1, y + 0.1 * slope))
        ang = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
        ax.text(x, y, text, rotation=ang, rotation_mode="anchor", fontsize=7.5,
                color=ps.INK_SECONDARY, ha="left", va="bottom", zorder=7)

    ax.set_xlim(0, 3)
    ax.set_ylim(-1, 3)
    along(1.55, -0.55 + 0.02, -1.0, "flat space")
    along(2.25, 1.125 + 0.02, 0.5, "speeding up ↑  slowing down ↓")
    ours = res["ours"]["params"]
    ax.plot(ours["Om"], ours["OL"], "*", color=ps.INK, ms=14, mec=ps.SURFACE, mew=1.0, zorder=9)
    ps.annotate(ax, "our universe", (ours["Om"], ours["OL"]), (0.62, 0.95), color=ps.INK,
                fontsize=8.5)
    for lbl, xy, xyt in (("Einstein–de Sitter", (1.0, 0.0), (1.04, -0.17)),
                         ("empty", (0.0, 0.0), (0.04, -0.17)),
                         ("de Sitter", (0.0, 1.0), (0.04, 0.88))):
        ax.plot(*xy, "o", color=ps.INK, ms=4, zorder=9, clip_on=False)
        ax.text(*xyt, lbl, fontsize=7.5, color=ps.INK_SECONDARY, zorder=9)
    ax.set_xlabel("Ωm   matter density, as a fraction of the critical density")
    ax.set_ylabel("ΩΛ   cosmological-constant density")
    ax.set_title("The map of possible universes: every point is one universe, every glyph its a(t)")
    ax.grid(False)
    # right: the highest redshift a bouncing universe can ever show
    a_star = A["edge_a_star"]
    om_c = A["edge_Om"]
    keep = om_c > 0
    ax2.plot(om_c[keep], 1 / a_star[keep] - 1, color=ps.INK, lw=1.8,
             label="with radiation (the CMB's photons)")
    a_nr = np.geomspace(1e-4, 0.9, 400)
    om_nr, _ = fr.double_root_curve(a_nr, 0.0)
    ax2.plot(om_nr, 1 / a_nr - 1, color=ps.INK_SECONDARY, lw=1.2, ls=(0, (4, 2)),
             label="ignoring radiation")
    zcap = s["max_bounce_z_radiation"]
    ps.annotate(ax2, f"radiation caps every bounce:\nnothing beyond z = {zcap:.1f}",
                (3e-5, zcap), (2e-5, 2.2), fontsize=7.5, color=ps.INK)
    for key, yoff, va in (("quasar", 0.84, "top"), ("galaxy", 1.1, "bottom"),
                          ("cmb", 1.1, "bottom")):
        o = cfg.DISTANT_OBJECTS[key]
        ax2.axhline(o["z"], color=ps.INK_SECONDARY, lw=1.0)
        ax2.text(2e-5, o["z"] * yoff, f"{o['label']}, z = {o['z']:.4g}", fontsize=7.5,
                 color=ps.INK_SECONDARY, va=va)
    ax2.axvspan(0.049, 3, color=ps.FATE_WASH["no big bang"], zorder=0, lw=0)
    ax2.text(0.058, 160, "atoms alone\n(Ωb ≈ 0.049)\nforbid any bounce\nbeyond z = "
             f"{s['max_bounce_z_atoms']:.1f}", fontsize=7.3, color=ps.INK_SECONDARY, va="center")
    zm = s["max_bounce_z_measured_Om"]
    ax2.plot(0.30, zm, "*", color=ps.INK, ms=11, mec=ps.SURFACE, zorder=6)
    ps.annotate(ax2, f"measured Ωm = 0.30:\nno bounce beyond z = {zm:.2f}", (0.30, zm), (0.06, 0.18),
                fontsize=7.3, color=ps.INK)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(1e-5, 3)
    ax2.set_ylim(0.1, 3000)
    ax2.set_xlabel("Ωm of the bouncing universe")
    ax2.set_ylabel("largest redshift it could ever show us")
    ax2.set_title("Did the universe have a beginning?")
    ax2.legend(loc="upper right", fontsize=7.5)
    ps.caption(fig, f"Left: the fate of every universe with matter and a cosmological constant "
               f"(radiation included, Ωr = {Or:.1e}). Lines: flat space, the edge between speeding up "
               "and slowing down today, the recollapse edge (Carroll, Press & Turner 1992 eq. 11) and "
               "the no-Big-Bang edge (eq. 12, here with radiation). Each glyph is that universe's "
               "a(t) from 30 Gyr ago to 40 Gyr ahead; its dot (today) sits on the universe it shows. "
               "Right: a universe that bounced at size a_min shows nothing beyond z = 1/a_min − 1, "
               "and a_min is smallest on the no-Big-Bang edge itself, so this curve is the largest "
               "redshift any bounce can show. We see quasars at z = 7.6, galaxies at z = 14.4 and "
               "the CMB at z = 1090: our universe had a Big Bang.", y=0.0)
    return ps.save(fig, FIG / "03_map_of_universes")


# ---------------------------------------------------------------------------
# 04  Three probes in the plane
# ---------------------------------------------------------------------------


def _combo_grid(name):
    A, s = _load(f"combo_{name}")
    OMg = np.repeat(A["Om"][:, None], A["Ok"].size, axis=1)
    return OMg, A["OL"], A["p"], A, s


def fig_three_probes():
    sn, _ = _load("sn")
    bao, _ = _load("bao")
    cmb, cs = _load("cmb_band_pr4")
    res = _results()
    fig, (ax, az) = plt.subplots(1, 2, figsize=(13, 6.3), gridspec_kw={"width_ratios": [1.25, 1]})
    # wide
    OMs, OLs = np.meshgrid(sn["Om"], sn["OL"], indexing="ij")
    filled(ax, OMs, OLs, sn["p"], ps.C_SN, zorder=2)
    OMb = np.repeat(bao["Om"][:, None], bao["OL"].shape[1], axis=1)
    filled(ax, OMb, bao["OL"], bao["p"], ps.C_BAO, zorder=3)
    OMc = np.repeat(cmb["Om"][:, None], cmb["OL"].shape[1], axis=1)
    filled(ax, OMc, cmb["OL"], cmb["p"], ps.C_CMB, zorder=4)
    X, Y, P, _, _ = _combo_grid("sn+bao+cmb")
    filled(ax, X, Y, P, ps.C_COMB, zorder=5, lines=False)
    _plane_lines(ax, 1.2, recollapse=True, nbb=True)
    ax.text(0.98, 0.0 + 0.03, "flat", fontsize=7.5, color=ps.INK_SECONDARY, rotation=-38)
    ax.text(1.0, 0.53, "accelerating ↑\ndecelerating ↓", fontsize=7.5, color=ps.INK_SECONDARY,
            rotation=22, ha="center")
    ax.text(0.06, 1.47, "no Big Bang", fontsize=7.5, color=ps.INK_SECONDARY)
    for txt, xy, col in (("supernovae\n(Pantheon+, D1)", (0.08, 0.42), ps.C_SN),
                         ("BAO\n(DESI DR2)", (0.30, 0.86), ps.C_BAO),
                         ("CMB\n(Planck)", (0.66, 0.53), ps.C_CMB)):
        ax.text(*xy, txt, color=col, fontsize=8.5, fontweight="bold", ha="center")
    ax.set_xlim(0, 1.2)
    ax.set_ylim(-0.2, 1.6)
    ax.set_xlabel("Ωm   matter density")
    ax.set_ylabel("ΩΛ   dark-energy density")
    ax.set_title("Three independent ways to weigh the universe")
    # zoom
    filled(az, OMb, bao["OL"], bao["p"], ps.C_BAO, zorder=2)
    filled(az, OMc, cmb["OL"], cmb["p"], ps.C_CMB, zorder=3)
    filled(az, OMs, OLs, sn["p"], ps.C_SN, zorder=1)
    for name, col in (("bao+cmb", ps.INK),):
        X, Y, P, _, _ = _combo_grid(name)
        filled(az, X, Y, P, ps.C_COMB, zorder=5, lines=True)
    X, Y, P, _, _ = _combo_grid("sn+bao+cmb")
    outline(az, X, Y, P, ps.INK, ls="-", lw=1.4, masses=(0.683,))
    try:
        ch = fd.load_desi_chain(f"desi-bao-all_{cfg.CMB_FULL}")
        Xd, Yd, Pd = chain_density(ch["omegam"], ch["omegal"], ch["weight"], (0.27, 0.34),
                                   (0.64, 0.73))
        outline(az, Xd, Yd, Pd, ps.MAGENTA, ls="--", lw=1.1)
        chc = fd.load_desi_chain(cfg.CMB_FULL)
        Xc, Yc, Pc = chain_density(chc["omegam"], chc["omegal"], chc["weight"], (0.2, 0.6),
                                   (0.4, 0.85), bins=140, smooth=2.0)
        outline(az, Xc, Yc, Pc, ps.C_CMB, ls=":", lw=1.2)
    except FileNotFoundError:
        pass
    om = np.linspace(0.2, 0.42, 50)
    az.plot(om, 1 - om, color=ps.INK_SECONDARY, lw=0.9, zorder=6)
    az.text(0.355, 0.648, "flat", fontsize=7.5, color=ps.INK_SECONDARY, rotation=-45)
    c = res["combos"]["sn+bao+cmb"]
    bc = res["combos"]["bao+cmb"]
    az.text(0.263, 0.757, f"BAO + CMB:  Ωk = {bc['Ok']['mean']:+.4f} ± {bc['Ok']['std']:.4f}\n"
            f"all three:   Ωk = {c['Ok']['mean']:+.4f} ± {c['Ok']['std']:.4f}",
            fontsize=8, va="top", color=ps.INK)
    handles = [Patch(color=ps.RAMPS[ps.C_BAO][1], label="DESI BAO alone"),
               Patch(color=ps.RAMPS[ps.C_CMB][1], label="CMB alone (compressed Planck PR4)"),
               Patch(color=ps.RAMPS[ps.C_SN][1], label="supernovae alone"),
               Patch(color=ps.RAMPS[ps.INK][1], label="BAO + CMB (ours)"),
               Line2D([], [], color=ps.INK, lw=1.4, label="all three, 68% (ours)"),
               Line2D([], [], color=ps.MAGENTA, lw=1.1, ls="--", label="DESI's BAO + full CMB chain"),
               Line2D([], [], color=ps.C_CMB, lw=1.2, ls=":", label="DESI's full-CMB-alone chain")]
    az.legend(handles=handles, loc="lower left", fontsize=7.3, frameon=True, facecolor=ps.SURFACE, edgecolor="none", framealpha=0.92)
    az.set_xlim(0.26, 0.37)
    az.set_ylim(0.62, 0.76)
    az.set_xlabel("Ωm")
    az.set_ylabel("ΩΛ")
    az.set_title("Where they cross: a flat universe")
    ps.caption(fig, "Credible regions (dark 68%, light 95%) with curvature left free. Supernovae "
               "measure how the expansion rate changed (D1's Pantheon+ likelihood); BAO measure "
               "distances in units of the sound horizon (DESI DR2, ruler size free); the CMB fixes "
               "the matter densities and the angle of its sound horizon, which leaves a thin band "
               "(the geometric degeneracy). Each alone allows curved universes; any two cross on "
               "the flat line. Dashed/dotted: DESI's own chains with the full Planck likelihood "
               "(incl. CMB lensing), for comparison.")
    fig.tight_layout()
    return ps.save(fig, FIG / "04_three_probes")


# ---------------------------------------------------------------------------
# 05  How flat is space?
# ---------------------------------------------------------------------------


def fig_curvature():
    res = _results()
    rows = [("supernovae alone", res["sn"]["Ok"], ps.C_SN),
            ("DESI BAO alone", res["bao"]["Ok"], ps.C_BAO),
            ("supernovae + BAO", res["combos"]["sn+bao"]["Ok"], ps.INK_SECONDARY),
            ("supernovae + CMB", res["combos"]["sn+cmb"]["Ok"], ps.INK_SECONDARY),
            ("BAO + CMB", res["combos"]["bao+cmb"]["Ok"], ps.INK),
            ("BAO + CMB, D3's own CMB fit", res["combos"].get("bao+cmb_d3", {}).get("Ok"), ps.INK),
            ("all three", res["combos"]["sn+bao+cmb"]["Ok"], ps.INK)]
    rows = [r for r in rows if r[1] is not None]
    ref = [("Planck 2018 alone (TT,TE,EE+lowE)", cfg.PLANCK18["Ok_TTTEEE_lowE"]),
           ("Planck 2018 + lensing", cfg.PLANCK18["Ok_lensing"]),
           ("Planck 2018 + lensing + BAO", cfg.PLANCK18["Ok_lensing_BAO"]),
           ("DESI DR2: full CMB alone", cfg.DESI_DR2_OMK["CMB (full Planck + ACT lensing)"]),
           ("DESI DR2: BAO + full CMB", cfg.DESI_DR2_OMK["DESI BAO + CMB (full)"] + (0.0011,))]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True)
    n = len(rows) + len(ref)
    ys = np.arange(n)[::-1]
    labels = []
    for ax in (a1, a2):
        for k, (lbl, s, col) in enumerate(rows):
            y = ys[k]
            ax.plot([s["lo95"], s["hi95"]], [y, y], color=col, lw=1.0, solid_capstyle="round")
            ax.plot([s["lo68"], s["hi68"]], [y, y], color=col, lw=3.0, solid_capstyle="round")
            ax.plot(s["mean"], y, "o", color=col, ms=6, mec=ps.SURFACE, mew=1.2, zorder=5)
        for k, (lbl, v) in enumerate(ref):
            y = ys[len(rows) + k]
            val, up = v[0], v[1]
            lo = v[2] if len(v) > 2 else v[1]
            ax.plot([val - lo, val + up], [y, y], color=ps.MUTED, lw=3.0, solid_capstyle="round")
            ax.plot(val, y, "D", color=ps.SURFACE, mec=ps.MUTED, ms=5.5, mew=1.3, zorder=5)
        ax.axvline(0, color=ps.INK_SECONDARY, lw=0.8, zorder=0)
    labels = [r[0] for r in rows] + [r[0] for r in ref]
    a1.set_yticks(ys)
    a1.set_yticklabels(labels, fontsize=8)
    a1.set_xlim(-0.35, 0.45)
    a2.set_xlim(-0.012, 0.012)
    a1.axvspan(-0.012, 0.012, color=ps.GRID, alpha=0.6, zorder=0)
    for ax in (a1, a2):
        ax.set_xlabel("Ωk   (negative: closed, like a sphere; positive: open, like a saddle)")
        ax.grid(axis="y", visible=False)
    a1.set_title("How flat is space? Each probe alone, and together")
    a2.set_title("Zoomed ×30 (the grey strip on the left)")
    ours = res["combos"]["sn+bao+cmb"]["Ok"]
    s = abs(ours["mean"]) + 2 * ours["std"]
    H0 = res["ours"]["params"]["H0"]
    R_gly = cfg.C_KM_S / H0 / np.sqrt(s) * 3.26156e-3
    y_free = ys[len(rows)]  # the "Planck alone" row: off-scale in the zoom, so empty there
    a2.text(0.0118, y_free, f"all three: |Ωk| < {s:.4f} (95%). If space is curved, its\n"
            f"radius exceeds c/(H0 √|Ωk|) ≈ {R_gly:.0f} billion light-years,\n"
            f"{R_gly / 46.3:.1f}× the radius of the observable universe", ha="right", va="center",
            fontsize=7.8, color=ps.INK)
    for ax in (a1, a2):
        ax.axhline(ys[len(rows) - 1] - 0.5, color=ps.AXIS, lw=0.8)
    a1.text(-0.345, ys[0] + 0.55, "our fits", fontsize=7.5, color=ps.MUTED, style="italic")
    a1.text(-0.345, ys[len(rows)] + 0.55, "published", fontsize=7.5, color=ps.MUTED, style="italic")
    ps.caption(fig, "Dots: our posterior means, thick bars 68%, thin bars 95% (our likelihoods, "
                    "curvature free). Grey diamonds: published values (Planck 2018 VI eqs. 46b, "
                    "47a-b; DESI DR2 Table V, whose 'CMB' adds CMB lensing, which our compressed "
                    "prior leaves out).")
    fig.tight_layout()
    return ps.save(fig, FIG / "05_curvature")


# ---------------------------------------------------------------------------
# 06  Ages across the plane, and the oldest stars
# ---------------------------------------------------------------------------


def fig_ages():
    A, s = _load("ages")
    M, ms = _load("map")
    res = _results()
    tmin = s["oldest_stars"]["t_min_gyr"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.4), gridspec_kw={"width_ratios": [1, 1, 0.95]},
                             layout="constrained")
    Om, OL = M["Om"], M["OL"]
    OMg, OLg = np.meshgrid(Om, OL, indexing="ij")
    sel = (OMg <= 1.6) & (OLg <= 2.0)
    levels = [6, 8, 10, 12, 14, 16, 20, 30]
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cmap = ListedColormap(ps.SEQ_BLUE + ["#0a2a54"])
    cmap.set_over("#0a2a54")
    bounds = [4, 6, 8, 10, 12, 14, 16, 20, 30, 40]
    norm = BoundaryNorm(bounds, cmap.N)
    for ax, key, title in ((axes[0], "cmb", "if H0 = {:.1f} (BAO + CMB)"),
                           (axes[1], "cepheids", "if H0 = {:.1f} (Cepheids, D1)")):
        age = A[f"age_{key}"]
        H0 = s[key]["H0"]
        agef = np.where(np.isfinite(age), age, np.nan).T
        mesh = ax.pcolormesh(Om, OL, agef, cmap=cmap, norm=norm, shading="auto", rasterized=True,
                             zorder=1)
        young = np.where(np.isfinite(age), age < tmin, False)
        ax.contourf(Om, OL, young.T.astype(float), levels=[0.5, 1.5], colors=[ps.SURFACE],
                    alpha=0.62, zorder=2)
        nb = np.where(np.isnan(age), 1.0, 0.0)
        ax.contourf(Om, OL, nb.T, levels=[0.5, 1.5], colors=[ps.FATE_WASH["no big bang"]], zorder=2)
        ax.contour(Om, OL, agef, levels=levels, colors=[ps.INK], linewidths=0.4, zorder=3,
                   alpha=0.6)
        ax.contour(Om, OL, np.where(np.isfinite(age), age, 0).T, levels=[tmin], colors=[ps.ORANGE],
                   linewidths=2.2, zorder=4)
        om = np.linspace(0, 1.6, 100)
        ax.plot(om, 1 - om, color=ps.INK_SECONDARY, lw=0.8, zorder=4)
        ours = res["ours"]["params"]
        ax.plot(ours["Om"], ours["OL"], "*", color=ps.INK, ms=13, mec=ps.SURFACE, mew=1.0, zorder=7)
        ax.plot(1, 0, "o", color=ps.SURFACE, mec=ps.INK, ms=5, zorder=7)
        ax.text(1.03, 0.04, f"Einstein–de Sitter\n{s[key]['eds_age_gyr']:.1f} Gyr", fontsize=7.5,
                color=ps.INK, zorder=7)
        ax.set_xlim(0, 1.6)
        ax.set_ylim(-0.6, 2.0)
        ax.set_title(title.format(H0))
        ax.set_xlabel("Ωm")
        ax.set_ylabel("ΩΛ")
        ax.grid(False)
        ax.text(0.05, 1.86, "no Big Bang", fontsize=7.5, color=ps.INK_SECONDARY, zorder=7)
        ax.text(0.98, -0.52, f"faded: younger than the oldest stars\n(orange line: {tmin:.1f} Gyr, B6)",
                fontsize=7.4, color=ps.INK, ha="right", zorder=7)
    cb = fig.colorbar(mesh, ax=[axes[0], axes[1]], extend="max", ticks=bounds, fraction=0.035,
                      pad=0.02)
    cb.set_label("age of the universe (Gyr)", color=ps.INK_SECONDARY)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_visible(False)
    # right: the age crisis
    ax = axes[2]
    om = np.linspace(0.05, 1.0, 200)
    t_flat = fr.age_grid(om, 1 - om, 0.0)
    t_open = fr.age_grid(om, 0 * om, 0.0)
    for H0, col, lab in ((s["cmb"]["H0"], ps.BLUE, "H0 = {:.1f} (CMB)"),
                         (s["cepheids"]["H0"], ps.ORANGE, "H0 = {:.1f} (Cepheids)"),
                         (cfg.AGE_CRISIS_1990S["H0_freedman94"][0], ps.INK_SECONDARY,
                          "H0 = {:.0f} (Freedman et al. 1994)")):
        tH = fr.hubble_time_gyr(H0)
        ax.plot(om, t_flat * tH, color=col, lw=1.8, label=lab.format(H0) + ", flat with Λ")
        ax.plot(om, t_open * tH, color=col, lw=1.0, ls=(0, (4, 2)))
    ax.axhspan(tmin - s["oldest_stars"]["t_min_err"], tmin + s["oldest_stars"]["t_min_err"],
               color=ps.lighten(ps.YELLOW, 0.6), zorder=0, lw=0)
    ax.axhline(tmin, color=ps.YELLOW, lw=1.4)
    ax.text(0.06, tmin + 0.2, f"oldest stars (B6): ≥ {tmin:.1f} ± {s['oldest_stars']['t_min_err']:.1f} Gyr",
            fontsize=7.5, color=ps.INK)
    bh = cfg.AGE_CRISIS_1990S["gc_age_bolte_hogan95"]
    ax.errorbar([0.97], [bh[0]], yerr=[bh[1]], fmt="s", color=ps.INK_SECONDARY, ms=4, capsize=2)
    ax.text(0.95, bh[0] + 0.3, "1995: clusters\n15.8 ± 2.1 Gyr", fontsize=7, ha="right",
            color=ps.INK_SECONDARY)
    ax.set_xlim(0.05, 1.02)
    ax.set_ylim(6, 20)
    ax.set_xlabel("Ωm")
    ax.set_ylabel("age of the universe (Gyr)")
    ax.set_title("The 1990s age crisis")
    ax.legend(loc="upper right", fontsize=7)
    ax.text(0.06, 6.4, "solid: flat with Λ; dashed: no Λ (open)", fontsize=7, color=ps.MUTED)
    ps.caption(fig, "Left, middle: the age of every Big Bang universe in the plane, for the expansion "
               "rate measured from the CMB + BAO and from the Cepheid distance ladder (D1). Orange "
               "line: universes exactly as old as the oldest stars (B6's globular clusters plus their "
               "formation time); faded: younger than their own stars, which is impossible; "
               "star: our universe. Right: before 1998 a matter-only universe (Ωm = 1) with "
               "H0 ≈ 70–80 was 8–9 Gyr old, younger than its clusters; dark energy resolved it.")
    return ps.save(fig, FIG / "06_ages")


# ---------------------------------------------------------------------------
# 07  Who ruled when
# ---------------------------------------------------------------------------


def fig_epochs():
    A, s = _load("ours")
    ep = s["epochs"]
    t = A["hist_t_gyr"] * 1e9  # years
    rad = A["hist_photons"] + A["hist_neutrinos"]
    mat = A["hist_dark matter"] + A["hist_atoms"]
    de = A["hist_dark energy"]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(12, 7.4), sharex=True,
                                  gridspec_kw={"height_ratios": [1.5, 1]})
    ax.stackplot(t, rad, mat, de, colors=[ps.YELLOW, ps.BLUE, ps.VIOLET],
                 edgecolor=ps.SURFACE, linewidth=0.0, alpha=0.9)
    for frac, name, x, y in ((rad, "radiation\n(light + neutrinos)", 3e2, 0.5),
                             (mat, "matter\n(dark + atoms)", 3e7, 0.5),
                             (de, "dark\nenergy", 4e10, 0.5)):
        ax.text(x, y, name, ha="center", va="center", fontsize=9, color=ps.SURFACE,
                fontweight="bold")
    for tt, lbl in ((ep["t_eq_yr"], f"matter overtakes radiation\n{ep['t_eq_yr'] / 1e3:,.0f} thousand "
                                    f"years (z = {ep['z_eq']:,.0f})"),
                    (ep["t_matter_lambda_gyr"] * 1e9, f"dark energy overtakes matter\n"
                                                       f"{ep['t_matter_lambda_gyr']:.1f} Gyr "
                                                       f"(z = {ep['z_matter_lambda']:.2f})")):
        ax.axvline(tt, color=ps.INK, lw=0.8)
        ax.text(tt * 1.25, 0.94, lbl, fontsize=7.6, va="top", color=ps.INK,
                bbox={"fc": ps.SURFACE, "ec": "none", "alpha": 0.85, "pad": 1.5})
    ax.set_ylim(0, 1)
    ax.set_ylabel("share of the total density")
    ax.set_xscale("log")
    ax.set_title("Who ruled when: radiation, then matter, then dark energy")
    ax.grid(False)
    # lower panel: the components separately, on a log scale
    for key, col, lbl in (("photons", ps.YELLOW, "photons"), ("neutrinos", ps.ORANGE, "neutrinos"),
                          ("dark matter", ps.BLUE, "dark matter"), ("atoms", ps.AQUA, "atoms"),
                          ("dark energy", ps.VIOLET, "dark energy")):
        ax2.plot(t, A[f"hist_{key}"], color=col, lw=1.6, label=lbl)
    ax2.set_yscale("log")
    ax2.set_ylim(1e-6, 1.5)
    ax2.set_ylabel("share (log scale)")
    ax2.set_xlabel("time since the Big Bang (years)")
    events = [(ep["t_star_yr"], f"the CMB is released\n{ep['t_star_yr'] / 1e3:,.0f} thousand years",
               "left"),
              (ep["t_acc_gyr"] * 1e9, f"expansion starts\nto speed up\n{ep['t_acc_gyr']:.1f} Gyr",
               "right"),
              (ep["age_gyr"] * 1e9, f"today\n{ep['age_gyr']:.2f} Gyr", "left")]
    for tt, lbl, side in events:
        ax2.axvline(tt, color=ps.INK_SECONDARY, lw=0.7, ls=(0, (3, 2)))
        ax2.text(tt * (1.15 if side == "left" else 0.87), 2e-6, lbl, fontsize=7.2,
                 color=ps.INK_SECONDARY, va="bottom", ha=side)
    ax2.legend(loc="lower left", fontsize=7.5, ncol=1)
    i_nu = np.argmin(np.abs(t - 2e9))
    ps.annotate(ax2, "the 0.06 eV neutrino slows down\nand starts to behave like matter",
                (t[i_nu], A["hist_neutrinos"][i_nu]), (3e5, 1.5e-3), fontsize=7.2)
    ax2.set_xlim(1.0, 1e11)
    ps.caption(fig, "Our universe (the combined fit, curvature free) run forward from 1 year after "
                    "the Big Bang to 100 billion years. Each ingredient dilutes at its own rate: "
                    "radiation as a⁻⁴, matter as a⁻³, a cosmological constant not at all, so they take "
                    "turns. One neutrino species carries 0.06 eV of mass, the minimum oscillation "
                    "experiments allow; it behaves like radiation until z ~ 100.")
    fig.tight_layout()
    return ps.save(fig, FIG / "07_epochs")


# ---------------------------------------------------------------------------
# 08  Our universe: past, future, and the expansion rate
# ---------------------------------------------------------------------------


def fig_ours():
    A, s = _load("ours")
    res = _results()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.6))
    tg = A["t_grid"]
    lo, hi = np.nanpercentile(A["curves"], [2.5, 97.5], axis=0)
    ax.fill_between(tg, lo, hi, color=ps.RAMPS[ps.BLUE][0], lw=0, label="95% range of the fit")
    ax.plot(A["best_t"], A["best_a"], color=ps.BLUE, lw=2.2, label="our universe")
    H0 = s["params"]["H0"]
    eds = fr.Universe(1.0, 0.0, H0=H0)
    t, a = eds.history(-14.5, 60, n=600)[0]
    ax.plot(t, a, color=ps.AQUA, lw=1.2, ls=(0, (4, 2)), label="Einstein–de Sitter, same H0")
    ax.axvline(0, color=ps.AXIS, lw=0.8)
    ep = s["epochs"]
    ax.plot(-ep["age_gyr"], 0, marker="|", color=ps.BLUE, ms=12, mew=2, clip_on=False)
    ps.annotate(ax, f"Big Bang, {ep['age_gyr']:.2f} ± {res['combos']['sn+bao+cmb']['age_gyr']['std']:.2f}"
                " Gyr ago", (-ep["age_gyr"], 0.02), (-13.5, 1.6))
    ps.annotate(ax, f"speeding up since {ep['lookback_acc_gyr']:.1f} Gyr ago", (-ep["lookback_acc_gyr"],
                float(np.interp(-ep["lookback_acc_gyr"], A["best_t"], A["best_a"]))), (-12, 2.4))
    fut = s["future"]
    ps.annotate(ax, f"twice today's size\nin {fut['doubling_time_gyr']:.1f} Gyr",
                (fut["doubling_time_gyr"], 2.0), (fut["doubling_time_gyr"] - 6, 4.2))
    ax.set_xlim(-15, 60)
    ax.set_ylim(0, 12)
    ax.set_xlabel("time from today (Gyr)")
    ax.set_ylabel("scale factor a")
    ax.set_title(f"Our universe: {ep['age_gyr']:.1f} billion years behind, no end ahead")
    ax.legend(loc="upper left", fontsize=8)
    # right: the expansion rate H(z) with data
    cc = fd.load_chronometers()
    z = np.linspace(0, 2.5, 300)
    p = s["params"]
    u = fr.Universe(p["Om"], p["OL"], Or=p["Or"], H0=p["H0"])
    ax2.plot(z, u.hubble(z) / (1 + z), color=ps.BLUE, lw=2.0, label="our universe")
    ax2.plot(z, fr.Universe(1.0, 0.0, H0=p["H0"]).hubble(z) / (1 + z), color=ps.AQUA, lw=1.2,
             ls=(0, (4, 2)), label="Einstein–de Sitter")
    ax2.plot(z, fr.Universe(0.0, 0.0, H0=p["H0"]).hubble(z) / (1 + z), color=ps.INK_SECONDARY,
             lw=1.0, ls=(0, (1, 2)), label="empty: coasting")
    ax2.errorbar(cc["z"], cc["H"] / (1 + cc["z"]), yerr=cc["err"] / (1 + cc["z"]), fmt="o",
                 color=ps.MUTED, ms=3.5, lw=0.8, capsize=0, label="cosmic chronometers (D8)")
    bao = fd.load_bao()
    rd = res["combos"]["sn+bao+cmb"].get("r_drag") or s["epochs"]["r_drag"]
    dh = bao["kind"] == "DH_over_rs"
    Hz = cfg.C_KM_S / (bao["value"][dh] * rd)
    eH = Hz * np.sqrt(np.diag(bao["cov"])[dh]) / bao["value"][dh]
    ax2.errorbar(bao["z"][dh], Hz / (1 + bao["z"][dh]), yerr=eH / (1 + bao["z"][dh]), fmt="s",
                 color=ps.INK, ms=4.5, lw=1.0, capsize=0, label="DESI DR2 BAO (D_H, our r_d)")
    za = ep["z_acc"]
    ya = float(u.hubble(za) / (1 + za))
    ps.annotate(ax2, f"slowest at z = {za:.2f}, {ep['lookback_acc_gyr']:.1f} Gyr ago:\n"
                "from here on the expansion speeds up", (za, ya), (0.95, 36), fontsize=7.8)
    ax2.set_xlim(0, 2.5)
    ax2.set_ylim(20, 110)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("expansion speed ȧ = H(z) / (1 + z)   (km/s/Mpc)")
    ax2.set_title("The universe's speedometer: slowing, then speeding up")
    ax2.legend(loc="upper left", fontsize=7.5)
    ps.caption(fig, f"Left: a(t) from our combined fit (H0 = {p['H0']:.2f}, Ωm = {p['Om']:.3f}, "
                    f"Ωk = {p['Ok']:+.4f}); the band holds 95% of 200 universes drawn from the "
                    "posterior. Right: ȧ falls with time (rises with z) while matter dominates and "
                    "rises again once dark energy takes over; the minimum is when the expansion "
                    "started to accelerate. Chronometer errors are the diagonal ones tabulated by "
                    "Moresco et al. (2022), drawn only; D8 fits them with their full covariance.")
    fig.tight_layout()
    return ps.save(fig, FIG / "08_our_universe")


# ---------------------------------------------------------------------------
# 09  Validation
# ---------------------------------------------------------------------------


def fig_validation():
    res = _results()
    bao, bs = _load("bao")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    ax = axes[0]
    OMb = np.repeat(bao["Om"][:, None], bao["OL"].shape[1], axis=1)
    filled(ax, OMb, bao["OL"], bao["p"], ps.C_BAO, lines=False)
    ch = fd.load_desi_chain("desi-bao-all")
    Xd, Yd, Pd = chain_density(ch["omegam"], ch["omegal"], ch["weight"], (0.24, 0.35), (0.52, 0.85))
    outline(ax, Xd, Yd, Pd, ps.INK, ls="--", lw=1.1)
    ax.set_xlim(0.25, 0.34)
    ax.set_ylim(0.56, 0.82)
    ax.set_xlabel("Ωm")
    ax.set_ylabel("ΩΛ")
    ax.set_title("BAO alone: ours (filled) vs DESI's chain (dashed)")
    sh = bs["shift_sigma"]
    ax.text(0.04, 0.04, "shifts of the means:\n" + "\n".join(f"{k}: {v:+.3f} σ" for k, v in sh.items()),
            transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY)
    # MCMC vs grid for BAO + CMB
    ax = axes[1]
    X, Y, P, _, _ = _combo_grid("bao+cmb")
    filled(ax, X, Y, P, ps.C_COMB, lines=False)
    mc, _ = _load("mcmc_bao_cmb")
    Xm, Ym, Pm = chain_density(mc["Om"], mc["OL"], np.ones(mc["Om"].size), (0.28, 0.33), (0.66, 0.72),
                               bins=90, smooth=1.5)
    outline(ax, Xm, Ym, Pm, ps.ORANGE, ls="-", lw=1.1)
    try:
        chd = fd.load_desi_chain(f"desi-bao-all_{cfg.CMB_FULL}")
        Xd, Yd, Pd = chain_density(chd["omegam"], chd["omegal"], chd["weight"], (0.28, 0.33),
                                   (0.66, 0.72), bins=90, smooth=1.5)
        outline(ax, Xd, Yd, Pd, ps.MAGENTA, ls="--", lw=1.1)
    except FileNotFoundError:
        pass
    ax.set_xlim(0.285, 0.325)
    ax.set_ylim(0.67, 0.715)
    ax.set_xlabel("Ωm")
    ax.set_ylabel("ΩΛ")
    ax.set_title("BAO + CMB: grid (filled), MCMC (orange), DESI (dashed)")
    # CAMB deviations
    ax = axes[2]
    rows = res["checks"]["camb"]
    keys = [("theta_star", "θ*"), ("r_drag", "r_drag"), ("D_M", "D_M(z)"), ("H", "H(z)")]
    x = np.arange(len(rows))
    for j, (k, lbl) in enumerate(keys):
        vals = np.abs([r[k] for r in rows])
        ax.plot(x + 0.15 * (j - 1.5), np.log10(np.maximum(vals, 1e-12)), "o",
                color=ps.SERIES[j] if hasattr(ps, "SERIES") else [ps.BLUE, ps.ORANGE, ps.AQUA,
                                                                    ps.YELLOW][j], label=lbl, ms=6)
    ax.axhline(np.log10(2.5e-4), color=ps.MUTED, lw=0.8)
    ax.text(-0.4, np.log10(2.5e-4) + 0.1, "θ* precision of Planck (2.5×10⁻⁴)", fontsize=7,
            color=ps.MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Ωk={r['omk']:+.2f}\nH0={r['H0']:.0f}" for r in rows], fontsize=7)
    ax.set_ylabel("log10 |ours / CAMB − 1|")
    ax.set_ylim(-9, -2.5)
    ax.set_title("Our curved-space background vs CAMB")
    ax.legend(loc="upper right", fontsize=7.5, ncol=2)
    ps.caption(fig, "Left: our DESI BAO likelihood with curvature free against DESI's released "
                    "MCMC chain for the same run. Middle: the grid with the physical densities "
                    "marginalised by the Laplace method against a direct MCMC of the same likelihood, "
                    "and DESI's chain with the full CMB likelihood (which adds CMB lensing). Right: "
                    "distances, rulers and θ* for six open and closed universes against CAMB.")
    fig.tight_layout()
    return ps.save(fig, FIG / "09_validation")


FIGURES = {"01": fig_universes, "02": fig_gallery, "03": fig_map, "04": fig_three_probes,
           "05": fig_curvature, "06": fig_ages, "07": fig_epochs, "08": fig_ours,
           "09": fig_validation}


def main(argv=None) -> int:
    which = (argv if argv is not None else sys.argv[1:]) or list(FIGURES)
    for key in which:
        out = FIGURES[key]()
        print("  wrote", ", ".join(str(p).rsplit("/", 1)[-1] for p in out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
