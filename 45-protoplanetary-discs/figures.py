"""Every figure of the project, from data/results.json and the cached images.

    python figures.py            # all figures into figures/
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

import alma_images as A  # noqa: E402
import discs as D  # noqa: E402
import geometry as G  # noqa: E402
import literature as L  # noqa: E402
import planets as P  # noqa: E402
import plotstyle as ps  # noqa: E402
from config import (CACHE_DIR, CO_CUBE_PATH, FIGURE_DIR, HUANG_SPIRAL_RANGES, KURTOVIC_GEOMETRY, PROFILE_DIR,  # noqa: E402
                    RESULTS_PATH, SNOWLINES_K, SOLAR_SYSTEM_AU, TARGETS, TEAGUE2019_HD163296)

ps.use_style()


def load_results() -> dict:
    return json.loads(RESULTS_PATH.read_text())


def asinh_stretch(img: np.ndarray, peak: float, soft: float = 0.03) -> np.ndarray:
    v = np.clip(img / peak, 0.0, None)
    return np.arcsinh(v / soft) / np.arcsinh(1.0 / soft)


def disc_view(key: str, res: dict, half_au: float):
    """Image of a disc in au around its centre (Gaia distance), plus its orientation."""
    d = res["distances"][key]["distance"]
    if key in ("AS205", "HTLup"):
        full = A.load_continuum(key)
        xc, yc, _, _ = D.locate(full)
        inc, pa = KURTOVIC_GEOMETRY[key]
    else:
        o = res["discs"][key]
        g = o["geometry"]
        full = A.load_continuum(key)
        xc, yc, inc, pa = g["x0"], g["y0"], g["inc"], g["pa"]
    half = min(half_au / d, 0.5 * (full.x.max() - full.x.min()) - 0.05)
    img = full.cutout(half, xc, yc)
    return img, d, xc, yc, inc, pa


def orbit_xy(a_au: float, inc: float, pa: float, n: int = 361):
    """A circular orbit of radius a seen at the disc's inclination and PA (sky offsets in au)."""
    th = np.linspace(-180, 180, n)
    return G.sky_from_polar(a_au, th, 0.0, 0.0, inc, pa)


def _peak_of(img: A.Image, xc: float, yc: float, r: float) -> float:
    xx, yy = img.grids()
    sel = np.hypot(xx - xc, yy - yc) < r
    return float(np.nanmax(img.data[sel]))


# ---------------------------------------------------------------------------
# 1. The gallery on one scale (the keeper)
# ---------------------------------------------------------------------------
def gallery_same_scale(res: dict) -> None:
    keys = [t.key for t in TARGETS] + ["HLTau"]
    sizes = {}
    for k in keys:
        if k in res["discs"]:
            sizes[k] = res["discs"][k]["r_dust_au"]
        else:
            sizes[k] = 80.0  # the two multiples: ordered by eye among the mid-sized discs
    order = sorted(keys, key=lambda k: -sizes[k])
    half_au = 280.0
    ncol, nrow = 7, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(17.5, 8.4), facecolor=ps.DARK_PANEL)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.08, wspace=0.03, hspace=0.12)
    for ax, key in zip(axes.flat, order):
        img, d, xc, yc, inc, pa = disc_view(key, res, half_au)
        x_au = (img.x - xc) * d
        y_au = (img.y - yc) * d
        peak = _peak_of(img, xc, yc, 0.6)
        ax.imshow(asinh_stretch(img.data, peak), origin="lower", cmap=ps.IMAGE_CMAP, vmin=0, vmax=1,
                  extent=[x_au[0], x_au[-1], y_au[0], y_au[-1]], interpolation="bilinear")
        for name, a in SOLAR_SYSTEM_AU.items():
            ox, oy = orbit_xy(a, inc, pa)
            ax.plot(ox, oy, color=ps.C_SOLAR, lw=0.8 if name != "Neptune" else 1.2, alpha=0.9)
        ax.set_xlim(half_au, -half_au)  # east to the left
        ax.set_ylim(-half_au, half_au)
        ax.set_facecolor("black")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        name = res["discs"][key]["name"] if key in res["discs"] else {"AS205": "AS 205", "HTLup": "HT Lup"}[key]
        star = L.dsharp1_stars().get(key)
        sub = f"{d:.0f} pc" + (f" · {star.mstar:.1f} M$_\\odot$" if star else "")
        ax.text(0.04, 0.96, name, transform=ax.transAxes, color="white", fontsize=9.5, fontweight="bold", va="top")
        ax.text(0.04, 0.86, sub, transform=ax.transAxes, color="#c3c2b7", fontsize=7.5, va="top")
    # scale bar in the last panel
    ax = axes.flat[-1]
    ax.plot([-100 + 150, 150], [-235, -235], color="white", lw=1.6)
    ax.text(100, -215, "100 au", color="white", fontsize=8, ha="center")
    fig.text(0.01, 0.965, "Planets being born: 21 young discs on one scale", color="white", fontsize=15,
             fontweight="bold", ha="left", va="top")
    fig.text(0.01, 0.925, "ALMA 1.25 mm images (DSHARP; HL Tau at 1.0 mm), every panel 560 au across. White "
                          "curves: the orbits of Jupiter, Saturn, Uranus and Neptune, tilted like each disc.",
             color="#c3c2b7", fontsize=9.5, ha="left", va="top")
    fig.text(0.01, 0.035, "Brightness on an asinh stretch, normalised per disc. Distances: Gaia DR3 (HL Tau: VLBA, "
                          "147 pc; Elias 27: mean of its cloud). Orientations: our own fits (AS 205, HT Lup: DSHARP IV). "
                          "Data: ALMA 2016.1.00484.L (Andrews et al. 2018) and 2011.0.00015.SV.",
             color=ps.MUTED, fontsize=7.5, ha="left", va="top")
    ps.save(fig, FIGURE_DIR / "gallery_same_scale", formats=("png",))


# ---------------------------------------------------------------------------
# 2. Each disc zoomed, with our gaps and rings
# ---------------------------------------------------------------------------
def gallery_annotated(res: dict) -> None:
    keys = [t.key for t in TARGETS if not t.multiple] + ["HLTau"]
    order = sorted(keys, key=lambda k: -res["discs"][k]["r_dust_au"])
    ncol, nrow = 5, 4
    fig, axes = plt.subplots(nrow, ncol, figsize=(15.5, 13.2), facecolor=ps.DARK_PANEL)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.05, wspace=0.04, hspace=0.12)
    for ax in axes.flat:
        ax.set_visible(False)
    for ax, key in zip(axes.flat, order):
        ax.set_visible(True)
        o = res["discs"][key]
        half_au = 1.25 * o["r_dust_au"]
        img, d, xc, yc, inc, pa = disc_view(key, res, half_au)
        x_au, y_au = (img.x - xc) * d, (img.y - yc) * d
        peak = _peak_of(img, xc, yc, 0.6)
        ax.imshow(asinh_stretch(img.data, peak, 0.05), origin="lower", cmap=ps.IMAGE_CMAP, vmin=0, vmax=1,
                  extent=[x_au[0], x_au[-1], y_au[0], y_au[-1]], interpolation="bilinear")
        for f in o["features"]:
            if not f["secure"]:
                continue
            ox, oy = orbit_xy(f["r_au"], inc, pa)
            col = ps.C_GAP if f["kind"] == "D" else ps.C_RING
            ax.plot(ox, oy, color=col, lw=0.9, alpha=0.9, ls=":" if f["kind"] == "D" else "-")
        ax.set_xlim(half_au, -half_au)
        ax.set_ylim(-half_au, half_au)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        # beam (lower left) and scale bar (lower right)
        bx, by = 0.82 * half_au, -0.82 * half_au
        ax.add_patch(Ellipse((bx, by), o["beam"][1] * d, o["beam"][0] * d, angle=-o["beam"][2], color="white"))
        bar = 10 if half_au < 60 else (20 if half_au < 120 else 50)
        ax.plot([-0.85 * half_au + bar, -0.85 * half_au], [by, by], color="white", lw=1.4)
        ax.text(-0.85 * half_au + bar / 2, by + 0.06 * half_au, f"{bar} au", color="white", fontsize=7, ha="center")
        n_sec = sum(f["secure"] for f in o["features"])
        ax.text(0.04, 0.96, o["name"], transform=ax.transAxes, color="white", fontsize=10, fontweight="bold", va="top")
        ax.text(0.04, 0.88, f"i = {inc:.0f}°, {n_sec} gaps/rings found", transform=ax.transAxes, color="#c3c2b7",
                fontsize=7.5, va="top")
    fig.text(0.01, 0.985, "The gaps and rings we find, disc by disc", color="white", fontsize=15, fontweight="bold",
             va="top")
    fig.text(0.01, 0.955, "Dotted blue: gaps (local minima of the deprojected profile, > 5σ); solid orange: rings. "
                          "Each panel zoomed to 1.25 × the disc's dust radius; ALMA's beam (resolution) in the lower left.",
             color="#c3c2b7", fontsize=9.5, va="top")
    ps.save(fig, FIGURE_DIR / "gallery_annotated", formats=("png",))


# ---------------------------------------------------------------------------
# 3. Radial profiles vs DSHARP's
# ---------------------------------------------------------------------------
def profiles_grid(res: dict) -> None:
    keys = [t.key for t in TARGETS if not t.multiple] + ["HLTau"]
    order = sorted(keys, key=lambda k: -res["discs"][k]["r_dust_au"])
    feats_pub = L.dsharp2_features()
    fig, axes = plt.subplots(5, 4, figsize=(15, 17))
    fig.subplots_adjust(hspace=0.45, wspace=0.25, top=0.925, bottom=0.04)
    for ax in axes.flat:
        ax.set_visible(False)
    for ax, key in zip(axes.flat, order):
        ax.set_visible(True)
        o = res["discs"][key]
        prof = D.load_profile(key)
        d_ds = o["d_dsharp"]
        r_au = prof.r * d_ds
        peak = np.nanmax(prof.intensity)
        rmax = 1.15 * o["r_dust_au_dsharp"]
        sel = r_au <= rmax
        floor = max(o["noise"] / 3, 1e-4 * peak)
        y = np.clip(prof.intensity / peak, floor / peak, None)
        ax.fill_between(r_au[sel], np.clip((prof.intensity - prof.error)[sel] / peak, floor / peak, None),
                        np.clip((prof.intensity + prof.error)[sel] / peak, floor / peak, None), color=ps.GRID, lw=0)
        ax.plot(r_au[sel], y[sel], color=ps.C_OURS, lw=1.2, label="ours")
        pf = PROFILE_DIR / f"{key}.profile.txt"
        if pf.exists():
            rel = np.loadtxt(pf)
            ss = rel[:, 0] <= rmax
            ax.plot(rel[ss, 0], np.clip(rel[ss, 2] / peak, floor / peak, None), color=ps.C_PUBLISHED, lw=0.9,
                    alpha=0.9, label="DSHARP II")
            # the released profiles mostly use 1-au bins, MY Lup's 2-au bins: compare at our radii
            pub_i = np.interp(r_au, rel[:, 0], rel[:, 2], left=np.nan, right=np.nan)
            inner = (r_au <= o["r_dust_au_dsharp"]) & np.isfinite(pub_i) & np.isfinite(prof.intensity) & (pub_i > 0)
            frac = np.median(np.abs(prof.intensity[inner] - pub_i[inner]) / pub_i[inner])
            ax.text(0.03, 0.04, f"median difference {frac * 100:.1f}%", transform=ax.transAxes, ha="left",
                    fontsize=7.5, color=ps.INK_SECONDARY)
        ax.set_yscale("log")
        ax.set_ylim(max(floor / peak, 1e-3), 1.6)
        ax.set_xlim(0, rmax)
        for f in o["features"]:
            col = ps.C_GAP if f["kind"] == "D" else ps.C_RING
            mk = "v" if f["kind"] == "D" else "^"
            ax.plot(f["r_au_dsharp"], 1.25, marker=mk, color=col, ms=6 if f["secure"] else 4,
                    mfc=col if f["secure"] else "none", mew=1.0, ls="none")
        for f in [f for f in feats_pub if f.key == key and f.r_au.value is not None]:
            ax.plot(f.r_au.value, 0.95, marker="|", color=ps.C_PUBLISHED, ms=7, mew=1.0)
        ax.set_title(o["name"], fontsize=10)
        ax.set_xlabel("radius (au, DSHARP distance)", fontsize=8)
        if key == order[0]:
            ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.1), fontsize=7.5)
    fig.suptitle("Deprojected brightness profiles: ours (black) and DSHARP's released ones (violet)", x=0.02, ha="left")
    fig.text(0.02, 0.952, "Intensity / peak, log scale. Triangles: our gaps (blue ▼) and rings (orange ▲), filled when "
                          "> 5σ; violet ticks: DSHARP II's features. Grey band: ±1σ (scatter / √number of beams). "
                          "Inclined discs averaged in ±20° wedges along the major axis, as in DSHARP II.",
             fontsize=8.5, color=ps.INK_SECONDARY, ha="left")
    ps.save(fig, FIGURE_DIR / "profiles", formats=("png",))


# ---------------------------------------------------------------------------
# 4. Orientations vs DSHARP II
# ---------------------------------------------------------------------------
def geometry_comparison(res: dict) -> None:
    rows = res["geometry"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.8))
    fig.subplots_adjust(wspace=0.28, top=0.82, bottom=0.14)
    for key, r in rows.items():
        o, p = r["ours"], r["published"]
        mk = "o" if o["method"] == "E" else "s"
        a1.errorbar(p["inc"], o["inc"], xerr=p["inc_err"] or 0, yerr=o["inc_err"] or 0, fmt=mk, color=ps.C_OURS,
                    ms=5, mfc=ps.C_OURS if o["method"] == "E" else "white", elinewidth=0.8, capsize=0)
        dpa, spa = r["d_pa"], np.hypot(o["pa_err"] or 0, p["pa_err"] or 0)
        a2.errorbar(p["inc"], dpa, yerr=spa, fmt=mk, color=ps.C_OURS, ms=5,
                    mfc=ps.C_OURS if o["method"] == "E" else "white", elinewidth=0.8, capsize=0)
        if abs(r["d_inc"]) > 2.0 or (abs(dpa) > 10 and spa < 200):
            a1.annotate(r["name"], (p["inc"], o["inc"]), xytext=(5, -10), textcoords="offset points", fontsize=7.5,
                        color=ps.INK_SECONDARY)
        if abs(dpa) > 6 and spa < 200:
            a2.annotate(r["name"], (p["inc"], dpa), xytext=(5, 3), textcoords="offset points", fontsize=7.5,
                        color=ps.INK_SECONDARY)
    lim = [0, 80]
    a1.plot(lim, lim, color=ps.AXIS, lw=1.0, zorder=0)
    a1.set_xlim(lim)
    a1.set_ylim(lim)
    a1.set_xlabel("inclination, DSHARP II / ALMA 2015 (deg)")
    a1.set_ylabel("inclination, ours (deg)")
    a1.set_title("Disc inclinations")
    a2.axhline(0, color=ps.AXIS, lw=1.0, zorder=0)
    a2.set_ylim(-60, 60)
    a2.set_xlim(lim)
    a2.set_xlabel("inclination, DSHARP II (deg)")
    a2.set_ylabel("PA ours − published (deg)")
    a2.set_title("Position angles: hard to pin down for face-on discs")
    a1.plot([], [], "o", color=ps.C_OURS, label="ellipse fits to rings (E)")
    a1.plot([], [], "s", color=ps.C_OURS, mfc="white", label="2-D Gaussian (G)")
    a1.legend(loc="upper left")
    meas = [r for r in rows.values() if (r["ours"]["inc_err"] or 99) < 10]
    n1 = sum(1 for r in meas if abs(r["d_inc_sigma"] or 0) <= 2)
    fig.suptitle(f"Our own disc orientations vs the published ones: {n1} of the {len(meas)} measurable inclinations "
                 "within 2σ (Sz 114, nearly face-on, is unconstrained)", x=0.02, ha="left")
    ps.save(fig, FIGURE_DIR / "geometry_comparison", formats=("png",))


# ---------------------------------------------------------------------------
# 5. Gap and ring radii, widths and depths vs DSHARP II
# ---------------------------------------------------------------------------
def features_comparison(res: dict) -> None:
    rows = res["features"]["rows"]
    widths = res["widths"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    fig.subplots_adjust(wspace=0.3, top=0.8, bottom=0.15)
    a1, a2, a3 = axes
    for r in rows:
        if not r["found"]:
            continue
        col = ps.C_GAP if r["kind"] == "D" else ps.C_RING
        a1.errorbar(r["r_pub_mas"], r["dr_mas"] / r["beam_mas"],
                    yerr=np.hypot(r["r_ours_err_mas"] or 0, r["r_pub_err_mas"] if np.isfinite(r["r_pub_err_mas"] or np.nan) else 0) / r["beam_mas"],
                    fmt="o" if r["method"] == "E" else ("s" if r["method"] in ("R", "-") else "D"), color=col,
                    mfc=col if r["method"] == "E" else "white", ms=4.5, elinewidth=0.7, capsize=0)
    a1.axhline(0, color=ps.AXIS, lw=1)
    a1.set_xscale("log")
    a1.set_xlabel("radius, DSHARP II (mas)")
    a1.set_ylabel("ours − published (beams)")
    a1.set_ylim(-1, 1)
    bm = res["features"]["by_method"]
    txt = "  ".join(f"{m if m != '-' else 'HL Tau'}: {v['found']}/{v['n']}" for m, v in bm.items())
    a1.set_title("Where the gaps and rings are")
    a1.text(0.02, 0.04, "recovered by method — " + txt, transform=a1.transAxes, fontsize=7.5, color=ps.INK_SECONDARY)
    for w in widths:
        if w.get("gap_width_pub") is not None and w["gap_width_au"] is not None and not w.get("gap_width_pub_upper"):
            a2.errorbar(w["gap_width_pub"], w["gap_width_au"], xerr=w.get("gap_width_pub_err") or 0,
                        yerr=w["gap_width_err"] if np.isfinite(w["gap_width_err"] or np.nan) else 0, fmt="o",
                        color=ps.C_GAP, ms=4.5, elinewidth=0.7)
        if w.get("ring_width_pub") is not None and w["ring_width_au"] is not None and np.isfinite(w["ring_width_au"]):
            a2.errorbar(w["ring_width_pub"], w["ring_width_au"], xerr=w.get("ring_width_pub_err") or 0,
                        yerr=w["ring_width_err"] if np.isfinite(w["ring_width_err"] or np.nan) else 0, fmt="o",
                        color=ps.C_RING, ms=4.5, elinewidth=0.7)
        if w.get("depth_pub") is not None:
            a3.errorbar(w["depth_pub"], w["depth"], xerr=w.get("depth_pub_err") or 0,
                        yerr=w["depth_err"] if np.isfinite(w["depth_err"] or np.nan) else 0, fmt="o", color=ps.C_GAP,
                        ms=4.5, elinewidth=0.7)
    for ax, lim in ((a2, (0, 25)), (a3, (0.005, 1.05))):
        ax.plot(lim, lim, color=ps.AXIS, lw=1, zorder=0)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
    a3.set_xscale("log")
    a3.set_yscale("log")
    a2.set_xlabel("width, DSHARP II (au)")
    a2.set_ylabel("width, ours (au)")
    a2.set_title("Widths of gaps (blue) and rings (orange)")
    a3.set_xlabel("depth I$_{gap}$/I$_{ring}$, DSHARP II")
    a3.set_ylabel("depth, ours")
    a3.set_title("Gap depths")
    a1.plot([], [], "o", color=ps.C_GAP, label="gap")
    a1.plot([], [], "o", color=ps.C_RING, label="ring")
    a1.plot([], [], "o", color=ps.INK_SECONDARY, label="DSHARP: ellipse fit")
    a1.plot([], [], "s", color=ps.INK_SECONDARY, mfc="white", label="DSHARP: profile extremum")
    a1.legend(loc="upper left", fontsize=7)
    fig.suptitle(f"Our gaps and rings vs DSHARP II: {res['features']['n_found']} of {res['features']['n_published']} "
                 "published features recovered", x=0.02, ha="left")
    ps.save(fig, FIGURE_DIR / "features_comparison", formats=("png",))


# ---------------------------------------------------------------------------
# 6. The unseen planets among the known exoplanets
# ---------------------------------------------------------------------------
def planets_mass_separation(res: dict) -> None:
    ex = res["exoplanets"]["by_method"]
    pl = [r for r in res["planets"]["rows"] if r.get("found")]
    fig, ax = plt.subplots(figsize=(9.5, 6.6))
    fig.subplots_adjust(top=0.85, bottom=0.11, right=0.97)
    imaging = ex.get("Imaging", {"a": [], "m": []})
    for meth, g in ex.items():
        if meth == "Imaging":
            continue
        ax.scatter(g["a"], g["m"], s=6, color=ps.C_CONTEXT, alpha=0.35, lw=0, zorder=1)
    ax.scatter(imaging["a"], imaging["m"], s=22, marker="s", facecolor="white", edgecolor=ps.INK_SECONDARY, lw=0.8,
               zorder=2)
    solar = {"Earth": (1.0, 3.003e-3), "Jupiter": (5.20, 1.0), "Saturn": (9.54, 0.2994), "Uranus": (19.19, 0.0457),
             "Neptune": (30.07, 0.0540)}
    for n, (a, m) in solar.items():
        ax.scatter(a, m, s=40, marker="*", color=ps.YELLOW, edgecolor=ps.INK, lw=0.5, zorder=4)
        ax.annotate(n, (a, m), xytext=(5, -9), textcoords="offset points", fontsize=7.5, color=ps.INK_SECONDARY)
    for r in pl:
        ms = r["masses_dr3"]
        pref = ms.get("DSD1_0.001")
        if pref is None:
            continue
        lo = min(v["m"] for k, v in ms.items() if k.startswith("DSD1"))
        hi = max(v["m"] for k, v in ms.items() if k.startswith("DSD1"))
        a = r["r_gap_au"]
        ax.plot([a, a], [lo, hi], color=ps.C_OURS, lw=1.0, zorder=5)
        ax.scatter(a, pref["m"], s=36, color=ps.C_OURS, zorder=6)
        ax.scatter(r["r_gap_pub"], r["masses"]["DSD1_0.001"]["pub"], s=34, marker="D", facecolor="none",
                   edgecolor=ps.C_PUBLISHED, lw=1.0, zorder=6)
        off = {("AS209", 9): (-58, 6), ("SR4", 11): (6, 6), ("Elias24", 57): (7, 4), ("AS209", 99): (7, -2),
               ("HD143006", 51): (7, -4), ("HD142666", 16): (-70, -3), ("HD163296", 48): (7, 2)}
        dx, dy = off.get((r["disc"], round(r["r_gap_pub"])), (6, 2))
        ax.annotate(f"{r['name']} {r['r_gap_au']:.0f} au", (a, pref["m"]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=7, color=ps.INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.005, 3000)
    ax.set_ylim(5e-4, 60)
    ax.set_xlabel("distance from the star (au)")
    ax.set_ylabel("planet mass (Jupiter masses; m sin i for velocity planets)")
    ax.scatter([], [], s=6, color=ps.C_CONTEXT, label=f"known exoplanets ({res['exoplanets']['n']}, NASA Exoplanet Archive)")
    ax.scatter([], [], s=22, marker="s", facecolor="white", edgecolor=ps.INK_SECONDARY, label="directly imaged")
    ax.scatter([], [], s=36, color=ps.C_OURS, label="planets that would carve the DSHARP gaps (ours; bar: α = 10$^{-4}$–10$^{-2}$)")
    ax.scatter([], [], s=34, marker="D", facecolor="none", edgecolor=ps.C_PUBLISHED, label="the same, DSHARP VII (Zhang et al. 2018)")
    ax.scatter([], [], s=40, marker="*", color=ps.YELLOW, edgecolor=ps.INK, lw=0.5, label="Solar System")
    ax.legend(loc="lower right", fontsize=7.5)
    fig.suptitle("If planets carve the gaps: Neptunes to Jupiters, 10-100 au out", x=0.02, ha="left")
    fig.text(0.02, 0.905, "Gap widths measured on our own profiles, turned into masses with DSHARP VII's simulations "
                          "(dust grains up to 0.1 mm, disc viscosity α = 10$^{-3}$). The bars span α = 10$^{-4}$ to 10$^{-2}$: "
                          "\nthe viscosity is unknown, so these masses are uncertain by factors of several, and a gap "
                          "need not hold a planet at all.", fontsize=8, color=ps.INK_SECONDARY, ha="left")
    ps.save(fig, FIGURE_DIR / "planets_mass_separation", formats=("png",))


def planet_mass_check(res: dict) -> None:
    pl = [r for r in res["planets"]["rows"] if r.get("found")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6))
    fig.subplots_adjust(wspace=0.3, top=0.8, bottom=0.15)
    for r in pl:
        a1.errorbar(r["delta_pub"], r["delta"], yerr=r["delta_err"], fmt="o", color=ps.C_GAP, ms=5, elinewidth=0.8)
        a1.annotate(f"{r['name']} {r['r_gap_pub']:.0f}", (r["delta_pub"], r["delta"]), xytext=(4, -9),
                    textcoords="offset points", fontsize=6.5, color=ps.INK_SECONDARY)
        m = r["masses"]["DSD1_0.001"]
        a2.errorbar(m["pub"], m["m"], yerr=[[m["m"] - m["lo"]], [m["hi"] - m["m"]]], fmt="o", color=ps.C_OURS, ms=5,
                    elinewidth=0.8)
    a1.plot([0.05, 0.75], [0.05, 0.75], color=ps.AXIS, lw=1, zorder=0)
    a1.set_xlim(0.05, 0.75)
    a1.set_ylim(0.05, 0.75)
    a1.set_xlabel("gap width Δ, DSHARP VII")
    a1.set_ylabel("gap width Δ, ours")
    a1.set_title("Gap widths Δ = (r$_{out}$ − r$_{in}$)/r$_{out}$")
    lim = (0.01, 40)
    a2.plot(lim, lim, color=ps.AXIS, lw=1, zorder=0)
    a2.set_xscale("log")
    a2.set_yscale("log")
    a2.set_xlim(lim)
    a2.set_ylim(lim)
    a2.set_xlabel("planet mass, DSHARP VII (M$_J$)")
    a2.set_ylabel("planet mass from our widths (M$_J$)")
    a2.set_title("Planet masses (α = 10$^{-3}$, DSD1); bars: fit scatter")
    fig.suptitle(f"From gap widths to planets: our numbers vs DSHARP VII (median ratio "
                 f"{res['planets']['median_ratio']:.2f}, rms {res['planets']['log_rms']:.2f} dex)", x=0.02, ha="left")
    ps.save(fig, FIGURE_DIR / "planet_mass_check", formats=("png",))


# ---------------------------------------------------------------------------
# 7. Snow lines
# ---------------------------------------------------------------------------
def snowlines(res: dict) -> None:
    rows = res["snowlines"]["rows"]
    discs = sorted({r["disc"] for r in rows}, key=lambda k: -L.dsharp1_stars()[k].lstar *
                   (res["distances"][k]["distance"] / L.dsharp1_stars()[k].distance.value) ** 2)
    fig, ax = plt.subplots(figsize=(9, 6.2))
    fig.subplots_adjust(left=0.16, top=0.86, bottom=0.11)
    shades = {"CO$_2$": ps.RED, "CO (on water ice)": ps.MUTED, "CO (pure ice)": ps.MUTED, "N$_2$ (pure ice)": ps.YELLOW}
    for i, (name, (lo, hi)) in enumerate(SNOWLINES_K.items()):
        ax.axvspan(lo, hi, color=shades[name], alpha=0.16, lw=0, zorder=0)
        ax.text(np.sqrt(lo * hi), len(discs) - 0.35 - (1.6 if name.startswith("N") else 0), name, ha="center",
                fontsize=7, color=ps.INK_SECONDARY, rotation=90, va="top")
    for i, k in enumerate(discs):
        for r in [r for r in rows if r["disc"] == k]:
            col = ps.C_GAP if r["kind"] == "D" else ps.C_RING
            ax.scatter(r["t_mid"], i, marker="v" if r["kind"] == "D" else "^", color=col, s=36, zorder=3)
    ax.set_yticks(range(len(discs)))
    lum = {k: L.dsharp1_stars()[k].lstar * (res["distances"][k]["distance"] / L.dsharp1_stars()[k].distance.value) ** 2
           for k in discs}
    ax.set_yticklabels([f"{res['discs'][k]['name']}  ({lum[k]:.1f} L$_\\odot$)" for k in discs], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(8, 150)
    ax.set_ylim(-0.7, len(discs) - 0.3)
    ax.set_xlabel("midplane temperature at the feature (K): T = (φ L$_*$ / 8π r$^2$ σ)$^{1/4}$, φ = 0.02")
    ax.scatter([], [], marker="v", color=ps.C_GAP, label="gap")
    ax.scatter([], [], marker="^", color=ps.C_RING, label="ring")
    ax.legend(loc="lower right")
    f = res["snowlines"]["fraction_in_range"]
    fig.suptitle("Are gaps ice lines? Their temperatures don't cluster", x=0.02, ha="left")
    nul = res["snowlines"]["uniform_log_t_expectation"]
    fig.text(0.02, 0.905, f"Shaded: where CO$_2$, CO and N$_2$ freeze out (DSHARP II Table 3). Of our {res['snowlines']['n_gaps']} "
                          "secure gaps, " + ", ".join(f"{v * 100:.0f}% fall in the {k} band (≈{nul[k] * 100:.0f}% expected "
                                                       "by chance)" for k, v in f.items() if v is not None) + ".",
             fontsize=8, color=ps.INK_SECONDARY, ha="left", wrap=True)
    ps.save(fig, FIGURE_DIR / "snowlines", formats=("png",))


# ---------------------------------------------------------------------------
# 8. HD 163296: weighing the star with its disc
# ---------------------------------------------------------------------------
def _kin_model(k: dict):
    import kinematics as K

    f = k["fit"]
    near = -1 if k["rotation_on_sky"] == "clockwise" else 1
    return K.DiscModel(f["x0_mas"] / 1e3, f["y0_mas"] / 1e3, f["inc"], f["pa"], f["mstar"], f["vsys"], f["z0_mas"] / 1e3,
                       f["phi"], f["z1_mas"] / 1e3, f["psi"], distance=f["distance"], near=near)


def hd163296_kinematics(res: dict) -> None:
    k = res["kinematics"]
    if "fit" not in k:
        return
    lc = np.load(CACHE_DIR / "kinematics" / "linecentre.npz")
    mdl = np.load(CACHE_DIR / "kinematics" / "model.npz")
    m = _kin_model(k)
    x, y = lc["x"], lc["y"]
    ext = [x[0], x[-1], y[0], y[-1]]
    v0, vmod, resid = lc["v0"], mdl["vmodel"], mdl["resid"]
    f = k["fit"]
    d = res["distances"]["HD163296"]["distance"]
    fig = plt.figure(figsize=(15, 9.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.75], hspace=0.32, wspace=0.45, top=0.88, bottom=0.07)
    # (a) the line-centre map with model isovelocity contours (a "spider diagram")
    a = fig.add_subplot(gs[0, 0])
    lim = 3.2
    im = a.imshow(v0, origin="lower", extent=ext, cmap=ps.VELOCITY_CMAP, vmin=f["vsys"] - 4.5, vmax=f["vsys"] + 4.5,
                  interpolation="nearest")
    levels = f["vsys"] + np.arange(-4.0, 4.01, 0.5)
    vm_masked = np.where(mdl["r"] < 3.0, vmod, np.nan)
    a.contour(x, y, vm_masked, levels=levels, colors=ps.INK, linewidths=0.5)
    a.set_xlim(lim, -lim)
    a.set_ylim(-lim, lim)
    a.set_xlabel("RA offset (arcsec)")
    a.set_ylabel("Dec offset (arcsec)")
    a.set_title("CO line centres + model isovelocity curves")
    a.grid(False)
    cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.02)
    cb.set_label("velocity (km/s, LSRK)")
    # (b) residuals with the continuum gaps and Pinte et al.'s kinks
    b = fig.add_subplot(gs[0, 1])
    im2 = b.imshow(resid * 1e3, origin="lower", extent=ext, cmap=ps.VELOCITY_CMAP, vmin=-150, vmax=150,
                   interpolation="nearest")
    o = res["discs"]["HD163296"]
    for feat in o["features"]:
        if feat["kind"] == "D" and feat["secure"]:
            th = np.linspace(-np.pi, np.pi, 361)
            xs, ys = m.sky(feat["r"], th, 0.0)
            b.plot(xs, ys, color=ps.INK, lw=0.6, ls=":")
    for name, kk in k["kinks"].items():
        sep, pa_k = kk["sep"], kk["pa"]
        b.plot(sep * np.sin(np.radians(pa_k)), sep * np.cos(np.radians(pa_k)), "o", mfc="none", mec=ps.INK, ms=12, mew=1.2)
        b.annotate(name.replace("HD163296 ", "kink "), (sep * np.sin(np.radians(pa_k)), sep * np.cos(np.radians(pa_k))),
                   xytext=(-8, 10), textcoords="offset points", fontsize=8)
    b.set_xlim(lim, -lim)
    b.set_ylim(-lim, lim)
    b.set_xlabel("RA offset (arcsec)")
    b.set_title("What rotation does not explain")
    b.grid(False)
    cb2 = fig.colorbar(im2, ax=b, fraction=0.046, pad=0.02)
    cb2.set_label("data − model (m/s)")
    # (c) the emission surface
    c = fig.add_subplot(gs[0, 2])
    rr = np.linspace(0.05, 4.0, 200)
    t = TEAGUE2019_HD163296
    z_t = t["z0_mas"][0] / 1e3 * rr ** t["phi"][0] - t["z1_mas"][0] / 1e3 * rr ** t["psi"][0]
    c.plot(rr * d, m.height(rr) * d, color=ps.C_OURS, lw=1.8, label="ours")
    c.plot(rr * d, np.clip(z_t, 0, None) * d, color=ps.C_PUBLISHED, lw=1.4, label="Teague et al. 2019")
    for zr in (0.1, 0.2, 0.3):
        c.plot(rr * d, zr * rr * d, color=ps.GRID, lw=0.8, zorder=0)
        xr_lab = min(3.3 * d, 95.0 / zr)
        c.text(xr_lab, zr * xr_lab + 1, f"z/r = {zr}", fontsize=7, color=ps.MUTED, ha="right", va="bottom")
    c.set_xlim(0, 3.5 * d)
    c.set_ylim(0, 100)
    c.set_xlabel("radius (au)")
    c.set_ylabel("height of the CO emitting layer (au)")
    c.set_title("The disc's CO 'skin' is flared")
    c.legend(loc="upper left")
    # (d) the rotation curve and its deviations
    rc = k["rotation_curve"]
    r_au = np.array(rc["r"]) * d
    vphi, vk, verr = np.array(rc["v_phi"]), np.array(rc["v_kep"]), np.array(rc["v_phi_err"])
    e = fig.add_subplot(gs[1, :2])
    e.axhline(0, color=ps.AXIS, lw=1)
    e.errorbar(r_au, (vphi / vk - 1) * 100, yerr=verr / vk * 100, fmt="o", color=ps.C_OURS, ms=4, elinewidth=0.8)
    for feat in o["features"]:
        if feat["secure"]:
            e.axvline(feat["r_au"], color=ps.C_GAP if feat["kind"] == "D" else ps.C_RING, lw=1.0,
                      ls=":" if feat["kind"] == "D" else "-", alpha=0.8)
    for pm in k.get("pressure_minima_au", []):
        e.annotate("", (pm, 0), xytext=(pm, -6.2), arrowprops=dict(arrowstyle="-|>", color=ps.C_OURS, lw=1.0))
        e.text(pm, -7.2, f"{pm:.0f} au", ha="center", fontsize=7.5, color=ps.C_OURS)
    for tr in (83.0, 137.0):
        e.plot(tr * 100.6 / 101.0, 6.8, marker="v", color=ps.C_PUBLISHED, ms=7)
    e.text(0.99, 0.94, "▼ planets inferred by Teague et al. 2018 (83, 137 au)\n↑ our pressure minima (v crosses Kepler upwards)",
           transform=e.transAxes, ha="right", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    e.set_xlim(0, 3.6 * d)
    e.set_ylim(-8, 8)
    e.set_xlabel("radius on the emission surface (au)")
    e.set_ylabel("v$_\\phi$ / v$_{Kepler}$ − 1 (%)")
    e.set_title("Rotation speed vs. Kepler: the gas speeds up and slows down around the dust gaps (blue) and rings (orange)")
    # (e) the numbers
    g = fig.add_subplot(gs[1, 2])
    g.axis("off")
    pub = t
    lines = [
        ("", "ours", "Teague+19"),
        ("M$_*$ (M$_\\odot$, at 101 pc)", f"{f['mstar']:.3f} ± {np.hypot(f['errors']['mstar'], k['mstar_systematic']):.3f}",
         f"{pub['mstar'][0]:.3f} ± {pub['mstar'][1]:.3f}"),
        ("inclination (°)", f"{f['inc']:.2f} ± {f['errors']['inc']:.2f}", f"{pub['inc_deg'][0]} ± {pub['inc_deg'][1]}"),
        ("PA, red side (°)", f"{f['pa']:.2f} ± {f['errors']['pa']:.2f}", f"{pub['pa_deg'][0]} ± {pub['pa_deg'][1]}"),
        ("v$_{sys}$ (km/s)", f"{f['vsys']:.3f} ± {f['errors']['vsys']:.3f}", f"{pub['vlsr_kms'][0]} ± {pub['vlsr_kms'][1]}"),
        ("z$_0$ (mas)", f"{f['z0_mas']:.0f} ± {f['errors']['z0']:.0f}", f"{pub['z0_mas'][0]:.0f} ± {pub['z0_mas'][1]:.0f}"),
        ("flaring φ", f"{f['phi']:.2f} ± {f['errors']['phi']:.2f}", f"{pub['phi'][0]:.2f} ± {pub['phi'][1]:.2f}"),
        ("near side", f"PA {k['near_side_pa']:.0f}° (NE)", "NE*"),
        (f"M$_*$ at {k['mstar_at_gaia']['distance']:.1f} pc (Gaia)", f"{k['mstar_at_gaia']['mstar']:.3f} M$_\\odot$", ""),
    ]
    for i, (lab, a1, a2) in enumerate(lines):
        yy = 0.95 - i * 0.105
        w = "bold" if i == 0 else "normal"
        g.text(0.0, yy, lab, fontsize=8.5, color=ps.INK_SECONDARY, transform=g.transAxes, fontweight=w)
        g.text(0.50, yy, a1, fontsize=8.5, color=ps.C_OURS, transform=g.transAxes, fontweight=w)
        g.text(0.86, yy, a2, fontsize=8.5, color=ps.C_PUBLISHED, transform=g.transAxes, fontweight=w)
    g.text(0.0, -0.06, f"Our mass error combines statistical (MCMC with a jitter term, ±{f['errors']['mstar']:.3f}) and\n"
                       f"systematic (radial range, inclination prior, surface taper: ±{k['mstar_systematic']:.3f}).\n"
                       "Teague et al. quote the statistical error only. M scales with the distance.\n"
                       "*not stated by Teague et al.; DSHARP II puts the far side to the SW.",
           fontsize=7.5, color=ps.MUTED, transform=g.transAxes)
    fig.suptitle("HD 163296 weighed by its own disc: a flared Keplerian disc fitted to the CO velocity field", x=0.02,
                 ha="left")
    fig.text(0.02, 0.925, f"DSHARP 12CO J=2-1 cube (0.32 km/s channels, 0.1\" beam); line centres from a parabola through "
                          f"each spectrum's peak; {f['n_points']} independent pixels (one per beam) between 0.4\" and 3\". "
                          f"Median |residual| {k['resid_mad_ms']:.0f} m/s.", fontsize=8.5, color=ps.INK_SECONDARY, ha="left")
    ps.save(fig, FIGURE_DIR / "hd163296_kinematics", formats=("png",))


def hd163296_channels(res: dict) -> None:
    """Channel maps around Pinte et al.'s two kinks, with the model's isovelocity curves."""
    k = res["kinematics"]
    if "fit" not in k or not CO_CUBE_PATH.exists():
        return
    cube = A.read_cube(CO_CUBE_PATH, "HD163296", 3.6)
    m = _kin_model(k)
    xx, yy = np.meshgrid(cube.x, cube.y)
    vmod, r_s = m.velocity(xx, yy)
    vmod = np.where(r_s < 3.4, vmod, np.nan)
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9.4), facecolor=ps.DARK_PANEL)
    fig.subplots_adjust(top=0.86, bottom=0.04, left=0.03, right=0.97, wspace=0.05, hspace=0.12)
    peak = np.nanpercentile(cube.data, 99.95)
    for row, name in enumerate(("HD163296 #1", "HD163296 #2")):
        kk = k["kinks"][name]
        ic = int(np.argmin(np.abs(cube.velocity - kk["v_planet"])))
        sep, pa_k = kk["sep"], kk["pa"]
        xk, yk = sep * np.sin(np.radians(pa_k)), sep * np.cos(np.radians(pa_k))
        dv_kink = 0.2 * kk["v_kep"] * np.sin(np.radians(m.inc))
        for col, j in enumerate((ic - 1, ic, ic + 1)):
            ax = axes[row, col]
            ax.imshow(np.clip(cube.data[j] / peak, 0, 1) ** 0.6, origin="lower", cmap=ps.IMAGE_CMAP,
                      extent=[cube.x[0], cube.x[-1], cube.y[0], cube.y[-1]], vmin=0, vmax=1)
            vc = cube.velocity[j]
            ax.contour(cube.x, cube.y, vmod, levels=[vc], colors="white", linewidths=1.0)
            ax.contour(cube.x, cube.y, vmod, levels=[vc - dv_kink, vc + dv_kink], colors="white", linewidths=0.6,
                       linestyles=":")
            ax.plot(xk, yk, "o", color="#5ce1e6", ms=5)
            ax.plot(*kk["surface_xy"], "o", mfc="none", mec="#5ce1e6", ms=13, mew=1.0)
            ax.text(0.03, 0.96, f"{vc:.2f} km/s", transform=ax.transAxes, color="white", fontsize=9, va="top")
            ax.set_xlim(3.4, -3.4)
            ax.set_ylim(-3.4, 3.4)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            if col == 0:
                ax.text(0.03, 0.06, f"{name.replace('HD163296 ', 'kink ')} (Pinte et al. 2020): {kk['sep']:.2f}\" at "
                                    f"PA {kk['pa']:.0f}° = {kk['r_mid_au']:.0f} au in our geometry (theirs: {kk['gap_au']:.0f} au)",
                        transform=ax.transAxes, color="#5ce1e6", fontsize=8)
    fig.text(0.03, 0.965, "Looking for a planet's wake in the gas: HD 163296's CO channel maps", color="white", fontsize=14,
             fontweight="bold", va="top")
    fig.text(0.03, 0.93, "Each panel is one 0.32 km/s channel: the gas at that line-of-sight velocity. Solid white: where our "
                         "Keplerian model puts that velocity on the emitting surface; dotted: ± 20% of the local Keplerian "
                         "speed (Pinte et al.'s gauge).\nCyan dot: the candidate planet in the midplane; circle: the same orbit "
                         "position on the CO surface. A planet's wake would bend the emission off the model curve there.",
             color="#c3c2b7", fontsize=8.5, va="top")
    ps.save(fig, FIGURE_DIR / "hd163296_channels", formats=("png",))


# ---------------------------------------------------------------------------
# 9. Spirals
# ---------------------------------------------------------------------------
def spirals_figure(res: dict) -> None:
    import spirals as S

    sp = res["spirals"]
    fig, axes = plt.subplots(3, 3, figsize=(15, 13.5))
    fig.subplots_adjust(top=0.9, bottom=0.07, hspace=0.55, wspace=0.28)
    for row, key in enumerate(("Elias27", "IMLup", "WaOph6")):
        o = res["discs"][key]
        g = o["geometry"]
        d = o["d_dsharp"]
        geom = G.DiscGeometry(g["x0"], g["y0"], g["inc"], g["pa"])
        img = A.load_continuum(key).cutout(o["half_width"], *o["centre_guess"])
        r_hi = HUANG_SPIRAL_RANGES[key][-1][1] / d
        pol = S.polar_image(img, geom, 1.25 * r_hi, img.beam_mean / 4, 1.0)
        s = sp[key]
        # (1) the image with the fitted arms
        a = axes[row, 0]
        xr = (img.x - g["x0"]) * d
        yr = (img.y - g["y0"]) * d
        peak = _peak_of(img, g["x0"], g["y0"], 0.5)
        a.imshow(asinh_stretch(img.data, peak, 0.05), origin="lower", cmap=ps.IMAGE_CMAP,
                 extent=[xr[0], xr[-1], yr[0], yr[-1]], vmin=0, vmax=1)
        for arm, col in zip(sorted(s["arms"]), (ps.AQUA, ps.AQUA, "#ffffff", "#ffffff")):
            ad = s["arms"][arm]
            th = np.linspace(ad["theta_range"][0], ad["theta_range"][1], 300)
            rr = ad["r0_au"] * np.exp(ad["b"] * np.radians(th))
            xs, ys = G.sky_from_polar(rr, th, 0, 0, g["inc"], g["pa"])
            a.plot(xs, ys, color="#5ce1e6" if arm.startswith("arm1") else "white", lw=1.2)
        lim = 1.1 * r_hi * d
        a.set_xlim(lim, -lim)
        a.set_ylim(-lim, lim)
        a.set_xticks([])
        a.set_yticks([])
        a.grid(False)
        a.set_title(f"{o['name']}: our fitted arms", fontsize=10)
        # (2) the residual map in (theta, r) with the traced points
        b = axes[row, 1]
        sel = (pol.r * d >= HUANG_SPIRAL_RANGES[key][0][0] * 0.8) & (pol.r * d <= r_hi * d * 1.05)
        resid = pol.residual[sel]
        vmax = np.nanpercentile(np.abs(resid), 98)
        b.imshow(resid, origin="lower", aspect="auto", cmap=ps.VELOCITY_CMAP, vmin=-vmax, vmax=vmax,
                 extent=[pol.theta[0], pol.theta[-1], pol.r[sel][0] * d, pol.r[sel][-1] * d])
        for arm in sorted(s["arms"]):
            ad = s["arms"][arm]
            th = np.array(ad["theta"])
            b.plot(((th + 180) % 360) - 180, ad["r_au"], ".", color=ps.INK, ms=4)
        b.set_yscale("log")
        b.set_xlabel("azimuth θ in the disc (deg)")
        b.set_ylabel("radius (au)")
        b.set_title("Brightness − azimuthal median; dots: arm crests", fontsize=9.5)
        b.grid(False)
        # (3) the log-spiral Fourier spectrum
        c = axes[row, 2]
        f = s["fourier"]
        p_grid, amp = np.array(f["p_grid"]), np.array(f["amp"])
        c.plot(p_grid, amp, color=ps.C_OURS, lw=1.4)
        c.axvline(f["p_max"], color=ps.C_OURS, lw=0.8, ls=":")
        for nm, pv in s["published"].items():
            pp = -2.0 / pv["b"]
            c.axvline(pp, color=ps.C_PUBLISHED, lw=0.9)
        c.set_xlabel("p (radial wavenumber in ln r)")
        c.set_ylabel("|A(m = 2, p)|")
        arms_txt = ", ".join(f"{k2}: {v['pitch']:.1f}±{v['pitch_err']:.1f}°" for k2, v in sorted(s["arms"].items()))
        pub_txt = ", ".join(f"{k2}: {v['pitch']:.1f}°" for k2, v in s["published"].items())
        c.set_title(f"Fourier transform: pitch {f['pitch']:.1f}°", fontsize=9.5)
        c.text(0.0, -0.27, f"traced arms: {arms_txt}\nDSHARP III: {pub_txt} (violet lines)",
               transform=c.transAxes, fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    fig.suptitle("Spiral arms in the dust: how tightly are they wound?", x=0.02, ha="left")
    fig.text(0.02, 0.935, "Pitch angle = angle between an arm and a circle; a logarithmic spiral r = R$_0$ e$^{bθ}$ has "
                          "tan(pitch) = |b|. Two independent measurements: tracing the arm crests (as DSHARP III), and the "
                          "2-D Fourier transform in (θ, ln r) used for galaxy arms.", fontsize=8.5, color=ps.INK_SECONDARY,
             ha="left")
    ps.save(fig, FIGURE_DIR / "spirals", formats=("png",))


# ---------------------------------------------------------------------------
# 10. HL Tau
# ---------------------------------------------------------------------------
def hltau_figure(res: dict) -> None:
    o = res["discs"]["HLTau"]
    g = o["geometry"]
    d = o["d_gaia"]
    img = A.load_continuum("HLTau").cutout(1.0, g["x0"], g["y0"])
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw=dict(width_ratios=[1, 1.5]))
    fig.subplots_adjust(top=0.84, bottom=0.12, wspace=0.18)
    xr, yr = (img.x - g["x0"]) * d, (img.y - g["y0"]) * d
    peak = _peak_of(img, g["x0"], g["y0"], 0.3)
    a.imshow(asinh_stretch(img.data, peak, 0.08), origin="lower", cmap=ps.IMAGE_CMAP, extent=[xr[0], xr[-1], yr[0], yr[-1]],
             vmin=0, vmax=1)
    a.set_xlim(130, -130)
    a.set_ylim(-130, 130)
    a.set_xlabel("au")
    a.set_ylabel("au")
    a.grid(False)
    a.set_title("HL Tau, 1.0 mm (ALMA 2014 Long Baseline Campaign)", fontsize=10)
    prof = D.load_profile("HLTau")
    r_au = prof.r * o["d_dsharp"]
    pk = np.nanmax(prof.intensity)
    b.fill_between(r_au, (prof.intensity - prof.error) / pk, (prof.intensity + prof.error) / pk, color=ps.GRID, lw=0)
    b.plot(r_au, prof.intensity / pk, color=ps.C_OURS, lw=1.4)
    for f in L.dsharp2_features():
        if f.key == "HLTau" and f.r_au.value:
            b.axvline(f.r_au.value, color=ps.C_PUBLISHED, lw=0.8, alpha=0.8, ls="-" if f.kind == "B" else ":")
    for f in o["features"]:
        b.plot(f["r_au_dsharp"], 1.08, marker="v" if f["kind"] == "D" else "^",
               color=ps.C_GAP if f["kind"] == "D" else ps.C_RING, ms=7 if f["secure"] else 5,
               mfc=None if f["secure"] else "none", ls="none")
    b.set_xlim(0, 130)
    b.set_ylim(0, 1.15)
    b.set_xlabel("radius (au, at 147 pc)")
    b.set_ylabel("intensity / peak (±20° wedges along the major axis)")
    gm = res["geometry"]["HLTau"]
    b.set_title(f"Our profile; violet: ALMA 2015's rings (B solid, D dotted). Our geometry i = {gm['ours']['inc']:.2f}°, "
                f"PA = {gm['ours']['pa']:.2f}° vs {gm['published']['inc']}°, {gm['published']['pa']}°", fontsize=9)
    fig.suptitle("HL Tau (2014): the first picture of a ringed planet-forming disc, re-measured", x=0.02, ha="left")
    ps.save(fig, FIGURE_DIR / "hltau", formats=("png",))


# ---------------------------------------------------------------------------
# 11. How a gap width becomes a planet mass
# ---------------------------------------------------------------------------
def gap_to_planet(res: dict) -> None:
    o = res["discs"]["AS209"]
    prof = D.load_profile("AS209")
    d = o["d_dsharp"]
    pair = next(p for p in o["pairs"] if p["gap"] == "D61")
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.5, 5.2))
    fig.subplots_adjust(top=0.82, bottom=0.13, wspace=0.25)
    r_au = prof.r * d
    a.plot(r_au, prof.intensity * 1e3, color=ps.C_OURS, lw=1.6)
    a.axhline(pair["i_mean"] * 1e3, color=ps.MUTED, lw=0.8)
    a.axvspan(pair["r_gap_in"] * d, pair["r_gap_out"] * d, color=ps.C_GAP, alpha=0.12, lw=0)
    a.axvspan(pair["r_gap_out"] * d, pair["r_ring_out"] * d, color=ps.C_RING, alpha=0.12, lw=0)
    a.plot(pair["gap_r"] * d, pair["i_gap"] * 1e3, "v", color=ps.C_GAP, ms=8)
    a.plot(pair["ring_r"] * d, pair["i_ring"] * 1e3, "^", color=ps.C_RING, ms=8)
    a.annotate("I$_{mean}$ = (I$_{gap}$ + I$_{ring}$)/2", (pair["r_gap_out"] * d, pair["i_mean"] * 1e3), xytext=(30, 25),
               textcoords="offset points", fontsize=8.5, arrowprops=dict(arrowstyle="-", color=ps.MUTED, lw=0.8))
    a.text(0.5 * (pair["r_gap_in"] + pair["r_gap_out"]) * d, 0.02, f"gap\n{pair['gap_width'] * d:.1f} au", ha="center",
           fontsize=8.5, color=ps.C_GAP, transform=a.get_xaxis_transform())
    a.text(0.5 * (pair["r_gap_out"] + pair["r_ring_out"]) * d, 0.02, f"ring\n{pair['ring_width'] * d:.1f} au",
           ha="center", fontsize=8.5, color=ps.C_RING, transform=a.get_xaxis_transform())
    a.set_xlim(40, 95)
    a.set_ylim(0, 0.3)
    a.set_xlabel("radius (au)")
    a.set_ylabel("intensity (mJy/beam)")
    a.set_title(f"AS 209's gap D61: width Δ = (r$_{{out}}$ − r$_{{in}}$)/r$_{{out}}$ = {pair['delta']:.3f}", fontsize=10)
    # the relation, for a 0.8 Msun star with h/r = 0.07
    deltas = np.linspace(0.1, 0.7, 200)
    for alpha, col in ((1e-4, ps.BLUE), (1e-3, ps.INK), (1e-2, ps.ORANGE)):
        ms = [P.planet_mass_mjup(dd, 0.07, alpha, 3.0, "DSD1", 1.0)[0] for dd in deltas]
        b.plot(deltas, ms, color=col, lw=1.6, label=f"α = {alpha:g}")
    pl = [r for r in res["planets"]["rows"] if r.get("found")]
    for r in pl:
        m = r["masses"]["DSD1_0.001"]["m"] / r["mstar"]  # per solar mass, for this illustration
        b.plot(r["delta"], m, "o", color=ps.C_OURS, ms=4, alpha=0.6)
    b.set_yscale("log")
    b.set_xlabel("gap width Δ")
    b.set_ylabel("planet mass (M$_J$ per M$_\\odot$ of star)")
    b.set_title("DSHARP VII's simulations: wider gaps need heavier planets (h/r = 0.07, St = 1.6×10$^{-3}$)",
                fontsize=10)
    b.legend(loc="upper left")
    fig.suptitle("From a gap to a planet: measure the width, ask the simulations what could carve it", x=0.02, ha="left")
    fig.text(0.02, 0.9, "Dots: the 11 DSHARP gaps we measured (each with its own disc temperature and grain size; "
                        "shown at their α = 10$^{-3}$ masses). More viscous discs refill gaps faster, so the same gap needs "
                        "a heavier planet.", fontsize=8.5, color=ps.INK_SECONDARY, ha="left")
    ps.save(fig, FIGURE_DIR / "gap_to_planet", formats=("png",))


if __name__ == "__main__":
    FIGURE_DIR.mkdir(exist_ok=True)
    r = load_results()
    gallery_same_scale(r)
    gallery_annotated(r)
    profiles_grid(r)
    geometry_comparison(r)
    features_comparison(r)
    planets_mass_separation(r)
    planet_mass_check(r)
    snowlines(r)
    hd163296_kinematics(r)
    hd163296_channels(r)
    spirals_figure(r)
    hltau_figure(r)
    gap_to_planet(r)
