"""Static figures (PNG + SVG) for the light-bending project.

    python figures.py            # every figure, from data/results.json (run analysis.py first)
    python figures.py 01 02      # only figures whose number starts with these
"""

from __future__ import annotations

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle

import config as C
import plotstyle as P

P.use_style()

HALO = [pe.withStroke(linewidth=2.6, foreground=P.SURFACE)]
DARK_HALO = [pe.withStroke(linewidth=2.4, foreground="#000000")]


def _out(stem: str) -> str:
    return str(C.FIGURE_DIR / stem)


def _results() -> dict:
    return json.loads(C.RESULTS_PATH.read_text())


# ---------------------------------------------------------------------------
# 01 The star field of 1919 on the real sky
# ---------------------------------------------------------------------------

#: The DSS2 colour image fetched by fetch_data.py: gnomonic (TAN), centred on
#: RA 66.05, Dec +21.85 (ICRS), 2100 x 1800 px of 0.002 deg, north up, east left
#: (checked: Gaia's 13 brightest stars land within 1.5 px).
SKY_CENTRE = (66.05, 21.85)
SKY_HALF_W, SKY_HALF_H = 2.1 * 60, 1.8 * 60          # arcmin


def _sky_arcmin(ra_deg, dec_deg):
    import stars as S

    xi, eta = S.standard_coords(np.deg2rad(np.asarray(ra_deg)), np.deg2rad(np.asarray(dec_deg)),
                                np.deg2rad(SKY_CENTRE[0]), np.deg2rad(SKY_CENTRE[1]))
    return np.rad2deg(xi) * 60, np.rad2deg(eta) * 60


#: label offsets (arcmin, east/north) for each Table I star in figure 01
FIELD_LABEL = {1: (-7, -8), 2: (-1, -9), 3: (-9, -5), 4: (-9, 5), 5: (-1, -9), 6: (-8, 5), 7: (-9, 4),
               8: (-9, -4), 9: (-9, -6), 10: (-8, 5), 11: (2, -9), 12: (2, 8), 13: (-9, 2)}


def fig_01_field(res: dict) -> list:
    """The 13 stars of Table I on a real image of the sky, the Sun hidden behind the
    Moon at Sobral's mid-totality, and Einstein's displacement of each star, enlarged."""
    from PIL import Image

    import stars as S

    img = np.asarray(Image.open(C.SKY_DIR / "hyades_dss2_color.jpg"))
    t1 = S.table1().set_index("no")
    g = S.gaia_table1()
    pos = S.gaia_at(g, S.MID_TOTALITY["sobral"])
    sx, sy = _sky_arcmin(pos.ra.deg, pos.dec.deg)
    sun = S.sun_at("sobral", S.MID_TOTALITY["sobral"])
    sunx, suny = _sky_arcmin(sun["ra_deg"], sun["dec_deg"])
    used_4in = [11, 5, 4, 3, 6, 10, 2]
    mag = 600.0          # arrow scale: 1" of deflection drawn as 600" = 10'

    fig, ax = plt.subplots(figsize=(9.6, 8.6))
    ax.imshow(img, extent=[SKY_HALF_W, -SKY_HALF_W, -SKY_HALF_H, SKY_HALF_H], origin="upper", interpolation="lanczos")
    ax.set_xlim(SKY_HALF_W, -SKY_HALF_W)
    ax.set_ylim(-SKY_HALF_H, SKY_HALF_H)
    ax.grid(False)
    # the Sun (behind the Moon) at mid-totality, Sobral
    ax.add_patch(Circle((float(sunx), float(suny)), sun["radius_arcmin"], facecolor="#000000",
                        edgecolor="#f2efe6", lw=1.2, zorder=3))
    ax.text(float(sunx), float(suny), "the Sun,\nhidden by\nthe Moon", ha="center", va="center", fontsize=7.5,
            color="#e9e6dc", zorder=4)
    for no in t1.index:
        x, y = float(sx[g.index.get_loc(no)]), float(sy[g.index.get_loc(no)])
        on4 = no in used_4in
        col = P.BLUE_300 if on4 else "#d7d4ca"
        ax.add_patch(Circle((x, y), 3.2, fill=False, ec=col, lw=1.4 if on4 else 0.9, zorder=5))
        dx, dy = t1.loc[no, "sobral_dx"] * mag / 60, t1.loc[no, "sobral_dy"] * mag / 60
        r0 = 3.4 / np.hypot(dx, dy)
        ax.annotate("", xy=(x + dx * (1 + r0), y + dy * (1 + r0)), xytext=(x + dx * r0, y + dy * r0),
                    arrowprops={"arrowstyle": "-|>", "color": "#ffd166", "lw": 1.3, "mutation_scale": 9}, zorder=6)
        ox, oy = FIELD_LABEL[no]
        ax.text(x + ox, y + oy, str(no), color=col, fontsize=9, fontweight="bold", ha="center", va="center",
                path_effects=DARK_HALO, zorder=7)
    # legend-like key, top-left on the dark sky
    ax.text(SKY_HALF_W - 4, SKY_HALF_H - 5,
            "○ blue: the 7 stars measured on the Sobral 4-inch plates\n○ grey: the other stars of the paper's Table I\n"
            "→ yellow: Einstein's outward shift of each star, drawn 600× longer\n    (1″ of deflection = 10′ on this map)",
            ha="left", va="top", fontsize=7.6, color="#e9e6dc", linespacing=1.45, zorder=8)
    ax.text(-SKY_HALF_W + 4, -SKY_HALF_H + 5, "N ↑   E ←", ha="right", va="bottom", fontsize=8, color="#e9e6dc")
    # a 30' scale bar, bottom-left
    x0, y0 = SKY_HALF_W - 6, -SKY_HALF_H + 8
    ax.plot([x0, x0 - 30], [y0, y0], color="#e9e6dc", lw=1.4)
    ax.text(x0, y0 + 2.5, "30′ (the Sun is 31.5′ across)", ha="left", va="bottom", fontsize=7.2, color="#e9e6dc")
    ax.set_xlabel("arcminutes east of the image centre (RA 4h 24.2m, Dec +21° 51′)")
    ax.set_ylabel("arcminutes north")
    ax.set_title("29 May 1919: the stars of the Hyades around the eclipsed Sun")
    # no-break spaces keep each "number name" entry on one caption line
    names = "; ".join(f"{n} {C.TABLE1_GAIA[n][1]}".replace(" ", "\u00a0") for n in t1.index)
    P.caption(fig, "Sky: DSS2 colour (POSS-II © Caltech, digitised by STScI; CDS hips2fits). Stars: Gaia DR3 "
                   "positions carried back to 29 May 1919; numbers "
                   f"are the paper's Table I ({names}; the paper's '66 Tauri' is 56 Tau). Sun: JPL Horizons, seen from Sobral "
                   "at mid-totality (12:01 UT). Arrows: the displacements printed in Table I (1.75″ at the limb).")
    return P.save(fig, _out("01_field_1919"))


# ---------------------------------------------------------------------------
# 02 Plate 1 of the paper, identified and dated
# ---------------------------------------------------------------------------

PLATE_CROP = (230, 310, 2370, 3060)       # x0, y0, x1, y1 in leaf 42 (full resolution)
PLATE_LABEL = {10: (-300, 110), 6: (-300, 60), 4: (230, -130), 3: (260, 110), 2: (60, -150), 5: (230, -150), 11: (-230, -150)}


def fig_02_plate1(res: dict) -> list:
    from PIL import Image

    import plate1 as P1

    r = res["plate1"]
    im = np.asarray(Image.open(C.PAGES_DIR / f"leaf{P1.PLATE1_LEAF:03d}.jpg").convert("L"))
    x0, y0, x1, y1 = PLATE_CROP
    fig = plt.figure(figsize=(10.4, 8.4))
    ax = fig.add_axes([0.02, 0.10, 0.58, 0.82])
    ax.imshow(im[y0:y1, x0:x1], cmap="gray", vmin=0, vmax=255, extent=[x0, x1, y1, y0], interpolation="lanczos")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_axis_off()
    for no, st in r["stars"].items():
        no = int(no)
        x, y = st["x_px"], st["y_px"]
        ox, oy = PLATE_LABEL[no]
        ax.add_patch(Circle((x, y), 34, fill=False, ec=P.BLUE, lw=1.3))
        ax.annotate(f"{no} {C.TABLE1_GAIA[no][1]}", (x, y), (x + ox, y + oy), fontsize=8.2, color=P.BLUE,
                    ha="center", va="center", fontweight="bold", path_effects=HALO,
                    arrowprops={"arrowstyle": "-", "color": P.BLUE, "lw": 0.9, "shrinkA": 3, "shrinkB": 8})
    m = r["moon"]
    ax.add_patch(Circle((m["cx"], m["cy"]), m["r_px"], fill=False, ec=P.ORANGE, lw=1.1, ls=(0, (4, 3))))
    k = r["scale_px_per_arcmin"]
    xb, yb = x0 + 90, y1 - 90
    ax.plot([xb, xb + 10 * k], [yb, yb], color=P.INK, lw=1.6)
    ax.text(xb + 5 * k, yb - 25, "10′", ha="center", va="bottom", fontsize=8, color=P.INK, path_effects=HALO)
    # compass from the fitted affine map: the print is not mirrored, but it is turned ~42 degrees
    M = np.array(r["affine"])[:2].T                      # (east, north) arcmin -> (right, down) px
    ox_, oy_, L = x1 - 170, y0 + 170, 110
    for vec, lab in ((M @ [0.0, 1.0], "N"), (M @ [1.0, 0.0], "E")):
        u = vec / np.hypot(*vec)
        ax.annotate("", (ox_ + L * u[0], oy_ + L * u[1]), (ox_, oy_),
                    arrowprops={"arrowstyle": "-|>", "color": P.INK, "lw": 1.1, "mutation_scale": 9})
        ax.text(ox_ + 1.35 * L * u[0], oy_ + 1.35 * L * u[1], lab, ha="center", va="center", fontsize=8.5,
                color=P.INK, path_effects=HALO)
    nx, ny = M @ [0.0, 1.0]
    turn = np.degrees(np.arctan2(nx, -ny))               # north, clockwise from the top of the page
    ax.set_title("Plate 1 of the 1920 paper: its stars identified, its moment found", loc="left", x=0.0)

    # right: when was the Moon where the plate shows it?
    bx = fig.add_axes([0.67, 0.50, 0.31, 0.38])
    import plate1 as P1b

    tr = P1b.moon_sun_track()
    sob = tr[tr["site"] == "sobral"].reset_index(drop=True)
    sol = np.array(r["affine"])
    mx, my = P1b.tangent_arcmin(np.deg2rad(sob["moon_ra"]), np.deg2rad(sob["moon_dec"]), *r["centre_rad"])
    mp = np.column_stack([mx, my, np.ones_like(mx)]) @ sol
    dist = np.hypot(mp[:, 0] - m["cx"], mp[:, 1] - m["cy"]) / k * 60          # arcsec
    bx.plot(sob["t"], dist, color=P.ORANGE, lw=1.6)
    t_star = pd.Timestamp("1919-05-29 " + r["dated_utc"])
    bx.axvline(t_star, color=P.INK, lw=0.9, ymax=0.9)          # stops below the exposure labels
    c2 = pd.Timestamp("1919-05-29 " + res["contacts"]["sobral"]["second_contact_utc"])
    c3 = pd.Timestamp("1919-05-29 " + res["contacts"]["sobral"]["third_contact_utc"])
    bx.axvspan(c2, c3, color=P.WASH, lw=0, zorder=0)
    ex = _sobral_4inch_exposures()
    if ex is not None:
        for _, e in ex.iterrows():
            bx.axvspan(e["start"], e["end"], color=P.BLUE_150, lw=0, zorder=0.5)
            bx.text(e["mid"], 162, str(e["ref"]), ha="center", va="bottom", fontsize=6.8, color=P.BLUE_550)
    bx.set_ylim(0, 175)
    bx.set_ylabel("Moon's predicted centre minus\nits centre on the plate (″)")
    bx.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    bx.set_xlim(pd.Timestamp("1919-05-29 11:57:30"), pd.Timestamp("1919-05-29 12:04:00"))
    bx.set_xlabel("29 May 1919, UT")
    em = r.get("exposure_match") or {}
    t_lab = pd.Timestamp("1919-05-29 " + r["dated_utc"]).round("s").strftime("%H:%M:%S")
    lab = f"best match:\n{t_lab} UT ± {r['sigma_t_s']:.0f} s"
    if em.get("plate"):
        lab += f"\n= 4-inch plate {em['plate']}"
    P.annotate(bx, lab, (t_star, 2), (t_star - pd.Timedelta(133, unit="s"), 118), va="center", ha="left")
    txt = (f"The seven marked stars match Gaia DR3 (carried back to 1919)\nto {r['rms_px'] / k * 60:.1f}″ rms; any other "
           f"assignment is ≥ {r['next_best_rms_px'] / k * 60:.0f}″ off.\n\nThe white disc is the Moon: no light reached the "
           f"negative there.\nIts edge (dashed orange) gives a radius of {m['r_arcmin']:.2f}′;\nJPL's Moon seen from Sobral: "
           f"{m['r_arcmin_horizons']:.2f}′.\n\nThe Moon slid across the Sun at 0.4′ a minute, so its\ncentre dates the exposure. "
           "Shaded: totality.\nBlue: the eight 4-inch exposures (p. 300), numbered as in\nthe paper's tables "
           "(VI, taken through cloud, shows no stars).")
    fig.text(0.67, 0.42, txt, ha="left", va="top", fontsize=7.8, color=P.INK_SECONDARY, linespacing=1.35)
    P.caption(fig, "Plate 1 of Dyson, Eddington & Davidson (1920), a half-tone of one of the Sobral 4-inch negatives "
                   "(Royal Society scan, Internet Archive philtrans06337895, Public Domain Mark). Star marks found by "
                   "machine; Moon and Sun: JPL Horizons, topocentric for Sobral. Compass: from the fit of the stars' "
                   f"positions; the print is not mirrored, but turned {turn:.0f}° clockwise from north-up.")
    return P.save(fig, _out("02_plate1"))


def _sobral_4inch_exposures() -> pd.DataFrame | None:
    """Start, middle and end (UT) of each 4-inch exposure on 29 May 1919 (p.300),
    if the exposure table has been transcribed; None otherwise."""
    path = C.TRANSCRIBED_DIR / "dyson1920_sobral_exposures.csv"
    if not path.exists():
        return None
    import plate1 as P1

    return P1.sobral_exposures(path, instrument="4-inch")


# ---------------------------------------------------------------------------
# 03 The stars moved: observed shifts against Einstein and Newton
# ---------------------------------------------------------------------------

def _read_tr(name: str) -> pd.DataFrame:
    return pd.read_csv(C.TRANSCRIBED_DIR / name, comment="#", skipinitialspace=True)


def _star_shifts() -> dict:
    """Observed and predicted displacements (arcsec) of each star, by station.
    Sobral 4-inch: p.309 (eclipse mean minus comparison mean, Table VIII), with the
    star-6 R.A. sign misprint corrected (+0.10 printed, -0.10 used everywhere else).
    Principe: Einstein's full displacement (Table I) plus the mean of the p.328
    residuals over the comparisons in which the star was measured."""
    import stars as S

    t1 = S.table1().set_index("no")
    p309 = _read_tr("dyson1920_p309_obs_calc.csv").set_index("star")
    p309.loc[6, "ra_observed"] = -0.10
    sob = pd.DataFrame({"x": t1.loc[p309.index, "x"], "y": t1.loc[p309.index, "y"],
                        "obs_dx": p309["ra_observed"], "obs_dy": p309["dec_observed"],
                        "gr_dx": p309["ra_calculated"], "gr_dy": p309["dec_calculated"]})
    # one-sigma of each star's radial shift: the scatter of its 7 eclipse and 7 comparison
    # plate values in Table VIII, projected away from the Sun (JPL's position)
    t8 = _read_tr("dyson1920_table8_residuals.csv")
    t8 = t8[t8["plate"] != "Mean"]
    sx, sy = json.loads(C.RESULTS_PATH.read_text())["frame"]["sites"]["sobral"]["sun_xy_horizons"]
    err = {}
    for star, g in t8.groupby("star"):
        ux, uy = t1.loc[star, "x"] - sx, t1.loc[star, "y"] - sy
        var = 0.0
        for _, gk in g.groupby("kind"):
            w = gk.pivot(index="plate", columns="coord", values="value_arcsec")
            rad = (w["RA"] * ux + w["Dec"] * uy) / np.hypot(ux, uy)
            var += rad.var(ddof=1) / len(rad)
        err[star] = float(np.sqrt(var))
    sob["err"] = pd.Series(err)
    r328 = _read_tr("dyson1920_p328_principe_residuals.csv").set_index("star")
    rx = r328[["x_G", "x_H", "x_D", "x_I"]].mean(axis=1)
    ry = r328[["y_G", "y_H", "y_D", "y_I"]].mean(axis=1)
    pr = pd.DataFrame({"x": t1.loc[r328.index, "x"], "y": t1.loc[r328.index, "y"],
                       "gr_dx": t1.loc[r328.index, "principe_dx"], "gr_dy": t1.loc[r328.index, "principe_dy"]})
    pr["obs_dx"], pr["obs_dy"] = pr["gr_dx"] + rx, pr["gr_dy"] + ry
    pr["n"] = r328[["y_G", "y_H", "y_D", "y_I"]].notna().sum(axis=1)
    return {"sobral": sob, "principe": pr}


def fig_03_star_shifts(res: dict) -> list:
    sh = _star_shifts()
    fr = res["frame"]["sites"]
    k = 0.50                      # Table I units (50') per arcsec of displacement: 1" drawn as 25'

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 6.1), gridspec_kw={"wspace": 0.08})
    for ax, (site, col, title, note) in zip(axes, (
            ("sobral", P.BLUE, "Sobral, 4-inch lens: 7 plates", "eclipse mean minus comparison mean (p.309)"),
            ("principe", P.ORANGE, "Príncipe: 2 plates, 4 comparisons", "Einstein + the p.328 residuals"))):
        d = sh[site]
        sx, sy = fr[site]["sun_xy_horizons"]
        rad = fr[site]["sun_radius_arcmin"] / C.UNIT_TABLE1_ARCMIN
        ax.add_patch(Circle((sx, sy), rad, facecolor=P.INK, edgecolor=P.INK, lw=0.8, zorder=2))
        for no, r in d.iterrows():
            ax.plot(r["x"], r["y"], marker="*", ms=9, color=P.INK_SECONDARY, mec=P.SURFACE, mew=0.5, zorder=3)
            gx, gy = r["x"] + k * r["gr_dx"], r["y"] + k * r["gr_dy"]
            nx, ny = r["x"] + k * r["gr_dx"] / 2, r["y"] + k * r["gr_dy"] / 2
            ox, oy = r["x"] + k * r["obs_dx"], r["y"] + k * r["obs_dy"]
            ax.annotate("", xy=(ox, oy), xytext=(r["x"], r["y"]),
                        arrowprops={"arrowstyle": "-|>", "mutation_scale": 10, "shrinkA": 0, "shrinkB": 0,
                                    "color": col, "lw": 2.2}, zorder=4)
            s1 = r["err"] if site == "sobral" else 0.22 / 0.6745 / np.sqrt(r["n"])     # Principe: as figure 04
            ax.add_patch(Circle((ox, oy), k * s1, fill=False, ec=col, lw=0.8, ls=(0, (2, 1.5)), zorder=4))
            ax.plot(gx, gy, "o", ms=7.5, mfc="none", mec=P.INK, mew=1.4, zorder=5)
            ax.plot(nx, ny, "+", ms=7, color=P.INK_SECONDARY, mew=1.3, zorder=5)
            # label on the side of the star away from its arrow; the close pair 3 and 4 sideways
            ux, uy = r["obs_dx"], r["obs_dy"]
            n = np.hypot(ux, uy) or 1.0
            lx, ly = {3: (0.09, 0.0), 4: (-0.09, 0.0)}.get(no, (-0.085 * ux / n, -0.085 * uy / n))
            ax.text(r["x"] + lx, r["y"] + ly, str(no), fontsize=8, color=P.INK_SECONDARY,
                    ha="center", va="center", path_effects=HALO, zorder=6)
        ax.set_aspect("equal")
        ax.set_xlim(1.78, -1.6)
        ax.set_ylim(-1.45, 1.75)
        ax.set_title(title, fontsize=10)
        ax.text(0.98, 0.98, f"observed: {note}", transform=ax.transAxes, fontsize=7.3, color=P.INK_SECONDARY,
                va="top", ha="right")
        ax.set_xlabel("east ← standard coordinate x (units of 50′)")
        ax.tick_params(labelsize=7.5)
    axes[0].set_ylabel("y, north (units of 50′)")
    axes[1].tick_params(labelleft=False)
    from matplotlib.lines import Line2D

    hs = [Line2D([], [], color=P.BLUE, lw=2.2, label="observed shift, Sobral 4-inch (dashed circle: 1σ)"),
          Line2D([], [], color=P.ORANGE, lw=2.2, label="observed shift, Príncipe (dashed circle: 1σ)"),
          Line2D([], [], ls="", marker="o", ms=7.5, mfc="none", mec=P.INK, mew=1.4, label="where Einstein puts it"),
          Line2D([], [], ls="", marker="+", ms=7, color=P.INK_SECONDARY, mew=1.3, label="where Newton puts it")]
    fig.legend(handles=hs, loc="upper center", ncol=4, fontsize=7.6, bbox_to_anchor=(0.5, 0.955), frameon=False)
    fig.suptitle("The stars moved outward, as Einstein said, by twice Newton's amount", y=1.0, fontsize=11.5,
                 fontweight="bold")
    P.caption(fig, "Each arrow is a star's displacement, drawn 1,500 times longer than on the sky (1″ = 25′), from its "
                   "1919 position. Black disc: the Sun (JPL). Sobral: the shifts after the plates' adopted scale and "
                   "orientation were removed (p.309; star 6's R.A. corrected from the printed +0.10″ to −0.10″, as Table "
                   "VIII and p.331 require). Príncipe: shifts implied by the residuals of the four comparisons (p.328). "
                   "Dashed circles, one sigma of the outward shift: at Sobral from each star's scatter over its 7 eclipse "
                   "and 7 comparison plates (Table VIII); at Príncipe the check plates' probable error of one measure "
                   "(0.22″, p.318), converted to one sigma (0.33″) and divided by √n, n being the number of comparisons "
                   "that measured the star (as if they were independent).")
    return P.save(fig, _out("03_star_shifts"))


# ---------------------------------------------------------------------------
# 04 The 1/r law: Diagram 2 of 1920, redrawn with both expeditions
# ---------------------------------------------------------------------------

def fig_04_radial_law(res: dict) -> list:
    sh = _star_shifts()
    fr = res["frame"]["sites"]
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    xmax = 0.62
    xs = np.linspace(0, xmax, 50)
    ax.plot(xs, C.DEFLECTION_LIMB_1919 * xs, **P.GR_LINE, label="Einstein: 1.75″ × (R☉/r)")
    ax.plot(xs, C.DEFLECTION_LIMB_1919 / 2 * xs, **P.NEWTON_LINE, label="Newtonian half")
    # Sobral 4-inch: p.331 radial observations; r from Table I and JPL's Sun
    p331 = _read_tr("dyson1920_p331_radial.csv").set_index("star")
    for site, col, mk, lab in (("sobral", P.BLUE, "o", "Sobral 4-inch, 1919 (p.331)"),
                               ("principe", P.ORANGE, "s", "Príncipe, 1919 (p.328)")):
        d = sh[site]
        sx, sy = fr[site]["sun_xy_horizons"]
        R = fr[site]["sun_radius_arcmin"] / C.UNIT_TABLE1_ARCMIN
        ux, uy = d["x"] - sx, d["y"] - sy
        r = np.hypot(ux, uy)
        if site == "sobral":
            radial = p331.loc[d.index, "observation_arcsec"].to_numpy(float)
            err = d["err"].to_numpy(float)
        else:
            radial = ((d["obs_dx"] * ux + d["obs_dy"] * uy) / r).to_numpy(float)
            # p.318: the check plates' p.e. of one measure, 0.22" (p.328 gives 0.21" for y on the
            # eclipse plates), as one sigma and over sqrt(n) comparisons, as if they were independent
            err = 0.22 / 0.6745 / np.sqrt(d["n"].to_numpy(float))
        ax.errorbar(R / r, radial, yerr=err, fmt=mk, color=col, mec=P.SURFACE, mew=0.6, ms=6.5, elinewidth=1.0,
                    capsize=0, zorder=4, label=lab)
    w = res.get("wallal")
    if w and w.get("stars"):
        ws = pd.DataFrame(w["stars"])
        ax.plot(ws["inv_r"], ws["dr"], ls="", marker="D", ms=2.4, color=P.AQUA, alpha=0.35, zorder=2,
                label=f"Wallal 1922, five-foot cameras:\n{len(ws)} stars (scale term added back)")
        edges = np.array([0.5, 0.9, 1.4, 1.9, 2.6, 3.5, 4.6, 5.9, 7.4, 9.5])        # degrees
        wt = 1 / ws["err"] ** 2
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (ws["r_deg"] >= lo) & (ws["r_deg"] < hi)
            if m.sum() < 3:
                continue
            mu = float(np.sum(wt[m] * ws["dr"][m]) / wt[m].sum())
            e = float(1 / np.sqrt(wt[m].sum()))
            x = float(np.sum(wt[m] * ws["inv_r"][m]) / wt[m].sum())
            ax.errorbar(x, mu, yerr=e, fmt="D", ms=6.5, color=P.AQUA, mec=P.SURFACE, mew=0.8, elinewidth=1.2,
                        capsize=0, zorder=5, label="… weighted means in distance bins" if lo == edges[0] else None)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.7, 1.4)                 # every Wallal star (-0.66" to +1.36") is on the plot
    ax.set_xlabel("R☉ / r: the Sun's radius divided by the star's distance from the Sun's centre")
    ax.set_ylabel("outward (radial) displacement of the star (″)")
    top = ax.secondary_xaxis("top", functions=(lambda v: v, lambda v: v))
    ticks_r = [10, 5, 3, 2]
    top.set_xticks([1 / t for t in ticks_r])
    top.set_xticklabels([f"{t} R☉" for t in ticks_r])
    top.tick_params(labelsize=7.5, colors=P.MUTED)
    # narrow enough to sit right of the one Wallal star in this corner (Rsun/r = 0.35, -0.47")
    ax.legend(loc="lower right", fontsize=7.6, frameon=True, facecolor=P.SURFACE, edgecolor="none", framealpha=0.85)
    ax.set_title("Twice Newton, falling off as 1/r: Diagram 2 of 1920, redrawn", pad=22)
    P.caption(fig, "Straight lines through the origin are the 1/r law: a star twice as far from the Sun is shifted "
                   "half as much. Distances from JPL's Sun at each station's mid-totality. Sobral error bars: each "
                   "star's scatter over its 7 eclipse and 7 comparison plates (Table VIII); Príncipe: the check plates' "
                   "probable error of one measure (0.22″, p.318), converted to one sigma (0.33″) and divided by √n, n "
                   "being the number of comparisons that measured the star (as if they were independent). Wallal: each "
                   "star's mean radial displacement (Table 6 of 1928) plus the plates' mean scale term (0.0098″ per "
                   "degree), so that Einstein's law is again a line through the origin.")
    return P.save(fig, _out("04_radial_law"))


# ---------------------------------------------------------------------------
# 05 Every determination, published and ours: let the reader judge
# ---------------------------------------------------------------------------

PE = 0.6745


def _forest_rows(res: dict) -> list:
    """(group, label, value, sigma or None, kind) for the summary figure; kind is
    'published' (as printed, probable errors converted to one sigma), 'harvey' or 'ours'."""
    rows = []
    s4 = res["sobral_4inch"]
    rows += [("Sobral 4-inch\n7 eclipse plates, 7 stars", "1920 paper", 1.98, 0.12 / PE, "published"),
             ("Sobral 4-inch\n7 eclipse plates, 7 stars", "1920 arithmetic, redone", s4["replication"]["limb"], None, "ours"),
             ("Sobral 4-inch\n7 eclipse plates, 7 stars", "least squares, RA + Dec (jackknife error)",
              s4["joint_fit"]["both"]["limb"], max(s4["joint_fit"]["both"]["limb_err"], s4["jackknife"]["err_plates"],
                                                   s4["jackknife"]["err_stars"]), "ours"),
             ("Sobral 4-inch\n7 eclipse plates, 7 stars", "… RA only", s4["joint_fit"]["RA"]["limb"], s4["joint_fit"]["RA"]["limb_err"], "ours"),
             ("Sobral 4-inch\n7 eclipse plates, 7 stars", "… Dec only", s4["joint_fit"]["Dec"]["limb"], s4["joint_fit"]["Dec"]["limb_err"], "ours"),
             ("Sobral 4-inch\n7 eclipse plates, 7 stars", "remeasured 1979 (5 plates)", 1.90, 0.11, "harvey")]
    a = res["astrographic"]
    fa = a["fits"]["scale_common_free"]
    g = "Sobral astrographic\n16 plates, Dec only"
    rows += [(g, "1920 paper (plates' own scale)", 0.93, None, "published"),
             (g, "1920 paper, refraction-only scale", 1.52, None, "published"),
             # the error rule of analysis.step_combined_1919: the larger of formal and plate jackknife
             (g, "least squares, one scale for all (jackknife error)", fa["limb"],
              max(fa["limb_err"], fa.get("limb_err_jackknife", 0.0)), "ours"),
             (g, "… scale fixed by refraction", a["fits"]["scale_fixed_refraction"]["limb"],
              a["fits"]["scale_fixed_refraction"]["limb_err"], "ours"),
             (g, "remeasured 1979 (9 plates reduced, RA + Dec)", 1.55, 0.34, "harvey")]
    if "principe" in res:
        pr = res["principe"]
        g = "Príncipe\n2 eclipse plates, 5 stars each"
        rows += [(g, "1920 paper", 1.61, 0.30 / PE, "published")]
        if pr.get("replication"):
            rows += [(g, "1920 arithmetic, redone", pr["replication"]["limb"], None, "ours")]
        if pr.get("fit"):
            rows += [(g, "least squares, check-plate scale", pr["fit"]["limb"], pr["fit"]["limb_err"], "ours"),
                     (g, "… scale left free", pr["fit_free_scale"]["limb"], pr["fit_free_scale"]["limb_err"], "ours")]
    comb = res.get("combined_1919")
    if comb:
        g = "1919, instruments combined\n(our least-squares fits)"
        c = comb["selection_1920"]
        rows.append((g, "the 1920 choice: 4-inch + Príncipe", c["limb"], c["limb_err"], "ours"))
        c = comb["all"]
        k = c["limb_err_scaled"] / c["limb_err"]
        rows.append((g, f"all three (they disagree, χ² = {c['chi2']:.0f} for {c['dof']};\nerror × {k:.1f})",
                     c["limb"], c["limb_err_scaled"], "ours"))
    if "wallal" in res:
        w = res["wallal"]
        # two subtitle lines, so that the text ends left of the Newtonian line
        g = "Wallal 1922 (Lick)\n15-foot: 92 stars\nfive-foot: 145 stars"
        rows += [(g, "1923: 15-foot cameras", 1.72, 0.11 / PE, "published"),
                 (g, "our fit, 15-foot stars (1923 Table II)", w["fit_15ft"]["D1"]["beta"], w["fit_15ft"]["D1"]["beta_err"], "ours"),
                 (g, "1928: five-foot cameras", 1.82, 0.15 / PE, "published"),
                 (g, "our fit, five-foot cameras", w["fit"]["limb"], w["fit"]["limb_err"], "ours"),
                 (g, "1928: both, combined", 1.75, 0.09 / PE, "published"),
                 (g, "our fits, combined", w["combined"]["beta"], w["combined"]["beta_err"], "ours")]
    return rows


def fig_05_forest(res: dict) -> list:
    rows = _forest_rows(res)
    groups = list(dict.fromkeys(r[0] for r in rows))
    n = len(rows) + 1.6 * len(groups)
    fig, ax = plt.subplots(figsize=(8.8, 0.31 * n + 1.9))
    gr = C.DEFLECTION_LIMB_1919
    ax.axvline(gr, **P.GR_LINE, zorder=1)
    ax.axvline(gr / 2, **P.NEWTON_LINE, zorder=1)
    style = {"published": dict(marker="s", color=P.INK_SECONDARY, mfc=P.SURFACE),
             "harvey": dict(marker="D", color=P.ORANGE, mfc=P.ORANGE),
             "ours": dict(marker="o", color=P.BLUE, mfc=P.BLUE)}
    y = 0.0
    yt, yl = [], []
    for gi, gname in enumerate(groups):
        head, _, sub = gname.partition("\n")
        y -= 0.6 if gi else 0.0
        ax.text(-0.015, y, head, transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=8.6, fontweight="bold", color=P.INK)
        ax.text(0.32, y, sub, ha="left", va="center", fontsize=7.2, color=P.MUTED, path_effects=HALO)
        y -= 1.0
        for g, lab, v, s, kind in [r for r in rows if r[0] == gname]:
            st = style[kind]
            if s is not None:
                ax.errorbar(v, y, xerr=s, fmt="none", ecolor=st["color"], elinewidth=1.3, capsize=0, zorder=3)
            ax.plot(v, y, ls="", marker=st["marker"], ms=6.5, color=st["color"], mfc=st["mfc"], mec=st["color"],
                    mew=1.3, zorder=4)
            txt = f"{v:.2f}" + (f" ± {s:.2f}" if s is not None else "")
            ax.text(2.72, y, txt, ha="right", va="center", fontsize=7.6, color=st["color"])
            yt.append(y)
            yl.append(lab)
            y -= 1.0
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=7.8)
    ax.tick_params(axis="y", length=0, pad=6)
    ax.set_xlim(0.3, 2.75)
    ax.set_ylim(y + 0.4, 0.9)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("deflection of starlight at the Sun's limb (″), with one-sigma errors")
    ax.text(gr + 0.02, 0.75, "Einstein 1.75″", ha="left", va="top", fontsize=8, color=P.INK)
    ax.text(gr / 2 + 0.02, 0.75, "Newtonian 0.87″", ha="left", va="top", fontsize=8, color=P.INK_SECONDARY)
    from matplotlib.lines import Line2D

    hs = [Line2D([], [], ls="", **{k: v for k, v in style[k2].items()}, ms=6.5, mew=1.3, label=lab)
          for k2, lab in (("published", "as published (probable errors × 1.483)"),
                          ("harvey", "plates remeasured in 1979 (Harvey)"),
                          ("ours", "this project, from the printed measures"))]
    ax.legend(handles=hs, loc="lower left", fontsize=7.4, bbox_to_anchor=(0.0, 1.0), ncol=3, frameon=False)
    ax.set_title("Every determination of the bending, then and now", pad=26)
    P.caption(fig, "Published errors before about 1950 are probable errors; here all are one sigma. The 1920 paper gave "
                   "no error for the astrographic plates. Our least-squares fits use the printed measures: 1920 (with "
                   "the two Table II misprints replaced by the values the 1920 reduction used), and 1923 and 1928 for "
                   "Wallal. Every plate's constants are fitted together with the deflection, except for the 15-foot "
                   "cameras, whose 1923 table gives only each star's mean over the four plate pairs. Where a row says "
                   "'jackknife', refitting with one plate left out at a time gave a larger error than the formal one; "
                   "the 1919 combinations use that error.")
    return P.save(fig, _out("05_every_determination"))


# ---------------------------------------------------------------------------
# 06 The astrographic plates: deflection and scale cannot be told apart
# ---------------------------------------------------------------------------

def fig_06_astrographic(res: dict) -> list:
    import eclipse1919 as E

    a = res["astrographic"]
    t9 = E.load_astro()["t9"]
    e = np.array(a["alpha_vs_e"]["e_rev"])
    lim = np.array(a["alpha_vs_e"]["limb"])
    gr = C.DEFLECTION_LIMB_1919
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    # the 1/r law's two predictions
    ax.axhline(gr, **P.GR_LINE)
    ax.axhline(gr / 2, **P.NEWTON_LINE)
    ax.text(0.05, gr + 0.04, "Einstein: 1.75″", ha="left", va="bottom", fontsize=8, color=P.INK)
    ax.text(-0.0135, gr / 2 - 0.05, "Newtonian half: 0.87″", ha="left", va="top", fontsize=8, color=P.INK_SECONDARY)
    # each plate's own solution (Table IX): deflection against the scale it found
    ax.scatter(t9["e_rev"].astype(float), t9["alpha_limb_arcsec"].astype(float), s=26, marker="D",
               color=P.AQUA, edgecolor=P.SURFACE, linewidth=0.6, zorder=4,
               label="the 18 solutions of Table IX (16 plates, each solved alone, 1920)")
    # what all 16 plates say, for every assumed scale (the p.311 combination)
    ax.plot(e, lim, color=P.AQUA, lw=2.0, label="all 16 plates, for an assumed scale change e")
    # the two choices of 1920
    for ee, lab, xy in ((0.082, "1920 choice: the plates' own\nmean scale e = 0.082 r → 0.93″", (0.092, 0.15)),
                        (0.020, "if only refraction and aberration\nchanged the scale (e = 0.020 r) → 1.52″", (0.026, 2.17))):
        y = float(np.interp(ee, e, lim))
        ax.plot([ee], [y], "o", color=P.INK, ms=6, mfc=P.SURFACE, mew=1.4, zorder=5)
        P.annotate(ax, lab, (ee, y), xy, va="center")
    # our joint fit with e free (one e for all plates), with its 1-sigma ellipse
    fc = a["fits"]["scale_common_free"]
    from matplotlib.patches import Ellipse

    cov = np.array([[fc["e_err_rev"] ** 2, fc["corr_alpha_e"] * fc["e_err_rev"] * fc["limb_err"]],
                    [fc["corr_alpha_e"] * fc["e_err_rev"] * fc["limb_err"], fc["limb_err"] ** 2]])
    w, v = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse((fc["e_rev"], fc["limb"]), 2 * np.sqrt(w[1]), 2 * np.sqrt(w[0]), angle=ang,
                         fill=False, ec=P.INK, lw=1.0, ls=(0, (3, 2)), zorder=4))
    # the scale needed for each theory
    for target, lab in ((gr, "Einstein's 1.75″\nneeds e = {} r"), (gr / 2, "the Newtonian\nhalf needs\ne = {} r")):
        ee = float(np.interp(target, lim[::-1], e[::-1]))
        ax.plot([ee, ee], [-0.3, target], color=P.MUTED, lw=0.7, ls=(0, (1, 2)))
        ax.text(ee + 0.0012, -0.22, lab.format(f"{ee:+.3f}".replace("-", "−")), ha="left", va="bottom", fontsize=7,
                color=P.INK_SECONDARY)
    ax.set_xlim(-0.015, 0.175)             # Table IX reaches e = 0.166 r (plate 7 against 18_3)
    ax.set_ylim(-0.25, 2.6)
    ax.set_xlabel("assumed change of scale between eclipse and comparison plates, e "
                  "(micrometer revolutions per 50′; 0.01 r = 0.12″ = 4 × 10⁻⁵)")
    ax.set_ylabel("deflection at the Sun's limb (″)")
    ax.set_title("Sobral astrographic: the answer depends on the plate scale you assume")
    ax.legend(loc="upper right", fontsize=7.6)
    corr = f"{fc['corr_alpha_e']:+.2f}".replace("-", "−")
    P.caption(fig, f"Declinations of the five bright stars on 16 eclipse plates (p.311) give α + 0.24 e for each plate; "
                   f"the line is their mean. Dashed ellipse: our least-squares fit with one common e (the Δχ² = 1 "
                   f"contour, whose extent along each axis is that parameter's one-sigma range): "
                   f"{fc['limb']:.2f} ± {fc['limb_err']:.2f}″ at e = {fc['e_rev']:.3f} ± {fc['e_err_rev']:.3f} r "
                   f"(correlation {corr}); refitting with one plate left out at a time gives "
                   f"± {fc['limb_err_jackknife']:.2f}″. Harvey (1979), remeasuring all 16 plates and reducing the 9 that show ten or more stars, in both coordinates "
                   "with the scale free on each, found 1.55 ± 0.34″.")
    return P.save(fig, _out("06_astrographic_scale"))


# ---------------------------------------------------------------------------
# 09 A century of gamma
# ---------------------------------------------------------------------------

FAMILY_STYLE = {"starlight": ("o", P.BLUE, "starlight bent at eclipses (and Hipparcos, 1997)"),
                "radio": ("s", P.ORANGE, "radio waves bent: interferometers and VLBI"),
                "delay": ("D", P.AQUA, "radar and spacecraft signals delayed")}


def fig_09_century(res: dict) -> list:
    """The one-sigma uncertainty of gamma from every independent experiment, 1919-2018,
    with the best-so-far ladder and this project's own re-analyses."""
    t = pd.DataFrame(res["history"]["timeline"])
    fr = pd.DataFrame(res["history"]["frontier"])
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    # best-so-far staircase (first-hand, unflagged results)
    xs, ys = [], []
    for i, r in fr.iterrows():
        nxt = fr["year"].iloc[i + 1] if i + 1 < len(fr) else 2026
        xs += [r["year"], nxt]
        ys += [r["sigma_gamma"], r["sigma_gamma"]]
    ax.plot(xs, ys, color=P.MUTED, lw=1.0, drawstyle="steps-post", zorder=1)
    for fam, (mk, col, lab) in FAMILY_STYLE.items():
        d = t[t["family"] == fam]
        hollow = d["flag"].astype(bool)
        ax.scatter(d.loc[~hollow, "year"], d.loc[~hollow, "sigma_gamma"], marker=mk, s=34, color=col,
                   edgecolor=P.SURFACE, linewidth=0.6, zorder=4, label=lab)
        ax.scatter(d.loc[hollow, "year"], d.loc[hollow, "sigma_gamma"], marker=mk, s=34, facecolor=P.SURFACE,
                   edgecolor=col, linewidth=1.1, zorder=4)
    # this project's re-analyses of the original data. The 1919 and 1922 stars are drawn 2.5 years
    # either side of their year: closer, their points touch the published points of those years
    shift = {1919: -2.5, 1922: 2.5}
    ours = [(yr + shift.get(yr, 0.0), sg, lab, off) for yr, sg, lab, off in _our_gamma_points(res)]
    if ours:
        ax.scatter([o[0] for o in ours], [o[1] for o in ours], marker="*", s=150, color=P.INK,
                   edgecolor=P.SURFACE, linewidth=0.8, zorder=6, label="this project, from the original measurements")
        for yr, sg, lab, off in ours:
            ha = "left" if off[0] >= 0 else "right"
            ax.annotate(lab, (yr, sg), xytext=off, textcoords="offset points", fontsize=7.4, color=P.INK, ha=ha,
                        va="center", path_effects=HALO,
                        arrowprops={"arrowstyle": "-", "color": P.INK_SECONDARY, "lw": 0.6, "shrinkA": 1, "shrinkB": 6,
                                    # leave from the text's corner nearest the star
                                    "relpos": (1.0 if ha == "right" else 0.0, 0.0 if off[1] > 0 else 1.0)}
                        if np.hypot(*off) > 15 else None)
    notes = [("E1919-sobral-4inch", "1919 eclipse,\nSobral 4-inch", (-14, 22), "right"),
             ("T1976-viking", "Viking landers\non Mars, 1976", (-10, 0), "right"),
             ("R1974-nrao-fomalont", "radio interferometers\nwatch quasars pass the Sun", (-8, -14), "right"),
             ("A1997-hipparcos", "Hipparcos", (8, 0), "left"),
             ("V2011-lambert", "VLBI, 2011", (8, 3), "left"),
             ("T2002-cassini", "Cassini, 2002:\n2.3 × 10⁻⁵", (-10, -6), "right")]
    for rid, txt, off, ha in notes:
        r = t[t["id"] == rid].iloc[0]
        ax.annotate(txt, (r["year"], r["sigma_gamma"]), xytext=off, textcoords="offset points", fontsize=7.4,
                    ha=ha, va="center", color=P.INK_SECONDARY,
                    arrowprops={"arrowstyle": "-", "color": P.MUTED, "lw": 0.6, "shrinkA": 1, "shrinkB": 4,
                                # leave from the text's corner nearest the point, not through the text
                                "relpos": (1.0 if ha == "right" else 0.0, 0.0 if off[1] > 0 else 1.0)}
                    if np.hypot(*off) > 15 else None)
    ax.set_yscale("log")
    ax.set_ylim(8e-6, 1.5)
    ax.set_xlim(1900, 2026)
    ax.set_xlabel("year of the observations")
    ax.set_ylabel("one-sigma uncertainty of γ\n(how tightly it is pinned to Einstein's γ = 1)")
    ax.set_title("A century of γ: from ±20% in 1919 to ±0.002% with Cassini")
    ax.legend(loc="lower left", fontsize=7.6)
    P.caption(fig, "One point per experiment (and one re-reduction), from the table verified in project N7 "
                   "(64-gravity-scorecard/literature/checks/gamma.csv; copy in transcribed/). Probable errors (the custom "
                   "before ~1950) converted to one sigma (× 1.483). Open symbols: re-reductions, numbers known only from "
                   "secondary sources, or results their authors flag as suspect. Grey staircase: the best first-hand "
                   "result so far. Our 1919 and 1922 refits are drawn 2.5 years before and after their year, clear of the "
                   "published points.")
    return P.save(fig, _out("08_century_of_gamma"))


def _our_gamma_points(res: dict) -> list:
    """(year, sigma_gamma, label, offset) for this project's own re-analyses."""
    out = []
    s4 = res.get("sobral_4inch", {})
    if s4.get("joint_fit", {}).get("both"):
        # the adopted error, as in analysis.step_combined_1919 and figure 05: the larger of the formal
        # error and the jackknives (over plates, the larger; see the README)
        jk = s4.get("jackknife", {})
        err = max(s4["joint_fit"]["both"]["limb_err"], jk.get("err_plates", 0.0), jk.get("err_stars", 0.0))
        out.append((1919, 2 * err / C.DEFLECTION_LIMB_GR, "our refit of the\nSobral 4-inch measures", (10, -26)))
    wal = res.get("wallal", {}).get("combined")
    if wal:
        out.append((1922, 2 * wal["beta_err"] / C.DEFLECTION_LIMB_GR, "our refit of\nWallal 1922", (12, 14)))
    cas = (res.get("cassini") or {}).get("summary") or {}
    if cas.get("gamma_err_total"):
        out.append((2002, cas["gamma_err_total"], "our fit of Cassini's\nraw Doppler data", (-12, 14)))
    return out


# ---------------------------------------------------------------------------
# GR's deflection across the sky; the Shapiro delay to Cassini
# ---------------------------------------------------------------------------

def fig_deflection_vs_elongation(regimes: dict) -> list:
    """GR's deflection from the Sun's limb to the anti-Sun, and where each kind of experiment
    looks. ``regimes``: label -> (min_deg, max_deg, label height or None, strip or None): a
    full-height band labelled at that data height, or a strip along the bottom starting at
    axes fraction ``strip``, labelled inside it (so overlapping ranges stay readable)."""
    import stars as S

    chi = np.logspace(np.log10(0.2667), np.log10(179.9), 600)
    gr = S.deflection_arcsec(np.deg2rad(chi), 1.0)
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.loglog(chi, gr, label="general relativity", **P.GR_LINE)
    ax.loglog(chi, gr / 2, label="Newtonian light corpuscles (half)", **P.NEWTON_LINE)
    ymin, ymax = 3e-4, 4.0
    for label, (lo, hi, ylab, strip) in regimes.items():
        if strip is None:          # a full-height band, label at data height ylab
            ax.axvspan(lo, hi, color=P.WASH, lw=0, zorder=0)
            ax.text(np.sqrt(lo * hi), ylab, label, ha="center", va="bottom", fontsize=7.5, color=P.INK_SECONDARY)
        else:                      # a strip along the bottom, from axes fraction `strip`, label inside it
            ax.axvspan(lo, hi, ymin=strip, ymax=strip + 0.05, color=P.WASH, lw=0, zorder=0.6)   # over the grid
            ax.text(np.sqrt(lo * hi), strip + 0.025, label, ha="center", va="center", fontsize=7.2,
                    color=P.INK_SECONDARY, transform=ax.get_xaxis_transform())
    limb = S.deflection_arcsec(np.deg2rad(0.2667), 1.0)
    ax.plot([0.2667], [limb], "o", color=P.INK, ms=5)
    P.annotate(ax, f"{limb:.2f}″ for a ray grazing the Sun", (0.2667, limb), (0.5, 2.6), va="center")
    d90 = S.deflection_arcsec(np.pi / 2, 1.0)
    ax.plot([90], [d90], "o", color=P.INK, ms=5)
    P.annotate(ax, f"{d90 * 1e3:.2f} milliarcseconds at 90° from the Sun:\n"
                   "a human hair seen from 4 km. Gaia corrects\nevery star it measures for this.",
               (90, d90), (1.2, 1.6e-3), va="center")
    ax.set_xlim(0.2, 200)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("angular distance of the star from the Sun's centre (degrees)")
    ax.set_ylabel("apparent outward shift of the star (arcsec)")
    ax.set_title("The Sun bends the light of every star in the sky")
    ax.legend(loc="upper right")
    P.caption(fig, "δ = (1+γ) GM/(c² d) cot(χ/2) for an observer at d = 1 au (γ = 1); "
                   "near the Sun this is the familiar 4GM/(c²b) = 1.75″ × (R☉/b).")
    return P.save(fig, _out("07_deflection_vs_elongation"))


def _sci(v: float) -> str:
    """2.2e-11 -> '2 × 10⁻¹¹' (one significant figure, Unicode superscripts)."""
    e = int(np.floor(np.log10(abs(v))))
    m = round(v / 10 ** e)
    if m == 10:
        m, e = 1, e + 1
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m} × 10{sup}"


def _cassini_series() -> pd.DataFrame | None:
    """The measured two-way frequency shift per time bin, from cassini_doppler.py (if it has
    run): the main data set only (the plasma-free X/X + Ka/Ka combination)."""
    path = C.CACHE_DIR / "cassini_measured_series.parquet"
    if not path.exists():
        return None
    d = pd.read_parquet(path)
    if "dataset" in d:
        main = d["dataset"].astype(str).str.startswith("main")
        d = d[main] if main.any() else d
    d = d.rename(columns={"date_utc": "date"})
    return d.reset_index(drop=True)


def fig_shapiro(sh: pd.DataFrame, summary: dict, cassini: dict | None = None) -> list:
    """Cassini 2002: the round-trip delay, and the frequency shift it causes."""
    ms_all = _cassini_series()
    if ms_all is not None:
        fig, axes = plt.subplots(4, 1, figsize=(7.4, 10.6), sharex=True,
                                 gridspec_kw={"height_ratios": [0.9, 1.15, 0.9, 0.75], "hspace": 0.13})
    else:
        fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.4), sharex=True,
                                 gridspec_kw={"height_ratios": [1.0, 1.15, 0.8], "hspace": 0.12})
    import shapiro

    for ax in axes:
        # SCE1 data: 6 June to the end of 5 July 2002 (PDS DATASET.CAT START/STOP_TIME)
        ax.axvspan(pd.Timestamp(shapiro.SCE1_START), pd.Timestamp(shapiro.SCE1_STOP), color=P.WASH, lw=0, zorder=0)
    ax = axes[0]
    ax.plot(sh["date"], sh["delay_2way_us"], color=P.BLUE, lw=1.8)
    i = int(np.argmin(sh["b_rsun"]))
    ax.plot(sh["date"].iloc[i], sh["delay_2way_us"].iloc[i], "o", color=P.BLUE, ms=6, mec=P.SURFACE, mew=1.2)
    P.annotate(ax, f"{sh['delay_2way_us'].iloc[i]:.0f} µs on 21 June:\nthe radio waves passed\n"
                   f"{summary['b_min_rsun']:.2f} solar radii from the Sun's centre",
               (sh["date"].iloc[i], sh["delay_2way_us"].iloc[i]),
               (sh["date"].iloc[i] + pd.Timedelta(6.5, unit="D"), sh["delay_2way_us"].iloc[i] - 18), va="top")
    ax.set_ylabel("extra round-trip time (µs)")
    ax.set_title("Cassini 2002: radio signals slowed by the Sun's gravity")
    ax.set_ylim(0, sh["delay_2way_us"].max() * 1.12)
    ax.text(pd.Timestamp("2002-06-06") + pd.Timedelta(12, unit="h"), 8, "Cassini tracked\n6 June – 5 July",
            fontsize=7.5, color=P.INK_SECONDARY, va="bottom")

    ax = axes[1]
    ax.axhline(0, color=P.AXIS, lw=0.8)
    dnu = (-sh["y_2way"] * 1e10).to_numpy()     # a growing delay lowers the received frequency
    ax.plot(sh["date"], dnu, color=P.ORANGE, lw=1.8, label="general relativity (γ = 1)", zorder=3)
    ms = ms_all
    if ms is not None:
        # the measured shift: Cassini's two-way Doppler with the fitted trajectory corrections,
        # plasma and other nuisance terms removed (cassini_doppler.py). Its sign convention is
        # matched to this panel's by comparing its own GR curve with ours.
        t_ms = pd.to_datetime(ms["date"])
        gr_here = np.interp(t_ms.astype("int64"), sh["date"].astype("int64"), sh["y_2way"].to_numpy())
        sgn = 1.0 if np.dot(gr_here, ms["y_gr"].to_numpy()) >= 0 else -1.0
        ax.errorbar(t_ms, -sgn * ms["y_measured"] * 1e10, yerr=ms["sigma"] * 1e10, fmt="o", ms=2.6,
                    color=P.BLUE, elinewidth=0.6, alpha=0.85, zorder=4,
                    label="measured by Cassini (our fit)")
        ax.legend(loc="lower right", fontsize=7.2, frameon=True, facecolor=P.SURFACE, edgecolor="none",
                  framealpha=0.85)            # over the grid and the shading's edge, no data
    ax.set_ylabel("two-way frequency shift  Δν/ν  (10⁻¹⁰)")
    j, k = int(np.argmin(dnu)), int(np.argmax(dnu))
    P.annotate(ax, "path closing in on the Sun:\nthe delay grows, so the\nsignal arrives slightly redder",
               (sh["date"].iloc[j], dnu[j]), (sh["date"].iloc[j] - pd.Timedelta(19, unit="D"), dnu[j] * 0.55),
               va="center")
    P.annotate(ax, "path moving away:\nthe delay shrinks, so bluer", (sh["date"].iloc[k], dnu[k]),
               (sh["date"].iloc[k] + pd.Timedelta(4, unit="D"), dnu[k] * 0.7), va="center")
    ax.text(0.01, 0.97, "Bertotti, Iess & Tortora (2003) fitted\nthis curve to Cassini's X- and Ka-band\n"
            "Doppler data: its amplitude, ∝ (1 + γ)/2,\nmatched relativity to 1.2 parts in\n"
            "100,000: γ − 1 = (2.1 ± 2.3) × 10⁻⁵.",
            transform=ax.transAxes, fontsize=7, color=P.INK_SECONDARY, va="top")

    if ms is not None:
        # measured minus GR (gamma = 1) in units of each point's error: the errors span three
        # decades (1e-14 far from the Sun, 2e-11 inside 3 solar radii), so plain residuals
        # would show only the few points nearest the Sun
        ax = axes[2]
        ax.axhline(0, color=P.AXIS, lw=0.8)
        z = -sgn * (ms["y_measured"] - ms["y_gr"]) / ms["sigma"]
        ax.plot(pd.to_datetime(ms["date"]), z, "o", ms=2.2, color=P.BLUE, alpha=0.7, mew=0)
        ax.set_ylabel("measured minus GR,\nin units of each\npoint's error")
        ax.set_ylim(-5.2, 5.2)
        sig = ms["sigma"].to_numpy()
        far, near = np.median(sig[ms["b_rsun"] > 20]), np.median(sig[ms["b_rsun"] < 3])
        cas = (cassini or {}).get("summary") or {}
        if cas.get("gamma_minus_1") is not None:
            g = f"{cas['gamma_minus_1'] * 1e4:+.1f}".replace("-", "−")
            ax.text(0.01, 0.95, f"our fit: γ − 1 = ({g} ± {cas['gamma_err_total'] * 1e4:.1f}) × 10⁻⁴;"
                    "  Bertotti et al. 2003: (+0.21 ± 0.23) × 10⁻⁴", transform=ax.transAxes, fontsize=7.2,
                    color=P.INK_SECONDARY, va="top", path_effects=HALO)
        ax.text(0.01, 0.05, f"{len(z)} points of 600 s; rms {np.sqrt(np.mean(z ** 2)):.2f}.  Typical error: "
                f"{_sci(far)} far from the Sun, {_sci(near)} within 3 solar radii",
                transform=ax.transAxes, fontsize=7.2, color=P.INK_SECONDARY, va="bottom", path_effects=HALO)
    ax = axes[-1]
    ax.plot(sh["date"], sh["b_rsun"], color=P.INK_SECONDARY, lw=1.4)
    ax.set_yscale("log")
    ax.set_ylabel("closest approach of the\nradio path (solar radii)")
    ax.axhline(1.0, color=P.MUTED, lw=0.8, ls=(0, (2, 2)))
    ax.text(pd.Timestamp("2002-07-09"), 1.07, "the Sun's surface", fontsize=7, color=P.MUTED, va="bottom")
    ax.set_ylim(0.8, 200)
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO, interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.set_xlim(pd.Timestamp("2002-05-24"), pd.Timestamp("2002-07-20"))
    P.caption(fig, "Geometry: JPL Horizons heliocentric vectors of Cassini (-82) and the Earth, with light time; "
                   "delay = (1+γ)GM/c³ ln[(r₁+r₂+R)/(r₁+r₂−R)] each way (γ = 1); "
                   "shaded: the Solar Conjunction Experiment (PDS volumes CORS_0021–0028).")
    return P.save(fig, _out("09_shapiro_cassini"))


# ---------------------------------------------------------------------------

#: Where each kind of experiment measures the bending: angular distance from the Sun (deg), then
#: (label height, None) for a full-height band or (None, axes fraction) for a strip along the bottom.
#: The 13 Table I stars of 1919 lie 0.35-1.56 deg from the Sun at Sobral (Principe measured six of
#: them, out to 1.51 deg); Wallal's five-foot stars of 1922 lie 0.55-10.3 deg from it.
REGIMES = {"eclipse photographs (1919: to 1.6°; Wallal 1922: to 10°)": (0.27, 10.3, None, 0.065),
           "radio waves from quasars": (2.5, 15.0, None, 0.01),
           "Hipparcos (1989-93)\nand Gaia (2014-25)\nscan here: 9.8 mas\nat 45°, 1.7 at 135°": (45.0, 135.0, 0.045, None)}


def fig_08_elongation(res: dict) -> list:
    return fig_deflection_vs_elongation(REGIMES)


def fig_10_shapiro(res: dict) -> list:
    import shapiro

    sh = shapiro.conjunction()
    return fig_shapiro(sh, shapiro.summary(sh), res.get("cassini"))


FIGURES = {"01": fig_01_field, "02": fig_02_plate1, "03": fig_03_star_shifts, "04": fig_04_radial_law,
           "05": fig_05_forest, "06": fig_06_astrographic, "07": fig_08_elongation, "08": fig_09_century,
           "09": fig_10_shapiro}


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    res = _results()
    for tag, fn in FIGURES.items():
        if argv and not any(tag.startswith(a) for a in argv):
            continue
        out = fn(res)
        print("  " + ", ".join(str(o).rsplit("/", 1)[-1] for o in out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
