"""Every static figure, drawn from the cleaned catalogue and the results dictionary.

Figures are numbered in reading order:

    01_world_map          where meteorites are seen to fall, and where they are found
    02_collection_areas   the dense collection areas, and Antarctica's stranding zones
    03_witnesses          witnessed falls follow people, country by country
    04_what_we_miss       recorded falls against the true rate at which they land
    05_what_is_found      what falls vs what each kind of search finds
    06_masses             how heavy the meteorites in each collection are
    07_history            how the catalogue grew

Every map uses an equal-area projection, so equal patches of map are equal areas
of ground and a density of dots can be read at face value.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

import analysis
import catalog
import geo
import plotstyle as ps
from config import FIGURE_DIR, MANIFEST_PATH, MILESTONES

ps.use_style()

FINDS_CMAP = LinearSegmentedColormap.from_list("finds", ps.BLUE_RAMP[3:])
LABEL_BOX = dict(boxstyle="round,pad=0.18", fc=ps.SURFACE, ec="none", alpha=0.82)


def _provenance(short: bool = False) -> str:
    try:
        m = json.loads(MANIFEST_PATH.read_text())["meteorites"]
        got = m["retrieved_utc"][:10]
    except (OSError, KeyError, ValueError):
        got = "?"
    if short:
        return (f"Data: Meteoritical Society catalogue (2013 snapshot via NASA Open Data, "
                f"retrieved {got}).")
    return (f"Data: The Meteoritical Society's catalogue of meteorites, 2013 snapshot distributed by "
            f"the NASA Open Data Portal (retrieved {got}); country outlines: Natural Earth.")


def _fmt(n) -> str:
    return f"{int(round(n)):,}"


def _mass(g: float) -> str:
    if g >= 1e6:
        return f"{g / 1e6:g} t"
    if g >= 1e3:
        return f"{g / 1e3:.3g} kg"
    if g >= 1:
        return f"{g:.3g} g"
    return f"{g:g} g"


def _title(ax, title: str, subtitle: str | None = None, *, lines: int = 1,
           fontsize: float = 8.2) -> None:
    """A left-aligned title with an optional subtitle tucked between it and the axes."""
    pad = 7 + (lines * fontsize * 1.3 + 3 if subtitle else 0)
    ax.set_title(title, loc="left", pad=pad)
    if subtitle:
        ax.annotate(subtitle, xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom", fontsize=fontsize,
                    color=ps.INK_SECONDARY, linespacing=1.3)


def _caption(fig, text: str) -> None:
    """Footnote placed just under the lowest element of the figure."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    y0 = min(a.get_tightbbox(r).y0 for a in fig.axes if a.get_visible())
    y = fig.transFigure.inverted().transform((0, y0))[1]
    fig.text(0.5, y - 0.012, text, ha="center", va="top", fontsize=7.5, color=ps.MUTED,
             wrap=True)


# --- Map drawing ------------------------------------------------------------------------


def _ring_path(xy: np.ndarray, *, ccw: bool) -> tuple[np.ndarray, list]:
    """A closed ring as path vertices + codes, wound the requested way (holes are
    wound opposite to exteriors so they stay unfilled)."""
    signed = geo.shoelace(xy[:, 0], xy[:, 1])
    if (signed > 0) != ccw:
        xy = xy[::-1]
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(xy) - 2) + [MplPath.CLOSEPOLY]
    return xy, codes


def _land_path(countries, proj) -> MplPath:
    verts, codes = [], []
    for c in countries:
        for ext, holes in c.polygons:
            for ring, ccw in [(ext, True)] + [(h, False) for h in holes]:
                x, y = proj(ring[:, 0], ring[:, 1])
                v, k = _ring_path(np.column_stack([x, y]), ccw=ccw)
                verts.append(v)
                codes.extend(k)
    return MplPath(np.vstack(verts), codes)


def draw_world(ax, countries, *, borders: bool = True):
    """Equal Earth base map: ocean, graticule, land."""
    ox, oy = geo.equal_earth(*geo.world_outline())
    ax.fill(ox, oy, color=ps.OCEAN, zorder=0, lw=0)
    for lon, lat in geo.graticule(30):
        gx, gy = geo.equal_earth(lon, lat)
        ax.plot(gx, gy, color=ps.GRID, lw=0.4, zorder=1)
    ax.add_patch(PathPatch(_land_path(countries, geo.equal_earth), facecolor=ps.LAND,
                           edgecolor=ps.COAST if borders else "none", lw=0.25, zorder=2))
    ax.plot(ox, oy, color=ps.AXIS, lw=0.8, zorder=6)
    ax.set_aspect("equal")
    ax.set_xlim(ox.min() * 1.005, ox.max() * 1.005)
    ax.set_ylim(oy.min() * 1.01, oy.max() * 1.01)
    ax.axis("off")


def _map_label(ax, text, lon, lat, dx, dy, *, ha="left", fontsize=7.5, proj=geo.equal_earth):
    x, y = proj(lon, lat)
    relx = {"left": 0.0, "right": 1.0, "center": 0.5}[ha]
    return ps.annotate(ax, text, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                       ha=ha, va="center", fontsize=fontsize, color=ps.INK_SECONDARY,
                       relpos=(relx, 0.5), zorder=8, bbox=LABEL_BOX)


# --- 01 World map ---------------------------------------------------------------------------


def figure_world_map(df, countries, results) -> list:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.5, 11.0))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01, hspace=0.13)

    # (a) Falls.
    draw_world(ax1, countries)
    f = df[df["valid"] & df["located"] & df["environment"].eq("fall")]
    x, y = geo.equal_earth(f["lon"], f["lat"])
    ax1.scatter(x, y, s=11, c=ps.C_FALL, edgecolors=ps.SURFACE, linewidths=0.45, zorder=5)
    _title(ax1, f"a   {_fmt(len(f))} witnessed falls: seen to fall, then picked up",
           "They trace where people live and look up: Europe, India, Japan, eastern China "
           "and the United States.", fontsize=8.8)

    # (b) Finds, as a density on equal-area hexagons.
    draw_world(ax2, countries)
    g = df[df["valid"] & df["located"] & df["environment"].isin(["antarctic", "desert", "chance"])]
    x, y = geo.equal_earth(g["lon"], g["lat"])
    ox, oy = geo.equal_earth(*geo.world_outline())
    hb = ax2.hexbin(x, y, gridsize=(150, 42), extent=(ox.min(), ox.max(), oy.min(), oy.max()),
                    norm=LogNorm(vmin=1, vmax=5000), mincnt=1, cmap=FINDS_CMAP,
                    linewidths=0.0, zorder=4)
    hex_area = abs(geo.shoelace(*hb.get_paths()[0].vertices[:-1].T)) * geo.EARTH_RADIUS_KM**2
    ca = results["collection_areas"]
    _title(ax2, f"b   {_fmt(len(g))} finds with a known find site",
           f"They trace where rocks survive and stand out: "
           f"{ca['antarctic_share_of_located_finds']['fraction']:.0%} come from Antarctica, "
           f"which is {ca['antarctica_share_of_land']:.0%} of the world's land.", fontsize=8.8)
    cax = ax2.inset_axes([0.07, 0.45, 0.17, 0.022])
    cb = fig.colorbar(hb, cax=cax, orientation="horizontal")
    cb.set_label(f"finds per hexagon\n(each {hex_area / 1e3:,.0f},000 km²)", fontsize=7,
                 color=ps.INK_SECONDARY, labelpad=3)
    cb.ax.tick_params(labelsize=6.5, length=2, color=ps.MUTED)
    cb.outline.set_linewidth(0)
    cax.xaxis.set_label_position("top")

    top = {t["area"]: t for t in ca["top"]}
    dens = {d["country"]: d for d in ca["all_country_density"]}
    nwa = top["Northwest Africa"]
    labels = [
        (f"Oman: {_fmt(dens['Oman']['finds'])} finds\n(Dhofar, Jiddat al Harasis, ...)",
         55, 19, 30, 22, "left"),
        (f"Libya: {_fmt(dens['Libya']['finds'])}\n(Dar al Gani, Hammadah al Hamra)",
         15, 27, 4, -44, "center"),
        (f"Northwest Africa: {_fmt(nwa['n'])} meteorites bought\nthrough dealers; only "
         f"{nwa['located']} have a find site", -6, 30.5, -40, -36, "right"),
        (f"Atacama, Chile: {_fmt(dens['Chile']['finds'])}", -70, -24.5, -30, -8, "right"),
        ("Nullarbor Plain, Australia", 128, -30.5, 34, -10, "left"),
        ("dry lakes of the\nUS Southwest", -116, 35, -46, -26, "right"),
        (f"Yamato & Asuka (Japanese teams): {_fmt(top['Yamato']['n'] + top['Asuka']['n'])}",
         34, -71.8, -14, 26, "right"),
        ("Transantarctic Mountains\n(US teams, since 1976)", 160, -80, -10, 36, "right"),
    ]
    for text, lon, lat, dx, dy, ha in labels:
        _map_label(ax2, text, lon, lat, dx, dy, ha=ha)

    dq = results["data_quality"]
    _caption(fig, f"Equal Earth projection: equal areas of map are equal areas of ground. "
             f"Not shown: {_fmt(dq['no_location'])} meteorites without usable coordinates "
             f"({_fmt(dq['coords_missing'])} blank, {_fmt(dq['coords_placeholder_0_0'])} at the "
             f"placeholder 0°, 0°, {dq['coords_zeroed_latitude']} Yamato meteorites with their "
             f"latitude zeroed). {_provenance()}")
    return ps.save(fig, FIGURE_DIR / "01_world_map")


# --- 02 Collection areas and Antarctica ------------------------------------------------------


def figure_collection_areas(df, countries, results) -> list:
    ca = results["collection_areas"]
    fig = plt.figure(figsize=(12.5, 6.4))
    ax1 = fig.add_axes([0.06, 0.08, 0.33, 0.8])
    ax2 = fig.add_axes([0.45, 0.0, 0.54, 0.88])
    ax2.set_anchor("N")

    # (a) The top collection areas.
    top = ca["top"][:16][::-1]
    ypos = np.arange(len(top))
    colors = [ps.ENV_COLORS[t["environment"]] for t in top]
    ax1.barh(ypos, [t["n"] for t in top], height=0.62, color=colors, zorder=3)
    for yy, t in zip(ypos, top):
        ax1.text(t["n"] + 90, yy, _fmt(t["n"]), va="center", fontsize=7.5,
                 color=ps.INK_SECONDARY)
    ax1.set_yticks(ypos, [t["area"] for t in top], fontsize=8)
    ax1.tick_params(axis="y", length=0, colors=ps.INK_SECONDARY)
    ax1.set_xlim(0, max(t["n"] for t in top) * 1.14)
    ax1.grid(axis="y", visible=False)
    ax1.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax1.set_xlabel("meteorites named after the area")
    _title(ax1, "a   The 16 richest collection areas",
           f"Each is one ice field or desert plain, searched again and again.\n"
           f"The top 10 hold {ca['top10_share_of_finds']:.0%} of all {_fmt(ca['n_finds'])} finds.",
           lines=2)
    handles = [plt.Rectangle((0, 0), 1, 1, color=ps.C_ANTARCTIC),
               plt.Rectangle((0, 0), 1, 1, color=ps.C_DESERT)]
    ax1.legend(handles, ["Antarctic ice", "hot desert"], loc="lower right", fontsize=8)

    # (b) Antarctica, south polar equal-area.
    proj = geo.south_polar
    ant = next(c for c in countries if c.name == "Antarctica")
    ax2.add_patch(PathPatch(_land_path([ant], proj), facecolor=ps.LAND, edgecolor=ps.COAST,
                            lw=0.4, zorder=2))
    for lat0 in (-60, -70, -80):
        lo = np.linspace(-180, 180, 361)
        ax2.plot(*proj(lo, np.full_like(lo, lat0)), color=ps.GRID, lw=0.5, zorder=1)
    for lon0 in range(-180, 180, 30):
        la = np.linspace(-90, -60, 50)
        ax2.plot(*proj(np.full_like(la, lon0), la), color=ps.GRID, lw=0.5, zorder=1)

    a = df[df["valid"] & df["environment"].eq("antarctic")]
    areas = (a[a["area"].notna()].groupby("area")
             .agg(n=("id", "size"),
                  lat=("lat", lambda s: s[a.loc[s.index, "located"]].median()),
                  lon=("lon", lambda s: s[a.loc[s.index, "located"]].median()))
             .dropna().sort_values("n", ascending=False))
    px, py = proj(areas["lon"], areas["lat"])
    ax2.scatter(px, py, s=areas["n"] * 0.16 + 4, c=ps.C_ANTARCTIC, alpha=0.55,
                edgecolors=ps.SURFACE, linewidths=0.8, zorder=4)
    areas["x"], areas["y"] = px, py
    unplaced = int(a["area"].notna().sum() - areas["n"].sum())

    # Areas elsewhere on the continent get individual labels...
    single = {"Yamato": (-12, 30, "right"), "Asuka": (-26, 10, "right"),
              "Grove Mountains": (18, 16, "left"), "LaPaz Icefield": (-34, 24, "right"),
              "Pecora Escarpment": (-40, 6, "right"), "Thiel Mountains": (-36, -14, "right")}
    for name, (dx, dy, ha) in single.items():
        r = areas.loc[name]
        _map_label(ax2, f"{name} {_fmt(r['n'])}", r["lon"], r["lat"], dx, dy, ha=ha,
                   fontsize=7.2, proj=proj)
    # ...while the crowded Transantarctic Mountains get a column of labels.
    ross = areas[(areas["lon"] > 140) | (areas["lon"] < -140)]
    ross = ross[ross["n"] >= 300].sort_values("y", ascending=False)
    col_x = 0.285
    ys = np.linspace(-0.075, -0.45, len(ross))
    for (name, r), yl in zip(ross.iterrows(), ys):
        ax2.annotate(f"{name} {_fmt(r['n'])}", xy=(r["x"], r["y"]), xytext=(col_x, yl),
                     textcoords="data", ha="left", va="center", fontsize=7.2,
                     color=ps.INK_SECONDARY, zorder=7, bbox=LABEL_BOX,
                     arrowprops=dict(arrowstyle="-", color=ps.INK_SECONDARY, lw=0.6,
                                     shrinkA=1, shrinkB=3, relpos=(0, 0.5)))
    ax2.text(col_x, -0.042, "Transantarctic Mountains:", fontsize=7.4, color=ps.INK,
             ha="left", va="center", zorder=7)

    x0, y0 = proj(0, -90)
    ax2.plot(x0, y0, marker="+", color=ps.INK_SECONDARY, ms=6, zorder=5)
    ax2.annotate("South Pole", (x0, y0), xytext=(5, 6), textcoords="offset points",
                 fontsize=7, color=ps.MUTED, ha="left")
    for text, lon, lat in [("Ross Sea", 180, -73), ("Weddell Sea", -45, -70),
                           ("0°", 0, -61), ("90° E", 90, -62.5), ("90° W", -90, -62.5)]:
        ax2.text(*proj(lon, lat), text, fontsize=7, color=ps.MUTED, ha="center", va="center")
    lim = 2 * np.sin(np.radians(29.2) / 2)
    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)
    ax2.set_aspect("equal")
    ax2.axis("off")
    _title(ax2, "b   Antarctica's meteorite stranding zones",
           "Ice carries buried meteorites towards the mountains; where the flow is blocked, "
           "wind strips\nthe ice away and strands them on the surface. Bubble area ∝ number of "
           "meteorites.", lines=2)
    _caption(fig, f"Lambert azimuthal equal-area projection centred on the South Pole; grid every "
             f"10° of latitude and 30° of longitude. {_fmt(unplaced)} meteorites from Antarctic "
             f"areas with no coordinates in the file (Larkman Nunatak, Roberts Massif, ...) "
             f"cannot be placed. {_provenance(short=True)}")
    return ps.save(fig, FIGURE_DIR / "02_collection_areas")


# --- 03 Witnessed falls follow people ------------------------------------------------------


def figure_witnesses(wtable, results) -> list:
    w = results["witness"]
    fit = w["fit"]
    t = wtable[wtable["population"].notna()].copy()
    fig, ax = plt.subplots(figsize=(9.4, 6.8))
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Individual countries: recessive, sized by area.
    hit = t[t["falls"] > 0]
    size = 8 + 70 * np.sqrt(hit["expected_if_all_found"] / t["expected_if_all_found"].max())
    ax.scatter(hit["density"], hit["efficiency"] * 100, s=size, c=ps.C_FALL, alpha=0.45,
               edgecolors=ps.SURFACE, linewidths=0.6, zorder=3)
    # Big countries with no recorded fall at all: their 95% upper limit.
    zero = t[(t["falls"] == 0) & (t["expected_if_all_found"] >= 1000)]
    ax.scatter(zero["density"], zero["eff_hi"] * 100, marker="v", s=26, facecolors="none",
               edgecolors=ps.INK_SECONDARY, linewidths=0.9, zorder=3)

    # Countries pooled by density: the sturdy points.
    for b in fit["binned"]:
        e = b["efficiency"] * 100
        lo, hi = (c * 100 for c in b["ci95"])
        ax.plot([b["density_typical"]] * 2, [lo, hi], color=ps.INK, lw=1.4, zorder=5)
        ax.scatter(b["density_typical"], e, s=58, marker="D", c=ps.INK, edgecolors=ps.SURFACE,
                   linewidths=1.2, zorder=6)

    d = np.logspace(-1, 3.4, 100)
    g_lo, g_hi = fit["gamma_ci95_bootstrap"]
    curve = lambda g: fit["eps0"] * (d / fit["d_ref"]) ** g * 100  # noqa: E731
    ax.plot(d, curve(fit["gamma"]), color=ps.C_FIT, lw=1.4, zorder=4)
    ax.fill_between(d, np.minimum(curve(g_lo), curve(g_hi)), np.maximum(curve(g_lo), curve(g_hi)),
                    color=ps.C_FIT, alpha=0.08, lw=0, zorder=1)

    labels = {
        "Japan": (12, 10), "China": (14, -12), "India": (14, -10),
        "United States of America": (-14, 12), "Canada": (10, 16), "Brazil": (12, -10),
        "Australia": (-8, -22), "Nigeria": (12, 10), "Czechia": (10, 6), "Mexico": (-10, 22),
        "Russia": (12, -6), "Dem. Rep. Congo": (-12, -12), "Egypt": (-16, 2),
    }
    short = {"United States of America": "USA", "Dem. Rep. Congo": "DR Congo"}
    for name, (dx, dy) in labels.items():
        row = t[t["country"] == name]
        if row.empty:
            continue
        r = row.iloc[0]
        if r["falls"] == 0 and r["expected_if_all_found"] < 1000:
            continue
        yv = (r["efficiency"] if r["falls"] > 0 else r["eff_hi"]) * 100
        txt = short.get(name, name) + (f" {int(r['falls'])}" if r["falls"] > 0 else ": none")
        ps.annotate(ax, txt, xy=(r["density"], yv), xytext=(dx, dy), textcoords="offset points",
                    fontsize=7.3, ha="left" if dx >= 0 else "right", va="center",
                    relpos=(0 if dx >= 0 else 1, 0.5), zorder=7)

    ax.set_xlim(0.2, 3000)
    ax.set_ylim(0.004, 5)
    ax.set_xlabel("people per km² (country average, 1960–2009)")
    ax.set_ylabel("falls of ≥ 100 g that were recorded (%)")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    handles = [
        Line2D([], [], ls="none", marker="o", ms=6, mfc=ps.C_FALL, mec=ps.SURFACE, alpha=0.6),
        Line2D([], [], ls="none", marker="D", ms=6, mfc=ps.INK, mec=ps.SURFACE),
        Line2D([], [], ls="none", marker="v", ms=6, mfc="none", mec=ps.INK_SECONDARY),
        Line2D([], [], color=ps.C_FIT, lw=1.4),
    ]
    ax.legend(handles, ["a country (area ∝ dot area; number = falls recorded)",
                        "countries pooled by density (95% interval)",
                        "large country with no recorded fall (95% upper limit)",
                        f"fit: recovery ∝ density^{fit['gamma']:.2f}"],
              loc="upper left", fontsize=7.6)
    dense = fit["binned"][-1]
    _title(ax, "Witnessed falls follow people",
           f"Recorded falls of ≥ 100 g, {w['years'][0]}–{w['years'][1]}, as a share of the "
           f"{w['reference_rate_per_million_km2_yr']:.0f} per million km² per year that land "
           f"(Halliday et al. 1989).\nWorldwide, 1 fall in {1 / w['global_efficiency']:,.0f} "
           f"is recorded; where more than {dense['density_lo']:.0f} people live per km², "
           f"1 in {1 / dense['efficiency']:,.0f}.", lines=2)
    _caption(fig, f"Poisson fit over {fit['n_countries']} countries: recovery grows as population "
             f"density to the power {fit['gamma']:.2f} (95% bootstrap range "
             f"{g_lo:.2f}–{g_hi:.2f}); where 100 people live per km², "
             f"{fit['eps0'] * 100:.2f}% of falls are recorded. Population: World Bank. "
             f"{_provenance(short=True)}")
    return ps.save(fig, FIGURE_DIR / "03_witnesses")


# --- 04 What we miss --------------------------------------------------------------------------


def figure_what_we_miss(results) -> list:
    rv = results["recovery_vs_mass"]
    m = np.asarray(rv["masses_g"])
    years = rv["years"][1] - rv["years"][0] + 1
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    ax.set_xscale("log")
    ax.set_yscale("log")
    grid = np.logspace(0, 7.5, 400)
    rec = np.array([(m >= g).sum() for g in grid]) / years
    ok = rec > 0
    ax.plot(grid[ok], rec[ok], color=ps.C_FALL, lw=2, drawstyle="steps-post", zorder=4)
    ax.annotate("falls recorded in the catalogue", (2e3, float(rec[np.searchsorted(grid, 2e3)])),
                xytext=(8, 10), textcoords="offset points", fontsize=8, color=ps.INK_SECONDARY)

    land = rv["land_area_ex_antarctica_km2"] / 1e6
    done = set()
    for r in rv["references"]:
        expected = r["expected_per_year_on_land"]
        marker = "o" if r["method"] == "camera network" else "s"
        ax.scatter(r["mass_g"], expected, s=46, marker=marker, c=ps.C_PUBLISHED,
                   edgecolors=ps.SURFACE, linewidths=0.8, zorder=5)
        if r["mass_g"] in done or r["method"] != "camera network":
            continue
        done.add(r["mass_g"])
        ax.annotate("", xy=(r["mass_g"], r["recorded_per_year"] * 1.12),
                    xytext=(r["mass_g"], expected / 1.12),
                    arrowprops=dict(arrowstyle="<->", color=ps.MUTED, lw=0.8), zorder=3)
        ax.text(r["mass_g"] * 1.18, np.sqrt(expected * r["recorded_per_year"]),
                f"1 in {1 / r['efficiency']:,.0f}", fontsize=8.5, color=ps.INK, va="center")

    notes = {
        (10, 222): ("Atacama desert stones,\n2-million-year average\n(Drouard et al. 2019)",
                    -12, 0, "right"),
        (10, 83): ("Canadian camera network\n(Halliday et al. 1989)", -12, 0, "right"),
        (100, 58): ("Halliday et al. 1989", 12, 9, "left"),
        (100, 48): ("European Fireball Network\n(Borovička et al. 2022)", 12, -14, "left"),
        (1000, 9): ("Halliday et al. 1989", 12, 0, "left"),
    }
    for r in rv["references"]:
        k = (r["mass_g"], r["rate"])
        if k in notes:
            text, dx, dy, ha = notes[k]
            ax.annotate(text, (r["mass_g"], r["expected_per_year_on_land"]), xytext=(dx, dy),
                        textcoords="offset points", fontsize=7.2, color=ps.INK_SECONDARY,
                        ha=ha, va="center")
    ax.set_xlim(0.12, 3e7)
    ax.set_ylim(0.01, 1e5)
    ax.set_xlabel("total mass of meteorites the fall delivers, at least")
    ax.set_ylabel("falls per year on all land outside Antarctica")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: {1: "1 g", 10: "10 g", 100: "100 g", 1e3: "1 kg", 1e4: "10 kg",
                      1e5: "100 kg", 1e6: "1 t", 1e7: "10 t"}.get(v, "")))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ref100 = next(r for r in rv["references"] if (r["mass_g"], r["rate"]) == (100, 58))
    unseen = float(f"{1 / ref100['efficiency']:.2g}")
    _title(ax, f"For every fall of 100 g or more that is recorded, about {unseen:,.0f} "
               f"land unseen",
           f"Orange: falls recorded per year, {rv['years'][0]}–{rv['years'][1]} (mass = recovered "
           f"mass). Violet: how often falls reach the ground, measured by\ncamera networks (●) "
           f"and by dating desert stones (■), scaled to the {land:.0f} million km² of land outside "
           f"Antarctica.", lines=2)
    _caption(fig, "The ocean, 71% of the Earth, receives more than twice as many again. Recovered "
             "mass is a lower limit on what fell, which moves the orange curve slightly left of "
             "the truth; the gap is three orders of magnitude either way. "
             + _provenance(short=True))
    return ps.save(fig, FIGURE_DIR / "04_what_we_miss")


# --- 05 What falls vs what is found ------------------------------------------------------------


def figure_what_is_found(results) -> list:
    comp = results["composition"]
    samples = list(analysis.SAMPLES)
    fig = plt.figure(figsize=(12.6, 5.2))
    ax1 = fig.add_axes([0.2, 0.22, 0.43, 0.64])
    ax2 = fig.add_axes([0.7, 0.22, 0.28, 0.64])

    ypos = np.arange(len(samples))[::-1]
    for yy, s in zip(ypos, samples):
        left = 0.0
        for fam in catalog.FAMILIES:
            frac = comp[s]["families"][fam]["fraction"]
            ax1.barh(yy, frac * 100, left=left * 100, height=0.6,
                     color=ps.FAMILY_COLORS[fam], edgecolor=ps.SURFACE, linewidth=1.0, zorder=3)
            if frac >= 0.045:
                ax1.text((left + frac / 2) * 100, yy, f"{frac:.0%}", ha="center", va="center",
                         fontsize=7.2, color="white" if fam in ("ordinary", "iron") else ps.INK,
                         zorder=4)
            left += frac
    names = [f"{analysis.SAMPLE_NAMES[s]}\n{_fmt(comp[s]['n'])} classified" for s in samples]
    ax1.set_yticks(ypos, names, fontsize=8)
    ax1.tick_params(axis="y", length=0, colors=ps.INK_SECONDARY)
    ax1.set_xlim(0, 100)
    ax1.grid(axis="y", visible=False)
    ax1.set_xlabel("share of the collection (%)")
    _title(ax1, "a   What each kind of search brings home")
    handles = [plt.Rectangle((0, 0), 1, 1, color=ps.FAMILY_COLORS[f]) for f in catalog.FAMILIES]
    ax1.legend(handles, [catalog.FAMILY_NAMES[f] for f in catalog.FAMILIES], ncol=3,
               loc="upper left", bbox_to_anchor=(0.0, -0.15), fontsize=7.8)

    # (b) Irons and stony-irons, on a log scale with 95% intervals.
    ref = comp["fall"]["iron_or_stony_iron"]["fraction"] * 100
    ax2.axvline(ref, color=ps.C_FALL, lw=1.0, zorder=1)
    ax2.text(ref * 1.08, len(samples) - 0.45, "what falls", color=ps.INK_SECONDARY,
             fontsize=7.5)
    for yy, s in zip(ypos, samples):
        p = comp[s]["iron_or_stony_iron"]
        v = p["fraction"] * 100
        lo, hi = (c * 100 for c in p["ci95"])
        ax2.plot([lo, hi], [yy, yy], color=ps.INK_SECONDARY, lw=1.2, zorder=2)
        ax2.scatter(v, yy, s=36, c=ps.FAMILY_COLORS["iron"], edgecolors=ps.SURFACE,
                    linewidths=1.2, zorder=3)
        ax2.text(hi * 1.15, yy, f"{v:.1f}%", va="center", fontsize=7.8, color=ps.INK)
    ax2.set_xscale("log")
    ax2.set_xlim(0.1, 100)
    ax2.set_ylim(-0.6, len(samples) - 0.2)
    ax2.set_yticks(ypos, [])
    ax2.tick_params(axis="y", length=0)
    ax2.grid(axis="y", visible=False)
    ax2.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}%"))
    ax2.set_xlabel("share that is metal (95% interval)")
    _title(ax2, "b   Irons and stony-irons")
    _caption(fig, "Witnessed falls are the fair sample of what reaches the ground. Iron meteorites "
             "are dense, magnetic, unlike any ordinary rock and outlast stones in wet climates, so "
             "people who stumble on them keep them; on Antarctic ice they may sink out of sight "
             "(Evatt et al. 2016); dealers get the rare stony types classified first. Counts are "
             "of stones: a shower that scattered many stones counts many times. "
             + _provenance(short=True))
    return ps.save(fig, FIGURE_DIR / "05_what_is_found")


# --- 06 Masses ------------------------------------------------------------------------------------


def figure_masses(df, results) -> list:
    ms = results["masses"]
    fig = plt.figure(figsize=(12.4, 5.5))
    ax1 = fig.add_axes([0.06, 0.12, 0.46, 0.74])
    ax2 = fig.add_axes([0.665, 0.12, 0.315, 0.74])

    v = df[df["valid"]]
    handles, names = [], []
    for env in ("antarctic", "desert", "fall", "chance"):
        m = np.sort(v.loc[v["environment"].eq(env), "mass_g"].dropna().to_numpy())
        cdf = np.arange(1, len(m) + 1) / len(m)
        (line,) = ax1.plot(m, cdf * 100, color=ps.ENV_COLORS[env], lw=2, zorder=3)
        med = ms["median_g"][env]["median"]
        ax1.scatter([med], [50], s=40, c=ps.ENV_COLORS[env], edgecolors=ps.SURFACE,
                    linewidths=1.2, zorder=4)
        handles.append(line)
        names.append(f"{analysis.ENV_NAMES[env]}: median {_mass(med)}")
    ax1.legend(handles, names, loc="lower right", fontsize=8)
    ax1.set_xscale("log")
    ax1.set_xlim(0.1, 1e8)
    ax1.set_ylim(0, 100)
    ax1.axhline(50, color=ps.AXIS, lw=0.6, zorder=1)
    ax1.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: _mass(x)))
    ax1.set_xlabel("recovered mass")
    ax1.set_ylabel("share of meteorites lighter than this (%)")
    _title(ax1, "a   Antarctic teams pick up pebbles; passers-by pick up boulders")

    heavy = ms["heaviest"][:10]
    stone = ms["heaviest_stones"][0]
    rows = heavy + [stone]
    ypos = np.arange(len(rows))[::-1].astype(float)
    ypos[-1] -= 0.6
    cols = [ps.FAMILY_COLORS["iron"]] * len(heavy) + [ps.FAMILY_COLORS["ordinary"]]
    ax2.barh(ypos, [r["mass_g"] / 1e6 for r in rows], height=0.62, color=cols, zorder=3)
    for yy, r in zip(ypos, rows):
        how = "fell" if r["fall"] == "Fell" else "found"
        ax2.text(r["mass_g"] / 1e6 + 0.8, yy, f"{r['mass_g'] / 1e6:g} t · {how} {int(r['year'])}",
                 va="center", fontsize=7.2, color=ps.INK_SECONDARY)
    names = [f"{r['name']} ({r.get('country') or '?'})".replace("United States of America", "USA")
             for r in rows]
    ax2.set_yticks(ypos, names, fontsize=7.8)
    ax2.tick_params(axis="y", length=0, colors=ps.INK_SECONDARY)
    ax2.grid(axis="y", visible=False)
    ax2.set_xlim(0, 88)
    ax2.set_xlabel("tonnes")
    _title(ax2, "b   The heavyweights are all iron")
    ax2.text(0.98, 0.02, "heaviest stone,\nfor scale", transform=ax2.transAxes, ha="right",
             va="bottom", fontsize=7.2, color=ps.MUTED)
    share = ms["mass_share_by_family"]["iron"]
    iron_falls = results["composition"]["fall"]["families"]["iron"]["fraction"]
    _caption(fig, f"Irons are {iron_falls:.1%} of witnessed falls but {share:.0%} of the "
             f"{ms['total_mass_kg'] / 1e3:,.0f} tonnes of meteorites in the catalogue. Mass is the "
             f"total recovered ('total known weight'); the curves in (a) are cumulative, with the "
             f"median marked. {_provenance(short=True)}")
    return ps.save(fig, FIGURE_DIR / "06_masses")


# --- 07 History --------------------------------------------------------------------------------


def figure_history(df, results) -> list:
    h = results["history"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8.4), sharex=True,
                                   gridspec_kw={"height_ratios": [1.5, 1], "hspace": 0.2})
    fig.subplots_adjust(top=0.93, bottom=0.07)
    years = np.arange(1750, 2014)
    per = {int(k): v for k, v in h["per_year"].items()}
    ends = {}
    for env in ("fall", "chance", "antarctic", "desert"):
        before = sum(v.get(env, 0) for y, v in per.items() if y < years[0])
        cum = before + np.cumsum([per.get(y, {}).get(env, 0) for y in years])
        ok = cum > 0
        ax1.plot(years[ok], cum[ok], color=ps.ENV_COLORS[env], lw=2, zorder=3)
        ends[env] = cum[-1]
    ax1.set_yscale("log")
    ax1.set_ylim(0.8, 8e4)
    ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    for env in ("antarctic", "desert", "chance", "fall"):
        ax1.annotate(f"{analysis.ENV_NAMES[env]}  {_fmt(ends[env])}", (2013, ends[env]),
                     xytext=(6, 0), textcoords="offset points", va="center", fontsize=8,
                     color=ps.INK_SECONDARY, annotation_clip=False)
    placement = {1803: ("right", 4e4, "top"), 1969: ("right", 4e4, "top"),
                 1976: ("left", 1.6, "bottom")}
    for yr, text in MILESTONES:
        ha, yy, va = placement.get(yr, ("right", 4e4, "top"))
        ax1.axvline(yr, color=ps.AXIS, lw=0.7, zorder=1)
        ax1.text(yr + (2 if ha == "left" else -2), yy, text, ha=ha, va=va, fontsize=7.2,
                 color=ps.INK_SECONDARY, bbox=LABEL_BOX, zorder=5)
    first = h["first_antarctic_find"]
    ax1.annotate(f"{first['name']}, {first['year']}:\nfirst meteorite found in Antarctica",
                 (first["year"], 1), xytext=(-14, 30), textcoords="offset points", ha="right",
                 fontsize=7.2, color=ps.INK_SECONDARY,
                 arrowprops=dict(arrowstyle="-", color=ps.INK_SECONDARY, lw=0.7))
    ax1.set_ylabel("meteorites known (cumulative)")
    _title(ax1, "How the catalogue grew",
           f"{h['share_of_finds_since_1969']:.0%} of all finds were made after 1969, when "
           f"Japanese glaciologists found meteorites concentrated on the ice.")

    dec = {int(k): v for k, v in h["falls_per_decade"].items() if 1750 <= int(k) <= 2000}
    ax2.bar(np.array(list(dec)) + 5, list(dec.values()), width=8, color=ps.C_FALL, zorder=3)
    for d0, n in dec.items():
        if d0 % 50 == 0 or d0 == 1930:
            ax2.text(d0 + 5, n + 2, str(n), ha="center", fontsize=7.2, color=ps.INK_SECONDARY)
    ax2.set_ylabel("witnessed falls per decade")
    ax2.set_xlim(1750, 2048)
    wp = h["world_population"]
    _title(ax2, "Witnessed falls: steady at about six a year since 1900",
           f"1800s: {h['falls_per_year_1800_1899']:.1f} per year; 1900–2009: "
           f"{h['falls_per_year_1900_2009']:.1f} per year. The 1960s gave {dec[1960]} and the "
           f"2000s {dec[2000]}, while the world's population grew {wp['ratio']:.1f}-fold.")
    _caption(fig, f"Year = year of fall or find (for Antarctic finds, the field season); "
             f"{_fmt(h['undated'])} meteorites without a year are not counted. The last few years "
             f"are incomplete: meteorites are catalogued only once classified, often years after "
             f"collection. {_provenance(short=True)}")
    return ps.save(fig, FIGURE_DIR / "07_history")


# --- All -------------------------------------------------------------------------------------------


def make_all(df, wtable, countries, results) -> list:
    written = []
    written += figure_world_map(df, countries, results)
    written += figure_collection_areas(df, countries, results)
    written += figure_witnesses(wtable, results)
    written += figure_what_we_miss(results)
    written += figure_what_is_found(results)
    written += figure_masses(df, results)
    written += figure_history(df, results)
    return written
