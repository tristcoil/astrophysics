"""The interactive 3D map: a self-contained HTML page you can rotate and fly into.

    python map3d.py            # after run.py has produced the catalogue

What is in it:

* Every star within 20 pc, named from SIMBAD in the hover text, plus a random
  one in six of the stars between 20 and 100 pc. The subsample is uniform --
  the same fraction of every type -- so the map still shows honestly that red
  dwarfs outnumber everything else. One legend entry per type: click to hide.
* Famous stars labelled; famous stars *beyond* the map as signposts on its
  edge, pointing the right way and giving their real distance.
* The two star clusters the analysis found, and (off by default) every known
  planet host.
* Camera presets: 3D, from above the Galactic plane, edge-on, 20 pc and 5 pc.

Star colour is the blackbody colour of each star's temperature, the same as in
the static figures; the colour bar is the key. The page embeds plotly.js, so it
works offline and can be shared as a single file.
"""

from __future__ import annotations

import html
import json

import numpy as np
import plotly.graph_objects as go

import plotstyle as ps
import stellar_physics as sp
from analysis import CLASS_LABELS
from config import MAP_HTML_PATH, MAX_DISTANCE_PC, NAMED_RADIUS_PC

#: One in this many stars beyond the named radius is drawn. Uniform across
#: types, so relative numbers stay true; sparse enough that the 100 pc ball
#: reads as a volume you can see into rather than a solid surface.
SUBSAMPLE = 6
SEED = 7

#: Temperature range spanned by the colour scale (log10 K).
LOG_T_MIN, LOG_T_MAX = np.log10(2300.0), np.log10(12000.0)

SCENE_BG = ps.PAGE
LAYER_ORDER = ["M", "K", "G", "F", "A", "B", "Giant", "WD", "L/T", "No colour"]
PLURAL = {"No colour": "Unclassified (no colour)"}


def _colourscale(n: int = 32) -> list:
    colours = sp.StarColours()
    lt = np.linspace(LOG_T_MIN, LOG_T_MAX, n)
    hexes = sp.srgb_to_hex(colours(10 ** lt))
    return [[float(i / (n - 1)), h] for i, h in enumerate(hexes)]


def _marker_size(m_g) -> np.ndarray:
    """Screen size in pixels, growing with luminosity (plotly sizes are pixels, not areas)."""
    m_g = np.where(np.isfinite(m_g), m_g, 14.0)
    return np.clip(8.0 - 0.5 * m_g, 1.0, 10.0)


def _star_label(row) -> str:
    """Rich hover text, used only for the stars within the named radius."""
    name = sp.tidy_simbad_name(str(row["name"])) if str(row["name"]).strip() else ""
    main = sp.tidy_simbad_name(str(row["simbad_id"])) if str(row["simbad_id"]).strip() else ""
    title = name or main or (f"Gaia DR3 {row['source_id']}" if row["origin"] == "Gaia DR3" else f"HIP {row['hip']}")
    lines = [f"<b>{html.escape(title)}</b>"]
    if main and main != title:
        lines.append(html.escape(main))
    spt = str(row["simbad_sptype"]).strip()
    if spt:
        lines.append(f"spectral type {html.escape(spt)} (SIMBAD)")
    src = f"Gaia DR3 {row['source_id']}" if row["origin"] == "Gaia DR3" else f"Hipparcos HIP {row['hip']}"
    lines.append(html.escape(src))
    return "<br>".join(lines) + "<br>"


def build_figure(cat, res) -> go.Figure:
    rng = np.random.default_rng(SEED)
    d = np.asarray(cat["dist_pc"], dtype=float)
    cls = np.asarray(cat["cls"]).astype(str)
    near = d <= NAMED_RADIUS_PC
    keep = near | (rng.random(len(d)) < 1.0 / SUBSAMPLE)
    teff = np.asarray(cat["display_teff"], dtype=float)
    logt = np.clip(np.log10(np.where(np.isfinite(teff), teff, 5800.0)), LOG_T_MIN, LOG_T_MAX)
    cscale = _colourscale()

    fig = go.Figure()
    for c in LAYER_ORDER:
        sel = np.flatnonzero(keep & (cls == c))
        if len(sel) == 0:
            continue
        text = np.array([""] * len(sel), dtype=object)
        for j, i in enumerate(sel):
            if near[i]:
                text[j] = _star_label(cat[i])
        # Sizes (whole pixels) and colours (256 steps along the temperature
        # scale) go in as bytes: indistinguishable on screen, a quarter the size.
        marker = dict(size=np.round(_marker_size(np.asarray(cat["m_g"], dtype=float)[sel])).astype(np.uint8),
                      opacity=0.5, line=dict(width=0))
        if c == "No colour":
            marker["color"] = ps.MUTED
        else:
            level = np.round(255 * (logt[sel] - LOG_T_MIN) / (LOG_T_MAX - LOG_T_MIN)).astype(np.uint8)
            marker.update(color=level, colorscale=cscale, cmin=0, cmax=255, showscale=False)
        label = PLURAL.get(c, CLASS_LABELS[c])
        n_all = int((cls == c).sum())
        fig.add_trace(go.Scatter3d(
            x=np.asarray(cat["x"])[sel].astype(np.float32), y=np.asarray(cat["y"])[sel].astype(np.float32),
            z=np.asarray(cat["z"])[sel].astype(np.float32), mode="markers", marker=marker,
            name=f"{label} ({n_all:,})", text=list(text),
            customdata=np.column_stack([d[sel], d[sel] * sp.LY_PER_PC]).astype(np.float32),
            hovertemplate=("%{text}" + f"<span style='color:{ps.INK_SECONDARY}'>{html.escape(label.rstrip('s'))}</span>"
                           "<br>%{customdata[0]:.2f} pc  ·  %{customdata[1]:.1f} light years<extra></extra>"),
        ))

    # The Sun, carrying the colour key so it never disappears with a layer.
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0], mode="markers+text", text=["Sun"], textposition="top center",
        textfont=dict(color=ps.INK, size=12),
        marker=dict(size=9, color=[np.log10(5772.0)], colorscale=cscale, cmin=LOG_T_MIN, cmax=LOG_T_MAX,
                    line=dict(width=0), showscale=True,
                    colorbar=dict(title=dict(text="surface<br>temperature", font=dict(color=ps.INK_SECONDARY, size=11)),
                                  tickvals=np.log10([3000, 4000, 5000, 5772, 7500, 10000]).tolist(),
                                  ticktext=["3,000 K", "4,000 K", "5,000 K", "Sun", "7,500 K", "10,000 K"],
                                  tickfont=dict(color=ps.INK_SECONDARY, size=10), len=0.45, thickness=12,
                                  x=0.0, xanchor="left", y=0.5, outlinewidth=0)),
        name="the Sun", hovertemplate="<b>The Sun</b><br>you are here<extra></extra>", showlegend=False))

    # Famous stars inside the map, labelled.
    inside = [f for f in res["famous"]["inside"] if f.get("xyz_pc")]
    fx, fy, fz = (np.array([f["xyz_pc"][k] for f in inside]) for k in range(3))
    hover = []
    for f in inside:
        src = f["source"] if f["source"] in ("Gaia DR3", "Hipparcos") else f"position from {f['source']}"
        extra = " (not in our catalogue)" if f.get("catalogue_index") is None else ""
        hover.append(f"<b>{html.escape(f['label'])}</b><br>{html.escape(f['note'])}<br>"
                     f"{f['dist_pc']:.2f} pc · {f['dist_pc'] * sp.LY_PER_PC:.1f} light years<br>"
                     f"{html.escape(str(f.get('sp_type') or ''))} · {html.escape(src)}{extra}")
    fig.add_trace(go.Scatter3d(
        x=fx, y=fy, z=fz, mode="markers+text", text=[f["label"] for f in inside], textposition="top center",
        textfont=dict(color=ps.INK_SECONDARY, size=11),
        marker=dict(size=7, symbol="circle-open", color=ps.INK, line=dict(width=1.5)),
        hovertext=hover, hovertemplate="%{hovertext}<extra></extra>", name="Famous stars (labelled)"))

    # Famous stars beyond the map: signposts on its edge, in their true directions.
    beyond = [f for f in res["famous"]["beyond"] if f.get("xyz_pc")]
    lx, ly, lz, tx, ty, tz, labels, hov = [], [], [], [], [], [], [], []
    for f in beyond:
        u = np.array(f["xyz_pc"]) / f["dist_pc"]
        a, b, t = u * MAX_DISTANCE_PC, u * (MAX_DISTANCE_PC + 12), u * (MAX_DISTANCE_PC + 17)
        lx += [a[0], b[0], None]; ly += [a[1], b[1], None]; lz += [a[2], b[2], None]
        tx.append(t[0]); ty.append(t[1]); tz.append(t[2])
        labels.append(f"{f['label']} → {f['dist_pc']:.0f} pc")
        hov.append(f"<b>{html.escape(f['label'])}</b><br>{html.escape(f['note'])}<br>really at "
                   f"{f['dist_pc']:.0f} pc ({f['dist_pc'] * sp.LY_PER_PC:.0f} ly), beyond this map<br>"
                   f"parallax: {html.escape(f['source'])}")
    fig.add_trace(go.Scatter3d(x=lx, y=ly, z=lz, mode="lines", line=dict(color=ps.MUTED, width=3),
                               hoverinfo="skip", name="Famous stars beyond the map", legendgroup="beyond"))
    fig.add_trace(go.Scatter3d(x=tx, y=ty, z=tz, mode="text", text=labels, textfont=dict(color=ps.MUTED, size=10),
                               hovertext=hov, hovertemplate="%{hovertext}<extra></extra>", showlegend=False,
                               legendgroup="beyond", name="Famous stars beyond the map"))

    # The clusters found by the analysis.
    peaks = [p for p in res["structure"]["peaks"] if p.get("match")]
    fig.add_trace(go.Scatter3d(
        x=[p["centre_xyz_pc"][0] for p in peaks], y=[p["centre_xyz_pc"][1] for p in peaks],
        z=[p["centre_xyz_pc"][2] for p in peaks], mode="markers+text",
        text=[p["match"]["name"] for p in peaks], textposition="bottom center",
        textfont=dict(color=ps.INK, size=12),
        marker=dict(size=26, symbol="circle-open", color=ps.INK_SECONDARY, line=dict(width=2)),
        hovertext=[(f"<b>{p['match']['name']}</b> star cluster<br>found by this analysis at {p['distance_pc']:.1f} pc"
                    f"<br>literature: {p['match']['literature_distance_pc']:.1f} pc<br>"
                    f"~{p['excess_within_mask_radius']:.0f} stars above the field within 10 pc") for p in peaks],
        hovertemplate="%{hovertext}<extra></extra>", name="Star clusters found"))

    # Known planet hosts (off until switched on in the legend).
    hosts = res["planet_hosts"]
    fig.add_trace(go.Scatter3d(
        x=[h["xyz_pc"][0] for h in hosts], y=[h["xyz_pc"][1] for h in hosts], z=[h["xyz_pc"][2] for h in hosts],
        mode="markers", marker=dict(size=4.5, symbol="diamond-open", color=ps.AQUA, line=dict(width=1.2)),
        hovertext=[(f"<b>{html.escape(h['host'])}</b><br>{h['n_planets']} known planet{'s' if h['n_planets'] > 1 else ''}: "
                    f"{html.escape(', '.join(h['planets'][:6]))}{' ...' if h['n_planets'] > 6 else ''}<br>"
                    f"found by {html.escape(', '.join(h['methods']))}<br>{h['dist_pc']:.1f} pc") for h in hosts],
        hovertemplate="%{hovertext}<extra></extra>", visible="legendonly",
        name=f"Known planet hosts ({len(hosts):,})"))

    # Orientation: distance rings in the Galactic plane, and the three axes.
    ring_x, ring_y, ring_z = [], [], []
    th = np.linspace(0, 2 * np.pi, 121)
    for r in (10, 20, 50, 100):
        ring_x += list(r * np.cos(th)) + [None]; ring_y += list(r * np.sin(th)) + [None]; ring_z += [0] * 121 + [None]
    fig.add_trace(go.Scatter3d(x=ring_x, y=ring_y, z=ring_z, mode="lines", line=dict(color=ps.AXIS, width=2),
                               hoverinfo="skip", name="Distance rings, 10-100 pc"))
    axes = [((125, 0, 0), "Galactic centre →<br>(8,200 pc away)"), ((0, 125, 0), "direction of<br>Galactic rotation"),
            ((0, 0, 125), "North Galactic Pole")]
    fig.add_trace(go.Scatter3d(
        x=sum(([0, a[0][0], None] for a in axes), []), y=sum(([0, a[0][1], None] for a in axes), []),
        z=sum(([0, a[0][2], None] for a in axes), []), mode="lines", line=dict(color=ps.MUTED, width=2),
        hoverinfo="skip", showlegend=False, legendgroup="rings"))
    fig.add_trace(go.Scatter3d(
        x=[a[0][0] * 1.06 for a in axes], y=[a[0][1] * 1.06 for a in axes], z=[a[0][2] * 1.06 for a in axes],
        mode="text", text=[a[1] for a in axes], textfont=dict(color=ps.MUTED, size=10), hoverinfo="skip",
        showlegend=False))

    menu = _camera_menu()
    axis = dict(showbackground=False, showgrid=True, gridcolor=ps.GRID, zeroline=False, showspikes=False,
                color=ps.MUTED, tickfont=dict(size=9, color=ps.MUTED), range=[-135, 135])
    fig.update_layout(
        paper_bgcolor=SCENE_BG, plot_bgcolor=SCENE_BG, font=dict(family="system-ui, -apple-system, 'Segoe UI', sans-serif",
                                                                color=ps.INK_SECONDARY),
        margin=dict(l=0, r=0, t=0, b=0), hoverlabel=dict(bgcolor=ps.SURFACE, bordercolor=ps.AXIS,
                                                         font=dict(color=ps.INK, size=12)),
        legend=dict(x=0.995, xanchor="right", y=0.98, bgcolor="rgba(13,13,13,0.7)", font=dict(size=11, color=ps.INK_SECONDARY),
                    itemsizing="constant", title=dict(text="click to show / hide", font=dict(size=10, color=ps.MUTED))),
        scene=dict(xaxis=dict(axis, title="X (pc)"), yaxis=dict(axis, title="Y (pc)"),
                   zaxis=dict(axis, title="Z (pc)"), aspectmode="cube", bgcolor=SCENE_BG,
                   camera=dict(eye=dict(x=1.35, y=1.1, z=0.55))),
        # Two copies of the camera presets: a button row for wide screens and a
        # dropdown for phones. The page script shows whichever fits.
        updatemenus=[dict(menu, type=kind, direction=direction, visible=visible)
                     for kind, direction, visible in (("buttons", "right", True), ("dropdown", "down", False))],
    )
    return fig


def _camera_menu() -> dict:
    return dict(
        x=0.01, xanchor="left", y=0.99, yanchor="top", showactive=True,
        bgcolor=ps.SURFACE, bordercolor=ps.AXIS, font=dict(color=ps.INK_SECONDARY, size=11), pad=dict(r=4, t=4),
        buttons=[
            dict(label="3D", method="relayout", args=[_view(135, (1.35, 1.1, 0.55))]),
            dict(label="From above", method="relayout", args=[_view(135, (0.0, 0.0, 2.1), up=(0, 1, 0))]),
            dict(label="Edge-on", method="relayout", args=[_view(135, (0.0, -2.1, 0.0))]),
            dict(label="Within 20 pc", method="relayout", args=[_view(24, (1.35, 1.1, 0.55))]),
            dict(label="Within 5 pc", method="relayout", args=[_view(6, (1.35, 1.1, 0.55))]),
        ])


def _view(half_range: float, eye, up=(0, 0, 1)) -> dict:
    r = [-half_range, half_range]
    return {"scene.xaxis.range": r, "scene.yaxis.range": r, "scene.zaxis.range": r,
            "scene.camera": {"eye": dict(zip("xyz", eye)), "up": dict(zip("xyz", up)),
                             "center": {"x": 0, "y": 0, "z": 0}}}


#: plotly.js version that plotly.py embeds; the artifact loads exactly this build.
PLOTLY_CDN = "https://cdnjs.cloudflare.com/ajax/libs/plotly.js/{version}/plotly.min.js"
FONTS_CSS = ("https://fonts.googleapis.com/css2?family=IBM+Plex+Sans+Condensed:wght@500;600"
             "&family=IBM+Plex+Sans:wght@400;500&display=swap")

TITLE = "Stars Within 100 Parsecs"

STYLE = """<title>{title}</title>
<style>
  /* One deliberate dark look: a star map is drawn on the night sky, and hot stars'
     true colours (pale blue-white) need a dark ground. Header band, full-bleed map, notes. */
  :root {{
    color-scheme: dark;
    --page: {page}; --surface: {surface}; --ink: {ink}; --ink-2: {ink2}; --muted: {muted}; --rule: {rule};
    --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
    --body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  }}
  * {{ box-sizing: border-box; }}
  html, body {{ background: var(--page); }}
  body {{ margin: 0; color: var(--ink-2); font: 15px/1.55 var(--body); }}
  .band {{ max-width: 1040px; margin: 0 auto; padding-inline: 16px; padding-block: 22px 10px; }}
  h1 {{ font: 600 clamp(1.45rem, 3.4vw, 2.1rem)/1.15 var(--display); color: var(--ink); margin: 0 0 8px;
        letter-spacing: 0.005em; text-wrap: balance; }}
  .lede {{ margin: 0; max-width: 68ch; }}
  .facts {{ display: flex; flex-wrap: wrap; gap: 10px 32px; margin: 16px 0 0; padding: 0; list-style: none; }}
  .facts li {{ font-size: 0.82rem; color: var(--muted); max-width: 30ch; min-width: 0; }}
  .facts b {{ display: block; font: 600 1.3rem/1.2 var(--display); color: var(--ink);
              font-variant-numeric: tabular-nums; }}
  #map {{ width: 100%; height: 80vh; min-height: 460px; }}
  .notes {{ font-size: 0.8rem; color: var(--muted); padding-block: 12px 28px; }}
  .notes p {{ margin: 0 0 8px; max-width: 92ch; }}
  @media (max-width: 720px) {{
    #map {{ height: 82vh; min-height: 460px; }}
    .facts {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 10px 18px; }}
  }}
</style>
"""

BODY = """<header class="band">
  <h1>Every star within 100 parsecs</h1>
  <p class="lede">Our corner of the Milky Way as mapped by ESA's Gaia mission: {n_total} stars within 326 light
     years of the Sun. Drag to rotate, scroll or pinch to zoom, hover or tap a star for its name. Each star is drawn
     in its real colour, set by its surface temperature.</p>
  <ul class="facts">
    <li><b>{m_share}%</b>of the stars within 20 pc are red dwarfs</li>
    <li><b>0</b>of those red dwarfs can be seen with the naked eye</li>
    <li><b>{hyades} pc</b>to the Hyades, the nearest star cluster, found blind in this data</li>
    <li><b>{n_planet_hosts}</b>stars here with known planets: switch them on in the legend</li>
  </ul>
</header>
<div id="map">{plot}</div>
<footer class="band notes">
  <p>Shown: every star within {named} pc, with its SIMBAD name, and a random one in {sub} of the stars beyond, the
     same fraction of every type so the mix stays true. Famous stars beyond 100 pc appear as signposts on the edge.
     Axes are Galactic coordinates centred on the Sun, in parsecs (1 pc = 3.26 light years).</p>
  <p>Data: ESA Gaia DR3, with spurious solutions removed by the Rybizki et al. (2022) classifier; stars too bright
     for Gaia from Hipparcos (van Leeuwen 2007); names and spectral types from SIMBAD (CDS, Strasbourg); planet
     hosts from the NASA Exoplanet Archive. Project B2 of the Astrophysics Project Atlas.</p>
</footer>
<script>
  /* At phone width, move the legend under the map and fold the view buttons into a menu. */
  (function () {{
    var gd = document.querySelector("#map .plotly-graph-div");
    if (!gd || !window.Plotly) return;
    var last = null;
    var sun = gd.data.findIndex(function (t) {{ return t.name === "the Sun"; }});
    function fit() {{
      var narrow = window.innerWidth < 720;
      if (narrow === last) return;
      last = narrow;
      window.Plotly.relayout(gd, narrow ? {{
        "legend.orientation": "h", "legend.x": 0, "legend.xanchor": "left", "legend.y": 0.27,
        "legend.yanchor": "top", "legend.font.size": 10, "legend.entrywidth": 0.5,
        "legend.entrywidthmode": "fraction", "scene.domain.y": [0.3, 1],
        "updatemenus[0].visible": false, "updatemenus[1].visible": true
      }} : {{
        "legend.orientation": "v", "legend.x": 0.995, "legend.xanchor": "right", "legend.y": 0.98,
        "legend.yanchor": "top", "legend.font.size": 11, "legend.entrywidth": 0,
        "legend.entrywidthmode": "pixels", "scene.domain.y": [0, 1],
        "updatemenus[0].visible": true, "updatemenus[1].visible": false
      }});
      if (sun >= 0) window.Plotly.restyle(gd, narrow
        ? {{ "marker.colorbar.len": 0.2, "marker.colorbar.y": 0.9, "marker.colorbar.thickness": 8,
             "marker.colorbar.tickfont.size": 8, "marker.colorbar.title.font.size": 9 }}
        : {{ "marker.colorbar.len": 0.45, "marker.colorbar.y": 0.5, "marker.colorbar.thickness": 12,
             "marker.colorbar.tickfont.size": 10, "marker.colorbar.title.font.size": 11 }}, [sun]);
    }}
    var t = null;
    window.addEventListener("resize", function () {{ clearTimeout(t); t = setTimeout(fit, 150); }});
    if (gd.layout) fit(); else gd.addEventListener("plotly_afterplot", fit, {{ once: true }});
  }})();
</script>
"""


def _page_parts(cat, res, plot: str) -> tuple[str, str]:
    hy = next((p for p in res["structure"]["peaks"] if (p.get("match") or {}).get("name") == "Hyades"), None)
    tokens = dict(page=ps.PAGE, surface=ps.SURFACE, ink=ps.INK, ink2=ps.INK_SECONDARY, muted=ps.MUTED, rule=ps.AXIS)
    style = STYLE.format(title=TITLE, **tokens)
    body = BODY.format(
        n_total=f"{res['sample']['n_total']:,}", m_share=f"{100 * res['census']['m_dwarf_share_of_stars_20pc']:.0f}",
        hyades=f"{hy['distance_pc']:.1f}" if hy else "?", n_planet_hosts=f"{len(res['planet_hosts']):,}",
        named=f"{NAMED_RADIUS_PC:g}", sub=SUBSAMPLE, plot=plot)
    return style, body


PLOT_CONFIG = {"displaylogo": False, "responsive": True, "modeBarButtonsToRemove": ["toImage"]}


def make_map(cat, res, *, path=MAP_HTML_PATH) -> str:
    """The self-contained page for the repository: one file, plotly.js embedded, works offline."""
    fig = build_figure(cat, res)
    plot = fig.to_html(full_html=False, include_plotlyjs=True, default_height="100%",
                       default_width="100%", config=PLOT_CONFIG)
    style, body = _page_parts(cat, res, plot)
    page = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
            f"{style}</head>\n<body>\n{body}</body>\n</html>\n")
    path.write_text(page, encoding="utf-8")
    return str(path)


def make_artifact_page(cat, res, path) -> str:
    """The same page shaped for publishing as a claude.ai Artifact.

    The Artifact host supplies the document skeleton, so this is body content
    only, and plotly.js loads from cdnjs (pinned to the exact version plotly.py
    embeds) instead of being inlined.
    """
    from plotly.offline import get_plotlyjs_version

    fig = build_figure(cat, res)
    plot = fig.to_html(full_html=False, include_plotlyjs=False, default_height="100%",
                       default_width="100%", config=PLOT_CONFIG)
    style, body = _page_parts(cat, res, plot)
    cdn = PLOTLY_CDN.format(version=get_plotlyjs_version())
    page = (f"{style}<link rel=\"stylesheet\" href=\"{FONTS_CSS}\">\n"
            f"<script src=\"{cdn}\"></script>\n{body}")
    path.write_text(page, encoding="utf-8")
    return str(path)


if __name__ == "__main__":
    import config as C
    from analysis import load_catalogue

    catalogue = load_catalogue()
    results = json.loads(C.RESULTS_PATH.read_text())
    out = make_map(catalogue, results)
    print(f"wrote {out} ({MAP_HTML_PATH.stat().st_size / 1e6:.1f} MB)")
