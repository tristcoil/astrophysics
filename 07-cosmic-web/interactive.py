"""The interactive 3D map: a self-contained HTML page.

135,000 galaxies of the volume-limited sample -- the fair census, complete to
z = 0.106 -- in true 3D comoving positions, rotatable and zoomable in any
browser. The page embeds plotly.js and the data, so it works offline and can
be shared as a single file. It follows the viewer's light/dark setting, with
colours re-stepped for each surface rather than inverted.

Seen all at once, 450 Mpc of galaxies stacked in projection wash the web out,
so the galaxies are also split into declination slabs a few degrees thick,
selectable from a menu: one slab at a time shows walls and voids the way the
static slice does, and can still be turned around in 3D.
"""

from __future__ import annotations

import json

import numpy as np

import cosmo_physics as cp
import plotstyle as ps
from config import FIGURE_DIR, SLICE_DEC_HALF_WIDTH, SLICE_RA_RANGE, VL_Z_MAX, VL_Z_MIN

THEMES = {
    "light": dict(surface=ps.SURFACE, ink=ps.INK, ink2=ps.INK_SECONDARY, muted=ps.MUTED,
                  red=ps.C_RED_GALAXY, blue=ps.C_BLUE_GALAXY, label=ps.INK,
                  slice="rgba(137,135,129,0.18)", menu_bg="#f0efec", menu_border="#c3c2b7"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                 red="#e66767", blue="#3987e5", label="#ffffff",
                 slice="rgba(195,194,183,0.14)", menu_bg="#2c2c2a", menu_border="#383835"),
}

#: Declination slabs, degrees. The second is exactly the slice of Fig. 01.
SLABS = [(-12.0, -SLICE_DEC_HALF_WIDTH), (-SLICE_DEC_HALF_WIDTH, SLICE_DEC_HALF_WIDTH),
         (SLICE_DEC_HALF_WIDTH, 12.0), (12.0, 22.0), (22.0, 32.0), (32.0, 42.0),
         (42.0, 52.0), (52.0, 72.0)]

#: Marker opacity with everything shown (dense) and with one slab shown.
OPACITY_ALL, OPACITY_ONE = 0.45, 0.9


def _slice_mesh(n=60):
    """A fan of triangles marking the equatorial slice of the static maps."""
    d = float(cp.comoving_distance(VL_Z_MAX))
    ra = np.radians(np.linspace(*SLICE_RA_RANGE, n))
    x = np.concatenate([[0.0], d * np.cos(ra)])
    y = np.concatenate([[0.0], d * np.sin(ra)])
    return x, y, np.zeros_like(x), np.zeros(n - 1, int), np.arange(1, n), np.arange(2, n + 1)


def _colorscale(th):
    return [[0.0, th["blue"]], [1.0, th["red"]]]


def cosmic_web_3d(res, art) -> str:
    import plotly.graph_objects as go

    L, D = THEMES["light"], THEMES["dark"]
    g, vl, red = art["g"], art["vl"], art["red"]
    rng = np.random.default_rng(11)

    fig = go.Figure()
    styles = []  # per trace, the colour settings for each theme

    def add(trace, light: dict, dark: dict):
        fig.add_trace(trace)
        styles.append({"light": light, "dark": dark})

    # Galaxies: one trace per declination slab, red and blue mixed in random
    # draw order so neither colour hides the other.
    slab_traces = []
    for lo, hi in SLABS:
        m = vl & (g["dec"] >= lo) & (g["dec"] < hi)
        idx = np.flatnonzero(m)
        idx = idx[rng.permutation(len(idx))]
        p = cp.cartesian(g["ra"][idx], g["dec"][idx], g["d_mpc"][idx]).astype(np.float32)
        label = ("Fig. 01's slice" if (lo, hi) == SLABS[1]
                 else f"declination {lo:+.0f}° to {hi:+.0f}°")
        slab_traces.append(len(fig.data))
        add(go.Scatter3d(x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="markers", name=label,
                         showlegend=False, hoverinfo="skip",
                         marker=dict(size=1.6, color=red[idx].astype(np.uint8), cmin=0, cmax=1,
                                     colorscale=_colorscale(L), showscale=False,
                                     opacity=OPACITY_ALL, line=dict(width=0))),
            {"marker.colorscale": [_colorscale(L)]}, {"marker.colorscale": [_colorscale(D)]})

    # Legend entries for the two classes (the points themselves are above).
    n_red, n_blue = int(red[vl].sum()), int((~red[vl]).sum())
    for key, text in (("red", f"red-sequence galaxies ({n_red:,}): old stars"),
                      ("blue", f"blue-cloud galaxies ({n_blue:,}): still forming stars")):
        add(go.Scatter3d(x=[None], y=[None], z=[None], mode="markers", name=text,
                         marker=dict(size=6, color=L[key]), hoverinfo="skip"),
            {"marker.color": L[key]}, {"marker.color": D[key]})

    x, y, z, i, j, k = _slice_mesh()
    add(go.Mesh3d(x=x, y=y, z=z, i=i, j=j, k=k, color=L["slice"], opacity=1.0,
                  name="plane of Fig. 01's slice", showlegend=True, visible="legendonly",
                  hoverinfo="skip"),
        {"color": L["slice"]}, {"color": D["slice"]})

    # Named clusters: the richest groups the friends-of-friends finder found,
    # labelled with their Abell catalogue names.
    named = [c for c in art["groups"] if c["abell"]][:10]
    cxyz = cp.cartesian([c["ra"] for c in named], [c["dec"] for c in named],
                        [c["distance_mpc"] for c in named])
    hover = [f"<b>{c['abell']}</b><br>{c['n']} bright member galaxies<br>"
             f"redshift {c['z']:.4f} · {c['distance_mpc']:,.0f} Mpc away<br>"
             f"velocity dispersion {c['sigma_v']:,.0f} km/s" for c in named]
    add(go.Scatter3d(x=cxyz[:, 0], y=cxyz[:, 1], z=cxyz[:, 2], mode="markers+text",
                     text=[c["abell"] for c in named], textposition="top center",
                     textfont=dict(color=L["label"], size=12),
                     marker=dict(size=8, color=L["label"], symbol="circle-open"),
                     name="richest clusters (hover for details)", hovertext=hover,
                     hovertemplate="%{hovertext}<extra></extra>"),
        {"textfont.color": L["label"], "marker.color": L["label"]},
        {"textfont.color": D["label"], "marker.color": D["label"]})
    add(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers+text", text=["Milky Way"],
                     textposition="bottom center", textfont=dict(color=L["label"], size=12),
                     marker=dict(size=5, color=L["label"]), name="us", showlegend=False,
                     hovertemplate="The Milky Way: where we observe from<extra></extra>"),
        {"textfont.color": L["label"], "marker.color": L["label"]},
        {"textfont.color": D["label"], "marker.color": D["label"]})

    # Menu: everything, or one slab at a time. Each button restyles only the
    # slab traces (the third argument), leaving legend, clusters and plane alone.
    buttons = []
    for name, show in [("the whole survey", slab_traces)] + [
            (fig.data[t].name, [t]) for t in slab_traces]:
        opacity = OPACITY_ALL if len(show) > 1 else OPACITY_ONE
        buttons.append(dict(label=f"show: {name}", method="restyle",
                            args=[{"visible": [t in show for t in slab_traces],
                                   "marker.opacity": [opacity] * len(slab_traces)},
                                  slab_traces]))

    axis = dict(visible=False, showbackground=False)
    fig.update_layout(
        height=820, margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor=L["surface"],
        font=dict(family="system-ui, -apple-system, 'Segoe UI', sans-serif", size=12,
                  color=L["ink"]),
        scene=dict(xaxis=axis, yaxis=axis, zaxis=axis, aspectmode="data",
                   bgcolor=L["surface"],
                   camera=dict(eye=dict(x=-0.35, y=0.3, z=1.1), up=dict(x=0, y=0, z=1))),
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(0,0,0,0)", font=dict(size=12),
                    itemsizing="constant"),
        hoverlabel=dict(bgcolor=L["surface"], font=dict(color=L["ink"], size=12)),
        updatemenus=[dict(buttons=buttons, direction="down", x=0.99, xanchor="right",
                          y=0.93, yanchor="top", bgcolor=L["menu_bg"],
                          bordercolor=L["menu_border"], font=dict(color=L["ink"], size=12),
                          showactive=True)],
    )
    theme_layout = {
        name: {"paper_bgcolor": th["surface"], "scene.bgcolor": th["surface"],
               "font.color": th["ink"], "hoverlabel.bgcolor": th["surface"],
               "hoverlabel.font.color": th["ink"], "updatemenus[0].bgcolor": th["menu_bg"],
               "updatemenus[0].bordercolor": th["menu_border"],
               "updatemenus[0].font.color": th["ink"]}
        for name, th in THEMES.items()}
    themes_js = json.dumps({name: {"layout": theme_layout[name],
                                   "traces": [s[name] for s in styles]} for name in THEMES})
    plot_html = fig.to_html(full_html=False, include_plotlyjs=True, div_id="web3d",
                            config={"displaylogo": False, "responsive": True})

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cosmic web in 3D</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']};
          --ink2: {L['ink2']}; --muted: {L['muted']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']};
            --ink2: {D['ink2']}; --muted: {D['muted']}; }}
  }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1180px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 80ch; }}
  .note {{ font-size: 0.85rem; color: var(--muted); }}
</style>
</head>
<body>
<main>
<h1>The cosmic web in 3D</h1>
<p>{int(vl.sum()):,} galaxies from the Sloan Digital Sky Survey, each at its position in space:
every galaxy brighter than absolute magnitude −20 (M<sub>r</sub> − 5 log h) between
{cp.comoving_distance(VL_Z_MIN):.0f} and {cp.comoving_distance(VL_Z_MAX):.0f} megaparsecs from the
Milky Way (300 million to 1.5 billion light-years). That is a complete census, so the density you
see is the real density of the universe, not an artefact of what the survey could detect.</p>
<p>Drag to rotate, scroll to zoom, hover the circled clusters for their details. Seen all at once
the structure is buried in projection, so use the menu at top right to show one slab of sky at a
time: walls and filaments of galaxies with empty voids between them appear, and rich clusters
are drawn out into spikes pointing back at the Milky Way — the "fingers of God", where the
orbital speeds of cluster galaxies are read as distance. Red galaxies concentrate in the clusters
and walls; blue ones are spread more evenly.</p>
{plot_html}
<p class="note">Data: SDSS DR17 Main Galaxy Sample via SkyServer; cluster names from the Abell,
Corwin &amp; Olowin (1989) catalogue. Positions are comoving distances for the Planck 2018
cosmology (H<sub>0</sub> = 67.7 km/s/Mpc). The big fan of galaxies is the survey's northern sky
region, the thin fans its three southern stripes; the empty spaces around them are where the
Milky Way's dust blocks the view, not empty space. Clusters found with a friends-of-friends group
finder ({res['groups']['n_groups_10plus']} groups with ten or more bright members). Made by run.py
in 07-cosmic-web.</p>
</main>
<script>
(function () {{
  var gd = document.getElementById("web3d");
  var THEMES = {themes_js};
  function apply(dark) {{
    var th = THEMES[dark ? "dark" : "light"];
    th.traces.forEach(function (s, i) {{ Plotly.restyle(gd, s, [i]); }});
    Plotly.relayout(gd, th.layout);
  }}
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  if (mq.matches) apply(true);
  mq.addEventListener("change", function (e) {{ apply(e.matches); }});
}})();
</script>
</body>
</html>
"""
    out = FIGURE_DIR / "cosmic_web_3d.html"
    out.write_text(page, encoding="utf-8")
    return str(out)
