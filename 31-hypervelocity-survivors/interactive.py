"""A self-contained, interactive page: the fastest stars in the Galaxy.

    figures/runaway_explorer.html

A 3D view of the Milky Way with the paths of all 17 known survivors of
thermonuclear supernovae, 3 Myr into the past and the future, and a time
slider that moves every star along its path; a sky map with Green's
supernova remnants and the traced-back explosion sites; and a table of every
star's numbers. Embeds plotly.js and its data, works offline, and follows the
reader's light/dark setting.
"""

from __future__ import annotations

import json

import numpy as np
from plotly.offline import get_plotlyjs

import plotstyle as ps
from config import FIGURE_DIR, R0_KPC, Z_SUN_KPC
from literature import CENSUS

THEMES = {
    "light": {"surface": ps.SURFACE, "ink": ps.INK, "ink2": ps.INK_SECONDARY, "muted": ps.MUTED,
              "grid": ps.GRID, "axis": ps.AXIS, "disc": "rgba(137,135,129,0.18)",
              "D6": ps.BLUE, "LP40-365": ps.ORANGE, "other": "#13895f", "sun": ps.YELLOW},
    "dark": {"surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#4a4a46", "disc": "rgba(195,194,183,0.14)",
             "D6": "#3987e5", "LP40-365": "#d95926", "other": "#1baf7a", "sun": "#eda100"},
}
CLASS_KEY = {"D6": "D6", "LP40-365": "LP40-365", "He-star donor": "other", "other": "other"}
CLASS_NAME = {"D6": "D⁶ star (surviving donor)", "LP40-365": "LP 40-365 star (Type Iax remnant?)",
              "He-star donor": "helium-star donor", "other": "oxygen white dwarf"}


def _r(a, nd=3):
    return np.round(np.asarray(a, float), nd).tolist()


def build_data(results, orbit_art, rem_art) -> dict:
    kin = results["kinematics"]
    stars = []
    for r in CENSUS:
        a = orbit_art[r.name]
        past, fut = a["past"], a["future"]          # (3, n, 151), t = 0 -> -3 / 0 -> +3 Myr
        # time axis -3 ... +3 Myr, median path over posterior draws
        path = np.concatenate([past[:, :, ::-1], fut[:, :, 1:]], axis=2)
        med = np.median(path, axis=1)
        k = kin[r.name]
        samples = [[_r(path[c, j, ::6], 2) for c in range(3)] for j in range(min(12, path.shape[1]))]
        tr = rem_art["tracks"][r.name]
        stars.append(dict(
            name=r.name, cls=CLASS_KEY[r.cls], clsname=CLASS_NAME[r.cls], disc=r.discovery,
            path=[_r(med[c], 3) for c in range(3)], samples=samples,
            sky_l=_r(((np.asarray(tr["lm"]) + 180) % 360) - 180, 2), sky_b=_r(tr["bm"], 2),
            rv=r.rv, d=[round(k["d"]["median"], 2), round(k["d"]["minus"], 2), round(k["d"]["plus"], 2)],
            vtot=round(k["v_tot"]["median"]), vej=[round(k["v_ej"]["median"]), round(k["v_ej"]["minus"]),
                                                   round(k["v_ej"]["plus"])],
            t=round(k["t_disk"]["median"], 2), punb=round(k["p_unbound"], 2),
            l=round(k["l"], 2), b=round(k["b"], 2), G=round(k["G"], 2)))
    g = rem_art["green"]
    rem = dict(name=[str(x) for x in g["name"]], l=_r(((np.asarray(g["l"]) + 180) % 360) - 180, 2),
               b=_r(g["b"], 2), r=_r(g["radius_deg"], 3))
    best = results["remnants"]["D6-2"]["best"]
    return dict(stars=stars, t=_r(np.linspace(-3, 3, 301), 2), t_sample_step=6, rem=rem,
                sun=[-R0_KPC, 0.0, Z_SUN_KPC], best=dict(name=best["name"], t=best["t_yr"]))


def build(results, orbit_art, rem_art) -> list[str]:
    data = build_data(results, orbit_art, rem_art)
    data_js = json.dumps(data, separators=(",", ":"))
    themes_js = json.dumps(THEMES)
    L, D = THEMES["light"], THEMES["dark"]
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Runaway explorer</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']}; --ink2: {L['ink2']};
          --muted: {L['muted']}; --grid: {L['grid']}; --accent: {L['D6']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']}; --ink2: {D['ink2']};
            --muted: {D['muted']}; --grid: {D['grid']}; --accent: {D['D6']}; }}
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1180px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 88ch; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 12px 24px; align-items: center; margin: 14px 0 6px; }}
  .control {{ display: flex; flex-direction: column; gap: 4px; min-width: 240px; flex: 1 1 320px; }}
  .control label, .sel label {{ font-size: 0.85rem; color: var(--muted); }}
  input[type=range] {{ width: 100%; accent-color: var(--accent); }}
  select {{ font: inherit; padding: 4px 6px; background: var(--surface); color: var(--ink);
           border: 1px solid var(--grid); border-radius: 6px; }}
  .toggle {{ display: flex; gap: 8px; align-items: center; font-size: 0.9rem; color: var(--ink2); }}
  .stats {{ font-size: 0.95rem; margin: 6px 0 10px; min-height: 1.5em; }}
  .charts {{ display: grid; grid-template-columns: 1fr; gap: 12px; }}
  .chart {{ width: 100%; height: 560px; }}
  #sky {{ height: 470px; }}
  @media (max-width: 700px) {{ .chart {{ height: 420px; }} #sky {{ height: 300px; }} }}
  table {{ border-collapse: collapse; font-size: 0.85rem; margin-top: 10px; width: 100%; }}
  th, td {{ border-bottom: 1px solid var(--grid); padding: 4px 8px; text-align: right; white-space: nowrap; }}
  th {{ white-space: normal; vertical-align: bottom; }}
  th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) {{ text-align: left; }}
  th {{ color: var(--muted); font-weight: 600; }}
  tr.current td {{ font-weight: 700; }}
  .table-wrap {{ overflow-x: auto; }}
  details {{ margin-top: 14px; }}
  summary {{ cursor: pointer; color: var(--ink2); }}
</style>
</head>
<body>
<main>
<h1>The fastest stars in the Galaxy</h1>
<p>When one white dwarf in a very tight pair explodes as a Type Ia supernova, its companion is released at
the orbital speed it had at that instant — 1,000 to 2,500 km/s — and leaves the Milky Way for good. These
are the 17 stars proposed so far as survivors of such thermonuclear explosions. Their paths are drawn
3 million years into the past (to where they were launched from the disc) and into the future. Move the
time slider to watch them fly; pick a star to follow it. Positions come from Gaia DR3 parallaxes and proper
motions plus published radial velocities; the faint lines show how uncertain the distances are.</p>
<div class="controls">
  <div class="control">
    <label for="time">time: <span id="tlabel"></span></label>
    <input type="range" id="time" min="0" max="300" step="1" value="150">
  </div>
  <div class="sel"><label for="star">follow a star</label><br>
    <select id="star"><option value="-1">all 17 runaways</option></select></div>
  <label class="toggle"><input type="checkbox" id="spread" checked> show distance uncertainty</label>
</div>
<div class="stats" id="stats" aria-live="polite"></div>
<div class="charts">
  <div id="space" class="chart" role="img" aria-label="3D view of the Milky Way with the runaways' paths"></div>
  <div id="sky" class="chart" role="img" aria-label="Sky map in Galactic coordinates with supernova remnants"></div>
</div>
<p style="font-size:0.85rem">Top: Galactocentric coordinates in kpc (the Sun is the yellow diamond at x = −8.1 kpc;
the grey disc is a 15-kpc guide, not data; drag to rotate, scroll to zoom). Bottom: today's sky in Galactic
coordinates (longitude increasing to the left, the Galactic Centre in the middle); open circles are the
310 supernova remnants of Green (2025), sized by their radius; each runaway's short line marks where an
explosion up to 300,000 years ago would appear today. Only D6-2's line ends inside a remnant
({data['best']['name']}, about {data['best']['t'] / 1e3:.0f},000 years ago).</p>
<details open>
<summary>Every star's numbers (table)</summary>
<div class="table-wrap"><table id="tbl"><thead><tr><th>star</th><th>class</th><th>G (mag)</th>
<th>RV (km/s)</th><th>distance (kpc)</th><th>speed (km/s)</th>
<th>launch speed (km/s)</th><th>flight time (Myr)</th><th>P(unbound)</th></tr></thead>
<tbody></tbody></table></div>
</details>
<p style="margin-top:14px;font-size:0.85rem">Distances: posterior medians (prior flat in absolute magnitude, launch
speed below 3,000 km/s; El-Badry et al. 2023's method). Launch speed: today's velocity minus the disc's rotation
where the path crosses the midplane. Flight time: negative = still heading toward the disc. Orbits in McMillan
(2017)'s potential with galpy. Project L5 of the Astrophysics Project Atlas (AI-generated; see README).</p>
</main>
<script>{get_plotlyjs()}</script>
<script>
(function () {{
  var D = {data_js};
  var T = {themes_js};
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  var slider = document.getElementById("time"), sel = document.getElementById("star"),
      spread = document.getElementById("spread");
  function theme() {{ return mq.matches ? T.dark : T.light; }}
  D.stars.forEach(function (s, i) {{
    var o = document.createElement("option"); o.value = i; o.textContent = s.name; sel.appendChild(o);
  }});
  function discMesh(t) {{
    var x = [], y = [], z = [], n = 64;
    x.push(0); y.push(0); z.push(0);
    for (var k = 0; k < n; k++) {{ var a = 2 * Math.PI * k / n; x.push(15 * Math.cos(a)); y.push(15 * Math.sin(a)); z.push(0); }}
    var i = [], j = [], kk = [];
    for (var k = 1; k <= n; k++) {{ i.push(0); j.push(k); kk.push(k % n + 1); }}
    return {{ type: "mesh3d", x: x, y: y, z: z, i: i, j: j, k: kk, color: t.disc, opacity: 0.35,
             hoverinfo: "skip", showlegend: false, flatshading: true }};
  }}
  function render() {{
    var t = theme(), it = +slider.value, tnow = D.t[it], pick = +sel.value;
    document.getElementById("tlabel").textContent =
      (tnow === 0 ? "today" : (Math.abs(tnow).toFixed(2) + " Myr " + (tnow < 0 ? "ago" : "from now")));
    var traces = [discMesh(t)];
    D.stars.forEach(function (s, i) {{
      var on = pick < 0 || pick === i, col = t[s.cls];
      if (spread.checked && on) {{
        s.samples.forEach(function (p) {{
          traces.push({{ type: "scatter3d", mode: "lines", x: p[0], y: p[1], z: p[2], hoverinfo: "skip",
                        showlegend: false, line: {{ color: col, width: 1 }}, opacity: 0.25 }});
        }});
      }}
      traces.push({{ type: "scatter3d", mode: "lines", x: s.path[0], y: s.path[1], z: s.path[2], name: s.name,
                    showlegend: false, line: {{ color: on ? col : t.muted, width: on ? 4 : 1.5 }},
                    opacity: on ? 1 : 0.35, hoverinfo: "skip" }});
      traces.push({{ type: "scatter3d", mode: "markers" + (pick === i ? "+text" : ""), x: [s.path[0][it]], y: [s.path[1][it]],
                    z: [s.path[2][it]], text: [s.name], textposition: "top center",
                    textfont: {{ color: t.ink2, size: 11 }}, showlegend: false,
                    marker: {{ size: on ? 5 : 3, color: on ? col : t.muted, line: {{ color: t.surface, width: 1 }} }},
                    hovertemplate: s.name + "<br>x %{{x:.1f}}, y %{{y:.1f}}, z %{{z:.1f}} kpc<extra></extra>" }});
    }});
    traces.push({{ type: "scatter3d", mode: "markers+text", x: [D.sun[0]], y: [D.sun[1]], z: [D.sun[2]], text: ["Sun"],
                  textposition: "bottom center", textfont: {{ color: t.ink2, size: 11 }}, showlegend: false,
                  marker: {{ size: 6, color: t.sun, symbol: "diamond" }}, hoverinfo: "skip" }});
    traces.push({{ type: "scatter3d", mode: "markers+text", x: [0], y: [0], z: [0], text: ["Galactic Centre"],
                  textposition: "bottom center", textfont: {{ color: t.ink2, size: 11 }}, showlegend: false,
                  marker: {{ size: 4, color: t.ink, symbol: "cross" }}, hoverinfo: "skip" }});
    var ax = function (title, c) {{ return {{ title: {{ text: title, font: {{ color: t.ink2 }} }}, range: [c - 13, c + 13],
             gridcolor: t.grid, zerolinecolor: t.axis, backgroundcolor: t.surface, showbackground: false,
             tickfont: {{ color: t.muted }} }}; }};
    Plotly.react("space", traces, {{
      paper_bgcolor: t.surface, font: {{ color: t.ink }}, margin: {{ l: 0, r: 0, t: 0, b: 0 }},
      scene: {{ xaxis: ax("x (kpc)", D.sun[0] + 3), yaxis: ax("y (kpc)", 0), zaxis: ax("z (kpc)", 0), aspectmode: "cube",
               camera: {{ eye: {{ x: -1.0, y: -1.0, z: 0.5 }} }} }},
      uirevision: "keep"
    }}, {{ displayModeBar: false, responsive: true }});

    // Galactic longitude increases to the LEFT on sky maps: plot -l, hover the real l
    var neg = function (a) {{ return a.map(function (x) {{ return -x; }}); }};
    var sk = [{{ type: "scattergeo", mode: "markers", lon: neg(D.rem.l), lat: D.rem.b, text: D.rem.name,
                marker: {{ size: D.rem.r.map(function (r) {{ return Math.max(3, Math.min(18, 3 + 6 * r)); }}),
                          color: "rgba(0,0,0,0)", line: {{ color: t.muted, width: 1 }} }},
                hovertemplate: "%{{text}}<extra>remnant</extra>", showlegend: false }}];
    D.stars.forEach(function (s, i) {{
      var on = pick < 0 || pick === i, col = t[s.cls];
      sk.push({{ type: "scattergeo", mode: "lines", lon: neg(s.sky_l), lat: s.sky_b, hoverinfo: "skip", showlegend: false,
                line: {{ color: on ? col : t.muted, width: on ? 2.5 : 1 }} }});
      sk.push({{ type: "scattergeo", mode: "markers", lon: [-s.sky_l[0]], lat: [s.sky_b[0]], showlegend: false,
                marker: {{ size: on ? 9 : 5, color: on ? col : t.muted, line: {{ color: t.surface, width: 1 }} }},
                customdata: [s.l], hovertemplate: s.name + "<br>l %{{customdata:.1f}}°, b %{{lat:.1f}}°<extra></extra>" }});
    }});
    Plotly.react("sky", sk, {{
      paper_bgcolor: t.surface, font: {{ color: t.ink }}, margin: {{ l: 0, r: 0, t: 0, b: 0 }},
      geo: {{ projection: {{ type: "mollweide", scale: 1 }}, showland: false, showcoastlines: false, showframe: true,
             showcountries: false, showlakes: false, showocean: false, framecolor: t.axis, bgcolor: t.surface,
             lonaxis: {{ range: [-180, 180], showgrid: true, gridcolor: t.grid, dtick: 60 }},
             lataxis: {{ range: [-90, 90], showgrid: true, gridcolor: t.grid, dtick: 30 }} }},
      uirevision: "keep"
    }}, {{ displayModeBar: false, responsive: true }});

    var st = document.getElementById("stats");
    if (pick >= 0) {{
      var s = D.stars[pick];
      st.innerHTML = "<b>" + s.name + "</b> \\u00b7 " + s.clsname + " \\u00b7 launched at <b>" + s.vej[0] +
        "</b> km/s (68%: " + (s.vej[0] - s.vej[1]) + "\\u2013" + (s.vej[0] + s.vej[2]) + ") \\u00b7 " + s.d[0] +
        " kpc away \\u00b7 found by " + s.disc;
    }} else {{
      st.innerHTML = "17 runaways: <b>8</b> D\\u2076 stars, <b>7</b> LP 40-365 stars, US 708 and J1240+6710";
    }}
    document.querySelectorAll("#tbl tbody tr").forEach(function (r, j) {{ r.className = j === pick ? "current" : ""; }});
  }}
  var tb = document.querySelector("#tbl tbody");
  D.stars.forEach(function (s) {{
    var tr = document.createElement("tr");
    var cells = [s.name, s.clsname, s.G.toFixed(2), s.rv, s.d[0].toFixed(2) + " (\\u2212" + s.d[1].toFixed(2) + " +" + s.d[2].toFixed(2) + ")",
                 s.vtot, s.vej[0] + " (\\u2212" + s.vej[1] + " +" + s.vej[2] + ")", s.t.toFixed(2), s.punb.toFixed(2)];
    cells.forEach(function (c) {{ var td = document.createElement("td"); td.textContent = c; tr.appendChild(td); }});
    tb.appendChild(tr);
  }});
  slider.addEventListener("input", render);
  sel.addEventListener("change", render);
  spread.addEventListener("change", render);
  mq.addEventListener("change", render);
  render();
}})();
</script>
</body>
</html>
"""
    out = FIGURE_DIR / "runaway_explorer.html"
    out.write_text(page, encoding="utf-8")
    return [str(out)]
