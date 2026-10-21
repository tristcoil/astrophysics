"""A self-contained, interactive page: the sausage without a velocity cut.

    figures/halo_explorer.html

Slide through iron abundance and watch the Milky Way's halo change shape:
below [M/H] ~ -1.7 the halo is a round cloud, around -1.3 it is a long thin
"sausage" (the debris of Gaia-Sausage-Enceladus), and towards -0.8 the
rotating thick disc takes over. No star was picked by its velocity -- only
by the iron abundance Gaia measured in its spectrum -- which is the
non-circular version of the discovery. A toggle shows which stars the usual
velocity cut would have kept.

The page embeds plotly.js and ~0.5 MB of data, works offline, follows the
reader's light/dark setting, and lists every window's numbers in a table.
"""

from __future__ import annotations

import json

import numpy as np
from plotly.offline import get_plotlyjs

import galaxy as gx
import halo as hl
import plotstyle as ps
from config import FIGURE_DIR, HALO_SPEED_KMS, REFERENCE
from fetch_data import GSPSPEC_PATH
from gd1 import _read

#: Metallicity windows: the same as the analysis (0.3 dex wide, every 0.1 dex).
WINDOW_HALF = hl.FINE_WINDOW_HALF
CENTRES = hl.FINE_CENTRES

THEMES = {
    "light": {"surface": ps.SURFACE, "ink": ps.INK, "ink2": ps.INK_SECONDARY, "muted": ps.MUTED,
              "grid": ps.GRID, "axis": ps.AXIS, "accent": ps.BLUE, "cut": ps.ORANGE,
              "ref": "rgba(235,104,52,0.16)"},
    "dark": {"surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#4a4a46", "accent": "#3987e5", "cut": "#d95926",
             "ref": "rgba(217,89,38,0.18)"},
}


def build_data() -> dict:
    mp = hl.enrich(_read(GSPSPEC_PATH), gx.AnalyticHelmi17(), gx.mcmillan17())
    s = hl.trusted_giants(mp)
    mh, vr, vt, vp = mp["mh"][s], mp["v_r"][s], mp["v_theta"][s], mp["v_phi"][s]
    passing = mp["speed_lsr"][s] > HALO_SPEED_KMS
    keep = mh < CENTRES[-1] + WINDOW_HALF
    windows = hl.fine_windows(mh, vr, vt, vp, passing, CENTRES, WINDOW_HALF)
    return {
        "mh": [round(float(x), 2) for x in mh[keep]],
        "vr": [int(round(float(x))) for x in vr[keep]],
        "vp": [int(round(float(x))) for x in vp[keep]],
        "cut": [int(x) for x in passing[keep]],
        "windows": windows,
        "half": WINDOW_HALF,
        "belokurov": REFERENCE["belokurov2018"],
        "n_total": int(keep.sum()),
    }


def write_explorer(data: dict | None = None) -> str:
    data = data or build_data()
    data_js = json.dumps(data, separators=(",", ":"))
    themes_js = json.dumps(THEMES)
    L, D = THEMES["light"], THEMES["dark"]
    start = int(np.argmin(np.abs(CENTRES + 1.25)))
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Halo explorer</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']}; --ink2: {L['ink2']};
          --muted: {L['muted']}; --grid: {L['grid']}; --accent: {L['accent']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']}; --ink2: {D['ink2']};
            --muted: {D['muted']}; --grid: {D['grid']}; --accent: {D['accent']}; }}
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1100px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 84ch; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 12px 24px; align-items: center; margin: 16px 0 6px; }}
  .control {{ display: flex; flex-direction: column; gap: 4px; min-width: 260px; flex: 1 1 320px; }}
  .control label {{ font-size: 0.85rem; color: var(--muted); }}
  input[type=range] {{ width: 100%; accent-color: var(--accent); }}
  .toggle {{ display: flex; gap: 8px; align-items: center; font-size: 0.9rem; color: var(--ink2); }}
  .stats {{ font-size: 0.95rem; margin: 6px 0 10px; }}
  .stats b {{ color: var(--ink); }}
  .charts {{ display: grid; grid-template-columns: minmax(0, 3fr) minmax(0, 2fr); gap: 16px; }}
  @media (max-width: 760px) {{ .charts {{ grid-template-columns: 1fr; }} }}
  .chart {{ width: 100%; height: 470px; }}
  table {{ border-collapse: collapse; font-size: 0.85rem; margin-top: 18px; width: 100%; }}
  th, td {{ border-bottom: 1px solid var(--grid); padding: 4px 8px; text-align: right; }}
  th {{ color: var(--muted); font-weight: 600; }}
  tr.current td {{ font-weight: 700; }}
  details {{ margin-top: 14px; }}
  summary {{ cursor: pointer; color: var(--ink2); }}
  .table-wrap {{ overflow-x: auto; }}
</style>
</head>
<body>
<main>
<h1>The sausage, without a velocity cut</h1>
<p>Each dot is a giant star within 3 kpc of the Sun, with its velocity towards (+) or away from (−)
the Galactic Centre across and its speed along the Galaxy's rotation up. Stars are chosen <b>only</b>
by the iron abundance Gaia measured in their spectra ([M/H]; −1 means a tenth of the Sun's), never by
how they move. Slide through iron abundance: the halo is a round cloud at the metal-poor end, becomes a
long thin "sausage" — the debris of the Gaia-Sausage-Enceladus galaxy — around [M/H] ≈ −1.3, and gives
way to the rotating thick disc (the blob at the top) towards −0.8.</p>
<div class="controls">
  <div class="control">
    <label for="mh">iron abundance window: <span id="mhlabel"></span></label>
    <input type="range" id="mh" min="0" max="{len(CENTRES) - 1}" step="1" value="{start}">
  </div>
  <label class="toggle"><input type="checkbox" id="cut"> highlight the stars a
    {HALO_SPEED_KMS:.0f} km/s velocity cut would keep</label>
</div>
<div class="stats" id="stats" aria-live="polite"></div>
<div class="charts">
  <div id="scatter" class="chart" role="img" aria-label="Velocity scatter plot of the selected stars"></div>
  <div id="beta" class="chart" role="img" aria-label="Halo anisotropy against iron abundance"></div>
</div>
<details>
<summary>Every window's numbers (table)</summary>
<div class="table-wrap"><table id="tbl"><thead><tr><th>[M/H] window</th><th>giants</th>
<th>halo β</th><th>halo share</th><th>σ_r, σ_θ, σ_φ (km/s)</th><th>velocity cut keeps</th></tr></thead>
<tbody></tbody></table></div>
</details>
<p style="margin-top:14px;font-size:0.85rem">β = 1 − (σ_θ² + σ_φ²)/(2σ_r²) of the halo component of a
disc + halo Gaussian mixture: 0 for random orbits, 1 for purely radial ones. Belokurov et al. (2018)
found β ≈ 0.9 above [Fe/H] ≈ −1.7 and 0.2–0.4 below (shaded). Data: Gaia DR3 radial velocities and
GSP-Spec metallicities (Recio-Blanco et al. 2023); project H3 of the Astrophysics Project Atlas.</p>
</main>
<script>{get_plotlyjs()}</script>
<script>
(function () {{
  var D = {data_js};
  var T = {themes_js};
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  var slider = document.getElementById("mh"), cut = document.getElementById("cut");
  function theme() {{ return mq.matches ? T.dark : T.light; }}
  function fmt(x) {{ return (x < 0 ? "\\u2212" : "") + Math.abs(x).toFixed(2); }}
  function axisStyle(t, title, range) {{
    return {{ title: {{ text: title, font: {{ color: t.ink2, size: 12 }} }}, range: range,
             gridcolor: t.grid, zerolinecolor: t.axis, linecolor: t.axis, tickfont: {{ color: t.muted }} }};
  }}
  function render() {{
    var t = theme(), i = +slider.value, w = D.windows[i];
    var lo = w.c - D.half, hi = w.c + D.half;
    document.getElementById("mhlabel").textContent = fmt(lo) + " to " + fmt(hi);
    var xs = [], ys = [], xc = [], yc = [];
    for (var k = 0; k < D.mh.length; k++) {{
      if (D.mh[k] >= lo && D.mh[k] < hi) {{
        if (cut.checked && D.cut[k]) {{ xc.push(D.vr[k]); yc.push(D.vp[k]); }}
        else {{ xs.push(D.vr[k]); ys.push(D.vp[k]); }}
      }}
    }}
    var base = cut.checked ? t.muted : t.accent;
    var traces = [{{ x: xs, y: ys, type: "scattergl", mode: "markers", name: cut.checked ? "removed by the cut" : "giants",
                    marker: {{ size: 4, color: base, opacity: 0.55 }}, hovertemplate: "v_r %{{x}} km/s<br>v_\\u03c6 %{{y}} km/s<extra></extra>" }}];
    if (cut.checked) traces.push({{ x: xc, y: yc, type: "scattergl", mode: "markers", name: "kept by the cut",
                                   marker: {{ size: 4, color: t.cut, opacity: 0.8 }},
                                   hovertemplate: "v_r %{{x}} km/s<br>v_\\u03c6 %{{y}} km/s<extra></extra>" }});
    Plotly.react("scatter", traces, {{
      paper_bgcolor: t.surface, plot_bgcolor: t.surface, font: {{ color: t.ink }},
      margin: {{ l: 60, r: 10, t: 10, b: 50 }}, showlegend: cut.checked,
      legend: {{ x: 0.01, y: 0.01, bgcolor: "rgba(0,0,0,0)", font: {{ color: t.ink2, size: 11 }} }},
      xaxis: axisStyle(t, "v_r: towards (+) / away from (\\u2212) the Galactic Centre (km/s)", [-450, 450]),
      yaxis: axisStyle(t, "v_\\u03c6: rotation (km/s)", [-350, 400])
    }}, {{ displayModeBar: false, responsive: true }});
    var bx = [], by = [];
    D.windows.forEach(function (r) {{ if (r.beta !== undefined) {{ bx.push(r.c); by.push(r.beta); }} }});
    var B = D.belokurov;
    Plotly.react("beta", [
      {{ x: [B.feh_transition, -0.7, -0.7, B.feh_transition], y: [0.85, 0.85, 0.95, 0.95], fill: "toself", mode: "lines",
         fillcolor: t.ref, line: {{ width: 0 }}, hoverinfo: "skip", showlegend: false, type: "scatter" }},
      {{ x: [-2.6, B.feh_transition, B.feh_transition, -2.6], y: [B.beta_metal_poor_halo[0], B.beta_metal_poor_halo[0],
         B.beta_metal_poor_halo[1], B.beta_metal_poor_halo[1]], fill: "toself", mode: "lines", fillcolor: t.ref, line: {{ width: 0 }},
         hoverinfo: "skip", showlegend: false, type: "scatter" }},
      {{ x: bx, y: by, type: "scatter", mode: "lines+markers", line: {{ color: t.accent, width: 2 }},
         marker: {{ size: 6, color: t.accent }}, hovertemplate: "[M/H] %{{x:.2f}}: \\u03b2 = %{{y:.2f}}<extra></extra>",
         showlegend: false }},
      {{ x: [w.c], y: [w.beta === undefined ? null : w.beta], type: "scatter", mode: "markers",
         marker: {{ size: 14, color: "rgba(0,0,0,0)", line: {{ color: t.ink, width: 2 }} }}, hoverinfo: "skip",
         showlegend: false }}
    ], {{
      paper_bgcolor: t.surface, plot_bgcolor: t.surface, font: {{ color: t.ink }},
      margin: {{ l: 55, r: 10, t: 10, b: 50 }},
      xaxis: axisStyle(t, "iron abundance [M/H]", [-2.65, -0.75]),
      yaxis: axisStyle(t, "halo anisotropy \\u03b2", [0, 1.02]),
      annotations: [{{ x: -2.55, y: 0.99, xanchor: "left", text: "shaded: Belokurov et al. 2018", showarrow: false,
                       font: {{ size: 11, color: t.ink2 }} }}]
    }}, {{ displayModeBar: false, responsive: true }});
    var s = "<b>" + w.n.toLocaleString() + "</b> giants";
    if (w.beta !== undefined) {{
      s += " \\u00b7 halo \\u03b2 = <b>" + w.beta.toFixed(2) + "</b> \\u00b7 halo share " + Math.round(100 * w.share) +
           "% \\u00b7 a velocity cut would keep " + Math.round(100 * w.comp) + "% of the halo";
    }} else {{ s += " \\u00b7 too few for a fit"; }}
    document.getElementById("stats").innerHTML = s;
    var rows = document.querySelectorAll("#tbl tbody tr");
    rows.forEach(function (r, j) {{ r.className = j === i ? "current" : ""; }});
  }}
  var tb = document.querySelector("#tbl tbody");
  D.windows.forEach(function (r) {{
    var tr = document.createElement("tr");
    var cells = [fmt(r.c - D.half) + " to " + fmt(r.c + D.half), r.n.toLocaleString(),
                 r.beta === undefined ? "\\u2013" : r.beta.toFixed(2),
                 r.share === undefined ? "\\u2013" : Math.round(100 * r.share) + "%",
                 r.sig === undefined ? "\\u2013" : r.sig.join(", "),
                 r.comp === undefined ? "\\u2013" : Math.round(100 * r.comp) + "%"];
    cells.forEach(function (c) {{ var td = document.createElement("td"); td.textContent = c; tr.appendChild(td); }});
    tb.appendChild(tr);
  }});
  slider.addEventListener("input", render);
  cut.addEventListener("change", render);
  mq.addEventListener("change", render);
  render();
}})();
</script>
</body>
</html>
"""
    out = FIGURE_DIR / "halo_explorer.html"
    out.write_text(page, encoding="utf-8")
    return str(out)


if __name__ == "__main__":
    import warnings

    warnings.filterwarnings("ignore")
    print(write_explorer())
