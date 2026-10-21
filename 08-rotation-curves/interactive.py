"""A self-contained, interactive rotation-curve explorer.

    figures/rotation_curve_explorer.html

Every SPARC galaxy, one at a time, with the stars' mass-to-light ratio on a
slider. It lets you ask, by hand, the question the whole project turns on: can
*any* amount of stars make the visible matter match the observed rotation?
Slide the ratio up until the outer curve is matched, and watch the inner curve
overshoot.

The page embeds plotly.js and ~200 kB of data, so it works offline and can be
shared as a single file. It follows the reader's light/dark setting (both
palettes are validated for colour-blind safety), and it carries a table of the
selected galaxy's numbers for anyone who cannot use the chart.
"""

from __future__ import annotations

import json

import numpy as np
from plotly.offline import get_plotlyjs

import analysis
import plotstyle as ps
from config import (
    FIGURE_DIR,
    MAX_DISK_TOLERANCE_SIGMA,
    MAX_QUALITY_FLAG,
    MIN_INCLINATION_DEG,
    SHOWCASE,
    SPARC_PAPER,
    UPSILON_BULGE,
    UPSILON_DISK,
    UPSILON_PLAUSIBLE,
)

#: Light values are the plotstyle palette; dark values are the same hues stepped
#: for a dark surface (validated separately against it).
THEMES = {
    "light": {
        "surface": ps.SURFACE, "ink": ps.INK, "ink2": ps.INK_SECONDARY, "muted": ps.MUTED,
        "grid": ps.GRID, "axis": ps.AXIS, "observed": ps.INK, "err": ps.MUTED,
        "visible": ps.BLUE, "gas": ps.AQUA, "disk": ps.YELLOW, "bulge": ps.MAGENTA,
        "gapfill": "rgba(235,104,52,0.14)", "overfill": "rgba(42,120,214,0.16)",
    },
    "dark": {
        "surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781",
        "grid": "#2c2c2a", "axis": "#4a4a46", "observed": "#ffffff", "err": "#898781",
        "visible": "#3987e5", "gas": "#199e70", "disk": "#c98500", "bulge": "#d55181",
        "gapfill": "rgba(217,89,38,0.13)", "overfill": "rgba(57,135,229,0.20)",
    },
}


def _round(a, nd=2):
    return [round(float(v), nd) for v in np.asarray(a)]


def _record(g, showcase: set[str]) -> dict:
    c = g.curve
    md = analysis.maximum_disk(g)
    need = md.upsilon_disk_needed if md and np.isfinite(md.upsilon_disk_needed) else None
    return {
        "name": g.name, "label": g.label, "type": g.type_name,
        "d": g.distance_mpc, "inc": g.inclination_deg, "q": g.quality,
        "in_sample": g.quality <= MAX_QUALITY_FLAG and g.inclination_deg >= MIN_INCLINATION_DEG,
        "showcase": g.name in showcase,
        "r": _round(c.r), "v": _round(c.v_obs), "e": _round(c.e_v_obs),
        "gas": _round(c.v_gas), "disk": _round(c.v_disk), "bul": _round(c.v_bulge),
        # Rounded down: at exactly the maximum the visible curve touches its
        # ceiling, so rounding up would tip it over.
        "ml_max": np.floor(md.upsilon_disk_max * 1000) / 1000 if md else None,
        "ml_need": round(need, 3) if need is not None else None,
    }


def write_explorer(galaxies: dict) -> str:
    """Write the explorer page; returns its path."""
    showcase = {s.name for s in SHOWCASE}
    records = [_record(g, showcase) for g in sorted(galaxies.values(), key=lambda g: g.name)]
    order = [s.name for s in SHOWCASE]
    data_js = json.dumps({"galaxies": records, "showcase": order}, separators=(",", ":"))
    settings = json.dumps({
        "ml_fiducial": UPSILON_DISK, "bulge_ratio": UPSILON_BULGE / UPSILON_DISK,
        "plausible": list(UPSILON_PLAUSIBLE), "tol": MAX_DISK_TOLERANCE_SIGMA,
        "start": "NGC3198",
    })
    themes_js = json.dumps(THEMES)
    L, D = THEMES["light"], THEMES["dark"]

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Rotation-curve explorer</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']};
          --ink2: {L['ink2']}; --muted: {L['muted']}; --grid: {L['grid']};
          --accent: {L['visible']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']};
            --ink2: {D['ink2']}; --muted: {D['muted']}; --grid: {D['grid']};
            --accent: {D['visible']}; }}
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1100px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 82ch; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 12px 20px; align-items: end;
              margin: 16px 0 10px; }}
  .control {{ display: flex; flex-direction: column; gap: 4px; min-width: 0; }}
  .control label {{ font-size: 0.8rem; color: var(--muted); }}
  select, button {{ font: inherit; color: var(--ink); background: var(--surface);
                   border: 1px solid var(--grid); border-radius: 6px; padding: 5px 9px; }}
  select {{ max-width: 100%; }}
  button {{ cursor: pointer; }}
  button:hover {{ border-color: var(--muted); }}
  .slider {{ display: flex; align-items: center; gap: 10px; }}
  input[type=range] {{ width: min(320px, 70vw); accent-color: var(--accent); }}
  output {{ font-variant-numeric: tabular-nums; min-width: 4.5ch; font-weight: 600; }}
  .buttons {{ display: flex; flex-wrap: wrap; gap: 6px; }}
  #readout {{ min-height: 4.6em; color: var(--ink); max-width: 90ch; }}
  #plot {{ width: 100%; }}
  details {{ margin-top: 12px; }}
  summary {{ cursor: pointer; color: var(--ink2); }}
  .table-wrap {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; font-size: 0.85rem; margin-top: 8px;
          font-variant-numeric: tabular-nums; }}
  th, td {{ padding: 3px 10px; text-align: right; border-bottom: 1px solid var(--grid); }}
  th {{ color: var(--muted); font-weight: 600; }}
  .note {{ font-size: 0.85rem; color: var(--muted); margin-top: 16px; }}
</style>
</head>
<body>
<main>
<h1>Rotation-curve explorer</h1>
<p>The observed orbital speeds of gas in {len(records)} real galaxies (dots), against what
their visible stars and gas can hold in orbit (blue). The one number nobody can measure
directly is how much mass the stars carry per unit of light: their <em>mass-to-light
ratio</em>, M/L. Stellar-population models put it near {UPSILON_DISK:g} in this infrared band.
Slide it, and try to make the blue curve match the dots everywhere.</p>
<div class="controls">
  <div class="control">
    <label for="galaxy">Galaxy</label>
    <select id="galaxy"></select>
  </div>
  <div class="control">
    <label for="ml">Stellar mass-to-light ratio, M/L (disc; bulge scales with it)</label>
    <div class="slider"><input type="range" id="ml" min="0" max="1000" step="1">
      <output id="ml-value" for="ml"></output></div>
  </div>
  <div class="control">
    <label>Jump to</label>
    <div class="buttons">
      <button type="button" id="b-fid">stellar populations</button>
      <button type="button" id="b-max">maximum disc</button>
      <button type="button" id="b-need">enough for the edge</button>
    </div>
  </div>
</div>
<p id="readout" aria-live="polite"></p>
<div id="plot"></div>
<details>
  <summary>Table of the selected galaxy's numbers</summary>
  <div class="table-wrap"><table id="table">
    <thead><tr><th>radius (kpc)</th><th>observed (km/s)</th><th>±</th>
      <th>all visible</th><th>gas</th><th>stellar disc</th><th>bulge</th></tr></thead>
    <tbody></tbody></table></div>
</details>
<p class="note">Data: SPARC ({SPARC_PAPER}). Speeds add in quadrature, so the blue curve is
not the sum of the thin ones. "Maximum disc" is the heaviest the stars can be before the
visible matter alone exceeds the observed speed by more than {MAX_DISK_TOLERANCE_SIGMA:.0f}σ
somewhere; "enough for the edge" is what it takes to reach the outermost observed speed (less
{MAX_DISK_TOLERANCE_SIGMA:.0f}σ). The gas mass is measured directly and is not adjustable.
Made by run.py in 08-rotation-curves.</p>
</main>
<script>{get_plotlyjs()}</script>
<script>
(function () {{
  "use strict";
  var DATA = {data_js.replace("</", "<\\/")};
  var SET = {settings};
  var THEMES = {themes_js};
  var LOG_MIN = Math.log10(0.05), LOG_MAX = Math.log10(50);
  var byName = {{}};
  DATA.galaxies.forEach(function (g) {{ byName[g.name] = g; }});
  var sel = document.getElementById("galaxy"), slider = document.getElementById("ml");
  var out = document.getElementById("ml-value"), readout = document.getElementById("readout");
  var tbody = document.querySelector("#table tbody"), gd = document.getElementById("plot");
  var mq = window.matchMedia("(prefers-color-scheme: dark)");

  function opt(g) {{
    var o = document.createElement("option");
    o.value = g.name;
    o.textContent = g.label + " (" + g.type + ")" + (g.in_sample ? "" : " — not in sample");
    return o;
  }}
  var grp = document.createElement("optgroup");
  grp.label = "Showcase galaxies";
  DATA.showcase.forEach(function (n) {{ grp.appendChild(opt(byName[n])); }});
  sel.appendChild(grp);
  grp = document.createElement("optgroup");
  grp.label = "All " + DATA.galaxies.length + " SPARC galaxies";
  DATA.galaxies.forEach(function (g) {{ grp.appendChild(opt(g)); }});
  sel.appendChild(grp);

  function sq(v) {{ return v * Math.abs(v); }}
  function root(x) {{ return x < 0 ? -Math.sqrt(-x) : Math.sqrt(x); }}
  // The buttons set an exact value; the slider can only land on its own steps
  // (0.3% apart), which would put "maximum disc" a hair over the maximum.
  var mlExact = null;
  function mlNow() {{
    return mlExact != null ? mlExact
      : Math.pow(10, LOG_MIN + (LOG_MAX - LOG_MIN) * slider.value / 1000);
  }}
  function setMl(ml) {{
    mlExact = Math.min(50, Math.max(0.05, ml));
    slider.value = Math.round(1000 * (Math.log10(mlExact) - LOG_MIN) / (LOG_MAX - LOG_MIN));
  }}
  function fmt(x) {{ return x == null ? "—" : (x < 10 ? x.toFixed(2) : x.toFixed(1)); }}

  function model(g, ml) {{
    var mb = ml * SET.bulge_ratio;
    return {{
      disk: g.disk.map(function (v) {{ return v * Math.sqrt(ml); }}),
      bul: g.bul.map(function (v) {{ return v * Math.sqrt(mb); }}),
      vbar: g.r.map(function (_, i) {{
        return root(sq(g.gas[i]) + ml * sq(g.disk[i]) + mb * sq(g.bul[i]));
      }})
    }};
  }}

  function render() {{
    var g = byName[sel.value], ml = mlNow(), m = model(g, ml);
    var th = THEMES[mq.matches ? "dark" : "light"];
    out.textContent = fmt(ml);
    var n = g.r.length, last = n - 1;
    var top = g.v.map(function (v, i) {{ return Math.max(v, m.vbar[i]); }});
    var rr = g.r.concat(g.r.slice().reverse());
    var overshoot = [];
    g.r.forEach(function (r, i) {{
      if (m.vbar[i] > (g.v[i] + SET.tol * g.e[i]) * 1.0005) overshoot.push(r);
    }});
    var hasBulge = g.bul.some(function (v) {{ return v > 0; }});
    var hv = function (name) {{ return name + ": %{{y:.1f}} km/s<extra></extra>"; }};
    var traces = [
      {{x: rr, y: top.concat(m.vbar.slice().reverse()), fill: "toself", fillcolor: th.gapfill,
        line: {{width: 0}}, hoverinfo: "skip", name: "the gap: unseen mass"}},
      {{x: rr, y: top.concat(g.v.slice().reverse()), fill: "toself", fillcolor: th.overfill,
        line: {{width: 0}}, hoverinfo: "skip", name: "visible matter overshoots",
        showlegend: m.vbar.some(function (v, i) {{ return v > g.v[i]; }})}},
      {{x: g.r, y: g.gas, mode: "lines", line: {{color: th.gas, width: 1.6}}, name: "gas",
        hovertemplate: hv("gas")}},
      {{x: g.r, y: m.disk, mode: "lines", line: {{color: th.disk, width: 1.6}},
        name: "stellar disc", hovertemplate: hv("stellar disc")}}
    ];
    if (hasBulge) traces.push({{x: g.r, y: m.bul, mode: "lines",
      line: {{color: th.bulge, width: 1.6}}, name: "bulge", hovertemplate: hv("bulge")}});
    traces.push({{x: g.r, y: m.vbar, mode: "lines", line: {{color: th.visible, width: 3}},
      name: "all visible matter", hovertemplate: hv("all visible")}});
    traces.push({{x: g.r, y: g.v, mode: "markers", name: "observed rotation",
      error_y: {{type: "data", array: g.e, color: th.err, thickness: 1, width: 0}},
      marker: {{color: th.observed, size: 8, line: {{color: th.surface, width: 1.5}}}},
      hovertemplate: "observed: %{{y:.1f}} ± %{{error_y.array:.1f}} km/s<extra></extra>"}});

    var ceiling = Math.max.apply(null, g.v.map(function (v, i) {{ return v + g.e[i]; }}));
    var ymax = Math.max(1.22 * ceiling, 1.05 * Math.max.apply(null, m.vbar));
    function axis(title, range) {{
      return {{title: {{text: title, font: {{color: th.ink2}}}}, range: range, zeroline: false,
        gridcolor: th.grid, linecolor: th.axis, tickcolor: th.axis, ticks: "outside",
        tickfont: {{color: th.muted}}}};
    }}
    var layout = {{
      paper_bgcolor: th.surface, plot_bgcolor: th.surface, height: 520,
      font: {{family: "system-ui, -apple-system, 'Segoe UI', sans-serif", size: 13, color: th.ink}},
      margin: {{l: 62, r: 12, t: 10, b: 54}},
      xaxis: axis("Distance from the galaxy's centre (kpc)", [0, g.r[last] * 1.04]),
      yaxis: axis("Orbital speed (km/s)", [0, ymax]),
      legend: {{orientation: "h", x: 0, y: 1.0, yanchor: "bottom", font: {{color: th.ink2}},
        traceorder: "reversed"}},
      hovermode: "x unified",
      hoverlabel: {{bgcolor: th.surface, bordercolor: th.axis, font: {{color: th.ink}}}}
    }};
    Plotly.react(gd, traces, layout, {{displaylogo: false, responsive: true}});

    var pull = sq(m.vbar[last]) / (g.v[last] * g.v[last]);
    var s = "<strong>" + g.label + "</strong> (" + g.type + ", " + g.d.toFixed(1) +
      " Mpc away, tilted " + g.inc.toFixed(0) + "° from face-on). At M/L = " + fmt(ml) +
      ", the visible matter supplies <strong>" + Math.round(100 * pull) +
      "%</strong> of the pull needed to hold the outermost orbit, " + g.r[last].toFixed(1) +
      " kpc out.";
    if (overshoot.length) {{
      s += " But it now spins the galaxy faster than observed (by more than " + SET.tol +
        "σ) at " + overshoot.length + " radi" + (overshoot.length > 1 ? "i" : "us") +
        ", from " + overshoot[0].toFixed(1) + " kpc: stars this heavy are ruled out.";
    }}
    s += " Stellar populations suggest M/L ≈ " + SET.plausible[0] + "–" + SET.plausible[1] +
      "; the inner curve allows at most " + fmt(g.ml_max) + "; the edge needs " +
      fmt(g.ml_need) + ".";
    if (!g.in_sample) s += " (Left out of the project's statistics: " + (g.q > 2 ?
      "its rotation curve is flagged as disturbed." : "it is too close to face-on.") + ")";
    readout.innerHTML = s;

    var rows = [];
    for (var i = 0; i < n; i++) {{
      rows.push("<tr><td>" + g.r[i].toFixed(2) + "</td><td>" + g.v[i].toFixed(1) + "</td><td>" +
        g.e[i].toFixed(1) + "</td><td>" + m.vbar[i].toFixed(1) + "</td><td>" + g.gas[i].toFixed(1) +
        "</td><td>" + m.disk[i].toFixed(1) + "</td><td>" + (hasBulge ? m.bul[i].toFixed(1) : "—") +
        "</td></tr>");
    }}
    tbody.innerHTML = rows.join("");
  }}

  sel.addEventListener("change", render);
  slider.addEventListener("input", function () {{ mlExact = null; render(); }});
  document.getElementById("b-fid").addEventListener("click", function () {{
    setMl(SET.ml_fiducial); render(); }});
  document.getElementById("b-max").addEventListener("click", function () {{
    var g = byName[sel.value]; if (g.ml_max) {{ setMl(g.ml_max); render(); }} }});
  document.getElementById("b-need").addEventListener("click", function () {{
    var g = byName[sel.value]; if (g.ml_need) {{ setMl(g.ml_need); render(); }} }});
  mq.addEventListener("change", render);

  var start = (location.hash || "").slice(1);
  sel.value = byName[start] ? start : SET.start;
  setMl(SET.ml_fiducial);
  render();
}})();
</script>
</body>
</html>
"""
    out = FIGURE_DIR / "rotation_curve_explorer.html"
    out.write_text(page, encoding="utf-8")
    return str(out)
