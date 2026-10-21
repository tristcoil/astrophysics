"""A self-contained, interactive world map of every located meteorite.

    figures/meteorite_map.html

Every meteorite with usable coordinates, on the same equal-area Equal Earth map
as the static figures, with its name, class, mass and year on hover. Filters
narrow it by how it was found, by family, by mass and by year, which is how
the patterns in the README can be checked by hand: show only irons and the
farm belts of North America light up; show only lunar and Martian meteorites
and the map empties to the deserts and Antarctica; show only falls before 1900
and you are looking at the population map of the 19th century.

The page embeds plotly.js, the data (~1.5 MB) and Natural Earth's 1:110m
outlines, so it works offline and can be shared as a single file. The map is
projected in the browser with the same formula as geo.equal_earth. It follows
the reader's light/dark setting and carries a summary table of whatever is
currently shown, for anyone who cannot use the map.
"""

from __future__ import annotations

import json

import numpy as np
from plotly.offline import get_plotlyjs

import catalog
import geo
import plotstyle as ps
from analysis import ENV_NAMES, ENVIRONMENTS
from config import COUNTRIES_110M, FIGURE_DIR

#: Light values are the plotstyle palette; dark values are the same hues stepped
#: for a dark surface (the atlas palette's validated dark steps).
THEMES = {
    "light": {"surface": ps.SURFACE, "ink": ps.INK, "ink2": ps.INK_SECONDARY, "muted": ps.MUTED,
              "grid": ps.GRID, "axis": ps.AXIS, "land": ps.LAND, "coast": ps.COAST,
              "fall": ps.ORANGE, "find": ps.BLUE, "ring": ps.SURFACE, "hit": ps.INK},
    "dark": {"surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#383835", "land": "#2a2a28", "coast": "#4a4a46",
             "fall": "#d95926", "find": "#3987e5", "ring": "#1a1a19", "hit": "#ffffff"},
}

FAMILY_OPTIONS = [
    ("all", "All families"),
    ("ordinary", "Ordinary chondrites"),
    ("carbonaceous", "Carbonaceous chondrites"),
    ("other_chondrite", "Enstatite & other chondrites"),
    ("achondrite", "Achondrites (all)"),
    ("group:Martian", "   Martian meteorites"),
    ("group:Lunar", "   Lunar meteorites"),
    ("group:HED", "   Vesta (howardites, eucrites, diogenites)"),
    ("stony_iron", "Stony-irons"),
    ("iron", "Irons"),
    ("unclassified", "Unclassified"),
]


def _outline_rings(countries, *, step: float = 0.0) -> list:
    """Every exterior and hole ring as [lon..., lat...] pairs rounded to 0.01 degree."""
    rings = []
    for c in countries:
        for ext, holes in c.polygons:
            for ring in [ext] + holes:
                rings.append([[round(float(x), 2), round(float(y), 2)] for x, y in ring])
    return rings


def _records(df) -> dict:
    """Column-oriented, compact data for the page."""
    d = df[df["valid"] & df["located"]].copy()
    d = d.sort_values("mass_g", ascending=True, na_position="first")  # heavy ones drawn last
    classes = sorted(d["recclass"].unique())
    cls_idx = {c: i for i, c in enumerate(classes)}
    groups = sorted(d["group"].unique())
    grp_idx = {g: i for i, g in enumerate(groups)}
    countries = sorted(x for x in d["country"].dropna().unique())
    cty_idx = {c: i for i, c in enumerate(countries)}
    env_codes = list(ENVIRONMENTS) + ["unknown"]
    fams = list(catalog.FAMILIES) + ["unclassified"]
    return {
        "name": d["name"].tolist(),
        "lat": [round(float(v), 3) for v in d["lat"]],
        "lon": [round(float(v), 3) for v in d["lon"]],
        "mass": [None if not np.isfinite(v) else float(f"{v:.4g}") for v in d["mass_g"]],
        "year": [None if not np.isfinite(v) else int(v) for v in d["year"]],
        "env": [env_codes.index(e) for e in d["environment"]],
        "fam": [fams.index(f) for f in d["family"]],
        "grp": [grp_idx[g] for g in d["group"]],
        "cls": [cls_idx[c] for c in d["recclass"]],
        "cty": [cty_idx.get(c, -1) if isinstance(c, str) else -1 for c in d["country"]],
        "classes": classes, "groups": groups, "countries": countries,
        "envs": env_codes, "env_names": [ENV_NAMES.get(e, "Unknown origin") for e in env_codes],
        "fams": fams, "fam_names": [catalog.FAMILY_NAMES[f] for f in fams],
    }


def write_map(df, results) -> str:
    """Write the interactive map page; returns its path."""
    countries = geo.load_countries(COUNTRIES_110M)
    data = _records(df)
    data["rings"] = _outline_rings(countries)
    data_js = json.dumps(data, separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")
    themes_js = json.dumps(THEMES)
    L, D = THEMES["light"], THEMES["dark"]
    census = results["census"]
    dq = results["data_quality"]
    n_located = len(data["name"])
    fam_opts = "\n".join(f'<option value="{v}">{t.replace("   ", "&nbsp;&nbsp;&nbsp;")}</option>'
                         for v, t in FAMILY_OPTIONS)
    env_boxes = "\n".join(
        f'<label class="chk"><input type="checkbox" data-env="{i}" checked> {ENV_NAMES[e]}</label>'
        for i, e in enumerate(ENVIRONMENTS))
    masses = [(0, "any"), (1, "1 g"), (10, "10 g"), (100, "100 g"), (1e3, "1 kg"), (1e4, "10 kg"),
              (1e5, "100 kg"), (1e6, "1 t")]
    mass_opts = "\n".join(f'<option value="{v:g}">{t}</option>' for v, t in masses)
    y0, y1 = census["year_range"]

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Meteorite landings map</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']}; --ink2: {L['ink2']};
          --muted: {L['muted']}; --grid: {L['grid']}; --fall: {L['fall']}; --find: {L['find']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']}; --ink2: {D['ink2']};
            --muted: {D['muted']}; --grid: {D['grid']}; --fall: {D['fall']}; --find: {D['find']}; }}
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1240px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 90ch; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 12px 22px; align-items: end; margin: 14px 0 8px; }}
  .control {{ display: flex; flex-direction: column; gap: 4px; min-width: 0; }}
  .control > label, .control > span {{ font-size: 0.8rem; color: var(--muted); }}
  .checks {{ display: flex; flex-wrap: wrap; gap: 2px 14px; }}
  .chk {{ font-size: 0.92rem; color: var(--ink); cursor: pointer; white-space: nowrap; }}
  select, input[type=number], input[type=search], button {{ font: inherit; color: var(--ink);
    background: var(--surface); border: 1px solid var(--grid); border-radius: 6px; padding: 5px 9px; }}
  input[type=number] {{ width: 6.2em; }}
  input[type=search] {{ width: min(260px, 80vw); }}
  button {{ cursor: pointer; }}
  button:hover, select:hover {{ border-color: var(--muted); }}
  .years {{ display: flex; gap: 6px; align-items: center; color: var(--muted); }}
  .key {{ display: flex; gap: 16px; font-size: 0.9rem; color: var(--ink2); margin: 4px 0 0; }}
  .dot {{ display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px;
          vertical-align: -1px; }}
  #readout {{ min-height: 3em; color: var(--ink); max-width: 100ch; margin-top: 8px; }}
  #map {{ width: 100%; }}
  details {{ margin-top: 12px; }}
  summary {{ cursor: pointer; color: var(--ink2); }}
  .table-wrap {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; font-size: 0.85rem; margin-top: 8px; font-variant-numeric: tabular-nums; }}
  th, td {{ padding: 3px 10px; text-align: right; border-bottom: 1px solid var(--grid); }}
  th:first-child, td:first-child {{ text-align: left; }}
  th {{ color: var(--muted); font-weight: 600; }}
  .note {{ font-size: 0.85rem; color: var(--muted); margin-top: 16px; }}
</style>
</head>
<body>
<main>
<h1>Meteorite landings: where they fall, and where we find them</h1>
<p>Every one of the {n_located:,} meteorites in the Meteoritical Society's catalogue (2013 snapshot)
that has usable coordinates. <strong>Falls</strong> were seen coming down and picked up soon after, so
they map where people are. <strong>Finds</strong> turned up later, anywhere from a Kansas field to an
Antarctic ice sheet, so they map where rocks survive and stand out. The map is equal-area, like the
static figures: equal patches of map are equal areas of ground. Hover for details; drag to zoom,
double-click to reset.</p>
<div class="controls">
  <div class="control"><span>How it reached the catalogue</span>
    <div class="checks">{env_boxes}</div></div>
  <div class="control"><label for="family">Family</label>
    <select id="family">{fam_opts}</select></div>
  <div class="control"><label for="mass">Recovered mass at least</label>
    <select id="mass">{mass_opts}</select></div>
  <div class="control"><span>Year of fall or find</span>
    <div class="years"><input type="number" id="y0" min="{y0}" max="{y1}" value="{y0}" aria-label="from year">
      to <input type="number" id="y1" min="{y0}" max="{y1}" value="{y1}" aria-label="to year"></div></div>
  <div class="control"><label for="search">Find a meteorite by name</label>
    <input type="search" id="search" placeholder="e.g. Allende, Allan Hills 84001" list="none"></div>
  <div class="control"><span>&nbsp;</span><button type="button" id="reset">Reset filters</button></div>
</div>
<div class="key"><span><span class="dot" style="background:var(--fall)"></span>seen to fall</span>
  <span><span class="dot" style="background:var(--find)"></span>found later</span></div>
<p id="readout" aria-live="polite"></p>
<div id="map"></div>
<details>
  <summary>Table: what is shown, by route and family</summary>
  <div class="table-wrap"><table id="table"><thead></thead><tbody></tbody></table></div>
</details>
<p class="note">Data: The Meteoritical Society, Meteoritical Bulletin Database, as distributed in
2013 by the NASA Open Data Portal ("Meteorite Landings"); country outlines: Natural Earth 1:110m.
Not shown: {dq['no_location']:,} meteorites without usable coordinates, among them
{dq['northwest_africa']['n'] - dq['northwest_africa']['located']:,} "Northwest Africa" meteorites
sold through dealers with no known find site, {dq['coords_placeholder_0_0']:,} placed at a dummy
0°,&nbsp;0° and {dq['coords_zeroed_latitude']} Antarctic Yamato meteorites whose latitude was zeroed
in the file; and one found on Mars by the Opportunity rover. Antarctic meteorites often share the
coordinates of their collection area, so thousands can sit on one dot. "Hot-desert collection areas"
are areas that name their meteorites with numbers (Dhofar 019); "individual finds" have names of
their own. Mass is the total recovered. Made by run.py in 14-meteorite-landings.</p>
</main>
<script>{get_plotlyjs()}</script>
<script>
(function () {{
  "use strict";
  var D = {data_js};
  var THEMES = {themes_js};
  var N = D.name.length;

  // Equal Earth (Savric, Patterson & Jenny 2018) -- the same formula as geo.equal_earth.
  var A1 = 1.340264, A2 = -0.081106, A3 = 0.000893, A4 = 0.003796, M = Math.sqrt(3) / 2;
  function ee(lon, lat) {{
    var th = Math.asin(M * Math.sin(lat * Math.PI / 180)), t2 = th * th, t6 = t2 * t2 * t2;
    var x = 2 * Math.sqrt(3) * (lon * Math.PI / 180) * Math.cos(th) /
            (3 * (A1 + 3 * A2 * t2 + t6 * (7 * A3 + 9 * A4 * t2)));
    return [x, th * (A1 + A2 * t2 + t6 * (A3 + A4 * t2))];
  }}
  var X = new Float64Array(N), Y = new Float64Array(N);
  for (var i = 0; i < N; i++) {{ var p = ee(D.lon[i], D.lat[i]); X[i] = p[0]; Y[i] = p[1]; }}

  function path(lines) {{
    var xs = [], ys = [];
    lines.forEach(function (ring) {{
      ring.forEach(function (q) {{ var p = ee(q[0], q[1]); xs.push(p[0]); ys.push(p[1]); }});
      xs.push(null); ys.push(null);
    }});
    return [xs, ys];
  }}
  var land = path(D.rings);
  var grat = [];
  for (var lo = -150; lo <= 150; lo += 30) {{
    var r = []; for (var la = -90; la <= 90; la += 2) r.push([lo, la]); grat.push(r); }}
  for (var la2 = -60; la2 <= 60; la2 += 30) {{
    var r2 = []; for (var lo2 = -180; lo2 <= 180; lo2 += 2) r2.push([lo2, la2]); grat.push(r2); }}
  var gratP = path(grat);
  var edge = [];
  for (var a = -90; a <= 90; a += 1) edge.push([180, a]);
  for (var b = 90; b >= -90; b -= 1) edge.push([-180, b]);
  edge.push([180, -90]);
  var edgeP = path([edge]);

  var fams = D.fams, famSel = document.getElementById("family");
  var massSel = document.getElementById("mass"), y0 = document.getElementById("y0"),
      y1 = document.getElementById("y1"), search = document.getElementById("search"),
      readout = document.getElementById("readout"), gd = document.getElementById("map");
  var envBoxes = Array.prototype.slice.call(document.querySelectorAll("input[data-env]"));
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  var picked = -1;
  var byName = {{}};
  for (var k = 0; k < N; k++) byName[D.name[k].toLowerCase()] = k;

  function fmtMass(g) {{
    if (g == null) return "unknown mass";
    if (g >= 1e6) return (g / 1e6).toLocaleString(undefined, {{maximumSignificantDigits: 3}}) + " t";
    if (g >= 1e3) return (g / 1e3).toLocaleString(undefined, {{maximumSignificantDigits: 3}}) + " kg";
    return g.toLocaleString(undefined, {{maximumSignificantDigits: 3}}) + " g";
  }}
  function describe(i) {{
    var fell = D.envs[D.env[i]] === "fall";
    return "<b>" + D.name[i] + "</b><br>" + D.classes[D.cls[i]] + " · " + D.fam_names[D.fam[i]] +
      "<br>" + fmtMass(D.mass[i]) + " · " + (fell ? "fell " : "found ") + (D.year[i] == null ? "(year unknown)" : D.year[i]) +
      "<br>" + D.env_names[D.env[i]] + (D.cty[i] >= 0 ? " · " + D.countries[D.cty[i]] : "") +
      "<br>" + D.lat[i].toFixed(2) + "°, " + D.lon[i].toFixed(2) + "°";
  }}

  function keep(i, f) {{
    if (!f.envs[D.env[i]]) return false;
    if (f.fam !== "all") {{
      if (f.fam.indexOf("group:") === 0) {{ if (D.groups[D.grp[i]] !== f.fam.slice(6)) return false; }}
      else if (fams[D.fam[i]] !== f.fam) return false;
    }}
    if (f.mass > 0 && !(D.mass[i] >= f.mass)) return false;
    if (f.y0 > {y0} || f.y1 < {y1}) {{
      if (D.year[i] == null || D.year[i] < f.y0 || D.year[i] > f.y1) return false;
    }}
    return true;
  }}
  function filters() {{
    var envs = {{}};
    envBoxes.forEach(function (b) {{ envs[+b.dataset.env] = b.checked; }});
    envs[D.envs.indexOf("unknown")] = envBoxes.every(function (b) {{ return b.checked; }});
    return {{envs: envs, fam: famSel.value, mass: +massSel.value,
             y0: +y0.value || {y0}, y1: +y1.value || {y1}}};
  }}

  function render(keepView) {{
    var th = THEMES[mq.matches ? "dark" : "light"], f = filters();
    var fx = [], fy = [], ft = [], gx = [], gy = [], gt = [], nf = 0, ng = 0;
    var counts = {{}};
    for (var i = 0; i < N; i++) {{
      if (!keep(i, f)) continue;
      var key = D.env[i] + "|" + D.fam[i]; counts[key] = (counts[key] || 0) + 1;
      if (D.envs[D.env[i]] === "fall") {{ fx.push(X[i]); fy.push(Y[i]); ft.push(i); nf++; }}
      else {{ gx.push(X[i]); gy.push(Y[i]); gt.push(i); ng++; }}
    }}
    var hover = function (idx) {{ return idx.map(describe); }};
    var traces = [
      {{x: gratP[0], y: gratP[1], mode: "lines", line: {{color: th.grid, width: 0.6}}, hoverinfo: "skip",
        showlegend: false}},
      {{x: land[0], y: land[1], mode: "lines", fill: "toself", fillcolor: th.land,
        line: {{color: th.coast, width: 0.5}}, hoverinfo: "skip", showlegend: false}},
      {{x: edgeP[0], y: edgeP[1], mode: "lines", line: {{color: th.axis, width: 1}}, hoverinfo: "skip",
        showlegend: false}},
      {{type: "scattergl", x: gx, y: gy, mode: "markers", name: "found later",
        marker: {{color: th.find, size: 5, opacity: 0.75, line: {{width: 0}}}},
        text: hover(gt), hovertemplate: "%{{text}}<extra></extra>"}},
      {{type: "scattergl", x: fx, y: fy, mode: "markers", name: "seen to fall",
        marker: {{color: th.fall, size: 8, line: {{color: th.ring, width: 1.2}}}},
        text: hover(ft), hovertemplate: "%{{text}}<extra></extra>"}}
    ];
    if (picked >= 0) {{
      traces.push({{type: "scattergl", x: [X[picked]], y: [Y[picked]], mode: "markers", name: "selected",
        marker: {{color: "rgba(0,0,0,0)", size: 22, line: {{color: th.hit, width: 2}}}},
        text: [describe(picked)], hovertemplate: "%{{text}}<extra></extra>"}});
    }}
    var ax = {{visible: false, showgrid: false, zeroline: false, fixedrange: false}};
    var layout = {{
      paper_bgcolor: th.surface, plot_bgcolor: th.surface, showlegend: false,
      height: Math.max(360, Math.min(720, gd.clientWidth * 0.52)),
      margin: {{l: 4, r: 4, t: 4, b: 4}}, dragmode: "zoom", hovermode: "closest",
      font: {{family: "system-ui, -apple-system, 'Segoe UI', sans-serif", color: th.ink}},
      xaxis: Object.assign({{range: [-2.75, 2.75]}}, ax),
      yaxis: Object.assign({{range: [-1.34, 1.34], scaleanchor: "x", scaleratio: 1}}, ax),
      hoverlabel: {{bgcolor: th.surface, bordercolor: th.axis, font: {{color: th.ink, size: 13}}}}
    }};
    if (keepView && gd.layout && gd.layout.xaxis) {{
      layout.xaxis.range = gd.layout.xaxis.range.slice();
      layout.yaxis.range = gd.layout.yaxis.range.slice();
    }}
    Plotly.react(gd, traces, layout, {{displaylogo: false, responsive: true, scrollZoom: true,
      modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"]}});
    var s = "Showing <strong>" + (nf + ng).toLocaleString() + "</strong> of " + N.toLocaleString() +
      " located meteorites: " + nf.toLocaleString() + " seen to fall, " + ng.toLocaleString() + " found later.";
    if (picked >= 0) s += "<br>Selected: " + describe(picked).replace(/<br>/g, " · ");
    readout.innerHTML = s;
    table(counts);
  }}

  function table(counts) {{
    var envOrder = [0, 1, 2, 3], famOrder = [0, 1, 2, 3, 4, 5, 6];
    var head = "<tr><th>Family</th>" + envOrder.map(function (e) {{ return "<th>" + D.env_names[e] + "</th>"; }}).join("") +
      "<th>Total</th></tr>";
    var rows = famOrder.map(function (fm) {{
      var tot = 0;
      var cells = envOrder.map(function (e) {{ var c = counts[e + "|" + fm] || 0; tot += c; return "<td>" + c.toLocaleString() + "</td>"; }});
      return "<tr><td>" + D.fam_names[fm] + "</td>" + cells.join("") + "<td>" + tot.toLocaleString() + "</td></tr>";
    }});
    document.querySelector("#table thead").innerHTML = head;
    document.querySelector("#table tbody").innerHTML = rows.join("");
  }}

  function find() {{
    var q = search.value.trim().toLowerCase();
    if (!q) {{ picked = -1; render(true); return; }}
    var i = byName[q];
    if (i == null) {{
      for (var j = 0; j < N; j++) if (D.name[j].toLowerCase().indexOf(q) === 0) {{ i = j; break; }}
    }}
    if (i == null) {{ readout.innerHTML = "No located meteorite is called “" + search.value + "”."; return; }}
    picked = i;
    render(false);
    Plotly.relayout(gd, {{"xaxis.range": [X[i] - 0.25, X[i] + 0.25], "yaxis.range": [Y[i] - 0.14, Y[i] + 0.14]}});
  }}

  envBoxes.concat([famSel, massSel, y0, y1]).forEach(function (el) {{
    el.addEventListener("change", function () {{ render(true); }});
  }});
  search.addEventListener("keydown", function (e) {{ if (e.key === "Enter") find(); }});
  search.addEventListener("search", find);
  document.getElementById("reset").addEventListener("click", function () {{
    envBoxes.forEach(function (b) {{ b.checked = true; }});
    famSel.value = "all"; massSel.value = "0"; y0.value = {y0}; y1.value = {y1};
    search.value = ""; picked = -1; render(false);
  }});
  mq.addEventListener("change", function () {{ render(true); }});
  window.addEventListener("resize", function () {{ render(true); }});
  render(false);
}})();
</script>
</body>
</html>
"""
    out = FIGURE_DIR / "meteorite_map.html"
    out.write_text(page, encoding="utf-8")
    return str(out)
