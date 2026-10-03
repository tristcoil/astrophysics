"""A self-contained, interactive page: the whole sunspot record, zoomable.

    figures/solar_cycle_explorer.html

Top: every daily sunspot number since 1818, the monthly means since 1749 and the
13-month smoothed curve. Bottom, on the same time axis: every sunspot group since
1874 as one dot at its mean position, sized by its largest area -- the butterfly
diagram, group by group. Zoom either panel and both follow, so you can go from 275
years down to a single month and see which spots made a given day's number.

The page embeds plotly.js and ~2 MB of data, so it works offline and can be shared
as a single file. It follows the reader's light/dark setting (both palettes pass
the colour-blind checks against their own surface), and carries a table of all 25
cycles for anyone who cannot use the chart.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from plotly.offline import get_plotlyjs

import plotstyle as ps
from config import FIGURE_DIR, SILSO_DOI, USAF_AREA_FACTOR_HATHAWAY

#: Light values are the plotstyle palette; dark values are the same hues stepped
#: for a dark surface (validated separately against it).
THEMES = {
    "light": {"surface": ps.SURFACE, "ink": ps.INK, "ink2": ps.INK_SECONDARY, "muted": ps.MUTED,
              "grid": ps.GRID, "axis": ps.AXIS, "smoothed": ps.BLUE, "groups": ps.SUN_RAMP[4],
              "daily": "rgba(137,135,129,0.35)", "monthly": "rgba(137,135,129,0.9)"},
    "dark": {"surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#383835", "smoothed": "#3987e5", "groups": "#d95926",
             "daily": "rgba(137,135,129,0.35)", "monthly": "rgba(195,194,183,0.75)"},
}

DAILY_START = "1818-01-01"
GROUP_EPOCH = "1874-01-01"


def _num(v, nd=1):
    return None if not np.isfinite(v) else round(float(v), nd)


def _payload(tables, results) -> dict:
    daily = tables["daily"]
    d = daily[daily["date"] >= DAILY_START]
    days = pd.date_range(DAILY_START, d["date"].max(), freq="D")
    sn = d.set_index("date")["sn"].reindex(days)
    m = tables["monthly"]

    s = tables["spots"]
    epoch = pd.Timestamp(GROUP_EPOCH)
    g = s.groupby("track").agg(first=("date", "min"), last=("date", "max"), lat=("lat", "mean"),
                               area=("area_adj", "max"), cycle=("cycle", "first"), era=("era", "first"),
                               group=("group", "first"))
    g = g.sort_values("first")
    t_mid = (g["first"] + (g["last"] - g["first"]) / 2 - epoch).dt.days.to_numpy()
    groups = {
        "d": t_mid.tolist(),
        "n": ((g["last"] - g["first"]).dt.days + 1).to_numpy().tolist(),
        "lat": np.round(g["lat"].to_numpy(), 1).tolist(),
        "a": np.round(g["area"].to_numpy()).astype(int).tolist(),
        "c": g["cycle"].astype(int).tolist(),
        "id": [(("RGO " if e == "RGO" else "NOAA ") + str(n)) for e, n in zip(g["era"], g["group"])],
    }
    cyc = tables["cycles"]
    cycles = [{"n": int(r.cycle), "min": r.min_ym, "min_sn": round(float(r.min_sn), 1), "max": r.max_ym,
               "max_sn": round(float(r.max_sn), 1), "rise": round(float(r.rise_yr), 1),
               "len": None if not np.isfinite(r.length_yr) else round(float(r.length_yr), 1),
               "t_max": round(float(r.t_max), 3)} for r in cyc.itertuples()]
    return {
        "daily_start": DAILY_START, "daily": [None if not np.isfinite(v) else int(v) for v in sn.to_numpy()],
        "monthly": [_num(v) for v in m["sn"].to_numpy()],
        "smoothed": [_num(v) for v in m["smoothed"].to_numpy()],
        "group_epoch": GROUP_EPOCH, "groups": groups, "cycles": cycles,
        "last_day": str(d["date"].max().date()),
    }


def _row(c: dict) -> str:
    length = "" if c["len"] is None else f"{c['len']:.1f}"
    return (f"<tr><td>{c['n']}</td><td>{c['min']}</td><td>{c['min_sn']:.1f}</td><td>{c['max']}</td>"
            f"<td>{c['max_sn']:.1f}</td><td>{c['rise']:.1f}</td><td>{length}</td></tr>")


def write_page(tables, results) -> str:
    data = _payload(tables, results)
    data_js = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    themes_js = json.dumps(THEMES)
    L, D = THEMES["light"], THEMES["dark"]
    n_groups = len(data["groups"]["d"])
    c25 = results["cycle25"]
    rows = "\n".join(_row(c) for c in data["cycles"])

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Sunspot cycle explorer</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']}; --ink2: {L['ink2']};
          --muted: {L['muted']}; --grid: {L['grid']}; --accent: {L['smoothed']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']}; --ink2: {D['ink2']};
            --muted: {D['muted']}; --grid: {D['grid']}; --accent: {D['smoothed']}; }}
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1200px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 88ch; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 8px 10px; align-items: center; margin: 14px 0 6px; }}
  .controls span {{ font-size: 0.8rem; color: var(--muted); margin-right: 4px; }}
  button {{ font: inherit; font-size: 0.9rem; color: var(--ink); background: var(--surface);
           border: 1px solid var(--grid); border-radius: 6px; padding: 4px 10px; cursor: pointer; }}
  button:hover, button:focus-visible {{ border-color: var(--muted); }}
  label.toggle {{ font-size: 0.9rem; color: var(--ink2); display: inline-flex; gap: 6px; align-items: center; }}
  #plot {{ width: 100%; height: min(78vh, 860px); min-height: 560px; }}
  details {{ margin-top: 14px; }}
  summary {{ cursor: pointer; color: var(--ink2); }}
  .table-wrap {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; font-size: 0.85rem; margin-top: 8px; font-variant-numeric: tabular-nums; }}
  th, td {{ padding: 3px 12px; text-align: right; border-bottom: 1px solid var(--grid); }}
  th {{ color: var(--muted); font-weight: 600; }}
  .note {{ font-size: 0.85rem; color: var(--muted); margin-top: 16px; }}
</style>
</head>
<body>
<main>
<h1>275 years of sunspots, zoomable</h1>
<p>Top: the international sunspot number &mdash; every day since 1818 (grey dots), monthly since 1749,
and the 13-month smoothed curve that defines the {len(data['cycles'])} numbered cycles. Bottom, on the
same time axis: all {n_groups:,} sunspot groups seen since 1874, one dot each at its mean position on the
Sun, sized by its largest area. Each cycle's spots start at high latitude and drift toward the equator:
the butterfly wings. Drag to zoom (both panels follow), double-click to reset; hover for details.</p>
<div class="controls" role="group" aria-label="Zoom presets">
  <span>Jump to</span>
  <button type="button" data-range="1749-01-01,{data['last_day']}">everything</button>
  <button type="button" data-range="1874-01-01,{data['last_day']}">since 1874</button>
  <button type="button" data-range="1790-01-01,1830-12-31">Dalton minimum</button>
  <button type="button" data-range="1954-01-01,1965-01-01">Cycle 19, the strongest</button>
  <button type="button" data-range="2007-01-01,{data['last_day']}">Cycles 24 and 25</button>
  <button type="button" data-range="2008-06-01,2010-06-01">the 2009 minimum</button>
  <label class="toggle"><input type="checkbox" id="show-daily" checked> daily values</label>
</div>
<div id="plot"></div>
<details>
  <summary>Table: all {len(data['cycles'])} cycles (13-month smoothed sunspot number)</summary>
  <div class="table-wrap"><table>
    <thead><tr><th>cycle</th><th>minimum</th><th>SN at minimum</th><th>maximum</th><th>SN at maximum</th>
      <th>rise (yr)</th><th>length (yr)</th></tr></thead>
    <tbody>
{rows}
    </tbody></table></div>
</details>
<p class="note">Cycle 25 peaked at {c25['observed_max']:.1f} in {c25['observed_max_date']}; the official 2019 forecast
was {c25['panel_max']:.0f} in {c25['panel_max_date']}. Data: sunspot number from WDC-SILSO, Royal Observatory of
Belgium, Brussels (version 2.0, {SILSO_DOI}), to {data['last_day']}; the last months are provisional. Sunspot groups:
Royal Greenwich Observatory 1874&ndash;1976 and USAF/NOAA from 1977, compiled by D. H. Hathaway
(solarcyclescience.com); USAF/NOAA areas &times;{USAF_AREA_FACTOR_HATHAWAY} to match RGO's scale. Made by run.py in
15-solar-cycles.</p>
</main>
<script>{get_plotlyjs()}</script>
<script>
(function () {{
  "use strict";
  var DATA = {data_js};
  var THEMES = {themes_js};
  var gd = document.getElementById("plot");
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  var DAY = 86400000;

  function iso(ms) {{ return new Date(ms).toISOString().slice(0, 10); }}
  var d0 = Date.parse(DATA.daily_start + "T00:00:00Z");
  var dx = [], dy = [];
  DATA.daily.forEach(function (v, i) {{ if (v !== null) {{ dx.push(iso(d0 + i * DAY)); dy.push(v); }} }});
  var mx = DATA.monthly.map(function (_, i) {{
    var y = 1749 + Math.floor(i / 12), m = i % 12 + 1;
    return y + "-" + (m < 10 ? "0" : "") + m + "-15";
  }});
  var g0 = Date.parse(DATA.group_epoch + "T00:00:00Z"), G = DATA.groups;
  var gx = G.d.map(function (d) {{ return iso(g0 + d * DAY); }});
  var gbase = G.a.map(function (a) {{ return Math.min(22, 0.55 * Math.sqrt(a)); }});
  var zoomK = 1;
  function gsizes(k) {{ return gbase.map(function (s) {{ return Math.max(2, s * k); }}); }}
  function spanYears(r) {{ return (Date.parse(r[1]) - Date.parse(r[0])) / (365.25 * DAY); }}
  // Dots shrink when zoomed out, so 150 years of groups do not merge into one blob.
  function kFor(r) {{ return Math.max(0.3, Math.min(1, Math.sqrt(18 / spanYears(r)))); }}
  // The top panel's height follows whatever is visible.
  function yMaxIn(r) {{
    var lo = r[0].slice(0, 10), hi = r[1].slice(0, 10), m = 0, i;
    var useDaily = document.getElementById("show-daily").checked;
    if (useDaily) for (i = 0; i < dx.length; i++) {{ if (dx[i] >= lo && dx[i] <= hi && dy[i] > m) m = dy[i]; }}
    for (i = 0; i < mx.length; i++) {{ if (mx[i] >= lo && mx[i] <= hi && DATA.monthly[i] > m) m = DATA.monthly[i]; }}
    return m > 0 ? m * 1.1 + 10 : 100;
  }}
  var gtext = G.id.map(function (id, i) {{
    return id + " · Cycle " + G.c[i] + "<br>" + (G.lat[i] >= 0 ? G.lat[i] + "°N" : -G.lat[i] + "°S") +
      " · largest area " + G.a[i] + " millionths · seen " + G.n[i] + " day" + (G.n[i] > 1 ? "s" : "");
  }});

  function build(th) {{
    var traces = [
      {{ type: "scattergl", mode: "markers", x: dx, y: dy, name: "daily", xaxis: "x", yaxis: "y",
        marker: {{ size: 2.5, color: th.daily }}, hovertemplate: "%{{x}}: %{{y}}<extra>daily</extra>",
        visible: document.getElementById("show-daily").checked }},
      {{ type: "scattergl", mode: "lines", x: mx, y: DATA.monthly, name: "monthly mean", xaxis: "x",
        yaxis: "y", line: {{ color: th.monthly, width: 1 }},
        hovertemplate: "%{{x|%b %Y}}: %{{y:.1f}}<extra>monthly mean</extra>" }},
      {{ type: "scattergl", mode: "lines", x: mx, y: DATA.smoothed, name: "13-month smoothed", xaxis: "x",
        yaxis: "y", line: {{ color: th.smoothed, width: 2.2 }},
        hovertemplate: "%{{x|%b %Y}}: %{{y:.1f}}<extra>13-month smoothed</extra>" }},
      {{ type: "scattergl", mode: "markers", x: gx, y: G.lat, name: "sunspot groups", xaxis: "x",
        yaxis: "y2", text: gtext, hovertemplate: "%{{x}}<br>%{{text}}<extra></extra>",
        marker: {{ size: gsizes(zoomK), color: th.groups, opacity: 0.5, line: {{ width: 0 }} }} }}
    ];
    var ann = DATA.cycles.map(function (c) {{
      var y = Math.floor(c.t_max), mo = Math.round((c.t_max - y) * 12 + 0.5);
      return {{ x: y + "-" + (mo < 10 ? "0" : "") + Math.min(mo, 12) + "-15", y: c.max_sn, xref: "x",
        yref: "y", text: String(c.n), showarrow: false, yshift: 10,
        font: {{ size: 11, color: th.ink2 }} }};
    }});
    var axis = {{ gridcolor: th.grid, linecolor: th.axis, zeroline: false, tickfont: {{ color: th.muted }} }};
    var tfont = {{ color: th.ink2, size: 13 }};
    var layout = {{
      paper_bgcolor: th.surface, plot_bgcolor: th.surface, font: {{ color: th.ink, family:
        "system-ui, -apple-system, Segoe UI, sans-serif" }},
      margin: {{ l: 64, r: 16, t: 16, b: 40 }}, showlegend: true,
      legend: {{ orientation: "h", x: 0, y: 1.06, font: {{ color: th.ink2 }} }},
      hovermode: "closest", dragmode: "zoom",
      xaxis: Object.assign({{ type: "date", domain: [0, 1], anchor: "y2",
        range: ["1745-01-01", DATA.last_day] }}, axis),
      yaxis: Object.assign({{ domain: [0.47, 1], title: {{ text: "sunspot number", font: tfont }}, rangemode: "tozero",
        fixedrange: false }}, axis),
      yaxis2: Object.assign({{ domain: [0, 0.42], title: {{ text: "latitude on the Sun", font: tfont }},
        range: [-50, 50], tickvals: [-40, -20, 0, 20, 40],
        ticktext: ["40°S", "20°S", "0°", "20°N", "40°N"] }}, axis),
      annotations: ann
    }};
    return {{ traces: traces, layout: layout }};
  }}

  function render() {{
    var th = mq.matches ? THEMES.dark : THEMES.light;
    var f = build(th);
    var keep = gd.layout && gd.layout.xaxis ? gd.layout.xaxis.range : null;
    if (keep) f.layout.xaxis.range = keep;
    Plotly.react(gd, f.traces, f.layout, {{ responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d"] }});
  }}
  var FULL = ["1745-01-01", DATA.last_day];
  zoomK = kFor(FULL);
  render();
  if (mq.addEventListener) mq.addEventListener("change", render);

  function applyZoom(r) {{
    zoomK = kFor(r);
    Plotly.restyle(gd, {{ "marker.size": [gsizes(zoomK)] }}, [3]);
    Plotly.relayout(gd, {{ "yaxis.range": [0, yMaxIn(r)] }});
  }}
  var busy = false;
  gd.on("plotly_relayout", function (ev) {{
    if (busy) return;
    var r = null;
    if (ev["xaxis.range[0]"] !== undefined) r = [String(ev["xaxis.range[0]"]), String(ev["xaxis.range[1]"])];
    else if (ev["xaxis.range"]) r = ev["xaxis.range"].map(String);
    else if (ev["xaxis.autorange"]) r = FULL;
    if (!r) return;
    busy = true;
    Promise.resolve(applyZoom(r)).then(function () {{ busy = false; }}, function () {{ busy = false; }});
    setTimeout(function () {{ busy = false; }}, 400);
  }});
  document.querySelectorAll("button[data-range]").forEach(function (b) {{
    b.addEventListener("click", function () {{
      Plotly.relayout(gd, {{ "xaxis.range": b.getAttribute("data-range").split(",") }});
    }});
  }});
  document.getElementById("show-daily").addEventListener("change", function (e) {{
    Plotly.restyle(gd, {{ visible: e.target.checked }}, [0]);
    applyZoom(gd.layout.xaxis.range.map(String));
  }});
}})();
</script>
</body>
</html>
"""
    out = FIGURE_DIR / "solar_cycle_explorer.html"
    out.write_text(page)
    return str(out)
