"""The interactive orbit viewer: one self-contained HTML page.

Every position in it is a Horizons sample -- nothing is propagated or
extrapolated in the browser. The page steps a clock through 1972-2050 and, for
each object, interpolates linearly between the two stored samples either side
of the current date.

To keep the page light, each object's track is thinned adaptively before it is
embedded: a sample is dropped only if straight-line interpolation between its
neighbours still passes within 0.3% of the object's distance from the Sun of
where Horizons actually put it. Straight cruise legs shrink to a handful of
points; Halley's hairpin at perihelion keeps many.

The page follows the viewer's light/dark setting, and carries a live table of
every object's distance and speed, so nothing can only be read by hovering.
"""

from __future__ import annotations

import json

import numpy as np

import horizons as hz
import orbit_physics as op
from config import (
    COMETS, DWARFS, INTERSTELLAR, PLANETS, PUBLISHED, SNAPSHOT_DATE, SPACECRAFT, VIEWER_PATH,
    VIEWER_SPAN, body,
)

#: Launch dates (UTC), from the Horizons headers and NASA mission pages.
LAUNCHES = {"pioneer10": "1972-03-03", "pioneer11": "1973-04-06", "voyager2": "1977-08-20",
            "voyager1": "1977-09-05", "new_horizons": "2006-01-19"}

#: Thinning tolerance: fraction of the distance from the Sun, plus a floor in au.
TOL_REL = 0.003
TOL_ABS = 0.0005


def _jd(date: str) -> float:
    from astropy.time import Time

    return float(Time(date, scale="tdb").jd)


T0, T1 = _jd(VIEWER_SPAN[0]), _jd(VIEWER_SPAN[1])


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------

def _series(b) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Full-resolution (jd, r, v) for a body over the viewer's span."""
    if b.group in ("planet", "dwarf"):
        d = hz.load(hz.orbit_key(b))
        jd, r, v = d["jd"], d["r"], d["v"]
    elif b.group == "interstellar":
        tr, pas = hz.load(hz.track_key(b)), hz.load(f"passage_{b.key}")
        before, after = tr["jd"] < pas["jd"][0], tr["jd"] > pas["jd"][-1]
        jd = np.r_[tr["jd"][before], pas["jd"], tr["jd"][after]]
        r = np.vstack([tr["r"][before], pas["r"], tr["r"][after]])
        v = np.vstack([tr["v"][before], pas["v"], tr["v"][after]])
    else:
        d = hz.load(hz.track_key(b))
        jd, r, v = d["jd"], d["r"], d["v"]
    # Keep one sample either side of the span so interpolation reaches the ends.
    lo = max(np.searchsorted(jd, T0) - 1, 0)
    hi = min(np.searchsorted(jd, T1) + 1, len(jd))
    return jd[lo:hi], r[lo:hi], v[lo:hi]


def _segment_error(t, xyz, i, j):
    """Worst distance between samples i..j and the straight line (in time) from i to j."""
    if j - i < 2:
        return 0.0
    f = ((t[i + 1:j] - t[i]) / (t[j] - t[i]))[:, None]
    interp = xyz[i] + f * (xyz[j] - xyz[i])
    return float(np.max(op.norm(xyz[i + 1:j] - interp)))


def thin(t, xyz) -> np.ndarray:
    """Indices of the samples to keep (greedy, with exponential then binary search)."""
    n = len(t)
    keep = [0]
    i = 0
    dist = op.norm(xyz)
    while i < n - 1:
        tol = lambda j: TOL_ABS + TOL_REL * min(dist[i], dist[j])  # noqa: E731
        good, step = i + 1, 1
        bad = None
        while True:
            j = min(i + step, n - 1)
            if _segment_error(t, xyz, i, j) <= tol(j):
                good = j
                if j == n - 1:
                    break
                step *= 2
            else:
                bad = j
                break
        if bad is not None:
            lo, hi = good, bad
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if _segment_error(t, xyz, i, mid) <= tol(mid):
                    lo = mid
                else:
                    hi = mid
            good = lo
        keep.append(good)
        i = good
    return np.array(keep)


def _body_payload(b) -> dict:
    jd, r, v = _series(b)
    idx = thin(jd, r)
    speed = op.to_kms(op.norm(v))
    out = {
        "key": b.key, "name": b.name, "group": b.group,
        "t": [round(float(x - T0), 3) for x in jd[idx]],
        "p": [round(float(x), 4) for x in r[idx].ravel()],
        "s": [round(float(x), 2) for x in speed[idx]],
    }
    if b.group in ("planet", "dwarf"):
        snap = _jd(SNAPSHOT_DATE)
        m = np.abs(jd - snap) <= b.period_days / 2
        loop = r[m]
        li = thin(jd[m], loop)
        loop = np.vstack([loop[li], loop[li][:1]])
        out["orbit"] = [round(float(x), 4) for x in loop.ravel()]
    return out


def _events(res) -> list[dict]:
    """Moments worth jumping to, all dated by the analysis or the headers."""
    ev = [{"date": d, "label": f"{body(k).name} launched", "view": "planets"}
          for k, d in LAUNCHES.items()]
    for f in res["spacecraft"]["flybys"]:
        ev.append({"date": f["closest_approach_utc"][:10],
                   "label": f"{body(f['craft']).name} flies past {f['planet'].capitalize()}",
                   "view": "planets"})
    k2 = res["kepler2"]
    ev.append({"date": k2["perihelion_1_utc"][:10], "label": "Halley's comet at perihelion",
               "view": "inner"})
    ev.append({"date": k2["aphelion_utc"][:10], "label": "Halley's comet at aphelion, 35 au out",
               "view": "planets"})
    hb = hz.load(hz.track_key(body("hale_bopp")))
    rn, rdot = op.radial_distance_and_rate(hb["r"], hb["v"])
    t_hb, q_hb = op.extrema(hb["jd"], rn, rdot, "min")[0]
    from analysis import tdb_to_utc_iso

    ev.append({"date": tdb_to_utc_iso(t_hb)[:10],
               "label": f"Comet Hale-Bopp at perihelion ({q_hb:.2f} au)", "view": "inner"})
    for key, o in res["interstellar"]["objects"].items():
        ev.append({"date": o["perihelion_utc"][:10],
                   "label": f"{o['name']} at perihelion", "view": "inner"})
    for key, h in PUBLISHED["heliopause"].items():
        ev.append({"date": h["date"], "label": f"{body(key).name} leaves the heliosphere",
                   "view": "craft"})
    ev.append({"date": "2015-07-14", "label": "New Horizons flies past Pluto", "view": "planets"})
    ev.append({"date": SNAPSHOT_DATE, "label": "Today (the project's snapshot date)",
               "view": "visitors"})
    ev.append({"date": res["light_day"]["crossing_utc"][:10],
               "label": "Voyager 1 one light-day from Earth", "view": "craft"})
    return sorted(ev, key=lambda e: e["date"])


def _days(date: str) -> float:
    return round(_jd(date) - T0, 3)


def build_payload(res) -> dict:
    bodies = PLANETS + DWARFS + COMETS + INTERSTELLAR + SPACECRAFT
    events = _events(res)
    for e in events:
        e["t"] = _days(e["date"])
    return {
        "t0_iso": VIEWER_SPAN[0], "span_days": round(T1 - T0, 3),
        "today": _days(SNAPSHOT_DATE), "today_iso": SNAPSHOT_DATE,
        "t0_jd": T0,
        "bodies": [_body_payload(b) for b in bodies],
        "events": events,
    }


# --------------------------------------------------------------------------
# The page
# --------------------------------------------------------------------------

PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Solar System Orbits</title>
<style>
  :root {
    color-scheme: light;
    --surface: #fcfcfb; --plane: #f9f9f7; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
    --grid: #e1e0d9; --axis: #c3c2b7; --ring: rgba(11,11,11,0.10);
    --craft: #2a78d6; --iso: #eb6834; --comet: #1baf7a;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      color-scheme: dark;
      --surface: #1a1a19; --plane: #0d0d0d; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
      --grid: #2c2c2a; --axis: #383835; --ring: rgba(255,255,255,0.10);
      --craft: #3987e5; --iso: #d95926; --comet: #199e70;
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --surface: #1a1a19; --plane: #0d0d0d; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --ring: rgba(255,255,255,0.10);
    --craft: #3987e5; --iso: #d95926; --comet: #199e70;
  }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--plane); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
  main { max-width: 1280px; margin: 0 auto; padding: 20px 16px 40px; }
  h1 { font-size: 1.35rem; margin: 0 0 4px; }
  p { margin: 0 0 10px; color: var(--ink2); max-width: 85ch; }
  .card { background: var(--surface); border: 1px solid var(--ring); border-radius: 10px; }
  .controls { display: flex; flex-wrap: wrap; gap: 10px 18px; align-items: center;
              padding: 12px 14px; margin: 12px 0 10px; }
  .date { font-size: 1.6rem; font-weight: 600; min-width: 9.5ch;
          font-variant-numeric: tabular-nums; }
  .slider { flex: 1 1 320px; display: flex; flex-direction: column; gap: 2px; min-width: 220px; }
  .slider input { width: 100%; accent-color: var(--craft); }
  .slider .ends { display: flex; justify-content: space-between; font-size: 0.75rem;
                  color: var(--muted); }
  button, select { font: inherit; font-size: 0.85rem; color: var(--ink); background: var(--surface);
                   border: 1px solid var(--axis); border-radius: 6px; padding: 4px 10px;
                   cursor: pointer; }
  button:hover, select:hover { border-color: var(--ink2); }
  button[aria-pressed="true"] { background: var(--ink); color: var(--surface);
                                border-color: var(--ink); }
  .group { display: flex; flex-wrap: wrap; gap: 6px; align-items: center; }
  .group .lbl { font-size: 0.75rem; color: var(--muted); margin-right: 2px; }
  .layout { display: grid; grid-template-columns: minmax(0, 1fr) 330px; gap: 12px; }
  #plot { width: 100%; height: min(74vh, 760px); min-height: 420px; }
  .readout { padding: 10px 12px; max-height: min(74vh, 760px); overflow: auto; }
  .readout h2 { font-size: 0.9rem; margin: 2px 0 6px; }
  table { border-collapse: collapse; width: 100%; font-size: 0.8rem; }
  th { text-align: left; font-weight: 600; color: var(--muted); font-size: 0.7rem;
       padding: 6px 4px 2px; text-transform: uppercase; letter-spacing: 0.04em; }
  td { padding: 2px 4px; border-top: 1px solid var(--grid); }
  td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
  .key { display: inline-block; width: 12px; height: 2px; vertical-align: middle;
         margin-right: 6px; border-radius: 1px; }
  .note { font-size: 0.8rem; color: var(--muted); margin-top: 10px; max-width: 110ch; }
  @media (max-width: 900px) {
    .layout { grid-template-columns: 1fr; }
    .readout { max-height: none; }
    #plot { height: 64vh; }
  }
</style>
</head>
<body>
<main>
<h1>Solar system orbits, 1972–2050</h1>
<p>Real positions from JPL Horizons for the planets, Pluto and Ceres, Halley's comet and
Hale-Bopp, the three interstellar objects, and the five spacecraft leaving the solar system.
Drag to rotate, scroll or pinch to zoom, hover for details. Every position is an ephemeris
sample (predictions after __TODAY__), not an animation.</p>

<div class="card controls">
  <div class="date" id="date" aria-live="polite"></div>
  <div class="slider">
    <input type="range" id="slider" min="0" max="__SPAN__" step="1" aria-label="Date">
    <div class="ends"><span>1972</span><span>1990</span><span>2010</span><span>2030</span><span>2050</span></div>
  </div>
  <div class="group">
    <button id="play" aria-pressed="false">Play</button>
    <select id="speed" aria-label="Playback speed">
      <option value="30">1 month per second</option>
      <option value="182.6">6 months per second</option>
      <option value="730.5" selected>2 years per second</option>
      <option value="1826">5 years per second</option>
    </select>
    <button id="today">Today</button>
  </div>
  <div class="group" role="group" aria-label="View">
    <span class="lbl">View</span>
    <button data-view="inner">Inner planets</button>
    <button data-view="planets">All planets</button>
    <button data-view="visitors">Visitors</button>
    <button data-view="craft">Spacecraft</button>
  </div>
  <div class="group">
    <select id="events" aria-label="Jump to a moment"><option value="">Jump to a moment…</option></select>
  </div>
</div>

<div class="layout">
  <div class="card"><div id="plot"></div></div>
  <div class="card readout">
    <h2>On <span id="date2"></span></h2>
    <table id="table"><tbody></tbody></table>
  </div>
</div>

<p class="note">Sun-centred coordinates, ecliptic of J2000 (x towards the March equinox), in au
(1 au = 149.6 million km). Planet orbits are drawn for one period around __TODAY__. Pale lines
are each object's whole path within 1972–2050; bold lines, its last five years. Positions between
stored samples are interpolated in straight lines; samples are kept so that this stays within
0.3% of the true distance. Through their planetary flybys the Voyager trajectories are
patched-conic reconstructions (JPL's own caveat), and every spacecraft position after its last
tracking data is a prediction. Made by run.py in 12-orbit-visualizer.</p>
</main>

<script>__PLOTLY__</script>
<script>
const DATA = __DATA__;
const VIEWS = {
  inner: { range: 2.2, eye: { x: 0.2, y: -0.95, z: 1.45 } },
  planets: { range: 36, eye: { x: 0.25, y: -0.95, z: 1.4 } },
  visitors: { range: 72, eye: { x: 0.3, y: -1.0, z: 1.35 } },
  craft: { range: 190, eye: { x: 0.9, y: -1.1, z: 0.9 } },
};
const GROUPS = [
  { id: "planet", name: "Planets", css: "--ink" },
  { id: "dwarf", name: "Dwarf planets", css: "--ink" },
  { id: "comet", name: "Comets", css: "--comet" },
  { id: "interstellar", name: "Interstellar objects", css: "--iso" },
  { id: "spacecraft", name: "Spacecraft", css: "--craft" },
];
const LEGEND = { planet: "Planets and dwarf planets", dwarf: null, comet: "Comets",
                 interstellar: "Interstellar objects", spacecraft: "Spacecraft" };

const bodies = DATA.bodies;
const gd = document.getElementById("plot");
const slider = document.getElementById("slider");
let state = { t: DATA.today, view: "visitors", playing: false, speed: 730.5, last: null };

function css(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}
function theme() {
  return { surface: css("--surface"), ink: css("--ink"), ink2: css("--ink2"),
           muted: css("--muted"), grid: css("--grid"), axis: css("--axis"),
           craft: css("--craft"), iso: css("--iso"), comet: css("--comet") };
}
function colorFor(group, th) {
  return { planet: th.ink, dwarf: th.ink, comet: th.comet, interstellar: th.iso,
           spacecraft: th.craft }[group];
}

// Where is body b at time t? Binary search its samples, then interpolate.
function locate(b, t) {
  const ts = b.t, n = ts.length;
  if (t < ts[0] || t > ts[n - 1]) return null;
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (ts[mid] <= t) lo = mid; else hi = mid; }
  const f = ts[hi] === ts[lo] ? 0 : (t - ts[lo]) / (ts[hi] - ts[lo]);
  const p = b.p, i = 3 * lo, j = 3 * hi;
  return { x: p[i] + f * (p[j] - p[i]), y: p[i + 1] + f * (p[j + 1] - p[i + 1]),
           z: p[i + 2] + f * (p[j + 2] - p[i + 2]), s: b.s[lo] + f * (b.s[hi] - b.s[lo]),
           idx: lo };
}

// Plotly draws 3D lines even outside the axis box, so vertices beyond the
// current view's cube are dropped (the line breaks there instead).
function inBox(px, py, pz) {
  const R = VIEWS[state.view].range;
  return Math.abs(px) <= R && Math.abs(py) <= R && Math.abs(pz) <= R;
}
function pushPoint(out, px, py, pz) {
  if (inBox(px, py, pz)) { out.x.push(px); out.y.push(py); out.z.push(pz); }
  else if (out.x.length && out.x[out.x.length - 1] !== null) { out.x.push(null); out.y.push(null); out.z.push(null); }
}
// Trails show the last five years of motion: enough to see direction and pace
// without old launch arcs cluttering the inner system decades later.
const TRAIL_DAYS = 1826;
function pathUpTo(list, t, full) {
  const out = { x: [], y: [], z: [] };
  for (const b of list) {
    const n = b.t.length;
    let start = 0, end = n;
    if (!full) {
      if (t < b.t[0]) continue;
      end = Math.min(n, (locate(b, t) || { idx: n - 1 }).idx + 1);
      const from = locate(b, Math.max(b.t[0], t - TRAIL_DAYS));
      if (from) { start = from.idx + 1; pushPoint(out, from.x, from.y, from.z); }
    }
    for (let k = start; k < end; k++) pushPoint(out, b.p[3 * k], b.p[3 * k + 1], b.p[3 * k + 2]);
    if (!full) { const now = locate(b, t); if (now) pushPoint(out, now.x, now.y, now.z); }
    out.x.push(null); out.y.push(null); out.z.push(null);
  }
  return out;
}

const byGroup = g => bodies.filter(b => b.group === g);
const moving = ["interstellar", "spacecraft"];

function staticTraces(th) {
  const traces = [];
  // Sun
  traces.push({ type: "scatter3d", mode: "markers+text", x: [0], y: [0], z: [0], text: ["Sun"],
    textposition: "top center", textfont: { color: th.ink2, size: 11 },
    marker: { size: 5, color: th.ink, line: { color: th.surface, width: 2 } },
    hovertemplate: "Sun<extra></extra>", showlegend: false, name: "Sun" });
  // Planet and dwarf-planet orbits
  const orb = { x: [], y: [], z: [] };
  for (const b of bodies.filter(b => b.orbit)) {
    for (let k = 0; k < b.orbit.length; k += 3) pushPoint(orb, b.orbit[k], b.orbit[k + 1], b.orbit[k + 2]);
    orb.x.push(null); orb.y.push(null); orb.z.push(null);
  }
  traces.push({ type: "scatter3d", mode: "lines", x: orb.x, y: orb.y, z: orb.z, hoverinfo: "skip",
    line: { color: th.axis, width: 2 }, showlegend: false, name: "orbits", legendgroup: "planet" });
  // Whole paths, pale
  for (const g of ["comet", "interstellar", "spacecraft"]) {
    const pth = pathUpTo(byGroup(g), 0, true);
    traces.push({ type: "scatter3d", mode: "lines", ...pth, hoverinfo: "skip", opacity: 0.35,
      line: { color: colorFor(g, th), width: 2 }, showlegend: false, legendgroup: g, name: g + " path" });
  }
  return traces;
}

function dynamicTraces(th, t) {
  const traces = [];
  for (const g of moving) {
    const pth = pathUpTo(byGroup(g), t, false);
    traces.push({ type: "scatter3d", mode: "lines", ...pth, hoverinfo: "skip",
      line: { color: colorFor(g, th), width: 4 }, showlegend: false, legendgroup: g, name: g + " trail" });
  }
  for (const g of ["planet", "dwarf", "comet", "interstellar", "spacecraft"]) {
    traces.push(markerTrace(g, th, t));
  }
  return traces;
}

function markerTrace(g, th, t) {
  const list = byGroup(g), R = VIEWS[state.view].range;
  const x = [], y = [], z = [], text = [], cd = [];
  for (const b of list) {
    const p = locate(b, t);
    if (!p) { x.push(null); y.push(null); z.push(null); text.push(""); cd.push([b.name, 0, 0]); continue; }
    const r = Math.hypot(p.x, p.y, p.z);
    const inside = inBox(p.x, p.y, p.z);
    x.push(inside ? p.x : null); y.push(inside ? p.y : null); z.push(inside ? p.z : null);
    // Name an object only where it is resolved at this scale.
    text.push(inside && r > 0.045 * R ? b.name : "");
    cd.push([b.name, r, p.s]);
  }
  return { type: "scatter3d", mode: "markers+text", x, y, z, text, customdata: cd,
    textposition: "top center", textfont: { color: th.ink2, size: 11 },
    marker: { size: g === "planet" ? 5 : 4.5, color: colorFor(g, th),
              symbol: g === "dwarf" ? "circle-open" : "circle",
              line: { color: th.surface, width: 1 } },
    hovertemplate: "<b>%{customdata[1]:.2f} au</b> from the Sun<br>%{customdata[2]:.1f} km/s<br>%{customdata[0]}<extra></extra>",
    name: LEGEND[g] || "", showlegend: LEGEND[g] !== null, legendgroup: g === "dwarf" ? "planet" : g };
}

function axis(th, title, R) {
  return { title: { text: title, font: { color: th.muted, size: 11 } }, range: [-R, R],
           showgrid: true, gridcolor: th.grid, zeroline: false, showbackground: false,
           showspikes: false, color: th.muted, tickfont: { color: th.muted, size: 10 },
           linecolor: th.axis };
}

function layout(th) {
  const v = VIEWS[state.view];
  return {
    paper_bgcolor: th.surface, plot_bgcolor: th.surface, margin: { l: 0, r: 0, t: 0, b: 0 },
    font: { family: "system-ui, -apple-system, 'Segoe UI', sans-serif", color: th.ink, size: 12 },
    showlegend: true, legend: { x: 0.01, y: 0.99, bgcolor: "rgba(0,0,0,0)",
                                font: { size: 11, color: th.ink2 } },
    hoverlabel: { bgcolor: th.surface, bordercolor: th.axis, font: { color: th.ink, size: 12 } },
    uirevision: state.view,
    scene: { aspectmode: "cube", bgcolor: th.surface,
             xaxis: axis(th, "x (au)", v.range), yaxis: axis(th, "y (au)", v.range),
             zaxis: axis(th, "z (au)", v.range),
             camera: { eye: v.eye, up: { x: 0, y: 0, z: 1 } } },
  };
}

let nStatic = 0;
function render(full) {
  const th = theme();
  if (full) {
    const st = staticTraces(th);
    nStatic = st.length;
    Plotly.react(gd, st.concat(dynamicTraces(th, state.t)), layout(th),
                 { displaylogo: false, responsive: true });
  } else {
    const dyn = dynamicTraces(th, state.t);
    const idx = dyn.map((_, i) => nStatic + i);
    Plotly.restyle(gd, { x: dyn.map(d => d.x), y: dyn.map(d => d.y), z: dyn.map(d => d.z),
                         text: dyn.map(d => d.text || null),
                         customdata: dyn.map(d => d.customdata || null) }, idx);
  }
  updateText();
}

// ---- dates and the readout table ----
const DAY_MS = 86400000;
const T0 = Date.parse(DATA.t0_iso + "T00:00:00Z");
function fmt(t, long) {
  const d = new Date(T0 + t * DAY_MS);
  return d.toLocaleDateString("en-GB", { day: "numeric", month: long ? "long" : "short",
                                         year: "numeric", timeZone: "UTC" });
}
const tbody = document.querySelector("#table tbody");
const cells = {};
function buildTable() {
  const th = document.createElement("tr");
  for (const [txt, cls] of [["Object", ""], ["From the Sun", "num"], ["Speed", "num"]]) {
    const h = document.createElement("th"); h.textContent = txt; if (cls) h.className = cls; th.appendChild(h);
  }
  tbody.appendChild(th);
  for (const g of GROUPS) {
    const list = byGroup(g.id);
    const hr = document.createElement("tr"), hc = document.createElement("th");
    hc.colSpan = 3; hc.textContent = g.name; hr.appendChild(hc); tbody.appendChild(hr);
    for (const b of list) {
      const tr = document.createElement("tr");
      const n = document.createElement("td");
      const key = document.createElement("span"); key.className = "key";
      key.style.background = `var(${g.css})`;
      n.appendChild(key); n.appendChild(document.createTextNode(b.name));
      const d = document.createElement("td"); d.className = "num";
      const s = document.createElement("td"); s.className = "num";
      tr.append(n, d, s); tbody.appendChild(tr);
      cells[b.key] = { d, s };
    }
  }
}
function updateText() {
  document.getElementById("date").textContent = fmt(state.t, false);
  document.getElementById("date2").textContent = fmt(state.t, true);
  slider.value = Math.round(state.t);
  for (const b of bodies) {
    const p = locate(b, state.t), c = cells[b.key];
    if (!p) { c.d.textContent = "—"; c.s.textContent = "—"; continue; }
    const r = Math.hypot(p.x, p.y, p.z);
    c.d.textContent = (r < 10 ? r.toFixed(3) : r.toFixed(1)) + " au";
    c.s.textContent = p.s.toFixed(1) + " km/s";
  }
}

// ---- controls ----
slider.addEventListener("input", () => { state.t = +slider.value; render(false); });
document.getElementById("today").addEventListener("click", () => { state.t = DATA.today; render(false); });
document.getElementById("speed").addEventListener("change", e => { state.speed = +e.target.value; });
const playBtn = document.getElementById("play");
playBtn.addEventListener("click", () => {
  state.playing = !state.playing;
  playBtn.setAttribute("aria-pressed", state.playing);
  playBtn.textContent = state.playing ? "Pause" : "Play";
  state.last = null;
  if (state.playing) {
    if (state.t >= DATA.span_days - 1) state.t = 0;
    requestAnimationFrame(tick);
  }
});
function tick(ts) {
  if (!state.playing) return;
  if (state.last !== null) {
    state.t = Math.min(DATA.span_days, state.t + state.speed * (ts - state.last) / 1000);
    render(false);
    if (state.t >= DATA.span_days) { playBtn.click(); return; }
  }
  state.last = ts;
  requestAnimationFrame(tick);
}
function setView(v) {
  state.view = v;
  for (const b of document.querySelectorAll("[data-view]")) b.setAttribute("aria-pressed", b.dataset.view === v);
  render(true);
}
for (const b of document.querySelectorAll("[data-view]")) b.addEventListener("click", () => setView(b.dataset.view));
const evSel = document.getElementById("events");
for (const e of DATA.events) {
  // Pioneer 10's trajectory starts the day after its launch: snap to the first day.
  if (e.t < -3 || e.t > DATA.span_days) continue;
  const o = document.createElement("option");
  o.value = Math.max(0, e.t); o.dataset.view = e.view;
  o.textContent = `${e.date} — ${e.label}`;
  evSel.appendChild(o);
}
evSel.addEventListener("change", () => {
  const o = evSel.selectedOptions[0];
  if (!o.value) return;
  state.t = +o.value;
  if (o.dataset.view && o.dataset.view !== state.view) setView(o.dataset.view); else render(false);
});
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => render(true));

// A moment can be linked: #date=1979-07-09&view=planets
function readHash() {
  const h = new URLSearchParams(location.hash.slice(1));
  const d = h.get("date"), v = h.get("view");
  if (d && !isNaN(Date.parse(d + "T00:00:00Z"))) {
    state.t = Math.max(0, Math.min(DATA.span_days, (Date.parse(d + "T00:00:00Z") - T0) / DAY_MS));
  }
  if (v && VIEWS[v]) state.view = v;
}
function writeHash() {
  const d = new Date(T0 + state.t * DAY_MS).toISOString().slice(0, 10);
  history.replaceState(null, "", `#date=${d}&view=${state.view}`);
}
for (const el of [slider, evSel]) el.addEventListener("change", writeHash);
document.getElementById("today").addEventListener("click", writeHash);
for (const b of document.querySelectorAll("[data-view]")) b.addEventListener("click", writeHash);
playBtn.addEventListener("click", () => { if (!state.playing) writeHash(); });

buildTable();
readHash();
setView(state.view);
</script>
</body>
</html>
"""


def build_viewer(res) -> str:
    from plotly.offline import get_plotlyjs

    payload = build_payload(res)
    data = json.dumps(payload, separators=(",", ":"))
    html = (PAGE.replace("__PLOTLY__", get_plotlyjs())
            .replace("__DATA__", data)
            .replace("__SPAN__", str(int(payload["span_days"])))
            .replace("__TODAY__", SNAPSHOT_DATE))
    VIEWER_PATH.write_text(html, encoding="utf-8")
    return str(VIEWER_PATH), len(data)
