"""The S2 orbit explorer: one self-contained HTML page, no external libraries.

A clock runs from 1992 (the first NTT image) to 2035, past the next pericentre
in 2034. At each moment the page shows where S2 is on the sky, its
line-of-sight velocity, its distance from the black hole and its speed, and a
table of the numbers -- so nothing can only be read by hovering.

Everything shown is our best-fitting joint (Keck + VLT) orbit with general
relativity, sampled evenly in *eccentric anomaly*, which automatically puts
most samples near pericentre where the star moves fastest; the page only
interpolates between those samples. The measured points appear as the clock
passes the night they were taken. After 2019 the curve is a prediction.

The page follows the viewer's light or dark setting.
"""

from __future__ import annotations

import json

import numpy as np

import orbit_physics as op
from config import VIEWER_PATH

T_START, T_END = 1992.0, 2035.0


def _track(p: dict) -> dict:
    """Model quantities on a grid that is even in eccentric anomaly."""
    n = 2 * np.pi / p["period"]
    m0, m1 = n * (T_START - p["t0"]), n * (T_END - p["t0"])
    e0, _ = op.solve_kepler(np.array([m0]), p["e"])
    e1, _ = op.solve_kepler(np.array([m1]), p["e"])
    big_e = np.linspace(e0[0], e1[0], 2400)
    t = p["t0"] + (big_e - p["e"] * np.sin(big_e)) / n
    q = {k: p[k] for k in ("gm6", "r0", "period", "t0", "e", "inc", "omega", "Omega", "f", "f_sp")}
    pred = op.predict(q, t, t, roemer=True)
    el = (p["gm6"], p["period"], p["t0"], p["e"], p["inc"], p["omega"], p["Omega"])
    t_em = op.roemer_emission_time(t, *el, f_sp=p["f_sp"])
    _, _, z, vx, vy, vz, r = op.orbit_state(t_em, *el, f_sp=p["f_sp"])
    speed = np.sqrt(vx**2 + vy**2 + vz**2)
    red = op.relativistic_redshift_kms(p["gm6"], r, speed)
    return {
        "t": np.round(t, 5).tolist(),
        "x": np.round(pred["x"], 3).tolist(),
        "y": np.round(pred["y"], 3).tolist(),
        "rv": np.round(pred["rv"] + p["vz0"], 1).tolist(),
        "r": np.round(r, 2).tolist(),
        "speed": np.round(speed, 1).tolist(),
        "redshift": np.round(red, 2).tolist(),
        "delay": np.round(z / op.C_AUYR * 365.25, 3).tolist(),
    }


def _data(art: dict, p: dict) -> dict:
    model = art["models"]["joint_gr"]
    out = {"astro": [], "rv": []}
    for fr in model.frames:
        dt = fr.data.t - fr.t_ref
        k = fr.key
        x = fr.data.x - p[f"x0_{k}"] - p[f"vx0_{k}"] * dt
        y = fr.data.y - p[f"y0_{k}"] - p[f"vy0_{k}"] * dt
        for i in range(fr.data.n):
            out["astro"].append([round(float(fr.data.t[i]), 4), round(float(x[i]), 2), round(float(y[i]), 2),
                                 round(float(fr.data.x_err[i]), 2), round(float(fr.data.y_err[i]), 2),
                                 k, str(fr.data.instrument[i])])
    rv = model.rv
    corr = np.where(rv.instrument == "NIRC2", p.get("nirc2_offset", 0.0), 0.0)
    for i in range(rv.n):
        team = "vlt" if rv.instrument[i] == "VLT" else "keck"
        out["rv"].append([round(float(rv.t[i]), 4), round(float(rv.v[i] + corr[i]), 1), float(rv.v_err[i]), team,
                          str(rv.instrument[i])])
    return out


def build_viewer(results: dict, art: dict):
    fit = art["fits"]["joint_gr"]
    p = dict(zip(fit["names"], fit["best_theta"])) | dict(art["models"]["joint_gr"].fixed)
    d = results["fits"]["joint_gr"]
    payload = {
        "track": _track(p),
        "data": _data(art, p),
        "meta": {
            "gm6": round(float(p["gm6"]), 4), "r0": round(float(p["r0"]), 4),
            "gm6_err": round(0.5 * (d["parameters"]["gm6"]["hi68"] - d["parameters"]["gm6"]["lo68"]), 4),
            "r0_err": round(0.5 * (d["parameters"]["r0"]["hi68"] - d["parameters"]["r0"]["lo68"]), 4),
            "t0": float(p["t0"]), "period": float(p["period"]),
            "rs_au": float(op.schwarzschild_radius_au(p["gm6"])),
            "light_hour_au": op.C_KMS * 3600 / op.AU_KM,
            "c_kms": op.C_KMS,
        },
    }
    blob = json.dumps(payload, separators=(",", ":"))
    html = PAGE.replace("__DATA__", blob)
    VIEWER_PATH.write_text(html, encoding="utf-8")
    return VIEWER_PATH, len(blob)


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S2 Orbit Explorer</title>
<style>
:root {
  --bg: #fcfcfb; --panel: #ffffff; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --keck: #2a78d6; --vlt: #eb6834; --gr: #1baf7a;
  --star: #0b0b0b; --orbit: #0b0b0b; --future: #c3c2b7;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
    --grid: #2c2c2a; --axis: #4a4945; --keck: #5a9be6; --vlt: #f08257; --gr: #34c48f;
    --star: #f2f1ec; --orbit: #f2f1ec; --future: #4a4945;
  }
}
:root[data-theme="dark"] {
  --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
  --grid: #2c2c2a; --axis: #4a4945; --keck: #5a9be6; --vlt: #f08257; --gr: #34c48f;
  --star: #f2f1ec; --orbit: #f2f1ec; --future: #4a4945;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
.wrap { max-width: 1180px; margin: 0 auto; padding: 18px 16px 28px; }
h1 { font-size: 20px; margin: 0 0 2px; }
.sub { color: var(--ink2); margin: 0 0 12px; max-width: 900px; }
.controls { display: flex; flex-wrap: wrap; gap: 10px 16px; align-items: center; margin: 8px 0 12px; }
.controls button, .controls select { font: inherit; color: var(--ink); background: var(--panel);
  border: 1px solid var(--axis); border-radius: 6px; padding: 4px 10px; cursor: pointer; }
.controls button[aria-pressed="true"] { border-color: var(--ink); }
.controls label { color: var(--ink2); display: inline-flex; gap: 6px; align-items: center; }
#slider { width: min(420px, 90vw); }
#date { font-variant-numeric: tabular-nums; min-width: 9.5em; font-weight: 600; }
.grid { display: grid; grid-template-columns: minmax(0, 1.05fr) minmax(0, 1fr); gap: 14px; }
@media (max-width: 860px) { .grid { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border: 1px solid var(--grid); border-radius: 8px; padding: 8px 10px; }
.card h2 { font-size: 13px; margin: 2px 0 4px; color: var(--ink2); font-weight: 600; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
svg text { fill: var(--muted); font-size: 11px; }
svg .lbl { fill: var(--ink2); }
.readout { width: 100%; border-collapse: collapse; margin-top: 12px; font-variant-numeric: tabular-nums; }
.readout td { padding: 3px 8px; border-bottom: 1px solid var(--grid); }
.readout td:first-child { color: var(--ink2); width: 42%; }
.legend { display: flex; flex-wrap: wrap; gap: 4px 14px; color: var(--ink2); font-size: 12px; margin: 2px 0 4px; }
.sw { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 5px; vertical-align: -1px; }
.note { color: var(--muted); font-size: 12px; margin-top: 14px; max-width: 1000px; }
#tip { position: fixed; pointer-events: none; background: var(--panel); color: var(--ink);
  border: 1px solid var(--axis); border-radius: 6px; padding: 4px 8px; font-size: 12px; display: none;
  box-shadow: 0 2px 8px rgba(0,0,0,.12); }
</style>
</head>
<body>
<div class="wrap">
<h1>S2 around Sagittarius A*</h1>
<p class="sub">A star on a 16-year orbit around the black hole at the centre of the Milky Way, from the positions and
velocities measured at Keck (UCLA) and the VLT (MPE). Drag the slider or press play; hover a chart to scrub
through time. The curve is our best-fitting orbit with general relativity; after 2019 it is a prediction.</p>
<div class="controls">
  <button id="play" aria-pressed="false">▶ Play</button>
  <input id="slider" type="range" min="1992" max="2035" step="0.001" value="2018.38" aria-label="date">
  <span id="date"></span>
  <label>speed <select id="speed"><option value="0.5">½ yr/s</option><option value="2" selected>2 yr/s</option>
    <option value="0.05">≈ 1 week/s</option></select></label>
  <button id="zoom" aria-pressed="false">Zoom to Sgr A*</button>
  <label><input type="checkbox" id="showKeck" checked> Keck</label>
  <label><input type="checkbox" id="showVlt" checked> VLT</label>
</div>
<div class="grid">
  <div class="card">
    <h2>On the sky (East to the left, as the sky looks)</h2>
    <div class="legend"><span><span class="sw" style="background:var(--keck)"></span>Keck positions</span>
      <span><span class="sw" style="background:var(--vlt)"></span>NTT/VLT positions</span>
      <span><span class="sw" style="background:var(--ink)"></span>S2 now</span>
      <span>● at the focus: Sgr A*</span></div>
    <svg id="sky" viewBox="0 0 560 560" role="img" aria-label="S2 position on the sky"></svg>
  </div>
  <div>
    <div class="card"><h2>Line-of-sight velocity (km/s; + = moving away from us)</h2><svg id="rvp" viewBox="0 0 560 190"></svg></div>
    <div class="card" style="margin-top:10px"><h2>Distance from Sgr A* (AU, log scale)</h2><svg id="rp" viewBox="0 0 560 150"></svg></div>
    <div class="card" style="margin-top:10px"><h2>Speed (km/s)</h2><svg id="vp" viewBox="0 0 560 150"></svg></div>
  </div>
</div>
<table class="readout" id="readout"></table>
<p class="note" id="note"></p>
</div>
<div id="tip"></div>
<script>
const D = __DATA__;
const T = D.track, M = D.meta;
const NS = "http://www.w3.org/2000/svg";
const $ = (id) => document.getElementById(id);
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
function lerpIndex(arr, v) { let lo = 0, hi = arr.length - 1;
  if (v <= arr[0]) return [0, 0]; if (v >= arr[hi]) return [hi, 0];
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (arr[mid] <= v) lo = mid; else hi = mid; }
  return [lo, (v - arr[lo]) / (arr[lo + 1] - arr[lo])]; }
function at(key, t) { const [i, f] = lerpIndex(T.t, t); const a = T[key];
  return i + 1 < a.length ? a[i] + f * (a[i + 1] - a[i]) : a[i]; }
function fmtDate(t) { const y = Math.floor(t); const d = new Date(Date.UTC(y, 0, 1) + (t - y) * 365.25 * 864e5);
  return d.toISOString().slice(0, 10); }
const tip = $("tip");
function showTip(ev, html) { tip.innerHTML = html; tip.style.display = "block";
  tip.style.left = Math.min(ev.clientX + 12, innerWidth - 220) + "px"; tip.style.top = (ev.clientY + 12) + "px"; }
function hideTip() { tip.style.display = "none"; }

// ---------- sky ----------
const sky = $("sky"); const SW = 560, SH = 560, PAD = 44;
const VIEWS = { full: { x: [80, -100], y: [-65, 200] }, zoom: { x: [26, -32], y: [-18, 26] } };
let view = "full";
function sx(x) { const v = VIEWS[view]; const w = SW - 2 * PAD, h = SH - 2 * PAD;
  const span = Math.max((v.x[0] - v.x[1]) / w, (v.y[1] - v.y[0]) / h);
  const cx = (v.x[0] + v.x[1]) / 2; return SW / 2 + (cx - x) / span; }
function sy(y) { const v = VIEWS[view]; const w = SW - 2 * PAD, h = SH - 2 * PAD;
  const span = Math.max((v.x[0] - v.x[1]) / w, (v.y[1] - v.y[0]) / h);
  const cy = (v.y[0] + v.y[1]) / 2; return SH / 2 - (y - cy) / span; }
let skyLayers = {};
function drawSky() {
  sky.innerHTML = "";
  const v = VIEWS[view];
  const clip = el("clipPath", { id: "clipsky" }, el("defs", {}, sky));
  el("rect", { x: 0, y: 0, width: SW, height: SH }, clip);
  const g = el("g", { "clip-path": "url(#clipsky)" }, sky);
  const step = view === "full" ? 50 : 10;
  for (let x = Math.ceil(v.x[1] / step) * step; x <= v.x[0]; x += step) {
    el("line", { x1: sx(x), x2: sx(x), y1: 0, y2: SH, stroke: "var(--grid)", "stroke-width": 1 }, g);
    el("text", { x: sx(x), y: SH - 4, "text-anchor": "middle" }, sky).textContent = x; }
  for (let y = Math.ceil(v.y[0] / step) * step; y <= v.y[1]; y += step) {
    el("line", { x1: 0, x2: SW, y1: sy(y), y2: sy(y), stroke: "var(--grid)", "stroke-width": 1 }, g);
    el("text", { x: 4, y: sy(y) - 3 }, sky).textContent = y; }
  el("text", { x: SW / 2, y: SH - 18, "text-anchor": "middle", class: "lbl" }, sky).textContent = "Δ right ascension (mas, East ←)";
  el("text", { x: 14, y: 16, class: "lbl" }, sky).textContent = "Δ declination (mas)";
  // one full orbit ending at the 2018 pericentre
  let d = ""; for (let i = 0; i < T.t.length; i++) { if (T.t[i] < M.t0 - M.period || T.t[i] > M.t0) continue;
    d += (d ? "L" : "M") + sx(T.x[i]).toFixed(1) + "," + sy(T.y[i]).toFixed(1); }
  el("path", { d, fill: "none", stroke: "var(--future)", "stroke-width": 1.2 }, g);
  skyLayers.trail = el("path", { fill: "none", stroke: "var(--orbit)", "stroke-width": 1.6 }, g);
  // scale bar: 500 AU on the whole orbit, 50 AU when zoomed; its West end sits near the right edge
  const au = view === "full" ? 500 : 50, bar = au / M.r0;
  const xw = v.x[1] + 0.08 * (v.x[0] - v.x[1]), xe = xw + bar;
  const byy = sy(v.y[0] + (v.y[1] - v.y[0]) * 0.06);
  el("line", { x1: sx(xe), x2: sx(xw), y1: byy, y2: byy, stroke: "var(--ink)", "stroke-width": 1.5 }, g);
  el("text", { x: (sx(xe) + sx(xw)) / 2, y: byy - 6, "text-anchor": "middle", class: "lbl" }, g).textContent =
    `${au} AU (${bar.toFixed(1)} mas)`;
  // data
  skyLayers.data = [];
  for (const a of D.data.astro) {
    const [t, x, y, ex, ey, team, inst] = a;
    const col = team === "keck" ? "var(--keck)" : "var(--vlt)";
    const grp = el("g", {}, g);
    el("line", { x1: sx(x - ex), x2: sx(x + ex), y1: sy(y), y2: sy(y), stroke: col, "stroke-opacity": 0.45 }, grp);
    el("line", { x1: sx(x), x2: sx(x), y1: sy(y - ey), y2: sy(y + ey), stroke: col, "stroke-opacity": 0.45 }, grp);
    const mk = team === "keck"
      ? el("circle", { cx: sx(x), cy: sy(y), r: 3.6, fill: col, stroke: "var(--panel)", "stroke-width": 1 }, grp)
      : el("rect", { x: sx(x) - 3, y: sy(y) - 3, width: 6, height: 6, fill: col, stroke: "var(--panel)", "stroke-width": 1 }, grp);
    mk.addEventListener("mousemove", (ev) => showTip(ev, `<b>${inst}</b> ${fmtDate(t)}<br>Δα ${x.toFixed(1)} ± ${ex} mas, Δδ ${y.toFixed(1)} ± ${ey} mas`));
    mk.addEventListener("mouseleave", hideTip);
    skyLayers.data.push({ t, team, grp });
  }
  el("circle", { cx: sx(0), cy: sy(0), r: view === "full" ? 5 : 6, fill: "var(--ink)", stroke: "var(--panel)", "stroke-width": 2 }, g);
  el("text", { x: sx(0) + 9, y: sy(0) + 15, class: "lbl" }, g).textContent = "Sgr A*";
  skyLayers.star = el("circle", { r: 7, fill: "none", stroke: "var(--ink)", "stroke-width": 2 }, g);
  skyLayers.dot = el("circle", { r: 3, fill: "var(--ink)" }, g);
}

// ---------- time-series panels ----------
function panel(svg, { key, ylog, ymin, ymax, ticks, data }) {
  const W = 560, H = +svg.getAttribute("viewBox").split(" ")[3], L = 50, R = 10, TP = 8, B = 22;
  const X = (t) => L + (t - 1992) / (2035 - 1992) * (W - L - R);
  const Y = ylog ? (v) => TP + (Math.log10(ymax) - Math.log10(v)) / (Math.log10(ymax) - Math.log10(ymin)) * (H - TP - B)
                 : (v) => TP + (ymax - v) / (ymax - ymin) * (H - TP - B);
  svg.innerHTML = "";
  for (const tk of ticks) { el("line", { x1: L, x2: W - R, y1: Y(tk), y2: Y(tk), stroke: "var(--grid)" }, svg);
    el("text", { x: L - 6, y: Y(tk) + 4, "text-anchor": "end" }, svg).textContent = tk.toLocaleString("en-US"); }
  for (let yr = 1995; yr <= 2035; yr += 5) { el("text", { x: X(yr), y: H - 6, "text-anchor": "middle" }, svg).textContent = yr; }
  let d = ""; for (let i = 0; i < T.t.length; i++) { const v = T[key][i]; if (ylog && v <= 0) continue;
    d += (d ? "L" : "M") + X(T.t[i]).toFixed(1) + "," + Y(Math.min(Math.max(v, ymin), ymax)).toFixed(1); }
  el("path", { d, fill: "none", stroke: "var(--orbit)", "stroke-width": 1.4 }, svg);
  el("line", { x1: X(2019.4), x2: X(2019.4), y1: TP, y2: H - B, stroke: "var(--axis)" }, svg);
  el("text", { x: X(2019.4) + 4, y: TP + 10 }, svg).textContent = "prediction →";
  const pts = [];
  if (data) for (const r of data) { const [t, v, e, team, inst] = r;
    const col = team === "keck" ? "var(--keck)" : "var(--vlt)";
    const grp = el("g", {}, svg);
    el("line", { x1: X(t), x2: X(t), y1: Y(v - e), y2: Y(v + e), stroke: col }, grp);
    const c = el("circle", { cx: X(t), cy: Y(v), r: 3, fill: col, stroke: "var(--panel)", "stroke-width": 1 }, grp);
    c.addEventListener("mousemove", (ev) => showTip(ev, `<b>${inst}</b> ${fmtDate(t)}<br>${v.toLocaleString("en-US")} ± ${e} km/s`));
    c.addEventListener("mouseleave", hideTip);
    pts.push({ t, team, grp }); }
  const cursor = el("line", { y1: TP, y2: H - B, stroke: "var(--ink)", "stroke-width": 1 }, svg);
  const hit = el("rect", { x: L, y: TP, width: W - L - R, height: H - TP - B, fill: "transparent" }, svg);
  hit.addEventListener("mousemove", (ev) => { const b = svg.getBoundingClientRect();
    const xx = (ev.clientX - b.left) / b.width * W; const t = 1992 + (xx - L) / (W - L - R) * (2035 - 1992);
    if (!playing) setTime(Math.min(2035, Math.max(1992, t))); });
  svg.insertBefore(hit, cursor);
  if (data) for (const p of pts) svg.appendChild(p.grp);
  return { cursor, X, pts };
}
const P = {};
function drawPanels() {
  P.rv = panel($("rvp"), { key: "rv", ymin: -2600, ymax: 4600, ticks: [-2000, 0, 2000, 4000], data: D.data.rv });
  P.r = panel($("rp"), { key: "r", ylog: true, ymin: 80, ymax: 3000, ticks: [100, 300, 1000, 3000] });
  P.v = panel($("vp"), { key: "speed", ymin: 0, ymax: 8200, ticks: [0, 2000, 4000, 6000, 8000] });
}

// ---------- state ----------
let now = 2018.38, playing = false, last = null;
const show = { keck: true, vlt: true };
function setTime(t) {
  now = t; $("slider").value = t; $("date").textContent = fmtDate(t);
  const x = at("x", t), y = at("y", t);
  skyLayers.star.setAttribute("cx", sx(x)); skyLayers.star.setAttribute("cy", sy(y));
  skyLayers.dot.setAttribute("cx", sx(x)); skyLayers.dot.setAttribute("cy", sy(y));
  let d = ""; const t0 = Math.max(1992, t - 16);
  for (let i = 0; i < T.t.length; i++) { if (T.t[i] < t0 || T.t[i] > t) continue;
    d += (d ? "L" : "M") + sx(T.x[i]).toFixed(1) + "," + sy(T.y[i]).toFixed(1); }
  d += (d ? "L" : "M") + sx(x).toFixed(1) + "," + sy(y).toFixed(1);
  skyLayers.trail.setAttribute("d", d);
  for (const p of skyLayers.data) p.grp.style.display = (p.t <= t && show[p.team]) ? "" : "none";
  for (const k of ["rv", "r", "v"]) { const c = P[k].cursor; c.setAttribute("x1", P[k].X(t)); c.setAttribute("x2", P[k].X(t)); }
  for (const p of P.rv.pts) p.grp.style.display = show[p.team] ? "" : "none";
  const r = at("r", t), v = at("speed", t);
  const rows = [
    ["Date", fmtDate(t) + (t > 2019.4 ? " (prediction)" : "")],
    ["Distance from Sgr A*", `${Math.round(r).toLocaleString("en-US")} AU = ${(r / M.light_hour_au).toFixed(1)} light-hours = ${Math.round(r / M.rs_au).toLocaleString("en-US")} Schwarzschild radii`],
    ["Speed", `${Math.round(v).toLocaleString("en-US")} km/s = ${(100 * v / M.c_kms).toFixed(2)}% of the speed of light`],
    ["Line-of-sight velocity", `${Math.round(at("rv", t)).toLocaleString("en-US")} km/s (${at("rv", t) > 0 ? "moving away" : "coming towards us"})`],
    ["Relativistic redshift", `${at("redshift", t).toFixed(1)} km/s (time dilation + gravity; Newton: 0)`],
    ["Light-travel delay across the orbit", `${at("delay", t).toFixed(2)} days (${at("delay", t) > 0 ? "behind" : "in front of"} the black hole)`],
    ["Position on the sky", `Δα ${x.toFixed(1)} mas, Δδ ${y.toFixed(1)} mas from Sgr A*`],
  ];
  $("readout").innerHTML = rows.map(([a, b]) => `<tr><td>${a}</td><td>${b}</td></tr>`).join("");
}
function frame(ts) { if (!playing) return; if (last !== null) { const dt = (ts - last) / 1000;
    let t = now + dt * +$("speed").value; if (t > 2035) t = 1992; setTime(t); }
  last = ts; requestAnimationFrame(frame); }
$("play").onclick = () => { playing = !playing; $("play").setAttribute("aria-pressed", playing);
  $("play").textContent = playing ? "❚❚ Pause" : "▶ Play"; last = null; if (playing) requestAnimationFrame(frame); };
$("slider").oninput = (e) => setTime(+e.target.value);
$("zoom").onclick = () => { view = view === "full" ? "zoom" : "full"; $("zoom").setAttribute("aria-pressed", view === "zoom");
  $("zoom").textContent = view === "zoom" ? "Whole orbit" : "Zoom to Sgr A*"; drawSky(); setTime(now); };
$("showKeck").onchange = (e) => { show.keck = e.target.checked; setTime(now); };
$("showVlt").onchange = (e) => { show.vlt = e.target.checked; setTime(now); };
$("note").textContent = `Model: our joint fit of both teams' data with general relativity — black hole mass ${M.gm6.toFixed(3)} ± ${M.gm6_err.toFixed(3)} million suns at ${M.r0.toFixed(3)} ± ${M.r0_err.toFixed(3)} kpc. ` +
  `Positions are drawn relative to the fitted position of the mass, each team's reference frame corrected for its own offset and drift. ` +
  `Data: Do et al. 2019, Science 365, 664 (Keck; velocities also from Subaru, Gemini and the VLT) and Gillessen et al. 2017, ApJ 837, 30 (NTT/VLT). ` +
  `Next pericentre: ${fmtDate(M.t0 + M.period)}.`;
drawSky(); drawPanels(); setTime(now);
</script>
</body>
</html>
"""
