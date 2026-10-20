"""The ring explorer: fit M87's black hole by hand. One self-contained HTML page, no libraries.

Sliders set a lopsided, blurred ring (an m = 1 'mG-ring': diameter, thickness,
lopsidedness, brightest direction, compact flux). The page recomputes, in the
browser, what the EHT would have measured from that ring on the real baselines
of a chosen night in April 2017, and draws it over the real data: visibility
amplitudes against baseline length, and closure phases on three triangles. The
chi^2 of the closure phases is shown as you drag. A second readout turns the
diameter into a mass.

The model's transform is V = F exp(-(pi w rho)^2 / 4 ln 2) [J0(pi d rho) +
2 i A J1(pi d rho) cos(phi_u - eta)], identical to models.MGRing with m = 1 (the
tests check this). The page follows the viewer's light or dark setting.
"""

from __future__ import annotations

import json

import numpy as np

import analysis as A
import uvdata as U
from config import FIGURE_DIR, M87_2017_DAYS, PUBLISHED, SHADOW_FACTOR_SCHWARZSCHILD, STATION_NAMES

VIEWER_PATH = FIGURE_DIR / "ring_explorer.html"
TRIANGLES = [("AA", "LM", "SM"), ("AA", "AZ", "LM"), ("AZ", "LM", "SM")]


def _night(day: str) -> dict:
    from figures import gmst_hours

    prep = U.prepare(A.m87_2017_path(day, "lo"))
    d = U.canonical_orientation(prep.avg)
    amp = U.debiased_amplitude(np.abs(d.vis), d.sigma)
    keep = ~d.intra_site()
    idx_map = -np.ones(len(d), dtype=int)
    idx_map[keep] = np.arange(keep.sum())
    vis = {"u": np.round(d.u[keep] / 1e9, 5).tolist(), "v": np.round(d.v[keep] / 1e9, 5).tolist(),
           "amp": np.round(amp[keep], 5).tolist(), "sig": np.round(d.sigma[keep], 5).tolist(),
           "bl": [f"{STATION_NAMES[a]}-{STATION_NAMES[b]}" for a, b in zip(d.t1[keep], d.t2[keep])]}
    # the minimal closure-phase set (diagonal errors here; the analysis uses the full covariance)
    cp = prep.cphase
    rows = []
    for k in range(len(cp)):
        nz = np.flatnonzero(cp.design[k])
        if np.all(idx_map[nz] >= 0):
            rows.append({"i": idx_map[nz].tolist(), "s": cp.design[k, nz].astype(int).tolist(),
                         "psi": round(float(np.degrees(cp.values[k])), 3), "sig": round(float(np.degrees(cp.sigma[k])), 3)})
    # closure phases on three named triangles, against sidereal time
    tris = []
    for tri in TRIANGLES:
        a, b, c = sorted(tri)
        pts = []
        for s in np.unique(d.scan):
            m = {(x, y): np.flatnonzero((d.scan == s) & (d.t1 == x) & (d.t2 == y)) for x, y in ((a, b), (b, c), (a, c))}
            if all(len(v) for v in m.values()):
                i, j, k = m[(a, b)][0], m[(b, c)][0], m[(a, c)][0]
                psi = np.angle(d.vis[i] * d.vis[j] * np.conj(d.vis[k]))
                sig = np.sqrt(sum((d.sigma[q] / abs(d.vis[q])) ** 2 for q in (i, j, k)))
                pts.append({"t": round(float(gmst_hours(prep.raw.mjd, d.time[i])), 4),
                            "psi": round(float(np.degrees(psi)), 2), "sig": round(float(np.degrees(sig)), 2),
                            "uv": [[round(float(d.u[q] / 1e9), 5), round(float(d.v[q] / 1e9), 5)] for q in (i, j, k)]})
        tris.append({"name": " – ".join(STATION_NAMES[x] for x in tri), "pts": pts})
    return {"label": M87_2017_DAYS[day], "vis": vis, "cp": rows, "tris": tris}


def _presets(res: dict) -> dict:
    out = {}
    for day in M87_2017_DAYS:
        f = res.get("m87_2017_fits", {}).get(f"m87_2017_{day}_lo", {}).get("fiducial")
        if f:
            out[day] = {"d": round(f["d"]["median"], 2), "w": round(f["w"]["median"], 2),
                        "A": round(f["asymmetry"]["median"], 3), "eta": round(f["pa_deg"]["median"], 1)}
    return out


def build(res: dict | None = None) -> str:
    res = res or {}
    data = {"nights": {day: _night(day) for day in M87_2017_DAYS}, "order": list(M87_2017_DAYS), "presets": _presets(res),
            "alpha_grmhd": PUBLISHED["alpha_xsring"], "alpha_shadow": SHADOW_FACTOR_SCHWARZSCHILD,
            "distance_mpc": PUBLISHED["m87_distance_mpc"][0]}
    html = PAGE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    VIEWER_PATH.write_text(html)
    return str(VIEWER_PATH)


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>M87 Ring Explorer</title>
<style>
:root {
  --bg: #fcfcfb; --panel: #ffffff; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --data: #2a78d6; --model: #eb6834; --gr: #1baf7a;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
    --grid: #2c2c2a; --axis: #4a4945; --data: #5a9be6; --model: #f08257; --gr: #34c48f;
  }
}
:root[data-theme="dark"] {
  --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
  --grid: #2c2c2a; --axis: #4a4945; --data: #5a9be6; --model: #f08257; --gr: #34c48f;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
.wrap { max-width: 1200px; margin: 0 auto; padding: 18px 16px 28px; }
h1 { font-size: 20px; margin: 0 0 2px; }
.sub { color: var(--ink2); margin: 0 0 12px; max-width: 920px; }
.controls { display: grid; grid-template-columns: repeat(auto-fill, minmax(360px, 1fr)); gap: 6px 26px; margin: 8px 0 12px; }
.controls label { color: var(--ink2); display: grid; grid-template-columns: 10.5em minmax(0, 1fr) 5.6em; gap: 10px; align-items: center; }
.controls input { width: 100%; }
.controls output { font-variant-numeric: tabular-nums; color: var(--ink); text-align: right; }
.bar { display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: center; margin-bottom: 8px; }
.bar button, .bar select { font: inherit; color: var(--ink); background: var(--panel);
  border: 1px solid var(--axis); border-radius: 6px; padding: 4px 10px; cursor: pointer; }
.grid { display: grid; grid-template-columns: minmax(0, 0.8fr) minmax(0, 1.2fr); gap: 14px; }
@media (max-width: 860px) { .grid { grid-template-columns: minmax(0, 1fr); } .controls label { grid-template-columns: 8em 1fr 4em; } }
.card { background: var(--panel); border: 1px solid var(--grid); border-radius: 8px; padding: 8px 10px; }
.card h2 { font-size: 13px; margin: 2px 0 4px; color: var(--ink2); font-weight: 600; }
canvas { display: block; width: 100%; height: auto; border-radius: 4px; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
svg text { fill: var(--muted); font-size: 11px; }
.readout { width: 100%; border-collapse: collapse; margin-top: 10px; font-variant-numeric: tabular-nums; }
.readout td { padding: 3px 8px; border-bottom: 1px solid var(--grid); }
.readout td:first-child { color: var(--ink2); width: 52%; }
.legend { color: var(--ink2); font-size: 12px; margin: 0 0 4px; }
.sw { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin: 0 4px 0 10px; vertical-align: -1px; }
.note { color: var(--muted); font-size: 12px; margin-top: 14px; max-width: 1000px; }
</style>
</head>
<body>
<div class="wrap">
<h1>Fit the black hole in M87 by hand</h1>
<p class="sub">The Event Horizon Telescope never took a picture: each pair of telescopes measured one Fourier component
of the sky. Shape a ring with the sliders and see, instantly, what the EHT would have recorded from it on the real
baselines of April 2017, drawn over what it did record. A ring about 45 µas across puts the dip in the amplitudes
where the data have it. Lopsidedness and its direction set the closure phases.</p>
<div class="bar">
  <label>Night <select id="night"></select></label>
  <button id="fit">Use our fitted ring</button>
  <button id="reset">Start from a disk-like guess</button>
</div>
<div class="controls">
  <label>Diameter d <input id="d" type="range" min="15" max="80" step="0.1" value="30"><output id="dO"></output></label>
  <label>Thickness (FWHM) <input id="w" type="range" min="0" max="35" step="0.1" value="20"><output id="wO"></output></label>
  <label>Lopsidedness |β₁| <input id="A" type="range" min="0" max="0.5" step="0.005" value="0"><output id="AO"></output></label>
  <label>Brightest side (PA) <input id="eta" type="range" min="0" max="360" step="1" value="180"><output id="etaO"></output></label>
  <label>Compact flux <input id="F" type="range" min="0.2" max="1.2" step="0.01" value="0.6"><output id="FO"></output></label>
</div>
<div class="grid">
  <div>
    <div class="card"><h2>The ring you made (east to the left, 128 µas across)</h2><canvas id="img" width="180" height="180" aria-label="model image" style="image-rendering:auto"></canvas></div>
    <div class="card" style="margin-top:10px"><h2>Its Fourier transform |V(u,v)|, and where the EHT sampled it</h2>
      <canvas id="uv" width="360" height="360" aria-label="Fourier amplitude with baseline tracks"></canvas></div>
  </div>
  <div>
    <div class="card"><h2>Visibility amplitude against baseline length (log scale)</h2>
      <div class="legend"><span class="sw" style="background:var(--data)"></span>EHT data (debiased, ±1σ)<span class="sw" style="background:var(--model)"></span>your ring, at the same (u, v)</div>
      <svg id="amp" viewBox="0 0 620 300"></svg></div>
    <div class="card" style="margin-top:10px"><h2>Closure phases (degrees) against sidereal time: station errors cancel</h2>
      <svg id="cp" viewBox="0 0 620 300"></svg></div>
    <table class="readout" id="readout"></table>
  </div>
</div>
<p class="note">Data: EHT 2017 M87 release 2019-D01-01 (HOPS pipeline, low band, scan-averaged; intra-site baselines
left out). The χ² uses the minimal set of closure phases with independent errors; the analysis behind the README
uses their full covariance and richer models (a second brightness mode, a central floor and two extra Gaussian
components), which is why "our fitted ring" does not reach χ²/N ≈ 1 here. Station gains are not modelled on this
page, so amplitudes from the LMT and SMT can sit off by tens of percent. Mass: θ_g = d / α, M = θ_g D c² / G with
D = 16.8 Mpc.</p>
</div>
<script>
const D = __DATA__;
const NS = "http://www.w3.org/2000/svg";
const $ = (id) => document.getElementById(id);
const RAD = Math.PI / 180 / 3600 / 1e6;  // radians per microarcsecond
function css(v) { return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }

// Bessel functions J0, J1 (Numerical Recipes rational approximations, |error| < 1e-7)
function j0(x) { const ax = Math.abs(x);
  if (ax < 8) { const y = x * x;
    return (57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * -184.9052456))))) /
           (57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y))))); }
  const z = 8 / ax, y = z * z, xx = ax - 0.785398164;
  return Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * (1 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6))))
    - z * Math.sin(xx) * (-0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7))))); }
function j1(x) { const ax = Math.abs(x);
  if (ax < 8) { const y = x * x;
    return x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * -30.16036606))))) /
           (144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))))); }
  const z = 8 / ax, y = z * z, xx = ax - 2.356194491;
  const a = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * (1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * -0.240337019e-6))))
    - z * Math.sin(xx) * (0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)))));
  return x < 0 ? -a : a; }

let P = {};
function readP() { for (const k of ["d", "w", "A", "eta", "F"]) P[k] = parseFloat($(k).value);
  $("dO").textContent = P.d.toFixed(1) + " µas"; $("wO").textContent = P.w.toFixed(1) + " µas";
  $("AO").textContent = P.A.toFixed(3); $("etaO").textContent = P.eta.toFixed(0) + "°"; $("FO").textContent = P.F.toFixed(2) + " Jy"; }

// model visibility at (u, v) in billions of wavelengths -> [re, im]
function vis(u, v) { const rho = Math.hypot(u, v) * 1e9, phi = Math.atan2(u, v);
  const z = Math.PI * P.d * RAD * rho, blur = Math.exp(-Math.pow(Math.PI * P.w * RAD * rho, 2) / (4 * Math.LN2));
  return [P.F * blur * j0(z), P.F * blur * 2 * P.A * j1(z) * Math.cos(phi - P.eta * Math.PI / 180)]; }

function drawImage() { const c = $("img"), g = c.getContext("2d"), n = c.width, fov = 128;
  const im = g.createImageData(n, n), vals = new Float32Array(n * n); let mx = 0;
  // exact image: the thin lopsided ring, sampled at 96 points, convolved with the Gaussian
  const sig = Math.max(P.w, 1.5) / 2.3548, r0 = P.d / 2, eta = P.eta * Math.PI / 180, K = 96;
  const rx = [], ry = [], rw = [];
  for (let k = 0; k < K; k++) { const th = (k + 0.5) * 2 * Math.PI / K;
    rx.push(r0 * Math.sin(th)); ry.push(r0 * Math.cos(th)); rw.push(1 + 2 * P.A * Math.cos(th - eta)); }
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const x = fov / 2 - (j + 0.5) * fov / n, y = fov / 2 - (i + 0.5) * fov / n;  // x east (left), y north (up)
    let val = 0;
    for (let k = 0; k < K; k++) { const dx = x - rx[k], dy = y - ry[k]; val += rw[k] * Math.exp(-(dx * dx + dy * dy) / (2 * sig * sig)); }
    vals[i * n + j] = val; if (val > mx) mx = val; }
  for (let k = 0; k < n * n; k++) { const t = Math.max(0, vals[k] / mx);
    // 'afmhot': black -> red -> yellow -> white
    im.data[4 * k] = 255 * Math.min(1, 2 * t); im.data[4 * k + 1] = 255 * Math.min(1, Math.max(0, 2 * t - 0.5));
    im.data[4 * k + 2] = 255 * Math.min(1, Math.max(0, 2 * t - 1)); im.data[4 * k + 3] = 255; }
  g.putImageData(im, 0, 0);
  g.strokeStyle = "#ffffff"; g.lineWidth = 1.5; g.beginPath(); g.moveTo(n * 0.06, n * 0.94); g.lineTo(n * 0.06 + n * 20 / fov, n * 0.94); g.stroke();
  g.fillStyle = "#ffffff"; g.font = "8px system-ui"; g.fillText("20 µas", n * 0.06, n * 0.90);
  g.fillText("N ↑  E ←", n * 0.74, n * 0.08); }

function drawUV() { const c = $("uv"), g = c.getContext("2d"), n = c.width, L = 9.5;
  const im = g.createImageData(n, n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const u = L - (j + 0.5) * 2 * L / n, v = L - (i + 0.5) * 2 * L / n, V = vis(u, v);
    const a = Math.hypot(V[0], V[1]) / P.F, t = Math.max(0, (Math.log10(Math.max(a, 1e-3)) + 3) / 3);
    const k = 4 * (i * n + j); im.data[k] = im.data[k + 1] = im.data[k + 2] = 255 * (0.12 + 0.88 * t); im.data[k + 3] = 255; }
  g.putImageData(im, 0, 0);
  const N = D.nights[$("night").value]; g.fillStyle = css("--data");
  for (let k = 0; k < N.vis.u.length; k++) for (const s of [1, -1]) {
    const px = (L - s * N.vis.u[k]) / (2 * L) * n, py = (L - s * N.vis.v[k]) / (2 * L) * n;
    g.beginPath(); g.arc(px, py, 2, 0, 2 * Math.PI); g.fill(); }
  g.fillStyle = "#ffffff"; g.font = "12px system-ui"; g.fillText("u, v within ±9.5 billion λ (grey: log |V|)", 8, 16); }

function axes(svg, x0, x1, y0, y1, xl, yl, logy) { svg.innerHTML = ""; const W = 620, H = 300, m = {l: 54, r: 10, t: 8, b: 34};
  const sx = (x) => m.l + (x - x0) / (x1 - x0) * (W - m.l - m.r);
  const sy = (y) => { const t = logy ? (Math.log10(y) - Math.log10(y0)) / (Math.log10(y1) - Math.log10(y0)) : (y - y0) / (y1 - y0);
    return H - m.b - t * (H - m.t - m.b); };
  const grid = css("--grid"), axis = css("--axis");
  const yt = logy ? [0.001, 0.01, 0.1, 1] : [-180, -90, 0, 90, 180];
  for (const t of yt) { if (t < y0 || t > y1) continue; el("line", {x1: m.l, x2: W - m.r, y1: sy(t), y2: sy(t), stroke: grid}, svg);
    el("text", {x: m.l - 6, y: sy(t) + 4, "text-anchor": "end"}, svg).textContent = t; }
  for (let t = Math.ceil(x0); t <= x1; t += (x1 - x0 > 12 ? 2 : 1)) { el("line", {x1: sx(t), x2: sx(t), y1: m.t, y2: H - m.b, stroke: grid}, svg);
    el("text", {x: sx(t), y: H - m.b + 15, "text-anchor": "middle"}, svg).textContent = t; }
  el("line", {x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: axis}, svg);
  el("text", {x: (W + m.l) / 2, y: H - 4, "text-anchor": "middle"}, svg).textContent = xl;
  el("text", {x: 12, y: H / 2, transform: `rotate(-90 12 ${H / 2})`, "text-anchor": "middle"}, svg).textContent = yl;
  return [sx, sy]; }

function drawAmp() { const N = D.nights[$("night").value], svg = $("amp");
  const [sx, sy] = axes(svg, 0, 9, 0.003, 2, "baseline length (billion wavelengths)", "amplitude (Jy)", true);
  const dc = css("--data"), mc = css("--model");
  for (let k = 0; k < N.vis.u.length; k++) { const r = Math.hypot(N.vis.u[k], N.vis.v[k]), a = N.vis.amp[k], s = N.vis.sig[k];
    if (a > 0) { el("line", {x1: sx(r), x2: sx(r), y1: sy(Math.max(a - s, 0.003)), y2: sy(Math.min(a + s, 2)), stroke: dc, "stroke-opacity": 0.5}, svg);
      el("circle", {cx: sx(r), cy: sy(a), r: 2.2, fill: dc}, svg); }
    const V = vis(N.vis.u[k], N.vis.v[k]), am = Math.max(Math.hypot(V[0], V[1]), 0.0031);
    el("circle", {cx: sx(r), cy: sy(am), r: 2.6, fill: "none", stroke: mc, "stroke-width": 1.2}, svg); }
  const nul = 2.404826 / (Math.PI * P.d * RAD) / 1e9;
  if (nul < 9) { el("line", {x1: sx(nul), x2: sx(nul), y1: 8, y2: 266, stroke: css("--gr"), "stroke-dasharray": "4 3"}, svg);
    el("text", {x: sx(nul) + 4, y: 20}, svg).textContent = `thin-ring null ${nul.toFixed(2)} Gλ`; } }

function closure(uv) { let ph = 0; const sg = [1, 1, -1];
  for (let q = 0; q < 3; q++) { const V = vis(uv[q][0], uv[q][1]); ph += sg[q] * Math.atan2(V[1], V[0]); }
  ph = Math.atan2(Math.sin(ph), Math.cos(ph)); return ph * 180 / Math.PI; }

function drawCP() { const N = D.nights[$("night").value], svg = $("cp");
  let t0 = 1e9, t1 = -1e9; for (const tr of N.tris) for (const p of tr.pts) { t0 = Math.min(t0, p.t); t1 = Math.max(t1, p.t); }
  const [sx, sy] = axes(svg, Math.floor(t0) - 0.2, Math.ceil(t1) + 0.2, -180, 180, "Greenwich sidereal time (hours)", "closure phase (°)", false);
  const cols = [css("--data"), css("--gr"), css("--ink2")];
  N.tris.forEach((tr, it) => { const c = cols[it];
    for (const p of tr.pts) { el("line", {x1: sx(p.t), x2: sx(p.t), y1: sy(Math.max(-180, p.psi - p.sig)), y2: sy(Math.min(180, p.psi + p.sig)), stroke: c, "stroke-opacity": 0.6}, svg);
      el("circle", {cx: sx(p.t), cy: sy(p.psi), r: 2.6, fill: c}, svg);
      el("rect", {x: sx(p.t) - 2.8, y: sy(closure(p.uv)) - 2.8, width: 5.6, height: 5.6, fill: "none", stroke: css("--model"), "stroke-width": 1.2}, svg); }
    el("circle", {cx: 70 + 175 * it, cy: 18, r: 4, fill: c}, svg);
    el("text", {x: 78 + 175 * it, y: 22}, svg).textContent = tr.name; });
  el("text", {x: 600, y: 22, "text-anchor": "end"}, svg).textContent = "□ your ring"; }

function readout() { const N = D.nights[$("night").value]; let chi = 0;
  const ph = N.vis.u.map((u, k) => { const V = vis(u, N.vis.v[k]); return Math.atan2(V[1], V[0]); });
  for (const r of N.cp) { let m = 0; r.i.forEach((i, q) => m += r.s[q] * ph[i]);
    let d = r.psi * Math.PI / 180 - m; d = Math.atan2(Math.sin(d), Math.cos(d)); chi += Math.pow(d / (r.sig * Math.PI / 180), 2); }
  const G = 6.6743e-11, c = 2.99792458e8, Msun = 1.98841e30, Mpc = 3.0856776e22;
  const mass = (alpha) => (P.d / alpha) * RAD * D.distance_mpc * Mpc * c * c / G / Msun / 1e9;
  const rows = [["closure-phase χ² / number of closure phases", `${chi.toFixed(0)} / ${N.cp.length} = ${(chi / N.cp.length).toFixed(2)}`],
    ["first null of a thin ring this size", `${(2.404826 / (Math.PI * P.d * RAD) / 1e9).toFixed(2)} billion λ`],
    [`mass if the ring is the shadow's edge (d = ${D.alpha_shadow.toFixed(2)} GM/c²D)`, `${mass(D.alpha_shadow).toFixed(2)} billion suns`],
    [`mass with the EHT's GRMHD calibration (d = ${D.alpha_grmhd} GM/c²D)`, `${mass(D.alpha_grmhd).toFixed(2)} billion suns`]];
  $("readout").innerHTML = rows.map(r => `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join(""); }

function update() { readP(); drawImage(); drawUV(); drawAmp(); drawCP(); readout(); }
function setP(p) { for (const k in p) $(k).value = p[k]; update(); }
for (const day of D.order) { const o = document.createElement("option"); o.value = day; o.textContent = D.nights[day].label + ", 2017"; $("night").appendChild(o); }
$("night").value = D.order[D.order.length - 1];  // April 11, the EHT's showcase night
for (const k of ["d", "w", "A", "eta", "F"]) $(k).addEventListener("input", update);
$("night").addEventListener("change", update);
$("fit").addEventListener("click", () => { const p = D.presets[$("night").value]; if (p) setP({d: p.d, w: p.w, A: p.A, eta: p.eta}); });
$("reset").addEventListener("click", () => setP({d: 30, w: 20, A: 0, eta: 180, F: 0.6}));
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", update);
update();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    from config import RESULTS_PATH

    r = json.loads(RESULTS_PATH.read_text()) if RESULTS_PATH.exists() else {}
    print(build(r))
