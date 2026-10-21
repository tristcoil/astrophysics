"""A self-contained, offline explorer: tilt a disc until it matches a galaxy.

    figures/velocity_field_explorer.html

Pick one of the four THINGS galaxies. The measured velocity field sits next
to a model disc whose inclination, orientation and systemic velocity you set
with sliders; the residual map and its rms update as you move them, so you
can try to beat the tilted-ring fit by hand (and see how a wrong inclination
leaves a four-leaf pattern and a wrong position angle a two-leaf one). Click
anywhere on the galaxy to see the 21-cm spectrum there, with the fitted
profiles and the velocity each recipe would report.

Plain HTML, canvas and SVG; no libraries; data embedded (~3 MB), so it works
offline and as a single shared file. Follows the reader's light/dark setting,
and carries the rotation-curve table for anyone who cannot use the maps.
"""

from __future__ import annotations

import base64
import json

import numpy as np

import plotstyle as P
import things
from config import FIGURE_DIR, GALAXIES, galaxy_by_name

THEMES = {
    "light": {"surface": P.SURFACE, "ink": P.INK, "ink2": P.INK_SECONDARY, "muted": P.MUTED, "grid": P.GRID,
              "axis": P.AXIS, "app": P.BLUE, "rec": P.RED, "model": P.ORANGE, "pub": P.VIOLET,
              "panel": "#f3f2ee"},
    "dark": {"surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a",
             "axis": "#4a4a46", "app": "#3987e5", "rec": "#e66767", "model": "#d95926", "pub": "#9085e9",
             "panel": "#242423"},
}

#: Diverging stops of plotstyle.velocity_cmap, light and dark midpoints.
CMAP_LIGHT = ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#f0efec", "#f9cfca", "#ee837e", "#c8302f", "#7a1717"]
CMAP_DARK = ["#5f9df0", "#3987e5", "#256abf", "#18406f", "#383835", "#6e2523", "#b73a37", "#e66767", "#f4a3a0"]

MAX_MAP_PIXELS = 280  # long side of the embedded velocity maps
MAX_SPECTRA = 4000  # per galaxy


def _b64(arr: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode("ascii")


def _galaxy_payload(gr) -> dict:
    import analysis as A

    g = galaxy_by_name(gr.name)
    with np.load(A.DERIVED_DIR / f"{gr.key}_maps.npz") as f:
        vf = f["velocity"]
    cube = things.load_cube(g)
    frame = A.frame_of(cube)
    tr = gr.primary
    rows, cols = np.nonzero(np.isfinite(vf))
    r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
    step = int(np.ceil(max(r1 - r0, c1 - c0) / MAX_MAP_PIXELS))
    sub = vf[r0:r1:step, c0:c1:step]
    enc = np.where(np.isfinite(sub), np.round((sub - tr.vsys) * 10), -32768).astype("<i2")
    # arcsec offsets (west, north) of the first kept pixel, relative to the dynamical centre
    x_first = (c0 - frame.cref) * frame.pixel - tr.centre[0]
    y_first = (r0 - frame.rref) * frame.pixel - tr.centre[1]

    maps = things.velocity_fields(g)
    ok = np.isfinite(vf)
    sp_step = max(1, int(np.ceil(np.sqrt(ok.sum() / MAX_SPECTRA))))
    grid = np.zeros_like(ok)
    grid[::sp_step, ::sp_step] = True
    sr, sc = np.nonzero(grid & ok)
    spectra = np.stack([np.asarray(cube.data[:, r, c], dtype=float) for r, c in zip(sr, sc)])
    smax = np.maximum(np.nanmax(np.abs(spectra), axis=1), 1e-9)
    q = np.clip(np.round(spectra / smax[:, None] * 127), -127, 127).astype(np.int8)
    params = np.column_stack([(sc - frame.cref) * frame.pixel - tr.centre[0],
                              (sr - frame.rref) * frame.pixel - tr.centre[1], smax * 1e3] +
                             [maps[k][sr, sc] for k in ("gh_amp", "gh_v0", "gh_sigma", "gh_h3", "gh_vpeak",
                                                        "gauss_amp", "gauss_v", "gauss_sigma", "peak_v",
                                                        "mom1_things")])
    params[:, [3, 8]] *= 1e3  # amplitudes to mJy/beam
    db = gr.deblok
    return {
        "key": gr.key, "name": gr.name, "distance": gr.distance_mpc, "kpcPerArcsec": gr.kpc_per_arcsec,
        "vsys": round(tr.vsys, 2), "pa": round(float(np.nanmedian(tr.pa)), 1),
        "incl": round(float(np.nanmedian(tr.incl)), 1), "beam": [round(b, 1) for b in gr.beam],
        "map": {"nx": int(sub.shape[1]), "ny": int(sub.shape[0]), "step": float(step * frame.pixel),
                "x0": float(x_first), "y0": float(y_first), "data": _b64(enc)},
        "curve": {"r": np.round(tr.radius, 2).tolist(), "v": np.round(tr.vrot, 2).tolist(),
                  "err": np.round(tr.vrot_err, 2).tolist(), "app": np.round(tr.vrot_app, 2).tolist(),
                  "rec": np.round(tr.vrot_rec, 2).tolist(), "pa": np.round(tr.pa, 2).tolist(),
                  "incl": np.round(tr.incl, 2).tolist()},
        "deblok": {"r": np.round(db.r_kpc / gr.kpc_per_arcsec, 2).tolist(), "v": np.round(db.v, 2).tolist()},
        "spec": {"v": np.round(cube.velocity, 2).tolist(), "n": int(len(sr)), "data": _b64(q.ravel()),
                 "params": np.round(np.nan_to_num(params, nan=-9999), 2).ravel().tolist()},
    }


def build(grs) -> str:
    payload = [_galaxy_payload(grs[g.key]) for g in GALAXIES if g.key in grs]
    L, D = THEMES["light"], THEMES["dark"]
    html = TEMPLATE
    for name, theme in (("L", L), ("D", D)):
        for k, v in theme.items():
            html = html.replace(f"%{name}_{k}%", v)
    html = html.replace("%CMAP_LIGHT%", json.dumps(CMAP_LIGHT)).replace("%CMAP_DARK%", json.dumps(CMAP_DARK))
    html = html.replace("%DATA%", json.dumps(payload, separators=(",", ":")))
    out = FIGURE_DIR / "velocity_field_explorer.html"
    out.write_text(html)
    return str(out)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Velocity Field Explorer</title>
<style>
  :root { color-scheme: light; --surface: %L_surface%; --ink: %L_ink%; --ink2: %L_ink2%; --muted: %L_muted%;
          --grid: %L_grid%; --axis: %L_axis%; --app: %L_app%; --rec: %L_rec%; --model: %L_model%;
          --pub: %L_pub%; --panel: %L_panel%; }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) { color-scheme: dark; --surface: %D_surface%; --ink: %D_ink%; --ink2: %D_ink2%;
          --muted: %D_muted%; --grid: %D_grid%; --axis: %D_axis%; --app: %D_app%; --rec: %D_rec%;
          --model: %D_model%; --pub: %D_pub%; --panel: %D_panel%; }
  }
  :root[data-theme="dark"] { color-scheme: dark; --surface: %D_surface%; --ink: %D_ink%; --ink2: %D_ink2%;
          --muted: %D_muted%; --grid: %D_grid%; --axis: %D_axis%; --app: %D_app%; --rec: %D_rec%;
          --model: %D_model%; --pub: %D_pub%; --panel: %D_panel%; }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--surface); color: var(--ink);
         font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
  main { max-width: 1240px; margin: 0 auto; padding: 20px 16px 40px; }
  h1 { font-size: 22px; margin: 0 0 4px; }
  p.lede { color: var(--ink2); margin: 0 0 14px; max-width: 900px; }
  .controls { display: flex; flex-wrap: wrap; gap: 10px 22px; align-items: center; padding: 12px 14px;
              background: var(--panel); border-radius: 8px; margin-bottom: 14px; }
  .controls label { display: flex; flex-direction: column; gap: 2px; font-size: 12px; color: var(--ink2); }
  .controls input[type=range] { width: 170px; }
  .controls .val { color: var(--ink); font-variant-numeric: tabular-nums; }
  select, button { font: inherit; color: var(--ink); background: var(--surface); border: 1px solid var(--axis);
                   border-radius: 6px; padding: 4px 8px; }
  button { cursor: pointer; }
  .score { margin-left: auto; text-align: right; }
  .score b { font-size: 20px; font-variant-numeric: tabular-nums; }
  .maps { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px; }
  .card { background: var(--panel); border-radius: 8px; padding: 10px; position: relative; min-width: 0; }
  .card h2 { font-size: 13px; margin: 0 0 6px; }
  canvas { width: 100%; height: auto; display: block; image-rendering: pixelated; cursor: crosshair; }
  .colorbar { display: flex; align-items: center; gap: 8px; font-size: 11px; color: var(--ink2); margin-top: 6px; }
  .colorbar .ramp { flex: 1; height: 10px; border-radius: 3px; }
  .lower { display: grid; grid-template-columns: minmax(0, 1.2fr) minmax(0, 1fr); gap: 12px; margin-top: 12px; }
  svg { width: 100%; height: auto; display: block; overflow: visible; }
  svg text { fill: var(--ink2); font-size: 11px; }
  .tip { position: fixed; pointer-events: none; background: var(--surface); color: var(--ink);
         border: 1px solid var(--axis); border-radius: 6px; padding: 6px 8px; font-size: 12px; display: none;
         box-shadow: 0 2px 8px rgba(0,0,0,.15); font-variant-numeric: tabular-nums; z-index: 5; }
  .legend { display: flex; flex-wrap: wrap; gap: 4px 14px; font-size: 11px; color: var(--ink2); margin-top: 4px; }
  .legend span::before { content: ""; display: inline-block; width: 14px; height: 3px; margin-right: 5px;
                         vertical-align: middle; background: var(--c); }
  details { margin-top: 14px; }
  table { border-collapse: collapse; font-size: 12px; font-variant-numeric: tabular-nums; }
  td, th { padding: 2px 10px; text-align: right; border-bottom: 1px solid var(--grid); }
  .note { color: var(--ink2); font-size: 12px; }
  @media (max-width: 860px) { .maps, .lower { grid-template-columns: minmax(0, 1fr); } .score { margin-left: 0; } }
</style>
</head>
<body>
<main>
<h1>Tilt a disc until it matches a galaxy</h1>
<p class="lede">The left map is the measured 21-cm velocity field of a galaxy (THINGS, VLA): blue gas comes toward us,
red moves away. The middle map is a flat disc in circular rotation, seen at the inclination and orientation you
choose. The right map is the difference. Our tilted-ring fit is the starting point; nudge the sliders and watch the
residuals: a wrong inclination leaves a four-leaf pattern, a wrong position angle a two-leaf one. Click the galaxy to
see its hydrogen spectrum there.</p>
<div class="controls">
  <label>Galaxy <select id="gal"></select></label>
  <label>inclination <span class="val" id="inclv"></span><input type="range" id="incl" min="20" max="88" step="0.5"></label>
  <label>position angle <span class="val" id="pav"></span><input type="range" id="pa" min="0" max="360" step="0.5"></label>
  <label>systemic velocity offset <span class="val" id="dvv"></span><input type="range" id="dv" min="-30" max="30" step="0.5"></label>
  <label>rotation curve <select id="curve"><option value="measured">our measured curve</option><option value="flat">flat curve (sliders)</option></select></label>
  <label>flat speed <span class="val" id="vflatv"></span><input type="range" id="vflat" min="10" max="320" step="1"></label>
  <label>turnover radius <span class="val" id="rturnv"></span><input type="range" id="rturn" min="5" max="400" step="1"></label>
  <button id="reset">back to our fit</button>
  <div class="score"><div class="note">rms residual</div><b id="rms">–</b> <span class="note">km/s</span><div class="note" id="best"></div></div>
</div>
<div class="maps">
  <div class="card"><h2>Measured velocity field</h2><canvas id="obs"></canvas></div>
  <div class="card"><h2>Your model disc</h2><canvas id="mod"></canvas></div>
  <div class="card"><h2>Measured minus model</h2><canvas id="res"></canvas></div>
</div>
<div class="colorbar"><span id="cblo"></span><div class="ramp" id="ramp"></div><span id="cbhi"></span>
  <span>line-of-sight velocity relative to systemic (km/s); residual panel: ±25 km/s</span></div>
<div class="lower">
  <div class="card"><h2>Rotation curve</h2><svg id="rc" viewBox="0 0 560 300"></svg>
    <div class="legend"><span style="--c: var(--ink)">measured, both halves</span><span style="--c: var(--app)">approaching half</span>
    <span style="--c: var(--rec)">receding half</span><span style="--c: var(--pub)">de Blok et al. 2008</span><span style="--c: var(--model)">your flat curve</span></div></div>
  <div class="card"><h2 id="spech">21-cm spectrum: click the galaxy</h2><svg id="sp" viewBox="0 0 460 300"></svg>
    <div class="note" id="specnote">Grey: the measured spectrum. Black: Gauss-Hermite fit. Blue: single Gaussian. Vertical lines: the velocity each recipe reports.</div></div>
</div>
<details><summary>Rotation curve as a table</summary><div id="tbl"></div></details>
<p class="note">Data: THINGS (Walter et al. 2008), natural-weighted cubes; rotation curves from our own Gauss-Hermite
velocity fields and tilted-ring fits (project J1 of the Astrophysics Project Atlas, AI-generated). The model disc here
has one inclination and position angle; the real fit lets them vary ring by ring.</p>
</main>
<div class="tip" id="tip"></div>
<script>
const DATA = %DATA%;
const CMAPS = {light: %CMAP_LIGHT%, dark: %CMAP_DARK%};
const $ = id => document.getElementById(id);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const dark = () => document.documentElement.dataset.theme === "dark" ||
  (document.documentElement.dataset.theme !== "light" && matchMedia("(prefers-color-scheme: dark)").matches);

function hex2rgb(h) { return [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16)); }
function lut() {
  const stops = (dark() ? CMAPS.dark : CMAPS.light).map(hex2rgb), out = [];
  for (let i = 0; i < 256; i++) {
    const t = i / 255 * (stops.length - 1), k = Math.min(Math.floor(t), stops.length - 2), f = t - k;
    out.push(stops[k].map((c, j) => Math.round(c + (stops[k + 1][j] - c) * f)));
  }
  return out;
}
function decode(b64, Type) {
  const bin = atob(b64), buf = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
  return new Type(buf.buffer);
}

let G, V, VM, VR, LUT, span, state = {};
DATA.forEach((g, i) => { const o = document.createElement("option"); o.value = i; o.textContent = g.name; $("gal").appendChild(o); });

function load(i) {
  G = DATA[i];
  const raw = decode(G.map.data, Int16Array);
  V = new Float32Array(raw.length);
  for (let k = 0; k < raw.length; k++) V[k] = raw[k] === -32768 ? NaN : raw[k] / 10;
  VM = new Float32Array(V.length); VR = new Float32Array(V.length);
  let m = 0; for (const v of V) if (v === v) m = Math.max(m, Math.abs(v));
  span = Math.ceil(m / 10) * 10;
  G.specData = decode(G.spec.data, Int8Array);
  for (const c of ["obs", "mod", "res"]) { $(c).width = G.map.nx; $(c).height = G.map.ny; }
  reset();
}
function reset() {
  const vmax = Math.max(...G.curve.v.filter(v => v === v));
  Object.assign(state, {incl: G.incl, pa: G.pa, dv: 0, vflat: Math.round(vmax), rturn: Math.round(G.curve.r[G.curve.r.length - 1] / 6)});
  for (const k of ["incl", "pa", "dv", "vflat", "rturn"]) $(k).value = state[k];
  $("curve").value = "measured";
  update(true);
}
function vcurve(r) {
  if ($("curve").value === "flat") return state.vflat * (1 - Math.exp(-r / state.rturn));
  const R = G.curve.r, Vc = G.curve.v;
  if (r <= R[0]) return Vc[0] * r / R[0];
  for (let k = 1; k < R.length; k++) if (r <= R[k]) {
    const a = Vc[k - 1], b = Vc[k];
    if (a !== a || b !== b) return NaN;
    return a + (b - a) * (r - R[k - 1]) / (R[k] - R[k - 1]);
  }
  return NaN;
}
function update(first) {
  $("inclv").textContent = state.incl.toFixed(1) + "°"; $("pav").textContent = state.pa.toFixed(1) + "°";
  $("dvv").textContent = (state.dv >= 0 ? "+" : "") + state.dv.toFixed(1) + " km/s";
  $("vflatv").textContent = state.vflat + " km/s"; $("rturnv").textContent = state.rturn + '"';
  const flat = $("curve").value === "flat";
  $("vflat").disabled = $("rturn").disabled = !flat;
  const pa = state.pa * Math.PI / 180, inc = state.incl * Math.PI / 180, si = Math.sin(inc), ci = Math.cos(inc);
  const {nx, ny, step, x0, y0} = G.map;
  let s2 = 0, n = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const k = j * nx + i, v = V[k];
    if (v !== v) { VM[k] = VR[k] = NaN; continue; }
    const x = x0 + i * step, y = y0 + j * step;
    const along = -x * Math.sin(pa) + y * Math.cos(pa), across = (-x * Math.cos(pa) - y * Math.sin(pa)) / ci;
    const r = Math.hypot(along, across), cosT = r > 0 ? along / r : 0;
    const m = state.dv + vcurve(r) * si * cosT;
    VM[k] = m; VR[k] = v - m;
    if (m === m) { s2 += (v - m) * (v - m); n++; }
  }
  const rms = Math.sqrt(s2 / Math.max(n, 1));
  $("rms").textContent = rms.toFixed(2);
  if (first) state.best = rms;
  $("best").textContent = `starting point, our fit's median tilt: ${state.best.toFixed(2)}`;
  draw("obs", V, span); draw("mod", VM, span); draw("res", VR, 25);
  drawCurve();
}
function draw(id, arr, lim) {
  const cv = $(id), ctx = cv.getContext("2d"), {nx, ny} = G.map, img = ctx.createImageData(nx, ny);
  const surf = hex2rgb(css("--surface"));
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    // canvas row 0 is the top, so rows flip to put north up; columns already run west to the right (east left)
    const k = j * nx + i, p = ((ny - 1 - j) * nx + i) * 4, v = arr[k];
    const c = v === v ? LUT[Math.max(0, Math.min(255, Math.round((v / lim + 1) * 127.5)))] : surf;
    img.data[p] = c[0]; img.data[p + 1] = c[1]; img.data[p + 2] = c[2]; img.data[p + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  $("cblo").textContent = -span; $("cbhi").textContent = "+" + span;
}
function svgEl(tag, attrs, parent) {
  const e = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}
function axes(svg, W, H, m, xr, yr, xlab, ylab) {
  svg.innerHTML = "";
  const X = x => m.l + (x - xr[0]) / (xr[1] - xr[0]) * (W - m.l - m.r);
  const Y = y => H - m.b - (y - yr[0]) / (yr[1] - yr[0]) * (H - m.t - m.b);
  const ticks = (a, b) => { const s = Math.pow(10, Math.floor(Math.log10((b - a) / 4))); const st = [1, 2, 5, 10].map(f => f * s).find(f => (b - a) / f <= 6); const out = []; for (let t = Math.ceil(a / st) * st; t <= b + 1e-9; t += st) out.push(+t.toFixed(6)); return out; };
  for (const t of ticks(...yr)) { svgEl("line", {x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), stroke: css("--grid")}, svg); svgEl("text", {x: m.l - 6, y: Y(t) + 4, "text-anchor": "end"}, svg).textContent = t; }
  for (const t of ticks(...xr)) { svgEl("text", {x: X(t), y: H - m.b + 16, "text-anchor": "middle"}, svg).textContent = t; }
  svgEl("line", {x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: css("--axis")}, svg);
  svgEl("text", {x: (W + m.l - m.r) / 2, y: H - 4, "text-anchor": "middle"}, svg).textContent = xlab;
  svgEl("text", {x: 12, y: (H - m.b + m.t) / 2, transform: `rotate(-90 12 ${(H - m.b + m.t) / 2})`, "text-anchor": "middle"}, svg).textContent = ylab;
  return {X, Y};
}
function path(svg, xs, ys, X, Y, colour, width) {
  let d = "", pen = false;
  xs.forEach((x, k) => { const y = ys[k]; if (y === y && y !== null) { d += (pen ? "L" : "M") + X(x).toFixed(1) + " " + Y(y).toFixed(1); pen = true; } else pen = false; });
  return svgEl("path", {d, fill: "none", stroke: colour, "stroke-width": width}, svg);
}
function drawCurve() {
  const svg = $("rc"), kpa = G.kpcPerArcsec, c = G.curve;
  const rk = c.r.map(r => r * kpa), top = Math.max(...c.v.filter(v => v === v).map((v, k) => v + c.err[k])) * 1.15;
  const {X, Y} = axes(svg, 560, 300, {l: 48, r: 12, t: 10, b: 36}, [0, Math.max(...rk) * 1.03], [0, top],
                      `radius (kpc, at ${G.distance} Mpc)`, "rotation speed (km/s)");
  path(svg, G.deblok.r.map(r => r * kpa), G.deblok.v, X, Y, css("--pub"), 1.2);
  path(svg, rk, c.app, X, Y, css("--app"), 1.2); path(svg, rk, c.rec, X, Y, css("--rec"), 1.2);
  rk.forEach((r, k) => { if (c.v[k] === c.v[k]) { svgEl("line", {x1: X(r), x2: X(r), y1: Y(c.v[k] - c.err[k]), y2: Y(c.v[k] + c.err[k]), stroke: css("--muted")}, svg); svgEl("circle", {cx: X(r), cy: Y(c.v[k]), r: 2.2, fill: css("--ink")}, svg); } });
  if ($("curve").value === "flat") { const xs = Array.from({length: 120}, (_, k) => k / 119 * c.r[c.r.length - 1]); path(svg, xs.map(x => x * kpa), xs.map(x => state.vflat * (1 - Math.exp(-x / state.rturn))), X, Y, css("--model"), 2); }
}
function gh(v, a, v0, s, h3) { const w = (v - v0) / s; return a * Math.exp(-0.5 * w * w) * (1 + h3 * (2 * Math.SQRT2 * w ** 3 - 3 * Math.SQRT2 * w) / Math.sqrt(6)); }
function showSpectrum(x, y) {
  const P = G.spec.params, n = G.spec.n, w = 13; let best = -1, bd = 1e30;
  for (let k = 0; k < n; k++) { const d = (P[k * w] - x) ** 2 + (P[k * w + 1] - y) ** 2; if (d < bd) { bd = d; best = k; } }
  if (best < 0 || Math.sqrt(bd) > 3 * G.map.step + 10) return;
  const p = P.slice(best * w, best * w + w), vel = G.spec.v, m = vel.length;
  const sp = Array.from({length: m}, (_, c) => G.specData[best * m + c] / 127 * p[2]);
  const vp = p[7], sig = p[5];
  const lo = vp - 6 * Math.max(sig, 5), hi = vp + 6 * Math.max(sig, 5);
  const top = Math.max(...sp) * 1.15, bot = Math.min(0, Math.min(...sp)) * 1.1;
  const svg = $("sp"), {X, Y} = axes(svg, 460, 300, {l: 48, r: 12, t: 10, b: 36}, [lo, hi], [bot, top], "heliocentric velocity (km/s)", "mJy/beam");
  let d = ""; vel.forEach((v, c) => { if (v < lo || v > hi) return; const half = (vel[1] - vel[0]) / 2; d += `${d ? "L" : "M"}${X(v - half).toFixed(1)} ${Y(sp[c]).toFixed(1)}L${X(v + half).toFixed(1)} ${Y(sp[c]).toFixed(1)}`; });
  svgEl("path", {d, fill: "none", stroke: css("--muted"), "stroke-width": 1.2}, svg);
  const xs = Array.from({length: 300}, (_, k) => lo + (hi - lo) * k / 299);
  path(svg, xs, xs.map(v => gh(v, p[3], p[4], p[5], p[6])), X, Y, css("--ink"), 1.8);
  path(svg, xs, xs.map(v => gh(v, p[8], p[9], p[10], 0)), X, Y, css("--app"), 1.1);
  const marks = [["Gauss-Hermite peak", p[7], css("--ink")], ["GH centre", p[4], css("--model")], ["Gaussian", p[9], css("--app")], ["brightest channel", p[11], css("--muted")], ["first moment", p[12], css("--pub")]];
  marks.forEach(([lab, v, col], k) => { if (v < -9000) return; svgEl("line", {x1: X(v), x2: X(v), y1: Y(bot), y2: Y(top * 0.98), stroke: col, "stroke-width": 1}, svg); });
  $("spech").textContent = `21-cm spectrum ${(-x / 60).toFixed(2)}′ east, ${(y / 60).toFixed(2)}′ north of the centre`;
  $("specnote").textContent = marks.filter(m => m[1] > -9000).map(m => `${m[0]} ${m[1].toFixed(1)}`).join(" · ") + " km/s (h3 = " + p[6].toFixed(2) + ")";
}
function canvasXY(ev, cv) {
  const r = cv.getBoundingClientRect(), {nx, ny, step, x0, y0} = G.map;
  const i = Math.floor((ev.clientX - r.left) / r.width * nx), j = ny - 1 - Math.floor((ev.clientY - r.top) / r.height * ny);
  return {i, j, x: x0 + i * step, y: y0 + j * step, k: j * nx + i};
}
for (const id of ["obs", "mod", "res"]) {
  const cv = $(id);
  cv.addEventListener("mousemove", ev => {
    const {i, j, x, y, k} = canvasXY(ev, cv), tip = $("tip");
    if (i < 0 || j < 0 || i >= G.map.nx || j >= G.map.ny || V[k] !== V[k]) { tip.style.display = "none"; return; }
    tip.innerHTML = `${(-x / 60).toFixed(2)}′ E, ${(y / 60).toFixed(2)}′ N<br>measured ${(V[k] + G.vsys).toFixed(1)} km/s<br>model ${(VM[k] + G.vsys).toFixed(1)}<br>residual ${VR[k].toFixed(1)}`;
    tip.style.display = "block"; tip.style.left = (ev.clientX + 14) + "px"; tip.style.top = (ev.clientY + 14) + "px";
  });
  cv.addEventListener("mouseleave", () => $("tip").style.display = "none");
  cv.addEventListener("click", ev => { const {x, y} = canvasXY(ev, cv); showSpectrum(x, y); });
}
for (const k of ["incl", "pa", "dv", "vflat", "rturn"]) $(k).addEventListener("input", () => { state[k] = +$(k).value; update(); });
$("curve").addEventListener("change", () => update());
$("reset").addEventListener("click", reset);
$("gal").addEventListener("change", () => { load(+$("gal").value); table(); });
function table() {
  const c = G.curve, kpa = G.kpcPerArcsec;
  let h = "<table><tr><th>R (arcsec)</th><th>R (kpc)</th><th>V (km/s)</th><th>±</th><th>approaching</th><th>receding</th><th>PA</th><th>incl.</th></tr>";
  c.r.forEach((r, k) => { h += `<tr><td>${r}</td><td>${(r * kpa).toFixed(2)}</td><td>${c.v[k]}</td><td>${c.err[k]}</td><td>${c.app[k]}</td><td>${c.rec[k]}</td><td>${c.pa[k]}</td><td>${c.incl[k]}</td></tr>`; });
  $("tbl").innerHTML = h + "</table>";
}
function restyle() { LUT = lut(); const st = (dark() ? CMAPS.dark : CMAPS.light); $("ramp").style.background = `linear-gradient(to right, ${st.join(",")})`; if (G) update(); }
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", restyle);
LUT = lut(); restyle(); load(0); table();
const hash = location.hash.slice(1).toUpperCase().replace(/\s/g, "");
DATA.forEach((g, i) => { if (g.key === hash) { $("gal").value = i; load(i); table(); } });
</script>
</body>
</html>
"""
