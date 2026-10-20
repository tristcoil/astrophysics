"""Build figures/contact_binary_explorer.html: an offline, dependency-free page
where you shape a contact binary and watch its light curve form.

The page carries a JavaScript port of roche.py's ray tracer (same potential,
same gravity and limb darkening, same Planck intensities), the Kepler mean
light curve of KIC 9832227 and the OGLE light curves of V1309 Sco's
progenitor, so a model can be laid over real stars.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import literature as L

PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Contact Binary Explorer</title>
<style>
:root {
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --blue: #2a78d6; --orange: #eb6834; --violet: #4a3aa7; --aqua: #1baf7a;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --surface: #161614; --panel: #201f1c; --ink: #f0efea; --ink2: #b9b7ae; --muted: #8f8d85;
    --grid: #2e2d29; --axis: #4a4943; --blue: #5598e7; --orange: #f08055; --violet: #8f82e0; --aqua: #34c590;
  }
}
:root[data-theme="dark"] {
  --surface: #161614; --panel: #201f1c; --ink: #f0efea; --ink2: #b9b7ae; --muted: #8f8d85;
  --grid: #2e2d29; --axis: #4a4943; --blue: #5598e7; --orange: #f08055; --violet: #8f82e0; --aqua: #34c590;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, "DejaVu Sans", sans-serif; }
main { max-width: 1120px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 22px; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 18px; max-width: 800px; }
.grid { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1.25fr); gap: 18px; }
@media (max-width: 820px) { .grid { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border-radius: 10px; padding: 14px; }
canvas { width: 100%; height: auto; display: block; border-radius: 6px; }
#star { background: #0b0b0b; aspect-ratio: 1 / 1; }
#lc { aspect-ratio: auto; }
.controls { display: grid; grid-template-columns: 1fr; gap: 10px; margin-top: 12px; }
label { display: grid; grid-template-columns: 150px 1fr 64px; align-items: center; gap: 10px; font-size: 14px; color: var(--ink2); }
@media (max-width: 480px) { label { grid-template-columns: 110px 1fr 56px; } }
input[type=range] { width: 100%; accent-color: var(--blue); }
output { font-variant-numeric: tabular-nums; color: var(--ink); text-align: right; }
.row { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-top: 10px; }
select, button { font: inherit; font-size: 14px; color: var(--ink); background: var(--surface); border: 1px solid var(--axis);
  border-radius: 6px; padding: 5px 10px; }
button { cursor: pointer; }
.facts { font-size: 13px; color: var(--ink2); margin-top: 10px; }
.facts b { color: var(--ink); font-weight: 600; }
.note { font-size: 13px; color: var(--muted); margin-top: 14px; max-width: 900px; }
.legend { display: flex; gap: 16px; flex-wrap: wrap; font-size: 13px; color: var(--ink2); margin-top: 6px; }
.sw { display: inline-block; width: 18px; height: 3px; vertical-align: middle; margin-right: 6px; border-radius: 2px; }
#status { font-size: 13px; color: var(--muted); }
</style>
</head>
<body>
<main>
<h1>Build a contact binary</h1>
<p class="lede">Two stars that share one envelope: choose their mass ratio, how far they overfill their Roche lobes,
and the angle we see them from. The image is ray-traced from the Roche potential; the light curve is the sum of the
image around the orbit. Lay it over real stars below.</p>
<div class="grid">
  <div class="card">
    <canvas id="star" width="360" height="360" aria-label="Rendered image of the binary"></canvas>
    <div class="row">
      <button id="play" type="button">Pause</button>
      <span id="status"></span>
    </div>
    <div class="facts" id="facts"></div>
  </div>
  <div class="card">
    <canvas id="lc" width="640" height="400" aria-label="Light curve"></canvas>
    <div class="legend">
      <span><span class="sw" style="background:var(--blue)"></span>model (this page)</span>
      <span id="leg-data"><span class="sw" style="background:var(--orange)"></span><span id="leg-data-text">data</span></span>
      <span><span class="sw" style="background:var(--ink2)"></span>current phase</span>
    </div>
    <div class="controls">
      <label>mass ratio q <input id="q" type="range" min="0.05" max="1" step="0.001" value="0.228"><output id="qv"></output></label>
      <label>fill-out f <input id="f" type="range" min="0" max="0.98" step="0.01" value="0.43"><output id="fv"></output></label>
      <label>inclination <input id="i" type="range" min="20" max="90" step="0.1" value="53.2"><output id="iv"></output></label>
      <label>T2 / T1 <input id="tr" type="range" min="0.8" max="1.2" step="0.001" value="1.021"><output id="trv"></output></label>
      <label>phase <input id="ph" type="range" min="0" max="1" step="0.005" value="0.12"><output id="phv"></output></label>
    </div>
    <div class="row">
      <select id="preset" aria-label="Preset">
        <option value="kic">KIC 9832227 (Molnar et al. 2017)</option>
        <option value="v1309">V1309 Sco progenitor, 2002 (Zhu et al. 2016)</option>
        <option value="twins">equal twins, deep contact</option>
        <option value="extreme">extreme mass ratio, q = 0.07</option>
      </select>
      <select id="data" aria-label="Compare with">
        <option value="kic">compare: KIC 9832227, Kepler 2009-13</option>
        <option value="v2002">compare: V1309 Sco, OGLE 2002</option>
        <option value="v2004">compare: V1309 Sco, OGLE 2004</option>
        <option value="v2006">compare: V1309 Sco, OGLE 2006</option>
        <option value="v2007">compare: V1309 Sco, OGLE 2007</option>
        <option value="none">no data</option>
      </select>
    </div>
  </div>
</div>
<p class="note">How it is computed: Kopal's Roche potential for a circular, synchronised orbit; the surface is the
equipotential with fill-out f between the inner (L1) and outer (L2) critical surfaces. Each pixel's ray is marched to the
surface, where temperature follows gravity (T ∝ g<sup>0.08</sup>, Lucy 1967) from each star's polar temperature
(5800 K for the heavier star), intensity is Planck's law at 640 nm (Kepler band) or 800 nm (OGLE I, chosen with the V1309 Sco data), with
quadratic limb darkening. No starspots, no reflection: the real stars' unequal maxima (spots) are the part the model
cannot reproduce. V1309 Sco preset: Zhu et al.'s q = 0.094 and f = 0.90; the inclination, 69°, is our match to the 2002
amplitude. The same code (roche.py) reproduces Kepler's KIC 9832227 amplitude to 0.002 mag with Molnar et al.'s parameters.
Data: Kepler long cadence (our 14-harmonic template); OGLE-III I band (Tylenda et al. 2011), folded on our exponential
ephemeris, magnitudes relative to each season's median.</p>
</main>
<script>
const DATA = __DATA__;
const $ = (id) => document.getElementById(id);
const HC_K = 1.438776877e-2, BETA = 0.32;

function omega(x, y, z, q) {
  const r1 = Math.sqrt(x * x + y * y + z * z), r2 = Math.sqrt((x - 1) * (x - 1) + y * y + z * z);
  return 1 / r1 + q * (1 / r2 - x) + 0.5 * (1 + q) * (x * x + y * y);
}
function grad(x, y, z, q) {
  const r1 = Math.sqrt(x * x + y * y + z * z), r2 = Math.sqrt((x - 1) * (x - 1) + y * y + z * z);
  const a = r1 * r1 * r1, b = r2 * r2 * r2;
  return [-x / a - q * (x - 1) / b - q + (1 + q) * x, -y / a - q * y / b + (1 + q) * y, -z / a - q * z / b];
}
function dOdx(x, q) {
  return -Math.sign(x) / (x * x) - q * Math.sign(x - 1) / ((x - 1) * (x - 1)) - q + (1 + q) * x;
}
function bisect(fn, lo, hi) {
  let flo = fn(lo);
  for (let k = 0; k < 100; k++) {
    const mid = 0.5 * (lo + hi), fm = fn(mid);
    if ((fm > 0) === (flo > 0)) { lo = mid; flo = fm; } else { hi = mid; }
  }
  return 0.5 * (lo + hi);
}
function lagrange(q) {
  const e = 1e-9;
  const x1 = bisect((x) => dOdx(x, q), e, 1 - e), x2 = bisect((x) => dOdx(x, q), 1 + e, 3), x3 = bisect((x) => dOdx(x, q), -3, -e);
  return { L1: x1, L2: x2, L3: x3, o1: omega(x1, 0, 0, q), o2: omega(x2, 0, 0, q), o3: omega(x3, 0, 0, q) };
}
function surfaceRadius(q, om, cx, ux, uy, uz, rmax) {
  // step outwards from the star centre to the first exit (the outer centrifugal region also has Omega > om)
  let lo = 1e-6, r = 0.004;
  while (r < rmax) {
    if (omega(cx + r * ux, r * uy, r * uz, q) < om) break;
    lo = r; r += 0.004;
  }
  if (r >= rmax) return rmax;
  let hi = r;
  for (let k = 0; k < 40; k++) {
    const mid = 0.5 * (lo + hi);
    if (omega(cx + mid * ux, mid * uy, mid * uz, q) > om) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}
function makeModel(q, f, t1, t2, lam) {
  const lp = lagrange(q);
  const om = lp.o1 - f * (lp.o1 - lp.o2);
  const xr = 1 + surfaceRadius(q, om, 1, 1, 0, 0, lp.L2 - 1);
  const xl = -surfaceRadius(q, om, 0, -1, 0, 0, -lp.L3);
  const cx = 0.5 * (xl + xr);
  // radius: the farthest surface point from (cx, 0, 0), sampled over directions from each star
  let rb = 0.5 * (xr - xl);
  for (const c of [0, 1]) {
    for (let a = 0; a < 24; a++) for (let b = 1; b < 12; b++) {
      const th = Math.PI * b / 12, phi = 2 * Math.PI * a / 24;
      const ux = Math.cos(th), uy = Math.sin(th) * Math.cos(phi), uz = Math.sin(th) * Math.sin(phi);
      const r = surfaceRadius(q, om, c, ux, uy, uz, 1.2);
      const px = c + r * ux - cx, py = r * uy, pz = r * uz;
      rb = Math.max(rb, Math.sqrt(px * px + py * py + pz * pz));
    }
  }
  rb += 0.03;
  const gp = (c) => { const r = surfaceRadius(q, om, c, 0, 0, 1, 1.5); const g = grad(c, 0, r, q); return Math.hypot(g[0], g[1], g[2]); };
  return { q, f, om, lp, cx, rb, t1, t2, lam: lam * 1e-9, g1: gp(0), g2: gp(1), ld: lam > 700 ? [0.55, 0.20] : [0.40, 0.26] };
}
function view(incl, phase) {
  const i = incl * Math.PI / 180, p = 2 * Math.PI * phase;
  const s = [Math.sin(i) * Math.cos(p), -Math.sin(i) * Math.sin(p), Math.cos(i)];
  let e2 = [-s[0] * s[2], -s[1] * s[2], 1 - s[2] * s[2]];
  const n2 = Math.hypot(e2[0], e2[1], e2[2]);
  e2 = n2 > 1e-12 ? e2.map((v) => v / n2) : [1, 0, 0];
  const e1 = [e2[1] * s[2] - e2[2] * s[1], e2[2] * s[0] - e2[0] * s[2], e2[0] * s[1] - e2[1] * s[0]];
  return { s, e1, e2 };
}
function render(M, incl, phase, npix, out) {
  const { s, e1, e2 } = view(incl, phase);
  const q = M.q, om = M.om, rb = M.rb, step = Math.max(0.004, 2 * rb / 260), c2 = HC_K / M.lam;
  let total = 0;
  for (let j = 0; j < npix; j++) {
    const v = -rb + (2 * rb) * (j + 0.5) / npix;
    for (let k = 0; k < npix; k++) {
      const u = -rb + (2 * rb) * (k + 0.5) / npix;
      const d2 = u * u + v * v;
      let val = 0;
      if (d2 < rb * rb) {
        const h = Math.sqrt(rb * rb - d2);
        const ox = M.cx + u * e1[0] + v * e2[0] + h * s[0], oy = u * e1[1] + v * e2[1] + h * s[1], oz = u * e1[2] + v * e2[2] + h * s[2];
        const nst = Math.floor(2 * h / step) + 2;
        let tp = 0, t = 0, hit = false;
        for (let n = 1; n < nst; n++) {
          t = Math.min(n * step, 2 * h);
          if (omega(ox - t * s[0], oy - t * s[1], oz - t * s[2], q) >= om) {
            let lo = tp, hi = t;
            for (let b = 0; b < 22; b++) { const mid = 0.5 * (lo + hi);
              if (omega(ox - mid * s[0], oy - mid * s[1], oz - mid * s[2], q) >= om) hi = mid; else lo = mid; }
            t = 0.5 * (lo + hi); hit = true; break;
          }
          tp = t;
        }
        if (hit) {
          const px = ox - t * s[0], py = oy - t * s[1], pz = oz - t * s[2];
          const g = grad(px, py, pz, q), gm = Math.hypot(g[0], g[1], g[2]);
          let mu = -(g[0] * s[0] + g[1] * s[1] + g[2] * s[2]) / gm;
          if (mu <= 0) mu = 1e-6;
          const temp = px < M.lp.L1 ? M.t1 * Math.pow(gm / M.g1, 0.25 * BETA) : M.t2 * Math.pow(gm / M.g2, 0.25 * BETA);
          const ld = 1 - M.ld[0] * (1 - mu) - M.ld[1] * (1 - mu) * (1 - mu);
          val = ld / (Math.exp(c2 / temp) - 1);
        }
      }
      if (out) out[j * npix + k] = val;
      total += val;
    }
  }
  return total * (2 * rb / npix) * (2 * rb / npix);
}

// ---- drawing ------------------------------------------------------------------
const starCv = $("star"), sctx = starCv.getContext("2d");
const lcCv = $("lc"), lctx = lcCv.getContext("2d");
const NP = 240;
let model = null, curve = null, imgMax = 1, playing = true, lastT = 0;
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

function starColour(x) {
  // dark -> warm white, the same ramp as the project's figures
  const stops = [[11, 11, 11], [61, 42, 28], [168, 117, 74], [242, 201, 154], [255, 246, 234]];
  const t = Math.max(0, Math.min(1, x)) * (stops.length - 1), k = Math.min(stops.length - 2, Math.floor(t)), f = t - k;
  return stops[k].map((a, n) => Math.round(a + f * (stops[k + 1][n] - a)));
}
function drawStar(phase) {
  const buf = new Float64Array(NP * NP);
  render(model, +$("i").value, phase, NP, buf);
  const im = sctx.createImageData(NP, NP);
  for (let j = 0; j < NP; j++) for (let k = 0; k < NP; k++) {
    const v = buf[(NP - 1 - j) * NP + k] / imgMax;
    const c = v > 0 ? starColour(v) : [11, 11, 11];
    const p = 4 * (j * NP + k);
    im.data[p] = c[0]; im.data[p + 1] = c[1]; im.data[p + 2] = c[2]; im.data[p + 3] = 255;
  }
  const off = document.createElement("canvas"); off.width = NP; off.height = NP;
  off.getContext("2d").putImageData(im, 0, 0);
  sctx.imageSmoothingEnabled = true;
  sctx.drawImage(off, 0, 0, starCv.width, starCv.height);
}
function computeCurve() {
  const n = 41, ph = [], fl = [];
  for (let k = 0; k < n; k++) { const p = 0.5 * k / (n - 1); ph.push(p); fl.push(render(model, +$("i").value, p, 84, null)); }
  const mx = Math.max(...fl);
  curve = { ph, mag: fl.map((f) => -2.5 * Math.log10(f / mx)) };
  // brightest image for the colour scale
  const buf = new Float64Array(NP * NP);
  render(model, +$("i").value, 0.25, NP, buf);
  imgMax = Math.max(...buf) || 1;
}
function curveAt(p) {
  p = ((p % 1) + 1) % 1; if (p > 0.5) p = 1 - p;
  const x = p / 0.5 * (curve.ph.length - 1), k = Math.min(curve.ph.length - 2, Math.floor(x)), f = x - k;
  return curve.mag[k] * (1 - f) + curve.mag[k + 1] * f;
}
function sizeCurveCanvas() {
  // draw in CSS pixels at the device's resolution, so text stays legible on phones
  const dpr = window.devicePixelRatio || 1, w = lcCv.clientWidth || 640, h = Math.round(w * 0.62);
  lcCv.width = Math.round(w * dpr); lcCv.height = Math.round(h * dpr);
  lctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return [w, h];
}
let LCW = 640, LCH = 400;
function drawCurve(phase) {
  const W = LCW, H = LCH, m = { l: 62, r: 12, t: 14, b: 40 };
  lctx.clearRect(0, 0, W, H);
  lctx.fillStyle = css("--panel"); lctx.fillRect(0, 0, W, H);
  const sel = $("data").value, d = DATA.curves[sel];
  let lo = 0, hi = Math.max(0.05, ...curve.mag);
  if (d) { hi = Math.max(hi, ...d.mag); lo = Math.min(lo, ...d.mag); }
  hi += 0.03; lo -= 0.03;
  const X = (p) => m.l + (W - m.l - m.r) * p / 2, Y = (v) => m.t + (H - m.t - m.b) * (v - lo) / (hi - lo);
  lctx.strokeStyle = css("--grid"); lctx.lineWidth = 1; lctx.font = "12px system-ui, sans-serif"; lctx.fillStyle = css("--muted");
  for (let p = 0; p <= 2.0001; p += 0.25) { lctx.beginPath(); lctx.moveTo(X(p), m.t); lctx.lineTo(X(p), H - m.b); lctx.stroke();
    if (Math.abs(p * 2 - Math.round(p * 2)) < 1e-9) lctx.fillText(p.toFixed(1), X(p) - 8, H - m.b + 16); }
  const stepv = (hi - lo) > 0.5 ? 0.1 : 0.05;
  lctx.textAlign = "right";
  for (let v = Math.ceil(lo / stepv) * stepv; v <= hi; v += stepv) { lctx.beginPath(); lctx.moveTo(m.l, Y(v)); lctx.lineTo(W - m.r, Y(v)); lctx.stroke();
    lctx.fillText(v.toFixed(2), m.l - 6, Y(v) + 4); }
  lctx.textAlign = "left";
  lctx.fillStyle = css("--ink2");
  lctx.fillText("orbital phase", W / 2 - 40, H - 8);
  lctx.save(); lctx.translate(13, H / 2 + 55); lctx.rotate(-Math.PI / 2); lctx.fillText("mag below maximum", 0, 0); lctx.restore();
  if (d) {
    lctx.fillStyle = css("--orange");
    for (let k = 0; k < d.ph.length; k++) for (const sft of [0, 1]) {
      lctx.beginPath(); lctx.arc(X(d.ph[k] + sft), Y(d.mag[k]), d.dense ? 1.2 : 2.4, 0, 2 * Math.PI); lctx.fill(); }
  }
  lctx.strokeStyle = css("--blue"); lctx.lineWidth = 2.2; lctx.beginPath();
  for (let n = 0; n <= 400; n++) { const p = 2 * n / 400, x = X(p), y = Y(curveAt(p)); if (n) lctx.lineTo(x, y); else lctx.moveTo(x, y); }
  lctx.stroke();
  lctx.strokeStyle = css("--ink2"); lctx.lineWidth = 1; lctx.beginPath(); lctx.moveTo(X(phase), m.t); lctx.lineTo(X(phase), H - m.b); lctx.stroke();
  lctx.fillStyle = css("--ink"); lctx.beginPath(); lctx.arc(X(phase), Y(curveAt(phase)), 5, 0, 2 * Math.PI); lctx.fill();
}
function facts() {
  const lp = model.lp, amp = Math.max(...curve.mag);
  const dmin = curveAt(0) - curveAt(0.5);
  $("facts").innerHTML = `L1 at <b>${lp.L1.toFixed(3)}</b>, L2 at <b>${lp.L2.toFixed(3)}</b> separations; ` +
    `amplitude <b>${amp.toFixed(3)} mag</b>; deeper minimum at phase ${dmin >= 0 ? "0" : "0.5"} ` +
    `(${Math.abs(dmin).toFixed(3)} mag deeper). Phase 0: the heavier star is eclipsed.`;
}
function setOut() {
  $("qv").value = (+$("q").value).toFixed(3); $("fv").value = (+$("f").value).toFixed(2);
  $("iv").value = (+$("i").value).toFixed(1) + "°"; $("trv").value = (+$("tr").value).toFixed(3); $("phv").value = (+$("ph").value).toFixed(2);
}
let pending = null;
function rebuild() {
  setOut();
  $("status").textContent = "computing...";
  clearTimeout(pending);
  pending = setTimeout(() => {
    const lam = $("data").value.startsWith("v") ? 800 : 640;
    model = makeModel(+$("q").value, +$("f").value, 5800, 5800 * +$("tr").value, lam);
    computeCurve(); facts(); redraw(); $("status").textContent = "";
  }, 60);
}
function redraw() { const p = +$("ph").value; drawStar(p); drawCurve(p); }
const PRESETS = DATA.presets;
function applyPreset(k) {
  const p = PRESETS[k];
  $("q").value = p.q; $("f").value = p.f; $("i").value = p.i; $("tr").value = p.tr;
  if (p.data) $("data").value = p.data;
  $("leg-data-text").textContent = DATA.curves[$("data").value] ? DATA.curves[$("data").value].label : "";
  rebuild();
}
for (const id of ["q", "f", "i", "tr"]) $(id).addEventListener("input", rebuild);
$("ph").addEventListener("input", () => { setOut(); redraw(); });
$("preset").addEventListener("change", (e) => applyPreset(e.target.value));
$("data").addEventListener("change", () => { $("leg-data-text").textContent = DATA.curves[$("data").value] ? DATA.curves[$("data").value].label : ""; rebuild(); });
$("play").addEventListener("click", () => { playing = !playing; $("play").textContent = playing ? "Pause" : "Play"; });
function tick(ts) {
  if (playing && model && ts - lastT > 90) {
    lastT = ts;
    $("ph").value = ((+$("ph").value + 0.01) % 1).toFixed(3); setOut(); redraw();
  }
  requestAnimationFrame(tick);
}
[LCW, LCH] = sizeCurveCanvas();
window.addEventListener("resize", () => { [LCW, LCH] = sizeCurveCanvas(); if (curve) drawCurve(+$("ph").value); });
applyPreset("kic");
requestAnimationFrame(tick);
</script>
</body>
</html>
"""


def build() -> str:
    tpl_npz = np.load(C.DERIVED_DIR / "roche_kic_model.npz")
    ph = tpl_npz["phase"]
    kep = -2.5 * np.log10(tpl_npz["kepler"] / tpl_npz["kepler"].max())
    curves = {"kic": {"label": "KIC 9832227, Kepler 2009-2013 (mean)", "ph": [round(float(x), 4) for x in ph[:-1]],
                      "mag": [round(float(x), 4) for x in kep[:-1]], "dense": True}}
    vc = list(np.load(C.DERIVED_DIR / "v1309_curves.npz", allow_pickle=True)["curves"])
    for c in vc:
        yr = int(c["year"])
        if yr not in (2002, 2004, 2006, 2007):
            continue
        p = np.array(c["phase"])
        dm = np.array(c["dmag"])
        fit = np.array(c["fit_dmag"])
        top = fit.min()
        curves[f"v{yr}"] = {"label": f"V1309 Sco, OGLE {yr} ({c['n']} points)", "ph": [round(float(x), 4) for x in p],
                            "mag": [round(float(x - top), 4) for x in dm], "dense": False}
    m = L.MOLNAR2017
    presets = {
        "kic": {"q": m["q"], "f": m["fillout"], "i": m["incl"], "tr": round(m["t2"] / m["t1"], 3), "data": "kic"},
        "v1309": {"q": 0.094, "f": 0.90, "i": 69.0, "tr": 1.0, "data": "v2002"},
        "twins": {"q": 1.0, "f": 0.6, "i": 85.0, "tr": 1.0, "data": "none"},
        "extreme": {"q": 0.07, "f": 0.8, "i": 80.0, "tr": 1.0, "data": "none"},
    }
    html = PAGE.replace("__DATA__", json.dumps({"curves": curves, "presets": presets}, separators=(",", ":")))
    out = C.FIGURE_DIR / "contact_binary_explorer.html"
    out.write_text(html)
    return str(out)


if __name__ == "__main__":
    print(build())
