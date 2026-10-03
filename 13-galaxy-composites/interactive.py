"""A hands-on page: make your own colour composite from the real data.

``write_explorer`` writes ``figures/composite_explorer.html``, one
self-contained file that works offline. It carries four galaxies' calibrated
g, r, i images (block-averaged to under 500 pixels, gzip-compressed and
inlined), and re-renders them in the browser as you change the stretch, the
white point and the filter-to-channel mapping, or switch to the eye's
colours. Hovering reads out a pixel's surface brightness, colour and
luminance, and whether an eye could see its colour at all.

The page does exactly the arithmetic ``composite.py`` does; the same numbers
(the colour matrix, the white-balance weights, the vision thresholds) are
written into it from the analysis, not typed in.
"""

from __future__ import annotations

import base64
import gzip
import json

import numpy as np

import color_physics as cp
from config import (
    BANDS,
    FIGURE_DIR,
    K_M_PHOTOPIC,
    LUPTON_Q,
    LUPTON_SATURATION,
    MESOPIC_LOWER_CD_M2,
    MESOPIC_UPPER_CD_M2,
    galaxy_by_key,
)

EXPLORER_GALAXIES = ("m51", "m81", "m104")
MAX_PX = 480

#: asinh companding of surface brightness into 12 bits (decoded in the page).
#: One step is 0.4% of the value for bright pixels and 0.0002 nanomaggies/arcsec^2
#: near zero -- both far below the images' own noise -- and 12-bit codes
#: compress much better than 16-bit ones, whose lowest bits are pure noise.
ENC_A, ENC_LO, ENC_HI, ENC_LEVELS = 0.05, -3.0, 20000.0, 4095


def _encode(img: np.ndarray) -> str:
    lo, hi = np.arcsinh(ENC_LO / ENC_A), np.arcsinh(ENC_HI / ENC_A)
    v = (np.arcsinh(np.clip(img, ENC_LO, ENC_HI) / ENC_A) - lo) / (hi - lo)
    code = np.round(v * ENC_LEVELS).astype("<u2")
    return base64.b64encode(gzip.compress(code.tobytes(), compresslevel=9)).decode("ascii")


def _block_mean(img: np.ndarray, k: int) -> np.ndarray:
    h, w = (img.shape[0] // k) * k, (img.shape[1] // k) * k
    return img[:h, :w].reshape(h // k, k, w // k, k).mean(axis=(1, 3))


def write_explorer(ref: dict, eye, frames: dict) -> str:
    """``frames``: key -> imaging.Frame for the galaxies in EXPLORER_GALAXIES."""
    galaxies = []
    for key in EXPLORER_GALAXIES:
        if key not in frames:
            continue
        fr = frames[key]
        g = galaxy_by_key(key)
        k = int(np.ceil(fr.shape[0] / MAX_PX))
        imgs = {b: _block_mean(fr.sb[b], k)[::-1] for b in BANDS}  # north up, row 0 at top
        n = imgs["r"].shape[0]
        galaxies.append({
            "key": key, "name": g.name, "common": g.common, "kind": g.kind, "why": g.why,
            "n": n, "pixscale": g.pixscale * k, "fov": g.fov_arcmin,
            "data": {b: _encode(imgs[b]) for b in BANDS},
        })
    sun_to_d65 = cp.M_XYZ_TO_SRGB @ cp.chromatic_adaptation(ref["sun_white_xyz"], cp.d65_white_xyz())
    temps, luv = cp.planckian_locus()
    pick = np.unique(np.linspace(0, len(temps) - 1, 400).astype(int))  # temps are log-spaced already
    payload = {
        "galaxies": galaxies,
        "enc": {"a": ENC_A, "lo": ENC_LO, "hi": ENC_HI, "levels": ENC_LEVELS},
        "weights": {name: ref["white_balance"][name]["weights"] for name in ("ours", "viewer", "flat")},
        "neutral": {name: ref["white_balance"][name]["neutral_color"] for name in ("ours", "viewer", "flat")},
        "matrix": np.asarray(eye.color.matrix).tolist(),  # fluxes (maggies) @ matrix -> XYZ
        "sun_to_d65_srgb": sun_to_d65.tolist(),
        "xyz_to_srgb": cp.M_XYZ_TO_SRGB.tolist(),
        "y_per_intensity": ref["y_per_intensity"],
        "lum_factor": K_M_PHOTOPIC * cp.ERG_S_CM2_IN_W_M2 / cp.ARCSEC2_SR,
        "mesopic": [MESOPIC_LOWER_CD_M2, MESOPIC_UPPER_CD_M2],
        "Q": LUPTON_Q, "saturation": LUPTON_SATURATION,
        "locus": {"T": temps[pick].tolist(), "u": luv[pick, 0].tolist(), "v": luv[pick, 1].tolist()},
    }
    html = TEMPLATE.replace("/*PAYLOAD*/null", json.dumps(payload, separators=(",", ":")))
    out = FIGURE_DIR / "composite_explorer.html"
    out.write_text(html, encoding="utf-8")
    return str(out)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Galaxy Composite Explorer</title>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --rule: #e1e0d9; --accent: #eb6834; --accent-ink: #b8461c; --well: #000000;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
    --rule: #2c2c2a; --accent: #d95926; --accent-ink: #f08a5d; --well: #000000;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
  --rule: #2c2c2a; --accent: #d95926; --accent-ink: #f08a5d; --well: #000000;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
.lede { color: var(--ink2); margin: 0 0 18px; max-width: 70ch; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 360px; gap: 22px; align-items: start; }
@media (max-width: 860px) { .layout { grid-template-columns: 1fr; } }
.stage { background: var(--well); border-radius: 6px; overflow: hidden; position: relative; }
canvas { display: block; width: 100%; height: auto; image-rendering: auto; cursor: crosshair; }
.readout { font-variant-numeric: tabular-nums; color: var(--ink2); font-size: .85rem; min-height: 3.2em;
  margin-top: 8px; }
.readout b { color: var(--ink); font-weight: 600; }
.chip { display: inline-block; width: .9em; height: .9em; border-radius: 2px; vertical-align: -1px;
  border: 1px solid var(--rule); margin-right: 4px; }
.panel { background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; padding: 14px 14px 6px; }
fieldset { border: 0; padding: 0; margin: 0 0 14px; }
legend { font-weight: 600; font-size: .9rem; margin-bottom: 6px; }
.seg { display: flex; flex-wrap: wrap; gap: 6px; }
.seg button { font: inherit; font-size: .82rem; padding: 5px 9px; border-radius: 5px; cursor: pointer;
  border: 1px solid var(--rule); background: transparent; color: var(--ink); }
.seg button[aria-pressed="true"] { border-color: var(--accent); color: var(--accent-ink);
  background: color-mix(in srgb, var(--accent) 12%, transparent); }
label.row { display: grid; grid-template-columns: 1fr auto; gap: 8px; font-size: .84rem; color: var(--ink2);
  margin: 6px 0 2px; }
input[type=range] { width: 100%; accent-color: var(--accent); }
select { font: inherit; font-size: .85rem; width: 100%; padding: 4px; background: var(--surface);
  color: var(--ink); border: 1px solid var(--rule); border-radius: 5px; }
.note { font-size: .8rem; color: var(--muted); margin: 4px 0 0; }
.about { color: var(--ink2); font-size: .9rem; margin-top: 14px; max-width: 75ch; }
.about h2 { font-size: 1rem; color: var(--ink); margin: 16px 0 4px; }
.theme { position: absolute; right: 16px; top: 18px; }
header { position: relative; }
</style>
</head>
<body>
<main>
<header>
  <h1>Make your own galaxy picture</h1>
  <p class="lede">Real, calibrated SDSS images in three filters. Every choice that turns them into a colour
  picture is a control on this page: the stretch, what counts as white, which filter goes to which colour,
  or none of that, and the colours a human eye would see instead.</p>
  <div class="theme seg"><button id="themeBtn" type="button" aria-label="Toggle light or dark theme">◐ theme</button></div>
</header>
<div class="layout">
  <div>
    <div class="stage"><canvas id="cv" width="480" height="480" aria-label="Galaxy image"></canvas></div>
    <div class="readout" id="readout" aria-live="polite">Hover over the picture to read out a pixel.</div>
    <div class="about" id="about"></div>
  </div>
  <div class="panel">
    <fieldset><legend>Galaxy</legend><div class="seg" id="galBtns"></div></fieldset>
    <fieldset><legend>Colour</legend>
      <div class="seg" id="modeBtns">
        <button type="button" data-mode="survey" aria-pressed="true">Filters → channels</button>
        <button type="button" data-mode="eye" aria-pressed="false">What an eye would see</button>
      </div>
      <p class="note" id="modeNote"></p>
    </fieldset>
    <fieldset id="wbSet"><legend>White point (what shows as grey)</legend>
      <div class="seg" id="wbBtns">
        <button type="button" data-wb="ours" aria-pressed="true">The Sun</button>
        <button type="button" data-wb="viewer" aria-pressed="false">Survey viewer</button>
        <button type="button" data-wb="flat" aria-pressed="false">Equal flux</button>
      </div>
      <label class="row" for="wg"><span>g weight (blue)</span><output id="wgo"></output></label>
      <input type="range" id="wg" min="-1" max="1" step="0.01">
      <label class="row" for="wi"><span>i weight (red)</span><output id="wio"></output></label>
      <input type="range" id="wi" min="-1" max="1" step="0.01">
      <p class="note" id="wbNote"></p>
    </fieldset>
    <fieldset id="mapSet"><legend>Filter → channel</legend>
      <select id="mapping">
        <option value="irg">i → red, r → green, g → blue (the usual)</option>
        <option value="gri">g → red, r → green, i → blue (reversed)</option>
        <option value="rig">r → red, i → green, g → blue</option>
      </select>
    </fieldset>
    <fieldset><legend>Stretch</legend>
      <div class="seg" id="stretchBtns">
        <button type="button" data-st="lupton" aria-pressed="true">Arcsinh, total light</button>
        <button type="button" data-st="perchan" aria-pressed="false">Arcsinh, per channel</button>
        <button type="button" data-st="linear" aria-pressed="false">Linear</button>
        <button type="button" data-st="log" aria-pressed="false">Log</button>
      </div>
      <label class="row" for="sat"><span>Saturation (surface brightness)</span><output id="sato"></output></label>
      <input type="range" id="sat" min="1" max="3.8" step="0.01">
      <label class="row" for="q"><span>Q (how much range is squeezed in)</span><output id="qo"></output></label>
      <input type="range" id="q" min="1" max="14" step="0.1">
      <p class="note">Saturation: the surface brightness shown as full brightness. Q: the arcsinh curve is linear
      below saturation / sinh(Q) and logarithmic above it.</p>
      <div class="seg"><button type="button" id="resetBtn">Reset to the project's settings</button></div>
    </fieldset>
  </div>
</div>
</main>
<script>
const P = /*PAYLOAD*/null;
const $ = (id) => document.getElementById(id);
const state = { gal: 0, mode: "survey", wb: "ours", st: "lupton", map: "irg",
  sat: P.saturation, Q: P.Q, w: {...P.weights.ours} };
const cache = {};

async function decode(b64) {
  const bin = Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));
  const ds = new DecompressionStream("gzip");
  const buf = await new Response(new Blob([bin]).stream().pipeThrough(ds)).arrayBuffer();
  const code = new Uint16Array(buf);
  const lo = Math.asinh(P.enc.lo / P.enc.a), hi = Math.asinh(P.enc.hi / P.enc.a);
  const out = new Float32Array(code.length);
  for (let k = 0; k < code.length; k++) out[k] = P.enc.a * Math.sinh(lo + (code[k] / P.enc.levels) * (hi - lo));
  return out;
}
async function galaxyData(i) {
  if (cache[i]) return cache[i];
  const g = P.galaxies[i];
  const [G, R, I] = await Promise.all(["g", "r", "i"].map((b) => decode(g.data[b])));
  cache[i] = { g: G, r: R, i: I };
  return cache[i];
}
const enc = (v) => (v <= 0.0031308 ? 12.92 * v : 1.055 * Math.pow(v, 1 / 2.4) - 0.055);
const dec = (v) => (v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4));
const asinhCurve = (x, s, sat) => Math.asinh(x / s) / Math.asinh(sat / s);

function channelsFor(map, sg, sr, si, w) {
  const vg = sg * w.g, vr = sr * w.r, vi = si * w.i;
  if (map === "gri") return [vg, vr, vi];
  if (map === "rig") return [vr, vi, vg];
  return [vi, vr, vg];
}

function render() {
  const g = P.galaxies[state.gal];
  const d = cache[state.gal];
  if (!d) return;
  const cv = $("cv");
  cv.width = g.n; cv.height = g.n;
  const ctx = cv.getContext("2d");
  const im = ctx.createImageData(g.n, g.n);
  const px = im.data;
  const sat = state.sat, s = sat / Math.sinh(state.Q);
  const M = P.matrix, A = P.sun_to_d65_srgb, yref = P.y_per_intensity;
  const floor = 0.02, lf = Math.log10(floor), ls = Math.log10(sat);
  for (let k = 0; k < d.r.length; k++) {
    let r = 0, gg = 0, b = 0;
    if (state.mode === "survey") {
      let [R, G, B] = channelsFor(state.map, d.g[k], d.r[k], d.i[k], state.w);
      if (state.st === "lupton") {
        const I = (R + G + B) / 3;
        const f = I > 0 ? asinhCurve(I, s, sat) / I : 0;
        R = Math.max(0, R * f); G = Math.max(0, G * f); B = Math.max(0, B * f);
        const m = Math.max(R, G, B); if (m > 1) { R /= m; G /= m; B /= m; }
      } else if (state.st === "perchan") {
        R = Math.min(1, Math.max(0, asinhCurve(Math.max(R, 0), s, sat)));
        G = Math.min(1, Math.max(0, asinhCurve(Math.max(G, 0), s, sat)));
        B = Math.min(1, Math.max(0, asinhCurve(Math.max(B, 0), s, sat)));
      } else if (state.st === "linear") {
        R = Math.min(1, Math.max(0, R / sat)); G = Math.min(1, Math.max(0, G / sat)); B = Math.min(1, Math.max(0, B / sat));
      } else {
        const L = (x) => Math.min(1, Math.max(0, (Math.log10(Math.max(x, floor)) - lf) / (ls - lf)));
        R = L(R); G = L(G); B = L(B);
      }
      r = R; gg = G; b = B;
    } else {
      const f0 = d.g[k] * 1e-9, f1 = d.r[k] * 1e-9, f2 = d.i[k] * 1e-9;
      const X = (f0 * M[0][0] + f1 * M[1][0] + f2 * M[2][0]) / yref;
      const Y = (f0 * M[0][1] + f1 * M[1][1] + f2 * M[2][1]) / yref;
      const Z = (f0 * M[0][2] + f1 * M[1][2] + f2 * M[2][2]) / yref;
      if (Y > 0) {
        let R = A[0][0] * X + A[0][1] * Y + A[0][2] * Z;
        let G = A[1][0] * X + A[1][1] * Y + A[1][2] * Z;
        let B = A[2][0] * X + A[2][1] * Y + A[2][2] * Z;
        const lo = Math.min(R, G, B);
        if (lo < 0) { const t = Y / (Y - lo); R = Y + t * (R - Y); G = Y + t * (G - Y); B = Y + t * (B - Y); }
        let target;
        if (state.st === "lupton" || state.st === "perchan") target = dec(Math.min(1, asinhCurve(Y, s, sat)));
        else if (state.st === "linear") target = dec(Math.min(1, Y / sat));
        else target = dec(Math.min(1, Math.max(0, (Math.log10(Math.max(Y, floor)) - lf) / (ls - lf))));
        const sc = target / Y;
        R = Math.max(0, R * sc); G = Math.max(0, G * sc); B = Math.max(0, B * sc);
        const m = Math.max(R, G, B); if (m > 1) { R /= m; G /= m; B /= m; }
        r = enc(R); gg = enc(G); b = enc(B);
      }
    }
    const o = 4 * k;
    px[o] = 255 * r; px[o + 1] = 255 * gg; px[o + 2] = 255 * b; px[o + 3] = 255;
  }
  ctx.putImageData(im, 0, 0);
}
let pending = false;
function schedule() { if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; render(); }); } }

function cct(X, Y, Z) {
  const dd = X + 15 * Y + 3 * Z, u = 4 * X / dd, v = 6 * Y / dd;
  let best = 0, bd = 1e9;
  for (let k = 0; k < P.locus.T.length; k++) {
    const e = Math.hypot(P.locus.u[k] - u, P.locus.v[k] - v);
    if (e < bd) { bd = e; best = k; }
  }
  return P.locus.T[best];
}
function readout(ev) {
  const g = P.galaxies[state.gal], d = cache[state.gal];
  if (!d) return;
  const rect = $("cv").getBoundingClientRect();
  const x = Math.floor((ev.clientX - rect.left) / rect.width * g.n);
  const y = Math.floor((ev.clientY - rect.top) / rect.height * g.n);
  if (x < 0 || y < 0 || x >= g.n || y >= g.n) return;
  const k = y * g.n + x, sg = d.g[k], sr = d.r[k], si = d.i[k];
  if (sr <= 0 || sg <= 0) { $("readout").innerHTML = "Sky: no measurable light here (noise)."; return; }
  const mu = 22.5 - 2.5 * Math.log10(sr), gr = -2.5 * Math.log10(sg / sr), ri = -2.5 * Math.log10(sr / Math.max(si, 1e-9));
  const M = P.matrix, f = [sg * 1e-9, sr * 1e-9, si * 1e-9];
  const X = f[0] * M[0][0] + f[1] * M[1][0] + f[2] * M[2][0];
  const Y = f[0] * M[0][1] + f[1] * M[1][1] + f[2] * M[2][1];
  const Z = f[0] * M[0][2] + f[1] * M[1][2] + f[2] * M[2][2];
  const L = Y * P.lum_factor;
  const regime = L < P.mesopic[0] ? "rods only, no colour vision" : L < P.mesopic[1] ? "mesopic: faint colour" : "photopic: full colour";
  const A = P.sun_to_d65_srgb;
  let rgb = [0, 1, 2].map((j) => A[j][0] * X + A[j][1] * Y + A[j][2] * Z);
  const lo = Math.min(...rgb); if (lo < 0) rgb = rgb.map((c) => c - lo);
  const mx = Math.max(...rgb); rgb = rgb.map((c) => Math.round(255 * enc(c / mx)));
  const temp = gr > -1 && gr < 2 ? cct(X, Y, Z) : NaN;
  $("readout").innerHTML =
    `<b>${mu.toFixed(2)}</b> r mag/arcsec² · g−r <b>${gr.toFixed(2)}</b> · r−i <b>${ri.toFixed(2)}</b><br>` +
    `eye colour <span class="chip" style="background:rgb(${rgb.join(",")})"></span>` +
    `${isFinite(temp) ? "≈ " + Math.round(temp / 50) * 50 + " K" : ""} · luminance <b>${L.toExponential(1)}</b> cd/m² — ${regime}`;
}

function pressed(groupId, attr, value) {
  for (const b of $(groupId).querySelectorAll("button")) b.setAttribute("aria-pressed", String(b.dataset[attr] === value));
}
function syncControls() {
  $("wg").value = Math.log10(state.w.g); $("wi").value = Math.log10(state.w.i);
  $("wgo").textContent = state.w.g.toFixed(2); $("wio").textContent = state.w.i.toFixed(2);
  $("sat").value = Math.log10(state.sat); $("q").value = state.Q;
  $("sato").textContent = `${state.sat.toFixed(0)} nmgy/arcsec² (${(22.5 - 2.5 * Math.log10(state.sat)).toFixed(1)} mag)`;
  $("qo").textContent = state.Q.toFixed(1);
  const gr = 2.5 * Math.log10(state.w.g / state.w.r), ri = 2.5 * Math.log10(state.w.r / state.w.i);
  $("wbNote").textContent = `Grey here means a source with g−r = ${gr.toFixed(2)}, r−i = ${ri.toFixed(2)}. ` +
    `The Sun has ${P.neutral.ours["g-r"].toFixed(2)} and ${P.neutral.ours["r-i"].toFixed(2)}.`;
  const eye = state.mode === "eye";
  $("wbSet").style.opacity = eye ? 0.45 : 1; $("mapSet").style.opacity = eye ? 0.45 : 1;
  for (const el of [$("wg"), $("wi"), $("mapping"), ...$("wbBtns").querySelectorAll("button")]) el.disabled = eye;
  $("modeNote").textContent = eye
    ? "Each pixel's g, r, i fluxes → the CIE XYZ an eye would register (a matrix fitted to 129 real galaxy spectra) → your screen, with sunlight as white. White point and filter mapping no longer apply."
    : "Each filter's image is weighted and placed in one screen channel, the way survey pictures are made.";
  const g = P.galaxies[state.gal];
  $("about").innerHTML = `<h2>${g.name} · ${g.common}</h2><p>${g.kind[0].toUpperCase() + g.kind.slice(1)}. ${g.why}</p>` +
    `<p>${g.fov}′ across, SDSS DR8 g, r, i; shown here block-averaged to ${g.pixscale.toFixed(1)}″ pixels.</p>`;
}

function build() {
  P.galaxies.forEach((g, i) => {
    const b = document.createElement("button");
    b.type = "button"; b.textContent = g.name; b.dataset.gal = String(i);
    b.setAttribute("aria-pressed", String(i === state.gal));
    b.onclick = async () => { state.gal = i; pressed("galBtns", "gal", String(i)); syncControls(); await galaxyData(i); schedule(); };
    $("galBtns").appendChild(b);
  });
  for (const b of $("modeBtns").querySelectorAll("button"))
    b.onclick = () => { state.mode = b.dataset.mode; pressed("modeBtns", "mode", state.mode); syncControls(); schedule(); };
  for (const b of $("wbBtns").querySelectorAll("button"))
    b.onclick = () => { state.wb = b.dataset.wb; state.w = {...P.weights[state.wb]}; pressed("wbBtns", "wb", state.wb); syncControls(); schedule(); };
  for (const b of $("stretchBtns").querySelectorAll("button"))
    b.onclick = () => { state.st = b.dataset.st; pressed("stretchBtns", "st", state.st); schedule(); };
  $("wg").oninput = (e) => { state.w.g = 10 ** +e.target.value; pressed("wbBtns", "wb", ""); syncControls(); schedule(); };
  $("wi").oninput = (e) => { state.w.i = 10 ** +e.target.value; pressed("wbBtns", "wb", ""); syncControls(); schedule(); };
  $("sat").oninput = (e) => { state.sat = 10 ** +e.target.value; syncControls(); schedule(); };
  $("q").oninput = (e) => { state.Q = +e.target.value; syncControls(); schedule(); };
  $("mapping").onchange = (e) => { state.map = e.target.value; schedule(); };
  $("resetBtn").onclick = () => {
    Object.assign(state, { mode: "survey", wb: "ours", st: "lupton", map: "irg", sat: P.saturation, Q: P.Q, w: {...P.weights.ours} });
    $("mapping").value = "irg"; pressed("modeBtns", "mode", "survey"); pressed("wbBtns", "wb", "ours");
    pressed("stretchBtns", "st", "lupton"); syncControls(); schedule();
  };
  $("cv").addEventListener("mousemove", readout);
  $("themeBtn").onclick = () => {
    const root = document.documentElement;
    const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    root.dataset.theme = dark ? "light" : "dark";
  };
}
build(); syncControls();
galaxyData(0).then(schedule);
</script>
</body>
</html>
"""
