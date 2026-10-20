"""The disc explorer: deproject a disc by hand, then weigh the planet in its gap.

One self-contained HTML page (no network, no libraries). Pick a disc; two sliders tilt
and turn a trial geometry. The page re-samples the real ALMA image along circles of that
trial disc and unrolls it into (azimuth, radius): when the inclination and position
angle are right, every ring becomes a straight horizontal stripe and the averaged
profile is at its sharpest. A "ring straightness" score -- the scatter of the brightness
around each circle, averaged over radius -- falls as you approach the right geometry.

The second panel turns a gap's width into a planet mass with DSHARP VII's fitted
relation, Delta = A K'^B with K' = q (h/r)^-0.18 alpha^-0.31 (the same formula and table
as planets.py; checked in a headless browser: AS 209's 9-au gap gives 2.058 M_J on the page
and 2.058 in Python at the same inputs), letting the viscosity alpha and the grain size be varied.

The unrolled image uses the same transform as geometry.sky_from_polar. The page follows
the viewer's light or dark setting.
"""

from __future__ import annotations

import base64
import json

import numpy as np
from matplotlib import colormaps

import alma_images as A
import literature as L
import planets as P
from config import FIGURE_DIR, ZHANG_WIDTH_FITS, ZHANG_WIDTH_FITS_DSD2

EXPLORER_PATH = FIGURE_DIR / "disc_explorer.html"
EXPLORER_DISCS = ("AS209", "HD163296", "HLTau", "Elias24", "GWLup", "DoAr25", "HD143006", "Elias27", "IMLup")
N_PIX = 256


def _disc_payload(key: str, r: dict) -> dict:
    o = r["discs"][key]
    g = o["geometry"]
    d = r["distances"][key]["distance"]
    half = 1.2 * o["r_dust_au"] / d  # arcsec
    img = A.load_continuum(key)
    # resample onto an N x N grid centred on our fitted centre (x east, y north)
    from geometry import sample_image

    t = np.linspace(-half, half, N_PIX)
    xx, yy = np.meshgrid(-t, t)  # column 0 is east (left), row 0 is south (bottom)
    v = sample_image(img, g["x0"] + xx, g["y0"] + yy, order=1)
    peak = np.nanpercentile(v, 99.9)
    s = np.arcsinh(np.clip(v / peak, 0, None) / 0.03) / np.arcsinh(1 / 0.03)
    u8 = np.clip(np.nan_to_num(s) * 255, 0, 255).astype(np.uint8)[::-1]  # row 0 = north (top) for the canvas
    pub = r["geometry"][key]["published"]
    star = L.dsharp1_stars().get(key)
    gaps = []
    for p in r["planets"]["rows"]:
        if p.get("found") and p["disc"] == key:
            gaps.append(dict(label=f"{p['our_gap']}{' (+D90, B97)' if p['merged'] else ''}", r_au=p["r_gap_au"],
                             delta=p["delta"], st=[x if x is not None else 0 for x in p["st_used"]],
                             pub=p["masses"]["DSD1_0.001"]["pub"]))
    lstar = None
    if star is not None:
        lstar = star.lstar * (d / star.distance.value) ** 2
    return dict(
        key=key, name=o["name"], d=d, half_arcsec=half, n=N_PIX, img=base64.b64encode(u8.tobytes()).decode(),
        ours=dict(inc=g["inc"], pa=g["pa"] % 180), dsharp=dict(inc=pub["inc"], pa=pub["pa"] % 180, src=pub["source"]),
        rdust_au=o["r_dust_au"], beam_au=float(np.sqrt(o["beam"][0] * o["beam"][1]) * d),
        features=[dict(kind=f["kind"], r_au=f["r_au"], label=f["label"]) for f in o["features"] if f["secure"]],
        mstar=star.mstar if star else None, lstar=lstar, gaps=gaps,
    )


def build(r: dict) -> None:
    discs = [_disc_payload(k, r) for k in EXPLORER_DISCS if k in r["discs"]]
    lut = (np.array([colormaps["inferno"](i / 255)[:3] for i in range(256)]) * 255).astype(int).ravel().tolist()
    fits = {"DSD1": {repr(k): v for k, v in ZHANG_WIDTH_FITS.items() if k != "gas"},
            "DSD2": {repr(k): v for k, v in ZHANG_WIDTH_FITS_DSD2.items() if v is not None}}
    consts = dict(G=6.67430e-8, kB=1.380649e-16, mH=1.6735575e-24, mu=2.35, sigma=5.670374419e-5, Lsun=3.828e33,
                  Msun=1.98840987e33, au=1.495978707e13, phi=0.02, MsunInMJ=P.M_SUN_IN_MJUP,
                  MJinMearth=317.828)
    data = json.dumps(dict(discs=discs, lut=lut, fits=fits, c=consts, stunit=5.23e-4))
    html = PAGE.replace("__DATA__", data)
    EXPLORER_PATH.write_text(html)


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Disc Explorer</title>
<style>
:root {
  --bg: #fcfcfb; --panel: #ffffff; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --ours: #0b0b0b; --pub: #4a3aa7; --gap: #2a78d6; --ring: #eb6834;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
    --grid: #2c2c2a; --axis: #4a4945; --ours: #f2f1ec; --pub: #9085e9; --gap: #3987e5; --ring: #d95926;
  }
}
:root[data-theme="dark"] {
  --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
  --grid: #2c2c2a; --axis: #4a4945; --ours: #f2f1ec; --pub: #9085e9; --gap: #3987e5; --ring: #d95926;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink); font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
.wrap { max-width: 1200px; margin: 0 auto; padding: 18px 16px 28px; }
h1 { font-size: 20px; margin: 0 0 2px; }
h2 { font-size: 15px; margin: 18px 0 4px; }
.sub { color: var(--ink2); margin: 0 0 12px; max-width: 940px; }
.bar { display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: center; margin: 6px 0 10px; }
.bar button, .bar select { font: inherit; color: var(--ink); background: var(--panel); border: 1px solid var(--axis);
  border-radius: 6px; padding: 4px 10px; cursor: pointer; }
.controls { display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 6px 26px; margin: 4px 0 10px; }
.controls label { color: var(--ink2); display: grid; grid-template-columns: 9.5em minmax(0, 1fr) 5.2em; gap: 10px; align-items: center; }
.controls input { width: 100%; }
.controls output { font-variant-numeric: tabular-nums; color: var(--ink); text-align: right; }
.grid { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1.25fr); gap: 14px; }
@media (max-width: 860px) { .grid { grid-template-columns: minmax(0, 1fr); } .controls label { grid-template-columns: 7.5em 1fr 4.4em; } }
.card { background: var(--panel); border: 1px solid var(--grid); border-radius: 8px; padding: 8px 10px; }
.card h3 { font-size: 13px; margin: 2px 0 6px; color: var(--ink2); font-weight: 600; }
canvas { display: block; width: 100%; height: auto; border-radius: 4px; image-rendering: auto; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
svg text { fill: var(--muted); font-size: 11px; }
.score { font-variant-numeric: tabular-nums; font-size: 22px; font-weight: 600; }
.readout { width: 100%; border-collapse: collapse; margin-top: 8px; font-variant-numeric: tabular-nums; }
.readout td { padding: 3px 8px; border-bottom: 1px solid var(--grid); }
.readout td:first-child { color: var(--ink2); width: 55%; }
.note { color: var(--muted); font-size: 12px; margin-top: 12px; max-width: 1000px; }
.legend { color: var(--ink2); font-size: 12px; }
.sw { display: inline-block; width: 18px; height: 0; border-top: 2px solid; margin: 0 4px 0 10px; vertical-align: 3px; }
</style>
</head>
<body>
<div class="wrap">
<h1>Unroll a planet-forming disc</h1>
<p class="sub">ALMA sees each disc tilted, as an ellipse. Tilt it back (inclination) and turn it (position angle) until
the rings are circles: then the unrolled image on the right shows every ring as a straight horizontal stripe, and the
brightness profile is at its sharpest. Real DSHARP and HL Tau images; the score below drops as the geometry improves.</p>
<div class="bar">
  <label class="legend" for="disc">disc</label> <select id="disc"></select>
  <button id="ours">our fit</button> <button id="dsharp">published</button> <button id="faceon">no deprojection</button>
</div>
<div class="controls">
  <label for="inc">inclination (°)<input id="inc" type="range" min="0" max="80" step="0.1"><output id="incv"></output></label>
  <label for="pa">position angle (°)<input id="pa" type="range" min="0" max="180" step="0.1"><output id="pav"></output></label>
</div>
<div class="grid">
  <div class="card"><h3>The sky as ALMA saw it, with circles of the trial disc</h3><canvas id="sky" width="512" height="512"
    aria-label="ALMA image of the selected disc"></canvas></div>
  <div class="card"><h3>Unrolled: azimuth across, radius up</h3><canvas id="polar" width="720" height="360"
    aria-label="The disc unrolled into azimuth and radius"></canvas>
    <svg id="prof" viewBox="0 0 720 180" aria-label="Azimuthally averaged brightness profile"></svg>
    <div class="legend">ring straightness (lower is better): <span class="score" id="score"></span>
      <span class="sw" style="border-color: var(--gap)"></span>our gaps <span class="sw" style="border-color: var(--ring)"></span>our rings</div>
  </div>
</div>
<h2>Weigh the planet that could carve a gap</h2>
<p class="sub" id="gapintro"></p>
<div class="bar"><label class="legend" for="gap">gap</label> <select id="gap"></select>
  <label class="legend" for="model">largest grains</label> <select id="model">
    <option value="DSD1">0.1 mm (DSHARP VII's preferred)</option><option value="1mm">1 mm</option><option value="DSD2">1 cm</option></select></div>
<div class="controls">
  <label for="delta">gap width Δ<input id="delta" type="range" min="0.10" max="0.70" step="0.001"><output id="deltav"></output></label>
  <label for="alpha">viscosity log α<input id="alpha" type="range" min="-4" max="-2" step="0.05" value="-3"><output id="alphav"></output></label>
</div>
<table class="readout" id="planet"></table>
<p class="note">Image brightness is shown on an asinh stretch, and the profile and score are computed from those
stretched values, so they are qualitative. Planet masses follow Zhang et al. 2018 (DSHARP VII): Δ = A K′<sup>B</sup>,
K′ = q (h/r)<sup>−0.18</sup> α<sup>−0.31</sup>, with A and B for the grains' Stokes number (DSHARP VII's own choice per
gap), h/r from a passively heated disc (T = (φ L<sub>*</sub>/8πr²σ)<sup>1/4</sup>, φ = 0.02). A gap is not proof of a
planet: ice lines, magnetic winds and dust growth can make rings too.</p>
</div>
<script>
const D = __DATA__;
const $ = (id) => document.getElementById(id);
const sel = $("disc"), inc = $("inc"), pa = $("pa");
const sky = $("sky").getContext("2d"), pol = $("polar").getContext("2d");
let cur = null, pix = null;
D.discs.forEach((d, i) => { const o = document.createElement("option"); o.value = i; o.textContent = d.name; sel.appendChild(o); });
function css(v) { return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }
function decode(b64) { const s = atob(b64); const a = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) a[i] = s.charCodeAt(i); return a; }
// disc (r [au], theta [deg]) -> sky offsets (au; x east, y north), as geometry.sky_from_polar
function sky_from_polar(r, th, incd, pad) {
  const p = pad * Math.PI / 180, c = Math.cos(incd * Math.PI / 180), t = th * Math.PI / 180;
  const xm = r * Math.cos(t) * c, ym = r * Math.sin(t);
  return [xm * Math.cos(p) + ym * Math.sin(p), -xm * Math.sin(p) + ym * Math.cos(p)];
}
function sample(x_au, y_au) {  // bilinear, image row 0 = north, column 0 = east
  const n = cur.n, half = cur.half_arcsec * cur.d;
  const fx = (half - x_au) / (2 * half) * (n - 1), fy = (half - y_au) / (2 * half) * (n - 1);
  if (fx < 0 || fy < 0 || fx > n - 1 || fy > n - 1) return NaN;
  const i = Math.floor(fx), j = Math.floor(fy), di = fx - i, dj = fy - j;
  const i1 = Math.min(i + 1, n - 1), j1 = Math.min(j + 1, n - 1);
  return (pix[j * n + i] * (1 - di) * (1 - dj) + pix[j * n + i1] * di * (1 - dj) + pix[j1 * n + i] * (1 - di) * dj + pix[j1 * n + i1] * di * dj);
}
function drawSky() {
  const n = cur.n, img = sky.createImageData(n, n);
  for (let k = 0; k < n * n; k++) { const v = pix[k]; img.data.set([D.lut[3 * v], D.lut[3 * v + 1], D.lut[3 * v + 2], 255], 4 * k); }
  const off = new OffscreenCanvas(n, n); off.getContext("2d").putImageData(img, 0, 0);
  sky.imageSmoothingEnabled = true; sky.drawImage(off, 0, 0, 512, 512);
  const half = cur.half_arcsec * cur.d, s = 256 / half;
  sky.strokeStyle = "rgba(255,255,255,0.75)"; sky.lineWidth = 1;
  const radii = cur.features.length ? cur.features.map(f => f.r_au) : [0.25, 0.5, 0.75, 1].map(f => f * cur.rdust_au);
  for (const r of radii) {
    sky.beginPath();
    for (let t = -180; t <= 180; t += 3) { const [x, y] = sky_from_polar(r, t, +inc.value, +pa.value); const X = 256 - x * s, Y = 256 - y * s; t === -180 ? sky.moveTo(X, Y) : sky.lineTo(X, Y); }
    sky.stroke();
  }
  const [mx, my] = sky_from_polar(1.15 * cur.rdust_au, 90, +inc.value, +pa.value);
  sky.setLineDash([4, 4]); sky.beginPath(); sky.moveTo(256 + mx * s, 256 + my * s); sky.lineTo(256 - mx * s, 256 - my * s); sky.stroke(); sky.setLineDash([]);
  sky.fillStyle = "white"; sky.font = "12px system-ui"; sky.fillText("N ↑  E ←", 10, 20);
  sky.fillText(`${(2 * half).toFixed(0)} au across`, 10, 500);
}
function drawPolar() {
  const W = 720, H = 360, nt = 360, nr = 180, rmax = 1.15 * cur.rdust_au;
  const img = pol.createImageData(nt, nr), prof = new Float64Array(nr), sd = new Float64Array(nr);
  for (let ir = 0; ir < nr; ir++) {
    const r = (ir + 0.5) / nr * rmax; let s1 = 0, s2 = 0, m = 0;
    for (let it = 0; it < nt; it++) {
      const [x, y] = sky_from_polar(r, it - 180, +inc.value, +pa.value);
      const v = sample(x, y), row = nr - 1 - ir;
      const q = isNaN(v) ? 0 : Math.round(v);
      img.data.set([D.lut[3 * q], D.lut[3 * q + 1], D.lut[3 * q + 2], 255], 4 * (row * nt + it));
      if (!isNaN(v)) { s1 += v; s2 += v * v; m++; }
    }
    prof[ir] = m ? s1 / m : NaN; sd[ir] = m ? Math.sqrt(Math.max(s2 / m - prof[ir] ** 2, 0)) : NaN;
  }
  const off = new OffscreenCanvas(nt, nr); off.getContext("2d").putImageData(img, 0, 0);
  pol.imageSmoothingEnabled = true; pol.drawImage(off, 0, 0, W, H);
  pol.fillStyle = "white"; pol.font = "11px system-ui";
  pol.fillText("θ = −180°", 6, H - 6); pol.fillText("+180°", W - 40, H - 6); pol.fillText(`r = ${rmax.toFixed(0)} au`, 6, 14);
  // score: azimuthal scatter relative to the mean, averaged where the disc is bright
  let num = 0, den = 0;
  for (let ir = 2; ir < nr; ir++) if (prof[ir] > 25) { num += sd[ir]; den += prof[ir]; }
  $("score").textContent = den ? (100 * num / den).toFixed(1) + "%" : "–";
  drawProfile(prof, rmax);
}
function drawProfile(prof, rmax) {
  const svg = $("prof"), W = 720, H = 180, L = 36, B = 22, nr = prof.length;
  const max = Math.max(...prof.filter(v => !isNaN(v)));
  let path = "";
  for (let i = 0; i < nr; i++) { const x = L + (i + 0.5) / nr * (W - L - 6), y = H - B - prof[i] / max * (H - B - 10); path += (i ? "L" : "M") + x.toFixed(1) + "," + y.toFixed(1); }
  let marks = "";
  for (const f of cur.features) { const x = L + f.r_au / rmax * (W - L - 6); if (x < W) marks += `<line x1="${x}" x2="${x}" y1="10" y2="${H - B}" stroke="${f.kind === "D" ? css("--gap") : css("--ring")}" stroke-width="1" ${f.kind === "D" ? 'stroke-dasharray="3 3"' : ""}/>`; }
  let ticks = "";
  for (let r = 0; r <= rmax; r += (rmax > 150 ? 50 : 20)) { const x = L + r / rmax * (W - L - 6); ticks += `<text x="${x}" y="${H - 6}" text-anchor="middle">${r}</text>`; }
  svg.innerHTML = `<line x1="${L}" x2="${W - 6}" y1="${H - B}" y2="${H - B}" stroke="${css("--axis")}"/>${marks}${ticks}
    <text x="${W - 6}" y="${H - 6 - 12}" text-anchor="end">radius (au)</text><text x="4" y="14">brightness</text>
    <path d="${path}" fill="none" stroke="${css("--ours")}" stroke-width="2"/>`;
}
function update() { $("incv").textContent = (+inc.value).toFixed(1); $("pav").textContent = (+pa.value).toFixed(1); drawSky(); drawPolar(); }
function setGeom(g) { inc.value = g.inc; pa.value = g.pa; update(); }
function loadDisc() {
  cur = D.discs[+sel.value]; pix = decode(cur.img); setGeom(cur.ours);
  const gs = $("gap"); gs.innerHTML = "";
  cur.gaps.forEach((g, i) => { const o = document.createElement("option"); o.value = i; o.textContent = `${g.label} at ${g.r_au.toFixed(0)} au`; gs.appendChild(o); });
  if (!cur.gaps.length) { const o = document.createElement("option"); o.textContent = "(none in DSHARP VII's list: try AS 209 or HD 163296)"; gs.appendChild(o); }
  $("gapintro").textContent = cur.mstar ? `${cur.name}: star of ${cur.mstar.toFixed(2)} M☉ and ${cur.lstar.toFixed(1)} L☉ (DSHARP I, luminosity rescaled to the Gaia DR3 distance).`
    : `${cur.name}: no stellar parameters in DSHARP I, so no planet masses here.`;
  if (cur.gaps.length) { $("delta").value = cur.gaps[0].delta; }
  planet();
}
function fitFor(stUnits, model) {
  const st = stUnits * D.stunit, table = model === "DSD2" ? D.fits.DSD2 : D.fits.DSD1;
  let best = null;
  for (const k in table) { const d = Math.abs(Math.log10(+k) - Math.log10(st)); if (!best || d < best.d) best = { d, k: +k, v: table[k] }; }
  return best;
}
function planet() {
  const t = $("planet"); $("deltav").textContent = (+$("delta").value).toFixed(3); $("alphav").textContent = (+$("alpha").value).toFixed(2);
  if (!cur.gaps.length || !cur.mstar) { t.innerHTML = ""; return; }
  const g = cur.gaps[+$("gap").value || 0], C = D.c, model = $("model").value;
  const stU = model === "DSD1" ? g.st[0] : model === "1mm" ? g.st[1] : g.st[2];
  if (!stU) { t.innerHTML = `<tr><td>DSHARP VII found no solution for this gap with these grains</td><td></td></tr>`; return; }
  const f = fitFor(stU, model), A = f.v[0], B = f.v[1];
  const r = g.r_au * C.au, T = Math.pow(C.phi * cur.lstar * C.Lsun / (8 * Math.PI * r * r * C.sigma), 0.25);
  const cs = Math.sqrt(C.kB * T / (C.mu * C.mH)), vk = Math.sqrt(C.G * cur.mstar * C.Msun / r), hr = cs / vk;
  const alpha = Math.pow(10, +$("alpha").value), delta = +$("delta").value;
  const q = Math.pow(delta / A, 1 / B) * Math.pow(hr, 0.18) * Math.pow(alpha, 0.31);
  const mj = q * cur.mstar * C.MsunInMJ, me = mj * C.MJinMearth;
  const like = me < 5 ? "a super-Earth" : me < 25 ? "Neptune/Uranus-like" : mj < 0.5 ? "Saturn-like" : mj < 2.5 ? "Jupiter-like" : "a super-Jupiter";
  t.innerHTML = `<tr><td>midplane temperature at ${g.r_au.toFixed(0)} au</td><td>${T.toFixed(0)} K</td></tr>
    <tr><td>disc aspect ratio h/r</td><td>${hr.toFixed(3)}</td></tr>
    <tr><td>Stokes number of the largest grains (DSHARP VII's choice)</td><td>${f.k.toExponential(2)}</td></tr>
    <tr><td>planet / star mass ratio q</td><td>${q.toExponential(2)}</td></tr>
    <tr><td><b>planet mass</b></td><td><b>${mj.toFixed(3)} M<sub>J</sub> = ${me.toFixed(0)} M<sub>⊕</sub></b> (${like})</td></tr>
    <tr><td>our measured width of this gap; DSHARP VII's mass for it (α = 10⁻³, 0.1 mm)</td><td>Δ = ${g.delta.toFixed(3)}; ${g.pub} M<sub>J</sub></td></tr>`;
}
sel.addEventListener("change", loadDisc);
[inc, pa].forEach(e => e.addEventListener("input", update));
$("ours").onclick = () => setGeom(cur.ours); $("dsharp").onclick = () => setGeom(cur.dsharp); $("faceon").onclick = () => setGeom({ inc: 0, pa: cur.ours.pa });
["gap", "model", "delta", "alpha"].forEach(id => $(id).addEventListener("input", () => { if (id === "gap" && cur.gaps.length) $("delta").value = cur.gaps[+$("gap").value].delta; planet(); }));
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", update);
loadDisc();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    build(json.loads((FIGURE_DIR.parent / "data" / "results.json").read_text()))
    print(EXPLORER_PATH)
