"""A self-contained HTML blink comparator: figures/remnant_explorer.html.

Two panels, no internet needed (images are embedded):

* **The Crab.** Blink the Gaia-registered images of 1951, 1991 and ~2012,
  overlay the 805 measured filament motions, or drag the "rewind" slider to
  run every filament backwards along its own motion and watch them converge
  on the expansion centre around 1106 -- while the pulsar walks back towards
  its birthplace.
* **Cas A in X-rays.** Step or play through all fifteen Chandra epochs,
  2000-2019, at a fixed stretch.
"""

from __future__ import annotations

import base64
import io
import json

import numpy as np
from PIL import Image

import config
import figures as F
import plotstyle as S

OUT = config.FIGURE_DIR / "remnant_explorer.html"


def _jpeg(arr01: np.ndarray, size: int, cmap=None, quality: int = 82) -> str:
    """0-1 array (row 0 = south) -> base64 JPEG, north up, east left."""
    a = np.clip(arr01[::-1], 0, 1)
    if cmap is not None:
        import matplotlib.pyplot as plt

        rgb = (plt.get_cmap(cmap)(a)[..., :3] * 255).astype(np.uint8)
        im = Image.fromarray(rgb)
    else:
        im = Image.fromarray((a * 255).astype(np.uint8)).convert("RGB")
    im = im.resize((size, size), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, format="JPEG", quality=quality)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def build(crab, crab_res: dict, vectors, frame, epochs, summ: dict) -> str:
    import casa_xray as X

    keys = ["crab_poss1_red", "crab_poss2_red", "crab_ps1_r"]
    disp = F._crab_display(crab, keys)
    crab_imgs = [{"label": f"{crab.regs[k].plate.label.split(' (')[0]}, {crab.regs[k].epoch:.1f}",
                  "src": _jpeg(disp[k], 540)} for k in keys]
    v = vectors
    vec = {"xi": np.round(v.xi, 2).tolist(), "eta": np.round(v.eta, 2).tolist(),
           "mx": np.round(v.mu_xi * 1000, 2).tolist(), "me": np.round(v.mu_eta * 1000, 2).tolist(),
           "tm": float(np.median(v.tmid))}
    m = crab_res["main"]
    p = crab_res["pulsar"]
    crab_meta = {"half": F.CRAB_HALF, "centre": m["centre"], "T": m["T"], "T_err": m["T_err_total"],
                 "pmra": p["pmra_mas_yr"], "pmdec": p["pmdec_mas_yr"], "martin": crab_res["martin2025_centre"]}

    cx, cy = frame.radec_to_box(X.CENTRE_RA, X.CENTRE_DEC)
    half = int(200 / config.CHANDRA_PIXEL_ARCSEC)
    sl = (slice(int(cy) - half, int(cy) + half), slice(int(cx) - half, int(cx) + half))
    fs = summ["regions"]["forward_shock"]["all"]
    rate = fs["rate_pct_yr"] / 100
    casa_imgs = []
    for ep in epochs:
        im = F._xray_norm(ep.counts)[sl]
        im = np.where(ep.footprint[sl], im, 0.0)
        casa_imgs.append({"label": f"{ep.year:.1f}", "src": _jpeg(im, 480, cmap=S.CMAP_XRAY),
                          "growth": 100 * rate * (ep.year - summ["t_ref"])})

    data = json.dumps({"crab": crab_imgs, "vec": vec, "crabMeta": crab_meta, "casa": casa_imgs,
                       "casaRate": fs["rate_pct_yr"]})
    html = TEMPLATE.replace("__DATA__", data)
    OUT.write_text(html)
    return str(OUT)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Rewind a Supernova Remnant</title>
<style>
:root { --bg:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --grid:#e1e0d9; --accent:#eb6834;
        --pulsar:#1baf7a; --pub:#4a3aa7; --panel:#f3f2ee; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#141413; --ink:#f0efe9; --ink2:#c3c2b7; --muted:#8f8d86; --grid:#2c2b28; --accent:#f07a4c;
  --pulsar:#2cc590; --pub:#8a7ee0; --panel:#1d1c1a; } }
:root[data-theme="dark"] { --bg:#141413; --ink:#f0efe9; --ink2:#c3c2b7; --muted:#8f8d86; --grid:#2c2b28;
  --accent:#f07a4c; --pulsar:#2cc590; --pub:#8a7ee0; --panel:#1d1c1a; }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--ink); font:15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1120px; margin: 0 auto; padding: 24px 16px 48px; }
h1 { font-size: 26px; margin: 0 0 6px; }
h2 { font-size: 19px; margin: 28px 0 6px; }
p.lede { color: var(--ink2); margin: 0 0 12px; max-width: 760px; }
.row { display:flex; gap:24px; flex-wrap:wrap; align-items:flex-start; }
.viewer { position:relative; width:min(540px, 100%); }
canvas { width:100%; height:auto; display:block; border-radius:6px; background:#000; }
.controls { flex:1 1 300px; min-width:260px; }
.btns { display:flex; gap:6px; flex-wrap:wrap; margin: 8px 0; }
button { font:inherit; font-size:14px; padding:6px 12px; border-radius:6px; border:1px solid var(--grid);
         background:var(--panel); color:var(--ink); cursor:pointer; }
button[aria-pressed="true"] { border-color: var(--accent); box-shadow: inset 0 0 0 1px var(--accent); }
label { display:block; margin: 10px 0 4px; color: var(--ink2); font-size:14px; }
input[type=range] { width:100%; accent-color: var(--accent); }
.readout { font-variant-numeric: tabular-nums; font-size: 22px; font-weight: 600; margin: 6px 0 2px; }
.small { color: var(--muted); font-size: 13px; }
.key span { display:inline-flex; align-items:center; gap:6px; margin-right:14px; font-size:13px; color:var(--ink2); }
.dot { width:10px; height:10px; border-radius:50%; display:inline-block; }
footer { margin-top: 32px; color: var(--muted); font-size: 12.5px; border-top: 1px solid var(--grid); padding-top: 10px; }
</style>
</head>
<body>
<main>
<h1>Rewind a supernova remnant</h1>
<p class="lede">Real images, registered on Gaia stars wound back to each date. Blink them to see the debris move,
then run every measured filament backwards to the explosion.</p>

<h2>The Crab Nebula, 1951&ndash;2012</h2>
<div class="row">
  <div class="viewer"><canvas id="crab" width="540" height="540" aria-label="Crab Nebula image viewer"></canvas></div>
  <div class="controls">
    <div class="btns" id="crabBtns" role="group" aria-label="Choose an image"></div>
    <div class="btns">
      <button id="blink" aria-pressed="false">Blink</button>
      <button id="vecs" aria-pressed="false">Show motions (x100 yr)</button>
    </div>
    <label for="rewind">Run the filaments back to the year:</label>
    <input id="rewind" type="range" min="1000" max="2012" step="1" value="2012">
    <div class="readout" id="year">2012</div>
    <div class="small" id="spread">Drag left to rewind. Each dot is one tracked filament moving back at its own measured speed.</div>
    <div class="key" style="margin-top:10px">
      <span><i class="dot" style="background:var(--accent)"></i>expansion centre (fit)</span>
      <span><i class="dot" style="background:var(--pulsar)"></i>Crab pulsar (Gaia DR3)</span>
      <span><i class="dot" style="background:var(--pub)"></i>Martin et al. 2025</span>
    </div>
    <p class="small" id="crabNote"></p>
  </div>
</div>

<h2>Cassiopeia A in X-rays, 2000&ndash;2019</h2>
<div class="row">
  <div class="viewer" style="width:min(480px,100%)"><canvas id="casa" width="480" height="480" aria-label="Cas A Chandra viewer"></canvas></div>
  <div class="controls">
    <label for="ep">Chandra epoch:</label>
    <input id="ep" type="range" min="0" max="0" step="1" value="0">
    <div class="readout" id="casaYear"></div>
    <div class="small" id="casaGrowth"></div>
    <div class="btns"><button id="play" aria-pressed="false">Play</button></div>
    <p class="small">Fifteen Chandra ACIS-S images at the same stretch (all energies). Watch the thin outer rim, the
    forward shock, creep outward; the bright ring of shocked debris inside changes as knots light up and fade.</p>
  </div>
</div>

<footer>Data: Digitized Sky Survey (POSS-I 1951, POSS-II 1991; STScI), Pan-STARRS1 (2009&ndash;2014), Gaia DR3,
Chandra Data Archive. Analysis and page: this atlas, project L3 (AI-generated; see the project README for methods and caveats).</footer>
</main>
<script>
const D = __DATA__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

// ---------- Crab ----------
const cv = document.getElementById('crab'), ctx = cv.getContext('2d');
const imgs = D.crab.map(o => { const im = new Image(); im.src = o.src; return im; });
let cur = 2, blinkTimer = null, showVec = false, year = 2012;
const H = D.crabMeta.half, S = cv.width / (2 * H);
const toC = (xi, eta) => [(H - xi) * S, (H - eta) * S];   // east to the left, north up
const btnBox = document.getElementById('crabBtns');
D.crab.forEach((o, i) => { const b = document.createElement('button'); b.textContent = o.label;
  b.onclick = () => { stopBlink(); cur = i; setYear(2012); draw(); }; btnBox.appendChild(b); });
function pulsarAt(y) { return [D.crabMeta.pmra * (y - 2016) / 1000, D.crabMeta.pmdec * (y - 2016) / 1000]; }
function cross(x, y, col, r) { ctx.strokeStyle = col; ctx.lineWidth = 2.5; ctx.beginPath();
  ctx.moveTo(x - r, y); ctx.lineTo(x + r, y); ctx.moveTo(x, y - r); ctx.lineTo(x, y + r); ctx.stroke(); }
function draw() {
  [...btnBox.children].forEach((b, i) => b.setAttribute('aria-pressed', String(i === cur && year === 2012)));
  ctx.clearRect(0, 0, cv.width, cv.height);
  const rewinding = year < 2012;
  ctx.globalAlpha = rewinding ? 0.35 : 1; ctx.drawImage(imgs[cur], 0, 0, cv.width, cv.height); ctx.globalAlpha = 1;
  const v = D.vec, n = v.xi.length;
  if (showVec && !rewinding) {
    ctx.strokeStyle = '#ffb347'; ctx.lineWidth = 1;
    for (let i = 0; i < n; i += 2) { const [x0, y0] = toC(v.xi[i], v.eta[i]);
      const [x1, y1] = toC(v.xi[i] + v.mx[i] * 0.1, v.eta[i] + v.me[i] * 0.1);
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke(); }
  }
  if (rewinding) {
    ctx.fillStyle = '#ffd27a'; let sx = 0, sy = 0, xs = [], ys = [];
    for (let i = 0; i < n; i++) { const dt = (v.tm - year) / 1000;
      const xi = v.xi[i] - v.mx[i] * dt, eta = v.eta[i] - v.me[i] * dt; xs.push(xi); ys.push(eta);
      const [x, y] = toC(xi, eta); ctx.fillRect(x - 1, y - 1, 2.2, 2.2); }
    const mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
    const rms = Math.sqrt(xs.reduce((a, x, i) => a + (x - mx) ** 2 + (ys[i] - my) ** 2, 0) / n);
    document.getElementById('spread').textContent = `Filaments within ${rms.toFixed(0)}" (rms) of each other.` +
      (Math.abs(year - D.crabMeta.T) < 15 ? ' This is where the fit puts the explosion.' : year < 1080 ? ' Too far: they have passed through each other.' : '');
  }
  const [px, py] = toC(...pulsarAt(rewinding ? year : 2016));
  ctx.fillStyle = css('--pulsar'); ctx.beginPath(); ctx.arc(px, py, 4.5, 0, 7); ctx.fill();
  const [cx, cy] = toC(...D.crabMeta.centre); cross(cx, cy, css('--accent'), 8);
  const [mx2, my2] = toC(...D.crabMeta.martin); ctx.strokeStyle = css('--pub'); ctx.lineWidth = 2;
  ctx.beginPath(); ctx.arc(mx2, my2, 4, 0, 7); ctx.stroke();
  ctx.fillStyle = '#fff'; ctx.font = '600 18px system-ui'; ctx.fillText(rewinding ? String(year) : D.crab[cur].label, 12, 26);
}
function stopBlink() { if (blinkTimer) clearInterval(blinkTimer); blinkTimer = null;
  document.getElementById('blink').setAttribute('aria-pressed', 'false'); }
document.getElementById('blink').onclick = e => {
  if (blinkTimer) { stopBlink(); return; }
  setYear(2012); e.target.setAttribute('aria-pressed', 'true');
  blinkTimer = setInterval(() => { cur = (cur + 1) % imgs.length; draw(); }, 800); };
document.getElementById('vecs').onclick = e => { showVec = !showVec; e.target.setAttribute('aria-pressed', String(showVec)); draw(); };
const rw = document.getElementById('rewind');
function setYear(y) { year = y; rw.value = y; document.getElementById('year').textContent = y; draw(); }
rw.oninput = () => { stopBlink(); setYear(+rw.value); };
document.getElementById('crabNote').textContent =
  `Best fit: the filaments left the centre in ${D.crabMeta.T.toFixed(0)} ± ${D.crabMeta.T_err.toFixed(0)} CE, ` +
  `about 50 years after the guest star of 1054: they have been sped up since, by the pulsar's wind.`;
imgs.forEach(im => im.onload = draw);

// ---------- Cas A ----------
const cc = document.getElementById('casa'), c2 = cc.getContext('2d');
const cimgs = D.casa.map(o => { const im = new Image(); im.src = o.src; return im; });
const ep = document.getElementById('ep'); ep.max = cimgs.length - 1;
let k = 0, playTimer = null;
function drawCasa() { c2.drawImage(cimgs[k], 0, 0, cc.width, cc.height);
  c2.fillStyle = '#fff'; c2.font = '600 18px system-ui'; c2.fillText(D.casa[k].label, 12, 26);
  document.getElementById('casaYear').textContent = D.casa[k].label;
  const g = D.casa[k].growth;
  document.getElementById('casaGrowth').textContent =
    `Forward shock ${g >= 0 ? 'larger' : 'smaller'} than in 2004 by ${Math.abs(g).toFixed(2)}% (fitted rate ${D.casaRate.toFixed(3)}% per year).`; }
ep.oninput = () => { k = +ep.value; drawCasa(); };
document.getElementById('play').onclick = e => {
  if (playTimer) { clearInterval(playTimer); playTimer = null; e.target.setAttribute('aria-pressed', 'false'); return; }
  e.target.setAttribute('aria-pressed', 'true');
  playTimer = setInterval(() => { k = (k + 1) % cimgs.length; ep.value = k; drawCasa(); }, 450); };
cimgs[0].onload = drawCasa;
</script>
</body>
</html>
"""
