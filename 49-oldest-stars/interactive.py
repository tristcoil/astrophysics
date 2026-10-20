"""An offline page to fit the clusters yourself: figures/isochrone_explorer.html.

Pick a cluster, then move the age, distance, reddening, [Fe/H] and colour-shift
sliders: the MIST isochrone ([alpha/Fe] = +0.4) moves over the HUGS stars, and a
small panel shows chi^2 against age for the other settings, as the real fit
computes it (in a simplified form: distance from each ridge point to the nearest
point of the isochrone, with the same error floors). Everything is embedded; no
network access is needed.
"""

from __future__ import annotations

import json

import numpy as np

import agefit as A
import config as C
from fiducial import K
from isochrones import Extinction, load_grid

AGES = np.round(np.arange(8.0, 16.01, 0.5), 2)
FEH_STEPS = np.round(np.arange(-0.3, 0.301, 0.1), 2)
N_STARS = 7000


def _cluster_payload(name: str, results: dict, art: dict, grid) -> dict:
    prep, inp = art["prepared"][name], art["inputs"][name]
    res = results["clusters"][name]
    r = prep["ridge"]
    win_lo, win_hi = r.m_to + 3.6, r.m_to - 4.5  # faint, bright limits (apparent F606W)
    rng = np.random.default_rng(1)
    sel = np.flatnonzero((prep["m606"] < win_lo) & (prep["m606"] > win_hi)
                         & (prep["c"] > r.c_to - 0.4) & (prep["c"] < r.c_to + 0.9))
    if len(sel) > N_STARS:
        sel = rng.choice(sel, N_STARS, replace=False)
    ext = Extinction(inp["feh"], 0.4)
    fehs = np.round(inp["feh"] + FEH_STEPS, 3)
    M, Cc, k606, k814 = [], [], [], []
    eep_keep = None
    blocks = []
    for feh in fehs:
        iso = grid.interpolate(float(feh), 0.4, AGES)
        blocks.append(iso)
        ok = np.isfinite(iso["F606W"]) & (iso["F606W"] < 7.5) & (iso["F606W"] > -1.5)
        cols = np.any(ok, axis=0)
        eep_keep = cols if eep_keep is None else (eep_keep | cols)
    for iso in blocks:
        m6 = iso["F606W"][:, eep_keep]
        m8 = iso["F814W"][:, eep_keep]
        lt, lg = iso["logTe"][:, eep_keep], iso["logg"][:, eep_keep]
        good = np.isfinite(m6) & np.isfinite(m8)
        av = C.R_V * max(inp["ebv"], 0.02)
        a6 = np.where(good, ext("F606W", np.where(good, lt, 3.8), np.where(good, lg, 4.0), av), np.nan)
        a8 = np.where(good, ext("F814W", np.where(good, lt, 3.8), np.where(good, lg, 4.0), av), np.nan)
        M.append(np.where(good, m6, np.nan))
        Cc.append(np.where(good, m6 - m8, np.nan))
        k606.append(a6 / max(inp["ebv"], 0.02))
        k814.append(a8 / max(inp["ebv"], 0.02))

    def pack(arr, nd=4):
        return [[[None if not np.isfinite(v) else round(float(v), nd) for v in row] for row in blk] for blk in arr]

    mist = res["models"]["MIST"]
    fit = A.ClusterFitter(grid, r, feh0=inp["feh"])
    return {
        "name": name,
        "stars": {"c": [round(float(v), 3) for v in prep["c"][sel]], "m": [round(float(v), 3) for v in prep["m606"][sel]]},
        "ridge": {"c": [round(float(v), 4) for v in fit.ridge.c], "m": [round(float(v), 4) for v in fit.ridge.m],
                  "sig": [round(float(v), 5) for v in np.sqrt(fit.sig2)],
                  "nx": [round(float(v), 4) for v in fit.ridge.nx], "ny": [round(float(v), 4) for v in fit.ridge.ny]},
        "feh_grid": [float(v) for v in fehs], "ages": [float(a) for a in AGES],
        "M": pack(M), "C": pack(Cc), "k606": pack(k606, 3), "k814": pack(k814, 3),
        "inputs": {"feh": inp["feh"], "dm": round(inp["dm"], 3), "dm_err": round(inp["dm_err"], 3), "ebv": inp["ebv"],
                   "dist_kpc": inp["dist_kpc"], "dm_geo": round(inp.get("dm_geo", inp["dm"]), 3)},
        "fit": {"age": round(mist["nominal_age"], 2), "dc": round(mist["nominal_dc"], 3),
                "combined": round(res["combined"]["age"], 2), "combined_err": round(res["combined"]["total_err"], 2)},
        "view": {"c0": round(r.c_to - 0.12, 3), "c1": round(r.c_to + 0.55, 3), "m0": round(r.m_to + 3.2, 2),
                 "m1": round(r.m_to - 4.2, 2)},
        "presets": _presets(name, inp, mist),
    }


def _presets(name, inp, mist):
    out = [{"label": "our MIST fit", "age": mist["nominal_age"], "dm": inp["dm"], "ebv": inp["ebv"],
            "feh": inp["feh"], "dc": mist["nominal_dc"]}]
    if "dm_geo" in inp:
        out.append({"label": "geometric distance", "dm": inp["dm_geo"]})
    if name == "M92":
        out.append({"label": "Ying et al. 2023 distance", "dm": C.YING2023_M92["dm"], "ebv": 0.01})
    return out


def build(results: dict, art: dict) -> str:
    grid = load_grid("MIST")
    payload = {name: _cluster_payload(name, results, art, grid) for name in C.CORE_CLUSTERS}
    payload["_meta"] = {"K": K, "floor_c": A.FLOOR_COLOUR, "floor_m": A.FLOOR_MAG,
                        "universe": results["cosmology"]["universe"].get("D3", {}).get("age", 13.8)}
    html = TEMPLATE.replace("__DATA__", json.dumps(payload, separators=(",", ":")))
    out = C.FIGURE_DIR / "isochrone_explorer.html"
    out.write_text(html)
    return str(out)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cluster Age Explorer</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --page: #f9f9f7; --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --stars: rgba(11,11,11,0.28); --ridge: #0b0b0b;
  --iso: #2a78d6; --band: #b7d3f6; --mark: #eb6834; --ring: rgba(11,11,11,0.10);
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --stars: rgba(255,255,255,0.30); --ridge: #ffffff;
    --iso: #3987e5; --band: #184f95; --mark: #d95926; --ring: rgba(255,255,255,0.10);
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --stars: rgba(255,255,255,0.30); --ridge: #ffffff;
  --iso: #3987e5; --band: #184f95; --mark: #d95926; --ring: rgba(255,255,255,0.10);
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1100px; margin: 0 auto; padding: 24px 16px 48px; }
h1 { font-size: 22px; margin: 0 0 4px; }
p.lede { color: var(--ink-2); margin: 0 0 18px; max-width: 760px; }
.tabs { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 14px; }
.tabs button { font: inherit; padding: 6px 12px; border-radius: 8px; border: 1px solid var(--axis);
               background: var(--surface); color: var(--ink); cursor: pointer; }
.tabs button[aria-pressed="true"] { border-color: var(--iso); box-shadow: 0 0 0 1px var(--iso) inset; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 320px; gap: 16px; }
@media (max-width: 820px) { .layout { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: 12px; padding: 14px; }
canvas { width: 100%; display: block; }
#cmd { aspect-ratio: 4 / 5; }
#chi { aspect-ratio: 16 / 7; margin-top: 8px; }
label { display: block; font-size: 13px; color: var(--ink-2); margin-top: 12px; }
label b { color: var(--ink); font-weight: 600; float: right; font-variant-numeric: tabular-nums; }
label.check { display: flex; gap: 8px; align-items: center; margin-top: 6px; }
input[type=range] { width: 100%; accent-color: var(--iso); }
.presets { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 14px; }
.presets button { font: inherit; font-size: 13px; padding: 4px 10px; border-radius: 6px;
                  border: 1px solid var(--axis); background: transparent; color: var(--ink); cursor: pointer; }
.readout { margin-top: 10px; font-size: 14px; }
.readout .big { font-size: 26px; font-weight: 600; }
.note { font-size: 13px; color: var(--ink-2); margin-top: 12px; }
.theme { float: right; font: inherit; font-size: 13px; padding: 4px 10px; border-radius: 6px;
         border: 1px solid var(--axis); background: var(--surface); color: var(--ink); cursor: pointer; }
</style>
</head>
<body>
<main>
<button class="theme" id="theme" type="button">Dark / light</button>
<h1>How old is this cluster? Fit it yourself</h1>
<p class="lede">Real Hubble photometry (HUGS) of four globular clusters with a MIST v2.5 isochrone on top.
Slide the age until the line runs through the turnoff and subgiant branch. Then change the distance and
watch the best age move: a cluster that is farther away has a brighter real turnoff, so it must be younger.</p>
<div class="tabs" id="tabs" role="group" aria-label="Cluster"></div>
<div class="layout">
  <div class="card">
    <canvas id="cmd" role="img" aria-label="Colour-magnitude diagram with isochrone"></canvas>
  </div>
  <div class="card">
    <div class="readout">age <span class="big" id="ageBig"></span> Gyr
      <div id="chiText" class="note"></div></div>
    <label for="age">age <b id="ageV"></b></label>
    <input type="range" id="age" min="8" max="16" step="0.05">
    <label for="dm">distance modulus (m − M)₀ <b id="dmV"></b></label>
    <input type="range" id="dm" step="0.005">
    <label for="ebv">reddening E(B−V) <b id="ebvV"></b></label>
    <input type="range" id="ebv" min="0" max="0.3" step="0.005">
    <label for="feh">[Fe/H] <b id="fehV"></b></label>
    <input type="range" id="feh" step="0.01">
    <label for="dc">colour shift Δc <b id="dcV"></b></label>
    <input type="range" id="dc" min="-0.08" max="0.08" step="0.002">
    <label class="check"><input type="checkbox" id="auto" checked> fit the colour shift automatically, as the real fit does</label>
    <div class="presets" id="presets"></div>
    <canvas id="chi" role="img" aria-label="chi-squared against age"></canvas>
    <div class="note">χ² of the ridge line (the dots) against isochrones of every age, other settings
      as chosen (with the colour shift re-fitted at each age when the box is ticked); the dot marks the best
      age. The project's real fit also varies the inputs 1,000 times and uses three isochrone sets.</div>
  </div>
</div>
<p class="note" id="foot"></p>
</main>
<script id="data" type="application/json">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById('data').textContent);
const META = D._meta; delete D._meta;
const names = Object.keys(D);
let cur = names[0];
const st = {};
const $ = id => document.getElementById(id);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

function bracket(arr, x) {
  let i = 0;
  while (i < arr.length - 2 && x > arr[i + 1]) i++;
  const t = Math.min(Math.max((x - arr[i]) / (arr[i + 1] - arr[i]), 0), 1);
  return [i, t];
}
function iso(cl, age, feh, dm, ebv, dc) {
  const [fi, ft] = bracket(cl.feh_grid, feh);
  const la = cl.ages.map(Math.log10);
  const [ai, at] = bracket(la, Math.log10(age));
  const n = cl.M[0][0].length, xs = [], ys = [];
  for (let k = 0; k < n; k++) {
    let M = 0, Cc = 0, k6 = 0, k8 = 0, ok = true;
    for (const [f, wf] of [[fi, 1 - ft], [fi + 1, ft]]) for (const [a, wa] of [[ai, 1 - at], [ai + 1, at]]) {
      const w = wf * wa; if (w === 0) continue;
      const m = cl.M[f][a][k], c = cl.C[f][a][k];
      if (m === null || c === null) { ok = false; break; }
      M += w * m; Cc += w * c; k6 += w * cl.k606[f][a][k]; k8 += w * cl.k814[f][a][k];
    }
    if (!ok) { xs.push(NaN); ys.push(NaN); continue; }
    ys.push(M + dm + k6 * ebv); xs.push(Cc + (k6 - k8) * ebv + dc);
  }
  return [xs, ys];
}
function chi2(cl, xs, ys, dc) {
  const K = META.K, r = cl.ridge; let s = 0; dc = dc || 0;
  for (let i = 0; i < r.c.length; i++) {
    const px = K * r.c[i], py = r.m[i]; let best = Infinity;
    for (let k = 0; k + 1 < xs.length; k++) {
      const ax = K * (xs[k] + dc), ay = ys[k], bx = K * (xs[k + 1] + dc) - ax, by = ys[k + 1] - ay;
      if (!isFinite(ax) || !isFinite(bx)) continue;
      const L = bx * bx + by * by; if (L === 0) continue;
      let u = ((px - ax) * bx + (py - ay) * by) / L; u = Math.max(0, Math.min(1, u));
      const dx = px - ax - u * bx, dy = py - ay - u * by, d = dx * dx + dy * dy;
      if (d < best) best = d;
    }
    s += Math.min(best, 0.09) / (r.sig[i] * r.sig[i]);
  }
  return s;
}
function bestDc(cl, xs, ys) {
  // the real fit lets one colour shift float at every age: golden-section search
  let a = -0.08, b = 0.08; const g = 0.6180339887;
  let c = b - g * (b - a), d = a + g * (b - a), fc = chi2(cl, xs, ys, c), fd = chi2(cl, xs, ys, d);
  for (let k = 0; k < 16; k++) {
    if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = chi2(cl, xs, ys, c); }
    else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = chi2(cl, xs, ys, d); }
  }
  return (a + b) / 2;
}
function setupCanvas(cv) {
  const dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); return [g, w, h];
}
function drawCMD(cl, xs, ys) {
  const [g, w, h] = setupCanvas($('cmd'));
  const v = cl.view, L = 58, R = 10, T = 10, B = 42;
  const X = c => L + (c - v.c0) / (v.c1 - v.c0) * (w - L - R);
  const Y = m => T + (m - v.m1) / (v.m0 - v.m1) * (h - T - B);
  g.fillStyle = css('--surface'); g.fillRect(0, 0, w, h);
  g.strokeStyle = css('--grid'); g.lineWidth = 1; g.fillStyle = css('--muted'); g.font = '12px system-ui';
  g.textAlign = 'center';
  for (let c = Math.ceil(v.c0 * 10) / 10; c <= v.c1; c += 0.1) {
    g.beginPath(); g.moveTo(X(c), T); g.lineTo(X(c), h - B); g.stroke(); g.fillText(c.toFixed(1), X(c), h - B + 15);
  }
  g.textAlign = 'right';
  for (let m = Math.ceil(v.m1); m <= v.m0; m += 1) {
    g.beginPath(); g.moveTo(L, Y(m)); g.lineTo(w - R, Y(m)); g.stroke(); g.fillText(m.toFixed(0), L - 6, Y(m) + 4);
  }
  g.textAlign = 'center';
  g.fillText('colour  F606W − F814W', L + (w - L - R) / 2, h - 6);
  g.save(); g.translate(13, T + (h - T - B) / 2); g.rotate(-Math.PI / 2); g.fillText('F606W (mag, fainter ↓)', 0, 0); g.restore();
  g.textAlign = 'left';
  g.fillStyle = css('--stars');
  const s = cl.stars;
  for (let i = 0; i < s.c.length; i++) g.fillRect(X(s.c[i]) - 0.75, Y(s.m[i]) - 0.75, 1.5, 1.5);
  g.fillStyle = css('--ridge');
  for (let i = 0; i < cl.ridge.c.length; i++) { g.beginPath(); g.arc(X(cl.ridge.c[i]), Y(cl.ridge.m[i]), 2, 0, 6.3); g.fill(); }
  g.strokeStyle = css('--iso'); g.lineWidth = 2; g.beginPath(); let pen = false;
  for (let k = 0; k < xs.length; k++) {
    if (!isFinite(xs[k])) { pen = false; continue; }
    if (!pen) { g.moveTo(X(xs[k]), Y(ys[k])); pen = true; } else g.lineTo(X(xs[k]), Y(ys[k]));
  }
  g.stroke();
}
function drawChi(cl) {
  const [g, w, h] = setupCanvas($('chi'));
  const ages = [], vals = [];
  for (let a = 8; a <= 16.001; a += 0.2) {
    const [xs, ys] = iso(cl, a, st.feh, st.dm, st.ebv, 0);
    const dc = st.auto ? bestDc(cl, xs, ys) : st.dc;
    ages.push(a); vals.push(chi2(cl, xs, ys, dc));
  }
  let k = 0; vals.forEach((v, i) => { if (v < vals[k]) k = i; });
  const lo = Math.min(...vals), hi = Math.min(Math.max(...vals), lo * 40 + 50);
  const L = 34, R = 8, T = 8, B = 22;
  const X = a => L + (a - 8) / 8 * (w - L - R), Y = v => T + (1 - (Math.min(v, hi) - lo) / (hi - lo)) * (h - T - B);
  g.fillStyle = css('--surface'); g.fillRect(0, 0, w, h);
  g.strokeStyle = css('--axis'); g.beginPath(); g.moveTo(L, h - B); g.lineTo(w - R, h - B); g.stroke();
  g.fillStyle = css('--muted'); g.font = '11px system-ui';
  for (let a = 8; a <= 16; a += 2) g.fillText(a + '', X(a) - 4, h - 6);
  const u = META.universe; g.strokeStyle = css('--ink-2'); g.beginPath(); g.moveTo(X(u), T); g.lineTo(X(u), h - B); g.stroke();
  g.fillText('universe', X(u) + 3, T + 10);
  g.strokeStyle = css('--iso'); g.lineWidth = 2; g.beginPath();
  ages.forEach((a, i) => i ? g.lineTo(X(a), Y(vals[i])) : g.moveTo(X(a), Y(vals[i]))); g.stroke();
  g.fillStyle = css('--mark'); g.beginPath(); g.arc(X(ages[k]), Y(vals[k]), 4, 0, 6.3); g.fill();
  g.strokeStyle = css('--ink'); g.beginPath(); g.moveTo(X(st.age), T); g.lineTo(X(st.age), h - B); g.stroke();
  return ages[k];
}
function render() {
  const cl = D[cur];
  for (const k of ['age', 'dm', 'ebv', 'feh', 'dc']) $(k).value = st[k];
  $('ageV').textContent = st.age.toFixed(2) + ' Gyr'; $('ageBig').textContent = st.age.toFixed(1);
  const dkpc = Math.pow(10, st.dm / 5 + 1) / 1000;
  $('dmV').textContent = st.dm.toFixed(3) + '  (' + dkpc.toFixed(2) + ' kpc)';
  $('ebvV').textContent = st.ebv.toFixed(3); $('fehV').textContent = st.feh.toFixed(2);
  let [xs, ys] = iso(cl, st.age, st.feh, st.dm, st.ebv, 0);
  if (st.auto) { st.dc = bestDc(cl, xs, ys); $('dc').value = st.dc; }
  $('dc').disabled = st.auto;
  $('dcV').textContent = (st.dc >= 0 ? '+' : '') + st.dc.toFixed(3) + (st.auto ? ' (fitted)' : '');
  const xd = xs.map(x => x + st.dc);
  drawCMD(cl, xd, ys);
  const best = drawChi(cl);
  $('chiText').textContent = 'χ² = ' + chi2(cl, xs, ys, st.dc).toFixed(1) + ' for ' + cl.ridge.c.length +
    ' ridge points · best age for these settings ≈ ' + best.toFixed(1) + ' Gyr';
  $('foot').textContent = cl.name + ': distance ' + cl.inputs.dist_kpc.toFixed(2) + ' kpc (Baumgardt & Vasiliev 2021), ' +
    'E(B−V) = ' + cl.inputs.ebv + ', [Fe/H] = ' + cl.inputs.feh + ' (Harris 2010). Project result (three isochrone sets): ' +
    cl.fit.combined + ' ± ' + cl.fit.combined_err + ' Gyr.';
}
function choose(name) {
  cur = name; const cl = D[name];
  document.querySelectorAll('#tabs button').forEach(b => b.setAttribute('aria-pressed', b.dataset.name === name));
  Object.assign(st, {age: cl.fit.age, dm: cl.inputs.dm, ebv: cl.inputs.ebv, feh: cl.inputs.feh, dc: cl.fit.dc});
  $('auto').checked = st.auto;
  $('dm').min = (cl.inputs.dm - 0.4).toFixed(3); $('dm').max = (cl.inputs.dm + 0.4).toFixed(3);
  $('feh').min = cl.feh_grid[0]; $('feh').max = cl.feh_grid[cl.feh_grid.length - 1];
  const p = $('presets'); p.innerHTML = '';
  for (const pr of cl.presets) {
    const b = document.createElement('button'); b.type = 'button'; b.textContent = pr.label;
    b.onclick = () => { for (const k of ['age', 'dm', 'ebv', 'feh', 'dc']) if (k in pr) st[k] = pr[k]; schedule(); };
    p.appendChild(b);
  }
  render();
}
for (const n of names) {
  const b = document.createElement('button'); b.type = 'button'; b.textContent = n; b.dataset.name = n;
  b.onclick = () => choose(n); $('tabs').appendChild(b);
}
let pending = false;
function schedule() { if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; render(); }); } }
for (const k of ['age', 'dm', 'ebv', 'feh', 'dc'])
  $(k).addEventListener('input', e => { st[k] = parseFloat(e.target.value); schedule(); });
st.auto = true;
$('auto').addEventListener('change', e => { st.auto = e.target.checked; schedule(); });
$('theme').onclick = () => {
  const root = document.documentElement;
  const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
  root.dataset.theme = dark ? 'light' : 'dark'; render();
};
window.addEventListener('resize', render);
choose(cur);
</script>
</body>
</html>
"""
