"""An offline page for measuring a time delay yourself: figures/delay_finder.html.

Pick a lens and one of its images, then drag the slider: the image's light
curve slides in time over the reference image's. The page computes, at every
trial shift, how badly the two curves disagree (the dispersion of their
differences, Pelt et al. 1996); where they agree best is the delay. A button
reveals our measurement and the published one, and for the H0LiCOW lenses
the page turns the shift into the H0 that lens would give.

Both curves have the telescope zero points and the slow microlensing drift
removed (from our joint-spline fit), so what is left is the quasar's own
flicker. The page is self-contained: data embedded as JSON, drawn as SVG,
light and dark themes, a table view.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import figures as F
import lightcurves as L

#: Light / dark colour steps of each image's hue (the repository palette).
COLORS = {"A": ["#2a78d6", "#3987e5"], "B": ["#eb6834", "#d95926"], "C": ["#1baf7a", "#199e70"],
          "D": ["#4a3aa7", "#9085e9"]}


def _payload(keys=None) -> dict:
    p = C.RESULTS / "h0.json"
    res = json.loads(p.read_text()) if p.exists() else {"ours": {}}
    lenses = []
    for key in keys or C.MEASURED:
        lc, m, fit, comb = F.final_model(key)
        ref = comb["ref"]
        y = {im: F._corrected(lc, m, fit, im) for im in lc.images}
        base = np.median(y[ref])
        t0 = float(np.floor(lc.t.min()))
        imgs = []
        for im in comb["others"]:
            tau, err = comb["taus"][im], comb["err"][im]
            pv, pem, pep, src = F.published_delay(key, ref, im)
            span = max(30.0, 4 * err)
            imgs.append(dict(name=im, ours=round(tau, 2), err=round(err, 2), pub=round(pv, 2),
                             pub_err=[pem, pep], pub_src=src, lo=round(tau - span, 1), hi=round(tau + span, 1)))
        main = max(imgs, key=lambda d: abs(d["ours"]))
        h0 = res["ours"].get(key)
        lenses.append(dict(
            key=key, name=C.LENSES[key]["name"], ref=ref, t0=t0, year0=float(L.decimal_year(t0)),
            t=np.round(lc.t - t0, 3).tolist(),
            y={im: np.round(y[im] - base, 4).tolist() for im in lc.images},
            e={im: np.round(lc.err[im] * np.array([m.err_scale.get((im, tel), 1.0) for tel in lc.tel]), 4).tolist()
               for im in lc.images},
            images=imgs, main=main["name"],
            h0=None if h0 is None else [round(v, 1) for v in h0],
            notes=" ".join(lc.notes),
            source=lc.source))
    return dict(colors=COLORS, lenses=lenses)


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Delay Finder</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --ref: #0b0b0b; --band: rgba(42,120,214,0.16);
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
    --grid: #2f2f2d; --axis: #4a4a46; --ref: #ffffff; --band: rgba(57,135,229,0.22);
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
  --grid: #2f2f2d; --axis: #4a4a46; --ref: #ffffff; --band: rgba(57,135,229,0.22);
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1060px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 16px; max-width: 800px; }
.row { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 10px; align-items: center; }
.row .tag { font-size: 0.82rem; color: var(--muted); min-width: 52px; }
button { font: inherit; font-size: 0.9rem; padding: 6px 12px; border-radius: 8px; cursor: pointer;
         border: 1px solid var(--axis); background: var(--surface); color: var(--ink); }
button[aria-pressed="true"] { background: var(--ink); color: var(--surface); border-color: var(--ink); }
button:focus-visible, input:focus-visible { outline: 2px solid #2a78d6; outline-offset: 2px; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 300px; gap: 16px; }
@media (max-width: 820px) { .layout { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border-radius: 12px; padding: 12px; min-width: 0; }
svg { width: 100%; height: auto; display: block; }
.ctl { margin-bottom: 14px; }
.ctl label { display: flex; justify-content: space-between; font-size: 0.88rem; color: var(--ink2); gap: 8px; }
.ctl label b { color: var(--ink); font-variant-numeric: tabular-nums; }
.ctl input[type=range] { width: 100%; }
.readout { font-size: 0.88rem; color: var(--ink2); border-top: 1px solid var(--grid); padding-top: 10px; }
.readout div { margin-bottom: 7px; }
.readout b { color: var(--ink); font-variant-numeric: tabular-nums; }
.big { font-size: 1.5rem; color: var(--ink); line-height: 1.2; font-variant-numeric: tabular-nums; }
.note { font-size: 0.86rem; color: var(--ink2); margin: 10px 2px 0; }
.legend { display: flex; flex-wrap: wrap; gap: 14px; font-size: 0.82rem; color: var(--ink2); margin: 6px 2px 0; }
.legend span { display: inline-flex; align-items: center; gap: 6px; white-space: nowrap; }
details { margin-top: 16px; font-size: 0.85rem; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
td, th { padding: 2px 10px; text-align: right; border-bottom: 1px solid var(--grid); }
.tablewrap { max-height: 280px; overflow: auto; }
.axis text { fill: var(--muted); font-size: 11px; }
.gridline { stroke: var(--grid); }
.lbl { fill: var(--ink2); font-size: 12px; }
footer { color: var(--muted); font-size: 0.78rem; margin-top: 18px; }
#tip { position: fixed; pointer-events: none; background: var(--panel); color: var(--ink); border: 1px solid var(--axis);
       border-radius: 6px; padding: 4px 8px; font-size: 0.8rem; display: none; font-variant-numeric: tabular-nums; }
</style>
</head>
<body>
<main>
  <h1>Find the delay</h1>
  <p class="lede">A galaxy in front of a distant quasar bends its light into several images. Each image's light took a
  different path, so the same flicker arrives at different times. Slide one image's light curve in time over the
  reference image until the wiggles line up: that shift is the time delay. Delays like these, together with a model
  of the lens galaxy's mass, measure the expansion rate of the universe without any distance ladder.</p>
  <div class="row" id="lenses" role="group" aria-label="lens"><span class="tag">lens</span></div>
  <div class="row" id="imgs" role="group" aria-label="image to slide"><span class="tag">image</span></div>
  <div class="layout">
    <div class="card">
      <svg id="lc" viewBox="0 0 720 300" role="img" aria-label="reference and shifted light curves"></svg>
      <div class="legend" id="legend"></div>
      <div class="ctl" style="margin:10px 2px 0"><label for="win">Window <b id="wv"></b></label>
        <input id="win" type="range" min="0" max="1" step="0.001" value="0"></div>
      <div class="row"><span class="tag">show</span>
        <button id="z2" aria-pressed="true">2 years</button><button id="z5" aria-pressed="false">5 years</button>
        <button id="zall" aria-pressed="false">everything</button></div>
      <svg id="ds" viewBox="0 0 720 190" role="img" aria-label="mismatch against shift"></svg>
      <p class="note" id="note"></p>
    </div>
    <div class="card">
      <div class="ctl"><label for="shift"><span>Slide image <span id="imname"></span> back by</span><b id="sv"></b></label>
        <input id="shift" type="range" step="0.1"></div>
      <div class="row">
        <button id="reveal" aria-pressed="false">show the measurements</button>
        <button id="zero">no shift</button>
      </div>
      <div class="readout" id="readout"></div>
    </div>
  </div>
  <details><summary>Table view of the data</summary><div class="tablewrap"><table id="tbl"></table></div></details>
  <footer>Light curves: COSMOGRAIL (Millon et al. 2020a; Bonvin et al. 2018, 2019; Eulaers et al. 2013; Courbin et al.
  2018), from VizieR. Telescope zero points and microlensing (a slow spline per image, relative to the reference) were
  removed with this project's joint fit; magnitudes relative to the reference image's median. Mismatch = variance of
  the magnitude difference between the shifted image and the reference, interpolated between reference points less
  than 12 days apart. The H0 readout scales the H0 this lens gives with our delays (H0LiCOW lens models) by
  our delay / your shift, the dependence D<sub>&Delta;t</sub> &prop; delay at fixed lens model.</footer>
</main>
<div id="tip"></div>
<script>
const D = __DATA__;
const NS = "http://www.w3.org/2000/svg";
let li = 0, ii = 0, span = 730, revealed = false;
const dark = () => {
  const forced = document.documentElement.dataset.theme;
  if (forced) return forced === "dark";
  return window.matchMedia("(prefers-color-scheme: dark)").matches;
};
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const col = k => (D.colors[k[0]] || D.colors.A)[dark() ? 1 : 0];
function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function ticks(lo, hi, n) {
  const sp = hi - lo, step0 = sp / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 5, 10].map(m => m * mag).find(s => sp / s <= n) || 10 * mag;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}
const fmt = (v, d) => v.toLocaleString("en-US", { maximumFractionDigits: d, minimumFractionDigits: d });
// mismatch of image `im` shifted back by `s` days against the reference
function mismatch(lens, im, s) {
  const t = lens.t, yr = lens.y[lens.ref], yi = lens.y[im];
  let n = 0, m1 = 0, m2 = 0, j = 0;
  for (let k = 0; k < t.length; k++) {
    const x = t[k] - s;
    while (j < t.length - 2 && t[j + 1] < x) j++;
    if (x < t[j] || x > t[j + 1] || t[j + 1] - t[j] > 12) continue;
    const yref = yr[j] + (yr[j + 1] - yr[j]) * (x - t[j]) / (t[j + 1] - t[j]);
    const d = yi[k] - yref; n++; m1 += d; m2 += d * d;
  }
  if (n < 20) return [NaN, n];
  const mean = m1 / n; return [m2 / n - mean * mean, n];
}
let curve = [];
function buttons() {
  const box = document.getElementById("lenses");
  D.lenses.forEach((l, i) => {
    const b = document.createElement("button"); b.textContent = l.name; b.setAttribute("aria-pressed", i === li);
    b.onclick = () => { li = i; ii = 0; [...box.querySelectorAll("button")].forEach((c, j) => c.setAttribute("aria-pressed", j === li)); setupLens(); };
    box.appendChild(b);
  });
}
function setupLens() {
  const l = D.lenses[li], box = document.getElementById("imgs");
  [...box.querySelectorAll("button")].forEach(b => b.remove());
  const mainIdx = l.images.findIndex(x => x.name === l.main);
  ii = mainIdx >= 0 ? mainIdx : 0;
  l.images.forEach((im, i) => {
    const b = document.createElement("button"); b.textContent = `${im.name} vs ${l.ref}`; b.setAttribute("aria-pressed", i === ii);
    b.onclick = () => { ii = i; [...box.querySelectorAll("button")].forEach((c, j) => c.setAttribute("aria-pressed", j === ii)); setupImage(); };
    box.appendChild(b);
  });
  document.getElementById("note").textContent = l.notes ? "Note: " + l.notes + "." : "";
  table(); setupImage();
}
const slider = document.getElementById("shift"), win = document.getElementById("win");
function setupImage() {
  const l = D.lenses[li], im = l.images[ii];
  slider.min = im.lo; slider.max = im.hi; slider.step = 0.1;
  slider.value = Math.round((im.lo + 0.25 * (im.hi - im.lo)) * 10) / 10;
  document.getElementById("imname").textContent = im.name;
  curve = [];
  for (let s = im.lo; s <= im.hi + 1e-9; s += 0.25) curve.push([s, mismatch(l, im.name, s)[0]]);
  draw();
}
function draw() {
  const l = D.lenses[li], im = l.images[ii], sh = +slider.value;
  const ref = css("--ref"), c = col(im.name), bg = css("--surface");
  // ---- light curves ----
  const svg = document.getElementById("lc"); svg.innerHTML = "";
  const W = 720, H = 300, Lm = 46, R = 12, T = 12, B = 36;
  const tmax = l.t[l.t.length - 1];
  const sp = Math.min(span, tmax + 40);
  const x0 = -20 + (+win.value) * Math.max(0, tmax + 20 - sp), x1 = x0 + sp;
  document.getElementById("wv").textContent = (l.year0 + x0 / 365.25).toFixed(1) + "–" + (l.year0 + x1 / 365.25).toFixed(1);
  const vis = [];
  l.t.forEach((t, k) => { if (t >= x0 && t <= x1) vis.push(l.y[l.ref][k]); if (t - sh >= x0 && t - sh <= x1) vis.push(l.y[im.name][k]); });
  let ymin = Math.min(...vis), ymax = Math.max(...vis); const pad = 0.08 * (ymax - ymin || 0.2); ymin -= pad; ymax += pad;
  const X = x => Lm + (x - x0) / (x1 - x0) * (W - Lm - R), Y = y => T + (y - ymin) / (ymax - ymin) * (H - T - B); // brighter = up
  const ax = el("g", { class: "axis" }, svg);
  const yrs = ticks(l.year0 + x0 / 365.25, l.year0 + x1 / 365.25, 7);
  for (const v of yrs) { const xx = X((v - l.year0) * 365.25);
    el("line", { x1: xx, x2: xx, y1: T, y2: H - B, class: "gridline" }, ax);
    el("text", { x: xx, y: H - B + 16, "text-anchor": "middle" }, ax).textContent = v; }
  for (const v of ticks(ymin, ymax, 4)) {
    el("line", { x1: Lm, x2: W - R, y1: Y(v), y2: Y(v), class: "gridline" }, ax);
    el("text", { x: Lm - 6, y: Y(v) + 4, "text-anchor": "end" }, ax).textContent = fmt(v, 1); }
  el("text", { x: (Lm + W - R) / 2, y: H - 4, "text-anchor": "middle", class: "lbl" }, svg).textContent = "year (arrival time of image " + l.ref + ")";
  el("text", { x: 12, y: T + (H - T - B) / 2, class: "lbl", transform: `rotate(-90 12 ${T + (H - T - B) / 2})`, "text-anchor": "middle" }, svg).textContent = "Δmag (fainter ↓)";
  let d = ""; l.t.forEach((t, k) => { if (t < x0 - 15 || t > x1 + 15) return;
    const gap = k === 0 || t - l.t[k - 1] > 20 || d === "";
    d += (gap ? "M" : "L") + X(t).toFixed(1) + " " + Y(l.y[l.ref][k]).toFixed(1); });
  el("path", { d: d, fill: "none", stroke: ref, "stroke-width": 1.4, "stroke-opacity": 0.85 }, svg);
  const g = el("g", {}, svg);
  l.t.forEach((t, k) => { const ts = t - sh; if (ts < x0 || ts > x1) return;
    const p = el("circle", { cx: X(ts), cy: Y(l.y[im.name][k]), r: 3, fill: c, stroke: bg, "stroke-width": 0.8 }, g);
    p.dataset.tip = `image ${im.name}: observed ${(l.year0 + t / 365.25).toFixed(2)}, shown at ${(l.year0 + ts / 365.25).toFixed(2)}; Δmag ${fmt(l.y[im.name][k], 3)}`; });
  document.getElementById("legend").innerHTML =
    `<span><svg width="22" height="10"><line x1="0" y1="5" x2="22" y2="5" stroke="${ref}" stroke-width="2"/></svg>image ${l.ref} (reference)</span>` +
    `<span><svg width="12" height="12"><circle cx="6" cy="6" r="4" fill="${c}"/></svg>image ${im.name}, moved back by the slider</span>`;
  // ---- mismatch curve ----
  const ds = document.getElementById("ds"); ds.innerHTML = "";
  const H2 = 190, B2 = 34, s0 = im.lo, s1 = im.hi;
  const good = curve.filter(p => isFinite(p[1]));
  const vmin = Math.min(...good.map(p => p[1])), vmax = Math.max(...good.map(p => p[1]));
  const X2 = x => Lm + (x - s0) / (s1 - s0) * (W - Lm - R), Y2 = v => T + (1 - (v - vmin) / (vmax - vmin || 1)) * (H2 - T - B2);
  const ax2 = el("g", { class: "axis" }, ds);
  for (const v of ticks(s0, s1, 8)) { el("line", { x1: X2(v), x2: X2(v), y1: T, y2: H2 - B2, class: "gridline" }, ax2);
    el("text", { x: X2(v), y: H2 - B2 + 16, "text-anchor": "middle" }, ax2).textContent = v; }
  el("text", { x: (Lm + W - R) / 2, y: H2 - 4, "text-anchor": "middle", class: "lbl" }, ds).textContent = "shift (days)";
  el("text", { x: Lm + 4, y: T + 12, class: "lbl" }, ds).textContent = "mismatch (lower = better match)";
  if (revealed) {
    el("rect", { x: X2(im.ours - im.err), width: Math.max(1, X2(im.ours + im.err) - X2(im.ours - im.err)), y: T, height: H2 - T - B2, fill: css("--band") }, ds);
    el("line", { x1: X2(im.pub), x2: X2(im.pub), y1: T, y2: H2 - B2, stroke: css("--muted"), "stroke-width": 1.5, "stroke-dasharray": "4 3" }, ds);
  }
  let d2 = ""; curve.forEach(p => { if (!isFinite(p[1])) { d2 += ""; return; } d2 += (d2 ? "L" : "M") + X2(p[0]).toFixed(1) + " " + Y2(p[1]).toFixed(1); });
  el("path", { d: d2, fill: "none", stroke: css("--ink2"), "stroke-width": 1.8 }, ds);
  const now = mismatch(l, im.name, sh);
  el("line", { x1: X2(sh), x2: X2(sh), y1: T, y2: H2 - B2, stroke: c, "stroke-width": 1.4 }, ds);
  if (isFinite(now[0])) el("circle", { cx: X2(sh), cy: Y2(now[0]), r: 5, fill: c, stroke: bg, "stroke-width": 1.5 }, ds);
  // ---- readout ----
  document.getElementById("sv").textContent = fmt(sh, 1) + " d";
  const best = good.reduce((a, p) => p[1] < a[1] ? p : a, good[0]);
  let h = `<div>mismatch at this shift<br><span class="big">${isFinite(now[0]) ? fmt(1000 * Math.sqrt(now[0]), 1) + " mmag" : "too little overlap"}</span></div>`;
  h += `<div>${now[1]} points of image ${im.name} overlap the reference here</div>`;
  if (revealed) {
    h += `<div>lowest mismatch on the page: <b>${fmt(best[0], 1)} d</b></div>`;
    h += `<div>our measurement: <b>${fmt(im.ours, 1)} ± ${fmt(im.err, 1)} d</b> (blue band; joint spline fit with microlensing)</div>`;
    h += `<div>published: <b>${fmt(im.pub, 1)} d</b> (−${fmt(im.pub_err[0], 1)}/+${fmt(im.pub_err[1], 1)}; ${im.pub_src}; dashed line)</div>`;
  } else {
    h += `<div>find the lowest point of the mismatch curve, then press <i>show the measurements</i>.</div>`;
  }
  if (l.h0 && im.name === l.main && sh * im.ours > 0) {
    const H0 = l.h0[0] * im.ours / sh;
    h += `<div>with the published lens model, this shift would mean <span class="big">H0 ≈ ${fmt(H0, 0)}</span> km/s/Mpc` +
         ` (our delay gives ${fmt(l.h0[0], 1)} −${fmt(l.h0[1], 1)}/+${fmt(l.h0[2], 1)})</div>`;
  } else if (l.h0 && im.name !== l.main) {
    h += `<div>H0 readout: choose image ${l.main}, this lens's longest delay.</div>`;
  } else if (!l.h0) {
    h += `<div>(no public H0LiCOW distance posterior for this lens; TDCOSMO uses it with kinematics)</div>`;
  }
  document.getElementById("readout").innerHTML = h;
}
function table() {
  const l = D.lenses[li];
  const ims = Object.keys(l.y);
  let h = `<tr><th>year</th>${ims.map(k => `<th>${k} (Δmag)</th>`).join("")}</tr>`;
  l.t.forEach((t, k) => { h += `<tr><td>${(l.year0 + t / 365.25).toFixed(3)}</td>${ims.map(m => `<td>${fmt(l.y[m][k], 3)}</td>`).join("")}</tr>`; });
  document.getElementById("tbl").innerHTML = h;
}
const zbtn = { z2: 730, z5: 1826, zall: 1e6 };
for (const id in zbtn) document.getElementById(id).onclick = () => { span = zbtn[id];
  for (const j in zbtn) document.getElementById(j).setAttribute("aria-pressed", j === id); draw(); };
slider.addEventListener("input", draw); win.addEventListener("input", draw);
document.getElementById("reveal").onclick = e => { revealed = !revealed; e.target.setAttribute("aria-pressed", revealed); draw(); };
document.getElementById("zero").onclick = () => { slider.value = Math.max(+slider.min, Math.min(+slider.max, 0)); draw(); };
const tip = document.getElementById("tip");
document.getElementById("lc").addEventListener("mousemove", e => { const t = e.target.dataset && e.target.dataset.tip;
  if (t) { tip.textContent = t; tip.style.display = "block"; tip.style.left = (e.clientX + 12) + "px"; tip.style.top = (e.clientY + 12) + "px"; }
  else tip.style.display = "none"; });
document.getElementById("lc").addEventListener("mouseleave", () => tip.style.display = "none");
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);
buttons(); setupLens();
</script>
</body>
</html>
"""


def write(keys=None, name="delay_finder.html"):
    data = _payload(keys)
    html = PAGE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    out = C.FIG / name
    out.write_text(html)
    return out


if __name__ == "__main__":
    p = write()
    print(p, f"{p.stat().st_size / 1e3:.0f} kB")
