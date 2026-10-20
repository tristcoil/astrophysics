"""An offline page for fitting the M-sigma relation yourself: figures/m_sigma_explorer.html.

Pick a compilation (van den Bosch 2016, Kormendy & Ho 2013, McConnell & Ma
2013), switch kinds of galaxies and measurements in or out, and the slope is
refitted in the browser with the same likelihood as the project's headline fit
(Gaussian scatter in log M, errors in both axes, the distribution of true sigma
modelled as one Gaussian, upper limits censored), and Laplace error bars.

The page is self-contained (no network, no libraries): data embedded as JSON,
drawn as SVG, light and dark themes, a table view of the selection.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import config as C
import literature as LIT


def _num(v, nd=3):
    return None if v is None or not np.isfinite(v) else round(float(v), nd)


def _rows(df: pd.DataFrame, kh_info: dict, src: str) -> list[dict]:
    out = []
    for _, r in df.iterrows():
        if not np.isfinite(r["logsig"]) or not np.isfinite(r["logM"]):
            continue
        info = kh_info.get(r["name"], {"kh": "unknown", "aside": False})
        if src == "kh":
            info = {"kh": {"elliptical": "elliptical", "classical": "classical", "pseudo": "pseudo"}[r["bulge"]],
                    "aside": bool(r["bulge"] in ("elliptical", "classical") and not r["kh13_fit"])}
        upper = bool(r.get("upper", False))
        lim = r.get("logM_limit", np.nan) if src == "vdb" else r["logM"]
        out.append({
            "n": r["name"], "s": _num(r["logsig"]), "es": _num(r.get("elogsig", 0.0)) or 0.0,
            "m": _num(r["logM"]), "em": _num(r.get("elogM", np.nan)) if not (src == "kh" and upper) else 0.0,
            "u": int(upper), "l": _num(lim) if upper else None, "t": r["method"], "kh": info["kh"],
            "x": int(info["aside"]), "d": _num(r.get("dist", np.nan), 2),
            "w": _num(r.get("weight", 1.0), 2) or 1.0,
            "r": str(r.get("ref", r.get("source_ref", r.get("mass_ref", "")))),
        })
    return out


def _payload(res: dict) -> dict:
    cat, fits = res["cat"], res["fits"]
    kh = cat["kh13"]
    kh_info = {}
    for _, r in kh.drop_duplicates("name", keep="first").iterrows():
        kh_info[r["name"]] = {"kh": r["bulge"], "aside": bool(r["bulge"] in ("elliptical", "classical")
                                                              and not r["kh13_fit"])}
    v = cat["vdb16"]
    vs = v[v["selected"] == 1]
    mm = cat["mm13"].copy()
    mm["ref"] = "McConnell & Ma 2013, ref. " + mm["mass_ref"].astype(str)
    samples = {
        "vdb": {"label": "van den Bosch 2016", "rows": _rows(vs, kh_info, "vdb"),
                "note": "230 galaxies in van den Bosch's fit; upper limits are the compilation's own (centre and "
                        "half-range, as he fitted them)."},
        "kh": {"label": "Kormendy & Ho 2013", "rows": _rows(kh, kh_info, "kh"),
               "note": "Their Tables 2 and 3, classifications included. Their own fit uses only ellipticals and "
                       "classical bulges, and sets aside mergers, monsters, gas masses without line widths and "
                       "the two largest black holes."},
        "mm": {"label": "McConnell & Ma 2013", "rows": _rows(mm, kh_info, "mm"),
               "note": "Their 72 galaxies (NGC 1399 twice, each at half weight). Bulge types are Kormendy & Ho's."},
    }
    published = {"vdb": LIT.VDB16["msigma"], "kh": LIT.KH13["msigma_eq3"], "mm": LIT.MM13["all"]}
    mw, m87 = C.MILKY_WAY, C.M87
    king = res["feedback"]["king_M_at_200"]
    return {
        "samples": samples,
        "published": {k: {"a": round(p["a"], 3), "b": p["b"], "eb": p["eb"], "eps": p["eps"]} for k, p in published.items()},
        "anchors": [
            {"n": "Milky Way (Sgr A*, from S2's orbit, project H1)", "s": float(np.log10(mw["sigma"])),
             "m": float(np.log10(mw["M"]))},
            {"n": "M87* (from its shadow, Event Horizon Telescope)", "s": float(np.log10(m87["sigma"])),
             "m": float(np.log10(m87["M"]))},
        ],
        "king": {"a": float(np.log10(king)), "b": 4.0},
        "headline": {"a": fits["vdb16"].a, "b": fits["vdb16"].b, "eb": fits["vdb16"].eb, "eps": fits["vdb16"].eps},
        "pivot": C.SIGMA_PIVOT,
    }


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>M-σ Explorer</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --wash: #ecebe6; --fit: #52514e;
  --c1: #2a78d6; --c2: #eb6834; --c3: #1baf7a; --ctx: #b9b8b1; --king: #1c5cab;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
    --grid: #2f2f2d; --axis: #4a4a46; --wash: #2a2a28; --fit: #c3c2b7;
    --c1: #3987e5; --c2: #d95926; --c3: #199e70; --ctx: #6e6d68; --king: #86b6ef;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
  --grid: #2f2f2d; --axis: #4a4a46; --wash: #2a2a28; --fit: #c3c2b7;
  --c1: #3987e5; --c2: #d95926; --c3: #199e70; --ctx: #6e6d68; --king: #86b6ef;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1100px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 14px; max-width: 800px; }
.row { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin-bottom: 10px; }
.row .lab { color: var(--muted); font-size: 0.85rem; margin-right: 2px; }
button { font: inherit; font-size: 0.88rem; padding: 5px 11px; border-radius: 8px; cursor: pointer;
         border: 1px solid var(--axis); background: var(--surface); color: var(--ink); }
button[aria-pressed="true"] { background: var(--ink); color: var(--surface); border-color: var(--ink); }
button:focus-visible, input:focus-visible { outline: 2px solid var(--c1); outline-offset: 2px; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 300px; gap: 16px; }
@media (max-width: 860px) { .layout { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border-radius: 12px; padding: 12px; min-width: 0; }
svg#chart { width: 100%; height: auto; display: block; }
fieldset { border: none; margin: 0 0 10px; padding: 0; }
legend { font-size: 0.8rem; color: var(--muted); margin-bottom: 4px; padding: 0; }
.check { font-size: 0.88rem; color: var(--ink2); display: flex; gap: 7px; align-items: center; margin: 2px 0; }
.readout { font-size: 0.88rem; color: var(--ink2); border-top: 1px solid var(--grid); padding-top: 10px; }
.readout .big { font-size: 1.6rem; color: var(--ink); font-weight: 600; line-height: 1.2; }
.readout b { color: var(--ink); }
.note { font-size: 0.85rem; color: var(--ink2); margin: 8px 2px 0; }
.legendrow { display: flex; flex-wrap: wrap; gap: 14px; font-size: 0.8rem; color: var(--ink2); margin: 6px 2px 0; }
.legendrow span { display: inline-flex; align-items: center; gap: 6px; white-space: nowrap; }
#tip { position: fixed; pointer-events: none; background: var(--surface); color: var(--ink);
       border: 1px solid var(--axis); border-radius: 6px; padding: 6px 9px; font-size: 0.8rem; display: none;
       box-shadow: 0 2px 8px rgba(0,0,0,0.15); max-width: 300px; }
#tip .v { font-weight: 600; }
#tip .s { color: var(--ink2); }
details { margin-top: 16px; font-size: 0.85rem; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; width: 100%; }
td, th { padding: 2px 8px; text-align: left; border-bottom: 1px solid var(--grid); white-space: nowrap; }
.tablewrap { max-height: 320px; overflow: auto; }
.axis text { fill: var(--muted); font-size: 11px; }
.axis line { stroke: var(--axis); }
.gridline { stroke: var(--grid); }
.lbl { fill: var(--ink2); font-size: 12px; }
footer { color: var(--muted); font-size: 0.78rem; margin-top: 18px; }
</style>
</head>
<body>
<main>
  <h1>Weigh the M-σ relation yourself</h1>
  <p class="lede">Every point is a galaxy whose central black hole has been weighed. Heavier black holes sit in
  galaxies whose stars move faster (larger velocity dispersion σ), but the <i>slope</i> of that trend depends on
  which galaxies you count. Switch kinds of galaxies in and out and the line is refitted here, in your browser.</p>
  <div class="row" role="group" aria-label="compilation"><span class="lab">Compilation</span><span id="srcs"></span></div>
  <div class="row" role="group" aria-label="colour by"><span class="lab">Colour by</span>
    <button id="cm" aria-pressed="true">how the mass was measured</button>
    <button id="cb" aria-pressed="false">bulge type (Kormendy &amp; Ho)</button>
    <button id="theme" style="margin-left:auto" aria-label="switch light or dark theme">theme</button>
  </div>
  <div class="layout">
    <div class="card">
      <svg id="chart" viewBox="0 0 760 520" role="img" aria-label="black hole mass against velocity dispersion"></svg>
      <div class="legendrow" id="legend"></div>
      <p class="note" id="note"></p>
    </div>
    <div class="card">
      <fieldset><legend>Measurements</legend>
        <label class="check"><input type="checkbox" id="i-stars" checked> stellar dynamics</label>
        <label class="check"><input type="checkbox" id="i-gas" checked> gas dynamics</label>
        <label class="check"><input type="checkbox" id="i-maser" checked> water masers</label>
        <label class="check"><input type="checkbox" id="i-reverb" checked> reverberation (light echoes)</label>
        <label class="check"><input type="checkbox" id="i-upper" checked> upper limits</label>
      </fieldset>
      <fieldset><legend>Galaxies (Kormendy &amp; Ho's classification)</legend>
        <label class="check"><input type="checkbox" id="i-pseudo" checked> pseudobulges</label>
        <label class="check"><input type="checkbox" id="i-aside" checked> ones Kormendy &amp; Ho set aside</label>
      </fieldset>
      <fieldset><legend>Lines</legend>
        <label class="check"><input type="checkbox" id="l-pub" checked> the paper's published fit</label>
        <label class="check"><input type="checkbox" id="l-king"> King's feedback prediction (slope 4)</label>
      </fieldset>
      <div class="row"><button id="kh-like">Kormendy &amp; Ho's choices</button><button id="all">everything</button></div>
      <div class="readout" id="readout" aria-live="polite"></div>
    </div>
  </div>
  <details><summary>Table view of the galaxies in the current fit</summary><div class="tablewrap"><table id="tbl"></table></div></details>
  <footer>Data: van den Bosch 2016 (ApJ 831, 134; author's GitHub table), Kormendy &amp; Ho 2013 (ARA&amp;A 51, 511;
  Tables 2-3 from the arXiv source), McConnell &amp; Ma 2013 (ApJ 764, 184; Table A1). Milky Way: this atlas's
  project H1; M87*: EHT Collaboration 2019. Fit: maximum likelihood with intrinsic scatter in log M, errors in
  both axes, the distribution of true σ modelled as a Gaussian (Kelly 2007), upper limits censored; ± from the
  curvature of the likelihood. Built by AI (Claude) for a hobby atlas; see README.</footer>
</main>
<div id="tip" role="tooltip"></div>
<script>
const D = __DATA__;
const NS = "http://www.w3.org/2000/svg";
const LP = Math.log10(D.pivot);
const state = { src: "vdb", color: "method" };
const INC = ["stars", "gas", "maser", "reverb", "upper", "pseudo", "aside"];
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const SHAPE = { stars: "circle", gas: "square", maser: "diamond", reverb: "triangle" };
const METHOD_NAME = { stars: "stellar dynamics", gas: "gas dynamics", maser: "water maser", reverb: "reverberation" };
const KH_NAME = { elliptical: "elliptical", classical: "classical bulge", pseudo: "pseudobulge", unknown: "not classified" };
function colorOf(r) {
  if (state.color === "method") return { stars: "--c1", gas: "--c2", maser: "--c3", reverb: "--ctx" }[r.t];
  return { elliptical: "--c1", classical: "--c2", pseudo: "--c3", unknown: "--ctx" }[r.kh];
}
function inc(id) { return document.getElementById("i-" + id).checked; }
function included(r) {
  if (!inc(r.t)) return false;
  if (r.u && !inc("upper")) return false;
  if (r.kh === "pseudo" && !inc("pseudo")) return false;
  if (r.x && !inc("aside")) return false;
  return true;
}

// ---------- statistics ----------
function logPhi(z) {          // log of the standard normal CDF, stable in both tails
  const x = Math.abs(z) / Math.SQRT2, t = 1 / (1 + 0.5 * x);
  const lerfc = Math.log(t) - x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))));
  return z < 0 ? Math.log(0.5) + lerfc : Math.log1p(-0.5 * Math.exp(lerfc));
}
function negLogL(p, d) {     // p = [a, b, eps, mu, ln tau]
  const [a, b, eps, mu, lt] = p;
  if (!(eps > 0 && eps < 3 && b > -20 && b < 30 && lt > Math.log(1e-3) && lt < Math.log(5))) return 1e30;
  const t2 = Math.exp(2 * lt);
  let s = 0;
  for (const r of d) {
    const vx = t2 + r.ex * r.ex, cxy = b * t2, dx = r.x - mu;
    const lpx = -0.5 * (dx * dx / vx + Math.log(2 * Math.PI * vx));
    const my = a + b * mu + cxy / vx * dx;
    const sy2 = b * b * t2 + eps * eps + r.ey * r.ey - cxy * cxy / vx;
    const z = (r.y - my) / Math.sqrt(sy2);
    const lpy = r.u ? logPhi(z) : -0.5 * (z * z + Math.log(2 * Math.PI * sy2));
    s += r.w * (lpx + lpy);
  }
  return -s;
}
function nelderMead(f, x0, steps, iters) {
  const n = x0.length;
  let simplex = [x0.slice()];
  for (let i = 0; i < n; i++) { const v = x0.slice(); v[i] += steps[i]; simplex.push(v); }
  let vals = simplex.map(f);
  for (let k = 0; k < iters; k++) {
    const idx = vals.map((v, i) => i).sort((i, j) => vals[i] - vals[j]);
    simplex = idx.map(i => simplex[i]); vals = idx.map(i => vals[i]);
    if (Math.abs(vals[n] - vals[0]) < 1e-10) break;
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += simplex[i][j] / n;
    const pt = (t) => c.map((cj, j) => cj + t * (simplex[n][j] - cj));
    const xr = pt(-1), fr = f(xr);
    if (fr < vals[0]) { const xe = pt(-2), fe = f(xe); if (fe < fr) { simplex[n] = xe; vals[n] = fe; } else { simplex[n] = xr; vals[n] = fr; } }
    else if (fr < vals[n - 1]) { simplex[n] = xr; vals[n] = fr; }
    else {
      const xc = fr < vals[n] ? pt(-0.5) : pt(0.5), fc = f(xc);
      if (fc < Math.min(fr, vals[n])) { simplex[n] = xc; vals[n] = fc; }
      else { for (let i = 1; i <= n; i++) { simplex[i] = simplex[i].map((v, j) => simplex[0][j] + 0.5 * (v - simplex[0][j])); vals[i] = f(simplex[i]); } }
    }
  }
  return simplex[0];
}
function invert(A) {
  const n = A.length, M = A.map((r, i) => r.concat(Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))));
  for (let c = 0; c < n; c++) {
    let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    const pv = M[c][c]; if (Math.abs(pv) < 1e-14) return null;
    for (let j = 0; j < 2 * n; j++) M[c][j] /= pv;
    for (let r = 0; r < n; r++) if (r !== c) { const f = M[r][c]; for (let j = 0; j < 2 * n; j++) M[r][j] -= f * M[c][j]; }
  }
  return M.map(r => r.slice(n));
}
function fit(rows) {
  const d = rows.map(r => ({ x: r.s - LP, ex: r.es || 0, y: r.m, ey: (r.em == null || !isFinite(r.em)) ? 0 : r.em, u: r.u, w: r.w }));
  if (d.filter(r => !r.u).length < 4) return null;
  const xs = d.map(r => r.x), mx = xs.reduce((a, b) => a + b, 0) / xs.length;
  const sx = Math.sqrt(xs.reduce((a, b) => a + (b - mx) ** 2, 0) / xs.length) || 0.1;
  const f = p => negLogL(p, d);
  let p = [8.3, 4.5, 0.4, mx, Math.log(sx)];
  for (let k = 0; k < 4; k++) p = nelderMead(f, p, [0.1, 0.5, 0.05, 0.05, 0.1], 3000);
  const h = [1e-3, 5e-3, 1e-3, 1e-3, 1e-3], n = 5, H = [];
  for (let i = 0; i < n; i++) { H.push([]); for (let j = 0; j < n; j++) {
    const e = (a, b) => { const q = p.slice(); q[i] += a * h[i]; q[j] += b * h[j]; return f(q); };
    H[i].push((e(1, 1) - e(1, -1) - e(-1, 1) + e(-1, -1)) / (4 * h[i] * h[j])); } }
  const cov = invert(H);
  const err = i => (cov && cov[i][i] > 0 ? Math.sqrt(cov[i][i]) : NaN);
  return { a: p[0], b: p[1], eps: p[2], ea: err(0), eb: err(1), eeps: err(2),
           n: d.length, nu: d.filter(r => r.u).length };
}

// ---------- drawing ----------
function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function marker(g, shape, x, y, r, fill, stroke, extra) {
  let e;
  if (shape === "circle") e = el("circle", { cx: x, cy: y, r: r }, g);
  else if (shape === "square") e = el("rect", { x: x - r * 0.9, y: y - r * 0.9, width: 1.8 * r, height: 1.8 * r }, g);
  else if (shape === "diamond") e = el("path", { d: `M${x},${y - 1.25 * r}L${x + 1.25 * r},${y}L${x},${y + 1.25 * r}L${x - 1.25 * r},${y}Z` }, g);
  else if (shape === "triangle") e = el("path", { d: `M${x},${y - 1.2 * r}L${x + 1.1 * r},${y + 0.8 * r}L${x - 1.1 * r},${y + 0.8 * r}Z` }, g);
  else if (shape === "caret") e = el("path", { d: `M${x - r},${y - 0.7 * r}L${x + r},${y - 0.7 * r}L${x},${y + 0.9 * r}Z` }, g);
  else {  // five-pointed star
    let d = "";
    for (let k = 0; k < 10; k++) {
      const ang = -Math.PI / 2 + k * Math.PI / 5, rr = k % 2 ? 0.45 * r : r;
      d += (k ? "L" : "M") + (x + rr * Math.cos(ang)).toFixed(2) + "," + (y + rr * Math.sin(ang)).toFixed(2);
    }
    e = el("path", { d: d + "Z" }, g);
  }
  e.style.fill = fill; e.style.stroke = stroke; e.style.strokeWidth = extra && extra.sw ? extra.sw : 1;
  if (extra && extra.op) e.style.opacity = extra.op;
  return e;
}
function setText(node, parts) {   // parts: [[text, className], ...] -- data goes in as text, never HTML
  node.textContent = "";
  for (const [t, c] of parts) { const s = document.createElement(c === "br" ? "br" : "span"); if (c !== "br") { s.textContent = t; if (c) s.className = c; } node.appendChild(s); }
}
const fmt = (v, n) => (isFinite(v) ? v.toFixed(n) : "n/a");

function draw() {
  const S = D.samples[state.src], rows = S.rows, svg = document.getElementById("chart");
  svg.textContent = "";
  const narrow = svg.clientWidth > 0 && svg.clientWidth < 600;
  const W = narrow ? 520 : 760, H = narrow ? 480 : 520, ml = 62, mr = 14, mt = 14, mb = 46;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const xmin = Math.log10(15), xmax = Math.log10(480), ymin = 2.7, ymax = 10.8;
  const X = s => ml + (s - xmin) / (xmax - xmin) * (W - ml - mr);
  const Y = m => mt + (ymax - m) / (ymax - ymin) * (H - mt - mb);
  const ax = el("g", { class: "axis" }, svg);
  for (let p = 3; p <= 10; p++) {
    el("line", { x1: ml, x2: W - mr, y1: Y(p), y2: Y(p), class: "gridline" }, ax);
    const t = el("text", { x: ml - 8, y: Y(p) + 4, "text-anchor": "end" }, ax); t.textContent = "1e" + p;
  }
  for (const s of [20, 30, 50, 70, 100, 150, 200, 300, 400]) {
    el("line", { x1: X(Math.log10(s)), x2: X(Math.log10(s)), y1: mt, y2: H - mb, class: "gridline" }, ax);
    const t = el("text", { x: X(Math.log10(s)), y: H - mb + 16, "text-anchor": "middle" }, ax); t.textContent = s;
  }
  el("line", { x1: ml, x2: W - mr, y1: H - mb, y2: H - mb }, ax);
  el("line", { x1: ml, x2: ml, y1: mt, y2: H - mb }, ax);
  const xl = el("text", { x: (ml + W - mr) / 2, y: H - 8, "text-anchor": "middle", class: "lbl" }, svg);
  xl.textContent = "stellar velocity dispersion σ (km/s)";
  const yl = el("text", { x: 14, y: (mt + H - mb) / 2, "text-anchor": "middle", class: "lbl",
                          transform: `rotate(-90 14 ${(mt + H - mb) / 2})` }, svg);
  yl.textContent = "black hole mass (solar masses)";
  const clip = el("clipPath", { id: "pc" }, el("defs", {}, svg));
  el("rect", { x: ml, y: mt, width: W - ml - mr, height: H - mt - mb }, clip);
  const g = el("g", { "clip-path": "url(#pc)" }, svg);

  const sel = rows.filter(included), F = fit(sel);
  const line = (a, b, attrs) => {
    const x1 = xmin, x2 = xmax, y1 = a + b * (x1 - LP), y2 = a + b * (x2 - LP);
    return el("line", Object.assign({ x1: X(x1), x2: X(x2), y1: Y(y1), y2: Y(y2) }, attrs), g);
  };
  if (F) {
    const band = el("path", { d: `M${X(xmin)},${Y(F.a + F.b * (xmin - LP) + F.eps)}L${X(xmax)},${Y(F.a + F.b * (xmax - LP) + F.eps)}` +
      `L${X(xmax)},${Y(F.a + F.b * (xmax - LP) - F.eps)}L${X(xmin)},${Y(F.a + F.b * (xmin - LP) - F.eps)}Z` }, g);
    band.style.fill = css("--wash");
  }
  if (document.getElementById("l-king").checked) {
    const k = line(D.king.a, D.king.b, {}); k.style.stroke = css("--king"); k.style.strokeWidth = 2;
  }
  if (document.getElementById("l-pub").checked) {
    const P = D.published[state.src];
    const k = line(P.a, P.b, {}); k.style.stroke = css("--ink2"); k.style.strokeWidth = 1.2; k.style.strokeDasharray = "6 4";
  }
  if (F) { const k = line(F.a, F.b, {}); k.style.stroke = css("--fit"); k.style.strokeWidth = 2.4; }

  const pts = [];
  const dots = el("g", {}, g);
  for (const r of rows) {
    const on = included(r), y = r.u ? (r.l != null ? r.l : r.m) : r.m;
    const col = on ? css(colorOf(r)) : "none", stroke = on ? css("--surface") : css("--ctx");
    const shape = r.u ? "caret" : SHAPE[r.t];
    marker(dots, shape, X(r.s), Y(y), on ? 4.4 : 3.6, r.u && on ? css("--muted") : col, stroke, { sw: on ? 0.9 : 1, op: on ? 1 : 0.7 });
    pts.push({ r, x: X(r.s), y: Y(y), on });
  }
  for (const a of D.anchors) {
    marker(dots, "star", X(a.s), Y(a.m), 8, css("--ink"), css("--surface"), { sw: 1 });
    pts.push({ r: { n: a.n, s: a.s, m: a.m, t: null, anchor: true }, x: X(a.s), y: Y(a.m), on: true });
  }
  // Nearest-point hover: the pointer only has to be closest, not on the dot.
  const hit = el("rect", { x: ml, y: mt, width: W - ml - mr, height: H - mt - mb, fill: "transparent" }, svg);
  const tip = document.getElementById("tip");
  hit.addEventListener("mousemove", ev => {
    const b = svg.getBoundingClientRect(), sx = (ev.clientX - b.left) * W / b.width, sy = (ev.clientY - b.top) * H / b.height;
    let best = null, bd = 1e9;
    for (const p of pts) { const d = (p.x - sx) ** 2 + (p.y - sy) ** 2; if (d < bd) { bd = d; best = p; } }
    if (!best || bd > 500) { tip.style.display = "none"; return; }
    const r = best.r, parts = [];
    if (r.anchor) {
      parts.push([`${(10 ** r.m / 1e6).toFixed(1)} million solar masses`, "v"], ["", "br"], [r.n, "s"], ["", "br"],
                 [`σ = ${(10 ** r.s).toFixed(0)} km/s; weighed without any galaxy model`, "s"]);
    } else {
      const m = r.u ? `below ${fmtM(r.l != null ? r.l : r.m)}` : fmtM(r.m);
      parts.push([`${m} solar masses`, "v"], ["", "br"], [r.n, ""], ["", "br"],
                 [`σ = ${(10 ** r.s).toFixed(0)} km/s · ${METHOD_NAME[r.t]} · ${KH_NAME[r.kh]}`, "s"]);
      if (r.d != null) parts.push(["", "br"], [`distance ${r.d} Mpc`, "s"]);
      if (r.r) parts.push(["", "br"], [`source: ${r.r}`, "s"]);
      if (!best.on) parts.push(["", "br"], ["(switched off: not in the fit)", "s"]);
    }
    setText(tip, parts);
    tip.style.display = "block"; tip.style.left = Math.min(ev.clientX + 14, window.innerWidth - 310) + "px"; tip.style.top = (ev.clientY + 12) + "px";
  });
  hit.addEventListener("mouseleave", () => { tip.style.display = "none"; });

  legend(); readout(F, sel); table(sel);
  document.getElementById("note").textContent = S.note;
}
function fmtM(lm) {
  const v = 10 ** lm;
  if (v >= 1e9) return (v / 1e9).toPrecision(2) + " billion";
  if (v >= 1e6) return (v / 1e6).toPrecision(2) + " million";
  return Math.round(v).toLocaleString("en-US");
}
function legend() {
  const L = document.getElementById("legend"); L.textContent = "";
  const items = state.color === "method"
    ? [["stars", "--c1", "circle", "stellar dynamics"], ["gas", "--c2", "square", "gas dynamics"], ["maser", "--c3", "diamond", "masers"], ["reverb", "--ctx", "triangle", "reverberation"]]
    : [["e", "--c1", "circle", "elliptical"], ["c", "--c2", "circle", "classical bulge"], ["p", "--c3", "circle", "pseudobulge"], ["u", "--ctx", "circle", "not classified"]];
  const add = (draw, text) => { const s = document.createElement("span"); const sv = document.createElementNS(NS, "svg");
    sv.setAttribute("width", "26"); sv.setAttribute("height", "12"); draw(sv); s.appendChild(sv); s.appendChild(document.createTextNode(text)); L.appendChild(s); };
  for (const [, c, sh, t] of items) add(sv => marker(sv, state.color === "method" ? sh : "circle", 13, 6, 4.4, css(c), css("--surface")), t);
  add(sv => marker(sv, "caret", 13, 6, 4.4, css("--muted"), css("--surface")), "upper limit");
  add(sv => marker(sv, "circle", 13, 6, 3.6, "none", css("--ctx")), "switched off");
  add(sv => { const l = el("line", { x1: 0, x2: 26, y1: 6, y2: 6 }, sv); l.style.stroke = css("--fit"); l.style.strokeWidth = 2.4; }, "fit to the selection");
  if (document.getElementById("l-pub").checked) add(sv => { const l = el("line", { x1: 0, x2: 26, y1: 6, y2: 6 }, sv); l.style.stroke = css("--ink2"); l.style.strokeWidth = 1.2; l.style.strokeDasharray = "6 4"; }, "published fit");
  if (document.getElementById("l-king").checked) add(sv => { const l = el("line", { x1: 0, x2: 26, y1: 6, y2: 6 }, sv); l.style.stroke = css("--king"); l.style.strokeWidth = 2; }, "King's prediction");
  add(sv => marker(sv, "star", 13, 6, 6, css("--ink"), css("--surface")), "Milky Way, M87*");
}
function readout(F, sel) {
  const R = document.getElementById("readout"); R.textContent = "";
  const P = D.published[state.src];
  const add = (t, c) => { const d = document.createElement("div"); d.textContent = t; if (c) d.className = c; R.appendChild(d); return d; };
  if (!F) { add("Too few masses to fit: switch something back on."); return; }
  add("slope of the fitted line");
  add(`${fmt(F.b, 2)} ± ${fmt(F.eb, 2)}`, "big");
  add(`log M = ${fmt(F.a, 2)} + ${fmt(F.b, 2)} log(σ / 200 km/s)`);
  add(`intrinsic scatter ${fmt(F.eps, 2)} dex (a factor of ${fmt(10 ** F.eps, 1)})`);
  add(`${F.n} galaxies, ${F.nu} of them upper limits`);
  const p = add(""); p.style.marginTop = "8px";
  p.textContent = `Published for this compilation: slope ${P.b} ± ${P.eb}, scatter ${P.eps} dex. ` +
    `Momentum-driven feedback (King 2003) predicts slope 4 and ${(10 ** D.king.a / 1e8).toFixed(1)} × 10⁸ solar masses at 200 km/s.`;
}
function table(sel) {
  const T = document.getElementById("tbl"); T.textContent = "";
  const hr = T.insertRow();
  for (const h of ["galaxy", "σ (km/s)", "log M", "method", "bulge type", "distance (Mpc)", "source"]) { const th = document.createElement("th"); th.textContent = h; hr.appendChild(th); }
  for (const r of sel) {
    const tr = T.insertRow();
    const vals = [r.n, (10 ** r.s).toFixed(0), (r.u ? "< " + (r.l != null ? r.l : r.m).toFixed(2) : r.m.toFixed(2)), METHOD_NAME[r.t], KH_NAME[r.kh], r.d == null ? "" : r.d, r.r];
    for (const v of vals) { const td = tr.insertCell(); td.textContent = v; }
  }
}
function setSrc(src) {
  state.src = src;
  document.querySelectorAll("#srcs button").forEach(b => b.setAttribute("aria-pressed", b.dataset.src === src));
  draw();
}
for (const k of ["vdb", "kh", "mm"]) {
  const b = document.createElement("button"); b.dataset.src = k; b.textContent = D.samples[k].label;
  b.style.marginRight = "6px"; b.addEventListener("click", () => setSrc(k)); document.getElementById("srcs").appendChild(b);
}
document.getElementById("cm").addEventListener("click", () => { state.color = "method"; document.getElementById("cm").setAttribute("aria-pressed", "true"); document.getElementById("cb").setAttribute("aria-pressed", "false"); draw(); });
document.getElementById("cb").addEventListener("click", () => { state.color = "bulge"; document.getElementById("cb").setAttribute("aria-pressed", "true"); document.getElementById("cm").setAttribute("aria-pressed", "false"); draw(); });
for (const id of INC) document.getElementById("i-" + id).addEventListener("change", draw);
for (const id of ["l-pub", "l-king"]) document.getElementById(id).addEventListener("change", draw);
document.getElementById("kh-like").addEventListener("click", () => {
  for (const [id, v] of [["stars", 1], ["gas", 1], ["maser", 1], ["reverb", 0], ["upper", 0], ["pseudo", 0], ["aside", 0]]) document.getElementById("i-" + id).checked = !!v;
  draw();
});
document.getElementById("all").addEventListener("click", () => { for (const id of INC) document.getElementById("i-" + id).checked = true; draw(); });
document.getElementById("theme").addEventListener("click", () => {
  const root = document.documentElement, dark = root.dataset.theme ? root.dataset.theme === "dark" : window.matchMedia("(prefers-color-scheme: dark)").matches;
  root.dataset.theme = dark ? "light" : "dark"; draw();
});
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);
let rt; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(draw, 120); });
setSrc("vdb");
</script>
</body>
</html>
"""


def write(res: dict):
    path = C.FIGURE_DIR / "m_sigma_explorer.html"
    payload = json.dumps(_payload(res), separators=(",", ":"), ensure_ascii=False)
    path.write_text(PAGE.replace("__DATA__", payload.replace("</", "<\\/")))
    return path
