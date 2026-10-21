"""Write figures/scorecard_explorer.html: the scorecard as an offline, interactive page.

One self-contained file (no external scripts or fonts). It shows the 46 tests as
precision vs strength of gravity (the same broken axis as figure 01) or as a
report card of pulls (figure 02). Hovering or tapping a test shows its reference,
quantity, value and source. There are family filters, a sortable table, and light
and dark themes.

Colours: blue, orange and aqua are validated categorical slots in both themes
(dataviz all-pairs check). The fourth family, galaxies and cosmology, is violet in
light mode, as in the static figures, and neutral grey in dark mode, where no
fourth hue clears the colour-vision checks. Every family also has its own marker
shape. The atlas's own measurements are stars in their family's colour with an
ink outline.
"""

from __future__ import annotations

import json
import math

import pandas as pd

import config as C

FAMILIES = {
    "equivalence": {"label": "equivalence principle", "shape": "circle"},
    "weak": {"label": "weak field (solar system)", "shape": "square"},
    "strong": {"label": "strong field and dynamics", "shape": "diamond"},
    "cosmology": {"label": "galaxies and cosmology", "shape": "triangle"},
}

#: label priority after the atlas's own points; the page places as many as fit without overlap
PRIORITY = ["MICROSCOPE", "LLR Earth–Moon", "Eöt-Wash", "atoms", "triple pulsar", "Cassini", "Mercury perihelion",
            "VLBI deflection", "Hulse–Taylor Ṗb", "double pulsar Ṗb", "growth index", "GP-B geodetic",
            "GP-B frame dragging", "Galileo clocks", "GP-A", "Pound–Snider", "S2 precession", "M87* shadow",
            "Sgr A* shadow", "GW250114 ringdown", "GW250114 inspiral", "antimatter", "spin precession", "LARES 2",
            "LLR geodetic", "Sun redshift", "S2 redshift", "J0348+0432", "double pulsar s", "IMR consistency",
            "ringdown damping", "E_G", "Σ0 (DESI)", "μ0 (DESI)", "galaxy γ"]

TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Gravity Test Scorecard</title>
<style>
:root {
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #6f6e69;
  --grid: #e8e7e1; --axis: #b9b8ad; --ring: rgba(11,11,11,0.10);
  --equivalence: #2a78d6; --weak: #eb6834; --strong: #1baf7a; --cosmology: #4a3aa7;
  --band1: #e4e3dc; --band2: #efeee9; --focus: #2a78d6;
  color-scheme: light;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #9a9890;
    --grid: #2a2a28; --axis: #4a4945; --ring: rgba(255,255,255,0.10);
    --equivalence: #3987e5; --weak: #d95926; --strong: #199e70; --cosmology: #c3c2b7;
    --band1: #33332f; --band2: #252523; --focus: #3987e5;
    color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #9a9890;
  --grid: #2a2a28; --axis: #4a4945; --ring: rgba(255,255,255,0.10);
  --equivalence: #3987e5; --weak: #d95926; --strong: #199e70; --cosmology: #c3c2b7;
  --band1: #33332f; --band2: #252523; --focus: #3987e5;
  color-scheme: dark;
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body { margin: 0; background: var(--page); color: var(--ink);
       font: 15px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width: 1120px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 12px; max-width: 72ch; }
.controls { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin: 10px 0 8px; }
.controls .sep { flex: 1 1 auto; }
.chip, .seg { font: inherit; font-size: 0.85rem; border: 1px solid var(--axis); background: var(--surface);
  color: var(--ink2); border-radius: 999px; padding: 4px 11px 4px 8px; display: inline-flex; gap: 6px; align-items: center; }
button.chip, button.seg { cursor: pointer; }
.seg { padding: 4px 11px; }
button[aria-pressed="true"] { color: var(--ink); border-color: var(--ink2); }
button.seg[aria-pressed="true"] { background: var(--band1); }
button.chip[aria-pressed="false"] { opacity: 0.5; text-decoration: line-through; }
.legend-star { border-style: dashed; }
button:focus-visible, .pt:focus-visible, th:focus-visible { outline: 2px solid var(--focus); outline-offset: 2px; }
.chip svg, .legend-star svg, td svg { width: 14px; height: 14px; flex: none; overflow: visible; }
#summary { color: var(--ink2); font-size: 0.88rem; margin: 2px 0 8px; min-height: 1.3em; }
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: 10px; padding: 8px 6px; }
#chartbox { width: 100%; }
#chart { display: block; width: 100%; height: auto; touch-action: manipulation; }
#chart .grid { stroke: var(--grid); stroke-width: 1; }
#chart .axis { stroke: var(--axis); stroke-width: 1; fill: none; }
#chart .zero { stroke: var(--ink); stroke-width: 1.5; }
#chart text { font-family: inherit; }
#chart .tick { fill: var(--muted); font-size: 11px; }
#chart .axislabel { fill: var(--ink2); font-size: 12px; }
#chart .lbl { fill: var(--ink2); font-size: 11px; }
#chart .lbl.atlas { fill: var(--ink); font-weight: 650; }
#chart .famhead { fill: var(--ink); font-size: 12px; font-weight: 650; }
#chart .ring { fill: none; stroke: var(--ink); stroke-width: 2; pointer-events: none; }
#chart .leader { stroke: var(--muted); stroke-width: 1; }
#chart .pt { outline: none; }
#tip { position: fixed; z-index: 10; width: max-content; max-width: min(340px, calc(100vw - 16px)); background: var(--surface);
  color: var(--ink); border: 1px solid var(--axis); border-radius: 8px; padding: 9px 11px; font-size: 0.82rem; line-height: 1.38;
  box-shadow: 0 4px 18px rgba(0,0,0,0.20); pointer-events: none; display: none; }
#tip .t { font-weight: 650; display: block; margin-bottom: 2px; }
#tip .k { color: var(--muted); display: block; margin-bottom: 4px; }
#tip .v { font-weight: 650; font-variant-numeric: tabular-nums; }
.note { color: var(--muted); font-size: 0.8rem; margin: 8px 0 0; max-width: 110ch; }
.tablewrap { overflow-x: auto; margin-top: 18px; }
table { width: 100%; border-collapse: collapse; font-size: 0.84rem; }
th, td { text-align: left; padding: 5px 6px; border-bottom: 1px solid var(--grid); vertical-align: top; }
th { color: var(--ink2); font-weight: 600; cursor: pointer; white-space: nowrap; user-select: none; }
td.num { font-variant-numeric: tabular-nums; white-space: nowrap; }
td.name { min-width: 15em; }
td.name svg { vertical-align: -2px; margin-right: 7px; }
td.ref { color: var(--ink2); min-width: 10em; }
</style>
</head>
<body>
<main>
<h1>Is Einstein still right?</h1>
<p class="lede">__N__ tests of general relativity, from laboratory balances to black holes. Each point is one test,
placed by the 1σ uncertainty of its measured deviation from GR (or of an effect GR forbids) against the strength of
the gravity it probes, GM/(rc²). Stars are measurements made in this atlas. Hover over or tap a point for its source.</p>
<div class="controls" id="filters"></div>
<p id="summary" aria-live="polite"></p>
<div class="card"><div id="chartbox"><svg id="chart" role="img" aria-label="Scorecard chart; the table below lists every test"></svg></div></div>
<p class="note">Values were checked against the primary papers (literature/checks/*.csv in project N7 of the atlas) and
converted to 1σ: 90% limits ÷ 1.645, 95% ÷ 1.96, probable errors × 1.4826, statistical and systematic errors added in
quadrature. GM/(rc²) is a characteristic value for each test, not a unique number. Overlapping points are nudged sideways;
the tooltips and the table give the exact values. In the report card, if GR is right and the errors are honest, about
68% of tests should land inside ±1σ and 95% inside ±2σ.</p>
<div class="tablewrap"><table id="table"></table></div>
<div id="tip" role="tooltip"></div>
</main>
<script>
const DATA = __DATA__;
const FAM = __FAMILIES__;
const PRIORITY = __PRIORITY__;
const ORDER = Object.keys(FAM);
const NS = "http://www.w3.org/2000/svg";
const state = { view: "precision", fams: new Set(ORDER), sort: { key: "x", dir: 1 } };
let live = { pts: [], W: 0, H: 0, mode: "points", ring: null };
let hovered = null;

function el(tag, attrs, parent, text) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (text !== undefined) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}
function esc(s) { return String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]); }
function sup(n) { const s = "⁰¹²³⁴⁵⁶⁷⁸⁹"; return (n < 0 ? "⁻" : "") + String(Math.abs(n)).split("").map(c => s[+c]).join(""); }
function minus(s) { return s.replace(/^-/, "−"); }
function fmt(v) {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 0.01 && a < 1000) return minus(String(+v.toPrecision(3)));
  let e = Math.floor(Math.log10(a)), mant = +(v / Math.pow(10, e)).toPrecision(3);
  if (Math.abs(mant) >= 10) { mant = +(mant / 10).toPrecision(3); e += 1; }
  return minus(String(mant)) + "×10" + sup(e);
}
function pullText(p, digits) { return (p >= 0 ? "+" : "−") + Math.abs(p).toFixed(digits) + "σ"; }
function shown() { return DATA.filter(d => state.fams.has(d.family)); }

// --- markers: one shape per family; the atlas's own points are stars in their family's colour ------------
function marker(parent, d, x, y) {
  const star = !!d.atlas;
  const a = { fill: `var(--${d.family})`, class: "mk", stroke: star ? "var(--ink)" : "var(--surface)", "stroke-width": star ? 1.2 : 1.5 };
  if (star) {
    const pts = [];
    for (let i = 0; i < 10; i++) {
      const ang = -Math.PI / 2 + i * Math.PI / 5, rr = i % 2 ? 4.2 : 9.2;
      pts.push((x + rr * Math.cos(ang)).toFixed(2) + "," + (y + rr * Math.sin(ang)).toFixed(2));
    }
    return el("polygon", { ...a, points: pts.join(" ") }, parent);
  }
  const s = FAM[d.family].shape;
  if (s === "circle") return el("circle", { ...a, cx: x, cy: y, r: 5.5 }, parent);
  if (s === "square") return el("rect", { ...a, x: x - 4.9, y: y - 4.9, width: 9.8, height: 9.8, rx: 1.5 }, parent);
  if (s === "diamond") return el("polygon", { ...a, points: `${x},${y - 6.8} ${x + 6.8},${y} ${x},${y + 6.8} ${x - 6.8},${y}` }, parent);
  return el("polygon", { ...a, points: `${x},${y - 6.6} ${x + 6.2},${y + 4.4} ${x - 6.2},${y + 4.4}` }, parent);
}
function glyph(d) {
  const s = document.createElementNS(NS, "svg");
  s.setAttribute("viewBox", "-10 -10 20 20"); s.setAttribute("aria-hidden", "true");
  marker(s, d, 0, 0);
  return s;
}

// --- tooltip -----------------------------------------------------------------------------------------------
const tip = document.getElementById("tip");
function tipHTML(d) {
  return `<span class="t">${esc(d.label)}</span>` +
    `<span class="k">${esc(FAM[d.family].label)} · ${d.year} · ${esc(d.ref)}</span>` +
    `${esc(d.quantity)}: <span class="v">${fmt(d.deviation)} ± ${fmt(d.sigma)}</span> (1σ)<br>` +
    `distance from GR: <span class="v">${pullText(d.pull, 2)}</span><br>` +
    `strength of gravity GM/(rc²) ≈ ${fmt(d.x)}` +
    (d.note ? `<span class="k" style="margin:4px 0 0">${esc(d.note)}</span>` : "") +
    (d.source ? `<span class="k" style="margin:4px 0 0">source: ${esc(d.source)}</span>` : "");
}
function placeTip(cx, cy) {
  tip.style.display = "block";
  const w = tip.offsetWidth, h = tip.offsetHeight, vw = document.documentElement.clientWidth, vh = window.innerHeight;
  let x = cx + 16, y = cy - h / 2;
  if (x + w > vw - 8) x = cx - 16 - w;
  if (x < 8) {                      // no room beside the point (phones): go below it, or above
    x = Math.max(8, Math.min(vw - w - 8, cx - w / 2));
    y = cy + 18 + h < vh - 8 ? cy + 18 : cy - 18 - h;
  }
  tip.style.left = x + "px";
  tip.style.top = Math.max(8, Math.min(vh - h - 8, y)) + "px";
}
function setHover(p) {
  hovered = p || null;
  if (live.ring) live.ring.setAttribute("visibility", p ? "visible" : "hidden");
  if (!p) { tip.style.display = "none"; return; }
  if (live.ring) { live.ring.setAttribute("cx", p.px); live.ring.setAttribute("cy", p.py); live.ring.setAttribute("r", p.rc + 4); }
  tip.innerHTML = tipHTML(p.d);
  const svg = document.getElementById("chart"), r = svg.getBoundingClientRect(), k = r.width / live.W;
  placeTip(r.left + p.px * k, r.top + p.py * k);
}
function nearest(ev) {
  const svg = document.getElementById("chart"), r = svg.getBoundingClientRect(), k = live.W / r.width;
  const mx = (ev.clientX - r.left) * k, my = (ev.clientY - r.top) * k;
  if (live.mode === "rows") return live.pts.find(p => Math.abs(p.py - my) <= live.rowH / 2 && mx > 0) || null;
  let best = null, bd = 22 * 22;
  for (const p of live.pts) { const d2 = (p.px - mx) ** 2 + (p.py - my) ** 2; if (d2 < bd) { bd = d2; best = p; } }
  return best;
}
function focusable(g, p) {
  g.setAttribute("tabindex", "0"); g.setAttribute("class", "pt");
  g.setAttribute("aria-label", `${p.d.label}: ${fmt(p.d.deviation)} ± ${fmt(p.d.sigma)}, ${pullText(p.d.pull, 1)} from GR`);
  g.addEventListener("focus", () => setHover(p));
  g.addEventListener("blur", () => setHover(null));
}

// --- chart 1: precision vs strength of gravity (broken y axis at 1e-3, as in figure 01) ----------------------
function overlaps(a, b) { return a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1; }
function drawPrecision(svg, W) {
  const narrow = W < 640;
  const m = { l: narrow ? 46 : 66, r: narrow ? 6 : 12, t: 10, b: 44 };
  const H = narrow ? 600 : 700, gap = 18;
  const ph = H - m.t - m.b - gap, topH = Math.round(ph * 0.55);
  const panels = [
    { y0: m.t, h: topH, lo: -3, hi: 0.4, ticks: [0, -1, -2, -3] },
    { y0: m.t + topH + gap, h: ph - topH, lo: -15.3, hi: -3, ticks: [-4, -6, -8, -10, -12, -14] },
  ];
  const xl = [-12.3, 0.7], pw = W - m.l - m.r;
  const X = lx => m.l + (lx - xl[0]) / (xl[1] - xl[0]) * pw;
  const panelOf = ly => (ly >= -3 ? panels[0] : panels[1]);
  const Y = ly => { const p = panelOf(ly); return p.y0 + (p.hi - ly) / (p.hi - p.lo) * p.h; };
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`); svg.setAttribute("width", W); svg.setAttribute("height", H);

  const axes = el("g", {}, svg);
  for (let e = -12; e <= 0; e += narrow ? 4 : 2) {
    for (const p of panels) el("line", { x1: X(e), x2: X(e), y1: p.y0, y2: p.y0 + p.h, class: "grid" }, axes);
    el("text", { x: X(e), y: H - m.b + 16, "text-anchor": "middle", class: "tick" }, axes, "10" + sup(e));
  }
  for (const p of panels) {
    for (const e of p.ticks) {
      el("line", { x1: m.l, x2: m.l + pw, y1: Y(e), y2: Y(e), class: "grid" }, axes);
      el("text", { x: m.l - 7, y: Y(e) + 4, "text-anchor": "end", class: "tick" }, axes, "10" + sup(e));
    }
    el("line", { x1: m.l, x2: m.l, y1: p.y0, y2: p.y0 + p.h, class: "axis" }, axes);
  }
  el("line", { x1: m.l, x2: m.l + pw, y1: H - m.b, y2: H - m.b, class: "axis" }, axes);
  for (const y of [panels[0].y0 + panels[0].h, panels[1].y0])           // break marks on the y axis
    el("line", { x1: m.l - 5, x2: m.l + 5, y1: y + 3, y2: y - 3, class: "axis", "stroke-width": 1.5 }, axes);
  el("text", { x: m.l + pw - 4, y: panels[1].y0 - 5, "text-anchor": "end", class: "tick" }, axes,
     narrow ? "scale changes at 10⁻³" : "the vertical scale changes at 10⁻³");
  el("text", { x: m.l + pw / 2, y: H - 8, "text-anchor": "middle", class: "axislabel" }, axes,
     "strength of gravity probed, GM/(rc²)" + (narrow ? "" : "   (weak → strong)"));
  const yc = m.t + ph / 2 + gap / 2;
  el("text", { x: 13, y: yc, "text-anchor": "middle", class: "axislabel", transform: `rotate(-90 13 ${yc})` }, axes,
     narrow ? "1σ uncertainty" : "1σ uncertainty of the deviation (smaller = more precise)");

  // points; nudge overlapping markers apart sideways (deterministic; the exact values are in the tooltip)
  const pts = shown().map(d => {
    const ly = Math.log10(d.sigma);
    return { d, px: X(Math.log10(d.x)), py: Y(ly), rc: d.atlas ? 8.5 : 6.3, panel: panelOf(ly) };
  });
  for (let it = 0; it < 120; it++) {
    let moved = false;
    for (let a = 0; a < pts.length; a++) for (let b = a + 1; b < pts.length; b++) {
      const p = pts[a], q = pts[b], need = p.rc + q.rc + 1.5, dy = q.py - p.py;
      if (Math.abs(dy) >= need) continue;
      const want = Math.sqrt(need * need - dy * dy), dx = q.px - p.px;
      if (Math.abs(dx) >= want - 1e-6) continue;
      const s = dx >= 0 ? 1 : -1, push = (want - Math.abs(dx)) / 2 + 0.01;
      p.px -= s * push; q.px += s * push; moved = true;
    }
    if (!moved) break;
  }
  const marks = el("g", {}, svg);
  for (const p of pts.filter(p => !p.d.atlas).concat(pts.filter(p => p.d.atlas))) {
    const g = el("g", {}, marks);
    marker(g, p.d, p.px, p.py);
    focusable(g, p);
  }

  // direct labels: the atlas's points first, then by PRIORITY, wherever they fit without touching anything
  const layer = el("g", {}, svg);
  const boxes = pts.map(p => ({ x0: p.px - p.rc, x1: p.px + p.rc, y0: p.py - p.rc, y1: p.py + p.rc, p }));
  const placed = [];
  const byShort = new Map(pts.map(p => [p.d.short, p]));
  const order = pts.filter(p => p.d.atlas).sort((a, b) => a.d.sigma - b.d.sigma).map(p => p.d.short)
    .concat(PRIORITY, pts.map(p => p.d.short));
  const done = new Set();
  const TOP = 10.5, BOT = 3;                         // ink extent of 11-px text above and below its baseline
  const fits = (b, p) => !(b.x0 < m.l + 3 || b.x1 > W - 2 || b.y0 < p.panel.y0 + 1 || b.y1 > p.panel.y0 + p.panel.h - 1 ||
                           boxes.some(o => o.p !== p && overlaps(b, o)) || placed.some(o => overlaps(b, o)));
  const clearPath = (p, x1, y1, x2, y2) => {         // a leader line must not cross other markers or labels
    for (let f = 0; f <= 1; f += 0.05) {
      const q = { x0: x1 + f * (x2 - x1) - 0.5, y0: y1 + f * (y2 - y1) - 0.5 };
      q.x1 = q.x0 + 1; q.y1 = q.y0 + 1;
      if (boxes.some(o => o.p !== p && overlaps(q, o)) || placed.some(o => overlaps(q, o))) return false;
    }
    return true;
  };
  for (const s of order) {
    const p = byShort.get(s);
    if (!p || done.has(s)) continue;
    done.add(s);
    const t = el("text", { class: "lbl" + (p.d.atlas ? " atlas" : ""), x: 0, y: 0 }, layer, s);
    const w = t.getComputedTextLength(), g = p.rc + 3;
    const R = p.px + g, L = p.px - g - w, Cx = p.px - w / 2, dR = p.px + g * 0.6, dL = p.px - g * 0.6 - w;
    const up = p.py - p.rc - 4, down = p.py + p.rc + TOP + 1.5;
    const cands = [[R, p.py + 4], [L, p.py + 4], [Cx, up], [Cx, down], [R, p.py + 13], [R, p.py - 5],
                   [L, p.py + 13], [L, p.py - 5], [dR, up + 2], [dR, down - 2], [dL, up + 2], [dL, down - 2]];
    let ok = false;
    for (const [x0, bl] of cands) {
      const b = { x0: x0 - 1, x1: x0 + w + 1, y0: bl - TOP, y1: bl + BOT };
      if (!fits(b, p)) continue;
      t.setAttribute("x", x0.toFixed(1)); t.setAttribute("y", bl.toFixed(1));
      placed.push(b); ok = true; break;
    }
    // the atlas's own points may sit a little further out, with a leader line back to the star
    for (const rad of p.d.atlas && !ok ? [24, 34, 46] : []) {
      for (const deg of [-35, 35, -145, 145, -65, 65, -115, 115, -90, 90]) {
        const a = deg * Math.PI / 180, c = Math.cos(a), sn = Math.sin(a);
        const ax = p.px + rad * c, ay = p.py + rad * sn, vertical = Math.abs(c) < 0.2;
        const x0 = vertical ? ax - w / 2 : c > 0 ? ax + 2 : ax - w - 2;
        const bl = vertical ? (sn < 0 ? ay - BOT - 1 : ay + TOP + 1) : ay + 4;
        const b = { x0: x0 - 1, x1: x0 + w + 1, y0: bl - TOP, y1: bl + BOT };
        const sx = p.px + (p.rc + 1) * c, sy = p.py + (p.rc + 1) * sn;
        if (!fits(b, p) || !clearPath(p, sx, sy, ax, ay)) continue;
        el("line", { x1: sx.toFixed(1), y1: sy.toFixed(1), x2: ax.toFixed(1), y2: ay.toFixed(1), class: "leader" }, layer);
        t.setAttribute("x", x0.toFixed(1)); t.setAttribute("y", bl.toFixed(1));
        placed.push(b); ok = true; break;
      }
      if (ok) break;
    }
    if (!ok) t.remove();
  }
  live = { pts, W, H, mode: "points", ring: el("circle", { class: "ring", r: 10, visibility: "hidden" }, svg) };
}

// --- chart 2: the report card (pulls) --------------------------------------------------------------------
function drawPulls(svg, W) {
  const narrow = W < 640;
  const rows = [];
  for (const f of ORDER) {
    if (!state.fams.has(f)) continue;
    const ds = DATA.filter(d => d.family === f).sort((a, b) => a.x - b.x);
    if (ds.length) { rows.push({ head: f }); ds.forEach(d => rows.push({ d })); }
  }
  const rowH = 21, headH = 26, labW = narrow ? Math.min(150, Math.round(W * 0.42)) : 230;
  const m = { l: labW + 14, r: narrow ? 46 : 64, t: 26, b: 30 };
  let y = m.t;
  for (const r of rows) { r.y = y; y += r.head ? headH : rowH; }
  const H = y + m.b, bottom = y;
  const X = v => m.l + (Math.max(-4.6, Math.min(4.6, v)) + 4.6) / 9.2 * (W - m.l - m.r);
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`); svg.setAttribute("width", W); svg.setAttribute("height", H);
  const g0 = el("g", {}, svg);
  el("rect", { x: X(-2), y: m.t, width: X(2) - X(-2), height: bottom - m.t, fill: "var(--band2)" }, g0);
  el("rect", { x: X(-1), y: m.t, width: X(1) - X(-1), height: bottom - m.t, fill: "var(--band1)" }, g0);
  for (let v = -4; v <= 4; v += narrow ? 2 : 1) {
    if (v !== 0) el("line", { x1: X(v), x2: X(v), y1: m.t, y2: bottom, class: "grid" }, g0);
    for (const ty of [m.t - 8, bottom + 16])
      el("text", { x: X(v), y: ty, "text-anchor": "middle", class: "tick" }, g0, v === 0 ? "GR" : (v > 0 ? "+" : "−") + Math.abs(v) + "σ");
  }
  el("line", { x1: X(0), x2: X(0), y1: m.t, y2: bottom, class: "zero" }, g0);
  const pts = [];
  for (const r of rows) {
    if (r.head) {
      el("text", { x: 4, y: r.y + headH - 8, class: "famhead" }, g0, FAM[r.head].label);
      continue;
    }
    const cy = r.y + rowH / 2, d = r.d;
    const t = el("text", { x: m.l - 10, y: cy + 4, "text-anchor": "end", class: "lbl" + (d.atlas ? " atlas" : "") }, g0, d.short);
    if (t.getComputedTextLength() > labW - 4) {
      let s = d.short;
      while (s.length > 4 && t.getComputedTextLength() > labW - 4) { s = s.slice(0, -1); t.textContent = s + "…"; }
    }
    const p = { d, px: X(d.pull), py: cy, rc: d.atlas ? 8.5 : 6.3 };
    const g = el("g", {}, g0);
    marker(g, d, p.px, p.py);
    focusable(g, p);
    el("text", { x: W - m.r + 8, y: cy + 4, class: "tick" }, g0, pullText(d.pull, 1));
    pts.push(p);
  }
  live = { pts, W, H, mode: "rows", rowH, ring: el("circle", { class: "ring", r: 10, visibility: "hidden" }, svg) };
}

function summary() {
  const ds = shown(), n = ds.length, s = document.getElementById("summary");
  if (!n) { s.textContent = "No tests shown; switch a family back on."; return; }
  const in1 = ds.filter(d => Math.abs(d.pull) <= 1).length, in2 = ds.filter(d => Math.abs(d.pull) <= 2).length;
  const worst = ds.reduce((a, b) => (Math.abs(b.pull) > Math.abs(a.pull) ? b : a));
  s.textContent = `${n} tests shown · ${in1} within ±1σ of GR (${Math.round(0.6827 * n)} expected by chance) · ` +
    `${in2} within ±2σ (${Math.round(0.9545 * n)} expected) · largest: ${worst.short}, ${pullText(worst.pull, 1)}`;
}

function render() {
  const svg = document.getElementById("chart"), W = Math.max(300, Math.floor(document.getElementById("chartbox").clientWidth));
  setHover(null);
  svg.innerHTML = "";
  (state.view === "precision" ? drawPrecision : drawPulls)(svg, W);
  summary();
}

// --- controls, table ---------------------------------------------------------------------------------------
function controls() {
  const box = document.getElementById("filters");
  for (const k of ORDER) {
    const b = document.createElement("button");
    b.className = "chip"; b.setAttribute("aria-pressed", "true"); b.title = "show or hide this family";
    b.appendChild(glyph({ family: k, atlas: "" }));
    b.appendChild(document.createTextNode(FAM[k].label));
    b.onclick = () => {
      if (state.fams.has(k)) state.fams.delete(k); else state.fams.add(k);
      b.setAttribute("aria-pressed", String(state.fams.has(k))); render(); table();
    };
    box.appendChild(b);
  }
  const st = document.createElement("span");
  st.className = "chip legend-star";
  st.appendChild(glyph({ family: "equivalence", atlas: "x" }));
  st.querySelector("polygon").setAttribute("fill", "var(--surface)");
  st.appendChild(document.createTextNode("measured in this atlas"));
  box.appendChild(st);
  const sp = document.createElement("span"); sp.className = "sep"; box.appendChild(sp);
  for (const [v, lab] of [["precision", "precision vs gravity"], ["pulls", "report card"]]) {
    const b = document.createElement("button");
    b.className = "seg"; b.textContent = lab; b.dataset.view = v; b.setAttribute("aria-pressed", String(state.view === v));
    b.onclick = () => {
      state.view = v;
      box.querySelectorAll("button[data-view]").forEach(x => x.setAttribute("aria-pressed", String(x === b)));
      render();
    };
    box.appendChild(b);
  }
  const th = document.createElement("button");
  th.className = "seg"; th.textContent = "light / dark"; th.title = "switch colour theme";
  th.onclick = () => {
    const r = document.documentElement;
    const dark = r.dataset.theme ? r.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    r.dataset.theme = dark ? "light" : "dark";
    try { localStorage.setItem("gravity-scorecard-theme", r.dataset.theme); } catch (e) {}
  };
  box.appendChild(th);
}
function table() {
  const cols = [["label", "test"], ["family", "family"], ["year", "year"], ["deviation", "deviation"], ["sigma", "1σ"],
                ["pull", "pull"], ["x", "GM/(rc²)"], ["ref", "reference"]];
  const t = document.getElementById("table");
  t.innerHTML = "";
  const head = t.createTHead().insertRow();
  for (const [k, lab] of cols) {
    const th = document.createElement("th");
    th.textContent = lab + (state.sort.key === k ? (state.sort.dir > 0 ? " ▲" : " ▼") : "");
    th.tabIndex = 0;
    th.setAttribute("aria-sort", state.sort.key === k ? (state.sort.dir > 0 ? "ascending" : "descending") : "none");
    const go = () => { state.sort = { key: k, dir: state.sort.key === k ? -state.sort.dir : 1 }; table(); };
    th.onclick = go; th.onkeydown = e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); go(); } };
    head.appendChild(th);
  }
  const key = state.sort.key, dir = state.sort.dir;
  const val = d => (key === "family" ? ORDER.indexOf(d.family) * 1e3 + Math.log10(d.x) : key === "pull" ? d.pull : d[key]);
  const rows = shown().slice().sort((a, b) => {
    const A = val(a), B = val(b);
    return (typeof A === "string" ? A.localeCompare(B) : A - B) * dir;
  });
  const body = t.createTBody();
  for (const d of rows) {
    const r = body.insertRow();
    const c0 = r.insertCell(); c0.className = "name"; c0.appendChild(glyph(d)); c0.appendChild(document.createTextNode(d.label));
    r.insertCell().textContent = FAM[d.family].label;
    r.insertCell().textContent = d.year;
    for (const v of [fmt(d.deviation), fmt(d.sigma), pullText(d.pull, 2), fmt(d.x)]) { const c = r.insertCell(); c.className = "num"; c.textContent = v; }
    const cr = r.insertCell(); cr.className = "ref"; cr.textContent = d.ref;
  }
}

try { const s = localStorage.getItem("gravity-scorecard-theme"); if (s === "light" || s === "dark") document.documentElement.dataset.theme = s; } catch (e) {}
controls(); render(); table();
const svgEl = document.getElementById("chart");
svgEl.addEventListener("pointermove", ev => { if (ev.pointerType !== "touch") setHover(nearest(ev)); });
svgEl.addEventListener("pointerdown", ev => setHover(nearest(ev)));
svgEl.addEventListener("pointerleave", ev => { if (ev.pointerType !== "touch" && !svgEl.contains(document.activeElement)) setHover(null); });
document.addEventListener("pointerdown", ev => { if (!svgEl.contains(ev.target)) setHover(null); });
document.addEventListener("keydown", ev => { if (ev.key === "Escape") setHover(null); });
window.addEventListener("scroll", () => setHover(null), { passive: true });
let lastW = document.getElementById("chartbox").clientWidth;
window.addEventListener("resize", () => {
  const w = document.getElementById("chartbox").clientWidth;
  if (w !== lastW) { lastW = w; render(); }
});
</script>
</body>
</html>
"""


def _text(v) -> str:
    return "" if v is None or (isinstance(v, float) and math.isnan(v)) else str(v)


def records(table: pd.DataFrame) -> list[dict]:
    order = {f: i for i, f in enumerate(FAMILIES)}
    t = table.assign(_o=table["family"].map(order)).sort_values(["_o", "x"]).reset_index(drop=True)
    return [{"label": _text(r.label), "short": _text(r.short), "family": r.family, "year": int(r.year),
             "ref": _text(r.ref), "quantity": _text(r.quantity), "deviation": float(r.deviation),
             "sigma": float(r.sigma), "pull": float(r.pull), "x": float(r.x), "atlas": _text(r.atlas),
             "note": _text(r.note), "source": _text(r.source)} for r in t.itertuples()]


def write(table: pd.DataFrame, path=None):
    recs = records(table)
    unknown = sorted(set(r["family"] for r in recs) - set(FAMILIES))
    if unknown:
        raise ValueError(f"families without a style: {unknown}")

    def js(o):                       # JSON that is also safe inside <script>
        return json.dumps(o, ensure_ascii=False).replace("</", "<\\/")

    html = (TEMPLATE.replace("__DATA__", js(recs)).replace("__FAMILIES__", js(FAMILIES))
            .replace("__PRIORITY__", js(PRIORITY)).replace("__N__", str(len(recs))))
    path = path or C.FIGURE_DIR / "scorecard_explorer.html"
    path.write_text(html)
    return path


if __name__ == "__main__":
    print(write(pd.read_csv(C.DERIVED_DIR / "scorecard_table.csv")))
