"""An offline page for folding the velocities yourself: figures/rv_explorer.html.

Pick a star, click any peak in its periodogram, and the velocities fold on that
period: the true period gives a clean curve, an alias gives a mess. The page
is self-contained (no network, no libraries): data embedded as JSON, drawn as
SVG, light and dark themes.

Datasets (each with the model's zero points removed):
  51 Peg      every velocity from four spectrographs
  55 Cnc e    Keck + Lick with planets b, c, d, f and the cycle term removed
  Proxima b   HARPS with the rotation signal removed
  GJ 581 'd'  HARPS with planets b, c and e removed: the activity signal
"""

from __future__ import annotations

import json

import numpy as np

import rv_physics as phys
from config import FIGURE_DIR

N_GRID = 5000  # periodogram points kept per star (max-pooled on a log-period grid)


def _pool(freqs, power, pmin, pmax, n=N_GRID):
    """Max-pool a fine periodogram onto n log-spaced period bins, so no peak
    is lost while the page stays small."""
    per = 1 / freqs
    edges = np.logspace(np.log10(pmin), np.log10(pmax), n + 1)
    idx = np.digitize(per, edges) - 1
    ok = (idx >= 0) & (idx < n)
    best = np.full(n, -1.0)
    bper = np.full(n, np.nan)
    order = np.argsort(power[ok])
    ii = idx[ok][order]
    best[ii] = power[ok][order]
    bper[ii] = per[ok][order]
    keep = best >= 0
    return np.round(bper[keep], 7).tolist(), np.round(best[keep], 4).tolist()


def _dataset(name, blurb, data, y, err, freqs, power, suggestions, model=None, unit_note=""):
    fam = [n.split()[0] for n in data.inst_names]
    t0 = float(np.floor(data.t.min()))
    per, pw = _pool(freqs, power, max(1 / freqs.max(), 0.5), min(1 / freqs.min(), 3e4))
    return {
        "name": name, "blurb": blurb, "t0": t0,
        "t": np.round(data.t - t0, 6).tolist(), "y": np.round(y, 3).tolist(), "e": np.round(err, 3).tolist(),
        "fam": [fam[k] for k in data.inst], "per": per, "pow": pw,
        "suggest": [{"P": float(p), "label": lab} for p, lab in suggestions],
        "model": model, "note": unit_note,
    }


def _model_curve(fit, j, P):
    """The fitted Keplerian of signal j over one orbit, for overlay."""
    Pa, K, e, w, tc = (a[0] for a in fit.model.planet_arrays(fit.theta))
    ph = np.linspace(0, 1, 201)
    v = fit.model.keplerian(fit.theta, tc[j] + ph * Pa[j], planets=[j])[0]
    return {"P": float(Pa[j]), "tc": float(tc[j]), "phase": np.round(ph, 4).tolist(), "v": np.round(v, 3).tolist()}


def build(results, arts, path=None):
    path = path or FIGURE_DIR / "rv_explorer.html"
    sets = []
    # 51 Peg: everything
    a = arts["51peg"]; fit = a["fit"]; d = a["data"]
    y = d.y - fit.model.offsets_and_trend(fit.theta)
    err = np.sqrt(d.err ** 2 + fit.model.jitter(fit.theta)[0, d.inst] ** 2)
    s0 = a["search"][0]
    P = fit.planets()[0]["P"]
    al = results["51peg"]["aliases_predicted"]
    sets.append(_dataset(
        "51 Peg b", "788 velocities of 51 Pegasi from ELODIE, Lick, Keck and HARPS, 1994–2018, each "
        "instrument's zero point removed. Click the tallest peak, then its echoes.",
        d, y, err, s0.freqs, s0.power,
        [(P, "the planet"), (al["sidereal_day"][0], "1-day alias"), (al["year"][0], "yearly alias"),
         (al["sidereal_day"][1], "other 1-day alias")], _model_curve(fit, 0, P)))
    # 55 Cnc e
    a = arts["55cnc"]; fit = a["fit"]; d = a["data"]
    y = d.y - fit.model.offsets_and_trend(fit.theta) - fit.model.keplerian(fit.theta, planets=[0, 1, 2, 4, 5])[0]
    err = fit.model.effective_errors(fit.theta)
    Pe = fit.planets()[3]["P"]
    alias = results["55cnc"]["e_step"]["best_period_if_search_stops_at_1d"]
    sets.append(_dataset(
        "55 Cnc e", "Keck and Lick velocities of 55 Cancri, 1989–2014, with planets b, c, d, f and the "
        "magnetic-cycle term removed. Two peaks one day apart: fold on each.",
        d, y, err, a["freqs"], a["no_e_power"], [(Pe, "the planet (0.74 d)"), (alias, "its alias (2.8 d)")],
        _model_curve(fit, 3, Pe)))
    # Proxima b
    a = arts["proxima"]; fit = a["fit"]; d = a["data"]
    y = fit.model.residuals(fit.theta) + fit.model.keplerian(fit.theta, planets=[0])[0]
    err = fit.model.effective_errors(fit.theta)
    s0 = a["search"][0]
    Pb = fit.planets()[0]["P"]
    sets.append(_dataset(
        "Proxima Cen b", "HARPS velocities of our nearest stellar neighbour, 2004–2019, with its 92-day "
        "rotation signal removed. A 1.2 m/s wobble under ~2 m/s of noise.",
        d, y, err, s0.freqs, s0.power, [(Pb, "the planet"), (1 / abs(1 / phys.SIDEREAL_DAY_D - 1 / Pb), "1-day alias")],
        _model_curve(fit, 0, Pb)))
    # GJ 581 'd'
    a = arts["gj581"]; three = a["three"]; d = a["data"]
    y = three.model.residuals(three.theta)
    err = three.model.effective_errors(three.theta)
    import rv_fit
    pw, _, _ = rv_fit.residual_periodogram(three, a["freqs"])
    Pd = results["gj581"]["signal_d"]["P"]
    sets.append(_dataset(
        "GJ 581 'd'", "HARPS velocities of Gliese 581, 2004–2012, with its three real planets removed. The "
        "66.8-day signal once published as a habitable-zone planet is the star's activity at half its "
        "rotation period.", d, y, err, a["freqs"], pw,
        [(Pd, "'planet d'"), (1 / (1 / Pd - 1 / phys.YEAR_D), "yearly alias (first published)"),
         (results["gj581"]["halpha_rotation"]["period"], "the star's rotation")]))
    payload = json.dumps({"sets": sets}, separators=(",", ":"))
    html = TEMPLATE.replace("__DATA__", payload)
    path.write_text(html, encoding="utf-8")
    return [str(path)]


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Fold It Yourself</title>
<style>
:root {
  --bg: #fcfcfb; --panel: #ffffff; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --signal: #2a78d6; --ctx: #d6d5ce; --model: #0b0b0b;
  --ELODIE: #4a3aa7; --Lick: #eb6834; --HIRES: #1baf7a; --HARPS: #2a78d6; --bin: #0b0b0b;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
    --grid: #2c2c2a; --axis: #4a4945; --signal: #3987e5; --ctx: #3a3a37; --model: #f2f1ec;
    --ELODIE: #9085e9; --Lick: #d95926; --HIRES: #199e70; --HARPS: #3987e5; --bin: #f2f1ec;
  }
}
:root[data-theme="dark"] {
  --bg: #141413; --panel: #1b1b1a; --ink: #f2f1ec; --ink2: #c3c2b7; --muted: #8f8d86;
  --grid: #2c2c2a; --axis: #4a4945; --signal: #3987e5; --ctx: #3a3a37; --model: #f2f1ec;
  --ELODIE: #9085e9; --Lick: #d95926; --HIRES: #199e70; --HARPS: #3987e5; --bin: #f2f1ec;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
.wrap { max-width: 1180px; margin: 0 auto; padding: 18px 16px 28px; }
h1 { font-size: 20px; margin: 0 0 2px; }
.sub { color: var(--ink2); margin: 0 0 10px; max-width: 920px; }
.controls { display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: center; margin: 8px 0 10px; }
.controls button { font: inherit; color: var(--ink); background: var(--panel); border: 1px solid var(--axis);
  border-radius: 6px; padding: 4px 10px; cursor: pointer; }
.controls button[aria-pressed="true"] { border-color: var(--ink); font-weight: 600; }
.controls label { color: var(--ink2); display: inline-flex; gap: 6px; align-items: center; }
#fine { width: min(360px, 80vw); }
.grid { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1fr); gap: 14px; }
@media (max-width: 860px) { .grid { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border: 1px solid var(--grid); border-radius: 8px; padding: 8px 10px; }
.card h2 { font-size: 13px; margin: 2px 0 4px; color: var(--ink2); font-weight: 600; }
.wide { grid-column: 1 / -1; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
svg text { fill: var(--muted); font-size: 11px; }
.readout { font-variant-numeric: tabular-nums; color: var(--ink); }
.readout b { font-weight: 600; }
.legend { display: flex; flex-wrap: wrap; gap: 4px 14px; color: var(--ink2); font-size: 12px; margin: 2px 0 4px; }
.sw { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 5px; vertical-align: -1px; }
.note { color: var(--muted); font-size: 12px; margin-top: 12px; max-width: 1000px; }
#blurb { color: var(--ink2); margin: 0 0 6px; }
#tip { position: fixed; pointer-events: none; background: var(--panel); color: var(--ink);
  border: 1px solid var(--axis); border-radius: 6px; padding: 4px 8px; font-size: 12px; display: none;
  box-shadow: 0 2px 8px rgba(0,0,0,.12); font-variant-numeric: tabular-nums; }
</style>
</head>
<body>
<div class="wrap">
<h1>Fold it yourself</h1>
<p class="sub">A planet makes its star's velocity repeat every orbit. Pick a star, then click a peak in the
periodogram (or a suggestion) to fold the velocities on that period. The right period stacks every orbit into one
clean curve; an alias, produced by observing only at night, leaves a mess.</p>
<div class="controls" id="stars" role="group" aria-label="dataset"></div>
<p id="blurb"></p>
<div class="controls" id="sugg" role="group" aria-label="suggested periods"></div>
<div class="grid">
  <div class="card wide">
    <h2>Periodogram: how much of the scatter a repeating signal at each period explains (click to fold)</h2>
    <svg id="pg" viewBox="0 0 1100 230" role="img" aria-label="periodogram"></svg>
  </div>
  <div class="card">
    <h2 id="foldTitle">Folded velocities</h2>
    <div class="legend" id="legend"></div>
    <svg id="fold" viewBox="0 0 540 330" role="img" aria-label="folded velocities"></svg>
    <div class="controls"><label>fine-tune <input id="fine" type="range" min="-1" max="1" step="0.001" value="0"></label></div>
  </div>
  <div class="card">
    <h2>The same velocities against time</h2>
    <svg id="ts" viewBox="0 0 540 330" role="img" aria-label="velocities against time"></svg>
    <p class="readout" id="readout"></p>
  </div>
</div>
<p class="note">Velocities from ELODIE (OHP archive), Lick (Fischer et al. 2014), Keck/HIRES (Butler et al. 2017)
and HARPS (RVBank, Trifonov et al. 2020), with the fitted zero point of each instrument removed. Large circles:
averages in 12 phase bins. "Explained" is the fraction of the variance a sine wave at the fold period removes (the
periodogram's power, recomputed live). Built by project A3 of the astrophysics atlas; everything here runs offline.</p>
</div>
<div id="tip"></div>
<script>
const D = __DATA__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const svgNS = "http://www.w3.org/2000/svg";
const el = (tag, attrs, parent) => { const e = document.createElementNS(svgNS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; };
const tip = document.getElementById("tip");
function showTip(ev, html) { tip.innerHTML = html; tip.style.display = "block";
  tip.style.left = (ev.clientX + 14) + "px"; tip.style.top = (ev.clientY + 12) + "px"; }
function hideTip() { tip.style.display = "none"; }
let cur = 0, P0 = null;

function niceTicks(lo, hi, n) { const span = hi - lo, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag; const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10)); return out; }

function axes(svg, W, H, m, xl, yl, xticks, yticks, sx, sy, xfmt, yfmt) {
  xticks.forEach(v => { const x = sx(v); el("line", {x1: x, x2: x, y1: m.t, y2: H - m.b, stroke: css("--grid")}, svg);
    const t = el("text", {x: x, y: H - m.b + 14, "text-anchor": "middle"}, svg); t.textContent = xfmt(v); });
  yticks.forEach(v => { const y = sy(v); el("line", {x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: css("--grid")}, svg);
    const t = el("text", {x: m.l - 6, y: y + 4, "text-anchor": "end"}, svg); t.textContent = yfmt(v); });
  el("line", {x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: css("--axis")}, svg);
  el("line", {x1: m.l, x2: m.l, y1: m.t, y2: H - m.b, stroke: css("--axis")}, svg);
  let t = el("text", {x: (m.l + W - m.r) / 2, y: H - 2, "text-anchor": "middle"}, svg); t.textContent = xl;
  t = el("text", {x: 12, y: (m.t + H - m.b) / 2, transform: `rotate(-90 12 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle"}, svg);
  t.textContent = yl;
}

function gls(t, y, e, P) { // power of a floating-mean sinusoid fit at period P
  const w = e.map(v => 1 / (v * v)), W = w.reduce((a, b) => a + b, 0);
  let Y = 0, C = 0, S = 0, YY = 0, YC = 0, YS = 0, CC = 0, SS = 0, CS = 0;
  for (let i = 0; i < t.length; i++) { const a = 2 * Math.PI * t[i] / P, c = Math.cos(a), s = Math.sin(a), wi = w[i] / W;
    Y += wi * y[i]; C += wi * c; S += wi * s; YY += wi * y[i] * y[i]; YC += wi * y[i] * c; YS += wi * y[i] * s;
    CC += wi * c * c; SS += wi * s * s; CS += wi * c * s; }
  YY -= Y * Y; YC -= Y * C; YS -= Y * S; CC -= C * C; SS -= S * S; CS -= C * S; const Dd = CC * SS - CS * CS;
  return (SS * YC * YC + CC * YS * YS - 2 * CS * YC * YS) / (YY * Dd);
}

function drawPG() {
  const d = D.sets[cur], svg = document.getElementById("pg"); svg.innerHTML = "";
  const W = 1100, H = 230, m = {l: 46, r: 12, t: 10, b: 32};
  const lx0 = Math.log10(Math.min(...d.per)), lx1 = Math.log10(Math.max(...d.per)), ymax = Math.max(...d.pow) * 1.1;
  const sx = p => m.l + (Math.log10(p) - lx0) / (lx1 - lx0) * (W - m.l - m.r), sy = v => H - m.b - v / ymax * (H - m.t - m.b);
  const xt = [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000].filter(v => v >= 10 ** lx0 && v <= 10 ** lx1);
  axes(svg, W, H, m, "period (days, log scale)", "power", xt, niceTicks(0, ymax, 4), sx, sy, v => v, v => v.toFixed(2));
  let dstr = ""; d.per.forEach((p, i) => { dstr += (i ? "L" : "M") + sx(p).toFixed(1) + "," + sy(d.pow[i]).toFixed(1); });
  el("path", {d: dstr, fill: "none", stroke: css("--signal"), "stroke-width": 1}, svg);
  if (P0) { const x = sx(P0); el("line", {x1: x, x2: x, y1: m.t, y2: H - m.b, stroke: css("--ink"), "stroke-width": 1.2}, svg); }
  const hit = el("rect", {x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent", style: "cursor:crosshair"}, svg);
  const toP = ev => { const r = svg.getBoundingClientRect(), x = (ev.clientX - r.left) / r.width * W;
    return 10 ** (lx0 + (x - m.l) / (W - m.l - m.r) * (lx1 - lx0)); };
  const nearest = P => { // snap to the highest pooled point within +-0.4% in log period
    let best = -1, bi = -1; d.per.forEach((p, i) => { if (Math.abs(Math.log(p / P)) < 0.004 && d.pow[i] > best) { best = d.pow[i]; bi = i; } });
    return bi; };
  hit.addEventListener("mousemove", ev => { const i = nearest(toP(ev)); if (i < 0) return hideTip();
    showTip(ev, `period <b>${d.per[i].toPrecision(6)} d</b><br>power ${d.pow[i].toFixed(3)}`); });
  hit.addEventListener("mouseleave", hideTip);
  hit.addEventListener("click", ev => { const i = nearest(toP(ev)); if (i >= 0) setP(d.per[i]); });
}

function famColor(f) { return css("--" + f) || css("--signal"); }

function drawFold() {
  const d = D.sets[cur], P = P0 * (1 + parseFloat(document.getElementById("fine").value) * 0.002);
  const svg = document.getElementById("fold"); svg.innerHTML = "";
  const W = 540, H = 330, m = {l: 46, r: 10, t: 10, b: 32};
  const ys = d.y.slice().sort((a, b) => a - b), lo = ys[Math.floor(0.01 * ys.length)], hi = ys[Math.floor(0.99 * ys.length) - 1];
  const pad = 0.15 * (hi - lo), y0 = lo - pad, y1 = hi + pad;
  const sx = p => m.l + p * (W - m.l - m.r), sy = v => H - m.b - (v - y0) / (y1 - y0) * (H - m.t - m.b);
  axes(svg, W, H, m, "phase (fraction of one orbit)", "velocity (m/s)", [0, 0.25, 0.5, 0.75, 1], niceTicks(y0, y1, 5), sx, sy,
       v => v, v => v);
  const ref = d.model ? d.model.tc - d.t0 : 0;
  const ph = d.t.map(t => ((t - ref) / P % 1 + 1) % 1);
  const bins = Array.from({length: 12}, () => [0, 0]);
  d.t.forEach((t, i) => { if (d.y[i] < y0 || d.y[i] > y1) return;
    const c = el("circle", {cx: sx(ph[i]), cy: sy(d.y[i]), r: 2.3, fill: famColor(d.fam[i]), "fill-opacity": 0.55}, svg);
    c.addEventListener("mousemove", ev => showTip(ev, `${d.fam[i]}<br>phase ${ph[i].toFixed(3)}<br>${d.y[i].toFixed(1)} ± ${d.e[i].toFixed(1)} m/s`));
    c.addEventListener("mouseleave", hideTip);
    const b = Math.min(11, Math.floor(ph[i] * 12)), w = 1 / (d.e[i] * d.e[i]); bins[b][0] += w * d.y[i]; bins[b][1] += w; });
  if (d.model && Math.abs(P / d.model.P - 1) < 2e-4) {
    let s = ""; d.model.phase.forEach((p, i) => { s += (i ? "L" : "M") + sx(p).toFixed(1) + "," + sy(d.model.v[i]).toFixed(1); });
    el("path", {d: s, fill: "none", stroke: css("--model"), "stroke-width": 1.6}, svg); }
  bins.forEach((b, k) => { if (b[1] > 0) el("circle", {cx: sx((k + 0.5) / 12), cy: sy(b[0] / b[1]), r: 4.5, fill: css("--bin"),
    stroke: css("--panel"), "stroke-width": 2}, svg); });
  const pw = gls(d.t, d.y, d.e, P);
  document.getElementById("foldTitle").textContent = `Folded on ${P.toPrecision(7)} days`;
  document.getElementById("readout").innerHTML =
    `fold period <b>${P.toPrecision(7)} d</b> · a sine wave at this period explains <b>${(100 * pw).toFixed(1)}%</b> of the variance` +
    (d.note ? `<br>${d.note}` : "");
}

function drawTS() {
  const d = D.sets[cur], svg = document.getElementById("ts"); svg.innerHTML = "";
  const W = 540, H = 330, m = {l: 46, r: 10, t: 10, b: 32};
  const ys = d.y.slice().sort((a, b) => a - b), lo = ys[Math.floor(0.01 * ys.length)], hi = ys[Math.floor(0.99 * ys.length) - 1];
  const pad = 0.15 * (hi - lo), y0 = lo - pad, y1 = hi + pad;
  const yr = d.t.map(t => 2000 + (t + d.t0 - 2451544.5) / 365.25), x0 = Math.floor(Math.min(...yr)), x1 = Math.ceil(Math.max(...yr));
  const sx = v => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r), sy = v => H - m.b - (v - y0) / (y1 - y0) * (H - m.t - m.b);
  axes(svg, W, H, m, "year", "velocity (m/s)", niceTicks(x0, x1, 6), niceTicks(y0, y1, 5), sx, sy, v => v, v => v);
  d.t.forEach((t, i) => { if (d.y[i] < y0 || d.y[i] > y1) return;
    el("circle", {cx: sx(yr[i]), cy: sy(d.y[i]), r: 2, fill: famColor(d.fam[i]), "fill-opacity": 0.6}, svg); });
}

function legend() {
  const d = D.sets[cur], fams = [...new Set(d.fam)], names = {ELODIE: "ELODIE", Lick: "Lick", HIRES: "Keck/HIRES", HARPS: "HARPS"};
  document.getElementById("legend").innerHTML = fams.map(f =>
    `<span><span class="sw" style="background:var(--${f})"></span>${names[f] || f}</span>`).join("") +
    `<span><span class="sw" style="background:var(--bin)"></span>phase-bin average</span>`;
}

function setP(P) { P0 = P; document.getElementById("fine").value = 0; drawPG(); drawFold();
  document.querySelectorAll("#sugg button").forEach(b => b.setAttribute("aria-pressed", Math.abs(+b.dataset.p / P - 1) < 1e-6)); }

function setStar(k) { cur = k; const d = D.sets[k];
  document.querySelectorAll("#stars button").forEach((b, i) => b.setAttribute("aria-pressed", i === k));
  document.getElementById("blurb").textContent = d.blurb;
  document.getElementById("sugg").innerHTML = "<span style='color:var(--ink2)'>try:</span>" + d.suggest.map(s =>
    `<button data-p="${s.P}">${s.label} · ${s.P.toPrecision(5)} d</button>`).join("");
  document.querySelectorAll("#sugg button").forEach(b => b.addEventListener("click", () => setP(+b.dataset.p)));
  legend(); drawTS(); setP(d.suggest[0].P); }

document.getElementById("stars").innerHTML = D.sets.map((s, i) => `<button>${s.name}</button>`).join("");
document.querySelectorAll("#stars button").forEach((b, i) => b.addEventListener("click", () => setStar(i)));
document.getElementById("fine").addEventListener("input", drawFold);
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { drawPG(); drawFold(); drawTS(); });
setStar(0);
</script>
</body>
</html>
"""
