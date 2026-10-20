"""A self-contained HTML explorer: turn the cosmology knobs, watch the counts.

Two sliders (Omega_m, sigma_8) or a click on the parameter plane choose a
universe; the page shows how many clusters eRASS1 should then have found, by
redshift and by brightness, against what it did find. Predictions are
precomputed on a grid with the published count rate-mass relation (the
consistency-scan model) and embedded, so the file works offline. No external
libraries: plain SVG drawn by a few lines of JavaScript.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import cosmology as CO
import model as MD

OM_STEPS = np.round(np.arange(0.20, 0.4201, 0.02), 3)
S8_STEPS = np.round(np.arange(0.70, 1.0001, 0.02), 3)


def _contours(grid):
    """68/95% contour polylines of the consistency-scan posterior."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = np.exp(grid["lnl"] - grid["lnl"].max())
    flat = np.sort(p.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    lev = [flat[np.searchsorted(cum, f)] for f in (0.95, 0.68)]
    fig, ax = plt.subplots()
    cs = ax.contour(grid["om"], grid["s8"], p.T, levels=lev)
    out = []
    for level_segs in cs.allsegs:
        out.append([np.round(s, 4).tolist() for s in level_segs])
    plt.close(fig)
    return out


def _payload(res) -> dict:
    cm = MD.CountsModel()
    pred = {}
    for om in OM_STEPS:
        for s8 in S8_STEPS:
            N = cm.expected(CO.Cosmology(Om=float(om), sigma8=float(s8)))
            pred[f"{om:.2f}_{s8:.2f}"] = [np.round(N.sum(0), 1).tolist(), np.round(N.sum(1), 1).tolist(),
                                          round(float(MD.loglike(N, res["counts"]["data"])), 2)]
    data = res["counts"]["data"]
    g = res["grid"]
    return {
        "om": OM_STEPS.tolist(), "s8": S8_STEPS.tolist(), "pred": pred,
        "data_z": np.round(data.sum(0), 1).tolist(), "data_cr": np.round(data.sum(1), 1).tolist(),
        "z_edges": np.round(res["counts"]["z_edges"], 3).tolist(),
        "cr_edges": np.round(np.exp(res["counts"]["cr_edges"]), 5).tolist(),
        "best": g["best"], "contours": _contours(g),
        "lnl_best": float(g["lnl"].max()),
        "planck": {"om": 0.3153, "s8": 0.8111}, "ghirardini": {"om": 0.29, "s8": 0.88},
    }


HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cluster Count Explorer</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --wash: #ecebe6;
  --c1: #2a78d6; --c2: #eb6834; --c3: #1baf7a; --ctx: #b9b8b1; --c1soft: #b7d3f6;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
    --grid: #33322f; --axis: #52514e; --wash: #2b2a28;
    --c1: #4c90e6; --c2: #f07d4f; --c3: #2fc28b; --ctx: #6f6e69; --c1soft: #1c3d66;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
  --grid: #33322f; --axis: #52514e; --wash: #2b2a28;
  --c1: #4c90e6; --c2: #f07d4f; --c3: #2fc28b; --ctx: #6f6e69; --c1soft: #1c3d66;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1120px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 20px; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 16px; max-width: 760px; }
.controls { display: flex; flex-wrap: wrap; gap: 18px 28px; align-items: center; background: var(--panel);
  border-radius: 10px; padding: 12px 16px; margin-bottom: 14px; }
.controls label { display: flex; gap: 10px; align-items: center; color: var(--ink2); }
.controls input[type=range] { width: 200px; accent-color: var(--c1); }
.val { font-variant-numeric: tabular-nums; color: var(--ink); min-width: 3.2em; display: inline-block; }
.tiles { display: flex; flex-wrap: wrap; gap: 10px; margin-bottom: 14px; }
.tile { background: var(--panel); border-radius: 10px; padding: 10px 14px; min-width: 150px; }
.tile .k { font-size: 12px; color: var(--muted); }
.tile .v { font-size: 22px; font-weight: 600; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 14px; }
.card { background: var(--panel); border-radius: 10px; padding: 12px 12px 6px; }
.card h2 { font-size: 14px; margin: 0 0 6px; }
svg { width: 100%; height: auto; display: block; }
.legend { display: flex; gap: 14px; flex-wrap: wrap; font-size: 12px; color: var(--ink2); margin: 2px 0 6px; }
.sw { display: inline-block; width: 14px; height: 3px; vertical-align: middle; margin-right: 5px; border-radius: 2px; }
.dot { display: inline-block; width: 8px; height: 8px; border-radius: 50%; vertical-align: middle; margin-right: 5px; }
#tip { position: fixed; pointer-events: none; background: var(--surface); color: var(--ink); border: 1px solid var(--axis);
  border-radius: 6px; padding: 6px 8px; font-size: 12px; display: none; box-shadow: 0 2px 8px rgba(0,0,0,.15); }
button { font: inherit; background: var(--surface); color: var(--ink); border: 1px solid var(--axis);
  border-radius: 6px; padding: 4px 10px; cursor: pointer; }
footer { color: var(--muted); font-size: 12px; margin-top: 16px; max-width: 900px; }
table { border-collapse: collapse; font-size: 12px; margin-top: 8px; font-variant-numeric: tabular-nums; }
td, th { padding: 2px 8px; border-bottom: 1px solid var(--grid); text-align: right; }
details { margin-top: 10px; color: var(--ink2); }
</style>
</head>
<body>
<main>
<h1>How clumpy is the universe? Count the clusters.</h1>
<p class="lede">eROSITA found about 5,000 galaxy clusters in half the sky. How many it <em>should</em> have found depends on
how much matter the universe holds (&Omega;<sub>m</sub>) and how clumpy it has become (&sigma;<sub>8</sub>). Move the
sliders, or click the parameter plane, and compare the prediction (blue) with the real counts (black).</p>
<div class="controls">
  <label>&Omega;<sub>m</sub> <input id="om" type="range" min="0" max="10" step="1"><span class="val" id="omv"></span></label>
  <label>&sigma;<sub>8</sub> <input id="s8" type="range" min="0" max="15" step="1"><span class="val" id="s8v"></span></label>
  <button id="best">best fit</button><button id="planck">Planck CMB</button>
  <button id="theme" aria-label="toggle colour theme">light / dark</button>
</div>
<div class="tiles">
  <div class="tile"><div class="k">clusters predicted</div><div class="v" id="npred"></div></div>
  <div class="tile"><div class="k">clusters observed (purity-weighted)</div><div class="v" id="nobs"></div></div>
  <div class="tile"><div class="k">S<sub>8</sub> = &sigma;<sub>8</sub>(&Omega;<sub>m</sub>/0.3)<sup>0.5</sup></div><div class="v" id="s8x"></div></div>
  <div class="tile"><div class="k">worse than best fit by &Delta;&chi;&sup2;</div><div class="v" id="dchi"></div></div>
</div>
<div class="grid">
  <div class="card"><h2>Parameter plane</h2>
    <div class="legend"><span><span class="sw" style="background:var(--c1)"></span>68 / 95% from the counts</span>
      <span><span class="dot" style="background:var(--c3)"></span>Planck CMB</span>
      <span><span class="dot" style="background:var(--c2)"></span>eRASS1 published</span></div>
    <svg id="plane" viewBox="0 0 360 300" role="img" aria-label="Omega_m versus sigma_8 plane"></svg></div>
  <div class="card"><h2>Clusters by redshift</h2>
    <div class="legend"><span><span class="sw" style="background:var(--c1)"></span>prediction</span>
      <span><span class="dot" style="background:var(--ink)"></span>eRASS1</span></div>
    <svg id="zplot" viewBox="0 0 360 300" role="img" aria-label="counts per redshift bin"></svg></div>
  <div class="card"><h2>Clusters by X-ray brightness</h2>
    <div class="legend"><span><span class="sw" style="background:var(--c1)"></span>prediction</span>
      <span><span class="dot" style="background:var(--ink)"></span>eRASS1</span></div>
    <svg id="crplot" viewBox="0 0 360 300" role="img" aria-label="counts per count-rate bin"></svg></div>
</div>
<details><summary>Table view (current settings)</summary><div id="tbl"></div></details>
<footer>Model: Tinker et al. (2008) mass function, the published eRASS1 count rate&ndash;mass relation (Ghirardini et al. 2024),
our reconstruction of the eRASS1 selection function and footprint; H<sub>0</sub> = 67.77, &Omega;<sub>b</sub> = 0.049,
n<sub>s</sub> = 0.96 fixed. Counts weighted by each cluster's probability of being real. Precomputed on a grid of
&Delta;&Omega;<sub>m</sub> = &Delta;&sigma;<sub>8</sub> = 0.02. AI-generated as part of a hobby astrophysics atlas (project D11).</footer>
</main>
<div id="tip"></div>
<script>
const D = __DATA__;
const $ = (id) => document.getElementById(id);
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const NS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
const tip = $("tip");
function showTip(ev, html) { tip.innerHTML = html; tip.style.display = "block";
  tip.style.left = Math.min(ev.clientX + 12, window.innerWidth - 180) + "px"; tip.style.top = (ev.clientY + 12) + "px"; }
function hideTip() { tip.style.display = "none"; }
let iom = D.om.indexOf(0.30), is8 = D.s8.indexOf(0.86);
if (iom < 0) iom = 5; if (is8 < 0) is8 = 8;
function key() { return D.om[iom].toFixed(2) + "_" + D.s8[is8].toFixed(2); }
const W = 360, H = 300, M = {l: 46, r: 12, t: 10, b: 38};
function axes(svg, xs, ys, xt, yt, xlab, ylab, xfmt, yfmt) {
  for (const v of yt) { const y = ys(v);
    el("line", {x1: M.l, x2: W - M.r, y1: y, y2: y, stroke: css("--grid"), "stroke-width": 1}, svg);
    el("text", {x: M.l - 6, y: y + 4, "text-anchor": "end", "font-size": 10, fill: css("--muted")}, svg).textContent = yfmt(v); }
  for (const v of xt) { const x = xs(v);
    el("text", {x: x, y: H - M.b + 14, "text-anchor": "middle", "font-size": 10, fill: css("--muted")}, svg).textContent = xfmt(v); }
  el("line", {x1: M.l, x2: W - M.r, y1: H - M.b, y2: H - M.b, stroke: css("--axis")}, svg);
  el("text", {x: (M.l + W - M.r) / 2, y: H - 6, "text-anchor": "middle", "font-size": 11, fill: css("--ink2")}, svg).textContent = xlab;
  const t = el("text", {x: 12, y: (M.t + H - M.b) / 2, "text-anchor": "middle", "font-size": 11, fill: css("--ink2"),
    transform: `rotate(-90 12 ${(M.t + H - M.b) / 2})`}, svg); t.textContent = ylab;
}
function histPlot(id, edges, model, data, logx, xlab) {
  const svg = $(id); svg.innerHTML = "";
  const lx = (v) => logx ? Math.log10(v) : v;
  const lo = logx ? Math.log10(0.005) : edges[0], hi = logx ? Math.log10(30) : edges[edges.length - 1];
  const xs = (v) => M.l + (lx(v) - lo) / (hi - lo) * (W - M.l - M.r);
  const ymax = Math.max(...model, ...data.map((d) => d + Math.sqrt(d))) * 1.08;
  const ys = (v) => H - M.b - v / ymax * (H - M.t - M.b);
  const step = Math.pow(10, Math.floor(Math.log10(ymax / 4)));
  const nice = [1, 2, 5, 10].map((m) => m * step).find((s) => ymax / s <= 6);
  const yt = []; for (let v = 0; v <= ymax; v += nice) yt.push(v);
  const xt = logx ? [0.01, 0.1, 1, 10] : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
  axes(svg, xs, ys, xt, yt, xlab, "clusters per bin", (v) => (logx ? (v >= 1 ? v.toFixed(0) : v.toString()) : v.toFixed(1)),
       (v) => v.toFixed(0));
  const cx = (i) => logx ? Math.sqrt(Math.max(edges[i], 0.005) * Math.min(edges[i + 1], 30))
                        : 0.5 * (edges[i] + edges[i + 1]);
  let d = "";
  for (let i = 0; i < model.length; i++) {
    const x0 = xs(Math.max(edges[i], logx ? 0.005 : edges[0])), x1 = xs(Math.min(edges[i + 1], logx ? 30 : edges[edges.length - 1]));
    d += (i ? "L" : "M") + x0 + "," + ys(model[i]) + "L" + x1 + "," + ys(model[i]);
  }
  el("path", {d: d, fill: "none", stroke: css("--c1"), "stroke-width": 2}, svg);
  for (let i = 0; i < data.length; i++) {
    const x = xs(cx(i)), e = Math.sqrt(Math.max(data[i], 0));
    el("line", {x1: x, x2: x, y1: ys(data[i] - e), y2: ys(data[i] + e), stroke: css("--ink"), "stroke-width": 1.2}, svg);
    el("circle", {cx: x, cy: ys(data[i]), r: 3.6, fill: css("--ink"), stroke: css("--surface"), "stroke-width": 1.5}, svg);
    const hit = el("rect", {x: xs(Math.max(edges[i], logx ? 0.005 : edges[0])), y: M.t,
      width: Math.max(xs(Math.min(edges[i + 1], logx ? 30 : edges[edges.length - 1])) - xs(Math.max(edges[i], logx ? 0.005 : edges[0])), 6),
      height: H - M.t - M.b, fill: "transparent"}, svg);
    const lab = logx ? `${edges[i].toPrecision(2)}&ndash;${edges[i + 1].toPrecision(2)} counts/s` : `z ${edges[i].toFixed(2)}&ndash;${edges[i + 1].toFixed(2)}`;
    hit.addEventListener("mousemove", (ev) => showTip(ev, `${lab}<br>observed ${data[i].toFixed(0)} &plusmn; ${e.toFixed(0)}<br>predicted ${model[i].toFixed(0)}`));
    hit.addEventListener("mouseleave", hideTip);
  }
}
const PL = {lo: 0.19, hi: 0.43, slo: 0.69, shi: 1.01};
function planePlot() {
  const svg = $("plane"); svg.innerHTML = "";
  const xs = (v) => M.l + (v - PL.lo) / (PL.hi - PL.lo) * (W - M.l - M.r);
  const ys = (v) => H - M.b - (v - PL.slo) / (PL.shi - PL.slo) * (H - M.t - M.b);
  axes(svg, xs, ys, [0.2, 0.25, 0.3, 0.35, 0.4], [0.7, 0.8, 0.9, 1.0], "Ωm (matter density)", "σ8 (clumpiness)",
       (v) => v.toFixed(2), (v) => v.toFixed(1));
  D.contours.forEach((segs, k) => segs.forEach((s) => {
    el("path", {d: "M" + s.map((p) => xs(p[0]) + "," + ys(p[1])).join("L"), fill: k ? css("--c1soft") : "none",
      "fill-opacity": k ? 0.9 : 0, stroke: css("--c1"), "stroke-width": k ? 1.6 : 1}, svg); }));
  for (const [p, c, lab] of [[D.planck, "--c3", "Planck CMB"], [D.ghirardini, "--c2", "eRASS1 published"]]) {
    el("circle", {cx: xs(p.om), cy: ys(p.s8), r: 5, fill: css(c), stroke: css("--surface"), "stroke-width": 2}, svg);
    el("text", {x: xs(p.om) + 8, y: ys(p.s8) + 4, "font-size": 10, fill: css("--ink2")}, svg).textContent = lab;
  }
  const x = xs(D.om[iom]), y = ys(D.s8[is8]);
  el("line", {x1: x, x2: x, y1: M.t, y2: H - M.b, stroke: css("--ink2"), "stroke-width": 0.8}, svg);
  el("line", {x1: M.l, x2: W - M.r, y1: y, y2: y, stroke: css("--ink2"), "stroke-width": 0.8}, svg);
  el("circle", {cx: x, cy: y, r: 5.5, fill: "none", stroke: css("--ink"), "stroke-width": 2}, svg);
  const hit = el("rect", {x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, fill: "transparent", cursor: "crosshair"}, svg);
  hit.addEventListener("click", (ev) => {
    const r = svg.getBoundingClientRect(), sx = W / r.width, sy = H / r.height;
    const om = PL.lo + ((ev.clientX - r.left) * sx - M.l) / (W - M.l - M.r) * (PL.hi - PL.lo);
    const s8 = PL.slo + (H - M.b - (ev.clientY - r.top) * sy) / (H - M.t - M.b) * (PL.shi - PL.slo);
    iom = nearest(D.om, om); is8 = nearest(D.s8, s8); update(); });
}
function nearest(arr, v) { let b = 0; arr.forEach((a, i) => { if (Math.abs(a - v) < Math.abs(arr[b] - v)) b = i; }); return b; }
function update() {
  $("om").value = iom; $("s8").value = is8; $("omv").textContent = D.om[iom].toFixed(2); $("s8v").textContent = D.s8[is8].toFixed(2);
  const [mz, mcr, lnl] = D.pred[key()];
  const ntot = mz.reduce((a, b) => a + b, 0), nobs = D.data_z.reduce((a, b) => a + b, 0);
  $("npred").textContent = Math.round(ntot).toLocaleString(); $("nobs").textContent = Math.round(nobs).toLocaleString();
  $("s8x").textContent = (D.s8[is8] * Math.sqrt(D.om[iom] / 0.3)).toFixed(3);
  $("dchi").textContent = Math.max(0, 2 * (D.lnl_best - lnl)).toFixed(0);
  histPlot("zplot", D.z_edges, mz, D.data_z, false, "photometric redshift");
  histPlot("crplot", D.cr_edges.map((v, i) => i === 0 ? 0.004 : (i === D.cr_edges.length - 1 ? 40 : v)), mcr, D.data_cr, true, "count rate (counts/s)");
  planePlot();
  let t = "<table><tr><th>z bin</th><th>observed</th><th>predicted</th></tr>";
  for (let i = 0; i < mz.length; i++) t += `<tr><td>${D.z_edges[i].toFixed(2)}&ndash;${D.z_edges[i + 1].toFixed(2)}</td><td>${D.data_z[i].toFixed(0)}</td><td>${mz[i].toFixed(0)}</td></tr>`;
  $("tbl").innerHTML = t + "</table>";
}
$("om").max = D.om.length - 1; $("s8").max = D.s8.length - 1;
$("om").addEventListener("input", (e) => { iom = +e.target.value; update(); });
$("s8").addEventListener("input", (e) => { is8 = +e.target.value; update(); });
$("best").addEventListener("click", () => { iom = nearest(D.om, D.best.Om); is8 = nearest(D.s8, D.best.sigma8); update(); });
$("planck").addEventListener("click", () => { iom = nearest(D.om, D.planck.om); is8 = nearest(D.s8, D.planck.s8); update(); });
$("theme").addEventListener("click", () => { const r = document.documentElement;
  const dark = r.dataset.theme ? r.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
  r.dataset.theme = dark ? "light" : "dark"; update(); });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", update);
update();
</script>
</body>
</html>
"""


def write(res) -> str:
    payload = _payload(res)
    path = C.FIGURE_DIR / "cluster_count_explorer.html"
    path.write_text(HTML.replace("__DATA__", json.dumps(payload, separators=(",", ":"))))
    return str(path)
