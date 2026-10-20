"""An offline page for building a supernova light curve yourself: figures/nickel_explorer.html.

Pick a supernova and drag three sliders: how much nickel-56 the explosion
made, how long light takes to diffuse out (Type Ia only), and when the gamma
rays start to escape (T0). The model is the one this project fits: the
56Ni -> 56Co -> 56Fe heating, Arnett's diffusion for the Type Ia peak, the
(T0/t)^2 gamma-ray leakage, positrons kept or let go. For the Type II
supernovae only the radioactive tail is modelled: their first months are
powered by the explosion's own heat, not by nickel.

The page is self-contained (no network, no libraries): data embedded as JSON,
drawn as SVG, light and dark themes, a table view of the data.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import config as C
import nickel as N


def _series(phase, lum, err, measured=None) -> dict:
    out = {"t": np.round(np.asarray(phase, float), 2).tolist(),
           "L": [float(f"{v:.4g}") for v in np.asarray(lum, float)],
           "e": [float(f"{v:.3g}") for v in np.asarray(err, float)]}
    out["m"] = (np.ones(len(out["t"]), bool) if measured is None else np.asarray(measured, bool)).astype(int).tolist()
    return out


def _payload(r: dict) -> dict:
    s = r["summary"]
    data = {}
    for sn in ("SN2011fe", "SN2014J"):
        tab = r[sn]["bolo"]["table"]
        af = r[sn]["fit"]["arnett_fit"]
        data[sn] = {"label": C.SUPERNOVAE[sn]["label"], "type": "Ia", "tmax": 360 if sn == "SN2011fe" else 140,
                    "ymin": 3e39 if sn == "SN2011fe" else 1e41, "series": _series(tab["phase"], tab["L"], tab["L_err"],
                                                                                tab["nir_source"] == "measured"),
                    "fit": {"m": round(af["m_ni"], 3), "tau": round(af["tau_m"], 2), "t0": round(af["t0"], 1)},
                    "tail_start": 0, "note": (
                        "Type Ia: a white dwarf blown apart. The peak is the nickel's heat diffusing out; "
                        "after it, the light follows the decay, then sinks as gamma rays escape."
                        if sn == "SN2011fe" else
                        "Type Ia behind thick dust in M82, dimmed 5 magnitudes in B. Its late points rise above "
                        "any model: background light and broad-filter dust effects, not physics.")}
    tab = r["SN2023ixf"]["bolo"]["table"]
    f = r["SN2023ixf"]["fit"]["tail_measured_nir"]
    data["SN2023ixf"] = {"label": "SN 2023ixf", "type": "II", "tmax": 365, "ymin": 1e40,
                         "series": _series(tab["phase"], tab["L"], tab["L_err"], tab["nir_source"] == "measured"),
                         "fit": {"m": round(f["m_ni"], 4), "tau": 0, "t0": round(f["t0"], 1)}, "tail_start": 90,
                         "note": "Type II: a red supergiant's core collapsed. The plateau is the explosion's own heat "
                                 "leaking out as hydrogen recombines; nickel only takes over after the drop near day 85."}
    sn87 = r["SN1987A"]["bolo"]
    pub = sn87["published"]
    late = sn87["published_late"][sn87["published_late"]["phase"] > pub["phase"].max()]
    both = pd.concat([pub[["phase", "L", "L_err"]], late[["phase", "L", "L_err"]]])
    f = r["SN1987A"]["fit"]["tail_with_escape"]
    data["SN1987A"] = {"label": "SN 1987A", "type": "II", "tmax": 1400, "ymin": 1e37,
                       "series": _series(both["phase"], both["L"], both["L_err"]),
                       "fit": {"m": round(f["m_ni"], 4), "tau": 0, "t0": round(f["t0"], 1)}, "tail_start": 120,
                       "note": "Type II from a compact blue supergiant: its shock heat faded fast, so nickel "
                               "powered even the slow rise to the peak at day 85. Past day ~800 cobalt-57, "
                               "titanium-44 and new dust take over (not modelled here)."}
    a, b = N.two_exponential_coefficients()
    return {"supernovae": data, "colors": {"SN2011fe": ["#2a78d6", "#3987e5"], "SN2014J": ["#4a3aa7", "#9085e9"],
                                            "SN1987A": ["#eb6834", "#d95926"], "SN2023ixf": ["#1baf7a", "#199e70"]},
            "nuclear": {"tauNi": N.TAU_NI, "tauCo": N.TAU_CO, "A": a, "B": b,
                        "posFrac": C.E_POSITRON_CO_MEV / (C.E_GAMMA_CO_MEV + C.E_POSITRON_CO_MEV)},
            "ztf_median": s["ztf"]["m_ni_median"]}


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Nickel Light Builder</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --model: #52514e;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
    --grid: #2f2f2d; --axis: #4a4a46; --model: #c3c2b7;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
  --grid: #2f2f2d; --axis: #4a4a46; --model: #c3c2b7;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1060px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 16px; max-width: 760px; }
.row { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 12px; }
button { font: inherit; font-size: 0.9rem; padding: 6px 12px; border-radius: 8px; cursor: pointer;
         border: 1px solid var(--axis); background: var(--surface); color: var(--ink); }
button[aria-pressed="true"] { background: var(--ink); color: var(--surface); border-color: var(--ink); }
button:focus-visible, input:focus-visible { outline: 2px solid #2a78d6; outline-offset: 2px; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 290px; gap: 16px; }
@media (max-width: 820px) { .layout { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border-radius: 12px; padding: 12px; }
svg { width: 100%; height: auto; display: block; }
.ctl { margin-bottom: 14px; }
.ctl label { display: flex; justify-content: space-between; font-size: 0.88rem; color: var(--ink2); }
.ctl label b { color: var(--ink); font-variant-numeric: tabular-nums; }
.ctl input[type=range] { width: 100%; }
.ctl.disabled { opacity: 0.45; }
.check { font-size: 0.88rem; color: var(--ink2); display: flex; gap: 6px; align-items: center; margin-bottom: 12px; }
.readout { font-size: 0.86rem; color: var(--ink2); border-top: 1px solid var(--grid); padding-top: 10px; }
.readout b { color: var(--ink); }
.note { font-size: 0.88rem; color: var(--ink2); margin: 10px 2px 0; }
.legend { display: flex; flex-wrap: wrap; gap: 14px; font-size: 0.82rem; color: var(--ink2); margin: 6px 2px 0; }
.legend span { display: inline-flex; align-items: center; gap: 6px; white-space: nowrap; }
#tip { position: fixed; pointer-events: none; background: var(--surface); color: var(--ink);
       border: 1px solid var(--axis); border-radius: 6px; padding: 5px 8px; font-size: 0.8rem; display: none;
       box-shadow: 0 2px 8px rgba(0,0,0,0.12); }
details { margin-top: 16px; font-size: 0.85rem; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
td, th { padding: 2px 10px; text-align: right; border-bottom: 1px solid var(--grid); }
.tablewrap { max-height: 280px; overflow: auto; }
.axis text { fill: var(--muted); font-size: 11px; }
.axis line, .axis path { stroke: var(--axis); }
.gridline { stroke: var(--grid); }
.lbl { fill: var(--ink2); font-size: 12px; }
footer { color: var(--muted); font-size: 0.78rem; margin-top: 18px; }
</style>
</head>
<body>
<main>
  <h1>Build a supernova's light from nickel</h1>
  <p class="lede">Explosions forge radioactive nickel-56, which decays to cobalt-56 and then to iron.
  The gamma rays from those decays heat the debris, and that heat is the light we see for months.
  Choose a supernova and set how much nickel it made: the model below is the decay chain,
  plus light diffusing out of the debris, plus gamma rays escaping once the debris thins.</p>
  <div class="row" id="sns" role="group" aria-label="supernova"></div>
  <div class="layout">
    <div class="card">
      <svg id="chart" viewBox="0 0 720 440" role="img" aria-label="luminosity against time, data and model"></svg>
      <div class="legend" id="legend"></div>
      <p class="note" id="note"></p>
    </div>
    <div class="card">
      <div class="ctl"><label for="m">Nickel-56 made <b id="mv"></b></label>
        <input id="m" type="range" min="-2.3" max="0.25" step="0.005"></div>
      <div class="ctl" id="taubox"><label for="tau">Diffusion time <b id="tauv"></b></label>
        <input id="tau" type="range" min="3" max="40" step="0.1"></div>
      <div class="ctl"><label for="t0">Gamma rays escape after T0 <b id="t0v"></b></label>
        <input id="t0" type="range" min="1" max="3.3" step="0.005"></div>
      <label class="check"><input id="pos" type="checkbox" checked> positrons stay in the debris</label>
      <div class="row">
        <button id="fit">project's fit</button>
        <button id="trap">trap every gamma ray</button>
      </div>
      <div class="readout" id="readout"></div>
    </div>
  </div>
  <details><summary>Table view of the data</summary><div class="tablewrap"><table id="tbl"></table></div></details>
  <footer>Data: Open Supernova Catalog photometry (SN 2011fe, SN 2014J), Hsu et al. 2025 and Singh et al. 2024
  (SN 2023ixf), published UVOIR of SN 1987A (Suntzeff et al. 1991, 1992; Bouchet et al. 1991), integrated
  over wavelength by this project. Open circles: near-infrared share estimated, not measured.
  Nuclear data: NNDC; energies per decay: Nadyozhin 1994.</footer>
</main>
<div id="tip" role="tooltip"></div>
<script>
const D = __DATA__;
const NS = "http://www.w3.org/2000/svg";
const nuc = D.nuclear;
let cur = "SN2011fe";
const st = {};

const dark = () => {
  const forced = document.documentElement.dataset.theme;
  if (forced) return forced === "dark";
  return window.matchMedia("(prefers-color-scheme: dark)").matches;
};
const color = sn => D.colors[sn][dark() ? 1 : 0];
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

// Heating per solar mass of nickel (erg/s), split into gamma rays and positrons.
function heating(t) {
  const total = nuc.A * Math.exp(-t / nuc.tauNi) + nuc.B * Math.exp(-t / nuc.tauCo);
  const co = nuc.B * (Math.exp(-t / nuc.tauCo) - Math.exp(-t / nuc.tauNi));
  const pos = co * nuc.posFrac;
  return { gamma: total - pos, pos: pos };
}
function deposit(t, T0, posTrapped) {
  const h = heating(Math.max(t, 1e-4));
  const f = T0 > 0 ? 1 - Math.exp(-((T0 / t) ** 2)) : 1;
  return h.gamma * f + (posTrapped ? h.pos : 0);
}
// Arnett (1982): L(t) = e^{-x^2} Int 2 t'/tau^2 e^{(t'/tau)^2} D(t') dt', stepped so no exponential
// overflows; after 8 diffusion times the light simply follows the heating.
function arnett(tmax, M, tau, T0, pos) {
  const out = [], dt = 0.05, tm2 = tau * tau, tSwitch = Math.min(tmax, 8 * tau);
  let L = 0, t = 0, d1 = M * deposit(1e-4, T0, pos);
  out.push([0, 0]);
  while (t < tSwitch) {
    const t2 = t + dt, d2 = M * deposit(t2, T0, pos);
    const damp = Math.exp(-(t2 * t2 - t * t) / tm2);
    L = L * damp + 0.5 * dt * (2 * t / tm2 * damp * d1 + 2 * t2 / tm2 * d2);
    t = t2; d1 = d2;
    if (Math.abs(t / 0.5 - Math.round(t / 0.5)) < 1e-6) out.push([t, L]);
  }
  for (let tt = Math.ceil(t) + 1; tt <= tmax; tt += 1) out.push([tt, M * deposit(tt, T0, pos)]);
  return out;
}

function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
const fmtL = v => { const p = Math.floor(Math.log10(v)); return (v / 10 ** p).toFixed(2) + "e" + p; };
const fmtM = v => v >= 0.1 ? v.toFixed(2) : v.toFixed(3);

function draw() {
  const sn = D.supernovae[cur], s = st[cur], svg = document.getElementById("chart");
  svg.innerHTML = "";
  // On narrow screens draw on a smaller canvas so the text stays legible.
  const narrow = svg.clientWidth > 0 && svg.clientWidth < 560;
  const W = narrow ? 460 : 720, H = narrow ? 400 : 440, ml = 62, mr = 12, mt = 14, mb = 44;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const ser = sn.series;
  const ymax = Math.max(...ser.L) * 4, ymin = sn.ymin, tmax = sn.tmax;
  const X = t => ml + (t / tmax) * (W - ml - mr);
  const Y = L => mt + (Math.log10(ymax) - Math.log10(Math.max(L, ymin * 0.5))) /
                 (Math.log10(ymax) - Math.log10(ymin)) * (H - mt - mb);
  const ax = el("g", { class: "axis" }, svg);
  for (let p = Math.ceil(Math.log10(ymin)); p <= Math.log10(ymax); p++) {
    el("line", { x1: ml, x2: W - mr, y1: Y(10 ** p), y2: Y(10 ** p), class: "gridline" }, ax);
    const tx = el("text", { x: ml - 8, y: Y(10 ** p) + 4, "text-anchor": "end" }, ax);
    tx.textContent = "1e" + p;
  }
  const step = tmax > 800 ? 200 : tmax > 200 ? 50 : 20;
  for (let t = 0; t <= tmax; t += step) {
    el("line", { x1: X(t), x2: X(t), y1: mt, y2: H - mb, class: "gridline" }, ax);
    const tx = el("text", { x: X(t), y: H - mb + 16, "text-anchor": "middle" }, ax);
    tx.textContent = t;
  }
  el("line", { x1: ml, x2: W - mr, y1: H - mb, y2: H - mb }, ax);
  el("line", { x1: ml, x2: ml, y1: mt, y2: H - mb }, ax);
  const xl = el("text", { x: (ml + W - mr) / 2, y: H - 8, "text-anchor": "middle", class: "lbl" }, svg);
  xl.textContent = sn.type === "Ia" ? "days since first light" : "days since explosion";
  const yl = el("text", { x: 14, y: (mt + H - mb) / 2, "text-anchor": "middle", class: "lbl",
                          transform: `rotate(-90 14 ${(mt + H - mb) / 2})` }, svg);
  yl.textContent = "luminosity (erg/s)";

  const clip = el("clipPath", { id: "plotclip" }, el("defs", {}, svg));
  el("rect", { x: ml, y: mt, width: W - ml - mr, height: H - mt - mb }, clip);
  const g = el("g", { "clip-path": "url(#plotclip)" }, svg);
  const c = color(cur);
  // Every gamma ray trapped (reference) and the model.
  const tStart = sn.tail_start;
  let pts = [];
  for (let t = Math.max(tStart, 0.5); t <= tmax; t += 1) pts.push([t, s.m * deposit(t, 0, true)]);
  el("path", { d: "M" + pts.map(p => X(p[0]).toFixed(1) + "," + Y(p[1]).toFixed(1)).join("L"), fill: "none",
               stroke: css("--model"), "stroke-width": 1.2, "stroke-dasharray": "5 4" }, g);
  let model;
  if (sn.type === "Ia") model = arnett(tmax, s.m, s.tau, s.t0, s.pos);
  else { model = []; for (let t = tStart; t <= tmax; t += 1) model.push([t, s.m * deposit(t, s.t0, s.pos)]); }
  model = model.filter(p => p[1] > 0);
  el("path", { d: "M" + model.map(p => X(p[0]).toFixed(1) + "," + Y(p[1]).toFixed(1)).join("L"), fill: "none",
               stroke: c, "stroke-width": 2.2 }, g);
  // Data.
  const dots = el("g", {}, g);
  ser.t.forEach((t, i) => {
    if (t > tmax) return;
    const open = !ser.m[i];
    el("circle", { cx: X(t), cy: Y(ser.L[i]), r: 3.2, fill: open ? css("--surface") : c, stroke: open ? c : css("--surface"),
                   "stroke-width": open ? 1.4 : 0.8 }, dots);
  });
  // Nearest-point hover over the whole plot.
  const hit = el("rect", { x: ml, y: mt, width: W - ml - mr, height: H - mt - mb, fill: "transparent" }, svg);
  const tip = document.getElementById("tip");
  hit.addEventListener("mousemove", ev => {
    const r = svg.getBoundingClientRect(), sx = (ev.clientX - r.left) * W / r.width, sy = (ev.clientY - r.top) * H / r.height;
    let best = -1, bd = 1e9;
    ser.t.forEach((t, i) => { if (t > tmax) return; const d = (X(t) - sx) ** 2 + (Y(ser.L[i]) - sy) ** 2; if (d < bd) { bd = d; best = i; } });
    if (best < 0 || bd > 900) { tip.style.display = "none"; return; }
    const t = ser.t[best], mod = sn.type === "Ia" || t >= tStart ? interp(model, t) : null;
    tip.innerHTML = `day ${t.toFixed(1)}<br>observed ${fmtL(ser.L[best])} erg/s` +
      (mod ? `<br>model ${fmtL(mod)} erg/s` : "") + (ser.m[best] ? "" : "<br><i>near-infrared estimated</i>");
    tip.style.display = "block"; tip.style.left = (ev.clientX + 14) + "px"; tip.style.top = (ev.clientY + 10) + "px";
  });
  hit.addEventListener("mouseleave", () => { tip.style.display = "none"; });

  // Readout.
  let rd = `<b>${fmtM(s.m)}</b> solar masses of nickel-56 heat the debris with <b>${fmtL(s.m * deposit(1, 0, true))}</b> erg/s on day 1.`;
  if (sn.type === "Ia") {
    let pk = model.reduce((a, b) => b[1] > a[1] ? b : a, [0, 0]);
    rd += `<br><br>Model peak: <b>${fmtL(pk[1])}</b> erg/s on day <b>${pk[0].toFixed(1)}</b>. ` +
          `Arnett's rule (peak light = heating at that moment): ${fmtL(s.m * deposit(pk[0], s.t0, s.pos))} erg/s.`;
  }
  const t1 = Math.min(tmax, 300);
  const kept = s.t0 > 0 ? 100 * (1 - Math.exp(-((s.t0 / t1) ** 2))) : 100;
  rd += `<br><br>At day ${t1}, ${kept.toFixed(0)}% of the gamma rays still heat the debris.`;
  rd += `<br><br>The project's fit: ${fmtM(sn.fit.m)} Msun, ` + (sn.type === "Ia" ? `diffusion ${sn.fit.tau} d, ` : "") +
        `T0 = ${sn.fit.t0} d. ZTF Type Ia median: ${D.ztf_median} Msun.`;
  document.getElementById("readout").innerHTML = rd;
  document.getElementById("note").textContent = sn.note;
  document.getElementById("legend").innerHTML =
    `<span><svg width="26" height="10"><line x1="0" y1="5" x2="26" y2="5" stroke="${c}" stroke-width="2.2"/></svg>model</span>` +
    `<span><svg width="26" height="10"><line x1="0" y1="5" x2="26" y2="5" stroke="${css("--model")}" stroke-width="1.2" stroke-dasharray="5 4"/></svg>heating, every gamma ray trapped</span>` +
    `<span><svg width="10" height="10"><circle cx="5" cy="5" r="4" fill="${c}"/></svg>${sn.label}</span>`;
}

function interp(arr, t) {
  for (let i = 1; i < arr.length; i++) if (arr[i][0] >= t) {
    const [t0, y0] = arr[i - 1], [t1, y1] = arr[i];
    return y0 + (y1 - y0) * (t - t0) / (t1 - t0);
  }
  return null;
}

function syncControls() {
  const s = st[cur], sn = D.supernovae[cur];
  document.getElementById("m").value = Math.log10(s.m);
  document.getElementById("tau").value = s.tau || 15;
  document.getElementById("t0").value = s.t0 > 0 ? Math.log10(s.t0) : 3.3;
  document.getElementById("pos").checked = s.pos;
  document.getElementById("taubox").classList.toggle("disabled", sn.type !== "Ia");
  document.getElementById("tau").disabled = sn.type !== "Ia";
  labels();
}
function labels() {
  const s = st[cur];
  document.getElementById("mv").textContent = fmtM(s.m) + " Msun";
  document.getElementById("tauv").textContent = D.supernovae[cur].type === "Ia" ? s.tau.toFixed(1) + " d" : "n/a";
  document.getElementById("t0v").textContent = s.t0 > 0 ? Math.round(s.t0) + " d" : "never";
}
function table() {
  const ser = D.supernovae[cur].series;
  let h = "<tr><th>day</th><th>luminosity (erg/s)</th><th>error</th><th>near-IR</th></tr>";
  ser.t.forEach((t, i) => { h += `<tr><td>${t.toFixed(1)}</td><td>${fmtL(ser.L[i])}</td><td>${fmtL(ser.e[i])}</td><td>${ser.m[i] ? "measured" : "estimated"}</td></tr>`; });
  document.getElementById("tbl").innerHTML = h;
}
function select(sn) {
  cur = sn;
  document.querySelectorAll("#sns button").forEach(b => b.setAttribute("aria-pressed", b.dataset.sn === sn));
  syncControls(); table(); draw();
}

for (const sn in D.supernovae) {
  const f = D.supernovae[sn].fit;
  st[sn] = { m: f.m, tau: f.tau || 15, t0: f.t0, pos: true };
  const b = document.createElement("button");
  b.textContent = D.supernovae[sn].label + " (Type " + D.supernovae[sn].type + ")";
  b.dataset.sn = sn;
  b.addEventListener("click", () => select(sn));
  document.getElementById("sns").appendChild(b);
}
document.getElementById("m").addEventListener("input", e => { st[cur].m = 10 ** +e.target.value; labels(); draw(); });
document.getElementById("tau").addEventListener("input", e => { st[cur].tau = +e.target.value; labels(); draw(); });
document.getElementById("t0").addEventListener("input", e => { st[cur].t0 = 10 ** +e.target.value; labels(); draw(); });
document.getElementById("pos").addEventListener("change", e => { st[cur].pos = e.target.checked; draw(); });
document.getElementById("fit").addEventListener("click", () => {
  const f = D.supernovae[cur].fit; st[cur] = { m: f.m, tau: f.tau || 15, t0: f.t0, pos: true }; syncControls(); draw();
});
document.getElementById("trap").addEventListener("click", () => { st[cur].t0 = 0; syncControls(); draw(); });
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);
let resizeTimer;
window.addEventListener("resize", () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(draw, 120); });
select("SN2011fe");
</script>
</body>
</html>
"""


def write(r: dict):
    path = C.FIGURE_DIR / "nickel_explorer.html"
    path.write_text(PAGE.replace("__DATA__", json.dumps(_payload(r), separators=(",", ":"))))
    return path
