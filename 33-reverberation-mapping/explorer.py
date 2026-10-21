"""An offline page for finding an echo yourself: figures/echo_finder.html.

Pick a light curve pair and drag the slider: the emission-line light curve
slides back in time over the continuum, and the cross-correlation at that
shift is read off the ICCF this project computed. Where the two line up best
is the lag; the page turns it into a size (light-days, AU) and, with the
published line width, a black hole mass.

The page is self-contained (no network, no libraries): data embedded as JSON,
drawn as SVG, light and dark themes, a table view of the data.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import iccf
import lightcurves as LC
import rl as RL

#: SDSS-RM quasars shown: RM078's H-beta (Shen et al. 2024's own showcase), then
#: for Mg II and C IV the lag with the highest correlation among those graded
#: 4-5 by Shen et al. and detected by both their pipeline and ours.


def sdssrm_picks(ours) -> list[tuple[int, str]]:
    t = RL.sdssrm_table().set_index("RMID")
    picks = [(78, "HB")]
    for line in ("MG2", "C4"):
        o = ours[(ours.line == line) & ours.ok & ours.detected].copy()
        o = o[[t.loc[r, f"{line}_grade"] >= 4 and t.loc[r, f"{line}_lag"] > 0 for r in o.RMID]]
        if len(o):
            picks.append((int(o.sort_values("rmax_near", ascending=False).iloc[0].RMID), line))
    return picks


def _std(y):
    y = np.asarray(y, float)
    return (y - y.mean()) / y.std()


def _r4(a):
    return [float(f"{v:.4g}") for v in np.asarray(a, float)]


def _entry(key, label, note, z, cont, line, lags, r, lag_obs, lag_err, sigma, sigma_ref, xlabel, t0=0.0,
           cont_label="continuum", line_label="line", detrend=False):
    yc = iccf.detrend(cont.t, cont.y) if detrend else cont.y
    yl = iccf.detrend(line.t, line.y) if detrend else line.y
    return dict(key=key, label=label, note=note, z=z, xlabel=xlabel,
                tc=_r4(cont.t - t0), yc=_r4(_std(yc)), tl=_r4(line.t - t0), yl=_r4(_std(yl)),
                lags=_r4(lags), r=[round(float(v), 3) for v in r], lag=round(float(lag_obs), 2),
                lag_err=[round(float(e), 2) for e in lag_err], sigma=sigma, sigma_ref=sigma_ref,
                cont_label=cont_label, line_label=line_label)


def _payload(res: dict) -> dict:
    out = []
    z1 = 1 + C.Z_NGC5548
    # NGC 5548 in 2014
    lc = res["storm"]["lcs"]
    raw = res["storm"]["raw"][("Hb", "F5100", "full")]
    m, lo, hi = raw.cent
    out.append(_entry("n5548_2014", "NGC 5548, 2014", "AGN STORM: 133 nights of spectra from the MDM Observatory. The "
                      "second half of the campaign is the 'BLR holiday': the line partly stopped following the "
                      "continuum, so a perfect match is impossible there.", C.Z_NGC5548, lc["F5100"], lc["Hb"],
                      raw.lags, raw.r, m, (lo, hi), 4278, "rms spectrum; Pei et al. 2017", "HJD − 2,450,000",
                      cont_label="5100 Å continuum", line_label="broad Hβ", detrend=True))
    # Selected AGN Watch seasons
    for s in res["watch"]["seasons"]:
        yr = int(round(s["year"] - 0.4))
        if yr not in (1989, 1998, 2005):
            continue
        w = res["watch"]["lcs"]
        c = w["F5100"].window(s["jd_start"], s["jd_end"])
        h = w["Hb"].window(s["jd_start"], s["jd_end"])
        rr = s["raw"]
        m, lo, hi = rr.cent
        notes = {1989: "AGN Watch year 1 (1988-89): the classic campaign that first measured a broad-line region's "
                       "size this way.",
                 1998: "1998 (AGN Watch year 11): NGC 5548 near its brightest; the longest lag it ever showed.",
                 2005: "2005: the faintest state ever monitored (Bentz et al. 2007): a short campaign and a small BLR."}
        out.append(_entry(f"n5548_{yr}", f"NGC 5548, {yr}", notes[yr], C.Z_NGC5548, c, h, rr.lags, rr.r, m, (lo, hi),
                          s["pub_sigma"][0], "rms spectrum; AGN Black Hole Mass Database", "JD − 2,400,000",
                          cont_label="5100 Å AGN continuum", line_label="broad Hβ"))
    # SDSS-RM quasars
    ours = res["sdssrm_ours"]
    t = RL.sdssrm_table().set_index("RMID")
    names = {"HB": ("hb", "broad Hβ"), "MG2": ("mg2", "broad Mg II"), "C4": ("c4", "broad C IV")}
    for rmid, line in sdssrm_picks(ours):
        row = ours[(ours.RMID == rmid) & (ours.line == line)]
        if row.empty or not bool(row.iloc[0].ok):
            continue
        row = row.iloc[0]
        cont = LC.sdssrm(rmid, "cont")
        ln = LC.sdssrm(rmid, names[line][0])
        lags = iccf.lag_grid(np.floor(row.search_lo / 2) * 2, np.ceil(row.search_hi / 2) * 2, C.SDSSRM_LAG_STEP)
        r, _ = iccf.ccf(cont.t, cont.y, ln.t, ln.y, lags)
        z = float(t.loc[rmid, "z"])
        sig = float(t.loc[rmid, f"{line}_sig_rms"])
        extra = ""
        if rmid == 78:
            extra = (f" For this quasar Shen et al.'s fiducial lag (PyROA) is {t.loc[rmid, 'HB_lag']:.0f} days observed: "
                     "its correlation peak is broad, so the centroid (used here) and the peak disagree. A lag is a "
                     "measurement with a method attached.")
        out.append(_entry(f"rm{rmid:03d}", f"SDSS-RM {rmid:03d} ({names[line][1].removeprefix('broad ')})",
                          f"A quasar at z = {z:.2f}: 11 years of imaging from four telescopes against 7 years of "
                          "SDSS spectra. Its lag of months is stretched by (1+z) in our clocks; seasonal gaps "
                          "make some shifts look good by accident (aliases)." + extra, z, cont, ln, lags, r, row.lag,
                          (row.lag_lo, row.lag_hi), round(sig) if np.isfinite(sig) else None,
                          "rms spectrum; Shen et al. 2024", "MJD", cont_label="continuum (g, r, i)",
                          line_label=names[line][1]))
    del z1
    return {"sets": out, "logf": res["f_dyn"]["mean"][0],
            "colors": {"cont": ["#2a78d6", "#3987e5"], "line": ["#eb6834", "#d95926"]},
            "consts": {"c": C.C_CGS, "G": C.G_CGS, "Msun": C.MSUN_G, "day": C.DAY_S, "au": C.AU_CM}}


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Echo Finder</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --model: #52514e; --mark: #0b0b0b;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
    --grid: #2f2f2d; --axis: #4a4a46; --model: #c3c2b7; --mark: #ffffff;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
  --grid: #2f2f2d; --axis: #4a4a46; --model: #c3c2b7; --mark: #ffffff;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1060px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 16px; max-width: 780px; }
.row { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 12px; }
button { font: inherit; font-size: 0.9rem; padding: 6px 12px; border-radius: 8px; cursor: pointer;
         border: 1px solid var(--axis); background: var(--surface); color: var(--ink); }
button[aria-pressed="true"] { background: var(--ink); color: var(--surface); border-color: var(--ink); }
button:focus-visible, input:focus-visible { outline: 2px solid #2a78d6; outline-offset: 2px; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 290px; gap: 16px; }
@media (max-width: 820px) { .layout { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border-radius: 12px; padding: 12px; min-width: 0; }
svg { width: 100%; height: auto; display: block; }
.ctl { margin-bottom: 14px; }
.ctl label { display: flex; justify-content: space-between; font-size: 0.88rem; color: var(--ink2); }
.ctl label b { color: var(--ink); font-variant-numeric: tabular-nums; }
.ctl input[type=range] { width: 100%; }
.readout { font-size: 0.88rem; color: var(--ink2); border-top: 1px solid var(--grid); padding-top: 10px; }
.readout div { margin-bottom: 6px; }
.readout b { color: var(--ink); font-variant-numeric: tabular-nums; }
.big { font-size: 1.6rem; color: var(--ink); line-height: 1.2; }
.note { font-size: 0.88rem; color: var(--ink2); margin: 10px 2px 0; }
.legend { display: flex; flex-wrap: wrap; gap: 14px; font-size: 0.82rem; color: var(--ink2); margin: 6px 2px 0; }
.legend span { display: inline-flex; align-items: center; gap: 6px; white-space: nowrap; }
details { margin-top: 16px; font-size: 0.85rem; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
td, th { padding: 2px 10px; text-align: right; border-bottom: 1px solid var(--grid); }
.tablewrap { max-height: 280px; overflow: auto; }
.axis text { fill: var(--muted); font-size: 11px; }
.axis line { stroke: var(--axis); }
.gridline { stroke: var(--grid); }
.lbl { fill: var(--ink2); font-size: 12px; }
footer { color: var(--muted); font-size: 0.78rem; margin-top: 18px; }
</style>
</head>
<body>
<main>
  <h1>Find the echo</h1>
  <p class="lede">The disk of hot gas around a black hole flickers. Clouds of gas farther out light up in the broad
  emission lines when the flicker reaches them, a little later. Slide the line's light curve back in time until it
  matches the continuum: the shift is the light-travel time to the clouds, and with the speed of the clouds
  (the line's width) it weighs the black hole.</p>
  <div class="row" id="sets" role="group" aria-label="light curves"></div>
  <div class="layout">
    <div class="card">
      <svg id="lc" viewBox="0 0 720 300" role="img" aria-label="continuum and shifted line light curves"></svg>
      <div class="legend" id="legend"></div>
      <svg id="cc" viewBox="0 0 720 200" role="img" aria-label="cross-correlation against shift"></svg>
      <p class="note" id="note"></p>
    </div>
    <div class="card">
      <div class="ctl"><label for="shift">Move the line back by <b id="sv"></b></label>
        <input id="shift" type="range" step="0.25"></div>
      <div class="row">
        <button id="best">jump to the measured lag</button>
        <button id="zero">no shift</button>
      </div>
      <div class="readout" id="readout"></div>
    </div>
  </div>
  <details><summary>Table view of the data</summary><div class="tablewrap"><table id="tbl"></table></div></details>
  <footer>Data: Pei et al. 2017 (NGC 5548, 2014); Peterson et al. 2013 (NGC 5548, 1988-2008); Shen et al. 2024
  (SDSS-RM). Both curves are scaled to zero mean and unit variance so their shapes can be compared. The correlation
  curve is this project's interpolated cross-correlation function; the measured lag is the median of 500-1000
  resamplings of the data. Mass = f c τ σ² / G with log f = __LOGF__ (dynamical models; uncertain by a factor of 2
  for any single AGN).</footer>
</main>
<script>
const D = __DATA__;
const NS = "http://www.w3.org/2000/svg";
let cur = 0;
const dark = () => {
  const forced = document.documentElement.dataset.theme;
  if (forced) return forced === "dark";
  return window.matchMedia("(prefers-color-scheme: dark)").matches;
};
const col = k => D.colors[k][dark() ? 1 : 0];
function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function ticks(lo, hi, n) {
  const span = hi - lo, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}
function interp(x, xs, ys) {
  if (x <= xs[0]) return ys[0]; if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
  let lo = 0, hi = xs.length - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] > x) hi = m; else lo = m; }
  return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
}
function fmt(v, d) { return v.toLocaleString("en-US", { maximumFractionDigits: d, minimumFractionDigits: d }); }
function buttons() {
  const box = document.getElementById("sets");
  D.sets.forEach((s, i) => {
    const b = document.createElement("button");
    b.textContent = s.label; b.setAttribute("aria-pressed", i === cur);
    b.onclick = () => { cur = i; [...box.children].forEach((c, j) => c.setAttribute("aria-pressed", j === cur)); setup(); };
    box.appendChild(b);
  });
}
const slider = document.getElementById("shift");
function setup() {
  const s = D.sets[cur];
  slider.min = Math.max(s.lags[0], -0.5 * (s.tl[s.tl.length - 1] - s.tl[0]));
  slider.max = Math.min(s.lags[s.lags.length - 1], 0.6 * (s.tl[s.tl.length - 1] - s.tl[0]) + Math.abs(s.lag));
  slider.step = (s.lags[1] - s.lags[0]) / 2;
  slider.value = 0;
  document.getElementById("note").textContent = s.note;
  document.getElementById("legend").innerHTML =
    (s.key.startsWith("rm")
      ? `<span><svg width="12" height="12"><circle cx="6" cy="6" r="2.5" fill="${col("cont")}"/></svg>${s.cont_label}</span>`
      : `<span><svg width="22" height="10"><line x1="0" y1="5" x2="22" y2="5" stroke="${col("cont")}" stroke-width="2"/></svg>${s.cont_label}</span>`) +
    `<span><svg width="12" height="12"><circle cx="6" cy="6" r="4" fill="${col("line")}"/></svg>${s.line_label}, moved back by the slider</span>`;
  table(); draw();
}
function draw() {
  const s = D.sets[cur], sh = +slider.value;
  // light curves
  const svg = document.getElementById("lc"); svg.innerHTML = "";
  const W = 720, H = 300, L = 46, R = 12, T = 12, B = 36;
  const xmin = Math.min(s.tc[0], s.tl[0] - Math.max(sh, 0)), xmax = Math.max(s.tc[s.tc.length - 1], s.tl[s.tl.length - 1] - Math.min(sh, 0));
  const ymin = -3.2, ymax = 3.2;
  const X = x => L + (x - xmin) / (xmax - xmin) * (W - L - R), Y = y => T + (ymax - y) / (ymax - ymin) * (H - T - B);
  const ax = el("g", { class: "axis" }, svg);
  for (const v of ticks(xmin, xmax, 7)) {
    el("line", { x1: X(v), x2: X(v), y1: T, y2: H - B, class: "gridline" }, ax);
    el("text", { x: X(v), y: H - B + 16, "text-anchor": "middle" }, ax).textContent = v;
  }
  for (const v of [-2, 0, 2]) {
    el("line", { x1: L, x2: W - R, y1: Y(v), y2: Y(v), class: "gridline" }, ax);
    el("text", { x: L - 6, y: Y(v) + 4, "text-anchor": "end" }, ax).textContent = v;
  }
  el("text", { x: (L + W - R) / 2, y: H - 4, "text-anchor": "middle", class: "lbl" }, svg).textContent = s.xlabel + " (days)";
  if (s.key.startsWith("rm")) {
    const gc = el("g", {}, svg);
    s.tc.forEach((t, i) => el("circle", { cx: X(t), cy: Y(s.yc[i]), r: 1.6, fill: col("cont") }, gc));
  } else {
    let d = ""; s.tc.forEach((t, i) => {
      const gap = i === 0 || t - s.tc[i - 1] > 30;
      d += (gap ? "M" : "L") + X(t).toFixed(1) + " " + Y(s.yc[i]).toFixed(1);
    });
    el("path", { d: d, fill: "none", stroke: col("cont"), "stroke-width": 1.6 }, svg);
  }
  const g = el("g", {}, svg);
  s.tl.forEach((t, i) => el("circle", { cx: X(t - sh), cy: Y(s.yl[i]), r: 2.6, fill: col("line"), stroke: getComputedStyle(document.body).backgroundColor, "stroke-width": 0.6 }, g));
  // correlation curve
  const cc = document.getElementById("cc"); cc.innerHTML = "";
  const H2 = 200, B2 = 34;
  const lx0 = +slider.min, lx1 = +slider.max;
  const X2 = x => L + (x - lx0) / (lx1 - lx0) * (W - L - R), Y2 = y => T + (1 - y) / 2 * (H2 - T - B2);
  const ax2 = el("g", { class: "axis" }, cc);
  for (const v of ticks(lx0, lx1, 8)) {
    el("line", { x1: X2(v), x2: X2(v), y1: T, y2: H2 - B2, class: "gridline" }, ax2);
    el("text", { x: X2(v), y: H2 - B2 + 16, "text-anchor": "middle" }, ax2).textContent = v;
  }
  for (const v of [-1, 0, 1]) {
    el("line", { x1: L, x2: W - R, y1: Y2(v), y2: Y2(v), class: "gridline" }, ax2);
    el("text", { x: L - 6, y: Y2(v) + 4, "text-anchor": "end" }, ax2).textContent = v;
  }
  el("text", { x: (L + W - R) / 2, y: H2 - 4, "text-anchor": "middle", class: "lbl" }, cc).textContent = "shift (observed days)";
  el("text", { x: L + 4, y: T + 12, class: "lbl" }, cc).textContent = "how well they match (correlation r)";
  let d2 = ""; s.lags.forEach((t, i) => { if (t >= lx0 && t <= lx1) d2 += (d2 ? "L" : "M") + X2(t).toFixed(1) + " " + Y2(s.r[i]).toFixed(1); });
  el("path", { d: d2, fill: "none", stroke: getComputedStyle(document.documentElement).getPropertyValue("--model"), "stroke-width": 1.8 }, cc);
  el("rect", { x: X2(s.lag - s.lag_err[0]), width: Math.max(1, X2(s.lag + s.lag_err[1]) - X2(s.lag - s.lag_err[0])), y: T, height: H2 - T - B2, fill: col("line"), opacity: 0.18 }, cc);
  const rnow = interp(sh, s.lags, s.r);
  el("line", { x1: X2(sh), x2: X2(sh), y1: T, y2: H2 - B2, stroke: getComputedStyle(document.documentElement).getPropertyValue("--mark"), "stroke-width": 1.2 }, cc);
  el("circle", { cx: X2(sh), cy: Y2(rnow), r: 4.5, fill: col("line"), stroke: getComputedStyle(document.documentElement).getPropertyValue("--mark"), "stroke-width": 1 }, cc);
  // readout
  const z1 = 1 + s.z, tauRest = sh / z1;
  document.getElementById("sv").textContent = fmt(sh, 2) + " d";
  const k = D.consts;
  let html = `<div>correlation at this shift<br><span class="big">r = ${fmt(rnow, 2)}</span></div>`;
  html += `<div>measured lag (median of resamplings): <b>${fmt(s.lag / z1, 1)} days</b> in the rest frame` +
          ` (${fmt(s.lag, 1)} observed; −${fmt(s.lag_err[0], 1)}/+${fmt(s.lag_err[1], 1)})</div>`;
  if (tauRest > 0) {
    const au = k.c * tauRest * k.day / k.au;
    html += `<div>your shift puts the clouds <b>${fmt(tauRest, 1)} light-days</b> out (${fmt(au, 0)} AU, ${fmt(au / 30.07, 0)}× Neptune's orbit)</div>`;
    if (s.sigma) {
      const vp = k.c * tauRest * k.day * Math.pow(s.sigma * 1e5, 2) / k.G / k.Msun;
      const m = vp * Math.pow(10, D.logf);
      html += `<div>with line width σ = <b>${fmt(s.sigma, 0)} km/s</b> (${s.sigma_ref}), black hole mass ` +
              `<b>${(m / 1e7).toFixed(1)} × 10⁷ suns</b></div>`;
    }
  } else {
    html += `<div>a negative shift would mean the line changed <i>before</i> the continuum: not an echo.</div>`;
  }
  document.getElementById("readout").innerHTML = html;
}
function table() {
  const s = D.sets[cur];
  let h = `<tr><th>time</th><th>${s.cont_label}</th><th>time</th><th>${s.line_label}</th></tr>`;
  const n = Math.max(s.tc.length, s.tl.length);
  for (let i = 0; i < n; i++) {
    h += `<tr><td>${s.tc[i] ?? ""}</td><td>${s.yc[i] ?? ""}</td><td>${s.tl[i] ?? ""}</td><td>${s.yl[i] ?? ""}</td></tr>`;
  }
  document.getElementById("tbl").innerHTML = h;
}
slider.addEventListener("input", draw);
document.getElementById("best").onclick = () => { slider.value = D.sets[cur].lag; draw(); };
document.getElementById("zero").onclick = () => { slider.value = 0; draw(); };
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", setup);
buttons(); setup();
</script>
</body>
</html>
"""


def write(res: dict):
    data = _payload(res)
    html = PAGE.replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__LOGF__", f"{data['logf']:.2f}")
    out = C.FIGURE_DIR / "echo_finder.html"
    out.write_text(html)
    return out
