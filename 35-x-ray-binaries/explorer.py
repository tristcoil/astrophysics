"""An offline page to weigh black holes yourself: figures/black_hole_explorer.html.

Two tools on one self-contained page (data embedded as JSON, drawn as SVG,
no network, light and dark themes):

1. Weigh a black hole: pick a system, drag the inclination, and watch its
   mass function turn into a mass. The inclination range measured from the
   companion's light curve is shaded; neutron stars cannot exceed the band
   at 2.1-3.2 solar masses.
2. Replay an outburst: step through MAXI J1820+070's 2018 outburst day by day
   on its light curve and its hardness-intensity diagram.
"""

from __future__ import annotations

import json

import numpy as np

from config import FIGURE_DIR


def _human(spec, unit: str = "") -> str:
    """('N', 51, 0.9) -> '51 ± 0.9°'; ('I', 37, 78) -> '37–78°'; None -> 'unknown'."""
    if spec is None:
        return "unknown"
    kind = spec[0]
    if kind == "N":
        return f"{spec[1]:g} ± {spec[2]:g}{unit}"
    if kind == "SN":
        return f"{spec[1]:g} +{spec[2]:g}/−{spec[3]:g}{unit}"
    if kind in ("U", "I"):
        return f"{spec[1]:g}–{spec[2]:g}{unit}"
    if kind == "I_ecl":
        return f"above {spec[1]:g}{unit} (no eclipses)"
    if kind == "max":
        return f"below {spec[1]:g}{unit}"
    return str(spec[1])


def _system_data(results, art) -> list[dict]:
    out = []
    draws = art["weighing"]["draws"]
    adopted = {a["name"]: a for a in art["weighing"]["adopted"]}
    for r in results["weighing"]["systems"]:
        d = draws[r["name"]]
        a = adopted[r["name"]]
        rec = {
            "name": r["name"], "common": r["common_name"], "sptype": r["sptype"],
            "porb_h": r["porb_h"][0], "k2": r["k2_kms"][0], "k2_err": r["k2_kms"][1],
            "fm": round(float(np.median(d["fm"])), 4), "q": round(float(np.median(d["q"])), 4),
            "q_txt": _human(a["q_spec"]), "i_txt": _human(a["i_spec"], "°"), "m1_pub": r["m1_published"],
            "refs": "; ".join(sorted(set(filter(None, [r["k2_ref"], r["q_ref"], r["i_ref"]])))),
        }
        if "incl" in d:
            lo, hi = np.quantile(d["incl"], [0.16, 0.84])
            rec["i_lo"], rec["i_hi"] = round(float(lo), 2), round(float(hi), 2)
            m = d["m1"]
            rec["m_med"] = round(float(np.median(m)), 2)
            rec["m_lo"], rec["m_hi"] = (round(float(x), 2) for x in np.quantile(m, [0.16, 0.84]))
        out.append(rec)
    return sorted(out, key=lambda x: x["fm"], reverse=True)


def _outburst_data(results, art) -> dict:
    df = art["outbursts"]["j1820"]
    tr = results["outbursts"]["j1820"]["transitions"]
    ok = np.isfinite(df["hr"]) & np.isfinite(df["crab_2-20"]) & (df["crab_2-20"] > 0)
    d = df[ok]
    return {
        "mjd": [float(x) for x in d["mjd"]],
        "flux": [round(float(x), 4) for x in d["crab_2-20"]],
        "hr": [round(float(x), 4) for x in d["hr"]],
        "bat": [None if not np.isfinite(x) else round(float(x), 4) for x in d["crab_bat"]],
        "h2s": tr["hard_to_soft_mjd"], "s2h": tr["soft_to_hard_mjd"],
    }


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Black Hole Explorer</title>
<style>
:root { color-scheme: light; --bg:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --grid:#e1e0d9;
  --axis:#c3c2b7; --gap:#efeee9; --ours:#eb6834; --pub:#4a3aa7; --curve:#2a78d6; --band:#cde2fb; --card:#ffffff; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme: dark; --bg:#1a1a19;
  --ink:#ffffff; --ink2:#c3c2b7; --muted:#9a998f; --grid:#33332f; --axis:#55544e; --gap:#262624; --ours:#d95926;
  --pub:#9085e9; --curve:#3987e5; --band:#1c3a5e; --card:#222220; } }
:root[data-theme="dark"] { color-scheme: dark; --bg:#1a1a19; --ink:#ffffff; --ink2:#c3c2b7; --muted:#9a998f;
  --grid:#33332f; --axis:#55544e; --gap:#262624; --ours:#d95926; --pub:#9085e9; --curve:#3987e5; --band:#1c3a5e;
  --card:#222220; }
body { margin:0; background:var(--bg); color:var(--ink); font:15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width:1060px; margin:0 auto; padding:24px 16px 48px; }
h1 { font-size:26px; margin:0 0 4px; } h2 { font-size:18px; margin:32px 0 6px; }
p.lede { color:var(--ink2); margin:0 0 16px; max-width:760px; }
.card { background:var(--card); border:1px solid var(--grid); border-radius:10px; padding:16px; }
.row { display:flex; gap:16px; flex-wrap:wrap; align-items:center; }
.controls label { color:var(--ink2); font-size:14px; }
select, button { font:inherit; color:var(--ink); background:var(--bg); border:1px solid var(--axis); border-radius:6px; padding:4px 8px; }
input[type=range] { width:260px; accent-color:var(--ours); }
.readout { display:grid; grid-template-columns:repeat(auto-fit, minmax(150px, 1fr)); gap:10px; margin:12px 0; }
.stat { border-left:3px solid var(--grid); padding:2px 10px; }
.stat .k { font-size:12.5px; color:var(--muted); } .stat .v { font-size:20px; font-weight:600; }
svg { width:100%; height:auto; display:block; }
svg text { fill:var(--ink2); font-size:12px; } .axis line, .axis path { stroke:var(--axis); }
.gridl { stroke:var(--grid); stroke-width:1; } .note { color:var(--ink2); font-size:13.5px; }
.two { display:grid; grid-template-columns:1fr 1fr; gap:16px; }
@media (max-width:760px) { .two { grid-template-columns:1fr; } input[type=range] { width:100%; } }
footer { color:var(--muted); font-size:12.5px; margin-top:32px; }
</style></head>
<body><main>
<h1>Black Hole Explorer</h1>
<p class="lede">The companion star's orbit gives the <b>mass function</b> f(M) = P K³ / 2πG, a strict lower limit on the
unseen mass. Turning it into a mass needs the orbit's tilt. Drag the inclination and watch the mass change. The real
measurements of 20 black-hole transients come from BlackCAT.</p>

<div class="card">
  <div class="row controls">
    <label>Black hole <select id="sys"></select></label>
    <label>Inclination <input id="inc" type="range" min="10" max="90" step="0.5" value="60"> <span id="incv"></span></label>
  </div>
  <div class="readout">
    <div class="stat"><div class="k">orbital period</div><div class="v" id="r-p"></div></div>
    <div class="stat"><div class="k">companion's speed K₂</div><div class="v" id="r-k"></div></div>
    <div class="stat"><div class="k">mass function f(M)</div><div class="v" id="r-f"></div></div>
    <div class="stat"><div class="k">mass at this inclination</div><div class="v" id="r-m"></div></div>
  </div>
  <svg id="plot" viewBox="0 0 900 380" role="img" aria-label="Black hole mass against inclination"></svg>
  <p class="note" id="verdict"></p>
</div>

<h2>Replay an outburst: MAXI J1820+070, 2018</h2>
<p class="lede">Daily X-ray fluxes from MAXI, in units of the Crab Nebula. The hardness is the ratio of the 4–10 keV
to the 2–4 keV flux. A Crab-like spectrum has hardness 1.</p>
<div class="card">
  <div class="row controls">
    <button id="play">▶ Play</button>
    <input id="day" type="range" min="0" value="0" step="1" style="width:420px">
    <span id="dayv" class="note"></span>
  </div>
  <div class="two">
    <svg id="lc" viewBox="0 0 480 300" role="img" aria-label="Light curve"></svg>
    <svg id="hid" viewBox="0 0 480 300" role="img" aria-label="Hardness-intensity diagram"></svg>
  </div>
  <p class="note" id="state"></p>
</div>
<footer>Data: BlackCAT online (Corral-Santana et al. 2016 and updates), MAXI/GSC (RIKEN, JAXA, MAXI team).
Masses here use each system's median mass ratio q; the project's figures propagate every uncertainty.
Made by project I4 of the Astrophysics Project Atlas (AI-generated).</footer>
</main>
<script>
const SYS = __SYSTEMS__;
const OB = __OUTBURST__;
const NS = "http://www.w3.org/2000/svg";
const $ = id => document.getElementById(id);
function el(tag, attrs, parent, text) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (text !== undefined) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const mass = (f, q, i) => f * (1 + q) ** 2 / Math.sin(i * Math.PI / 180) ** 3;

// ---------------------------------------------------------------- weighing
const sel = $("sys");
SYS.forEach((s, k) => { const o = document.createElement("option"); o.value = k;
  o.textContent = s.common ? `${s.name} (${s.common.split(" (")[0]})` : s.name; sel.appendChild(o); });
function drawWeigh() {
  const s = SYS[+sel.value], inc = +$("inc").value, svg = $("plot");
  svg.innerHTML = "";
  const W = 900, H = 380, L = 64, R = 20, T = 16, B = 44;
  const x = i => L + (i - 10) / 80 * (W - L - R);
  const y = m => T + (Math.log10(100) - Math.log10(Math.min(Math.max(m, 1), 100))) / 2 * (H - T - B);
  el("rect", {x: L, y: y(3.2), width: W - L - R, height: y(2.08) - y(3.2), fill: css("--gap")}, svg);
  el("text", {x: L + 6, y: y(2.6) + 4}, svg, "neutron stars stop here (2.1–3.2 M☉)");
  if (s.i_lo !== undefined) {
    el("rect", {x: x(s.i_lo), y: T, width: x(s.i_hi) - x(s.i_lo), height: H - T - B, fill: css("--band"), opacity: 0.6}, svg);
    el("text", {x: (x(s.i_lo) + x(s.i_hi)) / 2, y: T + 14, "text-anchor": "middle"}, svg, "measured inclination");
  }
  for (const m of [1, 2, 3, 5, 10, 20, 50, 100]) {
    el("line", {x1: L, x2: W - R, y1: y(m), y2: y(m), class: "gridl"}, svg);
    el("text", {x: L - 8, y: y(m) + 4, "text-anchor": "end"}, svg, m);
  }
  for (let i = 10; i <= 90; i += 10) {
    el("line", {x1: x(i), x2: x(i), y1: H - B, y2: H - B + 5, stroke: css("--axis")}, svg);
    el("text", {x: x(i), y: H - B + 18, "text-anchor": "middle"}, svg, i + "°");
  }
  el("text", {x: (L + W - R) / 2, y: H - 6, "text-anchor": "middle"}, svg, "inclination (90° = edge-on)");
  el("text", {x: 14, y: (T + H - B) / 2, transform: `rotate(-90 14 ${(T + H - B) / 2})`, "text-anchor": "middle"}, svg, "mass [M☉]");
  let d = "";
  for (let i = 10; i <= 90; i += 0.25) {
    const mm = mass(s.fm, s.q, i);
    if (mm > 100) continue;  // above the plot: start the curve where it enters
    d += (d ? "L" : "M") + x(i).toFixed(1) + "," + y(mm).toFixed(1);
  }
  el("path", {d, fill: "none", stroke: css("--curve"), "stroke-width": 2}, svg);
  const m = mass(s.fm, s.q, inc);
  el("circle", {cx: x(inc), cy: y(m), r: 7, fill: css("--ours"), stroke: css("--card"), "stroke-width": 2}, svg);
  $("incv").textContent = inc.toFixed(1) + "°";
  $("r-p").textContent = s.porb_h < 48 ? s.porb_h.toFixed(2) + " h" : (s.porb_h / 24).toFixed(2) + " d";
  $("r-k").textContent = s.k2.toFixed(0) + " km/s";
  $("r-f").textContent = s.fm.toFixed(2) + " M☉";
  $("r-m").textContent = m.toFixed(1) + " M☉";
  const lim = s.fm > 3.2 ? "Its mass function alone exceeds any possible neutron star: a black hole at any inclination." :
    `Its mass function (${s.fm.toFixed(2)} M☉) alone does not rule out a neutron star; the inclination and mass ratio do.`;
  const meas = s.m_med !== undefined ? ` With the measured inclination (${s.i_txt}) and q (${s.q_txt}), the mass is ` +
    `${s.m_med} M☉ (68%: ${s.m_lo}–${s.m_hi}). BlackCAT lists ${s.m1_pub}.` : ` Its inclination is not bounded from below, so only a minimum mass is known (BlackCAT: ${s.m1_pub}).`;
  $("verdict").textContent = lim + meas + (s.refs ? ` Sources: ${s.refs}.` : "");
}
sel.onchange = () => { const s = SYS[+sel.value]; if (s.i_lo !== undefined) $("inc").value = ((s.i_lo + s.i_hi) / 2).toFixed(1); drawWeigh(); };
$("inc").oninput = drawWeigh;

// ---------------------------------------------------------------- outburst
const N = OB.mjd.length; $("day").max = N - 1;
const t0 = OB.mjd[0];
function iso(mjd) { return new Date((mjd - 40587) * 86400000).toISOString().slice(0, 10); }
function axes(svg, xs, ys, xl, yl, xt, yt) {
  for (const v of yt) { el("line", {x1: xs(xt[0]), x2: xs(xt[xt.length - 1]), y1: ys(v), y2: ys(v), class: "gridl"}, svg);
    el("text", {x: 40, y: ys(v) + 4, "text-anchor": "end"}, svg, v); }
  for (const v of xt) el("text", {x: xs(v), y: 282, "text-anchor": "middle"}, svg, typeof v === "number" && v > 1000 ? iso(v).slice(5) : v);
  el("text", {x: 260, y: 298, "text-anchor": "middle"}, svg, xl);
  el("text", {x: 10, y: 140, transform: "rotate(-90 10 140)", "text-anchor": "middle"}, svg, yl);
}
function drawOB() {
  const k = +$("day").value, lc = $("lc"), hid = $("hid");
  lc.innerHTML = ""; hid.innerHTML = "";
  const lx = m => 48 + (m - t0) / (OB.mjd[N - 1] - t0) * 420, ly = f => 10 + (Math.log10(5) - Math.log10(Math.max(f, 0.02))) / Math.log10(250) * 256;
  axes(lc, lx, ly, "2018 (month-day)", "MAXI 2–20 keV [Crab]", [t0, t0 + 60, t0 + 120, t0 + 180, t0 + 240].map(v => Math.round(v)), [0.02, 0.1, 0.5, 1, 3]);
  let d = ""; OB.mjd.forEach((m, j) => d += (j ? "L" : "M") + lx(m).toFixed(1) + "," + ly(OB.flux[j]).toFixed(1));
  el("path", {d, fill: "none", stroke: css("--axis"), "stroke-width": 1.2}, lc);
  for (const tt of [OB.h2s, OB.s2h]) el("line", {x1: lx(tt), x2: lx(tt), y1: 10, y2: 266, stroke: css("--ours"), "stroke-width": 1}, lc);
  el("circle", {cx: lx(OB.mjd[k]), cy: ly(OB.flux[k]), r: 6, fill: css("--ours"), stroke: css("--card"), "stroke-width": 2}, lc);
  const hx = h => 48 + (Math.log10(h) - Math.log10(0.12)) / (Math.log10(2.5) - Math.log10(0.12)) * 420;
  axes(hid, hx, ly, "hardness (4–10)/(2–4) keV [Crab = 1]", "MAXI 2–20 keV [Crab]", [0.15, 0.3, 0.6, 1, 2], [0.02, 0.1, 0.5, 1, 3]);
  let p = ""; for (let j = 0; j <= k; j++) p += (j ? "L" : "M") + hx(OB.hr[j]).toFixed(1) + "," + ly(OB.flux[j]).toFixed(1);
  el("path", {d: p, fill: "none", stroke: css("--curve"), "stroke-width": 1.4, opacity: 0.8}, hid);
  el("circle", {cx: hx(OB.hr[k]), cy: ly(OB.flux[k]), r: 6, fill: css("--ours"), stroke: css("--card"), "stroke-width": 2}, hid);
  const m = OB.mjd[k];
  const st = m < OB.h2s - 4 ? "Hard state: a hot corona of electrons up-scatters the disc's light into hard X-rays; a steady radio jet." :
    m <= OB.h2s + 4 ? "The transition: within days the spectrum softens as the inner disc takes over; the jet switches off with a final ejection." :
    m < OB.s2h ? "Soft state: a thermal disc reaching close to the black hole dominates, cooling slowly as the flux falls." :
    "Back to the hard state at a flux several times lower than the one it left at (hysteresis), then fading towards quiescence.";
  $("dayv").textContent = `${iso(m)}  ·  ${OB.flux[k].toFixed(2)} Crab  ·  hardness ${OB.hr[k].toFixed(2)}`;
  $("state").textContent = st;
}
let timer = null;
$("day").oninput = drawOB;
$("play").onclick = () => {
  if (timer) { clearInterval(timer); timer = null; $("play").textContent = "▶ Play"; return; }
  $("play").textContent = "❚❚ Pause";
  timer = setInterval(() => { const s = $("day"); if (+s.value >= N - 1) s.value = 0; else s.value = +s.value + 1; drawOB(); }, 60);
};
sel.value = SYS.findIndex(s => s.name === "GS 2023+338"); sel.onchange(); drawOB();
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { drawWeigh(); drawOB(); });
</script>
</body></html>
"""


def build(results, art, path=None) -> list:
    path = path or FIGURE_DIR / "black_hole_explorer.html"
    html = (PAGE.replace("__SYSTEMS__", json.dumps(_system_data(results, art)))
                .replace("__OUTBURST__", json.dumps(_outburst_data(results, art))))
    path.write_text(html, encoding="utf-8")
    return [str(path)]
