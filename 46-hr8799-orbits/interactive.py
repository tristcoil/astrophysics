"""Offline, self-contained orbit explorer: figures/hr8799_orbit_explorer.html

    python interactive.py

A time slider (1990-2100) moves the four planets along the orbits fitted to
all the data (1998-2025); measured positions appear on the dates they were
taken; a spread of posterior orbits shows how well each future position is
known. Everything (positions on a time grid, the data) is embedded in the
page, so it works without a network connection.
"""

from __future__ import annotations

import json
import warnings

import numpy as np
from astropy.time import Time

import astrometry as A
import config as C
import fits
import plotstyle as S
import predict

warnings.filterwarnings("ignore")

YEARS = (1990.0, 2100.0)
STEP = 0.1          # years, median track
STEP_DRAWS = 0.5    # years, posterior draws
N_DRAWS = 24


def build(key: str = "ALL/near1248_coplanar", out_name: str = "hr8799_orbit_explorer.html") -> str:
    fit = fits.load(key)
    rng = np.random.default_rng(C.RANDOM_SEED)
    s = fit.equal[rng.choice(len(fit.equal), 400, replace=False)]
    yrs = np.round(np.arange(YEARS[0], YEARS[1] + 1e-9, STEP), 2)
    sky = predict.predicted_positions(fit.model, s, Time(yrs, format="decimalyear").mjd)
    yd = np.round(np.arange(YEARS[0], YEARS[1] + 1e-9, STEP_DRAWS), 2)
    sky_d = predict.predicted_positions(fit.model, s[:N_DRAWS], Time(yd, format="decimalyear").mjd)
    tab = fit.table()
    plx = float(np.median(tab.plx))
    # the drawn orbit outlines: the one posterior orbit whose 1998-2026 track
    # lies closest to the median track (a pointwise median of orbits with
    # different periods is not itself an orbit, and looks wiggly)
    span = (yrs >= 1998) & (yrs <= 2026)
    dist = sum(((sky[p][0][:, span] - np.median(sky[p][0][:, span], 0)) ** 2
                + (sky[p][1][:, span] - np.median(sky[p][1][:, span], 0)) ** 2).sum(1) for p in C.PLANETS)
    rep = int(np.argmin(dist))
    phys_rep = fit.model.physical(s[rep:rep + 1])
    planets = {}
    for p in C.PLANETS:
        ra, dec = sky[p]
        P = float(np.median(tab[f"P_{p}"]))
        a_, e_, inc_, aop_, pan_, tau_ = (float(np.ravel(v)[0]) for v in phys_rep["elements"][p])
        import kepler as K
        m_k = float(np.ravel(phys_rep["mstar"])[0]) + sum(float(np.ravel(phys_rep["masses"][q])[0])
                                                          for q in C.INNER_TO_OUTER[: C.INNER_TO_OUTER.index(p) + 1])
        P_rep = float(K.period_days(a_, m_k))
        t_full = Time(2000.0, format="decimalyear").mjd + np.linspace(0, P_rep, 361)
        ell = predict.predicted_positions(fit.model, s[rep:rep + 1], t_full)[p]
        planets[p] = {
            "track": np.round(np.stack([np.median(ra, 0), np.median(dec, 0)], 1), 1).tolist(),
            "sigma": np.round(np.sqrt(0.5 * (np.var(ra, 0) + np.var(dec, 0))), 1).tolist(),
            "draws": np.round(np.stack([sky_d[p][0], sky_d[p][1]], -1), 1).tolist(),
            "ellipse": np.round(np.stack([ell[0][0], ell[1][0]], 1), 1).tolist(),
            "period": round(P, 1), "a": round(float(np.median(tab[f"a_{p}"])), 1),
            "color": S.PLANET_COLOR[p], "color_dark": S.PLANET_COLOR_DARK[p],
        }
    data = A.load(C.ALL_SOURCES + ("AMI-N26",), planets=C.PLANETS)
    pts = [{"p": r.planet, "y": round(r.year, 3), "x": round(r.ra, 1), "d": round(r.dec, 1),
            "s": r.source, "date": r.date if isinstance(r.date, str) else ""} for _, r in data.iterrows()]
    payload = {"years": [YEARS[0], YEARS[1], STEP], "years_draws": [YEARS[0], YEARS[1], STEP_DRAWS],
               "planets": planets, "data": pts, "plx": round(plx, 3), "last_data": round(float(data.year.max()), 2),
               "fit": key}
    html = TEMPLATE.replace("/*DATA*/", json.dumps(payload, separators=(",", ":")))
    path = C.FIGURES / out_name
    path.write_text(html)
    return str(path)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HR 8799 orbit explorer</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --surface-2: #f3f2ee; --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --accent: #2a78d6;
  --pe: #2a78d6; --pd: #eb6834; --pc: #1baf7a; --pb: #eda100;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --surface-2: #242423; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #8f8e88;
    --grid: #2f2f2d; --axis: #4a4a47; --accent: #3987e5;
    --pe: #3987e5; --pd: #d95926; --pc: #199e70; --pb: #c98500;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --surface-2: #242423; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #8f8e88;
  --grid: #2f2f2d; --axis: #4a4a47; --accent: #3987e5;
  --pe: #3987e5; --pd: #d95926; --pc: #199e70; --pb: #c98500;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width: 980px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 22px; margin: 0 0 4px; }
.sub { color: var(--ink-2); margin: 0 0 16px; font-size: 14px; }
.controls { display: flex; flex-wrap: wrap; gap: 10px 16px; align-items: center; margin-bottom: 12px; }
button, select { font: inherit; color: var(--ink); background: var(--surface-2); border: 1px solid var(--axis);
  border-radius: 6px; padding: 6px 12px; cursor: pointer; }
button:focus-visible, select:focus-visible, input:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
input[type=range] { flex: 1 1 240px; accent-color: var(--accent); }
.year { font-size: 28px; font-weight: 600; min-width: 90px; font-variant-numeric: tabular-nums; }
.mode { color: var(--ink-2); font-size: 13px; }
.wrap { display: grid; grid-template-columns: minmax(0, 1fr) 260px; gap: 16px; }
@media (max-width: 760px) { .wrap { grid-template-columns: 1fr; } }
canvas { width: 100%; aspect-ratio: 1 / 1; background: var(--surface); border: 1px solid var(--grid); border-radius: 8px;
  touch-action: none; }
.panel { font-size: 13px; }
.panel h2 { font-size: 14px; margin: 0 0 6px; }
.row { display: grid; grid-template-columns: 18px 1fr; gap: 6px; padding: 6px 0; border-bottom: 1px solid var(--grid); }
.dot { width: 12px; height: 12px; border-radius: 50%; margin-top: 4px; }
.num { font-variant-numeric: tabular-nums; }
.note { color: var(--muted); font-size: 12px; margin-top: 12px; }
#tip { position: fixed; pointer-events: none; background: var(--surface-2); color: var(--ink); border: 1px solid var(--axis);
  border-radius: 6px; padding: 6px 8px; font-size: 12px; display: none; max-width: 240px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; margin-top: 8px; }
th, td { text-align: right; padding: 4px 6px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
details { margin-top: 14px; }
label { color: var(--ink-2); font-size: 13px; display: inline-flex; gap: 6px; align-items: center; }
</style>
</head>
<body>
<main>
  <h1>HR 8799: four planets orbiting their star</h1>
  <p class="sub">Orbits fitted to every precise measurement from 1998 to 2025. Drag the slider or press play.
    East is to the left and north up, as on the sky.</p>
  <div class="controls">
    <button id="play" aria-label="Play or pause">Play</button>
    <input id="t" type="range" min="1990" max="2100" step="0.1" value="2026.8" aria-label="Year">
    <span class="year num" id="yr">2026.8</span>
    <select id="speed" aria-label="Speed"><option value="1">1 yr/s</option><option value="4" selected>4 yr/s</option>
      <option value="12">12 yr/s</option></select>
    <label><input type="checkbox" id="showdata" checked> measured positions</label>
    <label><input type="checkbox" id="showcloud" checked> uncertainty</label>
    <button id="theme" aria-label="Toggle light or dark">Theme</button>
  </div>
  <div class="mode" id="mode"></div>
  <div class="wrap">
    <canvas id="c" width="900" height="900" role="img" aria-label="Sky plane with the four planets"></canvas>
    <div class="panel">
      <h2>Where the planets are</h2>
      <div id="rows"></div>
      <p class="note">Separation in milliarcseconds (mas) and in astronomical units at the Gaia distance
        (40.9 pc). ± is the 1σ spread of the orbits allowed by the data. The orbits come from the
        measured positions; the planets' masses barely affect them.</p>
    </div>
  </div>
  <details>
    <summary>Table view</summary>
    <table id="tbl"></table>
  </details>
  <p class="note" id="foot"></p>
</main>
<div id="tip" role="tooltip"></div>
<script>
const D = /*DATA*/;
const P = ["e", "d", "c", "b"];
const cv = document.getElementById("c"), ctx = cv.getContext("2d");
const css = () => getComputedStyle(document.documentElement);
let t = 2026.8, playing = false, last = null;
const lim = 1950;
const N = cv.width;
const sx = x => N / 2 - (x / lim) * (N / 2 - 30);     // east to the left
const sy = y => N / 2 - (y / lim) * (N / 2 - 30);
function idx(y, a) { return Math.max(0, Math.min(Math.round((y - a[0]) / a[2]), Math.round((a[1] - a[0]) / a[2]))); }
function col(p) { return css().getPropertyValue("--p" + p).trim(); }
function draw() {
  const c = css(); ctx.clearRect(0, 0, N, N);
  ctx.fillStyle = c.getPropertyValue("--surface"); ctx.fillRect(0, 0, N, N);
  // grid rings every 500 mas
  ctx.strokeStyle = c.getPropertyValue("--grid"); ctx.lineWidth = 1;
  for (let r = 500; r <= 2000; r += 500) { ctx.beginPath(); ctx.arc(N / 2, N / 2, r / lim * (N / 2 - 30), 0, 2 * Math.PI); ctx.stroke(); }
  ctx.fillStyle = c.getPropertyValue("--muted"); ctx.font = "20px system-ui";
  ctx.fillText("N", N / 2 - 6, 24); ctx.fillText("E", 10, N / 2 + 6);
  for (const p of P) {
    const pl = D.planets[p], e = pl.ellipse;
    ctx.strokeStyle = col(p); ctx.globalAlpha = 0.35; ctx.lineWidth = 2; ctx.beginPath();
    e.forEach((q, k) => k ? ctx.lineTo(sx(q[0]), sy(q[1])) : ctx.moveTo(sx(q[0]), sy(q[1]))); ctx.stroke();
    ctx.globalAlpha = 1;
  }
  if (document.getElementById("showdata").checked) {
    ctx.fillStyle = c.getPropertyValue("--ink");
    for (const d of D.data) if (d.y <= t) { ctx.beginPath(); ctx.arc(sx(d.x), sy(d.d), 3.2, 0, 2 * Math.PI); ctx.fill(); }
  }
  const i = idx(t, D.years), j = idx(t, D.years_draws);
  for (const p of P) {
    const pl = D.planets[p];
    if (document.getElementById("showcloud").checked) {
      ctx.fillStyle = col(p); ctx.globalAlpha = 0.35;
      for (const dr of pl.draws) { const q = dr[j]; ctx.beginPath(); ctx.arc(sx(q[0]), sy(q[1]), 2.4, 0, 2 * Math.PI); ctx.fill(); }
      ctx.globalAlpha = 1;
    }
    // trail: last 6 years
    ctx.strokeStyle = col(p); ctx.lineWidth = 4; ctx.lineCap = "round"; ctx.beginPath();
    const k0 = Math.max(0, i - 60);
    for (let k = k0; k <= i; k++) { const q = pl.track[k]; k === k0 ? ctx.moveTo(sx(q[0]), sy(q[1])) : ctx.lineTo(sx(q[0]), sy(q[1])); }
    ctx.stroke();
    const q = pl.track[i];
    ctx.fillStyle = c.getPropertyValue("--surface"); ctx.beginPath(); ctx.arc(sx(q[0]), sy(q[1]), 13, 0, 2 * Math.PI); ctx.fill();
    ctx.fillStyle = col(p); ctx.beginPath(); ctx.arc(sx(q[0]), sy(q[1]), 10, 0, 2 * Math.PI); ctx.fill();
    const r = Math.hypot(q[0], q[1]);
    ctx.fillStyle = c.getPropertyValue("--ink"); ctx.font = "bold 24px system-ui";
    ctx.fillText(p, sx(q[0] * (1 + 90 / r)) - 7, sy(q[1] * (1 + 90 / r)) + 8);
  }
  // star
  ctx.fillStyle = c.getPropertyValue("--ink"); ctx.font = "34px system-ui"; ctx.fillText("★", N / 2 - 15, N / 2 + 12);
  // scale bar 20 au
  const L = 20 * D.plx / lim * (N / 2 - 30);
  ctx.strokeStyle = c.getPropertyValue("--ink-2"); ctx.lineWidth = 3;
  ctx.beginPath(); ctx.moveTo(N - 40 - L, N - 34); ctx.lineTo(N - 40, N - 34); ctx.stroke();
  ctx.fillStyle = c.getPropertyValue("--ink-2"); ctx.font = "18px system-ui"; ctx.fillText("20 au", N - 40 - L / 2 - 22, N - 44);
  panel(i);
}
function panel(i) {
  document.getElementById("yr").textContent = t.toFixed(1);
  document.getElementById("mode").textContent = t <= D.last_data ? "Inside the span of the measurements: the orbit is fitted here."
    : "After the last measurement (Keck, Oct 2025): this is a prediction.";
  let html = "", tb = "<tr><th>planet</th><th>ΔRA (mas)</th><th>ΔDec (mas)</th><th>sep (mas)</th><th>± (mas)</th><th>PA (°)</th></tr>";
  for (const p of ["b", "c", "d", "e"]) {
    const pl = D.planets[p], q = pl.track[i], s = pl.sigma[i];
    const sep = Math.hypot(q[0], q[1]), pa = (Math.atan2(q[0], q[1]) * 180 / Math.PI + 360) % 360;
    html += `<div class="row"><span class="dot" style="background:${col(p)}"></span><div><b>${p}</b> &nbsp;
      <span class="num">${sep.toFixed(0)} ± ${s.toFixed(1)} mas = ${(sep / D.plx).toFixed(1)} au, PA ${pa.toFixed(1)}°</span><br>
      <span style="color:var(--muted)">orbit a ≈ ${pl.a} au, period ≈ ${pl.period.toFixed(0)} yr</span></div></div>`;
    tb += `<tr><td>${p}</td><td>${q[0].toFixed(1)}</td><td>${q[1].toFixed(1)}</td><td>${sep.toFixed(1)}</td><td>${s.toFixed(1)}</td><td>${pa.toFixed(2)}</td></tr>`;
  }
  document.getElementById("rows").innerHTML = html;
  document.getElementById("tbl").innerHTML = tb;
}
function tick(ts) {
  if (!playing) return;
  if (last !== null) { t += (ts - last) / 1000 * parseFloat(document.getElementById("speed").value); }
  last = ts;
  if (t > 2100) t = 1990;
  document.getElementById("t").value = t; draw(); requestAnimationFrame(tick);
}
document.getElementById("play").onclick = () => { playing = !playing; last = null;
  document.getElementById("play").textContent = playing ? "Pause" : "Play"; if (playing) requestAnimationFrame(tick); };
document.getElementById("t").oninput = e => { t = parseFloat(e.target.value); draw(); };
for (const id of ["showdata", "showcloud"]) document.getElementById(id).onchange = draw;
document.getElementById("theme").onclick = () => {
  const r = document.documentElement, dark = getComputedStyle(r).getPropertyValue("--surface").trim() === "#1a1a19";
  r.setAttribute("data-theme", dark ? "light" : "dark"); draw(); };
const tip = document.getElementById("tip");
cv.addEventListener("pointermove", ev => {
  if (!document.getElementById("showdata").checked) { tip.style.display = "none"; return; }
  const b = cv.getBoundingClientRect(), k = N / b.width, mx = (ev.clientX - b.left) * k, my = (ev.clientY - b.top) * k;
  let best = null, bd = 18 * k;
  for (const d of D.data) { if (d.y > t) continue; const dd = Math.hypot(sx(d.x) - mx, sy(d.d) - my); if (dd < bd) { bd = dd; best = d; } }
  if (best) { tip.style.display = "block"; tip.style.left = (ev.clientX + 12) + "px"; tip.style.top = (ev.clientY + 12) + "px";
    tip.innerHTML = `<b>planet ${best.p}</b> · ${best.s}<br>${best.date || best.y.toFixed(2)}<br>ΔRA ${best.x} mas, ΔDec ${best.d} mas`; }
  else tip.style.display = "none";
});
cv.addEventListener("pointerleave", () => tip.style.display = "none");
document.getElementById("foot").textContent = `Orbits: ${D.fit} fit (nested sampling, this project). Data: Keck/NIRC2, Gemini/GPI, VLT/SPHERE, LBT/LUCI, VLTI/GRAVITY, JWST, HST/NICMOS (see README). Measured positions from JWST/AMI (10-30 mas errors) are shown but were not fitted.`;
draw();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    print(build())
