"""Offline, self-contained explorer: figures/triple_star_explorer.html

    python interactive.py

Two views, everything embedded (works without a network connection):

* Eclipses that switch: pick SS Lac, HS Hya or V907 Sco and drag the year
  (1880-2300). Left: the eclipsing pair's orbit as seen from Earth that year
  (to scale: the two stars and their projected orbit), right: the tilt history
  with the plate and TESS measurements, and that year's light curve from the
  fitted model.
* The triples atlas: every inner/outer pair of the Multiple Star Catalog,
  outer against inner period, hover for the system.
"""

from __future__ import annotations

import json

import numpy as np

import config as C


def switch_data(name, n_draws=120):
    import figures as FG
    import switch_fit as F
    import switching as S
    draws, final = FG.posterior_draws(name, n=n_draws)
    m = S.Model(name)
    y0, y1 = 1880, (2300 if name != "V907 Sco" else 2100)
    yrs = np.arange(y0, y1 + 1, 1.0)
    jd = C.year_to_jd(yrs)
    tilts, J_d, m3_d, Pn_d = [], [], [], []
    for th in draws:
        g = m.geometry(m.unpack(th))
        tilts.append(FG.folded_tilt(np.degrees(m.inclination(g, jd))))
        J_d.append(np.degrees(g["J"]))
        m3_d.append(g["m3"])
        Pn_d.append(g["P_node_yr"])
    tilts = np.array(tilts)
    q = np.percentile(tilts, [16, 50, 84], axis=0)
    best = draws[0]
    g = m.geometry(m.unpack(best))
    sp = m.spec
    # orientation of the projected orbit each year: normal vector -> tilt and
    # position angle of the line of nodes (rotating about the line of sight
    # does not matter for eclipses, but it makes the drawing honest)
    h = g["h0"] + g["nu"] * (jd - sp.t_ref)
    eps, I = g["eps_in"], g["I"]
    nx = np.sin(eps) * np.sin(h)
    ny = np.cos(I) * (-np.sin(eps) * np.cos(h)) - np.sin(I) * np.cos(eps)
    nz = np.sin(I) * (-np.sin(eps) * np.cos(h)) + np.cos(I) * np.cos(eps)
    inc = np.degrees(np.arccos(np.clip(nz, -1, 1)))
    node = np.degrees(np.arctan2(nx, -ny))
    # light curves (plate band) every 5 years over one orbit, phase from the
    # primary conjunction
    lc_years = np.arange(y0, y1 + 1, 5.0)
    phases = np.linspace(-0.25, 0.75, 161)
    P = 2 * np.pi / g["n_lam"]
    lcs = []
    for yy in lc_years:
        tj = float(C.year_to_jd(yy))
        tc = F.GlobalFit.conjunction_time(_FakeFit(sp), g, tj)
        f = m.flux(g, tc + phases * P, band="B", L3=F.L3_B[name])
        lcs.append(np.round(-2.5 * np.log10(f), 4).tolist())
    fit = F.GlobalFit(name)
    eras = F.era_inclinations(name, best)
    tess = []
    rs = g["r1"] + g["r2"]
    for pr in F.profiles(name):
        x = np.array(pr["x"])
        c = (np.array(pr["chi2"]) - pr["chi2_min"]) / pr["scale"]
        ti = lambda xx: float(np.degrees(np.arccos(np.clip(xx * rs, 0, 1))))
        if pr["detected"]:
            ok = c < 1
            tess.append(dict(year=float(C.jd_to_year(pr["jd_mid"])), i=ti(pr["x_best"]), lo=ti(x[ok].max()),
                             hi=ti(x[ok].min()), det=True, sector=pr["sector"]))
        else:
            tess.append(dict(year=float(C.jd_to_year(pr["jd_mid"])), i=ti(x[c < 4].min()), det=False,
                             sector=pr["sector"]))
    ground = []
    gpath = F.CACHE / f"ground_profiles_{F._slug(name)}.json"
    for pr in (json.loads(gpath.read_text()) if gpath.exists() else []):
        x = np.array(pr["x"])
        c = (np.array(pr["chi2"]) - pr["chi2_min"]) / pr["scale"]
        ti = lambda xx: float(np.degrees(np.arccos(np.clip(xx * rs, 0, 1))))
        ok = c < 1
        measured = bool(pr["detected"] or x[ok].max() < x.max() - 0.05)
        ground.append(dict(year=float(C.jd_to_year(pr["jd_mid"])), label=pr["label"], det=measured,
                           i=ti(pr["x_best"]) if measured else ti(x[c < 4].min()),
                           lo=ti(x[ok].max()) if measured else None, hi=ti(x[ok].min()) if measured else None))
    return dict(name=name, final=bool(final), years=yrs.tolist(), q16=np.round(q[0], 2).tolist(),
                q50=np.round(q[1], 2).tolist(), q84=np.round(q[2], 2).tolist(),
                inc=np.round(inc, 2).tolist(), node=np.round(node, 2).tolist(),
                r1=float(g["r1"]), r2=float(g["r2"]), e=float(sp.e), omega_deg=float(np.degrees(g["g0"])),
                i_edge=float(np.degrees(np.arccos(rs))), lc_years=lc_years.tolist(), phases=phases.tolist(),
                lcs=lcs, eras=eras, tess=tess, ground=ground,
                P_node_yr=float(np.median(Pn_d)), J=float(np.median(J_d)), m3=float(np.median(m3_d)))


class _FakeFit:
    """Minimal stand-in so GlobalFit.conjunction_time can be used unbound."""

    def __init__(self, sp):
        self.sp = sp


def atlas_data():
    import architecture as A
    comp, sys_, orb, notes, pairs = A.msc_pairs()
    st = A.stability_table(pairs)
    margin = dict(zip(zip(st["wds"], st["inner"], st["outer"]), st["margin"]))
    out = []
    ok = np.isfinite(pairs["P_in"]) & np.isfinite(pairs["P_out"])
    for _, r in pairs[ok].iterrows():
        mg = margin.get((r["wds"], r["inner"], r["outer"]))
        out.append([r["wds"], round(float(np.log10(r["P_in"])), 3), round(float(np.log10(r["P_out"])), 3),
                    int(bool(r["inner_from_orbit"] and r["outer_from_orbit"])),
                    (r["ids"] or "")[:48], r["inner"], r["outer"],
                    None if mg is None or not np.isfinite(mg) else round(float(mg), 2)])
    return out


HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Triple star explorer</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --surface-2: #f2f1ed; --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --accent: #2a78d6; --violet: #4a3aa7;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --star2: #eb6834;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --surface-2: #242423; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #8f8e88;
    --grid: #2f2f2d; --axis: #4a4a47; --accent: #3987e5; --violet: #9085e9;
    --s1: #3987e5; --s2: #d95926; --s3: #199e70; --star2: #d95926;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --surface-2: #242423; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #8f8e88;
  --grid: #2f2f2d; --axis: #4a4a47; --accent: #3987e5; --violet: #9085e9;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --star2: #d95926;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width: 1080px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 22px; margin: 0 0 4px; }
h2 { font-size: 17px; margin: 26px 0 4px; }
.sub { color: var(--ink-2); margin: 0 0 14px; font-size: 14px; }
.controls { display: flex; flex-wrap: wrap; gap: 10px 16px; align-items: center; margin-bottom: 12px; }
button, select { font: inherit; color: var(--ink); background: var(--surface-2); border: 1px solid var(--axis);
  border-radius: 6px; padding: 6px 12px; cursor: pointer; }
button:focus-visible, select:focus-visible, input:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
input[type=range] { flex: 1 1 240px; accent-color: var(--accent); }
.year { font-size: 28px; font-weight: 600; min-width: 80px; font-variant-numeric: tabular-nums; }
.wrap { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1.25fr); gap: 16px; }
@media (max-width: 760px) { .wrap { grid-template-columns: 1fr; } }
canvas { width: 100%; background: var(--surface); border: 1px solid var(--grid); border-radius: 8px; touch-action: none; }
#sky { aspect-ratio: 1 / 1; }
#hist { aspect-ratio: 16 / 9; }
#lc { aspect-ratio: 16 / 7; margin-top: 10px; }
#atlas { aspect-ratio: 16 / 10; }
.status { font-size: 14px; margin: 6px 0 0; font-variant-numeric: tabular-nums; }
.note { color: var(--muted); font-size: 12px; margin-top: 10px; }
#tip { position: fixed; pointer-events: none; background: var(--surface-2); color: var(--ink); border: 1px solid var(--axis);
  border-radius: 6px; padding: 6px 8px; font-size: 12px; display: none; max-width: 280px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; margin-top: 8px; }
th, td { text-align: right; padding: 4px 6px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
details { margin-top: 12px; }
label { color: var(--ink-2); font-size: 13px; display: inline-flex; gap: 6px; align-items: center; }
</style>
</head>
<body>
<main>
  <h1>Triple stars: eclipses that switch on and off</h1>
  <p class="sub">A third star makes a close pair's orbit wobble like a spinning top. When the tilt passes
    edge-on, the pair eclipses; when it moves away, the eclipses stop. Drag the year.</p>
  <div class="controls">
    <select id="sys" aria-label="System"><option>SS Lac</option><option>HS Hya</option><option selected>V907 Sco</option></select>
    <button id="play" aria-label="Play or pause">Play</button>
    <input id="t" type="range" min="1880" max="2300" step="1" value="2026" aria-label="Year">
    <span class="year" id="yr">2026</span>
    <button id="theme" aria-label="Toggle light or dark">Theme</button>
  </div>
  <div class="wrap">
    <div>
      <canvas id="sky" width="700" height="700" role="img" aria-label="The pair's orbit as seen from Earth"></canvas>
      <p class="status" id="status"></p>
    </div>
    <div>
      <canvas id="hist" width="900" height="506" role="img" aria-label="Tilt of the orbit against year"></canvas>
      <canvas id="lc" width="900" height="394" role="img" aria-label="Light curve that year"></canvas>
    </div>
  </div>
  <p class="note" id="foot"></p>
  <details>
    <summary>Table view (measurements)</summary>
    <table id="tbl"></table>
  </details>
  <h2>The triples atlas</h2>
  <p class="sub">Every inner and outer orbit pair in Tokovinin's Multiple Star Catalog (2026). Hover a point.
    Lines: the Mardling-Aarseth stability limit for equal masses (circular and e = 0.6 outer orbits).</p>
  <div class="controls">
    <label><input type="checkbox" id="onlyorb"> only pairs with both orbits known</label>
  </div>
  <canvas id="atlas" width="1000" height="625" role="img" aria-label="Outer against inner period for every triple"></canvas>
</main>
<div id="tip" role="tooltip"></div>
<script>
const DATA = __DATA__;
const ATLAS = __ATLAS__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const SYSCOL = {"SS Lac": "--s1", "HS Hya": "--s2", "V907 Sco": "--s3"};
let sysName = "V907 Sco", year = 2026, playing = false, last = null;
const $ = id => document.getElementById(id);

function setupCanvas(cv) { const ctx = cv.getContext("2d"); return ctx; }

function axes(ctx, W, H, pad, xr, yr, xlab, ylab, xt, yt) {
  const X = x => pad.l + (x - xr[0]) / (xr[1] - xr[0]) * (W - pad.l - pad.r);
  const Y = y => H - pad.b - (y - yr[0]) / (yr[1] - yr[0]) * (H - pad.t - pad.b);
  ctx.strokeStyle = css("--grid"); ctx.lineWidth = 1; ctx.fillStyle = css("--muted");
  ctx.font = "24px system-ui"; ctx.textAlign = "center";
  for (const v of xt) { ctx.beginPath(); ctx.moveTo(X(v), Y(yr[0])); ctx.lineTo(X(v), Y(yr[1])); ctx.stroke();
    ctx.fillText(v, X(v), H - pad.b + 30); }
  ctx.textAlign = "right";
  for (const v of yt) { ctx.beginPath(); ctx.moveTo(X(xr[0]), Y(v)); ctx.lineTo(X(xr[1]), Y(v)); ctx.stroke();
    ctx.fillText(v, pad.l - 10, Y(v) + 8); }
  ctx.fillStyle = css("--ink-2"); ctx.textAlign = "center"; ctx.font = "26px system-ui";
  ctx.fillText(xlab, (W + pad.l - pad.r) / 2, H - 10);
  ctx.save(); ctx.translate(26, (H - pad.b + pad.t) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(ylab, 0, 0); ctx.restore();
  return [X, Y];
}

function interp(xs, ys, x) {
  if (x <= xs[0]) return ys[0]; if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
  let i = 0; while (xs[i + 1] < x) i++;
  const f = (x - xs[i]) / (xs[i + 1] - xs[i]); return ys[i] + f * (ys[i + 1] - ys[i]);
}

function drawSky() {
  const D = DATA[sysName], cv = $("sky"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height;
  ctx.fillStyle = css("--surface"); ctx.fillRect(0, 0, W, H);
  const inc = interp(D.years, D.inc, year) * Math.PI / 180, node = interp(D.years, D.node, year) * Math.PI / 180;
  const S = 0.36 * W, cx = W / 2, cy = H / 2;
  // projected orbit of the secondary around the primary (circular approximation for the drawing)
  const pts = [];
  for (let k = 0; k <= 240; k++) {
    const u = 2 * Math.PI * k / 240;
    const x = Math.cos(u), y = Math.sin(u) * Math.cos(inc), z = Math.sin(u) * Math.sin(inc);
    const X = x * Math.cos(node) - y * Math.sin(node), Y = x * Math.sin(node) + y * Math.cos(node);
    pts.push([cx + S * X, cy - S * Y, z]);
  }
  ctx.strokeStyle = css("--axis"); ctx.lineWidth = 2; ctx.beginPath();
  pts.forEach((p, i) => i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])); ctx.stroke();
  // animate the secondary along the orbit (a few seconds per orbit, for show)
  const u = (performance.now() / 3000) % 1 * 2 * Math.PI;
  const x = Math.cos(u), y = Math.sin(u) * Math.cos(inc), z = Math.sin(u) * Math.sin(inc);
  const X = x * Math.cos(node) - y * Math.sin(node), Y = x * Math.sin(node) + y * Math.cos(node);
  const p1 = [cx, cy, 0], p2 = [cx + S * X, cy - S * Y, z];
  const r1 = D.r1 * S, r2 = D.r2 * S;
  const order = p2[2] > 0 ? [[p1, r1, "--s1x"], [p2, r2, "--star2"]] : [[p2, r2, "--star2"], [p1, r1, "--s1x"]];
  for (const [p, r, c] of order) {
    ctx.fillStyle = c === "--star2" ? css("--star2") : css(SYSCOL[sysName]);
    ctx.strokeStyle = css("--surface"); ctx.lineWidth = 3;
    ctx.beginPath(); ctx.arc(p[0], p[1], r, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
  }
  ctx.fillStyle = css("--ink-2"); ctx.font = "26px system-ui"; ctx.textAlign = "left";
  ctx.fillText("the pair as seen from Earth (to scale)", 18, 36);
  const tilt = 90 - Math.abs(90 - inc * 180 / Math.PI);
  const ecl = tilt > D.i_edge;
  $("status").textContent = `${Math.round(year)}: tilt ${tilt.toFixed(1)}° — ` +
    (ecl ? "the stars eclipse each other every orbit" : "the stars miss each other: no eclipses");
}

function drawHist() {
  const D = DATA[sysName], cv = $("hist"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height;
  ctx.fillStyle = css("--surface"); ctx.fillRect(0, 0, W, H);
  const xr = [D.years[0], D.years[D.years.length - 1]];
  const xt = []; for (let v = Math.ceil(xr[0] / 50) * 50; v <= xr[1]; v += 50) xt.push(v);
  const pad = {l: 96, r: 40, t: 16, b: 74};
  const [X, Y] = axes(ctx, W, H, pad, xr, [60, 90.5], "year", "tilt (90° = edge-on)", xt, [60, 70, 80, 90]);
  const col = css(SYSCOL[sysName]);
  ctx.fillStyle = col; ctx.globalAlpha = 0.10; ctx.fillRect(X(xr[0]), Y(90), X(xr[1]) - X(xr[0]), Y(D.i_edge) - Y(90));
  ctx.globalAlpha = 0.28; ctx.beginPath();
  D.years.forEach((y, i) => i ? ctx.lineTo(X(y), Y(Math.max(60, D.q84[i]))) : ctx.moveTo(X(y), Y(Math.max(60, D.q84[i]))));
  for (let i = D.years.length - 1; i >= 0; i--) ctx.lineTo(X(D.years[i]), Y(Math.max(60, D.q16[i])));
  ctx.closePath(); ctx.fill(); ctx.globalAlpha = 1;
  ctx.strokeStyle = col; ctx.lineWidth = 3; ctx.beginPath();
  D.years.forEach((y, i) => { const v = Math.max(59.5, D.q50[i]); i ? ctx.lineTo(X(y), Y(v)) : ctx.moveTo(X(y), Y(v)); });
  ctx.stroke();
  ctx.fillStyle = css("--ink-2"); ctx.font = "22px system-ui"; ctx.textAlign = "right";
  ctx.fillText("eclipses in this band", X(xr[1]) - 8, Y(90) + 28);
  // measurements
  ctx.fillStyle = css("--ink"); ctx.strokeStyle = css("--ink"); ctx.lineWidth = 2;
  for (const e of D.eras) {
    const yc = e.year;   // an era's tilt refers to its median plate date
    if (e.dlnl_vs_none > 2) {
      ctx.beginPath(); ctx.moveTo(X(e.lo), Y(e.i_best)); ctx.lineTo(X(e.hi), Y(e.i_best));
      ctx.moveTo(X(yc), Y(e.i_lo)); ctx.lineTo(X(yc), Y(e.i_hi)); ctx.stroke();
      ctx.beginPath(); ctx.arc(X(yc), Y(e.i_best), 6, 0, 2 * Math.PI); ctx.fill();
    } else {
      ctx.beginPath(); ctx.moveTo(X(e.lo), Y(e.i_hi95)); ctx.lineTo(X(e.hi), Y(e.i_hi95)); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(X(yc) - 7, Y(e.i_hi95) - 2); ctx.lineTo(X(yc) + 7, Y(e.i_hi95) - 2);
      ctx.lineTo(X(yc), Y(e.i_hi95) + 10); ctx.closePath(); ctx.stroke();
    }
  }
  for (const s of D.ground) {
    ctx.lineWidth = 2;
    if (s.det) {
      ctx.beginPath(); ctx.moveTo(X(s.year), Y(s.lo)); ctx.lineTo(X(s.year), Y(s.hi)); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(X(s.year), Y(s.i) - 8); ctx.lineTo(X(s.year) + 8, Y(s.i));
      ctx.lineTo(X(s.year), Y(s.i) + 8); ctx.lineTo(X(s.year) - 8, Y(s.i)); ctx.closePath();
      ctx.fillStyle = css("--surface"); ctx.fill(); ctx.stroke(); ctx.fillStyle = css("--ink");
    } else { ctx.beginPath(); ctx.moveTo(X(s.year) - 6, Y(s.i) - 2); ctx.lineTo(X(s.year) + 6, Y(s.i) - 2);
      ctx.lineTo(X(s.year), Y(s.i) + 9); ctx.closePath(); ctx.stroke(); }
  }
  for (const s of D.tess) {
    if (s.det) { ctx.fillRect(X(s.year) - 6, Y(s.i) - 6, 12, 12); }
    else { ctx.beginPath(); ctx.moveTo(X(s.year) - 6, Y(s.i) - 2); ctx.lineTo(X(s.year) + 6, Y(s.i) - 2);
      ctx.lineTo(X(s.year), Y(s.i) + 9); ctx.closePath(); ctx.stroke(); }
  }
  ctx.strokeStyle = css("--accent"); ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(X(year), Y(60)); ctx.lineTo(X(year), Y(90.5)); ctx.stroke();
}

function drawLC() {
  const D = DATA[sysName], cv = $("lc"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height;
  ctx.fillStyle = css("--surface"); ctx.fillRect(0, 0, W, H);
  // nearest precomputed year
  let k = 0; for (let i = 0; i < D.lc_years.length; i++) if (Math.abs(D.lc_years[i] - year) < Math.abs(D.lc_years[k] - year)) k = i;
  const lc = D.lcs[k];
  const ymax = Math.max(0.3, Math.max(...D.lcs.map(a => Math.max(...a))) * 1.05);
  const pad = {l: 96, r: 24, t: 44, b: 74};
  const yt = [0, +(ymax / 2).toFixed(2), +ymax.toFixed(2)];
  const [X, Y] = axes(ctx, W, H, pad, [-0.25, 0.75], [ymax, -0.05], "orbital phase", "dimming (mag)",
    [0, 0.25, 0.5, 0.75], yt);
  ctx.strokeStyle = css(SYSCOL[sysName]); ctx.lineWidth = 3; ctx.beginPath();
  D.phases.forEach((p, i) => i ? ctx.lineTo(X(p), Y(lc[i])) : ctx.moveTo(X(p), Y(lc[i]))); ctx.stroke();
  ctx.fillStyle = css("--ink-2"); ctx.font = "24px system-ui"; ctx.textAlign = "left";
  ctx.fillText(`one orbit as seen in ${D.lc_years[k]} (model, blue-sensitive plates)`, pad.l, 30);
}

function table() {
  const D = DATA[sysName];
  let h = "<tr><th>Data</th><th>Years</th><th>Tilt (deg)</th><th>Note</th></tr>";
  for (const e of D.eras) h += `<tr><td>Harvard plates (${e.n})</td><td>${e.lo}-${e.hi}</td><td>` +
    (e.dlnl_vs_none > 2 ? `${e.i_best.toFixed(1)} (${e.i_lo.toFixed(1)}-${e.i_hi.toFixed(1)})` : `< ${e.i_hi95.toFixed(1)}`) +
    `</td><td>${e.dlnl_vs_none > 2 ? "eclipses seen" : "no eclipse"}</td></tr>`;
  for (const s of D.ground) h += `<tr><td>${s.label}</td><td>${s.year.toFixed(1)}</td><td>` +
    (s.det ? `${s.i.toFixed(1)} (${s.lo.toFixed(1)}-${s.hi.toFixed(1)})` : `< ${s.i.toFixed(1)}`) +
    `</td><td>${s.det ? "tilt measured" : "no eclipse"}</td></tr>`;
  for (const s of D.tess) h += `<tr><td>TESS sector ${s.sector}</td><td>${s.year.toFixed(2)}</td><td>` +
    (s.det ? s.i.toFixed(2) : `< ${s.i.toFixed(2)}`) + `</td><td>${s.det ? "eclipses seen" : "no eclipse"}</td></tr>`;
  $("tbl").innerHTML = h;
  $("foot").textContent = `${sysName} (posterior medians): mutual inclination ${D.J.toFixed(0)}°, third star ${D.m3.toFixed(2)} Msun, ` +
    `precession period ${Math.round(D.P_node_yr)} years. Band: 68% of fitted models. Dots: Harvard plate eras ` +
    `(at their median plate date); squares: TESS; open diamonds: Hipparcos and ASAS-3; triangles: no eclipse ` +
    `(upper limit).` + (D.final ? "" : " (preview)");
}

function drawAtlas(hover) {
  const cv = $("atlas"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height;
  ctx.fillStyle = css("--surface"); ctx.fillRect(0, 0, W, H);
  const pad = {l: 90, r: 16, t: 16, b: 74};
  const [X, Y] = axes(ctx, W, H, pad, [-1, 7.5], [0.5, 9.6], "inner period, log10 (days)", "outer period, log10 (days)",
    [-1, 0, 1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7, 8, 9]);
  const only = $("onlyorb").checked;
  for (const r of ATLAS) { if (only && !r[3]) continue;
    ctx.fillStyle = r[3] ? css("--s1") : css("--axis"); ctx.fillRect(X(r[1]) - 2, Y(r[2]) - 2, 4, 4); }
  ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.5;
  for (const off of [0.688, 1.527]) { ctx.beginPath(); ctx.moveTo(X(-1), Y(-1 + off)); ctx.lineTo(X(7.5), Y(7.5 + off)); ctx.stroke(); }
  if (hover) { ctx.strokeStyle = css("--ink"); ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(X(hover[1]), Y(hover[2]), 8, 0, 2 * Math.PI); ctx.stroke(); }
  return [X, Y];
}

function redraw() { drawSky(); drawHist(); drawLC(); }
function setSystem(n) {
  sysName = n; const D = DATA[n];
  $("t").min = D.years[0]; $("t").max = D.years[D.years.length - 1];
  if (year > D.years[D.years.length - 1]) { year = 2026; $("t").value = year; $("yr").textContent = year; }
  table(); redraw();
}
$("sys").onchange = e => setSystem(e.target.value);
$("t").oninput = e => { year = +e.target.value; $("yr").textContent = year; redraw(); };
$("play").onclick = () => { playing = !playing; $("play").textContent = playing ? "Pause" : "Play"; last = null; };
$("theme").onclick = () => { const r = document.documentElement;
  r.dataset.theme = (r.dataset.theme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light")) === "dark" ? "light" : "dark";
  redraw(); drawAtlas(); };
$("onlyorb").onchange = () => drawAtlas();
const tip = $("tip");
$("atlas").addEventListener("mousemove", ev => {
  const cv = $("atlas"), b = cv.getBoundingClientRect();
  const mx = (ev.clientX - b.left) * cv.width / b.width, my = (ev.clientY - b.top) * cv.height / b.height;
  const [X, Y] = drawAtlas(); let best = null, bd = 1e9; const only = $("onlyorb").checked;
  for (const r of ATLAS) { if (only && !r[3]) continue; const d = Math.hypot(X(r[1]) - mx, Y(r[2]) - my); if (d < bd) { bd = d; best = r; } }
  if (best && bd < 24) { drawAtlas(best); tip.style.display = "block"; tip.style.left = (ev.clientX + 14) + "px"; tip.style.top = (ev.clientY + 14) + "px";
    const pin = Math.pow(10, best[1]), pout = Math.pow(10, best[2]);
    const fmt = p => p < 1000 ? p.toFixed(p < 10 ? 2 : 0) + " d" : (p / 365.25 < 1e4 ? (p / 365.25).toFixed(0) + " yr" : (p / 365.25 / 1000).toFixed(0) + " kyr");
    tip.innerHTML = `<b>WDS ${best[0]}</b> ${best[4]}<br>inner ${best[5]}: ${fmt(pin)}<br>outer ${best[6]}: ${fmt(pout)}<br>` +
      (best[7] !== null ? `${best[7]}× the stability limit` : "masses unknown") + (best[3] ? "" : "<br>(a period estimated from separation)");
  } else { tip.style.display = "none"; }
});
$("atlas").addEventListener("mouseleave", () => { tip.style.display = "none"; drawAtlas(); });
function tick(ts) {
  if (playing) { if (last !== null) { year += (ts - last) / 1000 * 8; const D = DATA[sysName];
      if (year > D.years[D.years.length - 1]) year = D.years[0]; $("t").value = Math.round(year); $("yr").textContent = Math.round(year); }
    last = ts; }
  drawSky(); if (playing) { drawHist(); drawLC(); }
  requestAnimationFrame(tick);
}
setSystem(sysName); drawAtlas(); requestAnimationFrame(tick);
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { redraw(); drawAtlas(); });
</script>
</body>
</html>
"""


def build(out_name="triple_star_explorer.html"):
    data = {name: switch_data(name) for name in ("SS Lac", "HS Hya", "V907 Sco")}
    atlas = atlas_data()
    html = HTML.replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace(
        "__ATLAS__", json.dumps(atlas, separators=(",", ":")))
    path = C.FIGURES / out_name
    path.write_text(html)
    return str(path), path.stat().st_size


if __name__ == "__main__":
    print(build())
