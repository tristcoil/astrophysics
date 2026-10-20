"""The interactive orbit gallery: six famous binaries moving in real time, to one scale.

A self-contained HTML page (no external libraries; the orbital elements are
embedded and the Kepler equation is solved in the browser), so it works
offline as one file. Views: the true orbit face-on, or the orbit as seen from
Earth (north up, east left); motion of both stars about their centre of mass,
or the companion around the primary as visual observers measure it. It follows
the viewer's light/dark setting.
"""

from __future__ import annotations

import json

import config as C
import plotstyle as ps

THEMES = {
    "light": dict(surface=ps.SURFACE, ink=ps.INK, ink2=ps.INK_SECONDARY, muted=ps.MUTED, grid=ps.GRID,
                  a=ps.BLUE, b=ps.ORANGE, a_soft="#b7d3f6", b_soft="#f6c3ad", panel="#f5f4f0", border="#e1e0d9"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781", grid="#383835",
                 a="#3987e5", b="#d95926", a_soft="#1c4a80", b_soft="#7a3418", panel="#222220", border="#383835"),
}


def page_data(res) -> list:
    out = []
    for g in res["gallery"]:
        P, T, e, a, i, O, w = g["elements"]
        out.append({"name": g["name"], "wds": g["wds"], "disc": g["disc"].replace("  ", " "), "grade": g["grade"],
                    "ref": g["ref"], "P": P, "T": T, "e": e, "a_arcsec": a, "i": i, "Omega": O, "omega": w,
                    "a_au": g["a_au"], "fB": g["fB"], "MA": g["MA"], "MB": g["MB"], "Mtot": g["Mtot"],
                    "Mtot_err": g["Mtot_err"], "dist_pc": 1e3 / g["plx"], "plx_src": g["plx_src"],
                    "mass_src": g["mass_src"], "v1": g["v1"], "v2": g["v2"]})
    return out


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Binary Orbit Gallery</title>
<style>
  :root { color-scheme: light; --surface: __L_surface__; --ink: __L_ink__; --ink2: __L_ink2__; --muted: __L_muted__;
          --grid: __L_grid__; --a: __L_a__; --b: __L_b__; --a-soft: __L_a_soft__; --b-soft: __L_b_soft__;
          --panel: __L_panel__; --border: __L_border__; }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) { color-scheme: dark; --surface: __D_surface__; --ink: __D_ink__; --ink2: __D_ink2__;
          --muted: __D_muted__; --grid: __D_grid__; --a: __D_a__; --b: __D_b__; --a-soft: __D_a_soft__;
          --b-soft: __D_b_soft__; --panel: __D_panel__; --border: __D_border__; }
  }
  :root[data-theme="dark"] { color-scheme: dark; --surface: __D_surface__; --ink: __D_ink__; --ink2: __D_ink2__;
          --muted: __D_muted__; --grid: __D_grid__; --a: __D_a__; --b: __D_b__; --a-soft: __D_a_soft__;
          --b-soft: __D_b_soft__; --panel: __D_panel__; --border: __D_border__; }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
  main { max-width: 1180px; margin: 0 auto; padding: 20px 16px 40px; }
  h1 { font-size: 1.45rem; margin: 0 0 6px; }
  p.lede { margin: 0 0 14px; color: var(--ink2); max-width: 80ch; }
  .controls { display: flex; flex-wrap: wrap; gap: 10px 18px; align-items: center; padding: 10px 12px;
              background: var(--panel); border: 1px solid var(--border); border-radius: 8px; margin-bottom: 14px; }
  .seg { display: inline-flex; border: 1px solid var(--border); border-radius: 6px; overflow: hidden; }
  .seg button { background: transparent; color: var(--ink2); border: 0; padding: 6px 10px; font: inherit;
                font-size: 0.85rem; cursor: pointer; }
  .seg button[aria-pressed="true"] { background: var(--ink); color: var(--surface); }
  .ctl { display: inline-flex; align-items: center; gap: 8px; font-size: 0.85rem; color: var(--ink2); }
  button.play { background: var(--ink); color: var(--surface); border: 0; border-radius: 6px; padding: 6px 14px;
                font: inherit; font-size: 0.85rem; cursor: pointer; min-width: 74px; }
  #year { font-variant-numeric: tabular-nums; font-weight: 600; color: var(--ink); min-width: 5.5em; }
  input[type=range] { width: 130px; accent-color: var(--a); }
  .grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px; }
  @media (max-width: 900px) { .grid { grid-template-columns: repeat(2, minmax(0, 1fr)); } }
  @media (max-width: 560px) { .grid { grid-template-columns: 1fr; } }
  .card { border: 1px solid var(--border); border-radius: 8px; padding: 10px 10px 8px; background: var(--surface); }
  .card h2 { font-size: 1rem; margin: 0; }
  .card .sub { font-size: 0.78rem; color: var(--muted); margin: 0 0 4px; }
  canvas { width: 100%; height: auto; display: block; }
  .facts { font-size: 0.8rem; color: var(--ink2); display: grid; grid-template-columns: auto 1fr; gap: 1px 10px;
           margin-top: 4px; }
  .facts span:nth-child(odd) { color: var(--muted); }
  .legend { display: flex; gap: 16px; font-size: 0.8rem; color: var(--ink2); margin: 0 0 10px; flex-wrap: wrap; }
  .sw { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; vertical-align: -1px; }
  footer { margin-top: 18px; font-size: 0.78rem; color: var(--muted); max-width: 95ch; }
  .theme { margin-left: auto; }
</style>
</head>
<body>
<main>
  <h1>Six famous binaries, to one scale</h1>
  <p class="lede">Real orbits from the Sixth Catalog of Orbits of Visual Binary Stars, turned into astronomical units
  with each system's measured distance. Kepler's third law runs the clock: the wider the orbit and the lighter the
  stars, the slower the dance. The grey rings are the orbits of Jupiter, Saturn, Uranus and Neptune.</p>
  <div class="controls" role="group" aria-label="Animation controls">
    <button class="play" id="play" aria-label="Pause">Pause</button>
    <span class="ctl">Year <span id="year">2026</span></span>
    <label class="ctl">Speed <input type="range" id="speed" min="1" max="60" step="1" value="10">
      <span id="speedv">10 yr/s</span></label>
    <span class="seg" role="group" aria-label="View">
      <button data-k="view" data-v="face" aria-pressed="true">True shape</button>
      <button data-k="view" data-v="sky" aria-pressed="false">As seen from Earth</button>
    </span>
    <span class="seg" role="group" aria-label="Motion">
      <button data-k="frame" data-v="bary" aria-pressed="true">Both stars</button>
      <button data-k="frame" data-v="rel" aria-pressed="false">B around A</button>
    </span>
    <span class="seg" role="group" aria-label="Scale">
      <button data-k="scale" data-v="same" aria-pressed="true">One scale</button>
      <button data-k="scale" data-v="fit" aria-pressed="false">Fit each</button>
    </span>
    <span class="seg theme" role="group" aria-label="Theme">
      <button data-k="theme" data-v="auto" aria-pressed="true">Auto</button>
      <button data-k="theme" data-v="light" aria-pressed="false">Light</button>
      <button data-k="theme" data-v="dark" aria-pressed="false">Dark</button>
    </span>
  </div>
  <div class="legend"><span><i class="sw" style="background:var(--a)"></i>primary (A)</span>
    <span><i class="sw" style="background:var(--b)"></i>companion (B)</span>
    <span>+ centre of mass</span><span id="scalelabel"></span></div>
  <div class="grid" id="grid"></div>
  <footer>Orbital elements: ORB6 (Hartkopf, Mason &amp; Worley; GSU/USNO), October 2026; Sirius from this project's
  own fit to Bond et al.'s 1862-2016 measures. Distances and mass ratios as listed under each panel. "As seen from
  Earth" shows the projected orbit with north up and east to the left; the ascending node is ambiguous without radial
  velocities, so the sense of the tilt along the line of sight is not shown. Generated by
  50-binary-orbits/interactive.py.</footer>
</main>
<script>
const DATA = __DATA__;
const TAU = 2 * Math.PI, DEG = Math.PI / 180;
const state = { view: "face", frame: "bary", scale: "same", playing: true, speed: 10, year: 2026.0 };
const PLANETS = [5.20, 9.54, 19.19, 30.07];

function kepler(M, e) {
  M = ((M % TAU) + TAU) % TAU;
  let E = e < 0.8 ? M : Math.PI;
  for (let k = 0; k < 50; k++) {
    const d = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
    E -= d;
    if (Math.abs(d) < 1e-12) break;
  }
  return E;
}
// position of B relative to A (au): face-on (orbital plane, periastron along +x) or sky (east, north)
function relpos(s, t) {
  const E = kepler(TAU * (t - s.T) / s.P, s.e);
  const X = s.a_au * (Math.cos(E) - s.e), Y = s.a_au * Math.sqrt(1 - s.e * s.e) * Math.sin(E);
  if (state.view === "face") return [X, Y];
  const w = s.omega * DEG, O = s.Omega * DEG, i = s.i * DEG;
  const A = Math.cos(w) * Math.cos(O) - Math.sin(w) * Math.sin(O) * Math.cos(i);
  const B = Math.cos(w) * Math.sin(O) + Math.sin(w) * Math.cos(O) * Math.cos(i);
  const F = -Math.sin(w) * Math.cos(O) - Math.cos(w) * Math.sin(O) * Math.cos(i);
  const G = -Math.sin(w) * Math.sin(O) + Math.cos(w) * Math.cos(O) * Math.cos(i);
  const north = A * X + F * Y, east = B * X + G * Y;
  return [-east, north];   // east to the left on the screen
}
function positions(s, t) {
  const [x, y] = relpos(s, t);
  if (state.frame === "rel") return { A: [0, 0], B: [x, y] };
  return { A: [-s.fB * x, -s.fB * y], B: [(1 - s.fB) * x, (1 - s.fB) * y] };
}
function extent(s) {
  let m = 0;
  for (let k = 0; k < 360; k++) {
    const p = positions(s, s.T + s.P * k / 360);
    m = Math.max(m, Math.hypot(...p.A), Math.hypot(...p.B));
  }
  return m;
}
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const fmt = (x, d) => x.toLocaleString("en-US", { minimumFractionDigits: d, maximumFractionDigits: d });

const grid = document.getElementById("grid");
const cells = DATA.map((s) => {
  const card = document.createElement("div");
  card.className = "card";
  card.innerHTML = `<h2>${s.name}</h2><p class="sub">${s.wds} ${s.disc} · ORB6 grade ${s.grade} (${s.ref})</p>`;
  const cv = document.createElement("canvas");
  cv.setAttribute("role", "img");
  cv.setAttribute("aria-label", `${s.name}: orbit animation, period ${fmt(s.P, 1)} years`);
  card.appendChild(cv);
  const f = document.createElement("div");
  f.className = "facts";
  f.innerHTML = `<span>period</span><span>${fmt(s.P, 2)} yr</span>
    <span>semi-major axis</span><span>${fmt(s.a_arcsec, 3)}" = ${fmt(s.a_au, 1)} au</span>
    <span>eccentricity</span><span>${fmt(s.e, 3)}</span>
    <span>distance</span><span>${fmt(s.dist_pc, 2)} pc · ${s.plx_src}</span>
    <span>total mass</span><span>${fmt(s.Mtot, s.Mtot_err < 0.01 ? 3 : 2)} ± ${fmt(s.Mtot_err, s.Mtot_err < 0.01 ? 3 : 2)} M☉</span>
    <span>A + B</span><span>${fmt(s.MA, 2)} + ${fmt(s.MB, 2)} M☉ · ${s.mass_src}</span>`;
  card.appendChild(f);
  grid.appendChild(card);
  return { s, cv, trailA: [], trailB: [] };
});

function sizeCanvas(cv) {
  const r = window.devicePixelRatio || 1, w = cv.clientWidth || 300, h = Math.round(w * 0.78);
  if (cv.width !== Math.round(w * r) || cv.height !== Math.round(h * r)) {
    cv.width = Math.round(w * r); cv.height = Math.round(h * r); cv.style.height = h + "px";
  }
  return [w, h, r];
}
let SAME = 1;
function recomputeScale() {
  SAME = Math.max(...cells.map((c) => extent(c.s)));
  document.getElementById("scalelabel").textContent =
    state.scale === "same" ? `all panels: ${fmt(SAME, 0)} au from centre to edge` : "each panel scaled to fit";
}
function draw() {
  const ink = css("--ink"), muted = css("--muted"), gridc = css("--grid"), ca = css("--a"), cb = css("--b"),
        csa = css("--a-soft"), csb = css("--b-soft");
  for (const c of cells) {
    const [w, h, r] = sizeCanvas(c.cv);
    const ctx = c.cv.getContext("2d");
    ctx.setTransform(r, 0, 0, r, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const R = state.scale === "same" ? SAME : extent(c.s);
    const k = 0.46 * Math.min(w, h) / R, cx = w / 2, cy = h / 2;
    const X = (p) => cx + p[0] * k, Y = (p) => cy - p[1] * k;
    ctx.lineWidth = 1; ctx.strokeStyle = gridc;
    for (const pr of PLANETS) { ctx.beginPath(); ctx.arc(cx, cy, pr * k, 0, TAU); ctx.stroke(); }
    // full orbits
    for (const [who, col] of [["A", csa], ["B", csb]]) {
      ctx.beginPath(); ctx.strokeStyle = col; ctx.lineWidth = 1.2;
      for (let j = 0; j <= 360; j++) {
        const p = positions(c.s, c.s.T + c.s.P * j / 360)[who];
        j ? ctx.lineTo(X(p), Y(p)) : ctx.moveTo(X(p), Y(p));
      }
      ctx.stroke();
    }
    // trail over the last few percent of the period
    for (const [who, col] of [["A", ca], ["B", cb]]) {
      ctx.beginPath(); ctx.strokeStyle = col; ctx.lineWidth = 2;
      const span = Math.min(0.12 * c.s.P, 12);
      for (let j = 0; j <= 30; j++) {
        const p = positions(c.s, state.year - span * (1 - j / 30))[who];
        j ? ctx.lineTo(X(p), Y(p)) : ctx.moveTo(X(p), Y(p));
      }
      ctx.stroke();
    }
    ctx.strokeStyle = muted; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(cx - 5, cy); ctx.lineTo(cx + 5, cy); ctx.moveTo(cx, cy - 5); ctx.lineTo(cx, cy + 5); ctx.stroke();
    const pos = positions(c.s, state.year);
    const rad = (m) => 3 + 3.2 * Math.cbrt(Math.max(m, 0.05));
    for (const [who, col, m] of [["A", ca, c.s.MA], ["B", cb, c.s.MB]]) {
      ctx.beginPath(); ctx.fillStyle = col; ctx.strokeStyle = css("--surface"); ctx.lineWidth = 1.5;
      ctx.arc(X(pos[who]), Y(pos[who]), rad(m), 0, TAU); ctx.fill(); ctx.stroke();
    }
    if (state.view === "sky") {
      ctx.fillStyle = muted; ctx.font = "11px system-ui, sans-serif";
      ctx.fillText("N ↑   E ←", 8, 16);
    }
    // separation now
    const d = Math.hypot(pos.B[0] - pos.A[0], pos.B[1] - pos.A[1]);
    ctx.fillStyle = muted; ctx.font = "11px system-ui, sans-serif"; ctx.textAlign = "right";
    ctx.fillText(`${state.view === "sky" ? "projected " : ""}separation ${fmt(d, 1)} au`, w - 8, h - 8);
    ctx.textAlign = "left";
  }
  document.getElementById("year").textContent = String(Math.floor(state.year));
}
let last = null;
function tick(ts) {
  if (last !== null && state.playing) state.year += state.speed * (ts - last) / 1000;
  last = ts;
  draw();
  requestAnimationFrame(tick);
}
document.querySelectorAll(".seg button").forEach((b) => b.addEventListener("click", () => {
  const key = b.dataset.k, v = b.dataset.v;
  b.parentElement.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
  if (key === "theme") {
    if (v === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", v);
  } else { state[key] = v; recomputeScale(); }
  draw();
}));
const play = document.getElementById("play");
play.addEventListener("click", () => {
  state.playing = !state.playing;
  play.textContent = state.playing ? "Pause" : "Play";
  play.setAttribute("aria-label", state.playing ? "Pause" : "Play");
});
const sp = document.getElementById("speed");
sp.addEventListener("input", () => { state.speed = +sp.value; document.getElementById("speedv").textContent = `${sp.value} yr/s`; });
if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) { state.playing = false; play.textContent = "Play"; }
window.addEventListener("resize", draw);
recomputeScale();
requestAnimationFrame(tick);
</script>
</body>
</html>
"""


def write_page(res, path=C.FIGURE_DIR / "orbit_gallery.html"):
    html = TEMPLATE.replace("__DATA__", json.dumps(page_data(res), ensure_ascii=False))
    for mode, key in (("light", "L"), ("dark", "D")):
        for k, v in THEMES[mode].items():
            html = html.replace(f"__{key}_{k}__", v)
    path.write_text(html, encoding="utf-8")
    return path
