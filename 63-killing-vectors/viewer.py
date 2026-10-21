"""The Kerr orbit explorer: one self-contained HTML page, no external libraries.

Pick a black hole's spin and an orbit's closest and farthest distances and its
tilt; the page finds the three constants of motion that make that orbit,

* E   the energy, from the time-translation Killing vector,
* L_z the angular momentum about the spin axis, from the rotation Killing vector,
* Q   Carter's constant, from the Killing tensor (the hidden symmetry),

and then draws the orbit. With all three fixed, the motion separates (in
"Mino time" lambda, d tau = Sigma d lambda) into two independent oscillations,

    (dr/dlambda)^2  = R(r)  = [E (r^2 + a^2) - a L]^2 - Delta [r^2 + (L - a E)^2 + Q]
    (dz/dlambda)^2  = Z(z)  = Q (1 - z^2) - z^2 [a^2 (1 - E^2)(1 - z^2) + L^2],   z = cos(theta),

integrated in the page with a fourth-order Runge-Kutta step on the smooth
second-order forms r'' = R'(r)/2, z'' = Z'(z)/2 (no trouble at turning points).
The page shows R(r) and Z(z): the orbit lives where they are positive. Units
G = c = M = 1. Everything is computed in the browser; the page follows the
viewer's light or dark setting.
"""

from __future__ import annotations

import json

import config as C

PRESETS = [
    {"name": "Tilted orbit, fast spin", "a": 0.9, "rp": 4.0, "ra": 12.0, "inc": 45},
    {"name": "Zoom-whirl", "a": 0.9, "rp": 2.25, "ra": 30.0, "inc": 20},
    {"name": "Polar orbit", "a": 0.9, "rp": 6.0, "ra": 10.0, "inc": 90},
    {"name": "Retrograde", "a": 0.9, "rp": 9.5, "ra": 14.0, "inc": 150},
    {"name": "No spin (Schwarzschild)", "a": 0.0, "rp": 8.0, "ra": 20.0, "inc": 30},
    {"name": "Hugging the horizon", "a": 0.99, "rp": 1.6, "ra": 3.0, "inc": 8},
]

HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Kerr Orbit Explorer</title>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --wash: #ecebe6; --horizon: #0b0b0b;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --s1: #3987e5; --s2: #d95926; --s3: #199e70; --wash: #262624; --horizon: #000000;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --wash: #262624; --horizon: #000000;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
       font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 18px 16px 40px; }
h1 { font-size: 20px; margin: 0 0 4px; }
.lede { color: var(--ink2); margin: 0 0 14px; max-width: 900px; }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 12px; }
.grid { display: grid; grid-template-columns: minmax(0, 1.25fr) minmax(0, 1fr); gap: 12px; }
@media (max-width: 860px) { .grid { grid-template-columns: minmax(0, 1fr); } }
.controls { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 10px 16px; margin-bottom: 12px; }
@media (max-width: 860px) { .controls { grid-template-columns: repeat(2, minmax(0, 1fr)); } }
label { display: block; font-size: 12px; color: var(--ink2); }
label b { color: var(--ink); font-variant-numeric: tabular-nums; }
input[type=range] { width: 100%; accent-color: var(--s1); }
.presets { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 12px; }
button { font: inherit; font-size: 12.5px; color: var(--ink); background: var(--surface);
         border: 1px solid var(--border); border-radius: 999px; padding: 4px 11px; cursor: pointer; }
button:hover { background: var(--wash); }
button.on { border-color: var(--s1); }
canvas { width: 100%; height: auto; display: block; touch-action: none; cursor: grab; }
svg { width: 100%; height: auto; display: block; }
h2 { font-size: 13px; margin: 0 0 6px; color: var(--ink2); font-weight: 600; }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
td { padding: 3px 4px; border-bottom: 1px solid var(--grid); }
td.v { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.note { color: var(--muted); font-size: 12px; margin-top: 6px; }
.warn { color: var(--s2); font-size: 13px; min-height: 1.2em; }
.legend { display: flex; flex-wrap: wrap; gap: 4px 14px; font-size: 12px; color: var(--ink2); margin-top: 6px; }
.sw { display: inline-block; width: 18px; height: 3px; vertical-align: middle; margin-right: 5px; border-radius: 2px; }
.prose { color: var(--ink2); max-width: 900px; margin-top: 14px; }
.prose p { margin: 6px 0; }
</style>
</head>
<body>
<main>
<h1>Kerr orbit explorer</h1>
<p class="lede">An orbit around a spinning black hole is fixed by three conserved numbers: its energy <i>E</i> and its
angular momentum about the spin axis <i>L</i>, which come from the spacetime's two symmetries (Killing vectors),
and Carter's constant <i>Q</i>, which comes from a hidden symmetry (a Killing tensor). Choose the orbit's closest
and farthest distances and its tilt; the page solves for <i>E</i>, <i>L</i> and <i>Q</i> and lets the orbit run.
Drag the picture to turn it.</p>
<div class="card">
  <div class="presets" id="presets"></div>
  <div class="controls">
    <label>spin a = cJ/GM² <b id="la"></b><input id="a" type="range" min="0" max="0.998" step="0.001"></label>
    <label>closest distance r<sub>p</sub> <b id="lrp"></b><input id="rp" type="range" min="1.2" max="40" step="0.01"></label>
    <label>farthest distance r<sub>a</sub> <b id="lra"></b><input id="ra" type="range" min="1.2" max="60" step="0.01"></label>
    <label>tilt ι (cos ι = L/√(L²+Q)) <b id="linc"></b><input id="inc" type="range" min="0" max="180" step="1"></label>
  </div>
  <div class="warn" id="warn" role="status"></div>
</div>
<div class="grid" style="margin-top:12px">
  <div class="card">
    <h2>The orbit (distances in GM/c²)</h2>
    <canvas id="view" width="1200" height="900" aria-label="three-dimensional orbit view"></canvas>
    <div class="legend">
      <span><span class="sw" style="background:var(--s1)"></span>orbit</span>
      <span><span class="sw" style="background:var(--horizon)"></span>horizon</span>
      <span><span class="sw" style="background:var(--s2)"></span>ergosphere edge</span>
      <span><span class="sw" style="background:var(--s3)"></span>innermost stable circular orbit</span>
    </div>
    <div style="display:flex; gap:8px; margin-top:8px; align-items:center">
      <button id="play">Pause</button>
      <label style="flex:1">animation speed<input id="speed" type="range" min="0.2" max="5" step="0.1" value="1"></label>
    </div>
  </div>
  <div>
    <div class="card">
      <h2>Where the orbit can go: R(r)/r⁴ ≥ 0 (radial) and Z(z) ≥ 0 (up-down; z = −1 south pole, +1 north)</h2>
      <svg id="pot" viewBox="0 0 560 230" role="img" aria-label="radial potential"></svg>
      <svg id="zpot" viewBox="0 0 560 120" role="img" aria-label="polar potential" style="margin-top:6px"></svg>
    </div>
    <div class="card" style="margin-top:12px">
      <h2>The constants of motion (per unit mass)</h2>
      <table id="tab"></table>
      <div class="note">r is the Boyer-Lindquist radius. ISCO and photon-orbit radii are for equatorial orbits at this spin.</div>
    </div>
  </div>
</div>
<div class="prose">
<p><b>What to look for.</b> Without spin (Schwarzschild) the orbit stays in one plane, and that plane's ellipse turns
(the precession Einstein found for Mercury). With spin the plane itself swings round the spin axis (frame dragging,
the Lense-Thirring effect), so a tilted orbit slowly fills a doughnut-shaped region. It never fills it chaotically:
the third constant <i>Q</i> keeps the radial and up-down oscillations independent. "Zoom-whirl" orbits swing in
from far away, whirl round several times near the closest point and zoom out again. Sliders that ask for an orbit
that cannot exist (too close in) are refused: there is no stable orbit inside the innermost stable one, and the
page says why.</p>
<p>Computed in your browser from the equations above (G = c = M = 1). Part of project N6 (Killing vectors) of the
astrophysics atlas.</p>
</div>
</main>
<script>
const PRESETS = __PRESETS__;
const $ = id => document.getElementById(id);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

// ---------- Kerr functions (G = c = M = 1) ----------
function Rfun(r, a, E, L, Q) {
  const D = r*r - 2*r + a*a, P = E*(r*r + a*a) - a*L;
  return P*P - D*(r*r + (L - a*E)**2 + Q);
}
function dRfun(r, a, E, L, Q) {
  const P = E*(r*r + a*a) - a*L;
  return 4*E*r*P - (2*r - 2)*(r*r + (L - a*E)**2 + Q) - (r*r - 2*r + a*a)*2*r;
}
function Zfun(z, a, E, L, Q) { return Q*(1 - z*z) - z*z*(a*a*(1 - E*E)*(1 - z*z) + L*L); }
function dZfun(z, a, E, L, Q) { return -2*z*(Q + L*L + a*a*(1 - E*E)*(1 - 2*z*z)); }
function iscoPro(a) {
  const z1 = 1 + Math.cbrt(1 - a*a)*(Math.cbrt(1 + a) + Math.cbrt(1 - a));
  const z2 = Math.sqrt(3*a*a + z1*z1);
  return [3 + z2 - Math.sqrt((3 - z1)*(3 + z1 + 2*z2)), 3 + z2 + Math.sqrt((3 - z1)*(3 + z1 + 2*z2))];
}
function photonOrbits(a) { return [2*(1 + Math.cos(2/3*Math.acos(-a))), 2*(1 + Math.cos(2/3*Math.acos(a)))]; }

// Solve for (E, L, Q) given turning points rp, ra and the inclination iota,
// cos(iota) = L / sqrt(L^2 + Q): unknowns (E, Lt) with L = Lt cos(iota), Q = Lt^2 sin^2(iota).
// Same algorithm as orbits.kerr_constants in the project's Python code.
function circularKerr(r, a) {   // BPT 1972; a < 0 means retrograde
  const s = a >= 0 ? 1 : -1, aa = Math.abs(a);
  const den = Math.pow(r, 0.75)*Math.sqrt(Math.pow(r, 1.5) - 3*Math.sqrt(r) + 2*s*aa);
  return [(Math.pow(r, 1.5) - 2*Math.sqrt(r) + s*aa)/den, s*(r*r - 2*s*aa*Math.sqrt(r) + aa*aa)/den];
}
function constants(a, rp, ra, inc) {
  const ci = Math.cos(inc*Math.PI/180), si = Math.sin(inc*Math.PI/180);
  const p = 2*rp*ra/(rp + ra), e = (ra - rp)/(ra + rp);
  let E, Lt;
  if (p - 3 - e*e > 0.5) {
    E = Math.sqrt((p - 2 - 2*e)*(p - 2 + 2*e)/(p*(p - 3 - e*e))); Lt = p/Math.sqrt(p - 3 - e*e);
  } else {
    const c = circularKerr(0.5*(rp + ra), inc <= 90 ? a : -a);
    if (!isFinite(c[0]) || !isFinite(c[1])) return null;
    E = c[0]; Lt = Math.abs(c[1]);
  }
  const F = (E, Lt) => { const L = Lt*ci, Q = (Lt*si)**2; return [Rfun(rp, a, E, L, Q)/rp**4, Rfun(ra, a, E, L, Q)/ra**4]; };
  for (let it = 0; it < 100; it++) {
    const f = F(E, Lt), h = 1e-7;
    const fE = F(E + h, Lt), fL = F(E, Lt + h);
    const J = [[(fE[0]-f[0])/h, (fL[0]-f[0])/h], [(fE[1]-f[1])/h, (fL[1]-f[1])/h]];
    const det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    if (!isFinite(det) || Math.abs(det) < 1e-300) return null;
    const dE = ( J[1][1]*f[0] - J[0][1]*f[1])/det, dL = (-J[1][0]*f[0] + J[0][0]*f[1])/det;
    E -= dE; Lt -= dL;
    if (Math.abs(dE) < 1e-14 && Math.abs(dL) < 1e-13) break;
  }
  const L = Lt*ci, Q = (Lt*si)**2, f = F(E, Lt);
  if (!(E > 0 && E < 1 && Lt > 0 && Math.abs(f[0]) < 1e-9 && Math.abs(f[1]) < 1e-9)) return null;
  for (let k = 1; k < 40; k++) { const r = rp + (ra - rp)*k/40; if (Rfun(r, a, E, L, Q) <= 0) return null; }
  return {E, L, Q};
}

// Integrate in Mino time with RK4; returns arrays of t, x, y, z (Boyer-Lindquist to Cartesian), r, z
function integrate(a, c, rp, ra, nSteps) {
  const {E, L, Q} = c;
  let r = rp, vr = 0, z = 0, vz = Math.sqrt(Math.max(Zfun(0, a, E, L, Q), 0)), t = 0, ph = 0, tau = 0;
  const out = {t: [], x: [], y: [], zc: [], r: [], z: []};
  const dl = Math.min(0.02, 0.4/Math.sqrt(Math.abs(dRfun(rp, a, E, L, Q)) + Math.abs(dZfun(0.5, a, E, L, Q)) + 1)) ;
  const rates = (r, z) => {
    const D = r*r - 2*r + a*a, P = E*(r*r + a*a) - a*L, s2 = Math.max(1 - z*z, 1e-12);
    return [(r*r + a*a)/D*P - a*(a*E*s2 - L), a*P/D - a*E + L/s2, r*r + a*a*z*z];
  };
  for (let i = 0; i < nSteps; i++) {
    const f = (r, vr, z, vz) => [vr, dRfun(r, a, E, L, Q)/2, vz, dZfun(z, a, E, L, Q)/2];
    const k1 = f(r, vr, z, vz);
    const k2 = f(r + dl/2*k1[0], vr + dl/2*k1[1], z + dl/2*k1[2], vz + dl/2*k1[3]);
    const k3 = f(r + dl/2*k2[0], vr + dl/2*k2[1], z + dl/2*k2[2], vz + dl/2*k2[3]);
    const k4 = f(r + dl*k3[0], vr + dl*k3[1], z + dl*k3[2], vz + dl*k3[3]);
    const g1 = rates(r, z);
    r += dl/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    vr += dl/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    z += dl/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
    vz += dl/6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
    const g2 = rates(r, z);
    t += dl*(g1[0] + g2[0])/2; ph += dl*(g1[1] + g2[1])/2; tau += dl*(g1[2] + g2[2])/2;
    if (i % 3 === 0) {
      const s = Math.sqrt(Math.max(1 - z*z, 0)), rr = Math.sqrt(r*r + a*a);
      out.t.push(t); out.x.push(rr*s*Math.cos(ph)); out.y.push(rr*s*Math.sin(ph)); out.zc.push(r*z);
      out.r.push(r); out.z.push(z);
    }
    if (!isFinite(r) || r < 1 + Math.sqrt(1 - a*a)) break;
  }
  out.tau = tau;
  return out;
}

// ---------- drawing ----------
let state = {a: 0.9, rp: 4, ra: 12, inc: 45}, orbit = null, consts = null, idx = 0, playing = true;
let yaw = -0.6, pitch = 0.42, drag = null, speed = 1;
const cv = $("view"), ctx = cv.getContext("2d");

function project(x, y, z, scale) {
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const x1 = cy*x - sy*y, y1 = sy*x + cy*y;
  const yy = cp*y1 - sp*z, zz = sp*y1 + cp*z;
  return [cv.width/2 + x1*scale, cv.height/2 - zz*scale, yy];
}
function ring(rad, zc, scale, color, dash, width) {
  ctx.beginPath();
  for (let k = 0; k <= 160; k++) {
    const p = 2*Math.PI*k/160, q = project(rad*Math.cos(p), rad*Math.sin(p), zc, scale);
    k ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]);
  }
  ctx.setLineDash(dash || []); ctx.strokeStyle = color; ctx.lineWidth = width || 2; ctx.stroke(); ctx.setLineDash([]);
}
function meridian(fr, scale, color, width, fill) {
  // a curve r(theta) in the plane facing the viewer (phi = yaw-based), drawn as a closed outline
  const phi0 = -yaw;
  ctx.beginPath();
  for (let k = 0; k <= 200; k++) {
    const th = 2*Math.PI*k/200, r = fr(th), s = Math.sin(th);
    const R = Math.sqrt(r*r + state.a*state.a)*s;
    const q = project(R*Math.cos(phi0), R*Math.sin(phi0), r*Math.cos(th), scale);
    k ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]);
  }
  if (fill) { ctx.fillStyle = color; ctx.fill(); } else { ctx.strokeStyle = color; ctx.lineWidth = width; ctx.stroke(); }
}
function draw() {
  const W = cv.width, H = cv.height, a = state.a;
  ctx.fillStyle = css("--surface"); ctx.fillRect(0, 0, W, H);
  const extent = Math.max(state.ra*1.15, 6);
  const scale = 0.46*Math.min(W, H)/extent;
  for (let rr of [5, 10, 20, 40, 60]) if (rr < extent*1.1) ring(rr, 0, scale, css("--grid"), [], 1.5);
  // spin axis
  let p0 = project(0, 0, -extent*0.75, scale), p1 = project(0, 0, extent*0.75, scale);
  ctx.strokeStyle = css("--axis"); ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(p0[0], p0[1]); ctx.lineTo(p1[0], p1[1]); ctx.stroke();
  const rplus = 1 + Math.sqrt(1 - a*a);
  if (a > 0) meridian(th => 1 + Math.sqrt(1 - a*a*Math.cos(th)**2), scale, css("--s2"), 2, false);
  meridian(() => rplus, scale, css("--horizon"), 0, true);
  meridian(() => rplus, scale, css("--ink2"), 1.5, false);
  ring(iscoPro(a)[state.inc <= 90 ? 0 : 1], 0, scale, css("--s3"), [10, 8], 2);
  if (orbit) {
    const n = orbit.x.length, upto = Math.min(idx, n - 1), recent = Math.max(0, upto - 700);
    ctx.strokeStyle = css("--s1");
    for (const [from, alpha, width] of [[0, 0.18, 1.5], [recent, 0.95, 2.4]]) {
      ctx.globalAlpha = alpha; ctx.lineWidth = width; ctx.beginPath();
      for (let i = from; i <= upto; i++) {
        const q = project(orbit.x[i], orbit.y[i], orbit.zc[i], scale);
        i > from ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]);
      }
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
    const q = project(orbit.x[upto], orbit.y[upto], orbit.zc[upto], scale);
    ctx.beginPath(); ctx.arc(q[0], q[1], 9, 0, 2*Math.PI); ctx.fillStyle = css("--s1"); ctx.fill();
    ctx.lineWidth = 3; ctx.strokeStyle = css("--surface"); ctx.stroke();
  }
  ctx.fillStyle = css("--muted"); ctx.font = "26px system-ui, sans-serif";
  ctx.fillText("rings at r = 5, 10, 20, 40", 24, H - 24);
}
function svgPlot(el, xs, ys, xr, yr, extra) {
  const W = 560, H = el.viewBox.baseVal.height, l = 46, b = 26, t = 8, rr = 10;
  const X = v => l + (v - xr[0])/(xr[1] - xr[0])*(W - l - rr), Y = v => t + (yr[1] - v)/(yr[1] - yr[0])*(H - t - b);
  let s = `<rect x="0" y="0" width="${W}" height="${H}" fill="var(--surface)"/>`;
  s += `<line x1="${l}" x2="${W-rr}" y1="${Y(0)}" y2="${Y(0)}" stroke="var(--axis)" stroke-width="1"/>`;
  for (const tk of extra.xticks) s += `<text x="${X(tk)}" y="${H-8}" fill="var(--muted)" font-size="11" text-anchor="middle">${tk}</text>`;
  s += `<text x="${l-6}" y="${Y(yr[1])+10}" fill="var(--muted)" font-size="11" text-anchor="end">${extra.ylab}</text>`;
  let path = "", fillp = "";
  for (let i = 0; i < xs.length; i++) { path += (i ? "L" : "M") + X(xs[i]).toFixed(1) + "," + Y(Math.max(Math.min(ys[i], yr[1]), yr[0])).toFixed(1); }
  for (let i = 0; i < xs.length; i++) if (ys[i] > 0) fillp += `<line x1="${X(xs[i]).toFixed(1)}" x2="${X(xs[i]).toFixed(1)}" y1="${Y(0)}" y2="${Y(Math.min(ys[i], yr[1])).toFixed(1)}" stroke="var(--s1)" stroke-opacity="0.13" stroke-width="2"/>`;
  s += fillp + `<path d="${path}" fill="none" stroke="var(--s1)" stroke-width="2"/>`;
  for (const m of extra.marks || []) s += `<line x1="${X(m.x)}" x2="${X(m.x)}" y1="${t}" y2="${H-b}" stroke="${m.c}" stroke-width="1.5" ${m.d ? 'stroke-dasharray="4 3"' : ''}/><text x="${X(m.x)+3}" y="${t+11}" fill="var(--ink2)" font-size="10.5">${m.t}</text>`;
  if (extra.dot) s += `<circle cx="${X(extra.dot[0])}" cy="${Y(0)}" r="5" fill="var(--s1)" stroke="var(--surface)" stroke-width="2"/>`;
  s += `<text x="${W-rr}" y="${H-8}" fill="var(--muted)" font-size="11" text-anchor="end">${extra.xlab}</text>`;
  el.innerHTML = s;
}
function drawPotentials() {
  if (!consts) { $("pot").innerHTML = ""; $("zpot").innerHTML = ""; return; }
  const a = state.a, {E, L, Q} = consts, rplus = 1 + Math.sqrt(1 - a*a);
  const rmax = state.ra*1.3, xs = [], ys = [];
  for (let i = 0; i <= 400; i++) { const r = rplus + (rmax - rplus)*i/400; xs.push(r); ys.push(Rfun(r, a, E, L, Q)/r**4); }
  const ymax = Math.max(...ys.filter(v => v > 0), 1e-6)*1.25, ymin = -ymax*0.9;
  const ticks = []; const step = rmax > 40 ? 10 : rmax > 15 ? 5 : 2; for (let v = 0; v <= rmax; v += step) if (v >= rplus) ticks.push(v);
  const cur = orbit ? orbit.r[Math.min(idx, orbit.r.length - 1)] : null;
  svgPlot($("pot"), xs, ys, [rplus, rmax], [ymin, ymax], {xticks: ticks, xlab: "r (GM/c²)", ylab: "R/r⁴",
    marks: [{x: rplus, c: "var(--ink)", t: "horizon"}], dot: cur ? [cur] : null});
  const zs = [], zy = [];
  for (let i = 0; i <= 300; i++) { const z = -1 + 2*i/300; zs.push(z); zy.push(Zfun(z, a, E, L, Q)); }
  const zmax = Math.max(...zy.filter(v => v > 0), 1e-6)*1.25;
  const zc = orbit ? orbit.z[Math.min(idx, orbit.z.length - 1)] : null;
  svgPlot($("zpot"), zs, zy, [-1, 1], [-zmax*0.8, zmax], {xticks: [-1, -0.5, 0, 0.5, 1], xlab: "z = cos θ", ylab: "Z", dot: zc !== null ? [zc] : null});
}
function table() {
  const a = state.a, isco = iscoPro(a), ph = photonOrbits(a), rplus = 1 + Math.sqrt(1 - a*a);
  const rows = consts ? [
    ["energy E (units of mc²)", consts.E.toFixed(6)],
    ["binding energy 1 − E", (100*(1 - consts.E)).toFixed(3) + " %"],
    ["angular momentum about the spin axis L", consts.L.toFixed(4)],
    ["Carter's constant Q", consts.Q.toFixed(4)],
    ["Q + L² (total angular momentum², if a = 0)", (consts.Q + consts.L*consts.L).toFixed(4)],
  ] : [["no bound orbit with these settings", "—"]];
  rows.push(["horizon r₊", rplus.toFixed(3)], ["ISCO (prograde / retrograde)", isco[0].toFixed(3) + " / " + isco[1].toFixed(3)],
            ["photon orbit (prograde / retrograde)", ph[0].toFixed(3) + " / " + ph[1].toFixed(3)]);
  if (orbit) rows.push(["proper time shown (GM/c³)", orbit.tau.toFixed(0)]);
  $("tab").innerHTML = rows.map(r => `<tr><td>${r[0]}</td><td class="v">${r[1]}</td></tr>`).join("");
}
function recompute() {
  $("la").textContent = state.a.toFixed(3); $("lrp").textContent = state.rp.toFixed(2); $("lra").textContent = state.ra.toFixed(2);
  $("linc").textContent = state.inc + "°" + (state.inc > 90 ? " (retrograde)" : "");
  ["a", "rp", "ra", "inc"].forEach(k => { $(k).value = state[k]; });
  let msg = "";
  if (state.ra < state.rp) { msg = "The farthest distance must not be smaller than the closest."; consts = null; }
  else {
    consts = constants(state.a, state.rp, state.ra, state.inc);
    if (!consts) msg = "No bound orbit has these turning points: at this spin and tilt the closest approach is inside the region where orbits plunge (move r_p out, or r_a in).";
  }
  $("warn").textContent = msg;
  orbit = consts ? integrate(state.a, consts, state.rp, state.ra, 60000) : null;
  idx = 0; table(); drawPotentials(); draw();
}
function frame() {
  if (playing && orbit) {
    idx += Math.max(1, Math.round(6*speed));
    if (idx >= orbit.x.length) idx = 0;
    if (idx % 4 === 0) drawPotentials();
  }
  draw(); requestAnimationFrame(frame);
}
// ---------- wiring ----------
PRESETS.forEach((p, i) => {
  const b = document.createElement("button"); b.textContent = p.name;
  b.onclick = () => { state = {a: p.a, rp: p.rp, ra: p.ra, inc: p.inc}; recompute();
    document.querySelectorAll("#presets button").forEach(x => x.classList.remove("on")); b.classList.add("on"); };
  $("presets").appendChild(b);
});
["a", "rp", "ra", "inc"].forEach(k => $(k).addEventListener("input", e => { state[k] = parseFloat(e.target.value); recompute(); }));
$("play").onclick = () => { playing = !playing; $("play").textContent = playing ? "Pause" : "Play"; };
$("speed").addEventListener("input", e => { speed = parseFloat(e.target.value); });
cv.addEventListener("pointerdown", e => { drag = [e.clientX, e.clientY, yaw, pitch]; cv.setPointerCapture(e.pointerId); });
cv.addEventListener("pointermove", e => { if (!drag) return; yaw = drag[2] + (e.clientX - drag[0])*0.008;
  pitch = Math.max(-1.5, Math.min(1.5, drag[3] + (e.clientY - drag[1])*0.008)); });
cv.addEventListener("pointerup", () => { drag = null; });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { drawPotentials(); });
recompute(); document.querySelector("#presets button").classList.add("on"); requestAnimationFrame(frame);
</script>
</body>
</html>
"""


def make() -> str:
    html = HTML.replace("__PRESETS__", json.dumps(PRESETS))
    C.VIEWER_PATH.write_text(html)
    return str(C.VIEWER_PATH)


if __name__ == "__main__":
    print(make())
