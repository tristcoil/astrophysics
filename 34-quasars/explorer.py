"""An offline page: figures/quasar_time_machine.html.

Drag a redshift slider and watch, all from this project's own numbers:

* **The spectrum through the filters.** Our composite quasar spectrum moved to
  that redshift (with the mean intergalactic absorption of Madau 1995, the same
  model as figure 4), over the five SDSS filters, emission lines labelled.
* **Quasars among the stars.** Where that quasar lands in the u-g / g-r plane
  relative to real stars, and how complete the SDSS colour selection was there.
* **The cosmic clock.** The age of the universe in D3's cosmology, the density
  of luminous quasars then (our DR7 measurement), the black hole masses DR16Q
  sees at that redshift, and how many Eddington e-foldings a seed born at z = 30
  could have had.

Self-contained: data embedded as JSON, drawn as SVG, no network or libraries,
light and dark themes.
"""

from __future__ import annotations

import json

import numpy as np

import physics
from config import FIGURES, LUMINOUS_M_CUT
from photometry import load_filter

LINES = [(1033, "Lyβ"), (1216, "Lyα"), (1400, "Si IV"), (1549, "C IV"), (1909, "C III]"),
         (2798, "Mg II"), (3727, "[O II]"), (4861, "Hβ"), (5007, "[O III]"), (6563, "Hα")]


def _r(a, nd=4):
    return [None if not np.isfinite(v) else round(float(v), nd) for v in np.asarray(a, float)]


def build_payload(res: dict) -> dict:
    comp, ph, q, m, g = res["composite"], res["photometry"], res["qlf"], res["masses"], res["growth"]
    w = comp["wave"]
    keep = slice(0, None, 3)
    intr = ph["intrinsic"]
    filt = {}
    for b in "ugriz":
        fw, fr = load_filter(b)
        filt[b] = {"w": _r(fw, 0), "r": _r(fr / fr.max(), 3)}
    tr = ph["track"]
    rng = np.random.default_rng(3)
    st = ph["star_colours"]
    pick = rng.choice(len(st), size=min(2500, len(st)), replace=False)
    obs = ph["observed"]
    era = q["era"][LUMINOUS_M_CUT]
    pk = q["era_peak"]
    b = pk["bins"]
    x = np.log10(1 + b["z_mean"].to_numpy())
    y = np.log10(b["density"].to_numpy() * 1e9)   # per Gpc^3
    sig = b["err"].to_numpy() / b["density"].to_numpy() / np.log(10)
    c = np.polyfit(x, y, 2, w=1 / sig)
    zage = np.round(np.arange(0.0, 31.0, 0.02), 2)
    pz = m["per_z"]
    hz = g["hz"]
    return {
        "comp": {"w": _r(w[keep], 0), "f": _r(intr[keep] / np.nanmax(intr), 5)},
        "filters": filt,
        "lines": LINES,
        "track": {"z": _r(ph["zgrid"], 2), "ug": _r(tr["u-g"], 3), "gr": _r(tr["g-r"], 3),
                  "ri": _r(tr["r-i"], 3), "iz": _r(tr["i-z"], 3),
                  "dist": _r(ph["locus_distance_ugri"], 3), "distgriz": _r(ph["locus_distance_griz"], 3),
                  "sel": _r(ph["selfunc_marginal"], 3)},
        "stars": [[round(float(st[k, 1]), 3), round(float(st[k, 0]), 3)] for k in pick],
        "obs": {"z": _r(obs["u-g"]["z"], 3), "ug": _r(obs["u-g"]["median"], 3), "gr": _r(obs["g-r"]["median"], 3)},
        "era": {"z": _r(era["z_mean"], 3), "n": _r(era["density"] * 1e9, 3), "e": _r(era["err"] * 1e9, 3),
                "complete": [bool(v) for v in era["complete"]],
                "fit": _r(c, 6), "zmin": float(b["z_mean"].min()), "zmax": float(b["z_mean"].max()),
                "zpeak": round(pk["z_peak"], 3), "agepeak": round(pk["age_at_peak_gyr"], 3)},
        "ages": {"z": _r(zage, 2), "t": _r(physics.age_myr(zage), 1)},
        "mass": {"z": _r(0.5 * (pz["z_lo"] + pz["z_hi"]), 3), "p50": _r(pz["m_p50"], 2), "p95": _r(pz["m_p95"], 2),
                 "edd": _r(pz["edd_p50"], 2)},
        "highz": {"z": _r(hz["z"], 3), "m": _r(hz["log_mbh"], 2), "name": hz["name"].tolist()},
        "efold": round(g["efold_myr"], 2),
        "t30": round(float(physics.age_myr(30.0)), 1),
        "cut": LUMINOUS_M_CUT,
    }


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Quasar time machine</title>
<style>
:root{--bg:#f9f9f7;--panel:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;
--blue:#2a78d6;--orange:#eb6834;--violet:#4a3aa7;--star:#a9a79f;--filt:#ecebe5;--filtline:#c3c2b7;--warn:#d03b3b}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#0d0d0d;--panel:#1a1a19;--ink:#ffffff;--ink2:#c3c2b7;
--muted:#898781;--grid:#2c2c2a;--axis:#383835;--blue:#3987e5;--orange:#d95926;--violet:#9085e9;--star:#6f6d66;--filt:#262624;--filtline:#4a4944;--warn:#e66767}}
:root[data-theme="dark"]{--bg:#0d0d0d;--panel:#1a1a19;--ink:#ffffff;--ink2:#c3c2b7;--muted:#898781;--grid:#2c2c2a;--axis:#383835;
--blue:#3987e5;--orange:#d95926;--violet:#9085e9;--star:#6f6d66;--filt:#262624;--filtline:#4a4944;--warn:#e66767}
html,body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1180px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:6px 0 2px}p.lede{color:var(--ink2);margin:0 0 12px;max-width:900px}
.card{background:var(--panel);border-radius:10px;box-shadow:0 0 0 1px rgba(11,11,11,.08);padding:12px 14px;margin-bottom:14px}
:root[data-theme="dark"] .card{box-shadow:0 0 0 1px rgba(255,255,255,.10)}
.row{display:flex;flex-wrap:wrap;gap:14px}.row>.card{flex:1 1 360px;min-width:0}
h2{font-size:15px;margin:0 0 4px}.sub{color:var(--muted);font-size:12px;margin:0 0 6px}
svg{width:100%;height:auto;display:block}
.ctrl{display:flex;align-items:center;gap:12px;flex-wrap:wrap}
input[type=range]{flex:1 1 320px;accent-color:var(--orange)}
.zval{font-size:28px;font-weight:600;font-variant-numeric:tabular-nums;min-width:4.5em}
.readout{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:6px 18px;margin-top:8px}
.readout div{color:var(--ink2);font-size:13px}.readout b{color:var(--ink);font-weight:600;font-variant-numeric:tabular-nums}
.btns button{font:inherit;font-size:12px;background:none;border:1px solid var(--axis);color:var(--ink2);border-radius:6px;padding:3px 8px;cursor:pointer}
.btns button:focus-visible,input:focus-visible{outline:2px solid var(--blue);outline-offset:2px}
.note{color:var(--muted);font-size:12px;margin-top:6px}
text{fill:var(--ink2);font-size:11px;font-family:system-ui,-apple-system,"Segoe UI",sans-serif}
.tick{fill:var(--muted);font-size:10px}.lab{fill:var(--ink);font-size:11px}
.warn{color:var(--warn);font-weight:600}
</style></head><body><main>
<h1>Quasar time machine</h1>
<p class="lede">Move a typical quasar through cosmic time. The spectrum is our own average of 4,446 SDSS quasar spectra (project I3);
everything else on this page is measured in the same project. Drag the slider, or use the arrow keys.</p>
<div class="card"><div class="ctrl"><label for="z">Redshift</label>
<input id="z" type="range" min="0.10" max="7.64" step="0.02" value="2.70" aria-describedby="zread">
<span class="zval" id="zread">z = 2.70</span>
<span class="btns"><button data-z="0.3">0.3</button> <button data-z="1.0">1</button> <button data-z="2.0">2</button>
<button data-z="2.7">2.7</button> <button data-z="3.5">3.5</button> <button data-z="4.5">4.5</button>
<button data-z="6.2">6.2</button> <button data-z="7.64">7.64</button>
<button id="theme" aria-label="switch light/dark">◐</button></span></div>
<div class="readout" id="read"></div></div>
<div class="card"><h2>The light we receive</h2><p class="sub">Our composite spectrum shifted to this redshift (log flux), with the mean
absorption by intergalactic hydrogen (Madau 1995 model) and the five SDSS filters underneath.</p>
<svg id="spec" viewBox="0 0 1000 300" role="img" aria-label="redshifted quasar spectrum over the SDSS filters"></svg></div>
<div class="row">
<div class="card"><h2>Hiding among the stars</h2><p class="sub">u−g against g−r: grey dots are stars, blue dots the median colour of DR16Q
quasars per Δz = 0.05, the orange line our spectrum's track.</p>
<svg id="cc" viewBox="0 0 520 400" role="img" aria-label="colour-colour diagram"></svg></div>
<div class="card"><h2>The cosmic clock</h2><p class="sub">Luminous quasars (M<sub>i</sub>(z=2) &lt; −27.6) per Gpc³, SDSS DR7 uniform sample;
filled points complete, open points lower limits; orange: our fit.</p>
<svg id="clock" viewBox="0 0 520 300" role="img" aria-label="density of luminous quasars against cosmic age"></svg>
<svg id="sel" viewBox="0 0 520 124" role="img" aria-label="selection completeness against redshift"></svg></div>
</div>
<p class="note">Data: SDSS DR16Q (Lyke et al. 2020) with Wu &amp; Shen (2022) measurements; SDSS DR7 uniform sample (Shen et al. 2011);
selection function from Richards et al. (2006); stars from SDSS DR17 SkyServer; cosmology from project D3. Masses are single-epoch
estimates, uncertain by ~0.4 dex each.</p>
</main>
<script>
const D = __PAYLOAD__;
const $ = id => document.getElementById(id);
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function interp(x, xs, ys){ if (x <= xs[0]) return ys[0]; const n = xs.length; if (x >= xs[n-1]) return ys[n-1];
  let lo = 0, hi = n - 1; while (hi - lo > 1){ const mid = (lo + hi) >> 1; if (xs[mid] > x) hi = mid; else lo = mid; }
  const a = ys[lo], b = ys[hi]; if (a === null || b === null) return null; return a + (b - a) * (x - xs[lo]) / (xs[hi] - xs[lo]); }
function madau(lam, z){ // Madau (1995) mean IGM transmission, as in photometry.py
  let tau = 0; const zs1 = 1 + z;
  const S = [[1215.67,3.6e-3],[1025.72,1.7e-3],[972.54,1.2e-3],[949.74,9.3e-4]];
  for (const [lj, aj] of S){ if (lam < lj * zs1) tau += aj * Math.pow(lam / lj, 3.46); }
  if (lam < 911.75 * zs1){ const xc = lam / 911.75, xe = zs1;
    tau += 0.25*xc**3*(xe**0.46 - xc**0.46) + 9.4*xc**1.5*(xe**0.18 - xc**0.18) - 0.7*xc**3*(xc**-1.32 - xe**-1.32) - 0.023*(xe**1.68 - xc**1.68); }
  return Math.exp(-Math.min(Math.max(tau, 0), 50)); }
function age(z){ return interp(z, D.ages.z, D.ages.t); }
function inRange(x, xs, ys){ return (x < xs[0] || x > xs[xs.length - 1]) ? null : interp(x, xs, ys); }
function sc(d0, d1, r0, r1){ return v => r0 + (v - d0) * (r1 - r0) / (d1 - d0); }
function path(pts){ let s = "", pen = false; for (const p of pts){ if (p === null || !isFinite(p[0]) || !isFinite(p[1])){ pen = false; continue; }
  s += (pen ? "L" : "M") + p[0].toFixed(1) + "," + p[1].toFixed(1); pen = true; } return s; }
function fmt(v, d){ return v === null || !isFinite(v) ? "–" : v.toFixed(d); }

function drawSpec(z){
  const W = 1000, H = 300, L = 46, R = 12, T = 12, B = 34;
  const x = sc(3000, 11000, L, W - R), y = sc(-2.3, 0.15, H - B, T);
  let s = "";
  for (let lam = 3000; lam <= 11000; lam += 1000){ s += `<line x1="${x(lam)}" x2="${x(lam)}" y1="${T}" y2="${H-B}" stroke="${css('--grid')}" stroke-width="0.6"/>`;
    s += `<text class="tick" x="${x(lam)}" y="${H-B+14}" text-anchor="middle">${lam}</text>`; }
  s += `<text class="tick" x="${(L+W-R)/2}" y="${H-4}" text-anchor="middle">observed wavelength (Å)</text>`;
  const bands = ["u","g","r","i","z"];
  for (const b of bands){ const f = D.filters[b]; const pts = f.w.map((w, k) => [x(w), y(-2.3 + 0.9 * f.r[k])]);
    s += `<path d="${path(pts)}Z" fill="${css('--filt')}" stroke="${css('--filtline')}" stroke-width="0.8"/>`;
    const k = f.r.indexOf(Math.max(...f.r)); s += `<text class="lab" x="${x(f.w[k])}" y="${y(-2.3 + 0.92)}" text-anchor="middle">${b}</text>`; }
  const pts = []; let fmax = 0; const obs = [];
  for (let k = 0; k < D.comp.w.length; k++){ const lo = D.comp.w[k] * (1 + z); if (lo < 3000 || lo > 11000) continue;
    const f = D.comp.f[k] * madau(lo, z); obs.push([lo, f]); if (f > fmax) fmax = f; }
  for (const [lo, f] of obs){ pts.push([x(lo), y(Math.log10(Math.max(f / fmax, 1e-3)))]); }
  s += `<path d="${path(pts)}" fill="none" stroke="${css('--ink')}" stroke-width="1.2"/>`;
  for (const [lr, name] of D.lines){ const lo = lr * (1 + z); if (lo < 3050 || lo > 10950) continue;
    s += `<text class="lab" x="${x(lo)}" y="${T + 10}" text-anchor="middle">${name}</text><line x1="${x(lo)}" x2="${x(lo)}" y1="${T+14}" y2="${T+22}" stroke="${css('--muted')}"/>`; }
  const lyl = 911.75 * (1 + z);
  if (lyl > 3000 && lyl < 11000){ s += `<line x1="${x(lyl)}" x2="${x(lyl)}" y1="${T+26}" y2="${H-B}" stroke="${css('--muted')}" stroke-width="0.8"/>`;
    s += `<text class="tick" x="${x(lyl) + 4}" y="${T + 36}" text-anchor="start">◂ Lyman limit: bluer light absorbed</text>`; }
  for (let v = -2; v <= 0; v++) s += `<text class="tick" x="${L-6}" y="${y(v)+3}" text-anchor="end">${v === 0 ? "1" : "10" + (v === -1 ? "⁻¹" : "⁻²")}</text>`;
  $("spec").innerHTML = s; }

function drawCC(z){
  const W = 520, H = 400, L = 44, R = 12, T = 10, B = 34;
  const x = sc(-0.6, 2.0, L, W - R), y = sc(-0.6, 4.4, H - B, T);
  let s = "";
  for (let v = -0.5; v <= 2.0; v += 0.5) s += `<line x1="${x(v)}" x2="${x(v)}" y1="${T}" y2="${H-B}" stroke="${css('--grid')}" stroke-width="0.6"/><text class="tick" x="${x(v)}" y="${H-B+14}" text-anchor="middle">${v.toFixed(1)}</text>`;
  for (let v = 0; v <= 4; v += 1) s += `<line x1="${L}" x2="${W-R}" y1="${y(v)}" y2="${y(v)}" stroke="${css('--grid')}" stroke-width="0.6"/><text class="tick" x="${L-6}" y="${y(v)+3}" text-anchor="end">${v}</text>`;
  s += `<text class="tick" x="${(L+W-R)/2}" y="${H-4}" text-anchor="middle">g − r</text><text class="tick" x="12" y="${(T+H-B)/2}" transform="rotate(-90 12 ${(T+H-B)/2})" text-anchor="middle">u − g</text>`;
  const sc_ = css('--star'); for (const [gr, ug] of D.stars){ if (gr < -0.6 || gr > 2 || ug > 4.4) continue; s += `<circle cx="${x(gr).toFixed(1)}" cy="${y(ug).toFixed(1)}" r="1.1" fill="${sc_}"/>`; }
  const bl = css('--blue'); for (let k = 0; k < D.obs.z.length; k++){ const ug = D.obs.ug[k], gr = D.obs.gr[k]; if (ug === null || gr === null || ug > 4.4) continue;
    s += `<circle cx="${x(gr).toFixed(1)}" cy="${y(ug).toFixed(1)}" r="2.6" fill="${bl}"><title>DR16Q median, z = ${D.obs.z[k].toFixed(2)}</title></circle>`; }
  const tp = D.track.z.map((zz, k) => (D.track.ug[k] === null || D.track.gr[k] === null || D.track.ug[k] > 4.4) ? null : [x(D.track.gr[k]), y(D.track.ug[k])]);
  s += `<path d="${path(tp)}" fill="none" stroke="${css('--orange')}" stroke-width="1.6"/>`;
  const ug = inRange(z, D.track.z, D.track.ug), gr = inRange(z, D.track.z, D.track.gr);
  if (ug !== null && gr !== null && ug <= 4.4) s += `<circle cx="${x(gr)}" cy="${y(ug)}" r="7" fill="${css('--orange')}" stroke="${css('--panel')}" stroke-width="2.5"/>`;
  else s += `<text class="lab" x="${x(1.2)}" y="${y(4.1)}">u − g off the scale: the quasar has dropped out of the u band</text>`;
  $("cc").innerHTML = s; }

function drawClock(z){
  const W = 520, H = 300, L = 46, R = 12, T = 10, B = 34;
  const x = sc(0, 6.5, L, W - R), y = sc(-1.0, 2.5, H - B, T);
  let s = "";
  for (let v = 0; v <= 6; v++) s += `<line x1="${x(v)}" x2="${x(v)}" y1="${T}" y2="${H-B}" stroke="${css('--grid')}" stroke-width="0.6"/><text class="tick" x="${x(v)}" y="${H-B+14}" text-anchor="middle">${v}</text>`;
  for (const v of [-1, 0, 1, 2]) s += `<line x1="${L}" x2="${W-R}" y1="${y(v)}" y2="${y(v)}" stroke="${css('--grid')}" stroke-width="0.6"/><text class="tick" x="${L-6}" y="${y(v)+3}" text-anchor="end">${v === 0 ? "1" : (v === 1 ? "10" : (v === 2 ? "100" : "0.1"))}</text>`;
  s += `<text class="tick" x="${(L+W-R)/2}" y="${H-4}" text-anchor="middle">age of the universe (billion years)</text>`;
  const c = D.era.fit, fit = [];
  for (let zz = D.era.zmin; zz <= D.era.zmax + 1e-9; zz += 0.02){ const xx = Math.log10(1 + zz); fit.push([x(age(zz) / 1e3), y(c[0]*xx*xx + c[1]*xx + c[2])]); }
  s += `<path d="${path(fit)}" fill="none" stroke="${css('--orange')}" stroke-width="1.6"/>`;
  const bl = css('--blue');
  for (let k = 0; k < D.era.z.length; k++){ const n = D.era.n[k]; if (!(n > 0)) continue; const xx = x(age(D.era.z[k]) / 1e3), yy = y(Math.log10(n));
    s += D.era.complete[k] ? `<circle cx="${xx}" cy="${yy}" r="4" fill="${bl}"><title>z = ${D.era.z[k].toFixed(2)}: ${n.toFixed(1)} per Gpc³</title></circle>`
                           : `<circle cx="${xx}" cy="${yy}" r="4" fill="${css('--panel')}" stroke="${bl}" stroke-width="1.4"><title>z = ${D.era.z[k].toFixed(2)}: at least ${n.toFixed(1)} per Gpc³</title></circle>`; }
  const t = age(z) / 1e3, tc = Math.min(t, 6.5);
  s += `<line x1="${x(tc)}" x2="${x(tc)}" y1="${T}" y2="${H-B}" stroke="${css('--ink')}" stroke-width="1.2"/>`;
  s += `<text class="lab" x="${Math.min(x(tc) + 5, W - 160)}" y="${T + 12}">${t > 6.5 ? "z = " + z.toFixed(2) + " is later: " + t.toFixed(1) + " Gyr →" : "now: z = " + z.toFixed(2)}</text>`;
  $("clock").innerHTML = s;
  // selection completeness strip
  const W2 = 520, H2 = 124, B2 = 26; const xs = sc(0, 7.7, L, W2 - R), ys = sc(0, 1, H2 - B2, 26);
  let s2 = `<text class="tick" x="${L}" y="14">fraction of quasars the SDSS colour selection catches (Richards et al. 2006, i &lt; 19.1)</text>`;
  for (let v = 0; v <= 7; v++) s2 += `<text class="tick" x="${xs(v)}" y="${H2-B2+13}" text-anchor="middle">${v}</text>`;
  s2 += `<text class="tick" x="${(L+W2-R)/2}" y="${H2-2}" text-anchor="middle">redshift</text>`;
  s2 += `<line x1="${L}" x2="${W2-R}" y1="${ys(0)}" y2="${ys(0)}" stroke="${css('--axis')}"/><line x1="${L}" x2="${W2-R}" y1="${ys(1)}" y2="${ys(1)}" stroke="${css('--grid')}"/>`;
  s2 += `<path d="${path(D.track.z.map((zz, k) => D.track.sel[k] === null ? null : [xs(zz), ys(D.track.sel[k])]))}" fill="none" stroke="${css('--violet')}" stroke-width="1.6"/>`;
  s2 += `<line x1="${xs(z)}" x2="${xs(z)}" y1="22" y2="${H2-B2}" stroke="${css('--ink')}" stroke-width="1.2"/>`;
  $("sel").innerHTML = s2; }

function readout(z){
  const t = age(z), sel = inRange(z, D.track.z, D.track.sel), dist = inRange(z, D.track.z, D.track.dist);
  const ug = inRange(z, D.track.z, D.track.ug), gr = inRange(z, D.track.z, D.track.gr), ri = inRange(z, D.track.z, D.track.ri);
  let dens = "–";
  if (z >= D.era.zmin && z <= D.era.zmax){ const xx = Math.log10(1 + z), c = D.era.fit; dens = Math.pow(10, c[0]*xx*xx + c[1]*xx + c[2]).toFixed(0) + " per Gpc³"; }
  else if (z < D.era.zmin) dens = "not measured (too bright for SDSS's i = 15 limit this close)"; else dens = "not measured (too faint for the survey this far)";
  let mp = null, m95 = null, msrc = "DR16Q black holes here";
  if (z <= D.mass.z[D.mass.z.length - 1]){ mp = interp(z, D.mass.z, D.mass.p50); m95 = interp(z, D.mass.z, D.mass.p95); }
  else { const near = D.highz.z.map((zz, k) => [zz, D.highz.m[k]]).filter(([zz, mm]) => Math.abs(zz - z) < 0.3 && mm !== null).map(r => r[1]).sort((a, b) => a - b);
    msrc = "Known quasars within Δz = 0.3 (" + near.length + ")";
    if (near.length){ mp = near[Math.floor(near.length / 2)]; m95 = near[near.length - 1]; } }
  const dt = t - D.t30, ef = dt / D.efold;
  const seedlog = m95 === null ? null : m95 - ef / Math.LN10;
  let seed = "–";
  if (seedlog !== null && dt > 0) seed = seedlog < 1 ? "any stellar remnant: time is not the problem here" : "≥ 10<sup>" + seedlog.toFixed(1) + "</sup> M☉";
  const lya = 1215.67 * (1 + z);
  let band = "outside the filters"; for (const b of ["u","g","r","i","z"]){ const f = D.filters[b]; let k = f.r.indexOf(Math.max(...f.r)); const half = f.w.filter((w, j) => f.r[j] > 0.5); if (lya >= half[0] && lya <= half[half.length-1]) band = "in the " + b + " band"; }
  const inside = dist !== null && dist < 0.06 && z < 3.0;
  $("read").innerHTML =
    `<div>Universe's age: <b>${(t/1e3).toFixed(2)} billion years</b> (lookback ${(13.825 - t/1e3).toFixed(2)})</div>` +
    `<div>Lyman α lands at <b>${lya.toFixed(0)} Å</b>, ${band}</div>` +
    `<div>Colours: u−g <b>${fmt(ug, 2)}</b>, g−r <b>${fmt(gr, 2)}</b>, r−i <b>${fmt(ri, 2)}</b></div>` +
    `<div>SDSS colour selection catches <b>${sel === null ? "– (beyond its range)" : (100 * sel).toFixed(0) + "%"}</b>${inside ? ' <span class="warn">– inside the stellar locus</span>' : ""}</div>` +
    `<div>Luminous quasars: <b>${dens}</b> (peak z = ${D.era.zpeak.toFixed(2)})</div>` +
    `<div>${msrc}: median <b>10<sup>${fmt(mp, 1)}</sup> M☉</b>, ${z <= D.mass.z[D.mass.z.length - 1] ? "top 5% above" : "heaviest"} 10<sup>${fmt(m95, 1)}</sup></div>` +
    `<div>Since z = 30: <b>${dt > 0 ? dt.toFixed(0) + " Myr = " + ef.toFixed(1) + " e-folds" : "–"}</b> of Eddington-limited growth</div>` +
    `<div>Seed needed for the heaviest: <b>${seed}</b></div>`;
  $("zread").textContent = "z = " + z.toFixed(2); }

function draw(){ const z = parseFloat($("z").value); drawSpec(z); drawCC(z); drawClock(z); readout(z); }
$("z").addEventListener("input", draw);
document.querySelectorAll("button[data-z]").forEach(b => b.addEventListener("click", () => { $("z").value = b.dataset.z; draw(); }));
$("theme").addEventListener("click", () => { const cur = document.documentElement.getAttribute("data-theme");
  const dark = cur ? cur === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
  const nt = dark ? "light" : "dark"; document.documentElement.setAttribute("data-theme", nt); try { localStorage.setItem("i3theme", nt); } catch (e) {} draw(); });
try { const t = localStorage.getItem("i3theme"); if (t) document.documentElement.setAttribute("data-theme", t); } catch (e) {}
draw();
</script></body></html>
"""


def build(res: dict):
    payload = json.dumps(build_payload(res), separators=(",", ":"), allow_nan=False)
    out = FIGURES / "quasar_time_machine.html"
    out.write_text(HTML.replace("__PAYLOAD__", payload), encoding="utf-8")
    return out


if __name__ == "__main__":
    import pickle

    from config import RESULTS_PICKLE
    with open(RESULTS_PICKLE, "rb") as fh:
        res = pickle.load(fh)
    p = build(res)
    print(p, f"{p.stat().st_size / 1e3:.0f} kB")
