"""An offline page to see and hear stars ring: figures/star_sounds.html.

Pick a star, play its (sped-up) light curve, and drag the large separation
Dnu: the echelle diagram straightens into vertical ridges only at the right
value. The page is self-contained (no network, no libraries): spectra and
sounds embedded, drawn on canvas, light and dark themes, a table view of all
104 red giants.
"""

from __future__ import annotations

import base64
import json

import numpy as np
import pandas as pd

import config as C
import sonify as SO
import spectrum as S
import stars


def _b64_f32(a) -> str:
    return base64.b64encode(np.asarray(a, np.float32).tobytes()).decode("ascii")


def _entry(key, label, kind, nu, snr, numax, dnu, eps, sigma, facts, audio_b64, audio_info, n_sig=3.0):
    lo, hi = numax - n_sig * sigma, numax + n_sig * sigma
    m = (nu >= lo) & (nu <= hi)
    nu_m, s_m = nu[m], snr[m]
    df = float(nu_m[1] - nu_m[0])
    # a lightly smoothed copy for the spectrum panel (<= 4000 points)
    k = max(1, int(np.ceil(len(s_m) / 4000)))
    sm = S.smooth_boxcar(s_m, max(1, k))
    return {"key": key, "label": label, "kind": kind, "f0": float(nu_m[0]), "df": df,
            "snr": _b64_f32(s_m), "spec_f": _b64_f32(nu_m[::k]), "spec_s": _b64_f32(sm[::k]),
            "numax": float(numax), "dnu": float(dnu), "eps": float(eps), "sigma": float(sigma),
            "facts": facts, "audio": audio_b64, "audio_info": audio_info}


def _facts_giant(row, js):
    out = [("νmax", f"{js['numax']:.2f} ± {js['numax_err']:.2f} µHz"),
           ("Δν", f"{js['dnu']:.3f} ± {js['dnu_err']:.3f} µHz"),
           ("Radius", f"{row['R_self']:.1f} ± {row['R_self_err']:.1f} R☉"),
           ("Mass", f"{row['M_self']:.2f} ± {row['M_self_err']:.2f} M☉"),
           ("Luminosity", f"{row['L_self']:.0f} L☉"),
           ("Temperature (APOGEE)", f"{row['teff_ak']:.0f} K"),
           ("Distance (Gaia)", f"{row['distance_pc']:.0f} pc")]
    if row["dpi1_detected"]:
        burn = "hydrogen in a shell" if row["dpi1_hypothesis"] == "RGB" else "helium in the core"
        out.append(("Core (from ΔΠ₁)", f"ΔΠ₁ = {row['dpi1']:.1f} s: burning {burn}"))
    else:
        out.append(("Core (from ΔΠ₁)", "period spacing not measured"))
    out.append(("APOKASC-3 state", {"RGB": "red giant branch", "RC": "red clump"}.get(row["state_apokasc"],
                                                                                   row["state_apokasc"])))
    return out


def build(kics, res: dict, *, quiet: bool = False):
    df = pd.read_csv(C.DATA_DIR / "results_stars.csv", index_col="kic")
    audio = json.loads((C.AUDIO_DIR / "audio.json").read_text())
    entries = []
    from scipy.io import wavfile

    def b64(fname):
        return base64.b64encode((C.AUDIO_DIR / fname).read_bytes()).decode("ascii")

    for kic in kics:
        js, arr = stars.load_star(int(kic))
        row = df.loc[int(kic)]
        a = audio[str(int(kic))]
        kind = "RC" if row["state_apokasc"] == "RC" else "RGB"
        label = f"KIC {int(kic)}"
        entries.append(_entry(str(int(kic)), label, kind, arr["nu"], arr["psd"] / arr["bg"], js["numax"],
                              js["dnu"], js["epsilon"], js["envelope_fwhm"] / 2.3548,
                              _facts_giant(row, js), b64(a["file"]),
                              {"pitch": round(a["pitch_hz"]), "note": a["note"], "speedup": a["speedup"],
                               "days": round(a["data_days"])}))
    # the Sun and 16 Cyg A
    sun = json.loads((C.CACHE_DIR / "sun.json").read_text())
    z = np.load(C.CACHE_DIR / "sun.npz")
    y08 = sun["years"]["2008"]
    a = audio["sun"]
    entries.append(_entry("sun", "The Sun", "MS", z["y2008_nu"], z["y2008_psd"] / z["y2008_bg"], y08["numax"],
                          y08["dnu"], y08["epsilon"], y08["envelope_fwhm"] / 2.3548,
                          [("νmax (2008)", f"{y08['numax']:.0f} ± {y08['numax_err']:.0f} µHz"),
                           ("Δν (2008)", f"{y08['dnu']:.2f} ± {y08['dnu_err']:.2f} µHz"),
                           ("Radius, mass", "1 R☉, 1 M☉ by definition: the reference point"),
                           ("Data", "SOHO/VIRGO green sun photometer, 2008")],
                          b64(a["file"]), {"pitch": round(a["pitch_hz"]), "note": a["note"],
                                           "speedup": a["speedup"], "days": round(a["data_days"], 1)},
                          n_sig=2.0))
    ms = json.loads((C.CACHE_DIR / "sixteen_cyg.json").read_text())["16 Cyg A"]
    zc = np.load(C.CACHE_DIR / "sixteen_cyg.npz")
    a = audio["16 Cyg A"]
    sc = res["sixteen_cyg"]["16 Cyg A"]["scaling"]["self"]
    entries.append(_entry("16cyga", "16 Cyg A", "MS", zc["16_cyg_a_nu"], zc["16_cyg_a_psd"] / zc["16_cyg_a_bg"],
                          ms["numax"], ms["dnu"], ms["epsilon"], ms["envelope_fwhm"] / 2.3548,
                          [("νmax", f"{ms['numax']:.0f} ± {ms['numax_err']:.0f} µHz"),
                           ("Δν", f"{ms['dnu']:.2f} ± {ms['dnu_err']:.2f} µHz"),
                           ("Radius", f"{sc['radius']:.3f} ± {sc['radius_err']:.3f} R☉ (interferometry: 1.22 ± 0.02)"),
                           ("Mass", f"{sc['mass']:.2f} ± {sc['mass_err']:.2f} M☉ (models: 1.08 ± 0.02)"),
                           ("Data", "Kepler short cadence, Quarters 7-10")],
                          b64(a["file"]), {"pitch": round(a["pitch_hz"]), "note": a["note"],
                                           "speedup": a["speedup"], "days": round(a["data_days"], 1)},
                          n_sig=1.6))
    table = []
    for kic, r in df.sort_values("numax").iterrows():
        table.append([int(kic), round(r["numax"], 2), round(r["numax_ak"], 2), round(r["dnu"], 3),
                      round(r["dnu_ak"], 3), round(r["M_self"], 2), round(r["R_self"], 2),
                      round(r["M_ak"], 2), round(r["mass_ak"], 2),
                      ("" if not r["dpi1_detected"] else round(r["dpi1"], 1)),
                      r["state_apokasc"]])
    payload = {"stars": entries, "table": table, "speedup_giants": SO.SPEEDUP_GIANTS,
               "speedup_dwarfs": SO.SPEEDUP_DWARFS}
    html = PAGE.replace("__DATA__", json.dumps(payload, separators=(",", ":")))
    out = C.FIGURE_DIR / "star_sounds.html"
    out.write_text(html, encoding="utf-8")
    if not quiet:
        print(f"  {out} ({out.stat().st_size / 1e6:.1f} MB)")
    return out


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Star Sounds</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --s1: #2a78d6; --s2: #eb6834;
  --r0: #fcfcfb; --r1: #cde2fb; --r2: #86b6ef; --r3: #3987e5; --r4: #1c5cab; --r5: #0d366b;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
    --grid: #2f2f2d; --axis: #4a4a46; --s1: #3987e5; --s2: #d95926;
    --r0: #232321; --r1: #1b2c44; --r2: #184f95; --r3: #3987e5; --r4: #86b6ef; --r5: #e8f1fd;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --panel: #232321; --ink: #ffffff; --ink2: #c3c2b7; --muted: #8f8e88;
  --grid: #2f2f2d; --axis: #4a4a46; --s1: #3987e5; --s2: #d95926;
  --r0: #232321; --r1: #1b2c44; --r2: #184f95; --r3: #3987e5; --r4: #86b6ef; --r5: #e8f1fd;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1100px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 16px; max-width: 800px; }
.row { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 12px; }
button { font: inherit; font-size: 0.86rem; padding: 5px 11px; border-radius: 8px; cursor: pointer;
         border: 1px solid var(--axis); background: var(--surface); color: var(--ink); }
button[aria-pressed="true"] { background: var(--ink); color: var(--surface); border-color: var(--ink); }
button:focus-visible, input:focus-visible { outline: 2px solid var(--s1); outline-offset: 2px; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 300px; gap: 16px; }
@media (max-width: 860px) { .layout { grid-template-columns: minmax(0, 1fr); } }
.card { background: var(--panel); border-radius: 12px; padding: 12px; min-width: 0; position: relative; }
canvas { width: 100%; display: block; border-radius: 6px; }
.ctl { margin-bottom: 14px; }
.ctl label { display: flex; justify-content: space-between; font-size: 0.88rem; color: var(--ink2); }
.ctl label b { color: var(--ink); font-variant-numeric: tabular-nums; }
.ctl input[type=range] { width: 100%; }
.readout { font-size: 0.86rem; color: var(--ink2); border-top: 1px solid var(--grid); padding-top: 10px; }
.readout div { margin-bottom: 5px; }
.readout b { color: var(--ink); }
.cap { font-size: 0.8rem; color: var(--muted); margin: 4px 2px 10px; }
.note { font-size: 0.86rem; color: var(--ink2); margin: 8px 2px 0; }
.tip { position: absolute; pointer-events: none; background: var(--surface); color: var(--ink);
       border: 1px solid var(--axis); border-radius: 6px; padding: 3px 7px; font-size: 0.78rem;
       font-variant-numeric: tabular-nums; display: none; white-space: nowrap; }
audio { width: 100%; margin: 6px 0 4px; }
details { margin-top: 16px; font-size: 0.85rem; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
td, th { padding: 2px 9px; text-align: right; border-bottom: 1px solid var(--grid); }
th { color: var(--ink2); font-weight: 600; }
.tablewrap { max-height: 320px; overflow: auto; }
footer { color: var(--muted); font-size: 0.78rem; margin-top: 18px; }
</style>
</head>
<body>
<main>
  <h1>Star sounds</h1>
  <p class="lede">Stars ring like bells. Turbulent convection near the surface makes them vibrate with sound
  waves, which make their brightness flicker by a few parts in ten thousand. Kepler recorded that flicker for
  four years. Below are real Kepler light curves of red giants, sped up twenty million times so you can hear
  them: the bigger the star, the lower its note. The pattern of tones, spaced by the <b>large separation Δν</b>,
  measures the star's mean density. Drag Δν until the échelle diagram's ridges stand up straight.</p>
  <div class="row" id="picker" role="group" aria-label="stars"></div>
  <div class="layout">
    <div class="card">
      <canvas id="spec" height="190" aria-label="power spectrum around the oscillation hump"></canvas>
      <div class="cap" id="speccap"></div>
      <canvas id="ech" height="420" aria-label="echelle diagram"></canvas>
      <div class="cap">Échelle diagram: the spectrum cut into strips Δν wide and stacked. Stronger colour means more
      oscillation power. At the right Δν the modes of each degree line up vertically.</div>
      <div class="tip" id="tip"></div>
    </div>
    <div class="card">
      <div id="title" style="font-weight:600;font-size:1.05rem;margin-bottom:4px"></div>
      <audio id="player" controls preload="none"></audio>
      <div class="note" id="audionote"></div>
      <div class="ctl" style="margin-top:12px"><label for="dnu">Large separation Δν <b id="dnuv"></b></label>
        <input id="dnu" type="range" min="-30" max="30" step="0.25" value="0"></div>
      <div class="ctl"><label for="shift">Shift the strips <b id="shiftv"></b></label>
        <input id="shift" type="range" min="0" max="1" step="0.01" value="0.3"></div>
      <div class="row"><button id="snap">jump to the measured Δν</button></div>
      <div class="readout" id="facts"></div>
    </div>
  </div>
  <details><summary>Table view: all 104 red giants</summary><div class="tablewrap"><table id="tbl"></table></div></details>
  <footer>Data: Kepler long and short cadence (NASA, via MAST); SOHO/VIRGO (ESA/NASA). numax, Δν, masses and radii
  are this project's own measurements: masses and radii from the scaling relations calibrated on our own
  measurement of the Sun, with APOGEE temperatures. The sound is the measured light curve band-passed around the
  oscillation hump and normalised in loudness (loudness does not compare between stars); red giants are sped up
  2 × 10⁷ times, the Sun and 16 Cyg A 10⁵ times. Evolutionary states: APOKASC-3 (Pinsonneault et al. 2025).</footer>
</main>
<script>
const DATA = __DATA__;
const $ = id => document.getElementById(id);
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function f32(b64) { const b = atob(b64); const u = new Uint8Array(b.length);
  for (let i = 0; i < b.length; i++) u[i] = b.charCodeAt(i); return new Float32Array(u.buffer); }
const cache = {};
function star(i) { if (!cache[i]) { const s = DATA.stars[i];
  cache[i] = { ...s, snrA: f32(s.snr), sf: f32(s.spec_f), ss: f32(s.spec_s) }; } return cache[i]; }
let cur = 0, dnuNow = 0;
function hexToRgb(h) { h = h.replace('#', ''); if (h.length === 3) h = h.split('').map(c => c + c).join('');
  const n = parseInt(h, 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
function ramp() { return ['--r0', '--r1', '--r2', '--r3', '--r4', '--r5'].map(v => hexToRgb(css(v))); }
function rampColor(R, t) { t = Math.max(0, Math.min(1, t)) * (R.length - 1); const i = Math.min(R.length - 2, Math.floor(t));
  const f = t - i; return R[i].map((c, k) => Math.round(c + (R[i + 1][k] - c) * f)); }
function setupCanvas(c) { const r = window.devicePixelRatio || 1; const w = c.clientWidth, h = c.height;
  if (!c.dataset.h) c.dataset.h = h; const H = +c.dataset.h; c.width = w * r; c.height = H * r; c.style.height = H + 'px';
  const ctx = c.getContext('2d'); ctx.setTransform(r, 0, 0, r, 0, 0); return [ctx, w, H]; }
const PAD = { l: 58, r: 12, t: 10, b: 34 };
function axes(ctx, w, h, x0, x1, y0, y1, xlab, ylab) {
  ctx.strokeStyle = css('--axis'); ctx.lineWidth = 1; ctx.fillStyle = css('--muted'); ctx.font = '11px system-ui';
  ctx.beginPath(); ctx.moveTo(PAD.l, PAD.t); ctx.lineTo(PAD.l, h - PAD.b); ctx.lineTo(w - PAD.r, h - PAD.b); ctx.stroke();
  const nt = 6; ctx.textAlign = 'center';
  for (let i = 0; i <= nt; i++) { const v = x0 + (x1 - x0) * i / nt; const x = PAD.l + (w - PAD.l - PAD.r) * i / nt;
    ctx.fillText(v.toFixed(v < 10 ? 2 : (v < 100 ? 1 : 0)), x, h - PAD.b + 14); }
  ctx.textAlign = 'right';
  for (let i = 0; i <= 4; i++) { const v = y0 + (y1 - y0) * i / 4; const y = h - PAD.b - (h - PAD.t - PAD.b) * i / 4;
    ctx.fillText(v.toFixed(v < 10 ? 1 : 0), PAD.l - 5, y + 4); }
  ctx.textAlign = 'center'; ctx.fillStyle = css('--ink2'); ctx.fillText(xlab, (w + PAD.l) / 2, h - 4);
  ctx.save(); ctx.translate(11, PAD.t + (h - PAD.t - PAD.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(ylab, 0, 0);
  ctx.restore(); ctx.textAlign = 'left'; }
function drawSpec() {
  const s = star(cur); const c = $('spec'); const [ctx, w, h] = setupCanvas(c);
  ctx.clearRect(0, 0, w, h);
  const x0 = s.sf[0], x1 = s.sf[s.sf.length - 1]; let ymax = 0; for (const v of s.ss) ymax = Math.max(ymax, v);
  ymax *= 1.08; axes(ctx, w, h, x0, x1, 0, ymax, 'Frequency (µHz)', 'Power / background');
  const X = v => PAD.l + (v - x0) / (x1 - x0) * (w - PAD.l - PAD.r), Y = v => h - PAD.b - v / ymax * (h - PAD.t - PAD.b);
  // comb of l=0 modes at the current Dnu and shift
  const sh = +$('shift').value; ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1;
  const off = (s.eps - 0.3 + sh + 10) % 1;
  for (let n = Math.floor(x0 / dnuNow) - 1; n * dnuNow < x1 + dnuNow; n++) { const f = (n + off) * dnuNow + 0.3 * dnuNow;
    if (f > x0 && f < x1) { ctx.beginPath(); ctx.moveTo(X(f), PAD.t); ctx.lineTo(X(f), h - PAD.b); ctx.stroke(); } }
  ctx.strokeStyle = css('--s1'); ctx.lineWidth = 1.2; ctx.beginPath();
  for (let i = 0; i < s.sf.length; i++) { const x = X(s.sf[i]), y = Y(s.ss[i]); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }
  ctx.stroke();
  $('speccap').textContent = `Grey lines: a comb spaced by the Δν on the slider (${dnuNow.toFixed(3)} µHz). The hump peaks at νmax = ${s.numax.toFixed(s.numax < 100 ? 2 : 0)} µHz.`;
}
function drawEch() {
  const s = star(cur); const c = $('ech'); const [ctx, w, h] = setupCanvas(c);
  ctx.clearRect(0, 0, w, h); const R = ramp();
  const sh = +$('shift').value; const d = dnuNow; const off = ((s.eps - sh) % 1 + 1) % 1;
  const start = (Math.floor(s.f0 / d - off) + off) * d; const nx = 150;
  const nrow = Math.ceil((s.f0 + s.df * s.snrA.length - start) / d);
  const img = new Float32Array(nrow * nx), cnt = new Float32Array(nrow * nx);
  const sm = 2; // light smoothing
  for (let i = 0; i < s.snrA.length; i++) { const f = s.f0 + i * s.df; const rel = f - start; const r = Math.floor(rel / d);
    let col = Math.floor((rel - r * d) / d * nx); if (col >= nx) col = nx - 1; if (r < 0 || r >= nrow) continue;
    let v = 0; for (let k = -sm; k <= sm; k++) { const j = Math.min(s.snrA.length - 1, Math.max(0, i + k)); v += s.snrA[j]; }
    img[r * nx + col] += v / (2 * sm + 1); cnt[r * nx + col] += 1; }
  const vals = []; for (let i = 0; i < img.length; i++) { if (cnt[i]) { img[i] /= cnt[i]; vals.push(img[i]); } }
  vals.sort((a, b) => a - b); const vmax = vals[Math.floor(vals.length * 0.99)] || 1;
  const pw = w - PAD.l - PAD.r, ph = h - PAD.t - PAD.b; const cw = pw / nx, rh = ph / nrow;
  for (let r = 0; r < nrow; r++) for (let k = 0; k < nx; k++) { const t = img[r * nx + k] / vmax;
    const [cr, cg, cb] = rampColor(R, t); ctx.fillStyle = `rgb(${cr},${cg},${cb})`;
    ctx.fillRect(PAD.l + k * cw, h - PAD.b - (r + 1) * rh, cw + 0.6, rh + 0.6); }
  axes(ctx, w, h, 0, d, start, start + nrow * d, 'Frequency modulo Δν (µHz)', 'Frequency (µHz)');
  c._map = { start, d, nrow, nx, w, h };
}
function draw() { drawSpec(); drawEch(); }
function select(i) {
  cur = i; const s = star(i);
  document.querySelectorAll('#picker button').forEach((b, k) => b.setAttribute('aria-pressed', k === i));
  $('dnu').value = 0; dnuNow = s.dnu; $('shift').value = 0.3; update();
  $('title').textContent = s.label + (s.kind === 'RC' ? ' — red clump' : s.kind === 'RGB' ? ' — red giant branch' : '');
  $('player').src = 'data:audio/wav;base64,' + s.audio;
  const a = s.audio_info; $('audionote').innerHTML = `Sped up ${a.speedup.toExponential(0).replace('e+', ' × 10^')}: νmax becomes <b>${a.pitch} Hz</b> (${a.note}). Five seconds of sound = ${a.days} days of starlight.`;
  $('facts').innerHTML = s.facts.map(([k, v]) => `<div>${k}: <b>${v}</b></div>`).join('');
}
function update() { const s = star(cur); dnuNow = s.dnu * (1 + (+$('dnu').value) / 1000);
  $('dnuv').textContent = `${dnuNow.toFixed(s.dnu < 10 ? 3 : 2)} µHz (${(+$('dnu').value / 10).toFixed(1)}%)`;
  $('shiftv').textContent = (+$('shift').value).toFixed(2); draw(); }
DATA.stars.forEach((s, i) => { const b = document.createElement('button');
  b.textContent = `${s.label} · ${s.numax.toFixed(0)} µHz`; b.title = `νmax ${s.numax.toFixed(1)} µHz`;
  b.onclick = () => select(i); $('picker').appendChild(b); });
$('dnu').oninput = update; $('shift').oninput = update; $('snap').onclick = () => { $('dnu').value = 0; update(); };
window.addEventListener('resize', draw);
// hover readout on the echelle
$('ech').addEventListener('mousemove', e => { const c = $('ech'), m = c._map; if (!m) return;
  const rect = c.getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top;
  const pw = m.w - PAD.l - PAD.r, ph = m.h - PAD.t - PAD.b; if (x < PAD.l || x > m.w - PAD.r || y < PAD.t || y > m.h - PAD.b) { $('tip').style.display = 'none'; return; }
  const fx = (x - PAD.l) / pw * m.d, row = Math.floor((m.h - PAD.b - y) / ph * m.nrow); const f = m.start + row * m.d + fx;
  const t = $('tip'); t.style.display = 'block'; t.textContent = `${f.toFixed(2)} µHz`;
  t.style.left = (x + 18) + 'px'; t.style.top = (c.offsetTop + y - 10) + 'px'; });
$('ech').addEventListener('mouseleave', () => $('tip').style.display = 'none');
// table view
const head = ['KIC', 'νmax', 'νmax AK', 'Δν', 'Δν AK', 'M', 'R', 'M (AK scale)', 'M AK', 'ΔΠ₁ (s)', 'state AK'];
$('tbl').innerHTML = '<tr>' + head.map(h => `<th>${h}</th>`).join('') + '</tr>' +
  DATA.table.map(r => '<tr>' + r.map(v => `<td>${v}</td>`).join('') + '</tr>').join('');
if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', draw);
select(0);
</script>
</body>
</html>
"""
