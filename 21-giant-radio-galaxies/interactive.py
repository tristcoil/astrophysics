"""An offline, self-contained page for exploring the giants.

Pick a giant; blend its LOFAR radio emission over the optical sky; drag across
the image to measure any distance, in arcminutes and in megaparsecs at the
giant's distance. Everything (images and numbers) is embedded, so the page
works without a network connection.
"""

from __future__ import annotations

import base64
import io
import json

import numpy as np
from PIL import Image

import config
import cosmology
import figures
import radiomaps as rm
from config import FIGURE_DIR

SIZE_PX = 640


def _png(arr8: np.ndarray, mode: str) -> str:
    buf = io.BytesIO()
    Image.fromarray(arr8, mode=mode).save(buf, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def _jpeg(rgb8: np.ndarray, quality: int = 84) -> str:
    buf = io.BytesIO()
    Image.fromarray(rgb8, mode="RGB").save(buf, format="JPEG", quality=quality, optimize=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def _resize(arr: np.ndarray, n: int) -> np.ndarray:
    if arr.ndim == 2:
        return np.asarray(Image.fromarray(arr.astype(np.float32)).resize((n, n), Image.BILINEAR))
    return np.dstack([_resize(arr[..., i], n) for i in range(arr.shape[2])])


def giant_payload(g: config.Giant, row: dict, ext: rm.Extent) -> dict:
    opt, wopt = figures.legacy_image(g)
    img, sig = figures.display_map(g)
    rad = rm.reproject_to(img, wopt, opt.shape[:2])
    layer = figures.radio_layer(rad, sig)
    mask = figures.giant_mask_on(g, ext, wopt, opt.shape[:2])
    cra, cdec = figures._midpoint(row["adopted"]["ends"])
    half = 0.5 * row["adopted"]["las_arcmin"] * 1.25 + 0.8
    x0, x1, y0, y1 = figures._crop_box(wopt, opt.shape, cra, cdec, half)
    side = min(x1 - x0, y1 - y0)
    x1, y1 = x0 + side, y0 + side
    pix_arcmin = abs(wopt.pixel_scale_matrix[1, 1]) * 60 * side / SIZE_PX
    crop = lambda a: a[y0:y1, x0:x1][::-1]  # row 0 at the top, for the browser
    o8 = (np.clip(_resize(crop(opt), SIZE_PX), 0, 1) * 255).astype(np.uint8)
    lg = (np.clip(_resize(crop(layer * mask), SIZE_PX), 0, 1) * 255).astype(np.uint8)
    lo = (np.clip(_resize(crop(layer * ~mask), SIZE_PX), 0, 1) * 255).astype(np.uint8)
    hx, hy = wopt.world_to_pixel_values(g.ra, g.dec)
    host = [(hx - x0) / side * SIZE_PX, (y1 - 1 - hy) / side * SIZE_PX]
    return {
        "key": g.key, "name": g.name + (" (DR3)" if g.survey == "dr3" else ""), "host": g.host, "why": g.why,
        "z": g.z, "z_kind": g.z_kind,
        "z_ref": g.z_ref, "lp": row["adopted"]["lp_mpc"], "las": row["adopted"]["las_arcmin"],
        "beam": g.measure_beam,
        "kpc_per_arcmin": float(cosmology.kpc_per_arcmin(g.z, config.PLANCK18)),
        "lookback_gyr": float(cosmology.lookback_time(g.z, config.PLANCK18)),
        "pix_arcmin": pix_arcmin, "host_px": [float(host[0]), float(host[1])],
        "claims": [{"source": c[0], "las": c[1], "lp": c[2]} for c in g.claims],
        "optical": _jpeg(o8), "radio": _png(lg, "L"), "other": _png(lo, "L"),
    }


def write_explorer(measured: list[dict], art: dict) -> str:
    rows = {r["key"]: r for r in measured}
    giants = sorted([g for g in config.GIANTS if "gallery" in g.groups or "dr3" in g.groups],
                    key=lambda g: -rows[g.key]["adopted"]["lp_mpc"])
    payload = {"giants": [giant_payload(g, rows[g.key], figures.display_extent(g, art)) for g in giants],
               "local_group": {"mw_kpc": config.MILKY_WAY_DIAMETER_KPC, "m31_kpc": config.ANDROMEDA_DIAMETER_KPC,
                               "sep_kpc": config.ANDROMEDA_DISTANCE_KPC},
               "radio_rgb": [255, 122, 46]}
    html = TEMPLATE.replace("/*PAYLOAD*/null", json.dumps(payload, separators=(",", ":")))
    out = FIGURE_DIR / "giant_explorer.html"
    out.write_text(html, encoding="utf-8")
    return str(out)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Giant Radio Galaxies</title>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --rule: #e1e0d9; --accent: #eb6834; --accent-ink: #b8461c; --sel: #2a78d6; --well: #000000;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
    --rule: #2c2c2a; --accent: #d95926; --accent-ink: #f08a5d; --sel: #5598e7; --well: #000000;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
  --rule: #2c2c2a; --accent: #d95926; --accent-ink: #f08a5d; --sel: #5598e7; --well: #000000;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 20px 16px 40px; }
h1 { font-size: 1.45rem; margin: 0 0 4px; }
.lede { color: var(--ink2); margin: 0 0 16px; max-width: 75ch; }
.picker { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 14px; }
.picker button { font: inherit; font-size: .85rem; padding: 5px 10px; border-radius: 5px; cursor: pointer;
  background: var(--surface); color: var(--ink); border: 1px solid var(--rule); }
.picker button[aria-pressed="true"] { border-color: var(--sel); box-shadow: inset 0 0 0 1px var(--sel); }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 340px; gap: 22px; align-items: start; }
@media (max-width: 860px) { .layout { grid-template-columns: 1fr; } }
.stage { background: var(--well); border-radius: 6px; overflow: hidden; position: relative; }
canvas { display: block; width: 100%; height: auto; cursor: crosshair; touch-action: none; }
.panel { background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; padding: 14px; }
.panel h2 { font-size: 1.1rem; margin: 0 0 2px; }
.panel .why { color: var(--ink2); margin: 0 0 10px; font-size: .9rem; }
dl { display: grid; grid-template-columns: auto 1fr; gap: 3px 12px; margin: 0 0 12px; font-size: .88rem; }
dt { color: var(--muted); } dd { margin: 0; font-variant-numeric: tabular-nums; }
label { display: block; font-size: .88rem; margin: 8px 0 2px; color: var(--ink2); }
input[type=range] { width: 100%; }
.check { display: flex; align-items: center; gap: 6px; margin: 6px 0; font-size: .88rem; color: var(--ink2); }
.measure { font-variant-numeric: tabular-nums; font-size: .9rem; min-height: 2.6em; margin-top: 8px; color: var(--ink2); }
.measure b { color: var(--ink); }
.claims { font-size: .82rem; color: var(--ink2); margin: 0; padding-left: 18px; }
h3 { font-size: 1rem; margin: 26px 0 6px; }
.bars { width: 100%; }
.bars text { font: 12px system-ui, sans-serif; fill: var(--ink2); }
.bars .sel text { fill: var(--ink); font-weight: 600; }
.note { color: var(--muted); font-size: .8rem; margin-top: 18px; max-width: 90ch; }
</style>
</head>
<body>
<main>
<h1>Giant radio galaxies</h1>
<p class="lede">Jets from supermassive black holes that run millions of light-years beyond their galaxies. Orange: LOFAR's
view at 144 MHz (LoTSS DR2, or DR3 where DR2 has no images); background: the optical sky (DESI Legacy Surveys;
Pan-STARRS1 for NGC 2300). Drag across the image to measure a distance.</p>
<div class="picker" id="picker" role="group" aria-label="Choose a giant"></div>
<div class="layout">
  <div>
    <div class="stage"><canvas id="cv" width="640" height="640" aria-label="Radio emission over the optical sky"></canvas></div>
    <div class="measure" id="measure">Drag across the image to measure a distance.</div>
  </div>
  <div class="panel">
    <h2 id="name"></h2>
    <p class="why" id="why"></p>
    <dl id="facts"></dl>
    <label for="blend">Radio emission <span id="blendv"></span></label>
    <input id="blend" type="range" min="0" max="100" value="85">
    <div class="check"><input type="checkbox" id="others" checked><label for="others" style="margin:0">show unrelated radio sources</label></div>
    <div class="check"><input type="checkbox" id="mw" checked><label for="mw" style="margin:0">show the Milky Way to scale</label></div>
    <p style="font-size:.85rem;margin:12px 0 4px;color:var(--ink2)">Published lengths</p>
    <ul class="claims" id="claims"></ul>
  </div>
</div>
<h3>All of them at one scale</h3>
<svg class="bars" id="bars" role="img" aria-label="Projected lengths of the giants"></svg>
<p class="note">Lengths: our 3-sigma measurements on LoTSS DR2 (Planck 2018 cosmology), projected on the sky, so true 3D
lengths are longer. The Milky Way's disc is ~30 kpc across; Andromeda is 761 kpc away. Data: LoTSS DR2 (Shimwell et al.
2022) and DR3 (Shimwell et al. 2026), DESI Legacy Imaging Surveys DR10 (Dey et al. 2019), Pan-STARRS1 (Chambers et al. 2016). Part of the Astrophysics Project Atlas, project I18 (AI-generated).</p>
</main>
<script>
const P = /*PAYLOAD*/null;
const cv = document.getElementById("cv"), ctx = cv.getContext("2d");
const N = 640;
let cur = 0, drag = null, imgs = {};
function load(src) { return new Promise(r => { const i = new Image(); i.onload = () => r(i); i.src = src; }); }
function tint(img) {
  const c = document.createElement("canvas"); c.width = N; c.height = N;
  const x = c.getContext("2d"); x.drawImage(img, 0, 0, N, N);
  const d = x.getImageData(0, 0, N, N), a = d.data, [r, g, b] = P.radio_rgb;
  for (let i = 0; i < a.length; i += 4) { const v = a[i]; a[i] = r * v / 255; a[i + 1] = g * v / 255; a[i + 2] = b * v / 255; a[i + 3] = 255; }
  x.putImageData(d, 0, 0); return c;
}
async function prep(k) {
  if (imgs[k]) return imgs[k];
  const G = P.giants[k];
  const [o, r, t] = await Promise.all([load(G.optical), load(G.radio), load(G.other)]);
  return (imgs[k] = { o, r: tint(r), t: tint(t) });
}
function fmt(x, d) { return Number(x).toFixed(d); }
async function draw() {
  const G = P.giants[cur], I = await prep(cur);
  const f = document.getElementById("blend").value / 100;
  ctx.globalCompositeOperation = "source-over"; ctx.globalAlpha = 1;
  ctx.fillStyle = "#000"; ctx.fillRect(0, 0, N, N);
  ctx.globalAlpha = 0.85; ctx.drawImage(I.o, 0, 0, N, N);
  ctx.globalCompositeOperation = "screen";
  ctx.globalAlpha = f; ctx.drawImage(I.r, 0, 0);
  if (document.getElementById("others").checked) { ctx.globalAlpha = f * 0.45; ctx.drawImage(I.t, 0, 0); }
  ctx.globalCompositeOperation = "source-over"; ctx.globalAlpha = 1;
  // host
  ctx.strokeStyle = "rgba(255,255,255,.9)"; ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.arc(G.host_px[0], G.host_px[1], 9, 0, 2 * Math.PI); ctx.stroke();
  // scale bar + Milky Way
  const pxPerMpc = 1000 / G.kpc_per_arcmin / G.pix_arcmin;
  let bar = pxPerMpc, lab = "1 Mpc";
  if (bar > 0.55 * N) { bar = pxPerMpc / 2; lab = "0.5 Mpc"; }
  ctx.fillStyle = "#fff"; ctx.font = "13px system-ui, sans-serif";
  ctx.fillRect(24, N - 30, bar, 2); ctx.fillText(lab, 24 + bar / 2 - 20, N - 38);
  if (document.getElementById("mw").checked) {
    const w = Math.max(P.local_group.mw_kpc / 1000 * pxPerMpc, 3);
    ctx.beginPath(); ctx.ellipse(24 + bar + 40, N - 29, w / 2, Math.max(w * 0.15, 1.5), 0, 0, 2 * Math.PI); ctx.fill();
    ctx.fillText("Milky Way", 24 + bar + 50 + w / 2, N - 24);
  }
  if (drag && drag.b) {
    ctx.strokeStyle = "#7fd1ff"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(drag.a[0], drag.a[1]); ctx.lineTo(drag.b[0], drag.b[1]); ctx.stroke();
  }
}
function setGiant(k) {
  cur = k; drag = null;
  const G = P.giants[k];
  document.querySelectorAll("#picker button").forEach((b, i) => b.setAttribute("aria-pressed", i === k));
  document.getElementById("name").textContent = G.name;
  document.getElementById("why").textContent = G.why[0].toUpperCase() + G.why.slice(1) + ".";
  const ly = G.lp * 3.2616;
  document.getElementById("facts").innerHTML =
    `<dt>Length</dt><dd>${fmt(G.lp, 2)} Mpc (${fmt(ly, 1)} million light-years)</dd>` +
    `<dt>On the sky</dt><dd>${fmt(G.las, 1)}' (the full Moon: 31')</dd>` +
    `<dt>Redshift</dt><dd>${G.z} (${G.z_kind === "s" ? "spectroscopic" : "photometric"})</dd>` +
    `<dt>Light left it</dt><dd>${fmt(G.lookback_gyr, 2)} billion years ago</dd>` +
    `<dt>Host galaxy</dt><dd>${G.host}</dd>` +
    `<dt>Milky Ways end to end</dt><dd>${Math.round(G.lp * 1000 / P.local_group.mw_kpc)}</dd>`;
  document.getElementById("claims").innerHTML = G.claims.map(c =>
    `<li>${c.source}: ${c.lp} Mpc${c.las ? ` (${c.las}')` : ""}</li>`).join("");
  document.getElementById("measure").textContent = "Drag across the image to measure a distance.";
  draw(); bars();
}
function pos(e) { const r = cv.getBoundingClientRect(); return [(e.clientX - r.left) * N / r.width, (e.clientY - r.top) * N / r.height]; }
cv.addEventListener("pointerdown", e => { drag = { a: pos(e), b: null }; cv.setPointerCapture(e.pointerId); });
cv.addEventListener("pointermove", e => {
  if (!drag) return; drag.b = pos(e);
  const G = P.giants[cur];
  const d = Math.hypot(drag.b[0] - drag.a[0], drag.b[1] - drag.a[1]) * G.pix_arcmin;
  const mpc = d * G.kpc_per_arcmin / 1000;
  document.getElementById("measure").innerHTML = `<b>${fmt(d, 2)}'</b> on the sky = <b>${fmt(mpc, 2)} Mpc</b> at z = ${G.z} ` +
    `(${fmt(mpc * 3.2616, 1)} million light-years, ${Math.round(mpc * 1000 / P.local_group.mw_kpc)} Milky Ways)`;
  draw();
});
cv.addEventListener("pointerup", () => { if (drag && !drag.b) drag = null; });
document.getElementById("blend").addEventListener("input", e => { document.getElementById("blendv").textContent = e.target.value + "%"; draw(); });
["others", "mw"].forEach(id => document.getElementById(id).addEventListener("change", draw));
function bars() {
  const svg = document.getElementById("bars"), W = svg.clientWidth || 900, rowH = 22;
  const items = P.giants.map((G, i) => ({ name: G.name, v: G.lp, i }));
  items.push({ name: "Milky Way to Andromeda", v: P.local_group.sep_kpc / 1000, i: -1 });
  const max = Math.max(...items.map(d => d.v)), left = 190, right = 70;
  svg.setAttribute("viewBox", `0 0 ${W} ${items.length * rowH + 8}`); svg.setAttribute("height", items.length * rowH + 8);
  const ink = getComputedStyle(document.documentElement).getPropertyValue("--accent").trim();
  const sel = getComputedStyle(document.documentElement).getPropertyValue("--sel").trim();
  const mut = getComputedStyle(document.documentElement).getPropertyValue("--muted").trim();
  svg.innerHTML = items.map((d, k) => {
    const y = k * rowH + 4, w = (W - left - right) * d.v / max, isSel = d.i === cur;
    const col = d.i < 0 ? mut : (isSel ? sel : ink);
    return `<g class="${isSel ? "sel" : ""}" style="cursor:${d.i < 0 ? "default" : "pointer"}" data-i="${d.i}">` +
      `<text x="${left - 8}" y="${y + 14}" text-anchor="end">${d.name}</text>` +
      `<rect x="${left}" y="${y + 4}" width="${Math.max(w, 2)}" height="12" rx="3" fill="${col}"></rect>` +
      `<text x="${left + w + 6}" y="${y + 14}">${fmt(d.v, 2)} Mpc</text></g>`;
  }).join("");
  svg.querySelectorAll("g[data-i]").forEach(g => { const i = +g.dataset.i; if (i >= 0) g.addEventListener("click", () => setGiant(i)); });
}
const picker = document.getElementById("picker");
P.giants.forEach((G, i) => { const b = document.createElement("button"); b.textContent = G.name; b.onclick = () => setGiant(i); picker.appendChild(b); });
document.getElementById("blendv").textContent = "85%";
window.addEventListener("resize", bars);
setGiant(0);
</script>
</body>
</html>
"""
