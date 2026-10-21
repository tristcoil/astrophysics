"""A self-contained Einstein-ring explorer.

    figures/lens_explorer.html

Every lens we modelled, one at a time. Left: the Hubble image with the lens
galaxy's light removed, leaving the ring. Middle: what a lens model predicts.
Right: the source plane, where you can drag the background galaxy around behind
the lens and watch one image become two, then four, then a ring as it crosses
the caustics. Sliders set the lens's Einstein radius, flattening and
orientation; the page converts the Einstein radius into a mass and a velocity
dispersion as you go, and scores how well your model matches the real ring.

Buttons load our best fit, the published SLACS model, or the ring an
isothermal galaxy would make given the stars' measured speed.

Everything (data and code) is in one file, so it works offline. It follows the
reader's light/dark setting, and a table repeats the selected lens's numbers.
"""

from __future__ import annotations

import base64
import json

import numpy as np
from matplotlib import colormaps

import hstimages as H
import lensing as L
import plotstyle as ps
from config import FIGURE_DIR, MODEL_DIR, SHOWCASE

THEMES = {
    "light": {"surface": ps.SURFACE, "ink": ps.INK, "ink2": ps.INK_SECONDARY, "muted": ps.MUTED,
              "grid": ps.GRID, "axis": ps.AXIS, "ours": ps.ORANGE, "published": ps.VIOLET,
              "panel": "#0b0b0b"},
    "dark": {"surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#4a4a46", "ours": "#d95926", "published": "#9085e9",
             "panel": "#000000"},
}


def _b64(arr: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode("ascii")


def _record(row, npz, res) -> dict:
    from lensmodel import e_to_qphi

    data, var, fm = npz["data"], npz["var"], npz["fit_mask"]
    # Inside the fitted region, the joint fit's lens light; outside it, the
    # whole-box light model (the joint fit is not constrained out there).
    light = np.where(fm, npz["light"], npz["light_full"]) if "light_full" in npz else npz["light"]
    sig = float(np.sqrt(np.median(var[np.isfinite(var)])))
    resid = (data - light) / sig  # in units of the background noise
    # Beyond the fitted circle the light model suffers from the box edge (the
    # blurring cannot see light outside the cutout): show nothing there.
    n = data.shape[0]
    yy, xx = np.mgrid[:n, :n] - n // 2
    rr = np.hypot(xx, yy) * 0.04
    resid[(rr > 2.9) | (rr < 0.25)] = 0.0  # also the masked core, where no light model applies
    q8 = np.clip(np.round((resid + 4.0) * 4.0), 0, 255).astype(np.uint8)  # 0.25-sigma steps
    geo = H.model_geometry(npz)
    p = npz["params"]
    q_s, phi_s = e_to_qphi(p[7], p[8])
    north, parity = float(npz["north_angle"]), int(npz["parity"])
    rec = {
        "name": row.name, "n": int(data.shape[0]), "pix": 0.04,
        "img": _b64(q8), "mask": _b64(np.packbits(fm.ravel())),
        "north": north, "parity": parity,
        "fit": {"b": geo["b"], "q": geo["q"], "phi": geo["phi"], "x0": geo["x0"], "y0": geo["y0"],
                "xs": geo["xs"], "ys": geo["ys"], "rs": geo["r_s"], "ns": geo["n_s"],
                "qs": float(q_s), "phis": float(phi_s)},
        "zl": float(row.z_lens), "zs": float(row.z_source),
        "kpc": float(row.kpc_per_arcsec), "scr": float(row.sigma_cr),
        "dsdls": float(row.d_s / row.d_ls),
        "sigma": None if not np.isfinite(row.sigma_star) else float(row.sigma_star),
        "esigma": None if not np.isfinite(row.e_sigma_star) else float(row.e_sigma_star),
        "pub": None, "showcase": row.name in SHOWCASE,
        "chi2": float(res.get("chi2_red", np.nan)),
    }
    if np.isfinite(row.b_sie):
        rec["pub"] = {"b": float(row.b_sie), "q": float(row.q_sie),
                      "phi": north + parity * np.radians(float(row.pa_sie)), "src": "SLACS V"}
    elif np.isfinite(row.theta_e):
        rec["pub"] = {"b": float(row.theta_e), "q": None, "phi": None, "src": "SLACS IX"}
    return rec


def build(t, ours) -> list[str]:
    records = []
    table = ours.get("table")
    good = (set(table.loc[(table["failed"] == "") & (table["flag"] == ""), "name"])
            if table is not None and len(table) else set())
    for row in t.itertuples():
        if row.name not in good:
            continue
        npz, res = H.load_model(row.name)
        if npz is None:
            continue
        records.append(_record(row, npz, res))
    records.sort(key=lambda r: (not r["showcase"], r["name"]))
    lut = (np.array(colormaps["magma"](np.linspace(0, 1, 256)))[:, :3] * 255).astype(int).tolist()
    html = _page(records, lut)
    out = FIGURE_DIR / "lens_explorer.html"
    out.write_text(html)
    return [str(out)]


def _page(records, lut) -> str:
    Lt, D = THEMES["light"], THEMES["dark"]
    css_vars = lambda T: "; ".join(f"--{k}: {v}" for k, v in T.items())
    data_js = json.dumps(records, separators=(",", ":"))
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Einstein Ring Explorer</title>
<style>
  :root {{ color-scheme: light; {css_vars(Lt)}; }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{ color-scheme: dark; {css_vars(D)}; }}
  }}
  :root[data-theme="dark"] {{ color-scheme: dark; {css_vars(D)}; }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1180px; margin: 0 auto; padding: 18px 16px 40px; }}
  h1 {{ font-size: 22px; margin: 4px 0 2px; }}
  .lede {{ color: var(--ink2); margin: 0 0 14px; max-width: 900px; }}
  .bar {{ display: flex; flex-wrap: wrap; gap: 10px 14px; align-items: center; margin: 10px 0 14px; }}
  select, button {{ font: inherit; color: var(--ink); background: var(--surface);
         border: 1px solid var(--axis); border-radius: 6px; padding: 5px 10px; cursor: pointer; }}
  button:hover, select:hover {{ border-color: var(--ink2); }}
  .panels {{ display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 14px; }}
  .panel h2 {{ font-size: 14px; margin: 0 0 2px; }}
  .panel p {{ font-size: 12.5px; color: var(--ink2); margin: 0 0 6px; min-height: 2.6em; }}
  canvas {{ width: 100%; aspect-ratio: 1; background: var(--panel); border-radius: 6px;
           display: block; touch-action: none; }}
  #src {{ cursor: grab; }}
  .sliders {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 10px 18px;
             margin: 16px 0 8px; }}
  .sliders label {{ font-size: 13px; color: var(--ink2); display: block; }}
  .sliders output {{ color: var(--ink); font-variant-numeric: tabular-nums; }}
  input[type=range] {{ width: 100%; accent-color: var(--ours); }}
  .readout {{ display: grid; grid-template-columns: repeat(5, minmax(0, 1fr)); gap: 10px;
             margin: 12px 0; }}
  .tile {{ border: 1px solid var(--grid); border-radius: 8px; padding: 8px 10px; }}
  .tile .k {{ font-size: 12px; color: var(--ink2); }}
  .tile .v {{ font-size: 19px; font-variant-numeric: tabular-nums; }}
  .tile .s {{ font-size: 11.5px; color: var(--muted); }}
  table {{ border-collapse: collapse; font-size: 13px; margin-top: 8px; }}
  th, td {{ text-align: left; padding: 3px 12px 3px 0; border-bottom: 1px solid var(--grid); }}
  th {{ color: var(--ink2); font-weight: 600; }}
  .key {{ font-size: 12.5px; color: var(--ink2); }}
  .swatch {{ display: inline-block; width: 18px; height: 0; border-top: 2px dashed; vertical-align: middle;
             margin: 0 4px 0 10px; }}
  .note {{ font-size: 12.5px; color: var(--muted); margin-top: 18px; max-width: 900px; }}
  @media (max-width: 860px) {{
    .panels {{ grid-template-columns: 1fr; }}
    .sliders {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
    .readout {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
  }}
</style>
</head>
<body>
<main>
  <h1>Einstein ring explorer</h1>
  <p class="lede">A massive galaxy bends the light of a galaxy far behind it into arcs, or a ring.
  The ring's size weighs the lens. Pick a lens, then drag the background galaxy around in the
  right-hand panel and set the lens's mass and shape with the sliders. The middle panel shows what that
  lens would do; the left panel is what Hubble actually saw.</p>
  <div class="bar">
    <label>Lens <select id="pick" aria-label="Choose a lens"></select></label>
    <button id="ours">Our best fit</button>
    <button id="pub">Published model</button>
    <button id="iso">Isothermal guess from the stars' speed</button>
    <label class="key"><input type="checkbox" id="showpub" checked> show published critical curve</label>
  </div>
  <div class="panels">
    <div class="panel"><h2>What Hubble saw</h2>
      <p>F814W image minus our model of the lens galaxy's light. Dashed: your lens model's
      critical curve (white) and the published one (violet).</p>
      <canvas id="data" width="320" height="320" aria-label="Hubble image of the ring"></canvas></div>
    <div class="panel"><h2>What your lens predicts</h2>
      <p>The background galaxy, seen through the lens (same scale, 6.4&Prime; across).</p>
      <canvas id="model" width="320" height="320" aria-label="Predicted lensed image"></canvas></div>
    <div class="panel"><h2>Behind the lens: the source plane</h2>
      <p>Drag the galaxy (or use arrow keys). <span style="color:var(--ours)">Solid: caustic; dashed: cut.</span>
      Crossing them adds images.</p>
      <canvas id="src" width="320" height="320" tabindex="0" aria-label="Source plane: drag to move the background galaxy"></canvas></div>
  </div>
  <div class="sliders">
    <label>Einstein radius <output id="bv"></output>
      <input type="range" id="b" min="0.4" max="2.2" step="0.005"></label>
    <label>Axis ratio of the mass <output id="qv"></output>
      <input type="range" id="q" min="0.35" max="1" step="0.01"></label>
    <label>Orientation (° east of north) <output id="pav"></output>
      <input type="range" id="pa" min="0" max="180" step="1"></label>
    <label>Size of the background galaxy <output id="rsv"></output>
      <input type="range" id="rs" min="0.02" max="0.5" step="0.005"></label>
  </div>
  <div class="readout">
    <div class="tile"><div class="k">Mass inside the ring</div><div class="v" id="mass"></div><div class="s">from lensing alone</div></div>
    <div class="tile"><div class="k">Ring radius</div><div class="v" id="rkpc"></div><div class="s" id="rkly"></div></div>
    <div class="tile"><div class="k">Isothermal speed σ_SIE</div><div class="v" id="ssie"></div><div class="s" id="sstar"></div></div>
    <div class="tile"><div class="k">Images of the galaxy's centre</div><div class="v" id="nimg"></div><div class="s">1, 2, 3 or 4</div></div>
    <div class="tile"><div class="k">Match to Hubble's image</div><div class="v" id="match"></div><div class="s" id="matchbest"></div></div>
  </div>
  <table aria-label="Numbers for the selected lens"><thead><tr><th></th><th>Einstein radius</th><th>axis ratio</th><th>orientation</th><th>mass inside ring</th></tr></thead>
  <tbody id="tbl"></tbody></table>
  <p class="note">Data: Hubble ACS/WFC F814W images (MAST Hubble Advanced Products) of the SLACS lenses;
  redshifts and velocity dispersions from Bolton et al. (2008) and Auger et al. (2009). Lens models:
  singular isothermal ellipsoid + Sérsic source, fitted in this project. The match score compares the
  predicted image (scaled to its best brightness) with the data inside the fitted region; it ignores the
  telescope's blurring beyond a small smoothing, so 100% is out of reach. Distances assume
  H<sub>0</sub> = 70 km/s/Mpc, Ω<sub>m</sub> = 0.3.</p>
</main>
<script>
const LENSES = {data_js};
const LUT = {json.dumps(lut)};
const C_KMS = 299792.458, ARCSEC = Math.PI / 648000;
let lens = null, P = null, data = null, mask = null, bestMatch = null;

function decode(s) {{ const b = atob(s); const a = new Uint8Array(b.length);
  for (let i = 0; i < b.length; i++) a[i] = b.charCodeAt(i); return a; }}
function css(name) {{ return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }}

// Singular isothermal ellipsoid, intermediate-axis normalisation (Kormann et al. 1994)
function defl(x, y, p) {{
  const c = Math.cos(p.phi), s = Math.sin(p.phi);
  const dx = x - p.x0, dy = y - p.y0;
  const xr = c * dx + s * dy, yr = -s * dx + c * dy;
  const q = Math.min(p.q, 0.99999), qp = Math.sqrt(1 - q * q);
  const psi = Math.max(Math.sqrt(q * q * xr * xr + yr * yr), 1e-12);
  const sc = p.b * Math.sqrt(q) / qp;
  const ax = sc * Math.atan(qp * xr / psi);
  const t = Math.max(-1 + 1e-12, Math.min(1 - 1e-12, qp * yr / psi));
  const ay = sc * 0.5 * Math.log((1 + t) / (1 - t));
  return [c * ax - s * ay, s * ax + c * ay];
}}
function bn(n) {{ return 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n * n); }}
function sersic(bx, by, p) {{
  const c = Math.cos(p.phis), s = Math.sin(p.phis);
  const dx = bx - p.xs, dy = by - p.ys;
  const xr = c * dx + s * dy, yr = -s * dx + c * dy;
  const r = Math.sqrt(p.qs * xr * xr + yr * yr / p.qs);
  return Math.exp(-bn(p.ns) * (Math.pow(r / p.rs, 1 / p.ns) - 1));
}}
function coord(i, n) {{ return (i - Math.floor(n / 2)) * lens.pix; }}

function critCurve(p, npts) {{
  const out = []; const c = Math.cos(p.phi), s = Math.sin(p.phi);
  for (let k = 0; k <= npts; k++) {{ const t = 2 * Math.PI * k / npts;
    const xr = p.b / Math.sqrt(p.q) * Math.cos(t), yr = p.b * Math.sqrt(p.q) * Math.sin(t);
    out.push([p.x0 + c * xr - s * yr, p.y0 + s * xr + c * yr]); }}
  return out;
}}
function caustics(p) {{
  const cc = critCurve(p, 720).map(([x, y]) => {{ const a = defl(x, y, p); return [x - a[0], y - a[1]]; }});
  const cut = [];
  for (let k = 0; k <= 720; k++) {{ const t = 2 * Math.PI * k / 720;
    const x = p.x0 + 1e-7 * Math.cos(t), y = p.y0 + 1e-7 * Math.sin(t);
    const a = defl(x, y, p); cut.push([x - a[0], y - a[1]]); }}
  return [cc, cut];
}}
function inside(pt, poly) {{
  let c = false; for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {{
    const [xi, yi] = poly[i], [xj, yj] = poly[j];
    if (((yi > pt[1]) !== (yj > pt[1])) && (pt[0] < (xj - xi) * (pt[1] - yi) / (yj - yi) + xi)) c = !c; }}
  return c;
}}

function paint(ctx, vals, n, lo, hi, stretch) {{
  const img = ctx.createImageData(n, n);
  const f = v => Math.asinh(v / stretch);
  const a = f(lo), b = f(hi);
  for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) {{
    const v = vals[(n - 1 - r) * n + c];
    let u = Number.isFinite(v) ? (f(v) - a) / (b - a) : 0;
    u = Math.max(0, Math.min(1, u));
    const L = LUT[Math.round(u * 255)], o = 4 * (r * n + c);
    img.data[o] = L[0]; img.data[o + 1] = L[1]; img.data[o + 2] = L[2]; img.data[o + 3] = 255; }}
  return img;
}}
function drawImage(canvas, vals, n, hi) {{
  const off = document.createElement("canvas"); off.width = n; off.height = n;
  off.getContext("2d").putImageData(paint(off.getContext("2d"), vals, n, 0, hi, 3), 0, 0);
  const ctx = canvas.getContext("2d");
  ctx.imageSmoothingEnabled = true; ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.drawImage(off, 0, 0, canvas.width, canvas.height);
  return ctx;
}}
function toPx(canvas, x, y, half) {{
  return [(x + half) / (2 * half) * canvas.width, (half - y) / (2 * half) * canvas.height];
}}
function polyline(ctx, canvas, pts, half, color, dash, lw) {{
  ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = lw; ctx.setLineDash(dash); ctx.beginPath();
  pts.forEach(([x, y], k) => {{ const [u, v] = toPx(canvas, x, y, half); k ? ctx.lineTo(u, v) : ctx.moveTo(u, v); }});
  ctx.stroke(); ctx.restore();
}}

function modelImage(p) {{
  const n = lens.n, out = new Float32Array(n * n), ss = [-0.25, 0.25];
  for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) {{
    let acc = 0;
    for (const oy of ss) for (const ox of ss) {{
      const x = coord(c + ox, n), y = coord(r + oy, n);
      const a = defl(x, y, p); acc += sersic(x - a[0], y - a[1], p); }}
    out[r * n + c] = acc / 4; }}
  // light smoothing, standing in for Hubble's blur (FWHM ~2.5 pixels)
  const k = [0.25, 0.5, 0.25], tmp = new Float32Array(n * n), sm = new Float32Array(n * n);
  for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) {{ let s = 0;
    for (let d = -1; d <= 1; d++) s += k[d + 1] * out[r * n + Math.min(n - 1, Math.max(0, c + d))]; tmp[r * n + c] = s; }}
  for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) {{ let s = 0;
    for (let d = -1; d <= 1; d++) s += k[d + 1] * tmp[Math.min(n - 1, Math.max(0, r + d)) * n + c]; sm[r * n + c] = s; }}
  return sm;
}}
function match(m) {{
  let dm = 0, mm = 0, dd = 0;
  for (let i = 0; i < m.length; i++) if (mask[i]) {{ dm += data[i] * m[i]; mm += m[i] * m[i]; dd += data[i] * data[i]; }}
  const a = mm > 0 ? Math.max(dm / mm, 0) : 0;
  let rr = 0, n = 0; for (let i = 0; i < m.length; i++) if (mask[i]) {{ const r = data[i] - a * m[i]; rr += r * r; n++; }}
  // Share of the ring's signal (data power above the noise, which is 1 per pixel) that the model explains.
  return [Math.max(0, Math.min(1, 1 - (rr - n) / Math.max(dd - n, 1e-9))), a];
}}

function render() {{
  const n = lens.n, half = n / 2 * lens.pix, hi = lens.vmax;
  // data
  const dc = document.getElementById("data");
  const ctx = drawImage(dc, data, n, hi);
  polyline(ctx, dc, critCurve(P, 360), half, "#ffffff", [5, 4], 1.4);
  if (document.getElementById("showpub").checked && lens.pub && lens.pub.q !== null)
    polyline(ctx, dc, critCurve({{...lens.pub, x0: lens.fit.x0, y0: lens.fit.y0}}, 360), half, css("--published"), [3, 3], 1.4);
  else if (document.getElementById("showpub").checked && lens.pub) {{
    const circ = []; for (let k = 0; k <= 360; k++) {{ const t = 2 * Math.PI * k / 360;
      circ.push([lens.fit.x0 + lens.pub.b * Math.cos(t), lens.fit.y0 + lens.pub.b * Math.sin(t)]); }}
    polyline(ctx, dc, circ, half, css("--published"), [3, 3], 1.4); }}
  // model
  const m = modelImage(P);
  const [score, amp] = match(m);
  const scaled = m.map(v => v * amp);
  drawImage(document.getElementById("model"), scaled, n, hi);
  // source plane
  const sc = document.getElementById("src"), sctx = sc.getContext("2d"), sh = Math.max(0.6, 1.15 * P.b);
  const sn = 128, sv = new Float32Array(sn * sn);
  for (let r = 0; r < sn; r++) for (let c = 0; c < sn; c++)
    sv[r * sn + c] = sersic(P.x0 + (c - sn / 2 + 0.5) / sn * 2 * sh, P.y0 + (r - sn / 2 + 0.5) / sn * 2 * sh, P);
  const off = document.createElement("canvas"); off.width = sn; off.height = sn;
  off.getContext("2d").putImageData(paint(off.getContext("2d"), sv, sn, 0, 1.2, 0.3), 0, 0);
  sctx.imageSmoothingEnabled = true; sctx.drawImage(off, 0, 0, sc.width, sc.height);
  const [caus, cut] = caustics(P);
  const shift = pts => pts.map(([x, y]) => [x - P.x0, y - P.y0]);
  polyline(sctx, sc, shift(caus), sh, css("--ours"), [], 1.6);
  polyline(sctx, sc, shift(cut), sh, css("--ours"), [5, 4], 1.4);
  const [u, v] = toPx(sc, P.xs - P.x0, P.ys - P.y0, sh);
  sctx.strokeStyle = "#ffffff"; sctx.lineWidth = 1.2; sctx.beginPath();
  sctx.moveTo(u - 6, v); sctx.lineTo(u + 6, v); sctx.moveTo(u, v - 6); sctx.lineTo(u, v + 6); sctx.stroke();
  sctx.fillStyle = "#ffffff"; sctx.font = "11px system-ui"; sctx.fillText(`${{(2 * sh).toFixed(1)}}″ across`, 8, sc.height - 8);
  // images of the source centre
  const pt = [P.xs, P.ys], inC = inside(pt, caus), inCut = inside(pt, cut);
  const nimg = inC && inCut ? 4 : (inCut ? 2 : (inC ? 3 : 1));
  // readouts
  const rk = P.b * lens.kpc, mass = Math.PI * lens.scr * rk * rk;
  const ssie = C_KMS * Math.sqrt(P.b * ARCSEC * lens.dsdls / (4 * Math.PI));
  document.getElementById("mass").textContent = (mass / 1e11).toFixed(2) + " × 10¹¹ M☉";
  document.getElementById("rkpc").textContent = rk.toFixed(2) + " kpc";
  document.getElementById("rkly").textContent = (rk * 3.2616).toFixed(1) + " thousand light-years";
  document.getElementById("ssie").textContent = ssie.toFixed(0) + " km/s";
  document.getElementById("sstar").textContent = lens.sigma ? `stars: ${{lens.sigma}} ± ${{lens.esigma}} km/s` : "stars: not measured";
  document.getElementById("nimg").textContent = nimg;
  document.getElementById("match").textContent = (100 * score).toFixed(0) + "%";
  document.getElementById("matchbest").textContent = `our fit: ${{(100 * bestMatch).toFixed(0)}}%`;
  document.getElementById("bv").textContent = P.b.toFixed(3) + "″";
  document.getElementById("qv").textContent = P.q.toFixed(2);
  document.getElementById("pav").textContent = paOf(P.phi).toFixed(0) + "°";
  document.getElementById("rsv").textContent = P.rs.toFixed(3) + "″";
}}

function paOf(phi) {{ let pa = (phi - lens.north) * lens.parity * 180 / Math.PI; return ((pa % 180) + 180) % 180; }}
function phiOf(pa) {{ return lens.north + lens.parity * pa * Math.PI / 180; }}
function syncSliders() {{
  document.getElementById("b").value = P.b; document.getElementById("q").value = P.q;
  document.getElementById("pa").value = paOf(P.phi); document.getElementById("rs").value = P.rs;
}}
function table() {{
  const rows = [["Our fit", lens.fit.b, lens.fit.q, paOf(lens.fit.phi)]];
  if (lens.pub) rows.push([`Published (${{lens.pub.src}})`, lens.pub.b, lens.pub.q, lens.pub.phi === null ? null : paOf(lens.pub.phi)]);
  document.getElementById("tbl").innerHTML = rows.map(([lab, b, q, pa]) => {{
    const rk = b * lens.kpc, m = Math.PI * lens.scr * rk * rk;
    return `<tr><td>${{lab}}</td><td>${{b.toFixed(3)}}″ (${{rk.toFixed(2)}} kpc)</td><td>${{q === null ? "—" : q.toFixed(2)}}</td>` +
           `<td>${{pa === null ? "—" : pa.toFixed(0) + "°"}}</td><td>${{(m / 1e11).toFixed(2)}} × 10¹¹ M☉</td></tr>`; }}).join("");
}}
function load(name) {{
  lens = LENSES.find(l => l.name === name);
  const raw = decode(lens.img), bits = decode(lens.mask), n = lens.n;
  data = new Float32Array(n * n); mask = new Uint8Array(n * n);
  for (let i = 0; i < n * n; i++) {{ data[i] = raw[i] / 4 - 4; mask[i] = (bits[i >> 3] >> (7 - (i & 7))) & 1; }}
  const vals = []; for (let i = 0; i < n * n; i++) if (mask[i]) vals.push(data[i]);
  vals.sort((a, b) => a - b); lens.vmax = Math.max(vals[Math.floor(0.995 * (vals.length - 1))], 3);
  P = {{...lens.fit}};
  bestMatch = match(modelImage(P))[0];
  syncSliders(); table(); render();
}}

const pick = document.getElementById("pick");
LENSES.forEach(l => {{ const o = document.createElement("option"); o.value = l.name;
  o.textContent = l.name + (l.showcase ? " ★" : ""); pick.appendChild(o); }});
pick.onchange = () => {{ load(pick.value); history.replaceState(null, "", "#" + pick.value); }};
for (const id of ["b", "q", "pa", "rs"]) document.getElementById(id).oninput = e => {{
  const v = parseFloat(e.target.value);
  if (id === "pa") P.phi = phiOf(v); else P[id] = v; render(); }};
document.getElementById("showpub").onchange = render;
document.getElementById("ours").onclick = () => {{ P = {{...lens.fit}}; syncSliders(); render(); }};
document.getElementById("pub").onclick = () => {{ if (!lens.pub) return;
  P = {{...P, b: lens.pub.b}}; if (lens.pub.q !== null) {{ P.q = lens.pub.q; P.phi = lens.pub.phi; }} syncSliders(); render(); }};
document.getElementById("iso").onclick = () => {{ if (!lens.sigma) return;
  P = {{...P, b: 4 * Math.PI * Math.pow(lens.sigma / C_KMS, 2) / lens.dsdls / ARCSEC}}; syncSliders(); render(); }};
const sc = document.getElementById("src");
function drag(e) {{
  const rect = sc.getBoundingClientRect(), sh = Math.max(0.6, 1.15 * P.b);
  P.xs = P.x0 + ((e.clientX - rect.left) / rect.width * 2 - 1) * sh;
  P.ys = P.y0 + (1 - (e.clientY - rect.top) / rect.height * 2) * sh; render(); }}
sc.addEventListener("pointerdown", e => {{ sc.setPointerCapture(e.pointerId); sc.style.cursor = "grabbing"; drag(e); }});
sc.addEventListener("pointermove", e => {{ if (sc.hasPointerCapture(e.pointerId)) drag(e); }});
sc.addEventListener("pointerup", e => {{ sc.releasePointerCapture(e.pointerId); sc.style.cursor = "grab"; }});
sc.addEventListener("keydown", e => {{
  const d = 0.01 * (e.shiftKey ? 5 : 1), m = {{ArrowLeft: [-d, 0], ArrowRight: [d, 0], ArrowUp: [0, d], ArrowDown: [0, -d]}}[e.key];
  if (m) {{ e.preventDefault(); P.xs += m[0]; P.ys += m[1]; render(); }} }});
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
const start = decodeURIComponent(location.hash.slice(1));
const first = LENSES.find(l => l.name === start) ? start : LENSES[0].name;
pick.value = first; load(first);
</script>
</body>
</html>
"""
