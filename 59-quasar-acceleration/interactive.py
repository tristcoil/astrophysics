"""Offline explorer: refit the quasar drift yourself (figures/quasar_drift_explorer.html).

The page embeds Gaia's quasar proper motions summed per HEALPix pixel
(nside 16, ~3.7 deg), separately for 5 magnitude bins and 4 colour bins and
for three random tranches of the quasars (1%, 9% and 90%), so it can redo the
vector-spherical-harmonic fit in the browser for any combination:

* per pixel p: W_p = sum_k E_p^T (E_k C_k^-1 E_k^T) E_p (2 x 2) and
  b_p = sum_k E_p^T E_k C_k^-1 mu_k (2), in the galactic (east, north) frame
  at the pixel centre, i.e. each quasar's proper motion and covariance moved to
  the pixel centre through 3D;
* the real VSH basis up to l = 3 at every pixel centre, in galactic axes, from
  vsh.py (the same code as the main fit).

The baseline fit's outliers are left out. build() checks that the binned fit
reproduces the exact fit and records the difference in the page.
"""

from __future__ import annotations

import base64
import json
import math

import healpy as hp
import numpy as np
import pandas as pd

import analysis
import config as C
import vsh

NSIDE = 16
NSIDE_ARROWS = 4
G_EDGES = [(-np.inf, 18.0), (18.0, 19.0), (19.0, 20.0), (20.0, 20.5), (20.5, np.inf)]
G_LABELS = ["G ≤ 18", "18–19", "19–20", "20–20.5", "> 20.5"]
TRANCHES = [0.01, 0.09, 0.90]
A = np.array(C.A_G)


def _frames(l_deg, b_deg):
    l, b = np.radians(l_deg), np.radians(b_deg)
    u = np.column_stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])
    e = np.column_stack([-np.sin(l), np.cos(l), np.zeros_like(l)])
    n = np.column_stack([-np.sin(b) * np.cos(l), -np.sin(b) * np.sin(l), np.cos(b)])
    return u, e, n


def pixel_sums(df: pd.DataFrame, mask: np.ndarray) -> np.ndarray:
    """(npix, 6): W_ee, W_en, W_nn, b_e, b_n, count for the quasars in mask."""
    npix = hp.nside2npix(NSIDE)
    d = df[mask]
    s = analysis.make_sample(df, mask)
    w_aa, w_ad, w_dd = s.weights()
    u_icrs = analysis.unit_vectors(d["ra"].to_numpy(), d["dec"].to_numpy())
    e_ra, e_dec = analysis.tangent_basis(u_icrs)
    eg_ra, eg_dec = e_ra @ A.T, e_dec @ A.T
    pix = hp.ang2pix(NSIDE, d["l"].to_numpy(), d["b"].to_numpy(), lonlat=True)
    pl, pb = hp.pix2ang(NSIDE, np.arange(npix), lonlat=True)
    _, pe, pn = _frames(pl, pb)
    # quasar frame vectors projected on its pixel centre's (east, north)
    Pa = np.column_stack([np.einsum("ij,ij->i", eg_ra, pe[pix]), np.einsum("ij,ij->i", eg_ra, pn[pix])])
    Pd = np.column_stack([np.einsum("ij,ij->i", eg_dec, pe[pix]), np.einsum("ij,ij->i", eg_dec, pn[pix])])
    # W_c = Pa^T w_aa Pa + Pa^T w_ad Pd + Pd^T w_ad Pa + Pd^T w_dd Pd  (per quasar, 2x2)
    def outer(x, y, w):
        return np.einsum("ni,nj,n->nij", x, y, w)
    Wq = outer(Pa, Pa, w_aa) + outer(Pa, Pd, w_ad) + outer(Pd, Pa, w_ad) + outer(Pd, Pd, w_dd)
    bq = Pa * (w_aa * s.pmra + w_ad * s.pmdec)[:, None] + Pd * (w_ad * s.pmra + w_dd * s.pmdec)[:, None]
    out = np.zeros((npix, 6))
    for j, vals in enumerate((Wq[:, 0, 0], Wq[:, 0, 1], Wq[:, 1, 1], bq[:, 0], bq[:, 1], np.ones(len(d)))):
        out[:, j] = np.bincount(pix, weights=vals, minlength=npix)
    return out


def binned_fit(sums: np.ndarray, basis_e: np.ndarray, basis_n: np.ndarray, p: int):
    """The same normal equations the page solves (for the self-check)."""
    Be, Bn = basis_e[:, :p], basis_n[:, :p]
    W = sums
    N = (Be.T * W[:, 0]) @ Be + (Be.T * W[:, 1]) @ Bn + (Bn.T * W[:, 1]) @ Be + (Bn.T * W[:, 2]) @ Bn
    rhs = Be.T @ W[:, 3] + Bn.T @ W[:, 4]
    return np.linalg.solve(N, rhs), np.linalg.inv(N)


def _b64(arr: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(arr, dtype="<f4").tobytes()).decode()


def build(df: pd.DataFrame, results: dict, outlier: np.ndarray) -> dict:
    five = df["five_param"].to_numpy() & ~np.asarray(outlier)
    rng = np.random.default_rng(C.SEED + 99)
    tranche = np.searchsorted(np.cumsum(TRANCHES), rng.random(len(df)), side="right")
    tranche = np.minimum(tranche, len(TRANCHES) - 1)
    g = df["phot_g_mean_mag"].to_numpy()
    nu = df["nu_eff_used_in_astrometry"].to_numpy()
    nu_q = np.nanquantile(nu[five], [0, 0.25, 0.5, 0.75, 1.0])
    colour_bins = [(nu_q[i], nu_q[i + 1] + (1e-9 if i == 3 else 0)) for i in range(4)]
    colour_labels = ["reddest quarter", "red-middle", "blue-middle", "bluest quarter"]
    sets = []
    for axis, edges in (("mag", G_EDGES), ("colour", colour_bins)):
        col = g if axis == "mag" else nu
        for i, (lo, hi) in enumerate(edges):
            for t in range(len(TRANCHES)):
                m = five & (col > lo if axis == "mag" else col >= lo) & (col <= hi if axis == "mag" else col < hi) \
                    & (tranche == t)
                sets.append(pixel_sums(df, m))
    sums = np.stack(sets)  # (nsets, npix, 6)

    npix = hp.nside2npix(NSIDE)
    pl, pb = hp.pix2ang(NSIDE, np.arange(npix), lonlat=True)
    tl = vsh.terms(3)
    be, bn = vsh.basis(np.radians(pl), np.radians(pb), 3, tl)
    m_glide = vsh.l1_vector_matrix("S", tl)
    m_rot = vsh.l1_vector_matrix("T", tl)

    # self-check against the exact baseline fit (all magnitudes, all tranches)
    total = sums[: len(G_EDGES) * len(TRANCHES)].sum(axis=0)
    coef, cov = binned_fit(total, be, bn, vsh.n_params(3))
    g_binned = m_glide @ coef
    g_exact = np.array(results["baseline"]["scaled"]["glide_galactic"]["vec"])
    check = {"binned_glide_gal": g_binned.tolist(), "exact_glide_gal": g_exact.tolist(),
             "max_abs_diff_uas_yr": float(np.max(np.abs(g_binned - g_exact)))}

    al, ab = hp.pix2ang(NSIDE_ARROWS, np.arange(hp.nside2npix(NSIDE_ARROWS)), lonlat=True)
    base = results["baseline"]
    pred = results["predictions"]["constants"]
    payload = {
        "nside": NSIDE, "pix_l": np.round(pl, 4).tolist(), "pix_b": np.round(pb, 4).tolist(),
        "arrow_l": np.round(al, 3).tolist(), "arrow_b": np.round(ab, 3).tolist(),
        "n_terms": len(tl), "n_params": {str(l): vsh.n_params(l) for l in (1, 2, 3)},
        "basis_e": _b64(be), "basis_n": _b64(bn),
        "m_glide": m_glide.tolist(), "m_rot": m_rot.tolist(),
        "glide_idx": [i for i, t in enumerate(tl) if t.l == 1 and t.kind == "S"],
        "sums": _b64(sums.reshape(-1)), "n_sets": int(sums.shape[0]),
        "axes": {"mag": {"labels": G_LABELS, "first": 0},
                 "colour": {"labels": colour_labels, "first": len(G_EDGES) * len(TRANCHES),
                            "edges": [list(map(float, e)) for e in colour_bins]}},
        "n_tranches": len(TRANCHES), "tranche_fracs": TRANCHES,
        "uwe": base["scaled"]["uwe"],
        "sgra_lb": list(_sgra_lb()),
        "uas_to_m_s2": C.UAS_YR_TO_M_S2,
        "kms2_kpc_to_uas": C.KMS2_PER_KPC_TO_M_S2 / C.UAS_YR_TO_M_S2,
        "r0_h1": C.R0_H1_KPC, "r0_h1_err": C.R0_H1_ERR_KPC,
        "reference": {
            "ours": {"amp": base["bootstrap"]["glide_icrs"]["amplitude"]["value"],
                     "amp_err": base["bootstrap"]["glide_icrs"]["amplitude"]["err"],
                     "gX": base["bootstrap"]["glide_galactic"]["vec"][0],
                     "gX_err": base["bootstrap"]["glide_galactic"]["err"][0]},
            "klioner": {"amp": C.REFERENCE["klioner2021"]["amplitude_uas_yr"][0],
                        "amp_err": C.REFERENCE["klioner2021"]["amplitude_uas_yr"][1]},
            "pred_h1": {"g": pred["h1_sgra"]["g_uas_yr"]["median"], "sd": pred["h1_sgra"]["g_uas_yr"]["sd"]},
        },
        "check": check,
    }
    html = TEMPLATE.replace("__PAYLOAD__", json.dumps(payload, separators=(",", ":")))
    path = C.FIGURE_DIR / "quasar_drift_explorer.html"
    path.write_text(html, encoding="utf-8")
    return {"path": str(path), "bytes": path.stat().st_size, **check}


def _sgra_lb():
    ra, dec = math.radians(C.SGRA_RA_DEG), math.radians(C.SGRA_DEC_DEG)
    v = A @ np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(v[2]))


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Quasar Drift Explorer</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --panel: #f4f3ef; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --data: #0b0b0b; --model: #eb6834; --ours: #2a78d6;
  --star: #eda100; --chip: #ffffff; --chip-on: #2a78d6; --chip-on-ink: #ffffff; --wash: #dfe9f7;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --surface: #1a1a19; --panel: #232321; --ink: #f0efec; --ink2: #c3c2b7; --muted: #8f8e88;
    --grid: #383835; --axis: #5a5955; --data: #f0efec; --model: #d95926; --ours: #3987e5;
    --star: #c98500; --chip: #2a2a28; --chip-on: #3987e5; --chip-on-ink: #0b0b0b; --wash: #1f2c3d;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --panel: #232321; --ink: #f0efec; --ink2: #c3c2b7; --muted: #8f8e88;
  --grid: #383835; --axis: #5a5955; --data: #f0efec; --model: #d95926; --ours: #3987e5;
  --star: #c98500; --chip: #2a2a28; --chip-on: #3987e5; --chip-on-ink: #0b0b0b; --wash: #1f2c3d;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--surface); color: var(--ink);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 18px 16px 40px; }
h1 { font-size: 22px; margin: 4px 0 4px; }
p.lede { color: var(--ink2); margin: 0 0 14px; max-width: 900px; }
.controls { display: flex; flex-wrap: wrap; gap: 14px 22px; align-items: flex-start; background: var(--panel);
  border-radius: 10px; padding: 12px 14px; margin-bottom: 14px; }
.ctl { display: flex; flex-direction: column; gap: 6px; }
.ctl > span { font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; }
.chips { display: flex; flex-wrap: wrap; gap: 6px; }
button.chip { border: 1px solid var(--axis); background: var(--chip); color: var(--ink); border-radius: 999px;
  padding: 4px 11px; font: inherit; font-size: 13px; cursor: pointer; }
button.chip[aria-pressed="true"] { background: var(--chip-on); border-color: var(--chip-on); color: var(--chip-on-ink); }
button.chip:focus-visible { outline: 2px solid var(--ours); outline-offset: 2px; }
input[type=range] { width: 170px; accent-color: var(--ours); }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) 300px; gap: 16px; }
@media (max-width: 900px) { .layout { grid-template-columns: 1fr; } }
.card { background: var(--panel); border-radius: 10px; padding: 12px 14px; }
.big { font-size: 30px; font-weight: 700; line-height: 1.1; }
.sub { color: var(--ink2); font-size: 13px; }
.kv { display: grid; grid-template-columns: auto auto; gap: 3px 10px; font-size: 13px; margin-top: 8px; }
.kv div:nth-child(odd) { color: var(--muted); }
svg text { fill: var(--ink2); font: 11px system-ui, sans-serif; }
.legend { display: flex; gap: 16px; flex-wrap: wrap; font-size: 12px; color: var(--ink2); margin: 6px 2px 0; }
.legend i { display: inline-block; width: 22px; height: 3px; vertical-align: middle; margin-right: 6px; }
#tip { position: fixed; pointer-events: none; background: var(--ink); color: var(--surface); font-size: 12px;
  padding: 6px 8px; border-radius: 6px; opacity: 0; transition: opacity .1s; max-width: 260px; z-index: 5; }
.note { color: var(--muted); font-size: 12px; margin-top: 10px; }
.theme { position: absolute; top: 14px; right: 16px; }
details { margin-top: 12px; color: var(--ink2); font-size: 13px; }
table { border-collapse: collapse; font-size: 12px; margin-top: 6px; }
td, th { padding: 2px 8px; text-align: right; border-bottom: 1px solid var(--grid); }
</style>
</head>
<body>
<main>
<div style="position:relative">
<h1>Quasar drift explorer</h1>
<button class="chip theme" id="themeBtn" aria-label="Toggle dark mode">◐ theme</button>
</div>
<p class="lede">The Sun falls around the Milky Way, and its changing velocity tilts the light of every distant quasar a
few millionths of an arcsecond per year towards the Galactic Centre. Each Gaia quasar's error is ~90 times the signal,
so the pattern only appears when hundreds of thousands are averaged. Choose which quasars to use and the page refits
the pattern from Gaia DR3's proper motions.</p>

<div class="controls" role="group" aria-label="Fit controls">
  <div class="ctl"><span>Split quasars by</span><div class="chips" id="axisChips"></div></div>
  <div class="ctl"><span>Use these bins</span><div class="chips" id="binChips"></div></div>
  <div class="ctl"><span>Fraction of quasars</span><div class="chips" id="fracChips"></div></div>
  <div class="ctl"><span>VSH degrees fitted (l<sub>max</sub>)</span><div class="chips" id="lmaxChips"></div></div>
  <div class="ctl"><span>Arrow averaging radius: <b id="radVal"></b></span>
    <input type="range" id="radius" min="15" max="45" step="5" value="30" aria-label="Averaging radius in degrees"></div>
  <div class="ctl"><span>Show</span><div class="chips" id="showChips"></div></div>
  <div class="ctl"><span>Arrows and profile</span><div class="chips" id="cleanChips"></div></div>
</div>

<div class="layout">
  <div>
    <div class="card" style="padding:6px 6px 8px">
      <svg id="sky" viewBox="-300 -152 600 304" width="100%" role="img" aria-label="All-sky map of averaged quasar proper motions in galactic coordinates"></svg>
      <div class="legend">
        <span><i style="background:var(--data)"></i>averaged quasar proper motions (grey: fewer than half the usual quasars)</span>
        <span><i style="background:var(--model);height:5px"></i>fitted acceleration pattern (same averaging)</span>
        <span>★ Galactic Centre &nbsp; ⊕ fitted direction (1σ ellipse)</span>
        <span>arrow for 5 µas/yr: <svg width="44" height="10" style="vertical-align:middle"><line x1="2" y1="5" x2="38" y2="5" stroke="var(--data)" stroke-width="1.6" marker-end="url(#headInk)"/></svg></span>
      </div>
    </div>
    <div class="card" style="margin-top:14px">
      <svg id="prof" viewBox="0 0 640 230" width="100%" role="img" aria-label="Drift towards the Galactic Centre versus angular distance"></svg>
    </div>
  </div>
  <div>
    <div class="card" aria-live="polite">
      <div class="sub">Fitted drift (acceleration ÷ c)</div>
      <div class="big" id="amp">–</div>
      <div class="sub" id="ampSI"></div>
      <div class="kv" id="kv"></div>
    </div>
    <div class="card" style="margin-top:14px">
      <div class="sub">For comparison</div>
      <div class="kv" id="ref"></div>
    </div>
    <p class="note" id="checkNote"></p>
  </div>
</div>
<details><summary>What the page computes</summary>
<p>Gaia DR3's 1,215,942 five-parameter Gaia-CRF3 quasars, minus the 3,783 outliers of the full fit, are summed into
3,072 sky pixels (HEALPix nside 16) with their full 2×2 error weights, separately for each magnitude bin, colour bin and
three random tranches (1%, 9%, 90%). The page solves the weighted least-squares fit of vector spherical harmonics up to
the chosen degree (glide + rotation at l = 1, plus large-scale systematics at l = 2–3), exactly as the main analysis
but with each quasar placed at its pixel's centre. Errors are the fit's errors times the unit weight error of the full
fit. Arrows average all quasars within the chosen radius of each point of a coarser grid, so neighbouring arrows share
quasars. By default the arrows and the profile are drawn after subtracting the fit's own non-glide terms (the frame's
spin and, for l<sub>max</sub> ≥ 2, the larger-scale patterns of Gaia's systematic errors), which are bigger than the
acceleration itself; choose "raw" to see them. The profile below the map is the mean drift along (and across) the
great circle through the Galactic Centre in 10° rings.</p></details>
</main>
<div id="tip" role="tooltip"></div>
<script>
const P = __PAYLOAD__;
function f32(b64) { const s = atob(b64), u = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i); return new Float32Array(u.buffer); }
const BE = f32(P.basis_e), BN = f32(P.basis_n), SUMS = f32(P.sums);
const NT = P.n_terms, NPIX = P.pix_l.length, NSETS = P.n_sets, NTR = P.n_tranches;
const D2R = Math.PI / 180;
const state = { axis: "mag", bins: new Set([0, 1, 2, 3, 4]), frac: 2, lmax: 3, radius: 30, data: true, model: true, clean: true };

// ---------- geometry ----------
function wrap180(l) { return ((l + 180) % 360 + 360) % 360 - 180; }
function proj(l, b) { const lam = -wrap180(l) * D2R, ph = b * D2R, den = Math.sqrt(1 + Math.cos(ph) * Math.cos(lam / 2));
  return [100 * 2 * Math.SQRT2 * Math.cos(ph) * Math.sin(lam / 2) / den, -100 * Math.SQRT2 * Math.sin(ph) / den]; }
function unit(l, b) { l *= D2R; b *= D2R; return [Math.cos(b) * Math.cos(l), Math.cos(b) * Math.sin(l), Math.sin(b)]; }
function frames(l, b) { l *= D2R; b *= D2R; return { e: [-Math.sin(l), Math.cos(l), 0], n: [-Math.sin(b) * Math.cos(l), -Math.sin(b) * Math.sin(l), Math.cos(b)] }; }
function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function lonlat(v) { const r = Math.hypot(v[0], v[1], v[2]); return [((Math.atan2(v[1], v[0]) / D2R) % 360 + 360) % 360, Math.asin(v[2] / r) / D2R]; }
const PIXU = P.pix_l.map((l, i) => unit(l, P.pix_b[i])), PIXF = P.pix_l.map((l, i) => frames(l, P.pix_b[i]));

// ---------- data selection ----------
function selectedSums() {
  const out = new Float64Array(NPIX * 6);
  const ax = P.axes[state.axis], nb = ax.labels.length;
  for (let bin = 0; bin < nb; bin++) {
    if (!state.bins.has(bin)) continue;
    for (let t = 0; t <= state.frac; t++) {
      const set = ax.first + bin * NTR + t, off = set * NPIX * 6;
      for (let i = 0; i < NPIX * 6; i++) out[i] += SUMS[off + i];
    }
  }
  return out;
}

// ---------- linear algebra ----------
function cholSolve(A, b, n) { // returns {x, inv}
  const L = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = A[i * n + j]; for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
    if (i === j) { if (s <= 0) return null; L[i * n + i] = Math.sqrt(s); } else L[i * n + j] = s / L[j * n + j];
  }
  function solve(v) { const y = new Float64Array(n), x = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = v[i]; for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k]; y[i] = s / L[i * n + i]; }
    for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k]; x[i] = s / L[i * n + i]; }
    return x; }
  const x = solve(b), inv = new Float64Array(n * n);
  for (let j = 0; j < n; j++) { const e = new Float64Array(n); e[j] = 1; const c = solve(e); for (let i = 0; i < n; i++) inv[i * n + j] = c[i]; }
  return { x, inv };
}

function fitVSH(S) {
  const p = P.n_params[String(state.lmax)];
  const N = new Float64Array(p * p), r = new Float64Array(p);
  let nq = 0;
  for (let k = 0; k < NPIX; k++) {
    const o = k * 6, wee = S[o], wen = S[o + 1], wnn = S[o + 2], be = S[o + 3], bn = S[o + 4];
    nq += S[o + 5]; if (!(wee > 0)) continue;
    const ro = k * NT;
    for (let i = 0; i < p; i++) {
      const ei = BE[ro + i], ni = BN[ro + i];
      r[i] += ei * be + ni * bn;
      const ai = ei * wee + ni * wen, ci = ei * wen + ni * wnn;
      for (let j = 0; j <= i; j++) N[i * p + j] += ai * BE[ro + j] + ci * BN[ro + j];
    }
  }
  for (let i = 0; i < p; i++) for (let j = 0; j < i; j++) N[j * p + i] = N[i * p + j];
  const sol = cholSolve(N, r, p);
  if (!sol) return null;
  const u2 = P.uwe * P.uwe;
  function vec(M) { const v = [0, 0, 0], C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let a = 0; a < 3; a++) { for (let i = 0; i < p; i++) v[a] += M[a][i] * sol.x[i];
      for (let b = 0; b < 3; b++) { let s = 0; for (let i = 0; i < p; i++) for (let j = 0; j < p; j++) s += M[a][i] * sol.inv[i * p + j] * M[b][j]; C[a][b] = s * u2; } }
    return { v, C }; }
  return { glide: vec(P.m_glide), rot: vec(P.m_rot), nq, x: sol.x, p };
}

// Subtract the fit's own non-glide terms (frame spin, degree 2-3 patterns) from
// every pixel's weighted proper-motion sum: what is left is glide + noise + the
// smaller-scale patterns. Used for the pictures only, as in the README's keeper.
function cleanSums(S, res) {
  const out = Float64Array.from(S), keep = new Set(P.glide_idx);
  for (let k = 0; k < NPIX; k++) {
    const o = k * 6; if (!(S[o] > 0)) continue;
    let me = 0, mn = 0;
    for (let i = 0; i < res.p; i++) { if (keep.has(i)) continue; me += BE[k * NT + i] * res.x[i]; mn += BN[k * NT + i] * res.x[i]; }
    out[o + 3] -= S[o] * me + S[o + 1] * mn;
    out[o + 4] -= S[o + 1] * me + S[o + 2] * mn;
  }
  return out;
}

// ---------- averaging on the sky ----------
function capAverages(S, g) {
  const out = [], cr = Math.cos(state.radius * D2R);
  for (let c = 0; c < P.arrow_l.length; c++) {
    const uc = unit(P.arrow_l[c], P.arrow_b[c]), fc = frames(P.arrow_l[c], P.arrow_b[c]);
    let W = [0, 0, 0], b = [0, 0], bm = [0, 0], n = 0;
    for (let k = 0; k < NPIX; k++) {
      if (dot(uc, PIXU[k]) < cr) continue;
      const o = k * 6; if (!(S[o] > 0)) continue;
      const f = PIXF[k];
      // pixel frame vectors expressed in the centre frame
      const pe = [dot(f.e, fc.e), dot(f.e, fc.n)], pn = [dot(f.n, fc.e), dot(f.n, fc.n)];
      const wee = S[o], wen = S[o + 1], wnn = S[o + 2];
      // W_c += M^T W_p M with M = [pe; pn] (rows = pixel components)
      const m00 = pe[0], m01 = pe[1], m10 = pn[0], m11 = pn[1];
      const a00 = wee * m00 + wen * m10, a01 = wee * m01 + wen * m11, a10 = wen * m00 + wnn * m10, a11 = wen * m01 + wnn * m11;
      W[0] += m00 * a00 + m10 * a10; W[1] += m00 * a01 + m10 * a11; W[2] += m01 * a01 + m11 * a11;
      b[0] += m00 * S[o + 3] + m10 * S[o + 4]; b[1] += m01 * S[o + 3] + m11 * S[o + 4];
      if (g) { // glide model at the pixel centre, in the pixel frame, weighted
        const ge = dot(g, f.e), gn = dot(g, f.n), we = wee * ge + wen * gn, wn = wen * ge + wnn * gn;
        bm[0] += m00 * we + m10 * wn; bm[1] += m01 * we + m11 * wn; }
      n += S[o + 5];
    }
    const det = W[0] * W[2] - W[1] * W[1];
    if (!(det > 0) || n < 50) { out.push(null); continue; }
    const i00 = W[2] / det, i01 = -W[1] / det, i11 = W[0] / det;
    out.push({ l: P.arrow_l[c], b: P.arrow_b[c], n,
      e: i00 * b[0] + i01 * b[1], nn: i01 * b[0] + i11 * b[1], se: Math.sqrt(i00), sn: Math.sqrt(i11),
      me: i00 * bm[0] + i01 * bm[1], mn: i01 * bm[0] + i11 * bm[1] });
  }
  return out;
}

function profile(S, rot) {
  const gc = unit(P.sgra_lb[0], P.sgra_lb[1]), rings = [];
  for (let r = 0; r < 18; r++) rings.push({ W: [0, 0, 0], b: [0, 0], n: 0 });
  for (let k = 0; k < NPIX; k++) {
    const o = k * 6; if (!(S[o] > 0)) continue;
    const u = PIXU[k], c = Math.max(-1, Math.min(1, dot(u, gc))), th = Math.acos(c) / D2R;
    const r = Math.min(17, Math.floor(th / 10));
    let t = [gc[0] - c * u[0], gc[1] - c * u[1], gc[2] - c * u[2]]; const tn = Math.hypot(...t); if (tn < 1e-9) continue;
    t = t.map(x => x / tn);
    const f = PIXF[k], te = dot(t, f.e), tnn = dot(t, f.n), xe = -tnn, xn = te; // across = u x t
    const wee = S[o], wen = S[o + 1], wnn = S[o + 2];
    let be = S[o + 3], bn = S[o + 4];
    if (rot) { const w = rot; const mu = [w[1] * u[2] - w[2] * u[1], w[2] * u[0] - w[0] * u[2], w[0] * u[1] - w[1] * u[0]];
      const me = dot(mu, f.e), mn = dot(mu, f.n); be -= wee * me + wen * mn; bn -= wen * me + wnn * mn; }
    const R = rings[r];
    const a00 = wee * te + wen * tnn, a01 = wee * xe + wen * xn, a10 = wen * te + wnn * tnn, a11 = wen * xe + wnn * xn;
    R.W[0] += te * a00 + tnn * a10; R.W[1] += te * a01 + tnn * a11; R.W[2] += xe * a01 + xn * a11;
    R.b[0] += te * be + tnn * bn; R.b[1] += xe * be + xn * bn; R.n += S[o + 5];
  }
  return rings.map((R, i) => { const det = R.W[0] * R.W[2] - R.W[1] * R.W[1]; if (!(det > 0) || R.n < 500) return null;
    const i00 = R.W[2] / det, i01 = -R.W[1] / det, i11 = R.W[0] / det;
    return { th: i * 10 + 5, along: i00 * R.b[0] + i01 * R.b[1], ae: Math.sqrt(i00), across: i01 * R.b[0] + i11 * R.b[1], ce: Math.sqrt(i11), n: R.n }; });
}

// ---------- drawing ----------
const NS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
const tip = document.getElementById("tip");
function showTip(ev, html) { tip.innerHTML = html; tip.style.opacity = 1; tip.style.left = (ev.clientX + 14) + "px"; tip.style.top = (ev.clientY + 12) + "px"; }
function hideTip() { tip.style.opacity = 0; }
const ARROW_SCALE = 6.0; // svg units per uas/yr

function drawSky(avg, fitres) {
  const svg = document.getElementById("sky"); svg.innerHTML = "";
  const defs = el("defs", {}, svg);
  for (const [id, col] of [["headInk", "var(--data)"], ["headModel", "var(--model)"], ["headMuted", "var(--muted)"]]) {
    const m = el("marker", { id, viewBox: "0 0 10 10", refX: "8", refY: "5", markerWidth: "5", markerHeight: "5", orient: "auto-start-reverse" }, defs);
    el("path", { d: "M0,0 L10,5 L0,10 z", fill: col }, m); }
  const g0 = el("g", {}, svg);
  // graticule
  for (let l = -150; l <= 180; l += 30) { let d = ""; for (let b = -90; b <= 90; b += 3) { const [x, y] = proj(l === 180 ? 179.999 : l, b); d += (b === -90 ? "M" : "L") + x.toFixed(1) + "," + y.toFixed(1); } el("path", { d, fill: "none", stroke: "var(--grid)", "stroke-width": 0.6 }, g0); }
  for (const b of [-60, -30, 0, 30, 60]) { let d = ""; for (let l = -179.9; l <= 179.9; l += 3) { const [x, y] = proj(l, b); d += (l === -179.9 ? "M" : "L") + x.toFixed(1) + "," + y.toFixed(1); } el("path", { d, fill: "none", stroke: "var(--grid)", "stroke-width": 0.6 }, g0); }
  el("ellipse", { cx: 0, cy: 0, rx: 200 * Math.SQRT2, ry: 100 * Math.SQRT2, fill: "none", stroke: "var(--axis)" }, g0);
  for (const [l, t] of [[0, "l = 0°"], [90, "90°"], [270, "270°"]]) { const [x, y] = proj(l, 0); el("text", { x, y: y + 13, "text-anchor": "middle" }, g0).textContent = t; }
  // arrows
  function arrow(a, e, n, col, head, w, cls) {
    const amp = Math.hypot(e, n); if (!(amp > 0)) return null;
    const u = unit(a.l, a.b), f = frames(a.l, a.b), eps = 0.3 * D2R;
    const d = [(e * f.e[0] + n * f.n[0]) / amp, (e * f.e[1] + n * f.n[1]) / amp, (e * f.e[2] + n * f.n[2]) / amp];
    const u2 = [Math.cos(eps) * u[0] + Math.sin(eps) * d[0], Math.cos(eps) * u[1] + Math.sin(eps) * d[1], Math.cos(eps) * u[2] + Math.sin(eps) * d[2]];
    const [x0, y0] = proj(a.l, a.b), [lb1, bb1] = lonlat(u2), [x1, y1] = proj(lb1, bb1);
    let dx = x1 - x0, dy = y1 - y0; const nn = Math.hypot(dx, dy); if (!(nn > 0) || nn > 20) return null;
    dx /= nn; dy /= nn;
    return el("line", { x1: x0.toFixed(2), y1: y0.toFixed(2), x2: (x0 + dx * amp * ARROW_SCALE).toFixed(2), y2: (y0 + dy * amp * ARROW_SCALE).toFixed(2),
      stroke: col, "stroke-width": w, "marker-end": `url(#${head})`, class: cls }, svg);
  }
  for (const a of avg) { if (!a) continue;
    if (state.model) arrow(a, a.me, a.mn, "var(--model)", "headModel", 2.6, "m"); }
  const ns = avg.filter(a => a).map(a => a.n).sort((p, q) => p - q), nmed = ns[Math.floor(ns.length / 2)] || 0;
  for (const a of avg) { if (!a) continue;
    const sparse = a.n < 0.5 * nmed;
    if (state.data) { const ln = arrow(a, a.e, a.nn, sparse ? "var(--muted)" : "var(--data)", sparse ? "headMuted" : "headInk", sparse ? 1.0 : 1.3, "d");
      const [x0, y0] = proj(a.l, a.b);
      const hit = el("circle", { cx: x0, cy: y0, r: 9, fill: "transparent" }, svg);
      hit.addEventListener("mousemove", ev => showTip(ev, `(l, b) = (${wrap180(a.l).toFixed(0)}°, ${a.b.toFixed(0)}°)<br>${Math.round(a.n).toLocaleString()} quasars<br>east ${a.e.toFixed(1)} ± ${a.se.toFixed(1)}, north ${a.nn.toFixed(1)} ± ${a.sn.toFixed(1)} µas/yr<br>model: east ${a.me.toFixed(1)}, north ${a.mn.toFixed(1)}`));
      hit.addEventListener("mouseleave", hideTip); } }
  // Galactic Centre and fitted direction
  const [gx, gy] = proj(P.sgra_lb[0], P.sgra_lb[1]);
  el("text", { x: gx, y: gy + 5, "text-anchor": "middle", style: "font-size:18px;fill:var(--star);stroke:var(--ink);stroke-width:0.6" }, svg).textContent = "★";
  if (fitres) {
    const v = fitres.glide.v, [l, b] = lonlat(v);
    // 1-sigma ellipse via sampling the covariance in the tangent plane
    const C = fitres.glide.C, r = Math.hypot(...v), f = frames(l, b), u = unit(l, b);
    const ce = dot(f.e, mul(C, f.e)) / (r * r), cn = dot(f.n, mul(C, f.n)) / (r * r), cen = dot(f.e, mul(C, f.n)) / (r * r);
    const tr = ce + cn, dt = ce * cn - cen * cen, l1 = tr / 2 + Math.sqrt(Math.max(0, tr * tr / 4 - dt)), l2 = tr / 2 - Math.sqrt(Math.max(0, tr * tr / 4 - dt));
    const ang = 0.5 * Math.atan2(2 * cen, ce - cn);
    let d = "";
    for (let i = 0; i <= 72; i++) { const t = i / 72 * 2 * Math.PI, a = Math.sqrt(l1) * Math.cos(t), bb = Math.sqrt(Math.max(l2, 0)) * Math.sin(t);
      const de = a * Math.cos(ang) - bb * Math.sin(ang), dn = a * Math.sin(ang) + bb * Math.cos(ang);
      const w = [u[0] + de * f.e[0] + dn * f.n[0], u[1] + de * f.e[1] + dn * f.n[1], u[2] + de * f.e[2] + dn * f.n[2]];
      const [ll, bb2] = lonlat(w), [x, y] = proj(ll, bb2); d += (i ? "L" : "M") + x.toFixed(1) + "," + y.toFixed(1); }
    el("path", { d, fill: "var(--ours)", "fill-opacity": 0.18, stroke: "var(--ours)", "stroke-width": 1.2 }, svg);
    const [fx, fy] = proj(l, b);
    el("circle", { cx: fx, cy: fy, r: 3.5, fill: "var(--ours)", stroke: "var(--surface)", "stroke-width": 1 }, svg);
  }
}
function mul(C, v) { return [C[0][0] * v[0] + C[0][1] * v[1] + C[0][2] * v[2], C[1][0] * v[0] + C[1][1] * v[1] + C[1][2] * v[2], C[2][0] * v[0] + C[2][1] * v[1] + C[2][2] * v[2]]; }

function drawProfile(rings, gpar) {
  const svg = document.getElementById("prof"); svg.innerHTML = "";
  const X0 = 52, X1 = 620, Y0 = 20, Y1 = 190, ymin = -10, ymax = 14;
  const sx = th => X0 + (X1 - X0) * th / 180, sy = v => Y1 - (Y1 - Y0) * (v - ymin) / (ymax - ymin);
  for (let v = -10; v <= 14; v += 5) { el("line", { x1: X0, x2: X1, y1: sy(v), y2: sy(v), stroke: v === 0 ? "var(--axis)" : "var(--grid)" }, svg);
    el("text", { x: X0 - 6, y: sy(v) + 4, "text-anchor": "end" }, svg).textContent = v; }
  for (let th = 0; th <= 180; th += 30) el("text", { x: sx(th), y: Y1 + 16, "text-anchor": "middle" }, svg).textContent = th + "°";
  el("text", { x: (X0 + X1) / 2, y: Y1 + 34, "text-anchor": "middle" }, svg).textContent = "angular distance from the Galactic Centre";
  el("text", { x: X0, y: 12 }, svg).textContent = "drift towards the Galactic Centre (µas/yr): ● along, ○ sideways";
  let d = ""; for (let th = 0; th <= 180; th += 2) d += (th ? "L" : "M") + sx(th).toFixed(1) + "," + sy(gpar * Math.sin(th * D2R)).toFixed(1);
  el("path", { d, fill: "none", stroke: "var(--model)", "stroke-width": 2.4 }, svg);
  for (const r of rings) { if (!r) continue;
    const ya = Math.max(ymin, Math.min(ymax, r.along)), yc = Math.max(ymin, Math.min(ymax, r.across));
    el("line", { x1: sx(r.th) - 3, x2: sx(r.th) - 3, y1: sy(Math.min(ymax, r.along + r.ae)), y2: sy(Math.max(ymin, r.along - r.ae)), stroke: "var(--ours)", "stroke-width": 1.4 }, svg);
    const c = el("circle", { cx: sx(r.th) - 3, cy: sy(ya), r: 4, fill: "var(--ours)" }, svg);
    el("line", { x1: sx(r.th) + 3, x2: sx(r.th) + 3, y1: sy(Math.min(ymax, r.across + r.ce)), y2: sy(Math.max(ymin, r.across - r.ce)), stroke: "var(--muted)", "stroke-width": 1 }, svg);
    el("circle", { cx: sx(r.th) + 3, cy: sy(yc), r: 3, fill: "var(--surface)", stroke: "var(--muted)" }, svg);
    const hit = el("rect", { x: sx(r.th) - 14, y: Y0, width: 28, height: Y1 - Y0, fill: "transparent" }, svg);
    hit.addEventListener("mousemove", ev => showTip(ev, `${r.th - 5}–${r.th + 5}° from the centre<br>${Math.round(r.n).toLocaleString()} quasars<br>along ${r.along.toFixed(2)} ± ${r.ae.toFixed(2)}<br>sideways ${r.across.toFixed(2)} ± ${r.ce.toFixed(2)} µas/yr`));
    hit.addEventListener("mouseleave", hideTip); }
}

function fmt(x, d = 2) { return (x >= 0 ? "" : "−") + Math.abs(x).toFixed(d); }
function update() {
  document.getElementById("radVal").textContent = state.radius + "°";
  const S = selectedSums();
  const res = fitVSH(S);
  if (!res) { document.getElementById("amp").textContent = "–"; return; }
  const g = res.glide.v, C = res.glide.C, r = Math.hypot(...g);
  const sr = Math.sqrt(dot(g, mul(C, g))) / r;
  const [l, b] = lonlat(g);
  const gc = unit(P.sgra_lb[0], P.sgra_lb[1]), ang = Math.acos(Math.max(-1, Math.min(1, dot(g, gc) / r))) / D2R;
  document.getElementById("amp").textContent = `${r.toFixed(2)} ± ${sr.toFixed(2)} µas/yr`;
  document.getElementById("ampSI").textContent = `= (${(r * P.uas_to_m_s2 * 1e10).toFixed(2)} ± ${(sr * P.uas_to_m_s2 * 1e10).toFixed(2)}) × 10⁻¹⁰ m/s²`;
  const gx = g[0], sgx = Math.sqrt(C[0][0]);
  const vc = Math.sqrt(Math.max(gx, 0) / P.kms2_kpc_to_uas * P.r0_h1);
  const vce = vc * 0.5 * Math.hypot(sgx / Math.max(gx, 1e-9), P.r0_h1_err / P.r0_h1);
  const w = res.rot.v, sw = [0, 1, 2].map(i => Math.sqrt(res.rot.C[i][i]));
  document.getElementById("kv").innerHTML =
    `<div>quasars used</div><div>${Math.round(res.nq).toLocaleString()}</div>` +
    `<div>towards the centre g<sub>X</sub></div><div>${fmt(gx)} ± ${sgx.toFixed(2)}</div>` +
    `<div>along rotation g<sub>Y</sub></div><div>${fmt(g[1])} ± ${Math.sqrt(C[1][1]).toFixed(2)}</div>` +
    `<div>north g<sub>Z</sub></div><div>${fmt(g[2])} ± ${Math.sqrt(C[2][2]).toFixed(2)}</div>` +
    `<div>direction (l, b)</div><div>(${fmt(wrap180(l), 1)}°, ${fmt(b, 1)}°)</div>` +
    `<div>angle from Sgr A*</div><div>${ang.toFixed(1)}°</div>` +
    `<div>significance of g<sub>X</sub></div><div>${(gx / sgx).toFixed(1)} σ</div>` +
    `<div>circular speed √(g<sub>X</sub>c·R₀)</div><div>${vc.toFixed(0)} ± ${vce.toFixed(0)} km/s</div>` +
    `<div>frame rotation (gal.)</div><div>${w.map(x => fmt(x, 1)).join(", ")} µas/yr</div>`;
  const Sp = state.clean ? cleanSums(S, res) : S;
  const avg = capAverages(Sp, g);
  drawSky(avg, res);
  drawProfile(profile(Sp, state.clean ? null : res.rot.v), dot(g, gc));
}

function chips(id, items, isOn, onClick) {
  const box = document.getElementById(id); box.innerHTML = "";
  items.forEach((label, i) => { const b = document.createElement("button"); b.className = "chip"; b.textContent = label;
    b.setAttribute("aria-pressed", isOn(i) ? "true" : "false"); b.onclick = () => { onClick(i); renderControls(); update(); }; box.appendChild(b); });
}
function renderControls() {
  chips("axisChips", ["magnitude", "colour"], i => (i === 0) === (state.axis === "mag"), i => { state.axis = i === 0 ? "mag" : "colour";
    state.bins = new Set(P.axes[state.axis].labels.map((_, k) => k)); });
  chips("binChips", P.axes[state.axis].labels, i => state.bins.has(i), i => { if (state.bins.has(i)) { if (state.bins.size > 1) state.bins.delete(i); } else state.bins.add(i); });
  chips("fracChips", ["1%", "10%", "100%"], i => state.frac === i, i => { state.frac = i; });
  chips("lmaxChips", ["1", "2", "3"], i => state.lmax === i + 1, i => { state.lmax = i + 1; });
  chips("showChips", ["data arrows", "model arrows"], i => i === 0 ? state.data : state.model, i => { if (i === 0) state.data = !state.data; else state.model = !state.model; });
  chips("cleanChips", ["Gaia's patterns removed", "raw"], i => (i === 0) === state.clean, i => { state.clean = i === 0; });
}
document.getElementById("radius").addEventListener("input", e => { state.radius = +e.target.value; update(); });
document.getElementById("themeBtn").onclick = () => { const r = document.documentElement, dark = r.dataset.theme ? r.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches; r.dataset.theme = dark ? "light" : "dark"; };
const R = P.reference;
document.getElementById("ref").innerHTML =
  `<div>this project, all quasars</div><div>${R.ours.amp.toFixed(2)} ± ${R.ours.amp_err.toFixed(2)}</div>` +
  `<div>Klioner et al. 2021</div><div>${R.klioner.amp.toFixed(2)} ± ${R.klioner.amp_err.toFixed(2)}</div>` +
  `<div>predicted, H1's R₀ + Sgr A*</div><div>${R.pred_h1.g.toFixed(2)} ± ${R.pred_h1.sd.toFixed(2)}</div>`;
document.getElementById("checkNote").textContent = `Self-check: with every quasar, this page's binned fit gives a glide within ${P.check.max_abs_diff_uas_yr.toFixed(3)} µas/yr of the exact fit (l_max = 3). Bootstrap errors of the full fit are ~5% larger than the fit errors shown here.`;
renderControls(); update();
</script>
</body>
</html>
"""
