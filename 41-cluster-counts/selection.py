"""The X-ray selection function: which clusters did eROSITA detect?

eRASS1's own selection function is a Gaussian-process classifier trained on
simulations of the survey (Clerc et al. 2024). The model is not public, but
Ghirardini et al. (2024, Fig. 4) plot nine slices of it, the probability of
detection P(I | CR, z, T) against the unabsorbed count rate CR at three
redshifts and three exposure times, as a *vector* figure. We extract those
curves exactly from the PDF's drawing commands and fit them with a model that
is easy to interpolate:

    P(I) = P_max * Phi( ln(CR * a(N_H) / a(N_H,ref) * T / S_50) / s )

Phi is the normal cumulative distribution. S_50 is the number of photons
(count rate x exposure) at which half of all clusters are detected; it grows
slowly with exposure, because a longer stare also collects more background
to beat. s sets how gradual the transition is, and P_max < 1 means that
even bright, nearby clusters are sometimes missed on short exposures. At each
of the three redshifts, ln S_50, s and ln(1 - P_max) are linear in
ln(T / 150 s). Between redshifts we interpolate linearly; beyond them, the
default extrapolates linearly from the nearest pair, and an alternative
holds the parameters flat (a robustness test in the README).

Galactic gas absorbs soft X-rays, so the same unabsorbed count rate gives
fewer photons where N_H is high. a(N_H), the fraction of the count rate that
survives absorption, is calibrated on the catalogue itself: it lists both the
absorbed and the unabsorbed count rate of every cluster.

A second published figure (Clerc et al. 2024) shows the same model at 200 s
and 2000 s; it is a raster image, which we digitise to check the
interpolation in redshift and the extrapolation to long exposures.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.optimize import least_squares
from scipy.special import ndtr

import config as C

CURVES_CACHE = C.CACHE_DIR / "fig4_selection_curves.json"
Z_NODES = (0.15, 0.30, 0.50)
T_NODES = (50.0, 150.0, 500.0)
T0 = 150.0

# Colours and dash patterns of Ghirardini et al. Fig. 4 (from its legend).
_FIG4_T = {"rgb(17.254639%, 62.744141%, 17.254639%)": 50.0,
           "rgb(100%, 49.803162%, 5.490112%)": 150.0,
           "rgb(12.156677%, 46.665955%, 70.587158%)": 500.0}
_FIG4_Z = {None: 0.30, "11.1 4.8": 0.50, "3 4.95": 0.15}
# Axis calibration from the tick marks (PDF points): x major ticks at
# 10^-2, 10^-1, 10^0, 10^1 every 226.074; y ticks P = 0 ... 1 every 113.718.
_X0, _DX = 210.203125, 226.074219
_Y0, _DY = 112.1875, 680.777344 - 112.1875


# ---------------------------------------------------------------------------
# Absorption
# ---------------------------------------------------------------------------

#: a(N_H) = exp(-k (N_H / 1e20)^p), fitted to CR500 / CR500_NH0 of 5,190
#: cosmology-sample clusters against their HI4PI N_H (catalog.absorption_fit
#: recomputes it; rms 0.026 in ln a, no trend with redshift).
ABSORPTION = {"k": 0.05592, "p": 0.7982}


def absorbed_fraction(nh, k=None, p=None):
    k = ABSORPTION["k"] if k is None else k
    p = ABSORPTION["p"] if p is None else p
    return np.exp(-k * (np.asarray(nh, float) / 1e20) ** p)


# ---------------------------------------------------------------------------
# Reading the published curves
# ---------------------------------------------------------------------------


def extract_fig4_curves() -> dict:
    """The nine P(I | CR) curves of Ghirardini et al. Fig. 4, from the PDF.

    Returns {(z, T): (log10 CR array, P array)}. Uses poppler's pdftocairo to
    turn the PDF into SVG paths; the result is cached as JSON.
    """
    if CURVES_CACHE.exists():
        raw = json.loads(CURVES_CACHE.read_text())
        return {tuple(map(float, k.split("_"))): (np.array(v[0]), np.array(v[1])) for k, v in raw.items()}
    if shutil.which("pdftocairo") is None:
        raise RuntimeError("pdftocairo (poppler-utils) is needed once to read the vector figure")
    with tempfile.TemporaryDirectory() as tmp:
        svg = f"{tmp}/fig4.svg"
        subprocess.run(["pdftocairo", "-svg", str(C.RAW_DIR / "ghirardini2024_fig4_selection.pdf"), svg], check=True)
        text = open(svg).read()
    curves = {}
    for path in re.findall(r"<path[^>]*>", text):
        colour = re.search(r'stroke="(rgb\([^)]*\))"', path)
        if not colour or colour.group(1) not in _FIG4_T:
            continue
        xy = np.array(re.findall(r"([-\d.]+) ([-\d.]+)", re.search(r' d="([^"]*)"', path).group(1)), float)
        if len(xy) < 10:  # a legend swatch
            continue
        dash = re.search(r'stroke-dasharray="([^"]*)"', path)
        key = (_FIG4_Z[dash.group(1) if dash else None], _FIG4_T[colour.group(1)])
        curves[key] = (-2 + (xy[:, 0] - _X0) / _DX, (xy[:, 1] - _Y0) / _DY)
    if len(curves) != 9:
        raise RuntimeError(f"expected 9 curves in Fig. 4, found {len(curves)}")
    CURVES_CACHE.write_text(json.dumps({f"{k[0]}_{k[1]}": [v[0].tolist(), v[1].tolist()] for k, v in curves.items()}))
    return curves


def digitize_clerc() -> dict:
    """CR at P = 0.25, 0.5, 0.75 for the four curves of Clerc et al.'s figure.

    Raster image (600x450). Axes from the major ticks: 98.25 px per decade
    (10^-3 at x = 96), P = (371 - y) / 302. Black: 200 s, blue: 2000 s;
    solid z = 0.2, dashed z = 0.5 (the dashed curve of each pair lies to the
    left, higher redshift being easier to detect at fixed count rate).
    Model conditions: N_H = 3e20 cm^-2 and the nominal background.
    """
    import matplotlib.image as mpimg

    rgb = (mpimg.imread(C.RAW_DIR / "clerc2024_selection_200s_2000s.png")[..., :3] * 255).round()
    black = rgb.max(axis=2) < 70
    blue = (rgb[..., 0] < 70) & (rgb[..., 1] < 70) & (rgb[..., 2] > 200)
    black[290:, 340:] = blue[290:, 340:] = False  # legend
    out = {}
    for colour, mask, T in (("black", black, 200.0), ("blue", blue, 2000.0)):
        for P in (0.25, 0.5, 0.75):
            y = int(round(371 - 302 * P))
            xs = []
            for yy in (y - 1, y, y + 1):
                cols = np.where(mask[yy, 77:539])[0] + 77
                if cols.size:
                    runs = np.split(cols, np.where(np.diff(cols) > 2)[0] + 1)
                    xs.append([r.mean() for r in runs])
            # the two curves at this height: leftmost = z 0.5, rightmost = z 0.2
            pts = sorted(x for row in xs for x in row)
            if len(pts) < 2:
                continue
            groups = np.split(np.array(pts), np.where(np.diff(pts) > 4)[0] + 1)
            centres = [g.mean() for g in groups]
            for z, x in ((0.5, min(centres)), (0.2, max(centres))):
                out.setdefault((z, T), {})[P] = 10 ** (-3 + (x - 96.0) / 98.25)
    return out


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------


def _node_values(params, z):
    """Piecewise-linear interpolation of node parameters (3 nodes x 6) to z."""
    zn = np.array(Z_NODES)
    z = np.atleast_1d(np.asarray(z, float))
    lo = np.where(z < zn[1], 0, 1)
    w = (z - zn[lo]) / (zn[lo + 1] - zn[lo])
    return params[lo] + w[:, None] * (params[lo + 1] - params[lo])


def _p_detect(q, g, lncr, z, T):
    """Detection probability for node parameters q (3x6) and curvature g."""
    p = _node_values(q, z)
    x = np.log(T / T0)
    ln_s50 = p[:, 0] + p[:, 1] * x + g * _hinge(T)
    s = np.maximum(p[:, 2] + p[:, 3] * x, 0.15)
    pmax = 1 - np.exp(np.minimum(p[:, 4] + p[:, 5] * x, 0.0))
    return pmax * ndtr((lncr + np.log(T) - ln_s50) / s)


def _hinge(T):
    """ln(T / 500 s) above 500 s, else 0: where Fig. 4 stops constraining."""
    return np.maximum(np.log(np.asarray(T, float) / T_NODES[-1]), 0.0)


def fit_nodes(curves: dict | None = None, use_clerc: bool = True) -> dict:
    """Fit the selection model to the published curves.

    Ghirardini et al.'s nine vector curves (weight 1/0.005 in P, their
    reading precision) fix six parameters at each redshift node. Clerc et
    al.'s raster curves (weight 1/0.02) add the long-exposure behaviour
    through one extra slope g beyond 500 s, where Fig. 4 ends (ln S_50 gains
    g ln(T / 500 s)), shared by all redshifts. Their N_H of 3e20 is converted
    to the 3.5e20 reference.
    """
    curves = extract_fig4_curves() if curves is None else curves
    lncr, P, T, Z = [], [], [], []
    for (z, t), (lx, p) in curves.items():
        lncr.append(lx * np.log(10)); P.append(p); T.append(np.full(p.size, t)); Z.append(np.full(p.size, z))
    lncr, P, T, Z = (np.concatenate(v) for v in (lncr, P, T, Z))
    if use_clerc:
        shift = np.log(absorbed_fraction(3e20) / absorbed_fraction(C.NH_REF))
        pts = [(z, t, pp, np.log(cr) + shift) for (z, t), d in digitize_clerc().items() for pp, cr in d.items()]
        cz, ct, cp, clncr = (np.array(v) for v in zip(*pts))
    lo = [0, -1, 0.2, -1, -15, -6]
    hi = [8, 1, 2, 1, -0.5, 2]
    start = [np.log(65.0), 0.12, 0.68, -0.05, np.log(0.01), -1.5]

    def resid(v):
        q, g = v[:18].reshape(3, 6), (v[18] if use_clerc else 0.0)
        r = (_p_detect(q, g, lncr, Z, T) - P) / 0.005
        if use_clerc:
            r = np.concatenate([r, (_p_detect(q, g, clncr, cz, ct) - cp) / 0.02])
        return r

    nvar = 19 if use_clerc else 18
    sol = least_squares(resid, (start * 3 + [0.0])[:nvar],
                        bounds=((lo * 3 + [-1.0])[:nvar], (hi * 3 + [1.0])[:nvar]))
    q = sol.x[:18].reshape(3, 6)
    g = float(sol.x[18]) if use_clerc else 0.0
    rms = {}
    for z in Z_NODES:
        sel = Z == z
        d = _p_detect(q, g, lncr[sel], Z[sel], T[sel]) - P[sel]
        rms[z] = (float(np.sqrt(np.mean(d**2))), float(np.max(np.abs(d))))
    return {"z": np.array(Z_NODES), "params": q, "g": g, "rms": rms, "use_clerc": use_clerc}


class SelectionModel:
    """P(detected | unabsorbed count rate, redshift, exposure, N_H)."""

    def __init__(self, nodes: dict | None = None, extrapolation: str = "linear", s50_scale: float = 1.0):
        self.nodes = fit_nodes() if nodes is None else nodes
        self.extrapolation = extrapolation
        self.s50_scale = s50_scale  # robustness tests: shift the threshold
        self.a_ref = float(absorbed_fraction(C.NH_REF))

    def params_at(self, z):
        """Model parameters (a..f) at redshift(s) z, shape (len(z), 6)."""
        z = np.atleast_1d(np.asarray(z, float))
        zn, pn = self.nodes["z"], self.nodes["params"]
        zz = np.clip(z, zn[0], zn[-1]) if self.extrapolation == "flat" else z
        out = np.empty((z.size, pn.shape[1]))
        for j in range(pn.shape[1]):
            # piecewise linear, extended linearly beyond the end nodes
            out[:, j] = np.where(zz < zn[1],
                                 pn[0, j] + (zz - zn[0]) * (pn[1, j] - pn[0, j]) / (zn[1] - zn[0]),
                                 pn[1, j] + (zz - zn[1]) * (pn[2, j] - pn[1, j]) / (zn[2] - zn[1]))
        return out

    def __call__(self, cr, z, t, nh=None):
        """Detection probability; arguments broadcast against each other."""
        cr, z, t = np.broadcast_arrays(np.asarray(cr, float), np.asarray(z, float), np.asarray(t, float))
        eff = cr if nh is None else cr * absorbed_fraction(np.broadcast_to(nh, cr.shape)) / self.a_ref
        p = self.params_at(z.ravel())
        x = np.log(t.ravel() / T0)
        ln_s50 = p[:, 0] + p[:, 1] * x + self.nodes.get("g", 0.0) * _hinge(t.ravel()) + np.log(self.s50_scale)
        s = np.maximum(p[:, 2] + p[:, 3] * x, 0.15)
        pmax = 1 - np.exp(np.minimum(p[:, 4] + p[:, 5] * x, 0.0))
        val = pmax * ndtr((np.log(eff.ravel()) + np.log(t.ravel()) - ln_s50) / s)
        return val.reshape(cr.shape)

    def cr50(self, z, t, nh=None):
        """Count rate at which the detection probability is half of P_max."""
        p = self.params_at(np.atleast_1d(z))
        x = np.log(np.asarray(t, float) / T0)
        ln_s50 = p[:, 0] + p[:, 1] * x + self.nodes.get("g", 0.0) * _hinge(t) + np.log(self.s50_scale)
        cr = np.exp(ln_s50) / np.asarray(t, float)
        if nh is not None:
            cr = cr * self.a_ref / absorbed_fraction(nh)
        return cr


def validate(model: SelectionModel | None = None) -> dict:
    """How well the model reproduces Fig. 4 and predicts Clerc et al.'s curves."""
    model = SelectionModel() if model is None else model
    curves = extract_fig4_curves()
    fig4 = {}
    for (z, T), (lx, P) in curves.items():
        pred = model(10**lx, z, T)
        fig4[f"z={z}, T={T:.0f}s"] = float(np.max(np.abs(pred - P)))
    clerc = {}
    for (z, T), pts in digitize_clerc().items():
        for P, cr in pts.items():
            # invert the model at this probability, at Clerc's N_H = 3e20
            grid = np.geomspace(1e-3, 50, 4000)
            pm = model(grid, z, T, nh=3e20)
            cr_model = float(np.interp(P, pm, grid)) if pm.max() > P else np.nan
            clerc[f"z={z}, T={T:.0f}s, P={P}"] = {"published": cr, "model": cr_model,
                                                 "ratio": cr_model / cr}
    return {"fig4_max_abs_dP": fig4, "clerc": clerc, "node_rms": model.nodes["rms"]}
