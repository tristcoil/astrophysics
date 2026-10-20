"""The Edenhofer et al. (2024) map against a second one: Vergely, Lallement & Cox (2022).

Vergely et al. invert extinctions of 35 million stars (Gaia EDR3 + 2MASS,
inter-calibrated with spectroscopic ones) on a Cartesian grid with a 10 pc
correlation length (5 pc sampling): a different method, a different stellar
catalogue and a coarser resolution. Their density is A0 at 550 nm per pc; we
divide by 2.658 (ZGR23's curve at 552 nm) to get the map's ZGR23 units.

1. ``offset``: the shift between the two maps' mid-planes that best aligns
   them (tests our reading of the cube's header);
2. ``projections``: face-on views of both maps, the same box and stretch;
3. ``sightlines``: extinction integrated from the Sun to 100-1200 pc along
   3,072 directions in both maps (cf. Edenhofer et al. Fig. 7);
4. ``elongation``: how far the density stays correlated along the line of
   sight versus across it, as a function of distance. Distance errors smear
   clouds radially ("fingers of god"); a ratio near 1 means no smearing.
"""

from __future__ import annotations

import json

import healpy as hp
import numpy as np
from scipy.ndimage import map_coordinates

import catalogs
import config as C
import grids
from deredden import Extinction


def vergely_zgr23() -> dict:
    v = catalogs.vergely()
    return {**v, "cube": v["cube"] / C.ZGR23_TO_A550}


def _sample(cube, axes, pts):
    """Trilinear sample of a cube with uniform axes at points (N, 3); zero outside."""
    idx = [(pts[:, i] - axes[i][0]) / (axes[i][1] - axes[i][0]) for i in range(3)]
    return map_coordinates(cube, idx, order=1, mode="constant", cval=0.0)


def offset() -> dict:
    """Cross-correlate the two maps' face-on projections (|z| <= 400 pc), 5 pc cells."""
    from scipy.signal import fftconvolve
    v = vergely_zgr23()
    g = grids.build_main()
    xs = g["vol5_xs"]
    ours = g["vol5"][:, :, np.abs(g["vol5_zs"]) <= 400].sum(axis=2) * 5.0
    sel = (v["x"] >= xs[0] - 1e-6) & (v["x"] <= xs[-1] + 1e-6)
    vx = v["x"][sel]
    theirs = v["cube"][np.ix_(sel, sel, np.abs(v["z"]) <= 400)].sum(axis=2) * 5.0
    # interpolate theirs onto our cell centres (their nodes sit at multiples of 5 pc, ours at 2.5 + 5k)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    th = map_coordinates(theirs, [(X - vx[0]) / 5.0, (Y - vx[0]) / 5.0], order=1)
    R = np.hypot(X, Y)
    w = (R > 100) & (R < 1150)
    a = np.where(w, np.log10(ours + 1e-3), 0.0)
    b = np.where(w, np.log10(th + 1e-3), 0.0)
    a -= a[w].mean() * w
    b -= b[w].mean() * w
    cc = fftconvolve(a, b[::-1, ::-1], mode="same")
    i, j = np.unravel_index(np.argmax(cc), cc.shape)
    c0 = np.array(cc.shape) // 2
    # sub-pixel peak by a parabola through the maximum
    def sub(f_m, f_0, f_p):
        den = f_m - 2 * f_0 + f_p
        return 0.0 if den == 0 else 0.5 * (f_m - f_p) / den
    dx = (i - c0[0] + sub(cc[i - 1, j], cc[i, j], cc[i + 1, j])) * 5.0
    dy = (j - c0[1] + sub(cc[i, j - 1], cc[i, j], cc[i, j + 1])) * 5.0
    r = np.corrcoef(a[w], b[w])[0, 1]
    return {"shift_x_pc": float(dx), "shift_y_pc": float(dy), "log_corr": float(r),
            "note": "positive shift = Vergely's structures sit at larger x / y than ours"}


def projections() -> dict:
    """Face-on A_V (|z| <= 400 pc) of both maps on 5 pc cells, saved for the figure."""
    v = vergely_zgr23()
    g = grids.build_main()
    xs = g["vol5_xs"]
    ours = g["vol5"][:, :, np.abs(g["vol5_zs"]) <= 400].sum(axis=2) * 5.0
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    sel = np.abs(v["z"]) <= 400
    theirs_full = v["cube"][:, :, sel].sum(axis=2) * 5.0
    th = map_coordinates(theirs_full, [(X - v["x"][0]) / 5.0, (Y - v["y"][0]) / 5.0], order=1)
    np.savez_compressed(C.CACHE_DIR / "compare_projections.npz", xs=xs, ours=ours.astype(np.float32),
                        theirs=th.astype(np.float32))
    inside = np.hypot(X, Y) < 1240
    return {"total_ratio_vergely_over_ours": float(th[inside].sum() / ours[inside].sum()),
            "pearson_log": float(np.corrcoef(np.log10(ours[inside] + 1e-3), np.log10(th[inside] + 1e-3))[0, 1])}


def sightlines(nside: int = 16, step: float = 1.0) -> dict:
    """Integrated extinction from the Sun in both maps, every 50 pc from 100 to 1200 pc."""
    v = vergely_zgr23()
    ext = Extinction()
    l, b = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    dists = np.arange(100.0, 1201.0, 50.0)
    ours = np.array([ext.zgr23(l, b, np.full_like(l, d)) for d in dists])          # (nd, npix)
    u = np.stack([np.cos(np.radians(b)) * np.cos(np.radians(l)), np.cos(np.radians(b)) * np.sin(np.radians(l)),
                  np.sin(np.radians(b))], axis=1)
    r = np.arange(0.5 * step, dists[-1], step)
    theirs = np.empty_like(ours)
    inside_box = np.empty_like(ours, dtype=bool)
    for j in range(len(l)):
        pts = u[j][None, :] * r[:, None]
        rho = _sample(v["cube"], (v["x"], v["y"], v["z"]), pts)
        cum = np.cumsum(rho) * step
        theirs[:, j] = np.interp(dists, r, cum)
        far = dists[:, None] * u[j][None, :]
        inside_box[:, j] = (np.abs(far[:, 2]) <= 400) & (np.abs(far[:, 0]) <= 1500) & (np.abs(far[:, 1]) <= 1500)
    av_o, av_t = C.ZGR23_TO_AV * ours, C.ZGR23_TO_AV * theirs
    ok = inside_box & (av_o > 0.04)
    ratio = av_t[ok] / av_o[ok]
    bins = [(0.04, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 30.0)]
    by_level = [{"av_lo": lo, "av_hi": hi, "median_ratio": float(np.median(ratio[(av_o[ok] >= lo) & (av_o[ok] < hi)])),
                 "n": int(((av_o[ok] >= lo) & (av_o[ok] < hi)).sum())} for lo, hi in bins]
    np.savez_compressed(C.CACHE_DIR / "compare_sightlines.npz", l=l, b=b, d=dists, ours=av_o, theirs=av_t,
                        inside=inside_box)
    return {"n_pairs": int(ok.sum()), "median_ratio_vergely_over_ours": float(np.median(ratio)),
            "ratio_p16_p84": [float(np.quantile(ratio, 0.16)), float(np.quantile(ratio, 0.84))],
            "by_level": by_level}


def _corr_length(lags, corr):
    """Lag at which the correlation falls to 0.5 (linear interpolation), NaN if never."""
    below = np.flatnonzero(corr < 0.5)
    if below.size == 0:
        return np.nan
    k = below[0]
    if k == 0:
        return lags[0]
    return float(np.interp(0.5, [corr[k], corr[k - 1]], [lags[k], lags[k - 1]]))


def _highpass(cube: np.ndarray, cell: float, sigma_pc: float = 30.0) -> np.ndarray:
    """ln(density) minus its Gaussian-smoothed self: removes the bubble/disk-scale trends."""
    from scipy.ndimage import gaussian_filter
    ln = np.log(np.maximum(cube, 1e-6)).astype(np.float32)
    return ln - gaussian_filter(ln, sigma_pc / cell, mode="nearest")


def elongation(seed: int = 3, n: int = 60000) -> dict:
    """Correlation length of the small-scale density field along vs across the line of sight.

    Both maps on 5 pc cells (Edenhofer's block-averaged 2.5 pc grid,
    Vergely's native sampling); large scales removed by subtracting a 30 pc
    Gaussian smoothing of ln(density). Points are drawn in shells around the
    Sun (|z| < 150 pc); each is compared with the value a lag away along the
    radial and along the horizontal tangential direction.
    """
    rng = np.random.default_rng(seed)
    g = grids.build_main()
    v = vergely_zgr23()
    fields = {
        "edenhofer": (_highpass(g["vol5"], 5.0), (g["vol5_xs"], g["vol5_xs"], g["vol5_zs"]), g["vol5"]),
        "vergely": (_highpass(v["cube"], 5.0), (v["x"], v["y"], v["z"]), v["cube"]),
    }
    lags = np.array([2.5, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50, 60, 80, 100])
    shells = [(100, 200), (200, 350), (350, 500), (500, 700), (700, 900), (900, 1100)]
    out = {name: [] for name in fields}
    for lo, hi in shells:
        r = rng.uniform(lo, hi, n)
        phi = rng.uniform(0, 2 * np.pi, n)
        z = rng.uniform(-150, 150, n)
        rxy = np.sqrt(np.maximum(r**2 - z**2, 1.0))
        p = np.c_[rxy * np.cos(phi), rxy * np.sin(phi), z]
        rhat = p / np.linalg.norm(p, axis=1)[:, None]
        that = np.c_[-np.sin(phi), np.cos(phi), np.zeros(n)]
        for name, (hp_field, axes, raw) in fields.items():
            dense = _sample(raw, axes, p) > 1e-4          # only where there is some dust
            base = _sample(hp_field, axes, p)
            corr = {"radial": [], "tangential": []}
            for lag in lags:
                for key, uvec in (("radial", rhat), ("tangential", that)):
                    a = _sample(hp_field, axes, p + lag * uvec)
                    corr[key].append(float(np.corrcoef(base[dense], a[dense])[0, 1]))
            lr = _corr_length(lags, np.array(corr["radial"]))
            lt = _corr_length(lags, np.array(corr["tangential"]))
            out[name].append({"r_lo": lo, "r_hi": hi, "n_points": int(dense.sum()), "corr_len_radial_pc": lr,
                              "corr_len_tangential_pc": lt,
                              "elongation": lr / lt if np.isfinite(lr) and np.isfinite(lt) else float("nan"),
                              "corr_radial": corr["radial"], "corr_tangential": corr["tangential"]})
    out["lags_pc"] = lags.tolist()
    return out


def run() -> dict:
    res = {"offset": offset(), "projections": projections(), "sightlines": sightlines(), "elongation": elongation()}
    (C.CACHE_DIR / "compare.json").write_text(json.dumps(res, indent=1, default=float))
    return res


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=float))
