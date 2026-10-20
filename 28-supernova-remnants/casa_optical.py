"""Cas A's optical shell on Palomar plates (1954, 1991) and Pan-STARRS (~2012).

Cas A's optical knots are shocked clumps of oxygen- and sulphur-rich debris.
Each one lights up when the reverse shock hits it and fades within years to
decades, so the knots of 1954 are mostly not the knots of 2012. Tracking
small patches (as for the Crab) gave no coherent motions. What survives is
the shell as a whole: we find the single stretch about the centre (plus a
small shift) that best maps the earlier image onto the later one, using the
masked correlation of ``motions.py`` with stars excluded.

For a shell coasting since date T, a stretch s between dates t1 and t2 means
s = (t2 - T) / (t1 - T), so T = t1 - (t2 - t1) / (s - 1).

**Result: it does not work for Cas A on these images.** The correlation
peaks at s = 1.00-1.01 for every pair, red or near-infrared: what the plates
and Pan-STARRS have in common is structure that does not move (slow
circumstellar "flocculi" and foreground nebulosity), not the expanding shell,
whose knots are replaced rather than moved over 20-60 years. The same code on
the Crab finds s = 1.07, exactly what the filament vectors predict, so the
method itself is sound. Thorstensen et al. (2001) needed deep Palomar 5-m
plates and knot-by-knot identification.
"""

from __future__ import annotations

import json

import numpy as np
from scipy import ndimage, optimize

import astrometry as A
import config
import motions as M
import plates

#: The annulus holding the bright main shell (radius ~100-130") and the
#: fainter outer knots, about Thorstensen's centre (arcsec).
SHELL = (55.0, 175.0)
PAIRS = {
    "red 1954-1991": ("casa_poss1_red", "casa_poss2_red"),
    "red 1991-2012": ("casa_poss2_red", "casa_ps1_r"),
    "red 1954-2012": ("casa_poss1_red", "casa_ps1_r"),
    "near-IR 1991-2012": ("casa_poss2_ir", "casa_ps1_z"),
}
#: Thorstensen et al.'s outer-knot date: the stretch a coasting shell would show.
T_EXPECTED = 1671.3


def scan(p1: M.Prepared, p2: M.Prepared, grid: A.Grid, region: np.ndarray, centre=(0.0, 0.0),
         s_values=np.arange(0.98, 1.30, 0.005)) -> tuple[list, list]:
    """Correlation of image 2 with image 1 stretched by s about ``centre``
    (arcsec), for every s: the whole story in one curve."""
    out = []
    for s in s_values:
        w, ok = _warp(p1.image, p1.valid, grid, s, centre[0] * (1 - s), centre[1] * (1 - s))
        out.append(_ncc(w, p2.image, ok & p2.valid & region))
    return list(map(float, s_values)), out


def crab_control() -> dict:
    """The same global stretch on the Crab (1951 -> 2012), about its fitted
    centre: the method should find 1 + (t2 - t1) / (t1 - T) with T = 1106."""
    import crab as CR

    cr = CR.Crab()
    res = CR.analyse()
    a, b = CR.PAIRS[CR.MAIN]
    p1, p2 = cr.prepared(a), cr.prepared(b)
    s_vals, cc = scan(p1, p2, cr.grid, cr.region, centre=res["main"]["centre"],
                      s_values=np.arange(0.98, 1.15, 0.005))
    t1, t2 = cr.regs[a].epoch, cr.regs[b].epoch
    return {"s": s_vals, "ncc": cc, "s_best": s_vals[int(np.argmax(cc))],
            "s_expected": 1.0 + (t2 - t1) / (t1 - res["main"]["T"]), "t1": t1, "t2": t2}


def _warp(img: np.ndarray, valid: np.ndarray, grid: A.Grid, s: float, dxi: float, deta: float):
    """Image 1 as it would look stretched by s about the centre and shifted
    by (dxi, deta) arcsec: sample it at c + (p - c - d) / s."""
    n = grid.n
    c = (n - 1) / 2.0
    rr, cc = np.mgrid[:n, :n].astype(float)
    dc, dr = -dxi / grid.step, deta / grid.step  # east is -column
    src_r = c + (rr - c - dr) / s
    src_c = c + (cc - c - dc) / s
    out = ndimage.map_coordinates(img, [src_r, src_c], order=1, mode="constant", cval=0.0)
    ok = ndimage.map_coordinates(valid.astype(float), [src_r, src_c], order=1, mode="constant", cval=0.0) > 0.99
    return out, ok


def _ncc(a, b, m):
    a, b = a[m], b[m]
    a = a - a.mean()
    b = b - b.mean()
    return float(np.sum(a * b) / np.sqrt(np.sum(a * a) * np.sum(b * b)))


def fit_stretch(p1: M.Prepared, p2: M.Prepared, grid: A.Grid, region: np.ndarray) -> dict:
    """Best stretch and shift from image 1 to image 2 (maximum correlation)."""
    def score(p):
        w, ok = _warp(p1.image, p1.valid, grid, 1.0 + p[0], p[1], p[2])
        m = ok & p2.valid & region
        return -_ncc(w, p2.image, m)

    # Coarse scan of the stretch (shift 0), then a local polish of all three.
    trial = np.linspace(-0.02, 0.30, 65)
    vals = [score([e, 0.0, 0.0]) for e in trial]
    e0 = trial[int(np.argmin(vals))]
    best = optimize.minimize(score, [e0, 0.0, 0.0], method="Nelder-Mead",
                             options={"xatol": 1e-4, "fatol": 1e-6, "initial_simplex":
                                      [[e0, 0, 0], [e0 + 0.01, 0, 0], [e0, 1.0, 0], [e0, 0, 1.0]]})
    return {"s": 1.0 + float(best.x[0]), "dxi": float(best.x[1]), "deta": float(best.x[2]),
            "ncc": -float(best.fun), "scan_s": (1.0 + trial).tolist(), "scan_ncc": [-v for v in vals]}


def analyse(refresh: bool = False) -> dict:
    cache = config.DERIVED_DIR / "casa_optical.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())
    grid = plates.grid_for(config.CASA)
    gaia = A.read_gaia("casa")
    regs = plates.register_all("casa", grid=grid)
    xi, eta = grid.coords
    r = np.hypot(xi, eta)
    region = (r > SHELL[0]) & (r < SHELL[1])
    prep = {k: M.prepare(regs[k], grid, gaia, region, asinh_scale=50.0 if "ps1" in k else None)
            for k in sorted({k for pair in PAIRS.values() for k in pair})}
    out = {"plates": {k: {"label": regs[k].plate.label, "epoch": regs[k].epoch, "fwhm_arcsec": regs[k].fwhm_arcsec,
                          "solutions": regs[k].solutions} for k in prep}, "pairs": {}}
    # Uncertainty: jackknife over twelve 30-degree wedges of the shell (leave
    # one out, refit), which captures knots brightening or fading unevenly.
    pa = np.degrees(np.arctan2(xi, eta)) % 360.0
    wedge = (pa // 30).astype(int)
    for name, (a, b) in PAIRS.items():
        t1, t2 = regs[a].epoch, regs[b].epoch
        f = fit_stretch(prep[a], prep[b], grid, region)
        jack = np.array([fit_stretch(prep[a], prep[b], grid, region & (wedge != k))["s"] for k in range(12)])
        s = f["s"]
        s_err = float(np.sqrt(11.0 / 12.0 * np.sum((jack - jack.mean()) ** 2)))
        s_vals, cc = scan(prep[a], prep[b], grid, region)
        out["pairs"][name] = {**f, "s_err": s_err, "jackknife_s": jack.tolist(), "t1": t1, "t2": t2,
                              "s_expected": 1.0 + (t2 - t1) / (t1 - T_EXPECTED),
                              "scan_s": s_vals, "scan_ncc": cc}
        print(f"  {name}: best s = {s:.4f} +- {s_err:.4f} (a coasting shell from {T_EXPECTED} "
              f"would give {out['pairs'][name]['s_expected']:.3f})", flush=True)
    out["crab_control"] = crab_control()
    print(f"  Crab control: best s = {out['crab_control']['s_best']:.3f}, expected "
          f"{out['crab_control']['s_expected']:.3f}", flush=True)
    cache.write_text(json.dumps(out, indent=1, default=float))
    return out
