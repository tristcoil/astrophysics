"""Cluster photometry and the external inputs: distances, reddenings, metallicities.

* ``read_hugs``: the HUGS catalogue of one cluster (method 1), cached as .npz.
* ``quality_mask``: the selection of well-measured stars the HUGS team
  recommends (Nardiello et al. 2018, Sect. 3, after Milone et al. 2012): in
  magnitude bins, reject stars whose photometric scatter (RMS) lies more than
  3.5 sigma above the typical value, or whose PSF-fit quality (QFIT) lies 3.5 sigma
  below it, and keep only star-like shapes, |RADXS| < 0.15.
* ``members``: proper-motion membership probability >= 90%.
* ``differential_reddening``: dust in front of a cluster is patchy, so stars in
  one corner of the field can be dimmer and redder than in another. Following
  Milone et al. (2012), each star's offset from the main-sequence ridge line,
  measured along the reddening direction and averaged over its nearest
  well-measured neighbours, gives the local excess Delta E(B-V); subtracting it
  sharpens the sequences. Only *relative* reddening is removed here: the
  absolute E(B-V) of the cluster still comes from Harris.
* ``harris``, ``bv21``: Harris (2010 edition) and Baumgardt & Vasiliev (2021).
"""

from __future__ import annotations

import re

import numpy as np
from scipy.spatial import cKDTree

import config as C

BANDS = ("F275W", "F336W", "F438W", "F606W", "F814W")
_FIELDS = ("m", "rms", "qfit", "radxs", "nf", "ng")


def hugs_path(hugs_id: str):
    return C.HUGS_DIR / f"hlsp_hugs_hst_wfc3-uvis-acs-wfc_{hugs_id}_multi_v1_catalog-meth{C.HUGS_METHOD}.txt"


def read_hugs(hugs_id: str) -> dict:
    """All columns of one HUGS catalogue as a dict of arrays (cached)."""
    cache = C.DERIVED_DIR / f"hugs_{hugs_id}.npz"
    if cache.exists():
        z = np.load(cache)
        return {k: z[k] for k in z.files}
    import pandas as pd

    names = ["x", "y"] + [f"{f}_{b}" for b in BANDS for f in _FIELDS] + ["memb", "ra", "dec", "id", "iter"]
    df = pd.read_csv(hugs_path(hugs_id), comment="#", sep=r"\s+", header=None, names=names)
    out = {k: df[k].to_numpy() for k in names}
    out["id"] = out["id"].astype(str)  # identifiers like "R0000001": keep as text, not objects
    for b in BANDS:  # -99.9999 marks "not measured"
        bad = out[f"m_{b}"] < -50
        for f in ("m", "rms"):
            out[f"{f}_{b}"] = np.where(bad, np.nan, out[f"{f}_{b}"])
    np.savez_compressed(cache, **out)
    return out


def _clipped_stats(v: np.ndarray, nsig: float = 3.5, iters: int = 5) -> tuple[float, float]:
    keep = np.isfinite(v)
    for _ in range(iters):
        mu, sd = np.mean(v[keep]), np.std(v[keep])
        if sd == 0:
            break
        new = keep & (abs(v - mu) <= nsig * sd)
        if new.sum() == keep.sum():
            break
        keep = new
    return float(np.mean(v[keep])), float(np.std(v[keep]))


def quality_mask(cat: dict, bands=("F606W", "F814W"), nbins: int = 12, nsig: float = 3.5,
                 return_lines: bool = False):
    """Well-measured stars in every band of ``bands`` (Nardiello et al. 2018 recipe).

    Two adaptations: stars bright enough to saturate the long exposures are
    measured in a single short exposure (N_f = 1), where the RMS is 0 by
    construction, so the RMS cut applies to stars with N_f >= 2 only. And the
    shape cut |RADXS| < 0.15 also removes RADXS = -1, a value given to ~15% of
    stars, mostly in the crowded centre, whose colours scatter twice as much as
    the rest (blends).
    """
    ok = np.ones(len(cat["x"]), bool)
    lines = {}
    for b in bands:
        m, rms, q, sh = cat[f"m_{b}"], cat[f"rms_{b}"], cat[f"qfit_{b}"], cat[f"radxs_{b}"]
        nf = cat[f"nf_{b}"]
        good = np.isfinite(m) & (nf > 0) & (nf < 99)  # 99 marks saturated stars
        multi = good & (nf >= 2) & np.isfinite(rms)
        lo, hi = np.nanpercentile(m[good], [0.2, 99.8])
        edges = np.linspace(lo, hi, nbins + 1)
        c_r, rms_lim, c_q, q_lim = [], [], [], []
        for a, z in zip(edges[:-1], edges[1:]):
            s = good & (m >= a) & (m < z)
            if s.sum() >= 20:
                mu_q, sd_q = _clipped_stats(q[s], nsig)
                c_q.append(np.mean(m[s]))
                q_lim.append(mu_q - nsig * sd_q)
            s = multi & (m >= a) & (m < z)
            if s.sum() >= 20:
                mu_r, sd_r = _clipped_stats(rms[s], nsig)
                c_r.append(np.mean(m[s]))
                rms_lim.append(mu_r + nsig * sd_r)
        rl = np.interp(m, c_r, rms_lim)
        ql = np.interp(m, c_q, q_lim)
        rms_ok = np.where(nf >= 2, rms <= rl, True)
        ok &= good & rms_ok & (q >= ql) & (abs(sh) < 0.15)
        lines[b] = (np.array(c_r), np.array(rms_lim), np.array(c_q), np.array(q_lim))
    return (ok, lines) if return_lines else ok


def members(cat: dict, pmin: float = 90.0) -> np.ndarray:
    """Proper-motion members (probabilities are in per cent; -1 = not measured)."""
    return cat["memb"] >= pmin


# ---------------------------------------------------------------------------
# Differential reddening (Milone et al. 2012, Sect. 3)
# ---------------------------------------------------------------------------


def _ridge(x: np.ndarray, y: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Median of x in bins of y (sigma-clipped), as a smooth ridge line x(y)."""
    yc, xc = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (y >= a) & (y < b)
        if s.sum() < 15:
            continue
        v = x[s]
        med = np.median(v)
        mad = 1.4826 * np.median(abs(v - med))
        v = v[abs(v - med) < 3 * max(mad, 1e-3)]
        xc.append(np.median(v))
        yc.append(np.median(y[s]))
    return np.array(yc), np.array(xc)


def differential_reddening(mag: np.ndarray, col: np.ndarray, x: np.ndarray, y: np.ndarray,
                           ref: np.ndarray, a_mag: float, e_col: float, *,
                           n_neigh: int = 50, iterations: int = 3, mag_range=None) -> np.ndarray:
    """Local Delta E(B-V) for every star, from reference main-sequence stars.

    ``mag`` and ``col`` are the CMD (e.g. F606W and F606W - F814W); ``a_mag`` and
    ``e_col`` are A_mag / E(B-V) and E(col) / E(B-V), which fix the direction of
    the reddening vector. ``ref`` flags the reference stars (well-measured members
    on the main sequence below the turnoff, where the sequence runs at a large
    angle to the reddening vector). Returns Delta E(B-V) per star.
    """
    theta = np.arctan2(a_mag, e_col)  # angle of the reddening vector in (col, mag)
    unit = np.hypot(a_mag, e_col)  # its length per unit E(B-V)
    dE = np.zeros(len(mag))
    tree = cKDTree(np.column_stack([x[ref], y[ref]]))
    ref_idx = np.flatnonzero(ref)
    for _ in range(iterations):
        m_c = mag - a_mag * dE
        c_c = col - e_col * dE
        # rotate so the abscissa runs along the reddening vector
        absc = c_c * np.cos(theta) + m_c * np.sin(theta)
        ordi = -c_c * np.sin(theta) + m_c * np.cos(theta)
        lo, hi = (np.percentile(ordi[ref], [1, 99]) if mag_range is None else mag_range)
        edges = np.arange(lo, hi + 0.1, 0.1)
        yc, xc = _ridge(absc[ref], ordi[ref], edges)
        resid_ref = absc[ref] - np.interp(ordi[ref], yc, xc)
        # median offset of the nearest reference stars (excluding the star itself)
        dist, nn = tree.query(np.column_stack([x, y]), k=n_neigh + 1)
        self_hit = ref_idx[nn] == np.arange(len(mag))[:, None]
        nn_resid = np.where(self_hit, np.nan, resid_ref[nn])
        local = np.nanmedian(nn_resid[:, :n_neigh + 1], axis=1)
        dE = dE + local / unit
    return dE - np.median(dE[ref])


# ---------------------------------------------------------------------------
# External tables
# ---------------------------------------------------------------------------


def harris() -> dict:
    """Harris (2010 edition): {ID: {"name", "ra", "dec", "l", "b", "feh", "ebv", "v_hb", "dm_v"}}."""
    lines = (C.LIT_DIR / "mwgc.dat").read_text().splitlines()
    parts = [i for i, ln in enumerate(lines) if ln.strip().startswith("Part")]
    out = {}

    def num(s):
        s = s.strip()
        return float(s) if s else np.nan

    for ln in lines[parts[0]:parts[1]]:
        m = re.match(r"^\s(\S.{10})(.{12})\s+(\d\d \d\d [\d.]+)\s+([-+]\d\d \d\d [\d.]+)\s+([-\d.]+)\s+([-\d.]+)", ln)
        if not m:
            continue
        gid = m.group(1).strip()
        ra = [float(v) for v in m.group(3).split()]
        de = m.group(4).split()
        sign = -1 if de[0].startswith("-") else 1
        dec = sign * (abs(float(de[0])) + float(de[1]) / 60 + float(de[2]) / 3600)
        out[gid] = {"name": m.group(2).strip(), "ra": 15 * (ra[0] + ra[1] / 60 + ra[2] / 3600),
                    "dec": dec, "l": float(m.group(5)), "b": float(m.group(6))}
    for ln in lines[parts[1]:parts[2]]:
        gid = ln[:12].strip()
        if gid not in out or len(ln) < 40:
            continue
        out[gid].update({"feh": num(ln[12:18]), "ebv": num(ln[21:28]),
                         "v_hb": num(ln[28:34]), "dm_v": num(ln[34:40])})
    return out


def _pm(tex: str) -> tuple[float, float, float] | None:
    """Parse '$ 8.501^{+0.071}_{-0.070}$' or '$8.627 \\pm 0.473$' -> (value, minus, plus)."""
    tex = tex.replace("\\,", "").replace(" ", "")
    m = re.search(r"([\d.]+)\^\{\+([\d.]+)\}_\{-([\d.]+)\}", tex)
    if m:
        return float(m.group(1)), float(m.group(3)), float(m.group(2))
    m = re.search(r"([\d.]+)\\pm([\d.]+)", tex)
    if m:
        return float(m.group(1)), float(m.group(2)), float(m.group(2))
    return None


BV21_COLUMNS = ("gaia_parallax", "gaia_kinematic", "hst_kinematic", "subdwarf", "star_counts",
                "literature", "mean")


def bv21() -> dict:
    """Baumgardt & Vasiliev (2021) Table A? distances (kpc): {name: {column: (d, minus, plus)}}."""
    text = (C.LIT_DIR / "bv21_findis_tab.tex").read_text()
    out = {}
    for ln in text.splitlines():
        if "&" not in ln or ln.lstrip().startswith(("\\", "%", "Name", "&")):
            continue
        cells = [c.strip() for c in ln.split("\\\\")[0].split("&")]
        if len(cells) < 9:
            continue
        row = {col: _pm(cell) for col, cell in zip(BV21_COLUMNS, cells[1:8])}
        row["n_total"] = int(cells[8]) if cells[8].strip().isdigit() else None
        out[cells[0]] = row
    return out


def distance_modulus(entry: tuple[float, float, float]) -> tuple[float, float]:
    """(mu, sigma_mu) from a distance in kpc with asymmetric errors."""
    d, lo, hi = entry
    mu = 5 * np.log10(d * 100.0)
    sig = 5 / np.log(10) * 0.5 * (lo + hi) / d
    return float(mu), float(sig)


def geometric_distance(row: dict) -> tuple[float, float, float] | None:
    """Inverse-variance mean of the purely geometric/dynamical distances only.

    Gaia EDR3 parallax and the Gaia and HST kinematic distances (proper-motion
    dispersion against line-of-sight velocity dispersion) use no stellar-evolution
    model. The literature average B&V21 fold in includes isochrone fits of the
    CMD, which would make an age measured with that distance partly circular.
    """
    vals = [row[k] for k in ("gaia_parallax", "gaia_kinematic", "hst_kinematic") if row.get(k)]
    if not vals:
        return None
    d = np.array([v[0] for v in vals])
    s = np.array([0.5 * (v[1] + v[2]) for v in vals])
    w = 1 / s**2
    mean = float(np.sum(w * d) / np.sum(w))
    err = float(np.sqrt(1 / np.sum(w)))
    return mean, err, err


def hugs_to_ngc(hugs_id: str) -> str:
    """'ngc0104' -> 'NGC 104', the naming of Harris and Baumgardt & Vasiliev."""
    return f"NGC {int(hugs_id[3:])}"


def cluster_inputs(name: str) -> dict:
    """Everything the age fit needs to know about one cluster, before the data.

    ``name`` is one of the four core clusters ("M92", ...) or any HUGS id ("ngc0288").
    """
    info = C.CLUSTERS.get(name) or {"hugs": name, "harris": hugs_to_ngc(name), "bv21": hugs_to_ngc(name)}
    h = harris()[info["harris"]]
    b = bv21()[info["bv21"]]
    mu, mu_err = distance_modulus(b["mean"])
    geo = geometric_distance(b)
    out = {
        "name": name,
        "hugs": info["hugs"],
        "feh": h["feh"],
        "feh_err": C.FEH_SIGMA,
        "ebv": h["ebv"],
        "ebv_err": max(C.EBV_MIN_SIGMA, C.EBV_REL_SIGMA * h["ebv"]),
        "alpha": C.ALPHA_PRIOR[0],
        "alpha_err": C.ALPHA_PRIOR[1],
        "dist_kpc": b["mean"][0],
        "dist_err_kpc": 0.5 * (b["mean"][1] + b["mean"][2]),
        "dm": mu,
        "dm_err": float(np.hypot(mu_err, C.DM_SYS_FLOOR)),
        "dm_err_quoted": mu_err,
        "bv21": b,
        "harris": h,
    }
    if geo:
        out["dm_geo"], out["dm_geo_err"] = distance_modulus(geo)
        out["dist_geo_kpc"] = geo[0]
    return out
