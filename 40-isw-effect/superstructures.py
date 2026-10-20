"""Stacked ISW measurements on the three superstructure catalogues.

One filter bank holds every filter used anywhere: each catalogue at its
pre-registered filter size plus the exploratory grids. The same bank is applied
to the real sky, to every simulated sky, and to the ISW sky predicted from BOSS
galaxies, so data, noise and prediction are always measured identically.

Bank layout: ``index[(catalog_key, setting)]`` gives the bank rows of one stack,
where ``setting`` is ("R", radius_deg) for a fixed radius or ("alpha", a) for a
radius scaled to each structure (R = a * theta).
"""

from __future__ import annotations

import pickle
from dataclasses import dataclass

import numpy as np

import catalogs as K
import cmb
import stacking as S
from config import DERIVED_DIR, SCAN_ALPHA, SCAN_RADII_DEG

BANK_CACHE = DERIVED_DIR / "filter_bank.pkl"
#: Filters bigger than this are not built (their ring would span half the sky).
MAX_RADIUS_DEG = 30.0

GRANETT_KEYS = ("granett_voids", "granett_clusters")
NH14_KEYS = ("nh14_lrgdim_voids", "nh14_lrgdim_clusters",
             "nh14_lrgbright_voids", "nh14_lrgbright_clusters")
MAO_KEYS = ("mao_voids",)


@dataclass
class Layout:
    bank: S.FilterBank
    index: dict  # (catalog_key, setting) -> rows
    catalogs: dict  # catalog_key -> Catalog

    def rows(self, key, setting) -> np.ndarray:
        return self.index[(key, setting)]


def settings_for(key: str):
    if key in GRANETT_KEYS:
        return [("R", r) for r in SCAN_RADII_DEG]
    return [("alpha", a) for a in SCAN_ALPHA]


def build_layout(*, force: bool = False) -> Layout:
    if BANK_CACHE.exists() and not force:
        with open(BANK_CACHE, "rb") as fh:
            return pickle.load(fh)
    cats = K.all_catalogs()
    d = cmb.prepare_maps()
    masks = S.FastDegrader(d["mask"]).masks
    lon, lat, rad, index = [], [], [], {}
    for key in (*GRANETT_KEYS, *NH14_KEYS, *MAO_KEYS):
        c = cats[key]
        l, b = cmb.equatorial_to_galactic(c.ra, c.dec)
        for setting in settings_for(key):
            r = np.full(len(c), setting[1]) if setting[0] == "R" else setting[1] * c.theta
            start = len(rad)
            lon.extend(l)
            lat.extend(b)
            rad.extend(r)
            index[(key, setting)] = np.arange(start, len(rad))
    rad = np.asarray(rad)
    too_big = rad > MAX_RADIUS_DEG
    bank = S.FilterBank.build(np.asarray(lon), np.asarray(lat), np.where(too_big, 1.0, rad), masks)
    bank.valid[too_big] = False
    bank.radius = rad
    layout = Layout(bank, index, cats)
    with open(BANK_CACHE, "wb") as fh:
        pickle.dump(layout, fh)
    return layout


def stack(values, layout: Layout, key, setting, top_n: int | None = None):
    """Mean filtered temperature of one stack (data vector or sims array)."""
    rows = layout.rows(key, setting)
    if top_n is not None:
        rows = rows[:top_n]
    return S.stack_mean(values, rows)


def granett_combined(values, layout: Layout, radius: float, top_n: int | None = None):
    """Granett et al.'s combined statistic: (mean clusters - mean voids) / 2."""
    v = stack(values, layout, "granett_voids", ("R", radius), top_n)
    c = stack(values, layout, "granett_clusters", ("R", radius), top_n)
    return (c - v) / 2


# --------------------------------------------------------------------------
# Unblinded measurements
# --------------------------------------------------------------------------

PREREG_OF = {"granett_voids": "granett", "granett_clusters": "granett",
             "mao_voids": "mao", **{k: "nh14" for k in NH14_KEYS}}
PRIMARY = {"granett_voids": ("R", 4.0), "granett_clusters": ("R", 4.0), "mao_voids": ("alpha", 0.6),
           **{k: ("alpha", 0.6) for k in NH14_KEYS}}


def real_values(which: str = "smica", layout: Layout | None = None) -> np.ndarray:
    """Every filter of the bank applied to a real CMB map (only once unblinded)."""
    S.assert_unblinded("cmb", "granett", "nh14", "mao")
    layout = layout or build_layout()
    d = cmb.prepare_maps()
    deg = S.FastDegrader(d["mask"])
    return layout.bank.apply(deg(d[which]))


def predicted_values(tag: str = "fiducial", layout: Layout | None = None) -> np.ndarray:
    """Every filter applied to the ISW sky predicted from BOSS galaxies (CMB mask applied)."""
    import healpy as hp
    import boss

    path = DERIVED_DIR / f"pred_bank_{tag}.npy"
    if path.exists():
        return np.load(path)
    layout = layout or build_layout()
    m = hp.ud_grade(boss.predicted_isw(tag), 512).astype(np.float32)
    deg = S.FastDegrader(cmb.prepare_maps()["mask"])
    v = layout.bank.apply(deg(m))
    np.save(path, v)
    return v


def matched_filter(data, pred, sims, rows):
    """A_ISW = sum(p T / s^2) / sum(p^2 / s^2) over the given filter rows.

    Returns the estimate, its null distribution over the simulations (whose
    scatter is the error, correlations between overlapping filters included)
    and the S/N expected if LCDM is right (A = 1).
    """
    rows = rows[np.isfinite(data[rows]) & np.isfinite(pred[rows])]
    s2 = np.var(sims[:, rows], axis=0, ddof=1)
    w = pred[rows] / s2
    norm = np.sum(w * pred[rows])
    a = float(np.sum(w * data[rows]) / norm)
    a_sims = (sims[:, rows] @ w) / norm
    sigma = float(np.std(a_sims, ddof=1))
    return {"A": a, "sigma": sigma, "snr": a / sigma, "expected_snr_lcdm": 1.0 / sigma,
            "n": int(len(rows)), "A_sims": a_sims}


# --------------------------------------------------------------------------
# Look-elsewhere: what tuning R and N does to the significance
# --------------------------------------------------------------------------


def granett_scan(values, layout: Layout):
    """Combined Granett statistic on the (radius x top-N) grid: array (..., nR, nN)."""
    from config import SCAN_N

    out = np.stack([np.stack([granett_combined(values, layout, r, n) for n in SCAN_N], axis=-1)
                    for r in SCAN_RADII_DEG], axis=-2)
    return out


def look_elsewhere(data, sims, layout: Layout) -> dict:
    """Significance at the published choice versus the best choice on the grid.

    The S/N of each grid cell uses that cell's own simulated scatter. The
    global p-value asks how often a pure-noise sky reaches, somewhere on the
    grid, an S/N at least as large as the best one found in the data.
    """
    from config import SCAN_N

    d = granett_scan(data, layout)
    s = granett_scan(sims, layout)
    sd = s.std(axis=0, ddof=1)
    snr_d = d / sd
    snr_s = s / sd
    best = float(np.nanmax(snr_d))
    ib = np.unravel_index(np.nanargmax(snr_d), snr_d.shape)
    max_s = np.nanmax(snr_s.reshape(len(s), -1), axis=1)
    i4, i50 = SCAN_RADII_DEG.index(4.0), SCAN_N.index(50)
    fixed = float(snr_d[i4, i50])
    p_fixed = float((1 + np.sum(snr_s[:, i4, i50] >= fixed)) / (len(s) + 1))
    p_global = float((1 + np.sum(max_s >= best)) / (len(s) + 1))
    return {"snr_grid": snr_d, "value_grid": d, "sigma_grid": sd,
            "fixed_snr": fixed, "p_fixed": p_fixed,
            "best_snr": best, "best_radius": SCAN_RADII_DEG[ib[0]], "best_n": SCAN_N[ib[1]],
            "p_global": p_global, "max_snr_sims": max_s}


# --------------------------------------------------------------------------
# Errors from random positions on the real sky
# --------------------------------------------------------------------------


def random_positions(n: int, seed: int = 7):
    """Uniform random points in the BOSS footprint (north or south cap), Galactic (l, b)."""
    import healpy as hp
    import boss

    foot = boss.footprint()
    comp = foot["comp"].max(axis=0)
    pix = np.flatnonzero(comp > 0)
    rng = np.random.default_rng(seed)
    pick = rng.choice(pix, size=n, replace=True)
    ra, dec = hp.pix2ang(boss.BOSS_MASK_NSIDE, pick, lonlat=True)
    ra = ra + rng.uniform(-0.5, 0.5, n) * np.degrees(hp.nside2resol(boss.BOSS_MASK_NSIDE))
    l, b = cmb.equatorial_to_galactic(ra, dec)
    cap_north = dec > -15 + 0 * ra
    cap_north &= (ra > 90) & (ra < 290)
    return l, b, cap_north


def random_position_errors(key: str, setting, layout: Layout, *, n_pool: int | None = None,
                           n_sets: int = 5000, seed: int = 7) -> dict:
    """The catalogue's filters placed at random points of the real CMB map.

    Each random point gets the filter radius of a randomly chosen structure of
    the catalogue, and random catalogues of the same size are drawn from the
    pool. Positions are drawn in the same cap(s) as the catalogue's structures.
    """
    from config import N_RANDOM_POSITIONS

    S.assert_unblinded("cmb", PREREG_OF[key])
    n_pool = n_pool or N_RANDOM_POSITIONS
    cat = layout.catalogs[key]
    rows = layout.rows(key, setting)
    radii = layout.bank.radius[rows]
    l, b, north = random_positions(4 * n_pool, seed)
    cat_north = (cat.ra > 90) & (cat.ra < 290)
    if cat_north.all():
        sel = north
    elif (~cat_north).all():
        sel = ~north
    else:
        sel = np.ones_like(north)
    l, b = l[sel][:n_pool], b[sel][:n_pool]
    rng = np.random.default_rng(seed + 1)
    r = rng.choice(radii[np.isfinite(radii)], size=len(l))
    d = cmb.prepare_maps()
    deg = S.FastDegrader(d["mask"])
    bank = S.FilterBank.build(l, b, np.minimum(r, MAX_RADIUS_DEG), deg.masks)
    vals = bank.apply(deg(d["smica"]))
    vals = vals[np.isfinite(vals)]
    n = int(np.sum(np.isfinite(layout.bank.valid[rows]) & layout.bank.valid[rows]))
    means = np.array([rng.choice(vals, n, replace=False).mean() for _ in range(n_sets)])
    return {"single_rms": float(vals.std()), "set_mean": float(means.mean()),
            "set_sigma": float(means.std(ddof=1)), "n_pool": int(len(vals)), "n_per_set": n,
            "set_means": means}


# --------------------------------------------------------------------------
# Stacked images
# --------------------------------------------------------------------------

IMAGES_CACHE = DERIVED_DIR / "stack_images.npz"


def highpass(m, mask, lcut: int = 10, width: int = 4):
    """Remove multipoles below ~lcut (cosine taper over +-width), as Granett et al. 2015 did.

    Fluctuations larger than ~20 deg shift a whole stacked patch up or down and
    hide the structure-sized signal; the compensated filter ignores them anyway.
    """
    import healpy as hp

    lmax = 3 * hp.get_nside(m) - 1
    ell = np.arange(lmax + 1)
    fl = np.clip((ell - (lcut - width)) / (2 * width), 0, 1)
    fl = 0.5 - 0.5 * np.cos(np.pi * fl)
    alm = hp.map2alm(np.where(mask, m, 0.0), lmax=lmax)
    return hp.alm2map(hp.almxfl(alm, fl), hp.get_nside(m), lmax=lmax).astype(np.float32)


def stacked_images(*, force: bool = False) -> dict:
    """Mean CMB patches around the structures (real sky and BOSS-predicted ISW sky).

    Granett: 24 x 24 deg patches. Mao and NH14: patches of +-3 structure radii,
    so voids of different sizes line up.
    """
    if IMAGES_CACHE.exists() and not force:
        d = np.load(IMAGES_CACHE)
        return {k: d[k] for k in d.files}
    import healpy as hp
    import boss
    from config import STACK_IMAGE_DEG, STACK_IMAGE_PIX

    S.assert_unblinded("cmb", "granett", "nh14", "mao")
    layout = build_layout()
    d = cmb.prepare_maps()
    pred = hp.ud_grade(boss.predicted_isw("fiducial"), 512).astype(np.float32)
    skies = {"smica": highpass(d["smica"], d["mask"]), "pred": highpass(pred, d["mask"])}
    out = {}
    for key in (*GRANETT_KEYS, *NH14_KEYS, *MAO_KEYS):
        c = layout.catalogs[key]
        lon, lat = cmb.equatorial_to_galactic(c.ra, c.dec)
        signs = np.ones(len(c))
        scaled = not key.startswith("granett")
        for sky, m in skies.items():
            out[f"{key}__{sky}"] = S.stack_image(
                m, d["mask"], lon, lat, signs, size_deg=6.0 if scaled else STACK_IMAGE_DEG,
                npix=STACK_IMAGE_PIX, scale=c.theta if scaled else None)
    # individual Granett patches (coarser) for the interactive explorer
    for key in GRANETT_KEYS:
        c = layout.catalogs[key]
        lon, lat = cmb.equatorial_to_galactic(c.ra, c.dec)
        out[f"{key}__patches"] = np.stack([S.gnomonic_patch(d["smica"], d["mask"], lon[i], lat[i],
                                                            STACK_IMAGE_DEG, 49, 512)
                                           for i in range(len(c))]).astype(np.float32)
    np.savez(IMAGES_CACHE, **out)
    return out
