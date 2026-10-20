"""The Local Bubble: its dust walls traced in the Edenhofer et al. (2024) map.

We follow O'Neill et al. (2024, ApJ 973, 136) step by step, so their released
model is a check on our pipeline:

* every Nside-256 pixel centre is a line of sight; the density is sampled at
  1 pc steps from 69 to 650 pc (at a pixel centre the healpy weights reduce to
  the pixel itself, so this is the log-linear radial interpolation that
  dustmaps returns; O'Neill et al. truncated the map at 650 pc);
* the profile is smoothed with a Gaussian of sigma = 7 pc;
* peaks = local maxima with prominence > 2.04e-6 ZGR23/pc and a width of at
  least 3 pc at half prominence (scipy's find_peaks);
* the bubble's wall along that sightline is the *first* peak; its inner and
  outer edges are where the smoothed profile falls to the peak value minus
  half the prominence (and stop at a neighbouring peak's valley);
* the wall's gas: n = 1653 cm^-3 per ZGR23 mag/pc (A_G/N_H = 4e-22 cm^2 mag),
  mass = 1.37 m_p n dV summed over the unsmoothed samples between the edges.

Beyond the reproduction we fit a sphere to the wall (to compare with the
spherical expansion model of Zucker et al. 2022) and ask whether young star
clusters sit on the wall more often than older ones.
"""

from __future__ import annotations

import multiprocessing as mp

import healpy as hp
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths

import catalogs
import config as C
from dustmap import DustSphere, lbd_to_xyz

D_GRID = np.arange(69.0, 651.0, 1.0)     # pc, sampling of each line of sight
SIGMA_SMOOTH = 7.0                        # pc
P_MIN = 2.04e-6                           # ZGR23 mag/pc, minimum prominence
MIN_WIDTH = 3.0                           # pc, at half prominence
BATCH = 16384


N_OUT = 10


def _first_peak(prof_s: np.ndarray, prof: np.ndarray) -> tuple:
    """Wall properties along one smoothed profile (and its unsmoothed twin).

    Edges at half prominence (O'Neill's fiducial A'_0.5) and at 0.9 of the
    prominence (their generous A'_0.9, used for the interior volume).
    """
    peaks, props = find_peaks(prof_s, prominence=P_MIN, width=MIN_WIDTH, rel_height=0.5)
    if peaks.size == 0:
        return (np.nan,) * (N_OUT - 1) + (0,)
    p = peaks[0]
    out = [D_GRID[p]]
    for rel in (0.5, 0.9):
        if rel == 0.5:
            left, right = props["left_ips"][0], props["right_ips"][0]
        else:
            _, _, left, right = peak_widths(prof_s, peaks[:1], rel_height=rel,
                                            prominence_data=(props["prominences"][:1],
                                                             props["left_bases"][:1],
                                                             props["right_bases"][:1]))
            left, right = left[0], right[0]
        # a second peak inside the window: stop at the valley between them
        if peaks.size > 1 and peaks[1] < right:
            right = min(right, p + np.argmin(prof_s[p:peaks[1] + 1]))
        i0, i1 = int(np.ceil(left)), int(np.floor(right))
        seg = prof[i0:i1 + 1]
        d = D_GRID[i0:i1 + 1]
        out += [D_GRID[0] + left, D_GRID[0] + right, (seg * d**2).sum()]
        if rel == 0.5:
            out += [seg.max(), seg.sum()]          # peak density; extinction (dr = 1 pc)
    return tuple(out) + (len(peaks),)


def _batch(args) -> np.ndarray:
    pix, = args
    prof = _SPHERE.column_profiles(pix, D_GRID)                   # (n_d, n_pix)
    prof_s = gaussian_filter1d(prof, SIGMA_SMOOTH, axis=0)
    out = np.empty((pix.size, N_OUT))
    for j in range(pix.size):
        out[j] = _first_peak(prof_s[:, j], prof[:, j])
    return out


_SPHERE: DustSphere | None = None


def shell_model(force: bool = False, workers: int = 6) -> pd.DataFrame:
    """One row per Nside-256 sightline (NESTED order), cached as parquet."""
    path = C.CACHE_DIR / "bubble_shell.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    global _SPHERE
    _SPHERE = DustSphere("main")
    _ = _SPHERE.log_mean                                          # load before forking
    npix = _SPHERE.npix
    batches = [(np.arange(s, min(s + BATCH, npix)),) for s in range(0, npix, BATCH)]
    ctx = mp.get_context("fork")
    with ctx.Pool(workers, initializer=_one_thread) as pool:
        parts = pool.map(_batch, batches, chunksize=1)
    res = np.concatenate(parts)
    l, b = hp.pix2ang(_SPHERE.nside, np.arange(npix), nest=True, lonlat=True)
    df = pd.DataFrame(res, columns=["d_peak", "d_inner", "d_outer", "mass_w", "rho_peak", "a_zgr",
                                    "d_inner09", "d_outer09", "mass_w09", "n_peaks"])
    df["l"], df["b"] = l, b
    domega = 4 * np.pi / npix
    df["n_peak"] = C.ZGR23_TO_NH * df["rho_peak"]
    df["A_G"] = C.ZGR23_TO_AG * df["a_zgr"]
    to_msun = C.MU_H * C.M_PROTON_G * C.ZGR23_TO_NH * domega * C.PC_CM**3 / C.M_SUN_G
    df["mass"] = to_msun * df["mass_w"]
    df["mass09"] = to_msun * df["mass_w09"]
    df["thick"] = df["d_outer"] - df["d_inner"]
    df["x"], df["y"], df["z"] = lbd_to_xyz(l, b, df["d_peak"])
    df = df.drop(columns=["mass_w", "mass_w09"])
    df.to_parquet(path)
    _SPHERE = None
    return df


def _one_thread():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------

def compare_oneill(model: pd.DataFrame) -> dict:
    """Sightline-by-sightline comparison with O'Neill et al.'s released table."""
    ref = catalogs.oneill_shell()
    pix = hp.ang2pix(256, ref["l"].to_numpy(), ref["b"].to_numpy(), nest=True, lonlat=True)
    ours = model.iloc[pix].reset_index(drop=True)
    both = ref["d"].notna() & ours["d_peak"].notna()
    dd = (ours["d_peak"] - ref["d"])[both]
    stats = {
        "n_sightlines": int(len(ref)), "n_ref_with_peak": int(ref["d"].notna().sum()),
        "n_ours_with_peak": int(model["d_peak"].notna().sum()),
        "frac_same_peak_within_1pc": float((dd.abs() <= 1).mean()),
        "frac_within_5pc": float((dd.abs() <= 5).mean()),
        "median_dd": float(dd.median()),
        "d_inner_median_diff": float((ours["d_inner"] - ref["d_inner"])[both].median()),
        "d_outer_median_diff": float((ours["d_outer"] - ref["d_outer"])[both].median()),
    }
    for name, ser in (("ours", model["d_peak"]), ("oneill", ref["d"])):
        stats[f"d_peak_percentiles_{name}"] = [float(v) for v in np.nanpercentile(ser, [0, 2.28, 50, 97.72, 100])]
    for name, ser in (("ours", model["thick"]), ("oneill", ref["thick"])):
        stats[f"thick_percentiles_{name}"] = [float(v) for v in np.nanpercentile(ser, [0, 2.28, 50, 97.72, 100])]
    for name, ser in (("ours", model["n_peak"]), ("oneill", ref["npeak"])):
        stats[f"n_peak_percentiles_{name}"] = [float(v) for v in np.nanpercentile(ser, [0, 2.28, 50, 97.72, 100])]
    stats["mass_total_ours"] = float(np.nansum(model["mass"]))
    stats["mass_total_oneill_table"] = float(np.nansum(ref["mass"]))
    return stats


def sphere_fit(model: pd.DataFrame, rmax: float = 400.0) -> dict:
    """Least-squares sphere through the wall's peak points (sightlines weighted equally).

    Sightlines whose wall lies beyond ``rmax`` (the open 'chimney' toward the
    north) are left out, since a sphere cannot describe an opening.
    """
    from scipy.optimize import least_squares
    ok = model["d_peak"].notna() & (model["d_peak"] < rmax)
    P = model.loc[ok, ["x", "y", "z"]].to_numpy()

    def resid(p):
        return np.linalg.norm(P - p[:3], axis=1) - p[3]

    fit = least_squares(resid, x0=[0, 0, 0, 170.0], loss="soft_l1", f_scale=30.0)
    r = resid(fit.x)
    return {"centre_pc": [float(v) for v in fit.x[:3]], "radius_pc": float(fit.x[3]),
            "rms_pc": float(np.sqrt(np.mean(r**2))), "mad_pc": float(np.median(np.abs(r))),
            "n_sightlines": int(ok.sum()), "rmax_pc": rmax}


def interior_volume(model: pd.DataFrame) -> dict:
    """Volume inside the inner edge (incl. the 69-pc sphere), and its equivalent radius.

    O'Neill et al. define the interior with the generous A'_0.9 inner edge;
    sightlines without a wall (open to 650 pc) count to 650 pc.
    """
    npix = len(model)
    domega = 4 * np.pi / npix
    out = {"open_sightlines": int(model["d_peak"].isna().sum())}
    for col, tag in (("d_inner09", "a09"), ("d_inner", "a05")):
        d_in = model[col].fillna(650.0).to_numpy()
        vol = float(np.sum(domega * d_in**3 / 3))
        out[f"volume_pc3_{tag}"] = vol
        out[f"equivalent_radius_pc_{tag}"] = float((3 * vol / (4 * np.pi)) ** (1 / 3))
    out["mass_a09_msun"] = float(np.nansum(model["mass09"]))
    return out


def clusters_on_wall(model: pd.DataFrame, rmax: float = 400.0) -> pd.DataFrame:
    """Hunt & Reffert clusters within ``rmax`` with their offset from the wall along the sightline."""
    hr = catalogs.hunt_reffert()
    hr = hr[hr["reliable"] & (hr["dist50"] < rmax)].copy()
    pix = hp.ang2pix(256, hr["GLON"].to_numpy(), hr["GLAT"].to_numpy(), nest=True, lonlat=True)
    wall = model.iloc[pix]
    hr["d_wall"] = wall["d_peak"].to_numpy()
    hr["d_wall_inner"] = wall["d_inner"].to_numpy()
    hr["d_wall_outer"] = wall["d_outer"].to_numpy()
    hr["offset"] = hr["dist50"] - hr["d_wall"]
    return hr
