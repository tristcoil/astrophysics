"""Plate 1 of the 1920 paper: identify its stars, find the Moon, and date the exposure.

The paper's only photograph (Phil. Trans. A 220, Plate 1) is a half-tone of
"one of the negatives taken with the 4-inch lens at Sobral", with each
measured star marked by a pair of short horizontal dashes. The scan on the
Internet Archive (leaf 42 of philtrans06337895; Public Domain Mark) is good
enough to do three things:

1. find the seven dash pairs (a thin-line filter, then pairs of segments on the
   same row) and the star between each pair;
2. identify them: an affine map from the stars' Gaia DR3 positions, carried
   back to 1919 with their proper motions, onto the seven marks, tried for
   every assignment of the seven 4-inch stars (one assignment fits to a pixel;
   the next best is ~40 times worse);
3. fit a circle to the white disc (the Moon: no light reached the negative
   there) and compare its centre with JPL Horizons' Moon as seen from Sobral:
   the Moon slides across the Sun by ~0.4'/minute, so its offset dates the
   exposure to within seconds (its middle: the disc is the Moon's mean place
   over the 28-s exposure).
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

import config as C
import stars as S

PLATE1_LEAF = 42
#: The seven 4-inch stars (Table I numbers), in the order of the paper's tables.
STARS_4IN = [11, 5, 4, 3, 6, 10, 2]
#: Rough positions (px, full-resolution leaf 42) of the seven dash pairs and of
#: the Moon's disc, read by eye from a half-size view; refined by the code.
APPROX_MARKS = [(1600, 380), (1660, 720), (1440, 1352), (1376, 1420), (360, 1466), (780, 2836), (2230, 2836)]
APPROX_MOON = (902.0, 1830.0, 330.0)


def image() -> np.ndarray:
    from PIL import Image

    return np.asarray(Image.open(C.PAGES_DIR / f"leaf{PLATE1_LEAF:03d}.jpg").convert("L")).astype(float)


def find_marks(im: np.ndarray | None = None) -> pd.DataFrame:
    """Star positions (px) midway between the inner ends of each dash pair."""
    from scipy import ndimage as ndi

    im = image() if im is None else im
    ink = 255.0 - im
    along = ndi.uniform_filter1d(ink, 15, axis=1)
    resp = along - 0.5 * (np.roll(along, 3, axis=0) + np.roll(along, -3, axis=0))   # thin horizontal lines
    rows = []
    for k, (x0, y0) in enumerate(APPROX_MARKS):
        X0, Y0 = x0 - 170, y0 - 60
        r = resp[Y0:y0 + 60, X0:x0 + 170]
        lab, n = ndi.label(r > 0.4 * r.max(), structure=np.ones((3, 3)))
        segs = []
        for j in range(1, n + 1):
            yy, xx = np.nonzero(lab == j)
            if xx.max() - xx.min() >= 20:
                segs.append((X0 + xx.min(), X0 + xx.max(), Y0 + np.average(yy, weights=r[yy, xx])))
        best = None
        for a, b in itertools.combinations(segs, 2):
            left, right = (a, b) if a[0] < b[0] else (b, a)
            gap = right[0] - left[1]
            if abs(a[2] - b[2]) < 4 and 15 < gap < 110:
                score = abs(a[2] - b[2]) + abs(gap - 45) / 20
                if best is None or score < best[0]:
                    best = (score, left, right, gap)
        if best is None:
            raise RuntimeError(f"no dash pair near {(x0, y0)}")
        _, left, right, gap = best
        rows.append({"mark": k, "x_px": (left[1] + right[0]) / 2, "y_px": (left[2] + right[2]) / 2,
                     "dash_px": (left[1] - left[0] + right[1] - right[0]) / 2, "gap_px": gap})
    return pd.DataFrame(rows)


def tangent_arcmin(ra_rad, dec_rad, ra0: float, dec0: float) -> tuple[np.ndarray, np.ndarray]:
    xi, eta = S.standard_coords(np.asarray(ra_rad), np.asarray(dec_rad), ra0, dec0)
    return np.rad2deg(xi) * 60.0, np.rad2deg(eta) * 60.0


def identify(marks: pd.DataFrame, centre: tuple[float, float], utc: str = "1919-05-29T12:00:00") -> dict:
    """Best assignment of the seven 4-inch stars to the marks, by an affine fit
    from tangent-plane coordinates (arcmin, about ``centre`` = (ra, dec) rad)."""
    g = S.gaia_table1().loc[STARS_4IN]
    pos = S.gaia_at(g, utc)
    xs, ys = tangent_arcmin(pos.ra.rad, pos.dec.rad, *centre)
    A = np.column_stack([xs, ys, np.ones(len(xs))])
    P = marks[["x_px", "y_px"]].to_numpy()
    trials = []
    for perm in itertools.permutations(range(len(P))):
        B = P[list(perm)]
        sol, *_ = np.linalg.lstsq(A, B, rcond=None)
        trials.append((float(np.sqrt(np.mean((A @ sol - B) ** 2))), perm, sol))
    trials.sort(key=lambda t: t[0])
    rms, perm, sol = trials[0]
    resid = A @ sol - P[list(perm)]
    M = sol[:2].T
    return {"affine": sol, "rms_px": rms, "next_rms_px": trials[1][0],
            "scale_px_per_arcmin": float(np.sqrt(abs(np.linalg.det(M)))),
            # pixels run right/down and (xi, eta) east/north, so an unmirrored sky
            # image (east left, north up) has a positive determinant
            "mirrored": bool(np.linalg.det(M) < 0),
            "stars": pd.DataFrame({"no": STARS_4IN, "mark": list(perm),
                                   "x_px": P[list(perm), 0], "y_px": P[list(perm), 1],
                                   "resid_x_px": resid[:, 0], "resid_y_px": resid[:, 1]}).set_index("no")}


def fit_moon(im: np.ndarray | None = None) -> dict:
    """Circle fitted to the edge of the white disc (steepest bright-to-dark step
    along 360 radial profiles, outliers clipped)."""
    from scipy import ndimage as ndi
    from scipy.optimize import least_squares

    im = image() if im is None else im
    sm = ndi.gaussian_filter(im, 2.0)
    cx0, cy0, r0 = APPROX_MOON
    yy, xx = np.mgrid[0:im.shape[0], 0:im.shape[1]]
    inside = (sm > 200) & ((xx - cx0) ** 2 + (yy - cy0) ** 2 < (0.9 * r0) ** 2)
    cx0, cy0 = xx[inside].mean(), yy[inside].mean()
    pts = []
    for th in np.deg2rad(np.arange(0.0, 360.0, 1.0)):
        rr = np.arange(r0 - 60, r0 + 60, 0.5)
        x, y = cx0 + rr * np.cos(th), cy0 + rr * np.sin(th)
        grad = -np.gradient(ndi.map_coordinates(sm, [y, x], order=1))
        k = int(np.argmax(grad))
        if grad[k] > 3.0:
            pts.append((x[k], y[k]))
    pts = np.array(pts)

    def res(p, q):
        return np.hypot(q[:, 0] - p[0], q[:, 1] - p[1]) - p[2]

    sol = least_squares(res, [cx0, cy0, r0], args=(pts,), loss="soft_l1", f_scale=2.0)
    r = res(sol.x, pts)
    keep = np.abs(r) < 3 * 1.4826 * np.median(np.abs(r))
    sol = least_squares(res, sol.x, args=(pts[keep],))
    rr = sol.fun
    cov = np.linalg.inv(sol.jac.T @ sol.jac) * (rr @ rr) / (len(rr) - 3)
    return {"cx": sol.x[0], "cy": sol.x[1], "r_px": sol.x[2], "sigma": np.sqrt(np.diag(cov)),
            "rms_px": float(np.sqrt(np.mean(rr ** 2))), "n_edge": int(keep.sum())}


def moon_sun_track() -> pd.DataFrame:
    """Horizons' Moon and Sun seen from Sobral through totality (fetch_data.py)."""
    df = pd.read_csv(C.HORIZONS_DIR / "moon_sun_sobral_1919.csv")
    df["t"] = pd.to_datetime(df["utc"], format="%Y-%b-%d %H:%M:%S.%f")
    return df


def sobral_exposures(path=None, instrument: str = "4-inch") -> pd.DataFrame:
    """Start, middle and end (UT) of each exposure at Sobral on 29 May 1919 (p.300). The
    times were counted by a metronome from second contact, which the expedition assumed
    to be 11:58:18 UT; JPL Horizons puts it at 11:58:10 (see analysis.step_contacts)."""
    path = C.TRANSCRIBED_DIR / "dyson1920_sobral_exposures.csv" if path is None else path
    ex = pd.read_csv(path, comment="#", skipinitialspace=True)
    ex = ex[ex["instrument"] == instrument].copy()
    ex["start"] = pd.to_datetime(ex["utc_civil_derived"])
    ex["end"] = ex["start"] + pd.to_timedelta(ex["exposure_s"], unit="s")
    ex["mid"] = ex["start"] + (ex["end"] - ex["start"]) / 2
    ex["ref"] = ex["plate_label_tables2to8"].fillna(ex["ref_no"].astype(str)) if instrument == "4-inch" else ex["ref_no"]
    return ex[["ref_no", "ref", "start", "mid", "end", "exposure_s", "note"]].reset_index(drop=True)


def match_exposure(when: pd.Timestamp, sigma_t: float, clock_offset_s: float = 0.0) -> dict:
    """Which 4-inch exposure is the plate? The white disc is where the Moon covered the plate for
    the whole 28-s exposure (it moved ~4 px, 12", meanwhile), so its fitted centre is the Moon's
    mid-exposure position and the date (when, sigma_t) estimates the middle of an exposure. For
    each exposure that shows stars: z = (middle - date) / sigma_t, and a probability from Gaussian
    likelihoods with equal priors. The clock offset shifts the 1919 times (e.g. -8 s if second
    contact was at JPL's 11:58:10 rather than the assumed 11:58:18)."""
    ex = sobral_exposures()
    ex = ex[ex["ref"] != "VI"].copy()           # taken through cloud, shows no stars (p.300)
    z = ((ex["mid"] - when).dt.total_seconds() + clock_offset_s) / sigma_t
    like = np.exp(-0.5 * z ** 2)
    ex = ex.assign(z=z, prob=like / like.sum())
    best = ex.loc[ex["prob"].idxmax()]
    return {"plate": str(best["ref"]), "ref_no": int(best["ref_no"]), "prob": float(best["prob"]), "z": float(best["z"]),
            "start": best["start"].strftime("%H:%M:%S"), "mid": best["mid"].strftime("%H:%M:%S"),
            "end": best["end"].strftime("%H:%M:%S"),
            "z_by_plate": {str(r): float(q) for r, q in zip(ex["ref"], ex["z"])},
            "probs": {str(r): float(q) for r, q in zip(ex["ref"], ex["prob"]) if q > 1e-6}}


def date_exposure(ident: dict, moon: dict, centre: tuple[float, float], track: pd.DataFrame | None = None) -> dict:
    """When was the Moon where the plate shows it? Predicted Moon centre (px)
    against time, from the star solution; the time of closest approach to the
    measured centre, with the miss distance and the along-track speed."""
    track = moon_sun_track() if track is None else track
    sol = ident["affine"]

    def to_px(ra_deg, dec_deg):
        x, y = tangent_arcmin(np.deg2rad(ra_deg), np.deg2rad(dec_deg), *centre)
        return np.column_stack([x, y, np.ones_like(x)]) @ sol

    m = to_px(track["moon_ra"].to_numpy(), track["moon_dec"].to_numpy())
    s = to_px(track["sun_ra"].to_numpy(), track["sun_dec"].to_numpy())
    t = (track["t"] - track["t"].iloc[0]).dt.total_seconds().to_numpy()
    tt = np.linspace(t[0], t[-1], 20001)
    mx, my = np.interp(tt, t, m[:, 0]), np.interp(tt, t, m[:, 1])
    d = np.hypot(mx - moon["cx"], my - moon["cy"])
    k = int(np.argmin(d))
    speed = float(np.hypot(np.gradient(mx, tt), np.gradient(my, tt))[k])     # px/s
    when = track["t"].iloc[0] + pd.Timedelta(float(tt[k]), unit="s")
    radius_pred = float(np.interp(tt[k], t, track["moon_diam_arcsec"].to_numpy()) / 120.0) * ident["scale_px_per_arcmin"]
    return {"utc": when, "miss_px": float(d[k]), "speed_px_per_s": speed,
            "sigma_t_s": float(np.hypot(ident["rms_px"], d[k]) / speed),
            "moon_radius_pred_px": radius_pred, "sun_px": s, "moon_px": m, "t_track": track["t"]}
