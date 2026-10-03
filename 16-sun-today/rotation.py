"""The Sun's rotation, measured from a week of magnetograms.

The Sun is not solid: its equator turns once in about 25 days, higher latitudes
more slowly. Two measurements here:

1. **The rotation profile** (the magnetic-feature method of Snodgrass 1983 and
   Komm, Howard & Harvey 1993). Every magnetogram is remapped onto a grid of
   Carrington longitude and latitude: a frame that itself turns at 14.1844 deg/day,
   the IAU's reference rate. Anything rotating at exactly that rate stands still on
   the grid; faster gas drifts toward larger longitude. Cross-correlating each
   latitude band of the grid between magnetograms a day apart measures that drift
   for the small-scale magnetic network that covers the whole Sun.
2. **Active regions tracked across the disc**, each named by its NOAA number, with
   its own rotation rate and the magnetic polarity of its leading and following
   halves (Hale's law).
"""

from __future__ import annotations

import json
import warnings
from datetime import datetime

import astropy.units as u
import numpy as np
from scipy import ndimage

import config as C
import sdo
import solar_physics as sp

warnings.filterwarnings("ignore")

import sunpy.map  # noqa: E402
from sunpy.coordinates import frames, sun  # noqa: E402

#: Carrington latitude span of the remapped grids (degrees either side of the equator).
GRID_LAT_MAX = 64.0


def _time(m) -> datetime:
    return m.date.to_datetime().replace(tzinfo=None)


def radial_field_map(m):
    """B_r ~ B_los / mu: the line-of-sight field of vertical flux tubes, corrected for
    viewing angle, with the disc beyond ``ROT_MAX_DISK_ANGLE_DEG`` blanked."""
    r = sdo.radius_map(m)
    mu = sdo.mu_from_radius(r)
    br = np.asarray(m.data, dtype=float) / np.where(mu > 0, mu, np.nan)
    br[(mu < np.cos(np.radians(C.ROT_MAX_DISK_ANGLE_DEG))) | (r >= 1.0)] = np.nan
    return sunpy.map.Map(br, m.meta.copy()), mu


def carrington_header(m, lon_start: float, n_lon: int, n_lat: int) -> dict:
    """A plate-carree Carrington grid: pixel i covers longitude lon_start + i*step."""
    step = C.ROT_GRID_DEG
    obs = m.observer_coordinate
    return {
        "naxis": 2, "naxis1": n_lon, "naxis2": n_lat,
        "ctype1": "CRLN-CAR", "ctype2": "CRLT-CAR", "cunit1": "deg", "cunit2": "deg",
        "cdelt1": step, "cdelt2": step,
        "crval1": lon_start, "crpix1": 1.0,
        "crval2": 0.0, "crpix2": (n_lat + 1) / 2.0,
        "date-obs": m.date.isot,
        "hgln_obs": obs.lon.to_value(u.deg), "hglt_obs": obs.lat.to_value(u.deg),
        "dsun_obs": obs.radius.to_value(u.m), "rsun_ref": m.rsun_meters.to_value(u.m),
    }


def continuous_longitude(lon, ref: float):
    """Carrington longitude expressed within 180 degrees of ``ref``, so a week's worth
    of longitudes runs on without jumping from 0 back to 360."""
    return ref + (np.asarray(lon) - ref + 180.0) % 360.0 - 180.0


def remap(m, l0_ref: float) -> dict:
    """One magnetogram on a Carrington grid spanning the usable disc.

    Longitudes are kept continuous around ``l0_ref`` (the central meridian at mid
    week): the Sun's central longitude passes through 0/360 during some weeks, and
    grids on either side of that seam must still line up. Every grid starts on a
    multiple of the grid step, so grids from different times share pixel boundaries
    and can be compared by integer offsets.
    """
    br, _ = radial_field_map(m)
    l0 = float(continuous_longitude(sun.L0(m.date).to_value(u.deg), l0_ref))
    half = C.ROT_MAX_DISK_ANGLE_DEG + 3.0
    step = C.ROT_GRID_DEG
    start_idx = int(np.floor((l0 - half) / step))
    n_lon = int(np.ceil(2 * half / step)) + 1
    n_lat = int(round(2 * GRID_LAT_MAX / step)) + 1
    hdr = carrington_header(m, start_idx * step, n_lon, n_lat)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = br.reproject_to(sunpy.map.Map(np.zeros((n_lat, n_lon)), hdr).wcs)
    return {"time": _time(m), "start_idx": start_idx, "l0": l0,
            "data": np.asarray(out.data, dtype=np.float32)}


def prepare_for_correlation(grid: np.ndarray) -> np.ndarray:
    """Mask active regions and clip strong network so the correlation follows the
    small-scale magnetic network, as in Komm et al. (1993)."""
    g = np.where(np.isfinite(grid), grid, 0.0)
    smooth = ndimage.gaussian_filter(np.abs(g), C.ROT_AR_MASK_SMOOTH_PX)
    active = smooth > C.ROT_AR_MASK_GAUSS
    grow = int(round(C.ROT_AR_MASK_GROW_DEG / C.ROT_GRID_DEG))
    active = ndimage.binary_dilation(active, iterations=grow)
    out = np.clip(g, -C.ROT_FIELD_CLIP_GAUSS, C.ROT_FIELD_CLIP_GAUSS)
    out[active] = 0.0
    out[~np.isfinite(grid)] = np.nan
    return out


def band_rows(n_lat: int) -> list[tuple[float, slice]]:
    """(centre latitude, row slice) for each latitude band."""
    step = C.ROT_GRID_DEG
    lat = (np.arange(n_lat) - (n_lat - 1) / 2) * step
    out = []
    edges = C.ROT_BAND_EDGES_DEG
    for lo, hi in zip(edges[:-1], edges[1:]):
        rows = np.flatnonzero((lat >= lo) & (lat < hi))
        out.append((0.5 * (lo + hi), slice(rows[0], rows[-1] + 1)))
    return out


def correlate_band(a: np.ndarray, b: np.ndarray, max_shift: int):
    """Normalised cross-correlation of two aligned strips against longitude lag.

    ``a`` and ``b`` share columns (same Carrington longitudes); NaN marks pixels
    outside the usable disc. For lag k, a(x) is compared with b(x + k). Returns the
    sub-pixel lag of the peak (parabola through the top three) and the peak value.
    """
    va, vb = np.isfinite(a), np.isfinite(b)
    a0, b0 = np.where(va, a, 0.0), np.where(vb, b, 0.0)
    n = a.shape[1]
    lags = np.arange(-max_shift, max_shift + 1)
    cc = np.full(lags.size, np.nan)
    for i, k in enumerate(lags):
        if k >= 0:
            sa, sb = slice(0, n - k), slice(k, n)
        else:
            sa, sb = slice(-k, n), slice(0, n + k)
        w = va[:, sa] & vb[:, sb]
        if w.sum() < 2000:
            continue
        x, y = a0[:, sa][w], b0[:, sb][w]
        den = np.sqrt(np.sum(x * x) * np.sum(y * y))
        cc[i] = np.sum(x * y) / den if den > 0 else np.nan
    if not np.isfinite(cc).any():
        return np.nan, np.nan
    i = int(np.nanargmax(cc))
    if i == 0 or i == len(lags) - 1 or not np.isfinite(cc[i - 1:i + 2]).all():
        return np.nan, float(cc[i])
    y0, y1, y2 = cc[i - 1], cc[i], cc[i + 1]
    denom = y0 - 2 * y1 + y2
    frac = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    return float(lags[i] + frac), float(y1)


def rotation_profile(grids: list[dict]) -> dict:
    """Rotation rate against latitude from every pair of grids a day apart."""
    step = C.ROT_GRID_DEG
    max_shift = int(round(C.ROT_MAX_SHIFT_DEG / step))
    prepared = [prepare_for_correlation(g["data"]) for g in grids]
    n_lat = grids[0]["data"].shape[0]
    bands = band_rows(n_lat)
    pairs = []
    for i, gi in enumerate(grids):
        for j in range(i + 1, len(grids)):
            dt_h = (grids[j]["time"] - gi["time"]).total_seconds() / 3600
            if abs(dt_h - C.ROT_PAIR_SEPARATION_H) <= 1.6:
                pairs.append((i, j, dt_h))

    rates = {lat: [] for lat, _ in bands}
    peaks = {lat: [] for lat, _ in bands}
    pair_days = []
    for i, j, dt_h in pairs:
        gi, gj = grids[i], grids[j]
        # common longitude columns on the shared grid
        lo = max(gi["start_idx"], gj["start_idx"])
        hi = min(gi["start_idx"] + gi["data"].shape[1], gj["start_idx"] + gj["data"].shape[1])
        if hi - lo < 200:
            continue
        ai = prepared[i][:, lo - gi["start_idx"]:hi - gi["start_idx"]]
        bj = prepared[j][:, lo - gj["start_idx"]:hi - gj["start_idx"]]
        pair_days.append(gi["time"])
        for lat, rows in bands:
            lag, peak = correlate_band(ai[rows], bj[rows], max_shift)
            if np.isfinite(lag) and peak > 0.2:
                rates[lat].append(C.CARRINGTON_RATE_DEG_PER_DAY + lag * step / (dt_h / 24.0))
                peaks[lat].append(peak)

    # Pairs starting 3 h apart share most of their magnetic features, so they are not
    # independent; the independent samples are the separate days covered.
    span_days = (grids[-1]["time"] - grids[0]["time"]).total_seconds() / 86400
    n_indep = max(1, int(np.floor(span_days - C.ROT_PAIR_SEPARATION_H / 24.0)) + 1)
    out_bands = []
    for lat, _ in bands:
        v = np.array(rates[lat])
        if len(v) < 5:
            continue
        med = float(np.median(v))
        good = np.abs(v - med) < 5 * 1.4826 * np.median(np.abs(v - med)) + 1e-9
        vv = v[good]
        out_bands.append({
            "lat": lat, "omega": float(np.mean(vv)), "std": float(np.std(vv, ddof=1)),
            "sem": float(np.std(vv, ddof=1) / np.sqrt(min(n_indep, len(vv)))),
            "n_pairs": int(len(vv)), "n_rejected": int((~good).sum()),
            "median_peak_corr": float(np.median(np.array(peaks[lat])[good])),
        })
    lat = np.array([b["lat"] for b in out_bands])
    om = np.array([b["omega"] for b in out_bands])
    se = np.array([b["sem"] for b in out_bands])
    fits = {}
    for name, fit_c in (("ABC", True), ("AB", False)):
        if len(lat) < (5 if fit_c else 4):
            continue
        coef, cov, chi2, dof = sp.fit_rotation_law(lat, om, se, fit_c=fit_c)
        fits[name] = {"coef": coef.tolist(), "err": np.sqrt(np.diag(cov)).tolist(),
                      "cov": cov.tolist(), "chi2": chi2, "dof": dof}
    # Direct comparison with each published law, band by band (chi^2 with our errors)
    comparisons = {}
    for key, law in C.PUBLISHED_ROTATION.items():
        pred = sp.rotation_law(lat, law["A"], law["B"], law["C"])
        comparisons[key] = {
            "label": law["label"], "tracer": law["tracer"],
            "chi2": float(np.sum(((om - pred) / se) ** 2)), "n": int(len(lat)),
            "mean_diff": float(np.average(om - pred, weights=1 / se**2)),
            "mean_diff_err": float(1 / np.sqrt(np.sum(1 / se**2))),
        }
    eq = [b for b in out_bands if abs(b["lat"]) < 10]
    return {"bands": out_bands, "fits": fits, "comparisons": comparisons,
            "n_pairs": len(pair_days), "n_frames": len(grids), "n_independent_days": n_indep,
            "equator_mean": float(np.average([b["omega"] for b in eq],
                                              weights=[1 / b["sem"] ** 2 for b in eq])) if eq else None}


def synoptic_map(grids: list[dict], half_width_deg: float = 25.0) -> dict:
    """The week's magnetic field on one Carrington map.

    Each magnetogram contributes the strip within ``half_width_deg`` of its central
    meridian, where it sees the surface most squarely; overlapping strips are
    averaged. Longitude increases westward, the direction of rotation.
    """
    step = C.ROT_GRID_DEG
    lo = min(g["start_idx"] for g in grids)
    hi = max(g["start_idx"] + g["data"].shape[1] for g in grids)
    n_lat = grids[0]["data"].shape[0]
    total = np.zeros((n_lat, hi - lo))
    weight = np.zeros_like(total)
    for g in grids:
        cols = g["start_idx"] + np.arange(g["data"].shape[1])
        near = np.abs(cols * step - g["l0"]) <= half_width_deg
        d = g["data"][:, near]
        ok = np.isfinite(d)
        sl = cols[near] - lo
        total[:, sl] += np.where(ok, d, 0.0)
        weight[:, sl] += ok
    with np.errstate(invalid="ignore"):
        data = np.where(weight > 0, total / weight, np.nan)
    lon = (lo + np.arange(hi - lo)) * step
    lat = (np.arange(n_lat) - (n_lat - 1) / 2) * step
    keep = np.isfinite(data).any(axis=0)
    return {"lon": lon[keep], "lat": lat, "data": data[:, keep].astype(np.float32)}


# --- Active regions -----------------------------------------------------------------------------


def detect_regions(m) -> list[dict]:
    """Active regions on one magnetogram: patches of strong field, with their
    positions, areas, fluxes and the centroids of each polarity."""
    data = np.nan_to_num(np.asarray(m.data, dtype=float))
    r = sdo.radius_map(m)
    mu = sdo.mu_from_radius(r)
    usable = (r < 1.0) & (mu > np.cos(np.radians(75)))
    smooth = ndimage.gaussian_filter(np.abs(data), C.AR_DETECT_SMOOTH_PX)
    mask = (smooth > C.AR_DETECT_GAUSS) & usable
    groups, n = ndimage.label(ndimage.binary_dilation(mask, iterations=C.AR_GROUP_DILATE_PX))
    labels = np.where(mask, groups, 0)
    pix_mm = m.scale.axis1.to_value(u.arcsec / u.pix) * sdo.arcsec_to_mm(m)
    msh = sp.msh_per_pixel(pix_mm, np.where(usable, mu, 1.0))
    area_cm2 = (pix_mm * 1e8) ** 2 / np.where(usable, mu, 1.0)
    br = data / np.where(usable, mu, 1.0)
    out = []
    for k in range(1, n + 1):
        sel = labels == k
        area = float(msh[sel].sum())
        if area < C.AR_MIN_AREA_MSH:
            continue
        yy, xx = np.nonzero(sel)
        w = np.abs(br[sel])
        pos = br[sel] > 0
        neg = br[sel] < 0
        cent = {}
        for name, ww, sub in (("all", w, np.ones_like(pos)), ("pos", w, pos), ("neg", w, neg)):
            if ww[sub].sum() <= 0:
                cent[name] = (np.nan, np.nan)
                continue
            cent[name] = (float(np.average(xx[sub], weights=ww[sub])),
                          float(np.average(yy[sub], weights=ww[sub])))
        coords = {}
        for name, (x, y) in cent.items():
            if not np.isfinite(x):
                coords[name] = (np.nan, np.nan)
                continue
            c = m.pixel_to_world(x * u.pix, y * u.pix)
            car = c.transform_to(frames.HeliographicCarrington(obstime=m.date,
                                                               observer=m.observer_coordinate))
            coords[name] = (float(car.lon.to_value(u.deg)), float(car.lat.to_value(u.deg)))
        hpc = m.pixel_to_world(cent["all"][0] * u.pix, cent["all"][1] * u.pix)
        flux_pos = float(np.sum(br[sel][pos] * area_cm2[sel][pos]))
        flux_neg = float(np.sum(br[sel][neg] * area_cm2[sel][neg]))
        out.append({
            "time": m.date.isot, "area_msh": area,
            "lon": coords["all"][0], "lat": coords["all"][1],
            "pos_lon": coords["pos"][0], "pos_lat": coords["pos"][1],
            "neg_lon": coords["neg"][0], "neg_lat": coords["neg"][1],
            "hpc_x": float(hpc.Tx.to_value(u.arcsec)), "hpc_y": float(hpc.Ty.to_value(u.arcsec)),
            "flux_pos_mx": flux_pos, "flux_neg_mx": flux_neg,
            "mu": float(np.average(mu[sel], weights=w)),
        })
    return out


def _dlon(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def link_tracks(detections: list[list[dict]]) -> list[list[dict]]:
    """Chain detections into tracks: in Carrington coordinates a region barely moves
    (under 1 deg/day), so the nearest detection within ``AR_LINK_DEG`` on a later
    frame is the same region. A region may vanish for up to 12 hours (a missed
    frame, or a dip below the threshold) and still keep its track."""
    tracks: list[list[dict]] = []
    for frame in detections:
        taken = set()
        for d in sorted(frame, key=lambda d: -d["area_msh"]):
            t = datetime.fromisoformat(d["time"])
            best, best_dist = None, C.AR_LINK_DEG
            for k, tr in enumerate(tracks):
                if k in taken:
                    continue
                last = tr[-1]
                gap_h = (t - datetime.fromisoformat(last["time"])).total_seconds() / 3600
                if gap_h <= 0 or gap_h > 12.5:
                    continue
                dist = np.hypot(_dlon(d["lon"], last["lon"]) * np.cos(np.radians(d["lat"])),
                                d["lat"] - last["lat"])
                if dist < best_dist:
                    best, best_dist = k, dist
            if best is None:
                tracks.append([d])
                taken.add(len(tracks) - 1)
            else:
                tracks[best].append(d)
                taken.add(best)
    return tracks


def track_summary(track: list[dict], t_ref: datetime) -> dict:
    """Position, size, polarity and rotation rate of one tracked region."""
    t = np.array([(datetime.fromisoformat(d["time"]) - t_ref).total_seconds() / 86400
                  for d in track])
    lon = np.unwrap(np.radians([d["lon"] for d in track]))
    lon = np.degrees(lon)
    lat = np.array([d["lat"] for d in track])
    area = np.array([d["area_msh"] for d in track])
    span_h = (t[-1] - t[0]) * 24
    out = {
        "n": len(track), "span_h": float(span_h),
        "t_first": track[0]["time"], "t_last": track[-1]["time"],
        "lon": float(np.median(lon) % 360), "lat": float(np.median(lat)),
        "max_area_msh": float(area.max()),
        "hemisphere": "N" if np.median(lat) >= 0 else "S",
    }
    if span_h >= C.AR_MIN_TRACK_H and len(track) >= 6:
        (slope, icpt), cov = np.polyfit(t, lon, 1, cov=True)
        out["omega"] = float(C.CARRINGTON_RATE_DEG_PER_DAY + slope)
        out["omega_err"] = float(np.sqrt(cov[0, 0]))
    # polarity of the leading (westward = larger Carrington longitude) half, from the
    # frames where the region is largest
    big = sorted(track, key=lambda d: -d["area_msh"])[: max(3, len(track) // 3)]
    lead = []
    tilt = []
    for d in big:
        if not (np.isfinite(d["pos_lon"]) and np.isfinite(d["neg_lon"])):
            continue
        dl = _dlon(d["pos_lon"], d["neg_lon"]) * np.cos(np.radians(d["lat"]))
        lead.append(1 if dl > 0 else -1)
        # Joy's law tilt: positive when the leading polarity is nearer the equator
        lon_l, lat_l, lon_f, lat_f = ((d["pos_lon"], d["pos_lat"], d["neg_lon"], d["neg_lat"])
                                      if dl > 0 else
                                      (d["neg_lon"], d["neg_lat"], d["pos_lon"], d["pos_lat"]))
        dx = abs(_dlon(lon_l, lon_f)) * np.cos(np.radians(d["lat"]))
        dy = abs(lat_f) - abs(lat_l)
        tilt.append(float(np.degrees(np.arctan2(dy, dx))))
        out.setdefault("separation_deg", float(np.hypot(dx, lat_l - lat_f)))
    both = [d for d in big if np.isfinite(d["pos_lon"]) and np.isfinite(d["neg_lon"])]
    if both:
        ref = both[0]["pos_lon"]
        for key in ("pos_lon", "neg_lon"):
            out[key] = float(np.median([ref + _dlon(d[key], ref) for d in both]) % 360)
        for key in ("pos_lat", "neg_lat"):
            out[key] = float(np.median([d[key] for d in both]))
    if lead:
        out["leading_polarity"] = int(np.sign(np.sum(lead)) or lead[0])
        out["leading_consistency"] = float(np.mean(np.array(lead) == out["leading_polarity"]))
        out["tilt_deg"] = float(np.median(tilt))
        out["flux_balance"] = float(np.median([d["flux_pos_mx"] / -d["flux_neg_mx"]
                                               for d in big if d["flux_neg_mx"] < 0]))
    out["track"] = [{"time": d["time"], "lon": d["lon"], "lat": d["lat"],
                     "area_msh": d["area_msh"], "hpc_x": d["hpc_x"], "hpc_y": d["hpc_y"]}
                    for d in track]
    return out


def noaa_regions(path) -> list[dict]:
    """NOAA's region list: number, position (Carrington), area and class per day."""
    rows = json.loads(open(path).read())
    return [r for r in rows if r.get("carrington_longitude") is not None]


def name_tracks(tracks: list[dict], noaa: list[dict], t0: datetime, t1: datetime) -> None:
    """Attach the NOAA number whose reported position lies nearest each track."""
    cand = [r for r in noaa if t0.date().isoformat() <= r["observed_date"][:10] <= t1.date().isoformat()]
    for tr in tracks:
        best, best_d = None, C.AR_NOAA_MATCH_DEG
        for r in cand:
            d = np.hypot(_dlon(tr["lon"], r["carrington_longitude"]) * np.cos(np.radians(tr["lat"])),
                         tr["lat"] - r["latitude"])
            if d < best_d:
                best, best_d = r, d
        if best is not None:
            tr["noaa"] = int(best["region"])
            tr["noaa_offset_deg"] = float(best_d)
