"""The flare of the week, in eight AIA channels.

A flare heats loops of plasma to ~10 MK within minutes; then the loops cool, by
conduction into the chromosphere and by radiation. As they cool they pass through
the temperature bands of AIA's channels one after another, so each channel's light
curve peaks later than the one sensitive to hotter gas: the "cooling sequence"
(e.g. Viall & Klimchuk 2012). Measuring the peak time in each channel against the
temperature it sees traces the flare's plasma cooling from 10 MK to below 1 MK.
"""

from __future__ import annotations

import warnings

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from scipy import ndimage

import config as C
import sdo
import solar_physics as sp

warnings.filterwarnings("ignore")


def load_sequence(paths: C.EpochPaths, manifest: dict, wavelength: int) -> list:
    """All cut-outs of one channel, in DN/s, time-ordered."""
    rels = manifest["flare"]["frames"].get(str(wavelength), [])
    maps = [sdo.load_aia(paths.data / rel) for rel in rels]
    return sorted(maps, key=lambda m: m.date)


def aperture(m, x_arcsec: float, y_arcsec: float, radius_arcsec: float) -> np.ndarray:
    """Pixels within ``radius_arcsec`` of (x, y)."""
    c = m.world_to_pixel(SkyCoord(x_arcsec * u.arcsec, y_arcsec * u.arcsec,
                                  frame=m.coordinate_frame))
    ny, nx = m.data.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    scale = m.scale.axis1.to_value(u.arcsec / u.pix)
    return np.hypot(xx - c.x.value, yy - c.y.value) * scale <= radius_arcsec


def cube(paths: C.EpochPaths, manifest: dict, wavelength: int, begin: pd.Timestamp):
    """(times, frames, pre-flare median frame) for one channel, in DN/s."""
    maps = load_sequence(paths, manifest, wavelength)
    t = pd.DatetimeIndex([pd.Timestamp(m.date.isot, tz="UTC") for m in maps])
    data = np.array([m.data for m in maps])
    pre = np.median(data[t < begin], axis=0) if (t < begin).any() else data[0]
    return t, data, pre, maps


def arcade_mask(paths: C.EpochPaths, manifest: dict, *, level: float = 0.25) -> np.ndarray:
    """The flare's hot loops: pixels whose brightening in 131 A (Fe XXI, ~10 MK) or
    94 A (Fe XVIII, ~7 MK), at its maximum over the flare, exceeds ``level`` of the
    peak brightening near the flare site. Cooler channels are then measured in the
    same pixels, so every light curve follows the same loops as they cool."""
    loc = manifest["flare"]["location"]
    begin = pd.Timestamp(manifest["flare"]["selection"]["begin"])
    mask = None
    for wl in (131, 94):
        t, data, pre, maps = cube(paths, manifest, wl, begin)
        near = aperture(maps[0], loc["hpc_x_arcsec"], loc["hpc_y_arcsec"],
                        C.FLARE_APERTURE_ARCSEC)
        bright = ndimage.gaussian_filter((data - pre).max(axis=0), 1.0)
        m = (bright > level * bright[near].max()) & near
        mask = m if mask is None else (mask | m)
    return mask


def light_curves(paths: C.EpochPaths, manifest: dict, mask: np.ndarray) -> dict:
    """Brightening inside the arcade mask (DN/s, pre-flare image subtracted pixel by
    pixel), per channel."""
    begin = pd.Timestamp(manifest["flare"]["selection"]["begin"])
    out = {}
    for wl in C.FLARE_CHANNELS:
        t, data, pre, _ = cube(paths, manifest, wl, begin)
        excess = np.array([np.nansum((d - pre)[mask]) for d in data])
        base = float(np.nansum(pre[mask]))
        out[wl] = pd.DataFrame({"excess": excess, "relative": excess / base}, index=t)
        out[wl].attrs["background"] = base
    return out


def episode_end(xrs: pd.DataFrame, sel: dict, t_max: pd.Timestamp) -> pd.Timestamp:
    """When the flare's decay stops: the minimum of the 1-8 A flux (5-min running
    mean) between the peak and ``t_max``, if the flux then climbs again by more than
    10% (a new energy release in the region). Otherwise ``t_max``."""
    peak = pd.Timestamp(sel["peak"])
    seg = xrs["long"].loc[peak:t_max].rolling(5, center=True, min_periods=3).mean()
    i = seg.idxmin()
    after = seg.loc[i:]
    if after.max() > 1.10 * seg.loc[i]:
        return i
    return t_max


def peak_time(series: pd.Series) -> tuple[pd.Timestamp, float]:
    """Time of the maximum, refined by a parabola through the three samples around
    it (the cadence is 3 minutes; the peak rarely falls on a sample)."""
    y = series.to_numpy(dtype=float)
    t = (series.index - series.index[0]).total_seconds().to_numpy()
    i = int(np.nanargmax(y))
    if 0 < i < len(y) - 1:
        a, b, c = np.polyfit(t[i - 1:i + 2] - t[i], y[i - 1:i + 2], 2)
        if a < 0:
            dt = -b / (2 * a)
            if abs(dt) <= 1.5 * np.median(np.diff(t)):
                return series.index[i] + pd.Timedelta(seconds=float(dt)), float(c - b**2 / (4 * a))
    return series.index[i], float(y[i])


#: A cooling-phase peak must rise at least this fraction of the channel's full flare
#: amplitude above the post-impulsive dip to count as a peak.
MIN_COOLING_RISE = 0.15


def cooling_sequence(curves: dict, sel: dict, t_end: pd.Timestamp) -> list[dict]:
    """When each channel peaks, relative to the GOES 1-8 A peak.

    Three kinds of channel, three rules, all inside the flare's own episode
    (``t_end`` is where renewed heating begins):

    * the chromospheric channels (1600, 304 A) trace where the energy lands, so
      their peak is the maximum during the X-ray event itself;
    * the flare-hot channels (131 A for Fe XXI at 10 MK, 94 A for Fe XVIII at
      7 MK) peak while the plasma is hottest: their maximum over the episode;
    * the coronal channels (335, 211, 193, 171 A) also flicker during the impulsive
      phase, but their cooling-phase peak is the highest point after the dip that
      follows it. If that point is the end of the episode, or it rises less than
      ``MIN_COOLING_RISE`` of the channel's full flare amplitude above the dip (a
      wiggle, not a peak), the channel had not yet peaked when the episode ended and
      the delay is only a lower limit.
    """
    peak_t = pd.Timestamp(sel["peak"])
    begin, end = pd.Timestamp(sel["begin"]), pd.Timestamp(sel["end"])
    rows = []
    for wl, df in curves.items():
        ch = C.CHANNEL_BY_WAVELENGTH[wl]
        # 131 A answers to Fe XXI (10 MK) in a flare; 193 A's Fe XXIV needs 20 MK,
        # beyond this C-class flare, so it shows Fe XII (1.5 MK) cooling loops
        log_t = ch.log_t[1] if wl == 131 else ch.log_t[0]
        if wl in (1600, 304):
            kind, seg = "chromospheric", df["excess"].loc[begin:end]
        elif wl in (131, 94):
            kind, seg = "hot", df["excess"].loc[begin:t_end]
        else:
            kind = "coronal"
            gradual = df["excess"].loc[end:t_end]
            seg = gradual.loc[gradual.idxmin():]
        tp, value = peak_time(seg)
        lower = kind == "coronal" and (
            seg.index.get_loc(seg.idxmax()) >= len(seg) - 2
            or (seg.max() - seg.iloc[0]) < MIN_COOLING_RISE * df["excess"].max())
        if lower:
            tp = seg.index[-1]
        near_top = seg[seg >= 0.9 * seg.max()]
        rows.append({
            "wavelength": wl, "ion": ch.ions, "log_t": log_t, "kind": kind,
            "peak_time": tp.isoformat(),
            "delay_min": float((tp - peak_t).total_seconds() / 60),
            "lower_limit": bool(lower),
            "plateau_min": [float((near_top.index[0] - peak_t).total_seconds() / 60),
                            float((near_top.index[-1] - peak_t).total_seconds() / 60)],
        })
    return sorted(rows, key=lambda r: r["delay_min"])


def hot_area_cm2(paths: C.EpochPaths, manifest: dict, when: pd.Timestamp,
                 wavelength: int = 131, level: float = 0.5) -> dict:
    """Projected area of the hot flare loops: pixels whose brightening (frame nearest
    ``when`` minus the pre-flare median) exceeds ``level`` of its maximum near the
    flare site."""
    begin = pd.Timestamp(manifest["flare"]["selection"]["begin"])
    t, data, pre, maps = cube(paths, manifest, wavelength, begin)
    k = int(np.argmin(np.abs((t - when).total_seconds())))
    diff = data[k] - pre
    loc = manifest["flare"]["location"]
    near = aperture(maps[k], loc["hpc_x_arcsec"], loc["hpc_y_arcsec"], C.FLARE_APERTURE_ARCSEC)
    sel = near & (diff > level * np.nanmax(diff[near]))
    pix_cm = maps[k].scale.axis1.to_value(u.arcsec / u.pix) * sdo.arcsec_to_mm(maps[k]) * 1e8
    return {"area_cm2": float(sel.sum() * pix_cm**2), "n_pixels": int(sel.sum()),
            "frame_time": maps[k].date.isot, "channel": wavelength, "level": level}


def energetics(thermal: pd.DataFrame, area: dict, sel: dict) -> dict:
    """Density, thermal energy and cooling time from GOES T/EM and the AIA area.

    Volume ~ area^(3/2), the usual order-of-magnitude closure: a factor of a few
    either way.
    """
    thermal = thermal.loc[pd.Timestamp(sel["begin"]):pd.Timestamp(sel["end"])]
    k = thermal["em_cm3"].idxmax()
    t_k = float(thermal.loc[k, "temperature_MK"]) * 1e6
    em = float(thermal.loc[k, "em_cm3"])
    vol = area["area_cm2"] ** 1.5
    n_e = float(np.sqrt(em / vol))
    e_th = float(sp.thermal_energy_erg(t_k, em, vol))
    out = {
        "time_em_max": k.isoformat(), "temperature_MK": t_k / 1e6, "em_cm3": em,
        "volume_cm3": vol, "density_cm3": n_e, "thermal_energy_erg": e_th,
        "radiative_cooling_time_min": float(sp.radiative_cooling_time_s(t_k, n_e) / 60),
    }
    k_t = thermal["temperature_MK"].idxmax()
    out["time_t_max"] = k_t.isoformat()
    out["t_max_MK"] = float(thermal.loc[k_t, "temperature_MK"])
    return out
